"""Node 19 (v7.0, P4): Batch variability — DSD-quantile-resolved L2-L4.

Runs L2 → L3 → L4 across a set of representative bead radii sampled from
the L1 droplet-size distribution and reports distribution-level outputs
weighted by mass fraction in each bin.

Why this matters (per Scientific Advisor in v7.0 plan §3 / §4):
  Real microsphere batches have a DSD; chromatography pressure scales
  ~1/d_p² and film-mass-transfer ~1/d_p, so M3 predictions hand-fed a
  single d50 are wrong for any application that cares about packing.
  Batch variability is the largest scientific gap remaining for
  chromatography use.

Design:
  - **Backward-compatible**: ``PipelineOrchestrator.run_single`` is
    unchanged. ``run_batch`` is a new top-level convenience that calls
    ``run_single`` per quantile and aggregates the results.
  - **Mass-weighted aggregation**: for extensive properties (porosity,
    G_DN) we report the volume-weighted mean across quantiles; for
    intensive distributions we report p5/p50/p95 across quantiles
    (NOT inside-quantile MC — that's Node 18's UQ engine).
  - **Quantile choice**: default is 5 mass-volume quantile representatives
    (d10, d25, d50, d75, d90) per Doc 10 §3.3. Caller can override with
    a custom diameter list.

Out of scope for v7.0:
  - Full DSD discretisation (N_bins > 5) — current default is the
    minimum-viable resolution; v7.1 batch QC reporting will go finer.
  - Per-quantile UQ — calling Node 18's MC engine inside each quantile
    is straightforward but multiplies runtime by N_quantiles × N_samples.
    Defer to v7.1 when batch reporting is the primary use case.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, Sequence

import numpy as np

from ..datatypes import FullResult, RunContext, SimulationParameters
from .orchestrator import PipelineOrchestrator

logger = logging.getLogger(__name__)


@dataclass
class QuantileRun:
    """One representative-bead simulation result inside a batch.

    Audit N3 (v7.0.1, foot-gun callout): ``full_result.emulsification`` is
    SHARED BY REFERENCE with the base L1 run for performance reasons —
    every quantile sees the same DSD because L1 is run only once. So
    ``full_result.emulsification.d50`` does NOT match the
    ``representative_radius_m`` used to drive L2/L3/L4 in this quantile.

    Downstream consumers (M2 ingest, plotting, exports) MUST use
    ``representative_radius_m`` (or its diameter equivalent
    ``representative_diameter_m``) when they need the bead size that
    actually fed L2/L3/L4 for this quantile.
    """
    quantile: float           # e.g. 0.10 for d10
    representative_radius_m: float
    mass_fraction: float      # contribution of this bin to batch mass
    full_result: FullResult

    @property
    def representative_diameter_m(self) -> float:
        """Diameter equivalent of ``representative_radius_m`` (audit N3 helper).

        Use this instead of ``full_result.emulsification.d50`` for any
        downstream code that needs THIS quantile's bead size — the
        full_result.emulsification fields reflect the BASE L1 DSD, not
        the per-quantile representative.
        """
        return 2.0 * self.representative_radius_m


@dataclass
class BatchResult:
    """Aggregated batch-level result across DSD quantiles.

    The per-quantile FullResults are kept on the result for users who
    want to inspect specific bins; the headline statistics summarise
    the batch as a whole.
    """
    quantile_runs: list[QuantileRun] = field(default_factory=list)
    # Mass-weighted means
    mean_d32_m: float = 0.0
    mean_pore_m: float = 0.0
    mean_porosity: float = 0.0
    mean_G_DN_Pa: float = 0.0
    # Distribution-level p5/p50/p95 across quantiles (intensive properties)
    pore_p5_m: float = 0.0
    pore_p50_m: float = 0.0
    pore_p95_m: float = 0.0
    G_DN_p5_Pa: float = 0.0
    G_DN_p50_Pa: float = 0.0
    G_DN_p95_Pa: float = 0.0
    # Batch composition tracking
    n_quantiles: int = 0
    quantile_radii_m: list[float] = field(default_factory=list)
    quantile_mass_fractions: list[float] = field(default_factory=list)


def _representative_radii_from_dsd(
    full_result: FullResult,
    quantiles: Sequence[float],
) -> tuple[np.ndarray, np.ndarray]:
    """Extract representative diameters from the L1 PBE result.

    Uses the volume-weighted CDF over ``full_result.emulsification.d_bins``
    to pick the percentile diameters and returns the mass fraction in the
    bin centred on each chosen percentile. The mass fractions sum to 1
    (modulo numerical roundoff).
    """
    d = np.asarray(full_result.emulsification.d_bins, dtype=float)
    n = np.asarray(full_result.emulsification.n_d, dtype=float)
    # Volume-weighted density
    v = n * (d ** 3)
    if v.sum() <= 0:
        # Degenerate — fall back to d50 for every requested quantile.
        d50 = float(full_result.emulsification.d50)
        return (
            np.full(len(quantiles), d50 / 2.0),
            np.full(len(quantiles), 1.0 / len(quantiles)),
        )
    v = v / v.sum()
    cdf = np.cumsum(v)
    # For each requested quantile, find the bin where the CDF first
    # exceeds it. Mass fraction is approximated as the slice of the CDF
    # bracketing that quantile (uniform partition between requested points).
    radii = np.zeros(len(quantiles))
    for i, q in enumerate(quantiles):
        idx = int(np.searchsorted(cdf, q))
        idx = max(0, min(idx, d.size - 1))
        radii[i] = d[idx] / 2.0
    # Mass fraction per quantile bin: midpoint partition of [0,1].
    sorted_q = np.sort(quantiles)
    edges = np.concatenate(([0.0], (sorted_q[:-1] + sorted_q[1:]) / 2.0, [1.0]))
    mass = np.diff(edges)
    return radii, mass


def run_batch(
    params: SimulationParameters,
    quantiles: Sequence[float] = (0.10, 0.25, 0.50, 0.75, 0.90),
    db=None,
    output_dir: Optional[Path] = None,
    run_context: Optional[RunContext] = None,
    crosslinker_key: str = "genipin",
    l2_mode: str = "empirical",
) -> BatchResult:
    """Run the full pipeline on a DSD-quantile-resolved batch.

    Algorithm:
      1. Run L1 once to obtain the droplet-size distribution.
      2. Pick representative diameters at the requested quantiles.
      3. Re-run L2 → L3 → L4 with each representative bead radius
         (skip L1 — it would re-converge to the same DSD).
      4. Aggregate mass-weighted means + per-quantile percentiles.

    The current implementation calls ``run_single`` once per quantile
    for clarity; for performance a per-quantile L2-L4-only path could be
    added in v7.1 to avoid re-running L1 each time. The runtime cost is
    ``N_quantiles × T_run_single``.
    """
    quantiles = tuple(quantiles)
    if not quantiles:
        raise ValueError("At least one quantile must be requested.")
    if not all(0.0 < q < 1.0 for q in quantiles):
        raise ValueError("All quantiles must be in (0, 1).")
    # Audit N8 (v7.0.1): mass-fraction edges only make sense for sorted,
    # unique quantiles. Silently sort+dedupe rather than producing
    # negative or zero mass fractions for duplicate or unsorted input.
    quantiles = tuple(sorted(set(quantiles)))

    orch = PipelineOrchestrator(db=db, output_dir=output_dir)
    # First call gives us the DSD; we use its representative radii for
    # subsequent runs (each "run" sees the SAME L1 output but a different
    # bead-size assumption flowing into L2 onwards).
    base_run = orch.run_single(
        params, run_context=run_context,
        crosslinker_key=crosslinker_key, l2_mode=l2_mode,
    )
    radii, mass_fracs = _representative_radii_from_dsd(base_run, quantiles)

    quantile_runs: list[QuantileRun] = []
    for q, R, w in zip(quantiles, radii, mass_fracs):
        # Override d50 by mutating the orchestrator's per-quantile call:
        # the cleanest way is to re-run with phi_d unchanged but force the
        # representative R into L2/L3/L4 by using run_single with a
        # post-L1 patched FullResult. The orchestrator does not expose
        # that hook in v6.1, so we fall back to the simplest correct
        # approach: construct a per-quantile FullResult by re-using the
        # base L1 output and re-running L2-L4 with the override.
        per_q = _run_l2_l4_at_R(orch, params, db, base_run, R, run_context,
                                 crosslinker_key, l2_mode)
        quantile_runs.append(QuantileRun(
            quantile=float(q),
            representative_radius_m=float(R),
            mass_fraction=float(w),
            full_result=per_q,
        ))

    # Aggregate
    pores = np.array([r.full_result.gelation.pore_size_mean for r in quantile_runs])
    porosities = np.array([r.full_result.gelation.porosity for r in quantile_runs])
    g_dns = np.array([r.full_result.mechanical.G_DN for r in quantile_runs])
    weights = np.array([r.mass_fraction for r in quantile_runs])
    weights_norm = weights / weights.sum() if weights.sum() > 0 else np.ones_like(weights) / len(weights)

    return BatchResult(
        quantile_runs=quantile_runs,
        mean_d32_m=float(base_run.emulsification.d32),
        mean_pore_m=float(np.average(pores, weights=weights_norm)),
        mean_porosity=float(np.average(porosities, weights=weights_norm)),
        mean_G_DN_Pa=float(np.average(g_dns, weights=weights_norm)),
        pore_p5_m=float(np.percentile(pores, 5)),
        pore_p50_m=float(np.percentile(pores, 50)),
        pore_p95_m=float(np.percentile(pores, 95)),
        G_DN_p5_Pa=float(np.percentile(g_dns, 5)),
        G_DN_p50_Pa=float(np.percentile(g_dns, 50)),
        G_DN_p95_Pa=float(np.percentile(g_dns, 95)),
        n_quantiles=len(quantile_runs),
        quantile_radii_m=[float(r) for r in radii],
        quantile_mass_fractions=[float(w) for w in mass_fracs],
    )


def _run_l2_l4_at_R(
    orch: PipelineOrchestrator,
    params: SimulationParameters,
    db,
    base_run: FullResult,
    R_droplet: float,
    run_context,
    crosslinker_key: str,
    l2_mode: str,
) -> FullResult:
    """Re-run L2-L4 at a specific bead radius, reusing the base L1 output.

    Cheaper than calling run_single again because L1 PBE is the dominant
    cost in the default config. Returns a FullResult whose L1 fields are
    inherited from base_run (so downstream consumers see the unchanged
    DSD) but whose L2/L3/L4 reflect the specified R_droplet.
    """
    import copy as _copy
    from ..level2_gelation.solver import solve_gelation, solve_gelation_timing
    from ..level3_crosslinking.solver import solve_crosslinking
    from ..level4_mechanical.solver import solve_mechanical
    from ..datatypes import FullResult as _FR

    props = orch.db.update_for_conditions(
        T_oil=params.formulation.T_oil,
        c_agarose=params.formulation.c_agarose,
        c_chitosan=params.formulation.c_chitosan,
        c_span80=params.formulation.c_span80,
    )
    # Re-run L2 timing + pore at the new radius
    timing = solve_gelation_timing(params, props, R_droplet=R_droplet)
    gel = solve_gelation(params, props, R_droplet=R_droplet,
                          mode=l2_mode, timing=timing)
    xl = solve_crosslinking(
        params, props, R_droplet=R_droplet, porosity=gel.porosity,
        crosslinker_key=crosslinker_key,
    )
    mech = solve_mechanical(params, props, gel, xl, R_droplet=R_droplet)

    return _FR(
        parameters=base_run.parameters,
        emulsification=base_run.emulsification,
        gelation=gel,
        crosslinking=xl,
        mechanical=mech,
        gelation_timing=timing,
        run_report=base_run.run_report,  # share top-level evidence
    )
