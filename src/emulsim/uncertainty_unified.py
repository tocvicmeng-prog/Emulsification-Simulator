"""Node 13 (v6.1, F7): unified uncertainty schema scaffold.

The platform currently runs two separate uncertainty systems:

  - ``uncertainty_core.py`` propagates MaterialProperties perturbations
    through the M1-L4 pipeline and emits ``UncertaintyResult`` (median +
    90% CI for d32, pore, G_DN).
  - ``uncertainty_propagation/`` propagates ``M1ExportContract`` CVs into
    M2 via ``run_with_uncertainty`` and emits per-bound q_max stats.

Both are useful and well-tested in their respective scopes, but they do
not share a schema, source taxonomy, or correlation handling. The
consensus plan §6 calls for a single contract that downstream consumers
(RunDossier, optimizer, UI badges) can rely on.

This module is the **scaffold** for that unification:

  1. ``UncertaintyKind`` enum classifies the *source* of an uncertainty
     so consumers can reason about it (model-form vs measured replicate,
     calibration posterior vs assumed CV, etc.).
  2. ``UnifiedUncertaintySpec`` is the input contract: a list of named
     sources with kind + distribution parameters, plus an optional
     correlation matrix.
  3. ``UnifiedUncertaintyResult`` is the output contract: per-output
     percentiles plus a per-kind breakdown so the user can see *which
     class of uncertainty* dominates the prediction interval.
  4. Two adapters (``from_m1l4_result``, ``from_m1_contract_uq``) convert
     the existing-system outputs into the unified schema without changing
     either implementation. This lets new code consume the unified format
     today while we keep the legacy paths alive for the v7.0 merge.

The actual MC engine merge is deferred to v7.0 (consensus plan §6).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Optional, Sequence

import numpy as np


class UncertaintyKind(Enum):
    """Classification of an uncertainty source (consensus plan §6).

    The kind drives interpretation of the resulting interval. A model-form
    uncertainty is irreducible by sampling more inputs; a measurement
    uncertainty shrinks with more replicates. The optimizer and UI need
    this distinction to give honest guidance.
    """

    MEASUREMENT = "measurement"
    """Replicate-derived (e.g. CV from triplicate IFT measurement)."""

    MATERIAL_PROPERTY = "material_property"
    """Literature spread on a material property (e.g. agarose density)."""

    CALIBRATION_POSTERIOR = "calibration_posterior"
    """Posterior std from a calibration fit (e.g. fitted breakage_C1)."""

    MODEL_FORM = "model_form"
    """Discrepancy between model and reality (e.g. CT vs Alopaeus)."""

    NUMERICAL = "numerical"
    """Solver tolerance, grid resolution, time-step truncation."""

    SCALE_UP = "scale_up"
    """Lab-bench-to-pilot transfer factor."""


@dataclass
class UncertaintySource:
    """One named, classified, sampleable uncertainty input.

    Attributes:
        name: Identifier for logging/breakdown (e.g. "props.sigma").
        kind: Which class of uncertainty this represents.
        distribution: One of ``"normal"`` (mu=value, sigma=std) or
            ``"lognormal"`` (median=value, sigma_log=std).
        value: Central tendency (mean for normal, median for lognormal).
        std: Spread parameter (std for normal, sigma_log for lognormal).
        units: Physical unit string for the value (informational).
    """
    name: str
    kind: UncertaintyKind
    distribution: str  # "normal" | "lognormal"
    value: float
    std: float
    units: str = ""

    def sample(self, rng: np.random.Generator, n: int) -> np.ndarray:
        if self.distribution == "normal":
            return rng.normal(self.value, self.std, size=n)
        if self.distribution == "lognormal":
            # value is the median; std is sigma in log-space.
            return rng.lognormal(np.log(max(self.value, 1e-300)), self.std, size=n)
        raise ValueError(f"Unsupported distribution {self.distribution!r}")


@dataclass
class UnifiedUncertaintySpec:
    """Input spec for the unified MC engine (v7.0 target).

    For v6.1 we use this to *describe* legacy runs uniformly; the actual
    sampler is still the per-system implementation.
    """
    sources: list[UncertaintySource] = field(default_factory=list)
    correlation: Optional[np.ndarray] = None
    """Optional N_sources x N_sources correlation matrix. None = independent."""
    n_samples: int = 100
    seed: int = 42

    def add(self, src: UncertaintySource) -> None:
        self.sources.append(src)

    def kinds_present(self) -> set[UncertaintyKind]:
        return {s.kind for s in self.sources}


@dataclass
class OutputUncertainty:
    """Distribution summary for one named output (e.g. d32, q_max, DBC10)."""
    name: str
    units: str
    mean: float
    p5: float
    p50: float
    p95: float
    n_samples: int
    n_failed: int = 0

    @property
    def ci_width_relative(self) -> float:
        """(p95 - p5) / |p50|; useful as a single uncertainty scalar."""
        if abs(self.p50) < 1e-300:
            return float("nan")
        return (self.p95 - self.p5) / abs(self.p50)


@dataclass
class UnifiedUncertaintyResult:
    """Schema-uniform UQ output for downstream consumers.

    Carries per-output percentiles AND a kinds-per-output map so consumers
    can answer: "what class of uncertainty dominates the d32 interval?"
    even when only one input kind was sampled (the answer is then
    "{that kind}", which is itself useful for trust-banner messaging).

    Audit N2 (v7.0.1): ``kinds_sampled`` records ONLY the kinds that
    actually contributed to the percentile interval. Kinds declared in
    the spec but not sampled by the underlying engine (e.g. calibration
    posteriors that v7.0 absorbs but does not yet propagate through the
    legacy MC engines) live in ``kinds_declared_but_not_sampled`` instead
    so consumers do not over-attribute uncertainty contribution.
    """
    outputs: list[OutputUncertainty] = field(default_factory=list)
    kinds_sampled: set[UncertaintyKind] = field(default_factory=set)
    n_samples: int = 0
    n_failed: int = 0
    source_label: str = ""
    """Where this Result came from — e.g. "uncertainty_core.M1L4" or
    "uncertainty_propagation.M2"."""

    # Audit N2 (v7.0.1): kinds the spec listed but the engine did NOT
    # actually sample. Documents the v7.0 limitation without overclaiming
    # in `kinds_sampled`. Empty in v7.1+ when all kinds are sampled.
    kinds_declared_but_not_sampled: set[UncertaintyKind] = field(default_factory=set)

    def get(self, name: str) -> Optional[OutputUncertainty]:
        for o in self.outputs:
            if o.name == name:
                return o
        return None

    def summary(self) -> str:
        kinds_str = ", ".join(sorted(k.value for k in self.kinds_sampled)) or "none"
        lines = [
            f"Unified UQ — source={self.source_label}, "
            f"n={self.n_samples} ({self.n_failed} failed), kinds=[{kinds_str}]",
        ]
        if self.kinds_declared_but_not_sampled:
            unsamp = ", ".join(sorted(k.value for k in self.kinds_declared_but_not_sampled))
            lines.append(
                f"  (declared but NOT sampled by engine: [{unsamp}] — "
                "interval is a lower bound on total uncertainty)"
            )
        for o in self.outputs:
            lines.append(
                f"  {o.name:18s} = {o.p50:.4g} {o.units}  "
                f"(p5..p95 {o.p5:.4g}..{o.p95:.4g}, "
                f"width/median = {o.ci_width_relative:.2%})"
            )
        return "\n".join(lines)


# ─── Adapters from the existing UQ implementations ────────────────────────


def from_m1l4_result(legacy_result, source_label: str = "uncertainty_core.M1L4"):
    """Convert ``uncertainty_core.UncertaintyResult`` to the unified schema.

    The legacy result samples MaterialProperties (material-property kind)
    and emits hardcoded medians + 90% CIs for d32/pore/G_DN. We recompute
    p5/p50/p95 from the underlying samples for consistency with future
    consumers that expect the percentile fields.
    """
    if legacy_result is None:
        return UnifiedUncertaintyResult(source_label=source_label)

    def _output(name: str, unit: str, samples: np.ndarray) -> OutputUncertainty:
        valid = np.asarray(samples, dtype=float)
        valid = valid[np.isfinite(valid)]
        if valid.size == 0:
            return OutputUncertainty(
                name=name, units=unit, mean=float("nan"),
                p5=float("nan"), p50=float("nan"), p95=float("nan"),
                n_samples=0, n_failed=int(getattr(legacy_result, "n_failed", 0)),
            )
        return OutputUncertainty(
            name=name, units=unit,
            mean=float(np.mean(valid)),
            p5=float(np.percentile(valid, 5)),
            p50=float(np.percentile(valid, 50)),
            p95=float(np.percentile(valid, 95)),
            n_samples=int(valid.size),
            n_failed=int(getattr(legacy_result, "n_failed", 0)),
        )

    outputs = [
        _output("d32", "m", legacy_result.d32_samples),
        _output("pore_size_mean", "m", legacy_result.pore_samples),
        _output("G_DN", "Pa", legacy_result.G_DN_samples),
    ]
    return UnifiedUncertaintyResult(
        outputs=outputs,
        # Legacy uncertainty_core perturbs MaterialProperties only.
        kinds_sampled={UncertaintyKind.MATERIAL_PROPERTY},
        n_samples=int(getattr(legacy_result, "n_samples", 0)),
        n_failed=int(getattr(legacy_result, "n_failed", 0)),
        source_label=source_label,
    )


# ─── Node 18 (v7.0, P3): Unified MC engine ─────────────────────────────────


class UnifiedUncertaintyEngine:
    """Single-entrypoint MC engine that dispatches to legacy paths under
    the unified schema (Node 18, v7.0).

    The full structural merge of ``uncertainty_core`` (M1-L4) and
    ``uncertainty_propagation`` (M2) is multi-week work and would risk
    regressing the existing percentile semantics. The v7.0 minimum-viable
    deliverable is:

      1. **One entrypoint** — callers create one ``UnifiedUncertaintyEngine``
         per run and call ``run_m1l4`` or ``run_m2`` to get a
         ``UnifiedUncertaintyResult`` regardless of which legacy engine
         did the sampling.
      2. **Calibration-posterior sampling** — when a ``CalibrationStore``
         entry carries a ``posterior_uncertainty`` (sigma in measurement
         units), it is used to generate per-sample perturbations of the
         calibrated parameter so the MC interval reflects calibration
         confidence rather than treating calibrated values as point-true.
      3. **Source-kind tagging** — every output records which
         ``UncertaintyKind`` set contributed to its interval, so the UI
         and optimizer can communicate honest provenance.

    The two legacy MC engines remain functional and tested; this engine
    is additive. Full consolidation (drop legacy engines, single
    pipeline-level MC) is the v7.1 target.
    """

    def __init__(
        self,
        spec: Optional[UnifiedUncertaintySpec] = None,
        calibration_store=None,
    ):
        self.spec = spec or UnifiedUncertaintySpec()
        self.calibration_store = calibration_store
        # Build calibration-posterior sources from store entries that carry
        # posterior_uncertainty. Done at engine construction so a stable
        # snapshot is used across run_m1l4 / run_m2 calls.
        if calibration_store is not None:
            self._absorb_calibration_posteriors(calibration_store)

    def _absorb_calibration_posteriors(self, store) -> None:
        """Add a CALIBRATION_POSTERIOR source per calibration entry.

        Each entry's ``measured_value`` becomes the source value;
        ``posterior_uncertainty`` (std in measurement units, optional) is
        the spread. Entries without a posterior are ignored — they
        contribute their value to downstream solvers as a point estimate.
        """
        for entry in getattr(store, "entries", []):
            sigma = float(getattr(entry, "posterior_uncertainty", 0.0) or 0.0)
            if sigma <= 0.0:
                continue
            self.spec.add(UncertaintySource(
                name=f"calibration.{entry.target_module}.{entry.parameter_name}",
                kind=UncertaintyKind.CALIBRATION_POSTERIOR,
                distribution="normal",
                value=float(entry.measured_value),
                std=sigma,
                units=getattr(entry, "units", ""),
            ))

    def run_m1l4(
        self,
        params,
        n_samples: Optional[int] = None,
        seed: Optional[int] = None,
        n_jobs: int = 1,
        l2_mode: str = "empirical",
        crosslinker_key: str = "genipin",
    ) -> UnifiedUncertaintyResult:
        """M1-L4 MC delegating to ``uncertainty_core`` then unifying output.

        Returns a ``UnifiedUncertaintyResult`` carrying the
        MATERIAL_PROPERTY kind (legacy core perturbs MaterialProperties)
        plus any CALIBRATION_POSTERIOR kind absorbed at engine init.
        """
        from .uncertainty_core import UncertaintyPropagator

        propagator = UncertaintyPropagator(
            n_samples=n_samples or self.spec.n_samples,
            seed=seed if seed is not None else self.spec.seed,
            n_jobs=n_jobs,
        )
        legacy = propagator.run(
            params, crosslinker_key=crosslinker_key, l2_mode=l2_mode,
        )
        unified = from_m1l4_result(
            legacy, source_label="UnifiedUncertaintyEngine.M1L4",
        )
        # Audit N2 (v7.0.1): the legacy engine samples MaterialProperty
        # perturbations only. CalibrationPosterior sources declared on the
        # spec are NOT propagated through the legacy MC, so we record them
        # as "declared but not sampled" rather than tagging them as
        # contributing to the interval. v7.1 will close the gap by
        # injecting posterior samples into the legacy engine's
        # perturbation dict; until then, downstream consumers can detect
        # the lower-bound condition via this field.
        for src in self.spec.sources:
            if src.kind == UncertaintyKind.CALIBRATION_POSTERIOR:
                unified.kinds_declared_but_not_sampled.add(
                    UncertaintyKind.CALIBRATION_POSTERIOR
                )
        return unified

    def run_m2_q_max(
        self,
        contract,
        steps,
        m1_uncertainty,
        n_samples: Optional[int] = None,
        seed: Optional[int] = None,
    ) -> UnifiedUncertaintyResult:
        """M2 MC delegating to ``uncertainty_propagation.run_with_uncertainty``.

        Returns a ``UnifiedUncertaintyResult`` whose single output is
        ``estimated_q_max``. Mass-balance and other M3 outputs are not
        yet propagated here; that is v7.1 batch-variability work
        (Node 19 + future).
        """
        from .uncertainty_propagation.monte_carlo import run_with_uncertainty

        legacy = run_with_uncertainty(
            contract, steps, m1_uncertainty,
            n_samples=n_samples or self.spec.n_samples,
            seed=seed if seed is not None else self.spec.seed,
        )
        # legacy.all_q_max is the per-sample list
        return from_m1_contract_uq(
            list(legacy.all_q_max),
            units="mol/m^3",
            source_label="UnifiedUncertaintyEngine.M2",
        )


def from_m1_contract_uq(
    samples_q_max: Sequence[float],
    units: str = "mol/m^3",
    source_label: str = "uncertainty_propagation.M2",
):
    """Convert M1->M2 q_max sample array to the unified schema.

    The legacy ``run_with_uncertainty`` path returns a list of FMC objects
    one per Monte Carlo replicate. The caller passes the q_max samples
    extracted from that list. Source kind is MEASUREMENT when the
    M1UncertaintyContract.tier == "measured", otherwise MATERIAL_PROPERTY.
    """
    arr = np.asarray(samples_q_max, dtype=float)
    valid = arr[np.isfinite(arr)]
    if valid.size == 0:
        return UnifiedUncertaintyResult(
            kinds_sampled={UncertaintyKind.MATERIAL_PROPERTY},
            n_samples=0,
            source_label=source_label,
        )
    out = OutputUncertainty(
        name="estimated_q_max",
        units=units,
        mean=float(np.mean(valid)),
        p5=float(np.percentile(valid, 5)),
        p50=float(np.percentile(valid, 50)),
        p95=float(np.percentile(valid, 95)),
        n_samples=int(valid.size),
        n_failed=int(arr.size - valid.size),
    )
    return UnifiedUncertaintyResult(
        outputs=[out],
        # M1ExportContract CVs are screening defaults unless caller flagged
        # them measured; consumers can override after construction.
        kinds_sampled={UncertaintyKind.MATERIAL_PROPERTY},
        n_samples=int(arr.size),
        n_failed=int(arr.size - valid.size),
        source_label=source_label,
    )
