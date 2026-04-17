"""Unified Monte Carlo uncertainty quantification for EmulSim.

Node 30 (v7.1) consolidates the two legacy UQ engines
(``uncertainty_core`` — M1-L4 MaterialProperties MC, and
``uncertainty_propagation`` — M2 M1ExportContract MC) into a single MC
implementation that lives here. The public entry points are:

    from emulsim.uncertainty_unified import (
        UnifiedUncertaintyEngine,
        UnifiedUncertaintySpec,
        UnifiedUncertaintyResult,
        UncertaintySource,
        UncertaintyKind,
    )

The ``UnifiedUncertaintyEngine.run_m1l4`` method drives the full M1-L4
pipeline with:

1. **Default MaterialProperties perturbations** — the 10-factor
   perturbation scheme from v6.x/v7.0 (interfacial tension, viscosity,
   kinetic prefactors, L2 pore coefficients, model-form bridge + IPN
   coupling). RNG call order is preserved byte-for-byte from the
   v7.0.1 ``uncertainty_core.UncertaintyPropagator._generate_perturbations``
   so seed-identical output matches when no posterior sources are declared.

2. **Calibration-posterior injection (Audit N2 closure)** — when a
   ``CalibrationStore`` entry carries ``posterior_uncertainty > 0`` the
   engine samples per-MC-sample from ``Normal(measured_value, posterior)``
   and writes the draw onto the appropriate parameter object (KernelConfig
   for target_module=="L1", MaterialProperties otherwise). This closes
   the v7.0 gap where posteriors were absorbed into the spec but not
   actually propagated through the MC.

3. **Source-kind provenance** — every output records which
   ``UncertaintyKind`` set contributed to its interval. Posteriors
   declared on the spec but whose parameter_name does not dispatch to a
   known attribute are logged at DEBUG and surfaced via
   ``kinds_declared_but_not_sampled`` so consumers do not over-attribute.

Parallelism (Node 15) and the small-n-jobs auto-serial-fallback (Node 27
/ Audit N7) are preserved. PropertyDatabase is in-memory and process-safe
by construction.
"""

from __future__ import annotations

import copy
import logging
from dataclasses import dataclass, field
from enum import Enum
from typing import Optional

import numpy as np

from .datatypes import MaterialProperties, SimulationParameters
from .properties.database import PropertyDatabase

logger = logging.getLogger(__name__)


# ─── Schema ────────────────────────────────────────────────────────────────


class UncertaintyKind(Enum):
    """Classification of an uncertainty source (consensus plan §6)."""

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
    """One named, classified, sampleable uncertainty input."""

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
            return rng.lognormal(np.log(max(self.value, 1e-300)), self.std, size=n)
        raise ValueError(f"Unsupported distribution {self.distribution!r}")


@dataclass
class UnifiedUncertaintySpec:
    """Input spec for the unified MC engine."""

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
    raw_samples: Optional[np.ndarray] = None
    """Optional raw sample array. Populated by the engine for callers that
    need to run bespoke analyses (e.g. parallelism bit-identicality
    tests). May be None for downstream consumers that rebuild an
    ``OutputUncertainty`` from summary stats alone."""

    @property
    def ci_width_relative(self) -> float:
        """(p95 - p5) / |p50|; useful as a single uncertainty scalar."""
        if abs(self.p50) < 1e-300:
            return float("nan")
        return (self.p95 - self.p5) / abs(self.p50)


@dataclass
class UnifiedUncertaintyResult:
    """Schema-uniform UQ output for downstream consumers.

    ``kinds_sampled`` records ONLY the kinds that actually contributed to
    the percentile interval. Kinds declared on the spec but not
    dispatched by the engine (e.g. a posterior whose parameter_name does
    not match any attribute on MaterialProperties/KernelConfig) are
    surfaced via ``kinds_declared_but_not_sampled`` so downstream
    consumers can detect the lower-bound condition.
    """

    outputs: list[OutputUncertainty] = field(default_factory=list)
    kinds_sampled: set[UncertaintyKind] = field(default_factory=set)
    n_samples: int = 0
    n_failed: int = 0
    source_label: str = ""
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
            unsamp = ", ".join(
                sorted(k.value for k in self.kinds_declared_but_not_sampled)
            )
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


# ─── Default perturbation generator (byte-compat with v7.0.1) ──────────────


# Dispatch target for posterior injection: the calibration source name is
# "calibration.{target_module}.{parameter_name}". L1 posteriors target the
# KernelConfig (lives at params.emulsification.kernels); all other modules
# target MaterialProperties.
_KERNEL_MODULES = {"L1"}
_PROPS_MODULES = {"L2", "L3", "L4", "M2", "M3"}


def _generate_default_perturbations(rng: np.random.Generator) -> dict:
    """Generate the legacy per-sample MaterialProperties perturbation dict.

    RNG call order is byte-identical to the v7.0.1
    ``UncertaintyPropagator._generate_perturbations`` so seed-identical
    output matches when no posterior sources are declared.
    """
    return {
        # Parameter uncertainty (material properties)
        "sigma_factor": float(np.exp(rng.uniform(-0.26, 0.26))),        # IFT ±30%
        "mu_d_factor": float(np.exp(rng.uniform(-0.40, 0.40))),         # mu_d ±50%
        "k_xlink_factor": float(np.exp(rng.uniform(-0.34, 0.34))),      # k_xlink ±40%
        "G_prefactor_factor": float(np.exp(rng.uniform(-0.26, 0.26))),  # G_prefactor ±30%
        # Calibration uncertainty (empirical coefficients)
        "breakage_C3_factor": float(10 ** rng.uniform(-0.30, 0.30)),    # ±50%
        "pore_prefactor_factor": float(10 ** rng.uniform(-0.13, 0.13)), # ±30%
        "pore_exponent_offset": float(rng.uniform(-0.1, 0.1)),
        "confinement_alpha": float(rng.uniform(0.10, 0.25)),
        # Model-form uncertainty (structural assumptions)
        "f_bridge_factor": float(rng.uniform(0.25, 0.55)),
        "eta_coupling_offset": float(rng.uniform(-0.30, 0.0)),
    }


def _apply_default_perturbations(
    base_props: MaterialProperties, p: dict
) -> MaterialProperties:
    """Apply the default MaterialProperties perturbation factors."""
    q = copy.deepcopy(base_props)
    q.sigma *= p["sigma_factor"]
    q.mu_d *= p["mu_d_factor"]
    q.k_xlink_0 *= p["k_xlink_factor"]
    q.G_agarose_prefactor *= p["G_prefactor_factor"]
    q.f_bridge = p["f_bridge_factor"]
    q.eta_coupling = p["eta_coupling_offset"]
    q.breakage_C3 = q.breakage_C3 * p["breakage_C3_factor"]
    return q


def _apply_posterior_samples(
    props: MaterialProperties,
    params: SimulationParameters,
    posterior_samples: dict,
) -> set:
    """Apply sampled posterior values to the appropriate parameter object.

    Source names follow ``calibration.{target_module}.{parameter_name}``.
    L1 -> KernelConfig at ``params.emulsification.kernels``; L2-L4/M2-M3
    -> MaterialProperties. Malformed names, unknown modules, and unknown
    attributes are logged at DEBUG and skipped.

    Returns the set of UncertaintyKind values that actually reached a
    setattr (i.e. contributed to this sample's perturbation). When at
    least one posterior dispatch succeeded, the returned set contains
    ``CALIBRATION_POSTERIOR``; otherwise it is empty.
    """
    sampled_kinds: set[UncertaintyKind] = set()

    for name, value in posterior_samples.items():
        parts = name.split(".")
        if len(parts) != 3 or parts[0] != "calibration":
            logger.debug("Skipping malformed posterior source name %r", name)
            continue
        _, module, param = parts

        target: object | None = None
        if module in _KERNEL_MODULES:
            # KernelConfig is optional on EmulsificationParameters; create
            # a defaults instance on first posterior touch so the value
            # actually reaches the L1 solver rather than being swallowed
            # by the default-kernels fallback in the orchestrator.
            kernels = getattr(params.emulsification, "kernels", None)
            if kernels is None:
                from .datatypes import KernelConfig
                kernels = KernelConfig()
                params.emulsification.kernels = kernels
            target = kernels
        elif module in _PROPS_MODULES:
            target = props

        if target is None:
            logger.debug(
                "Posterior %r: no dispatch target for module %r", name, module
            )
            continue
        if not hasattr(target, param):
            logger.debug(
                "Posterior %r: target %s has no attribute %r",
                name, type(target).__name__, param,
            )
            continue

        setattr(target, param, value)
        sampled_kinds.add(UncertaintyKind.CALIBRATION_POSTERIOR)

    return sampled_kinds


# ─── Parallel worker (module-level so joblib loky can pickle it) ──────────


def _mc_one_sample(item, crosslinker_key, uv_intensity, l2_mode):
    """Run one full pipeline sample for the unified MC engine.

    Module-level (not a method) so joblib's loky backend can pickle it
    without dragging the engine instance across the process boundary.
    Returns (d32, pore, G_DN) on success, None on failure.

    Perturbations are pre-applied on the master process; this function
    just runs the deterministic pipeline given the perturbed props/params.
    """
    i, props_i, params_i, perturb = item
    try:
        # Local imports keep the worker bootstrap minimal — Numba JIT
        # compiles lazily on first call per process.
        from .level1_emulsification.solver import PBESolver
        from .level2_gelation.solver import solve_gelation
        from .level3_crosslinking.solver import solve_crosslinking
        from .level4_mechanical.solver import solve_mechanical

        pbe = PBESolver(
            n_bins=params_i.solver.l1_n_bins,
            d_min=params_i.solver.l1_d_min,
            d_max=params_i.solver.l1_d_max,
        )
        emul = pbe.solve(params_i, props_i)

        R = emul.d50 / 2.0
        gel = solve_gelation(params_i, props_i, R_droplet=R, mode=l2_mode)

        # Apply pore-model perturbations post-hoc (valid for empirical L2 only).
        pore_raw = gel.pore_size_mean
        pore_perturbed = (
            pore_raw
            * perturb["pore_prefactor_factor"]
            * np.exp(perturb["pore_exponent_offset"] * 0.5)
        )
        pore_conf_max = perturb["confinement_alpha"] * 2.0 * R
        if pore_conf_max > 0:
            pore_perturbed = min(pore_perturbed, pore_conf_max)
        gel.pore_size_mean = float(pore_perturbed)

        xl = solve_crosslinking(
            params_i, props_i, R_droplet=R, porosity=gel.porosity,
            crosslinker_key=crosslinker_key, uv_intensity=uv_intensity,
        )
        mech = solve_mechanical(params_i, props_i, gel, xl)
        return (emul.d32, gel.pore_size_mean, mech.G_DN)
    except Exception as exc:
        logger.debug("MC sample %d failed: %s", i, exc)
        return None


def _resolve_effective_n_jobs(n_jobs: int, n_prepared: int) -> int:
    """Audit N7 auto-serial-fallback for small workloads.

    joblib loky process-startup + Numba JIT cold-compile dominate when
    n_samples is small. Below 4 samples per worker, fall back to serial
    so the caller doesn't pay parallelism overhead for nothing.
    """
    if n_jobs == 1:
        return 1
    requested = n_jobs
    if requested == -1:
        import os
        requested = os.cpu_count() or 1
    min_samples_for_parallel = 4 * abs(requested)
    if n_prepared < min_samples_for_parallel:
        logger.info(
            "MC n_jobs=%d but n_samples=%d (<%d); falling back to serial — "
            "joblib startup + Numba JIT cold-compile would dominate.",
            n_jobs, n_prepared, min_samples_for_parallel,
        )
        return 1
    return n_jobs


def _build_output(
    name: str, unit: str, samples: list, n_failed: int
) -> OutputUncertainty:
    """Build an OutputUncertainty from a sample list, carrying raw array."""
    arr = np.asarray(samples, dtype=float)
    valid = arr[np.isfinite(arr)]
    if valid.size == 0:
        return OutputUncertainty(
            name=name, units=unit, mean=float("nan"),
            p5=float("nan"), p50=float("nan"), p95=float("nan"),
            n_samples=0, n_failed=int(n_failed),
            raw_samples=arr,
        )
    return OutputUncertainty(
        name=name, units=unit,
        mean=float(np.mean(valid)),
        p5=float(np.percentile(valid, 5)),
        p50=float(np.percentile(valid, 50)),
        p95=float(np.percentile(valid, 95)),
        n_samples=int(valid.size),
        n_failed=int(n_failed),
        raw_samples=arr,
    )


# ─── The unified engine ───────────────────────────────────────────────────


class UnifiedUncertaintyEngine:
    """Single-entrypoint MC engine for EmulSim uncertainty quantification.

    Call sites construct the engine with an optional spec + calibration
    store and invoke ``run_m1l4`` to execute the Monte Carlo. The engine
    is additive over the default MaterialProperties perturbations: any
    ``CalibrationStore`` entry with ``posterior_uncertainty > 0`` becomes
    an additional per-sample posterior draw that injects into the
    relevant parameter object at pipeline run time.
    """

    def __init__(
        self,
        spec: Optional[UnifiedUncertaintySpec] = None,
        calibration_store=None,
    ):
        self.spec = spec or UnifiedUncertaintySpec()
        self.calibration_store = calibration_store
        if calibration_store is not None:
            self._absorb_calibration_posteriors(calibration_store)

    def _absorb_calibration_posteriors(self, store) -> None:
        """Add a CALIBRATION_POSTERIOR source per calibration entry.

        Each entry's ``measured_value`` becomes the source value;
        ``posterior_uncertainty`` (std in measurement units) is the
        spread. Entries without a posterior are ignored — they contribute
        their value to downstream solvers as a point estimate via
        ``CalibrationStore.apply_to_model_params``.
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
        params: SimulationParameters,
        n_samples: Optional[int] = None,
        seed: Optional[int] = None,
        n_jobs: int = 1,
        l2_mode: str = "empirical",
        crosslinker_key: str = "genipin",
        uv_intensity: float = 0.0,
        db: Optional[PropertyDatabase] = None,
    ) -> UnifiedUncertaintyResult:
        """Run the full M1-L4 Monte Carlo and return a unified result.

        Parameters
        ----------
        params : SimulationParameters
            Base parameters. Deep-copied per sample so posterior
            injection into ``params.emulsification.kernels`` does not
            mutate the caller's object.
        n_samples : int, optional
            Number of MC samples. Defaults to ``self.spec.n_samples``.
        seed : int, optional
            Master RNG seed. Defaults to ``self.spec.seed``.
        n_jobs : int, default 1
            joblib workers. ``1`` keeps the serial loop and the
            byte-for-byte RNG determinism guarantee. ``-1`` uses every
            core. Auto-falls back to serial when ``n_samples < 4 * |n_jobs|``
            (Audit N7).
        l2_mode : str, default "empirical"
            L2 gelation solver mode ("empirical" | "ch_2d" | "ch_1d").
        crosslinker_key : str, default "genipin"
            Key into the crosslinker library passed to L3.
        uv_intensity : float, default 0.0
            UV intensity [mW/cm2] for UV-initiated crosslinkers.
        db : PropertyDatabase, optional
            Pre-built property database. Defaults to a fresh instance.
        """
        N = int(n_samples) if n_samples is not None else int(self.spec.n_samples)
        s = int(seed) if seed is not None else int(self.spec.seed)
        rng = np.random.default_rng(s)

        db = db or PropertyDatabase()
        base_props = db.update_for_conditions(
            T_oil=params.formulation.T_oil,
            c_agarose=params.formulation.c_agarose,
            c_chitosan=params.formulation.c_chitosan,
            c_span80=params.formulation.c_span80,
        )

        posterior_sources = [
            src for src in self.spec.sources
            if src.kind == UncertaintyKind.CALIBRATION_POSTERIOR
        ]

        prepared = []
        kinds_sampled_tracker: set[UncertaintyKind] = set()

        for i in range(N):
            # Byte-compat: default perturbations first, posteriors after.
            p = _generate_default_perturbations(rng)

            posterior_samples: dict = {}
            for src in posterior_sources:
                posterior_samples[src.name] = float(src.sample(rng, 1)[0])

            props_i = _apply_default_perturbations(base_props, p)
            params_i = copy.deepcopy(params)
            params_i.run_id = f"mc_{i:03d}"

            if posterior_samples:
                kinds_sampled_tracker.update(
                    _apply_posterior_samples(props_i, params_i, posterior_samples)
                )

            prepared.append((i, props_i, params_i, p))

        effective_n_jobs = _resolve_effective_n_jobs(n_jobs, len(prepared))
        if effective_n_jobs == 1:
            outputs = [
                _mc_one_sample(item, crosslinker_key, uv_intensity, l2_mode)
                for item in prepared
            ]
        else:
            from joblib import Parallel, delayed
            outputs = Parallel(n_jobs=effective_n_jobs, backend="loky")(
                delayed(_mc_one_sample)(
                    item, crosslinker_key, uv_intensity, l2_mode,
                ) for item in prepared
            )

        d32_samples: list[float] = []
        pore_samples: list[float] = []
        G_DN_samples: list[float] = []
        n_failed = 0
        for out in outputs:
            if out is None:
                n_failed += 1
                continue
            d32_i, pore_i, G_DN_i = out
            d32_samples.append(d32_i)
            pore_samples.append(pore_i)
            G_DN_samples.append(G_DN_i)

        if not d32_samples:
            raise RuntimeError("All MC samples failed")

        kinds_sampled: set[UncertaintyKind] = {UncertaintyKind.MATERIAL_PROPERTY}
        kinds_sampled.update(kinds_sampled_tracker)

        kinds_not_sampled: set[UncertaintyKind] = set()
        if posterior_sources and UncertaintyKind.CALIBRATION_POSTERIOR not in kinds_sampled:
            # Posteriors were declared but none of them dispatched to a known
            # attribute. Surface as lower-bound rather than silently dropping.
            kinds_not_sampled.add(UncertaintyKind.CALIBRATION_POSTERIOR)

        outputs_list = [
            _build_output("d32", "m", d32_samples, n_failed),
            _build_output("pore_size_mean", "m", pore_samples, n_failed),
            _build_output("G_DN", "Pa", G_DN_samples, n_failed),
        ]

        return UnifiedUncertaintyResult(
            outputs=outputs_list,
            kinds_sampled=kinds_sampled,
            n_samples=len(d32_samples),
            n_failed=n_failed,
            source_label="UnifiedUncertaintyEngine.M1L4",
            kinds_declared_but_not_sampled=kinds_not_sampled,
        )
