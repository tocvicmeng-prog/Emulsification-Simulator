"""Monte Carlo uncertainty propagation for simulation predictions."""

from __future__ import annotations

import copy
import logging
from dataclasses import dataclass
from typing import Optional

import numpy as np

from .datatypes import SimulationParameters, MaterialProperties, FullResult
from .pipeline.orchestrator import PipelineOrchestrator
from .properties.database import PropertyDatabase

logger = logging.getLogger(__name__)


@dataclass
class UncertaintyResult:
    """Uncertainty-quantified simulation results."""
    # Point estimates (median)
    d32_median: float        # [m]
    pore_median: float       # [m]
    G_DN_median: float       # [Pa]

    # Confidence intervals (5th-95th percentile)
    d32_ci: tuple[float, float]      # (low, high) [m]
    pore_ci: tuple[float, float]     # (low, high) [m]
    G_DN_ci: tuple[float, float]     # (low, high) [Pa]

    # All samples
    d32_samples: np.ndarray
    pore_samples: np.ndarray
    G_DN_samples: np.ndarray

    # Additional
    n_samples: int
    n_failed: int

    def summary(self) -> str:
        """Human-readable uncertainty summary."""
        lines = [
            "Uncertainty-Quantified Results:",
            f"  d32  = {self.d32_median*1e6:.2f} \u00b5m  "
            f"(90% CI: {self.d32_ci[0]*1e6:.2f} \u2013 {self.d32_ci[1]*1e6:.2f} \u00b5m)",
            f"  pore = {self.pore_median*1e9:.0f} nm  "
            f"(90% CI: {self.pore_ci[0]*1e9:.0f} \u2013 {self.pore_ci[1]*1e9:.0f} nm)",
            f"  G_DN = {self.G_DN_median/1000:.1f} kPa  "
            f"(90% CI: {self.G_DN_ci[0]/1000:.1f} \u2013 {self.G_DN_ci[1]/1000:.1f} kPa)",
            f"  ({self.n_samples} samples, {self.n_failed} failed)",
        ]
        return "\n".join(lines)


class UncertaintyPropagator:
    """Monte Carlo uncertainty propagation over material property uncertainties.

    Perturbation parameters are organised into three conceptually distinct
    uncertainty categories:

    Parameter uncertainty
        Stems from measurement or literature scatter in physical material
        properties that feed directly into the governing model equations
        (interfacial tension, viscosities, kinetic rate constants, modulus
        prefactors).  These parameters are propagated *pre-hoc* — they enter
        the solver before any physics is evaluated.

    Calibration uncertainty
        Reflects imprecision in the empirical coefficients that were fitted to
        experimental data.  Examples include the Alopaeus C3 viscous-correction
        constant, and the prefactor/exponent of the L2 pore-size correlation.
        These are also propagated pre-hoc where possible; the pore-model
        coefficients use a post-hoc approximation that is exact for the
        algebraic (empirical) L2 path (see note in ``run()``).

    Model-form uncertainty
        Accounts for structural assumptions embedded in the model that cannot
        be derived from first principles alone.  Examples are the bridge
        efficiency fraction for crosslinking and the IPN coupling coefficient.
        These are sampled from physically motivated uniform priors rather than
        being perturbed around a point estimate.
    """

    def __init__(self, n_samples: int = 20, seed: int = 42, n_jobs: int = 1):
        """Initialise an MC propagator.

        Parameters
        ----------
        n_samples : int
            Number of Monte Carlo samples.
        seed : int
            Master RNG seed for reproducibility.
        n_jobs : int
            Node 15 (v7.0, P5b): joblib parallelism. ``1`` (default) keeps
            the legacy serial loop. ``-1`` uses every available core.
            ``PropertyDatabase`` is in-memory and process-safe by
            construction (no disk cache, no shared mutable state), so
            joblib's ``"loky"`` backend is safe with default settings.
            Each worker re-imports the solver modules; first-call Numba JIT
            compile cost (Node 14) is paid per worker on cold start.
        """
        self.n_samples = n_samples
        self.rng = np.random.default_rng(seed)
        self.n_jobs = n_jobs

    def _generate_perturbations(self) -> list[dict]:
        """Generate LHS-like perturbation factors for uncertain parameters."""
        n = self.n_samples
        perturbations = []

        for i in range(n):
            p = {
                # ── Parameter uncertainty (propagated through model equations) ──
                'sigma_factor': np.exp(self.rng.uniform(-0.26, 0.26)),       # interfacial tension             +/-30%
                'mu_d_factor': np.exp(self.rng.uniform(-0.40, 0.40)),        # dispersed phase viscosity       +/-50%
                'k_xlink_factor': np.exp(self.rng.uniform(-0.34, 0.34)),     # crosslinking rate constant      +/-40%
                'G_prefactor_factor': np.exp(self.rng.uniform(-0.26, 0.26)), # agarose modulus prefactor       +/-30%

                # ── Calibration uncertainty (empirical model coefficients) ──
                'breakage_C3_factor': 10 ** self.rng.uniform(-0.30, 0.30),   # breakage viscous correction     +/-50%
                'pore_prefactor_factor': 10 ** self.rng.uniform(-0.13, 0.13), # L2 pore correlation prefactor  +/-30%
                'pore_exponent_offset': self.rng.uniform(-0.1, 0.1),          # L2 pore correlation exponent   on -0.7 exponent
                'confinement_alpha': self.rng.uniform(0.10, 0.25),            # confinement correction

                # ── Model-form uncertainty (structural assumptions) ──
                'f_bridge_factor': self.rng.uniform(0.25, 0.55),             # crosslinking bridge efficiency
                'eta_coupling_offset': self.rng.uniform(-0.30, 0.0),         # IPN coupling coefficient
            }
            perturbations.append(p)

        return perturbations

    def _apply_perturbation(self, props: MaterialProperties,
                             perturb: dict) -> MaterialProperties:
        """Apply a single perturbation to material properties."""
        p = copy.deepcopy(props)
        p.sigma *= perturb['sigma_factor']
        p.mu_d *= perturb['mu_d_factor']
        p.k_xlink_0 *= perturb['k_xlink_factor']
        p.G_agarose_prefactor *= perturb['G_prefactor_factor']
        p.f_bridge = perturb['f_bridge_factor']
        p.eta_coupling = perturb['eta_coupling_offset']
        return p

    def run(self, params: SimulationParameters,
            db: Optional[PropertyDatabase] = None,
            crosslinker_key: str = 'genipin',
            uv_intensity: float = 0.0,
            l2_mode: str = 'empirical') -> UncertaintyResult:
        """Run Monte Carlo uncertainty propagation.

        Parameters
        ----------
        params : SimulationParameters
        db : PropertyDatabase, optional
        crosslinker_key : str
            Key into CROSSLINKERS library to pass through to L3 solver.
        uv_intensity : float
            UV intensity [mW/cm2] for UV-initiated crosslinkers.
        l2_mode : str
            Gelation solver mode ('empirical', 'ch_2d', 'ch_1d').
        """
        db = db or PropertyDatabase()

        # Get base properties
        base_props = db.update_for_conditions(
            T_oil=params.formulation.T_oil,
            c_agarose=params.formulation.c_agarose,
            c_chitosan=params.formulation.c_chitosan,
            c_span80=params.formulation.c_span80,
        )

        perturbations = self._generate_perturbations()

        # Pre-compute perturbed prop snapshots so the worker function is a
        # pure mapping (i, props_i, params_i, perturb) -> sample triple. This
        # is the joblib boundary: keep RNG draws on the master process so the
        # results are bit-identical to the serial path (joblib does not
        # guarantee determinism if workers seed independently).
        prepared = []
        for i, perturb in enumerate(perturbations):
            props_i = self._apply_perturbation(base_props, perturb)
            # Apply L1 breakage_C3 perturbation here so workers see the
            # final perturbed props (legacy behaviour kept identical).
            props_i.breakage_C3 = props_i.breakage_C3 * perturb['breakage_C3_factor']
            params_i = copy.deepcopy(params)
            params_i.run_id = f"mc_{i:03d}"
            prepared.append((i, props_i, params_i, perturb))

        # Audit N7 (v7.0.1): joblib loky process-startup + Numba JIT
        # cold-compile dominate when n_samples is small. Auto-fallback to
        # serial when the workload is too small for parallelism to pay off.
        # Heuristic: at least 4 samples per worker. Caller can override
        # by passing n_jobs=1 explicitly (no fallback) or n_jobs=-1
        # (request all cores; only auto-fallbacks if n_samples < 4).
        effective_n_jobs = self.n_jobs
        if effective_n_jobs != 1:
            requested = effective_n_jobs
            if requested == -1:
                import os
                requested = os.cpu_count() or 1
            min_samples_for_parallel = 4 * abs(requested)
            if len(prepared) < min_samples_for_parallel:
                logger.info(
                    "MC n_jobs=%d but n_samples=%d (<%d); falling back to "
                    "serial — joblib startup + Numba JIT cold-compile would "
                    "dominate.", self.n_jobs, len(prepared),
                    min_samples_for_parallel,
                )
                effective_n_jobs = 1

        if effective_n_jobs == 1:
            # Serial path — preserves byte-for-byte legacy behaviour and
            # avoids joblib's process-fork overhead for tiny n_samples.
            outputs = [
                _mc_one_sample(item, crosslinker_key, uv_intensity, l2_mode)
                for item in prepared
            ]
        else:
            # Parallel path (Node 15). Loky backend uses cloudpickle to
            # transfer params; PropertyDatabase is rebuilt per worker via
            # solver imports.
            from joblib import Parallel, delayed
            outputs = Parallel(n_jobs=effective_n_jobs, backend="loky")(
                delayed(_mc_one_sample)(
                    item, crosslinker_key, uv_intensity, l2_mode,
                ) for item in prepared
            )

        d32_samples = []
        pore_samples = []
        G_DN_samples = []
        n_failed = 0
        for i, out in enumerate(outputs):
            if out is None:
                n_failed += 1
                continue
            d32_i, pore_i, G_DN_i = out
            d32_samples.append(d32_i)
            pore_samples.append(pore_i)
            G_DN_samples.append(G_DN_i)

        d32_arr = np.array(d32_samples)
        pore_arr = np.array(pore_samples)
        G_DN_arr = np.array(G_DN_samples)

        if len(d32_arr) == 0:
            raise RuntimeError("All MC samples failed")

        return UncertaintyResult(
            d32_median=float(np.median(d32_arr)),
            pore_median=float(np.median(pore_arr)),
            G_DN_median=float(np.median(G_DN_arr)),
            d32_ci=(float(np.percentile(d32_arr, 5)), float(np.percentile(d32_arr, 95))),
            pore_ci=(float(np.percentile(pore_arr, 5)), float(np.percentile(pore_arr, 95))),
            G_DN_ci=(float(np.percentile(G_DN_arr, 5)), float(np.percentile(G_DN_arr, 95))),
            d32_samples=d32_arr,
            pore_samples=pore_arr,
            G_DN_samples=G_DN_arr,
            n_samples=len(d32_arr),
            n_failed=n_failed,
        )


# ─── Node 15 (v7.0, P5b): module-level worker for joblib parallelism ──────


def _mc_one_sample(item, crosslinker_key, uv_intensity, l2_mode):
    """Run one full pipeline sample for the MC engine.

    Module-level (not a method) so joblib's loky backend can pickle it
    without dragging the UncertaintyPropagator instance across the process
    boundary. Returns (d32, pore, G_DN) on success, None on failure.

    The perturbations are pre-applied on the master process; this function
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
            * perturb['pore_prefactor_factor']
            * np.exp(perturb['pore_exponent_offset'] * 0.5)
        )
        pore_conf_max = perturb['confinement_alpha'] * 2.0 * R
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
