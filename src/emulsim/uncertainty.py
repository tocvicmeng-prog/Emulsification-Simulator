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
    """Monte Carlo uncertainty propagation over material property uncertainties."""

    def __init__(self, n_samples: int = 20, seed: int = 42):
        self.n_samples = n_samples
        self.rng = np.random.default_rng(seed)

    def _generate_perturbations(self) -> list[dict]:
        """Generate LHS-like perturbation factors for uncertain parameters."""
        n = self.n_samples
        perturbations = []

        for i in range(n):
            p = {
                # Multiplicative factors (lognormal-like via uniform on log scale)
                'sigma_factor': np.exp(self.rng.uniform(-0.26, 0.26)),      # +/-30%
                'mu_d_factor': np.exp(self.rng.uniform(-0.40, 0.40)),       # +/-50%
                'k_xlink_0_factor': np.exp(self.rng.uniform(-0.34, 0.34)), # +/-40%
                'G_prefactor_factor': np.exp(self.rng.uniform(-0.26, 0.26)), # +/-30%
                # Absolute values (uniform)
                'f_bridge': self.rng.uniform(0.25, 0.55),
                'eta_coupling': self.rng.uniform(-0.30, 0.0),
            }
            perturbations.append(p)

        return perturbations

    def _apply_perturbation(self, props: MaterialProperties,
                             perturb: dict) -> MaterialProperties:
        """Apply a single perturbation to material properties."""
        p = copy.deepcopy(props)
        p.sigma *= perturb['sigma_factor']
        p.mu_d *= perturb['mu_d_factor']
        p.k_xlink_0 *= perturb['k_xlink_0_factor']
        p.G_agarose_prefactor *= perturb['G_prefactor_factor']
        p.f_bridge = perturb['f_bridge']
        p.eta_coupling = perturb['eta_coupling']
        return p

    def run(self, params: SimulationParameters,
            db: Optional[PropertyDatabase] = None,
            crosslinker_key: str = 'genipin',
            uv_intensity: float = 0.0) -> UncertaintyResult:
        """Run Monte Carlo uncertainty propagation.

        Parameters
        ----------
        params : SimulationParameters
        db : PropertyDatabase, optional
        crosslinker_key : str
            Key into CROSSLINKERS library to pass through to L3 solver.
        uv_intensity : float
            UV intensity [mW/cm2] for UV-initiated crosslinkers.
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

        d32_samples = []
        pore_samples = []
        G_DN_samples = []
        n_failed = 0

        for i, perturb in enumerate(perturbations):
            props_i = self._apply_perturbation(base_props, perturb)

            # Override the database properties for this run
            params_i = copy.deepcopy(params)
            params_i.run_id = f"mc_{i:03d}"

            try:
                # Run pipeline with perturbed properties
                # Bypass the orchestrator's update_for_conditions
                # and inject perturbed props directly
                from .level1_emulsification.solver import PBESolver
                from .level2_gelation.solver import solve_gelation
                from .level3_crosslinking.solver import solve_crosslinking
                from .level4_mechanical.solver import solve_mechanical

                phi_d = params_i.formulation.phi_d

                # L1
                pbe = PBESolver(
                    n_bins=params_i.solver.l1_n_bins,
                    d_min=params_i.solver.l1_d_min,
                    d_max=params_i.solver.l1_d_max,
                )
                emul = pbe.solve(params_i, props_i, phi_d=phi_d)

                # L2
                R = emul.d50 / 2.0
                gel = solve_gelation(params_i, props_i, R_droplet=R, mode='empirical')

                # L3 -- Fix 8: pass crosslinker_key and uv_intensity through
                xl = solve_crosslinking(params_i, props_i, R_droplet=R,
                                        porosity=gel.porosity,
                                        crosslinker_key=crosslinker_key,
                                        uv_intensity=uv_intensity)

                # L4
                mech = solve_mechanical(params_i, props_i, gel, xl)

                d32_samples.append(emul.d32)
                pore_samples.append(gel.pore_size_mean)
                G_DN_samples.append(mech.G_DN)

            except Exception as e:
                logger.debug("MC sample %d failed: %s", i, e)
                n_failed += 1

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
