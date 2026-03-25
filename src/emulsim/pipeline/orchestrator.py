"""Pipeline orchestrator for sequential L1→L2→L3→L4 execution."""

from __future__ import annotations

import json
import logging
import time
from datetime import datetime
from pathlib import Path
from typing import Optional

from ..datatypes import (
    EmulsificationResult,
    GelationResult,
    MaterialProperties,
    SimulationParameters,
)
from ..properties.database import PropertyDatabase
from ..level1_emulsification.solver import PBESolver
from ..level2_gelation.solver import CahnHilliardSolver

logger = logging.getLogger(__name__)


class PipelineOrchestrator:
    """Orchestrates the sequential simulation pipeline.

    Currently implements Level 1 (emulsification). Levels 2-4 will be
    added in subsequent phases.
    """

    def __init__(self, db: Optional[PropertyDatabase] = None,
                 output_dir: Optional[Path] = None):
        self.db = db or PropertyDatabase()
        self.output_dir = Path(output_dir) if output_dir else Path("output")
        self.output_dir.mkdir(parents=True, exist_ok=True)

    def run_single(self, params: SimulationParameters,
                   phi_d: float = 0.05) -> dict:
        """Execute the pipeline for a single parameter set.

        Returns a dict with results from each level completed.
        """
        run_id = params.run_id or datetime.now().strftime("%Y%m%d_%H%M%S")
        run_dir = self.output_dir / f"run_{run_id}"
        run_dir.mkdir(parents=True, exist_ok=True)

        logger.info("Starting pipeline run: %s", run_id)
        results = {"run_id": run_id, "params": params}

        # Update material properties for actual conditions
        props = self.db.update_for_conditions(
            T_oil=params.formulation.T_oil,
            c_agarose=params.formulation.c_agarose,
            c_chitosan=params.formulation.c_chitosan,
            c_span80=params.formulation.c_span80,
        )

        # ── Level 1: Emulsification ──────────────────────────────────────
        logger.info("Level 1: Emulsification (RPM=%.0f)", params.emulsification.rpm)
        t0 = time.perf_counter()

        solver = PBESolver(
            n_bins=params.solver.l1_n_bins,
            d_min=params.solver.l1_d_min,
            d_max=params.solver.l1_d_max,
        )
        emul_result = solver.solve(params, props, phi_d=phi_d)
        dt1 = time.perf_counter() - t0

        results["emulsification"] = emul_result
        results["level1_time"] = dt1
        logger.info(
            "Level 1 complete: d32=%.2f µm, span=%.2f (%.1f s)",
            emul_result.d32 * 1e6, emul_result.span, dt1,
        )

        # ── Level 2: Gelation & Pore Formation ────────────────────────────
        R_droplet = emul_result.d32 / 2.0
        logger.info(
            "Level 2: Gelation (R=%.2f µm, cooling=%.3f K/s)",
            R_droplet * 1e6, params.formulation.cooling_rate,
        )
        t0 = time.perf_counter()

        ch_solver = CahnHilliardSolver(
            N_r=params.solver.l2_n_r,
            dt_initial=params.solver.l2_dt_initial,
            dt_max=params.solver.l2_dt_max,
            arrest_exponent=params.solver.l2_arrest_exponent,
        )
        gel_result = ch_solver.solve(params, props, R_droplet=R_droplet)
        dt2 = time.perf_counter() - t0

        results["gelation"] = gel_result
        results["level2_time"] = dt2
        logger.info(
            "Level 2 complete: pore=%.1f nm, porosity=%.2f, alpha=%.3f (%.1f s)",
            gel_result.pore_size_mean * 1e9, gel_result.porosity,
            gel_result.alpha_final, dt2,
        )

        # Save summary
        summary = {
            "run_id": run_id,
            "level1": {
                "d32_um": emul_result.d32 * 1e6,
                "d43_um": emul_result.d43 * 1e6,
                "d10_um": emul_result.d10 * 1e6,
                "d50_um": emul_result.d50 * 1e6,
                "d90_um": emul_result.d90 * 1e6,
                "span": emul_result.span,
                "converged": bool(emul_result.converged),
                "elapsed_s": dt1,
            },
            "level2": {
                "pore_size_mean_nm": gel_result.pore_size_mean * 1e9,
                "pore_size_std_nm": gel_result.pore_size_std * 1e9,
                "porosity": gel_result.porosity,
                "alpha_final": gel_result.alpha_final,
                "char_wavelength_nm": gel_result.char_wavelength * 1e9,
                "elapsed_s": dt2,
            },
        }
        with open(run_dir / "summary.json", "w") as f:
            json.dump(summary, f, indent=2)

        return results

    def run_rpm_sweep(self, rpms: list[float],
                      base_params: Optional[SimulationParameters] = None,
                      phi_d: float = 0.05) -> list[dict]:
        """Run Level 1 for a list of RPM values.

        Returns list of summary dicts with RPM and d32.
        """
        base_params = base_params or SimulationParameters()
        results = []

        props = self.db.update_for_conditions(
            T_oil=base_params.formulation.T_oil,
            c_agarose=base_params.formulation.c_agarose,
            c_chitosan=base_params.formulation.c_chitosan,
            c_span80=base_params.formulation.c_span80,
        )

        solver = PBESolver(
            n_bins=base_params.solver.l1_n_bins,
            d_min=base_params.solver.l1_d_min,
            d_max=base_params.solver.l1_d_max,
        )

        for rpm in rpms:
            import copy
            params = copy.deepcopy(base_params)
            params.emulsification.rpm = rpm

            t0 = time.perf_counter()
            emul = solver.solve(params, props, phi_d=phi_d)
            dt = time.perf_counter() - t0

            results.append({
                "rpm": rpm,
                "d32_um": emul.d32 * 1e6,
                "d50_um": emul.d50 * 1e6,
                "span": emul.span,
                "converged": emul.converged,
                "elapsed_s": dt,
            })
            logger.info("RPM=%d  d32=%.2f µm  (%.1f s)", rpm, emul.d32 * 1e6, dt)

        return results
