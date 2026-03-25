"""Pipeline orchestrator for sequential L1→L2→L3→L4 execution."""

from __future__ import annotations

import copy
import json
import logging
import time
from datetime import datetime
from pathlib import Path
from typing import Optional

from ..datatypes import (
    FullResult,
    MaterialProperties,
    SimulationParameters,
)
from ..properties.database import PropertyDatabase
from ..level1_emulsification.solver import PBESolver
from ..level2_gelation.solver import CahnHilliard2DSolver
from ..level3_crosslinking.solver import solve_crosslinking
from ..level4_mechanical.solver import solve_mechanical

logger = logging.getLogger(__name__)


class PipelineOrchestrator:
    """Orchestrates the full L1→L2→L3→L4 simulation pipeline."""

    def __init__(self, db: Optional[PropertyDatabase] = None,
                 output_dir: Optional[Path] = None):
        self.db = db or PropertyDatabase()
        self.output_dir = Path(output_dir) if output_dir else Path("output")
        self.output_dir.mkdir(parents=True, exist_ok=True)

    def run_single(self, params: SimulationParameters,
                   phi_d: float = None) -> FullResult:
        """Execute the full pipeline for a single parameter set."""
        if phi_d is None:
            phi_d = params.formulation.phi_d
        run_id = params.run_id or datetime.now().strftime("%Y%m%d_%H%M%S")
        run_dir = self.output_dir / f"run_{run_id}"
        run_dir.mkdir(parents=True, exist_ok=True)

        logger.info("Starting pipeline run: %s", run_id)
        timings = {}

        # Material properties interpolated to conditions
        props = self.db.update_for_conditions(
            T_oil=params.formulation.T_oil,
            c_agarose=params.formulation.c_agarose,
            c_chitosan=params.formulation.c_chitosan,
            c_span80=params.formulation.c_span80,
        )

        # ── Level 1: Emulsification ──────────────────────────────────────
        logger.info("L1: Emulsification (RPM=%.0f)", params.emulsification.rpm)
        t0 = time.perf_counter()
        pbe = PBESolver(
            n_bins=params.solver.l1_n_bins,
            d_min=params.solver.l1_d_min,
            d_max=params.solver.l1_d_max,
        )
        emul_result = pbe.solve(params, props, phi_d=phi_d)
        timings["L1"] = time.perf_counter() - t0
        logger.info("L1 done: d32=%.2f µm, span=%.2f (%.1fs)",
                     emul_result.d32 * 1e6, emul_result.span, timings["L1"])

        # ── Level 2: Gelation & Pore Formation ───────────────────────────
        # Use d50 (volume median) as the representative droplet size for L2.
        # d32 (Sauter mean) overweights large droplets; d50 better represents the
        # population center.  The full size distribution (d10–d90) is recorded in
        # the summary so downstream users can estimate pore-size spread via the
        # Cahn-Hilliard scaling λ* ~ sqrt(R).
        R_droplet = emul_result.d50 / 2.0
        logger.info("L2: Gelation (R=%.2f µm, d50-based)", R_droplet * 1e6)
        t0 = time.perf_counter()
        ch = CahnHilliard2DSolver(
            N_grid=params.solver.l2_n_grid,
            dt_initial=params.solver.l2_dt_initial,
            dt_max=params.solver.l2_dt_max,
            arrest_exponent=params.solver.l2_arrest_exponent,
        )
        gel_result = ch.solve(params, props, R_droplet=R_droplet)
        timings["L2"] = time.perf_counter() - t0
        logger.info("L2 done: pore=%.1f nm, porosity=%.2f (%.1fs)",
                     gel_result.pore_size_mean * 1e9, gel_result.porosity, timings["L2"])

        # ── Level 3: Crosslinking Kinetics ───────────────────────────────
        logger.info("L3: Crosslinking (genipin=%.1f mM, T=%.0f°C, t=%.0fh)",
                     params.formulation.c_genipin,
                     params.formulation.T_crosslink - 273.15,
                     params.formulation.t_crosslink / 3600)
        t0 = time.perf_counter()
        xlink_result = solve_crosslinking(
            params, props,
            R_droplet=R_droplet,
            porosity=gel_result.porosity,
        )
        timings["L3"] = time.perf_counter() - t0
        logger.info("L3 done: p=%.2f, G_chit=%.0f Pa, xi=%.1f nm (%.3fs)",
                     xlink_result.p_final, xlink_result.G_chitosan_final,
                     xlink_result.xi_final * 1e9, timings["L3"])

        # ── Level 4: Mechanical Properties ───────────────────────────────
        logger.info("L4: Mechanical property prediction")
        t0 = time.perf_counter()
        mech_result = solve_mechanical(params, props, gel_result, xlink_result)
        timings["L4"] = time.perf_counter() - t0
        logger.info("L4 done: G_DN=%.0f Pa, E*=%.0f Pa (%.4fs)",
                     mech_result.G_DN, mech_result.E_star, timings["L4"])

        total = sum(timings.values())
        logger.info("Pipeline complete: %.1fs total", total)

        # Save summary
        summary = {
            "run_id": run_id,
            "timings": timings,
            "total_time_s": total,
            "level1": {
                "d32_um": emul_result.d32 * 1e6,
                "d50_um": emul_result.d50 * 1e6,
                "span": emul_result.span,
                "converged": bool(emul_result.converged),
                "size_spread": {
                    "R_d10_um": emul_result.d10 * 0.5e6,
                    "R_d50_um": emul_result.d50 * 0.5e6,
                    "R_d90_um": emul_result.d90 * 0.5e6,
                    "note": "L2 run on d50; pore size scales ~sqrt(R) across distribution",
                },
            },
            "level2": {
                "pore_size_mean_nm": gel_result.pore_size_mean * 1e9,
                "porosity": gel_result.porosity,
                "alpha_final": gel_result.alpha_final,
            },
            "level3": {
                "p_final": xlink_result.p_final,
                "G_chitosan_Pa": xlink_result.G_chitosan_final,
                "xi_nm": xlink_result.xi_final * 1e9,
            },
            "level4": {
                "G_agarose_Pa": mech_result.G_agarose,
                "G_chitosan_Pa": mech_result.G_chitosan,
                "G_DN_Pa": mech_result.G_DN,
                "E_star_Pa": mech_result.E_star,
            },
        }
        with open(run_dir / "summary.json", "w") as f:
            json.dump(summary, f, indent=2)

        return FullResult(
            parameters=params,
            emulsification=emul_result,
            gelation=gel_result,
            crosslinking=xlink_result,
            mechanical=mech_result,
        )

    def run_rpm_sweep(self, rpms: list[float],
                      base_params: Optional[SimulationParameters] = None,
                      phi_d: float = None) -> list[dict]:
        """Run Level 1 for a list of RPM values."""
        base_params = base_params or SimulationParameters()
        if phi_d is None:
            phi_d = base_params.formulation.phi_d
        results = []

        props = self.db.update_for_conditions(
            T_oil=base_params.formulation.T_oil,
            c_agarose=base_params.formulation.c_agarose,
            c_chitosan=base_params.formulation.c_chitosan,
            c_span80=base_params.formulation.c_span80,
        )
        pbe = PBESolver(
            n_bins=base_params.solver.l1_n_bins,
            d_min=base_params.solver.l1_d_min,
            d_max=base_params.solver.l1_d_max,
        )

        for rpm in rpms:
            params = copy.deepcopy(base_params)
            params.emulsification.rpm = rpm
            t0 = time.perf_counter()
            emul = pbe.solve(params, props, phi_d=phi_d)
            dt = time.perf_counter() - t0
            results.append({
                "rpm": rpm, "d32_um": emul.d32 * 1e6,
                "d50_um": emul.d50 * 1e6, "span": emul.span,
                "converged": bool(emul.converged), "elapsed_s": dt,
            })
            logger.info("RPM=%d  d32=%.2f µm  (%.1fs)", rpm, emul.d32 * 1e6, dt)

        return results
