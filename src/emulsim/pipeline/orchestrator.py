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
    KernelConfig,
    M1ExportContract,
    MaterialProperties,
    RunContext,
    RunReport,
    SimulationParameters,
)
from ..trust import TrustAssessment, assess_trust
from ..properties.database import PropertyDatabase
from ..level1_emulsification.solver import PBESolver
from ..level2_gelation.solver import solve_gelation, solve_gelation_timing
from ..level3_crosslinking.solver import (
    solve_crosslinking,
    available_amine_concentration,
    available_hydroxyl_concentration,
)
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
                   phi_d: float = None,
                   l2_mode: str = 'empirical',
                   props_overrides: dict | None = None,
                   crosslinker_key: str = 'genipin',
                   uv_intensity: float = 0.0,
                   run_context: RunContext | None = None) -> FullResult:
        """Execute the full pipeline for a single parameter set.

        Parameters
        ----------
        l2_mode : str
            'empirical' (default) or 'ch_2d' for mechanistic phase-field.
        props_overrides : dict, optional
            Key-value pairs to override on the MaterialProperties object
            after it is built from the PropertyDatabase.  Used by the UI
            to inject reagent-library parameters (e.g. crosslinker kinetics,
            pre-computed IFT).
        run_context : RunContext, optional
            Cross-cutting run inputs (Node 7, v6.1). When provided with a
            populated ``calibration_store``, the orchestrator applies
            calibration overrides to the L1 KernelConfig and to
            MaterialProperties for L2-L4 BEFORE solving. ``None`` reproduces
            the v6.0 behaviour (no calibration injection).
        """
        errors = params.validate()
        if errors:
            raise ValueError("Invalid parameters:\n" + "\n".join(f"  - {e}" for e in errors))

        if phi_d is None:
            if params.emulsification.mode == "stirred_vessel":
                phi_d = params.formulation.phi_d_from_volumes
            else:
                phi_d = params.formulation.phi_d

        # Stirred-vessel: sync c_span80 from volumetric % if the TOML config
        # set c_span80_vol_pct but left c_span80 at the FormulationParameters
        # default. Compare against the volumetric-derived value to detect
        # whether the user explicitly set c_span80 to something different.
        if params.emulsification.mode == "stirred_vessel":
            c_from_vol = params.formulation.c_span80_from_vol_pct
            # If c_span80 does NOT match c_from_vol, the user set it explicitly
            # (e.g., from a calibrated run). Honour the explicit value.
            # If it's the legacy default (20.0) or matches the TOML c_span80
            # field that was loaded alongside vol_pct, overwrite with volumetric.
            if abs(params.formulation.c_span80 - c_from_vol) > 0.5:
                # User explicitly set c_span80 different from vol_pct derivation
                pass  # honour explicit value
            else:
                params.formulation.c_span80 = c_from_vol

        run_id = params.run_id or datetime.now().strftime("%Y%m%d_%H%M%S_%f")
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

        # Extract eta_is_custom before applying props overrides (it is not a
        # MaterialProperties field and must not be set on props).
        _eta_is_custom = bool(props_overrides.pop("eta_is_custom", False)) if props_overrides else False

        # Apply caller-supplied overrides (reagent library selections, etc.)
        if props_overrides:
            for k, v in props_overrides.items():
                setattr(props, k, v)

        # Node 7 (v6.1): apply calibration overrides from RunContext.
        # Order matters: caller props_overrides win for explicitly-set fields,
        # then calibration entries layer on top for any matching attribute.
        # The reasoning: explicit user/UI overrides represent a deliberate
        # choice for THIS run, while calibration entries represent the
        # current best evidence for the underlying constant.
        applied_calibrations: list[str] = []
        if run_context is not None and run_context.calibration_store is not None:
            store = run_context.calibration_store
            # Material properties pick up L2/L3/L4 calibrations (rate constants,
            # pore coefficients, IPN coupling). Each call rewrites only matching
            # attributes and returns a deep copy + override log.
            for module in ("L2", "L3", "L4"):
                if store.has_calibration_for(module):
                    props, overrides_log = store.apply_to_model_params(module, props)
                    applied_calibrations.extend(overrides_log)

            # L1 calibration targets KernelConfig. Build a working copy on
            # params.emulsification.kernels so the legacy `solve()` dispatch
            # (Node 3) honours it. Done after props because L1 may also read
            # props.breakage_C3 — but the dispatch wires KernelConfig as the
            # primary source.
            if store.has_calibration_for("L1"):
                # Audit N1 (v7.0.1): use a deep copy so the caller's
                # SimulationParameters instance is not silently mutated.
                # This matters for callers that reuse `params` across
                # multiple run_single() calls (e.g. batch_variability.run_batch,
                # parameter sweeps, optimisation campaigns) — without the copy
                # the SECOND call would see calibrated kernels left over
                # from the FIRST.
                source_kernels = (
                    params.emulsification.kernels
                    or KernelConfig.for_rotor_stator_legacy()
                )
                kernels_copy = copy.deepcopy(source_kernels)
                kernels_copy, overrides_log = store.apply_to_model_params(
                    "L1", kernels_copy,
                )
                applied_calibrations.extend(overrides_log)
                # Build a new EmulsificationParameters carrying the calibrated
                # kernels so we don't write back into the caller's params.
                params = copy.copy(params)
                params.emulsification = copy.copy(params.emulsification)
                params.emulsification.kernels = kernels_copy

            if applied_calibrations:
                logger.info(
                    "Applied %d calibration override(s) from RunContext.",
                    len(applied_calibrations),
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
        d_mode_um = getattr(emul_result, 'd_mode', 0.0) * 1e6
        logger.info("L1 done: d32=%.2f µm, d_mode=%.2f µm, span=%.2f (%.1fs)",
                     emul_result.d32 * 1e6, d_mode_um, emul_result.span, timings["L1"])

        # ── Level 2a: Gelation Timing ──────────────────────────────────
        # Use d50 (volume median) as the representative droplet size for L2.
        # d32 (Sauter mean) overweights large droplets; d50 better represents the
        # population center.  The full size distribution (d10–d90) is recorded in
        # the summary so downstream users can estimate pore-size spread via the
        # Cahn-Hilliard scaling λ* ~ sqrt(R).
        R_droplet = emul_result.d50 / 2.0
        logger.info("L2a: Gelation timing (R=%.2f um)", R_droplet * 1e6)
        timing_result = solve_gelation_timing(params, props, R_droplet=R_droplet)
        logger.info("L2a done: t_gel=%.1fs, alpha=%.3f, cool_eff=%.4f K/s",
                     timing_result.t_gel_onset, timing_result.alpha_final,
                     timing_result.cooling_rate_effective)

        # ── Level 2b: Gelation & Pore Formation ─────────────────────────
        # Node 8 (F8): pass the L2a timing result so the empirical pore model
        # can use the actual Avrami alpha_final instead of the hardcoded 0.999.
        # Mechanistic modes (ch_2d, ch_ternary) ignore the argument.
        logger.info("L2b: Gelation (R=%.2f µm, d50-based)", R_droplet * 1e6)
        t0 = time.perf_counter()
        gel_result = solve_gelation(
            params, props, R_droplet=R_droplet, mode=l2_mode, timing=timing_result,
        )
        timings["L2"] = time.perf_counter() - t0
        logger.info("L2 done: pore=%.1f nm, porosity=%.2f (%.1fs)",
                     gel_result.pore_size_mean * 1e9, gel_result.porosity, timings["L2"])

        # ── Level 3: Crosslinking Kinetics ───────────────────────────────
        logger.info("L3: Crosslinking (%s=%.1f mM, T=%.0f°C, t=%.0fh)",
                     crosslinker_key,
                     params.formulation.c_genipin,
                     params.formulation.T_crosslink - 273.15,
                     params.formulation.t_crosslink / 3600)
        t0 = time.perf_counter()
        xlink_result = solve_crosslinking(
            params, props,
            R_droplet=R_droplet,
            porosity=gel_result.porosity,
            crosslinker_key=crosslinker_key,
            uv_intensity=uv_intensity,
        )
        timings["L3"] = time.perf_counter() - t0
        logger.info("L3 done: p=%.2f, G_chit=%.0f Pa, xi=%.1f nm (%.3fs)",
                     xlink_result.p_final, xlink_result.G_chitosan_final,
                     xlink_result.xi_final * 1e9, timings["L3"])

        # ── Level 4: Mechanical Properties ───────────────────────────────
        logger.info("L4: Mechanical property prediction")
        t0 = time.perf_counter()
        mech_result = solve_mechanical(params, props, gel_result, xlink_result, R_droplet=R_droplet, eta_is_custom=_eta_is_custom)
        timings["L4"] = time.perf_counter() - t0
        logger.info("L4 done: G_DN=%.0f Pa, E*=%.0f Pa (%.4fs)",
                     mech_result.G_DN, mech_result.E_star, timings["L4"])

        total = sum(timings.values())
        logger.info("Pipeline complete: %.1fs total", total)

        # Trust assessment
        full_result = FullResult(
            parameters=params,
            emulsification=emul_result,
            gelation=gel_result,
            crosslinking=xlink_result,
            mechanical=mech_result,
            gelation_timing=timing_result,
        )
        trust = assess_trust(full_result, params, props, crosslinker_key=crosslinker_key, l2_mode=l2_mode)
        if not trust.trustworthy:
            logger.warning("TRUST GATE: %s", trust.summary())
        elif trust.warnings:
            logger.info("TRUST GATE: %s", trust.summary())

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
        summary["trust"] = {
            "level": trust.level,
            "trustworthy": trust.trustworthy,
            "warnings": trust.warnings,
            "blockers": trust.blockers,
        }
        with open(run_dir / "summary.json", "w") as f:
            json.dump(summary, f, indent=2)

        # Assemble RunReport (v6.1: model evidence provenance)
        _manifests = [
            r.model_manifest
            for r in (emul_result, gel_result, xlink_result, mech_result)
            if getattr(r, "model_manifest", None) is not None
        ]
        _diagnostics: dict = {"timings": timings, "total_time_s": total}
        if applied_calibrations:
            # Node 7: surface calibration provenance for downstream consumers
            # (UI badges, JSON export, optimizer Pareto labelling). Stored as
            # a list of human-readable override descriptions matching what
            # CalibrationStore.apply_to_model_params logs.
            _diagnostics["calibrations_applied"] = list(applied_calibrations)
            _diagnostics["calibration_count"] = len(applied_calibrations)
        run_report = RunReport(
            model_graph=_manifests,
            trust_level=trust.level,
            trust_warnings=list(trust.warnings),
            trust_blockers=list(trust.blockers),
            diagnostics=_diagnostics,
        )
        run_report.min_evidence_tier = run_report.compute_min_tier().value
        full_result.run_report = run_report

        return full_result

    def run_rpm_sweep(self, rpms: list[float],
                      base_params: Optional[SimulationParameters] = None,
                      phi_d: float = None) -> list[dict]:
        """Run Level 1 for a list of RPM values."""
        base_params = base_params or SimulationParameters()
        if phi_d is None:
            if base_params.emulsification.mode == "stirred_vessel":
                phi_d = base_params.formulation.phi_d_from_volumes
            else:
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


# ─── Amine-reactive crosslinker families ──────────────────────────────────────
# Crosslinkers in this set react with primary amines (NH2 on chitosan).
# Hydroxyl-reactive crosslinkers (DVS, ECH, citric acid) leave NH2 pools intact.
_AMINE_REACTIVE_FAMILIES = frozenset({"amine_covalent"})
_HYDROXYL_REACTIVE_FAMILIES = frozenset({"hydroxyl_covalent"})


def export_for_module2(result: FullResult,
                       trust: TrustAssessment,
                       crosslinker_key: str = "genipin",
                       props: Optional[MaterialProperties] = None) -> M1ExportContract:
    """Extract Module 2 inputs from Module 1 FullResult.

    Computes residual NH2 and OH concentrations from the crosslinking
    result and formulation parameters.  NH2 is only depleted when an
    amine-reactive crosslinker (e.g. genipin, glutaraldehyde) was used;
    hydroxyl-reactive crosslinkers (DVS, ECH) leave the NH2 pool intact.

    Parameters
    ----------
    result : FullResult
        Output of PipelineOrchestrator.run_single().
    trust : TrustAssessment
        Output of assess_trust() for the same run.
    crosslinker_key : str
        Name of the primary crosslinker (used for labelling only; chemistry
        family is read from result.crosslinking.network_metadata.solver_family).
    props : MaterialProperties, optional
        Material properties object containing DDA and M_GlcN.  If None,
        default values are used (DDA=0.90, M_GlcN=161.16).

    Returns
    -------
    M1ExportContract
        Fully populated stable interface object for Module 2.
    """
    params = result.parameters
    emul = result.emulsification
    gel = result.gelation
    xl = result.crosslinking
    mech = result.mechanical

    # Formulation source parameters
    c_agarose = params.formulation.c_agarose
    c_chitosan = params.formulation.c_chitosan

    # DDA and M_GlcN: prefer props (interpolated), fall back to defaults
    DDA = props.DDA if props is not None else 0.90
    M_GlcN = props.M_GlcN if props is not None else 161.16

    # ── Residual reactive groups ─────────────────────────────────────────
    # Total amine pool before any crosslinking
    nh2_total = available_amine_concentration(c_chitosan, DDA, M_GlcN)

    # Determine which reactive pool was consumed during primary crosslinking
    solver_family = (
        xl.network_metadata.solver_family
        if xl.network_metadata is not None
        else "amine_covalent"
    )

    if solver_family in _AMINE_REACTIVE_FAMILIES:
        # Amine-targeting crosslinker: NH2 is partially depleted by p_final
        nh2_residual = nh2_total * (1.0 - xl.p_final)
    else:
        # Hydroxyl-targeting or independent-network crosslinker: NH2 intact
        nh2_residual = nh2_total

    # Total hydroxyl pool (agarose OH, not consumed by amine crosslinkers)
    oh_total = available_hydroxyl_concentration(c_agarose)

    # ── Uncertainty notes (Node 0.3) ────────────────────────────────────
    notes = []
    if gel.model_tier in ("empirical_calibrated", "empirical_uncalibrated"):
        notes.append("Pore size from empirical correlation (+/-30%)")
    if mech.model_used == "phenomenological":
        notes.append("Modulus from phenomenological DN model (ranking only)")
    if xl.p_final < 0.05:
        notes.append("Low crosslinking conversion -- residual ACS estimate uncertain")
    uncertainty_notes = "; ".join(notes) if notes else "No special uncertainty flags"

    contract = M1ExportContract(
        # Geometry
        bead_radius=emul.d50 / 2.0,
        bead_d32=emul.d32,
        bead_d50=emul.d50,
        # Pore structure
        pore_size_mean=gel.pore_size_mean,
        pore_size_std=gel.pore_size_std,
        porosity=gel.porosity,
        l2_model_tier=gel.model_tier,
        # Network structure
        mesh_size_xi=xl.xi_final,
        p_final=xl.p_final,
        primary_crosslinker=crosslinker_key,
        # Residual reactive groups
        nh2_bulk_concentration=nh2_residual,
        oh_bulk_concentration=oh_total,
        # Mechanical
        G_DN=mech.G_DN,
        E_star=mech.E_star,
        model_used=mech.model_used,
        # Formulation
        c_agarose=c_agarose,
        c_chitosan=c_chitosan,
        DDA=DDA,
        # Trust
        trust_level=trust.level,
        trust_warnings=list(trust.warnings),
        uncertainty_notes=uncertainty_notes,
    )

    # Node 10 (F11): boundary unit/range check. Don't crash the pipeline on
    # violations — log as warnings so the user gets a clear diagnostic
    # without losing the run output. The orchestrator's trust gate already
    # downgrades suspicious results; this catches the more specific
    # "you sent this in the wrong unit" failure mode.
    _unit_violations = contract.validate_units()
    if _unit_violations:
        logger.warning(
            "M1ExportContract failed %d unit/range check(s):\n  %s",
            len(_unit_violations),
            "\n  ".join(_unit_violations),
        )
        # Surface in the contract's trust_warnings so downstream consumers
        # (M2 orchestrator, UI badges) see them too.
        contract.trust_warnings.extend(
            f"M1ExportContract unit check: {v}" for v in _unit_violations
        )

    return contract
