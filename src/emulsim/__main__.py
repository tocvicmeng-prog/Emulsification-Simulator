"""CLI entry point for the emulsification simulation."""

import argparse
import json
import logging
import sys
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(
        prog="emulsim",
        description="Multi-scale emulsification simulation for hydrogel microsphere preparation",
    )
    sub = parser.add_subparsers(dest="command")

    # run command
    run_p = sub.add_parser("run", help="Run the full L1->L2->L3->L4 pipeline")
    run_p.add_argument("config", nargs="?", default=None,
                       help="Path to TOML config file (default: configs/default.toml)")
    run_p.add_argument("--rpm", type=float, help="Override RPM")
    run_p.add_argument("--phi-d", type=float, help="Override dispersed phase volume fraction")
    run_p.add_argument(
        "--polymer-family", dest="polymer_family", default=None,
        choices=["agarose_chitosan", "alginate", "cellulose", "plga"],
        help="Override polymer family (Node F1-a Phase 2c). Routes the "
             "pipeline to the matching L2/L4 solver pair. Default: honour "
             "MaterialProperties from the config / property database.",
    )
    run_p.add_argument(
        "--gelant", default=None,
        choices=["cacl2_external", "gdl_caco3_internal"],
        help="Alginate gelant preset (Phase 2c polish). Applies the "
             "profile's effective Ca²⁺ concentration to "
             "formulation.c_Ca_bath using the current t_crosslink. Only "
             "meaningful with --polymer-family alginate.",
    )
    run_p.add_argument(
        "--cellulose-solvent", dest="cellulose_solvent", default=None,
        choices=["naoh_urea", "nmmo", "emim_ac", "dmac_licl"],
        help="Cellulose solvent-system preset (F1-b Phase 2 / 3). "
             "Patches MaterialProperties with the preset's χ, D, and "
             "modulus parameters. Only meaningful with "
             "--polymer-family cellulose.",
    )
    run_p.add_argument(
        "--plga-grade", dest="plga_grade", default=None,
        choices=["50_50", "75_25", "85_15", "pla"],
        help="PLGA grade preset (F1-c Phase 2). Patches "
             "MaterialProperties with the grade's D_DCM, phi_DCM_eq, "
             "G_glassy, and n_plga parameters. Only meaningful with "
             "--polymer-family plga.",
    )
    run_p.add_argument("--output", "-o", default="output", help="Output directory")
    run_p.add_argument("--quiet", "-q", action="store_true")

    # sweep command
    sweep_p = sub.add_parser("sweep", help="Run RPM sweep (Level 1 only)")
    sweep_p.add_argument("--rpm-min", type=float, default=3000)
    sweep_p.add_argument("--rpm-max", type=float, default=25000)
    sweep_p.add_argument("--rpm-steps", type=int, default=6)
    sweep_p.add_argument("config", nargs="?", default=None)

    # optimize command
    opt_p = sub.add_parser("optimize", help="Run Bayesian optimization campaign")
    opt_p.add_argument("--n-initial", type=int, default=15)
    opt_p.add_argument("--max-iter", type=int, default=200)
    opt_p.add_argument("--output", "-o", default="output/optimization")

    # design command (F3-c, v8.0 Phase 1): inverse-design BO with user targets
    des_p = sub.add_parser(
        "design",
        help="Inverse design: BO over user-specified (d32, pore, G_DN, Kav) targets",
    )
    des_p.add_argument("--d32", type=float, default=None,
                        help="Target d32 [m]. Omit to leave dimension out of spec.")
    des_p.add_argument("--d32-tol", type=float, default=None,
                        help="d32 tolerance [m]. Required if --d32 is set.")
    des_p.add_argument("--pore", type=float, default=None,
                        help="Target pore size [m].")
    des_p.add_argument("--pore-tol", type=float, default=None)
    des_p.add_argument("--G-DN", dest="G_DN", type=float, default=None,
                        help="Target shear modulus [Pa].")
    des_p.add_argument("--G-DN-log10-tol", dest="G_DN_log10_tol",
                        type=float, default=None,
                        help="log10-space tolerance for G_DN.")
    des_p.add_argument("--Kav", type=float, default=None,
                        help="Target distribution coefficient [-] (requires M3).")
    des_p.add_argument("--Kav-tol", dest="Kav_tol", type=float, default=None)
    des_p.add_argument(
        "--robust-variance-weight", dest="robust_variance_weight",
        type=float, default=0.0,
        help="F4-a: add lambda * std(obj) to mean(obj). 0 = nominal BO.",
    )
    des_p.add_argument(
        "--robust-n-samples", dest="robust_n_samples", type=int, default=5,
        help="Resamples per candidate for robust BO.",
    )
    des_p.add_argument(
        "--robust-cvar-alpha", dest="robust_cvar_alpha",
        type=float, default=0.0,
        help="F4-b: CVaR_α tail-risk aggregation. α ∈ (0, 1]. Takes "
             "precedence over --robust-variance-weight when both set. "
             "0 = nominal BO.",
    )
    des_p.add_argument("--n-initial", type=int, default=10)
    des_p.add_argument("--max-iter", type=int, default=50)
    des_p.add_argument("--output", "-o", default="output/design")

    # uncertainty command
    unc_p = sub.add_parser("uncertainty", help="Run Monte Carlo uncertainty propagation")
    unc_p.add_argument("config", nargs="?", default=None,
                        help="Path to TOML config file (default: configs/default.toml)")
    unc_p.add_argument("--n-samples", type=int, default=20,
                        help="Number of MC samples (default: 20)")
    unc_p.add_argument("--seed", type=int, default=42,
                        help="Random seed for reproducibility (default: 42)")
    unc_p.add_argument(
        "--n-jobs", type=int, default=1,
        help="Parallel MC workers (Node 15). 1 = serial; -1 = all cores. "
             "Note: process startup + Numba JIT cold-compile dominate when "
             "n_samples is small; use serial for n_samples < 8.",
    )
    unc_p.add_argument(
        "--engine", choices=["legacy", "unified"], default="unified",
        help="MC engine. 'unified' (default) runs the merged engine "
             "including any CalibrationStore posterior samples. "
             "'legacy' runs the merged engine with no posterior injection "
             "(byte-compat with v7.0.x uncertainty_core output for scripts "
             "that only want the default MaterialProperties perturbations).",
    )

    # info command
    sub.add_parser("info", help="Show default parameters and material properties")

    # ui command
    ui_p = sub.add_parser("ui", help="Launch the Streamlit web interface")
    ui_p.add_argument("--port", type=int, default=8501)

    # batch command (Node 24, audit N4): surface pipeline.batch_variability.run_batch
    batch_p = sub.add_parser(
        "batch",
        help="Run a DSD-quantile-resolved batch (L2-L4 across L1 size quantiles)",
    )
    batch_p.add_argument("config", nargs="?", default=None,
                          help="Path to TOML config file (default: configs/default.toml)")
    batch_p.add_argument(
        "--quantiles", default="0.10,0.25,0.50,0.75,0.90",
        help="Comma-separated quantiles in (0,1). Default: 0.10,0.25,0.50,0.75,0.90",
    )
    batch_p.add_argument("--output", "-o", default="output/batch",
                          help="Output directory for per-quantile FullResult dumps")
    batch_p.add_argument("--quiet", "-q", action="store_true")

    # dossier command (Node 24, audit N4): emit ProcessDossier JSON
    dossier_p = sub.add_parser(
        "dossier",
        help="Run the pipeline and write a ProcessDossier JSON artifact for reproducibility",
    )
    dossier_p.add_argument("config", nargs="?", default=None)
    dossier_p.add_argument("--output", "-o", default="output/dossier.json",
                            help="Path to write the dossier JSON")
    dossier_p.add_argument("--notes", default="",
                            help="Free-form notes attached to the dossier")
    dossier_p.add_argument("--quiet", "-q", action="store_true")

    # ingest command (Node 24, audit N4): AssayRecord -> CalibrationStore JSON
    ing_p = sub.add_parser(
        "ingest",
        help="Ingest wet-lab AssayRecord JSONs into a CalibrationStore-compatible fit JSON",
    )
    ing_p.add_argument(
        "module", choices=["L1"],
        help="Validation level to fit (only 'L1' wired in v7.0; L2/L3/L4/M2 stubs land in v7.1)",
    )
    ing_p.add_argument(
        "--assay-dir", default="data/validation/l1_dsd/assays",
        help="Directory containing AssayRecord JSON files",
    )
    ing_p.add_argument(
        "--output", "-o", default="data/validation/l1_dsd/fits/fit.json",
        help="Path to write the CalibrationStore-compatible JSON",
    )

    args = parser.parse_args()

    if args.command is None:
        parser.print_help()
        return

    # Setup logging
    level = logging.WARNING if getattr(args, 'quiet', False) else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
        datefmt="%H:%M:%S",
    )

    if args.command == "run":
        _cmd_run(args)
    elif args.command == "sweep":
        _cmd_sweep(args)
    elif args.command == "optimize":
        _cmd_optimize(args)
    elif args.command == "uncertainty":
        _cmd_uncertainty(args)
    elif args.command == "design":
        _cmd_design(args)
    elif args.command == "info":
        _cmd_info(args)
    elif args.command == "ui":
        _cmd_ui(args)
    elif args.command == "batch":
        _cmd_batch(args)
    elif args.command == "dossier":
        _cmd_dossier(args)
    elif args.command == "ingest":
        _cmd_ingest(args)


def _load_params(config_path):
    """Load params from TOML or use defaults."""
    from .datatypes import SimulationParameters
    if config_path:
        from .config import load_config
        return load_config(Path(config_path))
    return SimulationParameters()


def _cmd_run(args):
    from .pipeline.orchestrator import PipelineOrchestrator
    from .properties.database import PropertyDatabase

    params = _load_params(args.config)

    if args.rpm is not None:
        params.emulsification.rpm = args.rpm

    # phi_d: use CLI override if given, otherwise let orchestrator resolve
    # per mode (volumetric for stirred-vessel, formulation.phi_d for legacy)
    phi_d = args.phi_d if args.phi_d is not None else None

    # Node F1-a Phase 2c: --polymer-family overrides the default dispatch.
    props_overrides = None
    if getattr(args, "polymer_family", None) is not None:
        from .datatypes import PolymerFamily
        props_overrides = {"polymer_family": PolymerFamily(args.polymer_family)}

    # Node F1-a Phase 2c polish: --gelant applies a reagent-library preset
    # to formulation.c_Ca_bath (and resolves effective internal-release
    # concentration using the current t_crosslink).
    if getattr(args, "gelant", None) is not None:
        from .reagent_library_alginate import (
            GELANTS_ALGINATE,
            effective_bath_concentration,
        )
        profile = GELANTS_ALGINATE[args.gelant]
        t_end = params.formulation.t_crosslink
        params.formulation.c_Ca_bath = effective_bath_concentration(
            profile, t_end=t_end,
        )
        if not args.quiet:
            print(
                f"Gelant preset: {profile.name}  "
                f"(c_Ca_bath = {params.formulation.c_Ca_bath:.1f} mol/m³ "
                f"at t_crosslink = {t_end:.0f} s)"
            )

    # Node F1-b Phase 2: --cellulose-solvent folds a solvent-system
    # preset into props_overrides so that run_single receives the
    # correct χ, D, and modulus parameters.
    if getattr(args, "cellulose_solvent", None) is not None:
        from .properties.cellulose_defaults import CELLULOSE_SOLVENT_PRESETS
        preset = CELLULOSE_SOLVENT_PRESETS[args.cellulose_solvent]
        if props_overrides is None:
            props_overrides = {}
        props_overrides.update({
            "N_p_cellulose": preset.N_p,
            "chi_PS_cellulose": preset.chi_PS,
            "chi_PN_cellulose": preset.chi_PN,
            "chi_SN_cellulose": preset.chi_SN,
            "D_solvent_cellulose": preset.D_solvent,
            "D_nonsolvent_cellulose": preset.D_nonsolvent,
            "kappa_CH_cellulose": preset.kappa_CH,
            "K_cell_modulus": preset.K_cell,
            "alpha_cell_modulus": preset.alpha_cell,
        })
        if not args.quiet:
            print(
                f"Cellulose solvent preset: {preset.name}  "
                f"(N_p={preset.N_p}, χ_PS={preset.chi_PS}, "
                f"K_cell={preset.K_cell:.1e} Pa)"
            )

    # Node F1-c Phase 2: --plga-grade folds a PLGA grade preset into
    # props_overrides.
    if getattr(args, "plga_grade", None) is not None:
        from .properties.plga_defaults import PLGA_GRADE_PRESETS
        grade = PLGA_GRADE_PRESETS[args.plga_grade]
        if props_overrides is None:
            props_overrides = {}
        props_overrides.update({
            "D_DCM_plga": grade.D_DCM,
            "phi_DCM_eq": grade.phi_DCM_eq,
            "G_glassy_plga": grade.G_glassy,
            "n_plga_modulus": grade.n_plga,
        })
        if not args.quiet:
            print(
                f"PLGA grade preset: {grade.name}  "
                f"(D_DCM={grade.D_DCM:.1e} m²/s, "
                f"G_glassy={grade.G_glassy:.1e} Pa)"
            )

    orch = PipelineOrchestrator(output_dir=Path(args.output))
    result = orch.run_single(params, phi_d=phi_d, props_overrides=props_overrides)

    e = result.emulsification
    g = result.gelation
    x = result.crosslinking
    m = result.mechanical

    print()
    print("=== Simulation Results ===")
    print(f"  L1 Emulsification:  d32 = {e.d32*1e6:.2f} um   span = {e.span:.2f}")
    print(f"  L2 Gelation:        pore = {g.pore_size_mean*1e9:.1f} nm   porosity = {g.porosity:.3f}")
    print(f"  L3 Crosslinking:    p = {x.p_final:.3f}   G_chit = {x.G_chitosan_final:.0f} Pa")
    print(f"  L4 Mechanical:      G_DN = {m.G_DN:.0f} Pa   E* = {m.E_star:.0f} Pa")
    print()


def _cmd_sweep(args):
    import numpy as np
    from .pipeline.orchestrator import PipelineOrchestrator

    params = _load_params(args.config)
    if args.rpm_steps < 1:
        args.rpm_steps = 1
    rpms = np.linspace(args.rpm_min, args.rpm_max, args.rpm_steps).tolist()

    orch = PipelineOrchestrator()
    results = orch.run_rpm_sweep(rpms, base_params=params)

    print()
    print(f"{'RPM':>8s}  {'d32 (um)':>10s}  {'d50 (um)':>10s}  {'span':>6s}  {'time (s)':>8s}")
    print("-" * 50)
    for r in results:
        print(f"{r['rpm']:8.0f}  {r['d32_um']:10.2f}  {r['d50_um']:10.2f}  {r['span']:6.2f}  {r['elapsed_s']:8.1f}")
    print()


def _cmd_optimize(args):
    from .optimization.engine import OptimizationEngine
    from .optimization.analysis import pareto_summary, best_compromise

    engine = OptimizationEngine(
        n_initial=args.n_initial,
        max_iterations=args.max_iter,
        output_dir=Path(args.output),
    )
    state = engine.run()
    print()
    print(pareto_summary(state))
    idx = best_compromise(state)
    print(f"\nBest compromise: Pareto point #{idx+1}")
    print(f"  objectives = {state.pareto_Y[idx]}")
    print()


def _cmd_design(args):
    """Inverse-design BO with user-specified target spec (F3-c, v8.0).

    Builds a ``TargetSpec`` from CLI flags, constructs an
    ``OptimizationEngine(target_spec=...)``, and runs. The engine's
    legacy Pareto filter and trust-aware gating remain active; the
    only change is what "good" means.

    F4-a: --robust-variance-weight lambda turns on mean+lambda*std
    objective evaluation for robust BO.
    """
    from .optimization.engine import OptimizationEngine
    from .optimization.objectives import TargetSpec
    from .optimization.analysis import pareto_summary, best_compromise

    target = TargetSpec(
        d32_target=args.d32, d32_tol=args.d32_tol,
        pore_target=args.pore, pore_tol=args.pore_tol,
        G_DN_target=args.G_DN, G_DN_log10_tol=args.G_DN_log10_tol,
        Kav_target=args.Kav, Kav_tol=args.Kav_tol,
    )
    try:
        target.validate()
    except ValueError as exc:
        raise SystemExit(f"Invalid --design target: {exc}")

    print(f"Design targets: active dims = {target.active_dims()}")
    if args.robust_variance_weight > 0:
        print(
            f"Robust BO: lambda = {args.robust_variance_weight}, "
            f"n_samples/candidate = {args.robust_n_samples}"
        )

    engine = OptimizationEngine(
        n_initial=args.n_initial,
        max_iterations=args.max_iter,
        output_dir=Path(args.output),
        target_spec=target,
        robust_variance_weight=args.robust_variance_weight,
        robust_n_samples=args.robust_n_samples,
    )
    state = engine.run()
    print()
    print(pareto_summary(state))
    idx = best_compromise(state)
    print(f"\nBest compromise: Pareto point #{idx+1}")
    print(f"  objectives = {state.pareto_Y[idx]}")
    print()


def _cmd_uncertainty(args):
    """MC uncertainty propagation.

    Node 30 (v7.1): both ``--engine unified`` and ``--engine legacy`` run
    the single merged UnifiedUncertaintyEngine. ``unified`` includes any
    CalibrationStore posterior samples; ``legacy`` constructs a fresh
    engine with ``calibration_store=None`` for byte-compat with v7.0.x
    scripts that expect the default MaterialProperties perturbations
    only.
    """
    from .uncertainty_unified import UnifiedUncertaintyEngine

    params = _load_params(args.config)
    engine = UnifiedUncertaintyEngine()  # no calibration store wired in CLI yet
    label = "legacy" if args.engine == "legacy" else "unified"
    print(f"Running unified MC ({args.n_samples} samples, "
          f"n_jobs={args.n_jobs}, engine={label})...")
    result = engine.run_m1l4(
        params,
        n_samples=args.n_samples,
        seed=args.seed,
        n_jobs=args.n_jobs,
    )
    print()
    print(result.summary())
    print()


def _cmd_info(args):
    from .datatypes import SimulationParameters, MaterialProperties

    params = SimulationParameters()
    props = MaterialProperties()

    print()
    print("=== Default Simulation Parameters ===")
    print(f"  RPM:              {params.emulsification.rpm:.0f}")
    print(f"  Emul. time:       {params.emulsification.t_emulsification:.0f} s")
    print(f"  Agarose conc:     {params.formulation.c_agarose:.1f} kg/m3 ({params.formulation.c_agarose/10:.1f}% w/v)")
    print(f"  Chitosan conc:    {params.formulation.c_chitosan:.1f} kg/m3 ({params.formulation.c_chitosan/10:.1f}% w/v)")
    print(f"  Span-80 conc:     {params.formulation.c_span80:.1f} kg/m3 ({params.formulation.c_span80/10:.1f}% w/v)")
    print(f"  Oil temp:         {params.formulation.T_oil - 273.15:.0f} C")
    print(f"  Cooling rate:     {params.formulation.cooling_rate * 60:.1f} C/min")
    print(f"  Genipin conc:     {params.formulation.c_genipin:.1f} mM")
    print(f"  Crosslink temp:   {params.formulation.T_crosslink - 273.15:.0f} C")
    print(f"  Crosslink time:   {params.formulation.t_crosslink / 3600:.0f} h")
    print(f"  phi_d:            {params.formulation.phi_d}")
    print()
    print("=== Material Properties ===")
    print(f"  Oil density:      {props.rho_oil:.0f} kg/m3 (at 20 C)")
    print(f"  Oil viscosity:    {props.mu_oil*1000:.1f} mPa.s (at 90 C)")
    print(f"  IFT (Span-80):    {props.sigma*1000:.1f} mN/m")
    print(f"  T_gel:            {props.T_gel - 273.15:.0f} C")
    print(f"  Genipin Ea:       {props.E_a_xlink/1000:.0f} kJ/mol")
    print(f"  Bridge eff:       {props.f_bridge:.1%}")
    print(f"  IPN coupling:     eta = {props.eta_coupling}")
    print()


def _cmd_ui(args):
    import subprocess
    app_path = Path(__file__).parent / "visualization" / "app.py"
    print(f"Launching EmulSim UI on port {args.port}...")
    print(f"Open http://localhost:{args.port} in your browser")
    subprocess.run([
        sys.executable, "-m", "streamlit", "run", str(app_path),
        "--server.port", str(args.port),
        "--server.headless", "true",
    ])


# ─── Node 24 (v7.0.1, audit N4): batch / dossier / ingest CLI surfacing ────


def _parse_quantiles_arg(s: str) -> tuple[float, ...]:
    """Parse the --quantiles "0.10,0.25,0.50" string into a sorted tuple.

    Argparse error messages here mean the user typoed; raise a clear
    ValueError that argparse converts into a top-level exit message.
    """
    try:
        vals = [float(x) for x in s.split(",") if x.strip()]
    except ValueError as exc:
        raise ValueError(
            f"Could not parse --quantiles {s!r}: expected comma-separated "
            "floats in (0,1), e.g. '0.10,0.50,0.90'"
        ) from exc
    if not vals:
        raise ValueError("--quantiles must include at least one value.")
    if not all(0.0 < v < 1.0 for v in vals):
        raise ValueError(
            f"--quantiles must all lie in (0,1); got {vals!r}"
        )
    return tuple(vals)


def _cmd_batch(args):
    """Surface ``pipeline.batch_variability.run_batch`` on the CLI (audit N4)."""
    from .pipeline.batch_variability import run_batch

    quantiles = _parse_quantiles_arg(args.quantiles)
    params = _load_params(args.config)

    result = run_batch(
        params,
        quantiles=quantiles,
        output_dir=Path(args.output),
    )

    print()
    print("=== Batch Variability Results ===")
    print(f"  Quantiles:           {quantiles}")
    print(f"  Mass-weighted mean d32:  {result.mean_d32_m * 1e6:.2f} um")
    print(f"  Mass-weighted mean pore: {result.mean_pore_m * 1e9:.1f} nm")
    print(f"  Pore p5 / p50 / p95:     "
          f"{result.pore_p5_m*1e9:.1f} / "
          f"{result.pore_p50_m*1e9:.1f} / "
          f"{result.pore_p95_m*1e9:.1f} nm")
    print(f"  Mass-weighted mean G_DN: {result.mean_G_DN_Pa:.0f} Pa")
    print(f"  G_DN p5 / p50 / p95:     "
          f"{result.G_DN_p5_Pa:.0f} / "
          f"{result.G_DN_p50_Pa:.0f} / "
          f"{result.G_DN_p95_Pa:.0f} Pa")
    print()
    print("Per-quantile representative bead radii [um]:")
    for q, r, w in zip(quantiles, result.quantile_radii_m,
                        result.quantile_mass_fractions):
        print(f"  q={q:.3f}  R={r*1e6:.2f} um  mass_fraction={w:.3f}")
    print()


def _cmd_dossier(args):
    """Run the pipeline and write a ProcessDossier JSON (audit N4)."""
    from .pipeline.orchestrator import PipelineOrchestrator
    from .process_dossier import ProcessDossier

    params = _load_params(args.config)
    out_path = Path(args.output)
    run_dir = out_path.parent / "_dossier_runs"

    orch = PipelineOrchestrator(output_dir=run_dir)
    full_result = orch.run_single(params)

    dossier = ProcessDossier.from_run(
        full_result,
        notes=args.notes,
    )
    written = dossier.export_json(out_path)

    print()
    print(f"=== Process Dossier ===")
    print(f"  run_id:                {dossier.run_id}")
    print(f"  timestamp_utc:         {dossier.timestamp_utc}")
    print(f"  evidence_tier:         {full_result.run_report.min_evidence_tier}")
    print(f"  trust_level:           {full_result.run_report.trust_level}")
    print(f"  models recorded:       {len(full_result.run_report.model_graph)}")
    print(f"  calibrations applied:  "
          f"{full_result.run_report.diagnostics.get('calibration_count', 0)}")
    print(f"  dossier written to:    {written}")
    print()


def _cmd_ingest(args):
    """Wet-lab AssayRecord ingest -> CalibrationStore JSON (audit N4)."""
    from .calibration.fitters import (
        fit_l1_dsd_to_calibration_entries,
        load_assay_records,
        write_calibration_json,
    )

    if args.module != "L1":
        # Argparse already constrains choices, but be explicit for v7.1.
        print(f"Module {args.module} ingestion is not yet wired (v7.1).",
              file=sys.stderr)
        sys.exit(2)

    assay_dir = Path(args.assay_dir)
    if not assay_dir.exists():
        print(f"ERROR: --assay-dir {assay_dir} does not exist.", file=sys.stderr)
        print("Place AssayRecord JSONs there or pass --assay-dir.",
              file=sys.stderr)
        sys.exit(3)

    records = load_assay_records(assay_dir)
    if not records:
        print(f"WARNING: no AssayRecord JSONs found in {assay_dir}.")
        print("Nothing to fit; skipping write.")
        return

    entries = fit_l1_dsd_to_calibration_entries(records)
    if not entries:
        print(f"WARNING: fitter produced 0 CalibrationEntries from "
              f"{len(records)} records (no DROPLET_SIZE_DISTRIBUTION assays?).")
        return

    out_path = Path(args.output)
    fit_metadata = {
        "module": args.module,
        "assay_dir": str(assay_dir),
        "n_records_in": len(records),
        "n_entries_out": len(entries),
        "fitter": "stub_mean (v7.0)",
    }
    written = write_calibration_json(entries, out_path, fit_metadata=fit_metadata)
    print()
    print(f"=== AssayRecord Ingest ===")
    print(f"  module:               {args.module}")
    print(f"  records read:         {len(records)}")
    print(f"  calibration entries:  {len(entries)}")
    print(f"  fit JSON written to:  {written}")
    print(f"  metadata sidecar:     {written.with_suffix('.meta.json')}")
    print()
    print("Load with:")
    print(f"  store = CalibrationStore(); store.load_json({str(written)!r})")
    print(f"  ctx = RunContext(calibration_store=store)")
    print(f"  PipelineOrchestrator().run_single(params, run_context=ctx)")
    print()


if __name__ == "__main__":
    main()
