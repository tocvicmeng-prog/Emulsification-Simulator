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

    # uncertainty command
    unc_p = sub.add_parser("uncertainty", help="Run Monte Carlo uncertainty propagation")
    unc_p.add_argument("config", nargs="?", default=None,
                        help="Path to TOML config file (default: configs/default.toml)")
    unc_p.add_argument("--n-samples", type=int, default=20,
                        help="Number of MC samples (default: 20)")
    unc_p.add_argument("--seed", type=int, default=42,
                        help="Random seed for reproducibility (default: 42)")

    # info command
    sub.add_parser("info", help="Show default parameters and material properties")

    # ui command
    ui_p = sub.add_parser("ui", help="Launch the Streamlit web interface")
    ui_p.add_argument("--port", type=int, default=8501)

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
    elif args.command == "info":
        _cmd_info(args)
    elif args.command == "ui":
        _cmd_ui(args)


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

    phi_d = args.phi_d if args.phi_d is not None else params.formulation.phi_d

    orch = PipelineOrchestrator(output_dir=Path(args.output))
    result = orch.run_single(params, phi_d=phi_d)

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


def _cmd_uncertainty(args):
    from .uncertainty import UncertaintyPropagator

    params = _load_params(args.config)
    propagator = UncertaintyPropagator(n_samples=args.n_samples, seed=args.seed)

    print(f"Running Monte Carlo uncertainty propagation ({args.n_samples} samples)...")
    result = propagator.run(params)

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


if __name__ == "__main__":
    main()
