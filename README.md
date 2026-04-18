# EmulSim

[![CI](https://github.com/tocvicmeng-prog/Emulsification-Simulator/actions/workflows/ci.yml/badge.svg)](https://github.com/tocvicmeng-prog/Emulsification-Simulator/actions/workflows/ci.yml)

Multi-scale simulation of emulsification for macroporous double-network
polysaccharide hydrogel microsphere preparation.

EmulSim models the complete fabrication pipeline: from rotor-stator
emulsification of an agarose-chitosan aqueous phase in hot oil, through
thermally-induced gelation and spinodal decomposition, multi-mechanism
crosslinking of the chitosan network, to final double-network mechanical
properties. The target application is producing hydrogel microspheres with
controlled pore size for cell culture scaffolds, chromatography beads, and
drug delivery.

## Features

- **8 crosslinkers** (Genipin, Glutaraldehyde, EDC/NHS, PEGDA+UV, TPP, Epichlorohydrin, DVS, Citric Acid) with 4 kinetics models
- **6 surfactants** (Span-80, Span-60, Span-40, Span-85, PGPR, Lecithin) with Szyszkowski-Langmuir IFT model
- **Streamlit Web UI** with interactive reagent dropdowns, sliders, and Plotly dashboards
- **Bayesian optimization** (BoTorch) for multi-objective parameter tuning
- **Monte Carlo uncertainty propagation** for confidence intervals
- **Calibration protocol** -- 5-study wet-lab protocol to calibrate simulation constants (docs/04)

## Quick Start

### Install

```bash
git clone https://github.com/tocvicmeng-prog/emulsification-similation-for-microspheres-preparation.git
cd emulsification-similation-for-microspheres-preparation
pip install -e .

# With optimization support (BoTorch/PyTorch):
pip install -e ".[optimization]"

# Everything:
pip install -e ".[all]"
```

### Run

```bash
# Full pipeline with default parameters
python -m emulsim run

# Full pipeline with custom config
python -m emulsim run configs/default.toml

# Override specific parameters
python -m emulsim run --rpm 15000 --phi-d 0.08

# Show default parameters and material properties
python -m emulsim info

# Launch Streamlit web UI
python -m emulsim ui
```

### Example Output

The fast smoke config (`configs/fast_smoke.toml`, ~0.2 s) reliably produces:

```
=== Simulation Results ===
  L1 Emulsification:  d32 = 22.08 um   span = 1.04
  L2 Gelation:        pore = 180.9 nm  porosity = 0.871
  L3 Crosslinking:    p = 0.040        G_chit = 2062 Pa
  L4 Mechanical:      G_DN = 70766 Pa  E* = 257332 Pa
```

These are uncalibrated **semi-quantitative** outputs — the platform now ships
a `RunReport` with each result that carries the weakest model evidence tier
(see `docs/quickstart.md` for tier definitions and the new
`RunContext`/`CalibrationStore` calibration-injection workflow).

### RPM Sweep

```bash
python -m emulsim sweep --rpm-min 3000 --rpm-max 25000 --rpm-steps 6
```

### Bayesian Optimization

```bash
python -m emulsim optimize --n-initial 15 --max-iter 200
```

### Uncertainty Quantification

```bash
python -m emulsim uncertainty --n-samples 20
```

### Web UI

```bash
python -m emulsim ui
# Open http://localhost:8501
```

The UI provides:
- Dropdown selectors for crosslinker and surfactant with per-reagent defaults
- Per-constant Literature/Custom toggle with links to calibration protocol
- Interactive sliders for all process parameters
- Results dashboard with Plotly charts (size distribution, phase field, kinetics, Hertz contact, Kav)
- Trust assessment with 10 automated reliability checks
- Optimization assessment with actionable recommendations

## Pipeline Architecture

```
 Input Parameters (TOML config)
         |
         v
 +-------------------+
 | L1: Emulsification|   Population Balance Equation (PBE)
 | Rotor-stator      |   Breakage/coalescence kernels
 | droplet sizing    |   -> d32, d50, span, DSD
 +-------------------+
         |  d50
         v
 +-------------------+
 | L2: Gelation &    |   2D Cahn-Hilliard phase-field
 | Pore Formation    |   Flory-Huggins thermodynamics
 |                   |   Avrami gelation arrest
 |                   |   -> pore size, porosity
 +-------------------+
         |  porosity, R_droplet
         v
 +-------------------+
 | L3: Crosslinking  |   Multi-mechanism ODE kinetics
 | Kinetics          |   8 crosslinkers, Arrhenius rate law
 |                   |   -> conversion p, G_chitosan, mesh size
 +-------------------+
         |  G_chitosan, mesh size
         v
 +-------------------+
 | L4: Mechanical    |   Double-network (DN) model
 | Properties        |   IPN coupling, Hertz contact
 |                   |   -> G_DN, E*, sieving curve
 +-------------------+
```

## Configuration

All parameters are set via TOML files. See `configs/default.toml` for the
full reference configuration with comments, or `docs/configuration.md` for
a complete parameter reference with units and physical meaning.

Key sections:
- `[emulsification]` -- RPM, emulsification time, mixer geometry
- `[formulation]` -- polymer concentrations, temperatures, crosslinking conditions, phi_d
- `[solver]` -- numerical resolution for each level
- `[optimization]` -- Bayesian optimization campaign settings

## Python API

```python
from emulsim import run_pipeline, SimulationParameters

# Default parameters
result = run_pipeline()

# Custom parameters
params = SimulationParameters()
params.emulsification.rpm = 15000
params.formulation.phi_d = 0.08
result = run_pipeline(params)

print(f"d32 = {result.emulsification.d32 * 1e6:.1f} um")
print(f"pore = {result.gelation.pore_size_mean * 1e9:.0f} nm")
print(f"G_DN = {result.mechanical.G_DN:.0f} Pa")
```

## CLI Commands

| Command | Description |
|---------|-------------|
| `python -m emulsim run` | Run full 4-level pipeline |
| `python -m emulsim info` | Display parameters and material properties |
| `python -m emulsim sweep` | RPM parameter sweep |
| `python -m emulsim optimize` | BoTorch Bayesian optimization |
| `python -m emulsim uncertainty` | Monte Carlo uncertainty propagation |
| `python -m emulsim ui` | Launch Streamlit web interface |

## Project Structure

```
emulsification_sim/
  configs/
    default.toml              # Reference configuration
  data/
    properties.toml           # Material property database
  docs/
    quickstart.md             # Step-by-step first run guide
    configuration.md          # Full TOML parameter reference
    01_scientific_advisor_report.md
    02_computational_architecture.md
    03_scientific_review.md
    04_calibration_protocol.md
  src/emulsim/
    __init__.py               # Top-level API (run_pipeline, key types)
    __main__.py               # CLI entry point (python -m emulsim)
    config.py                 # TOML config loader
    datatypes.py              # Dataclasses for params, results
    pipeline/
      orchestrator.py         # L1->L2->L3->L4 sequencing
    level1_emulsification/    # PBE solver, breakage/coalescence kernels
    level2_gelation/          # Cahn-Hilliard 2D solver, pore analysis
    level3_crosslinking/      # Multi-mechanism crosslinking kinetics
    level4_mechanical/        # DN mechanical model
    properties/               # T-dependent viscosity, IFT, thermodynamics
    optimization/             # BoTorch multi-objective optimization
    visualization/            # Plotting utilities
    ui/                       # Streamlit web interface
  tests/
  notebooks/
  CHANGELOG.md
```

## Documentation

Detailed design documents are in `docs/`:

1. **Quickstart Guide** (`docs/quickstart.md`) --
   Step-by-step guide to your first simulation.
2. **Configuration Reference** (`docs/configuration.md`) --
   Full TOML parameter reference with units and physical meaning.
3. **Scientific Advisor Report** (`docs/01_scientific_advisor_report.md`) --
   Physical models, governing equations, and parameter choices.
4. **Computational Architecture** (`docs/02_computational_architecture.md`) --
   Software design, data flow, and solver implementations.
5. **Scientific Review** (`docs/03_scientific_review.md`) --
   Validation strategy, known limitations, and future directions.
6. **Calibration Protocol** (`docs/04_calibration_protocol.md`) --
   5-study wet-lab protocol (1081 lines) to calibrate simulation constants.

## Requirements

- Python >= 3.11
- NumPy, SciPy, h5py, Matplotlib
- Optional: BoTorch + PyTorch (for optimization), Plotly + Jupyter (for interactive viz)
- Optional: Streamlit (for web UI)

## License

See `LICENSE` for details.
