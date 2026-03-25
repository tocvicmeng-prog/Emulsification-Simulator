# EmulSim

Multi-scale simulation of emulsification for macroporous double-network
polysaccharide hydrogel microsphere preparation.

EmulSim models the complete fabrication pipeline: from rotor-stator
emulsification of an agarose-chitosan aqueous phase in hot oil, through
thermally-induced gelation and spinodal decomposition, genipin-mediated
crosslinking of the chitosan network, to final double-network mechanical
properties. The target application is producing hydrogel microspheres with
controlled pore size for cell culture scaffolds and drug delivery.

## Quick Start

### Install

```bash
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
```

### Example Output

```
=== Simulation Results ===
  L1 Emulsification:  d32 = 42.15 um   span = 1.23
  L2 Gelation:        pore = 78.3 nm    porosity = 0.412
  L3 Crosslinking:    p = 0.247         G_chit = 1850 Pa
  L4 Mechanical:      G_DN = 8420 Pa    E* = 25100 Pa
```

### RPM Sweep

```bash
python -m emulsim sweep --rpm-min 3000 --rpm-max 25000 --rpm-steps 6
```

### Bayesian Optimization

```bash
python -m emulsim optimize --n-initial 15 --max-iter 200
```

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
 | L3: Crosslinking  |   Genipin-chitosan ODE kinetics
 | Kinetics          |   Arrhenius rate law
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
full reference configuration with comments.

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
print(f"G_DN = {result.mechanical.G_DN:.0f} Pa")
```

## Project Structure

```
emulsification_sim/
  configs/
    default.toml              # Reference configuration
  data/
    properties.toml           # Material property database
  docs/
    01_scientific_advisor_report.md
    02_computational_architecture.md
    03_scientific_review.md
  src/emulsim/
    __init__.py               # Top-level API (run_pipeline, key types)
    __main__.py               # CLI entry point (python -m emulsim)
    config.py                 # TOML config loader
    datatypes.py              # Dataclasses for params, results
    pipeline/
      orchestrator.py         # L1->L2->L3->L4 sequencing
    level1_emulsification/    # PBE solver, breakage/coalescence kernels
    level2_gelation/          # Cahn-Hilliard 2D solver, pore analysis
    level3_crosslinking/      # Genipin-chitosan ODE solver
    level4_mechanical/        # DN mechanical model
    properties/               # T-dependent viscosity, IFT, thermodynamics
    optimization/             # BoTorch multi-objective optimization
    visualization/            # Plotting utilities
  tests/
  notebooks/
```

## Architecture Documentation

Detailed design documents are in `docs/`:

1. **Scientific Advisor Report** (`docs/01_scientific_advisor_report.md`) --
   Physical models, governing equations, and parameter choices.
2. **Computational Architecture** (`docs/02_computational_architecture.md`) --
   Software design, data flow, and solver implementations.
3. **Scientific Review** (`docs/03_scientific_review.md`) --
   Validation strategy, known limitations, and future directions.

## Requirements

- Python >= 3.11
- NumPy, SciPy, h5py, Matplotlib
- Optional: BoTorch + PyTorch (for optimization), Plotly + Jupyter (for interactive viz)

## License

See `LICENSE` for details.
