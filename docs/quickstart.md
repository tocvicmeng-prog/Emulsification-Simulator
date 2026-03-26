# Quickstart Guide

## Installation

```bash
git clone https://github.com/tocvicmeng-prog/emulsification-similation-for-microspheres-preparation.git
cd emulsification-similation-for-microspheres-preparation
pip install -e .
```

## Your First Simulation

### Option 1: Command Line
```bash
# Show default parameters
python -m emulsim info

# Run with defaults (Genipin + Span-80, 10000 RPM)
python -m emulsim run

# Run with custom RPM
python -m emulsim run --rpm 15000
```

### Option 2: Web UI
```bash
python -m emulsim ui
# Open http://localhost:8501
```

### Option 3: Python API
```python
from emulsim import run_pipeline
result = run_pipeline()
print(f"d32 = {result.emulsification.d32*1e6:.1f} µm")
print(f"pore = {result.gelation.pore_size_mean*1e9:.0f} nm")
print(f"G_DN = {result.mechanical.G_DN/1000:.1f} kPa")
```

## Editing Parameters

Edit `configs/default.toml` or pass values via the CLI/UI:

```toml
[emulsification]
rpm = 15000
t_emulsification = 60.0

[formulation]
c_agarose = 42.0      # kg/m³ (4.2% w/v)
c_chitosan = 18.0     # kg/m³ (1.8% w/v)
c_span80 = 20.0       # kg/m³ (2.0% w/v)
phi_d = 0.05           # dispersed phase volume fraction

[solver.level1]
n_bins = 20
```

## Choosing Reagents

The UI provides dropdown selectors for:
- **8 Crosslinkers**: Genipin, Glutaraldehyde, EDC/NHS, PEGDA+UV, TPP, Epichlorohydrin, DVS, Citric Acid
- **6 Surfactants**: Span-80, Span-60, Span-40, Span-85, PGPR, Lecithin

Each reagent has literature-sourced kinetic parameters that automatically feed into the simulation.

## Interpreting Results

| Output | What it means | Target range |
|--------|--------------|--------------|
| d32 | Sauter mean droplet diameter | 2-10 µm for chromatography beads |
| Pore size | Mean macropore diameter | 60-100 nm for protein SEC |
| G_DN | Double-network shear modulus | >10 kPa for column packing |
| Span | Size distribution width | <2.0 for uniformity |

## Running Optimization

```bash
python -m emulsim optimize --n-initial 15 --max-iter 100
```

## Uncertainty Quantification

```bash
python -m emulsim uncertainty --n-samples 20
```

## Calibration

See `docs/04_calibration_protocol.md` for a 5-study wet-lab protocol to calibrate the simulation constants against your specific materials.
