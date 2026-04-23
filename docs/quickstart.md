# Quickstart Guide

## Installation

```bash
git clone https://github.com/tocvicmeng-prog/Emulsification-Simulator.git
cd Emulsification-Simulator
pip install -e .
```

## Your First Simulation

### Option 1: Command Line
```bash
# Show default parameters
python -m emulsim info

# Fast smoke (~1 s) — useful for verifying the install and CI gates
python -m emulsim run configs/fast_smoke.toml --quiet

# Full default research run (~4 min, dominated by L2 phase-field)
python -m emulsim run

# Run with custom RPM
python -m emulsim run --rpm 15000
```

### Expected baseline output

`fast_smoke.toml` should reliably produce (within rounding):

```
=== Simulation Results ===
  L1 Emulsification:  d32 = 22.08 um   span = 1.04
  L2 Gelation:        pore = 180.9 nm  porosity = 0.871
  L3 Crosslinking:    p = 0.040        G_chit = 2062 Pa
  L4 Mechanical:      G_DN = 70766 Pa  E* = 257332 Pa
```

These are uncalibrated semi-quantitative outputs (see Calibration below).
The `tests/test_smoke.py` gate (`pytest -m smoke`) verifies these stay
inside generous sanity bounds across releases. If you see materially
different values, the model defaults have drifted.

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
- **9 Crosslinkers**: Genipin, Glutaraldehyde, EDC/NHS, PEGDA+UV, TPP, STMP, Epichlorohydrin, DVS, Citric Acid
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

## v7.0+ CLI commands

Beyond `run`/`sweep`/`optimize`/`uncertainty`/`info`/`ui`, v7.0.1 adds:

```bash
# Run L2-L4 across DSD quantiles (batch variability, Node 19)
python -m emulsim batch --quantiles 0.10,0.50,0.90 configs/default.toml

# Emit a ProcessDossier JSON for reproducibility (Node 16)
python -m emulsim dossier configs/default.toml --output dossier.json

# Ingest wet-lab AssayRecord JSONs into a CalibrationStore fit JSON (Node 20)
python -m emulsim ingest L1 --assay-dir data/validation/l1_dsd/assays \
    --output data/validation/l1_dsd/fits/fit.json

# Default uncertainty now uses the unified engine + parallel MC
python -m emulsim uncertainty --n-samples 50 --n-jobs 4
```

The default `uncertainty` engine is now `UnifiedUncertaintyEngine`
(consistent schema, calibration-store posterior absorption); pass
`--engine legacy` for v6.x byte-equivalent output.

## Calibration

See `docs/04_calibration_protocol.md` for a 5-study wet-lab protocol to calibrate the simulation constants against your specific materials.

> **Note on Node 8 (v6.1):** The L2 empirical pore-size formula is
> independent of `alpha_final`. Node 8 wired `solve_gelation_empirical`
> to receive the actual Avrami output via `timing=` and reflect it in
> `model_manifest.diagnostics["alpha_final_from_timing"]`, but
> the pore-size prediction itself depends only on concentration +
> cooling rate + bead radius. Users will not see different pore numbers
> after Node 8 — the change is honest metadata reporting, not a physics
> update.

Once you have calibration data, supply it via `RunContext` (Node 7):

```python
from emulsim.calibration.calibration_store import CalibrationStore
from emulsim.datatypes import RunContext
from emulsim.pipeline.orchestrator import PipelineOrchestrator

store = CalibrationStore()
store.load_json("my_calibration.json")
ctx = RunContext(calibration_store=store)
result = PipelineOrchestrator().run_single(params, run_context=ctx)
# result.run_report.diagnostics["calibrations_applied"] lists every override.
```

Calibration entries with `target_module="L1"` rewrite `KernelConfig` constants
(breakage_C1/C2/C3, coalescence_C4/C5); `"L2"`, `"L3"`, `"L4"` rewrite
`MaterialProperties` fields. Apply order is documented in
`pipeline/orchestrator.py`.

## Runtime expectations

| Config | Approximate runtime | When to use |
|---|---|---|
| `configs/fast_smoke.toml` | ~0.2 s | CI smoke gate, first install verification |
| `configs/default.toml` | ~4 minutes (n_grid=128 phase field) | Production research run |
| `configs/stirred_vessel.toml` | varies | Stirred-vessel mode comparison |

Set `[solver.level2].n_grid` lower (e.g. 32) to make the default config
finish in under a minute at the cost of pore-morphology resolution.

## Evidence tiers (v6.1)

Every `FullResult` now carries a `RunReport` with the weakest evidence
tier across L1-L4 + M2 + M3:

- `VALIDATED_QUANTITATIVE` — calibrated against your specific system
- `CALIBRATED_LOCAL` — calibrated against an analogous system
- `SEMI_QUANTITATIVE` — empirical model, not locally calibrated (default)
- `QUALITATIVE_TREND` — directional only (e.g. EDC/NHS approximate fallback)
- `UNSUPPORTED` — model not applicable to this chemistry/regime

The Bayesian optimizer (`python -m emulsim optimize`) excludes
`QUALITATIVE_TREND` and `UNSUPPORTED` candidates from the Pareto front by
default; each surviving Pareto point is labelled with its weakest tier in
`output/optimization/optimization_results.json`.
