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

### M1 — Fabrication (since v9.0 "Family-First")

- **Four polymer platforms**, each wired end-to-end with its own
  formulation page and physics model:
  - **Agarose + Chitosan** — thermal TIPS + optional covalent
    crosslinking (the original platform)
  - **Alginate** — ionotropic Ca²⁺ gelation, external CaCl₂ bath OR
    internal release (GDL + CaCO₃)
  - **Cellulose-NIPS** — non-solvent-induced phase separation
  - **PLGA** — solvent evaporation
- **Hardware Mode** inside M1 Emulsification: stirred-vessel (Stirrer A
  pitched-blade or Stirrer B rotor-stator) and rotor-stator legacy
- **8 crosslinkers** (Genipin, Glutaraldehyde, EDC/NHS, PEGDA+UV, TPP,
  Epichlorohydrin, DVS, Citric Acid) with 4 kinetics models
- **6 surfactants** (Span-80/60/40/85, PGPR, Lecithin) with
  Szyszkowski-Langmuir IFT model

### M2 — Functionalization

- **Surface chemistry workflows**: amine secondary crosslinking,
  hydroxyl activation (ECH/DVS), ligand and protein coupling, metal
  charging (IMAC), spacer-arm chemistry, washing, quenching
- **44 wet-lab functionalization protocols** shipped as a PDF
  (Appendix J) with full SDS-lite safety blocks
- Mechanistic two-step EDC/NHS kinetics with mass-balance enforcement

### M3 — Performance

- **Chromatography**: Lumped Rate Model (LRM) breakthrough simulation
  with Langmuir, HIC, and IMAC isotherms; gradient-aware elution
- **Catalysis**: packed-bed enzyme reactor with Michaelis-Menten kinetics,
  axial dispersion, internal effectiveness factor, enzyme deactivation

### Cross-cutting

- **Streamlit Web UI** with the Family-First M1 + M2 + M3 tabs and
  Plotly dashboards
- **Bayesian optimization** (BoTorch) for multi-objective parameter
  tuning, plus inverse design and CVaR-robust BO
- **Monte Carlo uncertainty propagation** for confidence intervals
- **Digital twin** (EnKF replay) and **lifetime model** for activity
  decay over operational cycles
- **One-click Windows installer** (Inno Setup, GPL-3.0 EULA shown on
  first page); see `installer/README.md`
- **Calibration protocol** — 5-study wet-lab protocol to calibrate
  simulation constants (`docs/04_calibration_protocol.md`)

## Quick Start

### Install

```bash
git clone https://github.com/tocvicmeng-prog/Emulsification-Simulator.git
cd Emulsification-Simulator
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
Emulsification-Simulator/
  configs/                    # TOML configs (default, fast_smoke, stirred_vessel, ...)
  data/
    properties.toml           # Material property database
    validation/               # Validation scaffold (l1_dsd, l2_pore, l3_kinetics, ...)
  docs/
    quickstart.md             # First-run guide
    configuration.md          # Full TOML parameter reference
    decisions/                # ADRs: Python policy, optimization stack pin, ...
    01_scientific_advisor_report.md
    02_computational_architecture.md
    INDEX.md                  # navigation map for docs/
    04_calibration_protocol.md
    INDEX.md                  # navigation map for docs/
    module2_history.md        # M2 architecture history + shipped audit findings
    ui_evolution.md           # Family-First UI history + shipped backend fixes
    19_ligand_protein_coupling_candidates.md
    20_linker_arm_candidates.md
    f1a_*, f1b_*, f1c_*, f2_*, f4b_*, f5_*  # per-platform protocols
    user_manual/              # First Edition + Appendix J PDFs (shipped in installer)
  installer/                  # Inno Setup installer sources (build_installer.bat, EmulSim.iss, ...)
  src/emulsim/
    __init__.py               # Top-level API (run_pipeline, key types)
    __main__.py               # CLI entry point (python -m emulsim)
    config.py                 # TOML config loader
    datatypes.py              # Dataclasses for params, results, manifests
    pipeline/                 # L1->L2->L3->L4 sequencing
    properties/               # T-dependent viscosity, IFT, thermodynamics
    level1_emulsification/    # PBE solver, breakage/coalescence kernels (numba JIT)
    level2_gelation/          # Cahn-Hilliard 2D, pore analysis, ionic Ca, NIPS, solvent evap
    level3_crosslinking/      # Multi-mechanism crosslinking kinetics
    level4_mechanical/        # Double-network mechanical model
    module2_functionalization/# Surface chemistry (EDC/NHS, ACS, modification steps, reagents)
    module3_performance/      # Chromatography (LRM, breakthrough, gradient elution),
                              # catalysis (packed-bed PFR), isotherms (Langmuir, HIC, IMAC),
                              # detection (UV, MS), hydrodynamics
    optimization/             # BoTorch multi-objective + inverse design + robust BO
    digital_twin/             # EnKF replay, schema
    lifetime/                 # Operational lifetime model
    calibration/              # CalibrationStore + run-context injection
    protocols/                # Mechanism data + protocol generator (Appendix J)
    visualization/            # Streamlit app, tabs (M1/M2/M3), plots, panels
  tests/                      # ~900 tests; smoke + slow markers
  release/                    # Installer .exe outputs (gitignored)
  CHANGELOG.md
  CLAUDE.md                   # Project instructions for Claude Code
  DESIGN.md                   # Visual design system (read before UI changes)
```

## Documentation

### Start here

- **Quickstart** (`docs/quickstart.md`) — your first simulation in five minutes
- **Configuration reference** (`docs/configuration.md`) — every TOML field with units and physical meaning
- **CHANGELOG.md** — what shipped in each release, in user-facing language

### Foundational scientific docs

- **Scientific Advisor Report** (`docs/01_scientific_advisor_report.md`) — physical models, governing equations, plus **Appendix A** consolidating audit findings (§A.1, §A.4, §A.5), design-principle rationale (§A.2), Cluster F platform science (§A.3), the EDC/NHS mechanism + rate constants (§A.6), and the full crosslinker library provenance (§A.7)
- **Computational Architecture** (`docs/02_computational_architecture.md`) — software design, data flow, solver implementations
- **Calibration Protocol** (`docs/04_calibration_protocol.md`) — 5-study wet-lab calibration plan
- **Documentation Index** (`docs/INDEX.md`) — navigation map for every file in `docs/`

### Architecture decisions (ADRs)

- `docs/decisions/ADR-001-python-version-policy.md` — why Python is pinned to `>=3.11,<3.13`
- `docs/decisions/ADR-002-optimization-stack-pin.md` — why botorch / gpytorch / torch are version-pinned

### Module / family deep dives

- `docs/f1a_alginate_protocol.md` — alginate ionic-Ca gelation
- `docs/f1b_cellulose_nips_protocol.md` — cellulose NIPS
- `docs/f1c_plga_protocol.md` — PLGA solvent evaporation
- `docs/f2_digital_twin_protocol.md` — digital twin / EnKF replay
- `docs/f4b_cvar_protocol.md` — CVaR-robust Bayesian optimization
- `docs/f5_md_ingest_protocol.md` — MARTINI MD parameter ingest
- `docs/module2_history.md` — Module 2 (functionalization) architecture, current reagent library, audit findings shipped
- `docs/ui_evolution.md` — Family-First UI architecture, v9.0 redesign backend-fix findings
- `docs/19_ligand_protein_coupling_candidates.md` — scientific candidate screen for M2 ligand/protein additions (CAS, MW, kinetic parameters, suitability tiers)
- `docs/20_linker_arm_candidates.md` — scientific candidate screen for M2 spacer arms (CAS, length, distal groups, compatibility matrix)

### History

The repository previously kept a full `docs/05_*` through `docs/35_*` set of
planning, audit, and remediation documents. Those were consolidated in the
2026-04-24 content audit: version-specific planning docs for superseded
releases (v5.8 / v5.9 / v6.0 / v7.x / v8.x) were removed; load-bearing
scientific content was preserved as Appendix A of
`docs/01_scientific_advisor_report.md`; M2 and UI iteration trails were
consolidated into `docs/module2_history.md` and `docs/ui_evolution.md`.
The pre-audit snapshot lives at tag `v9.2.2-pre-docs-audit`. For
user-facing "what shipped in each release", read `CHANGELOG.md`.

### Installer

- **`installer/README.md`** — how the Windows `.exe` is built and what it
  installs. End users don't read this; it's for maintainers cutting a release.

## Requirements

- Python `>=3.11,<3.13` (pinned in `pyproject.toml`; ADR-001 explains why
  newer Python versions are excluded for now)
- NumPy, SciPy, h5py, Matplotlib, Pydantic
- Optional: `[ui]` for Streamlit + Plotly, `[optimization]` for the
  pinned BoTorch / GpyTorch / PyTorch stack (see ADR-002), `[dev]` for
  pytest + numba + pytest-timeout

## License

EmulSim is distributed under the **GNU General Public License v3.0
(GPL-3.0)**. The full licence text is in `LICENSE` at the repo root and
is also available at <https://www.gnu.org/licenses/gpl-3.0.en.html>.

The intellectual property in this software, including the source code,
documentation, and accompanying assets, belongs to **Holocyte Pty Ltd**.

The Windows installer shows the same statement on its first page; see
`installer/LICENSE_AND_IP.txt` for the user-facing wording.
