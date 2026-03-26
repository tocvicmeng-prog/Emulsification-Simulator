# Changelog

## v0.1.0 (2026-03-26) — Initial Release

### Simulation Pipeline
- 4-level sequential pipeline: PBE emulsification → empirical gelation → multi-mechanism crosslinking → IPN mechanical properties
- 8 crosslinkers with 4 kinetics models (second-order amine/hydroxyl, UV dose, ionic instant)
- 6 surfactants with Szyszkowski-Langmuir IFT model
- Empirical pore-size model calibrated to literature (Pernodet 1997, Chen 2017)
- 2D Cahn-Hilliard phase-field solver available as advanced option

### Web UI (Streamlit)
- Interactive parameter input with sliders and dropdowns
- Reagent selection (crosslinker + surfactant) with per-reagent defaults
- Per-constant Literature/Custom toggle with calibration protocol links
- Results dashboard with Plotly charts (size distribution, phase field, kinetics, Hertz, Kav)
- Trust assessment with 10 automated reliability checks
- Optimization assessment with actionable recommendations

### CLI
- `python -m emulsim run` — full pipeline
- `python -m emulsim sweep` — RPM parameter sweep
- `python -m emulsim optimize` — BoTorch Bayesian optimization
- `python -m emulsim uncertainty` — Monte Carlo uncertainty propagation
- `python -m emulsim ui` — launch Streamlit web interface
- `python -m emulsim info` — display parameters and properties

### Documentation
- Scientific advisory report (docs/01)
- Computational architecture (docs/02)
- Scientific review with formula verification (docs/03)
- Calibration wet-lab protocol — 5 studies, 1081 lines (docs/04)
- Literature constants database with sources and DOIs
- Reagent library with 8 crosslinkers and 6 surfactants

### Quality Assurance
- 9 rounds of Codex (OpenAI) adversarial review — 63+ findings, all addressed
- Scientific Advisor review — 4 critical bugs fixed
- Dev-Orchestrator usability review — all priorities implemented
- Input validation, trust gates, uncertainty propagation
- 107+ unit tests
