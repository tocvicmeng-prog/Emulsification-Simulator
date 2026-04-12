# Changelog

## v6.0 (2026-04-12) — Calibration-Enabled Process Simulation

Transitions EmulSim from semi-quantitative chemistry simulator to calibration-enabled process simulation platform. All uncalibrated outputs remain semi-quantitative; calibrated outputs reflect user-supplied measurements.

### UI Restructure
- Split monolithic `app.py` (1480 lines) into modular tab architecture (7 UI files, orchestrator < 210 lines)
- `tabs/tab_m1.py`: M1 Fabrication tab (inputs, run, results, optimization, trust)
- `tabs/tab_m2.py`: M2 Functionalization tab (9 step types, 52 reagent profiles)
- `tabs/tab_m3.py`: M3 Performance tab (chromatography + catalysis)
- Sidebar panels for calibration, uncertainty, and lifetime frameworks

### Gradient-Aware LRM (H6)
- `solve_lrm()` accepts time-varying `ProcessState` via `gradient_program` + `equilibrium_adapter`
- Gradient values now mechanistically affect equilibrium during LRM time integration
- `run_gradient_elution()` auto-creates adapter for gradient-sensitive isotherms
- `gradient_sensitive` + `gradient_field` properties on SMA, HIC, IMAC, ProteinA, CompetitiveAffinity isotherms
- Fully backward compatible: existing callers unchanged

### Calibration Framework (v6.0-alpha)
- `CalibrationEntry` typed dataclass with units, target, validity domain (audit F2)
- `CalibrationStore` with JSON import/export, query, and `apply_to_fmc()` (audit F13)
- UI panel: JSON upload, manual entry, color-coded confidence display

### Uncertainty Propagation (v6.0-alpha)
- `M1UncertaintyContract` with 5 CVs and two tiers: measured (Tier 1) vs assumed (Tier 2, audit F4)
- `run_with_uncertainty()` Monte Carlo through M2 pipeline producing p5/p95 bounds on q_max
- UI panel: CV sliders, tier selection, sample count configuration

### Lifetime Projection (v6.0-rc)
- `LifetimeProjection` empirical first-order deactivation model (audit F6)
- `project_lifetime()` with cycles-to-80%/50% milestones
- UI panel: interactive Plotly decay curve, empirical confidence warning

### ProcessState (v6.0-beta)
- Typed `ProcessState` dataclass replacing loose dict for process conditions
- Carries salt, pH, imidazole, sugar competitor, temperature for multi-parameter isotherms
- `EquilibriumAdapter` dispatches by isotherm class name with ProcessState routing

### New Isotherms
- `HICIsotherm`: Salt-modulated Langmuir (K_eff = K_0 * exp(m * C_salt)), requires user calibration
- `CompetitiveAffinityIsotherm`: Generalized competitive binding for lectin elution (Con A, WGA)

### Quality
- 14/14 acceptance criteria from audit Section 7 verified passing
- 24 new integration tests (12 gradient LRM + 12 v6.0 end-to-end)
- All existing v5.9 workflows pass regression (280+ total tests, 0 failures)

---

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
