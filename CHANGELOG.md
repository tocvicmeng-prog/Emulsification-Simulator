# Changelog

## v7.0.1 (2026-04-17) — Audit remediation patch

Closes 8 of 10 findings from the post-Nodes-1-20 full-system audit. P0
ship-blockers fixed; v7.0 features now reachable from the CLI.

### P0 fixes (release blockers)
- **N1 (HIGH)** — `pipeline/orchestrator.py` no longer mutates the caller's
  `params.emulsification.kernels` in place when applying L1 calibration.
  Callers that reuse a `SimulationParameters` instance across multiple
  `run_single` calls (e.g. `batch_variability.run_batch`, parameter
  sweeps, optimisation campaigns) no longer see calibrated kernels leak
  between iterations. Regression test in `test_run_context.py`.
- **N2 (HIGH)** — `UnifiedUncertaintyEngine.run_m1l4` no longer claims to
  have sampled `CALIBRATION_POSTERIOR` when it has only absorbed the
  posterior into the spec. The new
  `UnifiedUncertaintyResult.kinds_declared_but_not_sampled` field
  records the v7.0 limitation honestly.

### P1 fixes (CLI surface — closes audit N4 + N5)
- **`python -m emulsim batch`** — surface
  `pipeline.batch_variability.run_batch` on the CLI. Pass `--quantiles`
  and `--output`; prints mass-weighted mean / per-quantile percentile
  table.
- **`python -m emulsim dossier`** — run the pipeline and emit a
  `ProcessDossier` JSON artifact for reproducibility. Records inputs,
  result summary, manifests, calibrations, environment.
- **`python -m emulsim ingest L1`** — ingest a directory of
  `AssayRecord` JSON files, run the L1 fitter, write a
  `CalibrationStore`-loadable fit JSON. v7.1 will add L2/L3/L4/M2.
- **`python -m emulsim uncertainty`** now defaults to the
  `UnifiedUncertaintyEngine` (Node 18) and exposes `--n-jobs` for
  Node 15's parallel MC. Pass `--engine legacy` for v6.x byte-equivalent
  output.
- **N3 follow-up** — `QuantileRun.representative_diameter_m` property
  added so downstream consumers don't accidentally read
  `full_result.emulsification.d50` (which is shared by reference across
  all per-quantile runs and reflects the BASE L1 DSD).

### P2 polish
- **N7** — `UncertaintyPropagator.run` auto-falls-back to serial when
  `n_samples < 4 × |n_jobs|`. Joblib startup + Numba JIT cold-compile
  dominate below this threshold.
- **N8** — `run_batch` silently sort+dedupes the `quantiles` argument.
  Duplicate or unsorted input no longer produces ill-defined mass
  fractions.

### P3 documentation
- **N6** — `INSTALL.md` documents the Numba JIT cache location and the
  `NUMBA_CACHE_DIR` environment-variable workaround for read-only
  Python installs (corporate, conda `--no-write-pkgs`,
  `pip install --user` on network shares).
- **N9** — Documenting that Node 8's L2 timing wiring was a metadata
  fix only; the empirical pore-size formula remains independent of
  `alpha_final`. The `model_manifest.diagnostics.alpha_final_from_timing`
  field now reflects the actual Avrami output instead of a hardcoded
  0.999, but pore predictions at typical conditions are unchanged.

### Tests
- 25 new tests across the patch (Nodes 22-29). 0 regressions.

---

## v7.0 (2026-04-17) — Engineering portion (Nodes 14-20)

Closes engineering items from the consensus v7.0 plan (doc 34 §9). F1
closure (kernel re-fit) remains gated on Study A wet-lab data.

### New modules
- `process_dossier.py` — `ProcessDossier` aggregator + JSON export
- `assay_record.py` — `AssayRecord` public data model with 12 `AssayKind` values
- `uncertainty_unified.py` — `UnifiedUncertaintyEngine` single entrypoint
- `pipeline/batch_variability.py` — `run_batch` over DSD quantiles
- `calibration/fitters.py` — stub L1 DSD fitter

### Performance
- Numba JIT for `breakage_rate_alopaeus`, `breakage_rate_coulaloglou`,
  `coalescence_rate_ct` matrix builder (5-10× on coalescence; matches
  NumPy to 1e-12 rtol).
- joblib parallel MC via `UncertaintyPropagator(n_jobs=-1)`.

### Calibration data scaffold
- `data/validation/{l1_dsd,l2_pore,l3_kinetics,l4_mechanics,m2_capacity}/`
  directory tree with JSON-Schema for L1 DSD assays.

---

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
