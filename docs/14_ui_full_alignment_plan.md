# EmulSim 3-Module UI Alignment Plan

**Date:** 2026-04-11
**Synthesised from:** Scientific Advisor, Dev-Orchestrator, Architect
**Scope:** Comprehensive UI redesign for architectural coherence across all 3 modules

---

## Executive Summary

Three specialist roles independently analyzed the UI alignment gap. Key consensus:

1. **Module 1 UI has 5 parameter defects** (DDA hardcoded, T_gel missing, phi_d uncapped,
   c_agarose/c_chitosan ranges too wide) and **4 missing validations**
2. **Modules 2 and 3 have ZERO UI** — entire backends are unwired
3. **Architecture recommendation: Single app.py dispatcher with extracted module files**
   (not multi-page Streamlit, not monolithic 2500-line app.py)
4. **~953 new LOC** across app.py modifications + 2 new plot files (~460 LOC)
5. **18 work nodes** in 6 waves, with clear dependency ordering
6. **Every numerical output needs a confidence label** (mechanistic/empirical/ranking-only)
7. **Trust must cascade monotonically** M1 → M2 → M3

---

## Architecture Decision

**Single app.py with extracted sidebar/results modules:**

```
src/emulsim/visualization/
  app.py              # ~400 LOC: config, state init, scope selector, run button, tab dispatch
  sidebar_m1.py       # ~300 LOC: existing M1 controls (extracted)
  sidebar_m2.py       # ~200 LOC: step builder, ACS preview
  sidebar_m3.py       # ~300 LOC: column, gradient, detection, catalysis
  results_m1.py       # ~250 LOC: tabs 1-5 (extracted)
  results_m2.py       # ~200 LOC: tab 6 (functionalization)
  results_m3.py       # ~250 LOC: tab 7 (chromatography/catalysis)
  plots.py            # existing, unchanged
  plots_m2.py         # NEW: ACS waterfall, surface area, modification timeline
  plots_m3.py         # NEW: chromatogram, breakthrough, M-M, eta, decay, ESI, pressure
```

---

## Sidebar Structure

```
[Pipeline Scope: M1 / M1+M2 / M1+M2+M3]

── Module 1: Fabrication ──        (always visible)
  Hardware Mode, Scientific Mode
  Emulsification, Formulation, Gelation, Crosslinking
  Material Constants, Optimization Targets

── Module 2: Functionalization ──  (visible when scope >= M1+M2)
  Step Builder (1-5 dynamic steps)
    Per step: type, reagent, concentration, T, time, pH
  ACS Budget Preview (real-time from M1 contract)

── Module 3: Performance ──        (visible when scope >= M1+M2+M3)
  Application Mode (Chromatography / Catalysis)
  Column/Reactor Geometry + pressure preview
  [Chromatography]: Sample, Gradient Designer, Isotherm, Detection
  [Catalysis]: V_max, K_m, D_eff, k_deact, substrate feed
```

---

## Results Tab Structure

```
[Dashboard | L1 | L2 | L3 | L4 | M2: Functionalization | M3: Performance]
                                    (conditional)          (conditional)
```

**Tab 6 (M2):** ACS waterfall, surface area breakdown, modification timeline,
step details table, updated G_DN/E_star, trust indicators

**Tab 7 (M3 Chromatography):** Chromatogram with gradient overlay, breakthrough
curve with DBC annotations, peak table, pressure gauge, mass balance indicator

**Tab 7 (M3 Catalysis):** Conversion vs time, Michaelis-Menten curve,
effectiveness factor diagnostic, activity decay, Thiele regime assessment

---

## 18 Work Nodes (6 Waves)

| Wave | Nodes | LOC | Key Deliverable |
|------|-------|-----|-----------------|
| 1 | U-M2.1 (scope selector), U-M1.1 (DDA+T_gel), U-M1.2 (range caps) | ~40 | Foundation controls |
| 2 | U-M2.2 (step builder), U-M3.1 (app mode selector) | ~88 | M2/M3 sidebar sections |
| 3 | U-M2.3 (ACS viz), U-M3.2 (column), U-M3.3 (sample), U-M3.9 (catalysis params) | ~220 | Parameter panels |
| 4 | U-M2.4 (M2 results), U-M3.4 (gradient), U-M3.5 (isotherm), U-M3.6 (detection) | ~180 | Complex controls |
| 5 | U-M3.7 (chrom results), U-M3.8 (catalysis results) | ~310 | Results + plots |
| 6 | U-INT.1 (pipeline wiring), U-INT.2 (validations), U-INT.3 (trust dashboard) | ~115 | Integration |
| **Total** | **18 nodes** | **~953** | |

---

## 13 New Plot Functions

**plots_m2.py (3 functions):**
- `plot_acs_waterfall()` — stacked bar of ACS states per modification step
- `plot_surface_area_comparison()` — external vs internal vs accessible
- `plot_modification_timeline()` — Gantt-style step timeline with conversions

**plots_m3.py (10 functions):**
- `plot_chromatogram()` — UV + gradient dual-axis overlay
- `plot_breakthrough_curve()` — C/C0 vs CV with DBC threshold annotations
- `plot_peak_table()` — formatted Plotly table
- `plot_gradient_preview()` — mini gradient shape
- `plot_michaelis_menten()` — intrinsic vs effective rate
- `plot_effectiveness_factor()` — eta vs Phi with regime zones
- `plot_activity_decay()` — exponential decay with half-life
- `plot_conversion_vs_time()` — substrate conversion
- `plot_esi_spectrum()` — ESI-MS charge envelope
- `plot_pressure_flow_curve()` — dP vs Q with safe zone

---

## Scientific Validity Requirements

### Confidence Labels (every output)
- **Mechanistic prediction**: calibrated parameters within validated domain
- **Empirical estimate**: fitted correlations, literature constants (+/- 20-30%)
- **Ranking only**: uncalibrated defaults, phenomenological models

### Trust Cascade
- M1 UNRELIABLE → M2 BLOCKER ("Fix M1 before proceeding")
- M1 CAUTION → M2 inherits CAUTION, displays wider uncertainty
- M2 trust = min(M1 trust, M2 gate results)
- M3 trust = min(M2 trust, M3 gate results)

### Parameter Validations (26 rules)
- Module 1: 5 range caps, 4 physical impossibility checks
- Module 2: 11 chemistry-specific validations (pH, pore exclusion, ACS conservation)
- Module 3: 11 transport/pressure/mass-balance checks

### Disclaimers per Module
- M1: "Uncalibrated results suitable for screening only"
- M2: "ACS estimates may differ from experimental titration by 2-5x"
- M3: "Performance predictions assume ideal packing and user-supplied isotherms"

---

## Pipeline Execution Flow

```
[Run Full Pipeline] →
  Phase 1: M1 pipeline (1-10s) → FullResult → assess_trust()
  Phase 2: export_for_module2() → ModificationOrchestrator.run() → FunctionalMicrosphere
  Phase 3: run_breakthrough() / run_gradient_elution() / solve_packed_bed() → Results
  Display: 7 tabs with trust indicators per module
```

**Invalidation cascade:** M1 param change → clear M2+M3 results. M2 change → clear M3.

---

## Session State Keys (15 new)

| Key | Type | Purpose |
|-----|------|---------|
| `module_scope` | str | Pipeline scope selector |
| `m1_export_contract` | M1ExportContract | M1→M2 bridge |
| `m2_steps` | list[ModificationStep] | User-built steps |
| `m2_result` | FunctionalMicrosphere | M2 output |
| `m3_app_mode` | str | Chromatography/Catalysis |
| `m3_result_bt` | BreakthroughResult | Breakthrough output |
| `m3_result_ge` | GradientElutionResult | Gradient elution output |
| `m3_result_cat` | CatalyticResult | Catalysis output |
| + trust, timing, hash keys | | |
