# EmulSim UI Final Implementation Plan

**Date:** 2026-04-11
**Input:** Original plan (doc 14) + 3rd-party audit (doc 15) + 3-role synthesis
**Status:** Implementation-ready — backend fixes verified, all 7 findings confirmed

---

## All 7 Audit Findings Confirmed by Code Inspection

| # | Finding | Severity | Backend Issue? | Fix Timing |
|---|---------|----------|---------------|------------|
| F1 | UI is M1-only, M2/M3 backends unwired | Critical | Integration gap | DURING UI work |
| F2 | No per-output confidence labels | Critical | Missing OutputMetadata schema | BEFORE (schema) |
| F3 | M2 UI implies broad chemistry coverage | Critical | 3 of 5 step types unimplemented | BEFORE (registry) |
| F4 | pH/temperature don't affect M2 rates | High | **YES: k0=0.0 bypasses Arrhenius** | BEFORE (~15 LOC) |
| F5 | Gradient not coupled to competitive Langmuir | High | **YES: gradient value discarded** | BEFORE (label+phase switch) |
| F6 | Mass balance >2% in some tests | High | Computation correct, exposure missing | BEFORE (quality enum) |
| F7 | M1 advice can recommend unsafe actions | High | DDA hardcoded, heuristic advice | BEFORE (DDA input) |

---

## Backend Fixes Required BEFORE UI Work

### BF-1: M2 Arrhenius Temperature Fix (F4)

**Problem:** `modification_steps.py` passes `k0=0.0` to `solve_second_order_consumption()`, 
bypassing Arrhenius. Temperature slider is inert.

**Fix:** Back-calculate `k0 = k_forward * exp(E_a/(R*T_ref))` from reagent profile, pass real k0.
```python
def _arrhenius_prefactor(k_ref, E_a, T_ref):
    return k_ref * math.exp(E_a / (8.314 * T_ref)) if E_a > 0 else 0.0
```
Apply in both `_solve_crosslinking_step` and `_solve_activation_step`.

**LOC:** ~15 | **Files:** modification_steps.py | **Tier:** Sonnet

### BF-2: Gradient Feed Phase Switching (F5)

**Problem:** `run_gradient_elution()` passes `C_feed.copy()` as inlet concentration at ALL times.
Gradient value is computed then discarded (`_ = gradient.value_at_time(t)`).

**Fix:** (a) Add load/elute phase switching — zero protein feed during gradient ramp.
(b) Add `gradient_affects_binding: bool` to `GradientElutionResult`.
(c) Add `gradient_sensitive` property to isotherm interface.

**LOC:** ~30 | **Files:** orchestrator.py, competitive_langmuir.py | **Tier:** Sonnet

### BF-3: Mass Balance Quality Enum (F6)

**Problem:** Mass balance error computed correctly but only emitted as `warnings.warn()`.

**Fix:** Add `MassBalanceQuality` enum (ACCEPTABLE/CAUTION/UNRELIABLE) and 
`classify_mass_balance()` function. Add to `LRMResult`.

**LOC:** ~20 | **Files:** lumped_rate.py | **Tier:** Sonnet

**Total backend fixes:** ~65 LOC, all Sonnet tier

---

## Final Build Order (7 Phases, 41 Work Nodes)

```
Pre-Phase: Backend Fixes (BF-1, BF-2, BF-3)         ~65 LOC
    ↓
Phase 0: UI Contract Layer (metadata, validators,     ~300 LOC
         state management, units)
    ↓
Phase 1: Module 1 Repairs (DDA, T_gel, phi_d caps,   ~195 LOC
         chemistry-specific crosslinker, confidence)
    ↓
Phase 2: Module 2 Minimal UI (4 workflows only,       ~440 LOC
         disabled LIGAND/PROTEIN/QUENCH)
    ↓
Phase 3: Module 3 Breakthrough UI (single-component,  ~350 LOC
         mass balance gates, pressure)
    ↓
Phase 4: Multi-Component + Gradient (SMA-coupled,     ~440 LOC
         gradient badge, peak metrics)
    ↓
Phase 5: Catalysis UI (gated on test reliability)     ~330 LOC
    ↓
Phase 6: Integration + Trust Dashboard                 ~240 LOC
```

**Total: ~2,360 LOC across 41 work nodes**

---

## Key Restrictions (from audit)

| Item | Action |
|------|--------|
| LIGAND_COUPLING, PROTEIN_COUPLING, QUENCHING | Show as "Planned" with disabled toggle |
| pH/temperature in M2 | Label as "validation metadata only" until Arrhenius fix deployed |
| Gradient overlay (competitive Langmuir) | Badge: "Gradient affects binding: NO" |
| Catalysis UI | Gated on test suite passing in <60s |
| Default isotherm parameters | Labeled "illustrative — user calibration required" |
| All numerical outputs | Must carry OutputMetadata with confidence label |

---

## Phase 0: UI Contract Layer (NEW — from audit)

| File | Purpose | LOC |
|------|---------|-----|
| `ui_model_metadata.py` | OutputMetadata dataclass, ModelBasis/ConfidenceLevel enums, pre-built metadata per output | ~80 |
| `ui_validators.py` | 26 validation rules (9 M1 + 11 M2 + 6 M3), hard blockers + soft warnings | ~120 |
| `ui_state.py` | SessionStateManager with M1→M2→M3 invalidation cascade, hash-based change detection | ~60 |
| `ui_units.py` | Unit conversion helpers (mol/m³→mM, Pa→bar, m→μm) | ~40 |

---

## Scientific Framing (all 3 roles agree)

Every output carries one of:
- **Mechanistic prediction** — validated governing equations + calibrated params
- **Empirical estimate** — fitted correlation within calibration range (+/-20-30%)
- **Ranking only** — uncalibrated defaults, relative ordering only
- **Not predicted** — placeholder, model not implemented

Trust cascades monotonically: M1 UNRELIABLE → M2 BLOCKER → M3 cannot run.
