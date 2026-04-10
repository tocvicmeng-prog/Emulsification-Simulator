# EmulSim UI Update Plan — Unified Synthesis

**Date:** 2026-04-11
**Synthesised from:** Scientific Advisor, Dev-Orchestrator, Architect
**Current UI:** `src/emulsim/visualization/app.py` (700-line monolith) + `plots.py` (166 lines)

---

## Executive Summary

Three specialist roles independently analyzed the Streamlit UI against the new M1-M4
backend capabilities and physical-world consistency. Key findings:

1. **ModelMode selector is completely missing** — the UI always defaults to `HYBRID_COUPLED`,
   making the Flory-Rehner affine model and mode-aware trust gates unreachable
2. **No real-time stoichiometry guidance** — users only learn about crosslinker-limited
   formulations after running the simulation
3. **Hashin-Shtrikman bounds are computed but never displayed**
4. **Several physically dangerous parameter ranges** exist (phi_d up to 0.75 in stirred
   vessel mode, chitosan up to 5%, agarose up to 10%)
5. **DDA (degree of deacetylation) is not user-adjustable** despite profoundly affecting
   crosslinking stoichiometry

---

## 1. Critical Physical Corrections (HIGH PRIORITY)

### 1.1 Parameter Range Corrections

| Parameter | Current Range | Recommended Range | Issue |
|-----------|--------------|-------------------|-------|
| phi_d (stirred vessel) | 0.10–0.75 via volumes | Cap at 0.55 | Phase inversion above ~0.60 (Salager 2000) |
| c_agarose | 1.0–10.0% | 1.0–8.0% | >8% exceeds empirical model calibration, nearly unpumpable |
| c_chitosan | 0.5–5.0% | 0.5–3.0% | >3% at acidic pH produces gel-like viscosity |
| T_crosslink | 0–120C | Crosslinker-dependent | 0C: negligible kinetics; >85C: gel re-melts |

### 1.2 Missing Physical Validations (Real-Time in Sidebar)

| Check | Condition | Display |
|-------|-----------|---------|
| Stoichiometric ratio | c_xlink / [NH2] | Traffic light: red <0.05, yellow 0.05-0.15, green >0.15 |
| T_crosslink < T_gel | Must hold | Error if violated (gel re-melts during crosslinking) |
| Total polymer > 10% | c_agarose + c_chitosan | Warning: extremely viscous |
| phi_d > 0.50 | Stirred vessel | Warning: phase inversion risk |

### 1.3 Missing User-Adjustable Parameters

| Parameter | Default | Range | Impact |
|-----------|---------|-------|--------|
| **DDA** (degree of deacetylation) | 0.90 | 0.70–0.95 | Directly controls [NH2], crosslinking capacity |
| **T_gel** (gelation temperature) | 38C | 25–42C | Low-melt vs standard agarose |

---

## 2. New Feature Integration (9 Work Nodes)

### Build Order and Parallelism

```
Wave 1 (independent, can parallelize):
  U8  Trust fix: pass l2_mode to assess_trust()     [Haiku, ~5 LOC]
  U1  Model Mode selector (sidebar radio)            [Haiku, ~15 LOC]
  U2  Adaptive PBE controls (collapsible expander)   [Haiku, ~20 LOC]
  U3  Adaptive PBE results display                   [Haiku, ~10 LOC]
  U4  Stoichiometry guidance widget                   [Sonnet, ~30 LOC]

Wave 2 (after U1 lands):
  U5  Per-chemistry eta display + override            [Sonnet, ~25 LOC]
  U6  Hashin-Shtrikman bounds visualization           [Sonnet, ~35 LOC]
  U7  Model used label + FR status                    [Haiku, ~15 LOC]

Wave 3 (polish):
  U9  Help text and mode documentation                [Haiku, ~15 LOC]
```

**Total: ~170 LOC modifications to existing files. No new files needed.**

---

### Node Specifications

#### U1: Model Mode Selector

**Placement:** Top of sidebar, after Hardware Mode radio
**Widget:** `st.sidebar.radio("Scientific Mode", [...])` mapping to `ModelMode` enum
**Mode-dependent visibility:**

| UI Element | Empirical | Hybrid | Mechanistic |
|------------|-----------|--------|-------------|
| Adaptive PBE controls | Hidden | Collapsed | Expanded |
| L2 CH-2D option | Available | Available | Pre-selected |
| Flory-Rehner viz | Hidden | Hidden | Visible |
| HS bounds detail | Caption | Error bars | Full chart |
| Trust W3 severity | Suppressed | Warning | Blocker |

**Files:** `app.py` (insert after line 76, modify param construction at lines 334/364)

#### U2: Adaptive PBE Controls

**Placement:** Collapsible expander at bottom of L1 section
**Widgets:** l1_t_max slider (60-1800s), l1_conv_tol slider (0.001-0.10), l1_max_extensions (0-5)
**Files:** `app.py` (insert after line 117)

#### U3: Adaptive PBE Results Display

**Placement:** L1 results tab
**Content:** `n_extensions`, `t_converged`, convergence status indicator
**Files:** `app.py` (modify lines 536-545)

#### U4: Stoichiometry Guidance Widget

**Placement:** Below crosslinker concentration slider in sidebar
**Logic:**
- Compute [NH2] from c_chitosan, DDA, M_GlcN (or [OH] for hydroxyl crosslinkers)
- Display ratio with color-coded traffic light
- Show recommended minimum from `recommended_crosslinker_concentration()`
- Chemistry-aware: routes to NH2 for amine, OH for hydroxyl crosslinkers
**Files:** `app.py` (insert after line 189)

#### U5: Per-Chemistry Eta Display

**Placement:** Replace global eta_coupling material constant with crosslinker-specific display
**Logic:** When crosslinker is selected, show its `eta_coupling_recommended` value
- Default: locked to library value (display-only)
- Advanced: allow override via Custom toggle
**Files:** `app.py` (modify material constants section)

#### U6: Hashin-Shtrikman Bounds Visualization

**Placement:** L4 tab + KPI metric card
**Components:**
- Error bars on G_DN metric card: `76.8 kPa [17.2 – 23.6]`
- New `plot_modulus_comparison()` bar chart in plots.py
- Shaded confidence band on Hertz contact plot
**Files:** `app.py` (modify L4 tab), `plots.py` (add 2 new functions)

#### U7: Model Used Label + Flory-Rehner Status

**Placement:** L4 tab header
**Content:** `st.caption(f"Model: {m.model_used}")` with tooltip
- In mechanistic mode with FR: show swelling equilibrium + network breakdown
**Files:** `app.py` (modify L4 tab), `plots.py` (add 2 new plot functions)

#### U8: Trust Fix

**What:** Pass `l2_mode` to `assess_trust()` call (currently missing)
**Files:** `app.py` (1-line fix at line 459)

#### U9: Help Text

**What:** Mode explanations, parameter tooltips, calibration cross-references
**Files:** `app.py` (static strings)

---

## 3. Progressive Disclosure Design (All 3 Roles Agree)

Map `ModelMode` to user experience complexity:

| Level | Mode | Sidebar | Results | Target User |
|-------|------|---------|---------|-------------|
| **Beginner** | Empirical Engineering | Essential params only, no advanced sections | KPI cards + dashboard, no HS detail | Formulation screening |
| **Standard** | Hybrid Coupled | All params + collapsed advanced | Full tabs, HS error bars, trust dashboard | Routine simulation |
| **Expert** | Mechanistic Research | Expanded advanced, chi/rho exposed, CH-2D default | FR visualization, model comparison, full trust | Research validation |

---

## 4. New Plot Functions (plots.py)

| Function | Content | When Shown |
|----------|---------|------------|
| `plot_modulus_comparison(m)` | Bar chart: G_agar, G_chit, G_DN with HS error bars | Always in L4 tab |
| `plot_hertz_with_bounds(m, R)` | Hertz curve + shaded HS confidence band | Always in L4 tab |
| `plot_swelling_equilibrium(...)` | phi_eq vs nu_e parametric curve | Mechanistic mode only |
| `plot_network_contribution(...)` | Stacked bar: affine G1+G2 vs phenomenological | Mechanistic mode only |

---

## 5. Scientific Advisor Key Recommendations

### Parameters That Should NOT Be User-Adjustable

- Per-chemistry eta values (locked to reagent library; display-only)
- Flory-Rehner chi1/chi2 in non-expert modes (expose only in Mechanistic)
- Dry polymer densities (rho_agarose=1640, rho_chitosan=1400)

### Parameters That SHOULD Be Added

| Priority | Parameter | Default | Physical Basis |
|----------|-----------|---------|----------------|
| HIGH | DDA | 0.90 | Controls [NH2], chitosan charge, crosslinking capacity |
| HIGH | T_gel | 38C (311.15K) | Varies with agarose type (low-melt: 25-30C) |
| MEDIUM | chi1 (agarose-water) | 0.50 | Only in Mechanistic mode, range 0.45-0.55 |
| MEDIUM | chi2 (chitosan-water) | 0.45 | Only in Mechanistic mode, range 0.30-0.60 |
| LOW | Pore model coefficients (A, alpha, beta) | 600nm, -0.7, -0.2 | Medium confidence, calibration-dependent |

### Physically Impossible Combinations to Block

1. phi_d > 0.60 — phase inversion
2. T_crosslink > T_gel_melting (~85C) — gel re-melts
3. T_oil < 65C with standard agarose — sol is gelling
4. Total polymer > 12% w/v — exceeds pumpability
5. UV crosslinker (PEGDA) with microspheres > 100um — UV cannot penetrate

---

## 6. Results Interpretation Enhancements

### Trust Gate Dashboard (Replace Flat List)

```
┌──────────┬──────────┬──────────┬──────────┐
│  L1: ●   │  L2: ●   │  L3: ●   │  L4: ●   │
│  PBE     │  Gelation │  Xlink   │  Mech    │
│  CAUTION │  OK       │  OK      │  CAUTION │
└──────────┴──────────┴──────────┴──────────┘
```

### Mechanical Results Display

```
G_DN = 76.8 kPa  (phenomenological)
  HS bounds: [17.2, 23.6] kPa
  Model: phenomenological (switch to Mechanistic for affine IPN)

--- In Mechanistic Mode ---
| Model                | G_DN (kPa) | Status    |
|----------------------|-----------|-----------|
| Phenomenological     | 76.8      | Default   |
| Flory-Rehner Affine  | 102.4     | Converged |
| HS Lower Bound       | 17.2      | Bound     |
| HS Upper Bound       | 23.6      | Bound     |
```

---

## 7. Implementation Estimates

| Wave | Nodes | Model Tier | LOC | Time |
|------|-------|-----------|-----|------|
| Wave 1 | U1, U2, U3, U4, U8 | 4×Haiku + 1×Sonnet | ~80 | Fast |
| Wave 2 | U5, U6, U7 | 2×Sonnet + 1×Haiku | ~75 | Medium |
| Wave 3 | U9 | 1×Haiku | ~15 | Fast |
| **Total** | **9 nodes** | | **~170** | |

**All changes are modifications to 2 existing files** (`app.py`, `plots.py`). No new files.
Backend is complete — the UI just needs to wire through existing capabilities.

---

## 8. Quality Gates

| Node | Verification |
|------|-------------|
| U1 | Toggle all 3 modes → verify params.model_mode changes. Mechanistic triggers FR in L4. |
| U2 | Change adaptive settings → verify in SolverSettings. Short t → extensions trigger. |
| U3 | Run sim → L1 tab shows t_converged and n_extensions. |
| U4 | Low c_genipin → red traffic light. High c_genipin → green. Change chitosan → updates. |
| U5 | Select DVS → shows eta=+0.10. Select genipin → shows eta=-0.15. |
| U6 | Run sim → L4 shows HS bounds. G_DN metric has range text. |
| U7 | Hybrid → "phenomenological". Mechanistic → "flory_rehner_affine". |
| U8 | Mechanistic + empirical L2 → mismatch warning fires. |
| U9 | Visual inspection of help text. |

**Integration test:** Run all 3 modes × 3 crosslinkers (genipin, DVS, TPP).
Verify mode flows through sidebar → params → orchestrator → trust → results display.
