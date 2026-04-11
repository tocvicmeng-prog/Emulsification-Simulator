# EmulSim Module 2 — Final Integration Plan (Post-Audit)

## Expanded Semi-Quantitative Candidate Library v5.7

**Date:** 2026-04-12
**Status:** Final — incorporates all audit findings from doc 23
**Prepared by:** Scientific Advisor + Architect + Dev-Orchestrator

---

## 1. Audit Finding Disposition

All 18 findings from the third-party audit (doc 23) addressed:

| ID | Severity | Finding | Disposition |
|---|---|---|---|
| **F1** | High | "Production-grade" claim overstated | **ACCEPT** — Rename to "expanded semi-quantitative candidate library" everywhere |
| **F2** | High | Spacer-as-multiplier is not mechanistic | **ACCEPT** — Keep Phase 1 multiplier, label as empirical accessibility correction |
| **F3** | High | Q reagent identity not explicit | **ACCEPT** — Add `reagent_identity="glycidyltrimethylammonium chloride"`, mark as functional approximation |
| **F4** | Med-High | CM chemistry ambiguous | **ACCEPT** — Use `reagent_identity="chloroacetic acid"`, document as CM-like weak cation exchanger |
| **F5** | High | NTA/IDA need metal-loading state | **ACCEPT** — Add `metal_ion`, `metal_loaded_fraction` fields to ReagentProfile for IMAC profiles |
| **F6** | High | Heparin should be macromolecule | **ACCEPT** — Set `is_macromolecule=True`, use ligand_accessible_area |
| **F7** | High | Streptavidin near-irreversible | **ACCEPT** — Cap effective stoichiometry at 2.5 (not 4), add `binding_model_hint="near_irreversible"` |
| **F8** | Med-High | Glutathione orientation | **ACCEPT** — Add `activity_retention=0.80` with uncertainty 0.15 (not 1.0) |
| **F9** | Med-High | HIC q_max not mappable | **ACCEPT** — Already `not_mapped`; confirm in code |
| **F10** | Medium | pH is validity gate not kinetic model | **ACCEPT** — Continue as-is, document limitation |
| **F11** | Positive | ACS accounting reliable | Acknowledged |
| **F12** | Positive | Solvers appropriate for Phase 1 | Acknowledged |
| **F13** | High | q_max area basis mismatch | **ACCEPT** — Add `ligand_density_area_basis` and `q_max_area_basis_note` to FMC |
| **F14** | High | charge_type mandatory | **ACCEPT** — Implement before adding Q/CM |
| **F15** | Med-High | New modes need binding_model_hint | **ACCEPT** — Add to FunctionalMediaContract |
| **F16** | Med-High | UI must expose confidence/hazards | **ACCEPT** — Display below reagent caption |
| **F17** | Medium | Spacer profiles must not be executable | **ACCEPT** — Filter by `reaction_type != "spacer"` |
| **F18** | Low-Med | Profile count inconsistent | **ACCEPT** — Canonical count: 14 existing + 8 coupling + 3 spacer = 25 total |

---

## 2. Revised Data Model Changes

### 2.1 ReagentProfile — 7 New Fields (expanded from 4 in doc 22)

```python
# ── Spacer arm (Phase 1 multiplier) ──
spacer_key: str = ""
spacer_length_angstrom: float = 0.0
spacer_activity_multiplier: float = 1.0

# ── Charge type for IEX (audit F14) ──
charge_type: str = ""                # "anion", "cation", ""

# ── IMAC metal state (audit F5) ──
metal_ion: str = ""                  # "Ni2+", "Co2+", "Cu2+", "Zn2+", ""
metal_loaded_fraction: float = 1.0   # [0,1] assumed metal loading; 1.0 = fully loaded

# ── Binding model hint for M3 (audit F15) ──
binding_model_hint: str = ""         # "charge_exchange", "metal_chelation",
                                     # "salt_promoted", "fc_affinity",
                                     # "gst_glutathione", "near_irreversible",
                                     # "mixed_mode", ""
```

### 2.2 FunctionalMediaContract — 3 New Fields (audit F13, F15)

```python
# ── Area basis (audit F13) ──
ligand_density_area_basis: str = ""  # "reagent_accessible", "ligand_accessible", "external"
q_max_area_basis_note: str = ""      # Human-readable note on q_max derivation

# ── Binding model hint (audit F15) ──
binding_model_hint: str = ""         # Passed from reagent profile for M3 routing
```

### 2.3 ACS / StepType Changes: NONE

No changes to ACSSiteType, ACSProfile, or ModificationStepType.

---

## 3. Revised Candidate Parameters

All candidates incorporate audit corrections:

### 3.1 Ligand Coupling (6 new profiles)

| Key | Reagent Identity (F3/F4) | Installed Ligand | charge_type (F14) | binding_model_hint (F15) | is_macro (F6) | activity_ret | confidence |
|---|---|---|---|---|---|---|---|
| `q_coupling` | Glycidyltrimethylammonium chloride | Q (quaternary ammonium) | **anion** | charge_exchange | false | 1.0 | semi_quantitative |
| `cm_coupling` | Chloroacetic acid + amino spacer | CM-like carboxymethyl | **cation** | charge_exchange | false | 1.0 | semi_quantitative |
| `nta_coupling` | Nitrilotriacetic acid | NTA chelator | — | metal_chelation | false | 1.0 | semi_quantitative |
| `butyl_coupling` | n-Butylamine | Butyl (HIC) | — | salt_promoted | false | 1.0 | semi_quantitative |
| `glutathione_coupling` | Glutathione (reduced, GSH) | Glutathione | — | gst_glutathione | false | **0.80** (F8) | semi_quantitative |
| `heparin_coupling` | Heparin sodium (porcine) | Heparin | — | mixed_mode | **true** (F6) | 1.0 | ranking_only |

**IMAC profiles (F5):** NTA and existing IDA get metal state fields:

| Key | metal_ion | metal_loaded_fraction | Notes |
|---|---|---|---|
| `nta_coupling` | Ni2+ | 1.0 | Assumed fully Ni-loaded; no leaching modeled |
| `ida_coupling` (existing) | Ni2+ | 1.0 | Backfill same assumption |

### 3.2 Protein Coupling (2 new profiles)

| Key | Reagent Identity | binding_model_hint (F15) | activity_ret | Effective Stoich (F7) | spacer | multiplier |
|---|---|---|---|---|---|---|
| `protein_ag_coupling` | Recombinant Protein A/G | fc_affinity | 0.55 | 2.0 | DADPA | 1.22 |
| `streptavidin_coupling` | Streptavidin (*S. avidinii*) | **near_irreversible** (F7) | 0.70 | **2.5** (F7, not 4) | DADPA | 1.22 |

**Streptavidin stoichiometry correction (F7):** Theoretical max = 4 biotin sites per tetramer. Practical effective stoichiometry = 2.5 due to steric occlusion of buried binding sites and orientation-dependent accessibility. Capped in FMC mapping.

### 3.3 Spacer Profiles (3 metadata-only profiles)

| Key | reaction_type | Notes |
|---|---|---|
| `dadpa_spacer` | **spacer** | NOT executable as step. Metadata only. |
| `aha_spacer` | **spacer** | NOT executable. AHA -COOH distal chemistry not modeled in Phase 1. |
| `dah_spacer` | **spacer** | NOT executable. Hydrophobicity/charge not modeled. |

### 3.4 Backfill Existing Profiles

| Existing Profile | New Fields Added |
|---|---|
| `deae_coupling` | `charge_type="anion"`, `binding_model_hint="charge_exchange"` |
| `sp_coupling` | `charge_type="cation"`, `binding_model_hint="charge_exchange"` |
| `ida_coupling` | `metal_ion="Ni2+"`, `metal_loaded_fraction=1.0`, `binding_model_hint="metal_chelation"` |
| `phenyl_coupling` | `binding_model_hint="salt_promoted"` |
| `protein_a_coupling` | `binding_model_hint="fc_affinity"` |
| `protein_g_coupling` | `binding_model_hint="fc_affinity"` |

---

## 4. Revised FunctionalMediaContract Logic

### 4.1 IEX Classification (F14 fix)

Replace string heuristic with:

```python
if rp.charge_type == "anion":
    ligand_type = "iex_anion"
elif rp.charge_type == "cation":
    ligand_type = "iex_cation"
```

### 4.2 q_max Mapping with Area Basis (F13 fix)

```python
# Determine area basis from reagent profile
if rp.is_macromolecule:
    area_per_particle = surface.ligand_accessible_area
    area_basis = "ligand_accessible"
else:
    area_per_particle = surface.reagent_accessible_area
    area_basis = "reagent_accessible"

# Specific surface area per bed volume
a_v = 6.0 * (1 - eps_bed) / d_particle  # external geometric
# Correction: use accessible area ratio if available
area_ratio = area_per_particle / (4 * math.pi * (d_particle/2)**2)  # accessible/external
a_v_accessible = a_v * area_ratio

q_max = functional_density * a_v_accessible * binding_stoich
```

### 4.3 Binding Stoichiometry Table (incorporating F7)

| ligand_type | binding_stoich | basis |
|---|---|---|
| iex_anion / iex_cation | 1.0 | Approximate; charge-dependent |
| imac | 1.0 | One His-tag per chelator |
| hic | N/A | `not_mapped` |
| affinity (Protein A/G) | 2.0 | Steric limit ~2 IgG per Protein A |
| biotin_affinity | **2.5** (not 4) | F7: steric occlusion reduces from theoretical 4 |
| gst_affinity | 1.0 | 1:1 GST:glutathione |
| heparin_affinity | 1.0 | Approximate; variable by target |

### 4.4 IMAC q_max Assumption Note (F5)

```python
if ligand_type == "imac":
    q_max_area_basis_note = (
        f"IMAC q_max assumes fully {rp.metal_ion}-loaded chelator "
        f"(metal_loaded_fraction={rp.metal_loaded_fraction:.0%}), "
        f"no metal leaching, no competing chelators. "
        f"Calibrate with batch uptake experiment."
    )
```

---

## 5. Revised UI Requirements (F16, F17)

### 5.1 Confidence Display (F16)

Below each reagent caption, add:

```python
st.caption(
    f"Confidence: {_profile.confidence_tier} | "
    f"Source: {_profile.calibration_source[:50]} | "
    f"Hazard: {_profile.hazard_class or 'low'}"
)
```

### 5.2 Spacer Dropdown Filtering (F17)

```python
# Filter spacer profiles from executable dropdowns
_reagent_options = {
    label: key for label, key in _all_options.items()
    if _REAGENT_PROFILES[key].reaction_type != "spacer"
}
```

Spacer selection appears as a separate optional field:

```python
if step_type in ("Ligand Coupling", "Protein Coupling"):
    _spacer_choice = st.selectbox(
        "Spacer Arm (optional)",
        ["None (Direct Coupling)", "DADPA (13 A)", "AHA (10 A)", "DAH (9 A)"],
        key=f"m2_spacer_{i}_{_sk}",
    )
```

---

## 6. Implementation Work Nodes

| # | Node | Tier | LOC | Depends On |
|---|---|---|---|---|
| WN-1 | Add 7 new ReagentProfile fields + backfill existing profiles | Opus | 60 | — |
| WN-2 | Add 6 ligand coupling profiles | Sonnet | 180 | WN-1 |
| WN-3 | Add 2 protein + 3 spacer profiles | Sonnet | 120 | WN-1 |
| WN-4 | Apply spacer multiplier + glutathione activity_ret in solvers | Sonnet | 15 | WN-1 |
| WN-5 | Extend FMC: charge_type routing, area basis, binding hints, IMAC notes | Opus | 80 | WN-1, WN-2 |
| WN-6 | UI: extend dropdowns, add spacer selectbox, display confidence/hazard | Sonnet | 60 | WN-2, WN-3 |
| WN-7 | Tests: profile completeness, charge_type, spacer filtering | Sonnet | 80 | WN-2, WN-3 |
| WN-8 | Tests: solver spacer multiplier, glutathione activity_ret cap | Sonnet | 40 | WN-4 |
| WN-9 | Tests: FMC q_max mapping, area basis, binding hints | Sonnet | 60 | WN-5 |
| WN-10 | Docstrings + rename "production-grade" → "semi-quantitative" | Haiku | 20 |All |
| **Total** | | | **~715** | |

### Parallelization

```
WN-1 → WN-2 + WN-3 (parallel) → WN-4 + WN-5 + WN-6 (parallel) → WN-7 + WN-8 + WN-9 (parallel) → WN-10
```

---

## 7. Revised Acceptance Criteria (from audit Section 10)

- [ ] 25 total profiles: 14 existing + 8 coupling + 3 spacer
- [ ] Every IEX profile has explicit `charge_type`; no string matching for anion/cation
- [ ] Spacer profiles filtered from executable dropdowns (`reaction_type != "spacer"`)
- [ ] `activity_retention * spacer_activity_multiplier` capped at 1.0
- [ ] IMAC profiles carry `metal_ion` and `metal_loaded_fraction`
- [ ] HIC q_max = `not_mapped`
- [ ] Streptavidin effective stoichiometry = 2.5 (not 4)
- [ ] Heparin has `is_macromolecule=True`
- [ ] Glutathione has `activity_retention=0.80`
- [ ] FMC includes `ligand_density_area_basis`, `q_max_area_basis_note`, `binding_model_hint`
- [ ] UI displays confidence_tier, calibration_source, hazard_class for every profile
- [ ] Existing 14 profiles pass regression (unchanged behavior)
- [ ] All existing 142 tests pass
- [ ] New tests cover: charge_type, spacer filtering, q_max area basis, IMAC assumptions, binding hints
- [ ] No document refers to "production-grade" — all say "semi-quantitative"

---

## 8. What is Deferred

### Phase 1.5 (near-term)
- `metal_charging` step type for IMAC
- `binding_model_hint` routing in M3 (M3 currently ignores it)
- `accessible_area_per_bed_volume` field on FMC
- Orientation/activity uncertainty metadata for proteins

### Phase 2 (future)
- SPACER_ARM as executable step type
- EDC/NHS chemistry for AHA carboxyl-terminated spacer
- pH-dependent rate scaling
- Salt-dependent HIC isotherms
- Reversible/irreversible affinity binding modes in M3

---

> **Disclaimer**: This is an expanded semi-quantitative candidate library. All kinetic parameters are order-of-magnitude estimates requiring experimental calibration. Spacer arm improvement factors are empirical corrections. q_max mappings are approximate and area-basis-dependent. Not suitable for quantitative wet-lab prediction without calibration against measured resin specifications.
