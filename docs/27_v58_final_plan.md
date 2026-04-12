# EmulSim Module 2 v5.8 — Final Plan (Post-Audit)

## Deferred Candidates + SM(PEG)n Heterobifunctional Crosslinker Path

**Date:** 2026-04-12
**Status:** Final — incorporates all 15 audit findings from doc 26
**Baseline:** v5.7 (25 profiles, 5 step types, 152 tests passing)

---

## 1. Audit Finding Disposition

| ID | Severity | Finding | Disposition |
|---|---|---|---|
| **F1** | Positive | SPACER_ARM correctly identified as new step type | Acknowledged — proceed |
| **F2** | High | Intermediate ACS profiles can contaminate FMC ligand density | **ACCEPT** — Add `profile_role` field; FMC scans only `profile_role="final_ligand"` |
| **F3** | High | Maleimide hydrolysis not modeled after immobilization | **ACCEPT** — Add `maleimide_decay_rate` on MALEIMIDE ACSProfile; apply first-order decay before protein coupling |
| **F4** | High | Protein-Cys coupling needs distinct profiles | **ACCEPT** — Add `protein_a_cys_coupling` etc. with `chemistry_class="maleimide_thiol"` |
| **F5** | High | Buffer compatibility missing | **ACCEPT** — Add `buffer_incompatibilities` field; UI warnings for Tris/DTT |
| **F6** | High | TMAE chemistry under-specified | **ACCEPT** — Defer TMAE to v5.9 until reactive precursor clarified |
| **F7** | High | Con A / WGA need lectin-specific assumptions | **ACCEPT** — Add `binding_model_hint="lectin_mannose"` / `"lectin_glcnac"`, cofactor metadata |
| **F8** | Med-High | Octyl HIC threshold needs unit conversion | **ACCEPT** — Implement density conversion with area basis |
| **F9** | Med-High | BDGE as activation loses spacer info | **ACCEPT** — Carry `spacer_length_angstrom=18` on the product EPOXIDE profile via modification result |
| **F10** | Med-High | Diamines can bridge/crosslink | **ACCEPT** — Add `distal_group_yield` field (default 0.7 for diamines) |
| **F11** | Med-High | Rule 5 too strict for accumulation | **ACCEPT** — Allow accumulation on same product if same chemistry; block only on conflicting chemistry |
| **F12** | Medium | Profile count inconsistent | **ACCEPT** — Canonical count: 25 + 6 Tier A + 3 Cys-protein + 1 EDA-arm + 3 diamine-arm + 4 SM(PEG)n = 42 |
| **F13** | Medium | "Zero risk" is not zero risk | **ACCEPT** — Rename to "low architecture risk"; add tests |
| **F14** | Medium | UI support underestimated | **ACCEPT** — Full UI spec for SPACER_ARM step type |
| **F15** | High | M3 isotherm coverage insufficient | **ACCEPT** — Add `m3_support_level` field on every profile |

---

## 2. Revised Scope: What Ships in v5.8

### 2.1 Phase 1 (v5.8.0-rc1): Low-Architecture-Risk Profiles

6 new profiles (TMAE deferred per F6):

| # | Candidate | Type | binding_model_hint | m3_support_level |
|---|---|---|---|---|
| 1 | Protein L | PROTEIN_COUPLING | kappa_light_chain_affinity | requires_user_calibration |
| 2 | Concanavalin A | PROTEIN_COUPLING | lectin_mannose_glucose | requires_user_calibration |
| 3 | Octyl | LIGAND_COUPLING | salt_promoted | not_mapped |
| 4 | WGA | PROTEIN_COUPLING | lectin_glcnac_sialic | requires_user_calibration |
| 5 | PEG-diamine Mn 600 | spacer (metadata) | — | — |
| 6 | BDGE | ACTIVATION | — | — |

**Excluded from Phase 1:** TMAE (F6: chemistry under-specified), EDA (needs SPACER_ARM)

### 2.2 Phase 2 (v5.8.0-rc2): SPACER_ARM Architecture + SM(PEG)n

New architecture:
- 2 new ACSSiteType: `AMINE_DISTAL`, `MALEIMIDE`
- 1 new ModificationStepType: `SPACER_ARM`
- 1 new solver: `_solve_spacer_arm_step()`

New profiles (11):
- 1 EDA spacer-arm
- 3 diamine spacer-arm variants (DADPA, DAH, PEG-diamine as executable)
- 4 SM(PEG)n heterobifunctional (SM(PEG)2, 4, 12, 24)
- 3 Protein-Cys coupling (Protein A-Cys, Protein G-Cys, generic Cys-protein)

### 2.3 Canonical Profile Count

| Category | v5.7 | +Phase 1 | +Phase 2 | v5.8 Total |
|---|---|---|---|---|
| Crosslinkers | 2 | 0 | 0 | 2 |
| Activators | 2 | +1 (BDGE) | 0 | 3 |
| Ligand coupling | 10 | +1 (Octyl) | 0 | 11 |
| Protein coupling | 4 | +3 (Protein L, Con A, WGA) | +3 (Cys variants) | 10 |
| Quenchers | 4 | 0 | 0 | 4 |
| Spacer metadata | 3 | +1 (PEG-diamine) | 0 | 4 |
| Spacer-arm (executable) | 0 | 0 | +4 (EDA, DADPA-arm, DAH-arm, PEG-arm) | 4 |
| SM(PEG)n | 0 | 0 | +4 | 4 |
| **Total** | **25** | **+6 = 31** | **+11 = 42** | **42** |

---

## 3. Revised Data Model (Incorporating Audit Corrections)

### 3.1 New ReagentProfile Fields (7 new, on top of v5.7's existing fields)

```python
# ── Profile role (audit F2) ──
profile_role: str = "final_ligand"      # "native", "activated", "spacer_intermediate",
                                         # "heterobifunctional_intermediate", "final_ligand",
                                         # "spacer_metadata", "quencher"

# ── M3 support level (audit F15) ──
m3_support_level: str = "mapped_estimated"  # "mapped_quantitative", "mapped_estimated",
                                             # "not_mapped", "requires_user_calibration"

# ── Spacer-arm distal yield (audit F10) ──
distal_group_yield: float = 1.0         # Fraction of consumed sites producing distal group [0,1]
                                         # <1.0 for diamines (bridging side reaction)

# ── Maleimide decay (audit F3) ──
maleimide_decay_rate: float = 0.0       # [1/s] first-order hydrolysis of immobilized maleimide

# ── Buffer incompatibilities (audit F5) ──
buffer_incompatibilities: str = ""      # Comma-separated: "Tris,glycine,DTT,beta-ME"

# ── Thiol-specific coupling (audit F4) ──
requires_reduced_thiol: bool = False    # True for maleimide-thiol coupling
thiol_accessibility_fraction: float = 1.0  # Fraction of protein Cys accessible [0,1]
```

### 3.2 New ACSSiteType Values

```python
AMINE_DISTAL = "amine_distal"       # -NH2 at spacer terminus
MALEIMIDE = "maleimide"             # Maleimide at crosslinker terminus
```

### 3.3 ACSProfile: profile_role Field (audit F2)

Each ACSProfile carries `profile_role` to distinguish intermediates from final ligands:

| profile_role | Created By | FMC Uses? |
|---|---|---|
| `"native"` | M1 initialization (AMINE_PRIMARY, HYDROXYL) | No |
| `"activated"` | ACTIVATION (EPOXIDE, VINYL_SULFONE) | No |
| `"spacer_intermediate"` | SPACER_ARM (AMINE_DISTAL) | No |
| `"heterobifunctional_intermediate"` | SPACER_ARM (MALEIMIDE) | No |
| `"final_ligand"` | LIGAND_COUPLING, PROTEIN_COUPLING | **Yes** |

**FMC contract logic (F2 fix):** `build_functional_media_contract()` scans only profiles where the last coupling step type was LIGAND_COUPLING or PROTEIN_COUPLING. Intermediate SPACER_ARM products are excluded from functional density calculation.

### 3.4 Maleimide Decay (audit F3)

Before `_solve_protein_coupling_step()` executes on a MALEIMIDE target, apply first-order decay:

```python
if target_profile.site_type == ACSSiteType.MALEIMIDE:
    k_decay = getattr(reagent_profile, 'maleimide_decay_rate', 0.0)
    if k_decay > 0 and step.time > 0:
        # Time between maleimide creation and protein coupling
        # ASSUMPTION: delay_time = step.time (user-specified reaction time)
        fraction_remaining = math.exp(-k_decay * step.time)
        sites_decayed = target_profile.remaining_activated * (1 - fraction_remaining)
        target_profile.hydrolyzed_sites += sites_decayed
```

Typical k_decay for maleimide at pH 7.0, 25C: ~1e-5 /s (half-life ~19 hours).
At pH 7.5: ~5e-5 /s (half-life ~4 hours). This is significant for overnight protocols.

### 3.5 Distal Group Yield (audit F10)

In `_solve_spacer_arm_step()`:

```python
distal_yield = getattr(reagent_profile, 'distal_group_yield', 1.0)
sites_created = sites_consumed * distal_yield
# Remainder goes to bridging/crosslinking (no useful product)
sites_bridged = sites_consumed * (1 - distal_yield)
target_profile.crosslinked_sites += sites_bridged  # bridging is a terminal state
```

Default yields:
- EDA: 0.60 (high bridging risk, short chain)
- DAH: 0.70 (moderate bridging)
- DADPA: 0.80 (internal amine reduces bridging)
- PEG-diamine: 0.90 (PEG flexibility favors monoattachment)

---

## 4. Revised Work Nodes

| # | Node | Tier | LOC | Dependencies |
|---|---|---|---|---|
| **WN-0** | 6 Tier A/B profiles + tests | Sonnet | 250 | None |
| **WN-1** | 7 new ReagentProfile fields + ACSSiteType extensions | Opus | 50 | None |
| **WN-2** | SPACER_ARM step type + solver with distal_yield + profile_role | Opus | 150 | WN-1 |
| **WN-3** | Maleimide decay in protein coupling solver | Opus | 40 | WN-2 |
| **WN-4** | Workflow validation Rules 5-6 + intermediate-aware ordering | Opus | 60 | WN-2 |
| **WN-5** | 4 spacer-arm + 4 SM(PEG)n + 3 Cys-protein profiles | Sonnet | 300 | WN-1 |
| **WN-6** | FMC: profile_role filtering + m3_support_level | Opus | 60 | WN-2, WN-5 |
| **WN-7** | UI: SPACER_ARM step, intermediate display, buffer warnings | Sonnet | 80 | WN-5 |
| **WN-8** | Integration tests: full SM(PEG)n path + conservation | Opus | 150 | All above |
| **WN-9** | Docs: canonical profile table, confidence labels | Haiku | 20 | WN-8 |
| **Total** | | | **~1,160** | |

### Dependency Graph

```
WN-0 ──────────────────────────────────> (independent, ship as v5.8.0-rc1)

WN-1 ──┬──> WN-2 ──┬──> WN-3 ──> WN-8 ──> WN-9
        │           ├──> WN-4 ──┘
        │           └──> WN-6 ──┘
        └──> WN-5 ──────────────┘
                    └──> WN-7 ──┘
```

### Phased Release

| Release | Nodes | Profiles | Tests |
|---|---|---|---|
| v5.8.0-rc1 | WN-0 | 25 → 31 | ~165 |
| v5.8.0-rc2 | WN-1 through WN-8 | 31 → 42 | ~200 |
| v5.8.0 | WN-9 | 42 (final) | ~200 |

---

## 5. Acceptance Criteria (Revised per Audit Section 8)

### Phase 1 (v5.8.0-rc1)

- [ ] 31 total profiles present; all instantiate without error
- [ ] Con A has `binding_model_hint="lectin_mannose_glucose"` and cofactor metadata
- [ ] WGA has `binding_model_hint="lectin_glcnac_sialic"`
- [ ] Octyl has `m3_support_level="not_mapped"`
- [ ] Protein L has `m3_support_level="requires_user_calibration"`
- [ ] BDGE activation carries `spacer_length_angstrom=18` downstream
- [ ] PEG-diamine is spacer metadata only (not executable)
- [ ] All 152 existing tests pass (regression)
- [ ] New profile tests pass

### Phase 2 (v5.8.0-rc2)

- [ ] SPACER_ARM consumes EPOXIDE and creates AMINE_DISTAL with `distal_group_yield`
- [ ] SPACER_ARM consumes AMINE_DISTAL and creates MALEIMIDE via SM(PEG)n
- [ ] Maleimide decay applies before protein-Cys coupling
- [ ] Protein-Cys profiles have `chemistry_class="maleimide_thiol"` and `requires_reduced_thiol=True`
- [ ] FMC ignores intermediate profiles (`profile_role != "final_ligand"`)
- [ ] Buffer incompatibility warnings displayed in UI for NHS and maleimide steps
- [ ] Rule 5: SPACER_ARM product_acs blocks conflicting duplicate creation
- [ ] Rule 6: SM(PEG)n requires prior AMINE_DISTAL in workflow
- [ ] End-to-end test: ECH → DADPA(SPACER_ARM) → SM(PEG)4(SPACER_ARM) → Protein-Cys → Quench
- [ ] Conservation holds across all intermediate and final ACS profiles
- [ ] All existing v5.7 workflows unchanged (backward compatibility)
- [ ] 42 profiles, 6 step types, ~200 tests passing

---

## 6. What is Deferred to v5.9+

| Item | Reason |
|---|---|
| TMAE | Reactive precursor chemistry under-specified (audit F6) |
| EDC/NHS activation for AHA -COOH terminus | Requires new activation chemistry, not in scope |
| Salt-dependent HIC isotherms | M3 extension, not M2 |
| Lectin-specific M3 elution models | M3 extension |
| pH-dependent rate scaling | Global M2 enhancement |
| Metal charging step for IMAC | Separate step type |
| Protein disulfide reduction step | Pre-treatment chemistry |

---

> **Disclaimer**: This is an expanded semi-quantitative simulation library. All kinetic parameters are order-of-magnitude estimates. SM(PEG)n coupling outcomes depend on buffer, pH, timing, and protein-specific thiol accessibility not fully captured in Phase 1. Calibrate against measured resin specifications before quantitative use.
