# M2 Candidates Integration Plan (v5.7)

## Consolidated Plan: Ligand, Protein, and Linker Arm Integration

**Date:** 2026-04-12
**Roles:** Scientific Advisor + Architect + Dev-Orchestrator
**Scope:** Expand Module 2 reagent library from 14 to 28 profiles; add spacer arm support

---

## PART 1 -- SCIENTIFIC ADVISOR

### 1.1 Candidate Triage: ADD NOW (v5.7) vs DEFER vs REJECT

#### ADD NOW -- Priority 1 (8 new reagent profiles)

| # | Candidate | Step Type | Key | Rationale |
|---|-----------|-----------|-----|-----------|
| 1 | **Q** (quaternary ammonium) | LIGAND_COUPLING | `q_coupling` | Completes IEX quartet: strong anion exchanger. Non-titratable charge at all pH. |
| 2 | **CM** (carboxymethyl) | LIGAND_COUPLING | `cm_coupling` | Completes IEX quartet: weak cation exchanger. |
| 3 | **NTA** (nitrilotriacetic acid) | LIGAND_COUPLING | `nta_coupling` | Industry-standard His-tag IMAC (Ni-NTA). Pairs with existing IDA. |
| 4 | **Butyl** (butylamine) | LIGAND_COUPLING | `butyl_coupling` | Lower-hydrophobicity HIC alternative to Phenyl. Covers practical HIC range. |
| 5 | **Glutathione** | LIGAND_COUPLING | `glutathione_coupling` | GST-tag affinity purification. Small molecule, full pore access. |
| 6 | **Heparin** | LIGAND_COUPLING | `heparin_coupling` | Dual-mode affinity + IEX. One of the most commercially important affinity resins. |
| 7 | **Protein A/G** (fusion) | PROTEIN_COUPLING | `protein_ag_coupling` | Broadest IgG subclass/species coverage. Fills gap between Protein A and G. |
| 8 | **Streptavidin** | PROTEIN_COUPLING | `streptavidin_coupling` | Biotin-tag capture. Kd ~10^-15 M. Widely used in bioprocessing. |

#### ADD NOW -- Priority 1 (3 spacer arm profiles)

| # | Spacer | Key | Length | Distal Group | Rationale |
|---|--------|-----|--------|--------------|-----------|
| 9 | **DADPA** | `dadpa_spacer` | 13 A (9 atoms) | -NH2 | EAH-Sepharose standard. Best for protein ligands. |
| 10 | **AHA** (6-aminohexanoic acid) | `aha_spacer` | 10 A (7 atoms) | -COOH | NHS-Sepharose HP standard. Enables EDC/NHS coupling path. |
| 11 | **DAH** (1,6-diaminohexane) | `dah_spacer` | 9 A (6 atoms) | -NH2 | AH-Sepharose standard. Simpler alternative to DADPA. |

#### DEFER to v5.8+ (Priority 2 -- adds coverage but not critical)

| Candidate | Reason for Deferral |
|-----------|-------------------|
| **Protein L** | Kappa light chain affinity -- niche (Fab/scFv). Add when users request. |
| **Concanavalin A** | Lectin affinity -- large tetramer (104 kDa). Needs long PEG spacer support first. |
| **Octyl** (octylamine) | Very hydrophobic HIC. Risk of irreversible binding. Butyl + Phenyl cover the range. |
| **WGA** (wheat germ agglutinin) | Complementary lectin. Low priority without Con A already in library. |
| **EDA** (ethylenediamine spacer) | Short spacer. Marginal benefit for small-molecule ligands that work fine direct. |
| **PEG-diamine Mn 600** | Long PEG spacer. Needed for Con A but that is itself deferred. |
| **BDGE** (bis-epoxide spacer) | Already inherent in ECH activation path. Redundant until two-stage activation modeled. |
| **TMAE** | Alternative Q route via DVS. Q already covers strong anion exchange via epoxide. |

#### REJECT (insufficient data or wrong chemistry)

| Candidate | Rejection Rationale |
|-----------|-------------------|
| **Cibacron Blue 3GA** | Requires direct triazine-OH coupling to unactivated agarose at 60C -- completely different activation chemistry from EPOXIDE/VS path. Would need a new activation step type AND new ACS site type (TRIAZINE). Not compatible with current M2 architecture. |
| **Reactive Red 120** | Same triazine-OH coupling problem as Cibacron Blue, plus niche application. |
| **Capto MMC / Capto Adhere** | Proprietary multi-step ligand synthesis. No public kinetic parameters. Cannot be calibrated. |
| **Anti-FLAG M2 antibody** | Very expensive, very low capacity. Not suitable for default simulation library. |
| **TED chelator** | Pentadentate -- very high specificity but very low capacity. Niche. |
| **Polymer brushes** (Dextran/PEI/Polyallylamine) | Fundamentally different "tentacle" functionalization. Requires grafting model, not stepwise coupling. Incompatible with current ACS state machine. |
| **SM(PEG)n crosslinkers** (P1-P4) | Heterobifunctional NHS-maleimide. Require amino-activated surface + thiol on protein. Not compatible with current epoxide activation path without intermediate spacer providing -NH2 terminus. Defer until amino-spacer path is validated. |
| **Lentil Lectin (LCA)** | Narrower specificity than Con A. If Con A is deferred, LCA has no case. |
| **Poly-glycine (Gly4)** | Peptide linker validated for rProtein A but requires custom synthesis and lacks commercial CAS/sourcing for general use. |

### 1.2 Spacer Arm Integration: Phase 1 Approach

#### Decision: Spacer-as-ReagentProfile-Attribute (NOT new step type)

**Phase 1 approach (v5.7):** Model spacer arms as **three new fields on ReagentProfile**:
- `spacer_key: str = ""` -- references a spacer profile (empty = direct coupling)
- `spacer_length_angstrom: float = 0.0` -- spacer length in angstroms
- `spacer_activity_multiplier: float = 1.0` -- multiplier on `activity_retention`

The spacer effect is applied as a **post-coupling multiplier** on `activity_retention` and `ligand_functional_sites`:

```
effective_activity = activity_retention * spacer_activity_multiplier
```

**Why NOT a new SPACER_ARM step type (Phase 2):**
1. **ACS constraint:** The ACS state model was just revised (v2 terminal states). Adding a new step type that creates intermediate ACS states (AMINE_DISTAL, CARBOXYL_DISTAL) would require new ACSSiteType enum values and new workflow routing in modification_steps.py. This violates the "no architectural changes to ACS" constraint.
2. **Complexity:** A SPACER_ARM step type requires a full new workflow solver (`_solve_spacer_step`), validation rules for spacer-coupling compatibility, and UI changes for two-step configuration.
3. **Diminishing returns:** For the 3 Priority 1 spacers, the improvement factor is 1.15-1.30x on activity_retention. A multiplier captures this to within the uncertainty band of the kinetic parameters themselves.

#### Scientific Validity of the Multiplier Simplification

The spacer-as-multiplier approach is **scientifically valid within the semi_quantitative confidence tier** because:

1. **The primary effect of a spacer is on activity retention**, not on coupling kinetics. The coupling reaction itself (epoxide ring-opening by amine) is chemically identical whether the amine is on the spacer terminus or on the protein directly.

2. **Published improvement factors** (1.05-1.40x) are empirical multipliers measured as ratio of active-ligand-with-spacer to active-ligand-without. This is exactly the abstraction we use.

3. **The uncertainty on activity_retention** (typically +/-0.15 for proteins) is comparable to or larger than the spacer improvement factor for short spacers. The multiplier is within the noise floor for small-molecule ligands (which show only 1.05-1.10x improvement).

4. **Limitation acknowledged:** This approach cannot model spacer-dependent changes in coupling kinetics (e.g., AHA-COOH terminus requiring EDC/NHS activation vs. DADPA-NH2 terminus coupling directly to epoxide). Phase 2 would address this with an explicit step.

#### Spacer-Ligand Compatibility Table (encoded in profiles)

| Spacer | Recommended For | Multiplier Range | Notes |
|--------|----------------|------------------|-------|
| None (direct) | IEX ligands (DEAE, Q, SP, CM), HIC (Phenyl, Butyl), IDA | 1.0 | Small molecules have full pore access |
| DADPA (13 A) | Protein A, G, A/G, L, Streptavidin | 1.15-1.30 | Standard for protein ligands |
| AHA (10 A) | Glutathione, Heparin | 1.10-1.20 | Provides -COOH for alternative coupling |
| DAH (9 A) | NTA, IDA (optional) | 1.05-1.10 | Simpler than DADPA, shorter chain |

### 1.3 FunctionalMediaContract: New Ligand Type Mappings

The current `ligand_type` field supports: `"iex_anion"`, `"iex_cation"`, `"affinity"`, `"imac"`, `"hic"`, `"none"`.

New candidates require these additional `ligand_type` values:

| New ligand_type | Triggered By | q_max Mapping |
|----------------|-------------|---------------|
| `"gst_affinity"` | Glutathione (`functional_mode="gst_affinity"`) | q_max = functional_density * a_v * 1.0 (1:1 GST:glutathione) |
| `"biotin_affinity"` | Streptavidin (`functional_mode="biotin_affinity"`) | q_max = functional_density * a_v * 4.0 (4 biotin sites per streptavidin tetramer) |
| `"heparin_affinity"` | Heparin (`functional_mode="heparin_affinity"`) | q_max = functional_density * a_v * 1.0 (approximate; variable stoichiometry) |

The existing `"affinity"` type with binding_stoich=2.0 for Protein A/G is retained. Protein A/G fusion uses the same mapping.

---

## PART 2 -- ARCHITECT

### 2.1 Data Model Changes

#### 2.1.1 ReagentProfile: New Fields (3 fields)

Add to `reagent_profiles.py` `ReagentProfile` dataclass:

```python
# -- Spacer arm support (Phase 1 multiplier model) --
spacer_key: str = ""                    # Reference to spacer profile key (empty = direct)
spacer_length_angstrom: float = 0.0     # Spacer length [angstrom]
spacer_activity_multiplier: float = 1.0 # Multiplier on activity_retention [>=1.0]
```

#### 2.1.2 ReagentProfile: New functional_mode Values

Current values: `"crosslinker"`, `"activator"`, `"iex_ligand"`, `"affinity_ligand"`, `"hic_ligand"`, `"imac_chelator"`, `"quencher"`

Add:
- `"gst_affinity"` -- for Glutathione
- `"biotin_affinity"` -- for Streptavidin
- `"heparin_affinity"` -- for Heparin
- `"spacer"` -- for spacer arm profiles (informational only; spacers are not independent steps in Phase 1)

#### 2.1.3 ACS Changes: NONE

No changes to `ACSSiteType` enum, `ACSProfile` dataclass, or conservation invariants.

#### 2.1.4 ModificationStepType: NONE

No new step types. The 5 existing types are sufficient for Phase 1.

#### 2.1.5 Spacer Profile Data Structure

Spacer arm profiles are stored in the same `REAGENT_PROFILES` dict but are NOT used as independent modification steps. They serve as lookup records for the `spacer_key` field on coupling profiles.

```python
# Spacer profiles -- metadata only, not executed as steps
"dadpa_spacer": ReagentProfile(
    name="DADPA (diaminodipropylamine spacer)",
    cas="56-18-8",
    reaction_type="spacer",       # new reaction_type value
    functional_mode="spacer",
    target_acs=ACSSiteType.EPOXIDE,
    product_acs=None,
    k_forward=0.0,  # not used in Phase 1
    E_a=0.0,
    stoichiometry=1.0,
    spacer_length_angstrom=13.0,
    spacer_activity_multiplier=1.22,  # midpoint of 1.15-1.30 range
    ...
)
```

### 2.2 Modification to Protein Coupling Solver

In `_solve_protein_coupling_step()` (modification_steps.py), apply spacer multiplier:

```python
# After computing activity_ret from reagent_profile:
activity_ret = getattr(reagent_profile, 'activity_retention', 1.0)
spacer_mult = getattr(reagent_profile, 'spacer_activity_multiplier', 1.0)
effective_activity = min(activity_ret * spacer_mult, 1.0)  # cap at 1.0
```

Same change in `_solve_ligand_coupling_step()` for completeness (though multiplier is 1.0 for small molecules by default).

### 2.3 FunctionalMediaContract Extensions

In `orchestrator.py`, `build_functional_media_contract()`:

1. Extend `_mode_map` to handle new `functional_mode` values:

```python
_mode_map = {
    "iex_ligand": ...,           # existing
    "affinity_ligand": "affinity",  # existing
    "imac_chelator": "imac",     # existing
    "hic_ligand": "hic",         # existing
    "gst_affinity": "gst_affinity",       # NEW
    "biotin_affinity": "biotin_affinity", # NEW
    "heparin_affinity": "heparin_affinity", # NEW
}
```

2. Add q_max mapping branches for new ligand types:

```python
elif ligand_type == "gst_affinity" and functional_density > 0:
    binding_stoich = 1.0  # 1 GST per glutathione
    q_max_est = functional_density * a_v * binding_stoich
    ...
elif ligand_type == "biotin_affinity" and functional_density > 0:
    binding_stoich = 4.0  # 4 biotin per streptavidin tetramer
    q_max_est = functional_density * a_v * binding_stoich
    ...
elif ligand_type == "heparin_affinity" and functional_density > 0:
    binding_stoich = 1.0  # approximate
    q_max_est = functional_density * a_v * binding_stoich
    ...
```

### 2.4 UI Changes (app.py)

Extend the reagent dropdown `_reagent_options` dictionaries:

**Ligand Coupling dropdown** (line ~865):
```python
_reagent_options = {
    "DEAE (Weak Anion Exchange)": "deae_coupling",
    "Q (Strong Anion Exchange)": "q_coupling",           # NEW
    "CM (Weak Cation Exchange)": "cm_coupling",           # NEW
    "Sulfopropyl (Strong Cation)": "sp_coupling",
    "IDA (IMAC Chelator)": "ida_coupling",
    "NTA (IMAC Chelator, His-tag)": "nta_coupling",       # NEW
    "Phenyl (HIC)": "phenyl_coupling",
    "Butyl (HIC, Mild)": "butyl_coupling",                 # NEW
    "Glutathione (GST-tag Affinity)": "glutathione_coupling", # NEW
    "Heparin (Affinity + IEX)": "heparin_coupling",       # NEW
}
```

**Protein Coupling dropdown** (line ~872):
```python
_reagent_options = {
    "Protein A (IgG Affinity)": "protein_a_coupling",
    "Protein G (IgG Broad Subclass)": "protein_g_coupling",
    "Protein A/G Fusion (Broadest IgG)": "protein_ag_coupling",  # NEW
    "Streptavidin (Biotin-tag)": "streptavidin_coupling",        # NEW
}
```

**Spacer selection (NEW UI element):**
Add an optional spacer selectbox that appears when a Ligand Coupling or Protein Coupling step is selected:

```python
_spacer_options = {
    "None (Direct Coupling)": "",
    "DADPA (13 A, EAH-standard)": "dadpa_spacer",
    "AHA (10 A, NHS-standard)": "aha_spacer",
    "DAH (9 A, AH-standard)": "dah_spacer",
}
```

The selected spacer key overrides the coupling profile's `spacer_activity_multiplier` at runtime.

### 2.5 File-Level Change List

| # | File | Changes | Est. LOC |
|---|------|---------|----------|
| 1 | `reagent_profiles.py` | Add 3 dataclass fields, 8 coupling profiles, 3 spacer profiles | +280 |
| 2 | `modification_steps.py` | Apply spacer_activity_multiplier in ligand + protein solvers (2 functions, ~5 lines each) | +12 |
| 3 | `orchestrator.py` | Extend _mode_map (3 entries), add 3 q_max branches, extend iex_anion/cation detection | +45 |
| 4 | `visualization/app.py` | Extend 2 dropdown dicts, add spacer selectbox, wire spacer_key to step construction | +40 |
| 5 | `tests/test_reagent_profiles.py` | Validate all 25 profiles have required metadata fields, test spacer multiplier bounds | +60 |
| 6 | `tests/test_modification_steps.py` | Test spacer multiplier applied in protein/ligand coupling, test new profiles execute | +80 |
| 7 | `tests/test_orchestrator.py` | Test new ligand_type mappings in FunctionalMediaContract, test q_max for new types | +50 |
| **Total** | | | **~567 LOC** |

---

## PART 3 -- DEV-ORCHESTRATOR

### 3.1 Work Node Breakdown

```
WN-1  [Opus]    Spacer field addition to ReagentProfile dataclass
WN-2  [Sonnet]  8 new coupling ReagentProfile definitions (Q, CM, NTA, Butyl, Glutathione, Heparin, Protein A/G, Streptavidin)
WN-3  [Sonnet]  3 spacer arm ReagentProfile definitions (DADPA, AHA, DAH)
WN-4  [Sonnet]  Apply spacer_activity_multiplier in _solve_ligand_coupling_step + _solve_protein_coupling_step
WN-5  [Sonnet]  Extend FunctionalMediaContract: _mode_map + q_max branches
WN-6  [Sonnet]  UI: extend dropdowns + add spacer selectbox in app.py
WN-7  [Sonnet]  Test: profile validation (all 25 profiles, metadata completeness)
WN-8  [Sonnet]  Test: coupling solver with spacer multiplier
WN-9  [Sonnet]  Test: FunctionalMediaContract new ligand type mappings
WN-10 [Haiku]   Docstring and comment updates across all changed files
```

### 3.2 Dependency Graph

```
WN-1 ─────┬──> WN-2 ──┬──> WN-4 ──> WN-8
           │           │
           ├──> WN-3 ──┘
           │
           └──────────────> WN-5 ──> WN-9
                            
WN-2 ──> WN-6 ──> WN-7
WN-3 ──┘

WN-10 depends on ALL of WN-1..WN-9
```

### 3.3 Parallelization Plan

| Phase | Parallel Nodes | Blocking Dependency |
|-------|---------------|-------------------|
| **Phase A** | WN-1 (Opus) | None -- must complete first (dataclass changes) |
| **Phase B** | WN-2, WN-3 (Sonnet x2) | WN-1 complete |
| **Phase C** | WN-4, WN-5, WN-6 (Sonnet x3) | WN-2 + WN-3 complete |
| **Phase D** | WN-7, WN-8, WN-9 (Sonnet x3) | WN-4 + WN-5 + WN-6 complete |
| **Phase E** | WN-10 (Haiku) | All prior complete |

**Critical path:** WN-1 -> WN-2 -> WN-4 -> WN-8 (4 serial phases)

### 3.4 Model Tier Assignments

| Tier | Nodes | Rationale |
|------|-------|-----------|
| **Opus** | WN-1 | Dataclass schema change requires understanding conservation invariants and downstream impact |
| **Sonnet** | WN-2 through WN-9 | Profile definitions, solver modifications, UI wiring, and test writing are well-scoped implementation tasks |
| **Haiku** | WN-10 | Docstring/comment updates are mechanical text changes |

### 3.5 Effort Estimate

| Phase | Nodes | Est. Time (per node) | Calendar Time |
|-------|-------|---------------------|---------------|
| A | 1 | 15 min | 15 min |
| B | 2 (parallel) | 20 min each | 20 min |
| C | 3 (parallel) | 15 min each | 15 min |
| D | 3 (parallel) | 20 min each | 20 min |
| E | 1 | 10 min | 10 min |
| **Total** | 10 | ~170 min aggregate | **~80 min wall-clock** |

### 3.6 Risk Assessment

| Risk | Severity | Likelihood | Mitigation |
|------|----------|------------|------------|
| **Spacer multiplier > 1.0 pushes activity_retention above 1.0** | Medium | Medium | Cap at `min(effective_activity, 1.0)` in solver. Add assertion in test. |
| **New functional_mode values break existing _mode_map logic** | Medium | Low | The `_mode_map.get(fm, "none")` fallback handles unknown modes gracefully. Add explicit entries. |
| **IEX anion/cation detection heuristic fails for Q/CM** | Medium | Medium | Current heuristic checks for "anion" or "deae" in installed_ligand string. Must add "q" to anion list and "cm" to cation list, or switch to a ligand-to-charge-type lookup table. **Recommend: replace string heuristic with explicit charge_type field on ReagentProfile.** |
| **Spacer profiles appear in step-type dropdowns** | Low | Medium | Filter spacer profiles out of step-type reagent options by checking `reaction_type != "spacer"`. |
| **New profiles missing required metadata fields** | Low | Low | WN-7 test validates all profiles have non-empty confidence_tier, calibration_source, and hazard_class. |
| **FunctionalMediaContract grows too many ligand_type values** | Low | Low | 3 new values is manageable. Consider enum if it grows beyond 10. |
| **ACS conservation violation from spacer multiplier** | None | None | Spacer multiplier only affects `ligand_functional_sites` (which is already a subset of `ligand_coupled_sites`). No conservation fields are changed. |

### 3.7 Validation Checklist (Definition of Done)

- [ ] All 25 reagent profiles (14 existing + 8 coupling + 3 spacer) present in `REAGENT_PROFILES`
- [ ] Every profile has non-empty `confidence_tier`, `calibration_source`, `hazard_class`
- [ ] Spacer multiplier capped at 1.0 for `effective_activity` in both coupling solvers
- [ ] Spacer profiles excluded from step-type reagent dropdowns in UI
- [ ] FunctionalMediaContract correctly maps all new `functional_mode` values to `ligand_type`
- [ ] IEX anion/cation detection works for Q (anion) and CM (cation) without string heuristic failure
- [ ] Existing 14 profiles unchanged (regression test)
- [ ] All existing tests pass (no regressions)
- [ ] New tests cover: profile completeness, spacer multiplier arithmetic, q_max mapping for 3 new affinity types

---

## APPENDIX A -- New Profile Parameters (Quick Reference)

### A.1 Ligand Coupling Profiles

| Key | k_forward | E_a | k_hydrol | pH_opt | T_default | t_default | MW | Spacer |
|-----|-----------|-----|----------|--------|-----------|-----------|-----|--------|
| `q_coupling` | 6e-5 | 50,000 | 1e-5 | 10.5 | 298.15 | 14,400 | 104 | None |
| `cm_coupling` | 3e-5 | 45,000 | 1e-5 | 10.5 | 298.15 | 21,600 | 94 | None |
| `nta_coupling` | 2e-5 | 45,000 | 1e-5 | 10.5 | 298.15 | 21,600 | 191 | None |
| `butyl_coupling` | 4e-5 | 48,000 | 5e-6 | 9.5 | 298.15 | 14,400 | 73 | None |
| `glutathione_coupling` | 3e-5 | 45,000 | 1e-5 | 9.5 | 298.15 | 14,400 | 307 | AHA |
| `heparin_coupling` | 1e-5 | 40,000 | 5e-6 | 9.5 | 298.15 | 28,800 | 14,000 | DADPA |

### A.2 Protein Coupling Profiles

| Key | k_forward | E_a | activity_ret | max_density | MW | Spacer | Multiplier |
|-----|-----------|-----|-------------|-------------|-----|--------|------------|
| `protein_ag_coupling` | 5e-7 | 25,000 | 0.55 | 2e-8 | 51,000 | DADPA | 1.22 |
| `streptavidin_coupling` | 4e-7 | 25,000 | 0.70 | 3e-8 | 53,000 | DADPA | 1.22 |

### A.3 Spacer Profiles

| Key | CAS | Length (A) | Distal Group | Multiplier | Hazard |
|-----|-----|-----------|-------------|------------|--------|
| `dadpa_spacer` | 56-18-8 | 13 | -NH2 | 1.22 | irritant |
| `aha_spacer` | 60-32-2 | 10 | -COOH | 1.15 | low_hazard |
| `dah_spacer` | 124-09-4 | 9 | -NH2 | 1.08 | irritant |

---

## APPENDIX B -- Architectural Decision Record

### ADR-1: Spacer as multiplier, not step type

**Status:** Accepted for v5.7 (Phase 1)
**Context:** Spacer arms scientifically operate as a separate coupling step between activation and ligand coupling. However, the ACS state model was just revised to v2 terminal states and adding a new step type would require new ACSSiteType values and solver routing.
**Decision:** Model spacer as a multiplier on activity_retention stored as a ReagentProfile attribute.
**Consequences:** Cannot model spacer-dependent coupling kinetics (e.g., COOH-terminus requiring EDC/NHS). Acceptable at semi_quantitative tier. Phase 2 can add SPACER_ARM step type if users need explicit spacer workflow simulation.

### ADR-2: Separate spacer profiles vs inline fields

**Status:** Accepted
**Context:** Spacer parameters could be stored directly on coupling profiles (inline) or as separate profiles referenced by key.
**Decision:** Both. Spacer profiles stored in REAGENT_PROFILES as metadata records. Coupling profiles reference them via `spacer_key` and carry the effective `spacer_activity_multiplier`.
**Consequences:** Spacer profiles must be filtered out of step-type dropdowns. The dual storage is slightly redundant but enables future Phase 2 where spacers become independent steps.

### ADR-3: Replace IEX charge heuristic with explicit field

**Status:** Proposed (implement in WN-5)
**Context:** Current `_mode_map` uses string matching on `installed_ligand` to distinguish anion vs cation IEX. Adding Q and CM makes this fragile.
**Decision:** Add `charge_type: str = ""` field to ReagentProfile with values `"anion"`, `"cation"`, `""`. Use this in `_mode_map` instead of string heuristic.
**Consequences:** One additional field on ReagentProfile. Must backfill on existing DEAE (anion) and SP (cation) profiles.

---

> **Disclaimer**: This integration plan is produced by an AI assistant acting in Scientific Advisor, Architect, and Dev-Orchestrator roles. All kinetic parameters are order-of-magnitude estimates requiring experimental calibration. Effort estimates assume familiarity with the codebase.
