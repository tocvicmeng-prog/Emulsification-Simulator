# EmulSim Module 2 v5.8 -- Deferred Candidates + SM(PEG)n Crosslinker Plan

**Date:** 2026-04-12
**Roles:** Scientific Advisor + Architect + Dev-Orchestrator
**Scope:** 8 deferred Priority 2 reagent profiles + SM(PEG)n heterobifunctional crosslinker path
**Baseline:** v5.7 (25 profiles, 5 step types, Phase 1 spacer-as-multiplier)

---

## PART 1 -- SCIENTIFIC ADVISOR

### 1.1 Deferred Candidate Scientific Parameters

#### 1.1.1 Protein L (kappa light chain affinity, 36 kDa)

- **Biology:** Protein L (Peptostreptococcus magnus) binds kappa light chains of immunoglobulins without interfering with antigen-binding site. Particularly useful for Fab, scFv, and single-domain antibodies that lack Fc regions.
- **MW:** 36 kDa monomer; r_h ~ 2.3 nm
- **Coupling chemistry:** Identical to Protein A/G -- epoxide-amine at 4C, 16h
- **k_forward:** 5e-7 m^3/(mol*s) (same order as Protein A -- diffusion-limited macromolecule)
- **E_a:** 25,000 J/mol
- **activity_retention:** 0.55 +/- 0.15 (orientation-sensitive; random coupling penalty)
- **max_surface_density:** 3e-8 mol/m^2 (smaller than Protein A, denser packing possible)
- **Spacer benefit:** DADPA (13A) recommended; multiplier 1.20
- **References:** Bjorck (1988) J. Immunol. 140:1194; Akerstrom & Bjorck (1989)

#### 1.1.2 Concanavalin A (lectin affinity, 104 kDa tetramer)

- **Biology:** Jack bean lectin, binds alpha-D-mannose and alpha-D-glucose residues on glycoproteins. Tetramer at pH > 7 (4 x 26 kDa), dimer below pH 6.
- **MW:** 104 kDa (tetramer); r_h ~ 4.0 nm
- **Coupling chemistry:** Epoxide-amine at 4C, requires long spacer (>25A) due to binding site geometry
- **k_forward:** 2e-7 m^3/(mol*s) (slower than Protein A due to larger size)
- **E_a:** 25,000 J/mol
- **activity_retention:** 0.40 +/- 0.20 (highly orientation-sensitive tetramer; some subunits occluded)
- **max_surface_density:** 1e-8 mol/m^2 (large footprint)
- **Spacer requirement:** PEG-diamine Mn 600 (35A) strongly recommended. With DADPA (13A) only, activity_retention drops to ~0.25. The spacer must present the Con A binding sites away from the matrix surface to accommodate the 104 kDa tetramer.
- **Optimal spacer length for Con A:** 25-40A. Below 15A, steric clash with the matrix reduces functional binding sites by ~50%. Published protocols (Pharmacia/Cytiva Con A Sepharose 4B) use 6-atom or longer spacer arms.
- **References:** Goldstein & Poretz (1986) "The Lectins"; Agrawal & Goldstein (1967) Biochim. Biophys. Acta 147:262

#### 1.1.3 Octyl (octylamine, strong HIC)

- **Chemistry:** n-Octylamine coupling via epoxide ring-opening. 8-carbon chain provides strongest hydrophobic interaction in the series (Phenyl < Butyl << Octyl).
- **MW:** 129 Da; r_h ~ 0.4 nm
- **k_forward:** 4e-5 m^3/(mol*s) (same class as butylamine)
- **E_a:** 48,000 J/mol
- **Coupling:** Standard epoxide-amine, pH 9.5, 25C, 4h
- **Concern:** Risk of irreversible protein binding at high ligand densities. The simulator should flag when Octyl density exceeds 20 umol/mL resin.
- **binding_model_hint:** "salt_promoted"
- **References:** Hjerten (1973) J. Chromatogr. 87:325

#### 1.1.4 WGA (wheat germ agglutinin, 36 kDa dimer)

- **Biology:** Lectin binding N-acetylglucosamine (GlcNAc) and sialic acid residues. Dimer of 18 kDa subunits.
- **MW:** 36 kDa (dimer); r_h ~ 2.3 nm
- **Coupling chemistry:** Epoxide-amine at 4C, 16h (same as protein coupling)
- **k_forward:** 4e-7 m^3/(mol*s)
- **E_a:** 25,000 J/mol
- **activity_retention:** 0.50 +/- 0.15
- **max_surface_density:** 3e-8 mol/m^2
- **Spacer:** DADPA (13A) sufficient; multiplier 1.20
- **References:** Nagata & Burger (1974) J. Biol. Chem. 249:3116

#### 1.1.5 EDA (ethylenediamine, 3A spacer)

- **Chemistry:** H2N-CH2-CH2-NH2 (MW 60). Shortest diamine spacer. Couples via one NH2 to epoxide, leaving distal -NH2 free.
- **Spacer length:** ~3A (2 carbon atoms)
- **k_forward:** 8e-5 m^3/(mol*s) (small, fast-diffusing diamine; higher than DADPA)
- **E_a:** 45,000 J/mol
- **Coupling:** pH 10-11, 25C, 4h
- **Role in v5.8:** Under Phase 1 (multiplier-only), EDA is a spacer profile (`reaction_type="spacer"`). Under Phase 2 (SPACER_ARM step), EDA becomes a SPACER_ARM reagent that consumes EPOXIDE and produces AMINE_DISTAL.
- **References:** Sundberg & Porath (1974)

#### 1.1.6 PEG-diamine Mn 600 (35A spacer)

- **Chemistry:** NH2-PEG-NH2 with average MW ~600 Da (~13 ethylene oxide units). Provides 35A flexible hydrophilic spacer.
- **k_forward:** 2e-5 m^3/(mol*s) (slower than EDA due to larger size and PEG drag)
- **E_a:** 45,000 J/mol
- **Coupling:** pH 10-11, 25C, 6h (longer time needed for PEG diffusion into pores)
- **CAS:** 929-59-9 (PEG-diamine, generic)
- **Critical for:** Con A coupling (see 1.1.2). Required spacer for any large (>50 kDa) lectin.
- **References:** Hearn & Bethell (1981) J. Chromatogr. 218:509

#### 1.1.7 BDGE (1,4-butanediol diglycidyl ether, 18A spacer)

- **Chemistry:** Bis-epoxide crosslinker. Reacts one epoxide with matrix -OH under alkaline conditions, leaving the second (distal) epoxide available for subsequent ligand coupling.
- **MW:** 202 Da
- **Spacer length:** ~18A (12-atom chain between epoxides)
- **k_forward:** 1.2e-5 m^3/(mol*s)
- **E_a:** 60,000 J/mol
- **Coupling:** pH 11-12, 25C, 4h (same conditions as ECH)
- **Hydrolysis:** ~5e-5 /s at pH 12 (slower than ECH because less ring strain)
- **Key distinction from ECH:** BDGE is functionally an ACTIVATION step that consumes HYDROXYL and produces EPOXIDE, but the product epoxide is presented on a long spacer arm. In the Phase 1 model this is indistinguishable from ECH activation. In Phase 2, it could be modeled as a combined ACTIVATION+SPACER step.
- **References:** Porath (1974) Methods Enzymol. 34:13; Sundberg & Porath (1974)

#### 1.1.8 TMAE (trimethylaminoethyl, strong anion via VS)

- **Chemistry:** 2-(Trimethylamino)ethyl chloride coupling via DVS-activated vinyl sulfone groups (Michael addition of TMAE amine to VS).
- **MW:** 138 Da; r_h ~ 0.4 nm
- **k_forward:** 3e-5 m^3/(mol*s) (vs-amine Michael addition)
- **E_a:** 40,000 J/mol
- **Coupling:** pH 8-9, 25C, 4h
- **Distinction from Q:** Q uses epoxide-amine path; TMAE uses VS-amine path. This gives the user a strong anion exchanger option for DVS-activated matrices where ECH activation is not available.
- **target_acs:** VINYL_SULFONE (not EPOXIDE)
- **chemistry_class:** "vs_amine"
- **charge_type:** "anion"
- **References:** Muller (1990) J. Chromatogr. 510:133

### 1.2 SM(PEG)n Scientific Validation

#### 1.2.1 The Amino-Spacer -> SM(PEG)n -> Cys-Protein Path

The complete reaction pathway is:

```
Step 1: ACTIVATION     OH + ECH -> EPOXIDE              (existing)
Step 2: SPACER_ARM     EPOXIDE + DADPA/DAH -> AMINE_DISTAL  (NEW)
Step 3: SPACER_ARM     AMINE_DISTAL + SM(PEG)n -> MALEIMIDE (NEW)
Step 4: PROTEIN_COUPLING  MALEIMIDE + Protein-Cys -> thioether (modified)
```

**Is this scientifically validated?** Yes, extensively:

1. **Thermo Fisher Scientific SM(PEG)n protocols** (Cat. 22101-22108): The canonical protocol is:
   - Activate surface with amine-reactive chemistry
   - React SM(PEG)n NHS-ester with surface amines (pH 7.2-7.5, 30 min, RT)
   - React maleimide terminus with protein thiol (pH 6.5-7.5, 2h, RT)

2. **Hermanson "Bioconjugate Techniques" (3rd ed., 2013), Chapter 5:** Describes the NHS-maleimide heterobifunctional crosslinker strategy for oriented protein immobilization via engineered Cys residues. The amino-spacer approach is explicitly recommended for solid supports lacking native amines.

3. **GE Healthcare Application Note 28-9078-88:** ECH-activated Sepharose -> diaminodipropylamine (DADPA) spacer -> SM(PEG)n -> oriented antibody coupling via hinge Cys after mild reduction.

4. **Karyakin et al. (2000) Anal. Chem. 72:3805:** Demonstrated the full ECH -> diamine -> NHS/maleimide -> enzyme-Cys path on agarose beads.

#### 1.2.2 Kinetics for SM(PEG)n Reactions

**NHS-ester + primary amine (Step 3a):**
- k_forward ~ 1e-3 m^3/(mol*s) at pH 7.4, 25C (Hermanson; NHS esters are fast)
- E_a ~ 35,000 J/mol
- Hydrolysis: k_hyd ~ 1e-3 /s at pH 7.4 (NHS ester half-life ~10 min in water; this is the dominant competing reaction)
- Optimal pH: 7.2-7.5 (above pH 8, hydrolysis dominates)
- Optimal time: 30 min (short! Must outpace hydrolysis)

**Maleimide + thiol (Step 3b -> Step 4):**
- k_forward ~ 5e-2 m^3/(mol*s) at pH 7.0 (maleimide-thiol is extremely fast, Michael addition)
- E_a ~ 20,000 J/mol
- Hydrolysis: k_hyd ~ 1e-5 /s at pH 7.0 (maleimide ring opening; slow)
- Optimal pH: 6.5-7.5 (above pH 7.5, maleimide reacts with amines non-specifically)
- Optimal time: 2h

#### 1.2.3 SM(PEG)n Variants

| Variant | PEG Units | Total Length (A) | MW (Da) | CAS |
|---------|-----------|------------------|---------|-----|
| SM(PEG)2 | 2 | ~18 | 425 | 1334179-85-1 |
| SM(PEG)4 | 4 | ~32 | 513 | 1229578-42-6 |
| SM(PEG)12 | 12 | ~60 | 865 | 1334179-86-2 |
| SM(PEG)24 | 24 | ~95 | 1393 | 1334179-87-3 |

#### 1.2.4 New ACSSiteType Values Required

| Value | Chemical Group | Created By | Consumed By |
|-------|---------------|------------|-------------|
| `AMINE_DISTAL` | -NH2 at spacer terminus | SPACER_ARM (DADPA/DAH/EDA on EPOXIDE) | SM(PEG)n NHS-ester coupling, or direct protein coupling |
| `MALEIMIDE` | Maleimide at crosslinker terminus | SPACER_ARM (SM(PEG)n on AMINE_DISTAL) | Thiol-protein coupling |

**CARBOXYL_DISTAL** is NOT needed in v5.8. The AHA spacer produces -COOH but that requires EDC/NHS activation, which is a separate chemistry path not in scope.

---

## PART 2 -- ARCHITECT

### 2.1 Core Design Decision: SPACER_ARM as New ModificationStepType

**Decision: YES, SPACER_ARM must be a new step type, not a subtype of ACTIVATION.**

Rationale:
1. ACTIVATION consumes a native surface group (HYDROXYL) and produces an activated group (EPOXIDE). It operates on the `activated_sites` counter.
2. SPACER_ARM consumes an activated or intermediate group and produces a NEW intermediate group on the SAME bead. The key semantic difference: the product is a new ACSSiteType entry in the acs_state dictionary, not a promotion within the same ACSProfile.
3. Modeling SPACER_ARM as ACTIVATION would require ACTIVATION to create new dictionary entries in acs_state, violating the current invariant that ACTIVATION only modifies the product_acs field within the existing profile hierarchy.
4. The ODE solver integration is identical (both route through `solve_second_order_consumption`), so code reuse is high.

### 2.2 ACS State Model Changes

#### 2.2.1 New ACSSiteType Enum Values

```python
class ACSSiteType(Enum):
    AMINE_PRIMARY = "amine_primary"
    HYDROXYL = "hydroxyl"
    CARBOXYL = "carboxyl"
    EPOXIDE = "epoxide"
    ALDEHYDE = "aldehyde"
    VINYL_SULFONE = "vinyl_sulfone"
    # v5.8 additions:
    AMINE_DISTAL = "amine_distal"     # -NH2 at spacer terminus
    MALEIMIDE = "maleimide"           # Maleimide at crosslinker terminus
```

#### 2.2.2 New ACSProfile Terminal State

**No new terminal states needed.** The existing terminal states (crosslinked, activated_consumed, hydrolyzed, ligand_coupled, blocked) are sufficient:

- AMINE_DISTAL sites consumed by SM(PEG)n -> `ligand_coupled_sites` on the AMINE_DISTAL profile (the "ligand" is the crosslinker)
- MALEIMIDE sites consumed by protein-Cys -> `ligand_coupled_sites` on the MALEIMIDE profile
- MALEIMIDE hydrolysis (ring opening) -> `hydrolyzed_sites` on the MALEIMIDE profile

The key insight is that each ACSProfile tracks ONE site type independently. When SPACER_ARM consumes EPOXIDE sites and creates AMINE_DISTAL sites, the EPOXIDE profile records `ligand_coupled_sites += consumed` and a NEW ACSProfile(site_type=AMINE_DISTAL) is inserted into acs_state with `accessible_sites = consumed * coupling_efficiency`.

#### 2.2.3 ACSProfile Initialization for Intermediate Types

New intermediate profiles (AMINE_DISTAL, MALEIMIDE) are NOT initialized from M1 export. They are created dynamically by SPACER_ARM steps:

```python
# In _solve_spacer_arm_step():
new_profile = ACSProfile(
    site_type=step.product_acs,
    total_sites=sites_created,
    accessible_sites=sites_created,  # All created sites are accessible
    activated_sites=sites_created,   # All are "activated" (ready for next step)
)
acs_state[step.product_acs] = new_profile
```

### 2.3 ModificationStepType Extension

```python
class ModificationStepType(Enum):
    SECONDARY_CROSSLINKING = "secondary_crosslinking"
    ACTIVATION = "activation"
    LIGAND_COUPLING = "ligand_coupling"
    PROTEIN_COUPLING = "protein_coupling"
    QUENCHING = "quenching"
    # v5.8:
    SPACER_ARM = "spacer_arm"
```

### 2.4 Workflow Validation Changes

The `_validate_workflow_ordering()` function needs these updates:

1. **New entry in `_STEP_ALLOWED_REACTION_TYPES`:**
   ```python
   ModificationStepType.SPACER_ARM: {"spacer_arm"},
   ```

2. **Rule 1 extension:** SPACER_ARM steps require either (a) existing target ACS in acs_state with remaining sites, or (b) a prior ACTIVATION or SPACER_ARM step that produces the target ACS type.

3. **Rule 5 (NEW):** SPACER_ARM product_acs must not already exist in acs_state (no overwriting existing profiles). This prevents ambiguous multi-source intermediate states.

4. **Rule 6 (NEW):** SM(PEG)n SPACER_ARM steps targeting AMINE_DISTAL must be preceded by a SPACER_ARM step producing AMINE_DISTAL (not just any amine source).

### 2.5 ReagentProfile Changes

New `reaction_type` values:
- `"spacer_arm"` -- for DADPA/DAH/EDA when used as SPACER_ARM steps (distinct from their Phase 1 `"spacer"` metadata role)
- `"heterobifunctional"` -- for SM(PEG)n profiles

New `chemistry_class` values:
- `"epoxide_amine_spacer"` -- diamine spacer coupling to epoxide (EPOXIDE -> AMINE_DISTAL)
- `"nhs_amine"` -- NHS-ester to amine (AMINE_DISTAL -> MALEIMIDE via SM(PEG)n)
- `"maleimide_thiol"` -- maleimide to thiol (MALEIMIDE -> thioether)

New `functional_mode` values:
- `"spacer"` (already exists conceptually)
- `"heterobifunctional_crosslinker"`

### 2.6 Solver: `_solve_spacer_arm_step()`

The new solver follows the same pattern as `_solve_activation_step()`:

1. Read remaining sites on target_acs (e.g., remaining EPOXIDE sites, or remaining AMINE_DISTAL sites)
2. Call `solve_second_order_consumption()` with spacer reagent kinetics
3. Compute sites consumed on the target profile
4. Create or update the product ACS profile in acs_state
5. The consumed sites on the source profile go to `ligand_coupled_sites` (they are "coupled" to the spacer)

**Critical difference from ACTIVATION:** The product is a new dictionary entry, not a field on the same profile. ACTIVATION updates `activated_sites` on a product profile that may already exist; SPACER_ARM creates a fresh ACSProfile.

### 2.7 File-Level Change Estimates

| File | Changes | LOC Est. |
|------|---------|----------|
| `acs.py` | Add AMINE_DISTAL, MALEIMIDE to ACSSiteType enum | +5 |
| `reagent_profiles.py` | Add 8 deferred profiles + 4 SM(PEG)n profiles + 3 spacer-arm profiles (DADPA/DAH/EDA as spacer_arm reaction_type) | +350 |
| `modification_steps.py` | Add SPACER_ARM to enum; add `_solve_spacer_arm_step()`; extend dispatch; add reaction_type "spacer_arm"/"heterobifunctional" to allowed types | +120 |
| `orchestrator.py` | Extend `_validate_workflow_ordering()` with Rules 5-6; extend `_STEP_ALLOWED_REACTION_TYPES` | +40 |
| `reactions.py` | No changes (reuses `solve_second_order_consumption`) | 0 |
| Tests: `test_modification_steps.py` | New tests for SPACER_ARM step, SM(PEG)n path, validation rules | +200 |
| Tests: `test_reagent_profiles.py` | Profile completeness checks for 12 new profiles | +60 |
| Tests: `test_integration_smpeg.py` | End-to-end: ECH -> DADPA -> SM(PEG)4 -> Protein-Cys | +100 |
| **Total** | | **~875** |

### 2.8 Existing Workflow Integration

The amino-spacer path integrates naturally with existing workflows:

**Existing workflow (v5.7):**
```
ECH(ACTIVATION) -> DEAE(LIGAND_COUPLING) -> Ethanolamine(QUENCHING)
```

**New SM(PEG)n workflow (v5.8):**
```
ECH(ACTIVATION) -> DADPA(SPACER_ARM) -> SM_PEG4(SPACER_ARM) -> Protein-Cys(PROTEIN_COUPLING) -> Ethanolamine(QUENCHING)
```

The orchestrator executes steps sequentially. The SPACER_ARM steps create intermediate ACS profiles that subsequent steps consume. No changes to the sequential execution model are needed.

**Backward compatibility:** All v5.7 workflows continue to work unchanged. The new step type is opt-in.

---

## PART 3 -- DEV-ORCHESTRATOR

### 3.1 Candidate Classification: Simple Profile vs Architecture-Dependent

#### Tier A: Simple profile additions (no arch change needed)

These candidates map directly onto existing step types with existing ACSSiteType targets:

| # | Candidate | Step Type | target_acs | Blocked By |
|---|-----------|-----------|------------|------------|
| 1 | Protein L | PROTEIN_COUPLING | EPOXIDE | Nothing |
| 3 | Octyl | LIGAND_COUPLING | EPOXIDE | Nothing |
| 4 | WGA | PROTEIN_COUPLING | EPOXIDE | Nothing |
| 8 | TMAE | LIGAND_COUPLING | VINYL_SULFONE | Nothing |

These 4 profiles can be added immediately with zero architecture changes.

#### Tier B: Simple profile + Phase 1 spacer (no arch change)

| # | Candidate | Role | Blocked By |
|---|-----------|------|------------|
| 2 | Concanavalin A | PROTEIN_COUPLING with spacer multiplier | PEG-diamine Mn 600 spacer profile (Tier B) |
| 6 | PEG-diamine Mn 600 | Spacer profile (reaction_type="spacer") | Nothing |

Con A can be added as a Phase 1 multiplier profile today, using PEG-diamine as a spacer reference. This works within the existing spacer-as-multiplier model. The Phase 2 SPACER_ARM version would improve fidelity later.

#### Tier C: Require SPACER_ARM step type (architecture change)

| # | Candidate | Role | Needs |
|---|-----------|------|-------|
| 5 | EDA | SPACER_ARM reagent | SPACER_ARM step type, AMINE_DISTAL site type |
| 7 | BDGE | ACTIVATION (Phase 1) / SPACER_ARM (Phase 2) | For Phase 1: simple ACTIVATION profile. For Phase 2: SPACER_ARM. |

EDA as a SPACER_ARM reagent is part of the SM(PEG)n capability. BDGE can be added as an ACTIVATION profile now (functionally equivalent to ECH with a longer spacer) and optionally upgraded to SPACER_ARM later.

#### SM(PEG)n profiles (fully dependent on architecture)

| Profile | PEG units | Needs |
|---------|-----------|-------|
| SM(PEG)2 | 2 | SPACER_ARM + AMINE_DISTAL + MALEIMIDE |
| SM(PEG)4 | 4 | SPACER_ARM + AMINE_DISTAL + MALEIMIDE |
| SM(PEG)12 | 12 | SPACER_ARM + AMINE_DISTAL + MALEIMIDE |
| SM(PEG)24 | 24 | SPACER_ARM + AMINE_DISTAL + MALEIMIDE |

### 3.2 Work Node Breakdown

```
WN-0: Tier A profiles (no dependencies)          Sonnet | 1h
  0.1  Add Protein L profile to reagent_profiles.py
  0.2  Add Octyl profile
  0.3  Add WGA profile
  0.4  Add TMAE profile (target_acs=VINYL_SULFONE)
  0.5  Add PEG-diamine Mn 600 spacer profile
  0.6  Add Concanavalin A profile (with PEG spacer multiplier)
  0.7  Add BDGE as ACTIVATION profile (HYDROXYL -> EPOXIDE, spacer_length=18)
  0.8  Unit tests for all 7 new profiles
  Deliverable: 7 new profiles (total 32), all existing tests pass

WN-1: ACSSiteType extension                      Sonnet | 0.5h
  1.1  Add AMINE_DISTAL, MALEIMIDE to ACSSiteType enum
  1.2  Update any downstream enum-dependent code (serialization, UI mapping)
  1.3  Unit tests for new enum values
  Deliverable: Extended enum, no functional change yet

WN-2: SPACER_ARM step type + solver              Opus | 2h
  2.1  Add SPACER_ARM to ModificationStepType enum
  2.2  Implement _solve_spacer_arm_step() in modification_steps.py
  2.3  Add dispatch case in solve_modification_step()
  2.4  Add "spacer_arm" to _STEP_ALLOWED_RTYPES
  2.5  Unit tests: SPACER_ARM consumes EPOXIDE, creates AMINE_DISTAL
  2.6  Unit tests: SPACER_ARM consumes AMINE_DISTAL, creates MALEIMIDE
  2.7  Conservation law tests: terminal-sum integrity after SPACER_ARM
  Depends on: WN-1
  Deliverable: Working SPACER_ARM solver with full test coverage

WN-3: Workflow validation extension               Opus | 1h
  3.1  Extend _validate_workflow_ordering() Rule 1 for SPACER_ARM
  3.2  Add Rule 5: product_acs must not pre-exist in acs_state
  3.3  Add Rule 6: SM(PEG)n requires prior AMINE_DISTAL producer
  3.4  Validation tests: legal and illegal SPACER_ARM orderings
  Depends on: WN-2
  Deliverable: Validated workflow rules for multi-step spacer paths

WN-4: Spacer-arm reagent profiles                Sonnet | 1.5h
  4.1  Add EDA spacer-arm profile (reaction_type="spacer_arm")
  4.2  Reclassify DADPA, DAH as dual-role: keep spacer metadata + add spacer_arm variant
       (e.g., dadpa_spacer_arm with reaction_type="spacer_arm")
  4.3  Add SM(PEG)2, SM(PEG)4, SM(PEG)12, SM(PEG)24 profiles
  4.4  Profile unit tests
  Depends on: WN-1 (needs AMINE_DISTAL, MALEIMIDE enum values)
  Deliverable: 7 new profiles (spacer-arm + SM(PEG)n variants)

WN-5: Integration testing                        Opus | 1.5h
  5.1  End-to-end test: ECH -> DADPA(SPACER_ARM) -> SM(PEG)4(SPACER_ARM) -> Cys-Protein
  5.2  End-to-end test: ECH -> EDA(SPACER_ARM) -> SM(PEG)2(SPACER_ARM) -> Cys-Protein
  5.3  Backward compat test: all v5.7 workflows still pass
  5.4  Conservation law audit: terminal sums across all ACS profiles
  5.5  Edge case: SM(PEG)n without prior spacer (should fail validation)
  5.6  Edge case: Double SPACER_ARM on same product_acs (should fail Rule 5)
  Depends on: WN-3, WN-4
  Deliverable: Full integration test suite; v5.8 release candidate

WN-6: Documentation + UI hints                   Sonnet | 0.5h
  6.1  Update REAGENT_PROFILES docstring with final count
  6.2  Add workflow examples to module docstring
  6.3  Update any UI reagent-selection filtering for new step type
  Depends on: WN-5
  Deliverable: Docs current, UI-ready
```

### 3.3 Dependency Graph

```
WN-0 ─────────────────────────────> (independent, can start immediately)

WN-1 ──> WN-2 ──> WN-3 ──> WN-5 ──> WN-6
              \              /
               WN-4 ────────
```

WN-0 is fully independent and can execute in parallel with the WN-1->WN-5 chain.

### 3.4 Model Tier Assignment

| Work Node | Model | Rationale |
|-----------|-------|-----------|
| WN-0 | Sonnet | Mechanical profile additions, pattern established |
| WN-1 | Sonnet | Enum extension, minimal logic |
| WN-2 | Opus | New solver with ODE integration, ACS state creation logic |
| WN-3 | Opus | Workflow validation rules require understanding full state machine |
| WN-4 | Sonnet | Profile additions following WN-0 pattern |
| WN-5 | Opus | Integration testing requires full system understanding |
| WN-6 | Sonnet | Documentation updates |

### 3.5 Total Effort Estimate

| Category | LOC | Time |
|----------|-----|------|
| Tier A profiles (WN-0) | ~280 | 1h |
| Enum extension (WN-1) | ~10 | 0.5h |
| SPACER_ARM solver (WN-2) | ~120 | 2h |
| Validation rules (WN-3) | ~50 | 1h |
| SM(PEG)n profiles (WN-4) | ~200 | 1.5h |
| Integration tests (WN-5) | ~200 | 1.5h |
| Docs + UI (WN-6) | ~15 | 0.5h |
| **Total** | **~875** | **~8h** |

### 3.6 Risk Register

| Risk | Impact | Mitigation |
|------|--------|------------|
| NHS-ester hydrolysis competes with coupling (k_hyd ~ k_coupling) | SM(PEG)n step may show low conversion at pH > 7.5 | Enforce pH < 7.5 validity window; warn user if conversion < 30% |
| MALEIMIDE -> amine cross-reaction above pH 7.5 | Non-specific coupling reduces site-directed benefit | Enforce pH < 7.5 for MALEIMIDE protein coupling |
| SPACER_ARM creating dynamic ACS profiles breaks serialization | FunctionalMicrosphere export may miss intermediate profiles | Ensure all ACS profiles (including dynamically created) are included in acs_profiles dict passed to FMC |
| Backward compat: old workflows referencing spacer profiles | Phase 1 spacer profiles (reaction_type="spacer") must remain filterable | Keep both variants: `dadpa_spacer` (Phase 1) and `dadpa_spacer_arm` (Phase 2) |

### 3.7 Profile Count Summary

| Version | Profiles | Step Types |
|---------|----------|------------|
| v5.7 (current) | 25 | 5 |
| v5.8 after WN-0 | 32 | 5 |
| v5.8 after WN-4 | 39 | 6 |
| Breakdown: 25 existing + 7 Tier A/B + 1 EDA-spacer-arm + 3 diamine-spacer-arm variants + 4 SM(PEG)n = 40 | | |

Note: The 3 existing spacer profiles (DADPA, AHA, DAH with reaction_type="spacer") are retained as-is. The spacer-arm variants (reaction_type="spacer_arm") are separate profile entries to maintain backward compatibility. Final count: 25 + 7 + 7 = 39 profiles (or 40 if BDGE is also given a SPACER_ARM variant in a later minor release).

---

## PART 4 -- IMPLEMENTATION SEQUENCE RECOMMENDATION

**Phase 1 (immediate, no arch change):** Execute WN-0. This adds 7 profiles with zero risk. Ship as v5.8.0-rc1.

**Phase 2 (SPACER_ARM architecture):** Execute WN-1 through WN-5 sequentially. Ship as v5.8.0-rc2.

**Phase 3 (polish):** Execute WN-6. Tag v5.8.0 release.

This two-phase approach means v5.8 can ship the deferred profiles immediately while the SPACER_ARM architecture is under development, reducing time-to-value for users who need Protein L, Octyl, WGA, TMAE, or Con A.
