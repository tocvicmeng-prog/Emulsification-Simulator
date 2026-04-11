# EmulSim Module 2 Expansion Work Plan

## Ligand Coupling, Protein Coupling, and Quenching Implementation

**Version:** 1.0
**Date:** 2026-04-12
**Status:** Approved by Scientific Advisor + Architect + Dev-Orchestrator
**Scope:** 3 new ModificationStepType workflows + 10 new reagent profiles + UI integration

---

## 1. Executive Summary

Module 2 currently implements 2 of 5 planned functionalization step types:
- **SECONDARY_CROSSLINKING** (amine-reactive: genipin, glutaraldehyde)
- **ACTIVATION** (hydroxyl-reactive: ECH -> epoxide, DVS -> vinyl sulfone)

This plan adds the remaining 3:
- **LIGAND_COUPLING** — small-molecule ligand immobilisation for IEX, IMAC, HIC
- **PROTEIN_COUPLING** — macromolecular immobilisation with steric blocking for affinity chromatography
- **QUENCHING** — terminal blocking of unreacted activated sites

Together, these complete the bead functionalization pipeline: Activation -> Coupling -> Quenching, enabling full chromatography media design simulation.

---

## 2. Scientific Assessment (Scientific Advisor)

### 2.1 Existing ODE Templates — Validation

| Template | Location | Intended Use | Verdict |
|---|---|---|---|
| Template 2: `_competitive_hydrolysis_rhs` | reactions.py:207-241 | Ligand coupling with competing hydrolysis | **CORRECT** — captures epoxide/VS coupling + aqueous hydrolysis competition |
| Template 3: `_steric_binding_rhs` | reactions.py:246-283 | Protein coupling with steric blocking | **CORRECT** — simplified RSA model; `max_sites` = steric jamming limit (~55% of geometric max) |
| `solve_second_order_consumption` | reactions.py:102-202 | Quenching (irreversible blocking) | **CORRECT** — high [reagent] drives >95% conversion |

### 2.2 Ligand Coupling — Scientific Design

**Chemistry:** Nucleophilic ring-opening (epoxide) or Michael addition (vinyl sulfone) by ligand amine/thiol groups, competing with aqueous hydrolysis.

**Priority ligands (4 profiles):**

| Ligand | Chromatography Mode | Target Site | CAS | MW (Da) | r_h (nm) | k_couple (m3/(mol*s)) | E_a (J/mol) | k_hydrol (1/s) | pH_opt | T (K) | t (s) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| DEAE | Weak anion exchange | EPOXIDE | 100-36-7 | 116 | ~0.4 | ~5e-5 | 50,000 | 1e-5 | 10-11 | 298.15 | 14,400 |
| IDA | IMAC | EPOXIDE | 142-73-4 | 133 | ~0.4 | ~2e-5 | 45,000 | 1e-5 | 10-11 | 298.15 | 21,600 |
| Phenylamine | HIC | EPOXIDE | 62-53-3 | 93 | ~0.3 | ~3e-5 | 48,000 | 5e-6 | 9-10 | 298.15 | 14,400 |
| Sulfopropyl | Strong cation exchange | EPOXIDE | 3680-02-2 | 122 | ~0.4 | ~4e-5 | 50,000 | 1e-5 | 10-11 | 298.15 | 14,400 |

**Key physics:**
- Small-molecule ligands use `reagent_accessible_area` (full pore access)
- Coupling yield = k_couple*[ligand] / (k_couple*[ligand] + k_hydrol)
- Product: `ligand_coupled_sites` = `ligand_functional_sites` (activity_retention = 1.0 for small molecules)

### 2.3 Protein Coupling — Scientific Design

**Chemistry:** Same activated site chemistry but with macromolecular steric constraints.

**Priority proteins (2 profiles):**

| Protein | Application | Target Site | MW (kDa) | r_h (nm) | k_couple (m3/(mol*s)) | E_a (J/mol) | Activity Retention | Max Surface Density (mol/m2) |
|---|---|---|---|---|---|---|---|---|
| Protein A | IgG affinity | EPOXIDE | 42 | 2.5 | ~5e-7 | 25,000 | 0.60 | ~2e-8 |
| Protein G | IgG affinity (broad) | EPOXIDE | 22 | 2.0 | ~5e-7 | 25,000 | 0.65 | ~3e-8 |

**Key physics:**
- Macromolecular ligands use `ligand_accessible_area` (pore-excluded)
- Steric blocking: f_steric = max(1 - coupled/max_sites, 0) per RSA
- `max_sites` = steric jamming limit (~55% of geometric maximum density)
- `ligand_functional_sites` = `ligand_coupled_sites` * `activity_retention`
- Activity loss from: random orientation, multi-point attachment, partial denaturation
- Coupling at 4 C (277 K) to preserve protein folding

### 2.4 Quenching — Scientific Design

**Chemistry:** High-concentration blocking reagent caps remaining activated sites.

**Priority reagents (4 profiles):**

| Reagent | Target Site | Product | CAS | k_quench (m3/(mol*s)) | E_a (J/mol) | Default [C] (mol/m3) | pH_opt | T (K) | t (s) |
|---|---|---|---|---|---|---|---|---|---|
| Ethanolamine | EPOXIDE | beta-hydroxyethylamine | 141-43-5 | ~1e-3 | 30,000 | 1000 | 8-9 | 298.15 | 7,200 |
| 2-Mercaptoethanol | VINYL_SULFONE | Thioether | 60-24-2 | ~5e-3 | 25,000 | 100 | 6-7 | 298.15 | 3,600 |
| NaBH4 | ALDEHYDE | Alcohol (-CH2OH) | 16940-66-2 | ~1e-1 | 15,000 | 50 | 7-8 | 298.15 | 1,800 |
| Acetic anhydride | AMINE_PRIMARY | Acetamide | 108-24-7 | ~5e-3 | 25,000 | 500 | 7-8 | 298.15 | 3,600 |

**Key physics:**
- High reagent concentrations (100-1000 mol/m3) drive >95% conversion
- ACS update: `blocked_sites += consumed`, no G_DN change, no product sites
- NaBH4 additionally stabilises Schiff base linkages (imine -> secondary amine)

### 2.5 Workflow Ordering Constraints

```
VALID:   Activation -> LigandCoupling -> Quenching
VALID:   Activation -> ProteinCoupling -> Quenching
VALID:   SecondaryCrosslinking -> Activation -> Coupling -> Quenching
INVALID: Coupling before Activation (no activated sites exist)
INVALID: Any step after Quenching (sites are blocked)
INVALID: Duplicate Quenching on same site type
```

**Enforcement rules for validator:**
1. LIGAND_COUPLING/PROTEIN_COUPLING requires `activated_sites > 0` on target ACS type
2. QUENCHING must be the last step (warn if followed by any step)
3. QUENCHING on a site type with `remaining_sites == 0` is a no-op (warn)

---

## 3. Architecture Design (Architect)

### 3.1 Module Dependency Graph

```
ReagentProfile (dataclass)          <- Add 5 new fields
    |
    v
reagent_profiles.py (REAGENT_PROFILES dict)   <- Add 10 new entries
    |
    v
reactions.py (ODE templates)        <- Wire Template 2 + Template 3 (already exist)
    |
    v
modification_steps.py               <- Add 3 new solver functions + dispatch
    |
    v
ui_validators.py                    <- Add M2 ordering validation rules
    |
    v
app.py (M2 tab)                     <- Add new step types to UI dropdowns
```

### 3.2 Data Structure Changes

#### ReagentProfile — New Fields (reagent_profiles.py)

```python
@dataclass
class ReagentProfile:
    # ... existing fields ...
    ligand_mw: float = 0.0            # [Da] molecular weight of ligand/protein
    ligand_r_h: float = 0.5e-9        # [m] hydrodynamic radius
    is_macromolecule: bool = False     # True -> use ligand_accessible_area
    activity_retention: float = 1.0   # Fraction retaining biological activity [0,1]
    max_surface_density: float = 0.0  # [mol/m2] steric jamming limit (0 = unlimited)
```

#### ACS Update Logic by Step Type

| Step Type | Consumes | Increments | G_DN Effect |
|---|---|---|---|
| SECONDARY_CROSSLINKING | target.consumed_sites | (none) | delta_G > 0 |
| ACTIVATION | target.consumed_sites | product.accessible_sites | 0 |
| LIGAND_COUPLING | target.consumed_sites | target.ligand_coupled_sites, target.ligand_functional_sites | 0 |
| PROTEIN_COUPLING | target.consumed_sites | target.ligand_coupled_sites, target.ligand_functional_sites (x activity_retention) | 0 |
| QUENCHING | target.consumed_sites | target.blocked_sites | 0 |

### 3.3 File Changes Summary

| File | Changes | Est. Lines |
|---|---|---|
| `reagent_profiles.py` | Add 5 fields to ReagentProfile + 10 new profiles | +180 |
| `modification_steps.py` | Add 3 solver functions + dispatch clauses | +200 |
| `reactions.py` | Wire `solve_competitive_hydrolysis()` + `solve_steric_binding()` wrappers | +80 |
| `ui_validators.py` | Add M2 ordering rules (3 new rules) | +40 |
| `app.py` | Add new step types + reagents to M2 tab dropdowns | +30 |
| `__init__.py` | Export new functions | +5 |
| **Total** | | **~535** |

---

## 4. Implementation Plan (Dev-Orchestrator)

### 4.1 Build Order

Modules must be built in dependency order. Each module goes through the full inner loop:
Protocol -> Implementation -> Audit -> Approval.

| Phase | Module | Depends On | Model Tier | Est. LOC |
|---|---|---|---|---|
| **A** | ReagentProfile expansion + 10 new profiles | None | Sonnet | 180 |
| **B** | reactions.py ODE wrappers | Phase A | Sonnet | 80 |
| **C** | `_solve_ligand_coupling_step()` | A, B | Sonnet | 80 |
| **D** | `_solve_protein_coupling_step()` | A, B | Opus (novel steric model) | 80 |
| **E** | `_solve_quenching_step()` | A | Haiku (simple) | 40 |
| **F** | Dispatch integration in `solve_modification_step()` | C, D, E | Sonnet | 30 |
| **G** | M2 ordering validation rules | F | Sonnet | 40 |
| **H** | UI integration (M2 tab dropdowns) | F, G | Sonnet | 30 |
| **I** | Integration testing (full M1->M2->M3 pipeline) | All | Opus | 50 |

### 4.2 Phase Details

#### Phase A: ReagentProfile Expansion

**Purpose:** Extend ReagentProfile dataclass with 5 new fields and add 10 new reagent profiles.

**Input:** Scientific Advisor reagent specifications (Section 2)
**Output:** Updated `reagent_profiles.py` with 14 total profiles (4 existing + 10 new)

**Test cases:**
- All 10 new profiles instantiate without error
- `is_macromolecule = True` only for Protein A, Protein G
- `activity_retention < 1.0` only for protein profiles
- `max_surface_density > 0` only for protein profiles
- All CAS numbers are non-empty strings

#### Phase B: ODE Wrapper Functions

**Purpose:** Create public solver functions that wrap Template 2 and Template 3.

**Functions:**
- `solve_competitive_hydrolysis(acs_conc, reagent_conc, k_couple, k_hydrol, stoich, time, temperature, E_a, k0)` -> `(conversion, reagent_remaining, hydrolysis_fraction)`
- `solve_steric_binding(acs_conc, ligand_conc, k_couple, max_sites, time, temperature, E_a, k0)` -> `(conversion, ligand_remaining)`

**Test cases:**
- Zero ligand concentration -> conversion = 0
- Very high k_couple, long time -> conversion approaches 1.0
- Steric model: conversion never exceeds max_sites/initial_sites
- Hydrolysis model: conversion < 1.0 even at infinite time (hydrolysis competes)

#### Phase C: Ligand Coupling Solver

**Purpose:** Implement `_solve_ligand_coupling_step()` using Template 2 ODE.

**ACS updates:**
1. Compute conversion via `solve_competitive_hydrolysis()`
2. `sites_consumed = conversion * target_profile.remaining_sites`
3. `target_profile.consumed_sites += sites_consumed`
4. `target_profile.ligand_coupled_sites += sites_consumed` (1:1 for small molecules)
5. `target_profile.ligand_functional_sites += sites_consumed * activity_retention`
6. `delta_G = 0` (no mechanical effect)

**Test cases:**
- DEAE coupling on epoxide: conversion > 0, ligand_coupled_sites > 0
- No activated sites -> conversion = 0
- Hydrolysis fraction increases at higher k_hydrol
- Conservation: consumed + remaining = original accessible

#### Phase D: Protein Coupling Solver

**Purpose:** Implement `_solve_protein_coupling_step()` using Template 3 ODE.

**ACS updates:**
1. Compute conversion via `solve_steric_binding()`
2. `sites_consumed = conversion * target_profile.remaining_sites`
3. `target_profile.consumed_sites += sites_consumed`
4. `target_profile.ligand_coupled_sites += sites_consumed`
5. `target_profile.ligand_functional_sites += sites_consumed * reagent.activity_retention`
6. `delta_G = 0`

**Steric limit enforcement:**
- `max_coupled = reagent.max_surface_density * surface_model.ligand_accessible_area`
- Coupling stops when `ligand_coupled_sites >= max_coupled`

**Test cases:**
- Protein A coupling: conversion > 0, functional < coupled (activity_retention = 0.6)
- Steric saturation: repeated coupling steps don't exceed max_surface_density
- Large protein (r_h = 5 nm) on small pores (50 nm): reduced accessible area
- Conservation checks pass

#### Phase E: Quenching Solver

**Purpose:** Implement `_solve_quenching_step()` using existing `solve_second_order_consumption()`.

**ACS updates:**
1. Compute conversion (expect >0.95 due to high [reagent])
2. `sites_consumed = conversion * target_profile.remaining_sites`
3. `target_profile.consumed_sites += sites_consumed`
4. `target_profile.blocked_sites += sites_consumed`
5. `delta_G = 0`

**Test cases:**
- Ethanolamine on epoxide: conversion > 0.95
- blocked_sites > 0 after quenching
- remaining_sites -> ~0 after quenching
- No G_DN change

#### Phase F: Dispatch Integration

**Purpose:** Add if/elif clauses in `solve_modification_step()` for the 3 new step types.

```python
if step.step_type == ModificationStepType.LIGAND_COUPLING:
    return _solve_ligand_coupling_step(step, acs_state, surface_model, reagent_profile)
elif step.step_type == ModificationStepType.PROTEIN_COUPLING:
    return _solve_protein_coupling_step(step, acs_state, surface_model, reagent_profile)
elif step.step_type == ModificationStepType.QUENCHING:
    return _solve_quenching_step(step, acs_state, surface_model, reagent_profile)
```

#### Phase G: M2 Ordering Validation

**Purpose:** Add 3 new validation rules to `ui_validators.py`.

| Rule | Check | Severity |
|---|---|---|
| M2-R12 | Coupling step requires activated_sites > 0 on target type | BLOCKER |
| M2-R13 | No steps after a Quenching step | WARNING |
| M2-R14 | No duplicate Quenching on same site type | WARNING |

#### Phase H: UI Integration

**Purpose:** Add new step types and reagent options to M2 tab in app.py.

Update the Chemistry selectbox to include:
```python
["Secondary Crosslinking", "Hydroxyl Activation",
 "Ligand Coupling", "Protein Coupling", "Quenching"]
```

Each selection dynamically populates the Reagent dropdown with the appropriate profiles.

#### Phase I: Integration Testing

**Purpose:** Run full M1->M2->M3 pipeline with new step types.

**Test workflows:**
1. `Activation(ECH) -> LigandCoupling(DEAE) -> Quenching(ethanolamine)` -> run M3 IEX breakthrough
2. `Activation(ECH) -> ProteinCoupling(ProteinA) -> Quenching(ethanolamine)` -> run M3 affinity breakthrough
3. `SecondaryCrosslinking(genipin) -> Activation(DVS) -> LigandCoupling(phenyl) -> Quenching(mercaptoethanol)` -> run M3 HIC
4. Invalid: `LigandCoupling before Activation` -> expect BLOCKER
5. Invalid: `Step after Quenching` -> expect WARNING

---

## 5. Risk Register

| Risk | Likelihood | Impact | Mitigation |
|---|---|---|---|
| Rate constants are order-of-magnitude estimates | HIGH | MEDIUM | Flag as "illustrative — user calibration required" in UI |
| Steric blocking model (linear RSA) is simplified | MEDIUM | LOW | Sufficient for ranking; note in trust tier as "semi_quantitative" |
| Protein activity retention varies widely in practice | HIGH | MEDIUM | Expose as user-adjustable parameter with literature default |
| Template 2/3 ODE stiffness at extreme parameters | LOW | HIGH | Use Radau solver (already used); add convergence check |
| UI widget key collisions with new dropdown options | LOW | HIGH | Prefix all new keys with `m2_` (already established) |
| M3 isotherm parameters don't match new ligand types | MEDIUM | MEDIUM | Document that M3 isotherm must be recalibrated for each ligand |

---

## 6. Milestone Definition

**Milestone:** M2 Complete Functionalization Pipeline

**Acceptance criteria:**
- [ ] All 5 ModificationStepTypes implemented and dispatching correctly
- [ ] All 14 reagent profiles (4 existing + 10 new) instantiate and validate
- [ ] ACS conservation checks pass for all step type combinations
- [ ] M2 ordering validation rules enforce Activation->Coupling->Quenching sequence
- [ ] UI M2 tab exposes all 5 step types with correct reagent dropdowns
- [ ] Full M1->M2->M3 integration test passes for IEX, affinity, and HIC workflows
- [ ] Trust tier labels correct: "semi_quantitative" for all new step types

---

## 7. Token Economy Estimate

| Phase | Model | Est. Tokens (input+output) |
|---|---|---|
| A: ReagentProfile expansion | Sonnet | ~4,000 |
| B: ODE wrappers | Sonnet | ~3,000 |
| C: Ligand coupling solver | Sonnet | ~4,000 |
| D: Protein coupling solver | Opus | ~5,000 |
| E: Quenching solver | Haiku | ~2,000 |
| F: Dispatch integration | Sonnet | ~2,000 |
| G: Validation rules | Sonnet | ~2,000 |
| H: UI integration | Sonnet | ~2,000 |
| I: Integration testing | Opus | ~4,000 |
| Audits (3x) | Opus | ~6,000 |
| **Total** | | **~34,000** |

---

> **Disclaimer**: Rate constants in this plan are order-of-magnitude estimates from
> general reaction chemistry principles. All values should be treated as illustrative
> defaults requiring experimental calibration for specific materials. This plan is
> provided for development purposes only and should be validated through appropriate
> testing before production use.
