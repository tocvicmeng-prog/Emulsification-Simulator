# EmulSim Second-Pass Review: Scientific Model Redesign and Module-by-Module Fix Plan

**Project:** Emulsification-Simulator  
**Date:** 2026-04-10  
**Purpose:** second-pass audit focused on:

1. redesigning the scientific model hierarchy
2. defining a concrete module-by-module fix plan

This document is a follow-on to:

- `docs/05_independent_scientific_audit_2026-04-10.md`

---

## 1. Executive Summary

The main problem in EmulSim is not that the individual formulas are all wrong. The main problem is that the **scientific hierarchy is incoherent**:

- Level 1 is presented as a mechanistic droplet-size model.
- Level 2 defaults to an empirical pore model that is largely independent of Level 1.
- Level 3 compresses chemically incompatible crosslinkers into one abstraction.
- Level 4 treats very different network constructions as if they were the same type of second network.

The result is a pipeline that is modular in software terms but not yet coherent in scientific terms.

### Core recommendation

EmulSim should be reorganized into **three clearly separated scientific operating modes**:

1. **Empirical Engineering Mode**
   - fast
   - calibration-driven
   - only modest mechanistic claims

2. **Hybrid Coupled Mode**
   - empirical where evidence exists
   - mechanistic only where coupling is defensible
   - preferred production mode

3. **Mechanistic Research Mode**
   - slower
   - explicitly exploratory
   - hypothesis-testing only

At present, the code mixes these modes implicitly. The redesign below makes them explicit.

---

## 2. Design Principles for a Scientifically Coherent EmulSim v2

## Principle 1: Each layer must have a single scientific role

Every module must answer one question only:

- L1: what droplet population is produced?
- L2: what internal microstructure forms inside a droplet of that size and thermal history?
- L3: what network topology and crosslink density emerge from a specific chemistry?
- L4: what bulk bead properties follow from that structure?

At present some layers do two jobs, some do none, and some do placeholders while presenting mechanistic authority.

## Principle 2: Empirical and mechanistic models must not be mixed silently

If a model is empirical:

- say so
- state its calibration range
- forbid interpretation outside that range

If a model is mechanistic:

- preserve state variables and scale couplings
- avoid hard-coded heuristics that override the mechanism

## Principle 3: Structural uncertainty must be explicit

The largest uncertainty in EmulSim is not parameter uncertainty. It is:

- whether the chosen hierarchy is even the right one
- whether a binary pore model is appropriate
- whether different crosslinkers can share one kinetic abstraction

That uncertainty must be surfaced, not buried.

## Principle 4: Local patch simulations must never be reused as full-bead geometry

If L2 uses a local microstructure patch for numerical feasibility, that patch must remain a local descriptor. It must not overwrite:

- bead radius
- particle contact scale
- packing mechanics scale

## Principle 5: Chemistry families must be separated

Amines, hydroxyls, ionic crosslinks, and independently polymerized PEG networks are not variations of one kinetic model. They are different scientific objects and must be treated as such.

---

## 3. Proposed Scientific Model Hierarchy v2

## 3.1 Top-Level Architecture

The revised pipeline should be:

```text
Input formulation + process conditions
        |
        v
L1: Emulsification population model
        |
        | outputs:
        | - DSD
        | - representative bead size classes
        | - interfacial area history
        v
L2a: Thermal / gelation trajectory model
        |
        | outputs:
        | - T(t) inside droplet class
        | - gelation onset / arrest timing
        v
L2b: Microstructure model
        |
        | outputs:
        | - macropore size distribution
        | - porosity
        | - morphology class / connectivity
        v
L3: Chemistry-specific network formation model
        |
        | outputs:
        | - crosslink density
        | - mesh size
        | - network type metadata
        v
L4: Property homogenization model
        |
        | outputs:
        | - modulus
        | - compression response
        | - permeability / partition metrics
```

The critical change is that Level 2 is split into:

- **L2a: thermal/gelation timing**
- **L2b: pore/microstructure formation**

This is scientifically cleaner than trying to make one solver represent everything.

---

## 3.2 Three Operating Modes

## Mode A: Empirical Engineering Mode

### Purpose

Fast screening once calibration data exist.

### Allowed model types

- L1: calibrated surrogate or semi-empirical PBE
- L2a: calibrated cooling/gelation timing
- L2b: empirical pore correlation
- L3: chemistry-family empirical kinetics
- L4: phenomenological modulus model

### Claims allowed

- trend screening
- calibrated interpolation
- design space ranking

### Claims not allowed

- first-principles pore prediction
- mechanistic explanation of polymer demixing

## Mode B: Hybrid Coupled Mode

### Purpose

Main operational mode for the project.

### Allowed model types

- L1: mechanistic or semi-mechanistic
- L2a: mechanistic thermal/gelation timing
- L2b: empirical but explicitly size-coupled
- L3: chemistry-family-specific reduced models
- L4: phenomenological homogenization

### Why this should be the default

This is the best balance between:

- runtime
- interpretability
- defensibility

## Mode C: Mechanistic Research Mode

### Purpose

Hypothesis testing for future publications and model development.

### Allowed model types

- L1: mechanistic PBE with calibrated kernels
- L2a/L2b: true multiphase field or ternary phase separation model
- L3: reaction-diffusion or spatially resolved chemistry
- L4: microstructure-aware constitutive model

### Usage note

This mode should be explicitly labeled:

- exploratory
- non-production
- calibration-dependent

---

## 4. What Each Current Layer Should Become

## 4.1 Level 1: Emulsification

### Current status

PBE framework exists, but legacy default behavior is not directionally robust and calibration is weak.

### v2 scientific role

L1 should become:

- a droplet population generator with explicit uncertainty and regime validity

### Required outputs

- `d10`, `d50`, `d90`, `d32`, `d_mode`
- class-based DSD
- representative droplet classes for downstream simulation
- regime flags:
  - inertial
  - viscous subrange
  - uncertain

### Required scientific improvements

- fit breakage and coalescence constants to real DSD measurements
- separate rotor-stator and stirred-vessel calibrations cleanly
- require monotonic sanity checks over RPM where physically expected
- expose kernel validity range, not just numerical output

---

## 4.2 Level 2a: Thermal / Gelation Timing

### Current status

Thermal handling is partly embedded in L1 and partly simplified in L2.

### v2 scientific role

L2a should answer:

- when does the droplet gel?
- how fast does gelation arrest coarsening?
- what is the relevant local cooling history?

### Why this split helps

Thermal history and gelation arrest are scientifically different from pore morphology. Separating them reduces conceptual overload and makes calibration easier.

### Recommended model content

- droplet cooling timescale
- Biot-number check for lumped vs radial cooling
- agarose gelation onset
- Avrami or better calibrated gelation arrest
- optional droplet-size dependence in heat transfer

### Minimal viable implementation

Keep the current cooling and Avrami logic, but move it into a dedicated Level 2a object that outputs:

- `T_history`
- `t_gel_onset`
- `alpha_final`
- `mobility_arrest_factor`

---

## 4.3 Level 2b: Pore / Microstructure Formation

### Current status

Two incompatible paradigms coexist:

- a fast empirical pore correlation
- a binary Cahn-Hilliard model that cannot represent agarose-chitosan demixing

### v2 scientific role

L2b should answer:

- what macropore size distribution forms?
- what morphology class results?
- how does that depend on droplet class and L2a arrest history?

### Recommended near-term approach

Do **not** attempt a full ternary mechanistic model immediately for production use.

Instead:

1. keep the empirical model
2. make it explicitly **size-coupled**
3. add morphology metadata
4. state calibration limits clearly

### Recommended long-term mechanistic approach

Implement a ternary or two-order-parameter model:

- agarose-rich field
- chitosan-rich field
- water implied by closure

This is the only way to mechanistically test the claimed polymer-polymer demixing hypothesis.

---

## 4.4 Level 3: Chemistry-Specific Network Formation

### Current status

Chemically distinct crosslinkers are over-compressed into one abstraction.

### v2 scientific role

L3 should answer:

- what network is being formed?
- where is it formed?
- what reactive groups limit it?
- is it covalent, ionic, reversible, or independent?

### Required redesign

Split L3 into four solver families:

1. `amine_covalent`
2. `hydroxyl_covalent`
3. `ionic_reversible`
4. `independent_network_polymerization`

Each family should expose:

- network target
- stoichiometric basis
- limiting reactant
- effective elastically active fraction
- diffusion-limited vs reaction-limited flag

### Why this matters

Without this split, downstream mechanics cannot distinguish:

- chitosan network strengthening
- agarose network strengthening
- PEGDA third-network formation
- weak ionic stabilization

---

## 4.5 Level 4: Property Homogenization

### Current status

Useful as a ranking layer, but scientifically overinterpretable.

### v2 scientific role

L4 should be explicitly a:

- homogenized property estimator

not a first-principles bead mechanics solver.

### Required inputs

- actual bead radius from L1
- pore descriptors from L2b
- network descriptors from L3
- morphology metadata

### Required outputs

- modulus estimate
- compression curve estimate
- partition/permeation estimate
- confidence / applicability flags

### Required conceptual correction

Never infer bead radius from a truncated L2 computational patch.

---

## 5. Concrete Module-by-Module Fix Plan

## 5.1 `src/emulsim/pipeline/orchestrator.py`

### Problems

- Default L2 path is uncoupled from L1 in substance.
- Single `solve_gelation()` call hides empirical vs mechanistic meaning.
- No explicit mode object for scientific operating mode.

### Required changes

1. Introduce a `model_mode` setting:
   - `empirical_engineering`
   - `hybrid_coupled`
   - `mechanistic_research`

2. Split current L2 call into:
   - `solve_gelation_timing()`
   - `solve_microstructure()`

3. Pass actual bead radius separately from any local microstructure domain size.

4. Allow L1 class-based sampling for L2b:
   - for example, simulate microstructure on `d10`, `d50`, `d90`
   - aggregate into bead-population summaries

### Acceptance criteria

- Changing RPM must change at least one downstream microstructure output in default coupled mode.
- `bead_radius_actual` and `microstructure_patch_size` are stored separately.

---

## 5.2 `src/emulsim/datatypes.py`

### Problems

- Current datatypes blur true state vs convenience fields.
- No place to represent model-family metadata or confidence structure.

### Required changes

Add explicit dataclasses for:

- `ModelMode`
- `GelationTimingResult`
- `MicrostructureResult`
- `NetworkTypeMetadata`
- `ApplicabilityFlags`

Add fields to results for:

- `actual_bead_radius`
- `local_microstructure_domain_size`
- `model_family`
- `calibration_range_ok`
- `structural_uncertainty_flag`

### Acceptance criteria

- Every result object states whether it came from empirical or mechanistic logic.
- No downstream solver has to guess what an upstream radius means.

---

## 5.3 `src/emulsim/properties/database.py`

### Problems

- Property file and hard-coded runtime values are not fully aligned.
- Some values that look database-driven are actually embedded in code elsewhere.

### Required changes

1. Make the property database genuinely authoritative.
2. Move hard-coded interfacial parameters, crosslinking defaults, and key constitutive constants into governed property/config data.
3. Add metadata fields:
   - source
   - calibration status
   - valid temperature range
   - valid composition range

### Acceptance criteria

- No scientifically important constant is silently overridden in code without metadata.
- User can inspect one source of truth for each calibrated parameter.

---

## 5.4 `src/emulsim/properties/interfacial.py`

### Problems

- Dynamic IFT exists but is mostly disconnected from practical prediction paths.
- Hard-coded `K_L = 0.75` conflicts with TOML metadata.

### Required changes

1. Parameterize the Szyszkowski-Langmuir model from the property database.
2. Add explicit switch:
   - `equilibrium_ift`
   - `dynamic_ift`

3. Use dynamic IFT only where breakup timescale justifies it.

### Acceptance criteria

- Documentation, TOML data, and runtime parameters agree.
- User can tell which IFT regime is being used.

---

## 5.5 `src/emulsim/properties/viscosity.py`

### Problems

- Current blend logic is simple and useful, but not yet enough for confidence in high-shear polymer emulsification.

### Required changes

1. Add a calibration slot for measured rheology curves.
2. Separate:
   - zero-shear viscosity
   - apparent viscosity at gap shear
   - viscosity-ratio validity flags

3. Add optional fitted Cross or Carreau-Yasuda parameters per formulation family.

### Acceptance criteria

- L1 can use calibrated apparent viscosity instead of generic defaults where data exist.
- Trust layer warns when viscosity model is extrapolated beyond calibration.

---

## 5.6 `src/emulsim/level1_emulsification/kernels.py`

### Problems

- Kernel choice is scientifically plausible, but default calibration is not reliable enough.
- Legacy default disables the viscous correction where it likely matters.

### Required changes

1. Formalize kernel families:
   - rotor-stator calibrated
   - stirred-vessel calibrated
   - exploratory research

2. Require calibrated `C1-C5` sets rather than one nominal default.
3. Add directionality sanity checks:
   - size vs RPM
   - size vs sigma
   - size vs phi_d

4. Track regime diagnostics explicitly:
   - `d / eta_K`
   - `mu_d / mu_c`
   - breakup regime label

### Acceptance criteria

- Default kernel set no longer produces clearly implausible RPM trends.
- Outputs include regime validity metadata.

---

## 5.7 `src/emulsim/level1_emulsification/solver.py`

### Problems

- Legacy and stirred-vessel modes are mixed into one large solver concept.
- Only a single representative size is passed downstream.

### Required changes

1. Split into:
   - `legacy_rotor_stator_solver`
   - `stirred_vessel_solver`

2. Add representative droplet class export:
   - `small`, `median`, `large`
   - or percentile-based classes

3. Add structured convergence diagnostics:
   - steady-state reached
   - oscillatory late state
   - regime mismatch

### Acceptance criteria

- Downstream modules can consume a population summary rather than only one `d50`.
- Solver summaries report scientific validity, not just numerical convergence.

---

## 5.8 `src/emulsim/level2_gelation/gelation.py`

### Problems

- Gelation kinetics are embedded in pore simulation logic rather than being a first-class layer.

### Required changes

Convert this module into Level 2a:

- `solve_gelation_timing()`
- `gelation_rate_constant()`
- `mobility_arrest_factor()`
- optional Biot-number and droplet-cooling corrections

### Acceptance criteria

- L2a outputs are reusable by both empirical and mechanistic microstructure models.

---

## 5.9 `src/emulsim/level2_gelation/solver.py`

### Problems

- Empirical and mechanistic models are mixed behind a single interface.
- Empirical model ignores droplet size in substance.
- Mechanistic model caps domain size and then passes it on incorrectly.

### Required changes

#### Empirical branch

1. Rename to something explicit, such as:
   - `solve_microstructure_empirical()`

2. Add real droplet-size coupling:
   - e.g. through local cooling time, gelation timing, or calibrated bead-size term

3. Include morphology outputs:
   - pore class
   - connectivity class
   - confidence band

#### Mechanistic branch

1. Rename to:
   - `solve_microstructure_phasefield()`

2. Store:
   - `actual_bead_radius`
   - `local_domain_size`

3. Prevent local domain size from being reused as bead radius.

4. Long-term: replace binary scalar model with ternary/two-order-parameter formulation.

### Acceptance criteria

- Default coupled mode responds to Level 1 size changes.
- Mechanistic mode no longer contaminates Level 4 bead radius.

---

## 5.10 `src/emulsim/level2_gelation/free_energy.py`

### Problems

- Present implementation is mathematically fine for a binary model, but scientifically incomplete for the stated chemistry.

### Required changes

Short term:

- keep as binary helper for exploratory mode only

Long term:

- replace with ternary free-energy library
- add polymer-polymer interaction parameter support
- allow separate agarose and chitosan order parameters

### Acceptance criteria

- No documentation claims polymer-polymer demixing while using a scalar binary model.

---

## 5.11 `src/emulsim/level2_gelation/pore_analysis.py`

### Problems

- Outputs are serviceable, but there is no morphology classification layer.

### Required changes

Add morphology descriptors:

- bicontinuous score
- anisotropy
- connectivity proxy
- chord-length skewness

These are needed because pore mean alone is not enough.

### Acceptance criteria

- L2b reports more than a single scalar pore size.

---

## 5.12 `src/emulsim/level3_crosslinking/solver.py`

### Problems

- Over-compressed chemistry abstraction
- single concentration field reused for all chemistries
- outputs not explicit about network identity

### Required changes

Split into family-specific solvers:

- `solve_amine_covalent()`
- `solve_hydroxyl_covalent()`
- `solve_ionic_reversible()`
- `solve_independent_network_polymerization()`

Introduce family-specific parameter containers:

- `c_amine_crosslinker`
- `c_hydroxyl_crosslinker`
- `c_ionic_crosslinker`
- `c_pegda`
- `uv_intensity`

Return metadata:

- `network_target = chitosan | agarose | independent | mixed`
- `bond_type = covalent | ionic | reversible`
- `is_true_second_network`

### Acceptance criteria

- Mechanical layer can distinguish what network was actually formed.
- UI and optimization no longer reuse `c_genipin` as universal chemistry concentration.

---

## 5.13 `src/emulsim/reagent_library.py`

### Problems

- Library is rich, but scientific metadata are not fully connected to model limitations.

### Required changes

Add fields such as:

- `solver_family`
- `network_target`
- `requires_alkaline_conditions`
- `reversible`
- `requires_diffusion_model`
- `confidence_level`

### Acceptance criteria

- Crosslinker selection automatically routes to the correct solver family.

---

## 5.14 `src/emulsim/level4_mechanical/solver.py`

### Problems

- Uses phenomenological modulus mixing without enough metadata.
- Can inherit incorrect radius from mechanistic L2.

### Required changes

1. Replace implicit radius inference with explicit bead radius input.
2. Rename current modulus model to indicate its status:
   - `phenomenological_dn_modulus()`

3. Add model branches:
   - `semi_ipn_estimator`
   - `triple_network_estimator`
   - `ionic_gel_estimator`

4. Make the output include:
   - property confidence score
   - constitutive model family

### Acceptance criteria

- Level 4 no longer assumes all L3 outputs mean the same network physics.

---

## 5.15 `src/emulsim/optimization/objectives.py`

### Problems

- Optimizer inherits scientific disconnects from upstream.

### Required changes

1. Prevent optimization runs in modes with broken coupling unless explicitly allowed.
2. Add penalty terms for:
   - structural uncertainty
   - regime mismatch
   - model-family incompatibility

3. Allow objective sets that differ by operating mode.

### Acceptance criteria

- Optimizer cannot silently optimize a scientifically incoherent pipeline.

---

## 5.16 `src/emulsim/uncertainty.py`

### Problems

- Only parameter uncertainty is propagated.
- Level 2 is forced to empirical mode.

### Required changes

Add two uncertainty classes:

1. `ParametricUncertaintyResult`
2. `StructuralUncertaintyResult`

Structural uncertainty should span:

- empirical vs mechanistic L2b
- alternative pore-coupling forms
- alternative L1 kernel sets
- alternative L4 constitutive closures

### Acceptance criteria

- Pore CI is no longer identically degenerate when structural uncertainty is material.
- Reports explicitly state whether they quantify parametric or structural uncertainty.

---

## 5.17 `src/emulsim/trust.py`

### Problems

- Trust checks are useful but still too permissive about structural issues.

### Required changes

Add blockers or warnings for:

- empirical L2 used in a supposedly fully mechanistic run
- local microstructure patch reused as bead radius
- crosslinker family mismatch
- chemistry outside calibration family
- uncertainty report that excludes structural uncertainty

### Acceptance criteria

- Trust output becomes a scientific validity gate, not just a numerical sanity checker.

---

## 6. Proposed Implementation Phases

## Phase 1: Scientific honesty and coupling repair

### Goal

Make the current code scientifically honest before making it more complex.

### Tasks

- add `model_mode`
- split L2 into timing and microstructure roles
- make empirical L2 size-coupled
- separate true bead radius from local patch size
- relabel Level 4 as phenomenological

### Expected result

The pipeline becomes scientifically interpretable even if still partly empirical.

## Phase 2: Chemistry-family refactor

### Goal

Remove false generality from Level 3.

### Tasks

- split solver families
- split parameter fields
- add network identity metadata
- update UI and optimization bindings

### Expected result

Crosslinker comparisons become scientifically meaningful.

## Phase 3: Calibration consolidation

### Goal

Make default mode trustworthy.

### Tasks

- calibrate L1 kernels
- calibrate empirical L2 size-coupled model
- calibrate L4 modulus estimator to compression data
- update property database governance

### Expected result

A usable `hybrid_coupled` production mode.

## Phase 4: Mechanistic research branch

### Goal

Develop the true research-grade microstructure model without destabilizing the production mode.

### Tasks

- ternary/two-field phase model
- local-to-global homogenization strategy
- reaction-diffusion extensions where justified

### Expected result

A separate research mode for hypothesis testing and future publications.

---

## 7. Recommended Immediate Deliverables

If development time is limited, the first five concrete deliverables should be:

1. Add `model_mode` and expose it in CLI/UI/config.
2. Refactor Level 2 into `timing` and `microstructure`.
3. Make empirical L2 explicitly depend on droplet size class.
4. Separate bead radius from local domain size in all result objects.
5. Split Level 3 into chemistry families and stop reusing `c_genipin` universally.

If only those five are done, the project will already become much more scientifically coherent.

---

## 8. Suggested New Repository Language

To avoid overstating the present system, the README and docs should describe the software as:

> A modular process-modeling platform combining calibrated empirical models and mechanistic submodels for emulsified hydrogel microsphere fabrication.

And not as:

> a fully predictive multiscale simulation of the complete fabrication pipeline

until the hierarchy is repaired and calibrated.

---

## 9. Final Recommendation

The project should now prioritize **scientific architecture refactoring** over adding more chemistries, UI features, or optimization options.

The correct next move is:

- first make the pipeline scientifically coherent,
- then calibrate it,
- then expand it.

In short:

**narrow the claims, repair the couplings, separate chemistry families, and only then deepen the mechanics.**

