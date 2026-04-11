# Audit Report: Module 2 Expansion Work Plan

Date: 2026-04-12

Audited plan: `docs/15_m2_expansion_workplan.md`

Audited implementation areas:

- `src/emulsim/module2_functionalization/`
- `src/emulsim/visualization/ui_validators.py`
- `src/emulsim/visualization/ui_model_metadata.py`
- `src/emulsim/visualization/app.py`
- `src/emulsim/module3_performance/`
- `tests/test_module2_acs.py`
- `tests/test_module2_workflows.py`
- `tests/test_ui_contract.py`
- `tests/test_module3_breakthrough.py`
- `tests/test_module3_multicomponent.py`

## 1. Executive Verdict

The Module 2 expansion work plan is directionally correct. Completing `Activation -> Coupling -> Quenching` is the right next scientific step if EmulSim is to simulate functional chromatography media rather than only primary hydrogel microspheres. The plan correctly recognizes that small-molecule ligand coupling, protein coupling, and quenching are distinct chemical workflows and that they must be tracked through ACS conservation.

However, the plan is not yet scientifically or computationally safe as written. Several proposed update rules are inconsistent with the current ACS state model, several chemistry assumptions are oversimplified for wetlab realism, and the M2-to-M3 bridge is under-specified. The most serious defect is the proposed quenching update: adding the same quenched sites to both `consumed_sites` and `blocked_sites` would violate the current conservation invariant `consumed_sites + blocked_sites <= accessible_sites` and can immediately make a valid ACS profile invalid at high quench conversion.

Overall assessment:

| Area | Rating | Audit conclusion |
|---|---:|---|
| Scientific direction | Good | Coupling/quenching are necessary for a real functionalization pipeline. |
| Current backend readiness | Moderate | ACS, surface area, activation, and secondary crosslinking exist; new workflows are still blocked. |
| Workplan scientific validity | Conditional | Valid only after correcting ACS bookkeeping, pH chemistry, reagent definitions, and M3 mapping. |
| Wetlab realism | Moderate-low | Useful screening defaults, but chemical identities, pH dependence, hazards, and calibration are under-specified. |
| Computational reliability | Moderate | Existing M2 tests pass; new ODE wrappers and solver diagnostics are not implemented. |
| M2->M3 consistency | Low | Functional ligand density is not yet mapped to binding capacity, isotherm constants, selectivity, or catalytic activity. |
| UI readiness | Low-medium | Current UI intentionally exposes only two M2 workflows and labels the other three as planned. |

Recommended decision:

- Proceed with the expansion only after revising the ACS update rules and chemistry registry.
- Do not expose ligand coupling, protein coupling, or quenching as runnable UI workflows until backend conservation tests, pH/temperature validity gates, and M2-to-M3 parameter mapping are implemented.
- Treat all new rate constants and functional ligand outputs as semi-quantitative and calibration-required.

## 2. Current Implementation Evidence

Current Module 2 implementation:

- `ModificationStepType` already includes `LIGAND_COUPLING`, `PROTEIN_COUPLING`, and `QUENCHING` as enum values in `src/emulsim/module2_functionalization/modification_steps.py:70`.
- `solve_modification_step()` still dispatches only `SECONDARY_CROSSLINKING` and `ACTIVATION`; other step types raise `KeyError` at `src/emulsim/module2_functionalization/modification_steps.py:201`.
- `ReagentProfile` currently has only the original fields and four profiles in `src/emulsim/module2_functionalization/reagent_profiles.py`; no ligand/protein/quench profiles exist yet.
- `_competitive_hydrolysis_rhs()` and `_steric_binding_rhs()` exist in `src/emulsim/module2_functionalization/reactions.py:207` and `src/emulsim/module2_functionalization/reactions.py:246`, but they are private templates without public wrappers.
- Existing M2 Arrhenius handling has already been improved: `_arrhenius_prefactor()` exists at `src/emulsim/module2_functionalization/modification_steps.py:46` and is used by the two current workflows.
- `ACSProfile.validate()` enforces `ligand_coupled_sites <= activated_sites`, `ligand_functional_sites <= ligand_coupled_sites`, and `consumed_sites + blocked_sites <= accessible_sites`.

Current UI and validators:

- The current M2 UI exposes only "Secondary Crosslinking" and "Hydroxyl Activation" at `src/emulsim/visualization/app.py:851`.
- The UI explicitly labels ligand coupling, protein coupling, and quenching as planned but not implemented at `src/emulsim/visualization/app.py:893`.
- `src/emulsim/visualization/ui_validators.py:51` supports only `secondary_crosslinking` and `activation`.
- `src/emulsim/visualization/ui_validators.py:235` blocks `ligand_coupling`, `protein_coupling`, and `quenching` with planned/unimplemented messages.
- `src/emulsim/visualization/ui_model_metadata.py:153` correctly warns that `LIGAND_COUPLING`, `PROTEIN_COUPLING`, and `QUENCHING` are not implemented.

Verification performed during this audit:

- `python -m pytest -q tests\test_module2_acs.py tests\test_module2_workflows.py tests\test_ui_contract.py`
- Result: 120 passed in 0.74 s, with one pytest cache warning unrelated to model behavior.
- `python -m pytest -q tests\test_module3_breakthrough.py tests\test_module3_multicomponent.py`
- Result: 95 passed in 63.13 s, with breakthrough mass-balance warnings above 2 percent in several tests.

## 3. Critical Findings

### F1. Proposed quenching update violates the current ACS conservation law

Severity: Critical

The workplan specifies quenching as:

- `target_profile.consumed_sites += sites_consumed`
- `target_profile.blocked_sites += sites_consumed`

This conflicts with the current `ACSProfile.remaining_sites` and `ACSProfile.validate()` definitions. In the current model:

```text
remaining_sites = accessible_sites - consumed_sites - blocked_sites
validate(): consumed_sites + blocked_sites <= accessible_sites
```

If a quench consumes 95 percent of remaining activated sites and the code adds those same sites to both `consumed_sites` and `blocked_sites`, the model records 190 percent site usage. This will either violate conservation or produce a physically impossible remaining-site count.

Scientific impact:

- Quenching would make high-conversion workflows fail conservation checks.
- Downstream ligand capacity would be artificially reduced by double counting.
- M3 capacity predictions would become nonphysical if based on the post-quench ACS state.

Required correction:

- Treat quenching as `blocked_sites += sites_blocked` only, not as both consumed and blocked.
- Alternatively redefine the ACS state model so every site has exactly one terminal state: `available`, `activated`, `coupled`, `hydrolyzed`, `blocked`, `crosslinked`, or `lost`. This is preferable for long-term extensibility.
- Add a regression test where 95 percent quenching leaves approximately 5 percent remaining sites and no conservation violation.

### F2. The plan uses "remaining sites" where it should use "remaining activated sites"

Severity: Critical

Coupling and quenching should operate on activated groups such as epoxide, vinyl sulfone, aldehyde, or other activated intermediates. The plan often states `target_profile.remaining_sites`, which is currently computed from `accessible_sites - consumed_sites - blocked_sites`, not from `activated_sites - coupled - hydrolyzed - blocked`.

Scientific impact:

- Coupling could consume unactivated accessible sites if the profile is not carefully constructed.
- Quenching could block all accessible sites, not only residual activated sites.
- Product ACS profiles from activation become ambiguous because `accessible_sites` and `activated_sites` are both set equal for product sites, while parent hydroxyl sites retain separate consumed state.

Required correction:

- Add explicit ACS state fields for activated-site accounting, or define helper properties:
  - `available_for_activation`
  - `available_for_coupling`
  - `available_for_quenching`
- For activated product profiles, compute coupling/quenching capacity from `activated_sites - ligand_coupled_sites - blocked_sites - hydrolyzed_sites`, not from generic `remaining_sites`.
- Add `hydrolyzed_sites` or `deactivated_sites` because competitive hydrolysis consumes activated groups without producing functional ligand.

### F3. Competitive hydrolysis template is not sufficient as a general ligand-coupling model

Severity: High

The plan describes Template 2 as correct for epoxide and vinyl sulfone ligand coupling with aqueous hydrolysis. The current private RHS has useful structure, but it is not yet a complete wetlab-realistic model.

Current limitations:

- The template does not include pH-dependent nucleophilicity or hydrolysis.
- The template does not include reagent/ligand charge state, pKa, buffer chemistry, ionic strength, or competing nucleophiles.
- The wrapper signature includes `stoich`, `temperature`, `E_a`, and `k0`, but the current RHS itself has no stoichiometry or Arrhenius logic.
- Vinyl sulfone and epoxide chemistries should not share identical hydrolysis assumptions.
- Hydrolysis should be tracked as a terminal inactive site state, not hidden only in reagent remaining fraction.

Scientific impact:

- pH 9, pH 10.5, and pH 12 can produce very different coupling and side-reaction outcomes.
- IDA, DEAE, phenylamine/aniline-type reagents, and sulfopropyl reagents have different nucleophilicity, charge state, solubility, and safety constraints.
- A single second-order template can be useful for ranking but is not sufficient for quantitative ligand density prediction.

Required correction:

- Implement pH correction factors for coupling and hydrolysis rates.
- Record hydrolyzed/deactivated sites explicitly in ACS state.
- Require each ligand profile to declare chemistry class: epoxide-amine, epoxide-thiol, vinyl-sulfone-amine, vinyl-sulfone-thiol, aldehyde-amine, NHS-amine, etc.
- Use chemistry-specific validity windows and do not run outside them without a warning or blocker.

### F4. Proposed ligand identities are chemically under-specified

Severity: High

The plan lists DEAE, IDA, phenylamine, and sulfopropyl as ligand profiles, but practical wetlab implementation requires the actual reactive reagent, leaving group, salt form, and attachment chemistry. A "ligand" is not always the reagent that reacts with epoxide or vinyl sulfone.

Examples:

- DEAE media are commonly prepared through a specific diethylaminoethylating reagent, not just an abstract DEAE group. The actual reagent identity, pKa, salt form, and leaving group matter.
- IDA creates an IMAC chelator only after metal charging. Coupling IDA alone does not produce an IMAC medium without a metal loading step such as Ni(II), Co(II), Cu(II), or Zn(II).
- Sulfopropyl strong cation-exchange groups require a defined sulfonating/sulfopropylating reagent. The plan should not treat "sulfopropyl" as a complete reagent without specifying the reactive form.
- Phenylamine/aniline is a weak nucleophile and toxic. HIC phenyl media often use different phenyl-bearing activated reagents; hydrophobic ligand density and spacer chemistry affect selectivity.

Scientific impact:

- The planned reagent registry may give users a false sense that wetlab protocols are fully specified.
- Wrong reagent identity leads to wrong kinetics, pH, hazard controls, quench chemistry, and final ligand functionality.

Required correction:

- Split each entry into `reagent_identity`, `installed_ligand`, `functional_mode`, and `reaction_chemistry`.
- Require verified CAS and molecular form for the actual reactive reagent, not only the final ligand group.
- Add metal charging as a separate Module 2 step for IMAC.
- Add charge-state/pKa metadata for ion-exchange ligands.

### F5. Protein coupling model is useful for ranking but too simplified for affinity-media design

Severity: High

The plan's protein coupling model uses a linear steric-blocking/RSA-like term and fixed activity retention. This is a reasonable first screening approximation, but it does not capture several dominant wetlab effects.

Missing effects:

- Protein orientation and active-site accessibility.
- Multipoint attachment and local denaturation.
- pH/ionic-strength effects on protein charge and surface adsorption.
- Diffusion into pores using pore-size distributions, not only mean pore diameter.
- Coupled protein leaching or irreversible deactivation.
- Spacer-arm effects between hydrogel surface and protein ligand.
- Competition between covalent coupling and nonspecific adsorption.

Scientific impact:

- Protein A or Protein G functional capacity may differ strongly from coupled mass.
- Activity retention cannot be a fixed universal number; it depends on chemistry, pH, time, temperature, ligand concentration, support hydrophilicity, and washing.
- Epoxide-protein coupling at 4 C may preserve folding but can be kinetically slow; the reactivity/stability tradeoff must be explicit.

Required correction:

- Label protein coupling as semi-quantitative/ranking-only unless calibrated.
- Make activity retention user-adjustable with a default uncertainty range.
- Add a spacer/orientation factor or explicitly mark orientation as not modeled.
- Use ligand-accessible area and solute radius from the actual protein profile, not the default 3 nm radius.
- Add protein-specific validity windows for pH, temperature, and contact time.

## 4. High-Severity Computational and Integration Findings

### F6. Unit consistency for steric limits needs explicit conversion

Severity: High

The plan states:

```text
max_coupled = reagent.max_surface_density * surface_model.ligand_accessible_area
```

This gives units of mol/particle. The current `_steric_binding_rhs()` template describes `max_sites` as a concentration-like quantity in the ODE state units. If wrappers pass mol/particle to an ODE state expressed in mol/m3, the steric ceiling will be wrong by a bead-volume factor.

Required correction:

- Define steric limit in one canonical unit.
- If ODE state uses mol/m3, convert `max_coupled_mol_per_particle / bead_volume`.
- If ACS state uses mol/particle, solve the coupling ODE in mol/particle consistently.
- Add unit tests that compare steric saturation for different bead sizes and confirm the expected area scaling.

### F7. M2-to-M3 mapping is not sufficient for IEX, IMAC, HIC, or affinity simulation

Severity: High

The plan's integration tests say new M2 workflows should run M3 IEX, affinity, and HIC workflows. However, M3 currently uses user-supplied or default isotherm parameters. Functional ligand density from M2 is not yet converted into:

- `q_max` for ion exchange, IMAC, HIC, or affinity.
- Binding affinity constants such as `K_L`, `K_protein`, or `K_imidazole`.
- Selectivity between target and impurities.
- Salt, pH, imidazole, or ammonium sulfate dependence.
- Active protein ligand capacity after orientation/activity losses.
- Metal loading and chelator occupancy for IMAC.

Scientific impact:

- M2 may predict ligand density, but M3 can still simulate a generic column unrelated to that ligand density.
- A passing M1->M2->M3 test can be computationally valid while scientifically disconnected.

Required correction:

- Add a `FunctionalLigandProfile` or equivalent object in `FunctionalMicrosphere`.
- Map ligand states into M3 parameters with explicit confidence labels.
- For IEX, derive maximum ionic capacity from functional charge density and pH/pKa.
- For IMAC, add metal charging efficiency and imidazole competition.
- For protein affinity, derive active binding capacity from functional protein ligand density and ligand stoichiometry.
- For HIC, require salt-dependent retention parameters; hydrophobic ligand density alone is not enough.

### F8. Validator expansion must be backend-level, not UI-only

Severity: High

The workplan adds ordering rules to `ui_validators.py`. This is necessary but not sufficient. Scientific validity should not depend on Streamlit UI validation. CLI, tests, notebooks, and future API calls can bypass UI validators.

Required correction:

- Add workflow validation in `ModificationOrchestrator.run()` or a Module 2 workflow validator.
- Enforce coupling-after-activation and no-step-after-quench before solving.
- Validate reagent target compatibility against the current ACS state.
- Validate that `step.target_acs` matches `reagent_profile.target_acs`.
- Validate that `step.product_acs` is either absent or matches the reagent profile.

Existing validator issue:

- `validate_m2_inputs()` currently checks `getattr(acs_state, "accessible", None)`, but `ACSProfile` uses `accessible_sites`. This means the accessible-site blocker may not fire for real `ACSProfile` objects unless this is corrected.

### F9. Public ODE wrappers need solver quality outputs, not only conversion

Severity: Medium-high

The planned wrappers return only conversion and remaining fractions. That is not enough for a scientific simulation system.

Required correction:

- Return a result dataclass containing conversion, terminal concentrations, hydrolysis fraction, solver success, solver message, time grid, mass/site balance residual, and warnings.
- Store solver diagnostics in `ModificationResult.notes` or structured fields.
- Gate trust if the ODE solver fails, if site balance error exceeds tolerance, or if conversion is clipped.

## 5. Wetlab and Physicochemical Validity Assessment

The plan's overall chemistry direction is realistic: hydrogel supports are often activated, coupled to ligands, quenched, washed, and packed into columns. The current simplification is acceptable for an early simulator if it is framed as screening.

Wetlab issues that must be added before production use:

- Hazard and handling metadata for ECH, DVS, glutaraldehyde, aniline/phenylamine, 2-mercaptoethanol, sodium borohydride, acetic anhydride, and alkaline activation conditions.
- Buffer compatibility and competing nucleophile effects.
- pH-dependent ionization of amine, thiol, carboxyl, IDA, DEAE, SP, protein residues, and target proteins.
- Hydrogel swelling or degradation at high pH, high salt, organic solvent exposure, and long reaction times.
- Washing efficiency, residual toxic reagent, and ligand leaching.
- Batch-to-batch DDA variability and uncertainty propagation from Module 1 ACS estimates.
- Pore-size distribution and protein exclusion, not only mean pore size.

The plan's disclaimer correctly states that rate constants are order-of-magnitude estimates. The UI and report outputs should go further: default rates should be labeled "illustrative defaults requiring calibration" wherever they are displayed, and absolute ligand densities should not be presented as validated wetlab yields.

## 6. Review of Proposed Implementation Phases

### Phase A: ReagentProfile expansion

Verdict: Accept with revision.

Required changes:

- Add fields for pH validity range, temperature validity range, solvent/buffer assumptions, reaction chemistry, hazard class, calibration source, and confidence tier.
- Do not rely on a single overloaded `ReagentProfile` for crosslinkers, activators, small ligands, proteins, and quenchers unless type-specific fields are validated.
- Verify CAS numbers and actual reactive reagent identity before treating profiles as chemistry references.

### Phase B: ODE wrappers

Verdict: Accept with revision.

Required changes:

- Include Arrhenius handling and pH correction in wrappers.
- Return structured solver diagnostics.
- Track hydrolyzed/deactivated ACS explicitly.
- Confirm stoichiometry is included in ligand depletion and site consumption.

### Phase C: Ligand coupling solver

Verdict: Scientifically plausible but incomplete.

Required changes:

- Couple only activated product profiles.
- Add `hydrolyzed_sites` or equivalent terminal inactive state.
- Use reagent-specific pH and hydrolysis models.
- Use actual reactive reagent concentration, not final ligand label alone.

### Phase D: Protein coupling solver

Verdict: Useful for ranking only.

Required changes:

- Enforce steric limit with correct units.
- Use ligand-accessible area with protein-specific radius.
- Make activity retention user-adjustable and uncertainty-bearing.
- Add orientation/spacer limitations or explicit "not modeled" warnings.

### Phase E: Quenching solver

Verdict: Must be corrected before implementation.

Required changes:

- Do not increment both `consumed_sites` and `blocked_sites` under the current ACS model.
- Quench only residual activated sites.
- Add quench-specific terminal states and conservation tests.
- Include side effects such as NaBH4 reduction of aldehyde/imine states only if those states exist in ACS bookkeeping.

### Phase F: Dispatch integration

Verdict: Accept only after backend validators exist.

Required changes:

- Do not dispatch based only on step enum.
- Validate reagent compatibility and state availability before solving.

### Phase G: M2 ordering validation

Verdict: Necessary but should not be UI-only.

Required changes:

- Put canonical ordering validation in Module 2 backend.
- Keep UI validators as user-facing preflight checks.

### Phase H: UI integration

Verdict: Premature as written.

Required changes:

- Keep new workflows disabled until backend tests pass.
- Show planned workflows as disabled with scientific rationale.
- Add confidence labels, rate-constant calibration warnings, and hazard disclaimers.

### Phase I: Integration testing

Verdict: Essential but incomplete.

Required changes:

- Test M2-to-M3 parameter mapping, not only that M3 runs.
- Include capacity monotonicity tests: higher functional ligand density should not reduce `q_max` unless pressure/pore-blockage penalties are explicitly modeled.
- Include invalid chemistry tests: IDA without metal charging should not behave as IMAC.
- Include quench conservation tests at >95 percent conversion.

## 7. Revised Implementation Strategy

### Step 1: Fix ACS state semantics first

Add explicit terminal states before adding new workflows:

- `hydrolyzed_sites` or `deactivated_sites`
- `crosslinked_sites` if crosslinking should be separated from generic consumed sites
- `blocked_sites`
- `ligand_coupled_sites`
- `ligand_functional_sites`

Define exactly one conservation equation and keep it invariant across all workflows.

### Step 2: Add a typed chemistry registry

Replace a flat reagent list with type-aware profiles:

- `CrosslinkerProfile`
- `ActivatorProfile`
- `SmallMoleculeLigandProfile`
- `ProteinLigandProfile`
- `QuencherProfile`
- `MetalChargingProfile` for IMAC

Each profile should include chemistry class, target ACS, product/terminal state, pH range, temperature range, default concentration, hydrodynamic radius, activity/charge parameters, hazard notes, calibration source, and model confidence.

### Step 3: Implement public ODE wrappers with diagnostics

Each solver wrapper should return structured diagnostics and balance errors. This is required for trust labels and for debugging wetlab-plausibility failures.

### Step 4: Implement backend workflow validation

The backend should block:

- Coupling without activated sites.
- Quenching without activated sites.
- Steps after terminal quenching unless explicitly allowed.
- Reagent-target mismatches.
- Protein coupling when ligand-accessible area is effectively external-only or pore-excluded beyond a threshold.

### Step 5: Implement one validated path before all paths

Recommended first production path:

```text
ECH activation -> ethanolamine quench
```

Recommended first functional ligand path:

```text
ECH activation -> one small-molecule ion-exchange ligand -> quench
```

Defer protein A/G coupling until steric, orientation, and activity-retention assumptions are explicitly labeled and tested.

### Step 6: Connect M2 to M3 through a functional media contract

Add a stable contract between Module 2 and Module 3:

- ligand type
- functional ligand density
- active ligand density
- accessible ligand area
- pore accessibility for target molecule
- estimated `q_max`
- confidence and calibration notes
- unsupported performance mappings

M3 should consume this contract rather than relying on unrelated default isotherm parameters.

## 8. Acceptance Criteria Needed Before Merge

Minimum scientific acceptance criteria:

- ACS conservation holds for every workflow, including >95 percent quenching.
- Hydrolyzed/deactivated activated sites are accounted for.
- pH and temperature are either mechanistically modeled or clearly labeled as validation metadata.
- Every new reagent profile has verified identity, target ACS, chemistry class, pH range, and confidence label.
- Protein coupling is labeled ranking-only unless calibrated.
- M2 functional ligand density maps to M3 `q_max` or is explicitly labeled "not mapped to process performance."
- UI does not allow unsupported workflows to run.

Minimum computational acceptance criteria:

- Unit tests for all 10 planned profiles.
- Unit tests for public ODE wrappers, including zero concentration, high concentration, hydrolysis-dominated, and steric-limited cases.
- Conservation tests after each step and full workflow.
- Backend workflow ordering tests independent of Streamlit.
- M1->M2->M3 integration tests that verify parameter propagation, not just successful execution.
- Solver diagnostics are captured and surfaced.

## 9. Final Scientific Judgment

The workplan is a valuable next-phase proposal, but it should be revised before implementation. The direction is correct, and the existing codebase has a strong foundation: `M1ExportContract`, `AccessibleSurfaceModel`, `ACSProfile`, Arrhenius-enabled current M2 workflows, UI validators, and Module 3 process models.

The plan's current weak points are ACS state semantics, wetlab chemistry specificity, pH dependence, steric/protein realism, and M2-to-M3 mapping. The quenching double-counting issue is a blocker because it directly contradicts the current conservation law. The M2-to-M3 mapping gap is the second major blocker because functional media design is only meaningful if functional ligand density changes the process model.

Recommended go/no-go decision:

- Go for a revised, phased implementation focused first on corrected ACS states and one small-molecule workflow.
- No-go for implementing all three new step types exactly as written.
- No-go for UI exposure of executable ligand/protein/quench workflows until conservation, validation, and M3 mapping tests pass.

