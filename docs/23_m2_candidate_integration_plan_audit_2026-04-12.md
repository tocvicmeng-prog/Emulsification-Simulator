# Audit Report: M2 Candidate Integration Plan

Date: 2026-04-12  
Auditor role: computational simulation scientist, wet-lab chemistry researcher, physicochemical model reviewer  
Scope: `docs/19_ligand_protein_coupling_candidates.md`, `docs/20_linker_arm_candidates.md`, `docs/21_m2_candidates_integration_plan.md`, `docs/22_m2_detailed_integration_plan.md`, and the current Module 2 implementation.

## 1. Executive Verdict

The M2 candidate integration plan is directionally sound and scientifically useful as a semi-quantitative expansion of the current functionalization simulator. It correctly recognizes that M2 must bridge M1-generated base microspheres to M3 process simulations by producing chemically meaningful ligand density, residual ACS state, accessibility, mechanical updates, and a `FunctionalMediaContract`.

However, the plan should not be described as "production-grade" without additional calibration, chemistry-specific states, and explicit limits. The proposed candidate set covers the major chromatography modes, but several planned profiles require clearer reagent identity, coupling orientation, charge state, metal-loading state, and M2-to-M3 capacity mapping. The current implementation is reliable for mass/accounting consistency and workflow execution, but the planned expansion remains semi-quantitative until calibrated against real resin specifications and wet-lab measurements.

Primary conclusion: accept the plan for Phase 1 only if it is explicitly labeled as an expanded semi-quantitative candidate library, with mandatory safeguards for charge typing, spacer interpretation, q_max area basis, IMAC metal loading, and UI confidence/hazard disclosure.

## 2. Sources Reviewed

Planning and study documents:

- `docs/19_ligand_protein_coupling_candidates.md`
- `docs/20_linker_arm_candidates.md`
- `docs/21_m2_candidates_integration_plan.md`
- `docs/22_m2_detailed_integration_plan.md`

Implementation areas reviewed:

- `src/emulsim/module2_functionalization/acs.py`
- `src/emulsim/module2_functionalization/reactions.py`
- `src/emulsim/module2_functionalization/modification_steps.py`
- `src/emulsim/module2_functionalization/reagent_profiles.py`
- `src/emulsim/module2_functionalization/orchestrator.py`
- `src/emulsim/visualization/app.py`
- `src/emulsim/visualization/ui_validators.py`

Tests run:

```powershell
python -m pytest -q tests\test_module2_acs.py tests\test_module2_workflows.py tests\test_ui_contract.py
```

Result: 130 passed, 1 non-functional pytest cache warning.

Previously observed M3 verification:

```powershell
python -m pytest -q tests\test_module3_breakthrough.py tests\test_module3_multicomponent.py
```

Result: 95 passed, with known mass-balance warning cases above the 2 percent threshold in some breakthrough scenarios.

## 3. Current System Structure and M2 Role

### 3.1 Current three-module system intent

The project is evolving into an integrated simulation chain:

1. M1: base microsphere fabrication by double emulsification.
2. M2: post-fabrication chemical modification, activation, crosslinking, ligand/protein immobilization, and quenching.
3. M3: process performance simulation in purification columns or catalytic beds.

M2 is therefore not simply a reagent selector. It is the scientific bridge that converts a fabricated microsphere into a functional chromatographic or catalytic medium.

### 3.2 Desired M2 function

M2 should accept the M1-derived base microsphere state and output:

- Updated ACS inventory, including residual, activated, consumed, hydrolyzed, coupled, functional, crosslinked, and blocked sites.
- Modified mechanical and transport properties after secondary crosslinking and functionalization.
- Ligand/protein density on a physically meaningful area basis.
- Accessibility-limited functional density, not just nominal coupling yield.
- A `FunctionalMediaContract` suitable for M3 adsorption, affinity, IMAC, HIC, IEX, or catalytic simulations.
- Confidence, calibration, hazard, and validity metadata so the UI and downstream modules do not imply false precision.

### 3.3 Current M2 implementation state

The current code is substantially more mature than the earlier two-step M2 concept:

- `acs.py` now implements an explicit terminal-state ACS model with `crosslinked_sites`, `activated_consumed_sites`, `hydrolyzed_sites`, `ligand_coupled_sites`, `ligand_functional_sites`, and `blocked_sites`.
- `modification_steps.py` dispatches five workflow types: secondary crosslinking, activation, ligand coupling, protein coupling, and quenching.
- `reactions.py` includes competitive coupling, steric coupling, and quenching solver wrappers with hydrolysis and site-balance reporting.
- `reagent_profiles.py` currently includes 14 profiles: 2 secondary crosslinkers, 2 activators, 4 ligand coupling profiles, 2 protein coupling profiles, and 4 quench profiles.
- `orchestrator.py` includes `FunctionalMicrosphere`, `FunctionalMediaContract`, workflow ordering validation, and initial q_max mapping.
- `visualization/app.py` exposes M2 secondary crosslinking, hydroxyl activation, ligand coupling, protein coupling, and quenching in the UI.

The current implementation is a valid foundation for the planned expansion. The main remaining issue is not software architecture; it is scientific parameterization, interpretation, and contract fidelity.

## 4. Assessment of the Candidate Expansion Plan

### 4.1 Overall candidate set

The proposed high-priority additions are appropriate for a chromatography-oriented platform:

- Q and CM fill missing strong anion and weak cation exchange coverage.
- NTA complements IDA for IMAC and is important for His-tag purification.
- Butyl complements Phenyl for HIC selectivity.
- Glutathione enables GST-tag affinity capture.
- Heparin enables mixed affinity/cation-exchange capture of heparin-binding proteins.
- Streptavidin enables biotin-tag affinity capture.
- Protein A/G fusion is a reasonable extension of the existing Protein A and Protein G profiles.

This set is scientifically relevant and commercially recognizable. It also maps well to the intended M3 process use cases.

### 4.2 Main concern: final ligand is not enough

For M2, each profile must define not only the final functional ligand but also the chemically plausible immobilization path. The plan sometimes treats final ligand identity as sufficient. That is risky because wet-lab immobilization outcomes depend on:

- Reactive precursor identity.
- Which functional group reacts with the activated microsphere.
- Orientation and multipoint attachment.
- Competing hydrolysis or side reactions.
- Whether the immobilized form preserves target-binding function.
- Whether the installed ligand introduces charge, hydrophobicity, swelling, or nonspecific binding.

The current plan is acceptable for ranking and exploratory simulation, but not for claiming predictive quantitative performance.

## 5. Scientific Validity Findings

### F1. The "production-grade library" claim is overstated

Severity: High  
Area: scientific validity, user trust  
Evidence: `docs/22_m2_detailed_integration_plan.md` proposes a "production-grade library" while also relying on estimated kinetics, activity multipliers, and simplified q_max mappings.

Assessment:

The plan is scientifically useful, but the data quality does not support a production-grade claim. Most kinetic constants, activity-retention values, steric limits, and q_max estimates are not currently calibrated to specific wet-lab protocols or commercial resin specifications. Profiles such as Protein A/G, heparin, streptavidin, and glutathione are highly sensitive to immobilization orientation and target molecule identity.

Recommendation:

Rename the planned library to "expanded semi-quantitative candidate library." Use `confidence_tier="semi_quantitative"` or `confidence_tier="ranking_only"` prominently in UI and reports. Reserve "production-grade" for calibrated profiles with documented source data and validation against measured resin capacity.

### F2. Spacer-as-multiplier is useful but not mechanistic

Severity: High  
Area: spacer chemistry, protein immobilization, process validity  
Evidence: `docs/20_linker_arm_candidates.md` proposes Phase 1 spacer fields and a multiplier on `activity_retention`; `docs/21_m2_candidates_integration_plan.md` and `docs/22_m2_detailed_integration_plan.md` adopt this approach.

Assessment:

The multiplier model is acceptable for Phase 1 ranking, especially for protein ligands. It captures the qualitative fact that spacer arms can reduce steric occlusion and improve retained activity. It does not simulate:

- Spacer coupling yield.
- Creation of a distal ACS state.
- Distal group chemistry, such as amine, carboxyl, thiol, or epoxide.
- Additional activation steps such as EDC/NHS after AHA.
- Multipoint attachment or crosslinking caused by diamines such as DADPA and DAH.
- Charge and nonspecific binding introduced by spacer amines.
- Hydrolytic instability or spacer cleavage.

Recommendation:

Keep the Phase 1 multiplier, but label it as an empirical accessibility/activity correction. Add `spacer_key`, `spacer_length_angstrom`, and `spacer_activity_multiplier` only as metadata-backed correction fields. Do not present spacer profiles as executable modification steps in the UI until a true `SPACER_ARM` step is implemented.

### F3. Q ligand chemistry needs explicit reagent identity

Severity: High  
Area: ion-exchange chemistry  

Assessment:

Q is a strong anion exchanger because the immobilized group carries a permanent positive quaternary ammonium charge. The plan must specify how the quaternary group is installed. A final "Q ligand" label is not sufficient because a pre-quaternized ammonium group is not generally modeled as the same nucleophilic epoxide-coupling species as DEAE or other amines.

Practical wet-lab routes may involve glycidyltrimethylammonium-type reagents, tertiary amine installation followed by quaternization, or commercial activated matrices with predesigned chemistries. These routes have different kinetics, side reactions, charge density, and hydrolysis behavior.

Recommendation:

Add explicit fields for `reagent_identity`, `reactive_group`, and `installed_ligand`. If the model uses a simplified Q profile, mark it as a functional approximation and not a mechanistic epoxide-amine coupling profile unless the actual reactive precursor supports that mechanism.

### F4. CM ligand chemistry is functionally plausible but chemically ambiguous

Severity: Medium-High  
Area: cation-exchange chemistry  

Assessment:

CM resins are commonly carboxymethylated polysaccharide matrices, often prepared through chloroacetic-acid chemistry on hydroxyl groups. If the plan models CM as epoxide coupling of an amino-acid-like carboxylated nucleophile, the installed group may function as a weak cation exchanger, but it is not the same chemistry as classical carboxymethyl agarose/cellulose formation.

Recommendation:

Distinguish "CM-like cation exchanger installed by amino-carboxyl ligand coupling" from "true carboxymethylation." If the intended profile is glycine or aminoacetic-acid coupling, name the precursor accordingly and document that the output is a carboxyl-bearing weak cation exchanger.

### F5. NTA and IDA require metal-loading state before M3 IMAC prediction

Severity: High  
Area: IMAC validity, M2->M3 contract  
Evidence: `docs/22_m2_detailed_integration_plan.md` explicitly treats NTA metal charging as implicit in Phase 1.

Assessment:

Uncharged IDA or NTA is not an IMAC resin. IMAC function depends on the chelated metal ion, metal loading fraction, ligand denticity, pH, buffer competition, imidazole, reducing agents, and metal leaching. Treating the NTA profile as implicitly charged can be acceptable for a simplified Phase 1 UI, but it is not sufficient for M3 process predictions where His-tag binding and imidazole elution are expected.

Recommendation:

Add at least a minimal `metal_ion`, `metal_loaded_fraction`, and `metal_loading_confidence` state to `FunctionalMediaContract` for IMAC profiles. If not implemented in Phase 1, M3 should mark IMAC q_max as estimated only under an explicit assumption such as "fully Ni2+-loaded chelator, no metal leaching, no competing chelators."

### F6. Heparin should be treated as a macromolecular ligand, not a small molecule

Severity: High  
Area: accessibility, steric exclusion, mixed-mode binding  
Evidence: `docs/19_ligand_protein_coupling_candidates.md` describes heparin as a polysaccharide but also says it is immobilized like a small molecule.

Assessment:

Heparin is a heterogeneous polyanionic polysaccharide with variable molecular weight, charge density, sulfation pattern, and multipoint attachment behavior. It can function as both affinity ligand and weak cation exchanger. Its installation and performance are not equivalent to a small ligand such as DEAE, IDA, or phenyl.

Recommendation:

Set `is_macromolecule=True` or otherwise route heparin through macromolecule accessibility constraints. Include `ligand_mw_distribution` or at least broad uncertainty metadata. Use heparin-specific q_max confidence as low or approximate because stoichiometry depends strongly on the target protein.

### F7. Streptavidin-biotin binding requires special process interpretation

Severity: High  
Area: affinity simulation, elution realism  
Evidence: `docs/19_ligand_protein_coupling_candidates.md` highlights streptavidin-biotin Kd around 10^-15 M; `docs/21_m2_candidates_integration_plan.md` proposes stoichiometry of 4.

Assessment:

Streptavidin-biotin binding is often effectively irreversible under mild chromatography conditions. A simple q_max estimate with stoichiometry 4 can overstate usable dynamic capacity. Steric constraints, orientation, tetramer accessibility, partial inactivation, and harsh elution requirements reduce practical performance.

Recommendation:

Mark streptavidin as `ranking_only` unless calibrated. In M3, represent the default behavior as near-irreversible capture or require the user to select a reversible Strep-tag-like binding system. Cap effective stoichiometry below 4 unless calibration data supports full accessibility.

### F8. Glutathione coupling requires orientation control

Severity: Medium-High  
Area: GST affinity chemistry  

Assessment:

Glutathione has multiple reactive groups and a specific recognition surface for GST. Random immobilization through the wrong group can reduce binding. Treating glutathione as a generic small ligand with full pore accessibility may overpredict functional density.

Recommendation:

Add orientation uncertainty to the profile. The profile should identify the coupling group and apply an activity-retention penalty unless a spacer and orientation-preserving chemistry are specified.

### F9. Butyl and Phenyl HIC cannot be converted to q_max by ligand density alone

Severity: Medium-High  
Area: HIC process validity  
Evidence: `docs/22_m2_detailed_integration_plan.md` correctly states HIC q_max should be `not_mapped`.

Assessment:

The plan is correct that HIC capacity depends heavily on salt concentration, protein hydrophobic patches, ligand type, temperature, pH, and additives. HIC density is useful as an M2 material descriptor, but M3 cannot infer a reliable q_max from ligand density alone.

Recommendation:

Keep HIC `q_max_confidence="not_mapped"` unless M3 has a salt-dependent adsorption isotherm for the selected target protein. The UI should not report HIC capacity as a numeric q_max by default.

### F10. Current pH handling is a validity gate, not a kinetic model

Severity: Medium  
Area: reaction kinetics  
Evidence: `reactions.py` has pH range warning checks, but pH does not yet mechanistically modify rate constants or ligand ionization.

Assessment:

The current implementation warns when pH is outside a valid range. That is useful but not equivalent to pH-dependent chemistry. Real coupling rates and final functionality depend on nucleophile protonation, epoxide/vinyl sulfone reactivity, hydrolysis, ligand charge state, and protein stability.

Recommendation:

Continue using pH validity warnings in Phase 1, but avoid claiming pH optimization. For future releases, add pH-dependent rate scaling for amines, thiols, carboxyl activation, and protein stability.

## 6. Computational and Model Reliability Findings

### F11. ACS accounting is now structurally reliable

Severity: Positive finding  
Area: mass/site conservation  

Assessment:

The current terminal-state ACS model is a major improvement. It prevents the earlier ambiguity between consumed, coupled, hydrolyzed, and blocked sites. This makes ligand density and residual ACS reporting more defensible.

Remaining limitation:

ACS density is still an abstract site inventory. It does not by itself guarantee physical pore accessibility, reaction penetration, or experimentally measurable ligand density. Those require calibration.

### F12. Competitive coupling and steric coupling solvers are appropriate Phase 1 templates

Severity: Positive finding with limitations  
Area: numerical modeling  

Assessment:

The competitive coupling solver with hydrolysis and the steric coupling solver for macromolecular ligands are appropriate templates. They provide better scientific structure than fixed-conversion formulas.

Remaining limitations:

- Diffusion-limited penetration is simplified.
- Reagent depletion is bulk-like and not spatially resolved.
- Protein immobilization orientation is not mechanistic.
- Multipoint attachment is not explicitly modeled.
- pH and ionic strength are not mechanistic modifiers.

The solvers are acceptable for semi-quantitative design screening, not predictive wet-lab optimization.

### F13. Current FunctionalMediaContract q_max mapping may mismatch area basis

Severity: High  
Area: M2->M3 consistency  
Evidence: `orchestrator.py` estimates q_max from functional density and packed-bed external particle surface area.

Assessment:

The plan uses:

```text
q_max [mol/m3 bed] = functional_ligand_density [mol/m2] * specific_surface_area [m2/m3 bed] * stoichiometry
specific_surface_area = 6 * (1 - eps_bed) / d_particle
```

This is dimensionally correct only if the ligand density is defined on the same surface area basis as the specific surface area. In porous microspheres, M2 ligand density may be based on internal and external accessible surface area, while `6*(1-eps)/dp` is an external geometric packed-particle area. Combining these without an area-basis correction can undercount or overcount capacity.

Recommendation:

Add explicit fields to `FunctionalMediaContract`:

- `ligand_density_area_basis`: external, internal_plus_external, reagent_accessible, ligand_accessible, or calibrated_bed.
- `accessible_area_per_bed_volume`.
- `q_max_area_basis_note`.

Use the M2 accessible area per particle and particle number per bed volume when possible, rather than only external particle area. If that is not available, label q_max as a rough estimate.

### F14. Explicit `charge_type` is mandatory for IEX

Severity: High  
Area: charge assignment, UI correctness, M3 mapping  
Evidence: `docs/21_m2_candidates_integration_plan.md` and `docs/22_m2_detailed_integration_plan.md` correctly propose `charge_type`; current code still uses string heuristics.

Assessment:

The current heuristic can classify DEAE but will not robustly classify Q, CM, future mixed-mode ligands, or custom names. IEX direction is a fundamental physical property: anion exchangers carry positive charge and bind anions; cation exchangers carry negative charge and bind cations.

Recommendation:

Implement `charge_type` before adding Q and CM. Do not expand IEX profiles while relying on string parsing.

### F15. New functional modes require stricter M3 contract behavior

Severity: Medium-High  
Area: M3 compatibility  
Evidence: `docs/21_m2_candidates_integration_plan.md` proposes `gst_affinity`, `biotin_affinity`, and `heparin_affinity`; current code maps only generic IEX, affinity, IMAC, and HIC.

Assessment:

New functional modes are scientifically justified, but they are not interchangeable:

- GST/glutathione is reversible and eluted with free glutathione.
- Streptavidin/biotin is effectively irreversible under mild conditions.
- Heparin is mixed-mode and target-dependent.
- Protein A/G has Fc-specific behavior and pH-sensitive elution.
- IMAC depends on metal and imidazole.

Recommendation:

`FunctionalMediaContract` should include both `ligand_type` and `binding_model_hint`, for example:

- `iex_charge_exchange`
- `imac_metal_chelation`
- `hic_salt_promoted`
- `protein_a_fc_affinity`
- `gst_glutathione_affinity`
- `biotin_streptavidin_near_irreversible`
- `heparin_mixed_mode`

This prevents M3 from treating all affinity ligands as one generic model.

## 7. UI and User-Workflow Findings

### F16. UI must expose confidence, assumptions, and hazards

Severity: Medium-High  
Area: user reliability  

Assessment:

The current UI exposes the main M2 workflow types and current reagent profiles. The expansion plan adds more profiles and spacer metadata. That will increase user-facing complexity and risk of false confidence.

Recommendation:

For every new profile, the UI should display:

- `confidence_tier`
- `calibration_source`
- `hazard_class`
- valid pH and temperature range
- functional mode
- whether q_max is mapped, approximate, or not mapped
- spacer assumption, if any

Do not hide these as developer-only metadata. They are essential scientific context.

### F17. Spacer profiles should not be selectable as independent reactions in Phase 1

Severity: Medium  
Area: UI correctness  
Evidence: `docs/21_m2_candidates_integration_plan.md` says spacer profiles are metadata and not independent steps.

Assessment:

If spacer profiles are placed in `REAGENT_PROFILES`, the UI must not automatically show them under ligand/protein coupling or quenching. A spacer selected as a normal reagent would imply an executable chemistry step that Phase 1 does not simulate.

Recommendation:

Filter `functional_mode="spacer"` from executable reagent dropdowns. Present spacer selection only as an optional field associated with compatible ligand/protein profiles.

### F18. Profile count bookkeeping is inconsistent across documents

Severity: Low-Medium  
Area: planning reliability  

Assessment:

The documents refer to counts such as "6-ligand + 2-protein", "14 existing + 8 coupling + 3 spacer", "25 profiles", and "22 coupling profiles + 3 spacer profiles." The current code contains 14 total profiles, not all of which are coupling profiles. The target count should be made unambiguous.

Recommendation:

Use one canonical target table:

- Existing total profiles: 14.
- Existing executable coupling/protein profiles: 6.
- New executable coupling/protein profiles: 8.
- New spacer metadata profiles: 3.
- Target total `REAGENT_PROFILES`: 25.

If "22 coupling profiles" means a different classification, define it explicitly.

## 8. Candidate-by-Candidate Audit

| Candidate | Scientific value | Main concern | Required plan change |
|---|---:|---|---|
| Q | High | Reactive precursor and quaternization route not explicit | Add actual reagent identity and `charge_type="anion"` |
| CM | High | Classical CM chemistry differs from amino-carboxyl ligand coupling | Rename/define precursor; add `charge_type="cation"` |
| NTA | High | IMAC requires metal charging and metal state | Add `metal_ion`, `metal_loaded_fraction`, and IMAC assumption notes |
| Butyl | Medium-High | HIC q_max not ligand-density-only; installed amine may add charge | Keep q_max `not_mapped`; document ligand chemistry |
| Glutathione | High | Orientation and multipoint reactions affect GST binding | Add coupling-site and activity-retention uncertainty |
| Heparin | High | Macromolecular, heterogeneous, mixed-mode | Treat as macromolecule; lower q_max confidence |
| Streptavidin | High | Near-irreversible binding; stoich 4 is theoretical maximum | Add near-irreversible M3 model hint and lower effective stoich |
| Protein A/G fusion | Medium | Commercial value lower than separate A/G unless calibrated | Keep ranking-only until validated |
| DADPA spacer | High | Can add charge and multipoint attachment | Use multiplier only; document amine charge/nonspecific binding |
| AHA spacer | High | COOH distal chemistry requires EDC/NHS step not modeled | Do not imply explicit AHA chemistry in Phase 1 |
| DAH spacer | Medium-High | Diamine side reactions and hydrophobicity | Keep as lower-confidence spacer metadata |

## 9. Recommended Corrections Before Implementation

### 9.1 Mandatory for Phase 1

1. Add `charge_type` to `ReagentProfile` and replace string-based IEX classification.
2. Add `spacer_key`, `spacer_length_angstrom`, and `spacer_activity_multiplier`, but keep spacer profiles non-executable.
3. Cap `activity_retention * spacer_activity_multiplier` at 1.0 and record the effective value.
4. Add explicit `functional_mode` mappings for `gst_affinity`, `biotin_affinity`, and `heparin_affinity`.
5. Keep HIC q_max as `not_mapped`.
6. Mark IMAC q_max as assumption-based unless metal loading state is represented.
7. Add q_max notes explaining the surface-area basis.
8. Update UI dropdowns to include new executable profiles and exclude spacer metadata profiles.
9. Display confidence tier, calibration source, hazard class, and q_max confidence in the UI.
10. Add tests for all new profile fields, spacer multiplier behavior, charge classification, FunctionalMediaContract mapping, and UI filtering.

### 9.2 Strongly recommended for Phase 1.5

1. Add `metal_ion` and `metal_loaded_fraction` for IMAC.
2. Add `binding_model_hint` to `FunctionalMediaContract`.
3. Add `ligand_density_area_basis` and `accessible_area_per_bed_volume`.
4. Add orientation/activity uncertainty metadata for protein, streptavidin, heparin, and glutathione profiles.
5. Add target-specific default assumptions for IgG, His-tag, GST-tag, biotin-tag, and heparin-binding proteins.

### 9.3 Recommended for Phase 2

1. Implement an explicit `SPACER_ARM` modification step.
2. Model spacer distal ACS states.
3. Add EDC/NHS chemistry for carboxyl-terminated spacers.
4. Add pH-dependent reaction kinetics and ligand ionization.
5. Add metal charging/regeneration and imidazole competition for IMAC.
6. Add salt-dependent HIC isotherms.
7. Add reversible/irreversible affinity binding modes for M3.

## 10. Revised Acceptance Criteria

The current acceptance criteria should be tightened as follows:

| Area | Required acceptance criterion |
|---|---|
| Profile completeness | 25 total profiles present: 14 existing, 8 new executable coupling/protein profiles, 3 spacer metadata profiles |
| Regression safety | Existing 14 profiles retain backward-compatible behavior unless deliberately changed |
| Charge typing | Every IEX profile has explicit `charge_type`; no IEX mapping depends on string matching |
| Spacer behavior | Spacer profiles cannot be executed as standalone steps in Phase 1 |
| Spacer multiplier | Effective activity retention is capped at 1.0 and tested |
| IMAC | NTA/IDA contracts state metal-loading assumptions or carry metal-loading fields |
| HIC | HIC contracts return q_max `not_mapped` unless a salt-dependent target model is available |
| Affinity modes | GST, biotin, heparin, Protein A/G, and Protein A/G fusion map to distinct model hints |
| q_max units | q_max mapping records the area basis and confidence |
| UI | UI displays confidence tier, calibration source, pH/temperature validity, hazards, and q_max confidence |
| Tests | Unit and integration tests cover profile completeness, charge mapping, spacer filtering, q_max branches, and full M1->M2->M3 contract generation |
| Documentation | Documentation states that Phase 1 is semi-quantitative and requires wet-lab calibration |

## 11. Calibration Requirements

Before using the expanded M2 library for quantitative predictions, at least the following wet-lab calibration data should be collected:

- Base microsphere size distribution, porosity, pore size distribution, swelling ratio, and accessible surface area.
- ACS density by chemical titration or spectroscopic derivatization.
- Activation density after ECH/DVS or other activation chemistry.
- Ligand density by elemental analysis, dye assay, UV/Vis assay, ICP-MS for metal chelators, or ligand-specific quantification.
- Residual reactive group density after quenching.
- Static binding capacity and dynamic binding capacity for representative targets.
- Pressure-flow behavior after functionalization.
- Breakthrough curves at multiple flow rates.
- Salt/pH/elution dependence for IEX, HIC, IMAC, and affinity modes.
- Protein activity retention after immobilization for Protein A/G, streptavidin, and enzyme/catalytic ligands.

Recommended calibration targets:

- Protein A/G: IgG static and dynamic binding capacity, mg/mL resin.
- IMAC: His-tag protein capacity and imidazole elution behavior, with metal leaching.
- IEX: BSA, lysozyme, or model protein binding across pH and conductivity.
- HIC: target protein binding across ammonium sulfate gradients.
- GST: GST-tag protein binding and glutathione elution.
- Streptavidin: biotinylated analyte capacity and elution feasibility.
- Heparin: representative heparin-binding protein capacity and salt elution.

## 12. Final Reliability Assessment

Current implementation reliability:

- Software workflow reliability: good for the implemented 14-profile M2 system.
- ACS/site conservation reliability: good, based on terminal-state accounting and passing tests.
- Numerical reaction-template reliability: acceptable for semi-quantitative screening.
- Scientific predictiveness: moderate for relative trends, low-to-moderate for absolute capacity without calibration.
- M2->M3 contract reliability: acceptable as a structural bridge, but q_max mapping requires area-basis and mode-specific corrections.

Plan reliability after proposed expansion:

- Reliable for exploratory design ranking if confidence tiers are shown.
- Not reliable as a quantitative wet-lab substitute without calibration.
- Not reliable for IMAC, HIC, streptavidin, heparin, or protein immobilization if the plan continues to treat these as simple ligand-density transformations.
- Suitable for the next development phase if the implementation includes explicit charge typing, non-executable spacer metadata, mode-specific M3 hints, and conservative q_max confidence labels.

## 13. Final Recommendation

Proceed with the candidate integration plan, but revise its framing and acceptance criteria before implementation. The plan should be implemented as an expanded semi-quantitative M2 library that improves coverage of chromatography chemistries and prepares the M2->M3 bridge for more realistic process modeling.

The highest-priority corrections are:

1. Replace IEX string heuristics with explicit `charge_type`.
2. Keep spacer profiles as metadata-only in Phase 1.
3. Add q_max area-basis metadata to `FunctionalMediaContract`.
4. Add IMAC metal-loading assumptions or fields.
5. Treat heparin and streptavidin as special high-uncertainty affinity cases.
6. Require UI disclosure of confidence tier, calibration status, hazards, and q_max confidence.
7. Rename the plan from "production-grade library" to "expanded semi-quantitative candidate library" until calibration is complete.

With these changes, the plan will be scientifically consistent enough for continued development and useful simulation-guided prioritization. Without these changes, the UI and M3 outputs may overstate predictive accuracy and mislead users about the real wet-lab performance of functionalized microspheres.
