# Audit Report: Module 2 v5.8 Deferred Candidates + SM(PEG)n Plan

Date: 2026-04-12  
Auditor role: computational simulation scientist, wet-lab chemistry researcher, physicochemical model reviewer  
Scope: `docs/25_v58_deferred_smpeg_plan.md` and the current Module 2 implementation.

## 1. Executive Verdict

The v5.8 plan is scientifically valuable and mostly consistent with the direction of the current M2 system. The deferred chromatography candidates are relevant, and the SM(PEG)n path is a legitimate wet-lab strategy for oriented immobilization of cysteine-bearing proteins on amine-functionalized supports.

The plan should proceed, but not as a simple profile expansion. The SM(PEG)n path is a true state-machine extension. It requires new intermediate ACS types, a new `SPACER_ARM` step, explicit distinction between intermediate spacer/crosslinker attachment and final functional ligand immobilization, maleimide decay handling, buffer/pH constraints, protein-Cys-specific profiles, and M3 model hints. Without these corrections, the UI and downstream process simulation may overstate real-world reliability.

Primary conclusion: accept the plan as a v5.8 architecture roadmap, but revise its acceptance criteria and risk register before implementation. The immediate Tier A/B profile additions are useful, but they are not "zero risk" because they affect UI behavior, q_max mapping, M3 routing, and scientific interpretation.

## 2. Sources Reviewed

Main plan:

- `docs/25_v58_deferred_smpeg_plan.md`

Relevant current implementation:

- `src/emulsim/module2_functionalization/acs.py`
- `src/emulsim/module2_functionalization/reactions.py`
- `src/emulsim/module2_functionalization/modification_steps.py`
- `src/emulsim/module2_functionalization/reagent_profiles.py`
- `src/emulsim/module2_functionalization/orchestrator.py`
- `src/emulsim/visualization/app.py`
- `src/emulsim/visualization/ui_validators.py`
- `src/emulsim/module3_performance/isotherms/imac.py`
- `src/emulsim/module3_performance/isotherms/protein_a.py`

Verification run:

```powershell
python -m pytest -q tests\test_module2_acs.py tests\test_module2_workflows.py tests\test_ui_contract.py
```

Result: 140 passed, 1 non-functional pytest cache warning.

## 3. Current M2 System Role and Structure

M2 is the chemical functionalization bridge between M1 microsphere fabrication and M3 column/process simulation. Its desired function is to convert a base hydrogel microsphere into a functional medium by tracking:

- Native ACS inventories from M1.
- Secondary crosslinking and mechanical reinforcement.
- Hydroxyl activation into epoxide or vinyl sulfone sites.
- Functional ligand or protein coupling.
- Hydrolysis and quenching losses.
- Surface accessibility and steric limitations.
- Functional ligand density passed to M3 through `FunctionalMediaContract`.
- Confidence, hazard, and calibration metadata.

The current implementation already supports five executable step types: secondary crosslinking, activation, ligand coupling, protein coupling, and quenching. It also has a terminal-state ACS accounting model and a v5.7-style reagent library with 25 current profiles, including Q, CM, NTA, butyl, glutathione, heparin, Protein A/G fusion, streptavidin, and three spacer metadata records.

The current code does not yet support:

- `SPACER_ARM` as a modification step.
- `AMINE_DISTAL` or `MALEIMIDE` ACS types.
- SM(PEG)n heterobifunctional linker profiles.
- Dynamic intermediate ACS creation outside activation.
- Protein-Cys-specific coupling profiles.
- UI selection for spacer-arm workflows.

This means the v5.8 SM(PEG)n route is not a minor library update. It is a controlled extension of the M2 state machine.

## 4. Scientific Validity Assessment

### 4.1 Deferred candidate profiles

The proposed deferred profiles are scientifically relevant:

- Protein L expands antibody-fragment and kappa light-chain affinity coverage.
- Concanavalin A and WGA add lectin affinity modes for glycoproteins.
- Octyl extends HIC strength beyond butyl and phenyl.
- TMAE adds a strong anion-exchange option for DVS-activated matrices.
- EDA, PEG-diamine, and BDGE extend spacer-arm options.

This is a sensible next step for a platform intended to simulate chromatographic resin development.

### 4.2 SM(PEG)n route

The planned route is chemically plausible:

```text
OH -> EPOXIDE -> AMINE_DISTAL -> MALEIMIDE -> protein-Cys thioether
```

This resembles established heterobifunctional crosslinker workflows: amine-terminated surface reacts with NHS ester, then maleimide reacts with a thiol-bearing protein. The major value is oriented or semi-oriented protein immobilization when the protein has an accessible engineered cysteine or controlled reduced disulfide.

However, the plan currently treats this route too cleanly. Real wet-lab performance depends on hydrolysis, buffer composition, reaction timing, thiol availability, maleimide ring opening, competing amines/thiols, protein stability, and orientation. These details must be represented as model assumptions or constraints.

## 5. Major Findings

### F1. `SPACER_ARM` is correctly identified as a new step type

Severity: Positive finding  
Area: architecture, ACS state semantics

The plan is correct that `SPACER_ARM` should not be forced into the existing activation step. Activation consumes a native group and creates an activated group. Spacer-arm chemistry consumes an existing activated or intermediate group and creates a new intermediate ACS type on the same bead. That semantic difference matters for ACS conservation and workflow validation.

Recommendation:

Proceed with a new `ModificationStepType.SPACER_ARM`, but implement it with explicit intermediate-role metadata so M3 does not mistake spacer/crosslinker intermediates for final functional ligands.

### F2. Intermediate ACS profiles can break M2-to-M3 ligand-density interpretation

Severity: High  
Area: FunctionalMediaContract, downstream reliability

The plan proposes recording spacer-arm consumption on the source profile as `ligand_coupled_sites` and creating a new product profile such as `AMINE_DISTAL` or `MALEIMIDE`. This is reasonable for local ACS conservation, but dangerous for `FunctionalMediaContract`.

Current contract-building scans ACS profiles for `ligand_coupled_sites` and `ligand_functional_sites`. If spacer intermediates use the same fields as final functional ligands, the contract can accidentally derive density from an intermediate rather than the final protein or ligand.

Recommendation:

Add one of these safeguards before implementing SM(PEG)n:

- `profile_role`: native, activated, spacer_intermediate, heterobifunctional_intermediate, final_ligand.
- `contributes_to_functional_media: bool`.
- A final-ligand pointer in `ModificationResult`.
- Contract logic that uses only the final executable `LIGAND_COUPLING` or `PROTEIN_COUPLING` product, never intermediate `SPACER_ARM` products.

### F3. Maleimide hydrolysis is not the same as soluble reagent hydrolysis

Severity: High  
Area: reaction kinetics, SM(PEG)n realism

The plan correctly notes NHS ester hydrolysis and maleimide ring opening. Current solvers can model hydrolysis of soluble reagent during a coupling step. They do not yet model time-dependent loss of immobilized maleimide ACS after it has been created.

This distinction is important. After SM(PEG)n creates a maleimide-terminated surface, delays, pH, temperature, and buffer exposure can hydrolyze maleimide before protein-Cys coupling. That loss should consume `MALEIMIDE` sites even if no protein is present.

Recommendation:

Add an immobilized-site decay path for `MALEIMIDE`, either as:

- a first-order intermediate decay during `SPACER_ARM` and `PROTEIN_COUPLING`,
- a dedicated aging/delay parameter between steps,
- or a `hydrolyzed_sites` update on the `MALEIMIDE` ACS profile before thiol coupling.

### F4. Protein-Cys coupling needs distinct profiles, not generic protein coupling

Severity: High  
Area: protein immobilization validity

The plan says the final step is `PROTEIN_COUPLING MALEIMIDE + Protein-Cys -> thioether`. Existing Protein A/G-style profiles are random amine-coupling models. Maleimide-thiol coupling is a different chemistry with different pH, kinetics, orientation, and protein-preparation requirements.

Recommendation:

Add dedicated profiles or fields for cysteine-directed coupling:

- `protein_a_cys_coupling`
- `protein_g_cys_coupling`
- `enzyme_cys_coupling` or generic `cys_protein_coupling`
- `reactive_residue="thiol"`
- `requires_reduced_thiol=True`
- `thiol_accessibility_fraction`
- `orientation_retention_multiplier`

Do not reuse generic epoxide-amine protein profiles on `MALEIMIDE` without changing their chemistry class and validity windows.

### F5. Buffer compatibility is missing from the SM(PEG)n risk register

Severity: High  
Area: wet-lab realism, UI safety

NHS ester and maleimide chemistries are strongly affected by buffer and additives. NHS ester reactions are incompatible with primary-amine buffers such as Tris and glycine. Maleimide-thiol coupling is affected by free thiols such as DTT, beta-mercaptoethanol, cysteine, and glutathione. EDTA, TCEP, reducing conditions, oxygen exposure, and protein disulfide reduction state also matter.

Recommendation:

Add buffer/additive compatibility metadata and UI warnings:

- Avoid primary-amine buffers during NHS coupling.
- Avoid free-thiol additives during maleimide coupling unless intentionally quenching.
- Warn for pH above 7.5 during maleimide coupling.
- Track reducing-agent assumptions for protein-Cys availability.
- Require short time between maleimide generation and protein coupling.

### F6. TMAE chemistry is under-specified and may be chemically inconsistent

Severity: High  
Area: ion-exchange chemistry

The plan describes 2-(trimethylamino)ethyl chloride coupling to vinyl sulfone by Michael addition of a "TMAE amine." A quaternary trimethylammonium group is not a normal neutral amine nucleophile. The actual reactive precursor and mechanism must be specified. Otherwise the model may assign a correct final strong-anion-exchanger function to an implausible coupling chemistry.

Recommendation:

Clarify the actual TMAE precursor and nucleophilic group. If the profile is a functional approximation rather than a mechanistic reaction, label it accordingly. The profile must include:

- `reagent_identity`
- actual reactive group
- final installed charged group
- `charge_type="anion"`
- `chemistry_class` consistent with the real reaction

### F7. Concanavalin A and WGA need lectin-specific process assumptions

Severity: High  
Area: affinity model validity, M3 fit

Con A and WGA are not generic affinity proteins. Con A binding depends on carbohydrate specificity, oligomeric state, pH, and divalent metal cofactors. WGA binds GlcNAc and sialic acid motifs and is also multivalent. Elution is typically sugar-competition driven, not simply pH-driven like Protein A.

Recommendation:

Add `binding_model_hint` values such as:

- `lectin_mannose_glucose_affinity`
- `lectin_glcnac_sialic_affinity`

Add metadata for:

- required cofactors for Con A, especially Ca2+/Mn2+ assumptions,
- oligomeric state and pH range,
- carbohydrate eluent species,
- target glycoprotein specificity,
- ranking-only confidence unless calibrated.

M3 should not route these through generic Protein A or generic Langmuir defaults without warning.

### F8. Octyl HIC threshold needs a unit-consistent density conversion

Severity: Medium-High  
Area: HIC realism, UI warning

The plan says the simulator should flag Octyl density above 20 umol/mL resin. Current M2 primarily tracks density in mol/m2 and M3 maps some modes to mol/m3 bed. A threshold in umol/mL resin requires a defined resin volume basis.

Recommendation:

Implement an explicit conversion:

```text
umol/mL resin = functional_density [mol/m2] * accessible_area_per_resin_volume [m2/m3] * 1000
```

Do not use external particle area alone for porous beads unless the warning is explicitly labeled as external-area based. HIC remains `q_max_confidence="not_mapped"` unless a salt-dependent HIC isotherm exists.

### F9. BDGE as a Phase 1 activation profile loses the spacer information

Severity: Medium-High  
Area: spacer realism, downstream coupling

The plan says BDGE can be added as an activation profile and is functionally indistinguishable from ECH in Phase 1. That is only true if the model cares only about epoxide count. It is not true if downstream coupling should benefit from an 18 A spacer.

Recommendation:

If BDGE is added as an activation profile in Phase 1, the product epoxide profile or modification result should carry `spacer_length_angstrom=18` and `spacer_origin="bdge_activation"`. Otherwise the UI may imply a long-spacer activation while downstream protein coupling behaves like direct ECH activation.

### F10. Diamine spacers can crosslink or deactivate sites

Severity: Medium-High  
Area: spacer chemistry

EDA, DAH, and DADPA do not guarantee one-end attachment with the other amine free. Depending on reagent excess, pH, and surface density, diamines can bridge two epoxides, react twice, introduce positive charge, or create nonspecific binding sites.

Recommendation:

For `SPACER_ARM` diamine profiles, add:

- `monofunctional_attachment_fraction`
- `bridging_fraction`
- `distal_group_yield`
- `nonspecific_charge_penalty`
- requirement for high soluble diamine excess to favor monoattachment

The created `AMINE_DISTAL` sites should be `sites_consumed * distal_group_yield`, not automatically equal to all consumed epoxide sites.

### F11. The current validation rules need intermediate-aware ordering

Severity: Medium-High  
Area: workflow reliability

The plan correctly adds workflow rules for `SPACER_ARM`, but Rule 5 says `product_acs` must not already exist. This prevents accidental overwrites, but it may be too strict for legitimate cumulative spacer creation from repeated or scaled steps.

Recommendation:

Use a stricter default in UI workflows, but support backend-safe accumulation when:

- the existing product profile was created by the same compatible upstream chemistry,
- the user explicitly chooses "append to existing intermediate sites,"
- or the workflow is a repeated-batch process.

At minimum, include tests for both rejected duplicates and intentionally allowed accumulation if the latter is supported.

### F12. Profile count bookkeeping is internally inconsistent

Severity: Medium  
Area: implementation planning

The document states:

- v5.7 baseline: 25 profiles.
- Scope: 8 deferred Priority 2 profiles plus SM(PEG)n path.
- WN-0 adds 7 profiles.
- WN-4 adds 7 profiles.
- Final count is described as 39, while one breakdown line says 40.

This is inconsistent. The likely intended counts are:

- Current v5.7: 25 profiles.
- WN-0 immediate profiles: Protein L, Octyl, WGA, TMAE, PEG-diamine metadata, Con A, BDGE = 7.
- EDA metadata/spacer-arm is deferred to architecture phase unless added separately.
- WN-4 spacer-arm/SM(PEG)n profiles: EDA spacer-arm, DADPA spacer-arm, DAH spacer-arm, possibly PEG-diamine spacer-arm, plus 4 SM(PEG)n = 7 or 8 depending on whether PEG-diamine is executable.

Recommendation:

Create one canonical target-count table before implementation. Do not use both 39 and 40 as possible release counts.

### F13. "Zero risk" profile additions are not zero risk

Severity: Medium  
Area: release reliability

WN-0 profile additions are lower-risk than architecture changes, but they still affect UI, validation, downstream contracts, and scientific interpretation. Con A, WGA, Protein L, Octyl, and TMAE need binding hints and UI warnings. BDGE needs spacer-state handling. PEG-diamine must not appear as an executable chemistry in Phase 1 unless `SPACER_ARM` exists.

Recommendation:

Rename WN-0 from "zero risk" to "low architecture risk." Add tests for:

- profile instantiation,
- UI filtering,
- M2 workflow compatibility,
- M2-to-M3 `FunctionalMediaContract` behavior,
- q_max confidence for unsupported modes.

### F14. UI support is underestimated

Severity: Medium  
Area: user reliability

The plan assigns UI hints only a small polish task. In practice, `SPACER_ARM` adds a new workflow class that users can mis-order easily. UI must represent intermediate site creation and consumption clearly.

Recommendation:

Add UI requirements:

- Step type dropdown includes `Spacer Arm` only after backend support exists.
- Valid reagent list is filtered by current available target ACS.
- UI displays created intermediate ACS after each step.
- UI blocks SM(PEG)n before `AMINE_DISTAL` exists.
- UI blocks protein-Cys coupling before `MALEIMIDE` exists.
- UI displays NHS/maleimide buffer and pH warnings.
- UI distinguishes Phase 1 spacer multiplier from Phase 2 executable spacer-arm chemistry.

### F15. M3 isotherm coverage does not yet satisfy all new ligands

Severity: High  
Area: system-level simulation validity

Current M3 has generic Langmuir/competitive models, Protein A pH-dependent affinity, and IMAC competition. It does not yet have dedicated lectin, Protein L, streptavidin-like, GST, heparin, or HIC salt-gradient models sufficient for all v5.8 candidates.

Recommendation:

Every new ligand should carry an M3 support level:

- `mapped_quantitative`: calibrated isotherm exists.
- `mapped_estimated`: basic q_max and generic isotherm only.
- `not_mapped`: material descriptor only.
- `requires_user_calibration`: target-specific binding constants required.

For Con A and WGA, default should be `requires_user_calibration` or `not_mapped` until lectin-specific target/eluent parameters exist.

## 6. Candidate-by-Candidate Audit

| Candidate | Scientific value | Main flaw or uncertainty | Recommended status |
|---|---:|---|---|
| Protein L | High | Kappa-subclass specificity and fragment formats need target metadata | Add as ranking-only affinity profile |
| Con A | High | Tetramer/dimer state, Ca2+/Mn2+, sugar elution, long spacer need explicit assumptions | Add only with lectin-specific warnings and PEG spacer metadata |
| Octyl | Medium-High | Strong HIC may irreversibly bind proteins; threshold needs unit conversion | Add as HIC descriptor, q_max not mapped |
| WGA | Medium-High | Glycan specificity and elution not generic affinity | Add as ranking-only lectin profile |
| EDA | Medium | High risk of bridging and charged nonspecific sites | Add only as spacer metadata until `SPACER_ARM` supports distal yield |
| PEG-diamine Mn 600 | High for large proteins | Average MW distribution and pore exclusion not captured | Add metadata; use macromolecule-accessibility warning for large ligands |
| BDGE | Medium-High | Phase 1 activation loses long-spacer state | Add only if spacer-origin is carried downstream |
| TMAE | High if chemically specified | Current reagent mechanism likely under-specified | Block implementation until reactive precursor is clarified |
| SM(PEG)2/4/12/24 | High | Requires intermediate ACS, maleimide decay, protein-Cys profiles, buffer constraints | Add only after `SPACER_ARM` architecture is implemented |

## 7. Revised Architecture Requirements

### 7.1 Required ACS and step additions

Implement:

- `ACSSiteType.AMINE_DISTAL`
- `ACSSiteType.MALEIMIDE`
- `ModificationStepType.SPACER_ARM`
- `_solve_spacer_arm_step()`
- intermediate-aware workflow validation

### 7.2 Required reagent metadata additions

Add or formalize:

- `reactive_group`
- `distal_group_yield`
- `profile_role`
- `contributes_to_functional_media`
- `buffer_incompatibilities`
- `requires_reduced_thiol`
- `thiol_accessibility_fraction`
- `maleimide_decay_rate`
- `orientation_retention_multiplier`

### 7.3 Required contract changes

`FunctionalMediaContract` should ignore intermediate spacer/crosslinker profiles unless they are final functional ligands. It should also carry:

- `binding_model_hint`
- `m3_support_level`
- `final_ligand_profile_key`
- `final_ligand_density_area_basis`
- `intermediate_profile_keys`
- `workflow_assumption_notes`

## 8. Revised Validation and Acceptance Criteria

Before v5.8 should be considered reliable, the following acceptance criteria should be added:

| Area | Required criterion |
|---|---|
| Baseline | Existing 25 v5.7 profiles and all current 5-step workflows remain backward-compatible |
| Profile count | One canonical final profile count is documented and tested |
| TMAE | Actual reactive precursor and chemistry class are clarified before profile activation |
| BDGE | Spacer origin and length are preserved downstream if used as activation |
| SPACER_ARM | Consumes target profile and creates a new intermediate profile without violating ACS conservation |
| Distal yield | Created intermediate sites are reduced by distal-group yield, not assumed equal to all consumed sites |
| Maleimide | Immobilized maleimide hydrolysis/aging can reduce available `MALEIMIDE` sites |
| Protein-Cys | Thiol-specific protein coupling profiles exist and are not confused with generic amine coupling |
| UI | UI blocks invalid SM(PEG)n ordering and displays intermediate ACS state |
| Buffer warnings | NHS and maleimide buffer/additive incompatibilities are surfaced to users |
| FMC | M2-to-M3 contract uses only final functional ligand density, not spacer intermediates |
| M3 support | Each new ligand declares whether M3 mapping is estimated, unsupported, or requires calibration |
| Tests | Integration tests cover ECH -> DADPA -> SM(PEG)n -> Protein-Cys and invalid ordering cases |

## 9. Reliability Assessment

Current M2 reliability:

- Reliable for the current v5.7 25-profile semi-quantitative workflows.
- Good ACS conservation for existing terminal states.
- Good regression status based on passing M2 tests.
- Not currently capable of executing SM(PEG)n spacer-arm workflows.

v5.8 plan reliability if implemented as written:

- Good conceptual architecture for adding `SPACER_ARM`.
- Moderate scientific reliability for deferred simple profiles.
- Low-to-moderate reliability for Con A, WGA, Protein L, and Octyl process predictions unless M3 target-specific binding assumptions are added.
- Low reliability for SM(PEG)n quantitative predictions unless maleimide decay, buffer constraints, thiol availability, and final-ligand contract logic are implemented.

v5.8 plan reliability after recommended corrections:

- Suitable for semi-quantitative workflow design and comparative ranking.
- Suitable for teaching and exploratory wet-lab planning.
- Still not a calibrated replacement for wet-lab optimization.
- Stronger as an M2-to-M3 bridge because unsupported process modes can be clearly marked instead of forced into generic isotherms.

## 10. Final Recommendation

Proceed with v5.8 in two phases, but revise the plan before implementation:

1. Ship low-architecture-risk deferred profiles only after profile-count, UI-filtering, and M3-support metadata are corrected.
2. Treat SM(PEG)n as a full architecture extension, not a profile-only addition.
3. Add `SPACER_ARM`, `AMINE_DISTAL`, and `MALEIMIDE`, but prevent intermediate profiles from contaminating final ligand-density calculations.
4. Add maleimide hydrolysis/aging and buffer compatibility rules before exposing SM(PEG)n in the UI.
5. Add protein-Cys-specific coupling profiles and do not reuse generic amine-coupled protein profiles without chemistry changes.
6. Mark lectin, HIC, and oriented-protein outputs as semi-quantitative or ranking-only unless calibrated to real resin data.

With these changes, the v5.8 plan is scientifically defensible and aligned with the system's intended simulation functions. Without these changes, the plan risks producing outputs that are chemically plausible in appearance but not reliable enough for real wet-lab or downstream process decision-making.
