# Audit Report: UI Full Alignment Plan for the EmulSim Three-Module System

Date: 2026-04-11

Audited plan: `docs/14_ui_full_alignment_plan.md`

Audited implementation areas:

- `src/emulsim/visualization/app.py`
- `src/emulsim/visualization/plots.py`
- `src/emulsim/pipeline/orchestrator.py`
- `src/emulsim/datatypes.py`
- `src/emulsim/module2_functionalization/`
- `src/emulsim/module3_performance/`
- `tests/test_module2_acs.py`
- `tests/test_module2_workflows.py`
- `tests/test_module3_breakthrough.py`
- `tests/test_module3_catalysis.py`
- `tests/test_module3_multicomponent.py`

## 1. Executive Verdict

The UI full-alignment plan is directionally correct and scientifically necessary. Its strongest decisions are the three-module pipeline scope selector, extraction of the monolithic Streamlit app into module-specific files, output-level confidence labels, and a trust cascade from Module 1 fabrication through Module 2 functionalization into Module 3 process performance.

However, the plan is not yet sufficient as a scientific UI specification. It describes a user interface that can expose Module 2 and Module 3 functions, but it does not fully specify the guardrails required to prevent users from interpreting provisional or semi-empirical model outputs as validated wetlab or process-scale predictions. The current backend contains real Module 2 and Module 3 packages, but the current UI remains Module 1-only. This creates a large integration gap: the planned UI is broadly aligned with the system goal, while the implemented UI does not yet satisfy the designed three-module simulation functions.

Overall assessment:

| Area | Rating | Audit conclusion |
|---|---:|---|
| Plan architecture | Good | Modular UI and pipeline scope design are appropriate. |
| Current UI coverage | Poor | Current Streamlit app exposes only Module 1/L1-L4 and no Module 2/3 workflows. |
| Scientific UI validity | Conditional | Valid if every output is labeled by model tier and invalid chemistry/process combinations are blocked. |
| Module 1 UI readiness | Moderate | Useful for screening but several input and guidance defects can mislead users. |
| Module 2 UI readiness | Moderate-low | Backend has ACS and two basic workflows, but UI must not imply broad chemistry coverage. |
| Module 3 chromatography UI readiness | Moderate | Breakthrough and gradient models exist, but gradient causality and mass balance warnings need explicit UI treatment. |
| Module 3 catalysis UI readiness | Low-medium | Model exists, but verification timed out and computational performance needs attention before interactive UI exposure. |
| Reliability for decision support | Low without guardrails | Suitable for education, ranking, and hypothesis generation; not ready for protocol-final or scale-up decisions. |

Primary recommendation: implement the plan as a contract-first scientific UI, not as a generic parameter dashboard. The UI must expose only the workflows that the backend can currently defend, show trust/confidence next to every computed number, block unsupported chemistries, display numerical residuals such as mass-balance error, and mark unvalidated outputs as "ranking only" or "semi-quantitative."

## 2. Evidence Summary

Current UI evidence:

- `src/emulsim/visualization/app.py` is still a monolithic Streamlit application with only Module 1 controls and L1-L4 result tabs.
- The result tab declaration at `src/emulsim/visualization/app.py:596` creates only five tabs: Dashboard, L1, L2, L3, and L4.
- The current UI has no imports or calls to `module2_functionalization`, `module3_performance`, `ModificationOrchestrator`, `run_breakthrough`, `run_gradient_elution`, or `solve_packed_bed`.
- `src/emulsim/visualization/plots.py` contains only Module 1/L1-L4 plotting functions and no `plots_m2.py` or `plots_m3.py`.
- The plan in `docs/14_ui_full_alignment_plan.md` correctly proposes `sidebar_m2.py`, `sidebar_m3.py`, `results_m2.py`, `results_m3.py`, `plots_m2.py`, and `plots_m3.py`, but these files are not yet present.

Backend evidence:

- `M1ExportContract` exists in `src/emulsim/datatypes.py:913` and contains geometry, pore structure, network structure, reactive group concentrations, mechanical properties, formulation, and trust metadata.
- `export_for_module2()` exists in `src/emulsim/pipeline/orchestrator.py:282` and computes residual NH2/OH concentrations and trust metadata for Module 2.
- Module 2 includes ACS profiles, surface-area models, reagent profiles, sequential modification workflows, and a `ModificationOrchestrator`.
- Module 3 includes packed-column hydrodynamics, single-component breakthrough, multi-component gradient elution, UV/MS/fluorescence/conductivity detection utilities, isotherms, and catalytic packed-bed simulation.

Verification evidence:

- `python -m pytest -q tests\test_module2_acs.py tests\test_module2_workflows.py tests\test_module3_breakthrough.py` passed: 76 tests passed in 34.89 s.
- The same run emitted LRM mass-balance warnings above the 2 percent threshold in several breakthrough tests, with examples around 2.48 percent, 3.60 percent, and 4.66 percent.
- `python -m pytest -q tests\test_module3_multicomponent.py` passed: 60 tests passed in 0.39 s.
- `python -m pytest -q tests\test_module3_catalysis.py` did not complete within 180 s. It progressed through most tests but did not produce a final pass/fail result in the available time.
- The full test suite collected 404 tests but did not complete within 180 s.

## 3. Current UI vs Planned UI

The current UI is not aligned with the three-module system. It remains a Module 1 fabrication UI with L1 emulsification, L2 gelation, L3 primary crosslinking, L4 mechanical outputs, optimization, and trust warnings. It does not yet expose Module 2 chemical modification/crosslinking or Module 3 chromatography/catalysis.

The plan correctly identifies the need for:

- A pipeline scope selector: M1, M1+M2, or M1+M2+M3.
- M2 controls for modification steps, ACS budget preview, and functionalized microsphere outputs.
- M3 controls for chromatography and catalysis.
- M2 plots: ACS waterfall, surface area, modification timeline.
- M3 plots: chromatograms, breakthrough curves, pressure diagnostics, enzyme effectiveness, and deactivation.
- Output confidence labels and a trust cascade.

The missing piece is a strict mapping from UI controls to defensible backend model contracts. A scientifically safe UI must not expose broad user choices simply because they are desirable from the initial product vision. It must expose only combinations for which the current solver has a clear mechanism, input range, validation state, and confidence level.

## 4. Critical and High-Severity Findings

### F1. The current UI is Module 1-only while the project now has partial Module 2/3 backends

Severity: Critical

The plan is written as an alignment plan for three modules, but the implemented UI does not yet call the available Module 2/3 backends. This means the current user experience cannot simulate functionalization, affinity chromatography, ion exchange, IMAC, protein ligand immobilization, catalytic beds, or process performance.

Scientific impact:

- Users cannot propagate fabrication uncertainty into downstream functionalization or column performance.
- The current UI cannot answer the central next-phase question: how preparation conditions affect functional microsphere performance in a realistic process environment.

Required solution:

- Implement a contract-first UI integration using `M1ExportContract` as the only M1-to-M2 input.
- Add M2/M3 tabs only after the corresponding result objects are created, validated, and confidence-labeled.
- Until then, show M2/M3 as "not yet connected" rather than hiding the gap.

### F2. The plan's confidence-label concept is correct, but must be applied to every individual output

Severity: Critical

The plan says every numerical output needs a confidence label. This is essential. The current implementation shows a global trust assessment at the bottom of the app, but individual outputs in the Dashboard, L1, L2, L3, and L4 tabs are not consistently labeled as mechanistic, empirical, ranking-only, semi-quantitative, or unvalidated.

Scientific impact:

- A user may treat a phenomenological modulus, empirical pore size, inferred ACS density, simulated chromatogram, or MS total-ion trace as equally reliable.
- This is especially dangerous for process development, where column scale-up decisions depend on pressure drop, binding capacity, mass transfer, and detector calibration.

Required solution:

- Add an output metadata object with fields such as `model_basis`, `confidence`, `calibration_required`, `validity_range`, `warnings`, and `recommended_use`.
- Render this metadata next to every key number and plot trace.
- Do not use the word "mechanistic" for a UI output unless the displayed value is produced by a model that actually uses the relevant physical or chemical mechanism.

### F3. Module 2 chemistry UI must not imply broad functionalization coverage

Severity: Critical

The current Module 2 backend supports a useful but limited set of workflows:

- Amine secondary crosslinking using genipin or glutaraldehyde.
- Hydroxyl activation using epichlorohydrin or divinyl sulfone.
- ACS state tracking for total, accessible, activated, consumed, blocked, coupled, and functional sites.

It does not yet implement a fully general coupling workflow for IMAC ligand immobilization, ion-exchange ligand installation, protein ligand coupling, quenching, washing, spacer chemistry, metal charging, ligand activity retention, or ligand leaching.

Scientific impact:

- A UI that allows arbitrary 1-5 step workflows could appear to support wetlab protocols that are not implemented.
- ACS conservation alone is not enough. Real protocols also require pH compatibility, reagent stability, hydrolysis, ligand steric exclusion, reaction selectivity, quench/wash losses, residual toxicity, and activity retention.

Required solution:

- Restrict the first UI release to named, implemented workflows.
- Label implemented workflows as "supported" and future workflows as "planned, disabled."
- Require a chemistry compatibility graph before users can chain arbitrary steps.
- For each step, show which ACS pool is consumed, which ACS pool is created, and what model terms are ignored.

### F4. Module 2 pH and temperature controls are currently scientifically misleading unless backend behavior changes

Severity: High

`ModificationStep` includes `temperature` and `ph`, and reagent profiles include `ph_optimum` and activation energies. However, the current workflow calls `solve_second_order_consumption()` with `k0=0.0`, which causes the solver to use `k_forward` directly and bypass Arrhenius temperature dependence. The `ph` field is not used in the reaction rate calculation. Hydrolysis is represented by a fixed first-order rate constant rather than a pH-dependent expression.

Scientific impact:

- If the UI allows users to vary pH and temperature and shows changed outputs only weakly or not at all, users may infer false mechanistic control.
- ECH and DVS activation chemistry is strongly pH-dependent in wetlab practice. Glutaraldehyde Schiff-base formation and genipin crosslinking also depend on pH, temperature, amine availability, and buffer composition.

Required solution:

- Either implement pH/temperature-dependent rate corrections before exposing those controls as mechanistic inputs, or label them as metadata/validation-only fields.
- Display an explicit warning: "Current Module 2 rates are reference-condition screening estimates; pH is used for compatibility warnings only unless pH-rate models are enabled."
- Bind `step.stoichiometry` to `reagent_profile.stoichiometry` by default and prevent accidental user mismatch unless advanced mode is enabled.

### F5. Module 3 gradient elution can be visually misleading because the current competitive Langmuir gradient is not causally coupled to binding

Severity: High

The gradient elution code records `gradient.value_at_time(t)`, but in the competitive Langmuir path the gradient value is effectively reserved for future SMA or pH-dependent models. The inlet concentration remains `C_feed.copy()` inside the RHS, and the gradient does not alter adsorption constants for the default competitive Langmuir model.

Scientific impact:

- A UI plot showing a gradient overlay on a chromatogram can imply that salt, pH, or imidazole causes elution.
- In the current default competitive Langmuir implementation, the gradient display may be diagnostic only, not a physical driver of elution.

Required solution:

- In the UI, disable "gradient-driven elution" wording for default competitive Langmuir.
- For true gradient elution, route ion-exchange cases through SMA or another isotherm where salt affects binding, and route IMAC cases through a model where imidazole concentration affects protein binding.
- Add a per-run badge: "Gradient affects binding: yes/no."

### F6. Module 3 process outputs must expose mass-balance and pressure checks as first-class UI results

Severity: High

The breakthrough model computes mass-balance error and pressure drop. Tests passed but emitted warnings when mass-balance error exceeded the 2 percent threshold. The UI plan mentions pressure and process plots, but it does not explicitly require a blocking or caution status when mass-balance error is too high.

Scientific impact:

- A simulated chromatogram with poor mass balance is not reliable for capacity, yield, or peak-shape interpretation.
- Pressure drop and bed compression affect column operability and can invalidate predicted chromatographic performance.

Required solution:

- Display mass-balance error beside every breakthrough or gradient result.
- Gate output trust: `<=2%` acceptable, `2-5%` caution, `>5%` unreliable unless justified.
- Show pressure drop, maximum safe flow rate, bed compression fraction, particle Reynolds number warning, and whether assumptions are violated.
- Do not treat pressure warnings as logs only; surface them in the main UI.

### F7. Current Module 1 UI advice can still recommend scientifically unsafe actions

Severity: High

The current UI recommendation logic can tell the user to increase RPM when droplet size is above target. Prior project audit identified nonphysical or unstable RPM trends in parts of the L1 PBE behavior. The current UI also hardcodes DDA, lacks a true gelation-temperature input, allows broad polymer concentration ranges, and has chemistry-agnostic crosslinker concentration controls.

Scientific impact:

- Recommendation text can be interpreted as process guidance.
- If the underlying model is nonmonotonic or outside its calibration regime, simple "increase/decrease" advice can be wrong.

Required solution:

- Replace deterministic recipe advice with sensitivity-tested advice.
- Only recommend changing a parameter if local perturbation simulations confirm the expected direction under the selected model.
- Label recommendations as "screening suggestions" unless validated against experiments.

## 5. Module 1 UI Scientific Audit

### Strengths

- The UI exposes key fabrication variables: agitation, time, dispersed phase fraction, polymer concentration, surfactant, temperature, cooling rate, crosslinker concentration, crosslinking time, and crosslinking temperature.
- It includes hardware mode distinctions for rotor-stator and stirred-vessel operation.
- It includes material constants and trust assessment, which are essential for a simulation platform intended to bridge process inputs and physicochemical outputs.

### Deficiencies

1. DDA is hardcoded in the UI. `src/emulsim/visualization/app.py:228` sets `_DDA = 0.90`, while Module 2 ACS and residual amine calculations depend directly on DDA. In real chitosan systems, DDA varies by supplier and batch and strongly affects NH2 density, crosslinking capacity, charge state, swelling, and ligand coupling capacity.

2. Gelation temperature is not exposed as a primary parameter. Agarose gelation is not governed only by oil temperature and cooling rate. Thermal history, gel point, droplet size, and heat-transfer conditions determine when phase separation and network locking occur.

3. Dispersed phase fraction is insufficiently guarded in stirred-vessel mode. The UI computes `phi_d` from oil and polymer solution volumes, allowing high internal phase fractions where coalescence, viscosity, impeller flooding, inversion risk, and PBE assumptions may fail.

4. Polymer concentration ranges are too broad for unrestricted quantitative prediction. Agarose at 1-10 percent w/v and chitosan at 0.5-5 percent w/v include regimes with major viscosity, mixing, dissolution, gelation, and mass-transfer differences. Without calibration, broad ranges should be screening-only.

5. Crosslinker UI uses a generic concentration slider from 0.5 to 500 mM. Realistic and safe concentration windows differ strongly among genipin, glutaraldehyde, ECH, DVS, TPP, PEGDA, and citric-acid chemistries. Toxicity, volatility, hydrolysis, pH, and safety constraints differ.

6. Stoichiometry guidance in the current UI is amine-centered. It computes crosslinker/NH2 guidance even when selected chemistry may be hydroxyl-reactive or otherwise not governed by chitosan primary amines.

7. Cahn-Hilliard L2 mode is exposed through grid options without enough runtime and reliability protection. The UI must prevent long blocking computations in interactive mode and must disclose when the empirical pore correlation is being used instead.

### Required UI corrections for Module 1

- Add DDA as a first-class, batch-specific input with a recommended range and source/calibration note.
- Add gelation temperature or gel point inputs and clarify how they feed L2.
- Add hard validation for `phi_d` regimes and flag likely inversion/coalescence zones.
- Split formulation ranges into validated, caution, and exploratory zones.
- Replace generic crosslinker concentration with chemistry-specific ranges, units, pH windows, safety notes, and ACS target type.
- Make every recommendation depend on a local sensitivity calculation rather than fixed heuristics.

## 6. Module 2 UI Scientific Audit

### What the plan gets right

The plan correctly recognizes ACS as the bridge between primary microspheres and functional microspheres. The ACS budget preview, ACS waterfall, surface-area breakdown, and modification timeline are scientifically appropriate UI primitives. They map well to the backend concepts of `M1ExportContract`, `AccessibleSurfaceModel`, `ACSProfile`, `ModificationStep`, `ModificationResult`, and `FunctionalMicrosphere`.

### Main scientific gaps

1. ACS density cannot be treated as a single scalar. It must be reported by site type, accessibility tier, solute size, pore/matrix accessibility, and uncertainty. Small-molecule reagent accessibility is not the same as protein ligand accessibility.

2. The surface-area model is currently a useful first approximation but not a validated pore-network model. It uses external area, porosity, mean pore diameter, tortuosity, and a steric accessibility function. Real hydrogels have pore-size distributions, swelling changes, polymer-rich domains, dead-end pores, and reaction-induced pore blockage.

3. The `EXTERNAL_ONLY` tier needs especially careful UI treatment. It is labeled unreliable, but current ACS initialization can make all bulk sites appear accessible in that degenerate tier. The UI should not present external-only ACS as a realistic coupling capacity without a severe warning.

4. Module 2 currently supports two workflow families, not a complete functionalization chemistry platform. Secondary crosslinking and OH activation are useful, but ligand coupling and protein coupling are not implemented as production workflows in the current Module 2 step solver.

5. Functional ligand performance is not determined only by the number of coupled ligands. Orientation, spacer length, steric occlusion, activity retention, metal charging efficiency, ionic strength, pH, ligand leaching, and nonspecific adsorption can dominate real-world performance.

### Required Module 2 UI behavior

- Build the step builder from a backend chemistry registry, not free-form UI choices.
- For each reagent, show target ACS, product ACS, supported step type, pH validity window, temperature validity window, solvent/buffer assumptions, hydrolysis/side-reaction status, toxicity note, and confidence.
- Disable `LIGAND_COUPLING`, `PROTEIN_COUPLING`, and `QUENCHING` as executable steps until production solvers exist.
- Show ACS conservation after every step.
- Show accessible area separately for reagent and ligand.
- Show the difference between total ACS, accessible ACS, activated ACS, consumed ACS, blocked ACS, coupled ligand, and functional ligand.
- Propagate Module 1 trust level and uncertainty notes into Module 2 outputs.
- Mark ligand density and functional capacity as "not predicted" unless actual ligand-coupling and activity-retention models are selected.

## 7. Module 3 UI Scientific Audit

### Chromatography

The current Module 3 chromatography backend has a credible skeleton for process simulation: column geometry, Kozeny-Carman pressure drop, lumped-rate transport, Langmuir/competitive Langmuir/SMA/IMAC isotherm components, UV detection, detector broadening, breakthrough metrics, peak metrics, and mass-balance diagnostics.

The UI plan is aligned with the intended process-science direction, but it needs stronger scientific constraints.

Required chromatography UI controls:

- Column diameter, bed height, bed volume, bed porosity, particle diameter, particle porosity, compression/modulus inputs inherited from M1/M2.
- Flow rate in both absolute units and column volumes per hour.
- Residence time, interstitial velocity, superficial velocity, Reynolds number, Peclet assumption, and pressure drop.
- Feed concentration, loading duration, wash/elution/strip durations, sample volume, and sample matrix composition.
- Isotherm model selection with parameter provenance: Langmuir, competitive Langmuir, IMAC competition, SMA for IEX.
- Calibration status for `q_max`, `K_L`, kinetic adsorption rate, diffusivity, and extinction coefficients.
- Detector type and calibration: UV, fluorescence, conductivity, MS/TIC/EIC.
- Mass-balance error and numerical solver status.

Required warnings:

- Default Langmuir and competitive Langmuir parameters are not universal resin properties.
- Gradient overlays are not necessarily mechanistic drivers unless the selected isotherm consumes the gradient variable.
- UV absorbance is quantitative only if extinction coefficient, path length, baseline, and linear range are valid.
- MS TIC/EIC outputs are semi-quantitative unless ionization efficiency and calibration are specified.
- Pressure drop is computed, but dynamic bed compression is not fully coupled back into chromatographic transport.

### Catalysis

The catalytic packed-bed model includes Michaelis-Menten kinetics, internal effectiveness factor, generalized Thiele modulus, axial dispersion, first-order deactivation, and transient packed-bed transport. This is a scientifically meaningful starting point.

However, UI readiness is lower than chromatography readiness because the catalytic test file did not complete within the timeout during this audit. The model also requires careful parameter interpretation.

Required catalysis UI controls:

- Enzyme loading or immobilized activity basis and conversion to `V_max`.
- `K_m`, `V_max`, `D_eff`, substrate concentration, stoichiometry, deactivation rate, temperature, pH, and buffer composition.
- Bed geometry, particle diameter, bed porosity, particle porosity, flow rate, residence time, pressure drop, and Reynolds/Peclet diagnostics.
- Output plots for substrate outlet, product outlet, conversion, activity history, effectiveness factor, Thiele modulus, productivity, and mass-balance error.

Required warnings:

- `V_max` and `K_m` are not transferable from free enzyme to immobilized enzyme without calibration.
- Effectiveness factor depends on internal diffusivity, particle size, and local substrate concentration.
- Temperature and pH effects are not complete unless kinetic and deactivation models include them explicitly.
- Product inhibition, substrate inhibition, enzyme leaching, multiphase effects, and non-1:1 stoichiometry are not covered by the current basic model.

## 8. UI Reliability Audit

### Positive reliability features in the plan

- Modular app decomposition should reduce maintenance risk.
- Session-state keys for M1, M2, and M3 results are appropriate.
- The trust cascade concept is necessary and should be preserved.
- The plan recognizes validation gates across modules.
- The planned M2/M3 plot files are the right separation of concerns.

### Reliability risks not fully resolved by the plan

1. Invalidation cascade needs to be strict. Any change in M1 inputs must invalidate M2 and M3 results. Any change in M2 steps must invalidate M3. UI state must not silently reuse stale downstream outputs.

2. Backend exceptions need user-facing scientific messages. A `KeyError` for unsupported reagent or a `ValueError` for missing ACS should become a clear UI block: "This chemistry cannot be simulated for the current microsphere because the required ACS is absent."

3. Long-running solvers need timeout, progress display, and safe cancellation. L2 Cahn-Hilliard and Module 3 catalysis are especially risky in interactive Streamlit.

4. Numerical quality metrics must be rendered, not hidden. Solver convergence, mass-balance error, grid size, time step, and warnings should be part of the result summary.

5. The UI must avoid stale confidence labels. Confidence should be recomputed from the actual model path and result diagnostics, not assigned from the selected tab.

6. Run reproducibility is under-specified. The UI should export a complete run manifest containing all inputs, model versions, parameter sources, confidence labels, warnings, and hashes.

## 9. Assessment Against the Initial Three-Module Requirement

### Module 1 requirement: double-emulsification microsphere preparation

Current UI satisfaction: Partial

The current UI supports the existing fabrication simulation and provides useful screening outputs. It needs repairs before it can serve as a reliable upstream source for Module 2. DDA, gelation temperature, phase fraction validity, polymer range validity, chemistry-specific crosslinker constraints, and sensitivity-based recommendations are required.

### Module 2 requirement: chemical modification and/or crosslinking of primary microspheres

Current UI satisfaction: Not implemented

Plan satisfaction: Partial

The plan addresses the right concepts: ACS density, crosslinking density, secondary crosslinking, ACS modification, functional ligand coupling, and ACS consumption tradeoffs. The backend currently supports ACS bookkeeping and two initial chemistry families. The UI plan must narrow the first release to supported chemistry and clearly distinguish simulated quantities from planned future chemistry.

### Module 3 requirement: process performance in purification and catalytic columns/beds

Current UI satisfaction: Not implemented

Plan satisfaction: Partial

The plan covers chromatography and catalysis at a high level. The backend includes meaningful first-generation models, but the UI must expose pressure, mass transfer, binding/isotherm calibration, detector calibration, mass balance, and solver reliability. It must also prevent gradient and detector plots from implying more mechanistic validity than the current selected backend provides.

## 10. Recommended Revised Implementation Strategy

### Phase 0: Scientific UI contract layer

Create a UI contract layer before building visual panels:

- `ui_model_metadata.py`: confidence labels, model basis, validity ranges, recommended use.
- `ui_validators.py`: module-specific hard blockers and warnings.
- `ui_state.py`: result hashes and invalidation cascade.
- `ui_units.py`: unit conversions and display formatting.

Acceptance criteria:

- Every displayed scalar and plot has confidence metadata.
- M1 input changes invalidate M2/M3.
- M2 input changes invalidate M3.
- Unsupported workflows are disabled, not hidden as if implemented.

### Phase 1: Repair Module 1 UI

Implement the plan's Module 1 fixes first:

- User-adjustable DDA.
- Gelation temperature or gel point handling.
- Hard `phi_d` validity gates.
- Validated/caution/exploratory ranges for polymer concentrations.
- Chemistry-specific crosslinker controls and stoichiometry.
- Sensitivity-tested recommendations.

Acceptance criteria:

- M1 export contract shows DDA source, ACS uncertainty, and trust metadata.
- Current M1 outputs clearly distinguish calibrated, empirical, and ranking-only values.

### Phase 2: Minimal Module 2 UI for supported workflows only

Expose only:

- Genipin secondary amine crosslinking.
- Glutaraldehyde secondary amine crosslinking.
- ECH hydroxyl activation.
- DVS hydroxyl activation.

Required M2 outputs:

- Surface-area breakdown.
- ACS state table by site type.
- ACS waterfall per step.
- Conversion per step.
- Remaining ACS.
- Updated modulus for secondary crosslinking.
- Conservation violations and uncertainty notes.

Acceptance criteria:

- Unsupported `LIGAND_COUPLING`, `PROTEIN_COUPLING`, and `QUENCHING` are visible as planned but disabled.
- pH and temperature are either mechanistically implemented or shown as validation metadata only.

### Phase 3: Minimal Module 3 chromatography UI

Expose single-component breakthrough first:

- Column geometry.
- Flow rate and CV/h.
- Langmuir parameters with calibration labels.
- UV detector parameters.
- Breakthrough curve and DBC metrics.
- Pressure drop, compression, mass-balance error, solver status.

Acceptance criteria:

- Mass-balance error is always displayed and gates trust.
- Pressure warnings are visible in the main result panel.
- Default isotherm parameters are labeled as illustrative unless user-calibrated.

### Phase 4: Multi-component and gradient UI

Expose gradient workflows only when the selected isotherm actually uses the gradient variable.

Acceptance criteria:

- Competitive Langmuir gradient overlay is labeled diagnostic only unless binding depends on the gradient.
- IMAC uses imidazole as a modeled competitor.
- IEX uses SMA or another salt-dependent isotherm.
- Peak purity/yield/resolution are labeled with mass-balance and calibration status.

### Phase 5: Catalytic packed-bed UI

Expose catalysis after performance profiling and test completion.

Acceptance criteria:

- Catalysis tests complete reliably in the target interactive environment.
- UI displays Thiele modulus, effectiveness factor, conversion, activity decay, productivity, and mass-balance error.
- Temperature and pH effects are not exposed as mechanistic unless implemented.

## 11. Final Scientific Judgment

The UI plan is scientifically valuable and should be implemented, but only with stricter model-aware constraints. The current plan can become a strong three-module interface if it treats the UI as a scientific instrument panel rather than a generic control form.

The most important change is to make the UI honest about model confidence. The project contains a mix of mechanistic equations, empirical correlations, semi-quantitative detector models, default illustrative constants, and unvalidated extrapolations. A professional wetlab/process-development user must be able to see that distinction immediately.

Current readiness:

- Module 1 UI: usable for screening after known repairs.
- Module 2 UI: ready for a restricted ACS and two-workflow prototype.
- Module 3 chromatography UI: ready for restricted breakthrough and cautiously labeled multi-component visualization.
- Module 3 catalysis UI: not ready for routine interactive use until tests complete and runtime is controlled.

Recommended go/no-go decision:

- Go for phased implementation with confidence labels, hard validators, and disabled unsupported workflows.
- No-go for a broad "full system" UI that presents all planned chemistry and process simulations as equally valid.
