# EmulSim Audit Report

Date: 2026-04-11

## Scope

This audit reviewed the EmulSim codebase as a scientific simulation system and as an experiment-support tool. The audit covered:

- Architecture and documentation claims in `README.md` and `docs/`
- Core implementation in `src/emulsim/`
- Test behavior in `tests/`
- Targeted runtime checks on the default pipeline and selected parameter sweeps

The goal was not to judge whether the code is well-organized in a software-only sense. The goal was to assess:

- Scientific validity
- Consistency with real physicochemical behavior
- Numerical reliability
- Practical usefulness for experimental planning
- Defects that can make outputs misleading even when the code executes cleanly

## Executive Summary

EmulSim is a serious and structured prototype, not a trivial mockup. It contains real mathematical models, a coherent L1-L4 pipeline, chemistry-specific branching, uncertainty propagation, and Bayesian optimization. It is clearly beyond a toy application.

However, in its current state it should **not** be treated as a quantitatively predictive simulator for laboratory decision-making without additional calibration and correction work.

The main reasons are:

1. The default Level 2 pore model is empirical, synthetic, and only weakly coupled to the upstream physics.
2. The default Level 1 pipeline shows a nonphysical RPM-to-droplet-size trend in direct audit runs.
3. The diffusion-limited Level 3 fallback can silently apply an amine/genipin-like model to the wrong chemistry for large droplets.
4. The mechanistic Level 2 solvers are too slow for routine use in the full pipeline and are therefore not the practical basis of the advertised workflow.
5. Optimization and uncertainty modules run, but part of their behavior is built on surrogate assumptions and post hoc perturbations rather than end-to-end mechanistic propagation.

## Overall Assessment

### Bottom line

- Good for: exploratory screening, software prototyping, UI-driven scenario comparison, rough trend ranking within a narrow calibrated regime
- Not yet good for: defensible quantitative prediction of real bead morphology, chemistry-specific crosslinking behavior under all supported chemistries, or automated optimization to drive experiments without human review

### Readiness judgement

| Aspect | Assessment |
|---|---|
| Software structure | Good |
| Numerical stability of default path | Moderate |
| Scientific fidelity of default path | Moderate to weak |
| Mechanistic completeness | Mixed |
| Predictive reliability | Limited |
| Experimental decision support readiness | Not yet sufficient |

## Audit Method

The audit used four kinds of evidence:

1. Code inspection of the main pipeline and all four levels.
2. Test execution.
3. Direct runtime experiments with the current default settings.
4. First-principles sanity checks against expected physicochemical behavior.

### Test observations

- `tests/test_level1_emulsification.py`: passed, but took about 64 s.
- `tests/test_level3_crosslinking.py`: passed.
- `tests/test_level4_mechanical.py`: passed.
- `tests/test_properties.py`: passed.
- `tests/test_data_layer.py`: passed.
- `tests/test_physics_layer.py`: passed.
- `tests/test_optimization.py`: passed, with BoTorch numerical warnings.
- `tests/test_trust_warning_fixes.py`: passed.
- `tests/test_level2_gelation.py`: did **not** finish within 10 minutes in audit runs.
- Full `python -m pytest -q`: did **not** finish within 10 minutes.

This matters scientifically. The Level 2 mechanistic path is not just computationally expensive in theory; it is operationally heavy enough to block end-to-end routine validation.

## Direct Runtime Findings

### Default pipeline

Default pipeline run produced approximately:

- `d32 = 2.81 um`
- `pore = 160.4 nm`
- `p_final = 0.0398`
- `G_DN = 70.8 kPa`

This is internally consistent in the code, but the interpretation is important:

- The pore value comes from the empirical Level 2 model, not from a resolved mechanistic phase-field prediction.
- The modulus is dominated by the agarose power-law contribution, not by a deeply resolved double-network mechanics model.

### RPM sweep

Audit runs with default settings gave:

| RPM | d32 (um) |
|---|---:|
| 3000 | 2.573 |
| 6000 | 1.627 |
| 10000 | 2.810 |
| 15000 | 3.703 |
| 25000 | 4.350 |

This is a major red flag. Once RPM exceeds 6000, the predicted droplet size increases instead of decreasing. For a rotor-stator emulsification model, that is the wrong qualitative trend unless some specific secondary mechanism has been clearly justified and validated. No such calibration evidence exists in the codebase.

### Cooling-rate sweep

Audit runs gave:

| Cooling rate (K/s) | Pore size (nm) |
|---|---:|
| 0.033 | 221.9 |
| 0.167 | 160.4 |
| 0.333 | 139.8 |

This trend is qualitatively plausible, but it is produced by the empirical pore correlation, not by the mechanistic CH solver in the default path.

### Crosslinking-time sweep

Audit runs gave:

| Crosslinking time | p_final |
|---|---:|
| 0.1 h | 0.00644 |
| 1 h | 0.03278 |
| 6 h | 0.03979 |
| 24 h | 0.03979 |
| 48 h | 0.03979 |

This behavior is chemically plausible for a crosslinker-limited system: the default formulation saturates by about 6 h, and increasing time further does not help.

### Mechanistic Level 2 timing

Audit timing results:

- Full pipeline with `l2_mode='empirical'`: about 5.1 s
- Full pipeline with `l2_mode='ch_1d'`: did not finish within 300 s
- Full pipeline with `l2_mode='ch_2d'` and `l2_n_grid=32`: did not finish within 300 s
- Full pipeline with `l2_mode='ch_ternary'` and `l2_n_grid=32`: about 8.4 s

Direct solver timing:

- `CahnHilliard2DSolver(N_grid=16)`: about 0.48 s, but with `grid_spacing = 93.75 nm`, `pore ~ 322 nm`, `char_wavelength ~ 600 nm`

This shows strong resolution sensitivity and confirms that the practical pipeline is using the empirical L2 path because the mechanistic route is too expensive.

## Primary Findings

### 1. High severity: Level 1 default pipeline gives nonphysical RPM behavior

**Finding**

The default Level 1 pipeline does not preserve the expected monotonic trend of smaller droplets at higher RPM. In audit runs, `d32` decreased from 3000 to 6000 RPM, then increased from 6000 to 25000 RPM.

**Why this matters**

This is not a small calibration issue. RPM is the dominant control variable in rotor-stator emulsification. If the sign of the trend is wrong in the default operating range, optimization, sensitivity analysis, and experiment recommendations become untrustworthy.

**Evidence**

- Kernel configuration explicitly states that `breakage_C3` was changed to restore monotonic RPM-to-`d32` behavior: [`src/emulsim/datatypes.py:312`](../src/emulsim/datatypes.py#L312)
- The current test suite no longer checks the PBE solver for monotonic RPM behavior; it only checks the Hinze scaling directly and explicitly avoids the PBE because of solver sensitivity: [`tests/test_level1_emulsification.py:180`](../tests/test_level1_emulsification.py#L180)

**Assessment**

The code is numerically runnable, but the default Level 1 calibration is not scientifically dependable as a predictive emulsification model.

### 2. High severity: The default Level 2 model is empirical, synthetic, and only weakly mechanistic

**Finding**

The pipeline default is `l2_mode='empirical'`: [`src/emulsim/pipeline/orchestrator.py:39`](../src/emulsim/pipeline/orchestrator.py#L39)

That empirical model:

- computes pore size from a hard-coded concentration/cooling-law correlation
- injects a synthetic log-normal pore distribution using a fixed RNG seed
- returns a uniform `phi_field`
- sets `alpha_final = 0.999` by assumption

**Evidence**

- Default mode selection: [`src/emulsim/pipeline/orchestrator.py:39`](../src/emulsim/pipeline/orchestrator.py#L39)
- Empirical pore law and confinement cap: [`src/emulsim/level2_gelation/solver.py:603`](../src/emulsim/level2_gelation/solver.py#L603)
- Hard-coded pore formula: [`src/emulsim/level2_gelation/solver.py:621`](../src/emulsim/level2_gelation/solver.py#L621)
- Synthetic pore distribution: [`src/emulsim/level2_gelation/solver.py:665`](../src/emulsim/level2_gelation/solver.py#L665)
- Hard-coded nearly complete gelation: [`src/emulsim/level2_gelation/solver.py:674`](../src/emulsim/level2_gelation/solver.py#L674)

**Why this matters**

The project is presented as a multiscale simulation platform, but the default pore result is not the output of the mechanistic Cahn-Hilliard model. It is a calibrated surrogate.

That is acceptable for engineering correlation mode, but it must be stated clearly:

- the pore result is not a first-principles prediction
- its uncertainty is dominated by calibration choice
- most upstream L1 details only influence L2 through the droplet-radius confinement cap, and often not even that

**Assessment**

The default L2 path is best interpreted as an empirical regression layer with synthetic morphology metadata, not as a faithful mechanistic pore-formation simulator.

### 3. High severity: The diffusion-limited Level 3 fallback can silently apply the wrong chemistry

**Finding**

When droplets are large enough and the computed Thiele modulus exceeds 1, `solve_crosslinking()` switches to `_solve_reaction_diffusion()` before chemistry dispatch.

That fallback:

- uses `props.k_xlink_0` and `props.E_a_xlink` rather than the selected crosslinker profile
- assumes amine-group chemistry on chitosan
- uses `c_genipin` as the crosslinker concentration field
- always tags the result as `amine_covalent`

**Evidence**

- Thiele-based fallback trigger uses base properties, not `xl`: [`src/emulsim/level3_crosslinking/solver.py:775`](../src/emulsim/level3_crosslinking/solver.py#L775)
- Reaction-diffusion solver hard-wires amine/genipin-style state variables: [`src/emulsim/level3_crosslinking/solver.py:329`](../src/emulsim/level3_crosslinking/solver.py#L329)
- Returned metadata is always `amine_covalent`: [`src/emulsim/level3_crosslinking/solver.py:790`](../src/emulsim/level3_crosslinking/solver.py#L790)

**Audit demonstration**

With `T_crosslink = 353.15 K` and `R_droplet = 500 um`, three chemically different crosslinkers:

- `dvs`
- `pegda_uv`
- `tpp`

all returned the same `p_final = 0.94176` and all were labeled `amine_covalent`.

That is a concrete implementation defect, not a philosophical modeling limitation.

**Why this matters**

This makes large-droplet predictions for non-amine chemistries physically wrong in a silent way. It is one of the most serious issues in the whole codebase.

### 4. High severity: Mechanistic Level 2 is not operationally viable in the full pipeline

**Finding**

The mechanistic CH solvers are present, but the 2D solver rebuilds a mobility-weighted sparse operator and solves a new sparse linear system on every timestep:

- mobility operator rebuilt every step: [`src/emulsim/level2_gelation/solver.py:447`](../src/emulsim/level2_gelation/solver.py#L447)
- matrix-matrix product and sparse solve every step: [`src/emulsim/level2_gelation/solver.py:468`](../src/emulsim/level2_gelation/solver.py#L468)

In the audit:

- `tests/test_level2_gelation.py` did not finish within 10 minutes
- full pipeline with `ch_1d` or `ch_2d` did not finish within 300 s

**Why this matters**

If the mechanistic branch is too slow to serve as the default analysis engine, then most user-facing results are effectively produced by the empirical branch. That narrows the scientifically defensible claims of the project.

### 5. Medium severity: The 2D CH solver imposes a capped square domain that weakens geometric fidelity

**Finding**

The 2D solver caps the domain size at `1.5e-6 m` regardless of actual droplet size:

- [`src/emulsim/level2_gelation/solver.py:341`](../src/emulsim/level2_gelation/solver.py#L341)

The ternary solver does the same:

- [`src/emulsim/level2_gelation/ternary_solver.py:68`](../src/emulsim/level2_gelation/ternary_solver.py#L68)

**Why this matters**

This helps grid resolution, but it weakens the physical meaning of confinement and geometry, especially for larger droplets. It turns the L2 mechanistic models into local morphology solvers rather than true whole-droplet structure solvers.

**Assessment**

Reasonable as a research simplification, but not consistent with strong claims about full-droplet pore prediction.

### 6. Medium severity: The fast ternary Level 2 path is exploratory, not production-grade

**Finding**

The ternary solver is scientifically attractive because it acknowledges agarose-chitosan-water demixing, but it uses:

- explicit Euler time stepping
- periodic boundary conditions
- a square Cartesian patch rather than a bead geometry

**Evidence**

- Explicit Euler and periodic BCs are stated directly: [`src/emulsim/level2_gelation/ternary_solver.py:9`](../src/emulsim/level2_gelation/ternary_solver.py#L9)
- Explicit update loop: [`src/emulsim/level2_gelation/ternary_solver.py:84`](../src/emulsim/level2_gelation/ternary_solver.py#L84)

**Why this matters**

This mode is useful as a hypothesis generator, but it is not yet a calibrated, geometry-faithful mechanistic production model.

### 7. Medium severity: Uncertainty propagation is partly ad hoc rather than end-to-end physical propagation

**Finding**

The uncertainty module perturbs some real parameters, but it also:

- names a factor `breakage_C1_factor` and then applies it to `breakage_C3`
- perturbs pore uncertainty post hoc after L2 rather than through the L2 model equations

**Evidence**

- Generated perturbations: [`src/emulsim/uncertainty.py:69`](../src/emulsim/uncertainty.py#L69)
- `breakage_C1_factor` applied to `breakage_C3`: [`src/emulsim/uncertainty.py:157`](../src/emulsim/uncertainty.py#L157)
- Post hoc pore perturbation: [`src/emulsim/uncertainty.py:172`](../src/emulsim/uncertainty.py#L172)

**Why this matters**

This is not invalid if described honestly as structural uncertainty injection, but it is not pure Monte Carlo propagation of a mechanistic simulator. The current naming and implementation blur that distinction.

### 8. Medium severity: Optimization uses a numerically discouraged acquisition function and incomplete feasibility checks

**Finding**

The optimization engine uses `qExpectedHypervolumeImprovement`:

- [`src/emulsim/optimization/engine.py:283`](../src/emulsim/optimization/engine.py#L283)

BoTorch emitted a runtime warning during audit that this acquisition has known numerical issues and recommends `qLogExpectedHypervolumeImprovement`.

Separately, the feasibility checks only enforce:

- span
- minimum modulus
- stirred-vessel modal droplet range

**Evidence**

- Constraint set: [`src/emulsim/optimization/objectives.py:78`](../src/emulsim/optimization/objectives.py#L78)

**Why this matters**

Optimization can still run, but it is not yet screening against several scientifically relevant failure modes, such as:

- trust-gate blockers
- chemistry-specific invalidity
- pore accessibility relative to mesh size
- L1 monotonicity/calibration failure

### 9. Medium severity: Mechanical prediction remains phenomenological in the default path

**Finding**

The default DN modulus is:

`G_DN = G_agarose + G_chitosan + eta * sqrt(G_agarose * G_chitosan)`

**Evidence**

- Phenomenological formula: [`src/emulsim/level4_mechanical/solver.py:41`](../src/emulsim/level4_mechanical/solver.py#L41)

**Why this matters**

This is acceptable as an engineering approximation, but it limits the meaning of absolute modulus predictions. In the default run, the final modulus is dominated by the agarose power-law term, while the crosslinked second network provides a much smaller increment.

That makes the pipeline better suited to ranking formulations than to predicting absolute compression behavior.

## Scientific Validity by Level

### Level 1: Emulsification

**Strengths**

- Uses a real PBE structure.
- Includes breakage and coalescence kernels.
- Tracks a full size distribution rather than only moments.

**Weaknesses**

- Current default behavior fails a basic qualitative trend check in audit runs.
- Kernel constants are not strongly system-calibrated.
- The scientific trust of the level is currently lower than the software polish suggests.

**Verdict**

Scientifically useful only after recalibration and verification of RPM dependence.

### Level 2: Gelation and pore formation

**Strengths**

- The codebase contains more than one modeling level.
- There is an empirical mode, a binary CH mode, and a ternary CH mode.

**Weaknesses**

- The default path is empirical rather than mechanistic.
- The binary CH mode is expensive enough to be operationally sidelined.
- The ternary mode is exploratory and not yet production-grade.

**Verdict**

Level 2 is the largest gap between the project’s mechanistic ambitions and its practical predictive basis.

### Level 3: Crosslinking

**Strengths**

- Multi-chemistry branching is implemented.
- The stoichiometric factor-of-two issue in amine consumption appears to be fixed.
- Crosslinking-time saturation behavior is plausible for a crosslinker-limited default formulation.

**Weaknesses**

- Large-droplet diffusion-limited fallback is chemistry-unsafe.
- Hydroxyl and other alternative chemistries are simplified and not consistently protected from misuse.

**Verdict**

Reasonable for small-droplet genipin-like screening; unsafe for broad chemistry claims until the fallback bug is fixed.

### Level 4: Mechanics

**Strengths**

- Clear separation of agarose, crosslinked-network, and combined modulus.
- Hertz and Ogston outputs are useful engineering observables.

**Weaknesses**

- Absolute modulus remains phenomenological in normal use.
- The link from network microstructure to bead mechanics is only partially mechanistic.

**Verdict**

Useful as a ranking layer, not yet a validated absolute mechanics predictor.

## Reliability and Effectiveness

### Reliability

Current reliability is **moderate for code execution** and **limited for scientific prediction**.

The system is reliable at:

- running
- producing structured outputs
- preserving internal data flow
- supporting UI-driven exploration

The system is not yet reliable at:

- guaranteeing physically correct L1 trends in the default setting
- giving chemistry-faithful large-droplet crosslinking results
- supporting routine mechanistic L2 use in reasonable time

### Effectiveness as a simulation tool

It is effective for:

- exploratory scenario comparison
- identifying obviously poor regimes
- exposing which variables the authors think matter
- giving users a structured experimentation interface

It is not yet effective for:

- high-confidence prediction of pore morphology from first principles
- chemistry-specific mechanistic design across all advertised crosslinkers
- optimization-driven experiment planning without expert review

## What Is Scientifically Defensible Today

The following claims are defensible if phrased carefully:

- The software provides a modular process-modeling framework for emulsification, gelation, crosslinking, and mechanics.
- The default workflow is best interpreted as a hybrid empirical-mechanistic engineering model.
- Crosslinking trends for genipin-like small-droplet cases are directionally useful.
- Cooling-rate trends in pore size are available as empirical estimates.

The following claims are **not** yet defensible without stronger evidence:

- The default pore prediction is mechanistic.
- All supported crosslinker chemistries are simulated with comparable fidelity.
- The optimized parameter sets are experimentally reliable by default.
- The simulator provides quantitatively predictive absolute pore and modulus values across regimes.

## Recommended Actions

### Immediate fixes

1. Fix the Level 3 reaction-diffusion fallback so it dispatches by chemistry before selecting the PDE model.
2. Recalibrate or redesign Level 1 so the default pipeline restores monotonic RPM-to-droplet-size behavior.
3. Clearly label default Level 2 output as empirical in CLI/UI/reporting paths.

### Near-term scientific improvements

4. Add calibration datasets and regression tests for `d32(RPM)`, pore size vs cooling rate, and modulus vs crosslinker loading.
5. Split uncertainty into:
   - parameter uncertainty
   - model-form uncertainty
   - post hoc empirical correlation uncertainty
6. Replace `qExpectedHypervolumeImprovement` with the numerically preferred log-EHVI variant.

### Medium-term model improvements

7. Make the ternary L2 model the real mechanistic research path, but with stable semi-implicit time stepping, realistic boundary conditions, and gelation arrest coupling.
8. Add chemistry-specific diffusion-limited crosslinking models or explicitly disable the PDE fallback for unsupported chemistries.
9. Add validation targets that enforce physically correct monotonicity and scale behavior in the actual pipeline, not just in isolated sub-model equations.

## Best Solution Strategy

The best way to cope with the findings is **not** to force every layer into a fully mechanistic model immediately. That would slow the project down and likely make reliability worse before calibration data exist.

The best solution is to deliberately split EmulSim into two scientifically honest operating tiers:

1. **Production Engineering Tier**
2. **Mechanistic Research Tier**

### Production Engineering Tier

This should become the default and should be presented explicitly as a calibrated engineering model.

It should use:

- calibrated L1 emulsification
- empirical or semi-empirical L2 pore model
- chemistry-safe L3 kinetics only in validated regimes
- phenomenological L4 mechanics

It should **not** claim:

- first-principles pore prediction
- chemistry-general mechanistic fidelity
- absolute predictive accuracy outside the calibration domain

The practical goal of this tier is:

- fast runtime
- monotonic physically sensible trends
- conservative trust warnings
- experiment ranking rather than absolute truth

### Mechanistic Research Tier

This should remain available, but it should be framed as slower hypothesis-testing infrastructure.

It should contain:

- ternary agarose-chitosan-water L2 model
- chemistry-specific diffusion-limited L3 solvers only where implemented
- stricter trust gates
- no silent fallback to empirical assumptions

The practical goal of this tier is:

- studying mechanism
- exploring failure of empirical assumptions
- generating new calibration hypotheses

This two-tier structure is already partly present in spirit, but the code and docs do not enforce the distinction strongly enough. That is the core strategic fix.

## Concrete Remediation Plan

### Phase 1: Make the current default path scientifically safe

This phase should happen before any major model expansion.

#### L1 emulsification

- Refit the default kernel constants against a small benchmark dataset of `d32` versus RPM for the intended rotor-stator hardware.
- Add a hard regression test on the **actual PBE output**, not just on `hinze_dmax()`.
- Remove any parameter choice whose only justification is “restores monotonicity” unless it also matches calibration data.

**Success criterion**

- For the calibrated baseline fluid system, `d32(RPM)` is monotonic decreasing over the intended operating window.
- The PBE output passes both trend tests and magnitude tests.

#### L2 pore model

- Keep the empirical model as default, but rename it in user-facing language to something like `empirical_calibrated`.
- Report its output as a calibrated estimate, not as a mechanistic pore simulation.
- Tie uncertainty bands to calibration residuals rather than synthetic post hoc variation alone.

**Success criterion**

- The default user cannot mistake the empirical L2 result for a first-principles pore prediction.
- The model gives stable, fast, calibrated pore estimates in the validated agarose/chitosan window.

#### L3 crosslinking

- Fix chemistry dispatch before the Thiele-based PDE branch.
- If a selected chemistry has no diffusion-limited solver, stop with a warning or fall back only with an explicit “unsupported in diffusion-limited regime” state.
- Do not silently run amine/genipin chemistry for DVS, PEGDA+UV, TPP, or other non-equivalent systems.

**Success criterion**

- Large-droplet runs no longer collapse distinct chemistries into the same amine result.
- Unsupported diffusion-limited chemistries fail loudly rather than fail silently.

### Phase 2: Strengthen scientific trust rather than adding more features

The next priority is not more chemistry options. It is trustworthy validation.

#### Add validation datasets

Minimum validation targets should include:

- `d32` vs RPM
- `d32` vs surfactant concentration
- pore size vs agarose concentration
- pore size vs cooling rate
- modulus vs agarose concentration
- modulus vs crosslinker concentration for at least one chemistry

Each dataset should produce:

- a calibration fit
- a holdout check
- an acceptance range

#### Upgrade trust logic

The trust system should block or penalize:

- non-monotonic L1 behavior in the declared calibration window
- unsupported chemistry/regime combinations
- use of empirical L2 in mechanistic mode without explicit acknowledgement
- optimization proposals outside validation bounds

**Success criterion**

- Trust output becomes regime-aware and validation-aware, not only rule-of-thumb aware.

### Phase 3: Rebuild mechanistic L2 around the ternary path

The binary CH solver is scientifically weaker than the ternary formulation for this material system, and it is too slow in its current form to justify being the main mechanistic branch.

The best long-term solution is:

- keep the empirical L2 as production default
- make the ternary L2 the only serious mechanistic development path
- de-emphasize the binary CH model unless it is needed as a reduced model

#### Required numerical changes

- replace explicit Euler in ternary L2 with a semi-implicit or convex-splitting scheme
- include gelation arrest coupling
- replace periodic square-box assumptions with geometry and boundary conditions that map to microspheres more honestly
- add runtime budgets and coarse/fine resolution modes

**Success criterion**

- mechanistic ternary L2 runs in a bounded research workflow with credible morphology trends and known numerical limits

### Phase 4: Make uncertainty and optimization scientifically honest

#### Uncertainty

Split uncertainty into three categories:

- parameter uncertainty
- calibration uncertainty
- model-form uncertainty

Do not describe post hoc pore perturbation as if it were equivalent to mechanistic propagation.

#### Optimization

- replace `qExpectedHypervolumeImprovement` with the more numerically robust log-EHVI variant
- include trust-gate and validation-domain checks in feasibility
- only optimize over chemistry/regime combinations supported by the validated tier

**Success criterion**

- optimization becomes a constrained search over trusted regions, not a search over the whole software state-space

## Recommended Scientific Positioning After Remediation

After the fixes above, the project should present itself as:

- a **hybrid calibrated process model** for production screening
- plus a **mechanistic research sandbox** for deeper morphology and chemistry studies

It should avoid presenting itself as a universally predictive first-principles simulator until:

- L1 is calibration-stable
- L2 mechanistic mode is both tractable and validated
- L3 diffusion-limited chemistry dispatch is correct
- end-to-end validation against real bead data exists

## Priority Order

If only a limited amount of development time is available, the best order is:

1. Fix Level 3 chemistry dispatch in diffusion-limited mode.
2. Fix and recalibrate Level 1 monotonic RPM behavior.
3. Reframe the default Level 2 path as calibrated empirical output everywhere in the product.
4. Add validation datasets and regression tests.
5. Upgrade uncertainty and optimization to respect validation bounds.
6. Invest in the ternary mechanistic L2 as the long-term research model.

This order gives the largest gain in scientific trust per unit engineering effort.

## Final Verdict

EmulSim is a credible research software prototype with real substance, but it is not yet a scientifically trustworthy quantitative simulator in its current default configuration.

Its biggest current strengths are architecture, modularity, and breadth.
Its biggest current weaknesses are:

- over-reliance on an empirical Level 2 default
- a serious chemistry-dispatch bug in Level 3 diffusion-limited mode
- weak trust in default Level 1 trend behavior
- mechanistic Level 2 paths that are too slow for routine practical use

### Practical recommendation

Use it today as a **screening and hypothesis-generation platform**.
Do not use it today as a **standalone decision engine for wet-lab parameter selection** without manual scientific review, explicit calibration to your own system, and targeted fixes to the issues listed above.
