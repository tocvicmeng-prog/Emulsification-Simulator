# EmulSim Architecture Modification Plan

**Document ID:** EMULSIM-ARCH-MOD-035  
**Date:** 2026-04-17  
**Scope:** Current `Emulsification-Simulator` source, tests, configs, and architecture documents.  
**Purpose:** Provide a practical architecture modification plan to repair, strengthen, and improve EmulSim without discarding its useful working core.

---

## 1. Executive Summary

EmulSim has grown from a four-level fabrication simulator into a broader hydrogel microsphere process platform:

1. **M1 fabrication:** emulsification, gelation/pore formation, crosslinking, and mechanical property prediction.
2. **M2 functionalization:** accessible chemical sites, activation, spacer arms, ligand/protein coupling, quenching, metal charging, washing, and secondary crosslinking.
3. **M3 performance:** chromatography hydrodynamics, breakthrough, gradient elution, multi-component isotherms, detection, and catalytic packed-bed screening.
4. **Cross-cutting systems:** calibration, uncertainty propagation, trust gates, optimization, protocol generation, lifetime projection, and Streamlit UI.

The platform should **not be rewritten from scratch**. The integrated M1 -> M2 -> M3 process chain and chemistry breadth are the core assets. The correct next architecture is an **evidence-centered, calibration-aware simulation architecture**: every output should carry its model basis, valid domain, calibration state, uncertainty source, diagnostics, and downstream trust consequences.

The code already contains important pieces of this architecture: `ModelEvidenceTier`, `ModelManifest`, `RunReport`, L1 dimensionless group utilities, L2 `empirical_uncalibrated` labeling, L3 Thiele diagnostics, L4 network classification, M2 ACS conservation, M3 mass-balance checks, and gradient-aware equilibrium routing. The main problem is incomplete wiring and several immediate defects.

Most important current defects:

- **Critical scientific defect:** L1 PBE default kernel constants still produce a nonphysical RPM trend. A runtime sweep produced `d32=8.66 um` at 5000 RPM and `d32=14.70 um` at 15000 RPM, so stronger agitation made larger droplets. The repository also marks the same RPM and interfacial-tension trend tests as expected failures.
- **Immediate runtime defect:** `python -m emulsim uncertainty` crashes because `__main__.py` imports `emulsim.uncertainty`, while the current implementations are `uncertainty_core.py` and `uncertainty_propagation/`.
- **Trust architecture gap:** trust-aware objective functions exist, but the optimizer still calls the non-trust-aware objective path.
- **Calibration gap:** `CalibrationStore.apply_to_model_params()` exists, but it is not injected into the M1 pipeline or evidence-tier upgrades.
- **UI/documentation drift:** UI metadata and docs still describe some implemented features as missing, and README/config/dataclass defaults are inconsistent.
- **UQ fragmentation:** M1-L4 uncertainty and M2 uncertainty exist as separate systems with different schemas.

Recommended priority: **repair scientific monotonicity and runtime reliability before adding new physics**. The next release should be a hardening release.

---

## 2. Evidence From Inspection

### 2.1 Source Areas Reviewed

Reviewed source and docs include:

- `src/emulsim/datatypes.py`
- `src/emulsim/pipeline/orchestrator.py`
- `src/emulsim/trust.py`
- `src/emulsim/level1_emulsification/`
- `src/emulsim/level2_gelation/`
- `src/emulsim/level3_crosslinking/`
- `src/emulsim/level4_mechanical/`
- `src/emulsim/module2_functionalization/`
- `src/emulsim/module3_performance/`
- `src/emulsim/calibration/`
- `src/emulsim/uncertainty_core.py`
- `src/emulsim/uncertainty_propagation/`
- `src/emulsim/optimization/`
- `src/emulsim/visualization/`
- `configs/`, `data/`, `tests/`, and relevant `docs/`.

Runtime checks performed:

```text
python -m emulsim info
python -m emulsim run --quiet
python -m emulsim sweep --rpm-min 5000 --rpm-max 15000 --rpm-steps 2
python -m emulsim uncertainty --n-samples 1
python -m pytest -q
```

Observed runtime state:

- `python -m emulsim info` succeeds.
- `python -m emulsim run --quiet` succeeds with dataclass defaults and produced:
  - L1 `d32 = 18.08 um`, `span = 0.85`
  - L2 `pore = 179.6 nm`, `porosity = 0.871`
  - L3 `p = 0.040`, `G_chit = 2062 Pa`
  - L4 `G_DN = 70766 Pa`, `E* = 257332 Pa`
- L1 RPM sweep shows nonphysical trend:
  - 5000 RPM -> `d32 = 8.66 um`
  - 15000 RPM -> `d32 = 14.70 um`
- `python -m emulsim uncertainty --n-samples 1` fails with `ModuleNotFoundError: No module named 'emulsim.uncertainty'`.
- Full test run collected 613 tests but exceeded 180 seconds before completion. Early progress showed many passes and two expected failures in `tests/test_l1_trends.py`.

---

## 3. Current Architecture

### 3.1 Current Data Flow

```text
SimulationParameters
        |
        v
PropertyDatabase.update_for_conditions()
        |
        v
MaterialProperties
        |
        v
M1 fabrication pipeline
  L1 PBESolver -> EmulsificationResult
  L2 timing + pore model -> GelationTimingResult, GelationResult
  L3 chemistry-dispatched crosslinking -> CrosslinkingResult
  L4 mechanics/accessibility -> MechanicalResult
        |
        v
FullResult + TrustAssessment + RunReport
        |
        v
M1ExportContract
        |
        v
M2 functionalization -> FunctionalMicrosphere
        |
        v
FunctionalMediaContract
        |
        v
M3 performance -> breakthrough, gradient, catalysis, detection
```

### 3.2 M1 Fabrication

**L1 emulsification** uses a fixed-pivot PBE solver with log-spaced size bins. It supports legacy rotor-stator and stirred-vessel modes, uses Alopaeus and Coulaloglou-Tavlarides style kernels, and includes energy dissipation, temperature-dependent stirred-vessel operation, adaptive PBE extension, size statistics, and dimensionless group utilities.

**L2 gelation and pore formation** has empirical, 1D Cahn-Hilliard, 2D Cahn-Hilliard, and ternary Cahn-Hilliard paths. The default empirical model is now correctly labeled `empirical_uncalibrated`, but it still hardcodes complete gelation and uses a synthetic pore distribution.

**L3 crosslinking** has chemistry-dispatched solvers for amine covalent, hydroxyl covalent, UV dose, ionic instant, and an approximate fallback for Michaelis-Menten-like EDC/NHS. It computes Thiele modulus and can use a reaction-diffusion PDE for amine-reactive large-particle cases.

**L4 mechanics** computes agarose modulus, crosslinked network modulus, DN/IPN phenomenological modulus, ionic and triple-network variants, Hashin-Shtrikman reference bounds, optional affine Flory-Rehner IPN mode, Hertz contact, and Ogston partitioning. It classifies the network using L3 metadata.

### 3.3 M2 Functionalization

M2 is built around stable module boundaries:

- `M1ExportContract`
- `AccessibleSurfaceModel`
- `ACSSiteType`
- `ACSProfile`
- `ModificationStep`
- `ModificationResult`
- `ModificationOrchestrator`
- `FunctionalMicrosphere`
- `FunctionalMediaContract`

This is the right direction. M2 should continue consuming `M1ExportContract` rather than internal M1 solver objects.

### 3.4 M3 Performance

M3 includes:

- `ColumnGeometry` with Kozeny-Carman pressure and compression screening.
- Lumped Rate Model breakthrough with finite-volume axial discretization.
- Langmuir, competitive Langmuir, SMA, HIC, IMAC, Protein A, and competitive affinity isotherms.
- `EquilibriumAdapter` for gradient-sensitive process-state variables.
- UV, fluorescence, conductivity, and MS detection models.
- Multi-component gradient elution with peak, yield, purity, resolution, pressure, and mass-balance outputs.
- Catalytic packed-bed screening.

M3 is useful for downstream process design, but most default isotherm parameters remain illustrative unless calibrated.

### 3.5 Cross-Cutting Systems

Existing cross-cutting components include:

- `TrustAssessment`
- `ModelEvidenceTier`, `ModelManifest`, and `RunReport`
- `CalibrationEntry` and `CalibrationStore`
- `UncertaintyPropagator` for M1-L4
- `M1UncertaintyContract` and `run_with_uncertainty()` for M2
- BoTorch optimization
- Streamlit UI modules and validators
- Protocol generation
- Lifetime projection

These are good foundations. They now need one common contract and enforcement path.

---

## 4. Strengths To Preserve

1. **Integrated process chain:** M1 -> M2 -> M3 captures fabrication, functionalization, and downstream performance.
2. **Chemistry-aware routing:** L3 and M2 do not treat all reagents as interchangeable.
3. **Stable data contracts:** `M1ExportContract` and `FunctionalMediaContract` decouple modules.
4. **Trust-gate culture:** the platform already warns about invalid assumptions, unsupported regimes, poor mass balance, pore limits, and low crosslinking.
5. **Calibration-first direction:** schema, stores, UI panels, and confidence fields already exist.
6. **Mechanistic upgrade paths:** L2 has empirical/CH/ternary paths, L3 has ODE/PDE paths, and M3 has simple/gradient-sensitive isotherm paths.
7. **Wet-lab realism:** profiles and comments recognize pH, toxicity, hydrolysis, reversibility, activity retention, sterics, and washing.
8. **Broad tests:** tests cover fabrication physics, M2 workflows, M3 chromatography/catalysis, UI contracts, calibration, protocols, and integration.

---

## 5. Critical Findings

### F1. L1 Droplet Physics Has A Nonphysical Default Trend

**Severity:** Critical  
**Subsystem:** `level1_emulsification`

The default L1 PBE predicts larger droplets at higher RPM in the checked sweep. This is the largest scientific risk because downstream pore, crosslinking, mechanics, functionalization, and chromatography predictions inherit bead size.

Repair:

- Build an L1 calibration dataset for DSD vs RPM, surfactant, viscosity, geometry, and phi_d.
- Fit geometry-specific breakage/coalescence and dissipation constants.
- Replace expected-fail monotonic tests with passing regression tests.
- Add a validity-domain downgrade when sign constraints fail.
- Add Hinze/viscous-Hinze sanity diagnostics as checks, not replacements for PBE.

### F2. CLI Uncertainty Command Is Broken

**Severity:** High  
**Subsystem:** CLI / uncertainty

`__main__.py` imports `from .uncertainty import UncertaintyPropagator`, but no `uncertainty.py` module exists.

Repair:

- Import `UncertaintyPropagator` from `uncertainty_core.py`, or add a compatibility module.
- Split CLI uncertainty modes if both M1-L4 and M2 uncertainty remain:
  - `emulsim uncertainty m1`
  - `emulsim uncertainty m2`
- Add a CLI smoke test.

### F3. Trust-Aware Optimization Is Not Used

**Severity:** High  
**Subsystem:** `optimization`

`compute_objectives_trust_aware()` exists, but `OptimizationEngine._evaluate()` calls `compute_objectives(result)`.

Repair:

- Make trust-aware objectives default.
- Add `trust_mode = "strict" | "penalize" | "off"`.
- Exclude `UNSUPPORTED` or blocker-level candidates from Pareto recommendations.
- Label Pareto candidates with weakest evidence tier and trust warnings.

### F4. Calibration Service Is Not A First-Class Parameter Provider

**Severity:** High  
**Subsystem:** calibration / all solvers

`CalibrationStore.apply_to_model_params()` exists, but the orchestrator does not accept or apply calibrated model parameters before solving.

Repair:

- Introduce a `ParameterProvider` or `RunContext`.
- Resolve each model parameter from user override, calibration store, config, property database, literature, or default.
- Record source, units, uncertainty, and valid domain.
- Upgrade evidence tier only when calibration applies inside its valid domain.

### F5. Evidence Metadata Is Present But Not Uniformly Enforced

**Severity:** High  
**Subsystem:** datatypes / solvers / UI / M2 / M3

M1 result classes have optional `model_manifest`, but M2 and M3 use string confidence fields or no manifest. `valid_domain` is often empty, and L1 dimensionless groups are not attached to result diagnostics.

Repair:

- Require every public result object to include a manifest.
- Populate `valid_domain` and `diagnostics` in solver entry points.
- Map legacy strings such as `ranking_only`, `semi_quantitative`, `requires_user_calibration`, and `mapped_estimated` to `ModelEvidenceTier`.
- Treat missing manifest as `SEMI_QUANTITATIVE` at best.

### F6. UI Metadata And Validators Drifted From Backend Reality

**Severity:** Medium-High  
**Subsystem:** `visualization`

UI metadata still says gradient competitive Langmuir does not affect binding, despite gradient-aware routing being implemented. `validate_m2_inputs()` lists five supported step types while the backend supports nine.

Repair:

- Generate UI evidence labels from backend manifests.
- Update M2 UI validation to mirror backend enums.
- Remove stale warnings that contradict implemented gradient-aware M3 behavior.
- Add UI contract tests comparing UI-supported values to backend-supported values.

### F7. UQ Is Split Into Two Systems

**Severity:** Medium-High  
**Subsystem:** uncertainty

`uncertainty_core.py` propagates M1-L4 material/model perturbations. `uncertainty_propagation/` propagates M1ExportContract uncertainty through M2. They do not share schema, source taxonomy, or calibration posterior handling.

Repair:

- Define shared `UncertaintySpec` and `UncertaintyResult`.
- Tag uncertainty sources as measurement, material property, calibration posterior, model form, numerical, or scale-up.
- Add optional correlation matrices.
- Sample calibrated model parameters from posterior distributions when available.

### F8. L2 Default Pore Model Is Honest But Still Weak For Quantitative Design

**Severity:** Medium-High  
**Subsystem:** L2 gelation

The empirical L2 path is now honestly labeled `empirical_uncalibrated`, but it still hardcodes `alpha_final = 0.999`, uses a synthetic log-normal pore distribution, and runs one representative bead size.

Repair:

- Replace hardcoded gelation completion with `solve_gelation_timing()` output.
- Add calibrated empirical coefficients for pore prefactor, concentration exponent, cooling exponent, chitosan factor, and confinement cap.
- Run L2 over L1 size quantiles or DSD-weighted bins.
- Use CH/ternary solvers to train a fast surrogate with explicit validity domain.

### F9. L3 Chemistry Fallbacks Need Evidence Consequences

**Severity:** Medium-High  
**Subsystem:** L3 crosslinking

EDC/NHS logs an approximate fallback to second-order amine kinetics. Trust gates warn that EDC/NHS requires carboxyl groups, but the solver can still produce a normal semi-quantitative result.

Repair:

- Add `ChemistryApplicability` checks before solving.
- Require EDC/NHS carboxyl-bearing chemistry or introduced COOH sites.
- Downgrade approximate fallbacks to `QUALITATIVE_TREND` or `UNSUPPORTED`.
- Include pH, ionic strength, temperature compatibility, diffusion regime, and optical penetration checks.

### F10. M3 Outputs Need Evidence And Calibration Inheritance

**Severity:** Medium  
**Subsystem:** M3 performance

`BreakthroughResult` and `GradientElutionResult` do not carry manifests or evidence tiers. M3 consumes FMC values but does not consistently inherit whether `q_max`, activity retention, and isotherm constants are estimated or calibrated.

Repair:

- Add `model_manifest`, `evidence_tier`, and `calibration_refs` to M3 result objects.
- If `estimated_q_max` is uncalibrated, M3 DBC/yield/purity cannot be stronger than semi-quantitative.
- Treat mass balance >5 percent as a blocker for decision-grade M3 outputs.
- Make `select_isotherm_from_fmc()` return both isotherm and manifest.

### F11. Unit Handling Is Mostly By Comments

**Severity:** Medium  
**Subsystem:** all modules

Dataclasses use comments such as `[kg/m3]`, `[mol/m3]`, and `[m]`, but interface-level unit enforcement is limited. UI validators use user-facing units that differ from internal SI.

Repair:

- Keep internal SI.
- Add targeted unit assertions at config load, UI conversion, M1ExportContract, FunctionalMediaContract, ColumnGeometry, and CalibrationEntry ingest.
- Add `UnitSpec` metadata for user-facing and file-facing fields.
- Do not add a full `pint` dependency unless a clear need appears.

### F12. Documentation And Defaults Are Diverging

**Severity:** Medium  
**Subsystem:** docs/config/README

README examples, dataclass defaults, and `configs/default.toml` do not fully agree. The inspected config run exceeded a 120 second timeout in the current environment.

Repair:

- Define one canonical default for CLI, README, quickstart, and tests.
- Add `configs/fast_smoke.toml` that finishes in under 30-60 seconds.
- Keep the realistic default config but document the expected runtime.
- Include expected baseline outputs in docs and smoke tests.

---

## 6. Target Architecture

### 6.1 Guiding Principle

EmulSim should become an evidence-carrying simulator:

> A numerical value is not a platform output unless it carries units, model basis, valid domain, calibration source, uncertainty source, diagnostics, and trust consequence.

### 6.2 Proposed Core Objects

#### RunContext

Carries cross-cutting inputs into a run:

```python
@dataclass
class RunContext:
    run_id: str
    model_mode: ModelMode
    calibration_store: CalibrationStore | None
    uncertainty_spec: UncertaintySpec | None
    parameter_provider: ParameterProvider
    random_seed: int
    output_policy: OutputPolicy
```

#### ParameterProvider

Resolves parameters with provenance:

```python
@dataclass
class ResolvedParameter:
    name: str
    value: float
    unit: str
    source_type: str       # user, calibration, config, database, literature, default
    source_ref: str
    uncertainty: float
    valid_domain: dict
```

#### Extended Model Manifest

Extend current `ModelManifest` rather than replacing it:

```python
@dataclass
class ModelRunManifest:
    model_name: str
    version: str
    evidence_tier: ModelEvidenceTier
    valid_domain: dict
    domain_evaluation: dict
    calibration_refs: list[str]
    assumptions: list[str]
    diagnostics: dict
    unit_checks: dict
    uncertainty_sources: list[str]
```

#### RunDossier

`RunDossier` should aggregate:

- input parameters,
- resolved parameters,
- model graph,
- result summaries,
- diagnostics,
- uncertainty summaries,
- trust warnings/blockers,
- calibration recommendations,
- environment and code version metadata,
- JSON/HDF5 export references.

### 6.3 Target Data Flow

```text
Config/UI/API input
        |
        v
Validated SimulationParameters + RunContext
        |
        v
ParameterProvider resolves model parameters with provenance
        |
        v
Solvers emit Result + ModelRunManifest + diagnostics
        |
        v
TrustEngine evaluates domains, conservation, evidence, and numerical quality
        |
        v
RunDossier
        |
        v
UI, CLI, optimizer, protocols, JSON/HDF5 exports
```

---

## 7. Modification Roadmap

### Phase 0: Immediate Repair And Release Stabilization

**Goal:** Make the current platform reliably runnable and prevent known-bad outputs from being recommended.

**Duration:** 2-4 engineering days.

Deliverables:

1. Fix CLI uncertainty import.
2. Add a fast smoke config and CLI smoke tests:
   - `python -m emulsim info`
   - `python -m emulsim run configs/fast_smoke.toml --quiet`
   - `python -m emulsim sweep ...`
   - `python -m emulsim uncertainty --n-samples 2`
3. Mark slow tests with `pytest.mark.slow`.
4. Switch optimization to trust-aware objectives by default.
5. Update README/quickstart/defaults.
6. Update stale UI metadata for gradient-aware M3 and M2 step types.
7. Add a release blocker stating that uncalibrated L1 PBE trends are not quantitative.

Acceptance criteria:

- CLI uncertainty command no longer crashes.
- Fast smoke suite completes under a defined time budget.
- Optimizer does not rank unsupported or blocker-level designs as Pareto candidates.
- Docs and default config agree on primary defaults and runtime expectations.

### Phase 1: Complete Evidence And Validity-Domain Wiring

**Goal:** Every public output carries evidence, valid domain, and diagnostics.

**Duration:** 1 week.

Deliverables:

1. Extend `ModelManifest` with version, domain evaluation, calibration refs, uncertainty sources, and unit checks.
2. Add manifests to M2 and M3 public result classes.
3. Add L1 dimensionless groups to manifest diagnostics.
4. Add L2 empirical calibration status and Avrami timing diagnostics.
5. Add L3 required-chemistry diagnostics and fallback flags.
6. Add M3 mass-balance quality, hydrodynamic validity, and isotherm calibration state.
7. Make UI evidence badges read from manifests, not parallel hardcoded tables.

Acceptance criteria:

- A complete M1 -> M2 -> M3 run can export a model graph.
- No numeric UI output lacks evidence tier.
- Missing manifest is a test failure for public result classes.

### Phase 2: L1 PBE Scientific Repair

**Goal:** Restore physically correct droplet-size trends and calibrate kernel parameters by equipment mode.

**Duration:** 1-2 weeks plus data collection.

Deliverables:

1. Add `data/validation/l1_dsd/` schema for RPM, geometry, viscosity, IFT, phi_d, and DSD.
2. Add fitting tools for breakage/coalescence and dissipation constants.
3. Add domain diagnostics: Re, We, Ca, viscosity ratio, Kolmogorov length, d32/eta_K, and residence/circulation time.
4. Add monotonic tests:
   - increasing RPM -> non-increasing d32,
   - increasing sigma -> non-decreasing d32,
   - increasing dispersed viscosity -> non-decreasing d32,
   - increasing surfactant reduces d32 until coalescence suppression saturates.
5. Retire current `xfail` trend tests after calibration.
6. Downgrade evidence when uncalibrated parameter regions violate signs.

Scientific note:

Do not force monotonicity by clipping the output. Diagnose whether the issue comes from breakage scale, coalescence scale, dissipation ratio, premix initialization, binning, convergence time, or average/max dissipation usage.

### Phase 3: L2 Pore Model Upgrade

**Goal:** Make pore predictions defensible for design rather than illustrative.

**Duration:** 1-2 weeks.

Deliverables:

1. Replace empirical `alpha_final=0.999` with actual gelation timing output.
2. Add calibrated empirical coefficients to `CalibrationStore`.
3. Add pore calibration dataset ingestion for cryo-SEM, tracer exclusion, SEC probe, or porosimetry data.
4. Run L2 over d10, d50, d90, or DSD-weighted bins.
5. Produce batch-level pore distribution and uncertainty.
6. Use CH/ternary solvers to generate mechanistic training data for a fast surrogate.

Acceptance criteria:

- Empirical L2 remains `SEMI_QUANTITATIVE` until local calibration exists.
- Calibrated L2 carries calibration reference and valid concentration/cooling domain.
- Batch pore distribution can be exported to M2/M3.

### Phase 4: L3/L4 Chemistry And Mechanics Hardening

**Goal:** Make chemistry applicability and mechanical predictions match wet-lab reality.

**Duration:** 1-2 weeks.

Deliverables:

1. Add `ChemistryApplicability` preflight checks.
2. Require carboxyl groups or introduced COOH functionality for EDC/NHS.
3. For hydroxyl crosslinkers, model or bound chitosan amine side reaction.
4. For UV crosslinking, include photoinitiator and optical attenuation before quantitative claims.
5. For ionic TPP, include salt/pH reversibility and buffer stability warnings.
6. Refine L4 network classification for covalent IPN, semi-IPN, PEG network, ionic reinforcement, and mixed hydroxyl/amine networks.
7. Add swelling and wash-equilibration diagnostics before M2/M3 handoff.

Acceptance criteria:

- Unsupported chemistry produces `UNSUPPORTED`.
- Approximate fallbacks produce `QUALITATIVE_TREND` and are excluded from trust-strict optimization.
- Mechanical result includes network class, model used, calibration status, and swelling assumptions.

### Phase 5: M2/M3 Process Simulation Hardening

**Goal:** Make downstream functionalization and performance predictions inherit upstream evidence and calibration state.

**Duration:** 1-2 weeks.

Deliverables:

1. Add manifests to:
   - `ModificationResult`,
   - `FunctionalMicrosphere`,
   - `FunctionalMediaContract`,
   - `BreakthroughResult`,
   - `GradientElutionResult`,
   - catalytic results.
2. Add calibration inheritance:
   - q_max from calibration -> `CALIBRATED_LOCAL` or `VALIDATED_QUANTITATIVE`,
   - q_max from ligand-density estimate -> `SEMI_QUANTITATIVE`,
   - activity-retention default -> warning,
   - uncalibrated HIC/lectin/affinity constants -> `requires_user_calibration`.
3. Return `(isotherm, manifest)` from isotherm selection.
4. Treat M3 mass balance as an output gate:
   - <=2 percent acceptable,
   - 2-5 percent caution,
   - >5 percent blocker.
5. Standardize `ProcessState` everywhere.
6. Add packed-bed validity diagnostics:
   - Re_p,
   - Peclet number,
   - pressure drop,
   - bed compression,
   - flow relative to crushing limit,
   - mass-transfer Biot/Thiele where relevant.

Acceptance criteria:

- M3 DBC/yield/purity outputs cannot be labeled stronger than their q_max/isotherm calibration.
- UI blocks M3 decision-grade reporting when mass balance is unreliable.
- M2/M3 tests verify evidence inheritance.

### Phase 6: Unified Calibration And Uncertainty

**Goal:** Replace fragmented UQ and calibration flows with one coherent architecture.

**Duration:** 2-3 weeks.

Deliverables:

1. Define `CalibrationFitResult` with fitted parameters, covariance/posterior samples, valid domain, residual diagnostics, and dataset reference.
2. Define a minimal `AssayRecord` schema for measurement type, sample identity, process conditions, units, replicate values, instrument/method, and notebook reference.
3. Define `UncertaintySpec` shared by M1, M2, and M3.
4. Let Monte Carlo sample:
   - measured input distributions,
   - calibrated posterior distributions,
   - model discrepancy priors,
   - numerical tolerance flags.
5. Propagate uncertainty through M1 -> M2 -> M3 for q_max, DBC, pressure, and selected process metrics.
6. Export uncertainty sources and confidence intervals to `RunDossier`.

Acceptance criteria:

- UQ output distinguishes measured uncertainty from assumed screening uncertainty.
- Calibrated parameters can contribute posterior uncertainty.
- M2/M3 uncertainty does not silently ignore upstream M1 uncertainty.

### Phase 7: Productization, UI, And Documentation

**Goal:** Make the platform reproducible, understandable, and hard to misuse.

**Duration:** 1-2 weeks.

Deliverables:

1. Add `RunDossier` JSON export in CLI and UI.
2. Add HDF5 or NPZ export for large arrays and histories.
3. Add canonical configs:
   - `configs/fast_smoke.toml`,
   - `configs/baseline_research.toml`,
   - `configs/stirred_vessel.toml`.
4. Split CLI uncertainty modes.
5. Generate protocol recommendations from trust gaps:
   - measure IFT,
   - run DSD RPM sweep,
   - measure q_max by static binding,
   - run DBC10 breakthrough.
6. Document dependency extras:
   - core,
   - ui,
   - optimization,
   - dev.
7. Add CI-style test commands:
   - fast unit,
   - slow mechanistic,
   - integration,
   - optional optimization.

Acceptance criteria:

- A new user can run a short example in under 60 seconds.
- A research user can run the full pipeline with a known expected time budget.
- UI, CLI, README, and docs describe the same evidence tiers and defaults.

---

## 8. Implementation Priority Table

| Priority | Work Item | Primary Files | Risk Reduced |
|---|---|---|---|
| P0 | Fix CLI uncertainty import | `src/emulsim/__main__.py` | Broken command |
| P0 | Add fast smoke config/tests | `configs/`, `tests/` | Release reliability |
| P0 | Make optimizer trust-aware | `optimization/engine.py`, `optimization/objectives.py` | Unsafe recommendations |
| P0 | Update stale UI metadata | `visualization/ui_model_metadata.py`, `visualization/ui_validators.py` | User confusion |
| P1 | Attach L1 dimensionless diagnostics to manifest | `level1_emulsification/solver.py`, `validation.py` | Domain misuse |
| P1 | Calibrate/repair L1 kernels | `level1_emulsification/`, `calibration/`, `data/validation/` | Scientific correctness |
| P1 | Generalize calibration injection | `pipeline/orchestrator.py`, `calibration/` | Evidence accuracy |
| P2 | Add manifests to M2/M3 outputs | `module2_functionalization/`, `module3_performance/` | Downstream trust |
| P2 | Unify UQ contracts | `uncertainty_core.py`, `uncertainty_propagation/` | Misleading intervals |
| P2 | L2 batch pore distribution and calibrated empirical model | `level2_gelation/` | Pore prediction quality |
| P3 | Chemistry applicability engine | `level3_crosslinking/`, `reagent_library.py`, `module2_functionalization/` | Wet-lab realism |
| P3 | RunDossier export | `datatypes.py`, `pipeline/orchestrator.py`, UI/CLI | Reproducibility |

---

## 9. Wet-Lab Calibration And Validation Plan

Architecture changes should be paired with a minimal experimental campaign. Without calibration, many outputs must remain semi-quantitative.

### Study A: Emulsification And DSD

Purpose: calibrate L1 PBE kernels and dissipation factors.

Measurements:

- RPM sweep at 5-7 speeds.
- Surfactant concentration sweep.
- Oil temperature sweep.
- DSD by laser diffraction or microscopy/image analysis.
- Interfacial tension by pendant drop tensiometry.
- Dispersed phase viscosity at relevant shear rates.

Parameters informed:

- `breakage_C1`, `breakage_C2`, `breakage_C3`.
- `coalescence_C4`, `coalescence_C5`.
- `dissipation_ratio`.
- IFT model parameters.
- Shear-thinning parameters.

### Study B: Gelation And Pore Formation

Purpose: calibrate L2 pore model.

Measurements:

- Agarose/chitosan concentration matrix.
- Cooling rate sweep.
- Cryo-SEM or confocal pore morphology.
- Porosity by tracer exclusion or wet/dry swelling.
- Gelation onset by rheology.

Parameters informed:

- empirical pore prefactor/exponents,
- chitosan factor,
- confinement cap,
- `T_gel`,
- Avrami constants,
- CH `kappa` and mobility scale.

### Study C: Crosslinking Kinetics

Purpose: calibrate L3 kinetics and bridge efficiency.

Measurements:

- Ninhydrin/TNBS residual amine assay vs time.
- Crosslinker concentration sweep.
- Temperature sweep.
- Rheology/compression vs crosslinking time.
- Wash leachate assay for residual crosslinker.

Parameters informed:

- `k_xlink_0`,
- `E_a_xlink`,
- `f_bridge`,
- stoichiometric ceiling,
- Thiele/reaction-diffusion validity.

### Study D: Mechanics And Swelling

Purpose: calibrate L4 mechanical models.

Measurements:

- Small-amplitude rheology for bulk gels.
- Bead compression/Hertz indentation.
- Packed-bed compression versus pressure.
- Swelling ratio in process buffers.

Parameters informed:

- agarose modulus prefactor/exponent,
- IPN coupling eta,
- Flory-Rehner chi proxies,
- bed compression safety limits.

### Study E: Functionalization And Capacity

Purpose: calibrate M2/M3 binding capacity and activity retention.

Measurements:

- Ligand density by colorimetric or elemental assay.
- Activity retention by binding assay.
- Static q_max.
- DBC5/DBC10 breakthrough.
- Gradient elution recovery and purity.
- Mass balance closure.

Parameters informed:

- ACS accessibility,
- reaction yields,
- q_max mapping,
- isotherm constants,
- M3 evidence tiers.

---

## 10. Architectural Rules To Adopt

1. Every public number carries evidence.
2. Every model declares its valid domain and evaluates the current run against it.
3. Every empirical model declares calibration provenance.
4. Every mechanistic model returns diagnostics.
5. Every fallback either preserves chemistry or downgrades/blocks evidence.
6. Every optimizer candidate passes trust gates before recommendation.
7. Every downstream result inherits upstream trust and uncertainty.
8. Every unit-operation model checks conservation or mass balance where applicable.
9. Every user-facing input has explicit units and conversion tests.
10. Every docs example is backed by a smoke test.

---

## 11. Proposed Next-Release Acceptance Criteria

The next release should be accepted only if:

- `python -m emulsim uncertainty --n-samples 2` runs without import error.
- Fast tests complete under a defined time budget.
- L1 trend tests for RPM and sigma are either passing after calibration or hard-blocked from quantitative recommendation.
- Optimizer uses trust-aware objectives by default.
- UI no longer advertises stale warnings for implemented gradient-sensitive binding.
- M2 UI validators support all backend step types or explicitly label unavailable UI paths.
- Every M1 result has a populated `ModelManifest`.
- M2/M3 outputs include evidence tier or are explicitly marked legacy/internal.
- `RunReport` or `RunDossier` export includes model graph, trust warnings, diagnostics, and weakest evidence tier.
- README, quickstart, and config defaults are mutually consistent.

---

## 12. Final Recommendation

The strongest architecture path is:

1. Repair reliability first: CLI uncertainty, fast test path, stale UI metadata, trust-aware optimizer.
2. Repair L1 physics before deeper optimization: the nonphysical RPM trend is the largest scientific risk.
3. Promote evidence metadata to a platform contract: UI, optimizer, protocols, and M3 should not consume unqualified numerical outputs.
4. Make calibration the parameter source of truth: model parameters should know whether they are measured, fitted, literature, estimated, or default.
5. Unify UQ across M1/M2/M3: uncertainty must travel downstream with source labels.
6. Then add advanced physics: L2 surrogate, richer L3 chemistry, and M3 process models should follow the evidence and calibration architecture.

EmulSim is already a valuable research platform. Its next step is not broader feature count; it is scientific defensibility, reproducibility, and trust-aware decision support.

