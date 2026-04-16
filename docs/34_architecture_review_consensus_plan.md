# Joint Architecture Review: Third-Party Redesign Report vs Current EmulSim

**Report:** EMULSIM-ARCH-CONSENSUS-001 Rev 1.0
**Date:** 2026-04-17
**Reviewers:** Scientific Advisor, Architect, Dev-Orchestrator (joint assessment)
**Subject:** EmulSim_Architecture_Redesign_Modification_Plan_2026-04-16.md

---

## 1. Executive Summary

Three perspectives (Scientific Advisor, Architect, Dev-Orchestrator) independently reviewed the third-party redesign report against the current EmulSim v6.0 codebase. We identify **14 consensus points** where the report's observations are scientifically valid, logically consistent, and technically feasible. We reject **4 proposals** as premature or architecturally inappropriate for the current project stage. We then synthesize these into a **6-sprint modification plan** that perfects, improves, strengthens, and repairs EmulSim while preserving its existing strengths.

**Key finding:** The third-party report is directionally correct and well-reasoned. Its core thesis -- that EmulSim needs stronger evidence contracts, not a rewrite -- aligns with the existing codebase's trajectory. However, several proposals (ProcessDossier, full unit validation, multi-tier UQ decomposition) represent scope that exceeds the current project's needs and should be deferred to v7.0+.

---

## 2. Method

### What was compared:
- **Third-party report**: 10 gaps identified, 10-phase modification plan, 10 architectural rules
- **Current codebase**: Every source file in src/emulsim/ (read in full during this session)
- **Assessment criteria**: Each proposal evaluated on three axes:
  1. **Scientifically valid** -- Does the proposal correctly identify a real scientific weakness?
  2. **Logically consistent** -- Is the proposed fix internally coherent and compatible with existing architecture?
  3. **Technically feasible** -- Can it be implemented without destabilizing the working system?

---

## 3. Current Architecture Factual Audit

Before comparing, here is what actually exists (verified by reading source):

| Subsystem | Current State | Evidence |
|---|---|---|
| **Trust gates** | 3-tier (TRUSTWORTHY/CAUTION/UNRELIABLE), 12+ checks including d32 bounds, convergence, pore range, gelation, stoichiometry, viscosity ratio, chemistry-specific warnings | `trust.py` lines 38-175 |
| **Model tier metadata** | L2 has `model_tier` ("empirical_calibrated"/"mechanistic"), L4 has `model_used` (5 variants), M1ExportContract carries both, L3 has `NetworkTypeMetadata` | `datatypes.py` lines 844, 887, 929, 943 |
| **Thiele modulus** | Computed, logged, warns at >0.3, triggers reaction-diffusion PDE at >1.0 for supported chemistries, raises NotImplementedError for unsupported | `level3_crosslinking/solver.py` lines 45-67, 205-227, 754-789 |
| **Chemistry dispatch** | Explicit per-chemistry routing (amine_covalent, hydroxyl_covalent, ionic_reversible, independent_network, uv_dose, michaelis_menten), no silent substitution for unsupported | `solver.py` lines 1-16, 840-940 |
| **Calibration** | CalibrationEntry + CalibrationStore, applies to FunctionalMediaContract, loads from JSON, logs overrides | `calibration/calibration_store.py` |
| **UQ** | Monte Carlo over material property uncertainties, 3 conceptual categories, produces 90% CIs for d32/pore/G_DN | `uncertainty_core.py` lines 1-60 |
| **L1 validation** | Mode-specific bounds, mass conservation check, steady-state check | `level1_emulsification/validation.py` |
| **Optimization** | BoTorch multi-objective (3 objectives), 7D parameter vector, Pareto fronts | `optimization/` |
| **Pipeline** | L1->L2->L3->L4 chain, exports M1ExportContract to M2 | `pipeline/orchestrator.py` |

---

## 4. Point-by-Point Assessment of Third-Party Gaps

### Gap 1: No Single Process Dossier
| Axis | Assessment |
|---|---|
| Scientifically valid | **Partially.** A richer experiment+evidence container is genuinely useful for reproducibility. |
| Logically consistent | **Yes** for the concept; **No** for immediate implementation -- the current SimulationParameters + FullResult + M1ExportContract already carry substantial provenance. |
| Technically feasible | **Deferred.** Creating ProcessDossier is a v7.0 feature. The current typed contracts are functional. |

**Verdict: AGREE in principle, DEFER to v7.0.** Add `RunReport` metadata to the orchestrator output (lightweight) as an interim step.

---

### Gap 2: Model Evidence Not Strong Enough
| Axis | Assessment |
|---|---|
| Scientifically valid | **Yes.** The model_tier/model_used fields exist but only on L2 and L4. L1 and L3 lack explicit evidence tier on their results. |
| Logically consistent | **Yes.** A uniform `ModelEvidenceTier` enum across all result objects would prevent confusion. |
| Technically feasible | **Yes.** Adding an enum + field to existing dataclasses is backward-compatible. |

**Verdict: STRONG AGREE. Priority 1.**

---

### Gap 3: L1 Droplet Physics Needs Trend Repair
| Axis | Assessment |
|---|---|
| Scientifically valid | **Yes.** The audit report documented a nonphysical RPM trend. Current validation checks bounds but does NOT enforce monotonic physical trends via regression tests. No dimensionless group report. |
| Logically consistent | **Yes.** Trend tests + dimensionless groups are complementary to existing bound checks. |
| Technically feasible | **Yes.** Add pytest-based monotonic trend tests and a dimensionless_groups() function. |

**Verdict: STRONG AGREE. Priority 2.**

---

### Gap 4: L2 Pore Model Bottleneck
| Axis | Assessment |
|---|---|
| Scientifically valid | **Yes.** The empirical pore model (`model_tier="empirical_calibrated"`) carries a misleading label -- it's empirical but NOT calibrated against the user's specific material lot. |
| Logically consistent | **Yes.** The 3-tier split (empirical/CH/surrogate) is the right architecture. The label should be `empirical_uncalibrated` by default. |
| Technically feasible | **Partially.** Renaming the label is trivial. Building a surrogate model is v7.0 scope. |

**Verdict: AGREE on labeling fix and CH test tiering. DEFER surrogate model.**

---

### Gap 5: Calibration Store Too Narrow
| Axis | Assessment |
|---|---|
| Scientifically valid | **Yes.** Current CalibrationStore only applies to FunctionalMediaContract. L1 kernel constants, L2 pore coefficients, and L3 kinetics cannot be calibrated through the store. |
| Logically consistent | **Yes.** A cross-module calibration service is the natural evolution. |
| Technically feasible | **Yes for v6.1.** Extend CalibrationEntry to target any model parameter. Keep the same JSON interface. |

**Verdict: AGREE. Priority 3.**

---

### Gap 6: Units Are Mostly Comments
| Axis | Assessment |
|---|---|
| Scientifically valid | **Yes** in principle. Unit errors are a real risk in mixed-unit systems. |
| Logically consistent | **Partially.** Full dimensional analysis (pint/unyt library) would add significant complexity and runtime overhead for marginal benefit in a simulation where all internal calculations use SI. |
| Technically feasible | **DEFER.** A full unit-validation layer is disproportionate for the current project. Instead, add explicit unit assertions at key interfaces (input validation, export contracts). |

**Verdict: PARTIALLY AGREE. Add targeted assertions, not a full unit framework.**

---

### Gap 7: UQ Still Approximate
| Axis | Assessment |
|---|---|
| Scientifically valid | **Yes.** The current MC samples independent priors. Correlated posteriors, model discrepancy, and numerical error are not decomposed. |
| Logically consistent | **Yes** for a production platform. **Premature** for current scope -- the existing MC UQ provides useful 90% CIs. |
| Technically feasible | **DEFER.** Full posterior UQ (measurement, material, calibration, model-form, numerical, scale-up as separate components) is a v7.0 research feature. |

**Verdict: PARTIALLY AGREE. Add UQ source tagging to existing MC, defer full decomposition.**

---

### Gap 8: M2/M3 Weak Calibration Integration
| Axis | Assessment |
|---|---|
| Scientifically valid | **Yes.** qmax, isotherm constants, and activity retention are default estimates, not calibrated values. The `confidence_tier` field exists but the UI doesn't block based on it. |
| Logically consistent | **Yes.** M3 outputs should inherit calibration status from M2 inputs. |
| Technically feasible | **Yes.** Add calibration provenance checks in M3 entry and label outputs accordingly. |

**Verdict: AGREE. Priority 4.**

---

### Gap 9: Optimization Can Outrun Trust
| Axis | Assessment |
|---|---|
| Scientifically valid | **Yes.** The optimizer searches freely across the parameter space without checking whether the model is valid in the proposed region. |
| Logically consistent | **Yes.** Trust penalties and hard blockers in the objective function are the correct approach. |
| Technically feasible | **Yes.** Add a trust-penalty term to the objective function and block recommendations where evidence tier is "unsupported". |

**Verdict: AGREE. Priority 5.**

---

### Gap 10: Test Suite Needs Fast/Slow Tiering
| Axis | Assessment |
|---|---|
| Scientifically valid | N/A (engineering concern, not scientific). |
| Logically consistent | **Yes.** The full suite takes >120s due to Hypothesis property tests and CH PDE solves. |
| Technically feasible | **Yes.** Use pytest markers (`@pytest.mark.slow`) and CI configuration. |

**Verdict: AGREE. Low priority but useful.**

---

## 5. Assessment of Third-Party Architectural Rules

| Rule | Current Status | Agree? | Action |
|---|---|---|---|
| 1. Every number must know where it came from | Partial (model_tier/model_used exist, not on all outputs) | **YES** | Add ModelEvidenceTier to all result objects |
| 2. Every model must declare where it is valid | Partial (trust gates check some ranges, not all) | **YES** | Add ValidationDomain to model manifests |
| 3. Every fallback must preserve chemistry or fail | **Already implemented** (L3 raises NotImplementedError for unsupported Thiele regimes) | YES (already done) | Verify coverage |
| 4. Every empirical model must carry calibration provenance | Partial (L2 has model_tier, calibration_store exists) | **YES** | Extend calibration provenance to all empirical models |
| 5. Every mechanistic model must return diagnostics | Partial (Thiele modulus logged but not in result object) | **YES** | Add diagnostics dict to result objects |
| 6. Every optimization candidate must satisfy trust gates | **NOT implemented** | **YES** | Add trust penalty to optimizer |
| 7. Every downstream prediction must inherit upstream uncertainty | **NOT implemented** (UQ runs L1-L4 only, M2/M3 separate) | **YES** | Thread UQ through M2/M3 pipeline |
| 8. Every unit operation must conserve quantity | Partial (L1 mass conservation, ACS conservation in M2) | YES | Add conservation checks to L3, L4, M3 |
| 9. Every assay should be ingestible as data | **NOT implemented** | DEFER | v7.0 assay ingestion framework |
| 10. Every UI panel should show evidence level | Partial (trust banner, confidence tier caption) | **YES** | Show evidence tier next to each numeric output |

---

## 6. Proposals Rejected or Deferred

| Proposal | Reason for Rejection/Deferral |
|---|---|
| **ProcessDossier** (Phase 1) | Architecturally sound but v7.0 scope. Current SimulationParameters + FullResult + M1ExportContract provide sufficient provenance for v6.x. Add lightweight RunReport instead. |
| **Full unit validation framework** (Gap 6) | Disproportionate complexity. All internal calculations are SI-consistent. Add targeted interface assertions instead. |
| **6-component UQ decomposition** (Gap 7) | Research-grade feature. Current MC provides useful CIs. Add source tagging but defer posterior decomposition. |
| **AssayRecord ingestion** (Rule 9) | Requires UI redesign and data format specification. v7.0 feature. |

---

## 7. Consensus Modification Plan

### Sprint 1: ModelEvidenceTier + Output Metadata (Priority 1)

**Goal:** Every result object carries evidence tier, model source, and validity domain.

**What to build:**

1. **`ModelEvidenceTier` enum** in `datatypes.py`:
   ```
   VALIDATED_QUANTITATIVE   — calibrated against experimental data for this system
   CALIBRATED_LOCAL         — calibrated against data from an analogous system
   SEMI_QUANTITATIVE        — empirical model, not locally calibrated
   QUALITATIVE_TREND        — directional predictions only, not numeric
   UNSUPPORTED              — model not applicable to this chemistry/regime
   ```

2. **`ModelManifest` dataclass**:
   ```
   model_name: str          — "L1.PBE.FixedPivot.AlopaeusCT"
   evidence_tier: ModelEvidenceTier
   valid_domain: dict       — {"Re": (100, 1e6), "We": (1, 1e4), ...}
   calibration_ref: str     — "" or "CalibrationEntry.id"
   assumptions: list[str]
   ```

3. **Attach to every result dataclass:**
   - `EmulsificationResult.model_manifest`
   - `GelationResult.model_manifest` (replaces bare `model_tier` string)
   - `CrosslinkingResult.model_manifest`
   - `MechanicalResult.model_manifest` (replaces bare `model_used` string)

4. **RunReport dataclass** in `pipeline/orchestrator.py`:
   ```
   inputs: SimulationParameters
   model_graph: list[ModelManifest]   — one per level
   trust: TrustAssessment
   diagnostics: dict                  — Thiele, conservation, convergence
   ```

5. **UI rendering:** Show evidence tier badge next to each output metric.

**Model tier:** Opus (novel scientific model architecture)
**Estimated LOC:** ~200 (datatypes + orchestrator + UI)
**Tests:** 15 new tests (tier assignment, manifest construction, RunReport assembly)

---

### Sprint 2: L1 Trend Tests + Dimensionless Groups (Priority 2)

**Goal:** L1 droplet predictions are physically monotonic and dimensionlessly characterized.

**What to build:**

1. **`dimensionless_groups()` function** in `level1_emulsification/validation.py`:
   - Reynolds number: Re = rho * N * D^2 / mu
   - Weber number: We = rho * N^2 * D^3 / sigma
   - Capillary number: Ca = mu * v_tip / sigma
   - Viscosity ratio: lambda = mu_d / mu_c
   - Kolmogorov length: eta = (nu^3 / epsilon)^0.25
   - Droplet/Kolmogorov ratio: d32 / eta
   - Return as a typed dict attached to EmulsificationResult

2. **Monotonic trend regression tests** in `tests/test_l1_trends.py`:
   - Increasing RPM -> d32 non-increasing (within calibrated domain)
   - Increasing sigma -> d32 non-decreasing
   - Increasing mu_d -> d32 non-decreasing (viscous suppression)
   - Increasing c_span80 -> d32 non-increasing (more surfactant)
   - These are sign-constraint tests, not exact-value tests

3. **Domain gate**: If Re or We is outside the kernel's validated range, evidence tier downgrades to QUALITATIVE_TREND.

**Model tier:** Sonnet (standard numerical + test writing)
**Estimated LOC:** ~150
**Tests:** 8 trend tests + 4 dimensionless group tests

---

### Sprint 3: Calibration Service Generalization (Priority 3)

**Goal:** CalibrationStore can calibrate any model parameter, not just FMC fields.

**What to build:**

1. **Extend `CalibrationEntry`** with:
   - `target_module`: "L1" | "L2" | "L3" | "L4" | "M2" | "M3"
   - `target_parameter`: "breakage_C1" | "pore_A" | "k_xlink_0" | etc.
   - `fit_method`: "manual" | "least_squares" | "bayesian"
   - `valid_domain`: dict of parameter ranges where calibration applies
   - `posterior_uncertainty`: float (standard deviation of fitted parameter)

2. **Extend `CalibrationStore.apply()`** to dispatch by target_module:
   - L1: override KernelConfig constants
   - L2: override pore_A, pore_alpha, pore_beta
   - L3: override k_xlink_0, E_a_xlink, f_bridge
   - M2: override reaction rate constants
   - M3: override isotherm parameters

3. **Evidence tier upgrade:** When a model parameter is calibrated, upgrade the level's evidence tier from SEMI_QUANTITATIVE to CALIBRATED_LOCAL.

**Model tier:** Sonnet
**Estimated LOC:** ~200
**Tests:** 12 tests (multi-module calibration application, evidence tier upgrade)

---

### Sprint 4: L2 Label Fix + L3 Diagnostics Surfacing (Priority 4)

**Goal:** Pore predictions are honestly labeled; crosslinking diagnostics are in result objects.

**What to build:**

1. **L2 label fix**: Change default `model_tier` from `"empirical_calibrated"` to `"empirical_uncalibrated"`. Only set to `"empirical_calibrated"` when CalibrationStore has pore calibration data for this system.

2. **L3 diagnostics in result**: Add to `CrosslinkingResult`:
   - `thiele_modulus: float`
   - `regime: str` ("reaction_limited" | "borderline" | "diffusion_limited")
   - `stoichiometric_ceiling: float` (maximum conversion given crosslinker/reactive-group ratio)
   - `residual_reactive_groups: float` (mol/m3 unreacted -NH2 or -OH)

3. **L4 network classification**: Add to `MechanicalResult`:
   - `network_type: str` ("true_IPN" | "semi_IPN" | "independent_network" | "ionic_reinforced")
   - Cannot label as "true_IPN" unless primary and secondary crosslinkers target different networks AND covalent coupling exists.

4. **Mark CH tests as slow**: Add `@pytest.mark.slow` to CH 2D solver tests. Add `pytest.ini` config for `--ignore-slow` in fast CI.

**Model tier:** Sonnet (L2/L3 changes), Haiku (test markers)
**Estimated LOC:** ~150
**Tests:** 10 tests

---

### Sprint 5: M2/M3 Evidence Inheritance + Trust-Aware Optimization (Priority 5)

**Goal:** Downstream outputs inherit calibration status; optimizer respects trust.

**What to build:**

1. **M2/M3 evidence inheritance**: When M3 consumes a FunctionalMediaContract, check:
   - Is qmax calibrated or estimated? -> Tag M3 breakthrough output accordingly
   - Is isotherm constant from literature or measured? -> Inherit evidence tier
   - Is activity_retention the profile default (0.60) or user-calibrated? -> Warn if default

2. **M3 output evidence**: Add `evidence_tier: ModelEvidenceTier` to M3 result objects (breakthrough, DBC, elution).

3. **Optimizer trust penalty**:
   - Add `trust_penalty(evidence_tiers: list[ModelEvidenceTier]) -> float` to optimization objective
   - UNSUPPORTED: hard block (do not evaluate)
   - QUALITATIVE_TREND: +100 penalty (effectively excluded from Pareto)
   - SEMI_QUANTITATIVE: +10 penalty (deprioritized)
   - CALIBRATED_LOCAL: +0 (no penalty)
   - VALIDATED_QUANTITATIVE: -5 bonus (preferred)

4. **Pareto output labeling**: Each Pareto candidate annotated with min evidence tier across its outputs.

**Model tier:** Sonnet
**Estimated LOC:** ~200
**Tests:** 15 tests

---

### Sprint 6: Conservation Checks + UI Evidence Display (Priority 6)

**Goal:** Conservation enforced everywhere; evidence visible in UI.

**What to build:**

1. **L3 conservation check**: Verify that consumed reactive groups + residual = initial. Add to `CrosslinkingResult.diagnostics`.

2. **M2 ACS conservation**: Already exists. Add formal mass balance check (total sites consumed + remaining = initial, across all step types).

3. **M3 mass balance**: Check that eluted + bound + flowthrough = loaded (within numerical tolerance).

4. **UI evidence badges**:
   - Next to every metric in the results dashboard, show a colored badge:
     - Green: VALIDATED_QUANTITATIVE or CALIBRATED_LOCAL
     - Yellow: SEMI_QUANTITATIVE
     - Orange: QUALITATIVE_TREND
     - Red: UNSUPPORTED
   - Expandable "Evidence Details" panel showing model manifest for each level
   - "What experiments would improve this?" recommendation based on which inputs are uncalibrated

5. **RunReport JSON export**: Add download button for the structured RunReport as JSON, suitable for lab notebook attachment.

**Model tier:** Sonnet (conservation), Haiku (UI badges)
**Estimated LOC:** ~250
**Tests:** 12 tests

---

## 8. Execution Summary

| Sprint | Deliverable | Model Tier | Est. LOC | Est. Tests | Priority |
|---|---|---|---|---|---|
| 1 | ModelEvidenceTier + ModelManifest + RunReport | Opus | 200 | 15 | **Critical** |
| 2 | L1 dimensionless groups + trend tests | Sonnet | 150 | 12 | **High** |
| 3 | Cross-module calibration service | Sonnet | 200 | 12 | **High** |
| 4 | L2 label fix + L3 diagnostics + L4 classification | Sonnet | 150 | 10 | **Medium** |
| 5 | M2/M3 evidence inheritance + trust-aware optimizer | Sonnet | 200 | 15 | **Medium** |
| 6 | Conservation checks + UI evidence display + RunReport export | Sonnet/Haiku | 250 | 12 | **Medium** |
| **Total** | | | **1150** | **76** | |

---

## 9. What This Plan Does NOT Do (Deferred to v7.0+)

| Feature | Why Deferred |
|---|---|
| ProcessDossier | Requires new data model spanning experiments, assays, evidence, targets. v7.0. |
| Full unit validation (pint/unyt) | Disproportionate complexity for SI-consistent internal code. |
| 6-component UQ decomposition | Research-grade feature requiring posterior inference infrastructure. |
| AssayRecord ingestion | Requires UI for data upload, format parsing, validation. |
| Surrogate L2 model | Requires training data from CH sweeps + experiments. |
| Target Product Profile | Useful but requires application-specific definition (chromo vs enzyme vs adsorbent). |
| Wet-lab measurement objects | Large data model expansion. |

---

## 10. Architectural Rules Adopted

From the third-party report's 10 rules, we adopt **8 immediately** and defer 2:

| # | Rule | Status |
|---|---|---|
| 1 | Every number must know where it came from | **ADOPTED** (Sprint 1) |
| 2 | Every model must declare where it is valid | **ADOPTED** (Sprint 1-2) |
| 3 | Every fallback must preserve chemistry or fail | **ALREADY IMPLEMENTED** |
| 4 | Every empirical model must carry calibration provenance | **ADOPTED** (Sprint 3) |
| 5 | Every mechanistic model must return diagnostics | **ADOPTED** (Sprint 4) |
| 6 | Every optimization candidate must satisfy trust gates | **ADOPTED** (Sprint 5) |
| 7 | Every downstream prediction must inherit upstream uncertainty | **ADOPTED** (Sprint 5) |
| 8 | Every unit operation must conserve quantity | **ADOPTED** (Sprint 6) |
| 9 | Every assay should be ingestible as data | **DEFERRED** (v7.0) |
| 10 | Every UI panel should show evidence level | **ADOPTED** (Sprint 6) |

---

## 11. Conclusion

The third-party redesign report is a well-reasoned, scientifically literate assessment. Its central thesis -- "rebuild the architecture around existing solvers, not from scratch" -- is correct. The 10 gaps identified are largely real, with the caveat that the current codebase already addresses Gaps 3 (partial), 4 (partial), and 10 (partially) through existing trust gates, chemistry dispatch, and validation checks.

The consensus plan extracts the highest-value, most feasible modifications into 6 sprints totaling ~1150 LOC and 76 tests. This strengthens EmulSim's scientific validity (evidence tiers, trend tests, conservation), logical consistency (uniform metadata, no misleading labels), and technical feasibility (backward-compatible extensions, no rewrite) while deferring research-grade features to v7.0.

**The three reviewers (Scientific Advisor, Architect, Dev-Orchestrator) unanimously agree on this plan.**

---

> **Disclaimer**: This architectural assessment and modification plan is provided for development guidance. All proposed changes should be validated through testing and peer review before production deployment. Scientific models require experimental calibration before being used for decision-making.
