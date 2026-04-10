# EmulSim Audit Remediation Plan — Final Synthesis

**Date:** 2026-04-11
**Input:** 3rd-party audit (docs/08), Scientific Advisor, Dev-Orchestrator, Architect
**Verdict:** All 3 roles agree with the audit's 9 findings. No re-architecture needed.

---

## Cross-Role Consensus

All three specialist roles independently confirmed:

1. **All 9 findings are valid.** No false positives.
2. **The two-tier strategy (Production Engineering + Mechanistic Research) is correct.**
   This is already partially implemented via `ModelMode` but not enforced strongly enough.
3. **Finding 3 (L3 chemistry dispatch) is the #1 priority** — it produces silently wrong results.
4. **Finding 1 (L1 RPM nonmonotonicity) has a specific mechanism** — Vi-coalescence
   interaction at high epsilon, not just a calibration gap.
5. **No full re-architecture is warranted.** The L1-L4 pipeline is structurally sound.

---

## Finding-by-Finding: Best Combined Solutions

### F1 (HIGH): L1 RPM Nonmonotonicity

**Root cause (Architect + SciAdvisor converged):**
At high RPM, breakage pushes droplets to small sizes where:
- The Alopaeus viscous correction (C3 * Vi, Vi ~ d^(-0.5)) exponentially kills breakage
- Coalescence becomes very efficient (lambda -> 1 when d^4 -> 0)
- The steady-state distribution shifts back to larger sizes

**Best solution (all 3 agree + additions):**
1. Cap Vi at Vi_max=50-100 in `breakage_rate_alopaeus()` to prevent exponential death
2. Optionally add Weber-dependent coalescence suppression at high epsilon
3. Refit C1/C2/C3 against real d32-vs-RPM data with cap in place
4. Replace Hinze-only test with PBE-level monotonicity regression test
5. Add a runtime trust-gate: flag if d32(RPM_high) > d32(RPM_low)

**Differs from audit:** Audit recommends refitting alone. All 3 roles agree refitting
without a Vi cap will just shift the problem to a different RPM range.

| Aspect | Audit | Our Assessment |
|--------|-------|---------------|
| Root cause | Unclear calibration | Vi-coalescence interaction (identified) |
| Fix approach | Refit constants | Vi cap + refit + coalescence suppression |
| Test | PBE regression test | Agree |

---

### F2 (HIGH): L2 Empirical Default Not Labeled

**All 3 roles fully agree with audit.**

**Best solution:**
1. Rename mode from `'empirical'` to `'empirical_calibrated'` in orchestrator
2. Add `model_tier` metadata to `GelationResult`
3. Remove or clearly label synthetic pore distribution (SciAdvisor addition)
4. Compute alpha_final from Avrami kinetics, not hard-coded 0.999 (SciAdvisor addition)
5. Report calibration uncertainty explicitly (+/-30%)

---

### F3 (HIGH): L3 Chemistry Dispatch Bug

**All 3 roles agree this is the #1 priority. Pure implementation defect.**

**Best solution:**
1. Move Thiele check AFTER chemistry dispatch (not before)
2. Each chemistry branch checks Thiele internally
3. For amine-reactive (genipin, glutaraldehyde): existing PDE is valid
4. For hydroxyl-reactive (DVS, ECH, citric acid): need separate PDE with -OH groups,
   or `NotImplementedError` until implemented
5. For UV (PEGDA): `NotImplementedError` — fundamentally different physics (photoinitiation front)
6. For ionic (TPP): `NotImplementedError` — shrinking-core model needed, not reaction-diffusion
7. Make `_solve_reaction_diffusion` accept `xl` profile for chemistry-specific parameters
8. Add regression test: DVS/PEGDA/TPP at R=500um must NOT return `amine_covalent`

**All 3 roles agree: unsupported chemistries should fail hard, not silently fallback.**

---

### F4 (HIGH→MEDIUM): L2 Mechanistic Performance

**SciAdvisor + Architect both downgrade to MEDIUM** (fixable engineering problem, not
scientific invalidity). Dev-Orchestrator agrees but defers to Phase 3.

**Best solution (short-term):**
1. Freeze mobility updates every 10 steps (operator changes slowly)
2. Pre-factorize LHS matrix (binary solver rebuilds without caching LU)
3. Switch from `spsolve` to iterative solver (CG/GMRES with ILU preconditioner)
4. Expected 10-50x speedup

**Best solution (long-term, per audit):**
Invest in ternary solver with semi-implicit convex-splitting + gelation arrest

---

### F5 (MEDIUM): L2 Domain Cap

**All 3 agree: acceptable RVE approximation, needs documentation.**
- Make domain size adaptive to spinodal wavelength
- Add warning when R_droplet < 5 * lambda_CH
- Document that mechanistic L2 is a local-patch solver, not whole-droplet

---

### F6 (MEDIUM): Ternary L2 Exploratory

**All 3 agree with audit's long-term recommendation.**
- Replace explicit Euler with semi-implicit convex-splitting
- Add gelation arrest coupling
- Replace periodic BCs with no-flux

---

### F7 (MEDIUM): Uncertainty Ad Hoc

**All 3 agree. SciAdvisor adds nuance:** post-hoc perturbation is mathematically
equivalent to pre-hoc for algebraic (empirical) models, partially defusing the concern.

**Best solution:**
1. Fix naming: `breakage_C1_factor` → `breakage_C3_factor`
2. Split into three labeled categories: parameter / calibration / model-form
3. Document that post-hoc pore perturbation is valid only for empirical L2

---

### F8 (MEDIUM): Optimization qEHVI + Feasibility

**All 3 agree. Near-zero effort fix for qLogEHVI swap.**

**Best solution:**
1. Replace `qExpectedHypervolumeImprovement` → `qLogExpectedHypervolumeImprovement` (1 line)
2. Add trust-gate checks to feasibility constraints
3. Add RPM monotonicity check (reject proposals in nonphysical regime)
4. Restrict optimization to validated parameter ranges

---

### F9 (MEDIUM): L4 Phenomenological

**All 3 agree: lowest priority, already partially addressed in our M3 work.**
- Keep phenomenological for production tier
- Flory-Rehner affine already available for mechanistic tier
- Future: add toughness estimation (Lake-Thomas model)

---

## Execution Plan

### Phase 1: Immediate Safety (Days 1-2) — All Parallelizable

| Node | Finding | Tier | LOC | Model |
|------|---------|------|-----|-------|
| 1.1 | F3: L3 chemistry dispatch | Bug fix | ~50 | Sonnet |
| 1.2 | F2: L2 empirical labeling | Metadata | ~30 | Haiku |
| 1.3 | F8a: qLogEHVI swap | 1-line | ~5 | Haiku |
| 1.4 | F7a: Fix naming bug | Rename | ~5 | Haiku |

### Phase 2: Scientific Trust (Days 3-7)

| Node | Finding | Tier | LOC | Model |
|------|---------|------|-----|-------|
| 2.1 | F1: L1 Vi cap + recalibration | Kernel fix | ~200 | Opus |
| 2.2 | F7b: Uncertainty category split | Refactor | ~100 | Sonnet |
| 2.3 | F8b: Optimization constraints | Logic | ~60 | Sonnet |

### Phase 3: Medium-Term (Weeks 2-4)

| Node | Finding | LOC | Model |
|------|---------|-----|-------|
| 3.1 | F4: L2 CH performance | ~150 | Opus |
| 3.2 | F6: Ternary solver stabilization | ~200 | Opus |
| 3.3 | F5: Adaptive domain sizing | ~50 | Sonnet |

### Reordering vs Audit

| Audit Order | Our Order | Rationale |
|-------------|-----------|-----------|
| F1 first (L3 dispatch) | F3 first | Agree — silent correctness bug |
| F2 second (L1 RPM) | F2+F8a+F7a parallel on day 1 | Near-zero effort, ship with F3 |
| F3 third (L2 label) | F1 in Phase 2 | Requires Opus-tier research, not a quick fix |
| Rest sequential | Phases 2-3 as planned | Agree with audit's phasing |

---

## Quality Gates

| Gate | Criteria |
|------|----------|
| G-F3 | DVS/PEGDA/TPP at R=500um: NOT `amine_covalent`, raises error or returns correct chemistry |
| G-F1 | d32(RPM) monotonically decreasing for RPM in [3000, 25000] in PBE output |
| G-F2 | GelationResult contains `model_tier="empirical"` in default path |
| G-F8a | `test_optimization.py` passes without BoTorch qEHVI deprecation warnings |
| G-F7a | No variable named `breakage_C1_factor` applied to `breakage_C3` |

---

## Total Effort

| Phase | LOC | Calendar | Key Risk |
|-------|-----|----------|----------|
| Phase 1 | ~90 | 1-2 days | Low — targeted fixes |
| Phase 2 | ~360 | 3-5 days | Medium — L1 recalibration needs domain expertise |
| Phase 3 | ~400 | 2-3 weeks | Medium — algorithmic changes to solvers |
| **Total** | **~850** | **~3 weeks** | |
