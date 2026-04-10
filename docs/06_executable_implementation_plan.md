# EmulSim Trust Warning Resolution -- Executable Implementation Plan

**Date:** 2026-04-11
**Prepared by:** Dev-Orchestrator (Opus 4.6)
**Source:** `docs/05_trust_warning_resolution_plan.md`
**Build order:** M1 -> M4 -> M2 -> M3 (Phases 1->2->3)

---

## 1. Task Stratification

### Overview

| Module | Description | Model Tier | Est. LOC | Context Budget (tokens) | Justification |
|--------|-------------|-----------|----------|------------------------|---------------|
| **M1** | PBE steady-state convergence | **Sonnet 4.6** | ~104 | ~9,800 | Standard algorithm (loop + convergence check), no novel science. Multiple branches but well-defined logic. MEDIUM complexity per Ref 02 heuristics (50-200 LOC, standard algorithms, domain-aware). |
| **M4** | Per-chemistry eta coupling | **Sonnet 4.6** | ~57 | ~5,700 | Pure data-flow plumbing -- reading existing fields and wiring them through. No new equations. MEDIUM complexity (multiple dispatch points, but straightforward). |
| **M2** | Crosslinker stoichiometry guidance | **Haiku** | ~24 | ~2,900 | Config change + simple arithmetic function + message enhancement. LOW complexity (<50 LOC, simple math, no domain expertise needed beyond the formula already in the plan). |
| **M3-P1** | Mode-aware warnings + metadata | **Sonnet 4.6** | ~14 | ~3,200 | Conditional logic in trust.py + dataclass fields. Requires understanding model modes. MEDIUM. |
| **M3-P2** | Hashin-Shtrikman bounds | **Sonnet 4.6** | ~40 | ~5,500 | Known composite-mechanics formula. Standard scientific coding -- equations are well-defined (Hashin 1963). No numerical solver. MEDIUM. |
| **M3-P3** | Flory-Rehner affine IPN model | **Opus 4.6** | ~210 | ~18,000 | Novel scientific computation: coupled nonlinear Flory-Rehner equilibrium + affine network theory. Requires deep domain knowledge (polymer physics), fsolve convergence management, robust fallback. HIGH complexity (>200 LOC, novel equations, deep domain expertise, numerical solver). |

### Token Budget Breakdown Per Module (from Ref 03 estimators)

```
Module   Protocol  Impl   Tests   Audit  Fix(x1)  Handoff  Subtotal  x1.5 Safety
M1       2500      416    250     1500   125      1000     5791      8687
M4       2000      228    137     1200   68       800      4433      6650
M2       1500      96     58      800    29       800      3283      4925
M3-P1    1500      56     34      1000   17       800      3407      5111
M3-P2    2000      160    96      1200   48       800      4304      6456
M3-P3    3000      840    504     2000   252      1200     7796      11694
```

---

## 2. Key Work Nodes

### M1: PBE Steady-State Convergence

#### M1-N1: Add adaptive convergence fields to datatypes

| Property | Value |
|----------|-------|
| **ID** | M1-N1 |
| **Description** | Add three fields to `SolverSettings` and one field to `EmulsificationResult` in `datatypes.py` |
| **Model tier** | Haiku |
| **Inputs** | Current `datatypes.py` (lines 534-552 for `SolverSettings`, lines 779-793 for `EmulsificationResult`) |
| **Outputs** | Modified `datatypes.py` |
| **Exact changes** | 1. In `SolverSettings` (after line 551, after `l3_atol`), add: `l1_t_max: float = 600.0` (max allowable L1 time [s]), `l1_conv_tol: float = 0.01` (d32 relative stability tolerance), `l1_max_extensions: int = 2` (max adaptive extension rounds). 2. In `EmulsificationResult` (after line 793, after `n_d_history`), add: `t_converged: Optional[float] = None` (time at which d32 stabilized, None if never). |
| **Acceptance criteria** | Fields exist with correct types and defaults. Existing code unaffected (all new fields have defaults). `SolverSettings()` and `EmulsificationResult(...)` construct without error. |
| **Est. LOC** | 4 |

#### M1-N2: Update default.toml

| Property | Value |
|----------|-------|
| **ID** | M1-N2 |
| **Description** | Change `t_emulsification` from 60.0 to 300.0 in `configs/default.toml` |
| **Model tier** | Haiku |
| **Inputs** | `configs/default.toml` line 7 |
| **Outputs** | Modified `configs/default.toml` |
| **Exact changes** | Line 7: `t_emulsification = 60.0` -> `t_emulsification = 300.0` with comment `# [s] (5 min -- allows PBE to reach steady state)` |
| **Acceptance criteria** | Config loads correctly. `params.emulsification.t_emulsification == 300.0` after `load_config()`. |
| **Est. LOC** | 1 |

#### M1-N3: Implement adaptive extension in `PBESolver.solve()`

| Property | Value |
|----------|-------|
| **ID** | M1-N3 |
| **Description** | Add adaptive time-extension loop to the legacy `solve()` method in `src/emulsim/level1_emulsification/solver.py` |
| **Model tier** | Sonnet 4.6 |
| **Inputs** | `solver.py` lines 207-317 (`solve()` method), `SolverSettings` with new fields from M1-N1 |
| **Outputs** | Modified `solve()` with extension loop |
| **Exact changes** | After the current `sol = solve_ivp(...)` block (line 283) and before `N_final` computation (line 288), wrap the convergence check + extension logic: (1) Compute convergence using `params.solver.l1_conv_tol` instead of hardcoded 0.01. (2) If not converged AND extension_count < `params.solver.l1_max_extensions` AND total_time < `params.solver.l1_t_max`: re-init from `sol.y[:, -1]`, extend by `t_emul / 2`, increment counter, re-integrate, append `sol.t` and `sol.y` histories. (3) Set `t_converged` to the time at which d32 first stabilized (within tolerance over last 10% window), or None. (4) Final `EmulsificationResult` constructor gets `t_converged=t_converged`. Volume conservation must be checked after each extension (assert `abs(total_vol_new / total_vol_init - 1.0) < 0.05`). (5) Log each extension: `logger.info("L1: extending by %.0fs (extension %d/%d), d32=%.2f um", ...)`. |
| **Acceptance criteria** | G-M1-1 through G-M1-5 from resolution plan. Specifically: default config converges; short t=5s triggers extension; t_history is monotonically increasing across extensions; volume conservation within 5%; W1 warning absent with defaults. |
| **Est. LOC** | ~40 |

#### M1-N4: Implement adaptive extension in `solve_stirred_vessel()`

| Property | Value |
|----------|-------|
| **ID** | M1-N4 |
| **Description** | Add identical adaptive extension logic to `solve_stirred_vessel()` |
| **Model tier** | Sonnet 4.6 |
| **Inputs** | `solver.py` lines 321-553 (`solve_stirred_vessel()`), M1-N3 as pattern |
| **Outputs** | Modified `solve_stirred_vessel()` |
| **Exact changes** | Mirror M1-N3 logic after the `sol = solve_ivp(...)` block (line 497). Key differences: (1) convergence tolerance from `params.solver.l1_conv_tol` (default 0.01) but the existing stirred-vessel code uses 0.05 -- keep 0.05 as the default for stirred-vessel mode by using `conv_tol = 0.05` locally (can be overridden later). (2) When extending, the temperature profile `T_profile` and `t_eval` must also be extended (append new segment). (3) The kernel cache (`g_cache`, `Q_cache`) does NOT need recomputation -- it is indexed by temperature, not time. |
| **Acceptance criteria** | Stirred-vessel mode also adaptively extends. Extension appends to `t_history` seamlessly. Temperature-dependent kernels interpolate correctly into the extended time range. |
| **Est. LOC** | ~35 |

#### M1-N5: Write convergence tests

| Property | Value |
|----------|-------|
| **ID** | M1-N5 |
| **Description** | Add tests for M1 convergence behavior in `tests/test_level1_emulsification.py` |
| **Model tier** | Sonnet 4.6 |
| **Inputs** | Existing test file, M1-N1 through M1-N4 completed |
| **Outputs** | New test functions in existing test file |
| **Exact changes** | Add functions: `test_default_config_converges()` -- load default.toml, run L1, assert `result.converged is True`. `test_short_time_triggers_extension()` -- set `t_emulsification=5.0`, run L1, assert `result.t_history[-1] > 5.0` (extension happened). `test_t_history_monotonic()` -- assert `np.all(np.diff(result.t_history) > 0)`. `test_volume_conservation()` -- assert `abs(result.total_volume_fraction / phi_d - 1.0) < 0.05`. `test_w1_absent_with_defaults()` -- run full pipeline with defaults, assert W1 text not in trust warnings. |
| **Acceptance criteria** | All 5 tests pass. No regressions in existing tests. |
| **Est. LOC** | ~25 |

---

### M4: Per-Chemistry Eta Coupling

#### M4-N1: Add `eta_coupling_recommended` to `NetworkTypeMetadata`

| Property | Value |
|----------|-------|
| **ID** | M4-N1 |
| **Description** | Add the `eta_coupling_recommended` field to `NetworkTypeMetadata` in `datatypes.py` |
| **Model tier** | Haiku |
| **Inputs** | `datatypes.py` lines 838-843 (`NetworkTypeMetadata`) |
| **Outputs** | Modified `NetworkTypeMetadata` |
| **Exact changes** | Add after line 843 (`is_true_second_network`): `eta_coupling_recommended: float = -0.15  # per-chemistry IPN coupling coefficient` |
| **Acceptance criteria** | Field exists with default -0.15. Backward compatible (default matches current universal value). |
| **Est. LOC** | 2 |

#### M4-N2: Wire eta from CrosslinkerProfile to NetworkTypeMetadata in L3 solver

| Property | Value |
|----------|-------|
| **ID** | M4-N2 |
| **Description** | At all 4 dispatch points in `solve_crosslinking()` where `NetworkTypeMetadata` is constructed, set `eta_coupling_recommended` from the `CrosslinkerProfile` |
| **Model tier** | Sonnet 4.6 |
| **Inputs** | `src/emulsim/level3_crosslinking/solver.py` lines 769-800 (dispatch metadata construction) |
| **Outputs** | Modified metadata construction at 4 points |
| **Exact changes** | At each `NetworkTypeMetadata(...)` construction in `solve_crosslinking()`: (1) **amine_covalent** (line 780-785): add `eta_coupling_recommended=xl.eta_coupling_recommended`. (2) **hydroxyl_covalent** (line 772-777): add `eta_coupling_recommended=xl.eta_coupling_recommended`. (3) **independent_network** (line 788-793): add `eta_coupling_recommended=xl.eta_coupling_recommended`. (4) **ionic_reversible** (line 796-800+): add `eta_coupling_recommended=xl.eta_coupling_recommended`. (5) Also in `_solve_reaction_diffusion()` (line 420-425): add `eta_coupling_recommended=xl.eta_coupling_recommended if xl else -0.15` -- but this path does not have access to `xl`. Instead, pass `xl` as a parameter to `_solve_reaction_diffusion()`, or set the metadata after return in the dispatcher. The simplest approach: in the dispatcher (line 763-766), after `rd_result` is returned, set `rd_result.network_metadata.eta_coupling_recommended = xl.eta_coupling_recommended` before returning. |
| **Acceptance criteria** | After calling `solve_crosslinking(..., crosslinker_key='dvs')`, result.network_metadata.eta_coupling_recommended == +0.10. Same for genipin -> -0.15, tpp -> -0.05, pegda_uv -> 0.0. |
| **Est. LOC** | ~10 |

#### M4-N3: Read per-chemistry eta in `select_modulus_model()` in L4 solver

| Property | Value |
|----------|-------|
| **ID** | M4-N3 |
| **Description** | Modify `select_modulus_model()` to use `network_metadata.eta_coupling_recommended` instead of hardcoded values |
| **Model tier** | Sonnet 4.6 |
| **Inputs** | `src/emulsim/level4_mechanical/solver.py` lines 74-103 (`select_modulus_model()`) |
| **Outputs** | Modified function |
| **Exact changes** | (1) In the `hydroxyl_covalent` branch (lines 95-99): replace `eta_coupling=+0.05` with `eta_coupling=getattr(network_metadata, 'eta_coupling_recommended', +0.05)`. Remove the dead `if hasattr(...)` block (lines 96-98). (2) In the `amine_covalent` default branch (lines 100-102): replace `eta_coupling` parameter with `eta_meta = getattr(network_metadata, 'eta_coupling_recommended', eta_coupling)` and pass `eta_meta` to `double_network_modulus()`. (3) In the `None` fallback (line 82): keep using `eta_coupling` parameter (no metadata available). (4) In `ionic_gel_modulus` branch: ionic model does not use eta, so no change needed. (5) In `independent_network` branch: `triple_network_modulus` uses `eta_13=0.0` by default; could use metadata, but PEGDA eta=0.0 matches, so no change needed for correctness. |
| **Acceptance criteria** | Genipin pipeline: G_DN unchanged (eta=-0.15 same as default). DVS pipeline: eta=+0.10 used, producing different G_DN than genipin. |
| **Est. LOC** | ~10 |

#### M4-N4: Suppress W4 trust warning when per-chemistry eta is active

| Property | Value |
|----------|-------|
| **ID** | M4-N4 |
| **Description** | Modify trust check #14 in `src/emulsim/trust.py` to only fire when the universal default is in use |
| **Model tier** | Haiku |
| **Inputs** | `trust.py` lines 177-183 (check #14) |
| **Outputs** | Modified condition |
| **Exact changes** | The current check fires when `props.eta_coupling == -0.15`. This is insufficient -- it should also check whether the pipeline actually used a per-chemistry value. The cleanest approach: add a `crosslinker_key` parameter (already present) and check whether the crosslinker has a non-default eta. Replace lines 177-183 with: (1) Look up xl from CROSSLINKERS. (2) If xl exists and xl.eta_coupling_recommended != -0.15, the pipeline will use the per-chemistry value, so skip the warning. (3) If xl is None or xl.eta_coupling_recommended == -0.15, keep the warning. Concretely: `if xl and xl.eta_coupling_recommended != -0.15: pass  # per-chemistry eta active, no warning` `else: warnings.append("IPN coupling coefficient (eta=-0.15) is the same default...")`. Note: `xl` is already fetched on line 138; reuse it. |
| **Acceptance criteria** | G-M4-5: W4 warning absent when crosslinker_key='dvs' or 'tpp' or 'pegda_uv'. W4 warning present when crosslinker_key='genipin' (because genipin eta IS -0.15). |
| **Est. LOC** | ~10 |

#### M4-N5: Write per-chemistry eta propagation tests

| Property | Value |
|----------|-------|
| **ID** | M4-N5 |
| **Description** | Add tests for eta propagation in `tests/test_level4_mechanical.py` and/or `tests/test_level3_crosslinking.py` |
| **Model tier** | Sonnet 4.6 |
| **Inputs** | M4-N1 through M4-N4 completed |
| **Outputs** | New test functions |
| **Exact changes** | Tests: `test_genipin_eta_backward_compatible()` -- run L3+L4 with genipin, assert eta=-0.15 used, G_DN matches pre-change value. `test_dvs_eta_positive()` -- run L3 with dvs, assert `result.network_metadata.eta_coupling_recommended == +0.10`. `test_tpp_eta()` -- assert eta=-0.05. `test_pegda_eta_zero()` -- assert eta=0.0. `test_gdn_differs_dvs_vs_genipin()` -- run both, assert different G_DN (different eta sign). `test_w4_suppressed_for_dvs()` -- run pipeline with dvs, assert W4 text absent from trust warnings. |
| **Acceptance criteria** | All 6 tests pass. G-M4-1 through G-M4-6. |
| **Est. LOC** | ~25 |

---

### M2: Crosslinker Stoichiometry Guidance

#### M2-N1: Update default c_genipin in config

| Property | Value |
|----------|-------|
| **ID** | M2-N1 |
| **Description** | Change `c_genipin` from 2.0 to 10.0 in `configs/default.toml` |
| **Model tier** | Haiku |
| **Inputs** | `configs/default.toml` line 21 |
| **Outputs** | Modified config |
| **Exact changes** | Line 21: `c_genipin = 2.0` -> `c_genipin = 10.0` with comment `# [mol/m3] (~10 mM, stoichiometric for target p=0.20)` |
| **Acceptance criteria** | Config loads. `params.formulation.c_genipin == 10.0`. |
| **Est. LOC** | 1 |

#### M2-N2: Add `recommended_crosslinker_concentration()` to L3 solver

| Property | Value |
|----------|-------|
| **ID** | M2-N2 |
| **Description** | Add a utility function to `src/emulsim/level3_crosslinking/solver.py` |
| **Model tier** | Haiku |
| **Inputs** | L3 solver, formulas from resolution plan |
| **Outputs** | New function |
| **Exact changes** | Add after `available_hydroxyl_concentration()` (after line 102): ```python def recommended_crosslinker_concentration(c_chitosan: float, DDA: float, M_GlcN: float, target_p: float = 0.20) -> float: """Minimum crosslinker concentration [mol/m3] for target crosslinking fraction. c_genipin_min = target_p * [NH2] / 2 Each crosslink bridge consumes 1 genipin and 2 NH2 groups. """ NH2 = available_amine_concentration(c_chitosan, DDA, M_GlcN) return target_p * NH2 / 2.0 ``` |
| **Acceptance criteria** | `recommended_crosslinker_concentration(18.0, 0.9, 161.16, 0.20)` returns approximately 10.05 mol/m3. `recommended_crosslinker_concentration(18.0, 0.9, 161.16, 0.30)` returns approximately 15.07 mol/m3. |
| **Est. LOC** | ~15 |

#### M2-N3: Enhance W2 trust warning message

| Property | Value |
|----------|-------|
| **ID** | M2-N3 |
| **Description** | Modify trust check #7 in `src/emulsim/trust.py` to include suggested minimum concentration |
| **Model tier** | Haiku |
| **Inputs** | `trust.py` lines 103-108 (check #7), M2-N2 function |
| **Outputs** | Enhanced warning message |
| **Exact changes** | Replace lines 104-108 with: (1) Compute `from .level3_crosslinking.solver import recommended_crosslinker_concentration`. (2) `c_min = recommended_crosslinker_concentration(params.formulation.c_chitosan, props.DDA, props.M_GlcN, target_p=0.05)`. (3) Warning message: `f"Crosslinker/NH2 ratio = {ratio:.3f} -- crosslinker-limited. Minimum c_genipin for ratio >= 0.05: {c_min:.1f} mol/m3 ({c_min:.1f} mM). Increasing crosslinking time will not help."` |
| **Acceptance criteria** | G-M2-1: default ratio >= 0.05 (because c_genipin=10 now). G-M2-2: W2 warning absent with defaults. G-M2-3: when warning fires (e.g., c_genipin=1.0), message includes suggested concentration. |
| **Est. LOC** | ~8 |

#### M2-N4: Write stoichiometry guidance tests

| Property | Value |
|----------|-------|
| **ID** | M2-N4 |
| **Description** | Add tests for M2 |
| **Model tier** | Haiku |
| **Inputs** | M2-N1 through M2-N3 completed |
| **Outputs** | New test functions in `tests/test_level3_crosslinking.py` |
| **Exact changes** | `test_default_ratio_above_threshold()` -- default config ratio >= 0.05. `test_w2_absent_with_defaults()` -- run pipeline, assert W2 text absent. `test_warning_includes_suggestion()` -- set c_genipin=1.0, assert "Minimum c_genipin" in warning text. `test_recommended_concentration()` -- call function with target_p=0.3, verify result is close to 15.07. |
| **Acceptance criteria** | G-M2-1 through G-M2-4. |
| **Est. LOC** | ~15 |

---

### M3: Predictive DN Modulus Model (Phased)

#### M3-P1-N1: Add model metadata fields to MechanicalResult

| Property | Value |
|----------|-------|
| **ID** | M3-P1-N1 |
| **Description** | Add `model_used`, `G_DN_lower`, `G_DN_upper` fields to `MechanicalResult` in `datatypes.py` |
| **Model tier** | Haiku |
| **Inputs** | `datatypes.py` lines 863-876 (`MechanicalResult`) |
| **Outputs** | Modified dataclass |
| **Exact changes** | After line 876 (`xi_mesh`), add: `model_used: str = "phenomenological"` (which model produced G_DN), `G_DN_lower: float = 0.0` (Hashin-Shtrikman lower bound), `G_DN_upper: float = 0.0` (Hashin-Shtrikman upper bound). |
| **Acceptance criteria** | Fields exist with defaults. Existing code unaffected. |
| **Est. LOC** | 4 |

#### M3-P1-N2: Mode-aware W3 trust warning

| Property | Value |
|----------|-------|
| **ID** | M3-P1-N2 |
| **Description** | Make W3 (phenomenological DN warning, trust check #12) conditional on model mode |
| **Model tier** | Sonnet 4.6 |
| **Inputs** | `trust.py` lines 163-167 (check #12) |
| **Outputs** | Modified check |
| **Exact changes** | Replace the unconditional `warnings.append(...)` on line 164 with conditional logic: (1) If `model_mode == 'empirical_engineering'`: suppress entirely (phenomenological is expected). (2) If `model_mode == 'hybrid_coupled'`: keep as INFO-level warning (current behavior). (3) If `model_mode == 'mechanistic_research'`: escalate to blocker if no affine model is active (check `m.model_used != 'flory_rehner_affine'`). Access `model_mode` from `params.model_mode`. Access `m.model_used` from `result.mechanical.model_used`. |
| **Acceptance criteria** | `empirical_engineering` mode: W3 text absent. `hybrid_coupled`: W3 present as warning. `mechanistic_research` with phenomenological model: W3 present as blocker. |
| **Est. LOC** | ~10 |

#### M3-P2-N1: Implement Hashin-Shtrikman bounds function

| Property | Value |
|----------|-------|
| **ID** | M3-P2-N1 |
| **Description** | Add `hashin_shtrikman_bounds()` to `src/emulsim/level4_mechanical/solver.py` |
| **Model tier** | Sonnet 4.6 |
| **Inputs** | Standard Hashin-Shtrikman (1963) composite mechanics formulas |
| **Outputs** | New function |
| **Exact changes** | Add after `triple_network_modulus()` (after line 71): ```python def hashin_shtrikman_bounds(G1: float, G2: float, phi1: float) -> tuple[float, float]: """Hashin-Shtrikman bounds for two-phase composite shear modulus [Pa]. G_lower = G_soft + phi_stiff / (1/(G_stiff - G_soft) + (1-phi_stiff)/(2*G_soft)) G_upper = G_stiff + (1-phi_stiff) / (1/(G_soft - G_stiff) + phi_stiff/(2*G_stiff)) Parameters ---------- G1, G2 : float Shear moduli of phases 1 and 2 [Pa]. phi1 : float Volume fraction of phase 1 [-]. Returns ------- G_lower, G_upper : float Lower and upper HS bounds [Pa]. """ phi2 = 1.0 - phi1 G_soft, G_stiff = sorted([G1, G2]) phi_stiff = phi1 if G1 >= G2 else phi2 phi_soft = 1.0 - phi_stiff if G_soft <= 0: return (0.0, G_stiff * phi_stiff) # Degenerate denom_lower = 1.0 / max(G_stiff - G_soft, 1e-30) + phi_soft / (2.0 * G_soft) G_lower = G_soft + phi_stiff / denom_lower denom_upper = 1.0 / max(G_soft - G_stiff, -1e30) + phi_stiff / (2.0 * G_stiff) # For upper bound, the denominator term is negative; use the standard form: G_upper = G_stiff + phi_soft / (1.0 / (G_soft - G_stiff) + phi_stiff / (2.0 * G_stiff)) return (max(G_lower, 0.0), max(G_upper, 0.0)) ``` Note: The exact Hashin-Shtrikman formulation must be verified against Hashin & Shtrikman (1963) during implementation. The /scientific-coder should validate with symmetric test cases. |
| **Acceptance criteria** | `hashin_shtrikman_bounds(G, G, 0.5)` returns `(G, G)` for any G > 0. `hashin_shtrikman_bounds(1000, 10000, 0.5)` returns `G_lower < G_upper`. Bounds contain the phenomenological G_DN value for standard formulations. |
| **Est. LOC** | ~30 |

#### M3-P2-N2: Populate bounds in `solve_mechanical()`

| Property | Value |
|----------|-------|
| **ID** | M3-P2-N2 |
| **Description** | Call `hashin_shtrikman_bounds()` in `solve_mechanical()` and populate `MechanicalResult` fields |
| **Model tier** | Sonnet 4.6 |
| **Inputs** | `solver.py` lines 176-247 (`solve_mechanical()`), M3-P2-N1 |
| **Outputs** | Modified `solve_mechanical()` |
| **Exact changes** | After G_DN computation (line 206) and before E_star (line 209), add: (1) Compute `phi1 = params.formulation.c_agarose / (params.formulation.c_agarose + params.formulation.c_chitosan)` (agarose volume fraction in polymer phase). (2) `G_lower, G_upper = hashin_shtrikman_bounds(G_agar, G_xlink, phi1)`. (3) In the `MechanicalResult(...)` constructor (line 236), add `model_used="phenomenological"`, `G_DN_lower=float(G_lower)`, `G_DN_upper=float(G_upper)`. |
| **Acceptance criteria** | G-M3-5: G_lower <= G_DN <= G_upper for standard formulations. Bounds are populated in every pipeline run. |
| **Est. LOC** | ~10 |

#### M3-P3-N1: Implement Flory-Rehner swelling equilibrium solver

| Property | Value |
|----------|-------|
| **ID** | M3-P3-N1 |
| **Description** | Add `flory_rehner_swelling()` to `src/emulsim/level4_mechanical/solver.py` |
| **Model tier** | Opus 4.6 |
| **Inputs** | Flory-Rehner theory, resolution plan specification |
| **Outputs** | New function |
| **Exact changes** | Add function that solves the coupled swelling equilibrium for two interpenetrating networks. Inputs: `nu_e1` (crosslink density of network 1), `nu_e2` (crosslink density of network 2), `chi1` (Flory-Huggins interaction parameter for network 1), `chi2` (for network 2), `phi1_0` (preparation volume fraction of network 1), `phi2_0` (preparation volume fraction of network 2). Outputs: `Q1, Q2` (equilibrium swelling ratios). Uses `scipy.optimize.fsolve` to solve the two coupled equations: `mu_mixing_1 + mu_elastic_1 = 0` and `mu_mixing_2 + mu_elastic_2 = 0` where `mu_mixing = RT * (ln(1-phi) + phi + chi * phi^2)` and `mu_elastic = RT * nu_e * V1 * (phi^(1/3) - phi/2)`. Must include: robust initial guess, convergence check, fallback to phenomenological on failure. |
| **Acceptance criteria** | Single-network limit: Q2 = 0 reproduces standard Flory-Rehner. Symmetric test: Q1 = Q2 when both networks identical. fsolve converges for chi in [0.4, 0.6] and nu_e in [1e22, 1e25]. |
| **Est. LOC** | ~60 |

#### M3-P3-N2: Implement affine double-network modulus

| Property | Value |
|----------|-------|
| **ID** | M3-P3-N2 |
| **Description** | Add `double_network_modulus_affine()` to L4 solver |
| **Model tier** | Opus 4.6 |
| **Inputs** | Affine network theory, M3-P3-N1 swelling solver |
| **Outputs** | New function |
| **Exact changes** | `G_DN = nu_e1 * kT * lambda_1^(-2/3) + nu_e2 * kT * lambda_2^(-2/3)` where `lambda_i = Q_i^(1/3)` is the linear stretch ratio from swelling. Inputs: `nu_e1, nu_e2, Q1, Q2, T`. Must check for physically reasonable outputs: `15 kPa < G_DN < 150 kPa` for standard formulations (warn if outside). |
| **Acceptance criteria** | G-M3-1: single-network recovery (G2=0 -> G_DN = G1 within 1%). G-M3-2: symmetric networks (G_DN >= G for identical networks). G-M3-3: literature validation within 2x for 4% agarose + 2% chitosan. |
| **Est. LOC** | ~80 |

#### M3-P3-N3: Update `select_modulus_model()` for mechanistic_research mode

| Property | Value |
|----------|-------|
| **ID** | M3-P3-N3 |
| **Description** | Route to affine model when `model_mode == MECHANISTIC_RESEARCH` |
| **Model tier** | Sonnet 4.6 |
| **Inputs** | M3-P3-N1, M3-P3-N2, L4 solver |
| **Outputs** | Modified `select_modulus_model()` or `solve_mechanical()` |
| **Exact changes** | In `solve_mechanical()`, before calling `select_modulus_model()`: check `params.model_mode`. If `MECHANISTIC_RESEARCH` and `family == 'amine_covalent'`: attempt `flory_rehner_swelling()` -> `double_network_modulus_affine()`. If fsolve fails, fall back to phenomenological and set `model_used = "phenomenological_fallback"`. If succeeds, set `model_used = "flory_rehner_affine"`. For other modes, keep current behavior with `model_used = "phenomenological"`. |
| **Acceptance criteria** | G-M3-4: fallback to phenomenological on fsolve failure. G-M3-6: empirical_engineering mode still uses phenomenological. |
| **Est. LOC** | ~20 |

#### M3-P3-N4: Write comprehensive M3 tests

| Property | Value |
|----------|-------|
| **ID** | M3-P3-N4 |
| **Description** | Add tests for all M3 phases |
| **Model tier** | Sonnet 4.6 |
| **Inputs** | All M3 nodes completed |
| **Outputs** | New test functions in `tests/test_level4_mechanical.py` |
| **Exact changes** | Tests: `test_hs_bounds_symmetric()` -- equal moduli -> bounds collapse. `test_hs_bounds_contain_gdn()` -- bounds bracket G_DN. `test_single_network_recovery()` -- G_DN(G1, 0) == G1 within 1%. `test_symmetric_networks()` -- G_DN(G, G) >= G. `test_affine_vs_literature()` -- 4% agarose + 2% chitosan within 2x of published value. `test_fsolve_fallback()` -- extreme params trigger fallback, no crash. `test_empirical_mode_uses_phenomenological()` -- model_used == "phenomenological". `test_mechanistic_mode_uses_affine()` -- model_used == "flory_rehner_affine". `test_w3_suppressed_empirical()` -- W3 absent in empirical mode. `test_bounds_populated()` -- G_DN_lower > 0 and G_DN_upper > 0 in every result. |
| **Acceptance criteria** | All 10 tests pass. G-M3-1 through G-M3-6. |
| **Est. LOC** | ~50 |

---

## 3. Multi-Agent Parallelism Map

### Dependency Graph

```
M1-N1 (datatypes)  ──> M1-N2 (config) ──> M1-N3 (solve) ──> M1-N4 (solve_stirred) ──> M1-N5 (tests)
                                              │
                                              │ (M1-N3 and M1-N4 can run in parallel
                                              │  if two agents are available)
                                              ▼
                                          M1-N5 (tests, depends on N3+N4)

M4-N1 (datatypes) ──> M4-N2 (L3 wiring) ──> M4-N3 (L4 wiring) ──> M4-N4 (trust) ──> M4-N5 (tests)

M2-N1 (config) ──> M2-N2 (function) ──> M2-N3 (trust msg) ──> M2-N4 (tests)

M3-P1-N1 (datatypes) ──> M3-P1-N2 (trust) ──> M3-P2-N1 (HS bounds) ──> M3-P2-N2 (populate) ──>
  M3-P3-N1 (Flory-Rehner) ──> M3-P3-N2 (affine modulus) ──> M3-P3-N3 (routing) ──> M3-P3-N4 (tests)
```

### Parallelism Opportunities

| Parallel Group | Nodes | Condition |
|---------------|-------|-----------|
| **PG-1** | M1-N1 + M4-N1 + M3-P1-N1 | All modify `datatypes.py` -- **CANNOT parallelize** (file conflict). Execute sequentially in one batch: M1-N1 then M4-N1 then M3-P1-N1. |
| **PG-2** | M1-N2 + M2-N1 | Both modify `configs/default.toml` -- **CANNOT parallelize**. Execute in one batch. |
| **PG-3** | M1-N3 + M1-N4 | Both modify `solver.py` (L1) but different methods -- could parallelize with careful merge, but **safer to serialize** (same file). |
| **PG-4** | M4-N2 + M4-N3 | Different files (L3 vs L4 solver) -- **CAN parallelize**. |
| **PG-5** | M1-N5 + M4-N5 + M2-N4 | Different test files -- **CAN parallelize** after their respective modules complete. |

### Recommended Single-Agent Execution Sequence

Since all datatype and config changes must serialize, and most changes touch shared files, the optimal single-agent sequence is:

```
Phase A: Datatypes + Config (all modules, one pass through each file)
  A1: M1-N1 + M4-N1 + M3-P1-N1  (datatypes.py -- 3 edits in one pass)
  A2: M1-N2 + M2-N1              (default.toml -- 2 edits in one pass)

Phase B: M1 Implementation
  B1: M1-N3  (L1 solver -- solve())
  B2: M1-N4  (L1 solver -- solve_stirred_vessel())

Phase C: M4 Implementation
  C1: M4-N2  (L3 solver wiring)
  C2: M4-N3  (L4 solver wiring)
  C3: M4-N4  (trust.py)

Phase D: M2 Implementation
  D1: M2-N2  (L3 solver -- new function)
  D2: M2-N3  (trust.py -- enhance warning)

Phase E: M3 Phases 1-2
  E1: M3-P1-N2  (trust.py -- mode-aware warning)
  E2: M3-P2-N1  (L4 solver -- HS bounds)
  E3: M3-P2-N2  (L4 solver -- populate bounds)

Phase F: M3 Phase 3 (Opus required)
  F1: M3-P3-N1  (Flory-Rehner swelling)
  F2: M3-P3-N2  (affine modulus)
  F3: M3-P3-N3  (model routing)

Phase G: All Tests
  G1: M1-N5  (L1 convergence tests)
  G2: M4-N5  (eta propagation tests)
  G3: M2-N4  (stoichiometry tests)
  G4: M3-P3-N4  (mechanical model tests)
```

### Collaboration Contracts

| Interface | Producer | Consumer | Contract |
|-----------|----------|----------|----------|
| `SolverSettings.l1_t_max/l1_conv_tol/l1_max_extensions` | M1-N1 | M1-N3, M1-N4 | Fields exist with correct defaults before solver edits begin |
| `EmulsificationResult.t_converged` | M1-N1 | M1-N3, M1-N4 | Field exists as `Optional[float] = None` |
| `NetworkTypeMetadata.eta_coupling_recommended` | M4-N1 | M4-N2, M4-N3 | Field exists with default -0.15 |
| `MechanicalResult.model_used/G_DN_lower/G_DN_upper` | M3-P1-N1 | M3-P2-N2, M3-P3-N3 | Fields exist with correct defaults |
| `recommended_crosslinker_concentration()` | M2-N2 | M2-N3 | Function importable from `emulsim.level3_crosslinking.solver` |
| `hashin_shtrikman_bounds()` | M3-P2-N1 | M3-P2-N2 | Function returns `(float, float)` tuple |
| `flory_rehner_swelling()` | M3-P3-N1 | M3-P3-N2 | Returns `(Q1, Q2)` floats; raises `RuntimeError` on convergence failure |
| `double_network_modulus_affine()` | M3-P3-N2 | M3-P3-N3 | Returns `float` G_DN; raises `RuntimeError` on invalid output |

---

## 4. Quality Gate Checklist

### M1: PBE Steady-State Convergence

| Gate ID | Type | Criterion | Test Method |
|---------|------|-----------|-------------|
| G-M1-1 | Functional | Default config (t=300s) produces `converged=True` | Unit test: load default.toml, run L1, assert converged |
| G-M1-2 | Functional | Short t=5s triggers adaptive extension, eventually converges | Unit test: set t_emulsification=5, run L1, assert t_history[-1] > 5 |
| G-M1-3 | Data integrity | `t_history` monotonically increasing after extension | Unit test: assert np.all(np.diff(result.t_history) > 0) |
| G-M1-4 | Physical | Volume conservation within 5% of phi_d | Unit test: assert abs(total_vol/phi_d - 1) < 0.05 |
| G-M1-5 | Integration | Trust warning W1 no longer fires with defaults | Integration test: run full pipeline, assert "PBE solver did not reach steady state" absent |
| G-M1-A1 | Architectural | Extension loop has max cap (no infinite loops) | Code review: l1_max_extensions enforced |
| G-M1-A2 | Architectural | Extension uses re-init from N_final (not restart from scratch) | Code review: verify N0 = sol.y[:, -1] |
| G-M1-A3 | Performance | Runtime with default config < 60s | Timing test (manual) |

### M4: Per-Chemistry Eta Coupling

| Gate ID | Type | Criterion | Test Method |
|---------|------|-----------|-------------|
| G-M4-1 | Backward compat | Genipin -> eta = -0.15 | Unit test |
| G-M4-2 | Functional | DVS -> eta = +0.10 | Unit test |
| G-M4-3 | Functional | TPP -> eta = -0.05 | Unit test |
| G-M4-4 | Functional | PEGDA -> eta = 0.0 | Unit test |
| G-M4-5 | Integration | W4 suppressed when per-chemistry eta active | Integration test |
| G-M4-6 | Physical | G_DN differs between DVS and genipin | Unit test: different eta sign -> different G_DN |
| G-M4-A1 | Architectural | NetworkTypeMetadata carries eta through L3->L4 cleanly | Code review: no hardcoded eta in L4 |
| G-M4-A2 | Architectural | Reaction-diffusion PDE path also gets eta wiring | Code review: verify line ~420 area in L3 solver |

### M2: Crosslinker Stoichiometry Guidance

| Gate ID | Type | Criterion | Test Method |
|---------|------|-----------|-------------|
| G-M2-1 | Functional | Default config ratio >= 0.05 | Unit test |
| G-M2-2 | Integration | W2 warning absent with defaults | Integration test |
| G-M2-3 | UX | Warning message includes suggested minimum concentration | String match test |
| G-M2-4 | Functional | `recommended_crosslinker_concentration(..., target_p=0.3)` -> ~15.07 | Unit test |
| G-M2-A1 | Backward compat | Old configs with explicit c_genipin=2.0 still work | Manual test |

### M3: Predictive DN Modulus Model

| Gate ID | Type | Criterion | Test Method |
|---------|------|-----------|-------------|
| G-M3-1 | Physical | Single network recovery: G_DN(G1, 0) = G1 +/- 1% | Unit test |
| G-M3-2 | Physical | Symmetric networks: G_DN(G, G) >= G | Unit test |
| G-M3-3 | Validation | Absolute prediction within 2x of literature for standard formulation | Unit test with known values |
| G-M3-4 | Robustness | Fallback to phenomenological when fsolve fails | Unit test with extreme params |
| G-M3-5 | Physical | G_lower <= G_DN <= G_upper (bounds containment) | Unit test |
| G-M3-6 | Modal | `empirical_engineering` mode still uses phenomenological | Unit test |
| G-M3-A1 | Architectural | Flory-Rehner solver has bounded iteration count | Code review |
| G-M3-A2 | Architectural | Sanity bounds enforced (15 kPa < G_DN < 150 kPa warning) | Code review |
| G-M3-A3 | Numerical | fsolve convergence verified (infodict check) | Code review |

---

## 5. Context Management Strategy

### Budget Allocation

| Phase | Cumulative Est. Tokens | Zone | Action |
|-------|----------------------|------|--------|
| Phase A (datatypes+config) | ~2,000 | GREEN | Proceed |
| Phase B (M1 impl) | ~8,000 | GREEN | Proceed |
| Phase C (M4 impl) | ~14,000 | GREEN | Proceed |
| Phase D (M2 impl) | ~17,000 | GREEN | Proceed |
| Phase E (M3 P1-P2) | ~24,000 | GREEN | Proceed |
| Phase F (M3 P3) | ~42,000 | GREEN/YELLOW | Pre-check: if context < 40%, compress before starting Phase F |
| Phase G (all tests) | ~52,000 | YELLOW | May need compression before final test phase |

### Compression Triggers

1. **Before Phase F (M3-P3):** If context consumed > 60%, execute Dialogue Compression Protocol. M3-P3 is the largest single unit (~18k tokens) and cannot be interrupted.

2. **Before Phase G (tests):** If context consumed > 70%, compress by summarizing all implementation details to one-line-per-module and carry forward only the test specs.

3. **Emergency:** If context hits RED (<30%) at any point, produce Milestone Handover with completed modules registered and next-module protocol pre-generated.

### Session Boundary Plan

If context is insufficient for all phases in one session:

- **Natural split point 1:** After Phase C (M1 + M4 complete). Produce handover with M2 + M3 protocols pre-generated.
- **Natural split point 2:** After Phase E (M1 + M4 + M2 + M3 P1-P2 complete). Produce handover with M3-P3 protocol pre-generated.
- **Each split:** Ensure all completed modules are tested before handover.

---

## 6. Module Registry Template

| # | Module | Sub-node | Version | Status | Model Tier | Fix Rounds | LOC | File(s) Modified |
|---|--------|----------|---------|--------|-----------|------------|-----|-----------------|
| 1 | M1 | M1-N1 | v1.0 | PENDING | Haiku | - | 4 | `src/emulsim/datatypes.py` |
| 2 | M1 | M1-N2 | v1.0 | PENDING | Haiku | - | 1 | `configs/default.toml` |
| 3 | M1 | M1-N3 | v1.0 | PENDING | Sonnet | - | 40 | `src/emulsim/level1_emulsification/solver.py` |
| 4 | M1 | M1-N4 | v1.0 | PENDING | Sonnet | - | 35 | `src/emulsim/level1_emulsification/solver.py` |
| 5 | M1 | M1-N5 | v1.0 | PENDING | Sonnet | - | 25 | `tests/test_level1_emulsification.py` |
| 6 | M4 | M4-N1 | v1.0 | PENDING | Haiku | - | 2 | `src/emulsim/datatypes.py` |
| 7 | M4 | M4-N2 | v1.0 | PENDING | Sonnet | - | 10 | `src/emulsim/level3_crosslinking/solver.py` |
| 8 | M4 | M4-N3 | v1.0 | PENDING | Sonnet | - | 10 | `src/emulsim/level4_mechanical/solver.py` |
| 9 | M4 | M4-N4 | v1.0 | PENDING | Haiku | - | 10 | `src/emulsim/trust.py` |
| 10 | M4 | M4-N5 | v1.0 | PENDING | Sonnet | - | 25 | `tests/test_level4_mechanical.py` |
| 11 | M2 | M2-N1 | v1.0 | PENDING | Haiku | - | 1 | `configs/default.toml` |
| 12 | M2 | M2-N2 | v1.0 | PENDING | Haiku | - | 15 | `src/emulsim/level3_crosslinking/solver.py` |
| 13 | M2 | M2-N3 | v1.0 | PENDING | Haiku | - | 8 | `src/emulsim/trust.py` |
| 14 | M2 | M2-N4 | v1.0 | PENDING | Haiku | - | 15 | `tests/test_level3_crosslinking.py` |
| 15 | M3 | M3-P1-N1 | v1.0 | PENDING | Haiku | - | 4 | `src/emulsim/datatypes.py` |
| 16 | M3 | M3-P1-N2 | v1.0 | PENDING | Sonnet | - | 10 | `src/emulsim/trust.py` |
| 17 | M3 | M3-P2-N1 | v1.0 | PENDING | Sonnet | - | 30 | `src/emulsim/level4_mechanical/solver.py` |
| 18 | M3 | M3-P2-N2 | v1.0 | PENDING | Sonnet | - | 10 | `src/emulsim/level4_mechanical/solver.py` |
| 19 | M3 | M3-P3-N1 | v1.0 | PENDING | **Opus** | - | 60 | `src/emulsim/level4_mechanical/solver.py` |
| 20 | M3 | M3-P3-N2 | v1.0 | PENDING | **Opus** | - | 80 | `src/emulsim/level4_mechanical/solver.py` |
| 21 | M3 | M3-P3-N3 | v1.0 | PENDING | Sonnet | - | 20 | `src/emulsim/level4_mechanical/solver.py` |
| 22 | M3 | M3-P3-N4 | v1.0 | PENDING | Sonnet | - | 50 | `tests/test_level4_mechanical.py` |

**Total nodes:** 22
**Total estimated LOC:** ~455
**Opus-required nodes:** 2 (M3-P3-N1, M3-P3-N2)
**Sonnet nodes:** 12
**Haiku nodes:** 8

---

## Appendix A: File Modification Summary

| File | Modules Touching It | Total Est. Changes |
|------|--------------------|--------------------|
| `src/emulsim/datatypes.py` | M1-N1, M4-N1, M3-P1-N1 | +10 lines (3 dataclass edits) |
| `configs/default.toml` | M1-N2, M2-N1 | 2 line changes |
| `src/emulsim/level1_emulsification/solver.py` | M1-N3, M1-N4 | +75 lines (2 method edits) |
| `src/emulsim/level3_crosslinking/solver.py` | M4-N2, M2-N2 | +25 lines (metadata wiring + new function) |
| `src/emulsim/level4_mechanical/solver.py` | M4-N3, M3-P2-N1, M3-P2-N2, M3-P3-N1, M3-P3-N2, M3-P3-N3 | +210 lines (1 edit + 4 new functions) |
| `src/emulsim/trust.py` | M4-N4, M2-N3, M3-P1-N2 | ~28 lines modified |
| `tests/test_level1_emulsification.py` | M1-N5 | +25 lines |
| `tests/test_level3_crosslinking.py` | M2-N4 | +15 lines |
| `tests/test_level4_mechanical.py` | M4-N5, M3-P3-N4 | +75 lines |

---

## Appendix B: Critical Line References

These are the exact locations where changes must be made, verified against the current codebase:

| Change | File | Line(s) | Current Content |
|--------|------|---------|-----------------|
| Add SolverSettings fields | `datatypes.py` | After 551 | `l3_atol: float = 1e-10` |
| Add EmulsificationResult.t_converged | `datatypes.py` | After 793 | `n_d_history: Optional[np.ndarray] = None` |
| Add NetworkTypeMetadata.eta_coupling_recommended | `datatypes.py` | After 843 | `is_true_second_network: bool = True` |
| Add MechanicalResult fields | `datatypes.py` | After 876 | `xi_mesh: float` |
| Change t_emulsification | `default.toml` | 7 | `t_emulsification = 60.0` |
| Change c_genipin | `default.toml` | 21 | `c_genipin = 2.0` |
| L1 solve() extension | `solver.py` (L1) | 276-317 | BDF solve_ivp + convergence check block |
| L1 solve_stirred() extension | `solver.py` (L1) | 490-553 | BDF solve_ivp + convergence check block |
| L3 metadata eta wiring | `solver.py` (L3) | 769-800 | 4 NetworkTypeMetadata constructors |
| L3 reaction-diffusion eta | `solver.py` (L3) | 420-425 | NetworkTypeMetadata without eta field |
| L4 select_modulus_model hydroxyl | `solver.py` (L4) | 95-99 | Hardcoded `eta_coupling=+0.05` |
| L4 select_modulus_model amine | `solver.py` (L4) | 100-102 | Uses `eta_coupling` parameter directly |
| Trust check #7 (W2) | `trust.py` | 103-108 | Simple ratio warning |
| Trust check #12 (W3) | `trust.py` | 163-167 | Unconditional append |
| Trust check #14 (W4) | `trust.py` | 177-183 | `props.eta_coupling == -0.15` check |

---

> **This plan is executable by /scientific-coder without ambiguity. Each node specifies exact file paths, line numbers, function names, and acceptance criteria. The build order (M1 -> M4 -> M2 -> M3) is preserved within the optimized execution sequence.**
