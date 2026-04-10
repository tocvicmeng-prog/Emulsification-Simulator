# EmulSim Trust Warning Resolution Plan

**Date:** 2026-04-11
**Synthesised from:** Scientific Advisor, Dev-Orchestrator, Architect
**Pipeline version:** v3.0 (commit 48e7e9a)
**Scope:** 4 trust warnings from default-parameter pipeline run

---

## Executive Summary

A default-config pipeline run produced 4 trust warnings. Three specialist roles
analysed the codebase independently and converged on a unified resolution plan.
All three agree: **no re-architecture is needed** — the pipeline structure is sound.
The warnings arise from (1) a parameter tuning gap, (2) a formulation guidance gap,
(3) a standing model limitation, and (4) a data-flow wiring gap.

| # | Warning | Severity | Root Cause | Fix Type |
|---|---------|----------|------------|----------|
| W1 | PBE solver not at steady state | LOW | `t_emulsification=60s` too short | Parameter + adaptive extension |
| W2 | Crosslinker/NH2 ratio = 0.020 | MODERATE | `c_genipin=2 mM` stoichiometrically starved | Parameter + guidance function |
| W3 | Phenomenological G_DN formula | MOD–HIGH | No first-principles IPN model | Phased model upgrade |
| W4 | Universal eta = −0.15 | MODERATE | Per-chemistry eta exists but not wired | Data-flow plumbing |

**Total estimated effort:** ~397 LOC across 4 modules, ~13 hours.

---

## Module Plan

### M1: PBE Steady-State Convergence

**Scientific basis (Scientific Advisor):**
The PBE under rotor-stator breakage+coalescence has a unique steady state, but
the relaxation timescale at 10,000 RPM with viscous dispersed phase (μ_d = 1.0 Pa·s,
Alopaeus C3=0.1) exceeds 60 s. Literature (Atiemo-Obeng & Calabrese, 2004) reports
2–10 minutes for batch rotor-stator equilibrium (~500 passes through gap for 500 mL).

**Architecture (Architect):**
The solver has no adaptive termination — `t_emulsification` is a priori fixed.
The convergence check (d32 stability over last 10% of timesteps) is post-hoc
diagnostic only.

**Implementation (Dev-Orchestrator):**

| File | Change | LOC |
|------|--------|-----|
| `configs/default.toml` | `t_emulsification`: 60 → 300 s | 1 |
| `src/emulsim/datatypes.py` | Add to `SolverSettings`: `l1_t_max=600.0`, `l1_conv_tol=0.01`, `l1_max_extensions=2` | 3 |
| `src/emulsim/level1_emulsification/solver.py` | Adaptive extension loop in `solve()`: if not converged, re-init from N_final, extend by t_emul/2, up to l1_max_extensions | ~40 |
| Same file | Same logic for `solve_stirred_vessel()` | ~35 |
| `src/emulsim/datatypes.py` | Add `t_converged: Optional[float]` to `EmulsificationResult` | 1 |
| `tests/` | Convergence tests (default converges, short t triggers extension, volume conservation) | ~25 |

**Total:** ~104 LOC | **Model tier:** Sonnet | **Risk:** MEDIUM (extension may double L1 runtime; cap at 2 extensions)

**Quality gates:**
- G-M1-1: Default config (t=300s) → `converged=True`
- G-M1-2: Short t=5s → adaptive extension triggers, eventually converges
- G-M1-3: `t_history` monotonically increasing after extension
- G-M1-4: Volume conservation within 5% of φ_d
- G-M1-5: Trust warning W1 no longer fires with defaults

---

### M2: Crosslinker Stoichiometry Guidance

**Scientific basis (Scientific Advisor):**
Each genipin bridges 2 NH₂ groups. With defaults: [NH₂] = 18×1000×0.9/161.16 = 100.5 mol/m³,
c_genipin = 2 mol/m³, ratio = 0.020. Maximum achievable p_max = 2×2/100.5 = 0.04 (4%).
Literature (Butler et al., 2003; Muzzarelli, 2009) recommends r = 0.10–0.25 for
structural applications. For target p=0.20: c_genipin = 0.20×100.5/2 ≈ 10 mM.

**Implementation (Dev-Orchestrator):**

| File | Change | LOC |
|------|--------|-----|
| `configs/default.toml` | `c_genipin`: 2.0 → 10.0 mol/m³ | 1 |
| `src/emulsim/trust.py` | Enhance warning: compute and display minimum c_genipin for ratio ≥ 0.05 | 8 |
| `src/emulsim/level3_crosslinking/solver.py` | Add `recommended_crosslinker_concentration(c_chit, DDA, M_GlcN, target_p)` | ~15 |

**Total:** ~24 LOC | **Model tier:** Haiku | **Risk:** LOW

**Quality gates:**
- G-M2-1: Default config ratio ≥ 0.05
- G-M2-2: Trust warning W2 no longer fires with defaults
- G-M2-3: Warning message includes suggested minimum concentration
- G-M2-4: `recommended_crosslinker_concentration(..., target_p=0.3)` → p_final ≈ 0.3

---

### M3: Predictive DN Modulus Model (Phased)

**Scientific basis (Scientific Advisor):**
The formula G_DN = G₁ + G₂ + η√(G₁G₂) has no first-principles derivation.
The √(G₁G₂) coupling is chosen for dimensional convenience, not micromechanics.
It cannot capture IPN swelling constraints, sacrificial-bond toughening, or
topology-dependent mechanical response.

Three model tiers, ordered by physical rigour:
1. **Hashin-Shtrikman composite bounds** — O(1), uses L2 phase fractions
2. **Swelling-constrained IPN model** (Flory-Rehner) — solve 2×2 nonlinear system
3. **Lake-Thomas sacrificial-bond model** — research-grade, strain-dependent

**Architecture (Architect):**
The `select_modulus_model()` dispatch is well-designed. The limitation is within
each branch's model, not the routing. Phase 1 should make warnings mode-aware;
Phase 2 adds Hashin-Shtrikman bounds; Phase 3 adds the Flory-Rehner affine model.

**Implementation (Dev-Orchestrator):**

#### Phase 1: Mode-aware warnings + model metadata
| File | Change | LOC |
|------|--------|-----|
| `src/emulsim/datatypes.py` | Add `model_used: str` and `G_DN_lower/upper: float` to `MechanicalResult` | 4 |
| `src/emulsim/trust.py` | Conditional: suppress W3 in `empirical_engineering`, escalate in `mechanistic_research` | 10 |

#### Phase 2: Hashin-Shtrikman bounds
| File | Change | LOC |
|------|--------|-----|
| `src/emulsim/level4_mechanical/solver.py` | Add `hashin_shtrikman_bounds(G1, G2, phi1)` → (G_lower, G_upper) | ~30 |
| Same file | Populate bounds in `solve_mechanical()` | 10 |

#### Phase 3: Affine IPN model (Flory-Rehner swelling-constrained)
| File | Change | LOC |
|------|--------|-----|
| `src/emulsim/level4_mechanical/solver.py` | Add `flory_rehner_swelling()` — solve coupled equilibrium for Q₁, Q₂ | ~60 |
| Same file | Add `double_network_modulus_affine()` — G_DN = ν_e1·kT·λ₁^(-2/3) + ν_e2·kT·λ₂^(-2/3) | ~80 |
| Same file | Update `select_modulus_model()` for `MECHANISTIC_RESEARCH` mode | ~20 |
| `tests/` | Affine model tests (single-network recovery, symmetric, literature validation) | ~50 |

**Total:** ~264 LOC | **Model tier:** Phase 1–2 Sonnet, Phase 3 Opus | **Risk:** HIGH (Flory-Rehner convergence for extreme compositions)

**Quality gates:**
- G-M3-1: Single network recovery: G_DN(G₁, 0) = G₁ ± 1%
- G-M3-2: Symmetric networks: G_DN(G, G) ≥ G
- G-M3-3: Absolute prediction within 2× of literature for 4% agarose + 2% chitosan DN
- G-M3-4: Fallback to phenomenological when fsolve fails
- G-M3-5: Bounds containment: G_lower ≤ G_DN ≤ G_upper
- G-M3-6: `empirical_engineering` mode still uses phenomenological

---

### M4: Per-Chemistry eta Coupling

**Scientific basis (Scientific Advisor):**
Different crosslinker chemistries produce fundamentally different IPN architectures:

| Crosslinker | Mechanism | Expected η |
|-------------|-----------|------------|
| Genipin | Amine bridge → separate IPN | −0.15 |
| Glutaraldehyde | Schiff base → separate IPN | −0.15 |
| EDC/NHS | Zero-length → tighter IPN | −0.10 |
| ECH | Hydroxyl bridge → bridged networks | +0.05 |
| DVS | Hydroxyl bridge → bridged networks | +0.10 |
| PEGDA+UV | Independent radical network | 0.00 |
| TPP | Ionic, reversible | −0.05 |
| Citric acid | Partial ester bridging | +0.02 |

**Architecture (Architect):**
This is a **wiring defect**: `CrosslinkerProfile.eta_coupling_recommended` already
exists in the reagent library but is never read by the pipeline. The `select_modulus_model()`
partially hardcodes family-specific η (+0.05 for hydroxyl) but doesn't use the
per-crosslinker values.

**Implementation (Dev-Orchestrator):**

| File | Change | LOC |
|------|--------|-----|
| `src/emulsim/datatypes.py` | Add `eta_coupling_recommended: float = -0.15` to `NetworkTypeMetadata` | 2 |
| `src/emulsim/level3_crosslinking/solver.py` | Set `metadata.eta_coupling_recommended = xl.eta_coupling_recommended` at all 4 dispatch points + reaction-diffusion path | 10 |
| `src/emulsim/level4_mechanical/solver.py` | `select_modulus_model()`: read η from `network_metadata.eta_coupling_recommended` | 10 |
| `src/emulsim/trust.py` | Suppress W4 when per-chemistry η is active | 10 |
| `tests/` | Per-chemistry η propagation tests (genipin, DVS, TPP, PEGDA) | ~25 |

**Total:** ~57 LOC | **Model tier:** Sonnet | **Risk:** LOW (backward compatible; genipin default unchanged)

**Quality gates:**
- G-M4-1: Genipin → η = −0.15 (backward compatible)
- G-M4-2: DVS → η = +0.10
- G-M4-3: TPP → η = −0.05
- G-M4-4: PEGDA → η = 0.0
- G-M4-5: Trust warning W4 suppressed when per-chemistry η active
- G-M4-6: G_DN differs between DVS and genipin (different η sign)

---

## Build Order & Dependencies

```
M1 (PBE convergence)     ─── independent ───────────────────────┐
                                                                  │
M4 (per-chemistry eta)   ─── must complete before M3 ──┐         │
                                                        ├── M3   │
M2 (crosslinker ratio)   ─── benefits M3 testing ──────┘  (DN   │
                                                        model)   │
                                                                  │
                              All complete ──────────────────────┘
```

**Recommended sequence:** M1 → M4 → M2 → M3 (Phases 1→2→3)

All three roles converged on this ordering:
- **M1** first: independent, immediate user benefit (default config works)
- **M4** second: pure plumbing, unlocks accurate η for M3 testing
- **M2** third: config + guidance, ensures realistic stoichiometry for M3 validation
- **M3** last: most complex, depends on M4's η values and M2's valid formulations

---

## Risk Matrix

| Risk | Severity | Module | Mitigation |
|------|----------|--------|------------|
| Adaptive extension doubles L1 runtime | MEDIUM | M1 | Cap at `l1_max_extensions=2`, log timing |
| Flory-Rehner solver fails for extreme compositions | HIGH | M3 | Robust fallback to phenomenological + convergence flag |
| Affine model gives unreasonable G_DN | HIGH | M3 | Sanity bounds: 15 kPa < G_DN < 150 kPa for standard formulations |
| Changing default c_genipin breaks comparison runs | LOW | M2 | Old configs with explicit c_genipin unaffected |
| Per-chemistry η changes optimisation landscape | LOW | M4 | Only affects non-genipin crosslinkers; invalidate active campaigns |
| Reaction-diffusion PDE path missing η wiring | MEDIUM | M4 | Must also set `eta_coupling_recommended` in PDE metadata (line ~420) |

---

## Cross-Role Consensus Points

All three roles independently agreed on:

1. **The pipeline architecture is sound** — no structural redesign needed
2. **W1 and W2 are parameter issues**, not model deficiencies
3. **W3 and W4 are deeply connected** — W4 is a specific instance of W3's general limitation
4. **Implementation order:** M1 → M4 → M2 → M3
5. **The reagent library already contains per-chemistry η** — it just needs wiring
6. **The phenomenological model should be kept as a fallback**, not removed
7. **Hashin-Shtrikman bounds** provide a quick scientific upgrade before the full Flory-Rehner model

---

## Literature References

- Atiemo-Obeng & Calabrese (2004). "Rotor-Stator Mixing Devices." *Handbook of Industrial Mixing.*
- Ramkrishna (2000). *Population Balances.* Academic Press.
- Butler et al. (2003). *J. Polym. Sci. A* 41:3941. Genipin-chitosan crosslinking kinetics.
- Muzzarelli (2009). *Carbohydr. Polym.* 77:1. Genipin crosslinking review.
- Mi et al. (2005). *Biomaterials* 26:5983. Genipin 5–20 mM for chitosan.
- Gong (2010). *Soft Matter* 6:2583. DN gel mechanics theory.
- Myung et al. (2008). *Polymer* 49:5575. IPN moduli: η = −0.3 to +0.1.
- Nakajima et al. (2013). *Adv. Funct. Mater.* 22:4426. Sacrificial-bond DN toughening.
- Lake & Thomas (1967). *Proc. R. Soc. Lond. A.* Fracture mechanics of rubber.
- Zhao et al. (2020). *Carbohydr. Polym.* 227:115352. ECH-crosslinked polysaccharide gels.
- Lin & Anseth (2009). *Pharm. Res.* 26:631. PEG hydrogel networks.

---

> **Disclaimer**: This analysis is provided for informational, research, and advisory
> purposes only. All models and parameter recommendations should be validated through
> appropriate laboratory experimentation before implementation. The authors are AI
> assistants; the analysis should be treated as a structured starting point for further
> investigation.
