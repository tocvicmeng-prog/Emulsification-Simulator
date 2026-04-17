# Node 31 — EDC/NHS Mechanistic Kinetic Model: /architect Protocol

**Prepared by:** /architect (via dev-orchestrator)
**Feeds:** /scientific-coder for Phase 2 implementation
**Date:** 2026-04-17
**Status:** Ready for Phase 2
**Scientific basis:** `docs/node31_edc_nhs_scientific_brief.md`

---

## 1. Purpose

Implement a two-step mechanistic EDC/NHS carbodiimide coupling kinetic
in M2 and expose it in L3 behind a surface-COOH availability gate.
Promote `edc_nhs_activation` from pseudo-single-step (ranking_only) to
a proper SEMI_QUANTITATIVE kinetic model. Close the Node 9 (F9) debt
item that forced L3 EDC/NHS to QUALITATIVE_TREND.

## 2. Interface specification

### 2.1 New dataclass `EdcNhsKinetics`

File: `src/emulsim/module2_functionalization/edc_nhs_kinetics.py` (new)

```python
@dataclass(frozen=True)
class EdcNhsKinetics:
    """Rate constants + activation energies for two-step EDC/NHS.

    All rate constants are at the reference temperature T_ref = 298.15 K;
    Arrhenius-extrapolated at runtime. Values default to literature
    medians from Hermanson 2013 / Wang 2011 / Cline & Hanna 1988.
    """
    # Step 1: EDC activation of COOH (pseudo-1st order when [EDC] >> [COOH])
    k1_0: float = 0.02             # [M^-1 s^-1] at 298 K, pH 5.5
    E_a_1: float = 40e3            # [J/mol]

    # Step 2: O-acylisourea + NHS -> NHS ester (productive stabilisation)
    k2_0: float = 0.5              # [M^-1 s^-1] at 298 K
    E_a_2: float = 50e3            # [J/mol]

    # Step 2b: O-acylisourea hydrolysis (non-productive)
    kh1_0: float = 1e-3            # [s^-1] at 298 K, pH 5.5
    E_a_h1: float = 60e3           # [J/mol]

    # Step 3: NHS ester hydrolysis (non-productive, pH-dependent)
    kh2_0: float = 3e-5            # [s^-1] at 298 K, pH 7.0
    E_a_h2: float = 55e3           # [J/mol]
    kh2_pH_slope: float = 1.5      # kh2 multiplies 10^(kh2_pH_slope * (pH - 7))

    # Step 4: NHS ester + NH2 -> amide (productive aminolysis)
    k3_0: float = 5.0              # [M^-1 s^-1] at 298 K (with [NH2] speciated)
    E_a_3: float = 50e3            # [J/mol]

    # Chitosan amine speciation
    pKa_amine: float = 6.4         # Strand et al. 2001 for chitosan
```

### 2.2 New solver `react_edc_nhs_two_step`

File: `src/emulsim/module2_functionalization/edc_nhs_kinetics.py`

```python
def react_edc_nhs_two_step(
    *,
    c_cooh_initial: float,        # [mol/m^3] initial COOH on matrix
    c_nh2_total: float,           # [mol/m^3] total NH2 sites
    c_edc_initial: float,         # [mol/m^3] EDC in solution
    c_nhs_initial: float,         # [mol/m^3] NHS in solution
    pH: float,
    T: float,                     # [K]
    time: float,                  # [s] reaction time
    kin: EdcNhsKinetics = EdcNhsKinetics(),
    rtol: float = 1e-6,
    atol: float = 1e-10,
) -> EdcNhsResult:
    """Integrate the four-variable ODE system from the scientific brief §2.

    Uses scipy.integrate.solve_ivp with the BDF method (same choice as
    the existing L3 ODE paths). When 4*k2*[NHS] + kh1 >> k1*[EDC] the
    o-acylisourea (A) reaches QSS in <1 s; the integrator handles this
    transparently but we optionally log a diagnostic.

    Returns
    -------
    EdcNhsResult with:
        p_final         — fraction of initial COOH converted to amide
        p_hydrolysed    — fraction lost to hydrolysis (C + acyl paths)
        p_residual_nhs_ester — NHS ester still present at end of run
        time_to_half    — [s] time to reach 0.5 * p_final (yields a
                          kinetic diagnostic for the trust banner)
        ode_diagnostics — dict with max_stiffness_ratio, n_steps, etc.
    """
```

### 2.3 Result dataclass `EdcNhsResult`

Same file. Frozen dataclass with the 5 fields above.

### 2.4 pH speciation helper

```python
def available_amine_fraction(pH: float, pKa: float = 6.4) -> float:
    """Deprotonated fraction of primary amines at given pH.

    [NH2] / [NH2 + NH3+] = 1 / (1 + 10^(pKa - pH))
    """
```

### 2.5 Arrhenius helper

Reuse the existing `arrhenius_rate_constant` in
`src/emulsim/level3_crosslinking/solver.py:72`. Do NOT duplicate.

## 3. Dispatch wiring

### 3.1 M2 side

File: `src/emulsim/module2_functionalization/reactions.py`

- Add a new branch in the step-kinetic dispatch for
  `profile.chemistry_class == "edc_nhs"` that calls
  `react_edc_nhs_two_step` and updates the `ACSState` accordingly:
  - `target_acs = CARBOXYL_DISTAL` → decremented by `p_final * c_cooh_initial`
  - `product_acs = AMIDE_BOND` (new enum value if not already present —
    verify in `acs.py`; otherwise reuse `NHS_ESTER` → `AMIDE` transition)
- Promote `edc_nhs_activation` profile's `confidence_tier` from
  `"ranking_only"` to `"semi_quantitative"` in `reagent_profiles.py:1563`.

### 3.2 L3 side (gated)

File: `src/emulsim/level3_crosslinking/solver.py`

Replace the `michaelis_menten` branch (lines 945–1000) with a gated
dispatch:

```python
elif xl.kinetics_model == 'michaelis_menten':
    # EDC/NHS requires COOH sites. Gate on surface COOH availability.
    surface_cooh = getattr(props, "surface_cooh_concentration", 0.0) or 0.0
    if surface_cooh <= 0:
        # Native matrix has no COOH — honest fallback, QUALITATIVE_TREND
        # (same behaviour as v7.0.1; preserved for safety).
        _fallback_used = True
        # ... existing fallback code ...
    else:
        # Run mechanistic model (M2 logic transplanted).
        from ..module2_functionalization.edc_nhs_kinetics import (
            react_edc_nhs_two_step, EdcNhsKinetics,
        )
        res = react_edc_nhs_two_step(
            c_cooh_initial=surface_cooh,
            c_nh2_total=available_amine_concentration(
                params.formulation.c_chitosan, props.DDA, props.M_GlcN,
            ),
            c_edc_initial=params.formulation.c_genipin,  # reused field
            c_nhs_initial=params.formulation.c_genipin,  # reused field
            pH=getattr(params.formulation, "pH", 7.0),
            T=params.formulation.T_crosslink,
            time=params.formulation.t_crosslink,
        )
        result = _build_crosslinking_result_from_edc_nhs(res, ...)
        metadata = NetworkTypeMetadata(
            solver_family="amine_covalent",
            network_target="edc_activated",
            bond_type="covalent",
            is_true_second_network=True,
            eta_coupling_recommended=xl.eta_coupling_recommended,
        )
        # Tier: SEMI_QUANTITATIVE (literature constants)
        _l3_tier = ModelEvidenceTier.SEMI_QUANTITATIVE
        _fallback_used = False
```

The field `surface_cooh_concentration` is a new optional attribute on
`MaterialProperties`. Initialise to `0.0` by default; populated by M2
when a prior succinylation / carboxyl-graft step has run.

### 3.3 New MaterialProperties field

File: `src/emulsim/datatypes.py` — add to the `MaterialProperties`
dataclass:

```python
surface_cooh_concentration: float = 0.0  # [mol/m^3] COOH sites on matrix
"""Surface carboxyl concentration introduced by prior M2 steps (e.g.
succinylation). Zero on native chitosan/agarose; non-zero gates L3
EDC/NHS to run the mechanistic path rather than fall back to
QUALITATIVE_TREND."""
```

Update `CalibrationStore.apply_to_model_params` documentation noting
this is L2/L3 target.

## 4. Algorithm (reference — see scientific brief §2 for full derivation)

Integrate four coupled ODEs from t=0 to t=time:

```
dC/dt  = -k1 * [EDC]_free * C + kh1 * A + kh2_eff * E
dA/dt  = +k1 * [EDC]_free * C - k2 * [NHS]_free * A - kh1 * A
dE/dt  = +k2 * [NHS]_free * A - kh2_eff * E - k3 * [NH2]_eff * E
dP/dt  = +k3 * [NH2]_eff * E
```

Coupled consumption balances:

```
d[EDC]_free/dt = -k1 * [EDC]_free * C
d[NHS]_free/dt = -k2 * [NHS]_free * A + k3 * [NH2]_eff * E
d[NH2]_total/dt = -k3 * [NH2]_eff * E
```

(where `[NH2]_eff = [NH2]_total * available_amine_fraction(pH, pKa)`)

Rate-constant evaluation:

```
k1 = k1_0 * exp(-E_a_1 / R * (1/T - 1/T_ref))
k2 = k2_0 * exp(-E_a_2 / R * (1/T - 1/T_ref))
kh1 = kh1_0 * exp(-E_a_h1 / R * (1/T - 1/T_ref))
kh2_eff = kh2_0 * exp(-E_a_h2 / R * (1/T - 1/T_ref))
           * 10^(kh2_pH_slope * (pH - 7))
k3 = k3_0 * exp(-E_a_3 / R * (1/T - 1/T_ref))
```

with `R = 8.314` J/(mol·K), `T_ref = 298.15` K.

Initial conditions: `C(0) = c_cooh_initial`, `A(0) = E(0) = P(0) = 0`,
`[EDC]_free(0) = c_edc_initial`, `[NHS]_free(0) = c_nhs_initial`,
`[NH2]_total(0) = c_nh2_total`.

Complexity: O(n_steps) where n_steps scales with the stiffness ratio.
Typical: 50–200 BDF steps for 1-hour reaction time.

## 5. Error handling

| Condition | Detection | Response |
|---|---|---|
| `c_cooh_initial <= 0` | Input validation | Raise `ValueError` with explicit "EDC/NHS requires surface COOH; set surface_cooh_concentration > 0 on MaterialProperties" |
| `T <= 0` or `not finite` | Input validation | Raise `ValueError` |
| `pH < 3 or pH > 10` | Input validation | Raise `ValueError` (model validity domain) |
| Integrator fails (returns `success=False`) | Check `sol.success` | Log warning, return `EdcNhsResult` with `p_final=nan` and diagnostic |
| Mass-conservation check `abs(C + A + E + P - c_cooh_initial) > 1% * c_cooh_initial` | Post-run verification | Log warning (ODE tolerance issue); do not fail |

## 6. Test cases

File: `tests/test_edc_nhs_kinetics.py` (new, ~150 LOC, ~8–10 tests)

### Core correctness
1. **Mass conservation**: after any integration, `C + A + E + P` equals
   `c_cooh_initial` within `1e-3 * c_cooh_initial`.
2. **Zero-EDC: no activation**: `c_edc_initial=0` → `p_final=0` exactly.
3. **Zero-NH2: no coupling** (but NHS ester can still form): `c_nh2_total=0` → `p_final=0`, `p_residual_nhs_ester > 0`.
4. **Zero-NHS: hydrolysis dominates**: `c_nhs_initial=0` → `p_final << p_hydrolysed`.

### Scientific trends
5. **pH plateau**: `p_final` at pH 6.5, 7.0, 7.5, 8.0 all within a 2×
   window; pH 5.0 and pH 9.0 drop noticeably (hydrolysis at high pH,
   low [NH2]_eff at low pH).
6. **Arrhenius scaling**: `p_final(T=310K) > p_final(T=298K) >
   p_final(T=277K)` for the same time budget (hydrolysis and coupling
   both accelerate but coupling wins at moderate T).
7. **Dose response**: doubling `c_edc_initial` with everything else
   fixed produces `p_final` closer to its asymptotic limit.

### Infrastructure
8. **pH speciation helper**: `available_amine_fraction(pKa, pKa) ==
   0.5` exactly; monotonic in pH.
9. **Input validation**: all the error cases above raise as declared.

### Integration with M2
10. **M2 profile uses two-step**: a smoke test that drives a
    ModificationStep with the `edc_nhs_activation` profile and
    confirms the kinetic path is exercised (by patching
    `react_edc_nhs_two_step` and checking it's called).

### Integration with L3 gate
11. **L3 native matrix: fallback**: `props.surface_cooh_concentration=0` →
    L3 still hits the QUALITATIVE_TREND fallback (no regression).
12. **L3 carboxylated matrix: mechanistic**: `surface_cooh_concentration=50`
    → L3 uses mechanistic path, manifest tier is SEMI_QUANTITATIVE.

## 7. Performance budget

| Metric | Target |
|---|---|
| `react_edc_nhs_two_step` single call (1-hour reaction) | ≤ 50 ms |
| Mass-balance error | ≤ 1e-3 relative |
| BDF integration steps | 50–500 typical |
| Test-suite addition (edc_nhs tests) | ≤ 5 s |
| LOC added | 300–450 (solver + profile + L3 gate + tests) |

## 8. Dependencies

- Upstream: scipy.integrate.solve_ivp (already in use; BDF method)
- Upstream: `arrhenius_rate_constant` in L3 solver (reuse)
- Downstream: M2 dispatch in `reactions.py`, L3 dispatch in
  `solver.py`, `MaterialProperties.surface_cooh_concentration` field
- No new external libraries

## 9. Logging

- `logger.info("EDC/NHS kinetic: T=%.1fK, pH=%.2f, t=%.0fs → p_final=%.3f", ...)` at run start/end
- `logger.debug("QSS diagnostic: k2*[NHS] = %.2e, kh1 = %.2e, k1*[EDC] = %.2e", ...)` for reproducibility
- `logger.warning("Mass conservation: C+A+E+P = %.6f, expected %.6f (rel err %.2%)", ...)` if balance slips

## 10. Risks + migration notes

1. **pH field on FormulationParameters may not exist.** Check
   `datatypes.py` and add `pH: float = 7.0` if absent.
2. **NHS and EDC currently share `c_genipin` as their concentration
   field.** This is a workaround (reuse of the crosslinker concentration
   slot). Cleaner would be a dedicated `FormulationParameters.c_edc`
   + `c_nhs`; keep this scope-creep for Node 31b.
3. **CrosslinkingResult.network_target="edc_activated"** is a new tag
   value. Verify downstream UI panels tolerate unknown tags gracefully.

## 11. Gate G1 (12-point completeness check)

| # | Criterion | Status |
|---|---|---|
| G1-01 | Purpose statement is one sentence, unambiguous | ✓ §1 |
| G1-02 | Every input parameter typed, shape, units, constraints, default | ✓ §2.2 |
| G1-03 | Every output parameter typed, shape, units, guarantees | ✓ §2.3 |
| G1-04 | Algorithm has plain-language + pseudocode | ✓ §4 (full derivation in brief §2) |
| G1-05 | Time and space complexity stated with justification | ✓ §4 |
| G1-06 | Numerical considerations (stability, precision) documented | ✓ §4 (stiffness, BDF, mass-balance check) |
| G1-07 | ≥3 unit test cases with concrete input → expected output | ✓ §6 (tests 1–4) |
| G1-08 | ≥1 boundary test per input constraint | ✓ §6 (tests 2–3, 8–9) |
| G1-09 | Error conditions with detection + response | ✓ §5 |
| G1-10 | Performance budget quantified | ✓ §7 |
| G1-11 | Upstream + downstream modules named | ✓ §8, §3 |
| G1-12 | Logging level + metrics specified | ✓ §9 |

**G1: PASSED.** Ready for Phase 2 (implementation, Sonnet tier —
standard scientific ODE integration; no novel numerical methods).

---

**Next session:** /scientific-coder reads `docs/node31_edc_nhs_scientific_brief.md` + this protocol, implements the 5 files per §2–§3, runs the 12 tests per §6, and hands off. Estimate: ~350–500 LOC, 1 fix round likely.
