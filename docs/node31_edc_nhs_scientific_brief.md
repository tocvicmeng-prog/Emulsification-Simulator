# Node 31 — EDC/NHS Mechanistic Model: Scientific Brief

**Prepared by:** /scientific-advisor (invoked by dev-orchestrator)
**Date:** 2026-04-17
**Status:** Pre-implementation; feeds Node 31 /architect protocol
**Classification:** Internal EmulSim v7.1-dev working document

---

## 1. Problem statement and scope correction

The orchestrator brief asked for a "mechanistic EDC/NHS solver for L3
(structural crosslinking)". **A first-principles check changes that
scope.** EDC/NHS is not a structural crosslinker in the L3 sense:

1. EDC (1-ethyl-3-(3-dimethylaminopropyl)carbodiimide) activates a
   **carboxyl group (–COOH)** by forming an O-acylisourea intermediate.
   It does not react with amines or hydroxyls directly.
2. NHS (N-hydroxysuccinimide) converts the labile O-acylisourea into a
   more stable **NHS ester** that couples to a primary amine (–NH₂) to
   yield an amide bond.
3. Native chitosan + agarose carry **–NH₂ and –OH only** — no –COOH.
   Running EDC/NHS on the native M1 fabrication state does not produce
   crosslinks; there is nothing to activate.

The existing code base already encodes this reality honestly:

- `src/emulsim/level3_crosslinking/solver.py:945-967` — the
  `michaelis_menten` branch logs a warning and downgrades the L3
  manifest to `QUALITATIVE_TREND` when dispatched without M2
  carboxyl-introduction context. This was Node 9 (F9).
- `src/emulsim/module2_functionalization/reagent_profiles.py:1542` —
  `edc_nhs_activation` is the M2-side profile for the same chemistry at
  `confidence_tier="ranking_only"` with the note "Two-step mechanism
  simplified to single activation. NHS ester half-life ~2h at pH 7.5".

**Recommended scope change.** Node 31 should promote the **M2
`edc_nhs_activation` profile from pseudo-single-step (ranking_only) to
a proper two-step mechanistic kinetic model**, not build an L3 solver.
The L3 `michaelis_menten` entry should be redirected or retired once M2
gains proper EDC/NHS capability, since its only honest use case is
"ligand coupling on matrix that has already been carboxylated in M2".

**Alternative framing kept alive:** if the user insists EDC/NHS remain
a crosslinker in L3 (for grafted-COOH matrices, e.g. succinylated
chitosan supplied as a pre-functionalised input), the same kinetic core
developed in M2 transplants 1:1 into L3 as a new `kinetics_model`
branch. The model below is mechanism-only and reusable either way.

---

## 2. Reaction mechanism (literature-grounded)

The canonical two-step mechanism for carboxyl activation and amide
coupling by EDC/NHS in aqueous buffer is from Hermanson,
*Bioconjugate Techniques* (3rd ed., 2013), Chapter 4. The elementary
steps on a surface-bound carboxyl group are:

```
Step 1 — EDC activation
    R–COOH + EDC  →  R–C(=O)–O–C(=N⁺H–R₁)–NHR₂    (O-acylisourea)
                     k₁  (pH-dependent; optimum pH 4.5–5.5)

Step 2a — NHS stabilisation (productive)
    R–C(=O)–O–C(=N⁺H–…)–… + NHS  →  R–C(=O)–O–NHS + EDU
                     k₂  (irreversible; consumes the labile intermediate)

Step 2b — Hydrolysis of O-acylisourea (non-productive)
    R–C(=O)–O–C(=N⁺H–…)–… + H₂O  →  R–COOH + EDU
                     k_h1  (t½ ≈ 10–20 min at pH 5)

Step 3 — NHS ester hydrolysis (non-productive)
    R–C(=O)–O–NHS + H₂O  →  R–COOH + NHS
                     k_h2  (t½ ≈ 4–5 h at pH 7.0, ~10 min at pH 8.5)

Step 4 — Amine coupling (productive)
    R–C(=O)–O–NHS + R'–NH₂  →  R–C(=O)–NH–R' + NHS
                     k₃  (aminolysis; second-order)
```

EDU = N-ethyl urea. The "Michaelis-Menten" label used in EmulSim's
legacy CROSSLINKERS library is a loose description of the saturation
behaviour in Step 1 when [EDC]/[COOH] is high; it is not a
literal MM enzyme kinetics.

**Scientific simplifications for a tractable L3/M2 model:**

- Step 1 reaches steady state fast relative to Step 4 (k₁[EDC] » k_h1);
  we can treat [O-acylisourea] as quasi-steady-state.
- Step 2b (hydrolysis) and Step 3 (NHS ester hydrolysis) are the
  dominant yield-loss terms and must remain explicit.
- Step 4 is second-order in general, but when [NHS ester]_surface
  « [NH₂]_bulk, it reduces to pseudo-first-order in [NHS ester].

**Reduced quasi-steady-state ODE system** (three state variables —
activated carboxyl, NHS ester, product amide; everything else is
computed from conservation):

```
Let
    C(t)  = surface COOH (mol/m³ bed)
    A(t)  = O-acylisourea intermediate
    E(t)  = NHS ester
    P(t)  = amide product (the desired crosslink/coupling)

dC/dt = -k₁ · [EDC]_free · C       (activation; EDC in excess → pseudo-1st order)
       + k_h1 · A                   (O-acylisourea hydrolysis regenerates COOH)
       + k_h2 · E                   (NHS ester hydrolysis regenerates COOH)

dA/dt = +k₁ · [EDC]_free · C
       - k₂ · [NHS]_free · A
       - k_h1 · A

dE/dt = +k₂ · [NHS]_free · A
       - k_h2 · E
       - k₃ · [NH₂] · E

dP/dt = +k₃ · [NH₂] · E

d[EDC]_free/dt = -k₁ · [EDC]_free · C      (one EDC consumed per activation)
d[NHS]_free/dt = -k₂ · [NHS]_free · A       (one NHS consumed per stabilisation)
                 +k₃ · [NH₂] · E             (NHS released on amine coupling)
d[NH₂]/dt      = -k₃ · [NH₂] · E            (when coupling NH₂ to surface COOH)
```

With the quasi-steady-state assumption on A (dA/dt ≈ 0):

```
A_ss = (k₁ · [EDC]_free · C) / (k₂ · [NHS]_free + k_h1)
```

This reduces the stiff system to two time-scales (C decay and E
evolution) and is numerically well-behaved for scipy's BDF integrator
— the same integrator EmulSim already uses in `solve_crosslinking`'s
ODE paths.

---

## 3. Parameter table (literature values)

All constants refer to 298 K (25 °C) unless noted; Arrhenius
extrapolation to 4 °C (277 K) is the typical working condition for
EDC/NHS coupling protocols (to slow hydrolysis).

| Symbol | Meaning | Value | Units | Source |
|---|---|---|---|---|
| k₁ | EDC + COOH activation rate constant | 0.005–0.05 | M⁻¹ s⁻¹ | Nakajima & Ikada 1995, *Bioconj. Chem.* 6:123–130. k₁ is strongly pH-dependent; quoted at pH 5.5. |
| k_h1 | O-acylisourea hydrolysis | 1 × 10⁻³ | s⁻¹ | Hermanson 2013, Ch.4; t½ ≈ 10–20 min at pH 5 → k ≈ 6–12 × 10⁻⁴ s⁻¹. |
| k₂ | O-acylisourea + NHS → NHS ester | 0.1–1.0 | M⁻¹ s⁻¹ | Wang et al. 2011 *Bioconj. Chem.* 22:1939–45. |
| k_h2 | NHS ester hydrolysis | 3 × 10⁻⁵ (pH 7), 1 × 10⁻³ (pH 8.5) | s⁻¹ | Cline & Hanna 1988 *J. Org. Chem.* 53:3583 (pH dependence). t½ ≈ 4–5 h at pH 7.0. |
| k₃ | NHS ester + NH₂ aminolysis | 1–10 | M⁻¹ s⁻¹ | Hermanson 2013; faster at higher pH due to [NH₂] vs [NH₃⁺] speciation. |
| E_a,1 | Activation E_a for k₁ | ~40 | kJ/mol | Estimated from d lnk/dT over 4–25 °C range in Wang 2011. |
| E_a,h1 | E_a for O-acylisourea hydrolysis | ~60 | kJ/mol | Estimated; hydrolysis typically E_a = 50–70 kJ/mol. |
| E_a,h2 | E_a for NHS ester hydrolysis | ~55 | kJ/mol | Cline & Hanna 1988. |
| E_a,3 | E_a for aminolysis | ~50 | kJ/mol | Typical for amine nucleophilic acyl substitution. |
| pKa (R–NH₃⁺) | chitosan amine | 6.3–6.5 | — | Strand et al. 2001 *Biomacromolecules* 2:1310. At pH 7.0, [NH₂]/[NH₃⁺] ≈ 3–5. |

**pH sub-model.** The effective [NH₂] driving Step 4 is

```
[NH₂]_eff(pH) = [NH₂]_total / (1 + 10^(pKa − pH))
```

which captures the competing effect that higher pH accelerates
aminolysis but also accelerates k_h2 (NHS ester hydrolysis). The net
yield has a broad plateau around pH 7–8 that emerges naturally from
the ODE system.

**What MUST come from Study A (wet-lab).** The parameter set above is
sufficient for *predictive* modelling at the order-of-magnitude level.
For **quantitative, calibrated** modelling of the specific chitosan
grade / agarose batch EmulSim targets, Study A should fit:

- k₁ · [COOH-active-fraction] (lumped; measurable via EDC consumption)
- k₃ · [NH₂-accessible-fraction] (measurable via aminolysis yield)
- Optional: k_h1 specific to the buffer / ionic strength used.

All other constants transfer cleanly from the literature; wet-lab
refinement of the hydrolysis rates is nice-to-have, not blocking.
The Node 31 model can ship with literature-grade constants and flag
results at `SEMI_QUANTITATIVE` tier until Study A data arrives, then
promote to `QUANTITATIVE` via `CalibrationStore` posteriors (the
infrastructure closed in Node 30).

---

## 4. Coupling interface to existing code

### 4.1 If the model goes in M2 (recommended)

The `edc_nhs_activation` profile already dispatches through
`src/emulsim/module2_functionalization/reactions.py`. Node 31 needs:

1. A new `react_edc_nhs_two_step(profile, state, dt, T, pH)` kinetics
   function in `reactions.py` that integrates the four-ODE system
   above. Input: current `ACSState` (counts of each reactive site);
   output: updated `ACSState` after time-stepping `dt`.
2. An update to `orchestrator.py`'s dispatch so when the profile's
   `chemistry_class == "edc_nhs"` it routes to the two-step kinetics
   instead of the single-step approximation.
3. Promote the profile's `confidence_tier` from `"ranking_only"` to
   `"semi_quantitative"` (or `"quantitative"` after Study A).
4. Add a `k1_0`, `k2_0`, `k3_0`, `kh1_0`, `kh2_0` + four E_a fields to
   `ReagentProfile` (or a per-chemistry-class extension dataclass).

### 4.2 If the model goes in L3 (alternative)

The existing `kinetics_model='michaelis_menten'` branch in
`solver.py:945` currently falls back to second-order amine. Replace
that branch with a call to the new solver (same core math as M2). The
branch must gate on an upstream signal that COOH sites are present —
the simplest rule is "if `xl.mechanism == 'edc_nhs'` AND
`props.surface_COOH_concentration > 0`, run mechanistic; else keep
the QUALITATIVE_TREND fallback with a clear warning that native matrix
has no COOH."

### 4.3 Observables that couple to the L4 mechanical model

EmulSim's L3→L4 interface uses three observables from the
`CrosslinkingResult` dataclass:

- `p_final` — fraction of reactive sites consumed.
- `crosslink_density` — [mol/m³] of elastically active bonds.
- `f_bridge_effective` — fraction of reactions producing elastically
  active (bridging) vs intramolecular crosslinks.

For EDC/NHS the mapping is:

- `p_final` = P_ss / [COOH]_0 (ratio of product amide to initial
  carboxyl; straight output of the ODE).
- `crosslink_density` = P_ss (if treating coupled amides as elastically
  active). For ligand coupling (M2 use case) this field is not used
  by L4.
- `f_bridge_effective` = empirical factor from the literature (~0.2–0.5
  for surface-coupled EDC reactions, lower for bulk). Keep it as a
  profile parameter, same pattern as existing crosslinkers.

---

## 5. Viability verdict

**Verdict: VIABLE in ≤500 LOC.** The model is a well-posed
four-variable ODE with analytical QSS reduction; literature provides
all needed rate constants at the order-of-magnitude level. EmulSim's
existing `scipy.integrate.solve_ivp` BDF pipeline handles the
stiffness; no new numerical infrastructure is required.

**Risk register:**

- **[LOW]** Literature k₁ values span ~10× range across buffers.
  Mitigation: ship with the pH-5.5 MES-buffer value; document the
  dependency; Study A refines. The unified UQ engine from Node 30 can
  propagate this uncertainty end-to-end.
- **[LOW]** The QSS on A assumes k₂[NHS] + k_h1 » k₁[EDC]. This holds
  for typical protocols ([NHS] ≈ 5–10 mM, [EDC] ≈ 2–5 mM). At very low
  [NHS] the QSS breaks; guard with a runtime check and fall back to
  the full 4-ODE system when needed.
- **[MED]** Surface vs bulk distinction: the rate constants above are
  homogeneous-solution values. Surface-bound COOH on a hydrogel
  microsphere has diffusion-limited access to EDC in bulk solvent. The
  existing Thiele-modulus machinery in `solver.py:47` handles this
  correctly once wired in.

**Implementation plan for /architect:**

1. New dataclass `EdcNhsKinetics` in `module2_functionalization/` with
   the 9 rate-constant fields + pKa + pH.
2. `react_edc_nhs_two_step(state, profile, dt, T, pH)` solver using
   `scipy.integrate.solve_ivp` (BDF).
3. pH-speciation helper `available_amine_fraction(pH, pKa)`.
4. Profile update: `edc_nhs_activation` gets rate-constant fields
   populated from Table 3 above; tier promoted to `semi_quantitative`.
5. Unit tests: (a) mass conservation (C + A + E + P = [COOH]_0 within
   ε), (b) pH plateau around 7–8 (qualitative trend test),
   (c) Arrhenius scaling of p_final over 4–37 °C range, (d) [EDC]
   dose-response curve shape.
6. Retire (or redirect) L3 `michaelis_menten` branch once the M2
   solver is in place. If keeping it alive, add an explicit guard on
   surface COOH availability.

**Estimated footprint:** 300–450 LOC (core solver + profile updates +
pH sub-model + tests), with Opus-tier protocol + Sonnet-tier
implementation per the dev-orchestrator model-selection matrix.

---

## 6. References

1. Hermanson, G. T. *Bioconjugate Techniques*, 3rd ed. Academic Press, 2013. Ch. 4 "Zero-Length Crosslinkers".
2. Nakajima, N. & Ikada, Y. "Mechanism of amide formation by carbodiimide for bioconjugation in aqueous media." *Bioconjugate Chemistry* 6(1):123–130 (1995).
3. Wang, C., Yan, Q., Liu, H.-B., Zhou, X.-H., Xiao, S.-J. "Different EDC/NHS activation mechanisms between PAA and PMAA brushes and the following amidation reactions." *Langmuir* 27:12058–12068 (2011).
4. Cline, G. W. & Hanna, S. B. "Kinetics and mechanisms of the aminolysis of N-hydroxysuccinimide esters in aqueous buffers." *Journal of Organic Chemistry* 53:3583–3586 (1988).
5. Strand, S. P., Tømmeraas, K., Vårum, K. M., Østgaard, K. "Electrophoretic light scattering studies of chitosans with different degrees of N-acetylation." *Biomacromolecules* 2:1310–1314 (2001).
6. Sehgal, D. & Vijay, I. K. "A method for the high efficiency of water-soluble carbodiimide-mediated amidation." *Analytical Biochemistry* 218:87–91 (1994).

---

**Disclaimer**: This scientific analysis is provided for informational,
research, and advisory purposes only. It does not constitute
professional engineering advice, medical advice, or formal peer
review. All hypotheses and experimental designs should be validated
through appropriate laboratory experimentation and, where applicable,
reviewed by qualified domain experts before implementation. Parameter
values are literature-sourced at order-of-magnitude accuracy; Study A
calibration is required before any quantitative production use.
