"""Mechanistic two-step EDC/NHS kinetic solver (Node 31, v7.1).

Integrates the four-variable ODE system derived in
``docs/node31_edc_nhs_scientific_brief.md`` §2. The state variables are

    C(t) — surface carboxyl (COOH)           [mol/m^3]
    A(t) — O-acylisourea intermediate        [mol/m^3]
    E(t) — NHS ester                         [mol/m^3]
    P(t) — product amide                     [mol/m^3]

and three coupled free-reagent / free-amine balances:

    [EDC]_free(t) — solution EDC             [mol/m^3]
    [NHS]_free(t) — solution NHS             [mol/m^3]
    [NH2]_total(t) — total surface amine     [mol/m^3]

Mass conservation at every t: C + A + E + P == c_cooh_initial (within
the integrator tolerance). Mass-balance slip is detected and reported
via ``EdcNhsResult.mass_balance_error``; values > 1 % trigger a WARNING
log but the result is still returned so callers can decide.

The mechanism reduces to a well-behaved stiff ODE handled by
scipy.integrate.solve_ivp with BDF. For the M2 activation use case the
caller passes ``c_nh2_total=0`` so step 4 (aminolysis) produces zero
amide; the useful output is the final NHS ester concentration. For the
L3 use case (COOH-grafted matrix) the caller passes
``c_nh2_total > 0`` and the useful output is ``p_final`` = P/c_cooh_initial.

See also: ``docs/node31_protocol.md`` for the architect's module
protocol and Gate G1 checklist.
"""

from __future__ import annotations

import logging
import math
from dataclasses import dataclass

import numpy as np
from scipy.integrate import solve_ivp

logger = logging.getLogger(__name__)


R_GAS = 8.314  # [J/(mol*K)]
T_REF = 298.15  # [K] reference temperature for all k_0 values


# ─── Parameter bundle (literature defaults from scientific brief §3) ────


@dataclass(frozen=True)
class EdcNhsKinetics:
    """Rate constants and activation energies for two-step EDC/NHS.

    Defaults are literature medians at T_ref = 298.15 K:

      - Hermanson 2013 *Bioconjugate Techniques* ch.4 (mechanism + k_h1, k_h2)
      - Nakajima & Ikada 1995 *Bioconj. Chem.* 6:123 (k_1)
      - Wang et al. 2011 *Langmuir* 27:12058 (k_2, E_a scaling)
      - Cline & Hanna 1988 *J. Org. Chem.* 53:3583 (k_h2 pH dependence)
      - Strand et al. 2001 *Biomacromolecules* 2:1310 (chitosan pKa)
    """
    # Step 1: EDC activation of COOH
    k1_0: float = 0.02             # [m^3 / (mol*s)] at 298 K, pH 5.5
    E_a_1: float = 40e3            # [J/mol]

    # Step 2: O-acylisourea + NHS -> NHS ester (productive stabilisation)
    k2_0: float = 0.5              # [m^3 / (mol*s)]
    E_a_2: float = 50e3            # [J/mol]

    # Step 2b: O-acylisourea hydrolysis (non-productive)
    kh1_0: float = 1e-3            # [1/s] at 298 K, pH 5.5
    E_a_h1: float = 60e3           # [J/mol]

    # Step 3: NHS ester hydrolysis (non-productive, pH-dependent)
    kh2_0: float = 3e-5            # [1/s] at 298 K, pH 7.0
    E_a_h2: float = 55e3           # [J/mol]
    kh2_pH_slope: float = 1.5
    """k_h2_effective = k_h2_0 * exp(-E_a_h2/R * (1/T - 1/T_ref)) *
    10^(kh2_pH_slope * (pH - 7)). Slope 1.5/pH-unit matches Cline &
    Hanna 1988."""

    # Step 4: NHS ester + NH2 -> amide (productive aminolysis)
    k3_0: float = 5.0              # [m^3 / (mol*s)]
    E_a_3: float = 50e3            # [J/mol]

    # Chitosan amine speciation
    pKa_amine: float = 6.4


# ─── Result dataclass ─────────────────────────────────────────────────────


@dataclass(frozen=True)
class EdcNhsResult:
    """Outcome of a single EDC/NHS integration run.

    Fractions are with respect to ``c_cooh_initial`` (the input COOH).
    """
    p_final: float
    """Fraction of initial COOH converted to amide bond (P/C0)."""

    p_hydrolysed: float
    """Fraction of initial COOH regenerated via hydrolysis pathways
    (O-acylisourea + NHS-ester hydrolysis; includes COOH that was
    re-activated and re-hydrolysed — this is 'lost to hydrolysis cycle'
    rather than committed amide)."""

    p_residual_nhs_ester: float
    """NHS ester still present at the end of the integration (E/C0).
    For M2 activation use case this is the active yield ready for
    subsequent amine coupling."""

    time_to_half: float
    """[s] Time to reach 0.5 * p_final. NaN when p_final is near zero."""

    mass_balance_error: float
    """|C + A + E + P - C0| / C0 at t=t_end. Target <= 1e-3."""

    solver_success: bool
    solver_message: str

    n_steps: int
    """Integrator step count (ode diagnostic)."""


# ─── pH speciation helper (thin wrapper over existing reactions.py) ───────


def available_amine_fraction(pH: float, pKa: float = 6.4) -> float:
    """Deprotonated (reactive) fraction of primary amines at given pH.

    [NH2] / ([NH2] + [NH3+]) = 1 / (1 + 10^(pKa - pH))

    Unlike the generic ``ph_rate_scaling`` in reactions.py (which uses a
    Hill coefficient and a disabled-if-pKa<=0 convention), this helper
    is a pure speciation function hard-wired to n_hill=1 for clarity at
    the EDC/NHS call site.
    """
    exponent = pKa - pH
    if exponent > 300:
        return 0.0
    if exponent < -300:
        return 1.0
    return 1.0 / (1.0 + 10.0 ** exponent)


# ─── Arrhenius helper (local copy to avoid cross-package import) ──────────


def _arrhenius(T: float, k0: float, E_a: float) -> float:
    """k(T) = k0 * exp(-E_a/R * (1/T - 1/T_ref)). Returns in same units as k0."""
    if T <= 0 or not math.isfinite(T):
        raise ValueError(f"T must be finite and positive, got {T}")
    return k0 * math.exp(-E_a / R_GAS * (1.0 / T - 1.0 / T_REF))


# ─── The solver ───────────────────────────────────────────────────────────


def react_edc_nhs_two_step(
    *,
    c_cooh_initial: float,
    c_nh2_total: float,
    c_edc_initial: float,
    c_nhs_initial: float,
    pH: float,
    T: float,
    time: float,
    kin: EdcNhsKinetics | None = None,
    rtol: float = 1e-6,
    atol: float = 1e-10,
) -> EdcNhsResult:
    """Integrate the four-ODE EDC/NHS system from t=0 to t=time.

    Parameters
    ----------
    c_cooh_initial : float
        [mol/m^3] initial COOH surface concentration. Must be > 0.
    c_nh2_total : float
        [mol/m^3] total surface amine concentration. Set to 0 for the
        M2 activation-only use case (no aminolysis occurs; NHS ester
        simply accumulates and hydrolyses).
    c_edc_initial, c_nhs_initial : float
        [mol/m^3] solution concentrations of EDC and NHS. Typical
        protocol: [EDC] = 2–10 mM (= 2–10 mol/m^3), [NHS] = 5–10 mM.
    pH : float
        Reaction pH. Valid domain 3.0–10.0; raises outside.
    T : float
        [K] reaction temperature. Typical: 277 K (4 °C) or 298 K (25 °C).
    time : float
        [s] total reaction time.
    kin : EdcNhsKinetics, optional
        Rate-constant bundle. Defaults to literature medians.
    rtol, atol : float
        scipy.integrate.solve_ivp tolerances. Defaults are conservative
        for the mass-balance invariant; looser tolerances may produce
        mass_balance_error > 1e-3.

    Returns
    -------
    EdcNhsResult
        Integration outcome; see dataclass for field semantics.

    Raises
    ------
    ValueError
        c_cooh_initial <= 0, T non-finite or non-positive, pH outside
        [3.0, 10.0].
    """
    if c_cooh_initial <= 0:
        raise ValueError(
            "EDC/NHS requires surface COOH: c_cooh_initial must be > 0. "
            "Did you mean to set MaterialProperties.surface_cooh_concentration "
            "on the pre-functionalised (e.g. succinylated) matrix?"
        )
    if not math.isfinite(T) or T <= 0:
        raise ValueError(f"T must be finite and positive, got {T}")
    if pH < 3.0 or pH > 10.0:
        raise ValueError(
            f"pH {pH} outside model validity domain [3.0, 10.0]. "
            "Rate constants are not parameterised outside this range."
        )
    if c_nh2_total < 0 or c_edc_initial < 0 or c_nhs_initial < 0:
        raise ValueError("Concentrations must be non-negative")
    if time < 0:
        raise ValueError("time must be non-negative")

    kin = kin or EdcNhsKinetics()

    # Temperature- and pH-corrected rate constants
    k1 = _arrhenius(T, kin.k1_0, kin.E_a_1)
    k2 = _arrhenius(T, kin.k2_0, kin.E_a_2)
    kh1 = _arrhenius(T, kin.kh1_0, kin.E_a_h1)
    kh2 = _arrhenius(T, kin.kh2_0, kin.E_a_h2) * (10.0 ** (kin.kh2_pH_slope * (pH - 7.0)))
    k3 = _arrhenius(T, kin.k3_0, kin.E_a_3)

    # QSS diagnostic (informational)
    qss_ratio = (k2 * c_nhs_initial + kh1) / max(k1 * c_edc_initial, 1e-300)
    logger.debug(
        "EDC/NHS rate constants @ T=%.1fK pH=%.2f: "
        "k1=%.3e k2=%.3e kh1=%.3e kh2=%.3e k3=%.3e (qss_ratio=%.2f)",
        T, pH, k1, k2, kh1, kh2, k3, qss_ratio,
    )

    amine_fraction = available_amine_fraction(pH, kin.pKa_amine)

    # State vector: [C, A, E, P, EDC_free, NHS_free, NH2_total]
    # (Instead of eliminating via conservation, keep all seven states
    # explicit so the RHS is symmetric and readable — the stiffness is
    # already handled by BDF and the cost is negligible.)
    y0 = np.array([
        c_cooh_initial,     # C
        0.0,                # A
        0.0,                # E
        0.0,                # P
        c_edc_initial,      # EDC_free
        c_nhs_initial,      # NHS_free
        c_nh2_total,        # NH2_total
    ], dtype=float)

    def rhs(_t, y):
        C, A, E, P, edc, nhs, nh2_tot = y
        nh2_eff = nh2_tot * amine_fraction

        # Activation step 1
        r_act = k1 * edc * C
        # Stabilisation step 2
        r_stab = k2 * nhs * A
        # Hydrolyses
        r_h1 = kh1 * A
        r_h2 = kh2 * E
        # Aminolysis step 4
        r_am = k3 * nh2_eff * E

        dC = -r_act + r_h1 + r_h2
        dA = +r_act - r_stab - r_h1
        dE = +r_stab - r_h2 - r_am
        dP = +r_am
        dEDC = -r_act
        dNHS = -r_stab + r_am  # NHS consumed in stabilisation, released in aminolysis
        dNH2 = -r_am

        return np.array([dC, dA, dE, dP, dEDC, dNHS, dNH2])

    # Solve. BDF handles the stiffness introduced by fast A <-> E
    # equilibration relative to slow P accumulation.
    sol = solve_ivp(
        rhs, (0.0, float(time)), y0,
        method="BDF", rtol=rtol, atol=atol,
        # Dense output would be nice but we only need t_end;
        # let the integrator pick its own step count.
    )

    # Extract final state
    if sol.y.shape[1] == 0:
        # solver never stepped; shouldn't happen but guard anyway
        C_f = c_cooh_initial
        A_f = E_f = P_f = 0.0
    else:
        C_f, A_f, E_f, P_f = float(sol.y[0, -1]), float(sol.y[1, -1]), float(sol.y[2, -1]), float(sol.y[3, -1])

    # Mass-balance check
    mb = (C_f + A_f + E_f + P_f) - c_cooh_initial
    mb_err = abs(mb) / max(c_cooh_initial, 1e-300)
    if mb_err > 1e-2 and sol.success:
        logger.warning(
            "EDC/NHS mass conservation slip: C+A+E+P=%.6e, C0=%.6e, rel err %.2e",
            C_f + A_f + E_f + P_f, c_cooh_initial, mb_err,
        )

    # Derived fractions
    p_final = max(P_f / c_cooh_initial, 0.0)
    # Hydrolysed = initial minus (still-intact + committed product)
    # = C0 - (C_f + A_f + E_f + P_f reached-product) — but the C/A/E
    # pool regenerates, so "hydrolysed" here is the total moles lost
    # across all hydrolysis events minus regenerated moles. The
    # invariant C+A+E+P=C0 holds; "hydrolysed pool" at end-time is
    # equivalent to C0 - (A_f + E_f + P_f) - C_f_never_activated, which
    # we approximate as the fraction that ended as C (regenerated) plus
    p_residual_nhs_ester = max(E_f / c_cooh_initial, 0.0)

    # p_hydrolysed: fraction ending as C (regenerated by hydrolysis) or
    # A (about-to-hydrolyse). Captures "not committed to product".
    p_hydrolysed = max((C_f + A_f) / c_cooh_initial, 0.0)

    # Time-to-half: simple linear interp on the P trajectory.
    time_to_half = float("nan")
    if p_final > 1e-6 and sol.y.shape[1] >= 2:
        P_traj = sol.y[3, :]
        t_traj = sol.t
        target = 0.5 * P_f
        idx = np.searchsorted(P_traj, target)
        if 0 < idx < len(P_traj):
            # linear interp between idx-1 and idx
            P0, P1 = P_traj[idx - 1], P_traj[idx]
            t0, t1 = t_traj[idx - 1], t_traj[idx]
            if P1 > P0:
                time_to_half = float(t0 + (target - P0) / (P1 - P0) * (t1 - t0))

    return EdcNhsResult(
        p_final=p_final,
        p_hydrolysed=p_hydrolysed,
        p_residual_nhs_ester=p_residual_nhs_ester,
        time_to_half=time_to_half,
        mass_balance_error=mb_err,
        solver_success=bool(sol.success),
        solver_message=str(sol.message),
        n_steps=int(sol.t.size),
    )
