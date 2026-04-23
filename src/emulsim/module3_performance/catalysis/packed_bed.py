"""Packed-bed catalytic reactor with axial dispersion.

Solves the transient plug-flow reactor (PFR) model for an immobilised-enzyme
packed bed with:
- Michaelis-Menten kinetics
- Internal effectiveness factor (generalised Thiele modulus, spherical geometry)
- First-order enzyme deactivation
- Axial dispersion (correlation: D_ax = u*d_p / (2*eps_b))

Governing equations
-------------------
eps_b * dS/dt = -u * dS/dz + D_ax * d^2S/dz^2
                - (1 - eps_b) * eta * V_max_eff(t) * S / (K_m + S)

dP/dt = (1 - eps_b) * eta * V_max_eff(t) * S / (K_m + S)  (1:1 stoichiometry)

where V_max_eff(t) = V_max * exp(-k_d * t)

Spatial discretisation: first-order upwind for advection, central difference
for dispersion, on a uniform grid of n_z finite-volume cells.
Time integration: scipy.integrate.solve_ivp with BDF method.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np
from scipy.integrate import solve_ivp

from ...datatypes import ModelEvidenceTier, ModelManifest
from .kinetics import (
    effectiveness_factor,
    generalized_thiele_modulus,
    michaelis_menten_rate,
)
from .deactivation import first_order_deactivation


# ---------------------------------------------------------------------------
# Result container
# ---------------------------------------------------------------------------

@dataclass
class CatalyticResult:
    """Results from a packed-bed enzyme reactor simulation.

    Attributes
    ----------
    time : np.ndarray, shape (n_t,)
        Time points [s].
    z : np.ndarray, shape (n_z,)
        Axial positions (cell centres) [m].
    S_outlet : np.ndarray, shape (n_t,)
        Substrate concentration at the outlet vs time [mol/m^3].
    P_outlet : np.ndarray, shape (n_t,)
        Product concentration at the outlet vs time [mol/m^3].
    conversion : float
        Steady-state (or final-time) conversion [-].
    effectiveness_factor : float
        Effectiveness factor at feed conditions [-].
    thiele_modulus : float
        Generalised Thiele modulus at feed conditions [-].
    productivity : float
        Volumetric productivity [mol/(m^3_bed * s)] at the final time.
    activity_history : np.ndarray, shape (n_t,)
        Remaining enzyme activity vs time [-].
    mass_balance_error : float
        Relative mass-balance error [-].
    """

    time: np.ndarray
    z: np.ndarray
    S_outlet: np.ndarray
    P_outlet: np.ndarray
    conversion: float
    effectiveness_factor: float
    thiele_modulus: float
    productivity: float
    activity_history: np.ndarray
    mass_balance_error: float
    # v6.1 (Node 5): evidence provenance — populated by solve_packed_bed.
    # Catalytic kinetics rely on user-supplied V_max/K_m and a Thiele-modulus
    # diagnostic; tier defaults to SEMI_QUANTITATIVE and is gated by mass
    # balance and Thiele regime checks.
    model_manifest: Optional[ModelManifest] = None


# ---------------------------------------------------------------------------
# Solver
# ---------------------------------------------------------------------------

def solve_packed_bed(
    bed_length: float,
    bed_diameter: float,
    particle_diameter: float,
    bed_porosity: float,
    particle_porosity: float,
    V_max: float,
    K_m: float,
    S_feed: float,
    flow_rate: float,
    D_eff: float = 1e-10,
    k_deact: float = 0.0,
    total_time: float = 3600.0,
    n_z: int = 50,
    D_molecular: float = 1e-9,
) -> CatalyticResult:
    """Solve the transient packed-bed enzyme reactor model.

    Parameters
    ----------
    bed_length : float
        Length of the packed bed [m].
    bed_diameter : float
        Inner diameter of the column [m].
    particle_diameter : float
        Diameter of the catalyst particles [m].
    bed_porosity : float
        Inter-particle (bed) void fraction [-], typically 0.3-0.5.
    particle_porosity : float
        Intra-particle porosity [-], typically 0.3-0.7.
    V_max : float
        Maximum reaction rate per unit bead volume [mol/(m^3*s)].
    K_m : float
        Michaelis constant [mol/m^3].
    S_feed : float
        Substrate feed concentration [mol/m^3].
    flow_rate : float
        Volumetric flow rate [m^3/s].
    D_eff : float, optional
        Effective intra-particle diffusivity [m^2/s].  Default 1e-10.
    k_deact : float, optional
        First-order deactivation rate constant [1/s].  Default 0 (no deactivation).
    total_time : float, optional
        Total simulation time [s].  Default 3600.
    n_z : int, optional
        Number of spatial grid cells.  Default 50.
    D_molecular : float, optional
        Molecular diffusivity in the mobile phase [m^2/s].  Default 1e-9.

    Returns
    -------
    CatalyticResult
        Simulation results.
    """
    # ---- Geometry ----------------------------------------------------------
    A_cross = np.pi / 4.0 * bed_diameter**2          # [m^2]
    u = flow_rate / (A_cross * bed_porosity)          # interstitial velocity [m/s]
    R_particle = particle_diameter / 2.0              # [m]
    bed_volume = A_cross * bed_length                 # [m^3]

    # ---- Spatial grid (uniform, cell-centred) ------------------------------
    dz = bed_length / n_z
    z = np.linspace(dz / 2.0, bed_length - dz / 2.0, n_z)

    # ---- Axial dispersion --------------------------------------------------
    # Correlation: D_ax = u * d_p / (2 * eps_b)  (low Pe regime approximation)
    D_ax = u * particle_diameter / (2.0 * bed_porosity)
    D_ax = max(D_ax, D_molecular * 0.1)  # floor to avoid zero dispersion

    # ---- Thiele modulus and effectiveness factor at feed conditions ---------
    phi_gen = generalized_thiele_modulus(R_particle, V_max, K_m, D_eff, S_feed)
    eta_feed = effectiveness_factor(phi_gen)

    # ---- Time integration setup --------------------------------------------
    # State vector: y = [S_1 .. S_Nz,  P_1 .. P_Nz]
    # Product tracked as accumulated product concentration in the mobile phase.
    y0 = np.zeros(2 * n_z)
    # Initial condition: bed is empty (S=0, P=0 everywhere)

    # Peclet-based time step guidance for output
    tau_residence = bed_length / u  # residence time [s]
    n_out = max(100, int(total_time / (tau_residence / 10.0)))
    n_out = min(n_out, 2000)
    t_eval = np.linspace(0.0, total_time, n_out)

    # ---- RHS function ------------------------------------------------------
    def rhs(t: float, y: np.ndarray) -> np.ndarray:
        S = y[:n_z]
        P = y[n_z:]
        dSdt = np.zeros(n_z)
        dPdt = np.zeros(n_z)

        # Enzyme activity at current time — scalar t, so a_t is a scalar.
        a_t = float(first_order_deactivation(t, k_deact))
        V_max_eff = V_max * a_t

        # Compute local effectiveness factor based on local S (vectorised)
        S_safe = np.maximum(S, 1e-30)
        F_local = K_m * S_safe - K_m**2 * np.log(1.0 + S_safe / K_m)
        F_local = np.maximum(F_local, 1e-30)
        phi_local = R_particle * np.sqrt(V_max_eff / (2.0 * D_eff * F_local))
        eta_local = effectiveness_factor(phi_local)

        # Reaction rate per unit bed volume
        # rate_bead = eta * V_max_eff * S / (K_m + S)  [per bead volume]
        # rate_bed  = (1 - eps_b) * rate_bead           [per bed volume]
        rate = (1.0 - bed_porosity) * eta_local * michaelis_menten_rate(
            S_safe, V_max_eff, K_m
        )

        # --- Advection (first-order upwind) ---------------------------------
        # Inlet boundary: S = S_feed (Dirichlet)
        S_face = np.empty(n_z + 1)
        S_face[0] = S_feed          # inlet ghost
        S_face[1:] = S              # upwind: face value = upstream cell value

        adv = -u * (S_face[1:] - S_face[:-1]) / dz  # upwind: (S_i - S_{i-1})/dz

        # Advection for product: P_feed = 0
        P_face = np.empty(n_z + 1)
        P_face[0] = 0.0
        P_face[1:] = P

        adv_P = -u * (P_face[1:] - P_face[:-1]) / dz

        # --- Dispersion (central difference) --------------------------------
        # Ghost cells: inlet S = S_feed, outlet dS/dz = 0 (Neumann)
        S_ext = np.empty(n_z + 2)
        S_ext[0] = S_feed                    # inlet Dirichlet (ghost)
        S_ext[1:-1] = S
        S_ext[-1] = S[-1]                    # outlet Neumann (zero gradient)

        disp = D_ax * (S_ext[2:] - 2.0 * S_ext[1:-1] + S_ext[:-2]) / dz**2

        # Dispersion for product (zero-gradient at both boundaries)
        P_ext = np.empty(n_z + 2)
        P_ext[0] = 0.0
        P_ext[1:-1] = P
        P_ext[-1] = P[-1]

        disp_P = D_ax * (P_ext[2:] - 2.0 * P_ext[1:-1] + P_ext[:-2]) / dz**2

        # --- Assemble -------------------------------------------------------
        dSdt = (adv + disp - rate) / bed_porosity
        dPdt = (adv_P + disp_P + rate) / bed_porosity

        return np.concatenate([dSdt, dPdt])

    # ---- Solve -------------------------------------------------------------
    # Method choice: LSODA auto-switches between stiff (BDF) and non-stiff
    # (Adams) and runs ~700× faster than fixed BDF on this advection-
    # dispersion-reaction problem (v9.1.1 issue #2 investigation, Apr 2026).
    # The PFR + Michaelis-Menten + dispersion system is mostly non-stiff at
    # the resolution used here; BDF was wasting work computing finite-
    # difference Jacobians that LSODA only does on demand.
    sol = solve_ivp(
        rhs,
        (0.0, total_time),
        y0,
        method="LSODA",
        t_eval=t_eval,
        rtol=1e-6,
        atol=1e-9,
        max_step=tau_residence / 2.0,
    )

    if not sol.success:
        raise RuntimeError(f"ODE solver failed: {sol.message}")

    # ---- Extract results ---------------------------------------------------
    S_all = sol.y[:n_z, :]   # (n_z, n_t)
    P_all = sol.y[n_z:, :]   # (n_z, n_t)

    S_outlet = S_all[-1, :]
    P_outlet = P_all[-1, :]

    # Steady-state conversion (use last time point)
    conversion = float(1.0 - S_outlet[-1] / S_feed) if S_feed > 0 else 0.0
    conversion = max(0.0, min(1.0, conversion))

    # Activity history
    activity = first_order_deactivation(sol.t, k_deact)

    # Productivity: moles of product per bed volume per second at final time
    productivity = float(P_outlet[-1] * flow_rate / bed_volume) if bed_volume > 0 else 0.0

    # Mass balance error
    # At steady state: S_in_rate = S_out_rate + P_out_rate
    S_in_rate = S_feed * flow_rate
    S_out_rate = S_outlet[-1] * flow_rate
    P_out_rate = P_outlet[-1] * flow_rate
    total_out = S_out_rate + P_out_rate
    if S_in_rate > 0:
        mass_balance_error = abs(S_in_rate - total_out) / S_in_rate
    else:
        mass_balance_error = 0.0

    manifest = _build_catalytic_manifest(
        thiele=float(phi_gen),
        eta=float(eta_feed),
        conversion=conversion,
        mass_balance_error=mass_balance_error,
        k_deact=k_deact,
    )

    return CatalyticResult(
        time=sol.t,
        z=z,
        S_outlet=S_outlet,
        P_outlet=P_outlet,
        conversion=conversion,
        effectiveness_factor=float(eta_feed),
        thiele_modulus=float(phi_gen),
        productivity=productivity,
        activity_history=np.asarray(activity),
        mass_balance_error=mass_balance_error,
        model_manifest=manifest,
    )


# ─── Node 5: catalytic-result manifest ────────────────────────────────────


# Catalytic mass-balance gates mirror the chromatography gates in the M3
# orchestrator: <=2% ok, 2-5% caution, >5% capped at QUALITATIVE_TREND.
_CAT_MB_CAUTION = 0.02
_CAT_MB_BLOCKER = 0.05

# Thiele-modulus regime cuts: phi <= 0.3 reaction-limited, > 1.0 strongly
# diffusion-limited (effectiveness factor << 1, predictions become sensitive
# to D_eff which the user usually cannot calibrate locally).
_THIELE_REACTION_LIMIT = 0.3
_THIELE_DIFFUSION_LIMIT = 1.0


def _build_catalytic_manifest(
    thiele: float,
    eta: float,
    conversion: float,
    mass_balance_error: float,
    k_deact: float,
) -> ModelManifest:
    """Build the CatalyticResult ModelManifest.

    Tier policy:
      * Default SEMI_QUANTITATIVE — V_max and K_m are user-supplied kinetics
        with no global calibration framework wired in v6.1.
      * Cap at QUALITATIVE_TREND when:
          - mass_balance_error > 5% (solver lost mass, rates not defendable);
          - phi > 1.0 AND eta < 0.5 (strongly diffusion-limited;
            effectiveness factor sensitive to D_eff which is rarely measured).
      * If conversion is essentially zero AND no deactivation was modelled,
        the run produced no chemistry — flag as UNSUPPORTED for downstream.
    """
    tier = ModelEvidenceTier.SEMI_QUANTITATIVE
    diagnostics: dict = {
        "thiele_modulus": float(thiele),
        "effectiveness_factor": float(eta),
        "conversion": float(conversion),
        "mass_balance_error": float(mass_balance_error),
        "deactivation_active": float(k_deact) > 0.0,
    }
    assumptions: list[str] = [
        "Michaelis-Menten kinetics with user-supplied V_max and K_m "
        "(no calibration framework wired in v6.1).",
    ]

    # Regime classification.
    if thiele <= _THIELE_REACTION_LIMIT:
        diagnostics["regime"] = "reaction_limited"
    elif thiele <= _THIELE_DIFFUSION_LIMIT:
        diagnostics["regime"] = "borderline"
    else:
        diagnostics["regime"] = "diffusion_limited"

    # Mass-balance gate.
    if mass_balance_error > _CAT_MB_BLOCKER:
        diagnostics["mass_balance_status"] = "blocker"
        tier = ModelEvidenceTier.QUALITATIVE_TREND
        assumptions.append(
            f"Mass balance error {mass_balance_error:.1%} exceeds {_CAT_MB_BLOCKER:.0%}; "
            "conversion and productivity are not decision-grade."
        )
    elif mass_balance_error > _CAT_MB_CAUTION:
        diagnostics["mass_balance_status"] = "caution"
    else:
        diagnostics["mass_balance_status"] = "ok"

    # Diffusion-limited regime gate.
    if thiele > _THIELE_DIFFUSION_LIMIT and eta < 0.5:
        # Tier is already SEMI; only weaken if not already capped.
        _ORDER = list(ModelEvidenceTier)
        capped = max(
            _ORDER.index(tier),
            _ORDER.index(ModelEvidenceTier.QUALITATIVE_TREND),
        )
        tier = _ORDER[capped]
        assumptions.append(
            f"Strongly diffusion-limited (Thiele={thiele:.2f}, eta={eta:.2f}); "
            "predictions are sensitive to D_eff which is typically uncalibrated."
        )

    # No-chemistry sentinel: if the solver produced essentially zero conversion
    # and no kinetic activity was simulated (no deactivation either), the run
    # carries no useful chemistry signal.
    if conversion <= 1e-6 and k_deact <= 0.0:
        tier = ModelEvidenceTier.UNSUPPORTED
        assumptions.append("No measurable conversion; treat as UNSUPPORTED.")

    return ModelManifest(
        model_name="M3.catalysis.packed_bed.MM",
        evidence_tier=tier,
        valid_domain={},
        calibration_ref="",
        assumptions=assumptions,
        diagnostics=diagnostics,
    )
