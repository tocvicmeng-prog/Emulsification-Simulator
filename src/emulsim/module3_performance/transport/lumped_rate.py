"""Lumped Rate Model (LRM) for chromatographic column transport.

Architecture: docs/13_module2_module3_final_implementation_plan.md, Phase C.

Governing equations (single component):

    eps_b * dC/dt = -u * dC/dz + D_ax * d2C/dz2
                    - (1 - eps_b) * (3/R_p) * k_f * (C - C_p)

    eps_p * dC_p/dt = (3/R_p) * k_f * (C - C_p)
                      - (1 - eps_p) * dq/dt

    dq/dt = k_ads * (q_eq(C_p) - q)

Discretization:
    - Finite volume in z: N_z cells (default 50)
    - Advection: first-order upwind (u * dC/dz)
    - Dispersion: central differences (D_ax * d2C/dz2)
    - BCs: Danckwerts inlet, zero-gradient outlet
    - Time: scipy solve_ivp with BDF (stiff system)
"""

from __future__ import annotations

import logging
import warnings
from dataclasses import dataclass, field

import numpy as np
from scipy.integrate import solve_ivp

from ..hydrodynamics import ColumnGeometry
from ..isotherms.langmuir import LangmuirIsotherm

logger = logging.getLogger(__name__)


# ─── Result Dataclass ─────────────────────────────────────────────────

@dataclass
class LRMResult:
    """Result container for a Lumped Rate Model simulation.

    Attributes:
        time: Time points [s], shape (N_t,).
        z: Cell-centre positions [m], shape (N_z,).
        C_outlet: Outlet concentration vs time [mol/m^3], shape (N_t,).
        q_average: Volume-averaged bound concentration [mol/m^3], shape (N_t,).
        mass_injected: Total mass injected [mol].
        mass_eluted: Total mass eluted [mol].
        mass_bound: Total mass currently bound in the column [mol].
        mass_balance_error: Relative mass balance error [-].
    """

    time: np.ndarray
    z: np.ndarray
    C_outlet: np.ndarray
    q_average: np.ndarray
    mass_injected: float
    mass_eluted: float
    mass_bound: float
    mass_balance_error: float


# ─── Correlation Functions ────────────────────────────────────────────

def _wilson_geankoplis_kf(
    u: float,
    dp: float,
    eps_b: float,
    D_m: float,
    mu: float = 1e-3,
    rho: float = 1000.0,
) -> float:
    """Film mass transfer coefficient from Wilson-Geankoplis correlation.

    Sh = (1.09 / eps_b) * (Re * Sc)^(1/3)
    k_f = Sh * D_m / dp

    Valid for Re_p < 55 and 0.0016 < Re*Sc < 55000.

    Args:
        u: Superficial velocity [m/s].
        dp: Particle diameter [m].
        eps_b: Bed porosity [-].
        D_m: Molecular diffusivity [m^2/s].
        mu: Dynamic viscosity [Pa.s].
        rho: Fluid density [kg/m^3].

    Returns:
        Film mass transfer coefficient k_f [m/s].
    """
    Re = rho * u * dp / mu
    Sc = mu / (rho * D_m)
    ReSc = Re * Sc

    # Guard against zero flow
    if ReSc < 1e-20:
        return D_m / dp  # diffusion limit

    Sh = (1.09 / eps_b) * ReSc ** (1.0 / 3.0)
    # Minimum Sherwood number = 2 (sphere in stagnant fluid)
    Sh = max(Sh, 2.0)
    return Sh * D_m / dp


def _axial_dispersion(u: float, dp: float, eps_b: float) -> float:
    """Axial dispersion coefficient from Peclet ~ 2 approximation.

    D_ax = u * dp / (2 * eps_b)

    Args:
        u: Superficial velocity [m/s].
        dp: Particle diameter [m].
        eps_b: Bed porosity [-].

    Returns:
        Axial dispersion coefficient D_ax [m^2/s].
    """
    if u < 1e-20:
        return 1e-10  # minimal dispersion for zero flow
    return u * dp / (2.0 * eps_b)


# ─── ODE RHS ─────────────────────────────────────────────────────────

def _build_rhs(
    n_z: int,
    dz: float,
    u: float,
    D_ax: float,
    eps_b: float,
    eps_p: float,
    R_p: float,
    k_f: float,
    k_ads: float,
    isotherm: LangmuirIsotherm,
    C_feed: float,
    feed_duration: float,
):
    """Build the RHS function for solve_ivp.

    State vector layout: y = [C_0..C_{Nz-1}, Cp_0..Cp_{Nz-1}, q_0..q_{Nz-1}]
    Total size: 3 * N_z.
    """
    mass_transfer_coeff = (3.0 / R_p) * k_f  # [1/s] * [1/m] -> combined

    def rhs(t: float, y: np.ndarray) -> np.ndarray:
        # Unpack state
        C = y[:n_z].copy()
        Cp = y[n_z:2 * n_z].copy()
        q = y[2 * n_z:3 * n_z].copy()

        # Positivity enforcement
        np.maximum(C, 0.0, out=C)
        np.maximum(Cp, 0.0, out=Cp)
        np.maximum(q, 0.0, out=q)

        # Inlet concentration (step function)
        C_in = C_feed if t <= feed_duration else 0.0

        # ── Advection: first-order upwind (u > 0) ──
        # dC/dz ~ (C[i] - C[i-1]) / dz
        dCdz = np.empty(n_z)
        # Danckwerts inlet: flux balance at z=0
        # u*C_in = u*C[0] - D_ax*(C[1]-C[0])/dz  (approx)
        # => ghost cell C_ghost = C_in (after rearranging with D_ax term in dispersion)
        dCdz[0] = (C[0] - C_in) / dz
        dCdz[1:] = (C[1:] - C[:-1]) / dz

        # ── Dispersion: central differences ──
        # d2C/dz2 ~ (C[i+1] - 2*C[i] + C[i-1]) / dz^2
        d2Cdz2 = np.empty(n_z)
        # Inlet ghost: C_ghost_in = C_in (Danckwerts)
        d2Cdz2[0] = (C[1] - 2.0 * C[0] + C_in) / dz ** 2
        # Interior
        d2Cdz2[1:-1] = (C[2:] - 2.0 * C[1:-1] + C[:-2]) / dz ** 2
        # Outlet: zero-gradient => C[Nz] = C[Nz-1]
        d2Cdz2[-1] = (C[-1] - 2.0 * C[-1] + C[-2]) / dz ** 2  # = (C[-2] - C[-1]) / dz^2

        # ── Isotherm ──
        q_eq = isotherm.equilibrium_loading(Cp)

        # ── Kinetics ──
        dqdt = k_ads * (q_eq - q)

        # ── Film mass transfer ──
        film_flux = mass_transfer_coeff * (C - Cp)  # [mol/m^3/s]

        # ── PDEs ──
        dCdt = (-u * dCdz + D_ax * d2Cdz2 - (1.0 - eps_b) * film_flux) / eps_b
        dCpdt = (film_flux - (1.0 - eps_p) * dqdt) / eps_p
        dqdt_out = dqdt

        return np.concatenate([dCdt, dCpdt, dqdt_out])

    return rhs


# ─── Solver ───────────────────────────────────────────────────────────

def solve_lrm(
    column: ColumnGeometry,
    isotherm: LangmuirIsotherm,
    C_feed: float,
    feed_duration: float,
    flow_rate: float,
    total_time: float,
    n_z: int = 50,
    D_molecular: float = 7e-11,
    k_ads: float = 100.0,
    mu: float = 1e-3,
    rho: float = 1000.0,
    rtol: float = 1e-6,
    atol: float = 1e-9,
) -> LRMResult:
    """Solve the Lumped Rate Model for a single-component breakthrough.

    Args:
        column: Packed column geometry.
        isotherm: Langmuir isotherm parameters.
        C_feed: Feed concentration [mol/m^3].
        feed_duration: Duration of feed step [s].
        flow_rate: Volumetric flow rate [m^3/s].
        total_time: Total simulation time [s].
        n_z: Number of finite-volume cells in z.
        D_molecular: Molecular diffusivity [m^2/s].
        k_ads: Adsorption rate constant [1/s].
        mu: Dynamic viscosity [Pa.s].
        rho: Fluid density [kg/m^3].
        rtol: Relative tolerance for ODE solver.
        atol: Absolute tolerance for ODE solver.

    Returns:
        LRMResult with breakthrough curves and mass balance.
    """
    # ── Derived quantities ──
    u = column.superficial_velocity(flow_rate)
    eps_b = column.bed_porosity
    eps_p = column.particle_porosity
    R_p = column.particle_radius
    dp = column.particle_diameter
    L = column.bed_height

    # Grid
    dz = L / n_z
    z = np.linspace(dz / 2.0, L - dz / 2.0, n_z)  # cell centres

    # CFL check
    dt_cfl = dz / u if u > 0 else 1e10
    if dt_cfl < 0.01:
        logger.warning(
            "CFL condition tight: dz/u = %.3e s. "
            "BDF solver handles this, but consider more cells.", dt_cfl
        )

    # Correlations
    k_f = _wilson_geankoplis_kf(u, dp, eps_b, D_molecular, mu, rho)
    D_ax = _axial_dispersion(u, dp, eps_b)

    logger.info(
        "LRM setup: u=%.3e m/s, k_f=%.3e m/s, D_ax=%.3e m^2/s, "
        "N_z=%d, dz=%.3e m",
        u, k_f, D_ax, n_z, dz,
    )

    # ── Build RHS ──
    rhs = _build_rhs(
        n_z=n_z,
        dz=dz,
        u=u,
        D_ax=D_ax,
        eps_b=eps_b,
        eps_p=eps_p,
        R_p=R_p,
        k_f=k_f,
        k_ads=k_ads,
        isotherm=isotherm,
        C_feed=C_feed,
        feed_duration=feed_duration,
    )

    # ── Initial conditions: column empty ──
    y0 = np.zeros(3 * n_z)

    # ── Time span ──
    # Dense output at ~200 points minimum
    n_eval = max(200, int(total_time / 0.5))
    t_eval = np.linspace(0.0, total_time, n_eval)

    # ── Solve ──
    sol = solve_ivp(
        rhs,
        t_span=(0.0, total_time),
        y0=y0,
        method="BDF",
        t_eval=t_eval,
        rtol=rtol,
        atol=atol,
        max_step=total_time / 20.0,
    )

    if not sol.success:
        raise RuntimeError(f"LRM solver failed: {sol.message}")

    # ── Extract results ──
    time = sol.t
    C_all = sol.y[:n_z, :]        # (N_z, N_t)
    Cp_all = sol.y[n_z:2*n_z, :]  # (N_z, N_t)
    q_all = sol.y[2*n_z:, :]      # (N_z, N_t)

    # Positivity clip on outputs
    C_all = np.maximum(C_all, 0.0)
    Cp_all = np.maximum(Cp_all, 0.0)
    q_all = np.maximum(q_all, 0.0)

    # Outlet = last cell
    C_outlet = C_all[-1, :]

    # Volume-averaged bound concentration
    q_average = np.mean(q_all, axis=0)

    # ── Mass balance ──
    A_cross = column.cross_section_area
    V_bed = column.bed_volume

    # Mass injected: integral of Q * C_feed over feed duration
    mass_injected = flow_rate * C_feed * feed_duration  # [mol]

    # Mass eluted: integral of Q * C_outlet(t) dt (trapezoidal)
    mass_eluted = float(np.trapezoid(C_outlet * flow_rate, time))  # [mol]

    # Mass currently in column:
    # Mobile phase in bed voids
    mass_mobile = float(np.mean(C_all[:, -1])) * eps_b * V_bed
    # Pore liquid
    mass_pore = float(np.mean(Cp_all[:, -1])) * eps_p * (1.0 - eps_b) * V_bed
    # Bound phase
    mass_bound = float(np.mean(q_all[:, -1])) * (1.0 - eps_p) * (1.0 - eps_b) * V_bed

    mass_in_column = mass_mobile + mass_pore + mass_bound

    # Relative error
    if mass_injected > 0:
        mass_balance_error = abs(mass_injected - mass_eluted - mass_in_column) / mass_injected
    else:
        mass_balance_error = 0.0

    if mass_balance_error > 0.02:
        warnings.warn(
            f"Mass balance error {mass_balance_error:.2%} exceeds 2% threshold. "
            f"Injected={mass_injected:.4e}, Eluted={mass_eluted:.4e}, "
            f"In column={mass_in_column:.4e}",
            stacklevel=2,
        )

    logger.info(
        "LRM solved: %d time steps, mass balance error = %.4f%%",
        len(time), mass_balance_error * 100,
    )

    return LRMResult(
        time=time,
        z=z,
        C_outlet=C_outlet,
        q_average=q_average,
        mass_injected=mass_injected,
        mass_eluted=mass_eluted,
        mass_bound=mass_bound,
        mass_balance_error=mass_balance_error,
    )
