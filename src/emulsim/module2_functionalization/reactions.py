"""Chemistry-agnostic sequential reaction engine for Module 2 functionalization.

Phase B: Minimal Module 2 — 2 Workflows.
Architecture: module2_module3_final_implementation_plan.md, Phase B.

Provides three ODE templates for surface modification chemistry:

Template 1 — Simple irreversible second-order (genipin, DVS, ECH crosslinking):
    d[ACS]/dt = -k * [reagent] * [ACS]
    d[reagent]/dt = -stoich * k * [reagent] * [ACS]

Template 2 — Competitive hydrolysis (NHS, CNBr — Phase E, template defined now):
    d[active]/dt = -k_couple * [active] * [ligand] - k_hydrol * [active]
    d[coupled]/dt = +k_couple * [active] * [ligand]
    d[hydrolyzed]/dt = +k_hydrol * [active]

Template 3 — Equilibrium binding with steric blocking (protein coupling):
    d[coupled]/dt = k_couple * [active_sites] * [ligand_free] * f_steric
    f_steric = max(1 - coupled/max_sites_per_area, 0)

Uses scipy.integrate.solve_ivp with Radau method (stiff-safe, same as L3).
"""

from __future__ import annotations

import logging
import math

import numpy as np
from scipy.integrate import solve_ivp

logger = logging.getLogger(__name__)

# ─── Constants ──────────────────────────────────────────────────────────

R_GAS = 8.314  # [J/(mol*K)]


# ─── Arrhenius helper ───────────────────────────────────────────────────

def arrhenius_rate_constant(T: float, k0: float, E_a: float) -> float:
    """Arrhenius rate constant.

    k(T) = k0 * exp(-E_a / (R*T))

    Args:
        T: Temperature [K].
        k0: Arrhenius prefactor [m^3/(mol*s)].
        E_a: Activation energy [J/mol].

    Returns:
        Rate constant at temperature T [m^3/(mol*s)].
    """
    if k0 <= 0 or T <= 0:
        return 0.0
    return k0 * math.exp(-E_a / (R_GAS * T))


# ─── Template 1: Simple irreversible second-order ──────────────────────

def _second_order_rhs(
    t: float,
    y: np.ndarray,
    k: float,
    stoichiometry: float,
    hydrolysis_rate: float,
) -> np.ndarray:
    """RHS of the second-order consumption ODE with optional hydrolysis.

    State variables y = [ACS, reagent]:
        ACS: available coupling site concentration [mol/m^3]
        reagent: reagent concentration [mol/m^3]

    d[ACS]/dt = -k * [reagent] * [ACS]
    d[reagent]/dt = -stoich * k * [reagent] * [ACS] - k_hydrol * [reagent]

    The hydrolysis term (-k_hydrol * [reagent]) models first-order
    degradation of the reagent competing with the coupling reaction.

    Args:
        t: Time [s] (unused, autonomous system).
        y: State vector [ACS, reagent].
        k: Forward rate constant [m^3/(mol*s)].
        stoichiometry: Moles reagent consumed per mole ACS consumed [-].
        hydrolysis_rate: First-order hydrolysis rate constant [1/s].

    Returns:
        dy/dt array of shape (2,).
    """
    acs = max(y[0], 0.0)
    reagent = max(y[1], 0.0)

    coupling_rate = k * reagent * acs
    hydrol_rate = hydrolysis_rate * reagent

    d_acs = -coupling_rate
    d_reagent = -stoichiometry * coupling_rate - hydrol_rate

    return np.array([d_acs, d_reagent])


def solve_second_order_consumption(
    acs_concentration: float,
    reagent_concentration: float,
    k_forward: float,
    stoichiometry: float,
    time: float,
    temperature: float,
    E_a: float = 0.0,
    k0: float = 0.0,
    hydrolysis_rate: float = 0.0,
) -> tuple[float, float]:
    """Solve second-order ACS consumption ODE.

    Uses scipy.integrate.solve_ivp with Radau method (implicit,
    A-stable, suitable for stiff kinetics).

    Args:
        acs_concentration: Initial ACS bulk concentration [mol/m^3].
        reagent_concentration: Initial reagent concentration [mol/m^3].
        k_forward: Forward rate constant [m^3/(mol*s)] at reference T.
            Used directly when E_a == 0 and k0 == 0.
        stoichiometry: Moles reagent consumed per mole ACS consumed [-].
        time: Reaction time [s].
        temperature: Reaction temperature [K].
        E_a: Activation energy [J/mol].  If > 0 and k0 > 0, Arrhenius
            rate is used instead of k_forward.
        k0: Arrhenius prefactor [m^3/(mol*s)].
        hydrolysis_rate: First-order hydrolysis rate constant [1/s].
            Models competing degradation of the reagent (e.g., ECH
            epoxide ring opening by water at alkaline pH).

    Returns:
        Tuple of (conversion, reagent_remaining_fraction):
            conversion: Fraction of initial ACS consumed, in [0, 1].
            reagent_remaining_fraction: Fraction of initial reagent
                remaining, in [0, 1].

    Raises:
        ValueError: If any input concentration is negative.
    """
    # --- Input validation ---
    if acs_concentration < 0:
        raise ValueError(f"acs_concentration must be >= 0, got {acs_concentration}")
    if reagent_concentration < 0:
        raise ValueError(f"reagent_concentration must be >= 0, got {reagent_concentration}")
    if time < 0:
        raise ValueError(f"time must be >= 0, got {time}")
    if temperature <= 0:
        raise ValueError(f"temperature must be > 0, got {temperature}")

    # Degenerate cases: no reaction possible
    if acs_concentration == 0 or reagent_concentration == 0 or time == 0:
        return 0.0, 1.0

    # --- Rate constant selection ---
    if E_a > 0 and k0 > 0:
        k = arrhenius_rate_constant(temperature, k0, E_a)
    else:
        k = k_forward

    if k <= 0:
        return 0.0, 1.0

    # --- ODE integration ---
    y0 = np.array([acs_concentration, reagent_concentration])

    sol = solve_ivp(
        fun=_second_order_rhs,
        t_span=(0.0, time),
        y0=y0,
        method="Radau",
        args=(k, stoichiometry, hydrolysis_rate),
        rtol=1e-8,
        atol=1e-12,
        max_step=time / 10.0,
    )

    if not sol.success:
        logger.warning(
            "ODE solver did not converge: %s. Returning partial result.",
            sol.message,
        )

    # --- Extract final state ---
    acs_final = max(sol.y[0, -1], 0.0)
    reagent_final = max(sol.y[1, -1], 0.0)

    # Conversion: fraction of ACS consumed
    conversion = 1.0 - acs_final / acs_concentration
    conversion = max(0.0, min(conversion, 1.0))

    # Reagent remaining fraction
    reagent_remaining = reagent_final / reagent_concentration
    reagent_remaining = max(0.0, min(reagent_remaining, 1.0))

    logger.debug(
        "solve_second_order: k=%.3e, conversion=%.4f, reagent_remaining=%.4f",
        k, conversion, reagent_remaining,
    )

    return conversion, reagent_remaining


# ─── Template 2: Competitive hydrolysis (Phase E, template only) ───────

def _competitive_hydrolysis_rhs(
    t: float,
    y: np.ndarray,
    k_couple: float,
    k_hydrol: float,
) -> np.ndarray:
    """RHS for competitive hydrolysis template.

    State variables y = [active, ligand, coupled, hydrolyzed]:
        active: activated site concentration [mol/m^3]
        ligand: free ligand concentration [mol/m^3]
        coupled: coupled product [mol/m^3]
        hydrolyzed: hydrolyzed (wasted) sites [mol/m^3]

    Args:
        t: Time [s].
        y: State vector [active, ligand, coupled, hydrolyzed].
        k_couple: Coupling rate constant [m^3/(mol*s)].
        k_hydrol: Hydrolysis rate constant [1/s].

    Returns:
        dy/dt array of shape (4,).
    """
    active = max(y[0], 0.0)
    ligand = max(y[1], 0.0)

    r_couple = k_couple * active * ligand
    r_hydrol = k_hydrol * active

    return np.array([
        -r_couple - r_hydrol,   # d[active]/dt
        -r_couple,              # d[ligand]/dt
        +r_couple,              # d[coupled]/dt
        +r_hydrol,              # d[hydrolyzed]/dt
    ])


# ─── Template 3: Equilibrium binding with steric blocking ─────────────

def _steric_binding_rhs(
    t: float,
    y: np.ndarray,
    k_couple: float,
    max_sites: float,
) -> np.ndarray:
    """RHS for equilibrium binding with steric blocking template.

    State variables y = [coupled, active, ligand_free]:
        coupled: coupled sites [mol/m^3]
        active: available active sites [mol/m^3]
        ligand_free: free ligand in solution [mol/m^3]

    Steric factor: f_steric = max(1 - coupled/max_sites, 0)

    Args:
        t: Time [s].
        y: State vector [coupled, active, ligand_free].
        k_couple: Coupling rate constant [m^3/(mol*s)].
        max_sites: Maximum coupling density [mol/m^3].

    Returns:
        dy/dt array of shape (3,).
    """
    coupled = max(y[0], 0.0)
    active = max(y[1], 0.0)
    ligand_free = max(y[2], 0.0)

    # Steric blocking factor: approaches 0 as surface saturates
    f_steric = max(1.0 - coupled / max_sites, 0.0) if max_sites > 0 else 0.0

    rate = k_couple * active * ligand_free * f_steric

    return np.array([
        +rate,   # d[coupled]/dt
        -rate,   # d[active]/dt
        -rate,   # d[ligand_free]/dt
    ])
