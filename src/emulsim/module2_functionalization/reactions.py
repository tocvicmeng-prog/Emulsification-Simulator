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
from dataclasses import dataclass, field

import numpy as np
from scipy.integrate import solve_ivp

logger = logging.getLogger(__name__)


# ─── CouplingResult Dataclass ──────────────────────────────────────────

@dataclass
class CouplingResult:
    """Structured output from coupling/quenching ODE solvers."""
    conversion: float = 0.0            # [0,1] fraction of target sites consumed
    sites_consumed: float = 0.0        # [mol/particle] absolute sites consumed
    sites_hydrolyzed: float = 0.0      # [mol/particle] lost to hydrolysis
    sites_coupled: float = 0.0         # [mol/particle] successfully coupled to ligand
    sites_blocked: float = 0.0         # [mol/particle] capped by quenching
    reagent_remaining_fraction: float = 0.0  # [0,1]
    hydrolysis_fraction: float = 0.0   # Fraction of consumed activated sites lost to hydrolysis
    solver_success: bool = True
    solver_message: str = ""
    site_balance_error: float = 0.0    # Residual of conservation check
    warnings: list[str] = field(default_factory=list)


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
    # Guard against overflow in both directions (Codex F2 fix)
    exponent = -E_a / (R_GAS * T)
    if exponent < -700:  # exp(-700) ~ 1e-304, underflow to zero
        return 0.0
    if exponent > 700:   # exp(700) ~ 1e304, overflow guard (negative E_a)
        return 1e300     # cap at large-but-finite rate
    return k0 * math.exp(exponent)


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


# ═════════════════════════════════════════════════════════════════════════
# Public ODE Wrappers (W4 — returns CouplingResult with diagnostics)
# ═════════════════════════════════════════════════════════════════════════


def solve_competitive_coupling(
    activated_sites_mol_per_particle: float,
    ligand_concentration: float,
    k_forward: float,
    E_a: float,
    temperature: float,
    time: float,
    bead_volume: float,
    hydrolysis_rate: float = 0.0,
    stoichiometry: float = 1.0,
    ph: float = 7.0,
    ph_min: float = 0.0,
    ph_max: float = 14.0,
) -> CouplingResult:
    """Solve ligand coupling with competitive hydrolysis (Template 2 wrapper).

    Converts mol/particle <-> mol/m^3 for the ODE, then back.
    Returns structured CouplingResult with hydrolysis tracking.

    Args:
        activated_sites_mol_per_particle: Remaining activated sites [mol/particle].
        ligand_concentration: Free ligand in solution [mol/m^3].
        k_forward: Forward rate constant [m^3/(mol*s)] at reference T.
        E_a: Activation energy [J/mol].
        temperature: Reaction temperature [K].
        time: Reaction time [s].
        bead_volume: Volume of one bead [m^3/particle].
        hydrolysis_rate: First-order hydrolysis rate constant [1/s].
        stoichiometry: Moles ligand per mole site consumed [-].
        ph: Reaction pH.
        ph_min: Minimum valid pH for this chemistry.
        ph_max: Maximum valid pH for this chemistry.

    Returns:
        CouplingResult with conversion, sites_coupled, sites_hydrolyzed, diagnostics.
    """
    warnings: list[str] = []

    # pH validity gate
    if ph < ph_min or ph > ph_max:
        warnings.append(
            f"pH {ph:.1f} outside valid range [{ph_min:.1f}, {ph_max:.1f}]; "
            f"rate constants may not apply."
        )

    # Degenerate cases
    if activated_sites_mol_per_particle <= 0 or ligand_concentration <= 0 or time <= 0:
        return CouplingResult(warnings=warnings)
    if bead_volume <= 0:
        return CouplingResult(solver_success=False, solver_message="bead_volume <= 0",
                              warnings=warnings)

    # Arrhenius rate constant
    # Back-calculate Arrhenius prefactor; clamp to prevent overflow
    _exp_ref = math.exp(max(-E_a / (R_GAS * 298.15), -700))
    k0 = min(k_forward / max(_exp_ref, 1e-30), 1e30) if E_a > 0 else 0.0
    k = arrhenius_rate_constant(temperature, k0, E_a) if E_a > 0 and k0 > 0 else k_forward

    if k <= 0:
        return CouplingResult(warnings=warnings)

    # Convert to bulk concentration [mol/m^3]
    active_conc = activated_sites_mol_per_particle / bead_volume

    # ODE initial state: [active, ligand, coupled, hydrolyzed]
    y0 = np.array([active_conc, ligand_concentration, 0.0, 0.0])

    sol = solve_ivp(
        fun=_competitive_hydrolysis_rhs,
        t_span=(0.0, time),
        y0=y0,
        method="Radau",
        args=(k, hydrolysis_rate),
        rtol=1e-8,
        atol=1e-12,
        max_step=time / 10.0,
    )

    if not sol.success:
        warnings.append(f"ODE solver warning: {sol.message}")

    # Extract final state [mol/m^3]
    active_final = max(sol.y[0, -1], 0.0)
    coupled_final = max(sol.y[2, -1], 0.0)
    hydrolyzed_final = max(sol.y[3, -1], 0.0)
    ligand_final = max(sol.y[1, -1], 0.0)

    # Convert back to mol/particle
    sites_coupled = coupled_final * bead_volume
    sites_hydrolyzed = hydrolyzed_final * bead_volume
    sites_consumed = sites_coupled + sites_hydrolyzed

    # Conversion
    conversion = sites_consumed / activated_sites_mol_per_particle if activated_sites_mol_per_particle > 0 else 0.0
    conversion = max(0.0, min(conversion, 1.0))

    # Hydrolysis fraction
    hydrol_frac = sites_hydrolyzed / sites_consumed if sites_consumed > 0 else 0.0

    # Reagent remaining
    reagent_rem = ligand_final / ligand_concentration if ligand_concentration > 0 else 0.0

    # Balance check: active_consumed should equal coupled + hydrolyzed
    active_consumed_conc = active_conc - active_final
    balance_error = abs(active_consumed_conc - coupled_final - hydrolyzed_final) / max(active_conc, 1e-30)

    return CouplingResult(
        conversion=conversion,
        sites_consumed=sites_consumed,
        sites_coupled=sites_coupled,
        sites_hydrolyzed=sites_hydrolyzed,
        reagent_remaining_fraction=max(0.0, min(reagent_rem, 1.0)),
        hydrolysis_fraction=hydrol_frac,
        solver_success=sol.success,
        solver_message=sol.message if not sol.success else "",
        site_balance_error=balance_error,
        warnings=warnings,
    )


def solve_steric_coupling(
    activated_sites_mol_per_particle: float,
    ligand_concentration: float,
    k_forward: float,
    E_a: float,
    temperature: float,
    time: float,
    bead_volume: float,
    max_coupled_mol_per_particle: float,
    ph: float = 7.0,
    ph_min: float = 0.0,
    ph_max: float = 14.0,
) -> CouplingResult:
    """Solve protein coupling with steric blocking (Template 3 wrapper).

    Uses RSA-like steric factor: f_steric = max(1 - coupled/max_sites, 0).
    Canonical unit = mol/particle; converts to mol/m^3 for ODE.

    Args:
        activated_sites_mol_per_particle: Remaining activated sites [mol/particle].
        ligand_concentration: Free protein in solution [mol/m^3].
        k_forward: Forward rate constant [m^3/(mol*s)] at reference T.
        E_a: Activation energy [J/mol].
        temperature: Reaction temperature [K].
        time: Reaction time [s].
        bead_volume: Volume of one bead [m^3/particle].
        max_coupled_mol_per_particle: Steric jamming limit [mol/particle].
        ph: Reaction pH.
        ph_min: Minimum valid pH.
        ph_max: Maximum valid pH.

    Returns:
        CouplingResult with conversion, sites_coupled, diagnostics.
    """
    warnings: list[str] = []

    if ph < ph_min or ph > ph_max:
        warnings.append(
            f"pH {ph:.1f} outside valid range [{ph_min:.1f}, {ph_max:.1f}]."
        )

    if activated_sites_mol_per_particle <= 0 or ligand_concentration <= 0 or time <= 0:
        return CouplingResult(warnings=warnings)
    if bead_volume <= 0:
        return CouplingResult(solver_success=False, solver_message="bead_volume <= 0",
                              warnings=warnings)

    # Arrhenius
    # Back-calculate Arrhenius prefactor; clamp to prevent overflow
    _exp_ref = math.exp(max(-E_a / (R_GAS * 298.15), -700))
    k0 = min(k_forward / max(_exp_ref, 1e-30), 1e30) if E_a > 0 else 0.0
    k = arrhenius_rate_constant(temperature, k0, E_a) if E_a > 0 and k0 > 0 else k_forward

    if k <= 0:
        return CouplingResult(warnings=warnings)

    # Convert to mol/m^3 (canonical ODE unit)
    active_conc = activated_sites_mol_per_particle / bead_volume
    max_sites_conc = max_coupled_mol_per_particle / bead_volume

    # ODE initial state: [coupled, active, ligand_free]
    y0 = np.array([0.0, active_conc, ligand_concentration])

    sol = solve_ivp(
        fun=_steric_binding_rhs,
        t_span=(0.0, time),
        y0=y0,
        method="Radau",
        args=(k, max_sites_conc),
        rtol=1e-8,
        atol=1e-12,
        max_step=time / 10.0,
    )

    if not sol.success:
        warnings.append(f"ODE solver warning: {sol.message}")

    coupled_final = max(sol.y[0, -1], 0.0)
    active_final = max(sol.y[1, -1], 0.0)
    ligand_final = max(sol.y[2, -1], 0.0)

    sites_coupled = coupled_final * bead_volume
    sites_coupled = min(sites_coupled, max_coupled_mol_per_particle)  # Codex P2-4: hard clamp
    sites_consumed = sites_coupled  # No hydrolysis in steric model

    conversion = sites_consumed / activated_sites_mol_per_particle if activated_sites_mol_per_particle > 0 else 0.0
    conversion = max(0.0, min(conversion, 1.0))

    reagent_rem = ligand_final / ligand_concentration if ligand_concentration > 0 else 0.0

    # Balance: active_consumed = coupled
    active_consumed = active_conc - active_final
    balance_error = abs(active_consumed - coupled_final) / max(active_conc, 1e-30)

    return CouplingResult(
        conversion=conversion,
        sites_consumed=sites_consumed,
        sites_coupled=sites_coupled,
        sites_hydrolyzed=0.0,
        reagent_remaining_fraction=max(0.0, min(reagent_rem, 1.0)),
        hydrolysis_fraction=0.0,
        solver_success=sol.success,
        solver_message=sol.message if not sol.success else "",
        site_balance_error=balance_error,
        warnings=warnings,
    )


def solve_quenching(
    activated_sites_mol_per_particle: float,
    reagent_concentration: float,
    k_forward: float,
    E_a: float,
    temperature: float,
    time: float,
    bead_volume: float,
    stoichiometry: float = 1.0,
    hydrolysis_rate: float = 0.0,
) -> CouplingResult:
    """Solve quenching reaction (Template 1 wrapper — high [reagent], fast blocking).

    Quenching uses the simple second-order ODE but returns sites_blocked
    instead of sites_coupled. No hydrolysis tracking needed (quench reagents
    are typically stable).

    Args:
        activated_sites_mol_per_particle: Remaining activated sites [mol/particle].
        reagent_concentration: Quenching reagent concentration [mol/m^3].
        k_forward: Forward rate constant [m^3/(mol*s)].
        E_a: Activation energy [J/mol].
        temperature: Reaction temperature [K].
        time: Reaction time [s].
        bead_volume: Volume of one bead [m^3/particle].
        stoichiometry: Moles reagent per mole site blocked [-].
        hydrolysis_rate: Reagent hydrolysis rate [1/s] (usually 0 for quench).

    Returns:
        CouplingResult with sites_blocked populated.
    """
    if activated_sites_mol_per_particle <= 0 or reagent_concentration <= 0 or time <= 0:
        return CouplingResult()
    if bead_volume <= 0:
        return CouplingResult(solver_success=False, solver_message="bead_volume <= 0")

    # Convert to mol/m^3
    active_conc = activated_sites_mol_per_particle / bead_volume

    # Arrhenius
    # Back-calculate Arrhenius prefactor; clamp to prevent overflow
    _exp_ref = math.exp(max(-E_a / (R_GAS * 298.15), -700))
    k0 = min(k_forward / max(_exp_ref, 1e-30), 1e30) if E_a > 0 else 0.0

    conversion, reagent_rem = solve_second_order_consumption(
        acs_concentration=active_conc,
        reagent_concentration=reagent_concentration,
        k_forward=k_forward,
        stoichiometry=stoichiometry,
        time=time,
        temperature=temperature,
        E_a=E_a,
        k0=k0,
        hydrolysis_rate=hydrolysis_rate,
    )

    sites_blocked = conversion * activated_sites_mol_per_particle

    return CouplingResult(
        conversion=conversion,
        sites_consumed=sites_blocked,
        sites_blocked=sites_blocked,
        sites_coupled=0.0,
        sites_hydrolyzed=0.0,
        reagent_remaining_fraction=reagent_rem,
        hydrolysis_fraction=0.0,
        solver_success=True,
    )
