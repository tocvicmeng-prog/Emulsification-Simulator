"""Thermal physics for the 'adjust cooling rate' suggestion's derivation page.

Inverse problem: given a user's target pore size and the current run's droplet
size + oil-phase properties, what cooling rate would hit the target?

The forward direction (pore ← cooling_rate) is encapsulated in the L2 Avrami
/ Cahn-Hilliard solver; this module re-implements the *scaling* relations
in a closed form so the derivation page can show the inversion step-by-step.

Governing chain (see JCP-EMULSIM-DERIV-001 §3.1):
    1. Biot number:      Bi = h*d/(2*k_oil)
    2. Lumped tau_th:    tau_th = rho*cp*d/(6*h)
    3. |dT/dt|:          dT/dt = (T_oil - T_bath)/tau_th
    4. Spinodal dwell:   dt_spin = 2*dT_band/|dT/dt|
    5. CH wavelength:    lambda* ~ sqrt(M*kappa) * dt_spin^(1/2)
    6. Pore ~ lambda*/2

Inverted step 5-6 at target pore gives the required dwell, then step 4 the
required cooling rate. Steps 1-2 check feasibility (lumped capacitance
valid; cooling rate achievable given the bath + convective coefficient).

Confidence tier:
    - Steps 1-3: VALIDATED (textbook lumped capacitance)
    - Step 4:    SEMI_QUANTITATIVE (scaling; dT_band is matrix-specific)
    - Steps 5-6: SEMI_QUANTITATIVE for l2_mode='ch_2d' or 'ch_ternary'
                 QUALITATIVE_TREND for l2_mode='empirical' (pore is a
                 correlation, not a first-principles derivation)

References:
    Incropera, DeWitt, Bergman, Lavine. Fundamentals of Heat and Mass Transfer.
    Cahn & Hilliard (1958) J. Chem. Phys. 28:258.
    Furukawa (1985) Adv. Phys. 34:703 (spinodal dynamics scaling).
"""

from __future__ import annotations

import math
from dataclasses import dataclass

# Default spinodal temperature band width — how many K around T_gel does
# the system spend "inside" the spinodal regime. Empirical; 2-5 K is typical
# for agarose in oil. Used as a scaling constant in the absence of
# matrix-specific data.
_DT_BAND_DEFAULT_K: float = 3.0

# Proportionality constant pulling the Cahn-Hilliard scaling into an absolute
# number. Calibrated so that at typical conditions (d=100 um, dT/dt=0.2 K/s,
# dT_band=3 K) the predicted pore sits in the 50-150 nm range observed for
# 4% agarose. Order-of-magnitude only.
_CH_SCALING_PREFACTOR: float = 3.0e-8  # [m / s^0.5]

# Tolerance band on the target pore for deriving min/max cooling rates.
_PORE_TOLERANCE: float = 0.10  # +/- 10%


@dataclass(frozen=True)
class CoolingRateTarget:
    """Output of cooling_rate_for_target_pore -- see suggestions.types.TargetRange."""

    nominal: float          # [K/s] best estimate
    min: float              # [K/s] slower bound (corresponds to pore + 10%)
    max: float              # [K/s] faster bound (corresponds to pore - 10%)
    biot_number: float      # dimensionless regime check
    tau_th: float           # [s] lumped thermal time constant
    dT_dt_effective: float  # [K/s] actual cooling rate from current inputs
    limited_by: str         # "ok" | "biot_invalid" | "bath_capacity" | "unachievable"
    confidence_tier: str    # "SEMI_QUANTITATIVE" | "QUALITATIVE_TREND"
    assumptions: tuple[str, ...]


def bead_biot_number(d: float, h_oil: float, k_oil: float) -> float:
    """Biot = h * (d/2) / k_oil. Bi < 0.1 -> lumped capacitance valid."""
    if k_oil <= 0:
        raise ValueError(f"k_oil must be positive, got {k_oil}")
    return h_oil * (d / 2.0) / k_oil


def lumped_cooling_time(d: float, rho_d: float, cp_d: float, h_oil: float) -> float:
    """Characteristic cooling time tau_th = rho*cp*V/(h*A) = rho*cp*d/(6*h).

    For a sphere V/A = d/6. Returns time to reach ~1/e of the initial
    temperature excess when cooled in a well-stirred bath.
    """
    if h_oil <= 0:
        raise ValueError(f"h_oil must be positive, got {h_oil}")
    return rho_d * cp_d * d / (6.0 * h_oil)


def effective_cooling_rate(T_oil: float, T_bath: float, tau_th: float) -> float:
    """|dT/dt| at the bead interior under lumped capacitance, evaluated at t=0."""
    if tau_th <= 0:
        raise ValueError(f"tau_th must be positive, got {tau_th}")
    return abs(T_oil - T_bath) / tau_th


def spinodal_dwell_time(cooling_rate: float, dT_band: float = _DT_BAND_DEFAULT_K) -> float:
    """Time spent inside the spinodal temperature band during linear cooling."""
    if cooling_rate <= 0:
        raise ValueError(f"cooling_rate must be positive, got {cooling_rate}")
    return 2.0 * dT_band / cooling_rate


def pore_from_dwell(dwell_time: float, prefactor: float = _CH_SCALING_PREFACTOR) -> float:
    """Cahn-Hilliard scaling: pore ~ sqrt(M*kappa) * dwell^(1/2).

    The (M*kappa)^(1/2) term is absorbed into `prefactor` for closed-form
    reporting. For fully mechanistic pore, see level2_gelation/spatial.py.
    """
    if dwell_time < 0:
        raise ValueError(f"dwell_time must be non-negative, got {dwell_time}")
    return prefactor * math.sqrt(dwell_time)


def dwell_for_target_pore(target_pore: float, prefactor: float = _CH_SCALING_PREFACTOR) -> float:
    """Inverse of pore_from_dwell."""
    if target_pore <= 0:
        raise ValueError(f"target_pore must be positive, got {target_pore}")
    return (target_pore / prefactor) ** 2


def cooling_rate_for_target_pore(
    *,
    target_pore: float,
    d_bead: float,
    T_oil: float,
    T_bath: float,
    rho_d: float,
    cp_d: float,
    h_oil: float,
    k_oil: float,
    l2_mode: str = "ch_2d",
    dT_band: float = _DT_BAND_DEFAULT_K,
    prefactor: float = _CH_SCALING_PREFACTOR,
    pore_tolerance: float = _PORE_TOLERANCE,
) -> CoolingRateTarget:
    """Derive the required |dT/dt| to land at target_pore, with a tolerance band.

    Parameters
    ----------
    target_pore : [m]
        User-entered target pore size.
    d_bead : [m]
        Bead diameter (usually d50 from L1).
    T_oil, T_bath : [K]
        Dispersed-phase temperature at emulsification and the bath.
    rho_d, cp_d, h_oil, k_oil : [SI]
        Dispersed-phase density, specific heat, convective coefficient
        and oil thermal conductivity.
    l2_mode : str
        'ch_2d' / 'ch_ternary' -> SEMI_QUANTITATIVE.
        'empirical' -> QUALITATIVE_TREND (the caller should refuse to
        display a number; we still compute one for internal use).

    Returns
    -------
    CoolingRateTarget : nominal + min/max + diagnostics.
    """
    # Forward check (current conditions)
    bi = bead_biot_number(d_bead, h_oil, k_oil)
    tau = lumped_cooling_time(d_bead, rho_d, cp_d, h_oil)
    dTdt_current = effective_cooling_rate(T_oil, T_bath, tau)

    # Tier
    tier = "SEMI_QUANTITATIVE" if l2_mode in ("ch_2d", "ch_ternary") else "QUALITATIVE_TREND"

    # Inverse: target pore -> required dwell -> required |dT/dt|
    dwell_target = dwell_for_target_pore(target_pore, prefactor=prefactor)
    rate_nominal = 2.0 * dT_band / dwell_target

    # +/-10% pore tolerance -> asymmetric rate band (rate ~ 1/dwell ~ 1/pore^2)
    pore_min = target_pore * (1.0 - pore_tolerance)
    pore_max = target_pore * (1.0 + pore_tolerance)
    rate_fast = 2.0 * dT_band / dwell_for_target_pore(pore_min, prefactor=prefactor)
    rate_slow = 2.0 * dT_band / dwell_for_target_pore(pore_max, prefactor=prefactor)

    # Feasibility diagnostic
    assumptions: list[str] = [
        f"Lumped capacitance (Bi={bi:.3f} < 0.1 required for validity)",
        f"Spinodal dwell band dT_band={dT_band:.1f} K",
        f"Cahn-Hilliard scaling prefactor {prefactor:.2e} m/s^0.5 (order-of-magnitude)",
    ]
    if bi >= 0.1:
        limited_by = "biot_invalid"
        assumptions.append(
            f"WARN: Bi={bi:.2f} exceeds 0.1 -> lumped capacitance is a rough "
            f"approximation; the internal gradient may dominate."
        )
    elif rate_nominal > dTdt_current * 10.0:
        limited_by = "bath_capacity"
        assumptions.append(
            f"Required rate {rate_nominal:.3f} K/s is >10x the rate achievable "
            f"with current (T_oil, T_bath, h_oil); consider a colder bath or "
            f"higher h_oil (e.g. pre-chilled oil with jacketed cooling)."
        )
    elif rate_nominal < dTdt_current / 10.0:
        limited_by = "operational"
        assumptions.append(
            f"Required rate {rate_nominal:.3f} K/s is <0.1x the current rate; "
            f"slowing the cooling that far typically needs a heated jacket "
            f"ramp rather than passive cooling."
        )
    else:
        limited_by = "ok"

    return CoolingRateTarget(
        nominal=rate_nominal,
        min=rate_slow,
        max=rate_fast,
        biot_number=bi,
        tau_th=tau,
        dT_dt_effective=dTdt_current,
        limited_by=limited_by,
        confidence_tier=tier,
        assumptions=tuple(assumptions),
    )


__all__ = [
    "CoolingRateTarget",
    "bead_biot_number",
    "lumped_cooling_time",
    "effective_cooling_rate",
    "spinodal_dwell_time",
    "pore_from_dwell",
    "dwell_for_target_pore",
    "cooling_rate_for_target_pore",
]
