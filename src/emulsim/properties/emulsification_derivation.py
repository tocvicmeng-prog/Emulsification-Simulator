"""Emulsification physics for increase_rpm / decrease_rpm suggestions.

Forward chain (Kolmogorov regime, stirred-vessel emulsification):
    We = rho_c * N^2 * D^3 / sigma
    d32/D = C_We * We^(-0.6) * (1 + beta * phi_d)                  [Sprow 1967]

Inverse: given target_d32, return required N (rpm/60), bounded by:
    - Kolmogorov microscale: d_min ~ (mu^3 / (rho * eps))^(1/4)
    - Transition to laminar: Re = N * D^2 * rho / mu > ~10^4

References:
    Calderbank (1958) Trans. IChemE 36:443.
    Sprow (1967) Chem. Eng. Sci. 22:435.
    Paul, Atiemo-Obeng, Kresta (2004) Handbook of Industrial Mixing.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

# Sprow d32/D = C_We * We^-0.6 * (1 + beta * phi_d). Values from the
# handbook reference above for standard Rushton impellers in low-viscosity
# aqueous/oil systems.
_C_WE: float = 0.06
_BETA: float = 9.0
_KOLMOGOROV_POWER_NUMBER: float = 5.0  # Rushton impeller Po in turbulent regime


@dataclass(frozen=True)
class RPMTarget:
    """Target RPM derived from a target droplet size."""

    nominal: float          # [rpm]
    min: float              # [rpm] slower bound (d32 + 20%)
    max: float              # [rpm] faster bound (d32 - 20%)
    kolmogorov_floor: float # [m] theoretical smallest droplet at nominal rpm
    reynolds: float         # Re at nominal
    limited_by: str         # "ok" | "kolmogorov" | "laminar"
    confidence_tier: str    # "SEMI_QUANTITATIVE"
    assumptions: tuple[str, ...]


def weber_number(N_rps: float, D: float, rho_c: float, sigma: float) -> float:
    """We = rho * N^2 * D^3 / sigma. N in rev/s, D impeller diameter [m]."""
    if sigma <= 0:
        raise ValueError(f"sigma must be positive, got {sigma}")
    return rho_c * (N_rps ** 2) * (D ** 3) / sigma


def d32_from_weber(We: float, D: float, phi_d: float) -> float:
    """Sprow correlation. d32/D = C_We * We^-0.6 * (1 + beta * phi_d)."""
    if We <= 0:
        raise ValueError(f"We must be positive, got {We}")
    return D * _C_WE * (We ** -0.6) * (1.0 + _BETA * phi_d)


def weber_for_target_d32(target_d32: float, D: float, phi_d: float) -> float:
    """Invert d32_from_weber -> We."""
    if target_d32 <= 0:
        raise ValueError(f"target_d32 must be positive, got {target_d32}")
    ratio = target_d32 / (D * _C_WE * (1.0 + _BETA * phi_d))
    return ratio ** (-1.0 / 0.6)


def kolmogorov_scale(N_rps: float, D: float, rho_c: float, mu_c: float) -> float:
    """eta_K = (mu^3 / (rho * eps))^(1/4); eps = Po * N^3 * D^5 / V (bulk avg)."""
    if mu_c <= 0:
        raise ValueError(f"mu_c must be positive, got {mu_c}")
    # Approximate bulk volume as (2*D)^3 = 8*D^3 (a standard rule-of-thumb
    # for a baffled stirred tank with T=D*T/D_ratio typical ~3).
    V = 8.0 * D ** 3
    eps = _KOLMOGOROV_POWER_NUMBER * rho_c * (N_rps ** 3) * (D ** 5) / V
    if eps <= 0:
        return float("inf")
    return ((mu_c ** 3) / (rho_c * eps)) ** 0.25


def rpm_from_N_rps(N_rps: float) -> float:
    return N_rps * 60.0


def N_rps_from_rpm(rpm: float) -> float:
    return rpm / 60.0


def rpm_for_target_d32(
    *,
    target_d32: float,
    D: float,
    rho_c: float,
    mu_c: float,
    sigma: float,
    phi_d: float,
    d32_tolerance: float = 0.20,
) -> RPMTarget:
    """Derive the required RPM to hit target_d32, with a tolerance band.

    Parameters
    ----------
    target_d32 : [m]
        Sauter mean diameter the user is targeting.
    D : [m]
        Impeller diameter.
    rho_c, mu_c : continuous (oil) phase density and viscosity.
    sigma : [N/m]
        Interfacial tension.
    phi_d : dispersed phase volume fraction.
    d32_tolerance : fractional band (default +/-20%).

    Returns
    -------
    RPMTarget with nominal + band + feasibility diagnostics.
    """
    # Inverse Sprow
    We_nominal = weber_for_target_d32(target_d32, D, phi_d)
    N_nominal_rps = math.sqrt(We_nominal * sigma / (rho_c * D ** 3))

    # Band: smaller d32 requires higher RPM
    d32_small = target_d32 * (1.0 - d32_tolerance)
    d32_large = target_d32 * (1.0 + d32_tolerance)
    We_fast = weber_for_target_d32(d32_small, D, phi_d)
    We_slow = weber_for_target_d32(d32_large, D, phi_d)
    N_fast_rps = math.sqrt(We_fast * sigma / (rho_c * D ** 3))
    N_slow_rps = math.sqrt(We_slow * sigma / (rho_c * D ** 3))

    # Feasibility checks at nominal
    eta_K = kolmogorov_scale(N_nominal_rps, D, rho_c, mu_c)
    Re = rho_c * N_nominal_rps * (D ** 2) / mu_c

    assumptions: list[str] = [
        "Sprow (1967) correlation d32/D = 0.06*We^-0.6*(1+9*phi_d); Kolmogorov regime.",
        "Rushton impeller, baffled tank, power number Po=5 assumed.",
        "Bulk dissipation rate averaged over 8*D^3 vessel volume.",
    ]

    if target_d32 < 3.0 * eta_K:
        limited_by = "kolmogorov"
        assumptions.append(
            f"WARN: target d32 ({target_d32*1e6:.1f} um) is within 3x the "
            f"Kolmogorov microscale eta_K ({eta_K*1e6:.2f} um) at the nominal "
            f"RPM -- breakage is viscosity-limited, further RPM gives "
            f"diminishing returns."
        )
    elif Re < 1e4:
        limited_by = "laminar"
        assumptions.append(
            f"WARN: Re={Re:.0f} is below the 1e4 turbulent transition; the "
            f"Kolmogorov-regime correlation overshoots actual breakage in "
            f"the transitional regime. Consider a smaller impeller (increase N)."
        )
    else:
        limited_by = "ok"

    return RPMTarget(
        nominal=rpm_from_N_rps(N_nominal_rps),
        min=rpm_from_N_rps(N_slow_rps),
        max=rpm_from_N_rps(N_fast_rps),
        kolmogorov_floor=eta_K,
        reynolds=Re,
        limited_by=limited_by,
        confidence_tier="SEMI_QUANTITATIVE",
        assumptions=tuple(assumptions),
    )


__all__ = [
    "RPMTarget",
    "weber_number",
    "d32_from_weber",
    "weber_for_target_d32",
    "kolmogorov_scale",
    "rpm_for_target_d32",
]
