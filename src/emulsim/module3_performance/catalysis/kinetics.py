"""Michaelis-Menten kinetics with effectiveness factor for immobilised enzymes.

Provides rate expressions and Thiele modulus calculations for spherical
catalyst particles in a packed-bed reactor.

Key equations
-------------
Michaelis-Menten:  v = V_max * S / (K_m + S)
Thiele modulus (first-order limit):
    Phi = R * sqrt(V_max / (K_m * D_eff))
Effectiveness factor (spherical, first-order):
    eta = (3/Phi) * (1/tanh(Phi) - 1/Phi)
Generalised Thiele modulus (full Michaelis-Menten):
    Phi_gen = R * sqrt(V_max / (2*D_eff*(K_m*S - K_m^2*ln(1+S/K_m))))
"""

from __future__ import annotations

import numpy as np


# ---------------------------------------------------------------------------
# Rate law
# ---------------------------------------------------------------------------

def michaelis_menten_rate(
    S: np.ndarray | float,
    V_max: float,
    K_m: float,
) -> np.ndarray | float:
    """Michaelis-Menten reaction rate.

    Parameters
    ----------
    S : array_like
        Substrate concentration [mol/m^3].  Negative values are clipped to 0.
    V_max : float
        Maximum reaction rate [mol/(m^3*s)].
    K_m : float
        Michaelis constant [mol/m^3].  Must be > 0.

    Returns
    -------
    v : same shape as *S*
        Reaction rate [mol/(m^3*s)].
    """
    if K_m <= 0:
        raise ValueError(f"K_m must be positive, got {K_m}")
    S_safe = np.maximum(S, 0.0)
    return V_max * S_safe / (K_m + S_safe)


# ---------------------------------------------------------------------------
# Thiele modulus — first-order limit
# ---------------------------------------------------------------------------

def thiele_modulus(
    R: float,
    V_max: float,
    K_m: float,
    D_eff: float,
) -> float:
    """First-order-limit Thiele modulus for a spherical catalyst particle.

    Parameters
    ----------
    R : float
        Particle radius [m].
    V_max : float
        Maximum reaction rate [mol/(m^3*s)].
    K_m : float
        Michaelis constant [mol/m^3].
    D_eff : float
        Effective intra-particle diffusivity [m^2/s].

    Returns
    -------
    phi : float
        Thiele modulus (dimensionless).
    """
    if D_eff <= 0:
        raise ValueError(f"D_eff must be positive, got {D_eff}")
    if K_m <= 0:
        raise ValueError(f"K_m must be positive, got {K_m}")
    return R * np.sqrt(V_max / (K_m * D_eff))


# ---------------------------------------------------------------------------
# Effectiveness factor — spherical geometry, first-order kinetics
# ---------------------------------------------------------------------------

def effectiveness_factor(phi: float | np.ndarray) -> float | np.ndarray:
    """Internal effectiveness factor for spherical geometry (first-order).

    eta = (3/phi) * (1/tanh(phi) - 1/phi)

    Handles limiting cases:
    - phi -> 0: eta -> 1  (no diffusion limitation)
    - phi -> inf: eta -> 3/phi

    Parameters
    ----------
    phi : float or array
        Thiele modulus (dimensionless).

    Returns
    -------
    eta : same shape as *phi*
        Effectiveness factor in [0, 1].
    """
    phi = np.asarray(phi, dtype=float)
    eta = np.ones_like(phi)

    # For very small phi, use Taylor expansion: eta ~ 1 - phi^2/15
    small = phi < 1e-6
    mid = ~small

    if np.any(mid):
        p = phi[mid]
        eta[mid] = (3.0 / p) * (1.0 / np.tanh(p) - 1.0 / p)

    # Clip to [0, 1] for safety
    eta = np.clip(eta, 0.0, 1.0)

    # Return scalar if input was scalar
    if eta.ndim == 0:
        return float(eta)
    return eta


# ---------------------------------------------------------------------------
# Generalised Thiele modulus — full Michaelis-Menten
# ---------------------------------------------------------------------------

def generalized_thiele_modulus(
    R: float,
    V_max: float,
    K_m: float,
    D_eff: float,
    S_bulk: float,
) -> float:
    """Generalised Thiele modulus for Michaelis-Menten kinetics.

    For a spherical particle with Michaelis-Menten kinetics, the generalised
    modulus accounts for the full saturation curve, not just the first-order
    limit.

    Phi_gen = R * sqrt(V_max / (2 * D_eff * F(S_bulk)))

    where F(S) = K_m * S - K_m^2 * ln(1 + S/K_m).

    When S_bulk -> 0 this reduces to the first-order Thiele modulus.

    Parameters
    ----------
    R : float
        Particle radius [m].
    V_max : float
        Maximum reaction rate [mol/(m^3*s)].
    K_m : float
        Michaelis constant [mol/m^3].
    D_eff : float
        Effective intra-particle diffusivity [m^2/s].
    S_bulk : float
        Bulk substrate concentration [mol/m^3].

    Returns
    -------
    phi_gen : float
        Generalised Thiele modulus (dimensionless).
    """
    if D_eff <= 0:
        raise ValueError(f"D_eff must be positive, got {D_eff}")
    if K_m <= 0:
        raise ValueError(f"K_m must be positive, got {K_m}")

    S_bulk = max(S_bulk, 1e-30)  # avoid log(1+0) = 0 in denominator

    F = K_m * S_bulk - K_m**2 * np.log(1.0 + S_bulk / K_m)

    # F -> S_bulk^2 / (2*K_m) for S_bulk << K_m, always positive
    F = max(F, 1e-30)

    phi_gen = R * np.sqrt(V_max / (2.0 * D_eff * F))
    return float(phi_gen)
