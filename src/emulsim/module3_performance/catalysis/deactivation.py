"""First-order enzyme deactivation models.

Provides time-dependent activity decay for immobilised enzymes, including
temperature dependence via an Arrhenius relation.

Key equations
-------------
Activity:   a(t) = exp(-k_d * t)
Half-life:  t_1/2 = ln(2) / k_d
Arrhenius:  k_d(T) = k_d_ref * exp(E_a_d/R * (1/T_ref - 1/T))
"""

from __future__ import annotations

import numpy as np

# Universal gas constant [J/(mol*K)]
_R_GAS = 8.314


# ---------------------------------------------------------------------------
# Activity decay
# ---------------------------------------------------------------------------

def first_order_deactivation(
    t: float | np.ndarray,
    k_d: float,
) -> float | np.ndarray:
    """First-order enzyme deactivation.

    Parameters
    ----------
    t : float or array
        Time [s].
    k_d : float
        Deactivation rate constant [1/s].  Must be >= 0.

    Returns
    -------
    a : same shape as *t*
        Activity fraction in [0, 1].
    """
    if k_d < 0:
        raise ValueError(f"k_d must be non-negative, got {k_d}")
    return np.exp(-k_d * np.asarray(t, dtype=float))


# ---------------------------------------------------------------------------
# Half-life
# ---------------------------------------------------------------------------

def half_life(k_d: float) -> float:
    """Compute the half-life from a first-order deactivation rate constant.

    Parameters
    ----------
    k_d : float
        Deactivation rate constant [1/s].  Must be > 0.

    Returns
    -------
    t_half : float
        Half-life [s].
    """
    if k_d <= 0:
        raise ValueError(f"k_d must be positive for finite half-life, got {k_d}")
    return np.log(2.0) / k_d


# ---------------------------------------------------------------------------
# Arrhenius temperature dependence
# ---------------------------------------------------------------------------

def arrhenius_deactivation_rate(
    T: float,
    k_d_ref: float,
    E_a_d: float,
    T_ref: float = 310.15,
) -> float:
    """Temperature-dependent deactivation rate via Arrhenius relation.

    k_d(T) = k_d_ref * exp(E_a_d / R * (1/T_ref - 1/T))

    Parameters
    ----------
    T : float
        Temperature [K].
    k_d_ref : float
        Reference deactivation rate constant at *T_ref* [1/s].
    E_a_d : float
        Activation energy for deactivation [J/mol].
    T_ref : float, optional
        Reference temperature [K].  Default 310.15 K (37 deg C).

    Returns
    -------
    k_d : float
        Deactivation rate constant at temperature *T* [1/s].
    """
    if T <= 0:
        raise ValueError(f"Temperature must be positive, got {T}")
    if k_d_ref < 0:
        raise ValueError(f"k_d_ref must be non-negative, got {k_d_ref}")

    return k_d_ref * np.exp(E_a_d / _R_GAS * (1.0 / T_ref - 1.0 / T))
