"""Avrami gelation kinetics and mobility arrest model."""

from __future__ import annotations
import numpy as np


def avrami_gelation(t_cool: float, k_gel: float, n_avrami: float) -> float:
    """Avrami gelation degree alpha(t).

    alpha = 1 - exp(-(k_gel * t_cool)^n)

    Parameters
    ----------
    t_cool : float
        Time since temperature crossed T_gel [s].
    k_gel : float
        Avrami rate constant [1/s].
    n_avrami : float
        Avrami exponent (2-3 for agarose).

    Returns
    -------
    float
        Gelation degree alpha in [0, 1].
    """
    if t_cool <= 0:
        return 0.0
    arg = (k_gel * t_cool) ** n_avrami
    return 1.0 - np.exp(-arg)


def gelation_rate_constant(T: float, T_gel: float, k_gel_0: float,
                           m: float = 1.0) -> float:
    """Temperature-dependent Avrami rate constant.

    k_gel(T) = k_gel_0 * (1 - T/T_gel)^m   for T < T_gel
             = 0                              for T >= T_gel

    Larger undercooling -> faster gelation.
    """
    if T >= T_gel:
        return 0.0
    return k_gel_0 * ((1.0 - T / T_gel) ** m)


def mobility(phi: np.ndarray, alpha: float, M_0: float,
             beta: float = 2.5) -> np.ndarray:
    """Composition- and gelation-dependent mobility for Cahn-Hilliard.

    M(phi, alpha) = M_0 * phi * (1-phi) * max(1-alpha, 0)^beta

    The phi*(1-phi) factor ensures mobility vanishes at pure phases.
    The (1-alpha)^beta factor arrests coarsening as gelation proceeds.
    """
    phi_c = np.clip(phi, 0.0, 1.0)
    gel_factor = max(1.0 - alpha, 0.0) ** beta
    return M_0 * phi_c * (1.0 - phi_c) * gel_factor


def cooling_temperature(t: float, T_init: float, T_ambient: float,
                        cooling_rate: float) -> float:
    """Bulk emulsion temperature during cooling [K].

    Linear cooling model: T(t) = T_init - cooling_rate * t
    Bounded below by T_ambient.

    Parameters
    ----------
    t : float
        Time since cooling started [s].
    T_init : float
        Initial temperature [K] (= T_oil).
    T_ambient : float
        Ambient/bath temperature [K]. Typically ~293 K (20 C).
    cooling_rate : float
        Cooling rate [K/s].
    """
    T = T_init - cooling_rate * t
    return max(T, T_ambient)
