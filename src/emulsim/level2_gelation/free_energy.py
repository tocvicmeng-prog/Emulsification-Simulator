"""Flory-Huggins free energy and chemical potential for Cahn-Hilliard solver.

Core Flory-Huggins functions are delegated to :mod:`emulsim.properties.thermodynamic`
(single source of truth). This module re-exports thin wrappers that accept ``chi``
as a required positional argument (the solver always pre-computes chi) and adds the
Eyre-splitting helper :func:`contractive_constant`.
"""

from __future__ import annotations
import numpy as np

from ..properties.thermodynamic import (
    flory_huggins_derivative as _fh_derivative,
    flory_huggins_second_derivative as _fh_second_derivative,
)


def flory_huggins_mu(phi: np.ndarray, T: float, chi: float,
                     N_p: float = 100.0, v0: float = 1.8e-29) -> np.ndarray:
    """Chemical potential mu = df/dphi (without gradient term).

    Delegates to :func:`emulsim.properties.thermodynamic.flory_huggins_derivative`.
    """
    return _fh_derivative(phi, T, N_p=N_p, chi=chi, v0=v0)


def flory_huggins_d2f(phi: np.ndarray, T: float, chi: float,
                      N_p: float = 100.0, v0: float = 1.8e-29) -> np.ndarray:
    """Second derivative d^2f/dphi^2.

    Delegates to :func:`emulsim.properties.thermodynamic.flory_huggins_second_derivative`.
    """
    return _fh_second_derivative(phi, T, N_p=N_p, chi=chi, v0=v0)


def contractive_constant(phi_range: tuple, T: float, chi: float,
                         N_p: float = 100.0, v0: float = 1.8e-29) -> float:
    """Compute the contractivity constant C for Eyre splitting.

    C must satisfy C > max|f''(phi)| over the expected phi range.
    The contractive (implicit) part is f_c'(phi) = -C*phi.
    The expansive (explicit) part is f_e'(phi) = mu(phi) + C*phi.
    """
    phi_test = np.linspace(max(phi_range[0], 0.01), min(phi_range[1], 0.99), 200)
    d2f = flory_huggins_d2f(phi_test, T, chi, N_p, v0)
    return 1.2 * np.max(np.abs(d2f))  # 20% safety margin
