"""Flory-Huggins free energy and chemical potential for Cahn-Hilliard solver."""

from __future__ import annotations
import numpy as np

K_BOLTZMANN = 1.38e-23  # J/K


def flory_huggins_mu(phi: np.ndarray, T: float, chi: float,
                     N_p: float = 100.0, v0: float = 1.8e-29) -> np.ndarray:
    """Chemical potential mu = df/dphi (without gradient term).

    mu = (kT/v0) [ln(phi)/Np + 1/Np - ln(1-phi) - 1 + chi*(1-2*phi)]
    """
    phi_c = np.clip(phi, 1e-12, 1.0 - 1e-12)
    pf = K_BOLTZMANN * T / v0
    return pf * (np.log(phi_c) / N_p + 1.0 / N_p - np.log(1.0 - phi_c) - 1.0 + chi * (1.0 - 2.0 * phi_c))


def flory_huggins_d2f(phi: np.ndarray, T: float, chi: float,
                      N_p: float = 100.0, v0: float = 1.8e-29) -> np.ndarray:
    """Second derivative d^2f/dphi^2.

    f'' = (kT/v0) [1/(Np*phi) + 1/(1-phi) - 2*chi]

    Negative in the spinodal region.
    """
    phi_c = np.clip(phi, 1e-12, 1.0 - 1e-12)
    pf = K_BOLTZMANN * T / v0
    return pf * (1.0 / (N_p * phi_c) + 1.0 / (1.0 - phi_c) - 2.0 * chi)


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
