"""Thermodynamic models for polymer solution free energy and phase behavior."""

from __future__ import annotations

import numpy as np

# Constants
R_GAS = 8.314       # J/(mol·K)
K_BOLTZMANN = 1.38e-23  # J/K
N_AVOGADRO = 6.022e23


def chi_flory_huggins(T: float, A: float = 150.0, B: float = 0.0) -> float:
    """Flory-Huggins interaction parameter χ(T).

    χ(T) = A/T + B

    For agarose-water: χ ≈ 0.497 at 25°C.
    A ≈ 150 K, B ≈ -0.006 gives χ(298) ≈ 0.497.

    Parameters
    ----------
    T : float
        Temperature [K].
    A : float
        Enthalpic contribution [K].
    B : float
        Entropic contribution [-].
    """
    return A / T + B


def flory_huggins_free_energy(phi: float | np.ndarray, T: float, N_p: float = 100.0,
                               chi: float | None = None, A: float = 150.0,
                               B: float = 0.0, v0: float = 1.8e-29) -> float | np.ndarray:
    """Flory-Huggins free energy density [J/m³].

    f(φ) = (k_B·T/v₀)·[φ·ln(φ)/N_p + (1-φ)·ln(1-φ) + χ·φ·(1-φ)]

    Parameters
    ----------
    phi : float or array
        Polymer volume fraction (0 < φ < 1).
    T : float
        Temperature [K].
    N_p : float
        Degree of polymerization of polymer.
    chi : float, optional
        Flory-Huggins χ. If None, computed from A, B, T.
    A, B : float
        χ(T) = A/T + B parameters.
    v0 : float
        Reference volume per lattice site [m³] (~molar volume of water / N_A).
    """
    if chi is None:
        chi = chi_flory_huggins(T, A, B)

    phi = np.asarray(phi, dtype=float)
    # Clip to avoid log(0)
    phi_c = np.clip(phi, 1e-15, 1.0 - 1e-15)

    prefactor = K_BOLTZMANN * T / v0
    f = prefactor * (
        phi_c * np.log(phi_c) / N_p
        + (1.0 - phi_c) * np.log(1.0 - phi_c)
        + chi * phi_c * (1.0 - phi_c)
    )
    return f


def flory_huggins_derivative(phi: float | np.ndarray, T: float, N_p: float = 100.0,
                              chi: float | None = None, A: float = 150.0,
                              B: float = 0.0, v0: float = 1.8e-29) -> float | np.ndarray:
    """First derivative of Flory-Huggins free energy df/dφ [J/m³].

    f'(φ) = (k_B·T/v₀)·[ln(φ)/N_p + 1/N_p - ln(1-φ) - 1 + χ·(1 - 2φ)]
    """
    if chi is None:
        chi = chi_flory_huggins(T, A, B)

    phi = np.asarray(phi, dtype=float)
    phi_c = np.clip(phi, 1e-15, 1.0 - 1e-15)

    prefactor = K_BOLTZMANN * T / v0
    return prefactor * (
        np.log(phi_c) / N_p + 1.0 / N_p
        - np.log(1.0 - phi_c) - 1.0
        + chi * (1.0 - 2.0 * phi_c)
    )


def flory_huggins_second_derivative(phi: float | np.ndarray, T: float, N_p: float = 100.0,
                                     chi: float | None = None, A: float = 150.0,
                                     B: float = 0.0, v0: float = 1.8e-29) -> float | np.ndarray:
    """Second derivative of Flory-Huggins free energy d²f/dφ² [J/m³].

    f''(φ) = (k_B·T/v₀)·[1/(N_p·φ) + 1/(1-φ) - 2χ]

    Negative in the spinodal region (unstable).
    """
    if chi is None:
        chi = chi_flory_huggins(T, A, B)

    phi = np.asarray(phi, dtype=float)
    phi_c = np.clip(phi, 1e-15, 1.0 - 1e-15)

    prefactor = K_BOLTZMANN * T / v0
    return prefactor * (
        1.0 / (N_p * phi_c)
        + 1.0 / (1.0 - phi_c)
        - 2.0 * chi
    )


def spinodal_boundary(T: float, N_p: float = 100.0, A: float = 150.0,
                      B: float = 0.0) -> tuple[float, float]:
    """Find spinodal boundary φ values where f''(φ) = 0.

    Returns (φ_s1, φ_s2) where φ_s1 < φ_s2, or (nan, nan) if no spinodal.
    """
    chi = chi_flory_huggins(T, A, B)

    # f''(φ) = 0  →  1/(N_p·φ) + 1/(1-φ) = 2χ
    # Solve numerically
    from scipy.optimize import brentq

    def f2(phi):
        return 1.0 / (N_p * phi) + 1.0 / (1.0 - phi) - 2.0 * chi

    # Check if spinodal exists: f''(φ) must be negative somewhere
    # Minimum of 1/(Np·φ) + 1/(1-φ) occurs at φ* = 1/(1 + sqrt(Np))
    phi_min = 1.0 / (1.0 + np.sqrt(N_p))
    if f2(phi_min) > 0:
        return (np.nan, np.nan)

    try:
        phi_s1 = brentq(f2, 1e-10, phi_min)
        phi_s2 = brentq(f2, phi_min, 1.0 - 1e-10)
        return (phi_s1, phi_s2)
    except ValueError:
        return (np.nan, np.nan)


def flory_rehner_crosslink_density(nu_2s: float, chi: float,
                                    V1: float = 1.8e-5) -> float:
    """Effective crosslink density from Flory-Rehner swelling theory [mol/m³].

    n_e = -[ln(1-ν₂,s) + ν₂,s + χ·ν₂,s²] / [V₁·(ν₂,s^(1/3) - ν₂,s/2)]

    Parameters
    ----------
    nu_2s : float
        Polymer volume fraction in swollen gel.
    chi : float
        Flory-Huggins interaction parameter.
    V1 : float
        Molar volume of solvent [m³/mol]. Water: 1.8e-5.
    """
    numerator = -(np.log(1.0 - nu_2s) + nu_2s + chi * nu_2s**2)
    denominator = V1 * (nu_2s**(1.0/3.0) - nu_2s / 2.0)
    if abs(denominator) < 1e-30:
        return 0.0
    return numerator / denominator


def mesh_size_canal_peppas(nu_2s: float, Mc: float) -> float:
    """Mesh size from Canal-Peppas equation [m].

    ξ = 0.071 · ν₂,s^(-1/3) · M_c^(1/2) [nm → converted to m]

    Parameters
    ----------
    nu_2s : float
        Polymer volume fraction in swollen state.
    Mc : float
        Molecular weight between crosslinks [g/mol].
    """
    if nu_2s <= 0 or Mc <= 0:
        return 0.0
    xi_nm = 0.071 * nu_2s**(-1.0/3.0) * np.sqrt(Mc)
    return xi_nm * 1e-9  # convert nm to m
