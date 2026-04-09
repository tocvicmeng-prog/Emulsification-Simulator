"""Ternary Flory-Huggins free energy for agarose-chitosan-water system.

Provides the free energy density and chemical potentials for a ternary
system with two conserved order parameters:
    phi_a = agarose volume fraction
    phi_c = chitosan volume fraction
    phi_s = 1 - phi_a - phi_c = solvent (water) volume fraction

The Flory-Huggins free energy density is:
    f = (kT/v0) * [phi_a*ln(phi_a)/N_a + phi_c*ln(phi_c)/N_c
                    + phi_s*ln(phi_s)
                    + chi_as*phi_a*phi_s + chi_cs*phi_c*phi_s
                    + chi_ac*phi_a*phi_c]
"""

from __future__ import annotations

import numpy as np


def ternary_flory_huggins(phi_a: np.ndarray, phi_c: np.ndarray,
                          T: float, kT_over_v0: float = 4.114e6,
                          N_a: float = 10.0, N_c: float = 5.0,
                          chi_as: float = 0.5, chi_cs: float = 0.6,
                          chi_ac: float = 1.0) -> np.ndarray:
    """Ternary Flory-Huggins free energy density [J/m^3].

    Parameters
    ----------
    phi_a, phi_c : np.ndarray
        Agarose and chitosan volume fractions.
    T : float
        Temperature [K] (unused; chi values passed directly).
    kT_over_v0 : float
        Thermal energy / monomer volume [J/m^3].
    N_a, N_c : float
        Degree of polymerisation for agarose and chitosan.
    chi_as, chi_cs : float
        Polymer-solvent Flory-Huggins interaction parameters.
    chi_ac : float
        Polymer-polymer interaction parameter.

    Returns
    -------
    np.ndarray
        Free energy density at each grid point [J/m^3].
    """
    eps = 1e-15
    phi_a = np.clip(phi_a, eps, 1.0 - eps)
    phi_c = np.clip(phi_c, eps, 1.0 - eps)
    phi_s = np.clip(1.0 - phi_a - phi_c, eps, 1.0 - eps)

    f = kT_over_v0 * (
        phi_a * np.log(phi_a) / N_a
        + phi_c * np.log(phi_c) / N_c
        + phi_s * np.log(phi_s)
        + chi_as * phi_a * phi_s
        + chi_cs * phi_c * phi_s
        + chi_ac * phi_a * phi_c
    )
    return f


def chemical_potential_a(phi_a: np.ndarray, phi_c: np.ndarray,
                         T: float, kT_over_v0: float = 4.114e6,
                         N_a: float = 10.0, N_c: float = 5.0,
                         chi_as: float = 0.5, chi_cs: float = 0.6,
                         chi_ac: float = 1.0) -> np.ndarray:
    r"""Chemical potential mu_a = df/d(phi_a) [J/m^3].

    Derived analytically from the ternary FH free energy with
    phi_s = 1 - phi_a - phi_c:

        mu_a = (kT/v0) * [(ln(phi_a)+1)/N_a - (ln(phi_s)+1)
                           + chi_as*(phi_s - phi_a)
                           - chi_cs*phi_c
                           + chi_ac*phi_c]
    """
    eps = 1e-15
    pa = np.clip(phi_a, eps, 1.0 - eps)
    pc = np.clip(phi_c, eps, 1.0 - eps)
    ps = np.clip(1.0 - pa - pc, eps, 1.0 - eps)

    mu = kT_over_v0 * (
        (np.log(pa) + 1.0) / N_a
        - (np.log(ps) + 1.0)
        + chi_as * (ps - pa)
        - chi_cs * pc
        + chi_ac * pc
    )
    return mu


def chemical_potential_c(phi_a: np.ndarray, phi_c: np.ndarray,
                         T: float, kT_over_v0: float = 4.114e6,
                         N_a: float = 10.0, N_c: float = 5.0,
                         chi_as: float = 0.5, chi_cs: float = 0.6,
                         chi_ac: float = 1.0) -> np.ndarray:
    r"""Chemical potential mu_c = df/d(phi_c) [J/m^3].

    Derived analytically from the ternary FH free energy with
    phi_s = 1 - phi_a - phi_c:

        mu_c = (kT/v0) * [(ln(phi_c)+1)/N_c - (ln(phi_s)+1)
                           + chi_cs*(phi_s - phi_c)
                           - chi_as*phi_a
                           + chi_ac*phi_a]
    """
    eps = 1e-15
    pa = np.clip(phi_a, eps, 1.0 - eps)
    pc = np.clip(phi_c, eps, 1.0 - eps)
    ps = np.clip(1.0 - pa - pc, eps, 1.0 - eps)

    mu = kT_over_v0 * (
        (np.log(pc) + 1.0) / N_c
        - (np.log(ps) + 1.0)
        + chi_cs * (ps - pc)
        - chi_as * pa
        + chi_ac * pa
    )
    return mu
