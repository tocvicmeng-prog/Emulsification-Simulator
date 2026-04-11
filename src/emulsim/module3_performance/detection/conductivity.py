"""Conductivity detection model for chromatographic monitoring.

Architecture: docs/13_module2_module3_final_implementation_plan.md, Phase E.

Electrical conductivity of an electrolyte solution is proportional to the
concentration and molar conductivity of each ionic species:

    kappa = sum_i( lambda_i * c_i )

where:
    kappa   : Electrical conductivity [S/m]
    lambda_i: Molar conductivity of species i [S*m^2/mol]
    c_i     : Concentration of species i [mol/m^3]

For a 1:1 electrolyte (e.g., NaCl) with one dominant salt:
    kappa = (lambda_cation + lambda_anion) * c_salt

Default molar conductivities at 25 C (infinite dilution):
    Na+ : 50.1e-4 S*m^2/mol
    Cl- : 76.3e-4 S*m^2/mol

These values give kappa in [S/m] when c is in [mol/m^3].

Conductivity is commonly reported as:
    - [mS/cm] in AKTA / NGC chromatography systems
    - Conversion: 1 S/m = 10 mS/cm

References:
    Atkins, P. & de Paula, J. (2010). Physical Chemistry, 9th Ed., Table 20B.1.
"""

from __future__ import annotations

import numpy as np


# Molar conductivities at infinite dilution, 25 C [S*m^2/mol]
_LAMBDA_NA = 50.1e-4   # Na+
_LAMBDA_CL = 76.3e-4   # Cl-
_LAMBDA_K  = 73.5e-4   # K+
_LAMBDA_SO4 = 80.0e-4  # SO4^2- (per equivalent)


def compute_conductivity(
    salt_concentration: float | np.ndarray,
    molar_conductivity_cation: float = _LAMBDA_NA,
    molar_conductivity_anion: float = _LAMBDA_CL,
) -> np.ndarray:
    """Compute solution conductivity from salt concentration.

    kappa = (lambda_cation + lambda_anion) * c_salt

    Args:
        salt_concentration: Salt concentration [mol/m^3]. Accepts scalar or array.
        molar_conductivity_cation: Cation molar conductivity [S*m^2/mol].
            Default 50.1e-4 is Na+ at 25 C.
        molar_conductivity_anion: Anion molar conductivity [S*m^2/mol].
            Default 76.3e-4 is Cl- at 25 C.

    Returns:
        Electrical conductivity [S/m], same shape as salt_concentration.
    """
    c = np.asarray(salt_concentration, dtype=float)
    c_safe = np.maximum(c, 0.0)
    kappa = (molar_conductivity_cation + molar_conductivity_anion) * c_safe
    return kappa


def conductivity_to_ms_per_cm(kappa_si: float | np.ndarray) -> np.ndarray:
    """Convert conductivity from SI units [S/m] to [mS/cm].

    1 S/m = 10 mS/cm

    Args:
        kappa_si: Conductivity [S/m].

    Returns:
        Conductivity [mS/cm].
    """
    return np.asarray(kappa_si, dtype=float) * 10.0


def conductivity_to_nacl_concentration(
    kappa_si: float | np.ndarray,
    molar_conductivity_cation: float = _LAMBDA_NA,
    molar_conductivity_anion: float = _LAMBDA_CL,
) -> np.ndarray:
    """Back-calculate NaCl concentration from measured conductivity.

    Inverse of compute_conductivity for a 1:1 electrolyte.

    Args:
        kappa_si: Measured conductivity [S/m].
        molar_conductivity_cation: Cation molar conductivity [S*m^2/mol].
        molar_conductivity_anion: Anion molar conductivity [S*m^2/mol].

    Returns:
        Salt concentration [mol/m^3].
    """
    kappa = np.asarray(kappa_si, dtype=float)
    lambda_total = molar_conductivity_cation + molar_conductivity_anion
    return np.maximum(kappa / lambda_total, 0.0)


def compute_chromatogram_conductivity(
    salt_profile: np.ndarray,
    molar_conductivity_cation: float = _LAMBDA_NA,
    molar_conductivity_anion: float = _LAMBDA_CL,
    report_units: str = "mS/cm",
) -> np.ndarray:
    """Compute conductivity chromatogram from a salt concentration profile.

    Convenience wrapper that converts units for display.

    Args:
        salt_profile: Salt concentration vs time [mol/m^3], shape (N_t,).
        molar_conductivity_cation: Cation molar conductivity [S*m^2/mol].
        molar_conductivity_anion: Anion molar conductivity [S*m^2/mol].
        report_units: Output units, "S/m" or "mS/cm" (default "mS/cm").

    Returns:
        Conductivity chromatogram in requested units, shape (N_t,).

    Raises:
        ValueError: If report_units is not recognized.
    """
    kappa = compute_conductivity(
        salt_profile, molar_conductivity_cation, molar_conductivity_anion
    )
    if report_units == "S/m":
        return kappa
    elif report_units == "mS/cm":
        return conductivity_to_ms_per_cm(kappa)
    else:
        raise ValueError(
            f"report_units must be 'S/m' or 'mS/cm', got '{report_units}'"
        )
