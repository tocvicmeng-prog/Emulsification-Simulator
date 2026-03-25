"""Interfacial tension models for Span-80 stabilized W/O emulsions."""

from __future__ import annotations

import numpy as np

# Constants
R_GAS = 8.314  # J/(mol·K)
M_SPAN80 = 428.6  # g/mol, molecular weight of Span-80


def sigma_bare(T: float) -> float:
    """Bare interfacial tension [N/m] for water/paraffin oil without surfactant.

    sigma_0(T) ~ 0.050 - 0.0001*(T - 293.15)

    Approximate linear decrease with temperature.
    """
    T_ref = 293.15  # 20 degC
    sigma_ref = 0.050  # N/m at 20 degC
    dsigma_dT = -1.0e-4    # N/(m·K)
    return max(sigma_ref + dsigma_dT * (T - T_ref), 0.001)


def interfacial_tension_span80(T: float, c_span80: float) -> float:
    """Interfacial tension [N/m] for Span-80 at oil/water interface.

    Szyszkowski-Langmuir model:
        sigma = sigma_0(T) - R*T*Gamma_inf*ln(1 + K_L*c_mol)

    Parameters
    ----------
    T : float
        Temperature [K].
    c_span80 : float
        Span-80 concentration in oil phase [kg/m3].

    Returns
    -------
    float
        Interfacial tension [N/m]. Clipped to minimum of 1e-4 N/m.
    """
    # Convert to molar concentration
    c_mol = c_span80 / M_SPAN80 * 1000.0  # mol/m3

    sigma_0 = sigma_bare(T)

    # Langmuir isotherm parameters for Span-80
    Gamma_inf = 3.5e-6   # mol/m2 (maximum surface excess)
    K_L = 0.75           # m3/mol (calibrated: ~5 mN/m at 2% w/v Span-80, 90°C)

    # Szyszkowski equation
    sigma = sigma_0 - R_GAS * T * Gamma_inf * np.log(1.0 + K_L * c_mol)

    # Physical lower bound
    return max(sigma, 1e-4)


def dynamic_interfacial_tension(T: float, c_span80: float, t: float,
                                 D_span80: float = 1e-10) -> float:
    """Time-dependent interfacial tension accounting for Span-80 adsorption kinetics.

    Simplified Ward-Tordai model for diffusion-controlled adsorption.
    At long times, approaches equilibrium value from interfacial_tension_span80.

    Parameters
    ----------
    T : float
        Temperature [K].
    c_span80 : float
        Bulk Span-80 concentration [kg/m3].
    t : float
        Time since interface creation [s].
    D_span80 : float
        Diffusion coefficient of Span-80 in oil [m2/s].

    Returns
    -------
    float
        Dynamic interfacial tension [N/m].
    """
    sigma_eq = interfacial_tension_span80(T, c_span80)
    sigma_0 = sigma_bare(T)

    # Characteristic adsorption time
    c_mol = c_span80 / M_SPAN80 * 1000.0
    Gamma_inf = 3.5e-6
    if c_mol > 0:
        tau_ads = (Gamma_inf ** 2) / (D_span80 * c_mol ** 2)
    else:
        return sigma_0

    # Exponential relaxation (simplified Ward-Tordai)
    sigma = sigma_eq + (sigma_0 - sigma_eq) * np.exp(-t / tau_ads)
    return sigma
