"""Dispersed phase viscosity models for agarose/chitosan polymer solutions."""

from __future__ import annotations

import numpy as np

# Physical constants
R_GAS = 8.314  # J/(mol·K)


def intrinsic_viscosity_agarose(M_w: float = 120000.0) -> float:
    """Intrinsic viscosity of agarose from Mark-Houwink relation.

    [eta] = K * M^a

    From Zhao et al. (2020): K = 0.07, a = 0.72
    Units: [eta] in mL/g when M in g/mol.

    Parameters
    ----------
    M_w : float
        Weight-average molecular weight [g/mol]. Default 120,000 for standard agarose.

    Returns
    -------
    float
        Intrinsic viscosity [mL/g].
    """
    K = 0.07
    a = 0.72
    return K * M_w ** a


def huggins_viscosity(eta_intrinsic: float, c: float, k_H: float = 0.4) -> float:
    """Solution specific viscosity from Huggins equation.

    eta_sp/c = [eta] + k_H*[eta]^2*c
    eta_sp = eta/eta_s - 1

    Valid for c*[eta] < ~1 (dilute/semi-dilute).

    Parameters
    ----------
    eta_intrinsic : float
        Intrinsic viscosity [mL/g].
    c : float
        Polymer concentration [g/mL].
    k_H : float
        Huggins coefficient (0.3-0.7 for most polymers).

    Returns
    -------
    float
        Specific viscosity eta_sp (dimensionless).
    """
    return c * (eta_intrinsic + k_H * eta_intrinsic**2 * c)


def martin_viscosity(eta_intrinsic: float, c: float, k_M: float = 0.28) -> float:
    """Solution specific viscosity from Martin equation.

    ln(eta_sp/c) = ln[eta] + k_M*[eta]*c

    Better for concentrated solutions (c*[eta] > 1).

    Parameters
    ----------
    eta_intrinsic : float
        Intrinsic viscosity [mL/g].
    c : float
        Polymer concentration [g/mL].
    k_M : float
        Martin coefficient.

    Returns
    -------
    float
        Specific viscosity eta_sp (dimensionless).
    """
    log_eta_sp_over_c = np.log(eta_intrinsic) + k_M * eta_intrinsic * c
    return c * np.exp(log_eta_sp_over_c)


def cross_model_correction(mu_0: float, shear_rate: float,
                           lambda_cross: float = 0.001,
                           m_cross: float = 0.6) -> float:
    """Cross model for shear-thinning viscosity.

    mu(gamma_dot) = mu_0 / (1 + (lambda * gamma_dot)^m)

    Parameters
    ----------
    mu_0 : float
        Zero-shear viscosity [Pa*s].
    shear_rate : float
        Shear rate [1/s]. Use gap_shear_rate from energy.py.
    lambda_cross : float
        Cross time constant [s].
    m_cross : float
        Cross exponent (0.5-0.9 for polymer solutions).

    Returns
    -------
    float
        Shear-rate-corrected viscosity [Pa*s].
    """
    return mu_0 / (1.0 + (lambda_cross * shear_rate) ** m_cross)


def water_viscosity(T: float) -> float:
    """Dynamic viscosity of water [Pa·s] at temperature T [K].

    VFT-type correlation valid 273-373 K.
    """
    # Simplified Arrhenius fit
    A = 2.414e-5  # Pa·s
    B = 247.8     # K
    C = 140.0     # K
    return A * 10.0 ** (B / (T - C))


def solution_viscosity(T: float, c_agarose: float, c_chitosan: float,
                       M_w_agarose: float = 120000.0) -> float:
    """Total dispersed phase viscosity [Pa·s] for agarose + chitosan solution.

    Uses Mark-Houwink for intrinsic viscosity, Martin equation for
    concentration dependence, and logarithmic mixing rule for blend.
    Temperature effect through solvent viscosity.

    Parameters
    ----------
    T : float
        Temperature [K]. Must be above gelation temperature.
    c_agarose : float
        Agarose concentration [kg/m3].
    c_chitosan : float
        Chitosan concentration [kg/m3].
    M_w_agarose : float
        Agarose molecular weight [g/mol].

    Returns
    -------
    float
        Dynamic viscosity [Pa·s].
    """
    mu_water = water_viscosity(T)

    # Convert kg/m3 to g/mL for viscosity equations
    c_agar_g_ml = c_agarose / 1000.0
    c_chit_g_ml = c_chitosan / 1000.0

    # Agarose contribution
    eta_intr_agar = intrinsic_viscosity_agarose(M_w_agarose)
    c_eta_agar = c_agar_g_ml * eta_intr_agar

    if c_eta_agar > 1.0:
        # Concentrated: use Martin equation
        eta_sp_agar = martin_viscosity(eta_intr_agar, c_agar_g_ml)
    else:
        # Dilute: use Huggins equation
        eta_sp_agar = huggins_viscosity(eta_intr_agar, c_agar_g_ml)

    mu_agarose = mu_water * (1.0 + eta_sp_agar)

    # Chitosan contribution
    # Intrinsic viscosity of chitosan ~5-15 dL/g (500-1500 mL/g) depending on MW and DD
    # Use typical value for high-DD chitosan (~300 kDa): [eta] ~ 800 mL/g
    eta_intr_chit = 800.0  # mL/g (typical high-MW, high-DD chitosan)
    c_eta_chit = c_chit_g_ml * eta_intr_chit

    if c_eta_chit > 1.0:
        eta_sp_chit = martin_viscosity(eta_intr_chit, c_chit_g_ml)
    else:
        eta_sp_chit = huggins_viscosity(eta_intr_chit, c_chit_g_ml)

    mu_chitosan = mu_water * (1.0 + eta_sp_chit)

    # Logarithmic mixing rule for polymer blend viscosity
    total = c_agarose + c_chitosan
    if total < 1e-10:
        return mu_water

    w_agar = c_agarose / total
    w_chit = c_chitosan / total

    log_mu = w_agar * np.log(mu_agarose) + w_chit * np.log(mu_chitosan)
    return np.exp(log_mu)
