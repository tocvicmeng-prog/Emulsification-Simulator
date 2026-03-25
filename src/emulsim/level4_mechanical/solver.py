"""Level 4: Mechanical property prediction for double-network microspheres.

Analytical models:
- Additive DN modulus: G_DN = G_agarose + G_chitosan
- Hertz contact for microsphere compression
- Ogston model for chromatographic Kav
"""

from __future__ import annotations

import numpy as np

from ..datatypes import (
    CrosslinkingResult,
    GelationResult,
    MaterialProperties,
    MechanicalResult,
    SimulationParameters,
)


def agarose_modulus(c_agarose: float, G_prefactor: float = 3000.0,
                    G_exponent: float = 2.2) -> float:
    """Agarose gel shear modulus [Pa] from power-law concentration dependence.

    G_agarose = G_prefactor · (c/c_ref)^G_exponent
    where c_ref = 10 kg/m³ (1% w/v).
    """
    c_ref = 10.0  # kg/m³ = 1% w/v
    if c_agarose <= 0:
        return 0.0
    return G_prefactor * (c_agarose / c_ref) ** G_exponent


def double_network_modulus(G_agarose: float, G_chitosan: float) -> float:
    """Additive double-network shear modulus [Pa].

    G_DN = G_agarose + G_chitosan

    Valid for semi-IPN where networks are independently load-bearing.
    """
    return G_agarose + G_chitosan


def effective_youngs_modulus(G: float, nu: float = 0.45) -> float:
    """Effective Young's modulus from shear modulus [Pa].

    E* = 2·G·(1+ν) for incompressible-like hydrogels (ν ≈ 0.45-0.5).
    """
    return 2.0 * G * (1.0 + nu)


def hertz_contact(E_star: float, R: float,
                  delta_max: float = None,
                  n_points: int = 100) -> tuple[np.ndarray, np.ndarray]:
    """Hertz contact model for sphere compression.

    F = (4/3) · E* · R^(1/2) · δ^(3/2)

    Parameters
    ----------
    E_star : float
        Effective Young's modulus [Pa].
    R : float
        Microsphere radius [m].
    delta_max : float, optional
        Maximum indentation [m]. Default: 10% of R.
    n_points : int
        Number of points in the force-displacement curve.

    Returns
    -------
    delta : np.ndarray
        Indentation [m].
    F : np.ndarray
        Force [N].
    """
    if delta_max is None:
        delta_max = 0.1 * R

    delta = np.linspace(0, delta_max, n_points)
    F = (4.0 / 3.0) * E_star * np.sqrt(R) * delta ** 1.5
    return delta, F


def ogston_kav(rh: np.ndarray, r_fiber: float, phi_fiber: float) -> np.ndarray:
    """Ogston partition coefficient for size-exclusion chromatography.

    Kav = exp(-π · (rh + r_f)² · n_f · L_f)

    Simplified Ogston model:
    Kav = exp(-phi_f · ((rh/r_f) + 1)²)

    Parameters
    ----------
    rh : np.ndarray
        Hydrodynamic radius of solute [m].
    r_fiber : float
        Gel fiber radius [m].
    phi_fiber : float
        Fiber volume fraction [-].
    """
    ratio = rh / r_fiber + 1.0
    return np.exp(-phi_fiber * ratio ** 2)


def solve_mechanical(params: SimulationParameters,
                     props: MaterialProperties,
                     gelation: GelationResult,
                     crosslinking: CrosslinkingResult) -> MechanicalResult:
    """Compute mechanical properties of the double-network microsphere.

    Parameters
    ----------
    params : SimulationParameters
    props : MaterialProperties
    gelation : GelationResult
        From Level 2.
    crosslinking : CrosslinkingResult
        From Level 3.
    """
    # Agarose network modulus
    G_agar = agarose_modulus(
        params.formulation.c_agarose,
        props.G_agarose_prefactor,
        props.G_agarose_exponent,
    )

    # Chitosan network modulus from crosslinking
    G_chit = crosslinking.G_chitosan_final

    # Double-network modulus
    G_DN = double_network_modulus(G_agar, G_chit)

    # Effective Young's modulus
    E_star = effective_youngs_modulus(G_DN)

    # Hertz contact curve
    R = gelation.r_grid[-1] + (gelation.r_grid[1] - gelation.r_grid[0]) / 2.0
    delta_arr, F_arr = hertz_contact(E_star, R)

    # Ogston Kav for a range of protein sizes
    rh_arr = np.logspace(np.log10(1e-9), np.log10(50e-9), 50)  # 1-50 nm
    phi_fiber = 1.0 - gelation.porosity  # polymer volume fraction
    Kav_arr = ogston_kav(rh_arr, props.r_fiber, phi_fiber)

    return MechanicalResult(
        G_agarose=float(G_agar),
        G_chitosan=float(G_chit),
        G_DN=float(G_DN),
        E_star=float(E_star),
        delta_array=delta_arr,
        F_array=F_arr,
        rh_array=rh_arr,
        Kav_array=Kav_arr,
        pore_size_mean=float(gelation.pore_size_mean),
        xi_mesh=float(crosslinking.xi_final),
    )
