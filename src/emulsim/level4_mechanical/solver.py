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


def double_network_modulus(G_agarose: float, G_chitosan: float,
                           eta_coupling: float = -0.15) -> float:
    """Phenomenological IPN modulus with coupling term [Pa].

    G_DN = G1 + G2 + eta * sqrt(G1 * G2)

    eta > 0: synergistic (DN gels with sacrificial bonds)
    eta < 0: antagonistic (mutual swelling constraint)
    eta = 0: simple additivity (upper bound)

    Default eta = -0.15 for sequential IPN without sacrificial bonds.
    """
    G_cross = eta_coupling * np.sqrt(max(G_agarose, 0) * max(G_chitosan, 0))
    return max(G_agarose + G_chitosan + G_cross, 0.0)


def ionic_gel_modulus(G_agarose: float, G_ionic: float,
                      f_ionic_strength: float = 0.3) -> float:
    """Modulus estimate for ionic (reversible) crosslinked gels [Pa].

    Ionic bonds are weaker and reversible. The ionic network provides
    only partial reinforcement, modeled as reduced additive contribution.
    """
    return max(G_agarose + f_ionic_strength * G_ionic, 0.0)


def triple_network_modulus(G_agarose: float, G_chitosan_inherent: float,
                            G_independent: float, eta_12: float = -0.15,
                            eta_13: float = 0.0) -> float:
    """Modulus estimate for triple-network (agarose + chitosan + PEGDA) [Pa].

    The independent network (PEGDA) adds linearly with no coupling to
    the agarose-chitosan IPN. eta_13 = 0 by default (no interaction).
    """
    G_ipn = double_network_modulus(G_agarose, G_chitosan_inherent, eta_12)
    G_cross_13 = eta_13 * np.sqrt(max(G_ipn, 0) * max(G_independent, 0))
    return max(G_ipn + G_independent + G_cross_13, 0.0)


def select_modulus_model(G_agarose: float, G_xlink: float,
                          network_metadata=None,
                          eta_coupling: float = -0.15) -> float:
    """Route to the appropriate modulus model based on network metadata.

    Falls back to phenomenological DN modulus if no metadata provided.
    """
    if network_metadata is None:
        return double_network_modulus(G_agarose, G_xlink, eta_coupling)

    family = getattr(network_metadata, 'solver_family', 'amine_covalent')

    if family == "ionic_reversible":
        return ionic_gel_modulus(G_agarose, G_xlink)
    elif family == "independent_network":
        # PEGDA forms a third network; agarose modulus is already the base,
        # G_xlink here is the independent PEGDA network modulus.
        # Use triple_network with zero chitosan contribution.
        return triple_network_modulus(G_agarose, 0.0, G_xlink, eta_13=0.0)
    elif family == "hydroxyl_covalent":
        # Hydroxyl crosslinkers bridge BOTH networks (synergistic coupling)
        eta_syn = getattr(network_metadata, 'eta_coupling_recommended', +0.05)
        if hasattr(network_metadata, 'eta_coupling_recommended'):
            # Not available on NetworkTypeMetadata, use from CrosslinkerProfile
            pass
        return double_network_modulus(G_agarose, G_xlink, eta_coupling=+0.05)
    else:
        # amine_covalent and default
        return double_network_modulus(G_agarose, G_xlink, eta_coupling)


def effective_youngs_modulus(G: float, nu: float = 0.45) -> float:
    """Reduced modulus for Hertz contact [Pa].

    E* = E / (1 - nu^2) = 2·G·(1+nu) / (1 - nu^2)

    The Hertz contact equation F = (4/3)·E*·sqrt(R)·delta^(3/2) requires
    the reduced modulus, not just Young's modulus.
    """
    E = 2.0 * G * (1.0 + nu)
    return E / (1.0 - nu**2)


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

    R = max(R, 1e-9)  # prevent sqrt(0)
    delta = np.linspace(0, delta_max, n_points)
    F = (4.0 / 3.0) * E_star * np.sqrt(R) * delta ** 1.5
    return delta, F


def ogston_kav(rh: np.ndarray, r_fiber: float, phi_fiber: float) -> np.ndarray:
    """Ogston partition coefficient for size-exclusion chromatography.

    Kav = exp(-phi_f * ((rh/r_f) + 1)^2)

    The fiber volume fraction ``phi_fiber`` should be derived from the
    actual polymer concentration divided by the dry polymer density
    (phi = c_polymer / rho_polymer), **not** from ``1 - porosity``.
    The thresholded porosity treats the entire polymer-rich phase as
    solid fiber, but that phase is still mostly water.

    Parameters
    ----------
    rh : np.ndarray
        Hydrodynamic radius of solute [m].
    r_fiber : float
        Gel fiber radius [m].
    phi_fiber : float
        Fiber (polymer) volume fraction [-].  Should be computed as
        (c_agarose + c_chitosan) / rho_polymer with rho_polymer ~ 1400 kg/m3.
    """
    ratio = rh / r_fiber + 1.0
    return np.exp(-phi_fiber * ratio ** 2)


def solve_mechanical(params: SimulationParameters,
                     props: MaterialProperties,
                     gelation: GelationResult,
                     crosslinking: CrosslinkingResult,
                     R_droplet: float = None) -> MechanicalResult:
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

    # Crosslinked network modulus from L3
    G_xlink = crosslinking.G_chitosan_final

    # Route to appropriate modulus model based on network metadata from L3.
    # If metadata is available, the model accounts for chemistry-specific
    # coupling (amine vs hydroxyl vs ionic vs independent PEGDA).
    network_meta = getattr(crosslinking, 'network_metadata', None)
    G_DN = select_modulus_model(G_agar, G_xlink, network_meta, props.eta_coupling)

    # Effective Young's modulus
    E_star = effective_youngs_modulus(G_DN)

    # Hertz contact curve
    # Bead radius: use actual R_droplet if provided, else fall back to
    # L_domain (legacy behavior with deprecation warning).
    if R_droplet is not None and R_droplet > 0:
        R = R_droplet
    else:
        import logging
        logging.getLogger(__name__).warning(
            "solve_mechanical: R_droplet not provided, falling back to "
            "L_domain/2 (may be truncated for large droplets)"
        )
        if gelation.L_domain > 0:
            R = gelation.L_domain / 2.0
        else:
            R = gelation.r_grid[-1] + (gelation.r_grid[1] - gelation.r_grid[0]) / 2.0
    delta_arr, F_arr = hertz_contact(E_star, R)

    # Ogston Kav for a range of protein sizes
    rh_arr = np.logspace(np.log10(1e-9), np.log10(50e-9), 50)  # 1-50 nm
    # Fiber volume fraction from actual polymer concentration
    # (not from thresholded porosity, which includes water in polymer-rich phase)
    rho_polymer = 1400.0  # kg/m³ (dry polymer density for agarose/chitosan)
    phi_fiber = max((params.formulation.c_agarose + params.formulation.c_chitosan) / rho_polymer, 0.0)
    Kav_arr = ogston_kav(rh_arr, props.r_fiber, phi_fiber)

    return MechanicalResult(
        G_agarose=float(G_agar),
        G_chitosan=float(G_xlink),
        G_DN=float(G_DN),
        E_star=float(E_star),
        delta_array=delta_arr,
        F_array=F_arr,
        rh_array=rh_arr,
        Kav_array=Kav_arr,
        pore_size_mean=float(gelation.pore_size_mean),
        xi_mesh=float(crosslinking.xi_final),
    )
