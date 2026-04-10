"""Level 4: Mechanical property prediction for double-network microspheres.

Analytical models:
- Additive DN modulus: G_DN = G_agarose + G_chitosan
- Hertz contact for microsphere compression
- Ogston model for chromatographic Kav
- Flory-Rehner affine IPN model (mechanistic_research mode)
"""

from __future__ import annotations

import logging

import numpy as np

from ..datatypes import (
    CrosslinkingResult,
    GelationResult,
    MaterialProperties,
    MechanicalResult,
    ModelMode,
    SimulationParameters,
)

logger = logging.getLogger(__name__)


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


def hashin_shtrikman_bounds(G1: float, G2: float, phi1: float) -> tuple[float, float]:
    """Hashin-Shtrikman bounds for two-phase composite shear modulus [Pa].

    Parameters
    ----------
    G1, G2 : float
        Shear moduli of phases 1 and 2 [Pa].
    phi1 : float
        Volume fraction of phase 1 [-].

    Returns
    -------
    G_lower, G_upper : float
        Lower and upper HS bounds [Pa].

    Notes
    -----
    Uses the Hashin-Shtrikman (1963) formulas for shear modulus bounds in
    incompressible two-phase composites:

        G_L = G_soft  + phi_stiff / (1/(G_stiff - G_soft) + 6*phi_soft/(5*G_soft*(4-5*nu)))
        G_U = G_stiff + phi_soft  / (1/(G_soft  - G_stiff) + 6*phi_stiff/(5*G_stiff*(4-5*nu)))

    For nearly incompressible materials (nu -> 0.5), (4 - 5*nu) -> 1.5, so the
    denominator factor simplifies to 6/(5*G*1.5) = 4/(5*G).  We use nu=0.5
    exactly (hydrogel limit), giving:

        G_L = G_soft  + phi_stiff / (1/(G_stiff - G_soft) + 2*phi_soft/(5*G_soft))

    which matches the standard incompressible HS shear modulus formula.
    """
    phi1 = float(np.clip(phi1, 0.0, 1.0))
    phi2 = 1.0 - phi1
    G1 = max(float(G1), 0.0)
    G2 = max(float(G2), 0.0)

    # Degenerate cases
    if G1 == G2:
        return G1, G1
    if phi1 == 0.0:
        return G2, G2
    if phi1 == 1.0:
        return G1, G1

    # Identify soft / stiff phase
    G_soft, G_stiff = (G1, G2) if G1 <= G2 else (G2, G1)
    phi_soft = phi2 if G1 <= G2 else phi1
    phi_stiff = 1.0 - phi_soft

    # Incompressible (nu=0.5) HS bounds
    # Lower bound: soft phase is matrix
    if G_soft > 0.0:
        denom_L = 1.0 / (G_stiff - G_soft) + 2.0 * phi_soft / (5.0 * G_soft)
        G_lower = G_soft + phi_stiff / denom_L
    else:
        # If G_soft = 0, lower bound is 0
        G_lower = 0.0

    # Upper bound: stiff phase is matrix
    # Guard: denominator can approach zero for large modulus contrasts
    denom_U = 1.0 / (G_soft - G_stiff) + 2.0 * phi_stiff / (5.0 * G_stiff)
    if abs(denom_U) < 1e-30:
        G_upper = G_stiff
    else:
        G_upper = G_stiff + phi_soft / denom_U
    # HS upper bound cannot exceed the modulus of the stiffest phase
    G_upper = max(min(G_upper, G_stiff), G_lower)

    return float(G_lower), float(G_upper)


# ── Flory-Rehner Affine IPN Model (Protocol §M3-P3) ──────────────────────


def flory_rehner_swelling(
    nu_e: float,
    phi_0: float,
    chi: float,
    T: float = 298.15,
) -> tuple[float, bool]:
    """Equilibrium polymer volume fraction from Flory-Rehner theory.

    Solves the Flory-Rehner equation for equilibrium swelling of a
    single crosslinked network:

        ln(1-phi) + phi + chi*phi^2 + nu_e*V_s*(phi^(1/3) - phi/2) = 0

    Parameters
    ----------
    nu_e : float
        Effective crosslink density [mol/m^3].
    phi_0 : float
        Reference (dry) polymer volume fraction [-].
    chi : float
        Flory-Huggins interaction parameter [-].
    T : float
        Temperature [K].

    Returns
    -------
    phi_eq : float
        Equilibrium swollen polymer volume fraction [-].
    converged : bool
        True if brentq found a root; False if fell back to phi_0.

    Notes
    -----
    Uses scipy.optimize.brentq for robust root-finding within the bracket
    [phi_0*0.01, 0.999].  If no root is found (e.g. very low crosslink density),
    falls back to phi_0 with a warning.
    """
    from scipy.optimize import brentq

    V_s = 18.015e-6  # [m^3/mol] molar volume of water

    def objective(phi: float) -> float:
        """Chemical potential balance: mu_mixing + mu_elastic = 0."""
        # Guard: phi must be strictly in (0, 1)
        if phi <= 0.0 or phi >= 1.0:
            return 1e10
        # Flory-Huggins mixing term + affine elastic term
        return (np.log(1.0 - phi) + phi + chi * phi**2
                + nu_e * V_s * (phi ** (1.0 / 3.0) - phi / 2.0))

    # Bracket: phi must be between a small fraction of phi_0 (very swollen)
    # and nearly pure polymer (0.999)
    phi_lo = max(phi_0 * 0.01, 1e-6)
    phi_hi = 0.999

    try:
        phi_eq = brentq(objective, phi_lo, phi_hi, xtol=1e-8, maxiter=200)
        return float(phi_eq), True
    except ValueError:
        # Bracket does not contain root; return reference volume fraction
        logger.warning(
            "Flory-Rehner: no root in [%.4f, %.4f] for nu_e=%.1f, chi=%.3f. "
            "Returning phi_0=%.4f.", phi_lo, phi_hi, nu_e, chi, phi_0,
        )
        return float(phi_0), False


def double_network_modulus_affine(
    nu_e1: float,
    nu_e2: float,
    phi_01: float,
    phi_02: float,
    chi1: float,
    chi2: float,
    T: float = 298.15,
) -> tuple[float, str]:
    """Affine IPN modulus from Flory-Rehner swelling-constrained model [Pa].

    Each network swells independently to equilibrium, then the IPN
    constrains both networks. The effective modulus is:

        G_DN = nu_e1 * N_A * kT * lambda_1^(-2/3) + nu_e2 * N_A * kT * lambda_2^(-2/3)

    where lambda_i = (phi_0i / phi_eq_i)^(1/3) is the swelling stretch ratio
    and phi_eq_i is the equilibrium polymer volume fraction from Flory-Rehner theory.

    Parameters
    ----------
    nu_e1 : float
        Effective crosslink density of network 1 (agarose) [mol/m^3].
    nu_e2 : float
        Effective crosslink density of network 2 (chitosan) [mol/m^3].
    phi_01 : float
        Reference polymer volume fraction of network 1 [-].
    phi_02 : float
        Reference polymer volume fraction of network 2 [-].
    chi1 : float
        Flory-Huggins parameter for network 1 [-].
        ASSUMPTION: chi1 ~ 0.50 for agarose-water (slightly poor solvent).
    chi2 : float
        Flory-Huggins parameter for network 2 [-].
        ASSUMPTION: chi2 ~ 0.45 for chitosan-water.
    T : float
        Temperature [K].

    Returns
    -------
    G_DN : float
        Double-network shear modulus [Pa].
    status : str
        "converged" if both Flory-Rehner solutions converged, else "fallback".

    Notes
    -----
    Physical constants:
    - k_B = 1.380649e-23 J/K  (exact, SI 2019)
    - N_A = 6.02214076e23 /mol (exact, SI 2019)
    - G_i = nu_e_i [mol/m^3] * N_A [/mol] * k_B [J/K] * T [K] * lambda_i^(-2/3)
          = nu_e_i * R * T * lambda_i^(-2/3)   where R = N_A * k_B = 8.314 J/(mol·K)
    """
    K_BOLTZMANN = 1.380649e-23  # [J/K]
    N_AVOGADRO = 6.02214076e23  # [1/mol]
    kT = K_BOLTZMANN * T  # [J]

    status = "converged"

    # --- Network 1: agarose ---
    if nu_e1 > 0 and phi_01 > 0:
        phi_eq1, fr_ok1 = flory_rehner_swelling(nu_e1, phi_01, chi1, T)
        if not fr_ok1:
            status = "fallback"
        if phi_eq1 > 0:
            lambda_1 = (phi_01 / phi_eq1) ** (1.0 / 3.0)
        else:
            lambda_1 = 1.0
            status = "fallback"
        # G1 = nu_e1 * N_A * kT * lambda_1^(-2/3)  [Pa]
        G1 = nu_e1 * N_AVOGADRO * kT * lambda_1 ** (-2.0 / 3.0)
    else:
        G1 = 0.0

    # --- Network 2: chitosan ---
    if nu_e2 > 0 and phi_02 > 0:
        phi_eq2, fr_ok2 = flory_rehner_swelling(nu_e2, phi_02, chi2, T)
        if not fr_ok2:
            status = "fallback"
        if phi_eq2 > 0:
            lambda_2 = (phi_02 / phi_eq2) ** (1.0 / 3.0)
        else:
            lambda_2 = 1.0
            status = "fallback"
        # G2 = nu_e2 * N_A * kT * lambda_2^(-2/3)  [Pa]
        G2 = nu_e2 * N_AVOGADRO * kT * lambda_2 ** (-2.0 / 3.0)
    else:
        G2 = 0.0

    G_DN = G1 + G2

    # Sanity: modulus cannot be negative
    if G_DN < 0:
        G_DN = 0.0
        status = "fallback"

    return float(G_DN), status


def select_modulus_model(G_agarose: float, G_xlink: float,
                          network_metadata=None,
                          eta_coupling: float = -0.15,
                          model_mode: str = "hybrid_coupled",
                          params: SimulationParameters | None = None,
                          crosslinking: CrosslinkingResult | None = None,
                          eta_is_custom: bool = False,
                          ) -> tuple[float, str]:
    """Route to the appropriate modulus model based on network metadata and mode.

    Returns (G_DN, model_label) where model_label identifies which model was used.

    In MECHANISTIC_RESEARCH mode with amine_covalent or hydroxyl_covalent
    chemistry, attempts the Flory-Rehner affine IPN model before falling
    back to the phenomenological formula.

    Falls back to phenomenological DN modulus if no metadata provided.
    Per-chemistry eta is read from network_metadata.eta_coupling_recommended
    when available; otherwise falls back to the eta_coupling parameter.
    If eta_is_custom is True, the eta_coupling parameter is used as-is
    (user explicitly set it via calibration).
    """
    if network_metadata is None:
        return double_network_modulus(G_agarose, G_xlink, eta_coupling), "phenomenological"

    # Extract per-chemistry eta from metadata, with fallback to parameter
    family = getattr(network_metadata, 'solver_family', 'amine_covalent')
    if eta_is_custom:
        eta_per_chem = eta_coupling  # User explicitly set via calibration
    else:
        eta_per_chem = getattr(network_metadata, 'eta_coupling_recommended', eta_coupling)

    if family == "ionic_reversible":
        return ionic_gel_modulus(G_agarose, G_xlink), "ionic_gel"
    elif family == "independent_network":
        return triple_network_modulus(G_agarose, 0.0, G_xlink, eta_13=eta_per_chem), "triple_network"

    # For amine_covalent and hydroxyl_covalent:
    # In MECHANISTIC_RESEARCH mode, attempt affine IPN model (Protocol §M3-P3-N3)
    if model_mode == "mechanistic_research" and crosslinking is not None and params is not None:
        try:
            N_AVOGADRO = 6.02214076e23  # [1/mol]

            # Network 2 (chitosan): crosslink density from L3 [1/m^3] -> [mol/m^3]
            nu_e2 = crosslinking.nu_e_final
            nu_e2_mol = nu_e2 / N_AVOGADRO

            # Network 1 (agarose): estimate nu_e from rubber elasticity G = nu_e * kT
            # ASSUMPTION: G_agarose arises from physical helix-bundle junctions, not
            # rubber-elastic crosslinks. Inverting via rubber elasticity gives a
            # fictitious crosslink density that is used only as a proxy for network
            # connectivity in the FR swelling calculation. This is a standard
            # approximation in IPN modelling (Drury & Mooney, 2003).
            T_xlink = params.formulation.T_crosslink
            kT = 1.380649e-23 * T_xlink
            nu_e1 = G_agarose / kT if kT > 0 else 0.0  # [1/m^3]
            nu_e1_mol = nu_e1 / N_AVOGADRO  # [mol/m^3]

            # Polymer volume fractions from concentration and dry density
            # ASSUMPTION: rho_agarose_dry ~ 1640 kg/m^3, rho_chitosan_dry ~ 1400 kg/m^3
            rho_agarose = 1640.0   # [kg/m^3]
            rho_chitosan = 1400.0  # [kg/m^3]
            phi_01 = params.formulation.c_agarose / rho_agarose
            phi_02 = params.formulation.c_chitosan / rho_chitosan

            # ASSUMPTION: chi1 ~ 0.50 for agarose-water, chi2 ~ 0.45 for chitosan-water
            chi1 = 0.50
            chi2 = 0.45

            G_affine, affine_status = double_network_modulus_affine(
                nu_e1_mol, nu_e2_mol, phi_01, phi_02, chi1, chi2, T_xlink,
            )

            if affine_status == "converged" and G_affine > 0:
                return G_affine, "flory_rehner_affine"
            else:
                logger.warning(
                    "Affine IPN model fallback: status=%s, G=%.1f Pa. "
                    "Using phenomenological.", affine_status, G_affine,
                )
                # Fall through to phenomenological with fallback label
                if family == "hydroxyl_covalent":
                    return double_network_modulus(G_agarose, G_xlink, eta_coupling=eta_per_chem), "phenomenological_fallback"
                else:
                    return double_network_modulus(G_agarose, G_xlink, eta_per_chem), "phenomenological_fallback"
        except Exception as exc:
            logger.warning(
                "Affine IPN model failed: %s. Using phenomenological.", exc,
            )

    # Phenomenological fallback (and default for hybrid_coupled / empirical_engineering)
    if family == "hydroxyl_covalent":
        return double_network_modulus(G_agarose, G_xlink, eta_coupling=eta_per_chem), "phenomenological"
    else:
        return double_network_modulus(G_agarose, G_xlink, eta_per_chem), "phenomenological"


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
                     R_droplet: float = None,
                     eta_is_custom: bool = False) -> MechanicalResult:
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
    # In MECHANISTIC_RESEARCH mode, attempts Flory-Rehner affine IPN model
    # before falling back to phenomenological (Protocol §M3-P3-N3).
    network_meta = getattr(crosslinking, 'network_metadata', None)
    model_mode_str = getattr(params.model_mode, 'value', 'hybrid_coupled')
    G_DN, model_label = select_modulus_model(
        G_agar, G_xlink, network_meta, props.eta_coupling,
        model_mode=model_mode_str,
        params=params,
        crosslinking=crosslinking,
        eta_is_custom=eta_is_custom,
    )

    # Hashin-Shtrikman bounds: use agarose volume fraction in the polymer phase
    c_agarose = params.formulation.c_agarose
    c_chitosan = params.formulation.c_chitosan
    _total_polymer = c_agarose + c_chitosan
    if _total_polymer > 0.0:
        phi1_hs = c_agarose / _total_polymer
    else:
        phi1_hs = 0.5
    G_DN_lower, G_DN_upper = hashin_shtrikman_bounds(G_agar, G_xlink, phi1_hs)

    # Effective Young's modulus
    E_star = effective_youngs_modulus(G_DN)

    # Hertz contact curve
    # Bead radius: use actual R_droplet if provided, else fall back to
    # L_domain (legacy behavior with deprecation warning).
    if R_droplet is not None and R_droplet > 0:
        R = R_droplet
    else:
        logger.warning(
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
        model_used=model_label,
        G_DN_lower=float(G_DN_lower),
        G_DN_upper=float(G_DN_upper),
    )
