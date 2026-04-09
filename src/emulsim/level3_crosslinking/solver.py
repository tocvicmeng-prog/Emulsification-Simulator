"""Level 3: Crosslinking kinetics solver with per-mechanism dispatch.

Models the evolution of crosslink density, mesh size, and modulus
during crosslinking of the polymer network within the already-gelled
microsphere.  Supports multiple crosslinker chemistries:

  - second_order (amine): genipin, glutaraldehyde, EDC/NHS
  - second_order (hydroxyl): epichlorohydrin, DVS, citric acid
  - uv_dose: PEGDA + UV photoinitiation
  - ionic_instant: TPP (sodium tripolyphosphate)

The Thiele modulus is computed from actual L1/L2 outputs (droplet radius,
porosity) to validate the ODE assumption.  When Phi << 1 the system is
reaction-limited and the simple ODE is valid; when Phi >> 1 a spatially
resolved reaction-diffusion PDE would be needed.
"""

from __future__ import annotations

import logging
import math
from typing import Optional

import numpy as np
from scipy.integrate import solve_ivp

from ..datatypes import (
    CrosslinkingResult,
    MaterialProperties,
    NetworkTypeMetadata,
    SimulationParameters,
)

logger = logging.getLogger(__name__)

# Constants
R_GAS = 8.314       # J/(mol*K)
K_BOLTZMANN = 1.38e-23  # J/K
N_AVOGADRO = 6.022e23

# Hydroxyl-reactive mechanisms (target agarose -OH network)
_HYDROXYL_MECHANISMS = frozenset({"hydroxyl", "michael_addition", "ester_bond"})


def compute_thiele_modulus(R: float, k: float, NH2: float,
                           D_genipin: float = 1e-10) -> float:
    """Thiele modulus for crosslinker diffusion into microsphere.

    Phi = R * sqrt(k * [reactive_groups] / D)

    Phi << 1: reaction-limited (uniform crosslinking, ODE valid)
    Phi >> 1: diffusion-limited (need reaction-diffusion PDE)

    Parameters
    ----------
    R : float
        Microsphere radius [m].
    k : float
        Second-order rate constant [m3/(mol*s)].
    NH2 : float
        Available reactive group concentration [mol/m3].
    D_genipin : float
        Effective diffusivity of crosslinker in gel [m2/s].
    """
    if D_genipin <= 0 or R <= 0:
        return 0.0
    return R * np.sqrt(k * NH2 / D_genipin)


def arrhenius_rate_constant(T: float, k0: float, E_a: float) -> float:
    """Arrhenius rate constant [m3/(mol*s)].

    k(T) = k0 * exp(-E_a / (R*T))
    """
    return k0 * np.exp(-E_a / (R_GAS * T))


# Keep old name as alias for backward compatibility
genipin_rate_constant = arrhenius_rate_constant


def available_amine_concentration(c_chitosan: float, DDA: float,
                                   M_GlcN: float) -> float:
    """Concentration of available primary amine groups [mol/m3].

    [NH2] = c_chitosan * DDA / M_GlcN
    where c_chitosan is in kg/m3, DDA is degree of deacetylation,
    and M_GlcN is molar mass of glucosamine unit in g/mol.
    """
    return c_chitosan * 1000.0 * DDA / M_GlcN  # kg/m3 -> g/m3, then /M -> mol/m3


def available_hydroxyl_concentration(c_agarose: float) -> float:
    """Concentration of available hydroxyl groups from agarose [mol/m3].

    Each agarose disaccharide repeat unit (~306 g/mol) has ~3 free OH groups.
    [OH] = c_agarose * 1000 * 3 / 306.0
    where c_agarose is in kg/m3.
    """
    M_agarose_repeat = 306.0  # g/mol for disaccharide repeat unit
    OH_per_repeat = 3.0
    return c_agarose * 1000.0 * OH_per_repeat / M_agarose_repeat


def crosslinking_odes(t: float, y: np.ndarray, k: float,
                      NH2_total: float) -> np.ndarray:
    """RHS of the crosslinking ODE system.

    State variables y = [X, G_free, NH2_free]:
        X: crosslink density [mol/m3]
        G_free: free crosslinker concentration [mol/m3]
        NH2_free: free reactive group concentration [mol/m3]

    dX/dt = k * G_free * NH2_free
    dG_free/dt = -k * G_free * NH2_free
    dNH2_free/dt = -2 * k * G_free * NH2_free

    Each crosslink bridge consumes 1 crosslinker and 2 reactive groups
    (one on each chain).
    """
    X, G_free, NH2_free = y
    # Ensure non-negative
    G_free = max(G_free, 0.0)
    NH2_free = max(NH2_free, 0.0)

    rate = k * G_free * NH2_free
    return np.array([rate, -rate, -2.0 * rate])


def crosslink_density_to_properties(
    X: float, c_polymer: float, T: float, DDA: float, M_GlcN: float,
    reactive_total: float, f_bridge: float = 0.4,
) -> tuple[float, float, float, float]:
    """Convert crosslink density to network properties.

    Parameters
    ----------
    c_polymer : float
        Concentration of the target polymer network [kg/m3].
    reactive_total : float
        Total initial reactive group concentration [mol/m3].
    f_bridge : float
        Fraction of crosslinker reactions that produce elastically
        active inter-chain crosslinks.

    Returns (nu_e, Mc, xi, G):
        nu_e: effective crosslink density [1/m3]
        Mc: molecular weight between crosslinks [g/mol]
        xi: mesh size [m]
        G: shear modulus of crosslinked network [Pa]
    """
    X_effective = f_bridge * X

    # Effective crosslink density (assuming tetrafunctional crosslinks)
    nu_e = X_effective * N_AVOGADRO

    # Molecular weight between crosslinks using effective conversion
    p_effective = X_effective / reactive_total if reactive_total > 0 else 0.0
    p_effective = min(p_effective, 1.0)

    # Mc = M_repeat / (2*p_eff)  for random crosslinking on a linear chain
    if p_effective > 1e-10:
        Mc = M_GlcN / (2.0 * p_effective)
    else:
        Mc = 1e10  # effectively infinite

    # Mesh size from Canal-Peppas
    rho_polymer = 1400.0  # kg/m3 for chitosan/agarose
    nu_2s = c_polymer / rho_polymer
    nu_2s = max(nu_2s, 0.001)

    if Mc > 0 and Mc < 1e8:
        xi_nm = 0.071 * nu_2s ** (-1.0 / 3.0) * np.sqrt(Mc)
        xi = xi_nm * 1e-9  # nm -> m
    else:
        xi = 1e-6  # 1 um default for uncrosslinked

    # Shear modulus from rubber elasticity: G = nu_e * kT
    G = nu_e * K_BOLTZMANN * T

    return nu_e, Mc, xi, G


def _thiele_check(R_droplet, porosity, k, reactive_total):
    """Run Thiele modulus diagnostic if droplet radius is available."""
    D_bare = 1e-10  # m2/s, small molecule in water
    if R_droplet is not None and R_droplet > 0:
        if porosity is not None and porosity > 0:
            D_eff = D_bare * porosity ** 2
        else:
            D_eff = D_bare

        thiele = compute_thiele_modulus(R_droplet, k, reactive_total, D_eff)
        logger.info("L3 Thiele modulus: %.4f  (R=%.2f um, D_eff=%.2e m2/s)",
                     thiele, R_droplet * 1e6, D_eff)

        if thiele > 1.0:
            logger.warning(
                "Thiele modulus %.2f >> 1: diffusion-limited regime. "
                "A reaction-diffusion PDE model is recommended.", thiele)
        elif thiele > 0.3:
            logger.warning(
                "Thiele modulus %.2f in borderline range (0.3-1.0). "
                "Uniform-crosslinking assumption may lose accuracy.", thiele)
    else:
        logger.debug("L3: R_droplet not supplied; skipping Thiele check.")


def _build_result_arrays(t_arr, X_arr, c_polymer, T, DDA, M_GlcN,
                         reactive_total, f_bridge):
    """Compute derived property arrays from crosslink density history."""
    nu_e_arr = np.zeros_like(t_arr)
    Mc_arr = np.zeros_like(t_arr)
    xi_arr = np.zeros_like(t_arr)
    G_arr = np.zeros_like(t_arr)

    for i in range(len(t_arr)):
        nu_e, Mc, xi, G = crosslink_density_to_properties(
            X_arr[i], c_polymer, T, DDA, M_GlcN, reactive_total, f_bridge,
        )
        nu_e_arr[i] = nu_e
        Mc_arr[i] = Mc
        xi_arr[i] = xi
        G_arr[i] = G

    # Fix 2: Each bridge consumes 2 reactive groups, so p = 2*X / reactive_total
    p_final = 2.0 * X_arr[-1] / reactive_total if reactive_total > 0 else 0.0

    return CrosslinkingResult(
        t_array=t_arr,
        X_array=X_arr,
        nu_e_array=nu_e_arr,
        Mc_array=Mc_arr,
        xi_array=xi_arr,
        G_chitosan_array=G_arr,
        p_final=float(min(p_final, 1.0)),
        nu_e_final=float(nu_e_arr[-1]),
        Mc_final=float(Mc_arr[-1]),
        xi_final=float(xi_arr[-1]),
        G_chitosan_final=float(G_arr[-1]),
    )


def _build_zero_result(params):
    """Return a CrosslinkingResult with all zeros (no crosslinking)."""
    t_xlink = params.formulation.t_crosslink
    t_arr = np.array([0.0, max(t_xlink, 1.0)])
    zeros = np.zeros_like(t_arr)
    xi_default = np.full_like(t_arr, 1e-6)  # 1 um uncrosslinked default
    Mc_default = np.full_like(t_arr, 1e10)
    return CrosslinkingResult(
        t_array=t_arr,
        X_array=zeros,
        nu_e_array=zeros,
        Mc_array=Mc_default,
        xi_array=xi_default,
        G_chitosan_array=zeros,
        p_final=0.0,
        nu_e_final=0.0,
        Mc_final=1e10,
        xi_final=1e-6,
        G_chitosan_final=0.0,
    )


# ═══════════════════════════════════════════════════════════════════════════
#  Per-mechanism solver implementations
# ═══════════════════════════════════════════════════════════════════════════

def _solve_second_order_amine(params, props, xl, R_droplet, porosity):
    """Second-order kinetics targeting chitosan -NH2 groups.

    Used by: genipin, glutaraldehyde, EDC/NHS.
    """
    T = params.formulation.T_crosslink
    t_xlink = params.formulation.t_crosslink
    c_crosslinker = max(params.formulation.c_genipin, 0.0)
    c_chitosan = params.formulation.c_chitosan

    if T <= 0:
        raise ValueError(f"Crosslinking temperature must be positive, got {T} K")
    if props.M_GlcN <= 0:
        raise ValueError("M_GlcN must be positive")

    # Fix 1: Prefer crosslinker profile values over props when xl is provided
    k0 = xl.k_xlink_0 if xl else props.k_xlink_0
    Ea = xl.E_a_xlink if xl else props.E_a_xlink
    fb = xl.f_bridge if xl else props.f_bridge

    # Fix 6: Clamp to non-negative
    k0 = max(k0, 0.0)
    fb = max(fb, 0.0)

    k = arrhenius_rate_constant(T, k0, Ea)
    NH2_total = max(available_amine_concentration(c_chitosan, props.DDA, props.M_GlcN), 0.0)

    _thiele_check(R_droplet, porosity, k, NH2_total)

    y0 = np.array([0.0, c_crosslinker, NH2_total])
    t_span = (0.0, t_xlink)
    t_eval = np.linspace(0, t_xlink, 200)

    sol = solve_ivp(
        crosslinking_odes, t_span, y0,
        args=(k, NH2_total),
        method=params.solver.l3_method,
        rtol=params.solver.l3_rtol,
        atol=params.solver.l3_atol,
        t_eval=t_eval,
    )
    if not sol.success:
        logger.warning("Crosslinking ODE solver did not converge: %s", sol.message)

    return _build_result_arrays(
        sol.t, sol.y[0], c_chitosan, T, props.DDA, props.M_GlcN,
        NH2_total, fb,
    )


def _solve_second_order_hydroxyl(params, props, xl, R_droplet, porosity):
    """Second-order kinetics targeting agarose -OH groups.

    Used by: epichlorohydrin (ECH), divinyl sulfone (DVS), citric acid.
    Same ODE structure as amine but reactive groups are agarose hydroxyl
    groups and the modulus contribution goes to the agarose network.

    NOTE: This model is simplified -- it treats all hydroxyl crosslinkers
    as bifunctional targeting agarose-OH only.  The actual chemistry of
    ECH (epoxide ring-opening) and citric acid (ester/amide formation)
    is more complex and may involve chitosan -NH2 as well.
    """
    T = params.formulation.T_crosslink
    t_xlink = params.formulation.t_crosslink
    c_crosslinker = max(params.formulation.c_genipin, 0.0)
    c_agarose = params.formulation.c_agarose

    if T <= 0:
        raise ValueError(f"Crosslinking temperature must be positive, got {T} K")

    # Fix 1: Prefer crosslinker profile values over props when xl is provided
    k0 = xl.k_xlink_0 if xl else props.k_xlink_0
    Ea = xl.E_a_xlink if xl else props.E_a_xlink
    fb = xl.f_bridge if xl else props.f_bridge

    # Fix 6: Clamp to non-negative
    k0 = max(k0, 0.0)
    fb = max(fb, 0.0)

    # Fix 3: Log simplification note for hydroxyl-reactive crosslinkers
    logger.info("L3 hydroxyl model: using simplified bifunctional agarose-OH model "
                "for '%s'. Actual chemistry (epoxide ring-opening / ester formation) "
                "is more complex and may also target chitosan -NH2.",
                xl.name if xl else "unknown")

    k = arrhenius_rate_constant(T, k0, Ea)
    OH_total = max(available_hydroxyl_concentration(c_agarose), 0.0)

    _thiele_check(R_droplet, porosity, k, OH_total)

    y0 = np.array([0.0, c_crosslinker, OH_total])
    t_span = (0.0, t_xlink)
    t_eval = np.linspace(0, t_xlink, 200)

    sol = solve_ivp(
        crosslinking_odes, t_span, y0,
        args=(k, OH_total),
        method=params.solver.l3_method,
        rtol=params.solver.l3_rtol,
        atol=params.solver.l3_atol,
        t_eval=t_eval,
    )
    if not sol.success:
        logger.warning("Crosslinking ODE solver did not converge: %s", sol.message)

    # Use agarose repeat unit MW for Mc calculation
    M_repeat = 306.0  # g/mol, agarose disaccharide
    return _build_result_arrays(
        sol.t, sol.y[0], c_agarose, T, props.DDA, M_repeat,
        OH_total, fb,
    )


def _solve_uv_dose(params, props, xl, uv_intensity, R_droplet, porosity):
    """UV-initiated radical polymerization kinetics.

    Used by: PEGDA + UV photoinitiation.

    Pseudo-first-order in PEGDA concentration:
        k_eff = k_uv * sqrt(I_uv)
        [PEGDA](t) = [PEGDA]_0 * exp(-k_eff * t)
        X(t) = [PEGDA]_0 * (1 - exp(-k_eff * t))

    PEGDA forms its own independent network (not chitosan, not agarose).
    """
    T = params.formulation.T_crosslink
    t_xlink = params.formulation.t_crosslink
    c_pegda = max(params.formulation.c_genipin, 0.0)  # reuse crosslinker conc field

    if T <= 0:
        raise ValueError(f"Crosslinking temperature must be positive, got {T} K")

    # Fix 1: Prefer crosslinker profile values over props when xl is provided
    k0 = xl.k_xlink_0 if xl else props.k_xlink_0
    Ea = xl.E_a_xlink if xl else props.E_a_xlink
    fb = xl.f_bridge if xl else props.f_bridge

    # Fix 6: Clamp to non-negative
    k0 = max(k0, 0.0)
    fb = max(fb, 0.0)

    # Fix 4: UV=0 should give zero conversion, not coerce to 0.001
    if uv_intensity <= 0:
        logger.info("L3 UV dose model: I_uv=0 -- no UV, returning zero conversion.")
        t_arr = np.linspace(0.0, t_xlink, 200)
        X_arr = np.zeros_like(t_arr)

        M_pegda_repeat = 700.0
        c_pegda_kg = c_pegda * M_pegda_repeat / 1000.0
        return _build_result_arrays(
            t_arr, X_arr, max(c_pegda_kg, 0.001), T, 1.0, M_pegda_repeat,
            c_pegda, fb,
        )

    # Fix 4: Check t_crosslink > 0
    if t_xlink <= 0:
        logger.warning("L3 UV dose model: t_crosslink <= 0, returning zero conversion.")
        return _build_zero_result(params)

    # k_eff = k_uv * sqrt(I_uv) where I_uv is in mW/cm2
    k_uv = arrhenius_rate_constant(T, k0, Ea)
    k_eff = k_uv * np.sqrt(uv_intensity)

    # Fix 4: Check k_eff is finite
    if not np.isfinite(k_eff):
        logger.warning("L3 UV dose model: k_eff is not finite (%.3e), clamping to 0.", k_eff)
        k_eff = 0.0

    logger.info("L3 UV dose model: k_uv=%.3e, I_uv=%.1f mW/cm2, k_eff=%.3e",
                k_uv, uv_intensity, k_eff)

    # Time array
    t_arr = np.linspace(0.0, t_xlink, 200)

    # First-order consumption: X(t) = c_pegda * (1 - exp(-k_eff * t))
    X_arr = c_pegda * (1.0 - np.exp(-k_eff * t_arr))

    # PEGDA network properties
    # Use PEGDA as its own polymer network; approximate repeat unit MW
    M_pegda_repeat = 700.0  # g/mol for PEGDA-700 (common)
    # Effective polymer concentration for mesh-size calculation
    # PEGDA is in the aqueous phase; use its mass concentration
    c_pegda_kg = c_pegda * M_pegda_repeat / 1000.0  # mol/m3 * g/mol / 1000 = kg/m3

    return _build_result_arrays(
        t_arr, X_arr, max(c_pegda_kg, 0.001), T, 1.0, M_pegda_repeat,
        c_pegda, fb,
    )


def _solve_ionic_instant(params, props, xl, R_droplet, porosity):
    """Instantaneous ionic gelation (equilibrium model).

    Used by: TPP (sodium tripolyphosphate).

    No ODE needed -- crosslinking is at equilibrium immediately.
    Degree of crosslinking depends on TPP:NH2 stoichiometric ratio.
    Ionic crosslinks are reversible and weaker than covalent bonds.
    """
    T = params.formulation.T_crosslink
    t_xlink = params.formulation.t_crosslink
    c_tpp = max(params.formulation.c_genipin, 0.0)  # reuse crosslinker conc field
    c_chitosan = params.formulation.c_chitosan

    if T <= 0:
        raise ValueError(f"Crosslinking temperature must be positive, got {T} K")
    if props.M_GlcN <= 0:
        raise ValueError("M_GlcN must be positive")

    NH2_total = max(available_amine_concentration(c_chitosan, props.DDA, props.M_GlcN), 0.0)

    # Fix 5: When NH2_total <= 0, return zero crosslinking result directly
    if NH2_total <= 1e-10:
        logger.info("L3 ionic model: NH2_total ~ 0, no reactive groups for ionic crosslinking.")
        return _build_zero_result(params)

    # Fix 1: Prefer crosslinker profile values over props when xl is provided
    fb = xl.f_bridge if xl else props.f_bridge

    # Fix 6: Clamp to non-negative
    fb = max(fb, 0.0)

    # TPP has 5 negative charges (P3O10^5-), each can bind one NH3+
    z_tpp = 5.0
    equivalents = c_tpp * z_tpp / max(NH2_total, 1e-10)
    p_ionic = min(equivalents, 1.0)

    # Ionic crosslinks are weaker than covalent -- stability factor
    f_ionic_strength = 0.3

    # Crosslink density at equilibrium
    X_final = p_ionic * NH2_total / 2.0

    logger.info("L3 ionic model: c_TPP=%.2f mM, NH2=%.1f mol/m3, "
                "equivalents=%.3f, p_ionic=%.3f",
                c_tpp, NH2_total, equivalents, p_ionic)

    # Build simple two-point time array (instantaneous equilibrium)
    t_arr = np.array([0.0, t_xlink])
    X_arr = np.array([X_final, X_final])  # already at equilibrium at t=0

    # Compute properties with ionic stability factor applied to f_bridge
    f_bridge_ionic = fb * f_ionic_strength

    nu_e_arr = np.zeros_like(t_arr)
    Mc_arr = np.zeros_like(t_arr)
    xi_arr = np.zeros_like(t_arr)
    G_arr = np.zeros_like(t_arr)

    for i in range(len(t_arr)):
        nu_e, Mc, xi, G = crosslink_density_to_properties(
            X_arr[i], c_chitosan, T, props.DDA, props.M_GlcN,
            NH2_total, f_bridge_ionic,
        )
        nu_e_arr[i] = nu_e
        Mc_arr[i] = Mc
        xi_arr[i] = xi
        G_arr[i] = G

    # Fix 2: Each bridge consumes 2 reactive groups, so p = 2*X / reactive_total
    p_final = 2.0 * X_final / NH2_total if NH2_total > 0 else 0.0

    return CrosslinkingResult(
        t_array=t_arr,
        X_array=X_arr,
        nu_e_array=nu_e_arr,
        Mc_array=Mc_arr,
        xi_array=xi_arr,
        G_chitosan_array=G_arr,
        p_final=float(min(p_final, 1.0)),
        nu_e_final=float(nu_e_arr[-1]),
        Mc_final=float(Mc_arr[-1]),
        xi_final=float(xi_arr[-1]),
        G_chitosan_final=float(G_arr[-1]),
    )


# ═══════════════════════════════════════════════════════════════════════════
#  Main dispatcher
# ═══════════════════════════════════════════════════════════════════════════

def solve_crosslinking(params: SimulationParameters,
                       props: MaterialProperties,
                       R_droplet: Optional[float] = None,
                       porosity: Optional[float] = None,
                       crosslinker_key: str = 'genipin',
                       uv_intensity: float = 0.0) -> CrosslinkingResult:
    """Dispatch to the appropriate kinetics model based on crosslinker.

    Parameters
    ----------
    params : SimulationParameters
        Contains crosslinking temperature, time, crosslinker concentration.
    props : MaterialProperties
        Contains kinetic constants, polymer properties.
    R_droplet : float, optional
        Microsphere radius [m] from L1 (d50/2).
    porosity : float, optional
        Gel porosity from L2.
    crosslinker_key : str
        Key into CROSSLINKERS library (default 'genipin' for backward
        compatibility).
    uv_intensity : float
        UV intensity [mW/cm2] for UV-initiated crosslinkers.
    """
    from ..reagent_library import CROSSLINKERS

    # Fix 6: Validate T_crosslink at the dispatcher level
    T = params.formulation.T_crosslink
    if not math.isfinite(T) or T <= 0:
        raise ValueError(f"T_crosslink must be finite and positive, got {T}")

    xl = CROSSLINKERS.get(crosslinker_key)

    # Fix 1: Unknown crosslinker key raises ValueError instead of silent fallthrough
    if xl is None:
        raise ValueError(f"Unknown crosslinker_key: {crosslinker_key!r}. "
                         f"Valid keys: {list(CROSSLINKERS.keys())}")

    # ── Dispatch to mechanism-specific solver ──────────────────────────
    if xl.kinetics_model == 'second_order':
        if xl.mechanism in _HYDROXYL_MECHANISMS:
            result = _solve_second_order_hydroxyl(params, props, xl, R_droplet, porosity)
            metadata = NetworkTypeMetadata(
                solver_family="hydroxyl_covalent",
                network_target="mixed",
                bond_type="covalent",
                is_true_second_network=True,
            )
        else:
            result = _solve_second_order_amine(params, props, xl, R_droplet, porosity)
            metadata = NetworkTypeMetadata(
                solver_family="amine_covalent",
                network_target="chitosan",
                bond_type="covalent",
                is_true_second_network=True,
            )
    elif xl.kinetics_model == 'uv_dose':
        result = _solve_uv_dose(params, props, xl, uv_intensity, R_droplet, porosity)
        metadata = NetworkTypeMetadata(
            solver_family="independent_network",
            network_target="independent",
            bond_type="covalent",
            is_true_second_network=True,
        )
    elif xl.kinetics_model == 'ionic_instant':
        result = _solve_ionic_instant(params, props, xl, R_droplet, porosity)
        metadata = NetworkTypeMetadata(
            solver_family="ionic_reversible",
            network_target="chitosan",
            bond_type="ionic",
            is_true_second_network=False,
        )
    elif xl.kinetics_model == 'michaelis_menten':
        # Michaelis-Menten crosslinkers (e.g. EDC/NHS) fall back to
        # second-order amine kinetics as an approximation until a
        # dedicated solver is implemented.
        logger.info("L3: michaelis_menten model not yet implemented; "
                    "falling back to second_order amine kinetics for '%s'. "
                    "NOTE: results are approximate -- EDC/NHS kinetics involve "
                    "competitive hydrolysis not captured by simple second-order model.",
                    crosslinker_key)
        result = _solve_second_order_amine(params, props, xl, R_droplet, porosity)
        metadata = NetworkTypeMetadata(
            solver_family="amine_covalent",
            network_target="chitosan",
            bond_type="covalent",
            is_true_second_network=True,
        )
    else:
        raise ValueError(f"Unknown kinetics model: {xl.kinetics_model}")

    # Attach per-chemistry network metadata to result
    result.network_metadata = metadata
    return result
