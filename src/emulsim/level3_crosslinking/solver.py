"""Level 3: Genipin-chitosan crosslinking kinetics ODE solver.

Models the evolution of crosslink density, mesh size, and modulus
during genipin crosslinking of the chitosan network within the
already-gelled microsphere.

The Thiele modulus is computed from actual L1/L2 outputs (droplet radius,
porosity) to validate the ODE assumption.  When Φ << 1 the system is
reaction-limited and the simple ODE is valid; when Φ >> 1 a spatially
resolved reaction-diffusion PDE would be needed.
"""

from __future__ import annotations

import logging
from typing import Optional

import numpy as np
from scipy.integrate import solve_ivp

from ..datatypes import CrosslinkingResult, MaterialProperties, SimulationParameters

logger = logging.getLogger(__name__)

# Constants
R_GAS = 8.314       # J/(mol·K)
K_BOLTZMANN = 1.38e-23  # J/K
N_AVOGADRO = 6.022e23


def compute_thiele_modulus(R: float, k: float, NH2: float,
                           D_genipin: float = 1e-10) -> float:
    """Thiele modulus for genipin diffusion into microsphere.

    Phi = R * sqrt(k * [NH2] / D_genipin)

    Phi << 1: reaction-limited (uniform crosslinking, ODE valid)
    Phi >> 1: diffusion-limited (need reaction-diffusion PDE)

    Parameters
    ----------
    R : float
        Microsphere radius [m].
    k : float
        Second-order rate constant [m3/(mol*s)].
    NH2 : float
        Available amine concentration [mol/m3].
    D_genipin : float
        Effective diffusivity of genipin in gel [m2/s].
    """
    if D_genipin <= 0 or R <= 0:
        return 0.0
    return R * np.sqrt(k * NH2 / D_genipin)


def genipin_rate_constant(T: float, k0: float, E_a: float) -> float:
    """Arrhenius rate constant for genipin-chitosan reaction [m³/(mol·s)].

    k(T) = k0 · exp(-E_a / (R·T))
    """
    return k0 * np.exp(-E_a / (R_GAS * T))


def available_amine_concentration(c_chitosan: float, DDA: float,
                                   M_GlcN: float) -> float:
    """Concentration of available primary amine groups [mol/m³].

    [NH2] = c_chitosan · DDA / M_GlcN
    where c_chitosan is in kg/m³, DDA is degree of deacetylation,
    and M_GlcN is molar mass of glucosamine unit in g/mol.
    """
    return c_chitosan * 1000.0 * DDA / M_GlcN  # kg/m³ → g/m³, then /M → mol/m³


def crosslinking_odes(t: float, y: np.ndarray, k: float,
                      NH2_total: float) -> np.ndarray:
    """RHS of the crosslinking ODE system.

    State variables y = [X, G_free, NH2_free]:
        X: crosslink density [mol/m³]
        G_free: free genipin concentration [mol/m³]
        NH2_free: free amine concentration [mol/m³]

    dX/dt = k · G_free · NH2_free
    dG_free/dt = -k · G_free · NH2_free
    dNH2_free/dt = -2 · k · G_free · NH2_free

    Each crosslink bridge consumes 1 genipin and 2 amines (one on each chain).
    """
    X, G_free, NH2_free = y
    # Ensure non-negative
    G_free = max(G_free, 0.0)
    NH2_free = max(NH2_free, 0.0)

    rate = k * G_free * NH2_free
    return np.array([rate, -rate, -2.0 * rate])


def crosslink_density_to_properties(
    X: float, c_chitosan: float, T: float, DDA: float, M_GlcN: float,
    NH2_total: float, f_bridge: float = 0.4,
) -> tuple[float, float, float, float]:
    """Convert crosslink density to network properties.

    Parameters
    ----------
    f_bridge : float
        Fraction of genipin-amine reactions that produce elastically
        active inter-chain crosslinks (not pendant modifications or
        intra-chain loops).  Literature range 0.3-0.5 for chitosan-
        genipin systems.

    Returns (nu_e, Mc, xi, G_chitosan):
        nu_e: effective crosslink density [1/m³]
        Mc: molecular weight between crosslinks [g/mol]
        xi: mesh size [m]
        G_chitosan: shear modulus of chitosan network [Pa]
    """
    # Only a fraction f_bridge of reacted genipin produces elastically
    # active inter-chain crosslinks (case 3).  The rest are pendant
    # modifications (case 1) or intra-chain loops (case 2).
    X_effective = f_bridge * X

    # Effective crosslink density (assuming tetrafunctional crosslinks)
    # nu_e = X_eff · N_A  (elastically active crosslinks per m³)
    nu_e = X_effective * N_AVOGADRO

    # Molecular weight between crosslinks using effective conversion
    p_effective = X_effective / NH2_total if NH2_total > 0 else 0.0
    p_effective = min(p_effective, 1.0)

    # Mc = M_GlcN / (2·p_eff)  for random crosslinking on a linear chain
    if p_effective > 1e-10:
        Mc = M_GlcN / (2.0 * p_effective)
    else:
        Mc = 1e10  # effectively infinite

    # Mesh size from Canal-Peppas
    # Estimate polymer volume fraction in gel
    rho_polymer = 1400.0  # kg/m³ for chitosan
    nu_2s = c_chitosan / rho_polymer
    nu_2s = max(nu_2s, 0.001)

    if Mc > 0 and Mc < 1e8:
        xi_nm = 0.071 * nu_2s ** (-1.0 / 3.0) * np.sqrt(Mc)
        xi = xi_nm * 1e-9  # nm → m
    else:
        xi = 1e-6  # 1 µm default for uncrosslinked

    # Shear modulus from rubber elasticity: G = nu_e · kT
    G_chitosan = nu_e * K_BOLTZMANN * T

    return nu_e, Mc, xi, G_chitosan


def solve_crosslinking(params: SimulationParameters,
                       props: MaterialProperties,
                       R_droplet: Optional[float] = None,
                       porosity: Optional[float] = None) -> CrosslinkingResult:
    """Solve the crosslinking kinetics ODE system.

    Parameters
    ----------
    params : SimulationParameters
        Contains crosslinking temperature, time, genipin concentration.
    props : MaterialProperties
        Contains kinetic constants, chitosan properties.
    R_droplet : float, optional
        Microsphere radius [m] from L1 (d50/2).  Used to compute the
        Thiele modulus and validate the lumped-ODE assumption.
    porosity : float, optional
        Gel porosity from L2.  Adjusts effective genipin diffusivity
        via D_eff = D_genipin * porosity**2 (tortuosity ~ 1/porosity).
    """
    T = params.formulation.T_crosslink
    t_xlink = params.formulation.t_crosslink
    c_genipin = params.formulation.c_genipin  # mol/m³
    c_chitosan = params.formulation.c_chitosan  # kg/m³

    # Rate constant
    k = genipin_rate_constant(T, props.k_xlink_0, props.E_a_xlink)

    # Available amine groups
    NH2_total = available_amine_concentration(c_chitosan, props.DDA, props.M_GlcN)

    # ── Thiele modulus check (couples L3 to L1/L2 outputs) ───────────
    D_genipin_bare = 1e-10  # m²/s, genipin in water
    if R_droplet is not None and R_droplet > 0:
        if porosity is not None and porosity > 0:
            # Effective diffusivity: D_eff = D_bare * porosity²
            # (tortuosity ~ 1/porosity for random porous media)
            D_eff = D_genipin_bare * porosity ** 2
        else:
            D_eff = D_genipin_bare

        thiele = compute_thiele_modulus(R_droplet, k, NH2_total, D_eff)
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

    # Initial conditions: X=0, G_free=c_genipin, NH2_free=NH2_total
    y0 = np.array([0.0, c_genipin, NH2_total])

    # Time span
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

    # Extract time histories
    t_arr = sol.t
    X_arr = sol.y[0]

    # Compute derived properties at each time point
    nu_e_arr = np.zeros_like(t_arr)
    Mc_arr = np.zeros_like(t_arr)
    xi_arr = np.zeros_like(t_arr)
    G_chit_arr = np.zeros_like(t_arr)

    for i in range(len(t_arr)):
        nu_e, Mc, xi, G_c = crosslink_density_to_properties(
            X_arr[i], c_chitosan, T, props.DDA, props.M_GlcN, NH2_total,
            props.f_bridge,
        )
        nu_e_arr[i] = nu_e
        Mc_arr[i] = Mc
        xi_arr[i] = xi
        G_chit_arr[i] = G_c

    # Final values
    p_final = X_arr[-1] / NH2_total if NH2_total > 0 else 0.0

    return CrosslinkingResult(
        t_array=t_arr,
        X_array=X_arr,
        nu_e_array=nu_e_arr,
        Mc_array=Mc_arr,
        xi_array=xi_arr,
        G_chitosan_array=G_chit_arr,
        p_final=float(min(p_final, 1.0)),
        nu_e_final=float(nu_e_arr[-1]),
        Mc_final=float(Mc_arr[-1]),
        xi_final=float(xi_arr[-1]),
        G_chitosan_final=float(G_chit_arr[-1]),
    )
