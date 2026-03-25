"""Level 3: Genipin-chitosan crosslinking kinetics ODE solver.

Models the evolution of crosslink density, mesh size, and modulus
during genipin crosslinking of the chitosan network within the
already-gelled microsphere.

Thiele modulus << 1 for 2 µm spheres → diffusion NOT limiting →
uniform crosslinking → simple ODE system (no spatial PDE needed).
"""

from __future__ import annotations

import numpy as np
from scipy.integrate import solve_ivp

from ..datatypes import CrosslinkingResult, MaterialProperties, SimulationParameters

# Constants
R_GAS = 8.314       # J/(mol·K)
K_BOLTZMANN = 1.38e-23  # J/K
N_AVOGADRO = 6.022e23


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
    dNH2_free/dt = -k · G_free · NH2_free

    Each crosslink consumes 1 genipin and 1 amine.
    """
    X, G_free, NH2_free = y
    # Ensure non-negative
    G_free = max(G_free, 0.0)
    NH2_free = max(NH2_free, 0.0)

    rate = k * G_free * NH2_free
    return np.array([rate, -rate, -rate])


def crosslink_density_to_properties(
    X: float, c_chitosan: float, T: float, DDA: float, M_GlcN: float,
    NH2_total: float,
) -> tuple[float, float, float, float]:
    """Convert crosslink density to network properties.

    Returns (nu_e, Mc, xi, G_chitosan):
        nu_e: effective crosslink density [1/m³]
        Mc: molecular weight between crosslinks [g/mol]
        xi: mesh size [m]
        G_chitosan: shear modulus of chitosan network [Pa]
    """
    # Crosslinking fraction
    p = X / NH2_total if NH2_total > 0 else 0.0
    p = min(p, 1.0)

    # Effective crosslink density (assuming tetrafunctional crosslinks)
    # nu_e = X · N_A  (crosslinks per m³)
    nu_e = X * N_AVOGADRO

    # Molecular weight between crosslinks
    # Mc = M_GlcN / (2·p)  for random crosslinking on a linear chain
    if p > 1e-10:
        Mc = M_GlcN / (2.0 * p)
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
                       props: MaterialProperties) -> CrosslinkingResult:
    """Solve the crosslinking kinetics ODE system.

    Parameters
    ----------
    params : SimulationParameters
        Contains crosslinking temperature, time, genipin concentration.
    props : MaterialProperties
        Contains kinetic constants, chitosan properties.
    """
    T = params.formulation.T_crosslink
    t_xlink = params.formulation.t_crosslink
    c_genipin = params.formulation.c_genipin  # mol/m³
    c_chitosan = params.formulation.c_chitosan  # kg/m³

    # Rate constant
    k = genipin_rate_constant(T, props.k_xlink_0, props.E_a_xlink)

    # Available amine groups
    NH2_total = available_amine_concentration(c_chitosan, props.DDA, props.M_GlcN)

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
