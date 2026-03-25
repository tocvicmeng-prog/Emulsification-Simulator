"""Core data structures for the emulsification simulation pipeline."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

import numpy as np


# ─── Simulation Parameters ───────────────────────────────────────────────

@dataclass
class MixerGeometry:
    """Rotor-stator mixer geometry."""
    rotor_diameter: float = 0.025       # [m] (25 mm default)
    stator_diameter: float = 0.026      # [m]
    gap_width: float = 0.0005           # [m] (0.5 mm)
    tank_volume: float = 0.0005         # [m³] (500 mL)
    power_number: float = 1.5           # [-]
    dissipation_ratio: float = 50.0     # ε_max / ε_avg


@dataclass
class EmulsificationParameters:
    """Process parameters for Level 1."""
    rpm: float = 10000.0                # [rev/min]
    t_emulsification: float = 600.0     # [s]
    mixer: MixerGeometry = field(default_factory=MixerGeometry)


@dataclass
class FormulationParameters:
    """Chemical formulation parameters."""
    c_agarose: float = 42.0             # [kg/m³] (4.2% w/v)
    c_chitosan: float = 18.0            # [kg/m³] (1.8% w/v)
    c_span80: float = 20.0              # [kg/m³] (~2% w/v)
    c_genipin: float = 2.0              # [mol/m³] (~2 mM)
    T_oil: float = 363.15              # [K] (90°C)
    cooling_rate: float = 0.167         # [K/s] (~10°C/min)
    T_crosslink: float = 310.15         # [K] (37°C)
    t_crosslink: float = 86400.0        # [s] (24 hours)
    phi_d: float = 0.05                 # [-] dispersed phase volume fraction

    @property
    def agarose_fraction(self) -> float:
        """Agarose mass fraction in polymer blend."""
        total = self.c_agarose + self.c_chitosan
        if total == 0:
            return 0.0
        return self.c_agarose / total

    @property
    def total_polymer(self) -> float:
        """Total polymer concentration [kg/m³]."""
        return self.c_agarose + self.c_chitosan


@dataclass
class SolverSettings:
    """Numerical solver settings."""
    # Level 1
    l1_n_bins: int = 50
    l1_d_min: float = 0.1e-6           # [m]
    l1_d_max: float = 500e-6           # [m]
    l1_rtol: float = 1e-6
    l1_atol: float = 1e-8
    # Level 2
    l2_n_r: int = 1000
    l2_n_grid: int = 128               # 2D grid side length (N×N)
    l2_dt_initial: float = 1e-4        # [s]
    l2_dt_max: float = 1.0             # [s]
    l2_arrest_exponent: float = 2.5
    # Level 3
    l3_method: str = "Radau"
    l3_rtol: float = 1e-8
    l3_atol: float = 1e-10


@dataclass
class SimulationParameters:
    """Top-level parameter container for the full pipeline."""
    emulsification: EmulsificationParameters = field(default_factory=EmulsificationParameters)
    formulation: FormulationParameters = field(default_factory=FormulationParameters)
    solver: SolverSettings = field(default_factory=SolverSettings)
    run_id: str = ""
    notes: str = ""

    def to_optimization_vector(self) -> np.ndarray:
        """Flatten to 7D vector for BoTorch."""
        return np.array([
            self.emulsification.rpm,
            self.formulation.c_span80,
            self.formulation.agarose_fraction,
            self.formulation.T_oil,
            self.formulation.cooling_rate,
            self.formulation.c_genipin,
            self.formulation.t_crosslink,
        ])

    @classmethod
    def from_optimization_vector(
        cls, x: np.ndarray, template: SimulationParameters
    ) -> SimulationParameters:
        """Reconstruct from 7D vector + template for fixed params."""
        import copy
        params = copy.deepcopy(template)
        params.emulsification.rpm = float(x[0])
        params.formulation.c_span80 = float(x[1])
        # Reconstruct agarose/chitosan from fraction and total
        frac = float(x[2])
        total = params.formulation.total_polymer
        params.formulation.c_agarose = frac * total
        params.formulation.c_chitosan = (1.0 - frac) * total
        params.formulation.T_oil = float(x[3])
        params.formulation.cooling_rate = float(x[4])
        params.formulation.c_genipin = float(x[5])
        params.formulation.t_crosslink = float(x[6])
        return params


# ─── Material Properties ─────────────────────────────────────────────────

@dataclass
class PropertyValue:
    """A single material property with metadata."""
    value: float
    unit: str
    uncertainty: float = 0.0
    source: str = ""
    T_ref: float = 298.15              # [K]


@dataclass
class MaterialProperties:
    """Aggregated material properties for the simulation."""
    # Oil phase
    rho_oil: float = 850.0              # [kg/m³] at 20°C reference (interpolated to T_oil by PropertyDatabase)
    mu_oil: float = 0.005               # [Pa·s] at 90°C reference (interpolated to T_oil by PropertyDatabase)

    # Aqueous / dispersed phase
    rho_aq: float = 1020.0              # [kg/m³]
    mu_d: float = 1.0                   # [Pa·s] dispersed phase viscosity at T_oil

    # Interfacial
    sigma: float = 5.0e-3               # [N/m] interfacial tension with Span-80

    # Thermodynamic
    chi_0: float = 0.497                # Flory-Huggins χ at reference T
    chi_T_coeffs: tuple = (505.0, -0.891)  # (A, B) for χ(T) = A/T + B; spinodal onset ~330 K (57°C)
    kappa_CH: float = 5.0e-12           # [J/m] Cahn-Hilliard gradient coefficient
    M_0: float = 1.0e-9                # [m⁵/(J·s)] bare mobility (calibrated for 50-100 nm coarsening)

    # Gelation
    T_gel: float = 311.15               # [K] (~38°C)
    k_gel_0: float = 1.0                # [1/s] Avrami rate prefactor
    n_avrami: float = 2.5               # [-]
    gel_arrest_exponent: float = 2.5    # β

    # Crosslinking
    k_xlink_0: float = 2806.0          # [m³/(mol·s)] calibrated: k(37°C)=5e-3 L/(mol·s)
    E_a_xlink: float = 52000.0         # [J/mol]
    DDA: float = 0.90                   # degree of deacetylation
    M_GlcN: float = 161.16             # [g/mol] glucosamine molar mass
    M_genipin: float = 226.23          # [g/mol]

    # Agarose gel
    G_agarose_prefactor: float = 3000.0  # [Pa] at 1% w/v
    G_agarose_exponent: float = 2.2      # power law exponent

    # Crosslinking bridge efficiency
    f_bridge: float = 0.4              # fraction of genipin reactions producing elastically active crosslinks

    # IPN coupling
    eta_coupling: float = -0.15        # IPN coupling coefficient

    # Network / pore
    r_fiber: float = 1.5e-9            # [m] agarose fiber radius


# ─── Result Structures ────────────────────────────────────────────────────

@dataclass
class EmulsificationResult:
    """Output of Level 1: PBE solver."""
    d_bins: np.ndarray                  # [m] pivot diameters (N_bins,)
    n_d: np.ndarray                     # [#/m³] number density (N_bins,)
    d32: float                          # [m] Sauter mean diameter
    d43: float                          # [m] volume-weighted mean
    d10: float                          # [m] 10th percentile
    d50: float                          # [m] median
    d90: float                          # [m] 90th percentile
    span: float                         # [-] (d90 - d10) / d50
    total_volume_fraction: float        # [-]
    converged: bool
    t_history: Optional[np.ndarray] = None
    n_d_history: Optional[np.ndarray] = None


@dataclass
class GelationResult:
    """Output of Level 2: Phase-field solver.

    Supports both 1D radial (legacy) and 2D Cartesian solvers.
    For 1D: r_grid is (N_r,), phi_field is (N_r,).
    For 2D: r_grid is (N,) coordinate array, phi_field is (N, N).
    """
    r_grid: np.ndarray                  # [m] (N_r,) or (N,) for 2D coords
    phi_field: np.ndarray               # [-] (N_r,) or (N, N)
    pore_size_mean: float               # [m]
    pore_size_std: float                # [m]
    pore_size_distribution: np.ndarray  # [m]
    porosity: float                     # [-]
    alpha_final: float                  # [-]
    char_wavelength: float              # [m]
    T_history: Optional[np.ndarray] = None
    phi_snapshots: Optional[np.ndarray] = None
    L_domain: float = 0.0              # [m] domain side length (2D solver)
    grid_spacing: float = 0.0          # [m] uniform grid spacing


@dataclass
class CrosslinkingResult:
    """Output of Level 3: ODE kinetics."""
    t_array: np.ndarray                 # [s] (N_t,)
    X_array: np.ndarray                 # [mol/m³] (N_t,)
    nu_e_array: np.ndarray              # [1/m³] (N_t,)
    Mc_array: np.ndarray                # [g/mol] (N_t,)
    xi_array: np.ndarray                # [m] (N_t,)
    G_chitosan_array: np.ndarray        # [Pa] (N_t,)
    p_final: float                      # [-]
    nu_e_final: float                   # [1/m³]
    Mc_final: float                     # [g/mol]
    xi_final: float                     # [m]
    G_chitosan_final: float             # [Pa]


@dataclass
class MechanicalResult:
    """Output of Level 4: Property prediction."""
    G_agarose: float                    # [Pa]
    G_chitosan: float                   # [Pa]
    G_DN: float                         # [Pa]
    E_star: float                       # [Pa]
    delta_array: np.ndarray             # [m]
    F_array: np.ndarray                 # [N]
    rh_array: np.ndarray                # [m]
    Kav_array: np.ndarray               # [-]
    pore_size_mean: float               # [m]
    xi_mesh: float                      # [m]


@dataclass
class FullResult:
    """Complete pipeline output."""
    parameters: SimulationParameters
    emulsification: EmulsificationResult
    gelation: GelationResult
    crosslinking: CrosslinkingResult
    mechanical: MechanicalResult

    def objective_vector(self) -> np.ndarray:
        """Compute the 3 objective values for optimization.

        Delegates to the canonical implementation in
        ``emulsim.optimization.objectives.compute_objectives``.
        """
        from emulsim.optimization.objectives import compute_objectives
        return compute_objectives(self)


@dataclass
class OptimizationState:
    """State of the Bayesian optimization campaign."""
    X_observed: np.ndarray              # (N_eval, 7)
    Y_observed: np.ndarray              # (N_eval, 3)
    pareto_X: np.ndarray
    pareto_Y: np.ndarray
    iteration: int
    hypervolume: float
    hypervolume_history: list = field(default_factory=list)
    converged: bool = False
    gp_state: Optional[dict] = None
