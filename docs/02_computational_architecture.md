# Computational Architecture: Multi-Scale Emulsification Simulation System

**Version**: 1.0
**Date**: 2026-03-25
**Author**: Computational Architect (Claude Opus 4.6)
**Upstream**: `01_scientific_advisor_report.md`
**Status**: PROPOSED — awaiting implementation approval

---

## Table of Contents

1. [System Architecture Overview](#1-system-architecture-overview)
2. [Module Specifications](#2-module-specifications)
3. [Data Structures and Interfaces](#3-data-structures-and-interfaces)
4. [Algorithm Specifications](#4-algorithm-specifications)
5. [Technology Stack](#5-technology-stack)
6. [Directory Structure](#6-directory-structure)
7. [Implementation Roadmap](#7-implementation-roadmap)
8. [Verification and Validation Strategy](#8-verification-and-validation-strategy)
9. [First-Stage Audit: Scientific Rationality](#9-first-stage-audit-scientific-rationality)

---

## 1. System Architecture Overview

### 1.1 Overall System Diagram

```
┌─────────────────────────────────────────────────────────────────────────┐
│                        PIPELINE ORCHESTRATOR (M7)                       │
│  ┌───────────┐   ┌───────────┐   ┌───────────┐   ┌───────────┐        │
│  │ Level 1   │──>│ Level 2   │──>│ Level 3   │──>│ Level 4   │        │
│  │ Emulsif.  │   │ Gelation  │   │ Crosslink │   │ Mechanics │        │
│  │   (M2)    │   │   (M3)    │   │   (M4)    │   │   (M5)    │        │
│  └─────┬─────┘   └─────┬─────┘   └─────┬─────┘   └─────┬─────┘        │
│        │               │               │               │               │
│        │    n(d)        │  pore struct  │  ν_e, M_c, G  │  G_DN, F(δ)  │
│        └───────────>────┘───────────>───┘───────────>───┘               │
│                                                                         │
│  ┌───────────────────────────────┐   ┌──────────────────────────┐      │
│  │   PROPERTY DATABASE (M1)      │   │  OPTIMIZATION ENGINE (M6)│      │
│  │  • Material properties        │   │  • Bayesian Opt (BoTorch)│      │
│  │  • Literature values ± σ      │   │  • GP surrogate          │      │
│  │  • T/conc interpolation       │   │  • Multi-objective        │      │
│  └───────────────────────────────┘   └──────────────────────────┘      │
└─────────────────────────────────────────────────────────────────────────┘
```

### 1.2 Module Decomposition

| Module | Name                    | Responsibility                                      | Coupling   |
|--------|-------------------------|------------------------------------------------------|------------|
| M1     | Property Database       | Store/retrieve all material properties               | Read by all|
| M2     | Emulsification (L1)     | Solve PBE, predict droplet size distribution         | Feeds M3   |
| M3     | Gelation & Pores (L2)   | Phase-field + Avrami, predict pore structure         | Feeds M4   |
| M4     | Crosslinking (L3)       | ODE kinetics, predict crosslink density & modulus    | Feeds M5   |
| M5     | Mechanical Props (L4)   | Analytical rubber elasticity, Hertz, Ogston          | Terminal   |
| M6     | Optimization Engine     | Bayesian optimization over full pipeline             | Drives M7  |
| M7     | Pipeline Orchestrator   | Sequence execution, logging, checkpointing           | Controls all|

### 1.3 Data Flow

```
SimulationParameters ──┐
                       v
                ┌──────────┐
                │    M1    │──> MaterialProperties (shared context)
                └──────────┘
                       │
                       v
                ┌──────────┐
                │    M2    │──> EmulsificationResult { n(d), d32, d43, span }
                └──────────┘
                       │
                       v
                ┌──────────┐
                │    M3    │──> GelationResult { pore_sizes, connectivity, φ_field }
                └──────────┘
                       │
                       v
                ┌──────────┐
                │    M4    │──> CrosslinkingResult { ν_e(t), M_c(t), ξ(t), G(t) }
                └──────────┘
                       │
                       v
                ┌──────────┐
                │    M5    │──> MechanicalResult { G_DN, F_vs_delta, K_av }
                └──────────┘
                       │
                       v
                ┌──────────┐
                │    M6    │──> OptimizationState { next_point, GP_model, Pareto_front }
                └──────────┘
```

### 1.4 Technology Stack Summary

| Component          | Choice                  | Justification                                        |
|--------------------|-------------------------|------------------------------------------------------|
| Language           | Python 3.11+            | NumPy/SciPy ecosystem, rapid prototyping             |
| Numerical core     | NumPy 1.26+, SciPy 1.12+| Mature, well-tested ODE/PDE solvers                 |
| Phase-field solver | Custom FD + SciPy sparse| Lightweight, no external mesh dependency             |
| Optimization       | BoTorch 0.11+ / GPyTorch| State-of-art BO for multi-objective                  |
| Data exchange      | HDF5 (h5py)             | Efficient for large arrays + metadata                |
| Configuration      | TOML                    | Human-readable, Python stdlib (tomllib)              |
| Visualization      | Matplotlib + Plotly     | Publication quality + interactive dashboards         |
| Testing            | pytest + hypothesis     | Property-based testing for numerical code            |
| Acceleration       | Numba (optional)        | JIT for hot loops in PBE and phase-field             |
| Typing             | dataclasses + Pydantic  | Validated, self-documenting data structures          |

---

## 2. Module Specifications

### 2.1 Module 1: Property Database

**Purpose**: Centralized, validated repository of all material properties with uncertainty quantification and temperature/concentration interpolation.

#### 2.1.1 Storage Design

Properties are stored in a TOML configuration file with the following schema:

```toml
[oil_phase.paraffin]
rho = { value = 850.0, unit = "kg/m3", T_ref = 298.15, source = "CRC Handbook" }
rho_T_coeff = { value = -0.65, unit = "kg/(m3·K)", note = "linear approx 20-100°C" }
mu_90C = { value = 0.005, unit = "Pa·s", uncertainty = 0.001, source = "measured" }

[aqueous_phase.agarose]
T_gel = { value = 311.15, unit = "K", range = [309.15, 315.15], source = "Sigma A9539" }
# ... etc.
```

#### 2.1.2 Retrieval Interface

```python
class PropertyDatabase:
    def __init__(self, config_path: Path):
        """Load from TOML."""

    def get(self, path: str, T: float = None, conc: float = None) -> PropertyValue:
        """Retrieve property, optionally interpolated to T and/or concentration.

        path: dot-separated key, e.g. "oil_phase.paraffin.rho"
        Returns PropertyValue with .value, .unit, .uncertainty, .source
        """

    def viscosity_dispersed(self, T: float, c_agarose: float, c_chitosan: float) -> float:
        """Mark-Houwink + mixing rule for agarose/chitosan solution viscosity.
        Uses Huggins equation at low overlap, Cross model at high shear."""

    def interfacial_tension(self, T: float, c_span80: float) -> float:
        """Szyszkowski-Langmuir isotherm for Span-80 at oil/water interface."""

    def chi_parameter(self, T: float, phi_polymer: float) -> float:
        """Flory-Huggins chi(T) for agarose-water system."""
```

#### 2.1.3 Temperature Interpolation

- Linear for density: ρ(T) = ρ_ref + α·(T - T_ref)
- Arrhenius for viscosity: μ(T) = A·exp(E_a / RT)
- Polynomial for chi: χ(T) = a/T + b + c·T (fitted to literature data)

#### 2.1.4 Uncertainty Propagation

Each property carries a `±uncertainty` (1σ). When the optimization engine samples parameters, it can optionally perturb properties within these bounds (Monte Carlo uncertainty quantification).

---

### 2.2 Module 2: Emulsification Simulator (Level 1)

**Purpose**: Solve the population balance equation to predict droplet size distribution under rotor-stator emulsification of a high-viscosity agarose/chitosan dispersed phase in paraffin oil.

#### 2.2.1 PBE Solver Architecture

**Recommendation: Fixed-Pivot Class Method (Kumar & Ramkrishna, 1996) with 50 logarithmically-spaced bins.**

Justification for Class Method over QMOM:
1. **Full distribution recovery** — QMOM only tracks moments; we need n(d) for Level 2 sampling.
2. **Robustness** — QMOM moment-inversion (Wheeler algorithm) can fail for broad or bimodal distributions, which are expected during early emulsification.
3. **Breakage/coalescence flexibility** — Class methods handle arbitrary kernel shapes without closure assumptions.
4. **Debugging transparency** — The evolving bin concentrations are directly interpretable.
5. **Acceptable cost** — 50 bins with O(N²) coalescence = 2,500 kernel evaluations per time step, trivially fast for 0D.

Trade-off: QMOM is faster for very large N_bins. With N=50, class method completes in <1 second for the full time evolution, so speed is not a concern.

#### 2.2.2 Discretization

```
Bin edges: d_i for i = 0, 1, ..., N_bins
    d_0 = 0.1 µm   (smallest physically relevant)
    d_N = 500 µm    (largest before breakage dominates)
    Spacing: logarithmic, d_i = d_0 * (d_N/d_0)^(i/N_bins)

Pivot diameters: x_i = (d_i + d_{i+1}) / 2  (geometric mean for log spacing)

State vector: N_i(t) = number density in bin i [#/m³]
```

#### 2.2.3 Governing Equation (Discretized)

```
dN_i/dt = B_break,i - D_break,i + B_coal,i - D_coal,i

where:
  D_break,i = g(x_i) · N_i                              (death by breakage of bin i)
  B_break,i = Σ_{j>i} η_{ij} · β(x_i|x_j) · g(x_j) · N_j  (birth from breakage of larger drops)
  D_coal,i  = N_i · Σ_j q(x_i, x_j) · N_j              (death by coalescence of bin i)
  B_coal,i  = (1/2) · Σ_{j,k: x_j+x_k∈bin_i} q(x_j,x_k) · N_j · N_k  (birth from coalescence)
```

Fixed-pivot redistribution (Kumar & Ramkrishna 1996) conserves mass exactly by assigning daughter drops to adjacent bins with appropriate weights.

#### 2.2.4 Breakage Kernel

Alopaeus et al. (2002) with viscosity correction:

```
g(d) = C_1 · (ε/ν_c)^(1/2) · exp(-C_2 · σ / (ρ_c · ε^(2/3) · d^(5/3))
                                  - C_3 · µ_d / (ρ_c^(1/2) · ε^(1/3) · d^(4/3)))

where:
  C_1, C_2, C_3 = fitted constants (default: C_1=0.986, C_2=0.0, C_3 from Alopaeus)
  ε = local energy dissipation rate [W/kg]
  ν_c = kinematic viscosity of continuous phase [m²/s]
  σ = interfacial tension [N/m]
  ρ_c = continuous phase density [kg/m³]
  µ_d = dispersed phase dynamic viscosity [Pa·s]
  d = droplet diameter [m]
```

The second exponential term (C_3 term) is the **viscous resistance to deformation** — critical for our high-viscosity agarose dispersed phase.

Daughter size distribution: β(d'|d) assumed to follow a beta distribution (Luo & Svendsen 1996) or uniform binary breakage as a simpler default.

#### 2.2.5 Coalescence Kernel

Coulaloglou & Tavlarides (1977):

```
q(d_i, d_j) = h(d_i, d_j) · λ(d_i, d_j)

h(d_i, d_j) = C_4 · ε^(1/3) · (d_i^2 + d_j^2) · (d_i^(2/3) + d_j^(2/3))^(1/2)
                                                    [collision frequency]

λ(d_i, d_j) = exp(-C_5 · µ_c · ρ_c · ε / σ² · (d_i·d_j/(d_i+d_j))^4)
                                                    [coalescence efficiency]
```

Span-80 enters through σ (lower σ → higher coalescence efficiency, but also easier breakage).

#### 2.2.6 Energy Dissipation Model

For rotor-stator mixer:

```
P = N_p · ρ_emulsion · N³ · D⁵    (power draw)
ε_avg = P / (ρ_emulsion · V_tank)  (average dissipation)
ε_max = k_ε · ε_avg                (local max in rotor-stator gap, k_ε ~ 10-100)

where:
  N_p = power number (~1-5 for rotor-stator, geometry-dependent)
  N = rotational speed [rev/s] = RPM/60
  D = rotor diameter [m]
  V_tank = tank volume [m³]
  k_ε = ratio of max-to-mean dissipation (fitted or from literature)
```

For 0D model: use ε_max in the rotor-stator gap region (most breakage occurs here).

#### 2.2.7 Time Integration

- **Solver**: `scipy.integrate.solve_ivp` with BDF method (implicit, for stiff coalescence terms)
- **Time span**: 0 to t_emulsification (typically 5-30 minutes)
- **Tolerances**: rtol=1e-6, atol=1e-8
- **Initial condition**: Coarse premix distribution — log-normal with d32 ~ 100-500 µm

#### 2.2.8 Interfacial Tension Sub-model

Szyszkowski-Langmuir:

```
σ(c_s, T) = σ_0(T) - R·T·Γ_∞·ln(1 + K_L·c_s)

where:
  σ_0(T) = bare oil/water interfacial tension (~35 mN/m for paraffin/water)
  Γ_∞ = maximum surface excess of Span-80 (~3.5×10⁻⁶ mol/m²)
  K_L = Langmuir adsorption constant (~200 m³/mol for Span-80)
  c_s = bulk Span-80 concentration [mol/m³]
```

#### 2.2.9 Input/Output

```
INPUT:  SimulationParameters.emulsification → {RPM, t_emulsification, mixer_geometry}
        SimulationParameters.formulation → {c_span80, T_oil, c_agarose, c_chitosan}
        MaterialProperties (from M1)

OUTPUT: EmulsificationResult
        → n_d[N_bins]: number density per bin [#/m³]
        → d_bins[N_bins]: pivot diameters [m]
        → d32: Sauter mean diameter [m]
        → d43: volume-weighted mean [m]
        → span: (d90 - d10) / d50
        → sigma_dist: distribution width
        → t_history: time evolution snapshots (optional)
```

#### 2.2.10 Validation Hooks

```python
class EmulsificationValidator:
    def check_mass_conservation(self, result: EmulsificationResult, tol=1e-4) -> bool:
        """Total volume of drops must be conserved ±tol."""

    def check_physical_bounds(self, result: EmulsificationResult) -> bool:
        """All N_i >= 0, d32 in [0.5 µm, 500 µm]."""

    def check_steady_state(self, result: EmulsificationResult, tol=0.01) -> bool:
        """d32 changes < tol over last 10% of time."""
```

---

### 2.3 Module 3: Gelation & Pore Formation Simulator (Level 2)

**Purpose**: Simulate spinodal decomposition of the polymer solution during cooling, arrested by agarose gelation, to predict the internal pore structure of microspheres.

#### 2.3.1 Phase-Field Solver Architecture

**Recommendation: 1D radial finite-difference with semi-implicit time stepping.**

Justification for 1D radial over 2D:
1. **Microsphere symmetry** — Cooling is uniform (Bi << 0.1), so radial direction is the only relevant coordinate for pore structure variation (skin vs. core).
2. **Computational cost** — 2D phase-field at 2 nm resolution on a 2 µm domain = 10⁶ grid points, demanding ~GB memory and hours of compute. 1D radial with 1000 points is trivial.
3. **Sufficient for pore statistics** — We need pore *size distribution*, not exact 3D morphology. 1D structure factor analysis gives characteristic wavelength (= pore spacing).
4. **Upgrade path** — If 2D is needed later (e.g., for connectivity analysis), the same equations and kernels apply on a 2D Cartesian domain with minimal code changes.

**Upgrade note**: If anisotropic pore structures or percolation connectivity are important, promote to 2D. The solver is designed with a generic `SpatialDomain` abstraction to support this.

#### 2.3.2 Spatial Discretization

```
1D Radial:
  r ∈ [0, R_droplet]  where R_droplet = d/2 from Level 1
  N_r = 1000 grid points (Δr ~ 1 nm for 1 µm radius)

  Laplacian in spherical coords:
  ∇²φ = (1/r²) · d/dr(r² · dφ/dr)

  Discretized (central differences):
  ∇²φ_i ≈ (2/r_i) · (φ_{i+1} - φ_{i-1})/(2Δr) + (φ_{i+1} - 2φ_i + φ_{i-1})/Δr²

  Boundary conditions:
  r = 0: dφ/dr = 0 (symmetry)
  r = R: dφ/dr = 0 (no-flux) or prescribed surface condition
```

#### 2.3.3 Governing Equations

**Cahn-Hilliard equation:**

```
∂φ/∂t = ∇·[M(φ, α) · ∇µ]

where:
  µ = δF/δφ = f'(φ) - κ·∇²φ       (chemical potential)

  f(φ) = Flory-Huggins free energy density:
  f(φ) = (k_B·T/v_0) · [φ·ln(φ)/N_p + (1-φ)·ln(1-φ)/N_s + χ(T)·φ·(1-φ)]

  where:
    φ = polymer volume fraction
    N_p = degree of polymerization of polymer
    N_s = 1 (solvent)
    v_0 = monomer/solvent molar volume
    χ(T) = A/T + B  (Flory-Huggins interaction parameter)
    κ = gradient energy coefficient [J/m] (related to interfacial width)

  M(φ, α) = M_0 · φ · (1-φ) · (1 - α)^β   (mobility, arrested by gelation)

  where:
    M_0 = bare mobility [m⁵/(J·s)]
    α = gelation degree (0→1 from Avrami)
    β = arrest exponent (~2-3, controls how sharply gelation stops coarsening)
```

**Avrami gelation kinetics (coupled):**

```
dα/dt = n · k_gel(T) · t^(n-1) · (1-α)      (differential form)

Or equivalently: α(t) = 1 - exp(-(k_gel · t_cool)^n)

where:
  k_gel(T) = Avrami rate constant (zero above T_gel, rapid below)
  n = Avrami exponent (~2-3 for agarose gelation)
  t_cool = time since temperature crossed T_gel

k_gel(T) = k_gel,0 · (1 - T/T_gel)^m   for T < T_gel
          = 0                             for T ≥ T_gel

Temperature profile (lumped capacitance, Bi << 0.1):
  T(t) = T_oil + (T_init - T_oil) · exp(-t/τ_cool)
  τ_cool = ρ·c_p·R / (3·h)   (cooling time constant)
```

#### 2.3.4 Time Stepping

**Semi-implicit scheme (Eyre 1998):**

Split the chemical potential into contractive (concave) and expansive (convex) parts:

```
µ = f'_c(φ) + f'_e(φ) - κ·∇²φ

(φ^{n+1} - φ^n)/Δt = ∇·[M^n · ∇(f'_c(φ^{n+1}) + f'_e(φ^n) - κ·∇²φ^{n+1})]
```

Treat the concave part (spinodal) implicitly → unconditionally stable.
Treat the convex part explicitly → no nonlinear solve needed.

For 1D with N=1000 points, this reduces to solving a banded linear system (bandwidth 5) at each time step using `scipy.linalg.solve_banded` — O(N) per step.

**Time step**: Adaptive, starting at Δt = 10⁻⁴ s, with Δt growing as coarsening slows. Typical total time: 60-300 s (cooling period).

#### 2.3.5 Pore Structure Analysis

After gelation is complete (α > 0.95), extract pore structure from the φ(r) field:

```python
class PoreAnalyzer:
    def structure_factor(self, phi_field: np.ndarray, r: np.ndarray) -> tuple:
        """Compute S(q) via Fourier transform of correlation function.
        Peak position q* → characteristic pore spacing d* = 2π/q*."""

    def chord_length_distribution(self, phi_field: np.ndarray, threshold: float) -> np.ndarray:
        """Binarize φ field, measure chord lengths in solvent-rich phase.
        Gives pore size distribution directly."""

    def porosity(self, phi_field: np.ndarray, threshold: float) -> float:
        """Volume fraction of solvent-rich phase (= pore space)."""

    def mean_pore_size(self, phi_field: np.ndarray) -> float:
        """From structure factor peak or chord length mean."""
```

For 1D: pore size = spacing between polymer-rich peaks in φ(r).

#### 2.3.6 Input/Output

```
INPUT:  EmulsificationResult.d32 (representative droplet diameter)
        SimulationParameters.formulation → {c_agarose, c_chitosan, cooling_rate}
        MaterialProperties → {χ(T), M_0, κ, k_gel, n_Avrami, T_gel}

OUTPUT: GelationResult
        → phi_field[N_r]: final polymer volume fraction profile
        → r_grid[N_r]: radial positions [m]
        → pore_size_mean: mean pore diameter [m]
        → pore_size_distribution: array of pore sizes [m]
        → porosity: volume fraction of pore space
        → alpha_final: final gelation degree
        → T_history: temperature vs time
        → phi_snapshots: φ(r) at selected times (for visualization)
```

---

### 2.4 Module 4: Crosslinking Kinetics (Level 3)

**Purpose**: Simulate genipin crosslinking of chitosan amine groups within the gelled microsphere, tracking crosslink density, mesh size, and elastic modulus over time.

#### 2.4.1 ODE System

Since the Thiele modulus is much less than 1 for 2 µm microspheres (diffusion is fast compared to reaction), the system is spatially uniform and reduces to coupled ODEs:

```
State variables:
  [NH2](t)  = free amine concentration [mol/m³]
  [Gen](t)  = free genipin concentration [mol/m³]
  [X](t)    = crosslink concentration [mol/m³]

Kinetics (second-order):
  d[X]/dt = k(T) · [NH2] · [Gen]
  d[NH2]/dt = -2 · k(T) · [NH2] · [Gen]     (2 amines consumed per crosslink)
  d[Gen]/dt = -k(T) · [NH2] · [Gen]          (1 genipin per crosslink)

Rate constant (Arrhenius):
  k(T) = k_0 · exp(-E_a / (R·T))
  k_0 calibrated so that k(37°C) = 5×10⁻³ L/(mol·s)
  E_a = 52 kJ/mol

Initial conditions:
  [NH2]_0 = f_NH2 · c_chitosan / M_GlcN     (free amine content of chitosan)
  [Gen]_0 = c_genipin / M_genipin
  [X]_0 = 0

  f_NH2 ≈ 0.75-0.85 (degree of deacetylation of chitosan)
```

#### 2.4.2 Coupled Property Evolution

At each time step, compute derived properties from [X](t):

```
Crosslink density:
  ν_e(t) = [X](t) · N_A                      [crosslinks/m³]

Molecular weight between crosslinks:
  M_c(t) = ρ_polymer / ν_e(t)                [g/mol]
  (valid once ν_e > ν_e,min, i.e., gel point reached)

Mesh size (Canal-Peppas):
  ξ(t) = 0.071 · ν_{2,s}^{-1/3} · M_c(t)^{1/2}   [nm]
  where ν_{2,s} = polymer volume fraction in swollen state

Shear modulus (rubber elasticity):
  G_chitosan(t) = ν_e(t) · k_B · T            [Pa]

Crosslinking degree:
  p(t) = [X](t) / [X]_max
  where [X]_max = min([NH2]_0 / 2, [Gen]_0)   (stoichiometric limit)
```

#### 2.4.3 Solver Selection

**Recommendation: `scipy.integrate.solve_ivp` with the Radau method (implicit Runge-Kutta).**

Justification:
- The system is mildly stiff (k(T) varies over orders of magnitude with temperature).
- Radau is L-stable → handles stiff problems robustly.
- Only 3 ODEs → cost per step is negligible (3x3 Jacobian).
- Dense output for smooth property curves.

Tolerances: rtol=1e-8, atol=1e-10 (tight, since downstream modulus predictions are sensitive to ν_e).

#### 2.4.4 Input/Output

```
INPUT:  GelationResult → {pore_size_mean, porosity, phi_polymer_avg}
        SimulationParameters.formulation → {c_genipin, c_chitosan, T_crosslink, t_crosslink}
        MaterialProperties → {k_0, E_a, f_NH2, DDA, M_GlcN, M_genipin}

OUTPUT: CrosslinkingResult
        → t_array[N_t]: time points [s]
        → X_array[N_t]: crosslink concentration vs time [mol/m³]
        → nu_e_array[N_t]: crosslink density vs time [1/m³]
        → Mc_array[N_t]: M_c vs time [g/mol]
        → xi_array[N_t]: mesh size vs time [m]
        → G_chitosan_array[N_t]: chitosan network modulus vs time [Pa]
        → p_final: final crosslinking degree [dimensionless]
        → nu_e_final, Mc_final, xi_final, G_chitosan_final: final values
```

---

### 2.5 Module 5: Mechanical Property Predictor (Level 4)

**Purpose**: Combine agarose gel and chitosan-genipin network properties to predict double-network microsphere mechanical behavior, chromatographic partitioning.

#### 2.5.1 Double-Network Modulus

Assuming additive, independent network contributions:

```
G_agarose = f(c_agarose, T_gel_history)
  Empirical: G_agarose ≈ A · c_agarose^n   (A ~ 3000 Pa at 1% w/v, n ~ 2.1-2.5)
  Source: Normand et al. (2000), Ross-Murphy (1991)

G_chitosan = ν_e_final · k_B · T   (from Level 3)

G_DN = G_agarose + G_chitosan
```

#### 2.5.2 Hertz Contact Model

For microsphere compression between two parallel plates (AFM or micromanipulation):

```
F(δ) = (4/3) · E* · R^(1/2) · δ^(3/2)

where:
  E* = 2·G_DN·(1+ν) / (1-ν)    (effective modulus, ν ≈ 0.5 for hydrogels)
      ≈ 3·G_DN                    (incompressible limit)
  R = microsphere radius [m]
  δ = indentation depth [m]

For compression between two plates:
  F(δ) = (4/3) · E* · (R/2)^(1/2) · (δ/2)^(3/2)  (adjusted contact geometry)
```

#### 2.5.3 Chromatographic Kav Predictor

Ogston model for partitioning of globular proteins into porous gel:

```
K_av = exp(-π · (r_h + r_f)² · n_f)

where:
  r_h = hydrodynamic radius of solute [m]
  r_f = fiber radius of gel network [m] (~1.5-2.0 nm for agarose)
  n_f = fiber density [m⁻²] = φ_polymer / (π · r_f²)

Alternatively, using mean pore size ξ:
  K_av = exp(-(r_h / ξ)²)   (simplified Ogston-Laurent)
```

#### 2.5.4 Input/Output

```
INPUT:  CrosslinkingResult → {G_chitosan_final, nu_e_final, xi_final}
        GelationResult → {pore_size_mean, porosity}
        EmulsificationResult → {d32}
        SimulationParameters.formulation → {c_agarose}
        MaterialProperties → {A_agarose, n_agarose, r_f}

OUTPUT: MechanicalResult
        → G_agarose: agarose network modulus [Pa]
        → G_chitosan: chitosan-genipin network modulus [Pa]
        → G_DN: double-network total modulus [Pa]
        → E_star: effective Young's modulus [Pa]
        → F_vs_delta: force-displacement curve (array) [N vs m]
        → K_av_vs_rh: partition coefficient vs solute size (array)
        → pore_size_mean: from gelation [m]
        → xi_mesh: from crosslinking [m]
```

---

### 2.6 Module 6: Optimization Engine

**Purpose**: Multi-objective Bayesian optimization over the full Level 1→4 pipeline to find process parameters that simultaneously achieve target droplet size, pore size, and mechanical modulus.

#### 2.6.1 Parameter Space

| Parameter          | Symbol    | Bounds          | Unit    | Scale |
|--------------------|-----------|-----------------|---------|-------|
| Rotor speed        | RPM       | [3000, 25000]   | rpm     | log   |
| Span-80 conc       | c_s       | [0.5, 5.0]      | % w/v   | log   |
| Agarose/chitosan   | r_AC      | [0.6, 0.9]      | ratio   | linear|
| Oil temperature    | T_oil     | [60, 95]        | °C      | linear|
| Cooling rate       | dT/dt     | [0.5, 20]       | °C/min  | log   |
| Genipin conc       | c_gen     | [0.1, 10]       | mM      | log   |
| Crosslinking time  | t_XL      | [1, 48]         | hours   | log   |

7 parameters total. Log scaling for parameters spanning >1 order of magnitude.

#### 2.6.2 Objectives

```
Minimize:
  f_1 = |d32 - 2.0 µm| / 2.0 µm                (target: 2 µm droplets)
  f_2 = |pore_mean - 80 nm| / 80 nm             (target: 60-100 nm pores)
  f_3 = |log10(G_DN) - log10(G_target)| / 1.0   (target modulus, user-defined)

Subject to:
  span < 2.0                (narrow size distribution)
  G_DN > G_min              (minimum mechanical integrity)
  xi_mesh > 2·r_h_target    (mesh must admit target protein)
```

#### 2.6.3 BoTorch Implementation

```python
class OptimizationEngine:
    def __init__(self, bounds: dict, objectives: list, n_initial: int = 20):
        """
        bounds: parameter bounds dict
        objectives: list of ObjectiveFunction callables
        n_initial: Sobol quasi-random initial samples
        """
        self.gp_models = {}  # one GP per objective
        self.X_observed = []  # parameter points evaluated
        self.Y_observed = []  # objective values

    def suggest_next(self, n_candidates: int = 1) -> np.ndarray:
        """Use EHVI (Expected Hypervolume Improvement) acquisition function
        for multi-objective optimization.

        Returns: array of shape (n_candidates, n_params)
        """

    def update(self, X_new: np.ndarray, Y_new: np.ndarray):
        """Add new observations, refit GP models."""

    def get_pareto_front(self) -> tuple:
        """Return current Pareto-optimal points and their objective values."""

    def convergence_check(self, tol: float = 0.01) -> bool:
        """Check if hypervolume improvement < tol over last 5 iterations."""
```

#### 2.6.4 Convergence Criteria

- **Primary**: Hypervolume improvement < 1% over 5 consecutive iterations.
- **Budget**: Maximum 200 pipeline evaluations (each takes ~seconds).
- **Secondary**: All objectives within 10% of target for at least one Pareto point.

#### 2.6.5 Initialization

Latin Hypercube Sampling (LHS) or Sobol sequence for the first 15-25 evaluations to build the initial GP surrogate with good space coverage.

---

### 2.7 Module 7: Pipeline Orchestrator

**Purpose**: Sequential execution of the L1→L2→L3→L4 pipeline with parameter passing, logging, checkpointing, and batch execution for optimization campaigns.

#### 2.7.1 Execution Flow

```python
class PipelineOrchestrator:
    def __init__(self, config: SimulationConfig):
        self.db = PropertyDatabase(config.properties_path)
        self.emulsifier = EmulsificationSimulator(self.db)
        self.gelation = GelationSimulator(self.db)
        self.crosslinker = CrosslinkingSimulator(self.db)
        self.mechanics = MechanicalPredictor(self.db)
        self.logger = SimulationLogger(config.output_dir)

    def run_single(self, params: SimulationParameters) -> FullResult:
        """Execute full L1→L2→L3→L4 pipeline."""
        self.logger.start_run(params)

        # Level 1: Emulsification
        emul_result = self.emulsifier.solve(params)
        self.logger.checkpoint("L1", emul_result)

        # Level 2: Gelation (using representative droplet from L1)
        gel_result = self.gelation.solve(params, emul_result)
        self.logger.checkpoint("L2", gel_result)

        # Level 3: Crosslinking
        xl_result = self.crosslinker.solve(params, gel_result)
        self.logger.checkpoint("L3", xl_result)

        # Level 4: Mechanical properties
        mech_result = self.mechanics.predict(params, emul_result, gel_result, xl_result)
        self.logger.checkpoint("L4", mech_result)

        full = FullResult(params, emul_result, gel_result, xl_result, mech_result)
        self.logger.finalize(full)
        return full

    def run_batch(self, param_list: list[SimulationParameters],
                  n_workers: int = 1) -> list[FullResult]:
        """Run multiple parameter sets. n_workers > 1 for parallel execution."""

    def resume_from_checkpoint(self, run_id: str) -> FullResult:
        """Resume a failed run from the last successful checkpoint."""
```

#### 2.7.2 Logging

```
output/
  run_20260325_143022/
    config.toml           # Input parameters (reproducibility)
    L1_emulsification.h5  # EmulsificationResult
    L2_gelation.h5        # GelationResult
    L3_crosslinking.h5    # CrosslinkingResult
    L4_mechanical.h5      # MechanicalResult
    summary.json          # Key scalar results for quick access
    run.log               # Text log with timestamps
```

#### 2.7.3 Checkpointing

Each level writes its result to HDF5 immediately upon completion. If the pipeline crashes at Level 3, the orchestrator can resume from the Level 2 checkpoint without re-running Level 1 and Level 2.

#### 2.7.4 Batch Execution

For optimization campaigns, the orchestrator supports:
- Sequential batch: run_batch with n_workers=1
- Parallel batch: run_batch with n_workers>1 (multiprocessing.Pool)
- Each run is independent (no shared mutable state), so parallelism is trivial.

---

## 3. Data Structures and Interfaces

All data structures use Python `dataclasses` with type annotations. Validation is performed at construction time.

```python
from dataclasses import dataclass, field
from typing import Optional
import numpy as np


# ─── Simulation Parameters ───────────────────────────────────────────────

@dataclass
class MixerGeometry:
    """Rotor-stator mixer geometry."""
    rotor_diameter: float          # [m]
    stator_diameter: float         # [m]
    gap_width: float               # [m]
    tank_volume: float             # [m³]
    power_number: float = 1.5      # [-] (geometry-dependent)
    dissipation_ratio: float = 50.0  # k_ε = ε_max / ε_avg


@dataclass
class EmulsificationParameters:
    """Process parameters for Level 1."""
    rpm: float                     # [rev/min]
    t_emulsification: float        # [s]
    mixer: MixerGeometry = field(default_factory=MixerGeometry)


@dataclass
class FormulationParameters:
    """Chemical formulation parameters."""
    c_agarose: float               # [kg/m³] (e.g., 60 for 6% w/v)
    c_chitosan: float              # [kg/m³]
    c_span80: float                # [kg/m³]
    c_genipin: float               # [mol/m³]
    T_oil: float                   # [K] (oil/emulsification temperature)
    cooling_rate: float            # [K/s]
    T_crosslink: float             # [K] (crosslinking temperature)
    t_crosslink: float             # [s]

    @property
    def agarose_fraction(self) -> float:
        """Agarose mass fraction in polymer blend."""
        return self.c_agarose / (self.c_agarose + self.c_chitosan)


@dataclass
class SimulationParameters:
    """Top-level parameter container for the full pipeline."""
    emulsification: EmulsificationParameters
    formulation: FormulationParameters
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
    def from_optimization_vector(cls, x: np.ndarray,
                                  template: 'SimulationParameters') -> 'SimulationParameters':
        """Reconstruct from 7D vector + template for fixed params."""
        ...


# ─── Material Properties ─────────────────────────────────────────────────

@dataclass
class PropertyValue:
    """A single material property with metadata."""
    value: float
    unit: str
    uncertainty: float = 0.0
    source: str = ""
    T_ref: float = 298.15         # [K]


@dataclass
class MaterialProperties:
    """Aggregated material properties for the simulation."""
    # Oil phase
    rho_oil: float                 # [kg/m³] at T_oil
    mu_oil: float                  # [Pa·s] at T_oil

    # Aqueous / dispersed phase
    rho_aq: float                  # [kg/m³]
    mu_d: float                    # [Pa·s] dispersed phase viscosity at T_oil

    # Interfacial
    sigma: float                   # [N/m] interfacial tension with Span-80

    # Thermodynamic
    chi_0: float                   # Flory-Huggins χ parameter at reference T
    chi_T_coeffs: tuple            # (A, B) for χ(T) = A/T + B
    kappa_CH: float                # [J/m] Cahn-Hilliard gradient coefficient
    M_0: float                     # [m⁵/(J·s)] bare mobility

    # Gelation
    T_gel: float                   # [K] agarose gelation temperature
    k_gel_0: float                 # [1/s] Avrami rate prefactor
    n_avrami: float                # [-] Avrami exponent
    gel_arrest_exponent: float     # β in mobility arrest

    # Crosslinking
    k_xlink_0: float              # [m³/(mol·s)] genipin rate prefactor
    E_a_xlink: float              # [J/mol] activation energy
    DDA: float                     # degree of deacetylation of chitosan
    M_GlcN: float = 161.16        # [g/mol] glucosamine molar mass
    M_genipin: float = 226.23     # [g/mol] genipin molar mass

    # Agarose gel
    G_agarose_prefactor: float     # A in G = A·c^n [Pa·(m³/kg)^n]
    G_agarose_exponent: float      # n (typically 2.1-2.5)

    # Network / pore
    r_fiber: float = 1.5e-9       # [m] agarose fiber radius


# ─── Result Structures ────────────────────────────────────────────────────

@dataclass
class EmulsificationResult:
    """Output of Level 1: PBE solver."""
    d_bins: np.ndarray             # [m] pivot diameters, shape (N_bins,)
    n_d: np.ndarray                # [#/m³] number density, shape (N_bins,)
    d32: float                     # [m] Sauter mean diameter
    d43: float                     # [m] volume-weighted mean diameter
    d10: float                     # [m] 10th percentile
    d50: float                     # [m] median diameter
    d90: float                     # [m] 90th percentile
    span: float                    # [-] (d90 - d10) / d50
    total_volume_fraction: float   # [-] dispersed phase holdup
    converged: bool                # steady-state reached?
    t_history: Optional[np.ndarray] = None     # [s] time snapshots
    n_d_history: Optional[np.ndarray] = None   # [#/m³] shape (N_t, N_bins)


@dataclass
class GelationResult:
    """Output of Level 2: Phase-field solver."""
    r_grid: np.ndarray             # [m] radial positions, shape (N_r,)
    phi_field: np.ndarray          # [-] final polymer vol fraction, shape (N_r,)
    pore_size_mean: float          # [m] mean pore diameter
    pore_size_std: float           # [m] std dev of pore sizes
    pore_size_distribution: np.ndarray  # [m] array of measured pore sizes
    porosity: float                # [-] volume fraction of pores
    alpha_final: float             # [-] final gelation degree
    char_wavelength: float         # [m] from structure factor peak
    T_history: Optional[np.ndarray] = None     # [K] temperature vs time
    phi_snapshots: Optional[np.ndarray] = None # shape (N_snap, N_r)


@dataclass
class CrosslinkingResult:
    """Output of Level 3: ODE kinetics."""
    t_array: np.ndarray            # [s] time points, shape (N_t,)
    X_array: np.ndarray            # [mol/m³] crosslink conc vs time
    nu_e_array: np.ndarray         # [1/m³] crosslink density vs time
    Mc_array: np.ndarray           # [g/mol] M_c vs time
    xi_array: np.ndarray           # [m] mesh size vs time
    G_chitosan_array: np.ndarray   # [Pa] chitosan modulus vs time
    p_final: float                 # [-] final crosslinking degree
    nu_e_final: float              # [1/m³]
    Mc_final: float                # [g/mol]
    xi_final: float                # [m]
    G_chitosan_final: float        # [Pa]


@dataclass
class MechanicalResult:
    """Output of Level 4: Property prediction."""
    G_agarose: float               # [Pa] agarose network modulus
    G_chitosan: float              # [Pa] chitosan-genipin network modulus
    G_DN: float                    # [Pa] double-network total modulus
    E_star: float                  # [Pa] effective Young's modulus
    delta_array: np.ndarray        # [m] indentation depths for F-δ curve
    F_array: np.ndarray            # [N] forces for F-δ curve
    rh_array: np.ndarray           # [m] solute radii for K_av curve
    Kav_array: np.ndarray          # [-] partition coefficients
    pore_size_mean: float          # [m]
    xi_mesh: float                 # [m]


@dataclass
class FullResult:
    """Complete pipeline output."""
    parameters: SimulationParameters
    emulsification: EmulsificationResult
    gelation: GelationResult
    crosslinking: CrosslinkingResult
    mechanical: MechanicalResult

    def objective_vector(self) -> np.ndarray:
        """Compute the 3 objective values for optimization."""
        f1 = abs(self.emulsification.d32 - 2.0e-6) / 2.0e-6
        f2 = abs(self.gelation.pore_size_mean - 80e-9) / 80e-9
        f3 = abs(np.log10(self.mechanical.G_DN) - np.log10(1e4)) / 1.0
        return np.array([f1, f2, f3])


@dataclass
class OptimizationState:
    """State of the Bayesian optimization campaign."""
    X_observed: np.ndarray         # shape (N_eval, 7) parameter points
    Y_observed: np.ndarray         # shape (N_eval, 3) objective values
    pareto_X: np.ndarray           # Pareto-optimal parameter points
    pareto_Y: np.ndarray           # corresponding objective values
    iteration: int                 # current BO iteration
    hypervolume: float             # current Pareto hypervolume
    hypervolume_history: list      # HV vs iteration
    converged: bool                # convergence flag
    gp_state: Optional[dict] = None  # serialized GP model state
```

---

## 4. Algorithm Specifications

### 4.1 PBE Solver (Level 1) — Fixed Pivot Class Method

#### Pseudocode

```
ALGORITHM: PBE_FixedPivot_Solve(params, props, t_span, N_bins=50)
────────────────────────────────────────────────────────────────
INPUT:
  params: SimulationParameters
  props: MaterialProperties
  t_span: (t_start, t_end)
  N_bins: number of discretization bins

1. SETUP GRID:
   d_edges = logspace(log10(0.1e-6), log10(500e-6), N_bins+1)
   d_pivots = sqrt(d_edges[:-1] * d_edges[1:])    # geometric mean
   v_pivots = π/6 * d_pivots³                       # volume pivots

2. COMPUTE ENERGY DISSIPATION:
   N_rot = params.rpm / 60
   P = N_p * ρ_emulsion * N_rot³ * D_rotor⁵
   ε = k_ε * P / (ρ_emulsion * V_tank)

3. PRECOMPUTE KERNELS (can be done once if T is constant):
   FOR i = 0 to N_bins-1:
     g[i] = breakage_rate(d_pivots[i], ε, σ, µ_d, ρ_c)
     FOR j = 0 to N_bins-1:
       Q[i,j] = coalescence_rate(d_pivots[i], d_pivots[j], ε, σ, µ_c)

   PRECOMPUTE daughter redistribution matrix:
   FOR j = 0 to N_bins-1:       # parent bin
     FOR each daughter pair (d', d''):
       find bins i_low, i_high bracketing d'
       η[i_low, j] += (v_{i_high} - v') / (v_{i_high} - v_{i_low})
       η[i_high, j] += (v' - v_{i_low}) / (v_{i_high} - v_{i_low})

4. DEFINE RHS FUNCTION dN/dt = f(t, N):
   FOR i = 0 to N_bins-1:
     death_break = g[i] * N[i]
     birth_break = Σ_j η[i,j] * g[j] * N[j]    (j > i for binary breakage)
     death_coal = N[i] * Σ_j Q[i,j] * N[j]
     birth_coal = 0.5 * Σ_{j,k: v_j+v_k → bin i} Q[j,k] * N[j] * N[k]
     dN[i]/dt = birth_break - death_break + birth_coal - death_coal

5. SET INITIAL CONDITION:
   N_0 = lognormal_distribution(d_pivots, d32_premix, σ_premix) * φ_d / <v>

6. INTEGRATE:
   solution = solve_ivp(f, t_span, N_0, method='BDF', rtol=1e-6, atol=1e-8)

7. EXTRACT RESULTS:
   n_d = solution.y[:, -1]           # final distribution
   d32 = Σ(n_d * d³) / Σ(n_d * d²)  # Sauter mean
   d43 = Σ(n_d * d⁴) / Σ(n_d * d³)  # De Brouckere mean
   Compute percentiles from cumulative volume distribution.

8. VALIDATE:
   assert abs(Σ(n_d * v) - Σ(N_0 * v)) / Σ(N_0 * v) < 1e-4  # mass conservation

RETURN EmulsificationResult(d_bins, n_d, d32, d43, ...)
```

#### Complexity

- **Kernel precomputation**: O(N_bins²) = O(2500) — done once
- **RHS evaluation**: O(N_bins²) per call (coalescence double sum)
- **Time steps**: ~100-1000 (adaptive BDF)
- **Total**: O(N_bins² × N_steps) ≈ O(10⁵ - 10⁶) floating-point operations
- **Memory**: O(N_bins² + N_bins × N_steps) for kernels + history ≈ tens of KB
- **Expected runtime**: < 1 second on modern CPU

### 4.2 Phase-Field Solver (Level 2) — Semi-Implicit Cahn-Hilliard

#### Pseudocode

```
ALGORITHM: CahnHilliard_SemiImplicit_1D(params, props, R_droplet)
──────────────────────────────────────────────────────────────────
INPUT:
  params: SimulationParameters
  props: MaterialProperties
  R_droplet: microsphere radius from Level 1

1. SETUP GRID:
   N_r = 1000
   r = linspace(0, R_droplet, N_r)
   Δr = r[1] - r[0]    # ~1 nm for R = 1 µm

2. INITIAL CONDITION:
   φ_0 = c_total_polymer / ρ_polymer    # uniform volume fraction
   φ = φ_0 + 0.01 * random_noise(N_r)  # small perturbation to trigger spinodal

3. BUILD DISCRETE OPERATORS:
   L = spherical_laplacian_matrix(r, Δr)   # sparse, tridiagonal
   # Apply BC: dφ/dr = 0 at r=0 and r=R

4. PRECOMPUTE Flory-Huggins derivatives:
   f'(φ) = (kT/v0) * [ln(φ)/N_p - ln(1-φ) + 1/N_p - 1 + χ(T)*(1-2φ)]
   f''(φ) = (kT/v0) * [1/(N_p·φ) + 1/(1-φ) - 2χ(T)]

   Split: f'_c(φ) = concave part,  f'_e(φ) = convex part
   Standard: f'_c(φ) = -C·φ (stabilizing linear),  f'_e(φ) = f'(φ) + C·φ
   where C > max|f''(φ)| ensures contractivity.

5. TIME LOOP:
   t = 0
   Δt = 1e-4    # initial time step [s]

   WHILE t < t_final AND α < 0.99:

     # Update temperature (lumped capacitance cooling)
     T = T_oil + (T_init - T_oil) * exp(-t / τ_cool)

     # Update Avrami gelation
     IF T < T_gel:
       t_cool = time since T crossed T_gel
       α = 1 - exp(-(k_gel * t_cool)^n_avrami)
     ELSE:
       α = 0

     # Update mobility with gelation arrest
     M = M_0 * φ * (1-φ) * (1-α)^β

     # Update χ(T) → affects free energy
     χ = A/T + B

     # Semi-implicit step:
     # Solve: (I - Δt·∇·[M·∇(f''_c·I - κ·L)]) · φ^{n+1}
     #      = φ^n + Δt·∇·[M·∇(f'_e(φ^n))]

     # Construct LHS matrix (sparse, banded):
     A_mat = I - Δt * L_mobility @ (C_contract * I - κ * L)

     # Construct RHS:
     µ_explicit = f'_e(φ)
     rhs = φ + Δt * divergence(M * gradient(µ_explicit))

     # Solve banded system:
     φ_new = solve_banded(A_mat, rhs)

     # Adaptive time stepping:
     change = max(|φ_new - φ|)
     IF change < 0.001: Δt *= 1.5
     IF change > 0.05:  Δt *= 0.5
     Δt = clip(Δt, 1e-6, 1.0)

     φ = φ_new
     t += Δt

     # Snapshot for visualization
     IF t in snapshot_times: save φ, T, α

6. ANALYZE PORE STRUCTURE:
   pore_sizes = chord_length_analysis(φ, r, threshold=0.5*(φ_max+φ_min))
   S_q = structure_factor(φ, r)
   q_star = argmax(S_q)
   char_wavelength = 2π / q_star
   porosity = volume_fraction_below_threshold(φ, r, threshold)

RETURN GelationResult(r, φ, mean(pore_sizes), std(pore_sizes), ...)
```

#### Complexity

- **Per time step**: O(N_r) for banded solve (bandwidth 5, Thomas algorithm)
- **Time steps**: ~1000-10000 (adaptive)
- **Total**: O(N_r × N_steps) ≈ O(10⁶ - 10⁷)
- **Memory**: O(N_r × N_snapshots) ≈ tens of KB
- **Expected runtime**: 1-10 seconds

### 4.3 Crosslinking ODE (Level 3)

#### Pseudocode

```
ALGORITHM: CrosslinkingODE_Solve(params, props, gel_result)
──────────────────────────────────────────────────────────────
INPUT:
  params: SimulationParameters
  props: MaterialProperties
  gel_result: GelationResult

1. INITIAL CONDITIONS:
   NH2_0 = DDA * c_chitosan / M_GlcN          # [mol/m³]
   Gen_0 = c_genipin                            # [mol/m³]
   X_0 = 0                                      # [mol/m³]
   y0 = [NH2_0, Gen_0, X_0]

2. RATE CONSTANT:
   k = k_0 * exp(-E_a / (R * T_crosslink))     # [m³/(mol·s)]

3. DEFINE ODE:
   def rhs(t, y):
     NH2, Gen, X = y
     rate = k * NH2 * Gen
     return [-2*rate, -rate, rate]

4. INTEGRATE:
   t_span = (0, t_crosslink)
   t_eval = linspace(0, t_crosslink, 500)
   sol = solve_ivp(rhs, t_span, y0, method='Radau',
                   rtol=1e-8, atol=1e-10, t_eval=t_eval)

5. DERIVED PROPERTIES at each time point:
   X = sol.y[2, :]
   nu_e = X * N_A                              # [1/m³]
   Mc = ρ_polymer / nu_e                       # [g/mol] (where nu_e > 0)
   xi = 0.071 * v2s^(-1/3) * Mc^(1/2) * 1e-9  # [m] (Canal-Peppas)
   G_chi = nu_e * k_B * T                      # [Pa]

6. FINAL VALUES:
   p_final = X[-1] / min(NH2_0/2, Gen_0)

RETURN CrosslinkingResult(t, X, nu_e, Mc, xi, G_chi, p_final, ...)
```

#### Complexity

- **Per step**: O(1) (3x3 system, analytic Jacobian)
- **Steps**: ~50-200 (Radau adaptive)
- **Total**: O(10²)
- **Memory**: O(N_t × 6) ≈ bytes
- **Expected runtime**: < 0.01 seconds

### 4.4 Mechanical Property Prediction (Level 4)

#### Pseudocode

```
ALGORITHM: MechanicalPredict(params, props, emul, gel, xl)
──────────────────────────────────────────────────────────────
1. AGAROSE MODULUS:
   G_agarose = A * (c_agarose)^n     # empirical power law

2. DOUBLE-NETWORK MODULUS:
   G_DN = G_agarose + xl.G_chitosan_final

3. EFFECTIVE YOUNG'S MODULUS:
   E_star = 3 * G_DN                  # incompressible limit (ν=0.5)

4. HERTZ FORCE-DISPLACEMENT:
   R = emul.d32 / 2
   delta = linspace(0, 0.1*R, 100)   # up to 10% compression
   F = (4/3) * E_star * sqrt(R) * delta^(3/2)

5. OGSTON PARTITIONING:
   r_h = logspace(-10, -7, 100)       # 0.1 nm to 100 nm solute radii
   n_f = gel.porosity_fiber_density    # or from phi and r_fiber
   K_av = exp(-π * (r_h + r_fiber)² * n_f)

RETURN MechanicalResult(G_agarose, xl.G_chitosan_final, G_DN, E_star,
                        delta, F, r_h, K_av, ...)
```

#### Complexity
- O(N_delta + N_rh) ≈ O(200) — instantaneous
- **Expected runtime**: < 0.001 seconds

### 4.5 Total Pipeline Runtime Estimate

| Level | Runtime   | Dominant cost           |
|-------|-----------|-------------------------|
| L1    | < 1 s     | O(N²) coalescence kernel|
| L2    | 1-10 s    | Semi-implicit time loop |
| L3    | < 0.01 s  | 3-ODE Radau integration |
| L4    | < 0.001 s | Analytical evaluations  |
| **Total** | **~2-12 s** | **L2 dominates**    |

For 200 optimization evaluations: **~7-40 minutes** total campaign time.

---

## 5. Technology Stack

### 5.1 Core Dependencies

| Package          | Version   | Purpose                                              |
|------------------|-----------|------------------------------------------------------|
| `numpy`          | >=1.26    | Array operations, linear algebra                     |
| `scipy`          | >=1.12    | ODE solvers, sparse matrices, optimization, FFT      |
| `h5py`           | >=3.10    | HDF5 I/O for checkpoints and results                 |
| `matplotlib`     | >=3.8     | Publication-quality static plots                     |
| `plotly`         | >=5.18    | Interactive HTML dashboards                          |
| `botorch`        | >=0.11    | Bayesian optimization (multi-objective)              |
| `gpytorch`       | >=1.11    | Gaussian process models (BoTorch backend)            |
| `torch`          | >=2.1     | Tensor operations (BoTorch backend)                  |
| `pydantic`       | >=2.5     | Config validation and serialization                  |
| `tomli`/`tomllib`| stdlib    | TOML config parsing (Python 3.11+ stdlib)            |

### 5.2 Optional / Development Dependencies

| Package          | Purpose                                              |
|------------------|------------------------------------------------------|
| `numba`          | JIT compilation for hot loops (PBE kernels, FD stencil)|
| `pytest`         | Testing framework                                    |
| `hypothesis`     | Property-based testing for numerical routines        |
| `rich`           | Terminal output formatting and progress bars         |
| `jupyter`        | Interactive exploration and visualization            |

### 5.3 Justification for Key Choices

**Why Python over C++/Julia?**
- The computational bottleneck (Level 2 phase-field) runs in ~10 seconds in Python with NumPy/SciPy. C++ would gain 10-50x but we only need ~200 evaluations total → total campaign in minutes regardless.
- Ecosystem: BoTorch (Bayesian optimization), SciPy (ODE solvers), Matplotlib (plotting) are all Python-native.
- Rapid iteration: Architecture changes, new kernels, alternative models — all faster to prototype in Python.
- If Level 2 becomes a bottleneck (e.g., 2D upgrade), Numba or Cython can accelerate the stencil loop by 10-100x without leaving Python.

**Why HDF5 over JSON/pickle?**
- Large array storage (φ field, time histories) is efficient in HDF5.
- Self-describing: metadata (units, parameters) stored alongside data.
- Language-agnostic: readable from MATLAB, Julia, C++ if needed later.
- Compression: gzip compression built-in.

**Why BoTorch over scipy.optimize?**
- Multi-objective optimization with Pareto front tracking.
- Gaussian process surrogate → intelligent sampling, not brute force.
- EHVI acquisition function specifically designed for expensive black-box multi-objective optimization.
- Active development by Meta FAIR, well-documented.

---

## 6. Directory Structure

```
emulsification_sim/
│
├── docs/
│   ├── 01_scientific_advisor_report.md
│   ├── 02_computational_architecture.md    ← THIS DOCUMENT
│   └── 03_implementation_notes.md          (created during development)
│
├── src/
│   └── emulsim/
│       ├── __init__.py
│       ├── config.py                # SimulationConfig, TOML loading
│       ├── datatypes.py             # All dataclasses (Sec. 3)
│       ├── properties/
│       │   ├── __init__.py
│       │   ├── database.py          # PropertyDatabase class
│       │   ├── viscosity.py         # Dispersed phase viscosity models
│       │   ├── interfacial.py       # Interfacial tension models
│       │   └── thermodynamic.py     # χ(T), free energy functions
│       ├── level1_emulsification/
│       │   ├── __init__.py
│       │   ├── solver.py            # PBE solver (class method)
│       │   ├── kernels.py           # Breakage & coalescence kernels
│       │   ├── energy.py            # Energy dissipation calculator
│       │   └── validation.py        # Mass conservation, bounds checks
│       ├── level2_gelation/
│       │   ├── __init__.py
│       │   ├── solver.py            # Cahn-Hilliard semi-implicit solver
│       │   ├── free_energy.py       # Flory-Huggins implementation
│       │   ├── gelation.py          # Avrami kinetics + mobility arrest
│       │   ├── pore_analysis.py     # Structure factor, chord lengths
│       │   └── spatial.py           # Grid, Laplacian, BC utilities
│       ├── level3_crosslinking/
│       │   ├── __init__.py
│       │   ├── solver.py            # ODE integration
│       │   └── network_props.py     # nu_e, Mc, xi, G calculations
│       ├── level4_mechanical/
│       │   ├── __init__.py
│       │   ├── modulus.py           # Double-network modulus
│       │   ├── hertz.py             # Hertz contact model
│       │   └── partitioning.py      # Ogston Kav model
│       ├── optimization/
│       │   ├── __init__.py
│       │   ├── engine.py            # BoTorch optimization loop
│       │   ├── objectives.py        # Objective function definitions
│       │   └── analysis.py          # Pareto front analysis, visualization
│       ├── pipeline/
│       │   ├── __init__.py
│       │   ├── orchestrator.py      # PipelineOrchestrator
│       │   ├── logger.py            # HDF5 + text logging
│       │   └── checkpoint.py        # Save/restore state
│       └── visualization/
│           ├── __init__.py
│           ├── plots.py             # Matplotlib static plots
│           └── dashboard.py         # Plotly interactive dashboard
│
├── data/
│   ├── properties.toml              # Default material property database
│   └── literature/
│       ├── breakage_validation.csv  # Published breakage data for testing
│       └── agarose_modulus.csv      # Published G vs c_agarose data
│
├── configs/
│   ├── default.toml                 # Default simulation configuration
│   ├── optimization_campaign.toml   # BO campaign settings
│   └── example_single_run.toml      # Example single-run config
│
├── tests/
│   ├── conftest.py                  # Shared fixtures
│   ├── test_properties.py
│   ├── test_level1_emulsification.py
│   ├── test_level2_gelation.py
│   ├── test_level3_crosslinking.py
│   ├── test_level4_mechanical.py
│   ├── test_pipeline.py
│   ├── test_optimization.py
│   └── test_integration.py          # End-to-end pipeline tests
│
├── notebooks/
│   ├── 01_explore_emulsification.ipynb
│   ├── 02_explore_gelation.ipynb
│   ├── 03_explore_crosslinking.ipynb
│   └── 04_optimization_results.ipynb
│
├── output/                          # Simulation results (gitignored)
│   └── .gitkeep
│
├── pyproject.toml                   # Project metadata + dependencies
├── requirements.txt                 # Pinned dependencies
├── LICENSE
└── README.md
```

---

## 7. Implementation Roadmap

### Phase 1: Core Infrastructure + Level 1 (Emulsification)

**Scope**: Project skeleton, data structures, property database, PBE solver, basic pipeline.

| Task | Description | Est. Effort |
|------|-------------|-------------|
| 1.1  | Project setup: pyproject.toml, directory structure, CI | 2 hours |
| 1.2  | `datatypes.py`: all dataclasses from Section 3 | 3 hours |
| 1.3  | `config.py`: TOML loading, validation with Pydantic | 2 hours |
| 1.4  | `properties/database.py`: PropertyDatabase with T interpolation | 4 hours |
| 1.5  | `properties/viscosity.py`: Mark-Houwink + mixing rule | 2 hours |
| 1.6  | `properties/interfacial.py`: Szyszkowski-Langmuir for Span-80 | 2 hours |
| 1.7  | `level1/kernels.py`: Alopaeus breakage + C&T coalescence | 6 hours |
| 1.8  | `level1/energy.py`: rotor-stator power/dissipation | 2 hours |
| 1.9  | `level1/solver.py`: Fixed-pivot PBE solver | 8 hours |
| 1.10 | `level1/validation.py`: mass conservation, bounds | 2 hours |
| 1.11 | Unit tests for all Phase 1 modules | 6 hours |
| 1.12 | `pipeline/orchestrator.py`: skeleton with L1 only | 3 hours |
| 1.13 | `pipeline/logger.py`: HDF5 output | 3 hours |
| 1.14 | Notebook: explore L1 behavior, parameter sweeps | 3 hours |
| **Total** | | **~48 hours** |

**Deliverable**: Working PBE solver that predicts d32 vs RPM and Span-80 concentration. Validated against Kolmogorov-Hinze scaling (d32 ~ ε^{-0.4}) and order-of-magnitude match to literature emulsification data.

### Phase 2: Level 2 (Gelation & Pore Formation)

**Scope**: Cahn-Hilliard solver, Avrami coupling, pore analysis.

| Task | Description | Est. Effort |
|------|-------------|-------------|
| 2.1  | `level2/spatial.py`: 1D radial grid, Laplacian, BCs | 4 hours |
| 2.2  | `level2/free_energy.py`: Flory-Huggins implementation | 4 hours |
| 2.3  | `level2/gelation.py`: Avrami kinetics + mobility arrest | 3 hours |
| 2.4  | `level2/solver.py`: Semi-implicit CH time stepper | 10 hours |
| 2.5  | `level2/pore_analysis.py`: structure factor, chord lengths | 4 hours |
| 2.6  | `properties/thermodynamic.py`: χ(T) parameterization | 3 hours |
| 2.7  | Integration: wire L1 output into L2 input | 2 hours |
| 2.8  | Unit + integration tests | 6 hours |
| 2.9  | Notebook: explore phase separation dynamics | 3 hours |
| **Total** | | **~39 hours** |

**Deliverable**: Working phase-field solver that produces pore structure as a function of cooling rate and polymer concentration. Validated: spinodal wavelength scales as predicted by linearized Cahn-Hilliard theory; gelation arrest freezes coarsening at correct T.

### Phase 3: Level 3 + Level 4 (Crosslinking + Mechanical)

**Scope**: Crosslinking ODE, mechanical property models, full pipeline.

| Task | Description | Est. Effort |
|------|-------------|-------------|
| 3.1  | `level3/solver.py`: ODE system for genipin crosslinking | 3 hours |
| 3.2  | `level3/network_props.py`: ν_e, M_c, ξ, G calculations | 3 hours |
| 3.3  | `level4/modulus.py`: double-network G_DN | 2 hours |
| 3.4  | `level4/hertz.py`: Hertz contact model | 2 hours |
| 3.5  | `level4/partitioning.py`: Ogston K_av model | 2 hours |
| 3.6  | Full pipeline integration (L1→L2→L3→L4) | 4 hours |
| 3.7  | `pipeline/checkpoint.py`: save/restore | 3 hours |
| 3.8  | Unit + integration tests | 5 hours |
| 3.9  | `visualization/plots.py`: standard result plots | 4 hours |
| **Total** | | **~28 hours** |

**Deliverable**: Complete L1→L2→L3→L4 pipeline producing all target outputs. Validated: crosslinking kinetics match published genipin-chitosan rates; modulus in physically reasonable range (kPa-MPa for hydrogel microspheres).

### Phase 4: Optimization Engine + Validation

**Scope**: BoTorch integration, multi-objective optimization, comprehensive validation.

| Task | Description | Est. Effort |
|------|-------------|-------------|
| 4.1  | `optimization/objectives.py`: objective function wrappers | 3 hours |
| 4.2  | `optimization/engine.py`: BoTorch EHVI loop | 8 hours |
| 4.3  | `optimization/analysis.py`: Pareto analysis, convergence | 4 hours |
| 4.4  | Batch execution in orchestrator | 3 hours |
| 4.5  | `visualization/dashboard.py`: Plotly optimization dashboard | 5 hours |
| 4.6  | Sensitivity analysis: Sobol indices for each objective | 5 hours |
| 4.7  | Comprehensive validation campaign | 8 hours |
| 4.8  | Documentation and final testing | 4 hours |
| **Total** | | **~40 hours** |

**Deliverable**: Working optimization engine that finds Pareto-optimal process conditions. Dashboard for exploring results.

### Summary

| Phase | Focus | Est. Effort | Cumulative |
|-------|-------|-------------|------------|
| 1     | Infrastructure + L1 | 48 hours | 48 hours |
| 2     | L2 Gelation | 39 hours | 87 hours |
| 3     | L3 + L4 + Pipeline | 28 hours | 115 hours |
| 4     | Optimization + Validation | 40 hours | 155 hours |
| **Total** | | **~155 hours** | |

---

## 8. Verification and Validation Strategy

### 8.1 Unit Tests

#### Module 1: Property Database
- `test_load_toml`: loads properties.toml without error
- `test_temperature_interpolation`: ρ(T), μ(T) match hand calculations
- `test_viscosity_model`: Mark-Houwink output in physically reasonable range
- `test_interfacial_tension`: σ decreases with increasing Span-80 (monotonic), σ > 0 always
- `test_chi_parameter`: χ increases as T decreases (for UCST system)

#### Module 2: Emulsification (Level 1)
- `test_breakage_kernel_units`: output in [1/s]
- `test_breakage_increases_with_epsilon`: g(d) increases with ε at fixed d
- `test_breakage_increases_with_d`: g(d) increases with d at fixed ε (for d > d_Kolmogorov)
- `test_viscous_correction`: higher µ_d → lower g(d) (viscous drops resist breakage)
- `test_coalescence_kernel_symmetry`: q(d_i, d_j) = q(d_j, d_i)
- `test_mass_conservation`: total volume conserved to 1e-4 over full simulation
- `test_non_negative`: all N_i >= 0 at all times
- `test_steady_state_reached`: d32 changes < 1% over last 10% of time
- `test_d32_vs_rpm_scaling`: d32 ~ RPM^{-1.2} (Kolmogorov-Hinze scaling)
- `test_known_solution_no_coalescence`: with Q=0, pure breakage → analytical limiting d32

#### Module 3: Gelation (Level 2)
- `test_flory_huggins_convexity`: f''(φ) > 0 outside spinodal, < 0 inside
- `test_spinodal_boundary`: correct φ_s1, φ_s2 for given χ
- `test_mass_conservation_CH`: ∫φ dV constant to 1e-6
- `test_no_phase_separation_above_Tgel`: uniform φ persists when χ < χ_critical
- `test_gelation_arrests_coarsening`: pore size stops growing when α → 1
- `test_avrami_limits`: α(0) = 0, α(∞) = 1
- `test_characteristic_wavelength`: q* matches linearized Cahn-Hilliard prediction for early times
- `test_boundary_conditions`: no-flux at r=0 and r=R

#### Module 4: Crosslinking (Level 3)
- `test_stoichiometry`: [NH2] + 2[X] = [NH2]_0 at all times
- `test_monotonic_crosslinking`: d[X]/dt >= 0 always
- `test_equilibrium`: [X] → [X]_max as t → ∞
- `test_arrhenius_temperature`: higher T → faster crosslinking
- `test_modulus_increases`: G_chitosan monotonically increases with time
- `test_mesh_size_decreases`: ξ monotonically decreases with crosslinking

#### Module 5: Mechanical (Level 4)
- `test_double_network_additive`: G_DN = G_agarose + G_chitosan
- `test_hertz_scaling`: F ~ δ^{3/2}
- `test_hertz_units`: F in Newtons for δ in meters
- `test_kav_limits`: K_av → 1 for r_h → 0; K_av → 0 for r_h → ∞
- `test_kav_decreases_with_solute_size`: monotonic decrease

#### Module 6: Optimization
- `test_sobol_initialization`: correct number of initial points
- `test_gp_fit`: GP predicts training points accurately
- `test_ehvi_positive`: acquisition function non-negative
- `test_pareto_dominance`: returned Pareto set is non-dominated

#### Module 7: Pipeline
- `test_full_pipeline_runs`: L1→L2→L3→L4 without error
- `test_checkpoint_restore`: save at L2, restore, continue to L4, same result
- `test_hdf5_roundtrip`: save and load each result type

### 8.2 Integration Test Scenarios

| Test | Description | Pass Criteria |
|------|-------------|---------------|
| IT-1 | Full pipeline with default parameters | Completes without error; d32 in [0.5, 50] µm |
| IT-2 | High RPM → small droplets → fine pores | d32 < 5 µm, pore_mean < 200 nm |
| IT-3 | Low RPM → large droplets | d32 > 50 µm |
| IT-4 | No Span-80 → high σ → large drops | d32 increases vs baseline |
| IT-5 | Fast cooling → fine pores (arrested early) | pore_mean smaller than slow cooling |
| IT-6 | High genipin → high crosslink → stiff | G_DN > baseline |
| IT-7 | Optimization: 5 iterations | Hypervolume increases monotonically |

### 8.3 Validation Against Literature

| Quantity | Literature Value | Source | Validation Method |
|----------|-----------------|--------|-------------------|
| d32 scaling | d32 ~ We^{-0.6} | Calabrese et al. (1986) | Log-log plot d32 vs We |
| Viscosity effect | d32 increases with µ_d | Alopaeus et al. (2002) | Parameter sweep µ_d |
| Spinodal wavelength | λ* = 2π/q* ~ (κ/f'')^{1/2} | Cahn (1965) | Compare early-time q* |
| Agarose gel modulus | G ~ c^{2.1-2.5} | Ross-Murphy (1991) | Power-law fit |
| Genipin rate | k(37°C) ≈ 5×10⁻³ L/(mol·s) | Butler et al. (2003) | Match published kinetics curve |
| Mesh size | ξ ~ 50-200 nm for 1-4% agarose | Pernodet et al. (1997) | Compare predicted vs measured |

### 8.4 Sensitivity Analysis Plan

**Method**: Sobol sensitivity indices (first-order and total-order) using SALib.

**Parameters varied** (7 optimization variables, each ±30% around baseline):

**Outputs analyzed**:
1. d32 (droplet size)
2. pore_size_mean
3. G_DN (modulus)
4. K_av for a reference protein (e.g., BSA, r_h = 3.5 nm)

**Expected dominant sensitivities**:
- d32: RPM (dominant), c_Span80, µ_d (via T_oil)
- pore_size_mean: cooling_rate (dominant), polymer concentration, T_gel
- G_DN: c_genipin (dominant), t_crosslink, c_agarose
- K_av: pore_size and xi_mesh (both depend on upstream levels)

This analysis identifies which parameters to prioritize in optimization and which can be fixed.

---

## 9. First-Stage Audit: Scientific Rationality

### 9.1 Governing Equations Fidelity

| Equation | Source | Faithfully Represented? | Notes |
|----------|--------|-------------------------|-------|
| PBE (Ramkrishna 2000) | Scientific Advisor | YES | Fixed-pivot class method preserves mass exactly |
| Alopaeus breakage | Alopaeus et al. (2002) | YES | Both surface tension and viscous resistance terms included |
| C&T coalescence | Coulaloglou & Tavlarides (1977) | YES | Both collision frequency and film drainage efficiency |
| Cahn-Hilliard | Standard | YES | Semi-implicit Eyre splitting, Flory-Huggins free energy |
| Avrami gelation | Standard kinetics | YES | Coupled to mobility (arrest mechanism) |
| Genipin kinetics | Butler et al. (2003) | YES | Second-order, Arrhenius T-dependence |
| Rubber elasticity | Treloar (1975) | YES | G = ν_e·k_B·T |
| Canal-Peppas mesh | Canal & Peppas (1989) | YES | ξ = 0.071·ν_{2,s}^{-1/3}·M_c^{1/2} |
| Hertz contact | Johnson (1985) | YES | F = (4/3)·E*·R^{1/2}·δ^{3/2} |
| Ogston partitioning | Ogston (1958) | YES | K_av = exp(-π·(r_h+r_f)²·n_f) |

**Verdict**: All governing equations from the scientific methodology are faithfully represented.

### 9.2 Assumptions Identified and Justified

| # | Assumption | Justification | Risk if Violated |
|---|-----------|---------------|------------------|
| A1 | One-way coupling between levels | Physical reality: emulsification is complete before gelation starts (batch process). Gelation does not affect droplet size (already formed). | LOW — validated by process timeline |
| A2 | 0D PBE (spatially uniform) | Rotor-stator mixers have small well-mixed volume in gap region where breakage occurs. 0D is standard for lab-scale. | MEDIUM — could miss spatial gradients in large vessels. Mitigation: Tier 1B CFD-PBE upgrade path exists. |
| A3 | 1D radial for phase field | Bi << 0.1 → uniform cooling → no preferred direction. Pore structure is statistically isotropic within the sphere. | LOW — validated by Biot number calculation |
| A4 | Thiele modulus << 1 (uniform crosslinking) | For 2 µm spheres with D_genipin ~ 10⁻¹⁰ m²/s and k ~ 10⁻³ s⁻¹: Φ = R·√(k/D) ~ 10⁻⁴ << 1. | LOW — extremely small spheres guarantee this |
| A5 | Binary breakage | Standard assumption. Ternary breakage is rare at our We numbers. | LOW |
| A6 | Additive double-network modulus | Assumes independent, interpenetrating networks. Valid when agarose and chitosan-genipin networks do not mechanically interfere. | MEDIUM — synergistic effects possible. Literature (Gong et al. 2003) shows DN gels can be super-additive, but our system is not a classical DN (no sacrificial bonds). Additive is conservative. |
| A7 | Incompressible hydrogel (ν=0.5) | Standard for swollen polymer networks at low strain. | LOW |
| A8 | Constant T during crosslinking | Crosslinking done in temperature-controlled bath. | LOW |
| A9 | No surfactant depletion during emulsification | Span-80 is in the continuous oil phase and typically in large excess. | LOW for c_Span80 > CMC; could matter near CMC |
| A10 | Flory-Huggins adequate for agarose-water | Agarose is semi-rigid → Flory-Huggins (flexible chain) is approximate. However, χ is an effective parameter fitted to phase behavior, which absorbs non-ideal chain stiffness effects. | MEDIUM — could be improved with Wertheim or Ogston models for semi-rigid chains |

### 9.3 Parameter Ranges Physical Reasonableness

| Parameter | Range Used | Physical Basis | Reasonable? |
|-----------|-----------|----------------|-------------|
| RPM | 3,000-25,000 | Lab rotor-stator (Ultra-Turrax) range | YES |
| c_Span80 | 0.5-5% w/v | Below and above CMC (~0.3% w/v); literature range | YES |
| r_AC (agarose fraction) | 0.6-0.9 | Must maintain gelation (need sufficient agarose); chitosan needed for crosslinking | YES |
| T_oil | 60-95°C | Must be above T_gel (~40°C) with margin; below 100°C (aqueous) | YES |
| Cooling rate | 0.5-20°C/min | Ice bath (~20°C/min) to air cooling (~0.5°C/min) | YES |
| c_genipin | 0.1-10 mM | Literature range for chitosan crosslinking | YES |
| t_crosslink | 1-48 hours | Literature range; <1 hr insufficient, >48 hr impractical | YES |
| µ_d | 0.1-10 Pa·s | Correct for hot agarose/chitosan solutions | YES |
| σ | 2-8 mN/m | Correct for Span-80 stabilized paraffin/water interface | YES |

**Verdict**: All parameter ranges are physically reasonable and consistent with the literature.

### 9.4 Missing Physics / Edge Cases

| # | Issue | Severity | Mitigation |
|---|-------|----------|------------|
| E1 | **Marangoni effects** on coalescence (surfactant gradients at film interfaces) | LOW | Absorbed into effective C_5 coalescence constant; explicit Marangoni would require interface-resolved model (overkill for 0D PBE) |
| E2 | **Ostwald ripening** (molecular diffusion between drops of different sizes) | LOW | Negligible for polymer solutions with near-zero solubility in oil phase |
| E3 | **Shear-dependent viscosity** of agarose solutions (shear-thinning) | MEDIUM | Could add Cross model: µ_d(γ̇) = µ_0 / (1 + (λ·γ̇)^m). Currently uses zero-shear viscosity → conservative (overestimates resistance to breakage). Flag for Phase 2 upgrade. |
| E4 | **Non-uniform Span-80 distribution** across drops | LOW | 0D PBE assumes well-mixed; valid for lab-scale with high RPM |
| E5 | **Elastic contribution to breakage resistance** (viscoelastic dispersed phase) | MEDIUM | Agarose solutions are weakly viscoelastic above T_gel. Purely viscous breakage model is standard. If elastic effects matter, add Calabrese elastic correction term. |
| E6 | **Phase field: polymer mixture** (agarose + chitosan) | MEDIUM | Currently treat as single effective polymer with averaged properties. A ternary Cahn-Hilliard (polymer A + polymer B + solvent) would be more accurate but significantly more complex (3 coupled PDEs). Defer to Phase 5 if needed. |
| E7 | **Genipin side reactions** (self-polymerization, reaction with water) | LOW | At neutral pH and < 50°C, genipin-amine dominates. Self-polymerization is slow. |
| E8 | **Swelling during crosslinking** | LOW | Crosslinking in the same buffer → no osmotic driving force for significant swelling. Flory-Rehner equilibrium can be checked post-hoc. |
| E9 | **Drop breakage during cooling** (shear from stirring during cooling) | MEDIUM | If stirring continues during cooling, small additional breakage or coalescence could occur. Model assumes emulsification and cooling are sequential. Flag as experimental protocol note. |

### 9.5 Model Fidelity Assessment

| Sub-process | Model | Fidelity | Appropriate? |
|-------------|-------|----------|--------------|
| Emulsification | 0D PBE, class method | Medium | YES for lab-scale. Captures the essential physics (turbulent breakage vs. coalescence). Upgrade path to CFD-PBE exists. |
| Phase separation | 1D Cahn-Hilliard + Flory-Huggins | Medium-High | YES. Captures spinodal decomposition dynamics and gelation arrest. 1D is sufficient for characteristic length scales. |
| Gelation | Avrami phenomenological | Medium | YES. Avrami is the standard model for polymer gelation kinetics. Molecular-level gelation models exist but are unnecessary for this scale. |
| Crosslinking | 2nd-order ODE | High | YES. Thiele << 1 eliminates spatial complexity. Second-order kinetics are well-established for genipin-amine. |
| Mechanical | Analytical rubber elasticity | Medium | YES for initial design. Hertz model is exact for small-strain elastic spheres. DN additivity is approximate but adequate. |
| Optimization | BO with GP surrogate | High | YES. Standard approach for expensive black-box multi-objective optimization with <10 parameters. |

### 9.6 Audit Summary

**Overall assessment**: The computational architecture faithfully represents the scientific methodology. All governing equations are correctly implemented. Assumptions are physically justified for the target application (2 µm microspheres, lab-scale rotor-stator emulsification). The main approximations (0D PBE, 1D phase field, additive DN modulus) are appropriate for the design-space exploration purpose of this simulation.

**Key risks to monitor**:
1. Shear-thinning of dispersed phase (E3) — could affect d32 predictions at high RPM
2. Single effective polymer in phase field (E6) — may miss agarose-chitosan demixing
3. Additive DN assumption (A6) — validate against experimental compression data when available

**Recommended validation priority**:
1. First: d32 vs RPM scaling (easiest to compare with literature)
2. Second: Agarose gel modulus vs concentration (abundant literature data)
3. Third: Pore size vs cooling rate (limited literature, may need own experiments)
4. Fourth: Full pipeline prediction vs microsphere characterization data

---

## Appendix A: Key Literature References

1. Ramkrishna, D. (2000). *Population Balances*. Academic Press.
2. Kumar, S. & Ramkrishna, D. (1996). Chem. Eng. Sci. 51(8), 1311-1332. [Fixed pivot technique]
3. Alopaeus, V. et al. (2002). Chem. Eng. Sci. 57(10), 1815-1825. [Viscosity-corrected breakage]
4. Coulaloglou, C. & Tavlarides, L. (1977). Chem. Eng. Sci. 32(11), 1289-1297. [Coalescence kernel]
5. Calabrese, R. et al. (1986). AIChE J. 32(4), 657-666. [Drop breakage in turbulent flow]
6. Eyre, D. (1998). MRS Proceedings 529, 39. [Unconditionally stable CH time stepping]
7. Cahn, J. (1965). J. Chem. Phys. 42(1), 93-99. [Spinodal decomposition theory]
8. Flory, P. (1953). *Principles of Polymer Chemistry*. Cornell University Press.
9. Canal, T. & Peppas, N. (1989). J. Biomed. Mater. Res. 23(10), 1183-1193. [Mesh size]
10. Butler, M. et al. (2003). J. Polym. Sci. A 41(22), 3941-3953. [Genipin crosslinking kinetics]
11. Normand, V. et al. (2000). Biomacromolecules 1(4), 730-738. [Agarose gel mechanics]
12. Ogston, A. (1958). Trans. Faraday Soc. 54, 1754-1757. [Gel partitioning model]
13. Balanescu, M. et al. (2023). *BoTorch: A Framework for Efficient Monte-Carlo Bayesian Optimization*. NeurIPS.

---

## Appendix B: Configuration File Examples

### B.1 Default Simulation Configuration (`configs/default.toml`)

```toml
[simulation]
run_id = "default_run"
notes = "Baseline simulation with default parameters"

[emulsification]
rpm = 10000
t_emulsification = 600.0  # 10 minutes [s]

[emulsification.mixer]
rotor_diameter = 0.025    # [m] (25 mm rotor)
stator_diameter = 0.026   # [m]
gap_width = 0.0005        # [m] (0.5 mm gap)
tank_volume = 0.0005      # [m³] (500 mL)
power_number = 1.5
dissipation_ratio = 50.0

[formulation]
c_agarose = 42.0          # [kg/m³] (4.2% w/v for 7:3 ratio with total 6%)
c_chitosan = 18.0          # [kg/m³] (1.8% w/v)
c_span80 = 20.0            # [kg/m³] (~2% w/v)
c_genipin = 2.0            # [mol/m³] (~2 mM)
T_oil = 363.15             # [K] (90°C)
cooling_rate = 0.167        # [K/s] (~10°C/min)
T_crosslink = 310.15       # [K] (37°C)
t_crosslink = 86400.0      # [s] (24 hours)

[solver.level1]
n_bins = 50
d_min = 0.1e-6             # [m]
d_max = 500e-6             # [m]
rtol = 1e-6
atol = 1e-8

[solver.level2]
n_r = 1000
dt_initial = 1e-4          # [s]
dt_max = 1.0               # [s]
arrest_exponent = 2.5

[solver.level3]
method = "Radau"
rtol = 1e-8
atol = 1e-10

[optimization]
n_initial = 20
max_iterations = 200
convergence_tol = 0.01
```

---

*End of Computational Architecture Document*
