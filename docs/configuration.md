# Configuration Reference

All simulation parameters are specified in TOML format. The default configuration
is in `configs/default.toml`. Parameters can also be overridden via CLI flags or
the Python API.

## `[simulation]`

| Parameter | Default | Type | Description |
|-----------|---------|------|-------------|
| `run_id` | `""` | string | Identifier for this simulation run |
| `notes` | `""` | string | Free-text notes |

## `[emulsification]`

Process parameters for Level 1 (rotor-stator emulsification).

| Parameter | Default | Unit | Description |
|-----------|---------|------|-------------|
| `rpm` | `10000` | rev/min | Rotor speed |
| `t_emulsification` | `60.0` | s | Emulsification duration (time to reach steady-state DSD) |

### `[emulsification.mixer]`

Rotor-stator mixer geometry.

| Parameter | Default | Unit | Description |
|-----------|---------|------|-------------|
| `rotor_diameter` | `0.025` | m | Rotor outer diameter (25 mm) |
| `stator_diameter` | `0.026` | m | Stator inner diameter |
| `gap_width` | `0.0005` | m | Rotor-stator gap (0.5 mm) |
| `tank_volume` | `0.0005` | m^3 | Vessel volume (500 mL) |
| `power_number` | `1.5` | - | Turbulent power number N_P |
| `dissipation_ratio` | `50.0` | - | Ratio of max to mean energy dissipation rate (epsilon_max / epsilon_avg) |

## `[formulation]`

Chemical formulation and process temperatures.

| Parameter | Default | Unit | Description |
|-----------|---------|------|-------------|
| `c_agarose` | `42.0` | kg/m^3 | Agarose concentration (4.2% w/v for 7:3 blend at 6% total) |
| `c_chitosan` | `18.0` | kg/m^3 | Chitosan concentration (1.8% w/v) |
| `c_span80` | `20.0` | kg/m^3 | Span-80 surfactant concentration (~2% w/v) |
| `c_genipin` | `2.0` | mol/m^3 | Genipin crosslinker concentration (~2 mM) |
| `T_oil` | `363.15` | K | Oil bath temperature during emulsification (90 deg C) |
| `cooling_rate` | `0.167` | K/s | Cooling rate after emulsification (~10 deg C/min) |
| `T_crosslink` | `310.15` | K | Crosslinking incubation temperature (37 deg C) |
| `t_crosslink` | `86400.0` | s | Crosslinking incubation time (24 hours) |
| `phi_d` | `0.05` | - | Dispersed phase volume fraction |

### Derived properties

- `agarose_fraction` = `c_agarose / (c_agarose + c_chitosan)` -- mass fraction of agarose in the polymer blend
- `total_polymer` = `c_agarose + c_chitosan` -- total polymer concentration [kg/m^3]

## `[solver.level1]`

Numerical settings for the Population Balance Equation (PBE) solver.

| Parameter | Default | Unit | Description |
|-----------|---------|------|-------------|
| `n_bins` | `20` | - | Number of size bins in the PBE discretization |
| `d_min` | `1e-6` | m | Minimum droplet diameter (1 um) |
| `d_max` | `500e-6` | m | Maximum droplet diameter (500 um). Must exceed premix d32 by >= 3 sigma |
| `rtol` | `1e-6` | - | Relative tolerance for ODE integrator |
| `atol` | `1e-8` | - | Absolute tolerance for ODE integrator |

## `[solver.level2]`

Numerical settings for Level 2 (gelation and phase-field solver).

| Parameter | Default | Unit | Description |
|-----------|---------|------|-------------|
| `n_r` | `1000` | - | Number of radial grid points (1D solver) |
| `n_grid` | `128` | - | 2D grid side length (N x N Cahn-Hilliard solver) |
| `dt_initial` | `1e-4` | s | Initial time step |
| `dt_max` | `1.0` | s | Maximum adaptive time step |
| `arrest_exponent` | `2.5` | - | Gelation arrest exponent beta: mobility ~ (1 - alpha)^beta |

## `[solver.level3]`

Numerical settings for Level 3 (crosslinking ODE kinetics).

| Parameter | Default | Unit | Description |
|-----------|---------|------|-------------|
| `method` | `"Radau"` | - | SciPy ODE solver method (Radau, BDF, RK45, etc.) |
| `rtol` | `1e-8` | - | Relative tolerance |
| `atol` | `1e-10` | - | Absolute tolerance |

## `[optimization]`

Settings for the BoTorch Bayesian optimization campaign.

| Parameter | Default | Unit | Description |
|-----------|---------|------|-------------|
| `n_initial` | `20` | - | Number of Sobol quasi-random initial evaluations |
| `max_iterations` | `200` | - | Maximum optimization iterations |
| `convergence_tol` | `0.01` | - | Relative hypervolume convergence tolerance |

## Material Properties

Material properties are loaded from `data/properties.toml` and are not normally
edited in the run configuration. Key defaults (from `MaterialProperties` dataclass):

| Property | Default | Unit | Description |
|----------|---------|------|-------------|
| `rho_oil` | `850.0` | kg/m^3 | Oil density at 20 deg C (interpolated to T_oil) |
| `mu_oil` | `0.005` | Pa s | Oil dynamic viscosity at 90 deg C |
| `rho_aq` | `1020.0` | kg/m^3 | Aqueous phase density |
| `mu_d` | `1.0` | Pa s | Dispersed phase viscosity at T_oil |
| `sigma` | `5.0e-3` | N/m | Interfacial tension with Span-80 |
| `chi_0` | `0.497` | - | Flory-Huggins interaction parameter at reference T |
| `kappa_CH` | `5.0e-12` | J/m | Cahn-Hilliard gradient energy coefficient |
| `M_0` | `1.0e-9` | m^5/(J s) | Bare mobility (calibrated for 50-100 nm coarsening) |
| `T_gel` | `311.15` | K | Gelation temperature (~38 deg C) |
| `k_gel_0` | `1.0` | 1/s | Avrami rate prefactor |
| `n_avrami` | `2.5` | - | Avrami exponent |
| `k_xlink_0` | `2806.0` | m^3/(mol s) | Crosslinking Arrhenius prefactor |
| `E_a_xlink` | `52000.0` | J/mol | Crosslinking activation energy |
| `DDA` | `0.90` | - | Chitosan degree of deacetylation |
| `G_agarose_prefactor` | `3000.0` | Pa | Agarose gel modulus prefactor at 1% w/v |
| `G_agarose_exponent` | `2.2` | - | Power-law exponent for agarose modulus vs concentration |
| `f_bridge` | `0.4` | - | Fraction of genipin reactions producing elastically active crosslinks |
| `eta_coupling` | `-0.15` | - | IPN coupling coefficient (negative = synergistic) |
| `breakage_C3` | `0.0` | - | Alopaeus viscous correction constant (0 = disabled) |
| `r_fiber` | `1.5e-9` | m | Agarose fiber radius |

## CLI Overrides

Common parameters can be overridden directly on the command line:

```bash
python -m emulsim run --rpm 15000 --phi-d 0.08
python -m emulsim run configs/custom.toml --rpm 12000
```

The precedence order is: CLI flags > TOML file > dataclass defaults.
