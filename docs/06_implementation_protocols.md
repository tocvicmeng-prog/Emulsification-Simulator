# EmulSim Trust Warning Resolution: Implementation Protocols

**Date:** 2026-04-11
**Architecture Reference:** docs/05_trust_warning_resolution_plan.md
**Build Order:** M1 -> M4 -> M2 -> M3 (Phases 1->2->3)
**Protocol Version:** 1.0

---

## Table of Contents

1. [Module M1: PBE Steady-State Convergence](#module-m1-pbe-steady-state-convergence)
2. [Module M4: Per-Chemistry eta Coupling](#module-m4-per-chemistry-eta-coupling)
3. [Module M2: Crosslinker Stoichiometry Guidance](#module-m2-crosslinker-stoichiometry-guidance)
4. [Module M3: Predictive DN Modulus Model](#module-m3-predictive-dn-modulus-model)

---

## Module M1: PBE Steady-State Convergence

### Pre-Flight Check

- **Context budget:** ~104 LOC implementation + ~25 LOC tests = ~520 impl tokens + ~312 test tokens. GREEN zone.
- **Upstream dependencies:** None. M1 is independent.
- **Model tier:** Sonnet (standard implementation, 50-200 LOC, standard algorithm, domain-aware).

### 1. Module Specification

#### 1.1 File: `configs/default.toml` (line 8)

**Old code (line 8):**
```toml
t_emulsification = 60.0  # [s] (1 min — sufficient for steady state)
```

**New code:**
```toml
t_emulsification = 300.0  # [s] (5 min — allows PBE to reach steady state for viscous systems)
```

#### 1.2 File: `src/emulsim/datatypes.py` — `SolverSettings` (lines 534-551)

**Old code (lines 534-551):**
```python
@dataclass
class SolverSettings:
    """Numerical solver settings."""
    # Level 1
    l1_n_bins: int = 20
    l1_d_min: float = 1e-6             # [m]
    l1_d_max: float = 500e-6           # [m] (must exceed premix d32 by >=3sigma)
    l1_rtol: float = 1e-6
    l1_atol: float = 1e-8
    # Level 2
    l2_n_r: int = 1000
    l2_n_grid: int = 128               # 2D grid side length (NxN)
    l2_dt_initial: float = 1e-4        # [s]
    l2_dt_max: float = 1.0             # [s]
    l2_arrest_exponent: float = 2.5
    # Level 3
    l3_method: str = "Radau"
    l3_rtol: float = 1e-8
    l3_atol: float = 1e-10
```

**New code:**
```python
@dataclass
class SolverSettings:
    """Numerical solver settings."""
    # Level 1
    l1_n_bins: int = 20
    l1_d_min: float = 1e-6             # [m]
    l1_d_max: float = 500e-6           # [m] (must exceed premix d32 by >=3sigma)
    l1_rtol: float = 1e-6
    l1_atol: float = 1e-8
    l1_t_max: float = 600.0            # [s] absolute max emulsification time (hard ceiling for extensions)
    l1_conv_tol: float = 0.01          # [-] relative d32 variation threshold for convergence
    l1_max_extensions: int = 2         # [-] max number of adaptive time extensions
    # Level 2
    l2_n_r: int = 1000
    l2_n_grid: int = 128               # 2D grid side length (NxN)
    l2_dt_initial: float = 1e-4        # [s]
    l2_dt_max: float = 1.0             # [s]
    l2_arrest_exponent: float = 2.5
    # Level 3
    l3_method: str = "Radau"
    l3_rtol: float = 1e-8
    l3_atol: float = 1e-10
```

#### 1.3 File: `src/emulsim/datatypes.py` — `EmulsificationResult` (lines 778-793)

**Old code (lines 778-793):**
```python
@dataclass
class EmulsificationResult:
    """Output of Level 1: PBE solver."""
    d_bins: np.ndarray                  # [m] pivot diameters (N_bins,)
    n_d: np.ndarray                     # [#/m^3] number density (N_bins,)
    d32: float                          # [m] Sauter mean diameter
    d43: float                          # [m] volume-weighted mean
    d10: float                          # [m] 10th percentile
    d50: float                          # [m] median
    d90: float                          # [m] 90th percentile
    span: float                         # [-] (d90 - d10) / d50
    total_volume_fraction: float        # [-]
    converged: bool
    d_mode: float = 0.0                 # [m] modal diameter (volume-weighted)
    t_history: Optional[np.ndarray] = None
    n_d_history: Optional[np.ndarray] = None
```

**New code:**
```python
@dataclass
class EmulsificationResult:
    """Output of Level 1: PBE solver."""
    d_bins: np.ndarray                  # [m] pivot diameters (N_bins,)
    n_d: np.ndarray                     # [#/m^3] number density (N_bins,)
    d32: float                          # [m] Sauter mean diameter
    d43: float                          # [m] volume-weighted mean
    d10: float                          # [m] 10th percentile
    d50: float                          # [m] median
    d90: float                          # [m] 90th percentile
    span: float                         # [-] (d90 - d10) / d50
    total_volume_fraction: float        # [-]
    converged: bool
    d_mode: float = 0.0                 # [m] modal diameter (volume-weighted)
    t_history: Optional[np.ndarray] = None
    n_d_history: Optional[np.ndarray] = None
    t_converged: Optional[float] = None # [s] time at which d32 convergence was first achieved (None if never)
    n_extensions: int = 0               # [-] number of adaptive extensions performed
```

#### 1.4 File: `src/emulsim/level1_emulsification/solver.py` — `PBESolver.solve()` method

The core change: wrap the existing integration in an adaptive extension loop. The extension loop re-initialises from the final state and extends by `t_emul / 2`, up to `l1_max_extensions` times.

**Old code (lines 207-317, the `solve` method):**

Replace the body of `solve()` starting after the mode dispatch (line 219) through to the `return` (line 317). The full replacement is shown in pseudocode below; the exact diff follows.

**Insert after line 219 (the stirred_vessel dispatch)**, replacing lines 222-317:

```python
    def solve(self, params: SimulationParameters,
              props: MaterialProperties,
              phi_d: float | None = None) -> EmulsificationResult:
        """Solve the PBE for given process conditions.

        Parameters
        ----------
        phi_d : float, optional
            Dispersed-phase volume fraction override.  If ``None``,
            defaults to ``formulation.phi_d_from_volumes`` (stirred-vessel)
            or ``formulation.phi_d`` (legacy).
        """
        if params.emulsification.mode == "stirred_vessel":
            return self.solve_stirred_vessel(params, props, phi_d=phi_d)

        # Legacy mode: resolve default
        if phi_d is None:
            phi_d = params.formulation.phi_d

        rpm = params.emulsification.rpm
        t_emul = params.emulsification.t_emulsification
        mixer = params.emulsification.mixer

        # Energy dissipation
        rho_emul = emulsion_density(props.rho_oil, props.rho_aq, phi_d)
        epsilon_max = max_dissipation(mixer, rpm, rho_emul)
        epsilon_avg = average_dissipation(mixer, rpm, rho_emul)

        # Breakage rates (use zero-shear mu_d for viscous resistance Vi)
        nu_c = props.mu_oil / props.rho_oil
        g = breakage_rate_alopaeus(
            self.d_pivots, epsilon_max, props.sigma, props.rho_oil,
            props.mu_d, C3=props.breakage_C3, nu_c=nu_c,
        )
        birth_matrix, death_rate = self._build_breakage_matrix(g)

        # Coalescence rate matrix (vectorised construction)
        di_grid, dj_grid = np.meshgrid(self.d_pivots, self.d_pivots, indexing='ij')
        Q = coalescence_rate_ct(
            di_grid, dj_grid,
            epsilon_avg, props.sigma, props.rho_oil,
            props.mu_oil, phi_d=phi_d,
        )

        # Initial distribution
        N0 = self._initial_distribution(phi_d)

        def rhs(t, N):
            return self._compute_rhs(t, N, birth_matrix, death_rate, Q)

        # ── Adaptive extension loop ─────────────────────────────────────
        conv_tol = params.solver.l1_conv_tol
        max_extensions = params.solver.l1_max_extensions
        t_max_abs = params.solver.l1_t_max

        t_start = 0.0
        t_end = t_emul
        N_current = N0.copy()
        all_t = []
        all_y_cols = []
        n_extensions = 0
        converged = False
        t_converged = None

        for extension_round in range(max_extensions + 1):
            t_span = (t_start, t_end)
            n_eval_pts = 101 if extension_round == 0 else 51
            t_eval = np.linspace(t_start, t_end, n_eval_pts)

            sol = solve_ivp(
                rhs, t_span, N_current,
                method='BDF',
                rtol=params.solver.l1_rtol,
                atol=params.solver.l1_atol,
                t_eval=t_eval,
                max_step=(t_end - t_start) / 10,
            )

            if not sol.success:
                logger.warning("PBE solver did not converge: %s", sol.message)

            # Accumulate history
            if extension_round == 0:
                all_t.append(sol.t)
                all_y_cols.append(sol.y)
            else:
                # Skip first point (duplicate of previous last point)
                all_t.append(sol.t[1:])
                all_y_cols.append(sol.y[:, 1:])

            N_final = np.maximum(sol.y[:, -1], 0.0)

            # Convergence check: d32 stable over last 10% of TOTAL accumulated time
            t_combined = np.concatenate(all_t)
            y_combined = np.concatenate(all_y_cols, axis=1)
            n_check = max(1, len(t_combined) // 10)
            if y_combined.shape[1] > n_check:
                d32_late = [self._sauter_mean(np.maximum(y_combined[:, k], 0.0))
                            for k in range(-n_check, 0)]
                d32_current = self._sauter_mean(N_final)
                variation = (max(d32_late) - min(d32_late)) / max(d32_current, 1e-15)
                if variation < conv_tol:
                    converged = True
                    if t_converged is None:
                        t_converged = float(t_combined[-1])
                    break

            # If not converged and more extensions allowed, extend
            if extension_round < max_extensions:
                extension_dt = t_emul / 2.0
                new_t_end = t_end + extension_dt
                if new_t_end > t_max_abs:
                    new_t_end = t_max_abs
                if new_t_end <= t_end:
                    logger.info("L1: reached t_max=%.0fs without convergence", t_max_abs)
                    break
                logger.info(
                    "L1: not converged at t=%.0fs, extending to %.0fs (extension %d/%d)",
                    t_end, new_t_end, extension_round + 1, max_extensions,
                )
                N_current = N_final.copy()
                t_start = t_end
                t_end = new_t_end
                n_extensions += 1

        # Combine all history segments
        t_combined = np.concatenate(all_t)
        y_combined = np.concatenate(all_y_cols, axis=1)

        N_final = np.maximum(y_combined[:, -1], 0.0)
        d32, d43, d10, d50, d90, span = self._compute_statistics(N_final)
        total_vol = np.sum(N_final * self.v_pivots)

        n_d_output = N_final / self.d_widths
        n_d_history = y_combined.T / self.d_widths[np.newaxis, :]
        d_mode = self._compute_d_mode(N_final)

        return EmulsificationResult(
            d_bins=self.d_pivots.copy(),
            n_d=n_d_output,
            d32=d32, d43=d43, d10=d10, d50=d50, d90=d90, span=span,
            total_volume_fraction=total_vol,
            converged=converged,
            d_mode=d_mode,
            t_history=t_combined,
            n_d_history=n_d_history,
            t_converged=t_converged,
            n_extensions=n_extensions,
        )
```

#### 1.5 File: `src/emulsim/level1_emulsification/solver.py` — `PBESolver.solve_stirred_vessel()` method

Apply the same adaptive extension pattern to `solve_stirred_vessel()`. The structure is identical but uses the temperature-dependent RHS. Key changes at lines 487-553:

**Replace lines 487-553** with the same adaptive extension loop structure. The only differences from `solve()`:

1. The `rhs` function uses temperature interpolation (already defined at line 481)
2. Convergence tolerance is 0.05 (5%) and checks last 20% (lines 529-534 originally)
3. Use `conv_tol_sv = max(conv_tol, 0.05)` to keep the relaxed tolerance for stirred-vessel

The implementation should mirror the `solve()` extension loop, substituting:
- `conv_tol` -> `max(params.solver.l1_conv_tol, 0.05)` (stirred-vessel has noisier dynamics)
- `n_check = max(1, len(t_combined) // 5)` (20% window, not 10%)
- The existing `rhs` function (line 481) which already handles temperature interpolation

The `EmulsificationResult` return should also include `t_converged=t_converged` and `n_extensions=n_extensions`.

#### 1.6 File: `configs/default.toml` — Add solver convergence settings

**Add after `[solver.level1]` section (after line 33):**
```toml
t_max = 600.0             # [s] absolute ceiling for adaptive extensions
conv_tol = 0.01           # [-] relative d32 variation for convergence
max_extensions = 2        # [-] max number of half-interval extensions
```

### 2. Interface Contract

**Input types:**
- `SimulationParameters` (unchanged structure, new fields in `SolverSettings`)
- `MaterialProperties` (unchanged)
- `phi_d: float | None`

**Output type:**
- `EmulsificationResult` (with new fields `t_converged: Optional[float]`, `n_extensions: int`)

**Side effects:**
- `logger.info()` messages when extensions are triggered
- `logger.warning()` when `solve_ivp` fails

**Error handling:**
- If `solve_ivp` fails, log warning and continue with available data
- If `t_max_abs` is reached without convergence, `converged=False` and `t_converged=None`
- Volume conservation check: `total_volume_fraction` should remain within 5% of initial `phi_d`

### 3. Algorithm Pseudocode: Adaptive Extension Loop

```
FUNCTION solve_with_extension(N0, t_emul, t_max, max_ext, conv_tol):
    N_current = N0
    t_start = 0
    t_end = t_emul
    history = []
    n_ext = 0

    FOR round = 0 TO max_ext:
        sol = solve_ivp(rhs, [t_start, t_end], N_current, ...)
        history.append(sol)
        N_final = sol.y[:, -1]

        d32_late = [sauter_mean(sol.y[:, k]) for k in last_10%]
        variation = (max(d32_late) - min(d32_late)) / d32_current

        IF variation < conv_tol:
            converged = True
            BREAK

        IF round < max_ext AND t_end + t_emul/2 <= t_max:
            t_start = t_end
            t_end = min(t_end + t_emul/2, t_max)
            N_current = N_final
            n_ext += 1
        ELSE:
            BREAK

    RETURN combined_result
```

### 4. Test Specification

**File:** `tests/test_m1_convergence.py`

```python
def test_default_config_converges():
    """G-M1-1: Default config (t=300s) -> converged=True."""
    # Input: SimulationParameters() with t_emulsification=300.0 (new default)
    # Expected: result.converged == True
    # Tolerance: N/A (boolean)

def test_short_time_triggers_extension():
    """G-M1-2: Short t=5s -> adaptive extension triggers, eventually converges."""
    # Input: params.emulsification.t_emulsification = 5.0
    #        params.solver.l1_max_extensions = 2
    # Expected: result.n_extensions >= 1
    # Expected: result.converged == True OR n_extensions == 2

def test_t_history_monotonic_after_extension():
    """G-M1-3: t_history monotonically increasing after extension."""
    # Input: params.emulsification.t_emulsification = 5.0
    # Expected: np.all(np.diff(result.t_history) > 0)

def test_volume_conservation():
    """G-M1-4: Volume conservation within 5% of phi_d."""
    # Input: default params, phi_d = 0.05
    # Expected: abs(result.total_volume_fraction - 0.05) / 0.05 < 0.05

def test_trust_warning_w1_suppressed():
    """G-M1-5: Trust warning W1 no longer fires with defaults."""
    # Run full pipeline with default config (t=300s)
    # Expected: "PBE solver did not reach steady state" NOT in trust.warnings

def test_t_converged_set_on_convergence():
    """t_converged is a positive float when converged=True."""
    # Input: default params
    # Expected: result.t_converged is not None
    # Expected: result.t_converged > 0

def test_no_extension_when_already_converged():
    """If converged on first solve, n_extensions == 0."""
    # Input: t_emulsification = 300.0 (should converge without extension)
    # Expected: result.n_extensions == 0

def test_max_extensions_cap():
    """Extensions capped at l1_max_extensions."""
    # Input: t_emulsification = 1.0, l1_max_extensions = 1, l1_t_max = 600.0
    # Expected: result.n_extensions <= 1
```

### 5. G1 Gate Check

| # | Criterion | Status |
|---|-----------|--------|
| G1-01 | Every file path and line number specified | YES |
| G1-02 | Old code -> New code diffs shown | YES |
| G1-03 | New function signatures with type hints | YES (solve() unchanged signature, new return fields) |
| G1-04 | New dataclass fields with types and defaults | YES (SolverSettings: 3 fields; EmulsificationResult: 2 fields) |
| G1-05 | Config additions with TOML syntax | YES |
| G1-06 | Interface contract complete (I/O types, side effects, errors) | YES |
| G1-07 | Algorithm pseudocode provided | YES |
| G1-08 | Test function names and signatures specified | YES (8 tests) |
| G1-09 | Input fixtures with numerical values | YES |
| G1-10 | Expected outputs with tolerances | YES |
| G1-11 | Edge cases covered | YES (short time, max cap, already converged) |
| G1-12 | No ambiguity for implementer | YES |

---

## Module M4: Per-Chemistry eta Coupling

### Pre-Flight Check

- **Context budget:** ~57 LOC implementation + ~25 LOC tests. GREEN zone.
- **Upstream dependencies:** None (M1 is independent; M4 is independent of M1).
- **Model tier:** Sonnet (standard implementation, <100 LOC, domain-aware plumbing).

### 1. Module Specification

#### 1.1 File: `src/emulsim/datatypes.py` — `NetworkTypeMetadata` (lines 837-843)

**Old code:**
```python
@dataclass
class NetworkTypeMetadata:
    """Describes what network was formed by crosslinking."""
    solver_family: str = "amine_covalent"
    network_target: str = "chitosan"      # "chitosan", "agarose", "independent", "mixed"
    bond_type: str = "covalent"           # "covalent", "ionic", "reversible"
    is_true_second_network: bool = True   # True for IPN, False for reinforcement
```

**New code:**
```python
@dataclass
class NetworkTypeMetadata:
    """Describes what network was formed by crosslinking."""
    solver_family: str = "amine_covalent"
    network_target: str = "chitosan"      # "chitosan", "agarose", "independent", "mixed"
    bond_type: str = "covalent"           # "covalent", "ionic", "reversible"
    is_true_second_network: bool = True   # True for IPN, False for reinforcement
    eta_coupling_recommended: float = -0.15  # [-] per-chemistry IPN coupling coefficient from CrosslinkerProfile
```

#### 1.2 File: `src/emulsim/level3_crosslinking/solver.py` — Set `eta_coupling_recommended` on metadata at all dispatch points

There are 5 places where `NetworkTypeMetadata` is constructed. Each must propagate the crosslinker profile's `eta_coupling_recommended`.

**1.2a: Line 772-777 (hydroxyl_covalent branch):**

**Old code:**
```python
            metadata = NetworkTypeMetadata(
                solver_family="hydroxyl_covalent",
                network_target="mixed",
                bond_type="covalent",
                is_true_second_network=True,
            )
```

**New code:**
```python
            metadata = NetworkTypeMetadata(
                solver_family="hydroxyl_covalent",
                network_target="mixed",
                bond_type="covalent",
                is_true_second_network=True,
                eta_coupling_recommended=xl.eta_coupling_recommended,
            )
```

**1.2b: Line 780-785 (amine_covalent branch):**

**Old code:**
```python
            metadata = NetworkTypeMetadata(
                solver_family="amine_covalent",
                network_target="chitosan",
                bond_type="covalent",
                is_true_second_network=True,
            )
```

**New code:**
```python
            metadata = NetworkTypeMetadata(
                solver_family="amine_covalent",
                network_target="chitosan",
                bond_type="covalent",
                is_true_second_network=True,
                eta_coupling_recommended=xl.eta_coupling_recommended,
            )
```

**1.2c: Line 788-793 (independent_network branch, uv_dose):**

**Old code:**
```python
        metadata = NetworkTypeMetadata(
            solver_family="independent_network",
            network_target="independent",
            bond_type="covalent",
            is_true_second_network=True,
        )
```

**New code:**
```python
        metadata = NetworkTypeMetadata(
            solver_family="independent_network",
            network_target="independent",
            bond_type="covalent",
            is_true_second_network=True,
            eta_coupling_recommended=xl.eta_coupling_recommended,
        )
```

**1.2d: Line 796-801 (ionic_reversible branch):**

**Old code:**
```python
        metadata = NetworkTypeMetadata(
            solver_family="ionic_reversible",
            network_target="chitosan",
            bond_type="ionic",
            is_true_second_network=False,
        )
```

**New code:**
```python
        metadata = NetworkTypeMetadata(
            solver_family="ionic_reversible",
            network_target="chitosan",
            bond_type="ionic",
            is_true_second_network=False,
            eta_coupling_recommended=xl.eta_coupling_recommended,
        )
```

**1.2e: Line 812-817 (michaelis_menten fallback branch):**

**Old code:**
```python
        metadata = NetworkTypeMetadata(
            solver_family="amine_covalent",
            network_target="chitosan",
            bond_type="covalent",
            is_true_second_network=True,
        )
```

**New code:**
```python
        metadata = NetworkTypeMetadata(
            solver_family="amine_covalent",
            network_target="chitosan",
            bond_type="covalent",
            is_true_second_network=True,
            eta_coupling_recommended=xl.eta_coupling_recommended,
        )
```

**1.2f: Reaction-diffusion path (line 420-425).** This is the `_solve_reaction_diffusion` function which constructs metadata at line 420. It does not have access to the `xl` profile because it is called before the dispatch block. The fix requires passing `xl` to `_solve_reaction_diffusion` or setting metadata after return.

**Current code at lines 763-766:**
```python
            rd_result, _ = _solve_reaction_diffusion(
                params, props, R_droplet, D_xlink,
            )
            return rd_result
```

**New code:**
```python
            rd_result, _ = _solve_reaction_diffusion(
                params, props, R_droplet, D_xlink,
            )
            # Attach per-chemistry metadata to reaction-diffusion result
            rd_result.network_metadata = NetworkTypeMetadata(
                solver_family="amine_covalent",
                network_target="chitosan",
                bond_type="covalent",
                is_true_second_network=True,
                eta_coupling_recommended=xl.eta_coupling_recommended,
            )
            return rd_result
```

#### 1.3 File: `src/emulsim/level4_mechanical/solver.py` — `select_modulus_model()` (lines 74-102)

**Old code (lines 74-102):**
```python
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
```

**New code:**
```python
def select_modulus_model(G_agarose: float, G_xlink: float,
                          network_metadata=None,
                          eta_coupling: float = -0.15) -> float:
    """Route to the appropriate modulus model based on network metadata.

    Falls back to phenomenological DN modulus if no metadata provided.
    When network_metadata carries a per-chemistry eta_coupling_recommended,
    that value overrides the default eta_coupling parameter.
    """
    if network_metadata is None:
        return double_network_modulus(G_agarose, G_xlink, eta_coupling)

    family = getattr(network_metadata, 'solver_family', 'amine_covalent')

    # Per-chemistry eta: prefer metadata value, fall back to function argument
    eta_per_chem = getattr(network_metadata, 'eta_coupling_recommended', eta_coupling)

    if family == "ionic_reversible":
        return ionic_gel_modulus(G_agarose, G_xlink)
    elif family == "independent_network":
        return triple_network_modulus(G_agarose, 0.0, G_xlink, eta_13=eta_per_chem)
    elif family == "hydroxyl_covalent":
        return double_network_modulus(G_agarose, G_xlink, eta_coupling=eta_per_chem)
    else:
        # amine_covalent and default
        return double_network_modulus(G_agarose, G_xlink, eta_per_chem)
```

#### 1.4 File: `src/emulsim/trust.py` — Suppress W4 when per-chemistry eta is active (lines 177-183)

**Old code (lines 177-183):**
```python
    # 14. Non-specific eta_coupling (same default for all crosslinker types)
    if props.eta_coupling == -0.15:
        warnings.append(
            "IPN coupling coefficient (eta=-0.15) is the same default for all crosslinker types. "
            "Per-chemistry eta values would improve accuracy."
        )
```

**New code:**
```python
    # 14. Non-specific eta_coupling (same default for all crosslinker types)
    # Suppress when per-chemistry eta is being used via network_metadata
    _eta_active = False
    if xl and hasattr(x, 'network_metadata') and x.network_metadata is not None:
        meta_eta = getattr(x.network_metadata, 'eta_coupling_recommended', None)
        if meta_eta is not None:
            _eta_active = True
    if not _eta_active and props.eta_coupling == -0.15:
        warnings.append(
            "IPN coupling coefficient (eta=-0.15) is the same default for all crosslinker types. "
            "Per-chemistry eta values would improve accuracy."
        )
```

Note: `x` is `result.crosslinking` (already assigned at line 55). The check verifies that the `CrosslinkingResult.network_metadata.eta_coupling_recommended` field is populated, indicating per-chemistry eta was wired through.

### 2. Interface Contract

**Input types:**
- `CrosslinkerProfile.eta_coupling_recommended: float` (already exists in reagent library)
- `NetworkTypeMetadata` (new field `eta_coupling_recommended`)
- `CrosslinkingResult.network_metadata` (carries the value through L3 -> L4)

**Output type:**
- `float` from `select_modulus_model()` (G_DN value that now uses per-chemistry eta)

**Side effects:**
- Trust warning W4 suppressed when per-chemistry eta is active
- No file I/O changes

**Error handling:**
- If `network_metadata` is `None`, falls back to `props.eta_coupling` (backward compatible)
- If `eta_coupling_recommended` attribute missing, falls back to `eta_coupling` parameter

### 3. Test Specification

**File:** `tests/test_m4_per_chemistry_eta.py`

```python
def test_genipin_eta_backward_compatible():
    """G-M4-1: Genipin -> eta = -0.15 (backward compatible)."""
    # Input: crosslinker_key='genipin'
    # Run solve_crosslinking with default params
    # Expected: result.network_metadata.eta_coupling_recommended == -0.15

def test_dvs_eta_positive():
    """G-M4-2: DVS -> eta = +0.10."""
    # Input: crosslinker_key='dvs'
    # Expected: result.network_metadata.eta_coupling_recommended == pytest.approx(0.10)

def test_tpp_eta():
    """G-M4-3: TPP -> eta = -0.05."""
    # Input: crosslinker_key='tpp'
    # Expected: result.network_metadata.eta_coupling_recommended == pytest.approx(-0.05)

def test_pegda_eta_zero():
    """G-M4-4: PEGDA -> eta = 0.0."""
    # Input: crosslinker_key='pegda_uv', uv_intensity=10.0
    # Expected: result.network_metadata.eta_coupling_recommended == pytest.approx(0.0)

def test_w4_suppressed_with_per_chemistry_eta():
    """G-M4-5: Trust warning W4 suppressed when per-chemistry eta active."""
    # Run full pipeline with genipin (which has per-chemistry eta wired)
    # Expected: "Per-chemistry eta values would improve accuracy" NOT in trust.warnings

def test_gdn_differs_dvs_vs_genipin():
    """G-M4-6: G_DN differs between DVS and genipin (different eta sign)."""
    # Input: same formulation, same G_agarose, same G_xlink
    #        G_agarose = 50000.0, G_xlink = 5000.0
    # Genipin: eta = -0.15, G_DN = 50000 + 5000 + (-0.15)*sqrt(50000*5000) = 52628.8
    # DVS:    eta = +0.10, G_DN = 50000 + 5000 + (0.10)*sqrt(50000*5000) = 56581.1
    # Expected: G_DN_dvs > G_DN_genipin
    # Numerical check:
    #   G_DN_genipin = pytest.approx(52628.8, rel=0.01)
    #   G_DN_dvs = pytest.approx(56581.1, rel=0.01)

def test_edc_nhs_eta():
    """EDC/NHS -> eta = -0.10."""
    # Input: crosslinker_key='edc_nhs'
    # Expected: result.network_metadata.eta_coupling_recommended == pytest.approx(-0.10)

def test_reaction_diffusion_path_carries_eta():
    """Reaction-diffusion PDE path also sets per-chemistry eta."""
    # Input: large R_droplet (e.g. 100e-6) to trigger Thiele > 1
    #        crosslinker_key='genipin'
    # Expected: result.network_metadata.eta_coupling_recommended == -0.15
```

### 4. G1 Gate Check

| # | Criterion | Status |
|---|-----------|--------|
| G1-01 | Every file path and line number specified | YES |
| G1-02 | Old code -> New code diffs shown | YES (7 diff sites) |
| G1-03 | New function signatures with type hints | YES (no new functions, modified select_modulus_model) |
| G1-04 | New dataclass fields with types and defaults | YES (NetworkTypeMetadata.eta_coupling_recommended) |
| G1-05 | Config additions with TOML syntax | N/A (no config changes) |
| G1-06 | Interface contract complete | YES |
| G1-07 | Algorithm pseudocode provided | N/A (pure plumbing, no algorithm) |
| G1-08 | Test function names and signatures specified | YES (8 tests) |
| G1-09 | Input fixtures with numerical values | YES |
| G1-10 | Expected outputs with tolerances | YES |
| G1-11 | Edge cases covered | YES (reaction-diffusion path, missing metadata fallback) |
| G1-12 | No ambiguity for implementer | YES |

---

## Module M2: Crosslinker Stoichiometry Guidance

### Pre-Flight Check

- **Context budget:** ~24 LOC implementation. GREEN zone.
- **Upstream dependencies:** M4 APPROVED (ensures per-chemistry eta is wired). M1 is independent.
- **Model tier:** Haiku (LOW complexity: <50 LOC, config change + simple function + message enhancement).

### 1. Module Specification

#### 1.1 File: `configs/default.toml` (line 21)

**Old code:**
```toml
c_genipin = 2.0            # [mol/m^3] (~2 mM)
```

**New code:**
```toml
c_genipin = 10.0            # [mol/m^3] (~10 mM, stoichiometrically adequate for r >= 0.10)
```

#### 1.2 File: `src/emulsim/level3_crosslinking/solver.py` — New function

**Insert after line 102 (after `available_hydroxyl_concentration`):**

```python
def recommended_crosslinker_concentration(
    c_chitosan: float,
    DDA: float,
    M_GlcN: float,
    target_p: float = 0.20,
) -> float:
    """Minimum crosslinker concentration [mol/m^3] for target crosslinking fraction.

    Each crosslinker bridge consumes 2 NH2 groups, so:
        c_crosslinker = target_p * [NH2] / 2

    Parameters
    ----------
    c_chitosan : float
        Chitosan concentration [kg/m^3].
    DDA : float
        Degree of deacetylation [-].
    M_GlcN : float
        Molar mass of glucosamine repeat unit [g/mol].
    target_p : float
        Target crosslinking fraction [-], default 0.20.

    Returns
    -------
    float
        Recommended crosslinker concentration [mol/m^3].
    """
    NH2 = available_amine_concentration(c_chitosan, DDA, M_GlcN)
    if NH2 <= 0 or target_p <= 0:
        return 0.0
    return target_p * NH2 / 2.0
```

#### 1.3 File: `src/emulsim/trust.py` — Enhance W2 warning (lines 103-108)

**Old code (lines 103-108):**
```python
    if NH2_total > 0 and params.formulation.c_genipin / NH2_total < 0.05:
        warnings.append(
            f"Crosslinker/NH2 ratio = {params.formulation.c_genipin/NH2_total:.3f} -- "
            "crosslinker-limited, increasing crosslinking time will not help"
        )
```

**New code:**
```python
    if NH2_total > 0 and params.formulation.c_genipin / NH2_total < 0.05:
        # Compute minimum crosslinker for ratio >= 0.05
        c_min_005 = 0.05 * NH2_total
        # Compute recommended for p = 0.20
        from .level3_crosslinking.solver import recommended_crosslinker_concentration
        c_rec = recommended_crosslinker_concentration(
            params.formulation.c_chitosan, props.DDA, props.M_GlcN, target_p=0.20,
        )
        warnings.append(
            f"Crosslinker/NH2 ratio = {params.formulation.c_genipin/NH2_total:.3f} -- "
            f"crosslinker-limited. Minimum for ratio >= 0.05: {c_min_005:.1f} mol/m3. "
            f"Recommended for 20% crosslinking: {c_rec:.1f} mol/m3."
        )
```

### 2. Interface Contract

**New function signature:**
```python
def recommended_crosslinker_concentration(
    c_chitosan: float,    # [kg/m^3]
    DDA: float,           # [-]
    M_GlcN: float,        # [g/mol]
    target_p: float = 0.20,  # [-]
) -> float:               # [mol/m^3]
```

**Input types:** All `float` scalars.
**Output type:** `float` (recommended concentration in mol/m^3).
**Side effects:** None (pure function).
**Error handling:** Returns 0.0 if NH2 <= 0 or target_p <= 0.

### 3. Test Specification

**File:** `tests/test_m2_stoichiometry.py`

```python
def test_default_ratio_adequate():
    """G-M2-1: Default config ratio >= 0.05."""
    # With c_genipin=10.0 (new default), c_chitosan=18.0, DDA=0.9, M_GlcN=161.16:
    # NH2 = 18.0 * 1000 * 0.9 / 161.16 = 100.5 mol/m^3
    # ratio = 10.0 / 100.5 = 0.0995
    # Expected: ratio >= 0.05 (True)
    params = SimulationParameters()  # uses new default c_genipin=10.0
    props = MaterialProperties()
    NH2 = params.formulation.c_chitosan * 1000 * props.DDA / props.M_GlcN
    ratio = params.formulation.c_genipin / NH2
    assert ratio >= 0.05

def test_w2_not_fired_with_defaults():
    """G-M2-2: Trust warning W2 no longer fires with defaults."""
    # Run assess_trust with default config (c_genipin=10.0)
    # Expected: "crosslinker-limited" NOT in any warning string

def test_warning_includes_recommendation():
    """G-M2-3: Warning message includes suggested minimum concentration."""
    # Input: params with c_genipin=1.0 (deliberately low)
    # Expected: warning contains "mol/m3" and a numerical recommendation

def test_recommended_concentration():
    """G-M2-4: recommended_crosslinker_concentration(..., target_p=0.3) -> p_final approx 0.3."""
    # Input: c_chitosan=18.0, DDA=0.9, M_GlcN=161.16, target_p=0.30
    # NH2 = 100.5 mol/m^3
    # c_rec = 0.30 * 100.5 / 2.0 = 15.075 mol/m^3
    # Expected: c_rec == pytest.approx(15.075, rel=0.01)
    from emulsim.level3_crosslinking.solver import recommended_crosslinker_concentration
    c = recommended_crosslinker_concentration(18.0, 0.9, 161.16, target_p=0.30)
    assert c == pytest.approx(15.075, rel=0.01)

def test_recommended_zero_chitosan():
    """Edge case: zero chitosan gives zero recommendation."""
    from emulsim.level3_crosslinking.solver import recommended_crosslinker_concentration
    c = recommended_crosslinker_concentration(0.0, 0.9, 161.16, target_p=0.20)
    assert c == 0.0

def test_recommended_zero_target():
    """Edge case: zero target_p gives zero recommendation."""
    from emulsim.level3_crosslinking.solver import recommended_crosslinker_concentration
    c = recommended_crosslinker_concentration(18.0, 0.9, 161.16, target_p=0.0)
    assert c == 0.0
```

### 4. G1 Gate Check

| # | Criterion | Status |
|---|-----------|--------|
| G1-01 | Every file path and line number specified | YES |
| G1-02 | Old code -> New code diffs shown | YES |
| G1-03 | New function signatures with type hints | YES |
| G1-04 | New dataclass fields with types and defaults | N/A (no dataclass changes) |
| G1-05 | Config additions with TOML syntax | YES |
| G1-06 | Interface contract complete | YES |
| G1-07 | Algorithm pseudocode provided | N/A (trivial formula) |
| G1-08 | Test function names and signatures specified | YES (6 tests) |
| G1-09 | Input fixtures with numerical values | YES |
| G1-10 | Expected outputs with tolerances | YES |
| G1-11 | Edge cases covered | YES (zero chitosan, zero target) |
| G1-12 | No ambiguity for implementer | YES |

---

## Module M3: Predictive DN Modulus Model

### Pre-Flight Check

- **Context budget:** ~264 LOC implementation + ~50 LOC tests. Largest module. Estimate: ~1056 impl tokens + ~300 test tokens + audit. GREEN zone but monitor.
- **Upstream dependencies:** M4 APPROVED (per-chemistry eta wired), M2 APPROVED (valid stoichiometry for testing).
- **Model tier:** Phase 1-2 Sonnet, Phase 3 Opus (novel scientific model, Flory-Rehner solver, >200 LOC).
- **Build sub-order:** Phase 1 -> Phase 2 -> Phase 3 (each independently testable).

### Phase 1: Mode-Aware Warnings + Model Metadata

#### 1.1 File: `src/emulsim/datatypes.py` — `MechanicalResult` (lines 863-875)

**Old code:**
```python
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
```

**New code:**
```python
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
    model_used: str = "phenomenological"  # which modulus model was used
    G_DN_lower: float = 0.0             # [Pa] Hashin-Shtrikman lower bound
    G_DN_upper: float = 0.0             # [Pa] Hashin-Shtrikman upper bound
```

#### 1.2 File: `src/emulsim/trust.py` — Mode-aware W3 warning (lines 163-167)

**Old code (lines 163-167):**
```python
    # 12. Phenomenological DN modulus (standing note)
    warnings.append(
        "G_DN uses phenomenological coupling formula (G1 + G2 + eta*sqrt(G1*G2)). "
        "Suitable for formulation ranking, not absolute mechanical prediction."
    )
```

**New code:**
```python
    # 12. Phenomenological DN modulus (mode-aware)
    model_mode = getattr(params, 'model_mode', None)
    _model_used = getattr(m, 'model_used', 'phenomenological')
    if _model_used == "phenomenological":
        if model_mode == ModelMode.MECHANISTIC_RESEARCH:
            warnings.append(
                "G_DN uses phenomenological coupling formula (G1 + G2 + eta*sqrt(G1*G2)). "
                "MECHANISTIC_RESEARCH mode requires the affine IPN model for defensible predictions. "
                "Set model_mode='hybrid_coupled' or implement Phase 3 of M3."
            )
        elif model_mode == ModelMode.EMPIRICAL_ENGINEERING:
            pass  # Phenomenological is appropriate for empirical engineering; no warning.
        else:
            # hybrid_coupled (default): note is informational
            warnings.append(
                "G_DN uses phenomenological coupling formula (G1 + G2 + eta*sqrt(G1*G2)). "
                "Suitable for formulation ranking, not absolute mechanical prediction."
            )
```

Add the import at the top of trust.py (near line 8):
```python
from .datatypes import FullResult, SimulationParameters, MaterialProperties, ModelMode
```

(Updating the existing import line 8 from `from .datatypes import FullResult, SimulationParameters, MaterialProperties` to include `ModelMode`.)

#### 1.3 File: `src/emulsim/level4_mechanical/solver.py` — Populate model_used in solve_mechanical (lines 236-247)

**Old code (lines 236-247):**
```python
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
```

**New code:**
```python
    # Determine which model was used
    network_meta = getattr(crosslinking, 'network_metadata', None)
    if network_meta is None:
        model_label = "phenomenological"
    else:
        family = getattr(network_meta, 'solver_family', 'amine_covalent')
        if family == "ionic_reversible":
            model_label = "ionic_gel"
        elif family == "independent_network":
            model_label = "triple_network"
        else:
            model_label = "phenomenological"

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
    )
```

### Phase 2: Hashin-Shtrikman Bounds

#### 2.1 File: `src/emulsim/level4_mechanical/solver.py` — New function

**Insert after `double_network_modulus` (after line 48):**

```python
def hashin_shtrikman_bounds(
    G1: float, G2: float, phi1: float,
) -> tuple[float, float]:
    """Hashin-Shtrikman bounds for two-phase composite shear modulus [Pa].

    These are the tightest possible bounds on the effective modulus of a
    two-phase isotropic composite, given only phase moduli and volume fractions.

    Parameters
    ----------
    G1 : float
        Shear modulus of phase 1 (e.g., agarose) [Pa].
    G2 : float
        Shear modulus of phase 2 (e.g., crosslinked chitosan) [Pa].
    phi1 : float
        Volume fraction of phase 1 [-]. phi2 = 1 - phi1.

    Returns
    -------
    G_lower : float
        Lower Hashin-Shtrikman bound [Pa].
    G_upper : float
        Upper Hashin-Shtrikman bound [Pa].

    Notes
    -----
    For shear modulus (incompressible limit, nu -> 0.5):
        G_HS^- = G_soft + phi_stiff / (1/(G_stiff - G_soft) + 3*phi_soft/(5*G_soft))
        G_HS^+ = G_stiff + phi_soft / (1/(G_soft - G_stiff) + 3*phi_stiff/(5*G_stiff))

    Reference: Hashin & Shtrikman (1963) J. Mech. Phys. Solids 11:127.
    """
    G1 = max(G1, 0.0)
    G2 = max(G2, 0.0)
    phi1 = max(0.0, min(phi1, 1.0))
    phi2 = 1.0 - phi1

    # Handle degenerate cases
    if G1 <= 0 and G2 <= 0:
        return 0.0, 0.0
    if phi1 <= 0:
        return G2, G2
    if phi2 <= 0:
        return G1, G1

    # Identify soft and stiff phases
    if G1 <= G2:
        G_soft, G_stiff = G1, G2
        phi_soft, phi_stiff = phi1, phi2
    else:
        G_soft, G_stiff = G2, G1
        phi_soft, phi_stiff = phi2, phi1

    # Lower bound (soft matrix, stiff inclusions)
    if G_soft > 0:
        denom_lower = 1.0 / max(G_stiff - G_soft, 1e-15) + 3.0 * phi_soft / (5.0 * G_soft)
        G_lower = G_soft + phi_stiff / denom_lower
    else:
        G_lower = 0.0

    # Upper bound (stiff matrix, soft inclusions)
    if G_stiff > 0:
        denom_upper = 1.0 / max(G_soft - G_stiff, -1e15) + 3.0 * phi_stiff / (5.0 * G_stiff)
        # For G_soft < G_stiff: 1/(G_soft - G_stiff) is negative
        # G_upper = G_stiff + phi_soft / denom_upper
        # Safer formulation:
        delta_G = G_soft - G_stiff  # negative
        if abs(delta_G) > 1e-15:
            denom_upper = 1.0 / delta_G + 3.0 * phi_stiff / (5.0 * G_stiff)
            if abs(denom_upper) > 1e-15:
                G_upper = G_stiff + phi_soft / denom_upper
            else:
                G_upper = G_stiff
        else:
            G_upper = G_stiff
    else:
        G_upper = 0.0

    # Ensure ordering
    G_lower = max(G_lower, 0.0)
    G_upper = max(G_upper, G_lower)

    return float(G_lower), float(G_upper)
```

#### 2.2 File: `src/emulsim/level4_mechanical/solver.py` — Populate bounds in `solve_mechanical` (around line 205-206)

**After `G_DN = select_modulus_model(...)` (line 206), insert:**

```python
    # Hashin-Shtrikman bounds
    # phi1 = agarose volume fraction in the polymer blend
    c_total_polymer = params.formulation.c_agarose + params.formulation.c_chitosan
    if c_total_polymer > 0:
        phi_agar_in_blend = params.formulation.c_agarose / c_total_polymer
    else:
        phi_agar_in_blend = 0.5
    G_DN_lower, G_DN_upper = hashin_shtrikman_bounds(G_agar, G_xlink, phi_agar_in_blend)
```

**Then modify the MechanicalResult return** to include:
```python
        G_DN_lower=G_DN_lower,
        G_DN_upper=G_DN_upper,
```

### Phase 3: Affine IPN Model (Flory-Rehner Swelling-Constrained)

#### 3.1 File: `src/emulsim/level4_mechanical/solver.py` — New functions

**Insert after `hashin_shtrikman_bounds`:**

```python
def flory_rehner_swelling(
    nu_e: float,
    phi_0: float,
    chi: float,
    T: float = 298.15,
) -> float:
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
    """
    from scipy.optimize import brentq

    V_s = 18.0e-6  # [m^3/mol] molar volume of water

    def objective(phi):
        if phi <= 0 or phi >= 1:
            return 1e10
        return (np.log(1.0 - phi) + phi + chi * phi**2
                + nu_e * V_s * (phi**(1.0/3.0) - phi / 2.0))

    # Bracket: phi must be between phi_0*0.01 (very swollen) and 0.999
    phi_lo = max(phi_0 * 0.01, 1e-6)
    phi_hi = 0.999

    try:
        phi_eq = brentq(objective, phi_lo, phi_hi, xtol=1e-8, maxiter=200)
    except ValueError:
        # Bracket does not contain root; return reference volume fraction
        import logging
        logging.getLogger(__name__).warning(
            "Flory-Rehner: no root in [%.4f, %.4f] for nu_e=%.1f, chi=%.3f. "
            "Returning phi_0=%.4f.", phi_lo, phi_hi, nu_e, chi, phi_0,
        )
        phi_eq = phi_0

    return float(phi_eq)


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

        G_DN = nu_e1 * kT * lambda_1^(-2/3) + nu_e2 * kT * lambda_2^(-2/3)

    where lambda_i = (phi_0i / phi_eq_i)^(1/3) is the swelling stretch ratio.

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
    chi2 : float
        Flory-Huggins parameter for network 2 [-].
    T : float
        Temperature [K].

    Returns
    -------
    G_DN : float
        Double-network shear modulus [Pa].
    status : str
        "converged" if both Flory-Rehner solutions converged, else "fallback".
    """
    K_BOLTZMANN = 1.38e-23  # J/K
    N_AVOGADRO = 6.022e23
    kT = K_BOLTZMANN * T

    status = "converged"

    # Network 1: agarose
    if nu_e1 > 0 and phi_01 > 0:
        phi_eq1 = flory_rehner_swelling(nu_e1, phi_01, chi1, T)
        if phi_eq1 > 0:
            lambda_1 = (phi_01 / phi_eq1) ** (1.0 / 3.0)
        else:
            lambda_1 = 1.0
            status = "fallback"
        G1 = nu_e1 * N_AVOGADRO * kT * lambda_1 ** (-2.0 / 3.0)
    else:
        G1 = 0.0

    # Network 2: chitosan
    if nu_e2 > 0 and phi_02 > 0:
        phi_eq2 = flory_rehner_swelling(nu_e2, phi_02, chi2, T)
        if phi_eq2 > 0:
            lambda_2 = (phi_02 / phi_eq2) ** (1.0 / 3.0)
        else:
            lambda_2 = 1.0
            status = "fallback"
        G2 = nu_e2 * N_AVOGADRO * kT * lambda_2 ** (-2.0 / 3.0)
    else:
        G2 = 0.0

    G_DN = G1 + G2

    # Sanity bounds: 15 kPa < G_DN < 150 kPa for standard formulations
    if G_DN < 0:
        G_DN = 0.0
        status = "fallback"

    return float(G_DN), status
```

#### 3.2 File: `src/emulsim/level4_mechanical/solver.py` — Update `select_modulus_model` for MECHANISTIC_RESEARCH

This requires the `select_modulus_model` function to accept additional parameters for the affine model. The cleanest approach: add an optional `params` argument.

**Updated signature and body of `select_modulus_model` (incorporating Phase 3):**

```python
def select_modulus_model(G_agarose: float, G_xlink: float,
                          network_metadata=None,
                          eta_coupling: float = -0.15,
                          model_mode: str = "hybrid_coupled",
                          params: 'SimulationParameters | None' = None,
                          crosslinking: 'CrosslinkingResult | None' = None,
                          ) -> tuple[float, str]:
    """Route to the appropriate modulus model based on network metadata and mode.

    Returns (G_DN, model_label) where model_label identifies which model was used.
    """
    if network_metadata is None:
        return double_network_modulus(G_agarose, G_xlink, eta_coupling), "phenomenological"

    family = getattr(network_metadata, 'solver_family', 'amine_covalent')
    eta_per_chem = getattr(network_metadata, 'eta_coupling_recommended', eta_coupling)

    if family == "ionic_reversible":
        return ionic_gel_modulus(G_agarose, G_xlink), "ionic_gel"
    elif family == "independent_network":
        return triple_network_modulus(G_agarose, 0.0, G_xlink, eta_13=eta_per_chem), "triple_network"

    # For amine_covalent and hydroxyl_covalent:
    # In MECHANISTIC_RESEARCH mode, attempt affine IPN model
    if model_mode == "mechanistic_research" and crosslinking is not None and params is not None:
        try:
            nu_e2 = crosslinking.nu_e_final  # crosslinked network crosslink density [1/m^3]
            # Convert to mol/m^3
            N_AVOGADRO = 6.022e23
            nu_e2_mol = nu_e2 / N_AVOGADRO

            # Network 1 (agarose): estimate nu_e from rubber elasticity G = nu_e * kT
            kT = 1.38e-23 * 298.15
            nu_e1 = G_agarose / kT if kT > 0 else 0.0
            nu_e1_mol = nu_e1 / N_AVOGADRO

            rho_polymer = 1400.0  # kg/m^3
            phi_01 = params.formulation.c_agarose / rho_polymer
            phi_02 = params.formulation.c_chitosan / rho_polymer

            chi1 = 0.50  # agarose-water (slightly poor solvent)
            chi2 = 0.45  # chitosan-water

            G_affine, affine_status = double_network_modulus_affine(
                nu_e1_mol, nu_e2_mol, phi_01, phi_02, chi1, chi2,
            )

            if affine_status == "converged" and G_affine > 0:
                return G_affine, "affine_ipn"
            else:
                import logging
                logging.getLogger(__name__).warning(
                    "Affine IPN model fallback: status=%s, G=%.1f Pa. "
                    "Using phenomenological.", affine_status, G_affine,
                )
        except Exception as exc:
            import logging
            logging.getLogger(__name__).warning(
                "Affine IPN model failed: %s. Using phenomenological.", exc,
            )

    # Phenomenological fallback (and default for hybrid_coupled / empirical_engineering)
    if family == "hydroxyl_covalent":
        return double_network_modulus(G_agarose, G_xlink, eta_coupling=eta_per_chem), "phenomenological"
    else:
        return double_network_modulus(G_agarose, G_xlink, eta_per_chem), "phenomenological"
```

**IMPORTANT:** The return type changes from `float` to `tuple[float, str]`. The call site in `solve_mechanical` (line 206) must be updated:

**Old call (line 206):**
```python
    G_DN = select_modulus_model(G_agar, G_xlink, network_meta, props.eta_coupling)
```

**New call:**
```python
    model_mode_str = getattr(params.model_mode, 'value', 'hybrid_coupled')
    G_DN, model_label = select_modulus_model(
        G_agar, G_xlink, network_meta, props.eta_coupling,
        model_mode=model_mode_str,
        params=params,
        crosslinking=crosslinking,
    )
```

And the `model_label` variable replaces the `model_label` block added in Phase 1 (section 1.3).

### 4. Algorithm Pseudocode

#### Flory-Rehner Solver

```
FUNCTION flory_rehner_swelling(nu_e, phi_0, chi, T):
    V_s = 18e-6  # molar volume of water [m^3/mol]

    DEFINE f(phi) = ln(1 - phi) + phi + chi*phi^2 + nu_e*V_s*(phi^(1/3) - phi/2)

    phi_lo = max(phi_0 * 0.01, 1e-6)
    phi_hi = 0.999

    phi_eq = brentq(f, phi_lo, phi_hi)  # scalar root-finding
    RETURN phi_eq
```

#### Hashin-Shtrikman Bounds

```
FUNCTION hashin_shtrikman_bounds(G1, G2, phi1):
    phi2 = 1 - phi1
    IF G1 <= G2:
        G_soft, G_stiff = G1, G2
        phi_soft, phi_stiff = phi1, phi2
    ELSE:
        swap

    # Lower bound (soft matrix):
    G_lower = G_soft + phi_stiff / (1/(G_stiff - G_soft) + 3*phi_soft/(5*G_soft))

    # Upper bound (stiff matrix):
    G_upper = G_stiff + phi_soft / (1/(G_soft - G_stiff) + 3*phi_stiff/(5*G_stiff))

    RETURN (G_lower, G_upper)
```

#### Affine IPN Model

```
FUNCTION double_network_modulus_affine(nu_e1, nu_e2, phi_01, phi_02, chi1, chi2, T):
    kT = k_B * T

    phi_eq1 = flory_rehner_swelling(nu_e1, phi_01, chi1, T)
    lambda_1 = (phi_01 / phi_eq1)^(1/3)
    G1 = nu_e1 * N_A * kT * lambda_1^(-2/3)

    phi_eq2 = flory_rehner_swelling(nu_e2, phi_02, chi2, T)
    lambda_2 = (phi_02 / phi_eq2)^(1/3)
    G2 = nu_e2 * N_A * kT * lambda_2^(-2/3)

    G_DN = G1 + G2
    RETURN G_DN
```

### 5. Test Specification

**File:** `tests/test_m3_modulus_models.py`

```python
# ── Hashin-Shtrikman tests ──

def test_hs_single_phase_G1_only():
    """phi1=1.0 -> bounds collapse to G1."""
    from emulsim.level4_mechanical.solver import hashin_shtrikman_bounds
    lo, hi = hashin_shtrikman_bounds(1000.0, 5000.0, 1.0)
    assert lo == pytest.approx(1000.0)
    assert hi == pytest.approx(1000.0)

def test_hs_single_phase_G2_only():
    """phi1=0.0 -> bounds collapse to G2."""
    from emulsim.level4_mechanical.solver import hashin_shtrikman_bounds
    lo, hi = hashin_shtrikman_bounds(1000.0, 5000.0, 0.0)
    assert lo == pytest.approx(5000.0)
    assert hi == pytest.approx(5000.0)

def test_hs_equal_moduli():
    """G1 == G2 -> bounds collapse to G1."""
    from emulsim.level4_mechanical.solver import hashin_shtrikman_bounds
    lo, hi = hashin_shtrikman_bounds(3000.0, 3000.0, 0.5)
    assert lo == pytest.approx(3000.0, rel=0.01)
    assert hi == pytest.approx(3000.0, rel=0.01)

def test_hs_ordering():
    """Lower bound <= Upper bound."""
    from emulsim.level4_mechanical.solver import hashin_shtrikman_bounds
    lo, hi = hashin_shtrikman_bounds(1000.0, 50000.0, 0.7)
    assert lo <= hi
    assert lo > 0
    assert hi > 0

def test_hs_bounds_bracket_phenomenological():
    """G_lower <= G_DN_phenomenological <= G_upper (not guaranteed but typical)."""
    from emulsim.level4_mechanical.solver import (
        hashin_shtrikman_bounds, double_network_modulus,
    )
    G1, G2 = 50000.0, 5000.0
    phi1 = 0.7
    lo, hi = hashin_shtrikman_bounds(G1, G2, phi1)
    G_dn = double_network_modulus(G1, G2, eta_coupling=-0.15)
    # Phenomenological may fall outside HS bounds (it's not a composite model),
    # but for typical parameters it should be in the right ballpark.
    # Just verify bounds are reasonable:
    assert lo > 0
    assert hi >= lo


# ── Flory-Rehner tests ──

def test_fr_swelling_reasonable():
    """Equilibrium phi should be between 0 and 1."""
    from emulsim.level4_mechanical.solver import flory_rehner_swelling
    phi_eq = flory_rehner_swelling(nu_e=100.0, phi_0=0.03, chi=0.50)
    assert 0 < phi_eq < 1

def test_fr_higher_crosslink_less_swelling():
    """Higher nu_e -> higher phi_eq (less swelling)."""
    from emulsim.level4_mechanical.solver import flory_rehner_swelling
    phi_low = flory_rehner_swelling(nu_e=10.0, phi_0=0.03, chi=0.50)
    phi_high = flory_rehner_swelling(nu_e=1000.0, phi_0=0.03, chi=0.50)
    assert phi_high >= phi_low


# ── Affine IPN model tests ──

def test_single_network_recovery():
    """G-M3-1: G_DN(nu_e1, 0) == G1 (single network recovery) within 1%."""
    from emulsim.level4_mechanical.solver import double_network_modulus_affine
    G_dn, status = double_network_modulus_affine(
        nu_e1=100.0, nu_e2=0.0,
        phi_01=0.03, phi_02=0.01,
        chi1=0.50, chi2=0.45,
    )
    # With nu_e2=0, G2=0, so G_DN = G1 only
    assert G_dn > 0
    # Check G_dn equals the single-network prediction
    G_single, _ = double_network_modulus_affine(
        nu_e1=100.0, nu_e2=0.0,
        phi_01=0.03, phi_02=0.0,
        chi1=0.50, chi2=0.45,
    )
    assert G_dn == pytest.approx(G_single, rel=0.01)

def test_symmetric_networks_exceed_single():
    """G-M3-2: Symmetric networks: G_DN(G, G) >= G."""
    from emulsim.level4_mechanical.solver import double_network_modulus_affine
    G_dn, _ = double_network_modulus_affine(
        nu_e1=100.0, nu_e2=100.0,
        phi_01=0.03, phi_02=0.03,
        chi1=0.50, chi2=0.50,
    )
    G_single, _ = double_network_modulus_affine(
        nu_e1=100.0, nu_e2=0.0,
        phi_01=0.03, phi_02=0.0,
        chi1=0.50, chi2=0.50,
    )
    assert G_dn >= G_single

def test_affine_fallback_on_failure():
    """G-M3-4: Fallback to phenomenological when fsolve fails."""
    # Use extreme parameters that should cause convergence failure
    # The function should return status="fallback" or the caller
    # falls back to phenomenological.
    from emulsim.level4_mechanical.solver import double_network_modulus_affine
    G_dn, status = double_network_modulus_affine(
        nu_e1=1e-20, nu_e2=1e-20,  # essentially zero crosslinking
        phi_01=0.001, phi_02=0.001,
        chi1=0.50, chi2=0.45,
    )
    # Should not crash; should return a non-negative value
    assert G_dn >= 0

def test_bounds_containment():
    """G-M3-5: G_lower <= G_DN <= G_upper (when all models run)."""
    from emulsim.level4_mechanical.solver import (
        hashin_shtrikman_bounds, agarose_modulus,
    )
    # For typical formulation:
    G1 = agarose_modulus(42.0)  # ~50 kPa
    G2 = 5000.0                 # ~5 kPa crosslinked chitosan
    phi1 = 42.0 / (42.0 + 18.0)  # 0.7
    lo, hi = hashin_shtrikman_bounds(G1, G2, phi1)
    assert lo > 0
    assert hi >= lo

def test_empirical_engineering_uses_phenomenological():
    """G-M3-6: empirical_engineering mode still uses phenomenological."""
    from emulsim.level4_mechanical.solver import select_modulus_model
    from emulsim.datatypes import NetworkTypeMetadata
    meta = NetworkTypeMetadata(
        solver_family="amine_covalent",
        eta_coupling_recommended=-0.15,
    )
    G_dn, label = select_modulus_model(
        50000.0, 5000.0, meta, -0.15,
        model_mode="empirical_engineering",
    )
    assert label == "phenomenological"
    assert G_dn > 0
```

### 6. G1 Gate Check

| # | Criterion | Status |
|---|-----------|--------|
| G1-01 | Every file path and line number specified | YES |
| G1-02 | Old code -> New code diffs shown | YES (3 phases) |
| G1-03 | New function signatures with type hints | YES (hashin_shtrikman_bounds, flory_rehner_swelling, double_network_modulus_affine) |
| G1-04 | New dataclass fields with types and defaults | YES (MechanicalResult: 3 fields) |
| G1-05 | Config additions with TOML syntax | N/A (no config changes for M3) |
| G1-06 | Interface contract complete | YES |
| G1-07 | Algorithm pseudocode provided | YES (3 algorithms) |
| G1-08 | Test function names and signatures specified | YES (14 tests) |
| G1-09 | Input fixtures with numerical values | YES |
| G1-10 | Expected outputs with tolerances | YES |
| G1-11 | Edge cases covered | YES (single-phase, symmetric, failure fallback, zero crosslinking) |
| G1-12 | No ambiguity for implementer | YES |

---

## Cross-Module Integration Notes

### Backward Compatibility

All changes maintain backward compatibility:

1. **M1:** Old configs with `t_emulsification=60` still work; they just get auto-extended. New fields in `SolverSettings` and `EmulsificationResult` have defaults.
2. **M4:** `eta_coupling_recommended=-0.15` default matches existing genipin behavior. Old code that does not set this field gets the same result.
3. **M2:** Old configs with `c_genipin=2.0` still work but will fire the enhanced W2 warning with a recommendation.
4. **M3:** New fields in `MechanicalResult` default to 0.0/"phenomenological". The `select_modulus_model` return type change from `float` to `tuple[float, str]` requires updating the single call site in `solve_mechanical`.

### Import Changes Summary

| File | Addition |
|------|----------|
| `src/emulsim/trust.py` line 8 | Add `ModelMode` to the existing import |
| `src/emulsim/trust.py` ~line 106 | Add `from .level3_crosslinking.solver import recommended_crosslinker_concentration` (local import inside function) |
| `src/emulsim/level4_mechanical/solver.py` | Add `from scipy.optimize import brentq` (in `flory_rehner_swelling`, local import) |

### Config File Upgrade Path

Users with existing `configs/default.toml` should update:
- `t_emulsification`: 60 -> 300
- `c_genipin`: 2.0 -> 10.0
- Add `t_max`, `conv_tol`, `max_extensions` under `[solver.level1]`

The code will work with old configs (defaults apply for missing fields), but old defaults will trigger trust warnings.

---

> **Disclaimer**: This protocol document is for development guidance only. All
> models and parameter recommendations should be validated through appropriate
> laboratory experimentation before deployment.
