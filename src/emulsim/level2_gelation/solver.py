"""Cahn-Hilliard phase-field solver for gelation and pore formation.

Provides two solver classes:
- ``CahnHilliardSolver``: legacy 1D radial solver (concentric shells only).
- ``CahnHilliard2DSolver``: 2D Cartesian solver with face-centred variable
  mobility.  This is the recommended default as it resolves true pore
  topology (bicontinuous structures from spinodal decomposition).

Both use semi-implicit Eyre time stepping coupled with Avrami gelation
kinetics for mobility arrest.
"""

from __future__ import annotations
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import factorized, spsolve

from ..datatypes import GelationResult, MaterialProperties, SimulationParameters
from .spatial import (
    create_radial_grid, build_laplacian_matrix,
    create_2d_grid, build_laplacian_2d, build_mobility_laplacian_2d,
)
from .free_energy import flory_huggins_mu, flory_huggins_d2f, contractive_constant
from .gelation import avrami_gelation, gelation_rate_constant, cooling_temperature, mobility
from .pore_analysis import (
    characteristic_wavelength, chord_length_distribution, compute_porosity,
    characteristic_wavelength_2d, chord_length_distribution_2d, compute_porosity_2d,
)


def _find_spinodal_temperature(A_chi: float, B_chi: float, N_p: float,
                               phi_0: float, v0: float = 1.8e-29) -> float:
    """Find temperature at which the system enters the spinodal region.

    Returns T_spinodal such that f''(phi_0, T_spinodal) = 0.
    Returns 0.0 if no spinodal crossing is found.
    """
    # f''(phi) = (kT/v0) * [1/(Np*phi) + 1/(1-phi) - 2*chi(T)]
    # f'' = 0 when chi(T) = 0.5 * [1/(Np*phi) + 1/(1-phi)]
    # chi(T) = A/T + B, so T = A / (chi_c - B)
    phi = np.clip(phi_0, 0.01, 0.99)
    chi_c_local = 0.5 * (1.0 / (N_p * phi) + 1.0 / (1.0 - phi))
    denom = chi_c_local - B_chi
    if denom <= 0:
        return 0.0
    return A_chi / denom


class CahnHilliardSolver:
    """1D radial Cahn-Hilliard solver with gelation arrest.

    Simulates spinodal decomposition of the polymer solution during
    cooling inside a single microsphere.
    """

    def __init__(self, N_r: int = 1000, dt_initial: float = 1e-4,
                 dt_max: float = 1.0, arrest_exponent: float = 2.5):
        self.N_r = N_r
        self.dt_initial = dt_initial
        self.dt_max = dt_max
        self.arrest_exponent = arrest_exponent

    def solve(self, params: SimulationParameters, props: MaterialProperties,
              R_droplet: float = 1.0e-6, rng_seed: int = 42) -> GelationResult:
        """Run the phase-field simulation for a single microsphere.

        Parameters
        ----------
        params : SimulationParameters
        props : MaterialProperties
        R_droplet : float
            Microsphere radius [m] (from Level 1 d32/2).
        rng_seed : int
            Random seed for initial perturbation.
        """
        # -- Setup ----------------------------------------------------------
        r, dr = create_radial_grid(R_droplet, self.N_r)
        L = build_laplacian_matrix(r, dr)
        I_mat = sparse.eye(self.N_r, format='csc')
        L_csc = L.tocsc()

        # Precompute L^2 once
        L_sq = (L @ L).tocsc()

        # Physical parameters
        T_init = params.formulation.T_oil
        T_ambient = 293.15  # 20 C room temperature
        cooling_rate = params.formulation.cooling_rate  # K/s
        T_gel = props.T_gel
        kappa = props.kappa_CH
        M_0 = props.M_0
        beta = self.arrest_exponent

        # Flory-Huggins parameters
        N_p = 100.0  # effective degree of polymerization
        v0 = 1.8e-29  # lattice site volume [m^3]
        A_chi, B_chi = props.chi_T_coeffs

        # Initial condition: uniform + small noise
        phi_0 = (params.formulation.c_agarose + params.formulation.c_chitosan) / 1400.0
        phi_0 = np.clip(phi_0, 0.02, 0.15)

        rng = np.random.default_rng(rng_seed)
        phi = phi_0 + 0.005 * rng.standard_normal(self.N_r)
        phi = np.clip(phi, 0.001, 0.999)

        # -- Determine time windows -----------------------------------------
        # Find when spinodal decomposition starts
        T_spinodal = _find_spinodal_temperature(A_chi, B_chi, N_p, phi_0, v0)

        # Determine total simulation time
        if cooling_rate > 0:
            t_final = (T_init - T_ambient) / cooling_rate * 1.5  # 50% extra
        else:
            t_final = 600.0

        # Fast-forward to just before spinodal onset (nothing happens before)
        t_start = 0.0
        if T_spinodal > T_ambient and cooling_rate > 0:
            t_spinodal = (T_init - T_spinodal) / cooling_rate
            t_start = max(0.0, t_spinodal - 5.0)  # start 5s before spinodal

        # -- Time stepping --------------------------------------------------
        t = t_start
        dt = self.dt_initial
        alpha = 0.0
        t_gel_cross = None

        # Check if gelation already started during fast-forward
        T_at_start = cooling_temperature(t_start, T_init, T_ambient, cooling_rate)
        if T_at_start < T_gel:
            t_gel_cross = max(0.0, (T_init - T_gel) / cooling_rate) if cooling_rate > 0 else 0.0
            t_cool = t_start - t_gel_cross
            if t_cool > 0:
                k_gel = gelation_rate_constant(T_at_start, T_gel, props.k_gel_0)
                alpha = avrami_gelation(t_cool, k_gel, props.n_avrami)

        # Storage for snapshots
        T_history = []
        alpha_history = []
        snap_times = np.linspace(t_start, t_final, 20)
        snap_idx = 0
        phi_snapshots = []

        step = 0
        max_steps = 200000

        # Cache for LHS factorization
        cached_factor = None
        cached_dm = -1.0
        cached_C = -1.0

        while t < t_final and step < max_steps:
            # Current temperature
            T = cooling_temperature(t, T_init, T_ambient, cooling_rate)

            # Chi parameter
            chi = A_chi / max(T, 100.0) + B_chi

            # Gelation kinetics
            if T < T_gel:
                if t_gel_cross is None:
                    t_gel_cross = t
                t_cool = t - t_gel_cross
                k_gel = gelation_rate_constant(T, T_gel, props.k_gel_0)
                alpha = avrami_gelation(t_cool, k_gel, props.n_avrami)

            # Mobility field
            M_field = mobility(phi, alpha, M_0, beta)
            M_avg = np.mean(M_field) if np.any(M_field > 0) else 0.0

            # If gelation is essentially complete, stop
            if alpha > 0.999:
                break

            # If mobility is negligible, stop
            if M_avg < 1e-30:
                break

            # Contractive constant for Eyre splitting
            C_contract = contractive_constant((0.01, 0.30), T, chi, N_p, v0)

            # Compute dt * M_avg product
            dm = dt * M_avg

            # Explicit chemical potential
            mu_explicit = flory_huggins_mu(phi, T, chi, N_p, v0) + C_contract * phi

            # RHS: phi + dm * L * mu_e
            rhs = phi + dm * (L @ mu_explicit)

            # LHS matrix: I + dm*C*L + dm*kappa*L^2
            # Use cached factorization if dm and C haven't changed significantly
            need_refactor = True
            if cached_dm > 0:
                rel_dm = abs(dm - cached_dm) / (cached_dm + 1e-50)
                rel_C = abs(C_contract - cached_C) / (cached_C + 1e-50)
                if rel_dm < 0.01 and rel_C < 0.01:
                    need_refactor = False

            if need_refactor:
                A_mat = I_mat + dm * C_contract * L_csc + dm * kappa * L_sq
                try:
                    cached_factor = factorized(A_mat)
                    cached_dm = dm
                    cached_C = C_contract
                except Exception:
                    cached_factor = None

            # Solve
            try:
                if cached_factor is not None:
                    phi_new = cached_factor(rhs)
                else:
                    A_mat = I_mat + dm * C_contract * L_csc + dm * kappa * L_sq
                    phi_new = spsolve(A_mat, rhs)
            except Exception:
                dt *= 0.25
                cached_factor = None
                cached_dm = -1.0
                if dt < 1e-10:
                    break
                continue

            # Clip to physical bounds
            phi_new = np.clip(phi_new, 1e-6, 1.0 - 1e-6)

            # Adaptive time stepping: target change_target per step
            change = np.max(np.abs(phi_new - phi))
            change_target = 0.02
            old_dt = dt
            if change > 0 and change < 0.2:
                # Scale dt to target change_target
                ratio = change_target / change
                ratio = np.clip(ratio, 0.25, 4.0)
                dt = np.clip(dt * ratio, 1e-10, self.dt_max)
            elif change >= 0.2:
                # Reject step: reduce dt significantly
                dt = max(dt * 0.1, 1e-10)
                cached_factor = None
                cached_dm = -1.0
                continue  # retry with smaller dt

            if abs(dt - old_dt) / (old_dt + 1e-50) > 0.05:
                cached_factor = None
                cached_dm = -1.0

            phi = phi_new
            t += dt
            step += 1

            # Record snapshot
            if snap_idx < len(snap_times) and t >= snap_times[snap_idx]:
                T_history.append((t, T))
                alpha_history.append(alpha)
                phi_snapshots.append(phi.copy())
                snap_idx += 1

        # -- Pore analysis --------------------------------------------------
        char_wl = characteristic_wavelength(phi, r)
        chords = chord_length_distribution(phi, r)
        pore_mean = np.mean(chords) if len(chords) > 0 else 0.0
        pore_std = np.std(chords) if len(chords) > 1 else 0.0
        porosity = compute_porosity(phi, r)

        T_hist_arr = np.array([T for _, T in T_history]) if T_history else np.array([0.0])

        return GelationResult(
            r_grid=r,
            phi_field=phi,
            pore_size_mean=pore_mean,
            pore_size_std=pore_std,
            pore_size_distribution=chords,
            porosity=porosity,
            alpha_final=alpha,
            char_wavelength=char_wl,
            T_history=T_hist_arr,
            phi_snapshots=np.array(phi_snapshots) if phi_snapshots else None,
        )


class CahnHilliard2DSolver:
    """2D Cartesian Cahn-Hilliard solver with face-centred variable mobility.

    Simulates spinodal decomposition on a square cross-section through the
    microsphere, resolving true 2D pore topology (bicontinuous structures).

    Key improvements over the 1D radial solver:
    - Captures pore connectivity and bicontinuous morphology.
    - Uses face-centred mobility averaging: M_{i+1/2,j} = 0.5*(M_{i,j} + M_{i+1,j})
      instead of collapsing to a single M_avg.
    """

    def __init__(self, N_grid: int = 128, dt_initial: float = 1e-4,
                 dt_max: float = 1.0, arrest_exponent: float = 2.5):
        self.N_grid = N_grid
        self.dt_initial = dt_initial
        self.dt_max = dt_max
        self.arrest_exponent = arrest_exponent

    def solve(self, params: SimulationParameters, props: MaterialProperties,
              R_droplet: float = 1.0e-6, rng_seed: int = 42) -> GelationResult:
        """Run the 2D phase-field simulation for a single microsphere.

        Parameters
        ----------
        params : SimulationParameters
        props : MaterialProperties
        R_droplet : float
            Microsphere radius [m] (from Level 1 d32/2).
        rng_seed : int
            Random seed for initial perturbation.
        """
        N = self.N_grid
        total = N * N

        # -- Setup ----------------------------------------------------------
        L_domain = 2.0 * R_droplet
        x, y, h = create_2d_grid(L_domain, N)
        L_std = build_laplacian_2d(N, h)      # standard (constant-coeff) Laplacian
        L_std_csc = L_std.tocsc()
        I_mat = sparse.eye(total, format='csc')

        # Physical parameters
        T_init = params.formulation.T_oil
        T_ambient = 293.15
        cooling_rate = params.formulation.cooling_rate
        T_gel = props.T_gel
        kappa = props.kappa_CH
        M_0 = props.M_0
        beta = self.arrest_exponent

        # Flory-Huggins parameters
        N_p = 100.0
        v0 = 1.8e-29
        A_chi, B_chi = props.chi_T_coeffs

        # Initial condition: uniform + small noise (on 2D grid)
        phi_0 = (params.formulation.c_agarose + params.formulation.c_chitosan) / 1400.0
        phi_0 = float(np.clip(phi_0, 0.02, 0.15))

        rng = np.random.default_rng(rng_seed)
        phi_2d = phi_0 + 0.005 * rng.standard_normal((N, N))
        phi_2d = np.clip(phi_2d, 0.001, 0.999)

        # Flatten for linear algebra (row-major)
        phi = phi_2d.ravel()

        # -- Determine time windows -----------------------------------------
        T_spinodal = _find_spinodal_temperature(A_chi, B_chi, N_p, phi_0, v0)

        if cooling_rate > 0:
            t_final = (T_init - T_ambient) / cooling_rate * 1.5
        else:
            t_final = 600.0

        t_start = 0.0
        if T_spinodal > T_ambient and cooling_rate > 0:
            t_spinodal_time = (T_init - T_spinodal) / cooling_rate
            t_start = max(0.0, t_spinodal_time - 5.0)

        # -- Time stepping --------------------------------------------------
        t = t_start
        dt = self.dt_initial
        alpha = 0.0
        t_gel_cross = None

        # Check if gelation already started during fast-forward
        T_at_start = cooling_temperature(t_start, T_init, T_ambient, cooling_rate)
        if T_at_start < T_gel:
            t_gel_cross = max(0.0, (T_init - T_gel) / cooling_rate) if cooling_rate > 0 else 0.0
            t_cool = t_start - t_gel_cross
            if t_cool > 0:
                k_gel = gelation_rate_constant(T_at_start, T_gel, props.k_gel_0)
                alpha = avrami_gelation(t_cool, k_gel, props.n_avrami)

        # Snapshot storage
        T_history = []
        alpha_history = []
        snap_times = np.linspace(t_start, t_final, 20)
        snap_idx = 0
        phi_snapshots = []

        step = 0
        max_steps = 200000

        while t < t_final and step < max_steps:
            # Current temperature
            T = cooling_temperature(t, T_init, T_ambient, cooling_rate)

            # Chi parameter
            chi = A_chi / max(T, 100.0) + B_chi

            # Gelation kinetics
            if T < T_gel:
                if t_gel_cross is None:
                    t_gel_cross = t
                t_cool = t - t_gel_cross
                k_gel = gelation_rate_constant(T, T_gel, props.k_gel_0)
                alpha = avrami_gelation(t_cool, k_gel, props.n_avrami)

            # If gelation is essentially complete, stop
            if alpha > 0.999:
                break

            # Mobility field: compute at every grid point
            M_flat = mobility(phi, alpha, M_0, beta)

            # If mobility is negligible, stop
            M_avg = np.mean(M_flat)
            if M_avg < 1e-30:
                break

            M_2d = M_flat.reshape(N, N)

            # Build mobility-weighted Laplacian: L_M = div(M * grad(.))
            L_M = build_mobility_laplacian_2d(M_2d, N, h)
            L_M_csc = L_M.tocsc()

            # Contractive constant for Eyre splitting
            C_contract = contractive_constant((0.01, 0.30), T, chi, N_p, v0)

            # Explicit chemical potential: mu_e = f'(phi) + C*phi
            mu_explicit = flory_huggins_mu(phi, T, chi, N_p, v0) + C_contract * phi

            # RHS: phi + dt * L_M @ mu_e
            rhs = phi + dt * (L_M @ mu_explicit)

            # LHS: I + dt*C*L_M + dt*kappa*(L_M @ L_std)
            # The implicit part handles -C*phi - kappa*laplacian(phi)
            LM_L = L_M_csc @ L_std_csc
            A_mat = I_mat + dt * C_contract * L_M_csc + dt * kappa * LM_L

            # Solve
            try:
                phi_new = spsolve(A_mat, rhs)
            except Exception:
                dt *= 0.25
                if dt < 1e-10:
                    break
                continue

            # Clip to physical bounds
            phi_new = np.clip(phi_new, 1e-6, 1.0 - 1e-6)

            # Adaptive time stepping
            change = np.max(np.abs(phi_new - phi))
            change_target = 0.02
            if change > 0 and change < 0.2:
                ratio = change_target / change
                ratio = np.clip(ratio, 0.25, 4.0)
                dt = np.clip(dt * ratio, 1e-10, self.dt_max)
            elif change >= 0.2:
                dt = max(dt * 0.1, 1e-10)
                continue  # retry with smaller dt

            phi = phi_new
            t += dt
            step += 1

            # Record snapshot
            if snap_idx < len(snap_times) and t >= snap_times[snap_idx]:
                T_history.append((t, T))
                alpha_history.append(alpha)
                phi_snapshots.append(phi.copy().reshape(N, N))
                snap_idx += 1

        # -- Pore analysis --------------------------------------------------
        phi_2d_final = phi.reshape(N, N)

        char_wl = characteristic_wavelength_2d(phi_2d_final, h)
        chords = chord_length_distribution_2d(phi_2d_final, h)
        pore_mean = float(np.mean(chords)) if len(chords) > 0 else 0.0
        pore_std = float(np.std(chords)) if len(chords) > 1 else 0.0
        porosity = compute_porosity_2d(phi_2d_final)

        T_hist_arr = np.array([T for _, T in T_history]) if T_history else np.array([0.0])

        return GelationResult(
            r_grid=x,  # 1D coordinate array
            phi_field=phi_2d_final,
            pore_size_mean=pore_mean,
            pore_size_std=pore_std,
            pore_size_distribution=chords,
            porosity=porosity,
            alpha_final=alpha,
            char_wavelength=char_wl,
            T_history=T_hist_arr,
            phi_snapshots=np.array(phi_snapshots) if phi_snapshots else None,
            L_domain=L_domain,
            grid_spacing=h,
        )


def solve_gelation(params: SimulationParameters, props: MaterialProperties,
                   R_droplet: float = 1.0e-6, use_2d: bool = True) -> GelationResult:
    """Convenience function for Level 2 gelation simulation.

    Parameters
    ----------
    params : SimulationParameters
    props : MaterialProperties
    R_droplet : float
        Microsphere radius [m].
    use_2d : bool
        If True (default), use the 2D Cartesian solver with face-centred
        variable mobility.  Set False to fall back to the legacy 1D radial
        solver.
    """
    if use_2d:
        solver = CahnHilliard2DSolver(
            N_grid=params.solver.l2_n_grid,
            dt_initial=params.solver.l2_dt_initial,
            dt_max=params.solver.l2_dt_max,
            arrest_exponent=params.solver.l2_arrest_exponent,
        )
    else:
        solver = CahnHilliardSolver(
            N_r=params.solver.l2_n_r,
            dt_initial=params.solver.l2_dt_initial,
            dt_max=params.solver.l2_dt_max,
            arrest_exponent=params.solver.l2_arrest_exponent,
        )
    return solver.solve(params, props, R_droplet)
