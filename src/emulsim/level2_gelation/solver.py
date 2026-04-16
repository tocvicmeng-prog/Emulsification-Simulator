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
import logging
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import factorized, spsolve, splu

from ..datatypes import GelationResult, GelationTimingResult, MaterialProperties, ModelEvidenceTier, ModelManifest, SimulationParameters
from .spatial import (
    create_radial_grid, build_laplacian_matrix,
    create_2d_grid, build_laplacian_2d, build_mobility_laplacian_2d,
)
from .free_energy import flory_huggins_mu, flory_huggins_d2f, contractive_constant
from .gelation import avrami_gelation, gelation_rate_constant, cooling_temperature, mobility
from .pore_analysis import (
    characteristic_wavelength, chord_length_distribution, compute_porosity,
    characteristic_wavelength_2d, chord_length_distribution_2d, compute_porosity_2d,
    morphology_descriptors,
)

logger = logging.getLogger(__name__)


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

        # Flory-Huggins parameters (same as 2D solver)
        N_p = 10.0   # effective DP for agarose helix aggregates
        v0 = 1.8e-29
        A_chi, B_chi = props.chi_T_coeffs

        # Effective hydrated volume fraction
        hydration_factor = 6.0
        phi_0_dry = (params.formulation.c_agarose + params.formulation.c_chitosan) / 1400.0
        phi_0 = np.clip(phi_0_dry * hydration_factor, 0.05, 0.45)

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
            if step == 0:
                # Precompute contractive constant once at the deepest expected quench
                chi_deep = A_chi / max(T_ambient, 200.0) + B_chi
                C_contract = contractive_constant(
                    (max(0.03, phi_0 * 0.3), min(0.97, phi_0 * 3.0)),
                    T_ambient, chi_deep, N_p, v0,
                )

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

            if not np.all(np.isfinite(phi_new)):
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

        # Morphology descriptors (1D field -> defaults)
        morph = morphology_descriptors(phi, dr)

        T_hist_arr = np.array([T for _, T in T_history]) if T_history else np.array([0.0])

        model_manifest = ModelManifest(
            model_name="L2.Pore.CahnHilliard2D",
            evidence_tier=ModelEvidenceTier.SEMI_QUANTITATIVE,
            assumptions=["Cahn-Hilliard free energy", "mobility arrest", "2D approximation"],
        )
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
            bicontinuous_score=morph['bicontinuous_score'],
            anisotropy=morph['anisotropy'],
            connectivity=morph['connectivity'],
            chord_skewness=morph['skewness'],
            model_tier="mechanistic",
            model_manifest=model_manifest,
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

        # -- Physical parameters (needed for adaptive domain sizing) ----------
        T_init = params.formulation.T_oil
        T_ambient = 293.15
        cooling_rate = params.formulation.cooling_rate
        T_gel = props.T_gel
        kappa = props.kappa_CH
        M_0 = props.M_0
        beta = self.arrest_exponent

        # Flory-Huggins parameters
        # N_p ~ 10: effective degree of polymerization for agarose helix aggregates
        # (not individual monomers — the structural unit is a hydrated fiber bundle)
        N_p = 10.0
        v0 = 1.8e-29
        A_chi, B_chi = props.chi_T_coeffs

        # Initial condition: effective hydrated volume fraction
        # Agarose helices bind ~5-10x their dry weight in water during gelation,
        # so the effective volume fraction is much larger than c_polymer/rho_polymer.
        # This places phi_0 in the bicontinuous spinodal regime (near phi_c ~ 0.24).
        hydration_factor = 6.0  # typical for agarose helix bundles
        phi_0_dry = (params.formulation.c_agarose + params.formulation.c_chitosan) / 1400.0
        phi_0 = float(np.clip(phi_0_dry * hydration_factor, 0.05, 0.45))

        # -- Setup ----------------------------------------------------------
        # Adaptive domain size: fit at least 5 Cahn-Hilliard spinodal wavelengths.
        # The CH wavelength lambda_CH = 2*pi*sqrt(2*kappa / f'') where f'' is the
        # Flory-Huggins second derivative at the mean composition.
        # This solver is a local-patch (RVE) representation, not the whole droplet.
        _chi_ref = A_chi / max(T_init, 100.0) + B_chi
        _f_pp = max(abs(1.0 / max(phi_0, 0.01) + 1.0 / max(1.0 - phi_0, 0.01) - 2.0 * _chi_ref), 1.0)
        lambda_CH = 2.0 * np.pi * np.sqrt(2.0 * kappa / _f_pp)
        # Domain must fit at least 5 wavelengths, with a minimum of 0.5 µm
        L_domain = min(2.0 * R_droplet, max(5.0 * lambda_CH, 0.5e-6))
        # Hard ceiling for computational feasibility
        L_domain = min(L_domain, 3.0e-6)
        if L_domain < 2.0 * R_droplet:
            logger.info("L2: domain %.1f nm < droplet %.1f nm (local-patch RVE mode)",
                        L_domain * 1e9, 2 * R_droplet * 1e9)
        if R_droplet < 5 * lambda_CH:
            logger.warning("L2: droplet R=%.0f nm is within 5x of spinodal wavelength "
                           "%.0f nm — confinement effects may be significant",
                           R_droplet * 1e9, lambda_CH * 1e9)
        x, y, h = create_2d_grid(L_domain, N)
        L_std = build_laplacian_2d(N, h)      # standard (constant-coeff) Laplacian
        L_std_csc = L_std.tocsc()
        I_mat = sparse.eye(total, format='csc')

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

        # -- Pre-loop: contractive constant (computed once at deepest quench) --
        chi_deep = A_chi / max(T_ambient, 200.0) + B_chi
        C_contract = contractive_constant(
            (max(0.03, phi_0 * 0.3), min(0.97, phi_0 * 3.0)),
            T_ambient, chi_deep, N_p, v0,
        )

        # -- Pre-loop: initial LU factorization ---------------------------------
        # Build mobility for the initial phi field
        M_flat_init = mobility(phi, 0.0, M_0, beta)
        M_2d_init = M_flat_init.reshape(N, N)
        L_M = build_mobility_laplacian_2d(M_2d_init, N, h)
        L_M_csc = L_M.tocsc()
        LM_L = L_M_csc @ L_std_csc
        A_mat = I_mat + dt * C_contract * L_M_csc + dt * kappa * LM_L
        try:
            A_factor = splu(A_mat.tocsc())
        except Exception:
            A_factor = None
        M_cached = M_2d_init.copy()

        REBUILD_INTERVAL = 10   # rebuild operator every N steps at most
        rebuild_counter = 0     # steps since last rebuild

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

            # Explicit chemical potential: mu_e = f'(phi) + C*phi
            mu_explicit = flory_huggins_mu(phi, T, chi, N_p, v0) + C_contract * phi

            # RHS: phi + dt * L_M @ mu_e   (uses current L_M, cached between rebuilds)
            rhs = phi + dt * (L_M @ mu_explicit)

            # -- Decide whether to rebuild the operator -----------------------
            rebuild_counter += 1
            need_rebuild = (rebuild_counter >= REBUILD_INTERVAL) or (A_factor is None)

            # Also rebuild if mobility has drifted significantly from cached value
            if not need_rebuild:
                M_avg_cached = float(np.mean(M_cached)) + 1e-30
                if np.max(np.abs(M_2d - M_cached)) / M_avg_cached > 0.05:
                    need_rebuild = True

            if need_rebuild:
                L_M = build_mobility_laplacian_2d(M_2d, N, h)
                L_M_csc = L_M.tocsc()
                LM_L = L_M_csc @ L_std_csc
                A_mat = I_mat + dt * C_contract * L_M_csc + dt * kappa * LM_L
                try:
                    A_factor = splu(A_mat.tocsc())
                except Exception:
                    A_factor = None
                M_cached = M_2d.copy()
                rebuild_counter = 0

            # Solve using cached LU factorization (or fall back to spsolve)
            try:
                if A_factor is not None:
                    phi_new = A_factor.solve(rhs)
                else:
                    phi_new = spsolve(A_mat, rhs)
            except Exception:
                dt *= 0.25
                # Force rebuild on next step with the new dt
                A_factor = None
                rebuild_counter = REBUILD_INTERVAL
                if dt < 1e-10:
                    break
                continue

            if not np.all(np.isfinite(phi_new)):
                dt *= 0.25
                # Force rebuild on next step with the new dt
                A_factor = None
                rebuild_counter = REBUILD_INTERVAL
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
                new_dt = np.clip(dt * ratio, 1e-10, self.dt_max)
                # If dt changes significantly, force operator rebuild next step
                if abs(new_dt - dt) / (dt + 1e-50) > 0.05:
                    A_factor = None
                    rebuild_counter = REBUILD_INTERVAL
                dt = new_dt
            elif change >= 0.2:
                dt = max(dt * 0.1, 1e-10)
                A_factor = None
                rebuild_counter = REBUILD_INTERVAL
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

        # Morphology descriptors
        morph = morphology_descriptors(phi_2d_final, h)

        T_hist_arr = np.array([T for _, T in T_history]) if T_history else np.array([0.0])

        model_manifest = ModelManifest(
            model_name="L2.Pore.CahnHilliard2D",
            evidence_tier=ModelEvidenceTier.SEMI_QUANTITATIVE,
            assumptions=["Cahn-Hilliard free energy", "mobility arrest", "2D approximation"],
        )
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
            bicontinuous_score=morph['bicontinuous_score'],
            anisotropy=morph['anisotropy'],
            connectivity=morph['connectivity'],
            chord_skewness=morph['skewness'],
            model_tier="mechanistic",
            model_manifest=model_manifest,
        )


def solve_gelation_timing(
    params: SimulationParameters,
    props: MaterialProperties,
    R_droplet: float = 50e-6,
) -> GelationTimingResult:
    """Level 2a: Compute gelation timing and arrest for a droplet.

    Separates the WHEN (timing) from the WHAT (microstructure).
    """
    T_oil = params.formulation.T_oil
    T_gel = props.T_gel
    cool_rate = params.formulation.cooling_rate

    # Effective cooling rate (size-dependent for small droplets)
    R_ref = 50e-6
    if R_droplet > 0 and R_droplet < R_ref:
        cool_rate_eff = cool_rate * (R_ref / R_droplet) ** 0.3
        cool_rate_eff = min(cool_rate_eff, cool_rate * 3.0)
    else:
        cool_rate_eff = cool_rate

    # Time to reach gelation temperature
    if T_oil > T_gel and cool_rate_eff > 0:
        t_gel_onset = (T_oil - T_gel) / cool_rate_eff
    else:
        t_gel_onset = 0.0  # already below T_gel

    # Temperature history (simple linear cooling for now)
    t_total = max(t_gel_onset * 2.0, 60.0)  # simulate past gel point
    t_arr = np.linspace(0, t_total, 200)
    T_arr = np.maximum(T_oil - cool_rate_eff * t_arr, 293.15)  # floor at 20C

    # Gelation fraction via simplified Avrami
    # alpha(t) = 1 - exp(-k_gel * (t - t_gel)^n) for t > t_gel
    k_gel = props.k_gel_0
    n_avrami = props.n_avrami
    alpha_arr = np.zeros_like(t_arr)
    mask = t_arr > t_gel_onset
    if mask.any():
        dt_gel = t_arr[mask] - t_gel_onset
        alpha_arr[mask] = 1.0 - np.exp(-k_gel * dt_gel ** n_avrami)

    alpha_final = float(alpha_arr[-1])

    # Mobility arrest factor: how much the gel network immobilises diffusion
    # ASSUMPTION: arrest ~ alpha^beta where beta is the arrest exponent
    beta = props.gel_arrest_exponent
    mobility_arrest_factor = max(1.0 - alpha_final ** beta, 0.01)

    return GelationTimingResult(
        T_history=T_arr,
        t_gel_onset=t_gel_onset,
        alpha_final=alpha_final,
        mobility_arrest_factor=mobility_arrest_factor,
        cooling_rate_effective=cool_rate_eff,
    )


def solve_gelation_empirical(params: SimulationParameters, props: MaterialProperties,
                             R_droplet: float = 1.0e-6,
                             timing: 'GelationTimingResult | None' = None) -> GelationResult:
    """Empirical pore-size model for agarose gels based on literature correlations.

    Uses the well-established relationships:
    - Pore diameter ~ c_agarose^(-0.5) (Pernodet et al. 1997, Narayanan et al. 2006)
    - Pore diameter decreases with faster cooling (Aymard et al. 2001)
    - Higher MW agarose → smaller pores (Zhao et al. 2020)

    This is the default L2 model. For mechanistic phase-field simulation,
    use solve_gelation(..., mode='ch_2d') instead.

    Node 8 (F8): When ``timing`` is provided, ``alpha_final`` is taken from
    ``timing.alpha_final`` (the Avrami output) instead of being hardcoded
    to 0.999. The hardcode was a v6.0 placeholder that decoupled the
    pore-size and gelation-completion outputs even though both are L2
    products. Backward-compat: ``timing=None`` reproduces the v6.0 behaviour.
    """
    c_agar = params.formulation.c_agarose  # kg/m³
    c_chit = params.formulation.c_chitosan
    cool_rate = params.formulation.cooling_rate  # K/s
    if cool_rate <= 0:
        cool_rate = 0.001  # fallback: very slow cooling

    # Empirical pore diameter [m] for agarose gels
    # d_pore = A * c^(-alpha) * (dT/dt)^(-beta)
    # Calibrated to: 2% agarose → ~300 nm, 4% → ~150 nm, 6% → ~80 nm
    # Cooling effect: 2°C/min → 1.3x, 10°C/min → 1.0x, 20°C/min → 0.8x
    c_pct = c_agar / 10.0  # convert to % w/v
    if c_pct < 0.5:
        c_pct = 0.5

    d_pore_base = 600e-9 * c_pct ** (-0.7)  # base pore size [m]

    # Cooling rate correction (normalised to 10°C/min = 0.167 K/s)
    cool_Cmin = cool_rate * 60.0  # K/s → °C/min
    cool_factor = (cool_Cmin / 10.0) ** (-0.2)  # mild effect
    cool_factor = np.clip(cool_factor, 0.5, 2.0)

    # Chitosan effect: presence of chitosan reduces pore size slightly
    # (phase separation between agarose and chitosan domains)
    chit_factor = 1.0 - 0.3 * (c_chit / (c_agar + c_chit + 1e-10))

    d_pore = d_pore_base * cool_factor * chit_factor

    # Confinement correction: spinodal wavelength truncated by sphere boundary
    # ASSUMPTION: max pore ~ 15% of bead diameter (Tanaka 2000, Physical Review E)
    d_pore_confinement_max = 0.15 * 2.0 * R_droplet
    if d_pore_confinement_max > 0:
        d_pore = min(d_pore, d_pore_confinement_max)

    # Size-dependent effective cooling: smaller droplets lose heat faster
    # (surface-to-volume ratio scales as 1/R), producing finer pores.
    # ASSUMPTION: mild power-law correction normalized to R_ref=50 um.
    R_ref = 50e-6  # [m] reference droplet radius
    if R_droplet > 0 and R_droplet < R_ref:
        size_cool_factor = (R_ref / R_droplet) ** 0.3  # faster effective cooling
        size_cool_factor = min(size_cool_factor, 3.0)   # cap at 3x
        d_pore *= size_cool_factor ** (-0.2)             # faster cooling → smaller pores

    d_pore = np.clip(d_pore, 10e-9, 2000e-9)

    # Porosity from Ogston model: porosity ≈ 1 - phi_fiber_effective
    rho_polymer = 1400.0
    phi_dry = (c_agar + c_chit) / rho_polymer
    porosity = 1.0 - phi_dry * 3.0  # hydrated fiber volume ~ 3x dry
    porosity = np.clip(porosity, 0.3, 0.98)

    # Pore size distribution (log-normal with CV ~ 30%)
    rng = np.random.default_rng(42)
    n_pores = 200
    pore_sizes = rng.lognormal(np.log(d_pore), 0.3, n_pores)

    # Build a synthetic 1D radial profile (for compatibility with L4)
    N_r = 100
    r = np.linspace(R_droplet * 0.005, R_droplet * 0.995, N_r)

    # Node 8 (F8): use the actual Avrami completion fraction from the
    # gelation timing solver when available. Falls back to 0.999 only when
    # the orchestrator did not pass a timing result (old callers, tests).
    T_gel = props.T_gel
    if timing is not None:
        alpha_final = float(timing.alpha_final)
        timing_diagnostics: dict = {
            "t_gel_onset_s": float(timing.t_gel_onset),
            "alpha_final_from_timing": alpha_final,
            "cooling_rate_effective_K_per_s": float(timing.cooling_rate_effective),
            "mobility_arrest_factor": float(timing.mobility_arrest_factor),
        }
    else:
        # Legacy fallback: assume complete gelation (placeholder pre-Node 8).
        alpha_final = 0.999
        timing_diagnostics = {
            "alpha_final_source": "hardcoded_legacy_fallback",
        }

    # Morphology descriptors: empirical model uses defaults (1D field)
    morph = morphology_descriptors(np.ones(N_r) * phi_dry, 0.0)

    _assumptions = [
        "Power-law pore-concentration scaling", "Pernodet AFM data",
    ]
    if timing is None:
        _assumptions.append(
            "alpha_final hardcoded to 0.999 (no timing result passed); "
            "callers should pass GelationTimingResult for evidence-grade output."
        )
    model_manifest = ModelManifest(
        model_name="L2.Pore.EmpiricalCorrelation",
        evidence_tier=ModelEvidenceTier.SEMI_QUANTITATIVE,
        assumptions=_assumptions,
        diagnostics=timing_diagnostics,
    )
    return GelationResult(
        r_grid=r,
        phi_field=np.ones(N_r) * phi_dry,  # uniform (empirical model)
        pore_size_mean=float(d_pore),
        pore_size_std=float(np.std(pore_sizes)),
        pore_size_distribution=pore_sizes,
        porosity=float(porosity),
        alpha_final=float(min(alpha_final, 0.999)),
        char_wavelength=float(d_pore * 1.5),  # approx periodicity
        T_history=np.array([T_gel]),
        phi_snapshots=None,
        L_domain=0.0,
        grid_spacing=0.0,
        bicontinuous_score=morph['bicontinuous_score'],
        anisotropy=morph['anisotropy'],
        connectivity=morph['connectivity'],
        chord_skewness=morph['skewness'],
        model_tier="empirical_uncalibrated",  # v6.1: honest label (set to "empirical_calibrated" only with user pore data)
        model_manifest=model_manifest,
    )


def solve_gelation(params: SimulationParameters, props: MaterialProperties,
                   R_droplet: float = 1.0e-6,
                   mode: str = 'empirical',
                   timing: 'GelationTimingResult | None' = None) -> GelationResult:
    """Level 2 gelation simulation.

    Parameters
    ----------
    params : SimulationParameters
    props : MaterialProperties
    R_droplet : float
        Microsphere radius [m].
    mode : str
        'empirical' (default) — literature-calibrated pore size correlation.
            Fast (<0.01s), well-calibrated for agarose systems.
        'ch_2d' — 2D Cahn-Hilliard phase-field simulation.
            Mechanistic but requires careful parameter calibration.
        'ch_1d' — 1D radial Cahn-Hilliard (legacy).
    timing : GelationTimingResult, optional
        Output of solve_gelation_timing for the same R_droplet. When provided
        in 'empirical' mode, alpha_final is taken from timing instead of
        the legacy hardcoded 0.999. Mechanistic modes ignore this argument
        because they compute their own alpha trajectory from the phase-field.
    """
    if mode == 'empirical':
        return solve_gelation_empirical(params, props, R_droplet, timing=timing)
    elif mode == 'ch_2d':
        solver = CahnHilliard2DSolver(
            N_grid=params.solver.l2_n_grid,
            dt_initial=params.solver.l2_dt_initial,
            dt_max=params.solver.l2_dt_max,
            arrest_exponent=params.solver.l2_arrest_exponent,
        )
        return solver.solve(params, props, R_droplet)
    elif mode == 'ch_ternary':
        from .ternary_solver import TernaryCahnHilliard2DSolver
        solver = TernaryCahnHilliard2DSolver()
        return solver.solve(params, props, R_droplet=R_droplet,
                            N=params.solver.l2_n_grid)
    elif mode == 'ch_1d':
        solver = CahnHilliardSolver(
            N_r=params.solver.l2_n_r,
            dt_initial=params.solver.l2_dt_initial,
            dt_max=params.solver.l2_dt_max,
            arrest_exponent=params.solver.l2_arrest_exponent,
        )
        return solver.solve(params, props, R_droplet)
    else:
        raise ValueError(f"Unknown gelation mode: {mode!r}. "
                         "Use 'empirical', 'ch_2d', 'ch_ternary', or 'ch_1d'.")
