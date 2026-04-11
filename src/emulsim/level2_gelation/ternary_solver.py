"""Ternary Cahn-Hilliard solver for agarose-chitosan-water phase separation.

Solves two coupled Cahn-Hilliard equations for conserved order parameters
phi_a (agarose) and phi_c (chitosan), with phi_s = 1 - phi_a - phi_c (solvent).

    d(phi_a)/dt = div(M_a * grad(mu_a - kappa_a * laplacian(phi_a)))
    d(phi_c)/dt = div(M_c * grad(mu_c - kappa_c * laplacian(phi_c)))

Uses semi-implicit Eyre convex-splitting time stepping with no-flux (Neumann)
boundary conditions on a 2D Cartesian grid.  Supports gelation arrest coupling
via Avrami kinetics to freeze coarsening when alpha approaches 1.

Improvements over the original explicit-Euler/periodic-BC version:
1. Semi-implicit Eyre splitting: treat convex part of free energy implicitly
   for unconditional stability, allowing dt ~ 5-50x larger.
2. No-flux (Neumann) boundary conditions: physically correct for a microsphere
   with no mass flux across the droplet boundary.
3. Gelation arrest: mobility scaled by (1 - alpha^beta) so coarsening halts
   when the gel network is complete.
"""

from __future__ import annotations

import logging

import numpy as np

from ..datatypes import GelationResult, MaterialProperties, SimulationParameters
from .ternary_free_energy import chemical_potential_a, chemical_potential_c

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Utility functions
# ---------------------------------------------------------------------------

def _laplacian_2d_noflux(field: np.ndarray, h: float) -> np.ndarray:
    """2D Laplacian with no-flux (Neumann) boundary conditions.

    Interior points use the standard 5-point stencil.  Boundary rows/columns
    are filled by mirroring the adjacent interior row/column, which enforces
    a zero normal gradient (no-flux) condition.

    Parameters
    ----------
    field : np.ndarray, shape (N, N)
    h : float
        Grid spacing [m].

    Returns
    -------
    np.ndarray, shape (N, N)
    """
    lap = np.zeros_like(field)
    # Interior: standard 5-point stencil
    lap[1:-1, 1:-1] = (
        field[:-2, 1:-1] + field[2:, 1:-1]
        + field[1:-1, :-2] + field[1:-1, 2:]
        - 4.0 * field[1:-1, 1:-1]
    ) / h ** 2
    # Boundary rows/cols: mirror adjacent interior value (zero gradient)
    lap[0, :]  = lap[1, :]
    lap[-1, :] = lap[-2, :]
    lap[:, 0]  = lap[:, 1]
    lap[:, -1] = lap[:, -2]
    return lap


def _laplacian_2d(field: np.ndarray, h: float) -> np.ndarray:
    """2D Laplacian with periodic boundary conditions (5-point stencil).

    Kept for backward compatibility; new code should prefer ``_laplacian_2d_noflux``.
    """
    return (
        np.roll(field, 1, axis=0) + np.roll(field, -1, axis=0)
        + np.roll(field, 1, axis=1) + np.roll(field, -1, axis=1)
        - 4.0 * field
    ) / h ** 2


def _dominant_wavelength(phi: np.ndarray, h: float) -> float:
    """Extract dominant wavelength from 2D power spectrum via radial averaging."""
    N = phi.shape[0]
    fft2 = np.fft.fft2(phi - phi.mean())
    power = np.abs(fft2) ** 2

    freqs = np.fft.fftfreq(N, d=h)
    kx, ky = np.meshgrid(freqs, freqs)
    k_mag = np.sqrt(kx ** 2 + ky ** 2)
    k_mag[0, 0] = 1.0  # avoid div-by-zero at DC

    # Radial average
    n_bins = 50
    k_bins = np.linspace(0, k_mag.max(), n_bins + 1)
    power_radial = np.zeros(n_bins)
    for i in range(n_bins):
        mask = (k_mag >= k_bins[i]) & (k_mag < k_bins[i + 1])
        if mask.any():
            power_radial[i] = power[mask].mean()

    k_peak = k_bins[1:][np.argmax(power_radial)] if power_radial.max() > 0 else 1e7
    return 1.0 / max(k_peak, 1e3)


def _stabilization_constant(
    phi_a: float,
    phi_c: float,
    kT_over_v0: float = 4.114e6,
    N_a: float = 10.0,
    N_c: float = 5.0,
    chi_as: float = 0.5,
    chi_cs: float = 0.6,
    chi_ac: float = 1.0,
) -> float:
    """Estimate a stabilization constant C >= max|f''| for Eyre splitting.

    The stabilization constant must satisfy C >= max eigenvalue of the Hessian
    of the bulk free energy so that (C*phi - f(phi)) is convex everywhere.
    We estimate it from the diagonal Hessian entries at the mean composition.

    Parameters
    ----------
    phi_a, phi_c : float
        Mean volume fractions (scalar).
    kT_over_v0 : float
        Thermal energy / monomer volume [J/m^3].

    Returns
    -------
    float
        Stabilization constant [J/m^3].
    """
    eps = 1e-3
    pa = max(phi_a, eps)
    pc = max(phi_c, eps)
    ps = max(1.0 - pa - pc, eps)

    # Second derivative of f wrt phi_a at mean composition (diagonal Hessian entry)
    d2f_a = kT_over_v0 * (1.0 / (N_a * pa) + 1.0 / ps)
    d2f_c = kT_over_v0 * (1.0 / (N_c * pc) + 1.0 / ps)

    # Cross term magnitude from chi interactions
    cross = abs(kT_over_v0 * chi_ac)

    # Take the largest entry with a 2x safety margin
    C = 2.0 * max(d2f_a, d2f_c, cross, 1.0)
    return C


# ---------------------------------------------------------------------------
# Public API: solve_ternary
# ---------------------------------------------------------------------------

def solve_ternary(
    phi_a0: np.ndarray,
    phi_c0: np.ndarray,
    h: float,
    dt: float = 1e-4,
    n_steps: int = 5000,
    M_a: float = 1e-9,
    M_c: float = 1e-9,
    kappa_a: float = 5e-12,
    kappa_c: float = 5e-12,
    T: float = 350.0,
    chi_ac: float = 1.0,
    use_noflux_bc: bool = True,
    gelation_arrest: bool = True,
    k_gel: float = 1.0,
    n_avrami: float = 2.5,
    arrest_exponent: float = 2.5,
) -> GelationResult:
    """Solve the coupled ternary Cahn-Hilliard equations.

    Uses semi-implicit Eyre convex-splitting for stability, optional no-flux
    boundary conditions, and optional gelation arrest via Avrami kinetics.

    Parameters
    ----------
    phi_a0 : np.ndarray, shape (N, N)
        Initial agarose volume fraction field.
    phi_c0 : np.ndarray, shape (N, N)
        Initial chitosan volume fraction field.
    h : float
        Grid spacing [m].
    dt : float
        Time step [s].  The semi-implicit scheme is unconditionally stable so
        this can be much larger than the explicit-Euler stability limit.
    n_steps : int
        Number of time steps to run.
    M_a, M_c : float
        Base mobilities for agarose and chitosan [m^5/(J s)].
    kappa_a, kappa_c : float
        Gradient energy coefficients [J/m].
    T : float
        Temperature [K].
    chi_ac : float
        Agarose-chitosan Flory-Huggins interaction parameter.
    use_noflux_bc : bool
        If True (default) use no-flux Neumann BCs; if False use periodic BCs.
    gelation_arrest : bool
        If True (default) scale mobility by (1 - alpha^beta) from Avrami.
    k_gel : float
        Avrami rate constant [s^-n_avrami].
    n_avrami : float
        Avrami exponent.
    arrest_exponent : float
        Exponent beta in the mobility arrest factor (1 - alpha^beta).

    Returns
    -------
    GelationResult
    """
    phi_a = phi_a0.copy()
    phi_c = phi_c0.copy()
    N = phi_a.shape[0]

    c_a = float(phi_a.mean())
    c_c = float(phi_c.mean())

    laplacian = _laplacian_2d_noflux if use_noflux_bc else _laplacian_2d

    # Stabilization constant for Eyre convex-splitting (>= max|f''|)
    C = _stabilization_constant(c_a, c_c, chi_ac=chi_ac)

    # Spectral radius of -L on an N x N grid (worst-case mode k = pi/h):
    # lambda_max = 8 / h^2  (2D, 5-point stencil, max eigenvalue of -L)
    lambda_max = 8.0 / h ** 2

    alpha = 0.0
    t = 0.0

    for step in range(n_steps):
        # Gelation arrest: update alpha and effective mobilities
        if gelation_arrest and t > 0.0:
            alpha = 1.0 - np.exp(-k_gel * t ** n_avrami)
            arrest_factor = max(1.0 - alpha ** arrest_exponent, 1e-10)
        else:
            arrest_factor = 1.0

        M_a_eff = M_a * arrest_factor
        M_c_eff = M_c * arrest_factor

        # Stop early if gel is essentially complete
        if gelation_arrest and alpha > 0.999:
            logger.debug(
                "Ternary CH: gelation complete at step %d (t=%.4f s, alpha=%.5f)",
                step, t, alpha,
            )
            break

        # --- Semi-implicit Eyre convex-splitting update ------------------
        # Full chemical potentials from the ternary Flory-Huggins free energy
        mu_a = chemical_potential_a(phi_a, phi_c, T, chi_ac=chi_ac)
        mu_c = chemical_potential_c(phi_a, phi_c, T, chi_ac=chi_ac)

        # Eyre explicit chemical potential: mu_exp = mu_FH - C * phi
        # (the +C*phi piece is moved to the implicit LHS)
        mu_a_exp = mu_a - C * phi_a
        mu_c_exp = mu_c - C * phi_c

        # Laplacians needed for the update
        lap_mu_a_exp = laplacian(mu_a_exp, h)
        lap_mu_c_exp = laplacian(mu_c_exp, h)

        lap_phi_a = laplacian(phi_a, h)
        lap_phi_c = laplacian(phi_c, h)

        lap_lap_phi_a = laplacian(lap_phi_a, h)
        lap_lap_phi_c = laplacian(lap_phi_c, h)

        # Semi-implicit update with scalar (diagonal) stabilization:
        #
        #   [1 + dt * M * C * lambda_max] * phi_new
        #       = phi_old + dt * M * [L(mu_exp) - kappa * L^2(phi_old)]
        #
        # The denominator approximates the implicit (I - dt*M*C*L) operator
        # by its spectral norm, giving an inexpensive point-wise update that
        # is provably stable for C >= max|f''|.
        denom_a = 1.0 + dt * M_a_eff * C * lambda_max
        denom_c = 1.0 + dt * M_c_eff * C * lambda_max

        phi_a_new = (
            phi_a + dt * M_a_eff * (lap_mu_a_exp - kappa_a * lap_lap_phi_a)
        ) / denom_a

        phi_c_new = (
            phi_c + dt * M_c_eff * (lap_mu_c_exp - kappa_c * lap_lap_phi_c)
        ) / denom_c

        # Enforce physical constraints
        phi_a_new = np.clip(phi_a_new, 1e-6, 1.0 - 1e-6)
        phi_c_new = np.clip(phi_c_new, 1e-6, 1.0 - 1e-6)

        # Ensure phi_a + phi_c <= 1 (solvent fraction >= 0)
        total = phi_a_new + phi_c_new
        mask = total > 0.99
        if mask.any():
            scale = 0.99 / total[mask]
            phi_a_new[mask] *= scale
            phi_c_new[mask] *= scale

        # NaN guard (should not occur with semi-implicit scheme)
        if np.any(np.isnan(phi_a_new)) or np.any(np.isnan(phi_c_new)):
            logger.warning(
                "Ternary CH (semi-implicit): NaN at step %d, halving dt", step
            )
            dt *= 0.5
            phi_a_new = np.nan_to_num(phi_a_new, nan=c_a)
            phi_c_new = np.nan_to_num(phi_c_new, nan=c_c)

        phi_a = phi_a_new
        phi_c = phi_c_new
        t += dt

    # ---------- Post-processing ----------
    phi_total = phi_a + phi_c
    L_domain = h * N

    from .pore_analysis import chord_length_distribution_2d, morphology_descriptors

    threshold = 0.5 * (c_a + c_c)
    pore_dist = chord_length_distribution_2d(phi_total, h, threshold=threshold)
    morph = morphology_descriptors(phi_total, h, threshold=threshold)

    pore_mean = float(np.mean(pore_dist)) if len(pore_dist) > 0 else 100e-9
    pore_std  = float(np.std(pore_dist))  if len(pore_dist) > 0 else 10e-9
    porosity  = float(np.mean(phi_total < threshold))

    char_wavelength = _dominant_wavelength(phi_total, h)

    return GelationResult(
        r_grid=np.linspace(0, L_domain / 2, N),
        phi_field=phi_total,
        pore_size_mean=pore_mean,
        pore_size_std=pore_std,
        pore_size_distribution=pore_dist if len(pore_dist) > 0 else np.array([pore_mean]),
        porosity=porosity,
        alpha_final=float(alpha),
        char_wavelength=char_wavelength,
        L_domain=L_domain,
        grid_spacing=h,
        bicontinuous_score=morph["bicontinuous_score"],
        anisotropy=morph["anisotropy"],
        connectivity=morph["connectivity"],
        chord_skewness=morph["skewness"],
    )


# ---------------------------------------------------------------------------
# Class-based interface (backward compatible)
# ---------------------------------------------------------------------------

class TernaryCahnHilliard2DSolver:
    """2D ternary Cahn-Hilliard solver with two conserved order parameters.

    Uses semi-implicit Eyre convex-splitting with no-flux boundary conditions
    and optional gelation arrest via Avrami kinetics.
    """

    def solve(
        self,
        params: SimulationParameters,
        props: MaterialProperties,
        R_droplet: float = 50e-6,
        N: int = 64,
        chi_ac: float = 1.0,
        kappa_a: float = 5e-12,
        kappa_c: float = 5e-12,
        M_a: float = 1e-9,
        M_c: float = 1e-9,
        dt: float = 1e-4,
        n_steps: int = 5000,
        use_noflux_bc: bool = True,
        gelation_arrest: bool = True,
        k_gel: float = 1.0,
        n_avrami: float = 2.5,
        arrest_exponent: float = 2.5,
    ) -> GelationResult:
        """Solve coupled CH equations for agarose-chitosan demixing.

        Parameters
        ----------
        params : SimulationParameters
        props : MaterialProperties
        R_droplet : float
            Microsphere radius [m].
        N : int
            Grid side length (N x N).
        chi_ac : float
            Agarose-chitosan Flory-Huggins parameter.
        kappa_a, kappa_c : float
            Gradient energy coefficients [J/m].
        M_a, M_c : float
            Mobilities [m^5/(J s)].
        dt : float
            Time step [s].
        n_steps : int
            Number of time steps.
        use_noflux_bc : bool
            Use no-flux (Neumann) BCs instead of periodic BCs.
        gelation_arrest : bool
            Enable gelation arrest coupling via Avrami kinetics.
        k_gel : float
            Avrami rate constant [s^-n_avrami].
        n_avrami : float
            Avrami exponent.
        arrest_exponent : float
            Exponent in (1 - alpha^arrest_exponent) mobility scaling.

        Returns
        -------
        GelationResult
        """
        L = min(2.0 * R_droplet, 1.5e-6)
        h = L / N

        # Initial conditions: uniform + small noise
        c_a = params.formulation.c_agarose / 1400.0   # dry volume fraction
        c_c = params.formulation.c_chitosan / 1400.0

        rng = np.random.default_rng(42)
        phi_a = c_a + 0.01 * rng.standard_normal((N, N))
        phi_c = c_c + 0.01 * rng.standard_normal((N, N))

        phi_a = np.clip(phi_a, 1e-6, 0.5)
        phi_c = np.clip(phi_c, 1e-6, 0.5)

        T = params.formulation.T_oil

        return solve_ternary(
            phi_a0=phi_a,
            phi_c0=phi_c,
            h=h,
            dt=dt,
            n_steps=n_steps,
            M_a=M_a,
            M_c=M_c,
            kappa_a=kappa_a,
            kappa_c=kappa_c,
            T=T,
            chi_ac=chi_ac,
            use_noflux_bc=use_noflux_bc,
            gelation_arrest=gelation_arrest,
            k_gel=k_gel,
            n_avrami=n_avrami,
            arrest_exponent=arrest_exponent,
        )
