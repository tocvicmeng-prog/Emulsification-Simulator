"""Ternary Cahn-Hilliard solver for agarose-chitosan-water phase separation.

Solves two coupled Cahn-Hilliard equations for conserved order parameters
phi_a (agarose) and phi_c (chitosan), with phi_s = 1 - phi_a - phi_c (solvent).

    d(phi_a)/dt = div(M_a * grad(mu_a - kappa_a * laplacian(phi_a)))
    d(phi_c)/dt = div(M_c * grad(mu_c - kappa_c * laplacian(phi_c)))

Uses explicit Euler time stepping with periodic boundary conditions on a
2D Cartesian grid.  Designed to capture polymer-polymer demixing that the
single-field CH solver cannot resolve.
"""

from __future__ import annotations

import logging

import numpy as np

from ..datatypes import GelationResult, MaterialProperties, SimulationParameters
from .ternary_free_energy import chemical_potential_a, chemical_potential_c

logger = logging.getLogger(__name__)


class TernaryCahnHilliard2DSolver:
    """2D ternary Cahn-Hilliard solver with two conserved order parameters."""

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

        # ---------- Time stepping (explicit Euler) ----------
        for step in range(n_steps):
            # Chemical potentials from FH free energy
            mu_a = chemical_potential_a(phi_a, phi_c, T)
            mu_c = chemical_potential_c(phi_a, phi_c, T)

            # Laplacians (periodic BC via np.roll)
            lap_phi_a = _laplacian_2d(phi_a, h)
            lap_phi_c = _laplacian_2d(phi_c, h)

            # Full chemical potential with gradient energy penalty
            mu_a_full = mu_a - kappa_a * lap_phi_a
            mu_c_full = mu_c - kappa_c * lap_phi_c

            # Diffusion step: d(phi)/dt = M * laplacian(mu_full)
            lap_mu_a = _laplacian_2d(mu_a_full, h)
            lap_mu_c = _laplacian_2d(mu_c_full, h)

            phi_a += dt * M_a * lap_mu_a
            phi_c += dt * M_c * lap_mu_c

            # Enforce physical constraints
            phi_a = np.clip(phi_a, 1e-6, 1.0 - 1e-6)
            phi_c = np.clip(phi_c, 1e-6, 1.0 - 1e-6)

            # Ensure phi_a + phi_c <= 1  (solvent fraction >= 0)
            total = phi_a + phi_c
            mask = total > 0.99
            if mask.any():
                scale = 0.99 / total[mask]
                phi_a[mask] *= scale
                phi_c[mask] *= scale

            # Adaptive stability check
            if np.any(np.isnan(phi_a)) or np.any(np.isnan(phi_c)):
                logger.warning("Ternary CH: NaN at step %d, halving dt", step)
                dt *= 0.5
                phi_a = np.nan_to_num(phi_a, nan=c_a)
                phi_c = np.nan_to_num(phi_c, nan=c_c)

        # ---------- Post-processing ----------
        phi_total = phi_a + phi_c

        # Pore analysis: pore = solvent-rich region (low phi_total)
        from .pore_analysis import chord_length_distribution_2d, morphology_descriptors

        threshold = 0.5 * (c_a + c_c)
        pore_dist = chord_length_distribution_2d(phi_total, h, threshold=threshold)
        morph = morphology_descriptors(phi_total, h, threshold=threshold)

        pore_mean = float(np.mean(pore_dist)) if len(pore_dist) > 0 else 100e-9
        pore_std = float(np.std(pore_dist)) if len(pore_dist) > 0 else 10e-9
        porosity = float(np.mean(phi_total < threshold))

        # Characteristic wavelength from FFT of phi_total
        char_wavelength = _dominant_wavelength(phi_total, h)

        return GelationResult(
            r_grid=np.linspace(0, L / 2, N),
            phi_field=phi_total,
            pore_size_mean=pore_mean,
            pore_size_std=pore_std,
            pore_size_distribution=pore_dist if len(pore_dist) > 0 else np.array([pore_mean]),
            porosity=porosity,
            alpha_final=1.0,
            char_wavelength=char_wavelength,
            L_domain=L,
            grid_spacing=h,
            bicontinuous_score=morph["bicontinuous_score"],
            anisotropy=morph["anisotropy"],
            connectivity=morph["connectivity"],
            chord_skewness=morph["skewness"],
        )


# ---------------------------------------------------------------------------
# Utility functions
# ---------------------------------------------------------------------------

def _laplacian_2d(field: np.ndarray, h: float) -> np.ndarray:
    """2D Laplacian with periodic boundary conditions (5-point stencil)."""
    return (
        np.roll(field, 1, axis=0) + np.roll(field, -1, axis=0)
        + np.roll(field, 1, axis=1) + np.roll(field, -1, axis=1)
        - 4.0 * field
    ) / h**2


def _dominant_wavelength(phi: np.ndarray, h: float) -> float:
    """Extract dominant wavelength from 2D power spectrum via radial averaging."""
    N = phi.shape[0]
    fft2 = np.fft.fft2(phi - phi.mean())
    power = np.abs(fft2) ** 2

    freqs = np.fft.fftfreq(N, d=h)
    kx, ky = np.meshgrid(freqs, freqs)
    k_mag = np.sqrt(kx**2 + ky**2)
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
