"""Spatial discretization utilities for 1D radial phase-field solver."""

from __future__ import annotations
import numpy as np
from scipy import sparse


def create_radial_grid(R: float, N_r: int) -> tuple[np.ndarray, float]:
    """Create 1D radial grid from r=0 to r=R.

    Uses uniform spacing. First point is at r=dr/2 (staggered to avoid r=0 singularity).
    Grid is built in normalised coordinates (0,1] then scaled to avoid
    floating-point cancellation at very small R.

    Returns (r, dr).
    """
    dr = R / N_r
    # Build in unit coordinates then scale
    r_unit = (np.arange(N_r) + 0.5) / N_r  # exact in [0.5/N, 1-0.5/N]
    r = r_unit * R
    return r, dr


def build_laplacian_matrix(r: np.ndarray, dr: float) -> sparse.csr_matrix:
    """Build sparse matrix for the spherical Laplacian operator.

    nabla^2 phi = (1/r^2) d/dr(r^2 dphi/dr)

    Conservative discretization on cell-centered grid with no-flux BCs.
    Internally works in normalised coordinates to avoid precision loss
    at very small R.  The returned matrix is in physical units (1/dr^2 scaling).
    """
    N = len(r)

    # Work in normalised coordinates: rhat = r / R  where R = r_max + dr/2
    R = r[-1] + dr / 2.0
    dr_hat = 1.0 / N  # normalised spacing (exact)

    # Normalised cell centres: (i+0.5)/N
    rhat = (np.arange(N) + 0.5) * dr_hat

    # Normalised face positions: i/N  (exact integers / N)
    rhat_face = np.arange(N + 1) * dr_hat

    # Face areas
    A_face = rhat_face ** 2

    # Cell-centre r^2
    rhat_sq = rhat ** 2

    # Coefficients in normalised space (dimensionless Laplacian, 1/dr_hat^2 = N^2)
    coeff_right = A_face[1:] / (dr_hat ** 2 * rhat_sq)
    coeff_left = A_face[:N] / (dr_hat ** 2 * rhat_sq)

    main_diag = -(coeff_right + coeff_left)
    upper_diag = coeff_right[:-1].copy()
    lower_diag = coeff_left[1:].copy()

    # BCs: coeff_left[0] = 0 (face at origin has zero area) — nothing to fold.
    # No-flux at r=R: fold coeff_right[-1] into main
    main_diag[-1] += coeff_right[-1]

    # Build in normalised coordinates → L_hat has units 1/dr_hat^2 = N^2
    # Physical Laplacian = L_hat / R^2
    L_hat = sparse.diags(
        [lower_diag, main_diag, upper_diag],
        [-1, 0, 1],
        shape=(N, N),
        format='csr',
    )
    return L_hat / (R ** 2)


def apply_gradient(phi: np.ndarray, r: np.ndarray, dr: float) -> np.ndarray:
    """Compute dphi/dr using central differences with no-flux BCs."""
    N = len(phi)
    grad = np.zeros(N)

    # Interior
    grad[1:-1] = (phi[2:] - phi[:-2]) / (2.0 * dr)

    # BCs: dphi/dr = 0 at both ends
    grad[0] = 0.0
    grad[-1] = 0.0

    return grad


def apply_divergence_flux(J: np.ndarray, r: np.ndarray, dr: float) -> np.ndarray:
    """Compute div(J) in spherical coordinates: (1/r^2) d(r^2 J)/dr."""
    N = len(J)
    div = np.zeros(N)

    r_sq_J = r ** 2 * J

    # Interior: central differences
    div[1:-1] = (r_sq_J[2:] - r_sq_J[:-2]) / (2.0 * dr * r[1:-1] ** 2)

    # Boundaries with no-flux
    if N > 1:
        div[0] = (r_sq_J[1]) / (dr * max(r[0] ** 2, 1e-30))
        div[-1] = (-r_sq_J[-2]) / (dr * r[-1] ** 2)

    return div
