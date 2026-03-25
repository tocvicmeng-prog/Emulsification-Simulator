"""Spatial discretization utilities for 1D radial and 2D Cartesian phase-field solvers."""

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


# ─── 2D Cartesian discretization ─────────────────────────────────────────


def create_2d_grid(L_domain: float, N: int) -> tuple[np.ndarray, np.ndarray, float]:
    """Create 2D Cartesian grid on [0, L_domain] x [0, L_domain].

    Cell-centred grid with uniform spacing h = L_domain / N.

    Returns (x, y, h) where x and y are 1D coordinate arrays of length N.
    """
    h = L_domain / N
    coords = (np.arange(N) + 0.5) * h
    return coords, coords, h


def build_laplacian_2d(N: int, h: float) -> sparse.csr_matrix:
    """Standard 5-point Laplacian on an N x N grid with Neumann (no-flux) BCs.

    Uses Kronecker product: L_2D = I (x) L_1D + L_1D (x) I
    where L_1D is the 1D second-derivative with Neumann BCs.

    Unknowns are in row-major order: k = i*N + j.

    Parameters
    ----------
    N : int
        Grid points per side.
    h : float
        Uniform grid spacing [m].

    Returns
    -------
    sparse.csr_matrix
        Shape (N*N, N*N).
    """
    # 1D second derivative with Neumann BCs
    e = np.ones(N)
    L1d = sparse.diags([e[:-1], -2.0 * e, e[:-1]], [-1, 0, 1],
                        shape=(N, N), format='lil')
    # Neumann BCs: ghost value equals interior value -> fold into diagonal
    L1d[0, 0] = -1.0
    L1d[-1, -1] = -1.0
    L1d = L1d.tocsr()

    I_N = sparse.eye(N, format='csr')

    L2d = sparse.kron(I_N, L1d, format='csr') + sparse.kron(L1d, I_N, format='csr')
    return L2d / (h * h)


def build_mobility_laplacian_2d(M_field: np.ndarray, N: int,
                                 h: float) -> sparse.csr_matrix:
    """Build div(M * grad(.)) operator with face-centred mobility averaging.

    This is the correct conservative discretization for variable mobility:

        div(M grad phi)_{i,j} = (1/h^2) * [
            M_{i+1/2,j}*(phi_{i+1,j} - phi_{i,j}) - M_{i-1/2,j}*(phi_{i,j} - phi_{i-1,j})
          + M_{i,j+1/2}*(phi_{i,j+1} - phi_{i,j}) - M_{i,j-1/2}*(phi_{i,j} - phi_{i,j-1})
        ]

    where M_{i+1/2,j} = 0.5*(M_{i,j} + M_{i+1,j}).

    No-flux BCs: faces at domain boundary have zero flux (no neighbour term).

    Parameters
    ----------
    M_field : np.ndarray
        Mobility at each grid point, shape (N, N).
    N : int
        Grid points per side.
    h : float
        Uniform grid spacing [m].

    Returns
    -------
    sparse.csr_matrix
        Shape (N*N, N*N), the mobility-weighted Laplacian operator.
    """
    total = N * N
    inv_h2 = 1.0 / (h * h)

    # Pre-compute face-averaged mobilities
    # x-direction faces: M_{i+1/2,j} for i=0..N-2, j=0..N-1
    Mx_face = 0.5 * (M_field[:-1, :] + M_field[1:, :])  # (N-1, N)
    # y-direction faces: M_{i,j+1/2} for i=0..N-1, j=0..N-2
    My_face = 0.5 * (M_field[:, :-1] + M_field[:, 1:])   # (N, N-1)

    # Vectorized assembly using COO format
    # Each interior face contributes to 4 matrix entries (2 off-diagonal, 2 diagonal)
    rows = []
    cols = []
    vals = []

    # x-direction faces (between rows i and i+1)
    for i in range(N - 1):
        j_arr = np.arange(N)
        k_lo = i * N + j_arr       # (i, j)
        k_hi = (i + 1) * N + j_arr  # (i+1, j)
        m_face = Mx_face[i, :] * inv_h2

        # Off-diagonal: k_lo -> k_hi and k_hi -> k_lo
        rows.extend(k_lo.tolist())
        cols.extend(k_hi.tolist())
        vals.extend(m_face.tolist())

        rows.extend(k_hi.tolist())
        cols.extend(k_lo.tolist())
        vals.extend(m_face.tolist())

        # Diagonal contributions: -m_face for both k_lo and k_hi
        rows.extend(k_lo.tolist())
        cols.extend(k_lo.tolist())
        vals.extend((-m_face).tolist())

        rows.extend(k_hi.tolist())
        cols.extend(k_hi.tolist())
        vals.extend((-m_face).tolist())

    # y-direction faces (between columns j and j+1)
    for j in range(N - 1):
        i_arr = np.arange(N)
        k_lo = i_arr * N + j       # (i, j)
        k_hi = i_arr * N + (j + 1)  # (i, j+1)
        m_face = My_face[:, j] * inv_h2

        # Off-diagonal
        rows.extend(k_lo.tolist())
        cols.extend(k_hi.tolist())
        vals.extend(m_face.tolist())

        rows.extend(k_hi.tolist())
        cols.extend(k_lo.tolist())
        vals.extend(m_face.tolist())

        # Diagonal
        rows.extend(k_lo.tolist())
        cols.extend(k_lo.tolist())
        vals.extend((-m_face).tolist())

        rows.extend(k_hi.tolist())
        cols.extend(k_hi.tolist())
        vals.extend((-m_face).tolist())

    L_M = sparse.coo_matrix((vals, (rows, cols)), shape=(total, total))
    return L_M.tocsr()
