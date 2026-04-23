"""Pore structure analysis from phase-field simulation results.

Supports both 1D radial and 2D Cartesian phase-field outputs.
"""

from __future__ import annotations
import numpy as np


def structure_factor_1d(phi: np.ndarray, r: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Compute 1D structure factor S(q) from radial composition profile.

    Uses the radial correlation function C(r) = <dphi(0)*dphi(r)> and its
    Fourier transform.

    Returns (q, S_q) arrays.
    """
    N = len(phi)
    dr = r[1] - r[0] if N > 1 else 1.0

    # Fluctuation
    dphi = phi - np.mean(phi)

    # Correlation function via FFT (1D autocorrelation)
    # Pad to avoid circular correlation artifacts
    n_pad = 2 * N
    dphi_pad = np.zeros(n_pad)
    dphi_pad[:N] = dphi

    ft = np.fft.rfft(dphi_pad)
    power = np.abs(ft) ** 2

    # Wavenumber array
    q = np.fft.rfftfreq(n_pad, d=dr) * 2 * np.pi
    S_q = power
    S_max = np.max(S_q)
    if S_max > 0:
        S_q = S_q / S_max

    return q, S_q


def characteristic_wavelength(phi: np.ndarray, r: np.ndarray) -> float:
    """Extract characteristic wavelength from structure factor peak.

    Returns lambda* = 2*pi/q* where q* is the peak wavenumber.
    Returns 0.0 if no clear peak is found.
    """
    q, S_q = structure_factor_1d(phi, r)

    if len(q) < 3:
        return 0.0

    # Skip q=0 (DC component)
    q_search = q[1:]
    S_search = S_q[1:]

    if len(S_search) == 0 or np.max(S_search) <= 0:
        return 0.0

    # Find peak
    idx_peak = np.argmax(S_search)
    q_star = q_search[idx_peak]

    if q_star > 0:
        return 2.0 * np.pi / q_star
    return 0.0


def chord_length_distribution(phi: np.ndarray, r: np.ndarray,
                              threshold: float | None = None) -> np.ndarray:
    """Compute chord lengths in the solvent-rich (pore) phase.

    Binarize the phi field at the threshold, then measure consecutive
    runs of pore-phase grid points.

    Parameters
    ----------
    phi : np.ndarray
        Polymer volume fraction profile.
    r : np.ndarray
        Radial positions [m].
    threshold : float, optional
        Binarization threshold. Default: midpoint between min and max phi.

    Returns
    -------
    np.ndarray
        Array of chord lengths [m]. Empty if no pores found.
    """
    if threshold is None:
        threshold = 0.5 * (np.min(phi) + np.max(phi))

    dr = r[1] - r[0] if len(r) > 1 else 0.0

    # Pore phase: phi < threshold
    is_pore = phi < threshold

    chords = []
    current_length = 0

    for i in range(len(is_pore)):
        if is_pore[i]:
            current_length += 1
        else:
            if current_length > 0:
                chords.append(current_length * dr)
                current_length = 0

    if current_length > 0:
        chords.append(current_length * dr)

    return np.array(chords) if chords else np.array([0.0])


def compute_porosity(phi: np.ndarray, r: np.ndarray,
                     threshold: float | None = None) -> float:
    """Compute porosity (volume fraction of pore space) in spherical geometry.

    Integrates the pore-phase volume using spherical shells.
    """
    if threshold is None:
        threshold = 0.5 * (np.min(phi) + np.max(phi))

    dr = r[1] - r[0] if len(r) > 1 else 0.0

    # Volume of each spherical shell
    shell_vol = 4.0 * np.pi * r ** 2 * dr
    total_vol = np.sum(shell_vol)

    if total_vol <= 0:
        return 0.0

    pore_vol = np.sum(shell_vol[phi < threshold])
    return pore_vol / total_vol


# ─── 2D Cartesian analysis ───────────────────────────────────────────────


def structure_factor_2d(phi_2d: np.ndarray, h: float) -> tuple[np.ndarray, np.ndarray]:
    """Compute radially averaged structure factor S(q) from a 2D composition field.

    Uses 2D FFT and radial binning of |FFT|^2.

    Parameters
    ----------
    phi_2d : np.ndarray
        Composition field, shape (N, N).
    h : float
        Grid spacing [m].

    Returns
    -------
    (q_bins, S_q) : tuple of 1D arrays
        Radially averaged structure factor.
    """
    Ny, Nx = phi_2d.shape
    dphi = phi_2d - np.mean(phi_2d)

    ft = np.fft.fft2(dphi)
    power = np.abs(ft) ** 2

    # Wavenumber grids
    qx = np.fft.fftfreq(Nx, d=h) * 2.0 * np.pi
    qy = np.fft.fftfreq(Ny, d=h) * 2.0 * np.pi
    QX, QY = np.meshgrid(qx, qy)
    Q_mag = np.sqrt(QX**2 + QY**2)

    # Radial binning
    q_max = np.pi / h  # Nyquist
    n_bins = min(Nx, Ny) // 2
    q_edges = np.linspace(0, q_max, n_bins + 1)
    q_bins = 0.5 * (q_edges[:-1] + q_edges[1:])
    S_q = np.zeros(n_bins)

    for ib in range(n_bins):
        mask = (Q_mag >= q_edges[ib]) & (Q_mag < q_edges[ib + 1])
        if np.any(mask):
            S_q[ib] = np.mean(power[mask])

    # Normalise
    S_max = np.max(S_q)
    if S_max > 0:
        S_q /= S_max

    return q_bins, S_q


def characteristic_wavelength_2d(phi_2d: np.ndarray, h: float) -> float:
    """Extract characteristic wavelength from 2D structure factor peak.

    Returns lambda* = 2*pi/q* where q* is the peak wavenumber.
    Returns 0.0 if no clear peak is found.
    """
    q, S_q = structure_factor_2d(phi_2d, h)

    if len(q) < 3:
        return 0.0

    # Skip DC-adjacent bins
    q_search = q[1:]
    S_search = S_q[1:]

    if len(S_search) == 0 or np.max(S_search) <= 0:
        return 0.0

    idx_peak = np.argmax(S_search)
    q_star = q_search[idx_peak]

    if q_star > 0:
        return 2.0 * np.pi / q_star
    return 0.0


def chord_length_distribution_2d(phi_2d: np.ndarray, h: float,
                                  threshold: float | None = None) -> np.ndarray:
    """Compute chord lengths from a 2D composition field.

    Scans all rows and all columns to collect pore-phase run lengths.

    Parameters
    ----------
    phi_2d : np.ndarray
        Composition field, shape (N, N).
    h : float
        Grid spacing [m].
    threshold : float, optional
        Binarization threshold. Default: midpoint of min/max phi.

    Returns
    -------
    np.ndarray
        Array of chord lengths [m].
    """
    if threshold is None:
        threshold = 0.5 * (np.min(phi_2d) + np.max(phi_2d))

    chords: list[float] = []

    # Scan rows
    for row in phi_2d:
        _collect_chords_1d(row < threshold, h, chords)

    # Scan columns
    for j in range(phi_2d.shape[1]):
        _collect_chords_1d(phi_2d[:, j] < threshold, h, chords)

    return np.array(chords) if chords else np.array([0.0])


def _collect_chords_1d(is_pore: np.ndarray, h: float, chords: list) -> None:
    """Helper: collect run lengths from a boolean array."""
    current = 0
    for val in is_pore:
        if val:
            current += 1
        else:
            if current > 0:
                chords.append(current * h)
                current = 0
    if current > 0:
        chords.append(current * h)


def compute_porosity_2d(phi_2d: np.ndarray, threshold: float | None = None) -> float:
    """Compute porosity (area fraction of pore space) from 2D field.

    For a uniform grid, this is simply the fraction of grid points
    below the threshold.
    """
    if threshold is None:
        threshold = 0.5 * (np.min(phi_2d) + np.max(phi_2d))

    return float(np.mean(phi_2d < threshold))


def morphology_descriptors(phi_field: np.ndarray, grid_spacing: float,
                            threshold: float = 0.5) -> dict:
    """Compute morphology descriptors from the phase-field.

    Parameters
    ----------
    phi_field : np.ndarray
        2D polymer volume fraction field (N x N).
    grid_spacing : float
        Grid spacing [m].
    threshold : float
        Threshold for pore vs polymer classification.

    Returns
    -------
    dict with keys:
        bicontinuous_score : float  0-1, how bicontinuous the structure is
        anisotropy : float          0-1, ratio of principal axis lengths
        connectivity : float        0-1, fraction of pore space connected
        skewness : float            skewness of chord length distribution
    """
    if phi_field.ndim != 2:
        # 1D field: return default descriptors
        return {
            'bicontinuous_score': 0.5,
            'anisotropy': 0.0,
            'connectivity': 1.0,
            'skewness': 0.0,
        }

    N = phi_field.shape[0]
    pore_mask = phi_field < threshold
    poly_mask = ~pore_mask

    # 1. Bicontinuous score: fraction of rows AND columns that contain
    #    both pore and polymer phases (indicates interpenetrating structure)
    row_both = sum(1 for i in range(N) if pore_mask[i].any() and poly_mask[i].any())
    col_both = sum(1 for j in range(N) if pore_mask[:, j].any() and poly_mask[:, j].any())
    bicontinuous_score = (row_both + col_both) / (2 * N) if N > 0 else 0.0

    # 2. Anisotropy: compare horizontal vs vertical chord lengths
    h_chords = []
    v_chords = []
    for i in range(N):
        row = pore_mask[i]
        chord = 0
        for j in range(N):
            if row[j]:
                chord += 1
            elif chord > 0:
                h_chords.append(chord * grid_spacing)
                chord = 0
        if chord > 0:
            h_chords.append(chord * grid_spacing)

    for j in range(N):
        col = pore_mask[:, j]
        chord = 0
        for i in range(N):
            if col[i]:
                chord += 1
            elif chord > 0:
                v_chords.append(chord * grid_spacing)
                chord = 0
        if chord > 0:
            v_chords.append(chord * grid_spacing)

    mean_h = float(np.mean(h_chords)) if h_chords else 1.0
    mean_v = float(np.mean(v_chords)) if v_chords else 1.0
    anisotropy = abs(mean_h - mean_v) / max(mean_h, mean_v) if max(mean_h, mean_v) > 0 else 0.0

    # 3. Connectivity: simple flood fill from top-left pore pixel
    #    fraction of total pore pixels reachable
    from scipy import ndimage
    if pore_mask.any():
        labeled, n_features = ndimage.label(pore_mask)
        if n_features > 0:
            largest_component = np.argmax(np.bincount(labeled.flat)[1:]) + 1
            connectivity = np.sum(labeled == largest_component) / np.sum(pore_mask)
        else:
            connectivity = 0.0
    else:
        connectivity = 0.0

    # 4. Skewness of all chord lengths
    all_chords = h_chords + v_chords
    if len(all_chords) > 2:
        from scipy.stats import skew
        skewness = float(skew(all_chords))
    else:
        skewness = 0.0

    return {
        'bicontinuous_score': float(bicontinuous_score),
        'anisotropy': float(anisotropy),
        'connectivity': float(connectivity),
        'skewness': float(skewness),
    }
