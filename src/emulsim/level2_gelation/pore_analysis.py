"""Pore structure analysis from phase-field simulation results."""

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
                              threshold: float = None) -> np.ndarray:
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
                     threshold: float = None) -> float:
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
