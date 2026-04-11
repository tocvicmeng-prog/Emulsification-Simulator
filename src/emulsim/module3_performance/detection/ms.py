"""Semi-quantitative mass spectrometry detection model.

Architecture: docs/13_module2_module3_final_implementation_plan.md, Phase E.

Provides two semi-quantitative MS functions:

1. ``simulate_esi_charge_envelope``: Generate the ESI-MS m/z spectrum
   for a protein of known molecular weight.  Uses an empirical charge-state
   distribution centred on z_avg ~ 0.0778 * MW^0.5 (MW in Da).

2. ``compute_tic``: Total Ion Chromatogram (TIC) as a concentration-
   proportional signal, useful for LC-MS detection of eluting peaks.

Physical basis:
    In positive-ion ESI, proteins acquire multiple protons (charge states).
    The observed m/z for charge state z is:
        m/z = (M + z * m_proton) / z
    where M = molecular mass [Da], m_proton = 1.00728 Da.

    The charge-state distribution is approximately Gaussian in z-space.
    The average charge z_avg follows an empirical power law:
        z_avg ~ 0.0778 * MW^0.5   (Kaltashov & Eyles, 2005)

References:
    Kaltashov, I.A. & Eyles, S.J. (2005). Studies of biomolecular
    conformations and conformational dynamics by mass spectrometry. Mass
    Spectrometry Reviews, 21(1), 37-71.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np


# Mass of a proton [Da]
_M_PROTON = 1.00728


# ─── ESI Charge Envelope ─────────────────────────────────────────────────────

@dataclass
class ESISpectrum:
    """ESI-MS spectrum for a single protein component.

    Attributes:
        mz: m/z values [Da/e], shape (n_points,).
        intensity: Relative intensity [-], shape (n_points,), normalised to max=1.
        molecular_weight: Input molecular weight [Da].
        z_avg: Average charge state [-].
        z_min: Minimum charge state represented.
        z_max: Maximum charge state represented.
    """

    mz: np.ndarray
    intensity: np.ndarray
    molecular_weight: float
    z_avg: float
    z_min: int
    z_max: int


def _charge_state_range(mw_da: float) -> tuple[int, int, float]:
    """Estimate the charge state range for a protein.

    Uses empirical power law: z_avg = 0.0778 * sqrt(MW_Da).
    Charge states typically span +/- 2*sigma_z around z_avg,
    where sigma_z ~ 0.4 * z_avg (empirical).

    Args:
        mw_da: Molecular weight [Da].

    Returns:
        Tuple of (z_min, z_max, z_avg).
    """
    z_avg = 0.0778 * np.sqrt(mw_da)
    sigma_z = max(0.4 * z_avg, 1.0)

    z_min = max(1, int(np.floor(z_avg - 3.0 * sigma_z)))
    z_max = max(z_min + 1, int(np.ceil(z_avg + 3.0 * sigma_z)))

    return z_min, z_max, z_avg


def simulate_esi_charge_envelope(
    molecular_weight: float,
    n_points: int = 100,
    peak_width_da: float | None = None,
) -> ESISpectrum:
    """Generate an ESI-MS m/z spectrum for a protein.

    Each charge state produces a Gaussian peak centred at:
        m/z_z = (MW + z * m_proton) / z

    The intensity at each charge state follows a Gaussian distribution
    in z-space centred on z_avg.

    Args:
        molecular_weight: Protein molecular weight [Da].
        n_points: Number of m/z points in the output spectrum.
        peak_width_da: Full-width at half-maximum of each charge-state
            peak in m/z space [Da].  Defaults to 0.1% of the mean m/z
            (representative of a high-resolution instrument).

    Returns:
        ESISpectrum containing the simulated m/z spectrum.

    Raises:
        ValueError: If molecular_weight is non-positive.

    Examples:
        >>> spec = simulate_esi_charge_envelope(50000.0)
        >>> assert spec.mz.min() > 0
        >>> assert spec.z_avg > 10
    """
    if molecular_weight <= 0:
        raise ValueError(f"molecular_weight must be positive, got {molecular_weight}")

    z_min, z_max, z_avg = _charge_state_range(molecular_weight)
    sigma_z = max(0.4 * z_avg, 1.0)

    # Charge states to model
    z_states = np.arange(z_min, z_max + 1, dtype=int)

    # m/z of each charge state centre
    mz_centres = (molecular_weight + z_states * _M_PROTON) / z_states.astype(float)

    # m/z range for the full spectrum
    mz_lo = mz_centres.min() * 0.95
    mz_hi = mz_centres.max() * 1.05
    mz_axis = np.linspace(mz_lo, mz_hi, n_points)

    # Peak width (sigma in m/z space)
    mean_mz = np.mean(mz_centres)
    if peak_width_da is None:
        peak_width_da = 0.001 * mean_mz  # 0.1% of mean m/z (FWHM)
    sigma_mz = peak_width_da / (2.0 * np.sqrt(2.0 * np.log(2.0)))

    # Charge-state amplitude weights (Gaussian in z)
    z_weights = np.exp(-0.5 * ((z_states - z_avg) / sigma_z) ** 2)
    z_weights /= z_weights.sum()

    # Build spectrum: sum of Gaussians
    spectrum = np.zeros(n_points)
    for z, mz_c, weight in zip(z_states, mz_centres, z_weights):
        peak = weight * np.exp(-0.5 * ((mz_axis - mz_c) / sigma_mz) ** 2)
        spectrum += peak

    # Normalise to max intensity = 1
    if spectrum.max() > 0:
        spectrum /= spectrum.max()

    return ESISpectrum(
        mz=mz_axis,
        intensity=spectrum,
        molecular_weight=molecular_weight,
        z_avg=z_avg,
        z_min=int(z_min),
        z_max=int(z_max),
    )


# ─── Total Ion Chromatogram ───────────────────────────────────────────────────

def compute_tic(
    C_outlet: np.ndarray,
    time: np.ndarray,
    sensitivity: float = 1.0,
) -> np.ndarray:
    """Compute a semi-quantitative Total Ion Chromatogram (TIC).

    In LC-MS, the TIC is the sum of all ion intensities at each time point.
    For a mixture, TIC ~ sum of concentrations (assuming similar ionisation
    efficiencies and detector response factors).

    This model provides a simple concentration-proportional signal:
        TIC(t) = sensitivity * sum_i(C_i_outlet(t))

    Args:
        C_outlet: Outlet concentration(s) vs time [mol/m^3].
            - For a single component: shape (N_t,).
            - For multi-component: shape (n_comp, N_t).
        time: Time array [s], shape (N_t,).
        sensitivity: Detector sensitivity factor [counts / (mol/m^3)].
            Arbitrary units; default 1.0.

    Returns:
        TIC signal [counts], shape (N_t,).
    """
    C = np.asarray(C_outlet, dtype=float)
    C_safe = np.maximum(C, 0.0)

    if C_safe.ndim == 1:
        # Single component
        tic = sensitivity * C_safe
    else:
        # Multi-component: sum over component axis (axis 0)
        tic = sensitivity * np.sum(C_safe, axis=0)

    return tic


def compute_extracted_ion_chromatogram(
    C_outlet: np.ndarray,
    time: np.ndarray,
    component_idx: int = 0,
    sensitivity: float = 1.0,
) -> np.ndarray:
    """Extracted Ion Chromatogram (EIC) for a single component.

    An EIC shows the signal from a specific m/z window, corresponding
    to a single molecular species.

    Args:
        C_outlet: Multi-component outlet concentrations [mol/m^3], shape (n_comp, N_t).
        time: Time array [s], shape (N_t,).
        component_idx: Index of the component to extract.
        sensitivity: Detector sensitivity factor.

    Returns:
        EIC signal [counts], shape (N_t,).
    """
    C = np.asarray(C_outlet, dtype=float)
    if C.ndim == 1:
        return compute_tic(C, time, sensitivity)

    C_comp = np.maximum(C[component_idx, :], 0.0)
    return sensitivity * C_comp
