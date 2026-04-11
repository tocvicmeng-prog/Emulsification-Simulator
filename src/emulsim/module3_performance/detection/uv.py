"""UV detection model: Beer-Lambert absorbance + detector broadening.

Architecture: docs/13_module2_module3_final_implementation_plan.md, Phase C.

Beer-Lambert law: A = epsilon * c * l

where:
    A       : absorbance [AU]
    epsilon : molar extinction coefficient [M^-1 cm^-1]  (converted to SI internally)
    c       : concentration [mol/m^3]
    l       : path length [m]

Output is in mAU (milli-absorbance units) for typical detector display.
"""

from __future__ import annotations

import numpy as np
from scipy.ndimage import gaussian_filter1d


def compute_uv_signal(
    C_outlet: np.ndarray,
    extinction_coeff: float = 36000.0,
    path_length: float = 0.01,
) -> np.ndarray:
    """Convert outlet concentration to UV absorbance signal.

    Args:
        C_outlet: Outlet concentration [mol/m^3], shape (N_t,).
        extinction_coeff: Molar extinction coefficient [1/(M*cm)].
            Default 36000 is typical for BSA at 280 nm.
        path_length: Detector flow cell path length [m].

    Returns:
        UV absorbance signal [mAU], shape (N_t,).
    """
    # Convert extinction_coeff from 1/(M*cm) to 1/(mol/m^3 * m):
    #   1 M = 1000 mol/m^3,  1 cm = 0.01 m
    #   epsilon_SI = epsilon / (1000 * 100) = epsilon / 1e5
    # But simpler: A = epsilon [1/(M*cm)] * c [M] * l [cm]
    #   c [M] = C_outlet [mol/m^3] / 1000
    #   l [cm] = path_length [m] * 100
    c_molar = C_outlet / 1000.0         # mol/m^3 -> mol/L (M)
    l_cm = path_length * 100.0          # m -> cm
    absorbance = extinction_coeff * c_molar * l_cm  # [AU]
    return absorbance * 1000.0          # [mAU]


def apply_detector_broadening(
    signal: np.ndarray,
    time: np.ndarray,
    sigma_detector: float = 1.0,
) -> np.ndarray:
    """Apply Gaussian convolution for extra-column band broadening.

    Args:
        signal: UV signal [mAU], shape (N_t,).
        time: Time array [s], shape (N_t,).
        sigma_detector: Detector band broadening time constant [s].

    Returns:
        Broadened UV signal [mAU], shape (N_t,).
    """
    if sigma_detector <= 0.0 or len(time) < 3:
        return signal.copy()

    # Convert sigma from time to index units
    dt = np.mean(np.diff(time))
    if dt <= 0:
        return signal.copy()

    sigma_idx = sigma_detector / dt
    if sigma_idx < 0.5:
        return signal.copy()

    return gaussian_filter1d(signal.astype(float), sigma=sigma_idx, mode="nearest")
