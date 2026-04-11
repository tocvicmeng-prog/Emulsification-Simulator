"""Fluorescence detection model for chromatographic monitoring.

Architecture: docs/13_module2_module3_final_implementation_plan.md, Phase E.

Fluorescence signal in the linear regime (low absorbance, A < 0.05 AU):

    F = Phi * epsilon * c * l

where:
    F       : Fluorescence signal [relative fluorescence units, RFU]
    Phi     : Quantum yield [-]   (0 < Phi <= 1)
    epsilon : Molar extinction coefficient at excitation wavelength [1/(M*cm)]
    c       : Concentration [mol/m^3]
    l       : Path length [m]

Note: This is the linear approximation valid when the inner filter effect
is negligible.  For concentrated solutions (A > 0.05), the signal
saturates and may decrease (inner filter effect), which is not modelled here.

Typical fluorophore parameters:
    - Fluorescein:    Phi=0.92, epsilon=76,900 M^-1 cm^-1 (at 490 nm ex)
    - Rhodamine B:    Phi=0.65, epsilon=106,000 M^-1 cm^-1 (at 554 nm ex)
    - mCherry:        Phi=0.22, epsilon=72,000 M^-1 cm^-1 (at 587 nm ex)
    - Intrinsic Trp:  Phi=0.13, epsilon=5,600 M^-1 cm^-1 (at 280 nm ex)
"""

from __future__ import annotations

import numpy as np


def compute_fluorescence_signal(
    C: float | np.ndarray,
    quantum_yield: float = 0.92,
    extinction_coeff: float = 76900.0,
    path_length: float = 0.01,
) -> np.ndarray:
    """Convert concentration to fluorescence signal in the linear regime.

    F = Phi * epsilon [1/(M*cm)] * c [M] * l [cm]

    Units are converted internally:
        c [M] = C [mol/m^3] / 1000
        l [cm] = path_length [m] * 100

    Args:
        C: Concentration [mol/m^3]. Accepts scalar or array.
        quantum_yield: Fluorescence quantum yield [-], 0 < Phi <= 1.
            Default 0.92 is for fluorescein.
        extinction_coeff: Molar extinction coefficient at excitation wavelength
            [1/(M*cm)].  Default 76,900 is for fluorescein at 490 nm.
        path_length: Detector flow cell optical path length [m].

    Returns:
        Fluorescence signal [RFU], same shape as C.

    Raises:
        ValueError: If quantum_yield is outside (0, 1].
    """
    if not 0 < quantum_yield <= 1.0:
        raise ValueError(
            f"quantum_yield must be in (0, 1], got {quantum_yield}"
        )
    if extinction_coeff <= 0:
        raise ValueError(
            f"extinction_coeff must be positive, got {extinction_coeff}"
        )
    if path_length <= 0:
        raise ValueError(
            f"path_length must be positive, got {path_length}"
        )

    C_arr = np.asarray(C, dtype=float)
    C_safe = np.maximum(C_arr, 0.0)

    c_molar = C_safe / 1000.0       # mol/m^3 -> mol/L (M)
    l_cm = path_length * 100.0      # m -> cm

    signal = quantum_yield * extinction_coeff * c_molar * l_cm
    return signal


def fluorescence_detection_limit(
    quantum_yield: float = 0.92,
    extinction_coeff: float = 76900.0,
    path_length: float = 0.01,
    min_detectable_rfu: float = 0.001,
) -> float:
    """Estimate the minimum detectable concentration.

    Args:
        quantum_yield: Fluorescence quantum yield [-].
        extinction_coeff: Molar extinction coefficient [1/(M*cm)].
        path_length: Path length [m].
        min_detectable_rfu: Minimum detectable signal [RFU].

    Returns:
        Minimum detectable concentration [mol/m^3].
    """
    l_cm = path_length * 100.0
    c_molar = min_detectable_rfu / (quantum_yield * extinction_coeff * l_cm)
    return c_molar * 1000.0  # mol/L -> mol/m^3
