"""Objective function definitions for multi-objective optimization.

Maps simulation FullResult → objective vector for BoTorch.
All objectives are MINIMIZED (lower is better).
"""

from __future__ import annotations

import numpy as np

from ..datatypes import FullResult


# Target values
TARGET_D32 = 2.0e-6      # 2 µm
TARGET_PORE = 80.0e-9    # 80 nm (centre of 60-100 nm range)
TARGET_LOG_G = 4.0        # log10(G_DN) = 4 → 10 kPa

# Constraint thresholds
MAX_SPAN = 2.0
MIN_G_DN = 1e3            # 1 kPa minimum


def compute_objectives(result: FullResult) -> np.ndarray:
    """Compute 3 objective values from a full pipeline result.

    f1 = |d32 - 2 µm| / 2 µm           (droplet size deviation)
    f2 = |pore - 80 nm| / 80 nm         (pore size deviation)
    f3 = |log10(G_DN) - 4.0|            (modulus deviation)

    All minimised.
    """
    d32 = result.emulsification.d32
    pore = result.gelation.pore_size_mean
    G_DN = max(result.mechanical.G_DN, 1.0)

    f1 = abs(d32 - TARGET_D32) / TARGET_D32
    f2 = abs(pore - TARGET_PORE) / TARGET_PORE
    f3 = abs(np.log10(G_DN) - TARGET_LOG_G)

    return np.array([f1, f2, f3])


def check_constraints(result: FullResult) -> tuple[bool, list[str]]:
    """Check optimisation constraints. Returns (feasible, violations)."""
    violations = []

    if result.emulsification.span > MAX_SPAN:
        violations.append(f"span={result.emulsification.span:.2f} > {MAX_SPAN}")

    if result.mechanical.G_DN < MIN_G_DN:
        violations.append(f"G_DN={result.mechanical.G_DN:.0f} Pa < {MIN_G_DN} Pa")

    return len(violations) == 0, violations


# Parameter bounds for optimisation (7-dimensional)
PARAM_NAMES = [
    "RPM", "c_span80", "agarose_frac", "T_oil",
    "cooling_rate", "c_genipin", "t_crosslink",
]

PARAM_BOUNDS = np.array([
    [3000.0,  25000.0],    # RPM
    [5.0,     50.0],       # c_span80 [kg/m³]
    [0.6,     0.9],        # agarose fraction
    [333.15,  368.15],     # T_oil [K] (60-95°C)
    [0.033,   0.333],      # cooling_rate [K/s] (2-20°C/min)
    [0.5,     10.0],       # c_genipin [mol/m³]
    [3600.0,  172800.0],   # t_crosslink [s] (1-48 h)
])

# Indices for log-scale parameters
LOG_SCALE_INDICES = [0, 1, 4, 5, 6]  # RPM, span80, cooling, genipin, t_xlink
