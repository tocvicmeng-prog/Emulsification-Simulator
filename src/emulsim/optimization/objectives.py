"""Objective function definitions for multi-objective optimization.

Maps simulation FullResult → objective vector for BoTorch.
All objectives are MINIMIZED (lower is better).

Supports two target sets:
  - Legacy rotor-stator: d32=2 µm, pore=80 nm
  - Stirred-vessel: d_mode=100 µm, pore=100 nm
"""

from __future__ import annotations

import numpy as np

from ..datatypes import FullResult


# ─── Legacy rotor-stator targets ─────────────────────────────────────────

TARGET_D32 = 2.0e-6      # 2 µm
TARGET_PORE = 80.0e-9    # 80 nm (centre of 60-100 nm range)
TARGET_LOG_G = 4.0        # log10(G_DN) = 4 → 10 kPa

# ─── Stirred-vessel targets ──────────────────────────────────────────────
# d_mode (volume-weighted modal diameter) is used instead of d32 for
# stirred-vessel mode because: (a) downstream gelation/crosslinking acts
# on individual droplets, so the most-probable size matters more than
# the surface-area-averaged d32; (b) d32 is dominated by fines in broad
# distributions, which is misleading for the stirred-vessel regime where
# the target is 75-150 um microspheres.

TARGET_D_MODE = 100.0e-6  # 100 µm modal peak
TARGET_PORE_SV = 100.0e-9 # 100 nm (centre of 80-120 nm range)

# ─── Constraint thresholds ───────────────────────────────────────────────

MAX_SPAN = 2.0
MIN_G_DN = 1e3            # 1 kPa minimum


def compute_objectives(result: FullResult, mode: str | None = None) -> np.ndarray:
    """Compute 3 objective values from a full pipeline result.

    Legacy mode:
        f1 = |d32 - 2 µm| / 2 µm
        f2 = |pore - 80 nm| / 80 nm
        f3 = |log10(G_DN) - 4.0|

    Stirred-vessel mode:
        f1 = |d_mode - 100 µm| / 100 µm
        f2 = |pore - 100 nm| / 100 nm
        f3 = |log10(G_DN) - 4.0|

    All minimised.

    If ``mode`` is None, it is inferred from the result's parameters.
    """
    if mode is None:
        mode = getattr(result.parameters.emulsification, 'mode', 'rotor_stator_legacy')

    pore = result.gelation.pore_size_mean
    G_DN = max(result.mechanical.G_DN, 1.0)

    if mode == "stirred_vessel":
        d_mode = getattr(result.emulsification, 'd_mode', 0.0)
        f1 = abs(d_mode - TARGET_D_MODE) / TARGET_D_MODE
        f2 = abs(pore - TARGET_PORE_SV) / TARGET_PORE_SV
    else:
        d32 = result.emulsification.d32
        f1 = abs(d32 - TARGET_D32) / TARGET_D32
        f2 = abs(pore - TARGET_PORE) / TARGET_PORE

    f3 = abs(np.log10(G_DN) - TARGET_LOG_G)

    return np.array([f1, f2, f3])


def check_constraints(result: FullResult, mode: str | None = None) -> tuple[bool, list[str]]:
    """Check optimisation constraints. Returns (feasible, violations)."""
    if mode is None:
        mode = getattr(result.parameters.emulsification, 'mode', 'rotor_stator_legacy')

    violations = []

    if result.emulsification.span > MAX_SPAN:
        violations.append(f"span={result.emulsification.span:.2f} > {MAX_SPAN}")

    if result.mechanical.G_DN < MIN_G_DN:
        violations.append(f"G_DN={result.mechanical.G_DN:.0f} Pa < {MIN_G_DN} Pa")

    # Stirred-vessel: check d_mode is in ideal range [75, 150] µm
    if mode == "stirred_vessel":
        d_mode = getattr(result.emulsification, 'd_mode', 0.0)
        if d_mode < 75e-6 or d_mode > 150e-6:
            violations.append(
                f"d_mode={d_mode*1e6:.1f} µm outside ideal range [75, 150] µm"
            )

    return len(violations) == 0, violations


# ─── Parameter bounds (legacy — rotor-stator) ────────────────────────────

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

# ─── Parameter bounds (stirred-vessel) ───────────────────────────────────

PARAM_NAMES_SV = [
    "RPM", "c_span80", "agarose_frac", "T_oil",
    "cooling_rate", "c_genipin", "t_crosslink",
]

PARAM_BOUNDS_SV = np.array([
    [800.0,   2000.0],     # RPM (Stirrer A: pitched blade)
    [5.0,     30.0],       # c_span80 [kg/m³] (1.5% v/v default)
    [0.6,     0.9],        # agarose fraction
    [333.15,  358.15],     # T_oil [K] (60-85°C, limited by heating)
    [0.005,   0.05],       # cooling_rate [K/s] (0.3-3°C/min, natural cooling)
    [0.5,     10.0],       # c_genipin [mol/m³]
    [3600.0,  172800.0],   # t_crosslink [s] (1-48 h)
])

PARAM_BOUNDS_SV_RS = np.array([
    [800.0,   9000.0],     # RPM (Stirrer B: small rotor-stator)
    [5.0,     30.0],       # c_span80 [kg/m³]
    [0.6,     0.9],        # agarose fraction
    [333.15,  358.15],     # T_oil [K]
    [0.005,   0.05],       # cooling_rate [K/s]
    [0.5,     10.0],       # c_genipin [mol/m³]
    [3600.0,  172800.0],   # t_crosslink [s]
])

LOG_SCALE_INDICES_SV = [0, 1, 4, 5, 6]


def get_param_bounds(mode: str = "rotor_stator_legacy",
                     stirrer_type: str = "pitched_blade") -> np.ndarray:
    """Return parameter bounds appropriate for the given mode and stirrer."""
    if mode == "stirred_vessel":
        if stirrer_type == "rotor_stator_small":
            return PARAM_BOUNDS_SV_RS.copy()
        return PARAM_BOUNDS_SV.copy()
    return PARAM_BOUNDS.copy()
