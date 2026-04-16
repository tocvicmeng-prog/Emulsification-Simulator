"""Validation utilities for the emulsification simulation.

Includes physical-bounds checks, steady-state checks, mass conservation,
and v6.1 dimensionless group calculations for domain validation.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

import numpy as np

from ..datatypes import EmulsificationResult


# ── Mode-specific constants ──────────────────────────────────────────────

# ASSUMPTION: rotor-stator targets d32 ~ 2 µm; stirred-vessel targets d_mode ~ 100 µm
_BOUNDS = {
    "rotor_stator_legacy": {
        "d32_min": 0.05e-6,       # 0.05 µm
        "d32_max": 1000e-6,       # 1000 µm
        "steady_state_tol": 0.01, # 1 % relative variation
        "steady_state_frac": 0.10,  # last 10 % of time history
    },
    "stirred_vessel": {
        "d32_min": 5e-6,          # 5 µm
        "d32_max": 2000e-6,       # 2 mm
        "steady_state_tol": 0.05, # 5 % — larger, noisier droplets
        "steady_state_frac": 0.20,  # last 20 % of time history
        # ASSUMPTION: modal diameter valid range for stirred-vessel
        "d_mode_min": 5e-6,       # 5 µm
        "d_mode_max": 500e-6,     # 500 µm
    },
}


def _get_bounds(mode: Optional[str] = None) -> dict:
    """Return physical bounds dict for the given mode.

    Falls back to rotor_stator_legacy when *mode* is ``None`` or unrecognised,
    preserving backward compatibility.
    """
    if mode is None or mode not in _BOUNDS:
        mode = "rotor_stator_legacy"
    return _BOUNDS[mode]


# ── Existing functions (signatures unchanged) ────────────────────────────

def check_mass_conservation(result: EmulsificationResult,
                            initial_phi_d: float = 0.05,
                            tol: float = 1e-3) -> tuple[bool, float]:
    """Check that total droplet volume is conserved.

    Parameters
    ----------
    result : EmulsificationResult
        Simulation result.
    initial_phi_d : float
        Expected total volume fraction (initial phi_d).
    tol : float
        Relative tolerance.

    Returns
    -------
    passed : bool
    relative_error : float
    """
    rel_err = abs(result.total_volume_fraction - initial_phi_d) / max(initial_phi_d, 1e-30)
    return rel_err < tol, rel_err


def check_non_negative(result: EmulsificationResult) -> bool:
    """Check all number densities are non-negative."""
    return np.all(result.n_d >= 0)


def check_physical_bounds(result: EmulsificationResult,
                          d32_min: float = 0.05e-6,
                          d32_max: float = 1000e-6) -> tuple[bool, str]:
    """Check that d32 and other statistics are physically reasonable.

    Returns (passed, message).
    """
    issues = []

    if result.d32 < d32_min:
        issues.append(f"d32={result.d32*1e6:.3f} µm below minimum {d32_min*1e6:.1f} µm")
    if result.d32 > d32_max:
        issues.append(f"d32={result.d32*1e6:.1f} µm above maximum {d32_max*1e6:.0f} µm")
    if result.span < 0:
        issues.append(f"span={result.span:.2f} is negative")
    if result.d10 > result.d50 or result.d50 > result.d90:
        issues.append("percentiles not monotonic: d10 > d50 or d50 > d90")

    if issues:
        return False, "; ".join(issues)
    return True, "OK"


def check_steady_state(result: EmulsificationResult,
                       tol: float = 0.01) -> tuple[bool, float]:
    """Check if the simulation reached steady state.

    Looks at d32 variation over the last 10% of time history.

    Returns (converged, relative_variation).
    """
    if result.n_d_history is None or len(result.n_d_history) < 2:
        return False, float('inf')

    # Compute d32 at last 10% of time steps
    n_steps = len(result.n_d_history)
    n_check = max(1, n_steps // 10)
    d32_vals = []
    for k in range(-n_check, 0):
        n_d_k = np.maximum(result.n_d_history[k] * (result.d_bins[1:] - result.d_bins[:-1]
                           if len(result.d_bins) > len(result.n_d_history[k]) else
                           np.ones_like(result.n_d_history[k])), 0.0)
        # Use raw n_d_history (already number density) to compute d32
        N_k = np.maximum(result.n_d_history[k], 0.0)
        d3 = np.sum(N_k * result.d_bins**3)
        d2 = np.sum(N_k * result.d_bins**2)
        d32_k = d3 / d2 if d2 > 0 else 0.0
        d32_vals.append(d32_k)

    if not d32_vals or max(d32_vals) == 0:
        return False, float('inf')

    variation = (max(d32_vals) - min(d32_vals)) / max(max(d32_vals), 1e-15)
    return variation < tol, variation


def validate_result(result: EmulsificationResult) -> dict:
    """Run all validation checks on an emulsification result.

    Returns dict with check names as keys and (passed, detail) as values.
    """
    checks = {}
    checks['non_negative'] = (check_non_negative(result), "")
    checks['physical_bounds'] = check_physical_bounds(result)
    checks['steady_state'] = check_steady_state(result)
    return checks


# ── Mode-aware validation functions ──────────────────────────────────────

def check_physical_bounds_mode(result: EmulsificationResult,
                               mode: str = "rotor_stator_legacy") -> tuple[bool, str]:
    """Mode-aware physical-bounds check.

    Uses pre-calibrated d32 ranges from ``_BOUNDS`` for the requested mode.
    For ``stirred_vessel`` mode, also validates the modal diameter (d_mode).

    Parameters
    ----------
    result : EmulsificationResult
        Simulation result.
    mode : str
        ``"rotor_stator_legacy"`` or ``"stirred_vessel"``.

    Returns
    -------
    passed : bool
    message : str
    """
    b = _get_bounds(mode)
    # Delegate d32 + basic checks to the original function with mode-specific limits
    passed, msg = check_physical_bounds(result,
                                        d32_min=b["d32_min"],
                                        d32_max=b["d32_max"])
    issues = [] if passed else [msg]

    # Additional d_mode check for stirred-vessel mode
    if mode == "stirred_vessel":
        d_mode = compute_d_mode(result)
        d_mode_min = b["d_mode_min"]
        d_mode_max = b["d_mode_max"]
        if d_mode < d_mode_min:
            issues.append(
                f"d_mode={d_mode*1e6:.1f} µm below minimum "
                f"{d_mode_min*1e6:.0f} µm")
        if d_mode > d_mode_max:
            issues.append(
                f"d_mode={d_mode*1e6:.1f} µm above maximum "
                f"{d_mode_max*1e6:.0f} µm")

    if issues:
        return False, "; ".join(issues)
    return True, "OK"


def check_steady_state_mode(result: EmulsificationResult,
                            mode: str = "rotor_stator_legacy") -> tuple[bool, float]:
    """Mode-aware steady-state check.

    For ``stirred_vessel`` mode the tolerance is relaxed to 5 % and the
    window is expanded to the last 20 % of time history (larger, noisier
    droplets converge more slowly).

    Parameters
    ----------
    result : EmulsificationResult
        Simulation result.
    mode : str
        ``"rotor_stator_legacy"`` or ``"stirred_vessel"``.

    Returns
    -------
    converged : bool
    relative_variation : float
    """
    b = _get_bounds(mode)
    tol = b["steady_state_tol"]
    frac = b["steady_state_frac"]

    if result.n_d_history is None or len(result.n_d_history) < 2:
        return False, float('inf')

    n_steps = len(result.n_d_history)
    n_check = max(1, int(n_steps * frac))
    d32_vals = []
    for k in range(-n_check, 0):
        N_k = np.maximum(result.n_d_history[k], 0.0)
        d3 = np.sum(N_k * result.d_bins**3)
        d2 = np.sum(N_k * result.d_bins**2)
        d32_k = d3 / d2 if d2 > 0 else 0.0
        d32_vals.append(d32_k)

    if not d32_vals or max(d32_vals) == 0:
        return False, float('inf')

    variation = (max(d32_vals) - min(d32_vals)) / max(max(d32_vals), 1e-15)
    return variation < tol, variation


def compute_d_mode(result: EmulsificationResult) -> float:
    """Compute the modal diameter (peak of volume-weighted distribution).

    ASSUMPTION: volume-weighted distribution is proportional to
    n_d * d^3; the mode is the bin with the highest volume density.

    Parameters
    ----------
    result : EmulsificationResult

    Returns
    -------
    d_mode : float
        Modal diameter [m].  Returns 0.0 if the distribution is empty.
    """
    vol_weighted = result.n_d * result.d_bins**3
    if np.sum(vol_weighted) == 0:
        return 0.0
    idx = np.argmax(vol_weighted)
    return float(result.d_bins[idx])


def check_d_mode(result: EmulsificationResult,
                 d_mode_min: float = 5e-6,
                 d_mode_max: float = 500e-6) -> tuple[bool, str]:
    """Standalone check that the modal diameter is within bounds.

    Primarily intended for stirred-vessel mode.

    Parameters
    ----------
    result : EmulsificationResult
    d_mode_min : float
        Minimum acceptable modal diameter [m].
    d_mode_max : float
        Maximum acceptable modal diameter [m].

    Returns
    -------
    passed : bool
    message : str
    """
    d_mode = compute_d_mode(result)
    issues = []
    if d_mode < d_mode_min:
        issues.append(
            f"d_mode={d_mode*1e6:.1f} µm below minimum {d_mode_min*1e6:.0f} µm")
    if d_mode > d_mode_max:
        issues.append(
            f"d_mode={d_mode*1e6:.1f} µm above maximum {d_mode_max*1e6:.0f} µm")
    if issues:
        return False, "; ".join(issues)
    return True, "OK"


def validate_result_mode(result: EmulsificationResult,
                         mode: str = "rotor_stator_legacy") -> dict:
    """Mode-aware validation suite.

    Runs all checks from :func:`validate_result` using mode-appropriate
    thresholds.  For ``stirred_vessel`` mode, an additional ``d_mode``
    check is included.

    Parameters
    ----------
    result : EmulsificationResult
    mode : str
        ``"rotor_stator_legacy"`` or ``"stirred_vessel"``.

    Returns
    -------
    dict
        Check names as keys, ``(passed, detail)`` as values.
    """
    checks: dict = {}
    checks['non_negative'] = (check_non_negative(result), "")
    checks['physical_bounds'] = check_physical_bounds_mode(result, mode=mode)
    checks['steady_state'] = check_steady_state_mode(result, mode=mode)

    if mode == "stirred_vessel":
        checks['d_mode'] = check_d_mode(result)

    return checks


# ── Dimensionless Groups (v6.1) ──────────────────────────────────────────

@dataclass
class DimensionlessGroups:
    """Dimensionless groups characterizing the emulsification regime.

    Used for domain validation and evidence-tier assignment.
    """
    Re: float = 0.0           # Reynolds number: rho_c * N * D^2 / mu_c
    We: float = 0.0           # Weber number: rho_c * N^2 * D^3 / sigma
    Ca: float = 0.0           # Capillary number: mu_c * v_tip / sigma
    viscosity_ratio: float = 0.0   # lambda = mu_d / mu_c
    eta_K: float = 0.0        # Kolmogorov length [m]: (nu^3 / epsilon)^0.25
    d32_over_eta_K: float = 0.0    # droplet/Kolmogorov ratio
    Np: float = 0.0           # Power number (from equipment)

    def as_dict(self) -> dict:
        """Return all groups as a plain dict for diagnostics/JSON."""
        return {
            "Re": self.Re,
            "We": self.We,
            "Ca": self.Ca,
            "viscosity_ratio": self.viscosity_ratio,
            "eta_K_m": self.eta_K,
            "d32_over_eta_K": self.d32_over_eta_K,
            "Np": self.Np,
        }


def compute_dimensionless_groups(
    rpm: float,
    D_impeller: float,
    rho_c: float,
    mu_c: float,
    sigma: float,
    mu_d: float,
    Np: float,
    V_tank: float,
    d32: float = 0.0,
) -> DimensionlessGroups:
    """Compute dimensionless groups for an emulsification system.

    Parameters
    ----------
    rpm : float
        Impeller speed [rev/min].
    D_impeller : float
        Impeller diameter [m].
    rho_c : float
        Continuous phase density [kg/m^3].
    mu_c : float
        Continuous phase viscosity [Pa.s].
    sigma : float
        Interfacial tension [N/m].
    mu_d : float
        Dispersed phase viscosity [Pa.s].
    Np : float
        Power number [-].
    V_tank : float
        Tank volume [m^3].
    d32 : float, optional
        Sauter mean diameter [m] (0 if not yet computed).

    Returns
    -------
    DimensionlessGroups
    """
    N = rpm / 60.0  # rev/s
    if N <= 0 or D_impeller <= 0 or rho_c <= 0:
        return DimensionlessGroups()

    nu_c = mu_c / rho_c if rho_c > 0 else 1e-6  # kinematic viscosity [m^2/s]
    v_tip = np.pi * D_impeller * N  # tip speed [m/s]

    # Reynolds number
    Re = rho_c * N * D_impeller**2 / mu_c if mu_c > 0 else 0.0

    # Weber number
    We = rho_c * N**2 * D_impeller**3 / sigma if sigma > 0 else 0.0

    # Capillary number
    Ca = mu_c * v_tip / sigma if sigma > 0 else 0.0

    # Viscosity ratio
    viscosity_ratio = mu_d / mu_c if mu_c > 0 else 0.0

    # Mean energy dissipation rate
    P = Np * rho_c * N**3 * D_impeller**5  # power [W]
    epsilon = P / (rho_c * V_tank) if V_tank > 0 and rho_c > 0 else 0.0

    # Kolmogorov length
    eta_K = (nu_c**3 / epsilon)**0.25 if epsilon > 0 else 0.0

    # d32 / Kolmogorov ratio
    d32_over_eta_K = d32 / eta_K if eta_K > 0 and d32 > 0 else 0.0

    return DimensionlessGroups(
        Re=Re,
        We=We,
        Ca=Ca,
        viscosity_ratio=viscosity_ratio,
        eta_K=eta_K,
        d32_over_eta_K=d32_over_eta_K,
        Np=Np,
    )
