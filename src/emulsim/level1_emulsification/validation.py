"""Validation utilities for the emulsification simulation."""

from __future__ import annotations

import numpy as np

from ..datatypes import EmulsificationResult


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
