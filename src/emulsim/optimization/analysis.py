"""Pareto front analysis and visualisation utilities."""

from __future__ import annotations

import numpy as np

from ..datatypes import OptimizationState
from .objectives import PARAM_NAMES, PARAM_BOUNDS, LOG_SCALE_INDICES


def pareto_summary(state: OptimizationState) -> str:
    """Generate a text summary of the Pareto front."""
    lines = [
        f"Optimisation Summary",
        f"  Total evaluations: {len(state.X_observed)}",
        f"  Pareto-optimal points: {len(state.pareto_X)}",
        f"  Final hypervolume: {state.hypervolume:.4f}",
        f"  Converged: {state.converged}",
        "",
        "Pareto Front:",
        f"  {'#':>3s}  {'f_d32':>8s}  {'f_pore':>8s}  {'f_G_DN':>8s}  {'RPM':>8s}  {'Span80':>8s}  {'AgarFrac':>8s}  {'T_oil_C':>8s}  {'Cool_C/m':>8s}  {'Genipin':>8s}  {'t_xlink_h':>8s}",
    ]

    for i in range(len(state.pareto_X)):
        x_ss = state.pareto_X[i]
        y = state.pareto_Y[i]
        x = x_ss.copy()
        for j in LOG_SCALE_INDICES:
            x[j] = 10.0 ** x[j]

        lines.append(
            f"  {i+1:3d}"
            f"  {y[0]:8.3f}"
            f"  {y[1]:8.3f}"
            f"  {y[2]:8.3f}"
            f"  {x[0]:8.0f}"
            f"  {x[1]:8.1f}"
            f"  {x[2]:8.3f}"
            f"  {x[3]-273.15:8.1f}"
            f"  {x[4]*60:8.2f}"
            f"  {x[5]:8.2f}"
            f"  {x[6]/3600:8.1f}"
        )

    return "\n".join(lines)


def best_compromise(state: OptimizationState) -> int:
    """Find the best compromise Pareto point (min sum of objectives)."""
    if len(state.pareto_Y) == 0:
        return 0  # fallback
    sums = np.sum(state.pareto_Y, axis=1)
    return int(np.argmin(sums))
