"""Chromatography breakthrough orchestrator with mass balance validation.

Architecture: docs/13_module2_module3_final_implementation_plan.md, Phase C.

Ties together hydrodynamics, isotherm, transport solver, and UV detection
into a single ``run_breakthrough`` call that returns a complete
BreakthroughResult including dynamic binding capacities.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

from .hydrodynamics import ColumnGeometry
from .isotherms.langmuir import LangmuirIsotherm
from .transport.lumped_rate import solve_lrm, LRMResult
from .detection.uv import compute_uv_signal, apply_detector_broadening

if TYPE_CHECKING:
    from emulsim.module2_functionalization.orchestrator import FunctionalMicrosphere

logger = logging.getLogger(__name__)


@dataclass
class BreakthroughResult:
    """Complete breakthrough simulation result.

    Attributes:
        time: Time array [s].
        uv_signal: UV absorbance signal [mAU].
        C_outlet: Outlet concentration [mol/m^3].
        dbc_5pct: Dynamic binding capacity at 5% breakthrough [mg/mL].
        dbc_10pct: Dynamic binding capacity at 10% breakthrough [mg/mL].
        dbc_50pct: Dynamic binding capacity at 50% breakthrough [mg/mL].
        pressure_drop: Column pressure drop [Pa].
        mass_balance_error: Relative mass balance error [-].
    """

    time: np.ndarray
    uv_signal: np.ndarray
    C_outlet: np.ndarray
    dbc_5pct: float
    dbc_10pct: float
    dbc_50pct: float
    pressure_drop: float
    mass_balance_error: float


def _compute_dbc(
    time: np.ndarray,
    C_outlet: np.ndarray,
    C_feed: float,
    flow_rate: float,
    column_volume: float,
    threshold_fraction: float,
) -> float:
    """Compute dynamic binding capacity at a given breakthrough threshold.

    DBC is defined as the mass of protein loaded per unit column volume
    at the time when the outlet concentration first reaches
    threshold_fraction * C_feed.

    Args:
        time: Time array [s].
        C_outlet: Outlet concentration [mol/m^3].
        C_feed: Feed concentration [mol/m^3].
        flow_rate: Volumetric flow rate [m^3/s].
        column_volume: Total bed volume [m^3].
        threshold_fraction: Breakthrough threshold (e.g., 0.05 for 5%).

    Returns:
        DBC [mol/m^3] (of column volume).
        Returns NaN if breakthrough is never reached.
    """
    C_threshold = threshold_fraction * C_feed

    # Find first time C_outlet >= C_threshold
    above = np.where(C_outlet >= C_threshold)[0]
    if len(above) == 0:
        # Breakthrough not reached — DBC is the total loaded
        t_bt = time[-1]
    else:
        idx = above[0]
        # Linear interpolation for more precise breakthrough time
        if idx > 0:
            C_lo, C_hi = C_outlet[idx - 1], C_outlet[idx]
            t_lo, t_hi = time[idx - 1], time[idx]
            if C_hi > C_lo:
                frac = (C_threshold - C_lo) / (C_hi - C_lo)
                t_bt = t_lo + frac * (t_hi - t_lo)
            else:
                t_bt = time[idx]
        else:
            t_bt = time[0]

    # Mass loaded up to breakthrough time (feed minus what passed through)
    mask = time <= t_bt
    if np.sum(mask) < 2:
        return 0.0

    t_bt_arr = time[mask]
    C_out_bt = C_outlet[mask]

    # Net mass captured = Q * integral(C_feed - C_outlet) dt
    mass_captured = flow_rate * float(np.trapezoid(C_feed - C_out_bt, t_bt_arr))
    mass_captured = max(mass_captured, 0.0)

    # DBC = mass_captured / column_volume [mol/m^3]
    return mass_captured / column_volume


def run_breakthrough(
    column: ColumnGeometry,
    microsphere=None,
    C_feed: float = 1.0,
    flow_rate: float = 1e-8,
    feed_duration: float = 600.0,
    total_time: float = 1200.0,
    extinction_coeff: float = 36000.0,
    sigma_detector: float = 1.0,
    isotherm: LangmuirIsotherm | None = None,
    n_z: int = 50,
    D_molecular: float = 7e-11,
    k_ads: float = 100.0,
) -> BreakthroughResult:
    """Run a single-component breakthrough simulation.

    If a FunctionalMicrosphere is provided, its mechanical properties
    (G_DN, E_star) and particle geometry override column defaults.

    Args:
        column: Column geometry and hydraulic parameters.
        microsphere: Optional FunctionalMicrosphere from Module 2.
        C_feed: Feed concentration [mol/m^3].
        flow_rate: Volumetric flow rate [m^3/s].
        feed_duration: Duration of feed/loading step [s].
        total_time: Total simulation time [s].
        extinction_coeff: Molar extinction coefficient [1/(M*cm)].
        sigma_detector: Detector broadening sigma [s].
        isotherm: Langmuir isotherm (default parameters if None).
        n_z: Number of axial finite-volume cells.
        D_molecular: Molecular diffusivity [m^2/s].
        k_ads: Adsorption rate constant [1/s].

    Returns:
        BreakthroughResult with DBC values, UV signal, and mass balance.
    """
    # ── Override column with microsphere properties if available ──
    if microsphere is not None:
        m1 = microsphere.m1_contract
        column = ColumnGeometry(
            diameter=column.diameter,
            bed_height=column.bed_height,
            particle_diameter=m1.bead_d50,
            bed_porosity=column.bed_porosity,
            particle_porosity=m1.porosity,
            G_DN=microsphere.G_DN_updated or m1.G_DN,
            E_star=microsphere.E_star_updated or m1.E_star,
        )

    # ── Validate flow rate ──
    flow_warnings = column.validate_flow_rate(flow_rate)
    for w in flow_warnings:
        logger.warning(w)

    # ── Default isotherm ──
    if isotherm is None:
        isotherm = LangmuirIsotherm()

    # ── Pressure drop ──
    pressure_drop = column.pressure_drop(flow_rate)

    # ── Solve LRM ──
    lrm_result = solve_lrm(
        column=column,
        isotherm=isotherm,
        C_feed=C_feed,
        feed_duration=feed_duration,
        flow_rate=flow_rate,
        total_time=total_time,
        n_z=n_z,
        D_molecular=D_molecular,
        k_ads=k_ads,
    )

    # ── UV detection ──
    uv_raw = compute_uv_signal(
        lrm_result.C_outlet,
        extinction_coeff=extinction_coeff,
    )
    uv_signal = apply_detector_broadening(uv_raw, lrm_result.time, sigma_detector)

    # ── Dynamic binding capacity ──
    V_col = column.bed_volume
    dbc_5 = _compute_dbc(
        lrm_result.time, lrm_result.C_outlet, C_feed, flow_rate, V_col, 0.05
    )
    dbc_10 = _compute_dbc(
        lrm_result.time, lrm_result.C_outlet, C_feed, flow_rate, V_col, 0.10
    )
    dbc_50 = _compute_dbc(
        lrm_result.time, lrm_result.C_outlet, C_feed, flow_rate, V_col, 0.50
    )

    return BreakthroughResult(
        time=lrm_result.time,
        uv_signal=uv_signal,
        C_outlet=lrm_result.C_outlet,
        dbc_5pct=dbc_5,
        dbc_10pct=dbc_10,
        dbc_50pct=dbc_50,
        pressure_drop=pressure_drop,
        mass_balance_error=lrm_result.mass_balance_error,
    )
