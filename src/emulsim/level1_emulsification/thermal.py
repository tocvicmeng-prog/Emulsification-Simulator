"""Thermal model for emulsification heating strategies.

Computes temperature T(t) during emulsification for three strategies:
  - ISOTHERMAL: constant temperature (legacy mode)
  - FLAT_PLATE: Newton's law of cooling (lumped capacitance)
  - HOT_WATER_JACKET: steady-state during emulsification, with optional
    post-emulsification cooling profile

Protocol §1: All temperatures in Kelvin [K], times in seconds [s].
"""

from __future__ import annotations

import math

import numpy as np

from ..datatypes import HeatingConfig, VesselGeometry


# ─── Constants ───────────────────────────────────────────────────────────

_LN20 = math.log(20.0)  # ≈ 2.996; used for 95% cooling convention
_REFERENCE_VOLUME = 0.0005  # [m³] — 500 mL calibration reference

# ASSUMPTION: jacket insulation factor — jacket cooling time constant is
# 1.5× that of an equivalent flat-plate setup (jacket retains heat longer).
_JACKET_INSULATION_FACTOR = 1.5


# ─── Core Functions ──────────────────────────────────────────────────────

def cooling_time_constant(
    cooldown_time: float,
    working_volume: float,
    reference_volume: float = _REFERENCE_VOLUME,
) -> float:
    """Compute Newton cooling time constant tau [s].

    Protocol §2: At t = cooldown_time, 95% of the temperature difference
    has decayed.  This gives:
        exp(-cooldown_time / tau) = 0.05
        tau = cooldown_time / ln(20)

    The time constant scales linearly with volume relative to the
    calibration reference (500 mL).

    Parameters
    ----------
    cooldown_time : float
        Time [s] for 95% cooling at the reference volume.
    working_volume : float
        Actual working volume [m³].
    reference_volume : float
        Calibration reference volume [m³] (default 500 mL = 0.0005 m³).

    Returns
    -------
    float
        Time constant tau [s].
    """
    # Protocol §2: tau = cooldown_time / ln(20), scaled by volume ratio
    tau_ref = cooldown_time / _LN20
    return tau_ref * (working_volume / reference_volume)


def flat_plate_temperature(
    t: np.ndarray,
    T_initial: float,
    T_final: float,
    tau: float,
) -> np.ndarray:
    """Newton cooling: T(t) = T_final + (T_initial - T_final) * exp(-t/tau).

    Parameters
    ----------
    t : np.ndarray
        Time points [s].
    T_initial : float
        Starting temperature [K].
    T_final : float
        Ambient / asymptotic temperature [K].
    tau : float
        Cooling time constant [s].

    Returns
    -------
    np.ndarray
        Temperature [K] at each time point, same shape as *t*.
    """
    return T_final + (T_initial - T_final) * np.exp(-t / tau)


def jacket_cooling_profile(
    t: np.ndarray,
    T_initial: float,
    T_final: float,
    cooldown_time: float,
    working_volume: float,
    reference_volume: float = _REFERENCE_VOLUME,
) -> np.ndarray:
    """Post-emulsification cooling profile for the hot-water jacket.

    After the jacket circulating pump is switched off, the vessel cools
    following the same Newton cooling model as flat-plate heating, but
    with a longer time constant (jacket insulates better).

    ASSUMPTION: jacket insulation factor = 1.5x flat-plate tau.

    Parameters
    ----------
    t : np.ndarray
        Time points [s] measured from the moment the jacket is turned off.
    T_initial : float
        Temperature [K] at the moment the jacket is turned off.
    T_final : float
        Ambient temperature [K].
    cooldown_time : float
        Calibration cooldown time [s] for 95% cooling at reference volume.
    working_volume : float
        Actual working volume [m³].
    reference_volume : float
        Calibration reference volume [m³] (default 500 mL).

    Returns
    -------
    np.ndarray
        Temperature [K] at each time point.
    """
    tau = cooling_time_constant(cooldown_time, working_volume, reference_volume)
    tau_jacket = tau * _JACKET_INSULATION_FACTOR
    return flat_plate_temperature(t, T_initial, T_final, tau_jacket)


# ─── Main Dispatcher ────────────────────────────────────────────────────

def temperature_profile(
    heating: HeatingConfig,
    vessel: VesselGeometry,
    t_array: np.ndarray,
) -> np.ndarray:
    """Compute T(t) [K] during emulsification.

    Dispatches to the appropriate thermal model based on the heating
    strategy specified in *heating*.

    Protocol §1-3: Supports ISOTHERMAL, FLAT_PLATE, and HOT_WATER_JACKET.

    Parameters
    ----------
    heating : HeatingConfig
        Heating/cooling configuration (strategy, temperatures, timings).
    vessel : VesselGeometry
        Vessel geometry (used for volume-dependent scaling).
    t_array : np.ndarray
        Time points [s] at which to evaluate temperature.

    Returns
    -------
    np.ndarray
        Temperature [K] at each time point, same shape as *t_array*.

    Raises
    ------
    ValueError
        If the heating strategy is not recognised.
    """
    # Compare by .value (string) to avoid Streamlit module-reload identity
    # mismatch where two HeatingStrategy enums from different import cycles
    # have the same value but fail == comparison.
    strategy_val = heating.strategy.value if hasattr(heating.strategy, 'value') else str(heating.strategy)

    # Protocol §1: ISOTHERMAL — constant temperature
    if strategy_val == "isothermal":
        return np.full_like(t_array, heating.T_initial, dtype=np.float64)

    # Protocol §2: FLAT_PLATE — Newton cooling (lumped capacitance)
    if strategy_val == "flat_plate":
        tau = cooling_time_constant(
            heating.cooldown_time,
            vessel.working_volume,
        )
        return flat_plate_temperature(
            t_array, heating.T_initial, heating.T_final, tau,
        )

    # Protocol §3: HOT_WATER_JACKET — steady state during emulsification
    # ASSUMPTION: jacket maintains T_initial throughout the emulsification
    # period via continuous hot-water circulation.
    if strategy_val == "hot_water_jacket":
        return np.full_like(t_array, heating.T_initial, dtype=np.float64)

    raise ValueError(f"Unknown heating strategy: {heating.strategy!r}")
