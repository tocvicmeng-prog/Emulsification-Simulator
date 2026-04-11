"""Gradient generator for chromatographic elution programs.

Architecture: docs/13_module2_module3_final_implementation_plan.md, Phase E.

Supports three gradient forms:
  - Linear: single ramp from start to end value.
  - Step:   instantaneous transitions between plateaus (modelled as 1 s ramps).
  - Multi-linear (GradientProgram): arbitrary piecewise-linear profile
    defined by (time, value) waypoints.

All gradient objects return a scalar value at any query time via
``value_at_time(t)``, which uses ``np.interp`` with constant extrapolation
at the boundaries.  This makes them safe to call inside ODE integrators.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np


# ─── Core Gradient Program ────────────────────────────────────────────────────

@dataclass
class GradientProgram:
    """Piecewise-linear gradient program defined by (time, value) waypoints.

    The gradient is linearly interpolated between waypoints and held constant
    outside the defined range (``np.interp`` boundary behaviour).

    Attributes:
        segments: List of (time_s, value) pairs defining the gradient profile.
            Must be sorted in ascending time order.

    Examples:
        A linear 0->1 M NaCl gradient over 600 s followed by a 300 s wash::

            g = GradientProgram(segments=[(0, 0.0), (600, 1.0), (900, 1.0)])
            g.value_at_time(300)   # 0.5
            g.value_at_time(750)   # 1.0
    """

    segments: list[tuple[float, float]] = field(
        default_factory=lambda: [(0.0, 0.0), (600.0, 1.0)]
    )

    def __post_init__(self) -> None:
        if len(self.segments) < 2:
            raise ValueError(
                "GradientProgram requires at least 2 waypoints (start, end)."
            )
        times = [s[0] for s in self.segments]
        if times != sorted(times):
            raise ValueError(
                "GradientProgram segments must be sorted in ascending time order."
            )

    @property
    def _times(self) -> np.ndarray:
        return np.array([s[0] for s in self.segments])

    @property
    def _values(self) -> np.ndarray:
        return np.array([s[1] for s in self.segments])

    def value_at_time(self, t: float | np.ndarray) -> float | np.ndarray:
        """Evaluate the gradient at time t by piecewise-linear interpolation.

        Extrapolates as constant outside the defined range.

        Args:
            t: Time [s], scalar or array.

        Returns:
            Gradient value (salt concentration, pH, etc.) at time t,
            same type as input.
        """
        scalar = np.isscalar(t)
        t_arr = np.asarray(t, dtype=float)
        result = np.interp(t_arr, self._times, self._values)
        return float(result) if scalar else result

    def total_duration(self) -> float:
        """Total gradient duration [s]."""
        return float(self._times[-1])

    def start_value(self) -> float:
        """Gradient value at t=0."""
        return float(self._values[0])

    def end_value(self) -> float:
        """Gradient value at end of program."""
        return float(self._values[-1])


# ─── Convenience Constructors ─────────────────────────────────────────────────

def make_linear_gradient(
    start_val: float,
    end_val: float,
    start_time: float = 0.0,
    end_time: float = 600.0,
) -> GradientProgram:
    """Create a simple linear gradient between two values.

    The gradient is constant (at start_val) before start_time and
    constant (at end_val) after end_time.

    Args:
        start_val: Initial gradient value (e.g., salt concentration [mol/m^3]).
        end_val: Final gradient value.
        start_time: Time at which the ramp begins [s].
        end_time: Time at which the ramp reaches end_val [s].

    Returns:
        GradientProgram for a simple linear ramp.

    Raises:
        ValueError: If end_time <= start_time.

    Examples:
        >>> g = make_linear_gradient(0.0, 1.0, end_time=300.0)
        >>> g.value_at_time(150.0)
        0.5
    """
    if end_time <= start_time:
        raise ValueError(
            f"end_time ({end_time}) must be greater than start_time ({start_time})."
        )
    segments: list[tuple[float, float]] = []
    if start_time > 0.0:
        segments.append((0.0, start_val))
    segments.append((start_time, start_val))
    segments.append((end_time, end_val))
    # Deduplicate if start_time == 0
    seen_times: set[float] = set()
    unique_segments: list[tuple[float, float]] = []
    for seg in segments:
        if seg[0] not in seen_times:
            unique_segments.append(seg)
            seen_times.add(seg[0])
    return GradientProgram(segments=unique_segments)


def make_step_gradient(
    steps: list[tuple[float, float]],
    ramp_duration: float = 1.0,
) -> GradientProgram:
    """Create a step gradient with short linear ramps between plateaus.

    Each step is specified as (time_s, value) where the gradient jumps
    to ``value`` at ``time_s``.  To avoid discontinuities (which can cause
    issues in ODE solvers), transitions are modelled as short linear ramps
    of duration ``ramp_duration`` seconds.

    Args:
        steps: List of (time_s, value) tuples defining each plateau.
            The first tuple sets the initial condition at t=0.
            Subsequent tuples define transitions.
        ramp_duration: Duration of each step transition ramp [s].

    Returns:
        GradientProgram with ramp-approximated steps.

    Raises:
        ValueError: If fewer than 2 steps are provided.

    Examples:
        A two-step gradient: 0 M for 0-300 s, then 0.5 M, then 1 M::

            g = make_step_gradient([(0, 0.0), (300, 0.5), (600, 1.0)])
            g.value_at_time(310)   # ~0.5
    """
    if len(steps) < 2:
        raise ValueError("make_step_gradient requires at least 2 steps.")

    # Sort by time
    steps_sorted = sorted(steps, key=lambda x: x[0])

    segments: list[tuple[float, float]] = []
    for i, (t, val) in enumerate(steps_sorted):
        if i == 0:
            segments.append((t, val))
        else:
            prev_val = steps_sorted[i - 1][1]
            # Add a point just before the ramp (at previous value)
            t_ramp_start = t
            t_ramp_end = t + ramp_duration
            # Hold previous value until ramp start
            segments.append((t_ramp_start, prev_val))
            # Ramp to new value
            segments.append((t_ramp_end, val))

    return GradientProgram(segments=segments)


def make_wash_gradient(
    load_value: float,
    wash_value: float,
    elute_value: float,
    load_duration: float = 300.0,
    wash_duration: float = 120.0,
    elute_duration: float = 300.0,
    strip_value: float | None = None,
    strip_duration: float = 60.0,
    ramp_duration: float = 1.0,
) -> GradientProgram:
    """Create a standard load-wash-elute-strip gradient program.

    This covers the typical 4-step affinity/IEX cycle:
      1. Load (isocratic at load_value)
      2. Wash  (isocratic at wash_value)
      3. Elute (isocratic at elute_value)
      4. Strip (optional, isocratic at strip_value)

    Args:
        load_value: Gradient value during loading.
        wash_value: Gradient value during wash.
        elute_value: Gradient value for elution.
        load_duration: Duration of load step [s].
        wash_duration: Duration of wash step [s].
        elute_duration: Duration of elution step [s].
        strip_value: Gradient value for strip step (None = omit strip).
        strip_duration: Duration of strip step [s].
        ramp_duration: Transition ramp duration [s].

    Returns:
        GradientProgram for load-wash-elute cycle.
    """
    t0 = 0.0
    t_wash = t0 + load_duration
    t_elute = t_wash + wash_duration
    t_end = t_elute + elute_duration

    steps: list[tuple[float, float]] = [
        (t0, load_value),
        (t_wash, wash_value),
        (t_elute, elute_value),
    ]

    if strip_value is not None:
        steps.append((t_end, strip_value))

    return make_step_gradient(steps, ramp_duration=ramp_duration)
