"""Node F2 Phase 1: trace-replay harness for the digital-twin EnKF.

Walks forward through a :class:`DigitalTwinTrace`, applying the
user-supplied ``state_transition`` between observation times and the
EnKF ``enkf_update`` whenever an observation arrives. Snapshots
(mean, std, full ensemble) are recorded at each observation time for
downstream plotting or validation.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, Optional

import numpy as np

from ..datatypes import ModelEvidenceTier, ModelManifest
from .enkf import enkf_update
from .schema import DigitalTwinTrace, Observation


@dataclass
class ReplayResult:
    """Outputs of :func:`run_replay`.

    Attributes
    ----------
    times : list[float]
        Observation times [s], in order.
    means : list[np.ndarray]
        Posterior ensemble mean at each observation time; one
        ``(n_state,)`` array per entry.
    stds : list[np.ndarray]
        Posterior ensemble standard deviation at each observation
        time; one ``(n_state,)`` array per entry.
    ensembles : list[np.ndarray]
        Full ``(N, n_state)`` posterior ensemble at each observation
        time. Stored only if ``store_ensembles=True`` was passed.
    final_ensemble : np.ndarray
        The posterior ensemble at the last observation (always
        retained for chaining into further analysis).
    model_manifest : ModelManifest
        Provenance record (SEMI_QUANTITATIVE).
    """

    times: list[float]
    means: list[np.ndarray]
    stds: list[np.ndarray]
    ensembles: list[np.ndarray]
    final_ensemble: np.ndarray
    model_manifest: ModelManifest = field(
        default_factory=lambda: ModelManifest(
            model_name="DigitalTwin.EnKFReplay",
            evidence_tier=ModelEvidenceTier.SEMI_QUANTITATIVE,
            assumptions=[],
            diagnostics={},
        )
    )


def run_replay(
    trace: DigitalTwinTrace,
    initial_ensemble: np.ndarray,
    state_transition: Callable[[np.ndarray, float], np.ndarray],
    observation_operator: Callable[[np.ndarray, str], np.ndarray],
    inflation: float = 1.02,
    seed: Optional[int] = 0,
    store_ensembles: bool = False,
) -> ReplayResult:
    """Sequentially replay a trace through a stochastic EnKF.

    Parameters
    ----------
    trace : DigitalTwinTrace
        Time-ordered observations + metadata.
    initial_ensemble : np.ndarray, shape (N, n_state)
        Prior ensemble at ``t = 0``.
    state_transition : callable
        ``f(x_ensemble, dt) -> x_ensemble`` — the physical model
        advance between observation times. ``dt`` is in seconds. May
        be the identity for a persistence / "stationary hidden
        parameter" model.
    observation_operator : callable
        ``h(x_ensemble, name) -> y_forecast`` — maps each ensemble
        member to a scalar forecast of the observable identified by
        ``name``. Must return shape ``(N,)``.
    inflation : float
        Multiplicative prior inflation passed to ``enkf_update``.
    seed : int, optional
        RNG seed for EnKF stochastic perturbations.
    store_ensembles : bool
        If True, retain the full ``(N, n_state)`` ensemble at every
        observation (memory-heavy for long traces with big N). If
        False, only the final ensemble is kept.

    Returns
    -------
    ReplayResult

    Raises
    ------
    ValueError
        If ``initial_ensemble`` is not 2-D or N < 2.
    """
    x_ens = np.asarray(initial_ensemble, dtype=float).copy()
    if x_ens.ndim != 2:
        raise ValueError(
            f"initial_ensemble must be 2-D; got shape {x_ens.shape}"
        )
    N, n_state = x_ens.shape
    if N < 2:
        raise ValueError(f"EnKF requires N >= 2, got {N}")

    rng = np.random.default_rng(seed)

    times: list[float] = []
    means: list[np.ndarray] = []
    stds: list[np.ndarray] = []
    ensembles: list[np.ndarray] = []

    t_prev = 0.0
    n_updates = 0
    for obs in trace.observations:
        dt = float(obs.t) - t_prev
        if dt < 0.0:
            raise ValueError(
                "DigitalTwinTrace observations must be time-ordered; "
                f"encountered dt = {dt} at obs {obs!r}"
            )
        # Advance the state ensemble
        if dt > 0.0:
            x_ens = state_transition(x_ens, dt)
            x_ens = np.asarray(x_ens, dtype=float)
            if x_ens.shape != (N, n_state):
                raise ValueError(
                    f"state_transition returned shape {x_ens.shape}; "
                    f"expected {(N, n_state)}"
                )

        # Forecast the observable
        y_fc = observation_operator(x_ens, obs.name)
        y_fc = np.asarray(y_fc, dtype=float).ravel()
        if y_fc.size != N:
            raise ValueError(
                f"observation_operator returned size {y_fc.size}; "
                f"expected ensemble size {N}"
            )

        # EnKF update
        x_ens = enkf_update(
            x_ens, y_fc, float(obs.value),
            R=float(obs.noise_std) ** 2,
            rng=rng, inflation=inflation,
        )
        n_updates += 1

        # Record snapshot
        times.append(float(obs.t))
        means.append(x_ens.mean(axis=0))
        stds.append(x_ens.std(axis=0, ddof=0))
        if store_ensembles:
            ensembles.append(x_ens.copy())

        t_prev = float(obs.t)

    manifest = ModelManifest(
        model_name="DigitalTwin.EnKFReplay",
        evidence_tier=ModelEvidenceTier.SEMI_QUANTITATIVE,
        assumptions=[
            "stochastic EnKF (Evensen 1994)",
            "scalar observations (Phase 1 limitation)",
            f"inflation = {inflation}",
            f"ensemble size N = {N}",
        ],
        diagnostics={
            "n_observations": len(trace.observations),
            "n_updates": n_updates,
            "trace_id": trace.trace_id,
            "final_mean": x_ens.mean(axis=0).tolist(),
            "final_std": x_ens.std(axis=0, ddof=0).tolist(),
        },
    )

    return ReplayResult(
        times=times,
        means=means,
        stds=stds,
        ensembles=ensembles,
        final_ensemble=x_ens,
        model_manifest=manifest,
    )
