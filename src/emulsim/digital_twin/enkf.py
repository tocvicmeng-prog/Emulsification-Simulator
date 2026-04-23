"""Node F2 Phase 1: stochastic Ensemble Kalman Filter update.

Phase 1 supports scalar observations only (the typical digital-twin
use case: d32 measured at one time, T measured at one time, …). Vector
observations with a matrix ``R`` come in Phase 2.

Algorithm (stochastic EnKF, Evensen 1994):

.. code-block:: text

    x_ens:         (N, n_state)  prior ensemble
    y_fc:          (N,)          forecast observable per member (= h(x_ens))
    x̄, H̄:         ensemble means
    X' = x_ens - x̄
    Y' = y_fc  - H̄
    P_xy = (X'.T @ Y') / (N - 1)      (n_state,)
    P_yy = (Y' @ Y')  / (N - 1)        (scalar)
    K    = P_xy / (P_yy + R)           (n_state,)
    y_pert^(i) = y_obs + ε^(i),  ε ~ N(0, R)
    x_post^(i) = x_ens^(i) + K · (y_pert^(i) - y_fc^(i))

Optional multiplicative prior inflation (Anderson & Anderson 1999)
inflates the deviation from the mean before the update to counter
sampling-driven ensemble collapse in long replays.
"""

from __future__ import annotations

from typing import Optional

import numpy as np


def enkf_update(
    x_ensemble: np.ndarray,
    y_forecast: np.ndarray,
    y_observed: float,
    R: float,
    rng: Optional[np.random.Generator] = None,
    inflation: float = 1.0,
) -> np.ndarray:
    """Stochastic EnKF scalar-observation update.

    Parameters
    ----------
    x_ensemble : np.ndarray, shape (N, n_state)
        Prior ensemble. ``N`` = ensemble size, ``n_state`` = state
        dimension.
    y_forecast : np.ndarray, shape (N,)
        Forecast observable per ensemble member (the output of
        ``h(x_ensemble[i])``).
    y_observed : float
        Scalar observation value.
    R : float
        Observation-noise variance (``noise_std**2``). Must be ≥ 0.
    rng : np.random.Generator, optional
        RNG for the stochastic perturbation. Defaults to a
        seed-independent default RNG.
    inflation : float
        Multiplicative prior inflation factor on the ensemble
        deviation from the mean. Applied **before** the Kalman update.
        Default 1.0 = no inflation.

    Returns
    -------
    np.ndarray, shape (N, n_state)
        Posterior ensemble.

    Raises
    ------
    ValueError
        If shapes are incompatible, N < 2, R < 0, or inflation < 0.
    """
    x_ens = np.asarray(x_ensemble, dtype=float)
    if x_ens.ndim != 2:
        raise ValueError(
            f"x_ensemble must be 2-D (N, n_state); got shape {x_ens.shape}"
        )
    N, n_state = x_ens.shape
    if N < 2:
        raise ValueError(f"EnKF requires N >= 2, got N = {N}")

    y_fc = np.asarray(y_forecast, dtype=float).ravel()
    if y_fc.size != N:
        raise ValueError(
            f"y_forecast size {y_fc.size} != ensemble size {N}"
        )

    if R < 0.0:
        raise ValueError(f"R must be >= 0, got {R}")
    if inflation < 0.0:
        raise ValueError(f"inflation must be >= 0, got {inflation}")

    if rng is None:
        rng = np.random.default_rng()

    # Apply multiplicative prior inflation
    if inflation != 1.0:
        x_mean = x_ens.mean(axis=0, keepdims=True)
        x_ens = x_mean + inflation * (x_ens - x_mean)
        # y_fc is a snapshot of h(x) — inflation on x doesn't
        # automatically re-forecast y. For a linear h the following
        # equivalent is correct:
        y_mean = y_fc.mean()
        y_fc = y_mean + inflation * (y_fc - y_mean)

    # Anomalies
    x_bar = x_ens.mean(axis=0)          # (n_state,)
    y_bar = float(y_fc.mean())           # scalar
    X_prime = x_ens - x_bar              # (N, n_state)
    Y_prime = y_fc - y_bar               # (N,)

    # Covariances (scalar-y simplification)
    P_xy = (X_prime.T @ Y_prime) / (N - 1)       # (n_state,)
    P_yy = float((Y_prime @ Y_prime) / (N - 1))   # scalar

    denom = P_yy + R
    if denom <= 0.0:
        # Degenerate: ensemble has collapsed AND R = 0. No update
        # possible. Return prior.
        raise ValueError(
            "EnKF update degenerate: P_yy + R <= 0 (ensemble collapse "
            "with zero observation noise)."
        )
    K = P_xy / denom                    # (n_state,)

    # Perturbed observations — np.asarray() normalises the rng.normal return
    # (which mypy types as float|ndarray) to a concrete ndarray.
    if R > 0.0:
        eps = np.asarray(rng.normal(loc=0.0, scale=np.sqrt(R), size=N))
    else:
        eps = np.zeros(N)
    y_pert = y_observed + eps            # (N,)

    # Posterior update
    innovation = y_pert - y_fc           # (N,)
    x_post = x_ens + np.outer(innovation, K)
    return x_post
