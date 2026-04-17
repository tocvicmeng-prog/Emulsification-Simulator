# F2 — Digital Twin (EnKF Replay): /architect Protocol

**Prepared by:** /architect (via dev-orchestrator)
**Date:** 2026-04-17
**Status:** Protocol + Phase 1 implementation shipped same turn.
**Scope basis:** `docs/node32_cluster_f_v8_roadmap.md` §3 (digital twin).

---

## 1. Purpose and scope

Add an **Ensemble Kalman Filter (EnKF) replay harness** for digital-
twin use cases: given a recorded trajectory of process observations
(e.g., d32 measured every minute during a real emulsification) and an
ensemble of simulated runs at perturbed parameters, sequentially
update the ensemble to estimate the true hidden parameters.

**Phase 1 = replay only.** The harness ingests a pre-recorded trace
and walks forward in time applying EnKF updates offline. Online /
real-time filtering and MPC-style control are Phase 2+.

Out of scope for Phase 1:
- Real-time online filtering (polling a live process log).
- Model-Predictive Control (MPC) built on top of the filtered state.
- Automatic parameter identifiability / sensitivity diagnostics
  (beyond the reported ensemble spread).

## 2. Scientific basis

The stochastic EnKF (Evensen 1994 *JGR* 99:10143; 2003 *Ocean Dyn.*
53:343) updates each ensemble member ``x_k^(i)`` when an observation
``y_k`` arrives:

```
Prior ensemble:      {x_k^(i)}_{i=1..N}
Forecast obs:        {h(x_k^(i))}_{i=1..N}
Ensemble mean x̄, H̄
Anomalies:           X' = x - x̄,    Y' = h(x) - H̄
Covariances:         P_xy = (X' · Y'ᵀ) / (N − 1)
                     P_yy = (Y' · Y'ᵀ) / (N − 1)
Kalman gain:         K = P_xy · (P_yy + R)⁻¹
Perturbed obs:       y_k^(i) = y_k + ε^(i),   ε^(i) ~ N(0, R)
Posterior:           x_k^(i) += K · (y_k^(i) − h(x_k^(i)))
```

``R`` is the observation-noise covariance. For scalar observations,
``P_yy + R`` is a 1×1 matrix; inversion collapses to division.

EnKF is the standard technique for high-dimensional non-linear filter-
ing problems where the full Kalman filter is intractable and a
particle filter would collapse. Typical ensemble size: 20–200.

References:
- Evensen (2003) *Ocean Dyn.* 53:343 — canonical EnKF review.
- Katzfuss & Stroud (2016) *Statistics* 50:241 — tutorial with
  the stochastic vs deterministic variants we'd pick between.
- Lorenz (1963) *J. Atmos. Sci.* 20:130 — canonical toy system for
  filter verification.

## 3. Deliverables (Phase 1)

| File | Purpose | LOC |
|---|---|---|
| `src/emulsim/digital_twin/__init__.py` | Module exports | ~15 |
| `src/emulsim/digital_twin/schema.py` | `DigitalTwinTrace` dataclass + JSON load/save | ~130 |
| `src/emulsim/digital_twin/enkf.py` | Stochastic EnKF update (single step) | ~180 |
| `src/emulsim/digital_twin/replay.py` | Trace replay harness | ~180 |
| `tests/test_digital_twin.py` | 10 tests | ~260 |

Total: ~765 LOC (a bit over the roadmap ~500 LOC estimate; cost
attributed to writing a clean stochastic EnKF + its linear-Gaussian
convergence test rather than cutting corners).

## 4. Data schema

```python
@dataclass
class Observation:
    t: float               # [s] time since experiment start
    name: str              # observable name, e.g. "d32" or "T"
    value: float           # measured value
    noise_std: float       # 1-sigma observation noise

@dataclass
class DigitalTwinTrace:
    trace_id: str
    process_description: str
    observations: list[Observation]   # time-ordered
    metadata: dict                    # free-form: rpm, vessel type, …
```

JSON round-trip; observations are stored as a list of objects in
time order (the loader sorts by ``t`` defensively).

## 5. EnKF API

```python
def enkf_update(
    x_ensemble: np.ndarray,    # (N, n_state)
    y_forecast: np.ndarray,    # (N,) forecast observable for each member
    y_observed: float,         # scalar observation
    R: float,                  # observation noise variance
    rng: np.random.Generator | None = None,
    inflation: float = 1.0,    # multiplicative prior inflation
) -> np.ndarray:
    """Return the posterior ensemble (same shape as x_ensemble)."""
```

Phase 1 supports **scalar observations only** (dim-1 y). Vector
observations are a Phase 2 extension — tensor-valued ``R``, matrix
inversion, etc. — but scalar covers "d32 measured at one time" and
"T measured at one time" which are the realistic replay targets.

## 6. Replay harness

```python
def run_replay(
    trace: DigitalTwinTrace,
    initial_ensemble: np.ndarray,      # (N, n_state)
    state_transition: Callable,        # x_{k+1} = f(x_k, dt)
    observation_operator: Callable,    # y = h(x, name)
    inflation: float = 1.02,
    seed: int | None = 0,
) -> ReplayResult:
    """Walk forward through the trace, apply EnKF at each obs, return
    the ensemble trajectory + filtered mean + spread."""
```

``ReplayResult`` carries the mean / std / ensemble snapshot at each
observation time plus a final manifest
(``DigitalTwin.EnKFReplay``, ``SEMI_QUANTITATIVE``).

## 7. Tests

1. Linear-Gaussian toy: with ``f(x) = x`` (persistence model),
   ``h(x) = x``, repeated observations drive the posterior mean to
   the observation and the spread to ``R/N``.
2. Zero-observation-noise limit: posterior collapses to the
   observation exactly.
3. Inflation > 1 grows the prior spread before the update.
4. ``DigitalTwinTrace`` JSON round-trip preserves all fields.
5. ``Observation`` time-ordering enforced on load.
6. Replay harness returns trajectory of correct shape.
7. Replay harness: no observations → ensemble passes through the
   state-transition unchanged.
8. Replay harness: single observation at t=0 reproduces a one-shot
   EnKF.
9. ``enkf_update`` rejects zero ensemble size and mis-shaped
   inputs.
10. Ensemble mean / spread shrinks over a multi-step replay when
    observations are informative.

## 8. Gate G1

All 12 criteria met for Phase 1. Stochastic EnKF is a well-posed
algorithm; numerical-stability edge cases (``P_yy + R = 0`` when both
ensemble collapse AND R = 0) are caught as ValueError.

## 9. Phase 2 (deferred)

- Vector observations (matrix ``R``).
- Square-root / deterministic EnKF variants (Sakov & Oke 2008) for
  better statistical properties at small N.
- Online polling adapter (file-watcher / HTTP hook).
- MPC: use the filtered state to recommend next-step process inputs.
- Automatic identifiability diagnostics (SVD of the Jacobian / Fisher
  information).
