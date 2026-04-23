# F4-b — CVaR Acquisition: /architect Protocol

**Prepared by:** /architect (via dev-orchestrator)
**Date:** 2026-04-17
**Status:** Protocol + implementation shipped in same turn.
**Scope basis:** `docs/01_scientific_advisor_report.md §A.3` §4 (robust BO).

---

## 1. Purpose

Add a Conditional-Value-at-Risk (CVaR) aggregation to the robust BO
layer in `OptimizationEngine`, complementing the existing mean-variance
weight (F4-a, `--robust-variance-weight λ`). CVaR penalises candidates
whose *worst* α-quantile of resamples is poor — useful when the
decision-maker cares about tail-risk more than average performance.

Trivial variant of F4-a by scope; the resample loop is reused
unchanged. Only the final aggregator differs.

## 2. Definitions

For a minimisation objective with resample values ``v₁, …, v_n`` sorted
ascending,

```
CVaR_α = mean of the top ⌈α·n⌉ values
       = average of the worst α-fraction
```

Example: with n = 10 samples and α = 0.3, CVaR₀.₃ = mean of the 3
largest (worst) values. α = 1.0 recovers ``mean``; α → 0 approaches
``max``.

Per-objective reduction: applied independently to each dimension of the
objective vector, since each dimension is scalar-minimised.

## 3. Deliverables

| File | Change | LOC |
|---|---|---|
| `src/emulsim/optimization/engine.py` | `robust_cvar_alpha` ctor param; CVaR branch in `_evaluate` | ~45 |
| `src/emulsim/__main__.py` | `--robust-cvar-alpha` CLI flag + routing | ~15 |
| `tests/test_f4b_cvar.py` | 6 tests | ~160 |

Total: ~220 LOC.

## 4. Algorithm

```
def cvar_per_objective(samples: np.ndarray, alpha: float) -> np.ndarray:
    # samples shape: (n_resamples, n_obj)
    n = samples.shape[0]
    k = max(1, int(np.ceil(alpha * n)))       # tail count
    sorted_asc = np.sort(samples, axis=0)      # small → large
    # For a minimisation problem, "worst" = largest, so take top-k.
    return sorted_asc[-k:].mean(axis=0)
```

## 5. Interface

```python
OptimizationEngine(
    ...,
    robust_variance_weight: float = 0.0,   # F4-a
    robust_n_samples: int = 0,
    robust_cvar_alpha: float = 0.0,        # F4-b (NEW)
)
```

Precedence: `robust_cvar_alpha > 0` takes precedence over
`robust_variance_weight`. Both 0 (default) → deterministic BO.

CLI:
```
python -m emulsim design --robust-cvar-alpha 0.3 --robust-n-samples 10 [...]
```

## 6. Tests

1. CVaR helper: monotone-increasing samples, α=0.3, n=10 → mean of top
   3 values.
2. CVaR(α = 1.0) recovers the sample mean.
3. CVaR decreases monotonically with α (at fixed samples).
4. Engine ctor rejects `robust_cvar_alpha` outside (0, 1].
5. Engine ctor rejects `robust_cvar_alpha > 0` with `robust_n_samples
   < 2`.
6. CVaR path takes precedence over mean-variance when both are set.

## 7. Gate G1

All 12 criteria met — CVaR is a pure-aggregator extension of the
existing F4-a infrastructure. No new failure modes.
