# ADR-002 — Optimization Stack Version Pin

**Status:** Accepted
**Date:** 2026-04-19
**Decision driver:** v9.1 health audit M4; mypy stub error in
`src/emulsim/optimization/engine.py:410` and broader risk of API drift across
botorch / gpytorch / torch.

## Context

The v9.0 `pyproject.toml` declared the `[optimization]` extra with loose
lower bounds:

```
optimization = ["botorch>=0.11", "gpytorch>=1.11", "torch>=2.1"]
```

The v9.1 health audit surfaced two specific drift symptoms:

1. **Mypy stub error** at `engine.py:410`: `qLogExpectedHypervolumeImprovement`
   declares `partitioning: NondominatedPartitioning`, but the runtime call site
   passes `FastNondominatedPartitioning`. Runtime works (verified by the M4
   smoke test); the stub is wrong because `FastNondominatedPartitioning`
   inherits from `FastPartitioning`, not from `NondominatedPartitioning`.
2. **Latitude in upgrade range** under the loose lower bounds: a fresh
   `pip install -e ".[optimization]"` could pull botorch 0.18+ at any time,
   at which point the stub-runtime divergence above could become a real
   API break (or vice versa, the stub could be fixed).

## Decision

Pin the optimization stack to compatible-release ranges allowing patch
updates only:

```
optimization = [
    "botorch~=0.17.2",
    "gpytorch~=1.15.2",
    "torch~=2.11.0",
]
```

`~=0.17.2` allows any 0.17.x but excludes 0.18+. Patch-only is the right
range because:

- BoTorch is pre-1.0 and minor bumps regularly rename or restructure
  acquisition modules (e.g., `logei` namespace).
- GpyTorch and Torch tend to keep API stable across patch releases but can
  shift compiled tensor semantics across minor releases.

## Consequences

### Positive

- **Reproducible installs.** Anyone running `pip install -e ".[optimization]"`
  gets the runtime stack the smoke test was verified against.
- **Mypy stub bug is contained.** The `# type: ignore[arg-type]` annotation
  at the partitioning= call site is anchored to a known-good botorch range.
- **CI matrix is meaningful.** Future M2 CI jobs that include `[optimization]`
  will exercise the same stack the developer ran.

### Negative

- **Security patches require manual bump.** When upstream botorch or torch
  release a security-driven minor (rare for these libraries), this ADR must
  be revisited and re-verified.
- **Downstream library compat.** If a future EmulSim feature depends on a
  newer botorch acquisition function, the pin must be lifted in the same
  PR with re-verification of the partitioning call.

### Neutral

- Numba is not in this extra and is unaffected.
- The `[dev]` and `[ui]` extras are not changed.

## Verification

The M4 verification harness (committed as
`tests/test_optimization_smoke.py`, marker `@pytest.mark.smoke`) exercises:

1. Construction of `FastNondominatedPartitioning(ref_point, Y)` with a tiny
   synthetic dataset.
2. `qLogExpectedHypervolumeImprovement(model, ref_point, partitioning, sampler)`
   instantiation with the FastNondominatedPartitioning passed in.
3. `optimize_acqf` end-to-end with `q=1, num_restarts=2, raw_samples=16`.

If any of these fails after a future bump, the pin is doing its job.

## References

- v9.1 architect blueprint, module M4
- `src/emulsim/optimization/engine.py:410` (the duck-typed call site)
- `tests/test_optimization_smoke.py` (the verification)
- ADR-001 (Python version policy) — companion pin for the dev environment
