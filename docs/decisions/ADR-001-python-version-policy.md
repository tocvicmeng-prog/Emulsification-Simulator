# ADR-001 — Python Version Policy

**Status:** Accepted
**Date:** 2026-04-19
**Decision driver:** v9.1 health audit; Py 3.14 incompatibility uncovered during M1 (test feedback loop) work.

## Context

EmulSim ships a one-click Windows installer that bundles a specific Python
runtime for end users. Until v9.0, `pyproject.toml` declared
`requires-python = ">=3.11"` with no upper bound, which permitted local dev
and CI to drift to whatever the latest CPython release happened to be.

During the v9.0 → v9.1 health audit (2026-04-19), running the test suite on
Python 3.14.3 produced two systemic failures:

1. **scipy `solve_ivp` BDF integration hangs** in many tests that exercise
   the PBE solver (level1 emulsification), the Cahn–Hilliard 1D/2D solvers
   (level2 gelation), and downstream pipeline tests that call into them.
   Symptoms: the process pegs a CPU at ~400 MB RSS with no progress for
   minutes. With `pytest-timeout` installed, these now fail loudly at the
   timeout instead of hanging the runner.
2. **`torch.jit.script` is explicitly unsupported on Python 3.14+** (per
   torch's own deprecation warning at site-packages/torch/jit/_script.py:1482).
   The optimization extra (`botorch`, `gpytorch`, `torch`) emits warnings
   and is at risk of silent miscompilation on future torch releases.

Numba support for Python 3.14 is also "best-effort"; the `.nbc` cache files
in `level1_emulsification/__pycache__/` show JIT compilation completed, but
cache-invalidation behaviour on 3.14 is unverified upstream.

The Windows installer ships Python 3.12 (matching the most recent stable
release at the time of v9.0). Customers who launch the installer therefore
run on a different Python than developers running on bleeding-edge CPython,
which is exactly the divergence that produced the v8.3.5 → v8.3.7 hotfix
cascade for installer-environment bugs.

## Decision

Pin the supported Python range to **`>=3.11,<3.13`** for v9.1 onward.

- **Lower bound (3.11):** matches the prior declared minimum and the typed
  syntax we depend on (PEP 654 exception groups, PEP 657 enhanced tracebacks).
  Dropping below 3.11 is not motivated.
- **Upper bound (`<3.13`):** excludes Python 3.13 and 3.14 from the supported
  matrix until upstream scientific deps (scipy BDF, numba, torch.jit, botorch)
  catch up. This is a conservative cap, not a permanent ceiling.

## Consequences

### Positive

- **Eliminates BDF hang storm.** The systemic test timeouts observed under
  3.14 are not present under 3.12 in upstream scipy issue trackers.
- **Aligns dev environment with installer environment.** The Python version
  developers test against is the same version end users run, removing a
  whole class of "works on my machine" install regressions.
- **Restores torch.jit support.** Removes the deprecation warning storm
  and the underlying risk of silent torch miscompilation.
- **Stable numba path.** Numba's officially-supported Python range covers
  3.11 and 3.12 cleanly.

### Negative

- **Developers on Python 3.14 must downgrade or use a per-project venv.**
  Recommended path:
  ```
  winget install Python.Python.3.12
  py -3.12 -m venv .venv
  .venv\Scripts\activate
  pip install -e ".[dev,ui]"
  ```
- **Future Python releases require an explicit decision.** Bumping to 3.13
  or 3.14 will require this ADR to be superseded; the bump must be gated on
  scipy / numba / torch upstream confirming support.

### Neutral

- The installer build pipeline is unchanged (already targets 3.12).
- Existing 3.11- and 3.12-compatible wheels for all transitive deps remain
  available on PyPI.

## Verification

Local verification on Python 3.12.10 in a fresh `.venv312` (2026-04-19):

- **Before (3.14.3):** Fast suite (`pytest -m "not slow"`) hangs at ~31%
  with no progress after 30+ minutes; ~150 tests reach the run phase.
- **After (3.12.10):** Fast suite progresses to ~49%; **283 tests pass**,
  1 fails, 2 xfail, 16 deselected, in 180s before the next hang.

So 3.12 is a large win but **not a silver bullet.** Two follow-ups remain:

1. **`test_solve_gelation_1d_fallback`** fails on 3.12 (was masked by the
   earlier hangs on 3.14). This is a real test or solver defect, not a
   version-pinning issue. Tracked separately.
2. **Packed-bed BDF tests** (`test_module3_catalysis::test_eta_in_range`
   and neighbours calling `solve_packed_bed`) still exceed the 60s
   timeout on 3.12. These need `@slow` markers or RHS conditioning fixes.
   Tracked separately.

This ADR is paired with the M2 CI scaffold (forthcoming), which will run
the test suite on `windows-latest` against the matrix `[3.11, 3.12]` on
every push and PR. A green CI signal is the operational verification of
this decision.

### Developer migration command

On Windows, from the project root:

```
winget install Python.Python.3.12 --silent
py -3.12 -m venv .venv312
.venv312\Scripts\activate
pip install -e ".[dev,ui]"
```

To include the optimization extra (torch / botorch / gpytorch), append
`,optimization` to the extras list. Note that `[optimization]` will be
pinned in a separate PR (module M4 in the v9.1 plan).

## References

- v9.1 architect blueprint, module M5
- Scientific Advisor root-cause memo (2026-04-19), Q4 (Python 3.14 risk
  assessment)
- pyproject.toml `[project] requires-python` change in the same commit
- torch deprecation warning: `site-packages/torch/jit/_script.py:1482`
