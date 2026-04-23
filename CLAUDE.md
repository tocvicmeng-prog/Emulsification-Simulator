# EmulSim — Claude Code instructions

## Design System
Always read `DESIGN.md` at the repo root before making any visual or UI
decisions. All font choices, colors, spacing, and aesthetic direction are
defined there. Do not deviate without explicit user approval.

In QA mode, flag any code that doesn't match DESIGN.md (e.g., hardcoded
colors outside the palette, fonts not in the approved list, purple
gradients, entrance animations, etc.).

## v9.0 Family-First UI contract

- Polymer family is the first M1 input (`tabs/m1/family_selector.py`). It
  drives all downstream rendering.
- Hardware Mode lives inside M1 Emulsification, NOT in Global Settings.
- Per SA §B matrix: never render crosslinker / cooling-rate / pore-model
  widgets for alginate / cellulose / PLGA. Never render chitosan+agarose %
  fields for any family other than AGAROSE_CHITOSAN.
- Always compare PolymerFamily members by `.value`, never by `is`. The
  Streamlit app reloads `emulsim.datatypes` on every rerun, minting a new
  enum class; identity comparisons silently break after the first rerun.
  See `RunReport.compute_min_tier` for the reference pattern.

## Known project quirks

- Repo path contains `文档` (Chinese). Any `print()` of an absolute path
  on Windows will crash with `UnicodeEncodeError` under cp1252. Scripts
  that print paths MUST do `sys.stdout.reconfigure(encoding="utf-8")` first
  (see `docs/user_manual/build_pdf.py` for the canonical pattern).

## v9.1.x decisions (read before changing the runtime stack)

- **Python pinned to `>=3.11,<3.13`** (ADR-001). Newer Python triggers
  scipy BDF + numba JIT cache issues. Don't widen the upper bound
  without re-running the test suite on the target Python.
- **Optimization stack pinned**: `botorch~=0.17.2`, `gpytorch~=1.15.2`,
  `torch~=2.11.0` (ADR-002). The smoke test
  `tests/test_optimization_smoke.py` exercises the duck-typed
  `FastNondominatedPartitioning` call site — it has a
  `# type: ignore[arg-type]` anchored to this pin range.
- **Solver method matrix** for scipy `solve_ivp`:
  - `module3_performance/catalysis/packed_bed.py` → LSODA (~700×
    faster than BDF on the non-stiff PFR + Michaelis-Menten problem)
  - `module3_performance/transport/lumped_rate.py::solve_lrm` → LSODA
    when no gradient adapter, BDF when both `gradient_program` and
    `equilibrium_adapter` are set (LSODA oscillates modes when the
    binding equilibrium is time-varying)
  - `module3_performance/orchestrator.py::run_gradient_elution` → BDF
    (always gradient path)
- **CI gates**: ruff must be 0; mypy must be 0 (cap lowered from 32
  to 0 in the v9.2.2 cleanup PR). PRs that add type errors fail CI;
  fix them in the same PR rather than raising the cap.
