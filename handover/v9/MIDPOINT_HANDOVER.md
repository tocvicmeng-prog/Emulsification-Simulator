# Midpoint Handover — EmulSim v9.0 "Family-First" UI redesign

**Date:** 2026-04-18
**Author:** Claude Code dev-orchestrator (Opus 4.7)
**Scope:** Half-way checkpoint after M0-M4 approved. Session-loss insurance.

## Status at handover

| Milestone | Status | User impact |
|---|---|---|
| M0 pre-flight | APPROVED | — |
| M1 scaffolding | APPROVED | — (invisible) |
| M2 Hardware Mode relocation | APPROVED | **User bug #2 CLOSED** |
| M3 Family selector | APPROVED | Family radio visible; 4 options shown, 3 are placeholders until M5-M7 |
| **M4 A+C extraction** | **APPROVED** | — (invisible; A+C pipeline behaviour preserved bit-for-bit) |
| M5 Alginate | pending | Closes bug #1 for alginate |
| M6 Cellulose | pending | Closes bug #1 for cellulose |
| M7 PLGA | pending | Closes bug #1 for PLGA |
| M8 Nav hide | pending | **Closes bug #3** |
| M9 TOML compat | pending | Regression protection |
| M10 Acceptance | pending | Final sign-off |

## What the next session inherits

### Files populated (functional, not stubs)

- `src/emulsim/visualization/nav.py` — STUB still (populates at M8)
- `src/emulsim/visualization/ui_links.py` — STUB still (populates at M8)
- `src/emulsim/visualization/tabs/m1/family_selector.py` — **populated (M3)**
- `src/emulsim/visualization/tabs/m1/hardware_section.py` — **populated (M2)** — only `render_hardware_mode_radio()` is live; the full `render_hardware_section()` is deferred
- `src/emulsim/visualization/tabs/m1/formulation_agarose_chitosan.py` — **populated (M4)**
- `src/emulsim/visualization/tabs/m1/formulation_alginate.py` — STUB (populates at M5)
- `src/emulsim/visualization/tabs/m1/formulation_cellulose.py` — STUB (populates at M6)
- `src/emulsim/visualization/tabs/m1/formulation_plga.py` — STUB (populates at M7)
- `src/emulsim/visualization/tabs/m1/crosslinking_section.py` — **populated (M4)**
- `src/emulsim/visualization/tabs/m1/targets_section.py` — **populated (M4)**
- `src/emulsim/visualization/tabs/m1/material_constants.py` — **populated (M4)**
- `src/emulsim/visualization/tabs/m1/runner.py` — STUB still (populates progressively in M5-M7)

### Files edited in live tree

- `src/emulsim/visualization/app.py` — Hardware Mode removed from sidebar; `is_stirred=None` passed to tab_m1
- `src/emulsim/visualization/tabs/tab_m1.py` — family selector at top; non-A+C shows info placeholder; M4 rendering now delegated to the four new modules

### Files intentionally unchanged

- `src/emulsim/pipeline/orchestrator.py` — backend complete, no change needed
- `src/emulsim/datatypes.py` — fields complete for all 4 families
- `src/emulsim/properties/*.py`, `reagent_library*.py` — authoritative backend
- `src/emulsim/visualization/panels/*.py` — family-agnostic, no change
- `src/emulsim/visualization/tabs/tab_m2.py`, `tab_m3.py` — out of scope (flagged in SA §E for v9.1)

## Verification at midpoint

- `pytest -m smoke -q` → **3/3 pass** (2.86s, same baseline)
- `python -c "from emulsim.pipeline.orchestrator import PipelineOrchestrator; from emulsim.config import load_config; from emulsim.properties.database import PropertyDatabase; orch=PipelineOrchestrator(db=PropertyDatabase()); r=orch.run_single(load_config('configs/fast_smoke.toml')); print(r.emulsification.d32*1e6)"` → **22.08 µm** (identical to v8.3.5 baseline)
- All 13 new modules import cleanly

## Session-state migration

Widget keys preserved (v8.x → v9.x compat): all 30 `m1_*` keys locked in `handover/v9/widget_keys.lock` are still honoured where the widget still exists.

New keys introduced by M2-M3: `m1v9_hardware_mode`, `m1v9_polymer_family` — namespaced to avoid collision with user-saved v8.x session state.

## How to resume (next session playbook)

1. Read: `handover/v9/MIDPOINT_HANDOVER.md` (this file), `MODULE_REGISTRY.md`, `widget_keys.lock`.
2. Confirm baseline: `pytest -m smoke -q` → 3/3 pass.
3. Start **M5 Alginate formulation**:
   - Populate `formulation_alginate.py::render_formulation_alginate`.
   - Add `FormulationParameters.polymer_family = PolymerFamily.ALGINATE` wiring in a dispatcher inside `tab_m1.py` (replace the info-banner short-circuit when family=ALGINATE with a real render call + run-button path).
   - The minimum runner needs to: (a) use the existing surfactant, T_oil, hardware context; (b) set `c_alginate`, `c_Ca_bath` (or internal-release params); (c) call `orch.run_single(params)` which already dispatches to `_run_alginate`.
   - Reference the gelant library at `src/emulsim/reagent_library_alginate.py::GELANTS_ALGINATE`.
   - Acceptance: AT-1, AT-2, SV-1 (see architect plan §7).
4. M6, M7, M8, M9 follow in order (see architect plan §4).
5. M10 runs the full acceptance suite + writes final handover.

## Risk register

- ✅ Streamlit version (1.55.0) — Option A viable for M8.
- ✅ Session-state collisions — namespaced `m1v9_*` keeps v8.x users safe.
- ⚠ SA §E out-of-scope items: family-aware trust (`trust.py`), M2/M3 tab compatibility. Flag at M10 as DONE_WITH_CONCERNS if not addressed.
- ⚠ `runner.py` stub — not yet used by tab_m1. M5 should decide whether to populate runner.py or keep the dispatch in tab_m1. Architect recommends moving to runner.py at M5 to set the pattern for M6/M7.
- ⚠ A small code duplication: `tab_m1.py` still imports `SURFACTANTS`, `CROSSLINKERS`, `ALL_CONSTANTS`, `available_amine_concentration`, `recommended_crosslinker_concentration` — now unused at top level. Cleanup candidate for M10.

## Status
**DONE** (midpoint checkpoint). User bug #2 shipped. Scientific drift reduced from 4-platforms-vs-1-UI to 4-platforms-vs-1-UI-plus-honest-placeholders. Zero regression on the A+C path. S3 (M0-M4) complete.
