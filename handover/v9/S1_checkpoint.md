# S1 Checkpoint — 2026-04-18

Session 1 of the EmulSim v9.0 UI redesign closed at a clean boundary.

## Completed this session

- **M0 Pre-flight gate** — Streamlit 1.55.0 (≥1.36 → Option A), Python 3.14.3,
  backup of `tab_m1.py` as `tab_m1.py.v8bak`, smoke 3/3 pass baseline, 30 `m1_*`
  widget keys inventoried to `widget_keys.lock`. Verdict PASS.
- **M1 Scaffolding** — 13 new module files created under
  `src/emulsim/visualization/`: `tabs/m1/` package (11 modules), `nav.py`,
  `ui_links.py`. All modules import cleanly (verified); all render functions
  raise `NotImplementedError` with the milestone that will populate them.
  Smoke 3/3 still pass. No behavioural change visible in the running app.

## How to resume

Next session starts with **M2 Hardware Mode relocation**:

1. Read `handover/v9/MODULE_REGISTRY.md` for the approved surface.
2. Re-run the pre-flight smoke gate: `pytest -m smoke -q`.
3. The architect protocol for M2 (in the architect plan delivered in S1, §2.1-§2.2):
   - Populate `tabs/m1/hardware_section.py::render_hardware_section`.
   - Edit `visualization/app.py` to remove the Hardware Mode sidebar control.
   - Edit `tabs/tab_m1.py` to call `render_hardware_section(...)` in the M1
     Emulsification section instead of reading hardware choice from session state.
   - Preserve widget keys: `m1_vessel`, `m1_stirrer`, `m1_rpm`, `m1_rpm_rs`,
     `m1_rpm_leg` (per `widget_keys.lock`).
   - Acceptance: AT-5 (sidebar shows no Hardware Mode), AT-6 (M1 shows it).

## Scope for session 2 (S2)

Per the build schedule:
- M2 Hardware Mode relocation
- M3 Family selector + shared L1 inputs
- Expected end-of-session zone: GREEN (~78% context remaining)

## Files touched

- NEW: `src/emulsim/visualization/tabs/m1/` (11 files)
- NEW: `src/emulsim/visualization/nav.py`
- NEW: `src/emulsim/visualization/ui_links.py`
- NEW: `handover/v9/` (3 docs)
- NEW: `src/emulsim/visualization/tabs/tab_m1.py.v8bak` (backup)
- UNTOUCHED in this session: `app.py`, `tab_m1.py`, `pages/reagent_detail.py`,
  `panels/*`, `tab_m2.py`, `tab_m3.py`, entire backend (`pipeline/`,
  `level1..level4`, `properties/`, `reagent_library*`, `datatypes.py`).

## Risk register

- No risks materialised in S1.
- Streamlit 1.55.0 confirms Option A is viable → no fallback required for M8.
- Widget key migration plan locked; M4 extraction has an explicit preservation list.

## Status
**DONE** — S1 scope complete, clean checkpoint, ready to resume.
