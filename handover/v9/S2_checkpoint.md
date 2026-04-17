# S2 Checkpoint — 2026-04-18

Session 2 closed after M2 + M3 approved. User bug #2 closed.

## Completed this session

- **M2 Hardware Mode relocation** — Removed `sim_mode` radio from `app.py` sidebar.
  Added `render_hardware_mode_radio()` to `tabs/m1/hardware_section.py`. Wired
  `tabs/tab_m1.py` to render the radio locally when `is_stirred=None` is passed
  from app.py. Widget key `m1v9_hardware_mode` (new, namespaced to avoid
  colliding with session state from v8.x users who had the sidebar selection).
  **Closes user bug #2.**

- **M3 Family selector** — `tabs/m1/family_selector.py` fully populated with
  `render_family_selector()`. 4-way radio (Agarose+Chitosan / Alginate /
  Cellulose (NIPS) / PLGA) at top of M1 tab, with scientific captions per
  family. When non-A+C is selected, tab_m1 short-circuits with an info banner
  telling the user which milestone will populate that family's UI. A+C path
  is unchanged (zero behavioural regression).

## Gates passed

- G1 Protocol: architect plan §2.1-2.2 — complete for M2, M3.
- G2 Implementation: all imports clean; `render_family_selector` and
  `render_hardware_mode_radio` both exposed at the module level.
- G3 Audit: no HIGH findings. Widget key namespace (`m1v9_*`) intentionally
  separate from `m1_*` to preserve session-state compat for existing users.
  Smoke tests 3/3 pass (same 2.77s baseline).

## User-visible changes at end of S2

- Sidebar simpler: Hardware Mode no longer present. Scientific Mode and v6.0
  Frameworks remain (genuinely global).
- M1 tab: a new "Polymer Family" section renders at the top, followed by
  "Emulsification (L1)" which now opens with the Hardware Mode radio.
- If user picks any of Alginate / Cellulose / PLGA, an info banner explains
  the family will be available in M5 / M6 / M7. Selecting A+C still lets the
  full existing pipeline run.

## Status of the 3 user bugs

| Bug | Description | Status after S2 |
|---|---|---|
| #1 | M1 only allows chitosan+agarose | **Partially addressed** — family selector is live and shows all 4 options, but selecting alginate/cellulose/PLGA shows a placeholder rather than a functional form. Closes fully at M7. |
| #2 | Hardware Mode belongs in M1, not Global Settings | **CLOSED** |
| #3 | App / reagent detail auto-list in top-left sidebar | Not yet addressed (closes at M8). |

## Scope for session 3 (S3)

Per the build schedule — **large milestone, ends with a handover**:
- M4 Agarose+Chitosan extraction: move ~600 LOC of A+C formulation,
  crosslinking, targets, and material_constants out of `tab_m1.py` into the
  dedicated modules (`formulation_agarose_chitosan.py`, `crosslinking_section.py`,
  `targets_section.py`, `material_constants.py`). Behaviour must stay bit-identical
  for the default TOML configs.

## Risk register updates

- Widget-key namespace collision: mitigated by using `m1v9_*` prefix for M2/M3
  new widgets. v8.x users' existing `m1_*` session state is preserved.
- Non-A+C placeholder: users who read the info banner will understand the gap;
  we decided against hiding the families from the selector because the backend
  already supports them and visibility builds trust with the scientific audit.

## Files touched in S2

- EDIT: `src/emulsim/visualization/app.py` (sidebar stripped, ~15 lines changed)
- EDIT: `src/emulsim/visualization/tabs/tab_m1.py` (local hardware radio + family selector wiring, ~30 lines changed)
- EDIT: `src/emulsim/visualization/tabs/m1/hardware_section.py` (stub → functional, ~45 new lines)
- EDIT: `src/emulsim/visualization/tabs/m1/family_selector.py` (stub → functional, ~60 new lines)

## Status
**DONE** — S2 scope complete, bug #2 closed, A+C backward-compatible.
