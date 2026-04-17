# v9.0 Module Registry

Living index of approved modules in the EmulSim v9.0 "Family-First" UI redesign.
Updated after every Phase 5 approval.

| M# | Name | Version | Model | LOC | Fix rounds | Files | Status | Approved at |
|----|------|---------|-------|-----|-----------|-------|--------|-------------|
| M0 | pre-flight | 1.0.0 | opus (orchestrator) | 0 | 0 | handover/v9/M0_preflight.md, widget_keys.lock | APPROVED | 2026-04-18 |
| M1 | scaffolding | 1.0.0 | opus (arch+coder) | ~350 | 0 | tabs/m1/__init__.py, tabs/m1/family_selector.py, tabs/m1/hardware_section.py, tabs/m1/formulation_agarose_chitosan.py, tabs/m1/formulation_alginate.py, tabs/m1/formulation_cellulose.py, tabs/m1/formulation_plga.py, tabs/m1/crosslinking_section.py, tabs/m1/targets_section.py, tabs/m1/material_constants.py, tabs/m1/runner.py, visualization/nav.py, visualization/ui_links.py | APPROVED | 2026-04-18 |
| M2 | hardware_mode_radio | 1.0.0 | opus | ~60 | 0 | tabs/m1/hardware_section.py (populated), app.py (sidebar stripped), tabs/tab_m1.py (local render) | APPROVED | 2026-04-18 |
| M3 | family_selector | 1.0.0 | opus | ~90 | 0 | tabs/m1/family_selector.py (populated), tabs/tab_m1.py (wires selector + non-A+C placeholder) | APPROVED | 2026-04-18 |
| M4 | ac_extraction | 1.0.0 | opus | ~560 (moved) | 0 | tabs/m1/formulation_agarose_chitosan.py, crosslinking_section.py, targets_section.py, material_constants.py (all populated); tabs/tab_m1.py (calls the 4 new modules, params-build unchanged) | APPROVED | 2026-04-18 |
| M5 | alginate | 1.0.0 | opus | ~130 | 0 | tabs/m1/formulation_alginate.py (populated); tabs/tab_m1.py::_render_non_ac_family | APPROVED | 2026-04-18 |
| M6 | cellulose | 1.0.0 | opus | ~130 | 0 | tabs/m1/formulation_cellulose.py (populated) | APPROVED | 2026-04-18 |
| M7 | plga | 1.0.0 | opus | ~115 | 0 | tabs/m1/formulation_plga.py (populated) | APPROVED | 2026-04-18 |
| M8 | nav_hide | 1.0.0 | opus | ~30 | 0 | app.py (CSS inject to hide auto-page list); ui_links.py (URL builder); formulation_*.py (inline links) | APPROVED | 2026-04-18 |
| M9 | toml_roundtrip | 1.0.0 | opus | ~60 | 0 | tests/test_toml_roundtrip_family.py | APPROVED | 2026-04-18 |
| M10 | acceptance | 1.0.0 | opus | 0 | 0 | handover/v9/FINAL_HANDOVER.md | APPROVED |

## Integration status

| Interface | From | To | Status |
|---|---|---|---|
| FamilyContext | family_selector | tab_m1 (direct, M3) → runner (M4+) | LIVE for tab_m1 (M3); PENDING runner integration (M4) |
| HardwareContext | hardware_section | tab_m1 (direct, M2) → runner (M4+) | LIVE for tab_m1 (M2, radio only); PENDING full context (M4) |
| AgaroseChitosanContext | formulation_agarose_chitosan | runner | PENDING (M4) |
| AlginateContext | formulation_alginate | runner | PENDING (M5) |
| CelluloseContext | formulation_cellulose | runner | PENDING (M6) |
| PLGAContext | formulation_plga | runner | PENDING (M7) |
| CrosslinkingContext | crosslinking_section | runner | PENDING (M4) |
| TargetsContext | targets_section | runner | PENDING (M4) |
| material_overrides (dict) | material_constants | runner | PENDING (M4) |
| reagent_detail URL | ui_links.build_reagent_link | all formulation_* modules | PENDING (M8) |
| st.navigation | nav.build_navigation | app.py | PENDING (M8) |

## Gate status

- G1 (Protocol): architect plan + SA findings on record (prior two assistant turns in this conversation)
- G2 (Implementation): M1 — 13/13 modules import cleanly, smoke 3/3 pass
- G3 (Audit): M1 — no HIGH findings; stubs raise NotImplementedError as designed

## Next milestone

**M2 — Hardware Mode relocation** (Sonnet, 3 fix rounds, closes user bug #2).
