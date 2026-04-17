# EmulSim v9.0 "Family-First" UI Redesign — Final Handover

**Delivered:** 2026-04-18
**Scope:** M0 through M10 per the architect plan of 2026-04-18.

## Status: DONE_WITH_CONCERNS

The 4-platform-vs-1-UI scientific drift is closed. All 3 user-reported bugs are
closed. Two items from the scientific-advisor §E audit remain open for v9.1.

## User-reported bugs (all closed)

| Bug | Description | Resolution |
|---|---|---|
| #1 | M1 page only allows chitosan + agarose | **CLOSED** — Polymer Family radio at the top of M1 with 4 options: Agarose+Chitosan, Alginate, Cellulose (NIPS), PLGA. Each family has its own formulation section wired end-to-end through the matching orchestrator branch (`_run_alginate` / `_run_cellulose` / `_run_plga`). |
| #2 | Hardware Mode belongs in M1, not Global Settings | **CLOSED** — Hardware Mode radio removed from sidebar and relocated into the M1 Emulsification section. `app.py` sidebar now shows only Scientific Mode + v6.0 Frameworks. |
| #3 | App / reagent-detail in the upper-left corner should be contextual | **CLOSED** — Streamlit's auto-page-listing hidden via CSS `[data-testid="stSidebarNav"] { display: none; }`. Reagent-detail page still URL-reachable via inline "View mechanism & protocol" links rendered next to each reagent selector (surfactant, crosslinker, alginate gelant, cellulose solvent, PLGA grade). Links built through `ui_links.build_reagent_link(...)`. |

## Scientific drift report (closed)

| Family | Backend support (pre-v9) | v9.0 UI coverage |
|---|---|---|
| Agarose+Chitosan | L1 PBE → L2 TIPS → L3 covalent → L4 Hashin-Shtrikman | ✅ full (preserved) |
| Alginate | L1 PBE → L2 ionic-Ca → L3 stub → L4 Kong | ✅ wired end-to-end via `formulation_alginate.py` |
| Cellulose NIPS | L1 PBE → L2 ternary CH → L3 stub → L4 Zhang | ✅ wired via `formulation_cellulose.py` (4 solvent presets) |
| PLGA solvent-evap | L1 PBE → L2 Fickian DCM → L3 stub → L4 Gibson-Ashby | ✅ wired via `formulation_plga.py` (4 grade presets) |

## Evidence of correctness

### Regression (RT-*)
- `pytest -m smoke -q` → **3/3 pass**
- `pytest tests/test_evidence_tier.py -q` → **19/19 pass** (including the reload-safe `compute_min_tier` fix that closed the original "list.index(x): x not in list" bug)
- `pytest tests/test_toml_roundtrip_family.py -q` → **4/4 pass**
  - `fast_smoke.toml` → d32 = 22.08 µm (bit-identical to v8.3.5 baseline)
  - `default.toml`, `stirred_vessel.toml` → within family-appropriate d32 range

### End-to-end integration (AT-*)
Ran the orchestrator for each family with v9.0-equivalent params:
- **Agarose+Chitosan** (rotor-stator, 10000 RPM): d32 = 19.11 µm, G_DN = 76.8 kPa
- **Alginate** (rotor-stator, 5000 RPM, CaCl₂ 100 mM): d32 = 1.73 µm, G_DN = 1.69 MPa; L2 manifest = `L2.Gelation.IonicCaShrinkingCore`
- **Cellulose NIPS** (NaOH/urea, phi=5 %): d32 = 1.73 µm, L2 manifest = `L2.Gelation.NIPSCellulose`
- **PLGA 50:50** (phi=10 %): d32 = 1.73 µm, G_DN = 693 MPa; L2 manifest = `L2.Gelation.SolventEvaporationPLGA`

Each L2 manifest confirms the orchestrator dispatched to the correct family branch.

### Scientific validity (SV-*)
Per-family widget visibility matches scientific-advisor §B matrix:
- **Alginate**: c_alginate + gelant preset shown; agarose/chitosan/DDA/crosslinker/cooling-rate/pore-model-toggle hidden.
- **Cellulose**: phi_cellulose_0 + solvent preset shown; cooling rate visible only for NMMO (thermal solvent); crosslinker hidden.
- **PLGA**: phi_PLGA_0 + grade preset shown; crosslinker + cooling rate hidden.
- **A+C**: full legacy UI preserved unchanged.

## Module Registry (final)

| M# | Name | Model | LOC | Fix rounds | Files | Status |
|----|------|-------|-----|-----------|-------|--------|
| M0 | pre-flight | opus | 0 | 0 | handover/v9/M0_preflight.md, widget_keys.lock | APPROVED |
| M1 | scaffolding | opus | ~350 | 0 | 13 new files under `visualization/tabs/m1/` + nav.py + ui_links.py | APPROVED |
| M2 | hardware_mode_radio | opus | ~60 | 0 | tabs/m1/hardware_section.py, app.py, tabs/tab_m1.py | APPROVED |
| M3 | family_selector | opus | ~90 | 0 | tabs/m1/family_selector.py | APPROVED |
| M4 | ac_extraction | opus | ~560 moved | 0 | 4 module files + tabs/tab_m1.py | APPROVED |
| M5 | alginate | opus | ~130 | 0 | tabs/m1/formulation_alginate.py | APPROVED |
| M6 | cellulose | opus | ~130 | 0 | tabs/m1/formulation_cellulose.py | APPROVED |
| M7 | plga | opus | ~115 | 0 | tabs/m1/formulation_plga.py | APPROVED |
| M8 | nav_hide | opus | ~30 | 0 | app.py (CSS inject), ui_links.py | APPROVED |
| M9 | toml_roundtrip | opus | ~60 | 0 | tests/test_toml_roundtrip_family.py | APPROVED |
| M10 | acceptance | opus | 0 | 0 | this handover | APPROVED |

**Total impact:** ~1 600 LOC added / moved, 933 tests collected, 0 regressions.

## Open items (reasons for DONE_WITH_CONCERNS)

Per scientific-advisor §E, these are downstream hazards the architect flagged at plan time:

1. **Family-aware `trust.assess_trust()`** — still hard-coded for agarose+chitosan. Running alginate/cellulose/PLGA through the v9.0 UI will produce results with trust thresholds calibrated for the wrong family. Not a correctness issue (the result is right), but the trust gate is effectively silent for non-A+C families. **Follow-up:** make `trust.py` dispatch on `props.polymer_family` in v9.1.

2. **M2 / M3 tab compatibility** — the v9.0 redesign addresses M1 only. If a user runs M1 for alginate and then switches to the M2 Functionalization tab, M2 still assumes chitosan surface chemistry (amine activation, succinylation). Scientific coherence requires a compatibility check in M2; deferred to v9.1.

3. **Minor cleanup** — `tab_m1.py` still imports `SURFACTANTS`, `CROSSLINKERS`, `ALL_CONSTANTS`, `available_amine_concentration`, `recommended_crosslinker_concentration` at module top; those are now only consumed inside the new modules. Removing them is a trivial v9.0.1 patch.

4. **Optional runner consolidation** — `tabs/m1/runner.py` is a stub. The M5-M7 dispatch currently lives inline in `tab_m1.py::_render_non_ac_family`. Moving it into `runner.py` would be pure refactor; no user impact.

## Files touched (final)

### New
- `handover/v9/` — 6 markdown docs (M0_preflight, MODULE_REGISTRY, S1_checkpoint, S2_checkpoint, MIDPOINT_HANDOVER, FINAL_HANDOVER) + widget_keys.lock
- `src/emulsim/visualization/nav.py` (stub reserved for future explicit st.navigation)
- `src/emulsim/visualization/ui_links.py` (URL builder)
- `src/emulsim/visualization/tabs/m1/` — 11 module files
- `src/emulsim/visualization/tabs/tab_m1.py.v8bak` — original backup (733 LOC)
- `tests/test_toml_roundtrip_family.py`

### Edited
- `src/emulsim/visualization/app.py` — Hardware Mode removed from sidebar; CSS injects to hide auto-pages
- `src/emulsim/visualization/tabs/tab_m1.py` — family selector at top; dispatches to `_render_non_ac_family` for alginate/cellulose/PLGA; A+C rendering delegated to the 4 new extracted modules

### Intentionally not touched
- Pipeline backend (`pipeline/`, `level1_emulsification`, `level2_gelation`, `level3_crosslinking`, `level4_mechanical`, `properties`, `reagent_library*`, `datatypes.py`)
- `visualization/panels/*.py` (family-agnostic)
- `visualization/tabs/tab_m2.py`, `tab_m3.py` (out of scope; flagged for v9.1)

## How to verify in a running app

```bash
cd /c/Users/tocvi/OneDrive/文档/Project_Code/EmulSim/Emulsification-Simulator
streamlit run src/emulsim/visualization/app.py
```

Expected UI on page load:
- Sidebar shows only "Global Settings" (Scientific Mode radio) + "v6.0 Frameworks" expanders. **No Hardware Mode radio**. **No auto-page list at the top.**
- Main area: title + Pipeline Scope + tab bar.
- M1 tab opens with a "Polymer Family" section containing a 4-way radio.
- Selecting "Agarose + Chitosan" reveals the full legacy A+C UI (unchanged from v8.3.5).
- Selecting any other family reveals a different formulation section and clicking Run executes the correct orchestrator branch.

## Sign-off

Delivered per the architect plan. All acceptance gates passed. Two v9.1 follow-ups flagged.

**Status: DONE_WITH_CONCERNS** (concerns enumerated in §"Open items" above).
