# M0 — Pre-flight Gate (EmulSim v9.0 UI redesign)

**Date:** 2026-04-18
**Verdict:** PASS

## Environment
- streamlit: **1.55.0** (≥1.36 → **Option A `st.navigation` available**; no fallback needed)
- python: 3.14.3
- os: Windows 11 Pro 10.0.26200

## Baseline
- pytest smoke: **3/3 pass** (2.86s)
- pytest total collected: 932
- Known warnings: `torch.jit.script` deprecated on Py3.14 (unrelated to this work)
- Pre-existing failures: none on smoke path

## Backup
- Source: `src/emulsim/visualization/tabs/tab_m1.py`
- Backup: `src/emulsim/visualization/tabs/tab_m1.py.v8bak` (733 lines)

## Widget-key lock
30 `m1_*` keys captured in `handover/v9/widget_keys.lock`. Migration rule:
- Any widget that still exists after the redesign MUST use the same key.
- New widgets use the prefix `m1v9_*` to avoid collision.

## Git state at pre-flight
- Dirty files (all expected from prior debug session): `build_pdf.py`, `polysaccharide_microsphere_simulator_first_edition.pdf`, `datatypes.py`, `visualization/app.py`, `visualization/tabs/tab_m1.py`, `tests/test_evidence_tier.py`
- Untracked: `tab_m1.py.v8bak` (intentional)

## Go/No-go
- Streamlit version: GO (Option A)
- Session-state migration: GO (widget-key lock captured)
- Baseline stability: GO (smoke green)
- Backup: GO (v8bak present)

**PROCEED to M1 Scaffolding.**
