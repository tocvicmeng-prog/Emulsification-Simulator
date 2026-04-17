# EmulSim 9.0.0 — Release Notes

**Release date:** 2026-04-18
**Platform:** Windows 11 x64 (Windows 10 x64 supported down to build 17763)
**Edition:** First Edition + Appendix J — Functionalisation Wet-Lab Protocols
**Licence:** GPL-3.0 — intellectual property: Holocyte Pty Ltd
**Upstream:** https://github.com/tocvicmeng-prog/Emulsification-Simulator

## What's new in v9.0.0 ("Family-First")

Major UI redesign. EmulSim was simulating four polymer platforms but
exposing only one to the user. v9.0 reshapes Module 1 around the
scientific truth.

### New in the UI

- **Polymer Family** radio at the top of Module 1. Four platforms, each
  with its own formulation section wired end-to-end:
  - **Agarose + Chitosan** — thermal TIPS + optional covalent crosslinking
  - **Alginate** — ionotropic Ca²⁺ gelation (CaCl₂ external bath OR
    GDL + CaCO₃ internal release with honoured k_release)
  - **Cellulose NIPS** — NaOH/urea, NMMO, EMIM-Ac, DMAc/LiCl presets
  - **PLGA** — 50:50, 75:25, 85:15, PLA grade presets
- **Hardware Mode** relocated from Global Settings into the M1 Emulsification
  section. It's an M1-local control (L1 PBE solver dispatch only).
- **Reagent-detail** auto-page nav hidden; page reachable only from inline
  "View mechanism & protocol" links next to each reagent selector.
- Scientific Instrument design system. Geist + Geist Mono + JetBrains Mono
  typography, slate neutrals, single teal accent. Dark-first.

### New in the documentation

- **Appendix J — Functionalisation Wet-Lab Protocols** (new appendix).
  44 reagent-level protocols with full SDS-lite safety blocks (GHS
  pictograms, H-codes, PPE, waste streams) written for users who have
  never done affinity chromatography. Covers hydroxyl activation, ligand
  coupling, protein coupling, spacer arms, IMAC metal charging, protein
  pretreatment, washing, quenching. Shipped as its own PDF download.

### Bug fixes (all back-ported from v8.3.x hotfix work)

- **compute_min_tier reload safety** — Streamlit's module reload produced
  a "list.index(x): x not in list" crash on every rerun. Fixed by
  comparing `ModelEvidenceTier` by `.value` instead of enum identity.
- **Family dispatch reload safety** — same enum-identity hazard caught in
  the v9 UI by `/codex review`; fixed at both dispatch points.
- **Alginate internal-release k_release** now actually reaches the solver
  via `effective_bath_concentration(profile, t_end)`.
- **build_pdf.py** Windows cp1252 crash when repo path contains non-ASCII
  characters; fixed with `sys.stdout.reconfigure(encoding="utf-8")`.

### Known limitations (deferred to v9.1)

- `trust.assess_trust()` thresholds are calibrated for Agarose + Chitosan
  only. Non-A+C families produce physically correct outputs but the
  trust gate is effectively silent for them.
- M2 / M3 tabs still assume chitosan surface chemistry. If you ran M1
  for alginate / cellulose / PLGA, do NOT rely on M2 / M3 output.

## What's in the installer

- EmulSim Python package (wheel, `emulsim-9.0.0-py3-none-any.whl`)
- User Manual First Edition (PDF)
- Appendix J — Functionalisation Protocols (PDF)
- 3 starter TOML configs (`default`, `fast_smoke`, `stirred_vessel`)
- Launch scripts (Web UI, CLI)
- Install / uninstall scripts

**Clean package.** No build artefacts, no `__pycache__`, no development
log files, no intermediate test output, no internal design / handover
documents. Runtime files only.

## Runtime prerequisites

- Windows 11 x64 (Windows 10 x64 ≥ 17763 also supported)
- Python 3.11 or newer on system PATH. The installer checks for Python
  and offers to open the download page if missing.

## One-click install

1. Download `EmulSim-9.0.0-Setup.exe` from the GitHub Releases page.
2. Double-click. Accept the EULA (GPL-3.0, Holocyte Pty Ltd).
3. Click Next through the wizard.
4. The installer runs `install.bat` automatically (creates a per-release
   Python venv, installs the wheel with UI + optimisation extras,
   smoke-tests the import).
5. Launch from the Start Menu (**EmulSim (Web UI)**) or the desktop icon.

## Changing from a prior version

- Uninstall the previous version first via Control Panel → Programs.
- Install v9.0.0.
- Session state between v8.3.x and v9.0 is preserved for the widget keys
  that still exist (`m1_rpm`, `m1_surfactant`, `m1_crosslinker`, etc.).

## Support

- Issues: https://github.com/tocvicmeng-prog/Emulsification-Simulator/issues
- Source (GPL-3.0): https://github.com/tocvicmeng-prog/Emulsification-Simulator
- Latest release: https://github.com/tocvicmeng-prog/Emulsification-Simulator/releases
