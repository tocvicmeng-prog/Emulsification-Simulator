# EmulSim 9.1.1 — Release Notes

**Release date:** 2026-04-19
**Platform:** Windows 11 x64 (Windows 10 x64 supported down to build 17763)
**Edition:** First Edition + Appendix J — Functionalisation Wet-Lab Protocols
**Licence:** GPL-3.0 — intellectual property: Holocyte Pty Ltd
**Upstream:** https://github.com/tocvicmeng-prog/Emulsification-Simulator
**Latest source code:** download from the "Releases" tab on GitHub

## What's new in v9.1.1

A backlog burndown release. No new simulator features; the work is in
performance, code quality, and continuous integration so future
releases stay green and ship faster.

### Performance

- **Packed-bed reactor and chromatography LRM solvers** switched to the
  scipy LSODA integrator (auto-switches between stiff and non-stiff
  methods). The packed-bed test workload that took 86 s under the
  previous fixed-BDF solver now finishes in 0.12 s — about 700x faster
  with identical scientific output. Gradient elution stays on BDF
  because the time-varying binding equilibrium is genuinely stiff.
- **Cahn-Hilliard 2D smoke test** runs on a tiny 16x16 grid with fast
  cooling so it returns in ~2 s instead of the prior ~600 s of
  simulated time. The full 32x32 scientific test suite is unchanged.

### Code quality

- Mypy type-error count cut from 71 to 32. The CI gate now refuses
  any pull request that adds new type errors past that baseline; the
  cap lowers as future cleanup lands.
- 17 dead-code variable assignments removed across the scientific
  solver tree.

### CI / packaging

- The CI installer-smoke job now compiles the actual `.exe` via
  Inno Setup, silently installs it to a temporary directory, and
  verifies the install tree before the wheel is even pip-installed.
  Catches the class of `.bat` parser / CRLF / Access Denied
  regressions that drove the v8.3.5–v8.3.7 hotfix cascade.

## What's new in v9.1.0 ("Health audit hardening")

The v9.1 series is a health-driven hardening pass driven by the v9.0
audit. Highlights from the v9.1.0 release that landed v9.1.1's
foundation:

- **`pytest-timeout` everywhere.** Default 120 s timeout in pytest so
  any future hang fails loudly within two minutes instead of stalling
  the runner indefinitely.
- **Python pinned to 3.11–3.12** via `requires-python = ">=3.11,<3.13"`.
  Documented in ADR-001. Earlier versions had `>=3.11` open-ended,
  which let dev environments drift to Python 3.14 where scipy and
  torch had compatibility issues.
- **Optimization stack pinned**: `botorch ~= 0.17.2`, `gpytorch ~= 1.15.2`,
  `torch ~= 2.11.0`. Documented in ADR-002. A paired smoke test
  exercises the BoTorch acquisition function so future drift breaks
  the build.
- **Ruff lint count reduced from 173 to 0.** Per-file-ignores
  introduced for the documented Streamlit reload pattern in
  `visualization/app.py`.
- **GitHub Actions CI** with three jobs: quick (3.11 + 3.12 matrix
  with ruff + mypy + fast pytest), smoke (minimal install + smoke
  marker), installer-smoke (the full .exe + silent install gate
  described above).

## What was in v9.0.0 ("Family-First UI redesign")

For background, the v9.0.0 baseline introduced:

- **Polymer Family** radio at the top of Module 1 — four platforms
  (Agarose+Chitosan, Alginate, Cellulose-NIPS, PLGA) each wired
  end-to-end with their own formulation section.
- **Hardware Mode** lives inside M1 Emulsification (Stirred Vessel
  vs Rotor-Stator legacy) — moved out of Global Settings.
- **Appendix J** — 44 wet-lab functionalisation protocols with
  full SDS-lite safety blocks, shipped as a second PDF alongside
  the First Edition manual.

## Installation

See `INSTALL.md` for the full Windows installation guide.
Briefly: install Python 3.11 or 3.12 from python.org (tick "Add
python.exe to PATH"), then run `install.bat` from the install
directory or let the installer's post-install step do it for you.

## Source code

The canonical upstream repository is
<https://github.com/tocvicmeng-prog/Emulsification-Simulator>.
Each tagged release is published under the "Releases" tab.
This release corresponds to tag `v9.1.1`.

## Licence

GPL-3.0. The full licence text ships as `LICENSE.txt` in the
install directory and is also at
<https://www.gnu.org/licenses/gpl-3.0.en.html>. The intellectual
property in this software belongs to Holocyte Pty Ltd; see
`LICENSE_AND_IP.txt` for the user agreement shown by the
installer.
