# EmulSim 9.2.0 — Release Notes

**Release date:** 2026-04-24
**Platform:** Windows 11 x64 (Windows 10 x64 supported down to build 17763)
**Edition:** First Edition + Appendix J — Functionalisation Wet-Lab Protocols
**Licence:** GPL-3.0 — intellectual property: Holocyte Pty Ltd
**Upstream:** https://github.com/tocvicmeng-prog/Emulsification-Simulator
**Latest source code:** download from the "Releases" tab on GitHub

## What's new in v9.2.0

Adds **hyperlinked derivation pages** for every optimization suggestion the
M1 tab produces after a run. Clicking the new [📊 derivation] link beside
a suggestion opens a dedicated page that reconstructs the full reasoning
chain: formulas, intermediate numbers, assumptions, confidence tier, and a
nominal-plus-band numeric target the user can dial toward.

### Why

Before v9.2.0 the "Optimization Assessment" produced a flat list like
*"Increase crosslinker concentration — G_DN 3 kPa below target 10 kPa"*
without explaining *why* or *by how much*. Users had to guess at the
magnitude of the change, cross-reference the manual, and hope their
chosen delta was in the right ballpark. The new derivation pages turn
each suggestion into a complete, auditable scientific argument.

### What the derivation page contains

Every page is built from three canonical sections:

1. **Derivation logic and formula pathway** — a step-by-step reasoning
   chain with LaTeX-rendered equations, the physical constants used,
   and the intermediate numbers computed from your inputs.
2. **Target value or range** — the nominal best-estimate plus a lower
   and upper bound that bracket the acceptable outcome. Units and
   feasibility-flag included.
3. **Assumptions and confidence** — an honest statement of what the
   derivation takes for granted, plus the confidence tier (VALIDATED,
   SEMI_QUANTITATIVE, or QUALITATIVE_TREND).

When the underlying model is only `QUALITATIVE_TREND` (e.g. the
empirical L2 pore correlation), the page **refuses to show a numeric
target** and explains why — matching the evidence tier rather than
giving the illusion of quantitative accuracy.

### Covered suggestions (v9.2.0, all five at once)

| Key | Physics used |
|---|---|
| `increase_rpm` / `decrease_rpm` | Sprow (1967) Weber-number correlation inverted for d32; Kolmogorov-scale and Reynolds feasibility checks. |
| `adjust_cooling_rate` | Lumped-capacitance heat transfer + Cahn-Hilliard spinodal-dwell scaling inverted for target pore size. |
| `increase_crosslinker` | Rubber elasticity G = ν_e·RT inverted for crosslink fraction, then stoichiometric relation for [crosslinker]. |
| `reduce_polymer` | Semi-dilute G ~ c² inverted for a uniform scaling factor α on both polymer concentrations. |

### Where it appears

- M1 tab → Optimization Assessment section: each suggestion now renders
  as a bullet ending in `[📊 derivation]`.
- Click the icon → `/suggestion_detail?key=...&ctx...` opens in the same
  window; a **← Back to Simulator** link returns you to the run.
- URLs are self-contained (full parameter round-trip) so you can
  bookmark a specific derivation or share it with a colleague.

### Under the hood

- New `src/emulsim/suggestions/` package (8 modules, ~750 LOC).
- New `src/emulsim/properties/*_derivation.py` physics helpers
  (3 modules, ~450 LOC).
- New Streamlit page `pages/suggestion_detail.py`.
- 55 new tests across 4 test files; all physics helpers include a
  forward-inverse round-trip property test.
- Zero solver-code changes; zero new dependencies.
- ruff 0 findings; mypy 32 errors (at the regression cap, zero added).

## What was in v9.1.2 ("STMP crosslinker")

Added Sodium Trimetaphosphate (STMP, Na₃P₃O₉, CAS 7785-84-4) as a new
crosslinker, with a first-principles scientific audit, a triggerable
two-phase wet-lab protocol, and homogeneity guardrails built into the
UI.

### Why STMP matters

STMP fills a genuine gap in the crosslinker library: a **food-grade,
covalent, dual-network, triggerable** crosslinker. Existing options
trade off against each other — TPP (ionic, reversible, chitosan-only),
genipin (covalent, slow, chitosan-only), ECH/DVS (covalent, toxic,
agarose-only). STMP crosslinks **both** networks (phosphate diester on
agarose -OH, phosphoramide on chitosan -NH₂) in the same bead, in a
single triggerable reaction.

**Do not confuse STMP with STPP.** Sodium *Tri*metaphosphate (STMP,
Na₃P₃O₉, CAS 7785-84-4) is the *cyclic* trimer used here — covalent,
alkaline pH, triggerable. Sodium *Tripolyphosphate* (STPP, Na₅P₃O₁₀,
CAS 7758-29-4) is the *linear* ionic crosslinker already in EmulSim's
`tpp` entry — different chemistry entirely. Always check the CAS on
every reagent bottle.

### What you got in v9.1.2

- New primary (L3) crosslinker `stmp` and secondary (M2) crosslinker
  `stmp_secondary`, both visible in the `AGAROSE_CHITOSAN` family UI.
- A 102-line Appendix J.1.7 protocol card: materials, safety, the
  three-phase procedure, QC acceptance criteria (ICP-OES P content,
  FTIR P=O and P-O-C bands, swelling and modulus changes), and
  troubleshooting.
- UI guardrails for the homogeneity window and the STMP/STPP
  disambiguation.
- Zero new solver code paths; STMP reuses the existing
  `mechanism="hydroxyl"` dispatch (same as ECH, DVS, citric acid).

## What was in v9.1.1 ("Backlog burndown")

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
This release corresponds to tag `v9.2.0`.

## Licence

GPL-3.0. The full licence text ships as `LICENSE.txt` in the
install directory and is also at
<https://www.gnu.org/licenses/gpl-3.0.en.html>. The intellectual
property in this software belongs to Holocyte Pty Ltd; see
`LICENSE_AND_IP.txt` for the user agreement shown by the
installer.
