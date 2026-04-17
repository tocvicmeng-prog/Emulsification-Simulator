# Changelog

## v8.3.4 — Per-user install by default (fixes Access Denied) (2026-04-17)

Hotfix for a second install-time failure reported after v8.3.3:

```
[install] Creating virtual environment at .venv\
Error: [WinError 5] Access is denied: 'C:\Program Files\EmulSim\.venv'
[install] ERROR: venv creation failed.
```

### Root cause

The v8.3.2 / v8.3.3 Inno Setup script used
`DefaultDirName={autopf}\EmulSim` with
`PrivilegesRequiredOverridesAllowed=dialog`. On a UAC-elevated
install Inno Setup placed files into `C:\Program Files\EmulSim`,
but the `[Run]` post-install step (`install.bat`) executes in the
user's non-elevated context. Non-admin cannot create `.venv\`
inside `C:\Program Files\...`, so venv creation fails.

### Fix (Inno Setup script)

- `DefaultDirName={userpf}\EmulSim` — per-user Program Files
  (`%LOCALAPPDATA%\Programs\EmulSim`), always user-writable.
- `PrivilegesRequiredOverridesAllowed=dialog` removed — user can no
  longer accidentally elevate to a location where the post-install
  step will fail.
- `UsedUserAreasWarning=no` — suppresses the Inno warning that
  would otherwise trigger for an all-user-area install script.

### Fix (install.bat)

- Venv-creation failure now prints an actionable diagnostic:
  "directory not writable → uninstall and reinstall v8.3.4+ per-user,
  or right-click install.bat → Run as administrator".
- No more silent exit 3.

### Migration note for existing admin installs

If v8.3.2 or v8.3.3 was installed into `C:\Program Files\EmulSim`:

1. Uninstall (Control Panel → Apps → EmulSim, or the Start-Menu
   "Uninstall EmulSim" shortcut).
2. Download `EmulSim-8.3.4-Setup.exe` from the GitHub Release.
3. Double-click — it installs into `%LOCALAPPDATA%\Programs\EmulSim`
   without UAC. The post-install step completes cleanly.

### Version bumps

- `pyproject.toml`, `src/emulsim/__init__.py`, installer script,
  build helper: 8.3.3 → 8.3.4.

### Artefacts

- `release/EmulSim-8.3.4-Setup.exe` (2.54 MB)
- `release/EmulSim-8.3.4-Windows-x64.zip` (563 KB)
- `dist/emulsim-8.3.4-py3-none-any.whl`

Wheel contents unchanged from v8.3.2.

---

## v8.3.3 — Self-healing launch scripts (2026-04-17)

Hotfix for a dead-end user experience on first run: if the installer's
post-install step was skipped or failed silently (e.g. because
Python 3.11+ was not on PATH at install time), the launcher batch
files previously printed only "Installation not found. Run install.bat
first." and exited, with no actionable guidance.

### Fixed

- `release/.../launch_ui.bat` and `release/.../launch_cli.bat`:
  **self-healing**. On missing `.venv`, they now
  1. report the exact expected path,
  2. probe for `python` on `PATH` and show the detected version,
  3. if Python is absent, print a hyperlink to
     `https://www.python.org/downloads/windows/` and abort cleanly
     with a press-any-key,
  4. if Python is present, offer to run `install.bat --no-test`
     automatically and then continue to the launch,
  5. if setup fails, show the error code and keep the window
     open so the user sees the cause.
- `release/.../install.bat`: always `pause` on completion so the
  user sees the success / failure message. Honours
  `NONINTERACTIVE=1` when invoked from automation. Explicit error
  message + pause on pip-upgrade failure (previously exited 4
  silently).

### Changed (version bumps)

- `pyproject.toml`: 8.3.2 → 8.3.3.
- `src/emulsim/__init__.py.__version__`: 8.3.2 → 8.3.3.
- `installer/EmulSim.iss`, `installer/build_installer.bat`: all
  `8.3.2` references updated to `8.3.3`.

### Artefacts

- `release/EmulSim-8.3.3-Setup.exe` (2.54 MB) — Inno Setup wizard.
- `release/EmulSim-8.3.3-Windows-x64.zip` (563 KB) — portable.
- `dist/emulsim-8.3.3-py3-none-any.whl` (~408 KB) — wheel.

All three are identical in wheel contents to v8.3.2; only the
launcher batch files changed. Users who already have a working
v8.3.2 install can just replace `launch_ui.bat` / `launch_cli.bat` /
`install.bat` with the v8.3.3 versions.

### Smoke verified

Fresh temp venv + `pip install emulsim-8.3.3-py3-none-any.whl` +
`import emulsim` — works end-to-end.

---

## v8.3.2 — One-click Windows 11 installer (.exe) (2026-04-17)

Ships a proper Windows installer wizard as
`release/EmulSim-8.3.2-Setup.exe` (2.54 MB), attached to the
existing v8.3.2 GitHub Release alongside the portable ZIP.

### Added

- `installer/EmulSim.iss` — Inno Setup 6 script defining the full
  wizard:
  1. **EULA page** (`LICENSE_AND_IP.txt`) declaring: intellectual
     property rights belong to Holocyte Pty Ltd; software licensed
     under GPL-3.0; canonical source at
     `github.com/tocvicmeng-prog/Emulsification-Simulator`.
  2. **Python presence check** with a clickable hyperlink to
     `https://www.python.org/downloads/windows/` if Python 3.11+ is
     not on PATH.
  3. **File layout** — wheel, configs, docs, launcher batch files,
     LICENSE, README, INSTALL, RELEASE_NOTES all extracted under a
     single install directory.
  4. **Shortcuts** — Start-Menu group with Web-UI, CLI, Manual
     PDF, and Uninstall entries; optional desktop shortcut.
  5. **Post-install hook** — runs the bundled `install.bat` which
     creates a self-contained `.venv` and pip-installs the wheel
     with `[ui,optimization]` extras, with a smoke-pipeline check.
  6. **Uninstaller** — purges `.venv` before removing files.
- `installer/LICENSE_AND_IP.txt` — the EULA text shown on the
  installer's first page.
- `installer/build_installer.bat` — four-step build helper
  (wheel, stage, locate ISCC, compile).
- `installer/README.md` — documentation of the installer build and
  runtime behaviour.

### Changed

- `.gitignore` — now also excludes `installer/stage/` (transient
  build directory rebuilt by `build_installer.bat`).

### GitHub Release (v8.3.2)

Two assets now attached:

| Asset | Size | Audience |
|---|---|---|
| `EmulSim-8.3.2-Setup.exe` | 2.54 MB | End users (one-click wizard installer) |
| `EmulSim-8.3.2-Windows-x64.zip` | 561 KB | Power users (portable, script-based install) |

---

## v8.3.2 — Clean Windows 11 x64 install package (2026-04-17)

Ships a self-contained, dev-artifact-free Windows 11 x64 install
bundle as `release/EmulSim-8.3.2-Windows-x64.zip` (0.55 MB
compressed, 14 files). A fresh Windows machine with Python 3.11+
installed can extract the zip and run `install.bat` to get a
fully working simulator — UI, CLI, and programmatic API — in a
self-contained `.venv\` that leaves system Python untouched.

### Version bumps

- `pyproject.toml`: 0.1.0 → 8.3.2 (caught up with feature releases).
- `src/emulsim/__init__.py.__version__`: 0.1.0 → 8.3.2.

### Build artefacts

- `dist/emulsim-8.3.2-py3-none-any.whl` — rebuilt wheel covering
  the full v8.3 feature set (four polymer platforms, inverse
  design, digital twin, MD ingest, Unicode-safe PDF manual).
- `dist/emulsim-8.3.2.tar.gz` — source distribution.

### Release tree (`release/EmulSim-8.3.2-Windows-x64/`)

| File | Purpose |
|---|---|
| `install.bat` | Create `.venv\`, install wheel with `[ui,optimization]` extras, verify import, run smoke pipeline. Accepts `--core` / `--no-opt` / `--no-test` flags. |
| `launch_ui.bat` | Start the Streamlit UI at `http://localhost:8501`. |
| `launch_cli.bat` | Open a Command Prompt with the venv activated and `emulsim` on PATH. |
| `uninstall.bat` | Confirm-and-delete the `.venv\`. |
| `README.txt` | One-page quickstart. |
| `INSTALL.md` | Detailed install + troubleshooting guide (7 sections). |
| `RELEASE_NOTES.md` | User-facing summary of what's in 8.3.2. |
| `LICENSE.txt` | Software licence. |
| `wheels/emulsim-8.3.2-py3-none-any.whl` | The wheel (408 KB). |
| `configs/{default,fast_smoke,stirred_vessel}.toml` | Example configs. |
| `docs/User_Manual_First_Edition.{pdf,md}` | First Edition manual. |

### "Clean" guarantees (validated at zip time)

The zip-builder refuses to create the archive if any of these are
present in the release tree:

- `__pycache__/`, `.pyc`, `.pyo`
- `.pytest_cache/`, `.mypy_cache/`, `.ruff_cache/`
- `.git/`, `.venv/`
- `build/`, `dist/`, `output/`
- `.log` files

Not shipped in the release (kept in the dev repo only):

- Source tree (`src/` — replaced by the wheel).
- Test suite (`tests/`).
- Internal design docs (`docs/f1a_*`, `docs/f1b_*`, `docs/f1c_*`,
  `docs/f2_*`, `docs/f4b_*`, `docs/f5_*`, `docs/node31_*`,
  `docs/node32_*`, `docs/node30_31_*`).
- Full dev `CHANGELOG.md` (condensed into `RELEASE_NOTES.md`).
- `skills/`, `.claude/` agent infrastructure.

### Smoke verification (performed before the ZIP was cut)

- Fresh temp venv + `pip install wheels/emulsim-8.3.2-py3-none-any.whl`:
  install succeeded.
- `import emulsim; emulsim.__version__ == '8.3.2'` — OK.
- `run_pipeline()` default: returns a FullResult with
  `d32 = 18.08 µm` — end-to-end L1 → L4 pipeline executes cleanly
  on a fresh install.

### Final archive

`release/EmulSim-8.3.2-Windows-x64.zip` — 561 KB compressed,
648 KB uncompressed, 14 files. Ready for distribution.

---

## v8.3.2 — First Edition PDF Unicode font fix (2026-04-17)

Fixes black-square ("tofu") rendering of scientific Unicode glyphs
(α, ²⁺, ⌈⌉, χ, ∇, ∂, μ, π, ≥, etc.) in the First Edition PDF.

### Root cause

reportlab's built-in Type-1 fonts (Helvetica / Courier) cover only
WinAnsi / Latin-1. Any glyph outside that band — superscripts, Greek
letters, ceiling / floor brackets, mathematical operators — was
rendered as a black filled square. Visible examples reported by the
user: `Ca²⁺-alginate`, `CVaR_α = (1 / ⌈α·N⌉)`.

### Fix

- `docs/user_manual/build_pdf.py` — register DejaVu Sans and DejaVu
  Sans Mono (shipped with matplotlib, each ~6 000 Unicode glyphs
  covering full Greek, super/subscripts, mathematical operators,
  arrows) as TTF fonts via `reportlab.pdfbase.pdfmetrics.registerFont`
  / `registerFontFamily` at module import time. All body, heading,
  code, bullet, caption, table-cell, and page-footer styles now
  reference `DejaVuSans` / `DejaVuSansMono` instead of
  `Helvetica` / `Courier`. Falls back gracefully to the Type-1
  fonts if DejaVu is missing.
- `docs/user_manual/polysaccharide_microsphere_simulator_first_edition.pdf`
  — rebuilt. File size 55 KB → 143 KB (DejaVu TTFs now embedded).

### Verification

- Font cmap coverage confirmed for 20 scientific codepoints
  (U+00B2, U+207A, U+207B, U+03B1, U+03C7, U+03BC, U+2308, U+2309,
  U+00B7, U+2192, U+2207, U+2202, U+222B, U+03C0, U+00B0, U+2265,
  U+2264, U+00B1, U+00D7, U+00D8): **20/20 present**.
- Round-trip text extraction via pypdf finds all user-flagged
  phrases literally in the rebuilt PDF:
  - `Ca²⁺-alginate` present
  - `⌈α·N⌉` present
  - α, χ, ∇, ∂, ≥, π, · present.

---

## v8.3.1 — First Edition user manual (2026-04-17)

Ships the Polysaccharide-Based Microsphere Emulsification Simulator
**First Edition** instruction manual as Markdown + PDF, with an
upper-right download button wired into the Streamlit UI. The
manual is written for first-time users (downstream-processing
technicians, junior researchers) who have no prior experience in
microsphere fabrication or downstream processing.

### Added

- `docs/user_manual/polysaccharide_microsphere_simulator_first_edition.md`
  — the authoritative instruction manual. Three-part structure:
  1. **Getting Started** — what the simulator does, workflow
     overview with ASCII chart, five-minute quickstart.
  2. **Platform Catalogue** — the four supported polymer families
     (agarose-chitosan, alginate, cellulose NIPS, PLGA), a
     crosslinker / gelant selection table, the EDC/NHS COOH
     warning.
  3. **Appendices A–I** — detailed input requirements (all
     parameter tables with units + ranges + defaults), process
     steps, essential pre-run checklist, 15-question FAQ,
     architectural ideas + working principles, chemical / physical
     principles, formulas + theorems, six standard wet-lab
     protocols (agarose-chitosan/genipin, Ca²⁺-alginate external &
     internal, cellulose NaOH/urea, EDC/NHS coupling, PLGA
     solvent-evaporation), and a 17-row troubleshooting table.
- `docs/user_manual/build_pdf.py` — compact Markdown-to-PDF
  renderer using reportlab. Handles the Markdown subset used in
  the manual (headings, paragraphs, ordered / unordered lists,
  GitHub tables, fenced code blocks, inline `**bold**` /
  `*italic*` / `` `code` `` with underscore-safe code-span
  extraction via placeholders). Run
  `python docs/user_manual/build_pdf.py` to rebuild.
- `docs/user_manual/polysaccharide_microsphere_simulator_first_edition.pdf`
  — the built artefact (~55 KB, A4, page-footer on every page).

### Changed

- `src/emulsim/visualization/app.py` — the page title row now uses a
  two-column layout with the title on the left and a
  **Manual (PDF)** download button in the upper-right corner.
  Button serves the PDF via `st.download_button` when the file
  exists; falls back to a caption telling the user to run the
  build script if the PDF is absent.

### Dependencies

- `reportlab` added (auto-installed into the user's pip environment
  during the build step). No new runtime requirement for users who
  don't need to regenerate the PDF.

### Tests

- Quick regression of 66 targeted tests (F4-b, F5, F2, PLGA Phase 2,
  cellulose Phase 2/3, alginate L4) pass with 0 regressions after
  the UI edit.

---

## v8.3.0-alpha — Cluster F finish: F4-b CVaR + F5 MD ingest + F2 digital twin (2026-04-17)

Closes the three remaining Cluster F workstreams from the Node 32
roadmap. With this release, **every workstream in Cluster F has a
Phase 1 shipment**.

### F4-b — CVaR robust BO

- `OptimizationEngine(robust_cvar_alpha=α)` — applies Conditional
  Value-at-Risk aggregation over resamples. When both
  `robust_variance_weight` and `robust_cvar_alpha` are set, CVaR
  takes precedence.
- `emulsim design --robust-cvar-alpha α` CLI flag.
- Algorithm: ``CVaR_α = mean of the worst ⌈α·n⌉ resamples per
  objective dimension``. α → 1 recovers the sample mean; α → 0
  approaches the worst-case sample.
- Validation: α ∈ [0, 1]; α > 0 requires `robust_n_samples >= 2`
  and a `target_spec`.
- `docs/f4b_cvar_protocol.md` — full /architect protocol.
- `tests/test_f4b_cvar.py` — 11 tests (math, engine validation,
  CLI, precedence, head-to-head vs mean-variance).

### F5 — MARTINI MD parameter ingest

- `src/emulsim/md_ingest.py` — `MartiniRecord` dataclass + JSON
  load / save + `apply_chi_to_props(props, record)` for cellulose
  χ fields.
- JSON schema: required `source / system_description / beads / chi /
  diagnostics`; optional `paper_doi / notes`; forward-compat
  unknown top-level keys preserved in `record.extra`.
- Current mapping: `polymer_solvent / polymer_nonsolvent /
  solvent_nonsolvent` → `chi_PS_cellulose / chi_PN_cellulose /
  chi_SN_cellulose` on `MaterialProperties`. Non-cellulose fields
  never modified. Missing χ sub-keys leave fields untouched.
- Validation: NaN / inf χ rejected at load; negative χ allowed
  (physically valid for attractive mixing).
- `data/validation/md/example_martini_cellulose.json` — reference
  fixture for tests and user-authoring template.
- `docs/f5_md_ingest_protocol.md` — full /architect protocol.
- `tests/test_md_ingest.py` — 11 tests (load, missing keys,
  extra-keys preservation, partial χ, apply to props, non-cellulose
  fields untouched, non-finite guards, negative χ, round-trip).

### F2 — Digital twin EnKF replay (Phase 1)

- `src/emulsim/digital_twin/enkf.py` — stochastic Ensemble Kalman
  Filter (Evensen 1994) `enkf_update(x, y_fc, y_obs, R, rng,
  inflation)`. Scalar observations only in Phase 1; multiplicative
  prior inflation optional.
- `src/emulsim/digital_twin/replay.py` — `run_replay(trace, x0,
  state_transition, observation_operator, ...)` walks forward
  through a `DigitalTwinTrace`, applies EnKF at each observation,
  returns `ReplayResult` with per-observation mean / std / optional
  full ensemble + a `DigitalTwin.EnKFReplay` SEMI_QUANTITATIVE
  manifest.
- `src/emulsim/digital_twin/schema.py` — `DigitalTwinTrace` +
  `Observation` dataclasses + JSON load / save (sorts observations
  by time on load).
- `src/emulsim/digital_twin/__init__.py` — module exports.
- `docs/f2_digital_twin_protocol.md` — full /architect protocol.
- `tests/test_digital_twin.py` — 11 tests (EnKF linear-Gaussian
  convergence, zero-noise collapse, inflation grows spread, EnKF
  input validation, trace round-trip + time-ordering on load,
  replay trajectory shape, empty-trace passthrough, multi-step
  spread shrinkage).

### Tests

- 33 new tests across F4-b (11) + F5 (11) + F2 (11). 218 targeted
  regression tests pass (PLGA Phase 1/2 + cellulose Phase 1/2/3 +
  internal-gelation + alginate 2a/b/c + EDC/NHS + UQ unified + CLI
  v7 + inverse-design + F4-b + F5 + F2) with 0 regressions.

### Footprint

- F4-b: ~220 LOC (engine edit + CLI + tests).
- F5: ~465 LOC (module + fixture + tests + docs).
- F2: ~895 LOC (schema + enkf + replay + __init__ + tests + docs).
- Total: ~1580 LOC added this turn.

### Cluster F status

All Cluster F workstreams from the Node 32 v8.0 roadmap are at
least Phase-1 shipped:

| Workstream | Status |
|---|---|
| F1-a Alginate | ✓ fully wired (v8.0-rc2) |
| F1-b Cellulose NIPS | ✓ fully wired (v8.1-beta) |
| F1-c PLGA | ✓ fully wired (v8.2-beta) |
| F2 Digital twin (EnKF replay) | ✓ Phase 1 shipped (v8.3-alpha) |
| F3-a Inverse design | ✓ complete (v8.0-alpha) |
| F3-b/c BO engine + CLI | ✓ complete (v8.0-alpha) |
| F4-a Robust BO (mean-variance) | ✓ complete (v8.0-alpha) |
| F4-b Robust BO (CVaR) | ✓ complete (v8.3-alpha) |
| F5 MD ingest (MARTINI) | ✓ Phase 1 shipped (v8.3-alpha) |

### Still deferred (Phase-2+ items, each needs fresh /architect kickoff)

- F2: vector observations (matrix R), square-root / deterministic
  EnKF variants, online polling adapter, MPC layer, identifiability
  diagnostics.
- F5: tabulated U(r) pair-potential ingestion, automatic
  MARTINI ↔ EmulSim bead-type mapping, CalibrationStore integration,
  reverse-direction emit.
- F4-b: automatic α selection (tail-risk auto-tune).
- PLGA moving-boundary ALE solver + Fujita `D(phi)`.
- v7.0 release still blocked on Study A wet-lab data.

---

## v8.2.0-beta — F1-c Phase 2: PLGA orchestrator + CLI + TOML (2026-04-17)

Closes F1-c. **All three Cluster F platforms (alginate, cellulose,
PLGA) are now fully wired end-to-end** — orchestrator dispatch, CLI
flags, and TOML config keys. F1 (multi-platform microsphere family)
is complete at the protocol-scope level.

### Added

- `PipelineOrchestrator._run_plga(...)` — mirrors `_run_cellulose` /
  `_run_alginate`. Applies `params.formulation.plga_grade` preset to
  MaterialProperties before the L2 solver runs. Skips L2a timing and
  L3 crosslinking (PLGA microspheres are glassy / physically
  entangled, not crosslinked); stubs `CrosslinkingResult`. Emits
  summary.json with `polymer_family = "plga"`, L2 diagnostics
  (`phi_plga_mean_final`, `t_vitrification_s`,
  `skin_thickness_proxy_m`, `R_shrunk_m`), and the L4 modulus.
- `run_single` branch: `props.polymer_family == PolymerFamily.PLGA`
  routes to `_run_plga`. Placed immediately after the cellulose
  branch for symmetry.
- `emulsim run --plga-grade {50_50 | 75_25 | 85_15 | pla}` CLI flag.
  Packs the grade's 4 PLGA-specific fields into `props_overrides`.
  Meaningful with `--polymer-family plga`; prints a one-line
  confirmation unless `--quiet`.
- TOML `[formulation].plga_grade = "..."` unpacks directly into the
  existing `FormulationParameters.plga_grade` field (shipped in
  F1-c Phase 1); orchestrator resolves at run time via
  `properties.plga_defaults.apply_preset`.

### Changed

- `orchestrator.py` — new imports (`solve_solvent_evaporation`,
  `solve_mechanical_plga`), new PLGA branch, new `_run_plga` method
  (~115 LOC).
- `__main__.py` — `--plga-grade` flag + `_cmd_run` hook to expand
  preset into `props_overrides`.
- `config.py` — no changes needed; TOML key unpacks naturally via
  the existing `plga_grade` field.

### Tests

- `tests/test_plga_phase2.py` — 12 tests:
  - Orchestrator dispatch (3): PLGA routes to `_run_plga`, summary.json
    records `polymer_family`, full pipeline reports
    SEMI_QUANTITATIVE end-to-end.
  - Preset application (2): orchestrator patches props (85:15 K_glassy
    = 1 × 10⁹ Pa verified); unknown grade raises `KeyError`.
  - TOML (2): `plga_grade` key unpacks; absent defaults to empty string.
  - CLI (2): all 4 choices in shipped parser source; argparse accepts
    the full flag invocation.
  - End-to-end sanity (3): non-zero G_DN; switching grade gives ≥
    1.5× modulus spread; L3 `p_final = 0` (stubbed).
- 185 targeted regression tests pass (PLGA Phase 1 + Phase 2 +
  cellulose Phases 1/2/3 + internal-gelation + alginate 2a/b/c +
  EDC/NHS + UQ unified + CLI v7 + inverse-design) with 0 regressions.

### Footprint

- New LOC: ~115 (orchestrator `_run_plga`) + ~30 (CLI) + ~340
  (tests) ≈ 485 LOC. Under the protocol's ~430 LOC estimate by a
  hair (simple wiring + trivial TOML).
- Cumulative F1 footprint across all three platforms: **~5000 LOC**.

### Cluster F status after v8.2.0-beta

| Platform | Programmatic | Orchestrator | CLI | TOML | Presets |
|---|---|---|---|---|---|
| Agarose-chitosan (default) | ✓ | ✓ | ✓ | ✓ | (built-in) |
| Alginate | ✓ | ✓ | ✓ | ✓ | 2 gelants |
| Cellulose NIPS | ✓ | ✓ | ✓ | ✓ | 4 solvents |
| PLGA solvent evap | ✓ | ✓ | ✓ | ✓ | 4 grades |

**F1 complete.** Other Cluster F workstreams remain un-started and
need their own /architect kickoffs:

- **F2 digital twin** (EnKF replay harness) — scoped in Node 32
  roadmap, protocol not drafted.
- **F4-b CVaR acquisition** — deferred v8.0 polish; trivial variant
  of F4-a once the resample strategy is finalised.
- **F5 MD parameter ingest** (MARTINI — ingest-only default) —
  scoped in Node 32 roadmap, protocol not drafted.

### Still deferred

- Moving-boundary ALE solver for PLGA shrinking-droplet correction
  (R5 from F1-c protocol).
- Fujita concentration-dependent `D(phi)` for late-stage PLGA
  evaporation dynamics.
- v7.0 release still blocked on Study A wet-lab data.

---

## v8.2.0-alpha — F1-c Phase 1: PLGA solvent-evaporation L2 + L4 + 4 grades (2026-04-17)

First shipment of Cluster F platform #3. Adds a PLGA /
solvent-evaporation L2 solver + Gibson-Ashby L4 modulus + full 4-grade
registry (PLGA 50:50 / 75:25 / 85:15 / PLA). Ships as
**programmatic API only**; orchestrator dispatch / CLI / TOML = F1-c
Phase 2.

### Added

- `docs/f1c_plga_protocol.md` — full /architect protocol doc (§1 scope,
  §2 mechanism + lit anchors, §4 algorithm, §5 4-grade parameter
  table, §6 16 test cases, §7 risks, §8 G1 12/12 for Phase 1,
  §9 execution plan). Matches the f1a / f1b protocol pattern.
- `src/emulsim/level2_gelation/solvent_evaporation.py` (~330 LOC):
  1D spherical Fickian DCM-depletion solver.
  - State: ``phi_DCM(r, t)`` single field; ``phi_PLGA = 1 − phi_DCM``
    algebraic.
  - Dirichlet sink at ``r = R`` (``phi_DCM_eq ≈ 0.005`` for DCM/water),
    symmetry at ``r = 0``.
  - BDF time integrator with dense-output vitrification-time probe.
  - Approximations: fixed droplet radius (moving-boundary ALE = Phase 2
    refinement); constant D (Fujita ``D(phi)`` = Phase 2). Both
    flagged in manifest assumptions.
  - Emits SEMI_QUANTITATIVE `GelationResult` tagged
    `L2.Gelation.SolventEvaporationPLGA` with
    `phi_plga_mean_final / phi_dcm_mean_final / t_vitrification /
    skin_thickness_proxy / core_porosity_proxy / R_shrunk_m`
    diagnostics.
- `src/emulsim/level4_mechanical/plga.py` (~130 LOC):
  `plga_modulus(phi_mean, G_glassy, n_plga)` Gibson-Ashby power law +
  `solve_mechanical_plga(...)` wrapper emitting
  `L4.Mechanical.PLGAGibsonAshby` SEMI_QUANTITATIVE
  `MechanicalResult` with `network_type="glassy_polymer"` and
  `model_used="plga_gibson_ashby"`.
- `src/emulsim/properties/plga_defaults.py` (~160 LOC):
  `PLGAGradeProfile` dataclass + `PLGA_GRADE_PRESETS` registry with
  **all four** grades populated (`50_50`, `75_25`, `85_15`, `pla`).
  Data sourced from Wang 1999 (D_DCM), Park 1998 (T_g, G_glassy),
  Freitas 2005 (process parameters). `apply_preset(props, grade)`
  helper mirrors the alginate / cellulose pattern. Phase 3 is
  effectively eliminated by shipping all 4 presets up front.
- `MaterialProperties` gains 4 PLGA-specific fields
  (`D_DCM_plga`, `phi_DCM_eq`, `G_glassy_plga`, `n_plga_modulus`).
- `FormulationParameters` gains `phi_PLGA_0` (initial polymer volume
  fraction in the droplet; default 0 = not-PLGA) and `plga_grade`
  (Phase 2 preset-selector field; default empty = skip).

### Tests

- `tests/test_plga_phase1.py` — 25 tests covering:
  - **Protocol §6 test 1**: monotone DCM depletion (with fixed
    transient-regime probe times)
  - **Protocol §6 test 2**: Dirichlet sink drives `phi_PLGA → 1` at
    long time
  - **Protocol §6 test 3**: early-regime √t front scaling (log-log
    slope 0.5 ± 0.15)
  - **Protocol §6 test 4**: 4× D_DCM gives earlier vitrification
  - **Protocol §6 test 5**: Gibson-Ashby `G ∝ phi^n` scaling (n ∈
    {1.5, 2.0, 2.5}); dense limit recovers `G_glassy`
  - **Protocol §6 test 6**: zero PLGA → UNSUPPORTED manifest + zero
    modulus
  - **Protocol §6 test 7**: L2 + L4 both SEMI_QUANTITATIVE; full
    diagnostic key presence
  - **Protocol §6 test 8**: `apply_preset` patches 4 PLGA fields;
    4 grades registered; physical-plausibility check on every grade
    (L_fraction, M_n, T_g, D, G_glassy, n ranges); switching grade
    gives ≥ 1.4× modulus spread
  - Edge cases: `plga_modulus` zero/negative inputs
  - Input validation: negative R, tiny grid, phi_0 out of [0, 1],
    negative time, non-positive D_DCM
  - Mass-conservation post-processing: `R_shrunk = R_0 · phi_0^(1/3)`
    at long time
- 173 targeted regression tests pass (PLGA Phase 1 + cellulose
  Phase 1/2/3 + internal-gelation + alginate 2a/b/c + EDC/NHS + UQ
  unified + CLI v7 + inverse-design) with 0 regressions.

### Footprint

- New LOC: ~330 (solver) + ~130 (L4) + ~160 (defaults) + ~25
  (datatypes) + ~390 (tests) ≈ 1035 LOC. Slightly over the protocol's
  ~910 LOC estimate because all 4 grade presets shipped in Phase 1
  instead of 3 being deferred to Phase 3.

### Still deferred

- **F1-c Phase 2** (~430 LOC, 1–2 sessions): orchestrator
  `_run_plga` branch (mirror `_run_cellulose`), `--polymer-family
  plga --plga-grade <name>` CLI surface, `[formulation].plga_grade`
  TOML key application, 12 integration tests.
- **F1-c Phase 3**: absorbed into Phase 1 (all 4 presets shipped).
- **Moving-boundary ALE solver** for shrinking-droplet correction
  (R5 from the protocol). Current fixed-R approximation reports
  `R_shrunk` as a post-processing diagnostic so users can plot the
  true final sphere.
- **Fujita concentration-dependent D**: v1 uses constant D; late-time
  (phi > 0.8) dynamics are order-of-magnitude-right, not quantitative.
- v7.0 release still blocked on Study A wet-lab data.

---

## v8.1.0-beta — F1-b Phases 2 + 3: cellulose orchestrator + all 4 solvents (2026-04-17)

Closes F1-b. Cellulose NIPS is now a first-class user-facing platform
(matches alginate surface area). All four solvent-system presets are
populated; the orchestrator dispatches `PolymerFamily.CELLULOSE`
through a dedicated `_run_cellulose` sub-pipeline; TOML and CLI flags
expose both family selection and solvent selection.

### Added

- `PipelineOrchestrator._run_cellulose(...)` — mirrors
  `_run_alginate`. Applies a solvent-system preset (if declared on
  `params.formulation.solvent_system`) to MaterialProperties before
  solving L2 NIPS. Skips L2a timing and L3 crosslinking (NIPS IS the
  gelation); stubs `CrosslinkingResult`. Emits summary.json with
  `polymer_family = "cellulose"`, `phi_mean_final`,
  `bicontinuous_score`, `demixing_index`, and the L4 modulus.
- `run_single` branch: `props.polymer_family == PolymerFamily.CELLULOSE`
  routes to `_run_cellulose`. Placed immediately after the alginate
  branch for symmetry.
- `FormulationParameters.solvent_system: str = ""` — TOML key
  `[formulation].solvent_system = "naoh_urea"` unpacks directly into
  this field, then the orchestrator resolves it at run time via
  `properties.cellulose_defaults.apply_preset`.
- `emulsim run --cellulose-solvent {naoh_urea | nmmo | emim_ac |
  dmac_licl}` CLI flag — packs the preset's 9 cellulose-specific
  fields into `props_overrides` before `run_single`. Meaningful with
  `--polymer-family cellulose`; prints a one-line confirmation
  (`Cellulose solvent preset: ...`) unless `--quiet`.
- Three new presets in `src/emulsim/properties/cellulose_defaults.py`:
  - **NMMO** (Lyocell, 80 wt% aq., T = 90 °C, higher N_p = 500,
    K_cell = 8 × 10⁵ Pa; Lenzing system).
  - **EMIM-Ac** (1-ethyl-3-methylimidazolium acetate, T = 80 °C,
    lowest χ_PS = 0.38; Swatloski 2002 IL system).
  - **DMAc/LiCl** (McCormick analytical system, T = 60 °C activation,
    K_cell = 6 × 10⁵ Pa). All values from the literature anchors in
    `docs/f1b_cellulose_nips_protocol.md` §5.

### Changed

- `orchestrator.py` — new imports (`solve_nips_cellulose`,
  `solve_mechanical_cellulose`), new CELLULOSE branch, new
  `_run_cellulose` method (~110 LOC).
- `__main__.py` — `--cellulose-solvent` flag + `_cmd_run` hook to
  expand preset into `props_overrides`.
- `config.py` — TOML `[formulation].solvent_system` unpacks naturally
  via the new `FormulationParameters.solvent_system` field. No
  special-case parsing; validation deferred to solver-time
  `apply_preset(...)`.

### Tests

- `tests/test_cellulose_phase2_phase3.py` — 14 tests:
  - Orchestrator dispatch (3): CELLULOSE routes to `_run_cellulose`,
    summary.json records polymer_family, full pipeline reports
    SEMI_QUANTITATIVE end-to-end.
  - TOML wiring (2): `solvent_system` key unpacks, absent key defaults
    to empty string.
  - Solvent preset application (2): orchestrator patches props so L4
    K_cell matches the NMMO preset (8 × 10⁵ Pa); unknown preset raises
    `KeyError`.
  - CLI (2): `--cellulose-solvent` flag in shipped parser, argparse
    accepts all 4 choices.
  - Registry (3): all 4 presets registered, each passes physical
    plausibility (χ_PN > χ_PS, D in bulk-water range, N_p in DP range,
    K_cell in 10⁴–10⁷ Pa band), water is the default non-solvent for
    all.
  - Diagnostics differentiation (1): switching preset changes G_DN by
    > 1.5× (spans real range).
  - Argparse rejection (1): unknown preset exits non-zero.
- 152 targeted regression tests pass (Phase 1 + Phase 2/3 +
  internal-gelation + alginate 2a/b/c + EDC/NHS + UQ unified + CLI v7
  + parallel MC + inverse-design) with 0 regressions.

### Footprint

- New LOC: ~110 (orchestrator `_run_cellulose`) + ~15 (CLI) + ~1
  (config.py, just the TOML doc comment) + ~160 (3 new presets) +
  ~330 (tests) ≈ 615 LOC. Cumulative F1-b footprint (Phases 1 + 2 +
  3): ~1525 LOC, a little under the 2000 LOC protocol budget because
  Phase 2 config wiring naturally unpacks via the single
  `solvent_system` field rather than needing a dedicated parser.

### Still deferred

- **v7.0 release** remains blocked on Study A wet-lab data.
- **F1-c PLGA solvent evaporation** still unscoped — a fresh
  /architect protocol is the natural next step if commercial
  prioritisation calls for it.

---

## v8.1.0-alpha — F1-b Phase 1: cellulose NIPS L2 + L4 + NaOH/urea (2026-04-17)

First shipment of Cluster F platform #2. Adds a cellulose /
non-solvent-induced phase separation (NIPS) L2 solver + L4 modulus +
NaOH/urea parameter preset. Ships as **programmatic API only** for now
— orchestrator dispatch / TOML config / CLI flags land in F1-b Phase 2.

### Added

- `src/emulsim/level2_gelation/nips_cellulose.py` (~380 LOC):
  1D spherical ternary Cahn-Hilliard + Fickian coupled-PDE solver.
  - State: `phi(r, t)` cellulose + `s(r, t)` solvent, `n = 1-phi-s`
    non-solvent (algebraic).
  - Flory-Huggins free energy with χ_PS / χ_PN / χ_SN.
  - Cahn-Hilliard gradient-energy regularisation on `mu_phi`.
  - Dirichlet bath BC at `r = R` (pure water by default), symmetry
    at `r = 0`.
  - 1 % noise on initial `phi` breaks spherical symmetry so spinodal
    decomposition can develop.
  - BDF time integration; clipped log arguments protect against
    spinodal excursions outside the physical simplex.
  - Emits SEMI_QUANTITATIVE `GelationResult` tagged
    `L2.Gelation.NIPSCellulose` with
    `phi_mean_final / s_mean_final / n_mean_final / phi_std_final /
    bicontinuous_score / demixing_index` diagnostics.
- `src/emulsim/level4_mechanical/cellulose.py` (~140 LOC):
  `cellulose_modulus(phi_mean, K_cell, alpha_cell)` power-law +
  `solve_mechanical_cellulose(...)` wrapper that emits
  `L4.Mechanical.CelluloseZhang2020` SEMI_QUANTITATIVE
  `MechanicalResult` with `network_type="physical_entangled"` and
  `model_used="cellulose_zhang2020"`.
- `src/emulsim/properties/cellulose_defaults.py` (~100 LOC):
  `CelluloseSolventPreset` dataclass + `CELLULOSE_SOLVENT_PRESETS`
  registry + `apply_preset(props, name)` helper. NaOH/urea preset
  (Zhang Lab Wuhan) populated from Lindman 2010 / Xu 2010 / Zhang
  2020. NMMO, EMIM-Ac, DMAc/LiCl stubs ship in F1-b Phase 3.
- `MaterialProperties` gains 9 cellulose-specific fields
  (`N_p_cellulose`, `chi_{PS,PN,SN}_cellulose`,
  `D_{solvent,nonsolvent}_cellulose`, `kappa_CH_cellulose`,
  `K_cell_modulus`, `alpha_cell_modulus`) — defaults match the
  NaOH/urea preset.
- `FormulationParameters.phi_cellulose_0` — initial cellulose volume
  fraction (default 0 = not-cellulose).

### Tests

- `tests/test_cellulose_nips_phase1.py` — 19 tests covering:
  - Protocol §6 test 2: ternary mass conservation (`phi + s + n = 1`)
  - Protocol §6 test 4: water-bath driven demixing (`phi_std` grows)
  - Protocol §6 test 7: modulus scaling `G ∝ phi^α`
  - Protocol §6 test 8: zero cellulose → zero gel / zero modulus
  - Protocol §6 test 10: SEMI_QUANTITATIVE manifests on both L2 and L4
  - NaOH/urea preset registry + `apply_preset` patching 9 fields
  - L4 modulus edge cases (zero phi / zero prefactor / negative phi)
  - Solver input validation (R ≤ 0, n_r < 8, phi_0 out of [0, 1],
    negative time / noise)
- 113 targeted regression tests pass (Phase 1 + internal-gelation +
  alginate Phase 2a/b/c + EDC/NHS + UQ unified + CLI v7) with 0
  regressions — the datatypes extensions did not perturb any existing
  consumers.

### Footprint

- New LOC: ~380 (solver) + ~140 (L4) + ~100 (defaults) + ~15
  (datatypes) + ~275 (tests) ≈ 910 LOC. On target for the protocol's
  Phase 1 ~900 LOC estimate.

### Still deferred (F1-b Phase 2 and 3)

- F1-b Phase 2: orchestrator `_run_cellulose` branch +
  `--polymer-family cellulose` CLI flag + `[formulation].solvent_system
  = "naoh_urea"` TOML key + 5 integration tests. ~500 LOC, 1–2
  sessions.
- F1-b Phase 3: NMMO / EMIM-Ac / DMAc/LiCl preset populations + 4
  solvent-dependence tests. ~350 LOC, 1 session.
- v7.0 release still blocked on Study A wet-lab data.

---

## v8.0.0-rc2 — Coupled GDL/CaCO₃ internal-release solver + F1-b protocol (2026-04-17)

Closes the last F1-a Phase 2c deferred item (the coupled
GDL/CaCO₃/alginate solver replacing the lumped-parameter exponential
approximation) and publishes the /architect protocol for F1-b
(cellulose NIPS), tee-ing up the next platform without committing to
implementation in this session.

### Added

- `solve_internal_gelation(params, props, *, R_droplet, C_CaCO3_0,
  L_GDL_0, k_hyd, k_diss, n_r, time, rtol, atol)` in
  `src/emulsim/level2_gelation/ionic_ca.py` — coupled ODE+PDE solver
  for homogeneous internal-release alginate gelation. State:
  - 3 spatially-uniform scalars (GDL, gluconic acid ≈ [H⁺], CaCO₃)
  - 3 radial fields (Ca²⁺, guluronate, egg-box crosslink density)
  - 6 coupled rate equations implementing GDL hydrolysis (Draget
    1997 k_hyd = 1.5 × 10⁻⁴ /s), CaCO₃ dissolution (Plummer 1978,
    Pokrovsky & Schott 2002 k_diss = 1 × 10⁻² m³/(mol·s)), and the
    existing Ca²⁺-guluronate egg-box binding.
  - No-flux outer BC (sealed droplet) + symmetry inner BC.
  - Emits SEMI_QUANTITATIVE `GelationResult` tagged
    `L2.Gelation.IonicCaInternalRelease` with `X_cov` homogeneity
    metric in diagnostics.
- `tests/test_internal_gelation.py` — 11 tests covering:
  - Schema + stoichiometric default for `L_GDL_0 = 2 × C_CaCO3_0`
  - First-order GDL hydrolysis decay (theory vs numerics within 2 %)
  - Monotone conversion
  - Zero-CaCO₃ and zero-alginate sanity
  - Ca²⁺ mass-balance under the no-flux BC
  - Homogeneity: internal-release CoV(X) < shrinking-core CoV(X) at
    matched Ca²⁺ budget and bead size (confirms the Draget 1997
    uniform-gel claim in simulation)
  - Input validation (negative R, negative CaCO₃, tiny grid)
- `docs/f1b_cellulose_nips_protocol.md` — full /architect protocol
  for the cellulose non-solvent-induced phase separation platform
  (F1-b). Includes: NIPS mechanism summary, 4-solvent parameter
  table (NaOH/urea, NMMO, EMIM-Ac, DMAc/LiCl), Cahn-Hilliard + ternary
  coupled-diffusion solver algorithm, 13 test cases, G1 gate status
  (10/12), and a 3-phase execution plan (4-6 fresh sessions total).

### Tests

- 11 new internal-gelation tests pass. 51 targeted regression
  (alginate Phase 2a + 2b + 2c + internal-release + CLI) with 0
  regressions.

### Still deferred

- **F1-b cellulose NIPS implementation** — protocol ready at
  `docs/f1b_cellulose_nips_protocol.md`; waits on fresh session(s)
  for the ~2000 LOC solver + tests.
- **F1-c PLGA solvent evaporation** — still un-scoped.
- **v7.0 release** — blocked on Study A wet-lab data.

---

## v8.0.0-rc1 — F1-a gelant preset wiring (2026-04-17)

Polish pass over v8.0.0-beta: the alginate reagent library shipped in
Phase 2c is now a first-class runtime input. Users can select a gelant
preset via `--gelant` on the CLI or `gelant = "..."` in the TOML
`[formulation]` section, and the simulator auto-wires the profile's
effective Ca²⁺ concentration into `FormulationParameters.c_Ca_bath`
using the current `t_crosslink` (static for external bath, saturating
exponential for internal release).

### Added

- `emulsim run --gelant {cacl2_external | gdl_caco3_internal}` CLI
  flag. When set, prints
  ``Gelant preset: <name> (c_Ca_bath = X mol/m³ at t_crosslink = Y s)``
  and overrides `formulation.c_Ca_bath` before the orchestrator runs.
- `[formulation].gelant = "<name>"` TOML key. Consumed by
  `load_config()` before unpacking the formulation section —
  `gelant` is a preset selector, not a persistent dataclass field.
  Unknown names raise `ValueError` with the list of available presets.
- 5 new tests in `tests/test_alginate_phase2c.py` (`TestGelantPreset`)
  covering external-bath static wiring, internal-release
  time-saturation, unknown-gelant rejection, and CLI argparse.

### Changed

- `_cmd_run` consults `GELANTS_ALGINATE` + `effective_bath_concentration`
  when `--gelant` is set, **after** `--polymer-family` has applied but
  **before** the orchestrator is instantiated.

### Tests

- 51 targeted regression pass (Phase 2c + Phase 2a/2b alginate + CLI),
  0 regressions. New Phase 2c total: 20 tests.

---

## v8.0.0-beta — F1-a Phase 2c: Alginate reagent library, TOML config, CLI (2026-04-17)

Closes the remaining three protocol §6 tests for the alginate platform
and exposes alginate as a first-class user-facing surface (CLI flag,
TOML config, reagent library). With Phase 2c in, alginate is no longer
a programmatic-API-only feature — end users can run
`python -m emulsim run --polymer-family alginate` against a TOML-defined
formulation and get the full L1 → L2-ionic-Ca → L4-Kong pipeline.

### Added

- `src/emulsim/reagent_library_alginate.py` — new
  `AlginateGelantProfile` dataclass and `GELANTS_ALGINATE` dict with two
  canonical entries:
  - **`cacl2_external`** — 100 mM CaCl₂ bath, shrinking-core mode,
    baseline for emulsification + drop-bath processes.
  - **`gdl_caco3_internal`** — glucono-δ-lactone + CaCO₃ in-situ
    release, lumped first-order release rate `k_release = 1.5e-4 s⁻¹`
    from Draget 1997.
  - `effective_bath_concentration(profile, t_end)` helper returns the
    static bath concentration for external mode and
    `C_source·(1 − exp(−k·t))` for internal mode.
- CLI `python -m emulsim run --polymer-family {agarose_chitosan |
  alginate | cellulose | plga}` flag routes `run_single` to the
  matching L2/L4 solver pair via `props_overrides={"polymer_family":
  ...}`.
- `tests/test_alginate_phase2c.py` — 15 tests covering:
  - **Protocol §6 test 3**: √t shrinking-core scaling
    (log-log slope of `X_mean` vs `t` is 0.5 ± 0.15 in the early
    diffusion-limited regime; R=1 mm, t=[10, 40, 160] s).
  - **Protocol §6 test 9**: TOML round-trip of
    `polymer_family = "alginate"` (and unknown-family rejection).
  - **Protocol §6 test 10**: L2 + L4 manifests both report
    SEMI_QUANTITATIVE when routed through the orchestrator, and
    RunReport's `min_evidence_tier` reflects that.
  - Alginate gelant library smoke (both modes, saturating release,
    unknown-mode ValueError).
  - CLI `--polymer-family alginate` argparse acceptance.

### Changed

- `src/emulsim/config.py` / `load_properties()` — accepts top-level
  scalar keys in addition to nested sections and coerces
  `polymer_family` strings to `PolymerFamily` enum members before
  constructing `MaterialProperties`. Unknown family names raise
  `ValueError` cleanly (no silent fall-through).
- `src/emulsim/__main__.py` / `_cmd_run` — `--polymer-family` override
  builds a `props_overrides` dict and passes it to
  `orchestrator.run_single`. Behaviour unchanged when the flag is
  omitted.

### Footprint

- New LOC: ~180 (reagent library) + ~20 (config) + ~15 (CLI) + ~270
  (tests) ≈ 485 LOC, within the protocol estimate for Phase 2c.
- Cumulative F1-a footprint (Phase 2a + 2b + 2c): ~1505 LOC of a
  projected ~1900 LOC; the ~400 LOC gap is cellulose-NIPS / PLGA
  scaffolding that was never in F1-a scope anyway.

### Tests

- 15 new Phase 2c tests; 115 targeted regression tests pass (alginate
  Phase 2a + 2b + 2c, EDC/NHS, UQ unified, UQ panel, inverse-design
  objectives + engine, parallel MC, CLI v7) with 0 regressions.

### Known limitations / still deferred

- Internal GDL/CaCO₃ mode uses the lumped-parameter
  `C_eff = C_source·(1 − exp(−k·t))` approximation; a fully coupled
  GDL + CaCO₃ + alginate solver is a Phase 3 follow-up if users
  request homogeneity predictions.
- Reagent library is surfaced as a module-level dict; a full
  `emulsim run --gelant cacl2_external` CLI flag that wires the
  profile into `FormulationParameters` automatically is a trivial
  v8.0-rc polish item.
- v7.0 release remains blocked on Study A wet-lab data for Node 21
  L1 PBE recalibration (unchanged from v7.0.1).

---

## v8.0.0-alpha — F1-a Phase 2b: Alginate L4 + orchestrator dispatch (2026-04-17)

Completes the functional alginate pipeline. `python`-level users can
now run a full L1 → L2 ionic-Ca → L4 Kong-2004 modulus pipeline for
alginate microspheres via `PipelineOrchestrator.run_single(params,
props_overrides={"polymer_family": PolymerFamily.ALGINATE, ...})`.

### Added

- `FormulationParameters.c_alginate` (kg/m³, default 0) +
  `FormulationParameters.c_Ca_bath` (mol/m³, default 100 mM CaCl₂).
  Zero `c_alginate` transparently falls back to the `c_agarose` slot
  for Phase 2a backward-compat.
- `src/emulsim/level4_mechanical/alginate.py` — Kong 2004 empirical
  modulus with `alginate_modulus(c, f_G, X_mean, K, n)` and
  `solve_mechanical_alginate(params, props, gelation, R_droplet=)`.
  Emits SEMI_QUANTITATIVE-tier `MechanicalResult` with
  `network_type="ionic_reinforced"` and `model_used="alginate_kong2004"`.
- `PipelineOrchestrator._run_alginate(...)` sub-pipeline: branches off
  `run_single` when `props.polymer_family == PolymerFamily.ALGINATE`.
  Skips L2a timing and L3 crosslinking (ionic gelation IS the
  crosslinking); stubs `CrosslinkingResult` to preserve the FullResult
  schema; records `polymer_family` in summary.json + RunReport
  diagnostics.
- `tests/test_alginate_l4_and_pipeline.py` — 7 tests covering:
  - `alginate_modulus` unit scaling in c² and f_G² (protocol §6 tests
    6, 7)
  - incomplete-gelation modulus reduction
  - `solve_mechanical_alginate` schema + zero-alginate edge case
  - full-pipeline orchestrator dispatch with `PolymerFamily.ALGINATE`
    (protocol §6 test 11)

### Tests

- 96/96 targeted regression tests pass in 88 s across F1-a Phase 2a/2b,
  F3 (inverse design + engine + CLI), F4-a, Node 30 / 31 / 30b, and
  the CLI contract. 0 regressions.

### Still deferred to F1-a Phase 2c / v8.0-beta

- Reagent library entries for CaCl₂ + internal GDL/CaCO₃ gelation.
- `config.py` TOML parser support for `polymer_family = "alginate"`
  (users currently set via `props_overrides` programmatically).
- The three remaining protocol §6 tests (§6 test 3 √t shrinking-core
  scaling, §6 test 9 TOML round-trip, §6 test 10 manifest-tier
  reporting from the orchestrator).
- CLI `python -m emulsim run --polymer-family alginate` surface.

### Footprint

- **Added:** ~380 LOC (L4 alginate + orchestrator branch + 7 tests) on
  top of Phase 2a's 640 LOC → cumulative F1-a footprint ~1020 LOC,
  about half of the ~1900 LOC roadmap estimate.

## v8.0.0-alpha — F1-a Phase 2a: Alginate ionic-Ca L2 solver (2026-04-17)

First non-chitosan-agarose platform lands. Shrinking-core Ca²⁺
diffusion + egg-box gelation gives EmulSim its first ionic-gelation
pipeline. Downstream L3 / L4 callers remain platform-agnostic —
the solver emits a standard `GelationResult`.

### Added

- `PolymerFamily` enum in `datatypes.py`: AGAROSE_CHITOSAN (default)
  / ALGINATE / CELLULOSE / PLGA. Drives future L2 dispatch.
- `MaterialProperties.polymer_family` field + alginate-specific
  defaults (`f_guluronate=0.5`, `D_Ca=1e-9 m²/s`, `k_bind_Ca=1e3
  M⁻²·s⁻¹`, `K_alg_modulus=30 kPa`, `n_alg_modulus=2.0`). Harmless
  for other families.
- `src/emulsim/level2_gelation/ionic_ca.py::solve_ionic_ca_gelation`:
  1D spherical finite-volume BDF solver for C(r,t) / G(r,t) / X(r,t)
  with second-order Ca²⁺ + 2 guluronate → egg-box junction binding.
  ~310 LOC. Returns a SEMI_QUANTITATIVE-tier `GelationResult`.
- `tests/test_alginate_ionic_ca.py` — 13 tests covering the
  PolymerFamily enum, guluronate-concentration helper, result
  schema, guluronate mass conservation (ε < 5 %), zero-Ca /
  zero-alginate edge cases, long-time conversion > 30 % at
  500 mM bath, and input validation.

### Deferred to F1-a Phase 2b

- L4 alginate modulus (`G_DN ∝ (c·f_G)² · X_mean / X_max`).
- Reagent library entries (CaCl₂, internal gelation with GDL + CaCO₃).
- Pipeline orchestrator dispatch by `polymer_family`.
- Config TOML parser support for `polymer_family = "alginate"`.
- Remaining tests from the protocol (√t scaling, modulus scaling,
  full-pipeline integration) — 6 tests deferred.
- Replace c_agarose-slot-as-alginate-proxy with a dedicated
  `FormulationParameters.c_alginate` field.

### Tests

- 77 targeted regression tests pass (F1-a + F3 + F4 + Node 30/31 +
  CLI) in 67 s; 0 regressions.

### Footprint

- **Added:** ~640 LOC (solver + 13 tests + datatypes edits). Matches
  the Phase 2a slice of the ~1900 LOC total projected in
  `docs/f1a_alginate_protocol.md`.

## v8.0.0-alpha — F3-b/c + F4-a: engine wiring + CLI + robust BO (2026-04-17)

Completes v8.0-alpha inverse-design surface. The Node F3-a objective
builders now have an engine-level accessor, a CLI, and a first robust
acquisition stacking mean-variance on top.

### Added

- **F3-b**: `OptimizationEngine(target_spec=...)` constructor param.
  When set, the engine uses `compute_inverse_design_objectives` and
  sizes its internal `REF_POINT` + failure-penalty arrays to
  `len(target_spec.active_dims())`. The 3-objective legacy mode is
  preserved as the default.
- **F3-c**: `python -m emulsim design --d32 ... --pore ... --G-DN ...
  --Kav ...` CLI subcommand with matching `--*-tol` flags.
  TargetSpec.validate() errors route to SystemExit with a clear
  message.
- **F4-a**: `--robust-variance-weight λ` flag + engine kwarg.
  Evaluates λ resamples per candidate and reports
  `mean(obj) + λ · std(obj)` per dimension. Requires `target_spec`
  at construction (robust BO is defined against user targets).
  Current resample strategy: ±1 %·k RPM jitter as a proxy — proper
  spec-driven MC resampling lands in F4-b.
- **F4-b CVaR** — protocol stub only: swap the mean+std layer for a
  CVaR quantile over resamples. Deferred to follow-up (trivial
  change to the same engine path once the resample strategy is
  finalised).

### Tests

- `tests/test_inverse_design_engine.py` — 9 tests covering constructor
  guards, `_n_obj` sizing, robust-BO configuration validation, CLI
  parser registration, and a mocked dispatch path.
- 76 tests pass across F3 + Node 30/31/30b + CLI surfaces; 0
  regressions.

### Deferred

- **F1-a alginate platform**: protocol-only this session at
  `docs/f1a_alginate_protocol.md`. ~1900 LOC projected across L2
  ionic-Ca solver, L4 alginate modulus, PolymerFamily dispatch,
  defaults, and 11 tests. Requires /scientific-advisor briefing at
  kickoff. 3-5 fresh sessions.

## v8.0.0-alpha — Node F3-a: Inverse-design TargetSpec objectives (2026-04-17)

First v8.0 node. Adds user-specified target matching to the
optimisation pipeline so BO can be run in "inverse design" mode
(given target specs, find optimal formulation) rather than the fixed
rotor-stator / stirred-vessel targets only.

### Added

- `TargetSpec` dataclass in `src/emulsim/optimization/objectives.py`:
  per-dimension target + tolerance pairs for d32 (or d_mode in
  stirred-vessel), pore size, G_DN (log10-distance), and Kav (M3
  distribution coefficient, optional). Dimensions are skipped when
  either the target or the tolerance is `None`, so users can target
  subsets. `TargetSpec.validate()` raises on empty spec or
  non-positive tolerance.
- `compute_inverse_design_objectives(result, target, trust_aware=True,
  mode=None)`: returns an objective vector sized to the active
  dimensions. Each component is the tolerance-normalised absolute
  distance (log10 for G_DN); trust penalty from Node 6 is added
  per-component when `trust_aware=True` so weak-evidence candidates
  still land above the engine REF_POINT.
- 12 unit tests in `tests/test_inverse_design_objectives.py` covering
  validate(), active_dims(), per-dimension distance math, trust
  penalty integration, Kav-missing fallback (inf), and stirred-vessel
  d_mode substitution.

### Not yet done (F3-b, F3-c)

- `OptimizationEngine.run(target_spec=...)` integration — the engine
  currently hard-wires `compute_objectives_trust_aware`. Switching the
  objective at runtime requires an engine-level accessor.
- CLI `python -m emulsim design --d32 2e-6 --pore 80e-9 ...` — wraps
  the above into a user surface.

### Footprint

- **Added:** ~200 LOC (TargetSpec + compute_inverse_design_objectives
  + tests). 67 tests across F3-a + UQ + EDC/NHS + panel + CLI pass;
  0 regressions.

## v7.1.0-dev — Node 32: Cluster F v8.0 roadmap (2026-04-17)

Architect-produced roadmap document (no code) at
`docs/node32_cluster_f_v8_roadmap.md`. Refines Doc 10 §4 into Node-level
deliverables for v8.0:

- **F1** Other microsphere platforms (alginate / cellulose NIPS /
  PLGA; alginate recommended first per smallest code delta)
- **F2** Digital twin (EnKF + online Bayesian + MPC; scoped to
  replay-only for v8.0 unless hardware partner emerges)
- **F3** Inverse design (constrained BO; leverages Node 30 UQ +
  Node 6 trust-aware evidence)
- **F4** Robust optimisation under uncertainty (mean-variance /
  CVaR acquisition stacked on F3)
- **F5** MD parameter estimation (MARTINI CG MD for χ, κ, M₀,
  f_bridge; ingest-only default scope)

Proposed v8.0 phasing:

1. Phase 1 (4 weeks): F3 + F4 — inverse design + robust optimisation
   on current platform; v8.0-alpha release.
2. Phase 2 (6 weeks): F1-a alginate; v8.0-beta.
3. Phase 3 (6 weeks): F5 ingest + F2 replay harness; v8.0 GA.

Hard entry criteria: v7.0 must ship (Study A wet-lab gate); CEO /
chief-economist / ip-auditor sign-off on commercial prioritisation
before first v8.0 node.

## v7.1.0-dev — Node 30b: Streamlit UQ panel migration (2026-04-17)

Closes the Node 30 deferral: the streamlit uncertainty panel now
builds a full `UnifiedUncertaintySpec` from user inputs instead of
showing an `st.info` placeholder. The built-in MaterialProperties
perturbations from `UnifiedUncertaintyEngine.run_m1l4` remain always-on;
the panel configures the *additional* spec-driven surface plus sampling
controls and surfaces a count of calibration-posterior sources that
will be absorbed from `st.session_state["_cal_store"]` at engine
construction time.

### Added

- `build_uncertainty_spec(n_samples, seed, custom_sources) ->
  UnifiedUncertaintySpec` — pure helper used by the panel. Invalid
  custom entries (blank name, std <= 0) are silently dropped;
  `n_samples < 1` raises.
- `count_store_posteriors(store) -> int` — tallies calibration entries
  with `posterior_uncertainty > 0` for the panel's status display.
- `CustomSourceInput` dataclass — typed bridge between streamlit
  widget state and the pure spec-builder.
- `tests/test_uncertainty_panel.py` — 12 unit tests of the spec
  builder, the posterior counter, and a panel-export smoke test.

### Changed

- `src/emulsim/visualization/panels/uncertainty.py` — rebuilt around
  the new helpers. UI exposes `n_samples`, `seed`, `n_jobs` (`1`, `2`,
  `4`, `-1`) in a three-column row plus an "Advanced" expander that
  lets the user add up to 10 custom `UncertaintySource` entries
  (name, kind, distribution `normal`/`lognormal`, value, std).
- Session-state surface: the panel persists the built spec at
  `st.session_state["_unc_spec"]` and the parallel-workers value at
  `st.session_state["_unc_n_jobs"]` for downstream run triggers.

### Tests

- 12 new panel tests pass.
- 146 tests across panel + UQ + EDC/NHS + CLI + UI contract surfaces
  pass; 0 regressions.

### Footprint

- **Added:** ~220 LOC (panel rewrite + tests). Net: −28 LOC was the
  placeholder; the migrated panel is ~190 LOC.

## v7.1.0-dev — Node 31: EDC/NHS mechanistic kinetic (2026-04-17)

Promotes EDC/NHS carbodiimide chemistry from QUALITATIVE_TREND
(Node 9 F9 fallback) to SEMI_QUANTITATIVE with a literature-grounded
two-step ODE model. The Hermanson 2013 / Wang 2011 / Cline & Hanna
1988 rate constants close the scientific debt item; Study A calibration
data can promote to QUANTITATIVE via the CalibrationStore posterior
machinery shipped in Node 30.

### Added

- `src/emulsim/module2_functionalization/edc_nhs_kinetics.py` — new
  mechanistic solver. Core: `react_edc_nhs_two_step(...)` integrates
  four ODEs (C → A → E → P) with competing O-acylisourea and
  NHS-ester hydrolyses, returning a structured `EdcNhsResult` with
  `p_final`, `p_hydrolysed`, `p_residual_nhs_ester`, `time_to_half`,
  mass-balance diagnostic, and solver diagnostics. `EdcNhsKinetics`
  dataclass carries the rate constants + activation energies; defaults
  are literature medians at T_ref=298 K. `available_amine_fraction(pH,
  pKa)` helper for chitosan amine speciation.
- `FormulationParameters.pH` field (default 7.0).
- `MaterialProperties.surface_cooh_concentration` field (default 0.0)
  — gates L3 EDC/NHS to run the mechanistic path when non-zero.
- `tests/test_edc_nhs_kinetics.py` — 18 tests covering mass
  conservation, edge cases, Arrhenius / pH / dose-response trends,
  input validation, and M2 + L3 integration.

### Changed

- `src/emulsim/module2_functionalization/modification_steps.py` —
  `_solve_activation_step` dispatches to the mechanistic ODE when
  `reagent_profile.chemistry_class == "edc_nhs"` (was generic
  single-step `solve_second_order_consumption`).
- `src/emulsim/module2_functionalization/reagent_profiles.py` — the
  `edc_nhs_activation` profile's `confidence_tier` promotes from
  `ranking_only` to `semi_quantitative`; `calibration_source` and
  `notes` updated to reference the mechanistic model.
- `src/emulsim/level3_crosslinking/solver.py` — the `michaelis_menten`
  branch now gates on `props.surface_cooh_concentration`. Native matrix
  (= 0) still falls back with QUALITATIVE_TREND (v7.0.1 behaviour
  preserved for safety); carboxylated matrix (> 0) runs the mechanistic
  ODE and ships SEMI_QUANTITATIVE.

### Kept deferred

- Dedicated `c_edc` / `c_nhs` concentration fields on
  `FormulationParameters` (Node 31 reuses `c_genipin`; cleanup = Node
  31b).
- pH-dependent kinetic constants beyond the k_h2 pH term (Node 31b).
- Study A calibration uptake for EDC/NHS-specific matrix chemistry
  (will arrive as `CalibrationStore` entries targeting M2 / L3).

### Tests

- 18 new EDC/NHS tests pass in <1 s.
- 157 tests across the EDC/NHS + UQ + CLI + M2 + L3 + batch +
  run-context surfaces pass in ≈38 s; 0 regressions.

### Footprint

- **Added:** ~420 LOC (solver + tests + L3 gate + M2 dispatch branch).
- **Touched:** 6 files.

## v7.1.0-dev — Node 30: Full UQ merge (2026-04-17)

Consolidates the two legacy Monte Carlo engines into a single
`UnifiedUncertaintyEngine` implementation and closes the Audit N2
calibration-posterior sampling gap left open by Node 18.

### Merged

- `uncertainty_core.py` (318 LOC) — deleted. The M1-L4
  `UncertaintyPropagator` logic lives in `uncertainty_unified.py`.
- `uncertainty_propagation/` package (216 LOC) — deleted. The M2-only
  `M1UncertaintyContract` / `run_with_uncertainty` path had no CLI
  surface and was unreachable outside the streamlit panel + two
  v6.0-era integration tests. An M2-specific UQ path can be rebuilt on
  top of the unified schema in v7.2 if user demand warrants it.
- `UnifiedUncertaintyEngine.run_m2_q_max`, `from_m1_contract_uq`,
  `from_m1l4_result` — deleted (dead adapters).

### Closed

- **Audit N2** (HIGH): `CalibrationStore` posteriors with
  `posterior_uncertainty > 0` now actually perturb the MC on each
  sample. The posterior draw is dispatched by `target_module` — L1
  posteriors land on `params.emulsification.kernels` (lazily
  instantiated if `None`), L2-L4/M2-M3 posteriors land on
  `MaterialProperties`. `result.kinds_sampled` honestly records
  `CALIBRATION_POSTERIOR` when a posterior dispatched to a real
  attribute; `result.kinds_declared_but_not_sampled` only contains
  `CALIBRATION_POSTERIOR` when EVERY posterior failed the dispatch
  (malformed name or unknown attribute).

### CLI

- `python -m emulsim uncertainty --engine {unified,legacy}` now routes
  both choices through the merged engine. `unified` includes
  posteriors; `legacy` runs with `calibration_store=None` for
  byte-compat with v7.0.x scripts that expected only the default
  MaterialProperties perturbations. The output schema is the unified
  summary in both cases — scripts parsing the legacy
  "Uncertainty-Quantified Results" header must migrate.

### Byte-compat

- `_generate_default_perturbations` preserves the exact RNG call order
  of v7.0.1 `UncertaintyPropagator._generate_perturbations`, so
  seed-identical output matches when no posterior sources are declared.
  Posterior draws come AFTER the default 10 draws per sample to avoid
  perturbing the default-only sequence.

### Deferred to Node 30b

- Streamlit UQ panel
  (`src/emulsim/visualization/panels/uncertainty.py`) now displays an
  info placeholder and returns `None`. Full migration to build a
  `UnifiedUncertaintySpec` from streamlit inputs is Node 30b. The CLI
  and programmatic API paths are fully functional in v7.1.

### Tests

- `tests/test_uncertainty_unified.py`: 14 tests. Rewrote the Node 23
  `test_n2_no_posterior_overclaim` as `test_posterior_now_actually_sampled`
  and `test_posterior_actually_perturbs_output` to verify the closure.
  Added `test_unknown_posterior_attribute_skipped` and
  `test_legacy_modules_are_gone`.
- `tests/test_parallel_mc.py`: 4 tests retargeted at the merged
  engine; parallel/serial bit-identicality invariant preserved via the
  new `OutputUncertainty.raw_samples` field.
- `tests/test_cli_v7.py::test_legacy_engine_byte_compat` → rewritten
  as `test_legacy_engine_flag_routes_through_unified`.
- `tests/test_v60_integration.py::TestUncertaintyIntegration` deleted
  (M2 MC path removed).

### Footprint

- **Deleted:** ~540 LOC (uncertainty_core.py + uncertainty_propagation/
  + dead adapters in uncertainty_unified.py).
- **Added:** ~200 LOC (merged sampler + posterior dispatch in
  uncertainty_unified.py, 5 new UQ tests).
- **Net:** ~−340 LOC. 52 tests passing in the targeted regression set
  (UQ + parallel + CLI + v6.0 integration + batch + run-context).

## v7.0.1 (2026-04-17) — Audit remediation patch

Closes 8 of 10 findings from the post-Nodes-1-20 full-system audit. P0
ship-blockers fixed; v7.0 features now reachable from the CLI.

### P0 fixes (release blockers)
- **N1 (HIGH)** — `pipeline/orchestrator.py` no longer mutates the caller's
  `params.emulsification.kernels` in place when applying L1 calibration.
  Callers that reuse a `SimulationParameters` instance across multiple
  `run_single` calls (e.g. `batch_variability.run_batch`, parameter
  sweeps, optimisation campaigns) no longer see calibrated kernels leak
  between iterations. Regression test in `test_run_context.py`.
- **N2 (HIGH)** — `UnifiedUncertaintyEngine.run_m1l4` no longer claims to
  have sampled `CALIBRATION_POSTERIOR` when it has only absorbed the
  posterior into the spec. The new
  `UnifiedUncertaintyResult.kinds_declared_but_not_sampled` field
  records the v7.0 limitation honestly.

### P1 fixes (CLI surface — closes audit N4 + N5)
- **`python -m emulsim batch`** — surface
  `pipeline.batch_variability.run_batch` on the CLI. Pass `--quantiles`
  and `--output`; prints mass-weighted mean / per-quantile percentile
  table.
- **`python -m emulsim dossier`** — run the pipeline and emit a
  `ProcessDossier` JSON artifact for reproducibility. Records inputs,
  result summary, manifests, calibrations, environment.
- **`python -m emulsim ingest L1`** — ingest a directory of
  `AssayRecord` JSON files, run the L1 fitter, write a
  `CalibrationStore`-loadable fit JSON. v7.1 will add L2/L3/L4/M2.
- **`python -m emulsim uncertainty`** now defaults to the
  `UnifiedUncertaintyEngine` (Node 18) and exposes `--n-jobs` for
  Node 15's parallel MC. Pass `--engine legacy` for v6.x byte-equivalent
  output.
- **N3 follow-up** — `QuantileRun.representative_diameter_m` property
  added so downstream consumers don't accidentally read
  `full_result.emulsification.d50` (which is shared by reference across
  all per-quantile runs and reflects the BASE L1 DSD).

### P2 polish
- **N7** — `UncertaintyPropagator.run` auto-falls-back to serial when
  `n_samples < 4 × |n_jobs|`. Joblib startup + Numba JIT cold-compile
  dominate below this threshold.
- **N8** — `run_batch` silently sort+dedupes the `quantiles` argument.
  Duplicate or unsorted input no longer produces ill-defined mass
  fractions.

### P3 documentation
- **N6** — `INSTALL.md` documents the Numba JIT cache location and the
  `NUMBA_CACHE_DIR` environment-variable workaround for read-only
  Python installs (corporate, conda `--no-write-pkgs`,
  `pip install --user` on network shares).
- **N9** — Documenting that Node 8's L2 timing wiring was a metadata
  fix only; the empirical pore-size formula remains independent of
  `alpha_final`. The `model_manifest.diagnostics.alpha_final_from_timing`
  field now reflects the actual Avrami output instead of a hardcoded
  0.999, but pore predictions at typical conditions are unchanged.

### Tests
- 25 new tests across the patch (Nodes 22-29). 0 regressions.

---

## v7.0 (2026-04-17) — Engineering portion (Nodes 14-20)

Closes engineering items from the consensus v7.0 plan (doc 34 §9). F1
closure (kernel re-fit) remains gated on Study A wet-lab data.

### New modules
- `process_dossier.py` — `ProcessDossier` aggregator + JSON export
- `assay_record.py` — `AssayRecord` public data model with 12 `AssayKind` values
- `uncertainty_unified.py` — `UnifiedUncertaintyEngine` single entrypoint
- `pipeline/batch_variability.py` — `run_batch` over DSD quantiles
- `calibration/fitters.py` — stub L1 DSD fitter

### Performance
- Numba JIT for `breakage_rate_alopaeus`, `breakage_rate_coulaloglou`,
  `coalescence_rate_ct` matrix builder (5-10× on coalescence; matches
  NumPy to 1e-12 rtol).
- joblib parallel MC via `UncertaintyPropagator(n_jobs=-1)`.

### Calibration data scaffold
- `data/validation/{l1_dsd,l2_pore,l3_kinetics,l4_mechanics,m2_capacity}/`
  directory tree with JSON-Schema for L1 DSD assays.

---

## v6.0 (2026-04-12) — Calibration-Enabled Process Simulation

Transitions EmulSim from semi-quantitative chemistry simulator to calibration-enabled process simulation platform. All uncalibrated outputs remain semi-quantitative; calibrated outputs reflect user-supplied measurements.

### UI Restructure
- Split monolithic `app.py` (1480 lines) into modular tab architecture (7 UI files, orchestrator < 210 lines)
- `tabs/tab_m1.py`: M1 Fabrication tab (inputs, run, results, optimization, trust)
- `tabs/tab_m2.py`: M2 Functionalization tab (9 step types, 52 reagent profiles)
- `tabs/tab_m3.py`: M3 Performance tab (chromatography + catalysis)
- Sidebar panels for calibration, uncertainty, and lifetime frameworks

### Gradient-Aware LRM (H6)
- `solve_lrm()` accepts time-varying `ProcessState` via `gradient_program` + `equilibrium_adapter`
- Gradient values now mechanistically affect equilibrium during LRM time integration
- `run_gradient_elution()` auto-creates adapter for gradient-sensitive isotherms
- `gradient_sensitive` + `gradient_field` properties on SMA, HIC, IMAC, ProteinA, CompetitiveAffinity isotherms
- Fully backward compatible: existing callers unchanged

### Calibration Framework (v6.0-alpha)
- `CalibrationEntry` typed dataclass with units, target, validity domain (audit F2)
- `CalibrationStore` with JSON import/export, query, and `apply_to_fmc()` (audit F13)
- UI panel: JSON upload, manual entry, color-coded confidence display

### Uncertainty Propagation (v6.0-alpha)
- `M1UncertaintyContract` with 5 CVs and two tiers: measured (Tier 1) vs assumed (Tier 2, audit F4)
- `run_with_uncertainty()` Monte Carlo through M2 pipeline producing p5/p95 bounds on q_max
- UI panel: CV sliders, tier selection, sample count configuration

### Lifetime Projection (v6.0-rc)
- `LifetimeProjection` empirical first-order deactivation model (audit F6)
- `project_lifetime()` with cycles-to-80%/50% milestones
- UI panel: interactive Plotly decay curve, empirical confidence warning

### ProcessState (v6.0-beta)
- Typed `ProcessState` dataclass replacing loose dict for process conditions
- Carries salt, pH, imidazole, sugar competitor, temperature for multi-parameter isotherms
- `EquilibriumAdapter` dispatches by isotherm class name with ProcessState routing

### New Isotherms
- `HICIsotherm`: Salt-modulated Langmuir (K_eff = K_0 * exp(m * C_salt)), requires user calibration
- `CompetitiveAffinityIsotherm`: Generalized competitive binding for lectin elution (Con A, WGA)

### Quality
- 14/14 acceptance criteria from audit Section 7 verified passing
- 24 new integration tests (12 gradient LRM + 12 v6.0 end-to-end)
- All existing v5.9 workflows pass regression (280+ total tests, 0 failures)

---

## v0.1.0 (2026-03-26) — Initial Release

### Simulation Pipeline
- 4-level sequential pipeline: PBE emulsification → empirical gelation → multi-mechanism crosslinking → IPN mechanical properties
- 8 crosslinkers with 4 kinetics models (second-order amine/hydroxyl, UV dose, ionic instant)
- 6 surfactants with Szyszkowski-Langmuir IFT model
- Empirical pore-size model calibrated to literature (Pernodet 1997, Chen 2017)
- 2D Cahn-Hilliard phase-field solver available as advanced option

### Web UI (Streamlit)
- Interactive parameter input with sliders and dropdowns
- Reagent selection (crosslinker + surfactant) with per-reagent defaults
- Per-constant Literature/Custom toggle with calibration protocol links
- Results dashboard with Plotly charts (size distribution, phase field, kinetics, Hertz, Kav)
- Trust assessment with 10 automated reliability checks
- Optimization assessment with actionable recommendations

### CLI
- `python -m emulsim run` — full pipeline
- `python -m emulsim sweep` — RPM parameter sweep
- `python -m emulsim optimize` — BoTorch Bayesian optimization
- `python -m emulsim uncertainty` — Monte Carlo uncertainty propagation
- `python -m emulsim ui` — launch Streamlit web interface
- `python -m emulsim info` — display parameters and properties

### Documentation
- Scientific advisory report (docs/01)
- Computational architecture (docs/02)
- Scientific review with formula verification (docs/03)
- Calibration wet-lab protocol — 5 studies, 1081 lines (docs/04)
- Literature constants database with sources and DOIs
- Reagent library with 8 crosslinkers and 6 surfactants

### Quality Assurance
- 9 rounds of Codex (OpenAI) adversarial review — 63+ findings, all addressed
- Scientific Advisor review — 4 critical bugs fixed
- Dev-Orchestrator usability review — all priorities implemented
- Input validation, trust gates, uncertainty propagation
- 107+ unit tests
