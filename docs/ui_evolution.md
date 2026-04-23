# UI Evolution — Family-First Architecture

> Consolidated from 4 UI design-iteration documents in the 2026-04-24
> content audit. Source files (`docs/07, 14, 15_ui, 16`) are removed;
> the backend-fix audit findings that drove the v9.0 UI redesign are
> captured below. Current UI implementation lives in
> `src/emulsim/visualization/`. Design system direction is in `DESIGN.md`.

## Current state (v9.2.2)

The v9.0 "Family-First" redesign is the current UI architecture:

- **Polymer family is the first M1 input.** Users select Agarose+Chitosan,
  Alginate, Cellulose-NIPS, or PLGA; all downstream M1 widgets are
  dispatched per family.
- **Hardware Mode lives inside M1 Emulsification** (stirred vessel
  Stirrer-A / Stirrer-B, rotor-stator legacy). Removed from Global
  Settings sidebar.
- **Per-output confidence labels.** Every result carries a
  `RunReport.evidence_tier` badge (VALIDATED_QUANTITATIVE /
  CALIBRATED_LOCAL / SEMI_QUANTITATIVE / QUALITATIVE_TREND / UNSUPPORTED)
  that the UI surfaces next to the numeric output.
- **Trust assessment panel** surfaces 10 automated reliability checks
  (non-monotonic L1 trends, chemistry/regime mismatch, empirical L2 in a
  "full-pipeline" claim, etc.).
- **Manual (PDF) download** in the upper-right corner; reagent detail
  pages reachable only via inline contextual links
  (`ui_links.build_reagent_link`).
- **M2 and M3 tabs** wired to their respective backends (not just M1).
- **Streamlit auto-page navigation** suppressed so users cannot
  accidentally land on uninitialized pages.

## Backend-fix audit findings addressed during the v9.0 UI redesign

The v9 redesign was gated on seven backend fixes. Preserved here as an
anchor for future contributors:

| ID | Finding | Severity | Backend fix shipped |
|---|---|---|---|
| F1 | M2/M3 backends unwired; UI was M1-only | Critical | M2 + M3 tabs now call their respective orchestrator paths |
| F2 | No per-output confidence labels | Critical | `OutputMetadata` schema + `evidence_tier` surfacing |
| F3 | M2 UI implied broad chemistry coverage that code didn't support | Critical | Reagent registry gates UI options by `confidence_tier` |
| F4 | pH/temperature didn't affect M2 rates (k0=0.0 bypassed Arrhenius) | High | Arrhenius fix; pH validity gates |
| F5 | Gradient value discarded; not coupled to competitive Langmuir | High | Gradient program wired through `equilibrium_adapter` |
| F6 | Mass balance >2 % in some tests, unexposed | High | Quality enum + backend-enforced tolerance |
| F7 | M1 could recommend unsafe process settings; DDA hardcoded | High | DDA as explicit user input; safety-gated advice |

## Critical reload pattern (v9 lesson preserved)

The Streamlit app reloads `emulsim.datatypes` on every rerun via
`importlib.reload`, minting a new `ModelEvidenceTier` class each time.
Enum members compare by identity → identity comparisons silently break
after the first rerun. **Always compare `PolymerFamily`,
`ModelEvidenceTier`, and similar enum values by `.value`, never by
`is`.** See `RunReport.compute_min_tier` for the reference pattern; it
has a dedicated regression test
(`tests/test_evidence_tier.py::test_min_tier_survives_reloaded_enum_class`).

## Design direction

`DESIGN.md` at the repo root captures the Scientific Instrument
aesthetic (Geist + Geist Mono + JetBrains Mono typography; slate
neutrals + single teal accent; dark-first; no gradients, no decorative
animations). Read before any visual change.

## Pointers

- Current UI code: `src/emulsim/visualization/` (panels, tabs/m1,
  tabs/m2, tabs/m3)
- Design system: `DESIGN.md`
- Claude-facing UI contract: `CLAUDE.md` §"v9.0 Family-First UI contract"
- Pre-audit snapshot of full design iteration trail: tag `v9.2.2-pre-docs-audit`
