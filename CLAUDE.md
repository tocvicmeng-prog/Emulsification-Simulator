# EmulSim — Claude Code instructions

## Design System
Always read `DESIGN.md` at the repo root before making any visual or UI
decisions. All font choices, colors, spacing, and aesthetic direction are
defined there. Do not deviate without explicit user approval.

In QA mode, flag any code that doesn't match DESIGN.md (e.g., hardcoded
colors outside the palette, fonts not in the approved list, purple
gradients, entrance animations, etc.).

## v9.0 Family-First UI contract

- Polymer family is the first M1 input (`tabs/m1/family_selector.py`). It
  drives all downstream rendering.
- Hardware Mode lives inside M1 Emulsification, NOT in Global Settings.
- Per SA §B matrix: never render crosslinker / cooling-rate / pore-model
  widgets for alginate / cellulose / PLGA. Never render chitosan+agarose %
  fields for any family other than AGAROSE_CHITOSAN.
- Always compare PolymerFamily members by `.value`, never by `is`. The
  Streamlit app reloads `emulsim.datatypes` on every rerun, minting a new
  enum class; identity comparisons silently break after the first rerun.
  See `RunReport.compute_min_tier` for the reference pattern.

## Known project quirks

- Repo path contains `文档` (Chinese). Any `print()` of an absolute path
  on Windows will crash with `UnicodeEncodeError` under cp1252. Scripts
  that print paths MUST do `sys.stdout.reconfigure(encoding="utf-8")` first
  (see `docs/user_manual/build_pdf.py` for the canonical pattern).
