# Design System — EmulSim

> Created 2026-04-18 via `/design-consultation` at the end of the v9.0 UI redesign.
> Read this before making any visual or UI decisions. Flag any code that deviates
> in QA mode.

## Product Context
- **What this is:** Multi-scale simulator for double-network polysaccharide microsphere fabrication (Emulsification → Gelation → Crosslinking → Mechanical Properties). Produced by Holocyte Pty Ltd / GMExpression.
- **Who it's for:** Chemistry / materials-science researchers, downstream-processing technicians, junior researchers, IP auditors, chief economists reviewing cost of goods.
- **Space:** Scientific instrument / simulation software. Adjacent category: Igor Pro, Origin, MATLAB, COMSOL, Bloomberg Terminal (for data density expectations).
- **Project type:** Streamlit-based desktop/web scientific instrument. NOT a marketing site. NOT a general-purpose SaaS dashboard.

## Aesthetic Direction
- **Direction:** Scientific Instrument (Industrial / Utilitarian base + editorial typography for documentation)
- **Decoration level:** Minimal. Typography + precise alignment + semantic color do all the work. No gradients, no decorative blobs, no 3-column icon grids.
- **Mood:** Calibrated. Measured. Serious. A lab instrument that respects the user's time and expertise.
- **North stars:** Igor Pro (data density), Linear (craft polish), Vercel dashboard (color restraint).
- **Anti-north-stars:** Anything that looks like a Stripe knockoff. Anything with a purple gradient as the primary accent. Anything that centers everything.

## Typography

Loaded via Google Fonts CDN in `app.py` head injection (see `src/emulsim/visualization/app.py`).

| Role | Font | Weight | Size | Rationale |
|---|---|---|---|---|
| Display / Hero (H1) | Geist Sans | 700 | 28-32px | Modern grotesque. Legible at any weight. Not on the overused list. |
| Section header (H2) | Geist Sans | 600 | 20-24px | |
| Sub-header (H3) | Geist Sans | 600 | 16-18px | |
| Body UI | Geist Sans | 400-500 | 14px | |
| Label / caption | Geist Sans | 500 | 12-13px | |
| Data / tables / metrics | Geist Mono | 400 | 14px | Tabular numerals enabled (`font-feature-settings: "tnum"`). Scientific data must align on the decimal point. |
| Equations / code / SMILES | JetBrains Mono | 400 | 13px | Wide Unicode support for µ, α, β, γ, χ, σ, Ω, arrow glyphs. |
| Manual / long-form docs | Geist Sans | 400 | 15px, 1.6 line-height | Editorial rhythm for the First Edition user manual. |

**Loading strategy:** Google Fonts `<link>` tag in page config. No self-hosting needed; repo stays lean.

**Do not use:** Inter, Roboto, Arial, Helvetica, Source Sans Pro (the Streamlit default), Open Sans, Montserrat, Poppins, Comic Sans, Papyrus.

## Color

Palette is deliberately **restrained**: one accent + cool neutrals. The product's primary chromatic output is the plotly scientific charts (DSD, phase field, breakthrough curves); brand colors must not compete with them.

### Neutrals (Slate)

| Token | Hex | Usage |
|---|---|---|
| `slate-50` | `#F8FAFC` | Light-mode page background |
| `slate-100` | `#F1F5F9` | Light-mode surface (cards, inputs) |
| `slate-200` | `#E2E8F0` | Light-mode border |
| `slate-400` | `#94A3B8` | Muted text, placeholder |
| `slate-600` | `#475569` | Secondary text |
| `slate-700` | `#334155` | Dark-mode surface (cards) |
| `slate-800` | `#1E293B` | Dark-mode surface-2 |
| `slate-900` | `#0F172A` | Dark-mode background |
| `slate-950` | `#020617` | Dark-mode deep background |

### Accent

| Token | Hex | Usage |
|---|---|---|
| `teal-400` | `#2DD4BF` | Dark-mode accent (primary buttons, links, focus ring) |
| `teal-500` | `#14B8A6` | Light-mode accent (primary buttons, links, focus ring) |
| `teal-600` | `#0D9488` | Hover state |

Rationale: teal is distinct from the default plotly palette (red/blue/green/orange), non-gendered, non-corporate. Reads as "measurement / sensor / scope" — appropriate for a simulator.

### Semantic (evidence tiers + trust gates)

These must map to `ModelEvidenceTier` enum in `datatypes.py`.

| Tier / state | Hex | Usage |
|---|---|---|
| `VALIDATED_QUANTITATIVE` | `#16A34A` (green-600) | Evidence badge "VALIDATED" |
| `CALIBRATED_LOCAL` | `#22C55E` (green-500) | Evidence badge "CALIBRATED" |
| `SEMI_QUANTITATIVE` | `#F59E0B` (amber-500) | Evidence badge "SEMI-QUANTITATIVE" |
| `QUALITATIVE_TREND` | `#F97316` (orange-500) | Evidence badge "QUALITATIVE TREND" |
| `UNSUPPORTED` | `#DC2626` (red-600) | Evidence badge "UNSUPPORTED" |
| `info` | `#0284C7` (sky-600) | Neutral informational banners |
| `success` | `#16A34A` | Success callouts (run completed) |
| `warning` | `#F59E0B` | Trust warnings (reagent outside calibration range) |
| `error` | `#DC2626` | Blockers (invalid inputs, simulation failed) |

### Dark mode

Dark mode is **primary**. Labs run dim, screens run dim, dark themes fatigue less. Light mode is available for print/publication screenshots.

- Background: `slate-900` (`#0F172A`)
- Surface: `slate-800` (`#1E293B`)
- Border: `slate-700` with 50% alpha
- Primary text: `slate-50`
- Muted text: `slate-400`
- Accent: `teal-400` (lower saturation for dark bg)
- Saturation reduced 10-15% on semantic colors vs. light mode.

## Spacing

- **Base unit:** 8px. All spacing is a multiple of 4px (half-unit for tight metric alignment).
- **Density:** Compact. Scientific data tooling values density over whitespace. Reject the "generous whitespace" default — that's a marketing-site pattern.

| Token | px | Usage |
|---|---|---|
| `space-0.5` | 2 | Inline separators |
| `space-1` | 4 | Tight padding, badge inner padding |
| `space-2` | 8 | Default gap between inline elements |
| `space-3` | 12 | Widget internal padding |
| `space-4` | 16 | Gap between sibling widgets |
| `space-6` | 24 | Gap between section subheads |
| `space-8` | 32 | Gap between top-level sections |
| `space-12` | 48 | Major structural breaks |

## Layout

- **Approach:** Grid-disciplined. Streamlit columns = fine. Enforce consistent column gutters.
- **Max content width:** fluid. Never constrain; scientists want to see their data table wide.
- **Border radius:** 4px uniform. No pill buttons, no 16px+ soft-corner cards. Tools, not toys.
- **Elevation:** minimal. One layer of elevation (surface on background), plus 1px border. Never drop shadows — noisy for data UIs.

## Motion

- **Approach:** Minimal-functional only.
- **Transitions:** 150ms ease-out on hover / focus state changes. Period.
- **Forbidden:** entrance animations, scroll-driven effects, spring/bouncy curves, page transitions, loading spinners with personality. Every animation is a moment when the user doesn't know if the result is final — in a simulator, that's trust leak.

## Component conventions

- **Buttons (primary):** teal-500 bg, white text, 4px radius, 8px × 16px padding. Hover: teal-600.
- **Buttons (secondary):** slate-100 bg (dark: slate-700), slate-900 text (dark: slate-50).
- **Inputs:** 1px slate-200 border (dark: slate-700), 4px radius, 8px × 12px padding. Focus: 2px teal-500 ring.
- **Evidence badges:** 2px × 8px padding, 2px radius, semantic-color bg with 10% alpha + same-color text at 100%. No uppercase.
- **Tables:** header row bold + slate-500 text + slate-50 bg (dark: slate-800). Zebra striping disabled. Tabular-nums on numeric cells.
- **Metric cards:** large number in Geist Mono 700 at 28px, label in Geist Sans 500 at 12px, delta chip below. Evidence badge inline with label.

## Decisions Log

| Date | Decision | Rationale |
|---|---|---|
| 2026-04-18 | Adopt Geist Sans + Mono over Inter / Source Sans Pro | Differentiates from every other scientific-SaaS UI; tabular nums by default in the mono; one family = coherence. |
| 2026-04-18 | Dark mode primary, not a toggle-afterthought | Target users work in labs with dimmed ambient light for long sessions. |
| 2026-04-18 | Teal as sole accent; neutralise everything else | The scientific plots are the colour performance. Brand colour must not compete. |
| 2026-04-18 | Reject "generous whitespace"; 8px base with compact density | Data tooling competes on information per viewport, not on breathing room. |
| 2026-04-18 | No entrance animations, ever | Every unfinished animation = "is the simulation still running?" uncertainty. |
