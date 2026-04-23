# EmulSim Documentation Index

> Navigation map for `docs/`. Maintained post-audit (2026-04-24; branch
> `docs/content-audit-2026-04`). Policy: every `.md` file in this
> directory is listed below. If you add a new doc, add an entry.

## User-facing

| File | Purpose |
|---|---|
| [`quickstart.md`](quickstart.md) | First simulation in five minutes |
| [`configuration.md`](configuration.md) | Every TOML field with units and physical meaning |
| [`04_calibration_protocol.md`](04_calibration_protocol.md) | 5-study wet-lab calibration plan |

## Scientific foundation

| File | Purpose |
|---|---|
| [`01_scientific_advisor_report.md`](01_scientific_advisor_report.md) | First-principles physics/chemistry decomposition (L1–L4). **Canonical scientific reference.** Includes appendices consolidated from older review/audit/brief docs in the 2026-04-24 content audit. |
| [`02_computational_architecture.md`](02_computational_architecture.md) | Software design, data flow, solver implementations |
| [`SA-EMULSIM-XL-001_Crosslinker_Evaluation.md`](SA-EMULSIM-XL-001_Crosslinker_Evaluation.md) | Scientific justification for the 9-crosslinker library (STMP provenance lives here) |

## Architecture Decision Records

| File | Purpose |
|---|---|
| [`decisions/ADR-001-python-version-policy.md`](decisions/ADR-001-python-version-policy.md) | Why Python is pinned to `>=3.11,<3.13` |
| [`decisions/ADR-002-optimization-stack-pin.md`](decisions/ADR-002-optimization-stack-pin.md) | Why botorch / gpytorch / torch are version-pinned |

## Platform / family deep-dives

| File | Purpose |
|---|---|
| [`f1a_alginate_protocol.md`](f1a_alginate_protocol.md) | Alginate ionic-Ca gelation |
| [`f1b_cellulose_nips_protocol.md`](f1b_cellulose_nips_protocol.md) | Cellulose NIPS |
| [`f1c_plga_protocol.md`](f1c_plga_protocol.md) | PLGA solvent evaporation |
| [`f2_digital_twin_protocol.md`](f2_digital_twin_protocol.md) | Digital twin / EnKF replay |
| [`f4b_cvar_protocol.md`](f4b_cvar_protocol.md) | CVaR-robust Bayesian optimization |
| [`f5_md_ingest_protocol.md`](f5_md_ingest_protocol.md) | MARTINI MD parameter ingest |

## User manual

| File | Purpose |
|---|---|
| [`user_manual/polysaccharide_microsphere_simulator_first_edition.md`](user_manual/polysaccharide_microsphere_simulator_first_edition.md) | First Edition — Getting Started, Platform Catalogue, Appendices A–I |
| [`user_manual/appendix_J_functionalization_protocols.md`](user_manual/appendix_J_functionalization_protocols.md) | 44 wet-lab functionalization protocols with SDS-lite safety blocks |

## Module / feature history

These files capture development history preserved for traceability. Current state lives in `src/` and CHANGELOG.

| File | Purpose |
|---|---|
| `module2_history.md` (created in PR-E) | Consolidated M2 expansion history (fold of 12 prior iteration docs) |
| `ui_evolution.md` (created in PR-E) | UI alignment / Family-First evolution |

## Policy

Version-specific planning documents for superseded releases (v5.x/v6.x/v7.x/v8.x) were removed in the 2026-04-24 content audit. The pre-audit snapshot lives at tag `v9.2.2-pre-docs-audit`. User-facing history is in `CHANGELOG.md`.
