# Module 2 (Functionalization) — Design History

> Consolidated from 10 M2 design-iteration documents in the 2026-04-24
> content audit. Source files (`docs/11, 12, 13, 15_m2, 17, 18, 21, 22,
> 23, 24`) are removed; their load-bearing architectural decisions are
> captured below. Current M2 implementation lives in
> `src/emulsim/module2_functionalization/`. Scientific candidate
> screens remain as standalone scientific references.

## Current state (v9.2.2)

M2 is wired end-to-end and ships with the Family-First UI (v9.0+):

- **Activation:** ECH (epichlorohydrin) → EPOXIDE; DVS (divinyl sulfone)
  → VINYL_SULFONE
- **Ligand coupling:** IEX (DEAE, SP), HIC (Phenyl), IMAC (IDA),
  plus the Priority-1 candidates from `docs/19_ligand_protein_coupling_candidates.md`
  staged for implementation (Q, CM, NTA, Butyl, Glutathione, Heparin)
- **Protein coupling:** Protein A, Protein G, plus Priority-1 candidates
  staged (Protein A/G fusion, Streptavidin, Protein L, Concanavalin A)
- **Spacer arms:** currently modelled as `spacer_length_angstrom`
  multiplier on activity retention (Phase 1 per `docs/20_linker_arm_candidates.md`);
  explicit SPACER_ARM step type is Phase 2
- **Quenching:** ethanolamine, mercaptoethanol
- **Mechanistic EDC/NHS:** shipped per Node 31
  (`module2_functionalization/edc_nhs_kinetics.py`), see Appendix A.6 of
  `01_scientific_advisor_report.md` for the full mechanism and rate-constant
  provenance
- **M2 → M3 bridge:** `FunctionalMediaContract` carries ligand identity,
  density, activity retention into M3

## Reagent library structure

Canonical reagent profiles live in
`src/emulsim/module2_functionalization/reagent_profiles.py`. Each profile
carries:

- `chemistry_class` (amine_covalent | hydroxyl_covalent | ionic_reversible | independent_network | edc_nhs | ...)
- `reagent_identity` (explicit compound identity — e.g. "chloroacetic acid" for CM; "glycidyltrimethylammonium chloride" for Q)
- `kinetics_model` (second_order | uv_dose | michaelis_menten | edc_nhs_two_step)
- `k0`, `E_a`, `f_bridge` (or equivalent per model)
- `activity_retention` with uncertainty (for protein coupling)
- `confidence_tier` (ranking_only | semi_quantitative | quantitative)
- `is_macromolecule` flag
- IMAC-specific: `metal_ion`, `metal_loaded_fraction`

## Audit findings addressed during M2 expansion

Ten audit findings from the v5.7 third-party audits were accepted and
shipped. Preserved here as an anchor for future contributors:

| ID | Original finding | How it shipped |
|---|---|---|
| F1 | Quenching double-counted consumed + blocked sites | Quenching increments `blocked_sites` only, not both |
| F2 | "remaining_sites" vs "remaining activated sites" ambiguous | Added `remaining_activated` property; coupling/quenching operate on activated profiles only |
| F3 | Template 2 (hydrolysis competition) needed pH correction | Added `hydrolyzed_sites` field; pH validity gates |
| F4 | Ligand identities under-specified | Split reagent identity from installed ligand; `chemistry_class` field |
| F5 | Protein coupling too simplified for production claim | Labelled as `ranking_only` tier; `activity_retention` user-adjustable with uncertainty |
| F6 | Unit inconsistency for steric limits | Canonical unit = mol/particle; convert in ODE wrapper |
| F7 | M2 → M3 mapping insufficient | `FunctionalMediaContract` bridging contract added |
| F8 | Validators were UI-only, not enforced in backend | Backend validation added in orchestrator.run() |
| F9 | ODE wrappers needed structured diagnostics | `CouplingResult` dataclass returns solver metadata |
| — | "Production-grade" claim overstated | Renamed to "expanded semi-quantitative candidate library" in user-facing surfaces |

## Staged expansions (recommended but not yet shipped)

Priority-1 additions from the scientific candidate screens (see
`docs/19_ligand_protein_coupling_candidates.md` for parameter tables):

- **Ligands:** Q (strong anion), CM (weak cation), NTA (IMAC), Butyl (HIC, mild), Glutathione (GST-tag), Heparin (affinity + IEX)
- **Proteins:** Protein A/G fusion (broadest IgG coverage), Streptavidin (biotin capture), Protein L (Fab/scFv), Concanavalin A (glycoprotein enrichment)
- **Spacers (Priority 1):** DADPA (13 Å, protein amine), AHA (10 Å, EDC/NHS coupling), DAH (9 Å, simple amine)

When these ship, CHANGELOG tier promotion from `ranking_only` → `semi_quantitative` requires at minimum literature rate-constant citation per profile.

## Wet-lab protocols

44 functionalization protocols for the shipped chemistry set live in
`docs/user_manual/appendix_J_functionalization_protocols.md` and are
bundled into the installer release as a PDF. Any new reagent profile
should ship with a matching appendix-J protocol entry.

## Pointers

- Scientific candidates: `docs/19_ligand_protein_coupling_candidates.md`, `docs/20_linker_arm_candidates.md`
- Mechanistic EDC/NHS derivation: `docs/01_scientific_advisor_report.md` §A.6
- Crosslinker provenance (L3 side, related): `docs/01_scientific_advisor_report.md` §A.7
- Wet-lab protocols: `docs/user_manual/appendix_J_functionalization_protocols.md`
- Current code: `src/emulsim/module2_functionalization/`
- Pre-audit snapshot of full design iteration trail: tag `v9.2.2-pre-docs-audit`
