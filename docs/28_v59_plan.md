# EmulSim Module 2 --- v5.9 Comprehensive Plan

**Date:** 2026-04-12
**Baseline:** v5.8 (42 profiles, 6 step types, 8 ACSSiteType values)
**Scope:** Address ALL deferred items from v5.7, v5.8, and audit reports (docs 23, 24, 26, 27)

---

## 1. Current Architecture Summary

| Component | Count | Key Files |
|-----------|-------|-----------|
| Reagent profiles | 42 | `reagent_profiles.py` (1382 lines) |
| Step types | 6 | `ModificationStepType` enum in `modification_steps.py` |
| ACS site types | 8 | `ACSSiteType` enum in `acs.py` |
| ODE templates | 3 | `reactions.py` (Templates 1-3 + quench) |
| M3 isotherms | 5 | `langmuir`, `competitive_langmuir`, `sma`, `imac`, `protein_a` |
| FMC fields | 18 | `FunctionalMediaContract` in `orchestrator.py` |

Validated workflows: SECONDARY_CROSSLINKING, ACTIVATION, LIGAND_COUPLING, PROTEIN_COUPLING, QUENCHING, SPACER_ARM.

---

## 2. Deferred Item Analysis (17 Items)

### Item D1: METAL_CHARGING Step Type for IMAC
**Source:** v5.7 doc 24, Phase 1.5

**SCIENTIFIC ADVISOR:**
- Ready for implementation. The chemistry is well-characterized: NTA/IDA chelators are loaded with Ni2+, Co2+, Cu2+, or Zn2+ from a metal salt solution (e.g., 50 mM NiSO4). The reaction is a rapid ligand-exchange equilibrium (seconds to minutes). Regeneration uses EDTA or imidazole stripping followed by re-charging.
- Data needed: Metal-chelator association constants (log K_NTA-Ni = 11.5, log K_IDA-Ni = 8.1); stripping efficiency with 50 mM EDTA (~99%).
- **Priority: MUST-HAVE for v5.9.** IMAC is the single most common His-tag purification chemistry and currently has a placeholder `metal_loaded_fraction=1.0`.

**ARCHITECT:**
- New `ModificationStepType.METAL_CHARGING` enum value.
- New solver `_solve_metal_charging_step()` in `modification_steps.py`: simple equilibrium binding model. No ODE needed --- metal loading is fast and reaches equilibrium.
- New fields on `ReagentProfile`: `metal_association_constant`, `stripping_agent`, `stripping_efficiency`.
- New reagent profiles: `nickel_charging`, `cobalt_charging`, `copper_charging`, `zinc_charging`, `edta_stripping`.
- Update `_STEP_ALLOWED_RTYPES` and `_STEP_ALLOWED_REACTION_TYPES` with `"metal_charging"` and `"metal_stripping"`.
- Update FMC builder to read actual `metal_loaded_fraction` from metal charging step result.
- LOC estimate: ~250

**DEV-ORCHESTRATOR:**
- Work node: **WN-D1**
- Model tier: Sonnet (straightforward equilibrium calc)
- Dependencies: None. Can start immediately.

---

### Item D2: binding_model_hint Routing in M3
**Source:** v5.7 doc 24, Phase 1.5

**SCIENTIFIC ADVISOR:**
- Ready for implementation. The `binding_model_hint` field is already propagated through FMC (line 109 of `orchestrator.py`) but M3's orchestrator (`module3_performance/orchestrator.py`) ignores it. M3 already has the right isotherm implementations (SMA for IEX, IMAC competitive Langmuir, Langmuir for affinity).
- The routing logic is a pure software concern --- no new science needed.
- **Priority: MUST-HAVE for v5.9.** Without routing, M3 defaults to Langmuir for everything, which is incorrect for IEX (should use SMA) and IMAC (should use competitive Langmuir with imidazole).

**ARCHITECT:**
- Add `_select_isotherm(fmc: FunctionalMediaContract)` dispatch function in M3 orchestrator.
- Routing map: `charge_exchange` -> SMA, `metal_chelation` -> IMACCompetitionIsotherm, `salt_promoted` -> Langmuir (HIC placeholder), `fc_affinity`/`gst_glutathione`/`near_irreversible`/`kappa_light_chain_affinity` -> Langmuir (affinity), `mixed_mode` -> SMA+Langmuir (warn: simplified).
- Auto-populate isotherm parameters from FMC: `q_max` from `estimated_q_max`, `Lambda` from `charge_density` for SMA.
- LOC estimate: ~180

**DEV-ORCHESTRATOR:**
- Work node: **WN-D2**
- Model tier: Sonnet
- Dependencies: None. Can run in parallel with WN-D1.

---

### Item D3: accessible_area_per_bed_volume on FMC
**Source:** v5.7 doc 24, Phase 1.5

**SCIENTIFIC ADVISOR:**
- Ready. This is a geometric calculation: `a_v_accessible = reagent_accessible_area * (6*(1-eps)/d_p) / A_external`. Currently `a_v` is computed inline in the FMC builder (line 193) using external geometric area only. The accessible version accounts for pore penetration.
- No new data needed; all inputs exist in `AccessibleSurfaceModel`.
- **Priority: MUST-HAVE for v5.9.** Corrects q_max estimates for porous beads (currently underestimates by ~2-5x for macroporous media).

**ARCHITECT:**
- Add `accessible_area_per_bed_volume: float` field to `FunctionalMediaContract`.
- Compute in `build_functional_media_contract()` using `surface_model.reagent_accessible_area`.
- Replace inline `a_v` calculation in q_max mapping with this field.
- LOC estimate: ~40

**DEV-ORCHESTRATOR:**
- Work node: **WN-D3**
- Model tier: Haiku (trivial field addition)
- Dependencies: None.

---

### Item D4: Orientation/Activity Uncertainty Metadata for Proteins
**Source:** v5.7 doc 24, Phase 1.5

**SCIENTIFIC ADVISOR:**
- Ready. `activity_retention_uncertainty` already exists on `ReagentProfile` (populated for all protein profiles, typically 0.10-0.15). Need to propagate this through `ModificationResult` and `FunctionalMediaContract` so M3 can produce confidence intervals.
- **Priority: Nice-to-have for v5.9.** The data exists; this is plumbing.

**ARCHITECT:**
- Add `activity_retention_uncertainty: float` to `FunctionalMediaContract`.
- Add `functional_density_uncertainty: float` to `FunctionalMediaContract`.
- Propagate from last protein coupling step's reagent profile through FMC builder.
- Add `q_max_lower` / `q_max_upper` bound fields to FMC.
- LOC estimate: ~60

**DEV-ORCHESTRATOR:**
- Work node: **WN-D4**
- Model tier: Haiku
- Dependencies: None.

---

### Item D5: EDC/NHS Chemistry for AHA -COOH Terminus
**Source:** v5.7 doc 24, Phase 2

**SCIENTIFIC ADVISOR:**
- Ready for implementation with caveats. EDC/NHS is the standard carboxyl activation method. The two-step mechanism is: (1) EDC activates -COOH to O-acylisourea intermediate (fast, k ~0.1 M^-1s^-1 at pH 5); (2) NHS converts to NHS-ester (stabilization step); (3) NHS-ester couples to primary amine (k ~1e-3 M^-1s^-1 at pH 7.5).
- Data needed: EDC hydrolysis rate (~0.003 /s at pH 5), NHS-ester hydrolysis rate (~1e-4 /s at pH 7.5). Well-characterized in Hermanson (2013).
- Requires new `ACSSiteType.NHS_ESTER` and `ACSSiteType.CARBOXYL` (CARBOXYL already exists in enum but has no M1 initialization path).
- **Priority: MUST-HAVE for v5.9.** Completes the AHA spacer pathway (currently noted as "EDC/NHS path not modeled" in aha_spacer notes).

**ARCHITECT:**
- New step type: `ModificationStepType.CARBOXYL_ACTIVATION` (or reuse ACTIVATION with CARBOXYL target).
- New ACSSiteType: `NHS_ESTER`.
- New reagent profiles: `edc_nhs_activation` (reaction_type="activation", target=CARBOXYL, product=NHS_ESTER), `edc_activation` (no NHS stabilization variant).
- New solver: Reuse Template 2 (competitive hydrolysis) for NHS-ester formation with EDC hydrolysis as competing pathway.
- Update AHA spacer notes to reference the new profiles.
- LOC estimate: ~200

**DEV-ORCHESTRATOR:**
- Work node: **WN-D5**
- Model tier: Opus (new chemistry pathway with multi-step mechanism)
- Dependencies: None for core implementation. Integrates with existing ACTIVATION workflow.

---

### Item D6: pH-Dependent Rate Scaling (Global M2 Enhancement)
**Source:** v5.7 doc 24, Phase 2

**SCIENTIFIC ADVISOR:**
- Ready with simplifications. Most M2 reactions have pH-dependent rate constants. Current implementation passes `ph` to coupling solvers but only uses it for validity window checks (ph_min/ph_max). A sigmoid scaling factor is the pragmatic approach:
  `k_eff = k_ref * sigmoid_ph_factor(ph, pKa, n_Hill)`
  where pKa and n_Hill are per-reagent parameters.
- Data needed: pKa values for nucleophilic groups (e.g., amine pKa ~10.5, thiol pKa ~8.3, hydroxyl pKa ~14). Hill coefficients can default to 1.
- **Priority: Nice-to-have for v5.9, MUST-HAVE for v6.0.** Improves quantitative accuracy but ranking is already correct with current validity windows.

**ARCHITECT:**
- Add `pKa_nucleophile: float` and `ph_hill_coefficient: float` fields to `ReagentProfile` (defaults: 0.0 = disabled).
- Add `_ph_rate_scaling(ph, pKa, n)` helper in `reactions.py`.
- Apply scaling in all four ODE templates when pKa > 0.
- Update all 42 reagent profiles with appropriate pKa values.
- LOC estimate: ~150 (code) + ~100 (profile updates)

**DEV-ORCHESTRATOR:**
- Work node: **WN-D6**
- Model tier: Sonnet
- Dependencies: None, but should be tested after WN-D5 (NHS hydrolysis is pH-sensitive).

---

### Item D7: Salt-Dependent HIC Isotherms (M3 Extension)
**Source:** v5.7 doc 24, Phase 2

**SCIENTIFIC ADVISOR:**
- NOT ready for implementation as a first-principles model. HIC retention depends on the Hofmeister series, and the relationship between ammonium sulfate concentration and hydrophobic interaction is empirically fitted, not mechanistically predicted from ligand density alone.
- A pragmatic placeholder: Langmuir with salt-modulated K_eq: `K_eff = K_0 * exp(m * C_salt)` where m is the molal surface tension increment. Literature values exist for phenyl/butyl/octyl with ammonium sulfate.
- **Priority: Defer to v6.0.** HIC is functional with Langmuir placeholder. Salt dependence requires user calibration data we cannot generate from M2 outputs alone.

**ARCHITECT (v6.0 sketch):**
- New `HICIsotherm` class in `isotherms/hic.py`.
- Fields: K_0, m_salt, q_max, salt_type.
- Integrate with binding_model_hint="salt_promoted" routing from WN-D2.
- LOC estimate: ~200

**DEV-ORCHESTRATOR:**
- Work node: **WN-D7** (deferred to v6.0)
- Model tier: Opus (new isotherm class with salt thermodynamics)
- Dependencies: WN-D2 (routing must exist first).

---

### Item D8: Reversible/Irreversible Affinity Binding Modes in M3
**Source:** v5.7 doc 24, Phase 2

**SCIENTIFIC ADVISOR:**
- Partially ready. Reversible binding is already modeled (Langmuir is inherently reversible). Irreversible binding (streptavidin-biotin, Kd ~10^-15 M) needs a one-way Langmuir: once bound, no desorption. This is a simple modification: set k_desorption = 0 in the LRM transport model.
- Also needed: "slow off-rate" regime for Protein A (Kd ~10^-8 M) where wash efficiency depends on residence time.
- **Priority: Nice-to-have for v5.9.** The streptavidin profile already has `binding_model_hint="near_irreversible"` but M3 treats it as reversible Langmuir.

**ARCHITECT:**
- Add `is_irreversible: bool` parameter to `LangmuirIsotherm`.
- When irreversible: modify `solve_lrm` to use one-way mass transfer (k_des = 0 after adsorption).
- Route from `binding_model_hint="near_irreversible"` in WN-D2.
- LOC estimate: ~80

**DEV-ORCHESTRATOR:**
- Work node: **WN-D8**
- Model tier: Sonnet
- Dependencies: WN-D2 (routing).

---

### Item D9: TMAE via DVS (Strong Anion Exchanger)
**Source:** v5.8 doc 27

**SCIENTIFIC ADVISOR:**
- NOT ready. TMAE coupling via DVS requires the vinyl sulfone to react with a tertiary amine (trimethylaminoethanol). The kinetics of VS-tertiary amine coupling are poorly characterized compared to VS-primary amine. The reactive precursor (2-(trimethylamino)ethanol or related quaternary amine nucleophile) has under-specified kinetic data.
- Would need: k_forward for VS + TMAE nucleophile, competitive hydrolysis rate of VS under TMAE coupling conditions. No reliable literature values found.
- **Priority: Defer to v6.0.** Requires experimental kinetic data or literature search. Q (quaternary amine via glycidyltrimethylammonium on epoxide) already covers the strong anion exchange use case.

**ARCHITECT (v6.0 sketch):**
- New reagent profile `tmae_vs_coupling` targeting `ACSSiteType.VINYL_SULFONE`.
- Requires extending coupling workflow to accept VS targets (currently only epoxide coupling profiles exist).
- LOC estimate: ~60 (profile) + ~30 (VS target support)

**DEV-ORCHESTRATOR:**
- Work node: **WN-D9** (deferred to v6.0)
- Dependencies: Literature kinetic data acquisition.

---

### Item D10: Lectin-Specific M3 Elution Models
**Source:** v5.8 doc 27

**SCIENTIFIC ADVISOR:**
- NOT ready for first-principles implementation. Lectin elution uses sugar competition (e.g., 0.5 M methyl-alpha-D-mannopyranoside for Con A, 0.5 M GlcNAc for WGA). This is mechanistically identical to IMAC imidazole competition --- competitive Langmuir between protein and sugar for lectin binding sites.
- Data needed: K_mannose-ConA (~1e3 M^-1), K_GlcNAc-WGA (~1e2 M^-1). These exist in lectin biochemistry literature but vary significantly with buffer conditions and Ca2+/Mn2+ cofactor concentration.
- **Priority: Defer to v6.0.** Con A and WGA profiles exist (profiles 30-31) with `m3_support_level="requires_user_calibration"`. The IMAC competitive isotherm template is reusable once sugar competition constants are populated.

**ARCHITECT (v6.0 sketch):**
- Parameterize `IMACCompetitionIsotherm` to accept generic "competitor" species (not just imidazole).
- Rename or create `CompetitiveAffinityIsotherm` base class.
- Add `sugar_competitor_K` field to lectin reagent profiles.
- LOC estimate: ~120

**DEV-ORCHESTRATOR:**
- Work node: **WN-D10** (deferred to v6.0)
- Dependencies: WN-D2 (routing), literature K values.

---

### Item D11: Protein Disulfide Reduction Step (Pre-Treatment)
**Source:** v5.8 doc 27

**SCIENTIFIC ADVISOR:**
- Ready for implementation. TCEP or DTT reduction of disulfide bonds exposes free thiols for maleimide-thiol coupling. TCEP is preferred (no interference with maleimide). The chemistry is well-characterized: TCEP reduces disulfides stoichiometrically at pH 7, room temperature, in minutes.
- This is a pre-treatment that modifies the protein, not the bead surface. It does not consume ACS sites. Instead, it determines the `thiol_accessibility_fraction` of the protein.
- Data needed: TCEP reduction efficiency (~95% for surface-exposed disulfides), reaction time (30 min typical).
- **Priority: MUST-HAVE for v5.9.** The `requires_reduced_thiol=True` flag exists on maleimide-thiol profiles but there is no mechanism to model the reduction step. Currently `thiol_accessibility_fraction` is a static value.

**ARCHITECT:**
- New `ModificationStepType.PROTEIN_PRETREATMENT` enum value.
- New solver `_solve_protein_pretreatment_step()`: simple first-order reduction model.
- New reagent profiles: `tcep_reduction`, `dtt_reduction`.
- The solver modifies a "protein state" object (not ACS state) that downstream PROTEIN_COUPLING reads for effective `thiol_accessibility_fraction`.
- Add `ProteinPretreatmentState` dataclass to track: `reduction_fraction`, `excess_reductant_removed`, `time_since_reduction` (for reoxidation kinetics).
- LOC estimate: ~200

**DEV-ORCHESTRATOR:**
- Work node: **WN-D11**
- Model tier: Sonnet
- Dependencies: None. Should complete before protein coupling integration testing.

---

### Item D12: Metal Loading/Regeneration and Imidazole Competition for IMAC
**Source:** Audit docs 23, 26

**SCIENTIFIC ADVISOR:**
- This overlaps significantly with D1 (metal charging) and D2 (binding_model_hint routing to IMAC isotherm). The imidazole competition part is already implemented in `isotherms/imac.py`.
- Remaining gap: Metal leaching model (Ni2+ release during operation, ~0.1-1% per cycle). This is a long-term stability concern (see D16).
- **Priority: Addressed by WN-D1 + WN-D2.** Leaching deferred to D16.

**ARCHITECT:**
- No separate work node. Covered by WN-D1 and WN-D2.

---

### Item D13: Target-Specific Default Assumptions
**Source:** Audit docs 23, 26

**SCIENTIFIC ADVISOR:**
- Ready for implementation as a data layer (no new physics). Provide default protein parameters for common purification targets:
  - **IgG:** MW 150 kDa, r_h 5.3 nm, pI 6-9, z_eff 3-8 (SMA)
  - **His-tag fusion (25-50 kDa):** 6xHis, K_protein ~1e4 m^3/mol on Ni-NTA
  - **GST-tag fusion (50-80 kDa):** K_GSH ~1e3 m^3/mol on glutathione
  - **Biotin-tag:** K_biotin ~1e10 m^3/mol on streptavidin (near-irreversible)
  - **Heparin-binding proteins:** Variable; AT-III as reference (K ~1e5 m^3/mol)
- **Priority: Nice-to-have for v5.9.** Improves usability but does not affect M2 accuracy.

**ARCHITECT:**
- New `target_protein_defaults.py` in module3_performance with `TARGET_PROTEIN_LIBRARY` dict.
- Each entry: `TargetProtein(name, mw, r_h, pI, z_eff, sigma_sma, K_affinity, tag_type, notes)`.
- M3 orchestrator uses these defaults when user does not specify target parameters.
- LOC estimate: ~150

**DEV-ORCHESTRATOR:**
- Work node: **WN-D13**
- Model tier: Haiku (data entry)
- Dependencies: WN-D2 (so M3 knows which isotherm to parameterize).

---

### Item D14: Batch-to-Batch DDA Variability and Uncertainty Propagation from M1
**Source:** Audit docs 23, 26

**SCIENTIFIC ADVISOR:**
- NOT ready for rigorous implementation. Requires M1 to provide uncertainty distributions on bead_d50, porosity, pore_size_mean, nh2_bulk_concentration, oh_bulk_concentration. Currently M1 provides point estimates only.
- A pragmatic approach: Add user-configurable `cv_bead_d50`, `cv_porosity`, etc. (coefficient of variation) to `M1ExportContract`, then propagate via Monte Carlo or analytical error propagation through M2.
- **Priority: Defer to v6.0.** Requires M1 interface changes and statistical framework.

**ARCHITECT (v6.0 sketch):**
- Add `uncertainty` sub-dataclass to `M1ExportContract`.
- Monte Carlo wrapper around `ModificationOrchestrator.run()` (N=100 samples).
- Report percentile bounds on FMC outputs.
- LOC estimate: ~400

**DEV-ORCHESTRATOR:**
- Work node: **WN-D14** (deferred to v6.0)
- Dependencies: M1 interface update (external).

---

### Item D15: Pore-Size Distribution (Not Just Mean)
**Source:** Audit docs 23, 26

**SCIENTIFIC ADVISOR:**
- NOT ready. Pore-size distribution (PSD) affects accessibility non-linearly: a bimodal PSD with 50% macropores (100 nm) and 50% micropores (5 nm) has very different protein accessibility than a unimodal 52.5 nm mean.
- Would need: PSD representation (log-normal parameters, or discrete bins), integration of accessibility fraction over the distribution.
- Current `AccessibleSurfaceModel` uses mean pore diameter with Renkin exclusion. This is adequate for unimodal distributions with CV < 0.5.
- **Priority: Defer to v6.0.** Requires PSD data from M1 (currently provides only mean).

**ARCHITECT (v6.0 sketch):**
- Add `pore_size_distribution: Optional[tuple[np.ndarray, np.ndarray]]` (diameters, fractions) to `M1ExportContract`.
- Integrate `f_accessible(r_h)` over PSD bins in `AccessibleSurfaceModel`.
- LOC estimate: ~200

**DEV-ORCHESTRATOR:**
- Work node: **WN-D15** (deferred to v6.0)
- Dependencies: M1 PSD data output.

---

### Item D16: Protein Leaching/Deactivation Kinetics
**Source:** Audit docs 23, 26

**SCIENTIFIC ADVISOR:**
- NOT ready for first-principles modeling. Protein leaching from affinity resins is governed by: (1) covalent bond stability (thioether > amide > ester), (2) alkaline CIP degradation, (3) conformational denaturation over cycles.
- Literature data: Protein A leaching ~10-50 ng/mL eluate per cycle (Cytiva data). Metal leaching from IMAC: 0.1-1% Ni per cycle.
- A simplified exponential decay model is feasible: `activity(cycle) = activity_0 * exp(-k_deact * n_cycles)` with `k_deact` as a per-chemistry constant.
- **Priority: Defer to v6.0.** Long-term stability is a process economics concern, not a media design concern. v5.9 focuses on fresh-media performance.

**ARCHITECT (v6.0 sketch):**
- Add `k_deactivation_per_cycle: float` to `ReagentProfile`.
- Post-processing function `predict_lifetime(fmc, n_cycles)` -> `LifetimeProjection`.
- LOC estimate: ~100

**DEV-ORCHESTRATOR:**
- Work node: **WN-D16** (deferred to v6.0)
- Dependencies: None.

---

### Item D17: Washing Efficiency and Residual Toxic Reagent Tracking
**Source:** Audit docs 23, 26

**SCIENTIFIC ADVISOR:**
- Partially ready. Washing is a simple mass-transfer problem: residual reagent in pores is removed by diffusion/convection during wash steps. The key equation is `C_residual = C_initial * exp(-D_eff * t / L^2)` where L is the pore depth.
- For toxic reagents (ECH, glutaraldehyde, DVS), regulatory limits exist: ECH < 1 ppm (EP), glutaraldehyde < 5 ppm.
- **Priority: Nice-to-have for v5.9.** Important for GMP compliance but does not affect chromatographic performance. Simple diffusion model is implementable.

**ARCHITECT:**
- New `ModificationStepType.WASHING` enum value.
- Simple diffusion-out model: `_solve_washing_step()` computes residual concentration.
- New reagent-level field: `regulatory_limit_ppm: float`, `log_Kow: float` (partition coefficient).
- Track `residual_toxic_concentration: dict[str, float]` on `FunctionalMicrosphere`.
- LOC estimate: ~180

**DEV-ORCHESTRATOR:**
- Work node: **WN-D17**
- Model tier: Sonnet
- Dependencies: None.

---

## 3. Priority Classification

### v5.9 MUST-HAVE (6 items, ~1110 LOC)

| ID | Item | LOC | Rationale |
|----|------|-----|-----------|
| D1 | Metal charging step type | 250 | IMAC is #1 His-tag method; placeholder metal_loaded_fraction=1.0 is wrong |
| D2 | binding_model_hint routing in M3 | 180 | M3 currently ignores M2 output; all isotherms default to Langmuir |
| D3 | accessible_area_per_bed_volume on FMC | 40 | q_max underestimated 2-5x for porous beads |
| D5 | EDC/NHS chemistry for AHA -COOH | 200 | Completes AHA spacer pathway (noted as "not modeled") |
| D11 | Protein disulfide reduction step | 200 | Maleimide-thiol profiles require reduced Cys; no mechanism exists |
| D17 | Washing / residual toxic tracking | 180 | GMP compliance; ECH/glutaraldehyde are toxic |

### v5.9 NICE-TO-HAVE (4 items, ~540 LOC)

| ID | Item | LOC | Rationale |
|----|------|-----|-----------|
| D4 | Activity uncertainty propagation | 60 | Data exists; plumbing only |
| D6 | pH-dependent rate scaling | 250 | Quantitative improvement; rankings already correct |
| D8 | Irreversible affinity binding | 80 | Simple flag for streptavidin-biotin |
| D13 | Target-specific default assumptions | 150 | Usability improvement |

### Deferred to v6.0 (5 items)

| ID | Item | Reason |
|----|------|--------|
| D7 | Salt-dependent HIC isotherms | Requires user calibration data; Langmuir placeholder is functional |
| D9 | TMAE via DVS | Kinetic data for VS-tertiary amine unavailable |
| D10 | Lectin-specific M3 elution models | Sugar competition constants poorly characterized |
| D14 | Batch-to-batch DDA variability | Requires M1 interface changes |
| D15 | Pore-size distribution | Requires M1 PSD data output |

### Deferred to v6.0+ (2 items)

| ID | Item | Reason |
|----|------|--------|
| D16 | Protein leaching/deactivation | Long-term stability, not fresh-media design |
| D12 | Metal loading detail (beyond D1) | Metal leaching per-cycle is v6.0 scope |

---

## 4. Work Node Breakdown and Dependency Graph

```
Phase 1 (Parallel, no dependencies):
  WN-D1  [Metal Charging]        ── Opus   ~250 LOC
  WN-D2  [M3 Hint Routing]       ── Sonnet ~180 LOC
  WN-D3  [a_v_accessible on FMC] ── Haiku  ~40  LOC
  WN-D11 [Protein Pretreatment]  ── Sonnet ~200 LOC
  WN-D17 [Washing Step]          ── Sonnet ~180 LOC

Phase 2 (Depends on Phase 1 completion):
  WN-D5  [EDC/NHS Chemistry]     ── Opus   ~200 LOC
    depends on: WN-D3 (needs corrected a_v for NHS-ester q_max mapping)

Phase 3 (Nice-to-have, parallel after Phase 1):
  WN-D4  [Uncertainty Propagation] ── Haiku ~60 LOC
  WN-D6  [pH Rate Scaling]         ── Sonnet ~250 LOC
    depends on: WN-D5 (NHS hydrolysis is pH-sensitive; test together)
  WN-D8  [Irreversible Binding]    ── Sonnet ~80 LOC
    depends on: WN-D2 (routing must exist)
  WN-D13 [Target Defaults]         ── Haiku ~150 LOC
    depends on: WN-D2 (M3 must know which isotherm to populate)
```

### Dependency Graph (ASCII)

```
                    +---------+
                    | WN-D1   |   Metal Charging
                    +---------+
                         |
                    (standalone)

  +---------+      +---------+      +---------+
  | WN-D3   |----->| WN-D5   |      | WN-D2   |
  | a_v FMC |      | EDC/NHS |      | M3 Route|
  +---------+      +---------+      +---------+
                        |                |
                        v                |-----> WN-D8  [Irreversible]
                   +---------+           |-----> WN-D13 [Defaults]
                   | WN-D6   |
                   | pH Scale|
                   +---------+

  +---------+      +---------+      +---------+
  | WN-D11  |      | WN-D17  |      | WN-D4   |
  | Protein  |      | Washing |      | Uncert. |
  | Pretreat|      |         |      |         |
  +---------+      +---------+      +---------+
  (standalone)     (standalone)     (standalone)
```

---

## 5. Detailed Implementation Specifications

### WN-D1: Metal Charging Step

**Files modified:**
- `modification_steps.py`: Add `METAL_CHARGING` to `ModificationStepType`, add `_solve_metal_charging_step()`.
- `reagent_profiles.py`: Add `metal_association_constant` field to `ReagentProfile`. Add 5 new profiles (`nickel_charging`, `cobalt_charging`, `copper_charging`, `zinc_charging`, `edta_stripping`).
- `orchestrator.py`: Update `_STEP_ALLOWED_RTYPES`, `_STEP_ALLOWED_REACTION_TYPES`. Update FMC builder to read metal_loaded_fraction from metal charging result.
- `acs.py`: No changes (metal state is not an ACS site type; it modifies the chelator profile's `metal_loaded_fraction`).

**Solver logic:**
```
metal_loaded_fraction = K_assoc * C_metal / (1 + K_assoc * C_metal)
```
This is a simple Langmuir equilibrium on the chelator sites. The result updates the `metal_loaded_fraction` on the IDA/NTA ACS profile (stored as a new field on `ACSProfile` or as a metadata dict).

**Test cases:**
1. NTA + 50 mM NiSO4 -> metal_loaded_fraction > 0.99
2. EDTA stripping -> metal_loaded_fraction < 0.01
3. Integration test: ECH -> IDA -> Ni-charging -> breakthrough with His-tag protein

---

### WN-D2: M3 Binding Model Hint Routing

**Files modified:**
- `module3_performance/orchestrator.py`: Add `_select_isotherm()` dispatch, modify `run_breakthrough()` and `run_gradient_elution()` to use FMC hint.

**Routing table:**

| binding_model_hint | Isotherm Class | Parameters from FMC |
|---|---|---|
| `charge_exchange` | `SMAIsotherm` | Lambda from charge_density, z/sigma from target defaults |
| `metal_chelation` | `IMACCompetitionIsotherm` | q_max from estimated_q_max, K from target defaults |
| `salt_promoted` | `LangmuirIsotherm` | q_max from estimated_q_max (placeholder) |
| `fc_affinity` | `LangmuirIsotherm` | q_max from estimated_q_max |
| `gst_glutathione` | `LangmuirIsotherm` | q_max, K ~1e3 |
| `near_irreversible` | `LangmuirIsotherm` + irreversible flag (WN-D8) | q_max |
| `kappa_light_chain_affinity` | `LangmuirIsotherm` | q_max |
| `lectin_*` | `LangmuirIsotherm` (placeholder) | q_max, warn: requires calibration |
| `mixed_mode` | `SMAIsotherm` | Lambda, warn: simplified |
| `""` (empty) | `LangmuirIsotherm` | q_max from estimated_q_max or default |

---

### WN-D5: EDC/NHS Chemistry

**Files modified:**
- `acs.py`: Add `NHS_ESTER` to `ACSSiteType`.
- `reagent_profiles.py`: Add `edc_nhs_activation` and `nhs_amine_coupling` profiles.
- `modification_steps.py`: Handle CARBOXYL -> NHS_ESTER activation through existing ACTIVATION workflow (no new step type needed; ACTIVATION already supports arbitrary target/product pairs).

**Chemistry:**
```
COOH + EDC -> O-acylisourea (fast, unstable)
O-acylisourea + NHS -> NHS-ester (stable intermediate)
NHS-ester + R-NH2 -> amide bond (coupling)
NHS-ester + H2O -> COOH + NHS (hydrolysis, competing)
```

Modeled as: One ACTIVATION step (CARBOXYL -> NHS_ESTER) with hydrolysis_rate for EDC, followed by standard LIGAND_COUPLING on NHS_ESTER target.

---

### WN-D11: Protein Disulfide Reduction

**Files modified:**
- `modification_steps.py`: Add `PROTEIN_PRETREATMENT` to `ModificationStepType`, add solver.
- `reagent_profiles.py`: Add `tcep_reduction` and `dtt_reduction` profiles.
- New dataclass `ProteinPretreatmentState` (in `modification_steps.py` or new file).

**Key design decision:** This step does NOT modify ACS state. It modifies a `protein_state` dict carried alongside `acs_state` through the orchestrator. The downstream PROTEIN_COUPLING solver reads `thiol_accessibility_fraction` from this state instead of from the static reagent profile.

---

### WN-D17: Washing Step

**Files modified:**
- `modification_steps.py`: Add `WASHING` to `ModificationStepType`, add `_solve_washing_step()`.
- `reagent_profiles.py`: Add `regulatory_limit_ppm` and `diffusion_coefficient` fields. Add `wash_buffer` profile.
- `orchestrator.py`: Track `residual_concentrations: dict[str, float]` on `FunctionalMicrosphere`.

**Solver logic:**
```
# Diffusion-out from pore volume
C_residual = C_initial * exp(-D_eff * t_wash / (L_pore^2 / pi^2))
# where L_pore = pore_depth ~ bead_radius, D_eff = D_free * porosity / tortuosity
# Convective term: multiply by (1 + Pe) where Pe = v * L / D_eff
passes_regulatory = all(C_residual[reagent] < regulatory_limit)
```

---

## 6. New Enum Values Summary

### ACSSiteType (acs.py): +1
- `NHS_ESTER = "nhs_ester"` (D5)

### ModificationStepType (modification_steps.py): +3
- `METAL_CHARGING = "metal_charging"` (D1)
- `PROTEIN_PRETREATMENT = "protein_pretreatment"` (D11)
- `WASHING = "washing"` (D17)

### New Reagent Profiles: +9
- `nickel_charging`, `cobalt_charging`, `copper_charging`, `zinc_charging` (D1)
- `edta_stripping` (D1)
- `edc_nhs_activation` (D5)
- `tcep_reduction`, `dtt_reduction` (D11)
- `wash_buffer` (D17)

**Post-v5.9 profile count: 51** (42 existing + 9 new)

---

## 7. Effort Estimates

### v5.9 Must-Have

| Work Node | LOC | Complexity | Estimated Hours |
|-----------|-----|------------|-----------------|
| WN-D1 | 250 | Medium | 6 |
| WN-D2 | 180 | Medium | 4 |
| WN-D3 | 40 | Low | 1 |
| WN-D5 | 200 | High | 6 |
| WN-D11 | 200 | Medium | 4 |
| WN-D17 | 180 | Medium | 4 |
| **Subtotal** | **1050** | | **25 hours** |

### v5.9 Nice-to-Have

| Work Node | LOC | Complexity | Estimated Hours |
|-----------|-----|------------|-----------------|
| WN-D4 | 60 | Low | 1 |
| WN-D6 | 250 | Medium | 5 |
| WN-D8 | 80 | Low | 2 |
| WN-D13 | 150 | Low | 3 |
| **Subtotal** | **540** | | **11 hours** |

### Integration Testing (all v5.9 items)

| Test Suite | Estimated Hours |
|------------|-----------------|
| Unit tests per WN | 8 |
| Integration: ECH -> IDA -> Ni -> breakthrough | 3 |
| Integration: ECH -> AHA -> EDC/NHS -> protein | 3 |
| Integration: ECH -> DADPA -> SM(PEG)4 -> TCEP -> Protein-Cys | 2 |
| Regression: all 42 existing profiles | 2 |
| **Subtotal** | **18 hours** |

### Total v5.9 Effort

| Category | Hours |
|----------|-------|
| Must-have implementation | 25 |
| Nice-to-have implementation | 11 |
| Testing | 18 |
| Documentation and audit prep | 4 |
| **Total** | **58 hours** |

---

## 8. Risk Register

| Risk | Impact | Mitigation |
|------|--------|------------|
| EDC/NHS two-step mechanism may need ODE Template 4 (sequential activation) | Medium | Fall back to pseudo-single-step (combined EDC+NHS as one reagent) if ODE stiffness is problematic |
| Metal charging equilibrium assumption may not hold at very low metal concentrations | Low | Add kinetic ODE as backup; equilibrium is valid above 1 mM metal salt |
| pH scaling (D6) may break existing test fixtures | Medium | Add `pKa_nucleophile=0.0` default (disabled) so all existing profiles are unaffected |
| PROTEIN_PRETREATMENT introduces a new state carrier alongside ACS | Medium | Keep design minimal: single `protein_state` dict, not a full parallel state machine |
| Washing model assumes well-mixed pore volume | Low | Adequate for ranking; flag as `semi_quantitative` confidence |

---

## 9. Acceptance Criteria for v5.9

1. **All 6 must-have items pass their integration test suites.**
2. **Profile count reaches 51.** All new profiles have non-zero kinetic parameters with literature citations.
3. **M3 routing is functional.** Given an FMC with `binding_model_hint="charge_exchange"`, M3 uses SMA (not Langmuir). Given `"metal_chelation"`, M3 uses IMACCompetitionIsotherm.
4. **IMAC end-to-end test passes.** ECH activation -> IDA coupling -> Ni charging -> FMC -> M3 breakthrough with His-tag protein using IMAC isotherm.
5. **EDC/NHS path validated.** ECH -> AHA spacer arm -> EDC/NHS activation -> amine coupling produces non-zero functional ligand density.
6. **Maleimide-thiol path extended.** TCEP reduction step -> maleimide-thiol coupling uses reduction-derived thiol_accessibility_fraction (not static default).
7. **Washing step tracks residual ECH.** After 3x column volume wash at standard flow, residual ECH < 1 ppm for the reference bead (d50=90 um, porosity=0.6).
8. **No regression.** All existing 42-profile test cases continue to pass.
9. **Conservation invariants hold.** All ACSProfile.validate() checks pass for all new workflows.
10. **FMC `accessible_area_per_bed_volume` is populated** and used in q_max calculations (replaces external-only a_v).

---

## 10. v6.0 Preview (Deferred Items Roadmap)

| Item | v6.0 Priority | Blocking Dependency |
|------|---------------|---------------------|
| D7: Salt-dependent HIC | High | User calibration framework |
| D9: TMAE via DVS | Medium | Kinetic data for VS-tertiary amine |
| D10: Lectin elution models | Medium | Sugar competition K values |
| D14: Batch-to-batch variability | High | M1 uncertainty interface |
| D15: Pore-size distribution | High | M1 PSD output |
| D16: Protein leaching/deactivation | Medium | Cycle-dependent kinetic data |
| D12: Metal leaching detail | Low | Subsumed by D16 |

---

*End of v5.9 plan. Document 28 of EmulSim Module 2 design series.*
