# EmulSim Modules 2 & 3: Final Implementation Plan

**Date:** 2026-04-11
**Synthesised from:** Scientific Advisor, Dev-Orchestrator, Architect
**Input:** Original design (doc 11) + 3rd-party audit (doc 12)
**Status:** Implementation-ready architecture

---

## Cross-Role Consensus

All 3 roles independently confirmed:

1. **All 10 audit findings are valid** — 7 HIGH, 3 MEDIUM, no false positives
2. **The audit's reframing is correct** — "validation-gated, tiered process simulator"
3. **The audit's Phase 0-E build order is strictly superior** to our original
4. **Module 2 is NOT a light L3 extension** — needs a sequential reaction-network simulator
5. **Surface area model is the #1 prerequisite** — without it, all ACS densities are undefined
6. **Start minimal (2 workflows, 1-component breakthrough), validate, then expand**

### Key Corrections to Original Design

| Original Claim | Audit + Our Assessment |
|----------------|----------------------|
| "Module 2 heavily reuses L3" | Module 2 needs competitive hydrolysis ODEs, steric blocking, and multi-step state tracking beyond L3 |
| "~4,650 LOC total" | ~5,270-8,500 LOC with proper surface area, pressure-drop, and transport tiering |
| "Build all libraries upfront" | Only 4 reagents for Phase B; rest deferred until core accounting validated |
| "LRM is sufficient" | Need 3-tier transport (EDM/LRM/GRM) with auto-selection via Biot number |
| "Pressure-drop is a trust gate mention" | Kozeny-Carman with compressibility feedback is a first-class solver |

---

## Final Build Order (6 Phases)

### Phase 0: Module 1 Dependency Stabilization (~170 LOC)
- `M1ExportContract` dataclass with trust-level tagging
- Uncertainty propagation metadata for downstream inheritance
- **Gate:** M1 outputs round-trip without data loss

### Phase A: Data Model + Surface Area Core (~750 LOC)
- `AccessibleSurfaceModel` (3 tiers: external_only, empirical_pore, morphology_based)
- `ACSProfile` with full hierarchy (total/accessible/activated/consumed/blocked/functional)
- `FunctionalMicrosphere`, `MobilePhaseComposition`, `ColumnGeometry`
- ACS conservation tests, surface-area validation
- **Gate:** Conservation invariant holds for 5 synthetic sequences; surface area within 20% of literature

### Phase B: Minimal Module 2 — 2 Workflows (~1,030 LOC)
- Chemistry-agnostic reaction engine with 3 ODE templates (simple 2nd-order, competitive hydrolysis, equilibrium binding)
- Workflow 1: Amine coupling on chitosan (genipin/glutaraldehyde secondary crosslinking)
- Workflow 2: Hydroxyl activation on agarose (ECH/DVS → epoxy/vinyl sulfone)
- Small-molecule ligand coupling, protein/enzyme coupling with activity model
- Only 4 reagent profiles (not 20)
- **Gate:** ACS conservation holds; protein coupling gives active < coupled; all 10 M2 trust gates fire

### Phase C: Minimal Module 3 Chromatography (~1,060 LOC)
- Kozeny-Carman pressure-drop with compressibility feedback
- Langmuir isotherm (single-component)
- Lumped Rate Model (LRM) transport: conservative FV, method-of-lines + BDF
- UV detection (Beer-Lambert + Gaussian broadening)
- Breakthrough curve extraction with DBC at 5%/10%/50%
- Mass balance validator
- **Gate:** Matches Thomas model analytical solution within 5%; mass balance < 1%; pressure within 10% of Ergun

### Phase D: Catalytic Packed Bed (~670 LOC)
- Plug-flow + axial dispersion (reuses Phase C FV engine)
- Michaelis-Menten with generalized Thiele modulus for effectiveness factor
- First-order enzyme deactivation
- **Gate:** Conversion matches analytical PFR within 5%; effectiveness factor within 2% of exact

### Phase E: Multi-Component + Gradients + Detection (~1,590 LOC)
- Competitive Langmuir, SMA (IEX), IMAC competition, pH-dependent Protein A
- Gradient generator (linear/step/multi-linear)
- Multi-component transport extension
- Fluorescence, conductivity, semi-quantitative MS detection
- Process metrics (yield, purity, resolution, productivity)
- Expanded reagent/ligand libraries
- Visualization (ACS waterfall, chromatograms, breakthrough overlays)
- Pipeline integration + UI tabs
- **Gate:** 2-component separation produces resolved peaks; all 20 trust gates active

---

## Key Architectural Decisions

| Decision | Rationale |
|----------|-----------|
| Single orchestrator extended (not separate) | Backward compatible; single FullResult |
| ACSProfile as state machine | Sequential chemistry demands full state tracking |
| 3-tier transport mirrors L2 tiering | Users apply same reasoning: empirical→mechanistic |
| ColumnGeometry owns pressure-drop | Every query checks mechanical feasibility |
| IsothermModel ABC with validate_conditions() | Each isotherm knows its validity envelope |
| DetectorModel splits ideal vs instrument response | Clean separation of physics vs artifacts |
| Surface area tier hooks into L2 morphology_descriptors() | Leverages existing phase-field infrastructure |

---

## Trust Gates: 20 Checks

**Module 2 (10):** ACS conservation (BLOCKER), accessibility ceiling (BLOCKER), activation ceiling, coupling ceiling, pH compatibility, pore exclusion (reagent), pore exclusion (ligand), excessive pore blockage, secondary XL capacity loss, residual reactive groups

**Module 3 (10):** Mechanical pressure limit (BLOCKER), solute pore exclusion, grid resolution, mass balance closure (BLOCKER), non-negative concentrations (BLOCKER), bed Reynolds number, Peclet number range, binding capacity vs ligand density, detector saturation, unsupported isotherm warning

---

## Total Effort

| Phase | LOC | Opus | Sonnet | Haiku |
|-------|-----|------|--------|-------|
| 0 | 170 | 0 | 3 | 1 |
| A | 750 | 2 | 4 | 0 |
| B | 1,030 | 2 | 5 | 1 |
| C | 1,060 | 1 | 5 | 1 |
| D | 670 | 2 | 2 | 0 |
| E | 1,590 | 2 | 9 | 3 |
| **Total** | **5,270** | **9** | **28** | **6** |

---

> **Strategic Principle:** Validate before extending, calibrate before predicting,
> be honest about what the model can and cannot do.
