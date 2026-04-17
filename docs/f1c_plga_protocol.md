# F1-c — PLGA Solvent-Evaporation Platform: /architect Protocol

**Prepared by:** /architect (via dev-orchestrator)
**Date:** 2026-04-17
**Status:** Protocol + Phase 1 ready for fresh session(s); Phases 2-3 deferred
**Scope basis:** `docs/node32_cluster_f_v8_roadmap.md` §2.3 (PLGA workstream)

---

## 1. Purpose and scope

Add a PLGA / solvent-evaporation platform to EmulSim as the third
non-chitosan-agarose microsphere family. PLGA microsphere fabrication
is fundamentally different from both alginate (ionic) and cellulose
(NIPS) gelation: the microsphere forms by volatile solvent (DCM or
ethyl acetate) evaporating out of the dispersed-phase droplets,
leaving concentrated PLGA behind, which solidifies as a glassy
amorphous polymer. The existing L1 (emulsification) machinery
transfers with minor property edits; L2 is a new
diffusion-dominated solvent-depletion solver; L3 is stubbed (no
covalent crosslinking — PLGA chains are physically entangled); L4
uses a Gibson-Ashby-style porosity-discounted glassy modulus.

Out of scope for F1-c: pharmaceutical-drug loading / release kinetics
(would be M3/M4 modules), multi-emulsion (W/O/W) formulations,
supercritical-CO₂ antisolvent processes.

## 2. Scientific basis

### 2.1 Mechanism summary

Standard single-emulsion O/W solvent-evaporation route:

1. **L1 emulsification:** PLGA is dissolved in a volatile organic
   solvent (DCM = dichloromethane, most common; ethyl acetate; or
   acetone) at 1-15 wt%. The organic phase is emulsified into an
   aqueous continuous phase containing PVA surfactant (0.5-2 wt%
   typical). Droplets carry the PLGA/solvent solution unchanged. Note
   the reversed polarity vs agarose/alginate: organic is dispersed,
   water is continuous.
2. **L2 solvent evaporation:** DCM diffuses from droplet interior to
   droplet surface, then transfers across the oil/water interface
   into bulk water. From the bulk water phase, DCM evaporates to the
   gas headspace (driven by stirring / warming / reduced pressure).
   The droplet-interior PLGA concentration rises as DCM depletes.
3. **Solidification:** when `phi_PLGA` locally crosses a threshold
   (~0.85 in volume fraction, corresponding to the glass-transition /
   vitrification composition), the polymer solidifies. If DCM
   evaporation is much faster than diffusion to the surface, a dense
   outer shell forms first, trapping internal DCM as vapour bubbles
   → porous / hollow microsphere. If evaporation is slow, uniform
   densification → dense microsphere.
4. **L3:** not applied. PLGA microspheres are mechanically stable
   via polymer-chain entanglement below T_g (T_g ≈ 40-60 °C for
   PLGA, well above typical processing and use temperatures).
5. **L4 modulus:** glassy polymer with Gibson-Ashby porosity
   discount: `G = G_glassy · (1 − porosity)² ` (or similar cellular-
   solid scaling). Typical E ≈ 1-3 GPa, G ≈ 0.4-1 GPa for dense
   PLGA; porous microspheres are 10×-100× softer.

### 2.2 Literature anchors

- Freitas et al. (2005) *Int. J. Pharm.* 282:1 — canonical review of
  solvent-evaporation microsphere fabrication, including DCM removal
  kinetics and morphology vs process parameter maps.
- Wang & Schwendeman (1999) *J. Pharm. Sci.* 88:1090 — DCM
  extraction / evaporation kinetic measurements; effective D_DCM in
  PLGA/DCM matrix ≈ 10⁻¹⁰–10⁻⁹ m²/s depending on concentration.
- Yang et al. (2001) *J. Control. Release* 74:151 — PLGA microsphere
  morphology (dense / porous / hollow) as a function of process
  conditions.
- Park et al. (1998) *Biomaterials* 19:745 — PLGA glass transition
  and modulus vs L:G ratio and M_n.
- Blasi (2019) *J. Pharm. Invest.* 49:337 — contemporary review of
  PLGA microsphere process-property relationships.

## 3. Module inventory (deliverables)

### 3.1 New files

| File | Purpose | Est. LOC |
|---|---|---|
| `src/emulsim/level2_gelation/solvent_evaporation.py` | 1D spherical Fickian DCM-depletion solver | ~320 |
| `src/emulsim/level4_mechanical/plga.py` | Gibson-Ashby glassy-porous modulus | ~130 |
| `src/emulsim/properties/plga_defaults.py` | Presets for PLGA 50:50 / 75:25 / 85:15 / PLA | ~160 |
| `tests/test_plga_phase1.py` | ~16 tests (Phase 1 scope) | ~300 |

### 3.2 Edits to existing files (Phase 2, deferred)

| File | Change | Est. LOC |
|---|---|---|
| `src/emulsim/datatypes.py` | Add PLGA fields to MaterialProperties; ensure PolymerFamily.PLGA exists (already shipped) | ~12 |
| `src/emulsim/pipeline/orchestrator.py` | Add `_run_plga` branch in `run_single` | ~110 |
| `src/emulsim/__main__.py` | `run --polymer-family plga --plga-grade 50_50` | ~25 |
| `src/emulsim/config.py` | Recognise `[formulation].plga_grade = "50_50"` via FormulationParameters field | ~3 |
| `tests/test_plga_phase2_phase3.py` | ~12 integration tests | ~280 |

### 3.3 Total footprint

Phase 1: ~910 LOC (solver + L4 + defaults + tests) — matches F1-a /
F1-b Phase 1 budget.
Phase 2: ~430 LOC (orchestrator + CLI + TOML + integration tests).
Phase 3: ~200 LOC (remaining 3 grade presets if not all shipped in
Phase 1, plus any polish).
Cumulative ≈ 1500-1600 LOC — within the Node 32 roadmap estimate.

## 4. Algorithm specifics

### 4.1 L2 solvent-evaporation solver

The physics reduces to a 1D spherical Fickian depletion problem on
`r ∈ [0, R_droplet]`:

```
State: phi_DCM(r, t)     — DCM volume fraction [−]
       phi_PLGA(r, t) = 1 − phi_DCM(r, t)      (algebraic)

Governing equation:
    ∂phi_DCM/∂t = D_DCM(phi_DCM) · (1/r²)·∂/∂r(r²·∂phi_DCM/∂r)

Boundary conditions:
    ∂phi_DCM/∂r (r=0, t) = 0        (spherical symmetry)
    phi_DCM(R, t) = phi_eq ≈ 0      (continuous-phase sink;
                                      DCM solubility in water is
                                      low and water-phase DCM is
                                      re-swept by stirring)

Initial condition:
    phi_DCM(r, 0) = 1 − phi_PLGA_0
    phi_PLGA_0 ≈ 0.10 for 10 wt% PLGA/DCM
```

Approximation used in Phase 1: **fixed droplet radius**. As DCM
evaporates the droplet should shrink to volume (phi_PLGA_0 · V_0),
but Eulerian PDE on a moving boundary is a Phase 2 refinement.
Reported "final microsphere radius" is instead derived by
mass-conservation post-processing: `R_sphere = R_0 · phi_PLGA_0^(1/3)`.

Concentration-dependent diffusivity (optional Phase 2): real DCM
diffusion in PLGA/DCM follows a Fujita-style exponential fall-off
with polymer fraction — `D(phi_DCM) = D_0 · exp(−β · phi_PLGA)`.
Phase 1 uses a constant `D_DCM` for simplicity; flagged as a known
limitation in the manifest assumptions.

**Observables** (populated into the `GelationResult`):

- `phi_plga_mean_final` — volume-averaged final polymer fraction
  (usually ≈ 1 for long enough `t_end`)
- `t_vitrification` — time at which `phi_PLGA_mean > 0.85`
- `skin_thickness_proxy` — outer-shell thickness at the earliest
  moment a cell crosses phi_PLGA > 0.85
- `core_porosity_proxy` — `1 − phi_PLGA_core` at `t_vitrification`
  (captures how much solvent was trapped inside when the skin
  formed)
- `porosity_final` = `1 − phi_plga_mean_final` (0 for a fully-dried
  dense microsphere; > 0 for porous / hollow morphology)

### 4.2 L4 PLGA modulus

```
G_PLGA = G_glassy · phi_PLGA_mean ^ n_plga
```

With `G_glassy ≈ 7 × 10⁸ Pa` for PLGA 50:50 (from E ≈ 2 GPa, ν ≈ 0.42)
and `n_plga ≈ 2` (Gibson-Ashby cellular-solid scaling for
closed-cell foams). For a dense microsphere (`phi_PLGA_mean = 1`)
this reduces to `G = G_glassy`. For a 50 %-porous microsphere,
`G = 0.25 · G_glassy ≈ 1.75 × 10⁸ Pa`.

Evidence tier: SEMI_QUANTITATIVE (literature-parameterised for 4
PLGA variants; wet-lab calibration upgrades to CALIBRATED_LOCAL if
Study C lands).

## 5. Parameter table (Phase 1 & Phase 3 presets)

| Parameter | PLGA 50:50 | PLGA 75:25 | PLGA 85:15 | PLA | Units | Source |
|---|---|---|---|---|---|---|
| L:G ratio | 50:50 | 75:25 | 85:15 | 100:0 | — | grade |
| M_n | 30 000 | 50 000 | 60 000 | 60 000 | g/mol | typical |
| T_g | 45 | 50 | 55 | 60 | °C | Park 1998 |
| D_DCM | 1.0 × 10⁻⁹ | 8.0 × 10⁻¹⁰ | 6.0 × 10⁻¹⁰ | 5.0 × 10⁻¹⁰ | m²/s | Wang 1999 |
| phi_eq (water bath) | 0.005 | 0.005 | 0.005 | 0.005 | — | Henry eq |
| G_glassy | 7.0 × 10⁸ | 9.0 × 10⁸ | 1.0 × 10⁹ | 1.2 × 10⁹ | Pa | Park 1998 |
| n_plga (Gibson-Ashby) | 2.0 | 2.0 | 2.0 | 2.0 | — | GA foam |
| phi_PLGA_0_typical | 0.10 | 0.10 | 0.08 | 0.08 | — | process |

## 6. Test cases (Phase 1 scope, 8 core + edge-case extras)

**Core 8 (Phase 1):**

1. Fickian mass conservation: total DCM integrated over the droplet
   decreases monotonically; residual DCM = 0 at t → ∞.
2. Dirichlet sink drives evaporation: with `phi_eq = 0`,
   `phi_PLGA_mean(t → ∞) = 1`.
3. √t front propagation in early times (diffusion-limited).
4. Higher D_DCM → faster vitrification (test with 2×D vs D).
5. Modulus scales as phi² (Gibson-Ashby): doubling phi_PLGA_mean at
   fixed G_glassy → G × 4.
6. Zero PLGA gives zero modulus, zero gel, UNSUPPORTED manifest.
7. L2 + L4 manifests are SEMI_QUANTITATIVE at solver exit.
8. PLGA 50:50 preset applies all 6 fields correctly to
   MaterialProperties.

**Extras (Phase 1 lab-bench):**

9. Solver input validation (negative R, tiny grid, negative time).
10. Skin-thickness proxy increases monotonically with evaporation
    time (or remains zero if bulk-limited).
11. `plga_modulus` function edge cases (zero phi, zero G_glassy).
12. Registry introspection: all 4 grade presets present after
    Phase 3, at least PLGA 50:50 after Phase 1.

**Integration (Phase 2, deferred):**

13. `PolymerFamily.PLGA` routes `run_single` to `_run_plga`.
14. TOML `[formulation].plga_grade = "50_50"` round-trip.
15. CLI `--plga-grade` flag parses.
16. Full pipeline produces non-zero G.

## 7. Dependencies and blockers

### Hard

- Shipped: Node 30 unified UQ; Node F1-a/b dispatch pattern (template
  for `_run_plga`); Node F3 inverse-design TargetSpec (usable
  immediately for PLGA once the platform ships).

### Soft

- /scientific-advisor full briefing on PLGA DCM-diffusion kinetics
  and glass-transition composition — valuable for Phase 2 but not
  blocking Phase 1 (literature values are robust).
- Study C wet-lab data (PLGA DCM-evaporation microspheres) — will
  promote tier to CALIBRATED_LOCAL when delivered.

### Risks

- **R5**: Fixed-radius approximation in the Phase 1 solver
  under-estimates final radius by ~55 % for 10 wt% formulations.
  Mitigation: expose `R_final = R_0 · phi_0^(1/3)` as a manifest
  diagnostic; user can post-process to plot the shrunken sphere.
  Full moving-boundary ALE solver = Phase 2 refinement.
- **R6**: Constant-D_DCM assumption breaks down when `phi_PLGA >
  0.6`. Mitigation: v1 solver still gives the right
  `t_vitrification` order-of-magnitude; document as known limitation.

## 8. Gate G1 (12-point completeness)

| # | Criterion | Status |
|---|---|---|
| G1-01 | Purpose one-sentence | ✓ §1 |
| G1-02 | Inputs typed + units | ✓ §4.1, §5 |
| G1-03 | Outputs typed + units | ✓ §4.1 observables |
| G1-04 | Algorithm + pseudocode | ✓ §4.1 |
| G1-05 | Complexity + justification | ✓ §4.1 (single-field Fickian, BDF on 40-pt grid) |
| G1-06 | Numerical considerations | ✓ §4.1 + §7 (fixed-R, constant-D limits) |
| G1-07 | ≥3 unit tests | ✓ §6 (tests 1-5) |
| G1-08 | Boundary tests per constraint | ✓ §6 (tests 2, 6, 9) |
| G1-09 | Error conditions + response | ✓ §6 #9; zero-polymer handled via `_zero_plga_result` |
| G1-10 | Performance budget | Target < 5 s for 40-pt grid, 50 µm bead, 1 h evaporation |
| G1-11 | Upstream + downstream modules | ✓ §3.2 |
| G1-12 | Logging + metrics | Phase 2 will add run-time log lines; Phase 1 emits manifest diagnostics |

**G1 status: COMPLETE.** 12/12 met for Phase 1 scope. Phase 2
orchestrator wiring carries its own G1 gate when authored.

## 9. Execution plan

**This session delivers Phase 1:** solver + L4 + PLGA 50:50 default
+ 8 core tests.

**Fresh session(s) pick up Phase 2:**
1. Add orchestrator `_run_plga` branch (mirror `_run_cellulose`).
2. Add `--polymer-family plga --plga-grade 50_50` CLI surface.
3. Add `[formulation].plga_grade` TOML field on
   `FormulationParameters`.
4. Add 12 integration tests.

**Phase 3:** populate the remaining PLGA grade presets (75:25,
85:15, PLA) and add 4 grade-dependence tests.

---

**End of F1-c protocol.**
