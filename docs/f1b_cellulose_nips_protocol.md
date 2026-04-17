# F1-b — Cellulose NIPS Platform: /architect Protocol

**Prepared by:** /architect (via dev-orchestrator)
**Date:** 2026-04-17
**Status:** Protocol only; implementation deferred to fresh session(s)
**Scope basis:** `docs/node32_cluster_f_v8_roadmap.md` §2.2

---

## 1. Purpose and scope

Add a cellulose / NIPS (non-solvent-induced phase separation) platform
to EmulSim as the second non-chitosan-agarose microsphere family.
Cellulose processing is mechanistically very different from both
agarose (thermal) and alginate (ionic) gelation: the gel forms by
polymer-solvent phase separation after the microsphere is exposed to a
non-solvent (water, typically). The existing L1 (emulsification)
machinery mostly transfers; L2 (gelation) requires a wholly new
Cahn–Hilliard-style coupled-diffusion solver; L3 is omitted by
default (NIPS is the gelation); L4 uses a different modulus law keyed
to polymer volume fraction.

Out of scope for F1-b: PLGA solvent evaporation (F1-c), cellulose
surface functionalisation, post-NIPS covalent crosslinking.

## 2. Scientific basis

### 2.1 Solvent systems considered

Four industrially-relevant cellulose solvent systems:

| Solvent | Typical `c_cell` | Non-solvent | Processability | Tier |
|---|---|---|---|---|
| NaOH/urea (7-9 wt% / 12 wt%) | 2-6 wt% | water | Easy; freezing-step pre-dissolution | Good v1 target |
| NMMO (N-methylmorpholine N-oxide, 50-80 wt%) | 5-15 wt% | water | Industrial Lyocell; T ≈ 90 °C | Good v1 target |
| EMIM-Ac (ionic liquid) | 3-12 wt% | water / ethanol | Research-scale; T ≈ 80 °C | Stretch |
| DMAc/LiCl | 3-10 wt% | water | Research; complex ternary | Skip v1 |

**F1-b v1 defaults to NaOH/urea** (widest literature base,
room-temperature processing, single-step freeze-dissolution).
Parameter table (§5) holds defaults for all four; the solver is
written to handle any ternary polymer/solvent/non-solvent system.

### 2.2 Mechanism summary

1. **Initial state:** cellulose dissolved in solvent (homogeneous, single phase).
2. **Emulsification (L1):** same PBE framework; cellulose-solvent
   droplets in paraffin oil (or silicone oil) with Span-80 or similar.
   Droplets carry the cellulose-solvent mixture unchanged — no phase
   separation yet.
3. **NIPS (L2):** droplets enter a non-solvent bath (water). Non-solvent
   diffuses in, solvent diffuses out. As non-solvent fraction rises,
   cellulose solubility drops. When the local composition crosses the
   binodal, phase separation starts; when it crosses the spinodal,
   spontaneous decomposition occurs. Final morphology:
   - Shallow quench (near binodal) → cellular / closed-pore.
   - Deep quench (past spinodal) → bicontinuous (preferred for
     chromatography beads — higher surface area, open porosity).
4. **L3:** not applied by default (NIPS gel is the end-state). Optional
   downstream covalent crosslinking (e.g., ECH, Node 9 library) plugs
   in via the existing L3 if the user selects it.
5. **L4 modulus:** empirical scaling
   `G = K_cell · phi_cell^alpha` where `alpha ≈ 2` in the semi-dilute
   regime (Rubinstein-Colby) and `alpha ≈ 9/4` in the entangled
   regime. `K_cell` depends on solvent and molecular weight.

### 2.3 Literature anchors

- Cuissinat & Navard (2006) *Macromol. Symp.* 244:1 — dissolution
  mechanisms in NaOH/urea, NMMO, IL.
- Xu et al. (2010) *Biomacromolecules* 11:1724 — cellulose/NaOH-urea
  gelation kinetics + NIPS morphologies.
- Lindman et al. (2010) *Phys. Chem. Chem. Phys.* 12:4369 — cellulose
  dissolution thermodynamics, Flory-Huggins χ parameters.
- Young et al. (2018) *Polymer* 134:76 — Cahn-Hilliard modelling of
  cellulose NIPS spherulite formation.
- Zhang et al. (2020) *Cellulose* 27:1071 — modulus vs phi_cell
  empirical scaling for regenerated cellulose aerogels/hydrogels.

## 3. Module inventory (deliverables)

### 3.1 New files

| File | Purpose | Est. LOC |
|---|---|---|
| `src/emulsim/level2_gelation/nips_cellulose.py` | Coupled-diffusion + Cahn-Hilliard NIPS solver | ~550 |
| `src/emulsim/level4_mechanical/cellulose.py` | phi^alpha modulus law | ~140 |
| `src/emulsim/properties/cellulose_defaults.py` | Defaults for 4 solvent systems | ~180 |
| `src/emulsim/reagent_library_cellulose.py` | SolventProfile + NonSolventProfile entries | ~220 |
| `tests/test_cellulose_nips.py` | ~18 tests | ~380 |

### 3.2 Edits to existing files

| File | Change | Est. LOC |
|---|---|---|
| `src/emulsim/datatypes.py` | Extend `PolymerFamily` use (enum value already exists); add `MaterialProperties` fields: `cellulose_solvent`, `chi_polymer_solvent`, `chi_polymer_nonsolvent`, `D_solvent`, `D_nonsolvent`, `K_cell`, `alpha_cell`, `M_n_cellulose` | ~40 |
| `src/emulsim/pipeline/orchestrator.py` | Add `PolymerFamily.CELLULOSE` branch in `run_single`: `_run_cellulose(...)` mirrors `_run_alginate` shape | ~110 |
| `src/emulsim/config.py` | Recognise `[formulation].solvent_system = "naoh_urea"` etc., `[formulation].nonsolvent_bath = "water"` | ~30 |
| `src/emulsim/__main__.py` | `run --polymer-family cellulose --solvent-system naoh_urea --nonsolvent water` | ~30 |

### 3.3 Total footprint

~1470 LOC new + ~210 LOC edits + ~380 LOC tests ≈ **~2060 LOC**,
somewhat over the Node 32 roadmap estimate (~1500-2000 LOC).
Practical plan: ship v1 with NaOH/urea only (skip the 4-system
defaults table for Phase 1), yielding ~1400 LOC + 260 LOC tests;
add the three remaining solvent presets in an F1-b.1 follow-up.

## 4. Algorithm specifics

### 4.1 L2 NIPS solver (core new work)

The 1D spherical coupled PDE system on `r ∈ [0, R_droplet]`:

```
State fields:
    phi(r, t)   — cellulose volume fraction [-]
    s(r, t)     — solvent (NaOH/urea, NMMO, …) volume fraction [-]
    ns(r, t) = 1 − phi − s   — non-solvent (water) volume fraction [-]

Coupled-diffusion equations (ternary):
    dphi/dt = ∇·(M_phi · ∇(mu_phi))
    ds/dt   = ∇·(M_s   · ∇(mu_s))

where mu_phi, mu_s are chemical potentials derived from
Flory-Huggins free energy:

    f(phi, s) = (phi/N_p)·ln(phi)
              + s·ln(s)
              + ns·ln(ns)
              + chi_PS · phi · s
              + chi_PN · phi · ns
              + chi_SN · s · ns

mu_phi = df/dphi − ns·df/dns       (Gibbs-Duhem)
mu_s   = df/ds   − ns·df/dns

The Cahn-Hilliard correction adds −κ·∇²phi to mu_phi (regularises the
interface; κ = 6σ·ξ / polymer_density where ξ is the interfacial width).

Boundary conditions:
    phi(R, t),  s(R, t)  :  Dirichlet matching bath composition
                             (phi = 0, s = 0, ns = 1 for pure water bath)
    dphi/dr(0, t) = ds/dr(0, t) = 0   :  symmetry
    phi(r, 0) = phi_0,  s(r, 0) = s_0,  ns(r, 0) = 0

Time integration: scipy.integrate.solve_ivp with BDF (stiff due to the
4th-order CH term after semi-discretisation). Radial grid of ~60
points; the CH term requires slightly more refinement than the Fickian
L2 solvers. Mobility M_phi, M_s are scalar constants at v1 (can be
composition-dependent in a follow-up).
```

**Observables for `GelationResult`:**

- `phi_cell_mean` — bulk-averaged cellulose volume fraction (high ≈ denser gel)
- `pore_size_mean` — from correlation-length analysis of phi(r) (FFT / autocorrelation)
- `porosity` = 1 − phi_cell_mean
- `bicontinuous_score` — from spatial-derivative histogram of phi
  (high if spinodal decomposition, low if cellular)
- `t_gel_complete` — time when `max|dphi/dt|` drops below 10⁻⁴ /s
  (composition frozen)

### 4.2 L4 cellulose modulus

```
G = K_cell · phi_cell^alpha_cell
```

With:
- `K_cell ≈ 0.5 × 10⁶ Pa` for NaOH/urea regenerated cellulose (Zhang 2020)
- `alpha_cell = 2.25` (entangled-regime default; use 2.0 for semi-dilute)
- `phi_cell_mean` from L2 output

Evidence tier: SEMI_QUANTITATIVE (literature-parameterised for 3
solvent systems; wet-lab calibration upgrades to CALIBRATED_LOCAL if
Study B lands).

## 5. Parameter table (NaOH/urea default)

Populate in `properties/cellulose_defaults.py`:

| Parameter | NaOH/urea | NMMO | EMIM-Ac | Units | Source |
|---|---|---|---|---|---|
| `phi_0` (initial cellulose vol frac) | 0.03-0.06 | 0.05-0.12 | 0.03-0.10 | — | process-set |
| `chi_PS` (χ polymer-solvent) | 0.45 | 0.40 | 0.38 | — | Lindman 2010 |
| `chi_PN` (χ polymer-nonsolvent = water) | 0.85 | 0.80 | 0.80 | — | Flory-Huggins fit |
| `chi_SN` (χ solvent-nonsolvent) | 0.30 | 0.20 | 0.25 | — | regular-solution |
| `D_solvent` (solvent self-diffusion in gel) | 5e-11 | 2e-11 | 1e-11 | m²/s | Xu 2010 |
| `D_nonsolvent` (water in gel) | 1e-10 | 8e-11 | 5e-11 | m²/s | Xu 2010 |
| `M_n_cellulose` | 60 000 | 80 000 | 60 000 | g/mol | batch-typical |
| `K_cell` | 5e5 | 8e5 | 4e5 | Pa | Zhang 2020 |
| `alpha_cell` | 2.25 | 2.25 | 2.25 | — | Rubinstein-Colby |
| `kappa_CH` (gradient energy coef) | 1e-17 | 1e-17 | 1e-17 | J·m⁻¹ | κ = 6σ·ξ order-of-magnitude |

All values carry SEMI_QUANTITATIVE until Study B' (cellulose wet-lab
campaign) promotes to QUANTITATIVE.

## 6. Test cases (≥11 required)

1. **Dispatch**: `MaterialProperties(polymer_family=CELLULOSE)` routes
   L2 to `nips_cellulose` solver via the orchestrator.
2. **Mass conservation**: ∫(phi + s + ns) dV = V_total at all times
   (the three fractions must sum to 1 everywhere).
3. **Initial-state stability**: for `chi_PS` below the binodal and no
   non-solvent exposure, phi(r, t) = phi_0 for all t (no phase separation).
4. **Water-bath driven demixing**: `ns(R) = 1` Dirichlet drives
   non-solvent ingress; phi_mean(t → large) > phi_0 (densification
   because cellulose precipitates locally).
5. **Binodal crossing produces morphology**: for a deep-quench recipe,
   the final phi(r) profile has bimodal spatial distribution (dense +
   sparse regions), quantified by `bicontinuous_score > 0.5`.
6. **Shallow quench → cellular** (bicontinuous_score < 0.3).
7. **Modulus scaling**: doubling `phi_cell_mean` at fixed `K_cell` and
   `alpha_cell = 2.0` → G × 4.
8. **Zero cellulose**: phi_0 = 0 → phi_cell_mean = 0 → G = 0.
9. **Config round-trip**: TOML with `polymer_family = "cellulose"` and
   `[formulation].solvent_system = "naoh_urea"` loads correctly.
10. **Evidence tier**: L2 manifest reports SEMI_QUANTITATIVE.
11. **Integration**: full pipeline `run_single` for a cellulose config
    produces a non-null `FullResult` with non-zero G, non-trivial
    `porosity` in [0.4, 0.95].
12. **Solvent-system selection**: switching from NaOH/urea to NMMO
    changes both L2 transport constants and L4 K_cell in manifest
    diagnostics.
13. **Regression**: existing alginate, agarose-chitosan, and EDC/NHS
    paths unaffected (smoke).

## 7. Dependencies and blockers

### Hard dependencies

- **Shipped**: Node 30 unified UQ; Node F1-a Phase 2a/2b/2c (the
  PolymerFamily dispatch pattern is the template for the cellulose
  branch); Node F3 inverse-design TargetSpec (usable immediately for
  cellulose once platform ships).
- **NOT yet done**: /scientific-advisor full briefing on cellulose NIPS
  mechanism + Flory-Huggins parameters + Cahn-Hilliard closure.
  **BLOCKING** for implementation kickoff.
- Study B' wet-lab data (cellulose NIPS) — SOFT; platform ships at
  SEMI_QUANTITATIVE without it, upgrades on delivery.

### Risks

- **R3** (Node 32 roadmap): Cahn-Hilliard numerics can be stiff and
  expensive. Mitigation: use semi-implicit Fourier-spectral rather
  than finite-volume if performance becomes a blocker (adds ~150 LOC
  and a `numpy.fft` dependency but is 5-10× faster).
- **R4**: ternary Flory-Huggins free energy has multiple local minima
  — solver may get stuck in metastable states for bad initial
  conditions. Mitigation: add a small noise term to `phi(r, 0)` (1%
  amplitude) to break symmetry, or run a short high-mobility
  "equilibration" phase at t < 0.

### Soft dependencies

- /chief-economist NPV: is cellulose chromatography media a
  commercially viable target vs alginate/PLGA/agarose?
- /ip-auditor FTO: NMMO cellulose processing has significant Lenzing
  patent overhang; NaOH/urea (Zhang Lab Wuhan) is more clearly free.

## 8. Gate G1 (12-point completeness)

| # | Criterion | Status |
|---|---|---|
| G1-01 | Purpose one-sentence | ✓ §1 |
| G1-02 | Inputs typed + units + constraints | ✓ §4.1, §5 |
| G1-03 | Outputs typed + units | ✓ §4.1 observables |
| G1-04 | Algorithm + pseudocode | ✓ §4.1 |
| G1-05 | Complexity + justification | ✓ §4.1 (60-pt grid, BDF + CH term) |
| G1-06 | Numerical considerations | ✓ §4.1, §7 R3 (stiffness + metastability) |
| G1-07 | ≥3 unit tests | ✓ §6 (tests 2-5) |
| G1-08 | Boundary tests per constraint | ✓ §6 (tests 2, 3, 8) |
| G1-09 | Error conditions + response | ✗ **Needs /scientific-advisor input — NaN trapping, solver divergence on deep quenches** |
| G1-10 | Performance budget | Partial — target < 30 s for 60-pt grid, 500 µm bead, 5 min NIPS |
| G1-11 | Upstream + downstream modules | ✓ §3.2 |
| G1-12 | Logging + metrics | ✗ Deferred to Phase 2 draft (`t_gel_complete`, `bicontinuous_score` log lines per timestep) |

**G1 status: PARTIAL.** 10/12 criteria met; items 09 and 12 require
/scientific-advisor consultation and Phase 2 authorship decisions.

## 9. Execution plan (fresh-session checklist)

1. Read `docs/node32_cluster_f_v8_roadmap.md` §2.2 + this protocol.
2. Invoke `/scientific-advisor` with §2 mechanism + §5 parameter table;
   request full literature review of cellulose NIPS thermodynamics
   (binodal / spinodal loci) and Cahn-Hilliard numerics choices.
3. Invoke `/chief-economist` (briefly) for cellulose-vs-alginate/PLGA
   NPV ranking.
4. Upgrade G1 to 12/12 with advisor input (items 09 and 12).
5. Execute Phase 2 in the order:
   - `datatypes.py` edits (MaterialProperties new fields)
   - `properties/cellulose_defaults.py` (NaOH/urea only for v1)
   - `level2_gelation/nips_cellulose.py` — **Opus-tier**, single
     largest file, expect 2 fix rounds
   - `level4_mechanical/cellulose.py`
   - `reagent_library_cellulose.py`
   - `pipeline/orchestrator.py` `_run_cellulose` branch (copy from
     `_run_alginate`)
   - `config.py` + `__main__.py` surface wiring
   - Tests in `tests/test_cellulose_nips.py`
6. Run full regression; expect ≥230 tests to pass (alginate +
   EDC/NHS + inverse-design + UQ + new cellulose tests).

### Phased split (recommended)

| Phase | Scope | LOC | Sessions |
|---|---|---|---|
| F1-b Phase 1 | L2 solver + L4 modulus + NaOH/urea defaults + schema + 8 core tests | ~900 | 2-3 |
| F1-b Phase 2 | Orchestrator dispatch + reagent library + config/CLI + 5 integration tests | ~500 | 1-2 |
| F1-b Phase 3 (polish) | NMMO + EMIM-Ac + DMAc/LiCl defaults + 4 extra tests | ~350 | 1 |

Total effort estimate: **4-6 fresh sessions**.

---

**End of F1-b protocol. Implementation deferred.**
