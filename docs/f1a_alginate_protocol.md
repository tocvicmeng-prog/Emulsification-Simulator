# F1-a — Alginate Platform: /architect Protocol + Handover

**Prepared by:** /architect (via dev-orchestrator)
**Date:** 2026-04-17
**Status:** Protocol only; implementation deferred to fresh session(s)
**Scope basis:** `docs/node32_cluster_f_v8_roadmap.md` §2.1

---

## 1. Purpose and scope

Add an alginate / Ca²⁺ ionic-gelation platform to EmulSim as the first
non-chitosan-agarose microsphere family. Alginate is the closest
migration of the five v8.0 F1 platforms: L1 (emulsification) reuses
most of the existing machinery; L2 (gelation) needs a new
diffusion-limited ionic-gelation solver; L3 and L4 are similar to the
chitosan/agarose path but with alginate-specific parameters.

Out of scope for F1-a: cellulose NIPS (F1-b), PLGA (F1-c), surface
functionalisation of alginate (v8.1+).

## 2. Scientific basis

**Mechanism summary** (requires /scientific-advisor full briefing at
implementation kickoff):

1. **L1:** Same PBE framework; alginate droplets in paraffin oil with
   appropriate surfactant (Span 80 / Tween 80 blend). Key differences:
   - `sigma` (IFT) ~ 3-7 mN/m, lower than agarose-chitosan
   - `mu_d` (alginate solution viscosity) strongly shear-thinning
     (Cross model)
   - No thermal gelation during emulsification (alginate is stable in
     hot oil); gelation happens only when Ca²⁺ is introduced.

2. **L2:** Diffusion-limited ionic gelation.
   - Ca²⁺ diffuses from an external bath (or internally-released
     source) into the alginate droplet.
   - Binding stoichiometry: ~2 carboxylate groups per Ca²⁺
     ("egg-box" junction zone).
   - Front propagation: shrinking-core model with r_gel(t) ∝ √(D·t).
   - Output observables: `porosity`, `pore_size_mean`, gelation time
     `t_gel`, residual un-gelled core fraction.

3. **L3:** Typically omitted for pure alginate (ionic gelation IS the
   crosslinking). For composite alginate systems (Ca²⁺-alginate +
   secondary covalent crosslinker) a separate L3 step can run; the
   existing generic L3 dispatch is reusable.

4. **L4:** Mechanical modulus from the alginate-guluronate (G-block)
   density. Empirical: `G_DN ∝ (c_alginate · f_G)²` where f_G is the
   guluronate fraction (typical 0.3-0.7 depending on source).

**Literature anchors** (for /scientific-advisor validation):
- Braccini & Pérez 2001 *Biomacromolecules* 2:1089 — egg-box model
- Kuo & Ma 2001 *Biomaterials* 22:511 — Ca²⁺-alginate mass transport
- Fernandez-Grajera et al. 2022 (review) — alginate microsphere
  fabrication parameters
- Kong et al. 2004 *Macromolecules* 37:6838 — modulus vs composition

## 3. Module inventory (deliverables)

### 3.1 New files

| File | Purpose | Est. LOC |
|---|---|---|
| `src/emulsim/level2_gelation/ionic_ca.py` | Shrinking-core diffusion-limited Ca²⁺ front solver | ~400 |
| `src/emulsim/level4_mechanical/alginate.py` | Modulus from G-block content | ~150 |
| `src/emulsim/properties/alginate_defaults.py` | Default material properties for alginate (IFT, viscosity curves, f_G, D_Ca, binding stoich) | ~150 |
| `src/emulsim/reagent_library_alginate.py` | Alginate crosslinker/source entries (CaCl₂, internal gelation with GDL + CaCO₃) | ~200 |
| `tests/test_alginate_platform.py` | ~15 tests covering L1 adaptation, L2 ionic gel, L4 modulus, integration | ~300 |

### 3.2 Edits to existing files

| File | Change | Est. LOC |
|---|---|---|
| `src/emulsim/datatypes.py` | Add `PolymerFamily` enum (AGAROSE_CHITOSAN | ALGINATE); `MaterialProperties.polymer_family`; alginate-specific fields | ~30 |
| `src/emulsim/pipeline/orchestrator.py` | Dispatch L2 by `polymer_family` | ~40 |
| `src/emulsim/level1_emulsification/solver.py` | No-op for alginate unless viscosity curve differs | ~0 |
| `src/emulsim/config.py` | TOML parser recognises `polymer_family = "alginate"` key | ~20 |

### 3.3 Total footprint estimate

~1500 LOC new + ~90 LOC edits + ~300 LOC tests ≈ **~1900 LOC**.
Matches the Node 32 roadmap estimate.

## 4. Algorithm specifics

### 4.1 L2 ionic-Ca front solver (core new work)

```
1D spherical shrinking-core in r ∈ [0, R_droplet]:

    Ca²⁺ concentration C_Ca(r, t):
        ∂C/∂t = D·(1/r²)·∂/∂r(r²·∂C/∂r) − 2·k_bind·C·G(r, t)
    Guluronate concentration G(r, t):
        ∂G/∂t = − k_bind·C·G
    Egg-box crosslink density X(r, t):
        dX/dt = ½·k_bind·C·G

Boundary conditions:
    C(R, t) = C_Ca_bath (external bath concentration)
    ∂C/∂r(0, t) = 0 (symmetry)
    G(r, 0) = G_0 (initial alginate concentration × f_G)
    C(r < R, 0) = 0
    X(r, 0) = 0

Reach criterion for "fully gelled":
    r such that G(r, t) < 0.01 · G_0 → gel front passed that r.
```

Scipy `solve_ivp` with method `BDF` or `LSODA` on a radial grid of
~50 points. Stiffness is moderate (diffusion + second-order reaction);
same integrator family as existing L3 reaction-diffusion solver.

**Observables for `GelationResult`:**

- `pore_size_mean` — correlated with residual G_0 fraction (low X =
  large pores between partial egg-box aggregates)
- `porosity` — polymer volume fraction
- `t_gel_complete` — time for front to reach r=0 (when reached; else
  "partial" flag)
- `egg_box_density_mean` — bulk average of X(r, t_final)

### 4.2 L4 alginate modulus

```
G_DN = K_alg · (c_alginate · f_G)^n_alg · (X_mean / X_max)
```

with `K_alg ≈ 30 kPa · (g/L)^-n_alg`, `n_alg ≈ 2.0` (Kong et al. 2004).
`X_max = c_alginate · f_G / 2` (stoichiometric maximum). The final
factor captures incomplete gelation.

## 5. Parameter table (alginate defaults)

Populate in `properties/alginate_defaults.py`:

| Parameter | Value | Units | Source |
|---|---|---|---|
| `sigma` | 5.0 × 10⁻³ | N/m | Typical Span-80 vs alginate |
| `mu_d` (zero-shear) | 0.5-5.0 | Pa·s | Concentration-dependent |
| `f_G` (guluronate fraction) | 0.3-0.7 | — | Batch-dependent; 0.5 default |
| `D_Ca` (Ca²⁺ in alginate gel) | 1.0 × 10⁻⁹ | m²/s | Kuo & Ma 2001 |
| `k_bind` (Ca²⁺ + 2 COO⁻) | ~10³ | M⁻²·s⁻¹ | Braccini & Pérez 2001 |
| `C_Ca_bath` (typical CaCl₂) | 100 | mol/m³ (100 mM) | Process standard |
| `K_alg` (modulus prefactor) | 30 × 10³ | Pa | Kong et al. 2004 |
| `n_alg` (modulus exponent) | 2.0 | — | Kong et al. 2004 |

All values require /scientific-advisor review at implementation
kickoff; flag at SEMI_QUANTITATIVE tier until Study A' (alginate
wet-lab campaign) promotes to QUANTITATIVE.

## 6. Test cases (11 required)

1. **Polymer family dispatch**: `MaterialProperties(polymer_family=ALGINATE)` routes L2 to the ionic solver.
2. **Mass conservation**: `∫ C_Ca dV` + 2·`∫ X dV` equals cumulative Ca²⁺ inflow.
3. **Shrinking-core front**: r_gel(t) ∝ √t (log-log slope 0.5 ± 0.1) in the diffusion-limited regime.
4. **Complete gel at infinite time**: `X(r, t → ∞) → X_max` for all r < R, for any non-zero C_Ca_bath.
5. **Zero Ca**: `C_Ca_bath = 0` → `X_mean = 0` everywhere, G_DN = 0.
6. **Modulus scaling**: doubling `c_alginate` at fixed f_G → G_DN roughly ×4 (n_alg=2.0).
7. **Guluronate fraction**: doubling f_G at fixed c_alginate → G_DN ~×4.
8. **Smoke test**: 100 µm droplet + 100 mM CaCl₂ bath + 1 h reaches ≥95 % gelation.
9. **Config round-trip**: TOML with `polymer_family = "alginate"` loads into `SimulationParameters` with the alginate L2 solver selected.
10. **Evidence tier**: alginate L2 manifest reports SEMI_QUANTITATIVE (literature-parameterised, not yet calibrated).
11. **Integration**: full pipeline `run_single` with an alginate config produces a non-null `FullResult` with non-zero G_DN.

## 7. Dependencies and blockers

### Hard dependencies
- Shipped: Node 30 unified UQ; Node 31 L3 chemistry expansion
  pattern; Node F3-a TargetSpec for inverse design across platforms.
- **NOT yet done**: /scientific-advisor briefing on alginate mechanism
  + parameter validation. **BLOCKING for implementation kickoff.**

### Soft dependencies
- /chief-economist NPV of alginate vs PLGA / cellulose — determines
  whether F1-a ships as v8.0 Phase 2 (planned) or gets rescheduled.
- /ip-auditor FTO review on alginate microsphere claims if commercial
  use is targeted.

### Risks (from Node 32 roadmap §6)

- **R2**: Alginate ionic gelation may diverge too much from current L2
  → full L2 rewrite risk. Mitigation: prototype `ionic_ca.py` as
  standalone (no shared state with existing L2) first.

## 8. Gate G1 (12-point completeness)

| # | Criterion | Status |
|---|---|---|
| G1-01 | Purpose one-sentence | ✓ §1 |
| G1-02 | Inputs typed + units + constraints | ✓ §4.1, §5 |
| G1-03 | Outputs typed + units | ✓ §4.1 "observables" |
| G1-04 | Algorithm + pseudocode | ✓ §4.1 |
| G1-05 | Complexity + justification | ✓ §4.1 (50-point radial grid, BDF stiff integrator) |
| G1-06 | Numerical considerations | ✓ §4.1 (stiffness, mass-conservation check) |
| G1-07 | ≥3 unit tests | ✓ §6 (tests 1-4) |
| G1-08 | Boundary tests per constraint | ✓ §6 (tests 2, 5) |
| G1-09 | Error conditions + response | ✗ **Needs /scientific-advisor input for realistic failure modes** |
| G1-10 | Performance budget | Partial — needs benchmarking target |
| G1-11 | Upstream + downstream modules | ✓ §3.2 |
| G1-12 | Logging + metrics | ✗ **Deferred to Phase 2 draft** |

**G1 status: PARTIAL.** 10/12 criteria met; items 09 and 12 require
/scientific-advisor consultation + Phase 2 authorship decisions
before full G1 pass.

## 9. Execution plan (fresh-session checklist)

A session picking up F1-a should:

1. Read `docs/node32_cluster_f_v8_roadmap.md` §2.1 + this protocol.
2. Invoke `/scientific-advisor` with the §2 "Mechanism summary" +
   §5 parameter table for validation. Request: full literature review
   of ionic-Ca gelation in alginate microspheres; confirm or refine
   all numerical values in §5 and §4.2 of this protocol.
3. Invoke `/chief-economist` (briefly) for alginate-vs-PLGA NPV
   ranking if not already decided.
4. Upgrade G1 to 12/12 by filling §8 items 09 and 12 with advisor's
   input.
5. Execute Phase 2 (implementation) in the order:
   - `datatypes.py` edits (PolymerFamily enum + field)
   - `properties/alginate_defaults.py`
   - `level2_gelation/ionic_ca.py` (~400 LOC, Opus-tier recommended)
   - `level4_mechanical/alginate.py`
   - `reagent_library_alginate.py`
   - `pipeline/orchestrator.py` dispatch edits
   - `config.py` TOML parser hook
   - 11 tests in `tests/test_alginate_platform.py`
6. Run full regression; expect ≥200 tests to pass (Node 30/31/F3/F4
   baseline + new alginate tests).

Total effort estimate: **3-5 sessions** (one for each of solver
implementation, dispatch wiring, test authoring, integration, and
polishing/docs).

---

**End of F1-a protocol. Implementation deferred.**
