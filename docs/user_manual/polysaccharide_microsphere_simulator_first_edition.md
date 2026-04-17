# Polysaccharide-Based Microsphere Emulsification Simulator — First Edition

**Edition:** 1.0 (v8.3 feature set)
**Date:** 2026-04-17
**Audience:** Downstream-processing technicians, junior R&D researchers,
first-time users of the EmulSim pipeline.

---

## Preface

This manual is written for someone who has never before fabricated a
polysaccharide microsphere or run a multi-scale process simulator.
It assumes college-level chemistry and basic Python literacy; no prior
experience with emulsification, hydrogel chemistry, or numerical
modelling is required.

The simulator's purpose is simple: **turn a written recipe into a
predicted product, before you touch the bench.** A recipe describes
what you put into the vessel (polymer, surfactant, crosslinker,
mixing conditions, temperature profile). A product is described by
measurable outputs (droplet size distribution, pore size, elastic
modulus, protein sieving coefficient). The simulator tells you what
to expect, with honest uncertainty bands, so your first wet-lab
attempt is an educated guess rather than a shot in the dark.

This manual is organised in three concentric layers. **Part I** is a
ten-minute read that teaches a first-time user to launch the UI,
enter a recipe, and interpret a result. **Part II** is a reference
for choosing between the four microsphere platforms the simulator
supports. **Part III** (appendices) is for users who want to
understand what the simulator is doing under the hood, verify the
science, or troubleshoot an unexpected result.

---

## Table of Contents

- **Part I — Getting Started**
  - 1. What the Simulator Does
  - 2. Workflow Overview
  - 3. Five-Minute Quickstart
- **Part II — Platform Catalogue**
  - 4. Agarose-Chitosan Double Network
  - 5. Alginate (Ca²⁺ Ionic Gelation)
  - 6. Cellulose NIPS (Non-Solvent-Induced Phase Separation)
  - 7. PLGA Solvent Evaporation
  - 8. Crosslinker / Gelant Selection
- **Part III — Appendices**
  - A. Detailed Input Requirements
  - B. Process Steps
  - C. Essential Input & Process Checklist
  - D. Frequently Asked Questions
  - E. Architectural Ideas & Working Principles
  - F. Chemical & Physical Principles
  - G. Formulas & Mathematical Theorems
  - H. Standard Wet-Lab Protocols
  - I. Troubleshooting Table

---

# Part I — Getting Started

## 1. What the Simulator Does

The EmulSim Polysaccharide Microsphere Simulator (below: "EmulSim")
is a multi-scale process simulator that predicts the outcome of
microsphere fabrication from first principles plus literature-derived
empirical laws. It covers four sequential stages, which it calls
**levels**:

| Level | Stage | Governs |
|---|---|---|
| L1 | Emulsification | Droplet size distribution |
| L2 | Gelation | Network formation, pore architecture |
| L3 | Crosslinking | Covalent or ionic crosslink density |
| L4 | Mechanical | Bulk modulus, protein sieving |

Not every level runs for every platform — for example, alginate skips
L3 because its ionic gelation **is** the crosslinking. The simulator
automatically selects the right sub-pipeline based on which polymer
family you pick.

### 1.1 What questions the simulator answers

- What will my droplet size distribution look like at a given RPM?
- What pore size and porosity will my gel have?
- Will my crosslinker reach the target crosslink density in the
  time allotted?
- What will the final bead modulus be?
- For chromatography: what size range of proteins will be sieved?
- How confident are we in each of the above?

### 1.2 What the simulator does not do

- It does not simulate downstream drug-release kinetics (that would
  be a Module 4 that is not yet built).
- It does not substitute for a regulatory tox / stability study.
- It cannot predict the behaviour of a novel polymer it has no
  parameters for. Supported polymer families in this edition are
  agarose-chitosan, alginate, cellulose, and PLGA.
- It cannot predict what will happen above the validation regime of
  its input parameters — it will warn you but not refuse to run.

### 1.3 Trust-tier vocabulary

Every result carries an **evidence tier** that tells you how much to
trust it. Users must learn these five tiers before reading any
result:

| Tier | Meaning | When to trust |
|---|---|---|
| VALIDATED_QUANTITATIVE | Predictions match wet-lab data within 10 % | Safe to use for process design |
| CALIBRATED_LOCAL | Fitted to local wet-lab data, not yet cross-validated | Safe for screening; validate before scale-up |
| SEMI_QUANTITATIVE | Literature-parameterised, plausible but un-calibrated | Rank-order is reliable; absolute values are indicative |
| QUALITATIVE_TREND | Directionally correct only | Use only for screening, not quantitative claims |
| UNSUPPORTED | Placeholder or degraded path | Do not use for decision-making |

The UI marks every reported value with its tier. In the first
edition of EmulSim, most outputs sit at **SEMI_QUANTITATIVE** —
literature-anchored but awaiting the Study A/B/C wet-lab
calibration campaigns that will promote them to VALIDATED status.

## 2. Workflow Overview

All runs follow the same five-step workflow regardless of platform.

```
 ┌─────────────────────┐
 │ 1. Choose Platform  │    agarose-chitosan  |  alginate
 └──────────┬──────────┘    cellulose NIPS    |  PLGA solvent-evap
            │
 ┌──────────▼──────────┐
 │ 2. Enter Formulation│    polymer wt%, surfactant, pH,
 └──────────┬──────────┘    temperature profile, t_gel, t_crosslink
            │
 ┌──────────▼──────────┐
 │ 3. Pick Crosslinker │    (or gelant / solvent system)
 └──────────┬──────────┘    genipin, Ca-bath, NaOH/urea, PLGA 50:50 …
            │
 ┌──────────▼──────────┐
 │ 4. Run Pipeline     │    L1 → L2 → L3 → L4
 └──────────┬──────────┘
            │
 ┌──────────▼──────────┐
 │ 5. Read Result      │    DSD, pore, G, Kav, trust tiers
 └─────────────────────┘
```

The UI implements exactly this flow. The CLI does it in one line.

## 3. Five-Minute Quickstart

### 3.1 Launching the UI

```bash
python -m emulsim ui
```

This opens a web browser at `http://localhost:8501` showing the
EmulSim dashboard. The upper-right corner of the page carries a
button labelled **Manual (PDF)** that downloads this document.

### 3.2 Entering a minimal recipe

The sidebar holds **Global Settings** (hardware mode, scientific
mode). The main panel has three tabs:

- **Module 1 — Emulsification & Particle Sizing**
- **Module 2 — Functionalisation & Crosslinking**
- **Module 3 — Mechanical & Chromatography**

A full run needs three numbers from Module 1 (RPM, emulsification
time, dispersed-phase fraction), two concentrations from Module 2
(polymer wt%, crosslinker concentration), and a "Run" click.
Everything else has a sensible default that the simulator prints
before running. First-time users should keep the defaults.

### 3.3 Choosing a platform (CLI example)

```bash
python -m emulsim run configs/default.toml \
    --polymer-family alginate \
    --gelant cacl2_external
```

### 3.4 Reading the output

The UI reports, for each level:

- L1: `d32` (Sauter mean diameter), `d50` (median), `span`
- L2: `pore_size_mean`, `porosity`, (platform-dependent) alpha
  or phi
- L3: `p_final` (crosslink conversion), `G_chit`
- L4: `G_DN` (effective shear modulus), `E_star` (reduced Young's
  modulus), `Kav` array (protein sieving)

Every value is tagged with its evidence tier and a small
provenance chip ("see Appendix E for model name"). **Before making
a recipe decision, confirm the tier of the specific value you are
using.**

### 3.5 Your first three runs

First-time users should do the following three runs on day one:

1. **Run the default recipe.** This is agarose-chitosan at 8000 RPM,
   genipin crosslinker, 24 h. Study the result.
2. **Double the RPM.** Observe how `d32` drops roughly as
   RPM⁻¹·². This teaches you the Kolmogorov-scaling intuition.
3. **Switch to alginate with `--gelant cacl2_external`.** Observe
   that L3 is stubbed (the pipeline tells you so) and that the
   modulus comes entirely from the L2 ionic-gel manifest.

After these three runs a new user is equipped to try their own
recipes.

---

# Part II — Platform Catalogue

The simulator supports four polymer families. Each has a different
gelation mechanism, different suitable crosslinkers, and different
typical applications. This section is the one-page "which do I
pick?" guide.

## 4. Agarose-Chitosan Double Network

**Default platform.** Suitable for: rigid chromatography beads,
cell-culture scaffolds, protein-resistant coatings.

**Mechanism:** hot agarose solution carrying dispersed chitosan is
emulsified in oil. On cooling below ~45 °C the agarose chains form
a thermally-gelled helix network. A subsequent amine-reactive
crosslinker (genipin, glutaraldehyde, epichlorohydrin) covalently
links chitosan into an interpenetrating second network. The
resulting double-network (DN) hydrogel is stiffer and more
chemically stable than either component alone.

**Typical recipe** (simulator defaults):

- Agarose: 4.2 wt% (42 kg/m³)
- Chitosan: 1.8 wt% (18 kg/m³)
- Span-80 surfactant: 2 wt%
- Oil: paraffin
- T_oil: 90 °C (during emulsification)
- Cooldown: ~10 °C/min
- Crosslinker: genipin 2 mM, 24 h at 37 °C

**Strengths:** rigid, well-validated, protein-compatible.
**Trade-offs:** slow crosslinking (24 h for genipin), requires
hot emulsification, blue colour from genipin.

## 5. Alginate (Ca²⁺ Ionic Gelation)

Suitable for: mild encapsulation (cells, probiotics, enzymes),
food-grade products, rapid prototyping.

**Mechanism:** alginate solution droplets contact Ca²⁺, which binds
pairs of guluronate ("G-block") residues in an "egg-box" junction,
immediately forming a gel. No covalent crosslinking.

**Sub-routes:**

- **External CaCl₂ bath** — fast (minutes), inhomogeneous gel
  (denser outer shell, softer core).
- **Internal GDL/CaCO₃ release** — slow (hours), homogeneous gel.

**Typical recipe:**

- Alginate: 2 wt% (20 kg/m³)
- f_guluronate: 0.5 (batch-dependent)
- CaCl₂ bath: 100 mM
- Process T: 25 °C (no thermal activation)
- t_gel: 30 min for 500 µm beads

**Strengths:** room-T processing, mild on biology, food-grade.
**Trade-offs:** ionic crosslinks dissolve in saline / pH extremes;
lower modulus (~10–50 kPa) vs covalent systems (~1 MPa).

## 6. Cellulose NIPS (Non-Solvent-Induced Phase Separation)

Suitable for: regenerated cellulose chromatography resins,
high-porosity scaffolds, sustainable materials.

**Mechanism:** cellulose is dissolved in a solvent (NaOH/urea,
NMMO, ionic liquid); droplets contact a non-solvent (water) which
diffuses in and precipitates cellulose via spinodal or binodal
phase separation. The resulting morphology depends on quench depth
— shallow quenches give cellular gels, deep quenches give
bicontinuous porous networks.

**Sub-routes (solvent systems, preset-selectable):**

- `naoh_urea` — aqueous, room-T, patent-clear (default)
- `nmmo` — Lyocell industrial, 90 °C, higher M_n
- `emim_ac` — ionic liquid, research scale
- `dmac_licl` — analytical-chemistry classic

**Typical recipe:**

- Cellulose: 5 wt% in NaOH/urea
- Non-solvent: water
- Process T: 25 °C
- t_gel: ≥10 min for 200 µm droplets

**Strengths:** bicontinuous porosity (open structure), cellulose's
excellent protein compatibility, no covalent crosslinker needed.
**Trade-offs:** solvent recovery dominates cost; the Cahn-Hilliard
L2 solver is the simulator's most mathematically demanding model.

## 7. PLGA Solvent Evaporation

Suitable for: drug-delivery microspheres, controlled-release
depot formulations, bioresorbable scaffolds.

**Mechanism:** PLGA is dissolved in a volatile organic solvent
(DCM); this is emulsified as droplets in an aqueous phase with
PVA surfactant. The DCM diffuses out into water and evaporates;
PLGA concentrates until it vitrifies into a glassy microsphere.
No crosslinking — mechanical stability comes from the glassy /
entangled state below T_g.

**Grade presets:** PLGA 50:50, 75:25, 85:15, PLA (100:0). Higher
lactide fraction → slower degradation, higher T_g, higher modulus.

**Typical recipe:**

- PLGA 50:50 dissolved in DCM: 10 wt%
- PVA aqueous phase: 1 wt%
- Process T: 25 °C (below all T_g)
- t_evaporation: 30 min to 4 h

**Strengths:** bioresorbable, FDA-familiar, tunable degradation.
**Trade-offs:** organic-solvent handling, potential for drug
degradation during evaporation, glassy mechanics (1 GPa+, not
squishy).

## 8. Crosslinker / Gelant Selection

| Chemistry | Appropriate for | Tier |
|---|---|---|
| Genipin (amine bridge) | Agarose-chitosan | SEMI_QUANTITATIVE |
| Glutaraldehyde | Agarose-chitosan | SEMI_QUANTITATIVE |
| EDC/NHS (amide bond) | COOH-bearing matrices only | SEMI_QUANTITATIVE with warning |
| Epichlorohydrin | Agarose + chitosan -OH | SEMI_QUANTITATIVE |
| TPP (ionic) | Chitosan (mild, reversible) | QUALITATIVE_TREND |
| Ca²⁺ (external bath) | Alginate | SEMI_QUANTITATIVE |
| GDL + CaCO₃ (internal) | Alginate (homogeneous) | SEMI_QUANTITATIVE |
| PEGDA + UV | PEG interpenetrating network | SEMI_QUANTITATIVE |
| Citric acid | Heat-cure ester bonds | QUALITATIVE_TREND |
| Divinyl Sulfone (DVS) | Agarose -OH | SEMI_QUANTITATIVE |

**Key warning on EDC/NHS:** this chemistry requires surface COOH
groups. Pure chitosan / pure agarose have essentially none. If you
select EDC/NHS without pre-modifying the matrix (e.g., carboxymethyl
chitosan or succinylation), the simulator runs a
QUALITATIVE_TREND-tier fallback and will not produce a trustworthy
modulus. The mechanistic EDC/NHS path activates only when
`MaterialProperties.surface_cooh_concentration > 0`.

---

# Part III — Appendices

## Appendix A. Detailed Input Requirements

### A.1 Global settings

| Field | Units | Typical range | Default | Notes |
|---|---|---|---|---|
| Hardware mode | — | rotor-stator / stirred-vessel | rotor-stator | Determines kernel parameters |
| Scientific mode | — | empirical / hybrid / mechanistic | hybrid | Hybrid is the default |
| Output directory | path | — | `./output` | Created if absent |

### A.2 Emulsification (Level 1)

| Field | Units | Range | Default | Role |
|---|---|---|---|---|
| `rpm` | rpm | 800–25 000 | 8 000 | Impeller speed |
| `t_emulsification` | s | 60–3 600 | 300 | Mixing duration |
| `phi_d` | — | 0.01–0.40 | 0.05 | Dispersed-phase volume fraction |
| `T_oil` | K | 298–373 | 363 | Oil-phase temperature |
| `cooling_rate` | K/s | 0.05–1.0 | 0.167 | Cool-down rate |
| `l1_n_bins` | — | 20–200 | 60 | Radial PBE bins |
| `l1_d_min` | m | 1e-7 | 1e-7 | Smallest tracked droplet |
| `l1_d_max` | m | 1e-3 | 1e-3 | Largest tracked droplet |

### A.3 Formulation

| Field | Units | Range | Default | Notes |
|---|---|---|---|---|
| `c_agarose` | kg/m³ | 10–80 | 42 | Primary for default platform |
| `c_chitosan` | kg/m³ | 5–30 | 18 | Primary for default platform |
| `c_alginate` | kg/m³ | 5–40 | 0 | F1-a platform |
| `phi_cellulose_0` | — | 0.02–0.15 | 0 | F1-b platform |
| `phi_PLGA_0` | — | 0.05–0.20 | 0 | F1-c platform |
| `c_span80` | kg/m³ | 5–40 | 20 | Surfactant |
| `c_genipin` | mol/m³ | 0.5–10 | 2 | Crosslinker (chitosan) |
| `c_Ca_bath` | mol/m³ | 50–500 | 100 | Alginate external bath |
| `t_crosslink` | s | 300–86 400 | 86 400 | Crosslinking / gelation time |
| `T_crosslink` | K | 277–353 | 310 | Crosslinking temperature |
| `pH` | — | 3–10 | 7 | EDC/NHS + future pH-dependent |

### A.4 Platform / preset selectors

| Field | TOML key | CLI flag | Allowed values |
|---|---|---|---|
| Polymer family | `[material].polymer_family` | `--polymer-family` | `agarose_chitosan`, `alginate`, `cellulose`, `plga` |
| Alginate gelant | `[formulation].gelant` | `--gelant` | `cacl2_external`, `gdl_caco3_internal` |
| Cellulose solvent | `[formulation].solvent_system` | `--cellulose-solvent` | `naoh_urea`, `nmmo`, `emim_ac`, `dmac_licl` |
| PLGA grade | `[formulation].plga_grade` | `--plga-grade` | `50_50`, `75_25`, `85_15`, `pla` |

### A.5 Material properties

Most users never need to edit MaterialProperties directly; the
simulator populates sensible defaults from the PropertyDatabase and
the selected preset. Power users who wish to override a single
parameter (e.g., raise `k_xlink_0` for a novel chitosan grade)
should pass it via `props_overrides` in the programmatic API or via
a `[material]` section in the properties TOML. The full field list
and units are in `src/emulsim/datatypes.py`.

## Appendix B. Process Steps

The typical first-time workflow, in order:

1. **Define the target.** Write down what you want the product to
   do. Example: "500 µm alginate beads with 70 % porosity for
   immobilised enzyme catalysis."
2. **Select the polymer family** based on the target (see §4-7).
3. **Pick the crosslinker or gelant** consistent with the family
   (see §8 and the table in Appendix A.4).
4. **Choose the hardware** (rotor-stator for small droplets,
   stirred vessel for larger). Set RPM and t_emulsification.
5. **Enter the formulation** (polymer wt%, surfactant, crosslinker
   concentration, pH, temperatures, times).
6. **Run the simulation.** Read the evidence tier of every output.
7. **Iterate.** Adjust RPM to hit a target droplet size; adjust
   crosslinker concentration to hit a target modulus; adjust
   emulsification time to achieve convergence (see FAQ 4).
8. **Export the ProcessDossier** (CLI `dossier` subcommand or UI
   "Export Dossier" button) and attach it to your notebook entry.
9. **Go to the bench.** Run the predicted recipe. Measure d32,
   pore, modulus.
10. **Return to the simulator.** Drop the measured values into
    `data/validation/<level>/assays/` as an AssayRecord. The next
    simulator run will consume them via `CalibrationStore` and
    promote tiers from SEMI_QUANTITATIVE toward CALIBRATED_LOCAL.

## Appendix C. Essential Input & Process Checklist

Before clicking Run, confirm every item below.

**Setup**
- [ ] Hardware mode matches your real-world vessel.
- [ ] Scientific mode is `hybrid` (first-time users keep the default).
- [ ] Output directory exists and is writable.

**Platform**
- [ ] Polymer family matches the recipe intent.
- [ ] Correct preset / gelant / grade selected for the family.
- [ ] Polymer concentration is populated (only one family at a time).

**Emulsification**
- [ ] RPM within the hardware's rated range.
- [ ] `t_emulsification` long enough for DSD convergence (typically
      ≥ 5 min; see FAQ 4).
- [ ] `phi_d` in the 0.01–0.40 range (outside this, PBE kernels
      lose accuracy).

**Formulation**
- [ ] Surfactant concentration is reasonable for the chosen
      polymer / oil combination.
- [ ] Temperatures are internally consistent (T_oil > T_gel for
      agarose; T_crosslink compatible with the crosslinker profile).
- [ ] Crosslinker concentration is not zero unless you intended a
      no-crosslinker run.
- [ ] For PLGA / cellulose: the CLI `--polymer-family` flag and the
      `phi_X_0` field are both set; leaving only one set is a
      common first-run error.

**Reading output**
- [ ] Check the evidence tier of each reported value.
- [ ] Look at the manifest diagnostics when a value looks
      surprising (UI: expand the "Provenance" chip).
- [ ] Export the ProcessDossier before closing the session.

## Appendix D. Frequently Asked Questions

**D1. I see a warning "trust blockers" — what does that mean?**
The simulator ran, but at least one result is below the trust
threshold for the scientific mode you chose. Options: drop to a less
strict scientific mode (Empirical), or accept the result at reduced
confidence. The UI shows the specific model name that blocked.

**D2. My d32 doesn't change when I change RPM.**
Two likely causes: (a) `t_emulsification` is too short — the DSD
hasn't converged (look for `converged = False` in the L1 output);
(b) you're at a breakage-coalescence equilibrium in which the two
rates cancel. Try raising RPM by 2×, or extending time to 10 min.

**D3. Can I simulate a stirred-tank bioreactor with cells?**
The emulsification physics is the same; the cell viability is not
modelled. The simulator will tell you what droplet size to expect;
you must check cell-compatibility of the process separately.

**D4. What does "DSD converged" mean?**
The population balance equation has reached a steady state between
breakage and coalescence. A non-converged DSD is typically too
broad and a poor predictor of the final product. Increase
`t_emulsification`.

**D5. Why does the alginate pipeline have no L3 output?**
Alginate's ionic gelation **is** its crosslinking. L3 is
intentionally stubbed for alginate; the modulus comes entirely
from L4. This is by design, not a bug.

**D6. Why is EDC/NHS giving a warning and a zero modulus?**
EDC/NHS requires surface COOH groups that pure chitosan does not
have. See §8 and Appendix H.5. Either pre-modify the matrix or
pick a different crosslinker.

**D7. How do I express uncertainty in my recipe inputs?**
Use `python -m emulsim uncertainty config.toml --n-samples 50`.
This runs Monte Carlo over the default perturbation spec; add
custom sources via the "Advanced" expander in the UI's
Uncertainty panel.

**D8. How is this different from a normal CFD simulation?**
CFD would solve Navier-Stokes at the droplet scale for hours to
days. EmulSim uses a **population balance equation** (PBE) at the
*droplet-size-distribution* scale, which is orders of magnitude
cheaper and gives the information you actually need (DSD, mean
diameters, span). Trade-off: we lose the per-droplet flow field.

**D9. Can I run this offline, without the UI?**
Yes. The CLI (`python -m emulsim run ...`) runs without a browser.
The UI is a convenience wrapper.

**D10. What is a ProcessDossier and why should I export it?**
A reproducibility artefact: a JSON bundle of every input, every
model manifest, every calibration applied, and every output. If
anyone (auditor, regulator, co-author) asks "how did you get this
number?", the ProcessDossier is your answer.

**D11. I want to add my own polymer chemistry. Is that possible?**
Yes, but it is an engineering task, not a UI task. See
`docs/f1a_alginate_protocol.md` and `docs/f1c_plga_protocol.md`
for the pattern: write an L2 solver, an L4 modulus, and a defaults
module, then wire them into the orchestrator.

**D12. Is my simulation reproducible?**
Yes, exactly reproducible given the same inputs and seed. The
random seed for Monte Carlo sampling is logged in the
ProcessDossier; re-running with that seed recovers the identical
output.

**D13. How big a bead can the simulator handle?**
The PBE tracks droplets from 100 nm to 1 mm. For beads > 1 mm,
use external CaCl₂ modes that are anyway process-limited by Ca²⁺
diffusion; the `solve_ionic_ca_gelation` handles up to ~5 mm with
no loss of accuracy (just longer solve time).

**D14. What's the fastest way to sweep a parameter?**
Use `python -m emulsim sweep` (L1 only) or the `design` subcommand
with a TargetSpec (full pipeline BO). Avoid naive loops over the
UI; the CLI is 10× faster.

**D15. Where do I report bugs?**
GitHub: `tocvicmeng-prog/Emulsification-Simulator`.

## Appendix E. Architectural Ideas & Working Principles

### E.1 The "levels" abstraction

EmulSim decomposes microsphere fabrication into four sequential
stages (**L1–L4**) plus two horizontal modules (**M1** export
contract, **M2** functionalisation, **M3** characterisation). Each
level has:

- A **solver** — a Python function that takes parameters and
  returns a result object.
- A **result schema** — a dataclass with typed, unit-carrying
  fields.
- A **model manifest** — a provenance record: which model was used,
  which assumptions it made, which evidence tier it carries, which
  diagnostics it produced.
- A **UI tab** (Modules 1–3, not directly Levels).

The **orchestrator** (`PipelineOrchestrator.run_single`) stitches
the levels together. It selects L2 and L4 based on the polymer
family, skips levels that don't apply (e.g. L3 for alginate /
cellulose / PLGA), and assembles the final `FullResult`.

### E.2 Multi-platform dispatch

All four polymer families share L1 (emulsification) — because
droplet formation is independent of chemistry. They diverge at L2:

```
L1 (PBE, universal)
    │
    ├─ if polymer_family == AGAROSE_CHITOSAN: L2-thermal-gelation
    ├─ if polymer_family == ALGINATE:          L2-ionic-Ca
    ├─ if polymer_family == CELLULOSE:         L2-NIPS
    └─ if polymer_family == PLGA:              L2-solvent-evaporation
        │
        ├─ (agarose-chitosan only): L3-crosslinking
        │
        L4 (family-specific modulus)
```

### E.3 Evidence & calibration

Every model carries an **evidence tier**. The RunReport computes the
minimum tier across the full pipeline and uses it to gate warnings.
A **CalibrationStore** lets wet-lab measurements promote tiers:
drop an AssayRecord JSON into `data/validation/<level>/assays/`,
run the fitter, and the next simulator run consumes the posterior
via `CalibrationStore.apply_to_model_params(target_module=...)`.

### E.4 Why a population balance equation, not CFD?

A full CFD simulation of 10⁷ droplets tumbling in a stirred
vessel is computationally prohibitive. The population balance
equation (PBE, §F.2) instead tracks the *number density function*
`n(v, t)` — how many droplets of each volume exist at each time.
Breakage and coalescence are represented as rates operating on
this distribution. The trade-off is that we lose information about
individual droplets (positions, instantaneous velocities) but we
retain everything we actually need: the DSD.

### E.5 Uncertainty propagation

`UnifiedUncertaintyEngine.run_m1l4` draws Monte Carlo samples from
user-specified perturbation distributions on material properties,
kernel coefficients, and formulation variables, pushes each sample
through the full pipeline in parallel, and reports the mean, std,
and raw samples of every observable. Calibration posteriors
(from `CalibrationStore`) are automatically included when
present.

### E.6 Inverse design

`emulsim design` takes a **TargetSpec** (desired d32, pore, G_DN,
K_av with tolerances) and runs Bayesian optimisation over the
process-parameter space, using a Gaussian-process surrogate and
Expected-Hypervolume-Improvement acquisition. Two robustness modes
are available: **mean-variance** (`--robust-variance-weight λ`)
penalises candidates whose resample std is large; **CVaR**
(`--robust-cvar-alpha α`) penalises candidates whose worst α-tail
is poor.

### E.7 Digital twin

`emulsim.digital_twin.run_replay` is an offline Ensemble Kalman
Filter harness. Given a recorded process trajectory (d32 measured
at several times, temperature log, etc.) and an ensemble of model
runs at perturbed parameters, it sequentially updates the ensemble
to estimate the hidden true parameters. Phase 1 is replay-only;
online filtering and model-predictive control are Phase 2.

## Appendix F. Chemical & Physical Principles

### F.1 Fluid mechanics of emulsification

Droplet breakup in a stirred vessel is driven by turbulent eddies
whose kinetic energy exceeds the interfacial-tension restoring force
of the droplet. The **Weber number**

```
We = ρ_c · ε^(2/3) · d^(5/3) / σ
```

compares inertial to capillary forces; breakage is active above a
critical `We ≈ 1.2` (Hinze 1955). The energy dissipation rate `ε`
in a stirred vessel scales with `N_p · n³ · D⁵ / V` where `N_p`
is the impeller power number, `n` is the rotational speed, `D` is
the impeller diameter, `V` is the vessel volume. This is what
makes d32 scale approximately as RPM⁻¹·² (the exact exponent
depends on the turbulence regime and the Weber-to-Kolmogorov
crossover, which the PBE resolves).

### F.2 Population balance equation

The discrete PBE in use (Hounslow 1988) tracks the number of
droplets `N_i` in `N_bins` geometric-spaced bins:

```
dN_i/dt = Σ [breakage_in − breakage_out + coalescence_in − coalescence_out]
```

- `breakage(v) = g(v)` — per-droplet breakage rate.
- `coalescence(v, v') = β(v, v') · h(v, v') · N(v) · N(v')` —
  second-order rate combining collision frequency `β` and
  coalescence efficiency `h`.

Kernels are parameterised; Node 30 unified UQ engine absorbs
calibration posteriors on the kernel coefficients automatically.

### F.3 Thermal gelation (agarose)

Agarose chains are random coils in solution above ~65 °C. On
cooling below the **gelation temperature** `T_gel ≈ 35–45 °C`,
helices form and aggregate into fibrils that cross-link into a
gel. The Avrami model captures the isothermal kinetics:

```
α(t) = 1 − exp[−(k·t)^n]
```

where `α` is the gel fraction, `k` is a rate constant, and `n` is
the Avrami exponent (3 for 3D heterogeneous nucleation).

### F.4 Ionic gelation (alginate)

Alginate is a copolymer of mannuronate (M) and guluronate (G)
residues. Ca²⁺ binds two G residues from adjacent chains in an
"egg-box" junction zone (Braccini & Pérez 2001). The L2 ionic-Ca
solver models the coupled mass transport:

```
dC/dt = D·∇²C − 2·k·C·G        (Ca²⁺ diffusion + binding)
dG/dt = −k·C·G                   (guluronate consumption)
dX/dt = +½·k·C·G                 (egg-box crosslink density)
```

### F.5 Non-solvent-induced phase separation (NIPS)

When a polymer solution contacts a non-solvent, the composition
crosses the binodal (metastable demixing) or spinodal (spontaneous
decomposition) of the ternary phase diagram. The **Cahn-Hilliard**
equation with a Flory-Huggins free energy captures the evolution:

```
∂φ/∂t = ∇ · (M · ∇(∂f/∂φ − κ ∇²φ))
```

`κ` is the gradient-energy coefficient (≈ 6·σ·ξ / ρ where σ is
interfacial tension and ξ is the interface width).

### F.6 Solvent evaporation (PLGA)

DCM diffuses out of the dispersed-phase droplet into the aqueous
continuous phase (perfect sink at droplet surface). The PLGA
volume fraction rises until it crosses the vitrification
threshold (`φ_PLGA ≈ 0.85`), at which point the local polymer
freezes into a glassy microsphere. The fixed-radius Fickian
approximation captures the kinetics adequately for `φ_PLGA < 0.8`.

### F.7 Rubber elasticity (covalent networks)

For a Gaussian-chain network with crosslink density `ν`,

```
G = ν · kT
```

(affine model, Flory-Rehner). For a double network, additive
superposition is a first approximation; non-trivial IPN coupling
coefficients (Node 10) refine this.

### F.8 Glassy polymer & Gibson-Ashby

Below its glass transition `T_g`, an amorphous polymer has a
modulus of ~1 GPa largely independent of crosslink density. For
a porous microsphere with polymer volume fraction `φ`,
Gibson-Ashby closed-cell scaling gives

```
G = G_glassy · φ^n
```

with `n ≈ 2` for closed-cell foams.

### F.9 Ogston sieving

The protein partition coefficient `K_av` between the pore space
and the bulk depends on the protein radius `r_h` and the fibre
volume fraction `φ_f`:

```
K_av = exp(−π · (r_h + r_f)² · L_f)       (Ogston 1958)
```

where `r_f` is the fibre radius and `L_f` is the fibre length per
volume. This controls chromatographic size exclusion.

## Appendix G. Formulas & Mathematical Theorems

### G.1 Kolmogorov microscale and Hinze limit

```
η_K = (ν³ / ε)^(1/4)            (Kolmogorov length)
d_max = C · ε^(−2/5) · σ^(3/5) · ρ^(−3/5)   (Hinze 1955)
```

### G.2 Discrete PBE (Hounslow 1988)

For geometric bin ratio 2:

```
dN_i/dt  = − g_i · N_i
           + Σ_j β_{ij} · S_{ij} · N_j · N_{i+1}
           − Σ_j (1 − β_{ij}) · C · N_i · N_j
           + (coalescence gain/loss terms)
```

Implementation: `src/emulsim/level1_emulsification/solver.py`.

### G.3 Avrami equation

```
α(t) = 1 − exp[−(k · t)^n]
```

### G.4 Kong 2004 alginate modulus

```
G_DN = K_alg · (c_alg · f_G)^n_alg · (X_mean / X_max)
X_max = ½ · c_alg · f_G / M_repeat
```

### G.5 Zhang 2020 cellulose modulus

```
G = K_cell · φ_cellulose^α_cell
```

with `K_cell ≈ 5 × 10⁵` Pa, `α_cell = 2.25` for NaOH/urea.

### G.6 Gibson-Ashby PLGA modulus

```
G = G_glassy · φ_PLGA^n_plga
```

with `G_glassy = 7 × 10⁸` Pa for PLGA 50:50, `n_plga = 2`.

### G.7 Ternary Flory-Huggins free energy

```
f/kT = (φ/N_p)·ln φ + s·ln s + n·ln n
     + χ_PS·φ·s + χ_PN·φ·n + χ_SN·s·n
```

where `n = 1 − φ − s` (volume conservation).

### G.8 Stochastic EnKF update (Evensen 1994)

```
K    = P_xy / (P_yy + R)
x_i^+ = x_i + K · (y_obs + ε_i − h(x_i))
```

with `ε_i ~ N(0, R)`. Implementation:
`src/emulsim/digital_twin/enkf.py`.

### G.9 Expected Hypervolume Improvement (EHVI)

The acquisition function for multi-objective BO (see
`optimization/engine.py`). For a reference point `r` and a
Pareto front `P`,

```
EHVI(x) = E_y[max(0, HV(P ∪ {y}) − HV(P))]
```

Approximated by MC integration over the GP posterior in our
implementation.

### G.10 CVaR aggregation (F4-b)

For a set of `N` resample objectives and tail fraction `α`,

```
CVaR_α = (1 / ⌈α·N⌉) · Σ_{i ∈ top_α(samples)} obj_i
```

Recovers `mean` at `α = 1` and approaches `max` as `α → 0`.

## Appendix H. Standard Wet-Lab Protocols

These are minimum-viable bench protocols that the simulator is
parameterised against. Each is a starting point; actual
laboratories will optimise around their local equipment and the
simulator's predictions.

### H.1 Agarose-chitosan double-network microspheres (genipin)

**Materials:** agarose (Type 1B), chitosan (medium M_w), paraffin
oil, Span-80, genipin, citric acid buffer pH 4.5.

1. Dissolve 4.2 g agarose in 100 mL buffer at 90 °C. In parallel,
   dissolve 1.8 g chitosan in 100 mL buffer at 60 °C. Combine
   above 80 °C to give the aqueous polymer phase.
2. Mix 200 mL paraffin oil with 4 g Span-80; heat to 90 °C in the
   emulsification vessel.
3. Add the aqueous phase to the oil under agitation (rotor-stator
   at the RPM recommended by the simulator) for the predicted
   `t_emulsification`.
4. Cool to 4 °C at ~10 °C/min; agarose gels below 45 °C.
5. Separate the beads by centrifugation; wash three times with
   warm (40 °C) buffer to remove residual oil and surfactant.
6. Resuspend in a 2 mM genipin solution in buffer at 37 °C;
   incubate 24 h with gentle stirring.
7. Wash thoroughly; store in buffer at 4 °C.

### H.2 Ca²⁺-alginate microspheres (external bath)

**Materials:** sodium alginate (medium G-content), paraffin oil,
Span-80, 100 mM CaCl₂ aqueous bath.

1. Dissolve 2 g sodium alginate in 100 mL deionised water.
2. Emulsify into 200 mL paraffin/Span-80 at room temperature under
   the simulator-recommended RPM for the predicted time.
3. Introduce the CaCl₂ bath by slow addition of 50 mL 100 mM
   CaCl₂ (aqueous) to the stirred emulsion. Alternatively, break
   the emulsion by centrifugation and drop the aqueous droplets
   directly into the Ca²⁺ bath.
4. Allow 30 min for gelation at room temperature (simulator
   predicts completion time based on bead radius).
5. Collect beads; wash three times with deionised water to remove
   excess Ca²⁺.
6. Store in 10 mM CaCl₂ at 4 °C.

### H.3 Ca²⁺-alginate (internal GDL/CaCO₃, homogeneous)

**Materials:** sodium alginate, CaCO₃ (<1 µm), glucono-δ-lactone
(GDL), paraffin oil, Span-80.

1. Dissolve 2 g sodium alginate in 100 mL water. Disperse 0.2 g
   CaCO₃ homogeneously; ultrasonicate if needed.
2. Immediately before emulsification, dissolve 0.7 g GDL in the
   aqueous phase (2:1 molar ratio GDL:Ca²⁺). Emulsify immediately
   — the clock is now running.
3. Emulsify for the simulator-recommended time at RT.
4. Maintain gentle agitation for 4 h while GDL hydrolyses and
   Ca²⁺ releases in situ.
5. Collect beads; wash with water; store at 4 °C.

### H.4 Cellulose NaOH/urea NIPS microspheres

**Materials:** cellulose (DP ~370, e.g., cotton linter),
7 wt% NaOH + 12 wt% urea aqueous solution (pre-chilled to
−12 °C), paraffin oil, Span-85 (more lipophilic for cold
emulsification), deionised water (non-solvent bath).

1. Dissolve 5 g cellulose in 100 mL pre-chilled NaOH/urea by one
   freeze-thaw cycle (−12 °C, 12 h; thaw with stirring at 5 °C).
2. Emulsify the cellulose solution into 200 mL paraffin/Span-85
   at 10 °C under the simulator-recommended RPM.
3. Slowly add 200 mL water to the emulsion (or break and drop
   into water bath); cellulose precipitates as demixing proceeds.
4. Wait 30–60 min at RT; wash the regenerated cellulose beads
   thoroughly with water to remove residual NaOH/urea.
5. Store in water at 4 °C, or freeze-dry for dry product.

### H.5 EDC/NHS activation and coupling (requires COOH)

**Pre-requisite:** the microsphere surface must carry COOH groups.
For chitosan, use carboxymethyl chitosan (CMC) or
succinyl-chitosan as the starting polymer. For agarose, pre-
oxidise with succinic anhydride.

**Materials:** MES buffer pH 5.5, EDC (1-ethyl-3-[3-
dimethylaminopropyl]carbodiimide), NHS (N-hydroxysuccinimide),
the ligand amine.

1. Suspend beads (~5 % w/v) in MES pH 5.5 at 4 °C.
2. Add EDC to 50 mM and NHS to 25 mM. Activate 30 min at 4 °C
   with gentle agitation.
3. Add ligand amine to 10 mM. React 2–4 h at 4 °C.
4. Wash beads in PBS to quench.
5. Store at 4 °C in PBS + 0.05 % sodium azide.

Simulator note: feed `surface_cooh_concentration` (mol/m³ of
grafted COOH) as a MaterialProperties override, not as a
formulation concentration.

### H.6 PLGA solvent-evaporation microspheres (O/W)

**Materials:** PLGA 50:50 (MW 30 kDa), dichloromethane (DCM),
1 % (w/v) PVA aqueous solution.

1. Dissolve 1 g PLGA in 10 mL DCM. (Organic phase: 10 wt%
   PLGA/DCM.)
2. Emulsify into 100 mL 1 % PVA aqueous phase under the
   simulator-recommended RPM at 25 °C for the predicted time.
3. Transfer the emulsion to a larger vessel with gentle stirring,
   open to air. Allow DCM to evaporate (typically 2–4 h).
4. Collect hardened beads by centrifugation; wash three times
   with water.
5. Freeze-dry for storage.

Safety: DCM is volatile and has low acute toxicity but is a
suspected carcinogen. Use in a fume hood with proper ventilation.

## Appendix I. Troubleshooting Table

| Symptom | Likely Cause | Fix |
|---|---|---|
| `d32` too large | RPM too low, or `t_emulsification` too short | Increase RPM 1.5×; extend time to 10 min |
| `d32` too small | RPM too high, or surfactant over-concentrated | Lower RPM; reduce surfactant to 1 % |
| DSD did not converge | `t_emulsification` too short | Extend to 10–15 min; check `converged` flag |
| Modulus zero | Crosslinker concentration 0, or wrong chemistry for the matrix | Re-verify crosslinker choice; for EDC/NHS confirm `surface_cooh_concentration > 0` |
| Pore size tiny (<10 nm) | Polymer concentration too high, or porosity model not applicable | Reduce polymer wt%; check platform-appropriate L2 model |
| Pore size huge (>10 µm) | Demixing too slow (cellulose) or incomplete gelation | Deeper quench; raise Ca²⁺ or crosslinker; extend gel time |
| Trust blockers warning | Model used is below threshold for current scientific mode | Drop scientific mode to `hybrid` or `empirical`; or upgrade with wet-lab calibration |
| `p_final = 0` for alginate | Expected behaviour — no L3 | Ignore; modulus comes from L2/L4 |
| EDC/NHS warning + QUALITATIVE_TREND | No surface COOH | Pre-modify matrix (CMC, succinylation) |
| Cellulose run times > 60 s | Cahn-Hilliard stiffness | Reduce `n_r` to 30; verify `phi_cellulose_0` is not out of range |
| PLGA `t_vitrification` unreported | Mean φ_PLGA never crossed 0.85 in the simulated time | Extend `t_crosslink`; confirm Dirichlet sink is effective |
| Calibration not applied | AssayRecord fit JSON missing / wrong target_module | Re-run `emulsim ingest`; check CalibrationStore loader logs |
| Reproducibility failure | RNG seed changed between runs | Pin `--seed` explicitly; check ProcessDossier for the used seed |
| Out-of-range parameter warning | Input outside model validation envelope | Re-check input against Appendix A.2/A.3 ranges; narrow the recipe or accept an elevated-uncertainty result |
| "Polymer family mismatch" error | Two competing polymer concentrations (e.g. c_agarose and phi_PLGA_0 both > 0) | Zero-out the unused family |
| UI tab blank | Session state corrupted | Refresh the browser; restart `python -m emulsim ui` |

---

**End of First Edition.**

*For updates, corrections, and contributions, visit the repository
at* `github.com/tocvicmeng-prog/Emulsification-Simulator`.
