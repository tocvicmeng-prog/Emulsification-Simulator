# Node 32 — Cluster F v8.0 Roadmap

**Prepared by:** /architect (via dev-orchestrator)
**Date:** 2026-04-17
**Status:** Roadmap document; scopes v8.0 workstreams. No code in this node.
**Supersedes / refines:** `docs/10_future_development_roadmap.md` §4

---

## 1. Purpose and scope

v8.0 is the "platform expansion + feedback-loop" release: it takes the
mature chitosan/agarose pipeline shipped in v7.x and turns EmulSim into
(a) a cross-platform microsphere simulator and (b) a closed-loop
design/optimisation/assimilation engine. This document scopes the five
workstreams Doc 10 §4 grouped as "Cluster F" into concrete Node-level
deliverables, identifies dependencies and entry criteria, and proposes
an execution order.

The roadmap is intentionally conservative on scope per item and
aggressive on dependency clarity — v8.0 will span multiple sessions and
the framework's value is in knowing what unblocks what before any code
is written.

## 2. Workstream inventory (one section per Cluster F item)

### 2.1 F1 — Other Microsphere Platforms

**Doc 10 §4.1.** Extend L1–L4 to platforms beyond chitosan/agarose.

**Scientific substance by platform:**

| Platform | Gelation mechanism | L2 solver reuse | New infrastructure |
|---|---|---|---|
| **Alginate (Ca²⁺)** | Diffusion-limited ionic gelation; shrinking-core | L2 needs new ionic-CH-like solver (concentration-dependent crosslink front) | CaCl₂ bath concentration, ion transport coefficients, Ca²⁺–COO⁻ binding isotherm |
| **PLGA (W/O/W + solvent extraction)** | Solvent–non-solvent phase inversion; not thermal TIPS | Fundamentally new L2: solvent extraction kinetics, polymer precipitation | Solvent partitioning coefficients, Hansen solubility parameters, extraction-rate kinetics |
| **Cellulose (NIPS)** | Non-solvent-induced phase separation; ternary CH applicable | Existing CH solver with ternary solvent/non-solvent/polymer composition | Three-component Flory-Huggins parameters, bath concentration dynamics |

**Complexity:** HIGH. Alginate is the closest migration (~4 weeks, Opus
tier architecture + Sonnet impl). PLGA is a ground-up new pipeline
(~8–12 weeks). Cellulose reuses much of the existing CH machinery
(~4–6 weeks).

**Entry criteria:**
- Scientific advisor review (mechanism + parameter tables per
  platform) — parallel to this document.
- Decision on which platform ships first. Recommendation: **alginate**
  (closest to current framework, smallest code delta, immediate
  commercial demand in chromatography).

**Dependencies:**
- Independent of F2–F5; can start immediately.
- Calls /scientific-advisor + /r&d-researcher for per-platform parameter tables.

**Footprint estimate (alginate only):** ~1500–2000 LOC across new
`level1_emulsification/alginate.py`, new `level2_gelation/ionic_ca.py`,
reagent library entries, validation data scaffold, and ~30 new tests.

**Nodes:** F1-a alginate (3–5 nodes), F1-b cellulose NIPS (2–3 nodes),
F1-c PLGA (5–8 nodes), each containing architect protocol + code +
tests + docs.

---

### 2.2 F2 — Digital Twin Capability

**Doc 10 §4.2.** Real-time process data assimilation via Ensemble
Kalman Filter (EnKF); online Bayesian parameter estimation; predictive
RPM / temperature / cooling-rate control.

**Scientific substance:**

- **State estimation:** EnKF over an ensemble of pipeline runs, each
  with its own MaterialProperties + KernelConfig perturbation. Update
  each ensemble member's state when live process measurements (d32,
  torque, turbidity) arrive.
- **Parameter estimation:** online Bayesian updates to the
  `CalibrationStore` posteriors as measurements accumulate. The Node
  30 posterior propagation machinery already supports this at a
  per-batch level; the digital twin turns that into a streaming loop.
- **Control handoff:** a lightweight MPC layer that chooses the next
  RPM / cooling-rate step to drive the predicted d32 trajectory toward
  the target specification.

**Complexity:** HIGH (novel control + estimation coupling). Requires
real or simulated live sensor data to test; without hardware-in-loop
or a sensor emulator, the value is limited to algorithmic
demonstration.

**Entry criteria:**
- F1 complete (or restricted to chitosan/agarose only).
- A defined sensor interface (Node 32b? or deferred to v8.1?) — either
  a live OPC-UA / MQTT / REST adapter or a replay harness.
- Scientific advisor + full-stack-engineer consultation on EnKF
  choice vs particle filter, ensemble size vs compute budget.

**Dependencies:**
- Strong dependency on **Node 30 unified UQ engine** — EnKF ensemble
  IS the UQ MC reframed as a time-series. Good reuse.
- Moderate dependency on F1 if applied to platforms other than
  chitosan/agarose.
- Independent of F3, F4, F5.

**Footprint estimate:** ~2500 LOC (EnKF core, sensor adapter, MPC
layer, streaming harness, 20+ tests). Opus tier architecture; Opus
implementation for the EnKF (novel numerics); Sonnet for sensor
adapter boilerplate.

**Nodes:** F2-a sensor harness + replay (2 nodes), F2-b EnKF state
estimation (3 nodes), F2-c online parameter update (2 nodes), F2-d MPC
control loop (3 nodes). Total ~10 nodes.

---

### 2.3 F3 — Inverse Design

**Doc 10 §4.3.** Given target specs (d32, pore, G_DN, Kav), find the
optimal formulation via constrained Bayesian Optimisation.

**Scientific substance:**

- Standard BO framework (GPyTorch / BoTorch) wrapping the existing
  pipeline.
- Objective: weighted distance between predicted (d32, pore, G_DN,
  Kav) and user-specified targets, under physical constraints
  (agarose–chitosan ratio, phi_d limits, feasible RPM range).
- Acquisition: expected hypervolume improvement for multi-objective
  target matching; expected improvement for single-objective.

**Complexity:** MEDIUM. The optimisation infrastructure is standard;
the hard work is defining physically meaningful constraints and a
well-conditioned objective that doesn't reward degenerate solutions
(e.g. zero-crosslinker giving a low G_DN trivially matching a "soft"
target).

**Entry criteria:**
- Existing optimiser (`pipeline.optimizer`) reviewed for reuse vs
  greenfield. Recommendation: reuse the trust-aware evidence gate
  from Node 6; wrap a new acquisition function around it.
- /scientific-advisor review of the constraint set to avoid physically
  impossible targets.

**Dependencies:**
- Depends on mature UQ (Node 30) for uncertainty-aware acquisition.
- Depends on F1 only if used for alt platforms.
- Independent of F2, F4, F5 conceptually; F4 is a natural extension.

**Footprint estimate:** ~800–1200 LOC. Opus architecture; Sonnet
implementation; Sonnet tests.

**Nodes:** F3-a constraint-aware acquisition (1 node), F3-b
multi-objective target matcher (1 node), F3-c CLI + UI surface (1
node). Total ~3 nodes.

---

### 2.4 F4 — Robust Optimisation Under Uncertainty

**Doc 10 §4.4.** Optimise E[f(x + δ)] rather than f(x) — find
formulations whose performance is insensitive to manufacturing
perturbations.

**Scientific substance:**

- For each candidate x in the BO loop, run a small MC ensemble over
  perturbations δ (input-uncertainty draw from the unified spec's
  MaterialProperty or custom sources).
- Acquisition on E[f] (mean of the objective) plus optional penalty
  on Var[f] for a "mean-variance" trade-off.
- Batch-variability framework from Node 19 already does the per-x
  ensemble work; F4 layers the BO decision on top.

**Complexity:** MEDIUM (builds on F3 + Node 19 machinery).

**Entry criteria:**
- F3 scaffold exists (F4 is an acquisition-function variant of F3).
- Agreement on the risk metric (E[f], E[f] + λ·Var[f], or CVaR).

**Dependencies:**
- Hard dependency on F3.
- Hard dependency on Node 19 (batch variability — already shipped).
- Hard dependency on Node 30 (unified UQ — already shipped).
- Independent of F1, F2, F5.

**Footprint estimate:** ~400–600 LOC incremental over F3. Sonnet
implementation.

**Nodes:** F4-a mean-variance acquisition (1 node), F4-b CVaR
acquisition (1 node). Total ~2 nodes.

---

### 2.5 F5 — MD Parameter Estimation

**Doc 10 §4.5.** Use MARTINI coarse-grained MD for `chi_ac` (Flory-Huggins
χ for agarose–chitosan), `kappa_CH` (gradient energy coefficient),
`M_0` (bare mobility), `f_bridge` (bridge efficiency).

**Scientific substance:**

- MARTINI mapping for agarose and chitosan (literature reference:
  Lopez et al. 2015 JCTC 11:2714 for polysaccharide CG params).
- Equilibrium MD on binary polymer mixtures to extract χ from
  radial distribution functions (Kirkwood-Buff integrals).
- Non-equilibrium MD for mobility `M_0` (forced-gradient method).
- Pulling MD for `f_bridge` (probability of crosslink bridging two
  chains vs looping back).
- Pipeline: MD run → postprocessing → `CalibrationStore` JSON →
  absorbed by Node 30 unified UQ engine.

**Complexity:** HIGH. MD infrastructure is external (GROMACS /
LAMMPS); the EmulSim-side work is a data-ingest pipeline, not an
MD engine. But the MD setup, equilibration, and analysis protocols
are non-trivial and need dedicated /scientific-coder + domain
consultation.

**Entry criteria:**
- Access to MD software + compute (GPU cluster or at least a
  dedicated workstation).
- /scientific-advisor consultation on MARTINI mapping validation for
  chitosan specifically (CG mappings for GlcN are published but vary).
- Infrastructure decision: does EmulSim launch MD jobs, or does it
  only ingest MD output JSONs produced externally?

**Dependencies:**
- Soft dependency on F1 (for platforms where χ is not already
  parameterised — alginate needs its own MD workstream).
- Strong dependency on Node 30 (CalibrationStore posterior absorption
  is the MD-to-pipeline handoff).
- Independent of F2, F3, F4.

**Footprint estimate (ingest-only):** ~600–800 LOC (MD-output parsers
for GROMACS/LAMMPS, RDF/Kirkwood-Buff math, `CalibrationEntry`
builders, ~15 tests). If we launch MD jobs: +1500 LOC (job dispatcher,
checkpoint-restart, result provenance).

**Nodes:** F5-a MD-output ingest for χ (2 nodes), F5-b mobility
extraction (1 node), F5-c f_bridge pulling analysis (1 node), F5-d
(optional) MD job dispatcher (3 nodes). Total 4–7 nodes.

---

## 3. Dependency graph

```
                 ┌──────────────────────┐
                 │ Node 30 UQ (shipped) │
                 └─────────┬────────────┘
                           │
   ┌──────────────┬────────┴────────┬──────────────┬──────────────┐
   │              │                 │              │              │
   ▼              ▼                 ▼              ▼              ▼
 [F1]          [F2]              [F3]           [F4 needs F3]   [F5]
Platforms   Digital Twin     Inverse Design     Robust Opt.     MD params
   │              │                 │              │              │
   │              │                 │              ▼              │
   └──────────────┴────┐        ┌───┴──────────────┘              │
                       │        │                                 │
                       ▼        ▼                                 │
                   [F2 benefits]                                  │
                   (EnKF uses UQ                                  │
                    machinery)                                    │
                                                                  │
   F5 feeds all (MD-calibrated constants -> MaterialProperties)   │
                                                                  │
   F1-alginate benefits from F5 (MARTINI for Ca-alginate binding) ┘
```

Critical observation: **F3 and F4 can be shipped before F1/F2/F5** —
they only need the current mature pipeline + unified UQ. Shipping them
first demonstrates commercial value (inverse design of chitosan/agarose
microspheres) without the long-pole platform work.

## 4. Proposed execution order

### v8.0 release (12–16 weeks)

**Phase 1 (weeks 1–4): Inverse Design foundation**
- F3-a, F3-b, F3-c — inverse design end-to-end for chitosan/agarose
  (current pipeline only). Demonstrates commercial value early.
- F4-a, F4-b — robust optimisation acquisition variants stacked on F3.
- **Exit:** v8.0-alpha ships with `python -m emulsim design --target …`
  CLI and streamlit UI.

**Phase 2 (weeks 5–10): First alt platform**
- F1-a alginate complete (diffusion-limited ionic gelation L2, reagent
  library entries, validation scaffold, ~5 nodes).
- F3/F4 re-validated on alginate to prove cross-platform inverse
  design works.
- **Exit:** v8.0-beta ships with alginate + chitosan/agarose platforms.

**Phase 3 (weeks 11–16): MD ingest + Digital Twin scaffold**
- F5-a, F5-b, F5-c ingest-only MD pipeline.
- F2-a sensor-replay harness (sensor adapter deferred to v8.1 if no
  hardware).
- **Exit:** v8.0.0 GA.

### v8.1 and beyond (post-v8.0)

- F1-b cellulose NIPS
- F1-c PLGA (if commercial demand)
- F2-b, F2-c, F2-d — full digital twin with live sensor integration
- F5-d MD job dispatcher (only if EmulSim grows MD infrastructure)

## 5. Entry criteria for v8.0 kickoff

Before any v8.0 node is implemented:

1. **v7.0 must ship** — current blocker is Study A wet-lab data for
   Node 21 F1 closure. v8.0 should not start until v7.0 releases so
   we don't stack validation debt.
2. **Scientific advisor reviews the roadmap** — particularly F1
   platform mechanisms, F2 EnKF choice, F5 MARTINI mappings.
3. **Chief economist** assesses NPV of each Cluster F workstream
   against wet-lab validation cost and time-to-market. Commercial
   prioritisation may reorder the phases above.
4. **ip-auditor** sweeps for F1/F3/F4 IP overlaps — inverse design
   and robust optimisation are patent-dense areas; FTO review before
   filing any commercial claims.

## 6. Risks and open questions

| # | Risk | Severity | Mitigation |
|---|---|---|---|
| R1 | Study A wet-lab data doesn't calibrate chitosan/agarose tightly enough for inverse design to find non-degenerate optima | HIGH | Ship F3 with UQ-aware acquisition (robust-to-calibration BO); F3 outputs flagged with evidence tier |
| R2 | Alginate ionic gelation diverges too much from current L2 — requires full L2 rewrite, not extension | MED | F1-a architecture review before code; prototype ionic solver as isolated module |
| R3 | MD parameter extraction requires compute budget the project doesn't have | MED | F5 ships ingest-only by default; external MD runs produce the JSONs |
| R4 | Digital twin (F2) needs hardware that doesn't exist yet | HIGH | Scope F2 to replay-only for v8.0; defer live integration to v8.1 |
| R5 | PLGA (F1-c) is a fundamentally different process (solvent extraction, not TIPS) — may require L1 rewrite too | MED | Defer to v8.1+; not in v8.0 scope |
| R6 | Inverse design rewards optima outside the model's validity domain | MED | F3 acquisition penalised by trust tier (from Node 6 trust-aware machinery); hard-clip at domain boundaries |

**Open questions for the user / project leadership:**

- Q1: Does the project have or plan to have MD compute for F5? If no,
  F5 scope shrinks to ingest-only.
- Q2: Is there a target hardware partner for F2 digital-twin? If no,
  F2 is algorithm-demo only for v8.0.
- Q3: Commercial priority: inverse design first (pull v8.0-alpha
  forward), or alginate first (pull v8.0-beta forward)?
- Q4: After Study A delivers, do we want to re-run calibration sweeps
  across F1 platforms, or assume each platform gets its own Study
  (A', A'', …)?

## 7. What this roadmap deliberately does NOT cover

- **v7.2 features.** Node 31b (dedicated c_edc / c_nhs fields),
  any further L3 chemistry additions, UI polish — all live outside
  Cluster F.
- **Internal refactors.** Anything that doesn't deliver new science
  or new workflow capability is v7.x maintenance, not v8.0.
- **Performance work** unless a Cluster F item specifically needs it
  (F2 digital twin probably does; nothing else does).
- **Infrastructure migrations.** Python version bumps, dependency
  pinning changes, CI overhauls — all separate from Cluster F.

## 8. Next actions for dev-orchestrator after Node 32

The dev-orchestrator's Node 32 ends with this document. The next
orchestrator session should:

1. Invite the user to review §6's open questions Q1–Q4.
2. Decide on the Phase 1 start item (F3 inverse design is the
   recommended default per §4).
3. Invoke /chief-economist for commercial prioritisation of F1
   platforms (alginate vs cellulose vs PLGA).
4. Invoke /scientific-advisor for per-workstream mechanism briefs
   when the first v8.0 node is scheduled.

**Suggested first v8.0 node:** F3-a — constraint-aware acquisition
wrapper around the existing `pipeline.optimizer`. Low-risk, demos the
inverse-design value proposition, and validates the v8.0 architecture
stance that UQ (Node 30) and trust-aware evidence (Node 6) are the
foundation we build on.

---

**End of Node 32 / Cluster F v8.0 roadmap.**
