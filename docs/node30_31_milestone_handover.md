# Milestone Handover: M30 — EmulSim v7.1-dev Node 30 complete, Node 31 ready for implementation

**Date:** 2026-04-17
**Session:** dev-orchestrator (Opus); originating session
`cd5306a2-97b9-4ce6-9bb9-27cc79361ba4`
**Project:** EmulSim — multi-scale hydrogel microsphere simulator
**Prepared by:** /architect (via dev-orchestrator)
**Classification:** Internal — Development Handover

---

## 1. Executive summary

Node 30 (Full UQ merge) shipped in this session: consolidated two
legacy MC engines into a single `UnifiedUncertaintyEngine`, closed
Audit N2 by actually propagating `CalibrationStore` posterior samples
into the MC (they were declared but never sampled in v7.0). Net
footprint −340 LOC, 52 targeted-regression tests passing, 0 new
regressions.

Node 31 (EDC/NHS mechanistic kinetic) completed its design phase this
session: a scientific brief from /scientific-advisor re-scoped the
work from "L3 solver" (scientifically inappropriate — EDC/NHS needs
COOH which native chitosan/agarose lacks) to "M2 mechanistic model
with gated L3 exposure". The /architect protocol passed Gate G1.
**Ready for /scientific-coder Phase 2 implementation** in a fresh
session.

Node 32 (Cluster F v8.0 roadmap) remains queued. Blocked on Node 31
completion.

## 2. Module registry

| # | Module / Node | Version | Status | Approved | Model Used | Fix Rounds | Lines | File Path |
|---|---|---|---|---|---|---|---|---|
| 30 | Full UQ merge | v7.1.0-dev | APPROVED | 2026-04-17 | Opus (design + impl) | 0 | −540 deleted / +200 added (net −340) | `src/emulsim/uncertainty_unified.py`, `tests/test_uncertainty_unified.py`, `tests/test_parallel_mc.py` |
| 31 | EDC/NHS kinetic | v7.1.0-dev | PROTOCOL-APPROVED | — | Opus (design) | — | ~350–500 projected | `docs/node31_edc_nhs_scientific_brief.md`, `docs/node31_protocol.md` |
| 32 | Cluster F roadmap | v8.0-scoping | QUEUED | — | Opus planned | — | doc only | — |

## 3. Integration status

| Interface | From Module | To Module | Status | Notes |
|---|---|---|---|---|
| `UnifiedUncertaintyResult` | `uncertainty_unified.py` | CLI, `ProcessDossier`, batch, optimizer | LIVE | `kinds_sampled` now honestly includes CALIBRATION_POSTERIOR when posteriors dispatched; `kinds_declared_but_not_sampled` tracks malformed/unknown posteriors |
| `CalibrationStore` posteriors | `calibration.calibration_store` | `UnifiedUncertaintyEngine` | LIVE | Audit N2 closed — posteriors actually perturb MC |
| `OutputUncertainty.raw_samples` | `uncertainty_unified.py` | Test harness, downstream analyses | LIVE | New field added to preserve bit-identicality tests |
| M2 `edc_nhs_activation` profile → mechanistic kinetic | `module2_functionalization` | — | PENDING (Node 31) | Currently routes to pseudo-single-step (ranking_only) |
| L3 `michaelis_menten` branch → gated mechanistic | `level3_crosslinking.solver` | — | PENDING (Node 31) | Currently falls back to second-order amine + QUALITATIVE_TREND |
| `MaterialProperties.surface_cooh_concentration` | `datatypes.py` | L3 gate | PENDING (Node 31) | New optional field default 0.0 |
| Streamlit uncertainty panel | `visualization.panels.uncertainty` | `visualization.app` | DEFERRED (Node 30b) | Shows st.info placeholder, returns None |

## 4. Current code inventory

### Node 30 (shipped)

| File | Purpose | Status | Lines |
|---|---|---|---|
| `src/emulsim/uncertainty_unified.py` | Unified MC engine, schema, sampler, posterior dispatch | APPROVED | 475 (was 420, +55 net after merge) |
| `src/emulsim/uncertainty_core.py` | DELETED | — | was 318 |
| `src/emulsim/uncertainty_propagation/` | DELETED (package) | — | was 216 |
| `src/emulsim/__main__.py` | `--engine {unified,legacy}` both route through merged engine | UPDATED | — |
| `src/emulsim/visualization/panels/uncertainty.py` | Placeholder (Node 30b) | UPDATED | 28 |
| `tests/test_uncertainty_unified.py` | Schema + engine tests incl Audit N2 closure | UPDATED | 261 |
| `tests/test_parallel_mc.py` | Parallel/serial bit-identicality on merged engine | UPDATED | 99 |
| `tests/test_cli_v7.py` | Legacy-flag test rewritten as unified-route test | UPDATED | — |
| `tests/test_v60_integration.py` | M2 UQ tests deleted (dead path) | UPDATED | — |
| `CHANGELOG.md` | v7.1.0-dev Node 30 entry prepended | UPDATED | — |

### Node 31 (protocol-approved, implementation pending)

| File | Purpose | Status | Lines |
|---|---|---|---|
| `docs/node31_edc_nhs_scientific_brief.md` | Mechanism, rate-constant table, viability | APPROVED (scientific-advisor) | 275 |
| `docs/node31_protocol.md` | Full module protocol; G1 passed | APPROVED (architect) | 340 |
| `src/emulsim/module2_functionalization/edc_nhs_kinetics.py` | New: kinetic core | PENDING | ~250 projected |
| `src/emulsim/module2_functionalization/reactions.py` | Dispatch branch for `chemistry_class="edc_nhs"` | PENDING | ~30 new |
| `src/emulsim/module2_functionalization/reagent_profiles.py` | Tier promotion for `edc_nhs_activation` | PENDING | ~5 edits |
| `src/emulsim/level3_crosslinking/solver.py` | Gated mechanistic branch replacing michaelis_menten fallback | PENDING | ~80 new / −50 deleted |
| `src/emulsim/datatypes.py` | New `MaterialProperties.surface_cooh_concentration` field | PENDING | ~3 new |
| `tests/test_edc_nhs_kinetics.py` | 10–12 tests per protocol §6 | PENDING | ~150 projected |

## 5. Architecture state

No high-level architectural changes since the Nodes 1–29 baseline. Two
notable boundary adjustments:

1. **UQ system reduced from 3-layer (two legacy engines + facade) to
   1-layer (merged engine)**. Public API (`UnifiedUncertaintyEngine`,
   `UnifiedUncertaintySpec`, `UnifiedUncertaintyResult`) unchanged;
   legacy `UncertaintyPropagator` and `run_with_uncertainty` entry
   points are gone.
2. **EDC/NHS scientific home clarified** (Node 31 design): the
   chemistry lives in M2 (functionalisation) with an optional L3
   surface-COOH gate. The L3-only framing assumed in the earlier
   roadmap (doc 10 §4) was a category error; Node 31 corrects it.

## 6. Design decisions log

| Decision | Rationale | Date | Alternatives considered |
|---|---|---|---|
| Node 30 Option B++ (merge uncertainty_core into uncertainty_unified + delete uncertainty_propagation) | M2 UQ path was unreachable from CLI; dead code. Merge preserves battle-tested `_mc_one_sample` worker; Audit N2 closure is the core value | 2026-04-17 | A (aggressive rewrite), Keep both engines |
| Panel migration deferred to Node 30b | Simulator-side merge is the scientifically important part; UI rebuild around UnifiedUncertaintySpec is mechanical and can land isolated | 2026-04-17 | Migrate in Node 30 (+150 LOC, streamlit visual QA), Delete panel |
| CLI `--engine legacy` retained as no-posterior synonym | Preserves byte-compat for scripts that want default MaterialProperties perturbations only; costs nothing | 2026-04-17 | Remove flag (breaks v7.0.x scripts) |
| Preserve exact RNG call order for default perturbations | Keeps numerical regression tests stable; posteriors draw AFTER default 10 to avoid perturbing the sequence | 2026-04-17 | Allow RNG reshuffle (cleaner code, breaks seed reproducibility) |
| Node 31 scope: M2 + L3 gated (not L3-only) | Scientifically correct — EDC on native matrix has nothing to activate; M2 is the chemistry's real home | 2026-04-17 | M2 only (loses L3 extensibility), L3 only (scientifically weak) |
| Node 31 tier target: SEMI_QUANTITATIVE with literature constants | Node 30 closed the calibration-posterior machinery so Study A can promote to QUANTITATIVE later without rework | 2026-04-17 | QUANTITATIVE blocked on Study A, SEMI+extra Arrhenius audit |

## 7. IP constraints encountered

No new IP constraints identified this session. EDC/NHS is an
80-year-old chemistry (patents long expired); all cited rate constants
are from peer-reviewed literature.

## 8. Open questions / unresolved issues

| # | Question | Awaiting | Priority |
|---|---|---|---|
| 1 | Pre-existing `tests/test_data_layer.py::TestKernelConfig::test_for_rotor_stator_legacy` failure (phi_d_correction True vs expected False) | Separate investigation — NOT caused by Node 30 | LOW |
| 2 | `FormulationParameters.pH` field may need adding for Node 31 | Verify in `datatypes.py` during Node 31 Phase 2 | MED |
| 3 | `AMIDE_BOND` ACS enum value may not exist | Verify in `acs.py` during Node 31 Phase 2; otherwise reuse NHS_ESTER transitions | MED |
| 4 | Node 30b panel migration (streamlit → UnifiedUncertaintySpec) | User prioritisation — post-Node 31 | LOW |
| 5 | Node 32 (Cluster F v8.0 roadmap) scope and depth | User direction after Node 31 ships | LOW |

## 9. Next module protocol

**See `docs/node31_protocol.md`** — complete, G1-passed, ready for
/scientific-coder. Estimated effort: 350–500 LOC, Sonnet tier,
1 fix round likely, ~25k tokens of implementation work.

Entry point for next session:

```
User: start Node 31 implementation
Claude: (reads docs/node31_protocol.md + docs/node31_edc_nhs_scientific_brief.md)
        (invokes /scientific-coder for Phase 2)
        ... implements 5 files ...
        ... runs tests, hands off ...
        (invokes /architect for Phase 3 audit)
```

## 10. Context compression summary

This session exited dev-orchestrator Phase 2 of Node 30 with Node 30
APPROVED, and completed Phase 0 + Phase 1 of Node 31 (scientific brief
+ protocol). Context budget at handover: estimated YELLOW zone (~40%
remaining). No active in-progress work; next session starts from a
clean slate.

**Approved work in this session (verbatim list):**
- Node 30 /architect Phase 0 (pre-flight, scope decisions)
- Node 30 /architect Phase 1 (protocol — compact form)
- Node 30 /scientific-coder Phase 2 (merged engine + CLI + panel guard + test updates)
- Node 30 /architect Phase 3 (audit — 52 tests pass, 0 regressions; APPROVED)
- Node 30 /architect Phase 5 (registration — CHANGELOG, memory, handover)
- Node 31 /scientific-advisor literature + mechanism brief
- Node 31 /architect Phase 1 (protocol — G1 passed)

**Compressed / dropped context:**
- Detailed grep/read output from file surveys (captured in scientific brief + protocol)
- Intermediate debugging of the 2 Node 30 test failures (fixed: lazy-create KernelConfig in posterior dispatch)
- Pytest progress output for partial regression runs (the key signal — 52/52 targeted tests passing — is preserved)

## 11. Model selection history

| Task | Model Used | Tokens (est.) | Outcome |
|---|---|---|---|
| Node 30 Phase 0 + Phase 1 (scope + compact protocol) | Opus | ~12k | Success; three scope questions yielded clean decisions |
| Node 30 Phase 2 (implementation) | Opus | ~18k | 1 fix needed (L1 posterior dispatch when kernels=None) — fixed on 2nd iteration |
| Node 30 Phase 3 (audit via regression tests) | Opus | ~5k | APPROVED after 2 test rewrites |
| Node 31 /scientific-advisor brief | Opus | ~10k | Key finding: scope correction to M2+L3 |
| Node 31 /architect protocol | Opus | ~8k | G1 passed |

**Total session tokens (est.):** ~55k
**Would-have-been all-Opus:** same (all phases here are Opus-mandated
per matrix — architecture design, full audit, protocol generation all
non-negotiably Opus; implementation was Opus because Node 30 was a
novel refactor touching UQ semantics).

## 12. Roadmap position

- **Current milestone:** v7.1.0-dev Nodes 30–32 (dev-orchestrator campaign initiated this session)
- **Modules completed this milestone:** 1 of 3 (Node 30)
- **Estimated remaining effort:**
  - Node 31: ~25k tokens implementation + audit, 1 fresh session
  - Node 32: ~10k tokens roadmap doc, part of same or next session
- **v7.0 release gate:** still blocked on Study A (Node 21 / F1 wet-lab data). Unchanged.
- **Process observations:**
  - The dev-orchestrator framework worked as designed — context compression discipline prompted this handover rather than a context-exhaustion failure.
  - Scientific brief before architect protocol caught a scope error (L3 vs M2 for EDC/NHS) that a rushed implementation would have missed.
  - Auto-test fallback in Node 30 posterior dispatch (handling None kernels) was a small trap that the test harness caught cleanly — validated the tight test-first loop.

---

*This handover is self-contained. A new dialogue can resume Node 31
implementation using only this document + `docs/node31_protocol.md`,
`docs/node31_edc_nhs_scientific_brief.md`, and the current repo state.*
