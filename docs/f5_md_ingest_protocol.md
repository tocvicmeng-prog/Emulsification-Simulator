# F5 — MD Parameter Ingest (MARTINI): /architect Protocol

**Prepared by:** /architect (via dev-orchestrator)
**Date:** 2026-04-17
**Status:** Protocol + Phase 1 implementation shipped same turn.
**Scope basis:** `docs/01_scientific_advisor_report.md §A.3` §5.

---

## 1. Purpose and scope

Expose coarse-grained MD (MARTINI) parameter tables as first-class
inputs to the EmulSim property pipeline. Users running MARTINI 3
simulations of polymer/solvent/surfactant systems produce:

- Effective χ (Flory-Huggins) parameters between bead types
- Effective interfacial tensions
- Pair-potential U(r) tables for specific bead pairs
- Radius-of-gyration / persistence-length estimates

F5 Phase 1 ingests a **JSON-serialised record** of these
MD-derived quantities and maps them onto existing EmulSim
MaterialProperties fields (primarily the χ parameters used by
cellulose NIPS). This is the **ingest-only default** agreed at the
Node 32 roadmap — we do NOT run MARTINI ourselves; we accept
user-supplied output.

Out of scope for F5 Phase 1:
- Running MARTINI or any MD package from within EmulSim.
- Pair-potential tabulated-force ingestion (beyond recording the
  values; there is no downstream consumer yet).
- Automatic MARTINI ↔ EmulSim bead-type mapping (user provides the
  mapping explicitly in the record).

## 2. Scientific basis

MARTINI is a coarse-grained force field (Marrink et al. 2007 *J.
Phys. Chem. B* 111:7812; Souza et al. 2021 *Nat. Methods* 18:382).
4-to-1 atom-to-bead mapping with a discrete bead-type library
(C1-C5 apolar, P1-P5 polar, etc.). Effective χ between two bead
types can be extracted from the Lennard-Jones ε matrix via the
standard Flory-Huggins mapping:

```
χ_ij ≈ z · (ε_ii + ε_jj − 2·ε_ij) / (2·kT)
```

where z is the bead-coordination number (~12 in MARTINI's
face-centred-cubic packing). For a polymer-solvent system with
many beads per chain, this gives an effective chain-averaged χ
suitable for dropping into a continuum Flory-Huggins model.

References:
- Marrink et al. (2007) *J. Phys. Chem. B* 111:7812 — canonical
  MARTINI 2 paper.
- Souza et al. (2021) *Nat. Methods* 18:382 — MARTINI 3.
- Lee & Larson (2009) *Macromolecules* 42:7528 — χ extraction
  for polymer blends.

## 3. Deliverables (Phase 1)

| File | Purpose | LOC |
|---|---|---|
| `src/emulsim/md_ingest.py` | `MartiniRecord` dataclass + `load_martini_record(path)` + `apply_chi_to_props(props, record)` | ~220 |
| `tests/test_md_ingest.py` | 10 tests | ~200 |
| `data/validation/md/example_martini_cellulose.json` | Reference record for test fixtures | ~30 |

Total: ~450 LOC.

## 4. Schema (JSON)

```json
{
  "source": "martini_3.0",
  "paper_doi": "10.1038/s41592-021-01098-3",
  "system_description": "cellulose / NaOH-urea / water (NaOH-urea lumped as P4)",
  "beads": {
    "polymer": "P3",
    "solvent": "P4",
    "nonsolvent": "P4"
  },
  "chi": {
    "polymer_solvent":    0.41,
    "polymer_nonsolvent": 0.88,
    "solvent_nonsolvent": 0.28
  },
  "diagnostics": {
    "n_chains":        64,
    "chain_length":    50,
    "simulation_ns":   500,
    "temperature_K":   298.15,
    "n_water_beads":   4096
  },
  "notes": "Run on MARTINI 3.0 force field; χ extracted via ε matrix → Flory-Huggins."
}
```

All top-level keys are required except `paper_doi` and `notes`.
Missing χ sub-keys are interpreted as "MD didn't measure this" and
left untouched on props; present χ values overwrite.

## 5. API

```python
from emulsim.md_ingest import (
    MartiniRecord,
    load_martini_record,
    apply_chi_to_props,
)

record = load_martini_record("path/to/file.json")
# record.chi is a dict {"polymer_solvent": 0.41, ...}

props = MaterialProperties()
apply_chi_to_props(props, record)
# props.chi_PS_cellulose, chi_PN_cellulose, chi_SN_cellulose are now
# overwritten by the MD-derived values.
```

Phase 2 (deferred): MaterialProperties adapter for non-cellulose
families (e.g., PLGA/solvent χ, alginate/water χ) — mostly a matter
of adding mapping keys to `apply_chi_to_props`.

## 6. Tests (Phase 1 scope)

1. Load a valid JSON file → populated `MartiniRecord`.
2. Missing required key → `ValueError`.
3. Extra unknown keys → preserved in `record.extra` (forward-compat).
4. χ sub-keys are optional (partial records accepted).
5. `apply_chi_to_props` patches cellulose χ fields correctly.
6. `apply_chi_to_props` leaves non-cellulose fields untouched.
7. χ with non-finite value (NaN, inf) → raises `ValueError` on apply.
8. Negative χ allowed (χ < 0 means "better than athermal", physically
   valid for attractive systems).
9. Round-trip: save → load → identical `MartiniRecord`.
10. Reference fixture file parses without error and contains three χ
    values.

## 7. Gate G1 (12-point completeness)

All 12 met for Phase 1 scope. Algorithm is a JSON parse + dict
merge; no new numerics; mapping is explicit and user-controlled.

## 8. Phase 2 (deferred)

- Tabulated U(r) pair-potential ingestion (needs a downstream
  consumer that's not yet designed).
- Automatic bead-type mapping from MARTINI → EmulSim polymer family.
- CalibrationStore integration so MD-derived χ propagates through
  the existing posterior-sampling infrastructure.
- Reverse direction: emit a MARTINI-shaped record FROM an EmulSim
  CalibrationStore fit (for cross-comparison with MD groups).
