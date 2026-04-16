# EmulSim Validation Datasets

**Purpose:** anchor the kernel constants, empirical correlations, and isotherm
parameters in the simulator against real wet-lab measurements so that the
evidence tier of model outputs can graduate from `SEMI_QUANTITATIVE` to
`CALIBRATED_LOCAL` (and ultimately `VALIDATED_QUANTITATIVE`).

This directory was scaffolded by Node 20 (v7.0, P1a). The actual data ingest
(Node 21, F1 closure) is gated on Study A wet-lab delivery.

## Directory map

| Subdirectory | Calibrates | Source assays (per docs/35 §9, docs/10 §5) |
|---|---|---|
| `l1_dsd/` | L1 PBE kernel constants (C1, C2, C3, C4, C5) | DROPLET_SIZE_DISTRIBUTION + INTERFACIAL_TENSION + DISPERSED_VISCOSITY |
| `l2_pore/` | L2 empirical pore-size coefficients | PORE_SIZE + POROSITY + GELATION_ONSET |
| `l3_kinetics/` | L3 rate constants (k_xlink_0, E_a, f_bridge) | CROSSLINK_CONVERSION |
| `l4_mechanics/` | L4 modulus prefactors + IPN coupling | BULK_MODULUS |
| `m2_capacity/` | M2 ligand density and activity retention | LIGAND_DENSITY + ACTIVITY_RETENTION + STATIC_BINDING_CAPACITY |

## File conventions

Each directory contains:

1. **`assays/`** — raw `AssayRecord` JSON files (Node 17 schema).
   Filename pattern: `<study_id>_<kind>.json`, e.g. `2026Q2_RPM_sweep.json`.

2. **`fits/`** — output of the calibration fitter (Node 21).
   Each fit produces a `CalibrationStore`-compatible JSON file:
   `fit_<dataset_id>_<timestamp>.json`. These can be loaded directly via
   `CalibrationStore.load_json()` and applied at run time via
   `RunContext(calibration_store=...)`.

3. **`schema.json`** — JSON-Schema for the AssayRecord variant accepted by
   that subdirectory's fitter (Node 21+, Sprint 1 of v7.0).

## Workflow

```
[Wet lab]                                      [Simulator]

  Run assay -- Sprint 1                          Use calibration -- v6.1+
       │                                              ▲
       ▼                                              │
  AssayRecord JSON (Node 17 schema) ─┐                │
                                     │                │
       ┌─────────────────────────────┘                │
       ▼                                              │
  data/validation/<level>/assays/                     │
       │                                              │
       ▼                                              │
  Fitter (Node 21) -- per-level kernel/coefficient    │
                       fitting + posterior std        │
       │                                              │
       ▼                                              │
  data/validation/<level>/fits/  -- CalibrationStore  │
                                    JSON files ──────┘
```

## Status (2026-04-17)

| Subdirectory | Assays present | Fits present | Tier achievable |
|---|---|---|---|
| `l1_dsd/` | 0 (gated on Study A) | — | SEMI_QUANTITATIVE only |
| `l2_pore/` | 0 | — | SEMI_QUANTITATIVE only |
| `l3_kinetics/` | 0 | — | SEMI_QUANTITATIVE only |
| `l4_mechanics/` | 0 | — | SEMI_QUANTITATIVE only |
| `m2_capacity/` | 0 | — | SEMI_QUANTITATIVE only |

The platform produces SEMI_QUANTITATIVE outputs on the default config. To
graduate to CALIBRATED_LOCAL, drop the Sprint-1 wet-lab outputs into the
appropriate `assays/` directory, run the fitter to produce a `fits/` JSON,
and pass `RunContext(calibration_store=CalibrationStore.load_json(...))`
to the orchestrator.
