# EmulSim Module 2 v5.9 — Final Plan (Post-Audit)

## Contract-First, Then Chemistry

**Date:** 2026-04-12
**Status:** Final — incorporates all 14 audit findings from doc 29
**Baseline:** v5.8 (42 profiles, 6 step types, 156 tests passing)
**Principle:** Fix the M2→M3 bridge FIRST, then expand chemistry.

---

## 1. Audit-Driven Restructuring

The auditor's central criticism: **"If M2 produces chemically detailed media but M3 still uses generic Langmuir defaults, the integrated simulator will appear more advanced than it actually is. Fix the contract and routing first, then expand chemistry."**

**Original plan (doc 28):** 6 must-haves in parallel, ~58h.
**Revised plan:** 5 sub-releases in dependency order, ~65h total.

### Disposition of All 14 Findings

| ID | Severity | Finding | Action |
|---|---|---|---|
| **F1** | Critical | M3 routing needs isotherm adapters, not just dispatch | **ACCEPT** — v5.9.0 core: adapters + transport interface |
| **F2** | High | accessible_area needs separate reagent/ligand bases | **ACCEPT** — Add both bases to FMC |
| **F3** | High | Metal state must live on material, not reagent profile | **ACCEPT** — Store on FunctionalMicrosphere after charging |
| **F4** | High | EDC/NHS needs executable AHA carboxyl-distal first | **ACCEPT** — Move EDC/NHS to v5.9.3 (after carboxyl-distal validated) |
| **F5** | Med-High | EDC/NHS pseudo-single-step is ranking_only | **ACCEPT** — Label accordingly |
| **F6** | High | Protein pretreatment needs typed dataclass | **ACCEPT** — ProteinPretreatmentState, not loose dict |
| **F7** | Med-High | Washing should be advisory screening only | **ACCEPT** — Scope to advisory model |
| **F8** | Medium | pH scaling pKa-only is incomplete | **ACCEPT** — Implement disabled by default |
| **F9** | Med-High | Irreversible affinity needs kinetic model | **ACCEPT** — dq/dt = k_ads * C * (q_max - q), not Langmuir flag |
| **F10** | Medium | Target defaults need confidence tiers | **ACCEPT** |
| **F11** | Medium | PSD/uncertainty limitations acknowledged | **ACCEPT** — Add capacity warnings |
| **F12** | High | UI plan incomplete | **ACCEPT** — Dedicated UI node per sub-release |
| **F13** | High | FMC needs coherent extension | **ACCEPT** — Design FMC v2 first in v5.9.0 |
| **F14** | Med-High | Scope too broad for single release | **ACCEPT** — Split into 5 sub-releases |

---

## 2. Revised Release Structure

| Sub-Release | Focus | Items | LOC | Depends On |
|---|---|---|---|---|
| **v5.9.0** | FMC v2 + M3 routing + accessible area | D2, D3, D4, D8, D13 | ~600 | Nothing (foundation) |
| **v5.9.1** | IMAC metal charging + state | D1 | ~300 | v5.9.0 (FMC v2) |
| **v5.9.2** | Protein pretreatment + maleimide completion | D11 | ~250 | v5.9.0 |
| **v5.9.3** | EDC/NHS + AHA carboxyl-distal path | D5 | ~300 | v5.9.0 + carboxyl-distal |
| **v5.9.4** | Washing advisory + pH scaling | D6, D17 | ~350 | v5.9.0 |
| **Total** | | 10 items | **~1,800** | |

### Why This Order

1. **v5.9.0 first** — without FMC v2 and M3 routing, new chemistry (metal charging, EDC/NHS) produces outputs M3 can't consume correctly.
2. **v5.9.1 after v5.9.0** — IMAC metal state needs FMC v2 to propagate `metal_loaded_fraction` to M3.
3. **v5.9.2 parallel with v5.9.1** — protein pretreatment is independent of IMAC but needs FMC v2 for uncertainty propagation.
4. **v5.9.3 after v5.9.0** — EDC/NHS needs the carboxyl-distal ACS path validated, which requires SPACER_ARM creating CARBOXYL_DISTAL.
5. **v5.9.4 last** — washing and pH scaling are enhancements that don't affect the core pipeline.

---

## 3. v5.9.0 — FMC v2 + M3 Routing (Foundation)

### 3.1 FunctionalMediaContract v2

New fields (consolidated per audit F13):

```python
# ── Accessible area (audit F2) ──
reagent_accessible_area_per_bed_volume: float = 0.0  # [m2/m3 bed]
ligand_accessible_area_per_bed_volume: float = 0.0   # [m2/m3 bed]
capacity_area_basis: str = ""  # "reagent_accessible" or "ligand_accessible"

# ── Uncertainty (audit D4) ──
activity_retention_uncertainty: float = 0.0
q_max_lower: float = 0.0
q_max_upper: float = 0.0

# ── Metal state (audit F3 prep) ──
metal_ion: str = ""
metal_loaded_fraction: float = 0.0

# ── M3 integration (audit F1) ──
m3_support_level: str = "not_mapped"
final_ligand_profile_key: str = ""
process_state_requirements: str = ""  # "salt_concentration", "imidazole", ""
residual_reagent_warnings: list[str] = field(default_factory=list)
```

### 3.2 M3 Isotherm Adapters (audit F1)

The critical insight: M3's LRM solver calls `isotherm.equilibrium_loading(C)` but SMA needs salt and IMAC needs imidazole. Solution: **EquilibriumAdapter pattern**.

```python
class EquilibriumAdapter:
    """Adapts multi-parameter isotherms to single-C interface for LRM solver."""
    def __init__(self, isotherm, process_state: dict):
        self._isotherm = isotherm
        self._state = process_state

    def equilibrium_loading(self, C):
        if isinstance(self._isotherm, SMAIsotherm):
            return self._isotherm.equilibrium_loading(C, self._state["salt_concentration"])
        elif isinstance(self._isotherm, IMACCompetitionIsotherm):
            return self._isotherm.equilibrium_loading(C, self._state.get("imidazole", 0.0))
        else:
            return self._isotherm.equilibrium_loading(C)
```

### 3.3 Irreversible Adsorption (audit F9)

NOT a Langmuir flag. Separate kinetic model:

```python
class IrreversibleAdsorptionIsotherm:
    """One-way binding: dq/dt = k_ads * C * (q_max - q), no desorption."""
    def __init__(self, q_max, k_ads):
        self.q_max = q_max
        self.k_ads = k_ads

    def equilibrium_loading(self, C):
        return self.q_max  # at equilibrium, all sites filled (irreversible)

    def kinetic_rate(self, C, q):
        return self.k_ads * C * max(self.q_max - q, 0.0)
```

### 3.4 Target Protein Defaults (D13)

```python
TARGET_PROTEIN_LIBRARY = {
    "IgG": TargetProtein(mw=150000, r_h=5.3e-9, pI=7.5, z_eff=5, K_affinity=1e8,
                          confidence="estimated", notes="Generic human IgG1"),
    "His6_50kDa": TargetProtein(mw=50000, r_h=2.8e-9, K_imac=1e4,
                                 confidence="estimated"),
    "GST_fusion": TargetProtein(mw=75000, r_h=3.5e-9, K_gst=1e3,
                                 confidence="estimated"),
    ...
}
```

### 3.5 Work Nodes for v5.9.0

| WN | Description | Tier | LOC |
|---|---|---|---|
| WN-0a | FMC v2 dataclass extension (all new fields) | Opus | 80 |
| WN-0b | FMC builder: accessible area per bed volume + uncertainty | Sonnet | 60 |
| WN-0c | EquilibriumAdapter + IrreversibleAdsorption in M3 | Opus | 150 |
| WN-0d | M3 routing dispatch: binding_model_hint → adapter | Opus | 120 |
| WN-0e | Target protein defaults library | Haiku | 80 |
| WN-0f | UI: expose SPACER_ARM + all v5.8 profiles in dropdowns | Sonnet | 60 |
| WN-0g | Tests: routing, adapter, area basis, uncertainty | Opus | 120 |

---

## 4. v5.9.1 — IMAC Metal Charging

### Key Change from Original Plan (audit F3)

Metal state stored on **FunctionalMicrosphere material state**, not on the reagent profile.

```python
@dataclass
class MetalChargingState:
    chelator_site_type: ACSSiteType  # EPOXIDE profile with IDA/NTA coupled
    metal_ion: str                   # "Ni2+", "Co2+", etc.
    metal_loaded_fraction: float     # [0,1] from equilibrium calc
    metal_capacity_density: float    # [mol/m2] metal per area
    charging_pH: float
    stripping_agent: str = ""        # "EDTA", "imidazole"
    stripped_fraction: float = 0.0
```

### Work Nodes

| WN | Description | Tier | LOC |
|---|---|---|---|
| WN-1a | METAL_CHARGING step type + equilibrium solver | Sonnet | 120 |
| WN-1b | MetalChargingState dataclass + orchestrator integration | Sonnet | 80 |
| WN-1c | 5 metal/stripping profiles (Ni, Co, Cu, Zn, EDTA) | Haiku | 100 |
| WN-1d | FMC builder reads metal state from material | Sonnet | 40 |
| WN-1e | Tests: charging, stripping, FMC propagation | Sonnet | 80 |

---

## 5. v5.9.2 — Protein Pretreatment

### Typed State (audit F6)

```python
@dataclass
class ProteinPretreatmentState:
    protein_key: str = ""
    free_thiol_fraction: float = 0.0
    activity_after_reduction: float = 1.0
    reductant_used: str = ""           # "TCEP", "DTT"
    excess_reductant_removed: bool = False
    time_since_reduction_s: float = 0.0
    reoxidation_rate: float = 1e-5     # [1/s] default slow reoxidation
    warnings: list[str] = field(default_factory=list)
```

### Work Nodes

| WN | Description | Tier | LOC |
|---|---|---|---|
| WN-2a | PROTEIN_PRETREATMENT step type + solver | Sonnet | 100 |
| WN-2b | ProteinPretreatmentState dataclass | Sonnet | 40 |
| WN-2c | Protein coupling reads thiol state from pretreatment | Sonnet | 30 |
| WN-2d | 2 profiles (TCEP, DTT) + validation rules | Sonnet | 60 |
| WN-2e | Tests: reduction, reductant removal check, coupling integration | Sonnet | 80 |

---

## 6. v5.9.3 — EDC/NHS + Carboxyl-Distal Path

### Prerequisite (audit F4)

Must first implement AHA spacer-arm creating CARBOXYL_DISTAL:

```
EPOXIDE → AHA(SPACER_ARM) → CARBOXYL_DISTAL → EDC/NHS(ACTIVATION) → NHS_ESTER → amine coupling
```

New ACSSiteType: `CARBOXYL_DISTAL` (distinct from native CARBOXYL to track provenance).

### Work Nodes

| WN | Description | Tier | LOC |
|---|---|---|---|
| WN-3a | Add CARBOXYL_DISTAL + NHS_ESTER to ACSSiteType | Sonnet | 10 |
| WN-3b | AHA spacer-arm profile (EPOXIDE → CARBOXYL_DISTAL) | Sonnet | 40 |
| WN-3c | EDC/NHS activation profile (CARBOXYL_DISTAL → NHS_ESTER) | Opus | 80 |
| WN-3d | NHS-ester hydrolysis modeling (Template 2) | Opus | 60 |
| WN-3e | Tests: full AHA → EDC/NHS → amine coupling path | Opus | 100 |

---

## 7. v5.9.4 — Washing Advisory + pH Scaling

### Washing (audit F7: advisory only)

Scoped to simple diffusion-out screening model. NOT GMP compliance.

```python
class WashingResult:
    residual_concentration: dict[str, float]  # reagent → [ppm]
    passes_screening: bool
    screening_basis: str = "diffusion_out_advisory"
    warnings: list[str]
```

### pH Scaling (audit F8: disabled by default)

```python
def _ph_rate_scaling(ph: float, pKa: float, n_hill: float = 1.0) -> float:
    """Sigmoid scaling for nucleophile ionization fraction."""
    if pKa <= 0:
        return 1.0  # disabled
    fraction_active = 1.0 / (1.0 + 10 ** (n_hill * (pKa - ph)))
    return fraction_active
```

### Work Nodes

| WN | Description | Tier | LOC |
|---|---|---|---|
| WN-4a | WASHING step type + advisory solver | Sonnet | 120 |
| WN-4b | wash_buffer profile + regulatory_limit_ppm field | Haiku | 40 |
| WN-4c | pH scaling helper + per-chemistry pKa fields | Sonnet | 120 |
| WN-4d | Tests: washing, pH scaling, residual tracking | Sonnet | 100 |

---

## 8. Complete Dependency Graph

```
v5.9.0 (FMC v2 + M3 routing)
  │
  ├──→ v5.9.1 (IMAC metal charging)
  │
  ├──→ v5.9.2 (Protein pretreatment)
  │
  ├──→ v5.9.3 (EDC/NHS + carboxyl-distal)
  │
  └──→ v5.9.4 (Washing + pH scaling)
```

All sub-releases depend on v5.9.0. Sub-releases 5.9.1-5.9.4 are independent of each other and can be implemented in any order after v5.9.0.

---

## 9. Post-v5.9 Metrics

| Metric | v5.8 | v5.9 |
|---|---|---|
| Profiles | 42 | **~56** (+5 metal, +2 reduction, +2 EDC/NHS, +1 wash, +4 target defaults) |
| Step types | 6 | **9** (+METAL_CHARGING, PROTEIN_PRETREATMENT, WASHING) |
| ACSSiteType | 8 | **10** (+CARBOXYL_DISTAL, NHS_ESTER) |
| FMC fields | ~20 | **~30** (v2 extension) |
| M3 isotherms | 5 | **6** (+IrreversibleAdsorption) |
| Tests | 156 | **~210** |

---

## 10. What Remains Deferred to v6.0+

| Item | Reason |
|---|---|
| D7: Salt-dependent HIC isotherms | Requires user calibration framework |
| D9: TMAE via DVS | Kinetic data unavailable |
| D10: Lectin-specific elution | Sugar competition constants needed |
| D14: Batch-to-batch variability | M1 interface changes |
| D15: Pore-size distribution | M1 PSD output |
| D16: Protein leaching/deactivation | Long-term stability modeling |
| D12: Metal leaching detail | Cycle-dependent kinetics |

---

## 11. Effort Estimate

| Sub-Release | Implementation | Testing | Total |
|---|---|---|---|
| v5.9.0 | 18h | 8h | 26h |
| v5.9.1 | 8h | 4h | 12h |
| v5.9.2 | 6h | 4h | 10h |
| v5.9.3 | 8h | 5h | 13h |
| v5.9.4 | 8h | 4h | 12h |
| **Total** | **48h** | **25h** | **~73h** |

---

> **Disclaimer**: This is an expanded semi-quantitative simulation library. All isotherms, kinetic constants, and process parameters are estimates requiring calibration. v5.9 improves the M2→M3 bridge fidelity but does not replace measured resin characterization or validated process models.
