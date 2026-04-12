# EmulSim v6.0 — Final Plan (Post-Audit)

## Calibration-Enabled M2/M3 Process Simulation

**Date:** 2026-04-12
**Status:** Final — incorporates all 15 audit findings from doc 32
**Baseline:** v5.9 (52 profiles, 9 step types, 10 ACSSiteType, 156 tests)
**Principle:** Harden v5.9 routing first, then build calibration/uncertainty/lifetime frameworks.

---

## 1. Audit Finding Disposition

| ID | Severity | Finding | Action |
|---|---|---|---|
| **F1** | Positive | Framework-first strategy correct | Acknowledged — proceed |
| **F2** | High | Calibration schema under-specified | **ACCEPT** — Typed schema with units, target, validity domain |
| **F3** | Critical | `uncertainty.py` naming conflict | **ACCEPT** — Rename to `uncertainty_propagation/` |
| **F4** | High | M1 default CVs need correlation handling | **ACCEPT** — Two-tier: measured (Tier 1) + assumed screening (Tier 2) |
| **F5** | Medium | Existing `pore_size_std` partially usable | **ACCEPT** — Use for log-normal PSD approximation |
| **F6** | Med-High | Lifetime must be calibrated empirical, not predictive | **ACCEPT** — Label as empirical projection |
| **F7** | High | HIC salt coefficient unit risk | **ACCEPT** — Calibration-first; require user K_0/m_salt |
| **F8** | Med-High | Lectin models need cofactor/oligomer metadata | **ACCEPT** — `requires_user_calibration` mandatory |
| **F9** | High | TMAE chemistry still unresolved | **ACCEPT** — Keep deferred to v7.0 |
| **F10** | Critical | M3 SMA adapter route currently fails | **ACCEPT** — v5.9 hardening gate BEFORE v6.0 |
| **F11** | High | Multi-component gradients need ProcessState dataclass | **ACCEPT** — First-class `ProcessState`, not loose dict |
| **F12** | High | UI must expose all 9 step types | **ACCEPT** — UI gate before public release |
| **F13** | Med-High | Calibration UI must be must-have | **ACCEPT** — Minimal import/inspect UI is must-have |
| **F14** | Medium | PSD only if M1 provides it | **ACCEPT** — Log-normal from `pore_size_std` first |
| **F15** | Medium | pH scaling can double-count | **ACCEPT** — Normalize: `k_eff(ph_optimum) = k_forward` |

---

## 2. Revised Release Structure

The audit's central requirement: **v5.9 hardening gate before any v6.0 framework work**.

| Sub-Release | Focus | Depends On |
|---|---|---|
| **v5.9.5** | Hardening: fix M3 adapter routing, array safety, UI gaps, uncertainty.py naming | Nothing |
| **v6.0-alpha** | Calibration framework + M1 uncertainty contract | v5.9.5 |
| **v6.0-beta** | ProcessState + multi-component gradients + HIC/IMAC validation | v6.0-alpha |
| **v6.0-rc** | Lifetime model + UI overhaul + lectin/PSD | v6.0-beta |
| **v6.0** | Integration testing + documentation | v6.0-rc |

---

## 3. v5.9.5 Hardening Gate (MANDATORY FIRST)

### Must Fix Before v6.0

| # | Issue | Fix |
|---|---|---|
| H1 | `SMAIsotherm` constructor incompatible with adapter | Fix `select_isotherm_from_fmc()` to use correct SMA constructor args |
| H2 | `IrreversibleAdsorptionIsotherm.equilibrium_loading()` not array-safe | Add `np.where(C > 0, q_max, 0.0)` for NumPy input |
| H3 | `run_breakthrough()` doesn't auto-select isotherm from FMC | Add optional `fmc` parameter to `run_breakthrough()` |
| H4 | UI missing METAL_CHARGING, PROTEIN_PRETREATMENT, WASHING | Add to Chemistry dropdown |
| H5 | `uncertainty.py` blocks future `uncertainty_propagation/` package | Rename existing module |
| H6 | Gradient elution doesn't pass process state to equilibrium | Wire gradient values into adapter's process_state |

**Estimated effort:** ~15h

---

## 4. v6.0-alpha: Calibration + Uncertainty Frameworks

### 4.1 Calibration Framework (F1)

**Typed schema (audit F2):**

```python
@dataclass
class CalibrationEntry:
    profile_key: str              # "protein_a_coupling"
    parameter_name: str           # "q_max", "K_L", "activity_retention"
    measured_value: float
    units: str                    # "mg/mL", "mol/m3", "fraction"
    target_molecule: str = ""     # "IgG1", "BSA", "His6-GFP"
    temperature_C: float = 25.0
    ph: float = 7.0
    salt_concentration_M: float = 0.0
    salt_type: str = ""
    measurement_type: str = ""    # "static_binding", "DBC10", "batch_uptake"
    confidence: str = "measured"  # "measured", "literature", "estimated"
    source_reference: str = ""
    replicates: int = 1
```

**CalibrationStore:** Load/save JSON, query by profile_key + conditions.
**apply_calibration():** Overrides FMC defaults; logs every override.
**UI:** Minimal import (JSON upload) + inspection panel showing active overrides.

### 4.2 M1 Uncertainty Interface (F2)

**Package:** `src/emulsim/uncertainty_propagation/` (renamed per audit F3)

**M1UncertaintyContract:**

```python
@dataclass
class M1UncertaintyContract:
    cv_bead_d50: float = 0.0     # Coefficient of variation [0,1]
    cv_porosity: float = 0.0
    cv_pore_size: float = 0.0
    cv_nh2_bulk: float = 0.0
    cv_oh_bulk: float = 0.0
    tier: str = "assumed"         # "measured" or "assumed"
    correlation_matrix: Optional[np.ndarray] = None  # 5x5 if measured
```

**Monte Carlo wrapper:**

```python
def run_with_uncertainty(contract, steps, n_samples=100, uncertainty=None):
    results = []
    for _ in range(n_samples):
        perturbed = _perturb_contract(contract, uncertainty)
        result = orchestrator.run(perturbed, steps)
        results.append(result)
    return UncertaintyResult(
        median=..., p5=..., p95=..., cv_q_max=...
    )
```

**Sanity constraints:** Reject samples with negative concentrations, porosity > 1, d50 < 0.

---

## 5. v6.0-beta: ProcessState + Gradients + HIC

### 5.1 ProcessState (audit F11)

```python
@dataclass
class ProcessState:
    """Time-varying process conditions for M3 simulations."""
    salt_concentration: float = 0.0     # [mol/m3]
    salt_type: str = ""
    ph: float = 7.0
    imidazole: float = 0.0              # [mol/m3]
    sugar_competitor: float = 0.0       # [mol/m3] (lectin elution)
    sugar_type: str = ""
    temperature: float = 298.15         # [K]
    conductivity: float = 0.0           # [mS/cm]
```

Replace all `process_state: dict` with `ProcessState` in adapters and routing.

### 5.2 HIC Isotherm (audit F7: calibration-first)

```python
class HICIsotherm:
    """Salt-modulated HIC binding: K_eff = K_0 * exp(m * C_salt)"""
    def __init__(self, q_max, K_0, m_salt, salt_type="ammonium_sulfate"):
        ...
    def equilibrium_loading(self, C, salt_concentration):
        K_eff = self.K_0 * math.exp(self.m_salt * salt_concentration)
        return self.q_max * K_eff * C / (1 + K_eff * C)
```

**Default:** `m3_support_level="requires_user_calibration"` unless user supplies K_0 and m_salt via calibration framework.

### 5.3 IMAC Imidazole Gradient Validation

End-to-end test: load Ni-NTA column → bind His-tag protein at 0 mM imidazole → gradient elute 0→500 mM → verify protein elutes between 100-300 mM.

---

## 6. v6.0-rc: Lifetime + UI + Lectin

### 6.1 Lifetime Model (audit F6: empirical)

```python
@dataclass
class LifetimeProjection:
    initial_capacity: float         # [mg/mL] or [mol/m3]
    k_deactivation: float           # [1/cycle]
    projected_cycles_to_80pct: int  # cycles until 80% of initial capacity
    assumption_notes: str           # CIP conditions, feed, storage
    confidence: str = "empirical"   # always empirical unless calibrated
```

### 6.2 UI Restructure (audit F12)

Split `app.py` (~1500 lines) into:

```
visualization/
    app.py              # Main entry + page config + sidebar
    tab_m1.py           # M1 Fabrication tab
    tab_m2.py           # M2 Functionalization tab
    tab_m3.py           # M3 Performance tab
    panels/
        calibration.py  # Calibration import/inspect
        uncertainty.py  # Uncertainty configuration
        lifetime.py     # Lifetime projection display
```

### 6.3 Lectin Competition (audit F8: calibration-required)

Reuse `CompetitiveAffinityIsotherm` (generalized from IMAC):

```python
class CompetitiveAffinityIsotherm:
    """Generic competitive binding: protein vs competitor for ligand sites."""
    def __init__(self, q_max, K_protein, K_competitor):
        ...
    def equilibrium_loading(self, C_protein, C_competitor):
        return q_max * K_protein * C_protein / (1 + K_protein * C_protein + K_competitor * C_competitor)
```

For Con A: `competitor = mannose`, `K_competitor ~ 1e3`.
For WGA: `competitor = GlcNAc`, `K_competitor ~ 1e2`.

---

## 7. What Defers to v7.0

| Item | Reason |
|---|---|
| D9: TMAE via DVS | Kinetic data unavailable; Q covers strong anion exchange |
| Mechanistic pH-temperature coupling | Research-grade; sigmoid scaling sufficient for v6.0 |
| Gradient optimization (automated method dev) | Requires objective function framework |
| GMP-validated washing compliance | Regulatory scope, not simulation scope |
| PDB-based protein orientation prediction | Research collaboration needed |
| Full correlated M1 uncertainty (>5 variables) | M1 interface expansion |

---

## 8. Effort Estimates

| Sub-Release | LOC | Hours |
|---|---|---|
| v5.9.5 Hardening | ~200 | 15h |
| v6.0-alpha (calibration + uncertainty) | ~800 | 35h |
| v6.0-beta (ProcessState + HIC + IMAC validation) | ~600 | 25h |
| v6.0-rc (lifetime + UI split + lectin) | ~900 | 35h |
| v6.0 (integration testing + docs) | ~200 | 10h |
| **Total** | **~2,700** | **~120h** |

---

## 9. Post-v6.0 Metrics (Projected)

| Metric | v5.9 | v6.0 |
|---|---|---|
| Profiles | 52 | ~55 |
| Step types | 9 | 9 (no new) |
| ACSSiteType | 10 | 10 (no new) |
| New packages | 0 | 3 (calibration/, uncertainty_propagation/, lifetime/) |
| New isotherm classes | 1 (Irreversible) | +2 (HIC, CompetitiveAffinity) |
| UI files | 1 (app.py) | 7 (modular) |
| ProcessState | dict | **Typed dataclass** |
| Tests | 156 | ~250 |

---

## 10. Acceptance Criteria (from audit Section 7)

- [ ] v5.9.5: `select_isotherm_from_fmc()` works for ALL binding_model_hint values
- [ ] v5.9.5: UI exposes all 9 step types
- [ ] v5.9.5: `uncertainty.py` renamed, no import conflicts
- [ ] v6.0-alpha: CalibrationStore loads JSON with typed schema validation
- [ ] v6.0-alpha: Calibration overrides visible in FMC and UI
- [ ] v6.0-alpha: Monte Carlo uncertainty produces p5/p95 bounds on q_max
- [ ] v6.0-beta: HIC isotherm requires user calibration (K_0, m_salt)
- [ ] v6.0-beta: IMAC imidazole gradient test: protein elutes at correct imidazole
- [ ] v6.0-beta: ProcessState is typed dataclass, not dict
- [ ] v6.0-rc: Lifetime projection labeled "empirical" with assumption notes
- [ ] v6.0-rc: UI split into >=5 modules, app.py < 300 lines
- [ ] v6.0-rc: Lectin models require_user_calibration
- [ ] v6.0: All existing v5.9 workflows pass regression
- [ ] v6.0: No profile uses `m3_support_level="mapped_quantitative"` without calibration

---

> **Disclaimer**: v6.0 transitions EmulSim from a semi-quantitative chemistry simulator to a calibration-enabled process simulation platform. All uncalibrated outputs remain semi-quantitative. Calibrated outputs reflect the quality of user-supplied measurements, not theoretical predictions. Not a substitute for measured resin characterization or regulatory validation.
