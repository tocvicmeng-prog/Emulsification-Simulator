# EmulSim Module 2 Expansion — Best & Final Implementation Plan

## Incorporating Third-Party Audit Findings (doc 17)

**Version:** 2.0 (Revised)
**Date:** 2026-04-12
**Status:** Final — addresses all 9 audit findings (F1-F9)
**Prepared by:** Scientific Advisor + Architect + Dev-Orchestrator

---

## 1. Audit Findings Disposition

| ID | Severity | Finding | Disposition | Resolution Phase |
|---|---|---|---|---|
| **F1** | Critical | Quenching double-counts consumed + blocked | **ACCEPT** — Quenching increments `blocked_sites` ONLY, not both | Phase 1 |
| **F2** | Critical | "remaining_sites" vs "remaining activated sites" | **ACCEPT** — Add `remaining_activated` property; coupling/quenching operate on activated profiles only | Phase 1 |
| **F3** | High | Template 2 needs pH correction + hydrolysis tracking | **ACCEPT** — Add `hydrolyzed_sites` field; add pH validity gates (not full mechanistic pH model) | Phase 1, 3 |
| **F4** | High | Ligand identities under-specified | **ACCEPT** — Split reagent identity from installed ligand; add chemistry_class field | Phase 2 |
| **F5** | High | Protein coupling too simplified for production | **ACCEPT** — Label as "ranking_only"; make activity_retention user-adjustable with uncertainty | Phase 5 |
| **F6** | High | Unit consistency for steric limits | **ACCEPT** — Canonical unit = mol/particle; convert in ODE wrapper | Phase 3 |
| **F7** | High | M2→M3 mapping insufficient | **ACCEPT** — Add FunctionalMediaContract bridging M2→M3 | Phase 7 |
| **F8** | High | Validators must be backend, not UI-only | **ACCEPT** — Add backend validation in orchestrator.run() | Phase 6 |
| **F9** | Med-High | ODE wrappers need structured diagnostics | **ACCEPT** — Return CouplingResult dataclass with solver metadata | Phase 3 |

---

## 2. Corrected ACS State Model (Phase 1 — Foundation)

### 2.1 Current Conservation Law

```
remaining_sites = accessible_sites - consumed_sites - blocked_sites
validate(): consumed_sites + blocked_sites <= accessible_sites
```

**Problem:** This conflates all consumption modes. Coupling, crosslinking, hydrolysis, and quenching all increment `consumed_sites`, but they have different physical meanings and downstream effects.

### 2.2 Revised Conservation Law

Every accessible site ends in exactly ONE terminal state:

```
accessible_sites = remaining_sites
                 + crosslinked_sites    (consumed by crosslinking — contributes to G_DN)
                 + hydrolyzed_sites     (lost to aqueous hydrolysis — waste)
                 + ligand_coupled_sites (coupled to functional ligand)
                 + blocked_sites        (capped by quenching reagent)
```

**New invariant:**
```python
terminal_sum = (self.crosslinked_sites + self.hydrolyzed_sites
                + self.ligand_coupled_sites + self.blocked_sites)
assert terminal_sum <= self.accessible_sites * _CONSERVATION_TOL
remaining_sites = accessible_sites - terminal_sum
```

### 2.3 New/Modified ACSProfile Fields

```python
@dataclass
class ACSProfile:
    site_type: ACSSiteType
    total_sites: float = 0.0           # Bulk including buried [mol/particle]
    accessible_sites: float = 0.0      # Surface + pore accessible [mol/particle]
    activated_sites: float = 0.0       # Activated by chemistry (product profiles)

    # ── Terminal states (mutually exclusive destinations) ──
    crosslinked_sites: float = 0.0     # Consumed by crosslinking (contributes to G_DN)
    hydrolyzed_sites: float = 0.0      # Lost to hydrolysis (waste) [NEW]
    ligand_coupled_sites: float = 0.0  # Coupled to ligand
    ligand_functional_sites: float = 0.0  # Subset of coupled retaining activity
    blocked_sites: float = 0.0         # Capped by quenching

    # ── Densities ──
    total_density: float = 0.0         # [mol/m^2]
    accessible_density: float = 0.0    # [mol/m^2]

    # ── Metadata ──
    accessibility_model: str = "empirical_pore"
    uncertainty_fraction: float = 0.1
```

**Backward compatibility:** The old `consumed_sites` field is REMOVED. Existing crosslinking and activation workflows are updated to use `crosslinked_sites` instead.

### 2.4 Updated remaining_sites Property

```python
@property
def remaining_sites(self) -> float:
    terminal = (self.crosslinked_sites + self.hydrolyzed_sites
                + self.ligand_coupled_sites + self.blocked_sites)
    return max(self.accessible_sites - terminal, 0.0)
```

### 2.5 Updated remaining_activated Property (NEW)

```python
@property
def remaining_activated(self) -> float:
    """Activated sites not yet consumed by coupling, hydrolysis, or quenching."""
    used = self.ligand_coupled_sites + self.hydrolyzed_sites + self.blocked_sites
    return max(self.activated_sites - used, 0.0)
```

**This resolves F1 and F2.** Coupling and quenching operate on `remaining_activated`, not `remaining_sites`.

### 2.6 Updated validate() Method

```python
def validate(self) -> list[str]:
    errors = []
    if self.accessible_sites > self.total_sites * _CONSERVATION_TOL:
        errors.append("accessible > total")
    if self.activated_sites > self.accessible_sites * _CONSERVATION_TOL:
        errors.append("activated > accessible")
    terminal = (self.crosslinked_sites + self.hydrolyzed_sites
                + self.ligand_coupled_sites + self.blocked_sites)
    if terminal > self.accessible_sites * _CONSERVATION_TOL:
        errors.append(f"terminal sum ({terminal:.2e}) > accessible ({self.accessible_sites:.2e})")
    if self.ligand_functional_sites > self.ligand_coupled_sites * _CONSERVATION_TOL:
        errors.append("ligand_functional > ligand_coupled")
    coupled_plus = self.ligand_coupled_sites + self.hydrolyzed_sites + self.blocked_sites
    if coupled_plus > self.activated_sites * _CONSERVATION_TOL:
        errors.append("coupled+hydrolyzed+blocked > activated (for product profiles)")
    # Non-negativity
    for field in [self.total_sites, self.accessible_sites, self.activated_sites,
                  self.crosslinked_sites, self.hydrolyzed_sites,
                  self.ligand_coupled_sites, self.ligand_functional_sites,
                  self.blocked_sites]:
        if field < 0:
            errors.append(f"negative field value: {field}")
    return errors
```

### 2.7 Migration of Existing Workflows

| Workflow | Old Update | New Update |
|---|---|---|
| SECONDARY_CROSSLINKING | `consumed_sites += sites` | `crosslinked_sites += sites` |
| ACTIVATION (target) | `consumed_sites += sites` | `crosslinked_sites += sites` (OH consumed by activation chemistry) |
| ACTIVATION (product) | `accessible/total/activated += sites` | No change |

---

## 3. Revised Reagent Registry (Phase 2)

### 3.1 Extended ReagentProfile

Per audit F4: separate reagent identity from installed functional group.

```python
@dataclass
class ReagentProfile:
    # ── Identity ──
    name: str                          # Display name
    cas: str                           # CAS number of ACTUAL reactive reagent
    reagent_identity: str              # Chemical name of reactive form
    installed_ligand: str              # Functional group after coupling
    functional_mode: str               # "crosslinker", "activator", "iex_ligand",
                                       # "affinity_ligand", "hic_ligand", "imac_chelator",
                                       # "quencher"

    # ── Chemistry ──
    reaction_type: str                 # "crosslinking", "activation", "coupling",
                                       # "protein_coupling", "blocking"
    chemistry_class: str               # "epoxide_amine", "epoxide_thiol",
                                       # "vs_amine", "vs_thiol", "aldehyde_amine",
                                       # "nhs_amine", "reduction", "acetylation"
    target_acs: ACSSiteType
    product_acs: Optional[ACSSiteType] = None

    # ── Kinetics ──
    k_forward: float = 0.0            # [m^3/(mol*s)] at reference T
    E_a: float = 0.0                  # [J/mol]
    stoichiometry: float = 1.0        # [-]
    hydrolysis_rate: float = 0.0      # [1/s]

    # ── Validity windows ──
    ph_optimum: float = 7.0
    ph_min: float = 4.0               # [NEW] blocker if pH outside range
    ph_max: float = 10.0              # [NEW]
    temperature_default: float = 298.15  # [K]
    temperature_min: float = 273.15   # [NEW]
    temperature_max: float = 353.15   # [NEW]
    time_default: float = 3600.0      # [s]

    # ── Macromolecule fields (protein coupling only) ──
    ligand_mw: float = 0.0            # [Da]
    ligand_r_h: float = 0.5e-9        # [m] hydrodynamic radius
    is_macromolecule: bool = False     # True -> use ligand_accessible_area
    activity_retention: float = 1.0   # [0,1] — user-adjustable
    activity_retention_uncertainty: float = 0.15  # [NEW] +/- range
    max_surface_density: float = 0.0  # [mol/m^2] steric jamming limit

    # ── Metadata ──
    confidence_tier: str = "semi_quantitative"  # [NEW] or "calibrated"
    calibration_source: str = ""      # [NEW] literature reference
    hazard_class: str = ""            # [NEW] GHS hazard category
    notes: str = ""
```

### 3.2 Ten New Profiles

**Ligand Coupling (4):**

| Key | Reagent Identity | Installed Ligand | Mode | Chemistry Class | Target | CAS |
|---|---|---|---|---|---|---|
| `deae_coupling` | 2-(Diethylamino)ethylamine | DEAE (weak anion) | iex_ligand | epoxide_amine | EPOXIDE | 100-36-7 |
| `ida_coupling` | Iminodiacetic acid | IDA (chelator) | imac_chelator | epoxide_amine | EPOXIDE | 142-73-4 |
| `phenyl_coupling` | 1,4-Butanediol diglycidyl ether + phenylamine | Phenyl (HIC) | hic_ligand | epoxide_amine | EPOXIDE | 62-53-3 |
| `sp_coupling` | 1,3-Propane sultone | Sulfopropyl (strong cation) | iex_ligand | epoxide_amine | EPOXIDE | 1120-71-4 |

**Protein Coupling (2):**

| Key | Protein | Mode | Chemistry Class | Target | MW (kDa) | r_h (nm) | Activity Ret. |
|---|---|---|---|---|---|---|---|
| `protein_a_coupling` | Protein A (rSPA) | affinity_ligand | epoxide_amine | EPOXIDE | 42 | 2.5 | 0.60 +/- 0.15 |
| `protein_g_coupling` | Protein G (rSPG) | affinity_ligand | epoxide_amine | EPOXIDE | 22 | 2.0 | 0.65 +/- 0.15 |

**Quenching (4):**

| Key | Reagent | Target | Chemistry Class | CAS | Default [C] (mol/m3) |
|---|---|---|---|---|---|
| `ethanolamine_quench` | Ethanolamine | EPOXIDE | epoxide_amine | 141-43-5 | 1000 |
| `mercaptoethanol_quench` | 2-Mercaptoethanol | VINYL_SULFONE | vs_thiol | 60-24-2 | 100 |
| `nabh4_quench` | Sodium borohydride | ALDEHYDE | reduction | 16940-66-2 | 50 |
| `acetic_anhydride_quench` | Acetic anhydride | AMINE_PRIMARY | acetylation | 108-24-7 | 500 |

---

## 4. ODE Wrappers with Structured Diagnostics (Phase 3)

### 4.1 CouplingResult Dataclass (Resolves F9)

```python
@dataclass
class CouplingResult:
    """Structured output from coupling/quenching ODE solvers."""
    conversion: float                  # [0,1] fraction of target sites consumed
    sites_consumed: float              # [mol/particle] absolute
    sites_hydrolyzed: float = 0.0      # [mol/particle] lost to hydrolysis
    sites_coupled: float = 0.0         # [mol/particle] successfully coupled
    sites_blocked: float = 0.0         # [mol/particle] quenched
    reagent_remaining_fraction: float = 0.0  # [0,1]
    solver_success: bool = True
    solver_message: str = ""
    site_balance_error: float = 0.0    # Residual of conservation check
    warnings: list[str] = field(default_factory=list)
```

### 4.2 Wrapper: solve_competitive_coupling() (Template 2)

Wraps `_competitive_hydrolysis_rhs` with:
- Arrhenius rate constant computation
- pH validity check (warn if outside profile's ph_min/ph_max)
- Returns CouplingResult with hydrolyzed_sites explicitly tracked
- Hydrolysis fraction = `y_final[3] / (y_final[2] + y_final[3])` from Template 2 states
- Unit conversion: mol/particle <-> mol/m3 via bead_volume

### 4.3 Wrapper: solve_steric_coupling() (Template 3)

Wraps `_steric_binding_rhs` with:
- **Unit fix (F6):** `max_sites` converted to mol/m3: `max_surface_density * ligand_accessible_area / bead_volume`
- Activity retention applied: `sites_functional = sites_coupled * activity_retention`
- Returns CouplingResult

### 4.4 Wrapper: solve_quenching() (existing Template 1)

Wraps `solve_second_order_consumption` with:
- Operates on `remaining_activated` (not `remaining_sites`) — fixes F2
- Updates ONLY `blocked_sites` (not `consumed_sites`) — fixes F1
- Returns CouplingResult with `sites_blocked = conversion * remaining_activated`

---

## 5. Three New Solver Functions (Phases 4-5)

### 5.1 _solve_ligand_coupling_step()

**Operates on:** Activated product profile (EPOXIDE, VINYL_SULFONE, ALDEHYDE)
**Uses:** `solve_competitive_coupling()` (Template 2 wrapper)
**ACS updates:**
```python
result = solve_competitive_coupling(...)
target_profile.ligand_coupled_sites += result.sites_coupled
target_profile.ligand_functional_sites += result.sites_coupled  # activity_retention = 1.0 for small molecules
target_profile.hydrolyzed_sites += result.sites_hydrolyzed
delta_G = 0.0
```

### 5.2 _solve_protein_coupling_step()

**Operates on:** Activated product profile (EPOXIDE, VINYL_SULFONE)
**Uses:** `solve_steric_coupling()` (Template 3 wrapper)
**ACS updates:**
```python
result = solve_steric_coupling(...)
target_profile.ligand_coupled_sites += result.sites_coupled
target_profile.ligand_functional_sites += result.sites_coupled * reagent.activity_retention
target_profile.hydrolyzed_sites += result.sites_hydrolyzed  # 0 for steric model
delta_G = 0.0
```
**Trust label:** `"ranking_only"` unless user provides calibrated activity_retention.

### 5.3 _solve_quenching_step()

**Operates on:** Activated product profile (EPOXIDE, VINYL_SULFONE, ALDEHYDE, AMINE_PRIMARY)
**Uses:** `solve_quenching()` (Template 1 wrapper)
**ACS updates:**
```python
result = solve_quenching(...)
target_profile.blocked_sites += result.sites_blocked
# NO update to crosslinked_sites, ligand_coupled_sites, or G_DN
delta_G = 0.0
```

---

## 6. Backend Workflow Validation (Phase 6 — Resolves F8)

### 6.1 Validation in ModificationOrchestrator.run()

Before dispatching each step, validate:

```python
def _validate_step(self, step, acs_state, step_index, all_steps):
    errors = []
    warnings = []

    # Rule 1: Coupling/quenching requires activated sites
    if step.step_type in (LIGAND_COUPLING, PROTEIN_COUPLING, QUENCHING):
        target = acs_state.get(step.target_acs)
        if target is None or target.remaining_activated <= 0:
            errors.append(f"Step {step_index+1}: {step.step_type.value} requires "
                         f"activated {step.target_acs.value} sites, but none available.")

    # Rule 2: No steps after quenching on same target
    prior_quench_targets = {s.target_acs for s in all_steps[:step_index]
                           if s.step_type == QUENCHING}
    if step.target_acs in prior_quench_targets:
        errors.append(f"Step {step_index+1}: {step.target_acs.value} was already quenched.")

    # Rule 3: Quenching should be last for its target type
    if step.step_type == QUENCHING:
        later_same_target = [s for s in all_steps[step_index+1:]
                            if s.target_acs == step.target_acs]
        if later_same_target:
            warnings.append(f"Step {step_index+1}: Quenching {step.target_acs.value} "
                          f"but {len(later_same_target)} later step(s) target the same site.")

    # Rule 4: Reagent-target compatibility
    # (reagent_profile.target_acs must match step.target_acs)

    # Rule 5: pH/temperature validity windows
    # (warn if step.ph outside reagent.ph_min/ph_max)

    return errors, warnings
```

### 6.2 UI Validators (Supplementary)

Keep existing UI validators in `ui_validators.py` as preflight UX checks. Backend validation is authoritative.

---

## 7. FunctionalMediaContract: M2→M3 Bridge (Phase 7 — Resolves F7)

### 7.1 New Dataclass

```python
@dataclass
class FunctionalMediaContract:
    """Stable interface between Module 2 and Module 3 for process simulation."""
    # ── From M1 (pass-through) ──
    bead_d50: float                    # [m]
    porosity: float                    # [-]
    pore_size_mean: float              # [m]

    # ── From M2 ──
    ligand_type: str                   # "iex_anion", "iex_cation", "affinity",
                                       # "imac", "hic", "none"
    installed_ligand: str              # e.g., "DEAE", "Protein A", "Phenyl"
    functional_ligand_density: float   # [mol/m^2] ligand_functional_sites / area
    total_coupled_density: float       # [mol/m^2] ligand_coupled_sites / area
    charge_density: float = 0.0        # [mol/m^2] for IEX (= functional_ligand_density)
    active_protein_density: float = 0.0  # [mol/m^2] for affinity
    G_DN_updated: float = 0.0         # [Pa]
    E_star_updated: float = 0.0       # [Pa]

    # ── Derived M3 parameter estimates ──
    estimated_q_max: float = 0.0       # [mol/m^3 bed] mapped from ligand density
    q_max_confidence: str = "not_mapped"  # "mapped_calibrated" | "mapped_estimated" | "not_mapped"
    q_max_mapping_notes: str = ""

    # ── Trust ──
    confidence_tier: str = "semi_quantitative"
    warnings: list[str] = field(default_factory=list)
```

### 7.2 Mapping Logic

For IEX:
```python
# Ionic capacity = functional_ligand_density * accessible_area_per_bed_volume
q_max_estimated = functional_ligand_density * specific_surface_area * (1 - bed_porosity)
```

For Affinity:
```python
# Active binding capacity = active_protein_density * ligand_accessible_area_per_bed_volume
# Stoichiometry: 1 Protein A binds 2 IgG Fc
q_max_estimated = active_protein_density * ligand_accessible_area_per_bed_vol * binding_stoichiometry
```

For HIC/IMAC: Label as `"not_mapped"` until salt/metal-dependent models are added.

### 7.3 M3 Integration

M3 should:
1. Check if `FunctionalMediaContract` exists in session_state
2. If `q_max_confidence != "not_mapped"`, use `estimated_q_max` as default (user can override)
3. Display confidence tier and mapping notes in M3 tab

---

## 8. Revised Implementation Phases

### Phase Order (Audit-Corrected)

| Phase | Module | Audit Finding Resolved | Model Tier | Est. LOC | Depends On |
|---|---|---|---|---|---|
| **1** | ACS state model revision | F1, F2 | Opus | 120 | — |
| **2** | ReagentProfile expansion + 10 profiles | F4 | Sonnet | 200 | — |
| **3** | ODE wrappers with diagnostics + unit fix | F3, F6, F9 | Opus | 120 | Phase 1 |
| **4** | Ligand coupling solver | F3 | Sonnet | 100 | Phases 1, 2, 3 |
| **5** | Protein coupling solver (ranking_only) | F5 | Opus | 100 | Phases 1, 2, 3 |
| **6** | Quenching solver + backend validation | F1, F8 | Sonnet | 120 | Phases 1, 2, 3 |
| **7** | FunctionalMediaContract (M2→M3 bridge) | F7 | Opus | 100 | Phases 4, 5, 6 |
| **8** | Dispatch integration + existing workflow migration | — | Sonnet | 60 | All above |
| **9** | UI integration (M2 tab expansion) | — | Sonnet | 50 | Phase 8 |
| **10** | Conservation + integration testing | All | Opus | 100 | All |
| **Total** | | | | **~1,070** | |

### Phase 1: ACS State Model Revision (FOUNDATION — Do First)

**Files:** `acs.py`
**Changes:**
1. Add `crosslinked_sites: float = 0.0` and `hydrolyzed_sites: float = 0.0`
2. Remove `consumed_sites` field
3. Update `remaining_sites` property to use terminal-sum formula
4. Add `remaining_activated` property
5. Update `validate()` with new invariant
6. Update `initialize_acs_from_m1()` — no changes needed (all terminal states start at 0)

**Migration:** Update `_solve_crosslinking_step` and `_solve_activation_step` to use `crosslinked_sites` instead of `consumed_sites`.

**Tests:**
- Conservation holds after crosslinking at 50% conversion
- Conservation holds after activation at 80% conversion
- `remaining_activated = activated_sites` for fresh product profile
- `remaining_activated = 0` after full coupling + quenching
- Negative field detection

### Audit-Recommended First Validated Path (Phase 6 Exit Criterion)

```
ECH activation → ethanolamine quench
```

This path exercises:
- Activation (existing, migrated to new ACS model)
- Quenching (new, operates on `remaining_activated`, updates `blocked_sites` only)
- Conservation invariant across the full sequence
- Backend validation (quenching requires activated sites)

### Audit-Recommended First Functional Ligand Path (Phase 8 Exit Criterion)

```
ECH activation → DEAE ligand coupling → ethanolamine quench
```

This exercises the full Activation → Coupling → Quenching sequence with:
- Hydrolysis tracking (competitive coupling)
- `ligand_coupled_sites` and `ligand_functional_sites` updated
- Conservation after all three steps
- `FunctionalMediaContract` generation with IEX q_max estimate

---

## 9. Acceptance Criteria (From Audit Section 8)

### Scientific Acceptance

- [ ] ACS conservation holds for ALL workflows including >95% quenching
- [ ] `hydrolyzed_sites` tracked explicitly for coupling workflows
- [ ] pH and temperature validated against reagent profile windows
- [ ] Every new reagent profile has verified CAS, chemistry_class, and confidence_tier
- [ ] Protein coupling labeled `"ranking_only"` with adjustable activity_retention
- [ ] M2 functional ligand density maps to M3 `estimated_q_max` for IEX and affinity
- [ ] UI does not allow unsupported workflows to run

### Computational Acceptance

- [ ] Unit tests for all 10 new profiles (instantiation + field validation)
- [ ] Unit tests for ODE wrappers (zero conc, high conc, hydrolysis-dominated, steric-limited)
- [ ] Conservation tests after each step type and full workflow sequences
- [ ] Backend ordering validation tests independent of Streamlit
- [ ] M1→M2→M3 integration tests with parameter propagation verification
- [ ] Solver diagnostics captured in CouplingResult and surfaced in ModificationResult.notes
- [ ] 95% quenching conservation regression test (audit F1 blocker)

---

## 10. What Is Explicitly Deferred

Per audit recommendation, the following are **out of scope** for this phase:

| Feature | Reason | When |
|---|---|---|
| Mechanistic pH model (pKa-dependent rates) | Over-engineering for screening tool | v6.0 if calibration data available |
| Metal charging step for IMAC | Requires new step type + metal inventory | v5.2 |
| Pore-size distribution (vs mean pore) | Requires L2 extension | v6.0 |
| Spacer-arm / orientation modeling for proteins | Insufficient published parameters | v6.0 |
| Protein leaching / deactivation kinetics | Long-term stability model | v6.0 |
| Buffer/ionic-strength effects | Complex multi-equilibrium | v6.0 |
| Salt-dependent HIC retention | Requires Solvophobic theory | v5.2 |
| Hazard management system | Not a simulation concern | External SDS database |

---

> **Disclaimer**: Rate constants are order-of-magnitude estimates from general reaction
> chemistry. All values are illustrative defaults requiring experimental calibration.
> Protein coupling outputs are ranking-only unless calibrated against specific materials.
