"""Modification step definitions and workflow solvers for Module 2.

Phase B: Minimal Module 2 — 2 Workflows.
Architecture: module2_module3_final_implementation_plan.md, Phase B.

Two validated workflows:

Workflow 1 — Amine secondary crosslinking:
    Target: AMINE_PRIMARY sites.
    Uses solve_second_order_consumption with reagent's k and E_a.
    Updates: crosslinked_sites increases, remaining decreases.
    If crosslinking: updates G_DN estimate (additive rubber elasticity).

Workflow 2 — Hydroxyl activation (ECH/DVS):
    Target: HYDROXYL sites.
    Product: EPOXIDE (ECH) or VINYL_SULFONE (DVS).
    Creates new ACSProfile entries for the product site type.
    Includes alkaline hydrolysis side reaction (k_hydrol for ECH ~1e-4 /s at pH 12).
"""

from __future__ import annotations

import copy
import logging
import math
from dataclasses import dataclass, field
from enum import Enum
from typing import Optional

from ..datatypes import ModelEvidenceTier, ModelManifest
from .acs import ACSSiteType, ACSProfile
from .reactions import (
    solve_second_order_consumption,
    solve_competitive_coupling,
    solve_steric_coupling,
    solve_quenching,
)
from .reagent_profiles import ReagentProfile
from .surface_area import AccessibleSurfaceModel

logger = logging.getLogger(__name__)

# ─── Constants ──────────────────────────────────────────────────────────

K_BOLTZMANN = 1.38e-23  # [J/K]
N_AVOGADRO = 6.022e23
R_GAS = 8.314  # J/(mol*K)


# ─── Arrhenius helpers ──────────────────────────────────────────────────

def _arrhenius_prefactor(k_ref: float, E_a: float, T_ref: float) -> float:
    """Back-calculate Arrhenius prefactor from reference-condition rate constant.

    k_ref = k0 * exp(-E_a/(R*T_ref))  =>  k0 = k_ref * exp(E_a/(R*T_ref))

    Args:
        k_ref: Rate constant at reference temperature [m^3/(mol*s)].
        E_a: Activation energy [J/mol].
        T_ref: Reference temperature [K].

    Returns:
        Arrhenius prefactor k0 [m^3/(mol*s)], or 0.0 if inputs are invalid.
    """
    if E_a <= 0 or T_ref <= 0 or k_ref <= 0:
        return 0.0
    # Guard overflow for large E_a (Codex R4 fix: continuous cap)
    exponent = E_a / (R_GAS * T_ref)
    if exponent > 700:
        return 1e20
    result = k_ref * math.exp(exponent)
    return min(result, 1e20)


# ─── Enums ──────────────────────────────────────────────────────────────

class ModificationStepType(Enum):
    """Type of surface modification step."""
    SECONDARY_CROSSLINKING = "secondary_crosslinking"
    ACTIVATION = "activation"
    LIGAND_COUPLING = "ligand_coupling"
    PROTEIN_COUPLING = "protein_coupling"
    QUENCHING = "quenching"
    SPACER_ARM = "spacer_arm"  # v5.8: consumes activated site, creates intermediate ACS
    METAL_CHARGING = "metal_charging"      # v5.9.1: loads metal onto IMAC chelator
    PROTEIN_PRETREATMENT = "protein_pretreatment"  # v5.9.2: reduces protein disulfides
    WASHING = "washing"                    # v5.9.4: removes residual reagents


# ─── Data classes ──────────────────────────────────────────────────────

@dataclass
class ModificationStep:
    """Specification for a single modification step.

    Attributes:
        step_type: Type of modification chemistry.
        reagent_key: Key into REAGENT_PROFILES dictionary.
        target_acs: ACS site type consumed by this step.
        product_acs: ACS site type produced (None for crosslinking/quenching).
        temperature: Reaction temperature [K].
        time: Reaction time [s].
        ph: Reaction pH [-].
        reagent_concentration: Reagent concentration [mol/m^3].
        stoichiometry: Moles reagent per mole ACS consumed [-].
        conversion: Populated after solving — fraction of target ACS consumed.
    """
    step_type: ModificationStepType
    reagent_key: str
    target_acs: ACSSiteType
    product_acs: Optional[ACSSiteType] = None
    temperature: float = 298.15       # [K]
    time: float = 3600.0              # [s]
    ph: float = 7.0
    reagent_concentration: float = 10.0  # [mol/m^3]
    stoichiometry: float = 1.0
    conversion: float = 0.0           # populated after solving


@dataclass
class ModificationResult:
    """Result of a single modification step.

    Attributes:
        step: The specification that was executed.
        acs_before: ACS state snapshot before this step.
        acs_after: ACS state snapshot after this step.
        conversion: Fraction of target ACS consumed, in [0, 1].
        delta_G_DN: Change in G_DN from crosslinking [Pa].
        notes: Diagnostic notes (trust warnings, etc.).
    """
    step: ModificationStep
    acs_before: dict[ACSSiteType, ACSProfile]
    acs_after: dict[ACSSiteType, ACSProfile]
    conversion: float = 0.0
    delta_G_DN: float = 0.0
    notes: str = ""
    # v6.1 (Node 4): evidence provenance — populated by solve_modification_step.
    # Tier defaults to SEMI_QUANTITATIVE for M2 chemistry (rate constants and
    # accessibility coefficients are literature/profile-derived, not locally
    # calibrated). Downgrade to QUALITATIVE_TREND when the reagent profile is
    # ranking-only or when conservation violations are detected.
    model_manifest: Optional[ModelManifest] = None


# ─── Core solver ────────────────────────────────────────────────────────

def solve_modification_step(
    step: ModificationStep,
    acs_state: dict[ACSSiteType, ACSProfile],
    surface_model: AccessibleSurfaceModel,
    reagent_profile: ReagentProfile,
) -> ModificationResult:
    """Execute one modification step, updating ACS state.

    Dispatches to the appropriate workflow based on step_type:
      - SECONDARY_CROSSLINKING: amine crosslinking workflow
      - ACTIVATION: hydroxyl activation workflow

    Args:
        step: Modification step specification.
        acs_state: Current ACS inventory (mutated in place).
        surface_model: Pre-computed surface area model.
        reagent_profile: Reagent kinetic profile.

    Returns:
        ModificationResult with before/after ACS snapshots and conversion.

    Raises:
        ValueError: If target ACS type is not present in acs_state.
        KeyError: If step_type is not implemented.
    """
    # --- Snapshot before ---
    acs_before = _snapshot_acs(acs_state)

    # --- Validate reagent compatibility BEFORE any early returns (Codex R3-F1) ---
    _STEP_ALLOWED_RTYPES = {
        ModificationStepType.SECONDARY_CROSSLINKING: {"crosslinking"},
        ModificationStepType.ACTIVATION: {"activation"},
        ModificationStepType.LIGAND_COUPLING: {"coupling"},
        ModificationStepType.PROTEIN_COUPLING: {"protein_coupling"},
        ModificationStepType.QUENCHING: {"blocking"},
        ModificationStepType.SPACER_ARM: {"spacer_arm", "heterobifunctional"},
        ModificationStepType.METAL_CHARGING: {"metal_charging", "metal_stripping"},
        ModificationStepType.PROTEIN_PRETREATMENT: {"protein_pretreatment"},
        ModificationStepType.WASHING: {"washing"},
    }
    allowed = _STEP_ALLOWED_RTYPES.get(step.step_type, set())
    if allowed and reagent_profile.reaction_type not in allowed:
        raise ValueError(
            f"Reagent '{step.reagent_key}' has reaction_type "
            f"'{reagent_profile.reaction_type}' incompatible with step type "
            f"'{step.step_type.value}'. Allowed: {sorted(allowed)}."
        )
    # Also check target_acs match (mirrors orchestrator Rule 3)
    if reagent_profile.target_acs != step.target_acs:
        raise ValueError(
            f"Reagent '{step.reagent_key}' targets "
            f"{reagent_profile.target_acs.value} but step targets "
            f"{step.target_acs.value}. Reagent-target mismatch."
        )

    # --- Validate target exists ---
    if step.target_acs not in acs_state:
        raise ValueError(
            f"Target ACS type {step.target_acs.value} not present in "
            f"current state. Available: {[s.value for s in acs_state.keys()]}"
        )

    target_profile = acs_state[step.target_acs]
    remaining = target_profile.remaining_sites

    if remaining <= 0:
        logger.warning(
            "No remaining %s sites for step %s. Skipping.",
            step.target_acs.value, step.reagent_key,
        )
        return ModificationResult(
            step=step,
            acs_before=acs_before,
            acs_after=_snapshot_acs(acs_state),
            conversion=0.0,
            notes="No remaining sites; step skipped.",
        )

    # --- Convert per-particle sites to bulk concentration ---
    V_bead = surface_model.bead_volume  # [m^3/particle]

    if V_bead <= 0:
        raise ValueError(f"Bead volume must be positive, got {V_bead}")

    acs_concentration = remaining / V_bead  # [mol/m^3]

    # --- Dispatch by step type ---
    if step.step_type == ModificationStepType.SECONDARY_CROSSLINKING:
        result = _solve_crosslinking_step(
            step, acs_state, target_profile, acs_concentration,
            reagent_profile, V_bead,
        )
    elif step.step_type == ModificationStepType.ACTIVATION:
        result = _solve_activation_step(
            step, acs_state, target_profile, acs_concentration,
            reagent_profile, V_bead,
        )
    elif step.step_type == ModificationStepType.LIGAND_COUPLING:
        result = _solve_ligand_coupling_step(
            step, acs_state, target_profile,
            reagent_profile, surface_model, V_bead,
        )
    elif step.step_type == ModificationStepType.PROTEIN_COUPLING:
        result = _solve_protein_coupling_step(
            step, acs_state, target_profile,
            reagent_profile, surface_model, V_bead,
        )
    elif step.step_type == ModificationStepType.QUENCHING:
        result = _solve_quenching_step(
            step, acs_state, target_profile,
            reagent_profile, V_bead,
        )
    elif step.step_type == ModificationStepType.SPACER_ARM:
        result = _solve_spacer_arm_step(
            step, acs_state, target_profile,
            reagent_profile, surface_model, V_bead,
        )
    elif step.step_type == ModificationStepType.METAL_CHARGING:
        result = _solve_metal_charging_step(
            step, acs_state, target_profile, reagent_profile,
        )
    elif step.step_type == ModificationStepType.PROTEIN_PRETREATMENT:
        result = _solve_protein_pretreatment_step(
            step, reagent_profile,
        )
    elif step.step_type == ModificationStepType.WASHING:
        result = _solve_washing_step(
            step, acs_state, surface_model, V_bead,
        )
    else:
        raise KeyError(
            f"Step type {step.step_type.value} not implemented. "
            f"Supported: all 9 ModificationStepType values."
        )

    # --- Snapshot after ---
    result.acs_before = acs_before
    result.acs_after = _snapshot_acs(acs_state)

    # --- Validate conservation after step ---
    conservation_violations: list[str] = []
    for profile in acs_state.values():
        violations = profile.validate()
        if violations:
            conservation_violations.extend(violations)
            result.notes += f" Conservation violations: {'; '.join(violations)}"
            logger.warning("ACS conservation violation: %s", violations)

    # --- Attach model manifest (Node 4: M2 evidence wiring) ---
    result.model_manifest = _build_step_manifest(
        step=step,
        reagent_profile=reagent_profile,
        result=result,
        conservation_violations=conservation_violations,
    )

    return result


def _build_step_manifest(
    step: ModificationStep,
    reagent_profile: ReagentProfile,
    result: ModificationResult,
    conservation_violations: list[str],
) -> ModelManifest:
    """Construct the per-step ModelManifest.

    Tier policy (Node 4):
      - UNSUPPORTED: solver dispatched but produced no chemistry (skipped).
      - QUALITATIVE_TREND: ranking-only reagent (affinity/biotin/heparin) OR
        any conservation violation detected after the step.
      - SEMI_QUANTITATIVE: literature/profile-derived rate constants without
        local calibration. This is the M2 default in v6.1; future calibration
        through CalibrationStore (Sprint 3 of consensus plan) can upgrade
        this to CALIBRATED_LOCAL or VALIDATED_QUANTITATIVE.

    Diagnostics carry the achieved conversion, delta_G_DN, and any conservation
    violations so RunReport.compute_min_tier() can roll the M2 evidence into a
    pipeline-wide tier without re-reading the reagent profile.
    """
    # Default to SEMI_QUANTITATIVE — the v6.1 default for M2 chemistry.
    tier = ModelEvidenceTier.SEMI_QUANTITATIVE

    # Ranking-only reagent classes cap evidence at QUALITATIVE_TREND. The
    # functional_mode strings come from reagent_profiles.py and indicate
    # binding modes whose q_max is intrinsically not quantitatively predictable
    # from ligand density alone (e.g., affinity stoichiometry depends on the
    # target protein, biotin/heparin involve cooperative binding effects).
    _ranking_modes = {"affinity_ligand", "biotin_affinity", "heparin_affinity"}
    if getattr(reagent_profile, "functional_mode", "") in _ranking_modes:
        tier = ModelEvidenceTier.QUALITATIVE_TREND

    # Conservation violations override everything: a step that broke ACS
    # bookkeeping cannot back a quantitative claim regardless of reagent class.
    if conservation_violations:
        tier = ModelEvidenceTier.QUALITATIVE_TREND

    # No chemistry happened — surface had nothing to react with.
    if result.conversion <= 0.0 and "skipped" in result.notes.lower():
        tier = ModelEvidenceTier.UNSUPPORTED

    diagnostics: dict = {
        "conversion": float(result.conversion),
        "delta_G_DN_Pa": float(result.delta_G_DN),
        "conservation_ok": len(conservation_violations) == 0,
    }
    if conservation_violations:
        diagnostics["conservation_violations"] = list(conservation_violations)

    assumptions: list[str] = []
    if getattr(reagent_profile, "k_ref", None) is not None:
        assumptions.append(
            f"Rate constant k_ref={reagent_profile.k_ref:g} from reagent profile "
            f"(literature value, not locally calibrated)."
        )
    if getattr(reagent_profile, "activity_retention", None) is not None:
        assumptions.append(
            f"activity_retention={reagent_profile.activity_retention:.2f} default "
            "(per-protein value not calibrated)."
        )

    valid_domain: dict = {}
    # Reagent profiles carry pH/T validity ranges in some cases; surface them
    # here if present so trust gates can warn when the step runs outside them.
    for attr, key in (("pH_range", "pH"), ("T_range_K", "T_K")):
        rng = getattr(reagent_profile, attr, None)
        if rng is not None:
            valid_domain[key] = tuple(rng) if isinstance(rng, (list, tuple)) else rng

    return ModelManifest(
        model_name=f"M2.{step.step_type.value}.{step.reagent_key}",
        evidence_tier=tier,
        valid_domain=valid_domain,
        calibration_ref="",
        assumptions=assumptions,
        diagnostics=diagnostics,
    )


# ─── Workflow 1: Amine secondary crosslinking ─────────────────────────

def _solve_crosslinking_step(
    step: ModificationStep,
    acs_state: dict[ACSSiteType, ACSProfile],
    target_profile: ACSProfile,
    acs_concentration: float,
    reagent_profile: ReagentProfile,
    V_bead: float,
) -> ModificationResult:
    """Amine secondary crosslinking workflow.

    Consumes AMINE_PRIMARY sites, increases crosslinked_sites count.
    Estimates delta_G_DN from additive rubber elasticity:
        delta_G = delta_nu_e * kB * T
    where delta_nu_e is the additional effective crosslink density.
    """
    k0 = _arrhenius_prefactor(
        reagent_profile.k_forward, reagent_profile.E_a, reagent_profile.temperature_default
    )
    conversion, reagent_remaining = solve_second_order_consumption(
        acs_concentration=acs_concentration,
        reagent_concentration=step.reagent_concentration,
        k_forward=reagent_profile.k_forward,
        stoichiometry=step.stoichiometry,
        time=step.time,
        temperature=step.temperature,
        E_a=reagent_profile.E_a,
        k0=k0,
        hydrolysis_rate=reagent_profile.hydrolysis_rate,
    )

    # Update ACS state — crosslinking consumes sites into crosslinked terminal state
    sites_consumed = conversion * target_profile.remaining_sites
    target_profile.crosslinked_sites += sites_consumed

    step.conversion = conversion

    # Estimate delta G_DN from rubber elasticity
    # Each crosslink bridge consumes 2 amine groups (stoich=0.5 means
    # 0.5 mol reagent per mol NH2, so 1 bridge per 2 NH2).
    # New crosslink density: sites_consumed / (2 * V_bead) * N_A [1/m^3]
    # G = nu_e * kB * T
    bridges_per_particle = sites_consumed / 2.0  # [mol/particle]
    nu_e_new = bridges_per_particle * N_AVOGADRO / V_bead  # [1/m^3]
    # Assume f_bridge ~ 0.4 for secondary crosslinking (same as L3)
    f_bridge = 0.4
    delta_G = f_bridge * nu_e_new * K_BOLTZMANN * step.temperature

    logger.info(
        "Crosslinking step '%s': conversion=%.4f, delta_G=%.1f Pa",
        step.reagent_key, conversion, delta_G,
    )

    return ModificationResult(
        step=step,
        acs_before={},  # populated by caller
        acs_after={},
        conversion=conversion,
        delta_G_DN=delta_G,
        notes=f"Crosslinking: {sites_consumed:.4e} mol NH2 consumed, "
              f"delta_G={delta_G:.1f} Pa",
    )


# ─── Workflow 2: Hydroxyl activation (ECH/DVS) ────────────────────────

def _solve_activation_step(
    step: ModificationStep,
    acs_state: dict[ACSSiteType, ACSProfile],
    target_profile: ACSProfile,
    acs_concentration: float,
    reagent_profile: ReagentProfile,
    V_bead: float,
) -> ModificationResult:
    """Hydroxyl activation workflow.

    Consumes HYDROXYL sites, produces EPOXIDE or VINYL_SULFONE sites.
    Includes hydrolysis side reaction for ECH (k_hydrol > 0).
    """
    k0 = _arrhenius_prefactor(
        reagent_profile.k_forward, reagent_profile.E_a, reagent_profile.temperature_default
    )
    conversion, reagent_remaining = solve_second_order_consumption(
        acs_concentration=acs_concentration,
        reagent_concentration=step.reagent_concentration,
        k_forward=reagent_profile.k_forward,
        stoichiometry=step.stoichiometry,
        time=step.time,
        temperature=step.temperature,
        E_a=reagent_profile.E_a,
        k0=k0,
        hydrolysis_rate=reagent_profile.hydrolysis_rate,
    )

    # Sites consumed from target (hydroxyl) — activation consumes into activated_consumed
    # (Codex P2-1: activation is NOT crosslinking, does not contribute to G_DN)
    # NOTE (Codex F4): The same sites appear as terminal on the parent (OH) profile
    # AND as new accessible/activated on the product (EPOXIDE) profile. This is by
    # design — each ACS profile tracks conservation independently. A global sum across
    # profiles will double-count, which is correct: the OH site is consumed, and a new
    # EPOXIDE site is created. They are chemically distinct inventory entries.
    sites_consumed = conversion * target_profile.remaining_sites
    target_profile.activated_consumed_sites += sites_consumed

    step.conversion = conversion

    # Compute how much reagent was lost to hydrolysis vs coupling
    # reagent consumed = initial - remaining
    reagent_consumed_total = step.reagent_concentration * (1.0 - reagent_remaining)
    # Of the reagent consumed, some went to coupling (stoich * sites_consumed / V_bead)
    # and the rest to hydrolysis.
    reagent_for_coupling = step.stoichiometry * sites_consumed / V_bead
    reagent_hydrolyzed = max(reagent_consumed_total - reagent_for_coupling, 0.0)

    # Product sites: coupling produces activated sites on the surface
    product_type = step.product_acs or reagent_profile.product_acs
    if product_type is not None:
        product_sites = sites_consumed  # 1:1 stoichiometry for activation

        if product_type in acs_state:
            # Add to existing product profile
            acs_state[product_type].accessible_sites += product_sites
            acs_state[product_type].total_sites += product_sites
            acs_state[product_type].activated_sites += product_sites
        else:
            # Create new product profile
            acs_state[product_type] = ACSProfile(
                site_type=product_type,
                total_sites=product_sites,
                accessible_sites=product_sites,
                activated_sites=product_sites,
                crosslinked_sites=0.0,
                activated_consumed_sites=0.0,
                hydrolyzed_sites=0.0,
                ligand_coupled_sites=0.0,
                ligand_functional_sites=0.0,
                blocked_sites=0.0,
                total_density=0.0,
                accessible_density=0.0,
                accessibility_model=target_profile.accessibility_model,
                uncertainty_fraction=target_profile.uncertainty_fraction,
            )

    hydrol_note = ""
    if reagent_profile.hydrolysis_rate > 0:
        hydrol_note = (
            f" Hydrolysis: ~{reagent_hydrolyzed:.4e} mol/m^3 reagent lost."
        )

    logger.info(
        "Activation step '%s': conversion=%.4f, product=%s%s",
        step.reagent_key, conversion,
        product_type.value if product_type else "none",
        hydrol_note,
    )

    return ModificationResult(
        step=step,
        acs_before={},
        acs_after={},
        conversion=conversion,
        delta_G_DN=0.0,
        notes=f"Activation: {sites_consumed:.4e} mol OH consumed, "
              f"produced {product_type.value if product_type else 'none'}.{hydrol_note}",
    )


# ─── Workflow 3: Ligand Coupling (W5) ────────────────────────────────

def _solve_ligand_coupling_step(
    step: ModificationStep,
    acs_state: dict[ACSSiteType, ACSProfile],
    target_profile: ACSProfile,
    reagent_profile: ReagentProfile,
    surface_model: AccessibleSurfaceModel,
    V_bead: float,
) -> ModificationResult:
    """Couple small-molecule ligand to activated sites with competitive hydrolysis.

    Operates on remaining_activated sites (not remaining_sites).
    Uses solve_competitive_coupling (Template 2 wrapper).
    Updates: ligand_coupled_sites, ligand_functional_sites, hydrolyzed_sites.
    No G_DN change.
    """
    remaining_act = target_profile.remaining_activated
    if remaining_act <= 0:
        step.conversion = 0.0
        return ModificationResult(
            step=step, acs_before={}, acs_after={},
            conversion=0.0, delta_G_DN=0.0,
            notes="Ligand coupling: no remaining activated sites.",
        )

    cr = solve_competitive_coupling(
        activated_sites_mol_per_particle=remaining_act,
        ligand_concentration=step.reagent_concentration,
        k_forward=reagent_profile.k_forward,
        E_a=reagent_profile.E_a,
        temperature=step.temperature,
        time=step.time,
        bead_volume=V_bead,
        hydrolysis_rate=reagent_profile.hydrolysis_rate,
        stoichiometry=step.stoichiometry,
        ph=step.ph,
        ph_min=getattr(reagent_profile, 'ph_min', 0.0),
        ph_max=getattr(reagent_profile, 'ph_max', 14.0),
    )

    # Update ACS — ligand coupling terminal states
    # Apply spacer multiplier on activity_retention (audit F2, capped at 1.0)
    _activity_ret = getattr(reagent_profile, 'activity_retention', 1.0)
    _spacer_mult = getattr(reagent_profile, 'spacer_activity_multiplier', 1.0)
    _effective_activity = min(_activity_ret * _spacer_mult, 1.0)
    target_profile.ligand_coupled_sites += cr.sites_coupled
    target_profile.ligand_functional_sites += cr.sites_coupled * _effective_activity
    target_profile.hydrolyzed_sites += cr.sites_hydrolyzed

    step.conversion = cr.conversion

    notes = (
        f"Ligand coupling '{step.reagent_key}': conv={cr.conversion:.4f}, "
        f"coupled={cr.sites_coupled:.4e}, hydrolyzed={cr.sites_hydrolyzed:.4e}"
    )
    if cr.warnings:
        notes += f" Warnings: {'; '.join(cr.warnings)}"

    logger.info(notes)

    return ModificationResult(
        step=step, acs_before={}, acs_after={},
        conversion=cr.conversion, delta_G_DN=0.0, notes=notes,
    )


# ─── Workflow 4: Protein Coupling (W6) ───────────────────────────────

def _solve_protein_coupling_step(
    step: ModificationStep,
    acs_state: dict[ACSSiteType, ACSProfile],
    target_profile: ACSProfile,
    reagent_profile: ReagentProfile,
    surface_model: AccessibleSurfaceModel,
    V_bead: float,
) -> ModificationResult:
    """Couple protein ligand to activated sites with steric blocking.

    Operates on remaining_activated sites. Uses solve_steric_coupling (Template 3).
    Steric limit from reagent.max_surface_density * ligand_accessible_area.
    Updates: ligand_coupled_sites, ligand_functional_sites (x activity_retention).
    Trust label: ranking_only.
    """
    # v5.8 WN-3: Maleimide decay — apply first-order hydrolysis before coupling
    # if the target is MALEIMIDE, immobilized maleimide rings open over time
    if target_profile.site_type == ACSSiteType.MALEIMIDE:
        _k_decay = getattr(reagent_profile, 'maleimide_decay_rate', 0.0)
        if _k_decay <= 0:
            # Try to get decay rate from the MALEIMIDE profile itself (set by SM(PEG)n)
            _k_decay = 1e-5  # ASSUMPTION: default maleimide hydrolysis at pH 7, 25C
        _avail_before = target_profile.remaining_activated
        if _avail_before > 0 and step.time > 0:
            _frac_remaining = math.exp(-_k_decay * step.time)
            _sites_decayed = _avail_before * (1.0 - _frac_remaining)
            target_profile.hydrolyzed_sites += _sites_decayed
            logger.info(
                "Maleimide decay: k=%.1e /s, t=%.0f s, decayed=%.4e mol/particle (%.1f%%)",
                _k_decay, step.time, _sites_decayed, (1 - _frac_remaining) * 100,
            )

    remaining_act = target_profile.remaining_activated
    if remaining_act <= 0:
        step.conversion = 0.0
        return ModificationResult(
            step=step, acs_before={}, acs_after={},
            conversion=0.0, delta_G_DN=0.0,
            notes="Protein coupling: no remaining activated sites (maleimide may have decayed).",
        )

    # Steric jamming limit in mol/particle
    max_density = getattr(reagent_profile, 'max_surface_density', 0.0)
    if max_density > 0 and surface_model.ligand_accessible_area > 0:
        max_coupled = max_density * surface_model.ligand_accessible_area
    else:
        max_coupled = remaining_act  # no steric limit if not specified

    # Account for already-coupled sites
    already_coupled = target_profile.ligand_coupled_sites
    max_coupled = max(max_coupled - already_coupled, 0.0)

    cr = solve_steric_coupling(
        activated_sites_mol_per_particle=remaining_act,
        ligand_concentration=step.reagent_concentration,
        k_forward=reagent_profile.k_forward,
        E_a=reagent_profile.E_a,
        temperature=step.temperature,
        time=step.time,
        bead_volume=V_bead,
        max_coupled_mol_per_particle=max_coupled,
        ph=step.ph,
        ph_min=getattr(reagent_profile, 'ph_min', 0.0),
        ph_max=getattr(reagent_profile, 'ph_max', 14.0),
    )

    # Update ACS — protein coupling with activity retention + spacer multiplier
    activity_ret = getattr(reagent_profile, 'activity_retention', 1.0)
    spacer_mult = getattr(reagent_profile, 'spacer_activity_multiplier', 1.0)
    effective_activity = min(activity_ret * spacer_mult, 1.0)  # cap at 1.0
    target_profile.ligand_coupled_sites += cr.sites_coupled
    target_profile.ligand_functional_sites += cr.sites_coupled * effective_activity

    step.conversion = cr.conversion

    notes = (
        f"Protein coupling '{step.reagent_key}': conv={cr.conversion:.4f}, "
        f"coupled={cr.sites_coupled:.4e}, functional={cr.sites_coupled * activity_ret:.4e} "
        f"(activity_retention={activity_ret:.2f}). RANKING_ONLY."
    )
    if cr.warnings:
        notes += f" Warnings: {'; '.join(cr.warnings)}"

    logger.info(notes)

    return ModificationResult(
        step=step, acs_before={}, acs_after={},
        conversion=cr.conversion, delta_G_DN=0.0, notes=notes,
    )


# ─── Workflow 5: Quenching (W7) ─────────────────────────────────────

def _solve_quenching_step(
    step: ModificationStep,
    acs_state: dict[ACSSiteType, ACSProfile],
    target_profile: ACSProfile,
    reagent_profile: ReagentProfile,
    V_bead: float,
) -> ModificationResult:
    """Block remaining sites with quenching reagent.

    For activated product profiles (EPOXIDE, VINYL_SULFONE, ALDEHYDE):
        operates on remaining_activated sites.
    For native profiles (AMINE_PRIMARY, HYDROXYL) where activated_sites == 0:
        operates on remaining_sites (Codex P1-2 fix — acetic anhydride on amines).

    Updates: blocked_sites ONLY (not crosslinked_sites — fixes audit F1).
    No G_DN change. No product sites created.
    """
    # Codex P1-2 + F5: Use remaining_activated for activated product profiles,
    # remaining_sites for explicitly native site types (AMINE_PRIMARY, HYDROXYL).
    # Other profiles with activated_sites==0 are likely malformed — warn and skip.
    _NATIVE_SITE_TYPES = {ACSSiteType.AMINE_PRIMARY, ACSSiteType.HYDROXYL}
    if target_profile.activated_sites > 0:
        available = target_profile.remaining_activated
    elif target_profile.site_type in _NATIVE_SITE_TYPES:
        available = target_profile.remaining_sites
    else:
        logger.warning(
            "Quenching targets %s with activated_sites=0 and non-native type. Skipping.",
            target_profile.site_type.value,
        )
        available = 0.0

    if available <= 0:
        step.conversion = 0.0
        return ModificationResult(
            step=step, acs_before={}, acs_after={},
            conversion=0.0, delta_G_DN=0.0,
            notes="Quenching: no remaining sites to block.",
        )

    cr = solve_quenching(
        activated_sites_mol_per_particle=available,
        reagent_concentration=step.reagent_concentration,
        k_forward=reagent_profile.k_forward,
        E_a=reagent_profile.E_a,
        temperature=step.temperature,
        time=step.time,
        bead_volume=V_bead,
        stoichiometry=step.stoichiometry,
        hydrolysis_rate=reagent_profile.hydrolysis_rate,
    )

    # Update ACS — quenching only increments blocked_sites (F1 fix)
    target_profile.blocked_sites += cr.sites_blocked

    step.conversion = cr.conversion

    notes = (
        f"Quenching '{step.reagent_key}': conv={cr.conversion:.4f}, "
        f"blocked={cr.sites_blocked:.4e}"
    )

    logger.info(notes)

    return ModificationResult(
        step=step, acs_before={}, acs_after={},
        conversion=cr.conversion, delta_G_DN=0.0, notes=notes,
    )


# ─── Workflow 6: Spacer Arm (v5.8 WN-2) ────────────────────────────

def _solve_spacer_arm_step(
    step: ModificationStep,
    acs_state: dict[ACSSiteType, ACSProfile],
    target_profile: ACSProfile,
    reagent_profile: ReagentProfile,
    surface_model: AccessibleSurfaceModel,
    V_bead: float,
) -> ModificationResult:
    """Couple spacer arm or heterobifunctional crosslinker to activated sites.

    Consumes activated sites on target profile and creates a NEW intermediate
    ACS profile for the product site type (e.g., EPOXIDE -> AMINE_DISTAL,
    or AMINE_DISTAL -> MALEIMIDE).

    Key differences from ACTIVATION:
    - Creates a new ACSProfile entry in acs_state (not just updating activated_sites)
    - Applies distal_group_yield: not all consumed sites produce useful distal groups
    - Bridged sites (1 - yield) go to crosslinked_sites on source profile
    - Product profile gets profile_role based on reagent type
    """
    remaining_act = target_profile.remaining_activated
    if remaining_act <= 0:
        step.conversion = 0.0
        return ModificationResult(
            step=step, acs_before={}, acs_after={},
            conversion=0.0, delta_G_DN=0.0,
            notes="Spacer arm: no remaining activated sites.",
        )

    # Arrhenius prefactor
    k0 = _arrhenius_prefactor(
        reagent_profile.k_forward, reagent_profile.E_a,
        reagent_profile.temperature_default,
    ) if reagent_profile.E_a > 0 else 0.0

    # ODE: same second-order consumption as crosslinking/activation
    acs_concentration = remaining_act / V_bead
    conversion, reagent_remaining = solve_second_order_consumption(
        acs_concentration=acs_concentration,
        reagent_concentration=step.reagent_concentration,
        k_forward=reagent_profile.k_forward,
        stoichiometry=step.stoichiometry,
        time=step.time,
        temperature=step.temperature,
        E_a=reagent_profile.E_a,
        k0=k0,
        hydrolysis_rate=reagent_profile.hydrolysis_rate,
    )

    sites_consumed = conversion * remaining_act
    step.conversion = conversion

    # Apply distal group yield (audit F10): not all sites produce useful product
    distal_yield = getattr(reagent_profile, 'distal_group_yield', 1.0)
    sites_created = sites_consumed * distal_yield
    sites_bridged = sites_consumed * (1.0 - distal_yield)

    # Update source profile: ALL consumed sites go to ligand_coupled (attached to spacer)
    # The distal yield reduction affects product creation, not source accounting
    target_profile.ligand_coupled_sites += sites_consumed

    # Determine product ACS type and profile role
    product_type = step.product_acs or getattr(reagent_profile, 'product_acs', None)
    if product_type is None:
        logger.warning("Spacer arm step has no product_acs; no intermediate created.")
        return ModificationResult(
            step=step, acs_before={}, acs_after={},
            conversion=conversion, delta_G_DN=0.0,
            notes=f"Spacer arm: {sites_consumed:.4e} consumed, no product_acs specified.",
        )

    # Determine profile role for the new intermediate
    _rt = getattr(reagent_profile, 'reaction_type', 'spacer_arm')
    if _rt == "heterobifunctional":
        _prod_role = "heterobifunctional_intermediate"
    else:
        _prod_role = "spacer_intermediate"

    # Create or update product ACS profile
    if product_type in acs_state:
        # Accumulate onto existing intermediate (allowed for compatible chemistry)
        acs_state[product_type].accessible_sites += sites_created
        acs_state[product_type].total_sites += sites_created
        acs_state[product_type].activated_sites += sites_created
    else:
        acs_state[product_type] = ACSProfile(
            site_type=product_type,
            total_sites=sites_created,
            accessible_sites=sites_created,
            activated_sites=sites_created,
            crosslinked_sites=0.0,
            activated_consumed_sites=0.0,
            hydrolyzed_sites=0.0,
            ligand_coupled_sites=0.0,
            ligand_functional_sites=0.0,
            blocked_sites=0.0,
            total_density=0.0,
            accessible_density=0.0,
            accessibility_model=target_profile.accessibility_model,
            uncertainty_fraction=target_profile.uncertainty_fraction,
        )

    notes = (
        f"Spacer arm '{step.reagent_key}': conv={conversion:.4f}, "
        f"consumed={sites_consumed:.4e}, created={sites_created:.4e} "
        f"(yield={distal_yield:.0%}), bridged={sites_bridged:.4e}, "
        f"product={product_type.value}, role={_prod_role}"
    )
    logger.info(notes)

    return ModificationResult(
        step=step, acs_before={}, acs_after={},
        conversion=conversion, delta_G_DN=0.0, notes=notes,
    )


# ─── Workflow 7: Metal Charging (v5.9.1) ──────────────���─────────────

def _solve_metal_charging_step(
    step: ModificationStep,
    acs_state: dict[ACSSiteType, ACSProfile],
    target_profile: ACSProfile,
    reagent_profile: ReagentProfile,
) -> ModificationResult:
    """Load or strip metal ions on IMAC chelator sites.

    Langmuir equilibrium: fraction = K * C / (1 + K * C)
    Does NOT consume ACS sites. Updates metal_loaded_fraction on target profile.
    For stripping (EDTA): sets metal_loaded_fraction to ~0.
    """
    K_assoc = getattr(reagent_profile, 'metal_association_constant', 1e10)
    C_metal = step.reagent_concentration  # [mol/m3]

    if reagent_profile.reaction_type == "metal_stripping":
        # EDTA strips metal: fraction -> ~0
        # ASSUMPTION: 50 mM EDTA at pH 7 strips >99% of Ni/Co/Cu
        stripping_eff = 0.99
        old_frac = target_profile.metal_loaded_fraction
        target_profile.metal_loaded_fraction = old_frac * (1.0 - stripping_eff)
        step.conversion = stripping_eff
        notes = f"Metal stripping: {old_frac:.3f} -> {target_profile.metal_loaded_fraction:.4f}"
    else:
        # Langmuir equilibrium loading
        if K_assoc > 0 and C_metal > 0:
            frac = K_assoc * C_metal / (1.0 + K_assoc * C_metal)
        else:
            frac = 0.0
        target_profile.metal_loaded_fraction = min(frac, 1.0)
        step.conversion = frac
        notes = (
            f"Metal charging '{step.reagent_key}': "
            f"K={K_assoc:.1e}, C={C_metal:.1f} mol/m3, "
            f"fraction={frac:.4f}"
        )

    logger.info(notes)
    return ModificationResult(
        step=step, acs_before={}, acs_after={},
        conversion=step.conversion, delta_G_DN=0.0, notes=notes,
    )


# ─── Workflow 8: Protein Pretreatment (v5.9.2) ──────────────────────

def _solve_protein_pretreatment_step(
    step: ModificationStep,
    reagent_profile: ReagentProfile,
) -> ModificationResult:
    """Reduce protein disulfides to expose free thiols for maleimide coupling.

    Does NOT modify ACS state. Returns a ModificationResult with pretreatment
    info in notes. The orchestrator stores the result for downstream coupling.

    Model: free_thiol_fraction = efficiency * (1 - exp(-k * t))
    """
    k_red = reagent_profile.k_forward  # reduction rate [1/s]
    t = step.time
    efficiency = getattr(reagent_profile, 'activity_retention', 0.95)
    # ASSUMPTION: activity_retention on pretreatment profile = reduction efficiency

    if k_red > 0 and t > 0:
        frac = efficiency * (1.0 - math.exp(-k_red * t))
    else:
        frac = 0.0

    step.conversion = frac

    # Check buffer incompatibilities
    _incompat = getattr(reagent_profile, 'buffer_incompatibilities', '')
    _warns = []
    if _incompat and "maleimide" in _incompat:
        _warns.append(f"WARNING: {step.reagent_key} is incompatible with maleimide. "
                      "Remove excess reductant before maleimide-thiol coupling.")

    notes = (
        f"Protein pretreatment '{step.reagent_key}': "
        f"free_thiol_fraction={frac:.3f}, efficiency={efficiency:.2f}"
    )
    if _warns:
        notes += " | " + "; ".join(_warns)

    logger.info(notes)
    return ModificationResult(
        step=step, acs_before={}, acs_after={},
        conversion=frac, delta_G_DN=0.0, notes=notes,
    )


# ─── Workflow 9: Washing (v5.9.4) ───────────────────────────────────

def _solve_washing_step(
    step: ModificationStep,
    acs_state: dict[ACSSiteType, ACSProfile],
    surface_model: AccessibleSurfaceModel,
    V_bead: float,
) -> ModificationResult:
    """Advisory washing model — simple diffusion-out from pore volume.

    C_residual = C_initial * exp(-D_eff * t / (R^2 / pi^2))
    Does NOT consume ACS sites. Reports residual concentration for screening.
    """
    R = surface_model.bead_radius
    D_eff = 1e-10  # ASSUMPTION: default diffusivity in gel [m2/s]
    t_wash = step.time
    C_initial = step.reagent_concentration  # [mol/m3] of wash buffer (proxy for volume)

    if R > 0 and t_wash > 0:
        # Characteristic diffusion time: R^2 / (pi^2 * D_eff)
        tau = R ** 2 / (math.pi ** 2 * D_eff)
        fraction_removed = 1.0 - math.exp(-t_wash / tau)
    else:
        fraction_removed = 0.0

    step.conversion = fraction_removed

    notes = (
        f"Washing: fraction_removed={fraction_removed:.4f}, "
        f"tau_diff={tau:.0f}s, t_wash={t_wash:.0f}s. "
        f"Advisory screening_only."
    )
    logger.info(notes)

    return ModificationResult(
        step=step, acs_before={}, acs_after={},
        conversion=fraction_removed, delta_G_DN=0.0, notes=notes,
    )


# ─── Helpers ───────────────────────────────────────────────────────────

def _snapshot_acs(
    acs_state: dict[ACSSiteType, ACSProfile],
) -> dict[ACSSiteType, ACSProfile]:
    """Deep copy of ACS state for before/after comparison."""
    return {k: copy.deepcopy(v) for k, v in acs_state.items()}
