"""Modification step definitions and workflow solvers for Module 2.

Phase B: Minimal Module 2 — 2 Workflows.
Architecture: module2_module3_final_implementation_plan.md, Phase B.

Two validated workflows:

Workflow 1 — Amine secondary crosslinking:
    Target: AMINE_PRIMARY sites.
    Uses solve_second_order_consumption with reagent's k and E_a.
    Updates: consumed_sites increases, remaining decreases.
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
from dataclasses import dataclass, field
from enum import Enum
from typing import Optional

from .acs import ACSSiteType, ACSProfile
from .reactions import solve_second_order_consumption
from .reagent_profiles import ReagentProfile
from .surface_area import AccessibleSurfaceModel

logger = logging.getLogger(__name__)

# ─── Constants ──────────────────────────────────────────────────────────

K_BOLTZMANN = 1.38e-23  # [J/K]
N_AVOGADRO = 6.022e23


# ─── Enums ──────────────────────────────────────────────────────────────

class ModificationStepType(Enum):
    """Type of surface modification step."""
    SECONDARY_CROSSLINKING = "secondary_crosslinking"
    ACTIVATION = "activation"
    LIGAND_COUPLING = "ligand_coupling"
    PROTEIN_COUPLING = "protein_coupling"
    QUENCHING = "quenching"


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
    # ACS remaining is in [mol/particle].
    # To get a bulk concentration [mol/m^3], we need volume per particle.
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
    else:
        raise KeyError(
            f"Step type {step.step_type.value} not implemented in Phase B. "
            f"Supported: SECONDARY_CROSSLINKING, ACTIVATION."
        )

    # --- Snapshot after ---
    result.acs_before = acs_before
    result.acs_after = _snapshot_acs(acs_state)

    # --- Validate conservation after step ---
    for profile in acs_state.values():
        violations = profile.validate()
        if violations:
            result.notes += f" Conservation violations: {'; '.join(violations)}"
            logger.warning("ACS conservation violation: %s", violations)

    return result


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

    Consumes AMINE_PRIMARY sites, increases consumed_sites count.
    Estimates delta_G_DN from additive rubber elasticity:
        delta_G = delta_nu_e * kB * T
    where delta_nu_e is the additional effective crosslink density.
    """
    conversion, reagent_remaining = solve_second_order_consumption(
        acs_concentration=acs_concentration,
        reagent_concentration=step.reagent_concentration,
        k_forward=reagent_profile.k_forward,
        stoichiometry=step.stoichiometry,
        time=step.time,
        temperature=step.temperature,
        E_a=reagent_profile.E_a,
        k0=0.0,  # Use k_forward directly (already at reference T)
        hydrolysis_rate=reagent_profile.hydrolysis_rate,
    )

    # Update ACS state
    sites_consumed = conversion * target_profile.remaining_sites
    target_profile.consumed_sites += sites_consumed

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
    conversion, reagent_remaining = solve_second_order_consumption(
        acs_concentration=acs_concentration,
        reagent_concentration=step.reagent_concentration,
        k_forward=reagent_profile.k_forward,
        stoichiometry=step.stoichiometry,
        time=step.time,
        temperature=step.temperature,
        E_a=reagent_profile.E_a,
        k0=0.0,  # Use k_forward directly
        hydrolysis_rate=reagent_profile.hydrolysis_rate,
    )

    # Sites consumed from target (hydroxyl)
    sites_consumed = conversion * target_profile.remaining_sites
    target_profile.consumed_sites += sites_consumed

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
                consumed_sites=0.0,
                blocked_sites=0.0,
                ligand_coupled_sites=0.0,
                ligand_functional_sites=0.0,
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


# ─── Helpers ───────────────────────────────────────────────────────────

def _snapshot_acs(
    acs_state: dict[ACSSiteType, ACSProfile],
) -> dict[ACSSiteType, ACSProfile]:
    """Deep copy of ACS state for before/after comparison."""
    return {k: copy.deepcopy(v) for k, v in acs_state.items()}
