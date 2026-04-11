"""Available Coupling Site (ACS) inventory for Module 2 functionalization.

Phase A+: Revised ACS state model with explicit terminal states.
Architecture: docs/18_m2_expansion_final_plan.md, Phase 1.

Provides ACSProfile — a state-machine dataclass tracking the full hierarchy
of reactive sites through sequential modification steps.

Conservation invariant (v2 — terminal-state model):
    accessible_sites = remaining_sites
                     + crosslinked_sites     (consumed by crosslinking)
                     + hydrolyzed_sites      (lost to aqueous hydrolysis)
                     + ligand_coupled_sites  (coupled to functional ligand)
                     + blocked_sites         (capped by quenching)

Every accessible site ends in exactly ONE terminal state.
The validate() method checks this and returns a list of violations.
"""

from __future__ import annotations

import logging
import math
from dataclasses import dataclass
from enum import Enum
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from emulsim.datatypes import M1ExportContract
    from emulsim.module2_functionalization.surface_area import AccessibleSurfaceModel

logger = logging.getLogger(__name__)


# ─── Constants ──────────────────────────────────────────────────────────

# Tolerance for floating-point conservation checks.
# 0.1% slack prevents false violations from rounding.
_CONSERVATION_TOL = 1.001


# ─── Enums ──────────────────────────────────────────────────────────────

class ACSSiteType(Enum):
    """Reactive coupling site types on hydrogel microspheres.

    AMINE_PRIMARY:   -NH2 on chitosan (from deacetylated glucosamine units).
    HYDROXYL:        -OH on agarose backbone.
    CARBOXYL:        -COOH (e.g., oxidised cellulose, carboxymethyl chitosan).
    EPOXIDE:         Epoxy groups introduced by ECH activation.
    ALDEHYDE:        -CHO from glutaraldehyde or periodate oxidation.
    VINYL_SULFONE:   VS groups from DVS activation.
    """
    AMINE_PRIMARY = "amine_primary"
    HYDROXYL = "hydroxyl"
    CARBOXYL = "carboxyl"
    EPOXIDE = "epoxide"
    ALDEHYDE = "aldehyde"
    VINYL_SULFONE = "vinyl_sulfone"


# ─── ACSProfile Dataclass ──────────────────────────────────────────────

@dataclass
class ACSProfile:
    """Available Coupling Site inventory for a single site type.

    Terminal-state conservation model (v2):
        accessible_sites = remaining_sites
                         + crosslinked_sites
                         + hydrolyzed_sites
                         + ligand_coupled_sites
                         + blocked_sites

    Every accessible site ends in exactly ONE terminal state.

    Hierarchy constraints:
        total_sites >= accessible_sites >= activated_sites
        ligand_functional_sites <= ligand_coupled_sites
        terminal_sum <= accessible_sites

    For activated product profiles (EPOXIDE, VINYL_SULFONE, ALDEHYDE):
        ligand_coupled + hydrolyzed + blocked <= activated_sites

    Attributes:
        site_type: Type of reactive group.
        total_sites: Total sites per particle (bulk, including buried) [mol/particle].
        accessible_sites: Sites accessible to reagent (surface + pore) [mol/particle].
        activated_sites: Sites activated by chemistry (epoxy, aldehyde, etc.) [mol/particle].
        crosslinked_sites: Sites consumed by crosslinking (contributes to G_DN) [mol/particle].
        hydrolyzed_sites: Sites lost to aqueous hydrolysis (waste) [mol/particle].
        ligand_coupled_sites: Ligand molecules coupled to activated sites [mol/particle].
        ligand_functional_sites: Coupled ligand retaining biological activity [mol/particle].
        blocked_sites: Sites intentionally blocked by quenching [mol/particle].
        total_density: Total site density on reagent-accessible area [mol/m^2].
        accessible_density: Accessible site density [mol/m^2].
        accessibility_model: Name of surface area tier used.
        uncertainty_fraction: Fractional uncertainty on site counts [-].
    """
    site_type: ACSSiteType
    total_sites: float = 0.0           # [mol/particle]
    accessible_sites: float = 0.0      # [mol/particle]
    activated_sites: float = 0.0       # [mol/particle]

    # ── Terminal states (mutually exclusive destinations) ──
    crosslinked_sites: float = 0.0     # [mol/particle] consumed by crosslinking (G_DN)
    activated_consumed_sites: float = 0.0  # [mol/particle] consumed by activation chemistry (not crosslinking)
    hydrolyzed_sites: float = 0.0      # [mol/particle] lost to hydrolysis
    ligand_coupled_sites: float = 0.0  # [mol/particle] coupled to ligand
    ligand_functional_sites: float = 0.0  # [mol/particle] subset retaining activity
    blocked_sites: float = 0.0         # [mol/particle] capped by quenching

    # Surface-normalised densities [mol/m^2]
    total_density: float = 0.0
    accessible_density: float = 0.0

    accessibility_model: str = "empirical_pore"
    uncertainty_fraction: float = 0.1

    @property
    def _terminal_sum(self) -> float:
        """Sum of all terminal-state site counts [mol/particle]."""
        return (self.crosslinked_sites + self.activated_consumed_sites
                + self.hydrolyzed_sites + self.ligand_coupled_sites
                + self.blocked_sites)

    @property
    def remaining_sites(self) -> float:
        """Sites still available for further chemistry [mol/particle].

        remaining = accessible - terminal_sum (clamped to >= 0).
        Terminal states: crosslinked + hydrolyzed + ligand_coupled + blocked.
        """
        return max(self.accessible_sites - self._terminal_sum, 0.0)

    @property
    def remaining_activated(self) -> float:
        """Activated sites not yet consumed by coupling, hydrolysis, or quenching.

        For activated product profiles (EPOXIDE, VINYL_SULFONE, ALDEHYDE):
        remaining_activated = activated - (ligand_coupled + hydrolyzed + blocked)
        """
        used = self.ligand_coupled_sites + self.hydrolyzed_sites + self.blocked_sites
        return max(self.activated_sites - used, 0.0)

    # Backward compatibility: consumed_sites as computed property
    @property
    def consumed_sites(self) -> float:
        """Backward-compatible view: total sites in any terminal state.

        Equivalent to the old consumed_sites + blocked_sites combined.
        Use terminal-state fields directly for new code.
        """
        return self._terminal_sum

    def validate(self) -> list[str]:
        """Conservation hierarchy checks (v2 terminal-state model).

        Returns:
            List of violation messages (empty = valid).
        """
        errors: list[str] = []
        name = self.site_type.value

        # Hierarchy: accessible <= total
        if self.accessible_sites > self.total_sites * _CONSERVATION_TOL:
            errors.append(f"{name}: accessible ({self.accessible_sites:.4e}) > total ({self.total_sites:.4e})")

        # Hierarchy: activated <= accessible
        if self.activated_sites > self.accessible_sites * _CONSERVATION_TOL:
            errors.append(f"{name}: activated ({self.activated_sites:.4e}) > accessible ({self.accessible_sites:.4e})")

        # Hierarchy: ligand_functional <= ligand_coupled
        if self.ligand_functional_sites > self.ligand_coupled_sites * _CONSERVATION_TOL:
            errors.append(f"{name}: ligand_functional ({self.ligand_functional_sites:.4e}) > ligand_coupled ({self.ligand_coupled_sites:.4e})")

        # Terminal-sum conservation: crosslinked + hydrolyzed + coupled + blocked <= accessible
        ts = self._terminal_sum
        if ts > self.accessible_sites * _CONSERVATION_TOL:
            errors.append(f"{name}: terminal_sum ({ts:.4e}) > accessible ({self.accessible_sites:.4e})")

        # For activated product profiles: ligand_coupled + hydrolyzed + blocked <= activated
        activated_used = self.ligand_coupled_sites + self.hydrolyzed_sites + self.blocked_sites
        if self.activated_sites > 0 and activated_used > self.activated_sites * _CONSERVATION_TOL:
            errors.append(
                f"{name}: ligand_coupled+hydrolyzed+blocked ({activated_used:.4e}) "
                f"> activated ({self.activated_sites:.4e})"
            )

        # Non-negativity
        for attr in ("total_sites", "accessible_sites", "activated_sites",
                      "crosslinked_sites", "activated_consumed_sites",
                      "hydrolyzed_sites",
                      "ligand_coupled_sites", "ligand_functional_sites",
                      "blocked_sites"):
            val = getattr(self, attr)
            if val < 0:
                errors.append(f"{name}: {attr} is negative ({val:.4e})")

        return errors


# ─── Initialization from M1 ───────────────────────────────────────────

def initialize_acs_from_m1(
    contract: M1ExportContract,
    surface_model: AccessibleSurfaceModel,
) -> dict[ACSSiteType, ACSProfile]:
    """Initialize ACS inventory from Module 1 export contract.

    Creates ACSProfile for each reactive group type present (NH2 and OH),
    using the surface area model to convert bulk concentrations to
    per-particle site counts and surface densities.

    Algorithm:
        1. Compute bead volume: V_bead = (4/3)*pi*R^3.
        2. Total sites per particle: n_total = C_bulk * V_bead.
        3. Accessible fraction = reagent_accessible_area / (external + internal).
           For EXTERNAL_ONLY tier, accessible = external_area / (external + internal)
           times total, but since internal=0, accessible = n_total * (ext/ext) = n_total
           for that degenerate case.
        4. Surface density = accessible_sites / reagent_accessible_area.

    Args:
        contract: M1ExportContract with nh2_bulk_concentration and oh_bulk_concentration.
        surface_model: Pre-computed AccessibleSurfaceModel (compute() must have been called).

    Returns:
        Dictionary mapping ACSSiteType to ACSProfile for AMINE_PRIMARY and HYDROXYL.

    Raises:
        ValueError: If surface model has not been computed (zero areas).
    """
    # --- Validate surface model was computed ---
    if surface_model.external_area <= 0:
        raise ValueError(
            "AccessibleSurfaceModel.compute() must be called before "
            "initialize_acs_from_m1(). external_area is 0."
        )

    R = surface_model.bead_radius
    V_bead = (4.0 / 3.0) * math.pi * R ** 3  # [m^3/particle]

    profiles: dict[ACSSiteType, ACSProfile] = {}

    # --- Build profiles for each reactive group ---
    _group_specs = [
        (ACSSiteType.AMINE_PRIMARY, contract.nh2_bulk_concentration),
        (ACSSiteType.HYDROXYL, contract.oh_bulk_concentration),
    ]

    for site_type, c_bulk in _group_specs:
        if c_bulk <= 0:
            logger.debug(
                "Skipping %s: bulk concentration is %.3e mol/m^3",
                site_type.value, c_bulk,
            )
            continue

        # Total sites per particle: n = C_bulk [mol/m^3] * V_bead [m^3]
        n_total = c_bulk * V_bead  # [mol/particle]

        # Accessible sites: scaled by ratio of reagent-accessible area to
        # total geometric area (external + internal).
        total_geometric = surface_model.external_area + surface_model.internal_geometric_area
        if total_geometric > 0:
            accessibility_ratio = surface_model.reagent_accessible_area / total_geometric
        else:
            accessibility_ratio = 1.0  # degenerate: no internal area
        n_accessible = n_total * accessibility_ratio  # [mol/particle]

        # Surface densities [mol/m^2]
        # Guard: division by zero if area is somehow zero (shouldn't happen after validation)
        total_density = (
            n_total / surface_model.reagent_accessible_area
            if surface_model.reagent_accessible_area > 0
            else 0.0
        )
        accessible_density = (
            n_accessible / surface_model.reagent_accessible_area
            if surface_model.reagent_accessible_area > 0
            else 0.0
        )

        profile = ACSProfile(
            site_type=site_type,
            total_sites=n_total,
            accessible_sites=n_accessible,
            activated_sites=0.0,
            crosslinked_sites=0.0,
            activated_consumed_sites=0.0,
            hydrolyzed_sites=0.0,
            ligand_coupled_sites=0.0,
            ligand_functional_sites=0.0,
            blocked_sites=0.0,
            total_density=total_density,
            accessible_density=accessible_density,
            accessibility_model=surface_model.tier.value,
            uncertainty_fraction=0.1 if surface_model.trust_level != "UNRELIABLE" else 0.3,
        )

        # Validate before storing
        violations = profile.validate()
        if violations:
            logger.warning(
                "ACS initialization violations for %s: %s",
                site_type.value, "; ".join(violations),
            )

        profiles[site_type] = profile

    logger.info(
        "Initialized ACS from M1: %d site types, V_bead=%.3e m^3",
        len(profiles), V_bead,
    )
    return profiles
