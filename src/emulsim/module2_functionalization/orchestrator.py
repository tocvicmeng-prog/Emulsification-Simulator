"""Module 2 functionalization orchestrator.

Phase B: Minimal Module 2 — 2 Workflows.
Architecture: module2_module3_final_implementation_plan.md, Phase B.

Provides FunctionalMicrosphere (complete description of a functionalized
microsphere) and ModificationOrchestrator (sequential execution of
modification steps with ACS tracking).
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import TYPE_CHECKING

from .acs import ACSSiteType, ACSProfile, initialize_acs_from_m1
from .modification_steps import (
    ModificationResult,
    ModificationStep,
    solve_modification_step,
)
from .reagent_profiles import REAGENT_PROFILES, ReagentProfile
from .surface_area import AccessibleSurfaceModel

if TYPE_CHECKING:
    from emulsim.datatypes import M1ExportContract

logger = logging.getLogger(__name__)


# ─── FunctionalMicrosphere ────────────────────────────────────────────

@dataclass
class FunctionalMicrosphere:
    """Complete description of a functionalized microsphere.

    Holds the M1 contract, surface model, current ACS state, and the
    full modification history with provenance.

    Attributes:
        m1_contract: Source Module 1 export contract.
        surface_model: Computed accessible surface area model.
        acs_profiles: Current ACS inventory (mutated by modification steps).
        modification_history: Ordered list of completed modification results.
        G_DN_updated: Updated double-network shear modulus [Pa] after
            secondary crosslinking.
        E_star_updated: Updated effective Young's modulus [Pa].
    """
    m1_contract: M1ExportContract
    surface_model: AccessibleSurfaceModel
    acs_profiles: dict[ACSSiteType, ACSProfile]
    modification_history: list[ModificationResult] = field(default_factory=list)
    G_DN_updated: float = 0.0
    E_star_updated: float = 0.0

    def validate(self) -> list[str]:
        """Validate ACS conservation across all profiles.

        Returns:
            List of violation messages (empty = all valid).
        """
        errors: list[str] = []
        for profile in self.acs_profiles.values():
            errors.extend(profile.validate())
        return errors


# ─── ModificationOrchestrator ─────────────────────────────────────────

class ModificationOrchestrator:
    """Execute sequential modification steps on a microsphere.

    Usage:
        orchestrator = ModificationOrchestrator()
        result = orchestrator.run(contract, steps)
        assert result.validate() == []
    """

    def run(
        self,
        contract: M1ExportContract,
        steps: list[ModificationStep],
    ) -> FunctionalMicrosphere:
        """Execute all modification steps sequentially, tracking ACS.

        Algorithm:
            1. Build AccessibleSurfaceModel from M1 contract.
            2. Initialize ACS inventory from contract + surface model.
            3. For each step: look up reagent, solve, record result.
            4. Accumulate G_DN updates from crosslinking steps.
            5. Return FunctionalMicrosphere with full provenance.

        Args:
            contract: M1ExportContract from Module 1.
            steps: Ordered list of modification steps to execute.

        Returns:
            FunctionalMicrosphere with updated ACS and modification history.

        Raises:
            KeyError: If a step references an unknown reagent_key.
            ValueError: If a step targets a missing ACS type.
        """
        # --- Build surface model ---
        surface_model = AccessibleSurfaceModel.from_m1_export(contract)

        # --- Initialize ACS from M1 ---
        acs_profiles = initialize_acs_from_m1(contract, surface_model)

        logger.info(
            "ModificationOrchestrator: %d steps, %d ACS types initialized.",
            len(steps), len(acs_profiles),
        )

        # --- Execute steps sequentially ---
        history: list[ModificationResult] = []
        total_delta_G = 0.0

        for i, step in enumerate(steps):
            # Look up reagent profile
            if step.reagent_key not in REAGENT_PROFILES:
                raise KeyError(
                    f"Step {i}: unknown reagent_key '{step.reagent_key}'. "
                    f"Available: {list(REAGENT_PROFILES.keys())}"
                )
            reagent_profile = REAGENT_PROFILES[step.reagent_key]

            logger.info(
                "Step %d/%d: %s with %s on %s",
                i + 1, len(steps), step.step_type.value,
                step.reagent_key, step.target_acs.value,
            )

            result = solve_modification_step(
                step=step,
                acs_state=acs_profiles,
                surface_model=surface_model,
                reagent_profile=reagent_profile,
            )
            history.append(result)
            total_delta_G += result.delta_G_DN

            logger.info(
                "Step %d: conversion=%.4f, delta_G=%.1f Pa",
                i + 1, result.conversion, result.delta_G_DN,
            )

        # --- Compute updated mechanical properties ---
        G_DN_base = contract.G_DN
        G_DN_updated = G_DN_base + total_delta_G
        # E* ~ 3*G for incompressible rubber (Poisson's ratio ~ 0.5)
        E_star_updated = 3.0 * G_DN_updated

        logger.info(
            "Orchestrator complete: G_DN %.1f -> %.1f Pa (delta=%.1f Pa)",
            G_DN_base, G_DN_updated, total_delta_G,
        )

        return FunctionalMicrosphere(
            m1_contract=contract,
            surface_model=surface_model,
            acs_profiles=acs_profiles,
            modification_history=history,
            G_DN_updated=G_DN_updated,
            E_star_updated=E_star_updated,
        )
