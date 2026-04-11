"""Module 2: Functionalization — surface chemistry and ligand coupling.

Phase A provides the data model and ACS conservation arithmetic:
  - AccessibleSurfaceModel: 3-tier surface area estimation.
  - ACSProfile: Available Coupling Site inventory with conservation checks.
  - initialize_acs_from_m1: Factory to seed ACS from Module 1 export contract.

Phase B adds 2 validated chemistry workflows:
  - Workflow 1: Amine secondary crosslinking (genipin, glutaraldehyde).
  - Workflow 2: Hydroxyl activation (ECH, DVS).
  - ModificationOrchestrator: Sequential step execution with ACS tracking.
  - FunctionalMicrosphere: Complete functionalized particle description.
  - ReagentProfile / REAGENT_PROFILES: 4 reagent kinetic profiles.
"""

from .acs import ACSSiteType, ACSProfile, initialize_acs_from_m1
from .modification_steps import (
    ModificationStep,
    ModificationStepType,
    ModificationResult,
    solve_modification_step,
)
from .orchestrator import FunctionalMicrosphere, ModificationOrchestrator
from .reagent_profiles import ReagentProfile, REAGENT_PROFILES
from .surface_area import AccessibleSurfaceModel, SurfaceAreaTier

__all__ = [
    # Phase A
    "ACSSiteType",
    "ACSProfile",
    "AccessibleSurfaceModel",
    "SurfaceAreaTier",
    "initialize_acs_from_m1",
    # Phase B
    "ModificationStep",
    "ModificationStepType",
    "ModificationResult",
    "solve_modification_step",
    "FunctionalMicrosphere",
    "ModificationOrchestrator",
    "ReagentProfile",
    "REAGENT_PROFILES",
]
