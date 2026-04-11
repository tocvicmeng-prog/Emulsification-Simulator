"""Module 2: Functionalization — surface chemistry and ligand coupling.

Phase A: ACS data model with terminal-state conservation (v2).
Phase B: 5 validated chemistry workflows:
  - Workflow 1: Amine secondary crosslinking (genipin, glutaraldehyde)
  - Workflow 2: Hydroxyl activation (ECH, DVS)
  - Workflow 3: Ligand coupling (DEAE, IDA, Phenyl, Sulfopropyl)
  - Workflow 4: Protein coupling (Protein A, Protein G) — ranking_only
  - Workflow 5: Quenching (ethanolamine, mercaptoethanol, NaBH4, acetic anhydride)

14 reagent profiles. Backend workflow validation.
FunctionalMediaContract bridges M2→M3 with ligand density → q_max mapping.
"""

from .acs import ACSSiteType, ACSProfile, initialize_acs_from_m1
from .modification_steps import (
    ModificationStep,
    ModificationStepType,
    ModificationResult,
    solve_modification_step,
)
from .orchestrator import (
    FunctionalMicrosphere,
    FunctionalMediaContract,
    ModificationOrchestrator,
    build_functional_media_contract,
)
from .reactions import CouplingResult
from .reagent_profiles import ReagentProfile, REAGENT_PROFILES
from .surface_area import AccessibleSurfaceModel, SurfaceAreaTier

__all__ = [
    # Data model
    "ACSSiteType",
    "ACSProfile",
    "AccessibleSurfaceModel",
    "SurfaceAreaTier",
    "initialize_acs_from_m1",
    # Workflows
    "ModificationStep",
    "ModificationStepType",
    "ModificationResult",
    "solve_modification_step",
    "CouplingResult",
    # Orchestration
    "FunctionalMicrosphere",
    "FunctionalMediaContract",
    "ModificationOrchestrator",
    "build_functional_media_contract",
    # Reagents
    "ReagentProfile",
    "REAGENT_PROFILES",
]
