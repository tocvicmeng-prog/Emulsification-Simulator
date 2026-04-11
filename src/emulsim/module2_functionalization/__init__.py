"""Module 2: Functionalization — surface chemistry and ligand coupling.

Phase A provides the data model and ACS conservation arithmetic:
  - AccessibleSurfaceModel: 3-tier surface area estimation.
  - ACSProfile: Available Coupling Site inventory with conservation checks.
  - initialize_acs_from_m1: Factory to seed ACS from Module 1 export contract.
"""

from .acs import ACSSiteType, ACSProfile, initialize_acs_from_m1
from .surface_area import AccessibleSurfaceModel, SurfaceAreaTier

__all__ = [
    "ACSSiteType",
    "ACSProfile",
    "AccessibleSurfaceModel",
    "SurfaceAreaTier",
    "initialize_acs_from_m1",
]
