"""3-tier internal surface area model bridging L2 pore morphology to M2 ACS.

Phase A: Data Model + ACS Surface-Area Core.
Architecture: module2_module3_final_implementation_plan.md, Phase A.

Provides AccessibleSurfaceModel which computes external, internal, and
accessibility-corrected surface areas for porous hydrogel microspheres.
These areas are the denominator for all ACS surface-density calculations.

Key formulas (cylindrical pore approximation):
  - External area:     A_ext = 4 * pi * R^2
  - Bead volume:       V_bead = (4/3) * pi * R^3
  - Specific surface:  S_v = 4 * epsilon / d_pore   [m^2/m^3]
  - Internal area:     A_int = S_v * V_bead
  - Accessible fraction for solute of hydrodynamic radius r_h:
        f = max((1 - (2*r_h / d_pore)^2), 0) / tortuosity
  - Reagent-accessible area:  A_ext + f_reagent * A_int
  - Ligand-accessible area:   A_ext + f_ligand  * A_int
"""

from __future__ import annotations

import logging
import math
from dataclasses import dataclass, field
from enum import Enum
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from emulsim.datatypes import M1ExportContract

logger = logging.getLogger(__name__)


# ─── Constants ──────────────────────────────────────────────────────────

_DEFAULT_REAGENT_RADIUS = 0.5e-9   # [m] small-molecule reagent (e.g. glutaraldehyde ~0.4 nm)
_DEFAULT_LIGAND_RADIUS = 3.0e-9    # [m] globular protein ligand (e.g. Protein A ~3 nm r_h)
_DEFAULT_TORTUOSITY = 1.5          # [-] typical for macroporous gels


# ─── Enums ──────────────────────────────────────────────────────────────

class SurfaceAreaTier(Enum):
    """Surface area estimation tier.

    EXTERNAL_ONLY:      Uses only the outer bead surface (4*pi*R^2).
                        Suitable when pore data is absent.
    EMPIRICAL_PORE:     Cylindrical-pore approximation using mean pore diameter
                        and porosity.  Default for most M1 exports.
    MORPHOLOGY_BASED:   Uses L2 morphology descriptors (chord lengths, connectivity)
                        for refined estimates.  Requires mechanistic L2 output.
    """
    EXTERNAL_ONLY = "external_only"
    EMPIRICAL_PORE = "empirical_pore"
    MORPHOLOGY_BASED = "morphology_based"


# ─── Trust level mapping ───────────────────────────────────────────────

_TIER_TRUST = {
    SurfaceAreaTier.EXTERNAL_ONLY: "UNRELIABLE",
    SurfaceAreaTier.EMPIRICAL_PORE: "CAUTION",
    SurfaceAreaTier.MORPHOLOGY_BASED: "TRUSTWORTHY",
}


# ─── Dataclass ──────────────────────────────────────────────────────────

@dataclass
class AccessibleSurfaceModel:
    """3-tier internal surface area model bridging L2 pore morphology to M2 ACS.

    Geometric inputs describe a single microsphere bead.
    All computed areas are per-particle [m^2/particle].

    Attributes:
        tier: Surface area estimation tier.
        bead_radius: Bead radius [m].
        porosity: Volume fraction of pore space [-].
        pore_diameter_mean: Mean pore diameter [m].
        xi_mesh: Crosslinked mesh size [m] (from L3).
        external_area: Outer bead surface area [m^2/particle].
        internal_geometric_area: Total internal pore area (geometric) [m^2/particle].
        reagent_accessible_area: Area accessible to small-molecule reagent [m^2/particle].
        ligand_accessible_area: Area accessible to protein/ligand [m^2/particle].
        trust_level: Inherited from tier; "UNRELIABLE", "CAUTION", or "TRUSTWORTHY".
    """
    tier: SurfaceAreaTier = SurfaceAreaTier.EMPIRICAL_PORE
    bead_radius: float = 50e-6          # [m]
    porosity: float = 0.7               # [-]
    pore_diameter_mean: float = 100e-9  # [m]
    xi_mesh: float = 20e-9              # [m]

    # Computed areas [m^2/particle] — populated by compute()
    external_area: float = 0.0
    internal_geometric_area: float = 0.0
    reagent_accessible_area: float = 0.0
    ligand_accessible_area: float = 0.0
    trust_level: str = "CAUTION"

    # Accessibility fractions (stored after compute for downstream reference)
    _f_reagent: float = field(default=0.0, repr=False)
    _f_ligand: float = field(default=0.0, repr=False)

    def compute(
        self,
        reagent_radius: float = _DEFAULT_REAGENT_RADIUS,
        ligand_radius: float = _DEFAULT_LIGAND_RADIUS,
        tortuosity: float = _DEFAULT_TORTUOSITY,
    ) -> None:
        """Populate area fields from geometric inputs.

        Args:
            reagent_radius: Hydrodynamic radius of small-molecule reagent [m].
            ligand_radius: Hydrodynamic radius of protein/ligand [m].
            tortuosity: Pore path tortuosity [-].  Must be >= 1.0.

        Raises:
            ValueError: If geometric inputs are non-physical.
        """
        # --- Input validation ---
        if self.bead_radius <= 0:
            raise ValueError(f"bead_radius must be positive, got {self.bead_radius}")
        if not (0.0 <= self.porosity < 1.0):
            raise ValueError(f"porosity must be in [0, 1), got {self.porosity}")
        if self.pore_diameter_mean <= 0 and self.tier != SurfaceAreaTier.EXTERNAL_ONLY:
            raise ValueError(
                f"pore_diameter_mean must be positive for tier {self.tier.value}, "
                f"got {self.pore_diameter_mean}"
            )
        if tortuosity < 1.0:
            raise ValueError(f"tortuosity must be >= 1.0, got {tortuosity}")

        R = self.bead_radius
        eps = self.porosity
        d_pore = self.pore_diameter_mean

        # --- External area: A_ext = 4 * pi * R^2 ---
        self.external_area = 4.0 * math.pi * R ** 2  # [m^2/particle]

        # --- Bead volume: V_bead = (4/3) * pi * R^3 ---
        V_bead = (4.0 / 3.0) * math.pi * R ** 3  # [m^3/particle]

        # --- Internal area depends on tier ---
        if self.tier == SurfaceAreaTier.EXTERNAL_ONLY:
            self.internal_geometric_area = 0.0
            self._f_reagent = 0.0
            self._f_ligand = 0.0
        else:
            # Specific surface area (cylindrical pore model): S_v = 4*eps/d_pore
            # Guard: d_pore > 0 already validated above for non-EXTERNAL tiers
            S_v = 4.0 * eps / d_pore  # [m^2/m^3]
            self.internal_geometric_area = S_v * V_bead  # [m^2/particle]

            # Accessibility fraction: f = max((1 - (2*r_h/d_pore)^2), 0) / tortuosity
            self._f_reagent = _accessibility_fraction(reagent_radius, d_pore, tortuosity)
            self._f_ligand = _accessibility_fraction(ligand_radius, d_pore, tortuosity)

        # --- Accessible areas ---
        self.reagent_accessible_area = (
            self.external_area + self._f_reagent * self.internal_geometric_area
        )
        self.ligand_accessible_area = (
            self.external_area + self._f_ligand * self.internal_geometric_area
        )

        # --- Trust level ---
        self.trust_level = _TIER_TRUST[self.tier]

        logger.debug(
            "AccessibleSurfaceModel.compute: tier=%s, A_ext=%.3e m^2, "
            "A_int=%.3e m^2, A_reagent=%.3e m^2, A_ligand=%.3e m^2",
            self.tier.value,
            self.external_area,
            self.internal_geometric_area,
            self.reagent_accessible_area,
            self.ligand_accessible_area,
        )

    @property
    def bead_volume(self) -> float:
        """Bead volume [m^3/particle]."""
        return (4.0 / 3.0) * math.pi * self.bead_radius ** 3

    @classmethod
    def from_m1_export(
        cls,
        contract: M1ExportContract,
        tier: SurfaceAreaTier = SurfaceAreaTier.EMPIRICAL_PORE,
        reagent_radius: float = _DEFAULT_REAGENT_RADIUS,
        ligand_radius: float = _DEFAULT_LIGAND_RADIUS,
        tortuosity: float = _DEFAULT_TORTUOSITY,
    ) -> AccessibleSurfaceModel:
        """Factory: build and compute from M1ExportContract.

        Args:
            contract: Stable M1 export interface.
            tier: Surface area estimation tier.
            reagent_radius: Small-molecule reagent radius [m].
            ligand_radius: Protein/ligand radius [m].
            tortuosity: Pore tortuosity [-].

        Returns:
            Fully computed AccessibleSurfaceModel instance.
        """
        # Auto-select tier if L2 model is empirical and user requests morphology_based
        effective_tier = tier
        if tier == SurfaceAreaTier.MORPHOLOGY_BASED:
            if contract.l2_model_tier != "mechanistic":
                logger.warning(
                    "MORPHOLOGY_BASED tier requested but L2 model is '%s'. "
                    "Falling back to EMPIRICAL_PORE.",
                    contract.l2_model_tier,
                )
                effective_tier = SurfaceAreaTier.EMPIRICAL_PORE

        model = cls(
            tier=effective_tier,
            bead_radius=contract.bead_radius,
            porosity=contract.porosity,
            pore_diameter_mean=contract.pore_size_mean,
            xi_mesh=contract.mesh_size_xi,
        )
        model.compute(
            reagent_radius=reagent_radius,
            ligand_radius=ligand_radius,
            tortuosity=tortuosity,
        )
        return model


# ─── Helper functions ───────────────────────────────────────────────────

def _accessibility_fraction(
    solute_radius: float,
    pore_diameter: float,
    tortuosity: float,
) -> float:
    """Compute accessible fraction of internal pore area for a solute.

    Uses steric exclusion model: pores smaller than 2 * solute_radius
    are inaccessible.  Remaining pores contribute with a quadratic
    penalty scaled by tortuosity.

    Args:
        solute_radius: Hydrodynamic radius of solute [m].
        pore_diameter: Mean pore diameter [m].
        tortuosity: Pore tortuosity [-].  Must be >= 1.0.

    Returns:
        Accessible fraction [-], in [0, 1].
    """
    if pore_diameter <= 0 or tortuosity <= 0:
        return 0.0
    ratio = 2.0 * solute_radius / pore_diameter
    if ratio >= 1.0:
        return 0.0
    # Quadratic steric exclusion: f = (1 - (2*r_h/d_pore)^2) / tortuosity
    f = (1.0 - ratio ** 2) / tortuosity
    return max(f, 0.0)


def _accessibility_fraction_psd(
    solute_radius: float,
    pore_diameter_mean: float,
    pore_diameter_std: float,
    tortuosity: float,
    n_bins: int = 20,
) -> float:
    """Compute PSD-weighted accessible fraction using log-normal distribution.

    v6.0: Uses pore_size_std from M1ExportContract for a log-normal
    pore-size distribution approximation. When std=0, falls back to
    mean-pore calculation.

    The log-normal is discretized into n_bins and each bin's accessibility
    is weighted by its probability density.

    Args:
        solute_radius: Hydrodynamic radius [m].
        pore_diameter_mean: Mean pore diameter [m].
        pore_diameter_std: Standard deviation of pore diameter [m].
        tortuosity: Pore tortuosity [-].
        n_bins: Number of log-normal bins.

    Returns:
        Weighted-average accessible fraction [-].
    """
    import numpy as _np

    if pore_diameter_std <= 0 or pore_diameter_mean <= 0:
        return _accessibility_fraction(solute_radius, pore_diameter_mean, tortuosity)

    # Log-normal parameters from mean and std
    cv = pore_diameter_std / pore_diameter_mean
    sigma_ln = _np.sqrt(_np.log(1 + cv ** 2))
    mu_ln = _np.log(pore_diameter_mean) - 0.5 * sigma_ln ** 2

    # Discretize: n_bins from mu-3sigma to mu+3sigma
    d_bins = _np.exp(_np.linspace(mu_ln - 3 * sigma_ln, mu_ln + 3 * sigma_ln, n_bins))
    # Log-normal PDF weights
    weights = _np.exp(-0.5 * ((_np.log(d_bins) - mu_ln) / sigma_ln) ** 2) / (d_bins * sigma_ln)
    weights /= weights.sum()  # normalize

    # Weighted-average accessibility
    f_total = 0.0
    for d_pore, w in zip(d_bins, weights):
        f_total += w * _accessibility_fraction(solute_radius, d_pore, tortuosity)

    return float(f_total)
