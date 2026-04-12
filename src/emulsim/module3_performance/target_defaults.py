"""Target protein default parameters for M3 process simulation.

v5.9.0: Provides default molecular properties for common purification
targets. These are used when the user does not specify target-specific
parameters in M3.

All values are ESTIMATES and should be treated as semi-quantitative
defaults. User calibration is required for quantitative prediction.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass
class TargetProtein:
    """Default molecular parameters for a purification target."""
    name: str
    mw: float                   # [Da]
    r_h: float                  # [m] hydrodynamic radius
    pI: float = 7.0             # isoelectric point
    z_eff: float = 5.0          # effective charge for SMA (at working pH)
    sigma_sma: float = 10.0     # SMA steric factor
    K_affinity: float = 1e6     # [m3/mol] generic affinity constant
    tag_type: str = ""          # "His6", "GST", "biotin", "Fc", ""
    confidence: str = "estimated"
    notes: str = ""


TARGET_PROTEIN_LIBRARY: dict[str, TargetProtein] = {

    "IgG": TargetProtein(
        name="Human IgG1 (generic)",
        mw=150000.0,
        r_h=5.3e-9,
        pI=7.5,
        z_eff=5.0,
        sigma_sma=50.0,
        K_affinity=1e8,
        tag_type="Fc",
        confidence="estimated",
        notes="Generic human IgG1; pI range 6-9 depending on clone",
    ),

    "His6_50kDa": TargetProtein(
        name="His6-tagged protein (50 kDa)",
        mw=50000.0,
        r_h=2.8e-9,
        K_affinity=1e4,
        tag_type="His6",
        confidence="estimated",
        notes="Generic His-tagged protein; K for Ni-NTA IMAC",
    ),

    "GST_fusion": TargetProtein(
        name="GST fusion protein (75 kDa)",
        mw=75000.0,
        r_h=3.5e-9,
        K_affinity=1e3,
        tag_type="GST",
        confidence="estimated",
        notes="GST-tagged; K for glutathione affinity",
    ),

    "biotin_tagged": TargetProtein(
        name="Biotinylated protein (50 kDa)",
        mw=50000.0,
        r_h=2.8e-9,
        K_affinity=1e10,
        tag_type="biotin",
        confidence="estimated",
        notes="Near-irreversible streptavidin-biotin; Kd ~10^-15 M",
    ),

    "BSA": TargetProtein(
        name="Bovine Serum Albumin (reference)",
        mw=66500.0,
        r_h=3.6e-9,
        pI=4.7,
        z_eff=-12.0,
        sigma_sma=20.0,
        K_affinity=0.0,
        confidence="literature",
        notes="Standard reference protein for IEX and SEC",
    ),

    "lysozyme": TargetProtein(
        name="Hen egg-white lysozyme (reference)",
        mw=14300.0,
        r_h=1.9e-9,
        pI=11.0,
        z_eff=8.0,
        sigma_sma=5.0,
        K_affinity=0.0,
        confidence="literature",
        notes="Standard reference for cation exchange",
    ),
}
