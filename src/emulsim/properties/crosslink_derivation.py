"""Crosslinker + polymer physics for increase_crosslinker / reduce_polymer.

Rubber elasticity (Gaussian network):
    G = nu_e * R * T
    nu_e = c_chain * p_final * f_bridge * DDA / M_GlcN       [chitosan]

Inverse (increase_crosslinker): given target_G_chit, solve for the
crosslinker concentration that lands p_final at the required value.
This wraps level3_crosslinking.solver.recommended_crosslinker_concentration
for the target-p computation.

Inverse (reduce_polymer): given target_G_DN < actual_G_DN, derive a
uniform scaling factor alpha such that c_agarose' = alpha*c_agarose and
c_chitosan' = alpha*c_chitosan produces the target G_DN. Assumes the
phenomenological double-network relation G_DN ~ alpha * (G_agarose_0 +
G_chitosan_0 + eta*sqrt(G_A*G_C)) -- a linear scaling in alpha.

Confidence tier: SEMI_QUANTITATIVE for both paths (rubber elasticity is
textbook; the alpha scaling is an approximation that will over- or
under-predict by tens of percent at extreme dilutions).

References:
    Treloar (1975) The Physics of Rubber Elasticity.
    Peppas & Merrill (1977) J. Polym. Sci. 15:739 (swollen network G).
    Flory & Rehner (1943) J. Chem. Phys. 11:521.
"""

from __future__ import annotations

from dataclasses import dataclass

from ..level3_crosslinking.solver import (
    available_amine_concentration,
    recommended_crosslinker_concentration,
)

_R_GAS: float = 8.314462618  # J/(mol.K)


@dataclass(frozen=True)
class CrosslinkerTarget:
    """Target crosslinker concentration to hit a G_chit target."""

    nominal: float          # [mol/m^3]
    min: float              # [mol/m^3] -15% G tolerance
    max: float              # [mol/m^3] +15% G tolerance
    required_p: float       # crosslink fraction that achieves target G
    stoichiometric_ceiling: float  # max p given [NH2]/2
    limited_by: str         # "ok" | "stoichiometry" | "kinetics"
    confidence_tier: str
    assumptions: tuple[str, ...]


@dataclass(frozen=True)
class PolymerScaleTarget:
    """Polymer concentration scaling factor alpha for a target G_DN."""

    alpha_nominal: float    # dimensionless multiplier on (c_agarose, c_chitosan)
    alpha_min: float        # -15% G
    alpha_max: float        # +15% G
    new_c_agarose: float    # [kg/m^3] alpha * c_agarose_current
    new_c_chitosan: float   # [kg/m^3] alpha * c_chitosan_current
    limited_by: str         # "ok" | "under_gelation" | "solubility"
    confidence_tier: str
    assumptions: tuple[str, ...]


def p_for_target_G_chit(
    target_G_chit: float,
    c_chitosan: float,
    DDA: float,
    M_GlcN: float,
    f_bridge: float,
    T: float,
) -> float:
    """Required crosslink fraction p given G_chit = nu_e * R * T."""
    if T <= 0:
        raise ValueError(f"T must be positive, got {T}")
    if c_chitosan <= 0 or f_bridge <= 0 or DDA <= 0:
        raise ValueError("c_chitosan, f_bridge, DDA must all be positive")
    # nu_e = target_G_chit / (R*T)
    nu_e = target_G_chit / (_R_GAS * T)
    # nu_e = c_chit * p * f_bridge * DDA / M_GlcN (kg/m^3 * 1000 g/kg -> g/m^3)
    nh2 = available_amine_concentration(c_chitosan, DDA, M_GlcN)
    # nu_e = nh2 * p * f_bridge (since nh2 already bakes in c_chit*DDA/M_GlcN)
    # Actually: nu_e (mol/m^3) = NH2 (mol/m^3) * p * f_bridge * 0.5
    # (0.5 because each bridge links 2 chains -> half as many "crosslink nodes"
    # per reacted amine site).
    denom = nh2 * f_bridge * 0.5
    if denom <= 0:
        return float("inf")
    return nu_e / denom


def crosslinker_conc_for_target_G(
    *,
    target_G_chit: float,
    c_chitosan: float,
    DDA: float,
    M_GlcN: float,
    f_bridge: float,
    T: float,
    G_tolerance: float = 0.15,
) -> CrosslinkerTarget:
    """Derive crosslinker concentration [mol/m^3] for target G_chit."""
    p_nom = p_for_target_G_chit(target_G_chit, c_chitosan, DDA, M_GlcN, f_bridge, T)
    p_min = p_for_target_G_chit(target_G_chit * (1.0 - G_tolerance),
                                 c_chitosan, DDA, M_GlcN, f_bridge, T)
    p_max = p_for_target_G_chit(target_G_chit * (1.0 + G_tolerance),
                                 c_chitosan, DDA, M_GlcN, f_bridge, T)

    c_nom = recommended_crosslinker_concentration(c_chitosan, DDA, M_GlcN, target_p=p_nom)
    c_min = recommended_crosslinker_concentration(c_chitosan, DDA, M_GlcN, target_p=p_min)
    c_max = recommended_crosslinker_concentration(c_chitosan, DDA, M_GlcN, target_p=p_max)

    # Stoichiometric ceiling: at p=1 every NH2 is reacted.
    ceiling = available_amine_concentration(c_chitosan, DDA, M_GlcN) / 2.0

    assumptions = [
        "Rubber elasticity: G_chit = nu_e * R * T (Gaussian affine network).",
        "Chain density nu_e = 0.5 * [NH2] * p * f_bridge.",
        "Stoichiometry: 1 bridge consumes 1 crosslinker + 2 NH2.",
        "Ignores pre-existing physical entanglements and chain-end defects.",
    ]
    if p_nom > 1.0:
        limited_by = "stoichiometry"
        assumptions.append(
            f"WARN: required p={p_nom:.2f} exceeds 1.0 -- target G_chit cannot be "
            f"reached by crosslinker alone at this polymer concentration; "
            f"consider also raising c_chitosan."
        )
    elif p_nom > 0.5:
        limited_by = "kinetics"
        assumptions.append(
            f"NOTE: required p={p_nom:.2f} is high; second-order kinetics slow "
            f"markedly near p=0.5. Realistic reaction times may need to roughly "
            f"double to hit this conversion."
        )
    else:
        limited_by = "ok"

    return CrosslinkerTarget(
        nominal=c_nom,
        min=min(c_min, c_max),
        max=max(c_min, c_max),
        required_p=p_nom,
        stoichiometric_ceiling=ceiling,
        limited_by=limited_by,
        confidence_tier="SEMI_QUANTITATIVE",
        assumptions=tuple(assumptions),
    )


def alpha_for_target_G_DN(
    *,
    G_DN_current: float,
    target_G_DN: float,
    c_agarose_current: float,
    c_chitosan_current: float,
    G_tolerance: float = 0.15,
) -> PolymerScaleTarget:
    """Scaling factor alpha on (c_agarose, c_chitosan) for target G_DN.

    Uses the phenomenological scaling G_DN ~ c^n with n ~ 2 for
    gel concentration-modulus relationship in the semi-dilute regime
    (Ferry, de Gennes). Inverting: alpha = (target_G_DN / G_DN_current)^(1/n).

    The relation saturates at low concentration near the gel point and is
    non-physical beyond solubility limits on the upper end.
    """
    if G_DN_current <= 0 or target_G_DN <= 0:
        raise ValueError("G_DN_current and target_G_DN must be positive")

    n_exponent = 2.0  # semi-dilute scaling
    alpha_nom = (target_G_DN / G_DN_current) ** (1.0 / n_exponent)
    alpha_lo = (target_G_DN * (1.0 - G_tolerance) / G_DN_current) ** (1.0 / n_exponent)
    alpha_hi = (target_G_DN * (1.0 + G_tolerance) / G_DN_current) ** (1.0 / n_exponent)

    new_cA = alpha_nom * c_agarose_current
    new_cC = alpha_nom * c_chitosan_current

    assumptions = [
        "Semi-dilute scaling G ~ c^2 (Ferry, de Gennes).",
        "Uniform alpha on both polymer concentrations -- ratio c_agarose/c_chitosan preserved.",
        "Valid above the gel threshold (typically c_agarose > 5 kg/m^3).",
    ]
    if new_cA < 5.0:
        limited_by = "under_gelation"
        assumptions.append(
            f"WARN: scaled c_agarose={new_cA:.1f} kg/m^3 is near or below the "
            f"gel threshold (~5 kg/m^3); the beads may fail to form a "
            f"continuous network at all."
        )
    elif new_cA > 100.0 or new_cC > 50.0:
        limited_by = "solubility"
        assumptions.append(
            f"WARN: scaled concentrations (agarose {new_cA:.1f}, chitosan "
            f"{new_cC:.1f}) exceed typical solubility limits for hot-mixed "
            f"double-network systems."
        )
    else:
        limited_by = "ok"

    return PolymerScaleTarget(
        alpha_nominal=alpha_nom,
        alpha_min=min(alpha_lo, alpha_hi),
        alpha_max=max(alpha_lo, alpha_hi),
        new_c_agarose=new_cA,
        new_c_chitosan=new_cC,
        limited_by=limited_by,
        confidence_tier="SEMI_QUANTITATIVE",
        assumptions=tuple(assumptions),
    )


# Re-export existing helpers so the suggestion modules can import from one place.
__all__ = [
    "CrosslinkerTarget",
    "PolymerScaleTarget",
    "p_for_target_G_chit",
    "crosslinker_conc_for_target_G",
    "alpha_for_target_G_DN",
    "available_amine_concentration",
    "recommended_crosslinker_concentration",
]
