"""Dataclasses shared by every suggestion module."""

from __future__ import annotations

from dataclasses import dataclass, field


@dataclass(frozen=True)
class SuggestionContext:
    """Runtime snapshot a suggestion needs to produce text + derive a target.

    Frozen so the round-trip URL codec cannot mutate it. Every field is a
    primitive (float, int, str, bool) or a dict of primitives -> serialization
    stays straightforward.
    """

    # Family routing
    family: str                         # PolymerFamily.value string

    # L1 outputs (droplet / bead size)
    d32_actual: float                   # [m]
    d50_actual: float                   # [m]

    # L2 outputs (pore / gelation)
    pore_actual: float                  # [m]
    l2_mode: str                        # "ch_2d" | "empirical" | "ch_ternary"
    cooling_rate_effective: float       # [K/s]

    # L3 outputs (crosslinking)
    p_final: float                      # crosslink conversion [0,1]

    # L4 outputs (mechanical)
    G_DN_actual: float                  # [Pa]

    # User targets
    target_d32: float                   # [m] (converted from µm in UI)
    target_pore: float                  # [m] (converted from nm in UI)
    target_G: float                     # [Pa] (converted from kPa in UI)

    # User inputs
    rpm: float
    T_oil: float                        # [K]
    cooling_rate_input: float           # [K/s]
    c_agarose: float                    # [kg/m^3]
    c_chitosan: float                   # [kg/m^3]
    c_crosslinker_mM: float
    crosslinker_key: str

    # Material snapshot (flattened primitives only)
    rho_oil: float                      # [kg/m^3]
    mu_oil: float                       # [Pa.s]
    rho_d: float                        # [kg/m^3] dispersed phase
    cp_d: float                         # [J/(kg.K)] dispersed phase
    k_oil: float                        # [W/(m.K)] oil thermal conductivity
    h_coeff: float                      # [W/(m^2.K)] convective coeff
    T_bath: float                       # [K] cooling bath temperature
    T_gel: float                        # [K] gelation onset
    DDA: float                          # chitosan degree of deacetylation
    M_GlcN: float                       # [g/mol] chitosan repeat unit MW
    f_bridge: float                     # crosslinker bridge efficiency
    impeller_D: float                   # [m] impeller diameter
    phi_d: float                        # dispersed phase volume fraction

    # Optional diagnostic
    run_id: str = ""                    # traceability back to a specific run


@dataclass(frozen=True)
class TargetRange:
    """Numeric target the user should aim for to fix a deviation.

    `nominal` is the best estimate for an exact hit; `min`/`max` bracket the
    acceptable band given the target-tolerance (default +/-10% pore, +/-20%
    d32, +/-15% G). `unit` is a human-readable unit string for display.
    `is_qualitative_only` is True when the underlying model is QUALITATIVE_TREND
    (pore from empirical model) — in that case the numeric fields should
    NOT be displayed and `qualitative_reason` carries the "why".
    """

    nominal: float                      # best estimate (0.0 if qualitative_only)
    min: float                          # lower bound of acceptable band
    max: float                          # upper bound of acceptable band
    unit: str                           # "K/s" | "rpm" | "mol/m^3" | "%"
    limited_by: str                     # diagnostic: "biot_number" | "laminar" | "stoichiometry" | ...
    confidence_tier: str                # "VALIDATED" | "SEMI_QUANTITATIVE" | "QUALITATIVE_TREND"
    assumptions: tuple[str, ...] = ()   # human-readable assumption list
    is_qualitative_only: bool = False   # True -> numeric target withheld
    qualitative_reason: str = ""        # explanation shown when is_qualitative_only


@dataclass(frozen=True)
class Suggestion:
    """A structured, hyperlinkable optimization suggestion."""

    key: str                            # URL key; must be in REGISTRY_KEYS
    display_text: str                   # surface rendering (M1 tab bullet)
    severity: str                       # "info" | "warning" | "error"
    context: SuggestionContext          # serialisable snapshot
    extras: dict[str, str] = field(default_factory=dict)
    # extras: per-suggestion knobs (e.g. {"direction": "increase"}) that the
    # derivation page consumes without polluting SuggestionContext.


__all__ = ["SuggestionContext", "TargetRange", "Suggestion"]
