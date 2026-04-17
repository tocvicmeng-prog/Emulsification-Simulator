"""Polymer-family selector (v9.0, milestone M3).

Renders a 4-way radio for AGAROSE_CHITOSAN / ALGINATE / CELLULOSE / PLGA
at the top of the M1 tab. The selected family becomes
FormulationParameters.polymer_family and drives every downstream
conditional-rendering decision in tab_m1 and its per-family modules.

Scientific caption is provided per option (per scientific-advisor §B)
so the user knows which L2 branch each choice dispatches to.

Milestone matrix:
    M3 (this file)  selector live; A+C formulation unchanged (still in tab_m1)
    M4              A+C formulation extracted to formulation_agarose_chitosan
    M5              alginate branch populated (formulation_alginate)
    M6              cellulose branch populated
    M7              PLGA branch populated
"""

from __future__ import annotations

from dataclasses import dataclass

import streamlit as st

from emulsim.datatypes import PolymerFamily


_FAMILY_DISPLAY = [
    ("Agarose + Chitosan", PolymerFamily.AGAROSE_CHITOSAN,
     "Thermal-TIPS gelation with optional covalent crosslinking "
     "(genipin / glutaraldehyde / EDC-NHS / others). Legacy platform."),
    ("Alginate", PolymerFamily.ALGINATE,
     "Ionotropic Ca²⁺ crosslinking (external CaCl₂ bath or internal "
     "GDL + CaCO₃ release). L3 crosslinking step does not apply."),
    ("Cellulose (NIPS)", PolymerFamily.CELLULOSE,
     "Non-solvent-induced phase separation in NaOH/urea, NMMO, "
     "EMIM-Ac, or DMAc/LiCl. L3 crosslinking step does not apply."),
    ("PLGA (solvent evaporation)", PolymerFamily.PLGA,
     "Glassy microspheres via DCM-depletion-driven solvent evaporation. "
     "L3 crosslinking step does not apply; grade presets 50:50/75:25/85:15/PLA."),
]


@dataclass
class FamilyContext:
    """Shared context emitted by the family selector."""

    family: PolymerFamily
    display_name: str


def render_family_selector(*, key: str = "m1v9_polymer_family") -> FamilyContext:
    """Render the polymer-family radio at the top of the M1 tab.

    Default is AGAROSE_CHITOSAN (legacy behaviour, session-state compatible).
    """
    st.subheader("Polymer Family")
    display_names = [row[0] for row in _FAMILY_DISPLAY]
    enums = [row[1] for row in _FAMILY_DISPLAY]
    helps = [row[2] for row in _FAMILY_DISPLAY]

    sel_name = st.radio(
        "Polymer family",
        display_names,
        index=0,
        horizontal=True,
        key=key,
        help="Selects the L2 gelation pathway and the set of scientifically "
             "applicable formulation inputs. Only fields that enter the "
             "chosen family's equations will be shown below.",
        label_visibility="collapsed",
    )
    idx = display_names.index(sel_name)
    family = enums[idx]
    st.caption(f"**{sel_name}** — {helps[idx]}")
    return FamilyContext(family=family, display_name=sel_name)
