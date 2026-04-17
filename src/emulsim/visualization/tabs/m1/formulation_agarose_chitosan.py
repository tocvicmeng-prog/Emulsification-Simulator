"""Agarose + Chitosan formulation (v9.0, milestone M4).

Renders the A+C-specific formulation inputs: surfactant, polymer %s,
oil T, cooling rate, pore-structure toggle (L2 model). Extracted from
tab_m1.py with widget keys preserved for session-state compatibility.

This module is invoked ONLY when family == AGAROSE_CHITOSAN.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import streamlit as st

from emulsim.reagent_library import SURFACTANTS


@dataclass
class AgaroseChitosanContext:
    c_agarose_pct: float
    c_chitosan_pct: float
    surfactant_key: str
    c_span80_pct: float        # % w/v (legacy) or v/v-derived (stirred)
    c_span80_vol_pct: float    # % v/v (stirred); 1.5 placeholder for legacy
    T_oil_C: float
    cooling_rate_Cmin: float
    l2_mode: str               # "empirical" | "ch_2d"
    grid_size: int
    surfactant: Any            # SurfactantProfile


def render_formulation_section(*, is_stirred: bool) -> AgaroseChitosanContext:
    """Render the left-column Formulation widgets.

    Widget keys preserved from v8.x:
        m1_surfactant, m1_c_agar, m1_c_chit, m1_span_vv, m1_span_wv,
        m1_T_oil, m1_T_oil_leg.
    """
    st.subheader("Formulation")

    surf_keys = list(SURFACTANTS.keys())
    surf_names = [SURFACTANTS[k].name for k in surf_keys]
    surf_sel_name = st.selectbox(
        "Surfactant", surf_names,
        index=surf_keys.index("span80"),
        help="Select surfactant for W/O emulsification",
        key="m1_surfactant",
    )
    surf_sel_key = surf_keys[surf_names.index(surf_sel_name)]
    surf = SURFACTANTS[surf_sel_key]
    st.caption(f"HLB={surf.hlb} | MW={surf.mw} g/mol | {surf.notes[:60]}")

    fc1, fc2 = st.columns(2)
    with fc1:
        c_agarose_pct = st.number_input("Agarose (% w/v)", 1.0, 10.0, 4.2, step=0.1, key="m1_c_agar")
    with fc2:
        c_chitosan_pct = st.number_input("Chitosan (% w/v)", 0.5, 5.0, 1.8, step=0.1, key="m1_c_chit")

    if is_stirred:
        c_span80_vol_pct = st.slider(
            "Span-80 in oil (% v/v)", 0.2, 5.0, 1.5, step=0.1, key="m1_span_vv",
        )
        c_span80_pct = c_span80_vol_pct * 986.0 / 1000.0
        T_oil_C = st.slider(
            "Paraffin Oil Temperature (C)", 65, 110, 80, step=1, key="m1_T_oil",
            help="Oil temperature at start of emulsification.",
        )
    else:
        c_span80_pct = st.slider(
            "Surfactant (% w/v)", 0.5, 5.0, 2.0, step=0.1, key="m1_span_wv",
        )
        c_span80_vol_pct = 1.5
        T_oil_C = st.slider(
            "Oil Temperature (C)", 60, 95, 90, key="m1_T_oil_leg",
        )

    # L2 cooling + pore model (right column in the legacy layout, but
    # logically part of A+C formulation so we include it here).
    st.subheader("Cooling & Gelation (L2)")
    if is_stirred:
        cooling_rate_Cmin = st.slider(
            "Cooling Rate (C/min)", 0.1, 15.0, 0.67, step=0.1,
            key="m1_cool_rate",
            help="Natural cooling: ~0.67 C/min for 500 mL beaker",
        )
    else:
        cooling_rate_Cmin = st.slider(
            "Cooling Rate (C/min)", 1.0, 20.0, 10.0, step=0.5, key="m1_cool_rate_leg",
        )
    l2_model = st.radio(
        "Pore Structure Model",
        ["Empirical (fast, calibrated)", "Cahn-Hilliard 2D (mechanistic)"],
        index=0, key="m1_l2_model",
        help="Empirical: literature-calibrated power law (~1 ms). "
             "CH 2D: phase-field spinodal decomposition (slower, requires parameter tuning).",
    )
    l2_mode = "empirical" if "Empirical" in l2_model else "ch_2d"
    grid_size = 64
    if l2_mode == "ch_2d":
        grid_size = st.select_slider(
            "Phase-Field Grid", [32, 64, 128], value=64, key="m1_grid",
        )

    return AgaroseChitosanContext(
        c_agarose_pct=float(c_agarose_pct),
        c_chitosan_pct=float(c_chitosan_pct),
        surfactant_key=surf_sel_key,
        c_span80_pct=float(c_span80_pct),
        c_span80_vol_pct=float(c_span80_vol_pct),
        T_oil_C=float(T_oil_C),
        cooling_rate_Cmin=float(cooling_rate_Cmin),
        l2_mode=l2_mode,
        grid_size=int(grid_size),
        surfactant=surf,
    )
