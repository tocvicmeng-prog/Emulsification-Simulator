"""Cellulose NIPS formulation (v9.0, milestone M6).

Ternary Cahn-Hilliard + Fickian solver via `solve_nips_cellulose`.
L2a timing, L3, cooling-rate model toggle all skipped. Cooling rate
visible only for NMMO (thermal solvent); other presets are isothermal
water-bath coagulation.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import streamlit as st

from emulsim.reagent_library import SURFACTANTS
from emulsim.properties.cellulose_defaults import CELLULOSE_SOLVENT_PRESETS


_THERMAL_SOLVENTS = {"nmmo"}  # requires hot coagulation bath


@dataclass
class CelluloseContext:
    phi_cellulose_0: float
    solvent_system: str
    cooling_rate_Cmin: float  # ignored for isothermal solvents
    surfactant_key: str
    c_span80_pct: float
    c_span80_vol_pct: float
    T_oil_C: float
    surfactant: Any


def render_formulation_cellulose(*, is_stirred: bool) -> CelluloseContext:
    """Render cellulose NIPS inputs."""
    st.subheader("Formulation — Cellulose (NIPS)")

    phi_cell_pct = st.slider(
        "Cellulose initial volume fraction (%)", 1.0, 20.0, 5.0, step=0.5,
        key="m1v9_phi_cell",
        help="Initial cellulose volume fraction in the droplet prior to non-solvent ingress.",
    )
    phi_cellulose_0 = phi_cell_pct / 100.0

    sol_keys = list(CELLULOSE_SOLVENT_PRESETS.keys())
    sol_display = {
        "naoh_urea": "NaOH / urea (aqueous, isothermal)",
        "nmmo": "NMMO (Lyocell, ~90 °C)",
        "emim_ac": "EMIM-Ac (ionic liquid)",
        "dmac_licl": "DMAc / LiCl",
    }
    sol_names = [sol_display.get(k, k) for k in sol_keys]
    sol_sel_display = st.selectbox(
        "Solvent system", sol_names, index=0,
        help="NIPS demixing kinetics + polymer-solvent χ parameters come from the preset.",
        key="m1v9_cell_solvent",
    )
    sol_sel_key = sol_keys[sol_names.index(sol_sel_display)]
    from emulsim.visualization.ui_links import build_reagent_link
    st.markdown(
        f"[View mechanism & protocol]({build_reagent_link(key=sol_sel_key, source='cellulose_solvents')})"
    )

    is_thermal = sol_sel_key in _THERMAL_SOLVENTS
    if is_thermal:
        cooling_rate_Cmin = st.slider(
            "Coagulation-bath cooling rate (°C/min)", 0.1, 15.0, 2.0, step=0.1,
            key="m1v9_cell_cool_rate",
            help="Applies to thermal solvents (NMMO) where demixing is triggered by quench.",
        )
    else:
        cooling_rate_Cmin = 0.0
        st.caption("Isothermal solvent — no cooling-rate input.")

    st.markdown("---")
    st.subheader("Surfactant")
    surf_keys = list(SURFACTANTS.keys())
    surf_names = [SURFACTANTS[k].name for k in surf_keys]
    surf_sel_name = st.selectbox(
        "Surfactant", surf_names, index=surf_keys.index("span80"),
        key="m1v9_cell_surf",
    )
    surf_sel_key = surf_keys[surf_names.index(surf_sel_name)]
    surf = SURFACTANTS[surf_sel_key]
    st.caption(f"HLB={surf.hlb} | {surf.notes[:60]}")

    if is_stirred:
        c_span80_vol_pct = st.slider(
            "Surfactant in oil (% v/v)", 0.2, 5.0, 1.5, step=0.1, key="m1v9_cell_span_vv",
        )
        c_span80_pct = c_span80_vol_pct * 986.0 / 1000.0
        T_oil_C = st.slider(
            "Continuous-phase T (°C)", 20, 100, 90 if is_thermal else 25,
            step=1, key="m1v9_cell_T_oil",
        )
    else:
        c_span80_pct = st.slider(
            "Surfactant (% w/v)", 0.5, 5.0, 2.0, step=0.1, key="m1v9_cell_span_wv",
        )
        c_span80_vol_pct = 1.5
        T_oil_C = st.slider(
            "Continuous-phase T (°C)", 20, 100, 90 if is_thermal else 25,
            step=1, key="m1v9_cell_T_oil_leg",
        )

    return CelluloseContext(
        phi_cellulose_0=float(phi_cellulose_0),
        solvent_system=sol_sel_key,
        cooling_rate_Cmin=float(cooling_rate_Cmin),
        surfactant_key=surf_sel_key,
        c_span80_pct=float(c_span80_pct),
        c_span80_vol_pct=float(c_span80_vol_pct),
        T_oil_C=float(T_oil_C),
        surfactant=surf,
    )
