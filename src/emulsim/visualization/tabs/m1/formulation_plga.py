"""PLGA formulation (v9.0, milestone M7).

Solvent-evaporation platform (DCM depletion). Orchestrator dispatches to
`_run_plga` which applies the grade preset and calls
`solve_solvent_evaporation`.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import streamlit as st

from emulsim.reagent_library import SURFACTANTS
from emulsim.properties.plga_defaults import PLGA_GRADE_PRESETS


@dataclass
class PLGAContext:
    phi_PLGA_0: float
    plga_grade: str
    surfactant_key: str
    c_span80_pct: float
    c_span80_vol_pct: float
    T_oil_C: float
    surfactant: Any


def render_formulation_plga(*, is_stirred: bool) -> PLGAContext:
    """Render PLGA inputs."""
    st.subheader("Formulation — PLGA (solvent evaporation)")

    grade_keys = list(PLGA_GRADE_PRESETS.keys())
    grade_display = {
        "50_50": "PLGA 50:50 (fast release, weeks)",
        "75_25": "PLGA 75:25 (moderate release, months)",
        "85_15": "PLGA 85:15 (slow release, 6-12 mo)",
        "pla": "PLA / PLLA (≥1 year, structural)",
    }
    grade_names = [grade_display.get(k, k) for k in grade_keys]
    grade_sel_name = st.selectbox(
        "PLGA grade", grade_names, index=0,
        help="Grade sets D_DCM, T_g, modulus, and degradation timescale via preset.",
        key="m1v9_plga_grade",
    )
    grade_sel_key = grade_keys[grade_names.index(grade_sel_name)]
    grade = PLGA_GRADE_PRESETS[grade_sel_key]
    st.caption(
        f"L-fraction={grade.L_fraction:.2f} | M_n={grade.M_n:.0f} g/mol | "
        f"T_g={grade.T_g_C:.0f}°C | D_DCM={grade.D_DCM:.1e} m²/s"
    )
    from emulsim.visualization.ui_links import build_reagent_link
    st.markdown(
        f"[View mechanism & protocol]({build_reagent_link(key=grade_sel_key, source='plga_grades')})"
    )

    phi_PLGA_pct = st.slider(
        "PLGA in DCM (% v/v)", 2.0, 30.0, float(grade.phi_PLGA_0_typical) * 100.0,
        step=0.5, key="m1v9_phi_plga",
        help="Initial PLGA volume fraction in the dispersed (DCM) phase. "
             "Higher values give larger, slower-to-form microspheres.",
    )
    phi_PLGA_0 = phi_PLGA_pct / 100.0

    st.markdown("---")
    st.subheader("Surfactant (PVA / Span / etc. in continuous phase)")
    surf_keys = list(SURFACTANTS.keys())
    surf_names = [SURFACTANTS[k].name for k in surf_keys]
    surf_sel_name = st.selectbox(
        "Surfactant", surf_names, index=surf_keys.index("span80"),
        key="m1v9_plga_surf",
    )
    surf_sel_key = surf_keys[surf_names.index(surf_sel_name)]
    surf = SURFACTANTS[surf_sel_key]
    st.caption(f"HLB={surf.hlb} | {surf.notes[:60]}")

    if is_stirred:
        c_span80_vol_pct = st.slider(
            "Surfactant in oil (% v/v)", 0.2, 5.0, 1.5, step=0.1, key="m1v9_plga_span_vv",
        )
        c_span80_pct = c_span80_vol_pct * 986.0 / 1000.0
        T_oil_C = st.slider("Continuous-phase T (°C)", 15, 40, 25, step=1, key="m1v9_plga_T_oil")
    else:
        c_span80_pct = st.slider(
            "Surfactant (% w/v)", 0.5, 5.0, 2.0, step=0.1, key="m1v9_plga_span_wv",
        )
        c_span80_vol_pct = 1.5
        T_oil_C = st.slider("Continuous-phase T (°C)", 15, 40, 25, step=1, key="m1v9_plga_T_oil_leg")

    return PLGAContext(
        phi_PLGA_0=float(phi_PLGA_0),
        plga_grade=grade_sel_key,
        surfactant_key=surf_sel_key,
        c_span80_pct=float(c_span80_pct),
        c_span80_vol_pct=float(c_span80_vol_pct),
        T_oil_C=float(T_oil_C),
        surfactant=surf,
    )
