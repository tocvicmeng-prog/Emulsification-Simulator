"""Alginate formulation (v9.0, milestone M5).

Renders ionotropic Ca²⁺ gelation inputs. L3 crosslinking does NOT apply
(ionic gelation IS the crosslinking); orchestrator dispatches to
`_run_alginate` which calls `solve_ionic_ca_gelation`.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import streamlit as st

from emulsim.reagent_library import SURFACTANTS
from emulsim.reagent_library_alginate import GELANTS_ALGINATE


@dataclass
class AlginateContext:
    c_alginate_kg_m3: float
    gelant_key: str
    c_Ca_bath_mM: float       # [mol/m³] used when mode=external_bath
    C_Ca_source_mM: float     # [mol/m³] used when mode=internal_release
    k_release_1_s: float      # [1/s]
    T_bath_C: float
    surfactant_key: str
    c_span80_pct: float
    c_span80_vol_pct: float
    T_oil_C: float
    surfactant: Any


def render_formulation_alginate(*, is_stirred: bool) -> AlginateContext:
    """Render alginate + gelant + surfactant inputs."""
    st.subheader("Formulation — Alginate")

    c_alginate_pct = st.number_input(
        "Alginate (% w/v)", 0.5, 5.0, 1.5, step=0.1, key="m1v9_c_alg",
        help="Typical sodium-alginate microsphere recipes: 1-3 % w/v.",
    )

    gelant_keys = list(GELANTS_ALGINATE.keys())
    gelant_names = [GELANTS_ALGINATE[k].name for k in gelant_keys]
    gel_sel_name = st.selectbox(
        "Gelant", gelant_names,
        index=0,
        help="External bath: pre-formed droplets fall into CaCl₂ (shrinking-core). "
             "Internal release: GDL+CaCO₃ dispersed in the alginate phase, homogeneous gel.",
        key="m1v9_alg_gelant",
    )
    gel_sel_key = gelant_keys[gelant_names.index(gel_sel_name)]
    gelant = GELANTS_ALGINATE[gel_sel_key]
    st.caption(f"Mode: {gelant.mode} | T_default={gelant.T_default-273.15:.0f}°C | Suitability: {gelant.suitability}/10")
    from emulsim.visualization.ui_links import build_reagent_link
    st.markdown(
        f"[View mechanism & protocol]({build_reagent_link(key=gel_sel_key, source='alginate_gelants', T_K=gelant.T_default, t_s=gelant.t_default, c_mM=gelant.C_Ca_bath)})"
    )

    if gelant.mode == "external_bath":
        c_Ca_bath_mM = st.slider(
            "External CaCl₂ bath (mM)", 20.0, 500.0, float(gelant.C_Ca_bath),
            step=5.0, key="m1v9_c_Ca_bath",
        )
        C_Ca_source_mM = 0.0
        k_release = 0.0
    else:
        c_Ca_bath_mM = 0.0
        C_Ca_source_mM = st.slider(
            "Internal Ca²⁺ source (mM, CaCO₃ equivalent)",
            5.0, 100.0, float(gelant.C_Ca_source), step=1.0, key="m1v9_C_Ca_source",
        )
        k_release = st.number_input(
            "GDL release rate k (1/s)", 1e-5, 1e-2,
            float(gelant.k_release), format="%.1e", key="m1v9_k_release",
        )

    T_bath_C = st.slider(
        "Bath temperature (°C)", 4, 80, int(gelant.T_default - 273.15), step=1,
        key="m1v9_T_bath",
    )

    st.markdown("---")
    st.subheader("Surfactant")
    surf_keys = list(SURFACTANTS.keys())
    surf_names = [SURFACTANTS[k].name for k in surf_keys]
    surf_sel_name = st.selectbox(
        "Surfactant", surf_names, index=surf_keys.index("span80"),
        key="m1v9_alg_surf",
    )
    surf_sel_key = surf_keys[surf_names.index(surf_sel_name)]
    surf = SURFACTANTS[surf_sel_key]
    st.caption(f"HLB={surf.hlb} | {surf.notes[:60]}")

    if is_stirred:
        c_span80_vol_pct = st.slider(
            "Surfactant in oil (% v/v)", 0.2, 5.0, 1.5, step=0.1, key="m1v9_alg_span_vv",
        )
        c_span80_pct = c_span80_vol_pct * 986.0 / 1000.0
        T_oil_C = st.slider("Continuous-phase T (°C)", 20, 60, 25, step=1, key="m1v9_alg_T_oil")
    else:
        c_span80_pct = st.slider(
            "Surfactant (% w/v)", 0.5, 5.0, 2.0, step=0.1, key="m1v9_alg_span_wv",
        )
        c_span80_vol_pct = 1.5
        T_oil_C = st.slider("Continuous-phase T (°C)", 20, 60, 25, step=1, key="m1v9_alg_T_oil_leg")

    return AlginateContext(
        c_alginate_kg_m3=float(c_alginate_pct * 10.0),
        gelant_key=gel_sel_key,
        c_Ca_bath_mM=float(c_Ca_bath_mM),
        C_Ca_source_mM=float(C_Ca_source_mM),
        k_release_1_s=float(k_release),
        T_bath_C=float(T_bath_C),
        surfactant_key=surf_sel_key,
        c_span80_pct=float(c_span80_pct),
        c_span80_vol_pct=float(c_span80_vol_pct),
        T_oil_C=float(T_oil_C),
        surfactant=surf,
    )
