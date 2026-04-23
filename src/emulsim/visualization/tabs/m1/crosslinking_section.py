"""Crosslinking section (v9.0, milestone M4).

Renders the L3 crosslinker selector + concentration + T + t + UV +
live NH2-ratio validation. Per SA §B, this section applies ONLY to
AGAROSE_CHITOSAN — alginate/cellulose/PLGA skip L3 in the orchestrator.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import streamlit as st

from emulsim.reagent_library import CROSSLINKERS
from emulsim.level3_crosslinking.solver import (
    available_amine_concentration,
    recommended_crosslinker_concentration,
)


@dataclass
class CrosslinkingContext:
    crosslinker_key: str
    c_genipin_mM: float
    T_xlink_C: float
    t_xlink_h: float
    uv_intensity: float
    crosslinker: Any  # CrosslinkerProfile


def render_crosslinking_section(
    *,
    c_chitosan_pct: float,
    DDA_out: list | None = None,
) -> CrosslinkingContext:
    """Render the crosslinker widgets. Caller must be in the AGAROSE_CHITOSAN branch.

    Parameters
    ----------
    c_chitosan_pct : float
        Chitosan % w/v from the formulation section. Used for the live
        crosslinker/NH2 ratio validation.
    DDA_out : list | None
        If provided, the DDA slider value is appended to this list so the
        caller can pass it to its MaterialProperties overrides. Kept as an
        out-parameter because the DDA widget is co-located with the
        crosslinker/NH2 ratio display, not in the formulation section.

    Widget keys preserved: m1_crosslinker, m1_c_xl, m1_DDA, m1_T_xl, m1_t_xl, m1_uv.
    """
    st.subheader("Crosslinking (L3)")

    xl_keys = list(CROSSLINKERS.keys())
    xl_names = [CROSSLINKERS[k].name for k in xl_keys]
    xl_sel_name = st.selectbox(
        "Crosslinker", xl_names,
        index=xl_keys.index("genipin"),
        help="Select crosslinking agent",
        key="m1_crosslinker",
    )
    xl_sel_key = xl_keys[xl_names.index(xl_sel_name)]
    xl = CROSSLINKERS[xl_sel_key]
    st.caption(f"{xl.mechanism} | k\u2080={xl.k_xlink_0:.1e} | Score: {xl.suitability}/10")
    if xl_sel_key == "stmp":
        st.info(
            "STMP: food-grade triggerable crosslinker. Homogeneous for bead "
            "radius d50/2 < 500 \u00b5m (Thiele modulus < 1). See Appendix J.1.7 "
            "for the two-phase wet-lab protocol. **Not the same as TPP (STPP).**"
        )

    c_genipin_mM = st.slider(
        "Crosslinker Concentration (mM)", 0.5, 500.0, 10.0, step=0.5, key="m1_c_xl",
    )
    st.markdown(
        f"[View mechanism & protocol](/reagent_detail"
        f"?key={xl_sel_key}&source=crosslinkers"
        f"&T={xl.T_crosslink_default}&t={xl.t_crosslink_default}"
        f"&c={c_genipin_mM}&pH=7.4)",
    )

    DDA = st.slider(
        "DDA (degree of deacetylation)", 0.50, 0.99, 0.85, step=0.01,
        key="m1_DDA",
        help="Chitosan degree of deacetylation. Affects NH2 density and crosslinking capacity.",
    )
    if DDA_out is not None:
        DDA_out.clear()
        DDA_out.append(DDA)

    M_GlcN = 161.16
    c_chit_kg = c_chitosan_pct * 10.0
    NH2 = available_amine_concentration(c_chit_kg, DDA, M_GlcN)
    if NH2 > 0:
        ratio = c_genipin_mM / NH2
        c_rec = recommended_crosslinker_concentration(c_chit_kg, DDA, M_GlcN, target_p=0.20)
        if ratio >= 0.10:
            st.success(f"Crosslinker/NH\u2082 = {ratio:.3f} \u2014 sufficient for p \u2265 0.20")
        elif ratio >= 0.05:
            st.warning(
                f"Crosslinker/NH\u2082 = {ratio:.3f} \u2014 may be limiting. "
                f"Recommend \u2265 {c_rec:.1f} mM for target p = 0.20"
            )
        else:
            st.error(
                f"Crosslinker/NH\u2082 = {ratio:.3f} \u2014 severely limiting! "
                f"Increase to at least {c_rec:.1f} mM for target p = 0.20"
            )

    T_xlink_default = int(xl.T_crosslink_default - 273.15)
    t_xlink_default_h = max(1, int(xl.t_crosslink_default / 3600))
    T_xlink_C = st.slider(
        "Crosslinking Temperature (\u00b0C)", 0, 120,
        min(max(T_xlink_default, 0), 120), key="m1_T_xl",
    )
    t_xlink_h = st.slider(
        "Crosslinking Time (hours)", 1, 48,
        min(t_xlink_default_h, 48), key="m1_t_xl",
    )

    if xl.kinetics_model == "uv_dose":
        uv_intensity = st.slider(
            "UV Intensity (mW/cm\u00b2)", 1.0, 100.0, 20.0, step=1.0, key="m1_uv",
        )
    else:
        uv_intensity = 0.0

    return CrosslinkingContext(
        crosslinker_key=xl_sel_key,
        c_genipin_mM=c_genipin_mM,
        T_xlink_C=float(T_xlink_C),
        t_xlink_h=float(t_xlink_h),
        uv_intensity=float(uv_intensity),
        crosslinker=xl,
    )
