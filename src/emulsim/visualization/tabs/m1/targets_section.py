"""Targets section (v9.0, milestone M4).

Family-aware optimization target inputs (d32 / d_mode / pore / G_DN).
Per SA §E.1 the physical meaning of pore_size differs per family, so the
help text is family-aware while the widget keys are preserved for
session-state compatibility.

M4 extracts the A+C targets from tab_m1.py. M5-M7 will add per-family
label overrides.
"""

from __future__ import annotations

from dataclasses import dataclass

import streamlit as st

from emulsim.datatypes import PolymerFamily


@dataclass
class TargetsContext:
    target_d32: float       # [um] volume-surface mean (legacy mode)
    target_d_mode: float    # [um] modal diameter (stirred-vessel mode)
    target_pore: float      # [nm] pore / mesh / free-volume scale
    target_G: float         # [kPa] target double-network shear modulus


_PORE_HELP = {
    PolymerFamily.AGAROSE_CHITOSAN: "Mean pore size from Cahn-Hilliard or empirical pore model.",
    PolymerFamily.ALGINATE: "Ionic-gel mesh size (gel-front depth at matched t).",
    PolymerFamily.CELLULOSE: "NIPS spinodal wavelength (characteristic pore).",
    PolymerFamily.PLGA: "Glassy polymer free-volume scale (SAXS-equivalent).",
}


def render_targets_section(*, family: PolymerFamily, is_stirred: bool) -> TargetsContext:
    """Render optimization targets.

    Keys preserved from v8.x:
        m1_tgt_d      (stirred vessel d_mode)
        m1_tgt_d32    (legacy d32)
        m1_tgt_pore   (stirred vessel pore)
        m1_tgt_pore_leg (legacy pore)
        m1_tgt_G      (G_DN)
    """
    st.subheader("Optimization Targets")
    pore_help = _PORE_HELP.get(family, "Characteristic pore / mesh size.")
    if is_stirred:
        target_d_mode = st.number_input(
            "Target d_mode (um)", 10.0, 500.0, 100.0, step=10.0,
            key="m1_tgt_d", help="Modal diameter of microspheres",
        )
        target_d32 = target_d_mode
        target_pore = st.number_input(
            "Target Pore Size (nm)", 10, 500, 100, step=10,
            key="m1_tgt_pore", help=pore_help,
        )
    else:
        target_d32 = st.number_input(
            "Target d32 (um)", 0.5, 50.0, 2.0, step=0.5, key="m1_tgt_d32",
        )
        target_d_mode = target_d32
        target_pore = st.number_input(
            "Target Pore Size (nm)", 10, 500, 80, step=10,
            key="m1_tgt_pore_leg", help=pore_help,
        )
    target_G = st.number_input(
        "Target G_DN (kPa)", 1.0, 500.0, 10.0, step=1.0, key="m1_tgt_G",
    )
    return TargetsContext(
        target_d32=target_d32,
        target_d_mode=target_d_mode,
        target_pore=target_pore,
        target_G=target_G,
    )
