"""Material constants expander (v9.0, milestone M4).

Family-aware visibility for calibrated material constants. Per SA §E.2:

    always shown (any family with L1 PBE + W/O emulsion)
        K_L, Gamma_inf  — Langmuir IFT parameters
        C3              — Alopaeus breakage viscous correction

    agarose+chitosan only
        eta_intr_chit   — chitosan intrinsic viscosity
        eta_coup        — IPN coupling coefficient

Returns a dict of overrides that the runner merges into MaterialProperties.
"""

from __future__ import annotations

from typing import Any, Callable

import numpy as np
import streamlit as st

from emulsim.datatypes import PolymerFamily
from emulsim.literature_constants import ALL_CONSTANTS


def render_material_constants(
    *,
    family: PolymerFamily,
    surfactant: Any,
    crosslinker: Any | None,
    const_input: Callable,
    T_oil_C: float,
    c_span80_pct: float,
) -> dict:
    """Render the Material Constants expander and return a props_overrides dict.

    Parameters
    ----------
    family : PolymerFamily
        Drives which constants are visible.
    surfactant : SurfactantProfile
    crosslinker : CrosslinkerProfile | None
        Only used for A+C family (to pre-fill eta_coupling_recommended).
    const_input : callable
        The _const_input helper defined in visualization/app.py.
    T_oil_C : float
        Oil temperature in Celsius (for on-the-fly IFT recomputation).
    c_span80_pct : float
        Surfactant concentration in % w/v.
    """
    overrides: dict = {}
    with st.expander("Material Constants", expanded=False):
        st.caption("For each constant: use literature or enter your own.")
        _mc = st

        _kl = ALL_CONSTANTS["K_L"]
        custom_K_L = const_input(
            _mc, "K_L (Span-80 adsorption)", "K_L",
            _kl.value, "m\u00b3/mol", "Santini 2007", 0.01, 10.0, 0.05, "K_L",
        )
        _gi = ALL_CONSTANTS["Gamma_inf"]
        custom_Gamma_inf = const_input(
            _mc, "\u0393\u221e (max surface excess)", "Gamma_inf",
            _gi.value * 1e6, "\u00d710\u207b\u2076 mol/m\u00b2",
            "Santini 2007", 1.0, 10.0, 0.1, "Gamma_inf", "%.1f",
        )
        _c3 = ALL_CONSTANTS["breakage_C3"]
        custom_C3 = const_input(
            _mc, "C3 (viscous breakage)", "C3",
            _c3.value, "-", "Alopaeus 2002", 0.0, 1.0, 0.05, "C3", "%.2f",
        )

        # A+C-only constants: intrinsic viscosity + IPN coupling
        custom_eta_chit = None
        custom_eta_coup = None
        eta_coup_source = "Literature"
        if family is PolymerFamily.AGAROSE_CHITOSAN:
            _ec = ALL_CONSTANTS["eta_intr_chit"]
            custom_eta_chit = const_input(
                _mc, "[\u03b7] chitosan", "eta_chit",
                _ec.value, "mL/g", "Rinaudo 1993", 100.0, 2000.0, 50.0,
                "eta_chit", "%.0f",
            )
            xl_eta = crosslinker.eta_coupling_recommended if crosslinker is not None else -0.15
            xl_label = crosslinker.name if crosslinker is not None else "crosslinker"
            st.caption(f"Per-chemistry \u03b7 for {xl_label}: {xl_eta:+.2f}")
            custom_eta_coup = const_input(
                _mc, "\u03b7 coupling (IPN)", "eta_coup",
                xl_eta, "-", f"{xl_label} library", -0.5, 0.5, 0.05,
                "eta_coup", "%.2f",
            )
            eta_coup_source = st.session_state.get("src_eta_coup", "Literature")

    # Build the props overrides dict
    overrides["breakage_C3"] = custom_C3
    if custom_eta_coup is not None:
        overrides["eta_coupling"] = custom_eta_coup
        overrides["eta_is_custom"] = (eta_coup_source == "Custom")
    if custom_eta_chit is not None:
        overrides["eta_intr_chit"] = custom_eta_chit

    # IFT from the Langmuir-Szyszkowski model, evaluated at the selected T
    # and surfactant concentration. Formula identical to the v8.x code path.
    R_gas = 8.314
    T_ift = T_oil_C + 273.15
    c_mol_ift = c_span80_pct * 10.0 / surfactant.mw * 1000.0
    sigma_0_T = max(surfactant.sigma_0_paraffin - 1.0e-4 * (T_ift - 293.15), 0.001)
    use_gamma = custom_Gamma_inf * 1e-6
    use_KL = custom_K_L
    sigma_calc = max(
        sigma_0_T - R_gas * T_ift * use_gamma * np.log(1 + use_KL * c_mol_ift),
        1e-4,
    )
    overrides["sigma"] = sigma_calc

    return overrides
