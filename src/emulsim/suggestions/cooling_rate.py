"""'adjust_cooling_rate' suggestion: generator + derivation-page renderer."""

from __future__ import annotations

from ..properties.thermal_derivation import (
    cooling_rate_for_target_pore,
)
from .generators import (
    DEVIATION_TRIGGER,
    cooling_rate_text,
    pore_deviation,
)
from .types import Suggestion, SuggestionContext, TargetRange


def generate(ctx: SuggestionContext) -> Suggestion | None:
    if pore_deviation(ctx) <= DEVIATION_TRIGGER:
        return None
    return Suggestion(
        key="adjust_cooling_rate",
        display_text=cooling_rate_text(ctx),
        severity="warning",
        context=ctx,
    )


def derive_target(ctx: SuggestionContext) -> TargetRange:
    result = cooling_rate_for_target_pore(
        target_pore=ctx.target_pore,
        d_bead=ctx.d50_actual,
        T_oil=ctx.T_oil,
        T_bath=ctx.T_bath,
        rho_d=ctx.rho_d,
        cp_d=ctx.cp_d,
        h_oil=ctx.h_coeff,
        k_oil=ctx.k_oil,
        l2_mode=ctx.l2_mode,
    )
    is_qualitative = result.confidence_tier == "QUALITATIVE_TREND"
    qual_reason = ""
    if is_qualitative:
        qual_reason = (
            f"L2 pore model is '{ctx.l2_mode}' (empirical correlation, not a "
            f"first-principles derivation). A numeric cooling-rate target would "
            f"imply a level of physical grounding the model does not support. "
            f"Direction-only guidance: slower cooling → larger pores; faster → "
            f"finer. Switch to mechanistic L2 mode ('ch_2d' or 'ch_ternary') in "
            f"Scientific Mode to unlock a numeric target."
        )
    return TargetRange(
        nominal=result.nominal if not is_qualitative else 0.0,
        min=result.min if not is_qualitative else 0.0,
        max=result.max if not is_qualitative else 0.0,
        unit="K/s",
        limited_by=result.limited_by,
        confidence_tier=result.confidence_tier,
        assumptions=result.assumptions,
        is_qualitative_only=is_qualitative,
        qualitative_reason=qual_reason,
    )


def render_derivation(ctx: SuggestionContext, target: TargetRange) -> None:
    """Streamlit-rendered three-section derivation page body."""
    import streamlit as st

    # Recompute the full CoolingRateTarget so we can expose intermediate
    # numbers on the page (tau_th, Bi, dT/dt_current).
    full = cooling_rate_for_target_pore(
        target_pore=ctx.target_pore,
        d_bead=ctx.d50_actual,
        T_oil=ctx.T_oil,
        T_bath=ctx.T_bath,
        rho_d=ctx.rho_d,
        cp_d=ctx.cp_d,
        h_oil=ctx.h_coeff,
        k_oil=ctx.k_oil,
        l2_mode=ctx.l2_mode,
    )

    st.markdown(
        "The suggestion comes from a four-step chain of physical scaling "
        "relations. Each step is reversible in closed form, which is why we "
        "can invert the chain at your target pore size to back-compute the "
        "required cooling rate."
    )

    st.subheader("Step 1 — Biot-number check (heat-transfer regime)")
    st.latex(r"\mathrm{Bi} = \dfrac{h \cdot (d/2)}{k_{\mathrm{oil}}}")
    st.markdown(
        f"- Bead diameter d = **{ctx.d50_actual*1e6:.1f} µm**\n"
        f"- Convective coefficient h = **{ctx.h_coeff:.0f} W/(m²·K)**\n"
        f"- Oil thermal conductivity k_oil = **{ctx.k_oil:.3f} W/(m·K)**\n"
        f"- → **Bi = {full.biot_number:.3f}**"
    )
    st.caption(
        "Bi < 0.1 means the bead temperature is effectively uniform — lumped "
        "capacitance applies. Bi > 1 means conduction into the bead interior "
        "is rate-limiting and this closed-form inversion is approximate."
    )

    st.subheader("Step 2 — Lumped thermal time constant")
    st.latex(r"\tau_{\mathrm{th}} = \dfrac{\rho \cdot c_p \cdot d}{6\,h}")
    st.markdown(
        f"- Dispersed-phase density ρ = **{ctx.rho_d:.0f} kg/m³**\n"
        f"- Specific heat c_p = **{ctx.cp_d:.0f} J/(kg·K)**\n"
        f"- → **τ_th = {full.tau_th:.1f} s**"
    )

    st.subheader("Step 3 — Current effective cooling rate")
    st.latex(r"\left| \dfrac{dT}{dt} \right| \approx \dfrac{T_{\mathrm{oil}} - T_{\mathrm{bath}}}{\tau_{\mathrm{th}}}")
    st.markdown(
        f"- T_oil = **{ctx.T_oil - 273.15:.1f} °C**\n"
        f"- T_bath = **{ctx.T_bath - 273.15:.1f} °C**\n"
        f"- → Current |dT/dt| = **{full.dT_dt_effective:.3f} K/s**\n"
        f"- Your selected cooling rate (input) = **{ctx.cooling_rate_input:.3f} K/s**"
    )

    st.subheader("Step 4 — Spinodal dwell and pore scaling")
    st.latex(r"\Delta t_{\mathrm{spin}} = \dfrac{2\,\Delta T_{\mathrm{band}}}{|dT/dt|}")
    st.latex(r"\mathrm{pore} \sim \sqrt{M \kappa} \cdot \Delta t_{\mathrm{spin}}^{1/2}")
    st.markdown(
        f"- Actual pore size (from this run) = **{ctx.pore_actual*1e9:.0f} nm**\n"
        f"- Target pore size = **{ctx.target_pore*1e9:.0f} nm**\n"
        f"- L2 model = `{ctx.l2_mode}`"
    )

    st.subheader("Step 5 — Inverted: target pore → required cooling rate")
    if target.is_qualitative_only:
        st.warning(
            "**Numeric cooling-rate target withheld.**\n\n" + target.qualitative_reason
        )
    else:
        st.latex(r"\left| \dfrac{dT}{dt} \right|_{\mathrm{req}} = \dfrac{2\,\Delta T_{\mathrm{band}}}{\left(\mathrm{pore}/\sqrt{M\kappa}\right)^{2}}")
        col1, col2, col3 = st.columns(3)
        col1.metric("Nominal |dT/dt|", f"{target.nominal:.3f} K/s")
        col2.metric("Lower bound", f"{target.min:.3f} K/s",
                    help="Corresponds to pore + 10% tolerance")
        col3.metric("Upper bound", f"{target.max:.3f} K/s",
                    help="Corresponds to pore − 10% tolerance")
