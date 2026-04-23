"""'increase_crosslinker' suggestion module."""

from __future__ import annotations

from ..properties.crosslink_derivation import (
    crosslinker_conc_for_target_G,
)
from .generators import (
    DEVIATION_TRIGGER,
    G_deviation,
    crosslinker_up_text,
)
from .types import Suggestion, SuggestionContext, TargetRange


def generate(ctx: SuggestionContext) -> Suggestion | None:
    if G_deviation(ctx) <= DEVIATION_TRIGGER:
        return None
    if ctx.G_DN_actual >= ctx.target_G:
        return None  # G is above target -> not the right fix
    return Suggestion(
        key="increase_crosslinker",
        display_text=crosslinker_up_text(ctx),
        severity="warning",
        context=ctx,
    )


def derive_target(ctx: SuggestionContext) -> TargetRange:
    # Crude split of G_DN into chitosan contribution: the phenomenological
    # G_DN is G_A + G_C + eta*sqrt(G_A*G_C). Treat G_target as the aggregate;
    # because chitosan carries the crosslinker-tunable portion, target G_chit
    # = target_G_DN - G_agarose_baseline. Use a conservative 30 % baseline.
    G_chit_target = 0.7 * ctx.target_G
    T_ref = 310.15  # 37 °C reference for rubber-elasticity evaluation
    result = crosslinker_conc_for_target_G(
        target_G_chit=G_chit_target,
        c_chitosan=ctx.c_chitosan,
        DDA=ctx.DDA,
        M_GlcN=ctx.M_GlcN,
        f_bridge=ctx.f_bridge,
        T=T_ref,
    )
    # Unit is mol/m^3; convert to mM for UI (1 mol/m^3 = 1 mM).
    return TargetRange(
        nominal=result.nominal,
        min=result.min,
        max=result.max,
        unit="mol/m³ (≈ mM)",
        limited_by=result.limited_by,
        confidence_tier=result.confidence_tier,
        assumptions=result.assumptions,
        is_qualitative_only=False,
        qualitative_reason="",
    )


def render_derivation(ctx: SuggestionContext, target: TargetRange) -> None:
    import streamlit as st

    from ..properties.crosslink_derivation import (
        available_amine_concentration,
        crosslinker_conc_for_target_G,
    )

    G_chit_target = 0.7 * ctx.target_G
    T_ref = 310.15
    full = crosslinker_conc_for_target_G(
        target_G_chit=G_chit_target,
        c_chitosan=ctx.c_chitosan,
        DDA=ctx.DDA,
        M_GlcN=ctx.M_GlcN,
        f_bridge=ctx.f_bridge,
        T=T_ref,
    )
    nh2 = available_amine_concentration(ctx.c_chitosan, ctx.DDA, ctx.M_GlcN)

    st.markdown(
        "The crosslinker suggestion comes from inverting the rubber-elasticity "
        "equation for the chitosan network against your target modulus."
    )

    st.subheader("Step 1 — Rubber elasticity: G → effective chain density")
    st.latex(r"G_{\mathrm{chit}} = \nu_e \cdot R \cdot T")
    st.markdown(
        f"- Target G_DN = **{ctx.target_G/1000:.1f} kPa**\n"
        f"- Chitosan-network target (70 % of G_DN) = **{G_chit_target/1000:.1f} kPa**\n"
        f"- Required ν_e = **{G_chit_target/(8.314*T_ref):.1f} mol/m³**"
    )

    st.subheader("Step 2 — Chain density → required conversion p")
    st.latex(r"\nu_e = \tfrac{1}{2} \cdot [\mathrm{NH}_2] \cdot p \cdot f_{\mathrm{bridge}}")
    st.markdown(
        f"- Available amines [NH₂] = **{nh2:.1f} mol/m³**\n"
        f"  (= c_chitosan × DDA / M_GlcN = {ctx.c_chitosan:.1f} × {ctx.DDA:.2f} / {ctx.M_GlcN:.2f})\n"
        f"- Bridge efficiency f_bridge = **{ctx.f_bridge:.2f}**\n"
        f"- Current conversion p (from run) = **{ctx.p_final:.1%}**\n"
        f"- Required conversion p_req = **{full.required_p:.1%}**"
    )

    st.subheader("Step 3 — Conversion → crosslinker concentration")
    st.latex(r"c_{\mathrm{xlink}} = \tfrac{1}{2} \cdot p \cdot [\mathrm{NH}_2]")
    st.markdown(
        f"- Current crosslinker = **{ctx.c_crosslinker_mM:.1f} mM**\n"
        f"- Stoichiometric ceiling (p=1) = **{full.stoichiometric_ceiling:.1f} mol/m³**"
    )

    st.subheader("Step 4 — Target concentration + band")
    col1, col2, col3 = st.columns(3)
    col1.metric("Nominal", f"{target.nominal:.1f} {target.unit}")
    col2.metric("Lower bound", f"{target.min:.1f}", help="-15 % G tolerance")
    col3.metric("Upper bound", f"{target.max:.1f}", help="+15 % G tolerance")
