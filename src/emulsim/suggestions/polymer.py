"""'reduce_polymer' suggestion module."""

from __future__ import annotations

from ..properties.crosslink_derivation import alpha_for_target_G_DN
from .generators import (
    DEVIATION_TRIGGER,
    G_deviation,
    polymer_down_text,
)
from .types import Suggestion, SuggestionContext, TargetRange


def generate(ctx: SuggestionContext) -> Suggestion | None:
    if G_deviation(ctx) <= DEVIATION_TRIGGER:
        return None
    if ctx.G_DN_actual <= ctx.target_G:
        return None  # G is below target -> crosslinker is the right fix
    return Suggestion(
        key="reduce_polymer",
        display_text=polymer_down_text(ctx),
        severity="warning",
        context=ctx,
    )


def derive_target(ctx: SuggestionContext) -> TargetRange:
    result = alpha_for_target_G_DN(
        G_DN_current=ctx.G_DN_actual,
        target_G_DN=ctx.target_G,
        c_agarose_current=ctx.c_agarose,
        c_chitosan_current=ctx.c_chitosan,
    )
    return TargetRange(
        nominal=result.alpha_nominal,
        min=result.alpha_min,
        max=result.alpha_max,
        unit="× (multiplier)",
        limited_by=result.limited_by,
        confidence_tier=result.confidence_tier,
        assumptions=result.assumptions,
        is_qualitative_only=False,
        qualitative_reason="",
    )


def render_derivation(ctx: SuggestionContext, target: TargetRange) -> None:
    import streamlit as st

    from ..properties.crosslink_derivation import alpha_for_target_G_DN

    full = alpha_for_target_G_DN(
        G_DN_current=ctx.G_DN_actual,
        target_G_DN=ctx.target_G,
        c_agarose_current=ctx.c_agarose,
        c_chitosan_current=ctx.c_chitosan,
    )

    st.markdown(
        "The polymer-reduction suggestion comes from inverting the "
        "semi-dilute modulus-vs-concentration scaling for hydrogel networks."
    )

    st.subheader("Step 1 — Semi-dilute scaling: G ∝ c²")
    st.latex(r"\dfrac{G_{\mathrm{DN}}^{\mathrm{target}}}{G_{\mathrm{DN}}^{\mathrm{actual}}} = \alpha^{2}")
    st.markdown(
        f"- Actual G_DN = **{ctx.G_DN_actual/1000:.1f} kPa**\n"
        f"- Target G_DN = **{ctx.target_G/1000:.1f} kPa**\n"
        f"- Ratio target/actual = **{ctx.target_G/ctx.G_DN_actual:.3f}**"
    )

    st.subheader("Step 2 — Solve for the scaling factor α")
    st.latex(r"\alpha = \sqrt{\dfrac{G_{\mathrm{DN}}^{\mathrm{target}}}{G_{\mathrm{DN}}^{\mathrm{actual}}}}")
    st.markdown(
        f"- α (nominal) = **{target.nominal:.3f}** (apply uniformly to both polymers)"
    )

    st.subheader("Step 3 — Apply α to current polymer concentrations")
    st.markdown(
        f"- Current c_agarose = **{ctx.c_agarose:.1f} kg/m³**  →  "
        f"new c_agarose = **{full.new_c_agarose:.1f} kg/m³**\n"
        f"- Current c_chitosan = **{ctx.c_chitosan:.1f} kg/m³**  →  "
        f"new c_chitosan = **{full.new_c_chitosan:.1f} kg/m³**"
    )

    st.subheader("Step 4 — Scaling-factor band (±15 % G tolerance)")
    col1, col2, col3 = st.columns(3)
    col1.metric("Nominal α", f"{target.nominal:.3f}")
    col2.metric("Lower bound α", f"{target.min:.3f}")
    col3.metric("Upper bound α", f"{target.max:.3f}")

    if full.limited_by != "ok":
        st.warning(f"Feasibility constraint: **{full.limited_by}**")
