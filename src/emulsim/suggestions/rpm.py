"""'increase_rpm' / 'decrease_rpm' suggestion module.

Single module, two branches: if d32 is too large -> increase_rpm; if d32
is too small -> decrease_rpm. Both share the same derivation page content
(the Weber-number inversion) but with different surface text.
"""

from __future__ import annotations

from ..properties.emulsification_derivation import (
    rpm_for_target_d32,
)
from .generators import (
    DEVIATION_TRIGGER,
    d_deviation,
    rpm_down_text,
    rpm_up_text,
)
from .types import Suggestion, SuggestionContext, TargetRange


def generate(ctx: SuggestionContext) -> Suggestion | None:
    if d_deviation(ctx) <= DEVIATION_TRIGGER:
        return None
    direction_up = ctx.d32_actual > ctx.target_d32
    if direction_up:
        return Suggestion(
            key="increase_rpm",
            display_text=rpm_up_text(ctx),
            severity="warning",
            context=ctx,
            extras={"direction": "increase"},
        )
    return Suggestion(
        key="decrease_rpm",
        display_text=rpm_down_text(ctx),
        severity="warning",
        context=ctx,
        extras={"direction": "decrease"},
    )


def _build_target(ctx: SuggestionContext) -> TargetRange:
    # Interfacial tension: pull from SurfactantProfile if we had the key;
    # fall back to 0.03 N/m (span-80 at 2 %, typical).
    sigma = 0.03
    result = rpm_for_target_d32(
        target_d32=ctx.target_d32,
        D=ctx.impeller_D,
        rho_c=ctx.rho_oil,
        mu_c=ctx.mu_oil,
        sigma=sigma,
        phi_d=ctx.phi_d,
    )
    return TargetRange(
        nominal=result.nominal,
        min=result.min,
        max=result.max,
        unit="rpm",
        limited_by=result.limited_by,
        confidence_tier=result.confidence_tier,
        assumptions=result.assumptions,
        is_qualitative_only=False,
        qualitative_reason="",
    )


def derive_target(ctx: SuggestionContext) -> TargetRange:
    return _build_target(ctx)


def render_derivation(ctx: SuggestionContext, target: TargetRange) -> None:
    import streamlit as st

    from ..properties.emulsification_derivation import (
        rpm_for_target_d32,
    )

    sigma = 0.03
    # Re-run so we can expose the intermediate numbers.
    full = rpm_for_target_d32(
        target_d32=ctx.target_d32,
        D=ctx.impeller_D,
        rho_c=ctx.rho_oil,
        mu_c=ctx.mu_oil,
        sigma=sigma,
        phi_d=ctx.phi_d,
    )

    st.markdown(
        "The RPM suggestion comes from inverting the Sprow (1967) droplet-size "
        "correlation for Kolmogorov-regime turbulent emulsification. The same "
        "chain is used whether we recommend an increase or a decrease — the "
        "direction is set by which side of the target your current d32 sits on."
    )

    st.subheader("Step 1 — Weber number (at your target d32)")
    st.latex(r"\dfrac{d_{32}}{D} = 0.06 \cdot \mathrm{We}^{-0.6} \cdot (1 + 9\,\phi_d)")
    st.latex(r"\mathrm{We} = \dfrac{\rho_c N^2 D^3}{\sigma}")
    st.markdown(
        f"- Impeller diameter D = **{ctx.impeller_D*1000:.0f} mm**\n"
        f"- Continuous-phase density ρ_c = **{ctx.rho_oil:.0f} kg/m³**\n"
        f"- Interfacial tension σ (assumed) = **{sigma*1000:.0f} mN/m**\n"
        f"- Dispersed-phase fraction φ_d = **{ctx.phi_d:.3f}**\n"
        f"- Target d32 = **{ctx.target_d32*1e6:.1f} µm**\n"
        f"- Actual d32 = **{ctx.d32_actual*1e6:.1f} µm**"
    )

    st.subheader("Step 2 — Solve for N (rev/s) from the target Weber")
    st.latex(r"N = \sqrt{\dfrac{\mathrm{We}_{\mathrm{req}} \cdot \sigma}{\rho_c\,D^3}}")
    st.markdown(
        f"- Current RPM = **{ctx.rpm:.0f}**\n"
        f"- Target RPM (nominal) = **{target.nominal:.0f}**\n"
    )

    st.subheader("Step 3 — Kolmogorov + laminar feasibility checks")
    st.latex(r"\eta_K = \left( \dfrac{\mu_c^3}{\rho_c \varepsilon} \right)^{1/4}")
    st.latex(r"\mathrm{Re} = \dfrac{\rho_c N D^2}{\mu_c}")
    st.markdown(
        f"- Kolmogorov microscale at nominal RPM = **{full.kolmogorov_floor*1e6:.2f} µm**\n"
        f"  (target d32 should be ≥ 3× this)\n"
        f"- Reynolds number at nominal RPM = **{full.reynolds:,.0f}**\n"
        f"  (must be > 10⁴ for the Kolmogorov regime correlation to hold)"
    )

    st.subheader("Step 4 — Target value + band")
    col1, col2, col3 = st.columns(3)
    col1.metric("Nominal RPM", f"{target.nominal:.0f}")
    col2.metric("Lower bound", f"{target.min:.0f}", help="Corresponds to d32 + 20%")
    col3.metric("Upper bound", f"{target.max:.0f}", help="Corresponds to d32 − 20%")
