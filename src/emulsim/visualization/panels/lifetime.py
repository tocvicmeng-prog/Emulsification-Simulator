"""Lifetime projection display panel for EmulSim UI.

v6.0-rc: Empirical resin lifetime projection (audit F6).
All outputs labeled "empirical" — not predictive without calibration.
"""

from __future__ import annotations

import numpy as np
import streamlit as st

from emulsim.lifetime import LifetimeProjection, project_lifetime


def render_lifetime_panel() -> LifetimeProjection | None:
    """Render the lifetime projection panel.

    Returns:
        LifetimeProjection if computed, else None.
    """
    with st.expander("\u23f3 Resin Lifetime Projection", expanded=False):
        st.caption(
            "Empirical first-order deactivation model: "
            "capacity(n) = initial \u00d7 exp(\u2212k \u00d7 n). "
            "Requires user-calibrated k from cycle studies."
        )

        _c1, _c2 = st.columns(2)
        with _c1:
            initial_cap = st.number_input(
                "Initial capacity (mol/m\u00b3)", 1.0, 500.0, 100.0,
                step=10.0, key="lt_init_cap",
                help="Fresh resin DBC or static binding capacity",
            )
            k_deact = st.number_input(
                "k_deactivation (1/cycle)", 0.0, 0.1, 0.005,
                step=0.001, format="%.4f", key="lt_k_deact",
                help="First-order deactivation constant. "
                     "Protein A: 0.002\u20130.01, IMAC Ni: 0.005\u20130.02",
            )
        with _c2:
            n_cycles = st.number_input(
                "Cycles to project", 10, 1000, 100,
                step=10, key="lt_n_cycles",
            )
            cip_desc = st.text_input(
                "CIP conditions assumed",
                "standard NaOH CIP",
                key="lt_cip",
                help="Document CIP protocol assumed for this projection",
            )

        # ── Compute ──
        if st.button("Project Lifetime", key="lt_run"):
            proj = project_lifetime(
                initial_capacity=initial_cap,
                k_deactivation=k_deact,
                n_cycles=n_cycles,
                cip_description=cip_desc,
            )
            st.session_state["_lt_result"] = proj

        # ── Display ──
        if "_lt_result" in st.session_state:
            proj = st.session_state["_lt_result"]

            _mc1, _mc2, _mc3 = st.columns(3)
            _mc1.metric(
                "Cycles to 80%",
                f"{proj.cycles_to_80pct}" if proj.cycles_to_80pct < 999999 else "\u221e",
            )
            _mc2.metric(
                "Cycles to 50%",
                f"{proj.cycles_to_50pct}" if proj.cycles_to_50pct < 999999 else "\u221e",
            )
            _mc3.metric(
                f"Capacity @ cycle {proj.n_cycles_queried}",
                f"{proj.capacity_at_n:.1f} mol/m\u00b3",
            )

            # ── Decay Curve ──
            if k_deact > 0:
                cycles = np.arange(0, n_cycles + 1)
                capacity = initial_cap * np.exp(-k_deact * cycles)

                import plotly.graph_objects as go
                fig = go.Figure()
                fig.add_trace(go.Scatter(
                    x=cycles, y=capacity,
                    mode="lines", name="Projected capacity",
                    line=dict(color="#1f77b4", width=2),
                ))
                # 80% line
                fig.add_hline(
                    y=initial_cap * 0.8, line_dash="dash",
                    line_color="orange",
                    annotation_text="80% threshold",
                )
                # 50% line
                fig.add_hline(
                    y=initial_cap * 0.5, line_dash="dash",
                    line_color="red",
                    annotation_text="50% threshold",
                )
                fig.update_layout(
                    title="Empirical Lifetime Projection",
                    xaxis_title="Cycle Number",
                    yaxis_title="Capacity (mol/m\u00b3)",
                    height=350,
                )
                st.plotly_chart(fig, use_container_width=True)

            # ── Confidence Label ──
            st.warning(
                f"\u26a0\ufe0f **{proj.confidence.upper()}** projection \u2014 "
                f"not predictive without calibrated k from real cycle data. "
                f"Assumes: {proj.assumption_notes}"
            )

            return proj

    return None
