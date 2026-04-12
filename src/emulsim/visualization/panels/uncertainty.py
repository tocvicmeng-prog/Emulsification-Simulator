"""Uncertainty configuration panel for EmulSim UI.

v6.0-rc: Exposes M1UncertaintyContract CVs and Monte Carlo settings.
Audit F4: Two-tier uncertainty — measured (Tier 1) vs assumed screening (Tier 2).
"""

from __future__ import annotations

import streamlit as st

from emulsim.uncertainty_propagation import M1UncertaintyContract


def render_uncertainty_panel() -> M1UncertaintyContract:
    """Render the uncertainty configuration panel.

    Returns:
        M1UncertaintyContract with user-configured CVs.
    """
    with st.expander("\U0001f4ca Uncertainty (M1)", expanded=False):
        st.caption(
            "Configure coefficient of variation (CV) for M1 fabrication parameters. "
            "Monte Carlo propagation runs N perturbed simulations through M2."
        )

        # ── Tier Selection ──
        tier = st.radio(
            "Uncertainty Tier",
            ["Assumed (screening defaults)", "Measured (user-supplied CVs)"],
            index=0,
            key="unc_tier",
            help="Assumed: literature-typical CVs for initial screening. "
                 "Measured: enter CVs from your replicate experiments.",
        )
        is_measured = "Measured" in tier

        if not is_measured:
            # Load screening defaults
            defaults = M1UncertaintyContract.screening_defaults()
            st.info(
                f"Using screening defaults: "
                f"d50 CV={defaults.cv_bead_d50:.0%}, "
                f"porosity CV={defaults.cv_porosity:.0%}, "
                f"pore CV={defaults.cv_pore_size:.0%}, "
                f"NH2 CV={defaults.cv_nh2_bulk:.0%}, "
                f"OH CV={defaults.cv_oh_bulk:.0%}"
            )
            contract = defaults
        else:
            # User-editable CVs
            _c1, _c2 = st.columns(2)
            with _c1:
                cv_d50 = st.slider(
                    "CV bead d50", 0.0, 0.30, 0.10, step=0.01,
                    key="unc_cv_d50",
                    help="Coefficient of variation on median bead diameter",
                )
                cv_porosity = st.slider(
                    "CV porosity", 0.0, 0.20, 0.05, step=0.01,
                    key="unc_cv_por",
                )
                cv_pore = st.slider(
                    "CV pore size", 0.0, 0.30, 0.10, step=0.01,
                    key="unc_cv_pore",
                )
            with _c2:
                cv_nh2 = st.slider(
                    "CV NH2 bulk", 0.0, 0.20, 0.08, step=0.01,
                    key="unc_cv_nh2",
                    help="DDA variability affects amine site density",
                )
                cv_oh = st.slider(
                    "CV OH bulk", 0.0, 0.20, 0.05, step=0.01,
                    key="unc_cv_oh",
                )

            contract = M1UncertaintyContract(
                cv_bead_d50=cv_d50,
                cv_porosity=cv_porosity,
                cv_pore_size=cv_pore,
                cv_nh2_bulk=cv_nh2,
                cv_oh_bulk=cv_oh,
                tier="measured",
            )

        # ── Monte Carlo Settings ──
        n_samples = st.number_input(
            "Monte Carlo samples", 10, 1000, 100, step=10,
            key="unc_n_samples",
            help="Number of perturbed M1->M2 runs. Higher = tighter bounds, slower.",
        )
        st.session_state["unc_n_samples_val"] = int(n_samples)

        # ── Display current config ──
        if contract.is_deterministic():
            st.warning("All CVs are zero \u2014 deterministic mode (no uncertainty propagation).")
        else:
            _tier_label = "Measured (Tier 1)" if is_measured else "Assumed (Tier 2)"
            st.caption(f"Tier: {_tier_label} | Samples: {n_samples}")

    return contract
