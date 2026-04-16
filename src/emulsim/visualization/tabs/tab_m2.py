"""M2 Functionalization tab — extracted from app.py for UI restructure.

v6.0: Renders the complete M2 Surface Functionalization tab including inputs,
run button, and results display. All widget keys preserved exactly (m2_* prefix).
Widget key pattern includes step_type AND reagent_key (Codex R3-F4 fix).
"""

from __future__ import annotations

import streamlit as st

from emulsim.properties.database import PropertyDatabase
from emulsim.trust import assess_trust


def render_tab_m2(tab_container, _smgr) -> None:
    """Render the M2 Functionalization tab inside the given Streamlit container.

    Args:
        tab_container: Streamlit tab container from st.tabs().
        _smgr: SessionStateManager instance.
    """
    with tab_container:
        st.header("Module 2: Surface Functionalization")

        # ── Upstream M1 Status Banner ────────────────────────────────────
        if "result" not in st.session_state:
            st.warning("\u26a0\ufe0f Module 1 has not been run yet. Run M1 first to provide upstream data.")
        else:
            _m1r = st.session_state["result"]
            _m1e = _m1r.emulsification
            _m1g = _m1r.gelation
            _m1m = _m1r.mechanical
            st.success(
                f"\u2705 M1 data available \u2014 d32={_m1e.d32*1e6:.2f} um | "
                f"porosity={_m1g.porosity:.1%} | G_DN={_m1m.G_DN/1000:.1f} kPa"
            )

        # ── M2 Inputs ───────────────────────────────────────────────────
        from emulsim.module2_functionalization.reagent_profiles import REAGENT_PROFILES as _REAGENT_PROFILES
        from emulsim.module2_functionalization import ModificationStep, ModificationStepType, ACSSiteType

        m2_steps = []
        n_m2_steps = st.number_input("Modification Steps", 0, 3, 1, key="m2_n_steps")
        for i in range(int(n_m2_steps)):
            with st.expander(f"Step {i + 1}", expanded=(i == 0)):
                step_type = st.selectbox(
                    "Chemistry",
                    ["Secondary Crosslinking", "Hydroxyl Activation", "Ligand Coupling", "Protein Coupling", "Spacer Arm", "Metal Charging", "Protein Pretreatment", "Washing", "Quenching"],
                    key=f"m2_type_{i}",
                )
                if "Crosslinking" in step_type:
                    _reagent_options = {
                        "Genipin": "genipin_secondary",
                        "Glutaraldehyde": "glutaraldehyde_secondary",
                    }
                elif "Activation" in step_type:
                    _reagent_options = {
                        "ECH (Epichlorohydrin)": "ech_activation",
                        "DVS (Divinyl Sulfone)": "dvs_activation",
                        "BDGE (18 A Long-Arm Epoxide)": "bdge_activation",
                    }
                elif step_type == "Ligand Coupling":
                    _reagent_options = {
                        "DEAE (Weak Anion Exchange)": "deae_coupling",
                        "Q (Strong Anion Exchange)": "q_coupling",
                        "CM (Weak Cation Exchange)": "cm_coupling",
                        "Sulfopropyl (Strong Cation)": "sp_coupling",
                        "IDA (IMAC Chelator)": "ida_coupling",
                        "NTA (IMAC His-tag)": "nta_coupling",
                        "Phenyl (HIC)": "phenyl_coupling",
                        "Butyl (HIC Mild)": "butyl_coupling",
                        "Octyl (HIC Strong)": "octyl_coupling",
                        "Glutathione (GST-tag)": "glutathione_coupling",
                        "Heparin (Affinity+IEX)": "heparin_coupling",
                    }
                elif step_type == "Protein Coupling":
                    _reagent_options = {
                        "Protein A (IgG Affinity)": "protein_a_coupling",
                        "Protein G (IgG Broad Subclass)": "protein_g_coupling",
                        "Protein A/G Fusion (Broadest IgG)": "protein_ag_coupling",
                        "Streptavidin (Biotin-tag)": "streptavidin_coupling",
                        "Protein L (Kappa Light Chain)": "protein_l_coupling",
                        "Con A (Lectin, Mannose)": "con_a_coupling",
                        "WGA (Lectin, GlcNAc)": "wga_coupling",
                        "Protein A-Cys (Oriented)": "protein_a_cys_coupling",
                        "Protein G-Cys (Oriented)": "protein_g_cys_coupling",
                        "Generic Cys-Protein": "generic_cys_protein_coupling",
                    }
                elif step_type == "Spacer Arm":
                    _reagent_options = {
                        "DADPA (13 A, Amine)": "dadpa_spacer_arm",
                        "DAH (9 A, Amine)": "dah_spacer_arm",
                        "EDA (3 A, Amine)": "eda_spacer_arm",
                        "PEG-diamine 600 (35 A)": "peg600_spacer_arm",
                        "SM(PEG)2 (18 A, NHS-Maleimide)": "sm_peg2",
                        "SM(PEG)4 (32 A, NHS-Maleimide)": "sm_peg4",
                        "SM(PEG)12 (60 A, NHS-Maleimide)": "sm_peg12",
                        "SM(PEG)24 (95 A, NHS-Maleimide)": "sm_peg24",
                    }
                elif step_type == "Metal Charging":
                    _reagent_options = {
                        "Nickel (Ni2+ for IMAC)": "nickel_charging",
                        "Cobalt (Co2+ for IMAC)": "cobalt_charging",
                        "Copper (Cu2+ for IMAC)": "copper_charging",
                        "Zinc (Zn2+ for IMAC)": "zinc_charging",
                        "EDTA Stripping": "edta_stripping",
                    }
                elif step_type == "Protein Pretreatment":
                    _reagent_options = {
                        "TCEP Reduction (maleimide-compatible)": "tcep_reduction",
                        "DTT Reduction (requires removal)": "dtt_reduction",
                    }
                elif step_type == "Washing":
                    _reagent_options = {
                        "Wash Buffer (advisory)": "wash_buffer",
                    }
                else:  # Quenching
                    _reagent_options = {
                        "Ethanolamine (Epoxide Quench)": "ethanolamine_quench",
                        "2-Mercaptoethanol (VS Quench)": "mercaptoethanol_quench",
                        "NaBH4 (Aldehyde Quench)": "nabh4_quench",
                        "Acetic Anhydride (Amine Quench)": "acetic_anhydride_quench",
                    }
                # Key includes step_type to reset selection when chemistry changes
                _reagent_sel_key = f"m2_reagent_{i}_{step_type.replace(' ', '_').lower()}"
                _reagent_label = st.selectbox("Reagent", list(_reagent_options.keys()), key=_reagent_sel_key)
                # Defensive: if Streamlit returns a stale label from session state, fall back to first option
                _reagent_key = _reagent_options.get(_reagent_label, list(_reagent_options.values())[0])

                _profile = _REAGENT_PROFILES.get(_reagent_key)
                if _profile is None:
                    _reagent_key = list(_reagent_options.values())[0]
                    _profile = _REAGENT_PROFILES[_reagent_key]
                st.caption(f"k={_profile.k_forward:.1e} | E_a={_profile.E_a / 1000:.0f} kJ/mol")
                # Audit F16: display confidence, calibration, hazard
                _conf = getattr(_profile, 'confidence_tier', 'semi_quantitative')
                _cal = getattr(_profile, 'calibration_source', '')[:60]
                _haz = getattr(_profile, 'hazard_class', '') or 'low'
                st.caption(f"Confidence: {_conf} | Hazard: {_haz}")
                st.markdown(
                    f"[View mechanism & protocol](/reagent_detail"
                    f"?key={_reagent_key}&source=reagent_profiles"
                    f"&T={_profile.temperature_default}&t={_profile.time_default}"
                    f"&c=10.0&pH={_profile.ph_optimum})",
                )

                # Spacer selectbox for coupling/protein steps
                if step_type in ("Ligand Coupling", "Protein Coupling"):
                    _spacer_options = {
                        "None (Direct Coupling)": "",
                        "DADPA (13 A, EAH-standard)": "dadpa_spacer",
                        "AHA (10 A, NHS-standard)": "aha_spacer",
                        "DAH (9 A, AH-standard)": "dah_spacer",
                    }
                    _spacer_label = st.selectbox(
                        "Spacer Arm (optional)",
                        list(_spacer_options.keys()),
                        key=f"m2_spacer_{i}_{step_type.replace(' ', '_').lower()}",
                    )
                    _selected_spacer = _spacer_options.get(_spacer_label, "")

                # Keys include step_type AND reagent so switching either resets defaults (Codex R3-F4)
                _wk = f"{step_type.replace(' ', '_').lower()}_{_reagent_key}"
                _conc = st.number_input("Concentration (mM)", 0.5, 200.0, 10.0, key=f"m2_conc_{i}_{_wk}")
                _temp_C = st.slider("Temperature (C)", 4, 80,
                                     int(_profile.temperature_default - 273.15), key=f"m2_temp_{i}_{_wk}")
                _time_h = st.number_input("Time (h)", 0.25, 48.0,
                                           float(_profile.time_default / 3600), key=f"m2_time_{i}_{_wk}")
                _ph = st.slider("pH", 3.0, 14.0, float(_profile.ph_optimum), step=0.5, key=f"m2_ph_{i}_{_wk}")
                st.markdown(
                    f"[View mechanism & protocol (with your parameters)](/reagent_detail"
                    f"?key={_reagent_key}&source=reagent_profiles"
                    f"&T={_temp_C + 273.15}&t={_time_h * 3600}"
                    f"&c={_conc}&pH={_ph})",
                )

                if step_type == "Protein Coupling":
                    st.caption("Ranking only \u2014 activity retention and steric limits are illustrative defaults.")

                _step_type_map = {
                    "Secondary Crosslinking": (ModificationStepType.SECONDARY_CROSSLINKING, ACSSiteType.AMINE_PRIMARY),
                    "Hydroxyl Activation": (ModificationStepType.ACTIVATION, ACSSiteType.HYDROXYL),
                    "Ligand Coupling": (ModificationStepType.LIGAND_COUPLING, None),
                    "Protein Coupling": (ModificationStepType.PROTEIN_COUPLING, None),
                    "Spacer Arm": (ModificationStepType.SPACER_ARM, None),
                    "Metal Charging": (ModificationStepType.METAL_CHARGING, None),
                    "Protein Pretreatment": (ModificationStepType.PROTEIN_PRETREATMENT, None),
                    "Washing": (ModificationStepType.WASHING, None),
                    "Quenching": (ModificationStepType.QUENCHING, None),
                }
                _step_type_enum, _default_target = _step_type_map[step_type]
                if _default_target is not None:
                    _target_acs = _default_target
                else:
                    _target_acs = _REAGENT_PROFILES[_reagent_key].target_acs
                # Set product_acs for activation steps (needed to create EPOXIDE/VS profiles)
                _product_acs = _REAGENT_PROFILES[_reagent_key].product_acs
                m2_steps.append(ModificationStep(
                    step_type=_step_type_enum,
                    reagent_key=_reagent_key,
                    target_acs=_target_acs,
                    product_acs=_product_acs,
                    temperature=_temp_C + 273.15,
                    time=_time_h * 3600,
                    ph=_ph,
                    reagent_concentration=_conc,
                ))

        # ── Run M2 Button ────────────────────────────────────────────────
        st.divider()
        _m2_can_run = "result" in st.session_state and len(m2_steps) > 0
        m2_run_btn = st.button("\u25b6 Run M2: Functionalization", type="primary",
                                use_container_width=True, disabled=not _m2_can_run)

        if m2_run_btn:
            _smgr.invalidate_downstream(from_module=2)
            with st.spinner("Running Module 2: Functionalization..."):
                from emulsim.pipeline.orchestrator import export_for_module2
                from emulsim.module2_functionalization import ModificationOrchestrator

                _m1_result = st.session_state["result"]
                _m1_params = st.session_state["params"]
                _m1_overrides = st.session_state.get("m1_overrides", {})
                _m1_xl_key = st.session_state.get("m1_xl_key", "genipin")
                _m1_l2_mode = st.session_state.get("m1_l2_mode", "empirical")

                db_trust = PropertyDatabase()
                props_trust_m2 = db_trust.update_for_conditions(
                    _m1_params.formulation.T_oil, _m1_params.formulation.c_agarose,
                    _m1_params.formulation.c_chitosan, _m1_params.formulation.c_span80,
                )
                for _k, _v in _m1_overrides.items():
                    if hasattr(props_trust_m2, _k):
                        setattr(props_trust_m2, _k, _v)
                _trust_m2 = assess_trust(_m1_result, _m1_params, props_trust_m2,
                                         crosslinker_key=_m1_xl_key, l2_mode=_m1_l2_mode)
                try:
                    _contract = export_for_module2(
                        _m1_result, _trust_m2,
                        crosslinker_key=_m1_xl_key,
                        props=props_trust_m2,
                    )
                    _m2_orch = ModificationOrchestrator()
                    _m2_result = _m2_orch.run(_contract, m2_steps)
                    st.session_state["m2_result"] = _m2_result
                    st.session_state["m1_contract"] = _contract
                except Exception as _m2_ex:
                    st.error(f"Module 2 failed: {_m2_ex}")
                    st.session_state.pop("m2_result", None)
                    st.stop()  # Don't rerun — show the error
            st.rerun()

        # ── M2 Results Display ───────────────────────────────────────────
        if "m2_result" in st.session_state:
            from emulsim.visualization.plots_m2 import plot_acs_waterfall, plot_surface_area_comparison
            _m2 = st.session_state["m2_result"]

            st.divider()
            st.header("\U0001f4ca M2 Results")

            _mc1, _mc2, _mc3 = st.columns(3)
            _mc1.metric("Steps Executed", len(_m2.modification_history))
            _mc2.metric("G_DN Updated", f"{_m2.G_DN_updated / 1000:.1f} kPa")
            _mc3.metric("E* Updated", f"{_m2.E_star_updated / 1000:.1f} kPa")

            st.divider()

            st.markdown("**ACS Site Inventory (remaining after all steps)**")
            for _site_type, _profile in _m2.acs_profiles.items():
                st.write(f"**{_site_type.value}**: remaining = {_profile.remaining_sites:.2e} mol/particle")

            st.divider()

            st.markdown("**Modification History**")
            for _i, _mr in enumerate(_m2.modification_history):
                _tgt = _mr.step.target_acs
                _sites_before = _mr.acs_before.get(_tgt)
                _consumed_str = ""
                if _sites_before is not None:
                    _consumed = _mr.conversion * _sites_before.remaining_sites
                    _consumed_str = f" | sites consumed: {_consumed:.2e} mol/particle"
                st.write(
                    f"Step {_i + 1}: {_mr.step.reagent_key} \u2014 "
                    f"conversion {_mr.conversion:.1%}{_consumed_str}"
                )

            st.divider()

            _fig_wf = plot_acs_waterfall(_m2.modification_history)
            if _fig_wf is not None:
                st.plotly_chart(_fig_wf, use_container_width=True)
            _fig_sa = plot_surface_area_comparison(_m2.surface_model)
            if _fig_sa is not None:
                st.plotly_chart(_fig_sa, use_container_width=True)

            st.caption(
                "[!] semi_quantitative \u2014 ACS inventory uses simplified site-density model. "
                "Default rate constants are illustrative \u2014 user calibration required."
            )
