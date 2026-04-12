"""M1 Fabrication tab — extracted from app.py for UI restructure.

v6.0: Renders the complete M1 Fabrication tab including inputs, run button,
results display, optimization assessment, trust assessment, and calibration.
All widget keys preserved exactly (m1_* prefix).
"""

from __future__ import annotations

import time
from pathlib import Path

import numpy as np
import streamlit as st

from emulsim.datatypes import (
    SimulationParameters, FormulationParameters,
    EmulsificationParameters, MixerGeometry, SolverSettings,
    VesselGeometry, StirrerGeometry, HeatingConfig, KernelConfig,
    ModelMode,
)
from emulsim.properties.database import PropertyDatabase
from emulsim.pipeline.orchestrator import PipelineOrchestrator
from emulsim.trust import assess_trust
from emulsim.reagent_library import SURFACTANTS, CROSSLINKERS
from emulsim.visualization.plots import (
    plot_droplet_size_distribution,
    plot_phase_field,
    plot_crosslinking_kinetics,
    plot_hertz_contact,
    plot_kav_curve,
    plot_results_dashboard,
    plot_modulus_comparison,
)
from emulsim.visualization.ui_validators import validate_m1_inputs as _validate_m1
from emulsim.literature_constants import ALL_CONSTANTS
from emulsim.level3_crosslinking.solver import (
    available_amine_concentration,
    recommended_crosslinker_concentration,
)


def render_tab_m1(
    tab_container,
    is_stirred: bool,
    model_mode_enum: ModelMode,
    _smgr,
    _const_input,
    _proto_sections: dict,
) -> None:
    """Render the M1 Fabrication tab inside the given Streamlit container.

    Args:
        tab_container: Streamlit tab container from st.tabs().
        is_stirred: Whether stirred vessel mode is selected.
        model_mode_enum: Scientific mode enum.
        _smgr: SessionStateManager instance.
        _const_input: Callable for rendering per-constant selector with protocol link.
        _proto_sections: Dict of calibration protocol sections.
    """
    with tab_container:
        st.header("Module 1: Fabrication Pipeline (L1→L2→L3→L4)")

        # ── M1 INPUT SECTION ─────────────────────────────────────────────────
        m1_col_left, m1_col_right = st.columns(2)

        # ── Left column: Emulsification + Formulation ──
        with m1_col_left:
            st.subheader("Emulsification (L1)")

            if is_stirred:
                vessel_choice = st.selectbox(
                    "Vessel", ["Glass Beaker (100 mm)", "Jacketed Vessel (92 mm)"],
                    help="Beaker: flat-plate heating. Jacketed: hot-water circulation.",
                    key="m1_vessel",
                )
                stirrer_choice = st.selectbox(
                    "Stirrer", ["Stirrer A - Pitched Blade (59 mm)", "Stirrer B - Rotor-Stator (32 mm)"],
                    key="m1_stirrer",
                )
                is_stirrer_A = "Pitched" in stirrer_choice

                if is_stirrer_A:
                    rpm = st.slider("Stirrer Speed (RPM)", 800, 2000, 1300, step=50, key="m1_rpm")
                else:
                    rpm = st.slider("Stirrer Speed (RPM)", 800, 9000, 1800, step=100, key="m1_rpm_rs")

                t_emul = st.number_input("Emulsification Time (min)", 1, 60, 10, key="m1_t_emul")

                st.caption("Working liquid volumes")
                v_oil_mL = st.slider("Oil + Span-80 (mL)", 100, 500, 300, step=10, key="m1_v_oil")
                v_poly_mL = st.slider("Polysaccharide solution (mL)", 50, 300, 200, step=10, key="m1_v_poly")
                total_mL = v_oil_mL + v_poly_mL
                phi_d = v_poly_mL / total_mL
                st.caption(f"Total: {total_mL} mL | phi_d = {phi_d:.2f}")

                if "Beaker" in vessel_choice:
                    heating_choice = "Flat Plate"
                    st.caption("Heating: flat-plate (150C -> 80C oil)")
                else:
                    heating_choice = "Hot Water Jacket"
                    st.caption("Heating: jacket (85C circulating water)")
            else:
                rpm = st.slider("Rotor Speed (RPM)", 3000, 25000, 10000, step=500, key="m1_rpm_leg")
                t_emul = st.number_input("Emulsification Time (min)", 1, 60, 10, key="m1_t_emul_leg")
                phi_d = st.slider("Dispersed Phase Fraction (phi_d)", 0.01, 0.30, 0.05, step=0.01, key="m1_phi_d")
                v_oil_mL = 300
                v_poly_mL = 200
                total_mL = 500

            if model_mode_enum != ModelMode.EMPIRICAL_ENGINEERING:
                with st.expander("Advanced PBE Settings",
                                  expanded=(model_mode_enum == ModelMode.MECHANISTIC_RESEARCH)):
                    l1_t_max = st.slider("Max emulsification time (s)", 60, 1800, 600, step=60,
                                          key="m1_t_max",
                                          help="Absolute ceiling for adaptive time extensions")
                    l1_conv_tol = st.slider("Convergence tolerance", 0.005, 0.10, 0.01, step=0.005,
                                             format="%.3f", key="m1_conv_tol",
                                             help="Relative d32 variation threshold for steady state")
                    l1_max_ext = st.number_input("Max extensions", 0, 5, 2, key="m1_max_ext",
                                                  help="Number of half-interval extensions if PBE not converged")
            else:
                l1_t_max = 600.0
                l1_conv_tol = 0.01
                l1_max_ext = 2

            st.subheader("Formulation")

            _surf_keys = list(SURFACTANTS.keys())
            _surf_names = [SURFACTANTS[k].name for k in _surf_keys]
            _surf_sel_name = st.selectbox(
                "Surfactant", _surf_names,
                index=_surf_keys.index("span80"),
                help="Select surfactant for W/O emulsification",
                key="m1_surfactant",
            )
            _surf_sel_key = _surf_keys[_surf_names.index(_surf_sel_name)]
            surf = SURFACTANTS[_surf_sel_key]
            st.caption(f"HLB={surf.hlb} | MW={surf.mw} g/mol | {surf.notes[:60]}")

            _fc1, _fc2 = st.columns(2)
            with _fc1:
                c_agarose_pct = st.number_input("Agarose (% w/v)", 1.0, 10.0, 4.2, step=0.1, key="m1_c_agar")
            with _fc2:
                c_chitosan_pct = st.number_input("Chitosan (% w/v)", 0.5, 5.0, 1.8, step=0.1, key="m1_c_chit")

            if is_stirred:
                c_span80_vol_pct = st.slider("Span-80 in oil (% v/v)", 0.2, 5.0, 1.5, step=0.1, key="m1_span_vv")
                c_span80_pct = c_span80_vol_pct * 986.0 / 1000.0
                T_oil_C = st.slider("Paraffin Oil Temperature (C)", 65, 110, 80, step=1, key="m1_T_oil",
                                     help="Oil temperature at start of emulsification.")
            else:
                c_span80_pct = st.slider("Surfactant (% w/v)", 0.5, 5.0, 2.0, step=0.1, key="m1_span_wv")
                c_span80_vol_pct = 1.5
                T_oil_C = st.slider("Oil Temperature (C)", 60, 95, 90, key="m1_T_oil_leg")

        # ── Right column: Gelation + Crosslinking + Material Constants + Targets ──
        with m1_col_right:
            st.subheader("Cooling & Gelation (L2)")
            if is_stirred:
                cooling_rate_Cmin = st.slider("Cooling Rate (C/min)", 0.1, 15.0, 0.67, step=0.1,
                                               key="m1_cool_rate",
                                               help="Natural cooling: ~0.67 C/min for 500 mL beaker")
            else:
                cooling_rate_Cmin = st.slider("Cooling Rate (C/min)", 1.0, 20.0, 10.0, step=0.5,
                                               key="m1_cool_rate_leg")
            l2_model = st.radio(
                "Pore Structure Model",
                ["Empirical (fast, calibrated)", "Cahn-Hilliard 2D (mechanistic)"],
                index=0, key="m1_l2_model",
                help="Empirical: literature-calibrated power law (~1 ms). "
                     "CH 2D: phase-field spinodal decomposition (slower, requires parameter tuning).",
            )
            l2_mode = "empirical" if "Empirical" in l2_model else "ch_2d"
            grid_size = 64
            if l2_mode == "ch_2d":
                grid_size = st.select_slider("Phase-Field Grid", [32, 64, 128], value=64, key="m1_grid")

            st.subheader("Crosslinking (L3)")

            _xl_keys = list(CROSSLINKERS.keys())
            _xl_names = [CROSSLINKERS[k].name for k in _xl_keys]
            _xl_sel_name = st.selectbox(
                "Crosslinker", _xl_names,
                index=_xl_keys.index("genipin"),
                help="Select crosslinking agent",
                key="m1_crosslinker",
            )
            _xl_sel_key = _xl_keys[_xl_names.index(_xl_sel_name)]
            xl = CROSSLINKERS[_xl_sel_key]
            st.caption(f"{xl.mechanism} | k\u2080={xl.k_xlink_0:.1e} | Score: {xl.suitability}/10")

            c_genipin_mM = st.slider("Crosslinker Concentration (mM)", 0.5, 500.0, 10.0, step=0.5,
                                      key="m1_c_xl")

            _DDA = st.slider("DDA (degree of deacetylation)", 0.50, 0.99, 0.85, step=0.01,
                              key="m1_DDA",
                              help="Chitosan degree of deacetylation. Affects NH2 density and crosslinking capacity.")
            _M_GlcN = 161.16
            _c_chit_kg = c_chitosan_pct * 10.0
            _NH2 = available_amine_concentration(_c_chit_kg, _DDA, _M_GlcN)
            if _NH2 > 0:
                _ratio = c_genipin_mM / _NH2
                _c_rec = recommended_crosslinker_concentration(_c_chit_kg, _DDA, _M_GlcN, target_p=0.20)
                if _ratio >= 0.10:
                    st.success(f"Crosslinker/NH\u2082 = {_ratio:.3f} \u2014 sufficient for p \u2265 0.20")
                elif _ratio >= 0.05:
                    st.warning(
                        f"Crosslinker/NH\u2082 = {_ratio:.3f} \u2014 may be limiting. "
                        f"Recommend \u2265 {_c_rec:.1f} mM for target p = 0.20"
                    )
                else:
                    st.error(
                        f"Crosslinker/NH\u2082 = {_ratio:.3f} \u2014 severely limiting! "
                        f"Increase to at least {_c_rec:.1f} mM for target p = 0.20"
                    )

            _T_xlink_default = int(xl.T_crosslink_default - 273.15)
            _t_xlink_default_h = max(1, int(xl.t_crosslink_default / 3600))
            T_xlink_C = st.slider("Crosslinking Temperature (\u00b0C)", 0, 120,
                                    min(max(_T_xlink_default, 0), 120), key="m1_T_xl")
            t_xlink_h = st.slider("Crosslinking Time (hours)", 1, 48,
                                    min(_t_xlink_default_h, 48), key="m1_t_xl")

            if xl.kinetics_model == 'uv_dose':
                uv_intensity = st.slider("UV Intensity (mW/cm\u00b2)", 1.0, 100.0, 20.0, step=1.0,
                                          key="m1_uv")
            else:
                uv_intensity = 0.0

            st.subheader("Optimization Targets")
            if is_stirred:
                target_d_mode = st.number_input("Target d_mode (um)", 10.0, 500.0, 100.0, step=10.0,
                                                 key="m1_tgt_d", help="Modal diameter of microspheres")
                target_d32 = target_d_mode
                target_pore = st.number_input("Target Pore Size (nm)", 10, 500, 100, step=10, key="m1_tgt_pore")
            else:
                target_d32 = st.number_input("Target d32 (um)", 0.5, 50.0, 2.0, step=0.5, key="m1_tgt_d32")
                target_d_mode = target_d32
                target_pore = st.number_input("Target Pore Size (nm)", 10, 500, 80, step=10, key="m1_tgt_pore_leg")
            target_G = st.number_input("Target G_DN (kPa)", 1.0, 500.0, 10.0, step=1.0, key="m1_tgt_G")

            # Material Constants
            with st.expander("Material Constants", expanded=False):
                st.caption("For each constant: use literature or enter your own.")
                _mc = st  # container for material constants

                _kl = ALL_CONSTANTS["K_L"]
                custom_K_L = _const_input(
                    _mc, "K_L (Span-80 adsorption)", "K_L",
                    _kl.value, "m\u00b3/mol", "Santini 2007", 0.01, 10.0, 0.05, "K_L",
                )
                _gi = ALL_CONSTANTS["Gamma_inf"]
                custom_Gamma_inf = _const_input(
                    _mc, "\u0393\u221e (max surface excess)", "Gamma_inf",
                    _gi.value * 1e6, "\u00d710\u207b\u2076 mol/m\u00b2", "Santini 2007", 1.0, 10.0, 0.1,
                    "Gamma_inf", "%.1f",
                )
                _ec = ALL_CONSTANTS["eta_intr_chit"]
                custom_eta_chit = _const_input(
                    _mc, "[\u03b7] chitosan", "eta_chit",
                    _ec.value, "mL/g", "Rinaudo 1993", 100.0, 2000.0, 50.0,
                    "eta_chit", "%.0f",
                )
                _c3 = ALL_CONSTANTS["breakage_C3"]
                custom_C3 = _const_input(
                    _mc, "C3 (viscous breakage)", "C3",
                    _c3.value, "-", "Alopaeus 2002", 0.0, 1.0, 0.05,
                    "C3", "%.2f",
                )
                _xl_eta = xl.eta_coupling_recommended
                st.caption(f"Per-chemistry \u03b7 for {xl.name}: {_xl_eta:+.2f}")
                custom_eta_coup = _const_input(
                    _mc, "\u03b7 coupling (IPN)", "eta_coup",
                    _xl_eta, "-", f"{xl.name} library", -0.5, 0.5, 0.05,
                    "eta_coup", "%.2f",
                )
                eta_coup_source = st.session_state.get("src_eta_coup", "Literature")

        # ── M1 Validation ────────────────────────────────────────────────────
        _m1_val = _validate_m1(
            rpm=float(rpm),
            phi_d=phi_d,
            c_agarose=c_agarose_pct,
            c_chitosan=c_chitosan_pct,
            dda=_DDA,
            crosslinker_key=_xl_sel_key,
            crosslinker_conc=float(c_genipin_mM),
            T_crosslink=float(T_xlink_C),
            T_oil=float(T_oil_C),
        )
        if _m1_val.blockers:
            for _blk in _m1_val.blockers:
                st.error(f"BLOCKER: {_blk}")
        if _m1_val.warnings:
            for _wrn in _m1_val.warnings:
                st.warning(_wrn)

        # ── Build Parameters ─────────────────────────────────────────────────
        if is_stirred:
            if "Beaker" in vessel_choice:
                _vessel = VesselGeometry.glass_beaker(working_volume=total_mL * 1e-6)
            else:
                _vessel = VesselGeometry.jacketed_vessel(working_volume=total_mL * 1e-6)

            if is_stirrer_A:
                _stirrer = StirrerGeometry.pitched_blade_A()
            else:
                _stirrer = StirrerGeometry.rotor_stator_B()

            if "Beaker" in vessel_choice:
                _heating = HeatingConfig.flat_plate()
            else:
                _heating = HeatingConfig.hot_water_jacket()
            _heating.T_initial = T_oil_C + 273.15
            _kernels = KernelConfig.for_stirrer_type(_stirrer.stirrer_type)

            params = SimulationParameters(
                model_mode=model_mode_enum,
                emulsification=EmulsificationParameters(
                    mode="stirred_vessel",
                    rpm=float(rpm),
                    t_emulsification=float(t_emul * 60),
                    vessel=_vessel,
                    stirrer=_stirrer,
                    heating=_heating,
                    kernels=_kernels,
                ),
                formulation=FormulationParameters(
                    c_agarose=c_agarose_pct * 10.0,
                    c_chitosan=c_chitosan_pct * 10.0,
                    c_span80=c_span80_pct * 10.0,
                    c_genipin=float(c_genipin_mM),
                    T_oil=T_oil_C + 273.15,
                    cooling_rate=cooling_rate_Cmin / 60.0,
                    T_crosslink=T_xlink_C + 273.15,
                    t_crosslink=float(t_xlink_h * 3600),
                    phi_d=phi_d,
                    c_span80_vol_pct=c_span80_vol_pct,
                    v_oil_span80_mL=float(v_oil_mL),
                    v_polysaccharide_mL=float(v_poly_mL),
                ),
                solver=SolverSettings(
                    l1_n_bins=40, l1_d_min=1e-6, l1_d_max=2000e-6,
                    l2_n_grid=grid_size,
                    l1_t_max=float(l1_t_max),
                    l1_conv_tol=float(l1_conv_tol),
                    l1_max_extensions=int(l1_max_ext),
                ),
            )
        else:
            params = SimulationParameters(
                model_mode=model_mode_enum,
                emulsification=EmulsificationParameters(
                    mode="rotor_stator_legacy",
                    rpm=float(rpm),
                    t_emulsification=float(t_emul * 60),
                    mixer=MixerGeometry(),
                ),
                formulation=FormulationParameters(
                    c_agarose=c_agarose_pct * 10.0,
                    c_chitosan=c_chitosan_pct * 10.0,
                    c_span80=c_span80_pct * 10.0,
                    c_genipin=float(c_genipin_mM),
                    T_oil=T_oil_C + 273.15,
                    cooling_rate=cooling_rate_Cmin / 60.0,
                    T_crosslink=T_xlink_C + 273.15,
                    t_crosslink=float(t_xlink_h * 3600),
                    phi_d=phi_d,
                ),
                solver=SolverSettings(
                    l2_n_grid=grid_size,
                    l1_t_max=float(l1_t_max),
                    l1_conv_tol=float(l1_conv_tol),
                    l1_max_extensions=int(l1_max_ext),
                ),
            )

        # Material constants overrides
        _custom_props_overrides = {
            "breakage_C3": custom_C3,
            "eta_coupling": custom_eta_coup,
            "eta_is_custom": (eta_coup_source == "Custom"),
        }
        _custom_props_overrides['k_xlink_0'] = xl.k_xlink_0
        _custom_props_overrides['E_a_xlink'] = xl.E_a_xlink
        _custom_props_overrides['f_bridge'] = xl.f_bridge
        _R_gas = 8.314
        _T_ift = T_oil_C + 273.15
        _c_mol_ift = c_span80_pct * 10.0 / surf.mw * 1000.0
        _sigma_0_T = max(surf.sigma_0_paraffin - 1.0e-4 * (_T_ift - 293.15), 0.001)
        _use_gamma = custom_Gamma_inf * 1e-6
        _use_KL = custom_K_L
        _sigma_calc = max(
            _sigma_0_T - _R_gas * _T_ift * _use_gamma * np.log(1 + _use_KL * _c_mol_ift),
            1e-4,
        )
        _custom_props_overrides['sigma'] = _sigma_calc
        _custom_props_overrides['eta_intr_chit'] = custom_eta_chit

        # ── Run M1 Button ────────────────────────────────────────────────────
        st.divider()
        m1_run_btn = st.button("\u25b6 Run M1: Fabrication Pipeline", type="primary",
                                use_container_width=True, disabled=bool(_m1_val.blockers))

        if m1_run_btn:
            _smgr.invalidate_downstream(from_module=1)
            with st.spinner("Running L1\u2192L2\u2192L3\u2192L4 pipeline..."):
                t_start = time.time()
                db = PropertyDatabase()
                orch = PipelineOrchestrator(db=db)
                progress = st.progress(0, text="Level 1: Emulsification (PBE solver)...")

                try:
                    result = orch.run_single(params, l2_mode=l2_mode,
                                             props_overrides=_custom_props_overrides,
                                             crosslinker_key=_xl_sel_key,
                                             uv_intensity=uv_intensity)
                except Exception as ex:
                    st.error(f"Simulation failed: {ex}")
                    st.stop()
                elapsed = time.time() - t_start
                progress.progress(100, text=f"M1 complete in {elapsed:.1f}s")

            st.session_state["result"] = result
            st.session_state["elapsed"] = elapsed
            st.session_state["params"] = params
            st.session_state["targets"] = (target_d32, target_pore, target_G)
            st.session_state["m1_overrides"] = _custom_props_overrides
            st.session_state["m1_xl_key"] = _xl_sel_key
            st.session_state["m1_l2_mode"] = l2_mode

            db_trust = PropertyDatabase()
            props_trust = db_trust.update_for_conditions(
                params.formulation.T_oil, params.formulation.c_agarose,
                params.formulation.c_chitosan, params.formulation.c_span80,
            )
            for _k, _v in _custom_props_overrides.items():
                if hasattr(props_trust, _k):
                    setattr(props_trust, _k, _v)
            st.session_state["trust"] = assess_trust(result, params, props_trust,
                                                        crosslinker_key=_xl_sel_key,
                                                        l2_mode=l2_mode)
            st.rerun()

        # ── M1 Results Display ───────────────────────────────────────────────
        if "result" in st.session_state:
            result = st.session_state["result"]
            elapsed = st.session_state["elapsed"]
            target_d32, target_pore, target_G = st.session_state["targets"]

            e = result.emulsification
            g = result.gelation
            x = result.crosslinking
            m = result.mechanical

            st.header("\U0001f4ca M1 Results")

            # KPI cards
            col1, col2, col3, col4 = st.columns(4)
            _d_mode = getattr(e, 'd_mode', 0.0)

            if is_stirred:
                d_primary_dev = abs(_d_mode * 1e6 - target_d_mode) / target_d_mode * 100
            else:
                d_primary_dev = abs(e.d32 * 1e6 - target_d32) / target_d32 * 100
            pore_dev = abs(g.pore_size_mean * 1e9 - target_pore) / target_pore * 100
            G_dev = abs(m.G_DN / 1000 - target_G) / target_G * 100

            if is_stirred:
                col1.metric("d_mode", f"{_d_mode*1e6:.1f} um",
                            delta=f"{d_primary_dev:.0f}% from target", delta_color="inverse")
            else:
                col1.metric("d32", f"{e.d32*1e6:.2f} um",
                            delta=f"{d_primary_dev:.0f}% from target", delta_color="inverse")
            col2.metric("Pore Size", f"{g.pore_size_mean*1e9:.1f} nm",
                        delta=f"{pore_dev:.0f}% from target", delta_color="inverse")
            _hs_lo = getattr(m, 'G_DN_lower', 0.0)
            _hs_hi = getattr(m, 'G_DN_upper', 0.0)
            if _hs_lo > 0 and _hs_hi > 0:
                col3.metric("G_DN", f"{m.G_DN/1000:.1f} kPa",
                            delta=f"Ref: [{_hs_lo/1000:.1f}, {_hs_hi/1000:.1f}] kPa (single-phase)",
                            delta_color="off")
            else:
                col3.metric("G_DN", f"{m.G_DN/1000:.1f} kPa",
                            delta=f"{G_dev:.0f}% from target", delta_color="inverse")
            col4.metric("Pipeline Time", f"{elapsed:.1f}s")

            col5, col6, col7, col8 = st.columns(4)
            col5.metric("Span", f"{e.span:.2f}")
            col6.metric("Porosity", f"{g.porosity:.1%}")
            col7.metric("Crosslink %", f"{x.p_final:.1%}")
            col8.metric("E*", f"{m.E_star/1000:.1f} kPa")

            st.divider()

            # ── L1-L4 Sub-tabs ───────────────────────────────────────────────
            _sub_labels = ["\U0001f4c8 Dashboard", "\U0001fab8 L1: Emulsification", "\U0001f9ca L2: Gelation",
                           "\U0001f517 L3: Crosslinking", "\U0001f4aa L4: Mechanical"]
            sub1, sub2, sub3, sub4, sub5 = st.tabs(_sub_labels)

            with sub1:
                st.plotly_chart(plot_results_dashboard(result), use_container_width=True)

            with sub2:
                st.subheader("Level 1: Emulsification -- Droplet Size Distribution")
                st.plotly_chart(plot_droplet_size_distribution(e), use_container_width=True)
                c1, c2, c3 = st.columns(3)
                c1.write(f"**d10** = {e.d10*1e6:.2f} um")
                c1.write(f"**d32** = {e.d32*1e6:.2f} um")
                if _d_mode > 0:
                    c1.write(f"**d_mode** = {_d_mode*1e6:.1f} um")
                c2.write(f"**d50** = {e.d50*1e6:.2f} um")
                c2.write(f"**d90** = {e.d90*1e6:.2f} um")
                c3.write(f"**d43** = {e.d43*1e6:.2f} um")
                c3.write(f"**Span** = {e.span:.3f}")
                _conv_icon = "\u2705 Yes" if e.converged else "\u26a0\ufe0f No (still evolving)"
                st.write(f"**Converged:** {_conv_icon}")
                _t_conv = getattr(e, 't_converged', None)
                _n_ext = getattr(e, 'n_extensions', 0)
                if _t_conv is not None:
                    st.write(f"Converged at t = {_t_conv:.1f} s")
                if _n_ext > 0:
                    st.write(f"Adaptive extensions: {_n_ext}")
                if is_stirred:
                    st.caption(
                        f"Mode: stirred-vessel | "
                        f"Vessel: {vessel_choice.split('(')[0].strip()} | "
                        f"Stirrer: {stirrer_choice.split('-')[0].strip()} | "
                        f"phi_d = {phi_d:.2f}"
                    )

            with sub3:
                st.subheader("Level 2: Gelation \u2014 Pore Structure")
                st.plotly_chart(plot_phase_field(g), use_container_width=True)
                c1, c2, c3 = st.columns(3)
                c1.write(f"**Pore size** = {g.pore_size_mean*1e9:.1f} nm")
                c2.write(f"**Porosity** = {g.porosity:.3f}")
                c3.write(f"**Gelation \u03b1** = {g.alpha_final:.4f}")
                if g.phi_field.ndim == 2:
                    st.write(f"Grid: {g.phi_field.shape[0]}\u00d7{g.phi_field.shape[1]}, "
                             f"spacing = {g.grid_spacing*1e9:.1f} nm, "
                             f"\u03c6 range: [{g.phi_field.min():.4f}, {g.phi_field.max():.4f}]")

            with sub4:
                st.subheader("Level 3: Crosslinking Kinetics")
                st.plotly_chart(plot_crosslinking_kinetics(x), use_container_width=True)
                c1, c2, c3 = st.columns(3)
                c1.write(f"**Crosslink fraction** = {x.p_final:.4f}")
                c1.write(f"**G_crosslinked** = {x.G_chitosan_final:.0f} Pa")
                c2.write(f"**Mesh size \u03be** = {x.xi_final*1e9:.1f} nm")
                c2.write(f"**M_c** = {x.Mc_final:.0f} g/mol")
                c3.write(f"**\u03bd_e** = {x.nu_e_final:.2e} /m\u00b3")

            with sub5:
                st.subheader("Level 4: Mechanical Properties")
                _model_label = getattr(m, 'model_used', 'phenomenological')
                st.caption(f"Model: {_model_label}")
                if _model_label == "flory_rehner_affine":
                    st.success("Flory-Rehner affine IPN model converged \u2014 mechanistic crosslink density used.")
                elif _model_label == "phenomenological":
                    if model_mode_enum == ModelMode.MECHANISTIC_RESEARCH:
                        st.info("Using phenomenological model. Switch to Mechanistic Research mode to enable Flory-Rehner affine IPN.")
                    else:
                        st.info("Phenomenological DN model (G1 + G2 + \u03b7\u221a(G1\u00b7G2)). Switch to Mechanistic Research for affine IPN.")

                left, right = st.columns(2)
                with left:
                    st.plotly_chart(plot_hertz_contact(m), use_container_width=True)
                with right:
                    st.plotly_chart(plot_kav_curve(m), use_container_width=True)
                st.plotly_chart(plot_modulus_comparison(m), use_container_width=True)

                c1, c2, c3 = st.columns(3)
                c1.write(f"**G_agarose** = {m.G_agarose:.0f} Pa ({m.G_agarose/1000:.1f} kPa)")
                c2.write(f"**G_crosslinked** = {m.G_chitosan:.0f} Pa ({m.G_chitosan/1000:.1f} kPa)")
                c3.write(f"**G_DN** = {m.G_DN:.0f} Pa ({m.G_DN/1000:.1f} kPa)")
                if _hs_lo > 0 and _hs_hi > 0:
                    st.write(f"**Single-phase reference (HS composite bounds, not IPN bounds):** [{_hs_lo/1000:.1f}, {_hs_hi/1000:.1f}] kPa")
                st.caption(f"Model: {_model_label}")

                if model_mode_enum == ModelMode.MECHANISTIC_RESEARCH:
                    st.write("**Model Comparison:**")
                    from emulsim.level4_mechanical.solver import double_network_modulus as _dnm
                    _eta_comp = getattr(x.network_metadata, 'eta_coupling_recommended', -0.15) if x.network_metadata else -0.15
                    _G_pheno = _dnm(m.G_agarose, m.G_chitosan, _eta_comp)
                    st.write(
                        f"Phenomenological: {_G_pheno/1000:.1f} kPa | "
                        f"{'Affine IPN' if _model_label == 'flory_rehner_affine' else _model_label}: {m.G_DN/1000:.1f} kPa"
                    )

            # ── Optimization Assessment ──────────────────────────────────────
            st.divider()
            st.header("\U0001f3af Optimization Assessment")

            if is_stirred:
                d_obj_val = _d_mode * 1e6
                d_obj_target = target_d_mode
                d_obj_label = "f_1 (d_mode deviation)"
                d_obj_help = f"|d_mode - {target_d_mode} um| / {target_d_mode} um"
            else:
                d_obj_val = e.d32 * 1e6
                d_obj_target = target_d32
                d_obj_label = "f_1 (d32 deviation)"
                d_obj_help = f"|d32 - {target_d32} um| / {target_d32} um"

            d_dev_obj = abs(d_obj_val - d_obj_target) / d_obj_target
            pore_dev_obj = abs(g.pore_size_mean * 1e9 - target_pore) / target_pore
            G_dev_obj = abs(np.log10(max(m.G_DN, 1)) - np.log10(target_G * 1000))
            obj = np.array([d_dev_obj, pore_dev_obj, G_dev_obj])

            st.write("**Objective Values** (lower = closer to target):")
            oc1, oc2, oc3 = st.columns(3)
            oc1.metric(d_obj_label, f"{obj[0]:.3f}", help=d_obj_help)
            oc2.metric("f_2 (pore deviation)", f"{obj[1]:.3f}",
                       help=f"|pore - {target_pore} nm| / {target_pore} nm")
            oc3.metric("f_3 (modulus deviation)", f"{obj[2]:.3f}",
                       help=f"|log10(G_DN) - log10({target_G*1000})|")

            overall = np.mean(obj)
            if overall < 0.3:
                st.success(f"\U0001f7e2 **Excellent match** (avg. deviation = {overall:.3f}). Parameters are near-optimal.")
            elif overall < 1.0:
                st.warning(f"\U0001f7e1 **Moderate match** (avg. deviation = {overall:.3f}). Consider optimization.")
            else:
                st.error(f"\U0001f534 **Poor match** (avg. deviation = {overall:.3f}). Significant parameter adjustment needed.")

            recs = []
            if obj[0] > 0.5:
                if d_obj_val > d_obj_target:
                    recs.append(f"**Increase RPM** (currently {rpm}) -- droplet size is {d_obj_val:.1f} um vs target {d_obj_target:.0f} um. "
                               f"Try higher RPM or increase Span-80 concentration.")
                else:
                    recs.append(f"**Decrease RPM** -- droplet size is {d_obj_val:.1f} um, smaller than target {d_obj_target:.0f} um.")
            if obj[1] > 0.5:
                recs.append(f"**Adjust cooling rate** \u2014 pore size ({g.pore_size_mean*1e9:.0f} nm) deviates from "
                            f"target ({target_pore} nm). Slower cooling \u2192 larger pores; faster \u2192 finer.")
            if obj[2] > 0.5:
                if m.G_DN / 1000 < target_G:
                    recs.append(f"**Increase crosslinker concentration** or crosslinking time \u2014 G_DN ({m.G_DN/1000:.1f} kPa) below "
                               f"target ({target_G} kPa). Current crosslinker is stoichiometry-limited at p={x.p_final:.1%}.")
                else:
                    recs.append(f"**Reduce polymer concentration** \u2014 G_DN ({m.G_DN/1000:.1f} kPa) exceeds target ({target_G} kPa).")
            if not recs:
                st.write("All objectives are within acceptable range. No adjustments needed.")
            else:
                for i, rec in enumerate(recs, 1):
                    st.write(f"{i}. {rec}")

            # ── Trust Assessment ─────────────────────────────────────────────
            if "trust" in st.session_state:
                st.divider()
                st.header("Trust Assessment")
                trust = st.session_state["trust"]
                _mode_desc = {
                    "empirical_engineering": "Empirical Engineering \u2014 trust warnings relaxed for screening",
                    "hybrid_coupled": "Hybrid Coupled \u2014 phenomenological models with trust warnings",
                    "mechanistic_research": "Mechanistic Research \u2014 strictest gates, Flory-Rehner when available",
                }
                st.caption(f"Mode: {_mode_desc.get(model_mode_enum.value, model_mode_enum.value)}")
                if trust.level == "TRUSTWORTHY":
                    st.success(f"**{trust.level}** -- All checks passed.")
                elif trust.level == "CAUTION":
                    st.warning(f"**{trust.level}** -- Some conditions are outside ideal range.")
                else:
                    st.error(f"**{trust.level}** -- Results should not be used for decisions.")
                if trust.blockers:
                    st.subheader("Blockers")
                    for b in trust.blockers:
                        st.error(b)
                if trust.warnings:
                    st.subheader("Warnings")
                    for w in trust.warnings:
                        st.warning(w)

            # ── Calibration Protocol ─────────────────────────────────────────
            st.divider()
            st.header("\U0001f4cb Calibration & Validation")
            cal_path = Path(__file__).resolve().parents[4] / "docs" / "04_calibration_protocol.md"
            if cal_path.exists():
                with st.expander("View Calibration Wet-Lab Protocol"):
                    st.markdown(cal_path.read_text(encoding="utf-8"))
            st.info(
                "The simulation uses literature-estimated constants that should be calibrated "
                "against your specific materials. See **docs/04_calibration_protocol.md** for "
                "a 5-study, ~30-preparation wet-lab protocol covering: interfacial tension "
                "(K_L, \u0393\u221e), chitosan viscosity (\u03b7_intr), breakage dynamics (C3), "
                "pore structure (empirical coefficients), and IPN mechanics (\u03b7_coupling)."
            )
