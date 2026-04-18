"""M3 Performance tab — extracted from app.py for UI restructure.

v6.0: Renders the complete M3 Performance Characterization tab including
inputs, run button, and results display. All widget keys preserved exactly
(m3_* prefix). Chromatography and catalysis modes supported.
"""

from __future__ import annotations

import numpy as np
import streamlit as st


def render_tab_m3(tab_container) -> None:
    """Render the M3 Performance tab inside the given Streamlit container.

    Args:
        tab_container: Streamlit tab container from st.tabs().
    """
    with tab_container:
        st.header("Module 3: Performance Characterization")

        # ── Upstream M2 Status Banner ────────────────────────────────────
        if "m2_result" not in st.session_state:
            st.warning("\u26a0\ufe0f Module 2 has not been run yet. Run M1 then M2 first to provide upstream data.")
        else:
            _m2r_banner = st.session_state["m2_result"]
            st.success(
                f"\u2705 M2 data available \u2014 G_DN={_m2r_banner.G_DN_updated/1000:.1f} kPa | "
                f"E*={_m2r_banner.E_star_updated/1000:.1f} kPa | "
                f"Steps={len(_m2r_banner.modification_history)}"
            )

        # ── M3 Inputs ───────────────────────────────────────────────────
        m3_app_mode = st.radio("Application", ["Chromatography", "Catalysis"], key="m3_mode")

        st.subheader("Column/Reactor")
        _m3c1, _m3c2 = st.columns(2)
        with _m3c1:
            col_diam_mm = st.number_input("Column I.D. (mm)", 1.0, 50.0, 10.0, key="m3_col_d")
            bed_height_cm = st.number_input("Bed height (cm)", 1.0, 30.0, 10.0, key="m3_bed_h")
        with _m3c2:
            bed_porosity = st.slider("Bed porosity", 0.25, 0.50, 0.38, step=0.01, key="m3_eps_b")
            flow_rate_mL = st.number_input("Flow rate (mL/min)", 0.01, 20.0, 1.0, step=0.1, key="m3_flow")

        # Pressure preview
        if "m2_result" in st.session_state:
            from emulsim.module3_performance.hydrodynamics import ColumnGeometry as _CG_preview
            m2r_prev = st.session_state["m2_result"]
            _preview_col = _CG_preview(
                diameter=col_diam_mm / 1000, bed_height=bed_height_cm / 100,
                particle_diameter=m2r_prev.m1_contract.bead_d50,
                bed_porosity=bed_porosity,
                G_DN=m2r_prev.G_DN_updated, E_star=m2r_prev.E_star_updated)
            _dP = _preview_col.pressure_drop(flow_rate_mL / 60e6)
            st.metric("Estimated Pressure Drop", f"{_dP / 1000:.1f} kPa")

        # Mode-specific inputs
        chrom_mode = "Breakthrough"
        C_feed_mg = 1.0
        feed_dur_min = 10.0
        total_time_min = 20.0
        q_max = 100.0
        K_L_m3 = 1000.0
        ext_coeff = 36000.0
        grad_start = 0.0
        grad_end = 500.0
        grad_dur_min = 10.0
        V_max = 1.0
        K_m = 1.0
        S_feed = 10.0
        D_eff = 1e-10
        k_deact = 0.0
        cat_time_h = 1.0

        if m3_app_mode == "Chromatography":
            chrom_mode = st.radio("Mode", ["Breakthrough", "Gradient Elution"], key="m3_chrom_mode")

            st.subheader("Feed")
            C_feed_mg = st.number_input("Feed conc. (mg/mL)", 0.01, 50.0, 1.0, key="m3_Cfeed")
            feed_dur_min = st.number_input("Feed duration (min)", 1.0, 60.0, 10.0, key="m3_feed_dur")
            total_time_min = st.number_input("Total time (min)", 5.0, 120.0, 20.0, key="m3_total_t")

            st.subheader("Isotherm")
            q_max = st.number_input("q_max (mol/m\u00b3)", 1.0, 500.0, 100.0, key="m3_qmax")
            K_L_m3 = st.number_input("K_L (m\u00b3/mol)", 10.0, 1e5, 1000.0, key="m3_KL")
            st.caption("Default isotherm parameters are illustrative \u2014 user calibration required.")

            ext_coeff = st.number_input(
                "\u03b5\u2082\u2088\u2080 (1/(M\u00b7cm))", 1000.0, 200000.0, 36000.0, key="m3_ext"
            )

            if chrom_mode == "Gradient Elution":
                st.subheader("Gradient")
                grad_start = st.number_input("Start (mM)", 0.0, 1000.0, 0.0, key="m3_grad_start")
                grad_end = st.number_input("End (mM)", 0.0, 1000.0, 500.0, key="m3_grad_end")
                grad_dur_min = st.number_input("Duration (min)", 1.0, 60.0, 10.0, key="m3_grad_dur")

        elif m3_app_mode == "Catalysis":
            st.subheader("Enzyme Kinetics")
            V_max = st.number_input("V_max (mol/(m\u00b3\u00b7s))", 0.001, 100.0, 1.0, key="m3_Vmax")
            K_m = st.number_input("K_m (mM)", 0.01, 100.0, 1.0, key="m3_Km")
            S_feed = st.number_input("Substrate feed (mM)", 0.1, 100.0, 10.0, key="m3_Sfeed")
            D_eff = st.number_input("D_eff (m\u00b2/s)", 1e-12, 1e-8, 1e-10, format="%.1e", key="m3_Deff")
            k_deact = st.number_input("k_d (1/s)", 0.0, 1e-3, 0.0, format="%.1e", key="m3_kd")
            cat_time_h = st.number_input("Sim. time (h)", 0.1, 48.0, 1.0, key="m3_cat_t")

        # ── Run M3 Button ────────────────────────────────────────────────
        st.divider()
        _m3_can_run = "m2_result" in st.session_state
        m3_run_btn = st.button("\u25b6 Run M3: Performance Simulation", type="primary",
                                use_container_width=True, disabled=not _m3_can_run)

        if m3_run_btn:
            with st.spinner("Running Module 3: Performance simulation..."):
                from emulsim.module3_performance.hydrodynamics import ColumnGeometry
                m2r = st.session_state["m2_result"]
                column = ColumnGeometry(
                    diameter=col_diam_mm / 1000,
                    bed_height=bed_height_cm / 100,
                    particle_diameter=m2r.m1_contract.bead_d50,
                    bed_porosity=bed_porosity,
                    particle_porosity=m2r.m1_contract.porosity,
                    G_DN=m2r.G_DN_updated,
                    E_star=m2r.E_star_updated,
                )
                st.session_state.pop("m3_result_bt", None)
                st.session_state.pop("m3_result_ge", None)
                st.session_state.pop("m3_result_cat", None)

                if m3_app_mode == "Chromatography":
                    from emulsim.module3_performance import (
                        run_breakthrough,
                        CompetitiveLangmuirIsotherm, run_gradient_elution,
                    )
                    from emulsim.module3_performance.gradient import make_linear_gradient
                    from emulsim.visualization.ui_validators import validate_m3_chromatography, validate_m3_result
                    MW_Da = 66500.0
                    C_feed_mol = C_feed_mg / (MW_Da * 1e-3)

                    _is_grad = (chrom_mode == "Gradient Elution")
                    _col_val = validate_m3_chromatography(
                        flow_rate=flow_rate_mL / 60e6, column=column,
                        isotherm_type="competitive_langmuir" if _is_grad else "langmuir",
                        gradient_enabled=_is_grad,
                    )
                    for _blk in _col_val.blockers:
                        st.warning(f"M3 input blocker: {_blk}")
                    for _wrn in _col_val.warnings:
                        st.info(f"M3 note: {_wrn}")

                    if not _col_val.blockers:
                        try:
                            if chrom_mode == "Breakthrough":
                                bt = run_breakthrough(
                                    column, microsphere=m2r,
                                    C_feed=C_feed_mol, flow_rate=flow_rate_mL / 60e6,
                                    feed_duration=feed_dur_min * 60,
                                    total_time=total_time_min * 60,
                                    extinction_coeff=ext_coeff,
                                )
                                st.session_state["m3_result_bt"] = bt
                                _res_val = validate_m3_result(
                                    mass_balance_error=bt.mass_balance_error,
                                    pressure_drop=bt.pressure_drop,
                                )
                                st.session_state["m3_result_val"] = _res_val
                            else:
                                gradient = make_linear_gradient(
                                    grad_start / 1000.0, grad_end / 1000.0, 0, grad_dur_min * 60
                                )
                                comp_iso = CompetitiveLangmuirIsotherm(
                                    q_max=np.array([q_max]),
                                    K_L=np.array([K_L_m3]),
                                )
                                ge = run_gradient_elution(
                                    column,
                                    C_feed=np.array([C_feed_mol]),
                                    gradient=gradient,
                                    flow_rate=flow_rate_mL / 60e6,
                                    total_time=total_time_min * 60,
                                    isotherm=comp_iso,
                                    feed_duration=feed_dur_min * 60,
                                )
                                st.session_state["m3_result_ge"] = ge
                        except Exception as _m3_ex:
                            st.error(f"Module 3 chromatography failed: {_m3_ex}")

                elif m3_app_mode == "Catalysis":
                    from emulsim.module3_performance.catalysis.packed_bed import solve_packed_bed
                    K_m_mol = K_m
                    S_feed_mol = S_feed
                    try:
                        cat = solve_packed_bed(
                            bed_length=bed_height_cm / 100,
                            bed_diameter=col_diam_mm / 1000,
                            particle_diameter=m2r.m1_contract.bead_d50,
                            bed_porosity=bed_porosity,
                            particle_porosity=m2r.m1_contract.porosity,
                            V_max=V_max, K_m=K_m_mol, S_feed=S_feed_mol,
                            flow_rate=flow_rate_mL / 60e6,
                            D_eff=D_eff, k_deact=k_deact,
                            total_time=cat_time_h * 3600,
                        )
                        st.session_state["m3_result_cat"] = cat
                    except Exception as _m3_ex:
                        st.error(f"Module 3 catalysis failed: {_m3_ex}")
            st.rerun()

        # ── M3 Results Display ───────────────────────────────────────────
        from emulsim.visualization.plots_m3 import (
            plot_chromatogram, plot_breakthrough_curve, plot_peak_table,
            plot_michaelis_menten, plot_effectiveness_factor, plot_activity_decay,
            plot_conversion_vs_time, plot_pressure_flow_curve,
        )

        _show_m3_bt = "m3_result_bt" in st.session_state
        _show_m3_ge = "m3_result_ge" in st.session_state
        _show_m3_cat = "m3_result_cat" in st.session_state

        if _show_m3_bt or _show_m3_ge or _show_m3_cat:
            st.divider()
            st.header("\U0001f4ca M3 Results")

            # Build sub-tabs for M3 results
            _m3_sub_labels = []
            if _show_m3_bt:
                _m3_sub_labels.append("\U0001f4ca Breakthrough")
            if _show_m3_ge:
                _m3_sub_labels.append("\U0001f4ca Gradient Elution")
            if _show_m3_cat:
                _m3_sub_labels.append("\u2697\ufe0f Catalysis")

            _m3_subs = st.tabs(_m3_sub_labels)
            _m3_idx = 0

            if _show_m3_bt:
                with _m3_subs[_m3_idx]:
                    from emulsim.visualization.ui_validators import validate_m3_result as _val_m3_res
                    _bt = st.session_state["m3_result_bt"]
                    st.subheader("Breakthrough Chromatography")

                    _mb_pct = abs(_bt.mass_balance_error) * 100.0
                    if _mb_pct > 5.0:
                        st.error(f"Mass balance error = {_mb_pct:.1f}% \u2014 results numerically unreliable.")
                    elif _mb_pct > 2.0:
                        st.warning(f"Mass balance error = {_mb_pct:.1f}% \u2014 treat with caution.")
                    else:
                        st.success(f"Mass balance error = {_mb_pct:.2f}% \u2014 acceptable.")

                    _bt_val = _val_m3_res(mass_balance_error=_bt.mass_balance_error,
                                          pressure_drop=_bt.pressure_drop)
                    if _bt_val.blockers:
                        for _blk in _bt_val.blockers:
                            st.error(f"Blocker: {_blk}")

                    _dbc_c1, _dbc_c2, _dbc_c3, _dbc_c4 = st.columns(4)
                    _dbc_c1.metric("DBC\u2085%", f"{_bt.dbc_5pct:.1f} mol/m\u00b3" if not np.isnan(_bt.dbc_5pct) else "N/A")
                    _dbc_c2.metric("DBC\u2081\u2080%", f"{_bt.dbc_10pct:.1f} mol/m\u00b3" if not np.isnan(_bt.dbc_10pct) else "N/A")
                    _dbc_c3.metric("DBC\u2085\u2080%", f"{_bt.dbc_50pct:.1f} mol/m\u00b3" if not np.isnan(_bt.dbc_50pct) else "N/A")
                    _dbc_c4.metric("Pressure drop", f"{_bt.pressure_drop / 1000:.1f} kPa")

                    st.divider()
                    _MW_Da_bt = 66500.0
                    _C_feed_mol_bt = C_feed_mg / (_MW_Da_bt * 1e-3)
                    st.plotly_chart(
                        plot_breakthrough_curve(
                            time=_bt.time, C_outlet=_bt.C_outlet, C_feed=_C_feed_mol_bt,
                            dbc_5=_bt.dbc_5pct, dbc_10=_bt.dbc_10pct, dbc_50=_bt.dbc_50pct,
                        ), use_container_width=True,
                    )
                    st.plotly_chart(
                        plot_chromatogram(time=_bt.time, uv_signal=_bt.uv_signal),
                        use_container_width=True,
                    )

                    if "m2_result" in st.session_state:
                        from emulsim.module3_performance.hydrodynamics import ColumnGeometry as _CG_bt
                        _m2r_bt = st.session_state["m2_result"]
                        _col_bt = _CG_bt(
                            diameter=col_diam_mm / 1000, bed_height=bed_height_cm / 100,
                            particle_diameter=_m2r_bt.m1_contract.bead_d50,
                            bed_porosity=bed_porosity,
                            particle_porosity=_m2r_bt.m1_contract.porosity,
                            G_DN=_m2r_bt.G_DN_updated, E_star=_m2r_bt.E_star_updated,
                        )
                        st.plotly_chart(plot_pressure_flow_curve(_col_bt), use_container_width=True)

                    st.caption(
                        "Mechanistic prediction: Lumped Rate Model (LRM) with Langmuir isotherm. "
                        "DBC values are model-based \u2014 calibrate isotherm parameters with batch uptake experiments."
                    )
                _m3_idx += 1

            if _show_m3_ge:
                with _m3_subs[_m3_idx]:
                    _ge = st.session_state["m3_result_ge"]
                    st.subheader("Gradient Elution Chromatography")

                    _grad_affects = getattr(_ge, "gradient_affects_binding", False)
                    if _grad_affects:
                        st.success("Gradient affects binding: YES - selected isotherm is gradient-sensitive.")
                    else:
                        st.info("Gradient affects binding: NO - plain competitive Langmuir is diagnostic/display only.")

                    _ge_time = getattr(_ge, "time", None)
                    _ge_uv = getattr(_ge, "uv_signal", None)
                    _ge_grad = getattr(_ge, "gradient_values", None)
                    if _ge_time is not None and _ge_uv is not None:
                        st.plotly_chart(
                            plot_chromatogram(
                                time=_ge_time, uv_signal=_ge_uv,
                                gradient_values=_ge_grad, gradient_affects_binding=_grad_affects,
                            ), use_container_width=True,
                        )

                    _ge_peaks = getattr(_ge, "peaks", [])
                    st.plotly_chart(plot_peak_table(_ge_peaks), use_container_width=True)

                    st.caption(
                        "Ranking only - gradient-sensitive adapters update binding during elution; "
                        "plain competitive Langmuir remains diagnostic. Quantitative yields require "
                        "isotherm calibration."
                    )
                _m3_idx += 1

            if _show_m3_cat:
                with _m3_subs[_m3_idx]:
                    _cat = st.session_state["m3_result_cat"]
                    st.subheader("Packed-Bed Catalytic Reactor")

                    _cat_c1, _cat_c2, _cat_c3, _cat_c4 = st.columns(4)
                    _cat_c1.metric("Final Conversion", f"{_cat.conversion:.1%}")
                    _cat_c2.metric("Effectiveness Factor \u03b7", f"{_cat.effectiveness_factor:.3f}")
                    _cat_c3.metric("Thiele Modulus \u03a6", f"{_cat.thiele_modulus:.2f}")
                    _cat_c4.metric("Productivity", f"{_cat.productivity:.3e} mol/(m\u00b3\u00b7s)")

                    _mb_cat_pct = abs(_cat.mass_balance_error) * 100.0
                    if _mb_cat_pct > 5.0:
                        st.error(f"Catalysis mass balance error = {_mb_cat_pct:.1f}% \u2014 results unreliable.")
                    elif _mb_cat_pct > 2.0:
                        st.warning(f"Catalysis mass balance error = {_mb_cat_pct:.1f}% \u2014 treat with caution.")
                    else:
                        st.success(f"Catalysis mass balance error = {_mb_cat_pct:.2f}% \u2014 acceptable.")

                    st.divider()
                    _time_h_cat = _cat.time / 3600.0
                    _S_in = float(_cat.S_outlet[0]) if len(_cat.S_outlet) > 0 else 1.0
                    if _S_in <= 0:
                        _S_in = S_feed
                    _conversion_arr = np.clip(1.0 - _cat.S_outlet / max(_S_in, 1e-12), 0.0, 1.0)
                    st.plotly_chart(
                        plot_conversion_vs_time(time_hours=_time_h_cat, conversion=_conversion_arr),
                        use_container_width=True,
                    )

                    _S_range_mm = np.linspace(0.01, max(S_feed * 3, 10.0), 200)
                    st.plotly_chart(
                        plot_michaelis_menten(
                            S_range=_S_range_mm, V_max=V_max, K_m=K_m,
                            eta=_cat.effectiveness_factor,
                        ), use_container_width=True,
                    )

                    _phi_range = np.logspace(-2, 2, 200)
                    st.plotly_chart(plot_effectiveness_factor(_phi_range), use_container_width=True)

                    if hasattr(_cat, "activity_history") and len(_cat.activity_history) > 0:
                        st.plotly_chart(
                            plot_activity_decay(time_hours=_time_h_cat, activity=_cat.activity_history),
                            use_container_width=True,
                        )

                    _phi_val = _cat.thiele_modulus
                    if _phi_val < 0.3:
                        st.success(f"Thiele modulus \u03a6 = {_phi_val:.2f} \u2014 kinetic regime (no diffusion limitation).")
                    elif _phi_val < 3.0:
                        st.warning(f"Thiele modulus \u03a6 = {_phi_val:.2f} \u2014 transition regime (partial diffusion limitation).")
                    else:
                        st.error(f"Thiele modulus \u03a6 = {_phi_val:.2f} \u2014 diffusion-limited regime. "
                                 "Effectiveness factor \u03b7 is significantly < 1. Consider smaller particles or higher D_eff.")

                    st.caption(
                        "Mechanistic prediction: Transient PFR with Michaelis-Menten kinetics, "
                        "Thiele modulus effectiveness factor, and first-order deactivation. "
                        "K_m and V_max require calibration against your specific enzyme."
                    )
