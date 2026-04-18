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
    ModelMode, ModelEvidenceTier,
)


def _evidence_badge(result_obj) -> str:
    """Return a colored evidence tier badge string for UI display."""
    manifest = getattr(result_obj, "model_manifest", None)
    if manifest is None:
        return ""
    tier = manifest.evidence_tier
    _COLORS = {
        ModelEvidenceTier.VALIDATED_QUANTITATIVE: ":green[VALIDATED]",
        ModelEvidenceTier.CALIBRATED_LOCAL: ":green[CALIBRATED]",
        ModelEvidenceTier.SEMI_QUANTITATIVE: ":orange[SEMI-QUANTITATIVE]",
        ModelEvidenceTier.QUALITATIVE_TREND: ":red[QUALITATIVE TREND]",
        ModelEvidenceTier.UNSUPPORTED: ":red[UNSUPPORTED]",
    }
    badge = _COLORS.get(tier, str(tier.value))
    model_name = manifest.model_name
    return f"Evidence: {badge} | Model: `{model_name}`"
from emulsim.properties.database import PropertyDatabase
from emulsim.pipeline.orchestrator import PipelineOrchestrator
from emulsim.trust import assess_trust
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


def _render_non_ac_family(*, tab_container, family, is_stirred_default, model_mode_enum, _smgr) -> None:
    """Render the M1 tab for alginate / cellulose / PLGA (v9.0 M5-M7).

    Shares Hardware Mode + shared L1 inputs (RPM, time, phi_d, v_oil/v_poly)
    with the A+C path, but dispatches to the family-specific formulation
    module and builds SimulationParameters with polymer_family set so the
    orchestrator routes through _run_alginate / _run_cellulose / _run_plga.
    """
    from emulsim.visualization.tabs.m1.hardware_section import render_hardware_mode_radio
    from emulsim.datatypes import PolymerFamily

    m1_col_left, m1_col_right = st.columns(2)
    with m1_col_left:
        st.subheader("Emulsification (L1)")
        is_stirred = render_hardware_mode_radio()

        if is_stirred:
            vessel_choice = st.selectbox(
                "Vessel", ["Glass Beaker (100 mm)", "Jacketed Vessel (92 mm)"],
                key="m1_vessel",
            )
            stirrer_choice = st.selectbox(
                "Stirrer", ["Stirrer A - Pitched Blade (59 mm)", "Stirrer B - Rotor-Stator (32 mm)"],
                key="m1_stirrer",
            )
            is_stirrer_A = "Pitched" in stirrer_choice
            rpm = st.slider(
                "Stirrer Speed (RPM)", 800, 2000 if is_stirrer_A else 9000,
                1300 if is_stirrer_A else 1800,
                step=50 if is_stirrer_A else 100,
                key="m1_rpm" if is_stirrer_A else "m1_rpm_rs",
            )
            t_emul = st.number_input("Emulsification Time (min)", 1, 60, 10, key="m1_t_emul")
            v_oil_mL = st.slider("Oil + surfactant (mL)", 100, 500, 300, step=10, key="m1_v_oil")
            v_poly_mL = st.slider("Dispersed phase (mL)", 50, 300, 200, step=10, key="m1_v_poly")
            total_mL = v_oil_mL + v_poly_mL
            phi_d = v_poly_mL / total_mL
            st.caption(f"Total: {total_mL} mL | phi_d = {phi_d:.2f}")
        else:
            rpm = st.slider("Rotor Speed (RPM)", 3000, 25000, 10000, step=500, key="m1_rpm_leg")
            t_emul = st.number_input("Emulsification Time (min)", 1, 60, 10, key="m1_t_emul_leg")
            phi_d = st.slider("Dispersed Phase Fraction (phi_d)", 0.01, 0.30, 0.05, step=0.01, key="m1_phi_d")
            v_oil_mL, v_poly_mL, total_mL = 300, 200, 500

    # Compare families by .value to survive importlib.reload of datatypes
    # (see tab_m1.py family-dispatch comment and RunReport.compute_min_tier).
    family_value = getattr(family, "value", family)
    with m1_col_right:
        if family_value == PolymerFamily.ALGINATE.value:
            from emulsim.visualization.tabs.m1.formulation_alginate import render_formulation_alginate
            fctx = render_formulation_alginate(is_stirred=is_stirred)
        elif family_value == PolymerFamily.CELLULOSE.value:
            from emulsim.visualization.tabs.m1.formulation_cellulose import render_formulation_cellulose
            fctx = render_formulation_cellulose(is_stirred=is_stirred)
        elif family_value == PolymerFamily.PLGA.value:
            from emulsim.visualization.tabs.m1.formulation_plga import render_formulation_plga
            fctx = render_formulation_plga(is_stirred=is_stirred)
        else:
            st.error(f"Unknown family: {family_value}")
            return

    st.divider()
    run_btn = st.button(
        "\u25b6 Run M1: Fabrication Pipeline", type="primary",
        use_container_width=True, key="m1v9_run_non_ac",
    )
    if not run_btn:
        return

    # ── Build SimulationParameters for the selected family ──────────────
    if is_stirred:
        if "Beaker" in vessel_choice:
            vessel = VesselGeometry.glass_beaker(working_volume=total_mL * 1e-6)
            heating = HeatingConfig.flat_plate()
        else:
            vessel = VesselGeometry.jacketed_vessel(working_volume=total_mL * 1e-6)
            heating = HeatingConfig.hot_water_jacket()
        stirrer = StirrerGeometry.pitched_blade_A() if is_stirrer_A else StirrerGeometry.rotor_stator_B()
        heating.T_initial = fctx.T_oil_C + 273.15
        kernels = KernelConfig.for_stirrer_type(stirrer.stirrer_type)
        emul = EmulsificationParameters(
            mode="stirred_vessel", rpm=float(rpm),
            t_emulsification=float(t_emul * 60),
            vessel=vessel, stirrer=stirrer, heating=heating, kernels=kernels,
        )
    else:
        emul = EmulsificationParameters(
            mode="rotor_stator_legacy", rpm=float(rpm),
            t_emulsification=float(t_emul * 60),
            mixer=MixerGeometry(),
        )

    formulation = FormulationParameters(
        c_span80=fctx.c_span80_pct * 10.0,
        T_oil=fctx.T_oil_C + 273.15,
        phi_d=phi_d,
        c_span80_vol_pct=fctx.c_span80_vol_pct,
        v_oil_span80_mL=float(v_oil_mL),
        v_polysaccharide_mL=float(v_poly_mL),
    )

    props_overrides: dict = {"polymer_family": family}
    if family_value == PolymerFamily.ALGINATE.value:
        formulation.c_alginate = fctx.c_alginate_kg_m3
        if fctx.c_Ca_bath_mM > 0:
            # External CaCl2 bath: direct concentration.
            formulation.c_Ca_bath = fctx.c_Ca_bath_mM
        else:
            # Internal release (GDL+CaCO3): lumped-parameter approximation
            # c_eff(t) = C_source * (1 - exp(-k*t)) per Draget 1997. Use the
            # gelant profile's recommended t_default as the process time.
            from emulsim.reagent_library_alginate import (
                GELANTS_ALGINATE,
                effective_bath_concentration,
            )
            profile = GELANTS_ALGINATE.get(fctx.gelant_key)
            if profile is not None and profile.mode == "internal_release":
                # Recompute with the UI-provided k_release (which overrides
                # profile.k_release when the user adjusts it) by building a
                # synthetic profile with the current widget values.
                from dataclasses import replace
                live_profile = replace(
                    profile,
                    C_Ca_source=fctx.C_Ca_source_mM,
                    k_release=fctx.k_release_1_s,
                )
                formulation.c_Ca_bath = effective_bath_concentration(
                    live_profile, profile.t_default,
                )
            else:
                formulation.c_Ca_bath = fctx.C_Ca_source_mM
    elif family_value == PolymerFamily.CELLULOSE.value:
        formulation.phi_cellulose_0 = fctx.phi_cellulose_0
        formulation.solvent_system = fctx.solvent_system
        formulation.cooling_rate = fctx.cooling_rate_Cmin / 60.0 if fctx.cooling_rate_Cmin > 0 else 0.0
    elif family_value == PolymerFamily.PLGA.value:
        formulation.phi_PLGA_0 = fctx.phi_PLGA_0
        formulation.plga_grade = fctx.plga_grade

    params = SimulationParameters(
        model_mode=model_mode_enum,
        emulsification=emul,
        formulation=formulation,
        solver=SolverSettings(l2_n_grid=64),
    )

    _smgr.invalidate_downstream(from_module=1)
    with st.spinner(f"Running L1→L2→L4 pipeline for {family.value}..."):
        t_start = time.time()
        db = PropertyDatabase()
        orch = PipelineOrchestrator(db=db)
        try:
            result = orch.run_single(
                params, l2_mode="empirical",
                props_overrides=props_overrides,
                crosslinker_key="genipin",  # ignored for non-A+C families
                uv_intensity=0.0,
            )
        except Exception as ex:
            st.error(f"Simulation failed: {ex}")
            st.exception(ex)
            return
        elapsed = time.time() - t_start

    st.success(f"Pipeline complete in {elapsed:.1f}s")

    # ── Family-neutral results display ────────────────────────────────
    e, g, m = result.emulsification, result.gelation, result.mechanical
    c1, c2, c3, c4 = st.columns(4)
    d_mode = getattr(e, "d_mode", 0.0)
    c1.metric("d_mode", f"{d_mode*1e6:.1f} µm" if d_mode > 0 else f"d32 {e.d32*1e6:.2f} µm")
    c2.metric("Pore size", f"{g.pore_size_mean*1e9:.1f} nm")
    c3.metric("Porosity", f"{g.porosity:.2f}")
    c4.metric("G (modulus)", f"{m.G_DN/1000:.1f} kPa")

    if getattr(e, "model_manifest", None):
        st.caption(f"L1 model: `{e.model_manifest.model_name}` ({e.model_manifest.evidence_tier.value})")
    if getattr(g, "model_manifest", None):
        st.caption(f"L2 model: `{g.model_manifest.model_name}` ({g.model_manifest.evidence_tier.value})")
    if getattr(m, "model_manifest", None):
        st.caption(f"L4 model: `{m.model_manifest.model_name}` ({m.model_manifest.evidence_tier.value})")

    # Save to session state so tabs M2/M3 can consume the result
    st.session_state["result"] = result
    st.session_state["elapsed"] = elapsed
    st.session_state["params"] = params


def render_tab_m1(
    tab_container,
    is_stirred: bool | None,
    model_mode_enum: ModelMode,
    _smgr,
    _const_input,
    _proto_sections: dict,
) -> None:
    """Render the M1 Fabrication tab inside the given Streamlit container.

    Args:
        tab_container: Streamlit tab container from st.tabs().
        is_stirred: Whether stirred vessel mode is selected. Pass ``None``
            (v9.0 default) to render the Hardware Mode radio locally at
            the top of the Emulsification section.
        model_mode_enum: Scientific mode enum.
        _smgr: SessionStateManager instance.
        _const_input: Callable for rendering per-constant selector with protocol link.
        _proto_sections: Dict of calibration protocol sections.
    """
    from emulsim.visualization.tabs.m1.hardware_section import render_hardware_mode_radio
    from emulsim.visualization.tabs.m1.family_selector import render_family_selector
    from emulsim.datatypes import PolymerFamily as _PF

    with tab_container:
        st.header("Module 1: Fabrication Pipeline (L1→L2→L3→L4)")

        # ── Polymer family (v9.0) — drives L2 dispatch and downstream rendering ──
        _family_ctx = render_family_selector()
        _family = _family_ctx.family

        # v9.0 M5-M7: dispatch non-A+C families to the family-specific runner.
        # Compare by .value (string), not enum identity — app.py reloads
        # emulsim.datatypes on every rerun, producing a new PolymerFamily
        # class. m1.family_selector is not in the reload list, so it returns
        # a stale-class enum member. Identity comparison would mis-classify
        # every family after the first rerun. See also
        # RunReport.compute_min_tier for the same hazard.
        if getattr(_family, "value", _family) != _PF.AGAROSE_CHITOSAN.value:
            _render_non_ac_family(
                tab_container=tab_container,
                family=_family,
                is_stirred_default=None,  # hardware_mode_radio already rendered in app.py sidebar v8 or by us in M2
                model_mode_enum=model_mode_enum,
                _smgr=_smgr,
            )
            return

        # ── M1 INPUT SECTION ─────────────────────────────────────────────────
        m1_col_left, m1_col_right = st.columns(2)

        # ── Left column: Emulsification + Formulation ──
        with m1_col_left:
            st.subheader("Emulsification (L1)")

            # v9.0: Hardware Mode relocated from Global Settings sidebar
            # into the M1 Emulsification section (see scientific-advisor
            # audit §C). Back-compat: if caller still passes a bool,
            # honour it; otherwise render the radio here.
            if is_stirred is None:
                is_stirred = render_hardware_mode_radio()

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
                    st.caption("Heating: flat-plate (150C -> 80C oil)")
                else:
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

            # v9.0 M4: formulation + cooling + pore-structure moved into
            # formulation_agarose_chitosan.render_formulation_section.
            from emulsim.visualization.tabs.m1.formulation_agarose_chitosan import (
                render_formulation_section as _render_ac_formulation,
            )
            _ac_ctx = _render_ac_formulation(is_stirred=is_stirred)
            c_agarose_pct = _ac_ctx.c_agarose_pct
            c_chitosan_pct = _ac_ctx.c_chitosan_pct
            _surf_sel_key = _ac_ctx.surfactant_key
            c_span80_pct = _ac_ctx.c_span80_pct
            c_span80_vol_pct = _ac_ctx.c_span80_vol_pct
            T_oil_C = _ac_ctx.T_oil_C
            cooling_rate_Cmin = _ac_ctx.cooling_rate_Cmin
            l2_mode = _ac_ctx.l2_mode
            grid_size = _ac_ctx.grid_size
            surf = _ac_ctx.surfactant

        # ── Right column: Crosslinking + Targets + Material Constants ──
        with m1_col_right:
            # v9.0 M4: crosslinking section moved into module.
            from emulsim.visualization.tabs.m1.crosslinking_section import (
                render_crosslinking_section as _render_crosslinking,
            )
            _DDA_out: list = []
            _xl_ctx = _render_crosslinking(
                c_chitosan_pct=c_chitosan_pct, DDA_out=_DDA_out,
            )
            _xl_sel_key = _xl_ctx.crosslinker_key
            c_genipin_mM = _xl_ctx.c_genipin_mM
            T_xlink_C = _xl_ctx.T_xlink_C
            t_xlink_h = _xl_ctx.t_xlink_h
            uv_intensity = _xl_ctx.uv_intensity
            xl = _xl_ctx.crosslinker
            _DDA = _DDA_out[0] if _DDA_out else 0.85

            # v9.0 M4: targets moved into module.
            from emulsim.visualization.tabs.m1.targets_section import (
                render_targets_section as _render_targets,
            )
            _tgt_ctx = _render_targets(family=_family, is_stirred=is_stirred)
            target_d32 = _tgt_ctx.target_d32
            target_d_mode = _tgt_ctx.target_d_mode
            target_pore = _tgt_ctx.target_pore
            target_G = _tgt_ctx.target_G

            # v9.0 M4: material constants moved into module (family-aware).
            from emulsim.visualization.tabs.m1.material_constants import (
                render_material_constants as _render_material_constants,
            )
            _mat_overrides = _render_material_constants(
                family=_family,
                surfactant=surf,
                crosslinker=xl,
                const_input=_const_input,
                T_oil_C=T_oil_C,
                c_span80_pct=c_span80_pct,
            )

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

        # Material constants overrides (v9.0 M4: sourced from material_constants module)
        _custom_props_overrides = dict(_mat_overrides)
        _custom_props_overrides['k_xlink_0'] = xl.k_xlink_0
        _custom_props_overrides['E_a_xlink'] = xl.E_a_xlink
        _custom_props_overrides['f_bridge'] = xl.f_bridge

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
                    st.exception(ex)
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
                _l1_badge = _evidence_badge(e)
                if _l1_badge:
                    st.caption(_l1_badge)
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
                _l2_badge = _evidence_badge(g)
                if _l2_badge:
                    st.caption(_l2_badge)

            with sub4:
                st.subheader("Level 3: Crosslinking Kinetics")
                st.plotly_chart(plot_crosslinking_kinetics(x), use_container_width=True)
                c1, c2, c3 = st.columns(3)
                c1.write(f"**Crosslink fraction** = {x.p_final:.4f}")
                c1.write(f"**G_crosslinked** = {x.G_chitosan_final:.0f} Pa")
                c2.write(f"**Mesh size \u03be** = {x.xi_final*1e9:.1f} nm")
                c2.write(f"**M_c** = {x.Mc_final:.0f} g/mol")
                c3.write(f"**\u03bd_e** = {x.nu_e_final:.2e} /m\u00b3")
                # v6.1: L3 diagnostics
                if getattr(x, 'regime', 'unknown') != 'unknown':
                    st.caption(f"Regime: {x.regime} | Thiele: {x.thiele_modulus:.2f} | Stoich. ceiling: {x.stoichiometric_ceiling:.2f}")
                _l3_badge = _evidence_badge(x)
                if _l3_badge:
                    st.caption(_l3_badge)

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
                # v6.1: network classification + evidence badge
                _ntype = getattr(m, 'network_type', 'unknown')
                if _ntype != 'unknown':
                    st.caption(f"Network type: {_ntype}")
                _l4_badge = _evidence_badge(m)
                if _l4_badge:
                    st.caption(_l4_badge)

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
