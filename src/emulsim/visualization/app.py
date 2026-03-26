"""Streamlit UI for the emulsification simulation system.

Launch: streamlit run src/emulsim/visualization/app.py
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import streamlit as st

# Ensure the package is importable and force-reload to pick up code changes
_root = Path(__file__).resolve().parents[3]
if str(_root / "src") not in sys.path:
    sys.path.insert(0, str(_root / "src"))

# Force reload of core modules to avoid Streamlit's import cache
import importlib
import emulsim.pipeline.orchestrator as _orch_mod
importlib.reload(_orch_mod)

from emulsim.datatypes import SimulationParameters, MaterialProperties, FormulationParameters, EmulsificationParameters, MixerGeometry, SolverSettings
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
)

# ─── Page Config ──────────────────────────────────────────────────────────

st.set_page_config(
    page_title="EmulSim — Microsphere Simulation",
    page_icon="🔬",
    layout="wide",
)

st.title("🔬 EmulSim — Double-Network Hydrogel Microsphere Simulator")
st.caption("Multi-scale simulation: Emulsification → Gelation → Crosslinking → Mechanical Properties")

# ─── Sidebar: Parameter Input ─────────────────────────────────────────────

st.sidebar.header("⚙️ Process Parameters")

st.sidebar.subheader("Emulsification (L1)")
rpm = st.sidebar.slider("Rotor Speed (RPM)", 3000, 25000, 10000, step=500)
t_emul = st.sidebar.number_input("Emulsification Time (min)", 1, 60, 10)
phi_d = st.sidebar.slider("Dispersed Phase Fraction (φ_d)", 0.01, 0.30, 0.05, step=0.01)

st.sidebar.subheader("Formulation")

# Surfactant selector
_surf_keys = list(SURFACTANTS.keys())
_surf_names = [SURFACTANTS[k].name for k in _surf_keys]
_surf_sel_name = st.sidebar.selectbox(
    "Surfactant",
    _surf_names,
    index=_surf_keys.index("span80"),
    help="Select surfactant for W/O emulsification",
)
_surf_sel_key = _surf_keys[_surf_names.index(_surf_sel_name)]
surf = SURFACTANTS[_surf_sel_key]
st.sidebar.caption(f"HLB={surf.hlb} | MW={surf.mw} g/mol | {surf.notes[:60]}")

col1, col2 = st.sidebar.columns(2)
with col1:
    c_agarose_pct = st.number_input("Agarose (% w/v)", 1.0, 10.0, 4.2, step=0.1)
with col2:
    c_chitosan_pct = st.number_input("Chitosan (% w/v)", 0.5, 5.0, 1.8, step=0.1)

c_span80_pct = st.sidebar.slider("Surfactant (% w/v)", 0.5, 5.0, 2.0, step=0.1)
T_oil_C = st.sidebar.slider("Oil Temperature (°C)", 60, 95, 90)

st.sidebar.subheader("Cooling & Gelation (L2)")
cooling_rate_Cmin = st.sidebar.slider("Cooling Rate (°C/min)", 1.0, 20.0, 10.0, step=0.5)
l2_model = st.sidebar.radio(
    "Pore Structure Model",
    ["Empirical (fast, calibrated)", "Cahn-Hilliard 2D (mechanistic)"],
    index=0,
    help="Empirical: literature-calibrated power law (~1 ms). "
         "CH 2D: phase-field spinodal decomposition (slower, requires parameter tuning).",
)
l2_mode = "empirical" if "Empirical" in l2_model else "ch_2d"
grid_size = 64
if l2_mode == "ch_2d":
    grid_size = st.sidebar.select_slider("Phase-Field Grid", [32, 64, 128], value=64)

st.sidebar.subheader("Crosslinking (L3)")

# Crosslinker selector
_xl_keys = list(CROSSLINKERS.keys())
_xl_names = [CROSSLINKERS[k].name for k in _xl_keys]
_xl_sel_name = st.sidebar.selectbox(
    "Crosslinker",
    _xl_names,
    index=_xl_keys.index("genipin"),
    help="Select crosslinking agent",
)
_xl_sel_key = _xl_keys[_xl_names.index(_xl_sel_name)]
xl = CROSSLINKERS[_xl_sel_key]
st.sidebar.caption(f"{xl.mechanism} | k\u2080={xl.k_xlink_0:.1e} | Score: {xl.suitability}/10")

c_genipin_mM = st.sidebar.slider("Crosslinker Concentration (mM)", 0.5, 10.0, 2.0, step=0.5)

# Adapt slider defaults and ranges from crosslinker profile
_T_xlink_default = int(xl.T_crosslink_default - 273.15)
_t_xlink_default_h = max(1, int(xl.t_crosslink_default / 3600))
T_xlink_C = st.sidebar.slider("Crosslinking Temperature (°C)", 0, 120,
                                min(max(_T_xlink_default, 0), 120))
t_xlink_h = st.sidebar.slider("Crosslinking Time (hours)", 1, 48,
                                min(_t_xlink_default_h, 48))

# UV intensity slider (only relevant for UV-initiated crosslinkers)
if xl.kinetics_model == 'uv_dose':
    uv_intensity = st.sidebar.slider("UV Intensity (mW/cm\u00b2)", 1.0, 100.0, 20.0, step=1.0)
else:
    uv_intensity = 0.0

st.sidebar.subheader("Optimization Targets")
target_d32 = st.sidebar.number_input("Target d32 (µm)", 0.5, 50.0, 2.0, step=0.5)
target_pore = st.sidebar.number_input("Target Pore Size (nm)", 10, 500, 80, step=10)
target_G = st.sidebar.number_input("Target G_DN (kPa)", 1.0, 500.0, 10.0, step=1.0)

# ─── Material Constants (per-constant Literature / Custom + protocol link) ─

st.sidebar.divider()
st.sidebar.subheader("Material Constants")
st.sidebar.caption("For each constant: use literature or enter your own. Click name to view calibration protocol.")

from emulsim.literature_constants import ALL_CONSTANTS

# Load the calibration protocol and extract per-study sections
_proto_path = Path(__file__).resolve().parents[3] / "docs" / "04_calibration_protocol.md"
_proto_sections = {}
if _proto_path.exists():
    _proto_text = _proto_path.read_text(encoding="utf-8")
    import re as _re
    # Extract each study section by heading
    _study_map = {
        "K_L":        ("Study 1 -- Interfacial Tension", "## 2. Study 1", "## 3. Study 2"),
        "Gamma_inf":  ("Study 1 -- Interfacial Tension", "## 2. Study 1", "## 3. Study 2"),
        "eta_chit":   ("Study 2 -- Chitosan Viscosity",  "## 3. Study 2", "## 4. Study 3"),
        "C3":         ("Study 3 -- Viscous Breakage",     "## 4. Study 3", "## 5. Study 4"),
        "pore":       ("Study 4 -- Pore Structure",       "## 5. Study 4", "## 6. Study 5"),
        "eta_coup":   ("Study 5 -- IPN Coupling",         "## 6. Study 5", "## 7. Inputting"),
    }
    for key, (title, start_marker, end_marker) in _study_map.items():
        s_idx = _proto_text.find(start_marker)
        e_idx = _proto_text.find(end_marker)
        if s_idx >= 0:
            section = _proto_text[s_idx:e_idx] if e_idx > s_idx else _proto_text[s_idx:]
            _proto_sections[key] = (title, section)


def _const_input(label, key, lit_val, unit, source_short, lo, hi, step,
                  proto_key, fmt="%.3f"):
    """Render per-constant selector with protocol link."""
    # Clickable protocol popover
    if proto_key in _proto_sections:
        title, section_md = _proto_sections[proto_key]
        with st.sidebar.popover(f"📋 {label}"):
            st.markdown(f"### Calibration Protocol: {title}")
            st.markdown(section_md[:3000])  # first 3000 chars to fit popover
            if len(section_md) > 3000:
                st.caption("... (see full protocol in docs/04_calibration_protocol.md)")
    else:
        st.sidebar.markdown(f"**{label}**")

    src = st.sidebar.radio(
        "Source",
        ["Literature", "Custom"],
        index=0, horizontal=True, key=f"src_{key}",
        label_visibility="collapsed",
    )
    if src == "Literature":
        st.sidebar.caption(f"  = {lit_val} {unit}  ({source_short})")
        return lit_val
    else:
        return st.sidebar.number_input(
            f"{label} ({unit})", lo, hi, float(lit_val), step=step,
            format=fmt, key=f"val_{key}",
        )

_kl = ALL_CONSTANTS["K_L"]
custom_K_L = _const_input(
    "K_L (Span-80 adsorption)", "K_L",
    _kl.value, "m³/mol", "Santini 2007", 0.01, 10.0, 0.05, "K_L",
)

_gi = ALL_CONSTANTS["Gamma_inf"]
custom_Gamma_inf = _const_input(
    "Γ∞ (max surface excess)", "Gamma_inf",
    _gi.value * 1e6, "×10⁻⁶ mol/m²", "Santini 2007", 1.0, 10.0, 0.1,
    "Gamma_inf", "%.1f",
)

_ec = ALL_CONSTANTS["eta_intr_chit"]
custom_eta_chit = _const_input(
    "[η] chitosan", "eta_chit",
    _ec.value, "mL/g", "Rinaudo 1993", 100.0, 2000.0, 50.0,
    "eta_chit", "%.0f",
)

_c3 = ALL_CONSTANTS["breakage_C3"]
custom_C3 = _const_input(
    "C3 (viscous breakage)", "C3",
    _c3.value, "-", "Alopaeus 2002", 0.0, 1.0, 0.05,
    "C3", "%.2f",
)

_et = ALL_CONSTANTS["eta_coupling"]
custom_eta_coup = _const_input(
    "η coupling (IPN)", "eta_coup",
    _et.value, "-", "Gong 2010 est.", -0.5, 0.5, 0.05,
    "eta_coup", "%.2f",
)

# ─── Build Parameters ─────────────────────────────────────────────────────

params = SimulationParameters(
    emulsification=EmulsificationParameters(
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
    solver=SolverSettings(l2_n_grid=grid_size),
)

# Apply custom/literature material constants
from emulsim.datatypes import MaterialProperties as _MP
_custom_props_overrides = {
    "breakage_C3": custom_C3,
    "eta_coupling": custom_eta_coup,
}

# ─── Build reagent-library overrides ─────────────────────────────────

# Crosslinker overrides (Level 3 kinetics)
_custom_props_overrides['k_xlink_0'] = xl.k_xlink_0
_custom_props_overrides['E_a_xlink'] = xl.E_a_xlink
_custom_props_overrides['f_bridge'] = xl.f_bridge

# Surfactant IFT override: use custom K_L/Gamma_inf if user set them,
# otherwise use the surfactant library values.
_R_gas = 8.314
_T_ift = T_oil_C + 273.15
_c_mol_ift = c_span80_pct * 10.0 / surf.mw * 1000.0  # kg/m3 -> mol/m3
_sigma_0_T = max(surf.sigma_0_paraffin - 1.0e-4 * (_T_ift - 293.15), 0.001)
_use_gamma = custom_Gamma_inf * 1e-6  # convert from display units (×10⁻⁶) to mol/m²
_use_KL = custom_K_L
_sigma_calc = max(
    _sigma_0_T - _R_gas * _T_ift * _use_gamma * np.log(1 + _use_KL * _c_mol_ift),
    1e-4,
)
_custom_props_overrides['sigma'] = _sigma_calc

# Chitosan intrinsic viscosity override — now a MaterialProperties field,
# so the PropertyDatabase will use it when computing dispersed phase viscosity.
_custom_props_overrides['eta_intr_chit'] = custom_eta_chit

# ─── Run Simulation ───────────────────────────────────────────────────────

st.divider()

run_btn = st.button("▶ Run Simulation", type="primary", use_container_width=True)

if run_btn:
    with st.spinner("Running L1→L2→L3→L4 pipeline..."):
        t_start = time.time()

        db = PropertyDatabase()
        orch = PipelineOrchestrator(db=db)

        # Progress tracking
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

        progress.progress(100, text=f"Complete in {elapsed:.1f}s")

    st.session_state["result"] = result
    st.session_state["elapsed"] = elapsed
    st.session_state["params"] = params
    st.session_state["targets"] = (target_d32, target_pore, target_G)

    # Compute trust assessment using the same overrides as the simulation
    db_trust = PropertyDatabase()
    props_trust = db_trust.update_for_conditions(
        params.formulation.T_oil, params.formulation.c_agarose,
        params.formulation.c_chitosan, params.formulation.c_span80,
    )
    for _k, _v in _custom_props_overrides.items():
        if hasattr(props_trust, _k):
            setattr(props_trust, _k, _v)
    st.session_state["trust"] = assess_trust(result, params, props_trust,
                                                crosslinker_key=_xl_sel_key)

# ─── Display Results ──────────────────────────────────────────────────────

if "result" in st.session_state:
    result = st.session_state["result"]
    elapsed = st.session_state["elapsed"]
    target_d32, target_pore, target_G = st.session_state["targets"]

    e = result.emulsification
    g = result.gelation
    x = result.crosslinking
    m = result.mechanical

    st.header("📊 Results Summary")

    # KPI cards
    col1, col2, col3, col4 = st.columns(4)

    d32_dev = abs(e.d32 * 1e6 - target_d32) / target_d32 * 100
    pore_dev = abs(g.pore_size_mean * 1e9 - target_pore) / target_pore * 100
    G_dev = abs(m.G_DN / 1000 - target_G) / target_G * 100

    col1.metric("d32", f"{e.d32*1e6:.2f} µm",
                delta=f"{d32_dev:.0f}% from target",
                delta_color="inverse")
    col2.metric("Pore Size", f"{g.pore_size_mean*1e9:.1f} nm",
                delta=f"{pore_dev:.0f}% from target",
                delta_color="inverse")
    col3.metric("G_DN", f"{m.G_DN/1000:.1f} kPa",
                delta=f"{G_dev:.0f}% from target",
                delta_color="inverse")
    col4.metric("Pipeline Time", f"{elapsed:.1f}s")

    # Additional metrics
    col5, col6, col7, col8 = st.columns(4)
    col5.metric("Span", f"{e.span:.2f}")
    col6.metric("Porosity", f"{g.porosity:.1%}")
    col7.metric("Crosslink %", f"{x.p_final:.1%}")
    col8.metric("E*", f"{m.E_star/1000:.1f} kPa")

    st.divider()

    # ── Detailed Results by Level ─────────────────────────────────────

    tab1, tab2, tab3, tab4, tab5 = st.tabs([
        "📈 Dashboard", "🫧 L1: Emulsification", "🧊 L2: Gelation",
        "🔗 L3: Crosslinking", "💪 L4: Mechanical",
    ])

    with tab1:
        st.plotly_chart(plot_results_dashboard(result), use_container_width=True)

    with tab2:
        st.subheader("Level 1: Emulsification — Droplet Size Distribution")
        st.plotly_chart(plot_droplet_size_distribution(e), use_container_width=True)

        c1, c2, c3 = st.columns(3)
        c1.write(f"**d10** = {e.d10*1e6:.2f} µm")
        c1.write(f"**d32** = {e.d32*1e6:.2f} µm")
        c2.write(f"**d50** = {e.d50*1e6:.2f} µm")
        c2.write(f"**d90** = {e.d90*1e6:.2f} µm")
        c3.write(f"**d43** = {e.d43*1e6:.2f} µm")
        c3.write(f"**Span** = {e.span:.3f}")
        st.write(f"Converged: {'✅' if e.converged else '❌'}")

    with tab3:
        st.subheader("Level 2: Gelation — Pore Structure")
        st.plotly_chart(plot_phase_field(g), use_container_width=True)

        c1, c2, c3 = st.columns(3)
        c1.write(f"**Pore size** = {g.pore_size_mean*1e9:.1f} nm")
        c2.write(f"**Porosity** = {g.porosity:.3f}")
        c3.write(f"**Gelation α** = {g.alpha_final:.4f}")

        if g.phi_field.ndim == 2:
            st.write(f"Grid: {g.phi_field.shape[0]}×{g.phi_field.shape[1]}, "
                     f"spacing = {g.grid_spacing*1e9:.1f} nm, "
                     f"φ range: [{g.phi_field.min():.4f}, {g.phi_field.max():.4f}]")

    with tab4:
        st.subheader("Level 3: Crosslinking Kinetics")
        st.plotly_chart(plot_crosslinking_kinetics(x), use_container_width=True)

        c1, c2, c3 = st.columns(3)
        c1.write(f"**Crosslink fraction** = {x.p_final:.4f}")
        c1.write(f"**G_crosslinked** = {x.G_chitosan_final:.0f} Pa")
        c2.write(f"**Mesh size ξ** = {x.xi_final*1e9:.1f} nm")
        c2.write(f"**M_c** = {x.Mc_final:.0f} g/mol")
        c3.write(f"**ν_e** = {x.nu_e_final:.2e} /m³")

    with tab5:
        st.subheader("Level 4: Mechanical Properties")

        left, right = st.columns(2)
        with left:
            st.plotly_chart(plot_hertz_contact(m), use_container_width=True)
        with right:
            st.plotly_chart(plot_kav_curve(m), use_container_width=True)

        c1, c2, c3 = st.columns(3)
        c1.write(f"**G_agarose** = {m.G_agarose:.0f} Pa ({m.G_agarose/1000:.1f} kPa)")
        c2.write(f"**G_crosslinked** = {m.G_chitosan:.0f} Pa ({m.G_chitosan/1000:.1f} kPa)")
        c3.write(f"**G_DN** = {m.G_DN:.0f} Pa ({m.G_DN/1000:.1f} kPa)")

    # ── Optimization Assessment ───────────────────────────────────────

    st.divider()
    st.header("🎯 Optimization Assessment")

    # Compute objectives inline using sidebar targets instead of hardcoded ones
    d32_dev_obj = abs(e.d32 * 1e6 - target_d32) / target_d32
    pore_dev_obj = abs(g.pore_size_mean * 1e9 - target_pore) / target_pore
    G_dev_obj = abs(np.log10(max(m.G_DN, 1)) - np.log10(target_G * 1000))
    obj = np.array([d32_dev_obj, pore_dev_obj, G_dev_obj])

    st.write("**Objective Values** (lower = closer to target):")
    oc1, oc2, oc3 = st.columns(3)
    oc1.metric("f_1 (d32 deviation)", f"{obj[0]:.3f}",
               help=f"|d32 - {target_d32} um| / {target_d32} um")
    oc2.metric("f_2 (pore deviation)", f"{obj[1]:.3f}",
               help=f"|pore - {target_pore} nm| / {target_pore} nm")
    oc3.metric("f_3 (modulus deviation)", f"{obj[2]:.3f}",
               help=f"|log10(G_DN) - log10({target_G*1000})|")

    overall = np.mean(obj)
    if overall < 0.3:
        st.success(f"🟢 **Excellent match** (avg. deviation = {overall:.3f}). Parameters are near-optimal.")
    elif overall < 1.0:
        st.warning(f"🟡 **Moderate match** (avg. deviation = {overall:.3f}). Consider optimization.")
    else:
        st.error(f"🔴 **Poor match** (avg. deviation = {overall:.3f}). Significant parameter adjustment needed.")

    # Recommendations
    st.subheader("Recommendations")
    recs = []
    if obj[0] > 0.5:
        if e.d32 * 1e6 > target_d32:
            recs.append(f"**Increase RPM** (currently {rpm}) — d32 is {e.d32*1e6:.1f} µm vs target {target_d32} µm. "
                       f"Try {rpm * 1.5:.0f}+ RPM or increase Span-80 concentration.")
        else:
            recs.append(f"**Decrease RPM** — d32 is {e.d32*1e6:.1f} µm, smaller than target {target_d32} µm.")

    if obj[1] > 0.5:
        recs.append(f"**Adjust cooling rate** — pore size ({g.pore_size_mean*1e9:.0f} nm) deviates from "
                    f"target ({target_pore} nm). Slower cooling → larger pores; faster → finer.")

    if obj[2] > 0.5:
        if m.G_DN / 1000 < target_G:
            recs.append(f"**Increase crosslinker concentration** or crosslinking time — G_DN ({m.G_DN/1000:.1f} kPa) below "
                       f"target ({target_G} kPa). Current crosslinker is stoichiometry-limited at p={x.p_final:.1%}.")
        else:
            recs.append(f"**Reduce polymer concentration** — G_DN ({m.G_DN/1000:.1f} kPa) exceeds target ({target_G} kPa).")

    if not recs:
        st.write("All objectives are within acceptable range. No adjustments needed.")
    else:
        for i, rec in enumerate(recs, 1):
            st.write(f"{i}. {rec}")

    # ── Trust Assessment ───────────────────────────────────────────────

    if "trust" in st.session_state:
        st.divider()
        st.header("Trust Assessment")

        trust = st.session_state["trust"]

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

    # ── Calibration Protocol Link ─────────────────────────────────────

    st.divider()
    st.header("📋 Calibration & Validation")
    cal_path = Path(__file__).resolve().parents[3] / "docs" / "04_calibration_protocol.md"
    if cal_path.exists():
        with st.expander("View Calibration Wet-Lab Protocol"):
            st.markdown(cal_path.read_text(encoding="utf-8"))
    st.info(
        "The simulation uses literature-estimated constants that should be calibrated "
        "against your specific materials. See **docs/04_calibration_protocol.md** for "
        "a 5-study, ~30-preparation wet-lab protocol covering: interfacial tension "
        "(K_L, Γ∞), chitosan viscosity (η_intr), breakage dynamics (C3), "
        "pore structure (empirical coefficients), and IPN mechanics (η_coupling)."
    )

# ─── Footer ───────────────────────────────────────────────────────────────

st.divider()
st.caption(
    "EmulSim v0.1.0 — Multi-scale simulation for double-network polysaccharide hydrogel microspheres. "
    "Pipeline: PBE (Alopaeus) → Cahn-Hilliard → Crosslinking kinetics → Rubber Elasticity + Ogston."
)
