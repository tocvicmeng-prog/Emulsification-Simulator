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

# Force reload of ALL core modules to avoid Streamlit's import cache.
# Order matters: reload dependencies before dependents.
import importlib
import emulsim.datatypes as _dt_mod; importlib.reload(_dt_mod)
import emulsim.level1_emulsification.thermal as _th_mod; importlib.reload(_th_mod)
import emulsim.level1_emulsification.energy as _en_mod; importlib.reload(_en_mod)
import emulsim.level1_emulsification.kernels as _kr_mod; importlib.reload(_kr_mod)
import emulsim.level1_emulsification.validation as _vl_mod; importlib.reload(_vl_mod)
import emulsim.level1_emulsification.solver as _sv_mod; importlib.reload(_sv_mod)
import emulsim.properties.database as _db_mod; importlib.reload(_db_mod)
import emulsim.pipeline.orchestrator as _orch_mod; importlib.reload(_orch_mod)
import emulsim.trust as _trust_mod; importlib.reload(_trust_mod)
import emulsim.level3_crosslinking.solver as _xl_mod; importlib.reload(_xl_mod)

from emulsim.datatypes import (
    SimulationParameters, MaterialProperties, FormulationParameters,
    EmulsificationParameters, MixerGeometry, SolverSettings,
    VesselGeometry, StirrerGeometry, HeatingConfig, KernelConfig,
    VesselType, StirrerType, HeatingStrategy, ModelMode,
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

# ─── Page Config ──────────────────────────────────────────────────────────

st.set_page_config(
    page_title="EmulSim — Microsphere Simulation",
    page_icon="🔬",
    layout="wide",
)

st.title("🔬 EmulSim — Double-Network Hydrogel Microsphere Simulator")
st.caption("Multi-scale simulation: Emulsification → Gelation → Crosslinking → Mechanical Properties")

# ─── Sidebar: Parameter Input ─────────────────────────────────────────────

st.sidebar.header("Process Parameters")

# ── Mode Selection ──
sim_mode = st.sidebar.radio(
    "Hardware Mode",
    ["Rotor-Stator (Legacy)", "Stirred Vessel (New)"],
    index=0,
    help="Legacy: 25 mm rotor-stator at 3k-25k RPM. "
         "Stirred Vessel: pitched-blade or small rotor-stator at 800-9000 RPM.",
)
is_stirred = "Stirred" in sim_mode

# ── Scientific Mode ──
model_mode = st.sidebar.radio(
    "Scientific Mode",
    ["Empirical Engineering", "Hybrid Coupled", "Mechanistic Research"],
    index=1,  # default: Hybrid Coupled
    help="Empirical: fast screening, suppresses model warnings. "
         "Hybrid: default, phenomenological DN model with trust warnings. "
         "Mechanistic: Flory-Rehner affine IPN model, strictest trust gates.",
)
_mode_map = {
    "Empirical Engineering": ModelMode.EMPIRICAL_ENGINEERING,
    "Hybrid Coupled": ModelMode.HYBRID_COUPLED,
    "Mechanistic Research": ModelMode.MECHANISTIC_RESEARCH,
}
model_mode_enum = _mode_map[model_mode]

st.sidebar.subheader("Emulsification (L1)")

if is_stirred:
    # ── Stirred-vessel controls ──
    vessel_choice = st.sidebar.selectbox(
        "Vessel", ["Glass Beaker (100 mm)", "Jacketed Vessel (92 mm)"],
        help="Beaker: flat-plate heating. Jacketed: hot-water circulation.",
    )
    stirrer_choice = st.sidebar.selectbox(
        "Stirrer", ["Stirrer A - Pitched Blade (59 mm)", "Stirrer B - Rotor-Stator (32 mm)"],
    )
    is_stirrer_A = "Pitched" in stirrer_choice

    if is_stirrer_A:
        rpm = st.sidebar.slider("Stirrer Speed (RPM)", 800, 2000, 1300, step=50)
    else:
        rpm = st.sidebar.slider("Stirrer Speed (RPM)", 800, 9000, 1800, step=100)

    t_emul = st.sidebar.number_input("Emulsification Time (min)", 1, 60, 10)

    # Volumes
    st.sidebar.caption("Working liquid volumes")
    v_oil_mL = st.sidebar.slider("Oil + Span-80 (mL)", 100, 500, 300, step=10)
    v_poly_mL = st.sidebar.slider("Polysaccharide solution (mL)", 50, 300, 200, step=10)
    total_mL = v_oil_mL + v_poly_mL
    phi_d = v_poly_mL / total_mL
    st.sidebar.caption(f"Total: {total_mL} mL | phi_d = {phi_d:.2f}")

    # Heating
    if "Beaker" in vessel_choice:
        heating_choice = "Flat Plate"
        st.sidebar.caption("Heating: flat-plate (150C -> 80C oil)")
    else:
        heating_choice = "Hot Water Jacket"
        st.sidebar.caption("Heating: jacket (85C circulating water)")

else:
    # ── Legacy rotor-stator controls ──
    rpm = st.sidebar.slider("Rotor Speed (RPM)", 3000, 25000, 10000, step=500)
    t_emul = st.sidebar.number_input("Emulsification Time (min)", 1, 60, 10)
    phi_d = st.sidebar.slider("Dispersed Phase Fraction (phi_d)", 0.01, 0.30, 0.05, step=0.01)
    v_oil_mL = 300
    v_poly_mL = 200
    total_mL = 500

if model_mode_enum != ModelMode.EMPIRICAL_ENGINEERING:
    with st.sidebar.expander("Advanced PBE Settings",
                              expanded=(model_mode_enum == ModelMode.MECHANISTIC_RESEARCH)):
        l1_t_max = st.slider("Max emulsification time (s)", 60, 1800, 600, step=60,
                              help="Absolute ceiling for adaptive time extensions")
        l1_conv_tol = st.slider("Convergence tolerance", 0.005, 0.10, 0.01, step=0.005,
                                 format="%.3f",
                                 help="Relative d32 variation threshold for steady state")
        l1_max_ext = st.number_input("Max extensions", 0, 5, 2,
                                      help="Number of half-interval extensions if PBE not converged")
else:
    l1_t_max = 600.0
    l1_conv_tol = 0.01
    l1_max_ext = 2

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

if is_stirred:
    c_span80_vol_pct = st.sidebar.slider("Span-80 in oil (% v/v)", 0.5, 5.0, 1.5, step=0.1)
    c_span80_pct = c_span80_vol_pct * 986.0 / 1000.0  # convert v/v to approx w/v
    T_oil_C = st.sidebar.slider("Paraffin Oil Temperature (C)", 65, 110, 80, step=1,
                                 help="Oil temperature at start of emulsification. "
                                      "Flat-plate heater can reach ~80 C; "
                                      "higher temps may need oil bath.")
else:
    c_span80_pct = st.sidebar.slider("Surfactant (% w/v)", 0.5, 5.0, 2.0, step=0.1)
    c_span80_vol_pct = 1.5
    T_oil_C = st.sidebar.slider("Oil Temperature (C)", 60, 95, 90)

st.sidebar.subheader("Cooling & Gelation (L2)")
if is_stirred:
    cooling_rate_Cmin = st.sidebar.slider("Cooling Rate (C/min)", 0.3, 3.0, 0.67, step=0.1,
                                           help="Natural cooling: ~0.67 C/min for 500 mL beaker")
else:
    cooling_rate_Cmin = st.sidebar.slider("Cooling Rate (C/min)", 1.0, 20.0, 10.0, step=0.5)
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

# ── Stoichiometry guidance ──
from emulsim.level3_crosslinking.solver import (
    available_amine_concentration,
    recommended_crosslinker_concentration,
)
_DDA = 0.90  # TODO: make user-adjustable in future
_M_GlcN = 161.16
_c_chit_kg = c_chitosan_pct * 10.0  # % -> kg/m3
_NH2 = available_amine_concentration(_c_chit_kg, _DDA, _M_GlcN)
if _NH2 > 0:
    _ratio = c_genipin_mM / _NH2
    _c_rec = recommended_crosslinker_concentration(_c_chit_kg, _DDA, _M_GlcN, target_p=0.20)
    if _ratio >= 0.10:
        st.sidebar.success(f"Crosslinker/NH\u2082 = {_ratio:.3f} \u2014 sufficient for p \u2265 0.20")
    elif _ratio >= 0.05:
        st.sidebar.warning(
            f"Crosslinker/NH\u2082 = {_ratio:.3f} \u2014 may be limiting. "
            f"Recommend \u2265 {_c_rec:.1f} mM for target p = 0.20"
        )
    else:
        st.sidebar.error(
            f"Crosslinker/NH\u2082 = {_ratio:.3f} \u2014 severely limiting! "
            f"Increase to at least {_c_rec:.1f} mM for target p = 0.20"
        )

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
if is_stirred:
    target_d_mode = st.sidebar.number_input("Target d_mode (um)", 10.0, 500.0, 100.0, step=10.0,
                                             help="Modal diameter of microspheres")
    target_d32 = target_d_mode  # for compatibility
    target_pore = st.sidebar.number_input("Target Pore Size (nm)", 10, 500, 100, step=10)
else:
    target_d32 = st.sidebar.number_input("Target d32 (um)", 0.5, 50.0, 2.0, step=0.5)
    target_d_mode = target_d32
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

# U5: Per-chemistry eta display — use crosslinker-specific recommended value
_xl_eta = xl.eta_coupling_recommended
st.sidebar.caption(f"Per-chemistry η for {xl.name}: {_xl_eta:+.2f}")
custom_eta_coup = _const_input(
    "η coupling (IPN)", "eta_coup",
    _xl_eta, "-", f"{xl.name} library", -0.5, 0.5, 0.05,
    "eta_coup", "%.2f",
)
# Detect whether user chose "Custom" for eta coupling (radio key is "src_eta_coup")
eta_coup_source = st.session_state.get("src_eta_coup", "Literature")

# ─── Build Parameters ─────────────────────────────────────────────────────

if is_stirred:
    # Build stirred-vessel equipment
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
    # Sync heating T_initial with user-selected oil temperature
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

# Apply custom/literature material constants
from emulsim.datatypes import MaterialProperties as _MP
_custom_props_overrides = {
    "breakage_C3": custom_C3,
    "eta_coupling": custom_eta_coup,
    "eta_is_custom": (eta_coup_source == "Custom"),
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
                                                crosslinker_key=_xl_sel_key,
                                                l2_mode=l2_mode)

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

    _d_mode = getattr(e, 'd_mode', 0.0)

    if is_stirred:
        d_primary_dev = abs(_d_mode * 1e6 - target_d_mode) / target_d_mode * 100
    else:
        d_primary_dev = abs(e.d32 * 1e6 - target_d32) / target_d32 * 100
    pore_dev = abs(g.pore_size_mean * 1e9 - target_pore) / target_pore * 100
    G_dev = abs(m.G_DN / 1000 - target_G) / target_G * 100

    if is_stirred:
        col1.metric("d_mode", f"{_d_mode*1e6:.1f} um",
                    delta=f"{d_primary_dev:.0f}% from target",
                    delta_color="inverse")
    else:
        col1.metric("d32", f"{e.d32*1e6:.2f} um",
                    delta=f"{d_primary_dev:.0f}% from target",
                    delta_color="inverse")
    col2.metric("Pore Size", f"{g.pore_size_mean*1e9:.1f} nm",
                delta=f"{pore_dev:.0f}% from target",
                delta_color="inverse")
    # U6: Show Hashin-Shtrikman bounds on G_DN metric if available
    _hs_lo = getattr(m, 'G_DN_lower', 0.0)
    _hs_hi = getattr(m, 'G_DN_upper', 0.0)
    if _hs_lo > 0 and _hs_hi > 0:
        col3.metric("G_DN", f"{m.G_DN/1000:.1f} kPa",
                    delta=f"Ref: [{_hs_lo/1000:.1f}, {_hs_hi/1000:.1f}] kPa (single-phase)",
                    delta_color="off")
    else:
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

        _conv_icon = "✅ Yes" if e.converged else "⚠️ No (still evolving)"
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

        # U7: Model used label + Flory-Rehner status
        _model_label = getattr(m, 'model_used', 'phenomenological')
        st.caption(f"Model: {_model_label}")
        if _model_label == "flory_rehner_affine":
            st.success("Flory-Rehner affine IPN model converged — mechanistic crosslink density used.")
        elif _model_label == "phenomenological":
            if model_mode_enum == ModelMode.MECHANISTIC_RESEARCH:
                st.info("Using phenomenological model. Switch to Mechanistic Research mode to enable Flory-Rehner affine IPN.")
            else:
                st.info("Phenomenological DN model (G1 + G2 + η√(G1·G2)). Switch to Mechanistic Research for affine IPN.")

        left, right = st.columns(2)
        with left:
            st.plotly_chart(plot_hertz_contact(m), use_container_width=True)
        with right:
            st.plotly_chart(plot_kav_curve(m), use_container_width=True)

        # U6c: Modulus comparison bar chart
        st.plotly_chart(plot_modulus_comparison(m), use_container_width=True)

        c1, c2, c3 = st.columns(3)
        c1.write(f"**G_agarose** = {m.G_agarose:.0f} Pa ({m.G_agarose/1000:.1f} kPa)")
        c2.write(f"**G_crosslinked** = {m.G_chitosan:.0f} Pa ({m.G_chitosan/1000:.1f} kPa)")
        c3.write(f"**G_DN** = {m.G_DN:.0f} Pa ({m.G_DN/1000:.1f} kPa)")

        # U6b: HS bounds display
        if _hs_lo > 0 and _hs_hi > 0:
            st.write(f"**Single-phase reference (HS composite bounds, not IPN bounds):** [{_hs_lo/1000:.1f}, {_hs_hi/1000:.1f}] kPa")
        st.caption(f"Model: {_model_label}")

        # U7: Model comparison in Mechanistic Research mode
        if model_mode_enum == ModelMode.MECHANISTIC_RESEARCH:
            st.write("**Model Comparison:**")
            from emulsim.level4_mechanical.solver import double_network_modulus as _dnm
            _eta_comp = getattr(x.network_metadata, 'eta_coupling_recommended', -0.15) if x.network_metadata else -0.15
            _G_pheno = _dnm(m.G_agarose, m.G_chitosan, _eta_comp)
            st.write(
                f"Phenomenological: {_G_pheno/1000:.1f} kPa | "
                f"{'Affine IPN' if _model_label == 'flory_rehner_affine' else _model_label}: {m.G_DN/1000:.1f} kPa"
            )

    # ── Optimization Assessment ───────────────────────────────────────

    st.divider()
    st.header("🎯 Optimization Assessment")

    # Compute objectives inline using sidebar targets
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
        st.success(f"🟢 **Excellent match** (avg. deviation = {overall:.3f}). Parameters are near-optimal.")
    elif overall < 1.0:
        st.warning(f"🟡 **Moderate match** (avg. deviation = {overall:.3f}). Consider optimization.")
    else:
        st.error(f"🔴 **Poor match** (avg. deviation = {overall:.3f}). Significant parameter adjustment needed.")

    # Recommendations
    st.subheader("Recommendations")
    recs = []
    if obj[0] > 0.5:
        if d_obj_val > d_obj_target:
            recs.append(f"**Increase RPM** (currently {rpm}) -- droplet size is {d_obj_val:.1f} um vs target {d_obj_target:.0f} um. "
                       f"Try higher RPM or increase Span-80 concentration.")
        else:
            recs.append(f"**Decrease RPM** -- droplet size is {d_obj_val:.1f} um, smaller than target {d_obj_target:.0f} um.")

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

        _mode_desc = {
            "empirical_engineering": "Empirical Engineering — trust warnings relaxed for screening",
            "hybrid_coupled": "Hybrid Coupled — phenomenological models with trust warnings",
            "mechanistic_research": "Mechanistic Research — strictest gates, Flory-Rehner when available",
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
    "EmulSim v3.1 -- Modular process-modeling platform combining calibrated empirical models "
    "and mechanistic submodels for emulsified hydrogel microsphere fabrication. "
    "L1: PBE emulsification (adaptive convergence) | L2: Empirical pore or Cahn-Hilliard 2D | "
    "L3: Chemistry-specific crosslinking (per-chemistry eta) | "
    "L4: Phenomenological + Flory-Rehner affine IPN + Hashin-Shtrikman bounds."
)
