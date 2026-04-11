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

# P1-1: Pipeline Scope selector
pipeline_scope = st.sidebar.radio(
    "Pipeline Scope",
    ["M1: Fabrication Only", "M1+M2: + Functionalization", "M1+M2+M3: Full Pipeline"],
    index=0,
    help="M1: microsphere fabrication. M1+M2: adds surface functionalization. M1+M2+M3: adds performance simulation."
)

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
    c_span80_vol_pct = st.sidebar.slider("Span-80 in oil (% v/v)", 0.2, 5.0, 1.5, step=0.1)
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
    cooling_rate_Cmin = st.sidebar.slider("Cooling Rate (C/min)", 0.1, 15.0, 0.67, step=0.1,
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

c_genipin_mM = st.sidebar.slider("Crosslinker Concentration (mM)", 0.5, 500.0, 10.0, step=0.5)

# ── Stoichiometry guidance ──
from emulsim.level3_crosslinking.solver import (
    available_amine_concentration,
    recommended_crosslinker_concentration,
)
# P1-2: DDA slider (replaces hardcoded 0.90)
_DDA = st.sidebar.slider("DDA (degree of deacetylation)", 0.50, 0.99, 0.85, step=0.01,
                          help="Chitosan degree of deacetylation. Affects NH2 density and crosslinking capacity.")
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

# P1-3: Validate M1 inputs and display blockers/warnings in sidebar
from emulsim.visualization.ui_validators import validate_m1_inputs as _validate_m1
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
        st.sidebar.error(f"BLOCKER: {_blk}")
if _m1_val.warnings:
    for _wrn in _m1_val.warnings:
        st.sidebar.warning(_wrn)

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

# ─── P2-1: Module 2 Sidebar Section ──────────────────────────────────────────

m2_steps = []
if "M2" in pipeline_scope:
    from emulsim.module2_functionalization.reagent_profiles import REAGENT_PROFILES as _REAGENT_PROFILES
    from emulsim.module2_functionalization import ModificationStep, ModificationStepType, ACSSiteType

    st.sidebar.divider()
    st.sidebar.header("Module 2: Functionalization")

    n_m2_steps = st.sidebar.number_input("Modification Steps", 0, 3, 1, key="m2_n_steps")
    for i in range(int(n_m2_steps)):
        with st.sidebar.expander(f"Step {i + 1}", expanded=(i == 0)):
            step_type = st.selectbox(
                "Chemistry",
                ["Secondary Crosslinking", "Hydroxyl Activation"],
                key=f"m2_type_{i}",
            )
            if "Crosslinking" in step_type:
                _reagent_options = {
                    "Genipin": "genipin_secondary",
                    "Glutaraldehyde": "glutaraldehyde_secondary",
                }
            else:
                _reagent_options = {
                    "ECH (Epichlorohydrin)": "ech_activation",
                    "DVS (Divinyl Sulfone)": "dvs_activation",
                }
            _reagent_label = st.selectbox("Reagent", list(_reagent_options.keys()), key=f"m2_reagent_{i}")
            _reagent_key = _reagent_options[_reagent_label]

            _profile = _REAGENT_PROFILES[_reagent_key]
            st.caption(f"k={_profile.k_forward:.1e} | E_a={_profile.E_a / 1000:.0f} kJ/mol")

            _conc = st.number_input("Concentration (mM)", 0.5, 200.0, 10.0, key=f"m2_conc_{i}")
            _temp_C = st.slider(
                "Temperature (C)",
                4, 80,
                int(_profile.temperature_default - 273.15),
                key=f"m2_temp_{i}",
            )
            _time_h = st.number_input(
                "Time (h)", 0.25, 48.0, float(_profile.time_default / 3600), key=f"m2_time_{i}"
            )
            _ph = st.slider("pH", 3.0, 14.0, float(_profile.ph_optimum), step=0.5, key=f"m2_ph_{i}")

            _target_acs = ACSSiteType.AMINE_PRIMARY if "Crosslinking" in step_type else ACSSiteType.HYDROXYL
            _step_type_enum = (
                ModificationStepType.SECONDARY_CROSSLINKING
                if "Crosslinking" in step_type
                else ModificationStepType.ACTIVATION
            )
            m2_steps.append(ModificationStep(
                step_type=_step_type_enum,
                reagent_key=_reagent_key,
                target_acs=_target_acs,
                temperature=_temp_C + 273.15,
                time=_time_h * 3600,
                ph=_ph,
                reagent_concentration=_conc,
            ))

    st.sidebar.caption(
        "Planned (not yet implemented): Ligand Coupling, Protein Coupling, Quenching"
    )

# ─── P3-5: Module 3 Sidebar Section ──────────────────────────────────────────

m3_app_mode = "Chromatography"
chrom_mode = "Breakthrough"
col_diam_mm = 10.0
bed_height_cm = 10.0
bed_porosity = 0.38
flow_rate_mL = 1.0
C_feed_mg = 1.0
feed_dur_min = 10.0
total_time_min = 20.0
q_max = 100.0
K_L = 1000.0
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

if "M3" in pipeline_scope:
    st.sidebar.divider()
    st.sidebar.header("Module 3: Performance")

    m3_app_mode = st.sidebar.radio("Application", ["Chromatography", "Catalysis"], key="m3_mode")

    # Column geometry (shared between modes)
    st.sidebar.subheader("Column/Reactor")
    col_diam_mm = st.sidebar.number_input("Column I.D. (mm)", 1.0, 50.0, 10.0, key="m3_col_d")
    bed_height_cm = st.sidebar.number_input("Bed height (cm)", 1.0, 30.0, 10.0, key="m3_bed_h")
    bed_porosity = st.sidebar.slider("Bed porosity", 0.25, 0.50, 0.38, step=0.01, key="m3_eps_b")
    flow_rate_mL = st.sidebar.number_input("Flow rate (mL/min)", 0.01, 20.0, 1.0, step=0.1, key="m3_flow")

    # Real-time pressure preview (if M2 result available)
    if "m2_result" in st.session_state:
        m2r = st.session_state["m2_result"]
        from emulsim.module3_performance.hydrodynamics import ColumnGeometry as _CG_preview
        _preview_col = _CG_preview(
            diameter=col_diam_mm / 1000, bed_height=bed_height_cm / 100,
            particle_diameter=m2r.m1_contract.bead_d50,
            bed_porosity=bed_porosity,
            G_DN=m2r.G_DN_updated, E_star=m2r.E_star_updated)
        _dP = _preview_col.pressure_drop(flow_rate_mL / 60e6)
        st.sidebar.metric("Pressure drop", f"{_dP / 1000:.1f} kPa")

    if m3_app_mode == "Chromatography":
        # Phase 4: Chromatography mode selector
        chrom_mode = st.sidebar.radio("Mode", ["Breakthrough", "Gradient Elution"], key="m3_chrom_mode")

        st.sidebar.subheader("Feed")
        C_feed_mg = st.sidebar.number_input("Feed conc. (mg/mL)", 0.01, 50.0, 1.0, key="m3_Cfeed")
        feed_dur_min = st.sidebar.number_input("Feed duration (min)", 1.0, 60.0, 10.0, key="m3_feed_dur")
        total_time_min = st.sidebar.number_input("Total time (min)", 5.0, 120.0, 20.0, key="m3_total_t")

        # Isotherm
        st.sidebar.subheader("Isotherm")
        q_max = st.sidebar.number_input("q_max (mol/m\u00b3)", 1.0, 500.0, 100.0, key="m3_qmax")
        K_L = st.sidebar.number_input("K_L (m\u00b3/mol)", 10.0, 1e5, 1000.0, key="m3_KL")
        st.sidebar.caption("Default isotherm parameters are illustrative — user calibration required.")

        # Detection
        ext_coeff = st.sidebar.number_input(
            "\u03b5\u2082\u2088\u2080 (1/(M\u00b7cm))", 1000.0, 200000.0, 36000.0, key="m3_ext"
        )

        # Phase 4: Gradient elution controls
        if chrom_mode == "Gradient Elution":
            st.sidebar.subheader("Gradient")
            grad_start = st.sidebar.number_input("Start (mM)", 0.0, 1000.0, 0.0, key="m3_grad_start")
            grad_end = st.sidebar.number_input("End (mM)", 0.0, 1000.0, 500.0, key="m3_grad_end")
            grad_dur_min = st.sidebar.number_input("Duration (min)", 1.0, 60.0, 10.0, key="m3_grad_dur")

    elif m3_app_mode == "Catalysis":
        # Phase 5: Catalysis UI
        st.sidebar.subheader("Enzyme Kinetics")
        V_max = st.sidebar.number_input("V_max (mol/(m\u00b3\u00b7s))", 0.001, 100.0, 1.0, key="m3_Vmax")
        K_m = st.sidebar.number_input("K_m (mM)", 0.01, 100.0, 1.0, key="m3_Km")
        S_feed = st.sidebar.number_input("Substrate feed (mM)", 0.1, 100.0, 10.0, key="m3_Sfeed")
        D_eff = st.sidebar.number_input("D_eff (m\u00b2/s)", 1e-12, 1e-8, 1e-10, format="%.1e", key="m3_Deff")
        k_deact = st.sidebar.number_input("k_d (1/s)", 0.0, 1e-3, 0.0, format="%.1e", key="m3_kd")
        cat_time_h = st.sidebar.number_input("Sim. time (h)", 0.1, 48.0, 1.0, key="m3_cat_t")

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

        progress.progress(40, text=f"M1 complete in {elapsed:.1f}s")

        # P2-2: Module 2 execution
        if "M2" in pipeline_scope and m2_steps:
            progress.progress(50, text="Module 2: Functionalization...")
            from emulsim.pipeline.orchestrator import export_for_module2
            from emulsim.module2_functionalization import ModificationOrchestrator

            db_trust = PropertyDatabase()
            props_trust_m2 = db_trust.update_for_conditions(
                params.formulation.T_oil, params.formulation.c_agarose,
                params.formulation.c_chitosan, params.formulation.c_span80,
            )
            for _k, _v in _custom_props_overrides.items():
                if hasattr(props_trust_m2, _k):
                    setattr(props_trust_m2, _k, _v)
            _trust_m2 = assess_trust(result, params, props_trust_m2,
                                     crosslinker_key=_xl_sel_key,
                                     l2_mode=l2_mode)
            try:
                _contract = export_for_module2(
                    result, _trust_m2,
                    crosslinker_key=_xl_sel_key,
                    props=props_trust_m2,
                )
                _m2_orch = ModificationOrchestrator()
                _m2_result = _m2_orch.run(_contract, m2_steps)
                st.session_state["m2_result"] = _m2_result
                st.session_state["m1_contract"] = _contract
            except Exception as _m2_ex:
                st.warning(f"Module 2 failed: {_m2_ex}")
                st.session_state.pop("m2_result", None)
        else:
            st.session_state.pop("m2_result", None)

        # ── Phase 6: Module 3 execution ──────────────────────────────────
        if "M3" in pipeline_scope and "m2_result" in st.session_state:
            progress.progress(75, text="Module 3: Performance simulation...")
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
                    LangmuirIsotherm,
                    run_breakthrough,
                    CompetitiveLangmuirIsotherm,
                    run_gradient_elution,
                )
                from emulsim.module3_performance.gradient import make_linear_gradient
                from emulsim.visualization.ui_validators import validate_m3_chromatography, validate_m3_result
                MW_Da = 66500.0  # default BSA
                C_feed_mol = C_feed_mg / (MW_Da * 1e-3)  # mg/mL → mol/m³

                _is_grad = (chrom_mode == "Gradient Elution")
                _col_val = validate_m3_chromatography(
                    flow_rate=flow_rate_mL / 60e6,
                    column=column,
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
                            isotherm = LangmuirIsotherm(q_max=q_max, K_L=K_L)
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
                                K_L=np.array([K_L]),
                                component_names=["Protein"],
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
                        st.warning(f"Module 3 chromatography failed: {_m3_ex}")

            elif m3_app_mode == "Catalysis":
                from emulsim.module3_performance.catalysis.packed_bed import solve_packed_bed
                # K_m and S_feed input as mM → convert to mol/m³
                K_m_mol = K_m  # mM = mol/m³ (since 1 mM = 1e-3 mol/L = 1 mol/m³)
                S_feed_mol = S_feed  # same unit
                try:
                    cat = solve_packed_bed(
                        bed_length=bed_height_cm / 100,
                        bed_diameter=col_diam_mm / 1000,
                        particle_diameter=m2r.m1_contract.bead_d50,
                        bed_porosity=bed_porosity,
                        particle_porosity=m2r.m1_contract.porosity,
                        V_max=V_max,
                        K_m=K_m_mol,
                        S_feed=S_feed_mol,
                        flow_rate=flow_rate_mL / 60e6,
                        D_eff=D_eff,
                        k_deact=k_deact,
                        total_time=cat_time_h * 3600,
                    )
                    st.session_state["m3_result_cat"] = cat
                except Exception as _m3_ex:
                    st.warning(f"Module 3 catalysis failed: {_m3_ex}")
        elif "M3" in pipeline_scope and "m2_result" not in st.session_state:
            st.warning("Module 3 requires Module 2 results. Enable M1+M2 scope and run with M2 steps defined.")

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

    _tab_labels = [
        "📈 Dashboard", "🫧 L1: Emulsification", "🧊 L2: Gelation",
        "🔗 L3: Crosslinking", "💪 L4: Mechanical",
    ]
    _show_m2_tab = "M2" in pipeline_scope and "m2_result" in st.session_state
    if _show_m2_tab:
        _tab_labels.append("🧪 M2: Functionalization")

    _show_m3_bt_tab = "m3_result_bt" in st.session_state
    _show_m3_ge_tab = "m3_result_ge" in st.session_state
    _show_m3_cat_tab = "m3_result_cat" in st.session_state
    if _show_m3_bt_tab:
        _tab_labels.append("📊 M3: Breakthrough")
    if _show_m3_ge_tab:
        _tab_labels.append("📊 M3: Gradient Elution")
    if _show_m3_cat_tab:
        _tab_labels.append("⚗️ M3: Catalysis")

    _all_tabs = st.tabs(_tab_labels)
    tab1, tab2, tab3, tab4, tab5 = _all_tabs[:5]
    _tab_idx = 5
    tab_m2 = _all_tabs[_tab_idx] if _show_m2_tab else None
    if _show_m2_tab:
        _tab_idx += 1
    tab_m3_bt = _all_tabs[_tab_idx] if _show_m3_bt_tab else None
    if _show_m3_bt_tab:
        _tab_idx += 1
    tab_m3_ge = _all_tabs[_tab_idx] if _show_m3_ge_tab else None
    if _show_m3_ge_tab:
        _tab_idx += 1
    tab_m3_cat = _all_tabs[_tab_idx] if _show_m3_cat_tab else None

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

    # ── P2-3: Module 2 Results Tab ────────────────────────────────────

    if _show_m2_tab and tab_m2 is not None:
        with tab_m2:
            from emulsim.visualization.plots_m2 import plot_acs_waterfall, plot_surface_area_comparison
            _m2 = st.session_state["m2_result"]
            st.subheader("Module 2: Surface Functionalization")

            # KPI row
            _mc1, _mc2, _mc3 = st.columns(3)
            _mc1.metric("Steps Executed", len(_m2.modification_history))
            _mc2.metric("G_DN Updated", f"{_m2.G_DN_updated / 1000:.1f} kPa")
            _mc3.metric("E* Updated", f"{_m2.E_star_updated / 1000:.1f} kPa")

            st.divider()

            # ACS site summary
            st.markdown("**ACS Site Inventory (remaining after all steps)**")
            for _site_type, _profile in _m2.acs_profiles.items():
                st.write(
                    f"**{_site_type.value}**: remaining = {_profile.remaining_sites:.2e} mol/particle"
                )

            st.divider()

            # Modification history
            st.markdown("**Modification History**")
            for _i, _mr in enumerate(_m2.modification_history):
                _tgt = _mr.step.target_acs
                _sites_before = _mr.acs_before.get(_tgt)
                _consumed_str = ""
                if _sites_before is not None:
                    _consumed = _mr.conversion * _sites_before.remaining_sites
                    _consumed_str = f" | sites consumed: {_consumed:.2e} mol/particle"
                st.write(
                    f"Step {_i + 1}: {_mr.step.reagent_key} — "
                    f"conversion {_mr.conversion:.1%}{_consumed_str}"
                )

            st.divider()

            # ACS waterfall chart
            _fig_wf = plot_acs_waterfall(_m2.modification_history)
            if _fig_wf is not None:
                st.plotly_chart(_fig_wf, use_container_width=True)

            # Surface area breakdown
            _fig_sa = plot_surface_area_comparison(_m2.surface_model)
            if _fig_sa is not None:
                st.plotly_chart(_fig_sa, use_container_width=True)

            st.caption(
                "[!] semi_quantitative — ACS inventory uses simplified site-density model. "
                "Default rate constants are illustrative — user calibration required."
            )

    # ── Phase 6: Module 3 Results Tabs ───────────────────────────────────

    from emulsim.visualization.plots_m3 import (
        plot_chromatogram,
        plot_breakthrough_curve,
        plot_peak_table,
        plot_michaelis_menten,
        plot_effectiveness_factor,
        plot_activity_decay,
        plot_conversion_vs_time,
        plot_pressure_flow_curve,
    )

    if _show_m3_bt_tab and tab_m3_bt is not None:
        with tab_m3_bt:
            from emulsim.visualization.ui_validators import validate_m3_result as _val_m3_res
            _bt = st.session_state["m3_result_bt"]
            st.subheader("Module 3: Breakthrough Chromatography")

            # Validation badges
            _bt_val = _val_m3_res(
                mass_balance_error=_bt.mass_balance_error,
                pressure_drop=_bt.pressure_drop,
            )
            _mb_pct = abs(_bt.mass_balance_error) * 100.0
            if _mb_pct > 5.0:
                st.error(f"Mass balance error = {_mb_pct:.1f}% — results numerically unreliable.")
            elif _mb_pct > 2.0:
                st.warning(f"Mass balance error = {_mb_pct:.1f}% — treat with caution.")
            else:
                st.success(f"Mass balance error = {_mb_pct:.2f}% — acceptable.")
            if _bt_val.blockers:
                for _blk in _bt_val.blockers:
                    st.error(f"Blocker: {_blk}")

            # DBC KPI row
            _dbc_c1, _dbc_c2, _dbc_c3, _dbc_c4 = st.columns(4)
            _dbc_c1.metric("DBC\u2085%", f"{_bt.dbc_5pct:.1f} mol/m\u00b3" if not np.isnan(_bt.dbc_5pct) else "N/A")
            _dbc_c2.metric("DBC\u2081\u2080%", f"{_bt.dbc_10pct:.1f} mol/m\u00b3" if not np.isnan(_bt.dbc_10pct) else "N/A")
            _dbc_c3.metric("DBC\u2085\u2080%", f"{_bt.dbc_50pct:.1f} mol/m\u00b3" if not np.isnan(_bt.dbc_50pct) else "N/A")
            _dbc_c4.metric("Pressure drop", f"{_bt.pressure_drop / 1000:.1f} kPa")

            st.divider()

            # Breakthrough curve
            _MW_Da_bt = 66500.0
            _C_feed_mol_bt = C_feed_mg / (_MW_Da_bt * 1e-3)
            st.plotly_chart(
                plot_breakthrough_curve(
                    time=_bt.time,
                    C_outlet=_bt.C_outlet,
                    C_feed=_C_feed_mol_bt,
                    dbc_5=_bt.dbc_5pct,
                    dbc_10=_bt.dbc_10pct,
                    dbc_50=_bt.dbc_50pct,
                ),
                use_container_width=True,
            )

            # UV chromatogram
            st.plotly_chart(
                plot_chromatogram(time=_bt.time, uv_signal=_bt.uv_signal),
                use_container_width=True,
            )

            # Pressure-flow curve (requires ColumnGeometry from session state)
            if "m2_result" in st.session_state:
                from emulsim.module3_performance.hydrodynamics import ColumnGeometry as _CG_bt
                _m2r_bt = st.session_state["m2_result"]
                _col_bt = _CG_bt(
                    diameter=col_diam_mm / 1000,
                    bed_height=bed_height_cm / 100,
                    particle_diameter=_m2r_bt.m1_contract.bead_d50,
                    bed_porosity=bed_porosity,
                    particle_porosity=_m2r_bt.m1_contract.porosity,
                    G_DN=_m2r_bt.G_DN_updated,
                    E_star=_m2r_bt.E_star_updated,
                )
                st.plotly_chart(plot_pressure_flow_curve(_col_bt), use_container_width=True)

            st.caption(
                "Mechanistic prediction: Lumped Rate Model (LRM) with Langmuir isotherm. "
                "DBC values are model-based — calibrate isotherm parameters with batch uptake experiments."
            )

    if _show_m3_ge_tab and tab_m3_ge is not None:
        with tab_m3_ge:
            _ge = st.session_state["m3_result_ge"]
            st.subheader("Module 3: Gradient Elution Chromatography")

            # Gradient affects binding badge
            _grad_affects = getattr(_ge, "gradient_affects_binding", False)
            if _grad_affects:
                st.success("Gradient affects binding: YES — competitive Langmuir mechanistically coupled.")
            else:
                st.info("Gradient affects binding: NO — gradient is diagnostic/display only (BF-2 not deployed).")

            # Chromatogram with gradient overlay
            _ge_time = getattr(_ge, "time", None)
            _ge_uv = getattr(_ge, "uv_signal", None)
            _ge_grad = getattr(_ge, "gradient_values", None)

            if _ge_time is not None and _ge_uv is not None:
                st.plotly_chart(
                    plot_chromatogram(
                        time=_ge_time,
                        uv_signal=_ge_uv,
                        gradient_values=_ge_grad,
                        gradient_affects_binding=_grad_affects,
                    ),
                    use_container_width=True,
                )

            # Peak table
            _ge_peaks = getattr(_ge, "peaks", [])
            st.plotly_chart(plot_peak_table(_ge_peaks), use_container_width=True)

            st.caption(
                "Ranking only — multi-component gradient elution with competitive Langmuir. "
                "Gradient salt concentration does not yet affect binding affinity (BF-2 pending). "
                "Elution order is indicative; quantitative yields require isotherm calibration."
            )

    if _show_m3_cat_tab and tab_m3_cat is not None:
        with tab_m3_cat:
            _cat = st.session_state["m3_result_cat"]
            st.subheader("Module 3: Packed-Bed Catalytic Reactor")

            # KPI row
            _cat_c1, _cat_c2, _cat_c3, _cat_c4 = st.columns(4)
            _cat_c1.metric("Final Conversion", f"{_cat.conversion:.1%}")
            _cat_c2.metric("Effectiveness Factor \u03b7", f"{_cat.effectiveness_factor:.3f}")
            _cat_c3.metric("Thiele Modulus \u03a6", f"{_cat.thiele_modulus:.2f}")
            _cat_c4.metric("Productivity", f"{_cat.productivity:.3e} mol/(m\u00b3\u00b7s)")

            _mb_cat_pct = abs(_cat.mass_balance_error) * 100.0
            if _mb_cat_pct > 5.0:
                st.error(f"Catalysis mass balance error = {_mb_cat_pct:.1f}% — results unreliable.")
            elif _mb_cat_pct > 2.0:
                st.warning(f"Catalysis mass balance error = {_mb_cat_pct:.1f}% — treat with caution.")
            else:
                st.success(f"Catalysis mass balance error = {_mb_cat_pct:.2f}% — acceptable.")

            st.divider()

            # Conversion vs time
            _time_h_cat = _cat.time / 3600.0
            _S_in = float(_cat.S_outlet[0]) if len(_cat.S_outlet) > 0 else 1.0
            if _S_in <= 0:
                _S_in = S_feed  # fallback to sidebar value (in mol/m³ = mM)
            _conversion_arr = np.clip(1.0 - _cat.S_outlet / max(_S_in, 1e-12), 0.0, 1.0)
            st.plotly_chart(
                plot_conversion_vs_time(time_hours=_time_h_cat, conversion=_conversion_arr),
                use_container_width=True,
            )

            # Michaelis-Menten curve
            _S_range_mm = np.linspace(0.01, max(S_feed * 3, 10.0), 200)
            st.plotly_chart(
                plot_michaelis_menten(
                    S_range=_S_range_mm,
                    V_max=V_max,
                    K_m=K_m,
                    eta=_cat.effectiveness_factor,
                ),
                use_container_width=True,
            )

            # Effectiveness factor vs Thiele modulus
            _phi_range = np.logspace(-2, 2, 200)
            st.plotly_chart(
                plot_effectiveness_factor(_phi_range),
                use_container_width=True,
            )

            # Activity decay
            if hasattr(_cat, "activity_history") and len(_cat.activity_history) > 0:
                st.plotly_chart(
                    plot_activity_decay(
                        time_hours=_time_h_cat,
                        activity=_cat.activity_history,
                    ),
                    use_container_width=True,
                )

            # Thiele regime diagnostic
            _phi_val = _cat.thiele_modulus
            if _phi_val < 0.3:
                st.success(f"Thiele modulus \u03a6 = {_phi_val:.2f} — kinetic regime (no diffusion limitation).")
            elif _phi_val < 3.0:
                st.warning(f"Thiele modulus \u03a6 = {_phi_val:.2f} — transition regime (partial diffusion limitation).")
            else:
                st.error(f"Thiele modulus \u03a6 = {_phi_val:.2f} — diffusion-limited regime. "
                         "Effectiveness factor \u03b7 is significantly < 1. Consider smaller particles or higher D_eff.")

            st.caption(
                "Mechanistic prediction: Transient PFR with Michaelis-Menten kinetics, "
                "Thiele modulus effectiveness factor, and first-order deactivation. "
                "K_m and V_max require calibration against your specific enzyme."
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
    "EmulSim v4.0 -- Modular process-modeling platform combining calibrated empirical models "
    "and mechanistic submodels for emulsified hydrogel microsphere fabrication. "
    "L1: PBE emulsification (adaptive convergence) | L2: Empirical pore or Cahn-Hilliard 2D | "
    "L3: Chemistry-specific crosslinking (per-chemistry eta) | "
    "L4: Phenomenological + Flory-Rehner affine IPN + Hashin-Shtrikman bounds | "
    "M2: Surface functionalization (secondary crosslinking, hydroxyl activation) | "
    "M3: Chromatography (breakthrough + gradient elution) + Packed-bed catalysis."
)
