"""Streamlit UI for the emulsification simulation system.

Launch: streamlit run src/emulsim/visualization/app.py

v6.0 — Modular tab layout with Pipeline Scope selector.
Tab rendering extracted to tabs/tab_m1.py, tab_m2.py, tab_m3.py.
This file is the thin orchestrator: page config, sidebar, tab routing, footer.
"""

from __future__ import annotations

import logging
import sys
from pathlib import Path

import streamlit as st

logger = logging.getLogger(__name__)

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

from emulsim.datatypes import ModelMode
from emulsim.visualization.ui_state import SessionStateManager
from emulsim.visualization.tabs import render_tab_m1, render_tab_m2, render_tab_m3
from emulsim.visualization.panels import (
    render_calibration_panel,
    render_uncertainty_panel,
    render_lifetime_panel,
)

# ─── Page Config ──────────────────────────────────────────────────────────

st.set_page_config(
    page_title="EmulSim — Microsphere Simulation",
    page_icon="\U0001f52c",
    layout="wide",
)

_title_col, _manual_col = st.columns([6, 1])
with _title_col:
    st.title("\U0001f52c EmulSim — Double-Network Hydrogel Microsphere Simulator")
    st.caption("Multi-scale simulation: Emulsification \u2192 Gelation \u2192 Crosslinking \u2192 Mechanical Properties")
with _manual_col:
    # Upper-right download button for the First Edition PDF manual.
    _manual_pdf = (
        _root
        / "docs" / "user_manual"
        / "polysaccharide_microsphere_simulator_first_edition.pdf"
    )
    if _manual_pdf.exists():
        with open(_manual_pdf, "rb") as _f:
            st.download_button(
                label="Manual (PDF)",
                data=_f.read(),
                file_name="EmulSim_First_Edition.pdf",
                mime="application/pdf",
                help="Polysaccharide-Based Microsphere Emulsification "
                     "Simulator — First Edition (instruction manual).",
                use_container_width=True,
            )
    else:
        st.caption(
            "Manual PDF not built yet — run "
            "`python docs/user_manual/build_pdf.py`"
        )

# ─── Session State Manager ────────────────────────────────────────────────

if "_state_mgr" not in st.session_state:
    st.session_state["_state_mgr"] = SessionStateManager()
    st.session_state["_state_mgr"].bind_store(st.session_state)
_smgr: SessionStateManager = st.session_state["_state_mgr"]

# ─── Sidebar: Global Settings Only ───────────────────────────────────────

st.sidebar.header("Global Settings")

sim_mode = st.sidebar.radio(
    "Hardware Mode",
    ["Rotor-Stator (Legacy)", "Stirred Vessel (New)"],
    index=0,
    help="Legacy: 25 mm rotor-stator at 3k-25k RPM. "
         "Stirred Vessel: pitched-blade or small rotor-stator at 800-9000 RPM.",
)
is_stirred = "Stirred" in sim_mode

model_mode = st.sidebar.radio(
    "Scientific Mode",
    ["Empirical Engineering", "Hybrid Coupled", "Mechanistic Research"],
    index=1,
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

# ─── Sidebar: v6.0 Framework Panels ──────────────────────────────────────

with st.sidebar:
    st.divider()
    st.subheader("v6.0 Frameworks")
    _cal_store = render_calibration_panel()
    _unc_contract = render_uncertainty_panel()
    _lt_proj = render_lifetime_panel()

# ─── Pipeline Scope (main area, above tabs) ──────────────────────────────

pipeline_scope = st.radio(
    "Pipeline Scope",
    ["M1: Fabrication Only", "M1\u2192M2: + Functionalization", "M1\u2192M2\u2192M3: Full Pipeline"],
    index=0,
    horizontal=True,
    key="pipeline_scope_selector",
    help="M1: microsphere fabrication. M1\u2192M2: adds surface functionalization. "
         "M1\u2192M2\u2192M3: adds performance simulation.",
)

# ─── Build Dynamic Tab List ──────────────────────────────────────────────

_module_tab_labels = ["\U0001f52c M1: Fabrication"]
if "M2" in pipeline_scope:
    _module_tab_labels.append("\U0001f9ea M2: Functionalization")
if "M3" in pipeline_scope:
    _module_tab_labels.append("\U0001f4ca M3: Performance")

_module_tabs = st.tabs(_module_tab_labels)
_tab_m1 = _module_tabs[0]
_tab_m2 = _module_tabs[1] if "M2" in pipeline_scope else None
_tab_m3 = _module_tabs[2] if "M3" in pipeline_scope and len(_module_tabs) > 2 else None

# ─── Calibration Protocol Loader (shared) ────────────────────────────────

_proto_path = Path(__file__).resolve().parents[3] / "docs" / "04_calibration_protocol.md"
_proto_sections: dict = {}
if _proto_path.exists():
    try:
        _proto_text = _proto_path.read_text(encoding="utf-8")
        _study_map = {
            "K_L":        ("Study 1 -- Interfacial Tension", "## 2. Study 1", "## 3. Study 2"),
            "Gamma_inf":  ("Study 1 -- Interfacial Tension", "## 2. Study 1", "## 3. Study 2"),
            "eta_chit":   ("Study 2 -- Chitosan Viscosity",  "## 3. Study 2", "## 4. Study 3"),
            "C3":         ("Study 3 -- Viscous Breakage",     "## 4. Study 3", "## 5. Study 4"),
            "pore":       ("Study 4 -- Pore Structure",       "## 5. Study 4", "## 6. Study 5"),
            "eta_coup":   ("Study 5 -- IPN Coupling",         "## 6. Study 5", "## 7. Inputting"),
        }
        # Marker discipline: when a section header in docs/04 is renamed, the
        # original silent-fallback behaviour (serve [start:] when end is missing,
        # skip entirely when start is missing) hid the drift. Preserve the
        # tolerant fallback for backward compatibility, but log per-key drift
        # so an operator notices stale markers in the console output.
        for key, (title, start_marker, end_marker) in _study_map.items():
            s_idx = _proto_text.find(start_marker)
            e_idx = _proto_text.find(end_marker)
            if s_idx >= 0:
                section = _proto_text[s_idx:e_idx] if e_idx > s_idx else _proto_text[s_idx:]
                _proto_sections[key] = (title, section)
                if e_idx <= s_idx:
                    logger.warning(
                        "Calibration protocol end marker %r not found for %s; "
                        "serving from start marker to end of document.",
                        end_marker, key,
                    )
            else:
                logger.warning(
                    "Calibration protocol start marker %r not found for %s; "
                    "popover disabled.",
                    start_marker, key,
                )
    except Exception as exc:  # pragma: no cover — defensive: file was read once at import
        logger.warning("Could not parse calibration protocol %s: %s", _proto_path, exc)
        _proto_sections = {}


def _const_input(container, label, key, lit_val, unit, source_short, lo, hi, step,
                  proto_key, fmt="%.3f"):
    """Render per-constant selector with protocol link inside a given container."""
    if proto_key in _proto_sections:
        title, section_md = _proto_sections[proto_key]
        with container.popover(f"\U0001f4cb {label}"):
            st.markdown(f"### Calibration Protocol: {title}")
            st.markdown(section_md[:3000])
            if len(section_md) > 3000:
                st.caption("... (see full protocol in docs/04_calibration_protocol.md)")
    else:
        container.markdown(f"**{label}**")

    src = container.radio(
        "Source",
        ["Literature", "Custom"],
        index=0, horizontal=True, key=f"src_{key}",
        label_visibility="collapsed",
    )
    if src == "Literature":
        container.caption(f"  = {lit_val} {unit}  ({source_short})")
        return lit_val
    else:
        return container.number_input(
            f"{label} ({unit})", lo, hi, float(lit_val), step=step,
            format=fmt, key=f"val_{key}",
        )


# ═════════════════════════════════════════════════════════════════════════
# TAB RENDERING
# ═════════════════════════════════════════════════════════════════════════

render_tab_m1(
    tab_container=_tab_m1,
    is_stirred=is_stirred,
    model_mode_enum=model_mode_enum,
    _smgr=_smgr,
    _const_input=_const_input,
    _proto_sections=_proto_sections,
)

if _tab_m2 is not None:
    render_tab_m2(tab_container=_tab_m2, _smgr=_smgr)

if _tab_m3 is not None:
    render_tab_m3(tab_container=_tab_m3)


# ─── Footer ───────────────────────────────────────────────────────────────

st.divider()
st.caption(
    "EmulSim v6.0 \u2014 Modular tab layout. "
    "L1: PBE emulsification (adaptive convergence) | L2: Empirical pore or Cahn-Hilliard 2D | "
    "L3: Chemistry-specific crosslinking (per-chemistry eta) | "
    "L4: Phenomenological + Flory-Rehner affine IPN + Hashin-Shtrikman bounds | "
    "M2: Surface functionalization (secondary crosslinking, hydroxyl activation) | "
    "M3: Chromatography (breakthrough + gradient elution) + Packed-bed catalysis."
)
