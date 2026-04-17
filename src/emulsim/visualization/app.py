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

# v9.0 M8 + design-consultation v1.0: hide Streamlit's auto-discovered pages/
# listing from the sidebar AND apply the EmulSim design system (DESIGN.md).
# Rationale (per scientific-advisor §D + design consultation): reagent_detail.py
# is meaningful only with specific query params, and the default Streamlit
# theme uses Source Sans Pro + a consumer-SaaS palette that fights the
# scientific-instrument positioning of EmulSim.
#
# IMPORTANT: use st.html (NOT st.markdown + unsafe_allow_html) for this block.
# st.markdown runs the content through a Markdown parser first, which treats
# "*" characters as italic markers. CSS attribute selectors like
# [class*="css"] get mangled: the markdown parser pairs two "*" together,
# drops both, and wraps the content between them in <em>…</em>. Result: the
# raw CSS shows up as literal text at the top of the page and none of the
# design system applies. st.html inserts the HTML verbatim, no markdown.
st.html(
    """
    <style>
        @import url('https://fonts.googleapis.com/css2?family=Geist:wght@400;500;600;700&family=Geist+Mono:wght@400;500&family=JetBrains+Mono:wght@400;500&display=swap');
    :root {
        --es-bg: #0F172A;
        --es-surface: #1E293B;
        --es-surface-2: #334155;
        --es-border: rgba(71, 85, 105, 0.5);
        --es-text: #F8FAFC;
        --es-text-muted: #94A3B8;
        --es-accent: #2DD4BF;
        --es-accent-hover: #5EEAD4;
        --es-success: #22C55E;
        --es-warning: #F59E0B;
        --es-error: #EF4444;
        --es-info: #38BDF8;
    }

    /* Hide auto-page nav (v9.0 M8) */
    [data-testid="stSidebarNav"] { display: none; }

    /* Hide the "Deploy" button from Streamlit's top-right toolbar.
     * EmulSim is a local instrument, not a Streamlit Cloud app; the
     * Deploy button was overlapping the Manual (PDF) download button
     * in the upper-right title column. Keep the running-state indicator
     * and hamburger "Rerun" menu.
     */
    [data-testid="stAppDeployButton"],
    [data-testid="stDeployButton"],
    .stAppDeployButton, .stDeployButton,
    [data-testid="stToolbar"] button[kind="headerNoPadding"]:first-child {
        display: none !important;
    }

    /* Keep a little breathing room below the top-right toolbar so the
     * Manual (PDF) button (rendered in the narrow right title column) is
     * never visually pinned under it even if a future Streamlit version
     * changes the toolbar positioning.
     */
    [data-testid="stDownloadButton"] { margin-top: 2.5rem; }

    /* ═══ Typography — Geist Sans / Mono + JetBrains Mono ══════════════ */
    /* Cascade the Geist family from the root and let descendants inherit,
     * EXCEPT icon and code elements. Do NOT use wildcard attribute
     * selectors like [class*="st-"] — they match every Streamlit wrapper
     * including the spans that carry Material Symbols icons, which makes
     * raw icon names (keyboard_double_arrow_left, expand_more, etc.) leak
     * as visible text because the icon font is no longer applied.
     */
    html, body, .stApp {
        font-family: "Geist", -apple-system, BlinkMacSystemFont, "Segoe UI",
                     system-ui, sans-serif;
        font-feature-settings: "ss01", "cv11";
    }
    /* Re-assert Material Symbols / Icons on Streamlit's icon spans, so the
     * chevrons, expanders, and toolbar icons actually render as glyphs.
     */
    .material-symbols-rounded, .material-symbols-outlined, .material-icons,
    [class*="material-symbols"], [class*="material-icons"],
    span[data-testid="stIconMaterial"] {
        font-family: "Material Symbols Rounded", "Material Symbols Outlined",
                     "Material Icons" !important;
    }
    /* Code and monospace surfaces */
    code, pre, kbd, samp, .stCode pre {
        font-family: "JetBrains Mono", "Geist Mono", Menlo, Monaco,
                     "Cascadia Code", monospace !important;
    }
    /* Numeric / data surfaces — Geist Mono with tabular-nums for decimal alignment */
    [data-testid="stMetricValue"], [data-testid="stDataFrame"] td,
    [data-testid="stTable"] td, [data-testid="stDataFrame"] th {
        font-family: "Geist Mono", "JetBrains Mono", Menlo, monospace !important;
        font-feature-settings: "tnum", "zero";
    }

    /* ═══ Spacing / density — compact, but leave room for slider bubbles */
    .block-container {
        padding-top: 2rem !important;
        padding-bottom: 3rem !important;
        max-width: none !important;
    }
    h1 { font-size: 1.75rem !important; font-weight: 700 !important;
         letter-spacing: -0.015em;
         margin-top: 0.5em !important;
         margin-bottom: 0.25em !important; }
    h2 { font-size: 1.25rem !important; font-weight: 600 !important;
         letter-spacing: -0.01em;
         margin-top: 1.75em !important;
         margin-bottom: 0.75em !important; }
    h3 { font-size: 1rem !important; font-weight: 600 !important;
         margin-top: 1.25em !important;
         margin-bottom: 0.5em !important; }

    /* Slider widgets render a value bubble ABOVE the thumb. Without
     * vertical clearance that bubble crashes into the next widget's
     * label. Give every slider container enough top room for the bubble
     * and bottom room so its own track doesn't kiss the next element.
     */
    [data-testid="stSlider"] {
        padding-top: 1.25rem !important;
        padding-bottom: 0.75rem !important;
    }

    /* Breathing room between every consecutive widget container */
    [data-testid="stElementContainer"] {
        margin-bottom: 0.25rem;
    }
    [data-testid="stVerticalBlock"] > [data-testid="stElementContainer"] + [data-testid="stElementContainer"] {
        margin-top: 0.1rem;
    }

    /* ═══ Buttons — 4px radius, teal accent, tool-like ═════════════════ */
    .stButton > button {
        border-radius: 4px !important;
        border: 1px solid var(--es-border) !important;
        font-weight: 500 !important;
        transition: background-color 150ms ease-out,
                    border-color 150ms ease-out !important;
    }
    .stButton > button[kind="primary"] {
        background-color: var(--es-accent) !important;
        color: var(--es-bg) !important;
        border-color: var(--es-accent) !important;
    }
    .stButton > button[kind="primary"]:hover {
        background-color: var(--es-accent-hover) !important;
        border-color: var(--es-accent-hover) !important;
    }

    /* ═══ Inputs — 4px radius, teal focus ring ═════════════════════════ */
    .stTextInput input, .stNumberInput input, .stSelectbox > div > div,
    .stTextArea textarea {
        border-radius: 4px !important;
        border: 1px solid var(--es-border) !important;
    }
    .stTextInput input:focus, .stNumberInput input:focus,
    .stTextArea textarea:focus {
        border-color: var(--es-accent) !important;
        box-shadow: 0 0 0 2px rgba(45, 212, 191, 0.25) !important;
    }

    /* ═══ Alerts — semantic colors matched to DESIGN.md ═══════════════ */
    [data-testid="stAlert"] { border-radius: 4px !important;
                              border-left-width: 3px !important; }

    /* ═══ Metric cards — give them some room to breathe ═══════════════ */
    [data-testid="stMetric"] {
        background: var(--es-surface);
        padding: 0.75rem 1rem !important;
        border-radius: 4px;
        border: 1px solid var(--es-border);
    }
    [data-testid="stMetricValue"] {
        font-size: 1.5rem !important;
        font-weight: 700 !important;
        letter-spacing: -0.01em;
    }
    [data-testid="stMetricLabel"] {
        font-size: 0.75rem !important;
        text-transform: uppercase;
        letter-spacing: 0.05em;
        color: var(--es-text-muted) !important;
    }

    /* ═══ Sidebar polish ══════════════════════════════════════════════ */
    [data-testid="stSidebar"] {
        background-color: var(--es-surface) !important;
        border-right: 1px solid var(--es-border);
    }
    [data-testid="stSidebar"] h2 { font-size: 0.875rem !important;
                                   text-transform: uppercase;
                                   letter-spacing: 0.06em;
                                   color: var(--es-text-muted) !important; }

    /* ═══ Data tables — tabular-nums, subtle header ════════════════════ */
    [data-testid="stDataFrame"] th {
        background-color: var(--es-surface-2) !important;
        color: var(--es-text) !important;
        font-weight: 600 !important;
        font-size: 0.75rem !important;
        text-transform: uppercase;
        letter-spacing: 0.05em;
    }

    /* ═══ Expanders — tighter, tool-like ══════════════════════════════ */
    [data-testid="stExpander"] details {
        border: 1px solid var(--es-border);
        border-radius: 4px;
    }
    [data-testid="stExpander"] summary { padding: 0.5rem 0.75rem !important; }

    /* ═══ Radio / checkbox — tighten spacing, teal accent on active ═══ */
    .stRadio > div { gap: 0.25rem !important; }
    input[type="radio"]:checked, input[type="checkbox"]:checked {
        accent-color: var(--es-accent) !important;
    }

    /* ═══ Tabs — align to the instrument aesthetic ═════════════════════ */
    [data-testid="stTabs"] button[role="tab"] {
        font-weight: 500 !important;
        padding: 0.5rem 1rem !important;
    }
    [data-testid="stTabs"] button[role="tab"][aria-selected="true"] {
        color: var(--es-accent) !important;
        border-bottom: 2px solid var(--es-accent) !important;
    }
    </style>
    """
)

_title_col, _manual_col = st.columns([6, 1])
with _title_col:
    st.title("\U0001f52c EmulSim — Double-Network Hydrogel Microsphere Simulator")
    st.caption("Multi-scale simulation: Emulsification \u2192 Gelation \u2192 Crosslinking \u2192 Mechanical Properties")
with _manual_col:
    # Upper-right download buttons for the printable PDFs. Two documents:
    #   1. Main user manual — First Edition. Principles, inputs, workflow.
    #   2. Appendix J — Functionalisation wet-lab protocols. 44 reagent
    #      protocols with full SDS-lite safety blocks for first-time users.
    _um_dir = _root / "docs" / "user_manual"
    _manual_pdf = _um_dir / "polysaccharide_microsphere_simulator_first_edition.pdf"
    _appendix_j_pdf = _um_dir / "appendix_J_functionalization_protocols.pdf"

    # Auto-build any missing target on first render. build_pdf.py iterates
    # its BUILD_TARGETS registry, so adding new PDFs only requires extending
    # that list + re-running.
    if not _manual_pdf.exists() or not _appendix_j_pdf.exists():
        _build_script = _um_dir / "build_pdf.py"
        if _build_script.exists():
            try:
                import runpy
                runpy.run_path(str(_build_script), run_name="__main__")
            except Exception as _build_ex:
                logger.warning("PDF auto-build failed: %s", _build_ex)

    if _manual_pdf.exists():
        with open(_manual_pdf, "rb") as _f:
            st.download_button(
                label="\U0001F4D8  Manual (PDF)",
                data=_f.read(),
                file_name="EmulSim_First_Edition.pdf",
                mime="application/pdf",
                help="Polysaccharide-Based Microsphere Emulsification "
                     "Simulator — First Edition (instruction manual).",
                width="stretch",
            )
    else:
        st.caption("Manual PDF missing — run `python docs/user_manual/build_pdf.py`")

    if _appendix_j_pdf.exists():
        with open(_appendix_j_pdf, "rb") as _f:
            st.download_button(
                label="\U0001F9EA  Appendix J (PDF)",
                data=_f.read(),
                file_name="EmulSim_Appendix_J_Functionalization.pdf",
                mime="application/pdf",
                help="Appendix J — Functionalisation Wet-Lab Protocols. "
                     "44 reagent-level protocols (hydroxyl activation, "
                     "ligand and protein coupling, spacer arms, IMAC metal "
                     "charging, pretreatment, washing, quenching) with "
                     "full SDS-lite safety blocks for first-time users.",
                width="stretch",
            )

# ─── Session State Manager ────────────────────────────────────────────────

if "_state_mgr" not in st.session_state:
    st.session_state["_state_mgr"] = SessionStateManager()
    st.session_state["_state_mgr"].bind_store(st.session_state)
_smgr: SessionStateManager = st.session_state["_state_mgr"]

# ─── Sidebar: Global Settings Only ───────────────────────────────────────
#
# v9.0 Family-First: Hardware Mode moved into M1 Emulsification section
# (scientific rationale in scientific-advisor audit §C — hardware gates
# the L1 PBE solver only, not downstream modules, so it belongs to M1).

st.sidebar.header("Global Settings")

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

# ─── Sidebar: Analysis Tools (as click-to-open overlays) ─────────────────
#
# Calibration data + MC uncertainty are genuinely cross-cutting (every run
# can consume them), so they stay accessible from the sidebar. Resin
# lifetime projection is M3-only, but leaving it here is acceptable until
# we complete the scope-move into the M3 tab (out of scope for this pass;
# scientific-advisor audit §E flagged it for v9.1).
#
# Design change: instead of three always-expanded expanders under the
# header "v6.0 Frameworks" (version numbers are not user-facing), each
# tool is now a popover button. Users click to open, configure, click
# anywhere else to dismiss. No permanent sidebar real estate consumed.

with st.sidebar:
    st.divider()
    st.caption("ANALYSIS TOOLS")
    with st.popover("\U0001F4C1  Calibration data", width="stretch"):
        _cal_store = render_calibration_panel()
    with st.popover("\U0001F4CA  Uncertainty (MC sampling)", width="stretch"):
        _unc_contract = render_uncertainty_panel()
    with st.popover("\u23F3  Resin lifetime projection", width="stretch"):
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
    is_stirred=None,  # v9.0: tab_m1 now renders Hardware Mode locally
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
