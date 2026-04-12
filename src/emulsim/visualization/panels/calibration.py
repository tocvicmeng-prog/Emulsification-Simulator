"""Calibration import/inspect panel for EmulSim UI.

v6.0-rc: Must-have per audit F13. Provides:
  - JSON file upload for CalibrationEntry data
  - Inspection table showing all loaded entries
  - Active override display when calibration is applied to FMC
"""

from __future__ import annotations

import json
import tempfile
from pathlib import Path

import streamlit as st

from emulsim.calibration import CalibrationEntry, CalibrationStore


def render_calibration_panel() -> CalibrationStore | None:
    """Render the calibration import/inspect panel.

    Placed in the sidebar or as an expander in the main area.
    Returns the CalibrationStore if entries are loaded, else None.

    Returns:
        CalibrationStore with loaded entries, or None if empty.
    """
    with st.expander("\U0001f4c1 Calibration Data", expanded=False):
        st.caption(
            "Import user-measured calibration data (JSON) to override "
            "default isotherm parameters. See docs/04_calibration_protocol.md."
        )

        # ── Initialize store in session state ──
        if "_cal_store" not in st.session_state:
            st.session_state["_cal_store"] = CalibrationStore()
        store: CalibrationStore = st.session_state["_cal_store"]

        # ── JSON Upload ──
        uploaded = st.file_uploader(
            "Upload calibration JSON",
            type=["json"],
            key="cal_json_upload",
            help="JSON array of CalibrationEntry objects. "
                 "See calibration/ package for schema.",
        )

        if uploaded is not None:
            try:
                data = json.loads(uploaded.read().decode("utf-8"))
                if not isinstance(data, list):
                    st.error("JSON must be an array of calibration entries.")
                else:
                    # Write to temp file and load via CalibrationStore
                    with tempfile.NamedTemporaryFile(
                        mode="w", suffix=".json", delete=False
                    ) as tmp:
                        json.dump(data, tmp)
                        tmp_path = tmp.name
                    new_store = CalibrationStore()
                    n_loaded = new_store.load_json(tmp_path)
                    st.session_state["_cal_store"] = new_store
                    store = new_store
                    st.success(f"Loaded {n_loaded} calibration entries.")
                    Path(tmp_path).unlink(missing_ok=True)
            except json.JSONDecodeError as e:
                st.error(f"Invalid JSON: {e}")
            except Exception as e:
                st.error(f"Failed to load calibration: {e}")

        # ── Manual Entry ──
        with st.popover("\u2795 Add Entry Manually"):
            _profile = st.text_input("Profile key", "protein_a_coupling",
                                      key="cal_man_profile")
            _param = st.text_input("Parameter", "estimated_q_max",
                                    key="cal_man_param")
            _value = st.number_input("Measured value", 0.0, 1e6, 100.0,
                                      key="cal_man_value")
            _units = st.text_input("Units", "mol/m3", key="cal_man_units")
            _molecule = st.text_input("Target molecule", "", key="cal_man_mol")
            _conf = st.selectbox("Confidence",
                                  ["measured", "literature", "estimated"],
                                  key="cal_man_conf")
            _mtype = st.selectbox("Measurement type",
                                   ["", "static_binding", "DBC10", "DBC5", "batch_uptake"],
                                   key="cal_man_mtype")
            _source = st.text_input("Source reference", "", key="cal_man_src")

            if st.button("Add", key="cal_man_add"):
                entry = CalibrationEntry(
                    profile_key=_profile,
                    parameter_name=_param,
                    measured_value=_value,
                    units=_units,
                    target_molecule=_molecule,
                    confidence=_conf,
                    measurement_type=_mtype,
                    source_reference=_source,
                )
                store.add(entry)
                st.success(f"Added: {_profile}.{_param} = {_value} {_units}")

        # ── Display Loaded Entries ──
        if len(store) > 0:
            st.markdown(f"**{len(store)} calibration entries loaded**")
            for i, entry in enumerate(store.entries):
                _conf_icon = {
                    "measured": "\U0001f7e2",
                    "literature": "\U0001f7e1",
                    "estimated": "\U0001f534",
                }.get(entry.confidence, "\u26aa")
                st.write(
                    f"{_conf_icon} **{entry.profile_key}**.{entry.parameter_name} "
                    f"= {entry.measured_value} {entry.units} "
                    f"({entry.confidence})"
                )
                if entry.target_molecule:
                    st.caption(f"  Molecule: {entry.target_molecule} | "
                               f"T={entry.temperature_C}\u00b0C | pH={entry.ph}")
        else:
            st.info("No calibration data loaded. Upload JSON or add entries manually.")

        # ── Clear Button ──
        if len(store) > 0:
            if st.button("Clear All Calibration Data", key="cal_clear"):
                st.session_state["_cal_store"] = CalibrationStore()
                st.rerun()

    return store if len(store) > 0 else None
