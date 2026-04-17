"""Explicit Streamlit navigation (v9.0, milestone M8).

Replaces Streamlit's default auto-listing of every file under pages/
with an explicit navigation configuration that:

    1. Exposes the main app at "/" (default page).
    2. Registers reagent_detail.py as a hidden page reachable via
       query-parameter links (/reagent_detail?key=...).
    3. Hides the auto-sidebar page list so "app" and "reagent detail"
       no longer appear in the top-left corner.

Requires streamlit ≥ 1.36 (Option A from the architect plan).
Pre-flight (M0) confirmed streamlit 1.55.0 is available.

Status: STUB (populated in M8).
"""

from __future__ import annotations

from typing import Any


def build_navigation() -> Any:
    """Build the st.navigation object for the v9.0 UI.

    Returns a streamlit.Navigation instance. Call .run() on it from app.py.
    """
    raise NotImplementedError("Populated in milestone M8")
