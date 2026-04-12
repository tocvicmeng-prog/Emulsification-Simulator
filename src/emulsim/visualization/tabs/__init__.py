"""Tab rendering modules for EmulSim Streamlit UI.

v6.0: Extracted from monolithic app.py into separate modules.
Each module provides a render_tab_*() function called from app.py.
"""

from .tab_m1 import render_tab_m1
from .tab_m2 import render_tab_m2
from .tab_m3 import render_tab_m3

__all__ = ["render_tab_m1", "render_tab_m2", "render_tab_m3"]
