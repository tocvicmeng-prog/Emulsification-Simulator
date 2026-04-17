"""Hardware section (v9.0, milestone M2).

Relocates the Hardware Mode selector from the Global Settings sidebar
into the M1 Emulsification section. Scientific rationale (per SA §C):
hardware choice (rotor-stator vs stirred vessel) gates the L1 PBE solver
only; L2/L3/L4 and M2/M3 do not depend on it. Placing it in Global
Settings falsely implied cross-module coupling.

M2 scope: render the Hardware Mode radio only. The vessel/stirrer/heating/
RPM widgets remain inside tab_m1.py for now; they will move into
HardwareContext via render_hardware_section in a later milestone.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import streamlit as st


@dataclass
class HardwareContext:
    """Hardware selections that feed EmulsificationParameters.

    Populated progressively: M2 sets only ``is_stirred``; later
    milestones add the full vessel/stirrer/heating/kernels chain.
    """

    is_stirred: bool
    vessel: Any = None
    stirrer: Any = None
    mixer: Any = None
    heating: Any = None
    kernels: Any = None
    rpm: float = 0.0


def render_hardware_mode_radio(*, key: str = "m1v9_hardware_mode") -> bool:
    """Render the Hardware Mode radio at the top of the M1 Emulsification section.

    Returns True for Stirred Vessel, False for Rotor-Stator (Legacy).
    Default: Rotor-Stator (Legacy) — matches v8.x behaviour.
    """
    sim_mode = st.radio(
        "Hardware Mode",
        ["Rotor-Stator (Legacy)", "Stirred Vessel (New)"],
        index=0,
        horizontal=True,
        help="Legacy: 25 mm rotor-stator at 3k-25k RPM. "
             "Stirred Vessel: pitched-blade or small rotor-stator at 800-9000 RPM. "
             "Hardware choice affects L1 PBE only; L2/L3/L4 and M2/M3 are hardware-agnostic.",
        key=key,
    )
    return "Stirred" in sim_mode


def render_hardware_section(*, t_emul_min: int = 10) -> HardwareContext:
    """Full hardware section — populated in M4 when tab_m1 is fully decomposed."""
    raise NotImplementedError("Populated in milestone M4; M2 only ships render_hardware_mode_radio")
