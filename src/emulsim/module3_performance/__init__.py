"""Module 3: Chromatographic Performance — breakthrough simulation.

Phase C provides minimal single-component Langmuir breakthrough:
  - ColumnGeometry: Kozeny-Carman pressure drop + compressibility.
  - LangmuirIsotherm: Single-component Langmuir equilibrium.
  - solve_lrm: Lumped Rate Model PDE solver (FV + BDF).
  - run_breakthrough: Full breakthrough orchestration with UV detection.
"""

from .hydrodynamics import ColumnGeometry
from .isotherms.langmuir import LangmuirIsotherm
from .transport.lumped_rate import LRMResult, solve_lrm
from .detection.uv import compute_uv_signal, apply_detector_broadening
from .orchestrator import BreakthroughResult, run_breakthrough

__all__ = [
    "ColumnGeometry",
    "LangmuirIsotherm",
    "LRMResult",
    "solve_lrm",
    "BreakthroughResult",
    "run_breakthrough",
    "compute_uv_signal",
    "apply_detector_broadening",
]
