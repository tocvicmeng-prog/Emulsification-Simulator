"""Module 3: Chromatographic Performance — breakthrough simulation.

Phase C provides minimal single-component Langmuir breakthrough:
  - ColumnGeometry: Kozeny-Carman pressure drop + compressibility.
  - LangmuirIsotherm: Single-component Langmuir equilibrium.
  - solve_lrm: Lumped Rate Model PDE solver (FV + BDF).
  - run_breakthrough: Full breakthrough orchestration with UV detection.

Phase E adds multi-component and gradient capabilities:
  - CompetitiveLangmuirIsotherm: Multi-component competitive Langmuir.
  - SMAIsotherm: Steric Mass Action for ion exchange.
  - IMACCompetitionIsotherm: IMAC with imidazole competition.
  - ProteinAIsotherm: pH-dependent affinity (Protein A / IgG).
  - GradientProgram: Piecewise-linear gradient generator.
  - make_linear_gradient, make_step_gradient, make_wash_gradient: Convenience.
  - compute_fluorescence_signal: Fluorescence detection.
  - compute_conductivity: Conductivity from salt concentration.
  - simulate_esi_charge_envelope: ESI-MS charge envelope simulation.
  - compute_tic: Total ion chromatogram.
  - run_gradient_elution: Multi-component gradient elution orchestration.
  - GradientElutionResult, PeakInfo: Result containers.
"""

# Phase C exports
from .hydrodynamics import ColumnGeometry
from .isotherms.langmuir import LangmuirIsotherm
from .transport.lumped_rate import LRMResult, solve_lrm
from .detection.uv import compute_uv_signal, apply_detector_broadening
from .orchestrator import BreakthroughResult, run_breakthrough

# Phase E: Isotherms
from .isotherms.competitive_langmuir import CompetitiveLangmuirIsotherm
from .isotherms.sma import SMAIsotherm
from .isotherms.imac import IMACCompetitionIsotherm
from .isotherms.protein_a import ProteinAIsotherm

# Phase E: Gradient
from .gradient import (
    GradientProgram,
    make_linear_gradient,
    make_step_gradient,
    make_wash_gradient,
)

# Phase E: Detection
from .detection.fluorescence import compute_fluorescence_signal, fluorescence_detection_limit
from .detection.conductivity import (
    compute_conductivity,
    conductivity_to_ms_per_cm,
    conductivity_to_nacl_concentration,
    compute_chromatogram_conductivity,
)
from .detection.ms import (
    ESISpectrum,
    simulate_esi_charge_envelope,
    compute_tic,
    compute_extracted_ion_chromatogram,
)

# Phase E: Gradient elution orchestrator
from .orchestrator import (
    PeakInfo,
    GradientElutionResult,
    run_gradient_elution,
)

__all__ = [
    # Phase C
    "ColumnGeometry",
    "LangmuirIsotherm",
    "LRMResult",
    "solve_lrm",
    "BreakthroughResult",
    "run_breakthrough",
    "compute_uv_signal",
    "apply_detector_broadening",
    # Phase E isotherms
    "CompetitiveLangmuirIsotherm",
    "SMAIsotherm",
    "IMACCompetitionIsotherm",
    "ProteinAIsotherm",
    # Phase E gradient
    "GradientProgram",
    "make_linear_gradient",
    "make_step_gradient",
    "make_wash_gradient",
    # Phase E detection
    "compute_fluorescence_signal",
    "fluorescence_detection_limit",
    "compute_conductivity",
    "conductivity_to_ms_per_cm",
    "conductivity_to_nacl_concentration",
    "compute_chromatogram_conductivity",
    "ESISpectrum",
    "simulate_esi_charge_envelope",
    "compute_tic",
    "compute_extracted_ion_chromatogram",
    # Phase E orchestrator
    "PeakInfo",
    "GradientElutionResult",
    "run_gradient_elution",
]
