"""EmulSim: Multi-scale emulsification simulation for hydrogel microsphere preparation."""

__version__ = "9.2.1"

from .datatypes import (
    SimulationParameters,
    MaterialProperties,
    EmulsificationResult,
    GelationTimingResult,
    GelationResult,
    CrosslinkingResult,
    MechanicalResult,
    FullResult,
)
from .pipeline.orchestrator import PipelineOrchestrator
from .properties.database import PropertyDatabase

__all__ = [
    "SimulationParameters",
    "MaterialProperties",
    "EmulsificationResult",
    "GelationTimingResult",
    "GelationResult",
    "CrosslinkingResult",
    "MechanicalResult",
    "FullResult",
    "PipelineOrchestrator",
    "PropertyDatabase",
    "run_pipeline",
    "__version__",
]


def run_pipeline(params: SimulationParameters | None = None, **kwargs) -> FullResult:
    """Run the full L1-L4 simulation pipeline.

    Convenience function for quick usage:
        from emulsim import run_pipeline
        result = run_pipeline()  # uses defaults
        result = run_pipeline(params)  # custom params
    """
    if params is None:
        params = SimulationParameters()
    orch = PipelineOrchestrator()
    return orch.run_single(params, **kwargs)
