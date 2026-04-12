"""Calibration framework for EmulSim v6.0.

Enables users to supply measured resin characterization data that
overrides semi-quantitative defaults in FunctionalMediaContract.
"""

from .calibration_data import CalibrationEntry
from .calibration_store import CalibrationStore

__all__ = ["CalibrationEntry", "CalibrationStore"]
