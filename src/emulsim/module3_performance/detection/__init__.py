"""Detection models for chromatographic signal processing."""

from .uv import compute_uv_signal, apply_detector_broadening

__all__ = ["compute_uv_signal", "apply_detector_broadening"]
