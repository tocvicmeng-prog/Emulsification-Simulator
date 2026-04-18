"""v9.0 M9 regression: TOML configs still load and run through the orchestrator.

Baseline: fast_smoke.toml produces d32 = 22.08 µm on v8.3.5. Any divergence
here means the v9.0 UI redesign has leaked into the backend pipeline.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from emulsim.config import load_config
from emulsim.pipeline.orchestrator import PipelineOrchestrator
from emulsim.properties.database import PropertyDatabase


_REPO = Path(__file__).resolve().parents[1]


@pytest.mark.slow
@pytest.mark.parametrize(
    "config_name,expected_d32_um_range",
    [
        ("fast_smoke.toml", (20.0, 25.0)),
        ("default.toml", (1.0, 50.0)),
        ("stirred_vessel.toml", (50.0, 300.0)),
    ],
)
def test_toml_config_loads_and_runs(config_name, expected_d32_um_range):
    """Each shipped TOML config loads, runs through orch.run_single,
    and produces a d32 in the family-appropriate range."""
    cfg_path = _REPO / "configs" / config_name
    if not cfg_path.exists():
        pytest.skip(f"{config_name} not present")
    params = load_config(str(cfg_path))
    orch = PipelineOrchestrator(db=PropertyDatabase())
    result = orch.run_single(params)
    d32_um = result.emulsification.d32 * 1e6
    lo, hi = expected_d32_um_range
    assert lo <= d32_um <= hi, (
        f"{config_name}: d32 = {d32_um:.2f} µm outside expected range "
        f"[{lo}, {hi}]. v9.0 UI redesign may have leaked into the backend."
    )


def test_fast_smoke_d32_preserved():
    """Tight check: fast_smoke.toml's d32 must stay at the v8.3.5 baseline."""
    cfg_path = _REPO / "configs" / "fast_smoke.toml"
    params = load_config(str(cfg_path))
    orch = PipelineOrchestrator(db=PropertyDatabase())
    result = orch.run_single(params)
    d32_um = result.emulsification.d32 * 1e6
    assert 21.5 <= d32_um <= 22.5, (
        f"fast_smoke d32 = {d32_um:.2f} µm; v8.3.5 baseline was 22.08 µm. "
        "The v9.0 UI redesign was supposed to be UI-only."
    )
