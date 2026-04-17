"""Tests for Node F3-b, F3-c, F4-a (v8.0 Phase 1): inverse-design BO
engine wiring + CLI + robust-BO constructor guards.

Heavy BO runs are out of scope here — this file validates:

- OptimizationEngine constructor accepts target_spec and sizes
  internal REF_POINT and penalty arrays accordingly.
- Invalid robust-BO configurations raise at construction.
- The design subcommand parser registers correctly and the
  _cmd_design handler routes a valid target spec into engine
  construction (mocked run).
"""

from __future__ import annotations

import sys
from pathlib import Path
from unittest.mock import patch, MagicMock

import numpy as np
import pytest

from emulsim.optimization.engine import OptimizationEngine
from emulsim.optimization.objectives import TargetSpec


class TestEngineTargetSpec:
    def test_default_uses_legacy_3d(self):
        engine = OptimizationEngine(
            n_initial=2, max_iterations=1,
            output_dir=Path("/tmp/_emulsim_f3b_test_legacy"),
        )
        assert engine.target_spec is None
        assert engine._n_obj == 3
        assert engine._ref_point.numel() == 3

    def test_target_spec_sizes_arrays(self):
        spec = TargetSpec(
            d32_target=2e-6, d32_tol=0.5e-6,
            pore_target=80e-9, pore_tol=10e-9,
        )
        engine = OptimizationEngine(
            n_initial=2, max_iterations=1,
            output_dir=Path("/tmp/_emulsim_f3b_test_spec"),
            target_spec=spec,
        )
        assert engine._n_obj == 2
        assert engine._ref_point.numel() == 2
        assert engine.target_spec is spec

    def test_invalid_target_spec_rejected_at_ctor(self):
        """Empty TargetSpec fails validate() at construction."""
        with pytest.raises(ValueError, match="no active dimension"):
            OptimizationEngine(
                n_initial=2, max_iterations=1,
                output_dir=Path("/tmp/_emulsim_f3b_test_empty"),
                target_spec=TargetSpec(),
            )

    def test_robust_requires_target_spec(self):
        """robust_variance_weight > 0 without target_spec raises."""
        with pytest.raises(ValueError, match="requires target_spec"):
            OptimizationEngine(
                n_initial=2, max_iterations=1,
                output_dir=Path("/tmp/_emulsim_f4_no_target"),
                robust_variance_weight=0.5,
            )

    def test_robust_needs_enough_samples(self):
        spec = TargetSpec(d32_target=2e-6, d32_tol=0.5e-6)
        with pytest.raises(ValueError, match="robust_n_samples"):
            OptimizationEngine(
                n_initial=2, max_iterations=1,
                output_dir=Path("/tmp/_emulsim_f4_few_samples"),
                target_spec=spec,
                robust_variance_weight=0.5,
                robust_n_samples=1,
            )

    def test_robust_happy_path_constructs(self):
        spec = TargetSpec(d32_target=2e-6, d32_tol=0.5e-6)
        engine = OptimizationEngine(
            n_initial=2, max_iterations=1,
            output_dir=Path("/tmp/_emulsim_f4_ok"),
            target_spec=spec,
            robust_variance_weight=0.25,
            robust_n_samples=3,
        )
        assert engine.robust_variance_weight == 0.25
        assert engine.robust_n_samples == 3


class TestDesignCLIParser:
    def test_parser_accepts_design_subcommand(self):
        """The argparse subparsers must register `design` with the key
        flags the F3-c protocol specifies."""
        # Lazy import to avoid running the CLI entry at import time
        from emulsim.__main__ import main
        with patch("sys.argv", ["emulsim", "design", "--help"]):
            with pytest.raises(SystemExit) as exc:
                main()
            # argparse exits 0 after printing help
            assert exc.value.code == 0

    def test_empty_target_spec_rejected_by_cli(self):
        """`emulsim design` with NO target flags must raise
        SystemExit with a clear error (TargetSpec.validate fails)."""
        from emulsim.__main__ import main
        with patch("sys.argv", ["emulsim", "design"]):
            with pytest.raises(SystemExit) as exc:
                main()
            # The handler raises SystemExit("Invalid --design target: ...")
            assert "Invalid" in str(exc.value) or exc.value.code != 0

    def test_design_dispatches_to_engine_with_target(self):
        """Providing --d32 + --d32-tol builds a TargetSpec and routes
        to OptimizationEngine. We mock engine.run so the test doesn't
        spin up a full BO campaign."""
        from emulsim.__main__ import main

        fake_state = MagicMock()
        fake_state.pareto_Y = np.zeros((1, 1))
        with patch(
            "emulsim.optimization.engine.OptimizationEngine.run",
            return_value=fake_state,
        ) as mock_run, patch(
            "emulsim.optimization.analysis.pareto_summary",
            return_value="ok",
        ), patch(
            "emulsim.optimization.analysis.best_compromise",
            return_value=0,
        ), patch("sys.argv", [
            "emulsim", "design",
            "--d32", "2e-6", "--d32-tol", "5e-7",
            "--n-initial", "2", "--max-iter", "1",
        ]):
            main()
        mock_run.assert_called_once()
