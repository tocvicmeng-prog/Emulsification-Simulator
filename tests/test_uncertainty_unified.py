"""Tests for Node 13 (v6.1, F7): unified UQ schema scaffold.

Acceptance for Node 13:
  1. UncertaintyKind enum covers the six classes from the consensus plan.
  2. UncertaintySource samples from the declared distribution.
  3. UnifiedUncertaintyResult percentile getters work for the canonical
     M1->L4 outputs (d32, pore, G_DN).
  4. The two adapters produce a UnifiedUncertaintyResult that preserves
     the legacy median/CI within numerical tolerance.
  5. The kinds_present helper correctly aggregates UncertaintyKind across
     a multi-source spec.
"""

from __future__ import annotations

import numpy as np
import pytest

from emulsim.uncertainty_unified import (
    OutputUncertainty,
    UncertaintyKind,
    UncertaintySource,
    UnifiedUncertaintyResult,
    UnifiedUncertaintySpec,
    from_m1_contract_uq,
    from_m1l4_result,
)


class TestUncertaintyKind:
    def test_consensus_plan_kinds_present(self):
        names = {k.value for k in UncertaintyKind}
        for required in (
            "measurement",
            "material_property",
            "calibration_posterior",
            "model_form",
            "numerical",
            "scale_up",
        ):
            assert required in names


class TestSourceSampling:
    def test_normal_sample_shape_and_mean(self):
        src = UncertaintySource(
            name="props.sigma", kind=UncertaintyKind.MATERIAL_PROPERTY,
            distribution="normal", value=5e-3, std=5e-4, units="N/m",
        )
        rng = np.random.default_rng(0)
        s = src.sample(rng, 5000)
        assert s.shape == (5000,)
        assert abs(np.mean(s) - 5e-3) < 5e-5

    def test_lognormal_positive(self):
        src = UncertaintySource(
            name="props.k_xlink_0", kind=UncertaintyKind.CALIBRATION_POSTERIOR,
            distribution="lognormal", value=1.33e4, std=0.30,
        )
        s = src.sample(np.random.default_rng(0), 1000)
        assert (s > 0).all()

    def test_unsupported_distribution_raises(self):
        src = UncertaintySource(
            name="x", kind=UncertaintyKind.NUMERICAL,
            distribution="bogus", value=1.0, std=0.1,
        )
        with pytest.raises(ValueError, match="bogus"):
            src.sample(np.random.default_rng(0), 10)


class TestSpecAggregation:
    def test_kinds_present_across_sources(self):
        spec = UnifiedUncertaintySpec()
        spec.add(UncertaintySource(
            "a", UncertaintyKind.MEASUREMENT, "normal", 1.0, 0.1,
        ))
        spec.add(UncertaintySource(
            "b", UncertaintyKind.MODEL_FORM, "normal", 1.0, 0.1,
        ))
        kinds = spec.kinds_present()
        assert UncertaintyKind.MEASUREMENT in kinds
        assert UncertaintyKind.MODEL_FORM in kinds
        assert UncertaintyKind.SCALE_UP not in kinds


class TestUnifiedResult:
    def test_get_returns_named_output(self):
        out = OutputUncertainty(
            name="d32", units="m", mean=2e-6, p5=1.5e-6, p50=2e-6, p95=2.5e-6,
            n_samples=100,
        )
        result = UnifiedUncertaintyResult(outputs=[out], n_samples=100)
        assert result.get("d32") is out
        assert result.get("nonexistent") is None

    def test_relative_ci_width(self):
        out = OutputUncertainty(
            name="d32", units="m", mean=2e-6, p5=1e-6, p50=2e-6, p95=3e-6,
            n_samples=100,
        )
        # (3e-6 - 1e-6) / 2e-6 = 1.0
        assert out.ci_width_relative == pytest.approx(1.0)


class TestAdapters:
    def test_from_m1l4_preserves_median(self):
        """Adapter median must match the median of the underlying sample array."""
        # Build a fake legacy UncertaintyResult-shaped object.
        class _Legacy:
            n_samples = 200
            n_failed = 0
            d32_samples = np.random.default_rng(1).normal(2e-6, 1e-7, 200)
            pore_samples = np.random.default_rng(2).normal(80e-9, 5e-9, 200)
            G_DN_samples = np.random.default_rng(3).normal(8000.0, 500.0, 200)
            d32_median = float(np.median(d32_samples))
            pore_median = float(np.median(pore_samples))
            G_DN_median = float(np.median(G_DN_samples))
            d32_ci = (float(np.percentile(d32_samples, 5)),
                      float(np.percentile(d32_samples, 95)))
            pore_ci = (float(np.percentile(pore_samples, 5)),
                       float(np.percentile(pore_samples, 95)))
            G_DN_ci = (float(np.percentile(G_DN_samples, 5)),
                       float(np.percentile(G_DN_samples, 95)))

        unified = from_m1l4_result(_Legacy(), source_label="test_legacy")

        assert unified.source_label == "test_legacy"
        assert unified.kinds_sampled == {UncertaintyKind.MATERIAL_PROPERTY}
        assert unified.n_samples == 200

        d32 = unified.get("d32")
        assert d32 is not None
        assert d32.p50 == pytest.approx(_Legacy.d32_median, rel=1e-9)
        assert d32.p5 == pytest.approx(_Legacy.d32_ci[0], rel=1e-9)
        assert d32.p95 == pytest.approx(_Legacy.d32_ci[1], rel=1e-9)

    def test_from_m1l4_handles_none(self):
        """None input -> empty unified result, no crash."""
        unified = from_m1l4_result(None)
        assert unified.outputs == []
        assert unified.n_samples == 0

    def test_from_m1_contract_uq_percentiles(self):
        samples = np.linspace(40.0, 60.0, 101)  # uniform 40..60
        unified = from_m1_contract_uq(samples, units="mol/m^3")
        out = unified.get("estimated_q_max")
        assert out is not None
        assert out.p50 == pytest.approx(50.0, abs=1e-9)
        assert out.p5 == pytest.approx(41.0, abs=1e-9)
        assert out.p95 == pytest.approx(59.0, abs=1e-9)
        assert out.n_samples == 101
        assert UncertaintyKind.MATERIAL_PROPERTY in unified.kinds_sampled

    def test_summary_runs(self):
        """summary() must produce a non-empty multi-line string."""
        unified = from_m1_contract_uq(np.array([1.0, 2.0, 3.0]))
        s = unified.summary()
        assert "estimated_q_max" in s
        assert s.count("\n") >= 1


# ─── Node 18 (v7.0, P3): Unified MC engine ────────────────────────────────


class TestUnifiedEngine:
    """Verifies the v7.0 single-entrypoint engine.

    The structural merge of legacy uncertainty_core + uncertainty_propagation
    is deferred to v7.1; v7.0 ships a single engine that delegates and
    unifies the output schema.
    """

    @pytest.fixture
    def smoke_params(self):
        from emulsim.config import load_config
        from pathlib import Path
        cfg = Path(__file__).resolve().parents[1] / "configs" / "fast_smoke.toml"
        if not cfg.exists():
            pytest.skip("fast_smoke.toml missing")
        return load_config(cfg)

    def test_run_m1l4_emits_unified_result(self, smoke_params):
        from emulsim.uncertainty_unified import (
            UnifiedUncertaintyEngine,
            UnifiedUncertaintyResult,
            UncertaintyKind,
        )
        engine = UnifiedUncertaintyEngine()
        result = engine.run_m1l4(smoke_params, n_samples=3, seed=42)
        assert isinstance(result, UnifiedUncertaintyResult)
        assert result.n_samples == 3
        assert UncertaintyKind.MATERIAL_PROPERTY in result.kinds_sampled
        assert result.get("d32") is not None
        assert result.get("pore_size_mean") is not None
        assert result.get("G_DN") is not None
        assert result.source_label.endswith(".M1L4")

    def test_calibration_posterior_absorbed(self):
        """Engine adds CALIBRATION_POSTERIOR sources from store entries."""
        from emulsim.calibration.calibration_data import CalibrationEntry
        from emulsim.calibration.calibration_store import CalibrationStore
        from emulsim.uncertainty_unified import (
            UnifiedUncertaintyEngine,
            UncertaintyKind,
        )

        store = CalibrationStore()
        # Entry WITHOUT posterior -> no source added
        store.add(CalibrationEntry(
            profile_key="x", parameter_name="breakage_C1",
            measured_value=0.986, units="-", confidence="medium",
            source_reference="ref", target_module="L1",
            posterior_uncertainty=0.0,
        ))
        # Entry WITH posterior -> source added
        store.add(CalibrationEntry(
            profile_key="x", parameter_name="coalescence_C5",
            measured_value=2.28e13, units="-", confidence="high",
            source_reference="ref", target_module="L1",
            posterior_uncertainty=5e12,
        ))

        engine = UnifiedUncertaintyEngine(calibration_store=store)
        kinds = engine.spec.kinds_present()
        assert UncertaintyKind.CALIBRATION_POSTERIOR in kinds
        # Only 1 source absorbed (the one with posterior > 0)
        post_sources = [s for s in engine.spec.sources
                        if s.kind == UncertaintyKind.CALIBRATION_POSTERIOR]
        assert len(post_sources) == 1
        assert post_sources[0].name.endswith("coalescence_C5")
        assert post_sources[0].std == pytest.approx(5e12)

    def test_n2_no_posterior_overclaim(self, smoke_params):
        """Audit N2 (v7.0.1): posterior NOT in kinds_sampled when not actually sampled.

        The v7.0 engine absorbs CalibrationStore posteriors into the spec
        but does not yet propagate them through the legacy MC engine. The
        result's ``kinds_sampled`` must therefore NOT advertise
        CALIBRATION_POSTERIOR; it must record the declared-but-not-sampled
        condition in ``kinds_declared_but_not_sampled`` instead so
        downstream UI/optimizer code does not over-attribute the interval.
        """
        from emulsim.calibration.calibration_data import CalibrationEntry
        from emulsim.calibration.calibration_store import CalibrationStore
        from emulsim.uncertainty_unified import (
            UnifiedUncertaintyEngine, UncertaintyKind,
        )

        store = CalibrationStore()
        store.add(CalibrationEntry(
            profile_key="x", parameter_name="coalescence_C5",
            measured_value=2.28e13, units="-", confidence="high",
            source_reference="ref", target_module="L1",
            posterior_uncertainty=5e12,
        ))
        engine = UnifiedUncertaintyEngine(calibration_store=store)
        result = engine.run_m1l4(smoke_params, n_samples=3, seed=42)

        # The legacy engine sampled MATERIAL_PROPERTY perturbations only.
        assert UncertaintyKind.MATERIAL_PROPERTY in result.kinds_sampled
        # The posterior was absorbed into the spec but NOT actually sampled.
        # That fact must be honestly recorded — kinds_sampled must NOT lie.
        assert UncertaintyKind.CALIBRATION_POSTERIOR not in result.kinds_sampled
        assert (UncertaintyKind.CALIBRATION_POSTERIOR
                in result.kinds_declared_but_not_sampled)
        # The summary string also surfaces the limitation
        assert "NOT sampled" in result.summary()
