"""Tests for the unified UQ schema + engine (Nodes 13, 18, 30).

Verifies:
  1. UncertaintyKind enum covers the six classes from the consensus plan.
  2. UncertaintySource samples from the declared distribution.
  3. UnifiedUncertaintyResult percentile getters work for the canonical
     M1->L4 outputs (d32, pore_size_mean, G_DN).
  4. ``kinds_present`` aggregates UncertaintyKind across a multi-source spec.
  5. UnifiedUncertaintyEngine.run_m1l4 returns a schema-conformant result.
  6. CalibrationStore posteriors are absorbed and actually perturb the
     MC (Audit N2 closed in Node 30).
  7. Malformed / unknown-attribute posterior sources are skipped at
     DEBUG log level and surface via kinds_declared_but_not_sampled.
  8. The removed legacy modules (uncertainty_core, uncertainty_propagation)
     are no longer importable.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from emulsim.uncertainty_unified import (
    OutputUncertainty,
    UncertaintyKind,
    UncertaintySource,
    UnifiedUncertaintyEngine,
    UnifiedUncertaintyResult,
    UnifiedUncertaintySpec,
)


# ─── Schema ────────────────────────────────────────────────────────────────


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

    def test_raw_samples_round_trip(self):
        arr = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
        out = OutputUncertainty(
            name="x", units="-",
            mean=float(np.mean(arr)),
            p5=float(np.percentile(arr, 5)),
            p50=float(np.median(arr)),
            p95=float(np.percentile(arr, 95)),
            n_samples=arr.size,
            raw_samples=arr,
        )
        np.testing.assert_array_equal(out.raw_samples, arr)


# ─── Engine ────────────────────────────────────────────────────────────────


@pytest.fixture(scope="module")
def smoke_params():
    from emulsim.config import load_config
    cfg = Path(__file__).resolve().parents[1] / "configs" / "fast_smoke.toml"
    if not cfg.exists():
        pytest.skip("fast_smoke.toml missing")
    return load_config(cfg)


class TestUnifiedEngine:
    def test_run_m1l4_emits_unified_result(self, smoke_params):
        engine = UnifiedUncertaintyEngine()
        result = engine.run_m1l4(smoke_params, n_samples=3, seed=42)
        assert isinstance(result, UnifiedUncertaintyResult)
        assert result.n_samples == 3
        assert UncertaintyKind.MATERIAL_PROPERTY in result.kinds_sampled
        assert result.get("d32") is not None
        assert result.get("pore_size_mean") is not None
        assert result.get("G_DN") is not None
        assert result.source_label.endswith(".M1L4")
        # Node 30: raw samples carried through for downstream analyses.
        d32 = result.get("d32")
        assert d32.raw_samples is not None
        assert d32.raw_samples.shape[0] >= 1

    def test_calibration_posterior_absorbed_into_spec(self):
        """Engine adds CALIBRATION_POSTERIOR sources from store entries."""
        from emulsim.calibration.calibration_data import CalibrationEntry
        from emulsim.calibration.calibration_store import CalibrationStore

        store = CalibrationStore()
        # Entry WITHOUT posterior -> no source added.
        store.add(CalibrationEntry(
            profile_key="x", parameter_name="breakage_C1",
            measured_value=0.986, units="-", confidence="medium",
            source_reference="ref", target_module="L1",
            posterior_uncertainty=0.0,
        ))
        # Entry WITH posterior -> source added.
        store.add(CalibrationEntry(
            profile_key="x", parameter_name="coalescence_C5",
            measured_value=2.28e13, units="-", confidence="high",
            source_reference="ref", target_module="L1",
            posterior_uncertainty=5e12,
        ))

        engine = UnifiedUncertaintyEngine(calibration_store=store)
        kinds = engine.spec.kinds_present()
        assert UncertaintyKind.CALIBRATION_POSTERIOR in kinds
        post_sources = [s for s in engine.spec.sources
                        if s.kind == UncertaintyKind.CALIBRATION_POSTERIOR]
        assert len(post_sources) == 1
        assert post_sources[0].name.endswith("coalescence_C5")
        assert post_sources[0].std == pytest.approx(5e12)

    def test_posterior_now_actually_sampled(self, smoke_params):
        """Node 30 closure of Audit N2: posteriors reach kinds_sampled.

        The v7.0 engine absorbed posteriors into the spec but never
        propagated them through the MC, recording them as
        ``kinds_declared_but_not_sampled``. Node 30 wires posterior
        injection into the sampler, so a well-formed L1 posterior with
        a known KernelConfig attribute must land in ``kinds_sampled``
        and NOT in the declared-but-not-sampled set.
        """
        from emulsim.calibration.calibration_data import CalibrationEntry
        from emulsim.calibration.calibration_store import CalibrationStore

        store = CalibrationStore()
        # coalescence_C5 is a real KernelConfig attribute; L1 dispatches
        # posteriors to params.emulsification.kernels.
        store.add(CalibrationEntry(
            profile_key="x", parameter_name="coalescence_C5",
            measured_value=2.28e13, units="-", confidence="high",
            source_reference="ref", target_module="L1",
            posterior_uncertainty=5e12,
        ))
        engine = UnifiedUncertaintyEngine(calibration_store=store)
        result = engine.run_m1l4(smoke_params, n_samples=3, seed=42)

        assert UncertaintyKind.MATERIAL_PROPERTY in result.kinds_sampled
        # Audit N2 CLOSED: posterior now actually perturbs the MC.
        assert UncertaintyKind.CALIBRATION_POSTERIOR in result.kinds_sampled
        assert (UncertaintyKind.CALIBRATION_POSTERIOR
                not in result.kinds_declared_but_not_sampled)
        assert "NOT sampled" not in result.summary()

    def test_posterior_actually_perturbs_output(self, smoke_params):
        """A real posterior on a sensitive kernel coefficient must widen
        the d32 interval relative to a no-posterior run with the same
        seed + n_samples. Uses a large std to make the effect robust
        against MC noise for small n."""
        from emulsim.calibration.calibration_data import CalibrationEntry
        from emulsim.calibration.calibration_store import CalibrationStore

        baseline = UnifiedUncertaintyEngine().run_m1l4(
            smoke_params, n_samples=12, seed=1,
        )
        baseline_width = baseline.get("d32").ci_width_relative

        store = CalibrationStore()
        # Big posterior std on breakage_C1 drives a wide d32 spread.
        store.add(CalibrationEntry(
            profile_key="x", parameter_name="breakage_C1",
            measured_value=0.986, units="-", confidence="high",
            source_reference="ref", target_module="L1",
            posterior_uncertainty=0.5,
        ))
        engine = UnifiedUncertaintyEngine(calibration_store=store)
        perturbed = engine.run_m1l4(smoke_params, n_samples=12, seed=1)
        perturbed_width = perturbed.get("d32").ci_width_relative

        # Posterior-enabled run widens the interval. Use a conservative
        # ratio (≥1.1×) to keep the test robust against small-N MC noise;
        # in practice the ratio is much larger for this std.
        assert perturbed_width > baseline_width * 1.1, (
            f"posterior did not widen d32 interval: "
            f"baseline={baseline_width:.4f} perturbed={perturbed_width:.4f}"
        )

    def test_unknown_posterior_attribute_skipped(self, smoke_params, caplog):
        """A posterior targeting a nonexistent attribute must NOT land in
        kinds_sampled, must be logged at DEBUG, and must appear in
        kinds_declared_but_not_sampled (sole posterior failed)."""
        import logging
        from emulsim.calibration.calibration_data import CalibrationEntry
        from emulsim.calibration.calibration_store import CalibrationStore

        caplog.set_level(logging.DEBUG, logger="emulsim.uncertainty_unified")
        store = CalibrationStore()
        store.add(CalibrationEntry(
            profile_key="x", parameter_name="not_a_real_field",
            measured_value=1.0, units="-", confidence="low",
            source_reference="ref", target_module="L1",
            posterior_uncertainty=0.1,
        ))
        engine = UnifiedUncertaintyEngine(calibration_store=store)
        result = engine.run_m1l4(smoke_params, n_samples=2, seed=0)

        assert (UncertaintyKind.CALIBRATION_POSTERIOR
                not in result.kinds_sampled)
        assert (UncertaintyKind.CALIBRATION_POSTERIOR
                in result.kinds_declared_but_not_sampled)
        # DEBUG log should mention the attribute issue.
        skip_msgs = [r for r in caplog.records
                     if "no attribute" in r.message
                     and "not_a_real_field" in r.message]
        assert len(skip_msgs) >= 1

    def test_legacy_modules_are_gone(self):
        """Node 30 removed uncertainty_core and uncertainty_propagation."""
        with pytest.raises(ImportError):
            import emulsim.uncertainty_core  # noqa: F401
        with pytest.raises(ImportError):
            import emulsim.uncertainty_propagation  # noqa: F401
