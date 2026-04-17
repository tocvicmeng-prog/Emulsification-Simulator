"""Tests for Node 30b: streamlit UQ panel spec-building helper.

The panel's streamlit rendering function is not unit-tested directly
(requires a streamlit runtime); its pure helper ``build_uncertainty_spec``
is. An import smoke-test verifies that ``render_uncertainty_panel`` is
a callable exported by the panel module for the app's import path.
"""

from __future__ import annotations

import pytest

from emulsim.uncertainty_unified import (
    UncertaintyKind,
    UnifiedUncertaintySpec,
)
from emulsim.visualization.panels.uncertainty import (
    CustomSourceInput,
    build_uncertainty_spec,
    count_store_posteriors,
)


class TestBuildSpec:
    def test_default_has_no_extra_sources(self):
        spec = build_uncertainty_spec(n_samples=100, seed=42)
        assert isinstance(spec, UnifiedUncertaintySpec)
        assert spec.n_samples == 100
        assert spec.seed == 42
        assert spec.sources == []

    def test_custom_source_appended(self):
        spec = build_uncertainty_spec(
            n_samples=50, seed=7,
            custom_sources=[CustomSourceInput(
                name="props.sigma",
                kind=UncertaintyKind.MATERIAL_PROPERTY,
                value=5e-3,
                std=5e-4,
            )],
        )
        assert len(spec.sources) == 1
        assert spec.sources[0].name == "props.sigma"
        assert spec.sources[0].kind == UncertaintyKind.MATERIAL_PROPERTY
        assert spec.sources[0].value == pytest.approx(5e-3)
        assert spec.sources[0].std == pytest.approx(5e-4)

    def test_blank_name_is_skipped(self):
        """User leaves name field empty → source skipped, not error."""
        spec = build_uncertainty_spec(
            n_samples=50, seed=0,
            custom_sources=[CustomSourceInput(
                name="   ",
                kind=UncertaintyKind.MEASUREMENT,
                value=1.0, std=0.1,
            )],
        )
        assert spec.sources == []

    def test_zero_std_is_skipped(self):
        """User enters std=0 → source skipped (no perturbation to add)."""
        spec = build_uncertainty_spec(
            n_samples=50, seed=0,
            custom_sources=[CustomSourceInput(
                name="placeholder",
                kind=UncertaintyKind.MEASUREMENT,
                value=1.0, std=0.0,
            )],
        )
        assert spec.sources == []

    def test_multiple_sources_preserve_order(self):
        inputs = [
            CustomSourceInput("a", UncertaintyKind.MEASUREMENT, 1.0, 0.1),
            CustomSourceInput("b", UncertaintyKind.MODEL_FORM, 2.0, 0.2),
            CustomSourceInput("c", UncertaintyKind.NUMERICAL, 3.0, 0.3),
        ]
        spec = build_uncertainty_spec(
            n_samples=100, seed=1, custom_sources=inputs,
        )
        assert [s.name for s in spec.sources] == ["a", "b", "c"]
        assert [s.kind for s in spec.sources] == [
            UncertaintyKind.MEASUREMENT,
            UncertaintyKind.MODEL_FORM,
            UncertaintyKind.NUMERICAL,
        ]

    def test_lognormal_distribution_preserved(self):
        spec = build_uncertainty_spec(
            n_samples=100, seed=1,
            custom_sources=[CustomSourceInput(
                name="x",
                kind=UncertaintyKind.MATERIAL_PROPERTY,
                value=1e-3, std=0.3,
                distribution="lognormal",
            )],
        )
        assert spec.sources[0].distribution == "lognormal"

    def test_n_samples_below_one_raises(self):
        with pytest.raises(ValueError, match=">= 1"):
            build_uncertainty_spec(n_samples=0, seed=0)


class TestCountStorePosteriors:
    def test_none_store_returns_zero(self):
        assert count_store_posteriors(None) == 0

    def test_empty_store_returns_zero(self):
        from emulsim.calibration.calibration_store import CalibrationStore
        assert count_store_posteriors(CalibrationStore()) == 0

    def test_zero_posterior_entries_not_counted(self):
        from emulsim.calibration.calibration_data import CalibrationEntry
        from emulsim.calibration.calibration_store import CalibrationStore

        store = CalibrationStore()
        store.add(CalibrationEntry(
            profile_key="x", parameter_name="breakage_C1",
            measured_value=0.986, units="-", confidence="medium",
            source_reference="ref", target_module="L1",
            posterior_uncertainty=0.0,
        ))
        assert count_store_posteriors(store) == 0

    def test_nonzero_posteriors_counted(self):
        from emulsim.calibration.calibration_data import CalibrationEntry
        from emulsim.calibration.calibration_store import CalibrationStore

        store = CalibrationStore()
        store.add(CalibrationEntry(
            profile_key="x", parameter_name="breakage_C1",
            measured_value=0.986, units="-", confidence="medium",
            source_reference="ref", target_module="L1",
            posterior_uncertainty=0.0,
        ))
        store.add(CalibrationEntry(
            profile_key="x", parameter_name="coalescence_C5",
            measured_value=2.28e13, units="-", confidence="high",
            source_reference="ref", target_module="L1",
            posterior_uncertainty=5e12,
        ))
        store.add(CalibrationEntry(
            profile_key="x", parameter_name="breakage_C2",
            measured_value=0.0115, units="-", confidence="low",
            source_reference="ref", target_module="L1",
            posterior_uncertainty=0.002,
        ))
        assert count_store_posteriors(store) == 2


class TestPanelExport:
    def test_render_uncertainty_panel_is_callable(self):
        """Import smoke test: the app's import path
        ``from emulsim.visualization.panels import render_uncertainty_panel``
        must still resolve to a callable after Node 30b."""
        from emulsim.visualization.panels import render_uncertainty_panel
        assert callable(render_uncertainty_panel)
