"""Tests for v6.1 ModelEvidenceTier, ModelManifest, and RunReport."""

from __future__ import annotations

import numpy as np
import pytest

from emulsim.datatypes import (
    ModelEvidenceTier,
    ModelManifest,
    RunReport,
    EmulsificationResult,
    GelationResult,
    CrosslinkingResult,
    MechanicalResult,
    FullResult,
    SimulationParameters,
)


# ── ModelEvidenceTier ────────────────────────────────────────────────────


class TestModelEvidenceTier:

    def test_tier_count(self):
        assert len(ModelEvidenceTier) == 5

    def test_tier_ordering(self):
        tiers = list(ModelEvidenceTier)
        assert tiers[0] == ModelEvidenceTier.VALIDATED_QUANTITATIVE
        assert tiers[-1] == ModelEvidenceTier.UNSUPPORTED

    def test_tier_values_are_strings(self):
        for tier in ModelEvidenceTier:
            assert isinstance(tier.value, str)

    def test_validated_is_strongest(self):
        tiers = list(ModelEvidenceTier)
        assert tiers.index(ModelEvidenceTier.VALIDATED_QUANTITATIVE) == 0

    def test_unsupported_is_weakest(self):
        tiers = list(ModelEvidenceTier)
        assert tiers.index(ModelEvidenceTier.UNSUPPORTED) == 4


# ── ModelManifest ────────────────────────────────────────────────────────


class TestModelManifest:

    def test_construction_minimal(self):
        m = ModelManifest(model_name="L1.PBE.FixedPivot")
        assert m.model_name == "L1.PBE.FixedPivot"
        assert m.evidence_tier == ModelEvidenceTier.SEMI_QUANTITATIVE

    def test_construction_full(self):
        m = ModelManifest(
            model_name="L2.Pore.CahnHilliard2D",
            evidence_tier=ModelEvidenceTier.CALIBRATED_LOCAL,
            valid_domain={"Re": (100, 1e6)},
            calibration_ref="CAL-001",
            assumptions=["2D approximation"],
            diagnostics={"converged": True},
        )
        assert m.evidence_tier == ModelEvidenceTier.CALIBRATED_LOCAL
        assert "Re" in m.valid_domain
        assert m.calibration_ref == "CAL-001"

    def test_default_lists_are_independent(self):
        m1 = ModelManifest(model_name="A")
        m2 = ModelManifest(model_name="B")
        m1.assumptions.append("test")
        assert "test" not in m2.assumptions


# ── RunReport ────────────────────────────────────────────────────────────


class TestRunReport:

    def test_empty_report(self):
        rr = RunReport()
        assert rr.compute_min_tier() == ModelEvidenceTier.UNSUPPORTED

    def test_single_model(self):
        m = ModelManifest(model_name="L1", evidence_tier=ModelEvidenceTier.SEMI_QUANTITATIVE)
        rr = RunReport(model_graph=[m])
        assert rr.compute_min_tier() == ModelEvidenceTier.SEMI_QUANTITATIVE

    def test_min_tier_picks_weakest(self):
        m1 = ModelManifest(model_name="L1", evidence_tier=ModelEvidenceTier.CALIBRATED_LOCAL)
        m2 = ModelManifest(model_name="L2", evidence_tier=ModelEvidenceTier.QUALITATIVE_TREND)
        rr = RunReport(model_graph=[m1, m2])
        assert rr.compute_min_tier() == ModelEvidenceTier.QUALITATIVE_TREND

    def test_all_validated(self):
        manifests = [
            ModelManifest(model_name=f"L{i}", evidence_tier=ModelEvidenceTier.VALIDATED_QUANTITATIVE)
            for i in range(4)
        ]
        rr = RunReport(model_graph=manifests)
        assert rr.compute_min_tier() == ModelEvidenceTier.VALIDATED_QUANTITATIVE

    def test_trust_fields(self):
        rr = RunReport(
            trust_level="CAUTION",
            trust_warnings=["d32 large"],
            trust_blockers=[],
        )
        assert rr.trust_level == "CAUTION"
        assert len(rr.trust_warnings) == 1


# ── Backward Compatibility ───────────────────────────────────────────────


class TestBackwardCompatibility:

    def test_emulsification_result_no_manifest(self):
        e = EmulsificationResult(
            d_bins=np.array([1e-5]),
            n_d=np.array([1e10]),
            d32=1e-5, d43=1e-5, d10=5e-6, d50=1e-5, d90=2e-5,
            span=1.5, total_volume_fraction=0.05, converged=True,
        )
        assert e.model_manifest is None

    def test_gelation_result_no_manifest(self):
        g = GelationResult(
            r_grid=np.array([0.0]),
            phi_field=np.array([0.5]),
            pore_size_mean=100e-9,
            pore_size_std=20e-9,
            pore_size_distribution=np.array([100e-9]),
            porosity=0.5,
            alpha_final=0.9,
            char_wavelength=200e-9,
        )
        assert g.model_manifest is None

    def test_crosslinking_result_no_manifest(self):
        x = CrosslinkingResult(
            t_array=np.array([0.0, 1.0]),
            X_array=np.array([0.0, 0.1]),
            nu_e_array=np.array([0.0, 1e20]),
            Mc_array=np.array([1e10, 1e4]),
            xi_array=np.array([1e-6, 50e-9]),
            G_chitosan_array=np.array([0.0, 1e4]),
            p_final=0.1,
            nu_e_final=1e20,
            Mc_final=1e4,
            xi_final=50e-9,
            G_chitosan_final=1e4,
        )
        assert x.model_manifest is None

    def test_mechanical_result_no_manifest(self):
        m = MechanicalResult(
            G_agarose=5e4, G_chitosan=1e4, G_DN=6e4, E_star=1.5e5,
            delta_array=np.array([1e-6]),
            F_array=np.array([1e-3]),
            rh_array=np.array([5e-9]),
            Kav_array=np.array([0.5]),
            pore_size_mean=100e-9,
            xi_mesh=50e-9,
        )
        assert m.model_manifest is None
        assert m.model_used == "phenomenological"

    def test_full_result_no_run_report(self):
        """FullResult can be created without run_report (backward compat)."""
        e = EmulsificationResult(
            d_bins=np.array([1e-5]), n_d=np.array([1e10]),
            d32=1e-5, d43=1e-5, d10=5e-6, d50=1e-5, d90=2e-5,
            span=1.5, total_volume_fraction=0.05, converged=True,
        )
        g = GelationResult(
            r_grid=np.array([0.0]), phi_field=np.array([0.5]),
            pore_size_mean=100e-9, pore_size_std=20e-9,
            pore_size_distribution=np.array([100e-9]),
            porosity=0.5, alpha_final=0.9, char_wavelength=200e-9,
        )
        x = CrosslinkingResult(
            t_array=np.array([0.0]), X_array=np.array([0.0]),
            nu_e_array=np.array([0.0]), Mc_array=np.array([1e10]),
            xi_array=np.array([1e-6]), G_chitosan_array=np.array([0.0]),
            p_final=0.0, nu_e_final=0.0, Mc_final=1e10, xi_final=1e-6,
            G_chitosan_final=0.0,
        )
        m = MechanicalResult(
            G_agarose=5e4, G_chitosan=0.0, G_DN=5e4, E_star=1.2e5,
            delta_array=np.array([1e-6]), F_array=np.array([1e-3]),
            rh_array=np.array([5e-9]), Kav_array=np.array([0.5]),
            pore_size_mean=100e-9, xi_mesh=1e-6,
        )
        fr = FullResult(
            parameters=SimulationParameters(),
            emulsification=e, gelation=g,
            crosslinking=x, mechanical=m,
        )
        assert fr.run_report is None
