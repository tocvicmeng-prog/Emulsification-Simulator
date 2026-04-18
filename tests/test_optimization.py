"""Tests for the optimisation engine."""

import numpy as np
import pytest
from pathlib import Path

from emulsim.datatypes import SimulationParameters
from emulsim.optimization.objectives import (
    PARAM_BOUNDS,
    PARAM_NAMES,
    LOG_SCALE_INDICES,
)


class TestObjectives:
    def test_param_names_count(self):
        assert len(PARAM_NAMES) == 7

    def test_bounds_shape(self):
        assert PARAM_BOUNDS.shape == (7, 2)

    def test_bounds_ordered(self):
        for i in range(7):
            assert PARAM_BOUNDS[i, 0] < PARAM_BOUNDS[i, 1]

    def test_log_scale_indices_valid(self):
        for idx in LOG_SCALE_INDICES:
            assert 0 <= idx < 7


class TestSearchSpace:
    def test_roundtrip(self):
        """Transform to search space and back should be identity."""
        from emulsim.optimization.engine import _to_search_space, _from_search_space
        x = np.array([10000, 20, 0.7, 363, 0.167, 2.0, 86400.0])
        x_ss = _to_search_space(x)
        x_back = _from_search_space(x_ss)
        np.testing.assert_allclose(x_back, x, rtol=1e-10)

    def test_log_transforms(self):
        from emulsim.optimization.engine import _to_search_space
        x = np.array([10000, 20, 0.7, 363, 0.167, 2.0, 86400.0])
        x_ss = _to_search_space(x)
        # RPM (index 0) should be log10(10000) = 4
        assert x_ss[0] == pytest.approx(4.0)
        # agarose_frac (index 2) should NOT be log-transformed
        assert x_ss[2] == pytest.approx(0.7)


class TestOptimizationEngine:
    def test_engine_initializes(self):
        from emulsim.optimization.engine import OptimizationEngine
        engine = OptimizationEngine(n_initial=3, max_iterations=2)
        assert engine.n_initial == 3
        assert engine.max_iterations == 2

    def test_evaluate_uses_trust_aware_objectives(self, monkeypatch):
        from emulsim.optimization import engine as opt_engine
        from emulsim.optimization.engine import OptimizationEngine

        sentinel_result = object()
        engine = OptimizationEngine(
            n_initial=1,
            max_iterations=0,
            output_dir=Path("output/test_optimization_trust"),
            template_params=SimulationParameters(),
        )

        class FakeOrchestrator:
            def run_single(self, params):
                return sentinel_result

        called = {}

        def fake_compute_objectives(result):
            called["result"] = result
            return np.array([1.0, 2.0, 3.0])

        monkeypatch.setattr(engine, "orchestrator", FakeOrchestrator())
        monkeypatch.setattr(opt_engine, "compute_objectives_trust_aware", fake_compute_objectives)
        monkeypatch.setattr(opt_engine, "check_constraints", lambda result: (True, []))

        objectives, result = engine._evaluate(np.array([10000, 20, 0.7, 363, 0.167, 2.0, 86400.0]))

        assert result is sentinel_result
        assert called["result"] is sentinel_result
        np.testing.assert_allclose(objectives, np.array([1.0, 2.0, 3.0]))

    @pytest.mark.slow
    def test_short_campaign(self):
        """Run a minimal 3+2 campaign to verify the engine works end-to-end."""
        from emulsim.optimization.engine import OptimizationEngine
        engine = OptimizationEngine(n_initial=3, max_iterations=2)
        state = engine.run()
        assert len(state.X_observed) == 5  # 3 init + 2 BO
        assert len(state.Y_observed) == 5
        assert state.pareto_X.shape[0] >= 1
        assert state.iteration >= 1


# ─── Node 6 (v6.1): Per-tier optimizer penalty ────────────────────────────


class TestNode6TrustPenalty:
    """Verifies tier-based optimizer penalty and Pareto exclusion semantics.

    Acceptance for Node 6:
      - trust_penalty_for_tier ordering matches the consensus-plan tier rank.
      - QUALITATIVE_TREND penalty puts max objective above the engine
        REF_POINT (5.0) so candidates drop off the Pareto front via
        feasible_mask in engine.run().
      - UNSUPPORTED penalty is even larger.
      - OptimizationState carries pareto_evidence_tiers aligned with pareto_X.
    """

    def test_penalty_ordering(self):
        """Penalties must be strictly weaker -> stronger from VALIDATED to UNSUPPORTED."""
        from emulsim.datatypes import ModelEvidenceTier
        from emulsim.optimization.objectives import trust_penalty_for_tier

        order = [
            ModelEvidenceTier.VALIDATED_QUANTITATIVE,
            ModelEvidenceTier.CALIBRATED_LOCAL,
            ModelEvidenceTier.SEMI_QUANTITATIVE,
            ModelEvidenceTier.QUALITATIVE_TREND,
            ModelEvidenceTier.UNSUPPORTED,
        ]
        penalties = [trust_penalty_for_tier(t) for t in order]
        # Strictly non-decreasing; ideally strictly increasing
        for prev, curr in zip(penalties, penalties[1:]):
            assert curr >= prev, f"Tier penalties out of order: {penalties}"
        # VALIDATED gives a (small) bonus
        assert penalties[0] < 0
        # SEMI is mild (<<1), QUALITATIVE crosses REF_POINT, UNSUPPORTED much larger
        assert penalties[2] < 1.0
        assert penalties[3] >= 5.0
        assert penalties[4] >= 10.0

    def test_qualitative_tier_excluded_from_pareto(self):
        """A QUALITATIVE_TREND candidate's max objective exceeds REF_POINT (5.0)."""
        from emulsim.datatypes import (
            ModelEvidenceTier,
        )
        from emulsim.optimization.objectives import compute_objectives_trust_aware

        # Build a minimal FullResult with a QUALITATIVE_TREND manifest and
        # near-target objectives so the only thing pushing max above REF_POINT
        # is the trust penalty itself.
        result = _make_full_result_with_tier(ModelEvidenceTier.QUALITATIVE_TREND)
        objectives = compute_objectives_trust_aware(result)
        assert objectives.max() > 5.0, (
            f"QUALITATIVE_TREND objectives {objectives} should exceed REF_POINT (5.0). "
            "If this fails, the tier penalty is too small to enforce Pareto exclusion."
        )

        # SEMI_QUANTITATIVE with the same near-target objectives stays below
        # REF_POINT and is therefore retained on the Pareto front.
        result_semi = _make_full_result_with_tier(ModelEvidenceTier.SEMI_QUANTITATIVE)
        objectives_semi = compute_objectives_trust_aware(result_semi)
        assert objectives_semi.max() <= 5.0, (
            f"SEMI_QUANTITATIVE objectives {objectives_semi} should be <= REF_POINT."
        )

    def test_optimization_state_carries_pareto_tiers(self):
        """OptimizationState.pareto_evidence_tiers is populated for surviving Pareto X."""
        from emulsim.datatypes import OptimizationState
        # Just verify the dataclass field exists and defaults to [] (engine
        # population is exercised by the slow test_short_campaign path; here
        # we only check the contract.)
        state = OptimizationState(
            X_observed=np.zeros((1, 7)),
            Y_observed=np.zeros((1, 3)),
            pareto_X=np.zeros((0, 7)),
            pareto_Y=np.zeros((0, 3)),
            iteration=0,
            hypervolume=0.0,
        )
        assert state.pareto_evidence_tiers == []


def _make_full_result_with_tier(tier):
    """Build a minimal FullResult whose run_report.compute_min_tier returns `tier`.

    Used by the Node 6 tests to drive compute_objectives_trust_aware without
    running the full pipeline. Base objectives end up near zero (d32 close to
    target etc.) so the assertion can isolate the tier penalty contribution.
    """
    from emulsim.datatypes import (
        EmulsificationResult, GelationResult, CrosslinkingResult,
        MechanicalResult, FullResult, ModelManifest, RunReport,
        SimulationParameters,
    )
    params = SimulationParameters()
    emul = EmulsificationResult(
        d_bins=np.array([1e-6, 2e-6, 4e-6]),
        n_d=np.array([1.0, 1.0, 1.0]),
        d32=2.0e-6, d43=2.0e-6,
        d10=1.5e-6, d50=2.0e-6, d90=2.5e-6,
        span=0.5,
        total_volume_fraction=0.05,
        converged=True,
    )
    gel = GelationResult(
        r_grid=np.array([0.0, 1e-6]),
        phi_field=np.array([0.5, 0.5]),
        pore_size_mean=80e-9,
        pore_size_std=10e-9,
        pore_size_distribution=np.array([80e-9]),
        porosity=0.7,
        alpha_final=1.0,
        char_wavelength=80e-9,
    )
    xl = CrosslinkingResult(
        t_array=np.array([0.0, 1.0]),
        X_array=np.array([0.0, 1.0]),
        nu_e_array=np.array([0.0, 1e23]),
        Mc_array=np.array([1e6, 1e3]),
        xi_array=np.array([1e-6, 1e-8]),
        G_chitosan_array=np.array([0.0, 5e3]),
        p_final=0.5, nu_e_final=1e23, Mc_final=1e3,
        xi_final=1e-8, G_chitosan_final=5e3,
    )
    mech = MechanicalResult(
        G_agarose=5e3, G_chitosan=5e3,
        G_DN=10e3,                    # log10(10e3) = 4.0 -> f3 = 0
        E_star=30e3,
        delta_array=np.array([0.0]),
        F_array=np.array([0.0]),
        rh_array=np.array([1e-9]),
        Kav_array=np.array([1.0]),
        pore_size_mean=80e-9,
        xi_mesh=1e-8,
    )
    full = FullResult(
        parameters=params, emulsification=emul,
        gelation=gel, crosslinking=xl, mechanical=mech,
    )
    full.run_report = RunReport(
        model_graph=[ModelManifest(model_name="synthetic", evidence_tier=tier)],
        trust_level="CAUTION",
    )
    return full
