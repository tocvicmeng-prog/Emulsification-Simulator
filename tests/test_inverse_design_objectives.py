"""Tests for Node F3-a (v8.0 Phase 1): inverse-design TargetSpec objectives.

Validates TargetSpec construction + validation, the
compute_inverse_design_objectives function's distance math, trust-penalty
integration, and the Kav fallback when the result lacks an M3 block.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
import pytest

from emulsim.datatypes import ModelEvidenceTier
from emulsim.optimization.objectives import (
    TargetSpec,
    compute_inverse_design_objectives,
    trust_penalty_for_tier,
)


# ─── Minimal result stub that walks the duck-typed attribute chain ────────


@dataclass
class _Emul:
    d32: float = 2.0e-6
    d_mode: float = 0.0
    span: float = 1.0


@dataclass
class _Gel:
    pore_size_mean: float = 80e-9


@dataclass
class _Mech:
    G_DN: float = 10_000.0


@dataclass
class _XL:
    xi_final: float = 1e-8


@dataclass
class _Form:
    c_agarose: float = 42.0
    c_chitosan: float = 18.0
    c_genipin: float = 2.0


@dataclass
class _EmulParams:
    mode: str = "rotor_stator_legacy"
    rpm: float = 10_000.0


@dataclass
class _Params:
    formulation: _Form
    emulsification: _EmulParams


@dataclass
class _RunReport:
    """Tier stub; compute_trust_penalty reads .compute_min_tier()."""
    tier: ModelEvidenceTier = ModelEvidenceTier.SEMI_QUANTITATIVE

    def compute_min_tier(self) -> ModelEvidenceTier:
        return self.tier


@dataclass
class _Result:
    emulsification: _Emul
    gelation: _Gel
    mechanical: _Mech
    crosslinking: _XL
    parameters: _Params
    run_report: _RunReport | None = None


def _result(
    *, d32=2.0e-6, pore=80e-9, G_DN=10_000.0, d_mode=0.0,
    mode="rotor_stator_legacy",
    tier: ModelEvidenceTier | None = ModelEvidenceTier.SEMI_QUANTITATIVE,
) -> _Result:
    return _Result(
        emulsification=_Emul(d32=d32, d_mode=d_mode),
        gelation=_Gel(pore_size_mean=pore),
        mechanical=_Mech(G_DN=G_DN),
        crosslinking=_XL(),
        parameters=_Params(formulation=_Form(), emulsification=_EmulParams(mode=mode)),
        run_report=_RunReport(tier=tier) if tier else None,
    )


# ─── TargetSpec ───────────────────────────────────────────────────────────


class TestTargetSpec:
    def test_empty_spec_fails_validate(self):
        with pytest.raises(ValueError, match="no active dimension"):
            TargetSpec().validate()

    def test_d32_only_is_valid(self):
        spec = TargetSpec(d32_target=2e-6, d32_tol=0.5e-6)
        spec.validate()
        assert spec.active_dims() == ["d32"]

    def test_all_dims_active(self):
        spec = TargetSpec(
            d32_target=2e-6, d32_tol=0.5e-6,
            pore_target=80e-9, pore_tol=10e-9,
            G_DN_target=10_000, G_DN_log10_tol=0.3,
            Kav_target=0.5, Kav_tol=0.05,
        )
        assert spec.active_dims() == ["d32", "pore", "G_DN", "Kav"]

    def test_non_positive_tol_raises(self):
        with pytest.raises(ValueError, match="d32_tol"):
            TargetSpec(d32_target=2e-6, d32_tol=0.0).validate()
        with pytest.raises(ValueError, match="pore_tol"):
            TargetSpec(pore_target=80e-9, pore_tol=-1e-9).validate()


# ─── compute_inverse_design_objectives ────────────────────────────────────


class TestInverseDesignObjective:
    def test_exact_target_is_zero_without_trust_penalty(self):
        spec = TargetSpec(
            d32_target=2e-6, d32_tol=0.5e-6,
            pore_target=80e-9, pore_tol=10e-9,
            G_DN_target=10_000, G_DN_log10_tol=0.3,
        )
        r = _result(d32=2e-6, pore=80e-9, G_DN=10_000.0)
        obj = compute_inverse_design_objectives(r, spec, trust_aware=False)
        np.testing.assert_allclose(obj, [0.0, 0.0, 0.0], atol=1e-10)

    def test_one_tolerance_away_is_one(self):
        spec = TargetSpec(d32_target=2e-6, d32_tol=0.5e-6)
        r = _result(d32=2.5e-6)  # exactly one tol above target
        obj = compute_inverse_design_objectives(r, spec, trust_aware=False)
        assert obj[0] == pytest.approx(1.0, rel=1e-9)

    def test_g_dn_uses_log10_distance(self):
        spec = TargetSpec(G_DN_target=10_000, G_DN_log10_tol=1.0)
        # G_DN=100,000 -> log10 diff = 1.0 -> normalized = 1.0
        r = _result(G_DN=100_000.0)
        obj = compute_inverse_design_objectives(r, spec, trust_aware=False)
        assert obj[0] == pytest.approx(1.0, rel=1e-9)

    def test_trust_penalty_added_per_component(self):
        spec = TargetSpec(
            d32_target=2e-6, d32_tol=0.5e-6,
            pore_target=80e-9, pore_tol=10e-9,
        )
        r = _result(
            d32=2e-6, pore=80e-9,
            tier=ModelEvidenceTier.QUALITATIVE_TREND,
        )
        obj = compute_inverse_design_objectives(r, spec, trust_aware=True)
        expected_penalty = trust_penalty_for_tier(
            ModelEvidenceTier.QUALITATIVE_TREND,
        )
        np.testing.assert_allclose(
            obj, [expected_penalty, expected_penalty], rtol=1e-9,
        )

    def test_trust_aware_false_skips_penalty(self):
        spec = TargetSpec(d32_target=2e-6, d32_tol=0.5e-6)
        r = _result(d32=2e-6, tier=ModelEvidenceTier.UNSUPPORTED)
        obj = compute_inverse_design_objectives(r, spec, trust_aware=False)
        assert obj[0] == pytest.approx(0.0, abs=1e-10)

    def test_kav_absent_gives_inf(self):
        """User asked for Kav but result has no m3 block -> +inf so the
        candidate is dropped by feasibility."""
        spec = TargetSpec(Kav_target=0.5, Kav_tol=0.05)
        r = _result()  # no m3 attribute
        obj = compute_inverse_design_objectives(r, spec, trust_aware=False)
        assert math.isinf(obj[0])

    def test_stirred_vessel_uses_d_mode_for_d_target(self):
        """In stirred-vessel mode the d32 target is compared against
        d_mode when it exists on the result."""
        spec = TargetSpec(d32_target=100e-6, d32_tol=10e-6)
        r = _result(
            d32=5e-6,        # should be ignored
            d_mode=100e-6,   # should match target -> 0
            mode="stirred_vessel",
        )
        obj = compute_inverse_design_objectives(r, spec, trust_aware=False)
        assert obj[0] == pytest.approx(0.0, abs=1e-10)

    def test_active_dims_drive_output_length(self):
        """Number of components in the objective vector = number of
        active dimensions in the TargetSpec."""
        spec = TargetSpec(
            d32_target=2e-6, d32_tol=0.5e-6,
            G_DN_target=10_000, G_DN_log10_tol=0.3,
        )
        r = _result()
        obj = compute_inverse_design_objectives(r, spec, trust_aware=False)
        assert obj.shape == (2,)
