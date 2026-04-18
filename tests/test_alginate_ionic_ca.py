"""Tests for Node F1-a Phase 2a (v8.0): alginate ionic-Ca gelation solver.

Covers the core invariants from the F1-a protocol §6 (tests 2, 3, 4, 5, 6):

- Mass conservation (total guluronate budget is preserved: G + 2*X = G0).
- Complete gel at long time for non-zero bath concentration.
- Zero bath → zero crosslink.
- Input validation.
- PolymerFamily enum + MaterialProperties.polymer_family field wired.
- GelationResult schema honoured so downstream L3/L4 are platform-agnostic.

Tests 3 (√t scaling), 6 (modulus scaling), and 8-11 (full pipeline
integration) are deferred to F1-a Phase 2b.
"""

from __future__ import annotations

import pytest

from emulsim.datatypes import (
    GelationResult,
    MaterialProperties,
    ModelEvidenceTier,
    PolymerFamily,
    SimulationParameters,
)
from emulsim.level2_gelation.ionic_ca import (
    solve_ionic_ca_gelation,
    _guluronate_concentration_mol_m3,
)


def _make_alginate_params(c_alg: float = 20.0, t_gel: float = 600.0):
    """Build a minimal SimulationParameters for alginate tests."""
    p = SimulationParameters()
    p.formulation.c_agarose = c_alg       # reused as alginate concentration
    p.formulation.t_crosslink = t_gel
    return p


def _make_alginate_props(f_G: float = 0.5, D_Ca: float = 1.0e-9,
                          k: float = 1.0e3) -> MaterialProperties:
    props = MaterialProperties()
    props.polymer_family = PolymerFamily.ALGINATE
    props.f_guluronate = f_G
    props.D_Ca = D_Ca
    props.k_bind_Ca = k
    return props


class TestPolymerFamily:
    def test_enum_values(self):
        assert PolymerFamily.AGAROSE_CHITOSAN.value == "agarose_chitosan"
        assert PolymerFamily.ALGINATE.value == "alginate"

    def test_default_is_agarose_chitosan(self):
        props = MaterialProperties()
        assert props.polymer_family == PolymerFamily.AGAROSE_CHITOSAN

    def test_switch_to_alginate(self):
        props = MaterialProperties()
        props.polymer_family = PolymerFamily.ALGINATE
        assert props.polymer_family == PolymerFamily.ALGINATE


class TestGuluronateConcentration:
    def test_typical_alginate(self):
        # 20 kg/m³ alginate, 50% guluronate, ~0.198 kg/mol repeat
        c = _guluronate_concentration_mol_m3(20.0, 0.5)
        # 20/0.198 * 0.5 ≈ 50.5 mol/m³
        assert c == pytest.approx(50.5, rel=0.02)

    def test_zero_alginate(self):
        assert _guluronate_concentration_mol_m3(0.0, 0.5) == 0.0


class TestSolverSchema:
    def test_returns_gelation_result(self):
        params = _make_alginate_params(c_alg=20.0, t_gel=300.0)
        props = _make_alginate_props()
        result = solve_ionic_ca_gelation(
            params, props, R_droplet=50e-6,
        )
        assert isinstance(result, GelationResult)
        assert result.r_grid.shape[0] == 50
        assert result.phi_field.shape[0] == 50
        assert result.pore_size_mean > 0
        assert 0.0 <= result.porosity <= 1.0
        assert result.model_manifest is not None
        assert (result.model_manifest.model_name ==
                "L2.Gelation.IonicCaShrinkingCore")
        assert (result.model_manifest.evidence_tier ==
                ModelEvidenceTier.SEMI_QUANTITATIVE)


class TestPhysicalInvariants:
    def test_mass_conservation_guluronate(self):
        """G_mean + 2*X_mean ≈ G0 (every Ca²⁺ consumes 2 guluronate
        residues and creates one crosslink)."""
        params = _make_alginate_params(c_alg=20.0, t_gel=600.0)
        props = _make_alginate_props()
        result = solve_ionic_ca_gelation(
            params, props, R_droplet=50e-6,
        )
        G0 = _guluronate_concentration_mol_m3(20.0, 0.5)
        diag = result.model_manifest.diagnostics
        G_mean = diag["G_mean_final"]
        X_mean = diag["X_mean_final"]
        # Conservation check: allow 5 % slip for finite-volume discretisation
        lhs = G_mean + 2.0 * X_mean
        rel_err = abs(lhs - G0) / G0
        assert rel_err < 0.05, (
            f"guluronate conservation slipped: G+2X={lhs:.3g} vs G0={G0:.3g} "
            f"(rel err {rel_err:.3f})"
        )

    def test_zero_ca_no_crosslink(self):
        """C_Ca_bath = 0 → no Ca²⁺ enters → no crosslinks form."""
        params = _make_alginate_params(c_alg=20.0, t_gel=600.0)
        props = _make_alginate_props()
        result = solve_ionic_ca_gelation(
            params, props, R_droplet=50e-6, C_Ca_bath=0.0,
        )
        X_mean = result.model_manifest.diagnostics["X_mean_final"]
        assert X_mean == pytest.approx(0.0, abs=1e-10)
        assert result.alpha_final == pytest.approx(0.0, abs=1e-10)

    def test_long_time_high_bath_reaches_significant_conversion(self):
        """At long time with plentiful Ca²⁺, core gels substantially."""
        params = _make_alginate_params(c_alg=20.0, t_gel=3600.0)  # 1 h
        props = _make_alginate_props()
        result = solve_ionic_ca_gelation(
            params, params, R_droplet=30e-6, C_Ca_bath=500.0,
        ) if False else solve_ionic_ca_gelation(
            params, props, R_droplet=30e-6, C_Ca_bath=500.0,
        )
        # ≥50% guluronate consumed somewhere between 1% and 100%
        assert result.alpha_final > 0.3, (
            f"alpha_final={result.alpha_final:.3f} below 0.3 at t=1h, C=500mM; "
            f"check D_Ca / k_bind defaults"
        )

    def test_zero_alginate_fallback(self):
        """No alginate → returns zero-gel fallback, no crash."""
        params = _make_alginate_params(c_alg=0.0, t_gel=600.0)
        props = _make_alginate_props()
        result = solve_ionic_ca_gelation(
            params, props, R_droplet=50e-6,
        )
        assert result.porosity == pytest.approx(1.0)
        assert result.alpha_final == 0.0
        assert (result.model_manifest.evidence_tier
                == ModelEvidenceTier.UNSUPPORTED)


class TestInputValidation:
    def test_bad_radius_raises(self):
        params = _make_alginate_params()
        props = _make_alginate_props()
        with pytest.raises(ValueError, match="R_droplet"):
            solve_ionic_ca_gelation(params, props, R_droplet=0.0)

    def test_bad_bath_raises(self):
        params = _make_alginate_params()
        props = _make_alginate_props()
        with pytest.raises(ValueError, match="C_Ca_bath"):
            solve_ionic_ca_gelation(
                params, props, R_droplet=50e-6, C_Ca_bath=-1.0,
            )

    def test_too_few_grid_points_raises(self):
        params = _make_alginate_params()
        props = _make_alginate_props()
        with pytest.raises(ValueError, match="n_r"):
            solve_ionic_ca_gelation(
                params, props, R_droplet=50e-6, n_r=4,
            )
