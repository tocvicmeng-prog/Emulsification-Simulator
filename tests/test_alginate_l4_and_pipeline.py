"""Tests for Node F1-a Phase 2b (v8.0): alginate L4 modulus + orchestrator
dispatch.

Closes protocol §6 tests 6, 7, 11:
  - Modulus scaling in c_alginate (test_l4_modulus_scales_with_c_sq)
  - Modulus scaling in f_guluronate (test_l4_modulus_scales_with_fg_sq)
  - Full-pipeline integration with PolymerFamily.ALGINATE
    (test_orchestrator_dispatch_runs_alginate)
"""

from __future__ import annotations

import numpy as np
import pytest

from emulsim.datatypes import (
    FullResult,
    MaterialProperties,
    MechanicalResult,
    ModelEvidenceTier,
    PolymerFamily,
    SimulationParameters,
)
from emulsim.level2_gelation.ionic_ca import solve_ionic_ca_gelation
from emulsim.level4_mechanical.alginate import (
    alginate_modulus,
    solve_mechanical_alginate,
)


def _alg_params(c_alg: float = 20.0, t_gel: float = 600.0,
                 C_Ca: float = 200.0) -> SimulationParameters:
    p = SimulationParameters()
    p.formulation.c_alginate = c_alg
    p.formulation.c_agarose = 0.0
    p.formulation.c_chitosan = 0.0
    p.formulation.c_Ca_bath = C_Ca
    p.formulation.t_crosslink = t_gel
    return p


def _alg_props(f_G: float = 0.5) -> MaterialProperties:
    props = MaterialProperties()
    props.polymer_family = PolymerFamily.ALGINATE
    props.f_guluronate = f_G
    return props


# ─── alginate_modulus unit tests ──────────────────────────────────────────


class TestAlginateModulusFn:
    def test_zero_inputs_give_zero(self):
        assert alginate_modulus(0.0, 0.5, 50.0, 30e3, 2.0) == 0.0
        assert alginate_modulus(20.0, 0.0, 50.0, 30e3, 2.0) == 0.0
        assert alginate_modulus(20.0, 0.5, 0.0, 30e3, 2.0) == 0.0

    def test_modulus_scales_as_c_alginate_squared(self):
        """For fixed X_mean/X_max=1, G ∝ (c·f_G)^n_alg=2 ⇒ doubling c → 4×."""
        G0 = _guluronate(20.0, 0.5)
        X_full = 0.5 * G0  # stoichiometric ceiling
        m1 = alginate_modulus(20.0, 0.5, X_full, 30e3, 2.0)

        G0_2 = _guluronate(40.0, 0.5)
        X_full_2 = 0.5 * G0_2
        m2 = alginate_modulus(40.0, 0.5, X_full_2, 30e3, 2.0)

        assert m2 == pytest.approx(4.0 * m1, rel=1e-6)

    def test_modulus_scales_as_fg_squared(self):
        """Fixed c; doubling f_G → 4× modulus (conversion stays 1.0)."""
        G0 = _guluronate(20.0, 0.3)
        m1 = alginate_modulus(20.0, 0.3, 0.5 * G0, 30e3, 2.0)

        G0_2 = _guluronate(20.0, 0.6)
        m2 = alginate_modulus(20.0, 0.6, 0.5 * G0_2, 30e3, 2.0)

        assert m2 == pytest.approx(4.0 * m1, rel=1e-6)

    def test_incomplete_gelation_reduces_modulus(self):
        """X_mean at half the stoichiometric ceiling halves the modulus."""
        G0 = _guluronate(20.0, 0.5)
        m_full = alginate_modulus(20.0, 0.5, 0.5 * G0, 30e3, 2.0)
        m_half = alginate_modulus(20.0, 0.5, 0.25 * G0, 30e3, 2.0)
        assert m_half == pytest.approx(0.5 * m_full, rel=1e-6)


def _guluronate(c_alg: float, f_G: float) -> float:
    """Helper to mirror the L2 solver's concentration calculation."""
    return (c_alg / 0.198) * f_G


# ─── solve_mechanical_alginate integration ────────────────────────────────


class TestSolveMechanicalAlginate:
    def test_returns_mechanical_result_with_alginate_manifest(self):
        params = _alg_params(c_alg=20.0, t_gel=600.0, C_Ca=200.0)
        props = _alg_props()
        gel = solve_ionic_ca_gelation(params, props, R_droplet=40e-6)
        mech = solve_mechanical_alginate(params, props, gel,
                                           R_droplet=40e-6)
        assert isinstance(mech, MechanicalResult)
        assert mech.model_used == "alginate_kong2004"
        assert mech.network_type == "ionic_reinforced"
        assert mech.G_agarose == 0.0
        assert mech.G_chitosan == 0.0
        assert mech.G_DN > 0.0
        assert mech.E_star > 0.0
        assert mech.model_manifest is not None
        assert (mech.model_manifest.model_name
                == "L4.Mechanical.AlginateKong2004")
        assert (mech.model_manifest.evidence_tier
                == ModelEvidenceTier.SEMI_QUANTITATIVE)

    def test_zero_alginate_gives_zero_modulus(self):
        params = _alg_params(c_alg=0.0)
        props = _alg_props()
        gel = solve_ionic_ca_gelation(params, props, R_droplet=40e-6)
        mech = solve_mechanical_alginate(params, props, gel,
                                           R_droplet=40e-6)
        assert mech.G_DN == 0.0


# ─── Orchestrator dispatch ────────────────────────────────────────────────


class TestOrchestratorDispatch:
    def test_orchestrator_dispatches_to_alginate(self, tmp_path):
        """PolymerFamily.ALGINATE routes run_single through the
        alginate sub-pipeline; the returned FullResult has the
        alginate L4 manifest and a stubbed L3 (no covalent
        crosslinking for pure alginate)."""
        from emulsim.pipeline.orchestrator import PipelineOrchestrator
        from emulsim.properties.database import PropertyDatabase

        # Build an alginate params / props pair. The PropertyDatabase
        # update_for_conditions returns a chitosan/agarose
        # MaterialProperties; we must force polymer_family back to
        # ALGINATE via props_overrides.
        params = _alg_params(c_alg=20.0, t_gel=600.0, C_Ca=200.0)
        # Set a moderate agarose value so the property DB doesn't
        # fail upstream lookups; we override to ALGINATE below.
        params.formulation.c_agarose = 30.0
        # Short emulsification so the test runs quickly.
        params.emulsification.t_emulsification = 2.0
        params.solver.l1_n_bins = 30

        db = PropertyDatabase()
        orch = PipelineOrchestrator(db=db, output_dir=tmp_path / "runs")

        result = orch.run_single(
            params,
            props_overrides={
                "polymer_family": PolymerFamily.ALGINATE,
                "f_guluronate": 0.5,
            },
        )

        assert isinstance(result, FullResult)
        assert result.mechanical.model_used == "alginate_kong2004"
        assert result.mechanical.network_type == "ionic_reinforced"
        # L3 stubbed, so p_final = 0.
        assert result.crosslinking.p_final == pytest.approx(0.0)
        # L2 used the ionic solver, so alpha_final reflects Ca²⁺
        # conversion (non-zero for C_Ca=200 mM).
        assert result.gelation.alpha_final > 0.0
        # Summary JSON should record polymer_family
        import json
        summaries = list((tmp_path / "runs").rglob("summary.json"))
        assert summaries, "orchestrator must emit summary.json"
        payload = json.loads(summaries[0].read_text())
        assert payload["polymer_family"] == "alginate"
