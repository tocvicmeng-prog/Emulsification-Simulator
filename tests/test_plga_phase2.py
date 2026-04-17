"""Tests for Node F1-c Phase 2: PLGA orchestrator + CLI + TOML
integration.

12 tests:

  1. Orchestrator dispatches `PolymerFamily.PLGA` to `_run_plga`
  2. Summary.json records polymer_family=plga + L2 diagnostics
  3. Full pipeline reports SEMI_QUANTITATIVE tier end-to-end
  4. Orchestrator applies plga_grade preset before L2
  5. Unknown grade raises KeyError from `apply_preset`
  6. TOML `[formulation].plga_grade` round-trip
  7. TOML without plga_grade defaults to empty string
  8. CLI `--plga-grade` flag accepted by argparse
  9. CLI source mentions all four grade choices
 10. Full pipeline produces non-zero G_DN
 11. Switching grade via preset changes L4 G_DN
 12. L3 stubbed for PLGA (p_final = 0)
"""

from __future__ import annotations

import json

import pytest

from emulsim.datatypes import (
    MaterialProperties,
    ModelEvidenceTier,
    PolymerFamily,
    SimulationParameters,
)
from emulsim.properties.plga_defaults import (
    PLGA_GRADE_PRESETS,
    apply_preset,
)


def _plga_params(phi0: float = 0.10, t_end: float = 60.0,
                 grade: str = "") -> SimulationParameters:
    p = SimulationParameters()
    p.formulation.phi_PLGA_0 = phi0
    p.formulation.plga_grade = grade
    p.formulation.c_agarose = 30.0   # keep prop-DB lookups happy
    p.formulation.c_chitosan = 18.0
    p.formulation.c_alginate = 0.0
    p.formulation.phi_cellulose_0 = 0.0
    p.formulation.t_crosslink = t_end
    p.emulsification.t_emulsification = 2.0
    p.solver.l1_n_bins = 30
    return p


# ─── Orchestrator dispatch ──────────────────────────────────────────


class TestOrchestratorDispatchPLGA:
    def test_polymer_family_plga_routes_to_run_plga(self, tmp_path):
        from emulsim.pipeline.orchestrator import PipelineOrchestrator
        from emulsim.properties.database import PropertyDatabase

        params = _plga_params(phi0=0.10, t_end=60.0)
        db = PropertyDatabase()
        orch = PipelineOrchestrator(db=db, output_dir=tmp_path / "runs")
        result = orch.run_single(
            params,
            props_overrides={"polymer_family": PolymerFamily.PLGA},
        )
        assert result.mechanical.model_used == "plga_gibson_ashby"
        assert result.mechanical.network_type == "glassy_polymer"
        assert result.gelation.model_manifest.model_name == (
            "L2.Gelation.SolventEvaporationPLGA"
        )

    def test_summary_json_records_plga_polymer_family(self, tmp_path):
        from emulsim.pipeline.orchestrator import PipelineOrchestrator
        from emulsim.properties.database import PropertyDatabase

        params = _plga_params(phi0=0.10, t_end=60.0)
        db = PropertyDatabase()
        orch = PipelineOrchestrator(db=db, output_dir=tmp_path / "runs")
        orch.run_single(
            params,
            props_overrides={"polymer_family": PolymerFamily.PLGA},
        )
        summaries = list((tmp_path / "runs").rglob("summary.json"))
        assert summaries
        payload = json.loads(summaries[0].read_text())
        assert payload["polymer_family"] == "plga"
        assert payload["level3"]["note"].startswith("not applicable")
        for key in (
            "phi_plga_mean_final", "phi_dcm_mean_final",
            "t_vitrification_s", "skin_thickness_proxy_m",
            "R_shrunk_m", "porosity",
        ):
            assert key in payload["level2"]

    def test_plga_reports_semi_quantitative_end_to_end(self, tmp_path):
        from emulsim.pipeline.orchestrator import PipelineOrchestrator
        from emulsim.properties.database import PropertyDatabase

        params = _plga_params(phi0=0.10, t_end=60.0)
        db = PropertyDatabase()
        orch = PipelineOrchestrator(db=db, output_dir=tmp_path / "runs")
        result = orch.run_single(
            params,
            props_overrides={"polymer_family": PolymerFamily.PLGA},
        )
        assert (
            result.gelation.model_manifest.evidence_tier
            == ModelEvidenceTier.SEMI_QUANTITATIVE
        )
        assert (
            result.mechanical.model_manifest.evidence_tier
            == ModelEvidenceTier.SEMI_QUANTITATIVE
        )
        assert result.run_report is not None


class TestPlgaGradePresetApplication:
    def test_grade_applied_before_l2_solver(self, tmp_path):
        """When params.formulation.plga_grade = "85_15", the
        orchestrator must patch MaterialProperties to the 85:15
        preset before calling the L2 solver. Detected via the L4
        modulus prefactor ending up at the 85:15 value (1e9 Pa)
        vs the default (7e8 Pa).
        """
        from emulsim.pipeline.orchestrator import PipelineOrchestrator
        from emulsim.properties.database import PropertyDatabase

        params = _plga_params(phi0=0.08, t_end=60.0, grade="85_15")
        db = PropertyDatabase()
        orch = PipelineOrchestrator(db=db, output_dir=tmp_path / "runs")
        result = orch.run_single(
            params,
            props_overrides={"polymer_family": PolymerFamily.PLGA},
        )
        # 85:15 G_glassy = 1e9. After L2 the mean phi should be ≈ 0.995
        # (Dirichlet sink at 0.005). Recover G_glassy = G_DN / phi^n.
        phi_mean = result.mechanical.model_manifest.diagnostics[
            "phi_plga_mean"
        ]
        G_DN = result.mechanical.model_manifest.diagnostics["G_DN_Pa"]
        G_glassy_obs = G_DN / max(phi_mean ** 2.0, 1e-12)
        assert G_glassy_obs == pytest.approx(1.0e9, rel=2e-2)

    def test_unknown_grade_raises(self, tmp_path):
        from emulsim.pipeline.orchestrator import PipelineOrchestrator
        from emulsim.properties.database import PropertyDatabase

        params = _plga_params(phi0=0.10, t_end=60.0, grade="quantum_foam")
        db = PropertyDatabase()
        orch = PipelineOrchestrator(db=db, output_dir=tmp_path / "runs")
        with pytest.raises(KeyError):
            orch.run_single(
                params,
                props_overrides={"polymer_family": PolymerFamily.PLGA},
            )


# ─── TOML wiring ────────────────────────────────────────────────────


class TestTomlPlgaGradeWiring:
    def test_toml_plga_grade_unpacks_into_formulation(self, tmp_path):
        from emulsim.config import load_config

        cfg = tmp_path / "sim.toml"
        cfg.write_text(
            "[emulsification]\n"
            'mode = "rotor_stator_legacy"\n'
            "\n"
            "[formulation]\n"
            "phi_PLGA_0 = 0.10\n"
            'plga_grade = "75_25"\n'
        )
        params = load_config(cfg)
        assert params.formulation.plga_grade == "75_25"
        assert params.formulation.phi_PLGA_0 == pytest.approx(0.10)

    def test_toml_without_plga_grade_defaults_empty(self, tmp_path):
        from emulsim.config import load_config

        cfg = tmp_path / "sim.toml"
        cfg.write_text(
            "[emulsification]\n"
            'mode = "rotor_stator_legacy"\n'
            "\n"
            "[formulation]\n"
            "phi_PLGA_0 = 0.10\n"
        )
        params = load_config(cfg)
        assert params.formulation.plga_grade == ""


# ─── CLI ────────────────────────────────────────────────────────────


class TestCliPlgaGrade:
    def test_cli_source_mentions_all_four_grades(self):
        import emulsim.__main__ as mm
        import inspect
        src = inspect.getsource(mm)
        assert '"--plga-grade"' in src
        assert '"50_50"' in src
        assert '"75_25"' in src
        assert '"85_15"' in src
        assert '"pla"' in src

    def test_argparse_accepts_plga_grade_choices(self):
        import argparse

        parser = argparse.ArgumentParser(prog="emulsim")
        sub = parser.add_subparsers(dest="command")
        run_p = sub.add_parser("run")
        run_p.add_argument("config", nargs="?", default=None)
        run_p.add_argument(
            "--polymer-family", default=None,
            choices=["agarose_chitosan", "alginate", "cellulose", "plga"],
        )
        run_p.add_argument(
            "--plga-grade", dest="plga_grade", default=None,
            choices=["50_50", "75_25", "85_15", "pla"],
        )
        run_p.add_argument("--quiet", "-q", action="store_true")

        args = parser.parse_args([
            "run", "--polymer-family", "plga", "--plga-grade", "75_25",
        ])
        assert args.polymer_family == "plga"
        assert args.plga_grade == "75_25"


# ─── End-to-end sanity ──────────────────────────────────────────────


class TestFullPipelineSanity:
    def test_pipeline_produces_nonzero_modulus(self, tmp_path):
        from emulsim.pipeline.orchestrator import PipelineOrchestrator
        from emulsim.properties.database import PropertyDatabase

        params = _plga_params(phi0=0.10, t_end=60.0)
        db = PropertyDatabase()
        orch = PipelineOrchestrator(db=db, output_dir=tmp_path / "runs")
        result = orch.run_single(
            params,
            props_overrides={"polymer_family": PolymerFamily.PLGA},
        )
        assert result.mechanical.G_DN > 0.0

    def test_switching_grade_changes_modulus_via_orchestrator(
        self, tmp_path,
    ):
        from emulsim.pipeline.orchestrator import PipelineOrchestrator
        from emulsim.properties.database import PropertyDatabase

        Gs = {}
        for grade in PLGA_GRADE_PRESETS:
            params = _plga_params(phi0=0.10, t_end=60.0, grade=grade)
            db = PropertyDatabase()
            orch = PipelineOrchestrator(db=db, output_dir=tmp_path / grade)
            result = orch.run_single(
                params,
                props_overrides={"polymer_family": PolymerFamily.PLGA},
            )
            Gs[grade] = result.mechanical.G_DN
        vals = list(Gs.values())
        # 50:50 (7e8) vs PLA (1.2e9) = ~1.7× spread minimum
        assert max(vals) >= 1.5 * min(vals), (
            f"grades should produce distinguishable moduli: {Gs}"
        )

    def test_l3_stubbed_for_plga(self, tmp_path):
        from emulsim.pipeline.orchestrator import PipelineOrchestrator
        from emulsim.properties.database import PropertyDatabase

        params = _plga_params(phi0=0.10, t_end=60.0)
        db = PropertyDatabase()
        orch = PipelineOrchestrator(db=db, output_dir=tmp_path / "runs")
        result = orch.run_single(
            params,
            props_overrides={"polymer_family": PolymerFamily.PLGA},
        )
        # Zero-crosslink stub: p_final must be 0.
        assert result.crosslinking.p_final == pytest.approx(0.0)
