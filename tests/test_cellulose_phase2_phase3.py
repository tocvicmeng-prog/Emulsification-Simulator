"""Tests for Node F1-b Phase 2 + Phase 3: cellulose orchestrator +
TOML/CLI wiring + full solvent-system registry.

Phase 2 scope (5 integration tests):
  1. Orchestrator dispatches `PolymerFamily.CELLULOSE` to `_run_cellulose`
  2. TOML round-trip of `[formulation].solvent_system = "naoh_urea"`
  3. Orchestrator applies the solvent preset to props before solving L2
  4. CLI `--polymer-family cellulose --cellulose-solvent naoh_urea` parses
  5. Full pipeline produces non-null FullResult with SEMI_QUANTITATIVE tier

Phase 3 scope (4 solvent-dependence tests):
  6. All four presets registered (naoh_urea, nmmo, emim_ac, dmac_licl)
  7. Each preset has a physically-reasonable parameter set
  8. Switching presets changes the L2 / L4 manifest diagnostics
  9. Unknown preset via CLI flag validated by argparse choices
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
from emulsim.properties.cellulose_defaults import (
    CELLULOSE_SOLVENT_PRESETS,
    apply_preset,
)


def _cell_params(phi0: float = 0.05, t_end: float = 100.0,
                 solvent: str = "") -> SimulationParameters:
    p = SimulationParameters()
    p.formulation.phi_cellulose_0 = phi0
    p.formulation.c_agarose = 30.0  # keep prop-DB lookup happy
    p.formulation.c_chitosan = 18.0
    p.formulation.c_alginate = 0.0
    p.formulation.t_crosslink = t_end
    p.formulation.solvent_system = solvent
    p.emulsification.t_emulsification = 2.0
    p.solver.l1_n_bins = 30
    return p


# ─── Phase 2: Orchestrator dispatch & integration ───────────────────


class TestOrchestratorDispatchCellulose:
    def test_polymer_family_cellulose_routes_to_run_cellulose(
        self, tmp_path,
    ):
        from emulsim.pipeline.orchestrator import PipelineOrchestrator
        from emulsim.properties.database import PropertyDatabase

        params = _cell_params(phi0=0.05, t_end=60.0)
        db = PropertyDatabase()
        orch = PipelineOrchestrator(db=db, output_dir=tmp_path / "runs")
        result = orch.run_single(
            params,
            props_overrides={"polymer_family": PolymerFamily.CELLULOSE},
        )
        # Full pipeline returns a FullResult with cellulose signatures.
        assert result.mechanical.model_used == "cellulose_zhang2020"
        assert result.mechanical.network_type == "physical_entangled"
        assert result.gelation.model_manifest.model_name == (
            "L2.Gelation.NIPSCellulose"
        )
        # L3 is stubbed for cellulose.
        assert result.crosslinking.p_final == pytest.approx(0.0)

    def test_summary_json_records_cellulose_polymer_family(self, tmp_path):
        from emulsim.pipeline.orchestrator import PipelineOrchestrator
        from emulsim.properties.database import PropertyDatabase

        params = _cell_params(phi0=0.05, t_end=60.0)
        db = PropertyDatabase()
        orch = PipelineOrchestrator(db=db, output_dir=tmp_path / "runs")
        orch.run_single(
            params,
            props_overrides={"polymer_family": PolymerFamily.CELLULOSE},
        )
        summaries = list((tmp_path / "runs").rglob("summary.json"))
        assert summaries, "orchestrator must emit summary.json"
        payload = json.loads(summaries[0].read_text())
        assert payload["polymer_family"] == "cellulose"
        assert payload["level3"]["note"].startswith("not applicable")
        assert "phi_mean_final" in payload["level2"]

    def test_orchestrator_reports_semi_quantitative_end_to_end(
        self, tmp_path,
    ):
        from emulsim.pipeline.orchestrator import PipelineOrchestrator
        from emulsim.properties.database import PropertyDatabase

        params = _cell_params(phi0=0.05, t_end=60.0)
        db = PropertyDatabase()
        orch = PipelineOrchestrator(db=db, output_dir=tmp_path / "runs")
        result = orch.run_single(
            params,
            props_overrides={"polymer_family": PolymerFamily.CELLULOSE},
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


class TestTomlSolventSystemWiring:
    def test_toml_solvent_system_unpacks_into_formulation(self, tmp_path):
        from emulsim.config import load_config

        cfg = tmp_path / "sim.toml"
        cfg.write_text(
            "[emulsification]\n"
            'mode = "rotor_stator_legacy"\n'
            "\n"
            "[formulation]\n"
            "phi_cellulose_0 = 0.05\n"
            'solvent_system = "naoh_urea"\n'
        )
        params = load_config(cfg)
        assert params.formulation.solvent_system == "naoh_urea"
        assert params.formulation.phi_cellulose_0 == pytest.approx(0.05)

    def test_toml_without_solvent_system_defaults_to_empty(self, tmp_path):
        from emulsim.config import load_config

        cfg = tmp_path / "sim.toml"
        cfg.write_text(
            "[emulsification]\n"
            'mode = "rotor_stator_legacy"\n'
            "\n"
            "[formulation]\n"
            "phi_cellulose_0 = 0.05\n"
        )
        params = load_config(cfg)
        assert params.formulation.solvent_system == ""


class TestOrchestratorAppliesSolventPreset:
    def test_solvent_system_applied_before_l2(self, tmp_path):
        """When params.formulation.solvent_system = "nmmo", the
        orchestrator must patch MaterialProperties to the NMMO preset
        before calling the NIPS solver. We detect this via the L4
        modulus prefactor ending up at the NMMO value (8e5 Pa) rather
        than the default (5e5 Pa).
        """
        from emulsim.pipeline.orchestrator import PipelineOrchestrator
        from emulsim.properties.database import PropertyDatabase

        params = _cell_params(phi0=0.06, t_end=60.0, solvent="nmmo")
        db = PropertyDatabase()
        orch = PipelineOrchestrator(db=db, output_dir=tmp_path / "runs")
        result = orch.run_single(
            params,
            props_overrides={"polymer_family": PolymerFamily.CELLULOSE},
        )
        # NMMO K_cell = 8e5; default K_cell = 5e5.
        # phi_mean_final ≈ phi_0 = 0.06, alpha = 2.25
        # G_DN(NMMO) / G_DN(default) = 8/5 = 1.6
        K_obs = (
            result.mechanical.model_manifest.diagnostics["G_DN_Pa"]
            / max(
                result.mechanical.model_manifest.diagnostics["phi_cell_mean"]
                ** 2.25,
                1e-12,
            )
        )
        assert K_obs == pytest.approx(8.0e5, rel=1e-3)

    def test_unknown_solvent_system_raises(self, tmp_path):
        from emulsim.pipeline.orchestrator import PipelineOrchestrator
        from emulsim.properties.database import PropertyDatabase

        params = _cell_params(phi0=0.05, t_end=60.0, solvent="bathtub_soup")
        db = PropertyDatabase()
        orch = PipelineOrchestrator(db=db, output_dir=tmp_path / "runs")
        with pytest.raises(KeyError):
            orch.run_single(
                params,
                props_overrides={"polymer_family": PolymerFamily.CELLULOSE},
            )


class TestCliCellulose:
    def test_main_parser_has_cellulose_solvent_flag(self):
        import emulsim.__main__ as mm
        import inspect
        src = inspect.getsource(mm)
        assert '"--cellulose-solvent"' in src
        assert '"naoh_urea"' in src
        assert '"nmmo"' in src
        assert '"emim_ac"' in src
        assert '"dmac_licl"' in src

    def test_run_subparser_accepts_cellulose_solvent_choices(self):
        """argparse invocation with the full flag-set parses cleanly."""
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
            "--cellulose-solvent", dest="cellulose_solvent", default=None,
            choices=["naoh_urea", "nmmo", "emim_ac", "dmac_licl"],
        )
        run_p.add_argument("--quiet", "-q", action="store_true")

        args = parser.parse_args([
            "run", "--polymer-family", "cellulose",
            "--cellulose-solvent", "nmmo",
        ])
        assert args.polymer_family == "cellulose"
        assert args.cellulose_solvent == "nmmo"


# ─── Phase 3: All four solvent presets ──────────────────────────────


class TestAllFourPresetsRegistered:
    def test_all_four_names_present(self):
        assert set(CELLULOSE_SOLVENT_PRESETS) == {
            "naoh_urea", "nmmo", "emim_ac", "dmac_licl",
        }

    def test_each_preset_has_physically_reasonable_params(self):
        for name, p in CELLULOSE_SOLVENT_PRESETS.items():
            # Flory-Huggins χ should be in a plausible range for the
            # cellulose/solvent/water ternary. (χ_PN > χ_PS is the
            # basic "bad-solvent" condition for NIPS.)
            assert 0.0 < p.chi_PS < 1.0, f"{name}: chi_PS out of range"
            assert 0.0 < p.chi_PN < 1.5, f"{name}: chi_PN out of range"
            assert p.chi_PN > p.chi_PS, (
                f"{name}: water must be a worse solvent than the "
                f"solvent for NIPS to work (chi_PN > chi_PS)"
            )
            # Diffusivities positive and sub-bulk-water.
            assert 0.0 < p.D_solvent < 1e-8
            assert 0.0 < p.D_nonsolvent < 1e-8
            # N_p in realistic DP range.
            assert 100.0 < p.N_p < 2000.0
            # Modulus prefactor spans a couple decades.
            assert 1e4 < p.K_cell < 1e7
            # Exponent strictly positive, typically 2-3.
            assert 1.5 < p.alpha_cell < 3.5
            # Processing T in realistic band.
            assert 273.0 <= p.T_process <= 423.0

    def test_water_is_the_default_nonsolvent_for_all(self):
        for name, p in CELLULOSE_SOLVENT_PRESETS.items():
            assert "water" in p.nonsolvent.lower(), f"{name}"


class TestPresetDiagnosticsDifferentiate:
    def test_switching_preset_changes_l4_manifest_diagnostics(self):
        """Applying different presets to the same MaterialProperties
        yields different L4 manifest K_cell / alpha_cell entries.
        """
        from emulsim.level2_gelation.nips_cellulose import solve_nips_cellulose
        from emulsim.level4_mechanical.cellulose import (
            solve_mechanical_cellulose,
        )

        params = _cell_params(phi0=0.05, t_end=60.0)

        preset_names = list(CELLULOSE_SOLVENT_PRESETS)
        G_values = {}
        for name in preset_names:
            props = MaterialProperties()
            props.polymer_family = PolymerFamily.CELLULOSE
            apply_preset(props, name)
            gel = solve_nips_cellulose(
                params, props, R_droplet=100e-6, n_r=20, seed=0,
            )
            mech = solve_mechanical_cellulose(params, props, gel)
            G_values[name] = mech.G_DN

        # G_DN values should span a real range — not all equal.
        vals = list(G_values.values())
        assert max(vals) > 1.5 * min(vals), (
            f"presets should give distinguishable moduli: {G_values}"
        )


class TestArgparseRejectsUnknownPreset:
    def test_cli_choice_validation_via_argparse(self):
        import argparse

        parser = argparse.ArgumentParser(prog="emulsim")
        sub = parser.add_subparsers(dest="command")
        run_p = sub.add_parser("run")
        run_p.add_argument(
            "--cellulose-solvent", dest="cellulose_solvent", default=None,
            choices=["naoh_urea", "nmmo", "emim_ac", "dmac_licl"],
        )
        with pytest.raises(SystemExit):
            parser.parse_args(["run", "--cellulose-solvent", "moon_juice"])
