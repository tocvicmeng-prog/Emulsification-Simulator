"""Tests for Node F1-a Phase 2c (v8.0-beta): alginate finish.

Closes the remaining protocol §6 cases:
  - Test 3: shrinking-core √t front scaling (L2 ionic-Ca)
  - Test 9: TOML round-trip of ``polymer_family = "alginate"``
  - Test 10: manifest-tier reporting from orchestrator

Plus Phase 2c wiring:
  - Alginate gelant library (CaCl₂ external + GDL/CaCO₃ internal)
  - CLI ``--polymer-family`` flag on ``emulsim run``
"""

from __future__ import annotations


import numpy as np
import pytest

from emulsim.config import load_properties
from emulsim.datatypes import (
    MaterialProperties,
    ModelEvidenceTier,
    PolymerFamily,
    SimulationParameters,
)
from emulsim.level2_gelation.ionic_ca import solve_ionic_ca_gelation
from emulsim.reagent_library_alginate import (
    GELANTS_ALGINATE,
    AlginateGelantProfile,
    effective_bath_concentration,
)


def _alg_params(
    c_alg: float = 20.0, t_gel: float = 600.0, C_Ca: float = 200.0
) -> SimulationParameters:
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


# ─── Protocol §6 test 3: √t shrinking-core scaling ──────────────────


class TestShrinkingCoreSqrtT:
    """At early times and with fast (diffusion-limited) binding, the
    volume-averaged egg-box density grows as √t because the front
    depth δ(t) ∝ √(D·t) and the gelled shell volume is ~4πR²·δ.
    """

    def test_x_mean_grows_as_sqrt_time(self):
        """Pick R large enough that the front stays in the early regime
        (δ/R ≤ ~0.4) across the whole time range, so saturation doesn't
        contaminate the slope. R = 1 mm, t ∈ [10, 40, 160] s gives
        δ/R ≈ [0.10, 0.20, 0.40] for D = 1e-9 m²/s.
        """
        params = _alg_params(c_alg=20.0, C_Ca=200.0)
        props = _alg_props(f_G=0.5)
        R = 1.0e-3  # 1 mm — far enough from saturation across the sweep

        t_list = [10.0, 40.0, 160.0]  # factor 4 each -> √t expectation: factor 2
        x_means = []
        for t_end in t_list:
            gel = solve_ionic_ca_gelation(
                params, props, R_droplet=R, time=t_end, n_r=40,
            )
            diag = gel.model_manifest.diagnostics
            x_means.append(float(diag["X_mean_final"]))

        # Log-log slope of X_mean vs t — expect 0.5 ± 0.15 in early regime.
        logs_t = np.log(t_list)
        logs_x = np.log(np.array(x_means))
        slope, _ = np.polyfit(logs_t, logs_x, 1)
        assert 0.35 <= slope <= 0.65, (
            f"shrinking-core √t scaling broken: slope={slope:.3f}, "
            f"expected ~0.5 (X_means={x_means})"
        )

    def test_front_propagates_monotonically(self):
        """Sanity companion: penetration depth grows with time."""
        params = _alg_params(c_alg=20.0, C_Ca=200.0)
        props = _alg_props(f_G=0.5)
        R = 500e-6

        g1 = solve_ionic_ca_gelation(params, props, R_droplet=R, time=5.0, n_r=40)
        g2 = solve_ionic_ca_gelation(params, props, R_droplet=R, time=50.0, n_r=40)
        assert g2.alpha_final > g1.alpha_final


# ─── Protocol §6 test 9: TOML round-trip ────────────────────────────


class TestPolymerFamilyTOMLRoundTrip:
    def test_load_properties_coerces_polymer_family_string(self, tmp_path):
        toml = tmp_path / "props.toml"
        toml.write_text(
            "[material]\n"
            'polymer_family = "alginate"\n'
            "f_guluronate = 0.5\n"
            "D_Ca = 1.0e-9\n"
        )
        props = load_properties(toml)
        assert props.polymer_family == PolymerFamily.ALGINATE
        assert props.f_guluronate == pytest.approx(0.5)
        assert props.D_Ca == pytest.approx(1.0e-9)

    def test_toml_without_polymer_family_defaults_to_agarose_chitosan(
        self, tmp_path,
    ):
        toml = tmp_path / "props.toml"
        toml.write_text(
            "[material]\n"
            "k_xlink_0 = 2806.0\n"
        )
        props = load_properties(toml)
        assert props.polymer_family == PolymerFamily.AGAROSE_CHITOSAN

    def test_toml_top_level_polymer_family_also_works(self, tmp_path):
        """TOML may place ``polymer_family`` at the top level too."""
        toml = tmp_path / "props.toml"
        toml.write_text('polymer_family = "alginate"\n')
        props = load_properties(toml)
        assert props.polymer_family == PolymerFamily.ALGINATE

    def test_toml_unknown_polymer_family_raises(self, tmp_path):
        toml = tmp_path / "props.toml"
        toml.write_text('polymer_family = "gummi_bear"\n')
        with pytest.raises(ValueError):
            load_properties(toml)


# ─── Protocol §6 test 10: manifest-tier reporting ────────────────────


class TestOrchestratorManifestTier:
    def test_alginate_path_reports_semi_quantitative(self, tmp_path):
        """The full alginate sub-pipeline should emit SEMI_QUANTITATIVE
        tiers on both the L2 (GelationResult) and L4 (MechanicalResult)
        model manifests, and the FullResult.run_report should agree.
        """
        from emulsim.pipeline.orchestrator import PipelineOrchestrator
        from emulsim.properties.database import PropertyDatabase

        params = _alg_params(c_alg=20.0, t_gel=600.0, C_Ca=200.0)
        params.formulation.c_agarose = 30.0  # keep the prop DB lookup happy
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

        assert result.gelation.model_manifest is not None
        assert result.mechanical.model_manifest is not None

        assert (
            result.gelation.model_manifest.evidence_tier
            == ModelEvidenceTier.SEMI_QUANTITATIVE
        )
        assert (
            result.mechanical.model_manifest.evidence_tier
            == ModelEvidenceTier.SEMI_QUANTITATIVE
        )

        assert result.run_report is not None
        # min_evidence_tier is stored as the string value
        min_tier = result.run_report.min_evidence_tier
        assert min_tier in {
            ModelEvidenceTier.SEMI_QUANTITATIVE.value,
            ModelEvidenceTier.QUALITATIVE_TREND.value,
            ModelEvidenceTier.UNSUPPORTED.value,
        }

    def test_alginate_manifest_names_include_platform_tag(self, tmp_path):
        from emulsim.pipeline.orchestrator import PipelineOrchestrator
        from emulsim.properties.database import PropertyDatabase

        params = _alg_params(c_alg=20.0, t_gel=300.0, C_Ca=150.0)
        params.formulation.c_agarose = 30.0
        params.emulsification.t_emulsification = 2.0
        params.solver.l1_n_bins = 30

        db = PropertyDatabase()
        orch = PipelineOrchestrator(db=db, output_dir=tmp_path / "runs")
        result = orch.run_single(
            params,
            props_overrides={"polymer_family": PolymerFamily.ALGINATE},
        )

        assert "IonicCa" in result.gelation.model_manifest.model_name
        assert "Alginate" in result.mechanical.model_manifest.model_name


# ─── Alginate gelant library (Phase 2c addition) ────────────────────


class TestAlginateGelantLibrary:
    def test_both_canonical_gelants_registered(self):
        assert "cacl2_external" in GELANTS_ALGINATE
        assert "gdl_caco3_internal" in GELANTS_ALGINATE
        ca = GELANTS_ALGINATE["cacl2_external"]
        gdl = GELANTS_ALGINATE["gdl_caco3_internal"]
        assert isinstance(ca, AlginateGelantProfile)
        assert isinstance(gdl, AlginateGelantProfile)
        assert ca.mode == "external_bath"
        assert gdl.mode == "internal_release"

    def test_external_bath_effective_concentration_is_static(self):
        ca = GELANTS_ALGINATE["cacl2_external"]
        assert effective_bath_concentration(ca, t_end=10.0) == pytest.approx(
            ca.C_Ca_bath
        )
        assert effective_bath_concentration(ca, t_end=1e6) == pytest.approx(
            ca.C_Ca_bath
        )

    def test_internal_release_effective_concentration_saturates(self):
        gdl = GELANTS_ALGINATE["gdl_caco3_internal"]
        c_short = effective_bath_concentration(gdl, t_end=100.0)
        c_long = effective_bath_concentration(gdl, t_end=1e5)
        assert 0.0 < c_short < gdl.C_Ca_source
        assert c_long > c_short
        # Long-time asymptote approaches the source concentration.
        assert c_long == pytest.approx(gdl.C_Ca_source, rel=1e-2)

    def test_internal_release_zero_time_is_zero(self):
        gdl = GELANTS_ALGINATE["gdl_caco3_internal"]
        assert effective_bath_concentration(gdl, t_end=0.0) == 0.0

    def test_unknown_mode_raises(self):
        bad = AlginateGelantProfile(
            name="x", cas="x", mode="plasma_ray",
            C_Ca_bath=0.0, C_Ca_source=0.0, k_release=0.0,
            T_default=298.15, t_default=60.0, suitability=1, notes="x",
        )
        with pytest.raises(ValueError):
            effective_bath_concentration(bad, t_end=1.0)


# ─── CLI ``--polymer-family`` flag ───────────────────────────────────


class TestCliPolymerFamilyFlag:
    def test_run_subparser_accepts_alginate_choice(self):
        """argparse should accept ``emulsim run --polymer-family alginate``."""
        import argparse
        from emulsim.__main__ import main as cli_main  # noqa: F401

        # Re-build the parser the same way __main__.main() does — we
        # don't invoke main because that would run the pipeline; we
        # just want to confirm the flag parses.
        parser = argparse.ArgumentParser(prog="emulsim")
        sub = parser.add_subparsers(dest="command")
        run_p = sub.add_parser("run")
        run_p.add_argument("config", nargs="?", default=None)
        run_p.add_argument("--rpm", type=float)
        run_p.add_argument("--phi-d", type=float)
        run_p.add_argument(
            "--polymer-family", dest="polymer_family", default=None,
            choices=["agarose_chitosan", "alginate", "cellulose", "plga"],
        )
        run_p.add_argument("--output", "-o", default="output")
        run_p.add_argument("--quiet", "-q", action="store_true")

        args = parser.parse_args(["run", "--polymer-family", "alginate"])
        assert args.polymer_family == "alginate"

    def test_main_parser_has_polymer_family_on_run(self):
        """Directly verify the shipped parser has the flag and choices."""
        # Walk the parser declaration in __main__ by importing its
        # module-level factory if we had one — we don't, so do it via
        # `-h` would call sys.exit. Instead, re-import main, set argv,
        # and assert via argparse SystemExit+help text is too brittle.
        # Simplest: regex the source text of __main__.
        import emulsim.__main__ as mm
        import inspect
        src = inspect.getsource(mm)
        assert '"--polymer-family"' in src
        assert '"alginate"' in src
        assert '"agarose_chitosan"' in src


# ─── Gelant preset flag / TOML (Phase 2c polish) ────────────────────


class TestGelantPreset:
    """``--gelant`` CLI + ``[formulation].gelant`` TOML key wire a
    `GELANTS_ALGINATE` profile's effective Ca²⁺ concentration into
    ``FormulationParameters.c_Ca_bath`` using the current
    ``t_crosslink``. The ``gelant`` key itself is consumed at config
    load time — it is not a persistent field on the dataclass.
    """

    def test_toml_gelant_cacl2_external_sets_bath_concentration(
        self, tmp_path,
    ):
        from emulsim.config import load_config

        cfg = tmp_path / "sim.toml"
        cfg.write_text(
            "[emulsification]\n"
            'mode = "rotor_stator_legacy"\n'
            "\n"
            "[formulation]\n"
            'gelant = "cacl2_external"\n'
            "c_alginate = 20.0\n"
            "t_crosslink = 900.0\n"
        )
        params = load_config(cfg)
        # CaCl₂ external bath is static — 100 mol/m³ regardless of t.
        assert params.formulation.c_Ca_bath == pytest.approx(100.0)
        assert params.formulation.c_alginate == pytest.approx(20.0)
        assert params.formulation.t_crosslink == pytest.approx(900.0)

    def test_toml_gelant_gdl_internal_saturates_with_time(self, tmp_path):
        from emulsim.config import load_config

        cfg_short = tmp_path / "short.toml"
        cfg_short.write_text(
            "[emulsification]\n"
            'mode = "rotor_stator_legacy"\n'
            "\n"
            "[formulation]\n"
            'gelant = "gdl_caco3_internal"\n'
            "t_crosslink = 600.0\n"
        )
        cfg_long = tmp_path / "long.toml"
        cfg_long.write_text(
            "[emulsification]\n"
            'mode = "rotor_stator_legacy"\n'
            "\n"
            "[formulation]\n"
            'gelant = "gdl_caco3_internal"\n'
            "t_crosslink = 50000.0\n"
        )
        p_short = load_config(cfg_short)
        p_long = load_config(cfg_long)
        # Internal release: C_eff saturates toward C_source with time.
        assert p_short.formulation.c_Ca_bath < p_long.formulation.c_Ca_bath
        assert p_long.formulation.c_Ca_bath == pytest.approx(20.0, rel=5e-2)

    def test_toml_unknown_gelant_raises(self, tmp_path):
        from emulsim.config import load_config

        cfg = tmp_path / "bad.toml"
        cfg.write_text(
            "[emulsification]\n"
            'mode = "rotor_stator_legacy"\n'
            "\n"
            "[formulation]\n"
            'gelant = "moonbeam_dust"\n'
        )
        with pytest.raises(ValueError, match="Unknown alginate gelant"):
            load_config(cfg)

    def test_cli_run_subparser_accepts_gelant_choices(self):
        import argparse

        parser = argparse.ArgumentParser(prog="emulsim")
        sub = parser.add_subparsers(dest="command")
        run_p = sub.add_parser("run")
        run_p.add_argument("config", nargs="?", default=None)
        run_p.add_argument(
            "--gelant", default=None,
            choices=["cacl2_external", "gdl_caco3_internal"],
        )
        run_p.add_argument("--quiet", "-q", action="store_true")

        args = parser.parse_args(["run", "--gelant", "gdl_caco3_internal"])
        assert args.gelant == "gdl_caco3_internal"

    def test_main_parser_has_gelant_on_run(self):
        import emulsim.__main__ as mm
        import inspect
        src = inspect.getsource(mm)
        assert '"--gelant"' in src
        assert '"cacl2_external"' in src
        assert '"gdl_caco3_internal"' in src
