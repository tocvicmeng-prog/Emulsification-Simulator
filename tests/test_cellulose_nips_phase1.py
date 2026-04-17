"""Tests for Node F1-b Phase 1 (v8.1-alpha): L2 NIPS solver + L4 modulus
+ NaOH/urea defaults.

Core 8 tests (per docs/f1b_cellulose_nips_protocol.md §6; Phase 1
subset — orchestrator dispatch / TOML integration land in Phase 2):

  1. Ternary mass conservation (phi + s + n = 1 everywhere, all times)
  2. Water-bath driven demixing (phi_std increases from initial)
  3. Zero cellulose → zero modulus + zero gel (both solvers agree)
  4. Modulus scaling (doubling phi → G scales as phi^alpha)
  5. Dispatch-ready manifest schema (SEMI_QUANTITATIVE, correct names)
  6. NaOH/urea preset patches MaterialProperties correctly
  7. L4 modulus function zero-input sanity
  8. Solver input validation (R_droplet <= 0, n_r < 8, phi_0 out of [0,1])
"""

from __future__ import annotations

import numpy as np
import pytest

from emulsim.datatypes import (
    MaterialProperties,
    ModelEvidenceTier,
    PolymerFamily,
    SimulationParameters,
)
from emulsim.level2_gelation.nips_cellulose import solve_nips_cellulose
from emulsim.level4_mechanical.cellulose import (
    cellulose_modulus,
    solve_mechanical_cellulose,
)
from emulsim.properties.cellulose_defaults import (
    CELLULOSE_SOLVENT_PRESETS,
    NAOH_UREA,
    apply_preset,
)


def _cell_params(phi0: float = 0.05, t_end: float = 300.0) -> SimulationParameters:
    p = SimulationParameters()
    p.formulation.phi_cellulose_0 = phi0
    p.formulation.c_agarose = 0.0
    p.formulation.c_chitosan = 0.0
    p.formulation.c_alginate = 0.0
    p.formulation.t_crosslink = t_end
    return p


def _cell_props() -> MaterialProperties:
    props = MaterialProperties()
    props.polymer_family = PolymerFamily.CELLULOSE
    apply_preset(props, "naoh_urea")
    return props


# ─── 1. Mass conservation ───────────────────────────────────────────


class TestTernaryMassConservation:
    def test_phi_plus_s_plus_n_equals_one_at_final_time(self):
        params = _cell_params(phi0=0.05, t_end=60.0)
        props = _cell_props()
        gel = solve_nips_cellulose(
            params, props, R_droplet=100e-6, n_r=20, seed=1,
        )
        diag = gel.model_manifest.diagnostics
        total = (
            diag["phi_mean_final"]
            + diag["s_mean_final"]
            + diag["n_mean_final"]
        )
        assert total == pytest.approx(1.0, rel=1e-2, abs=2e-2)


# ─── 2. Water-bath driven demixing ──────────────────────────────────


class TestDemixing:
    def test_phi_std_grows_when_water_bath_pulls_solvent_out(self):
        """With phi_bath = s_bath = 0 (pure water bath), the solvent
        is driven out of the droplet; cellulose concentrates in some
        regions and depletes in others. Phi spatial std must grow from
        its initial noise-only value.
        """
        params = _cell_params(phi0=0.05, t_end=600.0)
        props = _cell_props()
        gel = solve_nips_cellulose(
            params, props, R_droplet=200e-6, n_r=30, seed=0,
        )
        phi_std_final = gel.model_manifest.diagnostics["phi_std_final"]
        # Initial noise amplitude 1 % * phi_0 = 0.0005. Final must exceed this.
        assert phi_std_final > 0.0005


# ─── 3. Zero cellulose → zero gel + zero modulus ────────────────────


class TestZeroCellulose:
    def test_zero_phi_returns_no_gel_result(self):
        params = _cell_params(phi0=0.0, t_end=100.0)
        props = _cell_props()
        gel = solve_nips_cellulose(
            params, props, R_droplet=100e-6, n_r=20,
        )
        assert gel.alpha_final == 0.0
        assert gel.porosity == pytest.approx(1.0)
        assert gel.model_manifest.evidence_tier == ModelEvidenceTier.UNSUPPORTED

    def test_zero_phi_mean_gives_zero_modulus(self):
        params = _cell_params(phi0=0.0, t_end=100.0)
        props = _cell_props()
        gel = solve_nips_cellulose(
            params, props, R_droplet=100e-6, n_r=20,
        )
        mech = solve_mechanical_cellulose(params, props, gel)
        assert mech.G_DN == 0.0


# ─── 4. Modulus scaling with phi ────────────────────────────────────


class TestModulusScaling:
    def test_doubles_phi_scales_as_phi_alpha(self):
        """cellulose_modulus must respect the power-law exponent."""
        K = 5.0e5
        for alpha in (2.0, 2.25, 2.5):
            G1 = cellulose_modulus(0.05, K, alpha)
            G2 = cellulose_modulus(0.10, K, alpha)
            assert G2 / G1 == pytest.approx(2.0 ** alpha, rel=1e-6)

    def test_prefactor_scales_linearly(self):
        K = 1.0e5
        G = cellulose_modulus(0.1, K, 2.0)
        assert cellulose_modulus(0.1, 2.0 * K, 2.0) == pytest.approx(2.0 * G)


# ─── 5. Manifest schema (tier, model name) ──────────────────────────


class TestManifestSchema:
    def test_l2_nips_manifest_is_semi_quantitative(self):
        params = _cell_params(phi0=0.05, t_end=100.0)
        props = _cell_props()
        gel = solve_nips_cellulose(
            params, props, R_droplet=100e-6, n_r=20,
        )
        assert gel.model_manifest is not None
        assert gel.model_manifest.model_name == "L2.Gelation.NIPSCellulose"
        assert (
            gel.model_manifest.evidence_tier
            == ModelEvidenceTier.SEMI_QUANTITATIVE
        )
        assert "phi_mean_final" in gel.model_manifest.diagnostics
        assert "bicontinuous_score" in gel.model_manifest.diagnostics

    def test_l4_cellulose_manifest_schema(self):
        params = _cell_params(phi0=0.05, t_end=100.0)
        props = _cell_props()
        gel = solve_nips_cellulose(
            params, props, R_droplet=100e-6, n_r=20,
        )
        mech = solve_mechanical_cellulose(params, props, gel)
        assert mech.model_used == "cellulose_zhang2020"
        assert mech.network_type == "physical_entangled"
        assert mech.G_agarose == 0.0
        assert mech.G_chitosan == 0.0
        assert mech.model_manifest is not None
        assert (
            mech.model_manifest.model_name
            == "L4.Mechanical.CelluloseZhang2020"
        )
        assert (
            mech.model_manifest.evidence_tier
            == ModelEvidenceTier.SEMI_QUANTITATIVE
        )


# ─── 6. NaOH/urea preset ────────────────────────────────────────────


class TestCellulosePreset:
    def test_preset_registry_contains_naoh_urea(self):
        assert "naoh_urea" in CELLULOSE_SOLVENT_PRESETS
        assert CELLULOSE_SOLVENT_PRESETS["naoh_urea"] is NAOH_UREA

    def test_apply_preset_patches_all_nine_fields(self):
        props = MaterialProperties()
        apply_preset(props, "naoh_urea")
        assert props.N_p_cellulose == pytest.approx(NAOH_UREA.N_p)
        assert props.chi_PS_cellulose == pytest.approx(NAOH_UREA.chi_PS)
        assert props.chi_PN_cellulose == pytest.approx(NAOH_UREA.chi_PN)
        assert props.chi_SN_cellulose == pytest.approx(NAOH_UREA.chi_SN)
        assert props.D_solvent_cellulose == pytest.approx(NAOH_UREA.D_solvent)
        assert (
            props.D_nonsolvent_cellulose == pytest.approx(NAOH_UREA.D_nonsolvent)
        )
        assert props.kappa_CH_cellulose == pytest.approx(NAOH_UREA.kappa_CH)
        assert props.K_cell_modulus == pytest.approx(NAOH_UREA.K_cell)
        assert props.alpha_cell_modulus == pytest.approx(NAOH_UREA.alpha_cell)

    def test_apply_unknown_preset_raises(self):
        props = MaterialProperties()
        with pytest.raises(KeyError):
            apply_preset(props, "unicorn_saliva")


# ─── 7. L4 modulus zero-input sanity ────────────────────────────────


class TestModulusEdgeCases:
    def test_zero_phi_zero_modulus(self):
        assert cellulose_modulus(0.0, 5e5, 2.25) == 0.0

    def test_zero_prefactor_zero_modulus(self):
        assert cellulose_modulus(0.1, 0.0, 2.25) == 0.0

    def test_negative_phi_treated_as_zero(self):
        assert cellulose_modulus(-0.01, 5e5, 2.25) == 0.0


# ─── 8. Solver input validation ─────────────────────────────────────


class TestSolverValidation:
    def test_negative_R_raises(self):
        params = _cell_params()
        props = _cell_props()
        with pytest.raises(ValueError):
            solve_nips_cellulose(params, props, R_droplet=-1e-6)

    def test_tiny_grid_raises(self):
        params = _cell_params()
        props = _cell_props()
        with pytest.raises(ValueError):
            solve_nips_cellulose(params, props, R_droplet=100e-6, n_r=4)

    def test_phi0_above_one_raises(self):
        params = _cell_params(phi0=1.5)
        props = _cell_props()
        with pytest.raises(ValueError):
            solve_nips_cellulose(params, props, R_droplet=100e-6)

    def test_negative_time_raises(self):
        params = _cell_params()
        props = _cell_props()
        with pytest.raises(ValueError):
            solve_nips_cellulose(
                params, props, R_droplet=100e-6, time=-1.0,
            )

    def test_negative_noise_amp_raises(self):
        params = _cell_params()
        props = _cell_props()
        with pytest.raises(ValueError):
            solve_nips_cellulose(
                params, props, R_droplet=100e-6, noise_amp=-0.1,
            )
