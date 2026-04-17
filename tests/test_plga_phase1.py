"""Tests for Node F1-c Phase 1 (v8.2-alpha): L2 solvent-evaporation +
L4 PLGA modulus + 4-grade registry.

Core 8 tests (per docs/f1c_plga_protocol.md §6):

  1. Fickian mass conservation / monotone DCM depletion
  2. Dirichlet sink drives evaporation (phi_PLGA → 1 at t → ∞)
  3. √t early-time front propagation
  4. Higher D_DCM → faster vitrification
  5. Modulus scales as phi² (Gibson-Ashby)
  6. Zero PLGA → zero modulus + UNSUPPORTED manifest
  7. L2 + L4 manifests are SEMI_QUANTITATIVE at solver exit
  8. PLGA 50:50 preset applies all 4 fields correctly

Plus edge-case / input-validation extras.
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
from emulsim.level2_gelation.solvent_evaporation import (
    solve_solvent_evaporation,
)
from emulsim.level4_mechanical.plga import (
    plga_modulus,
    solve_mechanical_plga,
)
from emulsim.properties.plga_defaults import (
    PLGA_50_50,
    PLGA_GRADE_PRESETS,
    apply_preset,
)


def _plga_params(
    phi0: float = 0.10, t_end: float = 1800.0,
) -> SimulationParameters:
    p = SimulationParameters()
    p.formulation.phi_PLGA_0 = phi0
    p.formulation.c_agarose = 0.0
    p.formulation.c_chitosan = 0.0
    p.formulation.c_alginate = 0.0
    p.formulation.phi_cellulose_0 = 0.0
    p.formulation.t_crosslink = t_end
    return p


def _plga_props(grade: str = "50_50") -> MaterialProperties:
    props = MaterialProperties()
    props.polymer_family = PolymerFamily.PLGA
    apply_preset(props, grade)
    return props


# ─── Test 1: Monotone DCM depletion ─────────────────────────────────


class TestMonotoneDepletion:
    def test_phi_plga_grows_with_time(self):
        params = _plga_params(phi0=0.10, t_end=1.0)
        props = _plga_props()
        g_short = solve_solvent_evaporation(
            params, props, R_droplet=50e-6, n_r=20, time=1.0,
        )
        g_long = solve_solvent_evaporation(
            params, props, R_droplet=50e-6, n_r=20, time=100.0,
        )
        mean_short = g_short.model_manifest.diagnostics["phi_plga_mean_final"]
        mean_long = g_long.model_manifest.diagnostics["phi_plga_mean_final"]
        assert mean_long > mean_short
        # Initial mean should be phi_0 = 0.10; final long-time should be
        # approaching 1.0.
        assert mean_short >= 0.10 - 1e-6
        assert mean_long > 0.5

    def test_dcm_monotone_decrease(self):
        """DCM volume fraction should only ever go down, not up.
        Pick probe times well inside the transient regime (τ ≈ R²/D;
        for R = 100 µm, D = 1e-9 that's ~10 s, so t = [2, 5, 9] s
        gives a clearly-decreasing sequence).
        """
        params = _plga_params(phi0=0.10)
        props = _plga_props()
        R = 100e-6
        g1 = solve_solvent_evaporation(params, props, R_droplet=R, n_r=25, time=2.0)
        g2 = solve_solvent_evaporation(params, props, R_droplet=R, n_r=25, time=5.0)
        g3 = solve_solvent_evaporation(params, props, R_droplet=R, n_r=25, time=9.0)
        d1 = g1.model_manifest.diagnostics["phi_dcm_mean_final"]
        d2 = g2.model_manifest.diagnostics["phi_dcm_mean_final"]
        d3 = g3.model_manifest.diagnostics["phi_dcm_mean_final"]
        assert d1 > d2 > d3


# ─── Test 2: Dirichlet sink drives full evaporation ─────────────────


class TestEvaporationCompletes:
    def test_long_time_approaches_dense_microsphere(self):
        """With phi_eq ≈ 0 the long-time limit is phi_PLGA → 1."""
        params = _plga_params(phi0=0.10, t_end=3600.0)
        props = _plga_props()
        gel = solve_solvent_evaporation(
            params, props, R_droplet=30e-6, n_r=20,
        )
        phi_mean = gel.model_manifest.diagnostics["phi_plga_mean_final"]
        # For R=30 µm, D=1e-9, τ_diff ≈ R²/D ≈ 0.9 s. 3600 s = 4000 τ →
        # residual DCM essentially equals the Dirichlet equilibrium.
        assert phi_mean == pytest.approx(1.0 - 0.005, abs=1e-2)


# ─── Test 3: √t front propagation at early times ────────────────────


class TestSqrtTFront:
    def test_mean_dcm_depletion_scales_as_sqrt_t_early(self):
        """Early-regime (δ/R << 1): amount of DCM lost from the
        interior scales as √(D·t). Equivalently,
        phi_PLGA_mean - phi_0 ∝ √t. Pick a bead large enough that
        R²/D >> t for all probe times.
        """
        params = _plga_params(phi0=0.10)
        props = _plga_props()
        R = 500e-6  # τ_diff = R²/D = 250 s
        t_list = [2.0, 8.0, 32.0]  # well below τ_diff
        excess = []
        for t_k in t_list:
            gel = solve_solvent_evaporation(
                params, props, R_droplet=R, n_r=40, time=t_k,
            )
            mean_k = gel.model_manifest.diagnostics["phi_plga_mean_final"]
            excess.append(max(mean_k - 0.10, 1e-12))
        # log-log slope ≈ 0.5
        logs_t = np.log(t_list)
        logs_e = np.log(np.array(excess))
        slope, _ = np.polyfit(logs_t, logs_e, 1)
        assert 0.35 <= slope <= 0.65, (
            f"early-time √t scaling broken: slope={slope:.3f}, "
            f"excess={excess}"
        )


# ─── Test 4: D_DCM controls evaporation speed ───────────────────────


class TestDiffusivityControlsSpeed:
    def test_higher_D_gives_earlier_vitrification(self):
        """Pick R = 200 µm so τ_diff = R²/D ~ 40 s, and use t_end ~
        5·τ = 200 s so the vitrification-probe grid (t_end/200) is
        fine enough (~1 s) to resolve the 2× D speed-up.
        """
        params_fast = _plga_params(phi0=0.10, t_end=200.0)
        params_slow = _plga_params(phi0=0.10, t_end=200.0)
        props_fast = _plga_props()
        props_slow = _plga_props()
        # 4× diffusivity → expect ~2× speed-up in vitrification time.
        props_fast.D_DCM_plga = 4.0 * props_slow.D_DCM_plga

        g_fast = solve_solvent_evaporation(
            params_fast, props_fast, R_droplet=200e-6, n_r=30,
        )
        g_slow = solve_solvent_evaporation(
            params_slow, props_slow, R_droplet=200e-6, n_r=30,
        )
        t_fast = g_fast.model_manifest.diagnostics["t_vitrification"]
        t_slow = g_slow.model_manifest.diagnostics["t_vitrification"]
        assert t_fast < t_slow, (
            f"t_vitrification(4D)={t_fast:.2f} should be < "
            f"t_vitrification(D)={t_slow:.2f}"
        )


# ─── Test 5: Gibson-Ashby modulus scaling ───────────────────────────


class TestGibsonAshbyScaling:
    def test_doubles_phi_gives_phi_alpha_scaling(self):
        G_glassy = 7.0e8
        for n in (1.5, 2.0, 2.5):
            G_low = plga_modulus(0.3, G_glassy, n)
            G_high = plga_modulus(0.6, G_glassy, n)
            assert G_high / G_low == pytest.approx(2.0 ** n, rel=1e-6)

    def test_dense_limit_recovers_glassy_modulus(self):
        G_glassy = 7.0e8
        assert plga_modulus(1.0, G_glassy, 2.0) == pytest.approx(G_glassy)


# ─── Test 6: Zero PLGA input ────────────────────────────────────────


class TestZeroPLGA:
    def test_zero_phi0_returns_unsupported_gel(self):
        params = _plga_params(phi0=0.0, t_end=100.0)
        props = _plga_props()
        gel = solve_solvent_evaporation(
            params, props, R_droplet=50e-6, n_r=20,
        )
        assert gel.alpha_final == 0.0
        assert gel.porosity == pytest.approx(1.0)
        assert gel.model_manifest.evidence_tier == ModelEvidenceTier.UNSUPPORTED

    def test_zero_phi0_yields_zero_modulus(self):
        params = _plga_params(phi0=0.0, t_end=100.0)
        props = _plga_props()
        gel = solve_solvent_evaporation(
            params, props, R_droplet=50e-6, n_r=20,
        )
        mech = solve_mechanical_plga(params, props, gel)
        assert mech.G_DN == 0.0


# ─── Test 7: Manifest tiers + schema ────────────────────────────────


class TestManifestSchema:
    def test_l2_manifest_is_semi_quantitative(self):
        params = _plga_params(phi0=0.10, t_end=60.0)
        props = _plga_props()
        gel = solve_solvent_evaporation(
            params, props, R_droplet=50e-6, n_r=20,
        )
        assert gel.model_manifest is not None
        assert (
            gel.model_manifest.model_name
            == "L2.Gelation.SolventEvaporationPLGA"
        )
        assert (
            gel.model_manifest.evidence_tier
            == ModelEvidenceTier.SEMI_QUANTITATIVE
        )
        for k in (
            "phi_plga_mean_final", "phi_dcm_mean_final",
            "t_vitrification", "skin_thickness_proxy",
            "core_porosity_proxy", "R_shrunk_m",
        ):
            assert k in gel.model_manifest.diagnostics

    def test_l4_manifest_schema(self):
        params = _plga_params(phi0=0.10, t_end=60.0)
        props = _plga_props()
        gel = solve_solvent_evaporation(
            params, props, R_droplet=50e-6, n_r=20,
        )
        mech = solve_mechanical_plga(params, props, gel)
        assert mech.model_used == "plga_gibson_ashby"
        assert mech.network_type == "glassy_polymer"
        assert mech.G_agarose == 0.0
        assert mech.G_chitosan == 0.0
        assert mech.model_manifest is not None
        assert (
            mech.model_manifest.model_name
            == "L4.Mechanical.PLGAGibsonAshby"
        )
        assert (
            mech.model_manifest.evidence_tier
            == ModelEvidenceTier.SEMI_QUANTITATIVE
        )


# ─── Test 8: Grade preset application ───────────────────────────────


class TestGradePreset:
    def test_all_four_grades_registered(self):
        assert set(PLGA_GRADE_PRESETS) == {
            "50_50", "75_25", "85_15", "pla",
        }

    def test_apply_preset_patches_plga_fields(self):
        props = MaterialProperties()
        apply_preset(props, "50_50")
        assert props.D_DCM_plga == pytest.approx(PLGA_50_50.D_DCM)
        assert props.phi_DCM_eq == pytest.approx(PLGA_50_50.phi_DCM_eq)
        assert props.G_glassy_plga == pytest.approx(PLGA_50_50.G_glassy)
        assert props.n_plga_modulus == pytest.approx(PLGA_50_50.n_plga)

    def test_apply_unknown_grade_raises(self):
        props = MaterialProperties()
        with pytest.raises(KeyError):
            apply_preset(props, "100_0_pgla_mirror_world")

    def test_each_grade_has_physically_reasonable_params(self):
        for name, g in PLGA_GRADE_PRESETS.items():
            assert 0.0 <= g.L_fraction <= 1.0
            assert 5_000 < g.M_n < 500_000, f"{name}: M_n range"
            assert 20 <= g.T_g_C <= 70, f"{name}: T_g range"
            assert 1e-11 < g.D_DCM < 1e-8, f"{name}: D_DCM range"
            assert 1e7 < g.G_glassy < 1e10, f"{name}: G_glassy range"
            assert 1.5 < g.n_plga < 3.5, f"{name}: exponent range"
            assert 0 < g.phi_PLGA_0_typical < 0.5, (
                f"{name}: typical loading range"
            )

    def test_grade_switching_gives_different_moduli(self):
        """Different PLGA grades produce different G_glassy and
        therefore different G_DN for the same phi_mean."""
        params = _plga_params(phi0=0.10, t_end=60.0)
        Gs = {}
        for name in PLGA_GRADE_PRESETS:
            props = _plga_props(grade=name)
            gel = solve_solvent_evaporation(
                params, props, R_droplet=50e-6, n_r=20,
            )
            mech = solve_mechanical_plga(params, props, gel)
            Gs[name] = mech.G_DN
        vals = list(Gs.values())
        assert max(vals) >= 1.4 * min(vals), (
            f"grades must produce distinguishable moduli: {Gs}"
        )


# ─── Extras: edge cases + validation ────────────────────────────────


class TestModulusEdgeCases:
    def test_zero_phi_zero_modulus(self):
        assert plga_modulus(0.0, 7e8, 2.0) == 0.0

    def test_zero_prefactor_zero_modulus(self):
        assert plga_modulus(0.5, 0.0, 2.0) == 0.0

    def test_negative_phi_treated_as_zero(self):
        assert plga_modulus(-0.1, 7e8, 2.0) == 0.0


class TestSolverValidation:
    def test_negative_R_raises(self):
        params = _plga_params()
        props = _plga_props()
        with pytest.raises(ValueError):
            solve_solvent_evaporation(params, props, R_droplet=-1e-6)

    def test_tiny_grid_raises(self):
        params = _plga_params()
        props = _plga_props()
        with pytest.raises(ValueError):
            solve_solvent_evaporation(
                params, props, R_droplet=50e-6, n_r=4,
            )

    def test_phi0_above_one_raises(self):
        params = _plga_params(phi0=1.5)
        props = _plga_props()
        with pytest.raises(ValueError):
            solve_solvent_evaporation(params, props, R_droplet=50e-6)

    def test_negative_time_raises(self):
        params = _plga_params()
        props = _plga_props()
        with pytest.raises(ValueError):
            solve_solvent_evaporation(
                params, props, R_droplet=50e-6, time=-1.0,
            )

    def test_non_positive_D_raises(self):
        params = _plga_params()
        props = _plga_props()
        props.D_DCM_plga = 0.0
        with pytest.raises(ValueError):
            solve_solvent_evaporation(params, props, R_droplet=50e-6)


class TestMassBalancePostProcessing:
    def test_R_shrunk_follows_cube_root_of_phi_ratio(self):
        """For a dense microsphere (phi_mean → 1), R_shrunk should
        approach R_0 · phi_0^(1/3). Check at long time.
        """
        params = _plga_params(phi0=0.10, t_end=3600.0)
        props = _plga_props()
        R_0 = 30e-6
        gel = solve_solvent_evaporation(
            params, props, R_droplet=R_0, n_r=20,
        )
        R_shrunk = gel.model_manifest.diagnostics["R_shrunk_m"]
        R_expect = R_0 * (0.10 / 0.995) ** (1.0 / 3.0)
        assert R_shrunk == pytest.approx(R_expect, rel=5e-2)
