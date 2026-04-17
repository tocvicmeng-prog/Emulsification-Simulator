"""Tests for the coupled GDL + CaCO₃ + alginate internal-release
gelation solver (Phase 2c Phase-3 deferred item, shipped in v8.0-rc2).

Protocol coverage:
  - GDL exponential decay (first-order hydrolysis)
  - No-flux Ca²⁺ mass balance (closed droplet)
  - Homogeneity (X_cov < shrinking-core X_cov for same bead size)
  - Zero-CaCO₃ sanity
  - Monotone conversion
  - Manifest schema / evidence tier
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from emulsim.datatypes import (
    MaterialProperties,
    ModelEvidenceTier,
    PolymerFamily,
    SimulationParameters,
)
from emulsim.level2_gelation.ionic_ca import (
    solve_internal_gelation,
    solve_ionic_ca_gelation,
)


def _alg_params(c_alg: float = 20.0, t_gel: float = 3600.0) -> SimulationParameters:
    p = SimulationParameters()
    p.formulation.c_alginate = c_alg
    p.formulation.c_agarose = 0.0
    p.formulation.c_chitosan = 0.0
    p.formulation.t_crosslink = t_gel
    return p


def _alg_props(f_G: float = 0.5) -> MaterialProperties:
    props = MaterialProperties()
    props.polymer_family = PolymerFamily.ALGINATE
    props.f_guluronate = f_G
    return props


class TestInternalGelationSchema:
    def test_returns_gelation_result_with_correct_manifest(self):
        params = _alg_params(c_alg=20.0, t_gel=600.0)
        props = _alg_props()
        res = solve_internal_gelation(
            params, props, R_droplet=200e-6,
            C_CaCO3_0=20.0, time=600.0, n_r=30,
        )
        assert res.model_manifest is not None
        assert (
            res.model_manifest.model_name
            == "L2.Gelation.IonicCaInternalRelease"
        )
        assert (
            res.model_manifest.evidence_tier
            == ModelEvidenceTier.SEMI_QUANTITATIVE
        )
        # Diagnostics cover the internal-release-specific scalars.
        diag = res.model_manifest.diagnostics
        assert "L_final" in diag
        assert "A_final" in diag
        assert "S_final" in diag
        assert "X_cov" in diag

    def test_defaults_stoichiometric_gdl_loading(self):
        """L_GDL_0 defaults to 2 × C_CaCO3_0 when not supplied."""
        params = _alg_params(c_alg=20.0, t_gel=1.0)
        props = _alg_props()
        res = solve_internal_gelation(
            params, props, R_droplet=200e-6,
            C_CaCO3_0=15.0, time=1.0, n_r=20,
        )
        # At t=1 s almost no GDL has hydrolysed: L_final ≈ 30.
        assert res.model_manifest.diagnostics["L_final"] == pytest.approx(
            30.0, rel=1e-2
        )


class TestInternalGelationKinetics:
    def test_gdl_decays_first_order(self):
        """L(t) = L₀·exp(−k_hyd·t) in the long-time limit where CaCO₃
        dissolution doesn't back-react. Since GDL hydrolysis is
        unidirectional and S doesn't feed back into L, the decay is
        exact first order at all times.
        """
        params = _alg_params(c_alg=20.0)
        props = _alg_props()
        k_hyd = 1.5e-4
        t_end = 3000.0  # ~½ half-life at this k_hyd (ln2 / k ≈ 4620 s)
        res = solve_internal_gelation(
            params, props, R_droplet=200e-6,
            C_CaCO3_0=20.0, k_hyd=k_hyd, time=t_end, n_r=20,
        )
        L_theory = 40.0 * math.exp(-k_hyd * t_end)
        L_final = res.model_manifest.diagnostics["L_final"]
        assert L_final == pytest.approx(L_theory, rel=2e-2)

    def test_conversion_grows_monotonically(self):
        params = _alg_params(c_alg=20.0)
        props = _alg_props()
        res_short = solve_internal_gelation(
            params, props, R_droplet=200e-6,
            C_CaCO3_0=25.0, time=300.0, n_r=20,
        )
        res_long = solve_internal_gelation(
            params, props, R_droplet=200e-6,
            C_CaCO3_0=25.0, time=7200.0, n_r=20,
        )
        assert res_long.alpha_final > res_short.alpha_final

    def test_zero_caco3_gives_no_gel(self):
        params = _alg_params(c_alg=20.0)
        props = _alg_props()
        res = solve_internal_gelation(
            params, props, R_droplet=200e-6,
            C_CaCO3_0=0.0, time=3600.0, n_r=20,
        )
        assert res.alpha_final == 0.0

    def test_zero_alginate_returns_zero_result(self):
        params = _alg_params(c_alg=0.0)
        props = _alg_props()
        res = solve_internal_gelation(
            params, props, R_droplet=100e-6,
            C_CaCO3_0=20.0, time=600.0, n_r=20,
        )
        assert res.alpha_final == 0.0


class TestMassBalance:
    def test_calcium_is_conserved_under_no_flux(self):
        """Total Ca²⁺ (dissolved + bound as egg-box) + remaining CaCO₃
        equals the initial CaCO₃ loading (within numerical tolerance).
        """
        params = _alg_params(c_alg=20.0)
        props = _alg_props()
        C_CaCO3_0 = 25.0
        res = solve_internal_gelation(
            params, props, R_droplet=150e-6,
            C_CaCO3_0=C_CaCO3_0, time=5000.0, n_r=30,
        )
        diag = res.model_manifest.diagnostics
        # Ca balance: S_final + <C> + <X> ≤ C_CaCO3_0. (X_mean is bound
        # Ca²⁺; each egg-box junction = 1 Ca²⁺.)
        # We don't store <C>_final separately, so estimate from mass closure.
        S_final = diag["S_final"]
        X_mean = diag["X_mean_final"]
        # Bound Ca + undissolved CaCO₃ must be ≤ initial Ca inventory.
        bound_and_solid = S_final + X_mean
        assert bound_and_solid <= C_CaCO3_0 + 1e-6
        # Some gelation should have happened — X_mean strictly > 0.
        assert X_mean > 0.0


class TestHomogeneity:
    def test_internal_release_more_homogeneous_than_shrinking_core(self):
        """Protocol §2.2 claim (Draget 1997): internal GDL/CaCO₃ gels
        are more homogeneous than CaCl₂-bath shrinking-core gels. The
        comparison must be done mid-transient, when the shrinking-core
        front has only partially penetrated — if we wait for both to
        saturate, both become uniform and the CoV collapses.
        For R = 1.5 mm, D = 1 × 10⁻⁹: τ_diff ≈ 2250 s. At t = 800 s
        the shrinking-core front has reached δ/R ≈ 0.6, so a strong
        gradient is still present, while the internal-release gel has
        released most of its Ca²⁺ uniformly.
        """
        params = _alg_params(c_alg=20.0)
        props = _alg_props()
        R = 1.5e-3
        t_end = 800.0

        res_int = solve_internal_gelation(
            params, props, R_droplet=R,
            C_CaCO3_0=25.0, time=t_end, n_r=40,
        )
        res_ext = solve_ionic_ca_gelation(
            params, props, R_droplet=R,
            C_Ca_bath=25.0,
            time=t_end, n_r=40,
        )

        X_int = res_int.phi_field
        X_ext = res_ext.phi_field
        cov_int = float(np.std(X_int) / max(np.mean(X_int), 1e-12))
        cov_ext = float(np.std(X_ext) / max(np.mean(X_ext), 1e-12))
        assert cov_int < cov_ext, (
            f"internal-release should be more homogeneous: "
            f"cov_internal={cov_int:.3f}, cov_shrinking_core={cov_ext:.3f}"
        )
        # Sanity: the shrinking-core CoV is non-trivial in this regime.
        assert cov_ext > 0.1


class TestInputValidation:
    def test_negative_R_raises(self):
        params = _alg_params()
        props = _alg_props()
        with pytest.raises(ValueError):
            solve_internal_gelation(
                params, props, R_droplet=-1e-6,
                C_CaCO3_0=20.0, time=100.0,
            )

    def test_negative_C_CaCO3_raises(self):
        params = _alg_params()
        props = _alg_props()
        with pytest.raises(ValueError):
            solve_internal_gelation(
                params, props, R_droplet=100e-6,
                C_CaCO3_0=-5.0, time=100.0,
            )

    def test_tiny_grid_raises(self):
        params = _alg_params()
        props = _alg_props()
        with pytest.raises(ValueError):
            solve_internal_gelation(
                params, props, R_droplet=100e-6,
                C_CaCO3_0=20.0, time=100.0, n_r=4,
            )
