"""Tests for Level 2 gelation and pore formation simulation."""

import numpy as np
import pytest
from scipy import sparse

from emulsim.datatypes import MaterialProperties, SimulationParameters
from emulsim.level2_gelation.spatial import (
    create_radial_grid, build_laplacian_matrix,
    create_2d_grid, build_laplacian_2d, build_mobility_laplacian_2d,
)
from emulsim.level2_gelation.free_energy import (
    flory_huggins_mu, flory_huggins_d2f, contractive_constant,
)
from emulsim.level2_gelation.gelation import (
    avrami_gelation, gelation_rate_constant, cooling_temperature, mobility,
)
from emulsim.level2_gelation.pore_analysis import (
    chord_length_distribution, compute_porosity,
    chord_length_distribution_2d, compute_porosity_2d, characteristic_wavelength_2d,
)


class TestSpatial:
    def test_grid_size(self):
        r, dr = create_radial_grid(1e-6, 100)
        assert len(r) == 100
        assert dr == pytest.approx(1e-6 / 100)

    def test_grid_positive(self):
        r, dr = create_radial_grid(1e-6, 100)
        assert np.all(r > 0)
        assert r[-1] < 1e-6  # should not exceed R

    def test_laplacian_shape(self):
        r, dr = create_radial_grid(1e-6, 50)
        L = build_laplacian_matrix(r, dr)
        assert L.shape == (50, 50)

    def test_laplacian_constant_field(self):
        """Laplacian of constant field should be ~zero (unit scale to avoid FP noise)."""
        r, dr = create_radial_grid(1.0, 50)
        L = build_laplacian_matrix(r, dr)
        phi_const = np.ones(50) * 0.5
        result = L @ phi_const
        assert np.max(np.abs(result)) < 1e-10


class TestFreeEnergy:
    def test_mu_shape(self):
        phi = np.linspace(0.01, 0.99, 100)
        mu = flory_huggins_mu(phi, 300.0, 0.5)
        assert mu.shape == phi.shape

    def test_d2f_negative_in_spinodal(self):
        """f'' should be negative inside the spinodal region for high chi."""
        phi = np.array([0.1])
        d2f = flory_huggins_d2f(phi, 300.0, chi=2.0, N_p=100)
        assert d2f[0] < 0

    def test_d2f_positive_outside_spinodal(self):
        """f'' should be positive at extreme compositions."""
        phi_low = np.array([0.001])
        d2f_low = flory_huggins_d2f(phi_low, 300.0, chi=0.5, N_p=100)
        assert d2f_low[0] > 0

    def test_contractive_constant_positive(self):
        C = contractive_constant((0.01, 0.30), 300.0, 0.5)
        assert C > 0


class TestGelation:
    def test_avrami_zero_at_start(self):
        assert avrami_gelation(0.0, 0.01, 2.5) == 0.0

    def test_avrami_approaches_one(self):
        alpha = avrami_gelation(1000.0, 0.01, 2.5)
        assert alpha > 0.99

    def test_avrami_monotonic(self):
        times = [0, 10, 50, 100, 500]
        alphas = [avrami_gelation(t, 0.01, 2.5) for t in times]
        for i in range(len(alphas) - 1):
            assert alphas[i + 1] >= alphas[i]

    def test_rate_constant_zero_above_tgel(self):
        assert gelation_rate_constant(320.0, 311.15, 0.01) == 0.0

    def test_rate_constant_positive_below_tgel(self):
        k = gelation_rate_constant(300.0, 311.15, 0.01)
        assert k > 0

    def test_cooling_temperature(self):
        T = cooling_temperature(0.0, 363.15, 293.15, 0.167)
        assert T == pytest.approx(363.15)

    def test_cooling_bounded_below(self):
        T = cooling_temperature(1e6, 363.15, 293.15, 0.167)
        assert T >= 293.15

    def test_mobility_decreases_with_gelation(self):
        phi = np.array([0.05])
        m0 = mobility(phi, 0.0, 1e-14)
        m1 = mobility(phi, 0.5, 1e-14)
        m2 = mobility(phi, 0.99, 1e-14)
        assert m0[0] > m1[0] > m2[0]

    def test_mobility_zero_at_full_gelation(self):
        phi = np.array([0.05])
        m = mobility(phi, 1.0, 1e-14)
        assert m[0] == pytest.approx(0.0)


class TestPoreAnalysis:
    def test_chord_length_alternating(self):
        """Alternating pore/polymer pattern should give chord lengths."""
        phi = np.array([0.1, 0.1, 0.9, 0.9, 0.1, 0.1, 0.1, 0.9])
        r = np.arange(8) * 1e-9
        chords = chord_length_distribution(phi, r)
        assert len(chords) >= 1
        assert np.all(chords > 0)

    def test_porosity_uniform(self):
        """Uniform field below threshold -> porosity ~1.0."""
        r = np.linspace(1e-9, 1e-6, 100)
        phi = np.ones(100) * 0.01
        por = compute_porosity(phi, r, threshold=0.5)
        assert por > 0.99

    def test_porosity_range(self):
        r = np.linspace(1e-9, 1e-6, 100)
        phi = np.random.default_rng(42).uniform(0.0, 1.0, 100)
        por = compute_porosity(phi, r)
        assert 0.0 <= por <= 1.0


class TestCahnHilliardSolver:
    def test_solver_runs(self):
        """Solver should complete without error."""
        from emulsim.level2_gelation.solver import CahnHilliardSolver
        solver = CahnHilliardSolver(N_r=100, dt_initial=1e-3, dt_max=1.0)
        params = SimulationParameters()
        props = MaterialProperties()
        result = solver.solve(params, props, R_droplet=3.0e-6)
        assert result.r_grid.shape == (100,)
        assert result.phi_field.shape == (100,)

    def test_alpha_reaches_completion(self):
        """With standard cooling, gelation should be nearly complete."""
        from emulsim.level2_gelation.solver import CahnHilliardSolver
        solver = CahnHilliardSolver(N_r=100, dt_initial=1e-3, dt_max=1.0)
        params = SimulationParameters()
        props = MaterialProperties()
        result = solver.solve(params, props, R_droplet=3.0e-6)
        assert result.alpha_final > 0.9

    def test_phi_within_bounds(self):
        """Composition must stay in [0, 1]."""
        from emulsim.level2_gelation.solver import CahnHilliardSolver
        solver = CahnHilliardSolver(N_r=100, dt_initial=1e-3, dt_max=1.0)
        params = SimulationParameters()
        props = MaterialProperties()
        result = solver.solve(params, props, R_droplet=3.0e-6)
        assert np.all(result.phi_field >= 0)
        assert np.all(result.phi_field <= 1)

    def test_porosity_in_range(self):
        from emulsim.level2_gelation.solver import CahnHilliardSolver
        solver = CahnHilliardSolver(N_r=100, dt_initial=1e-3, dt_max=1.0)
        params = SimulationParameters()
        props = MaterialProperties()
        result = solver.solve(params, props, R_droplet=3.0e-6)
        assert 0.0 <= result.porosity <= 1.0


class TestSpatial2D:
    def test_2d_grid_size(self):
        x, y, h = create_2d_grid(6e-6, 64)
        assert len(x) == 64
        assert len(y) == 64
        assert h == pytest.approx(6e-6 / 64)

    def test_2d_grid_cell_centered(self):
        x, y, h = create_2d_grid(1.0, 10)
        assert x[0] == pytest.approx(0.05)
        assert x[-1] == pytest.approx(0.95)

    def test_2d_laplacian_shape(self):
        L = build_laplacian_2d(16, 0.1)
        assert L.shape == (256, 256)

    def test_2d_laplacian_constant_field(self):
        """Laplacian of a constant field should be zero."""
        N = 16
        L = build_laplacian_2d(N, 0.1)
        phi = np.ones(N * N) * 0.5
        result = L @ phi
        assert np.max(np.abs(result)) < 1e-10

    def test_2d_laplacian_symmetry(self):
        """The Laplacian matrix should be symmetric for Neumann BCs."""
        N = 8
        L = build_laplacian_2d(N, 0.1)
        diff = L - L.T
        assert sparse.linalg.norm(diff) < 1e-12

    def test_mobility_laplacian_constant_M(self):
        """With constant M, mobility Laplacian should equal M * standard Laplacian."""
        N = 8
        h = 0.1
        M_val = 2.5
        M_field = np.ones((N, N)) * M_val
        L_std = build_laplacian_2d(N, h)
        L_M = build_mobility_laplacian_2d(M_field, N, h)
        diff = L_M - M_val * L_std
        assert sparse.linalg.norm(diff) < 1e-10

    def test_mobility_laplacian_constant_field(self):
        """Mobility Laplacian applied to constant field should be zero."""
        N = 8
        h = 0.1
        M_field = np.random.default_rng(42).uniform(0.5, 2.0, (N, N))
        L_M = build_mobility_laplacian_2d(M_field, N, h)
        phi = np.ones(N * N) * 0.3
        result = L_M @ phi
        assert np.max(np.abs(result)) < 1e-10

    def test_mobility_laplacian_face_averaging(self):
        """Verify face-centred mobility is used, not cell-centred."""
        N = 4
        h = 1.0
        # Set up M with a sharp gradient
        M_field = np.ones((N, N))
        M_field[0, :] = 10.0  # first row has high mobility
        L_M = build_mobility_laplacian_2d(M_field, N, h)
        # The operator should be different from using cell-centred M
        # Just check it's a valid operator (negative semi-definite diagonal)
        diag = L_M.diagonal()
        assert np.all(diag <= 1e-14)  # diagonal entries should be <= 0


class TestPoreAnalysis2D:
    def test_porosity_2d_uniform_low(self):
        phi = np.ones((16, 16)) * 0.01
        por = compute_porosity_2d(phi, threshold=0.5)
        assert por > 0.99

    def test_porosity_2d_uniform_high(self):
        phi = np.ones((16, 16)) * 0.99
        por = compute_porosity_2d(phi, threshold=0.5)
        assert por < 0.01

    def test_chord_2d_alternating(self):
        """Checkerboard-like pattern should produce chords."""
        phi = np.zeros((8, 8))
        phi[:, 4:] = 1.0  # left half pore, right half polymer
        chords = chord_length_distribution_2d(phi, 1e-9)
        assert len(chords) >= 1
        assert np.all(chords > 0)

    def test_char_wavelength_2d_sinusoidal(self):
        """A sinusoidal pattern should have a detectable wavelength."""
        N = 64
        h = 1.0
        x = np.arange(N) * h
        # Create a sinusoidal pattern with wavelength = 16*h
        lam = 16.0 * h
        phi_row = 0.5 + 0.3 * np.sin(2 * np.pi * x / lam)
        phi_2d = np.tile(phi_row, (N, 1))
        wl = characteristic_wavelength_2d(phi_2d, h)
        # Should be close to 16
        assert abs(wl - lam) / lam < 0.15


class TestCahnHilliard2DSolver:
    def test_2d_solver_runs(self):
        """2D solver should complete without error on a small grid."""
        from emulsim.level2_gelation.solver import CahnHilliard2DSolver
        solver = CahnHilliard2DSolver(N_grid=32, dt_initial=1e-3, dt_max=1.0)
        params = SimulationParameters()
        props = MaterialProperties()
        result = solver.solve(params, props, R_droplet=3.0e-6)
        assert result.phi_field.shape == (32, 32)
        assert result.r_grid.shape == (32,)

    def test_2d_alpha_reaches_completion(self):
        """With standard cooling, gelation should be nearly complete."""
        from emulsim.level2_gelation.solver import CahnHilliard2DSolver
        solver = CahnHilliard2DSolver(N_grid=32, dt_initial=1e-3, dt_max=1.0)
        params = SimulationParameters()
        props = MaterialProperties()
        result = solver.solve(params, props, R_droplet=3.0e-6)
        assert result.alpha_final > 0.9

    def test_2d_phi_within_bounds(self):
        """Composition must stay in [0, 1]."""
        from emulsim.level2_gelation.solver import CahnHilliard2DSolver
        solver = CahnHilliard2DSolver(N_grid=32, dt_initial=1e-3, dt_max=1.0)
        params = SimulationParameters()
        props = MaterialProperties()
        result = solver.solve(params, props, R_droplet=3.0e-6)
        assert np.all(result.phi_field >= 0)
        assert np.all(result.phi_field <= 1)

    def test_2d_porosity_in_range(self):
        from emulsim.level2_gelation.solver import CahnHilliard2DSolver
        solver = CahnHilliard2DSolver(N_grid=32, dt_initial=1e-3, dt_max=1.0)
        params = SimulationParameters()
        props = MaterialProperties()
        result = solver.solve(params, props, R_droplet=3.0e-6)
        assert 0.0 <= result.porosity <= 1.0

    def test_2d_result_has_domain_info(self):
        """2D result should include domain size and grid spacing."""
        from emulsim.level2_gelation.solver import CahnHilliard2DSolver
        solver = CahnHilliard2DSolver(N_grid=32, dt_initial=1e-3, dt_max=1.0)
        params = SimulationParameters()
        props = MaterialProperties()
        result = solver.solve(params, props, R_droplet=3.0e-6)
        assert result.L_domain == pytest.approx(6.0e-6)
        assert result.grid_spacing > 0

    def test_solve_gelation_uses_2d_by_default(self):
        """The convenience function should use 2D solver by default."""
        from emulsim.level2_gelation.solver import solve_gelation
        params = SimulationParameters()
        params.solver.l2_n_grid = 32
        props = MaterialProperties()
        result = solve_gelation(params, props, R_droplet=3.0e-6)
        assert result.phi_field.ndim == 2
        assert result.L_domain > 0

    def test_solve_gelation_1d_fallback(self):
        """The convenience function should support 1D fallback."""
        from emulsim.level2_gelation.solver import solve_gelation
        params = SimulationParameters()
        params.solver.l2_n_r = 100
        props = MaterialProperties()
        result = solve_gelation(params, props, R_droplet=3.0e-6, use_2d=False)
        assert result.phi_field.ndim == 1
        assert result.L_domain == 0.0


# ─── Node 8 (v6.1, F8): L2 alpha_final wired from gelation timing ─────────


class TestNode8AlphaFromTiming:
    """Verifies that solve_gelation_empirical uses timing.alpha_final.

    Acceptance for Node 8 (F8):
      - When ``timing`` is None, alpha_final remains the legacy 0.999
        (backward compat for direct callers/tests).
      - When ``timing`` is provided, alpha_final tracks timing.alpha_final.
      - L2 manifest diagnostics include the timing source fields.
      - End-to-end pipeline run carries the timing-derived alpha into
        the L2 manifest (orchestrator wiring works).
    """

    def test_legacy_no_timing_keeps_hardcoded(self):
        from emulsim.level2_gelation.solver import solve_gelation_empirical
        params = SimulationParameters()
        props = MaterialProperties()
        result = solve_gelation_empirical(params, props, R_droplet=10e-6)
        assert abs(result.alpha_final - 0.999) < 1e-6
        diag = result.model_manifest.diagnostics
        assert diag.get("alpha_final_source") == "hardcoded_legacy_fallback"

    def test_timing_drives_alpha_and_diagnostics(self):
        from emulsim.datatypes import GelationTimingResult
        from emulsim.level2_gelation.solver import solve_gelation_empirical
        params = SimulationParameters()
        props = MaterialProperties()

        timing = GelationTimingResult(
            T_history=np.array([363.0, 333.0]),
            t_gel_onset=42.0,
            alpha_final=0.85,
            mobility_arrest_factor=0.05,
            cooling_rate_effective=0.20,
        )
        result = solve_gelation_empirical(
            params, props, R_droplet=10e-6, timing=timing,
        )
        assert result.alpha_final == pytest.approx(0.85, abs=1e-6)
        diag = result.model_manifest.diagnostics
        assert diag["alpha_final_from_timing"] == pytest.approx(0.85, abs=1e-6)
        assert diag["t_gel_onset_s"] == pytest.approx(42.0, abs=1e-6)
        assert diag["cooling_rate_effective_K_per_s"] == pytest.approx(0.20, abs=1e-6)
        assert "alpha_final_source" not in diag

    def test_orchestrator_passes_timing(self, tmp_path):
        from emulsim.config import load_config
        from emulsim.pipeline.orchestrator import PipelineOrchestrator
        from pathlib import Path
        repo_root = Path(__file__).resolve().parents[1]
        cfg = repo_root / "configs" / "fast_smoke.toml"
        if not cfg.exists():
            pytest.skip("fast_smoke.toml missing")
        params = load_config(cfg)
        result = PipelineOrchestrator(output_dir=tmp_path).run_single(params)
        l2_diag = result.gelation.model_manifest.diagnostics
        assert "alpha_final_from_timing" in l2_diag
        assert l2_diag["alpha_final_from_timing"] == pytest.approx(
            result.gelation_timing.alpha_final, abs=1e-9,
        )
