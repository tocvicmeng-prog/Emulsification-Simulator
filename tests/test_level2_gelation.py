"""Tests for Level 2 gelation and pore formation simulation."""

import numpy as np
import pytest

from emulsim.datatypes import MaterialProperties, SimulationParameters
from emulsim.level2_gelation.spatial import create_radial_grid, build_laplacian_matrix
from emulsim.level2_gelation.free_energy import (
    flory_huggins_mu, flory_huggins_d2f, contractive_constant,
)
from emulsim.level2_gelation.gelation import (
    avrami_gelation, gelation_rate_constant, cooling_temperature, mobility,
)
from emulsim.level2_gelation.pore_analysis import (
    chord_length_distribution, compute_porosity,
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
