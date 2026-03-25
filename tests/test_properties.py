"""Tests for material property models."""

import numpy as np
import pytest

from emulsim.properties.database import PropertyDatabase
from emulsim.properties.interfacial import interfacial_tension_span80, sigma_bare
from emulsim.properties.thermodynamic import (
    chi_flory_huggins,
    flory_huggins_free_energy,
    flory_huggins_second_derivative,
    mesh_size_canal_peppas,
    spinodal_boundary,
)
from emulsim.properties.viscosity import (
    intrinsic_viscosity_agarose,
    solution_viscosity,
    water_viscosity,
)


class TestViscosity:
    def test_water_viscosity_at_25C(self):
        mu = water_viscosity(298.15)
        assert 0.5e-3 < mu < 1.5e-3  # ~0.89 mPa·s

    def test_water_viscosity_decreases_with_T(self):
        assert water_viscosity(353.15) < water_viscosity(298.15)

    def test_intrinsic_viscosity_positive(self):
        eta = intrinsic_viscosity_agarose(120000.0)
        assert eta > 0

    def test_solution_viscosity_above_water(self):
        mu_sol = solution_viscosity(363.15, 42.0, 18.0)
        mu_water = water_viscosity(363.15)
        assert mu_sol > mu_water

    def test_solution_viscosity_increases_with_concentration(self):
        mu_low = solution_viscosity(363.15, 20.0, 10.0)
        mu_high = solution_viscosity(363.15, 42.0, 18.0)
        assert mu_high > mu_low

    def test_pure_water_when_no_polymer(self):
        mu = solution_viscosity(363.15, 0.0, 0.0)
        mu_water = water_viscosity(363.15)
        assert mu == pytest.approx(mu_water, rel=0.01)


class TestInterfacialTension:
    def test_bare_positive(self):
        assert sigma_bare(293.15) > 0
        assert sigma_bare(363.15) > 0

    def test_bare_decreases_with_T(self):
        assert sigma_bare(363.15) < sigma_bare(293.15)

    def test_span80_reduces_tension(self):
        sigma_no = sigma_bare(363.15)
        sigma_with = interfacial_tension_span80(363.15, 20.0)
        assert sigma_with < sigma_no

    def test_span80_monotonic_decrease(self):
        """More Span-80 → lower interfacial tension."""
        sigmas = [interfacial_tension_span80(363.15, c) for c in [1, 5, 10, 20, 30]]
        for i in range(len(sigmas) - 1):
            assert sigmas[i + 1] <= sigmas[i]

    def test_span80_positive(self):
        """Interfacial tension must be positive."""
        sigma = interfacial_tension_span80(363.15, 100.0)
        assert sigma > 0

    def test_span80_at_2pct(self):
        """At 2% w/v (~20 kg/m³), σ should be ~2-8 mN/m."""
        sigma = interfacial_tension_span80(363.15, 20.0)
        assert 0.001 < sigma < 0.010  # 1-10 mN/m


class TestThermodynamic:
    def test_chi_positive(self):
        chi = chi_flory_huggins(298.15)
        assert chi > 0

    def test_chi_increases_as_T_decreases(self):
        chi_hot = chi_flory_huggins(363.15)
        chi_cold = chi_flory_huggins(298.15)
        assert chi_cold > chi_hot

    def test_flory_huggins_convex_outside_spinodal(self):
        """f''(φ) > 0 at extreme compositions."""
        f2_low = flory_huggins_second_derivative(0.01, 298.15, N_p=100)
        f2_high = flory_huggins_second_derivative(0.99, 298.15, N_p=100)
        assert f2_low > 0
        assert f2_high > 0

    def test_spinodal_exists_at_high_chi(self):
        """With large enough χ, spinodal decomposition should occur."""
        phi_s1, phi_s2 = spinodal_boundary(200.0, N_p=100, A=150.0, B=0.0)
        # χ = 150/200 = 0.75, which is above critical for N_p=100
        assert not np.isnan(phi_s1)
        assert phi_s1 < phi_s2

    def test_mesh_size_positive(self):
        xi = mesh_size_canal_peppas(0.05, 50000.0)
        assert xi > 0

    def test_mesh_size_decreases_with_crosslinking(self):
        """Higher crosslink density (lower Mc) → smaller mesh."""
        xi_high_mc = mesh_size_canal_peppas(0.05, 100000.0)
        xi_low_mc = mesh_size_canal_peppas(0.05, 10000.0)
        assert xi_low_mc < xi_high_mc


class TestPropertyDatabase:
    def test_default_construction(self):
        db = PropertyDatabase()
        assert db.props.rho_oil == 850.0

    def test_oil_density_decreases_with_T(self):
        db = PropertyDatabase()
        rho_cold = db.oil_density(293.15)
        rho_hot = db.oil_density(363.15)
        assert rho_hot < rho_cold

    def test_update_for_conditions(self):
        db = PropertyDatabase()
        props = db.update_for_conditions(363.15, 42.0, 18.0, 20.0)
        assert props.sigma < 0.050  # reduced by surfactant
        assert props.mu_d > 0
        assert props.rho_oil < 850.0  # less dense when hot
