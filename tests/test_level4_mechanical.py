"""Tests for Level 4 mechanical property prediction."""

import numpy as np
import pytest

from emulsim.level4_mechanical.solver import (
    agarose_modulus,
    double_network_modulus,
    effective_youngs_modulus,
    hertz_contact,
    ogston_kav,
)


class TestAgaroseModulus:
    def test_positive(self):
        G = agarose_modulus(42.0)
        assert G > 0

    def test_increases_with_concentration(self):
        G1 = agarose_modulus(20.0)
        G2 = agarose_modulus(42.0)
        assert G2 > G1

    def test_zero_at_zero_conc(self):
        assert agarose_modulus(0.0) == 0.0


class TestDoubleNetwork:
    def test_coupling_default(self):
        """With default eta=-0.15, G_DN < G1 + G2 (antagonistic coupling)."""
        G_DN = double_network_modulus(1000.0, 500.0)
        # G_DN = 1000 + 500 + (-0.15)*sqrt(1000*500) ≈ 1393.9
        expected = 1000.0 + 500.0 + (-0.15) * np.sqrt(1000.0 * 500.0)
        assert G_DN == pytest.approx(expected)
        assert G_DN < 1500.0  # antagonistic → less than additive

    def test_additive_when_eta_zero(self):
        """eta=0 recovers simple additivity."""
        G_DN = double_network_modulus(1000.0, 500.0, eta_coupling=0.0)
        assert G_DN == pytest.approx(1500.0)

    def test_synergistic(self):
        """eta > 0 gives more than additive."""
        G_DN = double_network_modulus(1000.0, 500.0, eta_coupling=0.3)
        assert G_DN > 1500.0

    def test_non_negative(self):
        """Result is clamped to non-negative."""
        G_DN = double_network_modulus(1.0, 1.0, eta_coupling=-10.0)
        assert G_DN >= 0.0


class TestYoungsModulus:
    def test_relationship(self):
        """E* = 2G(1+nu)/(1-nu^2). For nu=0.5: E*=2G*1.5/0.75=4G."""
        G = 1000.0
        E_star = effective_youngs_modulus(G, nu=0.5)
        assert E_star == pytest.approx(4.0 * G)

    def test_reduced_modulus_formula(self):
        """Verify E* = E/(1-nu^2) where E = 2G(1+nu)."""
        G = 1000.0
        nu = 0.45
        E_star = effective_youngs_modulus(G, nu)
        E = 2.0 * G * (1.0 + nu)
        assert E_star == pytest.approx(E / (1.0 - nu**2))


class TestHertzContact:
    def test_zero_force_at_zero_indent(self):
        delta, F = hertz_contact(1e4, 1e-6)
        assert F[0] == pytest.approx(0.0)

    def test_force_increases(self):
        delta, F = hertz_contact(1e4, 1e-6)
        assert np.all(np.diff(F) >= 0)

    def test_shape(self):
        delta, F = hertz_contact(1e4, 1e-6, n_points=50)
        assert len(delta) == 50
        assert len(F) == 50


class TestOgstonKav:
    def test_kav_one_for_point_solute(self):
        """Infinitely small solute → Kav approaches 1."""
        rh = np.array([1e-12])
        Kav = ogston_kav(rh, 1.5e-9, 0.05)
        assert Kav[0] > 0.9

    def test_kav_decreases_with_size(self):
        rh = np.array([1e-9, 5e-9, 10e-9, 20e-9])
        Kav = ogston_kav(rh, 1.5e-9, 0.05)
        assert np.all(np.diff(Kav) < 0)

    def test_kav_in_range(self):
        rh = np.logspace(-9, -7.5, 20)
        Kav = ogston_kav(rh, 1.5e-9, 0.3)
        assert np.all(Kav >= 0)
        assert np.all(Kav <= 1)
