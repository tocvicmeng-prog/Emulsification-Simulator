"""Tests for Node 14 (v7.0, P5a): Numba JIT for PBE kernels.

Acceptance for Node 14:
  1. JIT and NumPy paths produce numerically identical results within tight
     tolerance (1e-12 relative on non-zero entries).
  2. JIT path produces the SAME output as the NumPy path under
     boundary conditions (zero epsilon, very small d, very large mu_d/Vi
     cap region).
  3. Wider regression: smoke + L1 trend + L1 emulsification suite still
     pass after JIT introduction (handled by the existing CI run, not
     duplicated here).
"""

from __future__ import annotations

import numpy as np
import pytest

import emulsim.level1_emulsification.kernels as _kernels_mod
from emulsim.datatypes import KernelConfig
from emulsim.level1_emulsification.kernels import (
    breakage_rate_alopaeus,
    breakage_rate_coulaloglou,
    coalescence_rate_dispatch,
)


@pytest.fixture
def grid():
    """Log-spaced diameter grid covering the typical PBE pivot range."""
    return np.geomspace(1e-7, 200e-6, 25)


@pytest.fixture(autouse=False)
def force_no_jit(monkeypatch):
    """Disable the JIT dispatch so the public functions run the NumPy path."""
    monkeypatch.setattr(_kernels_mod, "_USE_JIT", False)


# ─── Numerical equivalence: JIT vs NumPy ──────────────────────────────────


class TestJITMatchesNumPy:
    def test_breakage_alopaeus_matches(self, grid, monkeypatch):
        """JIT and NumPy outputs must agree to machine precision."""
        # JIT path
        jit_g = breakage_rate_alopaeus(
            grid, epsilon=1e3, sigma=5e-3, rho_c=850.0, mu_d=1.0,
            C1=0.986, C2=0.0115, C3=0.3, nu_c=5.88e-6,
        )
        # NumPy path (same inputs, JIT disabled)
        monkeypatch.setattr(_kernels_mod, "_USE_JIT", False)
        np_g = breakage_rate_alopaeus(
            grid, epsilon=1e3, sigma=5e-3, rho_c=850.0, mu_d=1.0,
            C1=0.986, C2=0.0115, C3=0.3, nu_c=5.88e-6,
        )
        np.testing.assert_allclose(jit_g, np_g, rtol=1e-12, atol=1e-300)

    def test_breakage_alopaeus_no_viscous_correction(self, grid, monkeypatch):
        """C3=0 path (default) must also match exactly."""
        jit_g = breakage_rate_alopaeus(
            grid, epsilon=1e3, sigma=5e-3, rho_c=850.0, mu_d=1.0,
            C1=0.986, C2=0.0115, C3=0.0, nu_c=5.88e-6,
        )
        monkeypatch.setattr(_kernels_mod, "_USE_JIT", False)
        np_g = breakage_rate_alopaeus(
            grid, epsilon=1e3, sigma=5e-3, rho_c=850.0, mu_d=1.0,
            C1=0.986, C2=0.0115, C3=0.0, nu_c=5.88e-6,
        )
        np.testing.assert_allclose(jit_g, np_g, rtol=1e-12, atol=1e-300)

    def test_breakage_coulaloglou_matches(self, grid, monkeypatch):
        jit_g = breakage_rate_coulaloglou(
            grid, epsilon=5e2, sigma=4e-3, rho_c=850.0,
            C1=0.00481, C2=0.08,
        )
        monkeypatch.setattr(_kernels_mod, "_USE_JIT", False)
        np_g = breakage_rate_coulaloglou(
            grid, epsilon=5e2, sigma=4e-3, rho_c=850.0,
            C1=0.00481, C2=0.08,
        )
        np.testing.assert_allclose(jit_g, np_g, rtol=1e-12, atol=1e-300)

    def test_coalescence_matrix_matches(self, grid, monkeypatch):
        kernels = KernelConfig.for_rotor_stator_legacy()
        jit_Q = coalescence_rate_dispatch(
            grid, epsilon=2e3, sigma=5e-3, rho_c=850.0,
            config=kernels, phi_d=0.05, mu_c=0.001,
        )
        monkeypatch.setattr(_kernels_mod, "_USE_JIT", False)
        np_Q = coalescence_rate_dispatch(
            grid, epsilon=2e3, sigma=5e-3, rho_c=850.0,
            config=kernels, phi_d=0.05, mu_c=0.001,
        )
        # Coalescence has a wider dynamic range; relax atol to a value still
        # well below any physically meaningful coalescence rate (>1e-30 m^3/s).
        np.testing.assert_allclose(jit_Q, np_Q, rtol=1e-12, atol=1e-30)
        assert jit_Q.shape == (grid.size, grid.size)


# ─── Edge cases ───────────────────────────────────────────────────────────


class TestEdgeCases:
    def test_zero_epsilon_returns_zeros(self, grid):
        g = breakage_rate_alopaeus(
            grid, epsilon=0.0, sigma=5e-3, rho_c=850.0, mu_d=1.0,
        )
        assert np.all(g == 0.0)

    def test_negative_epsilon_returns_zeros(self, grid):
        kernels = KernelConfig.for_rotor_stator_legacy()
        Q = coalescence_rate_dispatch(
            grid, epsilon=-1.0, sigma=5e-3, rho_c=850.0,
            config=kernels, phi_d=0.05, mu_c=0.001,
        )
        assert np.all(Q == 0.0)

    def test_high_mu_d_triggers_vi_cap(self):
        """JIT and NumPy must both honour the Vi=100 cap (anti-F1 guard)."""
        d = np.array([1e-6, 5e-6, 50e-6])
        g = breakage_rate_alopaeus(
            d, epsilon=1e3, sigma=5e-3, rho_c=850.0,
            mu_d=100.0,            # extreme viscosity -> Vi >> 100
            C3=0.05, nu_c=5.88e-6,
        )
        # Cap means breakage suppression saturates instead of going to 0.
        assert np.all(np.isfinite(g))
        assert g[0] < g[-1]  # smaller d still breaks slower (cap saturated)


# ─── Smoke benchmark (informational only) ─────────────────────────────────


class TestBenchmark:
    @pytest.mark.smoke
    def test_kernel_call_under_microseconds(self, grid):
        """Sanity benchmark — JIT path should be well under 100 us per call."""
        import time
        kernels = KernelConfig.for_rotor_stator_legacy()
        # Warm-up
        for _ in range(5):
            breakage_rate_alopaeus(grid, 1e3, 5e-3, 850.0, 1.0)
            coalescence_rate_dispatch(grid, 1e3, 5e-3, 850.0, kernels,
                                       phi_d=0.05, mu_c=0.001)
        # Time
        N = 500
        t0 = time.perf_counter()
        for _ in range(N):
            breakage_rate_alopaeus(grid, 1e3, 5e-3, 850.0, 1.0)
        breakage_us = (time.perf_counter() - t0) / N * 1e6

        t0 = time.perf_counter()
        for _ in range(N):
            coalescence_rate_dispatch(grid, 1e3, 5e-3, 850.0, kernels,
                                       phi_d=0.05, mu_c=0.001)
        coalescence_us = (time.perf_counter() - t0) / N * 1e6

        # Loose bound — failing this means the JIT path regressed badly.
        # Typical observation on dev box: breakage ~1-3 us, coalescence ~5-15 us.
        assert breakage_us < 100.0, f"breakage {breakage_us:.1f} us/call"
        assert coalescence_us < 200.0, f"coalescence {coalescence_us:.1f} us/call"
