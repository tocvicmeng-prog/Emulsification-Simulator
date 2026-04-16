"""Sprint 2: L1 dimensionless group calculations and monotonic trend tests.

Dimensionless group tests verify the computation is correct.
Monotonic trend tests verify that the PBE solver produces physically
correct directional behavior within its calibrated domain.
"""

from __future__ import annotations

import numpy as np
import pytest

from emulsim.level1_emulsification.validation import (
    DimensionlessGroups,
    compute_dimensionless_groups,
)


# ── DimensionlessGroups computation tests ────────────────────────────────


class TestDimensionlessGroups:

    def test_reynolds_number(self):
        """Re = rho * N * D^2 / mu. N=1000 rpm -> 16.67 rps, D=0.025 m, rho=850, mu=0.005."""
        g = compute_dimensionless_groups(
            rpm=1000, D_impeller=0.025, rho_c=850.0, mu_c=0.005,
            sigma=0.005, mu_d=1.0, Np=1.5, V_tank=5e-4,
        )
        # Re = 850 * 16.67 * 0.025^2 / 0.005 = 850 * 16.67 * 6.25e-4 / 0.005
        expected_Re = 850.0 * (1000 / 60.0) * 0.025**2 / 0.005
        assert abs(g.Re - expected_Re) / expected_Re < 1e-10

    def test_weber_number(self):
        """We = rho * N^2 * D^3 / sigma."""
        g = compute_dimensionless_groups(
            rpm=600, D_impeller=0.059, rho_c=850.0, mu_c=0.005,
            sigma=0.005, mu_d=0.1, Np=0.35, V_tank=5e-4,
        )
        N = 600 / 60.0
        expected_We = 850.0 * N**2 * 0.059**3 / 0.005
        assert abs(g.We - expected_We) / expected_We < 1e-10

    def test_capillary_number(self):
        """Ca = mu_c * v_tip / sigma."""
        g = compute_dimensionless_groups(
            rpm=1000, D_impeller=0.025, rho_c=850.0, mu_c=0.005,
            sigma=0.005, mu_d=1.0, Np=1.5, V_tank=5e-4,
        )
        v_tip = np.pi * 0.025 * (1000 / 60.0)
        expected_Ca = 0.005 * v_tip / 0.005
        assert abs(g.Ca - expected_Ca) / expected_Ca < 1e-10

    def test_viscosity_ratio(self):
        g = compute_dimensionless_groups(
            rpm=1000, D_impeller=0.025, rho_c=850.0, mu_c=0.005,
            sigma=0.005, mu_d=0.5, Np=1.5, V_tank=5e-4,
        )
        assert abs(g.viscosity_ratio - 100.0) < 1e-10

    def test_kolmogorov_length_positive(self):
        g = compute_dimensionless_groups(
            rpm=1000, D_impeller=0.025, rho_c=850.0, mu_c=0.005,
            sigma=0.005, mu_d=1.0, Np=1.5, V_tank=5e-4,
        )
        assert g.eta_K > 0

    def test_d32_ratio(self):
        g = compute_dimensionless_groups(
            rpm=1000, D_impeller=0.025, rho_c=850.0, mu_c=0.005,
            sigma=0.005, mu_d=1.0, Np=1.5, V_tank=5e-4,
            d32=50e-6,
        )
        assert g.d32_over_eta_K > 0

    def test_zero_rpm_returns_zeros(self):
        g = compute_dimensionless_groups(
            rpm=0, D_impeller=0.025, rho_c=850.0, mu_c=0.005,
            sigma=0.005, mu_d=1.0, Np=1.5, V_tank=5e-4,
        )
        assert g.Re == 0.0
        assert g.We == 0.0

    def test_as_dict(self):
        g = compute_dimensionless_groups(
            rpm=1000, D_impeller=0.025, rho_c=850.0, mu_c=0.005,
            sigma=0.005, mu_d=1.0, Np=1.5, V_tank=5e-4,
        )
        d = g.as_dict()
        assert "Re" in d
        assert "We" in d
        assert "eta_K_m" in d
        assert len(d) == 7


# ── Monotonic trend tests ────────────────────────────────────────────────
# These test physical sign constraints, not exact values.
# Each test runs the PBE solver at two parameter values and checks
# the expected direction of change.


class TestMonotonicTrends:
    """Physical monotonicity constraints for L1 PBE.

    These tests use the lightweight rotor-stator solver with short
    emulsification times to keep execution fast.
    """

    @staticmethod
    def _run_l1(rpm=10000, sigma=5e-3, mu_d=1.0, c_span80=20.0, t=10.0):
        """Run a quick L1 simulation and return d32."""
        from emulsim.datatypes import SimulationParameters, MaterialProperties
        from emulsim.level1_emulsification.solver import PBESolver

        params = SimulationParameters()
        params.emulsification.rpm = rpm
        params.emulsification.t_emulsification = t
        params.emulsification.mode = "rotor_stator_legacy"

        props = MaterialProperties()
        props.sigma = sigma
        props.mu_d = mu_d
        props.rho_oil = 850.0
        props.mu_oil = 0.005

        solver = PBESolver(n_bins=15, d_min=0.5e-6, d_max=200e-6)
        result = solver.solve(params, props, phi_d=0.05)
        return result.d32

    @pytest.mark.xfail(
        reason=(
            "F1 (doc 35): legacy rotor-stator PBE produces nonphysical RPM->d32 "
            "trend at mu_d=1 Pa.s, sigma=5 mN/m. Root cause is runaway CT "
            "coalescence at high droplet number density: high RPM produces many "
            "small droplets whose d_h^4 in the film-drainage damping collapses, "
            "letting coalescence pull d32 back up. Dispatch routing is now "
            "fixed (solver.py:251-274 honours KernelConfig), so the system is "
            "ready for kernel-constant calibration per docs/35 Phase 2 / "
            "Study A. Will pass after experimental DSD vs RPM data is fitted."
        )
    )
    def test_increasing_rpm_decreases_d32(self):
        """Higher agitation -> smaller droplets (within turbulent regime)."""
        d32_low = self._run_l1(rpm=5000)
        d32_high = self._run_l1(rpm=15000)
        assert d32_high <= d32_low, (
            f"d32 should decrease with RPM: d32(5k)={d32_low*1e6:.2f} um, "
            f"d32(15k)={d32_high*1e6:.2f} um"
        )

    @pytest.mark.xfail(
        reason=(
            "F1 (doc 35): legacy rotor-stator PBE produces nonphysical "
            "sigma->d32 trend with default uncalibrated CT/Alopaeus constants. "
            "Same calibration path as the RPM trend (Study A in docs/35 §9). "
            "Dispatch routing is now fixed; kernel constants need fitting."
        )
    )
    def test_increasing_sigma_increases_d32(self):
        """Higher IFT -> more resistance to breakage -> larger droplets."""
        d32_low_sigma = self._run_l1(sigma=2e-3)
        d32_high_sigma = self._run_l1(sigma=20e-3)
        assert d32_high_sigma >= d32_low_sigma, (
            f"d32 should increase with sigma: d32(2mN/m)={d32_low_sigma*1e6:.2f} um, "
            f"d32(20mN/m)={d32_high_sigma*1e6:.2f} um"
        )

    def test_increasing_viscosity_ratio_increases_d32(self):
        """Higher dispersed-phase viscosity -> harder to break -> larger droplets."""
        d32_low_mu = self._run_l1(mu_d=0.1)
        d32_high_mu = self._run_l1(mu_d=5.0)
        assert d32_high_mu >= d32_low_mu, (
            f"d32 should increase with mu_d: d32(0.1)={d32_low_mu*1e6:.2f} um, "
            f"d32(5.0)={d32_high_mu*1e6:.2f} um"
        )

    def test_d32_is_positive(self):
        """d32 should always be positive."""
        d32 = self._run_l1()
        assert d32 > 0
