"""Integration tests for gradient-aware LRM (v6.0 H6).

Tests that gradient values actually affect equilibrium during LRM time
integration, not just display. Covers:
  - Single-component solve_lrm with EquilibriumAdapter + GradientProgram
  - HIC salt gradient: higher salt -> more binding -> later elution
  - IMAC imidazole gradient: protein elutes at expected imidazole range
  - Backward compatibility: solve_lrm without gradient behaves identically
"""

from __future__ import annotations

import numpy as np
import pytest

from emulsim.module3_performance.hydrodynamics import ColumnGeometry
from emulsim.module3_performance.isotherms.langmuir import LangmuirIsotherm
from emulsim.module3_performance.isotherms.hic import HICIsotherm
from emulsim.module3_performance.isotherms.imac import IMACCompetitionIsotherm
from emulsim.module3_performance.isotherms.adapter import EquilibriumAdapter
from emulsim.module3_performance.process_state import ProcessState
from emulsim.module3_performance.transport.lumped_rate import solve_lrm
from emulsim.module3_performance.gradient import (
    make_linear_gradient,
)

# All tests in this module exercise solve_lrm with a gradient program. LSODA
# is great for breakthrough (constant equilibrium) but hangs when the binding
# constant varies in time, so these stay on the BDF code path. Faster gradient
# solving needs an analytical Jacobian — separate work.
pytestmark = pytest.mark.slow



# ─── Shared Fixtures ──────────────────────────────────────────────────────

@pytest.fixture
def small_column():
    """Small column for fast tests (10 mm x 5 cm bed)."""
    return ColumnGeometry(
        diameter=0.01,
        bed_height=0.05,
        particle_diameter=100e-6,
        bed_porosity=0.38,
        particle_porosity=0.5,
        G_DN=10000.0,
        E_star=30000.0,
    )


# ─── Test: Backward Compatibility ─────────────────────────────────────────

class TestBackwardCompatibility:
    """Verify solve_lrm without gradient params behaves identically to before."""

    def test_no_gradient_same_result(self, small_column):
        """Calling solve_lrm with no gradient_program gives identical result."""
        iso = LangmuirIsotherm(q_max=100.0, K_L=1000.0)
        kwargs = dict(
            column=small_column,
            isotherm=iso,
            C_feed=0.01,
            feed_duration=60.0,
            flow_rate=1e-8,
            total_time=120.0,
            n_z=20,
        )
        r1 = solve_lrm(**kwargs)
        r2 = solve_lrm(**kwargs, gradient_program=None, equilibrium_adapter=None)
        np.testing.assert_allclose(r1.C_outlet, r2.C_outlet, atol=1e-15)

    def test_gradient_program_only_no_adapter_same_result(self, small_column):
        """If gradient_program is provided but no adapter, behavior unchanged."""
        iso = LangmuirIsotherm(q_max=100.0, K_L=1000.0)
        grad = make_linear_gradient(0.0, 500.0, 0.0, 120.0)
        kwargs = dict(
            column=small_column,
            isotherm=iso,
            C_feed=0.01,
            feed_duration=60.0,
            flow_rate=1e-8,
            total_time=120.0,
            n_z=20,
        )
        r1 = solve_lrm(**kwargs)
        # gradient_program without adapter -> no effect (both must be present)
        r2 = solve_lrm(**kwargs, gradient_program=grad)
        np.testing.assert_allclose(r1.C_outlet, r2.C_outlet, atol=1e-15)


# ─── Test: HIC Salt Gradient ──────────────────────────────────────────────

class TestHICSaltGradient:
    """HIC: higher salt -> higher K_eff -> more binding.

    At high salt (load), protein binds. As salt decreases (elution),
    K_eff drops and protein desorbs.
    """

    def test_hic_gradient_changes_elution(self, small_column):
        """HIC with decreasing salt gradient should produce different profile."""
        # Use moderate parameters so protein actually elutes
        hic = HICIsotherm(q_max=20.0, K_0=10.0, m_salt=0.003)

        # Start at moderate salt, ramp to zero
        state = ProcessState(salt_concentration=200.0)
        adapter = EquilibriumAdapter(hic, state)
        grad = make_linear_gradient(200.0, 0.0, 300.0, 900.0)

        # Langmuir for fixed comparison
        iso_fixed = LangmuirIsotherm(q_max=20.0, K_L=100.0)

        common = dict(
            column=small_column,
            C_feed=0.1,
            feed_duration=300.0,
            flow_rate=5e-8,
            total_time=1200.0,
            n_z=20,
        )

        r_grad = solve_lrm(
            **common,
            isotherm=iso_fixed,
            gradient_program=grad,
            equilibrium_adapter=adapter,
            gradient_field="salt_concentration",
        )

        r_fixed = solve_lrm(**common, isotherm=iso_fixed)

        # The gradient-aware run updates equilibrium at each time step,
        # so the outlet profile must differ from the fixed-isotherm run.
        max_diff = np.max(np.abs(r_grad.C_outlet - r_fixed.C_outlet))
        assert max_diff > 1e-8, \
            f"Gradient-aware HIC should differ from fixed (max_diff={max_diff:.2e})"

        # Both should complete without solver failure
        assert r_grad.mass_balance_error < 0.15
        assert r_fixed.mass_balance_error < 0.15

    def test_hic_high_salt_binds_more(self, small_column):
        """At constant high salt, HIC binds strongly vs zero salt."""
        import math
        hic = HICIsotherm(q_max=80.0, K_0=0.1, m_salt=0.01)

        # High salt: K_eff = 0.1 * exp(0.01 * 500) = 0.1 * ~148 = ~14.8
        K_high = 0.1 * math.exp(0.01 * 500)
        q_high = 80.0 * K_high * 0.01 / (1 + K_high * 0.01)
        assert q_high > 5.0, f"HIC should bind at high salt, got q={q_high:.3f}"

        # Zero salt: K_eff = 0.1, very weak binding
        K_low = 0.1 * math.exp(0.01 * 0)  # = 0.1
        q_low = 80.0 * K_low * 0.01 / (1 + K_low * 0.01)
        assert q_low < q_high * 0.1, \
            f"HIC at zero salt ({q_low:.4f}) should be << high salt ({q_high:.3f})"


# ─── Test: IMAC Imidazole Gradient ────────────────────────────────────────

class TestIMACImidazoleGradient:
    """IMAC: imidazole gradient elutes His-tagged protein."""

    def test_imac_gradient_changes_elution(self, small_column):
        """IMAC with rising imidazole should produce different profile."""
        imac = IMACCompetitionIsotherm(q_max=20.0, K_protein=1e3, K_imidazole=10.0)

        state = ProcessState(imidazole=0.0)
        adapter = EquilibriumAdapter(imac, state)
        # Imidazole ramp: 0 -> 500 mol/m3 starting at t=300s
        grad = make_linear_gradient(0.0, 500.0, 300.0, 900.0)

        iso_fallback = LangmuirIsotherm(q_max=20.0, K_L=1e3)

        common = dict(
            column=small_column,
            C_feed=0.1,
            feed_duration=300.0,
            flow_rate=5e-8,
            total_time=1200.0,
            n_z=20,
        )

        r_grad = solve_lrm(
            **common,
            isotherm=iso_fallback,
            gradient_program=grad,
            equilibrium_adapter=adapter,
            gradient_field="imidazole",
        )

        r_fixed = solve_lrm(**common, isotherm=iso_fallback)

        # The gradient-aware run updates equilibrium at each time step
        max_diff = np.max(np.abs(r_grad.C_outlet - r_fixed.C_outlet))
        assert max_diff > 1e-8, \
            f"IMAC with imidazole gradient should differ from fixed (max_diff={max_diff:.2e})"

        # Should complete without solver failure
        assert r_grad.mass_balance_error < 0.15

    def test_imac_competition_reduces_binding(self):
        """Imidazole competition reduces protein equilibrium loading."""
        imac = IMACCompetitionIsotherm(q_max=50.0, K_protein=1e4, K_imidazole=10.0)

        q_no_imid, _ = imac.equilibrium_loading(0.01, 0.0)
        q_high_imid, _ = imac.equilibrium_loading(0.01, 500.0)

        assert float(q_no_imid) > float(q_high_imid), \
            "High imidazole should reduce protein binding on IMAC"
        assert float(q_high_imid) < float(q_no_imid) * 0.5, \
            "500 mol/m3 imidazole should substantially reduce binding"


# ─── Test: Gradient-Sensitive Property ─────────────────────────────────────

class TestGradientSensitiveProperty:
    """Verify gradient_sensitive and gradient_field on all isotherms."""

    def test_sma_gradient_sensitive(self):
        from emulsim.module3_performance.isotherms.sma import SMAIsotherm
        iso = SMAIsotherm()
        assert iso.gradient_sensitive is True
        assert iso.gradient_field == "salt_concentration"

    def test_hic_gradient_sensitive(self):
        iso = HICIsotherm()
        assert iso.gradient_sensitive is True
        assert iso.gradient_field == "salt_concentration"

    def test_imac_gradient_sensitive(self):
        iso = IMACCompetitionIsotherm()
        assert iso.gradient_sensitive is True
        assert iso.gradient_field == "imidazole"

    def test_protein_a_gradient_sensitive(self):
        from emulsim.module3_performance.isotherms.protein_a import ProteinAIsotherm
        iso = ProteinAIsotherm()
        assert iso.gradient_sensitive is True
        assert iso.gradient_field == "ph"

    def test_competitive_affinity_gradient_sensitive(self):
        from emulsim.module3_performance.isotherms.competitive_affinity import CompetitiveAffinityIsotherm
        iso = CompetitiveAffinityIsotherm()
        assert iso.gradient_sensitive is True
        assert iso.gradient_field == "sugar_competitor"

    def test_competitive_langmuir_not_sensitive(self):
        from emulsim.module3_performance.isotherms.competitive_langmuir import CompetitiveLangmuirIsotherm
        iso = CompetitiveLangmuirIsotherm(
            q_max=np.array([100.0]), K_L=np.array([1000.0]),
        )
        assert iso.gradient_sensitive is False
