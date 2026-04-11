"""Tests for Module 2 Phase A: ACS data model and surface area core.

Validates:
  - Surface area geometry (external, internal, accessibility)
  - ACS conservation hierarchy
  - M1 export contract integration
  - Synthetic multi-step modification sequences
"""

from __future__ import annotations

import math

import pytest

from emulsim.datatypes import M1ExportContract
from emulsim.module2_functionalization.acs import (
    ACSSiteType,
    ACSProfile,
    initialize_acs_from_m1,
)
from emulsim.module2_functionalization.surface_area import (
    AccessibleSurfaceModel,
    SurfaceAreaTier,
)


# ─── Fixtures ───────────────────────────────────────────────────────────

def _make_contract(
    bead_radius: float = 50e-6,
    porosity: float = 0.7,
    pore_size_mean: float = 100e-9,
    mesh_size_xi: float = 20e-9,
    nh2_bulk: float = 100.0,
    oh_bulk: float = 400.0,
    p_final: float = 0.5,
) -> M1ExportContract:
    """Helper: create a minimal M1ExportContract for testing."""
    return M1ExportContract(
        bead_radius=bead_radius,
        bead_d32=bead_radius * 2,
        bead_d50=bead_radius * 2,
        pore_size_mean=pore_size_mean,
        pore_size_std=pore_size_mean * 0.3,
        porosity=porosity,
        l2_model_tier="empirical_calibrated",
        mesh_size_xi=mesh_size_xi,
        p_final=p_final,
        primary_crosslinker="genipin",
        nh2_bulk_concentration=nh2_bulk,
        oh_bulk_concentration=oh_bulk,
        G_DN=5000.0,
        E_star=15000.0,
        model_used="phenomenological",
        c_agarose=42.0,
        c_chitosan=18.0,
        DDA=0.90,
        trust_level="CAUTION",
    )


# ─── Surface Area Tests ────────────────────────────────────────────────

class TestSurfaceAreaExternalOnly:
    """External-only tier: A = 4*pi*R^2, no internal contribution."""

    def test_external_area_value(self):
        R = 50e-6  # 50 um
        model = AccessibleSurfaceModel(
            tier=SurfaceAreaTier.EXTERNAL_ONLY,
            bead_radius=R,
        )
        model.compute()

        A_expected = 4.0 * math.pi * R ** 2
        assert math.isclose(model.external_area, A_expected, rel_tol=1e-10)

    def test_internal_area_is_zero(self):
        model = AccessibleSurfaceModel(
            tier=SurfaceAreaTier.EXTERNAL_ONLY,
            bead_radius=50e-6,
        )
        model.compute()
        assert model.internal_geometric_area == 0.0

    def test_accessible_equals_external(self):
        model = AccessibleSurfaceModel(
            tier=SurfaceAreaTier.EXTERNAL_ONLY,
            bead_radius=50e-6,
        )
        model.compute()
        assert model.reagent_accessible_area == model.external_area
        assert model.ligand_accessible_area == model.external_area

    def test_trust_level_unreliable(self):
        model = AccessibleSurfaceModel(tier=SurfaceAreaTier.EXTERNAL_ONLY, bead_radius=50e-6)
        model.compute()
        assert model.trust_level == "UNRELIABLE"


class TestSurfaceAreaEmpiricalPore:
    """Empirical pore tier: S_v = 4*eps/d_pore, total = A_ext + A_int."""

    def test_internal_area_formula(self):
        R = 50e-6
        eps = 0.7
        d_pore = 100e-9

        model = AccessibleSurfaceModel(
            tier=SurfaceAreaTier.EMPIRICAL_PORE,
            bead_radius=R,
            porosity=eps,
            pore_diameter_mean=d_pore,
        )
        model.compute()

        V_bead = (4.0 / 3.0) * math.pi * R ** 3
        S_v = 4.0 * eps / d_pore
        A_int_expected = S_v * V_bead

        assert math.isclose(model.internal_geometric_area, A_int_expected, rel_tol=1e-10)

    def test_external_area_still_correct(self):
        R = 50e-6
        model = AccessibleSurfaceModel(
            tier=SurfaceAreaTier.EMPIRICAL_PORE,
            bead_radius=R,
            porosity=0.7,
            pore_diameter_mean=100e-9,
        )
        model.compute()
        A_ext = 4.0 * math.pi * R ** 2
        assert math.isclose(model.external_area, A_ext, rel_tol=1e-10)

    def test_trust_level_caution(self):
        model = AccessibleSurfaceModel(
            tier=SurfaceAreaTier.EMPIRICAL_PORE,
            bead_radius=50e-6,
            porosity=0.7,
            pore_diameter_mean=100e-9,
        )
        model.compute()
        assert model.trust_level == "CAUTION"


class TestSurfaceAreaAccessibility:
    """Larger solutes have lower accessible area."""

    def test_reagent_more_accessible_than_ligand(self):
        model = AccessibleSurfaceModel(
            tier=SurfaceAreaTier.EMPIRICAL_PORE,
            bead_radius=50e-6,
            porosity=0.7,
            pore_diameter_mean=100e-9,
        )
        model.compute(reagent_radius=0.5e-9, ligand_radius=3e-9)
        assert model.reagent_accessible_area > model.ligand_accessible_area

    def test_very_large_solute_excluded(self):
        """Solute larger than pore diameter has no internal access."""
        model = AccessibleSurfaceModel(
            tier=SurfaceAreaTier.EMPIRICAL_PORE,
            bead_radius=50e-6,
            porosity=0.7,
            pore_diameter_mean=10e-9,  # 10 nm pores
        )
        # Ligand radius 10 nm => 2*r_h = 20 nm > d_pore = 10 nm
        model.compute(reagent_radius=0.5e-9, ligand_radius=10e-9)
        # Ligand-accessible area should be just external
        assert math.isclose(
            model.ligand_accessible_area, model.external_area, rel_tol=1e-10
        )

    def test_higher_tortuosity_reduces_access(self):
        model_low = AccessibleSurfaceModel(
            tier=SurfaceAreaTier.EMPIRICAL_PORE,
            bead_radius=50e-6,
            porosity=0.7,
            pore_diameter_mean=100e-9,
        )
        model_high = AccessibleSurfaceModel(
            tier=SurfaceAreaTier.EMPIRICAL_PORE,
            bead_radius=50e-6,
            porosity=0.7,
            pore_diameter_mean=100e-9,
        )
        model_low.compute(tortuosity=1.0)
        model_high.compute(tortuosity=3.0)
        assert model_low.reagent_accessible_area > model_high.reagent_accessible_area


class TestSurfaceAreaFromM1Export:
    """Factory from M1ExportContract populates correctly."""

    def test_factory_populates_areas(self):
        contract = _make_contract()
        model = AccessibleSurfaceModel.from_m1_export(contract)
        assert model.external_area > 0
        assert model.internal_geometric_area > 0
        assert model.reagent_accessible_area > model.external_area

    def test_factory_uses_contract_geometry(self):
        contract = _make_contract(bead_radius=25e-6, porosity=0.5, pore_size_mean=200e-9)
        model = AccessibleSurfaceModel.from_m1_export(contract)
        assert math.isclose(model.bead_radius, 25e-6)
        assert math.isclose(model.porosity, 0.5)
        assert math.isclose(model.pore_diameter_mean, 200e-9)

    def test_factory_morphology_fallback(self):
        """MORPHOLOGY_BASED falls back to EMPIRICAL_PORE for non-mechanistic L2."""
        contract = _make_contract()  # l2_model_tier = "empirical_calibrated"
        model = AccessibleSurfaceModel.from_m1_export(
            contract, tier=SurfaceAreaTier.MORPHOLOGY_BASED
        )
        assert model.tier == SurfaceAreaTier.EMPIRICAL_PORE


# ─── ACS Profile Tests ─────────────────────────────────────────────────

class TestACSProfileValidation:
    """ACS conservation hierarchy validation."""

    def test_validate_passes_valid_profile(self):
        profile = ACSProfile(
            site_type=ACSSiteType.AMINE_PRIMARY,
            total_sites=1e-9,
            accessible_sites=0.8e-9,
            activated_sites=0.5e-9,
            crosslinked_sites=0.2e-9,
            blocked_sites=0.1e-9,
            ligand_coupled_sites=0.15e-9,
            ligand_functional_sites=0.10e-9,
        )
        errors = profile.validate()
        assert errors == [], f"Unexpected violations: {errors}"

    def test_validate_catches_accessible_gt_total(self):
        profile = ACSProfile(
            site_type=ACSSiteType.AMINE_PRIMARY,
            total_sites=1e-9,
            accessible_sites=2e-9,  # violation
        )
        errors = profile.validate()
        assert any("accessible" in e and "total" in e for e in errors)

    def test_validate_catches_activated_gt_accessible(self):
        profile = ACSProfile(
            site_type=ACSSiteType.HYDROXYL,
            total_sites=1e-9,
            accessible_sites=0.8e-9,
            activated_sites=0.9e-9,  # violation
        )
        errors = profile.validate()
        assert any("activated" in e and "accessible" in e for e in errors)

    def test_validate_catches_ligand_coupled_gt_activated(self):
        profile = ACSProfile(
            site_type=ACSSiteType.AMINE_PRIMARY,
            total_sites=1e-9,
            accessible_sites=0.8e-9,
            activated_sites=0.5e-9,
            ligand_coupled_sites=0.6e-9,  # violation
        )
        errors = profile.validate()
        assert any("ligand_coupled" in e and "activated" in e for e in errors)

    def test_validate_catches_ligand_functional_gt_coupled(self):
        profile = ACSProfile(
            site_type=ACSSiteType.AMINE_PRIMARY,
            total_sites=1e-9,
            accessible_sites=0.8e-9,
            activated_sites=0.5e-9,
            ligand_coupled_sites=0.3e-9,
            ligand_functional_sites=0.4e-9,  # violation
        )
        errors = profile.validate()
        assert any("ligand_functional" in e and "ligand_coupled" in e for e in errors)

    def test_validate_catches_negative_sites(self):
        profile = ACSProfile(
            site_type=ACSSiteType.HYDROXYL,
            total_sites=-1e-9,  # violation
        )
        errors = profile.validate()
        assert any("negative" in e for e in errors)


class TestACSRemainingSites:
    """remaining = accessible - consumed - blocked."""

    def test_remaining_basic(self):
        profile = ACSProfile(
            site_type=ACSSiteType.AMINE_PRIMARY,
            accessible_sites=1e-9,
            crosslinked_sites=0.3e-9,
            blocked_sites=0.2e-9,
        )
        expected = 1e-9 - 0.3e-9 - 0.2e-9
        assert math.isclose(profile.remaining_sites, expected, rel_tol=1e-10)

    def test_remaining_clamped_to_zero(self):
        profile = ACSProfile(
            site_type=ACSSiteType.AMINE_PRIMARY,
            accessible_sites=1e-9,
            crosslinked_sites=0.6e-9,
            blocked_sites=0.6e-9,  # sum > accessible
        )
        assert profile.remaining_sites == 0.0


# ─── Initialize from M1 ────────────────────────────────────────────────

class TestInitializeACSFromM1:
    """Full initialization from M1 contract produces valid ACS profiles."""

    def test_produces_both_site_types(self):
        contract = _make_contract()
        surface = AccessibleSurfaceModel.from_m1_export(contract)
        profiles = initialize_acs_from_m1(contract, surface)

        assert ACSSiteType.AMINE_PRIMARY in profiles
        assert ACSSiteType.HYDROXYL in profiles

    def test_profiles_pass_validation(self):
        contract = _make_contract()
        surface = AccessibleSurfaceModel.from_m1_export(contract)
        profiles = initialize_acs_from_m1(contract, surface)

        for site_type, profile in profiles.items():
            errors = profile.validate()
            assert errors == [], f"{site_type.value} violations: {errors}"

    def test_total_sites_formula(self):
        """n_total = C_bulk * V_bead."""
        R = 50e-6
        nh2_bulk = 100.0  # mol/m^3
        contract = _make_contract(bead_radius=R, nh2_bulk=nh2_bulk)
        surface = AccessibleSurfaceModel.from_m1_export(contract)
        profiles = initialize_acs_from_m1(contract, surface)

        V_bead = (4.0 / 3.0) * math.pi * R ** 3
        expected_total = nh2_bulk * V_bead

        nh2 = profiles[ACSSiteType.AMINE_PRIMARY]
        assert math.isclose(nh2.total_sites, expected_total, rel_tol=1e-10)

    def test_accessible_less_equal_total(self):
        contract = _make_contract()
        surface = AccessibleSurfaceModel.from_m1_export(contract)
        profiles = initialize_acs_from_m1(contract, surface)

        for profile in profiles.values():
            assert profile.accessible_sites <= profile.total_sites * 1.001

    def test_surface_densities_positive(self):
        contract = _make_contract()
        surface = AccessibleSurfaceModel.from_m1_export(contract)
        profiles = initialize_acs_from_m1(contract, surface)

        for profile in profiles.values():
            assert profile.total_density > 0
            assert profile.accessible_density > 0

    def test_zero_concentration_skipped(self):
        contract = _make_contract(nh2_bulk=0.0)
        surface = AccessibleSurfaceModel.from_m1_export(contract)
        profiles = initialize_acs_from_m1(contract, surface)

        assert ACSSiteType.AMINE_PRIMARY not in profiles
        assert ACSSiteType.HYDROXYL in profiles

    def test_raises_on_uncomputed_surface(self):
        contract = _make_contract()
        surface = AccessibleSurfaceModel()  # not computed
        with pytest.raises(ValueError, match="compute"):
            initialize_acs_from_m1(contract, surface)


# ─── Conservation Through Synthetic Sequence ───────────────────────────

class TestACSConservationSyntheticSequence:
    """Synthetic 3-step modification preserves conservation.

    Step 1: Consume 30% of NH2 (secondary crosslinking with glutaraldehyde).
    Step 2: Activate 50% of OH (ECH epoxidation).
    Step 3: Block remaining activated OH sites.

    At every step, validate() must pass.
    """

    def test_three_step_conservation(self):
        # --- Initialize ---
        contract = _make_contract(nh2_bulk=100.0, oh_bulk=400.0)
        surface = AccessibleSurfaceModel.from_m1_export(contract)
        profiles = initialize_acs_from_m1(contract, surface)

        nh2 = profiles[ACSSiteType.AMINE_PRIMARY]
        oh = profiles[ACSSiteType.HYDROXYL]

        # --- Step 1: Consume 30% of accessible NH2 (secondary crosslinking) ---
        consumed_fraction = 0.30
        nh2.crosslinked_sites = consumed_fraction * nh2.accessible_sites

        errors = nh2.validate()
        assert errors == [], f"Step 1 NH2 violations: {errors}"
        # remaining + terminal_sum = accessible
        balance = nh2.remaining_sites + nh2._terminal_sum
        assert math.isclose(balance, nh2.accessible_sites, rel_tol=1e-10)

        # --- Step 2: Activate 50% of accessible OH (ECH epoxidation) ---
        activation_fraction = 0.50
        oh.activated_sites = activation_fraction * oh.accessible_sites

        errors = oh.validate()
        assert errors == [], f"Step 2 OH violations: {errors}"

        # --- Step 3: Block all remaining activated OH sites (quenching) ---
        oh.blocked_sites = oh.activated_sites  # quench all activated — only blocked_sites

        errors = oh.validate()
        assert errors == [], f"Step 3 OH violations: {errors}"
        # Terminal sum check: terminal_sum <= accessible
        assert oh._terminal_sum <= oh.accessible_sites * 1.001

    def test_ligand_coupling_sub_sequence(self):
        """After activation, couple ligand to 60% of activated sites,
        with 80% retaining function."""
        contract = _make_contract()
        surface = AccessibleSurfaceModel.from_m1_export(contract)
        profiles = initialize_acs_from_m1(contract, surface)
        oh = profiles[ACSSiteType.HYDROXYL]

        # Activate 50%
        oh.activated_sites = 0.50 * oh.accessible_sites

        # Couple ligand to 60% of activated
        oh.ligand_coupled_sites = 0.60 * oh.activated_sites

        # 80% of coupled ligand retains activity
        oh.ligand_functional_sites = 0.80 * oh.ligand_coupled_sites

        errors = oh.validate()
        assert errors == [], f"Ligand coupling violations: {errors}"

        # Functional <= coupled <= activated <= accessible <= total
        assert oh.ligand_functional_sites <= oh.ligand_coupled_sites
        assert oh.ligand_coupled_sites <= oh.activated_sites
        assert oh.activated_sites <= oh.accessible_sites
        assert oh.accessible_sites <= oh.total_sites


# ─── Physical Plausibility ──────────────────────────────────────────────

class TestSurfaceAreaValuesPhysical:
    """For 50um bead with 100nm pores and 70% porosity:
    A_ext ~ 3.14e-8 m^2, A_int ~ 1.47e-5 m^2.
    Internal >> external (as expected for macroporous beads)."""

    def test_orders_of_magnitude(self):
        R = 50e-6
        model = AccessibleSurfaceModel(
            tier=SurfaceAreaTier.EMPIRICAL_PORE,
            bead_radius=R,
            porosity=0.7,
            pore_diameter_mean=100e-9,
        )
        model.compute()

        # A_ext = 4*pi*(50e-6)^2 = 3.14e-8 m^2
        assert 3.0e-8 < model.external_area < 3.5e-8

        # S_v = 4*0.7/100e-9 = 2.8e7 m^2/m^3
        # V_bead = (4/3)*pi*(50e-6)^3 = 5.24e-13 m^3
        # A_int = 2.8e7 * 5.24e-13 = 1.47e-5 m^2
        assert 1.0e-5 < model.internal_geometric_area < 2.0e-5

        # Internal >> external (ratio ~ 470x)
        ratio = model.internal_geometric_area / model.external_area
        assert ratio > 100, f"Expected internal >> external, got ratio {ratio:.1f}"

    def test_bead_volume_property(self):
        R = 50e-6
        model = AccessibleSurfaceModel(bead_radius=R)
        V_expected = (4.0 / 3.0) * math.pi * R ** 3
        assert math.isclose(model.bead_volume, V_expected, rel_tol=1e-10)


# ─── Edge Cases ─────────────────────────────────────────────────────────

class TestEdgeCases:
    """Input validation edge cases."""

    def test_negative_radius_raises(self):
        model = AccessibleSurfaceModel(bead_radius=-1e-6)
        with pytest.raises(ValueError, match="bead_radius"):
            model.compute()

    def test_porosity_out_of_range_raises(self):
        model = AccessibleSurfaceModel(porosity=1.5)
        with pytest.raises(ValueError, match="porosity"):
            model.compute()

    def test_tortuosity_below_one_raises(self):
        model = AccessibleSurfaceModel(
            tier=SurfaceAreaTier.EMPIRICAL_PORE,
            bead_radius=50e-6,
            porosity=0.7,
            pore_diameter_mean=100e-9,
        )
        with pytest.raises(ValueError, match="tortuosity"):
            model.compute(tortuosity=0.5)
