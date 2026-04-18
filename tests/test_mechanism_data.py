"""Tests for emulsim.protocols.mechanism_data.

Covers:
  - MECHANISM_REGISTRY size and key validity
  - MechanismDescriptor field correctness for the genipin entry
  - ReactionStep field correctness for genipin step 1
  - get_mechanism() known key, fallback paths, and unknown key robustness
  - Full coverage: every CROSSLINKER and REAGENT_PROFILE key resolves
"""

from __future__ import annotations


from emulsim.protocols.mechanism_data import (
    MECHANISM_REGISTRY,
    MechanismDescriptor,
    ReactionStep,
    get_mechanism,
)
from emulsim.reagent_library import CROSSLINKERS
from emulsim.module2_functionalization.reagent_profiles import REAGENT_PROFILES


# ═══════════════════════════════════════════════════════════════════════════
#  TestMechanismRegistry
# ═══════════════════════════════════════════════════════════════════════════


class TestMechanismRegistry:
    """Tests covering the MECHANISM_REGISTRY dictionary itself."""

    def test_registry_has_all_manual_entries(self):
        """MECHANISM_REGISTRY must contain at least 22 manually authored entries."""
        assert len(MECHANISM_REGISTRY) >= 22

    def test_registry_keys_match_profiles(self):
        """Every key in MECHANISM_REGISTRY must exist in CROSSLINKERS or REAGENT_PROFILES."""
        all_known_keys = set(CROSSLINKERS) | set(REAGENT_PROFILES)
        for key in MECHANISM_REGISTRY:
            assert key in all_known_keys, (
                f"Registry key '{key}' not found in CROSSLINKERS or REAGENT_PROFILES"
            )

    def test_mechanism_descriptor_display_name(self):
        """Genipin entry must have a non-empty display_name."""
        descriptor = MECHANISM_REGISTRY["genipin"]
        assert descriptor.display_name != ""

    def test_mechanism_descriptor_overall_equation_latex(self):
        """Genipin entry must have a non-empty overall_equation_latex."""
        descriptor = MECHANISM_REGISTRY["genipin"]
        assert descriptor.overall_equation_latex != ""

    def test_mechanism_descriptor_mechanism_type(self):
        """Genipin mechanism_type must be 'nucleophilic_addition'."""
        descriptor = MECHANISM_REGISTRY["genipin"]
        assert descriptor.mechanism_type == "nucleophilic_addition"

    def test_mechanism_descriptor_steps_count(self):
        """Genipin entry must have exactly 3 reaction steps."""
        descriptor = MECHANISM_REGISTRY["genipin"]
        assert len(descriptor.steps) == 3

    def test_mechanism_descriptor_reversibility(self):
        """Genipin reversibility must be 'irreversible'."""
        descriptor = MECHANISM_REGISTRY["genipin"]
        assert descriptor.reversibility == "irreversible"


# ═══════════════════════════════════════════════════════════════════════════
#  TestReactionStep
# ═══════════════════════════════════════════════════════════════════════════


class TestReactionStep:
    """Tests covering ReactionStep fields for the genipin entry."""

    def _genipin_step1(self) -> ReactionStep:
        return MECHANISM_REGISTRY["genipin"].steps[0]

    def test_step_number(self):
        """Genipin step 1 must have step_number == 1."""
        assert self._genipin_step1().step_number == 1

    def test_step_description_not_empty(self):
        """Genipin step 1 description must be non-empty."""
        assert self._genipin_step1().description != ""

    def test_step_bond_formed_not_empty(self):
        """Genipin step 1 bond_formed must be non-empty."""
        assert self._genipin_step1().bond_formed != ""

    def test_step_equation_latex_not_empty(self):
        """Genipin step 1 equation_latex must be non-empty."""
        assert self._genipin_step1().equation_latex != ""


# ═══════════════════════════════════════════════════════════════════════════
#  TestGetMechanism
# ═══════════════════════════════════════════════════════════════════════════


class TestGetMechanism:
    """Tests covering the get_mechanism() function."""

    def test_get_mechanism_known_key(self):
        """get_mechanism('genipin') must return MechanismDescriptor with reagent_key == 'genipin'."""
        result = get_mechanism("genipin")
        assert isinstance(result, MechanismDescriptor)
        assert result.reagent_key == "genipin"

    def test_get_mechanism_fallback_crosslinker(self):
        """For a CROSSLINKER key absent from MECHANISM_REGISTRY, fallback must be valid."""
        # Find a key in CROSSLINKERS that is NOT in MECHANISM_REGISTRY.
        fallback_key = next(
            (k for k in CROSSLINKERS if k not in MECHANISM_REGISTRY),
            None,
        )
        if fallback_key is None:
            # All CROSSLINKER keys happen to be in the registry; test ultimate fallback instead.
            fallback_key = "totally_unknown_reagent_for_fallback_test_xyz"
        result = get_mechanism(fallback_key)
        assert isinstance(result, MechanismDescriptor)
        assert result.reagent_key == fallback_key

    def test_get_mechanism_fallback_reagent_profile(self):
        """get_mechanism('wash_buffer') should return a valid fallback MechanismDescriptor."""
        # wash_buffer is in REAGENT_PROFILES but not in MECHANISM_REGISTRY.
        assert "wash_buffer" not in MECHANISM_REGISTRY
        result = get_mechanism("wash_buffer")
        assert isinstance(result, MechanismDescriptor)
        assert result.display_name != ""

    def test_get_mechanism_fallback_reagent_profile_reagent_key(self):
        """Fallback MechanismDescriptor for wash_buffer must carry the correct reagent_key."""
        result = get_mechanism("wash_buffer")
        assert result.reagent_key == "wash_buffer"

    def test_get_mechanism_unknown_key_does_not_raise(self):
        """get_mechanism with an entirely unknown key must not raise."""
        result = get_mechanism("nonexistent_key_xyz")
        assert isinstance(result, MechanismDescriptor)

    def test_get_mechanism_unknown_key_returns_fallback(self):
        """get_mechanism with an unknown key must return a MechanismDescriptor (ultimate fallback)."""
        result = get_mechanism("nonexistent_key_xyz")
        # reagent_key should echo the input
        assert result.reagent_key == "nonexistent_key_xyz"

    def test_all_crosslinker_keys_have_mechanism(self):
        """Every key in CROSSLINKERS must resolve via get_mechanism without raising."""
        for key in CROSSLINKERS:
            result = get_mechanism(key)
            assert isinstance(result, MechanismDescriptor), (
                f"get_mechanism('{key}') did not return MechanismDescriptor"
            )

    def test_all_reagent_profile_keys_have_mechanism(self):
        """Every key in REAGENT_PROFILES must resolve via get_mechanism without raising."""
        for key in REAGENT_PROFILES:
            result = get_mechanism(key)
            assert isinstance(result, MechanismDescriptor), (
                f"get_mechanism('{key}') did not return MechanismDescriptor"
            )
