"""Shared test fixtures for emulsim tests."""

import pytest
import numpy as np

from emulsim.datatypes import (
    MaterialProperties,
    MixerGeometry,
    SimulationParameters,
)


@pytest.fixture
def default_params():
    """Default simulation parameters."""
    return SimulationParameters()


@pytest.fixture
def default_props():
    """Default material properties."""
    return MaterialProperties()


@pytest.fixture
def default_mixer():
    """Default mixer geometry."""
    return MixerGeometry()


@pytest.fixture
def rng():
    """Seeded random number generator for reproducible tests."""
    return np.random.default_rng(42)
