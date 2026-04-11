"""Steric Mass Action (SMA) isotherm for ion exchange chromatography.

Architecture: docs/13_module2_module3_final_implementation_plan.md, Phase E.

The SMA model (Brooks & Cramer, 1992) accounts for the steric shielding
effect of bound proteins on accessible ionic capacity:

    Lambda = sum_i( (z_i + sigma_i) * q_i ) + q_salt

At equilibrium for each protein component i:

    q_i = K_eq_i * C_i * (q_salt / C_salt)^z_i

where q_salt is the stationary-phase counter-ion concentration, found
by solving the ionic capacity constraint (fixed-point iteration).

References:
    Brooks, C.A. & Cramer, S.M. (1992). Steric mass-action ion exchange:
    Displacement profiles and induced salt gradients. AIChE J., 38(12), 1969.
"""

from __future__ import annotations

import logging
import warnings
from dataclasses import dataclass, field

import numpy as np

logger = logging.getLogger(__name__)

_MAX_ITER = 200
_TOL = 1e-10


@dataclass
class SMAIsotherm:
    """Steric Mass Action isotherm for ion exchange.

    Attributes:
        Lambda: Total ionic capacity of the resin [mol/m^3 solid].
        z: Characteristic charge per component [-], shape (n_comp,).
        sigma: Steric factor per component [-], shape (n_comp,).
        K_eq: Equilibrium constant per component [dimensionless], shape (n_comp,).
    """

    Lambda: float = 1000.0
    z: np.ndarray = field(default_factory=lambda: np.array([3.0, 5.0]))
    sigma: np.ndarray = field(default_factory=lambda: np.array([50.0, 80.0]))
    K_eq: np.ndarray = field(default_factory=lambda: np.array([1e-3, 1e-5]))

    def __post_init__(self) -> None:
        self.z = np.asarray(self.z, dtype=float)
        self.sigma = np.asarray(self.sigma, dtype=float)
        self.K_eq = np.asarray(self.K_eq, dtype=float)

        n = len(self.z)
        if not (len(self.sigma) == n == len(self.K_eq)):
            raise ValueError("z, sigma, K_eq must all have the same length.")

    @property
    def n_components(self) -> int:
        """Number of protein components (excludes salt)."""
        return len(self.z)

    def _q_salt_from_loading(self, q_proteins: np.ndarray) -> float:
        """Compute q_salt from the ionic capacity constraint.

        Lambda = sum_i( (z_i + sigma_i) * q_i ) + q_salt

        Args:
            q_proteins: Bound protein concentrations [mol/m^3 solid], shape (n_comp,).

        Returns:
            q_salt [mol/m^3 solid].
        """
        occupied = np.dot(self.z + self.sigma, q_proteins)
        return max(self.Lambda - occupied, 0.0)

    def equilibrium_loading(
        self,
        C: np.ndarray,
        salt_concentration: float,
    ) -> np.ndarray:
        """Compute equilibrium loading via fixed-point iteration.

        SMA equilibrium:
            q_i = K_eq_i * C_i * (q_salt / C_salt)^z_i
            q_salt = Lambda - sum_i( (z_i + sigma_i) * q_i )

        The two equations are coupled: q_i depends on q_salt which
        depends on q_i.  Solved by fixed-point iteration on q_salt.

        Args:
            C: Protein concentrations in mobile phase [mol/m^3], shape (n_comp,).
            salt_concentration: Salt concentration in mobile phase [mol/m^3].

        Returns:
            Equilibrium loadings [mol/m^3 solid], shape (n_comp,).
        """
        C = np.asarray(C, dtype=float)
        C_safe = np.maximum(C, 0.0)
        C_salt = max(salt_concentration, 1e-6)  # avoid division by zero

        # Initial guess: q_salt = Lambda (no protein bound yet)
        q_salt = float(self.Lambda)

        for iteration in range(_MAX_ITER):
            # Compute protein loadings from current q_salt
            ratio = q_salt / C_salt
            # Avoid numerical issues with very small ratio
            ratio = max(ratio, 1e-20)
            q_proteins = self.K_eq * C_safe * (ratio ** self.z)

            # Update q_salt from capacity constraint
            q_salt_new = self._q_salt_from_loading(q_proteins)
            q_salt_new = max(q_salt_new, 0.0)

            if abs(q_salt_new - q_salt) < _TOL * (1.0 + abs(q_salt)):
                q_salt = q_salt_new
                break

            # Damped update to improve convergence
            q_salt = 0.5 * q_salt + 0.5 * q_salt_new
        else:
            warnings.warn(
                f"SMA fixed-point iteration did not converge in {_MAX_ITER} "
                f"iterations (C_salt={C_salt:.3g}, C={C_safe}). "
                "Using last iterate.",
                stacklevel=2,
            )

        # Final protein loadings
        ratio = max(q_salt / C_salt, 1e-20)
        q_proteins = self.K_eq * C_safe * (ratio ** self.z)

        # Clip to physically feasible range
        q_proteins = np.maximum(q_proteins, 0.0)

        logger.debug(
            "SMA converged in %d iterations: q_salt=%.3g, q=%s",
            iteration + 1, q_salt, q_proteins,
        )

        return q_proteins

    def retention_factor(
        self,
        C_salt: float,
        component_idx: int = 0,
        C_protein: float = 1e-6,
    ) -> float:
        """Estimate retention factor k' for a single component.

        In SMA at low protein concentration:
            k' ~ K_eq_i * (Lambda / C_salt)^z_i * (1 - eps_b) / eps_b

        This simplified form (ignoring steric shielding at low load) is
        useful for gradient scouting.

        Args:
            C_salt: Salt concentration [mol/m^3].
            component_idx: Component index.
            C_protein: Protein concentration (for self-shielding estimate).

        Returns:
            Retention factor k' [-] (column-void-fraction-independent part).
        """
        C_salt = max(C_salt, 1e-6)
        ratio = self.Lambda / C_salt
        k_prime = self.K_eq[component_idx] * (ratio ** self.z[component_idx])
        return float(k_prime)

    def validate(self) -> list[str]:
        """Check parameter validity.

        Returns:
            List of violation messages (empty = OK).
        """
        errors: list[str] = []
        if self.Lambda <= 0:
            errors.append(f"Lambda must be positive, got {self.Lambda}")
        if np.any(self.z <= 0):
            errors.append(f"All characteristic charges z must be positive, got {self.z}")
        if np.any(self.sigma < 0):
            errors.append(f"Steric factors sigma must be non-negative, got {self.sigma}")
        if np.any(self.K_eq <= 0):
            errors.append(f"All K_eq must be positive, got {self.K_eq}")
        # Check that sigma doesn't prevent any binding
        min_sigma_sum = np.dot(self.z + self.sigma, np.ones(self.n_components))
        if min_sigma_sum >= self.Lambda:
            errors.append(
                f"sigma+z sum ({min_sigma_sum:.1f}) >= Lambda ({self.Lambda:.1f}); "
                "no capacity for binding."
            )
        return errors
