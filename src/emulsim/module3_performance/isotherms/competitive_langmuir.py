"""Competitive Langmuir isotherm for multi-component adsorption.

Architecture: docs/13_module2_module3_final_implementation_plan.md, Phase E.

Competitive Langmuir assumes all components compete for the same
adsorption sites with a shared denominator:

    q_i* = q_max_i * K_i * C_i / (1 + sum_j(K_j * C_j))

This reduces to single-component Langmuir when only one species is present.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np


@dataclass
class CompetitiveLangmuirIsotherm:
    """Multi-component competitive Langmuir isotherm.

    All components compete for the same adsorption sites.

    Attributes:
        q_max: Maximum binding capacity per component [mol/m^3 solid], shape (n_comp,).
        K_L: Langmuir equilibrium constant per component [m^3/mol], shape (n_comp,).
    """

    q_max: np.ndarray = field(default_factory=lambda: np.array([100.0, 80.0]))
    K_L: np.ndarray = field(default_factory=lambda: np.array([1e3, 5e2]))

    def __post_init__(self) -> None:
        self.q_max = np.asarray(self.q_max, dtype=float)
        self.K_L = np.asarray(self.K_L, dtype=float)
        if self.q_max.shape != self.K_L.shape:
            raise ValueError(
                f"q_max shape {self.q_max.shape} != K_L shape {self.K_L.shape}"
            )

    @property
    def n_components(self) -> int:
        """Number of components."""
        return len(self.q_max)

    def equilibrium_loading(self, C: np.ndarray) -> np.ndarray:
        """Compute equilibrium loading for all components.

        q_i = q_max_i * K_i * C_i / (1 + sum_j(K_j * C_j))

        Args:
            C: Mobile-phase concentrations [mol/m^3], shape (n_comp,) or (n_comp, N).
               For a grid of N spatial cells, shape is (n_comp, N).

        Returns:
            Equilibrium loading [mol/m^3 solid], same shape as C.
        """
        C = np.asarray(C, dtype=float)
        scalar_input = C.ndim == 1

        if scalar_input:
            # shape (n_comp,) -> treat as single point
            C_safe = np.maximum(C, 0.0)
            denominator = 1.0 + np.dot(self.K_L, C_safe)
            q = self.q_max * self.K_L * C_safe / denominator
        else:
            # shape (n_comp, N)
            C_safe = np.maximum(C, 0.0)
            # sum_j(K_j * C_j): shape (N,)
            denominator = 1.0 + np.einsum("i,ij->j", self.K_L, C_safe)
            # q_i = q_max_i * K_i * C_i / denom: shape (n_comp, N)
            q = (self.q_max[:, np.newaxis] * self.K_L[:, np.newaxis] * C_safe) / denominator[np.newaxis, :]

        return q

    def jacobian(self, C: np.ndarray) -> np.ndarray:
        """Compute dq_i/dC_j matrix (diagonal dominant for competitive systems).

        dq_i/dC_i = q_max_i * K_i * (1 + sum_{j!=i} K_j*C_j) / D^2
        dq_i/dC_j = -q_max_i * K_i * C_i * K_j / D^2   (j != i)

        where D = 1 + sum_k(K_k * C_k)

        Args:
            C: Mobile-phase concentrations [mol/m^3], shape (n_comp,).

        Returns:
            Jacobian matrix [m^3/m^3], shape (n_comp, n_comp).
        """
        C = np.asarray(C, dtype=float)
        C_safe = np.maximum(C, 0.0)
        n = self.n_components
        D = 1.0 + np.dot(self.K_L, C_safe)
        D2 = D ** 2

        J = np.zeros((n, n))
        for i in range(n):
            for j in range(n):
                if i == j:
                    # dq_i/dC_i
                    sum_other = np.dot(self.K_L, C_safe) - self.K_L[i] * C_safe[i]
                    J[i, j] = self.q_max[i] * self.K_L[i] * (1.0 + sum_other) / D2
                else:
                    # dq_i/dC_j
                    J[i, j] = -self.q_max[i] * self.K_L[i] * C_safe[i] * self.K_L[j] / D2

        return J

    def validate(self) -> list[str]:
        """Check parameter validity.

        Returns:
            List of violation messages (empty = OK).
        """
        errors: list[str] = []
        if np.any(self.q_max <= 0):
            errors.append(f"All q_max must be positive, got {self.q_max}")
        if np.any(self.K_L <= 0):
            errors.append(f"All K_L must be positive, got {self.K_L}")
        if self.n_components < 2:
            errors.append(
                "CompetitiveLangmuirIsotherm requires at least 2 components; "
                "use LangmuirIsotherm for single-component."
            )
        return errors
