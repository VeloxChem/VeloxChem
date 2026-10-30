# SPDX-License-Identifier: BSD-3-Clause
"""Small Cartesian potential result shared by evaluation and dynamics."""

from dataclasses import dataclass
import numpy as np


@dataclass(frozen=True)
class PotentialResult:
    """An owned energy (Hartree) and Cartesian gradient (Hartree/bohr).

    Atom-major gradients have shape (N, 3) and are never mass weighted.
    Results do not alias calculator storage or subsequent evaluations.
    """

    energy: float
    gradient: np.ndarray
    model_id: str

    def __post_init__(self):
        energy = float(self.energy)
        gradient = np.array(self.gradient, dtype=np.float64, copy=True)
        if gradient.ndim != 2 or gradient.shape[1] != 3 or not len(gradient):
            raise ValueError('Potential gradient must have shape (N, 3), N > 0.')
        if not np.isfinite(energy) or not np.all(np.isfinite(gradient)):
            raise ValueError('Potential energy and gradient must be finite.')
        gradient.setflags(write=False)
        object.__setattr__(self, 'energy', energy)
        object.__setattr__(self, 'gradient', gradient)

    @property
    def forces(self):
        """Physical Cartesian forces in Hartree/bohr."""
        return -self.gradient
