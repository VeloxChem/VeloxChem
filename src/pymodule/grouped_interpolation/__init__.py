"""Symmetry-aware grouped Modified Shepard interpolation.

This package is intentionally independent of the historical symmetry helpers.
It exposes small construction/runtime components and keeps the legacy flat MSI
path opt-in compatible through :class:`veloxchem.InterpolationDriver`.
"""

from .models import (
    CoordinateDefinition,
    GroupedCalculationRequest,
    GroupedCalculationResult,
    LocalSubset,
    PrimitiveRotor,
    SymmetryOperation,
)
from .motion_detection import PrimitiveRotorDetector
from .phase import estimate_collective_phase, principal_periodic_delta
from .structure_generation import RigidRotorStructureGenerator
from .runtime import GroupedRuntimeModel

__all__ = [
    "CoordinateDefinition",
    "GroupedCalculationRequest",
    "GroupedCalculationResult",
    "GroupedRuntimeModel",
    "LocalSubset",
    "PrimitiveRotor",
    "PrimitiveRotorDetector",
    "RigidRotorStructureGenerator",
    "SymmetryOperation",
    "estimate_collective_phase",
    "principal_periodic_delta",
]
