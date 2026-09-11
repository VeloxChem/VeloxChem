"""Periodic collective phase handling for redundant rotor torsions."""

from __future__ import annotations

import numpy as np

from .models import PrimitiveRotor


def principal_periodic_delta(delta):
    """Return angle differences in the half-open interval [-pi, pi)."""

    arr = np.asarray(delta, dtype=float)
    wrapped = (arr + np.pi) % (2.0 * np.pi) - np.pi
    if np.ndim(delta) == 0:
        return float(wrapped)
    return wrapped


def estimate_collective_phase(
    current_internal_coordinates,
    reference_internal_coordinates,
    rotor: PrimitiveRotor,
) -> float:
    """Estimate one collective rotor phase without multiplying redundant rows."""

    current = np.asarray(current_internal_coordinates, dtype=float)
    reference = np.asarray(reference_internal_coordinates, dtype=float)
    rows = np.asarray(rotor.signature_rows, dtype=int)
    signs = np.asarray(rotor.row_orientations, dtype=float)
    weights = 1.0 / np.square(np.asarray(rotor.signature_row_scales, dtype=float))
    deltas = principal_periodic_delta(current[rows] - reference[rows]) * signs
    sine = float(np.dot(weights, np.sin(deltas)))
    cosine = float(np.dot(weights, np.cos(deltas)))
    if abs(sine) < 1.0e-15 and abs(cosine) < 1.0e-15:
        return 0.0
    return float(np.arctan2(sine, cosine))
