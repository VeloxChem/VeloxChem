"""Projected Taylor models and local periodic Shepard factors."""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache

import numpy as np

from .phase import principal_periodic_delta


@lru_cache(maxsize=32)
def _torsion_masks(row_types):
    """Cache the row-kind masks; a chart is fixed for a model's lifetime."""

    kinds = np.asarray(row_types)
    return kinds == "proper_torsion", kinds == "improper_torsion"


@lru_cache(maxsize=32)
def _is_identity_projector(projector_rows, n_rows):
    return set(projector_rows) == set(range(n_rows))


def taylor_displacement(current, reference, row_types):
    """Return the persisted smooth displacement and its first derivative."""

    current = np.asarray(current, dtype=float)
    reference = np.asarray(reference, dtype=float)
    raw = current - reference
    displacement = raw.copy()
    chain = np.ones_like(raw)
    proper, improper = _torsion_masks(tuple(row_types))
    if proper.any():
        delta = principal_periodic_delta(raw[proper])
        # sin(delta) and its derivative are continuous across the
        # principal-angle branch at +/- pi.  The earlier half-angle chart
        # changed sign there and produced finite energy jumps in dynamics.
        displacement[proper] = np.sin(delta)
        chain[proper] = np.cos(delta)
    if improper.any():
        delta = principal_periodic_delta(raw[improper])
        cosine = np.maximum(np.cos(0.5 * delta), 1.0e-6)
        displacement[improper] = 2.0 * np.tan(0.5 * delta)
        chain[improper] = 1.0 / (cosine * cosine)
    return displacement, chain


def evaluate_projected_taylor(state, current_coordinates, row_types, projector_rows):
    """Evaluate a Taylor state and its internal gradient on explicit rows."""

    reference = np.asarray(state["internal_coordinates"], dtype=float)
    gradient = np.asarray(state["internal_gradient"], dtype=float)
    hessian = np.asarray(state["internal_hessian"], dtype=float)
    displacement, chain = taylor_displacement(
        current_coordinates, reference, row_types
    )
    if _is_identity_projector(tuple(projector_rows), displacement.size):
        # A projector that covers the whole chart is the identity.  Skipping
        # the fancy-index copies here matters: this is the inner loop of both
        # the additive assembly and its finite-difference derivative.
        projected_displacement = displacement
        projected_gradient = gradient
        h_d = hessian @ displacement
    else:
        rows = np.asarray(projector_rows, dtype=int)
        projected_displacement = np.zeros_like(displacement)
        projected_displacement[rows] = displacement[rows]
        projected_gradient = np.zeros_like(gradient)
        projected_gradient[rows] = gradient[rows]
        projected_hessian = np.zeros_like(hessian)
        projected_hessian[np.ix_(rows, rows)] = hessian[np.ix_(rows, rows)]
        h_d = projected_hessian @ projected_displacement
    energy = (
        float(state["energy"])
        + float(np.dot(projected_gradient, projected_displacement))
        + 0.5 * float(np.dot(projected_displacement, h_d))
    )
    internal_gradient = chain * (projected_gradient + h_d)
    return energy, internal_gradient


def periodic_signature_distance_and_gradient(
    current,
    reference,
    signature_rows,
    row_orientations,
    scales,
):
    """Periodic mean-square signature metric and its q-gradient."""

    current = np.asarray(current, dtype=float)
    reference = np.asarray(reference, dtype=float)
    rows = np.asarray(signature_rows, dtype=int)
    signs = np.asarray(row_orientations, dtype=float)
    scales = np.asarray(scales, dtype=float)
    if len(rows) == 0:
        raise ValueError("A factor requires at least one signature row.")
    delta = principal_periodic_delta(current[rows] - reference[rows]) * signs
    prefactor = 1.0 / (float(len(rows)) * np.square(scales))
    distance_squared = float(np.sum(2.0 * (1.0 - np.cos(delta)) * prefactor))
    gradient = np.zeros_like(current)
    gradient[rows] = 2.0 * np.sin(delta) * signs * prefactor
    return distance_squared, gradient


def local_shepard_weights(
    current,
    states,
    *,
    signature_rows,
    row_orientations,
    scales,
    confidence_radius,
    exponent_p,
    exponent_q,
    exact_center_tolerance,
):
    """Return normalized local weights and analytical q-gradients."""

    distances = []
    distance_gradients = []
    for state in states:
        d2, grad_d2 = periodic_signature_distance_and_gradient(
            current,
            state["internal_coordinates"],
            signature_rows,
            row_orientations,
            scales,
        )
        distances.append(d2)
        distance_gradients.append(grad_d2)
    distances = np.asarray(distances, dtype=float)
    distance_gradients = np.asarray(distance_gradients, dtype=float)
    centers = distances <= float(exact_center_tolerance) ** 2
    if np.any(centers):
        weights = centers.astype(float) / float(np.count_nonzero(centers))
        return weights, np.zeros_like(distance_gradients), distances

    radius2 = float(confidence_radius) ** 2
    if radius2 <= 0.0:
        raise ValueError("A factor confidence radius must be positive.")
    u = distances / radius2
    p = float(exponent_p)
    q = float(exponent_q)
    denominator = np.power(u, p) + np.power(u, q)
    raw = 1.0 / denominator
    # d(raw)/d(q) = -D'(u) / D(u)^2 * d(d2)/d(q) / rho^2
    denominator_prime = p * np.power(u, p - 1.0) + q * np.power(u, q - 1.0)
    raw_gradients = (
        -denominator_prime[:, None]
        / np.square(denominator)[:, None]
        * distance_gradients
        / radius2
    )
    scale = float(np.max(raw))
    if not np.isfinite(scale) or scale <= 0.0:
        nearest = int(np.argmin(distances))
        weights = np.zeros(len(states), dtype=float)
        weights[nearest] = 1.0
        return weights, np.zeros_like(distance_gradients), distances
    raw /= scale
    raw_gradients /= scale
    total = float(np.sum(raw))
    total_gradient = np.sum(raw_gradients, axis=0)
    weights = raw / total
    weight_gradients = (
        raw_gradients * total - raw[:, None] * total_gradient[None, :]
    ) / (total * total)
    return weights, weight_gradients, distances


def pair_shepard_weights(
    current,
    states,
    rotors,
    *,
    confidence_radius,
    exponent_p,
    exponent_q,
    exact_center_tolerance,
):
    """Two-motion Shepard weights with ``d_pair^2 = d_a^2 + d_b^2``."""

    distances = []
    distance_gradients = []
    for state in states:
        distance = 0.0
        gradient = np.zeros_like(np.asarray(current, dtype=float))
        for rotor in rotors:
            value, derivative = periodic_signature_distance_and_gradient(
                current,
                state["internal_coordinates"],
                rotor["signature_rows"],
                rotor["row_orientations"],
                rotor["signature_row_scales"],
            )
            distance += value
            gradient += derivative
        distances.append(distance)
        distance_gradients.append(gradient)
    distances = np.asarray(distances, dtype=float)
    distance_gradients = np.asarray(distance_gradients, dtype=float)
    centers = distances <= float(exact_center_tolerance) ** 2
    if np.any(centers):
        weights = centers.astype(float) / float(np.count_nonzero(centers))
        return weights, np.zeros_like(distance_gradients), distances
    radius2 = float(confidence_radius) ** 2
    if radius2 <= 0.0:
        raise ValueError("An interaction confidence radius must be positive.")
    u = distances / radius2
    p = float(exponent_p)
    q = float(exponent_q)
    denominator = np.power(u, p) + np.power(u, q)
    raw = 1.0 / denominator
    denominator_prime = p * np.power(u, p - 1.0) + q * np.power(u, q - 1.0)
    raw_gradients = (
        -denominator_prime[:, None]
        / np.square(denominator)[:, None]
        * distance_gradients
        / radius2
    )
    scale = float(np.max(raw))
    if not np.isfinite(scale) or scale <= 0.0:
        nearest = int(np.argmin(distances))
        weights = np.zeros(len(states), dtype=float)
        weights[nearest] = 1.0
        return weights, np.zeros_like(distance_gradients), distances
    raw /= scale
    raw_gradients /= scale
    total = float(np.sum(raw))
    total_gradient = np.sum(raw_gradients, axis=0)
    weights = raw / total
    weight_gradients = (
        raw_gradients * total - raw[:, None] * total_gradient[None, :]
    ) / (total * total)
    return weights, weight_gradients, distances


@dataclass
class LocalShepardFactor:
    subset: dict
    rotor: dict
    states: tuple[dict, ...]
    anchor_state: dict
    policy: dict
    row_types: tuple[str, ...]

    def evaluate(self, current_coordinates, *, include_virtual=True):
        states = self.states
        if not include_virtual:
            states = tuple(
                state for state in states if not state.get("is_virtual_image", False)
            )
        weights, weight_gradients, distances = local_shepard_weights(
            current_coordinates,
            states,
            signature_rows=self.rotor["signature_rows"],
            row_orientations=self.rotor["row_orientations"],
            scales=self.rotor["signature_row_scales"],
            confidence_radius=self.policy["confidence_radius"],
            exponent_p=self.policy["shepard_p"],
            exponent_q=self.policy["shepard_q"],
            exact_center_tolerance=self.policy["exact_center_tolerance"],
        )
        energies = []
        gradients = []
        projector_rows = self.subset["projector_rows"]
        for state in states:
            energy, gradient = evaluate_projected_taylor(
                state, current_coordinates, self.row_types, projector_rows
            )
            energies.append(energy)
            gradients.append(gradient)
        energies = np.asarray(energies, dtype=float)
        gradients = np.asarray(gradients, dtype=float)
        factor_energy = float(np.dot(weights, energies))
        factor_gradient = (
            np.tensordot(weights, gradients, axes=1)
            + np.tensordot(energies - factor_energy, weight_gradients, axes=1)
        )
        baseline_energy, baseline_gradient = evaluate_projected_taylor(
            self.anchor_state,
            current_coordinates,
            self.row_types,
            projector_rows,
        )
        return {
            "energy": factor_energy,
            "gradient": factor_gradient,
            "baseline_energy": baseline_energy,
            "baseline_gradient": baseline_gradient,
            "residual_energy": factor_energy - baseline_energy,
            "residual_gradient": factor_gradient - baseline_gradient,
            "weights": weights,
            "distances": distances,
        }


@dataclass
class PairInteractionFactor:
    """Interpolate an explicit inclusion-exclusion interaction residual."""

    subset: dict
    rotors: tuple[dict, dict]
    states: tuple[dict, ...]
    policy: dict
    row_types: tuple[str, ...]

    def evaluate(self, current_coordinates):
        weights, weight_gradients, distances = pair_shepard_weights(
            current_coordinates,
            self.states,
            self.rotors,
            confidence_radius=self.policy["confidence_radius"],
            exponent_p=self.policy["shepard_p"],
            exponent_q=self.policy["shepard_q"],
            exact_center_tolerance=self.policy["exact_center_tolerance"],
        )
        energies = []
        gradients = []
        for state in self.states:
            energy, gradient = evaluate_projected_taylor(
                state,
                current_coordinates,
                self.row_types,
                self.subset["projector_rows"],
            )
            energies.append(energy)
            gradients.append(gradient)
        energies = np.asarray(energies, dtype=float)
        gradients = np.asarray(gradients, dtype=float)
        energy = float(np.dot(weights, energies))
        gradient = (
            np.tensordot(weights, gradients, axes=1)
            + np.tensordot(energies - energy, weight_gradients, axes=1)
        )
        return {
            "residual_energy": energy,
            "residual_gradient": gradient,
            "weights": weights,
            "distances": distances,
        }
