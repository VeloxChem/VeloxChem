"""Numerical validation shared by grouped construction and user evaluation."""

from __future__ import annotations

import numpy as np

from .coordinates import build_interpolation_datapoint
from .structure_generation import clone_molecule_with_coordinates


def evaluate_runtime_at_coordinates(
    runtime,
    molecule_template,
    coordinates_bohr,
    z_matrix,
    settings,
):
    molecule = clone_molecule_with_coordinates(molecule_template, coordinates_bohr)
    point = build_interpolation_datapoint(
        z_matrix,
        coordinates_bohr,
        settings,
        masses=molecule.get_masses(),
        eq_bond_lengths=runtime.eq_bond_lengths,
    )
    return runtime.evaluate_coordinate(point, molecule)


def validate_training_points(runtime, molecule, z_matrix, settings, physical_results):
    records = []
    for result in physical_results:
        energy, gradient = evaluate_runtime_at_coordinates(
            runtime,
            molecule,
            result["cartesian_coordinates"],
            z_matrix,
            settings,
        )
        reference_gradient = np.asarray(result["cartesian_gradient"], dtype=float)
        error = np.asarray(gradient) - reference_gradient
        records.append(
            {
                "request_id": result["request_id"],
                "energy_error_hartree": float(energy - result["energy"]),
                "gradient_rms_error_hartree_per_bohr": float(np.sqrt(np.mean(error * error))),
                "maximum_force_error_hartree_per_bohr": float(np.max(np.abs(error))),
            }
        )
    return records


def finite_difference_runtime_gradient(
    runtime,
    molecule,
    z_matrix,
    settings,
    *,
    step=1.0e-4,
    dof_indices=None,
):
    coordinates = np.asarray(molecule.get_coordinates_in_bohr(), dtype=float)
    _, analytical = evaluate_runtime_at_coordinates(
        runtime, molecule, coordinates, z_matrix, settings
    )
    analytical_flat = analytical.reshape(-1)
    if dof_indices is None:
        dof_indices = tuple(range(min(12, coordinates.size)))
    errors = []
    for dof in dof_indices:
        atom, component = divmod(int(dof), 3)
        plus = coordinates.copy()
        minus = coordinates.copy()
        plus[atom, component] += step
        minus[atom, component] -= step
        energy_plus, _ = evaluate_runtime_at_coordinates(
            runtime, molecule, plus, z_matrix, settings
        )
        energy_minus, _ = evaluate_runtime_at_coordinates(
            runtime, molecule, minus, z_matrix, settings
        )
        finite_difference = (energy_plus - energy_minus) / (2.0 * step)
        errors.append(
            {
                "dof": int(dof),
                "analytical": float(analytical_flat[dof]),
                "finite_difference": float(finite_difference),
                "absolute_error": float(abs(finite_difference - analytical_flat[dof])),
            }
        )
    return {
        "step_bohr": float(step),
        "components": errors,
        "maximum_absolute_error": max(
            (entry["absolute_error"] for entry in errors), default=0.0
        ),
    }
