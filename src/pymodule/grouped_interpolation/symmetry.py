"""Discovery, validation, and application of methyl C3 operations."""

from __future__ import annotations

from dataclasses import replace
import hashlib
import itertools
import json

import numpy as np

from .models import GroupedCalculationResult, SymmetryOperation
from .phase import principal_periodic_delta
from .structure_generation import rotate_coordinates_about_axis


def _equivalent_row(candidate, row):
    candidate = tuple(candidate)
    row = tuple(row)
    if len(row) == 2:
        return candidate == row or candidate == row[::-1]
    if len(row) == 3:
        return candidate == row or candidate == row[::-1]
    if len(row) == 4:
        return candidate == row or candidate == row[::-1]
    return candidate == row


def _row_scales(coordinate_definition, row_permutation):
    """Chart-induced linear factor relating a row to its permutation partner.

    Relabelling equivalent atoms maps the *distance* of one row onto the
    distance of another, but the chart value is not always the distance itself.
    An equilibrium-scaled bond row stores ``r_eq - r_eq**2 / r`` with its own
    reference length, so two exchanged rows are related by a scale as well as an
    offset.  Inverse-distance and plain-distance charts, angles and torsions all
    have unit scale.
    """

    scales = np.ones(len(coordinate_definition.row_types), dtype=float)
    if coordinate_definition.bond_convention != "equilibrium_scaled_distance":
        return scales
    lengths = np.asarray(
        coordinate_definition.equilibrium_bond_lengths_bohr, dtype=float
    )
    bond_indices = [
        index
        for index, kind in enumerate(coordinate_definition.row_types)
        if kind == "bond"
    ]
    reference = {index: lengths[position] for position, index in enumerate(bond_indices)}
    for index in bond_indices:
        source = int(row_permutation[index])
        scales[index] = (reference[index] / reference[source]) ** 2
    return scales


def _discover_row_map(rows, atom_permutation, source_q, transformed_q, row_types):
    row_permutation = []
    row_signs = []
    for target_index, target_row in enumerate(rows):
        source_atoms = tuple(atom_permutation[atom] for atom in target_row)
        matches = [
            index for index, source_row in enumerate(rows)
            if len(source_row) == len(source_atoms) and _equivalent_row(source_atoms, source_row)
        ]
        if len(matches) != 1:
            raise ValueError(
                f"Coordinate chart is not closed under symmetry for row {target_row}; "
                f"found {len(matches)} matches."
            )
        source_index = matches[0]
        # Full reversal preserves the standard distance/angle/dihedral value.
        row_permutation.append(source_index)
        row_signs.append(1)
    return tuple(row_permutation), tuple(row_signs)


def _discover_row_offsets(
    row_permutation, row_signs, row_scales, source_q, transformed_q, row_types
):
    offsets = []
    for target_index, source_index in enumerate(row_permutation):
        raw_offset = (
            transformed_q[target_index]
            - row_signs[target_index] * row_scales[target_index] * source_q[source_index]
        )
        if row_types[target_index] in {"proper_torsion", "improper_torsion"}:
            raw_offset = principal_periodic_delta(raw_offset)
        offsets.append(float(raw_offset))
    return tuple(offsets)


def discover_cyclic_operations(
    rotor,
    coordinate_definition,
    anchor_coordinates,
    coordinate_value_function,
    atom_labels=None,
):
    """Build exact cyclic atom maps and measure rigid-rotation coverage."""

    if rotor.symmetry_order <= 1:
        return ()
    if len(rotor.equivalent_atom_groups) != 1:
        raise ValueError(
            f"Rotor {rotor.rotor_id} needs one explicit equivalent-atom group."
        )

    source_xyz = np.asarray(anchor_coordinates, dtype=float)
    source_q = np.asarray(coordinate_value_function(source_xyz), dtype=float)
    owned = tuple(rotor.equivalent_atom_groups[0])
    labels = None if atom_labels is None else tuple(str(value) for value in atom_labels)
    operations = []
    for step in (1, 2):
        phase_offset = step * 2.0 * np.pi / rotor.symmetry_order
        rotated = rotate_coordinates_about_axis(
            source_xyz,
            axis_atoms=rotor.axis_atoms,
            moving_atoms=rotor.moving_side_atoms,
            angle=phase_offset,
        )
        best = None
        for candidate in itertools.permutations(owned):
            permutation = np.arange(len(source_xyz), dtype=int)
            for target_atom, source_atom in zip(owned, candidate):
                permutation[target_atom] = source_atom
            transformed = source_xyz[permutation]
            errors = np.linalg.norm(transformed[list(owned)] - rotated[list(owned)], axis=1)
            maximum = float(np.max(errors))
            if best is None or maximum < best[0]:
                best = (maximum, permutation, transformed)
        geometry_error, atom_permutation, transformed_xyz = best
        exact_permutation = True
        if labels is not None:
            exact_permutation = all(
                labels[target] == labels[source]
                for target, source in enumerate(atom_permutation)
            )
        if not exact_permutation:
            raise ValueError(
                f"Symmetry operation for {rotor.rotor_id} exchanges unlike atoms."
            )
        transformed_q = np.asarray(
            coordinate_value_function(transformed_xyz), dtype=float
        )
        row_permutation, row_signs = _discover_row_map(
            coordinate_definition.rows,
            atom_permutation,
            source_q,
            transformed_q,
            coordinate_definition.row_types,
        )
        row_scales = _row_scales(coordinate_definition, row_permutation)
        row_offsets = _discover_row_offsets(
            row_permutation,
            row_signs,
            row_scales,
            source_q,
            transformed_q,
            coordinate_definition.row_types,
        )
        row_scales = tuple(float(value) for value in row_scales)
        payload = {
            "rotor_id": rotor.rotor_id,
            "step": step,
            "atom_permutation": atom_permutation.tolist(),
            "row_permutation": row_permutation,
            "row_signs": row_signs,
            "row_offsets": row_offsets,
            "row_scales": row_scales,
            "exact_atom_permutation": exact_permutation,
        }
        validation_hash = hashlib.sha256(
            json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
        ).hexdigest()
        operations.append(
            SymmetryOperation(
                operation_id=f"{rotor.rotor_id}.c3+{step}",
                source_rotor_id=rotor.rotor_id,
                symmetry_order=rotor.symmetry_order,
                atom_permutation=tuple(int(value) for value in atom_permutation),
                row_permutation=row_permutation,
                row_signs=row_signs,
                row_offsets=row_offsets,
                row_scales=row_scales,
                phase_offset=float(principal_periodic_delta(phase_offset)),
                # Element-preserving atom permutation plus chart closure is an
                # exact symmetry operation.  The geometry error only measures
                # how closely it covers an independently rigid-rotated state.
                geometry_validated=True,
                derivative_validated=exact_permutation,
                max_geometry_error=geometry_error,
                max_energy_error=0.0,
                max_gradient_error=0.0,
                max_hessian_error=0.0,
                validation_hash=validation_hash,
                exact_atom_permutation=exact_permutation,
                coverage_validated=False,
            )
        )
    return tuple(operations)


def _scales_of(operation, n_rows):
    """Row scales of an operation, tolerating models written before they existed."""

    if isinstance(operation, dict):
        values = operation.get("row_scales") or ()
    else:
        values = getattr(operation, "row_scales", ()) or ()
    if len(values) != n_rows:
        return np.ones(n_rows, dtype=float)
    return np.asarray(values, dtype=float)


def _permuted_cartesian_hessian(hessian, atom_permutation):
    dof_permutation = np.array(
        [3 * atom + component for atom in atom_permutation for component in range(3)],
        dtype=int,
    )
    return hessian[np.ix_(dof_permutation, dof_permutation)]


def transform_result(request, source_result, operation, measured_phase_signature):
    """Create a provenance-preserving virtual result without a QM calculation."""

    atom_permutation = np.asarray(operation.atom_permutation, dtype=int)
    row_permutation = np.asarray(operation.row_permutation, dtype=int)
    row_signs = np.asarray(operation.row_signs, dtype=float)
    row_scales = _scales_of(operation, row_permutation.size)
    coordinates = np.asarray(source_result["cartesian_coordinates"], dtype=float)[atom_permutation]
    cart_gradient = np.asarray(source_result["cartesian_gradient"], dtype=float)[atom_permutation]
    internal_coordinates = (
        row_signs
        * row_scales
        * np.asarray(source_result["internal_coordinates"], dtype=float)[row_permutation]
        + np.asarray(operation.row_offsets, dtype=float)
    )
    internal_gradient = (
        row_signs
        * np.asarray(source_result["internal_gradient"], dtype=float)[row_permutation]
        / row_scales
    )
    cart_hessian = source_result.get("cartesian_hessian")
    if cart_hessian is not None:
        cart_hessian = _permuted_cartesian_hessian(
            np.asarray(cart_hessian, dtype=float), atom_permutation
        )
    internal_hessian = source_result.get("internal_hessian")
    if internal_hessian is not None:
        source_hessian = np.asarray(internal_hessian, dtype=float)
        internal_hessian = (
            (row_signs / row_scales)[:, None]
            * source_hessian[np.ix_(row_permutation, row_permutation)]
            * (row_signs / row_scales)[None, :]
        )
    provenance = json.dumps(
        {
            "provider": "validated_symmetry_transform",
            "source_request_id": request.source_request_id,
            "symmetry_operation_id": operation.operation_id,
            "transformation_validation_hash": operation.validation_hash,
        },
        sort_keys=True,
    )
    return GroupedCalculationResult(
        request_id=request.request_id,
        converged=True,
        cartesian_coordinates=tuple(tuple(float(x) for x in row) for row in coordinates),
        energy=float(source_result["energy"]),
        cartesian_gradient=tuple(tuple(float(x) for x in row) for row in cart_gradient),
        cartesian_hessian=None
        if cart_hessian is None
        else tuple(tuple(float(x) for x in row) for row in cart_hessian),
        internal_coordinates=tuple(float(x) for x in internal_coordinates),
        internal_gradient=tuple(float(x) for x in internal_gradient),
        internal_hessian=None
        if internal_hessian is None
        else tuple(tuple(float(x) for x in row) for row in internal_hessian),
        measured_phase_signature=tuple(float(x) for x in measured_phase_signature),
        constraint_errors=(),
        coordinate_fingerprint=request.coordinate_fingerprint,
        source_point_label=source_result["source_point_label"],
        provenance=provenance,
    )


def transform_internal_taylor_state(source_state, operation):
    """Apply an exact symmetry map to an interaction-residual Taylor state."""

    row_permutation = np.asarray(operation["row_permutation"], dtype=int)
    row_signs = np.asarray(operation["row_signs"], dtype=float)
    row_scales = _scales_of(operation, row_permutation.size)
    transformed = dict(source_state)
    transformed["request_id"] = (
        f"{source_state['request_id']}::{operation['operation_id']}"
    )
    transformed["internal_coordinates"] = (
        row_signs
        * row_scales
        * np.asarray(source_state["internal_coordinates"], dtype=float)[
            row_permutation
        ]
        + np.asarray(operation["row_offsets"], dtype=float)
    )
    transformed["internal_gradient"] = (
        row_signs
        * np.asarray(source_state["internal_gradient"], dtype=float)[
            row_permutation
        ]
        / row_scales
    )
    source_hessian = np.asarray(source_state["internal_hessian"], dtype=float)
    transformed["internal_hessian"] = (
        (row_signs / row_scales)[:, None]
        * source_hessian[np.ix_(row_permutation, row_permutation)]
        * (row_signs / row_scales)[None, :]
    )
    transformed["is_virtual_image"] = True
    transformed["source_request_id"] = source_state["request_id"]
    transformed["symmetry_operation_id"] = operation["operation_id"]
    return transformed


def validate_operation(operation, source_result, physical_target_result, policy):
    """Measure rigid-phase coverage without vetoing an exact atom permutation."""

    atom_permutation = np.asarray(operation.atom_permutation, dtype=int)
    row_permutation = np.asarray(operation.row_permutation, dtype=int)
    row_signs = np.asarray(operation.row_signs, dtype=float)
    row_scales = _scales_of(operation, row_permutation.size)
    predicted_xyz = np.asarray(source_result["cartesian_coordinates"])[atom_permutation]
    target_xyz = np.asarray(physical_target_result["cartesian_coordinates"])
    geometry_error = float(np.max(np.linalg.norm(predicted_xyz - target_xyz, axis=1)))
    energy_error = abs(float(source_result["energy"]) - float(physical_target_result["energy"]))
    predicted_gradient = (
        row_signs
        * np.asarray(source_result["internal_gradient"])[row_permutation]
        / row_scales
    )
    gradient_error = float(
        np.max(np.abs(predicted_gradient - np.asarray(physical_target_result["internal_gradient"])))
    )
    hessian_error = 0.0
    source_hessian = source_result.get("internal_hessian")
    target_hessian = physical_target_result.get("internal_hessian")
    if source_hessian is not None and target_hessian is not None:
        predicted_hessian = (
            (row_signs / row_scales)[:, None]
            * np.asarray(source_hessian)[np.ix_(row_permutation, row_permutation)]
            * (row_signs / row_scales)[None, :]
        )
        hessian_error = float(
            np.max(np.abs(predicted_hessian - np.asarray(target_hessian)))
        )
    geometry_valid = geometry_error <= policy.geometry_symmetry_tolerance_bohr
    derivative_valid = (
        energy_error <= policy.energy_symmetry_tolerance_hartree
        and gradient_error <= policy.gradient_symmetry_tolerance_hartree_per_bohr
        and hessian_error <= policy.hessian_symmetry_tolerance_hartree_per_bohr2
    )
    coverage_validated = geometry_valid and derivative_valid
    return replace(
        operation,
        geometry_validated=(True if operation.exact_atom_permutation else geometry_valid),
        derivative_validated=(
            True if operation.exact_atom_permutation else derivative_valid
        ),
        max_geometry_error=geometry_error,
        max_energy_error=energy_error,
        max_gradient_error=gradient_error,
        max_hessian_error=hessian_error,
        coverage_validated=coverage_validated,
    )
