"""Orchestration of restartable rigid grouped-model construction."""

from __future__ import annotations

from dataclasses import asdict, replace
from itertools import combinations
import json
from pathlib import Path

import numpy as np

from .calculation_provider import GroupedCalculationProvider
from .coupling import (
    classify_rotor_pair,
    coupling_decision,
    order_pair_for_joint_generation,
    pair_nonadditivity,
)
from .coordinates import (
    build_coordinate_definition,
    build_interpolation_datapoint,
    symmetrize_equilibrium_bond_lengths,
)
from .models import GroupedCalculationRequest, GroupedModelPolicy
from .motion_detection import PrimitiveRotorDetector
from .phase import principal_periodic_delta
from .planning import (
    build_construction_plan,
    held_out_phases_for_rotor,
    rotor_images_cover_the_period,
    training_phases_for_rotor,
)
from .registry import GROUPED_SCHEMA_VERSION, GroupedModelRegistry
from .runtime import GroupedRuntimeModel
from .structure_generation import RigidRotorStructureGenerator, clone_molecule_with_coordinates
from .symmetry import (
    discover_cyclic_operations,
    transform_result,
    validate_operation,
)
from .validation import (
    evaluate_runtime_at_coordinates,
    finite_difference_runtime_gradient,
    validate_training_points,
)


def _driver_identity(driver):
    cls = driver.__class__
    method_id = f"{cls.__module__}.{cls.__name__}"
    version = getattr(driver, "version", None) or getattr(driver, "program_version", None)
    return method_id, "unknown" if version is None else str(version)


def _phase_key(degrees):
    """Round a phase onto the circle, so that 359.999999... is zero and not 360."""

    return round(float(degrees) % 360.0, 6) % 360.0


def _request_map(plan):
    return {request.request_id: request for request in plan.requests}


class GroupedModelConstructor:
    """Build a signed rigid-residual model through normal force-field setup."""

    def __init__(self, force_field_generator):
        self.generator = force_field_generator
        self.comm = force_field_generator.comm
        self.rank = force_field_generator.rank

    def _policy(self):
        generator = self.generator
        return GroupedModelPolicy(
            training_phases_degrees=tuple(
                float(value) for value in generator.grouped_training_phases_degrees
            ),
            held_out_phases_degrees=tuple(
                float(value) for value in generator.grouped_held_out_phases_degrees
            ),
            group_rotation_training_phases_degrees=tuple(
                float(value)
                for value in generator.grouped_group_rotation_training_phases_degrees
            ),
            group_rotation_held_out_phases_degrees=tuple(
                float(value)
                for value in generator.grouped_group_rotation_held_out_phases_degrees
            ),
            enforce_exact_permutation_symmetry=bool(
                generator.grouped_enforce_exact_permutation_symmetry
            ),
            construct_candidate_interaction_banks=bool(
                generator.grouped_run_validation_calculations
            ),
            shepard_p=float(generator.exponent_p),
            shepard_q=float(generator.exponent_q),
            confidence_radius=float(generator.grouped_confidence_radius),
            exact_center_tolerance=float(generator.grouped_exact_center_tolerance),
            phase_tolerance=float(generator.grouped_phase_tolerance),
            signature_scale=float(generator.grouped_signature_scale),
            geometry_symmetry_tolerance_bohr=float(generator.grouped_symmetry_geometry_tolerance_bohr),
            energy_symmetry_tolerance_hartree=float(generator.grouped_symmetry_energy_tolerance_hartree),
            gradient_symmetry_tolerance_hartree_per_bohr=float(generator.grouped_symmetry_gradient_tolerance),
            hessian_symmetry_tolerance_hartree_per_bohr2=float(generator.grouped_symmetry_hessian_tolerance),
            coupling_energy_threshold_hartree=float(generator.grouped_coupling_energy_threshold_hartree),
            coupling_gradient_rms_threshold=float(generator.grouped_coupling_gradient_rms_threshold),
            coupling_probe_phases_degrees=tuple(
                tuple(float(value) for value in pair)
                for pair in generator.grouped_coupling_probe_phases_degrees
            ),
            sidechain_basin_id=str(generator.grouped_sidechain_basin_id),
        )

    def construct(self, molecule, states_basis=None):
        generator = self.generator
        root = generator.roots_to_follow[0]
        filename = generator.imforcefieldfiles[root]
        settings = dict(generator.states_interpolation_settings[root])
        settings["construction_mode"] = "grouped_symmetry_aware"
        settings["grouped_model_id"] = getattr(generator, "grouped_model_id", None)
        z_matrix = generator.roots_z_matrix[root]
        policy = self._policy()

        eq_bond_lengths = np.array(
            [
                molecule.get_distance([int(a) + 1, int(b) + 1], "bohr")
                for a, b in z_matrix["bonds"]
            ],
            dtype=float,
        )
        coordinate_definition = build_coordinate_definition(
            z_matrix,
            use_inverse_bond_length=bool(generator.use_inverse_bond_length),
            use_eq_bond_length=bool(generator.use_eq_bond_length),
            use_cos_angle=bool(generator.use_cos_angle),
            use_mass_weight=bool(generator.use_mass_weight),
            eq_bond_lengths=eq_bond_lengths if generator.use_eq_bond_length else None,
        )

        def coordinate_values(coordinates):
            point = build_interpolation_datapoint(
                z_matrix,
                coordinates,
                settings,
                masses=molecule.get_masses(),
                eq_bond_lengths=eq_bond_lengths,
            )
            return point.internal_coordinates_values

        detector = PrimitiveRotorDetector(
            excluded_axis_bonds=getattr(generator, "grouped_excluded_axis_bonds", ()),
            detect_group_rotors=bool(
                getattr(generator, "grouped_detect_group_rotors", False)
            ),
        )
        detected_rotors = detector.detect(
            molecule, coordinate_definition, coordinate_values
        )
        if generator.use_eq_bond_length:
            # The chart has to commute with the permutations it is about to be
            # asked to trust, so equivalent bond rows share one reference length.
            # Detection only reads torsion rows, so the rotors carry over.
            symmetrized = symmetrize_equilibrium_bond_lengths(
                z_matrix, eq_bond_lengths, detected_rotors
            )
            if not np.allclose(symmetrized, eq_bond_lengths, rtol=0.0, atol=0.0):
                eq_bond_lengths = symmetrized
                coordinate_definition = build_coordinate_definition(
                    z_matrix,
                    use_inverse_bond_length=bool(generator.use_inverse_bond_length),
                    use_eq_bond_length=True,
                    use_cos_angle=bool(generator.use_cos_angle),
                    use_mass_weight=bool(generator.use_mass_weight),
                    eq_bond_lengths=eq_bond_lengths,
                )
                detected_rotors = detector.detect(
                    molecule, coordinate_definition, coordinate_values
                )
        detected_rotors = tuple(
            replace(
                rotor,
                signature_row_scales=(policy.signature_scale,) * len(rotor.signature_rows),
            )
            for rotor in detected_rotors
        )
        requested = getattr(generator, "requested_motion_ids", None)
        method_id, method_version = _driver_identity(generator.drivers["gs"][0])
        plan = build_construction_plan(
            coordinate_fingerprint=coordinate_definition.fingerprint,
            rotors=detected_rotors,
            n_internal_coordinates=len(coordinate_definition.rows),
            policy=policy,
            requested_motion_ids=requested,
            method_id=method_id,
            method_version=method_version,
            family_id=policy.sidechain_basin_id,
        )
        generator.grouped_model_id = plan.model_id
        settings["grouped_model_id"] = plan.model_id
        generator.states_interpolation_settings[root].update(
            {
                "construction_mode": "grouped_symmetry_aware",
                "grouped_model_id": plan.model_id,
            }
        )
        registry = GroupedModelRegistry(filename, plan.model_id)
        metadata = {
            "core_point_label": plan.anchor_request.request_id,
            "atom_labels": list(molecule.get_labels()),
            "charge": int(molecule.get_charge()),
            "multiplicity": int(molecule.get_multiplicity()),
            "basis_set_label": str(generator.gs_basis_set_label),
            "method_id": method_id,
            "method_version": method_version,
            "eq_bond_lengths_bohr": eq_bond_lengths.tolist(),
            "physical_cartesian_derivatives": True,
            "sidechain_basin_id": policy.sidechain_basin_id,
        }
        if self.rank == 0:
            status = registry.initialize(
                plan=plan,
                coordinate_definition=coordinate_definition,
                rotors=tuple(
                    rotor for rotor in detected_rotors
                    if rotor.rotor_id in plan.selected_rotor_ids
                ),
                policy=policy,
                metadata=metadata,
            )
        else:
            status = None
        status = self.comm.bcast(status, root=0)
        if status == "published":
            bundle = registry.load_bundle(require_published=True)
            validation = bundle["validation_json"]
            return {
                "n_datapoints": {root: validation["physical_point_count"]},
                "grouped_model_id": plan.model_id,
                "grouped_model_status": "published",
                "validation_report": validation,
            }

        selected_rotors = tuple(
            rotor for rotor in detected_rotors if rotor.rotor_id in plan.selected_rotor_ids
        )
        rotor_by_id = {rotor.rotor_id: rotor for rotor in selected_rotors}
        rotor_by_subset = {
            subset.subset_id: rotor_by_id[subset.canonical_subset_key[0]]
            for subset in plan.subsets
            if subset.role == "factor"
        }
        subset_by_id = {subset.subset_id: subset for subset in plan.subsets}
        anchor_q = coordinate_values(molecule.get_coordinates_in_bohr())
        provider = GroupedCalculationProvider(
            energy_function=generator._compute_energy,
            gradient_function=generator._compute_gradient,
            hessian_function=generator._compute_hessian,
            drivers=generator.drivers["gs"],
            z_matrix=z_matrix,
            interpolation_settings=settings,
            coordinate_fingerprint=coordinate_definition.fingerprint,
            eq_bond_lengths=eq_bond_lengths,
            anchor_internal_coordinates=anchor_q,
            rotor_by_subset=rotor_by_subset,
        )
        structure_generator = RigidRotorStructureGenerator(
            coordinate_values, phase_tolerance=policy.phase_tolerance
        )
        requests_by_id = _request_map(plan)
        request_molecules = {plan.anchor_request.request_id: molecule}
        for request in plan.requests:
            if request.is_anchor:
                continue
            subset = subset_by_id[request.subset_id]
            if subset.role == "factor":
                request_molecules[request.request_id] = structure_generator.generate(
                    molecule,
                    rotor_by_subset[request.subset_id],
                    request.phase_signature[0],
                )
                continue
            if subset.role == "interaction":
                rotor_a, rotor_b = (
                    rotor_by_id[rotor_id] for rotor_id in subset.canonical_subset_key
                )
                phase_by_id = dict(zip(
                    subset.canonical_subset_key, request.phase_signature
                ))
                outer, inner = order_pair_for_joint_generation(rotor_a, rotor_b)
                outer_molecule = structure_generator.generate(
                    molecule, outer, phase_by_id[outer.rotor_id]
                )
                request_molecules[request.request_id] = structure_generator.generate(
                    outer_molecule, inner, phase_by_id[inner.rotor_id]
                )
                continue
            raise ValueError(f"Unsupported request subset role: {subset.role}")

        basis = None
        if "XtbDriver" not in generator.drivers["gs"][0].__class__.__name__:
            from ..molecularbasis import MolecularBasis

            basis_label = generator.gs_basis_set_label if states_basis is None else states_basis.get("gs", generator.gs_basis_set_label)
            basis = MolecularBasis.read(molecule, basis_label)

        for request in plan.requests:
            if self.rank == 0:
                complete = registry.request_is_complete(request.request_id)
            else:
                complete = None
            complete = self.comm.bcast(complete, root=0)
            if complete:
                continue
            result = provider.compute(
                request, request_molecules[request.request_id], basis=basis, require_hessian=True
            )
            if not result.converged:
                raise RuntimeError(f"Grouped request did not converge: {request.request_id}")
            if result.constraint_errors and max(result.constraint_errors, default=0.0) > policy.phase_tolerance:
                raise RuntimeError(f"Grouped request failed phase validation: {request.request_id}")
            if self.rank == 0:
                registry.write_result(request, result)
            self.comm.barrier()

        # Symmetry is proposed from topology/chart closure, then checked against
        # independently calculated one-period rigid states.  An exact
        # element-preserving permutation is never vetoed by imperfect rigid
        # coverage; the latter is retained as a diagnostic.
        operations = []
        if self.rank == 0:
            anchor_result = registry.read_result(plan.anchor_request.request_id)
            for rotor in selected_rotors:
                if rotor.symmetry_order <= 1:
                    continue
                proposed = discover_cyclic_operations(
                    rotor,
                    coordinate_definition,
                    molecule.get_coordinates_in_bohr(),
                    coordinate_values,
                    atom_labels=molecule.get_labels(),
                )
                subset_id = f"factor.{rotor.rotor_id}"
                period_degrees = 360.0 / float(rotor.symmetry_order)
                physical_period = next(
                    request for request in plan.requests
                    if request.subset_id == subset_id
                    and request.purpose == "rigid_symmetry_coverage_validation"
                    and abs(
                        np.rad2deg(request.phase_signature[0]) - period_degrees
                    ) < 1.0e-8
                )
                result_period = registry.read_result(physical_period.request_id)
                operations.append(
                    validate_operation(
                        proposed[0], anchor_result, result_period, policy
                    )
                )
                operations.append(
                    validate_operation(
                        proposed[1], result_period, anchor_result, policy
                    )
                )
            registry.write_symmetry_registry(operations)
        operations = self.comm.bcast(operations if self.rank == 0 else None, root=0)

        if self.rank == 0:
            factor_registry = {
                subset_id: list(request_ids)
                for subset_id, request_ids in plan.factor_state_requests
            }
            for rotor in selected_rotors:
                if rotor.symmetry_order <= 1:
                    continue
                if not rotor_images_cover_the_period(rotor, policy):
                    # Every phase of this rotor is calculated, and its images sit
                    # off the rigid path by the group's own distortion, so adding
                    # them would only dilute physical states.
                    continue
                subset_id = f"factor.{rotor.rotor_id}"
                rotor_operations = sorted(
                    (
                        operation for operation in operations
                        if operation.source_rotor_id == rotor.rotor_id
                        and operation.validated
                    ),
                    key=lambda operation: operation.phase_offset,
                )
                source_ids = tuple(factor_registry[subset_id])
                existing_phases = set()
                for request_id in source_ids:
                    request = requests_by_id[request_id]
                    phase = 0.0 if request.is_anchor else np.rad2deg(
                        request.phase_signature[0]
                    )
                    existing_phases.add(_phase_key(phase))
                for source_request_id in source_ids:
                    source_request = requests_by_id[source_request_id]
                    base_degrees = 0.0 if source_request.is_anchor else np.rad2deg(
                        source_request.phase_signature[0]
                    )
                    for operation in rotor_operations:
                        target_degrees = (
                            base_degrees + np.rad2deg(operation.phase_offset)
                        ) % 360.0
                        target_key = _phase_key(target_degrees)
                        if target_key in existing_phases:
                            continue
                        request_id = (
                            f"{plan.model_id}.{rotor.rotor_id}.virtual_"
                            f"{int(round(target_degrees)) % 360:03d}"
                        )
                        if request_id in factor_registry[subset_id]:
                            continue
                        virtual_request = GroupedCalculationRequest(
                            request_id=request_id,
                            subset_id=subset_id,
                            purpose="exact_permutation_factor_state",
                            phase_signature=(float(principal_periodic_delta(
                                np.deg2rad(target_degrees))),),
                            source_geometry_id=source_request_id,
                            coordinate_fingerprint=coordinate_definition.fingerprint,
                            relaxation_policy_id=policy.relaxation_policy_id,
                            projector_policy_id=policy.projector_policy_id,
                            anchor_id=plan.anchor_request.request_id,
                            electronic_state_id="ground_state_0",
                            method_id=method_id,
                            method_version=method_version,
                            is_anchor=False,
                            is_virtual_image=True,
                            source_request_id=source_request_id,
                            symmetry_operation_id=operation.operation_id,
                        )
                        if not registry.request_is_complete(request_id):
                            virtual_result = transform_result(
                                virtual_request,
                                registry.read_result(source_request_id),
                                operation,
                                virtual_request.phase_signature,
                            )
                            registry.write_result(virtual_request, virtual_result)
                        factor_registry[subset_id].append(request_id)
                        existing_phases.add(target_key)
            registry.write_factor_registry(factor_registry)
        self.comm.barrier()

        runtime = GroupedRuntimeModel.from_hdf5(
            filename,
            model_id=plan.model_id,
            z_matrix=z_matrix,
            interpolation_settings=settings,
            require_published=False,
        )
        factor_training_ids = {
            request_id
            for _, request_ids in plan.factor_state_requests
            for request_id in request_ids
        }
        physical_results = [
            registry.read_result(request.request_id)
            for request in plan.requests
            if request.request_id in factor_training_ids
        ]
        training_records = validate_training_points(
            runtime, molecule, z_matrix, settings, physical_results
        )
        selected_dofs = tuple(
            3 * atom + component
            for rotor in selected_rotors
            for atom in rotor.owned_atoms[:1]
            for component in range(3)
        )
        held_out_records = []
        coupling_records = []
        pair_definitions = [
            classify_rotor_pair(rotor_a, rotor_b)
            for rotor_a, rotor_b in combinations(selected_rotors, 2)
        ]
        coupling_decisions = [
            {
                **pair,
                "decision": "NOT_PROBED",
                "interaction_bank_required": False,
                "number_of_probes": 0,
            }
            for pair in pair_definitions
        ]
        if generator.grouped_run_validation_calculations:
            for rotor in selected_rotors:
                subset_id = f"factor.{rotor.rotor_id}"
                for phase_degrees in held_out_phases_for_rotor(rotor, policy):
                    phase = float(np.deg2rad(phase_degrees))
                    validation_request = replace(
                        plan.anchor_request,
                        request_id=f"validation.{rotor.rotor_id}.{phase_degrees:.3f}",
                        subset_id=subset_id,
                        purpose="held_out_factor_validation",
                        phase_signature=(phase,),
                        is_anchor=False,
                    )
                    validation_molecule = structure_generator.generate(molecule, rotor, phase)
                    reference = provider.compute(
                        validation_request,
                        validation_molecule,
                        basis=basis,
                        require_hessian=False,
                    )
                    predicted_energy, predicted_gradient = evaluate_runtime_at_coordinates(
                        runtime,
                        molecule,
                        reference.cartesian_coordinates,
                        z_matrix,
                        settings,
                    )
                    gradient_error = predicted_gradient - np.asarray(reference.cartesian_gradient)
                    query_point = build_interpolation_datapoint(
                        z_matrix,
                        reference.cartesian_coordinates,
                        settings,
                        masses=molecule.get_masses(),
                        eq_bond_lengths=runtime.eq_bond_lengths,
                    )
                    no_virtual_energy, _ = runtime.evaluate_internal(
                        query_point.internal_coordinates_values,
                        include_virtual=False,
                    )
                    held_out_records.append(
                        {
                            "rotor_id": rotor.rotor_id,
                            "phase_degrees": float(phase_degrees),
                            "reference_energy_hartree": reference.energy,
                            "interpolated_energy_hartree": predicted_energy,
                            "energy_error_hartree": predicted_energy - reference.energy,
                            "gradient_rms_error_hartree_per_bohr": float(np.sqrt(np.mean(gradient_error**2))),
                            "maximum_force_error_hartree_per_bohr": float(np.max(np.abs(gradient_error))),
                            "error_decomposition_hartree": {
                                "assembled_model": predicted_energy - reference.energy,
                                "local_factor_without_virtualization": no_virtual_energy - reference.energy,
                                "symmetry_virtualization": predicted_energy - no_virtual_energy,
                                "taylor_truncation": None,
                            },
                        }
                    )
            coupling_decisions = []
            for (rotor_a, rotor_b), pair in zip(
                combinations(selected_rotors, 2), pair_definitions
            ):
                pair_records = []
                outer, inner = order_pair_for_joint_generation(rotor_a, rotor_b)
                for phase_a_deg, phase_b_deg in policy.coupling_probe_phases_degrees:
                    molecule_a = structure_generator.generate(
                        molecule, rotor_a, np.deg2rad(phase_a_deg)
                    )
                    molecule_b = structure_generator.generate(
                        molecule, rotor_b, np.deg2rad(phase_b_deg)
                    )
                    if outer.rotor_id == rotor_a.rotor_id:
                        joint_molecule = structure_generator.generate(
                            molecule_a, inner, np.deg2rad(phase_b_deg)
                        )
                    else:
                        joint_molecule = structure_generator.generate(
                            molecule_b, inner, np.deg2rad(phase_a_deg)
                        )
                    probe_request = replace(
                        plan.anchor_request,
                        request_id=(
                            f"coupling.{pair['pair_id']}."
                            f"{phase_a_deg:.0f}.{phase_b_deg:.0f}"
                        ),
                        subset_id=f"interaction.{pair['pair_id']}",
                        purpose="pair_coupling_probe",
                        phase_signature=(np.deg2rad(phase_a_deg), np.deg2rad(phase_b_deg)),
                        is_anchor=False,
                    )
                    reference_a = provider.compute(
                        replace(
                            probe_request,
                            request_id=f"{probe_request.request_id}.single_a",
                            subset_id=f"factor.{rotor_a.rotor_id}",
                            phase_signature=(np.deg2rad(phase_a_deg),),
                        ),
                        molecule_a,
                        basis=basis,
                        require_hessian=False,
                    )
                    reference_b = provider.compute(
                        replace(
                            probe_request,
                            request_id=f"{probe_request.request_id}.single_b",
                            subset_id=f"factor.{rotor_b.rotor_id}",
                            phase_signature=(np.deg2rad(phase_b_deg),),
                        ),
                        molecule_b,
                        basis=basis,
                        require_hessian=False,
                    )
                    reference_joint = provider.compute(
                        probe_request, joint_molecule, basis=basis, require_hessian=False
                    )
                    predicted_energy, predicted_gradient = evaluate_runtime_at_coordinates(
                        runtime,
                        molecule,
                        reference_joint.cartesian_coordinates,
                        z_matrix,
                        settings,
                    )
                    gradient_error = predicted_gradient - np.asarray(
                        reference_joint.cartesian_gradient
                    )
                    metrics = pair_nonadditivity(
                        {
                            "energy": reference_joint.energy,
                            "cartesian_gradient": reference_joint.cartesian_gradient,
                        },
                        {
                            "energy": reference_a.energy,
                            "cartesian_gradient": reference_a.cartesian_gradient,
                        },
                        {
                            "energy": reference_b.energy,
                            "cartesian_gradient": reference_b.cartesian_gradient,
                        },
                        registry.read_result(plan.anchor_request.request_id),
                    )
                    metrics.update(
                        {
                            "pair_id": pair["pair_id"],
                            "rotor_ids": list(pair["rotor_ids"]),
                            "relation": pair["relation"],
                            "phase_degrees": [phase_a_deg, phase_b_deg],
                            "assembled_model_energy_error_hartree": float(
                                predicted_energy - reference_joint.energy),
                            "assembled_model_gradient_rms_error_hartree_per_bohr": float(
                                np.sqrt(np.mean(gradient_error**2))),
                            "hessian_nonadditivity": None,
                        }
                    )
                    coupling_records.append(metrics)
                    pair_records.append(metrics)
                coupling_decisions.append(
                    coupling_decision(pair, pair_records, policy)
                )

        decisions_by_rotors = {
            tuple(decision["rotor_ids"]): decision
            for decision in coupling_decisions
        }
        selected_interactions = []
        selected_rotor_pairs = set()
        for subset_id, rotor_ids, records in plan.interaction_state_requests:
            decision = decisions_by_rotors.get(tuple(rotor_ids), {})
            if not decision.get("interaction_bank_required", False):
                continue
            selected_rotor_pairs.add(tuple(rotor_ids))
            selected_interactions.append(
                {
                    "subset_id": subset_id,
                    "pair_id": subset_id.removeprefix("interaction."),
                    "rotor_ids": list(rotor_ids),
                    "state_records": [
                        {
                            "joint_request_id": joint_id,
                            "state_a_request_id": state_a_id,
                            "state_b_request_id": state_b_id,
                        }
                        for joint_id, state_a_id, state_b_id in records
                    ],
                    "construction_identity": "joint_minus_additive_runtime_v1",
                }
            )
        coupling_decisions = [
            {
                **decision,
                "decision": (
                    "COUPLED_INTERACTION_INCLUDED"
                    if tuple(decision["rotor_ids"]) in selected_rotor_pairs
                    else decision["decision"]
                ),
                "interaction_bank_included": (
                    tuple(decision["rotor_ids"]) in selected_rotor_pairs
                ),
            }
            for decision in coupling_decisions
        ]
        if self.rank == 0:
            registry.write_interaction_registry(selected_interactions)
            registry.write_coupling_registry(coupling_decisions)
        self.comm.barrier()

        if selected_interactions:
            runtime = GroupedRuntimeModel.from_hdf5(
                filename,
                model_id=plan.model_id,
                z_matrix=z_matrix,
                interpolation_settings=settings,
                require_published=False,
            )
        finite_difference = finite_difference_runtime_gradient(
            runtime,
            molecule,
            z_matrix,
            settings,
            dof_indices=selected_dofs,
        )
        symmetry_sector_finite_difference = {
            "phase_degrees": None,
            "rotor_id": None,
            "step_bohr": 1.0e-4,
            "components": [],
            "maximum_absolute_error": 0.0,
        }
        symmetric_rotor = next(
            (rotor for rotor in selected_rotors if rotor.symmetry_order > 1),
            None,
        )
        if symmetric_rotor is not None:
            phase = -2.0 * np.pi / float(symmetric_rotor.symmetry_order)
            symmetry_molecule = structure_generator.generate(
                molecule, symmetric_rotor, phase
            )
            symmetry_sector_finite_difference = {
                "phase_degrees": float(np.rad2deg(phase) % 360.0),
                "rotor_id": symmetric_rotor.rotor_id,
                **finite_difference_runtime_gradient(
                    runtime,
                    symmetry_molecule,
                    z_matrix,
                    settings,
                    dof_indices=tuple(
                        3 * symmetric_rotor.owned_atoms[0] + component
                        for component in range(3)
                    ),
                ),
            }
        interaction_finite_difference = {
            "subset_id": None,
            "request_id": None,
            "step_bohr": 1.0e-4,
            "components": [],
            "maximum_absolute_error": 0.0,
        }
        if selected_interactions:
            interaction = selected_interactions[0]
            joint_record = next(
                record
                for record in interaction["state_records"]
                if record["joint_request_id"]
                not in {
                    record["state_a_request_id"],
                    record["state_b_request_id"],
                }
            )
            joint_result = registry.read_result(joint_record["joint_request_id"])
            joint_molecule = clone_molecule_with_coordinates(
                molecule, joint_result["cartesian_coordinates"]
            )
            first_pair_rotor = rotor_by_id[interaction["rotor_ids"][0]]
            interaction_finite_difference = {
                "subset_id": interaction["subset_id"],
                "request_id": joint_record["joint_request_id"],
                **finite_difference_runtime_gradient(
                    runtime,
                    joint_molecule,
                    z_matrix,
                    settings,
                    dof_indices=tuple(
                        3 * first_pair_rotor.owned_atoms[0] + component
                        for component in range(3)
                    ),
                ),
            }
        interaction_training_records = []
        for interaction in selected_interactions:
            for state_record in interaction["state_records"]:
                joint = registry.read_result(state_record["joint_request_id"])
                predicted_energy, predicted_gradient = evaluate_runtime_at_coordinates(
                    runtime,
                    molecule,
                    joint["cartesian_coordinates"],
                    z_matrix,
                    settings,
                )
                gradient_error = predicted_gradient - np.asarray(
                    joint["cartesian_gradient"], dtype=float
                )
                interaction_training_records.append(
                    {
                        "subset_id": interaction["subset_id"],
                        "request_id": state_record["joint_request_id"],
                        "energy_error_hartree": float(
                            predicted_energy - joint["energy"]
                        ),
                        "gradient_rms_error_hartree_per_bohr": float(
                            np.sqrt(np.mean(gradient_error**2))
                        ),
                        "maximum_force_error_hartree_per_bohr": float(
                            np.max(np.abs(gradient_error))
                        ),
                    }
                )

        point_index = registry.point_index()
        physical_count = len(plan.requests)
        virtual_count = len(point_index) - physical_count
        coupling_required = any(
            decision.get("interaction_bank_required", False)
            for decision in coupling_decisions
        )
        unresolved_coupling = any(
            decision.get("interaction_bank_required", False)
            and not decision.get("interaction_bank_included", False)
            for decision in coupling_decisions
        )
        max_training_energy = max(
            (abs(record["energy_error_hartree"]) for record in training_records), default=0.0
        )
        max_interaction_training_energy = max(
            (
                abs(record["energy_error_hartree"])
                for record in interaction_training_records
            ),
            default=0.0,
        )
        max_fd = max(
            finite_difference["maximum_absolute_error"],
            symmetry_sector_finite_difference["maximum_absolute_error"],
            interaction_finite_difference["maximum_absolute_error"],
        )
        publication_passed = (
            np.isfinite(max_training_energy)
            and max_training_energy <= 5.0e-8
            and np.isfinite(max_interaction_training_energy)
            and max_interaction_training_energy <= 5.0e-8
            and np.isfinite(max_fd)
            and max_fd <= 5.0e-4
            and not unresolved_coupling
        )
        phases_by_rotor = {
            rotor.rotor_id: list(training_phases_for_rotor(rotor, policy))
            for rotor in selected_rotors
        }
        held_out_by_rotor = {
            rotor.rotor_id: list(held_out_phases_for_rotor(rotor, policy))
            for rotor in selected_rotors
        }
        full_grid_count = int(np.prod([
            len(phases_by_rotor[rotor.rotor_id]) for rotor in selected_rotors
        ]))
        report = {
            "status": "passed" if publication_passed else "failed",
            "model_id": plan.model_id,
            "schema_version": GROUPED_SCHEMA_VERSION,
            "coordinate_fingerprint": coordinate_definition.fingerprint,
            "detected_rotors": [asdict(rotor) for rotor in detected_rotors],
            "selected_rotors": list(plan.selected_rotor_ids),
            "physical_point_count": physical_count,
            "virtual_point_count": virtual_count,
            "before_full_grid_point_count": full_grid_count,
            "after_grouped_physical_point_count": physical_count,
            "training_phases_degrees": list(policy.training_phases_degrees),
            "held_out_phases_degrees": list(policy.held_out_phases_degrees),
            "training_phases_by_rotor": phases_by_rotor,
            "held_out_phases_by_rotor": held_out_by_rotor,
            "training_errors": training_records,
            "interaction_training_errors": interaction_training_records,
            "held_out_errors": held_out_records,
            "symmetry_validation_errors": [
                {
                    "operation_id": operation.operation_id,
                    "geometry_validated": operation.geometry_validated,
                    "derivative_validated": operation.derivative_validated,
                    "exact_atom_permutation": operation.exact_atom_permutation,
                    "rigid_coverage_validated": operation.coverage_validated,
                    "max_geometry_error_bohr": operation.max_geometry_error,
                    "max_energy_error_hartree": operation.max_energy_error,
                    "max_gradient_error": operation.max_gradient_error,
                    "max_hessian_error": operation.max_hessian_error,
                    "validation_hash": operation.validation_hash,
                }
                for operation in operations
            ],
            "finite_difference_errors": finite_difference,
            "symmetry_sector_finite_difference_errors": (
                symmetry_sector_finite_difference
            ),
            "interaction_finite_difference_errors": (
                interaction_finite_difference
            ),
            "coupling_probe_results": coupling_records,
            "coupling_decisions": coupling_decisions,
            "coupling_decision": (
                "PAIR_INTERACTION_REQUIRED"
                if unresolved_coupling
                else (
                    "PAIR_INTERACTION_INCLUDED"
                    if coupling_required
                    else (
                        "DECOUPLED_MODEL_RETAINED"
                        if coupling_records
                        else "NOT_PROBED"
                    )
                )
            ),
            "sidechain_basin_id": policy.sidechain_basin_id,
            "publication_status": "published" if publication_passed else "staging",
            "content_hash": "",
            "limitations": [
                "rigid local factors only",
                "one canonical sidechain basin",
                "only nested-overlap interaction banks are constructed automatically",
            ],
        }
        if self.rank == 0:
            registry.update_metadata(
                {
                    "physical_point_count": physical_count,
                    "virtual_point_count": virtual_count,
                    "number_of_banks": len(plan.subsets),
                    "number_of_coupling_probes": len(coupling_records),
                    "number_of_selected_interaction_banks": len(selected_interactions),
                    "number_of_interaction_virtual_states": sum(
                        len(interaction["state_records"])
                        * sum(
                            1
                            for operation in operations
                            if operation.source_rotor_id
                            in interaction["rotor_ids"]
                            and operation.validated
                            and operation.exact_atom_permutation
                        )
                        for interaction in selected_interactions
                    ),
                }
            )
            registry.write_validation(report)
            if not publication_passed:
                raise RuntimeError(
                    f"Grouped model validation failed: max training energy error={max_training_energy:.3e}, "
                    f"max finite-difference error={max_fd:.3e}, "
                    f"unresolved pair interaction={unresolved_coupling}."
                )
            content_hash = registry.publish()
            report["content_hash"] = content_hash
            report_path = Path(
                generator.grouped_validation_report
                or f"{filename}.{plan.model_id}.validation.json"
            )
            report_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        else:
            report = None
        report = self.comm.bcast(report, root=0)
        return {
            "n_datapoints": {root: physical_count},
            "grouped_model_id": plan.model_id,
            "grouped_model_status": "published",
            "validation_report": report,
        }
