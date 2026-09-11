"""Deterministic grouped subset and physical-request planning."""

from __future__ import annotations

from dataclasses import asdict, dataclass
import hashlib
from itertools import combinations
import json

import numpy as np

from .models import (
    GroupedCalculationRequest,
    GroupedModelPolicy,
    LocalSubset,
    PrimitiveRotor,
)
from .coupling import classify_rotor_pair


@dataclass(frozen=True)
class GroupedConstructionPlan:
    model_id: str
    family_id: str
    anchor_request: GroupedCalculationRequest
    requests: tuple[GroupedCalculationRequest, ...]
    subsets: tuple[LocalSubset, ...]
    factor_state_requests: tuple[tuple[str, tuple[str, ...]], ...]
    selected_rotor_ids: tuple[str, ...]
    interaction_state_requests: tuple[
        tuple[str, tuple[str, str], tuple[tuple[str, str, str], ...]], ...
    ] = ()


def _phase_token(degrees: float) -> str:
    sign = "m" if degrees < 0.0 else "p"
    return f"{sign}{abs(float(degrees)):08.3f}".replace(".", "d")


def select_rotors(rotors, requested_motion_ids=None) -> tuple[PrimitiveRotor, ...]:
    registry = {rotor.rotor_id: rotor for rotor in rotors}
    if requested_motion_ids is None:
        selected = tuple(sorted(registry.values(), key=lambda rotor: rotor.rotor_id))
    else:
        missing = sorted(set(requested_motion_ids).difference(registry))
        if missing:
            raise ValueError(f"Requested grouped motions were not detected: {missing}")
        selected = tuple(registry[rotor_id] for rotor_id in requested_motion_ids)
    if not selected:
        raise ValueError("Grouped construction requires at least one validated primitive motion.")
    return selected


def rotor_images_cover_the_period(rotor, policy) -> bool:
    """Whether one sector plus cyclic images stands in for the full period.

    An exact atom permutation of a symmetric rotor leaves the energy unchanged
    and permutes the gradient and Hessian, so its images are genuine data and
    one sector is enough -- this is the default.

    The images relabel equivalent atoms while the phase grid rotates the group
    rigidly, and the two coincide only when the group's own geometry carries the
    symmetry.  A methoxy methyl does not: hyperconjugation shortens the anti C-H
    bond and narrows its O-C-H angle by several degrees.  Setting
    ``rigid_symmetry_defect_tolerance_bohr`` makes such a rotor fall back to a
    physically sampled full period, which only matters for a torsion that is
    driven rigidly rather than allowed to relax.
    """

    if not (policy.enforce_exact_permutation_symmetry and rotor.symmetry_order > 1):
        return False
    tolerance = policy.rigid_symmetry_defect_tolerance_bohr
    if tolerance is None:
        return True
    defect = float(getattr(rotor, "rigid_symmetry_defect_bohr", 0.0))
    return defect <= float(tolerance)


def training_phases_for_rotor(rotor, policy):
    """Return deterministic physical factor phases for one motion."""

    if rotor.motion_kind == "group_rotation":
        phases = policy.group_rotation_training_phases_degrees
    else:
        phases = policy.training_phases_degrees
    if not (policy.enforce_exact_permutation_symmetry and rotor.symmetry_order > 1):
        return tuple(float(value) for value in phases)
    period = 360.0 / float(rotor.symmetry_order)
    unique = []
    seen = set()
    for value in phases:
        reduced = float(value) % period
        if abs(reduced - period) < 1.0e-10 or abs(reduced) < 1.0e-10:
            reduced = 0.0
        key = round(reduced, 10)
        if key not in seen:
            seen.add(key)
            unique.append(reduced)
    if rotor_images_cover_the_period(rotor, policy):
        return tuple(unique)
    # The images are not on the sampled path, so replicate the sector grid over
    # the whole period and let every phase be a physical calculation.
    return tuple(
        sorted(
            value + sector * period
            for sector in range(int(rotor.symmetry_order))
            for value in unique
        )
    )


def held_out_phases_for_rotor(rotor, policy):
    if rotor.motion_kind == "group_rotation":
        return tuple(float(value) for value in policy.group_rotation_held_out_phases_degrees)
    return tuple(float(value) for value in policy.held_out_phases_degrees)


def build_construction_plan(
    *,
    coordinate_fingerprint: str,
    rotors,
    n_internal_coordinates: int,
    policy: GroupedModelPolicy,
    requested_motion_ids=None,
    method_id: str,
    method_version: str,
    electronic_state_id: str = "ground_state_0",
    family_id: str = "canonical_input_basin",
) -> GroupedConstructionPlan:
    selected = select_rotors(rotors, requested_motion_ids)
    all_active_rows = {
        int(row) for rotor in selected for row in rotor.torsion_rows
    }
    core_rows = set(range(int(n_internal_coordinates))).difference(all_active_rows)
    subsets = []
    for rotor in selected:
        active_rows = set(int(row) for row in rotor.torsion_rows)

        # active_rows still controls the phase/signature definition.
        # projector_rows controls which physical derivatives are retained.
        projector_rows = tuple(range(int(n_internal_coordinates)))
        environment_atoms = tuple(sorted(set(rotor.axis_atoms).difference(rotor.moving_side_atoms)))
        subsets.append(
            LocalSubset(
                subset_id=f"factor.{rotor.rotor_id}",
                canonical_subset_key=(rotor.rotor_id,),
                role="factor",
                parent_subset_ids=(),
                active_atoms=rotor.moving_side_atoms,
                environment_atoms=environment_atoms,
                relaxation_atoms=(),
                active_rows=tuple(sorted(active_rows)),
                response_rows=(),
                projector_rows=projector_rows,
                relaxation_policy_id=policy.relaxation_policy_id,
                projector_policy_id=policy.projector_policy_id,
                anchor_policy_id=policy.anchor_policy_id,
            )
        )

    identity_payload = {
        "coordinate_fingerprint": coordinate_fingerprint,
        "family_id": family_id,
        "rotors": [asdict(rotor) for rotor in selected],
        "policy": asdict(policy),
        "method_id": method_id,
        "method_version": method_version,
    }
    model_hash = hashlib.sha256(
        json.dumps(identity_payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest()[:16]
    model_id = f"grouped_rigid_{model_hash}"
    anchor_id = f"{model_id}.anchor"
    anchor_request = GroupedCalculationRequest(
        request_id=anchor_id,
        subset_id="core",
        purpose="canonical_anchor",
        phase_signature=(),
        source_geometry_id="input_geometry",
        coordinate_fingerprint=coordinate_fingerprint,
        relaxation_policy_id=policy.relaxation_policy_id,
        projector_policy_id="full_chart_v1",
        anchor_id=anchor_id,
        electronic_state_id=electronic_state_id,
        method_id=method_id,
        method_version=method_version,
        is_anchor=True,
        is_virtual_image=False,
        source_request_id=None,
        symmetry_operation_id=None,
    )
    requests = [anchor_request]
    factor_state_requests = []
    phase_request_ids = {}
    for rotor, subset in zip(selected, subsets):
        state_ids = []
        phase_ids = {}
        factor_phases = training_phases_for_rotor(rotor, policy)
        cyclic = policy.enforce_exact_permutation_symmetry and rotor.symmetry_order > 1
        period = 360.0 / float(rotor.symmetry_order) if cyclic else None
        for phase_degrees in factor_phases:
            phase = float(np.deg2rad(phase_degrees))
            if abs(float(phase_degrees)) < 1.0e-12:
                state_ids.append(anchor_id)
                phase_ids[0.0] = anchor_id
                continue
            # A phase sitting on a full period is the rigidly rotated partner of
            # the anchor.  It doubles as the reference that measures how well the
            # atom permutation reproduces that rotation, so it keeps the coverage
            # name -- but it is a factor state like any other, because the rigid
            # structure and the relabelled one only coincide for an internally
            # symmetric group.
            on_period = cyclic and abs(float(phase_degrees) % period) < 1.0e-10
            request_id = (
                f"{model_id}.{rotor.rotor_id}."
                f"{'coverage' if on_period else 'phase'}_{_phase_token(phase_degrees)}"
            )
            request = GroupedCalculationRequest(
                request_id=request_id,
                subset_id=subset.subset_id,
                purpose=(
                    "rigid_symmetry_coverage_validation"
                    if on_period
                    else "factor_training"
                ),
                phase_signature=(phase,),
                source_geometry_id=anchor_id,
                coordinate_fingerprint=coordinate_fingerprint,
                relaxation_policy_id=policy.relaxation_policy_id,
                projector_policy_id=policy.projector_policy_id,
                anchor_id=anchor_id,
                electronic_state_id=electronic_state_id,
                method_id=method_id,
                method_version=method_version,
                is_anchor=False,
                is_virtual_image=False,
                source_request_id=None,
                symmetry_operation_id=None,
            )
            requests.append(request)
            state_ids.append(request_id)
            phase_ids[round(float(phase_degrees) % 360.0, 10)] = request_id

        # The rigidly rotated one-period structure is always calculated: it is
        # the reference the symmetry validation is measured against.  Whether it
        # also joins the bank is a trade: it is a free physical point, but the
        # factor is then no longer closed under the rotor's own permutation
        # group and the model stops being exactly symmetric under it.
        if cyclic and round(period, 10) not in phase_ids:
            request_id = (
                f"{model_id}.{rotor.rotor_id}.coverage_"
                f"{_phase_token(period)}"
            )
            if policy.bank_includes_rigid_coverage_state:
                state_ids.append(request_id)
                phase_ids[round(period, 10)] = request_id
            requests.append(
                GroupedCalculationRequest(
                    request_id=request_id,
                    subset_id=subset.subset_id,
                    purpose="rigid_symmetry_coverage_validation",
                    phase_signature=(float(np.deg2rad(period)),),
                    source_geometry_id=anchor_id,
                    coordinate_fingerprint=coordinate_fingerprint,
                    relaxation_policy_id=policy.relaxation_policy_id,
                    projector_policy_id=policy.projector_policy_id,
                    anchor_id=anchor_id,
                    electronic_state_id=electronic_state_id,
                    method_id=method_id,
                    method_version=method_version,
                    is_anchor=False,
                    is_virtual_image=False,
                    source_request_id=None,
                    symmetry_operation_id=None,
                )
            )

        factor_state_requests.append((subset.subset_id, tuple(state_ids)))
        phase_request_ids[rotor.rotor_id] = phase_ids

    # Nested motions share moving atoms and receive a deterministic candidate
    # interaction bank.  Whether it becomes active is still decided from the
    # persisted nonadditivity probes.
    interaction_state_requests = []
    for rotor_a, rotor_b in combinations(selected, 2):
        pair = classify_rotor_pair(rotor_a, rotor_b)
        if (
            not policy.construct_candidate_interaction_banks
            or pair["relation"] != "nested_overlap"
        ):
            continue
        rotor_ids = tuple(sorted((rotor_a.rotor_id, rotor_b.rotor_id)))
        by_id = {rotor_a.rotor_id: rotor_a, rotor_b.rotor_id: rotor_b}
        first, second = (by_id[rotor_ids[0]], by_id[rotor_ids[1]])
        subset_id = f"interaction.{pair['pair_id']}"
        active_rows = tuple(sorted(set(first.torsion_rows).union(second.torsion_rows)))
        # active_rows defines the pair phase; the projector must still span the
        # whole chart, exactly as it does for a factor.  Restricting it to the
        # coupled torsions discards every other component of the residual
        # gradient, which leaves the assembled gradient wrong even at a
        # physical bank centre where the energy is reproduced exactly.
        projector_rows = tuple(range(int(n_internal_coordinates)))
        factor_parent_ids = tuple(f"factor.{rotor_id}" for rotor_id in rotor_ids)
        subsets.append(
            LocalSubset(
                subset_id=subset_id,
                canonical_subset_key=rotor_ids,
                role="interaction",
                parent_subset_ids=factor_parent_ids,
                active_atoms=tuple(sorted(
                    set(first.moving_side_atoms).union(second.moving_side_atoms)
                )),
                environment_atoms=tuple(sorted(
                    set(first.axis_atoms).union(second.axis_atoms)
                    .difference(first.moving_side_atoms)
                    .difference(second.moving_side_atoms)
                )),
                relaxation_atoms=(),
                active_rows=active_rows,
                response_rows=(),
                projector_rows=projector_rows,
                relaxation_policy_id=policy.relaxation_policy_id,
                projector_policy_id=policy.projector_policy_id,
                anchor_policy_id=policy.anchor_policy_id,
            )
        )
        records = []
        phases_first = training_phases_for_rotor(first, policy)
        phases_second = training_phases_for_rotor(second, policy)
        for phase_first in phases_first:
            for phase_second in phases_second:
                id_first = phase_request_ids[first.rotor_id][
                    round(float(phase_first) % 360.0, 10)
                ]
                id_second = phase_request_ids[second.rotor_id][
                    round(float(phase_second) % 360.0, 10)
                ]
                if abs(phase_first) < 1.0e-12 and abs(phase_second) < 1.0e-12:
                    joint_id = anchor_id
                elif abs(phase_second) < 1.0e-12:
                    joint_id = id_first
                elif abs(phase_first) < 1.0e-12:
                    joint_id = id_second
                else:
                    joint_id = (
                        f"{model_id}.{subset_id}.phase_"
                        f"{_phase_token(phase_first)}_"
                        f"{_phase_token(phase_second)}"
                    )
                    requests.append(
                        GroupedCalculationRequest(
                            request_id=joint_id,
                            subset_id=subset_id,
                            purpose="candidate_pair_interaction_training",
                            phase_signature=(
                                float(np.deg2rad(phase_first)),
                                float(np.deg2rad(phase_second)),
                            ),
                            source_geometry_id=anchor_id,
                            coordinate_fingerprint=coordinate_fingerprint,
                            relaxation_policy_id=policy.relaxation_policy_id,
                            projector_policy_id=policy.projector_policy_id,
                            anchor_id=anchor_id,
                            electronic_state_id=electronic_state_id,
                            method_id=method_id,
                            method_version=method_version,
                            is_anchor=False,
                            is_virtual_image=False,
                            source_request_id=None,
                            symmetry_operation_id=None,
                        )
                    )
                records.append((joint_id, id_first, id_second))
        interaction_state_requests.append(
            (subset_id, rotor_ids, tuple(records))
        )

    return GroupedConstructionPlan(
        model_id=model_id,
        family_id=family_id,
        anchor_request=anchor_request,
        requests=tuple(requests),
        subsets=tuple(subsets),
        factor_state_requests=tuple(factor_state_requests),
        selected_rotor_ids=tuple(rotor.rotor_id for rotor in selected),
        interaction_state_requests=tuple(interaction_state_requests),
    )
