"""Deterministic pair nonadditivity metrics and selection decisions."""

from __future__ import annotations

import hashlib
import json

import numpy as np


def classify_rotor_pair(rotor_a, rotor_b):
    """Classify two motions from canonical atom ownership, never raw names."""

    rotor_ids = tuple(sorted((rotor_a.rotor_id, rotor_b.rotor_id)))
    moving_a = set(rotor_a.moving_side_atoms)
    moving_b = set(rotor_b.moving_side_atoms)
    overlap = tuple(sorted(moving_a.intersection(moving_b)))
    if moving_a.issubset(moving_b) or moving_b.issubset(moving_a):
        relation = "nested_overlap"
    elif overlap:
        relation = "partial_overlap"
    else:
        relation = "disjoint"
    payload = {"rotor_ids": rotor_ids, "relation": relation}
    pair_hash = hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest()[:12]
    return {
        "pair_id": f"pair.{pair_hash}",
        "rotor_ids": rotor_ids,
        "relation": relation,
        "overlap_atoms": overlap,
        "probe_required": True,
    }


def order_pair_for_joint_generation(rotor_a, rotor_b):
    """Apply an outer group rotation before a motion nested inside it."""

    moving_a = set(rotor_a.moving_side_atoms)
    moving_b = set(rotor_b.moving_side_atoms)
    if moving_b < moving_a:
        return rotor_a, rotor_b
    if moving_a < moving_b:
        return rotor_b, rotor_a
    return tuple(sorted((rotor_a, rotor_b), key=lambda rotor: rotor.rotor_id))


def pair_nonadditivity(joint, state_a, state_b, anchor):
    delta_energy = (
        float(joint["energy"])
        - float(state_a["energy"])
        - float(state_b["energy"])
        + float(anchor["energy"])
    )
    delta_gradient = (
        np.asarray(joint["cartesian_gradient"], dtype=float)
        - np.asarray(state_a["cartesian_gradient"], dtype=float)
        - np.asarray(state_b["cartesian_gradient"], dtype=float)
        + np.asarray(anchor["cartesian_gradient"], dtype=float)
    )
    return {
        "energy_nonadditivity_hartree": delta_energy,
        "gradient_rms_nonadditivity_hartree_per_bohr": float(
            np.sqrt(np.mean(delta_gradient * delta_gradient))
        ),
        "maximum_gradient_nonadditivity_hartree_per_bohr": float(
            np.max(np.abs(delta_gradient))
        ),
    }


def pair_interaction_required(records, policy):
    return any(
        abs(record["energy_nonadditivity_hartree"])
        > policy.coupling_energy_threshold_hartree
        or record["gradient_rms_nonadditivity_hartree_per_bohr"]
        > policy.coupling_gradient_rms_threshold
        for record in records
    )


def coupling_decision(pair, records, policy):
    required = pair_interaction_required(records, policy)
    return {
        **pair,
        "decision": "COUPLED_INTERACTION_REQUIRED" if required else "DECOUPLED",
        "interaction_bank_required": bool(required),
        "maximum_absolute_energy_nonadditivity_hartree": max(
            (abs(record["energy_nonadditivity_hartree"]) for record in records),
            default=0.0,
        ),
        "maximum_gradient_rms_nonadditivity_hartree_per_bohr": max(
            (
                record["gradient_rms_nonadditivity_hartree_per_bohr"]
                for record in records
            ),
            default=0.0,
        ),
        "number_of_probes": len(records),
    }
