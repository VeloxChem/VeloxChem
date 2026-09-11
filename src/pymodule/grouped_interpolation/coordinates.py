"""Canonical coordinate-chart construction and conversion helpers."""

from __future__ import annotations

from dataclasses import asdict
import hashlib
import json

import numpy as np

from ..interpolationdatapoint import InterpolationDatapoint
from .models import CoordinateDefinition


COORDINATE_DEFINITION_VERSION = "grouped-coordinate-chart-v2"


def canonical_json(value) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _canonical_z_matrix(z_matrix) -> dict[str, list[tuple[int, ...]]]:
    required = ("bonds", "angles", "dihedrals", "impropers")
    return {
        key: [tuple(int(atom) for atom in row) for row in z_matrix.get(key, [])]
        for key in required
    }


def build_coordinate_definition(
    z_matrix,
    *,
    use_inverse_bond_length: bool,
    use_eq_bond_length: bool,
    use_cos_angle: bool,
    use_mass_weight: bool,
    eq_bond_lengths=None,
) -> CoordinateDefinition:
    """Build a content-addressed definition for the one accepted chart."""

    zmat = _canonical_z_matrix(z_matrix)
    rows = tuple(row for key in ("bonds", "angles", "dihedrals", "impropers") for row in zmat[key])
    row_types = (
        ("bond",) * len(zmat["bonds"])
        + ("angle",) * len(zmat["angles"])
        + ("proper_torsion",) * len(zmat["dihedrals"])
        + ("improper_torsion",) * len(zmat["impropers"])
    )
    if use_inverse_bond_length:
        bond_convention = "inverse_distance"
        bond_unit = "bohr^-1"
    elif use_eq_bond_length:
        bond_convention = "equilibrium_scaled_distance"
        bond_unit = "bohr"
    else:
        bond_convention = "distance"
        bond_unit = "bohr"

    units = tuple(
        bond_unit if kind == "bond" else "dimensionless" if kind == "angle" and use_cos_angle else "radian"
        for kind in row_types
    )
    periodicities = tuple(
        2.0 * np.pi if kind in {"proper_torsion", "improper_torsion"} else None
        for kind in row_types
    )
    sections = []
    start = 0
    for key in ("bonds", "angles", "dihedrals", "impropers"):
        end = start + len(zmat[key])
        sections.append((key, start, end))
        start = end

    equilibrium_lengths = ()
    if use_eq_bond_length:
        if eq_bond_lengths is None:
            raise ValueError(
                "Equilibrium bond lengths are part of an equilibrium-scaled coordinate chart."
            )
        equilibrium_lengths = tuple(float(value) for value in eq_bond_lengths)
        if len(equilibrium_lengths) != len(zmat["bonds"]):
            raise ValueError("Equilibrium bond-length count does not match the bond rows.")

    payload = {
        "rows": rows,
        "row_types": row_types,
        "units": units,
        "periodicities": periodicities,
        "torsion_convention": "sin_periodic_delta",
        "improper_convention": "principal_delta_then_2tan_half",
        "bond_convention": bond_convention,
        "angle_convention": "cosine" if use_cos_angle else "radian",
        "mass_weighting_convention": "cartesian_inverse_sqrt_mass" if use_mass_weight else "physical_cartesian",
        "equilibrium_bond_lengths_bohr": equilibrium_lengths,
        "version": COORDINATE_DEFINITION_VERSION,
        "sections": sections,
    }
    fingerprint = hashlib.sha256(canonical_json(payload).encode("utf-8")).hexdigest()
    return CoordinateDefinition(
        rows=rows,
        row_types=row_types,
        units=units,
        periodicities=periodicities,
        torsion_convention=payload["torsion_convention"],
        improper_convention=payload["improper_convention"],
        bond_convention=bond_convention,
        angle_convention=payload["angle_convention"],
        mass_weighting_convention=payload["mass_weighting_convention"],
        equilibrium_bond_lengths_bohr=equilibrium_lengths,
        version=COORDINATE_DEFINITION_VERSION,
        fingerprint=fingerprint,
        sections=tuple(sections),
    )


def symmetrize_equilibrium_bond_lengths(z_matrix, eq_bond_lengths, rotors):
    """Give bond rows exchanged by an exact rotor permutation one shared length.

    ``use_eq_bond_length`` stores ``r_eq - r_eq**2 / r`` per bond row, with the
    reference length read off the anchor.  The three C-H bonds of a methyl are
    never exactly equal there -- hyperconjugation alone splits them by more than
    a hundredth of a bohr -- so three rows that the C3 permutation exchanges get
    three different references.  The chart then fails to commute with the
    permutation: relabelling the hydrogens of a structure moves it in coordinate
    space even though nothing physical changed, and the exact symmetry image of
    a data point no longer lands where the chart says it should.  Averaging the
    reference over each permutation orbit restores ``q(pi X) = P q(X)`` exactly.
    """

    lengths = np.array(eq_bond_lengths, dtype=float)
    bonds = [tuple(int(atom) for atom in row) for row in z_matrix.get("bonds", [])]
    if lengths.size != len(bonds):
        raise ValueError("One equilibrium length per bond row is required.")
    index_of = {frozenset(bond): position for position, bond in enumerate(bonds)}

    parent = list(range(len(bonds)))

    def find(node):
        while parent[node] != node:
            parent[node] = parent[parent[node]]
            node = parent[node]
        return node

    def union(left, right):
        left, right = find(left), find(right)
        if left != right:
            parent[max(left, right)] = min(left, right)

    for rotor in rotors:
        if rotor.symmetry_order <= 1:
            continue
        for group in rotor.equivalent_atom_groups:
            group = tuple(int(atom) for atom in group)
            for shift in range(1, len(group)):
                mapping = {
                    group[position]: group[(position + shift) % len(group)]
                    for position in range(len(group))
                }
                for position, bond in enumerate(bonds):
                    image = frozenset(mapping.get(atom, atom) for atom in bond)
                    partner = index_of.get(image)
                    if partner is not None:
                        union(position, partner)

    orbits = {}
    for position in range(len(bonds)):
        orbits.setdefault(find(position), []).append(position)
    for members in orbits.values():
        if len(members) > 1:
            lengths[members] = float(np.mean(lengths[members]))
    return lengths


def coordinate_definition_to_dict(definition: CoordinateDefinition) -> dict:
    return asdict(definition)


def coordinate_definition_from_dict(payload: dict) -> CoordinateDefinition:
    converted = dict(payload)
    converted["rows"] = tuple(tuple(row) for row in converted["rows"])
    converted["row_types"] = tuple(converted["row_types"])
    converted["units"] = tuple(converted["units"])
    converted["periodicities"] = tuple(converted["periodicities"])
    converted["equilibrium_bond_lengths_bohr"] = tuple(
        converted.get("equilibrium_bond_lengths_bohr", ()))
    converted["sections"] = tuple(tuple(section) for section in converted["sections"])
    return CoordinateDefinition(**converted)


def build_interpolation_datapoint(
    z_matrix,
    coordinates_bohr,
    settings,
    *,
    masses=None,
    eq_bond_lengths=None,
) -> InterpolationDatapoint:
    """Create a chart-consistent coordinate object for one geometry."""

    point = InterpolationDatapoint(z_matrix)
    point.update_settings(settings)
    if settings.get("use_mass_weight", False):
        if masses is None:
            raise ValueError("Masses are required by the persisted mass-weighted chart.")
        point.inv_sqrt_masses = 1.0 / np.sqrt(np.repeat(np.asarray(masses, dtype=float), 3))
    point.eq_bond_lengths = None if eq_bond_lengths is None else np.asarray(eq_bond_lengths, dtype=float)
    point.reset_coordinates_impes_driver(np.asarray(coordinates_bohr, dtype=float))
    return point


def assert_coordinate_fingerprint(actual: str, expected: str, context: str) -> None:
    if actual != expected:
        raise ValueError(
            f"Coordinate fingerprint mismatch for {context}: expected {expected}, got {actual}."
        )
