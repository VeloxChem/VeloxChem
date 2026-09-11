"""Rigid local structure generation using an axis-angle rotation."""

from __future__ import annotations

import numpy as np

from ..molecule import Molecule
from .phase import estimate_collective_phase, principal_periodic_delta


def rotate_coordinates_about_axis(
    coordinates,
    *,
    axis_atoms: tuple[int, int],
    moving_atoms,
    angle: float,
):
    """Return a copy rotated by ``angle`` around the directed bond axis."""

    xyz = np.asarray(coordinates, dtype=float)
    rotated = xyz.copy()
    origin = xyz[int(axis_atoms[0])]
    axis = xyz[int(axis_atoms[1])] - origin
    norm = float(np.linalg.norm(axis))
    if norm < 1.0e-12:
        raise ValueError("Cannot rotate around a zero-length axis.")
    unit = axis / norm
    cosine = float(np.cos(angle))
    sine = float(np.sin(angle))
    for atom in moving_atoms:
        atom = int(atom)
        if atom in axis_atoms:
            continue
        vector = xyz[atom] - origin
        rotated[atom] = origin + (
            vector * cosine
            + np.cross(unit, vector) * sine
            + unit * np.dot(unit, vector) * (1.0 - cosine)
        )
    return rotated


def clone_molecule_with_coordinates(molecule, coordinates_bohr):
    clone = Molecule(
        molecule.get_labels(), np.asarray(coordinates_bohr, dtype=float), "bohr"
    )
    clone.set_charge(molecule.get_charge())
    clone.set_multiplicity(molecule.get_multiplicity())
    return clone


class RigidRotorStructureGenerator:
    """Generate chart-validated rigid rotor states from a canonical anchor."""

    def __init__(self, coordinate_value_function, phase_tolerance=2.0e-5):
        self.coordinate_value_function = coordinate_value_function
        self.phase_tolerance = float(phase_tolerance)

    def generate(self, anchor_molecule, rotor, target_phase) -> Molecule:
        anchor_coordinates = np.asarray(
            anchor_molecule.get_coordinates_in_bohr(), dtype=float
        )
        q_anchor = np.asarray(
            self.coordinate_value_function(anchor_coordinates), dtype=float
        )
        current_phase = estimate_collective_phase(q_anchor, q_anchor, rotor)
        increment = principal_periodic_delta(float(target_phase) - current_phase)
        generated_coordinates = rotate_coordinates_about_axis(
            anchor_coordinates,
            axis_atoms=rotor.axis_atoms,
            moving_atoms=rotor.moving_side_atoms,
            angle=increment,
        )
        q_generated = np.asarray(
            self.coordinate_value_function(generated_coordinates), dtype=float
        )
        measured = estimate_collective_phase(q_generated, q_anchor, rotor)
        error = abs(principal_periodic_delta(measured - float(target_phase)))
        if error > self.phase_tolerance:
            raise RuntimeError(
                f"Rigid generation for {rotor.rotor_id} missed its phase by "
                f"{error:.3e} rad (tolerance {self.phase_tolerance:.3e})."
            )
        return clone_molecule_with_coordinates(anchor_molecule, generated_coordinates)
