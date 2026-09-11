"""Narrow adapter around the existing VeloxChem E/G/H calculation path."""

from __future__ import annotations

import json

import numpy as np

from ..interpolationdatapoint import InterpolationDatapoint
from .models import GroupedCalculationResult


class GroupedCalculationProvider:
    """Reuse configured driver dispatch and expose physical Cartesian results."""

    def __init__(
        self,
        *,
        energy_function,
        gradient_function,
        hessian_function,
        drivers,
        z_matrix,
        interpolation_settings,
        coordinate_fingerprint,
        eq_bond_lengths,
        anchor_internal_coordinates,
        rotor_by_subset,
    ):
        self.energy_function = energy_function
        self.gradient_function = gradient_function
        self.hessian_function = hessian_function
        self.drivers = drivers
        self.z_matrix = z_matrix
        self.interpolation_settings = dict(interpolation_settings)
        self.coordinate_fingerprint = str(coordinate_fingerprint)
        self.eq_bond_lengths = np.asarray(eq_bond_lengths, dtype=float)
        self.anchor_internal_coordinates = np.asarray(
            anchor_internal_coordinates, dtype=float
        )
        self.rotor_by_subset = dict(rotor_by_subset)

    def compute(self, request, molecule, basis=None, require_hessian=True):
        energy_values, scf_results, rsp_results = self.energy_function(
            self.drivers[0], molecule, basis
        )
        gradient_values = self.gradient_function(
            self.drivers[1], molecule, basis, scf_results, rsp_results
        )
        hessian_values = None
        if require_hessian:
            hessian_values = self.hessian_function(self.drivers[2], molecule, basis)

        energy = float(np.asarray(energy_values).reshape(-1)[0])
        physical_gradient = np.asarray(gradient_values[0], dtype=float)
        physical_hessian = None
        if hessian_values is not None:
            physical_hessian = np.asarray(hessian_values[0], dtype=float).reshape(
                physical_gradient.size, physical_gradient.size
            )

        point = InterpolationDatapoint(self.z_matrix, atom_labels=molecule.get_labels())
        point.update_settings(self.interpolation_settings)
        point.cartesian_coordinates = np.asarray(
            molecule.get_coordinates_in_bohr(), dtype=float
        )
        point.eq_bond_lengths = self.eq_bond_lengths.copy()
        point.imp_int_coordinates = {
            "bonds": [], "angles": [], "dihedrals": [], "impropers": []
        }
        if self.interpolation_settings.get("use_mass_weight", False):
            inverse_sqrt_masses = 1.0 / np.sqrt(
                np.repeat(np.asarray(molecule.get_masses(), dtype=float), 3)
            )
            point.inv_sqrt_masses = inverse_sqrt_masses
            point.gradient = (
                inverse_sqrt_masses * physical_gradient.reshape(-1)
            ).reshape(physical_gradient.shape)
            if physical_hessian is not None:
                point.hessian = (
                    inverse_sqrt_masses[:, None]
                    * physical_hessian
                    * inverse_sqrt_masses[None, :]
                )
        else:
            point.gradient = physical_gradient.copy()
            point.hessian = None if physical_hessian is None else physical_hessian.copy()
        point.energy = energy
        if point.hessian is None:
            point.transform_gradient()
        else:
            point.transform_gradient_and_hessian()

        measured_phase = ()
        rotor = self.rotor_by_subset.get(request.subset_id)
        if rotor is not None:
            from .phase import estimate_collective_phase

            measured_phase = (
                estimate_collective_phase(
                    point.internal_coordinates_values,
                    self.anchor_internal_coordinates,
                    rotor,
                ),
            )
        target = np.asarray(request.phase_signature, dtype=float)
        measured = np.asarray(measured_phase, dtype=float)
        constraint_errors = ()
        if target.size and measured.size:
            from .phase import principal_periodic_delta

            constraint_errors = tuple(
                float(value)
                for value in np.abs(principal_periodic_delta(measured - target))
            )

        converged = bool(
            np.isfinite(energy)
            and np.all(np.isfinite(physical_gradient))
            and (physical_hessian is None or np.all(np.isfinite(physical_hessian)))
        )
        provenance = json.dumps(
            {
                "provider": "GroupedCalculationProvider",
                "physical_cartesian_derivatives": True,
                "mass_weighted_internal_transform": bool(
                    self.interpolation_settings.get("use_mass_weight", False)
                ),
                "electronic_structure_dispatch": "IMForceFieldGenerator",
            },
            sort_keys=True,
        )
        return GroupedCalculationResult(
            request_id=request.request_id,
            converged=converged,
            cartesian_coordinates=tuple(
                tuple(float(value) for value in row)
                for row in point.cartesian_coordinates
            ),
            energy=energy,
            cartesian_gradient=tuple(
                tuple(float(value) for value in row) for row in physical_gradient
            ),
            cartesian_hessian=None
            if physical_hessian is None
            else tuple(tuple(float(value) for value in row) for row in physical_hessian),
            internal_coordinates=tuple(
                float(value) for value in point.internal_coordinates_values
            ),
            internal_gradient=tuple(float(value) for value in point.internal_gradient),
            internal_hessian=None
            if point.internal_hessian is None
            else tuple(
                tuple(float(value) for value in row) for row in point.internal_hessian
            ),
            measured_phase_signature=measured_phase,
            constraint_errors=constraint_errors,
            coordinate_fingerprint=self.coordinate_fingerprint,
            source_point_label=request.request_id,
            provenance=provenance,
        )
