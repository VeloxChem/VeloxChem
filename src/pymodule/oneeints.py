#
#                              VELOXCHEM
#         ----------------------------------------------------
#                     An Electronic Structure Code
#
#  Copyright © 2018-2024 by VeloxChem developers. All rights reserved.
#
#  SPDX-License-Identifier: LGPL-3.0-or-later
#
#  This file is part of VeloxChem.
#
#  VeloxChem is free software: you can redistribute it and/or modify it under
#  the terms of the GNU Lesser General Public License as published by the Free
#  Software Foundation, either version 3 of the License, or (at your option)
#  any later version.
#
#  VeloxChem is distributed in the hope that it will be useful, but WITHOUT
#  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
#  License for more details.
#
#  You should have received a copy of the GNU Lesser General Public License
#  along with VeloxChem. If not, see <https://www.gnu.org/licenses/>.

import numpy as np

from .veloxchemlib import OverlapDriver
from .veloxchemlib import KineticEnergyDriver
from .veloxchemlib import NuclearPotentialDriver
from .veloxchemlib import NuclearPotentialGeom100Driver
from .veloxchemlib import ElectricDipoleMomentDriver
from .veloxchemlib import NuclearPotentialGeom200Driver
from .veloxchemlib import NuclearPotentialGeom101Driver
from .veloxchemlib import compute_linear_momentum_integrals
from .veloxchemlib import compute_angular_momentum_integrals
from .veloxchemlib import compute_electric_field_integrals
from .veloxchemlib import compute_electric_field_values
from .veloxchemlib import compute_electric_field_potential_gradient
from .veloxchemlib import compute_electric_field_fock_gradient
from .veloxchemlib import compute_electric_field_potential_gradient_for_mm
from .veloxchemlib import compute_electric_field_potential_hessian
from .matrices import Matrices


def compute_overlap_integrals(molecule, basis):
    """
    Computes overlap integrals.

    :param molecule:
        The molecule.
    :param basis:
        The molecular basis set.

    :return:
        The overlap integral matrix.
    """

    ovl_drv = OverlapDriver()
    ovl_mat = ovl_drv.compute(molecule, basis)

    return ovl_mat.to_numpy()


def compute_kinetic_energy_integrals(molecule, basis):
    """
    Computes kinetic energy integrals.

    :param molecule:
        The molecule.
    :param basis:
        The molecular basis set.

    :return:
        The kinetic energy integral matrix.
    """

    kin_drv = KineticEnergyDriver()
    kin_mat = kin_drv.compute(molecule, basis)

    return kin_mat.to_numpy()


def compute_nuclear_potential_integrals(molecule,
                                        basis,
                                        charges=None,
                                        coordinates=None):
    """
    Computes nuclear potential integrals.

    :param molecule:
        The molecule.
    :param basis:
        The molecular basis set.

    :return:
        The nuclear potential integral matrix.
    """

    npot_drv = NuclearPotentialDriver()

    if charges is None and coordinates is None:
        npot_mat = npot_drv.compute(molecule, basis)
    else:
        npot_mat = npot_drv.compute(molecule, basis, charges, coordinates)

    # Note: factor -1.0 for electron charge
    return -1.0 * npot_mat.to_numpy()


def compute_electric_dipole_integrals(molecule, basis, origin=(0.0, 0.0, 0.0)):
    """
    Computes electric dipole integrals.

    :param molecule:
        The molecule.
    :param basis:
        The molecular basis set.

    :return:
        A tuple containint the electric dipole integral matrices.
    """

    dip_drv = ElectricDipoleMomentDriver()
    dip_mats = dip_drv.compute(molecule, basis, list(origin))

    # Note: factor -1.0 for electron charge
    return tuple([
        -1.0 * dip_mats.matrix_to_numpy('X'),
        -1.0 * dip_mats.matrix_to_numpy('Y'),
        -1.0 * dip_mats.matrix_to_numpy('Z'),
    ])


def compute_nuclear_potential_gradient_bfs(molecule, basis, charges,
                                           coordinates, D):
    """
    Computes nuclear potential integrals contribution from point charges to
    molecular gradient.

    :param molecule:
        The molecule.
    :param basis:
        The molecular basis set.
    :param charges:
        The point charges.
    :param coordinates:
        The coordinates of point charges.

    :return:
        The nuclear potential integral contribution to molecular gradient.
    """

    natoms = molecule.number_of_atoms()

    grad = np.zeros((natoms, 3))

    npot_grad_100_drv = NuclearPotentialGeom100Driver()

    for iatom in range(natoms):
        gmats_100 = npot_grad_100_drv.compute(molecule, basis, iatom,
                                              coordinates, charges)

        for i, label in enumerate(['X', 'Y', 'Z']):
            gmat_100 = gmats_100.matrix_to_numpy(label)

            grad[iatom, i] += np.sum((gmat_100 + gmat_100.T) * D)

        gmats_100 = Matrices()

    # Note: factor -1.0 for electron charge
    return -1.0 * grad


def compute_electrostatic_potential_hessian(molecule, basis, mm_charges,
                                            mm_coordinates, density,
                                            qm_atom_index_i, qm_atom_index_j):

    hess = np.zeros((3, 3))

    i, j = qm_atom_index_i, qm_atom_index_j

    if i == j:
        npot_hess_200_drv = NuclearPotentialGeom200Driver()

        hmats_200 = npot_hess_200_drv.compute(molecule, basis, i,
                                              mm_coordinates, mm_charges)

        for x, label_x in enumerate('XYZ'):
            for y, label_y in enumerate('XYZ'):
                npot_label = label_x + label_y if x <= y else label_y + label_x
                npot_200_iixy = hmats_200.matrix_to_numpy(npot_label)
                hess[x, y] += 2.0 * (np.sum(density *
                                            (npot_200_iixy + npot_200_iixy.T)))

        hmats_200 = Matrices()

    npot_hess_101_drv = NuclearPotentialGeom101Driver()

    hmats_101 = npot_hess_101_drv.compute(molecule, basis, i, j, mm_coordinates,
                                          mm_charges)

    for x, label_x in enumerate('XYZ'):
        for y, label_y in enumerate('XYZ'):
            npot_xy_label = f'{label_x}_{label_y}'
            npot_101_ijxy = hmats_101.matrix_to_numpy(npot_xy_label)
            hess[x, y] += 2.0 * (np.sum(density *
                                        (npot_101_ijxy + npot_101_ijxy.T)))

    hmats_101 = Matrices()

    # Note: factor -1.0 for electron charge
    return -1.0 * hess


def compute_electrostatic_integrals_gradient(molecule, basis, mm_charges,
                                             mm_coordinates, qm_atom_index_i):

    naos = basis.get_dimensions_of_basis()

    ints_grad = np.zeros((3, naos, naos))

    i = qm_atom_index_i

    npot_grad_100_drv = NuclearPotentialGeom100Driver()

    gmats_100 = npot_grad_100_drv.compute(molecule, basis, i, mm_coordinates,
                                          mm_charges)

    for x, label in enumerate(['X', 'Y', 'Z']):
        gmat_100 = gmats_100.matrix_to_numpy(label)
        # Note: factor -1.0 for electron charge
        ints_grad[x] -= gmat_100 + gmat_100.T

    gmats_100 = Matrices()

    return ints_grad
