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
from .veloxchemlib import compute_linear_momentum_integrals
from .veloxchemlib import compute_angular_momentum_integrals
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
