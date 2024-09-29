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

from .veloxchemlib import NuclearPotentialDriver
from .veloxchemlib import ElectricDipoleMomentDriver


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
    return -1.0 * npot_mat.full_matrix().to_numpy()


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
        -1.0 * dip_mats.matrix('X').full_matrix().to_numpy(),
        -1.0 * dip_mats.matrix('Y').full_matrix().to_numpy(),
        -1.0 * dip_mats.matrix('Z').full_matrix().to_numpy(),
    ])
