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

from .veloxchemlib import DenseMatrix


def get_basis_function_indices_of_atoms(molecule, basis):
    """
    Gets AO indices of atoms.

    :param molecule:
        The molecule.
    :param basis:
        The AO basis set.

    :return:
        A list of list containing AO indices of atoms.
    """

    natoms = molecule.number_of_atoms()
    aoinds_atoms = [[] for atomidx in range(natoms)]

    max_angl = basis.max_angular_momentum()

    aoidx = 0

    for angl in range(max_angl + 1):
        indices = [[] for atomidx in range(natoms)]

        for s in range(-angl, angl + 1):

            for atomidx in range(natoms):
                indices[atomidx].append([])
                nao = basis.number_of_basis_functions([atomidx], angl)

                for i in range(nao):
                    indices[atomidx][-1].append(aoidx)
                    aoidx += 1

        for atomidx in range(natoms):

            reordered_indices = []

            if len(indices[atomidx]) == 3:
                #  1   2
                # -1   0
                #  0   1
                reordered_indices.append(indices[atomidx][2])
                reordered_indices.append(indices[atomidx][0])
                reordered_indices.append(indices[atomidx][1])
            else:
                for s in range(len(indices[atomidx])):
                    reordered_indices.append(indices[atomidx][s])

            aoinds_atoms[atomidx] += list(
                np.array(reordered_indices).T.reshape(-1))

    flat_inds = []
    for atomidx in range(natoms):
        flat_inds += aoinds_atoms[atomidx]

    return flat_inds


def ao_matrix_to_dalton(array, basis, molecule):

    bf_indices = get_basis_function_indices_of_atoms(molecule, basis)

    if isinstance(array, DenseMatrix):
        return DenseMatrix(array.to_numpy()[bf_indices, :][:, bf_indices])
    else:
        return array[bf_indices, :][:, bf_indices]


def ao_matrix_to_veloxchem(array, basis, molecule):

    bf_indices = get_basis_function_indices_of_atoms(molecule, basis)

    reverse_indices = [(x, i) for i, x in enumerate(bf_indices)]
    reverse_indices = [x[1] for x in sorted(reverse_indices)]

    if isinstance(array, DenseMatrix):
        return DenseMatrix(
            array.to_numpy()[reverse_indices, :][:, reverse_indices])
    else:
        return array[reverse_indices, :][:, reverse_indices]
