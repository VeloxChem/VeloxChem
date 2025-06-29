#
#                                   VELOXCHEM
#              ----------------------------------------------------
#                          An Electronic Structure Code
#
#  SPDX-License-Identifier: BSD-3-Clause
#
#  Copyright 2018-2025 VeloxChem developers
#
#  Redistribution and use in source and binary forms, with or without modification,
#  are permitted provided that the following conditions are met:
#
#  1. Redistributions of source code must retain the above copyright notice, this
#     list of conditions and the following disclaimer.
#  2. Redistributions in binary form must reproduce the above copyright notice,
#     this list of conditions and the following disclaimer in the documentation
#     and/or other materials provided with the distribution.
#  3. Neither the name of the copyright holder nor the names of its contributors
#     may be used to endorse or promote products derived from this software without
#     specific prior written permission.
#
#  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
#  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
#  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
#  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
#  FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
#  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
#  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
#  HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
#  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
#  OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

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
