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


class Seminario:
    """
    Class for calculating force constants using the Seminario method.
    """

    def __init__(self, hessian_matrix, coordinates):
        """
        Initializes the Seminaro object.
        """

        self.hessian_matrix = np.copy(hessian_matrix)
        self.coordinates = coordinates

    def get_reference(self):
        """
        Get string for reference paper.
        """

        return 'J. M. Seminario, Int. J. Quantum Chem. 60, 1271-1277 (1996)'

    def calculate_interatomic_k(self, atom_index_a, atom_index_b, vector_u):
        """
        Calculate the sum of eigenvectors multiplied by their respective
        eigenvalue for a given 3x3 submatrix of the Hessian matrix.
        """

        a, b = atom_index_a, atom_index_b

        submatrix = -1.0 * self.hessian_matrix[a * 3:(a + 1) * 3,
                                               b * 3:(b + 1) * 3]

        eigvals, eigvecs = np.linalg.eig(submatrix)
        eigvecs_T = np.transpose(eigvecs)

        # TODO: check for pairwise unstable interaction
        # - if one or more eigvals are negative
        # - if two of the three eigvals or eigvecs are complex

        return sum([
            eigvals[i] * abs(np.dot(vector_u, eigvecs_T[i])) for i in range(3)
        ]).real

    def unit_vector(self, atom_index_a, atom_index_b):
        """
        Calculate the unit vector defined by two atoms.
        """

        vec = (self.coordinates[atom_index_a] - self.coordinates[atom_index_b])
        return vec / np.linalg.norm(vec)

    def unit_normal_vector(self, atom_index_a, atom_index_b, atom_index_c):
        """
        Calculate the unit normal vector for the plane defined by three atoms.
        """

        vec_ab = self.unit_vector(atom_index_a, atom_index_b)
        vec_cb = self.unit_vector(atom_index_c, atom_index_b)

        cross_product = np.cross(vec_cb, vec_ab)
        return cross_product / np.linalg.norm(cross_product)

    def squared_distance(self, atom_index_a, atom_index_b):
        """
        Calculate the squared distance between two atoms.
        """

        return np.sum((self.coordinates[atom_index_a] -
                       self.coordinates[atom_index_b])**2)

    def calculate_bond_force_constant(self, atom_index_a, atom_index_b):
        """
        Calculate the bond force constant for a given pair of atom indices.
        """

        a, b = atom_index_a, atom_index_b

        vec = self.unit_vector(a, b)
        k_bond_1 = self.calculate_interatomic_k(a, b, vec)
        k_bond_1 = max(0.0, k_bond_1)

        vec = self.unit_vector(b, a)
        k_bond_2 = self.calculate_interatomic_k(b, a, vec)
        k_bond_2 = max(0.0, k_bond_2)

        return 0.5 * (k_bond_1 + k_bond_2)

    def calculate_angle_force_constant(self, atom_index_a, atom_index_b,
                                       atom_index_c):
        """
        Calculate the angle force constant for a given set of three atom
        indices.
        """

        a, b, c = atom_index_a, atom_index_b, atom_index_c

        normal_vector = self.unit_normal_vector(a, b, c)
        rab2 = self.squared_distance(a, b)
        rcb2 = self.squared_distance(c, b)

        perp_vector_a = np.cross(normal_vector, self.unit_vector(a, b))
        perp_vector_c = np.cross(self.unit_vector(c, b), normal_vector)

        denominator_a = rab2 * self.calculate_interatomic_k(a, b, perp_vector_a)
        denominator_c = rcb2 * self.calculate_interatomic_k(c, b, perp_vector_c)

        k_angle_1 = 1.0 / (1.0 / denominator_a + 1.0 / denominator_c)
        k_angle_1 = max(0.0, k_angle_1)

        perp_vector_a = np.cross(normal_vector, self.unit_vector(b, a))
        perp_vector_c = np.cross(self.unit_vector(b, c), normal_vector)

        denominator_a = rab2 * self.calculate_interatomic_k(b, a, perp_vector_a)
        denominator_c = rcb2 * self.calculate_interatomic_k(b, c, perp_vector_c)

        k_angle_2 = 1.0 / (1.0 / denominator_a + 1.0 / denominator_c)
        k_angle_2 = max(0.0, k_angle_2)

        return 0.5 * (k_angle_1 + k_angle_2)
