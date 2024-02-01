#
#                           VELOXCHEM 1.0-RC3
#         ----------------------------------------------------
#                     An Electronic Structure Code
#
#  Copyright © 2018-2022 by VeloxChem developers. All rights reserved.
#  Contact: https://veloxchem.org/contact
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
