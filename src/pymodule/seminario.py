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
        self.hessian_matrix = hessian_matrix
        self.coordinates = coordinates

    def calculate_eigenvalue_sum(self, atom_index_a, atom_index_b, direction_vector):
        """
        Calculate the sum of eigenvectors multiplied by their respective eigenvalue for a given 3x3 submatrix of the Hessian matrix.
        """
        submatrix = -self.hessian_matrix[3 * atom_index_a : 3 * (atom_index_a + 1), 3 * atom_index_b : 3 * (atom_index_b + 1)]
        eigenvalues, eigenvectors = np.linalg.eig(submatrix)
        eigenvectors = np.transpose(eigenvectors)

        eigenvalue_sum = sum(eigenvalues[i] * abs(np.dot(direction_vector, eigenvectors[i])) for i in range(3))
        return eigenvalue_sum

    def unit_vector(self, atom_index_a, atom_index_b):
        displacement = self.coordinates[atom_index_a] - self.coordinates[atom_index_b]
        return displacement / np.linalg.norm(displacement)
    
    def normal_unit_vector(self, atom_index_a, atom_index_b, atom_index_c):
        """
        Calculate the unit normal vector for the plane defined by three atoms.
        """
        vector_ab = self.unit_vector(atom_index_a, atom_index_b)
        vector_cb = self.unit_vector(atom_index_c, atom_index_b)
        cross_product = np.cross(vector_cb, vector_ab)
        return cross_product / np.linalg.norm(cross_product)

    def squared_distance(self, atom_index_a, atom_index_b):
        """
        Calculate the squared distance between two atoms.
        """
        return np.sum(np.square(self.coordinates[atom_index_a] - self.coordinates[atom_index_b]))

    def calculate_angle_force_constant(self, atom_index_a, atom_index_b, atom_index_c):
        """
        Calculate the angle force constant for a given set of three atom indices.
        """
        def calculate_partial_force_constant(i, j, k):
            normal_vector = self.normal_unit_vector(i, j, k)
            perp_vector_i = np.cross(normal_vector, self.unit_vector(i, j))
            perp_vector_k = np.cross(self.unit_vector(k, j), normal_vector)

            denominator1 = self.squared_distance(i, j) * self.calculate_eigenvalue_sum(i, j, perp_vector_i)
            denominator2 = self.squared_distance(k, j) * self.calculate_eigenvalue_sum(k, j, perp_vector_k)

            return 1 / (1 / denominator1 + 1 / denominator2).real
        
        return max(0, 0.5 * (calculate_partial_force_constant(atom_index_a, atom_index_b, atom_index_c) + calculate_partial_force_constant(atom_index_c, atom_index_b, atom_index_a)))
    
    def calculate_bond_force_constant(self, atom_index_a, atom_index_b):
        """
        Calculate the bond force constant for a given pair of atom indices.
        """
        def calculate_partial_force_constant(i, j):
            vector_ab = self.unit_vector(i, j)
            return self.calculate_eigenvalue_sum(i, j, vector_ab).real
        
        return max(0, 0.5 * (calculate_partial_force_constant(atom_index_a, atom_index_b) + calculate_partial_force_constant(atom_index_b, atom_index_a)))

    