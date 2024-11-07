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

from .errorhandler import assert_msg_critical


class Inversion:

    def __init__(self):

        self._axis = None

    def __str__(self):

        return 'Ci'

    def get_matrix(self):

        return -np.identity(3)


class Rotation:

    def __init__(self, axis, order=1, power=1):

        self._order = order
        self._power = power

        # normalize axis
        norm = np.linalg.norm(axis)
        assert_msg_critical(norm > 1e-8, 'Reflection: Invalid axis')
        self._axis = np.array(axis) / norm

    def __str__(self):

        if self._power == 1:
            return f'C{self._order}'
        else:
            return f'C{self._order}^{self._power}'

    def get_matrix(self):

        return rotation_matrix(self._axis,
                               (2.0 * np.pi / self._order) * self._power)


class Reflection:

    def __init__(self, axis):

        # normalize axis
        norm = np.linalg.norm(axis)
        assert_msg_critical(norm > 1e-8, 'Reflection: Invalid axis')
        self._axis = np.array(axis) / norm

    def __str__(self):

        return 'Cs'

    def get_matrix(self):

        return np.identity(3) - 2 * np.outer(self._axis, self._axis)


class ImproperRotation:

    def __init__(self, axis, order=1, power=1):

        self._order = order
        self._power = power

        # normalize axis
        norm = np.linalg.norm(axis)
        assert_msg_critical(norm > 1e-8, 'Reflection: Invalid axis')
        self._axis = np.array(axis) / norm

    def __str__(self):

        if self._power == 1:
            return f'S{self._order}'
        else:
            return f'S{self._order}^{self._power}'

    def get_matrix(self):

        # Build the rotation matrix
        rot_matrix = rotation_matrix(self._axis,
                                     (2 * np.pi / self._order) * self._power)

        # Build the reflexion matrix
        reflection = Reflection(self._axis)
        refl_matrix = reflection.get_matrix()

        return np.dot(rot_matrix, refl_matrix.T)


def rotation_matrix(axis, angle):
    """
    Build a rotation matrix to rotate of a given angle around a vector.

    :param axis:
        The normalized axis or vector around which the rotation is effectuated.
    :param angle:
        The angle to be rotated in radians.

    :return:
        The rotation matrix.
    """

    norm = np.linalg.norm(axis)
    assert_msg_critical(norm > 1e-8,
                        'SymmetryOperation.rotation_matrix: Invalid axis')

    # normalize axis
    u_vec = np.array(axis) / norm
    sin_theta = np.sin(angle)
    cos_theta = np.cos(angle)
    cos_term = 1 - cos_theta

    rot_matrix = [
        [
            u_vec[0] * u_vec[0] * cos_term + cos_theta,
            u_vec[0] * u_vec[1] * cos_term - u_vec[2] * sin_theta,
            u_vec[0] * u_vec[2] * cos_term + u_vec[1] * sin_theta,
        ],
        [
            u_vec[1] * u_vec[0] * cos_term + u_vec[2] * sin_theta,
            u_vec[1] * u_vec[1] * cos_term + cos_theta,
            u_vec[1] * u_vec[2] * cos_term - u_vec[0] * sin_theta,
        ],
        [
            u_vec[2] * u_vec[0] * cos_term - u_vec[1] * sin_theta,
            u_vec[2] * u_vec[1] * cos_term + u_vec[0] * sin_theta,
            u_vec[2] * u_vec[2] * cos_term + cos_theta,
        ],
    ]

    return np.array(rot_matrix)
