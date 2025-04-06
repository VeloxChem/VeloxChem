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
