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

from mpi4py import MPI
import math

from .veloxchemlib import assert_msg_critical as vlx_assert


def assert_msg_critical(condition, msg=''):
    """
    Asserts that the condition is true. Otherwise terminates the program with
    an error message.

    :param condition:
        The condition.
    :param msg:
        The error message.
    """
    if __debug__ and MPI.COMM_WORLD.Get_size() == 1:
        assert condition, msg
    else:
        vlx_assert(condition, msg)


def safe_arccos(val):
    """
    Safely uses math.acos and avoids the math domain error.

    :param val:
        The cosine value.

    :return:
        The angle in radian.
    """

    if abs(val) > 1.0:
        # avoid math domain error
        assert_msg_critical(
            abs(abs(val) - 1.0) < 1.0e-12, 'arccos: Invalid cosine value')
        cos_phi = 1.0 if val > 1.0 else -1.0
    else:
        cos_phi = val

    return math.acos(cos_phi)
