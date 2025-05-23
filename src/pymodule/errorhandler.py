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

from mpi4py import MPI
from sys import stderr
import numpy as np
import math


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
        if not condition:
            stderr.write(' **** Critical Error (process {}) **** {}\n'.format(
                MPI.COMM_WORLD.Get_rank(), msg))
            stderr.flush()
            MPI.COMM_WORLD.Abort()


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


def safe_arcsin(val):
    """
    Safely uses math.asin and avoids the math domain error.

    :param val:
        The sine value.

    :return:
        The angle in radian.
    """

    if abs(val) > 1.0:
        # avoid math domain error
        assert_msg_critical(
            abs(abs(val) - 1.0) < 1.0e-12, 'arcsin: Invalid sine value')
        sin_phi = 1.0 if val > 1.0 else -1.0
    else:
        sin_phi = val

    return math.asin(sin_phi)


def safe_solve(mat, b):
    """
    Safely solve Ax=b.
    """

    try:
        sol = np.linalg.solve(mat, b)
    except np.linalg.LinAlgError:
        sol = np.dot(np.linalg.pinv(mat), b)

    return sol
