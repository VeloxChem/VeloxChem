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

from .veloxchemlib import CubicGrid
from .errorhandler import assert_msg_critical


@staticmethod
def _CubicGrid_read_cube(cube_file):
    """
    Creates cubic grid from a cube file.

    :param cube_file:
        The name of the cube file.

    :return:
        The cubic grid.
    """

    f_cube = open(str(cube_file), 'r')

    f_cube.readline()
    f_cube.readline()

    nsteps = []
    stepsize = []

    content = f_cube.readline().split()
    is_mo = int(content[0]) < 0
    natoms = abs(int(content[0]))
    origin = [float(x) for x in content[1:4]]

    for d in range(3):
        content = f_cube.readline().split()
        nsteps.append(int(content[0]))
        stepsize.append(float(content[d + 1]))

    for a in range(natoms):
        f_cube.readline()

    if is_mo:
        f_cube.readline()

    values = []

    while True:
        line = f_cube.readline()
        if not line:
            break
        values += [float(x) for x in line.split()]

    f_cube.close()

    grid = CubicGrid(origin, stepsize, nsteps)
    grid.set_values(values)

    return grid


def _CubicGrid_compare(self, other_cubic_grid):
    """
    Compares self with another cubic grid.

    :param other_cubic_grid:
        The other cubic grid.

    :return:
        The maximum deviation.
    """

    vals_1 = self.values_to_numpy()
    vals_2 = other_cubic_grid.values_to_numpy()

    assert_msg_critical(
        vals_1.size == vals_2.size,
        'CubicGrid.compare: Inconsistent number of grid points')

    max_d_1 = np.max(np.abs(vals_1 - vals_2))
    max_d_2 = np.max(np.abs(vals_1 + vals_2))

    return min(max_d_1, max_d_2)


CubicGrid.read_cube = _CubicGrid_read_cube
CubicGrid.compare = _CubicGrid_compare
