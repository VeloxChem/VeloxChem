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
