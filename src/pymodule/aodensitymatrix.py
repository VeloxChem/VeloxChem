#
#                           VELOXCHEM 1.0-RC2
#         ----------------------------------------------------
#                     An Electronic Structure Code
#
#  Copyright © 2018-2021 by VeloxChem developers. All rights reserved.
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
import h5py

from .veloxchemlib import AODensityMatrix
from .veloxchemlib import denmat
from .errorhandler import assert_msg_critical


def _AODensityMatrix_write_hdf5(self, fname):
    """
    Writes AODensityMatrix to hdf5 file.

    :param fname:
        The name of the hdf5 file.
    """

    hf = h5py.File(fname, 'w')

    matrix_count = 0

    density_type = self.get_density_type()

    for density_id in range(self.number_of_density_matrices()):

        if density_type == denmat.rest:

            name = f'{matrix_count}_rest.alpha_{density_id}'
            array = self.alpha_to_numpy(density_id)
            hf.create_dataset(name, data=array, compression='gzip')
            matrix_count += 1

        elif density_type == denmat.unrest:

            name = f'{matrix_count}_unrest.alpha_{density_id}'
            array = self.alpha_to_numpy(density_id)
            hf.create_dataset(name, data=array, compression='gzip')
            matrix_count += 1

            name = f'{matrix_count}_unrest.beta_{density_id}'
            array = self.beta_to_numpy(density_id)
            hf.create_dataset(name, data=array, compression='gzip')
            matrix_count += 1

    hf.close()


@staticmethod
def _AODensityMatrix_read_hdf5(fname):
    """
    Reads AODensityMatrix from hdf5 file.

    :param fname:
        The name of the hdf5 file.

    :return:
        The AODensityMatrix.
    """

    dentype = {
        'rest.alpha': denmat.rest,
        'unrest.alpha': denmat.unrest,
        'unrest.beta': denmat.unrest,
    }

    hf = h5py.File(fname, 'r')

    matrix_id_and_keys = sorted([
        (int(key.split('_')[0]), key) for key in list(hf.keys())
    ])

    dens = [np.array(hf.get(key)) for matrix_id, key in matrix_id_and_keys]

    types = [
        dentype[key.split('_')[1]] for matrix_id, key in matrix_id_and_keys
    ]

    hf.close()

    assert_msg_critical(
        len(set(types)) == 1,
        'AODensityMatrix.read_hdf5: inconsistent density type')

    assert_msg_critical(types[0] in [denmat.rest, denmat.unrest],
                        'AODensityMatrix.read_hdf5: invalid density type')

    return AODensityMatrix(dens, types[0])


AODensityMatrix.write_hdf5 = _AODensityMatrix_write_hdf5
AODensityMatrix.read_hdf5 = _AODensityMatrix_read_hdf5
