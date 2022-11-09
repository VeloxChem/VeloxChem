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
import h5py

from .veloxchemlib import AOFockMatrix
from .veloxchemlib import fockmat
from .errorhandler import assert_msg_critical


def _AOFockMatrix_write_hdf5(self, fname):
    """
    Writes AOFockMatrix to hdf5 file.

    :param fname:
        The name of the hdf5 file.
    """

    focktype = {
        fockmat.restjk: 'restjk',
        fockmat.restjkx: 'restjkx',
        fockmat.restj: 'restj',
        fockmat.restk: 'restk',
        fockmat.restkx: 'restkx',
        fockmat.rgenjk: 'rgenjk',
        fockmat.rgenjkx: 'rgenjkx',
        fockmat.rgenj: 'rgenj',
        fockmat.rgenk: 'rgenk',
        fockmat.rgenkx: 'rgenkx',
        fockmat.unrestjk: 'unrestjk',
        fockmat.unrestjkx: 'unrestjkx',
        fockmat.unrestj: 'unrestj',
    }

    hf = h5py.File(fname, 'w')

    factors = []

    matrix_count = 0

    for fock_id in range(self.number_of_fock_matrices()):

        factors.append(self.get_scale_factor(fock_id, 'alpha'))

        fock_type_str = focktype[self.get_fock_type(fock_id,
                                                    'alpha')] + '.alpha'
        density_id = self.get_density_identifier(fock_id)
        name = f'{matrix_count}_{fock_type_str}_{density_id}'
        array = self.alpha_to_numpy(fock_id)
        hf.create_dataset(name, data=array, compression='gzip')

        matrix_count += 1

        if not self.is_closed_shell():

            factors.append(self.get_scale_factor(fock_id, 'beta'))

            fock_type_str = focktype[self.get_fock_type(fock_id,
                                                        'beta')] + '.beta'
            density_id = self.get_density_identifier(fock_id)
            name = f'{matrix_count}_{fock_type_str}_{density_id}'
            array = self.beta_to_numpy(fock_id)
            hf.create_dataset(name, data=array, compression='gzip')

            matrix_count += 1

    hf.create_dataset('factors', data=factors, compression='gzip')

    hf.close()


@staticmethod
def _AOFockMatrix_read_hdf5(fname):
    """
    Reads AOFockMatrix from hdf5 file.

    :param fname:
        The name of the hdf5 file.

    :return:
        The AOFockMatrix.
    """

    focktype = {
        'restjk': fockmat.restjk,
        'restjkx': fockmat.restjkx,
        'restj': fockmat.restj,
        'restk': fockmat.restk,
        'restkx': fockmat.restkx,
        'rgenjk': fockmat.rgenjk,
        'rgenjkx': fockmat.rgenjkx,
        'rgenj': fockmat.rgenj,
        'rgenk': fockmat.rgenk,
        'rgenkx': fockmat.rgenkx,
        'unrestjk': fockmat.unrestjk,
        'unrestjkx': fockmat.unrestjkx,
        'unrestj': fockmat.unrestj,
    }

    hf = h5py.File(fname, 'r')

    focks = []
    types = []
    factors = list(hf.get('factors'))
    density_ids = []

    matrix_id_and_keys = []
    for key in list(hf.keys()):
        if key != 'factors':
            matrix_id_and_keys.append((int(key.split('_')[0]), key))
    matrix_id_and_keys = sorted(matrix_id_and_keys)

    for i, key in matrix_id_and_keys:
        type_str, id_str = key.split('_')[1:]
        focks.append(np.array(hf.get(key)))
        types.append(focktype[type_str.split('.')[0]])
        density_ids.append(int(id_str))

    hf.close()

    for ftype in set(types):
        assert_msg_critical(ftype in list(focktype.values()),
                            'AOFockMatrix.read_hdf5: invalid Fock types')

    return AOFockMatrix(focks, types, factors, density_ids)


AOFockMatrix.write_hdf5 = _AOFockMatrix_write_hdf5
AOFockMatrix.read_hdf5 = _AOFockMatrix_read_hdf5
