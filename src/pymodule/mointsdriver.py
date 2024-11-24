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

from mpi4py import MPI
import numpy as np
import math
import sys

from .veloxchemlib import T4CScreener
from .veloxchemlib import mpi_master
from .fockdriver import FockDriver
from .outputstream import OutputStream
from .errorhandler import assert_msg_critical


class MOIntegralsDriver:
    """
    Implements MO integrals driver.

    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.

    Instance variable
        - eri_thresh: The electron repulsion integrals screening threshold.
        - num_matrices: Number of Fock matrices to be computed.
        - batch_size: Batch size for computing Fock matrices.
        - comm: The MPI communicator.
        - rank: The MPI rank.
        - nodes: Number of MPI processes.
        - ostream: The output stream.
    """

    def __init__(self, comm=None, ostream=None):
        """
        Initializes MO integrals driver  to default setup.
        """

        if comm is None:
            comm = MPI.COMM_WORLD

        if ostream is None:
            if comm.Get_rank() == mpi_master():
                ostream = OutputStream(sys.stdout)
            else:
                ostream = OutputStream(None)

        # screening scheme
        self.eri_thresh = 1.0e-12

        # Fock matrices
        self.num_matrices = 0
        self.batch_size = 3000

        # mpi information
        self.comm = comm
        self.rank = self.comm.Get_rank()
        self.nodes = self.comm.Get_size()

        # output stream
        self.ostream = ostream

    def compute_in_memory(self,
                          molecule,
                          basis,
                          mol_orbs,
                          moints_name,
                          moints_spin='aaaa'):
        """
        Performs in-memory MO integrals calculation for a molecule and a basis
        set.

        :param molecule:
            The molecule.
        :param basis:
            The AO basis set.
        :param mol_orbs:
            The molecular orbitals.
        :param moints_name:
            The name of MO integrals to be calculated, such as 'chem_oovv'.
        :param moints_spin:
            The spin of MO integrals to be calculated, such as 'aabb'.

        :return:
            The computed MO integrals.
        """

        mo_a = mol_orbs.alpha_to_numpy()
        mo_b = mol_orbs.beta_to_numpy()

        naos = mo_a.shape[0]

        nocc_a = molecule.number_of_alpha_electrons()
        nocc_b = molecule.number_of_beta_electrons()

        mo_occ_a = mo_a[:, :nocc_a].copy()
        mo_vir_a = mo_a[:, nocc_a:].copy()

        mo_occ_b = mo_b[:, :nocc_b].copy()
        mo_vir_b = mo_b[:, nocc_b:].copy()

        # get the type of MO integral, such as 'oovv'

        chem_notation = moints_name.lower().startswith('chem_')
        phys_notation = moints_name.lower().startswith('phys_')

        err_msg = 'MOIntegralsDriver.compute_in_memory: invalid moints_name'
        assert_msg_critical(chem_notation or phys_notation, err_msg)

        if chem_notation:
            moints_type = moints_name.lower().replace('chem_', '')
        elif phys_notation:
            moints_type = moints_name.lower().replace('phys_', '')

        assert_msg_critical(len(moints_type) == 4, err_msg)
        for x in moints_type:
            assert_msg_critical(x == 'o' or x == 'v', err_msg)

        # MO coefficients

        mo_coef_dict = {
            'oa': mo_occ_a,
            'va': mo_vir_a,
            'ob': mo_occ_b,
            'vb': mo_vir_b,
        }
        if chem_notation:
            mo_coefs = [
                mo_coef_dict[t + s] for t, s in zip(moints_type, moints_spin)
            ]
        elif phys_notation:
            mo_coefs = [
                mo_coef_dict[moints_type[0] + moints_spin[0]],
                mo_coef_dict[moints_type[2] + moints_spin[2]],
                mo_coef_dict[moints_type[1] + moints_spin[1]],
                mo_coef_dict[moints_type[3] + moints_spin[3]],
            ]

        # AO integrals

        t4c_drv = T4CScreener()
        t4c_drv.partition(basis, molecule, 'eri')

        thresh_int = int(-math.log10(self.eri_thresh))

        fock_drv = FockDriver()
        pqrs = fock_drv.compute_eri(t4c_drv, naos, thresh_int)

        # AO-to-MO integral transformation

        tqrs = np.einsum('pqrs,pt->tqrs', pqrs, mo_coefs[0], optimize=True)
        del pqrs
        turs = np.einsum('tqrs,qu->turs', tqrs, mo_coefs[1], optimize=True)
        del tqrs
        tuvs = np.einsum('turs,rv->tuvs', turs, mo_coefs[2], optimize=True)
        del turs
        tuvw = np.einsum('tuvs,sw->tuvw', tuvs, mo_coefs[3], optimize=True)
        del tuvs

        if chem_notation:
            return tuvw
        elif phys_notation:
            return tuvw.swapaxes(1, 2)
