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
import numpy as np
import sys

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

        fock_drv = FockDriver()
        pqrs = fock_drv.compute_eri(molecule, basis, self.eri_thresh)

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

    def compute_in_mem(self, molecule, basis, mol_orbs, mints_type):

        # for backward compatibility only
        return self.compute_in_memory(molecule, basis, mol_orbs,
                                      'phys_' + mints_type)
