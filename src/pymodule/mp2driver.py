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
from .veloxchemlib import mpi_master, mat_t
from .veloxchemlib import make_matrix
from .matrix import Matrix
from .fockdriver import FockDriver
from .molecularorbitals import MolecularOrbitals, molorb
from .outputstream import OutputStream
from .subcommunicators import SubCommunicators
from .sanitychecks import molecule_sanity_check
from .errorhandler import assert_msg_critical
from .inputparser import parse_input
from .batchsize import get_batch_size


class Mp2Driver:
    """
    Implements MP2 driver.

    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.

    Instance variable
        - e_mp2: The MP2 correlation energy.
        - comm: The MPI communicator.
        - rank: The MPI rank.
        - nodes: Number of MPI processes.
        - eri_thresh: The electron repulsion integrals screening threshold.
        - comm_size: The size of each subcommunicator.
        - ostream: The output stream.
    """

    def __init__(self, comm=None, ostream=None):
        """
        Initializes MP2 driver.
        """

        if comm is None:
            comm = MPI.COMM_WORLD

        if ostream is None:
            if comm.Get_rank() == mpi_master():
                ostream = OutputStream(sys.stdout)
            else:
                ostream = OutputStream(None)

        # mp2 energy
        self.e_mp2 = None

        # mpi information
        self.comm = comm
        self.rank = self.comm.Get_rank()
        self.nodes = self.comm.Get_size()

        # screening scheme and batch size for Fock build
        self.eri_thresh = 1.0e-12

        # size of subcommunicator
        self.comm_size = 1

        # output stream
        self.ostream = ostream

    def update_settings(self, mp2_dict, method_dict=None):
        """
        Updates settings in MP2 driver.

        :param mp2_dict:
            The dictionary of MP2 input.
        :param method_dict:
            The dictionary of method settings.
        """

        if method_dict is None:
            method_dict = {}

        mp2_keywords = {
            'eri_thresh': 'float',
            'comm_size': 'int',
        }

        parse_input(self, mp2_keywords, mp2_dict)

        if self.comm_size != 1:
            if self.nodes % self.comm_size != 0:
                self.comm_size = 1

        if 'xcfun' in method_dict:
            errmsg = 'Mp2Driver: The \'xcfun\' keyword is not supported in MP2 '
            errmsg += 'calculation.'
            if self.rank == mpi_master():
                assert_msg_critical(False, errmsg)

        if 'potfile' in method_dict:
            errmsg = 'Mp2Driver: The \'potfile\' keyword is not supported in '
            errmsg += 'MP2 calculation.'
            if self.rank == mpi_master():
                assert_msg_critical(False, errmsg)

        if 'electric_field' in method_dict:
            errmsg = 'Mp2Driver: The \'electric field\' keyword is not '
            errmsg += 'supported in MP2 calculation.'
            if self.rank == mpi_master():
                assert_msg_critical(False, errmsg)

    def compute(self, molecule, basis, scf_inp):
        """
        Performs MP2 calculation.

        :param molecule:
            The molecule.
        :param basis:
            The AO basis set.
        :param scf_inp:
            The input from SCF. In practice scf_inp could be scf_results, which
            is a dictionary returned by SCF driver, or mol_orbs, which is a
            molecular orbitals object from SCF driver.
        """

        if self.rank == mpi_master():

            # use case: mp2_drv.compute(molecule, basis, scf_results)

            if isinstance(scf_inp, dict):
                scf_results = scf_inp

                scf_type = scf_results['scf_type']

                mo_a = scf_results['C_alpha']
                ene_a = scf_results['E_alpha']
                occ_a = molecule.get_aufbau_alpha_occupation(ene_a.shape[0])

                if scf_type == 'restricted':
                    mol_orbs = MolecularOrbitals([mo_a], [ene_a], [occ_a],
                                                 molorb.rest)

                elif scf_type == 'unrestricted':
                    mo_b = scf_results['C_beta']
                    ene_b = scf_results['E_beta']
                    occ_b = molecule.get_aufbau_beta_occupation(ene_b.shape[0])
                    mol_orbs = MolecularOrbitals([mo_a, mo_b], [ene_a, ene_b],
                                                 [occ_a, occ_b], molorb.unrest)

                elif scf_type == 'restricted_openshell':
                    occ_b = molecule.get_aufbau_beta_occupation(ene_a.shape[0])
                    mol_orbs = MolecularOrbitals([mo_a], [ene_a],
                                                 [occ_a, occ_b],
                                                 molorb.restopen)

                else:
                    assert_msg_critical(False,
                                        'Mp2Driver.compute: Invalid SCF type')

            # use case: mp2_drv.compute(molecule, basis, mol_orbs)

            elif isinstance(scf_inp, MolecularOrbitals):
                mol_orbs = scf_inp

            # invalid use case

            else:
                assert_msg_critical(False, 'Mp2Driver.compute: Invalid input')

        else:
            mol_orbs = MolecularOrbitals()

        # run MP2

        return self.compute_distributed(molecule, basis, mol_orbs)

    def compute_distributed(self, molecule, basis, mol_orbs):
        """
        Performs MP2 calculation via distributed Fock builds.

        :param molecule:
            The molecule.
        :param basis:
            The AO basis set.
        :param mol_orbs:
            The molecular orbitals.
        """

        # synchronize mol_orbs

        if self.rank != mpi_master():
            mol_orbs = MolecularOrbitals()
        mol_orbs = mol_orbs.broadcast(self.comm, mpi_master())

        # sanity check

        molecule_sanity_check(molecule)

        assert_msg_critical(
            basis.get_dimensions_of_basis() == mol_orbs.number_of_aos(),
            'Mp2Driver.compute: Inconsistent number of AOs in basis set ' +
            'and molecular orbitals')

        assert_msg_critical(
            mol_orbs.get_orbitals_type() != molorb.restopen,
            'Mp2Driver.compute: Restricted open-shell MP2 not implemented')

        # screening

        thresh_int = int(-math.log10(self.eri_thresh))

        if self.rank == mpi_master():
            screening = T4CScreener()
            screening.partition(basis, molecule, 'eri')
        else:
            screening = None
        screening = self.comm.bcast(screening, root=mpi_master())

        # subcommunicators

        if self.rank == mpi_master():
            assert_msg_critical(self.nodes % self.comm_size == 0,
                                'MP2 driver: invalid size of subcommunicator')

        grps = [p // self.comm_size for p in range(self.nodes)]
        subcomm = SubCommunicators(self.comm, grps)
        local_comm = subcomm.local_comm
        cross_comm = subcomm.cross_comm

        cross_rank = cross_comm.Get_rank()
        cross_nodes = cross_comm.Get_size()

        local_master = (local_comm.Get_rank() == mpi_master())
        global_master = (self.rank == mpi_master())

        # prepare MO integrals

        # TODO: double check batch size
        batch_size = get_batch_size(None, 100, mol_orbs.number_of_aos(),
                                    self.comm)

        if local_master:
            e_mp2 = 0.0

            mol_orbs = mol_orbs.broadcast(cross_comm, mpi_master())

            nocc_a = molecule.number_of_alpha_electrons()
            nocc_b = molecule.number_of_beta_electrons()

            mo_a = mol_orbs.alpha_to_numpy()
            mo_b = mol_orbs.beta_to_numpy()

            mo_occ_a = mo_a[:, :nocc_a].copy()
            mo_vir_a = mo_a[:, nocc_a:].copy()

            mo_occ_b = mo_b[:, :nocc_b].copy()
            mo_vir_b = mo_b[:, nocc_b:].copy()

            orb_ene_a = mol_orbs.ea_to_numpy()
            orb_ene_b = mol_orbs.eb_to_numpy()

            evir_a = orb_ene_a[nocc_a:]
            evir_b = orb_ene_b[nocc_b:]

            evv_aa = evir_a.reshape(-1, 1) + evir_a
            evv_bb = evir_b.reshape(-1, 1) + evir_b
            evv_ab = evir_a.reshape(-1, 1) + evir_b

            # restricted

            if mol_orbs.get_orbitals_type() == molorb.rest:

                mo_ints_ids = [((i, j), 'aa')
                               for i in range(nocc_a)
                               for j in range(i + 1, nocc_a)]
                mo_ints_ids += [((i, i), 'aa') for i in range(nocc_a)]

            # unrestricted

            elif mol_orbs.get_orbitals_type() == molorb.unrest:

                mo_ints_ids = [((i, j), 'aa')
                               for i in range(nocc_a)
                               for j in range(i + 1, nocc_a)]
                mo_ints_ids += [((i, i), 'aa') for i in range(nocc_a)]

                mo_ints_ids += [((i, j), 'bb')
                                for i in range(nocc_b)
                                for j in range(i + 1, nocc_b)]
                mo_ints_ids += [((i, i), 'bb') for i in range(nocc_b)]

                mo_ints_ids += [
                    ((i, j), 'ab') for i in range(nocc_a) for j in range(nocc_b)
                ]

            self.print_header(len(mo_ints_ids), batch_size)
            valstr = 'Monitoring calculation on master node.'
            self.ostream.print_header(valstr.ljust(80))
            self.ostream.print_blank()

        else:
            mo_ints_ids = None

        mo_ints_ids = local_comm.bcast(mo_ints_ids, root=mpi_master())

        ave, res = divmod(len(mo_ints_ids), cross_nodes)
        count = [ave + 1 if i < res else ave for i in range(cross_nodes)]
        displ = [sum(count[:i]) for i in range(cross_nodes)]

        valstr = '{:d} Fock matrices '.format(count[cross_rank])
        valstr += 'will be processed on master node.'
        self.ostream.print_header(valstr.ljust(80))
        self.ostream.print_blank()
        self.ostream.flush()

        mo_ints_start = displ[cross_rank]
        mo_ints_end = mo_ints_start + count[cross_rank]

        if local_master:
            num_batches = count[cross_rank] // batch_size
            if count[cross_rank] % batch_size != 0:
                num_batches += 1
        else:
            num_batches = None
        num_batches = local_comm.bcast(num_batches, root=mpi_master())

        # compute MO integrals in batches

        for batch_ind in range(num_batches):

            batch_start = mo_ints_start + batch_ind * batch_size
            batch_end = min(batch_start + batch_size, mo_ints_end)
            batch_ids = mo_ints_ids[batch_start:batch_end]

            dks = []

            for (i, j), spin in batch_ids:

                if local_master:

                    if spin == 'aa':
                        mo_ij_aa = np.zeros((nocc_a, nocc_a))
                        mo_ij_aa[i, j] = 1.0
                        ao_dens = np.linalg.multi_dot(
                            [mo_occ_a, mo_ij_aa, mo_occ_a.T])

                    elif spin == 'bb':
                        mo_ij_bb = np.zeros((nocc_b, nocc_b))
                        mo_ij_bb[i, j] = 1.0
                        ao_dens = np.linalg.multi_dot(
                            [mo_occ_b, mo_ij_bb, mo_occ_b.T])

                    elif spin == 'ab':
                        mo_ij_ab = np.zeros((nocc_a, nocc_b))
                        mo_ij_ab[i, j] = 1.0
                        ao_dens = np.linalg.multi_dot(
                            [mo_occ_a, mo_ij_ab, mo_occ_b.T])

                else:
                    ao_dens = None

                ao_dens = local_comm.bcast(ao_dens, root=mpi_master())

                dks.append(ao_dens)

            fock_drv = FockDriver(local_comm)

            fock_type = 'k'
            exchange_scaling_factor = 1.0

            fock_arrays = []

            for idx in range(len(batch_ids)):
                den_mat_for_fock = make_matrix(basis, mat_t.general)
                den_mat_for_fock.set_values(dks[idx])

                fock_mat = fock_drv.compute(screening, den_mat_for_fock,
                                            fock_type, exchange_scaling_factor,
                                            0.0, thresh_int)

                fock_np = fock_mat.to_numpy()
                fock_mat = Matrix()

                fock_np = local_comm.reduce(fock_np, root=mpi_master())

                fock_arrays.append(fock_np)

            if local_master:

                # restricted

                if mol_orbs.get_orbitals_type() == molorb.rest:

                    for ind, ((i, j), spin) in enumerate(batch_ids):
                        f_aa = fock_arrays[ind]
                        vv = np.linalg.multi_dot([mo_vir_a.T, f_aa, mo_vir_a])
                        de = orb_ene_a[i] + orb_ene_a[j] - evv_aa
                        e_mp2 += np.sum(vv * (2.0 * vv - vv.T) / de)
                        if i != j:
                            e_mp2 += np.sum(vv.T * (2.0 * vv.T - vv) / de)

                # unrestricted

                elif mol_orbs.get_orbitals_type() == molorb.unrest:

                    for ind, ((i, j), spin) in enumerate(batch_ids):
                        if spin == 'aa':
                            f_aa = fock_arrays[ind]
                            vv = np.linalg.multi_dot(
                                [mo_vir_a.T, f_aa, mo_vir_a])
                            de = orb_ene_a[i] + orb_ene_a[j] - evv_aa
                            e_mp2 += 0.5 * np.sum(vv * (vv - vv.T) / de)
                            if i != j:
                                e_mp2 += 0.5 * np.sum(vv.T * (vv.T - vv) / de)

                        elif spin == 'bb':
                            f_bb = fock_arrays[ind]
                            vv = np.linalg.multi_dot(
                                [mo_vir_b.T, f_bb, mo_vir_b])
                            de = orb_ene_b[i] + orb_ene_b[j] - evv_bb
                            e_mp2 += 0.5 * np.sum(vv * (vv - vv.T) / de)
                            if i != j:
                                e_mp2 += 0.5 * np.sum(vv.T * (vv.T - vv) / de)

                        elif spin == 'ab':
                            f_ab = fock_arrays[ind]
                            vv = np.linalg.multi_dot(
                                [mo_vir_a.T, f_ab, mo_vir_b])
                            de = orb_ene_a[i] + orb_ene_b[j] - evv_ab
                            e_mp2 += np.sum(vv * vv / de)

        if local_master:
            e_mp2 = cross_comm.reduce(e_mp2, root=mpi_master())

        if global_master:
            self.e_mp2 = e_mp2
            mp2_str = '*** MP2 correlation energy: {:20.12f} a.u. '.format(
                self.e_mp2)
            self.ostream.print_header(mp2_str.ljust(80))
            self.ostream.print_blank()
            self.ostream.flush()

            return {'mp2_energy': self.e_mp2}
        else:
            return None

    def print_header(self, num_matrices, batch_size):
        """
        Prints header for the MP2 driver.

        :param num_matrices:
            The number of Fock matrices to be computed.
        :param batch_size:
            The batch size.
        """

        self.ostream.print_blank()
        self.ostream.print_header('MP2 Driver Setup')
        self.ostream.print_header(18 * '=')
        self.ostream.print_blank()

        str_width = 60
        cur_str = 'Number of Fock Matrices      : ' + str(num_matrices)
        self.ostream.print_header(cur_str.ljust(str_width))

        cur_str = 'Number of Subcommunicators   : '
        cur_str += str(self.nodes // self.comm_size)
        self.ostream.print_header(cur_str.ljust(str_width))

        cur_str = 'ERI Screening Threshold      : {:.1e}'.format(
            self.eri_thresh)
        self.ostream.print_header(cur_str.ljust(str_width))
        self.ostream.print_blank()
        self.ostream.flush()
