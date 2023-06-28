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
import numpy as np
import time as tm
import sys

from .veloxchemlib import ElectronRepulsionIntegralsDriver
from .veloxchemlib import mpi_master
from .veloxchemlib import denmat, fockmat, molorb
from .aodensitymatrix import AODensityMatrix
from .aofockmatrix import AOFockMatrix
from .molecularorbitals import MolecularOrbitals
from .outputstream import OutputStream
from .mointsdriver import MOIntegralsDriver
from .subcommunicators import SubCommunicators
from .qqscheme import get_qq_scheme
from .qqscheme import get_qq_type
from .sanitychecks import mp2_sanity_check
from .errorhandler import assert_msg_critical
from .inputparser import parse_input


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
        - qq_type: The electron repulsion integrals screening scheme.
        - eri_thresh: The electron repulsion integrals screening threshold.
        - batch_size: The number of Fock matrices in each batch.
        - comm_size: The size of each subcommunicator.
        - ostream: The output stream.
        - conventional: The flag for using conventional (in-memory) AO-to-MO
          integral transformation.
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
        self.qq_type = 'QQ_DEN'
        self.eri_thresh = 1.0e-12
        self.batch_size = 100

        # size of subcommunicator
        self.comm_size = 1

        # output stream
        self.ostream = ostream

        # use conventional (in-memory) AO-to-MO integral transformation?
        self.conventional = False

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
            'qq_type': 'str_upper',
            'eri_thresh': 'float',
            'batch_size': 'int',
            'comm_size': 'int',
            'conventional': 'bool',
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

    def compute(self, molecule, basis, scf_inp, scf_type=None):
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
        :param scf_type:
            The SCF type (restricted, unrestricted, restricted_openshell).
        """

        if self.rank == mpi_master():

            # use case: mp2_drv.compute(molecule, basis, scf_results)

            if isinstance(scf_inp, dict):
                scf_results = scf_inp

                if scf_type is not None:
                    assert_msg_critical(scf_type == scf_results['scf_type'],
                                        'Mp2Driver: Inconsistent SCF type')

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

            # use case: mp2_drv.compute(molecule, basis, mol_orbs, scf_type)

            elif isinstance(scf_inp, MolecularOrbitals):
                mol_orbs = scf_inp

            # invalid use case

            else:
                assert_msg_critical(False, 'Mp2Driver.compute: Invalid input')

        else:
            scf_type = None
            mol_orbs = MolecularOrbitals()

        scf_type = self.comm.bcast(scf_type, root=mpi_master())

        # run MP2

        if self.conventional:
            return self.compute_conventional(molecule, basis, mol_orbs,
                                             scf_type)
        else:
            return self.compute_distributed(molecule, basis, mol_orbs, scf_type)

    def compute_conventional(self, molecule, basis, mol_orbs, scf_type=None):
        """
        Performs conventional MP2 calculation.

        :param molecule:
            The molecule.
        :param basis:
            The AO basis set.
        :param mol_orbs:
            The molecular orbitals.
        :param scf_type:
            The SCF type (restricted, unrestricted, restricted_openshell).
        """

        # synchronize mol_orbs

        if self.rank != mpi_master():
            mol_orbs = MolecularOrbitals()
        mol_orbs.broadcast(self.rank, self.comm)

        # sanity check for MP2

        scf_type = mp2_sanity_check(molecule, basis, mol_orbs, scf_type)

        # compute MP2 in memory

        moints_drv = MOIntegralsDriver(self.comm, self.ostream)

        if self.rank == mpi_master():

            self.e_mp2 = 0.0

            if scf_type == 'restricted':

                orb_ene = mol_orbs.ea_to_numpy()
                nocc = molecule.number_of_alpha_electrons()
                eocc = orb_ene[:nocc]
                evir = orb_ene[nocc:]
                e_vv = evir.reshape(-1, 1) + evir

                phys_oovv = moints_drv.compute_in_memory(
                    molecule, basis, mol_orbs, 'phys_oovv')
                for i in range(phys_oovv.shape[0]):
                    for j in range(phys_oovv.shape[1]):
                        ab = phys_oovv[i, j, :, :]
                        denom = e_vv - eocc[i] - eocc[j]
                        self.e_mp2 -= np.sum(ab * (2.0 * ab - ab.T) / denom)

                mp2_str = '*** MP2 correlation energy: %20.12f a.u.' % self.e_mp2
                self.ostream.print_header(mp2_str.ljust(92))
                self.ostream.print_blank()

            elif scf_type == 'unrestricted':

                orb_ene_a = mol_orbs.ea_to_numpy()
                orb_ene_b = mol_orbs.eb_to_numpy()

                nocc_a = molecule.number_of_alpha_electrons()
                nocc_b = molecule.number_of_beta_electrons()

                eocc_a = orb_ene_a[:nocc_a]
                evir_a = orb_ene_a[nocc_a:]

                eocc_b = orb_ene_b[:nocc_b]
                evir_b = orb_ene_b[nocc_b:]

                e_vv_aa = evir_a.reshape(-1, 1) + evir_a
                e_vv_bb = evir_b.reshape(-1, 1) + evir_b
                e_vv_ab = evir_a.reshape(-1, 1) + evir_b

                phys_oovv_aaaa = moints_drv.compute_in_memory(
                    molecule, basis, mol_orbs, 'phys_oovv', 'aaaa')
                phys_oovv_bbbb = moints_drv.compute_in_memory(
                    molecule, basis, mol_orbs, 'phys_oovv', 'bbbb')
                phys_oovv_abab = moints_drv.compute_in_memory(
                    molecule, basis, mol_orbs, 'phys_oovv', 'abab')

                for i in range(phys_oovv_aaaa.shape[0]):
                    for j in range(phys_oovv_aaaa.shape[1]):
                        vv = phys_oovv_aaaa[i, j, :, :]
                        denom = e_vv_aa - eocc_a[i] - eocc_a[j]
                        self.e_mp2 -= 0.5 * np.sum(vv * (vv - vv.T) / denom)

                for i in range(phys_oovv_bbbb.shape[0]):
                    for j in range(phys_oovv_bbbb.shape[1]):
                        vv = phys_oovv_bbbb[i, j, :, :]
                        denom = e_vv_bb - eocc_b[i] - eocc_b[j]
                        self.e_mp2 -= 0.5 * np.sum(vv * (vv - vv.T) / denom)

                for i in range(phys_oovv_abab.shape[0]):
                    for j in range(phys_oovv_abab.shape[1]):
                        vv = phys_oovv_abab[i, j, :, :]
                        denom = e_vv_ab - eocc_a[i] - eocc_b[j]
                        self.e_mp2 -= np.sum(vv * vv / denom)

                mp2_str = '*** MP2 correlation energy: %20.12f a.u.' % self.e_mp2
                self.ostream.print_header(mp2_str.ljust(92))
                self.ostream.print_blank()

            elif scf_type == 'restricted_openshell':

                assert_msg_critical(
                    False, 'Restricted open-shell MP2 not implemented')

            return {'mp2_energy': self.e_mp2}
        else:
            return None

    def compute_distributed(self, molecule, basis, mol_orbs, scf_type=None):
        """
        Performs MP2 calculation via distributed Fock builds.

        :param molecule:
            The molecule.
        :param basis:
            The AO basis set.
        :param mol_orbs:
            The molecular orbitals.
        :param scf_type:
            The SCF type (restricted, unrestricted, restricted_openshell).
        """

        # synchronize mol_orbs

        if self.rank != mpi_master():
            mol_orbs = MolecularOrbitals()
        mol_orbs.broadcast(self.rank, self.comm)

        # sanity check for MP2

        scf_type = mp2_sanity_check(molecule, basis, mol_orbs, scf_type)

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

        # screening data

        eri_drv = ElectronRepulsionIntegralsDriver(local_comm)
        screening = eri_drv.compute(get_qq_scheme(self.qq_type),
                                    self.eri_thresh, molecule, basis)

        # prepare MO integrals

        if local_master:
            e_mp2 = 0.0

            mol_orbs.broadcast(cross_comm.Get_rank(), cross_comm)

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

            if scf_type == 'restricted':

                mo_ints_ids = [((i, j), 'aa')
                               for i in range(nocc_a)
                               for j in range(i + 1, nocc_a)]
                mo_ints_ids += [((i, i), 'aa') for i in range(nocc_a)]

            elif scf_type == 'unrestricted':

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

            self.print_header(len(mo_ints_ids))
            valstr = 'Monitoring calculation on master node.'
            self.ostream.print_header(valstr.ljust(80))
            self.ostream.print_blank()

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

            num_batches = count[cross_rank] // self.batch_size
            if count[cross_rank] % self.batch_size != 0:
                num_batches += 1
        else:
            num_batches = None

        num_batches = local_comm.bcast(num_batches, root=mpi_master())

        # compute MO integrals in batches

        batch_t0 = tm.time()

        for batch_ind in range(num_batches):

            if local_master:

                batch_start = mo_ints_start + batch_ind * self.batch_size
                batch_end = min(batch_start + self.batch_size, mo_ints_end)
                batch_ids = mo_ints_ids[batch_start:batch_end]

                dks = []

                for (i, j), spin in batch_ids:

                    if spin == 'aa':
                        mo_ij_aa = np.zeros((nocc_a, nocc_a))
                        mo_ij_aa[i, j] = 1.0
                        dks.append(
                            np.linalg.multi_dot(
                                [mo_occ_a, mo_ij_aa, mo_occ_a.T]))

                    elif spin == 'bb':
                        mo_ij_bb = np.zeros((nocc_b, nocc_b))
                        mo_ij_bb[i, j] = 1.0
                        dks.append(
                            np.linalg.multi_dot(
                                [mo_occ_b, mo_ij_bb, mo_occ_b.T]))

                    elif spin == 'ab':
                        mo_ij_ab = np.zeros((nocc_a, nocc_b))
                        mo_ij_ab[i, j] = 1.0
                        dks.append(
                            np.linalg.multi_dot(
                                [mo_occ_a, mo_ij_ab, mo_occ_b.T]))

                # Note: use restricted AODensityMatrix and Fock build
                dens = AODensityMatrix(dks, denmat.rest)

            else:
                dens = AODensityMatrix()

            dens.broadcast(local_comm.Get_rank(), local_comm)

            fock = AOFockMatrix(dens)
            for i in range(fock.number_of_fock_matrices()):
                fock.set_fock_type(fockmat.rgenk, i)

            eri_drv.compute(fock, dens, molecule, basis, screening)
            fock.reduce_sum(local_comm.Get_rank(), local_comm.Get_size(),
                            local_comm)

            if local_master:

                if scf_type == 'restricted':

                    for ind, ((i, j), spin) in enumerate(batch_ids):
                        f_aa = fock.alpha_to_numpy(ind)
                        vv = np.linalg.multi_dot([mo_vir_a.T, f_aa, mo_vir_a])
                        de = orb_ene_a[i] + orb_ene_a[j] - evv_aa
                        e_mp2 += np.sum(vv * (2.0 * vv - vv.T) / de)
                        if i != j:
                            e_mp2 += np.sum(vv.T * (2.0 * vv.T - vv) / de)

                elif scf_type == 'unrestricted':

                    for ind, ((i, j), spin) in enumerate(batch_ids):
                        if spin == 'aa':
                            f_aa = fock.alpha_to_numpy(ind)
                            vv = np.linalg.multi_dot(
                                [mo_vir_a.T, f_aa, mo_vir_a])
                            de = orb_ene_a[i] + orb_ene_a[j] - evv_aa
                            e_mp2 += 0.5 * np.sum(vv * (vv - vv.T) / de)
                            if i != j:
                                e_mp2 += 0.5 * np.sum(vv.T * (vv.T - vv) / de)

                        elif spin == 'bb':
                            # Note: use alpha_to_numpy due to restricted Fock
                            f_bb = fock.alpha_to_numpy(ind)
                            vv = np.linalg.multi_dot(
                                [mo_vir_b.T, f_bb, mo_vir_b])
                            de = orb_ene_b[i] + orb_ene_b[j] - evv_bb
                            e_mp2 += 0.5 * np.sum(vv * (vv - vv.T) / de)
                            if i != j:
                                e_mp2 += 0.5 * np.sum(vv.T * (vv.T - vv) / de)

                        elif spin == 'ab':
                            # Note: use alpha_to_numpy due to restricted Fock
                            f_ab = fock.alpha_to_numpy(ind)
                            vv = np.linalg.multi_dot(
                                [mo_vir_a.T, f_ab, mo_vir_b])
                            de = orb_ene_a[i] + orb_ene_b[j] - evv_ab
                            e_mp2 += np.sum(vv * vv / de)

            if global_master:
                valstr = '{:d} / {:d}'.format(batch_end - mo_ints_start,
                                              mo_ints_end - mo_ints_start)
                valstr += ' Fock matrices processed. Time: {:.2f} sec'.format(
                    tm.time() - batch_t0)
                self.ostream.print_header(valstr.ljust(80))
                self.ostream.print_blank()
                self.ostream.flush()

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

    def print_header(self, num_matrices):
        """
        Prints header for the MP2 driver.

        :param num_matrices:
            The number of Fock matrices to be computed.
        """

        self.ostream.print_blank()
        self.ostream.print_header('MP2 Driver Setup')
        self.ostream.print_header(18 * '=')
        self.ostream.print_blank()

        str_width = 60
        cur_str = 'Number of Fock Matrices      : ' + str(num_matrices)
        self.ostream.print_header(cur_str.ljust(str_width))
        cur_str = 'Size of Fock Matrices Batch  : ' + str(self.batch_size)
        self.ostream.print_header(cur_str.ljust(str_width))

        cur_str = 'Number of Subcommunicators   : '
        cur_str += str(self.nodes // self.comm_size)
        self.ostream.print_header(cur_str.ljust(str_width))

        cur_str = 'ERI Screening Scheme         : ' + get_qq_type(self.qq_type)
        self.ostream.print_header(cur_str.ljust(str_width))
        cur_str = 'ERI Screening Threshold      : {:.1e}'.format(
            self.eri_thresh)
        self.ostream.print_header(cur_str.ljust(str_width))
        self.ostream.print_blank()
        self.ostream.flush()
