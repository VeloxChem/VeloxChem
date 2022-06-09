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

from mpi4py import MPI
from pathlib import Path
import numpy as np
import time
import sys

from .veloxchemlib import mpi_master
from .outputstream import OutputStream
from .distributedarray import DistributedArray
from .cppsolver import ComplexResponse
from .linearsolver import LinearSolver
from .tpadriver import TPADriver
from .checkpoint import check_distributed_focks
from .checkpoint import read_distributed_focks
from .checkpoint import write_distributed_focks


class TPAReducedDriver(TPADriver):
    """
    Implements the reduced isotropic cubic response driver for two-photon
    absorption (TPA)

    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.
    """

    def __init__(self, comm=None, ostream=None):
        """
        Initializes the reduced isotropic cubic response driver for two-photon
        absorption (TPA)
        """

        if comm is None:
            comm = MPI.COMM_WORLD

        if ostream is None:
            if comm.Get_rank() == mpi_master():
                ostream = OutputStream(sys.stdout)
            else:
                ostream = OutputStream(None)

        super().__init__(comm, ostream)

    def update_settings(self, rsp_dict, method_dict=None):
        """
        Updates response and method settings for TPA

        :param rsp_dict:
            The dictionary of response input.
        :param method_dict:
            The dictionary of method settings.
        """

        if method_dict is None:
            method_dict = {}

        super().update_settings(rsp_dict, method_dict)

    def get_densities(self, wi, kX, mo, nocc):
        """
        Computes the compounded densities needed for the compounded Fock
        matrices F^{σ} used for the reduced iostropic cubic response function

        :param wi:
            A list of the frequencies
        :param kX:
            A dictonary with all the first-order response matrices
        :param mo:
            A matrix containing the MO coefficents
        :param nocc:
            Number of alpha electrons

        :return:
            A list of tranformed compounded densities
        """

        density_list = []

        D_mo = np.zeros((mo.shape[1], mo.shape[1]))
        D_mo[:nocc, :nocc] = np.eye(nocc)

        for w in wi:

            kx = kX['Nb'][('x', w)]
            ky = kX['Nb'][('y', w)]
            kz = kX['Nb'][('z', w)]

            # create the first order single indexed densiteies #

            Dx = self.commut(kx, D_mo)
            Dy = self.commut(ky, D_mo)
            Dz = self.commut(kz, D_mo)

            # create the first order two indexed densities #

            # σ terms #

            Dxx = self.commut(kx, Dx)
            Dyy = self.commut(ky, Dy)
            Dzz = self.commut(kz, Dz)

            D_sig_xx = 2 * (3 * Dxx + Dyy + Dzz)
            D_sig_yy = 2 * (Dxx + 3 * Dyy + Dzz)
            D_sig_zz = 2 * (Dxx + Dyy + 3 * Dzz)

            D_sig_xy = 2 * (self.commut(ky, Dx) + self.commut(kx, Dy))
            D_sig_xz = 2 * (self.commut(kx, Dz) + self.commut(kz, Dx))
            D_sig_yz = 2 * (self.commut(ky, Dz) + self.commut(kz, Dy))

            # density transformation from MO to AO basis

            D_sig_xx = np.linalg.multi_dot([mo, D_sig_xx, mo.T])
            D_sig_yy = np.linalg.multi_dot([mo, D_sig_yy, mo.T])
            D_sig_zz = np.linalg.multi_dot([mo, D_sig_zz, mo.T])

            D_sig_xy = np.linalg.multi_dot([mo, D_sig_xy, mo.T])
            D_sig_xz = np.linalg.multi_dot([mo, D_sig_xz, mo.T])
            D_sig_yz = np.linalg.multi_dot([mo, D_sig_yz, mo.T])

            density_list.append(D_sig_xx.real)
            density_list.append(D_sig_yy.real)
            density_list.append(D_sig_zz.real)
            density_list.append(D_sig_xy.real)
            density_list.append(D_sig_xz.real)
            density_list.append(D_sig_yz.real)

        return density_list

    def get_fock_dict(self, wi, density_list, F0_a, mo, molecule, ao_basis):
        """
        Computes the compounded Fock matrices F^{σ}  used for the reduced
        isotropic cubic response function

        :param wi:
            A list of the frequencies
        :param density_list:
            A list of tranformed compounded densities
        :param F0_a:
            The Fock matrix in MO basis
        :param mo:
            A matrix containing the MO coefficents
        :param molecule:
            The molecule
        :param ao_basis:
            The AO basis set

        :return:
            A dictonary of compounded first-order Fock-matrices
        """

        if self.rank == mpi_master():
            self.print_fock_header()

        keys = [
            'f_sig_xx',
            'f_sig_yy',
            'f_sig_zz',
            'f_sig_xy',
            'f_sig_xz',
            'f_sig_yz',
        ]

        if self.checkpoint_file is not None:
            fock_file = str(
                Path(self.checkpoint_file).with_suffix('.tpa_fock_1_red.h5'))
        else:
            fock_file = None

        if self.restart:
            if self.rank == mpi_master():
                self.restart = check_distributed_focks(fock_file, keys, wi)
            self.restart = self.comm.bcast(self.restart, mpi_master())

        if self.restart:
            focks = read_distributed_focks(fock_file, keys, wi, self.comm,
                                           self.ostream)
            focks['F0'] = F0_a
            return focks

        time_start_fock = time.time()
        dist_focks = self.comp_nlr_fock(mo, molecule, ao_basis, 'real', None,
                                        None, density_list, 'tpa')
        time_end_fock = time.time()

        total_time_fock = time_end_fock - time_start_fock
        self.print_fock_time(total_time_fock)

        focks = {'F0': F0_a}
        for key in keys:
            focks[key] = {}

        fock_index = 0
        for w in wi:
            for key in keys:
                focks[key][w] = DistributedArray(dist_focks.data[:, fock_index],
                                                 self.comm,
                                                 distribute=False)
                fock_index += 1

        write_distributed_focks(fock_file, focks, keys, wi, self.comm,
                                self.ostream)

        return focks

    def get_Nxy(self, w, d_a_mo, X, fock_dict, kX, nocc, norb, molecule,
                ao_basis, scf_tensors):
        """
        Computes all the second-order response vectors needed for the reduced
        isotropic cubic response computation

        :param w:
            A list of all the frequencies
        :param d_a_mo:
            The density matrix in MO basis
        :param X:
            Dipole integrals
        :param fock_dict:
            A dictonary containing all the Fock matricies
        :param kX:
            A dictonary containg all the response matricies
        :param nocc:
            The number of occupied orbitals
        :param norb:
            The number of total orbitals
        :param molecule:
            The molecule.
        :param basis:
            The AO basis.
        :param scf_tensors:
            The dictionary of tensors from converged SCF wavefunction.

        :return:
            A dictonary of Fock matrices from the subspace,second-order
            response vectors and second-order response matrices
        """

        # Get the second-order gradiants
        xy_dict = self.get_xy(d_a_mo, X, w, fock_dict, kX, nocc, norb)

        # Frequencies to compute
        if self.rank == mpi_master():
            freq = tuple([sum(x) for x in zip(w, w)])
        else:
            freq = None
        freq = self.comm.bcast(freq, root=mpi_master())

        N_total_drv = ComplexResponse(self.comm, self.ostream)
        N_total_drv.frequencies = freq

        cpp_keywords = {
            'damping', 'lindep_thresh', 'conv_thresh', 'max_iter', 'eri_thresh',
            'qq_type', 'timing', 'memory_profiling', 'batch_size', 'restart',
            'program_end_time'
        }

        for key in cpp_keywords:
            setattr(N_total_drv, key, getattr(self, key))

        if self.checkpoint_file is not None:
            N_total_drv.checkpoint_file = str(
                Path(self.checkpoint_file).with_suffix('.tpa_2_red.h5'))

        N_total_results = N_total_drv.compute(molecule, ao_basis, scf_tensors,
                                              xy_dict)

        kXY_dict = N_total_results['kappas']
        FXY_2_dict = N_total_results['focks']

        return (kXY_dict, FXY_2_dict)

    def get_xy(self, d_a_mo, X, wi, Fock, kX, nocc, norb):
        """
        Computes the compounded gradient vectors N^{σ}  used for the reduced
        isotropic cubic response function

        :param d_a_mo:
            The SCF density matrix in MO basis
        :param X:
            Dipole integrals
        :param wi:
            A list of the frequencies
        :param Fock:
            A dictonary containing all the Fock matricies
        :param kX:
            A dictonary with all the first-order response matrices
        :param nocc:
            The number of occupied orbitals
        :param norb:
            The number of total orbitals

        :return:
            A dictonary of compounded gradient vectors
        """

        xy_dict = {}

        for w in wi:

            vec_pack = np.array([
                Fock['Fb'][('x', w)].data,
                Fock['Fb'][('y', w)].data,
                Fock['Fb'][('z', w)].data,
                Fock['f_sig_xx'][w].data,
                Fock['f_sig_yy'][w].data,
                Fock['f_sig_zz'][w].data,
                Fock['f_sig_xy'][w].data,
                Fock['f_sig_xz'][w].data,
                Fock['f_sig_yz'][w].data,
            ]).T.copy()

            vec_pack = self.collect_vectors_in_columns(vec_pack)

            if self.rank != mpi_master():
                continue

            vec_pack = vec_pack.T.copy().reshape(-1, norb, norb)

            (f_x, f_y, f_z, f_sig_xx, f_sig_yy, f_sig_zz, f_sig_xy, f_sig_xz,
             f_sig_yz) = vec_pack

            mu_x = X['x']
            mu_y = X['y']
            mu_z = X['z']

            kx = kX['Nb'][('x', w)].T
            ky = kX['Nb'][('y', w)].T
            kz = kX['Nb'][('z', w)].T

            # REAL PART #
            F0 = Fock['F0']

            # BD σ gradients #

            xi_xx = self.xi(kx, kx, f_x, f_x, F0)
            xi_yy = self.xi(ky, ky, f_y, f_y, F0)
            xi_zz = self.xi(kz, kz, f_z, f_z, F0)

            x2_xx = self.x2_contract(kx.T, mu_x, d_a_mo, nocc, norb)
            x2_yy = self.x2_contract(ky.T, mu_y, d_a_mo, nocc, norb)
            x2_zz = self.x2_contract(kz.T, mu_z, d_a_mo, nocc, norb)

            key = (('N_sig_xx', w), 2 * w)
            mat = (3 * xi_xx + xi_yy + xi_zz + 0.5 * f_sig_xx).T
            xy_dict[key] = self.anti_sym(
                -LinearSolver.lrmat2vec(mat, nocc, norb))
            xy_dict[key] -= (3 * x2_xx + x2_yy + x2_zz)

            key = (('N_sig_yy', w), 2 * w)
            mat = (xi_xx + 3 * xi_yy + xi_zz + 0.5 * f_sig_yy).T
            xy_dict[key] = self.anti_sym(
                -LinearSolver.lrmat2vec(mat, nocc, norb))
            xy_dict[key] -= (x2_xx + 3 * x2_yy + x2_zz)

            key = (('N_sig_zz', w), 2 * w)
            mat = (xi_xx + xi_yy + 3 * xi_zz + 0.5 * f_sig_zz).T
            xy_dict[key] = self.anti_sym(
                -LinearSolver.lrmat2vec(mat, nocc, norb))
            xy_dict[key] -= (x2_xx + x2_yy + 3 * x2_zz)

            key = (('N_sig_xy', w), 2 * w)
            mat = (self.xi(ky, kx, f_y, f_x, F0) +
                   self.xi(kx, ky, f_x, f_y, F0) + 0.5 * f_sig_xy).T
            xy_dict[key] = self.anti_sym(
                -LinearSolver.lrmat2vec(mat, nocc, norb))
            xy_dict[key] -= self.x2_contract(ky.T, mu_x, d_a_mo, nocc, norb)
            xy_dict[key] -= self.x2_contract(kx.T, mu_y, d_a_mo, nocc, norb)

            key = (('N_sig_xz', w), 2 * w)
            mat = (self.xi(kz, kx, f_z, f_x, F0) +
                   self.xi(kx, kz, f_x, f_z, F0) + 0.5 * f_sig_xz).T
            xy_dict[key] = self.anti_sym(
                -LinearSolver.lrmat2vec(mat, nocc, norb))
            xy_dict[key] -= self.x2_contract(kz.T, mu_x, d_a_mo, nocc, norb)
            xy_dict[key] -= self.x2_contract(kx.T, mu_z, d_a_mo, nocc, norb)

            key = (('N_sig_yz', w), 2 * w)
            mat = (self.xi(kz, ky, f_z, f_y, F0) +
                   self.xi(ky, kz, f_y, f_z, F0) + 0.5 * f_sig_yz).T
            xy_dict[key] = self.anti_sym(
                -LinearSolver.lrmat2vec(mat, nocc, norb))
            xy_dict[key] -= self.x2_contract(kz.T, mu_y, d_a_mo, nocc, norb)
            xy_dict[key] -= self.x2_contract(ky.T, mu_z, d_a_mo, nocc, norb)

        return xy_dict

    def get_densities_II(self, wi, kX, kXY, mo, nocc):
        """
        Computes the compounded densities needed for the compounded
        second-order Fock matrices used for the reduced isotropic cubic response
        function. Note: All densities are 1/3 of those in the paper, and all
        the Fock matrices are later scaled by 3.

        :param wi:
            A list of the frequencies
        :param kX:
            A dictonary with all the first-order response matrices
        :param kXY:
            A dict of the two index response matrices
        :param mo:
            A matrix containing the MO coefficents
        :param nocc:
            Number of alpha electrons

        :return:
            A list of tranformed compounded densities
        """

        density_list = []

        D_mo = np.zeros((mo.shape[1], mo.shape[1]))
        D_mo[:nocc, :nocc] = np.eye(nocc)

        for w in wi:
            k_sig_xx = kXY[(('N_sig_xx', w), 2 * w)]
            k_sig_yy = kXY[(('N_sig_yy', w), 2 * w)]
            k_sig_zz = kXY[(('N_sig_zz', w), 2 * w)]

            k_sig_xy = kXY[(('N_sig_xy', w), 2 * w)]
            k_sig_xz = kXY[(('N_sig_xz', w), 2 * w)]
            k_sig_yz = kXY[(('N_sig_yz', w), 2 * w)]

            kx = kX['Nb'][('x', w)]
            ky = kX['Nb'][('y', w)]
            kz = kX['Nb'][('z', w)]

            kx_ = -kx.conj().T  # kX['Nc'][('x', -w)].T
            ky_ = -ky.conj().T  # kX['Nc'][('y', -w)].T
            kz_ = -kz.conj().T  # kX['Nc'][('z', -w)].T

            # SIGMA contributiatons #
            Dc_x_ = self.commut(kx_, D_mo)
            Dc_y_ = self.commut(ky_, D_mo)
            Dc_z_ = self.commut(kz_, D_mo)

            D_sig_xx = self.commut(k_sig_xx, D_mo)
            D_sig_yy = self.commut(k_sig_yy, D_mo)
            D_sig_zz = self.commut(k_sig_zz, D_mo)

            D_sig_xy = self.commut(k_sig_xy, D_mo)
            D_sig_xz = self.commut(k_sig_xz, D_mo)
            D_sig_yz = self.commut(k_sig_yz, D_mo)

            # x #
            Dx = self.commut(kx_, D_sig_xx)
            Dx += self.commut(k_sig_xx, Dc_x_)
            Dx += self.commut(ky_, D_sig_xy)
            Dx += self.commut(k_sig_xy, Dc_y_)

            Dx += self.commut(kz_, D_sig_xz)
            Dx += self.commut(k_sig_xz, Dc_z_)

            # y #
            Dy = self.commut(kx_, D_sig_xy)
            Dy += self.commut(k_sig_xy, Dc_x_)

            Dy += self.commut(ky_, D_sig_yy)
            Dy += self.commut(k_sig_yy, Dc_y_)

            Dy += self.commut(kz_, D_sig_yz)
            Dy += self.commut(k_sig_yz, Dc_z_)

            # z #
            Dz = self.commut(kx_, D_sig_xz)
            Dz += self.commut(k_sig_xz, Dc_x_)

            Dz += self.commut(ky_, D_sig_yz)
            Dz += self.commut(k_sig_yz, Dc_y_)

            Dz += self.commut(kz_, D_sig_zz)
            Dz += self.commut(k_sig_zz, Dc_z_)

            # density transformation from MO to AO basis

            Dx = np.linalg.multi_dot([mo, Dx, mo.T])
            Dy = np.linalg.multi_dot([mo, Dy, mo.T])
            Dz = np.linalg.multi_dot([mo, Dz, mo.T])

            density_list.append(Dx.real)
            density_list.append(Dx.imag)
            density_list.append(Dy.real)
            density_list.append(Dy.imag)
            density_list.append(Dz.real)
            density_list.append(Dz.imag)

        return density_list

    def get_fock_dict_II(self, wi, density_list, mo, molecule, ao_basis):
        """
        Computes the compounded second-order Fock matrices used for the
        isotropic cubic response function

        :param wi:
            A list of the frequencies
        :param density_list:
            A list of tranformed compounded densities
        :param mo:
            A matrix containing the MO coefficents
        :param molecule:
            The molecule
        :param ao_basis:
            The AO basis set

        :return:
            A dictonary of compounded second-order Fock-matrices
        """

        if self.rank == mpi_master():
            self.print_fock_header()

        keys = ['F123_x', 'F123_y', 'F123_z']

        if self.checkpoint_file is not None:
            fock_file = str(
                Path(self.checkpoint_file).with_suffix('.tpa_fock_2_red.h5'))
        else:
            fock_file = None

        if self.restart:
            if self.rank == mpi_master():
                self.restart = check_distributed_focks(fock_file, keys, wi)
            self.restart = self.comm.bcast(self.restart, mpi_master())

        if self.restart:
            return read_distributed_focks(fock_file, keys, wi, self.comm,
                                          self.ostream)

        time_start_fock = time.time()
        dist_focks = self.comp_nlr_fock(mo, molecule, ao_basis, 'real_and_imag',
                                        None, None, density_list, 'tpa')
        time_end_fock = time.time()

        total_time_fock = time_end_fock - time_start_fock
        self.print_fock_time(total_time_fock)

        focks = {}
        for key in keys:
            focks[key] = {}

        fock_index = 0
        for w in wi:
            for key in keys:
                focks[key][w] = DistributedArray(dist_focks.data[:, fock_index],
                                                 self.comm,
                                                 distribute=False)
                fock_index += 1

        write_distributed_focks(fock_file, focks, keys, wi, self.comm,
                                self.ostream)

        return focks

    def get_e3(self, wi, kX, kXY, fo, fo2, nocc, norb):
        """
        Contracts E[3]NxNyz for the isotropic cubic response function. Takes
        the Fock matrices from fock_dict and fock_dict_II and contracts them
        with the first and second-order response vectors.

        :param wi:
             A list of freqs
        :param kX:
            A dict of the single index response matricies
        :param kXY:
            A dict of the two index response matrices
        :param fo:
            A dictonary of transformed Fock matricies from fock_dict
        :param fo2:
            A dictonarty of transfromed Fock matricies from fock_dict_two
        :param nocc:
            The number of occupied orbitals
        :param norb:
            The total number of orbitals

        :return:
            A dictonary of compounded E[3] tensors for the isotropic reduced
            cubic response function for TPA
        """

        f_iso_x = {}
        f_iso_y = {}
        f_iso_z = {}

        for w in wi:

            vec_pack = np.array([
                fo['Fb'][('x', w)].data,
                fo['Fb'][('y', w)].data,
                fo['Fb'][('z', w)].data,
                fo2[(('N_sig_xx', w), 2 * w)].data,
                fo2[(('N_sig_yy', w), 2 * w)].data,
                fo2[(('N_sig_zz', w), 2 * w)].data,
                fo2[(('N_sig_xy', w), 2 * w)].data,
                fo2[(('N_sig_xz', w), 2 * w)].data,
                fo2[(('N_sig_yz', w), 2 * w)].data,
                fo2['F123_x'][w].data,
                fo2['F123_y'][w].data,
                fo2['F123_z'][w].data,
            ]).T.copy()

            vec_pack = self.collect_vectors_in_columns(vec_pack)

            if self.rank != mpi_master():
                continue

            vec_pack = vec_pack.T.copy().reshape(-1, norb, norb)

            (f_x, f_y, f_z, f_sig_xx, f_sig_yy, f_sig_zz, f_sig_xy, f_sig_xz,
             f_sig_yz, F123_x, F123_y, F123_z) = vec_pack

            f_x_ = f_x.T.conj()  # fo['Fc'][('x', -w)]
            f_y_ = f_y.T.conj()  # fo['Fc'][('y', -w)]
            f_z_ = f_z.T.conj()  # fo['Fc'][('z', -w)]

            f_sig_xx = f_sig_xx.T.conj()
            f_sig_yy = f_sig_yy.T.conj()
            f_sig_zz = f_sig_zz.T.conj()

            f_sig_xy = f_sig_xy.T.conj()
            f_sig_xz = f_sig_xz.T.conj()
            f_sig_yz = f_sig_yz.T.conj()

            F0_a = fo['F0']

            # Response #

            k_x = kX['Nb'][('x', w)].T
            k_y = kX['Nb'][('y', w)].T
            k_z = kX['Nb'][('z', w)].T

            k_x_ = -k_x.conj().T  # kX['Nc'][('x', -w)].T
            k_y_ = -k_y.conj().T  # kX['Nc'][('y', -w)].T
            k_z_ = -k_z.conj().T  # kX['Nc'][('z', -w)].T

            k_sig_xx = kXY[(('N_sig_xx', w), 2 * w)].T
            k_sig_yy = kXY[(('N_sig_yy', w), 2 * w)].T
            k_sig_zz = kXY[(('N_sig_zz', w), 2 * w)].T

            k_sig_xy = kXY[(('N_sig_xy', w), 2 * w)].T
            k_sig_xz = kXY[(('N_sig_xz', w), 2 * w)].T
            k_sig_yz = kXY[(('N_sig_yz', w), 2 * w)].T

            # x #

            zeta_sig_xx = self.xi(k_x_, k_sig_xx, f_x_, f_sig_xx, F0_a)
            zeta_sig_yy = self.xi(k_x_, k_sig_yy, f_x_, f_sig_yy, F0_a)
            zeta_sig_zz = self.xi(k_x_, k_sig_zz, f_x_, f_sig_zz, F0_a)

            zeta_sig_xy = self.xi(k_y_, k_sig_xy, f_y_, f_sig_xy, F0_a)
            zeta_sig_xz = self.xi(k_z_, k_sig_xz, f_z_, f_sig_xz, F0_a)

            X_terms = (zeta_sig_xx + zeta_sig_xy + zeta_sig_xz).T + (0.5 *
                                                                     F123_x).T
            ff_x = -2 * LinearSolver.lrmat2vec(X_terms, nocc, norb)
            ff_x = self.anti_sym(ff_x)
            f_iso_x[w] = ff_x

            # y #

            zeta_sig_yx = self.xi(k_x_, k_sig_xy, f_x_, f_sig_xy, F0_a)
            zeta_sig_yy = self.xi(k_y_, k_sig_yy, f_y_, f_sig_yy, F0_a)
            zeta_sig_yz = self.xi(k_z_, k_sig_yz, f_z_, f_sig_yz, F0_a)

            Y_terms = (zeta_sig_yx + zeta_sig_yy + zeta_sig_yz).T + (0.5 *
                                                                     F123_y).T
            ff_y = -2 * LinearSolver.lrmat2vec(Y_terms, nocc, norb)
            ff_y = self.anti_sym(ff_y)
            f_iso_y[w] = ff_y

            # z #

            zeta_sig_zx = self.xi(k_x_, k_sig_xz, f_x_, f_sig_xz, F0_a)
            zeta_sig_zy = self.xi(k_y_, k_sig_yz, f_y_, f_sig_yz, F0_a)
            zeta_sig_zz = self.xi(k_z_, k_sig_zz, f_z_, f_sig_zz, F0_a)

            Z_terms = (zeta_sig_zx + zeta_sig_zy + zeta_sig_zz).T + (0.5 *
                                                                     F123_z).T
            ff_z = -2 * LinearSolver.lrmat2vec(Z_terms, nocc, norb)
            ff_z = self.anti_sym(ff_z)
            f_iso_z[w] = ff_z

        return {'f_iso_x': f_iso_x, 'f_iso_y': f_iso_y, 'f_iso_z': f_iso_z}

    def get_other_terms(self, wi, track, X, kX, kXY, da, nocc, norb):
        """
        Computes the terms involving X[2],A[2] in the reduced isotropic cubic
        response function

        :param wi:
            A list containing all the frequencies
        :param track:
            A list that contains information about what γ components that are
            to be computed and which freqs
        :param X:
            A dictonray with all the property integral matricies
        :param kX:
            A dictonary with all the respone matricies
        :param kXY:
            A dictonary containing all the two-index response matricies
        :param da:
            The SCF density matrix in MO basis
        :param nocc:
            The number of occupied orbitals
        :param norb:
            The total number of orbitals

        :return:
            A dictonary of final X[2],A[2] contraction values
        """

        na_x2_nyz_dict = {}
        nx_a2_nyz_dict = {}

        comp_per_freq = len(track) // len(wi)

        inp_list = []

        for i in range(len(wi)):
            vals = track[i * comp_per_freq].split(',')
            w = float(vals[1])
            wa = float(vals[1])
            wb = float(vals[1])
            wc = float(vals[2])
            wd = float(vals[3])

            wbd = wb + wd

            for op_a in 'xyz':
                Na_ka = kX['Na'][(op_a, wa)]
                A = X[op_a]

                for op_c in 'xyz':
                    op_ac = op_a + op_c if op_a <= op_c else op_c + op_a

                    # BD
                    kbd = kXY[(('N_sig_' + op_ac, w), wbd)]
                    kc_kb = kX['Nb'][(op_c, -wc)]
                    C = X[op_c]

                    inp_list.append({
                        'flag': 'BD',
                        'freq': w,
                        'kbd': kbd,
                        'Na_ka': Na_ka,
                        'kc_kb': kc_kb,
                        'A': A,
                        'C': C,
                    })

        ave, res = divmod(len(inp_list), self.nodes)
        counts = [ave + 1 if p < res else ave for p in range(self.nodes)]
        starts = [sum(counts[:p]) for p in range(self.nodes)]
        ends = [sum(counts[:p + 1]) for p in range(self.nodes)]

        list_x2_a2 = [
            self.get_x2_a2(inp, da, nocc, norb)
            for inp in inp_list[starts[self.rank]:ends[self.rank]]
        ]

        list_x2_a2 = self.comm.gather(list_x2_a2, root=mpi_master())

        if self.rank == mpi_master():
            for terms in list_x2_a2:
                for term in terms:
                    key = term['key']
                    if key not in na_x2_nyz_dict:
                        na_x2_nyz_dict[key] = 0.0
                    if key not in nx_a2_nyz_dict:
                        nx_a2_nyz_dict[key] = 0.0
                    na_x2_nyz_dict[key] += term['x2']
                    nx_a2_nyz_dict[key] += term['a2']

            return {
                'NaX2Nyz': na_x2_nyz_dict,
                'NxA2Nyz': nx_a2_nyz_dict,
            }

        return None

    def print_results(self, freqs, gamma, comp, t4_dict, t3_dict, other_dict):
        """
        Prints the results from the reduced TPA calculation.

        :param freqs:
            List of frequencies
        :param gamma:
            A dictonary containing the reduced isotropic cubic response
            functions for TPA
        :param comp:
            List of gamma tensors components
        :param t4_dict:
            A dictonary containing the isotropic T[4] contractions (None for
            one-photon off-resonance TPA calculations)
        :param t3_dict:
            A dictonary containing the isotropic T[3] contractions for
            one-photon off-resonance TPA calculations
        :param other_dict:
            A dictonary containing the isotropic X[2] and A[2] contractions for
            one-photo off-resonance TPA calculations
        """

        NaX2Nyz = other_dict['NaX2Nyz']
        NxA2Nyz = other_dict['NxA2Nyz']

        self.ostream.print_blank()

        w_str = 'The Isotropic Average gamma Tensor and Its'
        self.ostream.print_header(w_str)
        w_str = 'Isotropic Contributions at Given Frequencies'
        self.ostream.print_header(w_str)
        self.ostream.print_header('=' * (len(w_str) + 2))
        self.ostream.print_blank()

        w_str = '*** Note: The reduced expression is an approximation to the  '
        self.ostream.print_header(w_str)
        w_str = '    second-order nonlinear hyperpolarizability (gamma) and is'
        self.ostream.print_header(w_str)
        w_str = '    intended for use in one-photon off-resonance regions.    '
        self.ostream.print_header(w_str)
        self.ostream.print_blank()

        for w in freqs:
            title = '{:<9s} {:>12s} {:>20s} {:>21s}'.format(
                'Component', 'Frequency', 'Real', 'Imaginary')
            width = len(title)
            self.ostream.print_header(title.ljust(width))
            self.ostream.print_header(('-' * len(title)).ljust(width))

            self.print_component('T3', w, t3_dict[w, -w, w], width)
            self.print_component('X2', w, NaX2Nyz[w, -w, w], width)
            self.print_component('A2', w, NxA2Nyz[w, -w, w], width)
            self.print_component('gamma', w, gamma[w, -w, w], width)

            self.ostream.print_blank()

        title = 'Reference: '
        title += 'K. Ahmadzadeh, M. Scott, M. Brand, O. Vahtras, X. Li, '
        self.ostream.print_header(title.ljust(width))
        title = 'Z. Rinkevicius, and P. Norman, '
        title += 'J. Chem. Phys. 154, 024111 (2021)'
        self.ostream.print_header(title.ljust(width))
        self.ostream.print_blank()

        self.ostream.print_blank()
