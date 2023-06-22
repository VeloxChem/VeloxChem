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
from pathlib import Path
import numpy as np
import time
import sys

from .veloxchemlib import mpi_master
from .outputstream import OutputStream
from .distributedarray import DistributedArray
from .cppsolver import ComplexResponse
from .linearsolver import LinearSolver
from .tpadriver import TpaDriver
from .checkpoint import check_distributed_focks
from .checkpoint import read_distributed_focks
from .checkpoint import write_distributed_focks


class TpaReducedDriver(TpaDriver):
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

    def get_densities(self, wi, Nx, mo, nocc, norb):
        """
        Computes the compounded densities needed for the compounded Fock
        matrices F^{σ} used for the reduced iostropic cubic response function

        :param wi:
            A list of the frequencies
        :param kX:
            A dictonary with all the first-order response vectors in
            distributed form
        :param mo:
            A matrix containing the MO coefficents
        :param nocc:
            Number of occupied orbitals
        :param norb:
            Number of orbitals

        :return:
            A list of tranformed compounded densities
        """

        distributed_density_1 = None
        distributed_density_2 = None

        for w in wi:

            nx = ComplexResponse.get_full_solution_vector(Nx[('x', w)])
            ny = ComplexResponse.get_full_solution_vector(Nx[('y', w)])
            nz = ComplexResponse.get_full_solution_vector(Nx[('z', w)])

            if self.rank == mpi_master():

                kx = self.complex_lrvec2mat(nx, nocc, norb)
                ky = self.complex_lrvec2mat(ny, nocc, norb)
                kz = self.complex_lrvec2mat(nz, nocc, norb)

                # create the first order single indexed densiteies #

                Dx = self.commut_mo_density(kx, nocc)
                Dy = self.commut_mo_density(ky, nocc)
                Dz = self.commut_mo_density(kz, nocc)

                # create the first order two indexed densities #

                # σ terms #

                Dxx = self.commut(kx, Dx)
                Dyy = self.commut(ky, Dy)
                Dzz = self.commut(kz, Dz)

                D_sig_xx = 6 * (3 * Dxx + Dyy + Dzz)
                D_sig_yy = 6 * (Dxx + 3 * Dyy + Dzz)
                D_sig_zz = 6 * (Dxx + Dyy + 3 * Dzz)

                D_sig_xy = 6 * (self.commut(ky, Dx) + self.commut(kx, Dy))
                D_sig_xz = 6 * (self.commut(kx, Dz) + self.commut(kz, Dx))
                D_sig_yz = 6 * (self.commut(ky, Dz) + self.commut(kz, Dy))

                # density transformation from MO to AO basis

                Dx = np.linalg.multi_dot([mo, Dx, mo.T])
                Dy = np.linalg.multi_dot([mo, Dy, mo.T])
                Dz = np.linalg.multi_dot([mo, Dz, mo.T])

                D_sig_xx = np.linalg.multi_dot([mo, D_sig_xx, mo.T])
                D_sig_yy = np.linalg.multi_dot([mo, D_sig_yy, mo.T])
                D_sig_zz = np.linalg.multi_dot([mo, D_sig_zz, mo.T])

                D_sig_xy = np.linalg.multi_dot([mo, D_sig_xy, mo.T])
                D_sig_xz = np.linalg.multi_dot([mo, D_sig_xz, mo.T])
                D_sig_yz = np.linalg.multi_dot([mo, D_sig_yz, mo.T])

                dist_den_1_freq = np.hstack((
                    Dx.real.reshape(-1, 1),
                    Dx.imag.reshape(-1, 1),
                    Dy.real.reshape(-1, 1),
                    Dy.imag.reshape(-1, 1),
                    Dz.real.reshape(-1, 1),
                    Dz.imag.reshape(-1, 1),
                ))

                dist_den_2_freq = np.hstack((
                    D_sig_xx.real.reshape(-1, 1),
                    D_sig_yy.real.reshape(-1, 1),
                    D_sig_zz.real.reshape(-1, 1),
                    D_sig_xy.real.reshape(-1, 1),
                    D_sig_xz.real.reshape(-1, 1),
                    D_sig_yz.real.reshape(-1, 1),
                ))

            else:

                dist_den_1_freq = None
                dist_den_2_freq = None

            dist_den_1_freq = DistributedArray(dist_den_1_freq, self.comm)
            dist_den_2_freq = DistributedArray(dist_den_2_freq, self.comm)

            if distributed_density_1 is None:
                distributed_density_1 = DistributedArray(dist_den_1_freq.data,
                                                         self.comm,
                                                         distribute=False)
            else:
                distributed_density_1.append(dist_den_1_freq, axis=1)

            if distributed_density_2 is None:
                distributed_density_2 = DistributedArray(dist_den_2_freq.data,
                                                         self.comm,
                                                         distribute=False)
            else:
                distributed_density_2.append(dist_den_2_freq, axis=1)

        return distributed_density_1, distributed_density_2, None

    def get_fock_dict(self, wi, density_list1, density_list2, density_list3,
                      F0_a, mo, molecule, ao_basis, dft_dict, profiler):
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
            self._print_fock_header()

        # generate key-frequency pairs

        key_freq_pairs = []

        for w in wi:
            for key in [
                    'f_sig_xx', 'f_sig_yy', 'f_sig_zz', 'f_sig_xy', 'f_sig_xz',
                    'f_sig_yz'
            ]:
                key_freq_pairs.append((key, w))

        # examine checkpoint file for distributed Focks

        if self.checkpoint_file is not None:
            fock_file = str(
                Path(self.checkpoint_file).with_suffix('.tpa_fock_1_red.h5'))
        else:
            fock_file = None

        if self.restart:
            if self.rank == mpi_master():
                self.restart = check_distributed_focks(fock_file,
                                                       key_freq_pairs)
            self.restart = self.comm.bcast(self.restart, mpi_master())

        # read or compute distributed Focks

        if self.restart:
            dist_focks = read_distributed_focks(fock_file, self.comm,
                                                self.ostream)
        else:
            time_start_fock = time.time()

            if self._dft:
                dist_focks = self._comp_nlr_fock(mo, molecule, ao_basis, 'real',
                                                 dft_dict, density_list1,
                                                 density_list2, None,
                                                 'redtpa_i', profiler)
            else:
                dist_focks = self._comp_nlr_fock(mo, molecule, ao_basis, 'real',
                                                 None, None, density_list2,
                                                 None, 'redtpa_i', profiler)

            self._print_fock_time(time.time() - time_start_fock)

            write_distributed_focks(fock_file, dist_focks, key_freq_pairs,
                                    self.comm, self.ostream)

        # assign distributed Focks to key-frequency pairs

        focks = {'F0': F0_a}

        for fock_index, (key, w) in enumerate(key_freq_pairs):
            if key not in focks:
                focks[key] = {}
            focks[key][w] = DistributedArray(dist_focks.data[:, fock_index],
                                             self.comm,
                                             distribute=False)

        return focks

    def get_Nxy(self, w, d_a_mo, X, fock_dict, Nx, nocc, norb, molecule,
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
        :param Nx:
            A dictonary containg all the response vectors in distributed form
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
        xy_dict = self.get_xy(d_a_mo, X, w, fock_dict, Nx, nocc, norb)

        # Frequencies to compute
        if self.rank == mpi_master():
            freq = tuple([sum(x) for x in zip(w, w)])
        else:
            freq = None
        freq = self.comm.bcast(freq, root=mpi_master())

        N_total_drv = ComplexResponse(self.comm, self.ostream)
        N_total_drv.frequencies = freq

        cpp_keywords = {
            'damping', 'norm_thresh', 'lindep_thresh', 'conv_thresh',
            'max_iter', 'eri_thresh', 'qq_type', 'timing', 'memory_profiling',
            'batch_size', 'restart', 'xcfun', 'grid_level', 'potfile',
            'electric_field', 'program_end_time'
        }

        for key in cpp_keywords:
            setattr(N_total_drv, key, getattr(self, key))

        if self.checkpoint_file is not None:
            N_total_drv.checkpoint_file = str(
                Path(self.checkpoint_file).with_suffix('.tpa_2_red.h5'))

        N_total_results = N_total_drv.compute(molecule, ao_basis, scf_tensors,
                                              xy_dict)

        self._is_converged = (self._is_converged and N_total_drv.is_converged)

        Nxy_dict = N_total_results['solutions']
        FXY_2_dict = N_total_results['focks']

        return (Nxy_dict, FXY_2_dict)

    def get_xy(self, d_a_mo, X, wi, Fock, Nx, nocc, norb):
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
        :param Nx:
            A dictonary with all the first-order response vector in distributed
            form
        :param nocc:
            The number of occupied orbitals
        :param norb:
            The number of total orbitals

        :return:
            A dictonary of compounded gradient vectors
        """

        xy_dict = {}

        one_third = 1.0 / 3.0

        for w in wi:

            vec_pack = np.array([
                Fock['Fb'][('x', w)].data,
                Fock['Fb'][('y', w)].data,
                Fock['Fb'][('z', w)].data,
                Fock['f_sig_xx'][w].data * one_third,
                Fock['f_sig_yy'][w].data * one_third,
                Fock['f_sig_zz'][w].data * one_third,
                Fock['f_sig_xy'][w].data * one_third,
                Fock['f_sig_xz'][w].data * one_third,
                Fock['f_sig_yz'][w].data * one_third,
            ]).T.copy()

            vec_pack = self._collect_vectors_in_columns(vec_pack)

            nx = ComplexResponse.get_full_solution_vector(Nx[('x', w)])
            ny = ComplexResponse.get_full_solution_vector(Nx[('y', w)])
            nz = ComplexResponse.get_full_solution_vector(Nx[('z', w)])

            if self.rank != mpi_master():
                continue

            vec_pack = vec_pack.T.copy().reshape(-1, norb, norb)

            (f_x, f_y, f_z, f_sig_xx, f_sig_yy, f_sig_zz, f_sig_xy, f_sig_xz,
             f_sig_yz) = vec_pack

            mu_x = X['x']
            mu_y = X['y']
            mu_z = X['z']

            kx = (self.complex_lrvec2mat(nx, nocc, norb)).T
            ky = (self.complex_lrvec2mat(ny, nocc, norb)).T
            kz = (self.complex_lrvec2mat(nz, nocc, norb)).T

            # REAL PART #
            F0 = Fock['F0']

            # BD σ gradients #

            xi_xx = self._xi(kx, kx, f_x, f_x, F0)
            xi_yy = self._xi(ky, ky, f_y, f_y, F0)
            xi_zz = self._xi(kz, kz, f_z, f_z, F0)

            x2_xx = self._x2_contract(kx.T, mu_x, d_a_mo, nocc, norb)
            x2_yy = self._x2_contract(ky.T, mu_y, d_a_mo, nocc, norb)
            x2_zz = self._x2_contract(kz.T, mu_z, d_a_mo, nocc, norb)

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
            mat = (self._xi(ky, kx, f_y, f_x, F0) +
                   self._xi(kx, ky, f_x, f_y, F0) + 0.5 * f_sig_xy).T
            xy_dict[key] = self.anti_sym(
                -LinearSolver.lrmat2vec(mat, nocc, norb))
            xy_dict[key] -= self._x2_contract(ky.T, mu_x, d_a_mo, nocc, norb)
            xy_dict[key] -= self._x2_contract(kx.T, mu_y, d_a_mo, nocc, norb)

            key = (('N_sig_xz', w), 2 * w)
            mat = (self._xi(kz, kx, f_z, f_x, F0) +
                   self._xi(kx, kz, f_x, f_z, F0) + 0.5 * f_sig_xz).T
            xy_dict[key] = self.anti_sym(
                -LinearSolver.lrmat2vec(mat, nocc, norb))
            xy_dict[key] -= self._x2_contract(kz.T, mu_x, d_a_mo, nocc, norb)
            xy_dict[key] -= self._x2_contract(kx.T, mu_z, d_a_mo, nocc, norb)

            key = (('N_sig_yz', w), 2 * w)
            mat = (self._xi(kz, ky, f_z, f_y, F0) +
                   self._xi(ky, kz, f_y, f_z, F0) + 0.5 * f_sig_yz).T
            xy_dict[key] = self.anti_sym(
                -LinearSolver.lrmat2vec(mat, nocc, norb))
            xy_dict[key] -= self._x2_contract(kz.T, mu_y, d_a_mo, nocc, norb)
            xy_dict[key] -= self._x2_contract(ky.T, mu_z, d_a_mo, nocc, norb)

        return xy_dict

    def get_densities_II(self, wi, Nx, Nxy, mo, nocc, norb):
        """
        Computes the compounded densities needed for the compounded
        second-order Fock matrices used for the reduced isotropic cubic response
        function. Note: All densities are 1/3 of those in the paper, and all
        the Fock matrices are later scaled by 3.

        :param wi:
            A list of the frequencies
        :param Nx:
            A dictonary with all the first-order response vectors in
            distributed form
        :param Nxy:
            A dict of the two index response vectors in distributed form
        :param mo:
            A matrix containing the MO coefficents
        :param nocc:
            Number of occupied orbitals
        :param norb:
            Number of orbitals

        :return:
            A list of tranformed compounded densities
        """

        distributed_density_1 = None
        distributed_density_2 = None

        for w in wi:

            nx = ComplexResponse.get_full_solution_vector(Nx[('x', w)])
            ny = ComplexResponse.get_full_solution_vector(Nx[('y', w)])
            nz = ComplexResponse.get_full_solution_vector(Nx[('z', w)])

            n_sig_xx = ComplexResponse.get_full_solution_vector(
                Nxy[(('N_sig_xx', w), 2 * w)])
            n_sig_yy = ComplexResponse.get_full_solution_vector(
                Nxy[(('N_sig_yy', w), 2 * w)])
            n_sig_zz = ComplexResponse.get_full_solution_vector(
                Nxy[(('N_sig_zz', w), 2 * w)])
            n_sig_xy = ComplexResponse.get_full_solution_vector(
                Nxy[(('N_sig_xy', w), 2 * w)])
            n_sig_xz = ComplexResponse.get_full_solution_vector(
                Nxy[(('N_sig_xz', w), 2 * w)])
            n_sig_yz = ComplexResponse.get_full_solution_vector(
                Nxy[(('N_sig_yz', w), 2 * w)])

            if self.rank == mpi_master():

                k_sig_xx = self.complex_lrvec2mat(n_sig_xx, nocc, norb)
                k_sig_yy = self.complex_lrvec2mat(n_sig_yy, nocc, norb)
                k_sig_zz = self.complex_lrvec2mat(n_sig_zz, nocc, norb)
                k_sig_xy = self.complex_lrvec2mat(n_sig_xy, nocc, norb)
                k_sig_xz = self.complex_lrvec2mat(n_sig_xz, nocc, norb)
                k_sig_yz = self.complex_lrvec2mat(n_sig_yz, nocc, norb)

                kx = self.complex_lrvec2mat(nx, nocc, norb)
                ky = self.complex_lrvec2mat(ny, nocc, norb)
                kz = self.complex_lrvec2mat(nz, nocc, norb)

                kx_ = -kx.conj().T  # kX['Nc'][('x', -w)].T
                ky_ = -ky.conj().T  # kX['Nc'][('y', -w)].T
                kz_ = -kz.conj().T  # kX['Nc'][('z', -w)].T

                # SIGMA contributiatons #
                Dc_x_ = self.commut_mo_density(kx_, nocc)
                Dc_y_ = self.commut_mo_density(ky_, nocc)
                Dc_z_ = self.commut_mo_density(kz_, nocc)

                D_sig_xx = self.commut_mo_density(k_sig_xx, nocc)
                D_sig_yy = self.commut_mo_density(k_sig_yy, nocc)
                D_sig_zz = self.commut_mo_density(k_sig_zz, nocc)

                D_sig_xy = self.commut_mo_density(k_sig_xy, nocc)
                D_sig_xz = self.commut_mo_density(k_sig_xz, nocc)
                D_sig_yz = self.commut_mo_density(k_sig_yz, nocc)

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

                Dc_x_ = np.linalg.multi_dot([mo, Dc_x_, mo.T])
                Dc_y_ = np.linalg.multi_dot([mo, Dc_y_, mo.T])
                Dc_z_ = np.linalg.multi_dot([mo, Dc_z_, mo.T])

                Dx = np.linalg.multi_dot([mo, Dx, mo.T])
                Dy = np.linalg.multi_dot([mo, Dy, mo.T])
                Dz = np.linalg.multi_dot([mo, Dz, mo.T])

                D_sig_xx = np.linalg.multi_dot([mo, D_sig_xx, mo.T])
                D_sig_yy = np.linalg.multi_dot([mo, D_sig_yy, mo.T])
                D_sig_zz = np.linalg.multi_dot([mo, D_sig_zz, mo.T])
                D_sig_xy = np.linalg.multi_dot([mo, D_sig_xy, mo.T])
                D_sig_xz = np.linalg.multi_dot([mo, D_sig_xz, mo.T])
                D_sig_yz = np.linalg.multi_dot([mo, D_sig_yz, mo.T])

                dist_den_1_freq = np.hstack(
                    (Dc_x_.real.reshape(-1, 1), Dc_x_.imag.reshape(-1, 1),
                     Dc_y_.real.reshape(-1, 1), Dc_y_.imag.reshape(-1, 1),
                     Dc_z_.real.reshape(-1, 1), Dc_z_.imag.reshape(-1, 1),
                     D_sig_xx.real.reshape(-1, 1), D_sig_xx.imag.reshape(-1, 1),
                     D_sig_yy.real.reshape(-1, 1), D_sig_yy.imag.reshape(-1, 1),
                     D_sig_zz.real.reshape(-1, 1), D_sig_zz.imag.reshape(-1, 1),
                     D_sig_xy.real.reshape(-1, 1), D_sig_xy.imag.reshape(-1, 1),
                     D_sig_xz.real.reshape(-1, 1), D_sig_xz.imag.reshape(-1, 1),
                     D_sig_yz.real.reshape(-1, 1), D_sig_yz.imag.reshape(-1,
                                                                         1)))

                dist_den_2_freq = np.hstack(
                    (Dx.real.reshape(-1, 1), Dx.imag.reshape(-1, 1),
                     Dy.real.reshape(-1, 1), Dy.imag.reshape(-1, 1),
                     Dz.real.reshape(-1, 1), Dz.imag.reshape(-1, 1)))
            else:
                dist_den_1_freq = None
                dist_den_2_freq = None

            dist_den_1_freq = DistributedArray(dist_den_1_freq, self.comm)
            dist_den_2_freq = DistributedArray(dist_den_2_freq, self.comm)

            if distributed_density_1 is None:
                distributed_density_1 = DistributedArray(dist_den_1_freq.data,
                                                         self.comm,
                                                         distribute=False)
            else:
                distributed_density_1.append(dist_den_1_freq, axis=1)

            if distributed_density_2 is None:
                distributed_density_2 = DistributedArray(dist_den_2_freq.data,
                                                         self.comm,
                                                         distribute=False)
            else:
                distributed_density_2.append(dist_den_2_freq, axis=1)

        return distributed_density_1, distributed_density_2

    def get_fock_dict_II(self, wi, density_list1, density_list2, mo, molecule,
                         ao_basis, dft_dict, profiler):
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
            self._print_fock_header()

        # generate key-frequnecy pairs

        key_freq_pairs = []

        for w in wi:
            for key in ['F123_x', 'F123_y', 'F123_z']:
                key_freq_pairs.append((key, w))

        # examine checkpoint file for distributed Focks

        if self.checkpoint_file is not None:
            fock_file = str(
                Path(self.checkpoint_file).with_suffix('.tpa_fock_2_red.h5'))
        else:
            fock_file = None

        if self.restart:
            if self.rank == mpi_master():
                self.restart = check_distributed_focks(fock_file,
                                                       key_freq_pairs)
            self.restart = self.comm.bcast(self.restart, mpi_master())

        # read or compute distributed Focks

        if self.restart:
            dist_focks = read_distributed_focks(fock_file, self.comm,
                                                self.ostream)
        else:
            time_start_fock = time.time()

            if self._dft:
                dist_focks = self._comp_nlr_fock(mo, molecule, ao_basis,
                                                 'real_and_imag', dft_dict,
                                                 density_list1, density_list2,
                                                 None, 'redtpa_ii', profiler)
            else:
                dist_focks = self._comp_nlr_fock(mo, molecule, ao_basis,
                                                 'real_and_imag', None, None,
                                                 density_list2, None,
                                                 'redtpa_ii', profiler)

            self._print_fock_time(time.time() - time_start_fock)

            write_distributed_focks(fock_file, dist_focks, key_freq_pairs,
                                    self.comm, self.ostream)

        # assign distributed Focks to key-frequency pairs

        focks = {}

        for fock_index, (key, w) in enumerate(key_freq_pairs):
            if key not in focks:
                focks[key] = {}
            focks[key][w] = DistributedArray(dist_focks.data[:, fock_index],
                                             self.comm,
                                             distribute=False)

        return focks

    def get_e3(self, wi, Nx, Nxy, fo, fo2, nocc, norb):
        """
        Contracts E[3]NxNyz for the isotropic cubic response function. Takes
        the Fock matrices from fock_dict and fock_dict_II and contracts them
        with the first and second-order response vectors.

        :param wi:
             A list of freqs
        :param Nx:
            A dict of the single index response vectors in distributed form
        :param Nxy:
            A dict of the two index response vectors in distributed form
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

            vec_pack = self._collect_vectors_in_columns(vec_pack)

            nx = ComplexResponse.get_full_solution_vector(Nx[('x', w)])
            ny = ComplexResponse.get_full_solution_vector(Nx[('y', w)])
            nz = ComplexResponse.get_full_solution_vector(Nx[('z', w)])

            n_sig_xx = ComplexResponse.get_full_solution_vector(
                Nxy[(('N_sig_xx', w), 2 * w)])
            n_sig_yy = ComplexResponse.get_full_solution_vector(
                Nxy[(('N_sig_yy', w), 2 * w)])
            n_sig_zz = ComplexResponse.get_full_solution_vector(
                Nxy[(('N_sig_zz', w), 2 * w)])
            n_sig_xy = ComplexResponse.get_full_solution_vector(
                Nxy[(('N_sig_xy', w), 2 * w)])
            n_sig_xz = ComplexResponse.get_full_solution_vector(
                Nxy[(('N_sig_xz', w), 2 * w)])
            n_sig_yz = ComplexResponse.get_full_solution_vector(
                Nxy[(('N_sig_yz', w), 2 * w)])

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

            k_x = (self.complex_lrvec2mat(nx, nocc, norb)).T
            k_y = (self.complex_lrvec2mat(ny, nocc, norb)).T
            k_z = (self.complex_lrvec2mat(nz, nocc, norb)).T

            k_x_ = -k_x.conj().T  # kX['Nc'][('x', -w)].T
            k_y_ = -k_y.conj().T  # kX['Nc'][('y', -w)].T
            k_z_ = -k_z.conj().T  # kX['Nc'][('z', -w)].T

            k_sig_xx = (self.complex_lrvec2mat(n_sig_xx, nocc, norb)).T
            k_sig_yy = (self.complex_lrvec2mat(n_sig_yy, nocc, norb)).T
            k_sig_zz = (self.complex_lrvec2mat(n_sig_zz, nocc, norb)).T
            k_sig_xy = (self.complex_lrvec2mat(n_sig_xy, nocc, norb)).T
            k_sig_xz = (self.complex_lrvec2mat(n_sig_xz, nocc, norb)).T
            k_sig_yz = (self.complex_lrvec2mat(n_sig_yz, nocc, norb)).T

            # x #

            zeta_sig_xx = self._xi(k_x_, k_sig_xx, f_x_, f_sig_xx, F0_a)
            zeta_sig_yy = self._xi(k_x_, k_sig_yy, f_x_, f_sig_yy, F0_a)
            zeta_sig_zz = self._xi(k_x_, k_sig_zz, f_x_, f_sig_zz, F0_a)

            zeta_sig_xy = self._xi(k_y_, k_sig_xy, f_y_, f_sig_xy, F0_a)
            zeta_sig_xz = self._xi(k_z_, k_sig_xz, f_z_, f_sig_xz, F0_a)

            X_terms = (zeta_sig_xx + zeta_sig_xy + zeta_sig_xz).T + (0.5 *
                                                                     F123_x).T
            ff_x = -2 * LinearSolver.lrmat2vec(X_terms, nocc, norb)
            ff_x = self.anti_sym(ff_x)
            f_iso_x[w] = ff_x

            # y #

            zeta_sig_yx = self._xi(k_x_, k_sig_xy, f_x_, f_sig_xy, F0_a)
            zeta_sig_yy = self._xi(k_y_, k_sig_yy, f_y_, f_sig_yy, F0_a)
            zeta_sig_yz = self._xi(k_z_, k_sig_yz, f_z_, f_sig_yz, F0_a)

            Y_terms = (zeta_sig_yx + zeta_sig_yy + zeta_sig_yz).T + (0.5 *
                                                                     F123_y).T
            ff_y = -2 * LinearSolver.lrmat2vec(Y_terms, nocc, norb)
            ff_y = self.anti_sym(ff_y)
            f_iso_y[w] = ff_y

            # z #

            zeta_sig_zx = self._xi(k_x_, k_sig_xz, f_x_, f_sig_xz, F0_a)
            zeta_sig_zy = self._xi(k_y_, k_sig_yz, f_y_, f_sig_yz, F0_a)
            zeta_sig_zz = self._xi(k_z_, k_sig_zz, f_z_, f_sig_zz, F0_a)

            Z_terms = (zeta_sig_zx + zeta_sig_zy + zeta_sig_zz).T + (0.5 *
                                                                     F123_z).T
            ff_z = -2 * LinearSolver.lrmat2vec(Z_terms, nocc, norb)
            ff_z = self.anti_sym(ff_z)
            f_iso_z[w] = ff_z

        return {'f_iso_x': f_iso_x, 'f_iso_y': f_iso_y, 'f_iso_z': f_iso_z}

    def get_other_terms(self, wi, track, X, Nx, Nxy, da, nocc, norb):
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
        :param Nx:
            A dictonary with all the respone vectors in distributed form
        :param Nxy:
            A dictonary containing all the two-index response vectors in
            distributed form
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
                Na = Nx[(op_a, wa)]
                A = X[op_a]

                for op_c in 'xyz':
                    op_ac = op_a + op_c if op_a <= op_c else op_c + op_a

                    # BD
                    Nbd = Nxy[(('N_sig_' + op_ac, w), wbd)]
                    Nc = Nx[(op_c, -wc)]
                    C = X[op_c]

                    inp_list.append({
                        'flag': 'BD',
                        'freq': w,
                        'Nbd': Nbd,
                        'Na': Na,
                        'Nc': Nc,
                        'A': A,
                        'C': C,
                    })

        list_x2_a2 = [self.get_x2_a2(inp, da, nocc, norb) for inp in inp_list]

        if self.rank == mpi_master():
            for terms in list_x2_a2:
                if terms['key'] not in na_x2_nyz_dict:
                    na_x2_nyz_dict[terms['key']] = 0.0
                if terms['key'] not in nx_a2_nyz_dict:
                    nx_a2_nyz_dict[terms['key']] = 0.0
                na_x2_nyz_dict[terms['key']] += terms['x2']
                nx_a2_nyz_dict[terms['key']] += terms['a2']

            return {
                'NaX2Nyz': na_x2_nyz_dict,
                'NxA2Nyz': nx_a2_nyz_dict,
            }

        return None

    def print_results(self, freqs, comp, result):
        """
        Prints the results from the reduced TPA calculation.

        :param freqs:
            List of frequencies
        :param comp:
            List of gamma tensors components
        :param result:
            A dictonary containing the isotropic gamma, T[4], T[3], X[3], A[3],
            X[2] and A[2] contractions.
        """

        gamma = result['gamma']

        t3_dict = result['t3_dict']

        NaX2Nyz = result['NaX2Nyz']
        NxA2Nyz = result['NxA2Nyz']

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

            self._print_component('T3', w, t3_dict[w, -w, w], width)
            self._print_component('X2', w, NaX2Nyz[w, -w, w], width)
            self._print_component('A2', w, NxA2Nyz[w, -w, w], width)
            self._print_component('gamma', w, gamma[w, -w, w], width)

            self.ostream.print_blank()

        title = 'Reference: '
        title += 'K. Ahmadzadeh, M. Scott, M. Brand, O. Vahtras, X. Li, '
        self.ostream.print_header(title.ljust(width))
        title = 'Z. Rinkevicius, and P. Norman, '
        title += 'J. Chem. Phys. 154, 024111 (2021)'
        self.ostream.print_header(title.ljust(width))
        self.ostream.print_blank()

        spectrum = self.get_spectrum(result, x_unit='au')
        self.print_spectrum(spectrum, width)

        self.ostream.print_blank()
