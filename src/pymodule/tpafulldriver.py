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


class TpaFullDriver(TpaDriver):
    """
    Implements the full isotropic cubic response driver for two-photon
    absorption (TPA)

    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.
    """

    def __init__(self, comm=None, ostream=None):
        """
        Initializes the full isotropic cubic response driver for two-photon
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
        matrices F^{σ},F^{λ+τ},F^{σλτ} used for the isotropic cubic response
        function. Note: All densities are 1/3 of those in the paper, and all
        the Fock matrices are later scaled by 3.

        :param wi:
            A list of the frequencies
        :param Nx:
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
        distributed_density_3 = None

        one_third = 1.0 / 3.0

        for w in wi:

            nx = ComplexResponse.get_full_solution_vector(Nx[('x', w)])
            ny = ComplexResponse.get_full_solution_vector(Nx[('y', w)])
            nz = ComplexResponse.get_full_solution_vector(Nx[('z', w)])

            if self.rank == mpi_master():

                kx = self.complex_lrvec2mat(nx, nocc, norb)
                ky = self.complex_lrvec2mat(ny, nocc, norb)
                kz = self.complex_lrvec2mat(nz, nocc, norb)

                kx_ = -kx.conj().T  # kX['Nc'][('x', -w)].T
                ky_ = -ky.conj().T  # kX['Nc'][('y', -w)].T
                kz_ = -kz.conj().T  # kX['Nc'][('z', -w)].T

                # create the first order single indexed densiteies #

                Dx = self.commut_mo_density(kx, nocc)
                Dy = self.commut_mo_density(ky, nocc)
                Dz = self.commut_mo_density(kz, nocc)

                Dx_ = self.commut_mo_density(kx_, nocc)
                Dy_ = self.commut_mo_density(ky_, nocc)
                Dz_ = self.commut_mo_density(kz_, nocc)

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

                # λ+τ terms #

                Dxx = self.commut(kx_, Dx) + self.commut(kx, Dx_)
                Dyy = self.commut(ky_, Dy) + self.commut(ky, Dy_)
                Dzz = self.commut(kz_, Dz) + self.commut(kz, Dz_)

                D_lamtau_xx = 6 * (3 * Dxx + Dyy + Dzz)
                D_lamtau_yy = 6 * (Dxx + 3 * Dyy + Dzz)
                D_lamtau_zz = 6 * (Dxx + Dyy + 3 * Dzz)

                D_lamtau_xy = 6 * (self.commut(ky_, Dx) + self.commut(kx_, Dy))
                D_lamtau_xy += 6 * (self.commut(ky, Dx_) + self.commut(kx, Dy_))

                D_lamtau_xz = 6 * (self.commut(kx_, Dz) + self.commut(kz_, Dx))
                D_lamtau_xz += 6 * (self.commut(kx, Dz_) + self.commut(kz, Dx_))

                D_lamtau_yz = 6 * (self.commut(ky_, Dz) + self.commut(kz_, Dy))
                D_lamtau_yz += 6 * (self.commut(ky, Dz_) + self.commut(kz, Dy_))

                # Create first order three indexed Densities #

                D_lam_sig_tau_x = (self.commut(kx_, one_third * D_sig_xx) +
                                   self.commut(ky_, one_third * D_sig_xy) +
                                   self.commut(kz_, one_third * D_sig_xz))

                D_lam_sig_tau_x += (self.commut(kx, one_third * D_lamtau_xx) +
                                    self.commut(ky, one_third * D_lamtau_xy) +
                                    self.commut(kz, one_third * D_lamtau_xz))

                D_lam_sig_tau_y = (self.commut(kx_, one_third * D_sig_xy) +
                                   self.commut(ky_, one_third * D_sig_yy) +
                                   self.commut(kz_, one_third * D_sig_yz))

                D_lam_sig_tau_y += (self.commut(kx, one_third * D_lamtau_xy) +
                                    self.commut(ky, one_third * D_lamtau_yy) +
                                    self.commut(kz, one_third * D_lamtau_yz))

                D_lam_sig_tau_z = (self.commut(kx_, one_third * D_sig_xz) +
                                   self.commut(ky_, one_third * D_sig_yz) +
                                   self.commut(kz_, one_third * D_sig_zz))

                D_lam_sig_tau_z += (self.commut(kx, one_third * D_lamtau_xz) +
                                    self.commut(ky, one_third * D_lamtau_yz) +
                                    self.commut(kz, one_third * D_lamtau_zz))

                # density transformation from MO to AO basis

                Dx = np.linalg.multi_dot([mo, Dx, mo.T])
                Dy = np.linalg.multi_dot([mo, Dy, mo.T])
                Dz = np.linalg.multi_dot([mo, Dz, mo.T])
                Dx_ = np.linalg.multi_dot([mo, Dx_, mo.T])
                Dy_ = np.linalg.multi_dot([mo, Dy_, mo.T])
                Dz_ = np.linalg.multi_dot([mo, Dz_, mo.T])

                D_sig_xx = np.linalg.multi_dot([mo, D_sig_xx, mo.T])
                D_sig_yy = np.linalg.multi_dot([mo, D_sig_yy, mo.T])
                D_sig_zz = np.linalg.multi_dot([mo, D_sig_zz, mo.T])
                D_sig_xy = np.linalg.multi_dot([mo, D_sig_xy, mo.T])
                D_sig_xz = np.linalg.multi_dot([mo, D_sig_xz, mo.T])
                D_sig_yz = np.linalg.multi_dot([mo, D_sig_yz, mo.T])

                D_lamtau_xx = np.linalg.multi_dot([mo, D_lamtau_xx, mo.T])
                D_lamtau_yy = np.linalg.multi_dot([mo, D_lamtau_yy, mo.T])
                D_lamtau_zz = np.linalg.multi_dot([mo, D_lamtau_zz, mo.T])
                D_lamtau_xy = np.linalg.multi_dot([mo, D_lamtau_xy, mo.T])
                D_lamtau_xz = np.linalg.multi_dot([mo, D_lamtau_xz, mo.T])
                D_lamtau_yz = np.linalg.multi_dot([mo, D_lamtau_yz, mo.T])

                D_lam_sig_tau_x = np.linalg.multi_dot(
                    [mo, D_lam_sig_tau_x, mo.T])
                D_lam_sig_tau_y = np.linalg.multi_dot(
                    [mo, D_lam_sig_tau_y, mo.T])
                D_lam_sig_tau_z = np.linalg.multi_dot(
                    [mo, D_lam_sig_tau_z, mo.T])

                # Reshape all the density matrices as vectors and combine them
                # all into one object

                dist_den_1_freq = np.hstack((
                    Dx.real.reshape(-1, 1),
                    Dx.imag.reshape(-1, 1),
                    Dy.real.reshape(-1, 1),
                    Dy.imag.reshape(-1, 1),
                    Dz.real.reshape(-1, 1),
                    Dz.imag.reshape(-1, 1),
                    Dx_.real.reshape(-1, 1),
                    Dx_.imag.reshape(-1, 1),
                    Dy_.real.reshape(-1, 1),
                    Dy_.imag.reshape(-1, 1),
                    Dz_.real.reshape(-1, 1),
                    Dz_.imag.reshape(-1, 1),
                ))

                dist_den_2_freq = np.hstack((
                    D_sig_xx.real.reshape(-1, 1),
                    D_sig_xx.imag.reshape(-1, 1),
                    D_sig_yy.real.reshape(-1, 1),
                    D_sig_yy.imag.reshape(-1, 1),
                    D_sig_zz.real.reshape(-1, 1),
                    D_sig_zz.imag.reshape(-1, 1),
                    D_sig_xy.real.reshape(-1, 1),
                    D_sig_xy.imag.reshape(-1, 1),
                    D_sig_xz.real.reshape(-1, 1),
                    D_sig_xz.imag.reshape(-1, 1),
                    D_sig_yz.real.reshape(-1, 1),
                    D_sig_yz.imag.reshape(-1, 1),
                    D_lamtau_xx.real.reshape(-1, 1),
                    D_lamtau_xx.imag.reshape(-1, 1),
                    D_lamtau_yy.real.reshape(-1, 1),
                    D_lamtau_yy.imag.reshape(-1, 1),
                    D_lamtau_zz.real.reshape(-1, 1),
                    D_lamtau_zz.imag.reshape(-1, 1),
                    D_lamtau_xy.real.reshape(-1, 1),
                    D_lamtau_xy.imag.reshape(-1, 1),
                    D_lamtau_xz.real.reshape(-1, 1),
                    D_lamtau_xz.imag.reshape(-1, 1),
                    D_lamtau_yz.real.reshape(-1, 1),
                    D_lamtau_yz.imag.reshape(-1, 1),
                ))

                dist_den_3_freq = np.hstack((
                    D_lam_sig_tau_x.real.reshape(-1, 1),
                    D_lam_sig_tau_x.imag.reshape(-1, 1),
                    D_lam_sig_tau_y.real.reshape(-1, 1),
                    D_lam_sig_tau_y.imag.reshape(-1, 1),
                    D_lam_sig_tau_z.real.reshape(-1, 1),
                    D_lam_sig_tau_z.imag.reshape(-1, 1),
                ))

            else:
                dist_den_1_freq = None
                dist_den_2_freq = None
                dist_den_3_freq = None

            dist_den_1_freq = DistributedArray(dist_den_1_freq, self.comm)
            dist_den_2_freq = DistributedArray(dist_den_2_freq, self.comm)
            dist_den_3_freq = DistributedArray(dist_den_3_freq, self.comm)

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

            if distributed_density_3 is None:
                distributed_density_3 = DistributedArray(dist_den_3_freq.data,
                                                         self.comm,
                                                         distribute=False)
            else:
                distributed_density_3.append(dist_den_3_freq, axis=1)

        return distributed_density_1, distributed_density_2, distributed_density_3

    def get_fock_dict(self, wi, density_list1, density_list2, density_list3,
                      F0_a, mo, molecule, ao_basis, dft_dict, profiler):
        """
        Computes the compounded Fock matrices F^{σ},F^{λ+τ},F^{σλτ} used for the
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

        for wb in wi:
            for key in [
                    'f_sig_xx', 'f_sig_yy', 'f_sig_zz', 'f_sig_xy', 'f_sig_xz',
                    'f_sig_yz', 'f_lamtau_xx', 'f_lamtau_yy', 'f_lamtau_zz',
                    'f_lamtau_xy', 'f_lamtau_xz', 'f_lamtau_yz'
            ]:
                key_freq_pairs.append((key, wb))

        for wb in wi:
            for key in ['F123_x', 'F123_y', 'F123_z']:
                key_freq_pairs.append((key, wb))

        # examine checkpoint file for distributed Focks

        if self.checkpoint_file is not None:
            fock_file = str(
                Path(self.checkpoint_file).with_suffix('.tpa_fock_1_full.h5'))
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
                                                 density_list3, 'tpa', profiler)
            else:
                density_list_23 = DistributedArray(density_list2.data,
                                                   self.comm,
                                                   distribute=False)
                density_list_23.append(density_list3, axis=1)
                dist_focks = self._comp_nlr_fock(mo, molecule, ao_basis,
                                                 'real_and_imag', None, None,
                                                 None, density_list_23, 'tpa',
                                                 profiler)

            self._print_fock_time(time.time() - time_start_fock)

            write_distributed_focks(fock_file, dist_focks, key_freq_pairs,
                                    self.comm, self.ostream)

        # assign distributed Focks to key-frequency pairs

        focks = {'F0': F0_a}

        for fock_index, (key, wb) in enumerate(key_freq_pairs):
            if key not in focks:
                focks[key] = {}
            focks[key][wb] = DistributedArray(dist_focks.data[:, fock_index],
                                              self.comm,
                                              distribute=False)

        return focks

    def get_e4(self, wi, Nx, fo, nocc, norb):
        """
        Contracts E[4]NxNyNz for the isotropic cubic response function. Takes
        the Fock matrices from fock_dict and contracts them with the response
        vectors.

        :param wi:
            A list of freqs
        :param Nx:
            A dict of the single index response vectors in distributed form
        :param fo:
            A dictonary of transformed Fock matricies from fock_dict
        :param nocc:
            The number of occupied orbitals
        :param norb:
            The total number of orbitals

        :return:
            A dictonary of compounded E[4] tensors for the isotropic cubic
            response function for TPA
        """

        f_iso_x = {}
        f_iso_y = {}
        f_iso_z = {}

        for w in wi:

            vec_pack = np.array([
                fo['Fb'][('x', w)].data,
                fo['Fb'][('y', w)].data,
                fo['Fb'][('z', w)].data,
                fo['f_lamtau_xx'][w].data,
                fo['f_lamtau_yy'][w].data,
                fo['f_lamtau_zz'][w].data,
                fo['f_lamtau_xy'][w].data,
                fo['f_lamtau_xz'][w].data,
                fo['f_lamtau_yz'][w].data,
                fo['f_sig_xx'][w].data,
                fo['f_sig_yy'][w].data,
                fo['f_sig_zz'][w].data,
                fo['f_sig_xy'][w].data,
                fo['f_sig_xz'][w].data,
                fo['f_sig_yz'][w].data,
                fo['F123_x'][w].data,
                fo['F123_y'][w].data,
                fo['F123_z'][w].data,
            ]).T.copy()

            vec_pack = self._collect_vectors_in_columns(vec_pack)

            nx = ComplexResponse.get_full_solution_vector(Nx[('x', w)])
            ny = ComplexResponse.get_full_solution_vector(Nx[('y', w)])
            nz = ComplexResponse.get_full_solution_vector(Nx[('z', w)])

            if self.rank != mpi_master():
                continue

            vec_pack = vec_pack.T.copy().reshape(-1, norb, norb)

            (Fx, Fy, Fz, f_lamtau_xx, f_lamtau_yy, f_lamtau_zz, f_lamtau_xy,
             f_lamtau_xz, f_lamtau_yz, f_sig_xx, f_sig_yy, f_sig_zz, f_sig_xy,
             f_sig_xz, f_sig_yz, f_x, f_y, f_z) = vec_pack

            Fx_ = Fx.T.conj()  # fo['Fc'][('x', -w)]
            Fy_ = Fy.T.conj()  # fo['Fc'][('y', -w)]
            Fz_ = Fz.T.conj()  # fo['Fc'][('z', -w)]

            F0 = fo['F0']

            # Get all the response matrices and Fock matrices

            kx = (self.complex_lrvec2mat(nx, nocc, norb)).T
            ky = (self.complex_lrvec2mat(ny, nocc, norb)).T
            kz = (self.complex_lrvec2mat(nz, nocc, norb)).T

            kx_ = -kx.conj().T  # kX['Nc'][('x', -w)].T
            ky_ = -ky.conj().T  # kX['Nc'][('y', -w)].T
            kz_ = -kz.conj().T  # kX['Nc'][('z', -w)].T

            # computes all the compounded Φ_αβ, see article, where small phi
            # here is defined as:
            #   φ(κa,κb,Fb,F0) = [κa,[κb,F0]+3Fb]
            #   Φ_αα^σ = φ(κα,κα,Fα,F0) + φ(κα,κα,Fα,F0) +
            #            Σ_{ρ}^{x,y,z}[φ(κρ,κρ,Fρ,F0)] (for α=β)
            #   Φ_αβ^σ = φ(κa,κa,Fb,F0) + φ(κb,κb,Fb,F0) (for α≠β)
            # For the Φ_{αβ}^{λ+τ} component see article.

            Pxx = self.phi(kx, kx, Fx, F0)
            Pyy = self.phi(ky, ky, Fy, F0)
            Pzz = self.phi(kz, kz, Fz, F0)

            Phi_sig_xx = 2 * (3 * Pxx + Pyy + Pzz)
            Phi_sig_yy = 2 * (Pxx + 3 * Pyy + Pzz)
            Phi_sig_zz = 2 * (Pxx + Pyy + 3 * Pzz)

            Phi_sig_xy = (2 * self.phi(kx, ky, Fy, F0) +
                          2 * self.phi(ky, kx, Fx, F0))

            Phi_sig_xz = (2 * self.phi(kx, kz, Fz, F0) +
                          2 * self.phi(kz, kx, Fx, F0))

            Phi_sig_yz = (2 * self.phi(ky, kz, Fz, F0) +
                          2 * self.phi(kz, ky, Fy, F0))

            Pxx = self.phi(kx, kx_, Fx_, F0) + self.phi(kx_, kx, Fx, F0)
            Pyy = self.phi(ky, ky_, Fy_, F0) + self.phi(ky_, ky, Fy, F0)
            Pzz = self.phi(kz, kz_, Fz_, F0) + self.phi(kz_, kz, Fz, F0)

            Phi_lamtau_xx = 2 * (3 * Pxx + Pyy + Pzz)
            Phi_lamtau_yy = 2 * (Pxx + 3 * Pyy + Pzz)
            Phi_lamtau_zz = 2 * (Pxx + Pyy + 3 * Pzz)

            Phi_lamtau_xy = 2 * (
                self.phi(kx_, ky, Fy, F0) + self.phi(ky, kx_, Fx_, F0) +
                self.phi(kx, ky_, Fy_, F0) + self.phi(ky_, kx, Fx, F0))

            Phi_lamtau_xz = 2 * (
                self.phi(kx_, kz, Fz, F0) + self.phi(kz, kx_, Fx_, F0) +
                self.phi(kx, kz_, Fz_, F0) + self.phi(kz_, kx, Fx, F0))

            Phi_lamtau_yz = 2 * (
                self.phi(ky_, kz, Fz, F0) + self.phi(kz, ky_, Fy_, F0) +
                self.phi(ky, kz_, Fz_, F0) + self.phi(kz_, ky, Fy, F0))

            # Computess all the elements of the Fock vector formed from the
            # E[4] contraction as E[4]NxNyNz = [f_{is} // f_{si}], for the
            # isotropic case, as derived in the article
            # The elements of f_{is} for each spatial component α is given by
            # an expression of the form
            # f_α = Σ_{β}^{x,y,z} [κ_{β}^{ω},Φ_{αβ}^{λ+τ}+f_{αβ}^{λ+τ}] +
            #       [κ_{β}^{-ω},Φ_{αβ}^{σ}+f_{αβ}^{σ}]

            # x

            # Creating the transformed total Fock matrices
            f_x += (self.commut(kx, Phi_lamtau_xx + f_lamtau_xx) +
                    self.commut(ky, Phi_lamtau_xy + f_lamtau_xy) +
                    self.commut(kz, Phi_lamtau_xz + f_lamtau_xz))
            f_x += (self.commut(kx_, Phi_sig_xx + f_sig_xx) +
                    self.commut(ky_, Phi_sig_xy + f_sig_xy) +
                    self.commut(kz_, Phi_sig_xz + f_sig_xz))

            # Taking the non redundant matrix elements {i,s} and forming the
            # anti-symmetric Fock vector
            f_x = -2. / 6 * LinearSolver.lrmat2vec(f_x.T, nocc, norb)
            f_x = self.anti_sym(f_x)
            f_iso_x[w] = f_x

            # y

            # Creating the transformed total Fock matrices
            f_y += (self.commut(kx, Phi_lamtau_xy + f_lamtau_xy) +
                    self.commut(ky, Phi_lamtau_yy + f_lamtau_yy) +
                    self.commut(kz, Phi_lamtau_yz + f_lamtau_yz))
            f_y += (self.commut(kx_, Phi_sig_xy + f_sig_xy) +
                    self.commut(ky_, Phi_sig_yy + f_sig_yy) +
                    self.commut(kz_, Phi_sig_yz + f_sig_yz))

            # Taking the non redundant matrix elements {i,s} and forming the
            # anti-symmetric Fock vector
            f_y = -2. / 6 * LinearSolver.lrmat2vec(f_y.T, nocc, norb)
            f_y = self.anti_sym(f_y)
            f_iso_y[w] = f_y

            # z

            # Creating the transformed total Fock matrices
            f_z += (self.commut(kx, Phi_lamtau_xz + f_lamtau_xz) +
                    self.commut(ky, Phi_lamtau_yz + f_lamtau_yz) +
                    self.commut(kz, Phi_lamtau_zz + f_lamtau_zz))
            f_z += (self.commut(kx_, Phi_sig_xz + f_sig_xz) +
                    self.commut(ky_, Phi_sig_yz + f_sig_yz) +
                    self.commut(kz_, Phi_sig_zz + f_sig_zz))

            # Taking the non redundant matrix elements {i,s} and forming the
            # anti-symmetric Fock vector
            f_z = -2. / 6 * LinearSolver.lrmat2vec(f_z.T, nocc, norb)
            f_z = self.anti_sym(f_z)
            f_iso_z[w] = f_z

        return {'f_iso_x': f_iso_x, 'f_iso_y': f_iso_y, 'f_iso_z': f_iso_z}

    def get_Nxy(self, w, d_a_mo, X, fock_dict, Nx, nocc, norb, molecule,
                ao_basis, scf_tensors):
        """
        Computes all the second-order response vectors needed for the isotropic
        cubic response computation

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

        # Get second-order gradients
        xy_dict = self.get_xy(d_a_mo, X, w, fock_dict, Nx, nocc, norb)

        # Frequencies to compute
        if self.rank == mpi_master():
            freq = tuple([0.0] + [sum(x) for x in zip(w, w)])
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
                Path(self.checkpoint_file).with_suffix('.tpa_2_full.h5'))

        # commutpute second-order response vectors

        N_total_results = N_total_drv.compute(molecule, ao_basis, scf_tensors,
                                              xy_dict)

        self._is_converged = (self._is_converged and N_total_drv.is_converged)

        Nxy_dict = N_total_results['solutions']
        FXY_2_dict = N_total_results['focks']

        return (Nxy_dict, FXY_2_dict)

    def get_xy(self, d_a_mo, X, wi, Fock, Nx, nocc, norb):
        """
        Computes the compounded gradient vectors N^{σ},N^{λ+τ} used for the
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
            A dictonary with all the first-order response vectors in
            distributed form
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
                Fock['f_lamtau_xx'][w].data * one_third,
                Fock['f_lamtau_yy'][w].data * one_third,
                Fock['f_lamtau_zz'][w].data * one_third,
                Fock['f_lamtau_xy'][w].data * one_third,
                Fock['f_lamtau_xz'][w].data * one_third,
                Fock['f_lamtau_yz'][w].data * one_third,
            ]).T.copy()

            vec_pack = self._collect_vectors_in_columns(vec_pack)

            nx = ComplexResponse.get_full_solution_vector(Nx[('x', w)])
            ny = ComplexResponse.get_full_solution_vector(Nx[('y', w)])
            nz = ComplexResponse.get_full_solution_vector(Nx[('z', w)])

            if self.rank != mpi_master():
                continue

            vec_pack = vec_pack.T.copy().reshape(-1, norb, norb)

            (f_x, f_y, f_z, f_sig_xx, f_sig_yy, f_sig_zz, f_sig_xy, f_sig_xz,
             f_sig_yz, f_lamtau_xx, f_lamtau_yy, f_lamtau_zz, f_lamtau_xy,
             f_lamtau_xz, f_lamtau_yz) = vec_pack

            f_x_ = f_x.T.conj()  # Fock['Fc'][('x', -w)]
            f_y_ = f_y.T.conj()  # Fock['Fc'][('y', -w)]
            f_z_ = f_z.T.conj()  # Fock['Fc'][('z', -w)]

            mu_x = X['x']
            mu_y = X['y']
            mu_z = X['z']

            kx = (self.complex_lrvec2mat(nx, nocc, norb)).T
            ky = (self.complex_lrvec2mat(ny, nocc, norb)).T
            kz = (self.complex_lrvec2mat(nz, nocc, norb)).T

            kx_ = -kx.conj().T  # kX['Nc'][('x', -w)].T
            ky_ = -ky.conj().T  # kX['Nc'][('y', -w)].T
            kz_ = -kz.conj().T  # kX['Nc'][('z', -w)].T

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

            # BC CD λ+τ gradients #

            xi_xx = self._xi(kx_, kx, f_x_, f_x, F0)
            xi_yy = self._xi(ky_, ky, f_y_, f_y, F0)
            xi_zz = self._xi(kz_, kz, f_z_, f_z, F0)

            x2_xx_ = self._x2_contract(kx_.T, mu_x, d_a_mo, nocc, norb)
            x2_yy_ = self._x2_contract(ky_.T, mu_y, d_a_mo, nocc, norb)
            x2_zz_ = self._x2_contract(kz_.T, mu_z, d_a_mo, nocc, norb)

            key = (('N_lamtau_xx', w), 0)
            mat = (6 * xi_xx + 2 * xi_yy + 2 * xi_zz + 0.5 * f_lamtau_xx).T
            xy_dict[key] = self.anti_sym(
                -LinearSolver.lrmat2vec(mat, nocc, norb))
            xy_dict[key] -= (3 * x2_xx_ + x2_yy_ + x2_zz_)
            xy_dict[key] -= (3 * x2_xx + x2_yy + x2_zz)

            key = (('N_lamtau_yy', w), 0)
            mat = (2 * xi_xx + 6 * xi_yy + 2 * xi_zz + 0.5 * f_lamtau_yy).T
            xy_dict[key] = self.anti_sym(
                -LinearSolver.lrmat2vec(mat, nocc, norb))
            xy_dict[key] -= (x2_xx_ + 3 * x2_yy_ + x2_zz_)
            xy_dict[key] -= (x2_xx + 3 * x2_yy + x2_zz)

            key = (('N_lamtau_zz', w), 0)
            mat = (2 * xi_xx + 2 * xi_yy + 6 * xi_zz + 0.5 * f_lamtau_zz).T
            xy_dict[key] = self.anti_sym(
                -LinearSolver.lrmat2vec(mat, nocc, norb))
            xy_dict[key] -= (x2_xx_ + x2_yy_ + 3 * x2_zz_)
            xy_dict[key] -= (x2_xx + x2_yy + 3 * x2_zz)

            key = (('N_lamtau_xy', w), 0)
            mat = (2 * self._xi(ky_, kx, f_y_, f_x, F0) +
                   2 * self._xi(kx_, ky, f_x_, f_y, F0) + 0.5 * f_lamtau_xy).T
            xy_dict[key] = self.anti_sym(
                -LinearSolver.lrmat2vec(mat, nocc, norb))
            xy_dict[key] -= self._x2_contract(ky.T, mu_x, d_a_mo, nocc, norb)
            xy_dict[key] -= self._x2_contract(kx.T, mu_y, d_a_mo, nocc, norb)
            xy_dict[key] -= self._x2_contract(ky_.T, mu_x, d_a_mo, nocc, norb)
            xy_dict[key] -= self._x2_contract(kx_.T, mu_y, d_a_mo, nocc, norb)

            key = (('N_lamtau_xz', w), 0)
            mat = (2 * self._xi(kz_, kx, f_z_, f_x, F0) +
                   2 * self._xi(kx_, kz, f_x_, f_z, F0) + 0.5 * f_lamtau_xz).T
            xy_dict[key] = self.anti_sym(
                -LinearSolver.lrmat2vec(mat, nocc, norb))
            xy_dict[key] -= self._x2_contract(kz.T, mu_x, d_a_mo, nocc, norb)
            xy_dict[key] -= self._x2_contract(kx.T, mu_z, d_a_mo, nocc, norb)
            xy_dict[key] -= self._x2_contract(kz_.T, mu_x, d_a_mo, nocc, norb)
            xy_dict[key] -= self._x2_contract(kx_.T, mu_z, d_a_mo, nocc, norb)

            key = (('N_lamtau_yz', w), 0)
            mat = (2 * self._xi(kz_, ky, f_z_, f_y, F0) +
                   2 * self._xi(ky_, kz, f_y_, f_z, F0) + 0.5 * f_lamtau_yz).T
            xy_dict[key] = self.anti_sym(
                -LinearSolver.lrmat2vec(mat, nocc, norb))
            xy_dict[key] -= self._x2_contract(kz.T, mu_y, d_a_mo, nocc, norb)
            xy_dict[key] -= self._x2_contract(ky.T, mu_z, d_a_mo, nocc, norb)
            xy_dict[key] -= self._x2_contract(kz_.T, mu_y, d_a_mo, nocc, norb)
            xy_dict[key] -= self._x2_contract(ky_.T, mu_z, d_a_mo, nocc, norb)

        return xy_dict

    def get_densities_II(self, wi, Nx, Nxy, mo, nocc, norb):
        """
        Computes the compounded densities needed for the compounded
        second-order Fock matrices used for the isotropic cubic response
        function

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

            n_lamtau_xx = ComplexResponse.get_full_solution_vector(
                Nxy[(('N_lamtau_xx', w), 0)])
            n_lamtau_yy = ComplexResponse.get_full_solution_vector(
                Nxy[(('N_lamtau_yy', w), 0)])
            n_lamtau_zz = ComplexResponse.get_full_solution_vector(
                Nxy[(('N_lamtau_zz', w), 0)])
            n_lamtau_xy = ComplexResponse.get_full_solution_vector(
                Nxy[(('N_lamtau_xy', w), 0)])
            n_lamtau_xz = ComplexResponse.get_full_solution_vector(
                Nxy[(('N_lamtau_xz', w), 0)])
            n_lamtau_yz = ComplexResponse.get_full_solution_vector(
                Nxy[(('N_lamtau_yz', w), 0)])

            if self.rank == mpi_master():

                k_sig_xx = self.complex_lrvec2mat(n_sig_xx, nocc, norb)
                k_sig_yy = self.complex_lrvec2mat(n_sig_yy, nocc, norb)
                k_sig_zz = self.complex_lrvec2mat(n_sig_zz, nocc, norb)
                k_sig_xy = self.complex_lrvec2mat(n_sig_xy, nocc, norb)
                k_sig_xz = self.complex_lrvec2mat(n_sig_xz, nocc, norb)
                k_sig_yz = self.complex_lrvec2mat(n_sig_yz, nocc, norb)

                k_lamtau_xx = self.complex_lrvec2mat(n_lamtau_xx, nocc, norb)
                k_lamtau_yy = self.complex_lrvec2mat(n_lamtau_yy, nocc, norb)
                k_lamtau_zz = self.complex_lrvec2mat(n_lamtau_zz, nocc, norb)
                k_lamtau_xy = self.complex_lrvec2mat(n_lamtau_xy, nocc, norb)
                k_lamtau_xz = self.complex_lrvec2mat(n_lamtau_xz, nocc, norb)
                k_lamtau_yz = self.complex_lrvec2mat(n_lamtau_yz, nocc, norb)

                kx = self.complex_lrvec2mat(nx, nocc, norb)
                ky = self.complex_lrvec2mat(ny, nocc, norb)
                kz = self.complex_lrvec2mat(nz, nocc, norb)

                kx_ = -kx.conj().T  # kX['Nc'][('x', -w)].T
                ky_ = -ky.conj().T  # kX['Nc'][('y', -w)].T
                kz_ = -kz.conj().T  # kX['Nc'][('z', -w)].T

                # SIGMA contributiatons #
                Dc_x = self.commut_mo_density(kx_, nocc)
                Dc_y = self.commut_mo_density(ky_, nocc)
                Dc_z = self.commut_mo_density(kz_, nocc)

                D_sig_xx = self.commut_mo_density(k_sig_xx, nocc)
                D_sig_yy = self.commut_mo_density(k_sig_yy, nocc)
                D_sig_zz = self.commut_mo_density(k_sig_zz, nocc)

                D_sig_xy = self.commut_mo_density(k_sig_xy, nocc)
                D_sig_xz = self.commut_mo_density(k_sig_xz, nocc)
                D_sig_yz = self.commut_mo_density(k_sig_yz, nocc)

                # LAM + TAU contributiatons #
                Db_x = self.commut_mo_density(kx, nocc)
                Db_y = self.commut_mo_density(ky, nocc)
                Db_z = self.commut_mo_density(kz, nocc)

                D_lamtau_xx = self.commut_mo_density(k_lamtau_xx, nocc)
                D_lamtau_yy = self.commut_mo_density(k_lamtau_yy, nocc)
                D_lamtau_zz = self.commut_mo_density(k_lamtau_zz, nocc)

                D_lamtau_xy = self.commut_mo_density(k_lamtau_xy, nocc)
                D_lamtau_xz = self.commut_mo_density(k_lamtau_xz, nocc)
                D_lamtau_yz = self.commut_mo_density(k_lamtau_yz, nocc)

                # x #
                Dx = self.commut(kx_, D_sig_xx)
                Dx += self.commut(k_sig_xx, Dc_x)
                Dx += self.commut(ky_, D_sig_xy)
                Dx += self.commut(k_sig_xy, Dc_y)

                Dx += self.commut(kz_, D_sig_xz)
                Dx += self.commut(k_sig_xz, Dc_z)

                Dx += self.commut(kx, D_lamtau_xx)
                Dx += self.commut(k_lamtau_xx, Db_x)

                Dx += self.commut(ky, D_lamtau_xy)
                Dx += self.commut(k_lamtau_xy, Db_y)

                Dx += self.commut(kz, D_lamtau_xz)
                Dx += self.commut(k_lamtau_xz, Db_z)

                # y #
                Dy = self.commut(kx_, D_sig_xy)
                Dy += self.commut(k_sig_xy, Dc_x)

                Dy += self.commut(ky_, D_sig_yy)
                Dy += self.commut(k_sig_yy, Dc_y)

                Dy += self.commut(kz_, D_sig_yz)
                Dy += self.commut(k_sig_yz, Dc_z)

                Dy += self.commut(kx, D_lamtau_xy)
                Dy += self.commut(k_lamtau_xy, Db_x)

                Dy += self.commut(ky, D_lamtau_yy)
                Dy += self.commut(k_lamtau_yy, Db_y)

                Dy += self.commut(kz, D_lamtau_yz)
                Dy += self.commut(k_lamtau_yz, Db_z)

                # z #
                Dz = self.commut(kx_, D_sig_xz)
                Dz += self.commut(k_sig_xz, Dc_x)

                Dz += self.commut(ky_, D_sig_yz)
                Dz += self.commut(k_sig_yz, Dc_y)

                Dz += self.commut(kz_, D_sig_zz)
                Dz += self.commut(k_sig_zz, Dc_z)

                Dz += self.commut(kx, D_lamtau_xz)
                Dz += self.commut(k_lamtau_xz, Db_x)

                Dz += self.commut(ky, D_lamtau_yz)
                Dz += self.commut(k_lamtau_yz, Db_y)

                Dz += self.commut(kz, D_lamtau_zz)
                Dz += self.commut(k_lamtau_zz, Db_z)

                # density transformation from MO to AO basis

                Dx = np.linalg.multi_dot([mo, Dx, mo.T])
                Dy = np.linalg.multi_dot([mo, Dy, mo.T])
                Dz = np.linalg.multi_dot([mo, Dz, mo.T])

                Db_x = np.linalg.multi_dot([mo, Db_x, mo.T])
                Db_y = np.linalg.multi_dot([mo, Db_y, mo.T])
                Db_z = np.linalg.multi_dot([mo, Db_z, mo.T])

                Dc_x = np.linalg.multi_dot([mo, Dc_x, mo.T])
                Dc_y = np.linalg.multi_dot([mo, Dc_y, mo.T])
                Dc_z = np.linalg.multi_dot([mo, Dc_z, mo.T])

                D_sig_xx = np.linalg.multi_dot([mo, D_sig_xx, mo.T])
                D_sig_yy = np.linalg.multi_dot([mo, D_sig_yy, mo.T])
                D_sig_zz = np.linalg.multi_dot([mo, D_sig_zz, mo.T])
                D_sig_xy = np.linalg.multi_dot([mo, D_sig_xy, mo.T])
                D_sig_xz = np.linalg.multi_dot([mo, D_sig_xz, mo.T])
                D_sig_yz = np.linalg.multi_dot([mo, D_sig_yz, mo.T])

                D_lamtau_xx = np.linalg.multi_dot([mo, D_lamtau_xx, mo.T])
                D_lamtau_yy = np.linalg.multi_dot([mo, D_lamtau_yy, mo.T])
                D_lamtau_zz = np.linalg.multi_dot([mo, D_lamtau_zz, mo.T])
                D_lamtau_xy = np.linalg.multi_dot([mo, D_lamtau_xy, mo.T])
                D_lamtau_xz = np.linalg.multi_dot([mo, D_lamtau_xz, mo.T])
                D_lamtau_yz = np.linalg.multi_dot([mo, D_lamtau_yz, mo.T])

                dist_den_1_freq = np.hstack((
                    Db_x.real.reshape(-1, 1),
                    Db_x.imag.reshape(-1, 1),
                    Db_y.real.reshape(-1, 1),
                    Db_y.imag.reshape(-1, 1),
                    Db_z.real.reshape(-1, 1),
                    Db_z.imag.reshape(-1, 1),
                    Dc_x.real.reshape(-1, 1),
                    Dc_x.imag.reshape(-1, 1),
                    Dc_y.real.reshape(-1, 1),
                    Dc_y.imag.reshape(-1, 1),
                    Dc_z.real.reshape(-1, 1),
                    Dc_z.imag.reshape(-1, 1),
                    D_sig_xx.real.reshape(-1, 1),
                    D_sig_xx.imag.reshape(-1, 1),
                    D_sig_yy.real.reshape(-1, 1),
                    D_sig_yy.imag.reshape(-1, 1),
                    D_sig_zz.real.reshape(-1, 1),
                    D_sig_zz.imag.reshape(-1, 1),
                    D_sig_xy.real.reshape(-1, 1),
                    D_sig_xy.imag.reshape(-1, 1),
                    D_sig_xz.real.reshape(-1, 1),
                    D_sig_xz.imag.reshape(-1, 1),
                    D_sig_yz.real.reshape(-1, 1),
                    D_sig_yz.imag.reshape(-1, 1),
                    D_lamtau_xx.real.reshape(-1, 1),
                    D_lamtau_xx.imag.reshape(-1, 1),
                    D_lamtau_yy.real.reshape(-1, 1),
                    D_lamtau_yy.imag.reshape(-1, 1),
                    D_lamtau_zz.real.reshape(-1, 1),
                    D_lamtau_zz.imag.reshape(-1, 1),
                    D_lamtau_xy.real.reshape(-1, 1),
                    D_lamtau_xy.imag.reshape(-1, 1),
                    D_lamtau_xz.real.reshape(-1, 1),
                    D_lamtau_xz.imag.reshape(-1, 1),
                    D_lamtau_yz.real.reshape(-1, 1),
                    D_lamtau_yz.imag.reshape(-1, 1),
                ))

                dist_den_2_freq = np.hstack((
                    Dx.real.reshape(-1, 1),
                    Dx.imag.reshape(-1, 1),
                    Dy.real.reshape(-1, 1),
                    Dy.imag.reshape(-1, 1),
                    Dz.real.reshape(-1, 1),
                    Dz.imag.reshape(-1, 1),
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

        # generate key-frequency paris

        key_freq_pairs = []
        for w in wi:
            for key in ['F123_x', 'F123_y', 'F123_z']:
                key_freq_pairs.append((key, w))

        # examine checkpoint file for distributed Focks

        if self.checkpoint_file is not None:
            fock_file = str(
                Path(self.checkpoint_file).with_suffix('.tpa_fock_2_full.h5'))
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
                                                 None, 'tpa_ii', profiler)
            else:
                dist_focks = self._comp_nlr_fock(mo, molecule, ao_basis,
                                                 'real_and_imag', None, None,
                                                 density_list2, None, 'tpa_ii',
                                                 profiler)

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
        Contracts E[3]

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
            A dictonary of compounded E[3] tensors for the isotropic cubic
            response function for TPA
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
                fo2[(('N_lamtau_xx', w), 0)].data,
                fo2[(('N_lamtau_yy', w), 0)].data,
                fo2[(('N_lamtau_zz', w), 0)].data,
                fo2[(('N_lamtau_xy', w), 0)].data,
                fo2[(('N_lamtau_xz', w), 0)].data,
                fo2[(('N_lamtau_yz', w), 0)].data,
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

            n_lamtau_xx = ComplexResponse.get_full_solution_vector(
                Nxy[(('N_lamtau_xx', w), 0)])
            n_lamtau_yy = ComplexResponse.get_full_solution_vector(
                Nxy[(('N_lamtau_yy', w), 0)])
            n_lamtau_zz = ComplexResponse.get_full_solution_vector(
                Nxy[(('N_lamtau_zz', w), 0)])
            n_lamtau_xy = ComplexResponse.get_full_solution_vector(
                Nxy[(('N_lamtau_xy', w), 0)])
            n_lamtau_xz = ComplexResponse.get_full_solution_vector(
                Nxy[(('N_lamtau_xz', w), 0)])
            n_lamtau_yz = ComplexResponse.get_full_solution_vector(
                Nxy[(('N_lamtau_yz', w), 0)])

            if self.rank != mpi_master():
                continue

            vec_pack = vec_pack.T.copy().reshape(-1, norb, norb)

            (f_x, f_y, f_z, f_sig_xx, f_sig_yy, f_sig_zz, f_sig_xy, f_sig_xz,
             f_sig_yz, f_lamtau_xx, f_lamtau_yy, f_lamtau_zz, f_lamtau_xy,
             f_lamtau_xz, f_lamtau_yz, F123_x, F123_y, F123_z) = vec_pack

            f_x_ = f_x.T.conj()  # fo['Fc'][('x', -w)]
            f_y_ = f_y.T.conj()  # fo['Fc'][('y', -w)]
            f_z_ = f_z.T.conj()  # fo['Fc'][('z', -w)]

            f_sig_xx = f_sig_xx.T.conj()
            f_sig_yy = f_sig_yy.T.conj()
            f_sig_zz = f_sig_zz.T.conj()

            f_sig_xy = f_sig_xy.T.conj()
            f_sig_xz = f_sig_xz.T.conj()
            f_sig_yz = f_sig_yz.T.conj()

            f_lamtau_xx = f_lamtau_xx.T.conj()
            f_lamtau_yy = f_lamtau_yy.T.conj()
            f_lamtau_zz = f_lamtau_zz.T.conj()

            f_lamtau_xy = f_lamtau_xy.T.conj()
            f_lamtau_xz = f_lamtau_xz.T.conj()
            f_lamtau_yz = f_lamtau_yz.T.conj()

            F0_a = fo['F0']

            # Response

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

            k_lamtau_xx = (self.complex_lrvec2mat(n_lamtau_xx, nocc, norb)).T
            k_lamtau_yy = (self.complex_lrvec2mat(n_lamtau_yy, nocc, norb)).T
            k_lamtau_zz = (self.complex_lrvec2mat(n_lamtau_zz, nocc, norb)).T
            k_lamtau_xy = (self.complex_lrvec2mat(n_lamtau_xy, nocc, norb)).T
            k_lamtau_xz = (self.complex_lrvec2mat(n_lamtau_xz, nocc, norb)).T
            k_lamtau_yz = (self.complex_lrvec2mat(n_lamtau_yz, nocc, norb)).T

            # Focks #

            # x

            zeta_sig_xx = self._xi(k_x_, k_sig_xx, f_x_, f_sig_xx, F0_a)
            zeta_sig_yy = self._xi(k_x_, k_sig_yy, f_x_, f_sig_yy, F0_a)
            zeta_sig_zz = self._xi(k_x_, k_sig_zz, f_x_, f_sig_zz, F0_a)

            zeta_sig_xy = self._xi(k_y_, k_sig_xy, f_y_, f_sig_xy, F0_a)
            zeta_sig_xz = self._xi(k_z_, k_sig_xz, f_z_, f_sig_xz, F0_a)

            zeta_lamtau_xx = self._xi(k_x, k_lamtau_xx, f_x, f_lamtau_xx, F0_a)
            zeta_lamtau_yy = self._xi(k_x, k_lamtau_yy, f_x, f_lamtau_yy, F0_a)
            zeta_lamtau_zz = self._xi(k_x, k_lamtau_zz, f_x, f_lamtau_zz, F0_a)

            zeta_lamtau_xy = self._xi(k_y, k_lamtau_xy, f_y, f_lamtau_xy, F0_a)
            zeta_lamtau_xz = self._xi(k_z, k_lamtau_xz, f_z, f_lamtau_xz, F0_a)

            X_terms = (zeta_sig_xx + zeta_sig_xy +
                       zeta_sig_xz).T + (zeta_lamtau_xx + zeta_lamtau_xy +
                                         zeta_lamtau_xz).T + (0.5 * F123_x).T
            Ff_x = -2 * LinearSolver.lrmat2vec(X_terms, nocc, norb)
            Ff_x = self.anti_sym(Ff_x)
            f_iso_x[w] = Ff_x

            # y

            zeta_sig_yx = self._xi(k_x_, k_sig_xy, f_x_, f_sig_xy, F0_a)
            zeta_sig_yy = self._xi(k_y_, k_sig_yy, f_y_, f_sig_yy, F0_a)
            zeta_sig_yz = self._xi(k_z_, k_sig_yz, f_z_, f_sig_yz, F0_a)

            zeta_lamtau_yx = self._xi(k_x, k_lamtau_xy, f_x, f_lamtau_xy, F0_a)
            zeta_lamtau_yy = self._xi(k_y, k_lamtau_yy, f_y, f_lamtau_yy, F0_a)
            zeta_lamtau_yz = self._xi(k_z, k_lamtau_yz, f_z, f_lamtau_yz, F0_a)

            Y_terms = (zeta_sig_yx + zeta_sig_yy +
                       zeta_sig_yz).T + (zeta_lamtau_yx + zeta_lamtau_yy +
                                         zeta_lamtau_yz).T + (0.5 * F123_y).T
            Ff_y = -2 * LinearSolver.lrmat2vec(Y_terms, nocc, norb)
            Ff_y = self.anti_sym(Ff_y)
            f_iso_y[w] = Ff_y

            # z

            zeta_sig_zx = self._xi(k_x_, k_sig_xz, f_x_, f_sig_xz, F0_a)
            zeta_sig_zy = self._xi(k_y_, k_sig_yz, f_y_, f_sig_yz, F0_a)
            zeta_sig_zz = self._xi(k_z_, k_sig_zz, f_z_, f_sig_zz, F0_a)

            zeta_lamtau_zx = self._xi(k_x, k_lamtau_xz, f_x, f_lamtau_xz, F0_a)
            zeta_lamtau_zy = self._xi(k_y, k_lamtau_yz, f_y, f_lamtau_yz, F0_a)
            zeta_lamtau_zz = self._xi(k_z, k_lamtau_zz, f_z, f_lamtau_zz, F0_a)

            Z_terms = (zeta_sig_zx + zeta_sig_zy +
                       zeta_sig_zz).T + (zeta_lamtau_zx + zeta_lamtau_zy +
                                         zeta_lamtau_zz).T + (0.5 * F123_z).T
            Ff_z = -2 * LinearSolver.lrmat2vec(Z_terms, nocc, norb)
            Ff_z = self.anti_sym(Ff_z)
            f_iso_z[w] = Ff_z

        return {'f_iso_x': f_iso_x, 'f_iso_y': f_iso_y, 'f_iso_z': f_iso_z}

    def get_other_terms(self, wi, track, X, Nx, Nxy, da, nocc, norb):
        """
        Computes the terms involving X[3],A[3],X[2],A[2] in the isotropic cubic
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

        na_a3_nx_ny_dict = {}
        na_x3_ny_nz_dict = {}
        na_x2_nyz_dict = {}
        nx_a2_nyz_dict = {}

        comp_per_freq = len(track) // len(wi)

        inp_list = []

        for j in range(len(wi)):
            for i in range(j * comp_per_freq, (j + 1) * comp_per_freq):
                comp_i = track[i]

                vals = comp_i.split(',')
                op_a, op_b, op_c, op_d = vals[0]
                w = float(vals[1])
                wa = float(vals[1])
                wb = float(vals[1])
                wc = float(vals[2])
                wd = float(vals[3])

                inp_list.append({
                    'freq': w,
                    'Na': Nx[(op_a, wa)],
                    'Nb': Nx[(op_b, wb)],
                    'Nc': Nx[(op_c, -wc)],
                    'Nd': Nx[(op_d, wd)],
                    'A': X[op_a],
                    'B': X[op_b],
                    'C': X[op_c],
                    'D': X[op_d],
                })

        list_x3_a3 = [self.get_x3_a3(inp, da, nocc, norb) for inp in inp_list]

        if self.rank == mpi_master():
            for terms in list_x3_a3:
                if terms['key'] not in na_x3_ny_nz_dict:
                    na_x3_ny_nz_dict[terms['key']] = 0.0
                if terms['key'] not in na_a3_nx_ny_dict:
                    na_a3_nx_ny_dict[terms['key']] = 0.0
                na_x3_ny_nz_dict[terms['key']] += terms['x3']
                na_a3_nx_ny_dict[terms['key']] += terms['a3']

        inp_list = []

        for i in range(len(wi)):
            vals = track[i * comp_per_freq].split(',')
            w = float(vals[1])
            wa = float(vals[1])
            wb = float(vals[1])
            wc = float(vals[2])
            wd = float(vals[3])

            wcd = 0
            wbd = wb + wd

            for op_a in 'xyz':
                Na = Nx[(op_a, wa)]
                A = X[op_a]

                for op_b in 'xyz':
                    op_ab = op_a + op_b if op_a <= op_b else op_b + op_a

                    # CD
                    Ncd = Nxy[(('N_lamtau_' + op_ab, w), wcd)]
                    Nb = Nx[(op_b, w)]
                    B = X[op_b]

                    inp_list.append({
                        'flag': 'CD',
                        'freq': w,
                        'Ncd': Ncd,
                        'Na': Na,
                        'Nb': Nb,
                        'A': A,
                        'B': B,
                    })

                    # BD
                    op_c = op_b
                    op_ac = op_ab
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
                'NaX3NyNz': na_x3_ny_nz_dict,
                'NaA3NxNy': na_a3_nx_ny_dict,
                'NaX2Nyz': na_x2_nyz_dict,
                'NxA2Nyz': nx_a2_nyz_dict,
            }

        return None

    def get_x3_a3(self, inp_dict, da, nocc, norb):
        """
        Computes X[3] and A[3] contributions.

        :param inp_dict:
            A dictionary containing input data for computing X[3] and A[3].
        :param da:
            The SCF density matrix in MO basis
        :param nocc:
            The number of occupied orbitals
        :param norb:
            The total number of orbitals

        :return:
            A dictionary containing frequencies, X[3] and A[3].
        """

        na_x3_ny_nz = 0.0
        na_a3_nx_ny = 0.0

        w = inp_dict['freq']

        Na = ComplexResponse.get_full_solution_vector(inp_dict['Na'])
        Nb = ComplexResponse.get_full_solution_vector(inp_dict['Nb'])
        Nc = ComplexResponse.get_full_solution_vector(inp_dict['Nc'])
        Nd = ComplexResponse.get_full_solution_vector(inp_dict['Nd'])

        if self.rank == mpi_master():
            kb = self.complex_lrvec2mat(Nb, nocc, norb)
            kc = -(self.complex_lrvec2mat(Nc, nocc, norb)).T.conj()
            kd = self.complex_lrvec2mat(Nd, nocc, norb)

            Nc = self.flip_yz(Nc)  # gets Nc from Nb

            A = inp_dict['A']
            B = inp_dict['B']
            C = inp_dict['C']
            D = inp_dict['D']

            # Na X[3]NyNz

            na_x3_ny_nz -= np.dot(Na.T,
                                  self._x3_contract(kc, kd, B, da, nocc, norb))
            na_x3_ny_nz -= np.dot(Na.T,
                                  self._x3_contract(kd, kc, B, da, nocc, norb))
            na_x3_ny_nz -= np.dot(Na.T,
                                  self._x3_contract(kd, kb, C, da, nocc, norb))
            na_x3_ny_nz -= np.dot(Na.T,
                                  self._x3_contract(kb, kd, C, da, nocc, norb))
            na_x3_ny_nz -= np.dot(Na.T,
                                  self._x3_contract(kb, kc, D, da, nocc, norb))
            na_x3_ny_nz -= np.dot(Na.T,
                                  self._x3_contract(kc, kb, D, da, nocc, norb))

            # NaA[3]NxNy

            na_a3_nx_ny += np.dot(self._a3_contract(kb, kc, A, da, nocc, norb),
                                  Nd)
            na_a3_nx_ny += np.dot(self._a3_contract(kb, kd, A, da, nocc, norb),
                                  Nc)
            na_a3_nx_ny += np.dot(self._a3_contract(kc, kb, A, da, nocc, norb),
                                  Nd)
            na_a3_nx_ny += np.dot(self._a3_contract(kc, kd, A, da, nocc, norb),
                                  Nb)
            na_a3_nx_ny += np.dot(self._a3_contract(kd, kb, A, da, nocc, norb),
                                  Nc)
            na_a3_nx_ny += np.dot(self._a3_contract(kd, kc, A, da, nocc, norb),
                                  Nb)

            return {
                'key': (w, -w, w),
                'x3': (1. / 15) * na_x3_ny_nz,
                'a3': (1. / 15) * na_a3_nx_ny,
            }
        else:
            return None

    def get_t4(self, wi, e4_dict, Nx, track, da, nocc, norb):
        """
        Computes the contraction of the E[4] tensor with that of the S[4] and
        R[4] tensors to return the contraction of T[4] as a dictonary of
        vectors. T[4]NxNyNz = (E^[4]-ω_1S^[4]-ω_1S^[4]-ω_3S^[4]-γiR^[4])

        :param wi:
            A list of all the freqs
        :param e4_dict:
            A dictonary of all the E[4] contraction
        :param Nx:
            A dictonray containng all the response vectors in distributed form
        :param track:
            A list containg information about all the γ components that are to
            be computed
        :param da:
            The SCF density matrix in MO basis
        :param nocc:
            The number of occupied orbitals
        :param norb:
            The total number of orbitals

        :return:
            A dictonary of final T[4] contraction values
        """

        T4 = {}
        S4, R4 = self.get_s4_and_r4(wi, Nx, track, da, nocc, norb)

        comp_per_freq = len(track) // len(wi)

        for i in range(len(wi)):
            vals = track[i * comp_per_freq].split(',')
            w = float(vals[1])
            ww = float(vals[1])

            na_x = ComplexResponse.get_full_solution_vector(Nx[('x', w)])
            na_y = ComplexResponse.get_full_solution_vector(Nx[('y', w)])
            na_z = ComplexResponse.get_full_solution_vector(Nx[('z', w)])

            if self.rank == mpi_master():
                t4val = (np.dot(na_x, e4_dict['f_iso_x'][ww] - S4[('x', ww)]) +
                         np.dot(na_y, e4_dict['f_iso_y'][ww] - S4[('y', ww)]) +
                         np.dot(na_z, e4_dict['f_iso_z'][ww] - S4[('z', ww)]))

                if self.damping > 0:
                    t4val += (R4[('x', ww)] + R4[('y', ww)] + R4[('z', ww)])

                T4[(ww, -ww, ww)] = -(1. / 15) * t4val

        if self.rank == mpi_master():
            return T4
        else:
            return None

    def get_s4_and_r4(self, wi, Nx, track, D0, nocc, norb):
        """
        Computes the S4 contractions

        :param wi:
            A list of all the freqs
        :param Nx:
            A dict with all the response vectors in distributed form
        :param track:
            A list containing information about all the components that are to
            be computed
        :param D0:
            The SCF density in MO basis
        :param nocc:
            The number of occupied obritals
        :param norb:
            The number of total orbitals

        :return:
            A dictonary of final S[4] contraction values
        """

        S4 = {}
        R4 = {}

        comp_per_freq = len(track) // len(wi)

        inp_list = []

        for j in range(len(wi)):
            vals = track[j * comp_per_freq].split(',')
            w = float(vals[1])
            w1 = float(vals[1])
            w2 = float(vals[2])
            w3 = float(vals[3])
            w_s = w1 + w2 + w3

            for i in range(j * comp_per_freq, (j + 1) * comp_per_freq):
                comp_i = track[i]
                op = comp_i[0]

                inp_dict = {
                    'w': w,
                    'w1': w1,
                    'w2': w2,
                    'w3': w3,
                    'op': op,
                    'Nb': Nx[(comp_i[1], w1)],
                    'Nc': Nx[(comp_i[2], -w2)],
                    'Nd': Nx[(comp_i[3], w3)],
                }

                if self.damping > 0:
                    inp_dict.update({
                        'Nb': Nx[(comp_i[1], w1)],
                        'Nc': Nx[(comp_i[2], -w2)],
                        'Nd': Nx[(comp_i[3], w3)],
                        'Na': Nx[(comp_i[0], w_s)],
                    })

                inp_list.append(inp_dict)

        list_s4_r4 = [
            self.get_s4_and_r4_terms(inp, D0, nocc, norb) for inp in inp_list
        ]

        if self.rank == mpi_master():
            local_s4_dict = {}
            local_r4_dict = {}

            for terms in list_s4_r4:
                if terms['s4_key'] not in local_s4_dict:
                    local_s4_dict[terms['s4_key']] = 0.0
                local_s4_dict[terms['s4_key']] += terms['s4']

                if terms['r4_key'] not in local_r4_dict:
                    local_r4_dict[terms['r4_key']] = 0.0
                local_r4_dict[terms['r4_key']] += terms['r4']

            list_s4_r4 = []

            for s4_key, r4_key in zip(list(local_s4_dict.keys()),
                                      list(local_r4_dict.keys())):
                list_s4_r4.append({
                    's4_key': s4_key,
                    's4': local_s4_dict[s4_key],
                    'r4_key': r4_key,
                    'r4': local_r4_dict[r4_key],
                })

            for terms in list_s4_r4:
                s4_key = terms['s4_key']
                r4_key = terms['r4_key']
                if s4_key not in S4:
                    S4[s4_key] = 0.0
                if r4_key not in R4:
                    R4[r4_key] = 0.0
                S4[s4_key] += terms['s4']
                R4[r4_key] += terms['r4']

            return S4, R4
        else:
            return None, None

    def get_s4_and_r4_terms(self, inp_dict, D0, nocc, norb):
        """
        Computes S[4] and R[4] contributions.

        :param inp_dict:
            A dictionary containing input data for computing S[4] and R[4].
        :param D0:
            The SCF density matrix in MO basis
        :param nocc:
            The number of occupied orbitals
        :param norb:
            The total number of orbitals

        :return:
            Dictionaries containing S[4] and R[4].
        """

        s4_term = 0.0
        r4_term = 0.0

        w = inp_dict['w']
        w1 = inp_dict['w1']
        w2 = inp_dict['w2']
        w3 = inp_dict['w3']
        op = inp_dict['op']

        s4_key = (op, w)
        r4_key = (op, w1)

        Nb = ComplexResponse.get_full_solution_vector(inp_dict['Nb'])
        Nc = ComplexResponse.get_full_solution_vector(inp_dict['Nc'])
        Nd = ComplexResponse.get_full_solution_vector(inp_dict['Nd'])

        if self.damping > 0:
            Na = ComplexResponse.get_full_solution_vector(inp_dict['Na'])

        if self.rank == mpi_master():
            kB = self.complex_lrvec2mat(Nb, nocc, norb)
            kC = -(self.complex_lrvec2mat(Nc, nocc, norb)).T.conj()
            kD = self.complex_lrvec2mat(Nd, nocc, norb)

            s4_term -= w1 * self._s4(kB, kC, kD, D0, nocc, norb)
            s4_term -= w2 * self._s4(kC, kB, kD, D0, nocc, norb)
            s4_term -= w3 * self._s4(kD, kB, kC, D0, nocc, norb)

            if self.damping > 0:
                kA = self.complex_lrvec2mat(Na, nocc, norb)

                Nc = self.flip_yz(Nc)
                Nb_h = self.flip_xy(Nb)
                Nc_h = self.flip_xy(Nc)
                Nd_h = self.flip_xy(Nd)

                r4_term += 1j * self.damping * np.dot(
                    Nd_h, self._s4_for_r4(kA.T, kB, kC, D0, nocc, norb))
                r4_term += 1j * self.damping * np.dot(
                    Nc_h, self._s4_for_r4(kA.T, kB, kD, D0, nocc, norb))
                r4_term += 1j * self.damping * np.dot(
                    Nd_h, self._s4_for_r4(kA.T, kC, kB, D0, nocc, norb))
                r4_term += 1j * self.damping * np.dot(
                    Nb_h, self._s4_for_r4(kA.T, kC, kD, D0, nocc, norb))
                r4_term += 1j * self.damping * np.dot(
                    Nc_h, self._s4_for_r4(kA.T, kD, kB, D0, nocc, norb))
                r4_term += 1j * self.damping * np.dot(
                    Nb_h, self._s4_for_r4(kA.T, kD, kC, D0, nocc, norb))

            return {
                's4_key': s4_key,
                'r4_key': r4_key,
                's4': s4_term,
                'r4': r4_term,
            }
        else:
            return None

    def print_results(self, freqs, gamma, comp, t4_dict, t3_dict, tpa_dict):
        """
        Prints the results from the TPA calculation.

        :param freqs:
            List of frequencies
        :param gamma:
            A dictonary containing the isotropic cubic response functions for
            TPA
        :param comp:
            List of gamma tensors components
        :param t4_dict:
            A dictonary containing the isotropic T[4] contractions
        :param t3_dict:
            A dictonary containing the isotropic T[3] contractions
        :param tpa_dict:
            A dictonary containing the isotropic X[3], A[3], X[2], A[2]
            contractions
        """

        NaX3NyNz = tpa_dict['NaX3NyNz']
        NaA3NxNy = tpa_dict['NaA3NxNy']
        NaX2Nyz = tpa_dict['NaX2Nyz']
        NxA2Nyz = tpa_dict['NxA2Nyz']

        self.ostream.print_blank()

        w_str = 'The Isotropic Average gamma Tensor and Its'
        self.ostream.print_header(w_str)
        w_str = 'Isotropic Contributions at Given Frequencies'
        self.ostream.print_header(w_str)
        self.ostream.print_header('=' * (len(w_str) + 2))
        self.ostream.print_blank()

        for w in freqs:
            title = '{:<9s} {:>12s} {:>20s} {:>21s}'.format(
                'Component', 'Frequency', 'Real', 'Imaginary')
            width = len(title)
            self.ostream.print_header(title.ljust(width))
            self.ostream.print_header(('-' * len(title)).ljust(width))

            self._print_component('T3', w, t3_dict[w, -w, w], width)
            self._print_component('T4', w, t4_dict[w, -w, w], width)
            self._print_component('X2', w, NaX2Nyz[w, -w, w], width)
            self._print_component('X3', w, NaX3NyNz[w, -w, w], width)
            self._print_component('A2', w, NxA2Nyz[w, -w, w], width)
            self._print_component('A3', w, NaA3NxNy[w, -w, w], width)
            self._print_component('gamma', w, gamma[w, -w, w], width)

            self.ostream.print_blank()

        title = 'Reference: '
        title += 'K. Ahmadzadeh, M. Scott, M. Brand, O. Vahtras, X. Li, '
        self.ostream.print_header(title.ljust(width))
        title = 'Z. Rinkevicius, and P. Norman, '
        title += 'J. Chem. Phys. 154, 024111 (2021)'
        self.ostream.print_header(title.ljust(width))
        self.ostream.print_blank()

        spectrum = self.get_spectrum(gamma)
        self.print_spectrum(spectrum, width)

        self.ostream.print_blank()
