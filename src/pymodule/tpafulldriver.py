import numpy as np
import time
import re

from .veloxchemlib import mpi_master
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

    def __init__(self, comm, ostream):
        """
        Initializes the full isotropic cubic response driver for two-photon
        absorption (TPA)
        """

        super().__init__(comm, ostream)

    def update_settings(self, rsp_dict, method_dict=None):
        """
        Updates response and method settings for TPA

        :param rsp_dict:
            The dictionary of response dict.
        :param method_dict:
            The dictionary of method rsp_dict.
        """

        if method_dict is None:
            method_dict = {}

        super().update_settings(rsp_dict, method_dict)

    def get_densities(self, wi, kX, S, D0, mo):
        """
        Computes the compounded densities needed for the compounded Fock
        matrics F^{σ},F^{λ+τ},F^{σλτ} used for the isotropic cubic response
        function. Note: All densities are 1/3 of those in the paper, and all
        the Fock matrices are later scaled by 3.

        :param wi:
            A list of the frequencies
        :param kX:
            A dictonary with all the first-order response matrices
        :param S:
            The overlap matrix
        :param D0:
            The SCF density matrix in AO basis
        :param mo:
            A matrix containing the MO coefficents

        :return:
            A list of tranformed compounded densities
        """

        density_list = []

        for w in wi:

            # convert response matrix to ao basis #

            kx = self.mo2ao(mo, kX['Nb'][('x', w)]).T
            ky = self.mo2ao(mo, kX['Nb'][('y', w)]).T
            kz = self.mo2ao(mo, kX['Nb'][('z', w)]).T

            kx_ = -kx.conj().T  # self.mo2ao(mo, kX['Nc'][('x', -w)]).T
            ky_ = -ky.conj().T  # self.mo2ao(mo, kX['Nc'][('y', -w)]).T
            kz_ = -kz.conj().T  # self.mo2ao(mo, kX['Nc'][('z', -w)]).T

            # create the first order single indexed densiteies #

            Dx = self.transform_dens(kx, D0, S)
            Dy = self.transform_dens(ky, D0, S)
            Dz = self.transform_dens(kz, D0, S)

            Dx_ = self.transform_dens(kx_, D0, S)
            Dy_ = self.transform_dens(ky_, D0, S)
            Dz_ = self.transform_dens(kz_, D0, S)

            # create the first order two indexed densities #

            # σ terms #

            Dxx = self.transform_dens(kx, Dx, S)
            Dyy = self.transform_dens(ky, Dy, S)
            Dzz = self.transform_dens(kz, Dz, S)

            D_sig_xx = 2 * (3 * Dxx + Dyy + Dzz)
            D_sig_yy = 2 * (Dxx + 3 * Dyy + Dzz)
            D_sig_zz = 2 * (Dxx + Dyy + 3 * Dzz)

            D_sig_xy = 2 * (self.transform_dens(ky, Dx, S) +
                            self.transform_dens(kx, Dy, S))
            D_sig_xz = 2 * (self.transform_dens(kx, Dz, S) +
                            self.transform_dens(kz, Dx, S))
            D_sig_yz = 2 * (self.transform_dens(ky, Dz, S) +
                            self.transform_dens(kz, Dy, S))

            # λ+τ terms #

            Dxx = (self.transform_dens(kx_, Dx, S) +
                   self.transform_dens(kx, Dx_, S))
            Dyy = (self.transform_dens(ky_, Dy, S) +
                   self.transform_dens(ky, Dy_, S))
            Dzz = (self.transform_dens(kz_, Dz, S) +
                   self.transform_dens(kz, Dz_, S))

            D_lamtau_xx = 2 * (3 * Dxx + Dyy + Dzz)
            D_lamtau_yy = 2 * (Dxx + 3 * Dyy + Dzz)
            D_lamtau_zz = 2 * (Dxx + Dyy + 3 * Dzz)

            D_lamtau_xy = 2 * (self.transform_dens(ky_, Dx, S) +
                               self.transform_dens(kx_, Dy, S))
            D_lamtau_xy += 2 * (self.transform_dens(ky, Dx_, S) +
                                self.transform_dens(kx, Dy_, S))

            D_lamtau_xz = 2 * (self.transform_dens(kx_, Dz, S) +
                               self.transform_dens(kz_, Dx, S))
            D_lamtau_xz += 2 * (self.transform_dens(kx, Dz_, S) +
                                self.transform_dens(kz, Dx_, S))

            D_lamtau_yz = 2 * (self.transform_dens(ky_, Dz, S) +
                               self.transform_dens(kz_, Dy, S))
            D_lamtau_yz += 2 * (self.transform_dens(ky, Dz_, S) +
                                self.transform_dens(kz, Dy_, S))

            # Create first order three indexed Densities #

            D_lam_sig_tau_x = (self.transform_dens(kx_, D_sig_xx, S) +
                               self.transform_dens(ky_, D_sig_xy, S) +
                               self.transform_dens(kz_, D_sig_xz, S))
            D_lam_sig_tau_x += (self.transform_dens(kx, D_lamtau_xx, S) +
                                self.transform_dens(ky, D_lamtau_xy, S) +
                                self.transform_dens(kz, D_lamtau_xz, S))

            D_lam_sig_tau_y = (self.transform_dens(kx_, D_sig_xy, S) +
                               self.transform_dens(ky_, D_sig_yy, S) +
                               self.transform_dens(kz_, D_sig_yz, S))
            D_lam_sig_tau_y += (self.transform_dens(kx, D_lamtau_xy, S) +
                                self.transform_dens(ky, D_lamtau_yy, S) +
                                self.transform_dens(kz, D_lamtau_yz, S))

            D_lam_sig_tau_z = (self.transform_dens(kx_, D_sig_xz, S) +
                               self.transform_dens(ky_, D_sig_yz, S) +
                               self.transform_dens(kz_, D_sig_zz, S))
            D_lam_sig_tau_z += (self.transform_dens(kx, D_lamtau_xz, S) +
                                self.transform_dens(ky, D_lamtau_yz, S) +
                                self.transform_dens(kz, D_lamtau_zz, S))

            density_list.append(D_sig_xx)
            density_list.append(D_sig_yy)
            density_list.append(D_sig_zz)
            density_list.append(D_sig_xy)
            density_list.append(D_sig_xz)
            density_list.append(D_sig_yz)

            density_list.append(D_lamtau_xx)
            density_list.append(D_lamtau_yy)
            density_list.append(D_lamtau_zz)
            density_list.append(D_lamtau_xy)
            density_list.append(D_lamtau_xz)
            density_list.append(D_lamtau_yz)

            density_list.append(D_lam_sig_tau_x)
            density_list.append(D_lam_sig_tau_y)
            density_list.append(D_lam_sig_tau_z)

        return density_list

    def get_fock_dict(self, wi, density_list, F0_a, mo, molecule, ao_basis):
        """
        Computes the compounded Fock matrics F^{σ},F^{λ+τ},F^{σλτ} used for the
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
            'f_lamtau_xx',
            'f_lamtau_yy',
            'f_lamtau_zz',
            'f_lamtau_xy',
            'f_lamtau_xz',
            'f_lamtau_yz',
            'F123_x',
            'F123_y',
            'F123_z',
        ]

        if self.checkpoint_file is not None:
            fock_file = re.sub(r'\.h5$', r'', self.checkpoint_file)
            fock_file += '_tpa_fock_1_full.h5'
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
        dist_focks = self.get_fock_r(mo, density_list, molecule, ao_basis,
                                     'real_and_imag')
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

    def get_e4(self, wi, kX, fo, nocc, norb):
        """
        Contracts E[4]NxNyNz for the isotropic cubic response function. Takes
        the Fock matrices from fock_dict and contracts them with the response
        vectors.

        :param wi:
            A list of freqs
        :param kX:
            A dict of the single index response matricies
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
                fo['f_lamtau_xx'][w].data * 3,
                fo['f_lamtau_yy'][w].data * 3,
                fo['f_lamtau_zz'][w].data * 3,
                fo['f_lamtau_xy'][w].data * 3,
                fo['f_lamtau_xz'][w].data * 3,
                fo['f_lamtau_yz'][w].data * 3,
                fo['f_sig_xx'][w].data * 3,
                fo['f_sig_yy'][w].data * 3,
                fo['f_sig_zz'][w].data * 3,
                fo['f_sig_xy'][w].data * 3,
                fo['f_sig_xz'][w].data * 3,
                fo['f_sig_yz'][w].data * 3,
                fo['F123_x'][w].data,
                fo['F123_y'][w].data,
                fo['F123_z'][w].data,
            ]).T.copy()

            vec_pack = self.collect_vectors_in_columns(vec_pack)

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

            kx = kX['Nb'][('x', w)].T
            ky = kX['Nb'][('y', w)].T
            kz = kX['Nb'][('z', w)].T

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

    def get_Nxy(self, w, d_a_mo, X, fock_dict, kX, nocc, norb, molecule,
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

        xy_dict = {}
        freq = None

        # Get second-order gradients
        xy_dict = self.get_xy(d_a_mo, X, w, fock_dict, kX, nocc, norb)

        if self.rank == mpi_master():

            # Frequencies to compute
            wbd = [sum(x) for x in zip(w, w)]

            freq = '0.0,'
            for i in range(len(wbd)):
                freq += str(wbd[i]) + ','

        freq = self.comm.bcast(freq, root=mpi_master())

        N_total_drv = ComplexResponse(self.comm, self.ostream)

        N_total_drv.update_settings({
            'frequencies': freq,
            'damping': self.damping,
            'conv_thresh': self.conv_thresh,
            'lindep_thresh': self.lindep_thresh,
            'max_iter': self.max_iter,
            'eri_thresh': self.eri_thresh,
            'qq_type': self.qq_type,
        })
        N_total_drv.timing = self.timing
        N_total_drv.memory_profiling = self.memory_profiling
        N_total_drv.batch_size = self.batch_size
        N_total_drv.restart = self.restart
        N_total_drv.program_start_time = self.program_start_time
        N_total_drv.maximum_hours = self.maximum_hours
        if self.checkpoint_file is not None:
            N_total_drv.checkpoint_file = re.sub(r'\.h5$', r'',
                                                 self.checkpoint_file)
            N_total_drv.checkpoint_file += '_tpa_2_full.h5'

        # commutpute second-order response vectors
        N_total_results = N_total_drv.compute(molecule, ao_basis, scf_tensors,
                                              xy_dict)

        Nxy_dict = N_total_results['solutions']
        kXY_dict = N_total_results['kappas']
        FXY_2_dict = N_total_results['focks']

        return (Nxy_dict, kXY_dict, FXY_2_dict)

    def get_xy(self, d_a_mo, X, wi, Fock, kX, nocc, norb):
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
                Fock['f_lamtau_xx'][w].data,
                Fock['f_lamtau_yy'][w].data,
                Fock['f_lamtau_zz'][w].data,
                Fock['f_lamtau_xy'][w].data,
                Fock['f_lamtau_xz'][w].data,
                Fock['f_lamtau_yz'][w].data,
            ]).T.copy()

            vec_pack = self.collect_vectors_in_columns(vec_pack)

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

            kx = kX['Nb'][('x', w)].T
            ky = kX['Nb'][('y', w)].T
            kz = kX['Nb'][('z', w)].T

            kx_ = -kx.conj().T  # kX['Nc'][('x', -w)].T
            ky_ = -ky.conj().T  # kX['Nc'][('y', -w)].T
            kz_ = -kz.conj().T  # kX['Nc'][('z', -w)].T

            F0 = Fock['F0']

            # BD σ gradients #

            xi_xx = self.xi(kx, kx, f_x, f_x, F0)
            xi_yy = self.xi(ky, ky, f_y, f_y, F0)
            xi_zz = self.xi(kz, kz, f_z, f_z, F0)

            key = (('N_sig_xx', w), 2 * w)
            mat = (3 * xi_xx + xi_yy + xi_zz + 0.5 * f_sig_xx).T
            xy_dict[key] = self.anti_sym(
                -2 * LinearSolver.lrmat2vec(mat, nocc, norb))
            xy_dict[key] -= 6 * self.x2_contract(kx.T, mu_x, d_a_mo, nocc, norb)
            xy_dict[key] -= 2 * self.x2_contract(ky.T, mu_y, d_a_mo, nocc, norb)
            xy_dict[key] -= 2 * self.x2_contract(kz.T, mu_z, d_a_mo, nocc, norb)

            key = (('N_sig_yy', w), 2 * w)
            mat = (xi_xx + 3 * xi_yy + xi_zz + 0.5 * f_sig_yy).T
            xy_dict[key] = self.anti_sym(
                -2 * LinearSolver.lrmat2vec(mat, nocc, norb))
            xy_dict[key] -= 2 * self.x2_contract(kx.T, mu_x, d_a_mo, nocc, norb)
            xy_dict[key] -= 6 * self.x2_contract(ky.T, mu_y, d_a_mo, nocc, norb)
            xy_dict[key] -= 2 * self.x2_contract(kz.T, mu_z, d_a_mo, nocc, norb)

            key = (('N_sig_zz', w), 2 * w)
            mat = (xi_xx + xi_yy + 3 * xi_zz + 0.5 * f_sig_zz).T
            xy_dict[key] = self.anti_sym(
                -2 * LinearSolver.lrmat2vec(mat, nocc, norb))
            xy_dict[key] -= 2 * self.x2_contract(kx.T, mu_x, d_a_mo, nocc, norb)
            xy_dict[key] -= 2 * self.x2_contract(ky.T, mu_y, d_a_mo, nocc, norb)
            xy_dict[key] -= 6 * self.x2_contract(kz.T, mu_z, d_a_mo, nocc, norb)

            key = (('N_sig_xy', w), 2 * w)
            mat = (self.xi(ky, kx, f_y, f_x, F0) +
                   self.xi(kx, ky, f_x, f_y, F0) + 0.5 * f_sig_xy).T
            xy_dict[key] = self.anti_sym(
                -2 * LinearSolver.lrmat2vec(mat, nocc, norb))
            xy_dict[key] -= 2 * self.x2_contract(ky.T, mu_x, d_a_mo, nocc, norb)
            xy_dict[key] -= 2 * self.x2_contract(kx.T, mu_y, d_a_mo, nocc, norb)

            key = (('N_sig_xz', w), 2 * w)
            mat = (self.xi(kz, kx, f_z, f_x, F0) +
                   self.xi(kx, kz, f_x, f_z, F0) + 0.5 * f_sig_xz).T
            xy_dict[key] = self.anti_sym(
                -2 * LinearSolver.lrmat2vec(mat, nocc, norb))
            xy_dict[key] -= 2 * self.x2_contract(kz.T, mu_x, d_a_mo, nocc, norb)
            xy_dict[key] -= 2 * self.x2_contract(kx.T, mu_z, d_a_mo, nocc, norb)

            key = (('N_sig_yz', w), 2 * w)
            mat = (self.xi(kz, ky, f_z, f_y, F0) +
                   self.xi(ky, kz, f_y, f_z, F0) + 0.5 * f_sig_yz).T
            xy_dict[key] = self.anti_sym(
                -2 * LinearSolver.lrmat2vec(mat, nocc, norb))
            xy_dict[key] -= 2 * self.x2_contract(kz.T, mu_y, d_a_mo, nocc, norb)
            xy_dict[key] -= 2 * self.x2_contract(ky.T, mu_z, d_a_mo, nocc, norb)

            # BC CD λ+τ gradients #

            xi_xx = self.xi(kx_, kx, f_x_, f_x, F0)
            xi_yy = self.xi(ky_, ky, f_y_, f_y, F0)
            xi_zz = self.xi(kz_, kz, f_z_, f_z, F0)

            key = (('N_lamtau_xx', w), 0)
            mat = (6 * xi_xx + 2 * xi_yy + 2 * xi_zz + 0.5 * f_lamtau_xx).T
            xy_dict[key] = self.anti_sym(
                -2 * LinearSolver.lrmat2vec(mat, nocc, norb))
            xy_dict[key] -= 6 * self.x2_contract(kx_.T, mu_x, d_a_mo, nocc,
                                                 norb)
            xy_dict[key] -= 2 * self.x2_contract(ky_.T, mu_y, d_a_mo, nocc,
                                                 norb)
            xy_dict[key] -= 2 * self.x2_contract(kz_.T, mu_z, d_a_mo, nocc,
                                                 norb)
            xy_dict[key] -= 6 * self.x2_contract(kx.T, mu_x, d_a_mo, nocc, norb)
            xy_dict[key] -= 2 * self.x2_contract(ky.T, mu_y, d_a_mo, nocc, norb)
            xy_dict[key] -= 2 * self.x2_contract(kz.T, mu_z, d_a_mo, nocc, norb)

            key = (('N_lamtau_yy', w), 0)
            mat = (2 * xi_xx + 6 * xi_yy + 2 * xi_zz + 0.5 * f_lamtau_yy).T
            xy_dict[key] = self.anti_sym(
                -2 * LinearSolver.lrmat2vec(mat, nocc, norb))
            xy_dict[key] -= 2 * self.x2_contract(kx_.T, mu_x, d_a_mo, nocc,
                                                 norb)
            xy_dict[key] -= 6 * self.x2_contract(ky_.T, mu_y, d_a_mo, nocc,
                                                 norb)
            xy_dict[key] -= 2 * self.x2_contract(kz_.T, mu_z, d_a_mo, nocc,
                                                 norb)
            xy_dict[key] -= 2 * self.x2_contract(kx.T, mu_x, d_a_mo, nocc, norb)
            xy_dict[key] -= 6 * self.x2_contract(ky.T, mu_y, d_a_mo, nocc, norb)
            xy_dict[key] -= 2 * self.x2_contract(kz.T, mu_z, d_a_mo, nocc, norb)

            key = (('N_lamtau_zz', w), 0)
            mat = (2 * xi_xx + 2 * xi_yy + 6 * xi_zz + 0.5 * f_lamtau_zz).T
            xy_dict[key] = self.anti_sym(
                -2 * LinearSolver.lrmat2vec(mat, nocc, norb))
            xy_dict[key] -= 2 * self.x2_contract(kx_.T, mu_x, d_a_mo, nocc,
                                                 norb)
            xy_dict[key] -= 2 * self.x2_contract(ky_.T, mu_y, d_a_mo, nocc,
                                                 norb)
            xy_dict[key] -= 6 * self.x2_contract(kz_.T, mu_z, d_a_mo, nocc,
                                                 norb)
            xy_dict[key] -= 2 * self.x2_contract(kx.T, mu_x, d_a_mo, nocc, norb)
            xy_dict[key] -= 2 * self.x2_contract(ky.T, mu_y, d_a_mo, nocc, norb)
            xy_dict[key] -= 6 * self.x2_contract(kz.T, mu_z, d_a_mo, nocc, norb)

            key = (('N_lamtau_xy', w), 0)
            mat = (2 * self.xi(ky_, kx, f_y_, f_x, F0) +
                   2 * self.xi(kx_, ky, f_x_, f_y, F0) + 0.5 * f_lamtau_xy).T
            xy_dict[key] = self.anti_sym(
                -2 * LinearSolver.lrmat2vec(mat, nocc, norb))
            xy_dict[key] -= 2 * self.x2_contract(ky.T, mu_x, d_a_mo, nocc, norb)
            xy_dict[key] -= 2 * self.x2_contract(kx.T, mu_y, d_a_mo, nocc, norb)
            xy_dict[key] -= 2 * self.x2_contract(ky_.T, mu_x, d_a_mo, nocc,
                                                 norb)
            xy_dict[key] -= 2 * self.x2_contract(kx_.T, mu_y, d_a_mo, nocc,
                                                 norb)

            key = (('N_lamtau_xz', w), 0)
            mat = (2 * self.xi(kz_, kx, f_z_, f_x, F0) +
                   2 * self.xi(kx_, kz, f_x_, f_z, F0) + 0.5 * f_lamtau_xz).T
            xy_dict[key] = self.anti_sym(
                -2 * LinearSolver.lrmat2vec(mat, nocc, norb))
            xy_dict[key] -= 2 * self.x2_contract(kz.T, mu_x, d_a_mo, nocc, norb)
            xy_dict[key] -= 2 * self.x2_contract(kx.T, mu_z, d_a_mo, nocc, norb)
            xy_dict[key] -= 2 * self.x2_contract(kz_.T, mu_x, d_a_mo, nocc,
                                                 norb)
            xy_dict[key] -= 2 * self.x2_contract(kx_.T, mu_z, d_a_mo, nocc,
                                                 norb)

            key = (('N_lamtau_yz', w), 0)
            mat = (2 * self.xi(kz_, ky, f_z_, f_y, F0) +
                   2 * self.xi(ky_, kz, f_y_, f_z, F0) + 0.5 * f_lamtau_yz).T
            xy_dict[key] = self.anti_sym(
                -2 * LinearSolver.lrmat2vec(mat, nocc, norb))
            xy_dict[key] -= 2 * self.x2_contract(kz.T, mu_y, d_a_mo, nocc, norb)
            xy_dict[key] -= 2 * self.x2_contract(ky.T, mu_z, d_a_mo, nocc, norb)
            xy_dict[key] -= 2 * self.x2_contract(kz_.T, mu_y, d_a_mo, nocc,
                                                 norb)
            xy_dict[key] -= 2 * self.x2_contract(ky_.T, mu_z, d_a_mo, nocc,
                                                 norb)

        return xy_dict

    def get_densities_II(self, wi, kX, kXY, S, D0, mo):
        """
        Computes the compounded densities needed for the compounded
        second-order Fock matrics used for the isotropic cubic response
        function

        :param wi:
            A list of the frequencies
        :param kX:
            A dictonary with all the first-order response matrices
        :param kXY:
            A dict of the two index response matrices
        :param S:
            The overlap matrix
        :param D0:
            The SCF density matrix in AO basis
        :param mo:
            A matrix containing the MO coefficents

        :return:
            A list of tranformed compounded densities
        """

        density_list = []

        for w in wi:
            k_sig_xx = self.mo2ao(mo, kXY[(('N_sig_xx', w), 2 * w)]).T
            k_sig_yy = self.mo2ao(mo, kXY[(('N_sig_yy', w), 2 * w)]).T
            k_sig_zz = self.mo2ao(mo, kXY[(('N_sig_zz', w), 2 * w)]).T

            k_sig_xy = self.mo2ao(mo, kXY[(('N_sig_xy', w), 2 * w)]).T
            k_sig_xz = self.mo2ao(mo, kXY[(('N_sig_xz', w), 2 * w)]).T
            k_sig_yz = self.mo2ao(mo, kXY[(('N_sig_yz', w), 2 * w)]).T

            k_lamtau_xx = self.mo2ao(mo, kXY[(('N_lamtau_xx', w), 0)]).T
            k_lamtau_yy = self.mo2ao(mo, kXY[(('N_lamtau_yy', w), 0)]).T
            k_lamtau_zz = self.mo2ao(mo, kXY[(('N_lamtau_zz', w), 0)]).T

            k_lamtau_xy = self.mo2ao(mo, kXY[(('N_lamtau_xy', w), 0)]).T
            k_lamtau_xz = self.mo2ao(mo, kXY[(('N_lamtau_xz', w), 0)]).T
            k_lamtau_yz = self.mo2ao(mo, kXY[(('N_lamtau_yz', w), 0)]).T

            kx = self.mo2ao(mo, kX['Nb'][('x', w)]).T
            ky = self.mo2ao(mo, kX['Nb'][('y', w)]).T
            kz = self.mo2ao(mo, kX['Nb'][('z', w)]).T

            kx_ = -kx.conj().T  # self.mo2ao(mo, kX['Nc'][('x', -w)]).T
            ky_ = -ky.conj().T  # self.mo2ao(mo, kX['Nc'][('y', -w)]).T
            kz_ = -kz.conj().T  # self.mo2ao(mo, kX['Nc'][('z', -w)]).T

            # SIGMA contributiatons #
            Dc_x_ = self.transform_dens(kx_, D0, S)
            Dc_y_ = self.transform_dens(ky_, D0, S)
            Dc_z_ = self.transform_dens(kz_, D0, S)

            D_sig_xx = self.transform_dens(k_sig_xx, D0, S)
            D_sig_yy = self.transform_dens(k_sig_yy, D0, S)
            D_sig_zz = self.transform_dens(k_sig_zz, D0, S)

            D_sig_xy = self.transform_dens(k_sig_xy, D0, S)
            D_sig_xz = self.transform_dens(k_sig_xz, D0, S)
            D_sig_yz = self.transform_dens(k_sig_yz, D0, S)

            # LAM + TAU contributiatons #
            Dc_x = self.transform_dens(kx, D0, S)
            Dc_y = self.transform_dens(ky, D0, S)
            Dc_z = self.transform_dens(kz, D0, S)

            D_lamtau_xx = self.transform_dens(k_lamtau_xx, D0, S)
            D_lamtau_yy = self.transform_dens(k_lamtau_yy, D0, S)
            D_lamtau_zz = self.transform_dens(k_lamtau_zz, D0, S)

            D_lamtau_xy = self.transform_dens(k_lamtau_xy, D0, S)
            D_lamtau_xz = self.transform_dens(k_lamtau_xz, D0, S)
            D_lamtau_yz = self.transform_dens(k_lamtau_yz, D0, S)

            # x #
            Dx = self.transform_dens(kx_, D_sig_xx, S)
            Dx += self.transform_dens(k_sig_xx, Dc_x_, S)
            Dx += self.transform_dens(ky_, D_sig_xy, S)
            Dx += self.transform_dens(k_sig_xy, Dc_y_, S)

            Dx += self.transform_dens(kz_, D_sig_xz, S)
            Dx += self.transform_dens(k_sig_xz, Dc_z_, S)

            Dx += self.transform_dens(kx, D_lamtau_xx, S)
            Dx += self.transform_dens(k_lamtau_xx, Dc_x, S)

            Dx += self.transform_dens(ky, D_lamtau_xy, S)
            Dx += self.transform_dens(k_lamtau_xy, Dc_y, S)

            Dx += self.transform_dens(kz, D_lamtau_xz, S)
            Dx += self.transform_dens(k_lamtau_xz, Dc_z, S)

            # y #
            Dy = self.transform_dens(kx_, D_sig_xy, S)
            Dy += self.transform_dens(k_sig_xy, Dc_x_, S)

            Dy += self.transform_dens(ky_, D_sig_yy, S)
            Dy += self.transform_dens(k_sig_yy, Dc_y_, S)

            Dy += self.transform_dens(kz_, D_sig_yz, S)
            Dy += self.transform_dens(k_sig_yz, Dc_z_, S)

            Dy += self.transform_dens(kx, D_lamtau_xy, S)
            Dy += self.transform_dens(k_lamtau_xy, Dc_x, S)

            Dy += self.transform_dens(ky, D_lamtau_yy, S)
            Dy += self.transform_dens(k_lamtau_yy, Dc_y, S)

            Dy += self.transform_dens(kz, D_lamtau_yz, S)
            Dy += self.transform_dens(k_lamtau_yz, Dc_z, S)

            # z #
            Dz = self.transform_dens(kx_, D_sig_xz, S)
            Dz += self.transform_dens(k_sig_xz, Dc_x_, S)

            Dz += self.transform_dens(ky_, D_sig_yz, S)
            Dz += self.transform_dens(k_sig_yz, Dc_y_, S)

            Dz += self.transform_dens(kz_, D_sig_zz, S)
            Dz += self.transform_dens(k_sig_zz, Dc_z_, S)

            Dz += self.transform_dens(kx, D_lamtau_xz, S)
            Dz += self.transform_dens(k_lamtau_xz, Dc_x, S)

            Dz += self.transform_dens(ky, D_lamtau_yz, S)
            Dz += self.transform_dens(k_lamtau_yz, Dc_y, S)

            Dz += self.transform_dens(kz, D_lamtau_zz, S)
            Dz += self.transform_dens(k_lamtau_zz, Dc_z, S)

            density_list.append(Dx)
            density_list.append(Dy)
            density_list.append(Dz)

        return density_list

    def get_fock_dict_II(self, wi, density_list, mo, molecule, ao_basis):
        """
        Computes the compounded second-order Fock matrics used for the
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
            fock_file = re.sub(r'\.h5$', r'', self.checkpoint_file)
            fock_file += '_tpa_fock_2_full.h5'
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
        dist_focks = self.get_fock_r(mo, density_list, molecule, ao_basis,
                                     'real_and_imag')
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
        Contracts E[3]

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

            vec_pack = self.collect_vectors_in_columns(vec_pack)

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

            k_lamtau_xx = kXY[(('N_lamtau_xx', w), 0)].T
            k_lamtau_yy = kXY[(('N_lamtau_yy', w), 0)].T
            k_lamtau_zz = kXY[(('N_lamtau_zz', w), 0)].T

            k_lamtau_xy = kXY[(('N_lamtau_xy', w), 0)].T
            k_lamtau_xz = kXY[(('N_lamtau_xz', w), 0)].T
            k_lamtau_yz = kXY[(('N_lamtau_yz', w), 0)].T

            # Focks #

            # x

            zeta_sig_xx = self.xi(k_x_, k_sig_xx, f_x_, f_sig_xx, F0_a)
            zeta_sig_yy = self.xi(k_x_, k_sig_yy, f_x_, f_sig_yy, F0_a)
            zeta_sig_zz = self.xi(k_x_, k_sig_zz, f_x_, f_sig_zz, F0_a)

            zeta_sig_xy = self.xi(k_y_, k_sig_xy, f_y_, f_sig_xy, F0_a)
            zeta_sig_xz = self.xi(k_z_, k_sig_xz, f_z_, f_sig_xz, F0_a)

            zeta_lamtau_xx = self.xi(k_x, k_lamtau_xx, f_x, f_lamtau_xx, F0_a)
            zeta_lamtau_yy = self.xi(k_x, k_lamtau_yy, f_x, f_lamtau_yy, F0_a)
            zeta_lamtau_zz = self.xi(k_x, k_lamtau_zz, f_x, f_lamtau_zz, F0_a)

            zeta_lamtau_xy = self.xi(k_y, k_lamtau_xy, f_y, f_lamtau_xy, F0_a)
            zeta_lamtau_xz = self.xi(k_z, k_lamtau_xz, f_z, f_lamtau_xz, F0_a)

            X_terms = (zeta_sig_xx + zeta_sig_xy +
                       zeta_sig_xz).T + (zeta_lamtau_xx + zeta_lamtau_xy +
                                         zeta_lamtau_xz).T + (0.5 * F123_x).T
            Ff_x = -2 * LinearSolver.lrmat2vec(X_terms, nocc, norb)
            Ff_x = self.anti_sym(Ff_x)
            f_iso_x[w] = Ff_x

            # y

            zeta_sig_yx = self.xi(k_x_, k_sig_xy, f_x_, f_sig_xy, F0_a)
            zeta_sig_yy = self.xi(k_y_, k_sig_yy, f_y_, f_sig_yy, F0_a)
            zeta_sig_yz = self.xi(k_z_, k_sig_yz, f_z_, f_sig_yz, F0_a)

            zeta_lamtau_yx = self.xi(k_x, k_lamtau_xy, f_x, f_lamtau_xy, F0_a)
            zeta_lamtau_yy = self.xi(k_y, k_lamtau_yy, f_y, f_lamtau_yy, F0_a)
            zeta_lamtau_yz = self.xi(k_z, k_lamtau_yz, f_z, f_lamtau_yz, F0_a)

            Y_terms = (zeta_sig_yx + zeta_sig_yy +
                       zeta_sig_yz).T + (zeta_lamtau_yx + zeta_lamtau_yy +
                                         zeta_lamtau_yz).T + (0.5 * F123_y).T
            Ff_y = -2 * LinearSolver.lrmat2vec(Y_terms, nocc, norb)
            Ff_y = self.anti_sym(Ff_y)
            f_iso_y[w] = Ff_y

            # z

            zeta_sig_zx = self.xi(k_x_, k_sig_xz, f_x_, f_sig_xz, F0_a)
            zeta_sig_zy = self.xi(k_y_, k_sig_yz, f_y_, f_sig_yz, F0_a)
            zeta_sig_zz = self.xi(k_z_, k_sig_zz, f_z_, f_sig_zz, F0_a)

            zeta_lamtau_zx = self.xi(k_x, k_lamtau_xz, f_x, f_lamtau_xz, F0_a)
            zeta_lamtau_zy = self.xi(k_y, k_lamtau_yz, f_y, f_lamtau_yz, F0_a)
            zeta_lamtau_zz = self.xi(k_z, k_lamtau_zz, f_z, f_lamtau_zz, F0_a)

            Z_terms = (zeta_sig_zx + zeta_sig_zy +
                       zeta_sig_zz).T + (zeta_lamtau_zx + zeta_lamtau_zy +
                                         zeta_lamtau_zz).T + (0.5 * F123_z).T
            Ff_z = -2 * LinearSolver.lrmat2vec(Z_terms, nocc, norb)
            Ff_z = self.anti_sym(Ff_z)
            f_iso_z[w] = Ff_z

        return {'f_iso_x': f_iso_x, 'f_iso_y': f_iso_y, 'f_iso_z': f_iso_z}

    def get_other_terms(self, wi, track, Nx, Nxy, X, kX, kXY, da, nocc, norb):
        """
        Computes the terms involving X[3],A[3],X[2],A[2] in the isotropic cubic
        response function

        :param wi:
            A list containing all the frequencies
        :param track:
            A list that contains information about what γ components that are
            to be computed and which freqs
        :param Nx:
            A dictonary containing all the single-index response vectors
        :param Nxy:
            A dictonary containing all the two-index response vectors
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
                    'Na': Nx['Na'][(op_a, wa)],
                    'Nb': Nx['Nb'][(op_b, wb)],
                    'Nc_Nb': Nx['Nb'][(op_c, -wc)],
                    'Nd': Nx['Nd'][(op_d, wd)],
                    'kb': kX['Nb'][(op_b, wb)],
                    'kc_kb': kX['Nb'][(op_c, -wc)],
                    'kd': kX['Nd'][(op_d, wd)],
                    'A': X[op_a],
                    'B': X[op_b],
                    'C': X[op_c],
                    'D': X[op_d],
                })

        ave, res = divmod(len(inp_list), self.nodes)
        counts = [ave + 1 if p < res else ave for p in range(self.nodes)]
        starts = [sum(counts[:p]) for p in range(self.nodes)]
        ends = [sum(counts[:p + 1]) for p in range(self.nodes)]

        list_x3_a3 = [
            self.get_x3_a3(inp, da, nocc, norb)
            for inp in inp_list[starts[self.rank]:ends[self.rank]]
        ]

        list_x3_a3 = self.comm.gather(list_x3_a3, root=mpi_master())

        if self.rank == mpi_master():
            for terms in list_x3_a3:
                for term in terms:
                    key = term['key']
                    if key not in na_x3_ny_nz_dict:
                        na_x3_ny_nz_dict[key] = 0.0
                    if key not in na_a3_nx_ny_dict:
                        na_a3_nx_ny_dict[key] = 0.0
                    na_x3_ny_nz_dict[key] += term['x3']
                    na_a3_nx_ny_dict[key] += term['a3']

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
                Na = Nx['Na'][(op_a, wa)]
                A = X[op_a]

                for op_b in 'xyz':
                    op_ab = op_a + op_b if op_a <= op_b else op_b + op_a

                    # CD
                    kcd = kXY[(('N_lamtau_' + op_ab, w), wcd)]
                    Ncd = Nxy[(('N_lamtau_' + op_ab, w), wcd)]
                    Nb = Nx['Nb'][(op_b, w)]
                    kb = kX['Nb'][(op_b, w)]
                    B = X[op_b]

                    inp_list.append({
                        'flag': 'CD',
                        'freq': w,
                        'kcd': kcd,
                        'Ncd': Ncd,
                        'Na': Na,
                        'Nb': Nb,
                        'kb': kb,
                        'A': A,
                        'B': B,
                    })

                    # BD
                    op_c = op_b
                    op_ac = op_ab
                    kbd = kXY[(('N_sig_' + op_ac, w), wbd)]
                    Nbd = Nxy[(('N_sig_' + op_ac, w), wbd)]
                    Nc_Nb = Nx['Nb'][(op_c, -wc)]
                    kc_kb = kX['Nb'][(op_c, -wc)]
                    C = X[op_c]

                    inp_list.append({
                        'flag': 'BD',
                        'freq': w,
                        'kbd': kbd,
                        'Nbd': Nbd,
                        'Na': Na,
                        'Nc_Nb': Nc_Nb,
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
        Na = inp_dict['Na']
        Nb = inp_dict['Nb']
        Nc = self.flip_yz(inp_dict['Nc_Nb'])  # gets Nc from Nb
        Nd = inp_dict['Nd']
        kb = inp_dict['kb']
        kc = -inp_dict['kc_kb'].T.conj()  # gets kc from kb
        kd = inp_dict['kd']
        A = inp_dict['A']
        B = inp_dict['B']
        C = inp_dict['C']
        D = inp_dict['D']

        # Na X[3]NyNz

        na_x3_ny_nz -= np.dot(Na.T, self.x3_contract(kc, kd, B, da, nocc, norb))
        na_x3_ny_nz -= np.dot(Na.T, self.x3_contract(kd, kc, B, da, nocc, norb))
        na_x3_ny_nz -= np.dot(Na.T, self.x3_contract(kd, kb, C, da, nocc, norb))
        na_x3_ny_nz -= np.dot(Na.T, self.x3_contract(kb, kd, C, da, nocc, norb))
        na_x3_ny_nz -= np.dot(Na.T, self.x3_contract(kb, kc, D, da, nocc, norb))
        na_x3_ny_nz -= np.dot(Na.T, self.x3_contract(kc, kb, D, da, nocc, norb))

        # NaA[3]NxNy

        na_a3_nx_ny += np.dot(self.a3_contract(kb, kc, A, da, nocc, norb), Nd)
        na_a3_nx_ny += np.dot(self.a3_contract(kb, kd, A, da, nocc, norb), Nc)
        na_a3_nx_ny += np.dot(self.a3_contract(kc, kb, A, da, nocc, norb), Nd)
        na_a3_nx_ny += np.dot(self.a3_contract(kc, kd, A, da, nocc, norb), Nb)
        na_a3_nx_ny += np.dot(self.a3_contract(kd, kb, A, da, nocc, norb), Nc)
        na_a3_nx_ny += np.dot(self.a3_contract(kd, kc, A, da, nocc, norb), Nb)

        return {
            'key': (w, -w, w),
            'x3': (1. / 15) * na_x3_ny_nz,
            'a3': (1. / 15) * na_a3_nx_ny,
        }

    def get_t4(self, wi, e4_dict, Nx, kX, track, da, nocc, norb):
        """
        Computes the contraction of the E[4] tensor with that of the S[4] and
        R[4] tensors to return the contraction of T[4] as a dictonary of
        vectors. T[4]NxNyNz = (E^[4]-ω_1S^[4]-ω_1S^[4]-ω_3S^[4]-γiR^[4])

        :param wi:
            A list of all the freqs
        :param e4_dict:
            A dictonary of all the E[4] contraction
        :param Nx:
            A dictonary with all the single index response vectors
        :param kX:
            A dictonray containng all the response matricies
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

        T4term = {}
        S4 = self.S4_dict(wi, kX, track, da, nocc, norb)

        if self.damping > 0:
            R4term = self.get_r4(wi, kX, Nx, track, da, nocc, norb)

        comp_per_freq = len(track) // len(wi)

        for i in range(len(wi)):
            vals = track[i * comp_per_freq].split(',')
            w = float(vals[1])
            ww = float(vals[1])

            t4term = (np.matmul(Nx['Na'][('x', w)],
                                e4_dict['f_iso_x'][ww] - S4[('x', ww)]) +
                      np.matmul(Nx['Na'][('y', w)],
                                e4_dict['f_iso_y'][ww] - S4[('y', ww)]) +
                      np.matmul(Nx['Na'][('z', w)],
                                e4_dict['f_iso_z'][ww] - S4[('z', ww)]))

            if self.damping > 0:
                t4term += (R4term[('x', ww)] + R4term[('y', ww)] +
                           R4term[('z', ww)])

            T4term[(ww, -ww, ww)] = -(1. / 15) * t4term

        return T4term

    def S4_dict(self, wi, kX, track, D0, nocc, norb):
        """
        Computes the S4 contractions

        :param wi:
            A list of all the freqs
        :param kX:
            A dict with all the response matricies in MO basis
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

        S4terms = {}
        comp_per_freq = len(track) // len(wi)

        for j in range(len(wi)):
            vals = track[j * comp_per_freq].split(',')
            w = float(vals[1])
            w1 = float(vals[1])
            w2 = float(vals[2])
            w3 = float(vals[3])

            S4_term_x = 0
            S4_term_y = 0
            S4_term_z = 0

            for i in range(j * comp_per_freq, (j + 1) * comp_per_freq):
                comp_i = track[i]

                kB = kX['Nb'][(comp_i[1], w1)]
                kC = -kX['Nb'][(comp_i[2], -w2)].T.conj()
                kD = kX['Nd'][(comp_i[3], w3)]

                if comp_i[0] in 'x':
                    S4_term_x += w1 * self.s4(kB, kC, kD, D0, nocc, norb)
                    S4_term_x += w2 * self.s4(kC, kB, kD, D0, nocc, norb)
                    S4_term_x += w3 * self.s4(kD, kB, kC, D0, nocc, norb)

                elif comp_i[0] in 'y':
                    S4_term_y += w1 * self.s4(kB, kC, kD, D0, nocc, norb)
                    S4_term_y += w2 * self.s4(kC, kB, kD, D0, nocc, norb)
                    S4_term_y += w3 * self.s4(kD, kB, kC, D0, nocc, norb)

                elif comp_i[0] == 'z':
                    S4_term_z += w1 * self.s4(kB, kC, kD, D0, nocc, norb)
                    S4_term_z += w2 * self.s4(kC, kB, kD, D0, nocc, norb)
                    S4_term_z += w3 * self.s4(kD, kB, kC, D0, nocc, norb)

            S4terms[('x', w)] = -S4_term_x
            S4terms[('y', w)] = -S4_term_y
            S4terms[('z', w)] = -S4_term_z

        return S4terms

    def s4(self, k1, k2, k3, D, nocc, norb):
        """
        Returns the contraction of S[4] for S[4] dict

        :param k1:
            A response matrix
        :param k2:
            A response matrix
        :param k3:
            A response matrix
        :param D:
            A density matrix
        :param nocc:
            The number of occupied orbtials
        :param norb:
            The number of total orbitals

        :return:
            The contraction of S[4] for S[4] dict
        """

        S4_123 = self.S4contract(k1, k2, k3, D, nocc, norb)
        S4_132 = self.S4contract(k1, k3, k2, D, nocc, norb)
        A = S4_123 + S4_132

        return A

    def S4contract(self, k1, k2, k3, D, nocc, norb):
        """
        Returns the contraction of S[4] for S[4] dict

        :param k1:
            A response matrix
        :param k2:
            A response matrix
        :param k3:
            A response matrix
        :param D:
            A density matrix
        :param nocc:
            The number of occupied orbtials
        :param norb:
            The number of total orbitals

        :return:
            The contraction of S[4] for S[4] dict
        """

        S4N1N2N3 = self.commut(self.commut(k3, self.commut(k2, k1)), D.T)
        S4N1N2N3 = [
            LinearSolver.lrmat2vec(S4N1N2N3.real, nocc, norb),
            LinearSolver.lrmat2vec(S4N1N2N3.imag, nocc, norb)
        ]
        S4N1N2N3_c = S4N1N2N3[0] + 1j * S4N1N2N3[1]
        return (2. / 6) * S4N1N2N3_c

    def get_r4(self, freqs, kX, Nx, track, d_a_mo, nocc, norb):
        """
        Returns a dict with all the R[4]NxNyNz contractions for the subsequent
        T[4] contraction

        :param freqs:
            A list of all the frequencies
        :param kX:
            A dictonary of all the first-order response matrices
        :param Nx:
            A dictonary of all the first-order response vectors
        :param track:
            A list of all the cubic response function components that are to be
            computed for the isotropic
        :param d_a_mo:
            The zeroth-order density in MO basis
        :param nocc:
            The number of occupied orbitals
        :param norb:
            The total number of orbitals

        :return:
            A dict with all the R[4]NxNyNz contractions for the subsequent
            T[4] contraction
        """

        R4terms = {}

        damp = self.damping
        comp_per_freq = len(track) // len(freqs)

        for j in range(len(freqs)):
            vals = track[j * comp_per_freq].split(',')
            w1 = float(vals[1])
            w2 = float(vals[2])
            w3 = float(vals[3])
            w_s = w1 + w2 + w3

            R4x = 0
            R4y = 0
            R4z = 0

            for i in range(j * comp_per_freq, (j + 1) * comp_per_freq):
                comp_i = track[i]

                # Na = Nx['Na'][(comp_i[0], w_s)]
                Nb = Nx['Nb'][(comp_i[1], w1)]
                Nc = self.flip_yz(Nx['Nb'][(comp_i[2], -w2)])  # gets Nc from Nb
                Nd = Nx['Nd'][(comp_i[3], w3)]
                kA = kX['Na'][(comp_i[0], w_s)]
                kB = kX['Nb'][(comp_i[1], w1)]
                kC = -kX['Nb'][(comp_i[2], -w2)].T.conj()  # gets kc from kb
                kD = kX['Nd'][(comp_i[3], w3)]

                Nb_h = self.flip_xy(Nb)
                Nc_h = self.flip_xy(Nc)
                Nd_h = self.flip_xy(Nd)

                if comp_i[0] == 'x':
                    R4x += -1j * damp * np.matmul(
                        Nd_h, self.s4_for_r4(kA.T, kB, kC, d_a_mo, nocc, norb))
                    R4x += -1j * damp * np.matmul(
                        Nc_h, self.s4_for_r4(kA.T, kB, kD, d_a_mo, nocc, norb))
                    R4x += -1j * damp * np.matmul(
                        Nd_h, self.s4_for_r4(kA.T, kC, kB, d_a_mo, nocc, norb))
                    R4x += -1j * damp * np.matmul(
                        Nb_h, self.s4_for_r4(kA.T, kC, kD, d_a_mo, nocc, norb))
                    R4x += -1j * damp * np.matmul(
                        Nc_h, self.s4_for_r4(kA.T, kD, kB, d_a_mo, nocc, norb))
                    R4x += -1j * damp * np.matmul(
                        Nb_h, self.s4_for_r4(kA.T, kD, kC, d_a_mo, nocc, norb))

                elif comp_i[0] == 'y':
                    R4y += -1j * damp * np.matmul(
                        Nd_h, self.s4_for_r4(kA.T, kB, kC, d_a_mo, nocc, norb))
                    R4y += -1j * damp * np.matmul(
                        Nc_h, self.s4_for_r4(kA.T, kB, kD, d_a_mo, nocc, norb))
                    R4y += -1j * damp * np.matmul(
                        Nd_h, self.s4_for_r4(kA.T, kC, kB, d_a_mo, nocc, norb))
                    R4y += -1j * damp * np.matmul(
                        Nb_h, self.s4_for_r4(kA.T, kC, kD, d_a_mo, nocc, norb))
                    R4y += -1j * damp * np.matmul(
                        Nc_h, self.s4_for_r4(kA.T, kD, kB, d_a_mo, nocc, norb))
                    R4y += -1j * damp * np.matmul(
                        Nb_h, self.s4_for_r4(kA.T, kD, kC, d_a_mo, nocc, norb))

                elif comp_i[0] == 'z':
                    R4z += -1j * damp * np.matmul(
                        Nd_h, self.s4_for_r4(kA.T, kB, kC, d_a_mo, nocc, norb))
                    R4z += -1j * damp * np.matmul(
                        Nc_h, self.s4_for_r4(kA.T, kB, kD, d_a_mo, nocc, norb))
                    R4z += -1j * damp * np.matmul(
                        Nd_h, self.s4_for_r4(kA.T, kC, kB, d_a_mo, nocc, norb))
                    R4z += -1j * damp * np.matmul(
                        Nb_h, self.s4_for_r4(kA.T, kC, kD, d_a_mo, nocc, norb))
                    R4z += -1j * damp * np.matmul(
                        Nc_h, self.s4_for_r4(kA.T, kD, kB, d_a_mo, nocc, norb))
                    R4z += -1j * damp * np.matmul(
                        Nb_h, self.s4_for_r4(kA.T, kD, kC, d_a_mo, nocc, norb))

            R4terms[('x', w1)] = -R4x
            R4terms[('y', w1)] = -R4y
            R4terms[('z', w1)] = -R4z

        return R4terms

    def s4_for_r4(self, k1, k2, k3, D, nocc, norb):
        """
        Returns the contraction of S[4] for the contraction of R[4]

        :param k1:
            A response matrix
        :param k2:
            A response matrix
        :param k3:
            A response matrix
        :param D:
            A density matrix
        :param nocc:
            The number of occupied orbtials
        :param norb:
            The number of total orbitals

        :return:
            The contraction of S[4] for the contraction of R[4]
        """

        S4_123 = self.S4contract(k1, k2, k3, D, nocc, norb)
        return S4_123

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

        width = 94

        w_str = 'Gamma tensor components computed per frequency'
        self.ostream.print_blank()
        self.ostream.print_header(w_str.ljust(width))
        self.ostream.print_blank()

        for a in range(len(comp) // len(freqs)):
            w_str = str(a + 1) + '. ' + str(comp[a].split(',')[0])
            self.ostream.print_header(w_str.ljust(width))

        self.ostream.print_blank()
        self.ostream.print_blank()

        w_str = 'Gamma Tensor Components at Given Frequencies'
        self.ostream.print_header(w_str.ljust(width))
        self.ostream.print_blank()

        for w in freqs:
            title = '{:<9s} {:>12s} {:>20s} {:>21s}'.format(
                'Component', 'Frequency', 'Real', 'Imaginary')
            self.ostream.print_header(title.ljust(width))
            self.ostream.print_header(('-' * len(title)).ljust(width))

            self.print_component('T3', w, t3_dict[w, -w, w], width)
            self.print_component('T4', w, t4_dict[w, -w, w], width)
            self.print_component('X2', w, NaX2Nyz[w, -w, w], width)
            self.print_component('X3', w, NaX3NyNz[w, -w, w], width)
            self.print_component('A2', w, NxA2Nyz[w, -w, w], width)
            self.print_component('A3', w, NaA3NxNy[w, -w, w], width)
            self.print_component('Gamma', w, gamma[w, -w, w], width)

            self.ostream.print_blank()

        self.ostream.print_blank()
