import numpy as np
import time
import re

from .veloxchemlib import mpi_master
from .cppsolver import ComplexResponse
from .linearsolver import LinearSolver
from .tpadriver import TpaDriver


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

            kx_ = self.mo2ao(mo, kX['Nc'][('x', -w)]).T
            ky_ = self.mo2ao(mo, kX['Nc'][('y', -w)]).T
            kz_ = self.mo2ao(mo, kX['Nc'][('z', -w)]).T
            # create the first order single indexed densiteies #
            Dx = self.transform_dens(kx, D0, S)
            Dy = self.transform_dens(ky, D0, S)
            Dz = self.transform_dens(kz, D0, S)

            Dx_ = self.transform_dens(kx_, D0, S)
            Dy_ = self.transform_dens(ky_, D0, S)
            Dz_ = self.transform_dens(kz_, D0, S)
            # create the first order two indexed densities #
            # σ terms #
            D_sig_xx = 2 * (3 * self.transform_dens(kx, Dx, S) +
                            self.transform_dens(ky, Dy, S) +
                            self.transform_dens(kz, Dz, S))
            D_sig_yy = 2 * (self.transform_dens(kx, Dx, S) +
                            3 * self.transform_dens(ky, Dy, S) +
                            self.transform_dens(kz, Dz, S))
            D_sig_zz = 2 * (self.transform_dens(kx, Dx, S) +
                            self.transform_dens(ky, Dy, S) +
                            3 * self.transform_dens(kz, Dz, S))

            D_sig_xy = 2 * (self.transform_dens(ky, Dx, S) +
                            self.transform_dens(kx, Dy, S))
            D_sig_xz = 2 * (self.transform_dens(kx, Dz, S) +
                            self.transform_dens(kz, Dx, S))
            D_sig_yz = 2 * (self.transform_dens(ky, Dz, S) +
                            self.transform_dens(kz, Dy, S))
            # λ+τ terms #
            D_lamtau_xx = 2 * (3 * self.transform_dens(kx_, Dx, S) +
                               self.transform_dens(ky_, Dy, S) +
                               self.transform_dens(kz_, Dz, S))
            D_lamtau_xx += 2 * (3 * self.transform_dens(kx, Dx_, S) +
                                self.transform_dens(ky, Dy_, S) +
                                self.transform_dens(kz, Dz_, S))

            D_lamtau_yy = 2 * (self.transform_dens(kx_, Dx, S) +
                               3 * self.transform_dens(ky_, Dy, S) +
                               self.transform_dens(kz_, Dz, S))
            D_lamtau_yy += 2 * (self.transform_dens(kx, Dx_, S) +
                                3 * self.transform_dens(ky, Dy_, S) +
                                self.transform_dens(kz, Dz_, S))

            D_lamtau_zz = 2 * (self.transform_dens(kx_, Dx, S) +
                               self.transform_dens(ky_, Dy, S) +
                               3 * self.transform_dens(kz_, Dz, S))
            D_lamtau_zz += 2 * (self.transform_dens(kx, Dx_, S) +
                                self.transform_dens(ky, Dy_, S) +
                                3 * self.transform_dens(kz, Dz_, S))

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

        time_start_fock = time.time()
        fock_list = self.get_fock_r(mo, density_list, molecule, ao_basis,
                                    'real_and_imag')
        time_end_fock = time.time()

        total_time_fock = time_end_fock - time_start_fock
        self.print_fock_time(total_time_fock)

        Fock = {}

        if self.rank == mpi_master():
            f_sig_xx = {}
            f_sig_yy = {}
            f_sig_zz = {}
            f_sig_xy = {}
            f_sig_xz = {}
            f_sig_yz = {}

            f_lamtau_xx = {}
            f_lamtau_yy = {}
            f_lamtau_zz = {}
            f_lamtau_xy = {}
            f_lamtau_xz = {}
            f_lamtau_yz = {}

            f_lam_sig_tau_x = {}
            f_lam_sig_tau_y = {}
            f_lam_sig_tau_z = {}

            for count, w in enumerate(wi):
                f_sig_xx[w] = fock_list[15 * count]
                f_sig_yy[w] = fock_list[15 * count + 1]
                f_sig_zz[w] = fock_list[15 * count + 2]
                f_sig_xy[w] = fock_list[15 * count + 3]
                f_sig_xz[w] = fock_list[15 * count + 4]
                f_sig_yz[w] = fock_list[15 * count + 5]

                f_lamtau_xx[w] = fock_list[15 * count + 6]
                f_lamtau_yy[w] = fock_list[15 * count + 7]
                f_lamtau_zz[w] = fock_list[15 * count + 8]
                f_lamtau_xy[w] = fock_list[15 * count + 9]
                f_lamtau_xz[w] = fock_list[15 * count + 10]
                f_lamtau_yz[w] = fock_list[15 * count + 11]

                f_lam_sig_tau_x[w] = fock_list[15 * count + 12]
                f_lam_sig_tau_y[w] = fock_list[15 * count + 13]
                f_lam_sig_tau_z[w] = fock_list[15 * count + 14]

            Fock = {
                'F0': F0_a,
                'f_sig_xx': f_sig_xx,
                'f_sig_yy': f_sig_yy,
                'f_sig_zz': f_sig_zz,
                'f_sig_xy': f_sig_xy,
                'f_sig_xz': f_sig_xz,
                'f_sig_yz': f_sig_yz,
                'f_lamtau_xx': f_lamtau_xx,
                'f_lamtau_yy': f_lamtau_yy,
                'f_lamtau_zz': f_lamtau_zz,
                'f_lamtau_xy': f_lamtau_xy,
                'f_lamtau_xz': f_lamtau_xz,
                'f_lamtau_yz': f_lamtau_yz,
                'F123_x': f_lam_sig_tau_x,
                'F123_y': f_lam_sig_tau_y,
                'F123_z': f_lam_sig_tau_z
            }

        return Fock

    def get_e4(self, wi, kX, fo, nocc, norb):
        """
        Contracts E[4]n_xNyNz for the isotropic cubic response function. Takes
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

        F0 = fo['F0']

        for w in wi:
            # Get all the response matrices and Fock matrices
            kx = kX['Nb'][('x', w)].T
            ky = kX['Nb'][('y', w)].T
            kz = kX['Nb'][('z', w)].T

            kx_ = kX['Nc'][('x', -w)].T
            ky_ = kX['Nc'][('y', -w)].T
            kz_ = kX['Nc'][('z', -w)].T

            Fx = fo['Fb'][('x', w)]
            Fy = fo['Fb'][('y', w)]
            Fz = fo['Fb'][('z', w)]

            Fx_ = fo['Fc'][('x', -w)]
            Fy_ = fo['Fc'][('y', -w)]
            Fz_ = fo['Fc'][('z', -w)]

            f_lamtau_xx = 3 * fo['f_lamtau_xx'][w]
            f_lamtau_yy = 3 * fo['f_lamtau_yy'][w]
            f_lamtau_zz = 3 * fo['f_lamtau_zz'][w]

            f_lamtau_xy = 3 * fo['f_lamtau_xy'][w]
            f_lamtau_xz = 3 * fo['f_lamtau_xz'][w]
            f_lamtau_yz = 3 * fo['f_lamtau_yz'][w]

            f_sig_xx = 3 * fo['f_sig_xx'][w]
            f_sig_yy = 3 * fo['f_sig_yy'][w]
            f_sig_zz = 3 * fo['f_sig_zz'][w]

            f_sig_xy = 3 * fo['f_sig_xy'][w]
            f_sig_xz = 3 * fo['f_sig_xz'][w]
            f_sig_yz = 3 * fo['f_sig_yz'][w]

            # computes all the compounded Φ_αβ, see article, where small phi
            # here is defined as:
            #   φ(κa,κb,Fb,F0) = [κa,[κb,F0]+3Fb]
            #   Φ_αα^σ = φ(κα,κα,Fα,F0) + φ(κα,κα,Fα,F0) +
            #            Σ_{ρ}^{x,y,z}[φ(κρ,κρ,Fρ,F0)] (for α=β)
            #   Φ_αβ^σ = φ(κa,κa,Fb,F0) + φ(κb,κb,Fb,F0) (for α≠β)
            #  For the Φ_{αβ}^{λ+τ} component see article.

            Phi_sig_xx = 2 * (3 * self.phi(kx, kx, Fx, F0) + self.phi(
                ky, ky, Fy, F0) + self.phi(kz, kz, Fz, F0))

            Phi_sig_yy = 2 * (self.phi(kx, kx, Fx, F0) +
                              3 * self.phi(ky, ky, Fy, F0) +
                              self.phi(kz, kz, Fz, F0))

            Phi_sig_zz = 2 * (self.phi(kx, kx, Fx, F0) + self.phi(
                ky, ky, Fy, F0) + 3 * self.phi(kz, kz, Fz, F0))

            Phi_sig_xy = (2 * self.phi(kx, ky, Fy, F0) +
                          2 * self.phi(ky, kx, Fx, F0))

            Phi_sig_xz = (2 * self.phi(kx, kz, Fz, F0) +
                          2 * self.phi(kz, kx, Fx, F0))

            Phi_sig_yz = (2 * self.phi(ky, kz, Fz, F0) +
                          2 * self.phi(kz, ky, Fy, F0))

            Phi_lamtau_xx = 2 * (
                3 * self.phi(kx, kx_, Fx_, F0) + 3 * self.phi(kx_, kx, Fx, F0) +
                self.phi(ky, ky_, Fy_, F0) + self.phi(ky_, ky, Fy, F0) +
                self.phi(kz, kz_, Fz_, F0) + self.phi(kz_, kz, Fz, F0))

            Phi_lamtau_yy = 2 * (
                self.phi(kx, kx_, Fx_, F0) + self.phi(kx_, kx, Fx, F0) +
                3 * self.phi(ky, ky_, Fy_, F0) + 3 * self.phi(ky_, ky, Fy, F0) +
                self.phi(kz, kz_, Fz_, F0) + self.phi(kz_, kz, Fz, F0))

            Phi_lamtau_zz = 2 * (
                self.phi(kx, kx_, Fx_, F0) + self.phi(kx_, kx, Fx, F0) +
                self.phi(ky, ky_, Fy_, F0) + self.phi(ky_, ky, Fy, F0) +
                3 * self.phi(kz, kz_, Fz_, F0) + 3 * self.phi(kz_, kz, Fz, F0))

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

            f_x = fo['F123_x'][w]

            f_x += (self.commut(kx, Phi_lamtau_xx + f_lamtau_xx) +
                    self.commut(ky, Phi_lamtau_xy + f_lamtau_xy) +
                    self.commut(kz, Phi_lamtau_xz + f_lamtau_xz))

            f_x += (self.commut(kx_, Phi_sig_xx + f_sig_xx) +
                    self.commut(ky_, Phi_sig_xy + f_sig_xy) +
                    self.commut(kz_, Phi_sig_xz + f_sig_xz))

            # Taking the non redundant matrix elements {i,s} and forming the
            # anti-symmetric Fock vector
            f_x = -2 * 1. / 6 * LinearSolver.lrmat2vec(f_x.T, nocc, norb)
            f_x = self.anti_sym(f_x)
            f_iso_x[w] = f_x

            # y

            # Creating the transformed total Fock matrices
            f_y = fo['F123_y'][w]

            f_y += (self.commut(kx, Phi_lamtau_xy + f_lamtau_xy) +
                    self.commut(ky, Phi_lamtau_yy + f_lamtau_yy) +
                    self.commut(kz, Phi_lamtau_yz + f_lamtau_yz))

            f_y += (self.commut(kx_, Phi_sig_xy + f_sig_xy) +
                    self.commut(ky_, Phi_sig_yy + f_sig_yy) +
                    self.commut(kz_, Phi_sig_yz + f_sig_yz))

            # Taking the non redundant matrix elements {i,s} and forming the
            # anti-symmetric Fock vector
            f_y = -2 * 1. / 6 * LinearSolver.lrmat2vec(f_y.T, nocc, norb)
            f_y = self.anti_sym(f_y)
            f_iso_y[w] = f_y

            # z

            # Creating the transformed total Fock matrices
            f_z = fo['F123_z'][w]

            f_z += (self.commut(kx, Phi_lamtau_xz + f_lamtau_xz) +
                    self.commut(ky, Phi_lamtau_yz + f_lamtau_yz) +
                    self.commut(kz, Phi_lamtau_zz + f_lamtau_zz))

            f_z += (self.commut(kx_, Phi_sig_xz + f_sig_xz) +
                    self.commut(ky_, Phi_sig_yz + f_sig_yz) +
                    self.commut(kz_, Phi_sig_zz + f_sig_zz))

            # Taking the non redundant matrix elements {i,s} and forming the
            # anti-symmetric Fock vector
            f_z = -2 * 1. / 6 * LinearSolver.lrmat2vec(f_z.T, nocc, norb)
            f_z = self.anti_sym(f_z)
            f_iso_z[w] = f_z

        return {'f_iso_x': f_iso_x, 'f_iso_y': f_iso_y, 'f_iso_z': f_iso_z}

    def get_n_xy(self, w, d_a_mo, X, fock_dict, kX, nocc, norb, molecule,
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

        if self.rank == mpi_master():

            # Get second-order gradients
            xy_dict = self.get_xy(d_a_mo, X, w, fock_dict, kX, nocc, norb)

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
        if self.checkpoint_file is not None:
            N_total_drv.checkpoint_file = re.sub(r'\.h5$', r'',
                                                 self.checkpoint_file)
            N_total_drv.checkpoint_file += '_tpa_2_full.h5'

        # commutpute second-order response vectors
        N_total_results = N_total_drv.compute(molecule, ao_basis, scf_tensors,
                                              xy_dict)

        if self.rank == mpi_master():
            n_xy_dict = N_total_results['solutions']
            kxy_dict = N_total_results['kappas']
            FXY_2_dict = N_total_results['focks']

            return (n_xy_dict, kxy_dict, FXY_2_dict, xy_dict)
        else:
            return (None, None, None, None)

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

        if self.rank == mpi_master():
            for w in wi:
                mu_x = X['x']
                mu_y = X['y']
                mu_z = X['z']

                kx = kX['Nb'][('x', w)].T
                ky = kX['Nb'][('y', w)].T
                kz = kX['Nb'][('z', w)].T
                kx_ = kX['Nc'][('x', -w)].T
                ky_ = kX['Nc'][('y', -w)].T
                kz_ = kX['Nc'][('z', -w)].T

                F0 = Fock['F0']
                f_x = Fock['Fb'][('x', w)]
                f_y = Fock['Fb'][('y', w)]
                f_z = Fock['Fb'][('z', w)]

                f_x_ = Fock['Fc'][('x', -w)]
                f_y_ = Fock['Fc'][('y', -w)]
                f_z_ = Fock['Fc'][('z', -w)]

                # BD σ gradients #

                key = (('N_sig_xx', w), 2 * w)
                mat = (3 * self.xi(kx, kx, f_x, f_x, F0) +
                       self.xi(ky, ky, f_y, f_y, F0) +
                       self.xi(kz, kz, f_z, f_z, F0) +
                       0.5 * Fock['f_sig_xx'][w]).T
                xy_dict[key] = self.anti_sym(
                    -2 * LinearSolver.lrmat2vec(mat, nocc, norb))
                xy_dict[key] -= 6 * self.x2_contract(kx.T, mu_x, d_a_mo, nocc,
                                                     norb)
                xy_dict[key] -= 2 * self.x2_contract(ky.T, mu_y, d_a_mo, nocc,
                                                     norb)
                xy_dict[key] -= 2 * self.x2_contract(kz.T, mu_z, d_a_mo, nocc,
                                                     norb)

                key = (('N_sig_yy', w), 2 * w)
                mat = (self.xi(kx, kx, f_x, f_x, F0) +
                       3 * self.xi(ky, ky, f_y, f_y, F0) +
                       self.xi(kz, kz, f_z, f_z, F0) +
                       0.5 * Fock['f_sig_yy'][w]).T
                xy_dict[key] = self.anti_sym(
                    -2 * LinearSolver.lrmat2vec(mat, nocc, norb))
                xy_dict[key] -= 2 * self.x2_contract(kx.T, mu_x, d_a_mo, nocc,
                                                     norb)
                xy_dict[key] -= 6 * self.x2_contract(ky.T, mu_y, d_a_mo, nocc,
                                                     norb)
                xy_dict[key] -= 2 * self.x2_contract(kz.T, mu_z, d_a_mo, nocc,
                                                     norb)

                key = (('N_sig_zz', w), 2 * w)
                mat = (self.xi(kx, kx, f_x, f_x, F0) +
                       self.xi(ky, ky, f_y, f_y, F0) +
                       3 * self.xi(kz, kz, f_z, f_z, F0) +
                       0.5 * Fock['f_sig_zz'][w]).T
                xy_dict[key] = self.anti_sym(
                    -2 * LinearSolver.lrmat2vec(mat, nocc, norb))
                xy_dict[key] -= 2 * self.x2_contract(kx.T, mu_x, d_a_mo, nocc,
                                                     norb)
                xy_dict[key] -= 2 * self.x2_contract(ky.T, mu_y, d_a_mo, nocc,
                                                     norb)
                xy_dict[key] -= 6 * self.x2_contract(kz.T, mu_z, d_a_mo, nocc,
                                                     norb)

                key = (('N_sig_xy', w), 2 * w)
                mat = (self.xi(ky, kx, f_y, f_x, F0) +
                       self.xi(kx, ky, f_x, f_y, F0) +
                       0.5 * Fock['f_sig_xy'][w]).T
                xy_dict[key] = self.anti_sym(
                    -2 * LinearSolver.lrmat2vec(mat, nocc, norb))
                xy_dict[key] -= 2 * self.x2_contract(ky.T, mu_x, d_a_mo, nocc,
                                                     norb)
                xy_dict[key] -= 2 * self.x2_contract(kx.T, mu_y, d_a_mo, nocc,
                                                     norb)

                key = (('N_sig_xz', w), 2 * w)
                mat = (self.xi(kz, kx, f_z, f_x, F0) +
                       self.xi(kx, kz, f_x, f_z, F0) +
                       0.5 * Fock['f_sig_xz'][w]).T
                xy_dict[key] = self.anti_sym(
                    -2 * LinearSolver.lrmat2vec(mat, nocc, norb))
                xy_dict[key] -= 2 * self.x2_contract(kz.T, mu_x, d_a_mo, nocc,
                                                     norb)
                xy_dict[key] -= 2 * self.x2_contract(kx.T, mu_z, d_a_mo, nocc,
                                                     norb)

                key = (('N_sig_yz', w), 2 * w)
                mat = (self.xi(kz, ky, f_z, f_y, F0) +
                       self.xi(ky, kz, f_y, f_z, F0) +
                       0.5 * Fock['f_sig_yz'][w]).T
                xy_dict[key] = self.anti_sym(
                    -2 * LinearSolver.lrmat2vec(mat, nocc, norb))
                xy_dict[key] -= 2 * self.x2_contract(kz.T, mu_y, d_a_mo, nocc,
                                                     norb)
                xy_dict[key] -= 2 * self.x2_contract(ky.T, mu_z, d_a_mo, nocc,
                                                     norb)

                # BC CD λ+τ gradients #

                key = (('N_lamtau_xx', w), 0)
                mat = (3 * 2 * self.xi(kx_, kx, f_x_, f_x, F0) +
                       2 * self.xi(ky_, ky, f_y_, f_y, F0) +
                       2 * self.xi(kz_, kz, f_z_, f_z, F0) +
                       0.5 * Fock['f_lamtau_xx'][w]).T
                xy_dict[key] = self.anti_sym(
                    -2 * LinearSolver.lrmat2vec(mat, nocc, norb))
                xy_dict[key] -= 6 * self.x2_contract(kx_.T, mu_x, d_a_mo, nocc,
                                                     norb)
                xy_dict[key] -= 2 * self.x2_contract(ky_.T, mu_y, d_a_mo, nocc,
                                                     norb)
                xy_dict[key] -= 2 * self.x2_contract(kz_.T, mu_z, d_a_mo, nocc,
                                                     norb)
                xy_dict[key] -= 6 * self.x2_contract(kx.T, mu_x, d_a_mo, nocc,
                                                     norb)
                xy_dict[key] -= 2 * self.x2_contract(ky.T, mu_y, d_a_mo, nocc,
                                                     norb)
                xy_dict[key] -= 2 * self.x2_contract(kz.T, mu_z, d_a_mo, nocc,
                                                     norb)

                key = (('N_lamtau_yy', w), 0)
                mat = (2 * self.xi(kx_, kx, f_x_, f_x, F0) +
                       3 * 2 * self.xi(ky_, ky, f_y_, f_y, F0) +
                       2 * self.xi(kz_, kz, f_z_, f_z, F0) +
                       0.5 * Fock['f_lamtau_yy'][w]).T
                xy_dict[key] = self.anti_sym(
                    -2 * LinearSolver.lrmat2vec(mat, nocc, norb))
                xy_dict[key] -= 2 * self.x2_contract(kx_.T, mu_x, d_a_mo, nocc,
                                                     norb)
                xy_dict[key] -= 6 * self.x2_contract(ky_.T, mu_y, d_a_mo, nocc,
                                                     norb)
                xy_dict[key] -= 2 * self.x2_contract(kz_.T, mu_z, d_a_mo, nocc,
                                                     norb)
                xy_dict[key] -= 2 * self.x2_contract(kx.T, mu_x, d_a_mo, nocc,
                                                     norb)
                xy_dict[key] -= 6 * self.x2_contract(ky.T, mu_y, d_a_mo, nocc,
                                                     norb)
                xy_dict[key] -= 2 * self.x2_contract(kz.T, mu_z, d_a_mo, nocc,
                                                     norb)

                key = (('N_lamtau_zz', w), 0)
                mat = (2 * self.xi(kx_, kx, f_x_, f_x, F0) +
                       2 * self.xi(ky_, ky, f_y_, f_y, F0) +
                       3 * 2 * self.xi(kz_, kz, f_z_, f_z, F0) +
                       0.5 * Fock['f_lamtau_zz'][w]).T
                xy_dict[key] = self.anti_sym(
                    -2 * LinearSolver.lrmat2vec(mat, nocc, norb))
                xy_dict[key] -= 2 * self.x2_contract(kx_.T, mu_x, d_a_mo, nocc,
                                                     norb)
                xy_dict[key] -= 2 * self.x2_contract(ky_.T, mu_y, d_a_mo, nocc,
                                                     norb)
                xy_dict[key] -= 6 * self.x2_contract(kz_.T, mu_z, d_a_mo, nocc,
                                                     norb)
                xy_dict[key] -= 2 * self.x2_contract(kx.T, mu_x, d_a_mo, nocc,
                                                     norb)
                xy_dict[key] -= 2 * self.x2_contract(ky.T, mu_y, d_a_mo, nocc,
                                                     norb)
                xy_dict[key] -= 6 * self.x2_contract(kz.T, mu_z, d_a_mo, nocc,
                                                     norb)

                key = (('N_lamtau_xy', w), 0)
                mat = (2 * self.xi(ky_, kx, f_y_, f_x, F0) +
                       2 * self.xi(kx_, ky, f_x_, f_y, F0) +
                       0.5 * Fock['f_lamtau_xy'][w]).T
                xy_dict[key] = self.anti_sym(
                    -2 * LinearSolver.lrmat2vec(mat, nocc, norb))
                xy_dict[key] -= 2 * self.x2_contract(ky.T, mu_x, d_a_mo, nocc,
                                                     norb)
                xy_dict[key] -= 2 * self.x2_contract(kx.T, mu_y, d_a_mo, nocc,
                                                     norb)
                xy_dict[key] -= 2 * self.x2_contract(ky_.T, mu_x, d_a_mo, nocc,
                                                     norb)
                xy_dict[key] -= 2 * self.x2_contract(kx_.T, mu_y, d_a_mo, nocc,
                                                     norb)

                key = (('N_lamtau_xz', w), 0)
                mat = (2 * self.xi(kz_, kx, f_z_, f_x, F0) +
                       2 * self.xi(kx_, kz, f_x_, f_z, F0) +
                       0.5 * Fock['f_lamtau_xz'][w]).T
                xy_dict[key] = self.anti_sym(
                    -2 * LinearSolver.lrmat2vec(mat, nocc, norb))
                xy_dict[key] -= 2 * self.x2_contract(kz.T, mu_x, d_a_mo, nocc,
                                                     norb)
                xy_dict[key] -= 2 * self.x2_contract(kx.T, mu_z, d_a_mo, nocc,
                                                     norb)
                xy_dict[key] -= 2 * self.x2_contract(kz_.T, mu_x, d_a_mo, nocc,
                                                     norb)
                xy_dict[key] -= 2 * self.x2_contract(kx_.T, mu_z, d_a_mo, nocc,
                                                     norb)

                key = (('N_lamtau_yz', w), 0)
                mat = (2 * self.xi(kz_, ky, f_z_, f_y, F0) +
                       2 * self.xi(ky_, kz, f_y_, f_z, F0) +
                       0.5 * Fock['f_lamtau_yz'][w]).T
                xy_dict[key] = self.anti_sym(
                    -2 * LinearSolver.lrmat2vec(mat, nocc, norb))
                xy_dict[key] -= 2 * self.x2_contract(kz.T, mu_y, d_a_mo, nocc,
                                                     norb)
                xy_dict[key] -= 2 * self.x2_contract(ky.T, mu_z, d_a_mo, nocc,
                                                     norb)
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

            kx_ = self.mo2ao(mo, kX['Nc'][('x', -w)]).T
            ky_ = self.mo2ao(mo, kX['Nc'][('y', -w)]).T
            kz_ = self.mo2ao(mo, kX['Nc'][('z', -w)]).T

            kx = self.mo2ao(mo, kX['Nb'][('x', w)]).T
            ky = self.mo2ao(mo, kX['Nb'][('y', w)]).T
            kz = self.mo2ao(mo, kX['Nb'][('z', w)]).T

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

        time_start_fock = time.time()
        fock_list = self.get_fock_r(mo, density_list, molecule, ao_basis,
                                    "real_and_imag")
        time_end_fock = time.time()

        total_time_fock = time_end_fock - time_start_fock
        self.print_fock_time(total_time_fock)

        fock_dict = {}

        if self.rank == mpi_master():
            F123_x = {}
            F123_y = {}
            F123_z = {}

            for count, w in enumerate(wi):
                F123_x[w] = fock_list[3 * count]
                F123_y[w] = fock_list[3 * count + 1]
                F123_z[w] = fock_list[3 * count + 2]

            fock_dict = {'F123_x': F123_x, 'F123_y': F123_y, 'F123_z': F123_z}

        return fock_dict

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

        F0_a = fo['F0']

        for w in wi:
            # Response
            k_x_ = kX['Nc'][('x', -w)].T
            k_y_ = kX['Nc'][('y', -w)].T
            k_z_ = kX['Nc'][('z', -w)].T

            k_x = kX['Nb'][('x', w)].T
            k_y = kX['Nb'][('y', w)].T
            k_z = kX['Nb'][('z', w)].T

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

            f_x_ = fo['Fc'][('x', -w)]
            f_y_ = fo['Fc'][('y', -w)]
            f_z_ = fo['Fc'][('z', -w)]

            f_x = fo['Fb'][('x', w)]
            f_y = fo['Fb'][('y', w)]
            f_z = fo['Fb'][('z', w)]

            f_sig_xx = np.conjugate(fo2[(('N_sig_xx', w), 2 * w)]).T
            f_sig_yy = np.conjugate(fo2[(('N_sig_yy', w), 2 * w)]).T
            f_sig_zz = np.conjugate(fo2[(('N_sig_zz', w), 2 * w)]).T

            f_sig_xy = np.conjugate(fo2[(('N_sig_xy', w), 2 * w)]).T
            f_sig_xz = np.conjugate(fo2[(('N_sig_xz', w), 2 * w)]).T
            f_sig_yz = np.conjugate(fo2[(('N_sig_yz', w), 2 * w)]).T

            f_lamtau_xx = np.conjugate(fo2[(('N_lamtau_xx', w), 0)]).T
            f_lamtau_yy = np.conjugate(fo2[(('N_lamtau_yy', w), 0)]).T
            f_lamtau_zz = np.conjugate(fo2[(('N_lamtau_zz', w), 0)]).T

            f_lamtau_xy = np.conjugate(fo2[(('N_lamtau_xy', w), 0)]).T
            f_lamtau_xz = np.conjugate(fo2[(('N_lamtau_xz', w), 0)]).T
            f_lamtau_yz = np.conjugate(fo2[(('N_lamtau_yz', w), 0)]).T

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

            X_terms = (zeta_sig_xx + zeta_sig_xy + zeta_sig_xz).T + (
                zeta_lamtau_xx + zeta_lamtau_xy +
                zeta_lamtau_xz).T + (0.5 * fo2['F123_x'][w]).T
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

            Y_terms = (zeta_sig_yx + zeta_sig_yy + zeta_sig_yz).T + (
                zeta_lamtau_yx + zeta_lamtau_yy +
                zeta_lamtau_yz).T + (0.5 * fo2['F123_y'][w]).T
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

            Z_terms = (zeta_sig_zx + zeta_sig_zy + zeta_sig_zz).T + (
                zeta_lamtau_zx + zeta_lamtau_zy +
                zeta_lamtau_zz).T + (0.5 * fo2['F123_z'][w]).T
            Ff_z = -2 * LinearSolver.lrmat2vec(Z_terms, nocc, norb)
            Ff_z = self.anti_sym(Ff_z)
            f_iso_z[w] = Ff_z

        return {'f_iso_x': f_iso_x, 'f_iso_y': f_iso_y, 'f_iso_z': f_iso_z}

    def get_other_terms(self, wi, track, n_x, n_xy, X, kX, kXY, da, nocc, norb):
        """
        Computes the terms involving X[3],A[3],X[2],A[2] in the isotropic cubic
        response function

        :param wi:
            A list containing all the frequencies
        :param track:
            A list that contains information about what γ components that are
            to be computed and which freqs
        :param n_x:
            A dictonary containing all the single-index response vectors
        :param n_xy:
            A dictonary containing all the two-index response vectors
        :param X:
            A dictonray with all the property integral matricies
        :param kX:
            A dictonary with all the respone matricies
        :param kXY:
            A dictonary containing all the two-index response matricies
        :param da:
            The SCF density matrix in MO bassi
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

        for j in range(len(wi)):
            na_x3_ny_nz = 0
            na_a3_nx_ny = 0

            for i in range(j * comp_per_freq, (j + 1) * comp_per_freq):
                comp_i = track[i]

                vals = comp_i.split(',')
                w1 = float(vals[1])
                w2 = float(vals[2])
                w3 = float(vals[3])

                Na = n_x['Na'][(comp_i[0], w1)]
                Nb = n_x['Nb'][(comp_i[1], w1)]

                Nc = n_x['Nc'][(comp_i[2], w2)]
                Nd = n_x['Nd'][(comp_i[3], w3)]

                kb = kX['Nb'][(comp_i[1], w1)]
                kc = kX['Nc'][(comp_i[2], w2)]
                kd = kX['Nd'][(comp_i[3], w3)]

                A = X[comp_i[0]]
                B = X[comp_i[1]]
                C = X[comp_i[2]]
                D = X[comp_i[3]]

                # Na X[3]NyNz

                na_x3_ny_nz += -np.matmul(
                    Na.T, self.x3_contract(kc, kd, B, da, nocc, norb))
                na_x3_ny_nz += -np.matmul(
                    Na.T, self.x3_contract(kd, kc, B, da, nocc, norb))
                na_x3_ny_nz += -np.matmul(
                    Na.T, self.x3_contract(kd, kb, C, da, nocc, norb))
                na_x3_ny_nz += -np.matmul(
                    Na.T, self.x3_contract(kb, kd, C, da, nocc, norb))
                na_x3_ny_nz += -np.matmul(
                    Na.T, self.x3_contract(kb, kc, D, da, nocc, norb))
                na_x3_ny_nz += -np.matmul(
                    Na.T, self.x3_contract(kc, kb, D, da, nocc, norb))

                # NaA[3]n_xNy

                na_a3_nx_ny += np.matmul(
                    self.a3_contract(kb, kc, A, da, nocc, norb), Nd)
                na_a3_nx_ny += np.matmul(
                    self.a3_contract(kb, kd, A, da, nocc, norb), Nc)
                na_a3_nx_ny += np.matmul(
                    self.a3_contract(kc, kb, A, da, nocc, norb), Nd)
                na_a3_nx_ny += np.matmul(
                    self.a3_contract(kc, kd, A, da, nocc, norb), Nb)
                na_a3_nx_ny += np.matmul(
                    self.a3_contract(kd, kb, A, da, nocc, norb), Nc)
                na_a3_nx_ny += np.matmul(
                    self.a3_contract(kd, kc, A, da, nocc, norb), Nb)

            na_a3_nx_ny_dict[(wi[j], -wi[j], wi[j])] = (1. / 15) * na_a3_nx_ny
            na_x3_ny_nz_dict[(wi[j], -wi[j], wi[j])] = (1. / 15) * na_x3_ny_nz

        for i in range(len(wi)):
            vals = track[i * comp_per_freq].split(',')
            w = float(vals[1])
            wa = float(vals[1])
            wb = float(vals[1])
            wc = float(vals[2])
            wd = float(vals[3])

            wcd = 0
            wbd = wb + wd

            na_x2_nyz = 0
            nx_a2_nyz = 0

            # CD
            # x
            kcd = kXY[(('N_lamtau_xx', w), wcd)]
            Ncd = n_xy[(('N_lamtau_xx', w), wcd)]

            Na = n_x['Na'][('x', wa)]
            Nb = n_x['Nb'][('x', w)]
            kb = kX['Nb'][('x', w)]
            A = X['x']
            B = X['x']
            na_x2_nyz += np.matmul(Na.T,
                                   self.x2_contract(kcd, B, da, nocc, norb))

            nx_a2_nyz += np.matmul(self.a2_contract(kb, A, da, nocc, norb), Ncd)
            nx_a2_nyz += np.matmul(self.a2_contract(kcd, A, da, nocc, norb), Nb)

            kcd = kXY[(('N_lamtau_xy', w), wcd)]
            Ncd = n_xy[(('N_lamtau_xy', w), wcd)]

            Nb = n_x['Nb'][('y', w)]
            kb = kX['Nb'][('y', w)]
            Na = n_x['Na'][('x', wa)]
            A = X['x']
            B = X['y']

            na_x2_nyz += np.matmul(Na.T,
                                   self.x2_contract(kcd, B, da, nocc, norb))

            nx_a2_nyz += np.matmul(self.a2_contract(kb, A, da, nocc, norb), Ncd)
            nx_a2_nyz += np.matmul(self.a2_contract(kcd, A, da, nocc, norb), Nb)

            kcd = kXY[(('N_lamtau_xz', w), wcd)]
            Ncd = n_xy[(('N_lamtau_xz', w), wcd)]
            Nb = n_x['Nb'][('z', w)]
            kb = kX['Nb'][('z', w)]
            Na = n_x['Na'][('x', wa)]
            A = X['x']
            B = X['z']

            na_x2_nyz += np.matmul(Na.T,
                                   self.x2_contract(kcd, B, da, nocc, norb))

            nx_a2_nyz += np.matmul(self.a2_contract(kb, A, da, nocc, norb), Ncd)
            nx_a2_nyz += np.matmul(self.a2_contract(kcd, A, da, nocc, norb), Nb)

            # y

            kcd = kXY[(('N_lamtau_xy', w), wcd)]
            Ncd = n_xy[(('N_lamtau_xy', w), wcd)]
            Nb = n_x['Nb'][('x', w)]
            kb = kX['Nb'][('x', w)]
            A = X['y']
            B = X['x']
            Na = n_x['Na'][('y', wa)]

            na_x2_nyz += np.matmul(Na.T,
                                   self.x2_contract(kcd, B, da, nocc, norb))

            nx_a2_nyz += np.matmul(self.a2_contract(kb, A, da, nocc, norb), Ncd)
            nx_a2_nyz += np.matmul(self.a2_contract(kcd, A, da, nocc, norb), Nb)

            kcd = kXY[(('N_lamtau_yy', w), wcd)]
            Ncd = n_xy[(('N_lamtau_yy', w), wcd)]
            Na = n_x['Na'][('y', wa)]
            Nb = n_x['Nb'][('y', w)]
            kb = kX['Nb'][('y', w)]
            A = X['y']
            B = X['y']

            na_x2_nyz += np.matmul(Na.T,
                                   self.x2_contract(kcd, B, da, nocc, norb))

            nx_a2_nyz += np.matmul(self.a2_contract(kb, A, da, nocc, norb), Ncd)
            nx_a2_nyz += np.matmul(self.a2_contract(kcd, A, da, nocc, norb), Nb)

            kcd = kXY[(('N_lamtau_yz', w), wcd)]
            Ncd = n_xy[(('N_lamtau_yz', w), wcd)]
            Nb = n_x['Nb'][('z', w)]
            kb = kX['Nb'][('z', w)]

            A = X['y']
            B = X['z']
            Na = n_x['Na'][('y', wa)]

            na_x2_nyz += np.matmul(Na.T,
                                   self.x2_contract(kcd, B, da, nocc, norb))

            nx_a2_nyz += np.matmul(self.a2_contract(kb, A, da, nocc, norb), Ncd)
            nx_a2_nyz += np.matmul(self.a2_contract(kcd, A, da, nocc, norb), Nb)

            # z

            kcd = kXY[(('N_lamtau_xz', w), wcd)]
            Ncd = n_xy[(('N_lamtau_xz', w), wcd)]
            Nb = n_x['Nb'][('x', w)]
            kb = kX['Nb'][('x', w)]
            Na = n_x['Na'][('z', wa)]
            A = X['z']
            B = X['x']

            na_x2_nyz += np.matmul(Na.T,
                                   self.x2_contract(kcd, B, da, nocc, norb))

            nx_a2_nyz += np.matmul(self.a2_contract(kb, A, da, nocc, norb), Ncd)
            nx_a2_nyz += np.matmul(self.a2_contract(kcd, A, da, nocc, norb), Nb)

            kcd = kXY[(('N_lamtau_yz', w), wcd)]
            Ncd = n_xy[(('N_lamtau_yz', w), wcd)]
            Na = n_x['Na'][('z', wa)]
            Nb = n_x['Nb'][('y', w)]
            kb = kX['Nb'][('y', w)]
            A = X['z']
            B = X['y']

            na_x2_nyz += np.matmul(Na.T,
                                   self.x2_contract(kcd, B, da, nocc, norb))

            nx_a2_nyz += np.matmul(self.a2_contract(kb, A, da, nocc, norb), Ncd)
            nx_a2_nyz += np.matmul(self.a2_contract(kcd, A, da, nocc, norb), Nb)

            kcd = kXY[(('N_lamtau_zz', w), wcd)]
            Ncd = n_xy[(('N_lamtau_zz', w), wcd)]
            Nb = n_x['Nb'][('z', w)]
            Na = n_x['Na'][('z', wa)]
            kb = kX['Nb'][('z', w)]
            A = X['z']
            B = X['z']

            na_x2_nyz += np.matmul(Na.T,
                                   self.x2_contract(kcd, B, da, nocc, norb))

            nx_a2_nyz += np.matmul(self.a2_contract(kb, A, da, nocc, norb), Ncd)
            nx_a2_nyz += np.matmul(self.a2_contract(kcd, A, da, nocc, norb), Nb)

            # BD

            kbd = kXY[(('N_sig_xx', w), wbd)]
            Nbd = n_xy[(('N_sig_xx', w), wbd)]
            Na = n_x['Na'][('x', wa)]
            Nc = n_x['Nc'][('x', wc)]
            kc = kX['Nc'][('x', wc)]
            A = X['x']
            C = X['x']

            na_x2_nyz += np.matmul(Na.T,
                                   self.x2_contract(kbd, C, da, nocc, norb))

            nx_a2_nyz += np.matmul(self.a2_contract(kc, A, da, nocc, norb), Nbd)
            nx_a2_nyz += np.matmul(self.a2_contract(kbd, A, da, nocc, norb), Nc)

            kbd = kXY[(('N_sig_xy', w), wbd)]
            Nbd = n_xy[(('N_sig_xy', w), wbd)]
            Na = n_x['Na'][('x', wa)]
            Nc = n_x['Nc'][('y', wc)]
            kc = kX['Nc'][('y', wc)]
            A = X['x']
            C = X['y']

            na_x2_nyz += np.matmul(Na.T,
                                   self.x2_contract(kbd, C, da, nocc, norb))

            nx_a2_nyz += np.matmul(self.a2_contract(kc, A, da, nocc, norb), Nbd)
            nx_a2_nyz += np.matmul(self.a2_contract(kbd, A, da, nocc, norb), Nc)

            kbd = kXY[(('N_sig_xz', w), wbd)]
            Nbd = n_xy[(('N_sig_xz', w), wbd)]
            Nc = n_x['Nc'][('z', wc)]
            kc = kX['Nc'][('z', wc)]
            Na = n_x['Na'][('x', wa)]
            A = X['x']
            C = X['z']

            na_x2_nyz += np.matmul(Na.T,
                                   self.x2_contract(kbd, C, da, nocc, norb))

            nx_a2_nyz += np.matmul(self.a2_contract(kc, A, da, nocc, norb), Nbd)
            nx_a2_nyz += np.matmul(self.a2_contract(kbd, A, da, nocc, norb), Nc)

            # y

            kbd = kXY[(('N_sig_xy', w), wbd)]
            Nbd = n_xy[(('N_sig_xy', w), wbd)]
            Na = n_x['Na'][('y', wa)]
            Nc = n_x['Nc'][('x', wc)]
            kc = kX['Nc'][('x', wc)]
            A = X['y']
            C = X['x']

            na_x2_nyz += np.matmul(Na.T,
                                   self.x2_contract(kbd, C, da, nocc, norb))

            nx_a2_nyz += np.matmul(self.a2_contract(kc, A, da, nocc, norb), Nbd)
            nx_a2_nyz += np.matmul(self.a2_contract(kbd, A, da, nocc, norb), Nc)

            kbd = kXY[(('N_sig_yy', w), wbd)]
            Nbd = n_xy[(('N_sig_yy', w), wbd)]
            Nc = n_x['Nc'][('y', wc)]
            Na = n_x['Na'][('y', wa)]
            kc = kX['Nc'][('y', wc)]
            A = X['y']
            C = X['y']

            na_x2_nyz += np.matmul(Na.T,
                                   self.x2_contract(kbd, C, da, nocc, norb))

            nx_a2_nyz += np.matmul(self.a2_contract(kc, A, da, nocc, norb), Nbd)
            nx_a2_nyz += np.matmul(self.a2_contract(kbd, A, da, nocc, norb), Nc)

            kbd = kXY[(('N_sig_yz', w), wbd)]
            Nbd = n_xy[(('N_sig_yz', w), wbd)]
            Nc = n_x['Nc'][('z', wc)]
            Na = n_x['Na'][('y', wa)]
            kc = kX['Nc'][('z', wc)]
            A = X['y']
            C = X['z']

            na_x2_nyz += np.matmul(Na.T,
                                   self.x2_contract(kbd, C, da, nocc, norb))

            nx_a2_nyz += np.matmul(self.a2_contract(kc, A, da, nocc, norb), Nbd)
            nx_a2_nyz += np.matmul(self.a2_contract(kbd, A, da, nocc, norb), Nc)

            # z

            kbd = kXY[(('N_sig_xz', w), wbd)]
            Nbd = n_xy[(('N_sig_xz', w), wbd)]
            Nc = n_x['Nc'][('x', wc)]
            kc = kX['Nc'][('x', wc)]
            Na = n_x['Na'][('z', wa)]
            A = X['z']
            C = X['x']

            na_x2_nyz += np.matmul(Na.T,
                                   self.x2_contract(kbd, C, da, nocc, norb))

            nx_a2_nyz += np.matmul(self.a2_contract(kc, A, da, nocc, norb), Nbd)
            nx_a2_nyz += np.matmul(self.a2_contract(kbd, A, da, nocc, norb), Nc)

            kbd = kXY[(('N_sig_yz', w), wbd)]
            Nbd = n_xy[(('N_sig_yz', w), wbd)]
            Nc = n_x['Nc'][('y', wc)]
            kc = kX['Nc'][('y', wc)]
            Na = n_x['Na'][('z', wa)]
            A = X['z']
            C = X['y']

            na_x2_nyz += np.matmul(Na.T,
                                   self.x2_contract(kbd, C, da, nocc, norb))

            nx_a2_nyz += np.matmul(self.a2_contract(kc, A, da, nocc, norb), Nbd)
            nx_a2_nyz += np.matmul(self.a2_contract(kbd, A, da, nocc, norb), Nc)

            kbd = kXY[(('N_sig_zz', w), wbd)]
            Nbd = n_xy[(('N_sig_zz', w), wbd)]
            Nc = n_x['Nc'][('z', wc)]
            Na = n_x['Na'][('z', wa)]
            kc = kX['Nc'][('z', wc)]
            A = X['z']
            C = X['z']

            na_x2_nyz += np.matmul(Na.T,
                                   self.x2_contract(kbd, C, da, nocc, norb))

            nx_a2_nyz += np.matmul(self.a2_contract(kc, A, da, nocc, norb), Nbd)
            nx_a2_nyz += np.matmul(self.a2_contract(kbd, A, da, nocc, norb), Nc)

            na_x2_nyz_dict[(w, -w, w)] = -(1. / 15) * na_x2_nyz
            nx_a2_nyz_dict[(w, -w, w)] = -(1. / 15) * nx_a2_nyz

        return {
            'NaX3NyNz': na_x3_ny_nz_dict,
            'NaA3NxNy': na_a3_nx_ny_dict,
            'NaX2Nyz': na_x2_nyz_dict,
            'NxA2Nyz': nx_a2_nyz_dict,
        }

    def get_t4(self, wi, e4_dict, n_x, kX, track, da, nocc, norb):
        """
        Computes the contraction of the E[4] tensor with that of the S[4] and
        R[4] tensors to return the contraction of T[4] as a dictonary of
        vectors. T[4]n_xNyNz = (E^[4]-ω_1S^[4]-ω_1S^[4]-ω_3S^[4]-γiR^[4])

        :param wi:
            A list of all the freqs
        :param e4_dict:
            A dictonary of all the E[4] contraction
        :param n_x:
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
            R4term = self.get_r4(wi, kX, n_x, track, da, nocc, norb)

        comp_per_freq = len(track) // len(wi)

        for i in range(len(wi)):
            vals = track[i * comp_per_freq].split(',')
            w = float(vals[1])
            ww = float(vals[1])

            t4term = (np.matmul(n_x['Na'][('x', w)],
                                e4_dict['f_iso_x'][ww] - S4[('x', ww)]) +
                      np.matmul(n_x['Na'][('y', w)],
                                e4_dict['f_iso_y'][ww] - S4[('y', ww)]) +
                      np.matmul(n_x['Na'][('z', w)],
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
                kC = kX['Nc'][(comp_i[2], w2)]
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

    def get_r4(self, freqs, kX, n_x, track, d_a_mo, nocc, norb):
        """
        Returns a dict with all the R[4]NxNyNz contractions for the subsequent
        T[4] contraction

        :param freqs:
            A list of all the frequencies
        :param kX:
            A dictonary of all the first-order response matrices
        :param n_x:
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

                # Na = n_x['Na'][(comp_i[0], w_s)]
                Nb = n_x['Nb'][(comp_i[1], w1)]
                Nc = n_x['Nc'][(comp_i[2], w2)]
                Nd = n_x['Nd'][(comp_i[3], w3)]
                kA = kX['Na'][(comp_i[0], w_s)]
                kB = kX['Nb'][(comp_i[1], w1)]
                kC = kX['Nc'][(comp_i[2], w2)]
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

        w_str = 'Gamma Tensor Components at Given Frequencies'
        self.ostream.print_blank()
        self.ostream.print_header(w_str.ljust(width))
        self.ostream.print_blank()

        for w in freqs:
            title = '{:<7s} {:>10s} {:>10s} {:>16s}'.format(
                'Contribution', 'Frequency', 'Real', 'Imaginary')
            self.ostream.print_header(title.ljust(width))
            self.ostream.print_header(('-' * len(title)).ljust(width))

            cont_label = "ΣNaT3NxNyz  {:10.4f}".format(w)

            w_str = '{:<15s} {:15.8f} {:13.8f}j'.format(cont_label,
                                                        t3_dict[w, -w, w].real,
                                                        t3_dict[w, -w, w].imag)
            self.ostream.print_header(w_str.ljust(width))

            cont_label = "ΣNaT4NxNyNz  {:9.4f}".format(w)

            w_str = '{:<15s} {:15.8f} {:13.8f}j'.format(cont_label,
                                                        t4_dict[w, -w, w].real,
                                                        t4_dict[w, -w, w].imag)
            self.ostream.print_header(w_str.ljust(width))

            cont_label = "ΣNaX2Nyz  {:12.4f}".format(w)

            w_str = '{:<15s} {:15.8f} {:13.8f}j'.format(cont_label,
                                                        NaX2Nyz[w, -w, w].real,
                                                        NaX2Nyz[w, -w, w].imag)
            self.ostream.print_header(w_str.ljust(width))

            cont_label = "ΣNaX3NyNz  {:11.4f}".format(w)

            w_str = '{:<15s} {:15.8f} {:13.8f}j'.format(cont_label,
                                                        NaX3NyNz[w, -w, w].real,
                                                        NaX3NyNz[w, -w, w].imag)
            self.ostream.print_header(w_str.ljust(width))

            cont_label = "ΣNxA2Nyz  {:12.4f}".format(w)

            w_str = '{:<15s} {:15.8f} {:13.8f}j'.format(cont_label,
                                                        NxA2Nyz[w, -w, w].real,
                                                        NxA2Nyz[w, -w, w].imag)
            self.ostream.print_header(w_str.ljust(width))

            cont_label = "ΣNaA3NxNy  {:11.4f}".format(w)

            w_str = '{:<15s} {:15.8f} {:13.8f}j'.format(cont_label,
                                                        NaA3NxNy[w, -w, w].real,
                                                        NaA3NxNy[w, -w, w].imag)
            self.ostream.print_header(w_str.ljust(width))

            cont_label = "Σ<<μ;μ,μ,μ>>  {:8.4f}".format(w)

            w_str = '{:<15s} {:15.8f} {:13.8f}j'.format(cont_label,
                                                        gamma[w, -w, w].real,
                                                        gamma[w, -w, w].imag)
            self.ostream.print_header(w_str.ljust(width))
            self.ostream.print_blank()

        self.ostream.print_blank()
