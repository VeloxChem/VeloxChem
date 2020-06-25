import numpy as np
import ctypes
import psutil
import time
import os

from .veloxchemlib import ElectronRepulsionIntegralsDriver
from .veloxchemlib import KineticEnergyIntegralsDriver
from .veloxchemlib import NuclearPotentialIntegralsDriver
from .veloxchemlib import mpi_master
from .veloxchemlib import denmat
from .veloxchemlib import fockmat
from .qqscheme import get_qq_scheme
from .cppsolver import ComplexResponse
from .linearsolver import LinearSolver
from .aofockmatrix import AOFockMatrix
from .aodensitymatrix import AODensityMatrix


class TPAdriver:
    """
    Computes the third order-gradients used for the isotropic cubic response
    function, also contains some other methods used for the evaluation of some
    relevant quantities for TPA calculations, such as A[3] and
    A[2] contractions

    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.
    """

    def __init__(self, comm, ostream):

        self.comm = comm
        self.rank = self.comm.Get_rank()
        self.nodes = self.comm.Get_size()

        self.ostream = ostream

        self.batch_size = None

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

        if 'damping' in rsp_dict:
            self.damping = float(rsp_dict['damping'])
        if 'lindep_thresh' in rsp_dict:
            self.lindep_thresh = float(rsp_dict['lindep_thresh'])
        if 'conv_thresh' in rsp_dict:
            self.conv_thresh = float(rsp_dict['conv_thresh'])
        if 'max_iter' in rsp_dict:
            self.max_iter = int(rsp_dict['max_iter'])
        if 'eri_thresh' in rsp_dict:
            self.eri_thresh = float(rsp_dict['eri_thresh'])
        if 'qq_type' in rsp_dict:
            self.qq_type = rsp_dict['qq_type']
        if 'batch_size' in rsp_dict:
            self.batch_size = int(rsp_dict['batch_size'])

    def get_e4(self, wi, kX, fo, nocc, norb):
        """
        Contracts E[4]n_xNyNz for the isotropic cubic response function. It
        takes the Fock matrices from fock_dict and contracts them with the
        response vectors.

        : param wi - A list of freqs
        : param kX - A dict of the single index response matricies
        : param fo - A dictonary of transformed Fock matricies from fock_dict
        : param nocc - The number of occupied orbitals
        : param norb - The total number of orbitals

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
            #   Φ_αα^σ = φ(κα,κα,Fα,F0) + φ(κα,κα,Fα,F0) + Σ_{ρ}^{x,y,z}[φ(κρ,κρ,Fρ,F0)] for α=β
            #   Φ_αβ^σ = φ(κa,κa,Fb,F0) + φ(κb,κb,Fb,F0) for α≠β
            #  For the Φ_{αβ}^{λ+τ} component see article.

            Phi_sig_xx = 2 * (3 * self.phi(kx, kx, Fx, F0) + self.phi(
                ky, ky, Fy, F0) + self.phi(kz, kz, Fz, F0))

            Phi_sig_yy = 2 * (self.phi(kx, kx, Fx, F0) +
                              3 * self.phi(ky, ky, Fy, F0) +
                              self.phi(kz, kz, Fz, F0))

            Phi_sig_zz = 2 * (self.phi(kx, kx, Fx, F0) + self.phi(
                ky, ky, Fy, F0) + 3 * self.phi(kz, kz, Fz, F0))

            Phi_sig_xy = 2 * self.phi(kx, ky, Fy, F0) + 2 * self.phi(
                ky, kx, Fx, F0)

            Phi_sig_xz = 2 * self.phi(kx, kz, Fz, F0) + 2 * self.phi(
                kz, kx, Fx, F0)

            Phi_sig_yz = 2 * self.phi(ky, kz, Fz, F0) + 2 * self.phi(
                kz, ky, Fy, F0)

            Phi_lamtau_xx = 3 * self.phi(kx, kx_, Fx_, F0) + 3 * self.phi(
                kx_, kx, Fx, F0) + self.phi(ky, ky_, Fy_, F0) + self.phi(
                    ky_, ky, Fy, F0) + self.phi(kz, kz_, Fz_, F0) + self.phi(
                        kz_, kz, Fz, F0)

            Phi_lamtau_xx = 2 * Phi_lamtau_xx

            Phi_lamtau_yy = self.phi(kx, kx_, Fx_, F0) + self.phi(
                kx_, kx, Fx,
                F0) + 3 * self.phi(ky, ky_, Fy_, F0) + 3 * self.phi(
                    ky_, ky, Fy, F0) + self.phi(kz, kz_, Fz_, F0) + self.phi(
                        kz_, kz, Fz, F0)

            Phi_lamtau_yy = 2 * Phi_lamtau_yy

            Phi_lamtau_zz = self.phi(kx, kx_, Fx_, F0) + self.phi(
                kx_, kx, Fx, F0) + self.phi(ky, ky_, Fy_, F0) + self.phi(
                    ky_, ky, Fy, F0) + 3 * self.phi(
                        kz, kz_, Fz_, F0) + 3 * self.phi(kz_, kz, Fz, F0)

            Phi_lamtau_zz = 2 * Phi_lamtau_zz

            Phi_lamtau_xy = self.phi(kx_, ky, Fy, F0) + self.phi(
                ky, kx_, Fx_, F0) + self.phi(kx, ky_, Fy_, F0) + self.phi(
                    ky_, kx, Fx, F0)

            Phi_lamtau_xy = 2 * Phi_lamtau_xy

            Phi_lamtau_xz = self.phi(kx_, kz, Fz, F0) + self.phi(
                kz, kx_, Fx_, F0) + self.phi(kx, kz_, Fz_, F0) + self.phi(
                    kz_, kx, Fx, F0)

            Phi_lamtau_xz = 2 * Phi_lamtau_xz

            Phi_lamtau_yz = self.phi(ky_, kz, Fz, F0) + self.phi(
                kz, ky_, Fy_, F0) + self.phi(ky, kz_, Fz_, F0) + self.phi(
                    kz_, ky, Fy, F0)
            Phi_lamtau_yz = 2 * Phi_lamtau_yz

            # Computess all the elements of the Fock vector formed from the
            # E[4] contraction as E[4]NxNyNz = [f_{is} // f_{si}], for the
            # isotropic case, as derived in the article
            # The elements of f_{is} for each spatial component α is given by
            # an expression of the form
            # f_α = Σ_{β}^{x,y,z} [κ_{β}^{ω},Φ_{αβ}^{λ+τ}+f_{αβ}^{λ+τ}] + [κ_{β}^{-ω},Φ_{αβ}^{σ}+f_{αβ}^{σ}]

            # x

            # Creating the transformed total Fock matrices

            f_x = fo['F123_x'][w]

            f_x += self.commut(kx, Phi_lamtau_xx + f_lamtau_xx) + self.commut(
                ky, Phi_lamtau_xy + f_lamtau_xy) + self.commut(
                    kz, Phi_lamtau_xz + f_lamtau_xz)

            f_x += self.commut(kx_, Phi_sig_xx + f_sig_xx) + self.commut(
                ky_, Phi_sig_xy + f_sig_xy) + self.commut(
                    kz_, Phi_sig_xz + f_sig_xz)

            # Taking the non redundant matrix elements {i,s} and forming the
            # anti-symmetric Fock vector
            f_x = -2 * 1 / 6 * LinearSolver.lrmat2vec(f_x.T, nocc, norb)
            f_x = np.array(self.anti_sym(f_x))
            f_iso_x.update({w: f_x})

            # y

            # Creating the transformed total Fock matrices
            f_y = fo['F123_y'][w]

            f_y += self.commut(kx, Phi_lamtau_xy + f_lamtau_xy) + self.commut(
                ky, Phi_lamtau_yy + f_lamtau_yy) + self.commut(
                    kz, Phi_lamtau_yz + f_lamtau_yz)

            f_y += self.commut(kx_, Phi_sig_xy + f_sig_xy) + self.commut(
                ky_, Phi_sig_yy + f_sig_yy) + self.commut(
                    kz_, Phi_sig_yz + f_sig_yz)

            # Taking the non redundant matrix elements {i,s} and forming the
            # anti-symmetric Fock vector
            f_y = -2 * 1 / 6 * LinearSolver.lrmat2vec(f_y.T, nocc, norb)
            f_y = np.array(self.anti_sym(f_y))
            f_iso_y.update({w: f_y})

            # z

            # Creating the transformed total Fock matrices
            f_z = fo['F123_z'][w]

            f_z += self.commut(kx, Phi_lamtau_xz + f_lamtau_xz) + self.commut(
                ky, Phi_lamtau_yz + f_lamtau_yz) + self.commut(
                    kz, Phi_lamtau_zz + f_lamtau_zz)

            f_z += self.commut(kx_, Phi_sig_xz + f_sig_xz) + self.commut(
                ky_, Phi_sig_yz + f_sig_yz) + self.commut(
                    kz_, Phi_sig_zz + f_sig_zz)

            # Taking the non redundant matrix elements {i,s} and forming the
            # anti-symmetric Fock vector
            f_z = -2 * 1 / 6 * LinearSolver.lrmat2vec(f_z.T, nocc, norb)
            f_z = np.array(self.anti_sym(f_z))
            f_iso_z.update({w: f_z})

        e4_dict = {'f_iso_x': f_iso_x, 'f_iso_y': f_iso_y, 'f_iso_z': f_iso_z}
        return e4_dict

    def get_e3_red(self, wi, kX, kXY, fo, fo2, nocc, norb):
        """
        This code contracts E[3]n_xNyz for the isotropic cubic response
        function. It takes the Fock matrices from fock_dict and fock_dict_II
        and contracts them with the first and second-order response vectors.

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

        F0_a = fo['F0']

        for w in wi:
            # Response #
            k_x_ = kX['Nc'][('x', -w)].T
            k_y_ = kX['Nc'][('y', -w)].T
            k_z_ = kX['Nc'][('z', -w)].T

            k_sig_xx = kXY[(('N_sig_xx', w), 2 * w)].T
            k_sig_yy = kXY[(('N_sig_yy', w), 2 * w)].T
            k_sig_zz = kXY[(('N_sig_zz', w), 2 * w)].T

            k_sig_xy = kXY[(('N_sig_xy', w), 2 * w)].T
            k_sig_xz = kXY[(('N_sig_xz', w), 2 * w)].T
            k_sig_yz = kXY[(('N_sig_yz', w), 2 * w)].T

            # Focks #
            f_x_ = fo['Fc'][('x', -w)]
            f_y_ = fo['Fc'][('y', -w)]
            f_z_ = fo['Fc'][('z', -w)]

            f_sig_xx = np.conjugate(fo2[(('N_sig_xx', w), 2 * w)]).T
            f_sig_yy = np.conjugate(fo2[(('N_sig_yy', w), 2 * w)]).T
            f_sig_zz = np.conjugate(fo2[(('N_sig_zz', w), 2 * w)]).T

            f_sig_xy = np.conjugate(fo2[(('N_sig_xy', w), 2 * w)]).T
            f_sig_xz = np.conjugate(fo2[(('N_sig_xz', w), 2 * w)]).T
            f_sig_yz = np.conjugate(fo2[(('N_sig_yz', w), 2 * w)]).T

            # x #
            zeta_sig_xx = self.xi(k_x_, k_sig_xx, f_x_, f_sig_xx, F0_a)
            zeta_sig_yy = self.xi(k_x_, k_sig_yy, f_x_, f_sig_yy, F0_a)
            zeta_sig_zz = self.xi(k_x_, k_sig_zz, f_x_, f_sig_zz, F0_a)

            zeta_sig_xy = self.xi(k_y_, k_sig_xy, f_y_, f_sig_xy, F0_a)
            zeta_sig_xz = self.xi(k_z_, k_sig_xz, f_z_, f_sig_xz, F0_a)

            X_terms = (zeta_sig_xx + zeta_sig_xy +
                       zeta_sig_xz).T + (0.5 * fo2['F123_x'][w]).T
            Ff_x = -2 * LinearSolver.lrmat2vec(X_terms, nocc, norb)
            Ff_x = np.array(self.anti_sym(Ff_x))
            f_iso_x.update({w: Ff_x})

            # y #
            zeta_sig_yx = self.xi(k_x_, k_sig_xy, f_x_, f_sig_xy, F0_a)
            zeta_sig_yy = self.xi(k_y_, k_sig_yy, f_y_, f_sig_yy, F0_a)
            zeta_sig_yz = self.xi(k_z_, k_sig_yz, f_z_, f_sig_yz, F0_a)

            Y_terms = (zeta_sig_yx + zeta_sig_yy +
                       zeta_sig_yz).T + (0.5 * fo2['F123_y'][w]).T
            Ff_y = -2 * LinearSolver.lrmat2vec(Y_terms, nocc, norb)
            Ff_y = np.array(self.anti_sym(Ff_y))
            f_iso_y.update({w: Ff_y})

            # z #
            zeta_sig_zx = self.xi(k_x_, k_sig_xz, f_x_, f_sig_xz, F0_a)
            zeta_sig_zy = self.xi(k_y_, k_sig_yz, f_y_, f_sig_yz, F0_a)
            zeta_sig_zz = self.xi(k_z_, k_sig_zz, f_z_, f_sig_zz, F0_a)

            Z_terms = (zeta_sig_zx + zeta_sig_zy +
                       zeta_sig_zz).T + (0.5 * fo2['F123_z'][w]).T
            Ff_z = -2 * LinearSolver.lrmat2vec(Z_terms, nocc, norb)
            Ff_z = np.array(self.anti_sym(Ff_z))
            f_iso_z.update({w: Ff_z})

        e3_dict = {'f_iso_x': f_iso_x, 'f_iso_y': f_iso_y, 'f_iso_z': f_iso_z}
        return e3_dict

    def get_e3(self, wi, kX, kXY, fo, fo2, nocc, norb):
        """
        This code contracts E[3]

        : param iso - A boolean that specifies if the progam is to compute the isotropic γ or a user specief combination of components
        : param wi - A list of freqs
        : param keys - A dict of lists of keys that give information about what components are present
        : param track - A list containing information about what γ components that are to be computed
        : param kX - A dict of the single index response matricies
        : param kXY - A dict of the two index response matrices
        : param fo - A dictonary of transformed Fock matricies from fock_dict
        : param fo2 - A dictonarty of transfromed Fock matricies from fock_dict_two
        : param nocc - The number of occupied orbitals
        : param norb - The total number of orbitals

        :return:
            A dictonary of compounded E[3] tensors for the isotropic cubic response function for TPA
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
            Ff_x = np.array(self.anti_sym(Ff_x))
            f_iso_x.update({w: Ff_x})

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
            Ff_y = np.array(self.anti_sym(Ff_y))
            f_iso_y.update({w: Ff_y})

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
            Ff_z = np.array(self.anti_sym(Ff_z))
            f_iso_z.update({w: Ff_z})

        e3_dict = {'f_iso_x': f_iso_x, 'f_iso_y': f_iso_y, 'f_iso_z': f_iso_z}
        return e3_dict

    def get_densities_red(self, wi, kX, S, D0, mo, nocc, norb):
        """
        Computes the compounded densities needed for the compounded Fock
        matrics F^{σ} used for the reduced iostropic cubic response function

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
        :param nocc:
            The number of occupied orbitals
        :param norb:
            The number of total orbitals

        :return:
            A list of tranformed compounded densities
        """

        density_list = []

        for w in wi:
            # convert response matrix to ao basis #
            kx = self.mo2ao(mo, kX['Nb'][('x', w)]).T
            ky = self.mo2ao(mo, kX['Nb'][('y', w)]).T
            kz = self.mo2ao(mo, kX['Nb'][('z', w)]).T

            # create the first order single indexed densiteies #
            Dx = self.transform_dens(kx, D0, S)
            Dy = self.transform_dens(ky, D0, S)
            Dz = self.transform_dens(kz, D0, S)

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

            # Create first order three indexed Densities #

            density_list.append(D_sig_xx.real)
            density_list.append(D_sig_yy.real)
            density_list.append(D_sig_zz.real)
            density_list.append(D_sig_xy.real)
            density_list.append(D_sig_xz.real)
            density_list.append(D_sig_yz.real)

        return density_list

    def get_densities(self, wi, kX, S, D0, mo, nocc, norb):
        """
        Computes the compounded densities needed for the compounded Fock
        matrics F^{σ},F^{λ+τ},F^{σλτ} used for the isotropic cubic response
        function

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
        :param nocc:
            The number of occupied orbitals
        :param norb:
            The number of total orbitals

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

            D_lam_sig_tau_x = self.transform_dens(
                kx_, D_sig_xx, S) + self.transform_dens(
                    ky_, D_sig_xy, S) + self.transform_dens(kz_, D_sig_xz, S)
            D_lam_sig_tau_x += self.transform_dens(
                kx, D_lamtau_xx, S) + self.transform_dens(
                    ky, D_lamtau_xy, S) + self.transform_dens(
                        kz, D_lamtau_xz, S)

            D_lam_sig_tau_y = self.transform_dens(
                kx_, D_sig_xy, S) + self.transform_dens(
                    ky_, D_sig_yy, S) + self.transform_dens(kz_, D_sig_yz, S)
            D_lam_sig_tau_y += self.transform_dens(
                kx, D_lamtau_xy, S) + self.transform_dens(
                    ky, D_lamtau_yy, S) + self.transform_dens(
                        kz, D_lamtau_yz, S)

            D_lam_sig_tau_z = self.transform_dens(
                kx_, D_sig_xz, S) + self.transform_dens(
                    ky_, D_sig_yz, S) + self.transform_dens(kz_, D_sig_zz, S)
            D_lam_sig_tau_z += self.transform_dens(
                kx, D_lamtau_xz, S) + self.transform_dens(
                    ky, D_lamtau_yz, S) + self.transform_dens(
                        kz, D_lamtau_zz, S)

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

    def get_fock_dict(self, wi, kX, density_list, D0, mo, molecule, ao_basis):
        """
        Computes the compounded Fock matrics F^{σ},F^{λ+τ},F^{σλτ} used for the
        isotropic cubic response function

        :param wi:
            A list of the frequencies
        :param kX:
            A dictonary with all the first-order response matrices
        :param D0:
            The SCF density matrix in AO basis
        :param mo:
            A matrix containing the MO coefficents
        :param molecule:
        :param ao_basis:

        :return:
            A dictonary of compounded first-order Fock-matrices

        """
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

        if self.rank == mpi_master():
            fock_num = 30 * len(wi)
            self.print_header(fock_num)
        time_start_fock = time.time()

        # computes both real and imaginary Fock matrices
        ff_AO = self.get_fock_r(mo, density_list, molecule, ao_basis, 1)
        # The zeroth order Fock matrix is computed
        F0_a = self.get_fock_r(mo, D0, molecule, ao_basis, 0)

        time_end_fock = time.time()
        total_time_fock = time_end_fock - time_start_fock
        self.print_time(total_time_fock)

        if self.rank == mpi_master():

            fock_list = ff_AO

            count = 0
            for w in wi:
                f_sig_xx.update({w: fock_list[15 * count]})
                f_sig_yy.update({w: fock_list[15 * count + 1]})
                f_sig_zz.update({w: fock_list[15 * count + 2]})
                f_sig_xy.update({w: fock_list[15 * count + 3]})
                f_sig_xz.update({w: fock_list[15 * count + 4]})
                f_sig_yz.update({w: fock_list[15 * count + 5]})

                f_lamtau_xx.update({w: fock_list[15 * count + 6]})
                f_lamtau_yy.update({w: fock_list[15 * count + 7]})
                f_lamtau_zz.update({w: fock_list[15 * count + 8]})
                f_lamtau_xy.update({w: fock_list[15 * count + 9]})
                f_lamtau_xz.update({w: fock_list[15 * count + 10]})
                f_lamtau_yz.update({w: fock_list[15 * count + 11]})

                f_lam_sig_tau_x.update({w: fock_list[15 * count + 12]})
                f_lam_sig_tau_y.update({w: fock_list[15 * count + 13]})
                f_lam_sig_tau_z.update({w: fock_list[15 * count + 14]})

                count += 1

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
        else:
            Fock = {}

        return Fock

    def get_fock_dict_red(self, wi, kX, density_list, D0, mo, molecule,
                          ao_basis):
        """
        Computes the compounded Fock matrics F^{σ}  used for the reduced
        isotropic cubic response function

        :param wi:
            A list of the frequencies
        :param kX:
            A dictonary with all the first-order response matrices
        :param D0:
            The SCF density matrix in AO basis
        :param mo:
            A matrix containing the MO coefficents
        :param molecule:
        :param ao_basis:

        :return:
            A dictonary of compounded first-order Fock-matrices
        """

        f_sig_xx = {}
        f_sig_yy = {}
        f_sig_zz = {}
        f_sig_xy = {}
        f_sig_xz = {}
        f_sig_yz = {}

        if self.rank == mpi_master():
            fock_num = 6 * len(wi)
            self.print_header(fock_num)
        time_start_fock = time.time()

        ff_AO = self.get_fock_r(mo, density_list, molecule, ao_basis, 2)
        F0_a = self.get_fock_r(mo, D0, molecule, ao_basis, 0)

        time_end_fock = time.time()
        total_time_fock = time_end_fock - time_start_fock
        self.print_time(total_time_fock)

        if self.rank == mpi_master():

            fock_list = ff_AO

            count = 0
            for w in wi:
                f_sig_xx.update({w: fock_list[6 * count]})
                f_sig_yy.update({w: fock_list[6 * count + 1]})
                f_sig_zz.update({w: fock_list[6 * count + 2]})
                f_sig_xy.update({w: fock_list[6 * count + 3]})
                f_sig_xz.update({w: fock_list[6 * count + 4]})
                f_sig_yz.update({w: fock_list[6 * count + 5]})

                count += 1

            Fock = {
                'F0': F0_a,
                'f_sig_xx': f_sig_xx,
                'f_sig_yy': f_sig_yy,
                'f_sig_zz': f_sig_zz,
                'f_sig_xy': f_sig_xy,
                'f_sig_xz': f_sig_xz,
                'f_sig_yz': f_sig_yz
            }
        else:
            Fock = {}

        return Fock

    def get_xy(self, d_a_mo, X, wi, Fock, kX, nocc, norb):
        """
        Computes the compounded gradient vectors N^{σ},N^{λ+τ} used for the
        isotropic cubic response function

        :param d_a_mo:
            The SCF density matrix in MO basis
        :param kX:
            A dictonary with all the first-order response matrices
        :param wi:
            A list of the frequencies
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
                xy_dict.update(
                    {
                        (('N_sig_xx', w), 2 * w):
                            self.anti_sym(
                                -2 * LinearSolver.lrmat2vec(
                                    (3 * self.xi(kx, kx, f_x, f_x, F0) +
                                     self.xi(ky, ky, f_y, f_y, F0) +
                                     self.xi(kz, kz, f_z, f_z, F0) +
                                     0.5 * Fock['f_sig_xx'][w]).T, nocc, norb))
                            - 3 * 2 *
                            self.x2_contract(kx.T, mu_x, d_a_mo, nocc, norb) -
                            2 * self.x2_contract(ky.T, mu_y, d_a_mo, nocc, norb)
                            -
                            2 * self.x2_contract(kz.T, mu_z, d_a_mo, nocc, norb)
                    })
                xy_dict.update(
                    {
                        (('N_sig_yy', w), 2 * w):
                            self.anti_sym(
                                -2 * LinearSolver.lrmat2vec(
                                    (self.xi(kx, kx, f_x, f_x, F0) +
                                     3 * self.xi(ky, ky, f_y, f_y, F0) +
                                     self.xi(kz, kz, f_z, f_z, F0) +
                                     0.5 * Fock['f_sig_yy'][w]).T, nocc, norb))
                            -
                            2 * self.x2_contract(kx.T, mu_x, d_a_mo, nocc, norb)
                            - 3 * 2 *
                            self.x2_contract(ky.T, mu_y, d_a_mo, nocc, norb) -
                            2 * self.x2_contract(kz.T, mu_z, d_a_mo, nocc, norb)
                    })
                xy_dict.update({
                    (('N_sig_zz', w), 2 * w):
                        self.anti_sym(-2 *
                                      LinearSolver.lrmat2vec(
                                          (self.xi(kx, kx, f_x, f_x, F0) +
                                           self.xi(ky, ky, f_y, f_y, F0) +
                                           3 * self.xi(kz, kz, f_z, f_z, F0) +
                                           0.5 * Fock['f_sig_zz'][w]).T, nocc,
                                          norb))
                        - 2 * self.x2_contract(kx.T, mu_x, d_a_mo, nocc, norb) -
                        2 * self.x2_contract(ky.T, mu_y, d_a_mo, nocc, norb) -
                        3 * 2 * self.x2_contract(kz.T, mu_z, d_a_mo, nocc, norb)
                })

                xy_dict.update({
                    (('N_sig_xy', w), 2 * w):
                        self.anti_sym(-2 * LinearSolver.lrmat2vec(
                            (self.xi(ky, kx, f_y, f_x, F0) +
                             self.xi(kx, ky, f_x, f_y, F0) +
                             0.5 * Fock['f_sig_xy'][w]).T, nocc, norb)) -
                        2 * self.x2_contract(ky.T, mu_x, d_a_mo, nocc, norb) -
                        2 * self.x2_contract(kx.T, mu_y, d_a_mo, nocc, norb)
                })
                xy_dict.update({
                    (('N_sig_xz', w), 2 * w):
                        self.anti_sym(-2 * LinearSolver.lrmat2vec(
                            (self.xi(kz, kx, f_z, f_x, F0) +
                             self.xi(kx, kz, f_x, f_z, F0) +
                             0.5 * Fock['f_sig_xz'][w]).T, nocc, norb)) -
                        2 * self.x2_contract(kz.T, mu_x, d_a_mo, nocc, norb) -
                        2 * self.x2_contract(kx.T, mu_z, d_a_mo, nocc, norb)
                })
                xy_dict.update({
                    (('N_sig_yz', w), 2 * w):
                        self.anti_sym(-2 * LinearSolver.lrmat2vec(
                            (self.xi(kz, ky, f_z, f_y, F0) +
                             self.xi(ky, kz, f_y, f_z, F0) +
                             0.5 * Fock['f_sig_yz'][w]).T, nocc, norb)) -
                        2 * self.x2_contract(kz.T, mu_y, d_a_mo, nocc, norb) -
                        2 * self.x2_contract(ky.T, mu_z, d_a_mo, nocc, norb)
                })

                # BC CD λ+τ gradients #
                xy_dict.update(
                    {
                        (('N_lamtau_xx', w), 0):
                            self.anti_sym(
                                -2 * LinearSolver.lrmat2vec(
                                    (3 * 2 * self.xi(kx_, kx, f_x_, f_x, F0) +
                                     2 * self.xi(ky_, ky, f_y_, f_y, F0) +
                                     2 * self.xi(kz_, kz, f_z_, f_z, F0) +
                                     0.5 * Fock['f_lamtau_xx'][w]).T, nocc,
                                    norb))
                            - 3 * 2 *
                            self.x2_contract(kx_.T, mu_x, d_a_mo, nocc, norb) -
                            2 *
                            self.x2_contract(ky_.T, mu_y, d_a_mo, nocc, norb) -
                            2 *
                            self.x2_contract(kz_.T, mu_z, d_a_mo, nocc, norb) -
                            3 * 2 *
                            self.x2_contract(kx.T, mu_x, d_a_mo, nocc, norb) -
                            2 * self.x2_contract(ky.T, mu_y, d_a_mo, nocc, norb)
                            -
                            2 * self.x2_contract(kz.T, mu_z, d_a_mo, nocc, norb)
                    })

                xy_dict.update(
                    {
                        (('N_lamtau_yy', w), 0):
                            self.anti_sym(
                                -2 * LinearSolver.lrmat2vec(
                                    (2 * self.xi(kx_, kx, f_x_, f_x, F0) +
                                     3 * 2 * self.xi(ky_, ky, f_y_, f_y, F0) +
                                     2 * self.xi(kz_, kz, f_z_, f_z, F0) +
                                     0.5 * Fock['f_lamtau_yy'][w]).T, nocc,
                                    norb))
                            - 2 * self.x2_contract(kx_.T, mu_x, d_a_mo,
                                                   nocc, norb) - 3 * 2 *
                            self.x2_contract(ky_.T, mu_y, d_a_mo, nocc, norb) -
                            2 *
                            self.x2_contract(kz_.T, mu_z, d_a_mo, nocc, norb) -
                            2 * self.x2_contract(kx.T, mu_x, d_a_mo, nocc, norb)
                            - 3 * 2 *
                            self.x2_contract(ky.T, mu_y, d_a_mo, nocc, norb) -
                            2 * self.x2_contract(kz.T, mu_z, d_a_mo, nocc, norb)
                    })

                xy_dict.update(
                    {
                        (('N_lamtau_zz', w), 0):
                            self.anti_sym(
                                -2 * LinearSolver.lrmat2vec(
                                    (2 * self.xi(kx_, kx, f_x_, f_x, F0) +
                                     2 * self.xi(ky_, ky, f_y_, f_y, F0) +
                                     3 * 2 * self.xi(kz_, kz, f_z_, f_z, F0) +
                                     0.5 * Fock['f_lamtau_zz'][w]).T, nocc,
                                    norb))
                            - 2 * self.x2_contract(kx_.T, mu_x, d_a_mo,
                                                   nocc, norb) - 2 *
                            self.x2_contract(ky_.T, mu_y, d_a_mo, nocc, norb) -
                            3 * 2 *
                            self.x2_contract(kz_.T, mu_z, d_a_mo, nocc, norb) -
                            2 * self.x2_contract(kx.T, mu_x, d_a_mo, nocc, norb)
                            -
                            2 * self.x2_contract(ky.T, mu_y, d_a_mo, nocc, norb)
                            - 3 * 2 *
                            self.x2_contract(kz.T, mu_z, d_a_mo, nocc, norb)
                    })

                xy_dict.update({
                    (('N_lamtau_xy', w), 0):
                        self.anti_sym(-2 * LinearSolver.lrmat2vec(
                            (2 * self.xi(ky_, kx, f_y_, f_x, F0) +
                             2 * self.xi(kx_, ky, f_x_, f_y, F0) +
                             0.5 * Fock['f_lamtau_xy'][w]).T, nocc, norb)) -
                        2 * self.x2_contract(ky.T, mu_x, d_a_mo, nocc, norb) -
                        2 * self.x2_contract(kx.T, mu_y, d_a_mo, nocc, norb) -
                        2 * self.x2_contract(ky_.T, mu_x, d_a_mo, nocc, norb) -
                        2 * self.x2_contract(kx_.T, mu_y, d_a_mo, nocc, norb)
                })
                xy_dict.update({
                    (('N_lamtau_xz', w), 0):
                        self.anti_sym(-2 * LinearSolver.lrmat2vec(
                            (2 * self.xi(kz_, kx, f_z_, f_x, F0) +
                             2 * self.xi(kx_, kz, f_x_, f_z, F0) +
                             0.5 * Fock['f_lamtau_xz'][w]).T, nocc, norb)) -
                        2 * self.x2_contract(kz.T, mu_x, d_a_mo, nocc, norb) -
                        2 * self.x2_contract(kx.T, mu_z, d_a_mo, nocc, norb) -
                        2 * self.x2_contract(kz_.T, mu_x, d_a_mo, nocc, norb) -
                        2 * self.x2_contract(kx_.T, mu_z, d_a_mo, nocc, norb)
                })
                xy_dict.update({
                    (('N_lamtau_yz', w), 0):
                        self.anti_sym(-2 * LinearSolver.lrmat2vec(
                            (2 * self.xi(kz_, ky, f_z_, f_y, F0) +
                             2 * self.xi(ky_, kz, f_y_, f_z, F0) +
                             0.5 * Fock['f_lamtau_yz'][w]).T, nocc, norb)) -
                        2 * self.x2_contract(kz.T, mu_y, d_a_mo, nocc, norb) -
                        2 * self.x2_contract(ky.T, mu_z, d_a_mo, nocc, norb) -
                        2 * self.x2_contract(kz_.T, mu_y, d_a_mo, nocc, norb) -
                        2 * self.x2_contract(ky_.T, mu_z, d_a_mo, nocc, norb)
                })
        else:
            xy_dict = {}
        return xy_dict

    def get_xy_red(self, d_a_mo, X, wi, Fock, kX, nocc, norb):
        """
        Computes the compounded gradient vectors N^{σ}  used for the reduced
        isotropic cubic response function

        :param d_a_mo:
            The SCF density matrix in MO basis
        :param kX:
            A dictonary with all the first-order response matrices
        :param wi:
            A list of the frequencies
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

                # REAL PART #
                F0 = Fock['F0']
                f_x = Fock['Fb'][('x', w)]
                f_y = Fock['Fb'][('y', w)]
                f_z = Fock['Fb'][('z', w)]

                # BD σ gradients #
                xy_dict.update(
                    {
                        (('N_sig_xx', w), 2 * w):
                            self.anti_sym(
                                -2 * LinearSolver.lrmat2vec(
                                    (3 * self.xi(kx, kx, f_x, f_x, F0) +
                                     self.xi(ky, ky, f_y, f_y, F0) +
                                     self.xi(kz, kz, f_z, f_z, F0) +
                                     0.5 * Fock['f_sig_xx'][w]).T, nocc, norb))
                            - 3 * 2 *
                            self.x2_contract(kx.T, mu_x, d_a_mo, nocc, norb) -
                            2 * self.x2_contract(ky.T, mu_y, d_a_mo, nocc, norb)
                            -
                            2 * self.x2_contract(kz.T, mu_z, d_a_mo, nocc, norb)
                    })
                xy_dict.update(
                    {
                        (('N_sig_yy', w), 2 * w):
                            self.anti_sym(
                                -2 * LinearSolver.lrmat2vec(
                                    (self.xi(kx, kx, f_x, f_x, F0) +
                                     3 * self.xi(ky, ky, f_y, f_y, F0) +
                                     self.xi(kz, kz, f_z, f_z, F0) +
                                     0.5 * Fock['f_sig_yy'][w]).T, nocc, norb))
                            -
                            2 * self.x2_contract(kx.T, mu_x, d_a_mo, nocc, norb)
                            - 3 * 2 *
                            self.x2_contract(ky.T, mu_y, d_a_mo, nocc, norb) -
                            2 * self.x2_contract(kz.T, mu_z, d_a_mo, nocc, norb)
                    })
                xy_dict.update({
                    (('N_sig_zz', w), 2 * w):
                        self.anti_sym(-2 *
                                      LinearSolver.lrmat2vec(
                                          (self.xi(kx, kx, f_x, f_x, F0) +
                                           self.xi(ky, ky, f_y, f_y, F0) +
                                           3 * self.xi(kz, kz, f_z, f_z, F0) +
                                           0.5 * Fock['f_sig_zz'][w]).T, nocc,
                                          norb))
                        - 2 * self.x2_contract(kx.T, mu_x, d_a_mo, nocc, norb) -
                        2 * self.x2_contract(ky.T, mu_y, d_a_mo, nocc, norb) -
                        3 * 2 * self.x2_contract(kz.T, mu_z, d_a_mo, nocc, norb)
                })

                xy_dict.update({
                    (('N_sig_xy', w), 2 * w):
                        self.anti_sym(-2 * LinearSolver.lrmat2vec(
                            (self.xi(ky, kx, f_y, f_x, F0) +
                             self.xi(kx, ky, f_x, f_y, F0) +
                             0.5 * Fock['f_sig_xy'][w]).T, nocc, norb)) -
                        2 * self.x2_contract(ky.T, mu_x, d_a_mo, nocc, norb) -
                        2 * self.x2_contract(kx.T, mu_y, d_a_mo, nocc, norb)
                })
                xy_dict.update({
                    (('N_sig_xz', w), 2 * w):
                        self.anti_sym(-2 * LinearSolver.lrmat2vec(
                            (self.xi(kz, kx, f_z, f_x, F0) +
                             self.xi(kx, kz, f_x, f_z, F0) +
                             0.5 * Fock['f_sig_xz'][w]).T, nocc, norb)) -
                        2 * self.x2_contract(kz.T, mu_x, d_a_mo, nocc, norb) -
                        2 * self.x2_contract(kx.T, mu_z, d_a_mo, nocc, norb)
                })
                xy_dict.update({
                    (('N_sig_yz', w), 2 * w):
                        self.anti_sym(-2 * LinearSolver.lrmat2vec(
                            (self.xi(kz, ky, f_z, f_y, F0) +
                             self.xi(ky, kz, f_y, f_z, F0) +
                             0.5 * Fock['f_sig_yz'][w]).T, nocc, norb)) -
                        2 * self.x2_contract(kz.T, mu_y, d_a_mo, nocc, norb) -
                        2 * self.x2_contract(ky.T, mu_z, d_a_mo, nocc, norb)
                })

        else:
            xy_dict = {}
        return xy_dict

    def other_red(self, iso, wi, track, n_x, n_xy, X, kX, kXY, da, mol_orbs,
                  task):
        """
        Computes the terms involving X[2],A[2] in the reduced isotropic cubic response function

        : param iso:
            A bolean value that states if its the isotrpoic gamma or a user defined case
        : param wi:
            A list containing all the frequencies
        : param keys:
            A dictonray or lists that are used to tell the program which elements are present
        : param track:
            A list that contains information about what γ components that are to be computed and which freqs
        : param n_x:
            A dictonary containing all the single-index response vectors
        : param n_xy:
            A dictonary containing all the two-index response vectors
        : param X:
            A dictonray with all the property integral matricies
        : param kX:
            A dictonary with all the respone matricies
        : param kXY:
            A dictonary containing all the two-index response matricies
        : param da:
            The SCF density matrix in MO bassi

        :return:
            A dictonary of final X[2],A[2] contraction values
        """

        na_x2_nyz_dict = {}
        nx_a2_nyz_dict = {}

        for i in range(len(wi)):
            wb = float(track[i * int(len(track) / len(wi))].split(",")[1])
            wa = float(track[i * int(len(track) / len(wi))].split(",")[1])
            wc = float(track[i * int(len(track) / len(wi))].split(",")[2])
            w = float(track[i * int(len(track) / len(wi))].split(",")[1])
            wd = float(track[i * int(len(track) / len(wi))].split(",")[3])

            wbd = wb + wd

            na_x2_nyz = 0
            nx_a2_nyz = 0

            # BD

            kbd = kXY[(('N_sig_xx', w), wbd)]
            Nbd = n_xy[(('N_sig_xx', w), wbd)]
            Na = n_x['Na'][('x', wa)]
            Nc = n_x['Nc'][('x', wc)]
            kc = kX['Nc'][('x', wc)]
            A = X['x']
            C = X['x']

            na_x2_nyz += Na.T @ self.x2_contract(kbd, C, da, mol_orbs, task)

            nx_a2_nyz += self.a2_contract(kc, A, da, mol_orbs, task) @ Nbd
            nx_a2_nyz += self.a2_contract(kbd, A, da, mol_orbs, task) @ Nc

            kbd = kXY[(('N_sig_xy', w), wbd)]
            Nbd = n_xy[(('N_sig_xy', w), wbd)]
            Na = n_x['Na'][('x', wa)]
            Nc = n_x['Nc'][('y', wc)]
            kc = kX['Nc'][('y', wc)]
            A = X['x']
            C = X['y']

            na_x2_nyz += Na.T @ self.x2_contract(kbd, C, da, mol_orbs, task)

            nx_a2_nyz += self.a2_contract(kc, A, da, mol_orbs, task) @ Nbd
            nx_a2_nyz += self.a2_contract(kbd, A, da, mol_orbs, task) @ Nc

            kbd = kXY[(('N_sig_xz', w), wbd)]
            Nbd = n_xy[(('N_sig_xz', w), wbd)]
            Nc = n_x['Nc'][('z', wc)]
            kc = kX['Nc'][('z', wc)]
            Na = n_x['Na'][('x', wa)]
            A = X['x']
            C = X['z']

            na_x2_nyz += Na.T @ self.x2_contract(kbd, C, da, mol_orbs, task)

            nx_a2_nyz += self.a2_contract(kc, A, da, mol_orbs, task) @ Nbd
            nx_a2_nyz += self.a2_contract(kbd, A, da, mol_orbs, task) @ Nc

            # y

            kbd = kXY[(('N_sig_xy', w), wbd)]
            Nbd = n_xy[(('N_sig_xy', w), wbd)]
            Na = n_x['Na'][('y', wa)]
            Nc = n_x['Nc'][('x', wc)]
            kc = kX['Nc'][('x', wc)]
            A = X['y']
            C = X['x']

            na_x2_nyz += Na.T @ self.x2_contract(kbd, C, da, mol_orbs, task)

            nx_a2_nyz += self.a2_contract(kc, A, da, mol_orbs, task) @ Nbd
            nx_a2_nyz += self.a2_contract(kbd, A, da, mol_orbs, task) @ Nc

            kbd = kXY[(('N_sig_yy', w), wbd)]
            Nbd = n_xy[(('N_sig_yy', w), wbd)]
            Nc = n_x['Nc'][('y', wc)]
            Na = n_x['Na'][('y', wa)]
            kc = kX['Nc'][('y', wc)]
            A = X['y']
            C = X['y']

            na_x2_nyz += Na.T @ self.x2_contract(kbd, C, da, mol_orbs, task)

            nx_a2_nyz += self.a2_contract(kc, A, da, mol_orbs, task) @ Nbd
            nx_a2_nyz += self.a2_contract(kbd, A, da, mol_orbs, task) @ Nc

            kbd = kXY[(('N_sig_yz', w), wbd)]
            Nbd = n_xy[(('N_sig_yz', w), wbd)]
            Nc = n_x['Nc'][('z', wc)]
            Na = n_x['Na'][('y', wa)]
            kc = kX['Nc'][('z', wc)]
            A = X['y']
            C = X['z']

            na_x2_nyz += Na.T @ self.x2_contract(kbd, C, da, mol_orbs, task)

            nx_a2_nyz += self.a2_contract(kc, A, da, mol_orbs, task) @ Nbd
            nx_a2_nyz += self.a2_contract(kbd, A, da, mol_orbs, task) @ Nc

            # z

            kbd = kXY[(('N_sig_xz', w), wbd)]
            Nbd = n_xy[(('N_sig_xz', w), wbd)]
            Nc = n_x['Nc'][('x', wc)]
            kc = kX['Nc'][('x', wc)]
            Na = n_x['Na'][('z', wa)]
            A = X['z']
            C = X['x']

            na_x2_nyz += Na.T @ self.x2_contract(kbd, C, da, mol_orbs, task)

            nx_a2_nyz += self.a2_contract(kc, A, da, mol_orbs, task) @ Nbd
            nx_a2_nyz += self.a2_contract(kbd, A, da, mol_orbs, task) @ Nc

            kbd = kXY[(('N_sig_yz', w), wbd)]
            Nbd = n_xy[(('N_sig_yz', w), wbd)]
            Nc = n_x['Nc'][('y', wc)]
            kc = kX['Nc'][('y', wc)]
            Na = n_x['Na'][('z', wa)]
            A = X['z']
            C = X['y']

            na_x2_nyz += Na.T @ self.x2_contract(kbd, C, da, mol_orbs, task)

            nx_a2_nyz += self.a2_contract(kc, A, da, mol_orbs, task) @ Nbd
            nx_a2_nyz += self.a2_contract(kbd, A, da, mol_orbs, task) @ Nc

            kbd = kXY[(('N_sig_zz', w), wbd)]
            Nbd = n_xy[(('N_sig_zz', w), wbd)]
            Nc = n_x['Nc'][('z', wc)]
            Na = n_x['Na'][('z', wa)]
            kc = kX['Nc'][('z', wc)]
            A = X['z']
            C = X['z']

            na_x2_nyz += Na.T @ self.x2_contract(kbd, C, da, mol_orbs, task)

            nx_a2_nyz += self.a2_contract(kc, A, da, mol_orbs, task) @ Nbd
            nx_a2_nyz += self.a2_contract(kbd, A, da, mol_orbs, task) @ Nc

            na_x2_nyz_dict.update({(w, -w, w): na_x2_nyz})
            nx_a2_nyz_dict.update({(w, -w, w): nx_a2_nyz})

        return na_x2_nyz_dict, nx_a2_nyz_dict

    def other(self, iso, wi, track, n_x, n_xy, X, kX, kXY, da, mol_orbs, task):
        """
        Computes the terms involving X[3],A[3],X[2],A[2] in the isotropic cubic
        response function

        : param iso:
            A bolean value that states if its the isotrpoic gamma or a user defined case
        : param wi:
            A list containing all the frequencies
        : param keys:
            A dictonray or lists that are used to tell the program which elements are present
        : param track:
            A list that contains information about what γ components that are to be computed and which freqs
        : param n_x:
            A dictonary containing all the single-index response vectors
        : param n_xy:
            A dictonary containing all the two-index response vectors
        : param X:
            A dictonray with all the property integral matricies
        : param kX:
            A dictonary with all the respone matricies
        : param kXY:
            A dictonary containing all the two-index response matricies
        : param da:
            The SCF density matrix in MO bassi

        :return:
            A dictonary of final X[2],A[2] contraction values
        """

        na_a3_nx_ny_dict = {}
        na_x3_ny_nz_dict = {}
        na_x2_nyz_dict = {}
        nx_a2_nyz_dict = {}

        count = 0
        for j in range(len(wi)):
            na_x3_ny_nz = 0
            na_a3_nx_ny = 0
            for i in range(int(len(track) * (1 / len(wi)))):
                i += count

                w1 = float(track[i].split(",")[1])
                w2 = float(track[i].split(",")[2])
                w3 = float(track[i].split(",")[3])

                Na = n_x['Na'][(track[i][0], w1)]
                Nb = n_x['Nb'][(track[i][1], w1)]

                #
                Nc = n_x['Nc'][(track[i][2], w2)]
                Nd = n_x['Nd'][(track[i][3], w3)]

                kb = kX['Nb'][(track[i][1], w1)]
                kc = kX['Nc'][(track[i][2], w2)]
                kd = kX['Nd'][(track[i][3], w3)]

                A = X[track[i][0]]
                B = X[track[i][1]]
                C = X[track[i][2]]
                D = X[track[i][3]]

                # Na X[3]NyNz

                na_x3_ny_nz += -Na.T @ self.x3_contract(kc, kd, B, da, mol_orbs,
                                                        task)
                na_x3_ny_nz += -Na.T @ self.x3_contract(kd, kc, B, da, mol_orbs,
                                                        task)
                na_x3_ny_nz += -Na.T @ self.x3_contract(kd, kb, C, da, mol_orbs,
                                                        task)
                na_x3_ny_nz += -Na.T @ self.x3_contract(kb, kd, C, da, mol_orbs,
                                                        task)
                na_x3_ny_nz += -Na.T @ self.x3_contract(kb, kc, D, da, mol_orbs,
                                                        task)
                na_x3_ny_nz += -Na.T @ self.x3_contract(kc, kb, D, da, mol_orbs,
                                                        task)

                # NaA[3]n_xNy

                na_a3_nx_ny += self.a3_contract(kb, kc, A, da, mol_orbs,
                                                task) @ Nd
                na_a3_nx_ny += self.a3_contract(kb, kd, A, da, mol_orbs,
                                                task) @ Nc
                na_a3_nx_ny += self.a3_contract(kc, kb, A, da, mol_orbs,
                                                task) @ Nd
                na_a3_nx_ny += self.a3_contract(kc, kd, A, da, mol_orbs,
                                                task) @ Nb
                na_a3_nx_ny += self.a3_contract(kd, kb, A, da, mol_orbs,
                                                task) @ Nc
                na_a3_nx_ny += self.a3_contract(kd, kc, A, da, mol_orbs,
                                                task) @ Nb

            count += int(len(track) / len(wi))

            na_a3_nx_ny_dict.update({(wi[j], -wi[j], wi[j]): na_a3_nx_ny})
            na_x3_ny_nz_dict.update({(wi[j], -wi[j], wi[j]): na_x3_ny_nz})

        for i in range(len(wi)):
            wb = float(track[i * int(len(track) / len(wi))].split(",")[1])
            wa = float(track[i * int(len(track) / len(wi))].split(",")[1])
            wc = float(track[i * int(len(track) / len(wi))].split(",")[2])
            w = float(track[i * int(len(track) / len(wi))].split(",")[1])
            wd = float(track[i * int(len(track) / len(wi))].split(",")[3])

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
            na_x2_nyz += Na.T @ self.x2_contract(kcd, B, da, mol_orbs, task)

            nx_a2_nyz += self.a2_contract(kb, A, da, mol_orbs, task) @ Ncd
            nx_a2_nyz += self.a2_contract(kcd, A, da, mol_orbs, task) @ Nb

            kcd = kXY[(('N_lamtau_xy', w), wcd)]
            Ncd = n_xy[(('N_lamtau_xy', w), wcd)]

            Nb = n_x['Nb'][('y', w)]
            kb = kX['Nb'][('y', w)]
            Na = n_x['Na'][('x', wa)]
            A = X['x']
            B = X['y']

            na_x2_nyz += Na.T @ self.x2_contract(kcd, B, da, mol_orbs, task)

            nx_a2_nyz += self.a2_contract(kb, A, da, mol_orbs, task) @ Ncd
            nx_a2_nyz += self.a2_contract(kcd, A, da, mol_orbs, task) @ Nb

            kcd = kXY[(('N_lamtau_xz', w), wcd)]
            Ncd = n_xy[(('N_lamtau_xz', w), wcd)]
            Nb = n_x['Nb'][('z', w)]
            kb = kX['Nb'][('z', w)]
            Na = n_x['Na'][('x', wa)]
            A = X['x']
            B = X['z']

            na_x2_nyz += Na.T @ self.x2_contract(kcd, B, da, mol_orbs, task)

            nx_a2_nyz += self.a2_contract(kb, A, da, mol_orbs, task) @ Ncd
            nx_a2_nyz += self.a2_contract(kcd, A, da, mol_orbs, task) @ Nb

            # y

            kcd = kXY[(('N_lamtau_xy', w), wcd)]
            Ncd = n_xy[(('N_lamtau_xy', w), wcd)]
            Nb = n_x['Nb'][('x', w)]
            kb = kX['Nb'][('x', w)]
            A = X['y']
            B = X['x']
            Na = n_x['Na'][('y', wa)]

            na_x2_nyz += Na.T @ self.x2_contract(kcd, B, da, mol_orbs, task)

            nx_a2_nyz += self.a2_contract(kb, A, da, mol_orbs, task) @ Ncd
            nx_a2_nyz += self.a2_contract(kcd, A, da, mol_orbs, task) @ Nb

            kcd = kXY[(('N_lamtau_yy', w), wcd)]
            Ncd = n_xy[(('N_lamtau_yy', w), wcd)]
            Na = n_x['Na'][('y', wa)]
            Nb = n_x['Nb'][('y', w)]
            kb = kX['Nb'][('y', w)]
            A = X['y']
            B = X['y']

            na_x2_nyz += Na.T @ self.x2_contract(kcd, B, da, mol_orbs, task)

            nx_a2_nyz += self.a2_contract(kb, A, da, mol_orbs, task) @ Ncd
            nx_a2_nyz += self.a2_contract(kcd, A, da, mol_orbs, task) @ Nb

            kcd = kXY[(('N_lamtau_yz', w), wcd)]
            Ncd = n_xy[(('N_lamtau_yz', w), wcd)]
            Nb = n_x['Nb'][('z', w)]
            kb = kX['Nb'][('z', w)]

            A = X['y']
            B = X['z']
            Na = n_x['Na'][('y', wa)]

            na_x2_nyz += Na.T @ self.x2_contract(kcd, B, da, mol_orbs, task)

            nx_a2_nyz += self.a2_contract(kb, A, da, mol_orbs, task) @ Ncd
            nx_a2_nyz += self.a2_contract(kcd, A, da, mol_orbs, task) @ Nb

            # z

            kcd = kXY[(('N_lamtau_xz', w), wcd)]
            Ncd = n_xy[(('N_lamtau_xz', w), wcd)]
            Nb = n_x['Nb'][('x', w)]
            kb = kX['Nb'][('x', w)]
            Na = n_x['Na'][('z', wa)]
            A = X['z']
            B = X['x']

            na_x2_nyz += Na.T @ self.x2_contract(kcd, B, da, mol_orbs, task)

            nx_a2_nyz += self.a2_contract(kb, A, da, mol_orbs, task) @ Ncd
            nx_a2_nyz += self.a2_contract(kcd, A, da, mol_orbs, task) @ Nb

            kcd = kXY[(('N_lamtau_yz', w), wcd)]
            Ncd = n_xy[(('N_lamtau_yz', w), wcd)]
            Na = n_x['Na'][('z', wa)]
            Nb = n_x['Nb'][('y', w)]
            kb = kX['Nb'][('y', w)]
            A = X['z']
            B = X['y']

            na_x2_nyz += Na.T @ self.x2_contract(kcd, B, da, mol_orbs, task)

            nx_a2_nyz += self.a2_contract(kb, A, da, mol_orbs, task) @ Ncd
            nx_a2_nyz += self.a2_contract(kcd, A, da, mol_orbs, task) @ Nb

            kcd = kXY[(('N_lamtau_zz', w), wcd)]
            Ncd = n_xy[(('N_lamtau_zz', w), wcd)]
            Nb = n_x['Nb'][('z', w)]
            Na = n_x['Na'][('z', wa)]
            kb = kX['Nb'][('z', w)]
            A = X['z']
            B = X['z']

            na_x2_nyz += Na.T @ self.x2_contract(kcd, B, da, mol_orbs, task)

            nx_a2_nyz += self.a2_contract(kb, A, da, mol_orbs, task) @ Ncd
            nx_a2_nyz += self.a2_contract(kcd, A, da, mol_orbs, task) @ Nb

            # BD

            kbd = kXY[(('N_sig_xx', w), wbd)]
            Nbd = n_xy[(('N_sig_xx', w), wbd)]
            Na = n_x['Na'][('x', wa)]
            Nc = n_x['Nc'][('x', wc)]
            kc = kX['Nc'][('x', wc)]
            A = X['x']
            C = X['x']

            na_x2_nyz += Na.T @ self.x2_contract(kbd, C, da, mol_orbs, task)

            nx_a2_nyz += self.a2_contract(kc, A, da, mol_orbs, task) @ Nbd
            nx_a2_nyz += self.a2_contract(kbd, A, da, mol_orbs, task) @ Nc

            kbd = kXY[(('N_sig_xy', w), wbd)]
            Nbd = n_xy[(('N_sig_xy', w), wbd)]
            Na = n_x['Na'][('x', wa)]
            Nc = n_x['Nc'][('y', wc)]
            kc = kX['Nc'][('y', wc)]
            A = X['x']
            C = X['y']

            na_x2_nyz += Na.T @ self.x2_contract(kbd, C, da, mol_orbs, task)

            nx_a2_nyz += self.a2_contract(kc, A, da, mol_orbs, task) @ Nbd
            nx_a2_nyz += self.a2_contract(kbd, A, da, mol_orbs, task) @ Nc

            kbd = kXY[(('N_sig_xz', w), wbd)]
            Nbd = n_xy[(('N_sig_xz', w), wbd)]
            Nc = n_x['Nc'][('z', wc)]
            kc = kX['Nc'][('z', wc)]
            Na = n_x['Na'][('x', wa)]
            A = X['x']
            C = X['z']

            na_x2_nyz += Na.T @ self.x2_contract(kbd, C, da, mol_orbs, task)

            nx_a2_nyz += self.a2_contract(kc, A, da, mol_orbs, task) @ Nbd
            nx_a2_nyz += self.a2_contract(kbd, A, da, mol_orbs, task) @ Nc

            # y

            kbd = kXY[(('N_sig_xy', w), wbd)]
            Nbd = n_xy[(('N_sig_xy', w), wbd)]
            Na = n_x['Na'][('y', wa)]
            Nc = n_x['Nc'][('x', wc)]
            kc = kX['Nc'][('x', wc)]
            A = X['y']
            C = X['x']

            na_x2_nyz += Na.T @ self.x2_contract(kbd, C, da, mol_orbs, task)

            nx_a2_nyz += self.a2_contract(kc, A, da, mol_orbs, task) @ Nbd
            nx_a2_nyz += self.a2_contract(kbd, A, da, mol_orbs, task) @ Nc

            kbd = kXY[(('N_sig_yy', w), wbd)]
            Nbd = n_xy[(('N_sig_yy', w), wbd)]
            Nc = n_x['Nc'][('y', wc)]
            Na = n_x['Na'][('y', wa)]
            kc = kX['Nc'][('y', wc)]
            A = X['y']
            C = X['y']

            na_x2_nyz += Na.T @ self.x2_contract(kbd, C, da, mol_orbs, task)

            nx_a2_nyz += self.a2_contract(kc, A, da, mol_orbs, task) @ Nbd
            nx_a2_nyz += self.a2_contract(kbd, A, da, mol_orbs, task) @ Nc

            kbd = kXY[(('N_sig_yz', w), wbd)]
            Nbd = n_xy[(('N_sig_yz', w), wbd)]
            Nc = n_x['Nc'][('z', wc)]
            Na = n_x['Na'][('y', wa)]
            kc = kX['Nc'][('z', wc)]
            A = X['y']
            C = X['z']

            na_x2_nyz += Na.T @ self.x2_contract(kbd, C, da, mol_orbs, task)

            nx_a2_nyz += self.a2_contract(kc, A, da, mol_orbs, task) @ Nbd
            nx_a2_nyz += self.a2_contract(kbd, A, da, mol_orbs, task) @ Nc

            # z

            kbd = kXY[(('N_sig_xz', w), wbd)]
            Nbd = n_xy[(('N_sig_xz', w), wbd)]
            Nc = n_x['Nc'][('x', wc)]
            kc = kX['Nc'][('x', wc)]
            Na = n_x['Na'][('z', wa)]
            A = X['z']
            C = X['x']

            na_x2_nyz += Na.T @ self.x2_contract(kbd, C, da, mol_orbs, task)

            nx_a2_nyz += self.a2_contract(kc, A, da, mol_orbs, task) @ Nbd
            nx_a2_nyz += self.a2_contract(kbd, A, da, mol_orbs, task) @ Nc

            kbd = kXY[(('N_sig_yz', w), wbd)]
            Nbd = n_xy[(('N_sig_yz', w), wbd)]
            Nc = n_x['Nc'][('y', wc)]
            kc = kX['Nc'][('y', wc)]
            Na = n_x['Na'][('z', wa)]
            A = X['z']
            C = X['y']

            na_x2_nyz += Na.T @ self.x2_contract(kbd, C, da, mol_orbs, task)

            nx_a2_nyz += self.a2_contract(kc, A, da, mol_orbs, task) @ Nbd
            nx_a2_nyz += self.a2_contract(kbd, A, da, mol_orbs, task) @ Nc

            kbd = kXY[(('N_sig_zz', w), wbd)]
            Nbd = n_xy[(('N_sig_zz', w), wbd)]
            Nc = n_x['Nc'][('z', wc)]
            Na = n_x['Na'][('z', wa)]
            kc = kX['Nc'][('z', wc)]
            A = X['z']
            C = X['z']

            na_x2_nyz += Na.T @ self.x2_contract(kbd, C, da, mol_orbs, task)

            nx_a2_nyz += self.a2_contract(kc, A, da, mol_orbs, task) @ Nbd
            nx_a2_nyz += self.a2_contract(kbd, A, da, mol_orbs, task) @ Nc

            na_x2_nyz_dict.update({(w, -w, w): na_x2_nyz})
            nx_a2_nyz_dict.update({(w, -w, w): nx_a2_nyz})

        return na_x3_ny_nz_dict, na_a3_nx_ny_dict, na_x2_nyz_dict, nx_a2_nyz_dict

    def get_t4(self, wi, e4_dict, n_x, Kx, track, da, mol_orbs, task):
        """
        Computes the contraction of the E[4] tensor with that of the S[4] and
        R[4] tensors to return the contraction of T[4] as a dictonary of
        vectors.
        T[4]n_xNyNz = (E^[4]-ω_1S^[4]-ω_1S^[4]-ω_3S^[4]-γiR^[4])

        : param wi:
            A list of all the freqs
        : param e4_dict:
            A dictonary of all the E[4] contraction
        : param n_x:
            A dictonary with all the single index response vectors
        : param kX:
            A dictonray containng all the response matricies
        : param track:
            A list containg information about all the γ components that are to be computed
        : param da:
            The SCF density matrix in MO basis

        :return:
            A dictonary of final T[4] contraction values
        """
        count = 0
        T4term = {}
        S4 = self.S4_dict(wi, Kx, track, da, mol_orbs, task)

        if self.damping > 0:
            R4term = self.get_r4(wi, Kx, n_x, track, da, mol_orbs, task)

        for i in range(len(wi)):

            w = float(track[i * int(len(track) / len(wi))].split(",")[1])
            ww = float(track[i * int(len(track) / len(wi))].split(",")[1])
            t4term = n_x['Na'][('x', w)] @ (e4_dict['f_iso_x'][ww] - S4[
                ('x', ww)]) + n_x['Na'][('y', w)] @ (
                    e4_dict['f_iso_y'][ww] - S4[('y', ww)]) + n_x['Na'][
                        ('z', w)] @ (e4_dict['f_iso_z'][ww] - S4[('z', ww)])

            if self.damping > 0:
                t4term += R4term[('x', ww)] + R4term[('y', ww)] + R4term[('z',
                                                                          ww)]

            T4term.update({(ww, -ww, ww): t4term})
            count += int(len(track) / len(wi))
        return T4term

    def S4_dict(self, wi, kX, track, D0, nocc, norb):
        """
        Computes the S4 contractions

        : param wi:
            A list of all the freqs
        : param kX:
            A dict with all the response matricies in MO basis
        : param track:
            A list containing information about all the components that are to be computed
        : param D0:
            The SCF density in MO basis
        : param nocc:
            The number of occupied obritals
        : param norb:
            The number of total orbitals

         :return:
                 A dictonary of final S[4] contraction values

        """

        S4terms = {}

        count = 0
        for j in range(len(wi)):
            w = float(track[j * int(len(track) / len(wi))].split(",")[1])
            w1 = float(track[j * int(len(track) / len(wi))].split(",")[1])
            w2 = float(track[j * int(len(track) / len(wi))].split(",")[2])
            w3 = float(track[j * int(len(track) / len(wi))].split(",")[3])

            S4_term_x = 0
            S4_term_y = 0
            S4_term_z = 0

            for i in range(int(len(track) * (1 / len(wi)))):
                i = i + count

                kB = kX['Nb'][(track[i][1], w1)]
                kC = kX['Nc'][(track[i][2], w2)]
                kD = kX['Nd'][(track[i][3], w3)]

                if track[i][0] in 'x':

                    S4_term_x += w1 * self.s4(kB, kC, kD, D0, nocc, norb)

                    S4_term_x += w2 * self.s4(kC, kB, kD, D0, nocc, norb)

                    S4_term_x += w3 * self.s4(kD, kB, kC, D0, nocc, norb)

                if track[i][0] in 'y':

                    S4_term_y += w1 * self.s4(kB, kC, kD, D0, nocc, norb)

                    S4_term_y += w2 * self.s4(kC, kB, kD, D0, nocc, norb)

                    S4_term_y += w3 * self.s4(kD, kB, kC, D0, nocc, norb)

                if track[i][0] in 'z':

                    S4_term_z += w1 * self.s4(kB, kC, kD, D0, nocc, norb)

                    S4_term_z += w2 * self.s4(kC, kB, kD, D0, nocc, norb)

                    S4_term_z += w3 * self.s4(kD, kB, kC, D0, nocc, norb)

            count += int(len(track) / len(wi))

            S4terms.update({('x', w): -S4_term_x})
            S4terms.update({('y', w): -S4_term_y})
            S4terms.update({('z', w): -S4_term_z})

        return S4terms

    def s4_for_r4(self, k1, k2, k3, D, nocc, norb):
        """
        Returns the contraction of S[4] for the contraction of R[4]
        :param k1,k2,k3:
            A response matrix
        :param D:
            A density matrix
        :param nocc:
            The number of occupied orbtials
        :param norb:
            The number of total orbitals
        """
        S4_123 = self.S4contract(k1, k2, k3, D, nocc, norb)
        return S4_123

    def s4(self, k1, k2, k3, D, nocc, norb):
        """

        This code returns the contraction of S[4] for S[4] dict

        :param k1,k2,k3:
            A response matrix
        :param D:
            A density matrix
        :param nocc:
            The number of occupied orbtials
        :param norb:
            The number of total orbitals

        """
        S4_123 = self.S4contract(k1, k2, k3, D, nocc, norb)
        S4_132 = self.S4contract(k1, k3, k2, D, nocc, norb)
        A = S4_123 + S4_132

        return A

    def S4contract(self, k1, k2, k3, D, nocc, norb):
        """

        This code returns the contraction of S[4] for S[4] dict

        :param k1,k2,k3:
            A response matrix
        :param D:
            A density matrix
        :param nocc:
            The number of occupied orbtials
        :param norb:
            The number of total orbitals

        """
        S4N1N2N3 = self.commut(self.commut(k3, self.commut(k2, k1)), D.T)
        S4N1N2N3 = [
            LinearSolver.lrmat2vec(S4N1N2N3.real, nocc, norb),
            LinearSolver.lrmat2vec(S4N1N2N3.imag, nocc, norb)
        ]
        S4N1N2N3_c = S4N1N2N3[0] + 1j * S4N1N2N3[1]
        return (2 / 6) * S4N1N2N3_c

    def flip_xy(self, X):
        # TODO: rewrite flip_xy
        NewXY = []
        for a in range(int(len(X) / 2), len(X), 1):
            NewXY.append(X[a])
        for a in range(0, int(len(X) / 2), 1):
            NewXY.append(X[a])
        NewXY = np.array(NewXY)
        return NewXY

    def flip_yz(self, X):
        # TODO: rewrite flip_yz
        NewXY = []
        for a in range(int(len(X) / 2), len(X), 1):
            NewXY.append(-X[a].real + 1j * X[a].imag)
        for a in range(0, int(len(X) / 2), 1):
            NewXY.append(-X[a].real + 1j * X[a].imag)
        NewXY = np.array(NewXY)
        return NewXY

    def transform_dens(self, k, D, S):
        """
        Creates the perturbed density

        : param k Response vector in matrix form in AO basis
        : param D The density that is to be perturbed in AO basis
        : param S Overlap matrix

        returns [k,D]
        """

        return k.T @ S @ D - D @ S @ k.T

    def mo2ao(self, mo, A):
        """
        Converts a matrix to atomic basis
        : param mo -  molecular orbital coefficent matrix
        : param A - The matrix in MO basis that is the converted to AO basis
        """

        return mo @ A @ mo.T

    def ao2mo(self, mo, A):
        """
        : param mo -  molecular orbital coefficent matrix
        : param A - The matrix in AO basis that is the converted to MO basis
        """

        return (mo.T @ A @ mo)

    def commut(self, A, B):
        """
        commutes two matricies A and B
        """
        return (A @ B - B @ A)

    def x2_contract(self, k, X, D, nocc, norb):
        """
        Contracts the generalized dipole gradient tensor of rank 2 with a second-order response matrix

        : param: k - Respose vector in matrix representation
        : param A - Property operator in matrix represatiation
        : param D - Density matrix
        : param nocc - Number of occupied orbitals
        : param norb - Number of total orbtials

        X[2]N1 = [[k1,X],D.T]

        :return :
            Returns a matrix

        """

        Xn_x = self.commut(self.commut(k, X), D.T)
        X2n_x_c = (LinearSolver.lrmat2vec(Xn_x.real, nocc, norb) +
                   1j * LinearSolver.lrmat2vec(Xn_x.imag, nocc, norb))
        return X2n_x_c

    def x3_contract(self, k1, k2, X, D, nocc, norb):
        """
        Contracts the generalized dipole gradient tensor of rank 3 with two first-order response matrices

        : param: k1,k2:
            First-order response matrix
        : param A:
            A dipole intergral  matrix
        : param D:
            Density matrix
        : param nocc:
            Number of occupied orbitals
        : param norb:
            Number of total orbtials

        X[3]N1N2 = -(1/2)[[k2,[k1,X]],D.T]

        :return :
            Returns a matrix

        """

        X3n_xNy = self.commut(self.commut(k2, self.commut(k1, X)), D.T)
        X3n_xNy = [
            LinearSolver.lrmat2vec(X3n_xNy.real, nocc, norb),
            LinearSolver.lrmat2vec(X3n_xNy.imag, nocc, norb)
        ]
        X3n_xNy_c = X3n_xNy[0] + 1j * X3n_xNy[1]
        return -(1 / 2) * X3n_xNy_c

    def a3_contract(self, k1, k2, A, D, nocc, norb):
        """
        Contracts the generalized dipole gradient tensor of rank 3 with two first-order response matrices

        : param: k1,k2:
            First-order response matrix
        : param A:
            A dipole intergral  matrix
        : param D:
            Density matrix
        : param nocc:
            Number of occupied orbitals
        : param norb:
            Number of total orbtials

        A[3]N1N2 = (1/6)[[k2,[k1,A]],D.T]

        :return :
            Returns a matrix

        """

        A3n_xNy = self.commut(self.commut(k2.T, self.commut(k1.T, A)), D.T)
        A3n_xNy = [
            LinearSolver.lrmat2vec(A3n_xNy.real, nocc, norb),
            LinearSolver.lrmat2vec(A3n_xNy.imag, nocc, norb)
        ]
        A3n_xNy_c = A3n_xNy[0] + 1j * A3n_xNy[1]
        return (1 / 6) * A3n_xNy_c

    def a2_contract(self, k, A, D, nocc, norb):
        """
        Contracts the generalized dipole gradient tensor of rank 2 with a second-order response matrix


        : param: k - Respose vector in matrix representation
        : param A - Property operator in matrix represatiation
        : param D - Density matrix
        : param nocc - Number of occupied orbitals
        : param norb - Number of total orbtials

        A[2]N1 = -(1 / 2)[[k1,X],D.T]

        :return :
            Returns a matrix

        """
        An_x = self.commut(self.commut(k.T, A), D.T)
        A2n_x = [
            LinearSolver.lrmat2vec(An_x.real, nocc, norb),
            LinearSolver.lrmat2vec(An_x.imag, nocc, norb)
        ]
        A2n_x_c = A2n_x[0] + 1j * A2n_x[1]
        return -(1 / 2) * A2n_x_c

    def get_fock_dict_II_red(self, wi, kX, density_list, D0, mo, molecule,
                             ao_basis):
        """
        Computes the compounded second-order Fock matrics used for the reduced isotropic cubic response function

        :param wi:
            A list of the frequencies
        :param kX:
            A dictonary with all the first-order response matrices
        :param D0:
            The SCF density matrix in AO basis
        :param mo:
            A matrix containing the MO coefficents
        :param molecule:
        :param ao_basis:

        :return:
            A dictonary of compounded second-order Fock-matrices

        """
        F123_x = {}
        F123_y = {}
        F123_z = {}

        if self.rank == mpi_master():
            fock_num = 6 * len(wi)
            self.print_header(fock_num)
        time_start_fock = time.time()

        ff_AO = self.get_fock_r(mo, density_list, molecule, ao_basis, 1)

        time_end_fock = time.time()
        total_time_fock = time_end_fock - time_start_fock
        self.print_time(total_time_fock)

        if self.rank == mpi_master():

            fock_list = ff_AO

            count = 0
            for w in wi:
                F123_x.update({w: fock_list[3 * count]})
                F123_y.update({w: fock_list[3 * count + 1]})
                F123_z.update({w: fock_list[3 * count + 2]})
                count += 1

            fock_dict = {'F123_x': F123_x, 'F123_y': F123_y, 'F123_z': F123_z}
        else:
            fock_dict = {}
        return fock_dict

    def get_fock_dict_II(self, wi, kX, density_list, D0, mo, molecule,
                         ao_basis):
        """
        Computes the compounded second-order Fock matrics used for the isotropic cubic response function

        :param wi:
            A list of the frequencies
        :param kX:
            A dictonary with all the first-order response matrices
        :param D0:
            The SCF density matrix in AO basis
        :param mo:
            A matrix containing the MO coefficents
        :param molecule:
        :param ao_basis:

        :return:
            A dictonary of compounded second-order Fock-matrices

        """

        F123_x = {}
        F123_y = {}
        F123_z = {}

        if self.rank == mpi_master():
            fock_num = 6 * len(wi)
            self.print_header(fock_num)
        time_start_fock = time.time()

        ff_AO = self.get_fock_r(mo, density_list, molecule, ao_basis, 1)

        time_end_fock = time.time()
        total_time_fock = time_end_fock - time_start_fock
        self.print_time(total_time_fock)

        if self.rank == mpi_master():

            fock_list = ff_AO

            count = 0
            for w in wi:
                F123_x.update({w: fock_list[3 * count]})
                F123_y.update({w: fock_list[3 * count + 1]})
                F123_z.update({w: fock_list[3 * count + 2]})
                count += 1

            fock_dict = {'F123_x': F123_x, 'F123_y': F123_y, 'F123_z': F123_z}
        else:
            fock_dict = {}
        return fock_dict

    def get_densities_II_red(self, wi, kX, kXY, S, D0, mo):
        """
        Computes the compounded densities needed for the compounded second-order Fock matrics used for the reduced isotropic cubic response function

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
        :param nocc:
            The number of occupied orbitals
        :param norb:
            The number of total orbitals

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

            kx_ = self.mo2ao(mo, kX['Nc'][('x', -w)]).T
            ky_ = self.mo2ao(mo, kX['Nc'][('y', -w)]).T
            kz_ = self.mo2ao(mo, kX['Nc'][('z', -w)]).T

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

            # x #

            Dx = self.transform_dens(kx_, D_sig_xx, S)
            Dx += self.transform_dens(k_sig_xx, Dc_x_, S)
            Dx += self.transform_dens(ky_, D_sig_xy, S)
            Dx += self.transform_dens(k_sig_xy, Dc_y_, S)

            Dx += self.transform_dens(kz_, D_sig_xz, S)
            Dx += self.transform_dens(k_sig_xz, Dc_z_, S)

            # y #
            Dy = self.transform_dens(kx_, D_sig_xy, S)
            Dy += self.transform_dens(k_sig_xy, Dc_x_, S)

            Dy += self.transform_dens(ky_, D_sig_yy, S)
            Dy += self.transform_dens(k_sig_yy, Dc_y_, S)

            Dy += self.transform_dens(kz_, D_sig_yz, S)
            Dy += self.transform_dens(k_sig_yz, Dc_z_, S)

            # z #
            Dz = self.transform_dens(kx_, D_sig_xz, S)
            Dz += self.transform_dens(k_sig_xz, Dc_x_, S)

            Dz += self.transform_dens(ky_, D_sig_yz, S)
            Dz += self.transform_dens(k_sig_yz, Dc_y_, S)

            Dz += self.transform_dens(kz_, D_sig_zz, S)
            Dz += self.transform_dens(k_sig_zz, Dc_z_, S)

            density_list.append(Dx)
            density_list.append(Dy)
            density_list.append(Dz)
        return density_list

    def get_densities_II(self, wi, kX, kXY, S, D0, mo):
        """
        Computes the compounded densities needed for the compounded second-order Fock matrics used for the isotropic cubic response function

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
        :param nocc:
            The number of occupied orbitals
        :param norb:
            The number of total orbitals

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

    def xi(self, kA, kB, Fa, Fb, F0):
        """
        Returns a matrix used for the E[4] contraction

        :param kA,kB:
            First-order response matrices
        :param Fa,Fb:
            First-order perturbed Fock matrices
        :param F0:
            SCF Fock matrix

        """
        xi = 0.5 * (self.commut(kA,
                                self.commut(kB, F0) + 2 * Fb) +
                    self.commut(kB,
                                self.commut(kA, F0) + 2 * Fa))
        return xi

    def phi(self, kA, kB, Fb, F0):
        """
        Returns a matrix used for the E[3] contraction

        :param kA,kB:
            First-order or Second-order response matrices
        :param Fa,Fb:
            First-order or Second-order perturbed Fock matrices
        :param F0:
            SCF Fock matrix
        """
        xi = self.commut(kA, self.commut(kB, F0) + 3 * Fb)
        return xi

    def main(self, Focks, iso, n_x, w, X, d_a_mo, kX, track, scf_tensors,
             molecule, ao_basis):
        """
        This code calls all the relevent functions to third-order isotropic gradient

        :param n_x:
            A dictonary containing all the single index response vectors
        :param w:
            A list of all the frequencies
        :param X:
            A dictonary of matricies containing all the dipole integrals
        :param d_a_mo:
            The SCF density in MO basis
        :param kX:
            A dictonary containing all the response matricies
        :param track:
            A list that contains all the information about which γ components
            and at what freqs they are to be computed
        """

        if self.rank == mpi_master():
            S = scf_tensors['S']
            D0 = scf_tensors['D'][0]
            mo = scf_tensors['C']
            nocc = molecule.number_of_alpha_electrons()
            norb = mo.shape[1]
        else:
            S = None
            D0 = None
            mo = None
            nocc = None
            norb = None

        # computing all compounded first-order densities
        if self.rank == mpi_master():
            density_list = self.get_densities(w, kX, S, D0, mo, nocc, norb)
            density_list_red = self.get_densities_red(w, kX, S, D0, mo, nocc,
                                                      norb)
        else:
            density_list = None
            density_list_red = None

        #  computing the compounded first-order Fock matrices
        fock_dict = self.get_fock_dict(w, kX, density_list, (D0, D0), mo,
                                       molecule, ao_basis)

        fock_dict_red = self.get_fock_dict_red(w, kX, density_list_red,
                                               (D0, D0), mo, molecule, ao_basis)

        if self.rank == mpi_master():
            fock_dict.update(Focks)
            fock_dict_red.update(Focks)
            e4_dict = self.get_e4(w, kX, fock_dict, nocc, norb)

        else:
            e4_dict = {}
            fock_dict = {}
            fock_dict_red = {}

        # computing all the compounded second-order response vectors and
        # extracting some of the second-order Fock matrices from the subspace
        n_xy_dict, kxy_dict, Focks_xy, XΥ_dict = self.get_n_xy(
            w, d_a_mo, X, fock_dict, kX, nocc, norb, molecule, ao_basis,
            scf_tensors)
        n_xy_dict_red, kxy_dict_red, Focks_xy_red, XΥ_dict_red = self.get_n_xy_red(
            w, d_a_mo, X, fock_dict_red, kX, nocc, norb, molecule, ao_basis,
            scf_tensors)

        # computing all second-order compounded densities based on the
        # second-order response vectors
        if self.rank == mpi_master():
            density_list_two = self.get_densities_II(w, kX, kxy_dict, S, D0, mo)
            density_list_two_red = self.get_densities_II_red(
                w, kX, kxy_dict_red, S, D0, mo)
        else:
            n_xy_dict = {}
            kxy_dict = {}
            Focks_xy = {}

            n_xy_dict_red = {}
            kxy_dict_red = {}
            Focks_xy_red = {}
            density_list_two = None
            density_list_two_red = None

        # computing the remaning second-order Fock matrices from the
        # second-order densities
        fock_dict_two = self.get_fock_dict_II(w, kX, density_list_two, D0, mo,
                                              molecule, ao_basis)
        fock_dict_two_red = self.get_fock_dict_II_red(w, kX,
                                                      density_list_two_red, D0,
                                                      mo, molecule, ao_basis)

        if self.rank == mpi_master():
            # Adding the Fock matrices extracted from the second-order response
            # vector subspace to the fock_dict's.
            fock_dict_two.update(Focks_xy)
            fock_dict_two_red.update(Focks_xy_red)

            # computing the compounded E[3] contractions for the isotropic
            # cubic response function
            e3_dict = self.get_e3(w, kX, kxy_dict, fock_dict, fock_dict_two,
                                  nocc, norb)
            e3_dict_red = self.get_e3_red(w, kX, kxy_dict_red, fock_dict_red,
                                          fock_dict_two_red, nocc, norb)

            # computing the X[3],A[3],X[2],A[2] contractions for the isotropic
            # cubic response function
            na_x3_ny_nz, na_a3_nx_ny, na_x2_nyz, nx_a2_nyz = self.other(
                iso, w, track, n_x, n_xy_dict, X, kX, kxy_dict, d_a_mo, nocc,
                norb)
            na_x2_nyz_red, nx_a2_nyz_red = self.other_red(
                iso, w, track, n_x, n_xy_dict_red, X, kX, kxy_dict_red, d_a_mo,
                nocc, norb)

        else:
            e3_dict = {}
            na_x3_ny_nz = None
            na_a3_nx_ny = None
            na_x2_nyz = None
            nx_a2_nyz = None
            e3_dict_red = {}
            na_x2_nyz_red = None
            nx_a2_nyz_red = None

        return (na_x3_ny_nz, na_a3_nx_ny, na_x2_nyz, nx_a2_nyz, e3_dict,
                e4_dict, na_x2_nyz_red, nx_a2_nyz_red, e3_dict_red)

    def get_t3(self, freqs, e3_dict, n_x, track):
        """
        Computes the T[3] contraction, for HF S[3] = 0, R[3] = 0 such that
        the T[3] contraction for the isotropic cubic response function in terms
        of compounded Fock matrices is given as:

                             [(ζ_{α}^{σσ} + ζ_{α}^{λλ+ττ} + f_{α}^{λσ,τ})_is]
        t3term = Σ_{α} N_{α} [(ζ_{α}^{σσ} + ζ_{α}^{λλ+ττ} + f_{α}^{λσ,τ})_si]

        For more details see article

        : param keys:
            Keys from the initial get_densitiesDict
        : param freqs:
            List of frequencies of the pertubations
        : param e3_dict:
            A dictonary that contains the contractions of E[3]
        : param n_x:
            A dictonary containing the response vectors n_x = (E[2]-wS[2])^-1 X[1]
        : param track:
            A list containing information about what tensor components that are being computed

        :return:
            A dictonary of the final values for the NaT[3]NxNyz contractions
        """

        count = 0
        t3_term = {}
        for i in range(len(freqs)):
            w = float(track[i * int(len(track) / len(freqs))].split(",")[1])

            t3term = -n_x['Na'][('x', w)] @ e3_dict['f_iso_x'][w] - n_x['Na'][
                ('y', w)] @ e3_dict['f_iso_y'][w] - n_x['Na'][
                    ('z', w)] @ e3_dict['f_iso_z'][w]

            t3_term.update({(w, -w, w): t3term})

            count += int(len(track) / len(freqs))
        return t3_term

    def get_r4(self, freqs, kX, n_x, track, d_a_mo, nocc, norb):
        """
        Returns a dict with all the R[4]NxNyNz contractions for the subsequent
        T[4] contraction

        : param freqs:
            A list of all the frequencies
        : param kX:
            A dictonary of all the first-order response matrices
        : param n_x:
            A dictonary of all the first-order response vectors
        : param track:
            A list of all the cubic response function components that are to be
            computed for the isotropic
        : param d_a_mo:
            The zeroth-order density in MO basis
        : param nocc:
            The number of occupied orbitals
        : param norb:
            The total number of orbitals

        """

        R4terms = {}
        count = 0

        damp = self.damping

        for j in range(len(freqs)):
            w1 = float(track[j * int(len(track) / len(freqs))].split(",")[1])
            w2 = float(track[j * int(len(track) / len(freqs))].split(",")[2])
            w3 = float(track[j * int(len(track) / len(freqs))].split(",")[3])
            w_s = w1 + w2 + w3

            R4x = 0
            R4y = 0
            R4z = 0

            for i in range(int(len(track) * (1 / len(freqs)))):
                i = i + count

                # Na = n_x['Na'][(track[i][0], w_s)]
                Nb = n_x['Nb'][(track[i][1], w1)]
                Nc = n_x['Nc'][(track[i][2], w2)]
                Nd = n_x['Nd'][(track[i][3], w3)]
                kA = kX['Na'][(track[i][0], w_s)]
                kB = kX['Nb'][(track[i][1], w1)]
                kC = kX['Nc'][(track[i][2], w2)]
                kD = kX['Nd'][(track[i][3], w3)]

                Nb_h = self.flip_xy(Nb)
                Nc_h = self.flip_xy(Nc)
                Nd_h = self.flip_xy(Nd)

                if track[i][0] in 'x':
                    R4x += -1j * damp * Nd_h @ self.s4_for_r4(
                        kA.T, kB, kC, d_a_mo, nocc, norb)
                    R4x += -1j * damp * Nc_h @ self.s4_for_r4(
                        kA.T, kB, kD, d_a_mo, nocc, norb)
                    R4x += -1j * damp * Nd_h @ self.s4_for_r4(
                        kA.T, kC, kB, d_a_mo, nocc, norb)
                    R4x += -1j * damp * Nb_h @ self.s4_for_r4(
                        kA.T, kC, kD, d_a_mo, nocc, norb)
                    R4x += -1j * damp * Nc_h @ self.s4_for_r4(
                        kA.T, kD, kB, d_a_mo, nocc, norb)
                    R4x += -1j * damp * Nb_h @ self.s4_for_r4(
                        kA.T, kD, kC, d_a_mo, nocc, norb)

                if track[i][0] in 'y':
                    R4y += -1j * damp * Nd_h @ self.s4_for_r4(
                        kA.T, kB, kC, d_a_mo, nocc, norb)
                    R4y += -1j * damp * Nc_h @ self.s4_for_r4(
                        kA.T, kB, kD, d_a_mo, nocc, norb)
                    R4y += -1j * damp * Nd_h @ self.s4_for_r4(
                        kA.T, kC, kB, d_a_mo, nocc, norb)
                    R4y += -1j * damp * Nb_h @ self.s4_for_r4(
                        kA.T, kC, kD, d_a_mo, nocc, norb)
                    R4y += -1j * damp * Nc_h @ self.s4_for_r4(
                        kA.T, kD, kB, d_a_mo, nocc, norb)
                    R4y += -1j * damp * Nb_h @ self.s4_for_r4(
                        kA.T, kD, kC, d_a_mo, nocc, norb)

                if track[i][0] in 'z':
                    R4z += -1j * damp * Nd_h @ self.s4_for_r4(
                        kA.T, kB, kC, d_a_mo, nocc, norb)
                    R4z += -1j * damp * Nc_h @ self.s4_for_r4(
                        kA.T, kB, kD, d_a_mo, nocc, norb)
                    R4z += -1j * damp * Nd_h @ self.s4_for_r4(
                        kA.T, kC, kB, d_a_mo, nocc, norb)
                    R4z += -1j * damp * Nb_h @ self.s4_for_r4(
                        kA.T, kC, kD, d_a_mo, nocc, norb)
                    R4z += -1j * damp * Nc_h @ self.s4_for_r4(
                        kA.T, kD, kB, d_a_mo, nocc, norb)
                    R4z += -1j * damp * Nb_h @ self.s4_for_r4(
                        kA.T, kD, kC, d_a_mo, nocc, norb)

            count += int(len(track) / len(freqs))
            R4terms.update({('x', w1): -R4x})
            R4terms.update({('y', w1): -R4y})
            R4terms.update({('z', w1): -R4z})

        return R4terms

    def anti_sym(self, Vec):
        """
        Returns an antisymetrized vector
        :param Vec:
            Vector to me anti-symetrized
        """
        N = len(Vec)
        NewVec = []
        Avec = []
        TotVec = []
        for a in range(int(N / 2)):
            NewVec.append(Vec[a])
        for a in range(int(N / 2), N):
            Avec.append(-Vec[a])
        for a in range(len(NewVec)):
            TotVec.append(NewVec[a])
        for a in range(len(Avec)):
            TotVec.append(Avec[a])
        return TotVec

    def get_n_xy(self, w, d_a_mo, X, fock_dict, kX, nocc, norb, molecule,
                 ao_basis, scf_tensors):
        """
        Computes all the second-order response vectors needed for the isotropic
        cubic response computation

        : param eri_tresh:
            Eri threshold
        : param w:
            A list of all the frequencies
        : param d_a_mo:
            The density matrix in MO basis
        : param X :
            Dipole integrals
        : param fock_dict:
            A dictonary containing all the Fock matricies
        : param kX :
            A dictonary containg all the response matricies
        : param nocc:
            The number of occupied orbitals
        : param norb:
            The number of total orbitals
        :return:
            A dictonary of Fock matrices from the subspace,second-order response vectors and second-order response matrices
        """

        N_total_drv = ComplexResponse(self.comm, self.ostream)

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

        N_total_drv.update_settings({
            'frequencies': freq,
            'damping': self.damping,
            'conv_thresh': self.conv_thresh,
            'lindep_thresh': self.lindep_thresh,
            'max_iter': self.max_iter,
            'eri_thresh': self.eri_thresh,
            'qq_type': self.qq_type,
        })
        if self.batch_size is not None:
            N_total_drv.update_settings({'batch_size': self.batch_size})

        # commutpute second-order response vectors
        N_total_Drv = N_total_drv.compute(molecule, ao_basis, scf_tensors,
                                          xy_dict)

        if self.rank == mpi_master():
            n_xy_dict = N_total_Drv['solutions']
            kxy_dict = N_total_Drv['kappas']
            FXY_2_dict = N_total_Drv['focks']

            return n_xy_dict, kxy_dict, FXY_2_dict, xy_dict
        else:
            return None, None, None, None

    def get_n_xy_red(self, w, d_a_mo, X, fock_dict, kX, nocc, norb, molecule,
                     ao_basis, scf_tensors):
        """
        Computes all the second-order response vectors needed for the reduced
        isotropic cubic response computation

        : param eri_tresh:
            Eri threshold
        : param w:
            A list of all the frequencies
        : param d_a_mo:
            The density matrix in MO basis
        : param X :
            Dipole integrals
        : param fock_dict:
            A dictonary containing all the Fock matricies
        : param kX :
            A dictonary containg all the response matricies
        : param nocc:
            The number of occupied orbitals
        : param norb:
            The number of total orbitals

        :return:
            A dictonary of Fock matrices from the subspace,second-order
            response vectors and second-order response matrices
        """

        xy_dict = {}
        freq = None

        N_total_drv_2 = ComplexResponse(self.comm, self.ostream)

        if self.rank == mpi_master():
            # Get the second-order gradiants
            xy_dict = self.get_xy_red(d_a_mo, X, w, fock_dict, kX, nocc, norb)

            wbd = [sum(x) for x in zip(w, w)]

            freq = ''
            for i in range(len(wbd)):
                freq += str(wbd[i]) + ','

        freq = self.comm.bcast(freq, root=mpi_master())

        N_total_drv_2.update_settings({
            'frequencies': freq,
            'damping': self.damping,
            'conv_thresh': self.conv_thresh,
            'lindep_thresh': self.lindep_thresh,
            'max_iter': self.max_iter,
            'eri_thresh': self.eri_thresh,
            'qq_type': self.qq_type,
        })
        if self.batch_size is not None:
            N_total_drv_2.update_settings({'batch_size': self.batch_size})

        N_total_Drv = N_total_drv_2.compute(molecule, ao_basis, scf_tensors,
                                            xy_dict)

        if self.rank == mpi_master():
            n_xy_dict = N_total_Drv['solutions']
            kxy_dict = N_total_Drv['kappas']
            FXY_2_dict = N_total_Drv['focks']

            return n_xy_dict, kxy_dict, FXY_2_dict, xy_dict
        else:
            return None, None, None, None

    # Fock code

    def get_fock_r(self, mo, D, molecule, ao_basis, rank):
        """
        Computes and returns a list of Fock matrices

        :param D:
            A list of densities
        :param molecule:
        :param ao_basis:
        :param rank:

        """
        if rank == 0:
            # computes the unperturbed Fock matrix from the SCF Density and
            # adds the one electron part
            fa = self.get_two_el_fock_mod_r(molecule, ao_basis, D)
            h = self.get_one_el_hamiltonian(molecule, ao_basis)

            if self.rank == mpi_master():
                fa = 0.5 * fa[0] + h
                return self.ao2mo(mo, fa)
            else:
                return None

        elif rank == 1:
            # computes complex Fock matrices (only two-eletron parts 2J-K)
            if self.rank == mpi_master():
                D_total = []
                for da in D:
                    D_total.append(da.real)
                    D_total.append(da.imag)
            else:
                D_total = None

            f_total = self.get_two_el_fock_mod_r(molecule, ao_basis, D_total)

            ff = []

            if self.rank == mpi_master():
                for i in range(len(f_total) // 2):
                    ff.append(
                        self.ao2mo(
                            mo,
                            0.5 * f_total[2 * i] + 0.5j * f_total[2 * i + 1]))
                return ff
            else:
                return None

        elif rank == 3:
            # computes real Fock Matrices (only two-eletron parts 2J-K)
            f_total = self.get_two_el_fock_mod_r(molecule, ao_basis, D)

            ff = []

            if self.rank == mpi_master():
                for i in range(len(f_total)):
                    ff.append(self.ao2mo(mo, 0.5j * f_total[i]))

                return ff
            else:
                return None

        else:
            # computes imaginary Fock matrices (only two-eletron parts 2J-K)
            f_total = self.get_two_el_fock_mod_r(molecule, ao_basis, D)

            ff = []

            if self.rank == mpi_master():
                for i in range(len(f_total)):
                    ff.append(self.ao2mo(mo, 0.5 * f_total[i]))

                return ff
            else:
                return None

    def get_two_el_fock_mod_r(self, molecule, ao_basis, *dabs):
        """
        Returns the two-electron part of the Fock matix 2J-K

        :param molecule:
        :param ao_basis:
        :param dabs:
            A list of densitiy matrices
        """

        eri_driver = ElectronRepulsionIntegralsDriver(self.comm)
        screening = eri_driver.compute(get_qq_scheme(self.qq_type),
                                       self.eri_thresh, molecule, ao_basis)

        # determine number of batches

        num_batches = 0

        total_mem = psutil.virtual_memory().total
        total_mem_list = self.comm.gather(total_mem, root=mpi_master())

        if self.rank == mpi_master():
            n_ao = dabs[0][0].shape[0]
            n_total = len(dabs[0])

            # check if master node has larger memory
            mem_adjust = 0.0
            if total_mem > min(total_mem_list):
                mem_adjust = total_mem - min(total_mem_list)

            # computes maximum batch size from available memory
            avail_mem = psutil.virtual_memory().available - mem_adjust
            mem_per_mat = n_ao**2 * ctypes.sizeof(ctypes.c_double)
            nthreads = int(os.environ['OMP_NUM_THREADS'])
            max_batch_size = int(avail_mem / mem_per_mat / (0.625 * nthreads))
            max_batch_size = max(1, max_batch_size)

            batch_size = self.batch_size
            if batch_size is None:
                batch_size = min(100, n_total, max_batch_size)

            # get number of batches
            num_batches = n_total // batch_size
            if n_total % batch_size != 0:
                num_batches += 1

        num_batches = self.comm.bcast(num_batches, root=mpi_master())

        # go through batches

        fabs = []

        if self.rank == mpi_master():
            batch_str = 'Processing Fock builds...'
            batch_str += ' (batch size: {:d})'.format(batch_size)
            self.ostream.print_info(batch_str)

        for batch_ind in range(num_batches):

            if self.rank == mpi_master():
                self.ostream.print_info('  batch {}/{}'.format(
                    batch_ind + 1, num_batches))
                self.ostream.flush()

            # form density matrices

            dts = []
            if self.rank == mpi_master():
                batch_start = batch_size * batch_ind
                batch_end = min(batch_start + batch_size, n_total)
                for dab in dabs[0][batch_start:batch_end]:
                    dt = 2 * dab
                    dts.append(dt)
                dens = AODensityMatrix(dts, denmat.rest)
            else:
                dens = AODensityMatrix()

            dens.broadcast(self.rank, self.comm)

            fock = AOFockMatrix(dens)
            for i in range(fock.number_of_fock_matrices()):
                fock.set_fock_type(fockmat.rgenjk, i)

            eri_driver.compute(fock, dens, molecule, ao_basis, screening)
            fock.reduce_sum(self.rank, self.nodes, self.comm)

            if self.rank == mpi_master():
                for i in range(fock.number_of_fock_matrices()):
                    fabs.append(fock.to_numpy(i).T)

        if self.rank == mpi_master():
            return tuple(fabs)
        else:
            return None

    def get_one_el_hamiltonian(self, molecule, ao_basis):
        """
        Returns the one electron part of the Fock matrix

        :param molecule:
        :param ao_basis:
        """

        kinetic_driver = KineticEnergyIntegralsDriver(self.comm)
        potential_driver = NuclearPotentialIntegralsDriver(self.comm)

        T = kinetic_driver.compute(molecule, ao_basis).to_numpy()
        V = potential_driver.compute(molecule, ao_basis).to_numpy()

        return T - V

    def print_header(self, fock_number):
        """
        Fock section setup header for output stream.
        :param fock_number:
           Total number of Fock matrices to compute
        """

        self.ostream.print_blank()
        self.ostream.print_header("Fock matrix computation")
        self.ostream.print_header(31 * "=")
        self.ostream.print_blank()

        width = 50

        cur_str = "Total number of Fock matrices to compute : " + str(
            fock_number)
        self.ostream.print_header(cur_str.ljust(width))

        self.ostream.flush()

    def print_time(self, time):
        """
        Prints time for Fock section
        :param time:
           Total time to compute Fock matrices
        """

        width = 50

        cur_str = "Total time for Fock matrices: " + str(time) + " sec "
        self.ostream.print_header(cur_str.ljust(width))

        self.ostream.flush()
