from .veloxchemlib import OverlapIntegralsDriver
from .veloxchemlib import ElectronRepulsionIntegralsDriver
from .veloxchemlib import KineticEnergyIntegralsDriver
from .veloxchemlib import NuclearPotentialIntegralsDriver
from .veloxchemlib import mpi_master
from .veloxchemlib import denmat, fockmat
from .veloxchemlib import ericut
import time
from .cppsolver import ComplexResponse
from .linearsolver import LinearSolver
from .aofockmatrix import AOFockMatrix
from .aodensitymatrix import AODensityMatrix
from mpi4py import MPI
import numpy as np


class tpa:

    def __init__(self):
        self.comm = MPI.COMM_WORLD
        self.rank = self.comm.Get_rank()
        self.size = self.comm.Get_size()
        self._fock = None
        self.XY_thresh = 1E-14
        self.num_thresh = 1E-12

    def E4(self, wi, kX, FO, nocc, norb):
        """ 
        This code contracts E[3] 

        : param iso - A boolean that specifies if the progam is to compute the isotropic γ or a user specief combination of components
        : param wi - A list of freqs
        : param keys - A dict of lists of keys that give information about what components are present
        : param track - A list containing information about what γ components that are to be computed
        : param kX - A dict of the single index response matricies 
        : param kXY - A dict of the two index response matrices
        : param FO - A dictonary of transformed Fock matricies from Fock_dict
        : param FO2 - A dictonarty of transfromed Fock matricies from Fock_dict_two
        : param nocc - The number of occupied orbitals
        : param norb - The total number of orbitals

        """

        F_iso_x = {}
        F_iso_y = {}
        F_iso_z = {}

        F0 = FO['F0']

        for w in wi:

            kx = kX['Nb'][('x', w)].T
            ky = kX['Nb'][('y', w)].T
            kz = kX['Nb'][('z', w)].T

            kx_ = kX['Nc'][('x', -w)].T
            ky_ = kX['Nc'][('y', -w)].T
            kz_ = kX['Nc'][('z', -w)].T

            Fx = FO['Fb'][('x', w)]
            Fy = FO['Fb'][('y', w)]
            Fz = FO['Fb'][('z', w)]

            Fx_ = FO['Fc'][('x', -w)]
            Fy_ = FO['Fc'][('y', -w)]
            Fz_ = FO['Fc'][('z', -w)]

            F_lamtau_xx = 3 * FO['F_lamtau_xx'][w]
            F_lamtau_yy = 3 * FO['F_lamtau_yy'][w]
            F_lamtau_zz = 3 * FO['F_lamtau_zz'][w]

            F_lamtau_xy = 3 * FO['F_lamtau_xy'][w]
            F_lamtau_xz = 3 * FO['F_lamtau_xz'][w]
            F_lamtau_yz = 3 * FO['F_lamtau_yz'][w]

            F_sig_xx = 3 * FO['F_sig_xx'][w]
            F_sig_yy = 3 * FO['F_sig_yy'][w]
            F_sig_zz = 3 * FO['F_sig_zz'][w]

            F_sig_xy = 3 * FO['F_sig_xy'][w]
            F_sig_xz = 3 * FO['F_sig_xz'][w]
            F_sig_yz = 3 * FO['F_sig_yz'][w]

            Phi_sig_xx = 2 * (3 * self.eps_E4(kx, kx, Fx, F0) + self.eps_E4(
                ky, ky, Fy, F0) + self.eps_E4(kz, kz, Fz, F0))
            Phi_sig_yy = 2 * (self.eps_E4(kx, kx, Fx, F0) +
                              3 * self.eps_E4(ky, ky, Fy, F0) +
                              self.eps_E4(kz, kz, Fz, F0))
            Phi_sig_zz = 2 * (self.eps_E4(kx, kx, Fx, F0) + self.eps_E4(
                ky, ky, Fy, F0) + 3 * self.eps_E4(kz, kz, Fz, F0))

            Phi_sig_xy = 2 * self.eps_E4(kx, ky, Fy, F0) + 2 * self.eps_E4(
                ky, kx, Fx, F0)
            Phi_sig_xz = 2 * self.eps_E4(kx, kz, Fz, F0) + 2 * self.eps_E4(
                kz, kx, Fx, F0)
            Phi_sig_yz = 2 * self.eps_E4(ky, kz, Fz, F0) + 2 * self.eps_E4(
                kz, ky, Fy, F0)

            Phi_lamtau_xx = 3 * self.eps_E4(kx, kx_, Fx_, F0) + 3 * self.eps_E4(
                kx_, kx, Fx, F0) + self.eps_E4(ky, ky_, Fy_, F0) + self.eps_E4(
                    ky_, ky, Fy, F0) + self.eps_E4(
                        kz, kz_, Fz_, F0) + self.eps_E4(kz_, kz, Fz, F0)
            Phi_lamtau_xx = 2 * Phi_lamtau_xx
            Phi_lamtau_yy = self.eps_E4(kx, kx_, Fx_, F0) + self.eps_E4(
                kx_, kx, Fx,
                F0) + 3 * self.eps_E4(ky, ky_, Fy_, F0) + 3 * self.eps_E4(
                    ky_, ky, Fy, F0) + self.eps_E4(
                        kz, kz_, Fz_, F0) + self.eps_E4(kz_, kz, Fz, F0)
            Phi_lamtau_yy = 2 * Phi_lamtau_yy
            Phi_lamtau_zz = self.eps_E4(kx, kx_, Fx_, F0) + self.eps_E4(
                kx_, kx, Fx, F0) + self.eps_E4(ky, ky_, Fy_, F0) + self.eps_E4(
                    ky_, ky, Fy, F0) + 3 * self.eps_E4(
                        kz, kz_, Fz_, F0) + 3 * self.eps_E4(kz_, kz, Fz, F0)
            Phi_lamtau_zz = 2 * Phi_lamtau_zz

            Phi_lamtau_xy = self.eps_E4(kx_, ky, Fy, F0) + self.eps_E4(
                ky, kx_, Fx_, F0) + self.eps_E4(kx, ky_, Fy_, F0) + self.eps_E4(
                    ky_, kx, Fx, F0)
            Phi_lamtau_xy = 2 * Phi_lamtau_xy
            Phi_lamtau_xz = self.eps_E4(kx_, kz, Fz, F0) + self.eps_E4(
                kz, kx_, Fx_, F0) + self.eps_E4(kx, kz_, Fz_, F0) + self.eps_E4(
                    kz_, kx, Fx, F0)
            Phi_lamtau_xz = 2 * Phi_lamtau_xz
            Phi_lamtau_yz = self.eps_E4(ky_, kz, Fz, F0) + self.eps_E4(
                kz, ky_, Fy_, F0) + self.eps_E4(ky, kz_, Fz_, F0) + self.eps_E4(
                    kz_, ky, Fy, F0)
            Phi_lamtau_yz = 2 * Phi_lamtau_yz

            ## x

            F_x = FO['F123_x'][w]
            F_x += self.Com(kx, Phi_lamtau_xx + F_lamtau_xx) + self.Com(
                ky, Phi_lamtau_xy + F_lamtau_xy) + self.Com(
                    kz, Phi_lamtau_xz + F_lamtau_xz)
            F_x += self.Com(kx_, Phi_sig_xx + F_sig_xx) + self.Com(
                ky_, Phi_sig_xy + F_sig_xy) + self.Com(kz_,
                                                       Phi_sig_xz + F_sig_xz)

            F_x = -2 * 1 / 6 * LinearSolver.lrmat2vec(F_x.T, nocc, norb)
            F_x = np.array(self.AntiSym(F_x))
            F_iso_x.update({w: F_x})

            ## y

            F_y = FO['F123_y'][w]
            F_y += self.Com(kx, Phi_lamtau_xy + F_lamtau_xy) + self.Com(
                ky, Phi_lamtau_yy + F_lamtau_yy) + self.Com(
                    kz, Phi_lamtau_yz + F_lamtau_yz)
            F_y += self.Com(kx_, Phi_sig_xy + F_sig_xy) + self.Com(
                ky_, Phi_sig_yy + F_sig_yy) + self.Com(kz_,
                                                       Phi_sig_yz + F_sig_yz)

            F_y = -2 * 1 / 6 * LinearSolver.lrmat2vec(F_y.T, nocc, norb)
            F_y = np.array(self.AntiSym(F_y))

            F_iso_y.update({w: F_y})

            ## z

            F_z = FO['F123_z'][w]
            F_z += self.Com(kx, Phi_lamtau_xz + F_lamtau_xz) + self.Com(
                ky, Phi_lamtau_yz + F_lamtau_yz) + self.Com(
                    kz, Phi_lamtau_zz + F_lamtau_zz)
            F_z += self.Com(kx_, Phi_sig_xz + F_sig_xz) + self.Com(
                ky_, Phi_sig_yz + F_sig_yz) + self.Com(kz_,
                                                       Phi_sig_zz + F_sig_zz)

            F_z = -2 * 1 / 6 * LinearSolver.lrmat2vec(F_z.T, nocc, norb)
            F_z = np.array(self.AntiSym(F_z))
            F_iso_z.update({w: F_z})

        E4_dict = {'F_iso_x': F_iso_x, 'F_iso_y': F_iso_y, 'F_iso_z': F_iso_z}
        return E4_dict

    def E3_red(self, wi, kX, kXY, FO, FO2, nocc, norb):
        """ 
        This code contracts E[3] 

        : param iso - A boolean that specifies if the progam is to compute the isotropic γ or a user specief combination of components
        : param wi - A list of freqs
        : param keys - A dict of lists of keys that give information about what components are present
        : param track - A list containing information about what γ components that are to be computed
        : param kX - A dict of the single index response matricies 
        : param kXY - A dict of the two index response matrices
        : param FO - A dictonary of transformed Fock matricies from Fock_dict
        : param FO2 - A dictonarty of transfromed Fock matricies from Fock_dict_two
        : param nocc - The number of occupied orbitals
        : param norb - The total number of orbitals

        """

        F_iso_x = {}
        F_iso_y = {}
        F_iso_z = {}

        F0_a = FO['F0']

        for w in wi:
            ### Response ####
            k_x_ = kX['Nc'][('x', -w)].T
            k_y_ = kX['Nc'][('y', -w)].T
            k_z_ = kX['Nc'][('z', -w)].T

            k_sig_xx = kXY[(('N_sig_xx', w), 2 * w)].T
            k_sig_yy = kXY[(('N_sig_yy', w), 2 * w)].T
            k_sig_zz = kXY[(('N_sig_zz', w), 2 * w)].T

            k_sig_xy = kXY[(('N_sig_xy', w), 2 * w)].T
            k_sig_xz = kXY[(('N_sig_xz', w), 2 * w)].T
            k_sig_yz = kXY[(('N_sig_yz', w), 2 * w)].T

            ### Focks #####

            F_x_ = FO['Fc'][('x', -w)]
            F_y_ = FO['Fc'][('y', -w)]
            F_z_ = FO['Fc'][('z', -w)]


            F_sig_xx = np.conjugate(FO2[(('N_sig_xx', w), 2 * w)]).T
            F_sig_yy = np.conjugate(FO2[(('N_sig_yy', w), 2 * w)]).T
            F_sig_zz = np.conjugate(FO2[(('N_sig_zz', w), 2 * w)]).T

            F_sig_xy = np.conjugate(FO2[(('N_sig_xy', w), 2 * w)]).T
            F_sig_xz = np.conjugate(FO2[(('N_sig_xz', w), 2 * w)]).T
            F_sig_yz = np.conjugate(FO2[(('N_sig_yz', w), 2 * w)]).T

            ## x

            zeta_sig_xx = self.eps(k_x_, k_sig_xx, F_x_, F_sig_xx, F0_a)
            zeta_sig_yy = self.eps(k_x_, k_sig_yy, F_x_, F_sig_yy, F0_a)
            zeta_sig_zz = self.eps(k_x_, k_sig_zz, F_x_, F_sig_zz, F0_a)

            zeta_sig_xy = self.eps(k_y_, k_sig_xy, F_y_, F_sig_xy, F0_a)
            zeta_sig_xz = self.eps(k_z_, k_sig_xz, F_z_, F_sig_xz, F0_a)

            X_terms = (zeta_sig_xx + zeta_sig_xy +
                       zeta_sig_xz).T + (0.5 * FO2['F123_x'][w]).T
            FF_x = -2 * LinearSolver.lrmat2vec(X_terms, nocc, norb)
            FF_x = np.array(self.AntiSym(FF_x))
            F_iso_x.update({w: FF_x})

            ## y

            zeta_sig_yx = self.eps(k_x_, k_sig_xy, F_x_, F_sig_xy, F0_a)
            zeta_sig_yy = self.eps(k_y_, k_sig_yy, F_y_, F_sig_yy, F0_a)
            zeta_sig_yz = self.eps(k_z_, k_sig_yz, F_z_, F_sig_yz, F0_a)

            Y_terms = (zeta_sig_yx + zeta_sig_yy +
                       zeta_sig_yz).T + (0.5 * FO2['F123_y'][w]).T
            FF_y = -2 * LinearSolver.lrmat2vec(Y_terms, nocc, norb)
            FF_y = np.array(self.AntiSym(FF_y))
            F_iso_y.update({w: FF_y})

            ## z

            zeta_sig_zx = self.eps(k_x_, k_sig_xz, F_x_, F_sig_xz, F0_a)
            zeta_sig_zy = self.eps(k_y_, k_sig_yz, F_y_, F_sig_yz, F0_a)
            zeta_sig_zz = self.eps(k_z_, k_sig_zz, F_z_, F_sig_zz, F0_a)

            Z_terms = (zeta_sig_zx + zeta_sig_zy +
                       zeta_sig_zz).T + (0.5 * FO2['F123_z'][w]).T
            FF_z = -2 * LinearSolver.lrmat2vec(Z_terms, nocc, norb)
            FF_z = np.array(self.AntiSym(FF_z))
            F_iso_z.update({w: FF_z})

        E3_dict = {'F_iso_x': F_iso_x, 'F_iso_y': F_iso_y, 'F_iso_z': F_iso_z}
        return E3_dict

    def E3(self, wi, kX, kXY, FO, FO2, nocc, norb):
        """ 
        This code contracts E[3] 

        : param iso - A boolean that specifies if the progam is to compute the isotropic γ or a user specief combination of components
        : param wi - A list of freqs
        : param keys - A dict of lists of keys that give information about what components are present
        : param track - A list containing information about what γ components that are to be computed
        : param kX - A dict of the single index response matricies 
        : param kXY - A dict of the two index response matrices
        : param FO - A dictonary of transformed Fock matricies from Fock_dict
        : param FO2 - A dictonarty of transfromed Fock matricies from Fock_dict_two
        : param nocc - The number of occupied orbitals
        : param norb - The total number of orbitals

        """

        F_iso_x = {}
        F_iso_y = {}
        F_iso_z = {}

        F0_a = FO['F0']

        for w in wi:
            ### Response ####
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

            ### Focks #####

            F_x_ = FO['Fc'][('x', -w)]
            F_y_ = FO['Fc'][('y', -w)]
            F_z_ = FO['Fc'][('z', -w)]

            F_x = FO['Fb'][('x', w)]
            F_y = FO['Fb'][('y', w)]
            F_z = FO['Fb'][('z', w)]

            F_sig_xx = np.conjugate(FO2[(('N_sig_xx', w), 2 * w)]).T
            F_sig_yy = np.conjugate(FO2[(('N_sig_yy', w), 2 * w)]).T
            F_sig_zz = np.conjugate(FO2[(('N_sig_zz', w), 2 * w)]).T

            F_sig_xy = np.conjugate(FO2[(('N_sig_xy', w), 2 * w)]).T
            F_sig_xz = np.conjugate(FO2[(('N_sig_xz', w), 2 * w)]).T
            F_sig_yz = np.conjugate(FO2[(('N_sig_yz', w), 2 * w)]).T

            F_lamtau_xx = np.conjugate(FO2[(('N_lamtau_xx', w), 0)]).T
            F_lamtau_yy = np.conjugate(FO2[(('N_lamtau_yy', w), 0)]).T
            F_lamtau_zz = np.conjugate(FO2[(('N_lamtau_zz', w), 0)]).T

            F_lamtau_xy = np.conjugate(FO2[(('N_lamtau_xy', w), 0)]).T
            F_lamtau_xz = np.conjugate(FO2[(('N_lamtau_xz', w), 0)]).T
            F_lamtau_yz = np.conjugate(FO2[(('N_lamtau_yz', w), 0)]).T

            ## x

            zeta_sig_xx = self.eps(k_x_, k_sig_xx, F_x_, F_sig_xx, F0_a)
            zeta_sig_yy = self.eps(k_x_, k_sig_yy, F_x_, F_sig_yy, F0_a)
            zeta_sig_zz = self.eps(k_x_, k_sig_zz, F_x_, F_sig_zz, F0_a)

            zeta_sig_xy = self.eps(k_y_, k_sig_xy, F_y_, F_sig_xy, F0_a)
            zeta_sig_xz = self.eps(k_z_, k_sig_xz, F_z_, F_sig_xz, F0_a)

            zeta_lamtau_xx = self.eps(k_x, k_lamtau_xx, F_x, F_lamtau_xx, F0_a)
            zeta_lamtau_yy = self.eps(k_x, k_lamtau_yy, F_x, F_lamtau_yy, F0_a)
            zeta_lamtau_zz = self.eps(k_x, k_lamtau_zz, F_x, F_lamtau_zz, F0_a)

            zeta_lamtau_xy = self.eps(k_y, k_lamtau_xy, F_y, F_lamtau_xy, F0_a)
            zeta_lamtau_xz = self.eps(k_z, k_lamtau_xz, F_z, F_lamtau_xz, F0_a)

            X_terms = (zeta_sig_xx + zeta_sig_xy + zeta_sig_xz).T + (
                zeta_lamtau_xx + zeta_lamtau_xy +
                zeta_lamtau_xz).T + (0.5 * FO2['F123_x'][w]).T
            FF_x = -2 * LinearSolver.lrmat2vec(X_terms, nocc, norb)
            FF_x = np.array(self.AntiSym(FF_x))
            F_iso_x.update({w: FF_x})

            ## y

            zeta_sig_yx = self.eps(k_x_, k_sig_xy, F_x_, F_sig_xy, F0_a)
            zeta_sig_yy = self.eps(k_y_, k_sig_yy, F_y_, F_sig_yy, F0_a)
            zeta_sig_yz = self.eps(k_z_, k_sig_yz, F_z_, F_sig_yz, F0_a)

            zeta_lamtau_yx = self.eps(k_x, k_lamtau_xy, F_x, F_lamtau_xy, F0_a)
            zeta_lamtau_yy = self.eps(k_y, k_lamtau_yy, F_y, F_lamtau_yy, F0_a)
            zeta_lamtau_yz = self.eps(k_z, k_lamtau_yz, F_z, F_lamtau_yz, F0_a)

            Y_terms = (zeta_sig_yx + zeta_sig_yy + zeta_sig_yz).T + (
                zeta_lamtau_yx + zeta_lamtau_yy +
                zeta_lamtau_yz).T + (0.5 * FO2['F123_y'][w]).T
            FF_y = -2 * LinearSolver.lrmat2vec(Y_terms, nocc, norb)
            FF_y = np.array(self.AntiSym(FF_y))
            F_iso_y.update({w: FF_y})

            ## z

            zeta_sig_zx = self.eps(k_x_, k_sig_xz, F_x_, F_sig_xz, F0_a)
            zeta_sig_zy = self.eps(k_y_, k_sig_yz, F_y_, F_sig_yz, F0_a)
            zeta_sig_zz = self.eps(k_z_, k_sig_zz, F_z_, F_sig_zz, F0_a)

            zeta_lamtau_zx = self.eps(k_x, k_lamtau_xz, F_x, F_lamtau_xz, F0_a)
            zeta_lamtau_zy = self.eps(k_y, k_lamtau_yz, F_y, F_lamtau_yz, F0_a)
            zeta_lamtau_zz = self.eps(k_z, k_lamtau_zz, F_z, F_lamtau_zz, F0_a)

            Z_terms = (zeta_sig_zx + zeta_sig_zy + zeta_sig_zz).T + (
                zeta_lamtau_zx + zeta_lamtau_zy +
                zeta_lamtau_zz).T + (0.5 * FO2['F123_z'][w]).T
            FF_z = -2 * LinearSolver.lrmat2vec(Z_terms, nocc, norb)
            FF_z = np.array(self.AntiSym(FF_z))
            F_iso_z.update({w: FF_z})

        E3_dict = {'F_iso_x': F_iso_x, 'F_iso_y': F_iso_y, 'F_iso_z': F_iso_z}
        return E3_dict

    def MakeAllDens_red(self, wi, kX, S, D0, mo, nocc, norb):
        D = {}
        Dens_list = []
        for w in wi:
            ### convert response matrix to ao basis ###
            kx = self.mo2ao(mo, kX['Nb'][('x', w)]).T
            ky = self.mo2ao(mo, kX['Nb'][('y', w)]).T
            kz = self.mo2ao(mo, kX['Nb'][('z', w)]).T

            ### create the first order single indexed densiteies ###
            Dx = self.MakeDens(kx, D0, S)
            Dy = self.MakeDens(ky, D0, S)
            Dz = self.MakeDens(kz, D0, S)

            ### create the first order two indexed densities ###
            ### σ terms ####
            D_sig_xx = 2 * (3 * self.MakeDens(kx, Dx, S) +
                            self.MakeDens(ky, Dy, S) + self.MakeDens(kz, Dz, S))
            D_sig_yy = 2 * (self.MakeDens(kx, Dx, S) + 3 *
                            self.MakeDens(ky, Dy, S) + self.MakeDens(kz, Dz, S))
            D_sig_zz = 2 * (self.MakeDens(kx, Dx, S) + self.MakeDens(ky, Dy, S)
                            + 3 * self.MakeDens(kz, Dz, S))

            D_sig_xy = 2 * (self.MakeDens(ky, Dx, S) + self.MakeDens(kx, Dy, S))
            D_sig_xz = 2 * (self.MakeDens(kx, Dz, S) + self.MakeDens(kz, Dx, S))
            D_sig_yz = 2 * (self.MakeDens(ky, Dz, S) + self.MakeDens(kz, Dy, S))

            ## Create first order three indexed Densities ###

            Dens_list.append(D_sig_xx.real)
            Dens_list.append(D_sig_yy.real)
            Dens_list.append(D_sig_zz.real)
            Dens_list.append(D_sig_xy.real)
            Dens_list.append(D_sig_xz.real)
            Dens_list.append(D_sig_yz.real)

        return Dens_list

    def MakeAllDens(self, wi, kX, S, D0, mo, nocc, norb):
        D = {}
        Dens_list = []
        for w in wi:
            ### convert response matrix to ao basis ###
            kx = self.mo2ao(mo, kX['Nb'][('x', w)]).T
            ky = self.mo2ao(mo, kX['Nb'][('y', w)]).T
            kz = self.mo2ao(mo, kX['Nb'][('z', w)]).T

            kx_ = self.mo2ao(mo, kX['Nc'][('x', -w)]).T
            ky_ = self.mo2ao(mo, kX['Nc'][('y', -w)]).T
            kz_ = self.mo2ao(mo, kX['Nc'][('z', -w)]).T
            ### create the first order single indexed densiteies ###
            Dx = self.MakeDens(kx, D0, S)
            Dy = self.MakeDens(ky, D0, S)
            Dz = self.MakeDens(kz, D0, S)

            Dx_ = self.MakeDens(kx_, D0, S)
            Dy_ = self.MakeDens(ky_, D0, S)
            Dz_ = self.MakeDens(kz_, D0, S)
            ### create the first order two indexed densities ###
            ### σ terms ####
            D_sig_xx = 2 * (3 * self.MakeDens(kx, Dx, S) +
                            self.MakeDens(ky, Dy, S) + self.MakeDens(kz, Dz, S))
            D_sig_yy = 2 * (self.MakeDens(kx, Dx, S) + 3 *
                            self.MakeDens(ky, Dy, S) + self.MakeDens(kz, Dz, S))
            D_sig_zz = 2 * (self.MakeDens(kx, Dx, S) + self.MakeDens(ky, Dy, S)
                            + 3 * self.MakeDens(kz, Dz, S))

            D_sig_xy = 2 * (self.MakeDens(ky, Dx, S) + self.MakeDens(kx, Dy, S))
            D_sig_xz = 2 * (self.MakeDens(kx, Dz, S) + self.MakeDens(kz, Dx, S))
            D_sig_yz = 2 * (self.MakeDens(ky, Dz, S) + self.MakeDens(kz, Dy, S))
            ### λ+τ terms ###
            D_lamtau_xx = 2 * (3 * self.MakeDens(kx_, Dx, S) + self.MakeDens(
                ky_, Dy, S) + self.MakeDens(kz_, Dz, S))
            D_lamtau_xx += 2 * (3 * self.MakeDens(kx, Dx_, S) + self.MakeDens(
                ky, Dy_, S) + self.MakeDens(kz, Dz_, S))

            D_lamtau_yy = 2 * (self.MakeDens(kx_, Dx, S) +
                               3 * self.MakeDens(ky_, Dy, S) +
                               self.MakeDens(kz_, Dz, S))
            D_lamtau_yy += 2 * (self.MakeDens(kx, Dx_, S) +
                                3 * self.MakeDens(ky, Dy_, S) +
                                self.MakeDens(kz, Dz_, S))

            D_lamtau_zz = 2 * (self.MakeDens(kx_, Dx, S) + self.MakeDens(
                ky_, Dy, S) + 3 * self.MakeDens(kz_, Dz, S))
            D_lamtau_zz += 2 * (self.MakeDens(kx, Dx_, S) + self.MakeDens(
                ky, Dy_, S) + 3 * self.MakeDens(kz, Dz_, S))

            D_lamtau_xy = 2 * (self.MakeDens(ky_, Dx, S) +
                               self.MakeDens(kx_, Dy, S))
            D_lamtau_xy += 2 * (self.MakeDens(ky, Dx_, S) +
                                self.MakeDens(kx, Dy_, S))

            D_lamtau_xz = 2 * (self.MakeDens(kx_, Dz, S) +
                               self.MakeDens(kz_, Dx, S))
            D_lamtau_xz += 2 * (self.MakeDens(kx, Dz_, S) +
                                self.MakeDens(kz, Dx_, S))

            D_lamtau_yz = 2 * (self.MakeDens(ky_, Dz, S) +
                               self.MakeDens(kz_, Dy, S))
            D_lamtau_yz += 2 * (self.MakeDens(ky, Dz_, S) +
                                self.MakeDens(kz, Dy_, S))

            ## Create first order three indexed Densities ###

            D_lam_sig_tau_x = self.MakeDens(kx_, D_sig_xx, S) + self.MakeDens(
                ky_, D_sig_xy, S) + self.MakeDens(kz_, D_sig_xz, S)
            D_lam_sig_tau_x += self.MakeDens(
                kx, D_lamtau_xx, S) + self.MakeDens(
                    ky, D_lamtau_xy, S) + self.MakeDens(kz, D_lamtau_xz, S)

            D_lam_sig_tau_y = self.MakeDens(kx_, D_sig_xy, S) + self.MakeDens(
                ky_, D_sig_yy, S) + self.MakeDens(kz_, D_sig_yz, S)
            D_lam_sig_tau_y += self.MakeDens(
                kx, D_lamtau_xy, S) + self.MakeDens(
                    ky, D_lamtau_yy, S) + self.MakeDens(kz, D_lamtau_yz, S)

            D_lam_sig_tau_z = self.MakeDens(kx_, D_sig_xz, S) + self.MakeDens(
                ky_, D_sig_yz, S) + self.MakeDens(kz_, D_sig_zz, S)
            D_lam_sig_tau_z += self.MakeDens(
                kx, D_lamtau_xz, S) + self.MakeDens(
                    ky, D_lamtau_yz, S) + self.MakeDens(kz, D_lamtau_zz, S)

            Dens_list.append(D_sig_xx)
            Dens_list.append(D_sig_yy)
            Dens_list.append(D_sig_zz)
            Dens_list.append(D_sig_xy)
            Dens_list.append(D_sig_xz)
            Dens_list.append(D_sig_yz)

            Dens_list.append(D_lamtau_xx)
            Dens_list.append(D_lamtau_yy)
            Dens_list.append(D_lamtau_zz)
            Dens_list.append(D_lamtau_xy)
            Dens_list.append(D_lamtau_xz)
            Dens_list.append(D_lamtau_yz)

            Dens_list.append(D_lam_sig_tau_x)
            Dens_list.append(D_lam_sig_tau_y)
            Dens_list.append(D_lam_sig_tau_z)

        return Dens_list

    def Fock_dict(self, wi, kX, Dens_list, S, D0, mo, molecule, ao_basis):
        F_sig_xx = {}
        F_sig_yy = {}
        F_sig_zz = {}
        F_sig_xy = {}
        F_sig_xz = {}
        F_sig_yz = {}

        F_lamtau_xx = {}
        F_lamtau_yy = {}
        F_lamtau_zz = {}
        F_lamtau_xy = {}
        F_lamtau_xz = {}
        F_lamtau_yz = {}

        F_lam_sig_tau_x = {}
        F_lam_sig_tau_y = {}
        F_lam_sig_tau_z = {}

        FF_AO = self.get_fock_r(Dens_list, molecule, ao_basis, 1)
        F0_a = self.get_fock_r(D0, molecule, ao_basis, 0)

        if self.rank == mpi_master():
            F0_a = self.ao2mo(mo, F0_a)

            Fock_list = []

            for i in range(len(FF_AO)):
                Fock_list.append(self.ao2mo(mo, FF_AO[i]))

            count = 0
            for w in wi:
                F_sig_xx.update({w: Fock_list[15 * count]})
                F_sig_yy.update({w: Fock_list[15 * count + 1]})
                F_sig_zz.update({w: Fock_list[15 * count + 2]})
                F_sig_xy.update({w: Fock_list[15 * count + 3]})
                F_sig_xz.update({w: Fock_list[15 * count + 4]})
                F_sig_yz.update({w: Fock_list[15 * count + 5]})

                F_lamtau_xx.update({w: Fock_list[15 * count + 6]})
                F_lamtau_yy.update({w: Fock_list[15 * count + 7]})
                F_lamtau_zz.update({w: Fock_list[15 * count + 8]})
                F_lamtau_xy.update({w: Fock_list[15 * count + 9]})
                F_lamtau_xz.update({w: Fock_list[15 * count + 10]})
                F_lamtau_yz.update({w: Fock_list[15 * count + 11]})

                F_lam_sig_tau_x.update({w: Fock_list[15 * count + 12]})
                F_lam_sig_tau_y.update({w: Fock_list[15 * count + 13]})
                F_lam_sig_tau_z.update({w: Fock_list[15 * count + 14]})

                count += 1

            Fock = {
                'F0': F0_a,
                'F_sig_xx': F_sig_xx,
                'F_sig_yy': F_sig_yy,
                'F_sig_zz': F_sig_zz,
                'F_sig_xy': F_sig_xy,
                'F_sig_xz': F_sig_xz,
                'F_sig_yz': F_sig_yz,
                'F_lamtau_xx': F_lamtau_xx,
                'F_lamtau_yy': F_lamtau_yy,
                'F_lamtau_zz': F_lamtau_zz,
                'F_lamtau_xy': F_lamtau_xy,
                'F_lamtau_xz': F_lamtau_xz,
                'F_lamtau_yz': F_lamtau_yz,
                'F123_x': F_lam_sig_tau_x,
                'F123_y': F_lam_sig_tau_y,
                'F123_z': F_lam_sig_tau_z
            }
        else:
            Fock = {}

        return Fock

    def Fock_dict_red(self, wi, kX, Dens_list, S, D0, mo, molecule, ao_basis):

        F_sig_xx = {}
        F_sig_yy = {}
        F_sig_zz = {}
        F_sig_xy = {}
        F_sig_xz = {}
        F_sig_yz = {}

        FF_AO = self.get_fock_r(Dens_list, molecule, ao_basis, 2)
        F0_a = self.get_fock_r(D0, molecule, ao_basis, 0)

        if self.rank == mpi_master():
            F0_a = self.ao2mo(mo, F0_a)

            Fock_list = []

            for i in range(len(FF_AO)):
                Fock_list.append(self.ao2mo(mo, FF_AO[i]))

            count = 0
            for w in wi:
                F_sig_xx.update({w: Fock_list[6 * count]})
                F_sig_yy.update({w: Fock_list[6 * count + 1]})
                F_sig_zz.update({w: Fock_list[6 * count + 2]})
                F_sig_xy.update({w: Fock_list[6 * count + 3]})
                F_sig_xz.update({w: Fock_list[6 * count + 4]})
                F_sig_yz.update({w: Fock_list[6 * count + 5]})

                count += 1

            Fock = {
                'F0': F0_a,
                'F_sig_xx': F_sig_xx,
                'F_sig_yy': F_sig_yy,
                'F_sig_zz': F_sig_zz,
                'F_sig_xy': F_sig_xy,
                'F_sig_xz': F_sig_xz,
                'F_sig_yz': F_sig_yz
            }
        else:
            Fock = {}

        return Fock

    def get_XY(self, d_a_mo, X, wi, Fock, kX, nocc, norb):
        XY_dict = {}
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

                ###### REAL PART #######
                F0 = Fock['F0']
                F_x = Fock['Fb'][('x', w)]
                F_y = Fock['Fb'][('y', w)]
                F_z = Fock['Fb'][('z', w)]

                F_x_ = Fock['Fc'][('x', -w)]
                F_y_ = Fock['Fc'][('y', -w)]
                F_z_ = Fock['Fc'][('z', -w)]

                ### BD σ gradients ####
                XY_dict.update({
                    (('N_sig_xx', w), 2 * w):
                        self.AntiSym(-2 *
                                     LinearSolver.lrmat2vec(
                                         (3 * self.eps(kx, kx, F_x, F_x, F0) +
                                          self.eps(ky, ky, F_y, F_y, F0) +
                                          self.eps(kz, kz, F_z, F_z, F0) + 0.5 *
                                          Fock['F_sig_xx'][w]).T, nocc, norb)) -
                        3 * 2 * self.X2contract(kx.T, mu_x, d_a_mo, nocc, norb)
                        - 2 * self.X2contract(ky.T, mu_y, d_a_mo, nocc, norb) -
                        2 * self.X2contract(kz.T, mu_z, d_a_mo, nocc, norb)
                })
                XY_dict.update({
                    (('N_sig_yy', w), 2 * w):
                        self.AntiSym(-2 *
                                     LinearSolver.lrmat2vec(
                                         (self.eps(kx, kx, F_x, F_x, F0) +
                                          3 * self.eps(ky, ky, F_y, F_y, F0) +
                                          self.eps(kz, kz, F_z, F_z, F0) + 0.5 *
                                          Fock['F_sig_yy'][w]).T, nocc, norb)) -
                        2 * self.X2contract(kx.T, mu_x, d_a_mo, nocc, norb) -
                        3 * 2 * self.X2contract(ky.T, mu_y, d_a_mo, nocc, norb)
                        - 2 * self.X2contract(kz.T, mu_z, d_a_mo, nocc, norb)
                })
                XY_dict.update({
                    (('N_sig_zz', w), 2 * w):
                        self.AntiSym(-2 *
                                     LinearSolver.lrmat2vec(
                                         (self.eps(kx, kx, F_x, F_x, F0) +
                                          self.eps(ky, ky, F_y, F_y, F0) +
                                          3 * self.eps(kz, kz, F_z, F_z, F0) +
                                          0.5 * Fock['F_sig_zz'][w]).T, nocc,
                                         norb))
                        - 2 * self.X2contract(kx.T, mu_x, d_a_mo, nocc, norb) -
                        2 * self.X2contract(ky.T, mu_y, d_a_mo, nocc, norb) -
                        3 * 2 * self.X2contract(kz.T, mu_z, d_a_mo, nocc, norb)
                })

                XY_dict.update({
                    (('N_sig_xy', w), 2 * w):
                        self.AntiSym(-2 * LinearSolver.lrmat2vec(
                            (self.eps(ky, kx, F_y, F_x, F0) +
                             self.eps(kx, ky, F_x, F_y, F0) +
                             0.5 * Fock['F_sig_xy'][w]).T, nocc, norb)) -
                        2 * self.X2contract(ky.T, mu_x, d_a_mo, nocc, norb) -
                        2 * self.X2contract(kx.T, mu_y, d_a_mo, nocc, norb)
                })
                XY_dict.update({
                    (('N_sig_xz', w), 2 * w):
                        self.AntiSym(-2 * LinearSolver.lrmat2vec(
                            (self.eps(kz, kx, F_z, F_x, F0) +
                             self.eps(kx, kz, F_x, F_z, F0) +
                             0.5 * Fock['F_sig_xz'][w]).T, nocc, norb)) -
                        2 * self.X2contract(kz.T, mu_x, d_a_mo, nocc, norb) -
                        2 * self.X2contract(kx.T, mu_z, d_a_mo, nocc, norb)
                })
                XY_dict.update({
                    (('N_sig_yz', w), 2 * w):
                        self.AntiSym(-2 * LinearSolver.lrmat2vec(
                            (self.eps(kz, ky, F_z, F_y, F0) +
                             self.eps(ky, kz, F_y, F_z, F0) +
                             0.5 * Fock['F_sig_yz'][w]).T, nocc, norb)) -
                        2 * self.X2contract(kz.T, mu_y, d_a_mo, nocc, norb) -
                        2 * self.X2contract(ky.T, mu_z, d_a_mo, nocc, norb)
                })

                ### BC CD λ+τ gradients ####
                XY_dict.update(
                    {
                        (('N_lamtau_xx', w), 0):
                            self.AntiSym(
                                -2 * LinearSolver.lrmat2vec(
                                    (3 * 2 * self.eps(kx_, kx, F_x_, F_x, F0) +
                                     2 * self.eps(ky_, ky, F_y_, F_y, F0) +
                                     2 * self.eps(kz_, kz, F_z_, F_z, F0) +
                                     0.5 * Fock['F_lamtau_xx'][w]).T, nocc,
                                    norb))
                            - 3 * 2 *
                            self.X2contract(kx_.T, mu_x, d_a_mo, nocc, norb) -
                            2 * self.X2contract(ky_.T, mu_y, d_a_mo, nocc, norb)
                            -
                            2 * self.X2contract(kz_.T, mu_z, d_a_mo, nocc, norb)
                            - 3 * 2 *
                            self.X2contract(kx.T, mu_x, d_a_mo, nocc, norb) -
                            2 * self.X2contract(ky.T, mu_y, d_a_mo, nocc, norb)
                            -
                            2 * self.X2contract(kz.T, mu_z, d_a_mo, nocc, norb)
                    })

                XY_dict.update(
                    {
                        (('N_lamtau_yy', w), 0):
                            self.AntiSym(
                                -2 * LinearSolver.lrmat2vec(
                                    (2 * self.eps(kx_, kx, F_x_, F_x, F0) +
                                     3 * 2 * self.eps(ky_, ky, F_y_, F_y, F0) +
                                     2 * self.eps(kz_, kz, F_z_, F_z, F0) +
                                     0.5 * Fock['F_lamtau_yy'][w]).T, nocc,
                                    norb))
                            - 2 * self.X2contract(kx_.T, mu_x, d_a_mo,
                                                  nocc, norb) - 3 * 2 *
                            self.X2contract(ky_.T, mu_y, d_a_mo, nocc, norb) -
                            2 * self.X2contract(kz_.T, mu_z, d_a_mo, nocc, norb)
                            -
                            2 * self.X2contract(kx.T, mu_x, d_a_mo, nocc, norb)
                            - 3 * 2 *
                            self.X2contract(ky.T, mu_y, d_a_mo, nocc, norb) -
                            2 * self.X2contract(kz.T, mu_z, d_a_mo, nocc, norb)
                    })

                XY_dict.update(
                    {
                        (('N_lamtau_zz', w), 0):
                            self.AntiSym(
                                -2 * LinearSolver.lrmat2vec(
                                    (2 * self.eps(kx_, kx, F_x_, F_x, F0) +
                                     2 * self.eps(ky_, ky, F_y_, F_y, F0) +
                                     3 * 2 * self.eps(kz_, kz, F_z_, F_z, F0) +
                                     0.5 * Fock['F_lamtau_zz'][w]).T, nocc,
                                    norb))
                            - 2 * self.X2contract(kx_.T, mu_x, d_a_mo,
                                                  nocc, norb) -
                            2 * self.X2contract(ky_.T, mu_y, d_a_mo, nocc, norb)
                            - 3 * 2 *
                            self.X2contract(kz_.T, mu_z, d_a_mo, nocc, norb) -
                            2 * self.X2contract(kx.T, mu_x, d_a_mo, nocc, norb)
                            -
                            2 * self.X2contract(ky.T, mu_y, d_a_mo, nocc, norb)
                            - 3 * 2 *
                            self.X2contract(kz.T, mu_z, d_a_mo, nocc, norb)
                    })

                XY_dict.update({
                    (('N_lamtau_xy', w), 0):
                        self.AntiSym(-2 * LinearSolver.lrmat2vec(
                            (2 * self.eps(ky_, kx, F_y_, F_x, F0) +
                             2 * self.eps(kx_, ky, F_x_, F_y, F0) +
                             0.5 * Fock['F_lamtau_xy'][w]).T, nocc, norb)) -
                        2 * self.X2contract(ky.T, mu_x, d_a_mo, nocc, norb) -
                        2 * self.X2contract(kx.T, mu_y, d_a_mo, nocc, norb) -
                        2 * self.X2contract(ky_.T, mu_x, d_a_mo, nocc, norb) -
                        2 * self.X2contract(kx_.T, mu_y, d_a_mo, nocc, norb)
                })
                XY_dict.update({
                    (('N_lamtau_xz', w), 0):
                        self.AntiSym(-2 * LinearSolver.lrmat2vec(
                            (2 * self.eps(kz_, kx, F_z_, F_x, F0) +
                             2 * self.eps(kx_, kz, F_x_, F_z, F0) +
                             0.5 * Fock['F_lamtau_xz'][w]).T, nocc, norb)) -
                        2 * self.X2contract(kz.T, mu_x, d_a_mo, nocc, norb) -
                        2 * self.X2contract(kx.T, mu_z, d_a_mo, nocc, norb) -
                        2 * self.X2contract(kz_.T, mu_x, d_a_mo, nocc, norb) -
                        2 * self.X2contract(kx_.T, mu_z, d_a_mo, nocc, norb)
                })
                XY_dict.update({
                    (('N_lamtau_yz', w), 0):
                        self.AntiSym(-2 * LinearSolver.lrmat2vec(
                            (2 * self.eps(kz_, ky, F_z_, F_y, F0) +
                             2 * self.eps(ky_, kz, F_y_, F_z, F0) +
                             0.5 * Fock['F_lamtau_yz'][w]).T, nocc, norb)) -
                        2 * self.X2contract(kz.T, mu_y, d_a_mo, nocc, norb) -
                        2 * self.X2contract(ky.T, mu_z, d_a_mo, nocc, norb) -
                        2 * self.X2contract(kz_.T, mu_y, d_a_mo, nocc, norb) -
                        2 * self.X2contract(ky_.T, mu_z, d_a_mo, nocc, norb)
                })
        else:
            XY_dict = {}
        return XY_dict

    def get_XY_red(self, d_a_mo, X, wi, Fock, kX, nocc, norb):
        XY_dict = {}
        if self.rank == mpi_master():
            for w in wi:
                mu_x = X['x']
                mu_y = X['y']
                mu_z = X['z']

                kx = kX['Nb'][('x', w)].T
                ky = kX['Nb'][('y', w)].T
                kz = kX['Nb'][('z', w)].T

                ###### REAL PART #######
                F0 = Fock['F0']
                F_x = Fock['Fb'][('x', w)]
                F_y = Fock['Fb'][('y', w)]
                F_z = Fock['Fb'][('z', w)]


                ### BD σ gradients ####
                XY_dict.update({
                    (('N_sig_xx', w), 2 * w):
                        self.AntiSym(-2 *
                                     LinearSolver.lrmat2vec(
                                         (3 * self.eps(kx, kx, F_x, F_x, F0) +
                                          self.eps(ky, ky, F_y, F_y, F0) +
                                          self.eps(kz, kz, F_z, F_z, F0) + 0.5 *
                                          Fock['F_sig_xx'][w]).T, nocc, norb)) -
                        3 * 2 * self.X2contract(kx.T, mu_x, d_a_mo, nocc, norb)
                        - 2 * self.X2contract(ky.T, mu_y, d_a_mo, nocc, norb) -
                        2 * self.X2contract(kz.T, mu_z, d_a_mo, nocc, norb)
                })
                XY_dict.update({
                    (('N_sig_yy', w), 2 * w):
                        self.AntiSym(-2 *
                                     LinearSolver.lrmat2vec(
                                         (self.eps(kx, kx, F_x, F_x, F0) +
                                          3 * self.eps(ky, ky, F_y, F_y, F0) +
                                          self.eps(kz, kz, F_z, F_z, F0) + 0.5 *
                                          Fock['F_sig_yy'][w]).T, nocc, norb)) -
                        2 * self.X2contract(kx.T, mu_x, d_a_mo, nocc, norb) -
                        3 * 2 * self.X2contract(ky.T, mu_y, d_a_mo, nocc, norb)
                        - 2 * self.X2contract(kz.T, mu_z, d_a_mo, nocc, norb)
                })
                XY_dict.update({
                    (('N_sig_zz', w), 2 * w):
                        self.AntiSym(-2 *
                                     LinearSolver.lrmat2vec(
                                         (self.eps(kx, kx, F_x, F_x, F0) +
                                          self.eps(ky, ky, F_y, F_y, F0) +
                                          3 * self.eps(kz, kz, F_z, F_z, F0) +
                                          0.5 * Fock['F_sig_zz'][w]).T, nocc,
                                         norb))
                        - 2 * self.X2contract(kx.T, mu_x, d_a_mo, nocc, norb) -
                        2 * self.X2contract(ky.T, mu_y, d_a_mo, nocc, norb) -
                        3 * 2 * self.X2contract(kz.T, mu_z, d_a_mo, nocc, norb)
                })

                XY_dict.update({
                    (('N_sig_xy', w), 2 * w):
                        self.AntiSym(-2 * LinearSolver.lrmat2vec(
                            (self.eps(ky, kx, F_y, F_x, F0) +
                             self.eps(kx, ky, F_x, F_y, F0) +
                             0.5 * Fock['F_sig_xy'][w]).T, nocc, norb)) -
                        2 * self.X2contract(ky.T, mu_x, d_a_mo, nocc, norb) -
                        2 * self.X2contract(kx.T, mu_y, d_a_mo, nocc, norb)
                })
                XY_dict.update({
                    (('N_sig_xz', w), 2 * w):
                        self.AntiSym(-2 * LinearSolver.lrmat2vec(
                            (self.eps(kz, kx, F_z, F_x, F0) +
                             self.eps(kx, kz, F_x, F_z, F0) +
                             0.5 * Fock['F_sig_xz'][w]).T, nocc, norb)) -
                        2 * self.X2contract(kz.T, mu_x, d_a_mo, nocc, norb) -
                        2 * self.X2contract(kx.T, mu_z, d_a_mo, nocc, norb)
                })
                XY_dict.update({
                    (('N_sig_yz', w), 2 * w):
                        self.AntiSym(-2 * LinearSolver.lrmat2vec(
                            (self.eps(kz, ky, F_z, F_y, F0) +
                             self.eps(ky, kz, F_y, F_z, F0) +
                             0.5 * Fock['F_sig_yz'][w]).T, nocc, norb)) -
                        2 * self.X2contract(kz.T, mu_y, d_a_mo, nocc, norb) -
                        2 * self.X2contract(ky.T, mu_z, d_a_mo, nocc, norb)
                })

        else:
            XY_dict = {}
        return XY_dict

    def other_red(self, iso, wi, track, Nx, Nxy, X, kX, kXY, da, mol_orbs,
                  task):
        """
        This code computes all the terms of form 
                                + Na(w1+w2+w3)[B[2]Ncd(w2,w3) + C[2]Nbd(w1,w3) + D[2]Nbc(w1,w2)]
                                - Na(w1+w2+3)[(P_23)B[3]Nc(w2)Nd(w3) + (P_13)C[3]Nb(w1)Nd(w3) + (P_12)D[3]Nb(w1)Nc(w2)] <--- other terms (computed with other method)

        : param iso - A bolean value that states if its the isotrpoic gamma or a user defined case
        : param wi - A list containing all the frequencies
        : param keys - A dictonray or lists that are used to tell the program which elements are present
        : param track - A list that contains information about what γ components that are to be computed and which freqs
        : param Nx -  A dictonary containing all the single-index response vectors
        : param Nxy - A dictonary containing all the two-index response vectors
        : param X - A dictonray with all the property integral matricies
        : param kX - A dictonary with all the respone matricies
        : param kXY - A dictonary containing all the two-index response matricies
        : param da - The SCF density matrix in MO bassi
        
        """

        START = time.time()

        NaX2Nyz_dict = {}
        NxA2Nyz_dict = {}

        count = 0

        for i in range(len(wi)):
            wb = float(track[i * int(len(track) / len(wi))].split(",")[1])
            wa = float(track[i * int(len(track) / len(wi))].split(",")[1])
            wc = float(track[i * int(len(track) / len(wi))].split(",")[2])
            w = float(track[i * int(len(track) / len(wi))].split(",")[1])
            wd = float(track[i * int(len(track) / len(wi))].split(",")[3])

            wcd = 0
            wbd = wb + wd
            wbc = 0

            NaX2Nyz = 0
            NxA2Nyz = 0

            ### BD

            kbd = kXY[(('N_sig_xx', w), wbd)]
            Nbd = Nxy[(('N_sig_xx', w), wbd)]
            Na = Nx['Na'][('x', wa)]
            Nc = Nx['Nc'][('x', wc)]
            kc = kX['Nc'][('x', wc)]
            A = X['x']
            C = X['x']

            NaX2Nyz += Na.T @ self.X2contract(kbd, C, da, mol_orbs, task)

            NxA2Nyz += self.A2contract(kc, A, da, mol_orbs, task) @ Nbd
            NxA2Nyz += self.A2contract(kbd, A, da, mol_orbs, task) @ Nc

            kbd = kXY[(('N_sig_xy', w), wbd)]
            Nbd = Nxy[(('N_sig_xy', w), wbd)]
            Na = Nx['Na'][('x', wa)]
            Nc = Nx['Nc'][('y', wc)]
            kc = kX['Nc'][('y', wc)]
            A = X['x']
            C = X['y']

            NaX2Nyz += Na.T @ self.X2contract(kbd, C, da, mol_orbs, task)

            NxA2Nyz += self.A2contract(kc, A, da, mol_orbs, task) @ Nbd
            NxA2Nyz += self.A2contract(kbd, A, da, mol_orbs, task) @ Nc

            kbd = kXY[(('N_sig_xz', w), wbd)]
            Nbd = Nxy[(('N_sig_xz', w), wbd)]
            Nc = Nx['Nc'][('z', wc)]
            kc = kX['Nc'][('z', wc)]
            Na = Nx['Na'][('x', wa)]
            A = X['x']
            C = X['z']

            NaX2Nyz += Na.T @ self.X2contract(kbd, C, da, mol_orbs, task)

            NxA2Nyz += self.A2contract(kc, A, da, mol_orbs, task) @ Nbd
            NxA2Nyz += self.A2contract(kbd, A, da, mol_orbs, task) @ Nc

            ### y

            kbd = kXY[(('N_sig_xy', w), wbd)]
            Nbd = Nxy[(('N_sig_xy', w), wbd)]
            Na = Nx['Na'][('y', wa)]
            Nc = Nx['Nc'][('x', wc)]
            kc = kX['Nc'][('x', wc)]
            A = X['y']
            C = X['x']

            NaX2Nyz += Na.T @ self.X2contract(kbd, C, da, mol_orbs, task)

            NxA2Nyz += self.A2contract(kc, A, da, mol_orbs, task) @ Nbd
            NxA2Nyz += self.A2contract(kbd, A, da, mol_orbs, task) @ Nc

            kbd = kXY[(('N_sig_yy', w), wbd)]
            Nbd = Nxy[(('N_sig_yy', w), wbd)]
            Nc = Nx['Nc'][('y', wc)]
            Na = Nx['Na'][('y', wa)]
            kc = kX['Nc'][('y', wc)]
            A = X['y']
            C = X['y']

            NaX2Nyz += Na.T @ self.X2contract(kbd, C, da, mol_orbs, task)

            NxA2Nyz += self.A2contract(kc, A, da, mol_orbs, task) @ Nbd
            NxA2Nyz += self.A2contract(kbd, A, da, mol_orbs, task) @ Nc

            kbd = kXY[(('N_sig_yz', w), wbd)]
            Nbd = Nxy[(('N_sig_yz', w), wbd)]
            Nc = Nx['Nc'][('z', wc)]
            Na = Nx['Na'][('y', wa)]
            kc = kX['Nc'][('z', wc)]
            A = X['y']
            C = X['z']

            NaX2Nyz += Na.T @ self.X2contract(kbd, C, da, mol_orbs, task)

            NxA2Nyz += self.A2contract(kc, A, da, mol_orbs, task) @ Nbd
            NxA2Nyz += self.A2contract(kbd, A, da, mol_orbs, task) @ Nc

            ### z

            kbd = kXY[(('N_sig_xz', w), wbd)]
            Nbd = Nxy[(('N_sig_xz', w), wbd)]
            Nc = Nx['Nc'][('x', wc)]
            kc = kX['Nc'][('x', wc)]
            Na = Nx['Na'][('z', wa)]
            A = X['z']
            C = X['x']

            NaX2Nyz += Na.T @ self.X2contract(kbd, C, da, mol_orbs, task)

            NxA2Nyz += self.A2contract(kc, A, da, mol_orbs, task) @ Nbd
            NxA2Nyz += self.A2contract(kbd, A, da, mol_orbs, task) @ Nc

            kbd = kXY[(('N_sig_yz', w), wbd)]
            Nbd = Nxy[(('N_sig_yz', w), wbd)]
            Nc = Nx['Nc'][('y', wc)]
            kc = kX['Nc'][('y', wc)]
            Na = Nx['Na'][('z', wa)]
            A = X['z']
            C = X['y']

            NaX2Nyz += Na.T @ self.X2contract(kbd, C, da, mol_orbs, task)

            NxA2Nyz += self.A2contract(kc, A, da, mol_orbs, task) @ Nbd
            NxA2Nyz += self.A2contract(kbd, A, da, mol_orbs, task) @ Nc

            kbd = kXY[(('N_sig_zz', w), wbd)]
            Nbd = Nxy[(('N_sig_zz', w), wbd)]
            Nc = Nx['Nc'][('z', wc)]
            Na = Nx['Na'][('z', wa)]
            kc = kX['Nc'][('z', wc)]
            A = X['z']
            C = X['z']

            NaX2Nyz += Na.T @ self.X2contract(kbd, C, da, mol_orbs, task)

            NxA2Nyz += self.A2contract(kc, A, da, mol_orbs, task) @ Nbd
            NxA2Nyz += self.A2contract(kbd, A, da, mol_orbs, task) @ Nc

            NaX2Nyz_dict.update({(w, -w, w): NaX2Nyz})
            NxA2Nyz_dict.update({(w, -w, w): NxA2Nyz})

        END = time.time()
        TIME_OTHER = END - START
        return NaX2Nyz_dict, NxA2Nyz_dict

    def other(self, iso, wi, track, Nx, Nxy, X, kX, kXY, da, mol_orbs, task):
        """
        This code computes all the terms of form 
                                + Na(w1+w2+w3)[B[2]Ncd(w2,w3) + C[2]Nbd(w1,w3) + D[2]Nbc(w1,w2)]
                                - Na(w1+w2+3)[(P_23)B[3]Nc(w2)Nd(w3) + (P_13)C[3]Nb(w1)Nd(w3) + (P_12)D[3]Nb(w1)Nc(w2)] <--- other terms (computed with other method)

        : param iso - A bolean value that states if its the isotrpoic gamma or a user defined case
        : param wi - A list containing all the frequencies
        : param keys - A dictonray or lists that are used to tell the program which elements are present
        : param track - A list that contains information about what γ components that are to be computed and which freqs
        : param Nx -  A dictonary containing all the single-index response vectors
        : param Nxy - A dictonary containing all the two-index response vectors
        : param X - A dictonray with all the property integral matricies
        : param kX - A dictonary with all the respone matricies
        : param kXY - A dictonary containing all the two-index response matricies
        : param da - The SCF density matrix in MO bassi
        
        """

        START = time.time()

        NaA3NxNy_dict = {}
        NaX3NyNz_dict = {}
        NaX2Nyz_dict = {}
        NxA2Nyz_dict = {}


        count = 0
        for j in range(len(wi)):
            NaX3NyNz = 0
            NaA3NxNy = 0
            for i in range(int(len(track) * (1 / len(wi)))):
                i += count

                w1 = float(track[i].split(",")[1])
                w2 = float(track[i].split(",")[2])
                w3 = float(track[i].split(",")[3])

                Na = Nx['Na'][(track[i][0], w1)]
                Nb = Nx['Nb'][(track[i][1], w1)]

                ####
                Nc = Nx['Nc'][(track[i][2], w2)]
                Nd = Nx['Nd'][(track[i][3], w3)]

                kb = kX['Nb'][(track[i][1], w1)]
                kc = kX['Nc'][(track[i][2], w2)]
                kd = kX['Nd'][(track[i][3], w3)]

                A = X[track[i][0]]
                B = X[track[i][1]]
                C = X[track[i][2]]
                D = X[track[i][3]]

                ## Na X[3]NyNz

                NaX3NyNz += -Na.T @ self.X3contract(kc, kd, B, da, mol_orbs,
                                                    task)
                NaX3NyNz += -Na.T @ self.X3contract(kd, kc, B, da, mol_orbs,
                                                    task)
                NaX3NyNz += -Na.T @ self.X3contract(kd, kb, C, da, mol_orbs,
                                                    task)
                NaX3NyNz += -Na.T @ self.X3contract(kb, kd, C, da, mol_orbs,
                                                    task)
                NaX3NyNz += -Na.T @ self.X3contract(kb, kc, D, da, mol_orbs,
                                                    task)
                NaX3NyNz += -Na.T @ self.X3contract(kc, kb, D, da, mol_orbs,
                                                    task)

                ### NaA[3]NxNy

                NaA3NxNy += self.A3contract(kb, kc, A, da, mol_orbs, task) @ Nd
                NaA3NxNy += self.A3contract(kb, kd, A, da, mol_orbs, task) @ Nc
                NaA3NxNy += self.A3contract(kc, kb, A, da, mol_orbs, task) @ Nd
                NaA3NxNy += self.A3contract(kc, kd, A, da, mol_orbs, task) @ Nb
                NaA3NxNy += self.A3contract(kd, kb, A, da, mol_orbs, task) @ Nc
                NaA3NxNy += self.A3contract(kd, kc, A, da, mol_orbs, task) @ Nb

            count += int(len(track) / len(wi))

            NaA3NxNy_dict.update({(wi[j], -wi[j], wi[j]): NaA3NxNy})
            NaX3NyNz_dict.update({(wi[j], -wi[j], wi[j]): NaX3NyNz})

        for i in range(len(wi)):
            wb = float(track[i * int(len(track) / len(wi))].split(",")[1])
            wa = float(track[i * int(len(track) / len(wi))].split(",")[1])
            wc = float(track[i * int(len(track) / len(wi))].split(",")[2])
            w = float(track[i * int(len(track) / len(wi))].split(",")[1])
            wd = float(track[i * int(len(track) / len(wi))].split(",")[3])

            wcd = 0
            wbd = wb + wd
            wbc = 0

            NaX2Nyz = 0
            NxA2Nyz = 0

            ## CD
            ## x
            kcd = kXY[(('N_lamtau_xx', w), wcd)]
            Ncd = Nxy[(('N_lamtau_xx', w), wcd)]

            Na = Nx['Na'][('x', wa)]
            Nb = Nx['Nb'][('x', w)]
            kb = kX['Nb'][('x', w)]
            A = X['x']
            B = X['x']
            NaX2Nyz += Na.T @ self.X2contract(kcd, B, da, mol_orbs, task)

            NxA2Nyz += self.A2contract(kb, A, da, mol_orbs, task) @ Ncd
            NxA2Nyz += self.A2contract(kcd, A, da, mol_orbs, task) @ Nb

            kcd = kXY[(('N_lamtau_xy', w), wcd)]
            Ncd = Nxy[(('N_lamtau_xy', w), wcd)]

            Nb = Nx['Nb'][('y', w)]
            kb = kX['Nb'][('y', w)]
            Na = Nx['Na'][('x', wa)]
            A = X['x']
            B = X['y']

            NaX2Nyz += Na.T @ self.X2contract(kcd, B, da, mol_orbs, task)

            NxA2Nyz += self.A2contract(kb, A, da, mol_orbs, task) @ Ncd
            NxA2Nyz += self.A2contract(kcd, A, da, mol_orbs, task) @ Nb

            kcd = kXY[(('N_lamtau_xz', w), wcd)]
            Ncd = Nxy[(('N_lamtau_xz', w), wcd)]
            Nb = Nx['Nb'][('z', w)]
            kb = kX['Nb'][('z', w)]
            Na = Nx['Na'][('x', wa)]
            A = X['x']
            B = X['z']

            NaX2Nyz += Na.T @ self.X2contract(kcd, B, da, mol_orbs, task)

            NxA2Nyz += self.A2contract(kb, A, da, mol_orbs, task) @ Ncd
            NxA2Nyz += self.A2contract(kcd, A, da, mol_orbs, task) @ Nb

            ### y

            kcd = kXY[(('N_lamtau_xy', w), wcd)]
            Ncd = Nxy[(('N_lamtau_xy', w), wcd)]
            Nb = Nx['Nb'][('x', w)]
            kb = kX['Nb'][('x', w)]
            A = X['y']
            B = X['x']
            Na = Nx['Na'][('y', wa)]

            NaX2Nyz += Na.T @ self.X2contract(kcd, B, da, mol_orbs, task)

            NxA2Nyz += self.A2contract(kb, A, da, mol_orbs, task) @ Ncd
            NxA2Nyz += self.A2contract(kcd, A, da, mol_orbs, task) @ Nb

            kcd = kXY[(('N_lamtau_yy', w), wcd)]
            Ncd = Nxy[(('N_lamtau_yy', w), wcd)]
            Na = Nx['Na'][('y', wa)]
            Nb = Nx['Nb'][('y', w)]
            kb = kX['Nb'][('y', w)]
            A = X['y']
            B = X['y']

            NaX2Nyz += Na.T @ self.X2contract(kcd, B, da, mol_orbs, task)

            NxA2Nyz += self.A2contract(kb, A, da, mol_orbs, task) @ Ncd
            NxA2Nyz += self.A2contract(kcd, A, da, mol_orbs, task) @ Nb

            kcd = kXY[(('N_lamtau_yz', w), wcd)]
            Ncd = Nxy[(('N_lamtau_yz', w), wcd)]
            Nb = Nx['Nb'][('z', w)]
            kb = kX['Nb'][('z', w)]

            A = X['y']
            B = X['z']
            Na = Nx['Na'][('y', wa)]

            NaX2Nyz += Na.T @ self.X2contract(kcd, B, da, mol_orbs, task)

            NxA2Nyz += self.A2contract(kb, A, da, mol_orbs, task) @ Ncd
            NxA2Nyz += self.A2contract(kcd, A, da, mol_orbs, task) @ Nb

            ### z

            kcd = kXY[(('N_lamtau_xz', w), wcd)]
            Ncd = Nxy[(('N_lamtau_xz', w), wcd)]
            Nb = Nx['Nb'][('x', w)]
            kb = kX['Nb'][('x', w)]
            Na = Nx['Na'][('z', wa)]
            A = X['z']
            B = X['x']

            NaX2Nyz += Na.T @ self.X2contract(kcd, B, da, mol_orbs, task)

            NxA2Nyz += self.A2contract(kb, A, da, mol_orbs, task) @ Ncd
            NxA2Nyz += self.A2contract(kcd, A, da, mol_orbs, task) @ Nb

            kcd = kXY[(('N_lamtau_yz', w), wcd)]
            Ncd = Nxy[(('N_lamtau_yz', w), wcd)]
            Na = Nx['Na'][('z', wa)]
            Nb = Nx['Nb'][('y', w)]
            kb = kX['Nb'][('y', w)]
            A = X['z']
            B = X['y']

            NaX2Nyz += Na.T @ self.X2contract(kcd, B, da, mol_orbs, task)

            NxA2Nyz += self.A2contract(kb, A, da, mol_orbs, task) @ Ncd
            NxA2Nyz += self.A2contract(kcd, A, da, mol_orbs, task) @ Nb

            kcd = kXY[(('N_lamtau_zz', w), wcd)]
            Ncd = Nxy[(('N_lamtau_zz', w), wcd)]
            Nb = Nx['Nb'][('z', w)]
            Na = Nx['Na'][('z', wa)]
            kb = kX['Nb'][('z', w)]
            A = X['z']
            B = X['z']

            NaX2Nyz += Na.T @ self.X2contract(kcd, B, da, mol_orbs, task)

            NxA2Nyz += self.A2contract(kb, A, da, mol_orbs, task) @ Ncd
            NxA2Nyz += self.A2contract(kcd, A, da, mol_orbs, task) @ Nb

            ### BD

            kbd = kXY[(('N_sig_xx', w), wbd)]
            Nbd = Nxy[(('N_sig_xx', w), wbd)]
            Na = Nx['Na'][('x', wa)]
            Nc = Nx['Nc'][('x', wc)]
            kc = kX['Nc'][('x', wc)]
            A = X['x']
            C = X['x']

            NaX2Nyz += Na.T @ self.X2contract(kbd, C, da, mol_orbs, task)

            NxA2Nyz += self.A2contract(kc, A, da, mol_orbs, task) @ Nbd
            NxA2Nyz += self.A2contract(kbd, A, da, mol_orbs, task) @ Nc

            kbd = kXY[(('N_sig_xy', w), wbd)]
            Nbd = Nxy[(('N_sig_xy', w), wbd)]
            Na = Nx['Na'][('x', wa)]
            Nc = Nx['Nc'][('y', wc)]
            kc = kX['Nc'][('y', wc)]
            A = X['x']
            C = X['y']

            NaX2Nyz += Na.T @ self.X2contract(kbd, C, da, mol_orbs, task)

            NxA2Nyz += self.A2contract(kc, A, da, mol_orbs, task) @ Nbd
            NxA2Nyz += self.A2contract(kbd, A, da, mol_orbs, task) @ Nc

            kbd = kXY[(('N_sig_xz', w), wbd)]
            Nbd = Nxy[(('N_sig_xz', w), wbd)]
            Nc = Nx['Nc'][('z', wc)]
            kc = kX['Nc'][('z', wc)]
            Na = Nx['Na'][('x', wa)]
            A = X['x']
            C = X['z']

            NaX2Nyz += Na.T @ self.X2contract(kbd, C, da, mol_orbs, task)

            NxA2Nyz += self.A2contract(kc, A, da, mol_orbs, task) @ Nbd
            NxA2Nyz += self.A2contract(kbd, A, da, mol_orbs, task) @ Nc

            ### y

            kbd = kXY[(('N_sig_xy', w), wbd)]
            Nbd = Nxy[(('N_sig_xy', w), wbd)]
            Na = Nx['Na'][('y', wa)]
            Nc = Nx['Nc'][('x', wc)]
            kc = kX['Nc'][('x', wc)]
            A = X['y']
            C = X['x']

            NaX2Nyz += Na.T @ self.X2contract(kbd, C, da, mol_orbs, task)

            NxA2Nyz += self.A2contract(kc, A, da, mol_orbs, task) @ Nbd
            NxA2Nyz += self.A2contract(kbd, A, da, mol_orbs, task) @ Nc

            kbd = kXY[(('N_sig_yy', w), wbd)]
            Nbd = Nxy[(('N_sig_yy', w), wbd)]
            Nc = Nx['Nc'][('y', wc)]
            Na = Nx['Na'][('y', wa)]
            kc = kX['Nc'][('y', wc)]
            A = X['y']
            C = X['y']

            NaX2Nyz += Na.T @ self.X2contract(kbd, C, da, mol_orbs, task)

            NxA2Nyz += self.A2contract(kc, A, da, mol_orbs, task) @ Nbd
            NxA2Nyz += self.A2contract(kbd, A, da, mol_orbs, task) @ Nc

            kbd = kXY[(('N_sig_yz', w), wbd)]
            Nbd = Nxy[(('N_sig_yz', w), wbd)]
            Nc = Nx['Nc'][('z', wc)]
            Na = Nx['Na'][('y', wa)]
            kc = kX['Nc'][('z', wc)]
            A = X['y']
            C = X['z']

            NaX2Nyz += Na.T @ self.X2contract(kbd, C, da, mol_orbs, task)

            NxA2Nyz += self.A2contract(kc, A, da, mol_orbs, task) @ Nbd
            NxA2Nyz += self.A2contract(kbd, A, da, mol_orbs, task) @ Nc

            ### z

            kbd = kXY[(('N_sig_xz', w), wbd)]
            Nbd = Nxy[(('N_sig_xz', w), wbd)]
            Nc = Nx['Nc'][('x', wc)]
            kc = kX['Nc'][('x', wc)]
            Na = Nx['Na'][('z', wa)]
            A = X['z']
            C = X['x']

            NaX2Nyz += Na.T @ self.X2contract(kbd, C, da, mol_orbs, task)

            NxA2Nyz += self.A2contract(kc, A, da, mol_orbs, task) @ Nbd
            NxA2Nyz += self.A2contract(kbd, A, da, mol_orbs, task) @ Nc

            kbd = kXY[(('N_sig_yz', w), wbd)]
            Nbd = Nxy[(('N_sig_yz', w), wbd)]
            Nc = Nx['Nc'][('y', wc)]
            kc = kX['Nc'][('y', wc)]
            Na = Nx['Na'][('z', wa)]
            A = X['z']
            C = X['y']

            NaX2Nyz += Na.T @ self.X2contract(kbd, C, da, mol_orbs, task)

            NxA2Nyz += self.A2contract(kc, A, da, mol_orbs, task) @ Nbd
            NxA2Nyz += self.A2contract(kbd, A, da, mol_orbs, task) @ Nc

            kbd = kXY[(('N_sig_zz', w), wbd)]
            Nbd = Nxy[(('N_sig_zz', w), wbd)]
            Nc = Nx['Nc'][('z', wc)]
            Na = Nx['Na'][('z', wa)]
            kc = kX['Nc'][('z', wc)]
            A = X['z']
            C = X['z']

            NaX2Nyz += Na.T @ self.X2contract(kbd, C, da, mol_orbs, task)

            NxA2Nyz += self.A2contract(kc, A, da, mol_orbs, task) @ Nbd
            NxA2Nyz += self.A2contract(kbd, A, da, mol_orbs, task) @ Nc

            NaX2Nyz_dict.update({(w, -w, w): NaX2Nyz})
            NxA2Nyz_dict.update({(w, -w, w): NxA2Nyz})

        END = time.time()
        TIME_OTHER = END - START
        return NaX3NyNz_dict, NaA3NxNy_dict, NaX2Nyz_dict, NxA2Nyz_dict

    def T4(self, damp, wi, E4_dict, Nx, Kx, track, da, mol_orbs, task):
        """ 
        This part of the code combines the contraction of the E[4] with that of S[4] and R[4] to return the contraction of T[4]

        : param damp - A scalar that is the damping parameter
        : param keys -  A dict of list of keys which specify which elements are present
        : param wi - A list of all the freqs
        : param E4_dict - A dictonary of all the E[4] contraction
        : param Nx - A dictonary with all the single index response vectors
        : param kX - A dictonray containng all the response matricies
        : param track - A list containg information about all the γ components that are to be computed
        : param da - The SCF density matrix in MO basis
        """
        count = 0
        T4term = {}
        S4 = self.S4_dict(wi, Kx, track, da, mol_orbs, task)

        if damp > 0:
            R4term = self.R4_dict(wi, damp, Kx, Nx, track, da, mol_orbs, task)

        for i in range(len(wi)):

            w = float(track[i * int(len(track) / len(wi))].split(",")[1])
            ww = float(track[i * int(len(track) / len(wi))].split(",")[1])
            t4term = Nx['Na'][('x', w)] @ (E4_dict['F_iso_x'][ww] - S4[
                ('x', ww)]) + Nx['Na'][('y', w)] @ (E4_dict['F_iso_y'][ww] - S4[
                    ('y', ww)]) + Nx['Na'][
                        ('z', w)] @ (E4_dict['F_iso_z'][ww] - S4[('z', ww)])

            if damp > 0:
                t4term += R4term[('x', ww)] + R4term[('y', ww)] + R4term[('z',
                                                                          ww)]

            T4term.update({(ww, -ww, ww): t4term})
            count += int(len(track) / len(wi))
        return T4term

    def S4_dict(self, wi, kX, track, D0, nocc, norb):
        """ 
        This performs the S4 contraction with permutations 

        : param wi - A list of all the freqs
        : param kX - A dict with all the response matricies in MO basis
        : param track - A list containing information about all the components that are to be computed
        : param D0 - The SCF density in MO basis
        : param nocc - The number of occupied obritals
        : param norb - The number of total orbitals

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

                    S4_term_x += w1 * self.S4NNN(kB, kC, kD, D0, nocc, norb)

                    S4_term_x += w2 * self.S4NNN(kC, kB, kD, D0, nocc, norb)

                    S4_term_x += w3 * self.S4NNN(kD, kB, kC, D0, nocc, norb)

                if track[i][0] in 'y':

                    S4_term_y += w1 * self.S4NNN(kB, kC, kD, D0, nocc, norb)

                    S4_term_y += w2 * self.S4NNN(kC, kB, kD, D0, nocc, norb)

                    S4_term_y += w3 * self.S4NNN(kD, kB, kC, D0, nocc, norb)

                if track[i][0] in 'z':

                    S4_term_z += w1 * self.S4NNN(kB, kC, kD, D0, nocc, norb)

                    S4_term_z += w2 * self.S4NNN(kC, kB, kD, D0, nocc, norb)

                    S4_term_z += w3 * self.S4NNN(kD, kB, kC, D0, nocc, norb)

            count += int(len(track) / len(wi))

            S4terms.update({('x', w): -S4_term_x})
            S4terms.update({('y', w): -S4_term_y})
            S4terms.update({('z', w): -S4_term_z})

        return S4terms

    def S4nnn(self, k1, k2, k3, D, nocc, norb):
        """ This code returns the contraction of S[4] for the contraction of R[4] """
        S4_123 = self.S4contract(k1, k2, k3, D, nocc, norb)
        return S4_123

    def S4NNN(self, k1, k2, k3, D, nocc, norb):
        """ This code returns the contraction of S[4] for S[4] dict"""
        S4_123 = self.S4contract(k1, k2, k3, D, nocc, norb)
        S4_132 = self.S4contract(k1, k3, k2, D, nocc, norb)
        A = S4_123 + S4_132
        return A

    def S4contract(self, k1, k2, k3, D, nocc, norb):
        S4N1N2N3 = self.Com(self.Com(k3, self.Com(k2, k1)), D.T)
        S4N1N2N3 = [
            LinearSolver.lrmat2vec(S4N1N2N3.real, nocc, norb),
            LinearSolver.lrmat2vec(S4N1N2N3.imag, nocc, norb)
        ]
        S4N1N2N3_c = S4N1N2N3[0] + 1j * S4N1N2N3[1]
        return (2 / 6) * S4N1N2N3_c

    def FlipXY(self, X):
        NewXY = []
        for a in range(int(len(X) / 2), len(X), 1):
            NewXY.append(X[a])
        for a in range(0, int(len(X) / 2), 1):
            NewXY.append(X[a])
        NewXY = np.array(NewXY)
        return NewXY

    def FlipYZ(self, X):
        NewXY = []
        for a in range(int(len(X) / 2), len(X), 1):
            NewXY.append(-X[a].real + 1j * X[a].imag)  ## Append (Y)
        for a in range(0, int(len(X) / 2), 1):
            NewXY.append(-X[a].real + 1j * X[a].imag)
        NewXY = np.array(NewXY)
        return NewXY

    def MakeDens(self, k, D, S):
        """ 
        Creates the perturbed density 

        : param k Response vector in matrix form 
        : param D The density that is to be perturbed
        : param S Overlap matrix

        returns [k,D]
        """

        da = D
        da_a = k.T @ S @ da - da @ S @ k.T
        return np.array(da_a)

    def mo2ao(self, mo, A):
        """ 
        Converts a matrix to atomic basis 
        : param mo -  molecular orbital coefficent matrix 
        : param A - The matrix in MO basis that is the converted to AO basis
        """

        return mo @ A @ mo.T

    def ao2mo_inv(self, mo, A):
        """
        : param mo - molecular orbital coefficent matrix
        """
        return np.linalg.inv(mo) @ A @ np.linalg.inv(mo.T)

    def ao2mo(self, mo, A):
        """
        : param mo -  molecular orbital coefficent matrix 
        : param A - The matrix in AO basis that is the converted to MO basis
        """

        return (mo.T @ A @ mo)

    def Com(self, A, B):
        """ 
        Commutes two matricies A and B
        """
        return (A @ B - B @ A)

    def X2contract(self, k, X, D, nocc, norb):
        """
        Contracts a property matrix with a response vector in their matrix form

        : param: k - Respose vector in matrix representation
        : param X - Property operator in matrix represatiation
        : param D - Density matrix
        : param nocc - Number of occupied orbitals
        : param norb - Number of total orbtials
        
        X[2]Nx = [[k,X],D.T]  <--- This is a matrix

        The function LinearSolver.lrmat2vec then uses the number of occupied and the total number orbitals to 
        extract the relevent components of matrix and rewrite it as a vector.

        returns the vector that results from the contraction of the rank 2 tensor

        """

        XNx = self.Com(self.Com(k, X), D.T)
        X2Nx = [
            LinearSolver.lrmat2vec(XNx.real, nocc, norb),
            LinearSolver.lrmat2vec(XNx.imag, nocc, norb)
        ]
        X2Nx_c = X2Nx[0] + 1j * X2Nx[1]
        return X2Nx_c

    def X3contract(self, k1, k2, X, D, nocc, norb):
        """
        Contracts a property tensor of rank 3 with two response vectors in their  matrix representation

        : param: k1 - Respose vector in matrix representation
        : param: k2 - Respose vector in matrix representation
        : param X - Property operator in matrix represatiation
        : param D - Density matrix
        : param nocc - Number of occupied orbitals
        : param norb - Number of total orbtials
        
        X[3]NxNy = -(1/2)[[k2,[k1,X]],D.T]  <--- This is a matrix

        The function LinearSolver.lrmat2vec then uses the number of occupied and the total number orbitals to 
        extract the relevent components of matrix and rewrite it as a vector.

        returns the vector that results from the contraction of the rank 3 tensor

        """

        X3NxNy = self.Com(self.Com(k2, self.Com(k1, X)), D.T)
        X3NxNy = [
            LinearSolver.lrmat2vec(X3NxNy.real, nocc, norb),
            LinearSolver.lrmat2vec(X3NxNy.imag, nocc, norb)
        ]
        X3NxNy_c = X3NxNy[0] + 1j * X3NxNy[1]
        return -(1 / 2) * X3NxNy_c

    def A3contract(self, k1, k2, A, D, nocc, norb):
        """
        Contracts a property tensor of rank 3 with two response vectors in their matrix representation

        : param: k1 - Respose vector in matrix representation
        : param: k2 - Respose vector in matrix representation
        : param A - Property operator in matrix represatiation
        : param D - Density matrix
        : param nocc - Number of occupied orbitals
        : param norb - Number of total orbtials
        
        A[3]NxNy = (1/6)[[k2,[k1,A]],D.T]  <--- This is a matrix

        The function LinearSolver.lrmat2vec then uses the number of occupied and the total number orbitals to 
        extract the relevent components of matrix and rewrite it as a vector.

        returns the vector that results from the contraction of the rank 3 tensor

        """

        A3NxNy = self.Com(self.Com(k2.T, self.Com(k1.T, A)), D.T)
        A3NxNy = [
            LinearSolver.lrmat2vec(A3NxNy.real, nocc, norb),
            LinearSolver.lrmat2vec(A3NxNy.imag, nocc, norb)
        ]
        A3NxNy_c = A3NxNy[0] + 1j * A3NxNy[1]
        return (1 / 6) * A3NxNy_c

    def A2contract(self, k, A, D, nocc, norb):
        """
        Contracts a property matrix with a response vector in their matrix form

        : param: k - Respose vector in matrix representation
        : param A - Property operator in matrix represatiation 
        : param D - Density matrix
        : param nocc - Number of occupied orbitals
        : param norb - Number of total orbtials
        
        A[2]Nx = (-1/2)*[[k,A],D.T]  <--- This is a matrix

        The function LinearSolver.lrmat2vec then uses the number of occupied and the total number orbitals to 
        extract the relevent components of matrix and rewrite it as a vector.

        returns the vector that results from the contraction of the rank 2 tensor

        """
        ANx = self.Com(self.Com(k.T, A), D.T)
        A2Nx = [
            LinearSolver.lrmat2vec(ANx.real, nocc, norb),
            LinearSolver.lrmat2vec(ANx.imag, nocc, norb)
        ]
        A2Nx_c = A2Nx[0] + 1j * A2Nx[1]
        return -(1 / 2) * A2Nx_c

    def Fock_dict_II_red(self, wi, kX, Dens_list, S, D0, mo, molecule,
                         ao_basis):

        F123_x = {}
        F123_y = {}
        F123_z = {}

        start = time.time()
        FF_AO = self.get_fock_r(Dens_list, molecule, ao_basis, 3)
        end = time.time()

        if self.rank == mpi_master():
            #print("Number of Focks for T[3] contraction")
            #print(len(Dens_list))
            Fock_list = []
            for i in range(len(FF_AO)):
                Fock_list.append(self.ao2mo(mo, FF_AO[i]))

            count = 0
            for w in wi:
                F123_x.update({w: Fock_list[3 * count]})
                F123_y.update({w: Fock_list[3 * count + 1]})
                F123_z.update({w: Fock_list[3 * count + 2]})
                count += 1

            Fock_dict = {'F123_x': F123_x, 'F123_y': F123_y, 'F123_z': F123_z}
        else:
            Fock_dict = {}
        return Fock_dict

    def Fock_dict_II(self, wi, kX, Dens_list, S, D0, mo, molecule, ao_basis):

        F123_x = {}
        F123_y = {}
        F123_z = {}

        start = time.time()
        FF_AO = self.get_fock_r(Dens_list, molecule, ao_basis, 1)
        end = time.time()

        if self.rank == mpi_master():
            print("Number of Focks for T[3] contraction")
            print(len(Dens_list))
            Fock_list = []
            for i in range(len(FF_AO)):
                Fock_list.append(self.ao2mo(mo, FF_AO[i]))

            count = 0
            for w in wi:
                F123_x.update({w: Fock_list[3 * count]})
                F123_y.update({w: Fock_list[3 * count + 1]})
                F123_z.update({w: Fock_list[3 * count + 2]})
                count += 1

            Fock_dict = {'F123_x': F123_x, 'F123_y': F123_y, 'F123_z': F123_z}
        else:
            Fock_dict = {}
        return Fock_dict

    def MakeAllDens_II_red(self, wi, kX, kXY, S, D0, mo):
        Dens_list = []
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


            ### SIGMA contributiatons ###
            Dc_x_ = self.MakeDens(kx_, D0, S)
            Dc_y_ = self.MakeDens(ky_, D0, S)
            Dc_z_ = self.MakeDens(kz_, D0, S)

            D_sig_xx = self.MakeDens(k_sig_xx, D0, S)
            D_sig_yy = self.MakeDens(k_sig_yy, D0, S)
            D_sig_zz = self.MakeDens(k_sig_zz, D0, S)

            D_sig_xy = self.MakeDens(k_sig_xy, D0, S)
            D_sig_xz = self.MakeDens(k_sig_xz, D0, S)
            D_sig_yz = self.MakeDens(k_sig_yz, D0, S)

            ### x ####

            Dx = self.MakeDens(kx_, D_sig_xx, S)
            Dx += self.MakeDens(k_sig_xx, Dc_x_, S)
            Dx += self.MakeDens(ky_, D_sig_xy, S)
            Dx += self.MakeDens(k_sig_xy, Dc_y_, S)

            Dx += self.MakeDens(kz_, D_sig_xz, S)
            Dx += self.MakeDens(k_sig_xz, Dc_z_, S)

            ## y ###
            Dy = self.MakeDens(kx_, D_sig_xy, S)
            Dy += self.MakeDens(k_sig_xy, Dc_x_, S)

            Dy += self.MakeDens(ky_, D_sig_yy, S)
            Dy += self.MakeDens(k_sig_yy, Dc_y_, S)

            Dy += self.MakeDens(kz_, D_sig_yz, S)
            Dy += self.MakeDens(k_sig_yz, Dc_z_, S)

            ## z ###
            Dz = self.MakeDens(kx_, D_sig_xz, S)
            Dz += self.MakeDens(k_sig_xz, Dc_x_, S)

            Dz += self.MakeDens(ky_, D_sig_yz, S)
            Dz += self.MakeDens(k_sig_yz, Dc_y_, S)

            Dz += self.MakeDens(kz_, D_sig_zz, S)
            Dz += self.MakeDens(k_sig_zz, Dc_z_, S)

            Dens_list.append(Dx.imag)
            Dens_list.append(Dy.imag)
            Dens_list.append(Dz.imag)
        return Dens_list

    def MakeAllDens_II(self, wi, kX, kXY, S, D0, mo):
        Dens_list = []
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

            ### SIGMA contributiatons ###
            Dc_x_ = self.MakeDens(kx_, D0, S)
            Dc_y_ = self.MakeDens(ky_, D0, S)
            Dc_z_ = self.MakeDens(kz_, D0, S)

            D_sig_xx = self.MakeDens(k_sig_xx, D0, S)
            D_sig_yy = self.MakeDens(k_sig_yy, D0, S)
            D_sig_zz = self.MakeDens(k_sig_zz, D0, S)

            D_sig_xy = self.MakeDens(k_sig_xy, D0, S)
            D_sig_xz = self.MakeDens(k_sig_xz, D0, S)
            D_sig_yz = self.MakeDens(k_sig_yz, D0, S)

            ### LAM + TAU contributiatons ###
            Dc_x = self.MakeDens(kx, D0, S)
            Dc_y = self.MakeDens(ky, D0, S)
            Dc_z = self.MakeDens(kz, D0, S)

            D_lamtau_xx = self.MakeDens(k_lamtau_xx, D0, S)
            D_lamtau_yy = self.MakeDens(k_lamtau_yy, D0, S)
            D_lamtau_zz = self.MakeDens(k_lamtau_zz, D0, S)

            D_lamtau_xy = self.MakeDens(k_lamtau_xy, D0, S)
            D_lamtau_xz = self.MakeDens(k_lamtau_xz, D0, S)
            D_lamtau_yz = self.MakeDens(k_lamtau_yz, D0, S)

            ### x ####

            Dx = self.MakeDens(kx_, D_sig_xx, S)
            Dx += self.MakeDens(k_sig_xx, Dc_x_, S)
            Dx += self.MakeDens(ky_, D_sig_xy, S)
            Dx += self.MakeDens(k_sig_xy, Dc_y_, S)

            Dx += self.MakeDens(kz_, D_sig_xz, S)
            Dx += self.MakeDens(k_sig_xz, Dc_z_, S)

            Dx += self.MakeDens(kx, D_lamtau_xx, S)
            Dx += self.MakeDens(k_lamtau_xx, Dc_x, S)

            Dx += self.MakeDens(ky, D_lamtau_xy, S)
            Dx += self.MakeDens(k_lamtau_xy, Dc_y, S)

            Dx += self.MakeDens(kz, D_lamtau_xz, S)
            Dx += self.MakeDens(k_lamtau_xz, Dc_z, S)

            ## y ###
            Dy = self.MakeDens(kx_, D_sig_xy, S)
            Dy += self.MakeDens(k_sig_xy, Dc_x_, S)

            Dy += self.MakeDens(ky_, D_sig_yy, S)
            Dy += self.MakeDens(k_sig_yy, Dc_y_, S)

            Dy += self.MakeDens(kz_, D_sig_yz, S)
            Dy += self.MakeDens(k_sig_yz, Dc_z_, S)

            Dy += self.MakeDens(kx, D_lamtau_xy, S)
            Dy += self.MakeDens(k_lamtau_xy, Dc_x, S)

            Dy += self.MakeDens(ky, D_lamtau_yy, S)
            Dy += self.MakeDens(k_lamtau_yy, Dc_y, S)

            Dy += self.MakeDens(kz, D_lamtau_yz, S)
            Dy += self.MakeDens(k_lamtau_yz, Dc_z, S)

            ## z ###
            Dz = self.MakeDens(kx_, D_sig_xz, S)
            Dz += self.MakeDens(k_sig_xz, Dc_x_, S)

            Dz += self.MakeDens(ky_, D_sig_yz, S)
            Dz += self.MakeDens(k_sig_yz, Dc_y_, S)

            Dz += self.MakeDens(kz_, D_sig_zz, S)
            Dz += self.MakeDens(k_sig_zz, Dc_z_, S)

            Dz += self.MakeDens(kx, D_lamtau_xz, S)
            Dz += self.MakeDens(k_lamtau_xz, Dc_x, S)

            Dz += self.MakeDens(ky, D_lamtau_yz, S)
            Dz += self.MakeDens(k_lamtau_yz, Dc_y, S)

            Dz += self.MakeDens(kz, D_lamtau_zz, S)
            Dz += self.MakeDens(k_lamtau_zz, Dc_z, S)

            Dens_list.append(Dx)
            Dens_list.append(Dy)
            Dens_list.append(Dz)
        return Dens_list

    def eps(self, kA, kB, Fa, Fb, F0):
        eps = 0.5 * (self.Com(kA,
                              self.Com(kB, F0) + 2 * Fb) +
                     self.Com(kB,
                              self.Com(kA, F0) + 2 * Fa))
        return eps

    def eps_E4(self, kA, kB, Fb, F0):
        eps = self.Com(kA, self.Com(kB, F0) + 3 * Fb)
        return eps

    def main(self, eri_thresh, conv_thresh, lindep_thresh, max_iter, Focks, iso,
             Nx, w, X, damp, d_a_mo, kX, track, S, D0, mo, nocc, norb,
             scf_tensors, molecule, ao_basis, comm, ostream):
        """ This code calls all the relevent functions to compute γ 
        : param iso - A boolean which states if its the isotropic case or if its a user specied components that are to be computed
        : param Nx - A dictonary containing all the single index response vectors
        : param w - a list of all the frequencies 
        : param X - A dict of matricies for all the property integrals
        : param damp - the damping constant
        : param d_a_mo - the SCF density in MO basis
        : param kX - A dictonary containing all the response matricies
        : param track - A list that contains all the information about which γ components and at what freqs they are to be computed
        : param S - The overlap matrix
        : param mo - MO coeff matrix
        : param nocc - The number of occupied orbitals
        : param norb - The number of total orbitals

        """

        if self.rank == mpi_master():
            Dens_list = self.MakeAllDens(w, kX, S, D0, mo, nocc, norb)
            Dens_list_red = self.MakeAllDens_red(w, kX, S, D0, mo, nocc, norb)
        else:
            Dens_list = None
            Dens_list_red = None
            Dens_dict = {}

        Fock_dict = self.Fock_dict(w, kX, Dens_list, S, (D0, D0), mo, molecule,
                                   ao_basis)
        Fock_dict_red = self.Fock_dict_red(w, kX, Dens_list_red, S, (D0, D0),
                                           mo, molecule, ao_basis)

        if self.rank == mpi_master():
            Fock_dict.update(Focks)
            Fock_dict_red.update(Focks)
            E4_dict = self.E4(w, kX, Fock_dict, nocc, norb)

        else:
            E4_dict = {}
            Fock_dict = {}
            Fock_dict_red = {}

        Nxy_dict, kXY_dict, Focks_xy, XΥ_dict = self.Nxy(
            eri_thresh, conv_thresh, lindep_thresh, max_iter, w, d_a_mo, damp,
            X, Fock_dict, kX, nocc, norb, molecule, ao_basis, scf_tensors, comm,
            ostream)
        Nxy_dict_red, kXY_dict_red, Focks_xy_red, XΥ_dict_red = self.Nxy_red(
            eri_thresh, conv_thresh, lindep_thresh, max_iter, w, d_a_mo, damp,
            X, Fock_dict_red, kX, nocc, norb, molecule, ao_basis, scf_tensors,
            comm, ostream)

        if self.rank == mpi_master():
            Dens_list_two = self.MakeAllDens_II(w, kX, kXY_dict, S, D0, mo)
            Dens_list_two_red = self.MakeAllDens_II_red(w, kX, kXY_dict_red, S,
                                                        D0, mo)
        else:
            Nxy_dict = {}
            kXY_dict = {}
            Focks_xy = {}

            Nxy_dict_red = {}
            kXY_dict_red = {}
            Focks_xy_red = {}
            Dens_list_two = None
            Dens_list_two_red = None

        Fock_dict_two = self.Fock_dict_II(w, kX, Dens_list_two, S, D0, mo,
                                          molecule, ao_basis)
        Fock_dict_two_red = self.Fock_dict_II_red(w, kX, Dens_list_two_red, S,
                                                  D0, mo, molecule, ao_basis)

        if self.rank == mpi_master():
            Fock_dict_two.update(Focks_xy)
            Fock_dict_two_red.update(Focks_xy_red)
            E3_dict = self.E3(w, kX, kXY_dict, Fock_dict, Fock_dict_two, nocc,
                              norb)
            E3_dict_red = self.E3_red(w, kX, kXY_dict_red, Fock_dict_red,
                                      Fock_dict_two_red, nocc, norb)

            NaX3NyNz, NaA3NxNy, NaX2Nyz, NxA2Nyz = self.other(
                iso, w, track, Nx, Nxy_dict, X, kX, kXY_dict, d_a_mo, nocc,
                norb)
            NaX2Nyz_red, NxA2Nyz_red = self.other_red(iso, w, track, Nx,
                                                      Nxy_dict_red, X, kX,
                                                      kXY_dict_red, d_a_mo,
                                                      nocc, norb)

        else:
            E3_dict = {}
            NaX3NyNz = None
            NaA3NxNy = None
            NaX2Nyz = None
            NxA2Nyz = None
            E3_dict_red = {}
            NaX2Nyz_red = None
            NxA2Nyz_red = None

        return NaX3NyNz, NaA3NxNy, NaX2Nyz, NxA2Nyz, E3_dict, E4_dict, NaX2Nyz_red, NxA2Nyz_red, E3_dict_red

    def T3(self, freqs, E3_dict, Nx, track):
        """
        Computes the  Na(x,w)(T[3]NbNcd + T[3]NcNbd + T[3]NdNbc) + Na(y,w)(T[3]NbNcd + T[3]NcNbd + T[3]NdNbc) + Na(z,w)(T[3]NbNcd + T[3]NcNbd + T[3]NdNbc)

        : param keys:
            Keys from the initial MakeAllDensDict
        : param freqs:
            List of frequencies of the pertubations
        : param E3_dict:
            A dictonary that contains the contractions of E[3]
        : param Nx:
            A dictonary containing the response vectors Nx = (E[2]-wS[2])^-1 X[1]
        : param track:
            A list containing information about what tensor components that are being computed
        """

        count = 0
        T3term = {}
        for i in range(len(freqs)):
            w = float(track[i * int(len(track) / len(freqs))].split(",")[1])

            t3term = -Nx['Na'][('x', w)] @ E3_dict['F_iso_x'][w] - Nx['Na'][
                ('y', w)] @ E3_dict['F_iso_y'][w] - Nx['Na'][
                    ('z', w)] @ E3_dict['F_iso_z'][w]


            T3term.update({(w, -w, w): t3term})

            count += int(len(track) / len(freqs))
        return T3term

    def R4_dict(self, freqs, damp, kX, Nx, track, d_a_mo, nocc, norb):
        """
        Returns a dict with R[4] contractions

        """

        R4terms = {}
        count = 0

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

                Na = Nx['Na'][(track[i][0], w_s)]
                Nb = Nx['Nb'][(track[i][1], w1)]
                Nc = Nx['Nc'][(track[i][2], w2)]
                Nd = Nx['Nd'][(track[i][3], w3)]
                kA = kX['Na'][(track[i][0], w_s)]
                kB = kX['Nb'][(track[i][1], w1)]
                kC = kX['Nc'][(track[i][2], w2)]
                kD = kX['Nd'][(track[i][3], w3)]

                Nb_h = self.FlipXY(Nb)
                Nc_h = self.FlipXY(Nc)
                Nd_h = self.FlipXY(Nd)

                if track[i][0] in 'x':
                    R4x += -1j * damp * Nd_h @ self.S4nnn(
                        kA.T, kB, kC, d_a_mo, nocc, norb)
                    R4x += -1j * damp * Nc_h @ self.S4nnn(
                        kA.T, kB, kD, d_a_mo, nocc, norb)
                    R4x += -1j * damp * Nd_h @ self.S4nnn(
                        kA.T, kC, kB, d_a_mo, nocc, norb)
                    R4x += -1j * damp * Nb_h @ self.S4nnn(
                        kA.T, kC, kD, d_a_mo, nocc, norb)
                    R4x += -1j * damp * Nc_h @ self.S4nnn(
                        kA.T, kD, kB, d_a_mo, nocc, norb)
                    R4x += -1j * damp * Nb_h @ self.S4nnn(
                        kA.T, kD, kC, d_a_mo, nocc, norb)

                if track[i][0] in 'y':
                    R4y += -1j * damp * Nd_h @ self.S4nnn(
                        kA.T, kB, kC, d_a_mo, nocc, norb)
                    R4y += -1j * damp * Nc_h @ self.S4nnn(
                        kA.T, kB, kD, d_a_mo, nocc, norb)
                    R4y += -1j * damp * Nd_h @ self.S4nnn(
                        kA.T, kC, kB, d_a_mo, nocc, norb)
                    R4y += -1j * damp * Nb_h @ self.S4nnn(
                        kA.T, kC, kD, d_a_mo, nocc, norb)
                    R4y += -1j * damp * Nc_h @ self.S4nnn(
                        kA.T, kD, kB, d_a_mo, nocc, norb)
                    R4y += -1j * damp * Nb_h @ self.S4nnn(
                        kA.T, kD, kC, d_a_mo, nocc, norb)

                if track[i][0] in 'z':
                    R4z += -1j * damp * Nd_h @ self.S4nnn(
                        kA.T, kB, kC, d_a_mo, nocc, norb)
                    R4z += -1j * damp * Nc_h @ self.S4nnn(
                        kA.T, kB, kD, d_a_mo, nocc, norb)
                    R4z += -1j * damp * Nd_h @ self.S4nnn(
                        kA.T, kC, kB, d_a_mo, nocc, norb)
                    R4z += -1j * damp * Nb_h @ self.S4nnn(
                        kA.T, kC, kD, d_a_mo, nocc, norb)
                    R4z += -1j * damp * Nc_h @ self.S4nnn(
                        kA.T, kD, kB, d_a_mo, nocc, norb)
                    R4z += -1j * damp * Nb_h @ self.S4nnn(
                        kA.T, kD, kC, d_a_mo, nocc, norb)

            count += int(len(track) / len(freqs))
            R4terms.update({('x', w1): -R4x})
            R4terms.update({('y', w1): -R4y})
            R4terms.update({('z', w1): -R4z})

        return R4terms

    def AntiSym(self, Vec):
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

    def Nxy(self, eri_thresh, conv_thresh, lindep_thresh, max_iter, w, d_a_mo,
            damp, X, Fock_dict, kX, nocc, norb, molecule, ao_basis, scf_tensors,
            comm, ostream):
        """ 
        This part of the code creates all the XY[1] gradient vectors that are to be sent to the response solver to get the two-index response vectors

        : param keys - a list containing all the frequencies that are present for the different components
        : param w - a list of all the frequencies
        : param d_a_mo - the density matrix in MO basis
        : param damp - the damping parameter
        : param X - dipole integrals
        : param Fock_dict - A dictonary containing all the Fock matricies
        : param track - A list containg all the gamma tensor components that are to be computed
        : param kX - A dictonary containg all the response matricies
        : param nocc - The number of occupied orbitals
        : param norb - The number of total orbitals


        """

        N_total_drv = ComplexResponse(comm, ostream)

        XY_dict = {}
        test_dict = {}
        freq = None

        if self.rank == mpi_master():

            XY_dict = self.get_XY(d_a_mo, X, w, Fock_dict, kX, nocc, norb)

            #for l in XY_dict.keys():
            #    test_dict.update({l: XY_dict[l].real})

            wbd = [sum(x) for x in zip(w, w)]
            freq_bd_s = zip(w, wbd)
            ff = list(freq_bd_s)

            freq = '0.0,'
            for i in range(len(wbd)):
                freq += str(wbd[i]) + ','

        freq = self.comm.bcast(freq, root=mpi_master())

        N_total_drv.update_settings({
            'frequencies': freq,
            'damping': damp,
            'eri_thresh': eri_thresh,
            'conv_thresh': conv_thresh,
            'lindep_thresh': lindep_thresh,
            'max_iter': max_iter,
        })

        start = time.time()

        N_total_Drv = N_total_drv.compute(molecule, ao_basis, scf_tensors,
                                          XY_dict)

        end = time.time()

        if self.rank == mpi_master():
            Nxy_dict = N_total_Drv['solutions']
            kXY_dict = N_total_Drv['kappas']
            FXY_2_dict = N_total_Drv['Focks']

            return Nxy_dict, kXY_dict, FXY_2_dict, XY_dict
        else:
            return None, None, None, None

    def Nxy_red(self, eri_thresh, conv_thresh, lindep_thresh, max_iter, w,
                d_a_mo, damp, X, Fock_dict, kX, nocc, norb, molecule, ao_basis,
                scf_tensors, comm, ostream):
        """ 
        This part of the code creates all the XY[1] gradient vectors that are to be sent to the response solver to get the two-index response vectors

        : param keys - a list containing all the frequencies that are present for the different components
        : param w - a list of all the frequencies
        : param d_a_mo - the density matrix in MO basis
        : param damp - the damping parameter
        : param X - dipole integrals
        : param Fock_dict - A dictonary containing all the Fock matricies
        : param track - A list containg all the gamma tensor components that are to be computed
        : param kX - A dictonary containg all the response matricies
        : param nocc - The number of occupied orbitals
        : param norb - The number of total orbitals


        """
        XY_dict = {}
        test_dict = {}
        freq = None

        N_total_drv_2 = ComplexResponse(comm, ostream)

        if self.rank == mpi_master():

            XY_dict = self.get_XY_red(d_a_mo, X, w, Fock_dict, kX, nocc, norb)

            #for l in XY_dict.keys():
            #    test_dict.update({l: XY_dict[l].real})

            wbd = [sum(x) for x in zip(w, w)]
            freq_bd_s = zip(w, wbd)
            ff = list(freq_bd_s)

            freq = ''
            for i in range(len(wbd)):
                freq += str(wbd[i]) + ','

        freq = self.comm.bcast(freq, root=mpi_master())

        N_total_drv_2.update_settings({
            'frequencies': freq,
            'damping': damp,
            'eri_thresh': eri_thresh,
            'conv_thresh': conv_thresh,
            'lindep_thresh': lindep_thresh,
            'max_iter': max_iter,
        })

        start = time.time()

        N_total_Drv = N_total_drv_2.compute(molecule, ao_basis, scf_tensors,
                                            XY_dict)

        end = time.time()

        if self.rank == mpi_master():
            Nxy_dict = N_total_Drv['solutions']
            kXY_dict = N_total_Drv['kappas']
            FXY_2_dict = N_total_Drv['Focks']
            time_Nxy = end - start

            return Nxy_dict, kXY_dict, FXY_2_dict, XY_dict
        else:
            return None, None, None, None


################################# FOCK CODE ##################################

    def get_fock_r(self, D, molecule, ao_basis, rank):
        if rank == 0:
            da = D
            fa = self.get_two_el_fock_mod_r(
                molecule,
                ao_basis,
                da,
            )
            h = self.get_one_el_hamiltonian(molecule, ao_basis)

            if self.rank == mpi_master():
                fa = 0.5 * fa[0]
                fa += h
                return fa
            else:
                return None

        elif rank == 1:

            D_total = []
            if self.rank == mpi_master():
                for da in D:
                    D_total.append(da.real)
                    D_total.append(da.imag)
            else:
                D_total = None

            F_total = self.get_two_el_fock_mod_r(molecule, ao_basis, D_total)

            FF = []

            if self.rank == mpi_master():
                for i in range(int(0.5 * len(F_total))):
                    if i == 0:
                        FF.append(0.5 * F_total[0] + 0.5j * F_total[1])
                    else:
                        FF.append(0.5 * F_total[2 * i] +
                                  0.5j * F_total[2 * i + 1])
                return FF
            else:
                return None

        elif rank == 3:
            F_total = self.get_two_el_fock_mod_r(molecule, ao_basis, D)

            FF = []

            if self.rank == mpi_master():
                for i in range(len(F_total)):
                    FF.append(0.5j * F_total[i])

                return FF
            else:
                return None

        else:

            F_total = self.get_two_el_fock_mod_r(molecule, ao_basis, D)

            FF = []

            if self.rank == mpi_master():
                for i in range(len(F_total)):
                    FF.append(0.5 * F_total[i])

                return FF
            else:
                return None

    def get_two_el_fock_mod_r(self, molecule, ao_basis, *dabs):
        dts = []
        if self.rank == mpi_master():
            for dab in dabs[0]:
                dt = 2 * dab
                dts.append(dt)
            dens = AODensityMatrix(dts, denmat.rest)
        else:
            dens = AODensityMatrix()

        dens.broadcast(self.rank, self.comm)
        fock = AOFockMatrix(dens)

        for i in range(0, AOFockMatrix(dens).number_of_fock_matrices()):
            fock.set_fock_type(fockmat.rgenjk, i)
            fock.set_fock_type(fockmat.rgenk, i + 1)

        eri_driver = ElectronRepulsionIntegralsDriver(self.comm)
        screening = eri_driver.compute(ericut.qqden, 1.0e-15, molecule,
                                       ao_basis)
        eri_driver.compute(ericut.qqden, 1.0e-15, molecule, ao_basis)
        eri_driver.compute(fock, dens, molecule, ao_basis, screening)
        fock.reduce_sum(self.rank, self.size, self.comm)
        fabs = []
        if self.rank == mpi_master():
            for i in range(0, len(dabs[0])):
                ft = fock.to_numpy(i).T
                fabs.append(ft)

            return tuple(fabs)
        else:
            return None

    def get_one_el_hamiltonian(self, molecule, ao_basis):
        kinetic_driver = KineticEnergyIntegralsDriver(self.comm)
        potential_driver = NuclearPotentialIntegralsDriver(self.comm)

        T = kinetic_driver.compute(molecule, ao_basis).to_numpy()
        V = potential_driver.compute(molecule, ao_basis).to_numpy()

        return T - V

    def print_properties(self, iso, w1, gamma, comp, T4, T3, NaX3NyNz, NaA3NxNy,
                         NaX2Nyz, NxA2Nyz):
        """
        Prints properties.

        :param props:
            The dictionary of properties.
        """
        width = 50

        self.ostream.print_blank()
        self.ostream.print_blank()
        w_str = "Gamma tensor components computed per frequency"
        self.ostream.print_blank()
        self.ostream.print_header(w_str.ljust(width))
        self.ostream.print_blank()
        count = 1
        for a in range(int(len(comp) * (1 / len(w1)))):
            w_str = str(count) + '. ' + str(comp[a].split(",")[0])
            self.ostream.print_header(w_str.ljust(width))
            count += 1

        self.ostream.print_blank()

        for w in w1:
            if iso is False:
                w_str = "ΣT4term =  {:.8f}".format(T4[w, -w, w])
            else:
                w_str = "ΣT4term =  {:.8f}".format(T4[w, -w, w] / 15)
            self.ostream.print_header(w_str.ljust(width))

            if iso is False:
                w_str = "ΣT3term =  {:.8f}".format(T3[w, -w, w])
            else:
                w_str = "ΣT3term =  {:.8f}".format(T3[w, -w, w] / 15)
            self.ostream.print_header(w_str.ljust(width))

            if iso is False:
                w_str = "ΣNaX3NyNz =  {:.8f}".format(NaX3NyNz[w, -w, w])
            else:
                w_str = "ΣNaX3NyNz =  {:.8f}".format(NaX3NyNz[w, -w, w] / 15)

            self.ostream.print_header(w_str.ljust(width))

            if iso is False:
                w_str = "ΣNaA3NxNy =  {:.8f}".format(NaA3NxNy[w, -w, w])
            else:
                w_str = "ΣNaA3NxNy =  {:.8f}".format(NaA3NxNy[w, -w, w] / 15)
            self.ostream.print_header(w_str.ljust(width))

            if iso is False:
                w_str = "ΣNaX2Nyz =  {:.8f}".format(NaX2Nyz[w, -w, w])
            else:
                w_str = "ΣNaX2Nyz =  {:.8f}".format(NaX2Nyz[w, -w, w] / 15)
            self.ostream.print_header(w_str.ljust(width))

            if iso is False:
                w_str = "ΣNxA2Nyz =  {:.8f}".format(NxA2Nyz[w, -w, w])
            else:
                w_str = "ΣNxA2Nyz =  {:.8f}".format(NxA2Nyz[w, -w, w] / 15)
            self.ostream.print_header(w_str.ljust(width))

            w_str = 'Σ<<μ;μ,μ,μ>>= {:.8f}, ω=({:.4f},{:.4f},{:.4f}), hω =({:.4f} eV)'.format(
                gamma[w, -w, w], w, -w, w, w * 54.437358841503915 / 2)
            self.ostream.print_header(w_str.ljust(width))
            self.ostream.print_header(('-' * len(w_str)).ljust(width))
            self.ostream.print_blank()
