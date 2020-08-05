import numpy as np
import time
import re

from .veloxchemlib import mpi_master
from .cppsolver import ComplexResponse
from .linearsolver import LinearSolver
from .tpadriver import TpaDriver


class TpaReducedDriver(TpaDriver):
    """
    Implements the reduced isotropic cubic response driver for two-photon
    absorption (TPA)

    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.
    """

    def __init__(self, comm, ostream):
        """
        Initializes the reduced isotropic cubic response driver for two-photon
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

    def get_fock_dict(self, wi, density_list, F0_a, mo, molecule, ao_basis):
        """
        Computes the compounded Fock matrics F^{σ}  used for the reduced
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
                                    'real')
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

            for count, w in enumerate(wi):
                f_sig_xx[w] = fock_list[6 * count]
                f_sig_yy[w] = fock_list[6 * count + 1]
                f_sig_zz[w] = fock_list[6 * count + 2]
                f_sig_xy[w] = fock_list[6 * count + 3]
                f_sig_xz[w] = fock_list[6 * count + 4]
                f_sig_yz[w] = fock_list[6 * count + 5]

            Fock = {
                'F0': F0_a,
                'f_sig_xx': f_sig_xx,
                'f_sig_yy': f_sig_yy,
                'f_sig_zz': f_sig_zz,
                'f_sig_xy': f_sig_xy,
                'f_sig_xz': f_sig_xz,
                'f_sig_yz': f_sig_yz
            }

        return Fock

    def get_n_xy(self, w, d_a_mo, X, fock_dict, kX, nocc, norb, molecule,
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

        xy_dict = {}
        freq = None

        if self.rank == mpi_master():
            # Get the second-order gradiants
            xy_dict = self.get_xy(d_a_mo, X, w, fock_dict, kX, nocc, norb)

            wbd = [sum(x) for x in zip(w, w)]

            freq = ''
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
            N_total_drv.checkpoint_file += '_tpa_2_red.h5'

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

        return xy_dict

    def get_densities_II(self, wi, kX, kXY, S, D0, mo):
        """
        Computes the compounded densities needed for the compounded
        second-order Fock matrics used for the reduced isotropic cubic response
        function. Note: All densities are 1/3 of those in the paper, and all
        the Fock matrices are later scaled by 3.

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
                                    'real_and_imag')
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
        Contracts E[3]n_xNyz for the isotropic cubic response function. Takes
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
            ff_x = -2 * LinearSolver.lrmat2vec(X_terms, nocc, norb)
            ff_x = self.anti_sym(ff_x)
            f_iso_x[w] = ff_x

            # y #
            zeta_sig_yx = self.xi(k_x_, k_sig_xy, f_x_, f_sig_xy, F0_a)
            zeta_sig_yy = self.xi(k_y_, k_sig_yy, f_y_, f_sig_yy, F0_a)
            zeta_sig_yz = self.xi(k_z_, k_sig_yz, f_z_, f_sig_yz, F0_a)

            Y_terms = (zeta_sig_yx + zeta_sig_yy +
                       zeta_sig_yz).T + (0.5 * fo2['F123_y'][w]).T
            ff_y = -2 * LinearSolver.lrmat2vec(Y_terms, nocc, norb)
            ff_y = self.anti_sym(ff_y)
            f_iso_y[w] = ff_y

            # z #
            zeta_sig_zx = self.xi(k_x_, k_sig_xz, f_x_, f_sig_xz, F0_a)
            zeta_sig_zy = self.xi(k_y_, k_sig_yz, f_y_, f_sig_yz, F0_a)
            zeta_sig_zz = self.xi(k_z_, k_sig_zz, f_z_, f_sig_zz, F0_a)

            Z_terms = (zeta_sig_zx + zeta_sig_zy +
                       zeta_sig_zz).T + (0.5 * fo2['F123_z'][w]).T
            ff_z = -2 * LinearSolver.lrmat2vec(Z_terms, nocc, norb)
            ff_z = self.anti_sym(ff_z)
            f_iso_z[w] = ff_z

        return {'f_iso_x': f_iso_x, 'f_iso_y': f_iso_y, 'f_iso_z': f_iso_z}

    def get_other_terms(self, wi, track, n_x, n_xy, X, kX, kXY, da, nocc, norb):
        """
        Computes the terms involving X[2],A[2] in the reduced isotropic cubic
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

        for i in range(len(wi)):
            vals = track[i * comp_per_freq].split(',')
            w = float(vals[1])
            wa = float(vals[1])
            wb = float(vals[1])
            wc = float(vals[2])
            wd = float(vals[3])

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
            'NaX2Nyz': na_x2_nyz_dict,
            'NxA2Nyz': nx_a2_nyz_dict,
        }

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

            cont_label = "ΣNaX2Nyz  {:12.4f}".format(w)

            w_str = '{:<15s} {:15.8f} {:13.8f}j'.format(cont_label,
                                                        NaX2Nyz[w, -w, w].real,
                                                        NaX2Nyz[w, -w, w].imag)
            self.ostream.print_header(w_str.ljust(width))

            cont_label = "ΣNxA2Nyz  {:12.4f}".format(w)

            w_str = '{:<15s} {:15.8f} {:13.8f}j'.format(cont_label,
                                                        NxA2Nyz[w, -w, w].real,
                                                        NxA2Nyz[w, -w, w].imag)
            self.ostream.print_header(w_str.ljust(width))

            # TODO: print Gamma, and perhaps T3/X2/A2

            cont_label = "Σ<<μ;μ,μ,μ>>  {:8.4f}".format(w)

            w_str = '{:<15s} {:15.8f} {:13.8f}j'.format(cont_label,
                                                        gamma[w, -w, w].real,
                                                        gamma[w, -w, w].imag)
            self.ostream.print_header(w_str.ljust(width))
            self.ostream.print_blank()

        self.ostream.print_blank()
