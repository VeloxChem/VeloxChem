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

from .veloxchemlib import ElectricDipoleIntegralsDriver
from .veloxchemlib import mpi_master, hartree_in_wavenumbers
from .profiler import Profiler
from .outputstream import OutputStream
from .cppsolver import ComplexResponse
from .linearsolver import LinearSolver
from .nonlinearsolver import NonlinearSolver
from .distributedarray import DistributedArray
from .errorhandler import assert_msg_critical
from .checkpoint import (check_distributed_focks, read_distributed_focks,
                         write_distributed_focks)


class CubicResponseDriver(NonlinearSolver):
    """
    Implements a general cubic response driver

    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.

    Instance variables
        - is_converged: The flag for convergence.
        - comp: The list of all the gamma tensor components
        - damping: The damping parameter.
        - conv_thresh: The convergence threshold for the solver.
        - max_iter: The maximum number of solver iterations.
    """

    def __init__(self, comm=None, ostream=None):
        """
        Initializes the cubic response driver.
        """

        if comm is None:
            comm = MPI.COMM_WORLD

        if ostream is None:
            if comm.Get_rank() == mpi_master():
                ostream = OutputStream(sys.stdout)
            else:
                ostream = OutputStream(None)

        super().__init__(comm, ostream)

        # cpp settings
        self.b_frequencies = (0,)
        self.c_frequencies = (0,)
        self.d_frequencies = (0,)
        self.comp = None
        self.damping = 1000.0 / hartree_in_wavenumbers()

        self.a_components = 'z'
        self.b_components = 'z'
        self.c_components = 'z'
        self.d_components = 'z'

        # input keywords
        self._input_keywords['response'].update({
            'b_frequencies': ('seq_range', 'B frequencies'),
            'c_frequencies': ('seq_range', 'C frequencies'),
            'd_frequencies': ('seq_range', 'D frequencies'),
            'damping': ('float', 'damping parameter'),
            'a_components': ('str_lower', 'Cartesian components of A operator'),
            'b_components': ('str_lower', 'Cartesian components of B operator'),
            'c_components': ('str_lower', 'Cartesian components of C operator'),
            'd_components': ('str_lower', 'Cartesian components of D operator'),
        })

    def update_settings(self, rsp_dict, method_dict=None):
        """
        Updates response and method settings

        :param rsp_dict:
            The dictionary of response input.
        :param method_dict:
            The dictionary of method settings.
        """

        if method_dict is None:
            method_dict = {}

        super().update_settings(rsp_dict, method_dict)

    def compute(self, molecule, ao_basis, scf_tensors):
        """
        Computes a cubic response function.

        :param molecule:
            The molecule.
        :param basis:
            The AO basis.
        :param scf_tensors:
            The dictionary of tensors from converged SCF wavefunction.

        :return:
              A dictonary containing the E[3], X[2], A[2] contractions
        """

        if self.norm_thresh is None:
            self.norm_thresh = self.conv_thresh * 1.0e-6
        if self.lindep_thresh is None:
            self.lindep_thresh = self.conv_thresh * 1.0e-6

        # double check SCF information
        self._check_scf_results(scf_tensors)

        # check dft setup
        self._dft_sanity_check_nonlinrsp()

        profiler = Profiler({
            'timing': False,
            'profiling': self.profiling,
            'memory_profiling': self.memory_profiling,
            'memory_tracing': self.memory_tracing,
        })

        if self.rank == mpi_master():
            self._print_header('Cubic Response Driver Setup')

        start_time = time.time()

        # sanity check
        nalpha = molecule.number_of_alpha_electrons()
        nbeta = molecule.number_of_beta_electrons()
        assert_msg_critical(
            nalpha == nbeta,
            'CubicResponseDriver: not implemented for unrestricted case')

        if self.rank == mpi_master():
            S = scf_tensors['S']
            da = scf_tensors['D_alpha']
            mo = scf_tensors['C_alpha']
            d_a_mo = np.linalg.multi_dot([mo.T, S, da, S, mo])
            norb = mo.shape[1]
        else:
            d_a_mo = None
            norb = None
        d_a_mo = self.comm.bcast(d_a_mo, root=mpi_master())
        norb = self.comm.bcast(norb, root=mpi_master())

        # Computing first-order gradient vectors
        dipole_drv = ElectricDipoleIntegralsDriver(self.comm)
        dipole_mats = dipole_drv.compute(molecule, ao_basis)

        operator = 'dipole'
        linear_solver = LinearSolver(self.comm, self.ostream)
        a_grad = linear_solver.get_complex_prop_grad(operator,
                                                     self.a_components,
                                                     molecule, ao_basis,
                                                     scf_tensors)
        b_grad = linear_solver.get_complex_prop_grad(operator,
                                                     self.b_components,
                                                     molecule, ao_basis,
                                                     scf_tensors)
        c_grad = linear_solver.get_complex_prop_grad(operator,
                                                     self.c_components,
                                                     molecule, ao_basis,
                                                     scf_tensors)
        d_grad = linear_solver.get_complex_prop_grad(operator,
                                                     self.d_components,
                                                     molecule, ao_basis,
                                                     scf_tensors)

        if self.rank == mpi_master():
            inv_sqrt_2 = 1.0 / np.sqrt(2.0)

            a_grad = list(a_grad)
            for ind in range(len(a_grad)):
                a_grad[ind] *= inv_sqrt_2

            b_grad = list(b_grad)
            for ind in range(len(b_grad)):
                b_grad[ind] *= inv_sqrt_2

            c_grad = list(c_grad)
            for ind in range(len(c_grad)):
                c_grad[ind] *= inv_sqrt_2

            d_grad = list(d_grad)
            for ind in range(len(d_grad)):
                d_grad[ind] *= inv_sqrt_2

        # Storing the dipole integral matrices used for the X[3],X[2],A[3] and
        # A[2] contractions in MO basis
        wa = [
            sum(x) for x in zip(
                self.b_frequencies,
                self.c_frequencies,
                self.d_frequencies,
            )
        ]

        freqtriples = [
            wl for wl in zip(
                self.b_frequencies,
                self.c_frequencies,
                self.d_frequencies,
            )
        ]

        ABCD = {}

        if self.rank == mpi_master():
            A = {(op, w): v for op, v in zip('A', a_grad) for w in wa}
            B = {
                (op, w): v for op, v in zip('B', b_grad)
                for w in self.b_frequencies
            }
            C = {
                (op, w): v for op, v in zip('C', c_grad)
                for w in self.c_frequencies
            }
            D = {
                (op, w): v for op, v in zip('D', d_grad)
                for w in self.d_frequencies
            }

            ABCD.update(A)
            ABCD.update(B)
            ABCD.update(C)
            ABCD.update(D)

            X = {
                'x': 2 * self.ao2mo(mo, dipole_mats.x_to_numpy()),
                'y': 2 * self.ao2mo(mo, dipole_mats.y_to_numpy()),
                'z': 2 * self.ao2mo(mo, dipole_mats.z_to_numpy())
            }
        else:
            X = None
            self.comp = None

        # Computing the first-order response vectors (3 per frequency)
        N_drv = ComplexResponse(self.comm, self.ostream)

        cpp_keywords = {
            'damping', 'norm_thresh', 'lindep_thresh', 'conv_thresh',
            'max_iter', 'eri_thresh', 'qq_type', 'timing', 'memory_profiling',
            'batch_size', 'restart', 'xcfun', 'grid_level', 'potfile',
            'electric_field', 'program_end_time'
        }

        for key in cpp_keywords:
            setattr(N_drv, key, getattr(self, key))

        if self.checkpoint_file is not None:
            N_drv.checkpoint_file = str(
                Path(self.checkpoint_file).with_suffix('.crf_1.h5'))

        N_results = N_drv.compute(molecule, ao_basis, scf_tensors, ABCD)

        self._is_converged = N_drv.is_converged

        kX = N_results['kappas']
        Focks = N_results['focks']

        profiler.check_memory_usage('1st CPP')

        cubic_dict = self.compute_cubic_components(Focks, freqtriples, X,
                                                   d_a_mo, kX, self.comp,
                                                   scf_tensors, molecule,
                                                   ao_basis, profiler)

        valstr = '*** Time spent in cubic response calculation: '
        valstr += '{:.2f} sec ***'.format(time.time() - start_time)
        self.ostream.print_header(valstr)
        self.ostream.print_blank()
        self.ostream.flush()

        profiler.end(self.ostream)

        return cubic_dict

    def compute_cubic_components(self, Focks, freqtriples, X, d_a_mo, kX, track,
                                 scf_tensors, molecule, ao_basis, profiler):
        """
        Computes all the relevent terms to compute a general cubic response function

        :param w:
            A list of all the frequencies
        :param X:
            A dictonary of matricies containing all the dipole integrals
        :param d_a_mo:
            The SCF density in MO basis
        :param kX:
            A dictonary containing all the response matricies
        :param scf_tensors:
            The dictionary of tensors from converged SCF wavefunction.
        :param molecule:
            The molecule.
        :param basis:
            The AO basis.
        :param profiler:
            The profiler.

        :return:

        """

        if self.rank == mpi_master():
            mo = scf_tensors['C_alpha']
            F0 = np.linalg.multi_dot([mo.T, scf_tensors['F_alpha'], mo])
            norb = mo.shape[1]
        else:
            mo = None
            F0 = None
            norb = None

        dft_dict = self._init_dft(molecule, scf_tensors)

        F0 = self.comm.bcast(F0, root=mpi_master())
        norb = self.comm.bcast(norb, root=mpi_master())

        nocc = molecule.number_of_alpha_electrons()

        # computing all compounded first-order densities
        if self.rank == mpi_master():
            density_list1, density_list2, density_list3 = self.get_densities(
                freqtriples, kX, mo, nocc)
        else:
            density_list1 = None
            density_list2 = None
            density_list3 = None

        profiler.check_memory_usage('1st densities')

        #  computing the compounded first-order Fock matrices
        fock_dict = self.get_fock_dict(freqtriples, density_list1,
                                       density_list2, density_list3, F0, mo,
                                       molecule, ao_basis, dft_dict)

        profiler.check_memory_usage('1st Focks')

        e4_dict, s4_dict, r4_dict = self.get_esr4(freqtriples, kX, fock_dict,
                                                  Focks, nocc, norb, d_a_mo)

        profiler.check_memory_usage('E[4]')

        k_xy, f_xy = self.get_nxy(freqtriples, kX, fock_dict, Focks, nocc, norb,
                                  d_a_mo, X, molecule, ao_basis, scf_tensors)

        profiler.check_memory_usage('2nd CPP')

        if self.rank == mpi_master():
            density_list1_two, density_list2_two = self.get_densities_II(
                freqtriples, kX, k_xy, mo, nocc)
        else:
            density_list1_two = None
            density_list2_two = None

        profiler.check_memory_usage('2nd densities')

        fock_dict_two = self.get_fock_dict_II(freqtriples, density_list1_two,
                                              density_list2_two, F0, mo,
                                              molecule, ao_basis, dft_dict)

        profiler.check_memory_usage('2nd Focks')

        e3_dict = self.get_e3(freqtriples, kX, k_xy, Focks, fock_dict_two, f_xy,
                              nocc, norb)

        profiler.check_memory_usage('E[3]')

        result = {}

        if self.rank == mpi_master():

            A = X[self.a_components]
            B = X[self.b_components]
            C = X[self.c_components]
            D = X[self.d_components]

            for (wb, wc, wd) in freqtriples:

                Na = self.complex_lrmat2vec(kX[('A', wb + wc + wd)], nocc, norb)
                Nb = self.complex_lrmat2vec(kX[('B', wb)], nocc, norb)
                Nc = self.complex_lrmat2vec(kX[('C', wc)], nocc, norb)
                Nd = self.complex_lrmat2vec(kX[('D', wd)], nocc, norb)

                Nbd = self.complex_lrmat2vec(k_xy[('BD', wb, wd), wb + wd],
                                             nocc, norb)
                Nbc = self.complex_lrmat2vec(k_xy[('BC', wb, wc), wb + wc],
                                             nocc, norb)
                Ncd = self.complex_lrmat2vec(k_xy[('CD', wc, wd), wc + wd],
                                             nocc, norb)

                NaE4NbNcNd = np.dot(Na, e4_dict[wb])
                NaS4NbNcNd = np.dot(Na, s4_dict[wb])
                NaR4NbNcNd = r4_dict[wb]

                NaE3NbNcd = np.dot(Na, e3_dict[(('E3'), (wb, wc, wd))])

                # X3 terms
                NaB3NcNd = np.dot(
                    Na.T,
                    self._x3_contract(kX[('C', wc)], kX[('D', wd)], B, d_a_mo,
                                      nocc, norb))
                NaB3NdNc = np.dot(
                    Na.T,
                    self._x3_contract(kX[('D', wd)], kX[('C', wc)], B, d_a_mo,
                                      nocc, norb))

                NaC3NbNd = np.dot(
                    Na.T,
                    self._x3_contract(kX[('B', wb)], kX[('D', wd)], C, d_a_mo,
                                      nocc, norb))
                NaC3NdNb = np.dot(
                    Na.T,
                    self._x3_contract(kX[('D', wd)], kX[('B', wb)], C, d_a_mo,
                                      nocc, norb))

                NaD3NbNc = np.dot(
                    Na.T,
                    self._x3_contract(kX[('B', wb)], kX[('C', wc)], D, d_a_mo,
                                      nocc, norb))
                NaD3NcNb = np.dot(
                    Na.T,
                    self._x3_contract(kX[('C', wc)], kX[('B', wb)], D, d_a_mo,
                                      nocc, norb))

                # X2 contraction
                NaB2Ncd = np.dot(
                    Na.T,
                    self._x2_contract(k_xy[('CD', wc, wd), wc + wd], B, d_a_mo,
                                      nocc, norb))
                NaC2Nbd = np.dot(
                    Na.T,
                    self._x2_contract(k_xy[('BD', wb, wd), wb + wd], C, d_a_mo,
                                      nocc, norb))
                NaD2Nbc = np.dot(
                    Na.T,
                    self._x2_contract(k_xy[('BC', wb, wc), wb + wc], D, d_a_mo,
                                      nocc, norb))

                # A3 contraction
                NdA3NbNc = np.dot(
                    self._a3_contract(kX[('B', wb)], kX[('C', wc)], A, d_a_mo,
                                      nocc, norb), Nd)
                NdA3NcNb = np.dot(
                    self._a3_contract(kX[('C', wc)], kX[('B', wb)], A, d_a_mo,
                                      nocc, norb), Nd)

                NbA3NcNd = np.dot(
                    self._a3_contract(kX[('C', wc)], kX[('D', wd)], A, d_a_mo,
                                      nocc, norb), Nb)
                NbA3NdNc = np.dot(
                    self._a3_contract(kX[('D', wd)], kX[('C', wc)], A, d_a_mo,
                                      nocc, norb), Nb)

                NcA3NbNd = np.dot(
                    self._a3_contract(kX[('B', wb)], kX[('D', wd)], A, d_a_mo,
                                      nocc, norb), Nc)
                NcA3NdNb = np.dot(
                    self._a3_contract(kX[('D', wd)], kX[('B', wb)], A, d_a_mo,
                                      nocc, norb), Nc)

                # A2 contraction
                NbA2Ncd = np.dot(
                    self._a2_contract(kX[('B', wb)], A, d_a_mo, nocc, norb),
                    Ncd)
                NcdA2Nb = np.dot(
                    self._a2_contract(k_xy[('CD', wc, wd), wc + wd], A, d_a_mo,
                                      nocc, norb), Nb)

                NcA2Nbd = np.dot(
                    self._a2_contract(kX[('C', wc)], A, d_a_mo, nocc, norb),
                    Nbd)
                NbdA2Nc = np.dot(
                    self._a2_contract(k_xy[('BD', wb, wd), wb + wd], A, d_a_mo,
                                      nocc, norb), Nc)

                NdA2Nbc = np.dot(
                    self._a2_contract(kX[('D', wd)], A, d_a_mo, nocc, norb),
                    Nbc)
                NbcA2Nd = np.dot(
                    self._a2_contract(k_xy[('BC', wb, wc), wb + wc], A, d_a_mo,
                                      nocc, norb), Nd)

                val_E3 = -(NaE3NbNcd)
                val_T4 = -(NaE4NbNcNd - NaS4NbNcNd - NaR4NbNcNd)
                val_X2 = NaB2Ncd + NaC2Nbd + NaD2Nbc
                val_X3 = NaB3NcNd + NaB3NdNc + NaC3NbNd + NaC3NdNb + NaD3NbNc + NaD3NcNb
                val_A2 = NbA2Ncd + NcdA2Nb + NcA2Nbd + NbdA2Nc + NdA2Nbc + NbcA2Nd
                val_A3 = -(NdA3NbNc + NdA3NcNb + NbA3NcNd + NbA3NdNc +
                           NcA3NbNd + NcA3NdNb)

                # Cubic response function
                gamma = val_T4 + val_E3 + val_X3 + val_A3 + val_X2 + val_A2

                self.ostream.print_blank()
                w_str = 'Cubic response function: << {};{},{},{} >>  ({},{},{})'.format(
                    self.a_components, self.b_components, self.c_components,
                    self.d_components, str(wb), str(wc), str(wd))
                self.ostream.print_header(w_str)
                self.ostream.print_header('=' * (len(w_str) + 2))
                self.ostream.print_blank()

                title = '{:<9s} {:>20s} {:>21s}'.format('Component', 'Real',
                                                        'Imaginary')
                width = len(title)
                self.ostream.print_header(title.ljust(width))
                self.ostream.print_header(('-' * len(title)).ljust(width))
                self._print_component('E3', val_E3, width)
                self._print_component('T4', val_T4, width)
                self._print_component('X2', val_X2, width)
                self._print_component('X3', val_X3, width)
                self._print_component('A2', val_A2, width)
                self._print_component('A3', val_A3, width)
                self._print_component('gamma', gamma, width)
                self.ostream.print_blank()

                result[('E3', wb, wc, wd)] = val_E3
                result[('T4', wb, wc, wd)] = val_T4
                result[('X3', wb, wc, wd)] = val_X3
                result[('X2', wb, wc, wd)] = val_X2
                result[('A3', wb, wc, wd)] = val_A3
                result[('A2', wb, wc, wd)] = val_A2
                result[('gamma', wb, wc, wd)] = gamma

        profiler.check_memory_usage('End of CRF')

        return result

    def get_e3(self, wi, kX, k_xy, fo, fo2, fo3, nocc, norb):
        """
        Contracts E[3] for CRF

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
            response function for CRF
        """

        e3vec = {}

        for (wb, wc, wd) in wi:

            vec_pack = np.array([
                fo[('B', wb)].data,
                fo[('C', wc)].data,
                fo[('D', wd)].data,
                fo2['F123'][(wb, wc, wd)].data,
                fo3[(('BC', wb, wc), wb + wc)].data,
                fo3[(('BD', wb, wd), wb + wd)].data,
                fo3[(('CD', wc, wd), wc + wd)].data,
            ]).T.copy()

            vec_pack = self._collect_vectors_in_columns(vec_pack)

            if self.rank != mpi_master():
                continue

            vec_pack = vec_pack.T.copy().reshape(-1, norb, norb)

            (fb, fc, fd, f123, fbc, fbd, fcd) = vec_pack

            fb = np.conjugate(fb).T
            fc = np.conjugate(fc).T
            fd = np.conjugate(fd).T
            fbc = np.conjugate(fbc).T
            fcd = np.conjugate(fcd).T
            fbd = np.conjugate(fbd).T

            F0_a = fo2['F0']

            # E3NbNcd

            kb = kX[('B', wb)].T
            kcd = k_xy[('CD', wc, wd), wc + wd].T

            xi_b_cd = self._xi(kb, kcd, fb, fcd, F0_a)

            # E3NcNbd

            kc = kX[('C', wc)].T
            kbd = k_xy[('BD', wb, wd), wb + wd].T

            xi_c_bd = self._xi(kc, kbd, fc, fbd, F0_a)

            # E3NdNbc

            kd = kX[('D', wd)].T
            kbc = k_xy[('BC', wb, wc), wb + wc].T

            xi_d_bc = self._xi(kd, kbc, fd, fbc, F0_a)

            e3fock = (xi_b_cd + xi_c_bd + xi_d_bc).T + (0.5 * f123).T

            e3vec[('E3', (wb, wc, wd))] = self.anti_sym(
                -2 * LinearSolver.lrmat2vec(e3fock, nocc, norb))

        return e3vec

    def get_densities(self, freqtriples, kX, mo, nocc):
        """
        Computes the  densities needed for the Fock matrices.

        :param freqtriples:
            A list of the frequency triples
        :param kX:
            A dictonary with all the first-order response matrices
        :param mo:
            A matrix containing the MO coefficents
        :param nocc:
            Number of occupied orbitals

        :return:
            A list of tranformed compounded densities
        """

        density_list1 = []
        density_list2 = []
        density_list3 = []

        for (wb, wc, wd) in freqtriples:

            # convert response matrix to ao basis #

            kb = kX[('B', wb)]
            kc = kX[('C', wc)]
            kd = kX[('D', wd)]

            # create the first order single indexed densiteies #

            Db = self.commut_mo_density(kb, nocc)
            Dc = self.commut_mo_density(kc, nocc)
            Dd = self.commut_mo_density(kd, nocc)

            # create the first order two indexed densities #

            Dbc = self.commut(kb, Dc)
            Dcb = self.commut(kc, Db)

            Dbd = self.commut(kb, Dd)
            Ddb = self.commut(kd, Db)

            Ddc = self.commut(kd, Dc)
            Dcd = self.commut(kc, Dd)

            # create the first order three indexed densities #

            Dbcd = self.commut(kb, Dcd)
            Dbdc = self.commut(kb, Ddc)

            Dcbd = self.commut(kc, Dbd)
            Dcdb = self.commut(kc, Ddb)

            Ddbc = self.commut(kd, Dbc)
            Ddcb = self.commut(kd, Dcb)

            # density transformation from MO to AO basis

            Db = np.linalg.multi_dot([mo, Db, mo.T])
            Dc = np.linalg.multi_dot([mo, Dc, mo.T])
            Dd = np.linalg.multi_dot([mo, Dd, mo.T])

            Dbc = np.linalg.multi_dot([mo, (Dbc + Dcb), mo.T])
            Dbd = np.linalg.multi_dot([mo, (Dbd + Ddb), mo.T])
            Dcd = np.linalg.multi_dot([mo, (Ddc + Dcd), mo.T])

            D123 = np.linalg.multi_dot(
                [mo, (Dbcd + Dbdc + Dcbd + Dcdb + Ddbc + Ddcb), mo.T])

            density_list1.append(Db.real)
            density_list1.append(Db.imag)
            density_list1.append(Dc.real)
            density_list1.append(Dc.imag)
            density_list1.append(Dd.real)
            density_list1.append(Dd.imag)

            density_list2.append(Dbc.real)
            density_list2.append(Dbc.imag)
            density_list2.append(Dbd.real)
            density_list2.append(Dbd.imag)
            density_list2.append(Dcd.real)
            density_list2.append(Dcd.imag)

            density_list3.append(D123.real)
            density_list3.append(D123.imag)

        return density_list1, density_list2, density_list3

    def get_densities_II(self, freqtriples, kX, k_xy, mo, nocc):
        """
        Computes the  densities needed for the Fock matrices.

        :param freqtriples:
            A list of the frequency triples
        :param kX:
            A dictonary with all the first-order response matrices
        :param k_xy:
            A dict of the two index response matrices
        :param mo:
            A matrix containing the MO coefficents
        :param nocc:
            Number of occupied orbitals

        :return:
            A list of tranformed compounded densities
        """

        density_list1 = []
        density_list2 = []

        for (wb, wc, wd) in freqtriples:

            # convert response matrix to ao basis #

            kb = kX[('B', wb)]
            kc = kX[('C', wc)]
            kd = kX[('D', wd)]

            kbc = k_xy[(('BC', wb, wc), wb + wc)]
            kbd = k_xy[(('BD', wb, wd), wb + wd)]
            kcd = k_xy[(('CD', wc, wd), wc + wd)]

            # create the first order single indexed densiteies #

            Db = self.commut_mo_density(kb, nocc)
            Dc = self.commut_mo_density(kc, nocc)
            Dd = self.commut_mo_density(kd, nocc)

            # create the second-order two indexed densities #

            Dbc = self.commut_mo_density(kbc, nocc)
            Dbd = self.commut_mo_density(kbd, nocc)
            Dcd = self.commut_mo_density(kcd, nocc)

            # create the second-order three indexed densities #

            Db_cd = self.commut(kb, Dcd)
            Dcd_b = self.commut(kcd, Db)

            Dc_bd = self.commut(kc, Dbd)
            Dbd_c = self.commut(kbd, Dc)

            Dd_bc = self.commut(kd, Dbc)
            Dbc_d = self.commut(kbc, Dd)

            # density transformation from MO to AO basis

            d_total = np.linalg.multi_dot(
                [mo, (Db_cd + Dcd_b + Dc_bd + Dbd_c + Dd_bc + Dbc_d), mo.T])

            Db = np.linalg.multi_dot([mo, Db, mo.T])
            Dc = np.linalg.multi_dot([mo, Dc, mo.T])
            Dd = np.linalg.multi_dot([mo, Dd, mo.T])
            Dbc = np.linalg.multi_dot([mo, Dbc, mo.T])
            Dbd = np.linalg.multi_dot([mo, Dbd, mo.T])
            Dcd = np.linalg.multi_dot([mo, Dcd, mo.T])

            density_list1.append(Db.real)
            density_list1.append(Db.imag)
            density_list1.append(Dc.real)
            density_list1.append(Dc.imag)
            density_list1.append(Dd.real)
            density_list1.append(Dd.imag)

            density_list1.append(Dbc.real)
            density_list1.append(Dbc.imag)
            density_list1.append(Dbd.real)
            density_list1.append(Dbd.imag)
            density_list1.append(Dcd.real)
            density_list1.append(Dcd.imag)

            density_list2.append(d_total.real)
            density_list2.append(d_total.imag)

        return density_list1, density_list2

    def get_fock_dict(self, wi, density_list1, density_list2, density_list3, F0,
                      mo, molecule, ao_basis, dft_dict):
        """
        Computes the Fock matrices for a cubic response function

        :param wi:
            A list of the frequencies
        :param density_list:
            A list of tranformed compounded densities
        :param F0:
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

        for (wb, wc, wd) in wi:
            for key in ['Fbc', 'Fbd', 'Fcd']:
                key_freq_pairs.append((key, wb))

        for (wb, wc, wd) in wi:
            for key in ['Fbcd']:
                key_freq_pairs.append((key, wb))

        # examine checkpoint file for distributed Focks

        if self.checkpoint_file is not None:
            fock_file = str(
                Path(self.checkpoint_file).with_suffix('.crf_fock_1.h5'))
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
                                                 density_list3, 'crf')
            else:
                if self.rank == mpi_master():
                    density_list_23 = density_list2 + density_list3
                else:
                    density_list_23 = None
                dist_focks = self._comp_nlr_fock(mo, molecule, ao_basis,
                                                 'real_and_imag', None, None,
                                                 None, density_list_23, 'crf')

            self._print_fock_time(time.time() - time_start_fock)

            write_distributed_focks(fock_file, dist_focks, key_freq_pairs,
                                    self.comm, self.ostream)

        # assign distributed Focks to key-frequency pairs

        focks = {'F0': F0}

        for fock_index, (key, wb) in enumerate(key_freq_pairs):
            if key not in focks:
                focks[key] = {}
            focks[key][wb] = DistributedArray(dist_focks.data[:, fock_index],
                                              self.comm,
                                              distribute=False)

        return focks

    def get_fock_dict_II(self, wi, density_list1, density_list2, F0, mo,
                         molecule, ao_basis, dft_dict):
        """
        Computes the Fock matrices for a cubic response function

        :param wi:
            A list of the frequencies
        :param density_list:
            A list of tranformed compounded densities
        :param F0:
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

        for (wb, wc, wd) in wi:
            for key in ['F123']:
                key_freq_pairs.append((key, (wb, wc, wd)))

        # examine checkpoint file for distributed Focks

        if self.checkpoint_file is not None:
            fock_file = str(
                Path(self.checkpoint_file).with_suffix('.crf_fock_2.h5'))
        else:
            fock_file = None

        if self.restart:
            if self.rank == mpi_master():
                self.restart = check_distributed_focks(fock_file,
                                                       key_freq_pairs)
            self.restart = self.comm.bcast(self.restart, mpi_master())

        # examine checkpoint file for distributed Focks

        if self.restart:
            dist_focks = read_distributed_focks(fock_file, self.comm,
                                                self.ostream)
        else:
            time_start_fock = time.time()

            if self._dft:
                dist_focks = self._comp_nlr_fock(mo, molecule, ao_basis,
                                                 'real_and_imag', dft_dict,
                                                 density_list1, density_list2,
                                                 None, 'crf_ii')
            else:
                dist_focks = self._comp_nlr_fock(mo, molecule, ao_basis,
                                                 'real_and_imag', None, None,
                                                 density_list2, None, 'crf_ii')

            self._print_fock_time(time.time() - time_start_fock)

            write_distributed_focks(fock_file, dist_focks, key_freq_pairs,
                                    self.comm, self.ostream)

        # assign distributed Focks to key-frequency pairs

        focks = {'F0': F0}

        for fock_index, (key, (wb, wc, wd)) in enumerate(key_freq_pairs):
            if key not in focks:
                focks[key] = {}
            focks[key][(wb, wc,
                        wd)] = DistributedArray(dist_focks.data[:, fock_index],
                                                self.comm,
                                                distribute=False)

        return focks

    def get_esr4(self, wi, kX, fo, fo2, nocc, norb, D0):
        """
        Contracts E[4], S[4], R[4]

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
            A dictonary of E[4], S[4], R[4] tensor contractions.
        """

        e4_vec = {}
        s4_vec = {}
        r4_vec = {}

        for (wb, wc, wd) in wi:

            vec_pack = np.array([
                fo2[('B', wb)].data,
                fo2[('C', wc)].data,
                fo2[('D', wd)].data,
                fo['Fbc'][wb].data,
                fo['Fbd'][wb].data,
                fo['Fcd'][wb].data,
                fo['Fbcd'][wb].data,
            ]).T.copy()

            vec_pack = self._collect_vectors_in_columns(vec_pack)

            if self.rank != mpi_master():
                continue

            vec_pack = vec_pack.T.copy().reshape(-1, norb, norb)

            (fb, fc, fd, fbc, fbd, fcd, fbcd) = vec_pack

            fb = np.conjugate(fb).T
            fc = np.conjugate(fc).T
            fd = np.conjugate(fd).T

            F0_a = fo['F0']

            # Response

            kb = kX[('B', wb)].T
            kc = kX[('C', wc)].T
            kd = kX[('D', wd)].T

            zi_bcd = self._zi(kb, kc, kd, fc, fd, fcd, F0_a)
            zi_cbd = self._zi(kc, kb, kd, fb, fd, fbd, F0_a)
            zi_dbc = self._zi(kd, kb, kc, fb, fc, fbc, F0_a)

            e4fock = (zi_bcd + zi_cbd + zi_dbc) + (fbcd)

            e4vec = 2. / 6 * self.anti_sym(
                LinearSolver.lrmat2vec(e4fock.T, nocc, norb))

            ka = kX[('A', (wb + wc + wd))]
            kb = kX[('B', wb)]
            kc = kX[('C', wc)]
            kd = kX[('D', wd)]

            s4_term = wb * self._s4(kb, kc, kd, D0, nocc, norb)
            s4_term += wc * self._s4(kc, kb, kd, D0, nocc, norb)
            s4_term += wd * self._s4(kd, kb, kc, D0, nocc, norb)

            Nb = self.complex_lrmat2vec(kb, nocc, norb)
            Nc = self.complex_lrmat2vec(kc, nocc, norb)
            Nd = self.complex_lrmat2vec(kd, nocc, norb)

            Nb_h = self.flip_xy(Nb)
            Nc_h = self.flip_xy(Nc)
            Nd_h = self.flip_xy(Nd)

            r4_term = 1j * self.damping * np.dot(
                Nd_h, self._s4_for_r4(ka.T, kb, kc, D0, nocc, norb))
            r4_term += 1j * self.damping * np.dot(
                Nc_h, self._s4_for_r4(ka.T, kb, kd, D0, nocc, norb))
            r4_term += 1j * self.damping * np.dot(
                Nd_h, self._s4_for_r4(ka.T, kc, kb, D0, nocc, norb))
            r4_term += 1j * self.damping * np.dot(
                Nb_h, self._s4_for_r4(ka.T, kc, kd, D0, nocc, norb))
            r4_term += 1j * self.damping * np.dot(
                Nc_h, self._s4_for_r4(ka.T, kd, kb, D0, nocc, norb))
            r4_term += 1j * self.damping * np.dot(
                Nb_h, self._s4_for_r4(ka.T, kd, kc, D0, nocc, norb))

            e4_vec[wb] = e4vec
            s4_vec[wb] = s4_term
            r4_vec[wb] = r4_term

        return e4_vec, s4_vec, r4_vec

    def get_nxy(self, wi, kX, fo, fo2, nocc, norb, d_a_mo, X, molecule,
                ao_basis, scf_tensors):
        """
        Computed NXY

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
            A dictonary of E[4], S[4], R[4] tensor contractions.
        """

        BC = {}
        CD = {}
        BD = {}

        XY = {}

        for (wb, wc, wd) in wi:

            vec_pack = np.array([
                fo2[('B', wb)].data,
                fo2[('C', wc)].data,
                fo2[('D', wd)].data,
                fo['Fbc'][wb].data,
                fo['Fbd'][wb].data,
                fo['Fcd'][wb].data,
            ]).T.copy()

            vec_pack = self._collect_vectors_in_columns(vec_pack)

            if self.rank != mpi_master():
                continue

            vec_pack = vec_pack.T.copy().reshape(-1, norb, norb)

            (fb, fc, fd, fbc, fbd, fcd) = vec_pack

            fb = np.conjugate(fb).T
            fc = np.conjugate(fc).T
            fd = np.conjugate(fd).T

            F0_a = fo['F0']

            # Response

            kb = kX[('B', wb)].T
            kc = kX[('C', wc)].T
            kd = kX[('D', wd)].T

            B = X[self.b_components]
            C = X[self.c_components]
            D = X[self.d_components]

            # BC

            xi = self._xi(kb, kc, fb, fc, F0_a)

            e3fock = xi.T + (0.5 * fbc).T
            E3NbNc = self.anti_sym(-LinearSolver.lrmat2vec(e3fock, nocc, norb))

            C2Nb = 0.5 * self._x2_contract(kX[('B', wb)], C, d_a_mo, nocc, norb)
            B2Nc = 0.5 * self._x2_contract(kX[('C', wc)], B, d_a_mo, nocc, norb)

            BC[(('BC', wb, wc), wb + wc)] = E3NbNc - C2Nb - B2Nc

            # BD

            xi = self._xi(kb, kd, fb, fd, F0_a)

            e3fock = xi.T + (0.5 * fbd).T
            E3NbNd = self.anti_sym(-LinearSolver.lrmat2vec(e3fock, nocc, norb))

            D2Nb = 0.5 * self._x2_contract(kX[('B', wb)], D, d_a_mo, nocc, norb)
            B2Nd = 0.5 * self._x2_contract(kX[('D', wd)], B, d_a_mo, nocc, norb)

            BD[(('BD', wb, wd), wb + wd)] = E3NbNd - D2Nb - B2Nd

            # CD

            xi = self._xi(kc, kd, fc, fd, F0_a)

            e3fock = xi.T + (0.5 * fcd).T
            E3NcNd = self.anti_sym(-LinearSolver.lrmat2vec(e3fock, nocc, norb))

            C2Nd = 0.5 * self._x2_contract(kX[('D', wd)], C, d_a_mo, nocc, norb)
            D2Nc = 0.5 * self._x2_contract(kX[('C', wc)], D, d_a_mo, nocc, norb)

            CD[(('CD', wc, wd), wc + wd)] = E3NcNd - C2Nd - D2Nc

            XY.update(BC)
            XY.update(BD)
            XY.update(CD)

        Nxy_drv = ComplexResponse(self.comm, self.ostream)

        cpp_keywords = {
            'damping', 'norm_thresh', 'lindep_thresh', 'conv_thresh',
            'max_iter', 'eri_thresh', 'qq_type', 'timing', 'memory_profiling',
            'batch_size', 'restart', 'xcfun', 'grid_level', 'potfile',
            'electric_field', 'program_end_time'
        }

        for key in cpp_keywords:
            setattr(Nxy_drv, key, getattr(self, key))

        if self.checkpoint_file is not None:
            Nxy_drv.checkpoint_file = str(
                Path(self.checkpoint_file).with_suffix('.crf_2.h5'))

        Nxy_results = Nxy_drv.compute(molecule, ao_basis, scf_tensors, XY)

        self._is_converged = (self._is_converged and Nxy_drv.is_converged)

        kX = Nxy_results['kappas']
        Focks = Nxy_results['focks']

        return kX, Focks

    def _print_component(self, label, value, width):
        """
        Prints response function components.

        :param label:
            The label
        :param freq:
            The frequency
        :param value:
            The complex value
        :param width:
            The width for the output
        """

        w_str = '{:<9s} {:20.8f} {:20.8f}j'.format(label, value.real,
                                                   value.imag)
        self.ostream.print_header(w_str.ljust(width))
