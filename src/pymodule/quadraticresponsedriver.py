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
from .veloxchemlib import mpi_master, hartree_in_wavenumber
from .profiler import Profiler
from .outputstream import OutputStream
from .cppsolver import ComplexResponse
from .linearsolver import LinearSolver
from .nonlinearsolver import NonlinearSolver
from .distributedarray import DistributedArray
from .sanitychecks import (molecule_sanity_check, scf_results_sanity_check,
                           dft_sanity_check)
from .errorhandler import assert_msg_critical
from .checkpoint import (check_distributed_focks, read_distributed_focks,
                         write_distributed_focks)


class QuadraticResponseDriver(NonlinearSolver):
    """
    Implements a general quadratic response driver.

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
        - a_component: Cartesian component of the A operator.
        - b_component: Cartesian component of the B operator.
        - c_component: Cartesian component of the C operator.
    """

    def __init__(self, comm=None, ostream=None):
        """
        Initializes the quadratic response driver.
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
        self.comp = None
        self.damping = 1000.0 / hartree_in_wavenumber()

        self.a_component = 'z'
        self.b_component = 'z'
        self.c_component = 'z'

        # input keywords
        self._input_keywords['response'].update({
            'b_frequencies': ('seq_range', 'B frequencies'),
            'c_frequencies': ('seq_range', 'C frequencies'),
            'damping': ('float', 'damping parameter'),
            'a_component': ('str_lower', 'Cartesian component of A operator'),
            'b_component': ('str_lower', 'Cartesian component of B operator'),
            'c_component': ('str_lower', 'Cartesian component of C operator'),
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
        Computes a quadratic response function.

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

        # check molecule
        molecule_sanity_check(molecule)

        # check SCF results
        scf_results_sanity_check(self, scf_tensors)

        # check dft setup
        dft_sanity_check(self, 'compute', 'nonlinear')

        profiler = Profiler({
            'timing': False,
            'profiling': self.profiling,
            'memory_profiling': self.memory_profiling,
            'memory_tracing': self.memory_tracing,
        })

        if self.rank == mpi_master():
            self._print_header('Quadratic Response Driver Setup')

        start_time = time.time()

        # sanity check
        nalpha = molecule.number_of_alpha_electrons()
        nbeta = molecule.number_of_beta_electrons()
        assert_msg_critical(
            nalpha == nbeta,
            'QuadaticResponseDriver: not implemented for unrestricted case')

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
        a_grad = linear_solver.get_complex_prop_grad(operator, self.a_component,
                                                     molecule, ao_basis,
                                                     scf_tensors)
        b_grad = linear_solver.get_complex_prop_grad(operator, self.b_component,
                                                     molecule, ao_basis,
                                                     scf_tensors)
        c_grad = linear_solver.get_complex_prop_grad(operator, self.c_component,
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

        # Storing the dipole integral matrices used for the X[2] and
        # A[2] contractions in MO basis
        wa = [sum(x) for x in zip(self.b_frequencies, self.c_frequencies)]

        freqpairs = [wl for wl in zip(self.b_frequencies, self.c_frequencies)]

        ABC = {}

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

            ABC.update(A)
            ABC.update(B)
            ABC.update(C)

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

        cpp_keywords = [
            'damping', 'norm_thresh', 'lindep_thresh', 'conv_thresh',
            'max_iter', 'eri_thresh', 'qq_type', 'timing', 'memory_profiling',
            'batch_size', 'restart', 'xcfun', 'grid_level', 'potfile',
            'electric_field', 'program_end_time'
        ]

        for key in cpp_keywords:
            setattr(N_drv, key, getattr(self, key))

        if self.checkpoint_file is not None:
            N_drv.checkpoint_file = str(
                Path(self.checkpoint_file).with_suffix('.qrf.h5'))

        N_results = N_drv.compute(molecule, ao_basis, scf_tensors, ABC)

        self._is_converged = N_drv.is_converged

        Nx = N_results['solutions']
        Focks = N_results['focks']

        profiler.check_memory_usage('CPP')

        quad_dict = self.compute_quad_components(Focks, freqpairs, X, d_a_mo,
                                                 Nx, self.comp, scf_tensors,
                                                 molecule, ao_basis, profiler)

        valstr = '*** Time spent in quadratic response calculation: '
        valstr += '{:.2f} sec ***'.format(time.time() - start_time)
        self.ostream.print_header(valstr)
        self.ostream.print_blank()
        self.ostream.flush()

        profiler.end(self.ostream)

        return quad_dict

    def compute_quad_components(self, Focks, freqpairs, X, d_a_mo, Nx, track,
                                scf_tensors, molecule, ao_basis, profiler):
        """
        Computes all the relevent terms to compute a general quadratic response function

        :param freqparis:
            A list of all the frequencies
        :param X:
            A dictonary of matricies containing all the property integrals
        :param d_a_mo:
            The SCF density in MO basis
        :param Nx:
            A dictonary containing all the response vectors in distributed form
        :param scf_tensors:
            The dictionary of tensors from converged SCF wavefunction.
        :param molecule:
            The molecule.
        :param basis:
            The AO basis.
        :param profiler:
            The profiler.

        :return:
            A dictionary containing all the relevent terms for quadratic response
        """

        if self.rank == mpi_master():
            mo = scf_tensors['C_alpha']
            F0 = np.linalg.multi_dot([mo.T, scf_tensors['F_alpha'], mo])
            norb = mo.shape[1]
        else:
            mo = None
            F0 = None
            norb = None
        F0 = self.comm.bcast(F0, root=mpi_master())
        norb = self.comm.bcast(norb, root=mpi_master())

        nocc = molecule.number_of_alpha_electrons()

        dft_dict = self._init_dft(molecule, scf_tensors)

        # computing all compounded first-order densities
        first_order_dens, second_order_dens = self.get_densities(
            freqpairs, Nx, mo, nocc, norb)

        profiler.check_memory_usage('Densities')

        #  computing the compounded first-order Fock matrices
        fock_dict = self.get_fock_dict(freqpairs, first_order_dens,
                                       second_order_dens, F0, mo, molecule,
                                       ao_basis, dft_dict)

        profiler.check_memory_usage('Focks')

        e3_dict = self.get_e3(freqpairs, Nx, fock_dict, Focks, nocc, norb)

        profiler.check_memory_usage('E[3]')

        result = {}

        for (wb, wc) in freqpairs:

            Na = ComplexResponse.get_full_solution_vector(Nx[('A', (wb + wc))])
            Nb = ComplexResponse.get_full_solution_vector(Nx[('B', wb)])
            Nc = ComplexResponse.get_full_solution_vector(Nx[('C', wc)])

            if self.rank == mpi_master():

                op_a = X[self.a_component]
                op_b = X[self.b_component]
                op_c = X[self.c_component]

                kb = self.complex_lrvec2mat(Nb, nocc, norb)
                kc = self.complex_lrvec2mat(Nc, nocc, norb)

                C2Nb = self._x2_contract(kb, op_c, d_a_mo, nocc, norb)
                B2Nc = self._x2_contract(kc, op_b, d_a_mo, nocc, norb)

                A2Nc = self._a2_contract(kc, op_a, d_a_mo, nocc, norb)
                A2Nb = self._a2_contract(kb, op_a, d_a_mo, nocc, norb)

                NaE3NbNc = np.dot(Na.T, e3_dict[wb])
                NaC2Nb = np.dot(Na.T, C2Nb)
                NaB2Nc = np.dot(Na.T, B2Nc)
                NbA2Nc = np.dot(Nb.T, A2Nc)
                NcA2Nb = np.dot(Nc.T, A2Nb)

                val_X2 = -(NaC2Nb + NaB2Nc)
                val_A2 = -(NbA2Nc + NcA2Nb)
                val_E3 = NaE3NbNc
                beta = val_E3 + val_A2 + val_X2

                self.ostream.print_blank()
                w_str = 'Quadratic response function: '
                w_str += '<< {};{},{} >>  ({},{})'.format(
                    self.a_component, self.b_component, self.c_component,
                    str(wb), str(wc))
                self.ostream.print_header(w_str)
                self.ostream.print_header('=' * (len(w_str) + 2))
                self.ostream.print_blank()

                title = '{:<9s} {:>20s} {:>21s}'.format('Component', 'Real',
                                                        'Imaginary')
                width = len(title)
                self.ostream.print_header(title.ljust(width))
                self.ostream.print_header(('-' * len(title)).ljust(width))
                self._print_component('X2', val_X2, width)
                self._print_component('A2', val_A2, width)
                self._print_component('E3', val_E3, width)
                self._print_component('beta', beta, width)
                self.ostream.print_blank()
                self.ostream.flush()

                result[(wb, wc)] = beta

        profiler.check_memory_usage('End of QRF')

        return result

    def get_densities(self, freqpairs, Nx, mo, nocc, norb):
        """
        Computes the densities needed for the perturbed Fock matrices.

        :param freqpairs:
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
            first_order_dens (first-order one-time tranformed compounded
            densities) and second_order_dens (first-order two-time tranformed
            compounded densities)
        """

        distributed_density_1 = None
        distributed_density_2 = None

        for (wb, wc) in freqpairs:

            Nb = ComplexResponse.get_full_solution_vector(Nx[('B', wb)])
            Nc = ComplexResponse.get_full_solution_vector(Nx[('C', wc)])

            if self.rank == mpi_master():

                kb = self.complex_lrvec2mat(Nb, nocc, norb)
                kc = self.complex_lrvec2mat(Nc, nocc, norb)

                # create the first order single indexed densiteies #

                Db = self.commut_mo_density(kb, nocc)
                Dc = self.commut_mo_density(kc, nocc)

                # create the first order two indexed densities #

                Dbc = self.commut(kb, Dc) + self.commut(kc, Db)

                # density transformation from MO to AO basis

                Db = np.linalg.multi_dot([mo, Db, mo.T])
                Dc = np.linalg.multi_dot([mo, Dc, mo.T])
                Dbc = np.linalg.multi_dot([mo, Dbc, mo.T])

                dist_den_1_freq = np.hstack((
                    Db.real.reshape(-1, 1),
                    Db.imag.reshape(-1, 1),
                    Dc.real.reshape(-1, 1),
                    Dc.imag.reshape(-1, 1),
                ))

                dist_den_2_freq = np.hstack((
                    Dbc.real.reshape(-1, 1),
                    Dbc.imag.reshape(-1, 1),
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

    def get_fock_dict(self,
                      wi,
                      first_order_dens,
                      second_order_dens,
                      F0,
                      mo,
                      molecule,
                      ao_basis,
                      dft_dict=None):
        """
        Computes the Fock matrices for a quadratic response function

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

        for (wb, wc) in wi:
            for key in ['FbcFcb']:
                key_freq_pairs.append((key, wb))

        # examine checkpoint file for distributed Focks

        if self.checkpoint_file is not None:
            fock_file = str(
                Path(self.checkpoint_file).with_suffix('.qrf_fock.h5'))
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

            dist_focks = self._comp_nlr_fock(mo, molecule, ao_basis,
                                             'real_and_imag', dft_dict,
                                             first_order_dens,
                                             second_order_dens, None, 'qrf')

            self._print_fock_time(time.time() - time_start_fock)

            write_distributed_focks(fock_file, dist_focks, key_freq_pairs,
                                    self.comm, self.ostream)

        focks = {'F0': F0}

        for fock_index, (key, wb) in enumerate(key_freq_pairs):
            if key not in focks:
                focks[key] = {}
            focks[key][wb] = DistributedArray(dist_focks.data[:, fock_index],
                                              self.comm,
                                              distribute=False)

        return focks

    def get_e3(self, wi, Nx, fo, fo2, nocc, norb):
        """
        Contracts E[3]

        :param wi:
            A list of freqs
        :param Nx:
            A dict of the single index response vectors in distributed form
        :param fo:
            A dictonary of transformed Fock matricies from fock_dict
        :param fo2:
            A dictonarty of transfromed Fock matricies from subspace of
            response solver
        :param nocc:
            The number of occupied orbitals
        :param norb:
            The total number of orbitals

        :return:
            A dictonary of compounded E[3] tensors for the isotropic cubic
            response function for QRF
        """

        e3vec = {}

        for (wb, wc) in wi:

            vec_pack = np.array([
                fo['FbcFcb'][wb].data,
                fo2[('B', wb)].data,
                fo2[('C', wc)].data,
            ]).T.copy()

            vec_pack = self._collect_vectors_in_columns(vec_pack)

            Nb = ComplexResponse.get_full_solution_vector(Nx[('B', wb)])
            Nc = ComplexResponse.get_full_solution_vector(Nx[('C', wc)])

            if self.rank != mpi_master():
                continue

            vec_pack = vec_pack.T.copy().reshape(-1, norb, norb)

            (fbcfcb, fb, fc) = vec_pack

            fb = np.conjugate(fb).T
            fc = np.conjugate(fc).T

            F0_a = fo['F0']

            # Response

            kb = (self.complex_lrvec2mat(Nb, nocc, norb)).T
            kc = (self.complex_lrvec2mat(Nc, nocc, norb)).T

            xi = self._xi(kb, kc, fb, fc, F0_a)

            e3fock = xi.T + 0.5 * fbcfcb.T

            e3vec[wb] = self.anti_sym(
                -2 * LinearSolver.lrmat2vec(e3fock, nocc, norb))

        return e3vec

    def _print_component(self, label, value, width):
        """
        Prints QRF component.

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
