#
#                                   VELOXCHEM
#              ----------------------------------------------------
#                          An Electronic Structure Code
#
#  SPDX-License-Identifier: BSD-3-Clause
#
#  Copyright 2018-2025 VeloxChem developers
#
#  Redistribution and use in source and binary forms, with or without modification,
#  are permitted provided that the following conditions are met:
#
#  1. Redistributions of source code must retain the above copyright notice, this
#     list of conditions and the following disclaimer.
#  2. Redistributions in binary form must reproduce the above copyright notice,
#     this list of conditions and the following disclaimer in the documentation
#     and/or other materials provided with the distribution.
#  3. Neither the name of the copyright holder nor the names of its contributors
#     may be used to endorse or promote products derived from this software without
#     specific prior written permission.
#
#  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
#  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
#  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
#  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
#  FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
#  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
#  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
#  HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
#  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
#  OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

from mpi4py import MPI
from pathlib import Path
import numpy as np
import time
import sys

from .oneeints import compute_electric_dipole_integrals
from .oneeints import compute_quadrupole_integrals
from .oneeints import compute_linear_momentum_integrals
from .oneeints import compute_angular_momentum_integrals
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
        self.damping = 1000.0 / hartree_in_wavenumber()

        self.a_operator = 'electric dipole'
        self.b_operator = 'electric dipole'
        self.c_operator = 'electric dipole'
        self.d_operator = 'electric dipole'

        self.a_component = None
        self.b_component = None
        self.c_component = None
        self.d_component = None

        # input keywords
        self._input_keywords['response'].update({
            'b_frequencies': ('seq_range', 'B frequencies'),
            'c_frequencies': ('seq_range', 'C frequencies'),
            'd_frequencies': ('seq_range', 'D frequencies'),
            'damping': ('float', 'damping parameter'),
            'a_operator': ('str_lower', 'A operator'),
            'b_operator': ('str_lower', 'B operator'),
            'c_operator': ('str_lower', 'C operator'),
            'd_operator': ('str_lower', 'D operator'),
            'a_component': ('str_lower', 'Cartesian component of A operator'),
            'b_component': ('str_lower', 'Cartesian component of B operator'),
            'c_component': ('str_lower', 'Cartesian component of C operator'),
            'd_component': ('str_lower', 'Cartesian component of D operator'),
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

        operator_mapping = {
            'dipole': 'electric_dipole',
            'electric dipole': 'electric_dipole',
            'quadrupole': 'electric_quadrupole',
            'electric quadrupole': 'electric_quadrupole',
            'linear momentum': 'linear_momentum',
            'angular momentum': 'angular_momentum',
            'magnetic dipole': 'magnetic_dipole',
        }

        if self.a_operator in operator_mapping:
            self._a_op_key = operator_mapping[self.a_operator]

        if self.b_operator in operator_mapping:
            self._b_op_key = operator_mapping[self.b_operator]

        if self.c_operator in operator_mapping:
            self._c_op_key = operator_mapping[self.c_operator]

        if self.d_operator in operator_mapping:
            self._d_op_key = operator_mapping[self.d_operator]

        # for backward compatibility
        if self.a_component is None and hasattr(self, 'a_components'):
            self.a_component = self.a_components

        if self.b_component is None and hasattr(self, 'b_components'):
            self.b_component = self.b_components

        if self.c_component is None and hasattr(self, 'c_components'):
            self.c_component = self.c_components

        if self.d_component is None and hasattr(self, 'd_components'):
            self.d_component = self.d_components

        # sanity check
        assert_msg_critical(
            self.a_component in ['x', 'y', 'z'],
            'CubicResponseDriver: Undefined or invalid a_component')

        assert_msg_critical(
            self.b_component in ['x', 'y', 'z'],
            'CubicResponseDriver: Undefined or invalid b_component')

        assert_msg_critical(
            self.c_component in ['x', 'y', 'z'],
            'CubicResponseDriver: Undefined or invalid c_component')

        assert_msg_critical(
            self.d_component in ['x', 'y', 'z'],
            'CubicResponseDriver: Undefined or invalid d_component')

        if self.norm_thresh is None:
            self.norm_thresh = self.conv_thresh * 1.0e-6
        if self.lindep_thresh is None:
            self.lindep_thresh = self.conv_thresh * 1.0e-6

        # check molecule
        molecule_sanity_check(molecule)

        # check SCF results
        scf_results_sanity_check(self, scf_tensors)

        # update checkpoint_file after scf_results_sanity_check
        if self.filename is not None and self.checkpoint_file is None:
            self.checkpoint_file = f'{self.filename}_rsp.h5'

        # check dft setup
        dft_sanity_check(self, 'compute', 'nonlinear')

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
        if self.rank == mpi_master():
            dipole_mats = compute_electric_dipole_integrals(
                molecule, ao_basis, [0.0, 0.0, 0.0])
            quadrupole_mats = compute_quadrupole_integrals(
                molecule, ao_basis, [0.0, 0.0, 0.0])
            linmom_mats = compute_linear_momentum_integrals(
                molecule, ao_basis)
            angmom_mats = compute_angular_momentum_integrals(
                molecule, ao_basis, [0.0, 0.0, 0.0])
        else:
            dipole_mats = None
            quadrupole_mats = None
            linmom_mats = None
            angmom_mats = None

        linear_solver = LinearSolver(self.comm, self.ostream)
        a_grad = linear_solver.get_complex_prop_grad(self._a_op_key, self.a_component,
                                                     molecule, ao_basis,
                                                     scf_tensors)
        b_grad = linear_solver.get_complex_prop_grad(self._b_op_key, self.b_component,
                                                     molecule, ao_basis,
                                                     scf_tensors)
        c_grad = linear_solver.get_complex_prop_grad(self._c_op_key, self.c_component,
                                                     molecule, ao_basis,
                                                     scf_tensors)
        d_grad = linear_solver.get_complex_prop_grad(self._d_op_key, self.d_component,
                                                     molecule, ao_basis,
                                                     scf_tensors)

        if self.rank == mpi_master():
            inv_sqrt_2 = 1.0 / np.sqrt(2.0)

            a_grad = list(a_grad)
            for ind in range(len(a_grad)):
                a_grad[ind] *= inv_sqrt_2
                # Note: nonliear response uses r instead of mu for dipole operator
                if self._a_op_key == 'electric_dipole':
                    a_grad[ind] *= -1.0

            b_grad = list(b_grad)
            for ind in range(len(b_grad)):
                b_grad[ind] *= inv_sqrt_2
                # Note: nonliear response uses r instead of mu for dipole operator
                if self._b_op_key == 'electric_dipole':
                    b_grad[ind] *= -1.0

            c_grad = list(c_grad)
            for ind in range(len(c_grad)):
                c_grad[ind] *= inv_sqrt_2
                # Note: nonliear response uses r instead of mu for dipole operator
                if self._c_op_key == 'electric_dipole':
                    c_grad[ind] *= -1.0

            d_grad = list(d_grad)
            for ind in range(len(d_grad)):
                d_grad[ind] *= inv_sqrt_2
                # Note: nonliear response uses r instead of mu for dipole operator
                if self._d_op_key == 'electric_dipole':
                    d_grad[ind] *= -1.0

        # Storing the dipole integral matrices used for the X[3],X[2],A[3] and
        # A[2] contractions in MO basis
        # note: use plain addition instead of sum
        wa = [
            x[0] + x[1] + x[2] for x in zip(
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

            quadrupole_r2 = (quadrupole_mats[0] + quadrupole_mats[3] +
                             quadrupole_mats[5])

            X = {
                # Note: nonliear response uses r instead of mu for dipole operator
                'electric_dipole': {
                    'x': 2 * self.ao2mo(mo, dipole_mats[0]) * (-1.0),
                    'y': 2 * self.ao2mo(mo, dipole_mats[1]) * (-1.0),
                    'z': 2 * self.ao2mo(mo, dipole_mats[2]) * (-1.0),
                },
                'electric_quadrupole': {
                    'xx': 2 * self.ao2mo(
                        mo, quadrupole_mats[0] * 3.0 - quadrupole_r2) * (-1.0),
                    'xy': 2 * self.ao2mo(mo, quadrupole_mats[1]) * (-1.0),
                    'xz': 2 * self.ao2mo(mo, quadrupole_mats[2]) * (-1.0),
                    'yy': 2 * self.ao2mo(
                        mo, quadrupole_mats[3] * 3.0 - quadrupole_r2) * (-1.0),
                    'yz': 2 * self.ao2mo(mo, quadrupole_mats[4]) * (-1.0),
                    'zz': 2 * self.ao2mo(
                        mo, quadrupole_mats[5] * 3.0 - quadrupole_r2) * (-1.0),
                },
                'linear_momentum': {
                    'x': 2 * self.ao2mo(mo, linmom_mats[0]) * (-1j),
                    'y': 2 * self.ao2mo(mo, linmom_mats[1]) * (-1j),
                    'z': 2 * self.ao2mo(mo, linmom_mats[2]) * (-1j),
                },
                'angular_momentum': {
                    'x': 2 * self.ao2mo(mo, angmom_mats[0]) * (-1j),
                    'y': 2 * self.ao2mo(mo, angmom_mats[1]) * (-1j),
                    'z': 2 * self.ao2mo(mo, angmom_mats[2]) * (-1j),
                },
                'magnetic_dipole': {
                    'x': 2 * self.ao2mo(mo, angmom_mats[0]) * (0.5j),
                    'y': 2 * self.ao2mo(mo, angmom_mats[1]) * (0.5j),
                    'z': 2 * self.ao2mo(mo, angmom_mats[2]) * (0.5j),
                },

            }
        else:
            X = None
            self.comp = None

        # Computing the first-order response vectors (3 per frequency)
        N_drv = ComplexResponse(self.comm, self.ostream)

        cpp_keywords = {
            'damping', 'norm_thresh', 'lindep_thresh', 'conv_thresh',
            'max_iter', 'eri_thresh', 'timing', 'memory_profiling',
            'batch_size', 'restart', 'xcfun', 'grid_level', 'potfile',
            'electric_field', 'program_end_time', '_debug', '_block_size_factor',
            'ri_coulomb'
        }

        for key in cpp_keywords:
            setattr(N_drv, key, getattr(self, key))

        if self.checkpoint_file is not None:
            fpath = Path(self.checkpoint_file)
            fpath = fpath.with_name(fpath.stem)
            N_drv.checkpoint_file = str(fpath) + '_crf_1.h5'

        N_results = N_drv.compute(molecule, ao_basis, scf_tensors, ABCD)

        self._is_converged = N_drv.is_converged

        Nx = N_results['solutions']
        Focks = N_results['focks']

        profiler.check_memory_usage('1st CPP')

        cubic_dict = self.compute_cubic_components(Focks, freqtriples, X,
                                                   d_a_mo, Nx, self.comp,
                                                   scf_tensors, molecule,
                                                   ao_basis, profiler)

        valstr = '*** Time spent in cubic response calculation: '
        valstr += '{:.2f} sec ***'.format(time.time() - start_time)
        self.ostream.print_header(valstr)
        self.ostream.print_blank()
        self.ostream.flush()

        profiler.end(self.ostream)

        return cubic_dict

    def compute_cubic_components(self, Focks, freqtriples, X, d_a_mo, Nx, track,
                                 scf_tensors, molecule, ao_basis, profiler):
        """
        Computes all the relevent terms to compute a general cubic response function

        :param w:
            A list of all the frequencies
        :param X:
            A dictonary of matricies containing all the dipole integrals
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

        eri_dict = self._init_eri(molecule, ao_basis)

        dft_dict = self._init_dft(molecule, scf_tensors)

        # computing all compounded first-order densities
        density_list1, density_list2, density_list3 = self.get_densities(
            freqtriples, Nx, mo, nocc, norb)

        profiler.check_memory_usage('1st densities')

        #  computing the compounded first-order Fock matrices
        fock_dict = self.get_fock_dict(freqtriples, density_list1,
                                       density_list2, density_list3, F0, mo,
                                       molecule, ao_basis, eri_dict, dft_dict)

        profiler.check_memory_usage('1st Focks')

        e4_dict, s4_dict, r4_dict = self.get_esr4(freqtriples, Nx, fock_dict,
                                                  Focks, nocc, norb, d_a_mo)

        profiler.check_memory_usage('E[4]')

        Nxy, f_xy = self.get_nxy(freqtriples, Nx, fock_dict, Focks, nocc, norb,
                                 d_a_mo, X, molecule, ao_basis, scf_tensors)

        profiler.check_memory_usage('2nd CPP')

        density_list1_two, density_list2_two = self.get_densities_II(
            freqtriples, Nx, Nxy, mo, nocc, norb)

        profiler.check_memory_usage('2nd densities')

        fock_dict_two = self.get_fock_dict_II(freqtriples, density_list1_two,
                                              density_list2_two, F0, mo,
                                              molecule, ao_basis, eri_dict,
                                              dft_dict)

        profiler.check_memory_usage('2nd Focks')

        e3_dict = self.get_e3(freqtriples, Nx, Nxy, Focks, fock_dict_two, f_xy,
                              nocc, norb)

        profiler.check_memory_usage('E[3]')

        result = {}

        for (wb, wc, wd) in freqtriples:

            Na = ComplexResponse.get_full_solution_vector(Nx[('A',
                                                              (wb + wc + wd))])
            Nb = ComplexResponse.get_full_solution_vector(Nx[('B', wb)])
            Nc = ComplexResponse.get_full_solution_vector(Nx[('C', wc)])
            Nd = ComplexResponse.get_full_solution_vector(Nx[('D', wd)])

            Nbc = ComplexResponse.get_full_solution_vector(Nxy[(('BC', wb, wc),
                                                                wb + wc)])
            Nbd = ComplexResponse.get_full_solution_vector(Nxy[(('BD', wb, wd),
                                                                wb + wd)])
            Ncd = ComplexResponse.get_full_solution_vector(Nxy[(('CD', wc, wd),
                                                                wc + wd)])

            if self.rank == mpi_master():

                A = X[self._a_op_key][self.a_component]
                B = X[self._b_op_key][self.b_component]
                C = X[self._c_op_key][self.c_component]
                D = X[self._d_op_key][self.d_component]

                kb = self.complex_lrvec2mat(Nb, nocc, norb)
                kc = self.complex_lrvec2mat(Nc, nocc, norb)
                kd = self.complex_lrvec2mat(Nd, nocc, norb)

                kbc = self.complex_lrvec2mat(Nbc, nocc, norb)
                kbd = self.complex_lrvec2mat(Nbd, nocc, norb)
                kcd = self.complex_lrvec2mat(Ncd, nocc, norb)

                NaE4NbNcNd = np.dot(Na, e4_dict[wb])
                NaS4NbNcNd = np.dot(Na, s4_dict[wb])
                NaR4NbNcNd = r4_dict[wb]

                NaE3NbNcd = np.dot(Na, e3_dict[('E3', (wb, wc, wd))])

                # X3 terms
                NaB3NcNd = np.dot(
                    Na.T, self._x3_contract(kc, kd, B, d_a_mo, nocc, norb))
                NaB3NdNc = np.dot(
                    Na.T, self._x3_contract(kd, kc, B, d_a_mo, nocc, norb))

                NaC3NbNd = np.dot(
                    Na.T, self._x3_contract(kb, kd, C, d_a_mo, nocc, norb))
                NaC3NdNb = np.dot(
                    Na.T, self._x3_contract(kd, kb, C, d_a_mo, nocc, norb))

                NaD3NbNc = np.dot(
                    Na.T, self._x3_contract(kb, kc, D, d_a_mo, nocc, norb))
                NaD3NcNb = np.dot(
                    Na.T, self._x3_contract(kc, kb, D, d_a_mo, nocc, norb))

                # X2 contraction
                NaB2Ncd = np.dot(Na.T,
                                 self._x2_contract(kcd, B, d_a_mo, nocc, norb))
                NaC2Nbd = np.dot(Na.T,
                                 self._x2_contract(kbd, C, d_a_mo, nocc, norb))
                NaD2Nbc = np.dot(Na.T,
                                 self._x2_contract(kbc, D, d_a_mo, nocc, norb))

                # A3 contraction
                NdA3NbNc = np.dot(
                    self._a3_contract(kb, kc, A, d_a_mo, nocc, norb), Nd)
                NdA3NcNb = np.dot(
                    self._a3_contract(kc, kb, A, d_a_mo, nocc, norb), Nd)

                NbA3NcNd = np.dot(
                    self._a3_contract(kc, kd, A, d_a_mo, nocc, norb), Nb)
                NbA3NdNc = np.dot(
                    self._a3_contract(kd, kc, A, d_a_mo, nocc, norb), Nb)

                NcA3NbNd = np.dot(
                    self._a3_contract(kb, kd, A, d_a_mo, nocc, norb), Nc)
                NcA3NdNb = np.dot(
                    self._a3_contract(kd, kb, A, d_a_mo, nocc, norb), Nc)

                # A2 contraction
                NbA2Ncd = np.dot(self._a2_contract(kb, A, d_a_mo, nocc, norb),
                                 Ncd)
                NcdA2Nb = np.dot(self._a2_contract(kcd, A, d_a_mo, nocc, norb),
                                 Nb)

                NcA2Nbd = np.dot(self._a2_contract(kc, A, d_a_mo, nocc, norb),
                                 Nbd)
                NbdA2Nc = np.dot(self._a2_contract(kbd, A, d_a_mo, nocc, norb),
                                 Nc)

                NdA2Nbc = np.dot(self._a2_contract(kd, A, d_a_mo, nocc, norb),
                                 Nbc)
                NbcA2Nd = np.dot(self._a2_contract(kbc, A, d_a_mo, nocc, norb),
                                 Nd)

                op_a_type = 'imag' if self.is_imag(self._a_op_key) else 'real'
                op_b_type = 'imag' if self.is_imag(self._b_op_key) else 'real'
                op_c_type = 'imag' if self.is_imag(self._c_op_key) else 'real'
                op_d_type = 'imag' if self.is_imag(self._d_op_key) else 'real'

                # flip sign for X3 and X2 terms
                if op_a_type == 'real':
                    if (op_b_type == op_c_type) and (op_c_type != op_d_type):
                        NaD3NbNc *= -1.0
                        NaD3NcNb *= -1.0
                        NaD2Nbc *= -1.0
                    elif (op_c_type == op_d_type) and (op_d_type != op_b_type):
                        NaB3NcNd *= -1.0
                        NaB3NdNc *= -1.0
                        NaB2Ncd *= -1.0
                    elif (op_d_type == op_b_type) and (op_b_type != op_c_type):
                        NaC3NbNd *= -1.0
                        NaC3NdNb *= -1.0
                        NaC2Nbd *= -1.0
                elif op_a_type == 'imag':
                    if (op_b_type == op_c_type) and (op_c_type != op_d_type):
                        if op_d_type == 'real':
                            NaD3NbNc *= -1.0
                            NaD3NcNb *= -1.0
                            NaD2Nbc *= -1.0
                        else:
                            NaB3NcNd *= -1.0
                            NaB3NdNc *= -1.0
                            NaB2Ncd *= -1.0
                            NaC3NbNd *= -1.0
                            NaC3NdNb *= -1.0
                            NaC2Nbd *= -1.0
                    elif (op_c_type == op_d_type) and (op_d_type != op_b_type):
                        if op_b_type == 'real':
                            NaB3NcNd *= -1.0
                            NaB3NdNc *= -1.0
                            NaB2Ncd *= -1.0
                        else:
                            NaC3NbNd *= -1.0
                            NaC3NdNb *= -1.0
                            NaC2Nbd *= -1.0
                            NaD3NbNc *= -1.0
                            NaD3NcNb *= -1.0
                            NaD2Nbc *= -1.0
                    elif (op_d_type == op_b_type) and (op_b_type != op_c_type):
                        if op_c_type == 'real':
                            NaC3NbNd *= -1.0
                            NaC3NdNb *= -1.0
                            NaC2Nbd *= -1.0
                        else:
                            NaD3NbNc *= -1.0
                            NaD3NcNb *= -1.0
                            NaD2Nbc *= -1.0
                            NaB3NcNd *= -1.0
                            NaB3NdNc *= -1.0
                            NaB2Ncd *= -1.0
                    elif (op_b_type == op_c_type) and (op_c_type == op_d_type):
                        if op_d_type == 'imag':
                            NaB3NcNd *= -1.0
                            NaB3NdNc *= -1.0
                            NaB2Ncd *= -1.0
                            NaC3NbNd *= -1.0
                            NaC3NdNb *= -1.0
                            NaC2Nbd *= -1.0
                            NaD3NbNc *= -1.0
                            NaD3NcNb *= -1.0
                            NaD2Nbc *= -1.0

                val_T4 = -(NaE4NbNcNd - NaS4NbNcNd - NaR4NbNcNd)
                val_E3 = -(NaE3NbNcd)
                val_X3 = NaB3NcNd + NaB3NdNc + NaC3NbNd + NaC3NdNb + NaD3NbNc + NaD3NcNb
                val_X2 = NaB2Ncd + NaC2Nbd + NaD2Nbc
                val_A3 = -(NdA3NbNc + NdA3NcNb + NbA3NcNd + NbA3NdNc + NcA3NbNd + NcA3NdNb)
                val_A2 = NbA2Ncd + NcdA2Nb + NcA2Nbd + NbdA2Nc + NdA2Nbc + NbcA2Nd

                # flip sign for T4 and E3
                if op_a_type == 'real':
                    if (op_c_type == op_d_type) and (op_c_type != op_b_type):
                        val_T4 *= -1.0
                        val_E3 *= -1.0
                    if op_b_type == 'imag':
                        val_T4 *= -1.0
                        val_E3 *= -1.0
                elif op_a_type == 'imag':
                    if (op_b_type != op_c_type) or (op_c_type != op_d_type):
                        val_T4 *= -1.0
                        val_E3 *= -1.0

                # flip sign for A3 and A2
                if op_a_type == 'real':
                    if (op_c_type == op_d_type) and (op_c_type != op_b_type):
                        val_A3 *= -1.0
                        val_A2 *= -1.0
                    if op_b_type == 'imag':
                        val_A3 *= -1.0
                        val_A2 *= -1.0
                elif op_a_type == 'imag':
                    if (op_b_type == op_c_type) and (op_c_type == op_d_type):
                        val_A3 *= -1.0
                        val_A2 *= -1.0

                # Cubic response function
                gamma = val_T4 + val_E3 + val_X3 + val_A3 + val_X2 + val_A2

                self.ostream.print_blank()
                w_str = 'Cubic response function: << {};{},{},{} >>  ({},{},{})'.format(
                    self.a_component, self.b_component, self.c_component,
                    self.d_component, str(wb), str(wc), str(wd))
                self.ostream.print_header(w_str)
                self.ostream.print_header('=' * (len(w_str) + 2))
                self.ostream.print_blank()

                title = '{:<9s} {:>20s} {:>21s}'.format('', 'Real', 'Imaginary')
                width = len(title)
                self.ostream.print_header(title.ljust(width))
                self.ostream.print_header(('-' * len(title)).ljust(width))
                self._print_component('CRF', gamma, width)
                self.ostream.print_blank()

                result[('crf', wb, wc, wd)] = gamma

                result['crf_terms'] = {
                    ('crf_T4_term', wb, wc, wd): val_T4,
                    ('crf_E3_term', wb, wc, wd): val_E3,
                    ('crf_X3_term', wb, wc, wd): val_X3,
                    ('crf_X2_term', wb, wc, wd): val_X2,
                    ('crf_A3_term', wb, wc, wd): val_A3,
                    ('crf_A2_term', wb, wc, wd): val_A2,
                }

        profiler.check_memory_usage('End of CRF')

        return result

    def get_e3(self, wi, Nx, Nxy, fo, fo2, fo3, nocc, norb):
        """
        Contracts E[3] for CRF

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

            Nb = ComplexResponse.get_full_solution_vector(Nx[('B', wb)])
            Nc = ComplexResponse.get_full_solution_vector(Nx[('C', wc)])
            Nd = ComplexResponse.get_full_solution_vector(Nx[('D', wd)])

            Nbc = ComplexResponse.get_full_solution_vector(Nxy[(('BC', wb, wc),
                                                                wb + wc)])
            Nbd = ComplexResponse.get_full_solution_vector(Nxy[(('BD', wb, wd),
                                                                wb + wd)])
            Ncd = ComplexResponse.get_full_solution_vector(Nxy[(('CD', wc, wd),
                                                                wc + wd)])

            if self.rank != mpi_master():
                continue

            kb = (self.complex_lrvec2mat(Nb, nocc, norb)).T
            kc = (self.complex_lrvec2mat(Nc, nocc, norb)).T
            kd = (self.complex_lrvec2mat(Nd, nocc, norb)).T

            kbc = (self.complex_lrvec2mat(Nbc, nocc, norb)).T
            kbd = (self.complex_lrvec2mat(Nbd, nocc, norb)).T
            kcd = (self.complex_lrvec2mat(Ncd, nocc, norb)).T

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

            xi_b_cd = self._xi(kb, kcd, fb, fcd, F0_a)

            # E3NcNbd

            xi_c_bd = self._xi(kc, kbd, fc, fbd, F0_a)

            # E3NdNbc

            xi_d_bc = self._xi(kd, kbc, fd, fbc, F0_a)

            e3fock = (xi_b_cd + xi_c_bd + xi_d_bc).T + (0.5 * f123).T

            e3vec[('E3', (wb, wc, wd))] = self.anti_sym(
                -2 * LinearSolver.lrmat2vec(e3fock, nocc, norb))

        return e3vec

    def get_densities(self, freqtriples, Nx, mo, nocc, norb):
        """
        Computes the  densities needed for the Fock matrices.

        :param freqtriples:
            A list of the frequency triples
        :param Nx:
            A dictonary with all the first-order response vectors in
            distributed form
        :param mo:
            A matrix containing the MO coefficents
        :param nocc:
            Number of occupied orbitals
        :param nocc:
            Number of orbitals

        :return:
            A list of tranformed compounded densities
        """

        distributed_density_1 = None
        distributed_density_2 = None
        distributed_density_3 = None

        for (wb, wc, wd) in freqtriples:

            # convert response matrix to ao basis #

            Nb = ComplexResponse.get_full_solution_vector(Nx[('B', wb)])
            Nc = ComplexResponse.get_full_solution_vector(Nx[('C', wc)])
            Nd = ComplexResponse.get_full_solution_vector(Nx[('D', wd)])

            if self.rank == mpi_master():

                kb = self.complex_lrvec2mat(Nb, nocc, norb)
                kc = self.complex_lrvec2mat(Nc, nocc, norb)
                kd = self.complex_lrvec2mat(Nd, nocc, norb)

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

                dist_den_1_freq = np.hstack((
                    Db.real.reshape(-1, 1),
                    Db.imag.reshape(-1, 1),
                    Dc.real.reshape(-1, 1),
                    Dc.imag.reshape(-1, 1),
                    Dd.real.reshape(-1, 1),
                    Dd.imag.reshape(-1, 1),
                ))

                dist_den_2_freq = np.hstack((
                    Dbc.real.reshape(-1, 1),
                    Dbc.imag.reshape(-1, 1),
                    Dbd.real.reshape(-1, 1),
                    Dbd.imag.reshape(-1, 1),
                    Dcd.real.reshape(-1, 1),
                    Dcd.imag.reshape(-1, 1),
                ))

                dist_den_3_freq = np.hstack((
                    D123.real.reshape(-1, 1),
                    D123.imag.reshape(-1, 1),
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

    def get_densities_II(self, freqtriples, Nx, Nxy, mo, nocc, norb):
        """
        Computes the  densities needed for the Fock matrices.

        :param freqtriples:
            A list of the frequency triples
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

        for (wb, wc, wd) in freqtriples:

            Nb = ComplexResponse.get_full_solution_vector(Nx[('B', wb)])
            Nc = ComplexResponse.get_full_solution_vector(Nx[('C', wc)])
            Nd = ComplexResponse.get_full_solution_vector(Nx[('D', wd)])

            Nbc = ComplexResponse.get_full_solution_vector(Nxy[(('BC', wb, wc),
                                                                wb + wc)])
            Nbd = ComplexResponse.get_full_solution_vector(Nxy[(('BD', wb, wd),
                                                                wb + wd)])
            Ncd = ComplexResponse.get_full_solution_vector(Nxy[(('CD', wc, wd),
                                                                wc + wd)])

            if self.rank == mpi_master():

                kb = self.complex_lrvec2mat(Nb, nocc, norb)
                kc = self.complex_lrvec2mat(Nc, nocc, norb)
                kd = self.complex_lrvec2mat(Nd, nocc, norb)

                kbc = self.complex_lrvec2mat(Nbc, nocc, norb)
                kbd = self.complex_lrvec2mat(Nbd, nocc, norb)
                kcd = self.complex_lrvec2mat(Ncd, nocc, norb)

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

                dist_den_1_freq = np.hstack((
                    Db.real.reshape(-1, 1),
                    Db.imag.reshape(-1, 1),
                    Dc.real.reshape(-1, 1),
                    Dc.imag.reshape(-1, 1),
                    Dd.real.reshape(-1, 1),
                    Dd.imag.reshape(-1, 1),
                    Dbc.real.reshape(-1, 1),
                    Dbc.imag.reshape(-1, 1),
                    Dbd.real.reshape(-1, 1),
                    Dbd.imag.reshape(-1, 1),
                    Dcd.real.reshape(-1, 1),
                    Dcd.imag.reshape(-1, 1),
                ))

                dist_den_2_freq = np.hstack((
                    d_total.real.reshape(-1, 1),
                    d_total.imag.reshape(-1, 1),
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

    def get_fock_dict(self, wi, density_list1, density_list2, density_list3, F0,
                      mo, molecule, ao_basis, eri_dict, dft_dict):
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
        :param eri_dict:
            The dictionary containing ERI information
        :param dft_dict:
            The dictionary containing DFT information

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
            fpath = Path(self.checkpoint_file)
            fpath = fpath.with_name(fpath.stem)
            fock_file = str(fpath) + '_crf_fock_1.h5'
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
                                                 'real_and_imag', eri_dict,
                                                 dft_dict, density_list1,
                                                 density_list2, density_list3,
                                                 'crf')
            else:
                density_list_23 = DistributedArray(density_list2.data,
                                                   self.comm,
                                                   distribute=False)
                density_list_23.append(density_list3, axis=1)
                dist_focks = self._comp_nlr_fock(mo, molecule, ao_basis,
                                                 'real_and_imag', eri_dict,
                                                 None, None, None,
                                                 density_list_23, 'crf')

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
                         molecule, ao_basis, eri_dict, dft_dict):
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
        :param eri_dict:
            The dictionary containing ERI information
        :param dft_dict:
            The dictionary containing DFT information

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
            fpath = Path(self.checkpoint_file)
            fpath = fpath.with_name(fpath.stem)
            fock_file = str(fpath) + '_crf_fock_2.h5'
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
                                                 'real_and_imag', eri_dict,
                                                 dft_dict, density_list1,
                                                 density_list2, None, 'crf_ii')
            else:
                dist_focks = self._comp_nlr_fock(mo, molecule, ao_basis,
                                                 'real_and_imag', eri_dict,
                                                 None, None, density_list2,
                                                 None, 'crf_ii')

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

    def get_esr4(self, wi, Nx, fo, fo2, nocc, norb, D0):
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

            Na = ComplexResponse.get_full_solution_vector(Nx[('A',
                                                              (wb + wc + wd))])
            Nb = ComplexResponse.get_full_solution_vector(Nx[('B', wb)])
            Nc = ComplexResponse.get_full_solution_vector(Nx[('C', wc)])
            Nd = ComplexResponse.get_full_solution_vector(Nx[('D', wd)])

            if self.rank != mpi_master():
                continue

            vec_pack = vec_pack.T.copy().reshape(-1, norb, norb)

            (fb, fc, fd, fbc, fbd, fcd, fbcd) = vec_pack

            fb = np.conjugate(fb).T
            fc = np.conjugate(fc).T
            fd = np.conjugate(fd).T

            F0_a = fo['F0']

            # Response

            kb = (self.complex_lrvec2mat(Nb, nocc, norb)).T
            kc = (self.complex_lrvec2mat(Nc, nocc, norb)).T
            kd = (self.complex_lrvec2mat(Nd, nocc, norb)).T

            zi_bcd = self._zi(kb, kc, kd, fc, fd, fcd, F0_a)
            zi_cbd = self._zi(kc, kb, kd, fb, fd, fbd, F0_a)
            zi_dbc = self._zi(kd, kb, kc, fb, fc, fbc, F0_a)

            e4fock = (zi_bcd + zi_cbd + zi_dbc) + (fbcd)

            e4vec = 2. / 6 * self.anti_sym(
                LinearSolver.lrmat2vec(e4fock.T, nocc, norb))

            ka = self.complex_lrvec2mat(Na, nocc, norb)
            kb = self.complex_lrvec2mat(Nb, nocc, norb)
            kc = self.complex_lrvec2mat(Nc, nocc, norb)
            kd = self.complex_lrvec2mat(Nd, nocc, norb)

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

    def get_nxy(self, wi, Nx, fo, fo2, nocc, norb, d_a_mo, X, molecule,
                ao_basis, scf_tensors):
        """
        Computed NXY

        :param wi:
            A list of freqs
        :param Nx:
            A dict of the single index response vectors in distributed form
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

            Nb = ComplexResponse.get_full_solution_vector(Nx[('B', wb)])
            Nc = ComplexResponse.get_full_solution_vector(Nx[('C', wc)])
            Nd = ComplexResponse.get_full_solution_vector(Nx[('D', wd)])

            if self.rank != mpi_master():
                continue

            vec_pack = vec_pack.T.copy().reshape(-1, norb, norb)

            (fb, fc, fd, fbc, fbd, fcd) = vec_pack

            fb = np.conjugate(fb).T
            fc = np.conjugate(fc).T
            fd = np.conjugate(fd).T

            F0_a = fo['F0']

            # Response

            kb = (self.complex_lrvec2mat(Nb, nocc, norb)).T
            kc = (self.complex_lrvec2mat(Nc, nocc, norb)).T
            kd = (self.complex_lrvec2mat(Nd, nocc, norb)).T

            B = X[self._b_op_key][self.b_component]
            C = X[self._c_op_key][self.c_component]
            D = X[self._d_op_key][self.d_component]

            op_b_type = 'imag' if self.is_imag(self._b_op_key) else 'real'
            op_c_type = 'imag' if self.is_imag(self._c_op_key) else 'real'
            op_d_type = 'imag' if self.is_imag(self._d_op_key) else 'real'

            # BC

            xi = self._xi(kb, kc, fb, fc, F0_a)

            e3fock = xi.T + (0.5 * fbc).T
            E3NbNc = self.anti_sym(-LinearSolver.lrmat2vec(e3fock, nocc, norb))

            C2Nb = 0.5 * self._x2_contract(kb.T, C, d_a_mo, nocc, norb)
            B2Nc = 0.5 * self._x2_contract(kc.T, B, d_a_mo, nocc, norb)

            # flip sign for BC term
            if op_b_type == 'imag':
                B2Nc *= -1.0
            if op_c_type == 'imag':
                C2Nb *= -1.0

            BC[(('BC', wb, wc), wb + wc)] = E3NbNc - C2Nb - B2Nc

            # BD

            xi = self._xi(kb, kd, fb, fd, F0_a)

            e3fock = xi.T + (0.5 * fbd).T
            E3NbNd = self.anti_sym(-LinearSolver.lrmat2vec(e3fock, nocc, norb))

            D2Nb = 0.5 * self._x2_contract(kb.T, D, d_a_mo, nocc, norb)
            B2Nd = 0.5 * self._x2_contract(kd.T, B, d_a_mo, nocc, norb)

            # flip sign for BD term
            if op_b_type == 'imag':
                B2Nd *= -1.0
            if op_d_type == 'imag':
                D2Nb *= -1.0

            BD[(('BD', wb, wd), wb + wd)] = E3NbNd - D2Nb - B2Nd

            # CD

            xi = self._xi(kc, kd, fc, fd, F0_a)

            e3fock = xi.T + (0.5 * fcd).T
            E3NcNd = self.anti_sym(-LinearSolver.lrmat2vec(e3fock, nocc, norb))

            C2Nd = 0.5 * self._x2_contract(kd.T, C, d_a_mo, nocc, norb)
            D2Nc = 0.5 * self._x2_contract(kc.T, D, d_a_mo, nocc, norb)

            # flip sign for CD term
            if op_c_type == 'imag':
                C2Nd *= -1.0
            if op_d_type == 'imag':
                D2Nc *= -1.0

            CD[(('CD', wc, wd), wc + wd)] = E3NcNd - C2Nd - D2Nc

            XY.update(BC)
            XY.update(BD)
            XY.update(CD)

        Nxy_drv = ComplexResponse(self.comm, self.ostream)

        cpp_keywords = {
            'damping', 'norm_thresh', 'lindep_thresh', 'conv_thresh',
            'max_iter', 'eri_thresh', 'timing', 'memory_profiling',
            'batch_size', 'restart', 'xcfun', 'grid_level', 'potfile',
            'electric_field', 'program_end_time', '_debug', '_block_size_factor',
            'ri_coulomb'
        }

        for key in cpp_keywords:
            setattr(Nxy_drv, key, getattr(self, key))

        if self.checkpoint_file is not None:
            fpath = Path(self.checkpoint_file)
            fpath = fpath.with_name(fpath.stem)
            Nxy_drv.checkpoint_file = str(fpath) + '_crf_2.h5'

        Nxy_results = Nxy_drv.compute(molecule, ao_basis, scf_tensors, XY)

        self._is_converged = (self._is_converged and Nxy_drv.is_converged)

        Nxy = Nxy_results['solutions']
        Focks = Nxy_results['focks']

        return Nxy, Focks

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
