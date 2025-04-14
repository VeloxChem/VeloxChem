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
from .veloxchemlib import (mpi_master, bohr_in_angstrom, hartree_in_ev,
                           hartree_in_inverse_nm, fine_structure_constant,
                           speed_of_light_in_vacuum_in_SI)
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
from .lreigensolver import LinearResponseEigenSolver
from .firstorderprop import FirstOrderProperties


class ThreePATransitionDriver(NonlinearSolver):
    """
    Implements a general quadratic response driver.

    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.

    Instance variables
        - is_converged: The flag for convergence.
        - lindep_thresh: The threshold for removing linear dependence in the
          trial vectors.
        - conv_thresh: The convergence threshold for the solver.
        - max_iter: The maximum number of solver iterations.
        - a_components: Cartesian components of the A operator.
        - b_components: Cartesian components of the B operator.
        - c_components: Cartesian components of the C operator.
    """

    def __init__(self, comm=None, ostream=None):
        """
        Initializes the three-photon absorption driver.
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
        self.damping = 0.0

        # tpa transition settings
        self.nstates = 3

        # input keywords
        self._input_keywords['response'].update({
            'nstates': ('int', 'number of excited states'),
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

    def compute(self, molecule, ao_basis, scf_results):
        """
        Computes a quadratic response function.

        :param molecule:
            The molecule.
        :param basis:
            The AO basis.
        :param scf_results:
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
        scf_results_sanity_check(self, scf_results)

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
            self.print_header()

        start_time = time.time()

        # sanity check
        nalpha = molecule.number_of_alpha_electrons()
        nbeta = molecule.number_of_beta_electrons()
        assert_msg_critical(
            nalpha == nbeta,
            'TpaTransitionDriver: not implemented for unrestricted case')

        if self.rank == mpi_master():
            S = scf_results['S']
            da = scf_results['D_alpha']
            mo = scf_results['C_alpha']
            d_a_mo = np.linalg.multi_dot([mo.T, S, da, S, mo])
            norb = mo.shape[1]
        else:
            d_a_mo = None
            norb = None
        d_a_mo = self.comm.bcast(d_a_mo, root=mpi_master())
        norb = self.comm.bcast(norb, root=mpi_master())

        # Computing first-order gradient vectors
        if self.rank == mpi_master():
            mu_dipole_mats = compute_electric_dipole_integrals(
                molecule, ao_basis)
            # Note: nonliear response uses r instead of mu for dipole operator
            dipole_mats = (mu_dipole_mats[0] * (-1.0),
                           mu_dipole_mats[1] * (-1.0),
                           mu_dipole_mats[2] * (-1.0))
        else:
            dipole_mats = tuple()

        operator = 'dipole'
        linear_solver = LinearSolver(self.comm, self.ostream)
        a_grad = linear_solver.get_complex_prop_grad(operator, 'xyz', molecule,
                                                     ao_basis, scf_results)

        b_grad = linear_solver.get_complex_prop_grad(operator, 'xyz', molecule,
                                                     ao_basis, scf_results)

        if self.rank == mpi_master():
            inv_sqrt_2 = 1.0 / np.sqrt(2.0)

            a_grad = list(a_grad)
            for ind in range(len(a_grad)):
                a_grad[ind] *= inv_sqrt_2
                # Note: nonliear response uses r instead of mu for dipole operator
                if operator == 'dipole':
                    a_grad[ind] *= -1.0

            b_grad = list(b_grad)
            for ind in range(len(b_grad)):
                b_grad[ind] *= inv_sqrt_2
                # Note: nonliear response uses r instead of mu for dipole operator
                if operator == 'dipole':
                    b_grad[ind] *= -1.0

        rpa_drv = LinearResponseEigenSolver(self.comm, self.ostream)
        rpa_drv.nonlinear = True

        rpa_keywords = [
            'nstates', 'norm_thresh', 'lindep_thresh', 'conv_thresh',
            'max_iter', 'eri_thresh', 'timing', 'memory_profiling',
            'batch_size', 'restart', 'xcfun', 'grid_level', 'potfile',
            'electric_field', 'program_end_time', '_debug', '_block_size_factor'
        ]

        for key in rpa_keywords:
            setattr(rpa_drv, key, getattr(self, key))

        if self.checkpoint_file is not None:
            fpath = Path(self.checkpoint_file)
            fpath = fpath.with_name(fpath.stem)
            rpa_drv.checkpoint_file = str(fpath) + '_3patrans_1_rpa.h5'

        rpa_results = rpa_drv.compute(molecule, ao_basis, scf_results)

        excitation_details = rpa_results['excitation_details']
        oscillator_strengths = rpa_results['oscillator_strengths']
        elec_trans_dipoles = rpa_results['electric_transition_dipoles']

        Xf = {}
        inv_sqrt_2 = 1.0 / np.sqrt(2.0)
        for i in range(self.nstates):
            Xf[i] = DistributedArray(
                -inv_sqrt_2 * rpa_results['eigenvectors_distributed'][i].data,
                self.comm,
                distribute=False)

        freqs = [-1/3 * a for a in rpa_results['eigenvalues']]
        freqs_for_response_vectors = freqs + [0.0]

        # print("freqs_for_response_vectors",freqs_for_response_vectors)

        # Storing the dipole integral matrices used for the X[2] and
        # A[2] contractions in MO basis

        B = {}

        if self.rank == mpi_master():
            B.update({
                (op, w): v for op, v in zip('xyz', b_grad)
                for w in freqs_for_response_vectors
            })

            X = {
                'x': 2 * self.ao2mo(mo, dipole_mats[0]),
                'y': 2 * self.ao2mo(mo, dipole_mats[1]),
                'z': 2 * self.ao2mo(mo, dipole_mats[2])
            }
        else:
            X = None

        # Computing the first-order response vectors (3 per frequency)
        N_drv = ComplexResponse(self.comm, self.ostream)

        cpp_keywords = {
            'damping', 'norm_thresh', 'lindep_thresh', 'conv_thresh',
            'max_iter', 'eri_thresh', 'timing', 'memory_profiling',
            'batch_size', 'restart', 'xcfun', 'grid_level', 'potfile',
            'electric_field', 'program_end_time', '_debug', '_block_size_factor'
        }

        for key in cpp_keywords:
            setattr(N_drv, key, getattr(self, key))

        if self.checkpoint_file is not None:
            fpath = Path(self.checkpoint_file)
            fpath = fpath.with_name(fpath.stem)
            N_drv.checkpoint_file = str(fpath) + '_3patrans_1_cpp.h5'

        N_results = N_drv.compute(molecule, ao_basis, scf_results, B)

        self._is_converged = N_drv.is_converged

        Nx = N_results['solutions']
        Focks = N_results['focks']
        Focks.update(rpa_results['focks'])

        profiler.check_memory_usage('CPP')

        ret_dict = self.compute_cubic_components(Focks, freqs, X, d_a_mo, Nx,
                                                scf_results, molecule, ao_basis,
                                                profiler, Xf)

        valstr = '*** Time spent in 3PA response calculation: '
        valstr += '{:.2f} sec ***'.format(time.time() - start_time)
        self.ostream.print_header(valstr)
        self.ostream.print_blank()
        self.ostream.flush()

        profiler.end(self.ostream)

        if self.rank == mpi_master():
            ret_dict.update({
                'oscillator_strengths': oscillator_strengths,
                'elec_trans_dipoles': elec_trans_dipoles,
                'excitation_details': excitation_details,
            })

        return ret_dict

    def compute_cubic_components(self, Focks, freqs, X, d_a_mo, Nx, scf_results,
                                molecule, ao_basis, profiler, Xf):
        """
        Computes all the relevent terms to compute a general quadratic response function

        :param freqparis:
            A list of all the frequencies
        :param X:
            A dictonary of matricies containing all the property integrals
        :param d_a_mo:
            The SCF density in MO basis
        :param kX:
            A dictonary containing all the response matricies
        :param scf_results:
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
            mo = scf_results['C_alpha']
            F0 = np.linalg.multi_dot([mo.T, scf_results['F_alpha'], mo])
            norb = mo.shape[1]
        else:
            mo = None
            F0 = None
            norb = None

        F0 = self.comm.bcast(F0, root=mpi_master())
        norb = self.comm.bcast(norb, root=mpi_master())

        nocc = molecule.number_of_alpha_electrons()

        eri_dict = self._init_eri(molecule, ao_basis)

        dft_dict = self._init_dft(molecule, scf_results)

        # computing all compounded first-order densities
        first_order_dens, second_order_dens, third_order_dens = self.get_densities(freqs, Nx, mo, nocc, norb, Xf)

        profiler.check_memory_usage('Densities')

        # computing the compounded first-order Fock matrices
        fock_dict = self.get_fock_dict(freqs, first_order_dens,
                                       second_order_dens, third_order_dens, F0, mo, molecule,
                                       ao_basis, eri_dict, dft_dict, profiler)

        profiler.check_memory_usage('Focks')

        e4_dict, s4_dict = self.get_es4(freqs, Nx, fock_dict, Focks, nocc, norb, d_a_mo,Xf)

        Nxy, f_xy = self.get_nxy(freqs, Nx, fock_dict, Focks, nocc, norb, d_a_mo, X, molecule, ao_basis, scf_results,Xf)

        # Needs tgo be fixed
        density_list1_two, density_list2_two = self.get_densities_II(freqs, Nx, Nxy, mo, nocc, norb,Xf)
        
        # Needs tgo be fixed
        fock_dict_two = self.get_fock_dict_II(freqs, density_list1_two,density_list2_two, F0, mo,molecule, ao_basis, eri_dict,dft_dict)

        e3_dict = self.get_e3(freqs, Nx, Nxy, Focks, fock_dict_two, f_xy, nocc, norb, Xf)

        profiler.check_memory_usage('E[3]')

        ret_dict = {}
        T_tensors = {}

        # Compute dipole vector
        scf_prop = FirstOrderProperties(self.comm, self.ostream)
        scf_prop.compute_scf_prop(molecule, ao_basis, scf_results)

        for w_ind, w in enumerate(freqs):
            for a_comp in 'xyz':
                for b_comp in 'xyz':
                    for c_comp in 'xyz':

                        Na_raw = ComplexResponse.get_full_solution_vector(Nx[(a_comp, w)])

                        Nb = ComplexResponse.get_full_solution_vector(Nx[(b_comp, w)])
                        Nc = ComplexResponse.get_full_solution_vector(Nx[(c_comp, w)])

                        Ncf = ComplexResponse.get_full_solution_vector(Nxy[(f"f{c_comp}", - 2 * w)])
                        Nbf = ComplexResponse.get_full_solution_vector(Nxy[(f"f{b_comp}", - 2 * w)])
                        try:
                            Nbc = ComplexResponse.get_full_solution_vector(Nxy[(f"{b_comp}{c_comp}", 2 * w)])
                        except KeyError:
                            Nbc = ComplexResponse.get_full_solution_vector(Nxy[(f"{c_comp}{b_comp}", 2 * w)])

                        Nf_raw = LinearResponseEigenSolver.get_full_solution_vector(Xf[w_ind])

                        if self.rank == mpi_master():

                            Na = self.flip_yz(Na_raw)
                            Nf = -Nf_raw

                            op_A = X[a_comp]
                            op_B = X[b_comp]
                            op_C = X[c_comp]

                            kb = LinearSolver.lrvec2mat(Nb, nocc, norb)
                            kc = LinearSolver.lrvec2mat(Nc, nocc, norb)
                            kf = LinearSolver.lrvec2mat(Nf, nocc, norb)
                            kbf = LinearSolver.lrvec2mat(Nbf, nocc, norb)
                            kcf = LinearSolver.lrvec2mat(Ncf, nocc, norb)
                            kbc = LinearSolver.lrvec2mat(Nbc, nocc, norb)

                            # A3 contraction
                            A3NbNcNf = np.dot(self._a3_contract(kb, kc, op_A, d_a_mo, nocc, norb), Nf)
                            A3NbNcNf += np.dot(self._a3_contract(kc, kb, op_A, d_a_mo, nocc, norb), Nf)

                            A3NbNcNf += np.dot(self._a3_contract(kf, kc, op_A, d_a_mo, nocc, norb), Nb)
                            A3NbNcNf += np.dot(self._a3_contract(kc, kf, op_A, d_a_mo, nocc, norb), Nb)

                            A3NbNcNf += np.dot(self._a3_contract(kb, kf, op_A, d_a_mo, nocc, norb), Nc)
                            A3NbNcNf += np.dot(self._a3_contract(kf, kb, op_A, d_a_mo, nocc, norb), Nc)

                            
                            # X3 contraction
                            NaB3NcNf = np.dot(Na.T, self._x3_contract(kc, kf, op_B, d_a_mo, nocc, norb))
                            NaB3NfNc = np.dot(Na.T, self._x3_contract(kf, kc, op_B, d_a_mo, nocc, norb))
                            NaC3NbNf = np.dot(Na.T, self._x3_contract(kb, kf, op_C, d_a_mo, nocc, norb))
                            NaC3NfNb = np.dot(Na.T, self._x3_contract(kf, kb, op_C, d_a_mo, nocc, norb))

                            # A2 contraction
                            NbA2Ncf = np.dot(self._a2_contract(kb, op_A, d_a_mo, nocc, norb), Ncf)
                            NbA2Ncf += np.dot(self._a2_contract(kcf, op_A, d_a_mo, nocc, norb), Nb)

                            NcA2Nbf = np.dot(self._a2_contract(kc, op_A, d_a_mo, nocc, norb), Nbf)
                            NcA2Nbf += np.dot(self._a2_contract(kbf, op_A, d_a_mo, nocc, norb), Nc)

                            NfA2Nbc = np.dot(self._a2_contract(-kf, op_A, d_a_mo, nocc, norb), Nbc)
                            NfA2Nbc += np.dot(self._a2_contract(kbc, op_A, d_a_mo, nocc, norb), -Nf)

                            # X2 contraction
                            NaB2Ncf = np.dot(Na.T,self._x2_contract(kcf, op_B, d_a_mo, nocc, norb))
                            NaC2Nbf = np.dot(Na.T,self._x2_contract(kbf, op_C, d_a_mo, nocc, norb))

                            # E3 contractions  
                            NaE3NbNcf = np.dot(Na, e3_dict[(f"{b_comp}{c_comp}", w)])

                            # T4 contraction
                            NaE4NbNcNf = np.dot(Na, e4_dict[(f"{b_comp}{c_comp}", w)])
                            NaS4NbNcNd = np.dot(Na, s4_dict[(f"{b_comp}{c_comp}", w)])

                            # Check signs of the contractions
                            T_abc = A3NbNcNf 
                            T_abc +=  -(NaB3NcNf + NaB3NfNc + NaC3NbNf + NaC3NfNb)
                            T_abc +=  NbA2Ncf + NcA2Nbf + NfA2Nbc
                            T_abc +=  NaB2Ncf + NaC2Nbf 
                            T_abc += - NaE3NbNcf 
                            T_abc +=  -(NaE4NbNcNf - NaS4NbNcNd)

                            T_tensors[(f"{a_comp}{b_comp}{c_comp}", round(float(w), 4))] = T_abc


        diagonalized_tensors = {}

        #if self.rank == mpi_master():
            #for key, tensor in M_tensors.items():
            #    eigvals, eigvecs = np.linalg.eigh(tensor)
            #    diagonalized_tensors[key] = np.diag(eigvals)

        self.ostream.print_blank()
        w_str = 'Summary of Three-photon Absorption'
        self.ostream.print_header(w_str)
        self.ostream.print_header('=' * (len(w_str) + 2))
        self.ostream.print_blank()

        # conversion factor for 3PA cross-sections in GM
        # a0 in cm, c in cm/s, gamma (0.1 eV) in au
        alpha = fine_structure_constant()
        a0 = bohr_in_angstrom() * 1.0e-8
        c = speed_of_light_in_vacuum_in_SI() * 100.0
        gamma = 0.1 / hartree_in_ev()
        # TODO: double check au2gm for 3PA
        au2gm = (4.0 * np.pi**3 * alpha * a0**8) / (3.0 * c **2  * gamma) * 1.0e+50

        if self.rank == mpi_master():
            tpa_strengths = {'linear': {}, 'circular': {}}
            tpa_cross_sections = {'linear': {}, 'circular': {}}

            for w in freqs:
                Df = 0.0
                Dg = 0.0
                for i in 'xyz':
                    for j in 'xyz':
                        for k in 'xyz':
                            Tiij = T_tensors[(f"{i}{i}{j}", round(float(w), 4))]
                            Tkkj = T_tensors[(f"{k}{k}{j}", round(float(w), 4))]
                            Tijk = T_tensors[(f"{i}{j}{k}", round(float(w), 4))]

                            Df += (Tiij * Tkkj).real
                            Dg += (Tijk * Tijk).real
                Df /= 35.0
                Dg /= 35.0

                D_linear = 2.0 * Dg + 3.0 * Df
                D_circular = 5.0 * Dg - 3.0 * Df

                tpa_strengths['linear'][w] = D_linear
                tpa_strengths['circular'][w] = D_circular

                # TODO: double check 3PA cross-sections

                tpa_cross_sections['linear'][w] = au2gm * w**3 * D_linear
                tpa_cross_sections['circular'][w] = au2gm * w**3 * D_circular

            profiler.check_memory_usage('End of 3PA')

        
            #a_target, b_target, c_target = 'z', 'z', 'y'

            #for key, value in debugg_dict.items():
            #    if isinstance(key, tuple) and len(key) == 3:
            #        abc, w, T = key
            #        if abc == f"{a_target}{b_target}{c_target}":
            #            print(f"{key}: {value}")

            ret_dict = {
                'photon_energies': [-w for w in freqs],
                'transition_moments': T_tensors,
                #'cross_sections': tpa_cross_sections,
                '3pa_strengths': tpa_strengths,
                'ground_state_dipole_moments':scf_prop.get_property('dipole moment')
            }


            self._print_results(ret_dict)

            return ret_dict
        else:
            return None, None, None, None


    def get_es4(self, wi, N, fo, fo2, nocc, norb, D0,Xf):
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

        for w_ind, w in enumerate(wi):

            vec_pack = np.array([
                fo['Fbfx'][w].data,
                fo['Fbfy'][w].data,
                fo['Fbfz'][w].data,
                fo['Fbc_xx'][w].data,
                fo['Fbc_yy'][w].data,
                fo['Fbc_zz'][w].data,
                fo['Fbc_xy'][w].data,
                fo['Fbc_xz'][w].data,
                fo['Fbc_yz'][w].data,
                fo['Fbcf_yz'][w].data,
                fo['Fbcf_xz'][w].data,
                fo['Fbcf_xy'][w].data,
                fo['Fbcf_xx'][w].data,
                fo['Fbcf_yy'][w].data,
                fo['Fbcf_zz'][w].data,
                fo2[('x', w)].data,
                fo2[('y', w)].data,
                fo2[('z', w)].data,
                fo2[w_ind].data,
            ]).T.copy()

            vec_pack = self._collect_vectors_in_columns(vec_pack)

            Nx = ComplexResponse.get_full_solution_vector(N[('x', w)])
            Ny = ComplexResponse.get_full_solution_vector(N[('y', w)])
            Nz = ComplexResponse.get_full_solution_vector(N[('z', w)])
            Nf = LinearResponseEigenSolver.get_full_solution_vector(Xf[w_ind])

            if self.rank != mpi_master():
                continue

            vec_pack = vec_pack.T.copy().reshape(-1, norb, norb)

            (Fbfx, 
             Fbfy, 
             Fbfz, 
             Fbc_xx, 
             Fbc_yy, 
             Fbc_zz, 
             Fbc_xy, 
             Fbc_xz, 
             Fbc_yz, 
             Fbcf_yz,
             Fbcf_xz,
             Fbcf_xy,
             Fbcf_xx,
             Fbcf_yy,
             Fbcf_zz,
             fx, 
             fy, 
             fz, 
             ff) = vec_pack

            # First-order Fock matrices
            fx = np.conjugate(fx).T
            fy = np.conjugate(fy).T
            fz = np.conjugate(fz).T
            ff = np.conjugate(ff).T * -1 / np.sqrt(2)

            F0_a = fo['F0']

            # Response

            kx = (self.complex_lrvec2mat(Nx, nocc, norb)).T
            ky = (self.complex_lrvec2mat(Ny, nocc, norb)).T
            kz = (self.complex_lrvec2mat(Nz, nocc, norb)).T
            kf = (self.complex_lrvec2mat(Nf, nocc, norb)).T

            # xx 
            zi_bcd_xx = self._zi(kx, kx, kf, fx, ff, Fbfx, F0_a)
            zi_cbd_xx = self._zi(kx, kx, kf, fx, ff, Fbfx, F0_a)
            zi_dbc_xx = self._zi(kf, kx, kx, fx, fx, Fbc_xx, F0_a)

            e4fock_xx = (zi_bcd_xx + zi_cbd_xx + zi_dbc_xx) + (Fbcf_xx)
            e4_vec[('xx',w)] = 2. / 6 * self.anti_sym(LinearSolver.lrmat2vec(e4fock_xx.T, nocc, norb))

            # yy 
            zi_bcd_yy = self._zi(ky, ky, kf, fy, ff, Fbfy, F0_a)
            zi_cbd_yy = self._zi(ky, ky, kf, fy, ff, Fbfy, F0_a)
            zi_dbc_yy = self._zi(kf, ky, ky, fy, fy, Fbc_yy, F0_a)

            e4fock_yy = (zi_bcd_yy + zi_cbd_yy + zi_dbc_yy) + (Fbcf_yy)
            e4_vec[('yy',w)] = 2. / 6 * self.anti_sym(LinearSolver.lrmat2vec(e4fock_yy.T, nocc, norb))

            # zz 
            zi_bcd_zz = self._zi(kz, kz, kf, fz, ff, Fbfz, F0_a)
            zi_cbd_zz = self._zi(kz, kz, kf, fz, ff, Fbfz, F0_a)
            zi_dbc_zz = self._zi(kf, kz, kz, fz, fz, Fbc_zz, F0_a)

            e4fock_zz = (zi_bcd_zz + zi_cbd_zz + zi_dbc_zz) + (Fbcf_zz)
            e4_vec[('zz',w)] = 2. / 6 * self.anti_sym(LinearSolver.lrmat2vec(e4fock_zz.T, nocc, norb))


            ## xy 
            zi_bcd_xy = self._zi(kx, ky, kf, fy, ff, Fbfy, F0_a)
            zi_cbd_xy = self._zi(ky, kx, kf, fx, ff, Fbfx, F0_a)
            zi_dbc_xy = self._zi(kf, kx, ky, fx, fy, Fbc_xy, F0_a)

            e4fock_xy = (zi_bcd_xy + zi_cbd_xy + zi_dbc_xy) + (Fbcf_xy)
            e4_vec[('xy',w)] = 2. / 6 * self.anti_sym(LinearSolver.lrmat2vec(e4fock_xy.T, nocc, norb))
            e4_vec[('yx',w)] = e4_vec[('xy',w)]

            # xz
            zi_bcd_xz = self._zi(kx, kz, kf, fz, ff, Fbfz, F0_a)
            zi_cbd_xz = self._zi(kz, kx, kf, fx, ff, Fbfx, F0_a)
            zi_dbc_xz = self._zi(kf, kx, kz, fx, fz, Fbc_xz, F0_a)

            e4fock_xz = (zi_bcd_xz + zi_cbd_xz + zi_dbc_xz) + (Fbcf_xz)
            e4_vec[('xz',w)] = 2. / 6 * self.anti_sym(LinearSolver.lrmat2vec(e4fock_xz.T, nocc, norb))
            e4_vec[('zx',w)] = e4_vec[('xz',w)]

            # yz
            zi_bcd_yz = self._zi(ky, kz, kf, fz, ff, Fbfz, F0_a)
            zi_cbd_yz = self._zi(kz, ky, kf, fy, ff, Fbfy, F0_a)
            zi_dbc_yz = self._zi(kf, ky, kz, fy, fz, Fbc_yz, F0_a)

            e4fock_yz = (zi_bcd_yz + zi_cbd_yz + zi_dbc_yz) + (Fbcf_yz)
            e4_vec[('yz',w)] = 2. / 6 * self.anti_sym(LinearSolver.lrmat2vec(e4fock_yz.T, nocc, norb))
            e4_vec[('zy',w)] = e4_vec[('yz',w)]

            # xx
            kx = kx.T 
            ky = ky.T 
            kz = kz.T 
            kf = kf.T 

            # xx
            s4_term_xx = w * self._s4(kx, kx, kf, D0, nocc, norb)
            s4_term_xx += w * self._s4(kx, kx, kf, D0, nocc, norb)
            s4_term_xx += - 3 * w * self._s4(kf, kx, kx, D0, nocc, norb)
            s4_vec[('xx',w)] = s4_term_xx

            # yy
            s4_term_yy = w * self._s4(ky, ky, kf, D0, nocc, norb)
            s4_term_yy += w * self._s4(ky, ky, kf, D0, nocc, norb)
            s4_term_yy += - 3 * w * self._s4(kf, ky, ky, D0, nocc, norb)
            s4_vec[('yy',w)] = s4_term_yy

            # zz
            s4_term_zz = w * self._s4(kz, kz, kf, D0, nocc, norb)
            s4_term_zz += w * self._s4(kz, kz, kf, D0, nocc, norb)
            s4_term_zz += - 3 * w * self._s4(kf, kz, kz, D0, nocc, norb)
            s4_vec[('zz',w)] = s4_term_zz

            # xy
            s4_term_xy = w * self._s4(kx, ky, kf, D0, nocc, norb)
            s4_term_xy += w * self._s4(ky, kx, kf, D0, nocc, norb)
            s4_term_xy += - 3 * w * self._s4(kf, kx, ky, D0, nocc, norb)
            s4_vec[('xy',w)] = s4_term_xy
            s4_vec[('yx',w)] = s4_term_xy

            # xz
            s4_term_xz = w * self._s4(kx, kz, kf, D0, nocc, norb)
            s4_term_xz += w * self._s4(kz, kx, kf, D0, nocc, norb)
            s4_term_xz += - 3 * w * self._s4(kf, kx, kz, D0, nocc, norb)
            s4_vec[('xz',w)] = s4_term_xz
            s4_vec[('zx',w)] = s4_term_xz

            # yz
            s4_term_yz = w * self._s4(ky, kz, kf, D0, nocc, norb)
            s4_term_yz += w * self._s4(kz, ky, kf, D0, nocc, norb)
            s4_term_yz += - 3 * w * self._s4(kf, ky, kz, D0, nocc, norb)
            s4_vec[('yz',w)] = s4_term_yz
            s4_vec[('zy',w)] = s4_term_yz


        return e4_vec, s4_vec


    def get_nxy(self, freqs, Nx, fo, fo2, nocc, norb, d_a_mo, X, molecule,ao_basis, scf_tensors,Xf):
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
        # We need N_xf, N_yf, N_zf

        BC = {}
        XY = {}

        for w_ind, w in enumerate(freqs):

            vec_pack = np.array([
                fo2[('x', w)].data,
                fo2[('y', w)].data,
                fo2[('z', w)].data,
                fo2[w_ind].data,
                fo['Fbfx'][w].data,
                fo['Fbfy' ][w].data,
                fo['Fbfz'][w].data,
                fo['Fbc_xx'][w].data,
                fo['Fbc_yy'][w].data,
                fo['Fbc_zz'][w].data,
                fo['Fbc_xy'][w].data,
                fo['Fbc_xz'][w].data,
                fo['Fbc_yz'][w].data,
            ]).T.copy()

            vec_pack = self._collect_vectors_in_columns(vec_pack)

            nx = ComplexResponse.get_full_solution_vector(Nx[('x', w)])
            ny = ComplexResponse.get_full_solution_vector(Nx[('y', w)])
            nz = ComplexResponse.get_full_solution_vector(Nx[('z', w)])

            Nf = LinearResponseEigenSolver.get_full_solution_vector(Xf[w_ind])


            if self.rank != mpi_master():
                continue

            vec_pack = vec_pack.T.copy().reshape(-1, norb, norb)


            # Fix Fock matrices 

            (fx, 
             fy, 
             fz,
             ff,
             Fbfx, 
             Fbfy, 
             Fbfz, 
             Fbc_xx, 
             Fbc_yy, 
             Fbc_zz, 
             Fbc_xy, 
             Fbc_xz, 
             Fbc_yz) = vec_pack

            fx = np.conjugate(fx).T
            fy = np.conjugate(fy).T
            fz = np.conjugate(fz).T
            ff = np.conjugate(ff).T * -1 / np.sqrt(2)

            F0_a = fo['F0']

            # Response

            kx = (self.complex_lrvec2mat(nx, nocc, norb)).T
            ky = (self.complex_lrvec2mat(ny, nocc, norb)).T
            kz = (self.complex_lrvec2mat(nz, nocc, norb)).T
            kf = LinearSolver.lrvec2mat(Nf, nocc, norb).T

            # N_xf

            xi_x = self._xi(kx, kf, fx, ff, F0_a)
            xi_y = self._xi(ky, kf, fy, ff, F0_a)
            xi_z = self._xi(kz, kf, fz, ff, F0_a)

            e3fock_x = xi_x.T + (0.5 * Fbfx).T
            e3fock_y = xi_y.T + (0.5 * Fbfy).T
            e3fock_z = xi_z.T + (0.5 * Fbfz).T

            E3NbNf_x = self.anti_sym(-LinearSolver.lrmat2vec(e3fock_x, nocc, norb))
            E3NbNf_y = self.anti_sym(-LinearSolver.lrmat2vec(e3fock_y, nocc, norb))
            E3NbNf_z = self.anti_sym(-LinearSolver.lrmat2vec(e3fock_z, nocc, norb))

            X2Nf = 0.5 * self._x2_contract(kf.T, X['x'], d_a_mo, nocc, norb)
            Y2Nf = 0.5 * self._x2_contract(kf.T, X['y'], d_a_mo, nocc, norb)
            Z2Nf = 0.5 * self._x2_contract(kf.T, X['z'], d_a_mo, nocc, norb)

            BC[('fx', - 2 * w)] = E3NbNf_x - X2Nf
            BC[('fy', - 2 * w)] = E3NbNf_y - Y2Nf
            BC[('fz', - 2 * w)] = E3NbNf_z - Z2Nf

            # E3 contractions 
            xi_xx = self._xi(kx, kx, fx, fx, F0_a)
            xi_yy = self._xi(ky, ky, fy, fy, F0_a)
            xi_zz = self._xi(kz, kz, fz, fz, F0_a)
            xi_xy = self._xi(kx, ky, fx, fy, F0_a)
            xi_xz = self._xi(kx, kz, fx, fz, F0_a)
            xi_yz = self._xi(ky, kz, fy, fz, F0_a)

            e3fock_xx = xi_xx.T + (0.5 * Fbc_xx).T
            e3fock_yy = xi_yy.T + (0.5 * Fbc_yy).T
            e3fock_zz = xi_zz.T + (0.5 * Fbc_zz).T
            e3fock_xy = xi_xy.T + (0.5 * Fbc_xy).T
            e3fock_xz = xi_xz.T + (0.5 * Fbc_xz).T
            e3fock_yz = xi_yz.T + (0.5 * Fbc_yz).T

            E3NbNf_xx = self.anti_sym(-LinearSolver.lrmat2vec(e3fock_xx, nocc, norb))
            E3NbNf_yy = self.anti_sym(-LinearSolver.lrmat2vec(e3fock_yy, nocc, norb))
            E3NbNf_zz = self.anti_sym(-LinearSolver.lrmat2vec(e3fock_zz, nocc, norb))
            E3NbNf_xy = self.anti_sym(-LinearSolver.lrmat2vec(e3fock_xy, nocc, norb))
            E3NbNf_xz = self.anti_sym(-LinearSolver.lrmat2vec(e3fock_xz, nocc, norb))
            E3NbNf_yz = self.anti_sym(-LinearSolver.lrmat2vec(e3fock_yz, nocc, norb))

            # X2 contractions 
            B2Nc_xx =  0.5 * self._x2_contract(kx.T, X['x'], d_a_mo, nocc, norb) +  0.5 * self._x2_contract(kx.T, X['x'], d_a_mo, nocc, norb)
            B2Nc_yy =  0.5 * self._x2_contract(ky.T, X['y'], d_a_mo, nocc, norb) + 0.5 * self._x2_contract(ky.T, X['y'], d_a_mo, nocc, norb)
            B2Nc_zz =  0.5 * self._x2_contract(kz.T, X['z'], d_a_mo, nocc, norb) + 0.5 * self._x2_contract(kz.T, X['z'], d_a_mo, nocc, norb)
            B2Nc_xy = 0.5 * self._x2_contract(kx.T, X['y'], d_a_mo, nocc, norb) + 0.5 * self._x2_contract(ky.T, X['x'], d_a_mo, nocc, norb) 
            B2Nc_xz = 0.5 * self._x2_contract(kx.T, X['z'], d_a_mo, nocc, norb) + 0.5 * self._x2_contract(kz.T, X['x'], d_a_mo, nocc, norb)
            B2Nc_yz = 0.5 * self._x2_contract(ky.T, X['z'], d_a_mo, nocc, norb) + 0.5 * self._x2_contract(kz.T, X['y'], d_a_mo, nocc, norb)

            BC[('xx', 2 * w)] = E3NbNf_xx - B2Nc_xx 
            BC[('yy', 2 * w)] = E3NbNf_yy - B2Nc_yy 
            BC[('zz', 2 * w)] = E3NbNf_zz - B2Nc_zz 
            BC[('xy', 2 * w)] = E3NbNf_xy - B2Nc_xy
            BC[('xz', 2 * w)] = E3NbNf_xz - B2Nc_xz
            BC[('yz', 2 * w)] = E3NbNf_yz - B2Nc_yz

            XY.update(BC)


        Nxy_drv = ComplexResponse(self.comm, self.ostream)

        cpp_keywords = {
            'damping', 'norm_thresh', 'lindep_thresh', 'conv_thresh',
            'max_iter', 'eri_thresh', 'timing', 'memory_profiling',
            'batch_size', 'restart', 'xcfun', 'grid_level', 'potfile',
            'electric_field', 'program_end_time', '_debug', '_block_size_factor'
        }

        for key in cpp_keywords:
            setattr(Nxy_drv, key, getattr(self, key))

        if self.checkpoint_file is not None:
            fpath = Path(self.checkpoint_file)
            fpath = fpath.with_name(fpath.stem)
            Nxy_drv.checkpoint_file = str(fpath) + '_3patrans_2.h5'

        Nxy_results = Nxy_drv.compute(molecule, ao_basis, scf_tensors, XY)

        self._is_converged = (self._is_converged and Nxy_drv.is_converged)

        Nxy = Nxy_results['solutions']
        Focks = Nxy_results['focks']

        return Nxy, Focks


    def get_densities_II(self, freq, N, N2, mo, nocc, norb,Xf):
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

        for w_ind, w in enumerate(freq):

            Nx = ComplexResponse.get_full_solution_vector(N[('x', w)])
            Ny = ComplexResponse.get_full_solution_vector(N[('y', w)])
            Nz = ComplexResponse.get_full_solution_vector(N[('z', w)])
            Nf = LinearResponseEigenSolver.get_full_solution_vector(Xf[w_ind])

            Nfx = ComplexResponse.get_full_solution_vector(N2[('fx', - 2 * w)])
            Nfy = ComplexResponse.get_full_solution_vector(N2[('fy', - 2 * w)])
            Nfz = ComplexResponse.get_full_solution_vector(N2[('fz', - 2 * w)])

            Nxx = ComplexResponse.get_full_solution_vector(N2[('xx', 2 * w)])
            Nyy = ComplexResponse.get_full_solution_vector(N2[('yy', 2 * w)])
            Nzz = ComplexResponse.get_full_solution_vector(N2[('zz', 2 * w)])
            Nxy = ComplexResponse.get_full_solution_vector(N2[('xy', 2 * w)])
            Nxz = ComplexResponse.get_full_solution_vector(N2[('xz', 2 * w)])
            Nyz = ComplexResponse.get_full_solution_vector(N2[('yz', 2 * w)])
            

            if self.rank == mpi_master():

                kx = self.complex_lrvec2mat(Nx, nocc, norb)
                ky = self.complex_lrvec2mat(Ny, nocc, norb)
                kz = self.complex_lrvec2mat(Nz, nocc, norb)
                kf = self.complex_lrvec2mat(Nf, nocc, norb)

                kfx = self.complex_lrvec2mat(Nfx, nocc, norb)
                kfy = self.complex_lrvec2mat(Nfy, nocc, norb)
                kfz = self.complex_lrvec2mat(Nfz, nocc, norb)

                kxx = self.complex_lrvec2mat(Nxx, nocc, norb)
                kyy = self.complex_lrvec2mat(Nyy, nocc, norb)
                kzz = self.complex_lrvec2mat(Nzz, nocc, norb)
                kxy = self.complex_lrvec2mat(Nxy, nocc, norb)
                kxz = self.complex_lrvec2mat(Nxz, nocc, norb)
                kyz = self.complex_lrvec2mat(Nyz, nocc, norb)
            
                # create the first order single indexed densiteies #

                Dx = self.commut_mo_density(kx, nocc)
                Dy = self.commut_mo_density(ky, nocc)
                Dz = self.commut_mo_density(kz, nocc)
                Df = self.commut_mo_density(kf, nocc)

                # create the second-order two indexed densities #

                Dfx = self.commut_mo_density(kfx, nocc)
                Dfy = self.commut_mo_density(kfy, nocc)
                Dfz = self.commut_mo_density(kfz, nocc)

                Dxx = self.commut_mo_density(kxx, nocc)
                Dyy = self.commut_mo_density(kyy, nocc)
                Dzz = self.commut_mo_density(kzz, nocc)
                Dxy = self.commut_mo_density(kxy, nocc)
                Dxz = self.commut_mo_density(kxz, nocc)
                Dyz = self.commut_mo_density(kyz, nocc)

                # create the second-order three indexed densities #

                Dx_fx = self.commut(kx, Dfx)
                Dfx_x = self.commut(kfx, Dx)

                Dy_fy = self.commut(ky, Dfy)
                Dfy_y = self.commut(kfy, Dy)

                Dz_fz = self.commut(kz, Dfz)
                Dfz_z = self.commut(kfz, Dz)

                # New terms

                Dx_fy = self.commut(kx, Dfy)
                Dfy_x = self.commut(kfy, Dx)

                Dx_fz = self.commut(kx, Dfz)
                Dfz_x = self.commut(kfz, Dx)

                Dy_fx = self.commut(ky, Dfx)
                Dfx_y = self.commut(kfx, Dy)

                Dy_fz = self.commut(ky, Dfz)
                Dfz_y = self.commut(kfz, Dy)

                Dz_fx = self.commut(kz, Dfx)
                Dfx_z = self.commut(kfx, Dz)

                Dz_fy = self.commut(kz, Dfy)
                Dfy_z = self.commut(kfy, Dz)

                ###

                Df_xx = self.commut(kf, Dxx)
                Dxx_f = self.commut(kxx, Df)

                Df_yy = self.commut(kf, Dyy)
                Dyy_f = self.commut(kyy, Df)

                Df_zz = self.commut(kf, Dzz)
                Dzz_f = self.commut(kzz, Df)

                Df_xy = self.commut(kf, Dxy)
                Dxy_f = self.commut(kxy, Df)

                Df_xz = self.commut(kf, Dxz)
                Dxz_f = self.commut(kxz, Df)

                Df_yz = self.commut(kf, Dyz)
                Dyz_f = self.commut(kyz, Df)
                
                # density transformation from MO to AO basis

                # xx
                # 0.5 * fxfx.T +  0.5 * fxfx.T  +   0.5 * ffxx.T
                Dfxx = np.linalg.multi_dot([mo, (2.0 * (Dx_fx + Dfx_x) + Df_xx + Dxx_f), mo.T])

                # yy 
                # 0.5 * fyfy.T +  0.5 * fyfy.T  +   0.5 * ffyy.T
                Dfyy = np.linalg.multi_dot([mo, (2.0 * (Dy_fy + Dfy_y) + Df_yy + Dyy_f), mo.T])

                # zz
                # 0.5 * fzfz.T +  0.5 * fzfz.T  +   0.5 * ffzz.T
                Dfzz = np.linalg.multi_dot([mo, (2.0 * (Dz_fz + Dfz_z) + Df_zz + Dzz_f), mo.T])

                #xy
                # 0.5 * fxfy.T +  0.5 * fyfx.T  +   0.5 * ffxy.T
                Dfxy = np.linalg.multi_dot([mo, (Dx_fy + Dfy_x) + (Dy_fx + Dfx_y) + (Df_xy + Dxy_f), mo.T])

                #xz
                # 0.5 * fxfz.T +  0.5 * fzfx.T  +   0.5 * ffxz.T
                Dfxz = np.linalg.multi_dot([mo, (Dx_fz + Dfz_x) + (Dz_fx + Dfx_z) + (Df_xz + Dxz_f), mo.T])

                #yz
                # 0.5 * fyfz.T +  0.5 * fzyf.T  +   0.5 * ffyz.T
                Dfyz = np.linalg.multi_dot([mo, (Dy_fz + Dfz_y) + (Dz_fy + Dfy_z) + (Df_yz + Dyz_f), mo.T])
                

                Dx = np.linalg.multi_dot([mo, Dx, mo.T])
                Dy = np.linalg.multi_dot([mo, Dy, mo.T])
                Dz = np.linalg.multi_dot([mo, Dz, mo.T])
                Df = np.linalg.multi_dot([mo, Df, mo.T])

                Dfx = np.linalg.multi_dot([mo, Dfx, mo.T])
                Dfy = np.linalg.multi_dot([mo, Dfy, mo.T])
                Dfz = np.linalg.multi_dot([mo, Dfz, mo.T])

                Dxx = np.linalg.multi_dot([mo, Dxx, mo.T])
                Dyy = np.linalg.multi_dot([mo, Dyy, mo.T])
                Dzz = np.linalg.multi_dot([mo, Dzz, mo.T])
                Dxy = np.linalg.multi_dot([mo, Dxy, mo.T])
                Dxz = np.linalg.multi_dot([mo, Dxz, mo.T])
                Dyz = np.linalg.multi_dot([mo, Dyz, mo.T])


                dist_den_1_freq = np.hstack((
                    Dx.real.reshape(-1, 1),
                    Dy.real.reshape(-1, 1),
                    Dz.real.reshape(-1, 1),
                    Df.real.reshape(-1, 1),
                    Dxx.real.reshape(-1, 1),
                    Dyy.real.reshape(-1, 1),
                    Dzz.real.reshape(-1, 1),
                    Dxy.real.reshape(-1, 1),
                    Dxz.real.reshape(-1, 1),
                    Dyz.real.reshape(-1, 1),
                    Dfx.real.reshape(-1, 1),
                    Dfy.real.reshape(-1, 1),
                    Dfz.real.reshape(-1, 1),
                ))

                dist_den_2_freq = np.hstack((
                    Dfxx.real.reshape(-1, 1),
                    Dfyy.real.reshape(-1, 1),
                    Dfzz.real.reshape(-1, 1),
                    Dfxy.real.reshape(-1, 1),
                    Dfxz.real.reshape(-1, 1),
                    Dfyz.real.reshape(-1, 1),
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


    def get_densities(self, freqs, Nx, mo, nocc, norb, Xf):
        """
        Computes the densities needed for the perturbed Fock matrices.

        :param freqs:
            A list of the frequencies
        :param kX:
            A dictonary with all the first-order response matrices
        :param mo:
            A matrix containing the MO coefficents
        :param nocc:
            Number of occupied orbitals

        :return:
            first_order_dens:
             A list of first-order one-time tranformed compounded densities
            second_order_dens:
             A list of first-order two-time tranformed compounded densities
        """

        distributed_density_1 = None
        distributed_density_2 = None
        distributed_density_3 = None


        for w_ind, w in enumerate(freqs):

            nx = ComplexResponse.get_full_solution_vector(Nx[('x', w)])
            ny = ComplexResponse.get_full_solution_vector(Nx[('y', w)])
            nz = ComplexResponse.get_full_solution_vector(Nx[('z', w)])

            Nf = LinearResponseEigenSolver.get_full_solution_vector(Xf[w_ind])

            if self.rank == mpi_master():

                kx = LinearSolver.lrvec2mat(nx, nocc, norb)
                ky = LinearSolver.lrvec2mat(ny, nocc, norb)
                kz = LinearSolver.lrvec2mat(nz, nocc, norb)

                kf = LinearSolver.lrvec2mat(Nf, nocc, norb)

                # create the first order single indexed densiteies #

                Dx = self.commut_mo_density(kx, nocc)
                Dy = self.commut_mo_density(ky, nocc)
                Dz = self.commut_mo_density(kz, nocc)
                Df = self.commut_mo_density(kf, nocc)

                # create the first order two indexed densities #

                Dbfx = self.commut(kx, Df) + self.commut(kf, Dx)
                Dbfy = self.commut(ky, Df) + self.commut(kf, Dy)
                Dbfz = self.commut(kz, Df) + self.commut(kf, Dz)
                Dbc_xx = self.commut(kx, Dx) + self.commut(kx, Dx)
                Dbc_yy = self.commut(ky, Dy) + self.commut(ky, Dy)
                Dbc_zz = self.commut(kz, Dz) + self.commut(kz, Dz)
                Dbc_xy = self.commut(kx, Dy) + self.commut(ky, Dx)
                Dbc_xz = self.commut(kx, Dz) + self.commut(kz, Dx)
                Dbc_yz = self.commut(ky, Dz) + self.commut(kz, Dy)

                # create the first order three indexed densities #
                
                # xx 
                Dbcf_xx = self.commut(kx, Dbfx) # bcd + bdc
                Dbcf_xx += self.commut(kx, Dbfx) # cbd + cdb
                Dbcf_xx += self.commut(kf, Dbc_xx) # dbc + dcb

                # yy
                Dbcf_yy = self.commut(ky, Dbfy) # bcd + bdc
                Dbcf_yy += self.commut(ky, Dbfy) # cbd + cdb
                Dbcf_yy += self.commut(kf, Dbc_yy) # dbc + dcb

                # zz
                Dbcf_zz = self.commut(kz, Dbfz) # bcd + bdc
                Dbcf_zz += self.commut(kz, Dbfz) # cbd + cdb
                Dbcf_zz += self.commut(kf, Dbc_zz) # dbc + dcb

                #  xy
                Dbcf_xy = self.commut(kx, Dbfy) # bcd + bdc
                Dbcf_xy += self.commut(ky, Dbfx) # cbd + cdb
                Dbcf_xy += self.commut(kf, Dbc_xy) # dbc + dcb

                #  xz
                Dbcf_xz = self.commut(kx, Dbfz) # bcd + bdc
                Dbcf_xz += self.commut(kz, Dbfx) # cbd + cdb
                Dbcf_xz += self.commut(kf, Dbc_xz) # dbc + dcb

                #  yz
                Dbcf_yz = self.commut(ky, Dbfz) # bcd + bdc
                Dbcf_yz += self.commut(kz, Dbfy) # cbd + cdb
                Dbcf_yz += self.commut(kf, Dbc_yz) # dbc + dcb

                # Density transformation from MO to AO basis

                # first-order one-index
                Dx = np.linalg.multi_dot([mo, Dx, mo.T])
                Dy = np.linalg.multi_dot([mo, Dy, mo.T])
                Dz = np.linalg.multi_dot([mo, Dz, mo.T])
                Df = np.linalg.multi_dot([mo, Df, mo.T])

                # first-order two-index
                Dbfx = np.linalg.multi_dot([mo, Dbfx, mo.T])
                Dbfy = np.linalg.multi_dot([mo, Dbfy, mo.T])
                Dbfz = np.linalg.multi_dot([mo, Dbfz, mo.T])
                Dbc_xx = np.linalg.multi_dot([mo, Dbc_xx, mo.T])
                Dbc_yy = np.linalg.multi_dot([mo, Dbc_yy, mo.T])
                Dbc_zz = np.linalg.multi_dot([mo, Dbc_zz, mo.T])
                Dbc_xy = np.linalg.multi_dot([mo, Dbc_xy, mo.T])
                Dbc_xz = np.linalg.multi_dot([mo, Dbc_xz, mo.T])
                Dbc_yz = np.linalg.multi_dot([mo, Dbc_yz, mo.T])


                # first-order three-index
                Dbcf_yz = np.linalg.multi_dot([mo, Dbcf_yz, mo.T])
                Dbcf_xz = np.linalg.multi_dot([mo, Dbcf_xz, mo.T])
                Dbcf_xy = np.linalg.multi_dot([mo, Dbcf_xy, mo.T])
                Dbcf_xx = np.linalg.multi_dot([mo, Dbcf_xx, mo.T])
                Dbcf_yy = np.linalg.multi_dot([mo, Dbcf_yy, mo.T])
                Dbcf_zz = np.linalg.multi_dot([mo, Dbcf_zz, mo.T])

           
                dist_den_1_freq = np.hstack((
                    Dx.real.reshape(-1, 1),
                    Dy.real.reshape(-1, 1),
                    Dz.real.reshape(-1, 1),
                    Df.real.reshape(-1, 1),
                ))

                dist_den_2_freq = np.hstack((
                    Dbfx.real.reshape(-1, 1),
                    Dbfy.real.reshape(-1, 1),
                    Dbfz.real.reshape(-1, 1),
                    Dbc_xx.real.reshape(-1, 1),
                    Dbc_yy.real.reshape(-1, 1),
                    Dbc_zz.real.reshape(-1, 1),
                    Dbc_xy.real.reshape(-1, 1),
                    Dbc_xz.real.reshape(-1, 1),
                    Dbc_yz.real.reshape(-1, 1),
                ))

                dist_den_3_freq = np.hstack((
                    Dbcf_xx.real.reshape(-1, 1),
                    Dbcf_yy.real.reshape(-1, 1),
                    Dbcf_zz.real.reshape(-1, 1),
                    Dbcf_xy.real.reshape(-1, 1),
                    Dbcf_xz.real.reshape(-1, 1),
                    Dbcf_yz.real.reshape(-1, 1),
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


    def get_fock_dict_II(self, wi, density_list1, density_list2, F0, mo,molecule, ao_basis, eri_dict, dft_dict):
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

        for w in wi:
            for key in ['Ffxx', 'Ffyy', 'Ffzz', 'Ffxy', 'Ffxz', 'Ffyz']: key_freq_pairs.append((key, w))

        # examine checkpoint file for distributed Focks

        if self.checkpoint_file is not None:
            fpath = Path(self.checkpoint_file)
            fpath = fpath.with_name(fpath.stem)
            fock_file = str(fpath) + '_3patrans_fock_2.h5'
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
                                                 'real', eri_dict,
                                                 dft_dict, density_list1,
                                                 density_list2, None, '3pa_ii')
            else:
                dist_focks = self._comp_nlr_fock(mo, molecule, ao_basis,
                                                 'real', eri_dict,
                                                 None, None, density_list2,
                                                 None, '3pa_ii')

            self._print_fock_time(time.time() - time_start_fock)

            write_distributed_focks(fock_file, dist_focks, key_freq_pairs,
                                    self.comm, self.ostream)

        # assign distributed Focks to key-frequency pairs

        focks = {'F0': F0}

        for fock_index, (key, w) in enumerate(key_freq_pairs):
            if key not in focks:
                focks[key] = {}
            focks[key][w] = DistributedArray(dist_focks.data[:, fock_index],self.comm,distribute=False)

        return focks

    def get_fock_dict(self,
                      wi,
                      first_order_dens,
                      second_order_dens,
                      third_order_dens,
                      F0,
                      mo,
                      molecule,
                      ao_basis,
                      eri_dict,
                      dft_dict=None,
                      profiler=None):
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

        for wb in wi:
            for key in ['Fbfx', 
                        'Fbfy', 
                        'Fbfz',
                        'Fbc_xx',
                        'Fbc_yy',
                        'Fbc_zz',
                        'Fbc_xy',
                        'Fbc_xz',
                        'Fbc_yz',
                        ]:
                
                key_freq_pairs.append((key, wb))


        for wb in wi:
            for key in ['Fbcf_xx',
                        'Fbcf_yy',
                        'Fbcf_zz',
                        'Fbcf_xy',
                        'Fbcf_xz',
                        'Fbcf_yz',
                        ]:
                
                key_freq_pairs.append((key, wb))

        # examine checkpoint for distributed Focks

        if self.checkpoint_file is not None:
            fpath = Path(self.checkpoint_file)
            fpath = fpath.with_name(fpath.stem)
            fock_file = str(fpath) + '_3patrans_fock_1.h5'
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
                                                 'real', eri_dict,
                                                 dft_dict, first_order_dens,
                                             second_order_dens, third_order_dens,
                                             '3pa',profiler)
            else:
                density_list_23 = DistributedArray(second_order_dens.data,
                                                   self.comm,
                                                   distribute=False)
                
                density_list_23.append(third_order_dens, axis=1)
                dist_focks = self._comp_nlr_fock(mo, molecule, ao_basis,
                                                 'real', eri_dict,
                                                 None, None, None,
                                                 density_list_23, '3pa')

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

    def get_e3(self, wi, Nx, N_xy, fo, fo2, fo3, nocc, norb, Xf):
        """
        Contracts E[3]

        :param wi:
            A list of freqs
        :param kX:
            A dict of the single index response matricies
        :param fo:
            A dictonary of transformed Fock matricies from fock_dict
        :param fo2:
            A dictonarty of transfromed Fock matricies from subspace of response solver
        :param nocc:
            The number of occupied orbitals
        :param norb:
            The total number of orbitals

        :return:
            A dictonary of compounded E[3] tensors for the isotropic cubic
            response function for QRF
        """

        e3vec = {}

        for w_ind, w in enumerate(wi):
            vec_pack = np.array([
                fo2['Ffxx'][w].data,
                fo2['Ffyy'][w].data,
                fo2['Ffzz'][w].data,
                fo2['Ffxy'][w].data,
                fo2['Ffxz'][w].data,
                fo2['Ffyz'][w].data,
                fo[('x', w)].data,
                fo[('y', w)].data,
                fo[('z', w)].data,
                fo3[('fx', - 2 * w)].data,
                fo3[('fy', - 2 * w)].data,
                fo3[('fz', - 2 * w)].data,
                fo3[('xx',  2 * w)].data,
                fo3[('yy',  2 * w)].data,
                fo3[('zz',  2 * w)].data,
                fo3[('xy',  2 * w)].data,
                fo3[('xz',  2 * w)].data,
                fo3[('yz',  2 * w)].data,
                fo[w_ind].data,
            ]).T.copy()

            vec_pack = self._collect_vectors_in_columns(vec_pack)

            # First-order response vectors
            Nbx = ComplexResponse.get_full_solution_vector(Nx[('x', w)])
            Nby = ComplexResponse.get_full_solution_vector(Nx[('y', w)])
            Nbz = ComplexResponse.get_full_solution_vector(Nx[('z', w)])
            Nf = LinearResponseEigenSolver.get_full_solution_vector(Xf[w_ind])

            # Second-order response vectors
            Nfx = ComplexResponse.get_full_solution_vector(N_xy[('fx', - 2 * w)])
            Nfy = ComplexResponse.get_full_solution_vector(N_xy[('fy', - 2 * w)])
            Nfz = ComplexResponse.get_full_solution_vector(N_xy[('fz', - 2 * w)])

            Nxx = ComplexResponse.get_full_solution_vector(N_xy[('xx', 2 * w)])
            Nyy = ComplexResponse.get_full_solution_vector(N_xy[('yy', 2 * w)])
            Nzz = ComplexResponse.get_full_solution_vector(N_xy[('zz', 2 * w)])
            Nxy = ComplexResponse.get_full_solution_vector(N_xy[('xy', 2 * w)])
            Nxz = ComplexResponse.get_full_solution_vector(N_xy[('xz', 2 * w)])
            Nyz = ComplexResponse.get_full_solution_vector(N_xy[('yz', 2 * w)])

            if self.rank != mpi_master():
                continue
            
            # Fock matrices
            vec_pack = vec_pack.T.copy().reshape(-1, norb, norb)

            (ffxx, ffyy, ffzz, ffxy, ffxz, ffyz,  fx, fy, fz, ffx, ffy, ffz, fxx,fyy,fzz,fxy,fxz,fyz, ff) = vec_pack

            fx = np.conjugate(fx).T
            fy = np.conjugate(fy).T
            fz = np.conjugate(fz).T

            # check these terms for factor of -1 / np.sqrt(2)
            ffx = np.conjugate(ffx).T 
            ffy = np.conjugate(ffy).T 
            ffz = np.conjugate(ffz).T 

            fxx = np.conjugate(fxx).T
            fyy = np.conjugate(fyy).T
            fzz = np.conjugate(fzz).T
            fxy = np.conjugate(fxy).T
            fxz = np.conjugate(fxz).T
            fyz = np.conjugate(fyz).T

            ff = np.conjugate(ff).T * -1 / np.sqrt(2)

            F0_a = fo2['F0']

            # Response matrices
            # First-order response matrices
            kx = (LinearSolver.lrvec2mat(Nbx, nocc, norb)).T
            ky = (LinearSolver.lrvec2mat(Nby, nocc, norb)).T
            kz = (LinearSolver.lrvec2mat(Nbz, nocc, norb)).T
            kf = (LinearSolver.lrvec2mat(Nf, nocc, norb)).T
            
            # Second-order response matrices
            kfx = (LinearSolver.lrvec2mat(Nfx, nocc, norb)).T 
            kfy = (LinearSolver.lrvec2mat(Nfy, nocc, norb)).T
            kfz = (LinearSolver.lrvec2mat(Nfz, nocc, norb)).T

            kxx = (LinearSolver.lrvec2mat(Nxx, nocc, norb)).T  
            kyy = (LinearSolver.lrvec2mat(Nyy, nocc, norb)).T
            kzz = (LinearSolver.lrvec2mat(Nzz, nocc, norb)).T
            kxy = (LinearSolver.lrvec2mat(Nxy, nocc, norb)).T
            kxz = (LinearSolver.lrvec2mat(Nxz, nocc, norb)).T
            kyz = (LinearSolver.lrvec2mat(Nyz, nocc, norb)).T

            # XF xi terms
            xi_xfx = self._xi(kx, kfx, fx, ffx, F0_a)
            xi_xfy = self._xi(kx, kfy, fx, ffy, F0_a)
            xi_xfz = self._xi(kx, kfz, fx, ffz, F0_a)
            xi_yfx = self._xi(ky, kfx, fy, ffx, F0_a)
            xi_yfy = self._xi(ky, kfy, fy, ffy, F0_a)
            xi_yfz = self._xi(ky, kfz, fy, ffz, F0_a)
            xi_zfx = self._xi(kz, kfx, fz, ffx, F0_a)
            xi_zfy = self._xi(kz, kfy, fz, ffy, F0_a)
            xi_zfz = self._xi(kz, kfz, fz, ffz, F0_a)

            # XX xi terms
            xi_fxx = self._xi(kf, kxx, ff, fxx, F0_a)
            xi_fyy = self._xi(kf, kyy, ff, fyy, F0_a)
            xi_fzz = self._xi(kf, kzz, ff, fzz, F0_a)
            xi_fxy = self._xi(kf, kxy, ff, fxy, F0_a)
            xi_fxz = self._xi(kf, kxz, ff, fxz, F0_a)
            xi_fyz = self._xi(kf, kyz, ff, fyz, F0_a)

            # xx 
            e3fock = xi_xfx.T  + xi_xfx.T + xi_fxx.T   +   0.5 * ffxx.T
            e3vec[('xx', w)] = self.anti_sym(-2 * LinearSolver.lrmat2vec(e3fock, nocc, norb))

            # yy 
            e3fock = xi_yfy.T  + xi_yfy.T + xi_fyy.T  +   0.5 * ffyy.T
            e3vec[('yy', w)] = self.anti_sym(-2 * LinearSolver.lrmat2vec(e3fock, nocc, norb))

            # zz
            e3fock = xi_zfz.T  + xi_zfz.T + xi_fzz.T  +   0.5 * ffzz.T
            e3vec[('zz', w)] = self.anti_sym(-2 * LinearSolver.lrmat2vec(e3fock, nocc, norb))

            # xy
            e3fock = xi_xfy.T  + xi_yfx.T + xi_fxy.T  +   0.5 * ffxy.T
            e3vec[('xy', w)] = self.anti_sym(-2 * LinearSolver.lrmat2vec(e3fock, nocc, norb))
            e3vec[('yx', w)] = self.anti_sym(-2 * LinearSolver.lrmat2vec(e3fock, nocc, norb))

            # xz
            e3fock = xi_xfz.T  + xi_zfx.T + xi_fxz.T  +   0.5 * ffxz.T
            e3vec[('xz', w)] = self.anti_sym(-2 * LinearSolver.lrmat2vec(e3fock, nocc, norb))
            e3vec[('zx', w)] = self.anti_sym(-2 * LinearSolver.lrmat2vec(e3fock, nocc, norb))

            # yz
            e3fock = xi_yfz.T  + xi_zfy.T + xi_fyz.T  +   0.5 * ffyz.T
            e3vec[('yz', w)] = self.anti_sym(-2 * LinearSolver.lrmat2vec(e3fock, nocc, norb))
            e3vec[('zy', w)] = self.anti_sym(-2 * LinearSolver.lrmat2vec(e3fock, nocc, norb))

        return e3vec

    def print_header(self):
        """
        Prints QRF setup header to output stream.
        """

        self.ostream.print_blank()

        title = 'Quadratic Response Driver Setup'
        self.ostream.print_header(title)
        self.ostream.print_header('=' * (len(title) + 2))
        self.ostream.print_blank()

        width = 50

        cur_str = 'ERI Screening Threshold         : {:.1e}'.format(
            self.eri_thresh)
        self.ostream.print_header(cur_str.ljust(width))
        cur_str = 'Convergance Threshold           : {:.1e}'.format(
            self.conv_thresh)
        self.ostream.print_header(cur_str.ljust(width))
        cur_str = 'Max. Number of Iterations       : {:d}'.format(self.max_iter)
        self.ostream.print_header(cur_str.ljust(width))
        self.ostream.print_header(cur_str.ljust(width))

        self.ostream.print_blank()
        self.ostream.flush()

    def _print_component(self, value, width):
        """
        Prints QRF component.

        :param value:
            The complex value
        :param width:
            The width for the output
        """

        w_str = '{:>12s}{:20.9f} {:20.9f} {:20.9f}'.format(
            'x', value[0][0].real, value[0][1].real, value[0][2].real)
        self.ostream.print_header(w_str.ljust(width))

        w_str = '{:>12s}{:20.9f} {:20.9f} {:20.9f}'.format(
            'y', value[1][0].real, value[1][1].real, value[1][2].real)
        self.ostream.print_header(w_str.ljust(width))

        w_str = '{:>12s}{:20.9f} {:20.9f} {:20.9f}'.format(
            'z', value[2][0].real, value[2][1].real, value[2][2].real)
        self.ostream.print_header(w_str.ljust(width))

    def _print_results(self, rsp_results):
        """
        Prints the results of 3PA calculation for a rank-3 tensor in three separate tables.

        :param rsp_results:
            A dictionary containing the results of the response calculation.
        """

        width = 82
        freqs = rsp_results['photon_energies']
        T_tensors = rsp_results['transition_moments']

        self.ostream.print_header('Three-Photon Absorption (3PA) Transition Moments (a.u.)')
        
        # Extract unique frequency values directly from stored tensor keys
        unique_freqs = sorted(set(w for _, w in T_tensors.keys()), reverse=True)

        xyz_triples = [f'{x}{y}{z}'
                       for i, x in enumerate('xyz')
                       for j, y in enumerate('xyz')
                       for k, z in enumerate('xyz')
                       if (i <= j and j <= k)]

        tensor_labels_groups = [xyz_triples[:5], xyz_triples[5:]]

        for tensor_labels in tensor_labels_groups:
            self.ostream.print_blank()

            # Create header with proper spacing
            header_format = '{:<7s} {:>11s}    ' + ' '.join(
                ['{:>11s}' for x in range(len(tensor_labels))])
            header_str = header_format.format('State', 'Ex.Energy', *tensor_labels)
            self.ostream.print_header(header_str)

            # Print a separator line with correct spacing
            self.ostream.print_header('-' * len(header_str))

            for w_ind, w in enumerate(unique_freqs):
                row_values = [f'  {w_ind + 1:<5d}', f'{- 3 * w * hartree_in_ev():>11.5f} eV']

                for label in tensor_labels:
                    tensor_key = (label, w)

                    if tensor_key in T_tensors:
                        row_values.append(f'{T_tensors[tensor_key].real:11.4f}')
                    else:
                        row_values.append(f'{"N/A":>11s}')

                self.ostream.print_header(' '.join(row_values))  # Print row data

        tpa_strengths = rsp_results['3pa_strengths']
        #tpa_cross_sections = rsp_results['cross_sections']

        self.ostream.print_blank()
        self.ostream.print_blank()

        #title = '3PA Strength and Cross-Section (Linear Polarization)'
        title = '3PA Strength (Linear Polarization)'
        self.ostream.print_header(title)
        self.ostream.print_header('-' * width)

        title = '{:<12s}{:>11s}{:>27s}'.format(
            'State', 'Ex.Energy', '3PA strength    ')
        self.ostream.print_header(title)
        self.ostream.print_header('-' * width)

        for w_ind, w in enumerate(freqs):
            exec_str = '  {:<4d}'.format(w_ind + 1)
            exec_str += '{:15.6f} eV'.format(3 * w * hartree_in_ev())
            exec_str += '{:20.6f} a.u.'.format(tpa_strengths['linear'][-w])
            #exec_str += '{:20.6f} GM'.format(tpa_cross_sections['linear'][-w])
            self.ostream.print_header(exec_str)

        self.ostream.print_blank()
        self.ostream.print_blank()

        #title = '3PA Strength and Cross-Section (Circular Polarization)'
        title = '3PA Strength (Circular Polarization)'
        self.ostream.print_header(title)
        self.ostream.print_header('-' * width)

        title = '{:<12s}{:>11s}{:>27s}'.format(
            'State', 'Ex.Energy', '3PA strength    ')
        self.ostream.print_header(title)
        self.ostream.print_header('-' * width)

        for w_ind, w in enumerate(freqs):
            exec_str = '  {:<4d}'.format(w_ind + 1)
            exec_str += '{:15.6f} eV'.format(3 * w * hartree_in_ev())
            exec_str += '{:20.6f} a.u.'.format(tpa_strengths['circular'][-w])
            #exec_str += '{:20.6f} GM'.format(tpa_cross_sections['circular'][-w])
            self.ostream.print_header(exec_str)
        self.ostream.print_blank()

        self.ostream.flush()
