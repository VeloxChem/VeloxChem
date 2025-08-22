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
from .veloxchemlib import mpi_master, hartree_in_ev
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


class DoubleResBetaDriver(NonlinearSolver):
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
        self.damping = 0.0

        # tpa transition settings

        self.initial_state = 3
        self.final_state = 3
        self.nstates = max(self.initial_state, self.final_state)



        # input keywords
        self._input_keywords['response'].update({
            'initial_state': ('int', 'index state'),
            'final_state': ('int', 'index of state'),
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
            'DoubleQuadResDriver: not implemented for unrestricted case')

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
            'nstates','norm_thresh', 'lindep_thresh', 'conv_thresh',
            'max_iter', 'eri_thresh', 'timing', 'memory_profiling',
            'batch_size', 'restart', 'xcfun', 'grid_level', 'potfile',
            'electric_field', 'program_end_time', '_debug', '_block_size_factor',
            'ri_coulomb'
        ]
        
        self.nstates = max(self.initial_state, self.final_state)

        for key in rpa_keywords:
            setattr(rpa_drv, key, getattr(self, key))

        if self.checkpoint_file is not None:
            fpath = Path(self.checkpoint_file)
            fpath = fpath.with_name(fpath.stem)
            rpa_drv.checkpoint_file = str(fpath) + '_doublequadres_rpa.h5'

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

        freqs = [rpa_results['eigenvalues'][self.final_state - 1] - rpa_results['eigenvalues'][self.initial_state - 1]]
        freqs_for_response_vectors = freqs + [0.0]

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
            'electric_field', 'program_end_time', '_debug', '_block_size_factor',
            'ri_coulomb'
        }

        for key in cpp_keywords:
            setattr(N_drv, key, getattr(self, key))

        if self.checkpoint_file is not None:
            fpath = Path(self.checkpoint_file)
            fpath = fpath.with_name(fpath.stem)
            N_drv.checkpoint_file = str(fpath) + '_doublequadres_cpp.h5'

        N_results = N_drv.compute(molecule, ao_basis, scf_results, B)

        self._is_converged = N_drv.is_converged

        Nx = N_results['solutions']
        Focks = N_results['focks']
        Focks.update(rpa_results['focks'])

        profiler.check_memory_usage('CPP')

        ret_dict = self.compute_quad_components(Focks, freqs, X, d_a_mo, Nx,
                                                scf_results, molecule, ao_basis,
                                                profiler, Xf)

        valstr = '*** Time spent in double residue of quadratic response calculation: '
        valstr += '{:.2f} sec ***'.format(time.time() - start_time)
        self.ostream.print_header(valstr)
        self.ostream.print_blank()
        self.ostream.flush()

        profiler.end(self.ostream)

        if self.rank == mpi_master():
            ret_dict.update({
                'oscillator_strengths': oscillator_strengths,
                'elec_trans_dipoles': elec_trans_dipoles,
                'excitation_details': excitation_details
            })

        return ret_dict


    def compute_quad_components(self, Focks, freqs, X, d_a_mo, Nx, scf_results,
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
        first_order_dens, second_order_dens = self.get_densities(
            freqs, Nx, mo, nocc, norb, Xf)

        profiler.check_memory_usage('Densities')

        # computing the compounded first-order Fock matrices
        fock_dict = self.get_fock_dict(freqs, first_order_dens,
                                       second_order_dens, F0, mo, molecule,
                                       ao_basis, eri_dict, dft_dict, profiler)

        profiler.check_memory_usage('Focks')

        e3_dict = self.get_e3(freqs, Nx, fock_dict, Focks, nocc, norb, Xf)

        profiler.check_memory_usage('E[3]')

        ret_dict = {}
        M_tensors = {}
        excited_state_dipole_moments = {}

        # Compute dipole vector
        scf_prop = FirstOrderProperties(self.comm, self.ostream)
        scf_prop.compute_scf_prop(molecule, ao_basis, scf_results)
    
        N0_x = ComplexResponse.get_full_solution_vector(Nx[('x', freqs[0])])
        N0_y = ComplexResponse.get_full_solution_vector(Nx[('y', freqs[0])])
        N0_z = ComplexResponse.get_full_solution_vector(Nx[('z', freqs[0])])


        Nf = LinearResponseEigenSolver.get_full_solution_vector(Xf[self.initial_state - 1])
        Ng_ = LinearResponseEigenSolver.get_full_solution_vector(Xf[self.final_state - 1])

        if self.rank == mpi_master():

            Nax0 = self.flip_yz(N0_x)
            Nay0 = self.flip_yz(N0_y)
            Naz0 = self.flip_yz(N0_z)

            Ng_ = self.flip_yz(Ng_)

            op_x = X['x']
            op_y = X['y']
            op_z = X['z']

            kf = LinearSolver.lrvec2mat(Nf, nocc, norb)
            kg_ = LinearSolver.lrvec2mat(Ng_, nocc, norb)

            # Double residue x-component

            ground_state_dipole_x = scf_prop.get_property('dipole moment')[0]

            A2Nc = self._a2_contract(kf, op_x, d_a_mo, nocc, norb)
            A2Nc_star = self._a2_contract(kg_, op_x, d_a_mo, nocc, norb)

            NaE3NbNc = np.dot(Nax0.T, e3_dict['fg'].real)

            Nc_starA2Nc = np.dot(Ng_.T, A2Nc)
            NcA2Nc_star = np.dot(Nf.T, A2Nc_star)

            val_A2 = -(Nc_starA2Nc + NcA2Nc_star)
            val_E3 = NaE3NbNc

            #excited_state_dipole_moments.update({('x',self.initial_state,self.final_state): (val_E3 + val_A2 + ground_state_dipole_x).real})
            excited_state_dipole_moments.update({('x',self.initial_state,self.final_state): (val_E3 + val_A2).real})

            # Double residue y-component

            ground_state_dipole_y = scf_prop.get_property('dipole moment')[1]

            A2Nc = self._a2_contract(kf, op_y, d_a_mo, nocc, norb)
            A2Nc_star = self._a2_contract(kg_, op_y, d_a_mo, nocc, norb)

            NaE3NbNc = np.dot(Nay0.T, e3_dict['fg'].real)

            Nc_starA2Nc = np.dot(Ng_.T, A2Nc)
            NcA2Nc_star = np.dot(Nf.T, A2Nc_star)

            val_A2 = -(Nc_starA2Nc + NcA2Nc_star)
            val_E3 = NaE3NbNc

            #excited_state_dipole_moments.update({('y',self.initial_state,self.final_state): (val_E3 + val_A2 + ground_state_dipole_y).real})
            excited_state_dipole_moments.update({('y',self.initial_state,self.final_state): (val_E3 + val_A2).real})

            # Double residue z-component

            ground_state_dipole_z = scf_prop.get_property('dipole moment')[2]

            A2Nc = self._a2_contract(kf, op_z, d_a_mo, nocc, norb)
            A2Nc_star = self._a2_contract(kg_, op_z, d_a_mo, nocc, norb)

            NaE3NbNc = np.dot(Naz0.T, e3_dict['fg'].real)

            Nc_starA2Nc = np.dot(Ng_.T, A2Nc)
            NcA2Nc_star = np.dot(Nf.T, A2Nc_star)

            val_A2 = -(Nc_starA2Nc + NcA2Nc_star)
            val_E3 = NaE3NbNc

            #excited_state_dipole_moments.update({('z',self.initial_state,self.final_state): (val_E3 + val_A2 + ground_state_dipole_z).real})
            excited_state_dipole_moments.update({('z',self.initial_state,self.final_state): (val_E3 + val_A2 ).real})
                

        self.ostream.print_blank()
        w_str = 'Summary of Double residue of Quadratic Response Function'
        self.ostream.print_header(w_str)
        self.ostream.print_header('=' * (len(w_str) + 2))
        self.ostream.print_blank()


        self._print_transition_dipoles("Transition dipole moments",excited_state_dipole_moments)
        if self.rank == mpi_master():

            profiler.check_memory_usage('End of QRF')

            ret_dict = {
                'photon_energies': [-w for w in freqs],
                'transition_moments': M_tensors,
                'excited_state_dipole_moments': excited_state_dipole_moments,
                'ground_state_dipole_moments':
                    scf_prop.get_property('dipole moment')
            }

            #self._print_results(ret_dict)

            return ret_dict
        else:
            return None


    def _print_transition_dipoles(self, title, trans_dipoles):
        """
        Prints transition dipole moments to output stream using bra-ket notation.

        :param title:
            The title to be printed to the output stream.
        :param trans_dipoles:
            Dictionary of transition dipole moments with keys ('x'/'y'/'z', initial_state, final_state)
            and values as the dipole moment components.
        """
        # Group dipole moments by (initial_state, final_state)
        dipole_dict = {}
        for (component, initial_state, final_state), value in trans_dipoles.items():
            key = (initial_state, final_state)
            if key not in dipole_dict:
                dipole_dict[key] = {'x': 0.0, 'y': 0.0, 'z': 0.0}
            dipole_dict[key][component.lower()] = value

        # Print header
        valstr = title
        self.ostream.print_header(valstr.ljust(80))
        self.ostream.print_header(('-' * len(valstr)).ljust(80))
        valstr = ' ' * 20 + '{:>12s}{:>12s}{:>12s}'.format('X', 'Y', 'Z')
        self.ostream.print_header(valstr.ljust(80))

        # Print dipole moments for each state transition using bra-ket notation
        for (initial_state, final_state), components in dipole_dict.items():
            label = f'<{initial_state}|mu|{final_state}> :'
            valstr = f'{label:<20}' + \
                    '{:>12.6f}{:>12.6f}{:>12.6f}'.format(
                        components['x'], components['y'], components['z']
                    )
            self.ostream.print_header(valstr.ljust(80))

        self.ostream.print_blank()
        self.ostream.flush()

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

        Nf = LinearResponseEigenSolver.get_full_solution_vector(Xf[self.initial_state - 1])
        Ng = LinearResponseEigenSolver.get_full_solution_vector(Xf[self.final_state - 1])

        if self.rank == mpi_master():

            kf = LinearSolver.lrvec2mat(Nf, nocc, norb).real
            kg_ = LinearSolver.lrvec2mat(self.flip_yz(Ng), nocc, norb).real

            # create the first order single indexed densiteies #

            Df = self.commut_mo_density(kf, nocc).real
            Dg_ = self.commut_mo_density(kg_, nocc).real

            # create the first order two indexed densities #

            Dfg_ = self.commut(kg_, Df) + self.commut(kf, Dg_)

            # Density transformation from MO to AO basis

            Df = np.linalg.multi_dot([mo, Df, mo.T]).real
            Dg_ = np.linalg.multi_dot([mo, Dg_, mo.T]).real

            Dfg_ = np.linalg.multi_dot([mo, Dfg_, mo.T]).real

            distributed_density_1_freq = np.hstack((
                Df.real.reshape(-1, 1),
                Df.imag.reshape(-1, 1),
                Dg_.real.reshape(-1, 1),
                Dg_.imag.reshape(-1, 1),
            ))

            distributed_density_2_freq = np.hstack((
                Dfg_.real.reshape(-1, 1),
                Dfg_.imag.reshape(-1, 1),
            ))

        else:
            distributed_density_1_freq = None
            distributed_density_2_freq = None

        distributed_density_1_freq = DistributedArray(
            distributed_density_1_freq, self.comm)
        distributed_density_2_freq = DistributedArray(
            distributed_density_2_freq, self.comm)

        if distributed_density_1 is None:
            distributed_density_1 = DistributedArray(
                distributed_density_1_freq.data,
                self.comm,
                distribute=False)
        else:
            distributed_density_1.append(distributed_density_1_freq, axis=1)

        if distributed_density_2 is None:
            distributed_density_2 = DistributedArray(
                distributed_density_2_freq.data,
                self.comm,
                distribute=False)
        else:
            distributed_density_2.append(distributed_density_2_freq, axis=1)

        return distributed_density_1, distributed_density_2

    def get_fock_dict(self,
                      wi,
                      first_order_dens,
                      second_order_dens,
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
            for key in ['F(fg_)']:
                key_freq_pairs.append((key, wb))

        # examine checkpoint for distributed Focks

        if self.checkpoint_file is not None:
            fpath = Path(self.checkpoint_file)
            fpath = fpath.with_name(fpath.stem)
            fock_file = str(fpath) + '_doublequadres_fock.h5'
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

            dist_focks = self._comp_nlr_fock(mo, molecule, ao_basis, 'real_and_imag',
                                             eri_dict, dft_dict,
                                             first_order_dens,
                                             second_order_dens, None,
                                             'qrf', profiler)

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

    def get_e3(self, wi, Nx, fo, fo2, nocc, norb, Xf):
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

    
        vec_pack = np.array([
            fo['F(fg_)'][wi[0]].data,
            fo2[self.initial_state - 1].data,
            fo2[self.final_state - 1].data,
        ]).T.copy()

        vec_pack = self._collect_vectors_in_columns(vec_pack)

        Nf = LinearResponseEigenSolver.get_full_solution_vector(Xf[self.initial_state - 1])
        Ng_ = LinearResponseEigenSolver.get_full_solution_vector(Xf[self.final_state - 1])

        if self.rank == mpi_master():
    
            vec_pack = vec_pack.T.copy().reshape(-1, norb, norb)

            (Ffg_, Ff, Fg_) = vec_pack

            ff = np.conjugate(Ff).T * -1 / np.sqrt(2)

            fg_star = np.conjugate(Fg_).T * -1 / np.sqrt(2)
            fg_star = np.conjugate(fg_star).T


            F0_a = fo['F0']

            # Response

            kf = (LinearSolver.lrvec2mat(Nf, nocc, norb)).T
            kg_star = (LinearSolver.lrvec2mat(self.flip_yz(Ng_), nocc, norb)).T

            # fg
            xi = self._xi(kg_star, kf, fg_star, ff, F0_a)
            e3fock = xi.T + 0.5 * Ffg_.T
            e3vec['fg'] = self.anti_sym(-2 * LinearSolver.lrmat2vec(e3fock, nocc, norb))


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

        w_str = '{:>12s}{:20.5f} {:20.5f} {:20.5f}'.format(
            'x', value[0][0].real, value[0][1].real, value[0][2].real)
        self.ostream.print_header(w_str.ljust(width))

        w_str = '{:>12s}{:20.5f} {:20.5f} {:20.5f}'.format(
            'y', value[1][0].real, value[1][1].real, value[1][2].real)
        self.ostream.print_header(w_str.ljust(width))

        w_str = '{:>12s}{:20.5f} {:20.5f} {:20.5f}'.format(
            'z', value[2][0].real, value[2][1].real, value[2][2].real)
        self.ostream.print_header(w_str.ljust(width))

    def _print_results(self, rsp_results):
        """
        Prints the results of TPA calculation.

        :param rsp_results:
            A dictonary containing the results of response calculation.
        """

        width = 92

        freqs = rsp_results['photon_energies']
        M_tensors = rsp_results['transition_moments']

        self.ostream.print_blank()
        title = 'Components of TPA Transition Moments (a.u.)'
        self.ostream.print_header(title)
        self.ostream.print_header('-' * width)

        title = '  {:<9s} {:>12s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s} '.format(
            'Ex. State', 'Ex. Energy', 'Sxx  ', 'Syy  ', 'Szz  ', 'Sxy  ',
            'Sxz  ', 'Syz  ')
        self.ostream.print_header(title.ljust(width))
        self.ostream.print_header('-' * width)

        for w_ind, w in enumerate(freqs):
            exec_str = '{:7d}   '.format(w_ind + 1)
            exec_str += '{:11.6f} eV'.format(w * hartree_in_ev())
            exec_str += '{:11.4f}'.format(M_tensors[-w][0, 0].real)
            exec_str += '{:11.4f}'.format(M_tensors[-w][1, 1].real)
            exec_str += '{:11.4f}'.format(M_tensors[-w][2, 2].real)
            exec_str += '{:11.4f}'.format(M_tensors[-w][0, 1].real)
            exec_str += '{:11.4f}'.format(M_tensors[-w][0, 2].real)
            exec_str += '{:11.4f}'.format(M_tensors[-w][1, 2].real)
            self.ostream.print_header(exec_str.ljust(width))
        self.ostream.print_blank()
        self.ostream.print_blank()

        tpa_strengths = rsp_results['tpa_strengths']
        tpa_cross_sections = rsp_results['cross_sections']

        title = 'TPA Strength and Cross-Section (Linear Polarization)'
        self.ostream.print_header(title)
        self.ostream.print_header('-' * width)

        title = '  {:<9s} {:>12s}{:>28s}{:>28s}'.format(
            'Ex. State', 'Ex. Energy', 'TPA strength    ',
            'TPA cross-section    ')
        self.ostream.print_header(title.ljust(width))
        self.ostream.print_header('-' * width)

        for w_ind, w in enumerate(freqs):
            exec_str = '{:7d}   '.format(w_ind + 1)
            exec_str += '{:11.6f} eV'.format(w * hartree_in_ev())
            exec_str += '{:20.6f} a.u.'.format(tpa_strengths['linear'][-w])
            exec_str += '{:20.6f} GM'.format(tpa_cross_sections['linear'][-w])
            self.ostream.print_header(exec_str.ljust(width))
        self.ostream.print_blank()
        self.ostream.print_blank()

        title = 'TPA Strength and Cross-Section (Circular Polarization)'
        self.ostream.print_header(title)
        self.ostream.print_header('-' * width)

        title = '  {:<9s} {:>12s}{:>28s}{:>28s}'.format(
            'Ex. State', 'Ex. Energy', 'TPA strength    ',
            'TPA cross-section    ')
        self.ostream.print_header(title.ljust(width))
        self.ostream.print_header('-' * width)

        for w_ind, w in enumerate(freqs):
            exec_str = '{:7d}   '.format(w_ind + 1)
            exec_str += '{:11.6f} eV'.format(w * hartree_in_ev())
            exec_str += '{:20.6f} a.u.'.format(tpa_strengths['circular'][-w])
            exec_str += '{:20.6f} GM'.format(tpa_cross_sections['circular'][-w])
            self.ostream.print_header(exec_str.ljust(width))
        self.ostream.print_blank()
        self.ostream.print_blank()

        self.ostream.flush()


