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
from .veloxchemlib import (mpi_master, bohr_in_angstroms, hartree_in_ev,
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


class TpaTransitionDriver(NonlinearSolver):
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
        dipole_drv = ElectricDipoleIntegralsDriver(self.comm)
        dipole_mats = dipole_drv.compute(molecule, ao_basis)

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

            b_grad = list(b_grad)
            for ind in range(len(b_grad)):
                b_grad[ind] *= inv_sqrt_2

        rpa_drv = LinearResponseEigenSolver(self.comm, self.ostream)
        rpa_drv.nonlinear = True

        rpa_keywords = [
            'nstates', 'norm_thresh', 'lindep_thresh', 'conv_thresh',
            'max_iter', 'eri_thresh', 'qq_type', 'timing', 'memory_profiling',
            'batch_size', 'restart', 'xcfun', 'grid_level', 'potfile',
            'electric_field', 'program_end_time'
        ]

        for key in rpa_keywords:
            setattr(rpa_drv, key, getattr(self, key))

        if self.checkpoint_file is not None:
            rpa_drv.checkpoint_file = str(
                Path(self.checkpoint_file).with_suffix('.tpatrans_rpa.h5'))

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

        freqs = [-0.5 * a for a in rpa_results['eigenvalues']]
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
                'x': 2 * self.ao2mo(mo, dipole_mats.x_to_numpy()),
                'y': 2 * self.ao2mo(mo, dipole_mats.y_to_numpy()),
                'z': 2 * self.ao2mo(mo, dipole_mats.z_to_numpy())
            }
        else:
            X = None

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
                Path(self.checkpoint_file).with_suffix('.tpatrans_cpp.h5'))

        N_results = N_drv.compute(molecule, ao_basis, scf_results, B)

        self._is_converged = N_drv.is_converged

        Nx = N_results['solutions']
        Focks = N_results['focks']
        Focks.update(rpa_results['focks'])

        profiler.check_memory_usage('CPP')

        ret_dict = self.compute_quad_components(Focks, freqs, X, d_a_mo, Nx,
                                                scf_results, molecule, ao_basis,
                                                profiler, Xf)

        valstr = '*** Time spent in quadratic response calculation: '
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

        dft_dict = self._init_dft(molecule, scf_results)

        # computing all compounded first-order densities
        first_order_dens, second_order_dens = self.get_densities(
            freqs, Nx, mo, nocc, norb, Xf)

        profiler.check_memory_usage('Densities')

        # computing the compounded first-order Fock matrices
        fock_dict = self.get_fock_dict(freqs, first_order_dens,
                                       second_order_dens, F0, mo, molecule,
                                       ao_basis, dft_dict, profiler)

        profiler.check_memory_usage('Focks')

        e3_dict = self.get_e3(freqs, Nx, fock_dict, Focks, nocc, norb, Xf)

        profiler.check_memory_usage('E[3]')

        ret_dict = {}
        M_tensors = {}
        excited_state_dipole_moments = {}

        # Compute dipole vector
        scf_prop = FirstOrderProperties(self.comm, self.ostream)
        scf_prop.compute_scf_prop(molecule, ao_basis, scf_results)

        for w_ind, w in enumerate(freqs):
            m = np.zeros((3, 3), dtype='complex128')

            N0_x = ComplexResponse.get_full_solution_vector(Nx[('x', 0)])
            N0_y = ComplexResponse.get_full_solution_vector(Nx[('y', 0)])
            N0_z = ComplexResponse.get_full_solution_vector(Nx[('z', 0)])

            N_x = ComplexResponse.get_full_solution_vector(Nx[('x', w)])
            N_y = ComplexResponse.get_full_solution_vector(Nx[('y', w)])
            N_z = ComplexResponse.get_full_solution_vector(Nx[('z', w)])

            Nc = LinearResponseEigenSolver.get_full_solution_vector(Xf[w_ind])

            if self.rank == mpi_master():

                Nax = self.flip_yz(N_x)
                Nay = self.flip_yz(N_y)
                Naz = self.flip_yz(N_z)

                Nax0 = self.flip_yz(N0_x)
                Nay0 = self.flip_yz(N0_y)
                Naz0 = self.flip_yz(N0_z)

                Nc_star = self.flip_yz(Nc)

                op_x = X['x']
                op_y = X['y']
                op_z = X['z']

                kx = LinearSolver.lrvec2mat(N_x, nocc, norb)
                ky = LinearSolver.lrvec2mat(N_y, nocc, norb)
                kz = LinearSolver.lrvec2mat(N_z, nocc, norb)
                kc = LinearSolver.lrvec2mat(Nc, nocc, norb)
                kc_ = LinearSolver.lrvec2mat(Nc_star, nocc, norb)

                # Double residue x-component

                ground_state_dipole_x = scf_prop.get_property(
                    'dipole moment')[0]

                A2Nc = self._a2_contract(kc, op_x, d_a_mo, nocc, norb)
                A2Nc_star = self._a2_contract(kc_, op_x, d_a_mo, nocc, norb)

                NaE3NbNc = np.dot(Nax0.T, e3_dict[('cc', w)].real)

                Nc_starA2Nc = np.dot(Nc_star.T, A2Nc)
                NcA2Nc_star = np.dot(Nc.T, A2Nc_star)

                val_A2 = -(Nc_starA2Nc + NcA2Nc_star)
                val_E3 = NaE3NbNc

                excited_state_dipole_moments.update({
                    (w, 'x'): (val_E3 + val_A2 + ground_state_dipole_x).real
                })

                # Double residue y-component

                ground_state_dipole_y = scf_prop.get_property(
                    'dipole moment')[1]

                A2Nc = self._a2_contract(kc, op_y, d_a_mo, nocc, norb)
                A2Nc_star = self._a2_contract(kc_, op_y, d_a_mo, nocc, norb)

                NaE3NbNc = np.dot(Nay0.T, e3_dict[('cc', w)].real)

                Nc_starA2Nc = np.dot(Nc_star.T, A2Nc)
                NcA2Nc_star = np.dot(Nc.T, A2Nc_star)

                val_A2 = -(Nc_starA2Nc + NcA2Nc_star)
                val_E3 = NaE3NbNc

                excited_state_dipole_moments.update({
                    (w, 'y'): (val_E3 + val_A2 + ground_state_dipole_y).real
                })

                # Double residue z-component

                ground_state_dipole_z = scf_prop.get_property(
                    'dipole moment')[2]

                A2Nc = self._a2_contract(kc, op_z, d_a_mo, nocc, norb)
                A2Nc_star = self._a2_contract(kc_, op_z, d_a_mo, nocc, norb)

                NaE3NbNc = np.dot(Naz0.T, e3_dict[('cc', w)].real)

                Nc_starA2Nc = np.dot(Nc_star.T, A2Nc)
                NcA2Nc_star = np.dot(Nc.T, A2Nc_star)

                val_A2 = -(Nc_starA2Nc + NcA2Nc_star)
                val_E3 = NaE3NbNc

                excited_state_dipole_moments.update({
                    (w, 'z'): (val_E3 + val_A2 + ground_state_dipole_z).real
                })

                # xx

                B2Nc = self._x2_contract(kc, op_x, d_a_mo, nocc, norb)
                A2Nc = self._a2_contract(kc, op_x, d_a_mo, nocc, norb)
                A2Nb = self._a2_contract(kx, op_x, d_a_mo, nocc, norb)

                NaE3NbNc = np.dot(Nax.T, e3_dict[('x', w)].real)
                NaB2Nc = np.dot(Nax.T, B2Nc)
                NbA2Nc = np.dot(N_x.T, A2Nc)
                NcA2Nb = np.dot(Nc.T, A2Nb)

                val_X2 = -(NaB2Nc)
                val_A2 = -(NbA2Nc + NcA2Nb)
                val_E3 = NaE3NbNc

                m[0, 0] = val_E3 + val_A2 + val_X2

                # yy

                B2Nc = self._x2_contract(kc, op_y, d_a_mo, nocc, norb)
                A2Nc = self._a2_contract(kc, op_y, d_a_mo, nocc, norb)
                A2Nb = self._a2_contract(ky, op_y, d_a_mo, nocc, norb)
                NaE3NbNc = np.dot(Nay.T, e3_dict[('y', w)].real)
                NaB2Nc = np.dot(Nay.T, B2Nc)
                NbA2Nc = np.dot(N_y.T, A2Nc)
                NcA2Nb = np.dot(Nc.T, A2Nb)

                val_X2 = -(NaB2Nc)
                val_A2 = -(NbA2Nc + NcA2Nb)
                val_E3 = NaE3NbNc

                m[1, 1] = val_E3 + val_A2 + val_X2

                # zz
                B2Nc = self._x2_contract(kc, op_z, d_a_mo, nocc, norb)
                A2Nc = self._a2_contract(kc, op_z, d_a_mo, nocc, norb)
                A2Nb = self._a2_contract(kz, op_z, d_a_mo, nocc, norb)
                NaE3NbNc = np.dot(Naz.T, e3_dict[('z', w)].real)
                NaB2Nc = np.dot(Naz.T, B2Nc)
                NbA2Nc = np.dot(N_z.T, A2Nc)
                NcA2Nb = np.dot(Nc.T, A2Nb)

                val_X2 = -(NaB2Nc)
                val_A2 = -(NbA2Nc + NcA2Nb)
                val_E3 = NaE3NbNc

                m[2, 2] = val_E3 + val_A2 + val_X2

                # A = x B = y
                B2Nc = self._x2_contract(kc, op_y, d_a_mo, nocc, norb)
                A2Nc = self._a2_contract(kc, op_x, d_a_mo, nocc, norb)
                A2Nb = self._a2_contract(ky, op_x, d_a_mo, nocc, norb)

                NaE3NbNc = np.dot(Nax.T, e3_dict[('y', w)].real)

                NaB2Nc = np.dot(Nax.T, B2Nc)
                NbA2Nc = np.dot(N_y.T, A2Nc)
                NcA2Nb = np.dot(Nc.T, A2Nb)

                val_X2 = -(NaB2Nc)
                val_A2 = -(NbA2Nc + NcA2Nb)
                val_E3 = NaE3NbNc

                m[0, 1] = val_E3 + val_A2 + val_X2
                m[1, 0] = val_E3 + val_A2 + val_X2

                # A = x B = z

                B2Nc = self._x2_contract(kc, op_z, d_a_mo, nocc, norb)
                A2Nc = self._a2_contract(kc, op_x, d_a_mo, nocc, norb)
                A2Nb = self._a2_contract(kz, op_x, d_a_mo, nocc, norb)

                NaE3NbNc = np.dot(Nax.T, e3_dict[('z', w)].real)

                NaB2Nc = np.dot(Nax.T, B2Nc)
                NbA2Nc = np.dot(N_z.T, A2Nc)
                NcA2Nb = np.dot(Nc.T, A2Nb)

                val_X2 = -(NaB2Nc)
                val_A2 = -(NbA2Nc + NcA2Nb)
                val_E3 = NaE3NbNc

                m[0, 2] = val_E3 + val_A2 + val_X2
                m[2, 0] = val_E3 + val_A2 + val_X2

                # yz

                B2Nc = self._x2_contract(kc, op_z, d_a_mo, nocc, norb)
                A2Nc = self._a2_contract(kc, op_y, d_a_mo, nocc, norb)
                A2Nb = self._a2_contract(kz, op_y, d_a_mo, nocc, norb)

                NaE3NbNc = np.dot(Nay.T, e3_dict[('z', w)].real)

                NaB2Nc = np.dot(Nay.T, B2Nc)
                NbA2Nc = np.dot(N_z.T, A2Nc)
                NcA2Nb = np.dot(Nc.T, A2Nb)

                val_X2 = -(NaB2Nc)
                val_A2 = -(NbA2Nc + NcA2Nb)
                val_E3 = NaE3NbNc

                m[1, 2] = val_E3 + val_A2 + val_X2
                m[2, 1] = val_E3 + val_A2 + val_X2

                M_tensors.update({w: m})

        diagonalized_tensors = {}

        if self.rank == mpi_master():
            for key, tensor in M_tensors.items():
                eigvals, eigvecs = np.linalg.eigh(tensor)
                diagonalized_tensors[key] = np.diag(eigvals)

        self.ostream.print_blank()
        w_str = 'Summary of Two-photon Absorption'
        self.ostream.print_header(w_str)
        self.ostream.print_header('=' * (len(w_str) + 2))
        self.ostream.print_blank()

        # conversion factor for TPA cross-sections in GM
        # a0 in cm, c in cm/s, gamma (0.1 eV) in au
        alpha = fine_structure_constant()
        a0 = bohr_in_angstroms() * 1.0e-8
        c = speed_of_light_in_vacuum_in_SI() * 100.0
        gamma = 0.1 / hartree_in_ev()
        au2gm = (8.0 * np.pi**2 * alpha * a0**5) / (c * gamma) * 1.0e+50

        if self.rank == mpi_master():
            tpa_strengths = {'linear': {}, 'circular': {}}
            tpa_cross_sections = {'linear': {}, 'circular': {}}

            for w in M_tensors.keys():
                Df = 0.0
                Dg = 0.0
                for i in range(3):
                    for j in range(3):
                        Mii = M_tensors[w][(i, i)]
                        Mij = M_tensors[w][(i, j)]
                        Mjj = M_tensors[w][(j, j)]
                        Df += (Mii * Mjj).real
                        Dg += (Mij * Mij).real
                Df /= 30.0
                Dg /= 30.0

                D_linear = 2.0 * Df + 4.0 * Dg
                D_circular = -2.0 * Df + 6.0 * Dg

                tpa_strengths['linear'][w] = D_linear
                tpa_strengths['circular'][w] = D_circular

                tpa_cross_sections['linear'][w] = au2gm * w**2 * D_linear
                tpa_cross_sections['circular'][w] = au2gm * w**2 * D_circular

            profiler.check_memory_usage('End of QRF')

            ret_dict = {
                'photon_energies': [-w for w in freqs],
                'transition_moments': M_tensors,
                'cross_sections': tpa_cross_sections,
                'tpa_strengths': tpa_strengths,
                'excited_state_dipole_moments': excited_state_dipole_moments,
                'ground_state_dipole_moments':
                    scf_prop.get_property('dipole moment')
            }

            self._print_results(ret_dict)

            return ret_dict
        else:
            return None

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

        for w_ind, w in enumerate(freqs):

            nx = ComplexResponse.get_full_solution_vector(Nx[('x', w)])
            ny = ComplexResponse.get_full_solution_vector(Nx[('y', w)])
            nz = ComplexResponse.get_full_solution_vector(Nx[('z', w)])

            Nc = LinearResponseEigenSolver.get_full_solution_vector(Xf[w_ind])

            if self.rank == mpi_master():

                kbx = LinearSolver.lrvec2mat(nx, nocc, norb)
                kby = LinearSolver.lrvec2mat(ny, nocc, norb)
                kbz = LinearSolver.lrvec2mat(nz, nocc, norb)

                kc = LinearSolver.lrvec2mat(Nc, nocc, norb)
                kc_ = LinearSolver.lrvec2mat(self.flip_yz(Nc), nocc, norb)

                # create the first order single indexed densiteies #

                Dbx = self.commut_mo_density(kbx, nocc)
                Dby = self.commut_mo_density(kby, nocc)
                Dbz = self.commut_mo_density(kbz, nocc)

                Dc = self.commut_mo_density(kc, nocc)
                Dc_ = self.commut_mo_density(kc_, nocc)

                # create the first order two indexed densities #

                Dbcx = self.commut(kbx, Dc) + self.commut(kc, Dbx)
                Dbcy = self.commut(kby, Dc) + self.commut(kc, Dby)
                Dbcz = self.commut(kbz, Dc) + self.commut(kc, Dbz)
                Dcc_ = self.commut(kc_, Dc) + self.commut(kc, Dc_)

                # Density transformation from MO to AO basis

                Dbx = np.linalg.multi_dot([mo, Dbx, mo.T])
                Dby = np.linalg.multi_dot([mo, Dby, mo.T])
                Dbz = np.linalg.multi_dot([mo, Dbz, mo.T])
                Dc = np.linalg.multi_dot([mo, Dc, mo.T])
                Dc_ = np.linalg.multi_dot([mo, Dc_, mo.T])

                Dbcx = np.linalg.multi_dot([mo, Dbcx, mo.T])
                Dbcy = np.linalg.multi_dot([mo, Dbcy, mo.T])
                Dbcz = np.linalg.multi_dot([mo, Dbcz, mo.T])

                Dcc_ = np.linalg.multi_dot([mo, Dcc_, mo.T])

                distributed_density_1_freq = np.hstack((
                    Dbx.real.reshape(-1, 1),
                    Dby.real.reshape(-1, 1),
                    Dbz.real.reshape(-1, 1),
                    Dc.reshape(-1, 1),
                    Dc_.reshape(-1, 1),
                ))

                distributed_density_2_freq = np.hstack((
                    Dbcx.real.reshape(-1, 1),
                    Dbcy.real.reshape(-1, 1),
                    Dbcz.real.reshape(-1, 1),
                    Dcc_.real.reshape(-1, 1),
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

        :return:
            A dictonary of compounded first-order Fock-matrices
        """

        if self.rank == mpi_master():
            self._print_fock_header()

        # generate key-frequency pairs

        key_freq_pairs = []

        for wb in wi:
            for key in ['F(bc)_x', 'F(bc)_y', 'F(bc)_z', 'F(cc)']:
                key_freq_pairs.append((key, wb))

        # examine checkpoint for distributed Focks

        if self.checkpoint_file is not None:
            fock_file = str(
                Path(self.checkpoint_file).with_suffix('.tpatrans_fock.h5'))
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

            dist_focks = self._comp_nlr_fock(mo, molecule, ao_basis, 'real',
                                             dft_dict, first_order_dens,
                                             second_order_dens, None,
                                             'tpa_quad', profiler)

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

        for w_ind, w in enumerate(wi):
            vec_pack = np.array([
                fo['F(bc)_x'][w].data,
                fo['F(bc)_y'][w].data,
                fo['F(bc)_z'][w].data,
                fo['F(cc)'][w].data,
                fo2[('x', w)].data,
                fo2[('y', w)].data,
                fo2[('z', w)].data,
                fo2[w_ind].data,
            ]).T.copy()

            vec_pack = self._collect_vectors_in_columns(vec_pack)

            Nbx = ComplexResponse.get_full_solution_vector(Nx[('x', w)])
            Nby = ComplexResponse.get_full_solution_vector(Nx[('y', w)])
            Nbz = ComplexResponse.get_full_solution_vector(Nx[('z', w)])

            Nc = LinearResponseEigenSolver.get_full_solution_vector(Xf[w_ind])

            if self.rank != mpi_master():
                continue

            vec_pack = vec_pack.T.copy().reshape(-1, norb, norb)

            (fbc_x, fbc_y, fbc_z, fcc, fbx, fby, fbz, fc) = vec_pack

            fbx = np.conjugate(fbx).T
            fby = np.conjugate(fby).T
            fbz = np.conjugate(fbz).T

            fc = np.conjugate(fc).T * -1 / np.sqrt(2)
            fc_star = np.conjugate(fc).T

            F0_a = fo['F0']

            # Response

            kbx = (LinearSolver.lrvec2mat(Nbx, nocc, norb)).T
            kby = (LinearSolver.lrvec2mat(Nby, nocc, norb)).T
            kbz = (LinearSolver.lrvec2mat(Nbz, nocc, norb)).T
            kc = (LinearSolver.lrvec2mat(Nc, nocc, norb)).T
            kc_star = (LinearSolver.lrvec2mat(self.flip_yz(Nc), nocc, norb)).T

            # x
            xi = self._xi(kbx, kc, fbx, fc, F0_a)
            e3fock = xi.T + 0.5 * fbc_x.T
            e3vec[('x', w)] = self.anti_sym(
                -2 * LinearSolver.lrmat2vec(e3fock, nocc, norb))

            # y
            xi = self._xi(kby, kc, fby, fc, F0_a)
            e3fock = xi.T + 0.5 * fbc_y.T
            e3vec[('y', w)] = self.anti_sym(
                -2 * LinearSolver.lrmat2vec(e3fock, nocc, norb))

            # z
            xi = self._xi(kbz, kc, fbz, fc, F0_a)
            e3fock = xi.T + 0.5 * fbc_z.T
            e3vec[('z', w)] = self.anti_sym(
                -2 * LinearSolver.lrmat2vec(e3fock, nocc, norb))

            # cc
            xi = self._xi(kc_star, kc, fc_star, fc, F0_a)
            e3fock = xi.T + 0.5 * fcc.T
            e3vec[('cc', w)] = self.anti_sym(
                -2 * LinearSolver.lrmat2vec(e3fock, nocc, norb))

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

    @staticmethod
    def get_spectrum(rsp_results, x_data, x_unit, b_value, b_unit):
        """
        Gets two-photon absorption spectrum.

        :param rsp_results:
            A dictonary containing the result of TPA transition calculation.
        :param x_data:
            The list or array of x values.
        :param x_unit:
            The unit of x values.
        :param b_value:
            The value of the broadening parameter.
        :param b_unit:
            The unit of the broadening parameter.

        :return:
            A dictionary containing photon energies and TPA cross-sections.
        """

        assert_msg_critical(
            x_unit.lower() in ['au', 'ev', 'nm'],
            'TpaTransitionDriver.get_spectrum: x_data should be au, ev or nm')

        assert_msg_critical(
            b_unit.lower() in ['au', 'ev'],
            'TpaTransitionDriver.get_spectrum: broadening parameter should be au or ev'
        )

        au2ev = hartree_in_ev()
        auxnm = 1.0 / hartree_in_inverse_nm()

        # conversion factor for TPA cross-sections in GM
        # * a0 in cm
        # * c in cm/s
        # * broadening parameter not included in au2gm
        alpha = fine_structure_constant()
        a0_in_cm = bohr_in_angstroms() * 1.0e-8
        c_in_cm_per_s = speed_of_light_in_vacuum_in_SI() * 100.0
        au2gm = (8.0 * np.pi**2 * alpha * a0_in_cm**5) / c_in_cm_per_s * 1.0e+50

        tpa_ene_au = []
        tpa_str = []

        for w, s in rsp_results['tpa_strengths']['linear'].items():
            tpa_ene_au.append(-w)
            tpa_str.append(s)

        spectrum = {}

        if x_unit.lower() == 'au':
            spectrum['x_label'] = 'Photon energy [a.u.]'
        elif x_unit.lower() == 'ev':
            spectrum['x_label'] = 'Photon energy [eV]'
        elif x_unit.lower() == 'nm':
            spectrum['x_label'] = 'Wavelength [nm]'

        spectrum['y_label'] = 'TPA cross-section [GM]'

        if x_unit.lower() == 'au':
            x_data_au = list(x_data)
        elif x_unit.lower() == 'ev':
            x_data_au = [x / au2ev for x in x_data]
        elif x_unit.lower() == 'nm':
            x_data_au = [auxnm / x for x in x_data]

        if b_unit.lower() == 'au':
            b_au = b_value
        elif b_unit.lower() == 'ev':
            b_au = b_value / au2ev

        y_data = []

        for x_au in x_data_au:
            y = 0.0
            for e, s in zip(tpa_ene_au, tpa_str):
                b_factor = b_au / ((e * 2 - x_au * 2)**2 + b_au**2)
                y += au2gm * b_factor * (s * x_au**2)
            y_data.append(y)

        spectrum['x_data'] = list(x_data)
        spectrum['y_data'] = y_data

        return spectrum
