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
from .firstorderprop import FirstOrderProperties
from .sanitychecks import scf_results_sanity_check, dft_sanity_check
from .errorhandler import assert_msg_critical
from .checkpoint import (check_distributed_focks, read_distributed_focks,
                         write_distributed_focks)


class ShgDriver(NonlinearSolver):
    """
    Implements a quadratic response driver for SHG calculations

    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.

    Instance variables
        - is_converged: The flag for convergence.
        - frequencies: The frequencies.
        - comp: The list of all the gamma tensor components
        - damping: The damping parameter.
        - conv_thresh: The convergence threshold for the solver.
        - max_iter: The maximum number of solver iterations.
        - a_operator: The A operator.
        - b_operator: The B operator.
        - c_operator: The C operator.
        - a_components: Cartesian components of the A operator.
        - b_components: Cartesian components of the B operator.
    """

    def __init__(self, comm=None, ostream=None):
        """
        Initializes the isotropic quadratic response driver for second harmonic
        generation (SHG).
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
        self.frequencies = (0,)
        self.comp = None
        self.damping = 1000.0 / hartree_in_wavenumbers()

        self.a_operator = 'dipole'
        self.b_operator = 'dipole'
        self.c_operator = 'dipole'
        self.a_components = 'xyz'
        self.b_components = 'xyz'

        self.shg_type = 'full'

        # input keywords
        self._input_keywords['response'].update({
            'frequencies': ('seq_range', 'frequencies'),
            'damping': ('float', 'damping parameter'),
            'a_operator': ('str_lower', 'A operator'),
            'b_operator': ('str_lower', 'B operator'),
            'c_operator': ('str_lower', 'C operator'),
            'shg_type': ('str_lower', 'full or reduced SHG calculation'),
        })

    def update_settings(self, rsp_dict, method_dict=None):
        """
        Updates response and method settings in SHG driver

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
        Computes the isotropic quadratic response function for second-harmonic
        generation.

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

        # check SCF results
        scf_results_sanity_check(self, scf_tensors)

        # check dft setup
        dft_sanity_check(self, 'compute', 'nonlinear')

        profiler = Profiler({
            'timing': self.timing,
            'profiling': self.profiling,
            'memory_profiling': self.memory_profiling,
            'memory_tracing': self.memory_tracing,
        })

        if self.rank == mpi_master():
            if self.shg_type == 'reduced':
                title = 'SHG Driver (Reduced) Setup'
            elif self.shg_type == 'full':
                title = 'SHG Driver Setup'
            self._print_header(title)

        start_time = time.time()

        # sanity check
        nalpha = molecule.number_of_alpha_electrons()
        nbeta = molecule.number_of_beta_electrons()
        assert_msg_critical(
            nalpha == nbeta,
            'SHG Driver: not implemented for unrestricted case')

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

        linear_solver = LinearSolver(self.comm, self.ostream)
        a_grad = linear_solver.get_complex_prop_grad(self.a_operator,
                                                     self.a_components,
                                                     molecule, ao_basis,
                                                     scf_tensors)
        b_grad = linear_solver.get_complex_prop_grad(self.b_operator,
                                                     self.b_components,
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

        # Storing the dipole integral matrices used for the X[3],X[2],A[3] and
        # A[2] contractions in MO basis
        wa = [sum(x) for x in zip(self.frequencies, self.frequencies)]

        freqpairs = [wl for wl in zip(self.frequencies, self.frequencies)]

        AB = {}

        if self.rank == mpi_master():
            A = {(op, w): v for op, v in zip('xyz', a_grad) for w in wa}
            B = {
                (op, w): v for op, v in zip('xyz', b_grad)
                for w in self.frequencies
            }

            AB.update(A)
            AB.update(B)

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
            'frequencies', 'damping', 'norm_thresh', 'lindep_thresh',
            'conv_thresh', 'max_iter', 'eri_thresh', 'qq_type', 'timing',
            'memory_profiling', 'batch_size', 'restart', 'xcfun', 'grid_level',
            'potfile', 'electric_field', 'program_end_time'
        ]

        for key in cpp_keywords:
            setattr(N_drv, key, getattr(self, key))

        if self.checkpoint_file is not None:
            N_drv.checkpoint_file = str(
                Path(self.checkpoint_file).with_suffix('.shg.h5'))

        N_results = N_drv.compute(molecule, ao_basis, scf_tensors, AB)

        self._is_converged = N_drv.is_converged

        Nx = N_results['solutions']
        Focks = N_results['focks']

        profiler.check_memory_usage('CPP')

        # Compute the isotropic parallel beta vector

        beta = self.compute_quad_components(Focks, freqpairs, X, d_a_mo, Nx,
                                            self.comp, scf_tensors, molecule,
                                            ao_basis, profiler)

        profiler.end(self.ostream)

        # Compute dipole vector
        scf_prop = FirstOrderProperties(self.comm, self.ostream)
        scf_prop.compute_scf_prop(molecule, ao_basis, scf_tensors)

        if self.rank == mpi_master():

            dip = scf_prop.get_property('dipole moment')
            # Norm of dipole
            dip_norm = np.linalg.norm(dip)

            # Compute isotropic beta along molecular dipole

            self.ostream.print_blank()
            w_str = 'Electronic dipole moment'
            self.ostream.print_header(w_str)
            self.ostream.print_header('=' * (len(w_str) + 2))

            self.ostream.print_blank()
            title = '{:<12s}{:>12s}{:>19s}{:>20s}'.format(
                'Component', 'Frequency', 'Real', 'Imaginary')
            width = len(title)
            self.ostream.print_header(title.ljust(width))
            self.ostream.print_header(('-' * len(title)).ljust(width))
            self._print_component('mu_x', 0, dip[0], width)
            self._print_component('mu_y', 0, dip[1], width)
            self._print_component('mu_z', 0, dip[2], width)
            self._print_component('|mu|', 0, dip_norm, width)
            self.ostream.print_blank()
            self.ostream.print_blank()

            w_str = 'SHG hyperpolarizability'
            self.ostream.print_header(w_str)
            self.ostream.print_header('=' * (len(w_str) + 2))
            self.ostream.print_blank()
            w_str = 'Observable quantity corresponding to EFISHG experiments with '
            self.ostream.print_header(w_str.ljust(width))
            w_str = 'parallel external fields:'
            self.ostream.print_header(w_str.ljust(width))
            self.ostream.print_blank()

            w_str = 'beta(vec)_i = 1/5*sum_{j=x,y,z}(beta_ijj + beta_jij + beta_jji)'
            self.ostream.print_header(w_str.ljust(width))
            self.ostream.print_blank()
            w_str = 'beta = sum_{i=x,y,z}(beta(vec)_i * mu_i / |mu|)'
            self.ostream.print_header(w_str.ljust(width))
            self.ostream.print_blank()

            beta_bar = {}
            for freq in beta:
                beta_bar[freq] = 0.0
                for a in range(len(beta[freq])):
                    beta_bar[freq] += 1.0 / dip_norm * dip[a] * beta[freq][a]

                self.ostream.print_header(title.ljust(width))
                self.ostream.print_header(('-' * len(title)).ljust(width))
                self._print_component('beta(vec)_x', freq, beta[freq][0], width)
                self._print_component('beta(vec)_y', freq, beta[freq][1], width)
                self._print_component('beta(vec)_z', freq, beta[freq][2], width)
                self._print_component('beta', freq, beta_bar[freq], width)
                self.ostream.print_blank()

            title = 'Reference: '
            title += 'K. Ahmadzadeh, X. Li, Z. Rinkevicius, P. Norman,'
            self.ostream.print_header(title.ljust(width))
            title = '           Electron. Struct. 2022, 4, 044004.'
            self.ostream.print_header(title.ljust(width))
            title = '           (DOI: 10.1088/2516-1075/aca859)'
            self.ostream.print_header(title.ljust(width))
            self.ostream.print_blank()

            self.ostream.print_blank()
            valstr = '*** Time spent in SHG calculation: '
            valstr += '{:.2f} sec ***'.format(time.time() - start_time)
            self.ostream.print_header(valstr)
            self.ostream.print_blank()
            self.ostream.flush()

            return {
                'dipole': dip,
                'beta': beta,
                'beta_bar': beta_bar,
            }

        else:
            return None

    def compute_quad_components(self, Focks, freqpairs, X, d_a_mo, Nx, track,
                                scf_tensors, molecule, ao_basis, profiler):
        """
        Computes all the relevent terms to compute the isotropic quadratic
        response function used for SHG.

        :param w:
            A list of all the frequencies
        :param X:
            A dictonary of matricies containing all the dipole integrals
        :param d_a_mo:
            The SCF density in MO basis
        :param Nx:
            A dictonary containing all the response vectors in distributed
            form.
        :param scf_tensors:
            The dictionary of tensors from converged SCF wavefunction.
        :param molecule:
            The molecule.
        :param basis:
            The AO basis.
        :param profiler:
            The profiler.

        :return:
            A dictionary containing all the relevent terms for SHG
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
                                       ao_basis, dft_dict, profiler)

        profiler.check_memory_usage('Focks')

        e3_dict = self.get_e3(freqpairs, Nx, fock_dict, Focks, nocc, norb)

        profiler.check_memory_usage('E[3]')

        beta = {}

        for (wb, wc) in freqpairs:

            Na = {
                'x': ComplexResponse.get_full_solution_vector(Nx[('x',
                                                                  (wb + wc))]),
                'y': ComplexResponse.get_full_solution_vector(Nx[('y',
                                                                  (wb + wc))]),
                'z': ComplexResponse.get_full_solution_vector(Nx[('z',
                                                                  (wb + wc))]),
            }

            Nb = {
                'x': ComplexResponse.get_full_solution_vector(Nx[('x', wb)]),
                'y': ComplexResponse.get_full_solution_vector(Nx[('y', wb)]),
                'z': ComplexResponse.get_full_solution_vector(Nx[('z', wb)]),
            }

            if self.rank == mpi_master():

                kX = {
                    ('x', wb): self.complex_lrvec2mat(Nb['x'], nocc, norb),
                    ('y', wb): self.complex_lrvec2mat(Nb['y'], nocc, norb),
                    ('z', wb): self.complex_lrvec2mat(Nb['z'], nocc, norb),
                }

                NaE3NbNc = {'x': 0.0, 'y': 0.0, 'z': 0.0}

                NaE3NbNc['x'] += np.dot(Na['x'].T, e3_dict[('sig_x', wb)])
                NaE3NbNc['x'] += 2 * np.dot(Na['y'].T, e3_dict[('lam_xy', wb)])
                NaE3NbNc['x'] += 2 * np.dot(Na['z'].T, e3_dict[('lam_xz', wb)])

                NaE3NbNc['y'] += np.dot(Na['y'].T, e3_dict[('sig_y', wb)])
                NaE3NbNc['y'] += 2 * np.dot(Na['z'].T, e3_dict[('lam_yz', wb)])
                NaE3NbNc['y'] += 2 * np.dot(Na['x'].T, e3_dict[('lam_xy', wb)])

                NaE3NbNc['z'] += np.dot(Na['z'].T, e3_dict[('sig_z', wb)])
                NaE3NbNc['z'] += 2 * np.dot(Na['x'].T, e3_dict[('lam_xz', wb)])
                NaE3NbNc['z'] += 2 * np.dot(Na['y'].T, e3_dict[('lam_yz', wb)])

                A2 = {'x': 0.0, 'y': 0.0, 'z': 0.0}
                X2 = {'x': 0.0, 'y': 0.0, 'z': 0.0}

                # pre-compute terms for A2 and X2
                a2_kX_X = {}
                x2_kX_X = {}
                for cart_1 in 'xyz':
                    for cart_2 in 'xyz':
                        a2_kX_X[cart_1 + cart_2] = self._a2_contract(
                            kX[(cart_1, wb)], X[cart_2], d_a_mo, nocc, norb)
                        x2_kX_X[cart_1 + cart_2] = self._x2_contract(
                            kX[(cart_1, wb)], X[cart_2], d_a_mo, nocc, norb)

                for eta in 'xyz':
                    for cart in 'xyz':
                        # A2 contractions
                        A2[cart] -= 2 * np.dot(Nb[eta].T, a2_kX_X[eta + cart])
                        A2[cart] -= 2 * np.dot(Nb[cart].T, a2_kX_X[eta + eta])
                        A2[cart] -= 2 * np.dot(Nb[eta].T, a2_kX_X[cart + eta])
                        # X2 contractions
                        X2[cart] -= 2 * np.dot(Na[eta].T, x2_kX_X[eta + cart])
                        X2[cart] -= 2 * np.dot(Na[cart].T, x2_kX_X[eta + eta])
                        X2[cart] -= 2 * np.dot(Na[eta].T, x2_kX_X[cart + eta])

                # β_i(ω) = -1/5 * (<< i;j,j >> + << j;i,j >> + << j;j,i >>),
                # for j in {x,y,z}

                beta[wb] = (
                    -1 / 5 * (NaE3NbNc['x'] + A2['x'] + X2['x']),
                    -1 / 5 * (NaE3NbNc['y'] + A2['y'] + X2['y']),
                    -1 / 5 * (NaE3NbNc['z'] + A2['z'] + X2['z']),
                )

        profiler.check_memory_usage('End of SHG')

        return beta

    def get_densities(self, freqpairs, Nx, mo, nocc, norb):
        """
        Computes the densities needed for the perturbed Fock matrices.

        :param freqpairs:
            A list of the frequency pairs
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
            first_order_dens:
             A list of first-order one-time tranformed compounded densities
            second_order_dens:
             A list of first-order two-time tranformed compounded densities
        """

        distributed_density_1 = None
        distributed_density_2 = None

        for (wb, wc) in freqpairs:

            nx = ComplexResponse.get_full_solution_vector(Nx[('x', wb)])
            ny = ComplexResponse.get_full_solution_vector(Nx[('y', wb)])
            nz = ComplexResponse.get_full_solution_vector(Nx[('z', wb)])

            if self.rank == mpi_master():

                k_x = self.complex_lrvec2mat(nx, nocc, norb)
                k_y = self.complex_lrvec2mat(ny, nocc, norb)
                k_z = self.complex_lrvec2mat(nz, nocc, norb)

                # create the first order single indexed densiteies #

                D_x = self.commut_mo_density(k_x, nocc)
                D_y = self.commut_mo_density(k_y, nocc)
                D_z = self.commut_mo_density(k_z, nocc)

                # create the first order two indexed densities #

                D_sig_x = 4 * self.commut(k_x, D_x)
                D_sig_x += 2 * (self.commut(k_x, D_x) + self.commut(k_y, D_y) +
                                self.commut(k_z, D_z))

                D_sig_y = 4 * self.commut(k_y, D_y)
                D_sig_y += 2 * (self.commut(k_x, D_x) + self.commut(k_y, D_y) +
                                self.commut(k_z, D_z))

                D_sig_z = 4 * self.commut(k_z, D_z)
                D_sig_z += 2 * (self.commut(k_x, D_x) + self.commut(k_y, D_y) +
                                self.commut(k_z, D_z))

                D_lam_xy = self.commut(k_x, D_y) + self.commut(k_y, D_x)
                D_lam_xz = self.commut(k_x, D_z) + self.commut(k_z, D_x)
                D_lam_yz = self.commut(k_y, D_z) + self.commut(k_z, D_y)

                # density transformation from MO to AO basis

                D_x = np.linalg.multi_dot([mo, D_x, mo.T])
                D_y = np.linalg.multi_dot([mo, D_y, mo.T])
                D_z = np.linalg.multi_dot([mo, D_z, mo.T])

                D_sig_x = np.linalg.multi_dot([mo, D_sig_x, mo.T])
                D_sig_y = np.linalg.multi_dot([mo, D_sig_y, mo.T])
                D_sig_z = np.linalg.multi_dot([mo, D_sig_z, mo.T])

                D_lam_xy = np.linalg.multi_dot([mo, D_lam_xy, mo.T])
                D_lam_xz = np.linalg.multi_dot([mo, D_lam_xz, mo.T])
                D_lam_yz = np.linalg.multi_dot([mo, D_lam_yz, mo.T])

                dist_den_1_freq = np.hstack((
                    D_x.real.reshape(-1, 1),
                    D_x.imag.reshape(-1, 1),
                    D_y.real.reshape(-1, 1),
                    D_y.imag.reshape(-1, 1),
                    D_z.real.reshape(-1, 1),
                    D_z.imag.reshape(-1, 1),
                ))

                if self.shg_type == 'reduced':
                    dist_den_2_freq = np.hstack((
                        D_sig_x.real.reshape(-1, 1),
                        D_sig_y.real.reshape(-1, 1),
                        D_sig_z.real.reshape(-1, 1),
                        D_lam_xy.real.reshape(-1, 1),
                        D_lam_xz.real.reshape(-1, 1),
                        D_lam_yz.real.reshape(-1, 1),
                    ))

                elif self.shg_type == 'full':
                    dist_den_2_freq = np.hstack((
                        D_sig_x.real.reshape(-1, 1),
                        D_sig_x.imag.reshape(-1, 1),
                        D_sig_y.real.reshape(-1, 1),
                        D_sig_y.imag.reshape(-1, 1),
                        D_sig_z.real.reshape(-1, 1),
                        D_sig_z.imag.reshape(-1, 1),
                        D_lam_xy.real.reshape(-1, 1),
                        D_lam_xy.imag.reshape(-1, 1),
                        D_lam_xz.real.reshape(-1, 1),
                        D_lam_xz.imag.reshape(-1, 1),
                        D_lam_yz.real.reshape(-1, 1),
                        D_lam_yz.imag.reshape(-1, 1),
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
                      dft_dict=None,
                      profiler=None):
        """
        Computes the compounded Fock matrices used for the
        isotropic quadratic response function used for SHG

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
            A dictonary of compounded first-order two-time Fock-matrices. For
            SHG there are 6 real and 6 imaginary Fock matrices.
        """

        if self.rank == mpi_master():
            self._print_fock_header()

        # generate key-frequency pairs

        key_freq_pairs = []

        for (wb, wc) in wi:
            for key in [
                    'F_sig_x', 'F_sig_y', 'F_sig_z', 'F_lam_xy', 'F_lam_xz',
                    'F_lam_yz'
            ]:
                key_freq_pairs.append((key, wb))

        # examine checkpoint for distributed Focks

        if self.checkpoint_file is not None:
            if self.shg_type == 'reduced':
                fock_suffix = '.shg_fock_red.h5'
            elif self.shg_type == 'full':
                fock_suffix = '.shg_fock_full.h5'
            fock_file = str(Path(self.checkpoint_file).with_suffix(fock_suffix))
        else:
            fock_file = None

        if self.restart:
            if self.rank == mpi_master():
                self.restart = check_distributed_focks(fock_file,
                                                       key_freq_pairs)
            self.restart = self.comm.bcast(self.restart, mpi_master())

        # examine checkpoint for distributed Focks

        if self.restart:
            dist_focks = read_distributed_focks(fock_file, self.comm,
                                                self.ostream)
        else:
            time_start_fock = time.time()

            if self.shg_type == 'reduced':
                dist_focks = self._comp_nlr_fock(mo, molecule, ao_basis, 'real',
                                                 dft_dict, first_order_dens,
                                                 second_order_dens, None,
                                                 'shg_red', profiler)
            elif self.shg_type == 'full':
                dist_focks = self._comp_nlr_fock(mo, molecule, ao_basis,
                                                 'real_and_imag', dft_dict,
                                                 first_order_dens,
                                                 second_order_dens, None, 'shg',
                                                 profiler)

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

    def get_e3(self, wi, Nx, fo, fo2, nocc, norb):
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
            response function for SHG
        """

        e3vec = {}

        for (wb, wc) in wi:

            vec_pack = np.array([
                fo2[('x', wb)].data,
                fo2[('y', wb)].data,
                fo2[('z', wb)].data,
                fo['F_sig_x'][wb].data,
                fo['F_sig_y'][wb].data,
                fo['F_sig_z'][wb].data,
                fo['F_lam_xy'][wb].data,
                fo['F_lam_xz'][wb].data,
                fo['F_lam_yz'][wb].data,
            ]).T.copy()

            vec_pack = self._collect_vectors_in_columns(vec_pack)

            nx = ComplexResponse.get_full_solution_vector(Nx[('x', wb)])
            ny = ComplexResponse.get_full_solution_vector(Nx[('y', wb)])
            nz = ComplexResponse.get_full_solution_vector(Nx[('z', wb)])

            if self.rank != mpi_master():
                continue

            vec_pack = vec_pack.T.copy().reshape(-1, norb, norb)

            (F_x, F_y, F_z, F_sig_x, F_sig_y, F_sig_z, F_lam_xy, F_lam_xz,
             F_lam_yz) = vec_pack

            F_x = np.conjugate(F_x).T
            F_y = np.conjugate(F_y).T
            F_z = np.conjugate(F_z).T

            F0_a = fo['F0']

            # Response

            k_x = (self.complex_lrvec2mat(nx, nocc, norb)).T
            k_y = (self.complex_lrvec2mat(ny, nocc, norb)).T
            k_z = (self.complex_lrvec2mat(nz, nocc, norb)).T

            # Make all Xi terms

            xi_sig_x = 3 * self._xi(k_x, k_x, F_x, F_x, F0_a)
            xi_sig_x += self._xi(k_y, k_y, F_y, F_y, F0_a)
            xi_sig_x += self._xi(k_z, k_z, F_z, F_z, F0_a)

            xi_sig_y = 3 * self._xi(k_y, k_y, F_y, F_y, F0_a)
            xi_sig_y += self._xi(k_z, k_z, F_z, F_z, F0_a)
            xi_sig_y += self._xi(k_x, k_x, F_x, F_x, F0_a)

            xi_sig_z = 3 * self._xi(k_z, k_z, F_z, F_z, F0_a)
            xi_sig_z += self._xi(k_x, k_x, F_x, F_x, F0_a)
            xi_sig_z += self._xi(k_y, k_y, F_y, F_y, F0_a)

            xi_lam_xy = self._xi(k_x, k_y, F_x, F_y, F0_a)
            xi_lam_xz = self._xi(k_x, k_z, F_x, F_z, F0_a)
            xi_lam_yz = self._xi(k_y, k_z, F_y, F_z, F0_a)

            # Store Transformed Fock Matrices

            e3fock_sig_x = xi_sig_x.T + (0.5 * F_sig_x).T
            e3fock_sig_y = xi_sig_y.T + (0.5 * F_sig_y).T
            e3fock_sig_z = xi_sig_z.T + (0.5 * F_sig_z).T

            e3fock_lam_xy = xi_lam_xy.T + (0.5 * F_lam_xy).T
            e3fock_lam_xz = xi_lam_xz.T + (0.5 * F_lam_xz).T
            e3fock_lam_yz = xi_lam_yz.T + (0.5 * F_lam_yz).T

            # Anti sym the Fock matrices and convert them to vectors

            e3vec[('sig_x', wb)] = self.anti_sym(
                -2 * LinearSolver.lrmat2vec(e3fock_sig_x, nocc, norb))
            e3vec[('sig_y', wb)] = self.anti_sym(
                -2 * LinearSolver.lrmat2vec(e3fock_sig_y, nocc, norb))
            e3vec[('sig_z', wb)] = self.anti_sym(
                -2 * LinearSolver.lrmat2vec(e3fock_sig_z, nocc, norb))

            e3vec[('lam_xy', wb)] = self.anti_sym(
                -2 * LinearSolver.lrmat2vec(e3fock_lam_xy, nocc, norb))
            e3vec[('lam_xz', wb)] = self.anti_sym(
                -2 * LinearSolver.lrmat2vec(e3fock_lam_xz, nocc, norb))
            e3vec[('lam_yz', wb)] = self.anti_sym(
                -2 * LinearSolver.lrmat2vec(e3fock_lam_yz, nocc, norb))

        return e3vec

    def _print_component(self, label, freq, value, width):
        """
        Prints SHG component.

        :param label:
            The label
        :param freq:
            The frequency
        :param value:
            The complex value
        :param width:
            The width for the output
        """

        w_str = '{:<12s}{:12.4f}{:19.8f}{:19.8f}j'.format(
            label, freq, value.real, value.imag)
        self.ostream.print_header(w_str.ljust(width))
