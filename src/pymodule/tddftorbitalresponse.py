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

import numpy as np
import time

from .veloxchemlib import mpi_master
from .veloxchemlib import XCIntegrator
from .profiler import Profiler
from .cphfsolver import CphfSolver
from .distributedarray import DistributedArray
from .dftutils import get_default_grid_level
from .inputparser import parse_input


class TddftOrbitalResponse(CphfSolver):
    """
    Implements orbital response Lagrange multipliers computation using a
    conjugate gradient scheme for the random phase approximation (RPA)
    level of theory.
    """

    def __init__(self, comm=None, ostream=None):
        """
        Initializes orbital response computation driver to default setup.

        :param comm:
            The MPI communicator.
        :param ostream:
            The output stream.
        """

        super().__init__(comm, ostream)

        self.orbrsp_type = 'tddftgrad'

        self.tamm_dancoff = False
        self.state_deriv_index = None

        self._input_keywords['orbitalresponse'].update({
            'tamm_dancoff': ('bool', 'whether RPA or TDA is calculated'),
            'state_deriv_index': ('seq_fixed_int', 'excited state of interest'),
        })

    def update_settings(self, orbrsp_dict, method_dict=None):
        """
        Updates response and method settings in orbital response computation
        driver.

        :param orbrsp_dict:
            The dictionary of orbital response settings.
        :param method_dict:
            The dictionary of method settings.
        """

        super().update_settings(orbrsp_dict, method_dict)

        orbrsp_keywords = {
            key: val[0]
            for key, val in self._input_keywords['orbitalresponse'].items()
        }

        parse_input(self, orbrsp_keywords, orbrsp_dict)

    def compute(self, molecule, basis, scf_tensors, rsp_results):
        """ Computes the lambda orbital response multipliers and
            excited state first order properties.

        :param molecule:
            The molecule.
        :param basis:
            The AO basis set.
        :param scf_tensors:
            The dictionary of tensors from converged SCF wavefunction.
        :param rsp_results:
            The dictionary of results from a converged linear response
            calculation.
        """

        super().compute(molecule, basis, scf_tensors, rsp_results)

    @staticmethod
    def get_full_solution_vector(solution):
        """
        Gets a full solution vector from the distributed solution.

        :param solution:
            The distributed solution as a tuple.

        :return:
            The full solution vector.
        """

        x_ger = solution.get_full_vector(0)
        x_ung = solution.get_full_vector(1)

        if solution.rank == mpi_master():
            x_ger_full = np.hstack((x_ger, x_ger))
            x_ung_full = np.hstack((x_ung, -x_ung))
            return x_ger_full + x_ung_full
        else:
            return None

    def compute_rhs(self, molecule, basis, scf_tensors, eri_dict, dft_dict,
                    pe_dict, rsp_results):
        """
        Computes the right-hand side (RHS) of the RPA orbital response equation
        including the necessary density matrices using molecular data.

        :param molecule:
            The molecule.
        :param basis:
            The AO basis set.
        :param scf_tensors:
            The dictionary of tensors from the converged SCF calculation.
        :param rsp_results:
            The results from a converged linear response calculation.

        :return:
            A dictionary containing the orbital-response RHS and
            unrelaxed one-particle density.
        """

        # TODO: replace with a sanity check?
        if 'eigenvectors' in rsp_results:
            self.tamm_dancoff = True

        profiler = Profiler({
            'timing': self.timing,
            'profiling': self.profiling,
            'memory_profiling': self.memory_profiling,
            'memory_tracing': self.memory_tracing,
        })

        profiler.set_timing_key('RHS')

        rhs_t0 = time.time()

        self.ostream.print_info(
            "Computing the right-hand side (RHS) of CPHF/CPKS equations...")
        self.ostream.print_blank()
        self.ostream.flush()

        profiler.start_timer('Prep')

        # Workflow:
        # 1) Construct the necessary density matrices
        # 2) Construct the RHS
        # 3) Construct the initial guess => in parent class
        # 4) Write the linear operator for matrix-vector product
        #    => in parent class
        # 5) Run the conjugate gradient => in parent class

        if self.rank == mpi_master():

            # 1) Calculate unrelaxed one-particle and transition density matrix
            ovlp = scf_tensors['S']
            mo = scf_tensors['C_alpha']

            nocc = molecule.number_of_alpha_electrons()
            mo_occ = mo[:, :nocc].copy()
            mo_vir = mo[:, nocc:].copy()
            nvir = mo_vir.shape[1]
            nao = mo.shape[0]

            # Take vectors of interest and convert to matrix form
            # n_states x nocc x nvir
            if self.state_deriv_index is None:
                # if no states are selected, calculate all
                # first excitd state is S1 (index starts at 1)
                self.state_deriv_index = list(
                    np.arange(1,
                              len(rsp_results['eigenvalues']) + 1))

            # number of degrees of freedon:
            dof = len(self.state_deriv_index)
            exc_vec = np.zeros((dof, nocc, nvir))
            deexc_vec = np.zeros((dof, nocc, nvir))
        else:
            self.state_deriv_index = None
            dof = None

        dof = self.comm.bcast(dof, root=mpi_master())
        self.state_deriv_index = self.comm.bcast(self.state_deriv_index,
                                                 root=mpi_master())

        for i in range(dof):
            ivec = self.state_deriv_index[i] - 1

            if self.tamm_dancoff:
                if self.rank == mpi_master():
                    exc_vec[i] = (rsp_results['eigenvectors'][:nocc * nvir,
                                                              ivec].reshape(
                                                                  nocc, nvir))
            else:
                eigenvector = self.get_full_solution_vector(
                    rsp_results['eigenvectors_distributed'][ivec])
                if self.rank == mpi_master():
                    exc_vec[i] = eigenvector[:nocc * nvir].reshape(nocc, nvir)
                    deexc_vec[i] = eigenvector[nocc * nvir:].reshape(nocc, nvir)

        if self.rank == mpi_master():
            # Construct plus/minus combinations of excitation
            # and de-excitation part
            x_plus_y = exc_vec + deexc_vec
            x_minus_y = exc_vec - deexc_vec

            # Transform the vectors to the AO basis
            x_plus_y_ao = np.zeros((dof, nao, nao))
            x_minus_y_ao = np.zeros((dof, nao, nao))
            dm_oo = np.zeros((dof, nocc, nocc))
            dm_vv = np.zeros((dof, nvir, nvir))
            unrel_dm_ao = np.zeros((dof, nao, nao))
            for i in range(dof):
                x_plus_y_ao[i] = np.linalg.multi_dot(
                    [mo_occ, x_plus_y[i], mo_vir.T])
                x_minus_y_ao[i] = np.linalg.multi_dot(
                    [mo_occ, x_minus_y[i], mo_vir.T])
                dm_oo[i] = -0.5 * (
                    np.linalg.multi_dot([x_plus_y[i], x_plus_y[i].T]) +
                    np.linalg.multi_dot([x_minus_y[i], x_minus_y[i].T]))
                dm_vv[i] = 0.5 * (
                    np.linalg.multi_dot([x_plus_y[i].T, x_plus_y[i]]) +
                    np.linalg.multi_dot([x_minus_y[i].T, x_minus_y[i]]))
                unrel_dm_ao[i] = (
                    np.linalg.multi_dot([mo_occ, dm_oo[i], mo_occ.T]) +
                    np.linalg.multi_dot([mo_vir, dm_vv[i], mo_vir.T]))

            # Make a list of unrel. DMs and excittion vectors:
            dm_ao_list = (list(unrel_dm_ao) + list(x_plus_y_ao) +
                          list(x_minus_y_ao))

            # 2) Construct the right-hand side

            if self._dft:
                # 3) Construct density matrices for E[3] term:
                # XCIntegrator expects a DM with real and imaginary part,
                # so we set the imaginary part to zero.
                perturbed_dm_ao_list = []
                zero_dm_ao_list = []

                # for each vector, we need to create a list with these elements:
                # x_minus_y_ao[i], 0*x_minus_y_ao[i],
                # x_minus_y_ao[i], 0*x_minus_y_ao[i];
                # and a list with 2 list of zeros for each 4 elements above.

                for s in range(dof):
                    perturbed_dm_ao_list.extend([
                        x_minus_y_ao[s], 0 * x_minus_y_ao[s], x_minus_y_ao[s],
                        0 * x_minus_y_ao[s]
                    ])

                    zero_dm_ao_list.extend(
                        [0 * x_minus_y_ao[s], 0 * x_minus_y_ao[s]])

        else:
            dof = None
            dm_ao_list = None

            if self._dft:
                perturbed_dm_ao_list = None
                zero_dm_ao_list = None

        dof = self.comm.bcast(dof, root=mpi_master())

        if self._dft:
            # TODO: bcast array by array
            perturbed_dm_ao_list = self.comm.bcast(perturbed_dm_ao_list,
                                                   root=mpi_master())
            # TODO: bcast array by array
            dm_ao_list = self.comm.bcast(dm_ao_list, root=mpi_master())
            # TODO: bcast array by array
            # TODO: consider only bcast size of zero dm
            zero_dm_ao_list = self.comm.bcast(zero_dm_ao_list,
                                              root=mpi_master())

        profiler.stop_timer('Prep')

        molgrid = dft_dict['molgrid']
        gs_density = dft_dict['gs_density']

        fock_ao_rhs = self._comp_lr_fock(dm_ao_list, molecule, basis, eri_dict,
                                         dft_dict, pe_dict, profiler)

        profiler.start_timer('Kxc')

        if self._dft:
            # Fock matrix for computing gxc
            fock_gxc_ao = []
            for i_mat in range(len(zero_dm_ao_list)):
                fock_gxc_ao.append(zero_dm_ao_list[i_mat].copy())
        else:
            fock_gxc_ao = None

        if self._dft:
            # Quadratic response routine for TDDFT E[3] term g^xc
            xc_drv = XCIntegrator()
            # molgrid.partition_grid_points()
            # molgrid.distribute_counts_and_displacements(self.rank,
            #                                    self.nodes, self.comm)
            xc_drv.integrate_kxc_fock(fock_gxc_ao, molecule, basis,
                                      perturbed_dm_ao_list, zero_dm_ao_list,
                                      gs_density, molgrid,
                                      self.xcfun.get_func_label(), "qrf")

            for idx in range(len(fock_gxc_ao)):
                fock_gxc_ao[idx] = self.comm.reduce(fock_gxc_ao[idx],
                                                    root=mpi_master())

        profiler.stop_timer('Kxc')

        profiler.start_timer('RHS_MO')

        dist_rhs_mo = []

        # Calculate the RHS and transform it to the MO basis
        if self.rank == mpi_master():
            # Extract the 1PDM contributions
            fock_ao_rhs_1pdm = np.zeros((dof, nao, nao))
            for ifock in range(dof):
                fock_ao_rhs_1pdm[ifock] = fock_ao_rhs[ifock]

            # Extract the excitation vector contributions
            fock_ao_rhs_x_plus_y = np.zeros((dof, nao, nao))
            fock_ao_rhs_x_minus_y = np.zeros((dof, nao, nao))
            for ifock in range(dof):
                fock_ao_rhs_x_plus_y[ifock] = fock_ao_rhs[dof + ifock]
                fock_ao_rhs_x_minus_y[ifock] = fock_ao_rhs[2 * dof + ifock]

        # Transform to MO basis:
        for i in range(dof):
            if self.rank == mpi_master():
                fmo_rhs_1pdm = 0.5 * np.linalg.multi_dot(
                    [mo_occ.T, fock_ao_rhs_1pdm[i], mo_vir])

                sdp_pds = 0.25 * (
                    np.linalg.multi_dot(
                        [ovlp, x_plus_y_ao[i], fock_ao_rhs_x_plus_y[i].T]) +
                    np.linalg.multi_dot(
                        [ovlp, x_minus_y_ao[i], fock_ao_rhs_x_minus_y[i].T]) -
                    np.linalg.multi_dot(
                        [fock_ao_rhs_x_plus_y[i].T, x_plus_y_ao[i], ovlp]) -
                    np.linalg.multi_dot(
                        [fock_ao_rhs_x_minus_y[i].T, x_minus_y_ao[i], ovlp]) -
                    np.linalg.multi_dot(
                        [ovlp, x_plus_y_ao[i], fock_ao_rhs_x_plus_y[i]]) +
                    np.linalg.multi_dot(
                        [ovlp, x_minus_y_ao[i], fock_ao_rhs_x_minus_y[i]]) +
                    np.linalg.multi_dot(
                        [fock_ao_rhs_x_plus_y[i], x_plus_y_ao[i], ovlp]) -
                    np.linalg.multi_dot(
                        [fock_ao_rhs_x_minus_y[i], x_minus_y_ao[i], ovlp]))

                rhs_mo_i = (fmo_rhs_1pdm +
                            np.linalg.multi_dot([mo_occ.T, sdp_pds, mo_vir]))

                # Add DFT E[3] contribution to the RHS:
                if self._dft:
                    gxc_mo_i = np.linalg.multi_dot(
                        [mo_occ.T, fock_gxc_ao[2 * i], mo_vir])
                    rhs_mo_i += 0.25 * gxc_mo_i

                rhs_mo_i = rhs_mo_i.reshape(nocc * nvir)
            else:
                rhs_mo_i = None
            dist_rhs_mo.append(DistributedArray(rhs_mo_i, self.comm))

        profiler.stop_timer('RHS_MO')

        profiler.print_timing(self.ostream)
        profiler.print_profiling_summary(self.ostream)

        self.ostream.print_info("RHS of CPHF/CPKS equations computed in " +
                                f"{time.time() - rhs_t0:.2f} sec.")
        self.ostream.print_blank()
        self.ostream.flush()

        if self.rank == mpi_master():
            return {
                'dist_cphf_rhs': dist_rhs_mo,
                'density_occ_occ': dm_oo,
                'density_vir_vir': dm_vv,
                'x_plus_y_ao': x_plus_y_ao,
                'x_minus_y_ao': x_minus_y_ao,
                'unrelaxed_density_ao': unrel_dm_ao,
                'fock_ao_rhs': fock_ao_rhs,
                'fock_gxc_ao': fock_gxc_ao,  # None if not DFT
            }
        else:
            return {
                'dist_cphf_rhs': dist_rhs_mo,
            }

    def compute_omega(self, molecule, basis, scf_tensors):
        """
        Calculates the omega Lagrange multipliers for the overlap matrix.

        :param molecule:
            The molecule.
        :param basis.
            The basis set.
        :param scf_tensors.
            The scf tensors.

        :return:
            a numpy array containing the Lagrange multipliers in AO basis.
        """

        profiler = Profiler({
            'timing': self.timing,
            'profiling': self.profiling,
            'memory_profiling': self.memory_profiling,
            'memory_tracing': self.memory_tracing,
        })

        profiler.set_timing_key('Omega')

        omega_t0 = time.time()

        self.ostream.print_info("Computing the omega Lagrange multipliers...")
        self.ostream.print_blank()
        self.ostream.flush()

        profiler.start_timer('Prep')

        # we don't want too much output about ERI/DFT/PE so this part is done
        # with mute/unmute
        self.ostream.mute()

        # ERI information
        eri_dict = self._init_eri(molecule, basis)
        # DFT information
        dft_dict = self._init_dft(molecule, scf_tensors, silent=True)
        # PE information
        pe_dict = self._init_pe(molecule, basis)

        self.ostream.unmute()

        profiler.stop_timer('Prep')

        profiler.start_timer('lambda')

        if self.rank == mpi_master():

            # Get overlap, MO coefficients from scf_tensors
            ovlp = scf_tensors['S']
            nocc = molecule.number_of_alpha_electrons()
            mo = scf_tensors['C_alpha']
            mo_occ = mo[:, :nocc]
            mo_vir = mo[:, nocc:]
            nocc = mo_occ.shape[1]
            nvir = mo_vir.shape[1]
            nao = mo_occ.shape[0]

            mo_energies = scf_tensors['E_alpha']
            eocc = mo_energies[:nocc]
            evir = mo_energies[nocc:]
            eo_diag = np.diag(eocc)
            ev_diag = np.diag(evir)

            x_plus_y_ao = self.cphf_results['x_plus_y_ao']
            x_minus_y_ao = self.cphf_results['x_minus_y_ao']
            dof = x_plus_y_ao.shape[0]

            fock_ao_rhs = self.cphf_results['fock_ao_rhs']

            # The density matrix; only alpha block;
            # Only works for the restricted case.
            D_occ = np.matmul(mo_occ, mo_occ.T)
            D_vir = np.matmul(mo_vir, mo_vir.T)

            # Because the excitation vector is not symmetric,
            # we need both the matrix (for omega_OO, and probably VO)
            # and its transpose (for omega_VV and OV).
            # This comes from the transformation of the 2PDM contribution
            # from MO to AO basis
            fock_ao_rhs_x_plus_y = np.zeros((dof, nao, nao))
            fock_ao_rhs_x_minus_y = np.zeros((dof, nao, nao))
            for ifock in range(dof):
                fock_ao_rhs_x_plus_y[ifock] = fock_ao_rhs[ifock + dof]
                fock_ao_rhs_x_minus_y[ifock] = fock_ao_rhs[ifock + 2 * dof]

            Fp1_vv = np.zeros((dof, nao, nao))
            Fm1_vv = np.zeros((dof, nao, nao))
            Fp2_vv = np.zeros((dof, nao, nao))
            Fm2_vv = np.zeros((dof, nao, nao))
            Fp1_oo = np.zeros((dof, nao, nao))
            Fm1_oo = np.zeros((dof, nao, nao))
            Fp2_oo = np.zeros((dof, nao, nao))
            Fm2_oo = np.zeros((dof, nao, nao))
            for s in range(dof):
                Fp1_vv[s] = 0.5 * np.linalg.multi_dot(
                    [fock_ao_rhs_x_plus_y[s].T, x_plus_y_ao[s], ovlp.T])
                Fm1_vv[s] = 0.5 * np.linalg.multi_dot(
                    [fock_ao_rhs_x_minus_y[s].T, x_minus_y_ao[s], ovlp.T])
                Fp2_vv[s] = 0.5 * np.linalg.multi_dot(
                    [fock_ao_rhs_x_plus_y[s], x_plus_y_ao[s], ovlp.T])
                Fm2_vv[s] = 0.5 * np.linalg.multi_dot(
                    [fock_ao_rhs_x_minus_y[s], x_minus_y_ao[s], ovlp.T])
                Fp1_oo[s] = 0.5 * np.linalg.multi_dot(
                    [fock_ao_rhs_x_plus_y[s], x_plus_y_ao[s].T, ovlp.T])
                Fm1_oo[s] = 0.5 * np.linalg.multi_dot(
                    [fock_ao_rhs_x_minus_y[s], x_minus_y_ao[s].T, ovlp.T])
                Fp2_oo[s] = 0.5 * np.linalg.multi_dot(
                    [fock_ao_rhs_x_plus_y[s].T, x_plus_y_ao[s].T, ovlp.T])
                Fm2_oo[s] = 0.5 * np.linalg.multi_dot(
                    [fock_ao_rhs_x_minus_y[s].T, x_minus_y_ao[s].T, ovlp.T])
        else:
            dof = None
        dof = self.comm.bcast(dof, root=mpi_master())

        # Construct fock_lambda (the lambda multipliers/cphf coefficients
        # contracted with the two-electron integrals)
        dist_cphf_ov = self.cphf_results['dist_cphf_ov']

        lambda_ao_list = []
        for x in range(dof):
            cphf_ov_x = dist_cphf_ov[x].get_full_vector(0)
            if self.rank == mpi_master():
                lambda_ao_list.append(
                    np.linalg.multi_dot(
                        [mo_occ,
                         cphf_ov_x.reshape(nocc, nvir), mo_vir.T]))

        profiler.stop_timer('lambda')

        fock_lambda = self._comp_lr_fock(lambda_ao_list, molecule, basis,
                                         eri_dict, dft_dict, pe_dict, profiler)

        profiler.start_timer('Other')

        if self.rank == mpi_master():
            # Compute the contributions from the relaxed 1PDM
            # to the omega Lagrange multipliers:
            fock_ao_lambda_np = np.zeros((dof, nao, nao))
            fock_ao_rhs_1pdm = np.zeros((dof, nao, nao))
            for ifock in range(dof):
                fock_ao_lambda_np[ifock] = fock_lambda[ifock]
                fock_ao_rhs_1pdm[ifock] = fock_ao_rhs[ifock]

            fmat = (fock_ao_lambda_np + fock_ao_lambda_np.transpose(0, 2, 1) +
                    0.5 * fock_ao_rhs_1pdm)

            omega_1pdm_2pdm_contribs = np.zeros((dof, nao, nao))
            for s in range(dof):
                omega_1pdm_2pdm_contribs[s] = 0.5 * (np.linalg.multi_dot([
                    D_vir,
                    (Fp1_vv[s] + Fm1_vv[s] - Fp2_vv[s] + Fm2_vv[s]), D_vir
                ]) + np.linalg.multi_dot([
                    D_occ,
                    (Fp1_vv[s] + Fm1_vv[s] - Fp2_vv[s] + Fm2_vv[s]), D_vir
                ]) + np.linalg.multi_dot([
                    D_occ,
                    (Fp1_vv[s] + Fm1_vv[s] - Fp2_vv[s] + Fm2_vv[s]), D_vir
                ]).T + np.linalg.multi_dot([
                    D_occ,
                    (Fp1_oo[s] + Fm1_oo[s] - Fp2_oo[s] + Fm2_oo[s]), D_occ
                ]) + 2.0 * np.linalg.multi_dot([D_occ, fmat[s], D_occ]))

            # Construct the energy-weighted one particle density matrix
            dm_oo = 0.5 * self.cphf_results['density_occ_occ']
            dm_vv = 0.5 * self.cphf_results['density_vir_vir']

            epsilon_dm_ao = np.zeros((dof, nao, nao))
            epsilon_lambda_ao = np.zeros((dof, nao, nao))

        for s in range(dof):
            cphf_ov_s = dist_cphf_ov[s].get_full_vector(0)

            if self.rank == mpi_master():
                epsilon_dm_ao[s] = np.linalg.multi_dot(
                    [mo_occ, eo_diag, dm_oo[s], mo_occ.T])
                epsilon_dm_ao[s] += np.linalg.multi_dot(
                    [mo_vir, ev_diag, dm_vv[s], mo_vir.T])
                epsilon_lambda_ao[s] = np.linalg.multi_dot(
                    [mo_occ, eo_diag,
                     cphf_ov_s.reshape(nocc, nvir), mo_vir.T])
                epsilon_dm_ao[s] += (epsilon_lambda_ao[s] +
                                     epsilon_lambda_ao[s].T)

        if self.rank == mpi_master():

            omega = -epsilon_dm_ao - omega_1pdm_2pdm_contribs

            fock_gxc_ao = self.cphf_results['fock_gxc_ao']

            if fock_gxc_ao is not None:
                factor = -0.25
                for ifock in range(dof):
                    fock_gxc_ao_np = fock_gxc_ao[2 * ifock]
                    omega[ifock] += factor * np.linalg.multi_dot(
                        [D_occ, fock_gxc_ao_np, D_occ])
        else:
            omega = None

        profiler.stop_timer('Other')

        profiler.print_timing(self.ostream)
        profiler.print_profiling_summary(self.ostream)

        self.ostream.print_info("The omega Lagrange multipliers computed in " +
                                f"{time.time() - omega_t0:.2f} sec.")
        self.ostream.print_blank()
        self.ostream.flush()

        return omega

    def print_cphf_header(self, title):
        """
        Prints information on the solver setup
        """

        self.ostream.print_blank()
        self.ostream.print_header('{:s} Setup'.format(title))
        self.ostream.print_header('=' * (len(title) + 8))
        self.ostream.print_blank()

        str_width = 70

        # print general info
        cur_str = 'Solver Type                     : '
        if self.use_subspace_solver:
            cur_str += 'Iterative Subspace Algorithm'
        else:
            cur_str += 'Conjugate Gradient'
        self.ostream.print_header(cur_str.ljust(str_width))

        cur_str = 'Max. Number of Iterations       : ' + str(self.max_iter)
        self.ostream.print_header(cur_str.ljust(str_width))
        cur_str = 'Convergence Threshold           : {:.1e}'.format(
            self.conv_thresh)
        self.ostream.print_header(cur_str.ljust(str_width))

        if self.ri_coulomb:
            cur_str = 'Resolution of the Identity      : RI-J'
            self.ostream.print_header(cur_str.ljust(str_width))

        if self._dft:
            cur_str = 'Exchange-Correlation Functional : '
            cur_str += self.xcfun.get_func_label().upper()
            self.ostream.print_header(cur_str.ljust(str_width))
            grid_level = (get_default_grid_level(self.xcfun)
                          if self.grid_level is None else self.grid_level)
            cur_str = 'Molecular Grid Level            : ' + str(grid_level)
            self.ostream.print_header(cur_str.ljust(str_width))

        self.ostream.print_blank()
        self.ostream.flush()
