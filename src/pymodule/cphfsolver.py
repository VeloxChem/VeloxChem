#
#                           VELOXCHEM 1.0-RC2
#         ----------------------------------------------------
#                     An Electronic Structure Code
#
#  Copyright © 2018-2021 by VeloxChem developers. All rights reserved.
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
#from pathlib import Path
import numpy as np
import time as tm
import sys

from .veloxchemlib import mpi_master
from .veloxchemlib import denmat
from .veloxchemlib import fockmat
from .veloxchemlib import AODensityMatrix
from .veloxchemlib import AOFockMatrix
from .veloxchemlib import ElectronRepulsionIntegralsDriver
from .outputstream import OutputStream
from .profiler import Profiler
from .distributedarray import DistributedArray
#from .signalhandler import SignalHandler
from .linearsolver import LinearSolver
from .errorhandler import assert_msg_critical
from .inputparser import parse_input
from .qqscheme import get_qq_scheme
#from .checkpoint import check_rsp_hdf5, create_hdf5, write_rsp_solution

# For PySCF integral derivatives
from .import_from_pyscf import overlap_deriv
from .import_from_pyscf import fock_deriv
from .import_from_pyscf import eri_deriv

class CphfSolver(LinearSolver):
    """
    Implements solver for the coupled-perturbed Hartree-Fock (CPHF) equations.

    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.

    ##Instance variables
    ##    - a_operator: The A operator.
    ##    - frequencies: The frequencies.
    """

    def __init__(self, comm=None, ostream=None):
        """
        Initializes CPHF solver to default setup.
        """

        if comm is None:
            comm = MPI.COMM_WORLD

        if ostream is None:
            ostream = OutputStream(sys.stdout)

        super().__init__(comm, ostream)


    def update_settings(self, orbrsp_dict, method_dict=None):
        """
        Updates response and method settings in CPHF solver.

        :param orbrsp_dict:
            The dictionary of CPHF (orbital response) settings.
        :param method_dict:
            The dictionary of method settings.
        """

        if method_dict is None:
            method_dict = {}

        super().update_settings(orbrsp_dict, method_dict)

##        rsp_keywords = {
##            'a_operator': 'str_lower',
##            'a_components': 'str_lower',
##            'b_operator': 'str_lower',
##            'b_components': 'str_lower',
##            'frequencies': 'seq_range',
##        }
##
##        parse_input(self, rsp_keywords, rsp_dict)


    def compute(self, molecule, basis, scf_tensors):
        """
        Performs CPHF calculation for a molecule and a basis set.

        :param molecule:
            The molecule.
        :param basis:
            The AO basis set.
        :param scf_tensors:
            The dictionary of tensors from converged SCF wavefunction.

        :return:
            To be determined. (Dict with RHS, CPHF coefficients, ...)
        """

        # TODO: what is needed everywhere, what only on the master node?
        density = scf_tensors['D_alpha']
        #overlap = self.scf_drv.scf_tensors['S']
        natm = molecule.number_of_atoms()
        mo = scf_tensors['C_alpha']
        nao = mo.shape[0]
        nocc = molecule.number_of_alpha_electrons()
        nvir = nao - nocc
        mo_occ = mo[:, :nocc]
        mo_vir = mo[:, nocc:]
        mo_energies = scf_tensors['E']
        eocc = mo_energies[:nocc]
        eoo = eocc.reshape(-1, 1) + eocc #ei+ej
        omega_ao = - np.linalg.multi_dot([mo_occ, np.diag(eocc), mo_occ.T])
        evir = mo_energies[nocc:]
        eov = eocc.reshape(-1, 1) - evir

        profiler = Profiler({
            'timing': self.timing,
            'profiling': self.profiling,
            'memory_profiling': self.memory_profiling,
            'memory_tracing': self.memory_tracing,
        })

        self.start_time = tm.time()

        # sanity check
        nalpha = molecule.number_of_alpha_electrons()
        nbeta = molecule.number_of_beta_electrons()
        assert_msg_critical(
            nalpha == nbeta,
            'CphfSolver: not implemented for unrestricted case')

        # preparing the CPHF RHS

        ovlp_deriv_ao = np.zeros((natm, 3, nao, nao))
        fock_deriv_ao = np.zeros((natm, 3, nao, nao))

        # import the integral derivatives
        for i in range(natm):
            ovlp_deriv_ao[i] = overlap_deriv(molecule, basis, i)
            fock_deriv_ao[i] = fock_deriv(molecule, basis, density, i)

        # transform integral derivatives to MO basis
        ovlp_deriv_ov = np.einsum('mi,xymn,na->xyia', mo_occ, ovlp_deriv_ao, mo_vir)
        ovlp_deriv_oo = np.einsum('mi,xymn,nj->xyij', mo_occ, ovlp_deriv_ao, mo_occ)
        fock_deriv_ov = np.einsum('mi,xymn,na->xyia', mo_occ, fock_deriv_ao, mo_vir)
        orben_ovlp_deriv_ov = np.einsum('i,xyia->xyia', eocc, ovlp_deriv_ov)

        # the oo part of the CPHF coefficients in AO basis,
        # transforming the oo overlap derivative back to AO basis (not equal to the initial one)
        uij_ao = np.einsum('mi,axij,nj->axmn', mo_occ, -0.5 * ovlp_deriv_oo, mo_occ).reshape((3*natm, nao, nao))
        uij_ao_list = list([uij_ao[x] for x in range(natm * 3)])

        # create AODensity and Fock matrix objects, contract with ERI
        ao_density_uij = AODensityMatrix(uij_ao_list, denmat.rest)
        fock_uij = AOFockMatrix(ao_density_uij)
        #fock_flag = fockmat.rgenjk
        #fock_uij.set_fock_type(fock_flag, 1)
        # ERI information
        eri_dict = self.init_eri(molecule, basis)
        # DFT information
        dft_dict = self.init_dft(molecule, scf_tensors)
        # PE information
        pe_dict = self.init_pe(molecule, basis)

        self.comp_lr_fock(fock_uij, ao_density_uij, molecule, basis, eri_dict, dft_dict, pe_dict, {})
        ##eri_drv = ElectronRepulsionIntegralsDriver(self.comm)
        ##screening = eri_drv.compute(get_qq_scheme(self.scf_drv.qq_type),
        ##                            self.scf_drv.eri_thresh, molecule, basis)
        ##eri_drv.compute(fock_uij, ao_density_uij, molecule, basis, screening)

        # TODO: how can this be done better?
        fock_uij_numpy = np.zeros((natm,3,nao,nao))
        for i in range(natm):
            for x in range(3):
                fock_uij_numpy[i,x] = fock_uij.to_numpy(3*i + x)

        # transform to MO basis
        fock_uij_mo = np.einsum('mi,axmn,nb->axib', mo_occ, fock_uij_numpy, mo_vir)

        # sum up the terms of the RHS
        cphf_rhs = (fock_deriv_ov - orben_ovlp_deriv_ov + 2 * fock_uij_mo).reshape(3*natm, nocc*nvir)

        # Initialize trial and sigma vectors
        self.dist_trials = None
        self.dist_sigmas = None

        # the preconditioner: 1 / (eocc - evir)
        precond = (1 / eov).reshape(nocc * nvir)

        # create list of distributed arrays for RHS
        dist_rhs = []
        # TODO: are these loops over 3*natm (all coordinates) really efficient?
        for k in range(3 * natm):
            dist_rhs.append(DistributedArray(cphf_rhs[k], self.comm))

        # setup and precondition trial vectors
        # TODO: remove renormalize=False
        dist_trials = self.setup_trials(molecule, precond, dist_rhs, renormalize=False)

        self.not_norm_dist_trials = dist_trials

        # construct the sigma (E*t) vectors
        self.build_sigmas(molecule, basis, scf_tensors, dist_trials, skip=True)

        dist_trials = self.setup_trials(molecule, precond, dist_rhs)

        # construct the sigma (E*t) vectors
        self.build_sigmas(molecule, basis, scf_tensors, dist_trials)

        # lists that will hold the solutions and residuals
        relative_residual_norm = np.ones((3*natm))
        solutions = list(np.zeros((3*natm)))
        residuals = list(np.zeros((3*natm)))

        # start iterations
        for iteration in range(self.max_iter):

            iter_start_time = tm.time()

            profiler.start_timer(iteration, 'ReducedSpace')

            # Orbital Hessian in reduced subspace
            orbhess_red = self.dist_trials.matmul_AtB(self.dist_sigmas)

            self.cur_iter = iteration
            num_vecs = self.dist_trials.shape(1)

            print(iteration, num_vecs)
            print()
            #print(self.dist_trials.data)
            #print()

            for x in range(3 * natm):
                # CPHF RHS in reduced space
                cphf_rhs_red = self.dist_trials.matmul_AtB(dist_rhs[x])

                # solve the equations exactly in the subspace
                u_red = np.linalg.solve(orbhess_red, cphf_rhs_red)

                # solution vector in full space
                u = self.dist_trials.matmul_AB_no_gather(u_red)

                # calculate residuals
                sigmas_u_red = self.dist_sigmas.matmul_AB_no_gather(u_red)
                residual = sigmas_u_red.data - dist_rhs[x].data

                # make distributed array out of residuals and current vectors
                dist_residual = DistributedArray(residual, self.comm, distribute=False)

                # calculate norms
                r_norm = np.sqrt(dist_residual.squared_norm(axis=0))
                u_norm = np.sqrt(u.squared_norm(axis=0))

                if u_norm != 0:
                    relative_residual_norm[x] = r_norm / u_norm
                else:
                    relative_residual_norm[x] = r_norm

                if relative_residual_norm[x] < self.conv_thresh:
                    print("Coordinate has converged (iteration, coordinate, shape): ", iteration, x, u.data.shape)
                    print("Solution:\n", u.data)
                    solutions[x] = u
                else:
                    residuals[x] = dist_residual


            # write to output
            if self.rank == mpi_master():
                self.ostream.print_info(
                    '{:d} trial vectors in reduced space'.format(num_vecs))
                self.ostream.print_blank()

                profiler.print_memory_subspace(
                    {
                        'dist_trials': self.dist_trials,
                        'dist_sigmas': self.dist_sigmas,
                        'precond': precond,
                        'solutions': solutions,
                        'residuals': residuals,
                    }, self.ostream)

                profiler.check_memory_usage(
                    'Iteration {:d} subspace'.format(iteration + 1))

                profiler.print_memory_tracing(self.ostream)

                ##self.print_iteration(relative_residual_norm, )

            profiler.stop_timer(iteration, 'ReducedSpace')

            # check convergence
            self.check_convergence(relative_residual_norm)

            if self.is_converged:
                break

            profiler.start_timer(iteration, 'Orthonorm.')

            # update trial vectors
            new_trials = self.setup_trials(molecule, precond, residuals, self.dist_trials)

            profiler.stop_timer(iteration, 'Orthonorm.')

            profiler.start_timer(iteration, 'FockBuild')

            # update sigma vectors
            self.build_sigmas(molecule, basis, scf_tensors, new_trials)

            profiler.stop_timer(iteration, 'FockBuild')
            profiler.check_memory_usage(
                'Iteration {:d} sigma build'.format(iteration + 1))

        # converged?
        if self.rank == mpi_master():
            self.print_convergence('Coupled-Perturbed Hartree-Fock')


        return {
            'cphf_ov': solutions,
            'cphf_rhs': cphf_rhs,
            'ovlp_deriv_oo': ovlp_deriv_oo,
            'fock_deriv_ao': fock_deriv_ao,
            'fock_uij': fock_uij_numpy,
        }



    def build_sigmas(self, molecule, basis, scf_tensors, dist_trials, skip=False):
        """
        Apply orbital Hessian matrix to a set of trial vectors.
        Appends sigma and trial vectors to member variable.

        :param molecule:
            The molecule.
        :param basis:
            The AO basis set.
        :param scf_tensors:
            The dictionary of tensors from converged SCF wavefunction.
        :param dist_trials:
            Distributed array of trial vectors
        """

        natm = molecule.number_of_atoms()
        mo = scf_tensors['C_alpha']
        nao = mo.shape[0]
        nocc = molecule.number_of_alpha_electrons()
        nvir = nao - nocc
        mo_occ = mo[:, :nocc]
        mo_vir = mo[:, nocc:]
        mo_energies = scf_tensors['E']
        eocc = mo_energies[:nocc]
        evir = mo_energies[nocc:]
        eov = eocc.reshape(-1, 1) - evir

        # ERI information
        eri_dict = self.init_eri(molecule, basis)
        # DFT information
        dft_dict = self.init_dft(molecule, scf_tensors)
        # PE information
        pe_dict = self.init_pe(molecule, basis)

        num_vecs = dist_trials.shape(1)
        vec_list = []
        # loop over columns / trial vectors
        for col in range(num_vecs):
            vec = dist_trials.get_full_vector(col).reshape(nocc, nvir)
            vec_ao = np.linalg.multi_dot([mo_occ, vec, mo_vir.T])
            vec_list.append(vec_ao)

        # TODO: remove variable
        if skip:
            self.vec_ao_list = vec_list
        
        # create density and Fock matrices, contract with two-electron integrals
        dens = AODensityMatrix(vec_list, denmat.rest)
        fock = AOFockMatrix(dens)

        fock_flag = fockmat.rgenjk
        # TODO: add dft
        for i in range(num_vecs):
            fock.set_fock_type(fock_flag, i)

        # TODO: replace empty dictionary by timing_dict
        self.comp_lr_fock(fock, dens, molecule, basis, eri_dict, dft_dict, pe_dict, {})

        # create sigma vectors
        sigmas = np.zeros((nocc * nvir, num_vecs))

        for ifock in range(num_vecs):
            fock_vec = fock.alpha_to_numpy(ifock)
            #TODO: remove this printout
            if skip and ifock == 0:
                print("Diagonal part in MO:\n")
                print(dist_trials.get_full_vector(ifock).reshape(nocc, nvir)) 
                print() 
                print("eov:\n")
                print(eov)

            cphf_mo = (- np.linalg.multi_dot([mo_occ.T, fock_vec, mo_vir])
                       - np.linalg.multi_dot([mo_vir.T, fock_vec, mo_occ]).T
                       + dist_trials.get_full_vector(ifock).reshape(nocc, nvir) * eov
                      )
            sigmas[:, ifock] = cphf_mo.reshape(nocc*nvir)

        dist_sigmas = DistributedArray(sigmas, self.comm)

        # append new sigma and trial vectors
        # TODO: create new function for this?

        if skip:
            # TODO: remove this variable
            self.sigma_init_guess = DistributedArray(dist_sigmas.data,
                                                self.comm,
                                                distribute=False)
        else:
            if self.dist_sigmas is None:
                self.dist_sigmas = DistributedArray(dist_sigmas.data,
                                                    self.comm,
                                                    distribute=False)
            else:
                self.dist_sigmas.append(dist_sigmas, axis=1)

            if self.dist_trials is None:
                self.dist_trials = DistributedArray(dist_trials.data,
                                                    self.comm,
                                                    distribute=False)
            else:
                self.dist_trials.append(dist_trials, axis=1)


    def setup_trials(self, molecule, precond, dist_rhs, old_trials=None, renormalize=True):
        """
        Set up trial vectors and apply preconditioner to them.

        :param molecule:
            The molecule.
        :param precond:
            The preconditioner.
        :param dist_rhs:
            The CPHF RHS as distributed array.
        :param old_trials:
            The previous trial vectors.
        :param renormalize:
            The flag for normalization.

        :return:
            The new set of orthonormalized trial vectors.
        """

        natm = molecule.number_of_atoms()

        # apply preconditioner to all right-hand sides and create distributed array for all trial vectors
        trials = []
        
        for k in range(3 * natm):
            v = DistributedArray(precond * dist_rhs[k].data, self.comm, distribute=False)
            norm = np.sqrt(v.squared_norm())
            
            if norm > 1e-10: # small_thresh of lrsolver
                trials.append(v.data[:])
            
        dist_trials = DistributedArray(np.array(trials).T, self.comm, distribute=False)

        # TODO: remove variable
        if old_trials is None:
            self.initial_guess = dist_trials.data
            

        if dist_trials.data.size == 0:
            dist_trials.data = np.zeros((dist_trials.shape(0), 0))

        if old_trials is not None:
            # t = t - (b (b.T t))
            bT_new_trials = old_trials.matmul_AtB_allreduce(dist_trials)
            dist_trials_proj = dist_trials.matmul_AB_no_gather(bT_new_trials)
            dist_trials.data -= dist_trials_proj.data 

        # remove linear dependencies and orthonormalize trial vectors
        if renormalize:
            if dist_trials.data.ndim > 0 and dist_trials.shape(0) > 0:
                dist_trials = self.remove_linear_dependence_half_distributed(
                                                dist_trials, self.lindep_thresh)
                dist_trials = self.orthogonalize_gram_schmidt_half_distributed(dist_trials)
                dist_trials = self.normalize_half_distributed(dist_trials)

        if self.rank == mpi_master():
            assert_msg_critical(
                dist_trials.data.size > 0,
                'CphfSolver: trial vectors are empty')

        return dist_trials


    def check_convergence(self, relative_residual_norm):
        """
        Checks convergence.

        :param relative_residual_norm:
            Relative residual norms.
        """

        if self.rank == mpi_master():
            max_residual = max(relative_residual_norm)
            if max_residual < self.conv_thresh:
                self.is_converged = True

        self.is_converged = self.comm.bcast(self.is_converged,
                                            root=mpi_master())


    # TODO: needs to be rewritten
    def print_iteration(self, relative_residual_norm):
        """
        Prints information of the iteration.

        :param relative_residual_norm:
            Relative residual norms.
        """

        width = 92
        output_header = '*** Iteration:   {} '.format(self.cur_iter + 1)
        output_header += '* Residuals (Max,Min): '
        output_header += '{:.2e} and {:.2e}'.format(
            max(relative_residual_norm.values()),
            min(relative_residual_norm.values()))
        self.ostream.print_header(output_header.ljust(width))
        self.ostream.print_blank()
        for op, freq, xv in xvs:
            ops_label = '<<{};{}>>_{:.4f}'.format(op, op, freq)
            rel_res = relative_residual_norm[(op, freq)]
            output_iter = '{:<15s}: {:15.8f} '.format(ops_label, -xv)
            output_iter += 'Residual Norm: {:.8f}'.format(rel_res)
            self.ostream.print_header(output_iter.ljust(width))
        self.ostream.print_blank()
        self.ostream.flush()



