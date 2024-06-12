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

# TODO: remove commented out code;
from mpi4py import MPI
import numpy as np
import time as tm
import sys

from .veloxchemlib import mpi_master
from .veloxchemlib import denmat
from .veloxchemlib import fockmat
from .veloxchemlib import AODensityMatrix
from .veloxchemlib import AOFockMatrix
from .veloxchemlib import ElectronRepulsionIntegralsDriver
from .veloxchemlib import XCMolecularHessian
from .outputstream import OutputStream
from .profiler import Profiler
from .distributedarray import DistributedArray
from .linearsolver import LinearSolver
from .errorhandler import assert_msg_critical
from .inputparser import parse_input
from .qqscheme import get_qq_scheme
from .batchsize import get_batch_size
from .batchsize import get_number_of_batches
from .dftutils import get_default_grid_level
from scipy.sparse import linalg

# For PySCF integral derivatives
from .import_from_pyscf import overlap_deriv
from .import_from_pyscf import fock_deriv
from .import_from_pyscf import vxc_deriv
from .import_from_pyscf import eri_deriv

class CphfSolver(LinearSolver):
    """
    Implements solver for the coupled-perturbed Hartree-Fock (CPHF) equations.

    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.

    Instance variables
        - use_subspace_solver: flag to use subspace solver
          instead of conjugate gradient.
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

        self.use_subspace_solver = True
        self.print_residuals = False
        self.max_iter = 25

        self._input_keywords['orbitalresponse'] = {
            'use_subspace_solver': ('bool', 'subspace or conjugate algorithm'),
            'print_residuals': ('bool', 'print iteration to output'),
            'max_iter': ('int', 'maximum number of iterations')
        }

    def update_settings(self, cphf_dict, method_dict=None):
        """
        Updates response and method settings in CPHF solver.

        :param cphf_dict:
            The dictionary of CPHF (orbital response) settings.
        :param method_dict:
            The dictionary of method settings.
        """

        super().update_settings(cphf_dict, method_dict)

        cphf_keywords = {
            key: val[0] for key, val in self._input_keywords['orbitalresponse'].items()
        }

        parse_input(self, cphf_keywords, cphf_dict)

    def compute(self, molecule, basis, scf_tensors, *args):
        """
        Performs CPHF calculation for a molecule and a basis set.

        :param molecule:
            The molecule.
        :param basis:
            The AO basis set.
        :param scf_tensors:
            The dictionary of tensors from converged SCF wavefunction.
        :param *args:
            Additional arguments, same as in compute_rhs.

        :return:
            A dictionary containing the RHS and solution (ov block)
            of the CPHF equations, the derivative of the AO Fock matrix,
            partial derivative of the overlap matrix (oo block),
            and an auxiliary Fock matrix (oo block of the CPHF coefficients
            contracted with the two-electron integrals).
        """

        self.profiler = Profiler({
            'timing': self.timing,
            'profiling': self.profiling,
            'memory_profiling': self.memory_profiling,
            'memory_tracing': self.memory_tracing,
        })

        if self.norm_thresh is None:
            self.norm_thresh = self.conv_thresh * 1.0e-6
        if self.lindep_thresh is None:
            self.lindep_thresh = self.conv_thresh * 1.0e-2

        # sanity check
        nalpha = molecule.number_of_alpha_electrons()
        nbeta = molecule.number_of_beta_electrons()
        assert_msg_critical(
            nalpha == nbeta,
            'CphfSolver: not implemented for unrestricted case')

        if self.rank == mpi_master():
            if self._dft:
                self.print_cphf_header('Coupled-Perturbed Kohn-Sham Solver')
            else:
                self.print_cphf_header('Coupled-Perturbed Hartree-Fock Solver')

        if self.use_subspace_solver:
            self.cphf_results = self.compute_subspace_solver(molecule, basis, scf_tensors, *args)
        else:
            self.cphf_results = self.compute_conjugate_gradient(molecule, basis, scf_tensors, *args)


    def compute_subspace_solver(self, molecule, basis, scf_tensors, *args):
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

        self.start_time = tm.time()

        # ERI information
        eri_dict = self._init_eri(molecule, basis)
        # DFT information
        dft_dict = self._init_dft(molecule, scf_tensors)
        # PE information
        pe_dict = self._init_pe(molecule, basis)
        # Timing information
        timing_dict = {}

        if self.rank == mpi_master():
            mo_energies = scf_tensors['E']
            # nmo is sometimes different than nao (because of linear
            # dependencies which get removed during SCF)
            nmo = mo_energies.shape[0]
            nocc = molecule.number_of_alpha_electrons()
            nvir = nmo - nocc
            nao = scf_tensors['C_alpha'].shape[0]
            eocc = mo_energies[:nocc]
            evir = mo_energies[nocc:]
            eov = eocc.reshape(-1, 1) - evir
        else:
            mo_energies = None
            eov = None
            nao = None

        mo_energies = self.comm.bcast(mo_energies, root=mpi_master())
        eov = self.comm.bcast(eov, root=mpi_master())
        nao = self.comm.bcast(nao, root=mpi_master())
        nmo = mo_energies.shape[0]
        nocc = molecule.number_of_alpha_electrons()
        natm =  molecule.number_of_atoms()
        nvir = nmo - nocc

        cphf_rhs_dict = self.compute_rhs(molecule, basis, scf_tensors, *args)

        if self.rank == mpi_master():
            # get rhs, find out how many degrees of freedom, and reshape
            cphf_rhs = cphf_rhs_dict['cphf_rhs'] #.reshape(3*natm, nocc*nvir)
            dof = cphf_rhs.shape[0]
            cphf_rhs = cphf_rhs.reshape(dof, nocc * nvir)
        else:
            cphf_rhs = None
            dof = None

        dof = self.comm.bcast(dof, root=mpi_master())

        # Initialize trial and sigma vectors
        self.dist_trials = None
        self.dist_sigmas = None

        # the preconditioner: 1 / (eocc - evir)
        precond = (1 / eov).reshape(nocc * nvir)
        dist_precond = DistributedArray(precond, self.comm)

        # create list of distributed arrays for RHS
        dist_rhs = []
        # TODO: are these loops over all degrees of freedom really efficient?
        for k in range(dof):
            if self.rank == mpi_master():
                cphf_rhs_k = cphf_rhs[k]
            else:
                cphf_rhs_k = None
            dist_rhs.append(DistributedArray(cphf_rhs_k, self.comm))

        # setup and precondition trial vectors
        dist_trials = self.setup_trials(molecule, dist_precond, dist_rhs)

        # construct the sigma (E*t) vectors
        self.build_sigmas(molecule, basis, scf_tensors, dist_trials, eri_dict,
                          dft_dict, pe_dict, timing_dict)

        # lists that will hold the solutions and residuals
        solutions = list(np.zeros((dof)))
        residuals = list(np.zeros((dof)))
        relative_residual_norm = np.ones((dof))

        # start iterations
        for iteration in range(self.max_iter):

            iter_start_time = tm.time()

            self.profiler.set_timing_key(f'Iteration {iteration+1}')
            self.profiler.start_timer('ReducedSpace')

            # Orbital Hessian in reduced subspace
            orbhess_red = self.dist_trials.matmul_AtB(self.dist_sigmas)

            self._cur_iter = iteration
            num_vecs = self.dist_trials.shape(1)

            for x in range(dof):
                if (iteration > 0 and
                        relative_residual_norm[x] < self.conv_thresh):
                    continue

                # CPHF RHS in reduced space
                cphf_rhs_red = self.dist_trials.matmul_AtB(dist_rhs[x])

                if self.rank == mpi_master():
                    # solve the equations exactly in the subspace
                    u_red = np.linalg.solve(orbhess_red, cphf_rhs_red)
                else:
                    u_red = None
                u_red = self.comm.bcast(u_red, root=mpi_master())

                # solution vector in full space
                u = self.dist_trials.matmul_AB_no_gather(u_red)

                # calculate residuals
                sigmas_u_red = self.dist_sigmas.matmul_AB_no_gather(u_red)
                residual = sigmas_u_red.data - dist_rhs[x].data

                # make distributed array out of residuals and current vectors
                dist_residual = DistributedArray(residual, self.comm,
                                                 distribute=False)

                # calculate norms
                r_norm = np.sqrt(dist_residual.squared_norm(axis=0))
                u_norm = np.sqrt(u.squared_norm(axis=0))

                if u_norm != 0:
                    relative_residual_norm[x] = r_norm / u_norm
                else:
                    relative_residual_norm[x] = r_norm

                if relative_residual_norm[x] < self.conv_thresh:
                    solutions[x] = u
                else:
                    residuals[x] = dist_residual

            # write to output
            if self.rank == mpi_master():
                self.ostream.print_info(
                    '{:d} trial vectors in reduced space'.format(num_vecs))
                self.ostream.print_blank()

                self.profiler.print_memory_subspace(
                    {
                        'dist_trials': self.dist_trials,
                        'dist_sigmas': self.dist_sigmas,
                        'precond': precond,
                        'solutions': solutions,
                        'residuals': residuals,
                    }, self.ostream)

                self.profiler.check_memory_usage(
                    'Iteration {:d} subspace'.format(iteration + 1))

                self.profiler.print_memory_tracing(self.ostream)

                if self.print_residuals:
                    self.print_iteration(relative_residual_norm, molecule)

            self.profiler.stop_timer('ReducedSpace')

            # check convergence
            self.check_convergence(relative_residual_norm)

            if self.is_converged:
                break

            self.profiler.start_timer('Orthonorm.')

            # update trial vectors
            new_trials = self.setup_trials(molecule, dist_precond,
                                           residuals, self.dist_trials)

            self.profiler.stop_timer('Orthonorm.')

            self.profiler.start_timer('FockBuild')

            # update sigma vectors
            self.build_sigmas(molecule, basis, scf_tensors, new_trials, eri_dict,
                              dft_dict, pe_dict, timing_dict)

            self.profiler.stop_timer('FockBuild')
            self.profiler.check_memory_usage(
                'Iteration {:d} sigma build'.format(iteration + 1))

        self.profiler.print_timing(self.ostream)
        self.profiler.print_profiling_summary(self.ostream)

        self.profiler.check_memory_usage('End of ComputeCPHF iterative subspace algorithm')
        self.profiler.print_memory_usage(self.ostream)

        # converged?
        if self.rank == mpi_master():
            if self._dft:
                self._print_convergence('Coupled-Perturbed Kohn-Sham')
            else:
                self._print_convergence('Coupled-Perturbed Hartree-Fock')

        # transform Distributed arrays into numpy arrays
        # TODO: decide whether using Distributed arrays in
        # ScfHessianDriver is more convenient/efficient etc.

            cphf_ov = np.zeros((dof, nocc, nvir))

        for i in range(dof):
            sol = solutions[i].get_full_vector(0)
            if self.rank == mpi_master():
                cphf_ov[i] = sol.reshape(nocc, nvir)

        # merge the rhs dict with the solution
        if self.rank == mpi_master():
            return {**cphf_rhs_dict,
                'cphf_ov': cphf_ov,
            }
        return None

    # To avoid init_eri, init_dft at each iteration, we do it once and pass it on to
    # build_sigmas.
    def build_sigmas(self, molecule, basis, scf_tensors, dist_trials, eri_dict,
                    dft_dict, pe_dict, timing_dict):
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

        if self.rank == mpi_master():
            natm = molecule.number_of_atoms()
            mo = scf_tensors['C_alpha']
            nao = mo.shape[0]
            mo_energies = scf_tensors['E']
            # nmo is sometimes different than nao (because of linear
            # dependencies which get removed during SCF)
            nmo = mo_energies.shape[0]
            nocc = molecule.number_of_alpha_electrons()
            nvir = nmo - nocc
            mo_occ = mo[:, :nocc]
            mo_vir = mo[:, nocc:]
            eocc = mo_energies[:nocc]
            evir = mo_energies[nocc:]
            eov = eocc.reshape(-1, 1) - evir
        else:
            nao = None

        num_vecs = dist_trials.shape(1)

        batch_size = get_batch_size(self.batch_size, num_vecs, nao, self.comm)
        num_batches = get_number_of_batches(num_vecs, batch_size, self.comm)

        if self.rank == mpi_master():
            batch_str = 'Processing Fock builds...'
            batch_str += ' (batch size: {:d})'.format(batch_size)
            self.ostream.print_info(batch_str)

        for batch_ind in range(num_batches):
            self.ostream.print_info('  batch {}/{}'.format(
                batch_ind + 1, num_batches))
            self.ostream.flush()

            batch_start = batch_size * batch_ind
            batch_end = min(batch_start + batch_size, num_vecs)

            if self.rank == mpi_master():
                vec_list = []
            # loop over columns / trial vectors
            for col in range(batch_start, batch_end):
                vec = dist_trials.get_full_vector(col)

                if self.rank == mpi_master():
                    vec = vec.reshape(nocc, nvir)
                    vec_ao = np.linalg.multi_dot([mo_occ, vec, mo_vir.T])
                    vec_list.append(vec_ao)

            if self.rank == mpi_master():
                # create density matrices
                dens = AODensityMatrix(vec_list, denmat.rest)
            else:
                dens = AODensityMatrix()

            dens.broadcast(self.rank, self.comm)

            # create Fock matrices and contract with two-electron integrals
            fock = AOFockMatrix(dens)

            self._comp_lr_fock(fock, dens, molecule, basis, eri_dict,
                              dft_dict, pe_dict, self.profiler)

            sigmas = None

            if self.rank == mpi_master():
                # create sigma vectors
                sigmas = np.zeros((nocc * nvir, num_vecs))

            for ifock in range(batch_start, batch_end):
                vec = dist_trials.get_full_vector(ifock)

                if self.rank == mpi_master():
                    fock_vec = fock.alpha_to_numpy(ifock)
                    cphf_mo = (
                        - np.linalg.multi_dot([mo_occ.T, fock_vec, mo_vir])
                        - np.linalg.multi_dot([mo_vir.T, fock_vec, mo_occ]).T
                        + vec.reshape(nocc,nvir) * eov
                      )
                    sigmas[:, ifock] = cphf_mo.reshape(nocc*nvir)

            dist_sigmas = DistributedArray(sigmas, self.comm)

            # append new sigma and trial vectors
            # TODO: create new function for this?
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


    def setup_trials(self, molecule, dist_precond, dist_rhs, old_trials=None,
                     renormalize=True):
        """
        Set up trial vectors and apply preconditioner to them.

        :param molecule:
            The molecule.
        :param dist_precond:
            The preconditioner as a distributed array.
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

        # apply preconditioner to all right-hand sides and create
        # distributed array for all trial vectors
        trials = []

        n = len(dist_rhs)

        for k in range(n):
            v = DistributedArray(dist_precond.data * dist_rhs[k].data,
                                 self.comm, distribute=False)
            norm = np.sqrt(v.squared_norm())

            if norm > 1e-10: # small_thresh of lrsolver
                trials.append(v.data[:])

        dist_trials = DistributedArray(np.array(trials).T,
                                       self.comm, distribute=False)

        if dist_trials.data.size == 0:
            dist_trials.data = np.zeros((dist_trials.shape(0), 0))

        if old_trials is not None:
            # t = t - (b (b.T t))
            bT_new_trials = old_trials.matmul_AtB_allreduce(dist_trials)
            dist_trials_proj = old_trials.matmul_AB_no_gather(bT_new_trials)
            dist_trials.data -= dist_trials_proj.data

        # remove linear dependencies and orthonormalize trial vectors
        if renormalize:
            if dist_trials.data.ndim > 0 and dist_trials.shape(0) > 0:
                dist_trials = self._remove_linear_dependence_half_size(
                                                dist_trials, self.lindep_thresh)
                dist_trials = (
                    self._orthogonalize_gram_schmidt_half_size(dist_trials)
                    )
                dist_trials = self._normalize_half_size(dist_trials)

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
                self._is_converged = True

        self._is_converged = self.comm.bcast(self.is_converged,
                                            root=mpi_master())


    def compute_conjugate_gradient(self, molecule, basis, scf_tensors, *args):
        """
        Computes the coupled-perturbed Hartree-Fock (CPHF) coefficients.

        :param molecule:
            The molecule.
        :param basis:
            The AO basis set.
        :param scf_tensors
            The SCF tensors.

        :return:
            A dictionary containing the RHS and solution (ov block)
            of the CPHF equations, the derivative of the AO Fock matrix,
            partial derivative of the overlap matrix (oo block),
            and an auxiliary Fock matrix (oo block of the CPHF coefficients
            contracted with the two-electron integrals).
        """

        self.start_time = tm.time()

        # ERI information
        eri_dict = self._init_eri(molecule, basis)
        # DFT information
        dft_dict = self._init_dft(molecule, scf_tensors)
        # PE information
        pe_dict = self._init_pe(molecule, basis)
        # Timing information
        timing_dict = {}

        self.profiler.set_timing_key('CPHF RHS')
        self.profiler.start_timer('CPHF RHS')

        if self.rank == mpi_master():
            mo_energies = scf_tensors['E']
            # nmo is sometimes different than nao (because of linear
            # dependencies which get removed during SCF)
            nmo = mo_energies.shape[0]
            nao = scf_tensors['C_alpha'].shape[0]
            nocc = molecule.number_of_alpha_electrons()
            nvir = nmo - nocc
            eocc = mo_energies[:nocc]
            evir = mo_energies[nocc:]
            eov = eocc.reshape(-1, 1) - evir
        else:
            mo_energies = None
            eov = None
            nao = None

        mo_energies = self.comm.bcast(mo_energies, root=mpi_master())
        eov = self.comm.bcast(eov, root=mpi_master())
        nao = self.comm.bcast(nao, root=mpi_master())
        nmo = mo_energies.shape[0]
        nocc = molecule.number_of_alpha_electrons()
        natm =  molecule.number_of_atoms()
        nvir = nmo - nocc

        cphf_rhs_dict = self.compute_rhs(molecule, basis, scf_tensors, *args)

        if self.rank == mpi_master():
            cphf_rhs = cphf_rhs_dict['cphf_rhs'] #.reshape(3*natm, nocc*nvir)
            dof = cphf_rhs.shape[0]
        else:
            cphf_rhs = None

        self.profiler.stop_timer('CPHF RHS')

        cphf_rhs = self.comm.bcast(cphf_rhs, root=mpi_master())

        # Solve the CPHF equations using conjugate gradient (cg)
        cphf_ov = self.solve_cphf_cg(molecule, basis, scf_tensors,
                                  cphf_rhs # TODO: possibly change the shape
                                  )

        self.profiler.print_timing(self.ostream)
        self.profiler.print_profiling_summary(self.ostream)

        self.profiler.check_memory_usage('End of ComputeCPHF CG algorithm')
        self.profiler.print_memory_usage(self.ostream)

        if self.rank == mpi_master():
            self.ostream.print_blank()
            if self._dft:
                self._print_convergence('Coupled-Perturbed Kohn-Sham')
            else:
                self._print_convergence('Coupled-Perturbed Hartree-Fock')

            # merge the rhs dict with the solution
            cphf_ov_dict = {**cphf_rhs_dict, 'cphf_ov': cphf_ov}

            return cphf_ov_dict

        return None

    def solve_cphf_cg(self, molecule, basis, scf_tensors, cphf_rhs):
        """
        Solves the CPHF equations using conjugate gradient
        for all atomic coordinates to obtain the ov block
        of the CPHF coefficients.

        :param molecule:
            The molecule.
        :param basis:
            The AO basis set.
        :param scf_tensors:
            The tensors from the converged SCF calculation.
        :param cphf_rhs:
            The right-hand side of the CPHF equations for all atomic coordinates.

        :returns:
            The ov block of the CPHF coefficients.
        """

        # ERI information
        eri_dict = self._init_eri(molecule, basis)
        # DFT information
        dft_dict = self._init_dft(molecule, scf_tensors)
        # PE information
        pe_dict = self._init_pe(molecule, basis)
        # Timing information
        timing_dict = {}

        nocc = molecule.number_of_alpha_electrons()
        natm = molecule.number_of_atoms()
        # degrees of freedom from rhs (can be different from number of atoms)
        dof = cphf_rhs.shape[0]

        if self.rank == mpi_master():
            mo = scf_tensors['C']
            nao = mo.shape[0]
            mo_energies = scf_tensors['E']

            mo_occ = mo[:, :nocc].copy()
            mo_vir = mo[:, nocc:].copy()
            nvir = mo_vir.shape[1]

            eocc = mo_energies[:nocc]
            evir = mo_energies[nocc:]
            eov = eocc.reshape(-1, 1) - evir
        else:
            nvir = None
            eov = None
        nvir = self.comm.bcast(nvir, root=mpi_master())
        eov = self.comm.bcast(eov, root=mpi_master())

        # Calculate the initial guess for the CPHF coefficients given by
        # the RHS divided by orbital-energy differences
        cphf_guess = cphf_rhs / eov

        if self.rank == mpi_master():
            # Create AODensityMatrix object from CPHF guess in AO
            cphf_ao = np.zeros((dof, nao, nao))
            cphf_ao = np.array([
                np.linalg.multi_dot([mo_occ, cphf_guess[x], mo_vir.T])
                for x in range(dof)
            ])

            cphf_ao_list = list([cphf_ao[x] for x in range(dof)])
            # create AODensityMatrix object
            ao_density_cphf = AODensityMatrix(cphf_ao_list, denmat.rest)
        else:
            ao_density_cphf = AODensityMatrix()
        ao_density_cphf.broadcast(self.rank, self.comm)


        # Create a Fock Matrix Object (initialized with zeros)
        fock_cphf = AOFockMatrix(ao_density_cphf)

        # Matrix-vector product of orbital Hessian with trial vector
        def cphf_matvec(v):
            """
            Function to carry out matrix multiplication of CPHF coefficient
            vector with orbital Hessian matrix.
            """

            self.profiler.set_timing_key('Iter ' + str(self._cur_iter) + 'CG')
            self.profiler.start_timer('Iter ' + str(self._cur_iter) + 'CG')

            # Create AODensityMatrix object from lambda in AO
            if self.rank == mpi_master():
                cphf_ao = np.zeros((dof, nao, nao))
                for i in range(dof):
                    cphf_ao[i] = np.linalg.multi_dot([
                        mo_occ, v.reshape(dof,nocc,nvir)[i], mo_vir.T
                        ])
                cphf_ao_list = list([cphf_ao[x] for x in range(dof)])
                ao_density_cphf = AODensityMatrix(cphf_ao_list, denmat.rest)
            else:
                ao_density_cphf = AODensityMatrix()
            ao_density_cphf.broadcast(self.rank, self.comm)

            self._comp_lr_fock(fock_cphf, ao_density_cphf, molecule,
                              basis, eri_dict, dft_dict, pe_dict, self.profiler)

            # Transform to MO basis (symmetrized w.r.t. occ. and virt.)
            # and add diagonal part
            if self.rank == mpi_master():
                cphf_mo = np.zeros((dof, nocc, nvir))
                for i in range(dof):
                    fock_cphf_numpy = fock_cphf.to_numpy(i)
                    cphf_mo[i] = ( - np.linalg.multi_dot([mo_occ.T,
                                                        fock_cphf_numpy,
                                                        mo_vir])
                                  - np.linalg.multi_dot([mo_vir.T,
                                                        fock_cphf_numpy,
                                                        mo_occ]).T
                                    )
                cphf_mo +=  v.reshape(dof, nocc, nvir) * eov
            else:
                cphf_mo = None

            cphf_mo = self.comm.bcast(cphf_mo, root=mpi_master())

            self.profiler.stop_timer('Iter ' + str(self._cur_iter) + 'CG')

            self.profiler.check_memory_usage(
                'CG Iteration {:d}'.format(self._cur_iter + 1))

            self.profiler.print_memory_tracing(self.ostream)

            # increase iteration counter every time this function is called
            self._cur_iter += 1

            if self.rank == mpi_master():
                if self.print_residuals:
                    residual_norms = np.zeros(dof)
                    for i in range(dof):
                        residual_norms[i] = np.linalg.norm(cphf_mo[i] - cphf_rhs[i])
                    self.print_iteration(residual_norms, molecule)

            return cphf_mo.reshape(dof * nocc * nvir)

        # Matrix-vector product for preconditioner using the
        # inverse of the diagonal (i.e. eocc - evir)
        def precond_matvec(v):
            """
            Function that defines the matrix-vector product
            required by the pre-conditioner for the conjugate gradient.
            It is an approximation for the inverse of matrix A in Ax = b.
            """
            current_v = v.reshape(dof, nocc, nvir)
            M_dot_v = current_v / eov

            return M_dot_v.reshape(dof * nocc * nvir)

        # 5) Define the linear operators and run conjugate gradient
        LinOp = linalg.LinearOperator((dof * nocc * nvir,
                                       dof * nocc * nvir),
                                       matvec=cphf_matvec)
        PrecondOp = linalg.LinearOperator((dof * nocc * nvir,
                                           dof * nocc * nvir),
                                           matvec=precond_matvec)

        b = cphf_rhs.reshape(dof * nocc * nvir)
        x0 = cphf_guess.reshape(dof * nocc * nvir)

        cphf_coefficients_ov, cg_conv = linalg.cg(A=LinOp,
                                                b=b,
                                                x0=x0,
                                                M=PrecondOp,
                                                rtol=self.conv_thresh,
                                                atol=0,
                                                maxiter=self.max_iter)

        self._is_converged = (cg_conv == 0)

        return cphf_coefficients_ov.reshape(dof, nocc, nvir)

    def compute_rhs(self, molecule, basis, scf_tensors, scf_drv):
        """
        Computes the right hand side for the CPHF equations for
        the analytical Hessian, all atomic coordinates.

        :param molecule:
            The molecule.
        :param basis:
            The AO basis set.
        :param scf_tensors:
            The tensors from the converged SCF calculation.
        :param scf_drv:
            The scf_driver

        :returns:
            The RHS of the CPHF equations.
        """

        return None

    def print_iteration_cg(self, residual_norm):
        """
        Prints information of the iteration.

        :param residual_norm:
            Residual norm.
        """
        width = 92
        output_header = '*** Iteration:   {:2d} '.format(self._cur_iter)
        output_header += '  * Residual Norm: '
        output_header += '{:.5e}'.format(residual_norm)
        self.ostream.print_header(output_header.ljust(width))
        self.ostream.flush()



    def print_iteration(self, relative_residual_norm, molecule):
        """
        Prints information of the iteration.

        :param relative_residual_norm:
            Relative residual norms.
        :param molecule:
            The molecule.
        """
        natm = molecule.number_of_atoms()
        atom_symbols = molecule.get_labels()

        width = 92
        output_header = '*** Iteration:   {} '.format(self._cur_iter + 1)
        output_header += '* Residuals (Max,Min): '
        output_header += '{:.2e} and {:.2e}'.format(
            max(relative_residual_norm),
            min(relative_residual_norm))
        self.ostream.print_header(output_header.ljust(width))
        self.ostream.print_blank()
        coord = 'xyz'

        # TODO: change power based on the convergence threshold
        if self.conv_thresh < 1e-4:
            power = int(str(self.conv_thresh).split("-")[1]) + 1
        else:
            power = 5
        residual_norm_string = 'Residual Norm: {:.%df}' % (power)

        for i in range(natm):
            for x in range(3):
                coord_label = '{:16s}'.format('   '+str(i+1)+' '
                                             +atom_symbols[i]+'('+coord[x]+')')
                rel_res = relative_residual_norm[3*i+x]
                output_iter = residual_norm_string.format(rel_res)
                self.ostream.print_line(coord_label + output_iter)
        self.ostream.print_blank()
        self.ostream.flush()

    def print_cphf_header(self, title):
        """
        Prints information on the solver setup
        """

        return None

