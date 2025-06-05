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
from datetime import datetime
import numpy as np
import time as tm
import sys
import os

# TODO import VisualizationDriver
from .veloxchemlib import MolecularGrid
from .veloxchemlib import AODensityMatrix, denmat
from .veloxchemlib import ScreeningData, GpuDevices
from .veloxchemlib import compute_fock_gpu
from .veloxchemlib import compute_electric_dipole_integrals_gpu
from .veloxchemlib import compute_linear_momentum_integrals_gpu
from .veloxchemlib import compute_angular_momentum_integrals_gpu
from .veloxchemlib import integrate_fxc_fock_gpu
from .veloxchemlib import mpi_master, hartree_in_ev, rotatory_strength_in_cgs
from .griddriver import GridDriver
from .distributedarray import DistributedArray
from .subcommunicators import SubCommunicators
from .molecularorbitals import MolecularOrbitals, molorb
from .sanitychecks import dft_sanity_check
from .errorhandler import assert_msg_critical
from .inputparser import (parse_input, print_keywords, print_attributes,
                          get_random_string_parallel)
from .dftutils import get_default_grid_level
from .checkpoint import write_rsp_hdf5


class LinearSolver:
    """
    Implements linear solver.

    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.

    Instance variables
        - eri_thresh: The electron repulsion integrals screening threshold.
        - batch_size: The batch size for computation of Fock matrices.
        - dft: The flag for running DFT.
        - grid_level: The accuracy level of DFT grid.
        - xcfun: The XC functional.
        - pe: The flag for running polarizable embedding calculation.
        - pe_options: The dictionary with options for polarizable embedding.
        - electric_field: The static electric field.
        - conv_thresh: The convergence threshold for the solver.
        - max_iter: The maximum number of solver iterations.
        - cur_iter: Index of the current iteration.
        - norm_thresh: The norm threshold for a vector to be considered a zero
          vector.
        - lindep_thresh: The threshold for removing linear dependence in the
          trial vectors.
        - is_converged: The flag for convergence.
        - comm: The MPI communicator.
        - rank: The MPI rank.
        - nodes: Number of MPI processes.
        - ostream: The output stream.
        - restart: The flag for restarting from checkpoint file.
        - checkpoint_file: The name of checkpoint file.
        - force_checkpoint: The flag for writing checkpoint every iteration.
        - timing: The flag for printing timing information.
        - profiling: The flag for printing profiling information.
        - memory_profiling: The flag for printing memory usage.
        - memory_tracing: The flag for tracing memory allocation.
        - program_end_time: The end time of the program.
        - filename: The filename.
        - dist_bger: The distributed gerade trial vectors.
        - dist_bung: The distributed ungerade trial vectors.
        - dist_e2bger: The distributed gerade sigma vectors.
        - dist_e2bung: The distributed ungerade sigma vectors.
        - nonlinear: The flag for running linear solver in nonlinear response.
        - dist_fock_ger: The distributed gerade Fock matrices in MO.
        - dist_fock_ung: The distributed ungerade Fock matrices in MO.
    """

    def __init__(self, comm, ostream):
        """
        Initializes linear solver to default setup.
        """

        # ERI settings
        self.eri_thresh = 1.0e-12
        self.pair_thresh = 1.0e-10
        self.density_thresh = 1.0e-10
        self.prelink_thresh = 5.0e-6

        # dft
        self.xcfun = None
        self.grid_level = None
        self._dft = False

        # polarizable embedding
        self.potfile = None
        self.pe_options = {}
        self._pe = False

        # static electric field
        self.electric_field = None

        # solver setup
        self.conv_thresh = 1.0e-4
        self.max_iter = 150
        self.norm_thresh = None
        self.lindep_thresh = None
        self._cur_iter = 0
        self._is_converged = False

        # mpi information
        self._comm = comm
        self._rank = self._comm.Get_rank()
        self._nodes = self._comm.Get_size()

        # output stream
        self._ostream = ostream

        # restart information
        self.restart = True
        self.checkpoint_file = None
        self.force_checkpoint = False
        self.save_solutions = True

        # timing and profiling
        self.timing = False
        self.profiling = False
        self.memory_profiling = False
        self.memory_tracing = False

        # verbosity of output (1-3)
        self.print_level = 2

        # program end time for graceful exit
        self.program_end_time = None

        # filename
        self.filename = None

        # distributed arrays
        self._dist_bger = None
        self._dist_bung = None
        self._dist_e2bger = None
        self._dist_e2bung = None

        # nonlinear flag and distributed Fock matrices
        self.nonlinear = False
        self._dist_fock_ger = None
        self._dist_fock_ung = None

        # serial ratio as in Amdahl's law for estimating parallel efficiency
        self.serial_ratio = 0.1

        # input keywords
        self._input_keywords = {
            'response': {
                'eri_thresh': ('float', 'ERI screening threshold'),
                'batch_size': ('int', 'batch size for Fock build'),
                'conv_thresh': ('float', 'convergence threshold'),
                'max_iter': ('int', 'maximum number of iterations'),
                'norm_thresh': ('float', 'norm threshold for adding vector'),
                'lindep_thresh': ('float', 'threshold for linear dependence'),
                'serial_ratio': ('float', 'serial ratio as in Amdahl\'s law'),
                'restart': ('bool', 'restart from checkpoint file'),
                'filename': ('str', 'base name of output files'),
                'checkpoint_file': ('str', 'name of checkpoint file'),
                'force_checkpoint':
                    ('bool', 'flag for writing checkpoint every iteration'),
                'save_solutions': ('bool', 'save solutions to file'),
                'timing': ('bool', 'print timing information'),
                'profiling': ('bool', 'print profiling information'),
                'memory_profiling': ('bool', 'print memory usage'),
                'memory_tracing': ('bool', 'trace memory allocation'),
                'print_level': ('int', 'verbosity of output (1-3)'),
            },
            'method_settings': {
                'xcfun': ('str_upper', 'exchange-correlation functional'),
                'grid_level': ('int', 'accuracy level of DFT grid'),
                'potfile': ('str', 'potential file for polarizable embedding'),
                'electric_field': ('seq_fixed', 'static electric field'),
            },
        }

    @property
    def comm(self):
        """
        Returns the MPI communicator.
        """

        return self._comm

    @property
    def rank(self):
        """
        Returns the MPI rank.
        """

        return self._rank

    @property
    def nodes(self):
        """
        Returns the number of MPI processes.
        """

        return self._nodes

    @property
    def nnodes(self):
        """
        Returns the number of MPI processes.
        """

        return self._nodes

    @property
    def ostream(self):
        """
        Returns the output stream.
        """

        return self._ostream

    @property
    def is_converged(self):
        """
        Returns whether linear solver is converged.
        """

        return self._is_converged

    def print_keywords(self):
        """
        Prints input keywords in linear solver.
        """

        print_keywords(self._input_keywords, self.ostream)

    def print_attributes(self):
        """
        Prints attributes in linear solver.
        """

        print_attributes(self._input_keywords, self.ostream)

    def update_settings(self, rsp_dict, method_dict=None):
        """
        Updates response and method settings in linear solver.

        :param rsp_dict:
            The dictionary of response input.
        :param method_dict:
            The dictionary of method settings.
        """

        if method_dict is None:
            method_dict = {}

        rsp_keywords = {
            key: val[0] for key, val in self._input_keywords['response'].items()
        }

        parse_input(self, rsp_keywords, rsp_dict)

        if 'program_end_time' in rsp_dict:
            self.program_end_time = rsp_dict['program_end_time']
        if 'filename' in rsp_dict:
            self.filename = rsp_dict['filename']
            if 'checkpoint_file' not in rsp_dict:
                self.checkpoint_file = f'{self.filename}.rsp.h5'

        method_keywords = {
            key: val[0]
            for key, val in self._input_keywords['method_settings'].items()
        }

        parse_input(self, method_keywords, method_dict)

        dft_sanity_check(self, 'update_settings')

        self._pe_sanity_check(method_dict)

        if self.electric_field is not None:
            assert_msg_critical(
                len(self.electric_field) == 3,
                'LinearSolver: Expecting 3 values in \'electric field\' input')
            assert_msg_critical(
                not self._pe,
                'LinearSolver: \'electric field\' input is incompatible ' +
                'with polarizable embedding')
            # disable restart of calculation with static electric field since
            # checkpoint file does not contain information about the electric
            # field
            self.restart = False

    def _pe_sanity_check(self, method_dict=None):
        """
        Checks PE settings and updates relevant attributes.

        :param method_dict:
            The dicitonary of method settings.
        """

        if method_dict:
            if 'pe_options' in method_dict:
                self.pe_options = dict(method_dict['pe_options'])
            else:
                self.pe_options = {}

        if self.potfile:
            self.pe_options['potfile'] = self.potfile

        self._pe = ('potfile' in self.pe_options)

    def _init_eri(self, molecule, basis):
        """
        Initializes ERI.

        :param molecule:
            The molecule.
        :param basis:
            The AO basis set.

        :return:
            The dictionary of ERI information.
        """

        num_gpus_per_node = self._get_num_gpus_per_node()

        # TODO: remove this screening
        screening = ScreeningData(molecule, basis, num_gpus_per_node,
                                  self.pair_thresh, self.density_thresh,
                                  self.comm.Get_rank(), self.comm.Get_size())

        return {
            'screening': screening,
        }

    def _init_dft(self, molecule, scf_tensors):
        """
        Initializes DFT.

        :param molecule:
            The molecule.
        :param scf_tensors:
            The dictionary of tensors from converged SCF wavefunction.

        :return:
            The dictionary of DFT information.
        """

        if self._dft:
            grid_drv = GridDriver(self.comm)
            grid_level = (get_default_grid_level(self.xcfun)
                          if self.grid_level is None else self.grid_level)
            grid_drv.set_level(grid_level)

            num_gpus_per_node = self._get_num_gpus_per_node()

            grid_t0 = tm.time()
            molgrid = grid_drv.generate(molecule, num_gpus_per_node)
            n_grid_points = molgrid.number_of_points()
            self.ostream.print_info(
                'Molecular grid with {0:d} points generated in {1:.2f} sec.'.
                format(n_grid_points,
                       tm.time() - grid_t0))
            self.ostream.print_blank()

            # TODO: create method for broadcasting large numpy array

            if self.rank == mpi_master():
                gs_density_np = scf_tensors['D_alpha']
                naos = gs_density_np.shape[0]
            else:
                naos = None
            naos = self.comm.bcast(naos, root=mpi_master())

            if self.rank != mpi_master():
                gs_density_np = np.zeros((naos, naos))

            self.comm.Bcast(gs_density_np, root=mpi_master())
            gs_density = AODensityMatrix([gs_density_np], denmat.rest)

            dft_func_label = self.xcfun.get_func_label().upper()
        else:
            molgrid = MolecularGrid()
            gs_density = AODensityMatrix()
            dft_func_label = 'HF'

        return {
            'molgrid': molgrid,
            'gs_density': gs_density,
            'dft_func_label': dft_func_label,
        }

    def _init_pe(self, molecule, basis):
        """
        Initializes polarizable embedding.

        :param molecule:
            The molecule.
        :param basis:
            The AO basis set.

        :return:
            The dictionary of polarizable embedding information.
        """

        if self._pe:
            if self.rank == mpi_master():
                potfile = self.pe_options['potfile']
                if not Path(potfile).is_file():
                    potfile = str(
                        Path(self._filename).parent / Path(potfile).name)
                with Path(potfile).open('r') as fh:
                    potfile_text = '\n'.join(fh.readlines())
            else:
                potfile_text = None
            potfile_text = self.comm.bcast(potfile_text, root=mpi_master())
        else:
            potfile_text = ''

        return {
            'potfile_text': potfile_text,
        }

    def _read_checkpoint(self, rsp_vector_labels):
        """
        Reads distributed arrays from checkpoint file.

        :param rsp_vector_labels:
            The list of labels of vectors.
        """

        dist_arrays = [
            DistributedArray.read_from_hdf5_file(self.checkpoint_file, label,
                                                 self.comm)
            for label in rsp_vector_labels
        ]

        if self.nonlinear:
            (self._dist_bger, self._dist_bung, self._dist_e2bger,
             self._dist_e2bung, self._dist_fock_ger,
             self._dist_fock_ung) = dist_arrays
        else:
            (self._dist_bger, self._dist_bung, self._dist_e2bger,
             self._dist_e2bung) = dist_arrays

        checkpoint_text = 'Restarting from checkpoint file: '
        checkpoint_text += self.checkpoint_file
        self.ostream.print_info(checkpoint_text)
        self.ostream.print_blank()

    def _append_trial_vectors(self, bger, bung):
        """
        Appends distributed trial vectors.

        :param bger:
            The distributed gerade trial vectors.
        :param bung:
            The distributed ungerade trial vectors.
        """

        if self._dist_bger is None:
            self._dist_bger = DistributedArray(bger.data,
                                               self.comm,
                                               distribute=False)
        else:
            self._dist_bger.append(bger, axis=1)

        if self._dist_bung is None:
            self._dist_bung = DistributedArray(bung.data,
                                               self.comm,
                                               distribute=False)
        else:
            self._dist_bung.append(bung, axis=1)

    def _append_sigma_vectors(self, e2bger, e2bung):
        """
        Appends distributed sigma (E2 b) vectors.

        :param e2bger:
            The distributed gerade sigma vectors.
        :param e2bung:
            The distributed ungerade sigma vectors.
        """

        if self._dist_e2bger is None:
            self._dist_e2bger = DistributedArray(e2bger.data,
                                                 self.comm,
                                                 distribute=False)
        else:
            self._dist_e2bger.append(e2bger, axis=1)

        if self._dist_e2bung is None:
            self._dist_e2bung = DistributedArray(e2bung.data,
                                                 self.comm,
                                                 distribute=False)
        else:
            self._dist_e2bung.append(e2bung, axis=1)

    def _append_fock_matrices(self, fock_ger, fock_ung):
        """
        Appends distributed Fock matrices in MO.

        :param fock_ger:
            The distributed gerade Fock matrices in MO.
        :param fock_ung:
            The distributed ungerade Fock matrices in MO.
        """

        if self._dist_fock_ger is None:
            self._dist_fock_ger = DistributedArray(fock_ger.data,
                                                   self.comm,
                                                   distribute=False)
        else:
            self._dist_fock_ger.append(fock_ger, axis=1)

        if self._dist_fock_ung is None:
            self._dist_fock_ung = DistributedArray(fock_ung.data,
                                                   self.comm,
                                                   distribute=False)
        else:
            self._dist_fock_ung.append(fock_ung, axis=1)

    def compute(self, molecule, basis, scf_tensors, v_grad=None):
        """
        Solves for the linear equations.

        :param molecule:
            The molecule.
        :param basis:
            The AO basis.
        :param scf_tensors:
            The dictionary of tensors from converged SCF wavefunction.
        :param v_grad:
            The gradients on the right-hand side. If not provided, v_grad will
            be computed for the B operator.

        :return:
            A dictionary containing response functions, solutions, etc.
        """

        return None

    def _e2n_half_size(self,
                       vecs_ger,
                       vecs_ung,
                       molecule,
                       basis,
                       scf_tensors,
                       eri_dict,
                       dft_dict,
                       pe_dict,
                       profiler=None):
        """
        Computes the E2 b matrix vector product.

        :param vecs_ger:
            The gerade trial vectors in half-size.
        :param vecs_ung:
            The ungerade trial vectors in half-size.
        :param molecule:
            The molecule.
        :param basis:
            The AO basis set.
        :param scf_tensors:
            The dictionary of tensors from converged SCF wavefunction.
        :param eri_dict:
            The dictionary containing ERI information.
        :param dft_dict:
            The dictionary containing DFT information.
        :param pe_dict:
            The dictionary containing PE information.
        :param profiler:
            The profiler.

        :return:
            The gerade and ungerade E2 b matrix vector product in half-size.
        """

        n_ger = vecs_ger.shape(1)
        n_ung = vecs_ung.shape(1)
        n_total = n_ger + n_ung

        # prepare molecular orbitals

        if self.rank == mpi_master():
            assert_msg_critical(
                vecs_ger.data.ndim == 2 and vecs_ung.data.ndim == 2,
                'LinearSolver._e2n_half_size: '
                'invalid shape of trial vectors')

            assert_msg_critical(
                vecs_ger.shape(0) == vecs_ung.shape(0),
                'LinearSolver._e2n_half_size: '
                'inconsistent shape of trial vectors')

            mo = scf_tensors['C_alpha']
            fa = scf_tensors['F_alpha']

            nocc = molecule.number_of_alpha_electrons()
            norb = mo.shape[1]

            if getattr(self, 'core_excitation', False):
                core_exc_orb_inds = list(range(self.num_core_orbitals)) + list(
                    range(nocc, norb))
                mo_core_exc = mo[:, core_exc_orb_inds]
                fa_mo = np.linalg.multi_dot([mo_core_exc.T, fa, mo_core_exc])
            else:
                fa_mo = np.linalg.multi_dot([mo.T, fa, mo])
        else:
            mo = None
            fa = None
            nocc = None
            norb = None
            fa_mo = None

        if self.rank == mpi_master():
            dt_and_subcomm_size = []

            for subcomm_size in range(1, self.nodes + 1):
                if self.nodes % subcomm_size != 0:
                    continue

                n_subcomms = self.nodes // subcomm_size

                ave, res = divmod(n_total, n_subcomms)
                counts = [ave + 1 if p < res else ave for p in range(n_subcomms)]

                time_per_fock = self.serial_ratio + (1 - self.serial_ratio) / subcomm_size
                dt = max(counts) * time_per_fock

                dt_and_subcomm_size.append((dt, subcomm_size))

            dt_and_subcomm_size.sort()
            subcomm_size = dt_and_subcomm_size[0][1]
        else:
            subcomm_size = None
        subcomm_size = self.comm.bcast(subcomm_size, root=mpi_master())

        grps = [p // subcomm_size for p in range(self.nodes)]
        subcomm = SubCommunicators(self.comm, grps)
        local_comm = subcomm.local_comm
        cross_comm = subcomm.cross_comm

        is_local_master = (local_comm.Get_rank() == mpi_master())
        is_global_master = (self.rank == mpi_master())

        subcomm_index = None
        local_master_ranks = None

        if is_local_master:
            mo = cross_comm.bcast(mo, root=mpi_master())
            fa = cross_comm.bcast(fa, root=mpi_master())
            nocc = cross_comm.bcast(nocc, root=mpi_master())
            norb = cross_comm.bcast(norb, root=mpi_master())
            fa_mo = cross_comm.bcast(fa_mo, root=mpi_master())

            subcomm_index = cross_comm.Get_rank()
            local_master_ranks = cross_comm.allgather(self.rank)

        subcomm_index = local_comm.bcast(subcomm_index, root=mpi_master())
        local_master_ranks = local_comm.bcast(local_master_ranks, root=mpi_master())

        batch_size = self.nodes // subcomm_size

        num_batches = n_total // batch_size
        if n_total % batch_size != 0:
            num_batches += 1

        num_gpus_per_node = self._get_num_gpus_per_node()

        screening_t0 = tm.time()

        local_screening = ScreeningData(molecule, basis, num_gpus_per_node,
                                        self.pair_thresh, self.density_thresh,
                                        local_comm.Get_rank(), local_comm.Get_size())

        if profiler is not None:
            profiler.add_timing_info('Screening', tm.time() - screening_t0)

        # go through batches

        if self.rank == mpi_master():
            batch_str = f'Processing {n_total} Fock builds'
            if batch_size > 1:
                batch_str += f' on {batch_size} subcommunicators'
            batch_str += '...'
            self.ostream.print_info(batch_str)
            self.ostream.flush()

        for batch_ind in range(num_batches):

            # form density matrices

            batch_start = batch_size * batch_ind
            batch_end = min(batch_start + batch_size, n_total)

            if is_local_master:
                dks = []
                kns = []

            symm_flags = [None for idx in range(len(local_master_ranks))]
            vec_list = [None for idx in range(len(local_master_ranks))]

            for idx, local_master_rank in enumerate(local_master_ranks):

                if idx + batch_start >= batch_end:
                    break

                col = idx + batch_start

                if col < n_ger:
                    v_ger = vecs_ger.get_full_vector(col, root=local_master_rank)
                    if self.rank == local_master_rank:
                        # full-size gerade trial vector
                        vec_list[idx] = np.hstack((v_ger, v_ger))
                    # Note: antisymmetric density matrix from gerade vector
                    # due to commut_mo_density
                    symm_flags[idx] = 'antisymm'
                else:
                    v_ung = vecs_ung.get_full_vector(col - n_ger, root=local_master_rank)
                    if self.rank == local_master_rank:
                        # full-size ungerade trial vector
                        vec_list[idx] = np.hstack((v_ung, -v_ung))
                    # Note: symmetric density matrix from ungerade vector
                    # due to commut_mo_density
                    symm_flags[idx] = 'symm'

            self.comm.barrier()

            if subcomm_index + batch_start < batch_end:

                local_master_rank = local_master_ranks[subcomm_index]

                symm_flag = symm_flags[subcomm_index]

                if is_local_master:
                    vec = vec_list[subcomm_index]

                    half_size = vec.shape[0] // 2

                    if getattr(self, 'core_excitation', False):
                        kn = self.lrvec2mat(vec, nocc, norb,
                                            self.num_core_orbitals)
                        dak = self.commut_mo_density(kn, nocc,
                                                     self.num_core_orbitals)
                        core_exc_orb_inds = list(range(
                            self.num_core_orbitals)) + list(range(nocc, norb))
                        mo_core_exc = mo[:, core_exc_orb_inds]
                        dak = np.linalg.multi_dot(
                            [mo_core_exc, dak, mo_core_exc.T])
                    else:
                        kn = self.lrvec2mat(vec, nocc, norb)
                        dak = self.commut_mo_density(kn, nocc)
                        dak = np.linalg.multi_dot([mo, dak, mo.T])

                    dks.append(dak)
                    kns.append(kn)

                if is_local_master:
                    den_mat_np = dks[0]
                    naos = den_mat_np.shape[0]
                else:
                    naos = None
                naos = local_comm.bcast(naos, root=mpi_master())

                if not is_local_master:
                    den_mat_np = np.zeros((naos, naos))

                local_comm.Bcast(den_mat_np, root=mpi_master())

                dens = AODensityMatrix([den_mat_np], denmat.rest)

                # form Fock matrices

                # Note: skipping Coulomb for antisymmetric density matrix
                coulomb_coef = 0.0 if symm_flag == 'antisymm' else 2.0
                fock_mat_local = self._comp_lr_fock(dens, molecule, basis, eri_dict,
                                                    dft_dict, pe_dict, coulomb_coef,
                                                    symm_flag, local_screening, local_comm,
                                                    profiler)

                if is_local_master:
                    fock_mat = np.zeros(fock_mat_local.shape)
                else:
                    fock_mat = None

                local_comm.Reduce(fock_mat_local,
                                  fock_mat,
                                  op=MPI.SUM,
                                  root=mpi_master())

                e2_ger = None
                e2_ung = None

                fock_ger = None
                fock_ung = None

                if is_local_master:

                    col = batch_start + subcomm_index
                    if col < n_ger:
                        batch_ger, batch_ung = 1, 0
                    else:
                        batch_ger, batch_ung = 0, 1

                    e2_ger = np.zeros((half_size, batch_ger))
                    e2_ung = np.zeros((half_size, batch_ung))

                    if self.nonlinear:
                        fock_ger = np.zeros((norb**2, batch_ger))
                        fock_ung = np.zeros((norb**2, batch_ung))

                    for ifock in range(batch_ger + batch_ung):
                        fak = fock_mat

                        if getattr(self, 'core_excitation', False):
                            core_exc_orb_inds = list(range(
                                self.num_core_orbitals)) + list(range(nocc, norb))
                            mo_core_exc = mo[:, core_exc_orb_inds]
                            fak_mo = np.linalg.multi_dot(
                                [mo_core_exc.T, fak, mo_core_exc])
                        else:
                            fak_mo = np.linalg.multi_dot([mo.T, fak, mo])

                        kfa_mo = self.commut(fa_mo.T, kns[ifock])

                        fat_mo = fak_mo + kfa_mo

                        if getattr(self, 'core_excitation', False):
                            gmo = -self.commut_mo_density(fat_mo, nocc,
                                                          self.num_core_orbitals)
                            gmo_vec_halfsize = self.lrmat2vec(
                                gmo, nocc, norb, self.num_core_orbitals)[:half_size]
                        else:
                            gmo = -self.commut_mo_density(fat_mo, nocc)
                            gmo_vec_halfsize = self.lrmat2vec(gmo, nocc,
                                                              norb)[:half_size]

                        if self.nonlinear:
                            # Note: fak_mo_vec uses full MO coefficients matrix
                            # since it is only for fock_ger/fock_ung in nonlinear
                            # response
                            fak_mo_vec = np.linalg.multi_dot([mo.T, fak,
                                                              mo]).reshape(norb**2)

                        if ifock < batch_ger:
                            e2_ger[:, ifock] = -gmo_vec_halfsize
                            if self.nonlinear:
                                fock_ger[:, ifock] = fak_mo_vec
                        else:
                            e2_ung[:, ifock - batch_ger] = -gmo_vec_halfsize
                            if self.nonlinear:
                                fock_ung[:, ifock - batch_ger] = fak_mo_vec

            self.comm.barrier()

            for idx, local_master_rank in enumerate(local_master_ranks):

                if idx + batch_start >= batch_end:
                    break

                vecs_e2_ger = DistributedArray(e2_ger, self.comm, root=local_master_rank)
                vecs_e2_ung = DistributedArray(e2_ung, self.comm, root=local_master_rank)
                self._append_sigma_vectors(vecs_e2_ger, vecs_e2_ung)

                if self.nonlinear:
                    dist_fock_ger = DistributedArray(fock_ger, self.comm, root=local_master_rank)
                    dist_fock_ung = DistributedArray(fock_ung, self.comm, root=local_master_rank)
                    self._append_fock_matrices(dist_fock_ger, dist_fock_ung)

        self._append_trial_vectors(vecs_ger, vecs_ung)

        self.ostream.print_blank()

    def _comp_lr_fock(self,
                      dens,
                      molecule,
                      basis,
                      eri_dict,
                      dft_dict,
                      pe_dict,
                      prefac_coulomb,
                      flag_exchange,
                      local_screening,
                      local_comm,
                      profiler=None):
        """
        Computes Fock/Fxc matrix (2e part) for linear response calculation.

        :param fock:
            The Fock matrix (2e part).
        :param dens:
            The density matrix.
        :param molecule:
            The molecule.
        :param basis:
            The basis set.
        :param eri_dict:
            The dictionary containing ERI information.
        :param dft_dict:
            The dictionary containing DFT information.
        :param pe_dict:
            The dictionary containing PE information.
        :param profiler:
            The profiler.
        """

        #screening = eri_dict['screening']

        molgrid = dft_dict['molgrid']
        gs_density = dft_dict['gs_density']

        t0 = tm.time()
        coulomb_timing = 0.0
        exchange_timing = 0.0

        fock_mat, fock_mat_erf_k = None, None

        if self._dft:

            if self.xcfun.is_hybrid():

                if self.xcfun.is_range_separated():
                    # range-separated hybrid
                    full_k_coef = (self.xcfun.get_rs_alpha() +
                                   self.xcfun.get_rs_beta())
                    erf_k_coef = -self.xcfun.get_rs_beta()
                    omega = self.xcfun.get_rs_omega()

                    fock_mat = compute_fock_gpu(molecule, basis, dens,
                                                prefac_coulomb, full_k_coef,
                                                0.0, flag_exchange,
                                                self.eri_thresh,
                                                self.prelink_thresh, local_screening,
                                                self.rank, self.nodes)

                    coulomb_timing += np.array([
                        float(dt.split()[0])
                        for dt in local_screening.get_coulomb_time()
                    ])

                    exchange_timing += np.array([
                        float(dt.split()[0])
                        for dt in local_screening.get_exchange_time()
                    ])

                    fock_mat_erf_k = compute_fock_gpu(molecule, basis, dens,
                                                      0.0, erf_k_coef, omega,
                                                      flag_exchange,
                                                      self.eri_thresh,
                                                      self.prelink_thresh,
                                                      local_screening,
                                                      self.rank, self.nodes)

                    coulomb_timing += np.array([
                        float(dt.split()[0])
                        for dt in local_screening.get_coulomb_time()
                    ])

                    exchange_timing += np.array([
                        float(dt.split()[0])
                        for dt in local_screening.get_exchange_time()
                    ])

                else:
                    # global hybrid
                    fock_mat = compute_fock_gpu(
                        molecule, basis, dens, prefac_coulomb,
                        self.xcfun.get_frac_exact_exchange(), 0.0,
                        flag_exchange, self.eri_thresh, self.prelink_thresh,
                        local_screening, self.rank, self.nodes)

                    coulomb_timing += np.array([
                        float(dt.split()[0])
                        for dt in local_screening.get_coulomb_time()
                    ])

                    exchange_timing += np.array([
                        float(dt.split()[0])
                        for dt in local_screening.get_exchange_time()
                    ])

            else:
                # pure DFT
                fock_mat = compute_fock_gpu(molecule, basis, dens,
                                            prefac_coulomb, 0.0, 0.0,
                                            flag_exchange, self.eri_thresh,
                                            self.prelink_thresh, local_screening,
                                            self.rank, self.nodes)

                coulomb_timing += np.array([
                    float(dt.split()[0]) for dt in local_screening.get_coulomb_time()
                ])

                exchange_timing += np.array([
                    float(dt.split()[0])
                    for dt in local_screening.get_exchange_time()
                ])

        else:
            # Hartree-Fock
            fock_mat = compute_fock_gpu(molecule, basis, dens, prefac_coulomb,
                                        1.0, 0.0, flag_exchange,
                                        self.eri_thresh, self.prelink_thresh,
                                        local_screening, self.rank, self.nodes)

            coulomb_timing += np.array(
                [float(dt.split()[0]) for dt in local_screening.get_coulomb_time()])

            exchange_timing += np.array(
                [float(dt.split()[0]) for dt in local_screening.get_exchange_time()])

        all_eri_timing = local_comm.allgather(coulomb_timing + exchange_timing)

        all_eri_timing = np.array(all_eri_timing).reshape(-1)

        max_eri_timing = np.max(all_eri_timing)

        if max_eri_timing > 0.0:
            eri_load_imb = 1.0 - np.sum(all_eri_timing) / (
                all_eri_timing.size * max_eri_timing)
        else:
            eri_load_imb = 0.0

        if profiler is not None:
            profiler.add_timing_info('FockERI', tm.time() - t0)
            profiler.add_timing_info('(loadimb)', eri_load_imb)

        if self._dft:
            t0 = tm.time()

            molgrid.re_distribute_counts_and_displacements(
                local_comm.Get_rank(), local_comm.Get_size())

            # TODO: enable meta-GGA
            # Note: skipping Fxc for antisymmetric density matrix
            if flag_exchange != 'antisymm':
                integrate_fxc_fock_gpu(fock_mat, molecule, basis, dens,
                                       gs_density, molgrid,
                                       self.xcfun.get_func_label(),
                                       local_screening.get_num_gpus_per_node(),
                                       self.rank, self.nodes)

            if profiler is not None:
                profiler.add_timing_info('FockXC', tm.time() - t0)

        fock_mat_local = fock_mat.to_numpy()
        if fock_mat_erf_k is not None:
            # add contribution from range-separation
            fock_mat_local += fock_mat_erf_k.to_numpy()

        return fock_mat_local

    def _write_checkpoint(self, molecule, basis, dft_dict, pe_dict, labels):
        """
        Writes checkpoint file.

        :param molecule:
            The molecule.
        :param basis:
            The basis set.
        :param dft_dict:
            The dictionary containing DFT information.
        :param pe_dict:
            The dictionary containing PE information.
        :param labels:
            The list of labels.
        """

        if self.checkpoint_file is None and self.filename is None:
            return

        if self.checkpoint_file is None and self.filename is not None:
            self.checkpoint_file = f'{self.filename}.rsp.h5'

        t0 = tm.time()

        if self.rank == mpi_master():
            success = write_rsp_hdf5(self.checkpoint_file, [], [], molecule,
                                     basis, dft_dict, pe_dict, self.ostream)
        else:
            success = False
        success = self.comm.bcast(success, root=mpi_master())

        if success:
            if self.nonlinear:
                dist_arrays = [
                    self._dist_bger, self._dist_bung, self._dist_e2bger,
                    self._dist_e2bung, self._dist_fock_ger, self._dist_fock_ung
                ]
            else:
                dist_arrays = [
                    self._dist_bger, self._dist_bung, self._dist_e2bger,
                    self._dist_e2bung
                ]

            for dist_array, label in zip(dist_arrays, labels):
                dist_array.append_to_hdf5_file(self.checkpoint_file, label)

            checkpoint_text = 'Time spent in writing checkpoint file: '
            checkpoint_text += f'{(tm.time() - t0):.2f} sec'
            self.ostream.print_info(checkpoint_text)
            self.ostream.print_blank()

    def _graceful_exit(self, molecule, basis, dft_dict, pe_dict, labels):
        """
        Gracefully exits the program.

        :param molecule:
            The molecule.
        :param basis:
            The basis set.
        :param dft_dict:
            The dictionary containing DFT information.
        :param pe_dict:
            The dictionary containing PE information.
        :param labels:
            The list of labels.

        :return:
            The return code.
        """

        self.ostream.print_blank()
        self.ostream.print_info('Preparing for a graceful termination...')
        self.ostream.print_blank()
        self.ostream.flush()

        self._write_checkpoint(molecule, basis, dft_dict, pe_dict, labels)

        self.ostream.print_info('...done.')
        self.ostream.print_blank()
        self.ostream.print_info('Exiting program.')
        self.ostream.print_blank()
        self.ostream.flush()

        sys.exit(0)

    def _need_graceful_exit(self, next_iter_in_hours):
        """
        Checks if a graceful exit is needed.

        :param: next_iter_in_hours:
            The predicted time for the next iteration in hours.

        :return:
            True if a graceful exit is needed, False otherwise.
        """

        if self.program_end_time is not None:
            remaining_hours = (self.program_end_time -
                               datetime.now()).total_seconds() / 3600
            # exit gracefully when the remaining time is not sufficient to
            # complete the next iteration (plus 25% to be on the safe side).
            if remaining_hours < next_iter_in_hours * 1.25:
                return True
        return False

    def _print_header(self, title, nstates=None, n_freqs=None, n_points=None):
        """
        Prints linear response solver setup header to output stream.

        :param title:
            The name of the solver.
        :param nstates:
            The number of excited states (TDA/RPA).
        :param n_freqs:
            The number of frequencies (LR/CPP).
        :param n_points:
            The number of integration points (C6).
        """

        self.ostream.print_blank()
        self.ostream.print_header('{:s} Setup'.format(title))
        self.ostream.print_header('=' * (len(title) + 8))
        self.ostream.print_blank()

        str_width = 60

        # print solver-specific info

        if nstates is not None:
            cur_str = 'Number of States                : ' + str(nstates)
            self.ostream.print_header(cur_str.ljust(str_width))
        if n_freqs is not None:
            cur_str = 'Number of Frequencies           : ' + str(n_freqs)
            self.ostream.print_header(cur_str.ljust(str_width))
        if n_points is not None:
            cur_str = 'Number of Integration Points    : ' + str(n_points)
            self.ostream.print_header(cur_str.ljust(str_width))

        # print general info

        cur_str = 'Max. Number of Iterations       : ' + str(self.max_iter)
        self.ostream.print_header(cur_str.ljust(str_width))
        cur_str = 'Convergence Threshold           : {:.1e}'.format(
            self.conv_thresh)
        self.ostream.print_header(cur_str.ljust(str_width))
        cur_str = 'ERI Screening Threshold         : {:.1e}'.format(
            self.eri_thresh)
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

    def _print_iteration(self, relative_residual_norm, xvs):

        return

    def _print_convergence(self, title):
        """
        Prints information after convergence.

        :param title:
            The name of the solver.
        """

        width = 92
        output_conv = '*** '
        if self._is_converged:
            output_conv += '{:s} converged'.format(title)
        else:
            output_conv += '{:s} NOT converged'.format(title)
        output_conv += ' in {:d} iterations. '.format(self._cur_iter + 1)
        output_conv += 'Time: {:.2f} sec'.format(tm.time() - self.start_time)
        self.ostream.print_header(output_conv.ljust(width))
        self.ostream.print_blank()
        self.ostream.print_blank()
        self.ostream.flush()

    def _check_convergence(self, relative_residual_norm):
        """
        Checks convergence.

        :param relative_residual_norm:
            Relative residual norms.
        """

        self._is_converged = False

        if self.rank == mpi_master():
            max_residual = max(relative_residual_norm.values())
            if max_residual < self.conv_thresh:
                self._is_converged = True

        self._is_converged = self.comm.bcast(self._is_converged,
                                             root=mpi_master())

    def _decomp_grad(self, grad):
        """
        Decomposes gradient into gerade and ungerade parts.

        :param grad:
            The gradient.

        :return:
            A tuple containing gerade and ungerade parts of gradient.
        """

        assert_msg_critical(grad.ndim == 1,
                            'LinearSolver.decomp_grad: Expecting a 1D array')

        assert_msg_critical(
            grad.shape[0] % 2 == 0,
            'LinearSolver.decomp_grad: size of array should be even')

        half_size = grad.shape[0] // 2

        # Note: grad may be complex
        grad_T = np.zeros_like(grad)
        grad_T[:half_size] = grad[half_size:]
        grad_T[half_size:] = grad[:half_size]

        ger = 0.5 * (grad + grad_T)[:half_size]
        ung = 0.5 * (grad - grad_T)[:half_size]

        return ger.T, ung.T

    def _get_precond(self, orb_ene, nocc, norb, w):

        return None

    def _preconditioning(self, precond, v_in):

        return None

    def _precond_trials(self, vectors, precond):
        """
        Applies preconditioner to trial vectors.

        :param vectors:
            The set of vectors.
        :param precond:
            The preconditioner.

        :return:
            The preconditioned trial vectors.
        """

        return None, None

    def _setup_trials(self,
                      vectors,
                      precond,
                      dist_bger=None,
                      dist_bung=None,
                      renormalize=True):
        """
        Computes orthonormalized trial vectors.

        :param vectors:
            The set of vectors.
        :param precond:
            The preconditioner.
        :param dist_bger:
            The distributed gerade subspace.
        :param dist_bung:
            The distributed ungerade subspace.
        :param renormalize:
            The flag for normalization.

        :return:
            The orthonormalized gerade and ungerade trial vectors.
        """

        dist_new_ger, dist_new_ung = self._precond_trials(vectors, precond)

        if dist_new_ger.data.size == 0:
            dist_new_ger.data = np.zeros((dist_new_ung.shape(0), 0))

        if dist_new_ung.data.size == 0:
            dist_new_ung.data = np.zeros((dist_new_ger.shape(0), 0))

        if dist_bger is not None:
            # t = t - (b (b.T t))
            bT_new_ger = dist_bger.matmul_AtB_allreduce(dist_new_ger, 2.0)
            dist_new_ger_proj = dist_bger.matmul_AB_no_gather(bT_new_ger)
            dist_new_ger.data -= dist_new_ger_proj.data

        if dist_bung is not None:
            # t = t - (b (b.T t))
            bT_new_ung = dist_bung.matmul_AtB_allreduce(dist_new_ung, 2.0)
            dist_new_ung_proj = dist_bung.matmul_AB_no_gather(bT_new_ung)
            dist_new_ung.data -= dist_new_ung_proj.data

        if renormalize:
            if dist_new_ger.data.ndim > 0 and dist_new_ger.shape(0) > 0:
                dist_new_ger = self._remove_linear_dependence_half_size(
                    dist_new_ger, self.lindep_thresh)
                dist_new_ger = self._orthogonalize_gram_schmidt_half_size(
                    dist_new_ger)
                dist_new_ger = self._normalize_half_size(dist_new_ger)

            if dist_new_ung.data.ndim > 0 and dist_new_ung.shape(0) > 0:
                dist_new_ung = self._remove_linear_dependence_half_size(
                    dist_new_ung, self.lindep_thresh)
                dist_new_ung = self._orthogonalize_gram_schmidt_half_size(
                    dist_new_ung)
                dist_new_ung = self._normalize_half_size(dist_new_ung)

        if self.rank == mpi_master():
            assert_msg_critical(
                dist_new_ger.data.size > 0 or dist_new_ung.data.size > 0,
                'LinearSolver: trial vectors are empty')

        return dist_new_ger, dist_new_ung

    def get_prop_grad(self, operator, components, molecule, basis, scf_tensors,
                      screening):
        """
        Computes property gradients for linear response equations.

        :param operator:
            The string for the operator.
        :param components:
            The string for Cartesian components.
        :param molecule:
            The molecule.
        :param basis:
            The AO basis set.
        :param scf_tensors:
            The dictionary of tensors from converged SCF wavefunction.

        :return:
            The property gradients.
        """

        # compute 1e integral

        assert_msg_critical(
            operator in [
                'dipole', 'electric dipole', 'electric_dipole',
                'linear_momentum', 'linear momentum', 'angular_momentum',
                'angular momentum', 'magnetic dipole', 'magnetic_dipole'
            ], f'LinearSolver.get_prop_grad: unsupported operator {operator}')

        if operator in ['dipole', 'electric dipole', 'electric_dipole']:
            mu_x, mu_y, mu_z = compute_electric_dipole_integrals_gpu(
                molecule, basis, [0.0, 0.0, 0.0], screening,
                self.rank, self.nodes)

            naos = mu_x.number_of_rows()

            if self.rank == mpi_master():
                mu_x_mat_np = np.zeros((naos, naos))
                mu_y_mat_np = np.zeros((naos, naos))
                mu_z_mat_np = np.zeros((naos, naos))
            else:
                mu_x_mat_np = None
                mu_y_mat_np = None
                mu_z_mat_np = None

            self.comm.Reduce(mu_x.to_numpy(),
                             mu_x_mat_np,
                             op=MPI.SUM,
                             root=mpi_master())
            self.comm.Reduce(mu_y.to_numpy(),
                             mu_y_mat_np,
                             op=MPI.SUM,
                             root=mpi_master())
            self.comm.Reduce(mu_z.to_numpy(),
                             mu_z_mat_np,
                             op=MPI.SUM,
                             root=mpi_master())

            if self.rank == mpi_master():
                integrals = (
                    mu_x_mat_np,
                    mu_y_mat_np,
                    mu_z_mat_np,
                )
            else:
                integrals = tuple()

        elif operator in ['linear_momentum', 'linear momentum']:
            mu_x, mu_y, mu_z = compute_linear_momentum_integrals_gpu(
                molecule, basis, screening, self.rank, self.nodes)

            naos = mu_x.number_of_rows()

            if self.rank == mpi_master():
                mu_x_mat_np = np.zeros((naos, naos))
                mu_y_mat_np = np.zeros((naos, naos))
                mu_z_mat_np = np.zeros((naos, naos))
            else:
                mu_x_mat_np = None
                mu_y_mat_np = None
                mu_z_mat_np = None

            self.comm.Reduce(mu_x.to_numpy(),
                             mu_x_mat_np,
                             op=MPI.SUM,
                             root=mpi_master())
            self.comm.Reduce(mu_y.to_numpy(),
                             mu_y_mat_np,
                             op=MPI.SUM,
                             root=mpi_master())
            self.comm.Reduce(mu_z.to_numpy(),
                             mu_z_mat_np,
                             op=MPI.SUM,
                             root=mpi_master())

            if self.rank == mpi_master():
                integrals = (
                    -1.0 * mu_x_mat_np,
                    -1.0 * mu_y_mat_np,
                    -1.0 * mu_z_mat_np,
                )
            else:
                integrals = tuple()

        elif operator in ['angular_momentum', 'angular momentum']:
            mu_x, mu_y, mu_z = compute_angular_momentum_integrals_gpu(
                molecule, basis, [0.0, 0.0, 0.0], screening,
                self.rank, self.nodes)

            naos = mu_x.number_of_rows()

            if self.rank == mpi_master():
                mu_x_mat_np = np.zeros((naos, naos))
                mu_y_mat_np = np.zeros((naos, naos))
                mu_z_mat_np = np.zeros((naos, naos))
            else:
                mu_x_mat_np = None
                mu_y_mat_np = None
                mu_z_mat_np = None

            self.comm.Reduce(mu_x.to_numpy(),
                             mu_x_mat_np,
                             op=MPI.SUM,
                             root=mpi_master())
            self.comm.Reduce(mu_y.to_numpy(),
                             mu_y_mat_np,
                             op=MPI.SUM,
                             root=mpi_master())
            self.comm.Reduce(mu_z.to_numpy(),
                             mu_z_mat_np,
                             op=MPI.SUM,
                             root=mpi_master())

            if self.rank == mpi_master():
                integrals = (
                    -1.0 * mu_x_mat_np,
                    -1.0 * mu_y_mat_np,
                    -1.0 * mu_z_mat_np,
                )
            else:
                integrals = tuple()

        elif operator in ['magnetic_dipole', 'magnetic dipole']:
            mu_x, mu_y, mu_z = compute_angular_momentum_integrals_gpu(
                molecule, basis, [0.0, 0.0, 0.0], screening,
                self.rank, self.nodes)

            naos = mu_x.number_of_rows()

            if self.rank == mpi_master():
                mu_x_mat_np = np.zeros((naos, naos))
                mu_y_mat_np = np.zeros((naos, naos))
                mu_z_mat_np = np.zeros((naos, naos))
            else:
                mu_x_mat_np = None
                mu_y_mat_np = None
                mu_z_mat_np = None

            self.comm.Reduce(mu_x.to_numpy(),
                             mu_x_mat_np,
                             op=MPI.SUM,
                             root=mpi_master())
            self.comm.Reduce(mu_y.to_numpy(),
                             mu_y_mat_np,
                             op=MPI.SUM,
                             root=mpi_master())
            self.comm.Reduce(mu_z.to_numpy(),
                             mu_z_mat_np,
                             op=MPI.SUM,
                             root=mpi_master())

            if self.rank == mpi_master():
                integrals = (
                    0.5 * mu_x_mat_np,
                    0.5 * mu_y_mat_np,
                    0.5 * mu_z_mat_np,
                )
            else:
                integrals = tuple()

        # compute right-hand side

        if self.rank == mpi_master():
            indices = {'x': 0, 'y': 1, 'z': 2}
            integral_comps = [integrals[indices[p]] for p in components]

            mo = scf_tensors['C_alpha']
            nocc = molecule.number_of_alpha_electrons()
            norb = mo.shape[1]

            factor = np.sqrt(2.0)

            if getattr(self, 'core_excitation', False):
                core_exc_orb_inds = list(range(self.num_core_orbitals)) + list(
                    range(nocc, norb))
                mo_core_exc = mo[:, core_exc_orb_inds]
                matrices = [
                    factor * (-1.0) * self.commut_mo_density(
                        np.linalg.multi_dot([mo_core_exc.T, P, mo_core_exc]),
                        nocc, self.num_core_orbitals) for P in integral_comps
                ]
                gradients = tuple(
                    self.lrmat2vec(m, nocc, norb, self.num_core_orbitals)
                    for m in matrices)
            else:
                matrices = [
                    factor * (-1.0) * self.commut_mo_density(
                        np.linalg.multi_dot([mo.T, P, mo]), nocc)
                    for P in integral_comps
                ]
                gradients = tuple(
                    self.lrmat2vec(m, nocc, norb) for m in matrices)

            return gradients

        else:
            return tuple()

    def get_complex_prop_grad(self, operator, components, molecule, basis,
                              scf_tensors, screening):
        """
        Computes complex property gradients for linear response equations.

        :param operator:
            The string for the operator.
        :param components:
            The string for Cartesian components.
        :param molecule:
            The molecule.
        :param basis:
            The AO basis set.
        :param scf_tensors:
            The dictionary of tensors from converged SCF wavefunction.

        :return:
            The complex property gradients.
        """

        # compute 1e integral

        assert_msg_critical(
            operator in [
                'dipole', 'electric dipole', 'electric_dipole',
                'linear_momentum', 'linear momentum', 'angular_momentum',
                'angular momentum', 'magnetic dipole', 'magnetic_dipole'
            ],
            f'LinearSolver.get_complex_prop_grad: unsupported operator {operator}'
        )

        if operator in ['dipole', 'electric dipole', 'electric_dipole']:
            mu_x, mu_y, mu_z = compute_electric_dipole_integrals_gpu(
                molecule, basis, [0.0, 0.0, 0.0], screening,
                self.rank, self.nodes)

            naos = mu_x.number_of_rows()

            if self.rank == mpi_master():
                mu_x_mat_np = np.zeros((naos, naos))
                mu_y_mat_np = np.zeros((naos, naos))
                mu_z_mat_np = np.zeros((naos, naos))
            else:
                mu_x_mat_np = None
                mu_y_mat_np = None
                mu_z_mat_np = None

            self.comm.Reduce(mu_x.to_numpy(),
                             mu_x_mat_np,
                             op=MPI.SUM,
                             root=mpi_master())
            self.comm.Reduce(mu_y.to_numpy(),
                             mu_y_mat_np,
                             op=MPI.SUM,
                             root=mpi_master())
            self.comm.Reduce(mu_z.to_numpy(),
                             mu_z_mat_np,
                             op=MPI.SUM,
                             root=mpi_master())

            if self.rank == mpi_master():
                integrals = (
                    mu_x_mat_np + 0j,
                    mu_y_mat_np + 0j,
                    mu_z_mat_np + 0j,
                )
            else:
                integrals = tuple()

        elif operator in ['linear_momentum', 'linear momentum']:
            mu_x, mu_y, mu_z = compute_linear_momentum_integrals_gpu(
                molecule, basis, screening, self.rank, self.nodes)

            naos = mu_x.number_of_rows()

            if self.rank == mpi_master():
                mu_x_mat_np = np.zeros((naos, naos))
                mu_y_mat_np = np.zeros((naos, naos))
                mu_z_mat_np = np.zeros((naos, naos))
            else:
                mu_x_mat_np = None
                mu_y_mat_np = None
                mu_z_mat_np = None

            self.comm.Reduce(mu_x.to_numpy(),
                             mu_x_mat_np,
                             op=MPI.SUM,
                             root=mpi_master())
            self.comm.Reduce(mu_y.to_numpy(),
                             mu_y_mat_np,
                             op=MPI.SUM,
                             root=mpi_master())
            self.comm.Reduce(mu_z.to_numpy(),
                             mu_z_mat_np,
                             op=MPI.SUM,
                             root=mpi_master())

            if self.rank == mpi_master():
                integrals = (
                    -1j * mu_x_mat_np,
                    -1j * mu_y_mat_np,
                    -1j * mu_z_mat_np,
                )
            else:
                integrals = tuple()

        elif operator in ['angular_momentum', 'angular momentum']:
            mu_x, mu_y, mu_z = compute_angular_momentum_integrals_gpu(
                molecule, basis, [0.0, 0.0, 0.0], screening,
                self.rank, self.nodes)

            naos = mu_x.number_of_rows()

            if self.rank == mpi_master():
                mu_x_mat_np = np.zeros((naos, naos))
                mu_y_mat_np = np.zeros((naos, naos))
                mu_z_mat_np = np.zeros((naos, naos))
            else:
                mu_x_mat_np = None
                mu_y_mat_np = None
                mu_z_mat_np = None

            self.comm.Reduce(mu_x.to_numpy(),
                             mu_x_mat_np,
                             op=MPI.SUM,
                             root=mpi_master())
            self.comm.Reduce(mu_y.to_numpy(),
                             mu_y_mat_np,
                             op=MPI.SUM,
                             root=mpi_master())
            self.comm.Reduce(mu_z.to_numpy(),
                             mu_z_mat_np,
                             op=MPI.SUM,
                             root=mpi_master())

            if self.rank == mpi_master():
                integrals = (
                    -1j * mu_x_mat_np,
                    -1j * mu_y_mat_np,
                    -1j * mu_z_mat_np,
                )
            else:
                integrals = tuple()

        elif operator in ['magnetic_dipole', 'magnetic dipole']:
            mu_x, mu_y, mu_z = compute_angular_momentum_integrals_gpu(
                molecule, basis, [0.0, 0.0, 0.0], screening,
                self.rank, self.nodes)

            naos = mu_x.number_of_rows()

            if self.rank == mpi_master():
                mu_x_mat_np = np.zeros((naos, naos))
                mu_y_mat_np = np.zeros((naos, naos))
                mu_z_mat_np = np.zeros((naos, naos))
            else:
                mu_x_mat_np = None
                mu_y_mat_np = None
                mu_z_mat_np = None

            self.comm.Reduce(mu_x.to_numpy(),
                             mu_x_mat_np,
                             op=MPI.SUM,
                             root=mpi_master())
            self.comm.Reduce(mu_y.to_numpy(),
                             mu_y_mat_np,
                             op=MPI.SUM,
                             root=mpi_master())
            self.comm.Reduce(mu_z.to_numpy(),
                             mu_z_mat_np,
                             op=MPI.SUM,
                             root=mpi_master())

            if self.rank == mpi_master():
                integrals = (
                    0.5j * mu_x_mat_np,
                    0.5j * mu_y_mat_np,
                    0.5j * mu_z_mat_np,
                )
            else:
                integrals = tuple()

        # compute right-hand side

        if self.rank == mpi_master():
            indices = {'x': 0, 'y': 1, 'z': 2}
            integral_comps = [integrals[indices[p]] for p in components]

            mo = scf_tensors['C_alpha']
            nocc = molecule.number_of_alpha_electrons()
            norb = mo.shape[1]

            factor = np.sqrt(2.0)
            matrices = [
                factor * (-1.0) *
                self.commut_mo_density(np.linalg.multi_dot([mo.T, P, mo]), nocc)
                for P in integral_comps
            ]

            gradients = tuple(self.lrmat2vec(m, nocc, norb) for m in matrices)
            return gradients

        else:
            return tuple()

    @staticmethod
    def commut_mo_density(A, nocc, num_core_orbitals=None):
        """
        Commutes matrix A and MO density

        :param A:
            Matrix A.
        :param nocc:
            Number of occupied orbitals.

        :return:
            A D_mo - D_mo A
        """

        # | 0    -A_ov |
        # | A_vo  0    |

        mat = np.zeros(A.shape, dtype=A.dtype)

        if num_core_orbitals is not None and num_core_orbitals > 0:
            mat[:num_core_orbitals,
                num_core_orbitals:] = -A[:num_core_orbitals, num_core_orbitals:]
            mat[num_core_orbitals:, :num_core_orbitals] = A[
                num_core_orbitals:, :num_core_orbitals]

        else:
            mat[:nocc, nocc:] = -A[:nocc, nocc:]
            mat[nocc:, :nocc] = A[nocc:, :nocc]

        return mat

    @staticmethod
    def commut(A, B):
        """
        Commutes two matricies A and B

        :param A:
            Matrix A.
        :param B:
            Matrix B.

        :return:
            AB - BA
        """

        return np.matmul(A, B) - np.matmul(B, A)

    @staticmethod
    def lrvec2mat(vec, nocc, norb, num_core_orbitals=None):
        """
        Converts vectors to matrices.

        :param vec:
            The vectors.
        :param nocc:
            Number of occupied orbitals.
        :param norb:
            Number of orbitals.

        :return:
            The matrices.
        """

        nvir = norb - nocc

        if num_core_orbitals is not None and num_core_orbitals > 0:
            n_ov = num_core_orbitals * nvir
            mat = np.zeros((num_core_orbitals + nvir, num_core_orbitals + nvir),
                           dtype=vec.dtype)
            # excitation and de-excitation
            mat[:num_core_orbitals, num_core_orbitals:] = vec[:n_ov].reshape(
                num_core_orbitals, nvir)
            mat[num_core_orbitals:, :num_core_orbitals] = vec[n_ov:].reshape(
                num_core_orbitals, nvir).T

        else:
            n_ov = nocc * nvir
            mat = np.zeros((norb, norb), dtype=vec.dtype)
            # excitation and de-excitation
            mat[:nocc, nocc:] = vec[:n_ov].reshape(nocc, nvir)
            mat[nocc:, :nocc] = vec[n_ov:].reshape(nocc, nvir).T

        return mat

    @staticmethod
    def lrmat2vec(mat, nocc, norb, num_core_orbitals=None):
        """
        Converts matrices to vectors.

        :param mat:
            The matrices.
        :param nocc:
            Number of occupied orbitals.
        :param norb:
            Number of orbitals.

        :return:
            The vectors.
        """

        nvir = norb - nocc

        if num_core_orbitals is not None and num_core_orbitals > 0:
            n_ov = num_core_orbitals * nvir
            vec = np.zeros(n_ov * 2, dtype=mat.dtype)
            # excitation and de-excitation
            vec[:n_ov] = mat[:num_core_orbitals,
                             num_core_orbitals:].reshape(n_ov)
            vec[n_ov:] = mat[num_core_orbitals:, :num_core_orbitals].T.reshape(
                n_ov)

        else:
            n_ov = nocc * nvir
            vec = np.zeros(n_ov * 2, dtype=mat.dtype)
            # excitation and de-excitation
            vec[:n_ov] = mat[:nocc, nocc:].reshape(n_ov)
            vec[n_ov:] = mat[nocc:, :nocc].T.reshape(n_ov)

        return vec

    @staticmethod
    def remove_linear_dependence(basis, threshold):
        """
        Removes linear dependence in a set of vectors.

        :param basis:
            The set of vectors.
        :param threshold:
            The threshold for removing linear dependence.

        :return:
            The new set of vectors.
        """

        Sb = np.matmul(basis.T, basis)
        l, T = np.linalg.eigh(Sb)
        b_norm = np.sqrt(Sb.diagonal())
        mask = l > b_norm * threshold
        return np.matmul(basis, T[:, mask])

    def _remove_linear_dependence_half_size(self, basis, threshold):
        """
        Removes linear dependence in a set of symmetrized vectors.

        :param basis:
            The set of upper parts of symmetrized vectors.
        :param threshold:
            The threshold for removing linear dependence.

        :return:
            The new set of vectors.
        """

        half_Sb = basis.matmul_AtB(basis)

        if self.rank == mpi_master():
            Sb = 2.0 * half_Sb
            l, T = np.linalg.eigh(Sb)
            b_norm = np.sqrt(Sb.diagonal())
            mask = l > b_norm * threshold
            Tmask = T[:, mask].copy()
        else:
            Tmask = None
        Tmask = self.comm.bcast(Tmask, root=mpi_master())

        return basis.matmul_AB_no_gather(Tmask)

    @staticmethod
    def orthogonalize_gram_schmidt(tvecs):
        """
        Applies modified Gram Schmidt orthogonalization to trial vectors.

        :param tvecs:
            The trial vectors.

        :return:
            The orthogonalized trial vectors.
        """

        if tvecs.shape[1] > 0:

            f = 1.0 / np.linalg.norm(tvecs[:, 0])
            tvecs[:, 0] *= f

            for i in range(1, tvecs.shape[1]):
                for j in range(i):
                    f = np.dot(tvecs[:, i], tvecs[:, j]) / np.dot(
                        tvecs[:, j], tvecs[:, j])
                    tvecs[:, i] -= f * tvecs[:, j]
                f = 1.0 / np.linalg.norm(tvecs[:, i])
                tvecs[:, i] *= f

        return tvecs

    def _orthogonalize_gram_schmidt_half_size(self, tvecs):
        """
        Applies modified Gram Schmidt orthogonalization to trial vectors.

        :param tvecs:
            The trial vectors.

        :return:
            The orthogonalized trial vectors.
        """

        invsqrt2 = 1.0 / np.sqrt(2.0)

        if tvecs.shape(1) > 0:

            n2 = tvecs.dot(0, tvecs, 0)
            f = invsqrt2 / np.sqrt(n2)
            tvecs.data[:, 0] *= f

            for i in range(1, tvecs.shape(1)):
                for j in range(i):
                    dot_ij = tvecs.dot(i, tvecs, j)
                    dot_jj = tvecs.dot(j, tvecs, j)
                    f = dot_ij / dot_jj
                    tvecs.data[:, i] -= f * tvecs.data[:, j]
                n2 = tvecs.dot(i, tvecs, i)
                f = invsqrt2 / np.sqrt(n2)
                tvecs.data[:, i] *= f

        return tvecs

    @staticmethod
    def normalize(vecs):
        """
        Normalizes vectors by dividing by vector norm.

        :param vecs:
            The vectors.

        :param Retruns:
            The normalized vectors.
        """

        if len(vecs.shape) != 1:
            for vec in range(vecs.shape[1]):
                invnorm = 1.0 / np.linalg.norm(vecs[:, vec])
                vecs[:, vec] *= invnorm
        else:
            invnorm = 1.0 / np.linalg.norm(vecs)
            vecs *= invnorm

        return vecs

    def _normalize_half_size(self, vecs):
        """
        Normalizes half-sized vectors by dividing by vector norm.

        :param vecs:
            The half-sized vectors.

        :param Retruns:
            The normalized vectors.
        """

        invsqrt2 = 1.0 / np.sqrt(2.0)

        invnorm = invsqrt2 / vecs.norm(axis=0)

        vecs.data *= invnorm

        return vecs

    @staticmethod
    def construct_ediag_sdiag_half(orb_ene, nocc, norb, num_core_orbitals=None):
        """
        Gets the upper half of E0 and S0 diagonal elements as arrays.

        :param orb_ene:
            Orbital energies.
        :param nocc:
            Number of occupied orbitals.
        :param norb:
            Number of orbitals.

        :return:
            The upper half of E0 and S0 diagonal elements as numpy arrays.
        """

        nvir = norb - nocc

        if num_core_orbitals is not None and num_core_orbitals > 0:
            n_ov = num_core_orbitals * nvir
            eocc = orb_ene[:num_core_orbitals]
        else:
            n_ov = nocc * nvir
            eocc = orb_ene[:nocc]

        evir = orb_ene[nocc:]

        ediag = (-eocc.reshape(-1, 1) + evir).reshape(n_ov)
        sdiag = np.ones(ediag.shape)

        return ediag, sdiag

    def get_nto(self, t_mat, mo_occ, mo_vir):
        """
        Gets the natural transition orbitals.

        :param t_mat:
            The excitation vector in matrix form (N_occ x N_virt).
        :param mo_occ:
            The MO coefficients of occupied orbitals.
        :param mo_vir:
            The MO coefficients of virtual orbitals.

        :return:
            The NTOs as a MolecularOrbitals object.
        """

        # SVD
        u_mat, s_diag, vh_mat = np.linalg.svd(t_mat, full_matrices=True)
        lam_diag = s_diag**2

        # holes in increasing order of lambda
        # particles in decreasing order of lambda
        nto_occ = np.flip(np.matmul(mo_occ, u_mat), axis=1)
        nto_vir = np.matmul(mo_vir, vh_mat.T)

        # NTOs including holes and particles
        nto_orbs = np.concatenate((nto_occ, nto_vir), axis=1)

        nto_lam = np.zeros(nto_orbs.shape[1])
        nocc = nto_occ.shape[1]
        for i_nto in range(lam_diag.size):
            nto_lam[nocc - 1 - i_nto] = -lam_diag[i_nto]
            nto_lam[nocc + i_nto] = lam_diag[i_nto]

        nto_mo = MolecularOrbitals.create_nto([nto_orbs], [nto_lam],
                                              molorb.rest)

        return nto_mo

    def write_nto_cubes(self,
                        cubic_grid,
                        molecule,
                        basis,
                        root,
                        nto_mo,
                        nto_pairs=None,
                        nto_thresh=0.1):
        """
        Writes cube files for natural transition orbitals.

        :param cubic_grid:
            The cubic grid.
        :param molecule:
            The molecule.
        :param basis:
            The AO basis set.
        :param root:
            The index of the root (0-based).
        :param nto_mo:
            The NTO coefficients (2D array).
        :param nto_pairs:
            The number of NTO pairs.
        :param nto_thresh:
            The threshold for writing NTO to cube file.

        :return:
            Values of lambda and names of cube files.
        """

        filenames = []

        if self.filename is not None:
            base_fname = self.filename
        else:
            name_string = get_random_string_parallel(self.comm)
            base_fname = 'vlx_' + name_string

        vis_drv = VisualizationDriver(self.comm)

        if getattr(self, 'core_excitation', False):
            nocc = self.num_core_orbitals
        else:
            nocc = molecule.number_of_alpha_electrons()
        nvir = nto_mo.number_of_mos() - nocc
        lam_diag = nto_mo.occa_to_numpy()[nocc:nocc + min(nocc, nvir)]

        for i_nto in range(lam_diag.size):
            if lam_diag[i_nto] < nto_thresh:
                continue

            if (nto_pairs is not None) and (i_nto >= nto_pairs):
                continue

            self.ostream.print_info('  lambda: {:.4f}'.format(lam_diag[i_nto]))

            # hole
            ind_occ = nocc - i_nto - 1
            vis_drv.compute(cubic_grid, molecule, basis, nto_mo, ind_occ,
                            'alpha')

            if self.rank == mpi_master():
                occ_cube_name = '{:s}_S{:d}_NTO_H{:d}.cube'.format(
                    base_fname, root + 1, i_nto + 1)
                vis_drv.write_data(occ_cube_name, cubic_grid, molecule, 'nto',
                                   ind_occ, 'alpha')
                filenames.append(occ_cube_name)

                self.ostream.print_info(
                    '    Cube file (hole)     : {:s}'.format(occ_cube_name))
                self.ostream.flush()

            # electron
            ind_vir = nocc + i_nto
            vis_drv.compute(cubic_grid, molecule, basis, nto_mo, ind_vir,
                            'alpha')

            if self.rank == mpi_master():
                vir_cube_name = '{:s}_S{:d}_NTO_P{:d}.cube'.format(
                    base_fname, root + 1, i_nto + 1)
                vis_drv.write_data(vir_cube_name, cubic_grid, molecule, 'nto',
                                   ind_vir, 'alpha')
                filenames.append(vir_cube_name)

                self.ostream.print_info(
                    '    Cube file (particle) : {:s}'.format(vir_cube_name))
                self.ostream.flush()

        self.ostream.print_blank()

        return lam_diag, filenames

    def get_detach_attach_densities(self, z_mat, y_mat, mo_occ, mo_vir):
        """
        Gets the detachment and attachment densities.

        :param z_mat:
            The excitation vector in matrix form (N_occ x N_virt).
        :param y_mat:
            The de-excitation vector in matrix form (N_occ x N_virt).
        :param mo_occ:
            The MO coefficients of occupied orbitals.
        :param mo_vir:
            The MO coefficients of virtual orbitals.

        :return:
            The detachment and attachment densities.
        """

        dens_D = -np.linalg.multi_dot([mo_occ, z_mat, z_mat.T, mo_occ.T])
        dens_A = np.linalg.multi_dot([mo_vir, z_mat.T, z_mat, mo_vir.T])

        if y_mat is not None:
            dens_D += np.linalg.multi_dot([mo_occ, y_mat, y_mat.T, mo_occ.T])
            dens_A -= np.linalg.multi_dot([mo_vir, y_mat.T, y_mat, mo_vir.T])

        return dens_D, dens_A

    def write_detach_attach_cubes(self, cubic_grid, molecule, basis, root,
                                  dens_DA):
        """
        Writes cube files for detachment and attachment densities.

        :param cubic_grid:
            The cubic grid.
        :param molecule:
            The molecule.
        :param basis:
            The AO basis set.
        :param root:
            The index of the root (0-based).
        :param dens_DA:
            The AODensityMatrix object containing detachment and attachment
            densities.

        :return:
            The list containing the names of the cube files.
        """

        filenames = []

        if self.filename is not None:
            base_fname = self.filename
        else:
            name_string = get_random_string_parallel(self.comm)
            base_fname = 'vlx_' + name_string

        vis_drv = VisualizationDriver(self.comm)

        vis_drv.compute(cubic_grid, molecule, basis, dens_DA, 0, 'alpha')

        if self.rank == mpi_master():
            detach_cube_name = '{:s}_S{:d}_detach.cube'.format(
                base_fname, root + 1)
            vis_drv.write_data(detach_cube_name, cubic_grid, molecule,
                               'detachment', 0, 'alpha')
            filenames.append(detach_cube_name)

            self.ostream.print_info(
                '  Cube file (detachment) : {:s}'.format(detach_cube_name))
            self.ostream.flush()

        vis_drv.compute(cubic_grid, molecule, basis, dens_DA, 1, 'alpha')

        if self.rank == mpi_master():
            attach_cube_name = '{:s}_S{:d}_attach.cube'.format(
                base_fname, root + 1)
            vis_drv.write_data(attach_cube_name, cubic_grid, molecule,
                               'attachment', 1, 'alpha')
            filenames.append(attach_cube_name)

            self.ostream.print_info(
                '  Cube file (attachment) : {:s}'.format(attach_cube_name))
            self.ostream.flush()

        self.ostream.print_blank()

        return filenames

    def get_excitation_details(self, eigvec, nocc, nvir, coef_thresh=0.2):
        """
        Get excitation details.

        :param eigvec:
            The eigenvectors.
        :param nocc:
            The number of occupied molecular orbitals.
        :param nvir:
            The number of virtual molecular orbitals.
        :param coef_thresh:
            The threshold for including transitions.

        :return:
            The excitation details as a list of strings.
        """

        n_ov = nocc * nvir
        assert_msg_critical(
            eigvec.size == n_ov or eigvec.size == n_ov * 2,
            'LinearSolver.get_excitation_details: Inconsistent size')

        excitations = []
        de_excitations = []

        for i in range(nocc):
            if getattr(self, 'core_excitation', False):
                homo_str = f'core_{i+1}'
            else:
                homo_str = 'HOMO' if i == nocc - 1 else f'HOMO-{nocc-1-i}'

            for a in range(nvir):
                lumo_str = 'LUMO' if a == 0 else f'LUMO+{a}'

                ia = i * nvir + a

                exc_coef = eigvec[ia]
                if abs(exc_coef) > coef_thresh:
                    excitations.append((
                        abs(exc_coef),
                        f'{homo_str:<8s} -> {lumo_str:<8s} {exc_coef:10.4f}',
                    ))

                if eigvec.size == n_ov * 2:
                    de_exc_coef = eigvec[n_ov + ia]
                    if abs(de_exc_coef) > coef_thresh:
                        de_excitations.append((
                            abs(de_exc_coef),
                            f'{homo_str:<8s} <- {lumo_str:<8s} {de_exc_coef:10.4f}',
                        ))

        excitation_details = []
        for exc in sorted(excitations, reverse=True):
            excitation_details.append(exc[1])
        for de_exc in sorted(de_excitations, reverse=True):
            excitation_details.append(de_exc[1])

        return excitation_details

    def _print_transition_dipoles(self, title, trans_dipoles):
        """
        Prints transition dipole moments to output stream.

        :param title:
            The title to be printed to the output stream.
        :param trans_dipoles:
            The transition dipole moments.
        """

        spin_str = 'S'

        valstr = title
        self.ostream.print_header(valstr.ljust(92))
        self.ostream.print_header(('-' * len(valstr)).ljust(92))
        valstr = '                     '
        valstr += '{:>13s}{:>13s}{:>13s}'.format('X', 'Y', 'Z')
        self.ostream.print_header(valstr.ljust(92))
        for s, r in enumerate(trans_dipoles):
            valstr = 'Excited State {:>5s}: '.format(spin_str + str(s + 1))
            valstr += '{:13.6f}{:13.6f}{:13.6f}'.format(r[0], r[1], r[2])
            self.ostream.print_header(valstr.ljust(92))
        self.ostream.print_blank()

    def _print_absorption(self, title, results):
        """
        Prints absorption to output stream.

        :param title:
            The title to be printed to the output stream.
        :param results:
            The dictionary containing response results.
        """

        spin_str = 'S'

        valstr = title
        self.ostream.print_header(valstr.ljust(92))
        self.ostream.print_header(('-' * len(valstr)).ljust(92))
        for s, e in enumerate(results['eigenvalues']):
            valstr = 'Excited State {:>5s}: '.format(spin_str + str(s + 1))
            valstr += '{:15.8f} a.u. '.format(e)
            valstr += '{:12.5f} eV'.format(e * hartree_in_ev())
            f = results['oscillator_strengths'][s]
            valstr += '    Osc.Str. {:9.4f}'.format(f)
            self.ostream.print_header(valstr.ljust(92))
        self.ostream.print_blank()

    def _print_ecd(self, title, results):
        """
        Prints electronic circular dichroism to output stream.

        :param title:
            The title to be printed to the output stream.
        :param results:
            The dictionary containing response results.
        """

        spin_str = 'S'

        valstr = title
        self.ostream.print_header(valstr.ljust(92))
        self.ostream.print_header(('-' * len(valstr)).ljust(92))
        for s, R in enumerate(results['rotatory_strengths']):
            valstr = 'Excited State {:>5s}: '.format(spin_str + str(s + 1))
            valstr += '    Rot.Str. '
            valstr += f'{(R / rotatory_strength_in_cgs()):13.6f} a.u.'
            valstr += f'{R:11.4f} [10**(-40) cgs]'
            self.ostream.print_header(valstr.ljust(92))
        self.ostream.print_blank()

    def _print_excitation_details(self, title, results):
        """
        Prints excitation details.

        :param title:
            The title to be printed to the output stream.
        :param results:
            The dictionary containing response results.
        """

        self.ostream.print_header(title.ljust(92))
        self.ostream.print_blank()

        nstates = results['eigenvalues'].size
        excitation_details = results['excitation_details']

        for s in range(nstates):
            valstr = 'Excited state {}'.format(s + 1)
            self.ostream.print_header(valstr.ljust(92))
            self.ostream.print_header(('-' * len(valstr)).ljust(92))

            for exc_str in excitation_details[s]:
                self.ostream.print_header(exc_str.ljust(92))
            self.ostream.print_blank()

        self.ostream.flush()

    def _print_excited_state_absorption(self, title, results):
        """
        Prints excited state absorption.

        :param title:
            The title to be printed to the output stream.
        :param results:
            The dictionary containing response results.
        """

        self.ostream.print_blank()
        self.ostream.print_header(title.ljust(92))
        self.ostream.print_blank()

        for excitation in results['esa_results']:
            from_state = excitation['from_state']
            to_state = excitation['to_state']
            ene_in_ev = excitation['excitation_energy'] * hartree_in_ev()
            osc_str = excitation['oscillator_strength']
            trans_dipole = excitation['transition_dipole']

            valstr = f'{from_state:<4s}-> {to_state:<4s}'
            valstr += f' {ene_in_ev:10.5f} eV'
            valstr += f'     Osc.Str. {osc_str:8.4f}'
            valstr += '    Trans.Dipole '
            for i in range(3):
                valstr += f'{trans_dipole[i]:8.4f}'
            self.ostream.print_header(valstr.ljust(92))

        self.ostream.print_blank()
        self.ostream.flush()

    def _get_num_gpus_per_node(self):
        """
        Gets number of GPUs per MPI process.

        :return:
            The number of GPUs per MPI process.
        """

        if 'VLX_NUM_GPUS_PER_NODE' in os.environ:
            num_gpus_per_node = int(os.environ['VLX_NUM_GPUS_PER_NODE'])
        else:
            devices = GpuDevices()
            num_gpus_per_node = devices.get_number_devices()
            if 'SLURM_NTASKS_PER_NODE' in os.environ:
                num_gpus_per_node //= int(os.environ['SLURM_NTASKS_PER_NODE'])

        assert_msg_critical(
            num_gpus_per_node > 0,
            'LinearSolver: Invalid number of GPUs per MPI process')

        return num_gpus_per_node
