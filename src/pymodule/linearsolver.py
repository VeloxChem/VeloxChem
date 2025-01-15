#
#                              VELOXCHEM
#         ----------------------------------------------------
#                     An Electronic Structure Code
#
#  Copyright © 2018-2024 by VeloxChem developers. All rights reserved.
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

from datetime import datetime
import numpy as np
import time as tm
import math
import sys

from .veloxchemlib import T4CScreener
from .veloxchemlib import MolecularGrid, XCIntegrator
from .veloxchemlib import mpi_master, hartree_in_ev
from .veloxchemlib import rotatory_strength_in_cgs
from .veloxchemlib import make_matrix, mat_t
from .matrix import Matrix
from .distributedarray import DistributedArray
from .subcommunicators import SubCommunicators
from .fockdriver import FockDriver
from .griddriver import GridDriver
from .molecularorbitals import MolecularOrbitals, molorb
from .visualizationdriver import VisualizationDriver
from .profiler import Profiler
from .oneeints import (compute_electric_dipole_integrals,
                       compute_linear_momentum_integrals,
                       compute_angular_momentum_integrals)
from .sanitychecks import dft_sanity_check, pe_sanity_check
from .errorhandler import assert_msg_critical
from .inputparser import (parse_input, print_keywords, print_attributes,
                          get_random_string_parallel)
from .dftutils import get_default_grid_level, print_libxc_reference
from .checkpoint import write_rsp_hdf5
from .batchsize import get_batch_size
from .batchsize import get_number_of_batches


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
        self.eri_thresh = 1.0e-15
        self.batch_size = None

        # dft
        self.xcfun = None
        self.grid_level = None
        self._dft = False

        # polarizable embedding
        self.potfile = None
        self.pe_options = {}
        self._pe = False
        self.embedding_options = None
        self._embedding_drv = None

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

        self._debug = False
        self._block_size_factor = 8
        self._xcfun_ldstaging = 1024

        # serial ratio as in Amdahl's law for estimating parallel efficiency
        self.serial_ratio = 0.05
        self.use_subcomms = False

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
                'use_subcomms': ('bool', 'use subcommunicators for Fock build'),
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
                '_debug': ('bool', 'print debug info'),
                '_block_size_factor': ('int', 'block size factor for ERI'),
                '_xcfun_ldstaging': ('int', 'max batch size for DFT grid'),
            },
            'method_settings': {
                'xcfun': ('str_upper', 'exchange-correlation functional'),
                'grid_level': ('int', 'accuracy level of DFT grid'),
                'potfile': ('str', 'potential file for polarizable embedding'),
                'embedding_options':
                    ('dict', 'dictionary to set the embedding options'),
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

        pe_sanity_check(self, method_dict)

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

        self._print_mem_debug_info('before screener')
        if self.rank == mpi_master():
            screening = T4CScreener()
            screening.partition(basis, molecule, 'eri')
        else:
            screening = None
        self._print_mem_debug_info('after  screener')
        screening = self.comm.bcast(screening, root=mpi_master())
        self._print_mem_debug_info('after  bcast screener')

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
            print_libxc_reference(self.xcfun, self.ostream)

            grid_drv = GridDriver(self.comm)
            grid_level = (get_default_grid_level(self.xcfun)
                          if self.grid_level is None else self.grid_level)
            grid_drv.set_level(grid_level)

            grid_t0 = tm.time()
            molgrid = grid_drv.generate(molecule, self._xcfun_ldstaging)
            n_grid_points = molgrid.number_of_points()
            self.ostream.print_info(
                'Molecular grid with {0:d} points generated in {1:.2f} sec.'.
                format(n_grid_points,
                       tm.time() - grid_t0))
            self.ostream.print_blank()

            if self.rank == mpi_master():
                # Note: make gs_density a tuple
                gs_density = (scf_tensors['D_alpha'].copy(),)
            else:
                gs_density = None
            gs_density = self.comm.bcast(gs_density, root=mpi_master())

            dft_func_label = self.xcfun.get_func_label().upper()
        else:
            molgrid = MolecularGrid()
            gs_density = None
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
            assert_msg_critical(
                self.embedding_options['settings']['embedding_method'] == 'PE',
                'PolarizableEmbedding: Invalid embedding_method. Only PE is supported.'
            )

            from .embedding import PolarizableEmbeddingLRS

            self._embedding_drv = PolarizableEmbeddingLRS(
                molecule=molecule,
                ao_basis=basis,
                options=self.embedding_options,
                comm=self.comm)

            # TODO: print PyFraME info

            pot_info = 'Reading polarizable embedding potential: {}'.format(
                self.pe_options['potfile'])
            self.ostream.print_info(pot_info)
            self.ostream.print_blank()

            with open(str(self.pe_options['potfile']), 'r') as f_pot:
                potfile_text = '\n'.join(f_pot.readlines())
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

    def _print_mem_debug_info(self, label):
        """
        Prints memory debug information.

        :param label:
            The label of memory debug information.
        """

        if self._debug:
            profiler = Profiler()
            self.ostream.print_info(f'==DEBUG==   available memory {label}: ' +
                                    profiler.get_available_memory())
            self.ostream.flush()

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

        if self.use_subcomms:
            self._e2n_half_size_subcomms(vecs_ger, vecs_ung, molecule, basis,
                                         scf_tensors, eri_dict, dft_dict,
                                         pe_dict, profiler)
        else:
            self._e2n_half_size_single_comm(vecs_ger, vecs_ung, molecule, basis,
                                            scf_tensors, eri_dict, dft_dict,
                                            pe_dict, profiler)

    def _e2n_half_size_subcomms(self,
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

        n_general = min(n_ger, n_ung)
        n_extra_ger = n_ger - n_general
        n_extra_ung = n_ung - n_general

        n_total = n_general + n_extra_ger + n_extra_ung

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

        # determine subcomm size

        if self.rank == mpi_master():
            dt_and_subcomm_size = []

            for subcomm_size in range(1, self.nodes + 1):
                if self.nodes % subcomm_size != 0:
                    continue

                n_subcomms = self.nodes // subcomm_size

                ave, res = divmod(n_total, n_subcomms)
                counts = [
                    ave + 1 if p < res else ave for p in range(n_subcomms)
                ]

                time_per_fock = self.serial_ratio + (
                    1 - self.serial_ratio) / subcomm_size
                dt = max(counts) * time_per_fock

                dt_and_subcomm_size.append((dt, subcomm_size))

            # find subcomm_size with smallest dt
            dt_and_subcomm_size.sort()
            subcomm_size = dt_and_subcomm_size[0][1]
        else:
            subcomm_size = None
        subcomm_size = self.comm.bcast(subcomm_size, root=mpi_master())

        # create subcomms

        grps = [p // subcomm_size for p in range(self.nodes)]
        subcomm = SubCommunicators(self.comm, grps)
        local_comm = subcomm.local_comm
        cross_comm = subcomm.cross_comm

        is_local_master = (local_comm.Get_rank() == mpi_master())

        # make data available on local master processes
        # and determine subcomm indices

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
        local_master_ranks = local_comm.bcast(local_master_ranks,
                                              root=mpi_master())

        # batch_size corresponds to the number of subcomms
        batch_size = self.nodes // subcomm_size

        num_batches = n_total // batch_size
        if n_total % batch_size != 0:
            num_batches += 1

        # go through batches

        if self.rank == mpi_master():
            batch_str = f'Processing {n_total} Fock builds'
            if batch_size > 1:
                batch_str += f' on {batch_size} subcommunicators'
            batch_str += '...'
            self.ostream.print_info(batch_str)
            self.ostream.flush()

        if self._debug:
            self.ostream.print_info(
                '==DEBUG== batch_size: {}'.format(batch_size))
            self.ostream.print_blank()
            self.ostream.flush()

        for batch_ind in range(num_batches):

            if self._debug:
                self.ostream.print_info('==DEBUG== batch {}/{}'.format(
                    batch_ind + 1, num_batches))
                self.ostream.flush()

            # form density matrices

            batch_start = batch_size * batch_ind
            batch_end = min(batch_start + batch_size, n_total)

            if is_local_master:
                dks = []
                kns = []

            vec_list = [None for idx in range(len(local_master_ranks))]

            e2_ger, e2_ung = None, None
            fock_ger, fock_ung = None, None

            for idx, local_master_rank in enumerate(local_master_ranks):

                if idx + batch_start >= batch_end:
                    break

                col = idx + batch_start

                # determine vec_type

                if col < n_general:
                    vec_type = 'general'
                elif n_extra_ger > 0:
                    vec_type = 'gerade'
                elif n_extra_ung > 0:
                    vec_type = 'ungerade'

                # form full-size vec

                if vec_type == 'general':
                    v_ger = vecs_ger.get_full_vector(col,
                                                     root=local_master_rank)
                    v_ung = vecs_ung.get_full_vector(col,
                                                     root=local_master_rank)
                    if self.rank == local_master_rank:
                        # full-size gerade trial vector
                        vec_list[idx] = np.hstack((v_ger, v_ger))
                        vec_list[idx] += np.hstack((v_ung, -v_ung))

                elif vec_type == 'gerade':
                    v_ger = vecs_ger.get_full_vector(col,
                                                     root=local_master_rank)
                    if self.rank == local_master_rank:
                        # full-size gerade trial vector
                        vec_list[idx] = np.hstack((v_ger, v_ger))

                elif vec_type == 'ungerade':
                    v_ung = vecs_ung.get_full_vector(col,
                                                     root=local_master_rank)
                    if self.rank == local_master_rank:
                        # full-size ungerade trial vector
                        vec_list[idx] = np.hstack((v_ung, -v_ung))

            self.comm.barrier()

            if subcomm_index + batch_start < batch_end:

                local_master_rank = local_master_ranks[subcomm_index]

                if is_local_master:
                    vec = vec_list[subcomm_index]

                    half_size = vec.shape[0] // 2

                    # build density

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
                else:
                    dks = None

                # form Fock matrices

                self._print_mem_debug_info('before Fock build')

                fock = self._comp_lr_fock(dks, molecule, basis, eri_dict,
                                          dft_dict, pe_dict, profiler,
                                          local_comm)

                self._print_mem_debug_info('after  Fock build')

                if is_local_master:
                    raw_fock_ger = []
                    raw_fock_ung = []

                    raw_kns_ger = []
                    raw_kns_ung = []

                    col = subcomm_index + batch_start

                    if col < n_general:
                        raw_fock_ger.append(0.5 * (fock[0] - fock[0].T))
                        raw_fock_ung.append(0.5 * (fock[0] + fock[0].T))
                        raw_kns_ger.append(0.5 * (kns[0] + kns[0].T))
                        raw_kns_ung.append(0.5 * (kns[0] - kns[0].T))

                    elif n_extra_ger > 0:
                        raw_fock_ger.append(0.5 * (fock[0] - fock[0].T))
                        raw_kns_ger.append(0.5 * (kns[0] + kns[0].T))

                    elif n_extra_ung > 0:
                        raw_fock_ung.append(0.5 * (fock[0] + fock[0].T))
                        raw_kns_ung.append(0.5 * (kns[0] - kns[0].T))

                    fock = raw_fock_ger + raw_fock_ung
                    kns = raw_kns_ger + raw_kns_ung

                    batch_ger = len(raw_fock_ger)
                    batch_ung = len(raw_fock_ung)

                if is_local_master:

                    e2_ger = np.zeros((half_size, batch_ger))
                    e2_ung = np.zeros((half_size, batch_ung))

                    if self.nonlinear:
                        fock_ger = np.zeros((norb**2, batch_ger))
                        fock_ung = np.zeros((norb**2, batch_ung))

                    for ifock in range(batch_ger + batch_ung):
                        fak = fock[ifock]

                        if getattr(self, 'core_excitation', False):
                            core_exc_orb_inds = list(
                                range(self.num_core_orbitals)) + list(
                                    range(nocc, norb))
                            mo_core_exc = mo[:, core_exc_orb_inds]
                            fak_mo = np.linalg.multi_dot(
                                [mo_core_exc.T, fak, mo_core_exc])
                        else:
                            fak_mo = np.linalg.multi_dot([mo.T, fak, mo])

                        kfa_mo = self.commut(fa_mo.T, kns[ifock])

                        fat_mo = fak_mo + kfa_mo

                        if getattr(self, 'core_excitation', False):
                            gmo = -self.commut_mo_density(
                                fat_mo, nocc, self.num_core_orbitals)
                            gmo_vec_halfsize = self.lrmat2vec(
                                gmo, nocc, norb,
                                self.num_core_orbitals)[:half_size]
                        else:
                            gmo = -self.commut_mo_density(fat_mo, nocc)
                            gmo_vec_halfsize = self.lrmat2vec(gmo, nocc,
                                                              norb)[:half_size]

                        # Note: fak_mo_vec uses full MO coefficients matrix since
                        # it is only for fock_ger/fock_ung in nonlinear response
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

                vecs_e2_ger = DistributedArray(e2_ger,
                                               self.comm,
                                               root=local_master_rank)
                vecs_e2_ung = DistributedArray(e2_ung,
                                               self.comm,
                                               root=local_master_rank)
                self._append_sigma_vectors(vecs_e2_ger, vecs_e2_ung)

                if self.nonlinear:
                    dist_fock_ger = DistributedArray(fock_ger,
                                                     self.comm,
                                                     root=local_master_rank)
                    dist_fock_ung = DistributedArray(fock_ung,
                                                     self.comm,
                                                     root=local_master_rank)
                    self._append_fock_matrices(dist_fock_ger, dist_fock_ung)

        self._append_trial_vectors(vecs_ger, vecs_ung)

        self.ostream.print_blank()

    def _e2n_half_size_single_comm(self,
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

        n_general = min(n_ger, n_ung)
        n_extra_ger = n_ger - n_general
        n_extra_ung = n_ung - n_general

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

        # determine number of batches

        n_ao = mo.shape[0] if self.rank == mpi_master() else None

        n_total = n_general + n_extra_ger + n_extra_ung

        batch_size = get_batch_size(self.batch_size, n_total, n_ao, self.comm)
        num_batches = get_number_of_batches(n_total, batch_size, self.comm)

        # go through batches

        if self.rank == mpi_master():
            batch_str = f'Processing {n_total} Fock build'
            if n_total > 1:
                batch_str += 's'
            batch_str += '...'
            self.ostream.print_info(batch_str)
            self.ostream.flush()

        if self._debug:
            self.ostream.print_info(
                '==DEBUG== batch_size: {}'.format(batch_size))
            self.ostream.print_blank()
            self.ostream.flush()

        for batch_ind in range(num_batches):

            if self._debug:
                self.ostream.print_info('==DEBUG== batch {}/{}'.format(
                    batch_ind + 1, num_batches))
                self.ostream.flush()

            # form density matrices

            batch_start = batch_size * batch_ind
            batch_end = min(batch_start + batch_size, n_total)

            if self.rank == mpi_master():
                dks = []
                kns = []
            else:
                dks = None

            for col in range(batch_start, batch_end):

                # determine vec_type

                if col < n_general:
                    vec_type = 'general'
                elif n_extra_ger > 0:
                    vec_type = 'gerade'
                elif n_extra_ung > 0:
                    vec_type = 'ungerade'

                # form full-size vec

                if vec_type == 'general':
                    v_ger = vecs_ger.get_full_vector(col)
                    v_ung = vecs_ung.get_full_vector(col)
                    if self.rank == mpi_master():
                        # full-size trial vector
                        vec = np.hstack((v_ger, v_ger))
                        vec += np.hstack((v_ung, -v_ung))

                elif vec_type == 'gerade':
                    v_ger = vecs_ger.get_full_vector(col)
                    if self.rank == mpi_master():
                        # full-size gerade trial vector
                        vec = np.hstack((v_ger, v_ger))

                elif vec_type == 'ungerade':
                    v_ung = vecs_ung.get_full_vector(col)
                    if self.rank == mpi_master():
                        # full-size ungerade trial vector
                        vec = np.hstack((v_ung, -v_ung))

                # build density

                if self.rank == mpi_master():
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

            # form Fock matrices

            self._print_mem_debug_info('before Fock build')

            fock = self._comp_lr_fock(dks, molecule, basis, eri_dict, dft_dict,
                                      pe_dict, profiler)

            self._print_mem_debug_info('after  Fock build')

            if self.rank == mpi_master():
                raw_fock_ger = []
                raw_fock_ung = []

                raw_kns_ger = []
                raw_kns_ung = []

                for col in range(batch_start, batch_end):
                    ifock = col - batch_start

                    if col < n_general:
                        raw_fock_ger.append(0.5 * (fock[ifock] - fock[ifock].T))
                        raw_fock_ung.append(0.5 * (fock[ifock] + fock[ifock].T))
                        raw_kns_ger.append(0.5 * (kns[ifock] + kns[ifock].T))
                        raw_kns_ung.append(0.5 * (kns[ifock] - kns[ifock].T))

                    elif n_extra_ger > 0:
                        raw_fock_ger.append(0.5 * (fock[ifock] - fock[ifock].T))
                        raw_kns_ger.append(0.5 * (kns[ifock] + kns[ifock].T))

                    elif n_extra_ung > 0:
                        raw_fock_ung.append(0.5 * (fock[ifock] + fock[ifock].T))
                        raw_kns_ung.append(0.5 * (kns[ifock] - kns[ifock].T))

                fock = raw_fock_ger + raw_fock_ung
                kns = raw_kns_ger + raw_kns_ung

                batch_ger = len(raw_fock_ger)
                batch_ung = len(raw_fock_ung)

            e2_ger = None
            e2_ung = None

            fock_ger = None
            fock_ung = None

            if self.rank == mpi_master():

                e2_ger = np.zeros((half_size, batch_ger))
                e2_ung = np.zeros((half_size, batch_ung))

                if self.nonlinear:
                    fock_ger = np.zeros((norb**2, batch_ger))
                    fock_ung = np.zeros((norb**2, batch_ung))

                for ifock in range(batch_ger + batch_ung):
                    fak = fock[ifock]

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

                    # Note: fak_mo_vec uses full MO coefficients matrix since
                    # it is only for fock_ger/fock_ung in nonlinear response
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

            vecs_e2_ger = DistributedArray(e2_ger, self.comm)
            vecs_e2_ung = DistributedArray(e2_ung, self.comm)
            self._append_sigma_vectors(vecs_e2_ger, vecs_e2_ung)

            if self.nonlinear:
                dist_fock_ger = DistributedArray(fock_ger, self.comm)
                dist_fock_ung = DistributedArray(fock_ung, self.comm)
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
                      profiler=None,
                      comm=None):
        """
        Computes Fock/Fxc matrix (2e part) for linear response calculation.

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

        :return:
            The Fock matrix (2e part).
        """

        if comm is None:
            comm = self.comm
            redistribute_xc_molgrid = False
        else:
            # use subcommunicators
            redistribute_xc_molgrid = True

        comm_rank = comm.Get_rank()

        screening = eri_dict['screening']

        molgrid = dft_dict['molgrid']
        gs_density = dft_dict['gs_density']

        if comm_rank == mpi_master():
            num_densities = len(dens)
        else:
            num_densities = None
        num_densities = comm.bcast(num_densities, root=mpi_master())

        if comm_rank != mpi_master():
            dens = [None for idx in range(num_densities)]

        for idx in range(num_densities):
            dens[idx] = comm.bcast(dens[idx], root=mpi_master())

        thresh_int = int(-math.log10(self.eri_thresh))

        t0 = tm.time()

        fock_drv = FockDriver(comm)
        fock_drv._set_block_size_factor(self._block_size_factor)

        # determine fock_type and exchange_scaling_factor
        fock_type = '2jk'
        exchange_scaling_factor = 1.0
        if self._dft:
            if self.xcfun.is_hybrid():
                fock_type = '2jkx'
                exchange_scaling_factor = self.xcfun.get_frac_exact_exchange()
            else:
                fock_type = 'j'
                exchange_scaling_factor = 0.0

        # further determine exchange_scaling_factor, erf_k_coef and omega
        need_omega = (self._dft and self.xcfun.is_range_separated())
        if need_omega:
            exchange_scaling_factor = (self.xcfun.get_rs_alpha() +
                                       self.xcfun.get_rs_beta())
            erf_k_coef = -self.xcfun.get_rs_beta()
            omega = self.xcfun.get_rs_omega()
        else:
            erf_k_coef, omega = None, None

        fock_arrays = []

        for idx in range(num_densities):
            den_mat_for_fock = make_matrix(basis, mat_t.general)
            den_mat_for_fock.set_values(dens[idx])

            self._print_mem_debug_info('before restgen Fock build')
            fock_mat = fock_drv.compute(screening, den_mat_for_fock, fock_type,
                                        exchange_scaling_factor, 0.0,
                                        thresh_int)
            self._print_mem_debug_info('after  restgen Fock build')

            fock_np = fock_mat.to_numpy()
            fock_mat = Matrix()

            if fock_type == 'j':
                # for pure functional
                fock_np *= 2.0

            if need_omega:
                # for range-separated functional
                self._print_mem_debug_info('before restgen erf Fock build')
                fock_mat = fock_drv.compute(screening, den_mat_for_fock,
                                            'kx_rs', erf_k_coef, omega,
                                            thresh_int)
                self._print_mem_debug_info('after  restgen erf Fock build')

                fock_np -= fock_mat.to_numpy()
                fock_mat = Matrix()

            fock_arrays.append(fock_np)

        if profiler is not None:
            profiler.add_timing_info('FockERI', tm.time() - t0)

        if self._dft:
            t0 = tm.time()

            if redistribute_xc_molgrid:
                molgrid.re_distribute_counts_and_displacements(
                    comm.Get_rank(), comm.Get_size())

            xc_drv = XCIntegrator()
            xc_drv.integrate_fxc_fock(fock_arrays, molecule, basis, dens,
                                      gs_density, molgrid, self.xcfun)

            if profiler is not None:
                profiler.add_timing_info('FockXC', tm.time() - t0)

        if self._pe:
            t0 = tm.time()
            for idx in range(num_densities):
                # Note: only closed shell density for now
                dm = dens[idx] * 2.0
                V_emb = self._embedding_drv.compute_pe_contributions(
                    density_matrix=dm)
                if comm_rank == mpi_master():
                    fock_arrays[idx] += V_emb

            if profiler is not None:
                profiler.add_timing_info('FockPE', tm.time() - t0)

        for idx in range(len(fock_arrays)):
            fock_arrays[idx] = comm.reduce(fock_arrays[idx], root=mpi_master())

        if comm_rank == mpi_master():
            return fock_arrays
        else:
            return None

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

    def get_prop_grad(self, operator, components, molecule, basis, scf_tensors):
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
            if self.rank == mpi_master():
                dipole_mats = compute_electric_dipole_integrals(
                    molecule, basis, [0.0, 0.0, 0.0])
                integrals = tuple(dipole_mats)
            else:
                integrals = tuple()

        elif operator in ['linear_momentum', 'linear momentum']:
            if self.rank == mpi_master():
                linmom_mats = compute_linear_momentum_integrals(molecule, basis)
                integrals = (
                    -1.0 * linmom_mats[0],
                    -1.0 * linmom_mats[1],
                    -1.0 * linmom_mats[2],
                )
            else:
                integrals = tuple()

        elif operator in ['angular_momentum', 'angular momentum']:
            if self.rank == mpi_master():
                angmom_mats = compute_angular_momentum_integrals(
                    molecule, basis, [0.0, 0.0, 0.0])
                integrals = (
                    -1.0 * angmom_mats[0],
                    -1.0 * angmom_mats[1],
                    -1.0 * angmom_mats[2],
                )
            else:
                integrals = tuple()

        elif operator in ['magnetic_dipole', 'magnetic dipole']:
            if self.rank == mpi_master():
                angmom_mats = compute_angular_momentum_integrals(
                    molecule, basis, [0.0, 0.0, 0.0])
                integrals = (
                    0.5 * angmom_mats[0],
                    0.5 * angmom_mats[1],
                    0.5 * angmom_mats[2],
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
                        np.linalg.multi_dot([mo.T, P.T, mo]), nocc)
                    for P in integral_comps
                ]
                gradients = tuple(
                    self.lrmat2vec(m, nocc, norb) for m in matrices)

            return gradients

        else:
            return tuple()

    def get_complex_prop_grad(self, operator, components, molecule, basis,
                              scf_tensors):
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
            if self.rank == mpi_master():
                dipole_mats = compute_electric_dipole_integrals(
                    molecule, basis, [0.0, 0.0, 0.0])
                integrals = (
                    dipole_mats[0] + 0j,
                    dipole_mats[1] + 0j,
                    dipole_mats[2] + 0j,
                )
            else:
                integrals = tuple()

        elif operator in ['linear_momentum', 'linear momentum']:
            if self.rank == mpi_master():
                linmom_mats = compute_linear_momentum_integrals(molecule, basis)
                integrals = (
                    -1j * linmom_mats[0],
                    -1j * linmom_mats[1],
                    -1j * linmom_mats[2],
                )
            else:
                integrals = tuple()

        elif operator in ['angular_momentum', 'angular momentum']:
            if self.rank == mpi_master():
                angmom_mats = compute_angular_momentum_integrals(
                    molecule, basis, [0.0, 0.0, 0.0])
                integrals = (
                    -1j * angmom_mats[0],
                    -1j * angmom_mats[1],
                    -1j * angmom_mats[2],
                )
            else:
                integrals = tuple()

        elif operator in ['magnetic_dipole', 'magnetic dipole']:
            if self.rank == mpi_master():
                angmom_mats = compute_angular_momentum_integrals(
                    molecule, basis, [0.0, 0.0, 0.0])
                integrals = (
                    0.5j * angmom_mats[0],
                    0.5j * angmom_mats[1],
                    0.5j * angmom_mats[2],
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
                factor * (-1.0) * self.commut_mo_density(
                    np.linalg.multi_dot([mo.T, P.T, mo]), nocc)
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

            nto_coefs = nto_mo.alpha_to_numpy()

            # hole
            ind_occ = nocc - i_nto - 1
            vis_drv.compute(cubic_grid, molecule, basis, nto_coefs, ind_occ,
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
            vis_drv.compute(cubic_grid, molecule, basis, nto_coefs, ind_vir,
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
                homo_str = f'core_{i + 1}'
            else:
                homo_str = 'HOMO' if i == nocc - 1 else f'HOMO-{nocc - 1 - i}'

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
