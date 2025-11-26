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

from datetime import datetime
import numpy as np
import time as tm
import math
import sys

from .veloxchemlib import T4CScreener
from .veloxchemlib import MolecularGrid, XCIntegrator
from .veloxchemlib import (mpi_master, rotatory_strength_in_cgs, hartree_in_ev,
                           hartree_in_inverse_nm, fine_structure_constant,
                           extinction_coefficient_from_beta, avogadro_constant,
                           bohr_in_angstrom, hartree_in_wavenumber)
from .veloxchemlib import make_matrix, mat_t
from .matrix import Matrix
from .distributedarray import DistributedArray
from .subcommunicators import SubCommunicators
from .rifockdriver import RIFockDriver
from .fockdriver import FockDriver
from .griddriver import GridDriver
from .molecularorbitals import MolecularOrbitals, molorb
from .visualizationdriver import VisualizationDriver
from .oneeints import (compute_electric_dipole_integrals,
                       compute_quadrupole_integrals,
                       compute_linear_momentum_integrals,
                       compute_angular_momentum_integrals)
from .sanitychecks import (dft_sanity_check, pe_sanity_check,
                           solvation_model_sanity_check)
from .errorhandler import assert_msg_critical
from .inputparser import (parse_input, print_keywords, print_attributes,
                          get_random_string_parallel)
from .dftutils import get_default_grid_level, print_xc_reference
from .checkpoint import write_rsp_hdf5
from .batchsize import get_batch_size
from .batchsize import get_number_of_batches
from .cpcmdriver import CpcmDriver

try:
    import matplotlib.pyplot as plt
    import matplotlib.lines as mlines
except ImportError:
    pass


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

        # RI-J
        self.ri_coulomb = False
        self.ri_auxiliary_basis = 'def2-universal-jfit'
        self._ri_drv = None

        # dft
        self.xcfun = None
        self.grid_level = None
        self._dft = False

        # polarizable embedding
        self.potfile = None
        self.pe_options = {}
        self._pe = False
        self.embedding = None
        self._embedding_drv = None

        # static electric field
        self.electric_field = None

        # point charges
        self.point_charges = None

        # solvation model
        self.solvation_model = None
        self.non_equilibrium_solv = True

        # C-PCM setup
        self._cpcm = False
        self.cpcm_drv = None
        self.cpcm_epsilon = 78.39
        self.cpcm_optical_epsilon = 1.777849
        self.cpcm_grid_per_sphere = (194, 110)
        self.cpcm_cg_thresh = 1.0e-8
        self.cpcm_x = 0
        self.cpcm_custom_vdw_radii = None

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
        self.print_level = 1

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

        # group label used to save the response results in a checkpoint file
        self.group_label = 'rsp'

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
                'non_equilibrium_solv':
                    ('bool',
                     'toggle use of non-equilibrium solvation for response'),
            },
            'method_settings': {
                'ri_coulomb': ('bool', 'use RI-J approximation'),
                'ri_auxiliary_basis': ('str', 'RI-J auxiliary basis set'),
                'xcfun': ('str_upper', 'exchange-correlation functional'),
                'grid_level': ('int', 'accuracy level of DFT grid'),
                'potfile': ('str', 'potential file for polarizable embedding'),
                'electric_field': ('seq_fixed', 'static electric field'),
                'solvation_model': ('str', 'solvation model'),
                'cpcm_grid_per_sphere':
                    ('seq_fixed_int', 'number of C-PCM grid points per sphere'),
                'cpcm_cg_thresh':
                    ('float', 'threshold for solving C-PCM charges'),
                'cpcm_epsilon':
                    ('float', 'dielectric constant of solvent (C-PCM)'),
                'cpcm_optical_epsilon':
                    ('float', 'optical dielectric constant of solvent (C-PCM)'),
                'cpcm_x': ('float', 'parameter for scaling function (C-PCM)'),
                'cpcm_custom_vdw_radii':
                    ('seq_fixed_str', 'custom vdw radii for C-PCM'),
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
                self.checkpoint_file = f'{self.filename}_rsp.h5'

        method_keywords = {
            key: val[0]
            for key, val in self._input_keywords['method_settings'].items()
        }

        parse_input(self, method_keywords, method_dict)

        dft_sanity_check(self, 'update_settings')

        pe_sanity_check(self, method_dict)

        solvation_model_sanity_check(self)

        assert_msg_critical(not self._smd,
                            'Cannot use SMD in response calculation')

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

        if self.rank == mpi_master():
            screening = T4CScreener()
            screening.partition(basis, molecule, 'eri')
        else:
            screening = None
        screening = self.comm.bcast(screening, root=mpi_master())

        if self.ri_coulomb:
            self._ri_drv = RIFockDriver(self.comm, self.ostream)
            self._ri_drv.prepare_buffers(molecule,
                                         basis,
                                         self.ri_auxiliary_basis,
                                         verbose=True)

        return {
            'screening': screening,
        }

    def _init_dft(self, molecule, scf_tensors, silent=False):
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
            if not silent:
                print_xc_reference(self.xcfun, self.ostream)

            grid_drv = GridDriver(self.comm)
            grid_level = (get_default_grid_level(self.xcfun)
                          if self.grid_level is None else self.grid_level)
            grid_drv.set_level(grid_level)

            grid_t0 = tm.time()
            molgrid = grid_drv.generate(molecule, self._xcfun_ldstaging)
            n_grid_points = molgrid.number_of_points()
            if (not silent) and (self.print_level > 1):
                self.ostream.print_info(
                    'Molecular grid with {0:d} points generated in {1:.2f} sec.'
                    .format(n_grid_points,
                            tm.time() - grid_t0))
                self.ostream.print_blank()

            if self.rank == mpi_master():
                # Note: make gs_density a tuple
                if scf_tensors['scf_type'] == 'restricted':
                    gs_density = (scf_tensors['D_alpha'].copy(),)
                else:
                    gs_density = (scf_tensors['D_alpha'].copy(),
                                  scf_tensors['D_beta'].copy())
            else:
                gs_density = None
            # TODO: bcast D_alpha and D_beta separately
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

    def _init_pe(self, molecule, basis, silent=False):
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
                self.embedding['settings']['embedding_method'] == 'PE',
                'PolarizableEmbedding: Invalid embedding_method. Only PE is supported.'
            )

            from .embedding import PolarizableEmbeddingLRS

            self._embedding_drv = PolarizableEmbeddingLRS(
                molecule=molecule,
                ao_basis=basis,
                options=self.embedding,
                comm=self.comm)

            # TODO: print PyFraME info

            if not silent:
                pot_info = 'Reading polarizable embedding: {}'.format(
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

    def _init_cpcm(self, molecule):
        """
        Initializes C-PCM.

        :param molecule:
            The molecule.
        """

        # C-PCM setup
        if self._cpcm:

            self.cpcm_drv = CpcmDriver(self.comm, self.ostream)
            self.cpcm_drv.print_cpcm_info()

            self.cpcm_drv.grid_per_sphere = self.cpcm_grid_per_sphere
            self.cpcm_drv.epsilon = self.cpcm_epsilon
            self.cpcm_drv.optical_epsilon = self.cpcm_optical_epsilon
            self.cpcm_drv.x = self.cpcm_x
            self.cpcm_drv.custom_vdw_radii = self.cpcm_custom_vdw_radii

            cpcm_grid_t0 = tm.time()

            self.cpcm_drv.init(molecule, do_nuclear=False)

            if self.print_level > 1:
                self.ostream.print_info(
                    f'C-PCM grid with {self.cpcm_drv._cpcm_grid.shape[0]} points '
                    + f'generated in {tm.time() - cpcm_grid_t0:.2f} sec.')
                self.ostream.print_blank()
                self.ostream.flush()

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
                       profiler=None,
                       method_type='restricted'):

        if self.use_subcomms and self.ri_coulomb:
            self.use_subcomms = False
            warn_msg = 'Use of subcomms is disabled for RI-J.'
            self.ostream.print_warning(warn_msg)
            self.ostream.print_blank()
            self.ostream.flush()

        # TODO: enable subcomms for RI-J

        if self.use_subcomms:
            if method_type == 'restricted':
                self._e2n_half_size_subcomms(vecs_ger, vecs_ung, molecule,
                                             basis, scf_tensors, eri_dict,
                                             dft_dict, pe_dict, profiler)
            else:
                # TODO: enable subcomms for unrestricted
                assert_msg_critical(
                    False, 'LinearSolver._e2n_half_size: '
                    'Cannot use subcomms for unrestricted case')
        else:
            if method_type == 'restricted':
                self._e2n_half_size_single_comm(vecs_ger, vecs_ung, molecule,
                                                basis, scf_tensors, eri_dict,
                                                dft_dict, pe_dict, profiler)
            else:
                self._e2n_half_size_single_comm_unrestricted(
                    vecs_ger, vecs_ung, molecule, basis, scf_tensors, eri_dict,
                    dft_dict, pe_dict, profiler)

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

        prep_t0 = tm.time()

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

                # make sure that number of subcomms does not exceed number
                # of trial vectors
                if n_subcomms > n_total:
                    continue

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

        if profiler is not None:
            profiler.add_timing_info('PreProc', tm.time() - prep_t0)

        vecs_e2_ger_data = None
        vecs_e2_ung_data = None

        # go through batches

        if self.rank == mpi_master() and self.print_level > 1:
            batch_str = f'Processing {n_total} Fock build'
            if n_total > 1:
                batch_str += 's'
            if batch_size > 1:
                batch_str += f' on {batch_size} subcommunicators'
            batch_str += '...'
            self.ostream.print_info(batch_str)
            self.ostream.print_blank()
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

            # TODO: use Alltoallv

            prep_t0 = tm.time()

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

            if profiler is not None:
                profiler.add_timing_info('PreProc', tm.time() - prep_t0)

            if subcomm_index + batch_start < batch_end:

                prep_t0 = tm.time()

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

                if profiler is not None:
                    profiler.add_timing_info('PreProc', tm.time() - prep_t0)

                # form Fock matrices

                fock = self._comp_lr_fock(dks, molecule, basis, eri_dict,
                                          dft_dict, pe_dict, profiler,
                                          local_comm)

                if profiler is not None:
                    # only increment FockCount on local master
                    if is_local_master:
                        profiler.add_timing_info('_FockCount_', 1)

                prep_t0 = tm.time()

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

                if profiler is not None:
                    profiler.add_timing_info('PostProc', tm.time() - prep_t0)

            self.comm.barrier()

            prep_t0 = tm.time()

            # TODO: use Alltoallv

            local_e2_ger_data = None
            local_e2_ung_data = None

            for idx, local_master_rank in enumerate(local_master_ranks):

                if idx + batch_start >= batch_end:
                    break

                dist_e2_ger = DistributedArray(e2_ger,
                                               self.comm,
                                               root=local_master_rank)
                dist_e2_ung = DistributedArray(e2_ung,
                                               self.comm,
                                               root=local_master_rank)

                # Note: accumulate per-rank data to avoid incremental appending

                if local_e2_ger_data is None:
                    local_e2_ger_data = dist_e2_ger.data.copy()
                else:
                    local_e2_ger_data = np.hstack(
                        (local_e2_ger_data, dist_e2_ger.data))

                if local_e2_ung_data is None:
                    local_e2_ung_data = dist_e2_ung.data.copy()
                else:
                    local_e2_ung_data = np.hstack(
                        (local_e2_ung_data, dist_e2_ung.data))

                if self.nonlinear:
                    dist_fock_ger = DistributedArray(fock_ger,
                                                     self.comm,
                                                     root=local_master_rank)
                    dist_fock_ung = DistributedArray(fock_ung,
                                                     self.comm,
                                                     root=local_master_rank)
                    # TODO: avoid incremental append
                    self._append_fock_matrices(dist_fock_ger, dist_fock_ung)

            # Note: accumulate per-batch data to avoid incremental appending

            if vecs_e2_ger_data is None:
                vecs_e2_ger_data = local_e2_ger_data.copy()
            else:
                vecs_e2_ger_data = np.hstack(
                    (vecs_e2_ger_data, local_e2_ger_data))

            if vecs_e2_ung_data is None:
                vecs_e2_ung_data = local_e2_ung_data.copy()
            else:
                vecs_e2_ung_data = np.hstack(
                    (vecs_e2_ung_data, local_e2_ung_data))

            if profiler is not None:
                profiler.add_timing_info('PostProc', tm.time() - prep_t0)

        prep_t0 = tm.time()

        # Note: append sigma and trial vectors only once

        vecs_e2_ger = DistributedArray(vecs_e2_ger_data,
                                       self.comm,
                                       distribute=False)

        vecs_e2_ung = DistributedArray(vecs_e2_ung_data,
                                       self.comm,
                                       distribute=False)

        self._append_sigma_vectors(vecs_e2_ger, vecs_e2_ung)

        self._append_trial_vectors(vecs_ger, vecs_ung)

        if profiler is not None:
            profiler.add_timing_info('AppendVec', tm.time() - prep_t0)

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

        prep_t0 = tm.time()

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

        if profiler is not None:
            profiler.add_timing_info('PreProc', tm.time() - prep_t0)

        vecs_e2_ger_data = None
        vecs_e2_ung_data = None

        # go through batches

        if self.rank == mpi_master() and self.print_level > 1:
            batch_str = f'Processing {n_total} Fock build'
            if n_total > 1:
                batch_str += 's'
            batch_str += '...'
            self.ostream.print_info(batch_str)
            self.ostream.print_blank()
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

            prep_t0 = tm.time()

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

            if profiler is not None:
                profiler.add_timing_info('PreProc', tm.time() - prep_t0)

            # form Fock matrices

            fock = self._comp_lr_fock(dks, molecule, basis, eri_dict, dft_dict,
                                      pe_dict, profiler)

            if profiler is not None:
                # only increment FockCount on master rank
                if self.rank == mpi_master():
                    profiler.add_timing_info('_FockCount_', len(dks))

            prep_t0 = tm.time()

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

            dist_e2_ger = DistributedArray(e2_ger, self.comm)
            dist_e2_ung = DistributedArray(e2_ung, self.comm)

            if vecs_e2_ger_data is None:
                vecs_e2_ger_data = dist_e2_ger.data.copy()
            else:
                vecs_e2_ger_data = np.hstack(
                    (vecs_e2_ger_data, dist_e2_ger.data))

            if vecs_e2_ung_data is None:
                vecs_e2_ung_data = dist_e2_ung.data.copy()
            else:
                vecs_e2_ung_data = np.hstack(
                    (vecs_e2_ung_data, dist_e2_ung.data))

            if self.nonlinear:
                dist_fock_ger = DistributedArray(fock_ger, self.comm)
                dist_fock_ung = DistributedArray(fock_ung, self.comm)
                # TODO: avoid incremental append
                self._append_fock_matrices(dist_fock_ger, dist_fock_ung)

            if profiler is not None:
                profiler.add_timing_info('PostProc', tm.time() - prep_t0)

        prep_t0 = tm.time()

        # Note: append sigma and trial vectors only once

        vecs_e2_ger = DistributedArray(vecs_e2_ger_data,
                                       self.comm,
                                       distribute=False)

        vecs_e2_ung = DistributedArray(vecs_e2_ung_data,
                                       self.comm,
                                       distribute=False)

        self._append_sigma_vectors(vecs_e2_ger, vecs_e2_ung)

        self._append_trial_vectors(vecs_ger, vecs_ung)

        if profiler is not None:
            profiler.add_timing_info('AppendVec', tm.time() - prep_t0)

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
            if self.ri_coulomb:
                assert_msg_critical(
                    fock_type == 'j',
                    'LinearSolver: RI is only applicable to pure DFT functional'
                )

            if self.ri_coulomb and fock_type == 'j':
                # symmetrize density for RI-J
                den_mat_for_ri_j = make_matrix(basis, mat_t.symmetric)
                den_mat_for_ri_j.set_values(0.5 * (dens[idx] + dens[idx].T))

                fock_mat = self._ri_drv.compute(den_mat_for_ri_j, 'j')
            else:
                den_mat_for_fock = make_matrix(basis, mat_t.general)
                den_mat_for_fock.set_values(dens[idx])

                fock_mat = fock_drv.compute(screening, den_mat_for_fock,
                                            fock_type, exchange_scaling_factor,
                                            0.0, thresh_int)

            fock_np = fock_mat.to_numpy()
            fock_mat = Matrix()

            if fock_type == 'j':
                # for pure functional
                fock_np *= 2.0

            if need_omega:
                # for range-separated functional
                fock_mat = fock_drv.compute(screening, den_mat_for_fock,
                                            'kx_rs', erf_k_coef, omega,
                                            thresh_int)

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

        if self._cpcm:

            t0 = tm.time()

            for idx in range(num_densities):
                Fock_sol = self.cpcm_drv.compute_response_fock(
                    molecule, basis, dens[idx] * 2.0, self.cpcm_cg_thresh,
                    self.non_equilibrium_solv)

                if comm_rank == mpi_master():
                    fock_arrays[idx] += Fock_sol

            if profiler is not None:
                profiler.add_timing_info('FockCPCM', tm.time() - t0)

        for idx in range(len(fock_arrays)):
            fock_arrays[idx] = comm.reduce(fock_arrays[idx], root=mpi_master())

        if comm_rank == mpi_master():
            return fock_arrays
        else:
            return None

    def _comp_lr_fock_unrestricted(self,
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

        dens_a, dens_b = dens

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
            num_densities = len(dens_a)
        else:
            num_densities = None
        num_densities = comm.bcast(num_densities, root=mpi_master())

        if comm_rank != mpi_master():
            dens_a = [None for idx in range(num_densities)]
            dens_b = [None for idx in range(num_densities)]

        for idx in range(num_densities):
            dens_a[idx] = comm.bcast(dens_a[idx], root=mpi_master())
            dens_b[idx] = comm.bcast(dens_b[idx], root=mpi_master())

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
            if self.ri_coulomb:
                assert_msg_critical(
                    fock_type == 'j',
                    'LinearSolver: RI is only applicable to pure DFT functional'
                )

            if self.ri_coulomb and fock_type == 'j':
                # TODO: double check
                # symmetrize density for RI-J
                den_mat_for_ri_j = make_matrix(basis, mat_t.symmetric)
                dens_Jab = dens_a[idx] + dens_b[idx]
                den_mat_for_ri_j.set_values(0.5 * (dens_Jab + dens_Jab.T))

                fock_mat = self._ri_drv.compute(den_mat_for_ri_j, 'j')
            else:
                # for now we calculate Ka, Kb and Jab separately for open-shell
                den_mat_for_Ka = make_matrix(basis, mat_t.general)
                den_mat_for_Ka.set_values(dens_a[idx])

                fock_mat = fock_drv.compute(screening, den_mat_for_Ka, 'kx',
                                            exchange_scaling_factor, 0.0,
                                            thresh_int)
                K_a_np = fock_mat.to_numpy()

                den_mat_for_Kb = make_matrix(basis, mat_t.general)
                den_mat_for_Kb.set_values(dens_b[idx])

                fock_mat = fock_drv.compute(screening, den_mat_for_Kb, 'kx',
                                            exchange_scaling_factor, 0.0,
                                            thresh_int)

                K_b_np = fock_mat.to_numpy()

                den_mat_for_Jab = make_matrix(basis, mat_t.general)
                den_mat_for_Jab.set_values(dens_a[idx] + dens_b[idx])

                fock_mat = fock_drv.compute(screening, den_mat_for_Jab, 'j',
                                            exchange_scaling_factor, 0.0,
                                            thresh_int)

                J_ab_np = fock_mat.to_numpy()

                fock_mat_a_np = J_ab_np - K_a_np
                fock_mat_b_np = J_ab_np - K_b_np

            if need_omega:
                # for range-separated functional
                fock_mat = fock_drv.compute(screening, den_mat_for_Ka, 'kx_rs',
                                            erf_k_coef, omega, thresh_int)

                fock_mat_a_np -= fock_mat.to_numpy()

                fock_mat = fock_drv.compute(screening, den_mat_for_Kb, 'kx_rs',
                                            erf_k_coef, omega, thresh_int)

                fock_mat_b_np -= fock_mat.to_numpy()

            den_mat_for_Ka = Matrix()
            den_mat_for_Kb = Matrix()
            den_mat_for_Jab = Matrix()
            fock_mat = Matrix()

            fock_arrays.append(fock_mat_a_np)
            fock_arrays.append(fock_mat_b_np)

        if profiler is not None:
            profiler.add_timing_info('FockERI', tm.time() - t0)

        if self._dft:
            # TODO: unrestricted
            t0 = tm.time()

            if redistribute_xc_molgrid:
                molgrid.re_distribute_counts_and_displacements(
                    comm.Get_rank(), comm.Get_size())

            dens_a_and_b = []
            for da, db in zip(dens_a, dens_b):
                dens_a_and_b.append(da)
                dens_a_and_b.append(db)

            xc_drv = XCIntegrator()
            xc_drv.integrate_fxc_fock(fock_arrays, molecule, basis,
                                      dens_a_and_b, gs_density, molgrid,
                                      self.xcfun)

            if profiler is not None:
                profiler.add_timing_info('FockXC', tm.time() - t0)

        if self._pe:
            t0 = tm.time()

            for idx in range(num_densities):
                dm = dens_a[idx] + dens_b[idx]
                V_emb = self._embedding_drv.compute_pe_contributions(
                    density_matrix=dm)

                if comm_rank == mpi_master():
                    fock_arrays[idx * 2 + 0] += V_emb
                    fock_arrays[idx * 2 + 1] += V_emb

            if profiler is not None:
                profiler.add_timing_info('FockPE', tm.time() - t0)

        if self._cpcm:
            t0 = tm.time()

            for idx in range(num_densities):
                Fock_sol = self.cpcm_drv.compute_response_fock(
                    molecule, basis, dens_a[idx] + dens_b[idx],
                    self.cpcm_cg_thresh, self.non_equilibrium_solv)

                if comm_rank == mpi_master():
                    fock_arrays[idx * 2 + 0] += Fock_sol
                    fock_arrays[idx * 2 + 1] += Fock_sol

            if profiler is not None:
                profiler.add_timing_info('FockCPCM', tm.time() - t0)

        for idx in range(num_densities):
            fock_arrays[idx * 2 + 0] = comm.reduce(fock_arrays[idx * 2 + 0],
                                                   root=mpi_master())
            fock_arrays[idx * 2 + 1] = comm.reduce(fock_arrays[idx * 2 + 1],
                                                   root=mpi_master())

        if comm_rank == mpi_master():
            return fock_arrays
        else:
            return None

    def _e2n_half_size_single_comm_unrestricted(self,
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

        prep_t0 = tm.time()

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

            mo_a = scf_tensors['C_alpha']
            mo_b = scf_tensors['C_beta']

            fa = scf_tensors['F_alpha']
            fb = scf_tensors['F_beta']

            nocc_a = molecule.number_of_alpha_electrons()
            nocc_b = molecule.number_of_beta_electrons()

            norb = mo_a.shape[1]

            if getattr(self, 'core_excitation', False):
                core_exc_orb_inds_a = list(range(
                    self.num_core_orbitals)) + list(range(nocc_a, norb))
                core_exc_orb_inds_b = list(range(
                    self.num_core_orbitals)) + list(range(nocc_b, norb))
                mo_core_exc_a = mo_a[:, core_exc_orb_inds_a]
                mo_core_exc_b = mo_b[:, core_exc_orb_inds_b]
                fa_mo = np.linalg.multi_dot(
                    [mo_core_exc_a.T, fa, mo_core_exc_a])
                fb_mo = np.linalg.multi_dot(
                    [mo_core_exc_b.T, fb, mo_core_exc_b])
            else:
                fa_mo = np.linalg.multi_dot([mo_a.T, fa, mo_a])
                fb_mo = np.linalg.multi_dot([mo_b.T, fb, mo_b])

        else:
            nocc_a = None
            norb = None

        nocc_a, norb = self.comm.bcast((nocc_a, norb), root=mpi_master())

        # determine number of batches

        n_ao = mo_a.shape[0] if self.rank == mpi_master() else None

        n_total = n_general + n_extra_ger + n_extra_ung

        batch_size = get_batch_size(self.batch_size, n_total, n_ao, self.comm)
        num_batches = get_number_of_batches(n_total, batch_size, self.comm)

        if profiler is not None:
            profiler.add_timing_info('PreProc', tm.time() - prep_t0)

        vecs_e2_ger_data = None
        vecs_e2_ung_data = None

        # go through batches

        if self.rank == mpi_master() and self.print_level > 1:
            batch_str = f'Processing {n_total} Fock build'
            if n_total > 1:
                batch_str += 's'
            batch_str += '...'
            self.ostream.print_info(batch_str)
            self.ostream.print_blank()
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
                dks_a = []
                kns_a = []
                dks_b = []
                kns_b = []
            else:
                dks_a = None
                kns_a = None
                dks_b = None
                kns_b = None

            prep_t0 = tm.time()

            for col in range(batch_start, batch_end):

                # determine vec_type

                if col < n_general:
                    vec_type = 'general'
                elif n_extra_ger > 0:
                    vec_type = 'gerade'
                elif n_extra_ung > 0:
                    vec_type = 'ungerade'

                # form full-size vec

                if getattr(self, 'core_excitation', False):
                    n_ov_a = self.num_core_orbitals * (norb - nocc_a)
                else:
                    n_ov_a = nocc_a * (norb - nocc_a)

                if vec_type == 'general':
                    v_ger = vecs_ger.get_full_vector(col)
                    v_ung = vecs_ung.get_full_vector(col)
                    if self.rank == mpi_master():
                        v_ger_a = v_ger[:n_ov_a]
                        v_ger_b = v_ger[n_ov_a:]
                        v_ung_a = v_ung[:n_ov_a]
                        v_ung_b = v_ung[n_ov_a:]
                        # full-size trial vector
                        vec_a = np.hstack((v_ger_a, v_ger_a))
                        vec_a += np.hstack((v_ung_a, -v_ung_a))
                        vec_b = np.hstack((v_ger_b, v_ger_b))
                        vec_b += np.hstack((v_ung_b, -v_ung_b))

                elif vec_type == 'gerade':
                    v_ger = vecs_ger.get_full_vector(col)
                    if self.rank == mpi_master():
                        v_ger_a = v_ger[:n_ov_a]
                        v_ger_b = v_ger[n_ov_a:]
                        # full-size gerade trial vector
                        vec_a = np.hstack((v_ger_a, v_ger_a))
                        vec_b = np.hstack((v_ger_b, v_ger_b))

                elif vec_type == 'ungerade':
                    v_ung = vecs_ung.get_full_vector(col)
                    if self.rank == mpi_master():
                        v_ung_a = v_ung[:n_ov_a]
                        v_ung_b = v_ung[n_ov_a:]
                        # full-size ungerade trial vector
                        vec_a = np.hstack((v_ung_a, -v_ung_a))
                        vec_b = np.hstack((v_ung_b, -v_ung_b))

                # build density

                if self.rank == mpi_master():
                    half_size_a = vec_a.shape[0] // 2
                    half_size_b = vec_b.shape[0] // 2

                    if getattr(self, 'core_excitation', False):
                        kn_a = self.lrvec2mat(vec_a, nocc_a, norb,
                                              self.num_core_orbitals)
                        kn_b = self.lrvec2mat(vec_b, nocc_b, norb,
                                              self.num_core_orbitals)
                        dak = self.commut_mo_density(kn_a, nocc_a,
                                                     self.num_core_orbitals)
                        dbk = self.commut_mo_density(kn_b, nocc_b,
                                                     self.num_core_orbitals)
                        core_exc_orb_inds_a = list(range(
                            self.num_core_orbitals)) + list(range(nocc_a, norb))
                        core_exc_orb_inds_b = list(range(
                            self.num_core_orbitals)) + list(range(nocc_b, norb))
                        mo_core_exc_a = mo_a[:, core_exc_orb_inds_a]
                        mo_core_exc_b = mo_b[:, core_exc_orb_inds_b]
                        dak = np.linalg.multi_dot(
                            [mo_core_exc_a, dak, mo_core_exc_a.T])
                        dbk = np.linalg.multi_dot(
                            [mo_core_exc_b, dbk, mo_core_exc_b.T])
                    else:
                        kn_a = self.lrvec2mat(vec_a, nocc_a, norb)
                        kn_b = self.lrvec2mat(vec_b, nocc_b, norb)
                        dak = self.commut_mo_density(kn_a, nocc_a)
                        dbk = self.commut_mo_density(kn_b, nocc_b)
                        dak = np.linalg.multi_dot([mo_a, dak, mo_a.T])
                        dbk = np.linalg.multi_dot([mo_b, dbk, mo_b.T])

                    dks_a.append(dak)
                    dks_b.append(dbk)
                    kns_a.append(kn_a)
                    kns_b.append(kn_b)

            if profiler is not None:
                profiler.add_timing_info('PreProc', tm.time() - prep_t0)

            # form Fock matrices

            fock = self._comp_lr_fock_unrestricted(
                (dks_a, dks_b), molecule, basis, eri_dict, dft_dict, pe_dict,
                profiler)

            if profiler is not None:
                # only increment FockCount on master rank
                if self.rank == mpi_master():
                    # Note: Jab/Ka/Kb -> 3 Fock builds for openshell
                    profiler.add_timing_info('_FockCount_', len(dks_a) * 3)

            prep_t0 = tm.time()

            if self.rank == mpi_master():
                raw_fock_ger_a = []
                raw_fock_ung_a = []
                raw_fock_ger_b = []
                raw_fock_ung_b = []

                raw_kns_ger_a = []
                raw_kns_ung_a = []
                raw_kns_ger_b = []
                raw_kns_ung_b = []

                for col in range(batch_start, batch_end):
                    ifock = col - batch_start
                    fock_a_np = fock[ifock * 2 + 0]
                    fock_b_np = fock[ifock * 2 + 1]

                    if col < n_general:
                        raw_fock_ger_a.append(0.5 * (fock_a_np - fock_a_np.T))
                        raw_fock_ung_a.append(0.5 * (fock_a_np + fock_a_np.T))
                        raw_kns_ger_a.append(0.5 *
                                             (kns_a[ifock] + kns_a[ifock].T))
                        raw_kns_ung_a.append(0.5 *
                                             (kns_a[ifock] - kns_a[ifock].T))
                        raw_fock_ger_b.append(0.5 * (fock_b_np - fock_b_np.T))
                        raw_fock_ung_b.append(0.5 * (fock_b_np + fock_b_np.T))
                        raw_kns_ger_b.append(0.5 *
                                             (kns_b[ifock] + kns_b[ifock].T))
                        raw_kns_ung_b.append(0.5 *
                                             (kns_b[ifock] - kns_b[ifock].T))

                    elif n_extra_ger > 0:
                        raw_fock_ger_a.append(0.5 * (fock_a_np - fock_a_np.T))
                        raw_kns_ger_a.append(0.5 *
                                             (kns_a[ifock] + kns_a[ifock].T))
                        raw_fock_ger_b.append(0.5 * (fock_b_np - fock_b_np.T))
                        raw_kns_ger_b.append(0.5 *
                                             (kns_b[ifock] + kns_b[ifock].T))

                    elif n_extra_ung > 0:
                        raw_fock_ung_a.append(0.5 * (fock_a_np + fock_a_np.T))
                        raw_kns_ung_a.append(0.5 *
                                             (kns_a[ifock] - kns_a[ifock].T))
                        raw_fock_ung_b.append(0.5 * (fock_b_np + fock_b_np.T))
                        raw_kns_ung_b.append(0.5 *
                                             (kns_b[ifock] - kns_b[ifock].T))

                fock_a = raw_fock_ger_a + raw_fock_ung_a
                kns_a = raw_kns_ger_a + raw_kns_ung_a
                fock_b = raw_fock_ger_b + raw_fock_ung_b
                kns_b = raw_kns_ger_b + raw_kns_ung_b

                batch_ger = len(raw_fock_ger_a)
                batch_ung = len(raw_fock_ung_a)

            e2_ger_a = None
            e2_ung_a = None
            e2_ger_b = None
            e2_ung_b = None

            fock_ger_a = None
            fock_ung_a = None
            fock_ger_b = None
            fock_ung_b = None

            if self.rank == mpi_master():

                e2_ger_a = np.zeros((half_size_a, batch_ger))
                e2_ung_a = np.zeros((half_size_a, batch_ung))
                e2_ger_b = np.zeros((half_size_b, batch_ger))
                e2_ung_b = np.zeros((half_size_b, batch_ung))

                if self.nonlinear:
                    fock_ger_a = np.zeros((norb**2, batch_ger))
                    fock_ung_a = np.zeros((norb**2, batch_ung))

                for ifock in range(batch_ger + batch_ung):
                    fak = fock_a[ifock]
                    fbk = fock_b[ifock]

                    if getattr(self, 'core_excitation', False):
                        core_exc_orb_inds_a = list(range(
                            self.num_core_orbitals)) + list(range(nocc_a, norb))
                        core_exc_orb_inds_b = list(range(
                            self.num_core_orbitals)) + list(range(nocc_b, norb))
                        mo_core_exc_a = mo_a[:, core_exc_orb_inds_a]
                        mo_core_exc_b = mo_b[:, core_exc_orb_inds_b]
                        fak_mo = np.linalg.multi_dot(
                            [mo_core_exc_a.T, fak, mo_core_exc_a])
                        fbk_mo = np.linalg.multi_dot(
                            [mo_core_exc_b.T, fbk, mo_core_exc_b])
                    else:
                        fak_mo = np.linalg.multi_dot([mo_a.T, fak, mo_a])
                        fbk_mo = np.linalg.multi_dot([mo_b.T, fbk, mo_b])

                    kfa_mo = self.commut(fa_mo.T, kns_a[ifock])
                    kfb_mo = self.commut(fb_mo.T, kns_b[ifock])

                    fat_mo = fak_mo + kfa_mo
                    fbt_mo = fbk_mo + kfb_mo

                    if getattr(self, 'core_excitation', False):
                        gmo_a = -self.commut_mo_density(fat_mo, nocc_a,
                                                        self.num_core_orbitals)
                        gmo_b = -self.commut_mo_density(fbt_mo, nocc_b,
                                                        self.num_core_orbitals)
                        gmo_vec_halfsize_a = self.lrmat2vec(
                            gmo_a, nocc_a, norb,
                            self.num_core_orbitals)[:half_size_a]
                        gmo_vec_halfsize_b = self.lrmat2vec(
                            gmo_b, nocc_b, norb,
                            self.num_core_orbitals)[:half_size_b]
                    else:
                        gmo_a = -self.commut_mo_density(fat_mo, nocc_a)
                        gmo_b = -self.commut_mo_density(fbt_mo, nocc_b)
                        gmo_vec_halfsize_a = self.lrmat2vec(
                            gmo_a, nocc_a, norb)[:half_size_a]
                        gmo_vec_halfsize_b = self.lrmat2vec(
                            gmo_b, nocc_b, norb)[:half_size_b]

                    # Note: fak_mo_vec uses full MO coefficients matrix since
                    # it is only for fock_ger/fock_ung in nonlinear response
                    fak_mo_vec = np.linalg.multi_dot([mo_a.T, fak,
                                                      mo_a]).reshape(norb**2)
                    fbk_mo_vec = np.linalg.multi_dot([mo_b.T, fbk,
                                                      mo_b]).reshape(norb**2)

                    if ifock < batch_ger:
                        e2_ger_a[:, ifock] = -gmo_vec_halfsize_a
                        e2_ger_b[:, ifock] = -gmo_vec_halfsize_b
                        if self.nonlinear:
                            fock_ger_a[:, ifock] = fak_mo_vec
                            fock_ger_b[:, ifock] = fbk_mo_vec
                    else:
                        e2_ung_a[:, ifock - batch_ger] = -gmo_vec_halfsize_a
                        e2_ung_b[:, ifock - batch_ger] = -gmo_vec_halfsize_b
                        if self.nonlinear:
                            fock_ung_a[:, ifock - batch_ger] = fak_mo_vec
                            fock_ung_b[:, ifock - batch_ger] = fbk_mo_vec

                # put alpha and beta together
                e2_ger = np.vstack((e2_ger_a, e2_ger_b))
                e2_ung = np.vstack((e2_ung_a, e2_ung_b))

            else:
                e2_ger = None
                e2_ung = None

            dist_e2_ger = DistributedArray(e2_ger, self.comm)
            dist_e2_ung = DistributedArray(e2_ung, self.comm)

            if vecs_e2_ger_data is None:
                vecs_e2_ger_data = dist_e2_ger.data.copy()
            else:
                vecs_e2_ger_data = np.hstack(
                    (vecs_e2_ger_data, dist_e2_ger.data))

            if vecs_e2_ung_data is None:
                vecs_e2_ung_data = dist_e2_ung.data.copy()
            else:
                vecs_e2_ung_data = np.hstack(
                    (vecs_e2_ung_data, dist_e2_ung.data))

            if self.nonlinear:
                # TODO: unrestricted for nonlinear
                # dist_fock_ger = DistributedArray(fock_ger, self.comm)
                # dist_fock_ung = DistributedArray(fock_ung, self.comm)
                # TODO: avoid incremental append
                # self._append_fock_matrices(dist_fock_ger, dist_fock_ung)
                pass

            if profiler is not None:
                profiler.add_timing_info('PostProc', tm.time() - prep_t0)

        prep_t0 = tm.time()

        # Note: append sigma and trial vectors only once

        vecs_e2_ger = DistributedArray(vecs_e2_ger_data,
                                       self.comm,
                                       distribute=False)

        vecs_e2_ung = DistributedArray(vecs_e2_ung_data,
                                       self.comm,
                                       distribute=False)

        self._append_sigma_vectors(vecs_e2_ger, vecs_e2_ung)

        self._append_trial_vectors(vecs_ger, vecs_ung)

        if profiler is not None:
            profiler.add_timing_info('AppendVec', tm.time() - prep_t0)

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
            self.checkpoint_file = f'{self.filename}_rsp.h5'

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

        need_exit = False

        if self.program_end_time is not None:
            remaining_hours = (self.program_end_time -
                               datetime.now()).total_seconds() / 3600
            # exit gracefully when the remaining time is not sufficient to
            # complete the next iteration (plus 25% to be on the safe side).
            if remaining_hours < next_iter_in_hours * 1.25:
                need_exit = True

        need_exit = self.comm.bcast(need_exit, root=mpi_master())

        return need_exit

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

        if self._cpcm:
            cur_str = 'Solvation Model                 : '
            cur_str += 'C-PCM'
            self.ostream.print_header(cur_str.ljust(str_width))
            cur_str = 'C-PCM Points per Hydrogen Sphere: '
            cur_str += f'{self.cpcm_grid_per_sphere[1]}'
            self.ostream.print_header(cur_str.ljust(str_width))
            cur_str = 'C-PCM Points per non-H Sphere   : '
            cur_str += f'{self.cpcm_grid_per_sphere[0]}'
            self.ostream.print_header(cur_str.ljust(str_width))
            cur_str = 'Non-Equilibrium solvation       : '
            cur_str += f'{self.non_equilibrium_solv}'
            self.ostream.print_header(cur_str.ljust(str_width))
            cur_str = 'C-PCM Dielectric Constant       : '
            cur_str += f'{self.cpcm_epsilon}'
            self.ostream.print_header(cur_str.ljust(str_width))
            if self.non_equilibrium_solv:
                cur_str = 'C-PCM Optical Dielectric Const. : '
                cur_str += f'{self.cpcm_optical_epsilon}'
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

    def get_prop_grad(self,
                      operator,
                      components,
                      molecule,
                      basis,
                      scf_tensors,
                      spin='alpha'):
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
                'dipole', 'electric dipole', 'electric_dipole', 'quadrupole',
                'electric quadrupole', 'electric_quadrupole', 'linear_momentum',
                'linear momentum', 'angular_momentum', 'angular momentum',
                'magnetic dipole', 'magnetic_dipole'
            ], f'LinearSolver.get_prop_grad: unsupported operator {operator}')

        integrals = None

        if operator in ['dipole', 'electric dipole', 'electric_dipole']:
            if self.rank == mpi_master():
                dipole_mats = compute_electric_dipole_integrals(
                    molecule, basis, [0.0, 0.0, 0.0])
                integrals = {
                    'x': dipole_mats[0],
                    'y': dipole_mats[1],
                    'z': dipole_mats[2],
                }

        elif operator in [
                'quadrupole', 'electric quadrupole', 'electric_quadrupole'
        ]:
            if self.rank == mpi_master():
                quadrupole_mats = compute_quadrupole_integrals(
                    molecule, basis, [0.0, 0.0, 0.0])
                integrals = {
                    'xx': quadrupole_mats[0],
                    'xy': quadrupole_mats[1],
                    'xz': quadrupole_mats[2],
                    'yy': quadrupole_mats[3],
                    'yz': quadrupole_mats[4],
                    'zz': quadrupole_mats[5],
                }

        elif operator in ['linear_momentum', 'linear momentum']:
            if self.rank == mpi_master():
                linmom_mats = compute_linear_momentum_integrals(molecule, basis)
                integrals = {
                    'x': -1.0 * linmom_mats[0],
                    'y': -1.0 * linmom_mats[1],
                    'z': -1.0 * linmom_mats[2],
                }

        elif operator in ['angular_momentum', 'angular momentum']:
            if self.rank == mpi_master():
                angmom_mats = compute_angular_momentum_integrals(
                    molecule, basis, [0.0, 0.0, 0.0])
                integrals = {
                    'x': -1.0 * angmom_mats[0],
                    'y': -1.0 * angmom_mats[1],
                    'z': -1.0 * angmom_mats[2],
                }

        elif operator in ['magnetic_dipole', 'magnetic dipole']:
            if self.rank == mpi_master():
                angmom_mats = compute_angular_momentum_integrals(
                    molecule, basis, [0.0, 0.0, 0.0])
                integrals = {
                    'x': 0.5 * angmom_mats[0],
                    'y': 0.5 * angmom_mats[1],
                    'z': 0.5 * angmom_mats[2],
                }

        # compute right-hand side

        if self.rank == mpi_master():
            integral_comps = [integrals[p] for p in components]

            if spin == 'alpha':
                mo = scf_tensors['C_alpha']
                nocc = molecule.number_of_alpha_electrons()
            elif spin == 'beta':
                mo = scf_tensors['C_beta']
                nocc = molecule.number_of_beta_electrons()
            norb = mo.shape[1]

            factor = np.sqrt(2.0)

            if getattr(self, 'core_excitation', False):
                core_exc_orb_inds = list(range(self.num_core_orbitals)) + list(
                    range(nocc, norb))
                mo_core_exc = mo[:, core_exc_orb_inds]
                matrices = [
                    factor * (-1.0) * self.commut_mo_density(
                        np.linalg.multi_dot([mo_core_exc.T, P.T, mo_core_exc]),
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
                'dipole', 'electric dipole', 'electric_dipole', 'quadrupole',
                'electric quadrupole', 'electric_quadrupole', 'linear_momentum',
                'linear momentum', 'angular_momentum', 'angular momentum',
                'magnetic dipole', 'magnetic_dipole'
            ],
            f'LinearSolver.get_complex_prop_grad: unsupported operator {operator}'
        )

        integrals = None

        if operator in ['dipole', 'electric dipole', 'electric_dipole']:
            if self.rank == mpi_master():
                dipole_mats = compute_electric_dipole_integrals(
                    molecule, basis, [0.0, 0.0, 0.0])
                integrals = {
                    'x': dipole_mats[0] + 0j,
                    'y': dipole_mats[1] + 0j,
                    'z': dipole_mats[2] + 0j,
                }

        elif operator in [
                'quadrupole', 'electric quadrupole', 'electric_quadrupole'
        ]:
            if self.rank == mpi_master():
                quadrupole_mats = compute_quadrupole_integrals(
                    molecule, basis, [0.0, 0.0, 0.0])
                integrals = {
                    'xx': quadrupole_mats[0] + 0j,
                    'xy': quadrupole_mats[1] + 0j,
                    'xz': quadrupole_mats[2] + 0j,
                    'yy': quadrupole_mats[3] + 0j,
                    'yz': quadrupole_mats[4] + 0j,
                    'zz': quadrupole_mats[5] + 0j,
                }

        elif operator in ['linear_momentum', 'linear momentum']:
            if self.rank == mpi_master():
                linmom_mats = compute_linear_momentum_integrals(molecule, basis)
                integrals = {
                    'x': -1j * linmom_mats[0],
                    'y': -1j * linmom_mats[1],
                    'z': -1j * linmom_mats[2],
                }

        elif operator in ['angular_momentum', 'angular momentum']:
            if self.rank == mpi_master():
                angmom_mats = compute_angular_momentum_integrals(
                    molecule, basis, [0.0, 0.0, 0.0])
                integrals = {
                    'x': -1j * angmom_mats[0],
                    'y': -1j * angmom_mats[1],
                    'z': -1j * angmom_mats[2],
                }

        elif operator in ['magnetic_dipole', 'magnetic dipole']:
            if self.rank == mpi_master():
                angmom_mats = compute_angular_momentum_integrals(
                    molecule, basis, [0.0, 0.0, 0.0])
                integrals = {
                    'x': 0.5j * angmom_mats[0],
                    'y': 0.5j * angmom_mats[1],
                    'z': 0.5j * angmom_mats[2],
                }

        # compute right-hand side

        if self.rank == mpi_master():
            integral_comps = [integrals[p] for p in components]

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

    def get_nto_unrestricted(self, t_mat, mo_occ, mo_vir):
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

        t_mat_a, t_mat_b = t_mat
        mo_occ_a, mo_occ_b = mo_occ
        mo_vir_a, mo_vir_b = mo_vir

        # SVD
        u_mat_a, s_diag_a, vh_mat_a = np.linalg.svd(t_mat_a, full_matrices=True)
        lam_diag_a = s_diag_a**2

        u_mat_b, s_diag_b, vh_mat_b = np.linalg.svd(t_mat_b, full_matrices=True)
        lam_diag_b = s_diag_b**2

        # holes in increasing order of lambda
        # particles in decreasing order of lambda
        nto_occ_a = np.flip(np.matmul(mo_occ_a, u_mat_a), axis=1)
        nto_vir_a = np.matmul(mo_vir_a, vh_mat_a.T)

        nto_occ_b = np.flip(np.matmul(mo_occ_b, u_mat_b), axis=1)
        nto_vir_b = np.matmul(mo_vir_b, vh_mat_b.T)

        # NTOs including holes and particles
        nto_orbs_a = np.concatenate((nto_occ_a, nto_vir_a), axis=1)
        nto_orbs_b = np.concatenate((nto_occ_b, nto_vir_b), axis=1)

        nto_lam_a = np.zeros(nto_orbs_a.shape[1])
        nocc_a = nto_occ_a.shape[1]
        for i_nto in range(lam_diag_a.size):
            nto_lam_a[nocc_a - 1 - i_nto] = -lam_diag_a[i_nto]
            nto_lam_a[nocc_a + i_nto] = lam_diag_a[i_nto]

        nto_lam_b = np.zeros(nto_orbs_b.shape[1])
        nocc_b = nto_occ_b.shape[1]
        for i_nto in range(lam_diag_b.size):
            nto_lam_b[nocc_b - 1 - i_nto] = -lam_diag_b[i_nto]
            nto_lam_b[nocc_b + i_nto] = lam_diag_b[i_nto]

        nto_mo = MolecularOrbitals.create_nto([nto_orbs_a, nto_orbs_b],
                                              [nto_lam_a, nto_lam_b],
                                              molorb.unrest)

        return nto_mo

    def write_nto_cubes(self,
                        cubic_grid,
                        molecule,
                        basis,
                        root,
                        nto_mo,
                        nto_pairs=None,
                        nto_thresh=0.1,
                        nto_spin=''):
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
            if nto_spin == '' or nto_spin == 'alpha':
                nocc = molecule.number_of_alpha_electrons()
            elif nto_spin == 'beta':
                nocc = molecule.number_of_beta_electrons()
        nvir = nto_mo.number_of_mos() - nocc

        if nto_spin == '' or nto_spin == 'alpha':
            lam_diag = nto_mo.occa_to_numpy()[nocc:nocc + min(nocc, nvir)]
            nto_coefs = nto_mo.alpha_to_numpy()
            cube_spin_str = 'alpha'
        elif nto_spin == 'beta':
            lam_diag = nto_mo.occb_to_numpy()[nocc:nocc + min(nocc, nvir)]
            nto_coefs = nto_mo.beta_to_numpy()
            cube_spin_str = 'beta'

        for i_nto in range(lam_diag.size):
            if lam_diag[i_nto] < nto_thresh:
                continue

            if (nto_pairs is not None) and (i_nto >= nto_pairs):
                continue

            self.ostream.print_info('  lambda: {:.4f}'.format(lam_diag[i_nto]))

            # hole
            ind_occ = nocc - i_nto - 1
            vis_drv.compute(cubic_grid, molecule, basis, nto_coefs, ind_occ)

            if self.rank == mpi_master():
                occ_cube_name = '{:s}_S{:d}_NTO_H{:d}{:s}.cube'.format(
                    base_fname, root + 1, i_nto + 1, nto_spin[:1])
                vis_drv.write_data(occ_cube_name, cubic_grid, molecule, 'nto',
                                   ind_occ, cube_spin_str)
                filenames.append(occ_cube_name)

                self.ostream.print_info(
                    '    Cube file (hole)     : {:s}'.format(occ_cube_name))
                self.ostream.flush()

            # electron
            ind_vir = nocc + i_nto
            vis_drv.compute(cubic_grid, molecule, basis, nto_coefs, ind_vir)

            if self.rank == mpi_master():
                vir_cube_name = '{:s}_S{:d}_NTO_P{:d}{:s}.cube'.format(
                    base_fname, root + 1, i_nto + 1, nto_spin[:1])
                vis_drv.write_data(vir_cube_name, cubic_grid, molecule, 'nto',
                                   ind_vir, cube_spin_str)
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

        if getattr(self, 'core_excitation', False):
            n_ov = self.num_core_orbitals * nvir
        else:
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

    def get_excitation_details_unrestricted(self,
                                            eigvec,
                                            nocc,
                                            nvir,
                                            coef_thresh=0.2):
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

        eigvec_a, eigvec_b = eigvec
        nocc_a, nocc_b = nocc
        nvir_a, nvir_b = nvir

        if getattr(self, 'core_excitation', False):
            n_ov_a = self.num_core_orbitals * nvir_a
            n_ov_b = self.num_core_orbitals * nvir_b
        else:
            n_ov_a = nocc_a * nvir_a
            n_ov_b = nocc_b * nvir_b

        assert_msg_critical(
            eigvec_a.size == n_ov_a or eigvec_a.size == n_ov_a * 2,
            'LinearSolver.get_excitation_details: Inconsistent size')

        assert_msg_critical(
            eigvec_b.size == n_ov_b or eigvec_b.size == n_ov_b * 2,
            'LinearSolver.get_excitation_details: Inconsistent size')

        excitations = []
        de_excitations = []

        for nocc, nvir, eigvec, spin in [(nocc_a, nvir_a, eigvec_a, 'a'),
                                         (nocc_b, nvir_b, eigvec_b, 'b')]:
            for i in range(nocc):
                if getattr(self, 'core_excitation', False):
                    homo_str = f'core_{i + 1}'
                else:
                    homo_str = 'HOMO' if i == nocc - 1 else f'HOMO-{nocc - 1 - i}'
                homo_str += f'({spin})'

                for a in range(nvir):
                    lumo_str = 'LUMO' if a == 0 else f'LUMO+{a}'
                    lumo_str += f'({spin})'

                    ia = i * nvir + a

                    exc_coef = eigvec[ia]
                    if abs(exc_coef) > coef_thresh:
                        excitations.append((
                            abs(exc_coef),
                            f'{homo_str:<12s} -> {lumo_str:<12s} {exc_coef:10.4f}',
                        ))

                    if eigvec.size == nocc * nvir * 2:
                        de_exc_coef = eigvec[nocc * nvir + ia]
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
        self.ostream.flush()

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
        self.ostream.flush()

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
        self.ostream.flush()

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

    @staticmethod
    def is_imag(op):
        """
        Checks if an operator is imaginary.

        :return:
            True if operator is imaginary, False otherwise
        """

        return op in [
            'linear momentum',
            'linear_momentum',
            'angular momentum',
            'angular_momentum',
            'magnetic dipole',
            'magnetic_dipole',
        ]

    @staticmethod
    def is_quadrupole(op):
        """
        Checks if an operator is quadrupole.

        :param op:
            The operator.

        :return:
            True if operator is quadrupole, False otherwise
        """

        return op in [
            'quadrupole',
            'electric quadrupole',
            'electric_quadrupole',
        ]

    def is_valid_component(self, comp, op):
        """
        Checks if a component is valid.

        :param comp:
            The component.
        :param op:
            The operator.

        :return:
            True if component is valid, False otherwise
        """

        if self.is_quadrupole(op):
            return comp in ['xx', 'xy', 'xz', 'yy', 'yz', 'zz']
        else:
            return comp in ['x', 'y', 'z']

    def plot_xas(self,
                 rsp_results,
                 broadening_type="lorentzian",
                 broadening_value=0.1239842,
                 ax=None):
        """
        Plot the X-ray absorption spectrum from the response calculation.

        :param rsp_results:
            The dictionary containing the linear response results.
        :param broadening_type:
            The type of broadening to use. Either 'lorentzian' or 'gaussian'.
        :param broadening_value:
            The broadening value in eV.
        :param ax:
            The matplotlib axis to plot on.
        """

        assert_msg_critical('matplotlib' in sys.modules,
                            'matplotlib is required.')

        au2ev = hartree_in_ev()
        ev2au = 1.0 / au2ev

        # initialize the plot
        if ax is None:
            fig, ax = plt.subplots(figsize=(8, 5))

        ax.set_xlabel('Photon energy [eV]')
        ax.set_ylabel(r'$\epsilon$ [L mol$^{-1}$ cm$^{-1}$]')

        ax.set_title("Absorption Spectrum")

        x = (rsp_results['eigenvalues'])
        y = rsp_results['oscillator_strengths']
        xmin = max(0.0, min(x) - 0.03)
        xmax = max(x) + 0.03
        xstep = 0.0001

        ax2 = ax.twinx()

        for i in np.arange(len(rsp_results['eigenvalues'])):
            ax2.plot(
                [
                    (rsp_results['eigenvalues'][i] * au2ev),
                    (rsp_results['eigenvalues'][i] * au2ev),
                ],
                [0.0, rsp_results['oscillator_strengths'][i]],
                alpha=0.7,
                linewidth=2,
                color="darkcyan",
            )

        c = 1.0 / fine_structure_constant()
        NA = avogadro_constant()
        a_0 = bohr_in_angstrom() * 1.0e-10

        if broadening_type.lower() == "lorentzian":
            xi, yi = self.lorentzian_absorption(x, y, xmin, xmax, xstep,
                                                broadening_value * ev2au)

        elif broadening_type.lower() == "gaussian":
            xi, yi = self.gaussian_absorption(x, y, xmin, xmax, xstep,
                                              broadening_value * ev2au)

        sigma = (2 * np.pi * np.pi * xi * yi) / c
        sigma_m2 = sigma * a_0**2
        sigma_cm2 = sigma_m2 * 10**4
        epsilon = sigma_cm2 * NA / (np.log(10) * 10**3)
        ax.plot(xi * au2ev,
                epsilon,
                color="black",
                alpha=0.9,
                linewidth=2.5)

        legend_bars = mlines.Line2D([], [],
                                    color='darkcyan',
                                    alpha=0.7,
                                    linewidth=2,
                                    label='Oscillator strength')
        label_spectrum = f'{broadening_type.capitalize()} '
        label_spectrum += f'broadening ({broadening_value:.3f} eV)'
        legend_spectrum = mlines.Line2D([], [],
                                        color='black',
                                        linestyle='-',
                                        linewidth=2.5,
                                        label=label_spectrum)
        ax2.legend(handles=[legend_bars, legend_spectrum],
                   frameon=False,
                   borderaxespad=0.,
                   loc='center left',
                   bbox_to_anchor=(1.15, 0.5))
        ax2.set_ylim(0, max(abs(rsp_results['oscillator_strengths'])) * 1.1)
        ax.set_ylim(0, max(epsilon) * 1.1)
        ax.set_ylim(bottom=0)
        ax2.set_ylim(bottom=0)
        ax2.set_ylabel("Oscillator strength")
        ax.set_xlim(xmin * au2ev, xmax * au2ev)

    def plot_uv_vis(self,
                    rsp_results,
                    broadening_type="lorentzian",
                    broadening_value=(1000.0 / hartree_in_wavenumber() *
                                      hartree_in_ev()),
                    ax=None):
        """
        Plot the UV spectrum from the response calculation.

        :param rsp_results:
            The dictionary containing linear response results.
        :param broadening_type:
            The type of broadening to use. Either 'lorentzian' or 'gaussian'.
        :param broadening_value:
            The broadening value in eV.
        :param ax:
            The matplotlib axis to plot on.
        """

        assert_msg_critical('matplotlib' in sys.modules,
                            'matplotlib is required.')

        ev_x_nm = hartree_in_ev() / hartree_in_inverse_nm()
        au2ev = hartree_in_ev()
        ev2au = 1.0 / au2ev

        # initialize the plot
        if ax is None:
            fig, ax = plt.subplots(figsize=(8, 5))

        ax.set_xlabel('Wavelength [nm]')
        ax.set_ylabel(r'$\epsilon$ [L mol$^{-1}$ cm$^{-1}$]')

        ax.set_title("Absorption Spectrum")

        x = (rsp_results['eigenvalues'])
        y = rsp_results['oscillator_strengths']
        xmin = max(0.0, min(x) - 0.03)
        xmax = max(x) + 0.03
        xstep = 0.0001

        ax2 = ax.twinx()

        for i in np.arange(len(rsp_results['eigenvalues'])):
            ax2.plot(
                [
                    ev_x_nm / (rsp_results['eigenvalues'][i] * au2ev),
                    ev_x_nm / (rsp_results['eigenvalues'][i] * au2ev),
                ],
                [0.0, rsp_results['oscillator_strengths'][i]],
                alpha=0.7,
                linewidth=2,
                color="darkcyan",
            )

        c = 1.0 / fine_structure_constant()
        NA = avogadro_constant()
        a_0 = bohr_in_angstrom() * 1.0e-10

        if broadening_type.lower() == "lorentzian":
            xi, yi = self.lorentzian_absorption(x, y, xmin, xmax, xstep,
                                                broadening_value * ev2au)

        elif broadening_type.lower() == "gaussian":
            xi, yi = self.gaussian_absorption(x, y, xmin, xmax, xstep,
                                              broadening_value * ev2au)

        sigma = (2 * np.pi * np.pi * xi * yi) / c
        sigma_m2 = sigma * a_0**2
        sigma_cm2 = sigma_m2 * 10**4
        epsilon = sigma_cm2 * NA / (np.log(10) * 10**3)
        ax.plot(ev_x_nm / (xi * au2ev),
                epsilon,
                color="black",
                alpha=0.9,
                linewidth=2.5)

        legend_bars = mlines.Line2D([], [],
                                    color='darkcyan',
                                    alpha=0.7,
                                    linewidth=2,
                                    label='Oscillator strength')
        label_spectrum = f'{broadening_type.capitalize()} '
        label_spectrum += f'broadening ({broadening_value:.3f} eV)'
        legend_spectrum = mlines.Line2D([], [],
                                        color='black',
                                        linestyle='-',
                                        linewidth=2.5,
                                        label=label_spectrum)
        ax2.legend(handles=[legend_bars, legend_spectrum],
                   frameon=False,
                   borderaxespad=0.,
                   loc='center left',
                   bbox_to_anchor=(1.15, 0.5))
        ax2.set_ylim(0, max(abs(rsp_results['oscillator_strengths'])) * 1.1)
        ax.set_ylim(0, max(epsilon) * 1.1)
        ax.set_ylim(bottom=0)
        ax2.set_ylim(bottom=0)
        ax2.set_ylabel("Oscillator strength")
        ax.set_xlim(ev_x_nm / (xmax * au2ev), ev_x_nm / (xmin * au2ev))

    def plot_xcd(self,
                 rsp_results,
                 broadening_type="lorentzian",
                 broadening_value=0.1239842,
                 ax=None):
        """
        Plot the X-ray CD spectrum from the response calculation.

        :param rsp_results:
            The dictionary containing linear response results.
        :param broadening_type:
            The type of broadening to use. Either 'lorentzian' or 'gaussian'.
        :param broadening_value:
            The broadening value in eV.
        :param ax:
            The matplotlib axis to plot on.
        """

        assert_msg_critical('matplotlib' in sys.modules,
                            'matplotlib is required.')

        au2ev = hartree_in_ev()

        # initialize the plot
        if ax is None:
            fig, ax = plt.subplots(figsize=(8, 5))

        ax.set_xlabel("Photon energy [eV]")
        ax.set_title("ECD Spectrum")
        ax.set_ylabel(r'$\Delta \epsilon$ [L mol$^{-1}$ cm$^{-1}$]')

        ax2 = ax.twinx()
        ax2.set_ylabel('Rotatory strength [10$^{-40}$ cgs]')

        for i in np.arange(len(rsp_results["eigenvalues"])):
            ax2.plot(
                [
                    rsp_results["eigenvalues"][i] * au2ev,
                    rsp_results["eigenvalues"][i] * au2ev,
                ],
                [0.0, rsp_results["rotatory_strengths"][i]],
                alpha=0.7,
                linewidth=2,
                color="darkcyan",
            )
        ax2.set_ylim(-max(abs(rsp_results["rotatory_strengths"])) * 1.1,
                     max(abs(rsp_results["rotatory_strengths"])) * 1.1)

        ax.axhline(y=0,
                   marker=',',
                   color='k',
                   linestyle='-.',
                   markersize=0,
                   linewidth=0.2)

        x = (rsp_results["eigenvalues"]) * au2ev
        y = rsp_results["rotatory_strengths"]
        xmin = max(0.0, min(x) - 0.8)
        xmax = max(x) + 0.8
        xstep = 0.003

        if broadening_type.lower() == "lorentzian":
            xi, yi = self.lorentzian_ecd(x, y, xmin, xmax, xstep,
                                         broadening_value)

        elif broadening_type.lower() == "gaussian":
            xi, yi = self.gaussian_ecd(x, y, xmin, xmax, xstep,
                                       broadening_value)

        # denorm_factor is roughly 22.96 * PI
        denorm_factor = (rotatory_strength_in_cgs() /
                         (extinction_coefficient_from_beta() / 3.0))
        yi = (yi * xi) / denorm_factor

        ax.set_ylim(-max(abs(yi)) * 1.1, max(abs(yi)) * 1.1)

        ax.plot(xi, yi, color="black", alpha=0.9, linewidth=2.5)
        ax.set_xlim(xmin, xmax)

        # include a legend for the bar and for the broadened spectrum
        legend_bars = mlines.Line2D([], [],
                                    color='darkcyan',
                                    alpha=0.7,
                                    linewidth=2,
                                    label='Rotatory strength')
        label_spectrum = f'{broadening_type.capitalize()} '
        label_spectrum += f'broadening ({broadening_value:.3f} eV)'
        legend_spectrum = mlines.Line2D([], [],
                                        color='black',
                                        linestyle='-',
                                        linewidth=2.5,
                                        label=label_spectrum)
        ax.legend(handles=[legend_bars, legend_spectrum],
                  frameon=False,
                  borderaxespad=0.,
                  loc='center left',
                  bbox_to_anchor=(1.15, 0.5))

    def plot_ecd(self,
                 rsp_results,
                 broadening_type="lorentzian",
                 broadening_value=(1000.0 / hartree_in_wavenumber() *
                                   hartree_in_ev()),
                 ax=None):
        """
        Plot the ECD spectrum from the response calculation.

        :param rsp_results:
            The dictionary containing linear response results.
        :param broadening_type:
            The type of broadening to use. Either 'lorentzian' or 'gaussian'.
        :param broadening_value:
            The broadening value in eV.
        :param ax:
            The matplotlib axis to plot on.
        """

        assert_msg_critical('matplotlib' in sys.modules,
                            'matplotlib is required.')

        ev_x_nm = hartree_in_ev() / hartree_in_inverse_nm()
        au2ev = hartree_in_ev()

        # initialize the plot
        if ax is None:
            fig, ax = plt.subplots(figsize=(8, 5))

        ax.set_xlabel("Wavelength [nm]")
        ax.set_title("ECD Spectrum")
        ax.set_ylabel(r'$\Delta \epsilon$ [L mol$^{-1}$ cm$^{-1}$]')

        ax2 = ax.twinx()
        ax2.set_ylabel('Rotatory strength [10$^{-40}$ cgs]')

        for i in np.arange(len(rsp_results["eigenvalues"])):
            ax2.plot(
                [
                    ev_x_nm / (rsp_results["eigenvalues"][i] * au2ev),
                    ev_x_nm / (rsp_results["eigenvalues"][i] * au2ev),
                ],
                [0.0, rsp_results["rotatory_strengths"][i]],
                alpha=0.7,
                linewidth=2,
                color="darkcyan",
            )
        ax2.set_ylim(-max(abs(rsp_results["rotatory_strengths"])) * 1.1,
                     max(abs(rsp_results["rotatory_strengths"])) * 1.1)

        ax.axhline(y=0,
                   marker=',',
                   color='k',
                   linestyle='-.',
                   markersize=0,
                   linewidth=0.2)

        x = (rsp_results["eigenvalues"]) * au2ev
        y = rsp_results["rotatory_strengths"]
        xmin = max(0.0, min(x) - 0.8)
        xmax = max(x) + 0.8
        xstep = 0.003

        if broadening_type.lower() == "lorentzian":
            xi, yi = self.lorentzian_ecd(x, y, xmin, xmax, xstep,
                                         broadening_value)

        elif broadening_type.lower() == "gaussian":
            xi, yi = self.gaussian_ecd(x, y, xmin, xmax, xstep,
                                       broadening_value)

        # denorm_factor is roughly 22.96 * PI
        denorm_factor = (rotatory_strength_in_cgs() /
                         (extinction_coefficient_from_beta() / 3.0))
        yi = (yi * xi) / denorm_factor

        ax.set_ylim(-max(abs(yi)) * 1.1, max(abs(yi)) * 1.1)

        ax.plot(ev_x_nm / xi, yi, color="black", alpha=0.9, linewidth=2.5)
        ax.set_xlim(ev_x_nm / xmax, ev_x_nm / xmin)

        # include a legend for the bar and for the broadened spectrum
        legend_bars = mlines.Line2D([], [],
                                    color='darkcyan',
                                    alpha=0.7,
                                    linewidth=2,
                                    label='Rotatory strength')
        label_spectrum = f'{broadening_type.capitalize()} '
        label_spectrum += f'broadening ({broadening_value:.3f} eV)'
        legend_spectrum = mlines.Line2D([], [],
                                        color='black',
                                        linestyle='-',
                                        linewidth=2.5,
                                        label=label_spectrum)
        ax.legend(handles=[legend_bars, legend_spectrum],
                  frameon=False,
                  borderaxespad=0.,
                  loc='center left',
                  bbox_to_anchor=(1.15, 0.5))

    @staticmethod
    def lorentzian_absorption(x, y, xmin, xmax, xstep, gamma):
        xi = np.arange(xmin, xmax, xstep)
        yi = np.zeros(len(xi))
        for i in range(len(xi)):
            for k in range(len(x)):
                yi[i] = yi[i] + y[k] / x[k] * gamma / (
                    (xi[i] - x[k])**2 + gamma**2)
        yi = yi / np.pi
        return xi, yi

    @staticmethod
    def gaussian_absorption(x, y, xmin, xmax, xstep, sigma):
        xi = np.arange(xmin, xmax, xstep)
        yi = np.zeros(len(xi))
        for i in range(len(xi)):
            for k in range(len(x)):
                yi[i] = yi[i] + y[k] / x[k] * np.exp(-((xi[i] - x[k])**2) /
                                                     (2 * sigma**2))
        yi = yi / (sigma * np.sqrt(2 * np.pi))
        return xi, yi

    @staticmethod
    def lorentzian_ecd(x, y, xmin, xmax, xstep, gamma):
        xi = np.arange(xmin, xmax, xstep)
        yi = np.zeros(len(xi))
        for i in range(len(xi)):
            for k in range(len(x)):
                yi[i] = yi[i] + y[k] * (gamma) / ((xi[i] - x[k])**2 +
                                                  (gamma)**2)
        return xi, yi

    @staticmethod
    def gaussian_ecd(x, y, xmin, xmax, xstep, sigma):
        xi = np.arange(xmin, xmax, xstep)
        yi = np.zeros(len(xi))
        for i in range(len(xi)):
            for k in range(len(x)):
                yi[i] = yi[i] + y[k] * np.exp(-((xi[i] - x[k])**2) /
                                              (2 * sigma**2))
        yi = np.pi * yi / (sigma * np.sqrt(2 * np.pi))
        return xi, yi

    def plot(self,
             rsp_results,
             broadening_type="lorentzian",
             broadening_value=(1000.0 / hartree_in_wavenumber() *
                               hartree_in_ev()),
             plot_type="electronic"):
        """
        Plot the absorption or ECD spectrum from the response calculation.

        :param rsp_results:
            The dictionary containing linear response results.
        :param broadening_type:
            The type of broadening to use. Either 'lorentzian' or 'gaussian'.
        :param broadening_value:
            The broadening value in eV.
        :param plot_type:
            The type of plot to generate. Either 'uv', 'xas', 'ecd', 'xcd', or 'electronic'.
        """

        assert_msg_critical('matplotlib' in sys.modules,
                            'matplotlib is required.')

        if plot_type.lower() in ["uv", "uv-vis", "uv_vis"]:
            self.plot_uv_vis(rsp_results,
                             broadening_type=broadening_type,
                             broadening_value=broadening_value)
        if plot_type.lower() in ["xas"]:
            self.plot_xas(rsp_results,
                          broadening_type=broadening_type,
                          broadening_value=broadening_value)

        elif plot_type.lower() == "ecd":
            self.plot_ecd(rsp_results,
                          broadening_type=broadening_type,
                          broadening_value=broadening_value)
        elif plot_type.lower() == "xcd":
            self.plot_xcd(rsp_results,
                          broadening_type=broadening_type,
                          broadening_value=broadening_value)

        elif plot_type.lower() == "electronic":
            fig, axs = plt.subplots(2, 1, figsize=(8, 10))
            # Increase the height space between subplots
            fig.subplots_adjust(hspace=0.3)
            if self.core_excitation:
                self.plot_xas(rsp_results,
                              broadening_type=broadening_type,
                              broadening_value=broadening_value,
                              ax=axs[0])
                self.plot_xcd(rsp_results,
                              broadening_type=broadening_type,
                              broadening_value=broadening_value,
                              ax=axs[1])
            else:
                self.plot_uv_vis(rsp_results,
                                 broadening_type=broadening_type,
                                 broadening_value=broadening_value,
                                 ax=axs[0])
                self.plot_ecd(rsp_results,
                              broadening_type=broadening_type,
                              broadening_value=broadening_value,
                              ax=axs[1])

        else:
            assert_msg_critical(False, 'Invalid plot type')

        plt.show()
