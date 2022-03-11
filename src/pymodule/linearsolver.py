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

from pathlib import Path
from datetime import datetime
import numpy as np
import time as tm
import sys

from .veloxchemlib import ElectronRepulsionIntegralsDriver
from .veloxchemlib import ElectricDipoleIntegralsDriver
from .veloxchemlib import LinearMomentumIntegralsDriver
from .veloxchemlib import AngularMomentumIntegralsDriver
from .veloxchemlib import AODensityMatrix
from .veloxchemlib import AOFockMatrix
from .veloxchemlib import DenseMatrix
from .veloxchemlib import GridDriver
from .veloxchemlib import XCIntegrator
from .veloxchemlib import MolecularGrid
from .veloxchemlib import mpi_master
from .veloxchemlib import denmat
from .veloxchemlib import fockmat
from .veloxchemlib import molorb
from .veloxchemlib import parse_xc_func
from .distributedarray import DistributedArray
from .subcommunicators import SubCommunicators
from .molecularorbitals import MolecularOrbitals
from .visualizationdriver import VisualizationDriver
from .errorhandler import assert_msg_critical
from .inputparser import parse_input, print_keywords, get_datetime_string
from .qqscheme import get_qq_scheme
from .qqscheme import get_qq_type
from .checkpoint import write_rsp_hdf5
from .batchsize import get_batch_size
from .batchsize import get_number_of_batches


class LinearSolver:
    """
    Implements linear solver.

    :param comm:
        The MPI communicator.
    :param use_split_comm:
        The flag for using split communicators.

    Instance variables
        - eri_thresh: The electron repulsion integrals screening threshold.
        - qq_type: The electron repulsion integrals screening scheme.
        - batch_size: The batch size for computation of Fock matrices.
        - dft: The flag for running DFT.
        - grid_level: The accuracy level of DFT grid.
        - xcfun: The XC functional.
        - pe: The flag for running polarizable embedding calculation.
        - pe_options: The dictionary with options for polarizable embedding.
        - use_split_comm: The flag for using split communicators.
        - split_comm_ratio: The list of ratios for split communicators.
        - electric_field: The static electric field.
        - conv_thresh: The convergence threshold for the solver.
        - max_iter: The maximum number of solver iterations.
        - cur_iter: Index of the current iteration.
        - small_thresh: The norm threshold for a vector to be considered a zero
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
        self.qq_type = 'QQ_DEN'
        self.batch_size = None

        # dft
        self.dft = False
        self.grid_level = 4
        self.xcfun = None

        # polarizable embedding
        self.pe = False
        self.pe_options = {}

        # split communicators
        self.use_split_comm = False
        self.split_comm_ratio = None

        # static electric field
        self.electric_field = None

        # solver setup
        self.conv_thresh = 1.0e-4
        self.max_iter = 150
        self.cur_iter = 0
        self.small_thresh = 1.0e-10
        self.lindep_thresh = 1.0e-6
        self.is_converged = False

        # mpi information
        self.comm = comm
        self.rank = self.comm.Get_rank()
        self.nodes = self.comm.Get_size()

        # output stream
        self.ostream = ostream

        # restart information
        self.restart = True
        self.checkpoint_file = None

        # timing and profiling
        self.timing = False
        self.profiling = False
        self.memory_profiling = False
        self.memory_tracing = False

        # program end time for graceful exit
        self.program_end_time = None

        # filename
        self.filename = f'veloxchem_rsp_{get_datetime_string()}'

        # distributed arrays
        self.dist_bger = None
        self.dist_bung = None
        self.dist_e2bger = None
        self.dist_e2bung = None

        # nonlinear flag and distributed Fock matrices
        self.nonlinear = False
        self.dist_fock_ger = None
        self.dist_fock_ung = None

        # input keywords
        self.input_keywords = {
            'response': {
                'eri_thresh': ('float', 'ERI screening threshold'),
                'qq_type': ('str_upper', 'ERI screening scheme'),
                'batch_size': ('int', 'batch size for Fock build'),
                'conv_thresh': ('float', 'convergence threshold'),
                'max_iter': ('int', 'maximum number of iterations'),
                'lindep_thresh': ('float', 'threshold for linear dependence'),
                'restart': ('bool', 'restart from checkpoint file'),
                'checkpoint_file': ('str', 'name of checkpoint file'),
                'timing': ('bool', 'print timing information'),
                'profiling': ('bool', 'print profiling information'),
                'memory_profiling': ('bool', 'print memory usage'),
                'memory_tracing': ('bool', 'trace memory allocation'),
            },
            'method_settings': {
                'dft': ('bool', 'use DFT'),
                'xcfun': ('str_upper', 'exchange-correlation functional'),
                'grid_level': ('int', 'accuracy level of DFT grid'),
                'pe': ('bool', 'use polarizable embedding'),
                'potfile': ('str', 'potential file for polarizable embedding'),
                'electric_field': ('seq_fixed', 'static electric field'),
                'use_split_comm': ('bool', 'use split communicators'),
            },
        }

    def print_keywords(self):
        """
        Prints input keywords in linear solver.
        """

        print_keywords(self.input_keywords, self.ostream)

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
            key: val[0] for key, val in self.input_keywords['response'].items()
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
            for key, val in self.input_keywords['method_settings'].items()
        }

        parse_input(self, method_keywords, method_dict)

        if 'xcfun' in method_dict:
            if 'dft' not in method_dict:
                self.dft = True
            self.xcfun = parse_xc_func(method_dict['xcfun'].upper())
            assert_msg_critical(not self.xcfun.is_undefined(),
                                'LinearSolver: Undefined XC functional')

        if 'pe_options' not in method_dict:
            method_dict['pe_options'] = {}

        if ('potfile' in method_dict) or method_dict['pe_options']:
            if 'pe' not in method_dict:
                self.pe = True

        if self.pe:
            from .polembed import PolEmbed

            if ('potfile' in method_dict and
                    'potfile' not in method_dict['pe_options']):
                method_dict['pe_options']['potfile'] = method_dict['potfile']
            assert_msg_critical('potfile' in method_dict['pe_options'],
                                'LinearSolver: No potential file defined')
            self.pe_options = dict(method_dict['pe_options'])

            cppe_potfile = None
            if self.rank == mpi_master():
                potfile = self.pe_options['potfile']
                if not Path(potfile).is_file():
                    potfile = str(
                        Path(self.filename).parent / Path(potfile).name)
                cppe_potfile = PolEmbed.write_cppe_potfile(potfile)
            cppe_potfile = self.comm.bcast(cppe_potfile, root=mpi_master())
            self.pe_options['potfile'] = cppe_potfile

        if self.electric_field is not None:
            assert_msg_critical(
                len(self.electric_field) == 3,
                'LinearSolver: Expecting 3 values in \'electric field\' input')
            assert_msg_critical(
                not self.pe,
                'LinearSolver: \'electric field\' input is incompatible ' +
                'with polarizable embedding')
            # disable restart of calculation with static electric field since
            # checkpoint file does not contain information about the electric
            # field
            self.restart = False

    def init_eri(self, molecule, basis):
        """
        Initializes ERI.

        :param molecule:
            The molecule.
        :param basis:
            The AO basis set.

        :return:
            The dictionary of ERI information.
        """

        if self.use_split_comm:
            self.use_split_comm = ((self.dft or self.pe) and self.nodes >= 8)

        if self.use_split_comm:
            screening = None
            valstr = 'ERI'
            if self.dft:
                valstr += '/DFT'
            if self.pe:
                valstr += '/PE'
            self.ostream.print_info(
                'Using sub-communicators for {}.'.format(valstr))
            self.ostream.print_blank()
        else:
            eri_drv = ElectronRepulsionIntegralsDriver(self.comm)
            screening = eri_drv.compute(get_qq_scheme(self.qq_type),
                                        self.eri_thresh, molecule, basis)

        return {
            'screening': screening,
        }

    def init_dft(self, molecule, scf_tensors):
        """
        Initializes DFT.

        :param molecule:
            The molecule.
        :param scf_tensors:
            The dictionary of tensors from converged SCF wavefunction.

        :return:
            The dictionary of DFT information.
        """

        if self.dft:
            # check dft setup
            assert_msg_critical(self.xcfun is not None,
                                'LinearSolver: Undefined XC functional')
            if isinstance(self.xcfun, str):
                self.xcfun = parse_xc_func(self.xcfun.upper())
                assert_msg_critical(not self.xcfun.is_undefined(),
                                    'LinearSolver: Undefined XC functional')

            # generate integration grid
            grid_drv = GridDriver(self.comm)
            grid_drv.set_level(self.grid_level)

            grid_t0 = tm.time()
            molgrid = grid_drv.generate(molecule)
            n_grid_points = molgrid.number_of_points()
            self.ostream.print_info(
                'Molecular grid with {0:d} points generated in {1:.2f} sec.'.
                format(n_grid_points,
                       tm.time() - grid_t0))
            self.ostream.print_blank()

            if self.rank == mpi_master():
                gs_density = AODensityMatrix([scf_tensors['D_alpha']],
                                             denmat.rest)
            else:
                gs_density = AODensityMatrix()
            gs_density.broadcast(self.rank, self.comm)

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

    def init_pe(self, molecule, basis):
        """
        Initializes polarizable embedding.

        :param molecule:
            The molecule.
        :param basis:
            The AO basis set.

        :return:
            The dictionary of polarizable embedding information.
        """

        # set up polarizable embedding
        if self.pe:
            from .polembed import PolEmbed
            pe_drv = PolEmbed(molecule, basis, self.pe_options, self.comm)
            V_es = pe_drv.compute_multipole_potential_integrals()

            cppe_info = 'Using CPPE {} for polarizable embedding.'.format(
                pe_drv.get_cppe_version())
            self.ostream.print_info(cppe_info)
            self.ostream.print_blank()

            pot_info = "Reading polarizable embedding potential: {}".format(
                self.pe_options['potfile'])
            self.ostream.print_info(pot_info)
            self.ostream.print_blank()

            with open(str(self.pe_options['potfile']), 'r') as f_pot:
                potfile_text = '\n'.join(f_pot.readlines())
        else:
            pe_drv = None
            V_es = None
            potfile_text = ''

        return {
            'pe_drv': pe_drv,
            'V_es': V_es,
            'potfile_text': potfile_text,
        }

    def read_checkpoint(self, rsp_vector_labels):
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
            (self.dist_bger, self.dist_bung, self.dist_e2bger, self.dist_e2bung,
             self.dist_fock_ger, self.dist_fock_ung) = dist_arrays
        else:
            (self.dist_bger, self.dist_bung, self.dist_e2bger,
             self.dist_e2bung) = dist_arrays

        checkpoint_text = 'Restarting from checkpoint file: '
        checkpoint_text += self.checkpoint_file
        self.ostream.print_info(checkpoint_text)
        self.ostream.print_blank()

    def append_trial_vectors(self, bger, bung):
        """
        Appends distributed trial vectors.

        :param bger:
            The distributed gerade trial vectors.
        :param bung:
            The distributed ungerade trial vectors.
        """

        if self.dist_bger is None:
            self.dist_bger = DistributedArray(bger.data,
                                              self.comm,
                                              distribute=False)
        else:
            self.dist_bger.append(bger, axis=1)

        if self.dist_bung is None:
            self.dist_bung = DistributedArray(bung.data,
                                              self.comm,
                                              distribute=False)
        else:
            self.dist_bung.append(bung, axis=1)

    def append_sigma_vectors(self, e2bger, e2bung):
        """
        Appends distributed sigma (E2 b) vectors.

        :param e2bger:
            The distributed gerade sigma vectors.
        :param e2bung:
            The distributed ungerade sigma vectors.
        """

        if self.dist_e2bger is None:
            self.dist_e2bger = DistributedArray(e2bger.data,
                                                self.comm,
                                                distribute=False)
        else:
            self.dist_e2bger.append(e2bger, axis=1)

        if self.dist_e2bung is None:
            self.dist_e2bung = DistributedArray(e2bung.data,
                                                self.comm,
                                                distribute=False)
        else:
            self.dist_e2bung.append(e2bung, axis=1)

    def append_fock_matrices(self, fock_ger, fock_ung):
        """
        Appends distributed Fock matrices in MO.

        :param fock_ger:
            The distributed gerade Fock matrices in MO.
        :param fock_ung:
            The distributed ungerade Fock matrices in MO.
        """

        if self.dist_fock_ger is None:
            self.dist_fock_ger = DistributedArray(fock_ger.data,
                                                  self.comm,
                                                  distribute=False)
        else:
            self.dist_fock_ger.append(fock_ger, axis=1)

        if self.dist_fock_ung is None:
            self.dist_fock_ung = DistributedArray(fock_ung.data,
                                                  self.comm,
                                                  distribute=False)
        else:
            self.dist_fock_ung.append(fock_ung, axis=1)

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

    def e2n_half_size(self,
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
                'LinearSolver.e2n_half_size: '
                'invalid shape of trial vectors')

            assert_msg_critical(
                vecs_ger.shape(0) == vecs_ung.shape(0),
                'LinearSolver.e2n_half_size: '
                'inconsistent shape of trial vectors')

            mo = scf_tensors['C_alpha']
            S = scf_tensors['S']
            da = scf_tensors['D_alpha']
            db = scf_tensors['D_beta']
            fa = scf_tensors['F_alpha']
            # fb = scf_tensors['F_beta']

            nocc = molecule.number_of_alpha_electrons()
            norb = mo.shape[1]

        # determine number of batches

        n_ao = mo.shape[0] if self.rank == mpi_master() else None

        batch_size = get_batch_size(self.batch_size, n_total, n_ao, self.comm)
        num_batches = get_number_of_batches(n_total, batch_size, self.comm)

        # go through batches

        if self.rank == mpi_master():
            batch_str = 'Processing Fock builds...'
            batch_str += ' (batch size: {:d})'.format(batch_size)
            self.ostream.print_info(batch_str)

        for batch_ind in range(num_batches):

            self.ostream.print_info('  batch {}/{}'.format(
                batch_ind + 1, num_batches))
            self.ostream.flush()

            # form density matrices

            batch_start = batch_size * batch_ind
            batch_end = min(batch_start + batch_size, n_total)

            if self.rank == mpi_master():
                dks = []
                kns = []

            for col in range(batch_start, batch_end):
                if col < n_ger:
                    v_ger = vecs_ger.get_full_vector(col)
                    if self.rank == mpi_master():
                        # full-size gerade trial vector
                        vec = np.hstack((v_ger, v_ger))
                else:
                    v_ung = vecs_ung.get_full_vector(col - n_ger)
                    if self.rank == mpi_master():
                        # full-size ungerade trial vector
                        vec = np.hstack((v_ung, -v_ung))

                if self.rank == mpi_master():
                    half_size = vec.shape[0] // 2

                    kN = self.lrvec2mat(vec, nocc, norb).T
                    kn = np.linalg.multi_dot([mo, kN, mo.T])

                    dak = np.linalg.multi_dot([kn.T, S, da])
                    dak -= np.linalg.multi_dot([da, S, kn.T])
                    # dbk = kn.T @ S @ db - db @ S @ kn.T

                    dks.append(dak)
                    kns.append(kn)

            if self.rank == mpi_master():
                dens = AODensityMatrix(dks, denmat.rest)
            else:
                dens = AODensityMatrix()
            dens.broadcast(self.rank, self.comm)

            # form Fock matrices

            fock = AOFockMatrix(dens)

            self.comp_lr_fock(fock, dens, molecule, basis, eri_dict, dft_dict,
                              pe_dict, profiler)

            e2_ger = None
            e2_ung = None

            fock_ger = None
            fock_ung = None

            if self.rank == mpi_master():

                batch_ger, batch_ung = 0, 0
                for col in range(batch_start, batch_end):
                    if col < n_ger:
                        batch_ger += 1
                    else:
                        batch_ung += 1

                e2_ger = np.zeros((half_size, batch_ger))
                e2_ung = np.zeros((half_size, batch_ung))

                if self.nonlinear:
                    fock_ger = np.zeros((norb**2, batch_ger))
                    fock_ung = np.zeros((norb**2, batch_ung))

                for ifock in range(batch_ger + batch_ung):
                    fak = fock.alpha_to_numpy(ifock).T
                    # fbk = fock.beta_to_numpy(ifock).T

                    kn = kns[ifock]

                    kfa = (np.linalg.multi_dot([S, kn, fa]) -
                           np.linalg.multi_dot([fa, kn, S]))
                    # kfb = S @ kn @ fb - fb @ kn @ S

                    fat = fak + kfa
                    fbt = fat
                    # fbt = fbk + kfb

                    gao = np.matmul(S,
                                    np.matmul(da, fat.T) + np.matmul(db, fbt.T))
                    gao -= np.matmul(
                        np.matmul(fat.T, da) + np.matmul(fbt.T, db), S)
                    gao *= 0.5

                    gmo = np.linalg.multi_dot([mo.T, gao, mo])

                    gmo_vec_halfsize = self.lrmat2vec(gmo, nocc,
                                                      norb)[:half_size]
                    fak_T_mo_vec = np.linalg.multi_dot([mo.T, fak.T,
                                                        mo]).reshape(norb**2)

                    if ifock < batch_ger:
                        e2_ger[:, ifock] = -gmo_vec_halfsize
                        if self.nonlinear:
                            fock_ger[:, ifock] = fak_T_mo_vec
                    else:
                        e2_ung[:, ifock - batch_ger] = -gmo_vec_halfsize
                        if self.nonlinear:
                            fock_ung[:, ifock - batch_ger] = fak_T_mo_vec

            vecs_e2_ger = DistributedArray(e2_ger, self.comm)
            vecs_e2_ung = DistributedArray(e2_ung, self.comm)
            self.append_sigma_vectors(vecs_e2_ger, vecs_e2_ung)

            if self.nonlinear:
                dist_fock_ger = DistributedArray(fock_ger, self.comm)
                dist_fock_ung = DistributedArray(fock_ung, self.comm)
                self.append_fock_matrices(dist_fock_ger, dist_fock_ung)

        self.append_trial_vectors(vecs_ger, vecs_ung)

        self.ostream.print_blank()

    def comp_lr_fock(self,
                     fock,
                     dens,
                     molecule,
                     basis,
                     eri_dict,
                     dft_dict,
                     pe_dict,
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

        screening = eri_dict['screening']

        molgrid = dft_dict['molgrid']
        gs_density = dft_dict['gs_density']

        V_es = pe_dict['V_es']
        pe_drv = pe_dict['pe_drv']

        # set flags for Fock matrices

        fock_flag = fockmat.rgenjk
        if self.dft:
            if self.xcfun.is_hybrid():
                fock_flag = fockmat.rgenjkx
                fact_xc = self.xcfun.get_frac_exact_exchange()
                for i in range(fock.number_of_fock_matrices()):
                    fock.set_scale_factor(fact_xc, i)
            else:
                fock_flag = fockmat.rgenj
        for i in range(fock.number_of_fock_matrices()):
            fock.set_fock_type(fock_flag, i)

        # calculate Fock on subcommunicators

        if self.use_split_comm:
            self.comp_lr_fock_split_comm(fock, dens, molecule, basis, eri_dict,
                                         dft_dict, pe_dict)

        else:
            t0 = tm.time()
            eri_drv = ElectronRepulsionIntegralsDriver(self.comm)
            eri_drv.compute(fock, dens, molecule, basis, screening)
            if profiler is not None:
                profiler.add_timing_info('ERI', tm.time() - t0)

            if self.dft:
                t0 = tm.time()
                if not self.xcfun.is_hybrid():
                    for ifock in range(fock.number_of_fock_matrices()):
                        fock.scale(2.0, ifock)
                xc_drv = XCIntegrator(self.comm)
                molgrid.distribute(self.rank, self.nodes, self.comm)
                xc_drv.integrate(fock, dens, gs_density, molecule, basis,
                                 molgrid, self.xcfun.get_func_label())
                if profiler is not None:
                    profiler.add_timing_info('DFT', tm.time() - t0)

            if self.pe:
                t0 = tm.time()
                pe_drv.V_es = V_es.copy()
                for ifock in range(fock.number_of_fock_matrices()):
                    dm = dens.alpha_to_numpy(ifock) + dens.beta_to_numpy(ifock)
                    e_pe, V_pe = pe_drv.get_pe_contribution(dm, elec_only=True)
                    if self.rank == mpi_master():
                        fock.add_matrix(DenseMatrix(V_pe), ifock)
                if profiler is not None:
                    profiler.add_timing_info('PE', tm.time() - t0)

            fock.reduce_sum(self.rank, self.nodes, self.comm)

    def comp_lr_fock_split_comm(self, fock, dens, molecule, basis, eri_dict,
                                dft_dict, pe_dict):
        """
        Computes linear response Fock/Fxc matrix on split communicators.

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
        """

        molgrid = dft_dict['molgrid']
        gs_density = dft_dict['gs_density']

        V_es = pe_dict['V_es']
        pe_drv = pe_dict['pe_drv']

        if self.split_comm_ratio is None:
            if self.dft and self.pe:
                self.split_comm_ratio = [0.34, 0.33, 0.33]
            elif self.dft:
                self.split_comm_ratio = [0.5, 0.5, 0.0]
            elif self.pe:
                self.split_comm_ratio = [0.5, 0.0, 0.5]
            else:
                self.split_comm_ratio = [1.0, 0.0, 0.0]

        if self.dft:
            dft_nodes = int(float(self.nodes) * self.split_comm_ratio[1] + 0.5)
            dft_nodes = max(1, dft_nodes)
        else:
            dft_nodes = 0

        if self.pe:
            pe_nodes = int(float(self.nodes) * self.split_comm_ratio[2] + 0.5)
            pe_nodes = max(1, pe_nodes)
        else:
            pe_nodes = 0

        eri_nodes = max(1, self.nodes - dft_nodes - pe_nodes)

        if eri_nodes == max(eri_nodes, dft_nodes, pe_nodes):
            eri_nodes = self.nodes - dft_nodes - pe_nodes
        elif dft_nodes == max(eri_nodes, dft_nodes, pe_nodes):
            dft_nodes = self.nodes - eri_nodes - pe_nodes
        else:
            pe_nodes = self.nodes - eri_nodes - dft_nodes

        node_grps = [0] * eri_nodes + [1] * dft_nodes + [2] * pe_nodes
        eri_comm = (node_grps[self.rank] == 0)
        dft_comm = (node_grps[self.rank] == 1)
        pe_comm = (node_grps[self.rank] == 2)

        subcomms = SubCommunicators(self.comm, node_grps)
        local_comm = subcomms.local_comm
        cross_comm = subcomms.cross_comm

        # reset molecular grid for DFT and V_es for PE
        if self.rank != mpi_master():
            molgrid = MolecularGrid()
            V_es = np.zeros(0)
        if self.dft:
            if local_comm.Get_rank() == mpi_master():
                molgrid.broadcast(cross_comm.Get_rank(), cross_comm)
        if self.pe:
            if local_comm.Get_rank() == mpi_master():
                V_es = cross_comm.bcast(V_es, root=mpi_master())

        t0 = tm.time()

        # calculate Fock on ERI nodes
        if eri_comm:
            eri_drv = ElectronRepulsionIntegralsDriver(local_comm)
            local_screening = eri_drv.compute(get_qq_scheme(self.qq_type),
                                              self.eri_thresh, molecule, basis)
            eri_drv.compute(fock, dens, molecule, basis, local_screening)
            if self.dft and not self.xcfun.is_hybrid():
                for ifock in range(fock.number_of_fock_matrices()):
                    fock.scale(2.0, ifock)

        # calculate Fxc on DFT nodes
        if dft_comm:
            xc_drv = XCIntegrator(local_comm)
            molgrid.distribute(local_comm.Get_rank(), local_comm.Get_size(),
                               local_comm)
            xc_drv.integrate(fock, dens, gs_density, molecule, basis, molgrid,
                             self.xcfun.get_func_label())

        # calculate e_pe and V_pe on PE nodes
        if pe_comm:
            from .polembed import PolEmbed
            pe_drv = PolEmbed(molecule, basis, self.pe_options, local_comm)
            pe_drv.V_es = V_es.copy()
            for ifock in range(fock.number_of_fock_matrices()):
                dm = dens.alpha_to_numpy(ifock) + dens.beta_to_numpy(ifock)
                e_pe, V_pe = pe_drv.get_pe_contribution(dm, elec_only=True)
                if local_comm.Get_rank() == mpi_master():
                    fock.add_matrix(DenseMatrix(V_pe), ifock)

        dt = tm.time() - t0

        # collect Fock on master node
        fock.reduce_sum(self.rank, self.nodes, self.comm)

        if local_comm.Get_rank() == mpi_master():
            dt = cross_comm.gather(dt, root=mpi_master())

        if self.rank == mpi_master():
            time_eri = dt[0] * eri_nodes
            time_dft = 0.0
            if self.dft:
                time_dft = dt[1] * dft_nodes
            time_pe = 0.0
            if self.pe:
                pe_root = 2 if self.dft else 1
                time_pe = dt[pe_root] * pe_nodes
            time_sum = time_eri + time_dft + time_pe
            self.split_comm_ratio = [
                time_eri / time_sum,
                time_dft / time_sum,
                time_pe / time_sum,
            ]
        self.split_comm_ratio = self.comm.bcast(self.split_comm_ratio,
                                                root=mpi_master())

    def write_checkpoint(self, molecule, basis, dft_dict, pe_dict, labels):
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

        if self.checkpoint_file is None:
            return

        if self.rank == mpi_master():
            success = write_rsp_hdf5(self.checkpoint_file, [], [], molecule,
                                     basis, dft_dict, pe_dict, self.ostream)
        else:
            success = False
        success = self.comm.bcast(success, root=mpi_master())

        if success:
            if self.nonlinear:
                dist_arrays = [
                    self.dist_bger, self.dist_bung, self.dist_e2bger,
                    self.dist_e2bung, self.dist_fock_ger, self.dist_fock_ung
                ]
            else:
                dist_arrays = [
                    self.dist_bger, self.dist_bung, self.dist_e2bger,
                    self.dist_e2bung
                ]

            for dist_array, label in zip(dist_arrays, labels):
                dist_array.append_to_hdf5_file(self.checkpoint_file, label)

    def graceful_exit(self, molecule, basis, dft_dict, pe_dict, labels):
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

        self.write_checkpoint(molecule, basis, dft_dict, pe_dict, labels)

        self.ostream.print_info('...done.')
        self.ostream.print_blank()
        self.ostream.print_info('Exiting program.')
        self.ostream.print_blank()
        self.ostream.flush()

        sys.exit(0)

    def need_graceful_exit(self, next_iter_in_hours):
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

    def print_header(self, title, nstates=None, n_freqs=None, n_points=None):
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

        cur_str = 'ERI Screening Scheme            : ' + get_qq_type(
            self.qq_type)
        self.ostream.print_header(cur_str.ljust(str_width))
        cur_str = 'ERI Screening Threshold         : {:.1e}'.format(
            self.eri_thresh)
        self.ostream.print_header(cur_str.ljust(str_width))
        if self.batch_size is not None:
            cur_str = 'Batch Size of Fock Matrices     : {:d}'.format(
                self.batch_size)
            self.ostream.print_header(cur_str.ljust(str_width))

        if self.dft:
            cur_str = 'Exchange-Correlation Functional : '
            cur_str += self.xcfun.get_func_label().upper()
            self.ostream.print_header(cur_str.ljust(str_width))
            cur_str = 'Molecular Grid Level            : ' + str(
                self.grid_level)
            self.ostream.print_header(cur_str.ljust(str_width))

        self.ostream.print_blank()
        self.ostream.flush()

    def print_iteration(self, relative_residual_norm, xvs):

        return

    def print_convergence(self, title):
        """
        Prints information after convergence.

        :param title:
            The name of the solver.
        """

        width = 92
        output_conv = '*** '
        if self.is_converged:
            output_conv += '{:s} converged'.format(title)
        else:
            output_conv += '{:s} NOT converged'.format(title)
        output_conv += ' in {:d} iterations. '.format(self.cur_iter + 1)
        output_conv += 'Time: {:.2f} sec'.format(tm.time() - self.start_time)
        self.ostream.print_header(output_conv.ljust(width))
        self.ostream.print_blank()
        self.ostream.print_blank()
        self.ostream.flush()

    def check_convergence(self, relative_residual_norm):
        """
        Checks convergence.

        :param relative_residual_norm:
            Relative residual norms.
        """

        if self.rank == mpi_master():
            max_residual = max(relative_residual_norm.values())
            if max_residual < self.conv_thresh:
                self.is_converged = True

        self.is_converged = self.comm.bcast(self.is_converged,
                                            root=mpi_master())

    def decomp_grad(self, grad):
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

    def get_precond(self, orb_ene, nocc, norb, w):

        return None

    def preconditioning(self, precond, v_in):

        return None

    def precond_trials(self, vectors, precond):
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

    def setup_trials(self,
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

        dist_new_ger, dist_new_ung = self.precond_trials(vectors, precond)

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
                dist_new_ger = self.remove_linear_dependence_half_distributed(
                    dist_new_ger, self.lindep_thresh)
                dist_new_ger = self.orthogonalize_gram_schmidt_half_distributed(
                    dist_new_ger)
                dist_new_ger = self.normalize_half_distributed(dist_new_ger)

            if dist_new_ung.data.ndim > 0 and dist_new_ung.shape(0) > 0:
                dist_new_ung = self.remove_linear_dependence_half_distributed(
                    dist_new_ung, self.lindep_thresh)
                dist_new_ung = self.orthogonalize_gram_schmidt_half_distributed(
                    dist_new_ung)
                dist_new_ung = self.normalize_half_distributed(dist_new_ung)

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
            dipole_drv = ElectricDipoleIntegralsDriver(self.comm)
            dipole_mats = dipole_drv.compute(molecule, basis)

            if self.rank == mpi_master():
                integrals = (dipole_mats.x_to_numpy(), dipole_mats.y_to_numpy(),
                             dipole_mats.z_to_numpy())
            else:
                integrals = tuple()

        elif operator in ['linear_momentum', 'linear momentum']:
            linmom_drv = LinearMomentumIntegralsDriver(self.comm)
            linmom_mats = linmom_drv.compute(molecule, basis)

            if self.rank == mpi_master():
                integrals = (linmom_mats.x_to_numpy(), linmom_mats.y_to_numpy(),
                             linmom_mats.z_to_numpy())
            else:
                integrals = tuple()

        elif operator in ['angular_momentum', 'angular momentum']:
            angmom_drv = AngularMomentumIntegralsDriver(self.comm)
            angmom_mats = angmom_drv.compute(molecule, basis)

            if self.rank == mpi_master():
                integrals = (angmom_mats.x_to_numpy(), angmom_mats.y_to_numpy(),
                             angmom_mats.z_to_numpy())
            else:
                integrals = tuple()

        elif operator in ['magnetic_dipole', 'magnetic dipole']:
            angmom_drv = AngularMomentumIntegralsDriver(self.comm)
            angmom_mats = angmom_drv.compute(molecule, basis)

            if self.rank == mpi_master():
                integrals = (-0.5 * angmom_mats.x_to_numpy(),
                             -0.5 * angmom_mats.y_to_numpy(),
                             -0.5 * angmom_mats.z_to_numpy())
            else:
                integrals = tuple()

        # compute right-hand side

        if self.rank == mpi_master():
            indices = {'x': 0, 'y': 1, 'z': 2}
            integral_comps = [integrals[indices[p]] for p in components]

            mo = scf_tensors['C_alpha']
            S = scf_tensors['S']
            D = scf_tensors['D_alpha'] + scf_tensors['D_beta']

            nocc = molecule.number_of_alpha_electrons()
            norb = mo.shape[1]

            factor = 0.5 * np.sqrt(2.0)
            matrices = tuple(factor * np.linalg.multi_dot([
                mo.T,
                (np.linalg.multi_dot([S, D, P.T]) -
                 np.linalg.multi_dot([P.T, D, S])), mo
            ]) for P in integral_comps)

            gradients = tuple(self.lrmat2vec(m, nocc, norb) for m in matrices)
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
            dipole_drv = ElectricDipoleIntegralsDriver(self.comm)
            dipole_mats = dipole_drv.compute(molecule, basis)

            if self.rank == mpi_master():
                integrals = (dipole_mats.x_to_numpy() + 0j,
                             dipole_mats.y_to_numpy() + 0j,
                             dipole_mats.z_to_numpy() + 0j)
            else:
                integrals = tuple()

        elif operator in ['linear_momentum', 'linear momentum']:
            linmom_drv = LinearMomentumIntegralsDriver(self.comm)
            linmom_mats = linmom_drv.compute(molecule, basis)

            if self.rank == mpi_master():
                integrals = (-1j * linmom_mats.x_to_numpy(),
                             -1j * linmom_mats.y_to_numpy(),
                             -1j * linmom_mats.z_to_numpy())
            else:
                integrals = tuple()

        elif operator in ['angular_momentum', 'angular momentum']:
            angmom_drv = AngularMomentumIntegralsDriver(self.comm)
            angmom_mats = angmom_drv.compute(molecule, basis)

            if self.rank == mpi_master():
                integrals = (-1j * angmom_mats.x_to_numpy(),
                             -1j * angmom_mats.y_to_numpy(),
                             -1j * angmom_mats.z_to_numpy())
            else:
                integrals = tuple()

        elif operator in ['magnetic_dipole', 'magnetic dipole']:
            angmom_drv = AngularMomentumIntegralsDriver(self.comm)
            angmom_mats = angmom_drv.compute(molecule, basis)

            if self.rank == mpi_master():
                integrals = (0.5j * angmom_mats.x_to_numpy(),
                             0.5j * angmom_mats.y_to_numpy(),
                             0.5j * angmom_mats.z_to_numpy())
            else:
                integrals = tuple()

        # compute right-hand side

        if self.rank == mpi_master():
            indices = {'x': 0, 'y': 1, 'z': 2}
            integral_comps = [integrals[indices[p]] for p in components]

            mo = scf_tensors['C_alpha']
            S = scf_tensors['S']
            D = scf_tensors['D_alpha'] + scf_tensors['D_beta']

            nocc = molecule.number_of_alpha_electrons()
            norb = mo.shape[1]

            factor = 0.5 * np.sqrt(2.0)
            matrices = tuple(factor * np.linalg.multi_dot([
                mo.T,
                (np.linalg.multi_dot([S, D, P.conj().T]) -
                 np.linalg.multi_dot([P.conj().T, D, S])), mo
            ]) for P in integral_comps)

            gradients = tuple(self.lrmat2vec(m, nocc, norb) for m in matrices)
            return gradients

        else:
            return tuple()

    @staticmethod
    def lrvec2mat(vec, nocc, norb):
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
        n_ov = nocc * nvir
        mat = np.zeros((norb, norb), dtype=vec.dtype)

        # excitation and de-excitation
        mat[:nocc, nocc:] = vec[:n_ov].reshape(nocc, nvir)
        mat[nocc:, :nocc] = vec[n_ov:].reshape(nocc, nvir).T

        return mat

    @staticmethod
    def lrmat2vec(mat, nocc, norb):
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

    def remove_linear_dependence_half_distributed(self, basis, threshold):
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

    def orthogonalize_gram_schmidt_half_distributed(self, tvecs):
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

    def normalize_half_distributed(self, vecs):
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
    def construct_ed_sd_half(orb_ene, nocc, norb):
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
            The lambda values (1D array) and the NTO coefficients (2D array).
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

        nto_lam = np.concatenate((np.flip(lam_diag),lam_diag))
        nto_ener = np.zeros(nto_orbs.shape[1])
        nto_mo = MolecularOrbitals([nto_orbs], [nto_ener], [nto_lam], molorb.rest)

        return lam_diag, nto_mo

    def write_nto_cubes(self,
                        cubic_grid,
                        molecule,
                        basis,
                        root,
                        lam_diag,
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
        :param lam_diag:
            The lambda values (1D array).
        :param nto_mo:
            The NTO coefficients (2D array).
        :param nto_pairs:
            The number of NTO pairs.
        :param nto_thresh:
            The threshold for writing NTO to cube file.

        :return:
            The list containing the names of the cube files.
        """

        filenames = []

        vis_drv = VisualizationDriver(self.comm)

        nocc = molecule.number_of_alpha_electrons()

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
                    self.filename, root + 1, i_nto + 1)
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
                    self.filename, root + 1, i_nto + 1)
                vis_drv.write_data(vir_cube_name, cubic_grid, molecule, 'nto',
                                   ind_vir, 'alpha')
                filenames.append(vir_cube_name)

                self.ostream.print_info(
                    '    Cube file (particle) : {:s}'.format(vir_cube_name))
                self.ostream.flush()

        self.ostream.print_blank()

        return filenames

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

        vis_drv = VisualizationDriver(self.comm)

        vis_drv.compute(cubic_grid, molecule, basis, dens_DA, 0, 'alpha')

        if self.rank == mpi_master():
            detach_cube_name = '{:s}_S{:d}_detach.cube'.format(
                self.filename, root + 1)
            vis_drv.write_data(detach_cube_name, cubic_grid, molecule,
                               'detachment', 0, 'alpha')
            filenames.append(detach_cube_name)

            self.ostream.print_info(
                '  Cube file (detachment) : {:s}'.format(detach_cube_name))
            self.ostream.flush()

        vis_drv.compute(cubic_grid, molecule, basis, dens_DA, 1, 'alpha')

        if self.rank == mpi_master():
            attach_cube_name = '{:s}_S{:d}_attach.cube'.format(
                self.filename, root + 1)
            vis_drv.write_data(attach_cube_name, cubic_grid, molecule,
                               'attachment', 1, 'alpha')
            filenames.append(attach_cube_name)

            self.ostream.print_info(
                '  Cube file (attachment) : {:s}'.format(attach_cube_name))
            self.ostream.flush()

        self.ostream.print_blank()

        return filenames
