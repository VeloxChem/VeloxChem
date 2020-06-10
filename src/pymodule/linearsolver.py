import numpy as np
import time as tm
import itertools
import ctypes
import psutil
import sys
import os

from .veloxchemlib import ElectronRepulsionIntegralsDriver
from .veloxchemlib import ElectricDipoleIntegralsDriver
from .veloxchemlib import LinearMomentumIntegralsDriver
from .veloxchemlib import AngularMomentumIntegralsDriver
from .veloxchemlib import AODensityMatrix
from .veloxchemlib import AOFockMatrix
from .veloxchemlib import DenseMatrix
from .veloxchemlib import ExcitationVector
from .veloxchemlib import GridDriver
from .veloxchemlib import XCFunctional
from .veloxchemlib import XCIntegrator
from .veloxchemlib import MolecularGrid
from .veloxchemlib import mpi_master
from .veloxchemlib import denmat
from .veloxchemlib import fockmat
from .veloxchemlib import szblock
from .veloxchemlib import parse_xc_func
from .distributedarray import DistributedArray
from .subcommunicators import SubCommunicators
from .errorhandler import assert_msg_critical
from .qqscheme import get_qq_scheme
from .qqscheme import get_qq_type
from .checkpoint import write_rsp_hdf5


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
        - program_start_time: The start time of the program.
        - maximum_hours: The timelimit in hours.
        - timing: The flag for printing timing information.
        - profiling: The flag for printing profiling information.
        - memory_profiling: The flag for printing memory usage.
        - memory_tracing: The flag for tracing memory allocation.
        - dist_bger: The distributed gerade trial vectors.
        - dist_bung: The distributed ungerade trial vectors.
        - dist_e2bger: The distributed gerade sigma vectors.
        - dist_e2bung: The distributed ungerade sigma vectors.
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
        self.xcfun = XCFunctional()

        # polarizable embedding
        self.pe = False
        self.pe_options = {}

        # split communicators
        self.use_split_comm = False
        self.split_comm_ratio = None

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

        # information for graceful exit
        self.program_start_time = None
        self.maximum_hours = None

        # timing and profiling
        self.timing = False
        self.profiling = False
        self.memory_profiling = False
        self.memory_tracing = False

        self.dist_bger = None
        self.dist_bung = None
        self.dist_e2bger = None
        self.dist_e2bung = None

    def update_settings(self, rsp_dict, method_dict={}):
        """
        Updates response and method settings in linear solver.

        :param rsp_dict:
            The dictionary of response dict.
        :param method_dict:
            The dictionary of method rsp_dict.
        """

        if 'eri_thresh' in rsp_dict:
            self.eri_thresh = float(rsp_dict['eri_thresh'])
        if 'qq_type' in rsp_dict:
            self.qq_type = rsp_dict['qq_type']
        if 'batch_size' in rsp_dict:
            self.batch_size = int(rsp_dict['batch_size'])

        if 'dft' in method_dict:
            key = method_dict['dft'].lower()
            self.dft = True if key == 'yes' else False
        if 'grid_level' in method_dict:
            self.grid_level = int(method_dict['grid_level'])
        if 'xcfun' in method_dict:
            if 'dft' not in method_dict:
                self.dft = True
            self.xcfun = parse_xc_func(method_dict['xcfun'].upper())
            assert_msg_critical(not self.xcfun.is_undefined(),
                                'Undefined XC functional')

        if 'pe_options' not in method_dict:
            method_dict['pe_options'] = {}

        if 'pe' in method_dict:
            key = method_dict['pe'].lower()
            self.pe = True if key == 'yes' else False
        else:
            if ('potfile' in method_dict) or method_dict['pe_options']:
                self.pe = True

        if self.pe:
            if ('potfile' in method_dict and
                    'potfile' not in method_dict['pe_options']):
                method_dict['pe_options']['potfile'] = method_dict['potfile']
            assert_msg_critical('potfile' in method_dict['pe_options'],
                                'SCF driver: No potential file defined')
            self.pe_options = dict(method_dict['pe_options'])

        if 'use_split_comm' in method_dict:
            key = method_dict['use_split_comm'].lower()
            self.use_split_comm = True if key == 'yes' else False

        if 'conv_thresh' in rsp_dict:
            self.conv_thresh = float(rsp_dict['conv_thresh'])
        if 'max_iter' in rsp_dict:
            self.max_iter = int(rsp_dict['max_iter'])
        if 'lindep_thresh' in rsp_dict:
            self.lindep_thresh = float(rsp_dict['lindep_thresh'])

        if 'restart' in rsp_dict:
            key = rsp_dict['restart'].lower()
            self.restart = True if key == 'yes' else False
        if 'checkpoint_file' in rsp_dict:
            self.checkpoint_file = rsp_dict['checkpoint_file']

        if 'program_start_time' in rsp_dict:
            self.program_start_time = rsp_dict['program_start_time']
        if 'maximum_hours' in rsp_dict:
            self.maximum_hours = rsp_dict['maximum_hours']

        if 'timing' in rsp_dict:
            key = rsp_dict['timing'].lower()
            self.timing = True if key in ['yes', 'y'] else False
        if 'profiling' in rsp_dict:
            key = rsp_dict['profiling'].lower()
            self.profiling = True if key in ['yes', 'y'] else False

        if 'memory_profiling' in rsp_dict:
            key = rsp_dict['memory_profiling'].lower()
            self.memory_profiling = True if key in ['yes', 'y'] else False
        if 'memory_tracing' in rsp_dict:
            key = rsp_dict['memory_tracing'].lower()
            self.memory_tracing = True if key in ['yes', 'y'] else False

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

        # generate integration grid
        if self.dft:
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
                gs_density = AODensityMatrix([scf_tensors['D'][0]], denmat.rest)
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
            pe_drv = PolEmbed(molecule, basis, self.comm, self.pe_options)
            V_es = pe_drv.compute_multipole_potential_integrals()

            pot_info = "Reading polarizable embedding potential: {}".format(
                self.pe_options['potfile'])
            self.ostream.print_info(pot_info)
            self.ostream.print_blank()

            with open(self.pe_options['potfile'], 'r') as f_pot:
                potfile_text = os.linesep.join(f_pot.readlines())
        else:
            pe_drv = None
            V_es = None
            potfile_text = ''

        return {
            'pe_drv': pe_drv,
            'V_es': V_es,
            'potfile_text': potfile_text,
        }

    def read_vectors(self, rsp_vector_labels):
        """
        Reads vectors from checkpoint file.

        :param rsp_vector_labels:
            The list of labels of response vectors.
        """

        self.dist_bger, self.dist_bung, self.dist_e2bger, self.dist_e2bung = [
            DistributedArray.read_from_hdf5_file(self.checkpoint_file, label,
                                                 self.comm)
            for label in rsp_vector_labels
        ]
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

    def compute(self, molecule, basis, scf_tensors, v1=None):
        """
        Solves for the linear equations.

        :param molecule:
            The molecule.
        :param basis:
            The AO basis.
        :param scf_tensors:
            The dictionary of tensors from converged SCF wavefunction.
        :param v1:
            The gradients on the right-hand side. If not provided, v1 will be
            computed for the B operator.

        :return:
            A dictionary containing response functions, solutions, etc.
        """

        return None

    def e2n_half_size(self, vecs_ger, vecs_ung, molecule, basis, scf_tensors,
                      eri_dict, dft_dict, pe_dict, timing_dict):
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
        :param timing_dict:
            The dictionary containing timing information.

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
                'LinearResponse.e2n_half_size: '
                'invalid shape of trial vectors')

            assert_msg_critical(
                vecs_ger.shape(0) == vecs_ung.shape(0),
                'LinearResponse.e2n_half_size: '
                'inconsistent shape of trial vectors')

            mo = scf_tensors['C']
            S = scf_tensors['S']
            da, db = scf_tensors['D']
            fa, fb = scf_tensors['F']

            nocc = molecule.number_of_alpha_electrons()
            norb = mo.shape[1]

        # determine number of batches

        num_batches = 0
        batch_size = None

        total_mem = psutil.virtual_memory().total
        total_mem_list = self.comm.gather(total_mem, root=mpi_master())

        if self.rank == mpi_master():
            # check if master node has larger memory
            mem_adjust = 0.0
            if total_mem > min(total_mem_list):
                mem_adjust = total_mem - min(total_mem_list)

            # compute maximum batch size from available memory
            avail_mem = psutil.virtual_memory().available - mem_adjust
            mem_per_mat = mo.shape[0]**2 * ctypes.sizeof(ctypes.c_double)
            nthreads = int(os.environ['OMP_NUM_THREADS'])
            max_batch_size = int(avail_mem / mem_per_mat / (0.625 * nthreads))
            max_batch_size = max(1, max_batch_size)

            # set batch size
            batch_size = self.batch_size
            if batch_size is None:
                batch_size = min(100, n_total, max_batch_size)

            # get number of batches
            num_batches = n_total // batch_size
            if n_total % batch_size != 0:
                num_batches += 1
        batch_size = self.comm.bcast(batch_size, root=mpi_master())
        num_batches = self.comm.bcast(num_batches, root=mpi_master())

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
                              pe_dict, timing_dict)

            if self.rank == mpi_master():

                batch_ger, batch_ung = 0, 0
                for col in range(batch_start, batch_end):
                    if col < n_ger:
                        batch_ger += 1
                    else:
                        batch_ung += 1

                e2_ger = np.zeros((half_size, batch_ger))
                e2_ung = np.zeros((half_size, batch_ung))


                F_ger = []
                F_ung = []

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

                    gmo = np.linalg.multi_dot([mo.T, gao, mo])

                    if ifock < batch_ger:
                        e2_ger[:, ifock] = -self.lrmat2vec(gmo, nocc,
                                                           norb)[:half_size]
                        F_ger.append(np.linalg.multi_dot([mo.T, fak.T, mo]))

                    else:
                        e2_ung[:, ifock - batch_ger] = -self.lrmat2vec(
                            gmo, nocc, norb)[:half_size]
                        F_ung.append(np.linalg.multi_dot([mo.T, fak.T, mo]))

            else:
                e2_ger = None
                e2_ung = None

            vecs_e2_ger = DistributedArray(e2_ger, self.comm)
            vecs_e2_ung = DistributedArray(e2_ung, self.comm)

            #F_e2_ger = DistributedArray(F_ger,self.comm)
            #F_e2_ung = DistributedArray(F_ung,self.comm)

            self.append_sigma_vectors(vecs_e2_ger, vecs_e2_ung)

        self.append_trial_vectors(vecs_ger, vecs_ung)

        self.ostream.print_blank()
        if self.rank == mpi_master():
            return F_ger,F_ung
        else:
            return None,None

    def comp_lr_fock(self, fock, dens, molecule, basis, eri_dict, dft_dict,
                     pe_dict, timing_dict):
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
        :param timing_dict:
            The dictionary containing timing information.
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
            if timing_dict is not None:
                timing_dict['ERI'] = tm.time() - t0

            if self.dft:
                t0 = tm.time()
                if not self.xcfun.is_hybrid():
                    for ifock in range(fock.number_of_fock_matrices()):
                        fock.scale(2.0, ifock)
                xc_drv = XCIntegrator(self.comm)
                molgrid.distribute(self.rank, self.nodes, self.comm)
                xc_drv.integrate(fock, dens, gs_density, molecule, basis,
                                 molgrid, self.xcfun.get_func_label())
                if timing_dict is not None:
                    timing_dict['DFT'] = tm.time() - t0

            if self.pe:
                t0 = tm.time()
                pe_drv.V_es = V_es.copy()
                for ifock in range(fock.number_of_fock_matrices()):
                    dm = dens.alpha_to_numpy(ifock) + dens.beta_to_numpy(ifock)
                    e_pe, V_pe = pe_drv.get_pe_contribution(dm, elec_only=True)
                    if self.rank == mpi_master():
                        fock.add_matrix(DenseMatrix(V_pe), ifock)
                if timing_dict is not None:
                    timing_dict['PE'] = tm.time() - t0

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
            pe_drv = PolEmbed(molecule, basis, local_comm, self.pe_options)
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

        if self.maximum_hours is not None:
            remaining_hours = (self.maximum_hours -
                               (tm.time() - self.program_start_time) / 3600)
            # exit gracefully when the remaining time is not sufficient to
            # complete the next iteration (plus 25% to be on the safe side).
            if remaining_hours < next_iter_in_hours * 1.25:
                return True
        return False

    def print_header(self, title, nstates=None, n_points=None):
        """
        Prints linear response solver setup header to output stream.

        :param title:
            The name of the solver.
        :param nstates:
            The number of excited states.
        """

        self.ostream.print_blank()
        self.ostream.print_header('{:s} Setup'.format(title))
        self.ostream.print_header('=' * (len(title) + 8))
        self.ostream.print_blank()

        str_width = 60

        if nstates is not None:
            cur_str = 'Number of States                : ' + str(nstates)
            self.ostream.print_header(cur_str.ljust(str_width))

        if n_points is not None:
            cur_str = 'Number of integration points    : ' + str(n_points)
            self.ostream.print_header(cur_str.ljust(str_width))

        cur_str = 'Max. Number of Iterations       : ' + str(self.max_iter)
        self.ostream.print_header(cur_str.ljust(str_width))
        cur_str = 'Convergence Threshold           : ' + \
            '{:.1e}'.format(self.conv_thresh)
        self.ostream.print_header(cur_str.ljust(str_width))

        cur_str = 'ERI Screening Scheme            : ' + get_qq_type(
            self.qq_type)
        self.ostream.print_header(cur_str.ljust(str_width))
        cur_str = 'ERI Screening Threshold         : ' + \
            '{:.1e}'.format(self.eri_thresh)
        self.ostream.print_header(cur_str.ljust(str_width))
        if self.batch_size is not None:
            cur_str = 'Batch Size of Fock Matrices     : ' + \
                '{:d}'.format(self.batch_size)
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

    def initial_guess(self, v1, precond):

        return None

    def decomp_grad(self, grad):
        """
        Decomposes gradient into gerade and ungerade parts.

        :param grad:
            The gradient.

        :return:
            A tuple containing gerade and ungerade parts of gradient.
        """

        assert_msg_critical(grad.ndim == 1, 'decomp_grad: Expecting a 1D array')

        assert_msg_critical(grad.shape[0] % 2 == 0,
                            'decomp_grad: size of array should be even')

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

    def get_rhs(self, operator, components, molecule, basis, scf_tensors):
        """
        Creates right-hand side of linear response equations.

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
            The right-hand sides (gradients).
        """

        # compute 1e integral

        assert_msg_critical(
            operator in [
                'dipole', 'linear_momentum', 'linear momentum',
                'angular_momentum', 'angular momentum'
            ], 'get_rhs: unsupported operator {}'.format(operator))

        if operator == 'dipole':
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

        # compute right-hand side

        if self.rank == mpi_master():
            indices = {'x': 0, 'y': 1, 'z': 2}
            integral_comps = [integrals[indices[p]] for p in components]

            mo = scf_tensors['C']
            S = scf_tensors['S']
            D = scf_tensors['D'][0] + scf_tensors['D'][1]

            nocc = molecule.number_of_alpha_electrons()
            norb = mo.shape[1]

            matrices = tuple(
                np.linalg.multi_dot([
                    mo.T,
                    (np.linalg.multi_dot([S, D, P.T]) -
                     np.linalg.multi_dot([P.T, D, S])), mo
                ]) for P in integral_comps)

            gradients = tuple(self.lrmat2vec(m, nocc, norb) for m in matrices)
            return gradients

        else:
            return tuple()

    def get_complex_rhs(self, operator, components, molecule, basis,
                        scf_tensors):
        """
        Creates right-hand side of linear response equations.

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
            The right-hand sides (gradients).
        """

        # compute 1e integral

        assert_msg_critical(
            operator in [
                'dipole', 'linear_momentum', 'linear momentum',
                'angular_momentum', 'angular momentum'
            ], 'get_rhs: unsupported operator {}'.format(operator))

        if operator == 'dipole':
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

        # compute right-hand side

        if self.rank == mpi_master():
            indices = {'x': 0, 'y': 1, 'z': 2}
            integral_comps = [integrals[indices[p]] for p in components]

            mo = scf_tensors['C']
            S = scf_tensors['S']
            D = scf_tensors['D'][0] + scf_tensors['D'][1]

            nocc = molecule.number_of_alpha_electrons()
            norb = mo.shape[1]

            matrices = tuple(
                np.linalg.multi_dot([
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

        zlen = len(vec) // 2
        z, y = vec[:zlen], vec[zlen:]

        xv = ExcitationVector(szblock.aa, 0, nocc, nocc, norb, True)
        xv.set_yzcoefficients(z, y)

        kz = xv.get_zmatrix()
        ky = xv.get_ymatrix()

        rows = kz.number_of_rows() + ky.number_of_rows()
        cols = kz.number_of_columns() + ky.number_of_columns()

        kzy = np.zeros((rows, cols))
        kzy[:kz.number_of_rows(), ky.number_of_columns():] = kz.to_numpy()
        kzy[kz.number_of_rows():, :ky.number_of_columns()] = ky.to_numpy()

        return kzy

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

        xv = ExcitationVector(szblock.aa, 0, nocc, nocc, norb, True)
        excitations = list(
            itertools.product(xv.bra_unique_indexes(), xv.ket_unique_indexes()))

        z = [mat[i, j] for i, j in excitations]
        y = [mat[j, i] for i, j in excitations]
        return np.array(z + y)

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

    @staticmethod
    def remove_linear_dependence_half(basis, threshold):
        """
        Removes linear dependence in a set of symmetrized vectors.

        :param basis:
            The set of upper parts of symmetrized vectors.
        :param threshold:
            The threshold for removing linear dependence.

        :return:
            The new set of vectors.
        """

        Sb = 2 * np.matmul(basis.T, basis)
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

    @staticmethod
    def orthogonalize_gram_schmidt_half(tvecs):
        """
        Applies modified Gram Schmidt orthogonalization to trial vectors.

        :param tvecs:
            The trial vectors.

        :return:
            The orthogonalized trial vectors.
        """

        invsqrt2 = 1.0 / np.sqrt(2.0)

        if tvecs.shape[1] > 0:

            f = invsqrt2 / np.linalg.norm(tvecs[:, 0])
            tvecs[:, 0] *= f

            for i in range(1, tvecs.shape[1]):
                for j in range(i):
                    f = np.dot(tvecs[:, i], tvecs[:, j]) / np.dot(
                        tvecs[:, j], tvecs[:, j])
                    tvecs[:, i] -= f * tvecs[:, j]
                f = invsqrt2 / np.linalg.norm(tvecs[:, i])
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

    @staticmethod
    def normalize_half(vecs):
        """
        Normalizes half-sized vectors by dividing by vector norm.

        :param vecs:
            The half-sized vectors.

        :param Retruns:
            The normalized vectors.
        """

        invsqrt2 = 1.0 / np.sqrt(2.0)

        if len(vecs.shape) != 1:
            for vec in range(vecs.shape[1]):
                invnorm = invsqrt2 / np.linalg.norm(vecs[:, vec])
                vecs[:, vec] *= invnorm
        else:
            invnorm = invsqrt2 / np.linalg.norm(vecs)
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
    def construct_ed_sd(orb_ene, nocc, norb):
        """
        Gets the E0 and S0 diagonal elements as arrays.

        :param orb_ene:
            Orbital energies.
        :param nocc:
            Number of occupied orbitals.
        :param norb:
            Number of orbitals.

        :return:
            The E0 and S0 diagonal elements as numpy arrays.
        """

        xv = ExcitationVector(szblock.aa, 0, nocc, nocc, norb, True)
        excitations = list(
            itertools.product(xv.bra_unique_indexes(), xv.ket_unique_indexes()))

        z = [2.0 * (orb_ene[j] - orb_ene[i]) for i, j in excitations]
        ediag = np.array(z + z)

        lz = len(excitations)
        sdiag = 2.0 * np.ones(2 * lz)
        sdiag[lz:] = -2.0

        return ediag, sdiag

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

        xv = ExcitationVector(szblock.aa, 0, nocc, nocc, norb, True)
        excitations = list(
            itertools.product(xv.bra_unique_indexes(), xv.ket_unique_indexes()))

        z = [2.0 * (orb_ene[j] - orb_ene[i]) for i, j in excitations]
        ediag = np.array(z)

        lz = len(excitations)
        sdiag = 2.0 * np.ones(lz)

        return ediag, sdiag
