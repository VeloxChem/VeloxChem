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
import numpy as np
import time as tm

from .veloxchemlib import ElectronRepulsionIntegralsDriver
from .veloxchemlib import AODensityMatrix
from .veloxchemlib import AOFockMatrix
from .veloxchemlib import GridDriver
from .veloxchemlib import XCFunctional
from .veloxchemlib import MolecularGrid
from .veloxchemlib import mpi_master
from .veloxchemlib import denmat
from .veloxchemlib import fockmat
from .linearsolver import LinearSolver
from .distributedarray import DistributedArray
from .errorhandler import assert_msg_critical
from .inputparser import parse_input, print_keywords, get_datetime_string
from .qqscheme import get_qq_scheme
from .batchsize import get_batch_size
from .batchsize import get_number_of_batches
from .veloxchemlib import parse_xc_func
from .veloxchemlib import XCIntegrator

class NonLinearSolver:
    """
    Implements nonlinear solver.

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
        - filename: The filename.
        - program_end_time: The end time of the program.
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

        # static electric field
        self.electric_field = None

        # solver setup
        self.conv_thresh = 1.0e-4
        self.max_iter = 150
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

        # input keywords
        self.input_keywords = {
            'response': {
                'eri_thresh': ('float', 'ERI screening threshold'),
                'qq_type': ('str_upper', 'ERI screening scheme'),
                'batch_size': ('int', 'batch size for Fock build'),
                'max_iter': ('int', 'maximum number of iterations'),
                'conv_thresh': ('float', 'convergence threshold'),
                'lindep_thresh': ('float', 'threshold for linear dependence'),
                'restart': ('bool', 'restart from checkpoint file'),
                'checkpoint_file': ('str', 'name of checkpoint file'),
                'timing': ('bool', 'print timing information'),
                'profiling': ('bool', 'print profiling information'),
                'memory_profiling': ('bool', 'print memory usage'),
                'memory_tracing': ('bool', 'trace memory allocation'),
            },
        }

    def print_keywords(self):
        """
        Prints input keywords in nonlinear solver.
        """

        print_keywords(self.input_keywords, self.ostream)

    def update_settings(self, rsp_dict, method_dict=None):
        """
        Updates response and method settings in nonlinear solver.

        :param rsp_dict:
            The dictionary of response dict.
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

        if 'xcfun' in method_dict:
            if 'dft' not in method_dict:
                self.dft = True
            self.xcfun = parse_xc_func(method_dict['xcfun'].upper())

            assert_msg_critical(not self.xcfun.is_undefined(),
                                'Nonlinear solver: Undefined XC functional')

        if 'potfile' in method_dict:
            errmsg = 'NonLinearSolver: The \'potfile\' keyword is not supported '
            errmsg += 'in nonlinear response calculation.'
            if self.rank == mpi_master():
                assert_msg_critical(False, errmsg)

        if 'electric_field' in method_dict:
            errmsg = 'NonLinearSolver: The \'electric field\' keyword is not '
            errmsg += 'supported in nonlinear response calculation.'
            if self.rank == mpi_master():
                assert_msg_critical(False, errmsg)

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

    def compute(self, molecule, basis, scf_tensors, v1=None):
        """
        Solves for the nonlinear response functions.

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

    def comp_nlr_fock(self, mo, molecule, ao_basis, fock_flag, dft_dict = None, d_dft_1 = None, d_dft_2 = None, mode=None):
        """-
        Computes and returns a list of Fock matrices 

        :param mo:
            The MO coefficients
        :param D:
            A list of densities
        :param molecule:
            The molecule
        :param ao_basis:
            The AO basis set
        :param fock_flag:
            The type of Fock matrices

        :return:
            A list of Fock matrices
        """

        if fock_flag == 'real_and_imag':

            if d_dft_1 and d_dft_2:
                f_total = self.comp_two_el_int(mo, molecule, ao_basis,dft_dict,d_dft_1,d_dft_2,mode)
            if d_dft_1 and not d_dft_2 :
                f_total = self.comp_two_el_int(mo, molecule, ao_basis,dft_dict,d_dft_1,mode)

            nrows = f_total.data.shape[0]
            half_ncols = f_total.data.shape[1] // 2
            ff_data = np.zeros((nrows, half_ncols), dtype=np.complex128)
            for i in range(half_ncols):
                ff_data[:, i] = (f_total.data[:, 2 * i] +
                                 1j * f_total.data[:, 2 * i + 1])
            return DistributedArray(ff_data, self.comm, distribute=False)

        elif fock_flag == 'real':
            return self.comp_two_el_int(mo, molecule, ao_basis, d_hf)

        else:
            return None

    def comp_two_el_int(self, mo, molecule, ao_basis , dft_dict = None, dens_1 = None, dens_2 = None,mode=None):
        """
        Returns the two-electron part of the Fock matix in MO basis

        :param mo:
            The MO coefficients
        :param molecule:
            The molecule
        :param ao_basis:
            The AO basis set
        :param dabs:
            A list of densitiy matrices

        :return:
            A tuple containing the two-electron part of the Fock matix (in MO
            basis)
        """

        eri_driver = ElectronRepulsionIntegralsDriver(self.comm)
        screening = eri_driver.compute(get_qq_scheme(self.qq_type),
                                       self.eri_thresh, molecule, ao_basis)

        # determine number of batches

        if self.rank == mpi_master():
            n_total = len(dens_2)
            n_ao = dens_2[0].shape[0]
            norb = mo.shape[1]
        else:
            n_total = None
            n_ao = None
        
        batch_size = get_batch_size(self.batch_size, n_total, n_ao, self.comm)
        num_batches = get_number_of_batches(n_total, batch_size, self.comm)

        # go through batches

        dist_fabs = None

        if self.rank == mpi_master():
            batch_str = 'Processing Fock builds...'
            batch_str += ' (batch size: {:d})'.format(batch_size)
            self.ostream.print_info(batch_str)

        for batch_ind in range(num_batches):

            if self.rank == mpi_master():
                self.ostream.print_info('  batch {}/{}'.format(
                    batch_ind + 1, num_batches))
                self.ostream.flush()

            # form density matrices

            if self.rank == mpi_master():
                batch_start = batch_size * batch_ind
                batch_end = min(batch_start + batch_size, n_total)

                if dens_1:
                    dts1 = [
                    np.ascontiguousarray(dab)
                    for dab in dens_1 ]
                    dens1 = AODensityMatrix(dts1, denmat.rest)

                if dens_2:
                    dts2 = [
                    np.ascontiguousarray(dab)
                    for dab in dens_2 ]
                    dens2 = AODensityMatrix(dts2, denmat.rest)
                
            else:
                dens1 = AODensityMatrix()
                dens2 = AODensityMatrix()
            
            if dens_1:
                dens1.broadcast(self.rank, self.comm)
            if dens_2:
                dens2.broadcast(self.rank, self.comm)

            fock = AOFockMatrix(dens2)
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

            eri_driver.compute(fock, dens2, molecule, ao_basis, screening)
            if self.dft and not self.xcfun.is_hybrid():
                for ifock in range(fock.number_of_fock_matrices()):
                    fock.scale(2.0, ifock)
            if self.dft: 
                xc_drv = XCIntegrator(self.comm)
                molgrid = dft_dict['molgrid']
                gs_density = dft_dict['gs_density']

                molgrid.distribute(self.rank, self.nodes, self.comm)
                xc_drv.integrate(fock, dens1, dens2, gs_density, molecule, ao_basis,
                                 molgrid, self.xcfun.get_func_label(),mode)

            fock.reduce_sum(self.rank, self.nodes, self.comm)

            if self.rank == mpi_master():
                nfocks = fock.number_of_fock_matrices()
                fock_mo = np.zeros((norb**2, nfocks))
                for i in range(nfocks):
                    fock_mo[:, i] = self.ao2mo(mo,
                                               fock.to_numpy(i).T).reshape(-1)
            else:
                fock_mo = None

            dist_fock_mo = DistributedArray(fock_mo, self.comm)

            if dist_fabs is None:
                dist_fabs = DistributedArray(dist_fock_mo.data,
                                             self.comm,
                                             distribute=False)
            else:
                dist_fabs.append(dist_fock_mo, axis=1)

        self.ostream.print_blank()

        return dist_fabs


    def comp_nlr_fock_cubic(self, mo, D, molecule, ao_basis, fock_flag):
        """-
        Computes and returns a list of Fock matrices.

        :param mo:
            The MO coefficients
        :param D:
            A list of densities
        :param molecule:
            The molecule
        :param ao_basis:
            The AO basis set
        :param fock_flag:
            The type of Fock matrices

        :return:
            A list of Fock matrices
        """

        if fock_flag == 'real_and_imag':
            if self.rank == mpi_master():
                D_total = []
                for da in D:
                    D_total.append(da.real)
                    D_total.append(da.imag)
            else:
                D_total = None

            f_total = self.comp_two_el_int_cubic(mo, molecule, ao_basis, D_total)

            nrows = f_total.data.shape[0]
            half_ncols = f_total.data.shape[1] // 2
            ff_data = np.zeros((nrows, half_ncols), dtype=np.complex128)
            for i in range(half_ncols):
                ff_data[:, i] = (f_total.data[:, 2 * i] +
                                 1j * f_total.data[:, 2 * i + 1])
            return DistributedArray(ff_data, self.comm, distribute=False)

        elif fock_flag == 'real':
            return self.comp_two_el_int_cubic(mo, molecule, ao_basis, D)

        else:
            return None

    def comp_two_el_int_cubic(self, mo, molecule, ao_basis, dabs):
        """
        Returns the two-electron part of the Fock matix in MO basis

        :param mo:
            The MO coefficients
        :param molecule:
            The molecule
        :param ao_basis:
            The AO basis set
        :param dabs:
            A list of densitiy matrices

        :return:
            A tuple containing the two-electron part of the Fock matix (in MO
            basis)
        """

        eri_driver = ElectronRepulsionIntegralsDriver(self.comm)
        screening = eri_driver.compute(get_qq_scheme(self.qq_type),
                                       self.eri_thresh, molecule, ao_basis)

        # determine number of batches

        if self.rank == mpi_master():
            n_total = len(dabs)
            n_ao = dabs[0].shape[0]
            norb = mo.shape[1]
        else:
            n_total = None
            n_ao = None

        batch_size = get_batch_size(self.batch_size, n_total, n_ao, self.comm)
        num_batches = get_number_of_batches(n_total, batch_size, self.comm)

        # go through batches

        dist_fabs = None

        if self.rank == mpi_master():
            batch_str = 'Processing Fock builds...'
            batch_str += ' (batch size: {:d})'.format(batch_size)
            self.ostream.print_info(batch_str)

        for batch_ind in range(num_batches):

            if self.rank == mpi_master():
                self.ostream.print_info('  batch {}/{}'.format(
                    batch_ind + 1, num_batches))
                self.ostream.flush()

            # form density matrices

            if self.rank == mpi_master():
                batch_start = batch_size * batch_ind
                batch_end = min(batch_start + batch_size, n_total)
                dts = [
                    np.ascontiguousarray(dab)
                    for dab in dabs[batch_start:batch_end]
                ]
                dens = AODensityMatrix(dts, denmat.rest)
            else:
                dens = AODensityMatrix()

            dens.broadcast(self.rank, self.comm)

            fock = AOFockMatrix(dens)
            for i in range(fock.number_of_fock_matrices()):
                fock.set_fock_type(fockmat.rgenjk, i)

            eri_driver.compute(fock, dens, molecule, ao_basis, screening)
            fock.reduce_sum(self.rank, self.nodes, self.comm)

            if self.rank == mpi_master():
                nfocks = fock.number_of_fock_matrices()
                fock_mo = np.zeros((norb**2, nfocks))
                for i in range(nfocks):
                    fock_mo[:, i] = self.ao2mo(mo,
                                               fock.to_numpy(i).T).reshape(-1)
            else:
                fock_mo = None

            dist_fock_mo = DistributedArray(fock_mo, self.comm)

            if dist_fabs is None:
                dist_fabs = DistributedArray(dist_fock_mo.data,
                                             self.comm,
                                             distribute=False)
            else:
                dist_fabs.append(dist_fock_mo, axis=1)

        self.ostream.print_blank()

        return dist_fabs

    def flip_yz(self, X):
        """
        This method takes a first-order response vector with a given sign of
        the frequency and returns the first-order response vector with reversed
        frequency argument.

        :param X:
            A response vector N(ω,x) = (Z,-Y^*)

        :return:
            A response vector with reversed optical frequency N(-ω,x) =
            (Y,-Z^*)
        """

        if X.ndim == 1:
            new_yz = np.zeros_like(X)
            half_len = X.shape[0] // 2
            new_yz[:half_len] = -X.real[half_len:] + 1j * X.imag[half_len:]
            new_yz[half_len:] = -X.real[:half_len] + 1j * X.imag[:half_len]
            return new_yz

        return None

    def print_fock_header(self):
        """
        Prints header for Fock computation
        """

        title = 'Fock Matrix Computation'
        self.ostream.print_header(title)
        self.ostream.print_header('=' * (len(title) + 2))
        self.ostream.print_blank()

    def print_fock_time(self, time):
        """
        Prints time for Fock computation

        :param time:
            Total time to compute Fock matrices
        """

        cur_str = 'Time spent in Fock matrices: {:.2f} sec'.format(time)
        self.ostream.print_info(cur_str)
        self.ostream.print_blank()
        self.ostream.flush()

    def print_component(self, label, value, width):
        """
        Prints response function components.

        :param label:
            The label
        :param freq:
            The frequency
        :param value:
            The complex value
        :param width:
            The width for the output
        """

        w_str = '{:<9s} {:20.8f} {:20.8f}j'.format(label, value.real,
                                                   value.imag)
        self.ostream.print_header(w_str.ljust(width))

    def s4(self, k1, k2, k3, D, nocc, norb):
        """
        Used for the contraction of the S[4] tensor with three response vectors

        :param k1:
            A response matrix
        :param k2:
            A response matrix
        :param k3:
            A response matrix
        :param D:
            A density matrix
        :param nocc:
            The number of occupied orbtials
        :param norb:
            The number of total orbitals

        :return:
            The contraction of S[4] for S[4] dict
        """

        S4_123 = self.s4_contract(k1, k2, k3, D, nocc, norb)
        S4_132 = self.s4_contract(k1, k3, k2, D, nocc, norb)

        return S4_123 + S4_132

    def s4_contract(self, k1, k2, k3, D, nocc, norb):
        """
        Returns the contraction of the S[4] tensor 

        :param k1:
            A response matrix
        :param k2:
            A response matrix
        :param k3:
            A response matrix
        :param D:
            A density matrix
        :param nocc:
            The number of occupied orbtials
        :param norb:
            The number of total orbitals

        :return:
            The contraction of S[4] for S[4] dict
        """

        S4N1N2N3 = self.commut(self.commut(k3, self.commut(k2, k1)), D.T)
        S4N1N2N3_c = self.complex_lrmat2vec(S4N1N2N3, nocc, norb)
        return (2. / 6) * S4N1N2N3_c

    def flip_xy(self, X):
        """
        Swaps upper and lower parts of a response vector.

        :param X:
            A response vector v = (Z,-Y^*)

        :return:
            A response vector of the form v' = (-Y^*,Z)
        """

        if X.ndim == 1:
            new_xy = np.zeros_like(X)
            half_len = X.shape[0] // 2
            new_xy[:half_len] = X[half_len:]
            new_xy[half_len:] = X[:half_len]
            return new_xy

        return None

    def s4_for_r4(self, k1, k2, k3, D, nocc, norb):
        """
        Returns the contraction of S[4] for the contraction of R[4]

        :param k1:
            A response matrix
        :param k2:
            A response matrix
        :param k3:
            A response matrix
        :param D:
            A density matrix
        :param nocc:
            The number of occupied orbtials
        :param norb:
            The number of total orbitals

        :return:
            The contraction of S[4] for the contraction of R[4]
        """

        S4_123 = self.s4_contract(k1, k2, k3, D, nocc, norb)

        return S4_123

    def collect_vectors_in_columns(self, sendbuf):
        """
        Collects vectors into 2d array (column-wise).

        :param sendbuf:
            The 2d array containing the vector segments in columns.

        :return:
            A 2d array containing the full vectors in columns.
        """

        counts = self.comm.gather(sendbuf.size, root=mpi_master())
        if self.rank == mpi_master():
            displacements = [sum(counts[:p]) for p in range(self.nodes)]
            recvbuf = np.zeros(sum(counts), dtype=sendbuf.dtype).reshape(
                -1, sendbuf.shape[1])
        else:
            displacements = None
            recvbuf = None

        if sendbuf.dtype == np.float64:
            mpi_data_type = MPI.DOUBLE
        elif sendbuf.dtype == np.complex128:
            mpi_data_type = MPI.C_DOUBLE_COMPLEX

        self.comm.Gatherv(sendbuf,
                          [recvbuf, counts, displacements, mpi_data_type],
                          root=mpi_master())

        return recvbuf

    def transform_dens(self, k, D, S):
        """
        Creates the perturbed density

        :param k:
            Response vector in matrix form in AO basis
        :param D:
            The density that is to be perturbed in AO basis
        :param S:
            Overlap matrix

        :return:
            [k,D]
        """

        return (np.linalg.multi_dot([k, S, D]) - np.linalg.multi_dot([D, S, k]))

    def mo2ao(self, mo, A):
        """
        Transform a matrix to atomic basis

        :param mo:
            molecular orbital coefficent matrix
        :param A:
            The matrix in MO basis that is the converted to AO basis

        :return:
            The matrix in AO basis
        """

        return np.linalg.multi_dot([mo, A, mo.T])

    def ao2mo(self, mo, A):
        """
        Transform a matrix to molecular basis

        :param mo:
            molecular orbital coefficent matrix
        :param A:
            The matrix in AO basis that is the converted to MO basis

        :return:
            The matrix in MO basis
        """

        return np.linalg.multi_dot([mo.T, A, mo])

    def commut(self, A, B):
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

    def x3_contract(self, k1, k2, X, D, nocc, norb):
        """
        Contracts the generalized dipole gradient tensor of rank 3 with two
        first-order response matrices. X[3]N1N2 = (1/2)[[k2,[k1,X]],D.T]

        :param: k1:
            First-order response matrix
        :param: k2:
            First-order response matrix
        :param X:
            Dipole intergral matrix
        :param D:
            Density matrix
        :param nocc:
            Number of occupied orbitals
        :param norb:
            Number of total orbtials

        :return:
            Returns a matrix
        """

        X3NxNy = self.commut(self.commut(k2, self.commut(k1, X)), D.T)
        X3NxNy_c = self.complex_lrmat2vec(X3NxNy, nocc, norb)
        return (1. / 2) * X3NxNy_c

    def x2_contract(self, k, X, D, nocc, norb):
        """
        Contracts the generalized dipole gradient tensor of rank 2 with a
        second-order response matrix. X[2]N1 = [[k1,X],D.T]

        :param: k:
            Respose vector in matrix representation
        :param X:
            Property operator in matrix represatiation
        :param D:
            Density matrix
        :param nocc:
            Number of occupied orbitals
        :param norb:
            Number of total orbtials

        :return:
            Returns a matrix
        """

        XNx = self.commut(self.commut(k, X), D.T)
        X2Nx_c = self.complex_lrmat2vec(XNx, nocc, norb)
        return X2Nx_c

    def a2_contract(self, k, A, D, nocc, norb):
        """
        Contracts the generalized dipole gradient tensor of rank 2 with a
        second-order response matrix. A[2]N1 = -(1 / 2)[[k1,X],D.T]

        # Note that the sign needs further investigation.

        :param: k:
            Respose vector in matrix representation
        :param A:
            Property operator in matrix represatiation
        :param D:
            Density matrix
        :param nocc:
            Number of occupied orbitals
        :param norb:
            Number of total orbtials

        :return:
            Returns a matrix
        """

        ANx = self.commut(self.commut(k.T, A), D.T)
        A2Nx_c = self.complex_lrmat2vec(ANx, nocc, norb)
        return -(1. / 2) * A2Nx_c

    def a3_contract(self, k1, k2, A, D, nocc, norb):
        """
        Contracts the generalized dipole gradient tensor of rank 3 with two
        first-order response matrices. A[3]N1N2 = -(1/6)[[k2,[k1,A]],D.T]

        :param: k1:
            First-order response matrix
        :param: k2:
            First-order response matrix
        :param A:
            A dipole intergral matrix
        :param D:
            Density matrix
        :param nocc:
            Number of occupied orbitals
        :param norb:
            Number of total orbtials

        :return:
            Returns a matrix
        """

        A3NxNy = self.commut(self.commut(k2.T, self.commut(k1.T, A)), D.T)
        A3NxNy_c = self.complex_lrmat2vec(A3NxNy, nocc, norb)
        return -(1. / 6) * A3NxNy_c

    def zi(self, kB, kC, kD, Fc, Fd, Fbc, Fcb, F0):
        """
        Returns a matrix used for the E[4] contraction

        :param kA:
            First-order response matrix
        :param kB:
            First-order response matrix
        :param Fa:
            First-order perturbed Fock matrix
        :param Fb:
            First-order perturbed Fock matrix
        :param F0:
            SCF Fock matrix

        :return:
            Returns a matrix
        """
        M1 = self.commut(kC, self.commut(kD, F0) + 3 * Fd)
        M2 = self.commut(kD, self.commut(kC, F0) + 3 * Fc)

        return (self.commut(kB, M1 + M2 + 3 * (Fbc + Fcb)))

    def anti_sym(self, vec):
        """
        Returns an antisymetrized vector

        :param vec:
            The vector to be anti-symetrized

        :return:
            An antisymetrized vector
        """

        if vec.ndim == 1:
            new_vec = np.zeros_like(vec)
            half_len = vec.shape[0] // 2
            new_vec[:half_len] = vec[:half_len]
            new_vec[half_len:] = -vec[half_len:]
            return new_vec

        return None

    def xi(self, kA, kB, Fa, Fb, F0):
        """
        Returns a matrix used for the E[3] contraction

        :param kA:
            First-order response matrix
        :param kB:
            First-order response matrix
        :param Fa:
            First-order perturbed Fock matrix
        :param Fb:
            First-order perturbed Fock matrix
        :param F0:
            SCF Fock matrix

        :return:
            Returns a matrix
        """

        return 0.5 * (self.commut(kA,
                                  self.commut(kB, F0) + 2 * Fb) +
                      self.commut(kB,
                                  self.commut(kA, F0) + 2 * Fa))

    def complex_lrmat2vec(self, mat, nocc, norb):
        """
        Converts complex matrix to vector.

        :param mat:
            The complex matrix.
        :param nocc:
            Number of occupied orbitals.
        :param norb:
            Number of orbitals.

        :return:
            The complex vector.
        """

        return (LinearSolver.lrmat2vec(mat.real, nocc, norb) +
                1j * LinearSolver.lrmat2vec(mat.imag, nocc, norb))
