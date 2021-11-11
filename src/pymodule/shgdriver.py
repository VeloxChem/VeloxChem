from mpi4py import MPI
from pathlib import Path
import numpy as np
import time
import sys
import re

from .veloxchemlib import (ElectricDipoleIntegralsDriver, mpi_master)
from .profiler import Profiler
from .outputstream import OutputStream
from .cppsolver import ComplexResponse
from .linearsolver import LinearSolver
from .nonlinearsolver import NonLinearSolver
from .distributedarray import DistributedArray
from .scffirstorderprop import ScfFirstOrderProperties
from .errorhandler import assert_msg_critical
from .checkpoint import (check_distributed_focks, read_distributed_focks,
                         write_distributed_focks)
from .inputparser import parse_input


class SHGDriver(NonLinearSolver):
    """
    Implements a quadratic response driver for SHG calculations

    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.

    Instance variables
        - is_converged: The flag for convergence.
        - eri_thresh: The electron repulsion integrals screening threshold.
        - qq_type: The electron repulsion integrals screening scheme.
        - batch_size: The batch size for computation of Fock matrices.
        - frequencies: The frequencies.
        - comp: The list of all the gamma tensor components
        - damping: The damping parameter.
        - lindep_thresh: The threshold for removing linear dependence in the
          trial vectors.
        - conv_thresh: The convergence threshold for the solver.
        - max_iter: The maximum number of solver iterations.
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
    """

    def __init__(self, comm=None, ostream=None):
        """
        Initializes the isotropic quadratic response driver for second harmonic
        generation (SHG).
        """

        if comm is None:
            comm = MPI.COMM_WORLD

        if ostream is None:
            ostream = OutputStream(sys.stdout)

        super().__init__(comm, ostream)

        self.is_converged = False

        # ERI settings
        self.eri_thresh = 1.0e-15
        self.qq_type = 'QQ_DEN'
        self.batch_size = None

        # cpp settings
        self.frequencies = (0,)
        self.comp = None
        self.damping = 0.004556335294880438
        self.lindep_thresh = 1.0e-10
        self.conv_thresh = 1.0e-4
        self.max_iter = 50
        self.a_component = 'xyz'
        self.b_component = 'xyz'
        self.a_operator = 'dipole'
        self.b_operator = 'dipole'
        self.c_operator = 'dipole'

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

    def update_settings(self, rsp_dict, method_dict=None):
        """
        Updates response and method settings in SHG driver

        :param rsp_dict:
            The dictionary of response dict.
        :param method_dict:
            The dictionary of method rsp_dict.
        """

        if method_dict is None:
            method_dict = {}

        rsp_keywords = {
            'frequencies': 'seq_range',
            'damping': 'float',
            'a_operator': 'str',
            'b_operator': 'str',
            'c_operator': 'str',
            'eri_thresh': 'float',
            'qq_type': 'str_upper',
            'batch_size': 'int',
            'max_iter': 'int',
            'conv_thresh': 'float',
            'lindep_thresh': 'float',
            'restart': 'bool',
            'checkpoint_file': 'str',
            'timing': 'bool',
            'profiling': 'bool',
            'memory_profiling': 'bool',
            'memory_tracing': 'bool',
        }

        parse_input(self, rsp_keywords, rsp_dict)

        if 'program_start_time' in rsp_dict:
            self.program_start_time = rsp_dict['program_start_time']
        if 'maximum_hours' in rsp_dict:
            self.maximum_hours = rsp_dict['maximum_hours']

        if 'xcfun' in method_dict:
            errmsg = 'ShgDriver: The \'xcfun\' keyword is not supported in SHG '
            errmsg += 'calculation.'
            if self.rank == mpi_master():
                assert_msg_critical(False, errmsg)

        if 'potfile' in method_dict:
            errmsg = 'ShgDriver: The \'potfile\' keyword is not supported in '
            errmsg += 'SHG calculation.'
            if self.rank == mpi_master():
                assert_msg_critical(False, errmsg)

        if 'electric_field' in method_dict:
            errmsg = 'ShgDriver: The \'electric field\' keyword is not '
            errmsg += 'supported in SHG calculation.'
            if self.rank == mpi_master():
                assert_msg_critical(False, errmsg)

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

        profiler = Profiler({
            'timing': False,
            'profiling': self.profiling,
            'memory_profiling': self.memory_profiling,
            'memory_tracing': self.memory_tracing,
        })

        if self.rank == mpi_master():
            self.print_header()

        start_time = time.time()

        # sanity check
        nalpha = molecule.number_of_alpha_electrons()
        nbeta = molecule.number_of_beta_electrons()
        assert_msg_critical(
            nalpha == nbeta,
            'SHG Driver: not implemented for unrestricted case')

        if self.rank == mpi_master():
            S = scf_tensors['S']
            da = scf_tensors['D'][0]
            mo = scf_tensors['C']
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
        a_rhs = linear_solver.get_complex_prop_grad(self.a_operator,
                                                    self.a_component, molecule,
                                                    ao_basis, scf_tensors)
        b_rhs = linear_solver.get_complex_prop_grad(self.b_operator,
                                                    self.b_component, molecule,
                                                    ao_basis, scf_tensors)

        if self.rank == mpi_master():
            inv_sqrt_2 = 1.0 / np.sqrt(2.0)
            b_rhs = list(b_rhs)
            for ind in range(len(b_rhs)):
                b_rhs[ind] *= inv_sqrt_2
            a_rhs = list(a_rhs)
            for ind in range(len(a_rhs)):
                a_rhs[ind] *= inv_sqrt_2

        # Storing the dipole integral matrices used for the X[3],X[2],A[3] and
        # A[2] contractions in MO basis
        wa = [sum(x) for x in zip(self.frequencies, self.frequencies)]

        freqpairs = [wl for wl in zip(self.frequencies, self.frequencies)]

        AB = {}

        if self.rank == mpi_master():
            A = {(op, w): v for op, v in zip('xyz', a_rhs) for w in wa}
            B = {(op, w): v for op, v in zip('xyz', b_rhs)
                 for w in self.frequencies}

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

        N_drv.update_settings({
            'damping': self.damping,
            'lindep_thresh': self.lindep_thresh,
            'conv_thresh': self.conv_thresh,
            'max_iter': self.max_iter,
            'eri_thresh': self.eri_thresh,
            'qq_type': self.qq_type,
        })

        N_drv.timing = self.timing
        N_drv.memory_profiling = self.memory_profiling
        N_drv.batch_size = self.batch_size
        N_drv.restart = self.restart
        N_drv.program_start_time = self.program_start_time
        N_drv.maximum_hours = self.maximum_hours
        if self.checkpoint_file is not None:
            N_drv.checkpoint_file = re.sub(r'\.h5$', r'', self.checkpoint_file)
            N_drv.checkpoint_file += '_quadb_1.h5'

        N_results = N_drv.compute(molecule, ao_basis, scf_tensors, AB)

        kX = N_results['kappas']
        Focks = N_results['focks']

        profiler.check_memory_usage('CPP')

        # Compute the isotropic parallel beta vector

        beta = self.compute_quad_components(Focks, freqpairs, X, d_a_mo, kX,
                                            self.comp, scf_tensors, molecule,
                                            ao_basis, profiler)

        # Compute dipole vector
        scf_prop = ScfFirstOrderProperties(self.comm, self.ostream)
        scf_prop.compute(molecule, ao_basis, scf_tensors)

        if self.rank == mpi_master():

            dip = list(scf_prop.get_property('dipole moment'))
            # Norm of dipole
            dip_norm = np.linalg.norm(scf_prop.get_property('dipole moment'))

            # Compute isotropic beta along molecular dipole

            self.ostream.print_blank()
            w_str = 'Electronic dipole moment: '
            self.ostream.print_header(w_str)
            self.ostream.print_header('=' * (len(w_str) + 2))
            title = '{:<9s} {:>12s} {:>20s} {:>21s}'.format(
                'Component', 'Frequency', 'Real', 'Imaginary')
            width = len(title)
            self.ostream.print_header(title.ljust(width))
            self.ostream.print_header(('-' * len(title)).ljust(width))
            self.print_component('μ_x', 0, dip[0], width)
            self.print_component('μ_y', 0, dip[1], width)
            self.print_component('μ_z', 0, dip[2], width)
            self.print_component('|μ|', 0, dip_norm, width)
            self.ostream.print_blank()
            w_str = 'Averaged first-order hyperpolarizability: '
            self.ostream.print_header(w_str)
            self.ostream.print_header('=' * (len(w_str) + 2))

            # TODO: look into beta_bar

            beta_bar = {}

            for key in beta.keys():
                betaa = 0
                for a in range(len(beta[key])):
                    betaa += 1 / dip_norm * dip[a] * beta[key][a]

                self.ostream.print_blank()

                self.ostream.print_header(title.ljust(width))
                self.ostream.print_header(('-' * len(title)).ljust(width))
                self.print_component('β_x', key, beta[key][0], width)
                self.print_component('β_y', key, beta[key][1], width)
                self.print_component('β_z', key, beta[key][2], width)
                self.print_component('β ||', key, 1 / 5 * betaa, width)
                beta_bar = {key: betaa}

            self.ostream.print_blank()
            valstr = '*** Time spent in quadratic response calculation: '
            valstr += '{:.2f} sec ***'.format(time.time() - start_time)
            self.ostream.print_header(valstr)
            self.ostream.print_blank()
            self.ostream.flush()

        profiler.end(self.ostream)

        self.is_converged = True

        return beta

    def compute_quad_components(self, Focks, freqpairs, X, d_a_mo, kX, track,
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
        :param kX:
            A dictonary containing all the response matricies
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
            S = scf_tensors['S']
            D0 = scf_tensors['D'][0]
            mo = scf_tensors['C']
            F0 = np.linalg.multi_dot([mo.T, scf_tensors['F'][0], mo])
            norb = mo.shape[1]
        else:
            S = None
            D0 = None
            mo = None
            F0 = None
            norb = None
        F0 = self.comm.bcast(F0, root=mpi_master())
        norb = self.comm.bcast(norb, root=mpi_master())

        nocc = molecule.number_of_alpha_electrons()

        # computing all compounded first-order densities
        if self.rank == mpi_master():
            density_list = self.get_densities(freqpairs, kX, S, D0, mo)
        else:
            density_list = None

        profiler.check_memory_usage('1st densities')

        #  computing the compounded first-order Fock matrices
        fock_dict = self.get_fock_dict(freqpairs, density_list, F0, mo,
                                       molecule, ao_basis)

        e3_dict = self.get_e3(freqpairs, kX, fock_dict, Focks, nocc, norb)

        beta = {}

        Na = {}
        Nb = {}

        if self.rank == mpi_master():
            for (wb, wc) in freqpairs:

                Na = {
                    'x': (LinearSolver.lrmat2vec(kX[('x', wb + wc)].real, nocc,
                                                 norb) +
                          1j * LinearSolver.lrmat2vec(kX[
                              ('x', wb + wc)].imag, nocc, norb)),
                    'y': (LinearSolver.lrmat2vec(kX[('y', wb + wc)].real, nocc,
                                                 norb) +
                          1j * LinearSolver.lrmat2vec(kX[
                              ('y', wb + wc)].imag, nocc, norb)),
                    'z': (LinearSolver.lrmat2vec(kX[('z', wb + wc)].real, nocc,
                                                 norb) +
                          1j * LinearSolver.lrmat2vec(kX[
                              ('z', wb + wc)].imag, nocc, norb))
                }

                Nb = {
                    'x': (
                        LinearSolver.lrmat2vec(kX[('x', wb)].real, nocc, norb) +
                        1j *
                        LinearSolver.lrmat2vec(kX[('x', wb)].imag, nocc, norb)),
                    'y': (
                        LinearSolver.lrmat2vec(kX[('y', wb)].real, nocc, norb) +
                        1j *
                        LinearSolver.lrmat2vec(kX[('y', wb)].imag, nocc, norb)),
                    'z': (
                        LinearSolver.lrmat2vec(kX[('z', wb)].real, nocc, norb) +
                        1j *
                        LinearSolver.lrmat2vec(kX[('z', wb)].imag, nocc, norb))
                }

                NaE3NbNc_x = np.dot(Na['x'].T, e3_dict[
                    ('sig_x', wb)]) + 2 * np.dot(Na['y'].T, e3_dict[
                        ('lam_xy', wb)]) + 2 * np.dot(Na['z'].T, e3_dict[
                            ('lam_xz', wb)])

                NaE3NbNc_y = np.dot(Na['y'].T, e3_dict[
                    ('sig_y', wb)]) + 2 * np.dot(Na['x'].T, e3_dict[
                        ('lam_xy', wb)]) + 2 * np.dot(Na['z'].T, e3_dict[
                            ('lam_yz', wb)])

                NaE3NbNc_z = np.dot(Na['z'].T, e3_dict[
                    ('sig_z', wb)]) + 2 * np.dot(Na['y'].T, e3_dict[
                        ('lam_yz', wb)]) + 2 * np.dot(Na['x'].T, e3_dict[
                            ('lam_xz', wb)])

                A2_x = 0
                A2_y = 0
                A2_z = 0

                X2_x = 0
                X2_y = 0
                X2_z = 0

                for eta in 'xyz':

                    # A2 contractions

                    A2_x -= np.dot(
                        Nb[eta].T,
                        self.a2_contract(kX[
                            (eta, wb)], X['x'], d_a_mo, nocc, norb)) + np.dot(
                                Nb[eta].T,
                                self.a2_contract(kX[(eta, wb)], X['x'], d_a_mo,
                                                 nocc, norb))
                    A2_x -= np.dot(
                        Nb['x'].T,
                        self.a2_contract(kX[
                            (eta, wb)], X[eta], d_a_mo, nocc, norb)) + np.dot(
                                Nb[eta].T,
                                self.a2_contract(kX[('x', wb)], X[eta], d_a_mo,
                                                 nocc, norb))
                    A2_x -= np.dot(
                        Nb[eta].T,
                        self.a2_contract(kX[
                            ('x', wb)], X[eta], d_a_mo, nocc, norb)) + np.dot(
                                Nb['x'].T,
                                self.a2_contract(kX[(eta, wb)], X[eta], d_a_mo,
                                                 nocc, norb))

                    A2_y -= np.dot(
                        Nb[eta].T,
                        self.a2_contract(kX[
                            (eta, wb)], X['y'], d_a_mo, nocc, norb)) + np.dot(
                                Nb[eta].T,
                                self.a2_contract(kX[(eta, wb)], X['y'], d_a_mo,
                                                 nocc, norb))
                    A2_y -= np.dot(
                        Nb['y'].T,
                        self.a2_contract(kX[
                            (eta, wb)], X[eta], d_a_mo, nocc, norb)) + np.dot(
                                Nb[eta].T,
                                self.a2_contract(kX[('y', wb)], X[eta], d_a_mo,
                                                 nocc, norb))
                    A2_y -= np.dot(
                        Nb[eta].T,
                        self.a2_contract(kX[
                            ('y', wb)], X[eta], d_a_mo, nocc, norb)) + np.dot(
                                Nb['y'].T,
                                self.a2_contract(kX[(eta, wb)], X[eta], d_a_mo,
                                                 nocc, norb))

                    A2_z -= np.dot(
                        Nb[eta].T,
                        self.a2_contract(kX[
                            (eta, wb)], X['z'], d_a_mo, nocc, norb)) + np.dot(
                                Nb[eta].T,
                                self.a2_contract(kX[(eta, wb)], X['z'], d_a_mo,
                                                 nocc, norb))
                    A2_z -= np.dot(
                        Nb['z'].T,
                        self.a2_contract(kX[
                            (eta, wb)], X[eta], d_a_mo, nocc, norb)) + np.dot(
                                Nb[eta].T,
                                self.a2_contract(kX[('z', wb)], X[eta], d_a_mo,
                                                 nocc, norb))
                    A2_z -= np.dot(
                        Nb[eta].T,
                        self.a2_contract(kX[
                            ('z', wb)], X[eta], d_a_mo, nocc, norb)) + np.dot(
                                Nb['z'].T,
                                self.a2_contract(kX[(eta, wb)], X[eta], d_a_mo,
                                                 nocc, norb))

                    # X2 contractions

                    X2_x -= 2 * np.dot(
                        Na['x'].T,
                        self.x2_contract(kX[(eta, wb)], X[eta], d_a_mo, nocc,
                                         norb))
                    X2_x -= 2 * np.dot(
                        Na[eta].T,
                        self.x2_contract(kX[(eta, wb)], X['x'], d_a_mo, nocc,
                                         norb))
                    X2_x -= 2 * np.dot(
                        Na[eta].T,
                        self.x2_contract(kX[('x', wb)], X[eta], d_a_mo, nocc,
                                         norb))

                    X2_y -= 2 * np.dot(
                        Na['y'].T,
                        self.x2_contract(kX[(eta, wb)], X[eta], d_a_mo, nocc,
                                         norb))
                    X2_y -= 2 * np.dot(
                        Na[eta].T,
                        self.x2_contract(kX[(eta, wb)], X['y'], d_a_mo, nocc,
                                         norb))
                    X2_y -= 2 * np.dot(
                        Na[eta].T,
                        self.x2_contract(kX[('y', wb)], X[eta], d_a_mo, nocc,
                                         norb))

                    X2_z -= 2 * np.dot(
                        Na['z'].T,
                        self.x2_contract(kX[(eta, wb)], X[eta], d_a_mo, nocc,
                                         norb))
                    X2_z -= 2 * np.dot(
                        Na[eta].T,
                        self.x2_contract(kX[(eta, wb)], X['z'], d_a_mo, nocc,
                                         norb))
                    X2_z -= 2 * np.dot(
                        Na[eta].T,
                        self.x2_contract(kX[('z', wb)], X[eta], d_a_mo, nocc,
                                         norb))

                beta.update({
                    wb: (NaE3NbNc_x + A2_x + X2_x, NaE3NbNc_y + A2_y + X2_y,
                         NaE3NbNc_z + A2_z + X2_z)
                })

        profiler.check_memory_usage('End of SHG')

        return beta

    def get_densities(self, freqpairs, kX, S, D0, mo):
        """

        :param freqpairs:
            A list of the frequencies
        :param kX:
            A dictonary with all the first-order response matrices
        :param S:
            The overlap matrix
        :param D0:
            The SCF density matrix in AO basis
        :param mo:
            A matrix containing the MO coefficents

        :return:
            A list of tranformed compounded densities
        """

        density_list = []

        for (wb, wc) in freqpairs:

            # convert response matrix to ao basis #

            k_x = self.mo2ao(mo, kX[('x', wb)])
            k_y = self.mo2ao(mo, kX[('y', wb)])
            k_z = self.mo2ao(mo, kX[('z', wb)])

            # create the first order single indexed densiteies #

            D_x = self.transform_dens(k_x, D0, S)
            D_y = self.transform_dens(k_y, D0, S)
            D_z = self.transform_dens(k_z, D0, S)

            # create the first order two indexed densities #

            D_sig_x = 4 * self.transform_dens(k_x, D_x, S) + 2 * (
                self.transform_dens(k_x, D_x, S) + self.transform_dens(
                    k_y, D_y, S) + self.transform_dens(k_z, D_z, S))
            D_sig_y = 4 * self.transform_dens(k_y, D_y, S) + 2 * (
                self.transform_dens(k_x, D_x, S) + self.transform_dens(
                    k_y, D_y, S) + self.transform_dens(k_z, D_z, S))
            D_sig_z = 4 * self.transform_dens(k_z, D_z, S) + 2 * (
                self.transform_dens(k_x, D_x, S) + self.transform_dens(
                    k_y, D_y, S) + self.transform_dens(k_z, D_z, S))

            D_lam_xy = self.transform_dens(k_x, D_y, S) + self.transform_dens(
                k_y, D_x, S)
            D_lam_xz = self.transform_dens(k_x, D_z, S) + self.transform_dens(
                k_z, D_x, S)
            D_lam_yz = self.transform_dens(k_y, D_z, S) + self.transform_dens(
                k_z, D_y, S)

            density_list.append(D_sig_x)
            density_list.append(D_sig_y)
            density_list.append(D_sig_z)
            density_list.append(D_lam_xy)
            density_list.append(D_lam_xz)
            density_list.append(D_lam_yz)

        return density_list

    def get_fock_dict(self, wi, density_list, F0, mo, molecule, ao_basis):
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
            A dictonary of compounded first-order Fock-matrices
        """

        if self.rank == mpi_master():
            self.print_fock_header()

        ww = []
        for (wb, wc) in wi:
            ww.append(wb)

        keys = [
            'F_sig_x', 'F_sig_y', 'F_sig_z', 'F_lam_xy', 'F_lam_xz', 'F_lam_yz'
        ]

        if self.checkpoint_file is not None:
            fock_file = str(
                Path(self.checkpoint_file).with_suffix('.shg_fock.h5'))
        else:
            fock_file = None

        if self.restart:
            if self.rank == mpi_master():
                self.restart = check_distributed_focks(fock_file, keys, ww)
            self.restart = self.comm.bcast(self.restart, mpi_master())

        if self.restart:
            focks = read_distributed_focks(fock_file, keys, ww, self.comm,
                                           self.ostream)
            focks['F0'] = F0
            return focks

        time_start_fock = time.time()
        dist_focks = self.comp_nlr_fock(mo, density_list, molecule, ao_basis,
                                        'real_and_imag')
        time_end_fock = time.time()

        total_time_fock = time_end_fock - time_start_fock
        self.print_fock_time(total_time_fock)

        focks = {'F0': F0}
        for key in keys:
            focks[key] = {}

        fock_index = 0

        for (wb, wc) in wi:
            for key in keys:
                focks[key][wb] = DistributedArray(dist_focks.data[:,
                                                                  fock_index],
                                                  self.comm,
                                                  distribute=False)
                fock_index += 1

        write_distributed_focks(fock_file, focks, keys, ww, self.comm,
                                self.ostream)

        return focks

    def get_e3(self, wi, kX, fo, fo2, nocc, norb):
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
                fo2[('x', wb)].data, fo2[('y', wb)].data, fo2[('z', wb)].data,
                fo['F_sig_x'][wb].data, fo['F_sig_y'][wb].data,
                fo['F_sig_z'][wb].data, fo['F_lam_xy'][wb].data,
                fo['F_lam_xz'][wb].data, fo['F_lam_yz'][wb].data
            ]).T.copy()

            vec_pack = self.collect_vectors_in_columns(vec_pack)

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

            k_x = kX[('x', wb)].T

            k_y = kX[('y', wb)].T

            k_z = kX[('z', wb)].T

            # Make all Xi terms

            xi_sig_x = 3 * self.xi(k_x, k_x, F_x, F_x, F0_a) + self.xi(
                k_y, k_y, F_y, F_y, F0_a) + self.xi(k_z, k_z, F_z, F_z, F0_a)

            xi_sig_y = 3 * self.xi(k_y, k_y, F_y, F_y, F0_a) + self.xi(
                k_z, k_z, F_z, F_z, F0_a) + self.xi(k_z, k_z, F_z, F_z, F0_a)

            xi_sig_z = 3 * self.xi(k_z, k_z, F_z, F_z, F0_a) + self.xi(
                k_y, k_y, F_y, F_y, F0_a) + self.xi(k_x, k_x, F_x, F_x, F0_a)

            xi_lam_xy = self.xi(k_x, k_y, F_x, F_y, F0_a)

            xi_lam_xz = self.xi(k_x, k_z, F_x, F_z, F0_a)

            xi_lam_yz = self.xi(k_y, k_z, F_y, F_z, F0_a)

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

    def print_header(self):
        """
        Prints SHG setup header to output stream.
        """

        self.ostream.print_blank()

        title = 'SHG Driver Setup'
        self.ostream.print_header(title)
        self.ostream.print_header('=' * (len(title) + 2))
        self.ostream.print_blank()

        width = 50

        cur_str = 'ERI Screening Threshold         : {:.1e}'.format(
            self.eri_thresh)
        self.ostream.print_header(cur_str.ljust(width))
        cur_str = 'Convergance Threshold           : {:.1e}'.format(
            self.conv_thresh)
        self.ostream.print_header(cur_str.ljust(width))
        cur_str = 'Max. Number of Iterations       : {:d}'.format(self.max_iter)
        self.ostream.print_header(cur_str.ljust(width))
        cur_str = 'Damping Parameter               : {:.6e}'.format(
            self.damping)
        self.ostream.print_header(cur_str.ljust(width))

        self.ostream.print_blank()
        self.ostream.flush()

    def print_component(self, label, freq, value, width):
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

        w_str = '{:<9s} {:12.4f} {:20.8f} {:20.8f}j'.format(
            label, freq, value.real, value.imag)
        self.ostream.print_header(w_str.ljust(width))
