from mpi4py import MPI
import numpy as np
import time
import re
import time as tm

from .veloxchemlib import ElectronRepulsionIntegralsDriver
from .veloxchemlib import ElectricDipoleIntegralsDriver
from .veloxchemlib import denmat
from .veloxchemlib import fockmat
from .veloxchemlib import mpi_master
from .qqscheme import get_qq_scheme
from .profiler import Profiler
from .cppsolver import ComplexResponse
from .linearsolver import LinearSolver
from .aofockmatrix import AOFockMatrix
from .aodensitymatrix import AODensityMatrix
from .distributedarray import DistributedArray
from .inputparser import parse_seq_range
from .errorhandler import assert_msg_critical
from .batchsize import get_batch_size
from .batchsize import get_number_of_batches

from .veloxchemlib import XCFunctional
from .veloxchemlib import XCIntegrator
from .veloxchemlib import parse_xc_func
from .veloxchemlib import GridDriver

from .checkpoint import check_distributed_focks
from .checkpoint import read_distributed_focks
from .checkpoint import write_distributed_focks


from .inputparser import parse_input

class CubicResponseDriver:
    """
    Implements a general quadratic response driver 

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

    def __init__(self, comm, ostream):
        """
        Initializes the quadratic response driver
        """

        self.is_converged = False

        # ERI settings
        self.eri_thresh = 1.0e-15
        self.qq_type = 'QQ_DEN'
        self.batch_size = None

        # cpp settings
        self.b_frequencies = (0,)
        self.c_frequencies = (0,)
        self.d_frequencies = (0,)
        self.comp = None
        self.damping = 0.004556335294880438
        self.lindep_thresh = 1.0e-10
        self.conv_thresh = 1.0e-4
        self.max_iter = 50
        self.a_component = 'z'
        self.b_component = 'z'
        self.c_component = 'z'
        self.d_component = 'z'

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
        Updates response and method settings

        :param rsp_dict:
            The dictionary of response dict.
        :param method_dict:
            The dictionary of method rsp_dict.
        """

        if method_dict is None:
            method_dict = {}


        rsp_keywords = {
            'b_frequencies': 'seq_range',
            'c_frequencies': 'seq_range',
            'd_frequencies': 'seq_range',
            'a_component': 'str',
            'b_component': 'str',
            'c_component': 'str',
            'd_component': 'str',
            'damping': 'float',
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
            errmsg = 'CrfDriver: The \'xcfun\' keyword is not supported in Crf '
            errmsg += 'calculation.'
            if self.rank == mpi_master():
                assert_msg_critical(False, errmsg)

        if 'potfile' in method_dict:
            errmsg = 'CrfDriver: The \'potfile\' keyword is not supported in '
            errmsg += 'Crf calculation.'
            if self.rank == mpi_master():
                assert_msg_critical(False, errmsg)

        if 'electric_field' in method_dict:
            errmsg = 'CrfDriver: The \'electric field\' keyword is not '
            errmsg += 'supported in Crf calculation.'
            if self.rank == mpi_master():
                assert_msg_critical(False, errmsg)


    def compute(self, molecule, ao_basis, scf_tensors):
        """
        Computes a quadratic response function 

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
        assert_msg_critical(nalpha == nbeta,
                            'Cubic response driver: not implemented for unrestricted case')

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

        operator = 'dipole'

        linear_solver = LinearSolver(self.comm, self.ostream)

        a_rhs = linear_solver.get_complex_prop_grad(operator, self.a_component,
                                                    molecule, ao_basis,
                                                    scf_tensors)

        b_rhs = linear_solver.get_complex_prop_grad(operator, self.b_component,
                                                    molecule, ao_basis,
                                                    scf_tensors)

        c_rhs = linear_solver.get_complex_prop_grad(operator, self.c_component,
                                                    molecule, ao_basis,
                                                    scf_tensors)

        d_rhs = linear_solver.get_complex_prop_grad(operator, self.d_component,
                                                    molecule, ao_basis,
                                                    scf_tensors)

        if self.rank == mpi_master():
            inv_sqrt_2 = 1.0 / np.sqrt(2.0)
            b_rhs = list(b_rhs)
            for ind in range(len(b_rhs)):
                b_rhs[ind] *= inv_sqrt_2
            a_rhs = list(a_rhs)
            for ind in range(len(a_rhs)):
                a_rhs[ind] *= inv_sqrt_2
            c_rhs = list(c_rhs)
            for ind in range(len(c_rhs)):
                c_rhs[ind] *= inv_sqrt_2
            d_rhs = list(d_rhs)
            for ind in range(len(d_rhs)):
                d_rhs[ind] *= inv_sqrt_2

        # Storing the dipole integral matrices used for the X[3],X[2],A[3] and
        # A[2] contractions in MO basis
        wa = [sum(x) for x in zip(self.b_frequencies, self.c_frequencies,self.d_frequencies)]

        freqpairs = [wl for wl in zip(self.b_frequencies, self.c_frequencies,self.d_frequencies)]

        if self.rank == mpi_master():
            A = {(op, w): v for op, v in zip('A', a_rhs) for w in wa}
            B = {(op, w): v for op, v in zip('B', b_rhs) for w in self.b_frequencies}
            C = {(op, w): v for op, v in zip('C', c_rhs) for w in self.c_frequencies}
            D = {(op, w): v for op, v in zip('D', d_rhs) for w in self.d_frequencies}

            X = {
                'x': 2 * self.ao2mo(mo, dipole_mats.x_to_numpy()),
                'y': 2 * self.ao2mo(mo, dipole_mats.y_to_numpy()),
                'z': 2 * self.ao2mo(mo, dipole_mats.z_to_numpy())
            }

        else:
            v1 = None
            X = None
            self.comp = None

        ABCD = {}
        ABCD.update(A)
        ABCD.update(B)
        ABCD.update(C)
        ABCD.update(D)

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
            N_drv.checkpoint_file += '_quada_1.h5'

        N_results = N_drv.compute(molecule, ao_basis, scf_tensors, ABCD)

        kX = {}
        Focks = {}

        kX = N_results['kappas']

        Focks = N_results['focks']

        profiler.check_memory_usage('CPP')

        cubic_dict = self.compute_cubic_components(Focks, freqpairs, X,
                                               d_a_mo, kX, self.comp,
                                               scf_tensors, molecule, ao_basis,
                                               profiler)

        valstr = '*** Time spent in cubic response calculation: {:.2f} sec ***'.format(
            time.time() - start_time)
        self.ostream.print_header(valstr)
        self.ostream.print_blank()
        self.ostream.flush()

        profiler.end(self.ostream)

        self.is_converged = True

        return cubic_dict

    def compute_cubic_components(self, Focks, freqpairs, X, d_a_mo, kX, track,
                               scf_tensors, molecule, ao_basis, profiler):
        """
        Computes all the relevent terms to compute a general quadratic response function

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
        fock_dict = self.get_fock_dict(freqpairs, density_list, F0, mo, molecule,
                                       ao_basis)

        e4_dict,s4_dict,r4_dict = self.get_esr4(freqpairs, kX, fock_dict,Focks, nocc,norb,d_a_mo)


        k_xy,f_xy = self.get_nxy(freqpairs, kX, fock_dict,Focks, nocc, norb,d_a_mo,X,molecule,ao_basis,scf_tensors)

        if self.rank == mpi_master():
            density_list_ii = self.get_densities_ii(freqpairs, kX,k_xy, S, D0, mo)
        else:
            density_list_ii = None

        fock_dict_ii = self.get_fock_dict_ii(freqpairs, density_list_ii, F0, mo, molecule,
                                       ao_basis)

        e3_dict = self.get_e3(freqpairs, kX, k_xy, Focks,fock_dict_ii,f_xy, nocc, norb)

        A = X[self.a_component]
        B = X[self.b_component]
        C = X[self.c_component]
        D = X[self.d_component]
        
        result = {}

        for (wb,wc,wd) in freqpairs:

            Na = (LinearSolver.lrmat2vec(kX[('A',wb+wc+wd)].real, nocc, norb) +
                1j * LinearSolver.lrmat2vec(kX[('A',wb+wc+wd)].imag, nocc, norb))

            Nb = (LinearSolver.lrmat2vec(kX[('B',wb)].real, nocc, norb) +
                1j * LinearSolver.lrmat2vec(kX[('B',wb)].imag, nocc, norb))

            Nc = (LinearSolver.lrmat2vec(kX[('C',wc)].real, nocc, norb) +
                1j * LinearSolver.lrmat2vec(kX[('C',wc)].imag, nocc, norb))
            
            Nd = (LinearSolver.lrmat2vec(kX[('D',wd)].real, nocc, norb) +
                1j * LinearSolver.lrmat2vec(kX[('D',wd)].imag, nocc, norb))

            Nbd = (LinearSolver.lrmat2vec(k_xy[('BD', wb,wd),wb+wd].real, nocc, norb) +
                   1j * LinearSolver.lrmat2vec(k_xy[('BD', wb,wd),wb+wd].imag, nocc, norb))

            Nbc = (LinearSolver.lrmat2vec(k_xy[('BC', wb,wc),wb+wc].real, nocc, norb) +
                   1j * LinearSolver.lrmat2vec(k_xy[('BC', wb,wc),wb+wc].imag, nocc, norb))

            Ncd = (LinearSolver.lrmat2vec(k_xy[('CD', wc,wd),wc+wd].real, nocc, norb) +
                   1j * LinearSolver.lrmat2vec(k_xy[('CD', wc,wd),wc+wd].imag, nocc, norb))

            NaE4NbNcNd = np.dot(Na, e4_dict[wb])

            NaS4NbNcNd = np.dot(Na, s4_dict[wb])

            NaR4NbNcNd =  r4_dict[wb]

            NaE3NbNcd = np.dot(Na,e3_dict[(('E3NbNcd'),(wb,wc,wd))])
            NaE3NcNbd = np.dot(Na,e3_dict[(('E3NcNbd'),(wb,wc,wd))])
            NaE3NdNbc = np.dot(Na,e3_dict[(('E3NdNbc'),(wb,wc,wd))])

            # X3 terms
            NaB3NcNd = np.dot(Na.T, self.x3_contract(kX[('C',wc)], kX[('D',wd)], B, d_a_mo, nocc, norb))
            NaB3NdNc = np.dot(Na.T, self.x3_contract(kX[('D',wd)], kX[('C',wc)], B, d_a_mo, nocc, norb))

            NaC3NbNd = np.dot(Na.T, self.x3_contract(kX[('B',wb)], kX[('D',wd)], C, d_a_mo, nocc, norb))
            NaC3NdNb = np.dot(Na.T, self.x3_contract(kX[('D',wd)], kX[('B',wb)], C, d_a_mo, nocc, norb))

            NaD3NbNc = np.dot(Na.T, self.x3_contract(kX[('B',wb)], kX[('C',wc)], D, d_a_mo, nocc, norb))
            NaD3NcNb = np.dot(Na.T, self.x3_contract(kX[('C',wc)], kX[('B',wb)], D, d_a_mo, nocc, norb))

            # X2 contraction
            NaB2Ncd = np.dot(Na.T, self.x2_contract(k_xy[('CD', wc,wd),wc+wd], B, d_a_mo, nocc, norb))
            NaC2Nbd = np.dot(Na.T, self.x2_contract(k_xy[('BD', wb,wd),wb+wd], C, d_a_mo, nocc, norb))
            NaD2Nbc = np.dot(Na.T, self.x2_contract(k_xy[('BC', wb,wc),wb+wc], D, d_a_mo, nocc, norb))

            # A3 contraction  
            NdA3NbNc = np.dot(self.a3_contract(kX[('B',wb)], kX[('C',wc)], A, d_a_mo, nocc, norb), Nd)
            NdA3NcNb = np.dot(self.a3_contract(kX[('C',wc)], kX[('B',wb)], A, d_a_mo, nocc, norb), Nd)

            NbA3NcNd = np.dot(self.a3_contract(kX[('C',wc)], kX[('D',wd)], A, d_a_mo, nocc, norb), Nb)
            NbA3NdNc = np.dot(self.a3_contract(kX[('D',wd)], kX[('C',wc)], A, d_a_mo, nocc, norb), Nb)

            NcA3NbNd = np.dot(self.a3_contract(kX[('B',wb)], kX[('D',wd)], A, d_a_mo, nocc, norb), Nc)
            NcA3NdNb = np.dot(self.a3_contract(kX[('D',wd)], kX[('B',wb)], A, d_a_mo, nocc, norb), Nc)

            # A2 contraction 
            NbA2Ncd = np.dot(self.a2_contract(kX[('B',wb)], A, d_a_mo, nocc, norb), Ncd)
            NcdA2Nb = np.dot(self.a2_contract(k_xy[('CD', wc,wd),wc+wd], A, d_a_mo, nocc, norb), Nb)

            NcA2Nbd = np.dot(self.a2_contract(kX[('C',wc)], A, d_a_mo, nocc, norb), Nbd)
            NbdA2Nc = np.dot(self.a2_contract(k_xy[('BD', wb,wd),wb+wd], A, d_a_mo, nocc, norb), Nc)

            NdA2Nbc = np.dot(self.a2_contract(kX[('D',wd)], A, d_a_mo, nocc, norb), Nbc)
            NbcA2Nd = np.dot(self.a2_contract(k_xy[('BC', wb,wc),wb+wc], A, d_a_mo, nocc, norb), Nd)

            # Cubic response function
            Gamma = -(NaE4NbNcNd-NaS4NbNcNd-NaR4NbNcNd) -(NaE3NbNcd+NaE3NcNbd+NaE3NdNbc) 
            Gamma +=  NaB3NcNd+NaB3NdNc +NaC3NbNd+NaC3NdNb +NaD3NbNc +NaD3NcNb
            Gamma += -(NdA3NbNc + NdA3NcNb + NbA3NcNd + NbA3NdNc + NcA3NbNd +NcA3NdNb)
            Gamma += NaB2Ncd + NaC2Nbd + NaD2Nbc
            Gamma += NbA2Ncd + NcdA2Nb + NcA2Nbd + NbdA2Nc + NdA2Nbc + NbcA2Nd

            self.ostream.print_blank()
            w_str = 'Cubic response function: ' + '<< ' + str(self.a_component) +';' + str(self.b_component) + ',' + str(self.c_component) +  ',' + str(self.d_component) + ' >> ' + ' ('+str(wb) + ' ,' + str(wc) + ',' + str(wd) + ')'
            self.ostream.print_header(w_str)
            self.ostream.print_header('=' * (len(w_str) + 2))
            self.ostream.print_blank()
            title = '{:<9s}  {:>20s} {:>21s}'.format(
                'Component', 'Real', 'Imaginary')
            width = len(title)
            self.ostream.print_header(title.ljust(width))
            self.ostream.print_header(('-' * len(title)).ljust(width))
            self.print_component('E3', -(NaE3NbNcd+NaE3NcNbd+NaE3NdNbc), width)
            self.print_component('T4', -(NaE4NbNcNd-NaS4NbNcNd-NaR4NbNcNd), width)
            self.print_component('X2', NaB2Ncd + NaC2Nbd + NaD2Nbc, width)
            self.print_component('X3', NaB3NcNd+NaB3NdNc +NaC3NbNd+NaC3NdNb +NaD3NbNc +NaD3NcNb , width)
            self.print_component('A2', NbA2Ncd + NcdA2Nb + NcA2Nbd + NbdA2Nc + NdA2Nbc + NbcA2Nd , width)
            self.print_component('A3', -(NdA3NbNc + NdA3NcNb + NbA3NcNd + NbA3NdNc + NcA3NbNd +NcA3NdNb)   , width)
            self.print_component('γ', Gamma , width)
            self.ostream.print_blank()
            result.update({('E3',wb,wc,wd): -(NaE3NbNcd+NaE3NcNbd+NaE3NdNbc)})
            result.update({('T4',wb,wc,wd): -(NaE4NbNcNd-NaS4NbNcNd-NaR4NbNcNd)})
            result.update({('X3',wb,wc,wd): (NaB3NcNd+NaB3NdNc +NaC3NbNd+NaC3NdNb +NaD3NbNc +NaD3NcNb)})
            result.update({('X2',wb,wc,wd): NaB2Ncd + NaC2Nbd + NaD2Nbc})
            result.update({('A3',wb,wc,wd): -(NdA3NbNc + NdA3NcNb + NbA3NcNd + NbA3NdNc + NcA3NbNd +NcA3NdNb)})
            result.update({('A2',wb,wc,wd): NbA2Ncd + NcdA2Nb + NcA2Nbd + NbdA2Nc + NdA2Nbc + NbcA2Nd})
        

        profiler.check_memory_usage('End of QRF')

        return result

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
        A3NxNy_c = (LinearSolver.lrmat2vec(A3NxNy.real, nocc, norb) +
                    1j * LinearSolver.lrmat2vec(A3NxNy.imag, nocc, norb))
        return -(1. / 6) * A3NxNy_c


    def get_e3(self, wi, kX,k_xy, fo,fo2,fo3, nocc, norb):
        """
        Contracts E[3] for CRF

        :param wi:
            A list of freqs
        :param kX:
            A dict of the single index response matricies
        :param kXY:
            A dict of the two index response matrices
        :param fo:
            A dictonary of transformed Fock matricies from fock_dict
        :param fo2:
            A dictonarty of transfromed Fock matricies from fock_dict_two
        :param nocc:
            The number of occupied orbitals
        :param norb:
            The total number of orbitals

        :return:
            A dictonary of compounded E[3] tensors for the isotropic cubic
            response function for QRF
        """

        e3vec = {}

        for (wb,wc,wd) in wi:

            vec_pack = np.array([
                fo[('B', wb)].data,
                fo[('C', wc)].data,
                fo[('D', wd)].data,
                fo2['Fb_cd'][(wb,wc,wd)].data,
                fo2['Fcd_b'][(wb,wc,wd)].data,
                fo2['Fc_bd'][(wb,wc,wd)].data,
                fo2['Fbd_c'][(wb,wc,wd)].data,
                fo2['Fd_bc'][(wb,wc,wd)].data,
                fo2['Fbc_d'][(wb,wc,wd)].data,
                fo3[(('BC',wb,wc),wb+wc)].data,
                fo3[(('BD',wb,wd),wb+wd)].data,
                fo3[(('CD',wc,wd),wc+wd)].data,
            ]).T.copy()

            vec_pack = self.collect_vectors_in_columns(vec_pack)

            if self.rank != mpi_master():
                continue

            vec_pack = vec_pack.T.copy().reshape(-1, norb, norb)

            (fb, fc,fd,fb_cd,fcd_b,fc_bd,fbd_c,fd_bc,fbc_d,fbc,fbd,fcd) = vec_pack

            fb = np.conjugate(fb).T
            fc = np.conjugate(fc).T
            fd = np.conjugate(fd).T
            fbc = np.conjugate(fbc).T
            fcd = np.conjugate(fcd).T
            fbd = np.conjugate(fbd).T

            F0_a = fo2['F0']

            # E3NbNcd

            kb = kX[('B', wb)].T
            kcd = k_xy[('CD', wc,wd),wc+wd].T

            xi = self.xi(kb, kcd, fb, fcd, F0_a)

            e3fock = xi.T + (0.5 * fb_cd + 0.5*fcd_b).T
            e3vec[('E3NbNcd',(wb,wc,wd))] = self.anti_sym(-2 * LinearSolver.lrmat2vec(e3fock, nocc, norb))

            # E3NcNbd

            kc = kX[('C', wc)].T
            kbd = k_xy[('BD', wb,wd),wb+wd].T

            xi = self.xi(kc, kbd, fc, fbd, F0_a)

            e3fock = xi.T + (0.5 * fc_bd + 0.5*fbd_c).T
            e3vec[('E3NcNbd',(wb,wc,wd))] = self.anti_sym(-2 * LinearSolver.lrmat2vec(e3fock, nocc, norb))

            # E3NdNbc

            kd = kX[('D', wd)].T
            kbc = k_xy[('BC', wb,wc),wb+wc].T

            xi = self.xi(kd, kbc, fd, fbc, F0_a)

            e3fock = xi.T + (0.5 * fd_bc + 0.5*fbc_d).T
            e3vec[('E3NdNbc',(wb,wc,wd))] = self.anti_sym(-2 * LinearSolver.lrmat2vec(e3fock, nocc, norb))

        return e3vec

    def get_densities(self, freqpairs, kX, S, D0, mo):
        """
        Computes the  densities needed for the  Fock
        matrics 

        :param wi:
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

        for (wb,wc,wd) in freqpairs:

                # convert response matrix to ao basis #

                kb = self.mo2ao(mo, kX[('B', wb)])
                kc = self.mo2ao(mo, kX[('C', wc)])
                kd = self.mo2ao(mo, kX[('D', wd)])

                # create the first order single indexed densiteies #

                Db = self.transform_dens(kb, D0, S)
                Dc = self.transform_dens(kc, D0, S)
                Dd = self.transform_dens(kd, D0, S)

                # create the first order two indexed densities #

                Dbc = self.transform_dens(kb, Dc, S)
                Dcb = self.transform_dens(kc, Db, S)

                Dbd = self.transform_dens(kb, Dd, S)
                Ddb = self.transform_dens(kd, Db, S)

                Ddc = self.transform_dens(kd, Dc, S)
                Dcd = self.transform_dens(kc, Dd, S)

                # create the first order three indexed densities #

                Dbcd = self.transform_dens(kb, Dcd, S)
                Dbdc = self.transform_dens(kb, Ddc, S)

                Dcbd = self.transform_dens(kc, Dbd, S)
                Dcdb = self.transform_dens(kc, Ddb, S)

                Ddbc = self.transform_dens(kd, Dbc, S)
                Ddcb = self.transform_dens(kd, Dcb, S)

                density_list.append(Dbc)
                density_list.append(Dcb)
                density_list.append(Dbd)
                density_list.append(Ddb)
                density_list.append(Ddc)
                density_list.append(Dcd)

                density_list.append(Dbcd)
                density_list.append(Dbdc)
                density_list.append(Dcbd)
                density_list.append(Dcdb)
                density_list.append(Ddbc)
                density_list.append(Ddcb)

        return density_list


    def get_densities_ii(self, freqpairs, kX, k_xy, S, D0, mo):
        """
        Computes the  densities needed for the  Fock
        matrics 

        :param wi:
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

        for (wb,wc,wd) in freqpairs:

                # convert response matrix to ao basis #

                kb = self.mo2ao(mo, kX[('B', wb)])
                kc = self.mo2ao(mo, kX[('C', wc)])
                kd = self.mo2ao(mo, kX[('D', wd)])

                kbc = self.mo2ao(mo, k_xy[(('BC', wb,wc),wb+wc)])
                kbd = self.mo2ao(mo, k_xy[(('BD', wb,wd),wb+wd)])
                kcd = self.mo2ao(mo, k_xy[(('CD', wc,wd),wc+wd)])

                # create the first order single indexed densiteies #

                Db = self.transform_dens(kb, D0, S)
                Dc = self.transform_dens(kc, D0, S)
                Dd = self.transform_dens(kd, D0, S)

                # create the second-order two indexed densities #

                Dbc = self.transform_dens(kbc, D0, S)
                Dbd = self.transform_dens(kbd, D0, S)
                Dcd = self.transform_dens(kcd, D0, S)

                # create the second-order three indexed densities #

                Db_cd = self.transform_dens(kb, Dcd, S)
                Dcd_b = self.transform_dens(kcd, Db, S)

                Dc_bd = self.transform_dens(kc, Dbd, S)
                Dbd_c = self.transform_dens(kbd, Dc, S)

                Dd_bc = self.transform_dens(kd, Dbc, S)
                Dbc_d = self.transform_dens(kbc, Dd, S)

                density_list.append(Db_cd)
                density_list.append(Dcd_b)
                density_list.append(Dc_bd)
                density_list.append(Dbd_c)
                density_list.append(Dd_bc)
                density_list.append(Dbc_d)

        return density_list


    def get_fock_dict(self, wi, density_list, F0, mo, molecule, ao_basis):
        """
        Computes the Fock matrices for a quadratic response function

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

        keys = [
            'Fbc',
            'Fcb',
            'Fbd',
            'Fdb',
            'Fdc',
            'Fcd',
            'Fbcd',
            'Fbdc',
            'Fcbd',
            'Fcdb',
            'Fdbc',
            'Fdcb',
        ]

        if self.checkpoint_file is not None:
            fock_file = str(
                Path(self.checkpoint_file).with_suffix('.tpa_fock_1_full.h5'))
        else:
            fock_file = None

        if self.restart:
            if self.rank == mpi_master():
                self.restart = check_distributed_focks(fock_file, keys, wi)
            self.restart = self.comm.bcast(self.restart, mpi_master())

        if self.restart:
            focks = read_distributed_focks(fock_file, keys, wi, self.comm,
                                           self.ostream)
            focks['F0'] = F0
            return focks

        time_start_fock = time.time()
        dist_focks = self.get_fock(mo, density_list, molecule, ao_basis,
                                     'real_and_imag')
        time_end_fock = time.time()

        total_time_fock = time_end_fock - time_start_fock
        self.print_fock_time(total_time_fock)
        
        focks = {'F0': F0}
        for key in keys:
            focks[key] = {}

        fock_index = 0
        for (wb,wc,wd) in wi:
            for key in keys:
                focks[key][wb] = DistributedArray(dist_focks.data[:, fock_index],
                                                 self.comm,
                                                 distribute=False)
                fock_index += 1

        write_distributed_focks(fock_file, focks, keys, wi, self.comm,
                                self.ostream)

        return focks

    def get_fock_dict_ii(self, wi, density_list, F0, mo, molecule, ao_basis):
        """
        Computes the Fock matrices for a quadratic response function

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

        keys = [
            'Fb_cd',
            'Fcd_b',
            'Fc_bd',
            'Fbd_c',
            'Fd_bc',
            'Fbc_d',
        ]

        if self.checkpoint_file is not None:
            fock_file = str(
                Path(self.checkpoint_file).with_suffix('.crf_fock_2_full.h5'))
        else:
            fock_file = None

        if self.restart:
            if self.rank == mpi_master():
                self.restart = check_distributed_focks(fock_file, keys, wi)
            self.restart = self.comm.bcast(self.restart, mpi_master())

        if self.restart:
            focks = read_distributed_focks(fock_file, keys, wi, self.comm,
                                           self.ostream)
            focks['F0'] = F0
            return focks

        time_start_fock = time.time()
        dist_focks = self.get_fock(mo, density_list, molecule, ao_basis,
                                     'real_and_imag')
        time_end_fock = time.time()

        total_time_fock = time_end_fock - time_start_fock
        self.print_fock_time(total_time_fock)
        
        focks = {'F0': F0}
        for key in keys:
            focks[key] = {}

        fock_index = 0
        for (wb,wc,wd) in wi:
            for key in keys:
                focks[key][(wb,wc,wd)] = DistributedArray(dist_focks.data[:, fock_index],
                                                 self.comm,
                                                 distribute=False)
                fock_index += 1

        write_distributed_focks(fock_file, focks, keys, wi, self.comm,
                                self.ostream)

        return focks



    def get_esr4(self, wi, kX, fo,fo2, nocc, norb,D0):
        """
        Contracts E[4], S[4], R[4]

        :param wi:
            A list of freqs
        :param kX:
            A dict of the single index response matricies
        :param kXY:
            A dict of the two index response matrices
        :param fo:
            A dictonary of transformed Fock matricies from fock_dict
        :param fo2:
            A dictonarty of transfromed Fock matricies from fock_dict_two
        :param nocc:
            The number of occupied orbitals
        :param norb:
            The total number of orbitals

        :return:
            A dictonary of E[4], S[4], R[4] tensor contractions 
        """

        e4_vec = {}
        s4_vec = {}
        r4_vec = {}

        for (wb,wc,wd) in wi:

            vec_pack = np.array([
                fo2[('B', wb)].data,
                fo2[('C', wc)].data,
                fo2[('D', wd)].data,
                fo['Fbc'][wb].data,
                fo['Fcb'][wb].data,
                fo['Fbd'][wb].data,
                fo['Fdb'][wb].data,
                fo['Fcd'][wb].data,
                fo['Fdc'][wb].data,
                fo['Fbcd'][wb].data,
                fo['Fbdc'][wb].data,
                fo['Fcbd'][wb].data,
                fo['Fcdb'][wb].data,
                fo['Fdbc'][wb].data,
                fo['Fdcb'][wb].data,
            ]).T.copy()

            vec_pack = self.collect_vectors_in_columns(vec_pack)

            if self.rank != mpi_master():
                continue

            vec_pack = vec_pack.T.copy().reshape(-1, norb, norb)

            (fb, fc, fd, fbc, fcb,fbd,fdb,fcd,fdc,fbcd,fbdc,fcbd,fcdb,fdbc,fdcb) = vec_pack

            fb = np.conjugate(fb).T
            fc = np.conjugate(fc).T
            fd = np.conjugate(fd).T

            F0_a = fo['F0']

            # Response

            kb = kX[('B', wb)].T

            kc = kX[('C', wc)].T

            kd = kX[('D', wd)].T


            zi_bcd = self.zi(kb, kc, kd, fc, fd, fcd, fdc, F0_a)

            zi_cbd = self.zi(kc, kb, kd, fb, fd, fbd, fdb, F0_a)

            zi_dbc = self.zi(kd, kb, kc, fb, fc, fbc, fcb, F0_a)

            e4fock = (zi_bcd + zi_cbd + zi_dbc) + ( fbcd + fbdc + fcbd + fcdb + fdbc + fdcb )

            e4vec = 2./ 6 * self.anti_sym(LinearSolver.lrmat2vec(e4fock.T, nocc, norb))

            ka = kX[('A', (wb+wc+wd))]
            kb = kX[('B', wb)]
            kc = kX[('C', wc)]
            kd = kX[('D', wd)]

            s4_term =  wb * self.s4(kb, kc, kd, D0, nocc, norb)

            s4_term += wc * self.s4(kc, kb, kd, D0, nocc, norb)

            s4_term += wd * self.s4(kd, kb, kc, D0, nocc, norb)

            Nb = (LinearSolver.lrmat2vec(kb.real, nocc, norb) +
                  1j * LinearSolver.lrmat2vec(kb.imag, nocc, norb))

            Nc = (LinearSolver.lrmat2vec(kc.real, nocc, norb) +
                     1j * LinearSolver.lrmat2vec(kc.imag, nocc, norb))

            Nd = (LinearSolver.lrmat2vec(kd.real, nocc, norb) +
                  1j * LinearSolver.lrmat2vec(kd.imag, nocc, norb))

            Nb_h = self.flip_xy(Nb)

            Nc_h = self.flip_xy(Nc)

            Nd_h = self.flip_xy(Nd)

            r4_term = - 1j * self.damping * np.dot(Nd_h, self.s4_for_r4(ka.T, kb, kc, D0, nocc, norb))

            r4_term += -1j * self.damping * np.dot(Nc_h, self.s4_for_r4(ka.T, kb, kd, D0, nocc, norb))
                
            r4_term += -1j * self.damping * np.dot(Nd_h, self.s4_for_r4(ka.T, kc, kb, D0, nocc, norb))

            r4_term += -1j * self.damping * np.dot(Nb_h, self.s4_for_r4(ka.T, kc, kd, D0, nocc, norb))

            r4_term += -1j * self.damping * np.dot(Nc_h, self.s4_for_r4(ka.T, kd, kb, D0, nocc, norb))

            r4_term += -1j * self.damping * np.dot(Nb_h, self.s4_for_r4(ka.T, kd, kc, D0, nocc, norb))

            e4_vec.update({wb: e4vec})

            s4_vec.update({wb: s4_term})

            r4_vec.update({wb: r4_term})

        return e4_vec, s4_vec, r4_vec

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

    def get_nxy(self, wi, kX, fo,fo2, nocc, norb,d_a_mo,X,molecule,ao_basis,scf_tensors):
        """
        Computed NXY

        :param wi:
            A list of freqs
        :param kX:
            A dict of the single index response matricies
        :param kXY:
            A dict of the two index response matrices
        :param fo:
            A dictonary of transformed Fock matricies from fock_dict
        :param fo2:
            A dictonarty of transfromed Fock matricies from fock_dict_two
        :param nocc:
            The number of occupied orbitals
        :param norb:
            The total number of orbitals

        :return:
            A dictonary of E[4], S[4], R[4] tensor contractions 
        """

        BC = {}
        CD = {}
        BD = {}

        XY = {}

        for (wb,wc,wd) in wi:

            vec_pack = np.array([
                fo2[('B', wb)].data,
                fo2[('C', wc)].data,
                fo2[('D', wd)].data,
                fo['Fbc'][wb].data,
                fo['Fcb'][wb].data,
                fo['Fbd'][wb].data,
                fo['Fdb'][wb].data,
                fo['Fcd'][wb].data,
                fo['Fdc'][wb].data,
                fo['Fbcd'][wb].data,
                fo['Fbdc'][wb].data,
                fo['Fcbd'][wb].data,
                fo['Fcdb'][wb].data,
                fo['Fdbc'][wb].data,
                fo['Fdcb'][wb].data,
            ]).T.copy()

            vec_pack = self.collect_vectors_in_columns(vec_pack)

            if self.rank != mpi_master():
                continue

            vec_pack = vec_pack.T.copy().reshape(-1, norb, norb)

            (fb, fc, fd, fbc, fcb,fbd,fdb,fcd,fdc,fbcd,fbdc,fcbd,fcdb,fdbc,fdcb) = vec_pack

            fb = np.conjugate(fb).T
            fc = np.conjugate(fc).T
            fd = np.conjugate(fd).T

            F0_a = fo['F0']

            # Response

            kb = kX[('B', wb)].T

            kc = kX[('C', wc)].T

            kd = kX[('D', wd)].T

            B = X[self.b_component]
            C = X[self.c_component]
            D = X[self.d_component]

            Nb = (LinearSolver.lrmat2vec(kX[('B',wb)].real, nocc, norb) +
                1j * LinearSolver.lrmat2vec(kX[('B',wb)].imag, nocc, norb))

            Nc = (LinearSolver.lrmat2vec(kX[('C',wc)].real, nocc, norb) +
                1j * LinearSolver.lrmat2vec(kX[('C',wc)].imag, nocc, norb))
            
            Nd = (LinearSolver.lrmat2vec(kX[('D',wd)].real, nocc, norb) +
                1j * LinearSolver.lrmat2vec(kX[('D',wd)].imag, nocc, norb))
            

        # BC

            xi = self.xi(kb, kc, fb, fc, F0_a)

            e3fock = xi.T + (0.5 * fbc + 0.5*fcb).T
            E3NbNc = self.anti_sym(-LinearSolver.lrmat2vec(e3fock, nocc, norb))

            C2Nb = 0.5*self.x2_contract(kX[('B',wb)], C, d_a_mo, nocc, norb)
            B2Nc = 0.5*self.x2_contract(kX[('C',wc)], B, d_a_mo, nocc, norb)

            BC.update({(('BC',wb,wc),wb+wc): (E3NbNc - C2Nb  - B2Nc)})

        
        # BD

            xi = self.xi(kb, kd, fb, fd, F0_a)

            e3fock = xi.T + (0.5*fbd + 0.5*fdb).T
            E3NbNd = self.anti_sym(- LinearSolver.lrmat2vec(e3fock, nocc, norb))

            D2Nb = 0.5*self.x2_contract(kX[('B',wb)], D, d_a_mo, nocc, norb)
            B2Nd = 0.5*self.x2_contract(kX[('D',wd)], B, d_a_mo, nocc, norb)

            BD.update({ (('BD',wb,wd),wb+wd): (E3NbNd - D2Nb - B2Nd)})



        # CD

            xi = self.xi(kc, kd, fc, fd, F0_a)

            e3fock = xi.T + (0.5 * fcd + 0.5*fdc).T
            E3NcNd = self.anti_sym(-  LinearSolver.lrmat2vec(e3fock, nocc, norb))

            C2Nd = 0.5*self.x2_contract(kX[('D',wd)], C, d_a_mo, nocc, norb)
            D2Nc = 0.5*self.x2_contract(kX[('C',wc)], D, d_a_mo, nocc, norb)

            CD.update({ (('CD',wc,wd),wc+wd): (E3NcNd - C2Nd  - D2Nc)})

            XY.update(BC)
            XY.update(BD)
            XY.update(CD)

        Nxy_drv = ComplexResponse(self.comm, self.ostream)

        Nxy_drv.update_settings({
            'damping': self.damping,
            'lindep_thresh': self.lindep_thresh,
            'conv_thresh': self.conv_thresh,
            'max_iter': self.max_iter,
            'eri_thresh': self.eri_thresh,
            'qq_type': self.qq_type,
        })

        Nxy_drv.timing = self.timing
        Nxy_drv.memory_profiling = self.memory_profiling
        Nxy_drv.batch_size = self.batch_size
        Nxy_drv.restart = self.restart
        Nxy_drv.program_start_time = self.program_start_time
        Nxy_drv.maximum_hours = self.maximum_hours
        if self.checkpoint_file is not None:
            Nxy_drv.checkpoint_file = re.sub(r'\.h5$', r'', self.checkpoint_file)
            Nxy_drv.checkpoint_file += '_crf_2.h5'

        Nxy_results = Nxy_drv.compute(molecule, ao_basis, scf_tensors, XY)

        kX = {}
        Focks = {}

        kX = Nxy_results['kappas']

        Focks = Nxy_results['focks']


        return kX,Focks



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

    def print_results(self, freqs, gamma, comp, t4_dict, t3_dict, tpa_dict):
        """
        Prints the results from the TPA calculation.

        :param freqs:
            List of frequencies
        :param gamma:
            A dictonary containing the isotropic cubic response functions for
            TPA
        :param comp:
            List of gamma tensors components
        :param t4_dict:
            A dictonary containing the isotropic T[4] contractions
        :param t3_dict:
            A dictonary containing the isotropic T[3] contractions
        :param tpa_dict:
            A dictonary containing the isotropic X[3], A[3], X[2], A[2]
            contractions
        """

        return None

    def print_header(self):
        """
        Prints TPA setup header to output stream.
        """

        self.ostream.print_blank()

        title = 'Two-Photon Absorbtion Driver Setup'
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

        return (np.linalg.multi_dot([k, S, D]) -
                np.linalg.multi_dot([D, S, k]))

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
        X3NxNy_c = (LinearSolver.lrmat2vec(X3NxNy.real, nocc, norb) +
                    1j * LinearSolver.lrmat2vec(X3NxNy.imag, nocc, norb))
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
        X2Nx_c = (LinearSolver.lrmat2vec(XNx.real, nocc, norb) +
                  1j * LinearSolver.lrmat2vec(XNx.imag, nocc, norb))
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
        A2Nx_c = (LinearSolver.lrmat2vec(ANx.real, nocc, norb) +
                  1j * LinearSolver.lrmat2vec(ANx.imag, nocc, norb))
        return -(1. / 2) * A2Nx_c

    
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
        
        return (self.commut(kB, M1 + M2 + 3*(Fbc+ Fcb) ) )


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


    def get_fock(self, mo, D, molecule, ao_basis, fock_flag):
        """
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
            if self.rank == mpi_master():
                D_total = []
                for da in D:
                    D_total.append(da.real)
                    D_total.append(da.imag)
            else:
                D_total = None

            f_total = self.get_two_el_fock(mo, molecule, ao_basis,
                                                 D_total)

            nrows = f_total.data.shape[0]
            half_ncols = f_total.data.shape[1] // 2
            ff_data = np.zeros((nrows, half_ncols), dtype=np.complex128)
            for i in range(half_ncols):
                ff_data[:, i] = (f_total.data[:, 2 * i] +
                                 1j * f_total.data[:, 2 * i + 1])
            return DistributedArray(ff_data, self.comm, distribute=False)

        elif fock_flag == 'real':
            return self.get_two_el_fock(mo, molecule, ao_basis, D)

        else:
            return None

    def get_two_el_fock(self, mo, molecule, ao_basis, dabs):
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
        Prints TPA component.

        :param label:
            The label
        :param freq:
            The frequency
        :param value:
            The complex value
        :param width:
            The width for the output
        """

        w_str = '{:<9s} {:20.8f} {:20.8f}j'.format(
            label, value.real, value.imag)
        self.ostream.print_header(w_str.ljust(width))

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

    def s4(self, k1, k2, k3, D, nocc, norb):
        """
        Returns the contraction of S[4] for S[4] dict

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
        Returns the contraction of S[4] for S[4] dict

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
        S4N1N2N3_c = (LinearSolver.lrmat2vec(S4N1N2N3.real, nocc, norb) +
                      1j * LinearSolver.lrmat2vec(S4N1N2N3.imag, nocc, norb))
        return (2. / 6) * S4N1N2N3_c

    def flip_xy(self, X):
        """
        Swaps upper and lower parts of a response vector. This is used when
        rewriting the R^[4] tensor contraction in terms of S^[4] tensor
        contractions.

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