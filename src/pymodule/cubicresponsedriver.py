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
        self.dcomponent = 'z'

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

        if 'grid_level' in method_dict:
            self.grid_level = int(method_dict['grid_level'])

        if 'xcfun' in method_dict:
            if 'dft' not in method_dict:
                self.dft = True
            self.xcfun = parse_xc_func(method_dict['xcfun'].upper())
            assert_msg_critical(not self.xcfun.is_undefined(),
                                'Response solver: Undefined XC functional')

        if 'b_frequencies' in rsp_dict:
            self.b_frequencies = parse_seq_range(rsp_dict['b_frequencies'])

        if 'c_frequencies' in rsp_dict:
            self.c_frequencies = parse_seq_range(rsp_dict['c_frequencies'])
        
        if 'd_frequencies' in rsp_dict:
            self.d_frequencies = parse_seq_range(rsp_dict['d_frequencies'])

        if 'a_component' in rsp_dict:
            self.a_component = rsp_dict['a_component']
        
        if 'b_component' in rsp_dict:
            self.b_component = rsp_dict['b_component']

        if 'c_component' in rsp_dict:
            self.c_component = rsp_dict['c_component']
        
        if 'd_component' in rsp_dict:
            self.dcomponent = rsp_dict['d_component']

        if 'damping' in rsp_dict:
            self.damping = float(rsp_dict['damping'])

        if 'eri_thresh' in rsp_dict:
            self.eri_thresh = float(rsp_dict['eri_thresh'])
        if 'qq_type' in rsp_dict:
            self.qq_type = rsp_dict['qq_type']
        if 'batch_size' in rsp_dict:
            self.batch_size = int(rsp_dict['batch_size'])

        if 'max_iter' in rsp_dict:
            self.max_iter = int(rsp_dict['max_iter'])
        if 'conv_thresh' in rsp_dict:
            self.conv_thresh = float(rsp_dict['conv_thresh'])
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

        d_rhs = linear_solver.get_complex_prop_grad(operator, self.dcomponent,
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


        # Computing the first-order response vectors (3 per frequency)
        Na_drv = ComplexResponse(self.comm, self.ostream)

        Na_drv.update_settings({
            'frequencies': wa,
            'damping': self.damping,
            'lindep_thresh': self.lindep_thresh,
            'conv_thresh': self.conv_thresh,
            'max_iter': self.max_iter,
            'eri_thresh': self.eri_thresh,
            'qq_type': self.qq_type,
        })

        Na_drv.timing = self.timing
        Na_drv.memory_profiling = self.memory_profiling
        Na_drv.batch_size = self.batch_size
        Na_drv.restart = self.restart
        Na_drv.program_start_time = self.program_start_time
        Na_drv.maximum_hours = self.maximum_hours
        if self.checkpoint_file is not None:
            Na_drv.checkpoint_file = re.sub(r'\.h5$', r'', self.checkpoint_file)
            Na_drv.checkpoint_file += '_quada_1.h5'

        Na_results = Na_drv.compute(molecule, ao_basis, scf_tensors, A)

        Nb_drv = ComplexResponse(self.comm, self.ostream)

        Nb_drv.update_settings({
            'frequencies': self.b_frequencies,
            'damping': self.damping,
            'lindep_thresh': self.lindep_thresh,
            'conv_thresh': self.conv_thresh,
            'max_iter': self.max_iter,
            'eri_thresh': self.eri_thresh,
            'qq_type': self.qq_type,
        })

        Nb_drv.timing = self.timing
        Nb_drv.memory_profiling = self.memory_profiling
        Nb_drv.batch_size = self.batch_size
        Nb_drv.restart = self.restart
        Nb_drv.program_start_time = self.program_start_time
        Nb_drv.maximum_hours = self.maximum_hours
        if self.checkpoint_file is not None:
            Nb_drv.checkpoint_file = re.sub(r'\.h5$', r'', self.checkpoint_file)
            Nb_drv.checkpoint_file += '_quadb_1.h5'

        Nb_results = Nb_drv.compute(molecule, ao_basis, scf_tensors, B)

        Nc_drv = ComplexResponse(self.comm, self.ostream)

        Nc_drv.update_settings({
            'frequencies': self.c_frequencies,
            'damping': self.damping,
            'lindep_thresh': self.lindep_thresh,
            'conv_thresh': self.conv_thresh,
            'max_iter': self.max_iter,
            'eri_thresh': self.eri_thresh,
            'qq_type': self.qq_type,
        })

        Nc_drv.timing = self.timing
        Nc_drv.memory_profiling = self.memory_profiling
        Nc_drv.batch_size = self.batch_size
        Nc_drv.restart = self.restart
        Nc_drv.program_start_time = self.program_start_time
        Nc_drv.maximum_hours = self.maximum_hours
        if self.checkpoint_file is not None:
            Nc_drv.checkpoint_file = re.sub(r'\.h5$', r'', self.checkpoint_file)
            Nc_drv.checkpoint_file += '_quadc_1.h5'

        Nc_results = Nc_drv.compute(molecule, ao_basis, scf_tensors, C)

        Nd_drv = ComplexResponse(self.comm, self.ostream)

        Nd_drv.update_settings({
            'frequencies': self.c_frequencies,
            'damping': self.damping,
            'lindep_thresh': self.lindep_thresh,
            'conv_thresh': self.conv_thresh,
            'max_iter': self.max_iter,
            'eri_thresh': self.eri_thresh,
            'qq_type': self.qq_type,
        })

        Nd_drv.timing = self.timing
        Nd_drv.memory_profiling = self.memory_profiling
        Nd_drv.batch_size = self.batch_size
        Nd_drv.restart = self.restart
        Nd_drv.program_start_time = self.program_start_time
        Nd_drv.maximum_hours = self.maximum_hours
        if self.checkpoint_file is not None:
            Nd_drv.checkpoint_file = re.sub(r'\.h5$', r'', self.checkpoint_file)
            Nd_drv.checkpoint_file += '_quadc_1.h5'

        Nd_results = Nd_drv.compute(molecule, ao_basis, scf_tensors, D)

        kX = {}
        Focks = {}

        kX = Nb_results['kappas']

        kX.update(Na_results['kappas'])

        kX.update(Nc_results['kappas'])

        kX.update(Nd_results['kappas'])


        Focks = Nb_results['focks']

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
        :param track:
            A list that contains all the information about which γ components
            and at what freqs they are to be computed
        :param scf_tensors:
            The dictionary of tensors from converged SCF wavefunction.
        :param molecule:
            The molecule.
        :param basis:
            The AO basis.
        :param profiler:
            The profiler.

        :return:
            A dictionary containing all the relevent terms to third-order
            isotropic gradient
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

        e4_dict,s4_dict,r4_dict = self.get_e4(freqpairs, kX, fock_dict, nocc,norb,d_a_mo)

        op_a = X[self.a_component]
        op_b = X[self.b_component]
        op_c = X[self.c_component]
        op_d = X[self.dcomponent]
        
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
            

            NaE4NbNcNd = np.dot(Na, e4_dict[wb])

            NaS4NbNcNd = np.dot(Na, s4_dict[wb])

            NaR4NbNcNd =  r4_dict[wb]

            self.ostream.print_blank()
            w_str = 'Cubic response function: ' + '<< ' + str(self.a_component) +';' + str(self.b_component) + ',' + str(self.c_component) +  ',' + str(self.dcomponent) + ' >> '
            self.ostream.print_header(w_str)
            self.ostream.print_header('=' * (len(w_str) + 2))
            self.ostream.print_blank()
            title = '{:<9s} {:>12s} {:>20s} {:>21s}'.format(
                'Component', 'Frequency', 'Real', 'Imaginary')
            width = len(title)
            self.ostream.print_header(title.ljust(width))
            self.ostream.print_header(('-' * len(title)).ljust(width))
            self.print_component('T4',wb, NaE4NbNcNd-NaS4NbNcNd+NaR4NbNcNd, width)
            result.update({wb: NaE4NbNcNd+NaS4NbNcNd-NaR4NbNcNd})
        

        profiler.check_memory_usage('End of QRF')

        return result

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

                density_list.append(Db)
                density_list.append(Dc)
                density_list.append(Dd)

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
            'Fb',
            'Fc',
            'Fd',
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

    def get_e4(self, wi, kX, fo, nocc, norb,D0):
        """
        Contracts E[4]

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
            response function for TPA
        """

        e4_vec = {}
        s4_vec = {}
        r4_vec = {}

        for (wb,wc,wd) in wi:

            vec_pack = np.array([
                fo['Fb'][wb].data,
                fo['Fc'][wb].data,
                fo['Fd'][wb].data,
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

    def xi(self, kA, kB, Fa, Fb, F0):
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

        return 0.5 * (self.commut(kA,
                                  self.commut(kB, F0) + 2 * Fb) +
                      self.commut(kB,
                                  self.commut(kA, F0) + 2 * Fa))

    
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

    def print_component(self, label, freq, value, width):
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

        w_str = '{:<9s} {:12.4f} {:20.8f} {:20.8f}j'.format(
            label, freq, value.real, value.imag)
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