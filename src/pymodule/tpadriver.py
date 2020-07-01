import numpy as np
import ctypes
import psutil
import time
import os
import re

from .veloxchemlib import ElectronRepulsionIntegralsDriver
from .veloxchemlib import KineticEnergyIntegralsDriver
from .veloxchemlib import NuclearPotentialIntegralsDriver
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
from .inputparser import parse_frequencies
from .errorhandler import assert_msg_critical


class TpaDriver:
    """
    Implements the isotropic cubic response driver for two-photon absorption
    (TPA)

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
        - timing: The flag for printing timing information.
        - profiling: The flag for printing profiling information.
        - memory_profiling: The flag for printing memory usage.
        - memory_tracing: The flag for tracing memory allocation.
    """

    def __init__(self, comm, ostream):
        """
        Initializes the isotropic cubic response driver for two-photon
        absorption (TPA)
        """

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

    def update_settings(self, rsp_dict, method_dict=None):
        """
        Updates response and method settings in TPA driver

        :param rsp_dict:
            The dictionary of response dict.
        :param method_dict:
            The dictionary of method rsp_dict.
        """

        if method_dict is None:
            method_dict = {}

        if 'frequencies' in rsp_dict:
            self.frequencies = parse_frequencies(rsp_dict['frequencies'])
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
        Computes the isotropic cubic response function for two-photon
        absorption

        :param molecule:
            The molecule.
        :param basis:
            The AO basis.
        :param scf_tensors:
            The dictionary of tensors from converged SCF wavefunction.

        :return:
              A dictonary containing the isotropic T[4], T[3], X[3], A[3],
              X[2], A[2] contractions and the isotropic cubic response
              functions for TPA
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
                            'TpaDriver: not implemented for unrestricted case')

        if self.rank == mpi_master():
            S = scf_tensors['S']
            da = scf_tensors['D'][0]
            mo = scf_tensors['C']
            d_a_mo = np.linalg.multi_dot([mo.T, S, da, S, mo])
            nocc = molecule.number_of_alpha_electrons()
            norb = mo.shape[1]
        else:
            d_a_mo = None

        # Computing first-order gradient vectors
        dipole_drv = ElectricDipoleIntegralsDriver(self.comm)
        dipole_mats = dipole_drv.compute(molecule, ao_basis)

        operator = 'dipole'
        component = 'xyz'

        linear_solver = LinearSolver(self.comm, self.ostream)
        b_rhs = linear_solver.get_complex_rhs(operator, component, molecule,
                                              ao_basis, scf_tensors)

        # Storing the dipole integral matrices used for the X[3],X[2],A[3] and
        # A[2] contractions in MO basis
        if self.rank == mpi_master():
            v1 = {(op, w): v for op, v in zip(component, b_rhs)
                  for w in self.frequencies}
            X = {
                'x': 2 * self.ao2mo(mo, dipole_mats.x_to_numpy()),
                'y': 2 * self.ao2mo(mo, dipole_mats.y_to_numpy()),
                'z': 2 * self.ao2mo(mo, dipole_mats.z_to_numpy())
            }
            self.comp = self.get_comp(self.frequencies)
        else:
            v1 = None
            X = None
            self.comp = None

        # Computing the first-order response vectors (3 per frequency)
        Nb_drv = ComplexResponse(self.comm, self.ostream)

        Nb_drv.update_settings({
            'frequencies': self.frequencies,
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
        if self.checkpoint_file is not None:
            Nb_drv.checkpoint_file = re.sub(r'\.h5$', r'', self.checkpoint_file)
            Nb_drv.checkpoint_file += '_tpa_1.h5'

        Nb_results = Nb_drv.compute(molecule, ao_basis, scf_tensors, v1)

        Nx = {}
        kX = {}
        Focks = {}

        if self.rank == mpi_master():
            Nx['Nb'] = Nb_results['solutions']
            kX['Nb'] = Nb_results['kappas']
            Focks['Fb'] = Nb_results['focks']

            Nx['Nc'] = {}
            kX['Nc'] = {}
            Focks['Fc'] = {}
            Focks['Fd'] = {}

            # The first-order response vectors with negative frequency are
            # obtained from the first-order response vectors with positive
            # frequency by using flip_zy, see article.

            for (op, w) in Nx['Nb']:
                Nx['Nc'][(op, -w)] = self.flip_yz(Nx['Nb'][(op, w)])

                # Creating the response matrix for the negative first-order
                # response vectors

                kX['Nc'][(op, -w)] = (
                    LinearSolver.lrvec2mat(Nx['Nc'][(op, -w)].real, nocc,
                                           norb) + 1j *
                    LinearSolver.lrvec2mat(Nx['Nc'][(op, -w)].imag, nocc, norb))

                # The first-order Fock matrices with positive and negative
                # wuencies are each other complex conjugates

                Focks['Fc'][(op, -w)] = Focks['Fb'][(op, w)]
                Focks['Fd'][(op, w)] = np.conjugate(Focks['Fb'][(op, w)]).T

            # For cubic-response with all operators being the dipole μ Nb=Na=Nd
            # Likewise, Fb=Fd

            Focks['Fb'].update(Focks['Fd'])

            Nx['Na'] = Nx['Nb']
            kX['Na'] = kX['Nb']

            Nx['Nd'] = Nx['Nb']
            kX['Nd'] = kX['Nb']

        profiler.check_memory_usage('1st CPP')

        # Computing the third-order gradient and also the contractions of
        # A[3] and A[2] which formally are not part of the third-order gradient
        # but which are used for the cubic response function

        tpa_dict = self.main(Focks, Nx, self.frequencies, X, d_a_mo, kX,
                             self.comp, scf_tensors, molecule, ao_basis,
                             profiler)

        valstr = '*** Time spent in TPA calculation: {:.2f} sec ***'.format(
            time.time() - start_time)
        self.ostream.print_header(valstr)
        self.ostream.print_blank()
        self.ostream.flush()

        profiler.end(self.ostream)

        self.is_converged = True

        return tpa_dict

    def main(self, Focks, Nx, w, X, d_a_mo, kX, track, scf_tensors, molecule,
             ao_basis, profiler):
        """
        Computes all the relevent terms to third-order isotropic gradient

        :param Nx:
            A dictonary containing all the single index response vectors
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
            nocc = molecule.number_of_alpha_electrons()
            norb = mo.shape[1]
        else:
            D0 = None
            mo = None
            nocc = None
            norb = None

        # computing all compounded first-order densities
        if self.rank == mpi_master():
            density_list = self.get_densities(w, kX, S, D0, mo)
        else:
            density_list = None

        profiler.check_memory_usage('1st densities')

        #  computing the compounded first-order Fock matrices
        fock_dict = self.get_fock_dict(w, density_list, (D0, D0), mo, molecule,
                                       ao_basis)

        if self.rank == mpi_master():
            fock_dict.update(Focks)
            e4_dict = self.get_e4(w, kX, fock_dict, nocc, norb)

        profiler.check_memory_usage('1st Focks')

        # computing all the compounded second-order response vectors and
        # extracting some of the second-order Fock matrices from the subspace
        (n_xy_dict, kxy_dict, Focks_xy,
         XΥ_dict) = self.get_n_xy(w, d_a_mo, X, fock_dict, kX, nocc, norb,
                                  molecule, ao_basis, scf_tensors)

        profiler.check_memory_usage('2nd CPP')

        # computing all second-order compounded densities based on the
        # second-order response vectors
        if self.rank == mpi_master():
            density_list_two = self.get_densities_II(w, kX, kxy_dict, S, D0, mo)
        else:
            density_list_two = None

        profiler.check_memory_usage('2nd densities')

        # computing the remaning second-order Fock matrices from the
        # second-order densities
        fock_dict_two = self.get_fock_dict_II(w, density_list_two, mo, molecule,
                                              ao_basis)

        profiler.check_memory_usage('2nd Focks')

        if self.rank == mpi_master():
            # Adding the Fock matrices extracted from the second-order response
            # vector subspace to the fock_dict's.
            fock_dict_two.update(Focks_xy)

            # computing the compounded E[3] contractions for the isotropic
            # cubic response function
            e3_dict = self.get_e3(w, kX, kxy_dict, fock_dict, fock_dict_two,
                                  nocc, norb)

            # computing the X[3],A[3],X[2],A[2] contractions for the isotropic
            # cubic response function
            other_dict = self.other(w, track, Nx, n_xy_dict, X, kX, kxy_dict,
                                    d_a_mo, nocc, norb)

        profiler.check_memory_usage('Other terms')

        # Combining all the terms to evaluate the iso-tropic cubic response
        # function. For TPA Full and reduced, see article

        result = {}

        if self.rank == mpi_master():
            t4_dict = self.get_t4(self.frequencies, e4_dict, Nx, kX, self.comp,
                                  d_a_mo, nocc, norb)
            t3_dict = self.get_t3(self.frequencies, e3_dict, Nx, self.comp)

            gamma = {}

            for w in self.frequencies:
                sum_val = t3_dict[(w, -w, w)]
                if t4_dict is not None:
                    sum_val += t4_dict[(w, -w, w)]
                for key, val in other_dict.items():
                    sum_val += val[(w, -w, w)]
                gamma[(w, -w, w)] = 1. / 15 * sum_val

            self.print_results(self.frequencies, gamma, self.comp, t4_dict,
                               t3_dict, other_dict)

            result.update(other_dict)

            result.update({
                't4_dict': t4_dict,
                't3_dict': t3_dict,
                'gamma': gamma,
                'w': self.frequencies,
            })

        profiler.check_memory_usage('End of TPA')

        return result

    def get_t4(self, wi, e4_dict, n_x, kX, track, da, nocc, norb):
        """
        Computes the contraction of the E[4] tensor with that of the S[4] and
        R[4] tensors to return the contraction of T[4] as a dictonary of
        vectors. T[4]n_xNyNz = (E^[4]-ω_1S^[4]-ω_1S^[4]-ω_3S^[4]-γiR^[4])

        :param wi:
            A list of all the freqs
        :param e4_dict:
            A dictonary of all the E[4] contraction
        :param n_x:
            A dictonary with all the single index response vectors
        :param kX:
            A dictonray containng all the response matricies
        :param track:
            A list containg information about all the γ components that are to
            be computed
        :param da:
            The SCF density matrix in MO basis
        :param nocc:
            The number of occupied orbitals
        :param norb:
            The total number of orbitals

        :return:
            A dictonary of final T[4] contraction values
        """

        return None

    def get_t3(self, freqs, e3_dict, n_x, track):
        """
        Computes the T[3] contraction, for HF S[3] = 0, R[3] = 0 such that
        the T[3] contraction for the isotropic cubic response function in terms
        of compounded Fock matrices is given as:

                             [(ζ_{α}^{σσ} + ζ_{α}^{λλ+ττ} + f_{α}^{λσ,τ})_is]
        t3term = Σ_{α} N_{α} [(ζ_{α}^{σσ} + ζ_{α}^{λλ+ττ} + f_{α}^{λσ,τ})_si]

        For more details see article

        :param freqs:
            List of frequencies of the pertubations
        :param e3_dict:
            A dictonary that contains the contractions of E[3]
        :param n_x:
            A dictonary containing the response vectors n_x = (E[2]-wS[2])^-1
            X[1]
        :param track:
            A list containing information about what tensor components that are
            being computed

        :return:
            A dictonary of the final values for the NaT[3]NxNyz contractions
        """

        t3_term = {}

        for i in range(len(freqs)):
            w = float(track[i * (len(track) // len(freqs))].split(",")[1])

            t3term = -(np.matmul(n_x['Na'][('x', w)], e3_dict['f_iso_x'][w]) +
                       np.matmul(n_x['Na'][('y', w)], e3_dict['f_iso_y'][w]) +
                       np.matmul(n_x['Na'][('z', w)], e3_dict['f_iso_z'][w]))

            t3_term[(w, -w, w)] = t3term

        return t3_term

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

    def get_comp(self, freqs):
        """
        Makes a list of all the gamma tensor components that are to be computed
        for printing purposes and for the contraction of X[3],X[2],A[3],A[2]

        :param freqs:
            A list of all the frequencies for the TPA calculation

        :return:
            A list of gamma tensors components inlcuded in the isotropic cubic
            response with their corresponding frequencies
        """

        comp_iso = []
        spat_A = 'xyz'

        for w in freqs:
            w_key = '{},{},{}'.format(w, -w, w)
            for b in spat_A:
                for a in spat_A:
                    aabb = '{}{}{}{}'.format(a, a, b, b)
                    abab = '{}{}{}{}'.format(a, b, a, b)
                    abba = '{}{}{}{}'.format(a, b, b, a)
                    comp_iso.append(aabb + ',' + w_key)
                    comp_iso.append(abab + ',' + w_key)
                    comp_iso.append(abba + ',' + w_key)

        return sorted(comp_iso, key=comp_iso.index)

    def get_e4(self, wi, kX, fo, nocc, norb):
        """
        Contracts E[4]n_xNyNz for the isotropic cubic response function. Takes
        the Fock matrices from fock_dict and contracts them with the response
        vectors.

        :param wi:
            A list of freqs
        :param kX:
            A dict of the single index response matricies
        :param fo:
            A dictonary of transformed Fock matricies from fock_dict
        :param nocc:
            The number of occupied orbitals
        :param norb:
            The total number of orbitals

        :return:
            A dictonary of compounded E[4] tensors for the isotropic cubic
            response function for TPA
        """

        return None

    def flip_xy(self, X):

        # TODO: add docstring for flip_xy

        if X.ndim == 1:
            new_xy = np.zeros_like(X)
            half_len = X.shape[0] // 2
            new_xy[:half_len] = X[half_len:]
            new_xy[half_len:] = X[:half_len]
            return new_xy

        return None

    def flip_yz(self, X):

        # TODO: add docstring for flip_yz

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

        return (np.linalg.multi_dot([k.T, S, D]) -
                np.linalg.multi_dot([D, S, k.T]))

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

        Xn_x = self.commut(self.commut(k, X), D.T)
        X2n_x_c = (LinearSolver.lrmat2vec(Xn_x.real, nocc, norb) +
                   1j * LinearSolver.lrmat2vec(Xn_x.imag, nocc, norb))
        return X2n_x_c

    def x3_contract(self, k1, k2, X, D, nocc, norb):
        """
        Contracts the generalized dipole gradient tensor of rank 3 with two
        first-order response matrices. X[3]N1N2 = -(1/2)[[k2,[k1,X]],D.T]

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

        X3n_xNy = self.commut(self.commut(k2, self.commut(k1, X)), D.T)
        X3n_xNy = [
            LinearSolver.lrmat2vec(X3n_xNy.real, nocc, norb),
            LinearSolver.lrmat2vec(X3n_xNy.imag, nocc, norb)
        ]
        X3n_xNy_c = X3n_xNy[0] + 1j * X3n_xNy[1]
        return -(1. / 2) * X3n_xNy_c

    def a3_contract(self, k1, k2, A, D, nocc, norb):
        """
        Contracts the generalized dipole gradient tensor of rank 3 with two
        first-order response matrices. A[3]N1N2 = (1/6)[[k2,[k1,A]],D.T]

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

        A3n_xNy = self.commut(self.commut(k2.T, self.commut(k1.T, A)), D.T)
        A3n_xNy = [
            LinearSolver.lrmat2vec(A3n_xNy.real, nocc, norb),
            LinearSolver.lrmat2vec(A3n_xNy.imag, nocc, norb)
        ]
        A3n_xNy_c = A3n_xNy[0] + 1j * A3n_xNy[1]
        return (1. / 6) * A3n_xNy_c

    def a2_contract(self, k, A, D, nocc, norb):
        """
        Contracts the generalized dipole gradient tensor of rank 2 with a
        second-order response matrix. A[2]N1 = -(1 / 2)[[k1,X],D.T]

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

        An_x = self.commut(self.commut(k.T, A), D.T)
        A2n_x = [
            LinearSolver.lrmat2vec(An_x.real, nocc, norb),
            LinearSolver.lrmat2vec(An_x.imag, nocc, norb)
        ]
        A2n_x_c = A2n_x[0] + 1j * A2n_x[1]
        return -(1. / 2) * A2n_x_c

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

    def phi(self, kA, kB, Fb, F0):
        """
        Returns a matrix used for the E[3] contraction

        :param kA:
            First-order or Second-order response matrix
        :param kB:
            First-order or Second-order response matrix
        :param Fa:
            First-order or Second-order perturbed Fock matrix
        :param Fb:
            First-order or Second-order perturbed Fock matrix
        :param F0:
            SCF Fock matrix

        :return:
            Returns a matrix
        """

        return self.commut(kA, self.commut(kB, F0) + 3 * Fb)

    def anti_sym(self, vec):
        """
        Returns an antisymetrized vector

        :param vec:
            The vector to be anti-symetrized

        :return:
            An antisymetrized vector
        """

        # TODO: look into "anti_sym"

        if vec.ndim == 1:
            new_vec = np.zeros_like(vec)
            half_len = vec.shape[0] // 2
            new_vec[:half_len] = vec[:half_len]
            new_vec[half_len:] = -vec[half_len:]
            return new_vec

        return None

    def get_fock_r(self, mo, D, molecule, ao_basis, rank):
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
        :param rank:
            ...

        :return:
            A list of Fock matrices
        """

        # TODO: look into "rank"

        # TODO: move AO-to-MO transformation into get_two_el_fock_mod_r

        if rank == 0:
            # computes the unperturbed Fock matrix from the SCF Density and
            # adds the one electron part
            fa = self.get_two_el_fock_mod_r(molecule, ao_basis, D)
            h = self.get_one_el_hamiltonian(molecule, ao_basis)

            if self.rank == mpi_master():
                fa = 0.5 * fa[0] + h
                return self.ao2mo(mo, fa)
            else:
                return None

        elif rank == 1:
            # computes complex Fock matrices (only two-eletron parts 2J-K)
            if self.rank == mpi_master():
                D_total = []
                for da in D:
                    D_total.append(da.real)
                    D_total.append(da.imag)
            else:
                D_total = None

            f_total = self.get_two_el_fock_mod_r(molecule, ao_basis, D_total)

            ff = []

            if self.rank == mpi_master():
                for i in range(len(f_total) // 2):
                    ff.append(
                        self.ao2mo(
                            mo,
                            0.5 * f_total[2 * i] + 0.5j * f_total[2 * i + 1]))
                return ff
            else:
                return None

        elif rank == 3:
            # computes real Fock Matrices (only two-eletron parts 2J-K)
            f_total = self.get_two_el_fock_mod_r(molecule, ao_basis, D)

            ff = []

            if self.rank == mpi_master():
                for i in range(len(f_total)):
                    ff.append(self.ao2mo(mo, 0.5j * f_total[i]))

                return ff
            else:
                return None

        else:
            # computes imaginary Fock matrices (only two-eletron parts 2J-K)
            f_total = self.get_two_el_fock_mod_r(molecule, ao_basis, D)

            ff = []

            if self.rank == mpi_master():
                for i in range(len(f_total)):
                    ff.append(self.ao2mo(mo, 0.5 * f_total[i]))

                return ff
            else:
                return None

    def get_two_el_fock_mod_r(self, molecule, ao_basis, dabs):
        """
        Returns the two-electron part of the Fock matix 2J-K

        :param molecule:
            The molecule
        :param ao_basis:
            The AO basis set
        :param dabs:
            A list of densitiy matrices

        :return:
            A tuple containing the two-electron part of the Fock matix
        """

        # TODO: make a function to determine batch size

        eri_driver = ElectronRepulsionIntegralsDriver(self.comm)
        screening = eri_driver.compute(get_qq_scheme(self.qq_type),
                                       self.eri_thresh, molecule, ao_basis)

        # determine number of batches

        num_batches = 0

        total_mem = psutil.virtual_memory().total
        total_mem_list = self.comm.gather(total_mem, root=mpi_master())

        if self.rank == mpi_master():
            n_ao = dabs[0].shape[0]
            n_total = len(dabs)

            # check if master node has larger memory
            mem_adjust = 0.0
            if total_mem > min(total_mem_list):
                mem_adjust = total_mem - min(total_mem_list)

            # computes maximum batch size from available memory
            avail_mem = psutil.virtual_memory().available - mem_adjust
            mem_per_mat = n_ao**2 * ctypes.sizeof(ctypes.c_double)
            nthreads = int(os.environ['OMP_NUM_THREADS'])
            max_batch_size = int(avail_mem / mem_per_mat / (0.625 * nthreads))
            max_batch_size = max(1, max_batch_size)

            batch_size = self.batch_size
            if batch_size is None:
                batch_size = min(100, n_total, max_batch_size)

            # get number of batches
            num_batches = n_total // batch_size
            if n_total % batch_size != 0:
                num_batches += 1

        num_batches = self.comm.bcast(num_batches, root=mpi_master())

        # go through batches

        fabs = []

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

            dts = []
            if self.rank == mpi_master():
                batch_start = batch_size * batch_ind
                batch_end = min(batch_start + batch_size, n_total)
                for dab in dabs[batch_start:batch_end]:
                    dt = 2 * dab
                    dts.append(dt)
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
                for i in range(fock.number_of_fock_matrices()):
                    fabs.append(fock.to_numpy(i).T)

        self.ostream.print_blank()

        if self.rank == mpi_master():
            return tuple(fabs)
        else:
            return None

    def get_one_el_hamiltonian(self, molecule, ao_basis):
        """
        Returns the one electron part of the Fock matrix

        :param molecule:
            The molecule
        :param ao_basis:
            The AO basis set

        :return:
            The one electron part of the Fock matrix
        """

        kinetic_driver = KineticEnergyIntegralsDriver(self.comm)
        potential_driver = NuclearPotentialIntegralsDriver(self.comm)

        T = kinetic_driver.compute(molecule, ao_basis).to_numpy()
        V = potential_driver.compute(molecule, ao_basis).to_numpy()

        return T - V

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
