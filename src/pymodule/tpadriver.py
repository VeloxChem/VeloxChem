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
from pathlib import Path
import numpy as np
import time

from .veloxchemlib import ElectronRepulsionIntegralsDriver
from .veloxchemlib import ElectricDipoleIntegralsDriver
from .veloxchemlib import denmat, fockmat
from .veloxchemlib import mpi_master, hartree_in_wavenumbers
from .qqscheme import get_qq_scheme
from .profiler import Profiler
from .cppsolver import ComplexResponse
from .linearsolver import LinearSolver
from .aofockmatrix import AOFockMatrix
from .aodensitymatrix import AODensityMatrix
from .distributedarray import DistributedArray
from .errorhandler import assert_msg_critical
from .inputparser import parse_input
from .batchsize import get_batch_size
from .batchsize import get_number_of_batches


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
        - program_start_time: The start time of the program.
        - maximum_hours: The timelimit in hours.
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
        self.damping = 1000.0 / hartree_in_wavenumbers()
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
        Updates response and method settings in TPA driver

        :param rsp_dict:
            The dictionary of response dict.
        :param method_dict:
            The dictionary of method settings.
        """

        if method_dict is None:
            method_dict = {}

        rsp_keywords = {
            'frequencies': 'seq_range',
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
            errmsg = 'TpaDriver: The \'xcfun\' keyword is not supported in TPA '
            errmsg += 'calculation.'
            if self.rank == mpi_master():
                assert_msg_critical(False, errmsg)

        if 'potfile' in method_dict:
            errmsg = 'TpaDriver: The \'potfile\' keyword is not supported in '
            errmsg += 'TPA calculation.'
            if self.rank == mpi_master():
                assert_msg_critical(False, errmsg)

        if 'electric_field' in method_dict:
            errmsg = 'TpaDriver: The \'electric field\' keyword is not '
            errmsg += 'supported in TPA calculation.'
            if self.rank == mpi_master():
                assert_msg_critical(False, errmsg)

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
            da = scf_tensors['D_alpha']
            mo = scf_tensors['C_alpha']
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
        component = 'xyz'

        linear_solver = LinearSolver(self.comm, self.ostream)
        b_rhs = linear_solver.get_complex_prop_grad(operator, component,
                                                    molecule, ao_basis,
                                                    scf_tensors)

        # This is a workaround for the sqrt(2) factor in the property gradient
        if self.rank == mpi_master():
            inv_sqrt_2 = 1.0 / np.sqrt(2.0)
            b_rhs = list(b_rhs)
            for ind in range(len(b_rhs)):
                b_rhs[ind] *= inv_sqrt_2

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
        X = self.comm.bcast(X, root=mpi_master())
        self.comp = self.comm.bcast(self.comp, root=mpi_master())

        # Computing the first-order response vectors (3 per frequency)
        Nb_drv = ComplexResponse(self.comm, self.ostream)

        cpp_keywords = {
            'frequencies', 'damping', 'lindep_thresh', 'conv_thresh',
            'max_iter', 'eri_thresh', 'qq_type', 'timing', 'memory_profiling',
            'batch_size', 'restart', 'program_start_time', 'maximum_hours'
        }

        for key in cpp_keywords:
            setattr(Nb_drv, key, getattr(self, key))

        if self.checkpoint_file is not None:
            Nb_drv.checkpoint_file = str(
                Path(self.checkpoint_file).with_suffix('.tpa_1.h5'))

        Nb_results = Nb_drv.compute(molecule, ao_basis, scf_tensors, v1)

        kX = {}
        Focks = {}

        kX['Nb'] = Nb_results['kappas']
        Focks['Fb'] = Nb_results['focks']
        Focks['Fd'] = {}

        # The first-order response vectors with negative frequency are
        # obtained from the first-order response vectors with positive
        # frequency by using flip_zy, see article.

        for (op, w) in kX['Nb']:

            # The first-order Fock matrices with positive and negative
            # frequencies are each other complex conjugates

            Fb_op_w = Focks['Fb'][(op, w)].get_full_vector()

            if self.rank == mpi_master():
                Fd_op_w = Fb_op_w.reshape(norb, norb).T.conj().reshape(-1)
            else:
                Fd_op_w = None

            Focks['Fd'][(op, w)] = DistributedArray(Fd_op_w, self.comm)

        # For cubic-response with all operators being the dipole μ Nb=Na=Nd
        # Likewise, Fb=Fd

        Focks['Fb'].update(Focks['Fd'])

        kX['Na'] = kX['Nb']
        kX['Nd'] = kX['Nb']

        profiler.check_memory_usage('1st CPP')

        # Computing the third-order gradient and also the contractions of
        # A[3] and A[2] which formally are not part of the third-order gradient
        # but which are used for the cubic response function

        tpa_dict = self.compute_tpa_components(Focks, self.frequencies, X,
                                               d_a_mo, kX, self.comp,
                                               scf_tensors, molecule, ao_basis,
                                               profiler)

        valstr = '*** Time spent in TPA calculation: {:.2f} sec ***'.format(
            time.time() - start_time)
        self.ostream.print_header(valstr)
        self.ostream.print_blank()
        self.ostream.flush()

        profiler.end(self.ostream)

        self.is_converged = True

        return tpa_dict

    def compute_tpa_components(self, Focks, w, X, d_a_mo, kX, track,
                               scf_tensors, molecule, ao_basis, profiler):
        """
        Computes all the relevent terms to third-order isotropic gradient

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
            D0 = scf_tensors['D_alpha']
            mo = scf_tensors['C_alpha']
            F0 = np.linalg.multi_dot([mo.T, scf_tensors['F_alpha'], mo])
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
            density_list = self.get_densities(w, kX, S, D0, mo)
        else:
            density_list = None

        profiler.check_memory_usage('1st densities')

        #  computing the compounded first-order Fock matrices
        fock_dict = self.get_fock_dict(w, density_list, F0, mo, molecule,
                                       ao_basis)

        profiler.check_memory_usage('1st Focks')

        fock_dict.update(Focks)
        e4_dict = self.get_e4(w, kX, fock_dict, nocc, norb)

        profiler.check_memory_usage('E[4]')

        # computing all the compounded second-order response vectors and
        # extracting some of the second-order Fock matrices from the subspace
        (kXY_dict, Focks_xy) = self.get_Nxy(w, d_a_mo, X, fock_dict, kX, nocc,
                                            norb, molecule, ao_basis,
                                            scf_tensors)

        profiler.check_memory_usage('2nd CPP')

        # computing all second-order compounded densities based on the
        # second-order response vectors
        if self.rank == mpi_master():
            density_list_two = self.get_densities_II(w, kX, kXY_dict, S, D0, mo)
        else:
            density_list_two = None

        profiler.check_memory_usage('2nd densities')

        # computing the remaning second-order Fock matrices from the
        # second-order densities
        fock_dict_two = self.get_fock_dict_II(w, density_list_two, mo, molecule,
                                              ao_basis)

        profiler.check_memory_usage('2nd Focks')

        # Adding the Fock matrices extracted from the second-order response
        # vector subspace to the fock_dict's.
        fock_dict_two.update(Focks_xy)

        # computing the compounded E[3] contractions for the isotropic
        # cubic response function
        e3_dict = self.get_e3(w, kX, kXY_dict, fock_dict, fock_dict_two, nocc,
                              norb)

        profiler.check_memory_usage('E[3]')

        # computing the X[3],A[3],X[2],A[2] contractions for the isotropic
        # cubic response function
        other_dict = self.get_other_terms(w, track, X, kX, kXY_dict, d_a_mo,
                                          nocc, norb)

        profiler.check_memory_usage('X[3],A[3],X[2],A[2]')

        # Combining all the terms to evaluate the iso-tropic cubic response
        # function. For TPA Full and reduced, see article

        t4_dict = self.get_t4(self.frequencies, e4_dict, kX, self.comp, d_a_mo,
                              nocc, norb)
        if self.rank == mpi_master():
            t3_dict = self.get_t3(self.frequencies, e3_dict, kX, self.comp,
                                  nocc, norb)

        profiler.check_memory_usage('T[4],T[3]')

        result = {}

        if self.rank == mpi_master():
            gamma = {}

            for w in self.frequencies:
                sum_val = t3_dict[(w, -w, w)]
                if t4_dict is not None:
                    sum_val += t4_dict[(w, -w, w)]
                for key, val in other_dict.items():
                    sum_val += val[(w, -w, w)]
                gamma[(w, -w, w)] = sum_val

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

    def get_densities(self, wi, kX, S, D0, mo):
        """
        Computes the compounded densities needed for the compounded Fock
        matrics F^{σ},F^{λ+τ},F^{σλτ} used for the isotropic cubic response
        function

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

        return None

    def get_fock_dict(self, wi, density_list, F0, mo, molecule, ao_basis):
        """
        Computes the compounded Fock matrics F^{σ},F^{λ+τ},F^{σλτ} used for the
        isotropic cubic response function

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

        return None

    def get_e4(self, wi, kX, fo, nocc, norb):
        """
        Contracts E[4]NxNyNz for the isotropic cubic response function. Takes
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

    def get_Nxy(self, w, d_a_mo, X, fock_dict, kX, nocc, norb, molecule,
                ao_basis, scf_tensors):
        """
        Computes all the second-order response vectors needed for the isotropic
        cubic response computation

        :param w:
            A list of all the frequencies
        :param d_a_mo:
            The density matrix in MO basis
        :param X:
            Dipole integrals
        :param fock_dict:
            A dictonary containing all the Fock matricies
        :param kX:
            A dictonary containg all the response matricies
        :param nocc:
            The number of occupied orbitals
        :param norb:
            The number of total orbitals
        :param molecule:
            The molecule.
        :param basis:
            The AO basis.
        :param scf_tensors:
            The dictionary of tensors from converged SCF wavefunction.

        :return:
            A dictonary of Fock matrices from the subspace,second-order
            response vectors and second-order response matrices
        """

        return None

    def get_densities_II(self, wi, kX, kXY, S, D0, mo):
        """
        Computes the compounded densities needed for the compounded
        second-order Fock matrics used for the isotropic cubic response
        function

        :param wi:
            A list of the frequencies
        :param kX:
            A dictonary with all the first-order response matrices
        :param kXY:
            A dict of the two index response matrices
        :param S:
            The overlap matrix
        :param D0:
            The SCF density matrix in AO basis
        :param mo:
            A matrix containing the MO coefficents

        :return:
            A list of tranformed compounded densities
        """

        return None

    def get_fock_dict_II(self, wi, density_list, mo, molecule, ao_basis):
        """
        Computes the compounded second-order Fock matrics used for the
        isotropic cubic response function

        :param wi:
            A list of the frequencies
        :param density_list:
            A list of tranformed compounded densities
        :param mo:
            A matrix containing the MO coefficents
        :param molecule:
            The molecule
        :param ao_basis:
            The AO basis set

        :return:
            A dictonary of compounded second-order Fock-matrices
        """

        return None

    def get_e3(self, wi, kX, kXY, fo, fo2, nocc, norb):
        """
        Contracts E[3]

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

        return None

    def get_other_terms(self, wi, track, X, kX, kXY, da, nocc, norb):
        """
        Computes the terms involving X[3],A[3],X[2],A[2] in the isotropic cubic
        response function

        :param wi:
            A list containing all the frequencies
        :param track:
            A list that contains information about what γ components that are
            to be computed and which freqs
        :param X:
            A dictonray with all the property integral matricies
        :param kX:
            A dictonary with all the respone matricies
        :param kXY:
            A dictonary containing all the two-index response matricies
        :param da:
            The SCF density matrix in MO basis
        :param nocc:
            The number of occupied orbitals
        :param norb:
            The total number of orbitals

        :return:
            A dictonary of final X[2],A[2] contraction values
        """

        return None

    def get_t4(self, wi, e4_dict, kX, track, da, nocc, norb):
        """
        Computes the contraction of the E[4] tensor with that of the S[4] and
        R[4] tensors to return the contraction of T[4] as a dictonary of
        vectors. T[4]NxNyNz = (E^[4]-ω_1S^[4]-ω_1S^[4]-ω_3S^[4]-γiR^[4])

        :param wi:
            A list of all the freqs
        :param e4_dict:
            A dictonary of all the E[4] contraction
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

    def get_t3(self, freqs, e3_dict, kX, track, nocc, norb):
        """
        Computes the T[3] contraction, for HF S[3] = 0, R[3] = 0 such that
        the T[3] contraction for the isotropic cubic response function in terms
        of compounded Fock matrices is given as:

                             [(ζ_{α}^{σσ} + ζ_{α}^{λλ+ττ} + f_{α}^{λσ,τ})_is]
        t3term = Σ_{α} N_{α} [(ζ_{α}^{σσ} + ζ_{α}^{λλ+ττ} + f_{α}^{λσ,τ})_si]

        For more details see article

        :param freqs:
            List of frequencies of the pertubations
        :param e3_dict:
            A dictonary that contains the contractions of E[3]
        :param kX:
            A dictonray containng all the response matricies
        :param track:
            A list containing information about what tensor components that are
            being computed
        :param nocc:
            The number of occupied orbitals
        :param norb:
            The total number of orbitals

        :return:
            A dictonary of the final values for the NaT[3]NxNyz contractions
        """

        t3_term = {}

        for i in range(len(freqs)):
            w = float(track[i * (len(track) // len(freqs))].split(",")[1])

            ka_x = kX['Na'][('x', w)]
            ka_y = kX['Na'][('y', w)]
            ka_z = kX['Na'][('z', w)]

            na_x = (LinearSolver.lrmat2vec(ka_x.real, nocc, norb) +
                    1j * LinearSolver.lrmat2vec(ka_x.imag, nocc, norb))
            na_y = (LinearSolver.lrmat2vec(ka_y.real, nocc, norb) +
                    1j * LinearSolver.lrmat2vec(ka_y.imag, nocc, norb))
            na_z = (LinearSolver.lrmat2vec(ka_z.real, nocc, norb) +
                    1j * LinearSolver.lrmat2vec(ka_z.imag, nocc, norb))

            t3term = (np.dot(na_x, e3_dict['f_iso_x'][w]) +
                      np.dot(na_y, e3_dict['f_iso_y'][w]) +
                      np.dot(na_z, e3_dict['f_iso_z'][w]))

            t3_term[(w, -w, w)] = 1. / 15 * t3term

        return t3_term

    def get_x3_a3(self, inp_dict, da, nocc, norb):
        """
        Computes X[3] and A[3] contributions.

        :param inp_dict:
            A dictionary containing input data for computing X[3] and A[3].
        :param da:
            The SCF density matrix in MO basis
        :param nocc:
            The number of occupied orbitals
        :param norb:
            The total number of orbitals

        :return:
            A dictionary containing frequencies, X[3] and A[3].
        """

        return None

    def get_x2_a2(self, inp_dict, da, nocc, norb):
        """
        Computes X[2] and A[2] contributions.

        :param inp_dict:
            A dictionary containing input data for computing X[2] and A[2].
        :param da:
            The SCF density matrix in MO basis
        :param nocc:
            The number of occupied orbitals
        :param norb:
            The total number of orbitals

        :return:
            A dictionary containing frequencies, X[2] and A[2].
        """

        na_x2_nyz = 0.0
        nx_a2_nyz = 0.0

        w = inp_dict['freq']
        ka_na = inp_dict['Na_ka']
        A = inp_dict['A']

        Na = (LinearSolver.lrmat2vec(ka_na.real, nocc, norb) +
              1j * LinearSolver.lrmat2vec(ka_na.imag, nocc, norb))

        if inp_dict['flag'] == 'CD':
            kcd = inp_dict['kcd']
            kb = inp_dict['kb']
            B = inp_dict['B']

            Ncd = (LinearSolver.lrmat2vec(kcd.real, nocc, norb) +
                   1j * LinearSolver.lrmat2vec(kcd.imag, nocc, norb))
            Nb = (LinearSolver.lrmat2vec(kb.real, nocc, norb) +
                  1j * LinearSolver.lrmat2vec(kb.imag, nocc, norb))

            na_x2_nyz += np.dot(Na.T, self.x2_contract(kcd, B, da, nocc, norb))
            nx_a2_nyz += np.dot(self.a2_contract(kb, A, da, nocc, norb), Ncd)
            nx_a2_nyz += np.dot(self.a2_contract(kcd, A, da, nocc, norb), Nb)

        elif inp_dict['flag'] == 'BD':
            kbd = inp_dict['kbd']
            kc = -inp_dict['kc_kb'].T.conj()  # gets kc from kb
            C = inp_dict['C']

            Nbd = (LinearSolver.lrmat2vec(kbd.real, nocc, norb) +
                   1j * LinearSolver.lrmat2vec(kbd.imag, nocc, norb))
            Nc = (LinearSolver.lrmat2vec(kc.real, nocc, norb) +
                  1j * LinearSolver.lrmat2vec(kc.imag, nocc, norb))

            na_x2_nyz += np.dot(Na.T, self.x2_contract(kbd, C, da, nocc, norb))
            nx_a2_nyz += np.dot(self.a2_contract(kc, A, da, nocc, norb), Nbd)
            nx_a2_nyz += np.dot(self.a2_contract(kbd, A, da, nocc, norb), Nc)

        return {
            'key': (w, -w, w),
            'x2': -(1. / 15) * na_x2_nyz,
            'a2': -(1. / 15) * nx_a2_nyz,
        }

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

        XNx = self.commut(self.commut(k, X), D.T)
        X2Nx_c = (LinearSolver.lrmat2vec(XNx.real, nocc, norb) +
                  1j * LinearSolver.lrmat2vec(XNx.imag, nocc, norb))
        return X2Nx_c

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

        if vec.ndim == 1:
            new_vec = np.zeros_like(vec)
            half_len = vec.shape[0] // 2
            new_vec[:half_len] = vec[:half_len]
            new_vec[half_len:] = -vec[half_len:]
            return new_vec

        return None

    def get_fock_r(self, mo, D, molecule, ao_basis, fock_flag):
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

            f_total = self.get_two_el_fock_mod_r(mo, molecule, ao_basis,
                                                 D_total)

            nrows = f_total.data.shape[0]
            half_ncols = f_total.data.shape[1] // 2
            ff_data = np.zeros((nrows, half_ncols), dtype=np.complex128)
            for i in range(half_ncols):
                ff_data[:, i] = (f_total.data[:, 2 * i] +
                                 1j * f_total.data[:, 2 * i + 1])
            return DistributedArray(ff_data, self.comm, distribute=False)

        elif fock_flag == 'real':
            return self.get_two_el_fock_mod_r(mo, molecule, ao_basis, D)

        else:
            return None

    def get_two_el_fock_mod_r(self, mo, molecule, ao_basis, dabs):
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
