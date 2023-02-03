#
#                           VELOXCHEM 1.0-RC3
#         ----------------------------------------------------
#                     An Electronic Structure Code
#
#  Copyright © 2018-2022 by VeloxChem developers. All rights reserved.
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
import numpy as np
import time

from .veloxchemlib import ElectricDipoleIntegralsDriver
from .veloxchemlib import mpi_master, hartree_in_wavenumbers
from .profiler import Profiler
from .cppsolver import ComplexResponse
from .linearsolver import LinearSolver
from .nonlinearsolver import NonlinearSolver
from .distributedarray import DistributedArray
from .errorhandler import assert_msg_critical


class TpaDriver(NonlinearSolver):
    """
    Implements the isotropic cubic response driver for two-photon absorption
    (TPA)

    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.

    Instance variables
        - is_converged: The flag for convergence.
        - frequencies: The frequencies.
        - comp: The list of all the gamma tensor components
        - damping: The damping parameter.
        - lindep_thresh: The threshold for removing linear dependence in the
          trial vectors.
        - conv_thresh: The convergence threshold for the solver.
        - max_iter: The maximum number of solver iterations.
    """

    def __init__(self, comm, ostream):
        """
        Initializes the isotropic cubic response driver for two-photon
        absorption (TPA)
        """

        super().__init__(comm, ostream)

        # cpp settings
        self.frequencies = (0,)
        self.comp = None
        self.damping = 1000.0 / hartree_in_wavenumbers()
        self.lindep_thresh = 1.0e-10
        self.conv_thresh = 1.0e-4
        self.max_iter = 50

        # input keywords
        self._input_keywords['response'].update({
            'frequencies': ('seq_range', 'frequencies'),
            'damping': ('float', 'damping parameter'),
        })

    def update_settings(self, rsp_dict, method_dict=None):
        """
        Updates response and method settings in TPA driver

        :param rsp_dict:
            The dictionary of response input.
        :param method_dict:
            The dictionary of method settings.
        """

        if method_dict is None:
            method_dict = {}

        super().update_settings(rsp_dict, method_dict)

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

        # double check SCF information
        self._check_scf_results(scf_tensors)

        # check dft setup
        self._dft_sanity_check_nonlinrsp()

        profiler = Profiler({
            'timing': self.timing,
            'profiling': self.profiling,
            'memory_profiling': self.memory_profiling,
            'memory_tracing': self.memory_tracing,
        })

        if self.rank == mpi_master():
            self.print_header()

        start_time = time.time()

        dft_dict = self._init_dft(molecule, scf_tensors)

        # sanity check
        nalpha = molecule.number_of_alpha_electrons()
        nbeta = molecule.number_of_beta_electrons()
        assert_msg_critical(
            nalpha == nbeta,
            'TPA Driver: not implemented for unrestricted case')

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
        b_grad = linear_solver.get_complex_prop_grad(operator, component,
                                                     molecule, ao_basis,
                                                     scf_tensors)

        # This is a workaround for the sqrt(2) factor in the property gradient
        if self.rank == mpi_master():
            inv_sqrt_2 = 1.0 / np.sqrt(2.0)
            b_grad = list(b_grad)
            for ind in range(len(b_grad)):
                b_grad[ind] *= inv_sqrt_2

        # Storing the dipole integral matrices used for the X[3],X[2],A[3] and
        # A[2] contractions in MO basis
        if self.rank == mpi_master():
            v_grad = {(op, w): v for op, v in zip(component, b_grad)
                      for w in self.frequencies}
            X = {
                'x': 2 * self.ao2mo(mo, dipole_mats.x_to_numpy()),
                'y': 2 * self.ao2mo(mo, dipole_mats.y_to_numpy()),
                'z': 2 * self.ao2mo(mo, dipole_mats.z_to_numpy())
            }
            self.comp = self.get_comp(self.frequencies)
        else:
            v_grad = None
            X = None
            self.comp = None
        X = self.comm.bcast(X, root=mpi_master())
        self.comp = self.comm.bcast(self.comp, root=mpi_master())

        # Computing the first-order response vectors (3 per frequency)
        Nb_drv = ComplexResponse(self.comm, self.ostream)

        cpp_keywords = {
            'frequencies', 'damping', 'lindep_thresh', 'conv_thresh',
            'max_iter', 'eri_thresh', 'qq_type', 'timing', 'memory_profiling',
            'batch_size', 'restart', 'xcfun', 'grid_level', 'potfile',
            'electric_field', 'program_end_time'
        }

        for key in cpp_keywords:
            setattr(Nb_drv, key, getattr(self, key))

        if self.checkpoint_file is not None:
            Nb_drv.checkpoint_file = str(
                Path(self.checkpoint_file).with_suffix('.tpa_1.h5'))

        Nb_results = Nb_drv.compute(molecule, ao_basis, scf_tensors, v_grad)

        self._is_converged = Nb_drv.is_converged

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
                                               profiler, dft_dict)

        valstr = '*** Time spent in TPA calculation: {:.2f} sec ***'.format(
            time.time() - start_time)
        self.ostream.print_header(valstr)
        self.ostream.print_blank()
        self.ostream.flush()

        profiler.end(self.ostream)

        return tpa_dict

    def compute_tpa_components(self, Focks, w, X, d_a_mo, kX, track,
                               scf_tensors, molecule, ao_basis, profiler,
                               dft_dict):
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
            mo = scf_tensors['C_alpha']
            F0 = np.linalg.multi_dot([mo.T, scf_tensors['F_alpha'], mo])
            norb = mo.shape[1]
        else:
            mo = None
            F0 = None
            norb = None
        F0 = self.comm.bcast(F0, root=mpi_master())
        norb = self.comm.bcast(norb, root=mpi_master())

        nocc = molecule.number_of_alpha_electrons()

        # computing all compounded first-order densities
        if self.rank == mpi_master():
            density_list1, density_list2, density_list3 = self.get_densities(
                w, kX, mo, nocc)
        else:
            density_list1 = None
            density_list2 = None
            density_list3 = None

        profiler.check_memory_usage('1st densities')

        fock_profiler = Profiler({
            'timing': self.timing,
            'profiling': self.profiling,
            'memory_profiling': self.memory_profiling,
            'memory_tracing': self.memory_tracing,
        })

        #  computing the compounded first-order Fock matrices
        fock_dict = self.get_fock_dict(w, density_list1, density_list2,
                                       density_list3, F0, mo, molecule,
                                       ao_basis, dft_dict, fock_profiler)

        fock_profiler.end(self.ostream)

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
            density_list_two1, density_list_two2 = self.get_densities_II(
                w, kX, kXY_dict, mo, nocc)
        else:
            density_list_two1 = None
            density_list_two2 = None

        profiler.check_memory_usage('2nd densities')

        fock_profiler_two = Profiler({
            'timing': self.timing,
            'profiling': self.profiling,
            'memory_profiling': self.memory_profiling,
            'memory_tracing': self.memory_tracing,
        })

        # computing the remaning second-order Fock matrices from the
        # second-order densities
        fock_dict_two = self.get_fock_dict_II(w, density_list_two1,
                                              density_list_two2, mo, molecule,
                                              ao_basis, dft_dict,
                                              fock_profiler_two)

        fock_profiler_two.end(self.ostream)

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
                'w': self.frequencies
            })

        profiler.check_memory_usage('End of TPA')

        return result

    def get_densities(self, wi, kX, mo, nocc):
        """
        Computes the compounded densities needed for the compounded Fock
        matrices F^{σ},F^{λ+τ},F^{σλτ} used for the isotropic cubic response
        function

        :param wi:
            A list of the frequencies
        :param kX:
            A dictonary with all the first-order response matrices
        :param mo:
            A matrix containing the MO coefficents
        :param nocc:
            Number of alpha electrons

        :return:
            A list of tranformed compounded densities
        """

        return None

    def get_fock_dict(self, wi, density_list, F0, mo, molecule, ao_basis,
                      dft_dict, profiler):
        """
        Computes the compounded Fock matrices F^{σ},F^{λ+τ},F^{σλτ} used for the
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

    def get_densities_II(self, wi, kX, kXY, mo, nocc):
        """
        Computes the compounded densities needed for the compounded
        second-order Fock matrices used for the isotropic cubic response
        function

        :param wi:
            A list of the frequencies
        :param kX:
            A dictonary with all the first-order response matrices
        :param kXY:
            A dict of the two index response matrices
        :param mo:
            A matrix containing the MO coefficents
        :param nocc:
            Number of alpha electrons

        :return:
            A list of tranformed compounded densities
        """

        return None

    def get_fock_dict_II(self, wi, density_list, mo, molecule, ao_basis,
                         dft_dict, profiler):
        """
        Computes the compounded second-order Fock matrices used for the
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

            na_x = self.complex_lrmat2vec(ka_x, nocc, norb)
            na_y = self.complex_lrmat2vec(ka_y, nocc, norb)
            na_z = self.complex_lrmat2vec(ka_z, nocc, norb)

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

        Na = self.complex_lrmat2vec(ka_na, nocc, norb)

        if inp_dict['flag'] == 'CD':
            kcd = inp_dict['kcd']
            kb = inp_dict['kb']
            B = inp_dict['B']

            Ncd = self.complex_lrmat2vec(kcd, nocc, norb)
            Nb = self.complex_lrmat2vec(kb, nocc, norb)

            na_x2_nyz += np.dot(Na.T, self._x2_contract(kcd, B, da, nocc, norb))
            nx_a2_nyz += np.dot(self._a2_contract(kb, A, da, nocc, norb), Ncd)
            nx_a2_nyz += np.dot(self._a2_contract(kcd, A, da, nocc, norb), Nb)

        elif inp_dict['flag'] == 'BD':
            kbd = inp_dict['kbd']
            kc = -inp_dict['kc_kb'].T.conj()  # gets kc from kb
            C = inp_dict['C']

            Nbd = self.complex_lrmat2vec(kbd, nocc, norb)
            Nc = self.complex_lrmat2vec(kc, nocc, norb)

            na_x2_nyz += np.dot(Na.T, self._x2_contract(kbd, C, da, nocc, norb))
            nx_a2_nyz += np.dot(self._a2_contract(kc, A, da, nocc, norb), Nbd)
            nx_a2_nyz += np.dot(self._a2_contract(kbd, A, da, nocc, norb), Nc)

        return {
            'key': (w, -w, w),
            'x2': -(1. / 15) * na_x2_nyz,
            'a2': -(1. / 15) * nx_a2_nyz,
        }

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

    def _print_component(self, label, freq, value, width):
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
