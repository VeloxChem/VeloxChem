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

from pathlib import Path
import numpy as np
import time

from .oneeints import compute_electric_dipole_integrals
from .veloxchemlib import (mpi_master, bohr_in_angstrom, hartree_in_ev,
                           hartree_in_inverse_nm, hartree_in_wavenumber,
                           fine_structure_constant,
                           speed_of_light_in_vacuum_in_SI)
from .profiler import Profiler
from .cppsolver import ComplexResponse
from .linearsolver import LinearSolver
from .nonlinearsolver import NonlinearSolver
from .distributedarray import DistributedArray
from .sanitychecks import (molecule_sanity_check, scf_results_sanity_check,
                           dft_sanity_check)
from .errorhandler import assert_msg_critical
from .checkpoint import check_distributed_focks
from .checkpoint import read_distributed_focks
from .checkpoint import write_distributed_focks


class ThgRedDriver(NonlinearSolver):
    """
    Implements the isotropic cubic response driver for Third-harmonic gerneration
    (thg)

    # vlxtag: RHF, thg, CR

    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.

    Instance variables
        - is_converged: The flag for convergence.
        - frequencies: The frequencies.
        - comp: The list of all the gamma tensor components
        - damping: The damping parameter.
        - conv_thresh: The convergence threshold for the solver.
        - max_iter: The maximum number of solver iterations.
    """

    def __init__(self, comm, ostream):
        """
        Initializes the isotropic cubic response driver for two-photon
        absorption (thg)
        """

        super().__init__(comm, ostream)

        # cpp settings
        self.frequencies = (0,)
        self.comp = None
        self.damping = 1000.0 / hartree_in_wavenumber()

        # input keywords
        self._input_keywords['response'].update({
            'frequencies': ('seq_range', 'frequencies'),
            'damping': ('float', 'damping parameter'),
        })

    def update_settings(self, rsp_dict, method_dict=None):
        """
        Updates response and method settings in thg driver

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
            A dictonary containing the isotropic T[4], T[3], X[3], A[3], X[2],
            A[2] contractions and the isotropic cubic response functions for
            thg.
        """

        if self.norm_thresh is None:
            self.norm_thresh = self.conv_thresh * 1.0e-6
        if self.lindep_thresh is None:
            self.lindep_thresh = self.conv_thresh * 1.0e-6

        # check molecule
        molecule_sanity_check(molecule)

        # check SCF results
        scf_results_sanity_check(self, scf_tensors)

        # update checkpoint_file after scf_results_sanity_check
        if self.filename is not None and self.checkpoint_file is None:
            self.checkpoint_file = f'{self.filename}_rsp.h5'

        # check dft setup
        dft_sanity_check(self, 'compute', 'nonlinear')

        profiler = Profiler({
            'timing': self.timing,
            'profiling': self.profiling,
            'memory_profiling': self.memory_profiling,
            'memory_tracing': self.memory_tracing,
        })

        if self.rank == mpi_master():
            self._print_header('Third-Harmonic Generation Driver Setup')

        start_time = time.time()

        eri_dict = self._init_eri(molecule, ao_basis)

        dft_dict = self._init_dft(molecule, scf_tensors)

        # sanity check
        nalpha = molecule.number_of_alpha_electrons()
        nbeta = molecule.number_of_beta_electrons()
        assert_msg_critical(
            nalpha == nbeta,
            'thg Driver: not implemented for unrestricted case')

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
        if self.rank == mpi_master():
            mu_dipole_mats = compute_electric_dipole_integrals(
                molecule, ao_basis)
            # Note: nonliear response uses r instead of mu for dipole operator
            dipole_mats = (mu_dipole_mats[0] * (-1.0),
                           mu_dipole_mats[1] * (-1.0),
                           mu_dipole_mats[2] * (-1.0))
        else:
            dipole_mats = tuple()

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
                # Note: nonliear response uses r instead of mu for dipole operator
                if operator == 'dipole':
                    b_grad[ind] *= -1.0

        # Storing the dipole integral matrices used for the X[3],X[2],A[3] and
        # A[2] contractions in MO basis
        if self.rank == mpi_master():
            v_grad = {
                (op, w_): v
                for op, v in zip(component, b_grad)
                for w in self.frequencies
                for w_ in (w, 3 * w)
            }
            X = {
                'x': 2 * self.ao2mo(mo, dipole_mats[0]),
                'y': 2 * self.ao2mo(mo, dipole_mats[1]),
                'z': 2 * self.ao2mo(mo, dipole_mats[2])
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
            'frequencies', 'damping', 'norm_thresh', 'lindep_thresh',
            'conv_thresh', 'max_iter', 'eri_thresh', 'timing',
            'memory_profiling', 'batch_size', 'restart', 'xcfun', 'grid_level',
            'potfile', 'electric_field', 'program_end_time', '_debug',
            '_block_size_factor', 'ri_coulomb'
        }

        for key in cpp_keywords:
            setattr(Nb_drv, key, getattr(self, key))

        if self.checkpoint_file is not None:
            fpath = Path(self.checkpoint_file)
            fpath = fpath.with_name(fpath.stem)
            Nb_drv.checkpoint_file = str(fpath) + '_tpa_1.h5'

        Nb_results = Nb_drv.compute(molecule, ao_basis, scf_tensors, v_grad)

        self._is_converged = Nb_drv.is_converged

        Nx = Nb_results['solutions']
        Focks = {'Fb': Nb_results['focks']}

        # The first-order response vectors with negative frequency are
        # obtained from the first-order response vectors with positive
        # frequency by using flip_zy, see article.

        for (op, w) in Focks['Fb']:

            # The first-order Fock matrices with positive and negative
            # frequencies are each other complex conjugates

            Fb_op_w = Focks['Fb'][(op, w)].get_full_vector()

            if self.rank == mpi_master():
                Fd_op_w = Fb_op_w.reshape(norb, norb).T.conj().reshape(-1)
            else:
                Fd_op_w = None

            Focks['Fb'][(op, w)] = DistributedArray(Fd_op_w, self.comm)

        profiler.check_memory_usage('1st CPP')

        # Computing the third-order gradient and also the contractions of
        # A[3] and A[2] which formally are not part of the third-order gradient
        # but which are used for the cubic response function

        tpa_dict = self.compute_thg_components(Focks, self.frequencies, X,
                                               d_a_mo, Nx, self.comp,
                                               scf_tensors, molecule, ao_basis,
                                               profiler, eri_dict, dft_dict)

        valstr = '*** Time spent in thg calculation: {:.2f} sec ***'.format(
            time.time() - start_time)
        self.ostream.print_header(valstr)
        self.ostream.print_blank()
        self.ostream.flush()

        profiler.end(self.ostream)

        return tpa_dict

    def compute_thg_components(self, Focks, w, X, d_a_mo, Nx, track,
                               scf_tensors, molecule, ao_basis, profiler,
                               eri_dict, dft_dict):
        """
        Computes all the relevent terms to third-order isotropic gradient

        :param w:
            A list of all the frequencies
        :param X:
            A dictonary of matricies containing all the dipole integrals
        :param d_a_mo:
            The SCF density in MO basis
        :param Nx:
            A dictonary containing all the response vectors in distributed form
        :param track:
            A list that contains all the information about which γ components
            and at what freqs they are to be computed
        :param scf_tensors:
            The dictionary of tensors from converged SCF wavefunction.
        :param molecule:
            The molecule.
        :param basis:
            The AO basis.
        :param eri_dict:
            The dictionary containing ERI information
        :param dft_dict:
            The dictionary containing DFT information

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
        density_list1, density_list2, density_list3 = self.get_densities(w, Nx, mo, nocc, norb)

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
                                       ao_basis, eri_dict, dft_dict,
                                       fock_profiler)

        fock_profiler.end(self.ostream)

        profiler.check_memory_usage('1st Focks')

        fock_dict.update(Focks)
        e4_dict = self.get_e4(w, Nx, fock_dict, nocc, norb)

        profiler.check_memory_usage('E[4]')

        # computing all the compounded second-order response vectors and
        # extracting some of the second-order Fock matrices from the subspace
        (Nxy_dict, Focks_xy) = self.get_Nxy(w, d_a_mo, X, fock_dict, Nx, nocc,
                                            norb, molecule, ao_basis,
                                            scf_tensors)

        profiler.check_memory_usage('2nd CPP')

        # computing all second-order compounded densities based on the
        # second-order response vectors
        density_list_two1, density_list_two2 = self.get_densities_II(w, Nx, Nxy_dict, mo, nocc, norb)

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
                                              ao_basis, eri_dict, dft_dict,
                                              fock_profiler_two)

        fock_profiler_two.end(self.ostream)

        profiler.check_memory_usage('2nd Focks')

        # Adding the Fock matrices extracted from the second-order response
        # vector subspace to the fock_dict's.
        fock_dict_two.update(Focks_xy)

        # computing the compounded E[3] contractions for the isotropic
        # cubic response function
        e3_dict = self.get_e3(w, Nx, Nxy_dict, fock_dict, fock_dict_two, nocc,norb)

        profiler.check_memory_usage('E[3]')

        # computing the X[3],A[3],X[2],A[2] contractions for the isotropic
        # cubic response function
        other_dict = self.get_other_terms(w, track, X, Nx, Nxy_dict, d_a_mo, nocc, norb)

        profiler.check_memory_usage('X[3],A[3],X[2],A[2]')

        # Combining all the terms to evaluate the iso-tropic cubic response
        # function. For thg Full and reduced, see article

        t4_dict = self.get_t4(self.frequencies, e4_dict, Nx, self.comp, d_a_mo,nocc, norb)
        t3_dict = self.get_t3(self.frequencies, e3_dict, Nx, self.comp, nocc,norb)

        profiler.check_memory_usage('T[4],T[3]')

        ret_dict = {}

        if self.rank == mpi_master():
            THG = {}

            for w in self.frequencies:
                sum_val = t3_dict[(w, w, w)]
                if t4_dict is not None:
                    sum_val += t4_dict[(w, w, w)]
                for key, val in other_dict.items():
                    sum_val += val[(w, w, w)]
                THG[(w, w, w)] = sum_val.real

            ret_dict.update(other_dict)

            ret_dict.update({
                't4_dict': t4_dict,
                't3_dict': t3_dict,
                'other': other_dict,
                'THG': THG,
                'frequencies': list(self.frequencies),
            })


            self._print_results2(ret_dict)



            
        
        

        profiler.check_memory_usage('End of thg')

        return ret_dict

    def get_densities(self, wi, Nx, mo, nocc, norb):
        """
        Computes the compounded densities needed for the compounded Fock
        matrices F^{σ},F^{λ+τ},F^{σλτ} used for the isotropic cubic response
        function. Note: All densities are 1/3 of those in the paper, and all
        the Fock matrices are later scaled by 3.

        :param wi:
            A list of the frequencies
        :param Nx:
            A dictonary with all the first-order response vectors in
            distributed form
        :param mo:
            A matrix containing the MO coefficents
        :param nocc:
            Number of occupied orbitals
        :param norb:
            Number of orbitals

        :return:
            A list of tranformed compounded densities
        """

        distributed_density_1 = None
        distributed_density_2 = None
        distributed_density_3 = None

        one_third = 1.0 / 3.0

        for w in wi:

            nx = ComplexResponse.get_full_solution_vector(Nx[('x', w)])
            ny = ComplexResponse.get_full_solution_vector(Nx[('y', w)])
            nz = ComplexResponse.get_full_solution_vector(Nx[('z', w)])

            if self.rank == mpi_master():

                kx = self.complex_lrvec2mat(nx, nocc, norb)
                ky = self.complex_lrvec2mat(ny, nocc, norb)
                kz = self.complex_lrvec2mat(nz, nocc, norb)

                # create the first order single indexed densiteies #

                Dx = self.commut_mo_density(kx, nocc)
                Dy = self.commut_mo_density(ky, nocc)
                Dz = self.commut_mo_density(kz, nocc)

                # create the first order two indexed densities #


                # σ terms #

                Dxx = self.commut(kx, Dx)
                Dyy = self.commut(ky, Dy)
                Dzz = self.commut(kz, Dz)

                D_sig_xx = 6 * (3 * Dxx + Dyy + Dzz)
                D_sig_yy = 6 * (Dxx + 3 * Dyy + Dzz)
                D_sig_zz = 6 * (Dxx + Dyy + 3 * Dzz)
                D_sig_xy = 6 * (self.commut(ky, Dx) + self.commut(kx, Dy))
                D_sig_xz = 6 * (self.commut(kx, Dz) + self.commut(kz, Dx))
                D_sig_yz = 6 * (self.commut(ky, Dz) + self.commut(kz, Dy))

                # λ+τ terms #

                Dxx = self.commut(kx, Dx) + self.commut(kx, Dx)
                Dyy = self.commut(ky, Dy) + self.commut(ky, Dy)
                Dzz = self.commut(kz, Dz) + self.commut(kz, Dz)

                D_lamtau_xx = 2.0 * D_sig_xx
                D_lamtau_yy = 2.0 * D_sig_yy
                D_lamtau_zz = 2.0 * D_sig_zz
                D_lamtau_xy = 2.0 * D_sig_xy 
                D_lamtau_xz = 2.0 * D_sig_xz 
                D_lamtau_yz = 2.0 * D_sig_yz 

                # Create first order three indexed Densities #

                D_lam_sig_tau_x = (self.commut(kx, one_third * D_sig_xx) +
                                   self.commut(ky, one_third * D_sig_xy) +
                                   self.commut(kz, one_third * D_sig_xz))

                D_lam_sig_tau_x += (self.commut(kx, one_third * D_lamtau_xx) +
                                    self.commut(ky, one_third * D_lamtau_xy) +
                                    self.commut(kz, one_third * D_lamtau_xz))

                D_lam_sig_tau_y = (self.commut(kx, one_third * D_sig_xy) +
                                   self.commut(ky, one_third * D_sig_yy) +
                                   self.commut(kz, one_third * D_sig_yz))

                D_lam_sig_tau_y += (self.commut(kx, one_third * D_lamtau_xy) +
                                    self.commut(ky, one_third * D_lamtau_yy) +
                                    self.commut(kz, one_third * D_lamtau_yz))

                D_lam_sig_tau_z = (self.commut(kx, one_third * D_sig_xz) +
                                   self.commut(ky, one_third * D_sig_yz) +
                                   self.commut(kz, one_third * D_sig_zz))

                D_lam_sig_tau_z += (self.commut(kx, one_third * D_lamtau_xz) +
                                    self.commut(ky, one_third * D_lamtau_yz) +
                                    self.commut(kz, one_third * D_lamtau_zz))
                
                # density transformation from MO to AO basis

                Dx = np.linalg.multi_dot([mo, Dx, mo.T])
                Dy = np.linalg.multi_dot([mo, Dy, mo.T])
                Dz = np.linalg.multi_dot([mo, Dz, mo.T])

                D_sig_xx = np.linalg.multi_dot([mo, D_sig_xx, mo.T])
                D_sig_yy = np.linalg.multi_dot([mo, D_sig_yy, mo.T])
                D_sig_zz = np.linalg.multi_dot([mo, D_sig_zz, mo.T])
                D_sig_xy = np.linalg.multi_dot([mo, D_sig_xy, mo.T])
                D_sig_xz = np.linalg.multi_dot([mo, D_sig_xz, mo.T])
                D_sig_yz = np.linalg.multi_dot([mo, D_sig_yz, mo.T])

                D_lam_sig_tau_x = np.linalg.multi_dot([mo, D_lam_sig_tau_x, mo.T])
                D_lam_sig_tau_y = np.linalg.multi_dot([mo, D_lam_sig_tau_y, mo.T])
                D_lam_sig_tau_z = np.linalg.multi_dot([mo, D_lam_sig_tau_z, mo.T])

                # Reshape all the density matrices as vectors and combine them
                # all into one object

                dist_den_1_freq = np.hstack((
                    Dx.real.reshape(-1, 1),
                    Dy.real.reshape(-1, 1),
                    Dz.real.reshape(-1, 1),
                ))

                dist_den_2_freq = np.hstack((
                    D_sig_xx.real.reshape(-1, 1),
                    D_sig_yy.real.reshape(-1, 1),
                    D_sig_zz.real.reshape(-1, 1),
                    D_sig_xy.real.reshape(-1, 1),
                    D_sig_xz.real.reshape(-1, 1),
                    D_sig_yz.real.reshape(-1, 1),
                ))

                dist_den_3_freq = np.hstack((
                    D_lam_sig_tau_x.real.reshape(-1, 1),
                    D_lam_sig_tau_y.real.reshape(-1, 1),
                    D_lam_sig_tau_z.real.reshape(-1, 1),
                ))

            else:
                dist_den_1_freq = None
                dist_den_2_freq = None
                dist_den_3_freq = None

            dist_den_1_freq = DistributedArray(dist_den_1_freq, self.comm)
            dist_den_2_freq = DistributedArray(dist_den_2_freq, self.comm)
            dist_den_3_freq = DistributedArray(dist_den_3_freq, self.comm)

            if distributed_density_1 is None:
                distributed_density_1 = DistributedArray(dist_den_1_freq.data,
                                                         self.comm,
                                                         distribute=False)
            else:
                distributed_density_1.append(dist_den_1_freq, axis=1)

            if distributed_density_2 is None:
                distributed_density_2 = DistributedArray(dist_den_2_freq.data,
                                                         self.comm,
                                                         distribute=False)
            else:
                distributed_density_2.append(dist_den_2_freq, axis=1)

            if distributed_density_3 is None:
                distributed_density_3 = DistributedArray(dist_den_3_freq.data,
                                                         self.comm,
                                                         distribute=False)
            else:
                distributed_density_3.append(dist_den_3_freq, axis=1)

        return distributed_density_1, distributed_density_2, distributed_density_3


    def get_fock_dict(self, wi, density_list1, density_list2, density_list3,
                      F0_a, mo, molecule, ao_basis, eri_dict, dft_dict,
                      profiler):
        """
        Computes the compounded Fock matrices F^{σ},F^{λ+τ},F^{σλτ} used for the
        isotropic cubic response function

        :param wi:
            A list of the frequencies
        :param density_list:
            A list of tranformed compounded densities
        :param F0_a:
            The Fock matrix in MO basis
        :param mo:
            A matrix containing the MO coefficents
        :param molecule:
            The molecule
        :param ao_basis:
            The AO basis set
        :param eri_dict:
            The dictionary containing ERI information
        :param dft_dict:
            The dictionary containing DFT information

        :return:
            A dictonary of compounded first-order Fock-matrices
        """

        if self.rank == mpi_master():
            self._print_fock_header()

        # generate key-frequency pairs

        key_freq_pairs = []

        for wb in wi:
            for key in [
                    'f_sig_xx', 'f_sig_yy', 'f_sig_zz', 'f_sig_xy', 'f_sig_xz',
                    'f_sig_yz']:
                key_freq_pairs.append((key, wb))

        for wb in wi:
            for key in ['F123_x', 'F123_y', 'F123_z']:
                key_freq_pairs.append((key, wb))

        # examine checkpoint file for distributed Focks

        if self.checkpoint_file is not None:
            fpath = Path(self.checkpoint_file)
            fpath = fpath.with_name(fpath.stem)
            fock_file = str(fpath) + '_thg_red_fock_1_full.h5'
        else:
            fock_file = None

        if self.restart:
            if self.rank == mpi_master():
                self.restart = check_distributed_focks(fock_file,key_freq_pairs)
            self.restart = self.comm.bcast(self.restart, mpi_master())

        # read or compute distributed Focks

        if self.restart:
            dist_focks = read_distributed_focks(fock_file, self.comm,self.ostream)
        else:
            time_start_fock = time.time()

            if self._dft:
                dist_focks = self._comp_nlr_fock(mo, molecule, ao_basis,
                                                 'real', eri_dict,
                                                 dft_dict, density_list1,
                                                 density_list2, density_list3,
                                                 'thgred', profiler)
            else:
                density_list_23 = DistributedArray(density_list2.data,
                                                   self.comm,
                                                   distribute=False)
                density_list_23.append(density_list3, axis=1)
                dist_focks = self._comp_nlr_fock(mo, molecule, ao_basis,
                                                 'real', eri_dict,
                                                 None, None, None,
                                                 density_list_23, 'thgred',
                                                 profiler)

            self._print_fock_time(time.time() - time_start_fock)

            write_distributed_focks(fock_file, dist_focks, key_freq_pairs,
                                    self.comm, self.ostream)

        # assign distributed Focks to key-frequency pairs

        focks = {'F0': F0_a}

        for fock_index, (key, wb) in enumerate(key_freq_pairs):
            if key not in focks:
                focks[key] = {}
            focks[key][wb] = DistributedArray(dist_focks.data[:, fock_index],
                                              self.comm,
                                              distribute=False)

        return focks

    def get_e4(self, wi, Nx, fo, nocc, norb):
        """
        Contracts E[4]NxNyNz for the isotropic cubic response function. Takes
        the Fock matrices from fock_dict and contracts them with the response
        vectors.

        :param wi:
            A list of freqs
        :param Nx:
            A dict of the single index response vectors in distributed form
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

        f_iso_x = {}
        f_iso_y = {}
        f_iso_z = {}

        for w in wi:

            vec_pack = np.array([
                fo['Fb'][('x', w)].data,
                fo['Fb'][('y', w)].data,
                fo['Fb'][('z', w)].data,
                fo['f_sig_xx'][w].data,
                fo['f_sig_yy'][w].data,
                fo['f_sig_zz'][w].data,
                fo['f_sig_xy'][w].data,
                fo['f_sig_xz'][w].data,
                fo['f_sig_yz'][w].data,
                fo['F123_x'][w].data,
                fo['F123_y'][w].data,
                fo['F123_z'][w].data,
            ]).T.copy()

            vec_pack = self._collect_vectors_in_columns(vec_pack)

            nx = ComplexResponse.get_full_solution_vector(Nx[('x', w)])
            ny = ComplexResponse.get_full_solution_vector(Nx[('y', w)])
            nz = ComplexResponse.get_full_solution_vector(Nx[('z', w)])

            if self.rank != mpi_master():
                continue

            vec_pack = vec_pack.T.copy().reshape(-1, norb, norb)

            (Fx, Fy, Fz, f_sig_xx, f_sig_yy, f_sig_zz, f_sig_xy,
             f_sig_xz, f_sig_yz, f_x, f_y, f_z) = vec_pack


            F0 = fo['F0']

            # Get all the response matrices and Fock matrices

            kx = (self.complex_lrvec2mat(nx, nocc, norb)).T
            ky = (self.complex_lrvec2mat(ny, nocc, norb)).T
            kz = (self.complex_lrvec2mat(nz, nocc, norb)).T


            # computes all the compounded Φ_αβ, see article, where small phi
            # here is defined as:
            #   φ(κa,κb,Fb,F0) = [κa,[κb,F0]+3Fb]
            #   Φ_αα^σ = φ(κα,κα,Fα,F0) + φ(κα,κα,Fα,F0) +
            #            Σ_{ρ}^{x,y,z}[φ(κρ,κρ,Fρ,F0)] (for α=β)
            #   Φ_αβ^σ = φ(κa,κa,Fb,F0) + φ(κb,κb,Fb,F0) (for α≠β)
            # For the Φ_{αβ}^{λ+τ} component see article.

            Pxx = self.phi(kx, kx, Fx, F0)
            Pyy = self.phi(ky, ky, Fy, F0)
            Pzz = self.phi(kz, kz, Fz, F0)

            Phi_sig_xx = 2 * (3 * Pxx + Pyy + Pzz)
            Phi_sig_yy = 2 * (Pxx + 3 * Pyy + Pzz)
            Phi_sig_zz = 2 * (Pxx + Pyy + 3 * Pzz)

            Phi_sig_xy = (2 * self.phi(kx, ky, Fy, F0) +
                          2 * self.phi(ky, kx, Fx, F0))

            Phi_sig_xz = (2 * self.phi(kx, kz, Fz, F0) +
                          2 * self.phi(kz, kx, Fx, F0))

            Phi_sig_yz = (2 * self.phi(ky, kz, Fz, F0) +
                          2 * self.phi(kz, ky, Fy, F0))

            Pxx = self.phi(kx, kx, Fx, F0) + self.phi(kx, kx, Fx, F0)
            Pyy = self.phi(ky, ky, Fy, F0) + self.phi(ky, ky, Fy, F0)
            Pzz = self.phi(kz, kz, Fz, F0) + self.phi(kz, kz, Fz, F0)

            Phi_lamtau_xx = 2 * (3 * Pxx + Pyy + Pzz)
            Phi_lamtau_yy = 2 * (Pxx + 3 * Pyy + Pzz)
            Phi_lamtau_zz = 2 * (Pxx + Pyy + 3 * Pzz)

            Phi_lamtau_xy = 2 * (
                self.phi(kx, ky, Fy, F0) + self.phi(ky, kx, Fx, F0) +
                self.phi(kx, ky, Fy, F0) + self.phi(ky, kx, Fx, F0))

            Phi_lamtau_xz = 2 * (
                self.phi(kx, kz, Fz, F0) + self.phi(kz, kx, Fx, F0) +
                self.phi(kx, kz, Fz, F0) + self.phi(kz, kx, Fx, F0))

            Phi_lamtau_yz = 2 * (
                self.phi(ky, kz, Fz, F0) + self.phi(kz, ky, Fy, F0) +
                self.phi(ky, kz, Fz, F0) + self.phi(kz, ky, Fy, F0))

            # Computess all the elements of the Fock vector formed from the
            # E[4] contraction as E[4]NxNyNz = [f_{is} // f_{si}], for the
            # isotropic case, as derived in the article
            # The elements of f_{is} for each spatial component α is given by
            # an expression of the form
            # f_α = Σ_{β}^{x,y,z} [κ_{β}^{ω},Φ_{αβ}^{λ+τ}+f_{αβ}^{λ+τ}] +
            #       [κ_{β}^{-ω},Φ_{αβ}^{σ}+f_{αβ}^{σ}]

            # x
                

            # Creating the transformed total Fock matrices
            f_x += (self.commut(kx, Phi_sig_xx + Phi_lamtau_xx +  3.0 * f_sig_xx) +
                    self.commut(ky, Phi_sig_xy + Phi_lamtau_xy +  3.0 * f_sig_xy) +
                    self.commut(kz, Phi_sig_xz + Phi_lamtau_xz +  3.0 * f_sig_xz))
            
            # Taking the non redundant matrix elements {i,s} and forming the
            # anti-symmetric Fock vector
            f_x = -2. / 6 * LinearSolver.lrmat2vec(f_x.T, nocc, norb)
            f_x = self.anti_sym(f_x)
            f_iso_x[w] = f_x

            # y

            # Creating the transformed total Fock matrices
            f_y += (self.commut(kx, Phi_sig_xy + Phi_lamtau_xy + 3.0 * f_sig_xy) +
                    self.commut(ky, Phi_sig_yy + Phi_lamtau_yy + 3.0 * f_sig_yy) +
                    self.commut(kz, Phi_sig_yz + Phi_lamtau_yz + 3.0 * f_sig_yz))
            
            # Taking the non redundant matrix elements {i,s} and forming the
            # anti-symmetric Fock vector
            f_y = -2. / 6 * LinearSolver.lrmat2vec(f_y.T, nocc, norb)
            f_y = self.anti_sym(f_y)
            f_iso_y[w] = f_y

            # z

            # Creating the transformed total Fock matrices
            f_z += (self.commut(kx, Phi_sig_xz + Phi_lamtau_xz + 3.0 * f_sig_xz) +
                    self.commut(ky, Phi_sig_yz + Phi_lamtau_yz + 3.0 * f_sig_yz) +
                    self.commut(kz, Phi_sig_zz + Phi_lamtau_zz + 3.0 * f_sig_zz))
            
            # Taking the non redundant matrix elements {i,s} and forming the
            # anti-symmetric Fock vector
            f_z = -2. / 6 * LinearSolver.lrmat2vec(f_z.T, nocc, norb)
            f_z = self.anti_sym(f_z)
            f_iso_z[w] = f_z

        return {'f_iso_x': f_iso_x, 'f_iso_y': f_iso_y, 'f_iso_z': f_iso_z}

    def get_Nxy(self, w, d_a_mo, X, fock_dict, Nx, nocc, norb, molecule,
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
        :param Nx:
            A dictonary containg all the response vectors in distributed form
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

        # Get second-order gradients
        xy_dict = self.get_xy(d_a_mo, X, w, fock_dict, Nx, nocc, norb)

        # Frequencies to compute
        if self.rank == mpi_master():
            # note: use plain addition instead of sum
            freq = tuple([0.0] + [x[0] + x[1] for x in zip(w, w)])
        else:
            freq = None

        freq = self.comm.bcast(freq, root=mpi_master())

        N_total_drv = ComplexResponse(self.comm, self.ostream)
        N_total_drv.frequencies = freq

        cpp_keywords = {
            'damping', 'norm_thresh', 'lindep_thresh', 'conv_thresh',
            'max_iter', 'eri_thresh', 'timing', 'memory_profiling',
            'batch_size', 'restart', 'xcfun', 'grid_level', 'potfile',
            'electric_field', 'program_end_time', '_debug', '_block_size_factor',
            'ri_coulomb'
        }

        for key in cpp_keywords:
            setattr(N_total_drv, key, getattr(self, key))

        if self.checkpoint_file is not None:
            fpath = Path(self.checkpoint_file)
            fpath = fpath.with_name(fpath.stem)
            N_total_drv.checkpoint_file = str(fpath) + '_thg_2_full.h5'

        # commutpute second-order response vectors

        N_total_results = N_total_drv.compute(molecule, ao_basis, scf_tensors,
                                              xy_dict)

        self._is_converged = (self._is_converged and N_total_drv.is_converged)

        Nxy_dict = N_total_results['solutions']
        FXY_2_dict = N_total_results['focks']

        return (Nxy_dict, FXY_2_dict)

    def get_densities_II(self, wi, Nx, Nxy, mo, nocc, norb):
        """
        Computes the compounded densities needed for the compounded
        second-order Fock matrices used for the isotropic cubic response
        function

        :param wi:
            A list of the frequencies
        :param Nx:
            A dictonary with all the first-order response vectors in
            distributed form
        :param Nxy:
            A dict of the two index response vectors in distributed form
        :param mo:
            A matrix containing the MO coefficents
        :param nocc:
            Number of occupied orbitals
        :param norb:
            Number of orbitals

        :return:
            A list of tranformed compounded densities
        """

        distributed_density_1 = None
        distributed_density_2 = None

        for w in wi:

            nx = ComplexResponse.get_full_solution_vector(Nx[('x', w)])
            ny = ComplexResponse.get_full_solution_vector(Nx[('y', w)])
            nz = ComplexResponse.get_full_solution_vector(Nx[('z', w)])

            n_sig_xx = ComplexResponse.get_full_solution_vector(Nxy[(('N_sig_xx', w), 2 * w)])
            n_sig_yy = ComplexResponse.get_full_solution_vector(Nxy[(('N_sig_yy', w), 2 * w)])
            n_sig_zz = ComplexResponse.get_full_solution_vector(Nxy[(('N_sig_zz', w), 2 * w)])
            n_sig_xy = ComplexResponse.get_full_solution_vector(Nxy[(('N_sig_xy', w), 2 * w)])
            n_sig_xz = ComplexResponse.get_full_solution_vector(Nxy[(('N_sig_xz', w), 2 * w)])
            n_sig_yz = ComplexResponse.get_full_solution_vector(Nxy[(('N_sig_yz', w), 2 * w)])


            if self.rank == mpi_master():

                k_sig_xx = self.complex_lrvec2mat(n_sig_xx, nocc, norb)
                k_sig_yy = self.complex_lrvec2mat(n_sig_yy, nocc, norb)
                k_sig_zz = self.complex_lrvec2mat(n_sig_zz, nocc, norb)
                k_sig_xy = self.complex_lrvec2mat(n_sig_xy, nocc, norb)
                k_sig_xz = self.complex_lrvec2mat(n_sig_xz, nocc, norb)
                k_sig_yz = self.complex_lrvec2mat(n_sig_yz, nocc, norb)


                kx = self.complex_lrvec2mat(nx, nocc, norb)
                ky = self.complex_lrvec2mat(ny, nocc, norb)
                kz = self.complex_lrvec2mat(nz, nocc, norb)


                # SIGMA contributiatons #
                Dc_x = self.commut_mo_density(kx, nocc)
                Dc_y = self.commut_mo_density(ky, nocc)
                Dc_z = self.commut_mo_density(kz, nocc)

                D_sig_xx = self.commut_mo_density(k_sig_xx, nocc)
                D_sig_yy = self.commut_mo_density(k_sig_yy, nocc)
                D_sig_zz = self.commut_mo_density(k_sig_zz, nocc)
                D_sig_xy = self.commut_mo_density(k_sig_xy, nocc)
                D_sig_xz = self.commut_mo_density(k_sig_xz, nocc)
                D_sig_yz = self.commut_mo_density(k_sig_yz, nocc)


                # x #
                Dx =  3.0 * self.commut(kx, D_sig_xx)
                Dx += 3.0 * self.commut(k_sig_xx, Dc_x)
                Dx += 3.0 * self.commut(ky, D_sig_xy)
                Dx += 3.0 * self.commut(k_sig_xy, Dc_y)
                Dx += 3.0 * self.commut(kz, D_sig_xz)
                Dx += 3.0 * self.commut(k_sig_xz, Dc_z)


                # y #
                Dy =  3.0 * self.commut(kx, D_sig_xy)
                Dy += 3.0 * self.commut(k_sig_xy, Dc_x)
                Dy += 3.0 * self.commut(ky, D_sig_yy)
                Dy += 3.0 * self.commut(k_sig_yy, Dc_y)
                Dy += 3.0 * self.commut(kz, D_sig_yz)
                Dy += 3.0 * self.commut(k_sig_yz, Dc_z)


                # z #
                Dz =  3.0 * self.commut(kx, D_sig_xz)
                Dz += 3.0 * self.commut(k_sig_xz, Dc_x)
                Dz += 3.0 * self.commut(ky, D_sig_yz)
                Dz += 3.0 * self.commut(k_sig_yz, Dc_y)
                Dz += 3.0 * self.commut(kz, D_sig_zz)
                Dz += 3.0 * self.commut(k_sig_zz, Dc_z)


                # density transformation from MO to AO basis

                Dx = np.linalg.multi_dot([mo, Dx, mo.T])
                Dy = np.linalg.multi_dot([mo, Dy, mo.T])
                Dz = np.linalg.multi_dot([mo, Dz, mo.T])

                Dc_x = np.linalg.multi_dot([mo, Dc_x, mo.T])
                Dc_y = np.linalg.multi_dot([mo, Dc_y, mo.T])
                Dc_z = np.linalg.multi_dot([mo, Dc_z, mo.T])

                D_sig_xx = np.linalg.multi_dot([mo, D_sig_xx, mo.T])
                D_sig_yy = np.linalg.multi_dot([mo, D_sig_yy, mo.T])
                D_sig_zz = np.linalg.multi_dot([mo, D_sig_zz, mo.T])
                D_sig_xy = np.linalg.multi_dot([mo, D_sig_xy, mo.T])
                D_sig_xz = np.linalg.multi_dot([mo, D_sig_xz, mo.T])
                D_sig_yz = np.linalg.multi_dot([mo, D_sig_yz, mo.T])


                dist_den_1_freq = np.hstack((
                    Dc_x.real.reshape(-1, 1),
                    Dc_y.real.reshape(-1, 1),
                    Dc_z.real.reshape(-1, 1),
                    D_sig_xx.real.reshape(-1, 1),
                    D_sig_yy.real.reshape(-1, 1),
                    D_sig_zz.real.reshape(-1, 1),
                    D_sig_xy.real.reshape(-1, 1),
                    D_sig_xz.real.reshape(-1, 1),
                    D_sig_yz.real.reshape(-1, 1),
                ))

                dist_den_2_freq = np.hstack((
                    Dx.real.reshape(-1, 1),
                    Dy.real.reshape(-1, 1),
                    Dz.real.reshape(-1, 1),
                ))

            else:

                dist_den_1_freq = None
                dist_den_2_freq = None

            dist_den_1_freq = DistributedArray(dist_den_1_freq, self.comm)
            dist_den_2_freq = DistributedArray(dist_den_2_freq, self.comm)

            if distributed_density_1 is None:
                distributed_density_1 = DistributedArray(dist_den_1_freq.data,
                                                         self.comm,
                                                         distribute=False)
            else:
                distributed_density_1.append(dist_den_1_freq, axis=1)

            if distributed_density_2 is None:
                distributed_density_2 = DistributedArray(dist_den_2_freq.data,
                                                         self.comm,
                                                         distribute=False)
            else:
                distributed_density_2.append(dist_den_2_freq, axis=1)

        return distributed_density_1, distributed_density_2


    def get_xy(self, d_a_mo, X, wi, Fock, Nx, nocc, norb):
        """
        Computes the compounded gradient vectors N^{σ},N^{λ+τ} used for the
        isotropic cubic response function

        :param d_a_mo:
            The SCF density matrix in MO basis
        :param X:
            Dipole integrals
        :param wi:
            A list of the frequencies
        :param Fock:
            A dictonary containing all the Fock matricies
        :param Nx:
            A dictonary with all the first-order response vectors in
            distributed form
        :param nocc:
            The number of occupied orbitals
        :param norb:
            The number of total orbitals

        :return:
            A dictonary of compounded gradient vectors
        """

        xy_dict = {}

        one_third = 1.0 / 3.0

        for w in wi:

            vec_pack = np.array([
                Fock['Fb'][('x', w)].data,
                Fock['Fb'][('y', w)].data,
                Fock['Fb'][('z', w)].data,
                Fock['f_sig_xx'][w].data * one_third,
                Fock['f_sig_yy'][w].data * one_third,
                Fock['f_sig_zz'][w].data * one_third,
                Fock['f_sig_xy'][w].data * one_third,
                Fock['f_sig_xz'][w].data * one_third,
                Fock['f_sig_yz'][w].data * one_third
            ]).T.copy()

            vec_pack = self._collect_vectors_in_columns(vec_pack)

            nx = ComplexResponse.get_full_solution_vector(Nx[('x', w)])
            ny = ComplexResponse.get_full_solution_vector(Nx[('y', w)])
            nz = ComplexResponse.get_full_solution_vector(Nx[('z', w)])

            if self.rank != mpi_master():
                continue

            vec_pack = vec_pack.T.copy().reshape(-1, norb, norb)

            (f_x, f_y, f_z, f_sig_xx, f_sig_yy, f_sig_zz, f_sig_xy, f_sig_xz,
             f_sig_yz) = vec_pack


            mu_x = X['x']
            mu_y = X['y']
            mu_z = X['z']

            kx = (self.complex_lrvec2mat(nx, nocc, norb)).T
            ky = (self.complex_lrvec2mat(ny, nocc, norb)).T
            kz = (self.complex_lrvec2mat(nz, nocc, norb)).T


            F0 = Fock['F0']

            # BD σ gradients #

            xi_xx = self._xi(kx, kx, f_x, f_x, F0)
            xi_yy = self._xi(ky, ky, f_y, f_y, F0)
            xi_zz = self._xi(kz, kz, f_z, f_z, F0)

            x2_xx = self._x2_contract(kx.T, mu_x, d_a_mo, nocc, norb)
            x2_yy = self._x2_contract(ky.T, mu_y, d_a_mo, nocc, norb)
            x2_zz = self._x2_contract(kz.T, mu_z, d_a_mo, nocc, norb)

            key = (('N_sig_xx', w), 2 * w)
            mat = (3 * xi_xx + xi_yy + xi_zz + 0.5 * f_sig_xx).T
            xy_dict[key] = self.anti_sym(
                -LinearSolver.lrmat2vec(mat, nocc, norb))
            xy_dict[key] -= (3 * x2_xx + x2_yy + x2_zz)

            key = (('N_sig_yy', w), 2 * w)
            mat = (xi_xx + 3 * xi_yy + xi_zz + 0.5 * f_sig_yy).T
            xy_dict[key] = self.anti_sym(
                -LinearSolver.lrmat2vec(mat, nocc, norb))
            xy_dict[key] -= (x2_xx + 3 * x2_yy + x2_zz)

            key = (('N_sig_zz', w), 2 * w)
            mat = (xi_xx + xi_yy + 3 * xi_zz + 0.5 * f_sig_zz).T
            xy_dict[key] = self.anti_sym(
                -LinearSolver.lrmat2vec(mat, nocc, norb))
            xy_dict[key] -= (x2_xx + x2_yy + 3 * x2_zz)

            key = (('N_sig_xy', w), 2 * w)
            mat = (self._xi(ky, kx, f_y, f_x, F0) +
                   self._xi(kx, ky, f_x, f_y, F0) + 0.5 * f_sig_xy).T
            xy_dict[key] = self.anti_sym(
                -LinearSolver.lrmat2vec(mat, nocc, norb))
            xy_dict[key] -= self._x2_contract(ky.T, mu_x, d_a_mo, nocc, norb)
            xy_dict[key] -= self._x2_contract(kx.T, mu_y, d_a_mo, nocc, norb)

            key = (('N_sig_xz', w), 2 * w)
            mat = (self._xi(kz, kx, f_z, f_x, F0) +
                   self._xi(kx, kz, f_x, f_z, F0) + 0.5 * f_sig_xz).T
            xy_dict[key] = self.anti_sym(
                -LinearSolver.lrmat2vec(mat, nocc, norb))
            xy_dict[key] -= self._x2_contract(kz.T, mu_x, d_a_mo, nocc, norb)
            xy_dict[key] -= self._x2_contract(kx.T, mu_z, d_a_mo, nocc, norb)

            key = (('N_sig_yz', w), 2 * w)
            mat = (self._xi(kz, ky, f_z, f_y, F0) +
                   self._xi(ky, kz, f_y, f_z, F0) + 0.5 * f_sig_yz).T
            xy_dict[key] = self.anti_sym(
                -LinearSolver.lrmat2vec(mat, nocc, norb))
            xy_dict[key] -= self._x2_contract(kz.T, mu_y, d_a_mo, nocc, norb)
            xy_dict[key] -= self._x2_contract(ky.T, mu_z, d_a_mo, nocc, norb)

        return xy_dict



    def get_fock_dict_II(self, wi, density_list1, density_list2, mo, molecule,
                         ao_basis, eri_dict, dft_dict, profiler):
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
        :param eri_dict:
            The dictionary containing ERI information
        :param dft_dict:
            The dictionary containing DFT information

        :return:
            A dictonary of compounded second-order Fock-matrices
        """

        if self.rank == mpi_master():
            self._print_fock_header()

        # generate key-frequency paris

        key_freq_pairs = []
        for w in wi:
            for key in ['F123_x', 'F123_y', 'F123_z']:
                key_freq_pairs.append((key, w))

        # examine checkpoint file for distributed Focks

        if self.checkpoint_file is not None:
            fpath = Path(self.checkpoint_file)
            fpath = fpath.with_name(fpath.stem)
            fock_file = str(fpath) + '_thg_red_fock_2_full.h5'
        else:
            fock_file = None

        if self.restart:
            if self.rank == mpi_master():
                self.restart = check_distributed_focks(fock_file,
                                                       key_freq_pairs)
            self.restart = self.comm.bcast(self.restart, mpi_master())

        # read or compute distributed Focks

        if self.restart:
            dist_focks = read_distributed_focks(fock_file, self.comm,
                                                self.ostream)
        else:
            time_start_fock = time.time()

            if self._dft:
                dist_focks = self._comp_nlr_fock(mo, molecule, ao_basis,
                                                 'real', eri_dict,
                                                 dft_dict, density_list1,
                                                 density_list2, None, 'thgred_ii',
                                                 profiler)
            else:
                dist_focks = self._comp_nlr_fock(mo, molecule, ao_basis,
                                                 'real', eri_dict,
                                                 None, None, density_list2,
                                                 None, 'thgred_ii', profiler)

            self._print_fock_time(time.time() - time_start_fock)

            write_distributed_focks(fock_file, dist_focks, key_freq_pairs,
                                    self.comm, self.ostream)

        # assign distributed Focks to key-frequency pairs

        focks = {}

        for fock_index, (key, w) in enumerate(key_freq_pairs):
            if key not in focks:
                focks[key] = {}
            focks[key][w] = DistributedArray(dist_focks.data[:, fock_index],
                                             self.comm,
                                             distribute=False)

        return focks

    def get_e3(self, wi, Nx, Nxy, fo, fo2, nocc, norb):
        """
        Contracts E[3]

        :param wi:
            A list of freqs
        :param Nx:
            A dict of the single index response vectors in distributed form
        :param Nxy:
            A dict of the two index response vectors in distributed form
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

        f_iso_x = {}
        f_iso_y = {}
        f_iso_z = {}

        for w in wi:

            vec_pack = np.array([
                fo['Fb'][('x', w)].data,
                fo['Fb'][('y', w)].data,
                fo['Fb'][('z', w)].data,
                fo2[(('N_sig_xx', w), 2 * w)].data,
                fo2[(('N_sig_yy', w), 2 * w)].data,
                fo2[(('N_sig_zz', w), 2 * w)].data,
                fo2[(('N_sig_xy', w), 2 * w)].data,
                fo2[(('N_sig_xz', w), 2 * w)].data,
                fo2[(('N_sig_yz', w), 2 * w)].data,
                fo2['F123_x'][w].data,
                fo2['F123_y'][w].data,
                fo2['F123_z'][w].data,
            ]).T.copy()

            vec_pack = self._collect_vectors_in_columns(vec_pack)

            nx = ComplexResponse.get_full_solution_vector(Nx[('x', w)])
            ny = ComplexResponse.get_full_solution_vector(Nx[('y', w)])
            nz = ComplexResponse.get_full_solution_vector(Nx[('z', w)])

            n_sig_xx = ComplexResponse.get_full_solution_vector(Nxy[(('N_sig_xx', w), 2 * w)])
            n_sig_yy = ComplexResponse.get_full_solution_vector(Nxy[(('N_sig_yy', w), 2 * w)])
            n_sig_zz = ComplexResponse.get_full_solution_vector(Nxy[(('N_sig_zz', w), 2 * w)])
            n_sig_xy = ComplexResponse.get_full_solution_vector(Nxy[(('N_sig_xy', w), 2 * w)])
            n_sig_xz = ComplexResponse.get_full_solution_vector(Nxy[(('N_sig_xz', w), 2 * w)])
            n_sig_yz = ComplexResponse.get_full_solution_vector(Nxy[(('N_sig_yz', w), 2 * w)])


            if self.rank != mpi_master():
                continue

            vec_pack = vec_pack.T.copy().reshape(-1, norb, norb)

            (f_x, f_y, f_z, f_sig_xx, f_sig_yy, f_sig_zz, f_sig_xy, f_sig_xz, f_sig_yz, F123_x, F123_y, F123_z) = vec_pack

            f_sig_xx = f_sig_xx.T.conj()
            f_sig_yy = f_sig_yy.T.conj()
            f_sig_zz = f_sig_zz.T.conj()
            f_sig_xy = f_sig_xy.T.conj()
            f_sig_xz = f_sig_xz.T.conj()
            f_sig_yz = f_sig_yz.T.conj()


            F0_a = fo['F0']

            # Response

            k_x = (self.complex_lrvec2mat(nx, nocc, norb)).T
            k_y = (self.complex_lrvec2mat(ny, nocc, norb)).T
            k_z = (self.complex_lrvec2mat(nz, nocc, norb)).T

            k_sig_xx = (self.complex_lrvec2mat(n_sig_xx, nocc, norb)).T
            k_sig_yy = (self.complex_lrvec2mat(n_sig_yy, nocc, norb)).T
            k_sig_zz = (self.complex_lrvec2mat(n_sig_zz, nocc, norb)).T
            k_sig_xy = (self.complex_lrvec2mat(n_sig_xy, nocc, norb)).T
            k_sig_xz = (self.complex_lrvec2mat(n_sig_xz, nocc, norb)).T
            k_sig_yz = (self.complex_lrvec2mat(n_sig_yz, nocc, norb)).T

            # Focks #

            # x

            zeta_sig_xx = self._xi(k_x, k_sig_xx, f_x, f_sig_xx, F0_a) + self._xi(k_x, 2.0 * k_sig_xx, f_x, 2.0 * f_sig_xx, F0_a) 
            zeta_sig_yy = self._xi(k_x, k_sig_yy, f_x, f_sig_yy, F0_a) + self._xi(k_x, 2.0 * k_sig_yy, f_x, 2.0 * f_sig_yy, F0_a)
            zeta_sig_zz = self._xi(k_x, k_sig_zz, f_x, f_sig_zz, F0_a) + self._xi(k_x, 2.0 * k_sig_zz, f_x, 2.0 * f_sig_zz, F0_a)
            zeta_sig_xy = self._xi(k_y, k_sig_xy, f_y, f_sig_xy, F0_a) + self._xi(k_y, 2.0 * k_sig_xy, f_y, 2.0 * f_sig_xy, F0_a)
            zeta_sig_xz = self._xi(k_z, k_sig_xz, f_z, f_sig_xz, F0_a) + self._xi(k_z, 2.0 * k_sig_xz, f_z, 2.0 * f_sig_xz, F0_a)


            X_terms = (zeta_sig_xx + zeta_sig_xy + zeta_sig_xz).T + (0.5 * F123_x).T
            Ff_x = -2 * LinearSolver.lrmat2vec(X_terms, nocc, norb)
            Ff_x = self.anti_sym(Ff_x)
            f_iso_x[w] = Ff_x

            # y

            zeta_sig_yx = self._xi(k_x, k_sig_xy, f_x, f_sig_xy, F0_a) + self._xi(k_x, 2.0 * k_sig_xy, f_x, 2.0 * f_sig_xy, F0_a)
            zeta_sig_yy = self._xi(k_y, k_sig_yy, f_y, f_sig_yy, F0_a) + self._xi(k_y, 2.0 * k_sig_yy, f_y, 2.0 * f_sig_yy, F0_a)
            zeta_sig_yz = self._xi(k_z, k_sig_yz, f_z, f_sig_yz, F0_a) + self._xi(k_z, 2.0 * k_sig_yz, f_z, 2.0 * f_sig_yz, F0_a)


            Y_terms = (zeta_sig_yx + zeta_sig_yy + zeta_sig_yz).T + (0.5 * F123_y).T
            Ff_y = -2 * LinearSolver.lrmat2vec(Y_terms, nocc, norb)
            Ff_y = self.anti_sym(Ff_y)
            f_iso_y[w] = Ff_y

            # z
            zeta_sig_zx = self._xi(k_x, k_sig_xz, f_x, f_sig_xz, F0_a) + self._xi(k_x, 2.0 * k_sig_xz, f_x, 2.0 * f_sig_xz, F0_a) 
            zeta_sig_zy = self._xi(k_y, k_sig_yz, f_y, f_sig_yz, F0_a) + self._xi(k_y, 2.0 * k_sig_yz, f_y, 2.0 * f_sig_yz, F0_a)
            zeta_sig_zz = self._xi(k_z, k_sig_zz, f_z, f_sig_zz, F0_a) + self._xi(k_z, 2.0 * k_sig_zz, f_z, 2.0 * f_sig_zz, F0_a)


            Z_terms = (zeta_sig_zx + zeta_sig_zy + zeta_sig_zz).T + (0.5 * F123_z).T
            Ff_z = -2 * LinearSolver.lrmat2vec(Z_terms, nocc, norb)
            Ff_z = self.anti_sym(Ff_z)
            f_iso_z[w] = Ff_z

        return {'f_iso_x': f_iso_x, 'f_iso_y': f_iso_y, 'f_iso_z': f_iso_z}

    def get_other_terms(self, wi, track, X, Nx, Nxy, da, nocc, norb):
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
        :param Nx:
            A dictonary with all the respone vectors in distributed form
        :param Nxy:
            A dictonary containing all the two-index response vectors in
            distributed form
        :param da:
            The SCF density matrix in MO basis
        :param nocc:
            The number of occupied orbitals
        :param norb:
            The total number of orbitals

        :return:
            A dictonary of final X[2],A[2] contraction values
        """

        na_a3_nx_ny_dict = {}
        na_x3_ny_nz_dict = {}
        na_x2_nyz_dict = {}
        nx_a2_nyz_dict = {}

        comp_per_freq = len(track) // len(wi)

        inp_list = []

        for j in range(len(wi)):
            for i in range(j * comp_per_freq, (j + 1) * comp_per_freq):
                comp_i = track[i]

                vals = comp_i.split(',')
                op_a, op_b, op_c, op_d = vals[0]
                w = float(vals[1])
                wa = 3 * float(vals[1])
                wb = float(vals[1])
                wc = float(vals[2])
                wd = float(vals[3])

                inp_list.append({
                    'freq': w,
                    'Na': Nx[(op_a, wa)],
                    'Nb': Nx[(op_b, wb)],
                    'Nc': Nx[(op_c, wc)],
                    'Nd': Nx[(op_d, wd)],
                    'A': X[op_a],
                    'B': X[op_b],
                    'C': X[op_c],
                    'D': X[op_d],
                })

        list_x3_a3 = [self.get_x3_a3(inp, da, nocc, norb) for inp in inp_list]

        if self.rank == mpi_master():
            for terms in list_x3_a3:
                if terms['key'] not in na_x3_ny_nz_dict:
                    na_x3_ny_nz_dict[terms['key']] = 0.0
                if terms['key'] not in na_a3_nx_ny_dict:
                    na_a3_nx_ny_dict[terms['key']] = 0.0
                na_x3_ny_nz_dict[terms['key']] += terms['x3']
                na_a3_nx_ny_dict[terms['key']] += terms['a3']

        inp_list = []

        for i in range(len(wi)):
            vals = track[i * comp_per_freq].split(',')
            w = float(vals[1])
            wa = 3 * float(vals[1])
            wb = float(vals[1])
            wc = float(vals[2])
            wd = float(vals[3])

            wcd = wb + wd
            wbd = wb + wd

            for op_a in 'xyz':
                Na = Nx[(op_a, wa)]
                A = X[op_a]

                for op_b in 'xyz':
                    op_ab = op_a + op_b if op_a <= op_b else op_b + op_a

                    # CD
                    Ncd = Nxy[(('N_sig_' + op_ab, w), wcd)]
                    Nb = Nx[(op_b, w)]
                    B = X[op_b]

                    inp_list.append({
                        'flag': 'CD',
                        'freq': w,
                        'Ncd': Ncd,
                        'Na': Na,
                        'Nb': Nb,
                        'A': A,
                        'B': B,
                    })

                    # BD
                    op_c = op_b
                    op_ac = op_ab
                    Nbd = Nxy[(('N_sig_' + op_ac, w), wbd)]
                    Nc = Nx[(op_c, wc)]
                    C = X[op_c]

                    inp_list.append({
                        'flag': 'BD',
                        'freq': w,
                        'Nbd': Nbd,
                        'Na': Na,
                        'Nc': Nc,
                        'A': A,
                        'C': C,
                    })

        list_x2_a2 = [self.get_x2_a2(inp, da, nocc, norb) for inp in inp_list]

        if self.rank == mpi_master():
            for terms in list_x2_a2:
                if terms['key'] not in na_x2_nyz_dict:
                    na_x2_nyz_dict[terms['key']] = 0.0
                if terms['key'] not in nx_a2_nyz_dict:
                    nx_a2_nyz_dict[terms['key']] = 0.0
                na_x2_nyz_dict[terms['key']] += terms['x2']
                nx_a2_nyz_dict[terms['key']] += terms['a2']

            return {
                'NaX3NyNz': na_x3_ny_nz_dict,
                'NaA3NxNy': na_a3_nx_ny_dict,
                'NaX2Nyz': na_x2_nyz_dict,
                'NxA2Nyz': nx_a2_nyz_dict,
            }


        return None



    def get_s4_and_r4(self, wi, Nx, track, D0, nocc, norb):
        """
        Computes the S4 contractions

        :param wi:
            A list of all the freqs
        :param Nx:
            A dict with all the response vectors in distributed form
        :param track:
            A list containing information about all the components that are to
            be computed
        :param D0:
            The SCF density in MO basis
        :param nocc:
            The number of occupied obritals
        :param norb:
            The number of total orbitals

        :return:
            A dictonary of final S[4] contraction values
        """

        S4 = {}
        R4 = {}

        comp_per_freq = len(track) // len(wi)

        inp_list = []

        for j in range(len(wi)):
            vals = track[j * comp_per_freq].split(',')
            
            w = float(vals[1])
            w1 = float(vals[1])
            w2 = float(vals[2])
            w3 = float(vals[3])
            w_s = w1 + w2 + w3

            for i in range(j * comp_per_freq, (j + 1) * comp_per_freq):
                comp_i = track[i]
                op = comp_i[0]

                inp_dict = {
                    'w': w,
                    'w1': w1,
                    'w2': w2,
                    'w3': w3,
                    'op': op,
                    'Nb': Nx[(comp_i[1], w1)],
                    'Nc': Nx[(comp_i[2], w2)],
                    'Nd': Nx[(comp_i[3], w3)],
                }

                if self.damping > 0:
                    inp_dict.update({
                        'Nb': Nx[(comp_i[1], w1)],
                        'Nc': Nx[(comp_i[2], w2)],
                        'Nd': Nx[(comp_i[3], w3)],
                        'Na': Nx[(comp_i[0], w_s)],
                    })

                inp_list.append(inp_dict)

        list_s4_r4 = [
            self.get_s4_and_r4_terms(inp, D0, nocc, norb) for inp in inp_list
        ]

        if self.rank == mpi_master():
            local_s4_dict = {}
            local_r4_dict = {}

            for terms in list_s4_r4:
                if terms['s4_key'] not in local_s4_dict:
                    local_s4_dict[terms['s4_key']] = 0.0
                local_s4_dict[terms['s4_key']] += terms['s4']

                if terms['r4_key'] not in local_r4_dict:
                    local_r4_dict[terms['r4_key']] = 0.0
                local_r4_dict[terms['r4_key']] += terms['r4']

            list_s4_r4 = []

            for s4_key, r4_key in zip(list(local_s4_dict.keys()),
                                      list(local_r4_dict.keys())):
                list_s4_r4.append({
                    's4_key': s4_key,
                    's4': local_s4_dict[s4_key],
                    'r4_key': r4_key,
                    'r4': local_r4_dict[r4_key],
                })

            for terms in list_s4_r4:
                s4_key = terms['s4_key']
                r4_key = terms['r4_key']
                if s4_key not in S4:
                    S4[s4_key] = 0.0
                if r4_key not in R4:
                    R4[r4_key] = 0.0
                S4[s4_key] += terms['s4']
                R4[r4_key] += terms['r4']

            return S4, R4
        else:
            return None, None



    def get_s4_and_r4_terms(self, inp_dict, D0, nocc, norb):
        """
        Computes S[4] and R[4] contributions.

        :param inp_dict:
            A dictionary containing input data for computing S[4] and R[4].
        :param D0:
            The SCF density matrix in MO basis
        :param nocc:
            The number of occupied orbitals
        :param norb:
            The total number of orbitals

        :return:
            Dictionaries containing S[4] and R[4].
        """

        s4_term = 0.0
        r4_term = 0.0

        w = inp_dict['w']
        w1 = inp_dict['w1']
        w2 = inp_dict['w2']
        w3 = inp_dict['w3']
        op = inp_dict['op']

        s4_key = (op, w)
        r4_key = (op, w1)

        Nb = ComplexResponse.get_full_solution_vector(inp_dict['Nb'])
        Nc = ComplexResponse.get_full_solution_vector(inp_dict['Nc'])
        Nd = ComplexResponse.get_full_solution_vector(inp_dict['Nd'])

        if self.damping > 0:
            Na = ComplexResponse.get_full_solution_vector(inp_dict['Na'])

        if self.rank == mpi_master():
            kB = self.complex_lrvec2mat(Nb, nocc, norb)
            kC = self.complex_lrvec2mat(Nc, nocc, norb)
            kD = self.complex_lrvec2mat(Nd, nocc, norb)

            s4_term -= w1 * self._s4(kB, kC, kD, D0, nocc, norb)
            s4_term -= w2 * self._s4(kC, kB, kD, D0, nocc, norb)
            s4_term -= w3 * self._s4(kD, kB, kC, D0, nocc, norb)

            if self.damping > 0:
                kA = self.complex_lrvec2mat(Na, nocc, norb)

                Nb_h = self.flip_xy(Nb)
                Nc_h = self.flip_xy(Nc)
                Nd_h = self.flip_xy(Nd)

                r4_term += 1j * self.damping * np.dot(Nd_h, self._s4_for_r4(kA.T, kB, kC, D0, nocc, norb))
                r4_term += 1j * self.damping * np.dot(Nc_h, self._s4_for_r4(kA.T, kB, kD, D0, nocc, norb))
                r4_term += 1j * self.damping * np.dot(Nd_h, self._s4_for_r4(kA.T, kC, kB, D0, nocc, norb))
                r4_term += 1j * self.damping * np.dot(Nb_h, self._s4_for_r4(kA.T, kC, kD, D0, nocc, norb))
                r4_term += 1j * self.damping * np.dot(Nc_h, self._s4_for_r4(kA.T, kD, kB, D0, nocc, norb))
                r4_term += 1j * self.damping * np.dot(Nb_h, self._s4_for_r4(kA.T, kD, kC, D0, nocc, norb))

            return {
                's4_key': s4_key,
                'r4_key': r4_key,
                's4': s4_term,
                'r4': r4_term,
            }
        else:
            return None

    def get_t4(self, wi, e4_dict, Nx, track, da, nocc, norb):
        """
        Computes the contraction of the E[4] tensor with that of the S[4] and
        R[4] tensors to return the contraction of T[4] as a dictonary of
        vectors. T[4]NxNyNz = (E^[4]-ω_1S^[4]-ω_1S^[4]-ω_3S^[4]-γiR^[4])

        :param wi:
            A list of all the freqs
        :param e4_dict:
            A dictonary of all the E[4] contraction
        :param Nx:
            A dictonray containng all the response vectors in distributed form
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

        T4 = {}
        S4, R4 = self.get_s4_and_r4(wi, Nx, track, da, nocc, norb)

        comp_per_freq = len(track) // len(wi)

        for i in range(len(wi)):
            vals = track[i * comp_per_freq].split(',')
            w = float(vals[1])
            ww = float(vals[1])

            na_x = ComplexResponse.get_full_solution_vector(Nx[('x', 3 * w)])
            na_y = ComplexResponse.get_full_solution_vector(Nx[('y', 3 * w)])
            na_z = ComplexResponse.get_full_solution_vector(Nx[('z', 3 * w)])

            if self.rank == mpi_master():
                t4val = (np.dot(na_x, e4_dict['f_iso_x'][ww] - S4[('x', ww)]) +
                         np.dot(na_y, e4_dict['f_iso_y'][ww] - S4[('y', ww)]) +
                         np.dot(na_z, e4_dict['f_iso_z'][ww] - S4[('z', ww)]))

                if self.damping > 0:
                    t4val += (R4[('x', ww)] + R4[('y', ww)] + R4[('z', ww)])

                T4[(ww, ww, ww)] = -(1. / 15) * t4val

        if self.rank == mpi_master():
            return T4
        else:
            return None

    def get_t3(self, freqs, e3_dict, Nx, track, nocc, norb):
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
        :param Nx:
            A dictonray containng all the response vectors in distributed form
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

            na_x = ComplexResponse.get_full_solution_vector(Nx[('x', 3.0 * w)])
            na_y = ComplexResponse.get_full_solution_vector(Nx[('y', 3.0 * w)])
            na_z = ComplexResponse.get_full_solution_vector(Nx[('z', 3.0 * w)])

            if self.rank == mpi_master():

                t3val = (np.dot(na_x, e3_dict['f_iso_x'][w]) +
                         np.dot(na_y, e3_dict['f_iso_y'][w]) +
                         np.dot(na_z, e3_dict['f_iso_z'][w]))

                t3_term[(w, w, w)] = 1. / 15 * t3val

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

        na_x3_ny_nz = 0.0
        na_a3_nx_ny = 0.0

        w = inp_dict['freq']

        Na = ComplexResponse.get_full_solution_vector(inp_dict['Na'])
        Nb = ComplexResponse.get_full_solution_vector(inp_dict['Nb'])
        Nc = ComplexResponse.get_full_solution_vector(inp_dict['Nc'])
        Nd = ComplexResponse.get_full_solution_vector(inp_dict['Nd'])

        if self.rank == mpi_master():
            kb = self.complex_lrvec2mat(Nb, nocc, norb)
            kc = self.complex_lrvec2mat(Nc, nocc, norb)
            kd = self.complex_lrvec2mat(Nd, nocc, norb)

            A = inp_dict['A']
            B = inp_dict['B']
            C = inp_dict['C']
            D = inp_dict['D']

            # Na X[3]NyNz

            na_x3_ny_nz -= np.dot(Na.T, self._x3_contract(kc, kd, B, da, nocc, norb))
            na_x3_ny_nz -= np.dot(Na.T, self._x3_contract(kd, kc, B, da, nocc, norb))
            na_x3_ny_nz -= np.dot(Na.T, self._x3_contract(kd, kb, C, da, nocc, norb))
            na_x3_ny_nz -= np.dot(Na.T, self._x3_contract(kb, kd, C, da, nocc, norb))
            na_x3_ny_nz -= np.dot(Na.T, self._x3_contract(kb, kc, D, da, nocc, norb))
            na_x3_ny_nz -= np.dot(Na.T, self._x3_contract(kc, kb, D, da, nocc, norb))

            # NaA[3]NxNy

            na_a3_nx_ny += np.dot(self._a3_contract(kb, kc, A, da, nocc, norb),Nd)
            na_a3_nx_ny += np.dot(self._a3_contract(kb, kd, A, da, nocc, norb),Nc)
            na_a3_nx_ny += np.dot(self._a3_contract(kc, kb, A, da, nocc, norb),Nd)
            na_a3_nx_ny += np.dot(self._a3_contract(kc, kd, A, da, nocc, norb),Nb)
            na_a3_nx_ny += np.dot(self._a3_contract(kd, kb, A, da, nocc, norb),Nc)
            na_a3_nx_ny += np.dot(self._a3_contract(kd, kc, A, da, nocc, norb),Nb)

            return {
                'key': (w, w, w),
                'x3': (1. / 15) * na_x3_ny_nz,
                'a3': (1. / 15) * na_a3_nx_ny,
            }
        else:
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
        A = inp_dict['A']

        Na = ComplexResponse.get_full_solution_vector(inp_dict['Na'])

        if inp_dict['flag'] == 'CD':
            Ncd = 2.0  * ComplexResponse.get_full_solution_vector(inp_dict['Ncd'])
            Nb = ComplexResponse.get_full_solution_vector(inp_dict['Nb'])

        elif inp_dict['flag'] == 'BD':
            Nbd = ComplexResponse.get_full_solution_vector(inp_dict['Nbd'])
            Nc = ComplexResponse.get_full_solution_vector(inp_dict['Nc'])

        if self.rank == mpi_master():

            if inp_dict['flag'] == 'CD':
                kcd = self.complex_lrvec2mat(Ncd, nocc, norb)
                kb = self.complex_lrvec2mat(Nb, nocc, norb)
                B = inp_dict['B']

                na_x2_nyz += np.dot(Na.T,self._x2_contract(kcd, B, da, nocc, norb))
                nx_a2_nyz += np.dot(self._a2_contract(kb, A, da, nocc, norb), Ncd)
                nx_a2_nyz += np.dot(self._a2_contract(kcd, A, da, nocc, norb),Nb)

            elif inp_dict['flag'] == 'BD':
                kbd = self.complex_lrvec2mat(Nbd, nocc, norb)
                kc = self.complex_lrvec2mat(Nc, nocc, norb)
                C = inp_dict['C']

                na_x2_nyz += np.dot(Na.T, self._x2_contract(kbd, C, da, nocc, norb))
                nx_a2_nyz += np.dot(self._a2_contract(kc, A, da, nocc, norb), Nbd)
                nx_a2_nyz += np.dot(self._a2_contract(kbd, A, da, nocc, norb), Nc)

            return {
                'key': (w, w, w),
                'x2': -(1. / 15) * na_x2_nyz,
                'a2': -(1. / 15) * nx_a2_nyz,
            }

        else:
            return {}


    def _print_results(self, rsp_results):
        """
        Prints the results from the TPA calculation.

        :param rsp_results:
            A dictionary containing the results of response calculation.
        """

        freqs = rsp_results['frequencies']
        component = 'THG'

        self.ostream.print_blank()

        w_str = 'Isotropic Average of Re(gamma(-3w; w, w, w)) for Thid-Harmonic Generation (THG) intensities'
        self.ostream.print_header(w_str)
        self.ostream.print_header('=' * (len(w_str) + 2))
        self.ostream.print_blank()

        # Updated header with only two columns
        title = '{:>12s} {:>25s}'.format('Frequency', 'THG intensity (a.u.)')
        width = len(title)
        self.ostream.print_header(title.ljust(width))
        self.ostream.print_header(('-' * len(title)).ljust(width))

        tensor = rsp_results.get(component)
        if tensor is not None:
            for w in freqs:
                key = (w, w, w)
                val = tensor.get(key)
                if val is not None:
                    intensity = val.real if hasattr(val, 'real') else val
                    line = '{:>12.4f} {:>25.8f}'.format(w, intensity)
                    self.ostream.print_header(line.ljust(width))
                else:
                    self.ostream.print_header(f"{w:>12.4f} {'N/A':>25}")

        self.ostream.print_blank()
        self.ostream.print_header('Reference:'.ljust(width))
        self.ostream.print_blank()
        self.ostream.flush()



    def _print_results2(self, rsp_results):
        """
        Prints the results from the TPA calculation.

        :param rsp_results:
            A dictionary containing the results of response calculation.
        """

        freqs = rsp_results['frequencies']
        components = ['THG', 't4_dict', 't3_dict', 'NaX3NyNz', 'NaA3NxNy', 'NaX2Nyz', 'NxA2Nyz']   # list of components to print

        self.ostream.print_blank()

        # Title section
        w_str = 'Isotropic Average of Re(gamma(-3w; w, w, w)) for Third-Harmonic Generation (THG) intensities'
        self.ostream.print_header(w_str)
        self.ostream.print_header('=' * (len(w_str) + 2))
        self.ostream.print_blank()

        for component in components:
            tensor = rsp_results.get(component)
            if tensor is None:
                continue

            # Header for this component
            component_title = f'Component: {component}'
            self.ostream.print_header(component_title)
            self.ostream.print_header('-' * len(component_title))
            self.ostream.print_blank()

            # Column headers
            title = '{:>12s} {:>25s}'.format('Frequency', f'{component} value (a.u.)')
            width = len(title)
            self.ostream.print_header(title.ljust(width))
            self.ostream.print_header(('-' * len(title)).ljust(width))

            # Values
            for w in freqs:
                key = (w, w, w)
                val = tensor.get(key)
                if val is not None:
                    value = val.real if hasattr(val, 'real') else val
                    line = '{:>12.4f} {:>25.8f}'.format(w, value)
                    self.ostream.print_header(line.ljust(width))
                else:
                    self.ostream.print_header(f"{w:>12.4f} {'N/A':>25}")

            self.ostream.print_blank()

        # Reference footer
        self.ostream.print_header('Reference:'.ljust(width))
        self.ostream.print_blank()
        self.ostream.flush()

    def get_comp(self, freqs):
        """
        Makes a list of all the gamma tensor components that are to be computed
        for printing purposes and for the contraction of X[3],X[2],A[3],A[2]

        :param freqs:
            A list of all the frequencies for the thg calculation

        :return:
            A list of gamma tensors components inlcuded in the isotropic cubic
            response with their corresponding frequencies
        """

        comp_iso = []
        spat_A = 'xyz'

        for w in freqs:
            w_key = '{},{},{}'.format(w, w, w)
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
        Prints thg component.

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

    @staticmethod
    def get_spectrum(rsp_results, x_unit):
        """
        Gets Third-harmonic gerneration spectrum.

        :param rsp_results:
            A dictonary containing the results of response calculation.
        :param x_unit:
            The unit of x-axis.

        :return:
            A dictionary containing photon energies and thg cross-sections.
        """

        assert_msg_critical(
            x_unit.lower() in ['au', 'ev', 'nm'],
            'TpaDriver.get_spectrum: x_unit should be au, ev or nm')

        au2ev = hartree_in_ev()
        auxnm = 1.0 / hartree_in_inverse_nm()

        # conversion factor for thg cross-sections in GM
        # * a0 in cm
        # * c in cm/s
        # * broadening parameter not included in au2gm
        alpha = fine_structure_constant()
        a0_in_cm = bohr_in_angstrom() * 1.0e-8
        c_in_cm_per_s = speed_of_light_in_vacuum_in_SI() * 100.0
        au2gm = (8.0 * np.pi**2 * alpha * a0_in_cm**5) / c_in_cm_per_s * 1.0e+50

        gamma = rsp_results['gamma']

        spectrum = {'x_data': [], 'y_data': []}

        if x_unit.lower() == 'au':
            spectrum['x_label'] = 'Photon energy [a.u.]'
        elif x_unit.lower() == 'ev':
            spectrum['x_label'] = 'Photon energy [eV]'
        elif x_unit.lower() == 'nm':
            spectrum['x_label'] = 'Wavelength [nm]'

        spectrum['y_label'] = 'thg cross-section [GM]'

        freqs = rsp_results['frequencies']

        for w in freqs:
            if w == 0.0:
                continue

            if x_unit.lower() == 'au':
                spectrum['x_data'].append(w)
            elif x_unit.lower() == 'ev':
                spectrum['x_data'].append(au2ev * w)
            elif x_unit.lower() == 'nm':
                spectrum['x_data'].append(auxnm / w)

            cross_section_in_GM = gamma[(w, -w, w)].imag * w**2 * au2gm

            spectrum['y_data'].append(cross_section_in_GM)

        return spectrum

    def _print_spectrum(self, spectrum, width):
        """
        Prints Third-harmonic gerneration spectrum.

        :param spectrum:
            The spectrum.
        :param width:
            The width of the output.
        """

        self.ostream.print_blank()

        title = 'Third-harmonic gerneration Spectrum'
        self.ostream.print_header(title)
        self.ostream.print_header('=' * (len(title) + 2))
        self.ostream.print_blank()

        if len(self.frequencies) == 1 and self.frequencies[0] == 0.0:
            text = '*** No Third-harmonic gerneration spectrum at zero frequency.'
            self.ostream.print_header(text.ljust(width))
            self.ostream.print_blank()
            return

        assert_msg_critical(
            '[a.u.]' in spectrum['x_label'],
            'TpaDriver._print_spectrum: In valid unit in x_label')
        assert_msg_critical(
            '[GM]' in spectrum['y_label'],
            'TpaDriver._print_spectrum: In valid unit in y_label')

        title = '{:<20s}{:<20s}{:>15s}'.format('Frequency[a.u.]',
                                               'Frequency[eV]',
                                               'thg cross-section[GM]')
        self.ostream.print_header(title.ljust(width))
        self.ostream.print_header(('-' * len(title)).ljust(width))

        for w, cross_section in zip(spectrum['x_data'], spectrum['y_data']):
            output = '{:<20.4f}{:<20.5f}{:>13.8f}'.format(
                w, w * hartree_in_ev(), cross_section)
            self.ostream.print_header(output.ljust(width))

        self.ostream.print_blank()
        self.ostream.flush()
