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

from mpi4py import MPI
import numpy as np
import time
import math

from .veloxchemlib import (OverlapGeom100Driver, KineticEnergyGeom100Driver,
                           NuclearPotentialGeom100Driver,
                           NuclearPotentialGeom010Driver, FockGeom1000Driver)
from .veloxchemlib import RIFockGradDriver
from .veloxchemlib import T4CScreener
from .veloxchemlib import XCMolecularGradient
from .veloxchemlib import mpi_master, mat_t
from .veloxchemlib import make_matrix
from .matrices import Matrices
from .molecularbasis import MolecularBasis
from .tdaeigensolver import TdaEigenSolver
from .tddftorbitalresponse import TddftOrbitalResponse
from .gradientdriver import GradientDriver
from .scfgradientdriver import ScfGradientDriver
from .firstorderprop import FirstOrderProperties
from .errorhandler import assert_msg_critical
from .inputparser import parse_input
from .sanitychecks import (molecule_sanity_check, scf_results_sanity_check,
                           dft_sanity_check)


class TddftGradientDriver(GradientDriver):
    """
    Implements the analytic gradient driver for excited states at the
    Tamm-Dancoff approximation (TDA) and random phase approximation (RPA)
    level based on an SCF ground state (both HF and DFT).

    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.

    Instance variables
        - gradient: The gradient.
        - tamm_dancoff: Flag if Tamm-Dancoff approximation is employed.
        - state_deriv_index: The excited state of interest.
        - do_first_order_prop: Controls the printout of first-order properties.
        - relaxed_dipole_moment: The relaxed excited-state dipole moment.
        - delta_h: The displacement for finite difference.
        - do_four_point: Flag for four-point finite difference.
    """

    # TODO: add response driver here? Save scf_drv as an instance variable?
    def __init__(self, scf_drv):
        """
        Initializes gradient driver.

        :param scf_drv:
            The Scf driver.
        """

        super().__init__(scf_drv.comm, scf_drv.ostream)
        self._debug = scf_drv._debug
        self._scf_drv = scf_drv

        self._rsp_results = None

        self._block_size_factor = 4

        self._xcfun_ldstaging = scf_drv._xcfun_ldstaging

        self.timing = scf_drv.timing

        self.flag = 'RPA Gradient Driver'

        # flag on whether RPA or TDA is calculated
        self.tamm_dancoff = False

        # excited state information; if it is set to None,
        # all available states will be calculated.
        self.state_deriv_index = None

        self.do_first_order_prop = False
        self.relaxed_dipole_moment = None
        self.unrelaxed_dipole_moment = None

        # option dictionaries from input
        # TODO: cleanup
        self.method_dict = {}
        self.orbrsp_dict = {}
        self.grad_dict = {}

        self._input_keywords['gradient'].update({
            'tamm_dancoff': ('bool', 'whether RPA or TDA is calculated'),
            'state_deriv_index': ('seq_fixed_int', 'excited state information'),
            'do_first_order_prop': ('bool', 'do first-order property'),
        })

    def update_settings(self,
                        grad_dict,
                        rsp_dict,
                        orbrsp_dict=None,
                        method_dict=None):
        """
        Updates settings in gradient driver.

        :param grad_dict:
            The input dictionary of gradient settings group.
        :param rsp_dict:
            The input dictionary of response settings  group.
        :param orbrsp_dict:
            The input dictionary of orbital response settings group.
        :param method_dict:
            The input dicitonary of method settings group.
        """

        if method_dict is None:
            method_dict = {}

        if orbrsp_dict is None:
            orbrsp_dict = {}

        super().update_settings(grad_dict, method_dict)

        grad_keywords = {
            key: val[0] for key, val in self._input_keywords['gradient'].items()
        }

        if 'tamm_dancoff' in rsp_dict.keys():
            grad_dict['tamm_dancoff'] = rsp_dict['tamm_dancoff']
            orbrsp_dict['tamm_dancoff'] = rsp_dict['tamm_dancoff']

        if 'do_first_order_prop' in grad_dict.keys():
            orbrsp_dict['do_first_order_prop'] = grad_dict[
                'do_first_order_prop']

        parse_input(self, grad_keywords, grad_dict)

        if self.tamm_dancoff:
            self.flag = 'TDA Gradient Driver'

        if self.state_deriv_index is not None:
            orbrsp_dict['state_deriv_index'] = self.state_deriv_index

        self.grad_dict = dict(grad_dict)
        self.rsp_dict = dict(rsp_dict)
        self.method_dict = dict(method_dict)
        self.orbrsp_dict = dict(orbrsp_dict)

    def compute(self, molecule, basis, scf_drv, rsp_drv, rsp_results):
        """
        Performs calculation of analytical or numerical gradient.

        :param molecule:
            The molecule.
        :param basis:
            The AO basis set.
        :param rsp_drv:
            The RPA or TDA driver.
        :param rsp_results:
            The results from the RPA or TDA calculation.
        """

        # sanity checks
        molecule_sanity_check(molecule)
        scf_results_sanity_check(self, self._scf_drv.scf_tensors)
        dft_sanity_check(self, 'compute')

        # TODO: replace with a sanity check?
        if isinstance(rsp_drv, TdaEigenSolver):
            self.tamm_dancoff = True

        if self.rank == mpi_master():
            if self.tamm_dancoff:
                self.flag = 'TDA Gradient Driver'

            all_states = list(np.arange(1, len(rsp_results['eigenvalues']) + 1))
            if self.state_deriv_index is not None:

                # make self.state_deriv_index a tuple in case it is set as an integer
                # TODO: reconsider how state_deriv_index is used
                if isinstance(self.state_deriv_index, int):
                    self.state_deriv_index = (self.state_deriv_index,)

                error_message = 'TdscfGradientDriver: some of the '
                error_message += 'selected states have not been calculated.'
                assert_msg_critical(
                    set(self.state_deriv_index).issubset(all_states),
                    error_message)
            else:
                self.state_deriv_index = all_states

            if self.numerical:
                self.print_header(self.state_deriv_index[:1])
            else:
                self.print_header(self.state_deriv_index)

        self.state_deriv_index = self.comm.bcast(self.state_deriv_index,
                                                 root=mpi_master())

        start_time = time.time()

        # NOTE: the numerical gradient is calculated for the first state only.
        if self.numerical:
            assert_msg_critical(not self.unrelaxed,
                            'TddftGradientDriver: Numerical unrelaxed gradient not available')
            scf_drv.ostream.mute()
            rsp_drv.ostream.mute()
            self.compute_numerical(molecule, basis, scf_drv, rsp_drv,
                                   rsp_results)
            rsp_drv.ostream.unmute()
            scf_drv.ostream.unmute()
        else:
            self.compute_analytical(molecule, basis, rsp_results)

        if self.rank == mpi_master():
            self.print_geometry(molecule)
            if self.numerical:
                self.print_gradient(molecule, self.state_deriv_index[:1])
            else:
                self.print_gradient(molecule, self.state_deriv_index)

            valstr = '*** Time spent in gradient calculation: '
            valstr += '{:.2f} sec ***'.format(time.time() - start_time)
            self.ostream.print_header(valstr)
            self.ostream.print_blank()
            self.ostream.flush()

    def compute_analytical(self, molecule, basis, rsp_results):
        """
        Performs calculation of analytical gradient.

        :param molecule:
            The molecule.
        :param basis:
            The AO basis set.
        :param rsp_results:
            The results of the RPA or TDA calculation.
        """

        scf_tensors = self._scf_drv.scf_tensors
        if self._rsp_results is None:
            self._rsp_results = rsp_results

        self.ostream.print_info('Computing orbital response...')
        self.ostream.print_blank()
        self.ostream.flush()

        # compute orbital response
        orbrsp_drv = TddftOrbitalResponse(self.comm, self.ostream)
        orbrsp_drv.update_settings(self.orbrsp_dict, self.method_dict)
        # TODO: also check other options
        if 'state_deriv_index' not in self.orbrsp_dict:
            orbrsp_drv.state_deriv_index = self.state_deriv_index
        if 'timing' not in self.orbrsp_dict:
            orbrsp_drv.timing = self.timing
        if 'filename' not in self.orbrsp_dict:
            orbrsp_drv.filename = self.filename
        orbrsp_drv.compute(molecule, basis, scf_tensors, self._rsp_results)

        omega_ao = orbrsp_drv.compute_omega(molecule, basis, scf_tensors)

        grad_timing = {
            'Relaxed_density': 0.0,
            'Ground_state_grad': 0.0,
            'Kinetic_energy_grad': 0.0,
            'Nuclear_potential_grad': 0.0,
            'Overlap_grad': 0.0,
            'Fock_prep': 0.0,
            'Fock_grad': 0.0,
            'Vxc_grad': 0.0,
            'Fxc_grad': 0.0,
            'Kxc_grad': 0.0,
        }

        t0 = time.time()

        orbrsp_results = orbrsp_drv.cphf_results

        natm = molecule.number_of_atoms()

        if self.rank == mpi_master():
            # only alpha part
            gs_dm = scf_tensors['D_alpha']
            nocc = molecule.number_of_alpha_electrons()
            mo = scf_tensors['C_alpha']
            mo_occ = mo[:, :nocc]
            mo_vir = mo[:, nocc:]
            nocc = mo_occ.shape[1]
            nvir = mo_vir.shape[1]

            # TODO: check variable names and make sure they are consistent
            # with cphfsolver.
            # spin summation already included
            x_plus_y_ao = orbrsp_results['x_plus_y_ao']
            x_minus_y_ao = orbrsp_results['x_minus_y_ao']

            dof = x_plus_y_ao.shape[0]
        else:
            dof = None
        dof = self.comm.bcast(dof, root=mpi_master())

        dist_cphf_ov = orbrsp_results['dist_cphf_ov']

        cphf_ao = []
        for x in range(dof):
            cphf_ov_x = dist_cphf_ov[x].get_full_vector(0)

            if self.rank == mpi_master():
                # CPHF/CPKS coefficients (lambda Lagrange multipliers)
                cphf_ao.append(
                    np.linalg.multi_dot(
                        [mo_occ,
                         cphf_ov_x.reshape(nocc, nvir), mo_vir.T]))

        if self.rank == mpi_master():
            cphf_ao = np.array(cphf_ao)
            unrelaxed_density_ao = orbrsp_results['unrelaxed_density_ao']
            relaxed_density_ao = (unrelaxed_density_ao + 2.0 * cphf_ao +
                                  2.0 * cphf_ao.transpose(0, 2, 1))
        else:
            gs_dm = None
            relaxed_density_ao = None
            x_plus_y_ao = None
            x_minus_y_ao = None

        gs_dm = self.comm.bcast(gs_dm, root=mpi_master())
        relaxed_density_ao = self.comm.bcast(relaxed_density_ao,
                                             root=mpi_master())
        x_plus_y_ao = self.comm.bcast(x_plus_y_ao, root=mpi_master())
        x_minus_y_ao = self.comm.bcast(x_minus_y_ao, root=mpi_master())
        omega_ao = self.comm.bcast(omega_ao, root=mpi_master())

        grad_timing['Relaxed_density'] += time.time() - t0

        # ground state gradient

        self.ostream.print_info('Computing ground-state gradient...')
        self.ostream.print_blank()
        self.ostream.flush()

        t0 = time.time()

        gs_grad_drv = ScfGradientDriver(self._scf_drv)
        gs_grad_drv.update_settings(self.grad_dict, self.method_dict)

        if self.unrelaxed:
            gs_grad_drv.unrelaxed = True

        gs_grad_drv.ostream.mute()
        gs_grad_drv.compute(molecule, basis, scf_tensors)
        gs_grad_drv.ostream.unmute()

        gs_grad = gs_grad_drv.get_gradient()

        grad_timing['Ground_state_grad'] += time.time() - t0

        self.gradient = np.zeros((dof, natm, 3))

        self.ostream.print_info('Computing excited-state gradient...')
        self.ostream.print_blank()
        self.ostream.flush()

        local_atoms = molecule.partition_atoms(self.comm)

        # kinetic energy contribution to gradient

        t0 = time.time()

        kin_grad_drv = KineticEnergyGeom100Driver()

        for iatom in local_atoms:
            gmats = kin_grad_drv.compute(molecule, basis, iatom)

            for i, label in enumerate(['X', 'Y', 'Z']):
                gmat = gmats.matrix_to_numpy(label)
                for s in range(dof):
                    # Sum of alpha + beta already in relaxed_density_ao
                    if self.unrelaxed:
                        self.gradient[s, iatom, i] += np.sum(
                            (gmat + gmat.T) * unrelaxed_density_ao[s])
                    else:
                        self.gradient[s, iatom, i] += np.sum(
                            (gmat + gmat.T) * relaxed_density_ao[s])

            gmats = Matrices()

        grad_timing['Kinetic_energy_grad'] += time.time() - t0

        # nuclear potential contribution to gradient

        t0 = time.time()

        npot_grad_100_drv = NuclearPotentialGeom100Driver()
        npot_grad_010_drv = NuclearPotentialGeom010Driver()

        for iatom in local_atoms:
            gmats_100 = npot_grad_100_drv.compute(molecule, basis, iatom)
            gmats_010 = npot_grad_010_drv.compute(molecule, basis, iatom)

            for i, label in enumerate(['X', 'Y', 'Z']):
                gmat_100 = gmats_100.matrix_to_numpy(label)
                gmat_010 = gmats_010.matrix_to_numpy(label)

                # TODO: move minus sign into function call (such as in oneints)
                for s in range(dof):
                    # summation of alpha and beta already included
                    # in relaxed_density_ao
                    if self.unrelaxed:
                        self.gradient[s, iatom, i] -= np.sum(
                            (gmat_100 + gmat_100.T) * unrelaxed_density_ao[s])
                        self.gradient[s, iatom,
                                      i] -= np.sum(gmat_010 * unrelaxed_density_ao[s])
                    else:
                        self.gradient[s, iatom, i] -= np.sum(
                            (gmat_100 + gmat_100.T) * relaxed_density_ao[s])
                        self.gradient[s, iatom,
                                      i] -= np.sum(gmat_010 * relaxed_density_ao[s])

            gmats_100 = Matrices()
            gmats_010 = Matrices()

        grad_timing['Nuclear_potential_grad'] += time.time() - t0

        if not self.unrelaxed:
            # orbital response contribution to gradient

            t0 = time.time()

            ovl_grad_drv = OverlapGeom100Driver()

            for iatom in local_atoms:
                gmats = ovl_grad_drv.compute(molecule, basis, iatom)

                for i, label in enumerate(['X', 'Y', 'Z']):
                    gmat = gmats.matrix_to_numpy(label)
                    for s in range(dof):
                        self.gradient[s, iatom, i] += 2.0 * np.sum(
                            (gmat + gmat.T) * omega_ao[s])
                gmats = Matrices()

            grad_timing['Overlap_grad'] += time.time() - t0

        # ERI contribution to gradient

        t0 = time.time()

        # determine the Fock type and exchange scaling factor
        if self._dft:
            if self.xcfun.is_hybrid():
                exchange_scaling_factor = self.xcfun.get_frac_exact_exchange()
                fock_type = '2jkx'
            else:
                exchange_scaling_factor = 0.0
                fock_type = 'j'
        else:
            exchange_scaling_factor = 1.0
            fock_type = '2jk'

        need_omega = (self._dft and self.xcfun.is_range_separated())
        if need_omega:
            exchange_scaling_factor = (self.xcfun.get_rs_alpha() +
                                       self.xcfun.get_rs_beta())
            erf_k_coef = -self.xcfun.get_rs_beta()
            omega = self.xcfun.get_rs_omega()
        else:
            erf_k_coef, omega = None, None

        fock_grad_drv = FockGeom1000Driver()
        fock_grad_drv._set_block_size_factor(self._block_size_factor)

        screener = T4CScreener()
        screener.partition(basis, molecule, 'eri')

        grad_timing['Fock_prep'] += time.time() - t0

        t0 = time.time()

        thresh_int = int(-math.log10(self._scf_drv.eri_thresh))

        den_mat_for_fock_gs = make_matrix(basis, mat_t.symmetric)
        den_mat_for_fock_gs.set_values(gs_dm)

        if self._scf_drv.ri_coulomb:

            sym_den_mat_for_fock_rel = make_matrix(basis, mat_t.symmetric)

            sym_den_mat_for_fock_xmy = make_matrix(basis, mat_t.symmetric)
            sym_den_mat_for_fock_xmy_p_xmyT = make_matrix(
                basis, mat_t.symmetric)

            assert_msg_critical(
                basis.get_label().lower().startswith('def2-'),
                'TddftGradientDriver: Invalid basis set for RI-J')

            if self.rank == mpi_master():
                basis_ri_j = MolecularBasis.read(
                    molecule, self._scf_drv.ri_auxiliary_basis)
            else:
                basis_ri_j = None
            basis_ri_j = self.comm.bcast(basis_ri_j, root=mpi_master())

            self._scf_drv._ri_drv.prepare_buffers(molecule,
                                                  basis,
                                                  basis_ri_j,
                                                  verbose=False)

            ri_grad_drv = RIFockGradDriver()

            ri_gvec_gs = self._scf_drv._ri_drv.compute_bq_vector(
                den_mat_for_fock_gs)

            for idx in range(dof):

                if self.unrelaxed:
                    sym_den_mat_for_fock_rel.set_values(
                        self.get_sym_mat(unrelaxed_density_ao[idx]))
                else:
                    sym_den_mat_for_fock_rel.set_values(
                        self.get_sym_mat(relaxed_density_ao[idx]))

                sym_den_mat_for_fock_xmy.set_values(
                    self.get_sym_mat(x_minus_y_ao[idx]))
                sym_den_mat_for_fock_xmy_p_xmyT.set_values(
                    self.get_sym_mat(x_minus_y_ao[idx] + x_minus_y_ao[idx].T))

                ri_gvec_rel = self._scf_drv._ri_drv.compute_bq_vector(
                    sym_den_mat_for_fock_rel)
                ri_gvec_xmy = self._scf_drv._ri_drv.compute_bq_vector(
                    sym_den_mat_for_fock_xmy)
                ri_gvec_xmy_2 = self._scf_drv._ri_drv.compute_bq_vector(
                    sym_den_mat_for_fock_xmy_p_xmyT)

                for iatom in local_atoms:

                    atomgrad_rel = ri_grad_drv.direct_compute(
                        screener, basis, basis_ri_j, molecule, ri_gvec_rel,
                        ri_gvec_gs, sym_den_mat_for_fock_rel,
                        den_mat_for_fock_gs, iatom, thresh_int)

                    atomgrad_xmy = ri_grad_drv.direct_compute(
                        screener, basis, basis_ri_j, molecule, ri_gvec_xmy_2,
                        ri_gvec_xmy, sym_den_mat_for_fock_xmy_p_xmyT,
                        sym_den_mat_for_fock_xmy, iatom, thresh_int)

                    # Note: RI gradient from direct_compute does NOT contain factor of 2
                    self.gradient[idx, iatom, :] += np.array(
                        atomgrad_rel.coordinates()) * 2.0
                    self.gradient[idx, iatom, :] += 0.5 * np.array(
                        atomgrad_xmy.coordinates()) * 2.0

        else:

            den_mat_for_fock_rel = make_matrix(basis, mat_t.general)

            den_mat_for_fock_xpy = make_matrix(basis, mat_t.general)
            den_mat_for_fock_xpy_m_xpyT = make_matrix(basis, mat_t.general)

            den_mat_for_fock_xmy = make_matrix(basis, mat_t.general)
            den_mat_for_fock_xmy_p_xmyT = make_matrix(basis, mat_t.general)

            factor = 2.0 if fock_type == 'j' else 1.0

            for iatom in local_atoms:

                screener_atom = T4CScreener()
                screener_atom.partition_atom(basis, molecule, 'eri', iatom)

                for idx in range(dof):

                    if self.unrelaxed:
                        den_mat_for_fock_rel.set_values(unrelaxed_density_ao[idx])
                    else:
                        den_mat_for_fock_rel.set_values(relaxed_density_ao[idx])

                    den_mat_for_fock_xpy.set_values(x_plus_y_ao[idx])
                    den_mat_for_fock_xpy_m_xpyT.set_values(x_plus_y_ao[idx] -
                                                           x_plus_y_ao[idx].T)

                    den_mat_for_fock_xmy.set_values(x_minus_y_ao[idx])
                    den_mat_for_fock_xmy_p_xmyT.set_values(x_minus_y_ao[idx] +
                                                           x_minus_y_ao[idx].T)

                    atomgrad_rel = fock_grad_drv.compute(
                        basis, screener_atom, screener, den_mat_for_fock_gs,
                        den_mat_for_fock_rel, iatom, fock_type,
                        exchange_scaling_factor, 0.0, thresh_int)

                    atomgrad_xpy = fock_grad_drv.compute(
                        basis, screener_atom, screener, den_mat_for_fock_xpy,
                        den_mat_for_fock_xpy_m_xpyT, iatom, fock_type,
                        exchange_scaling_factor, 0.0, thresh_int)

                    atomgrad_xmy = fock_grad_drv.compute(
                        basis, screener_atom, screener, den_mat_for_fock_xmy,
                        den_mat_for_fock_xmy_p_xmyT, iatom, fock_type,
                        exchange_scaling_factor, 0.0, thresh_int)

                    self.gradient[idx,
                                  iatom, :] += np.array(atomgrad_rel) * factor
                    self.gradient[
                        idx, iatom, :] += 0.5 * np.array(atomgrad_xpy) * factor
                    self.gradient[
                        idx, iatom, :] += 0.5 * np.array(atomgrad_xmy) * factor

                    if need_omega:
                        # for range-separated functional
                        atomgrad_rel_rs = fock_grad_drv.compute(
                            basis, screener_atom, screener, den_mat_for_fock_gs,
                            den_mat_for_fock_rel, iatom, 'kx_rs', erf_k_coef,
                            omega, thresh_int)

                        atomgrad_xpy_rs = fock_grad_drv.compute(
                            basis, screener_atom, screener,
                            den_mat_for_fock_xpy, den_mat_for_fock_xpy_m_xpyT,
                            iatom, 'kx_rs', erf_k_coef, omega, thresh_int)

                        atomgrad_xmy_rs = fock_grad_drv.compute(
                            basis, screener_atom, screener,
                            den_mat_for_fock_xmy, den_mat_for_fock_xmy_p_xmyT,
                            iatom, 'kx_rs', erf_k_coef, omega, thresh_int)

                        self.gradient[idx,
                                      iatom, :] -= np.array(atomgrad_rel_rs)
                        self.gradient[
                            idx, iatom, :] -= 0.5 * np.array(atomgrad_xpy_rs)
                        self.gradient[
                            idx, iatom, :] -= 0.5 * np.array(atomgrad_xmy_rs)

        grad_timing['Fock_grad'] += time.time() - t0

        # TODO: enable multiple DMs for DFT to avoid for-loops.
        if self._dft:
            xcfun_label = self._scf_drv.xcfun.get_func_label()

            xcgrad_drv = XCMolecularGradient()
            mol_grid = self._scf_drv._mol_grid

            for s in range(dof):
                if self.rank == mpi_master():
                    if self.unrelaxed:
                        rhow_dm = 0.5 * unrelaxed_density_ao[s]
                    else:
                        rhow_dm = 0.5 * relaxed_density_ao[s]
                    rhow_dm_sym = 0.5 * (rhow_dm + rhow_dm.T)
                    xmy_sym = 0.5 * (x_minus_y_ao[s] + x_minus_y_ao[s].T)
                else:
                    rhow_dm_sym = None
                    xmy_sym = None

                rhow_dm_sym = self.comm.bcast(rhow_dm_sym, root=mpi_master())
                xmy_sym = self.comm.bcast(xmy_sym, root=mpi_master())

                xc_grad_t0 = time.time()

                tddft_xcgrad = xcgrad_drv.integrate_vxc_gradient(
                    molecule, basis, [rhow_dm_sym], [gs_dm], mol_grid,
                    xcfun_label)

                grad_timing['Vxc_grad'] += time.time() - xc_grad_t0
                xc_grad_t0 = time.time()

                tddft_xcgrad += xcgrad_drv.integrate_fxc_gradient(
                    molecule, basis, [rhow_dm_sym], [gs_dm], [gs_dm], mol_grid,
                    xcfun_label)

                tddft_xcgrad += xcgrad_drv.integrate_fxc_gradient(
                    molecule, basis, [xmy_sym], [xmy_sym], [gs_dm], mol_grid,
                    xcfun_label)

                grad_timing['Fxc_grad'] += time.time() - xc_grad_t0
                xc_grad_t0 = time.time()

                tddft_xcgrad += xcgrad_drv.integrate_kxc_gradient(
                    molecule, basis, [xmy_sym], [xmy_sym], [gs_dm], mol_grid,
                    xcfun_label)

                grad_timing['Kxc_grad'] += time.time() - xc_grad_t0

                self.gradient[s] += tddft_xcgrad

        if self.rank == mpi_master():
            for s in range(dof):
                self.gradient[s] += gs_grad

        if self.do_first_order_prop:
            if self.rank == mpi_master():
                unrelaxed_total_density = (scf_tensors['D_alpha'] +
                                           scf_tensors['D_beta'] +
                                           unrelaxed_density_ao)
                relaxed_total_density = (scf_tensors['D_alpha'] +
                                         scf_tensors['D_beta'] +
                                         relaxed_density_ao)
            else:
                unrelaxed_total_density = None
                relaxed_total_density = None

            self.unrelaxed_dipole_moment = self.compute_dipole_moment(
                molecule,
                basis,
                unrelaxed_total_density,
                dipole_moment_type='Unrelaxed')
            self.relaxed_dipole_moment = self.compute_dipole_moment(
                molecule,
                basis,
                relaxed_total_density,
                dipole_moment_type='Relaxed')

        self.gradient = self.comm.allreduce(self.gradient, op=MPI.SUM)

        if self.timing and self.rank == mpi_master():
            self.ostream.print_info('Gradient timing decomposition')
            for key, val in grad_timing.items():
                self.ostream.print_info(f'    {key:<25s}:  {val:.2f} sec')
            self.ostream.print_blank()

    def compute_dipole_moment(self,
                              molecule,
                              basis,
                              density,
                              dipole_moment_type="Relaxed"):
        """ Computes the dipole moment based on the density.

        :param molecule:
            The molecule.
        :param basis:
            The basis set.
        :param density:
            The total density matrix (alpha + beta).
        :param dipole_moment_type:
            The type of dipole moment (Relaxed or Unrelaxed).

        """

        first_order_prop = FirstOrderProperties(self.comm, self.ostream)
        first_order_prop.compute(molecule, basis, density)

        if self.rank == mpi_master():
            if self.tamm_dancoff:
                method = 'TDA'
            else:
                method = 'RPA'
            title = method + ' ' + dipole_moment_type + ' Dipole Moment(s) '
            first_order_prop.print_properties(molecule, title,
                                              self.state_deriv_index)
            self.ostream.print_blank()
            dipole_moment = first_order_prop.get_property('dipole moment')
            return dipole_moment
        else:
            return None

    def compute_energy(self, molecule, basis, scf_drv, rsp_drv, rsp_results):
        """
        Computes the energy at the current position.

        :param molecule:
            The molecule.
        :param basis:
            The basis set.
        :param scf_drv:
            The SCF driver.
        :param rsp_drv:
            The linear response driver.
        :param rsp_results:
            The dictionary containing response results.
        """

        if self.rank == mpi_master():
            if isinstance(self.state_deriv_index, int):
                state_deriv_index = self.state_deriv_index - 1
            else:
                state_deriv_index = self.state_deriv_index[0] - 1
        else:
            state_deriv_index = None
        state_deriv_index = self.comm.bcast(state_deriv_index,
                                            root=mpi_master())

        if self.numerical:
            # disable restarting scf for numerical gradient
            scf_drv.restart = False
        else:
            # always try restarting scf for analytical gradient
            scf_drv.restart = True
        scf_results = scf_drv.compute(molecule, basis)
        assert_msg_critical(scf_drv.is_converged,
                            'TddftGradientDriver: SCF did not converge')
        self._scf_drv = scf_drv

        # response should not be restarted
        rsp_drv.restart = False
        rsp_results = rsp_drv.compute(molecule, basis, scf_results)
        assert_msg_critical(rsp_drv.is_converged,
                            'TddftGradientDriver: response did not converge')
        self._rsp_results = rsp_results

        if self.rank == mpi_master():
            scf_ene = scf_results['scf_energy']
            exc_ene = rsp_results['eigenvalues'][state_deriv_index]
            total_ene = scf_ene + exc_ene
        else:
            total_ene = None
        total_ene = self.comm.bcast(total_ene, root=mpi_master())

        return total_ene

    def compute_numerical_dipole(self,
                                 molecule,
                                 ao_basis,
                                 scf_drv,
                                 rsp_drv,
                                 field_strength=1e-5):
        """
        Performs calculation of numerical dipole moment at RPA or TDA level.

        :param molecule:
            The molecule.
        :param ao_basis:
            The AO basis set.
        :param scf_drv:
            The SCF driver.
        :param rsp_drv:
            The RPA or TDA driver.
        :param field_strength:
            The strength of the external electric field.
        :return:
            The electric dipole moment vector.
        """

        scf_drv.ostream.mute()
        rsp_drv.ostream.mute()

        # numerical dipole moment of the excited states
        n_states = len(self.state_deriv_index)
        dipole_moment = np.zeros((n_states, 3))
        field = [0.0, 0.0, 0.0]

        for s in range(n_states):

            self.ostream.unmute()
            self.ostream.print_info(
                f'Processing excited state {s + 1}/{n_states}...')
            self.ostream.flush()
            self.ostream.mute()

            for i in range(3):
                field[i] = field_strength
                scf_drv.electric_field = field
                scf_drv.compute(molecule, ao_basis)
                scf_tensors = scf_drv.scf_tensors
                rsp_drv._is_converged = False  # only needed for RPA
                rsp_results = rsp_drv.compute(molecule, ao_basis, scf_tensors)
                exc_en_plus = rsp_results['eigenvalues'][s]
                e_plus = scf_drv.get_scf_energy() + exc_en_plus

                field[i] = -field_strength
                scf_drv.compute(molecule, ao_basis)
                rsp_drv._is_converged = False
                rsp_results = rsp_drv.compute(molecule, ao_basis,
                                              scf_drv.scf_tensors)
                exc_en_minus = rsp_results['eigenvalues'][s]
                e_minus = scf_drv.get_scf_energy() + exc_en_minus

                field[i] = 0.0
                dipole_moment[s, i] = (-(e_plus - e_minus) /
                                       (2.0 * field_strength))

        scf_drv.ostream.unmute()
        rsp_drv.ostream.unmute()

        return dipole_moment

    @staticmethod
    def get_sym_mat(array):

        return 0.5 * (array + array.T)
