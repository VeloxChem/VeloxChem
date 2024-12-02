#
#                              VELOXCHEM
#         ----------------------------------------------------
#                     An Electronic Structure Code
#
#  Copyright © 2018-2024 by VeloxChem developers. All rights reserved.
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
from copy import deepcopy
from os import environ
import numpy as np
import time
import math

from .veloxchemlib import (OverlapGeom100Driver, KineticEnergyGeom100Driver,
                           NuclearPotentialGeom100Driver,
                           NuclearPotentialGeom010Driver, FockGeom1000Driver)
from .veloxchemlib import mpi_master, mat_t
from .veloxchemlib import denmat
from .veloxchemlib import XCIntegrator
from .veloxchemlib import T4CScreener
from .veloxchemlib import partition_atoms, make_matrix
from .matrices import Matrices
from .profiler import Profiler
from .tddftorbitalresponse import TddftOrbitalResponse
from .molecule import Molecule
from .lrsolver import LinearResponseSolver
from .gradientdriver import GradientDriver
from .scfgradientdriver import ScfGradientDriver
from .errorhandler import assert_msg_critical
from .inputparser import (parse_input, parse_seq_fixed)
from .sanitychecks import dft_sanity_check

# For PySCF integral derivatives
from .import_from_pyscf import overlap_deriv
from .import_from_pyscf import hcore_deriv
from .import_from_pyscf import eri_deriv

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
        self._block_size_factor = scf_drv._block_size_factor
        self._scf_drv = scf_drv

        self.flag = 'RPA Gradient Driver'

        # flag on whether RPA or TDA is calculated
        self.tamm_dancoff = False

        # excited state information; if it is set to None, 
        # all available states will be calculated.
        self.state_deriv_index = None

        self.do_first_order_prop = False
        self.relaxed_dipole_moment = None

        # TODO: remove relaxed_dipole_moment (not input variable)
        self._input_keywords['gradient'].update({
            'tamm_dancoff': ('bool', 'whether RPA or TDA is calculated'),
            'state_deriv_index': ('seq_fixed_int', 'excited state information'),
            'do_first_order_prop': ('bool', 'do first-order property'),
            'relaxed_dipole_moment': ( 'float','relaxed excited-state dipole moment'),
            }
        )

    def update_settings(self, grad_dict, rsp_dict, orbrsp_dict=None, method_dict=None):
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
            orbrsp_dict['do_first_order_prop'] = grad_dict['do_first_order_prop']

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
        dft_sanity_check(self, 'compute')

        if self.rank == mpi_master():
            all_states = list(np.arange(1, len(rsp_results['eigenvalues'])+1))
            if self.state_deriv_index is not None:
                error_message =  'TdscfGradientDriver: some of the '
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

        start_time = time.time()

        # NOTE: the numerical gradient is calculated for the first state only.
        if self.numerical:
            scf_drv.ostream.mute()
            self.compute_numerical(molecule, basis, scf_drv, rsp_drv, rsp_results)
            scf_drv.ostream.unmute()
        else:
            self.compute_analytical(molecule, basis, scf_drv, rsp_results)

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

    def compute_analytical_vlx(self, molecule, basis, rsp_results):
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

        # compute orbital response
        orbrsp_drv = TddftOrbitalResponse(self.comm, self.ostream)
        orbrsp_drv.update_settings(self.orbrsp_dict, self.method_dict)
        orbrsp_drv.compute(molecule, basis, scf_tensors,
                           rsp_results)

        orbrsp_results = orbrsp_drv.cphf_results
        omega_ao = orbrsp_drv.compute_omega(molecule, basis, scf_tensors)

        if self.rank == mpi_master():
            # only alpha part
            gs_dm = scf_tensors['D_alpha']
            nocc = molecule.number_of_alpha_electrons()
            natm = molecule.number_of_atoms()
            mo = scf_tensors['C_alpha']
            mo_occ = mo[:, :nocc]
            mo_vir = mo[:, nocc:]
            nocc = mo_occ.shape[1]
            nvir = mo_vir.shape[1]
            nao = mo_occ.shape[0]

            # TODO: check variable names and make sure they are consistent
            # with cphfsolver.
            # spin summation already included
            x_plus_y_ao = orbrsp_results['x_plus_y_ao']
            x_minus_y_ao = orbrsp_results['x_minus_y_ao']

            # CPHF/CPKS coefficients (lambda Lagrange multipliers)
            cphf_ov = orbrsp_results['cphf_ov']
            unrelaxed_density_ao = orbrsp_results['unrelaxed_density_ao']
            dof = x_plus_y_ao.shape[0]
            cphf_ao = np.array([
                np.linalg.multi_dot([mo_occ, cphf_ov[x], mo_vir.T])
                for x in range(dof)
            ])
            relaxed_density_ao = ( unrelaxed_density_ao + 2.0 * cphf_ao
                        + 2.0 * cphf_ao.transpose(0,2,1) )
        else:
            dof = None
            relaxed_density_ao = None

        dof = self.comm.bcast(dof, root=mpi_master())
        relaxed_density_ao = self.comm.bcast(relaxed_density_ao,
                                             root=mpi_master())

        # ground state gradient
        t1 = time.time()
        gs_grad_drv = ScfGradientDriver(self._scf_drv)
        gs_grad_drv.update_settings(self.grad_dict, self.method_dict)

        gs_grad_drv.ostream.mute()
        gs_grad_drv.compute(molecule, basis, scf_tensors)
        gs_grad_drv.ostream.unmute()
        t2 = time.time()

        if self.rank == mpi_master():
            if gs_grad_drv.numerical:
                gs_gradient_type = "Numerical GS gradient"
            else:
                gs_gradient_type = "Analytical GS gradient"
            self.ostream.print_info(gs_gradient_type
                                    + ' computed in'
                                    + ' {:.2f} sec.'.format(t2 - t1))
            self.ostream.print_blank()
            self.ostream.flush()

            gs_gradient = gs_grad_drv.get_gradient()

        self.gradient = np.zeros((dof, natm, 3))

        local_atoms = partition_atoms(natm, self.rank, self.nodes)

        # kinetic energy contribution to gradient

        kin_grad_drv = KineticEnergyGeom100Driver()

        self._print_debug_info('before kin_grad')

        for iatom in local_atoms:
            gmats = kin_grad_drv.compute(molecule, basis, iatom)

            for i, label in enumerate(['X', 'Y', 'Z']):
                gmat = gmats.matrix_to_numpy(label)
                for s in range(dof):
                    # Sum of alpha + beta already in relaxed_density_ao
                    self.gradient[s, iatom, i] += np.sum((gmat + gmat.T) 
                                                 * relaxed_density_ao[s])

            gmats = Matrices()

        self._print_debug_info('after  kin_grad')

        # nuclear potential contribution to gradient

        self._print_debug_info('before npot_grad')

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
                    self.gradient[s, iatom, i] -=  np.sum(
                        (gmat_100 + gmat_100.T) * relaxed_density_ao[s])
                    self.gradient[s, iatom, i] -=  np.sum(gmat_010 
                                                * relaxed_density_ao[s])

            gmats_100 = Matrices()
            gmats_010 = Matrices()

        self._print_debug_info('after  npot_grad')
        
        # orbital contribution to gradient

        self._print_debug_info('before ovl_grad')

        ovl_grad_drv = OverlapGeom100Driver()

        for iatom in local_atoms:
            gmats = ovl_grad_drv.compute(molecule, basis, iatom)

            for i, label in enumerate(['X', 'Y', 'Z']):
                gmat = gmats.matrix_to_numpy(label)
                for s in range(dof):
                    self.gradient[s, iatom, i] += 2.0 * np.sum(
                        (gmat + gmat.T) * omega_ao[s])
            gmats = Matrices()

        self._print_debug_info('after  ovl_grad')

        # ERI contribution to gradient

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

        # TODO: range-separated Fock
        need_omega = (self._dft and self.xcfun.is_range_separated())
        if need_omega:
            assert_msg_critical(
                False, 'ScfGradientDriver: Not implemented for' +
                ' range-separated functional')

        fock_timing = {
            'Screening': 0.0,
            'FockGrad': 0.0,
            }

        self._print_debug_info('before fock_grad')
        fock_grad_drv = FockGeom1000Driver()

        extra_factor = self._get_extra_block_size_factor(
            basis.get_dimensions_of_basis())
        fock_grad_drv._set_block_size_factor(self._block_size_factor *
                                             extra_factor)
        t0 = time.time()

        screener = T4CScreener()
        screener.partition(basis, molecule, 'eri')

        fock_timing['Screening'] += time.time() - t0

        thresh_int = int(-math.log10(self._scf_drv.eri_thresh))

        factor = 2.0 if fock_type == 'j' else 1.0

        for iatom in local_atoms:

            screener_atom = T4CScreener()
            screener_atom.partition_atom(basis, molecule, 'eri', iatom)

            t0 = time.time()

            for idx in range(dof):
                den_mat_for_fock_rel = make_matrix(basis, mat_t.symmetric)
                den_mat_for_fock_rel.set_values(relaxed_density_ao[idx])
                den_mat_for_fock_rel2 = make_matrix(basis, mat_t.general)
                den_mat_for_fock_rel2.set_values(relaxed_density_ao[idx])

                atomgrad = fock_grad_drv.compute(basis, screener_atom, screener,
                                             den_mat_for_fock_rel,
                                             den_mat_for_fock_rel2, iatom,
                                             fock_type, exchange_scaling_factor,
                                             0.0, thresh_int)
                
                self.gradient[idx, iatom, :] += np.array(atomgrad) * factor

        fock_timing['FockGrad'] += time.time() - t0

        self._print_debug_info('after  fock_grad')


    # TODO: remove this routine once compute_analytical_vlx
    # confirmed to be correct.
    def compute_analytical(self, molecule, basis, scf_drv, rsp_results):
        """
        Performs calculation of analytical gradient.

        :param molecule:
            The molecule.
        :param basis:
            The AO basis set.
        :param scf_drv:
            The SCF driver.
        :param rsp_results:
            The results of the RPA or TDA calculation.
        """

        # SCF results
        scf_tensors = scf_drv.scf_tensors

        # compute orbital response
        orbrsp_drv = TddftOrbitalResponse(self.comm, self.ostream)
        orbrsp_drv.update_settings(self.orbrsp_dict, self.method_dict)
        orbrsp_drv.compute(molecule, basis, scf_tensors,
                           rsp_results)

        orbrsp_results = orbrsp_drv.cphf_results
        omega_ao = orbrsp_drv.compute_omega(molecule, basis, scf_tensors)

        if self.rank == mpi_master():
            # only alpha part
            gs_dm = scf_drv.scf_tensors['D_alpha']
            nocc = molecule.number_of_alpha_electrons()
            natm = molecule.number_of_atoms()
            mo = scf_tensors['C_alpha']
            mo_occ = mo[:, :nocc]
            mo_vir = mo[:, nocc:]
            nocc = mo_occ.shape[1]
            nvir = mo_vir.shape[1]
            nao = mo_occ.shape[0]

            # TODO: check variable names and make sure they are consistent
            # with cphfsolver.
            # spin summation already included
            x_plus_y_ao = orbrsp_results['x_plus_y_ao']
            x_minus_y_ao = orbrsp_results['x_minus_y_ao']

            # CPHF/CPKS coefficients (lambda Lagrange multipliers)
            cphf_ov = orbrsp_results['cphf_ov']
            unrelaxed_density_ao = orbrsp_results['unrelaxed_density_ao']
            dof = x_plus_y_ao.shape[0]
            cphf_ao = np.array([
                np.linalg.multi_dot([mo_occ, cphf_ov[x], mo_vir.T])
                for x in range(dof)
            ])
            relaxed_density_ao = ( unrelaxed_density_ao + 2.0 * cphf_ao
                        + 2.0 * cphf_ao.transpose(0,2,1) )
        else:
            dof = None

        dof = self.comm.bcast(dof, root=mpi_master())

        # ground state gradient
        t1 = time.time()
        gs_grad_drv = ScfGradientDriver(scf_drv)
        gs_grad_drv.update_settings(self.grad_dict, self.method_dict)

        gs_grad_drv.ostream.mute()
        gs_grad_drv.compute(molecule, basis, scf_drv.scf_tensors)
        gs_grad_drv.ostream.unmute()
        t2 = time.time()

        if self.rank == mpi_master():
            if gs_grad_drv.numerical:
                gs_gradient_type = "Numerical GS gradient"
            else:
                gs_gradient_type = "Analytical GS gradient"
            self.ostream.print_info(gs_gradient_type
                                    + ' computed in'
                                    + ' {:.2f} sec.'.format(t2 - t1))
            self.ostream.print_blank()
            self.ostream.flush()

            self.gradient = np.zeros((dof, natm, 3))
            self.gradient += gs_grad_drv.get_gradient()


            self.partial_contrib = np.zeros((dof, natm, 3))

            # loop over atoms and contract integral derivatives
            # with density matrices
            # add the corresponding contribution to the gradient
            for i in range(natm):
                t0 = time.time()
                # taking integral derivatives from pyscf
                d_ovlp = overlap_deriv(molecule, basis, i)
                d_hcore = hcore_deriv(molecule, basis, i)
                d_eri = eri_deriv(molecule, basis, i)

                if self._dft:
                    if self.xcfun.is_hybrid():
                        frac_K = self.xcfun.get_frac_exact_exchange()
                    else:
                        frac_K = 0.0
                else:
                    frac_K = 1.0

                for s in range(dof):
                    for x in range(3):
                        self.gradient[s,i,x] += np.linalg.multi_dot([
                            relaxed_density_ao[s].reshape(nao**2),
                            (d_hcore[x].T).reshape(nao**2)
                        ])
                        self.partial_contrib[s,i,x] += np.linalg.multi_dot([
                            relaxed_density_ao[s].reshape(nao**2),
                            (d_hcore[x].T).reshape(nao**2)
                        ])

                        self.gradient[s,i,x] += np.linalg.multi_dot([
                            2.0 * omega_ao[s].reshape(nao**2),
                            (d_ovlp[x].T).reshape(nao**2)
                        ])
                        self.partial_contrib[s,i,x] += np.linalg.multi_dot([
                            2.0 * omega_ao[s].reshape(nao**2),
                            (d_ovlp[x].T).reshape(nao**2)
                        ])
                        self.gradient[s,i,x] += 2.0 * np.linalg.multi_dot([
                            gs_dm.reshape(nao**2),
                            d_eri[x].transpose(1,0,3,2).reshape(nao**2,nao**2),
                            relaxed_density_ao[s].reshape(nao**2)
                        ])
                        self.partial_contrib[s,i,x] += 2.0 * np.linalg.multi_dot([
                            gs_dm.reshape(nao**2),
                            d_eri[x].transpose(1,0,3,2).reshape(nao**2,nao**2),
                            relaxed_density_ao[s].reshape(nao**2)
                        ])
                        self.gradient[s,i,x] += -1.0 * frac_K * np.linalg.multi_dot([
                           gs_dm.reshape(nao**2),
                           d_eri[x].transpose(3,0,2,1).reshape(nao**2,nao**2),
                           relaxed_density_ao[s].reshape(nao**2)
                        ])
                        self.partial_contrib[s,i,x] += -1.0 * frac_K * np.linalg.multi_dot([
                           gs_dm.reshape(nao**2),
                           d_eri[x].transpose(3,0,2,1).reshape(nao**2,nao**2),
                           relaxed_density_ao[s].reshape(nao**2)
                        ])
                        self.gradient[s,i,x] += np.linalg.multi_dot([
                            x_plus_y_ao[s].reshape(nao**2),
                            d_eri[x].transpose(3,2,0,1).reshape(nao**2,nao**2),
                            (x_plus_y_ao[s] - x_plus_y_ao[s].T).reshape(nao**2)
                        ])
                        self.gradient[s,i,x] += - 0.5 * frac_K * np.linalg.multi_dot([
                            x_plus_y_ao[s].reshape(nao**2),
                            d_eri[x].transpose(1,2,0,3).reshape(nao**2,nao**2),
                            (x_plus_y_ao[s] - x_plus_y_ao[s].T).reshape(nao**2)
                        ]) 
                        self.gradient[s,i,x] += np.linalg.multi_dot([
                            x_minus_y_ao[s].reshape(nao**2),
                            d_eri[x].transpose(3,2,0,1).reshape(nao**2,nao**2),
                            (x_minus_y_ao[s] + x_minus_y_ao[s].T).reshape(nao**2)
                        ])
                        self.gradient[s,i,x] -= 0.5 * frac_K * np.linalg.multi_dot([
                            x_minus_y_ao[s].reshape(nao**2),
                            d_eri[x].transpose(1,2,0,3).reshape(nao**2,nao**2),
                            (x_minus_y_ao[s] + x_minus_y_ao[s].T).reshape(nao**2)
                        ])
                        

        # TODO: enable multiple DMs for DFT to avoid for-loops.
        if self._dft:
            xcfun_label = scf_drv.xcfun.get_func_label()

            for s in range(dof):
                if self.rank == mpi_master():
                    gs_dm = scf_drv.scf_tensors['D_alpha']

                    rhow_dm = 0.5 * relaxed_density_ao[s]
                    rhow_dm_sym = 0.5 * (rhow_dm + rhow_dm.T)

                    x_minus_y_sym = 0.5 * (x_minus_y_ao[s] + x_minus_y_ao[s].T)
                else:
                    gs_dm = None
                    rhow_dm_sym = None
                    x_minus_y_sym = None

                gs_dm = self.comm.bcast(gs_dm, root=mpi_master())
                rhow_dm_sym = self.comm.bcast(rhow_dm_sym, root=mpi_master())
                x_minus_y_sym = self.comm.bcast(x_minus_y_sym, root=mpi_master())

                tddft_xcgrad = self.grad_tddft_xc_contrib(molecule, basis,
                                               [rhow_dm_sym], [x_minus_y_sym],
                                               [gs_dm], xcfun_label)

                if self.rank == mpi_master():
                    self.gradient[s] += tddft_xcgrad

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
        :param state_deriv_index:
            The index of the excited state of interest.
        """
        if self.rank == mpi_master():
            if isinstance(self.state_deriv_index, int):
                # Python numbering starts at 0
                state_deriv_index = self.state_deriv_index - 1
            else:
                state_deriv_index = self.state_deriv_index[0] - 1
        scf_drv.restart = False
        scf_results = scf_drv.compute(molecule, basis)
        assert_msg_critical(scf_drv.is_converged,
                            'TddftGradientDriver: SCF did not converge')

        rsp_drv.restart = False
        rsp_results = rsp_drv.compute(molecule, basis, scf_results)
        assert_msg_critical(rsp_drv.is_converged,
                            'TddftGradientDriver: response did not converge')

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

        # numerical dipole moment of the excited states
        n_states = len(self.state_deriv_index)
        dipole_moment = np.zeros((n_states, 3))
        field = [0.0, 0.0, 0.0]

        for s in range(n_states):
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
                dipole_moment[s, i] = ( -(e_plus - e_minus)
                                         / (2.0 * field_strength) )

        scf_drv.ostream.unmute()

        return dipole_moment

    # TODO: this routine is used by both ScfGradientDriver and
    # TddftGradientDriver. It can be moved to the parent class.
    def _print_debug_info(self, label):
        """
        Prints debug information.

        :param label:
            The label of debug information.
        """

        if self._debug:
            profiler = Profiler()
            self.ostream.print_info(f'==DEBUG==   available memory {label}: ' +
                                    profiler.get_available_memory())
            self.ostream.flush()

    # TODO: this routine is used by both ScfGradientDriver and
    # TddftGradientDriver. It can be moved to the parent class.
    def _get_extra_block_size_factor(self, naos):

        total_cores = self.nodes * int(environ['OMP_NUM_THREADS'])

        if total_cores >= 2048:
            if naos >= 4500:
                extra_factor = 4
            else:
                extra_factor = 2

        elif total_cores >= 1024:
            extra_factor = 2

        else:
            extra_factor = 1

        return extra_factor
