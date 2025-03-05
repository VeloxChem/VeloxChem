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

import numpy as np
import time as tm
import math

from .veloxchemlib import T4CScreener
from .veloxchemlib import XCMolecularHessian
from .veloxchemlib import mpi_master
from .veloxchemlib import make_matrix, mat_t, partition_atoms
from .veloxchemlib import (OverlapGeom100Driver, KineticEnergyGeom100Driver,
                           NuclearPotentialGeom100Driver,
                           NuclearPotentialGeom010Driver, FockGeom1000Driver)
from .matrices import Matrices
from .profiler import Profiler
from .distributedarray import DistributedArray
from .cphfsolver import CphfSolver
from .errorhandler import assert_msg_critical
from .dftutils import get_default_grid_level
from .batchsize import get_batch_size


class HessianOrbitalResponse(CphfSolver):
    """
    Implements solver for the coupled-perturbed Hartree-Fock (CPHF) equations
    for the Hessian orbital response

    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.

    Instance variables
        - use_subspace_solver: flag to use subspace solver
          instead of conjugate gradient.
    """

    def __init__(self, comm=None, ostream=None):
        """
        Initializes orbital response driver to default setup.

        """

        super().__init__(comm, ostream)

        self.orbrsp_type = 'hessian'
        self._embedding_hess_drv = None

    def update_settings(self, cphf_dict, method_dict=None):
        """
        Updates response and method settings in CPHF solver.

        :param cphf_dict:
            The dictionary of CPHF (orbital response) settings.
        :param method_dict:
            The dictionary of method settings.
        """

        super().update_settings(cphf_dict, method_dict)

    def compute_rhs(self, molecule, basis, scf_tensors, eri_dict, dft_dict,
                    pe_dict):
        """
        Computes the right hand side for the CPHF equations for
        the analytical Hessian, all atomic coordinates.

        :param molecule:
            The molecule.
        :param basis:
            The AO basis set.
        :param scf_tensors:
            The tensors from the converged SCF calculation.

        :returns:
            The RHS of the CPHF equations.
        """

        profiler = Profiler({
            'timing': self.timing,
            'profiling': self.profiling,
            'memory_profiling': self.memory_profiling,
            'memory_tracing': self.memory_tracing,
        })

        profiler.set_timing_key('RHS')

        natm = molecule.number_of_atoms()
        nocc = molecule.number_of_alpha_electrons()

        if self.rank == mpi_master():
            density = scf_tensors['D_alpha']
            eocc = scf_tensors['E_alpha'][:nocc]
            mo = scf_tensors['C_alpha']
            point_charges = scf_tensors.get('point_charges', None)
            qm_vdw_params = scf_tensors.get('qm_vdw_params', None)
        else:
            density = None
            eocc = None
            mo = None
            point_charges = None
            qm_vdw_params = None

        density = self.comm.bcast(density, root=mpi_master())
        eocc = self.comm.bcast(eocc, root=mpi_master())
        mo = self.comm.bcast(mo, root=mpi_master())
        point_charges = self.comm.bcast(point_charges, root=mpi_master())
        qm_vdw_params = self.comm.bcast(qm_vdw_params, root=mpi_master())

        # TODO: double check dft_dict['gs_density']
        mol_grid = dft_dict['molgrid']
        gs_density = [density]

        # MO coefficients
        mo_occ = mo[:, :nocc].copy()
        mo_vir = mo[:, nocc:].copy()
        nvir = mo_vir.shape[1]

        omega_ao = -1.0 * np.linalg.multi_dot([mo_occ, np.diag(eocc), mo_occ.T])

        # partition atoms for parallellisation
        # TODO: use partition_atoms in e.g. scfgradientdriver
        local_atoms = partition_atoms(natm, self.rank, self.nodes)

        atom_idx_rank = [(iatom, self.rank) for iatom in local_atoms]
        gathered_atom_idx_rank = self.comm.gather(atom_idx_rank)
        if self.rank == mpi_master():
            all_atom_idx_rank = []
            for local_list in gathered_atom_idx_rank:
                for elem in local_list:
                    all_atom_idx_rank.append(elem)
            all_atom_idx_rank = sorted(all_atom_idx_rank)
        else:
            all_atom_idx_rank = None
        all_atom_idx_rank = self.comm.bcast(all_atom_idx_rank,
                                            root=mpi_master())

        # preparing the CPHF RHS

        t0 = tm.time()

        # overlap derivatives

        profiler.start_timer('dOvlp')

        ovlp_deriv_ao_dict = {
            (iatom, x): None for iatom in range(natm) for x in range(3)
        }

        ovlp_grad_drv = OverlapGeom100Driver()

        for iatom in local_atoms:

            gmats = ovlp_grad_drv.compute(molecule, basis, iatom)

            for x, label in enumerate(['X', 'Y', 'Z']):
                gmat = gmats.matrix_to_numpy(label)
                ovlp_deriv_ao_dict[(iatom, x)] = gmat + gmat.T

            gmats = Matrices()

        dist_ovlp_deriv_ao = {}

        for iatom, root_rank in all_atom_idx_rank:
            for x in range(3):
                key_ix = (iatom, x)

                dist_ovlp_deriv_ao[key_ix] = DistributedArray(
                    ovlp_deriv_ao_dict[key_ix], self.comm, root=root_rank)

                ovlp_deriv_ao_dict[key_ix] = None

        ovlp_deriv_ao_dict.clear()

        profiler.stop_timer('dOvlp')

        # Fock derivatives

        profiler.start_timer('dFock')

        fock_deriv_ao_dict = {
            (iatom, x): None for iatom in range(natm) for x in range(3)
        }

        for iatom in local_atoms:

            fock_deriv_ao_i = self._compute_fmat_deriv(molecule, basis, density,
                                                       iatom, eri_dict,
                                                       point_charges,
                                                       qm_vdw_params)

            fock_deriv_ao_dict[(iatom, 0)] = fock_deriv_ao_i[0]
            fock_deriv_ao_dict[(iatom, 1)] = fock_deriv_ao_i[1]
            fock_deriv_ao_dict[(iatom, 2)] = fock_deriv_ao_i[2]

        dist_fock_deriv_ao = {}

        for iatom, root_rank in all_atom_idx_rank:
            for x in range(3):
                key_ix = (iatom, x)

                dist_fock_deriv_ao[key_ix] = DistributedArray(
                    fock_deriv_ao_dict[key_ix], self.comm, root=root_rank)

                fock_deriv_ao_dict[key_ix] = None

        fock_deriv_ao_dict.clear()

        profiler.stop_timer('dFock')

        # XC derivatives

        # Note: Parallelization of DFT integration is done over grid points
        # instead of atoms.

        profiler.start_timer('dXC')

        if self._dft:
            xc_mol_hess = XCMolecularHessian()

            naos = basis.get_dimensions_of_basis()

            batch_size = get_batch_size(None, natm * 3, naos, self.comm)
            batch_size = batch_size // 3

            num_batches = natm // batch_size
            if natm % batch_size != 0:
                num_batches += 1

            for batch_ind in range(num_batches):

                batch_start = batch_ind * batch_size
                batch_end = min(batch_start + batch_size, natm)

                atom_list = list(range(batch_start, batch_end))

                vxc_deriv_batch = xc_mol_hess.integrate_vxc_fock_gradient(
                    molecule, basis, gs_density, mol_grid,
                    self.xcfun.get_func_label(), atom_list)

                for vecind, (iatom, root_rank) in enumerate(
                        all_atom_idx_rank[batch_start:batch_end]):

                    for x in range(3):
                        vxc_deriv_ix = self.comm.reduce(
                            vxc_deriv_batch[vecind * 3 + x])

                        dist_vxc_deriv_ix = DistributedArray(
                            vxc_deriv_ix, self.comm)

                        key_ix = (iatom, x)
                        dist_fock_deriv_ao[
                            key_ix].data += dist_vxc_deriv_ix.data

        profiler.stop_timer('dXC')

        t1 = tm.time()

        if self.rank == mpi_master():
            self.ostream.print_info("CPHF/CPKS integral derivatives computed" +
                                    ' in {:.2f} sec.'.format(t1 - t0))
            self.ostream.print_blank()
            self.ostream.flush()

        hessian_first_integral_derivatives = np.zeros((natm, 3, natm, 3))

        profiler.start_timer('1stHess')

        for iatom, root_rank_i in all_atom_idx_rank:
            for x in range(3):
                key_ix = (iatom, x)

                ovlp_deriv_ao_ix = dist_ovlp_deriv_ao[key_ix].get_full_matrix()
                fock_deriv_ao_ix = dist_fock_deriv_ao[key_ix].get_full_matrix()

                if self.rank == mpi_master():
                    DFD_ix = np.linalg.multi_dot(
                        [density, fock_deriv_ao_ix, density])
                    DSD_ix = np.linalg.multi_dot(
                        [density, ovlp_deriv_ao_ix, density])
                    OmegaSD_ix = np.linalg.multi_dot(
                        [omega_ao, ovlp_deriv_ao_ix, density])
                else:
                    DFD_ix = None
                    DSD_ix = None
                    OmegaSD_ix = None

                dist_DFD_ix = DistributedArray(DFD_ix, self.comm)
                dist_DSD_ix = DistributedArray(DSD_ix, self.comm)
                dist_OmegaSD_ix = DistributedArray(OmegaSD_ix, self.comm)

                for jatom, root_rank_j in all_atom_idx_rank:
                    if jatom < iatom:
                        continue
                    for y in range(3):
                        key_jy = (jatom, y)

                        Fix_Sjy = np.dot(
                            dist_DFD_ix.data.reshape(-1),
                            dist_ovlp_deriv_ao[key_jy].data.reshape(-1))

                        Fjy_Six = np.dot(
                            dist_DSD_ix.data.reshape(-1),
                            dist_fock_deriv_ao[key_jy].data.reshape(-1))

                        Six_Sjy = 2.0 * np.dot(
                            dist_OmegaSD_ix.data.reshape(-1),
                            dist_ovlp_deriv_ao[key_jy].data.reshape(-1))

                        hess_ijxy = -2.0 * (Fix_Sjy + Fjy_Six + Six_Sjy)
                        hessian_first_integral_derivatives[iatom, x, jatom,
                                                           y] += hess_ijxy
                        if iatom != jatom:
                            hessian_first_integral_derivatives[jatom, y, iatom,
                                                               x] += hess_ijxy

        hessian_first_integral_derivatives = self.comm.reduce(
            hessian_first_integral_derivatives, root=mpi_master())

        profiler.stop_timer('1stHess')

        dist_cphf_rhs = []

        hessian_eri_overlap = np.zeros((natm, 3, natm, 3))

        for iatom, root_rank in all_atom_idx_rank:

            uij_t0 = tm.time()

            # the oo part of the CPHF coefficients in AO basis
            if self.rank == mpi_master():
                uij_ao_list = []
            else:
                uij_ao_list = None

            for x in range(3):
                key_ix = (iatom, x)
                ovlp_deriv_ao_ix = dist_ovlp_deriv_ao[key_ix].get_full_matrix()

                if self.rank == mpi_master():
                    uij_ao_list.append(
                        np.linalg.multi_dot(
                            [density, -0.5 * ovlp_deriv_ao_ix, density]))

            profiler.add_timing_info('UijAO', tm.time() - uij_t0)

            # create AODensity and Fock matrix objects, contract with ERI
            fock_uij = self._comp_lr_fock(uij_ao_list, molecule, basis,
                                          eri_dict, dft_dict, pe_dict, profiler)

            if self.rank == mpi_master():
                uij_ao_list.clear()

            uij_t0 = tm.time()

            for x in range(3):
                if self.rank == mpi_master():
                    DFD_ix = np.linalg.multi_dot(
                        [density, fock_uij[x], density])
                else:
                    DFD_ix = None

                dist_DFD_ix = DistributedArray(DFD_ix, self.comm)

                # Note: hessian_eri_overlap, i.e. inner product of P_P_Six_fock_ao and
                # P_P_Sjy, is obtained from fock_uij and uij_ao_jy in hessian orbital
                # response

                for jatom, root_rank_j in all_atom_idx_rank:
                    if jatom < iatom:
                        continue

                    for y in range(3):
                        key_jy = (jatom, y)

                        # 2.0 * (np.dot(DFD_ix_list[x], (-0.5 * ovlp_deriv_ao_jy)))
                        hess_ijxy = -np.dot(
                            dist_DFD_ix.data.reshape(-1),
                            dist_ovlp_deriv_ao[key_jy].data.reshape(-1))

                        hess_ijxy *= 4.0

                        hessian_eri_overlap[iatom, x, jatom, y] += hess_ijxy
                        if iatom != jatom:
                            hessian_eri_overlap[jatom, y, iatom, x] += hess_ijxy

            profiler.add_timing_info('UijDot', tm.time() - uij_t0)

            uij_t0 = tm.time()

            for x in range(3):

                # form dist_fock_uij_ov_2 from mpi_master

                if self.rank == mpi_master():
                    # transform to MO basis
                    fock_uij_ov_2 = 2.0 * np.linalg.multi_dot(
                        [mo_occ.T, fock_uij[x], mo_vir])
                    fock_uij_ov_2 = fock_uij_ov_2.reshape(nocc * nvir)
                else:
                    fock_uij_ov_2 = None

                dist_fock_uij_ov_2 = DistributedArray(fock_uij_ov_2, self.comm)

                # form dist_fock_deriv_ov_ix and dist_orben_ovlp_deriv_ov_ix
                # from root_rank

                key_ix = (iatom, x)
                ovlp_deriv_ao_ix = dist_ovlp_deriv_ao[key_ix].get_full_matrix(
                    root=root_rank)
                fock_deriv_ao_ix = dist_fock_deriv_ao[key_ix].get_full_matrix(
                    root=root_rank)

                # TODO: do this on mpi_master
                if self.rank == root_rank:

                    fock_deriv_ov_dict_ix = np.linalg.multi_dot(
                        [mo_occ.T, fock_deriv_ao_ix,
                         mo_vir]).reshape(nocc * nvir)

                    ovlp_deriv_ov_ix = np.linalg.multi_dot(
                        [mo_occ.T, ovlp_deriv_ao_ix, mo_vir])

                    orben_ovlp_deriv_ov_ix = (eocc.reshape(-1, 1) *
                                              ovlp_deriv_ov_ix).reshape(nocc *
                                                                        nvir)

                else:
                    fock_deriv_ov_dict_ix = None
                    orben_ovlp_deriv_ov_ix = None

                dist_fock_deriv_ov_ix = DistributedArray(fock_deriv_ov_dict_ix,
                                                         self.comm,
                                                         root=root_rank)

                dist_orben_ovlp_deriv_ov_ix = DistributedArray(
                    orben_ovlp_deriv_ov_ix, self.comm, root=root_rank)

                dist_cphf_rhs_ix_data = (dist_fock_uij_ov_2.data +
                                         dist_fock_deriv_ov_ix.data -
                                         dist_orben_ovlp_deriv_ov_ix.data)

                dist_cphf_rhs_ix = DistributedArray(dist_cphf_rhs_ix_data,
                                                    self.comm,
                                                    distribute=False)

                dist_cphf_rhs.append(dist_cphf_rhs_ix)

            profiler.add_timing_info('distRHS', tm.time() - uij_t0)

        hessian_eri_overlap = self.comm.reduce(hessian_eri_overlap,
                                               root=mpi_master())

        t2 = tm.time()

        if self.rank == mpi_master():
            self.ostream.print_info('CPHF/CPKS right-hand side computed' +
                                    ' in {:.2f} sec.'.format(t2 - t1))
            self.ostream.print_blank()
            self.ostream.flush()

        profiler.print_timing(self.ostream)
        profiler.print_profiling_summary(self.ostream)

        return {
            'dist_cphf_rhs': dist_cphf_rhs,
            # 'ovlp_deriv_oo': ovlp_deriv_oo,
            # 'fock_deriv_ao': fock_deriv_ao,
            # 'fock_uij': fock_uij,
            'hessian_first_integral_derivatives': hessian_first_integral_derivatives,
            'hessian_eri_overlap': hessian_eri_overlap,
        }

    def _compute_fmat_deriv(self,
                            molecule,
                            basis,
                            density,
                            i,
                            eri_dict,
                            point_charges=None,
                            qm_vdw_params=None):
        """
        Computes the derivative of the Fock matrix with respect
        to the coordinates of atom i.

        :param molecule:
            The molecule.
        :param basis:
            The basis set.
        :param density:
            The density matrix in AO basis.
        :param i:
            The atom index.
        :param eri_dict:
            The dictionary containing ERI information.

        :return fmat_deriv:
            The derivative of the Fock matrix wrt. atom i.
        """

        # number of atomic orbitals
        nao = basis.get_dimensions_of_basis()

        # initialize fmat variable
        fmat_deriv = np.zeros((3, nao, nao))

        # kinetic integral gadient
        kin_grad_drv = KineticEnergyGeom100Driver()

        gmats_kin = kin_grad_drv.compute(molecule, basis, i)

        for x, label in enumerate(['X', 'Y', 'Z']):
            gmat_kin = gmats_kin.matrix_to_numpy(label)
            fmat_deriv[x] += gmat_kin + gmat_kin.T

        gmats_kin = Matrices()

        # nuclear potential integral gradients
        npot_grad_100_drv = NuclearPotentialGeom100Driver()
        npot_grad_010_drv = NuclearPotentialGeom010Driver()

        gmats_npot_100 = npot_grad_100_drv.compute(molecule, basis, i)
        gmats_npot_010 = npot_grad_010_drv.compute(molecule, basis, i)

        for x, label in enumerate(['X', 'Y', 'Z']):
            gmat_npot_100 = gmats_npot_100.matrix_to_numpy(label)
            gmat_npot_010 = gmats_npot_010.matrix_to_numpy(label)
            fmat_deriv[x] -= gmat_npot_100 + gmat_npot_100.T + gmat_npot_010

        gmats_npot_100 = Matrices()
        gmats_npot_010 = Matrices()

        # point charges contribution
        if point_charges is not None:
            npoints = point_charges.shape[1]

            mm_coords = []
            mm_charges = []
            for p in range(npoints):
                xyz_p = point_charges[:3, p]
                chg_p = point_charges[3, p]
                mm_coords.append(xyz_p.copy())
                mm_charges.append(chg_p)

            gmats_100 = npot_grad_100_drv.compute(molecule, basis, i, mm_coords,
                                                  mm_charges)

            for x, label in enumerate(['X', 'Y', 'Z']):
                gmat_100 = gmats_100.matrix_to_numpy(label)
                fmat_deriv[x] -= gmat_100 + gmat_100.T

            gmats_100 = Matrices()

        if self._embedding_hess_drv is not None:
            pe_fock_grad_contr = self._embedding_hess_drv.compute_pe_fock_gradient_contributions(i=i)
            for x in range(3):
                fmat_deriv[x] += pe_fock_grad_contr[x]

        if self._dft:
            if self.xcfun.is_hybrid():
                fock_type = '2jkx'
                exchange_scaling_factor = self.xcfun.get_frac_exact_exchange()
            else:
                fock_type = 'j'
                exchange_scaling_factor = 0.0
        else:
            exchange_scaling_factor = 1.0
            fock_type = "2jk"

        # TODO: range-separated Fock
        need_omega = (self._dft and self.xcfun.is_range_separated())
        if need_omega:
            assert_msg_critical(
                False, 'HessianOrbitalResponse: Not implemented for' +
                ' range-separated functional')

        den_mat_for_fock = make_matrix(basis, mat_t.symmetric)
        den_mat_for_fock.set_values(density)

        # ERI threshold
        thresh_int = int(-math.log10(self.eri_thresh))

        # screening
        screener = eri_dict['screening']

        screener_atom = T4CScreener()
        screener_atom.partition_atom(basis, molecule, 'eri', i)

        # Fock gradient
        fock_grad_drv = FockGeom1000Driver()
        gmats_eri = fock_grad_drv.compute(basis, screener_atom, screener,
                                          den_mat_for_fock, i, fock_type,
                                          exchange_scaling_factor, 0.0,
                                          thresh_int)

        # scaling of Fock gradient for non-hybrid functionals
        factor = 2.0 if fock_type == 'j' else 1.0

        # calculate gradient contributions
        for x, label in enumerate(['X', 'Y', 'Z']):
            gmat_eri = gmats_eri.matrix_to_numpy(label)
            fmat_deriv[x] += gmat_eri * factor

        gmats_eri = Matrices()

        return fmat_deriv

    def print_cphf_header(self, title):
        """
        Prints information on the solver setup
        """

        self.ostream.print_blank()
        self.ostream.print_header('{:s} Setup'.format(title))
        self.ostream.print_header('=' * (len(title) + 8))
        self.ostream.print_blank()

        str_width = 70

        # print general info
        cur_str = 'Solver Type                     : '
        if self.use_subspace_solver:
            cur_str += 'Iterative Subspace Algorithm'
        else:
            cur_str += 'Conjugate Gradient'
        self.ostream.print_header(cur_str.ljust(str_width))

        cur_str = 'Max. Number of Iterations       : ' + str(self.max_iter)
        self.ostream.print_header(cur_str.ljust(str_width))
        cur_str = 'Convergence Threshold           : {:.1e}'.format(
            self.conv_thresh)
        self.ostream.print_header(cur_str.ljust(str_width))

        if self._dft:
            cur_str = 'Exchange-Correlation Functional : '
            cur_str += self.xcfun.get_func_label().upper()
            self.ostream.print_header(cur_str.ljust(str_width))
            grid_level = (get_default_grid_level(self.xcfun)
                          if self.grid_level is None else self.grid_level)
            cur_str = 'Molecular Grid Level            : ' + str(grid_level)
            self.ostream.print_header(cur_str.ljust(str_width))

        self.ostream.print_blank()
        self.ostream.flush()
