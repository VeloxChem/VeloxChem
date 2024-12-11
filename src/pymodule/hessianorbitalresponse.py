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

# TODO: remove commented out code;
from mpi4py import MPI
from pathlib import Path
import numpy as np
import time as tm
import sys
import math

from .veloxchemlib import AODensityMatrix
from .veloxchemlib import T4CScreener
from .veloxchemlib import XCMolecularHessian
from .veloxchemlib import mpi_master, denmat
from .veloxchemlib import make_matrix, mat_t, partition_atoms
from .veloxchemlib import (OverlapGeom100Driver, KineticEnergyGeom100Driver,
                           NuclearPotentialGeom100Driver,
                           NuclearPotentialGeom010Driver, FockGeom1000Driver)
from .outputstream import OutputStream
from .matrices import Matrices
from .distributedarray import DistributedArray
from .cphfsolver import CphfSolver
from .errorhandler import assert_msg_critical
from .inputparser import parse_input
from .batchsize import get_batch_size
from .batchsize import get_number_of_batches
from .dftutils import get_default_grid_level


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

    def update_settings(self, cphf_dict, method_dict=None):
        """
        Updates response and method settings in CPHF solver.

        :param cphf_dict:
            The dictionary of CPHF (orbital response) settings.
        :param method_dict:
            The dictionary of method settings.
        """

        super().update_settings(cphf_dict, method_dict)

    def compute_rhs(self, molecule, basis, scf_tensors, eri_dict, dft_dict, pe_dict):
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

        self.profiler.start_timer('RHS')

        natm = molecule.number_of_atoms()
        nocc = molecule.number_of_alpha_electrons()
        nao = basis.get_dimensions_of_basis()

        if self.rank == mpi_master():
            density = scf_tensors['D_alpha']
            eocc = scf_tensors['E_alpha'][:nocc]
            mo = scf_tensors['C_alpha']
        else:
            density = None
            eocc = None
            mo = None

        density = self.comm.bcast(density, root=mpi_master())
        eocc = self.comm.bcast(eocc, root=mpi_master())
        mo = self.comm.bcast(mo, root=mpi_master())

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
        all_atom_idx_rank = self.comm.bcast(all_atom_idx_rank, root=mpi_master())

        # preparing the CPHF RHS

        ovlp_deriv_ao_dict = {(iatom, x): None
                              for iatom in range(natm)
                              for x in range(3)}

        fock_deriv_ao_dict = {(iatom, x): None
                              for iatom in range(natm)
                              for x in range(3)}

        # TODO: double check the use of profiler
        #self.profiler.set_timing_key('derivs')
        self.profiler.start_timer('derivs')

        t0 = tm.time()

        ovlp_grad_drv = OverlapGeom100Driver()

        for iatom in local_atoms:

            # compute overlap gradient integrals matrices
            gmats = ovlp_grad_drv.compute(molecule, basis, iatom)

            for x, label in enumerate(['X', 'Y', 'Z']):
                gmat = gmats.matrix_to_numpy(label)
                ovlp_deriv_ao_dict[(iatom, x)] = gmat + gmat.T

            gmats = Matrices()

            # compute full Fock integrals matrix
            fock_deriv_ao_i = self._compute_fmat_deriv(
                molecule, basis, density, iatom, eri_dict)

            fock_deriv_ao_dict[(iatom, 0)] = fock_deriv_ao_i[0]
            fock_deriv_ao_dict[(iatom, 1)] = fock_deriv_ao_i[1]
            fock_deriv_ao_dict[(iatom, 2)] = fock_deriv_ao_i[2]

        # Note: Parallelization of DFT integration is done over grid points
        # instead of atoms.

        if self._dft:
            xc_mol_hess = XCMolecularHessian()

            for iatom, root_rank in all_atom_idx_rank:
                vxc_deriv_i = xc_mol_hess.integrate_vxc_fock_gradient(
                    molecule, basis, gs_density, mol_grid,
                    self.xcfun.get_func_label(), iatom)

                for x in range(3):
                    vxc_deriv_ix = self.comm.reduce(vxc_deriv_i[x],
                                                    root=root_rank)
                    if self.rank == root_rank:
                        fock_deriv_ao_dict[(iatom, x)] += vxc_deriv_ix

        t1 = tm.time()

        if self.rank == mpi_master():
            self.ostream.print_info("CPHF/CPKS integral derivatives computed"
                                    + ' in {:.2f} sec.'.format(t1 - t0))
            self.ostream.print_blank()
            self.ostream.flush()
            
        self.profiler.stop_timer('derivs')

        if self.rank == mpi_master():
            hessian_first_integral_derivatives = np.zeros((natm, natm, 3, 3))
        else:
            hessian_first_integral_derivatives = None

        for iatom, root_rank_i in all_atom_idx_rank:
            for x in range(3):

                ovlp_deriv_ao_ix = self.comm.bcast(
                        ovlp_deriv_ao_dict[(iatom,x)], root=root_rank_i)

                fock_deriv_ao_ix = self.comm.bcast(
                        fock_deriv_ao_dict[(iatom,x)], root=root_rank_i)

                for jatom, root_rank_j in all_atom_idx_rank:
                    if jatom < iatom:
                        continue
                    for y in range(3):

                        ovlp_deriv_ao_jy = self.comm.bcast(
                                ovlp_deriv_ao_dict[(jatom,y)], root=root_rank_j)

                        fock_deriv_ao_jy = self.comm.bcast(
                                fock_deriv_ao_dict[(jatom,y)], root=root_rank_j)

                        if self.rank == mpi_master():

                            Fix_Sjy = np.sum(
                                np.matmul(density, fock_deriv_ao_ix) *
                                np.matmul(density, ovlp_deriv_ao_jy).T)

                            Fjy_Six = np.sum(
                                np.matmul(density, fock_deriv_ao_jy) *
                                np.matmul(density, ovlp_deriv_ao_ix).T)

                            Six_Sjy = 2.0 * np.sum(
                                np.matmul(omega_ao, ovlp_deriv_ao_ix) *
                                np.matmul(density, ovlp_deriv_ao_jy).T)

                            hess_ijxy = -2.0 * (Fix_Sjy + Fjy_Six + Six_Sjy)
                            hessian_first_integral_derivatives[iatom,jatom,x,y] += hess_ijxy
                            if iatom != jatom:
                                hessian_first_integral_derivatives[jatom,iatom,y,x] += hess_ijxy

        dist_cphf_rhs = []

        if self.rank == mpi_master():
            hessian_eri_overlap = np.zeros((natm, natm, 3, 3))
        else:
            hessian_eri_overlap = None

        for iatom, root_rank in all_atom_idx_rank:

            # the oo part of the CPHF coefficients in AO basis
            uij_ao_list = []
            for x in range(3):
                ovlp_deriv_ao_ix = self.comm.bcast(
                        ovlp_deriv_ao_dict[(iatom,x)], root=root_rank)
                uij_ao_list.append(np.linalg.multi_dot([
                   density, -0.5 * ovlp_deriv_ao_ix, density
                ]))

            # create AODensity and Fock matrix objects, contract with ERI
            fock_uij = self._comp_lr_fock(uij_ao_list, molecule, basis, eri_dict,
                                          dft_dict, pe_dict, self.profiler)

            # Note: hessian_eri_overlap, i.e. inner product of P_P_Six_fock_ao and
            # P_P_Sjy, is obtained from fock_uij and uij_ao_jy in hessian orbital
            # response

            for jatom, root_rank_j in all_atom_idx_rank:
                if jatom < iatom:
                    continue

                hess_ij = np.zeros((3, 3))

                for y in range(3):
                    ovlp_deriv_ao_jy = self.comm.bcast(
                            ovlp_deriv_ao_dict[(jatom,y)], root=root_rank_j)
                    uij_ao_jy = np.linalg.multi_dot([
                       density, -0.5 * ovlp_deriv_ao_jy, density
                    ])

                    if self.rank == mpi_master():
                        for x in range(3):
                            # xmn,ymn->xy
                            hess_ij[x, y] = 2.0 * (np.sum(
                                fock_uij[x] * uij_ao_jy))

                if self.rank == mpi_master():
                    hessian_eri_overlap[iatom,jatom,:,:] += 4.0 * hess_ij
                    if iatom != jatom:
                        hessian_eri_overlap[jatom,iatom,:,:] += 4.0 * hess_ij.T

            for x in range(3):

                # form dist_fock_uij_ov_2 from mpi_master

                if self.rank == mpi_master():
                    # transform to MO basis
                    fock_uij_ov_2 = 2.0 * np.linalg.multi_dot([
                        mo_occ.T, fock_uij[x], mo_vir])
                    fock_uij_ov_2 = fock_uij_ov_2.reshape(nocc * nvir)
                else:
                    fock_uij_ov_2 = None

                dist_fock_uij_ov_2 = DistributedArray(fock_uij_ov_2, self.comm)

                # form dist_fock_deriv_ov_ix and dist_orben_ovlp_deriv_ov_ix
                # from root_rank

                if self.rank == root_rank:

                    fock_deriv_ov_dict_ix = np.linalg.multi_dot([
                        mo_occ.T, fock_deriv_ao_dict[(iatom,x)], mo_vir
                    ]).reshape(nocc * nvir)

                    ovlp_deriv_ov_ix = np.linalg.multi_dot([
                        mo_occ.T, ovlp_deriv_ao_dict[(iatom,x)], mo_vir
                    ])

                    orben_ovlp_deriv_ov_ix = (
                        eocc.reshape(-1, 1) * ovlp_deriv_ov_ix
                    ).reshape(nocc * nvir)

                else:
                    fock_deriv_ov_dict_ix = None
                    orben_ovlp_deriv_ov_ix = None

                dist_fock_deriv_ov_ix = DistributedArray(
                        fock_deriv_ov_dict_ix,
                        self.comm,
                        root=root_rank)

                dist_orben_ovlp_deriv_ov_ix = DistributedArray(
                        orben_ovlp_deriv_ov_ix,
                        self.comm,
                        root=root_rank)

                dist_cphf_rhs_ix_data = (dist_fock_uij_ov_2.data +
                                         dist_fock_deriv_ov_ix.data -
                                         dist_orben_ovlp_deriv_ov_ix.data)

                dist_cphf_rhs_ix = DistributedArray(
                        dist_cphf_rhs_ix_data,
                        self.comm,
                        distribute=False)

                dist_cphf_rhs.append(dist_cphf_rhs_ix)

        ovlp_deriv_ao_dict.clear()
        fock_deriv_ao_dict.clear()

        t2 = tm.time() 

        if self.rank == mpi_master():
            self.ostream.print_info('CPHF/CPKS right-hand side computed' +
                                     ' in {:.2f} sec.'.format(t2 - t1))
            self.ostream.print_blank()
            self.ostream.flush()

        self.profiler.stop_timer('RHS')

        return {
            'dist_cphf_rhs': dist_cphf_rhs,
            # 'ovlp_deriv_oo': ovlp_deriv_oo,
            # 'fock_deriv_ao': fock_deriv_ao,
            # 'fock_uij': fock_uij,
            'hessian_first_integral_derivatives': hessian_first_integral_derivatives,
            'hessian_eri_overlap': hessian_eri_overlap,
        }

    def _compute_fmat_deriv(self, molecule, basis, density, i, eri_dict):
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

        self.ostream.print_blank()
        self.ostream.flush()
