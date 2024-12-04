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
import numpy as np
import time as tm
import sys
import math

from .veloxchemlib import mpi_master
from .veloxchemlib import denmat
from .veloxchemlib import AODensityMatrix
#from .veloxchemlib import ElectronRepulsionIntegralsDriver
from .veloxchemlib import XCMolecularHessian
from .veloxchemlib import make_matrix, mat_t, partition_atoms
from .veloxchemlib import T4CScreener
from .outputstream import OutputStream
from .matrices import Matrices
from .distributedarray import DistributedArray
from .cphfsolver import CphfSolver
from .errorhandler import assert_msg_critical
from .inputparser import parse_input
from .batchsize import get_batch_size
from .batchsize import get_number_of_batches
from .dftutils import get_default_grid_level
from scipy.sparse import linalg

# For PySCF integral derivatives
from .import_from_pyscf import overlap_deriv
from .import_from_pyscf import fock_deriv
from .import_from_pyscf import vxc_deriv
from .import_from_pyscf import eri_deriv
from .import_from_pyscf import hcore_deriv

# import veloxchem integrals
from .veloxchemlib import (OverlapGeom100Driver, KineticEnergyGeom100Driver,
NuclearPotentialGeom100Driver,
NuclearPotentialGeom010Driver, FockGeom1000Driver)

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

    def compute_rhs(self, molecule, basis, scf_tensors, scf_drv):
        """
        Computes the right hand side for the CPHF equations for
        the analytical Hessian, all atomic coordinates.

        :param molecule:
            The molecule.
        :param basis:
            The AO basis set.
        :param scf_tensors:
            The tensors from the converged SCF calculation.
        :param scf_drv:
            The scf_driver

        :returns:
            The RHS of the CPHF equations.
        """

        # DFT information
        dft_dict = self._init_dft(molecule, scf_tensors)
        # ERI information
        eri_dict = self._init_eri(molecule, basis)
        # PE information
        pe_dict = self._init_pe(molecule, basis)

        self.profiler.start_timer('RHS')

        # number of atoms
        natm = molecule.number_of_atoms()
        nao = basis.get_dimensions_of_basis()
        nocc = molecule.number_of_alpha_electrons()

        if self.rank == mpi_master():
            density = scf_tensors['D_alpha']
            mo = scf_tensors['C_alpha']
            mo_energies = scf_tensors['E_alpha']
            eocc = mo_energies[:nocc]
            tmp_eocc = eocc.reshape(-1,1)
            # DFT
            gs_density = dft_dict['gs_density']
        else:
            density = None
            mol_grid = None
            gs_density = None
            tmp_eocc = None
            mo = None

        # DFT grid
        mol_grid = dft_dict['molgrid'] 

        density = self.comm.bcast(density, root=mpi_master())
        gs_density = self.comm.bcast(gs_density, root=mpi_master())
        tmp_eocc = self.comm.bcast(tmp_eocc, root=mpi_master())
        mo = self.comm.bcast(mo, root=mpi_master())

        # MO coefficients
        mo_occ = mo[:, :nocc]
        mo_vir = mo[:, nocc:]
        nvir = mo_vir.shape[1]

        # partition atoms for parallellisation
        local_atoms = partition_atoms(natm, self.rank, self.nodes)

        # preparing the CPHF RHS
        ovlp_deriv_ao = np.zeros((natm, 3, nao, nao))
        fock_deriv_ao = np.zeros((natm, 3, nao, nao))

        # import the integral derivatives
        # For some reason the commented-out line results in a
        # profiler error which I don't understand. Not sure why.
        # It works if it's commented out. TODO: remove?
        #self.profiler.set_timing_key('derivs')
        self.profiler.start_timer('derivs')

        t0 = tm.time()

        if scf_drv._dft: 
            xc_mol_hess = XCMolecularHessian()

        # overlap gradient driver:
        ovlp_grad_drv = OverlapGeom100Driver()

        for iatom in local_atoms:

            # compute overlap gradient integrals matrices
            gmats = ovlp_grad_drv.compute(molecule, basis, iatom)

            for x, label in enumerate(['X', 'Y', 'Z']):
                gmat = gmats.matrix_to_numpy(label)
                ovlp_deriv_ao[iatom, x] = gmat + gmat.T

            gmat = Matrices()

            # compute full Fock integrals  matrix
            fock_deriv_ao[iatom] = self._compute_fmat_deriv(molecule, basis,
                                                        scf_drv, density, iatom)

            if scf_drv._dft:
                vxc_deriv_i = xc_mol_hess.integrate_vxc_fock_gradient(
                                    molecule, basis, gs_density, mol_grid,
                                    scf_drv.xcfun.get_func_label(), iatom)
                vxc_deriv_i = self.comm.reduce(vxc_deriv_i, root=mpi_master())

                fock_deriv_ao[iatom] += vxc_deriv_i

        ovlp_deriv_ao = self.comm.allreduce(ovlp_deriv_ao, op=MPI.SUM)
        fock_deriv_ao = self.comm.allreduce(fock_deriv_ao, op=MPI.SUM)

        t1 = tm.time()

        if self.rank == mpi_master():
            self.ostream.print_info("CPHF/CPKS import of integral derivatives"
                                    + ' took'
                                    + ' {:.2f} sec.'.format(t1 - t0))
            self.ostream.print_blank()
            self.ostream.flush()
            
        self.profiler.stop_timer('derivs')

        ovlp_deriv_ov = np.zeros((natm, 3, nocc, nvir))
        ovlp_deriv_oo = np.zeros((natm, 3, nocc, nocc))
        fock_deriv_ov = np.zeros((natm, 3, nocc, nvir))
        orben_ovlp_deriv_ov = np.zeros((natm, 3, nocc, nvir))

        # transform integral derivatives to MO basis
        for iatom in local_atoms:
            for x in range(3):
                ovlp_deriv_ov[iatom,x] = np.linalg.multi_dot([
                    mo_occ.T, ovlp_deriv_ao[iatom,x], mo_vir
                ])
                ovlp_deriv_oo[iatom,x] = np.linalg.multi_dot([
                    mo_occ.T, ovlp_deriv_ao[iatom,x], mo_occ
                ])
                fock_deriv_ov[iatom,x] = np.linalg.multi_dot([
                    mo_occ.T, fock_deriv_ao[iatom,x], mo_vir
                ])
                orben_ovlp_deriv_ov[iatom,x] = np.multiply(
                    tmp_eocc, ovlp_deriv_ov[iatom,x]
                ) 

        # the oo part of the CPHF coefficients in AO basis,
        # transforming the oo overlap derivative back to AO basis
        # (not equal to the initial one)
        uij_ao = np.zeros((natm, 3, nao, nao))
        
        for iatom in local_atoms:
            for x in range(3):
                uij_ao[iatom,x] = np.linalg.multi_dot([
                   mo_occ, -0.5 * ovlp_deriv_oo[iatom,x], mo_occ.T 
                ])

        uij_ao = uij_ao.reshape(3*natm, nao, nao)
        uij_ao = self.comm.reduce(uij_ao, root=mpi_master())
                   
        if self.rank == mpi_master():
            num_uij = len(uij_ao)
        else:
            num_uij = None
        num_uij = self.comm.bcast(num_uij, root=mpi_master())

        if self.rank == mpi_master():
            uij_ao_list = list([uij_ao[x] for x in range(natm * 3)])
        else:
            uij_ao_list = [None for idx in range(num_uij)]

        for idx in range(num_uij):
            uij_ao_list[idx] = self.comm.bcast(uij_ao_list[idx], root=mpi_master())

        # create AODensity and Fock matrix objects, contract with ERI
        fock_uij = self._comp_lr_fock(uij_ao_list, molecule, basis,
                           eri_dict, dft_dict, pe_dict, self.profiler)

        # _comp_lr_fock only returns value to master
        fock_uij = self.comm.bcast(fock_uij, root=mpi_master())
       
        fock_uij_numpy = np.zeros((natm, 3, nao, nao))
        fock_uij_mo = np.zeros((natm, 3, nocc, nvir))

        for iatom in local_atoms:
            for x in range(3):
                tmp_uij_ao = fock_uij[3*iatom + x]
                # transform to MO basis
                fock_uij_mo[iatom,x] = np.linalg.multi_dot([
                    mo_occ.T, tmp_uij_ao, mo_vir
                    ])
                fock_uij_numpy[iatom,x] = tmp_uij_ao
                
        # sum up the terms of the RHS
        cphf_rhs = (fock_deriv_ov - orben_ovlp_deriv_ov
                        + 2 * fock_uij_mo)

        fock_uij_numpy = self.comm.reduce(fock_uij_numpy, root=mpi_master())
        cphf_rhs = self.comm.reduce(cphf_rhs, root=mpi_master())

        t2 = tm.time() 

        if self.rank == mpi_master():
            self.ostream.print_info('CPHF/CPKS RHS computed in' +
                                     ' {:.2f} sec.'.format(t2 - t1))
            self.ostream.print_blank()
            self.ostream.flush()

        self.profiler.stop_timer('RHS')

        if self.rank == mpi_master():
            return {
                'cphf_rhs': cphf_rhs.reshape(natm*3, nocc, nvir),
                'ovlp_deriv_oo': ovlp_deriv_oo,
                'fock_deriv_ao': fock_deriv_ao,
                'fock_uij': fock_uij_numpy,
            }
        else:
            return {}

    def _compute_fmat_deriv(self, molecule, basis, scf_drv, density, i):
        """
        Computes the derivative of the Fock matrix with respect
        to the coordinates of atom i.

        :param molecule:
            The molecule.
        :param basis:
            The basis set.
        :param density:
            The density matrix in AO basis.
        :param scf_drv:
            The SCF Driver.
        :param i:
            The atom index.

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

        # nuclear potential integral gradients
        npot_grad_100_drv = NuclearPotentialGeom100Driver()
        npot_grad_010_drv = NuclearPotentialGeom010Driver()
        gmats_npot_100 = npot_grad_100_drv.compute(molecule, basis, i)
        gmats_npot_010 = npot_grad_010_drv.compute(molecule, basis, i)

        # for Fock matrix gradient
        exchange_scaling_factor = 1.0
        fock_type = "2jk"

        if scf_drv._dft:
            # TODO: range-separated Fock
            if scf_drv.xcfun.is_hybrid():
                fock_type = '2jkx'
                exchange_scaling_factor = scf_drv.xcfun.get_frac_exact_exchange()
            else:
                fock_type = 'j'
                exchange_scaling_factor = 0.0
        
        den_mat_for_fock = make_matrix(basis, mat_t.symmetric)
        den_mat_for_fock.set_values(density)

        # ERI threshold
        thresh_int = int(-math.log10(scf_drv.eri_thresh))

        # screening
        screener = T4CScreener()
        screener.partition(basis, molecule, 'eri')
        screener_atom = T4CScreener()
        screener_atom.partition_atom(basis, molecule, 'eri', i)

        # Fock gradient
        fock_grad_drv = FockGeom1000Driver()
        gmats_eri = fock_grad_drv.compute(basis, screener_atom, screener,
                                      den_mat_for_fock, i, fock_type,
                                      exchange_scaling_factor, 0.0,
                                      thresh_int)

        # scaling of ERI gradient for non-hybrid functionals
        factor = 2.0 if fock_type == 'j' else 1.0

        # calculate gradient contributions
        for x, label in enumerate(['X', 'Y', 'Z']):
            gmat_kin = gmats_kin.matrix_to_numpy(label)
            gmat_npot_100 = gmats_npot_100.matrix_to_numpy(label)
            gmat_npot_010 = gmats_npot_010.matrix_to_numpy(label)
            gmat_eri = gmats_eri.matrix_to_numpy(label)

            fmat_deriv[x] += gmat_kin + gmat_kin.T
            fmat_deriv[x] -= gmat_npot_100 + gmat_npot_100.T + gmat_npot_010
            fmat_deriv[x] += gmat_eri * factor

        gmats_kin = Matrices()
        gmats_npot_100 = Matrices()
        gmats_npot_010 = Matrices()
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
