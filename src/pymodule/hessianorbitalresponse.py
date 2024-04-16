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

# TODO: remove commented out code;
from mpi4py import MPI
import numpy as np
import time as tm
import sys

from .veloxchemlib import mpi_master
from .veloxchemlib import denmat
from .veloxchemlib import fockmat
from .veloxchemlib import AODensityMatrix
from .veloxchemlib import AOFockMatrix
from .veloxchemlib import ElectronRepulsionIntegralsDriver
from .veloxchemlib import XCMolecularHessian
from .outputstream import OutputStream
from .profiler import Profiler
from .distributedarray import DistributedArray
from .cphfsolver import CphfSolver
from .errorhandler import assert_msg_critical
from .inputparser import parse_input
from .qqscheme import get_qq_scheme
from .batchsize import get_batch_size
from .batchsize import get_number_of_batches
from .dftutils import get_default_grid_level
from scipy.sparse import linalg

# For PySCF integral derivatives
from .import_from_pyscf import overlap_deriv
from .import_from_pyscf import fock_deriv
from .import_from_pyscf import vxc_deriv
from .import_from_pyscf import eri_deriv

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

        self.profiler = Profiler({
            'timing': self.timing,
            'profiling': self.profiling,
            'memory_profiling': self.memory_profiling,
            'memory_tracing': self.memory_tracing,
        })

        # DFT information
        dft_dict = self._init_dft(molecule, scf_tensors)
        # ERI information
        eri_dict = self._init_eri(molecule, basis)
        # PE information
        pe_dict = self._init_pe(molecule, basis)

        self.profiler.start_timer('RHS')
        natm = molecule.number_of_atoms()

        if self.rank == mpi_master():

            density = scf_tensors['D_alpha']
            mo = scf_tensors['C_alpha']
            mo_energies = scf_tensors['E']
            nao = mo.shape[0]

            # nmo is sometimes different than nao (because of linear
            # dependencies which get removed during SCF)
            nmo = mo_energies.shape[0]
            nocc = molecule.number_of_alpha_electrons()
            nvir = nmo - nocc
            mo_occ = mo[:, :nocc]
            mo_vir = mo[:, nocc:]
            eocc = mo_energies[:nocc]
            eoo = eocc.reshape(-1, 1) + eocc # ei + ej
            omega_ao = - np.linalg.multi_dot([mo_occ, np.diag(eocc), mo_occ.T])
            evir = mo_energies[nocc:]
            eov = eocc.reshape(-1, 1) - evir

            # preparing the CPHF RHS
            ovlp_deriv_ao = np.zeros((natm, 3, nao, nao))
            fock_deriv_ao = np.zeros((natm, 3, nao, nao))
        else:
            density = None

        density = self.comm.bcast(density, root=mpi_master())
        mol_grid = dft_dict['molgrid'] 
        gs_density = dft_dict['gs_density']

        # import the integral derivatives
        self.profiler.set_timing_key('derivs')
        self.profiler.start_timer('derivs')

        t0 = tm.time()
        # TODO: test if XCMolecularHessian obejct outside 
        # for-loop gives the correct Hessian
        if scf_drv._dft: 
            xc_mol_hess = XCMolecularHessian()
        for i in range(natm):
            if scf_drv._dft:
                vxc_deriv_i = xc_mol_hess.integrate_vxc_fock_gradient(
                                    molecule, basis, gs_density, mol_grid,
                                    scf_drv.xcfun.get_func_label(), i)
                vxc_deriv_i = self.comm.reduce(vxc_deriv_i, root=mpi_master())
            if self.rank == mpi_master():
                ovlp_deriv_ao[i] = overlap_deriv(molecule, basis, i)
                fock_deriv_ao[i] = fock_deriv(molecule, basis, density,
                                                i, scf_drv)
                if scf_drv._dft:
                    fock_deriv_ao[i] += vxc_deriv_i
        t1 = tm.time()

        if self.rank == mpi_master():
            self.ostream.print_info("CPHF/CPKS import of integral derivatives"
                                    + ' took'
                                    + ' {:.2f} sec.'.format(t1 - t0))
            self.ostream.print_blank()
            self.ostream.flush()
            
        self.profiler.stop_timer('derivs')

        if self.rank == mpi_master():
            # transform integral derivatives to MO basis
            ovlp_deriv_ov = np.zeros((natm, 3, nocc, nvir))
            ovlp_deriv_oo = np.zeros((natm, 3, nocc, nocc))
            fock_deriv_ov = np.zeros((natm, 3, nocc, nvir))
            orben_ovlp_deriv_ov = np.zeros((natm, 3, nocc, nvir))
            tmp_eocc = eocc.reshape(-1,1)

            for a in range(natm):
                for x in range(3):
                    ovlp_deriv_ov[a,x] = np.linalg.multi_dot([
                        mo_occ.T, ovlp_deriv_ao[a,x], mo_vir
                    ])
                    ovlp_deriv_oo[a,x] = np.linalg.multi_dot([
                        mo_occ.T, ovlp_deriv_ao[a,x], mo_occ
                    ])
                    fock_deriv_ov[a,x] = np.linalg.multi_dot([
                        mo_occ.T, fock_deriv_ao[a,x], mo_vir
                    ])
                    orben_ovlp_deriv_ov[a,x] = np.multiply(
                        tmp_eocc, ovlp_deriv_ov[a,x]
                    ) 

            # the oo part of the CPHF coefficients in AO basis,
            # transforming the oo overlap derivative back to AO basis
            # (not equal to the initial one)
            uij_ao = np.zeros((natm, 3, nao, nao))
            for a in range(natm):
                for x in range(3):
                    uij_ao[a,x] = np.linalg.multi_dot([
                       mo_occ, -0.5 * ovlp_deriv_oo[a,x], mo_occ.T 
                    ])
            uij_ao = uij_ao.reshape(3*natm, nao, nao)
                       
            uij_ao_list = list([uij_ao[x] for x in range(natm * 3)])

            # create AODensity and Fock matrix objects, contract with ERI
            ao_density_uij = AODensityMatrix(uij_ao_list, denmat.rest)
        else:
            ao_density_uij = AODensityMatrix()

        ao_density_uij.broadcast(self.rank, self.comm)

        fock_uij = AOFockMatrix(ao_density_uij)

        self._comp_lr_fock(fock_uij, ao_density_uij, molecule, basis,
                           eri_dict, dft_dict, pe_dict, self.profiler)
       
        t2 = tm.time() 
        if self.rank == mpi_master():
            self.ostream.print_info('CPHF/CPKS RHS computed in' +
                                     ' {:.2f} sec.'.format(t2 - t1))
            self.ostream.print_blank()
            self.ostream.flush()

            # TODO: how can this be done better?
            fock_uij_numpy = np.zeros((natm, 3, nao, nao))
            fock_uij_mo = np.zeros((natm, 3, nocc, nvir))
            for i in range(natm):
                for x in range(3):
                    fock_uij_numpy[i,x] = fock_uij.to_numpy(3*i + x)
                    # transform to MO basis
                    fock_uij_mo[i,x] = np.linalg.multi_dot([
                        mo_occ.T, fock_uij_numpy[i, x], mo_vir
                        ])
                    
            # sum up the terms of the RHS
            cphf_rhs = (fock_deriv_ov - orben_ovlp_deriv_ov
                        + 2 * fock_uij_mo)

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
