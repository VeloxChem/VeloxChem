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
import numpy as np
import time as tm
import sys

from .molecule import Molecule
from .gradientdriver import GradientDriver
from .hessiandriver import HessianDriver
from .scfgradientdriver import ScfGradientDriver
from .outputstream import OutputStream
from .firstorderprop import FirstOrderProperties
from .orbitalresponse import OrbitalResponse
from .lrsolver import LinearResponseSolver
from .profiler import Profiler
from .qqscheme import get_qq_scheme
from .veloxchemlib import mpi_master
from .veloxchemlib import ElectronRepulsionIntegralsDriver
from .veloxchemlib import denmat
from .veloxchemlib import fockmat
from .veloxchemlib import AODensityMatrix
from .veloxchemlib import AOFockMatrix
from .veloxchemlib import ElectricDipoleIntegralsDriver
from scipy.sparse import linalg

# For PySCF integral derivatives
from .import_from_pyscf import overlap_deriv
from .import_from_pyscf import fock_deriv
from .import_from_pyscf import eri_deriv
from .import_from_pyscf import overlap_second_deriv
from .import_from_pyscf import hcore_second_deriv
from .import_from_pyscf import eri_second_deriv


class ScfHessianDriver(HessianDriver):
    """
    Implements SCF Hessian driver.

    :param scf_drv:
        The SCF driver.
    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.

    Instance variables
        - scf_drv: The SCF driver.
        - do_raman: Additionally calculate Raman intensities
                (at significantly higher computational cost).
    """

    def __init__(self, scf_drv, comm=None, ostream=None):
        """
        Initializes SCF Hessian driver.
        """

        super().__init__(comm, ostream)

        self.flag = 'SCF Hessian Driver'
        self.scf_drv = scf_drv
        self.do_raman = False
        self.pople = False

        # Solver setup
        self.conv_thresh = 1.0e-4
        self.max_iter = 50
        self.iter_count = 0
        self.is_converged = False


    def update_settings(self, method_dict, freq_dict=None, cphf_dict=None):
        """
        Updates settings in ScfHessianDriver.

        :param method_dict:
            The input dictionary of method settings group.
        :param freq_dict:
            The input dictionary of Hessian/frequency settings group.
        :param cphf_dict:
            The input dictionary of CPHF (orbital response) settings.
        """

        super().update_settings(method_dict, freq_dict)

        if freq_dict is None:
            freq_dict = {}

        # Settings for orbital response module
        if cphf_dict is None:
            cphf_dict = {}
        self.cphf_dict = dict(cphf_dict)

        if 'conv_thresh' in cphf_dict:
            self.conv_thresh = float(cphf_dict['conv_thresh'])

        if 'max_iter' in cphf_dict:
            self.max_iter = int(cphf_dict['max_iter']) 

        # check if Raman intensities are to be calculated (numerically)
        if 'do_raman' in freq_dict:
            key = freq_dict['do_raman'].lower()
            self.do_raman = True if key in ['yes', 'y'] else False

        # The electronic energy
        self.elec_energy = self.scf_drv.get_scf_energy()

        # TODO: maybe this should be solved differently
        self.profiler = Profiler({
            'timing': self.timing,
            'profiling': self.profiling,
            'memory_profiling': self.memory_profiling,
            'memory_tracing': self.memory_tracing,
        })


    def compute(self, molecule, ao_basis, min_basis=None):
        """
        Computes the analytical or numerical nuclear Hessian.

        :param molecule:
            The molecule.
        :param ao_basis:
            The AO basis set.
        :param min_basis:
            The minimal AO basis set.
        """
        self.print_header()

        start_time = tm.time()

        if self.numerical:
            self.compute_numerical(molecule, ao_basis, min_basis)
        else:
            self.compute_analytical(molecule, ao_basis)

        # print Hessian
        if self.do_print_hessian is True:
            self.print_geometry(molecule)
            self.ostream.print_blank()
            self.print_hessian(molecule)

        valstr = '*** Time spent in Hessian calculation: '
        valstr += '{:.2f} sec ***'.format(tm.time() - start_time)
        self.ostream.print_header(valstr)
        self.ostream.print_blank()
        self.ostream.print_blank()
        self.ostream.flush()

    def compute_analytical(self, molecule, ao_basis):
        """
        Computes the analytical nuclear Hessian.
        So far only for restricted Hartree-Fock with PySCF integral derivatives...

        :param molecule:
            The molecule.
        :param ao_basis:
            The AO basis set.
        """

        density = self.scf_drv.scf_tensors['D_alpha']
        natm = molecule.number_of_atoms()
        mo = self.scf_drv.scf_tensors['C_alpha']
        nao = mo.shape[0]
        nocc = molecule.number_of_alpha_electrons()
        mo_occ = mo[:, :nocc]
        mo_vir = mo[:, nocc:]
        mo_energies = self.scf_drv.scf_tensors['E']
        eocc = mo_energies[:nocc]
        eoo = eocc.reshape(-1, 1) + eocc #ei+ej
        omega_ao = - np.linalg.multi_dot([mo_occ, np.diag(eocc), mo_occ.T])

        # Solve the CPHF equations
        cphf_solution_dict = self.compute_cphf(molecule, ao_basis)
        cphf_ov = cphf_solution_dict['cphf_ov']
        ovlp_deriv_oo = cphf_solution_dict['ovlp_deriv_oo']

        # Calculate the perturbed density matrix
        # cphf_oo = -0.5*ovlp_deriv_oo
        perturbed_density = ( - np.einsum('mj,xyij,ni->xymn',
                                        mo_occ, ovlp_deriv_oo, mo_occ)
                              + np.einsum('ma,xyia,ni->xymn',
                                        mo_vir, cphf_ov, mo_occ)
                              + np.einsum('mi,xyia,na->xymn',
                                        mo_occ, cphf_ov, mo_vir)
                            )

        # Parts related to first-order integral derivatives
        if self.pople:
            fock_uij = cphf_solution_dict['fock_uij']
            fock_deriv_ao = cphf_solution_dict['fock_deriv_ao']
            fock_deriv_oo = np.einsum('mi,xymn,nj->xyij', mo_occ, fock_deriv_ao, mo_occ)
            # (ei+ej)S^\chi_ij
            orben_ovlp_deriv_oo = np.einsum('ij,xyij->xyij', eoo, ovlp_deriv_oo)
            hessian_first_order_derivatives = self.compute_pople(molecule, ao_basis,
                                             -0.5 * ovlp_deriv_oo, cphf_ov, fock_uij,
                                              fock_deriv_oo, orben_ovlp_deriv_oo, perturbed_density)
        else:
            cphf_rhs = cphf_solution_dict['cphf_rhs']
            hessian_first_order_derivatives = self.compute_furche(molecule, ao_basis,
                                                         cphf_rhs, -0.5 * ovlp_deriv_oo, cphf_ov)

        # Parts related to second-order integral derivatives
        hessian_2nd_order_derivatives = np.zeros((natm, natm, 3, 3))
        for i in range(natm):
            # do only upper triangular matrix
            for j in range(i, natm):
                # Get integral second-order derivatives
                ovlp_2nd_deriv_ii, ovlp_2nd_deriv_ij = overlap_second_deriv(molecule, ao_basis, i, j)
                hcore_2nd_deriv_ij = hcore_second_deriv(molecule, ao_basis, i, j)
                eri_2nd_deriv_ii, eri_2nd_deriv_ij = eri_second_deriv(molecule, ao_basis, i, j)
                if i == j:
                    # Add diagonal (same atom) contributions, 2S + 2J - K
                    hessian_2nd_order_derivatives[i,i] += 2*np.einsum('mn,xymn->xy', omega_ao, ovlp_2nd_deriv_ii)
                    hessian_2nd_order_derivatives[i,i] += 2*np.einsum('mn,kl,xymnkl->xy', density, density, eri_2nd_deriv_ii)
                    hessian_2nd_order_derivatives[i,i] -= np.einsum('mk,nl,xymnkl->xy', density, density, eri_2nd_deriv_ii)
                # Add non-diagonal contributions, 2S + 2J - K + 2h
                hessian_2nd_order_derivatives[i,j] += 2*np.einsum('mn,kl,xymnkl->xy', density, density, eri_2nd_deriv_ij)
                hessian_2nd_order_derivatives[i,j] -= np.einsum('mk,nl,xymnkl->xy', density, density, eri_2nd_deriv_ij)
                hessian_2nd_order_derivatives[i,j] += 2*np.einsum('mn,xymn->xy', omega_ao, ovlp_2nd_deriv_ij)
                hessian_2nd_order_derivatives[i,j] += 2*np.einsum('mn,xymn->xy', density, hcore_2nd_deriv_ij)

            # lower triangle is transpose of the upper part
            for j in range(i):
                hessian_2nd_order_derivatives[i,j] += hessian_2nd_order_derivatives[j,i].T

        # Nuclear-nuclear repulsion contribution
        hessian_nuclear_nuclear = self.hess_nuc_contrib(molecule)

        # Sum up the terms and reshape for final Hessian
        self.hessian = ( hessian_first_order_derivatives + hessian_2nd_order_derivatives
                       + hessian_nuclear_nuclear ).transpose(0,2,1,3).reshape(3*natm, 3*natm)

        # Calculate the gradient of the dipole moment, needed for IR intensities
        self.compute_dipole_gradient(molecule, ao_basis, perturbed_density)


    def compute_pople(self, molecule, ao_basis, cphf_oo, cphf_ov, fock_uij,
                      fock_deriv_oo, orben_ovlp_deriv_oo, perturbed_density):
        """
        Computes the analytical nuclear Hessian the Pople way.
        Int. J. Quantum Chem. Quantum Chem. Symp. 13, 225-241 (1979).
        DOI: 10.1002/qua.560160825

        :param molecule:
            The molecule.
        :param ao_basis:
            The AO basis set.
        :param cphf_oo:
            The oo block of the CPHF coefficients.
        :param cphf_ov:
            The ov block of the CPHF coefficients.
        :param fock_uij:
            The auxiliary Fock matrix constructed
            using the oo block of the CPHF coefficients
            and two electron integrals.
        :param fock_deriv_oo:
            The oo block of the derivative of the Fock matrix
            with respect to nuclear coordinates
        :param orben_ovlp_deriv_oo:
            The oo block of the derivative of the overlap matrix
            with respect to nuclear coordinates, multiplied with
            orbital energies (ei+ej)S^\chi_ij
        :param perturbed_density:
            The perturbed density matrix.
        """

        natm = molecule.number_of_atoms()
        nocc = molecule.number_of_alpha_electrons()
        mo = self.scf_drv.scf_tensors['C']
        mo_occ = mo[:, :nocc].copy()
        mo_vir = mo[:, nocc:].copy()
        nao = mo.shape[0]
        density = self.scf_drv.scf_tensors['D_alpha']
        mo_energies = self.scf_drv.scf_tensors['E']
        eocc = mo_energies[:nocc]
        eo_diag = np.diag(eocc)
        epsilon_dm_ao = - np.linalg.multi_dot([mo_occ, eo_diag, mo_occ.T])

        # Construct the perturbed density matrix, and perturbed omega
        # TODO: consider if using the transpose makes the
        # computation faster; consider using cphf coefficients in AO
        # to compute the perturbed density matrix.
        ##perturbed_density = ( 2 * np.einsum('mj,xyij,ni->xymn',
        ##                                mo_occ, cphf_oo, mo_occ)
        ##                      + np.einsum('ma,xyia,ni->xymn',
        ##                                mo_vir, cphf_ov, mo_occ)
        ##                      + np.einsum('mi,xyia,na->xymn',
        ##                                mo_occ, cphf_ov, mo_vir)
        ##                    )
        orben_perturbed_density = ( np.einsum('i,mj,xyij,ni->xymn',
                                            eocc, mo_occ, cphf_oo, mo_occ)
                                  + np.einsum('i,mi,xyij,nj->xymn',
                                            eocc, mo_occ, cphf_oo, mo_occ)
                                  + np.einsum('i,ma,xyia,ni->xymn',
                                            eocc, mo_vir, cphf_ov, mo_occ)
                                  +np.einsum('i,mi,xyia,na->xymn',
                                            eocc, mo_occ, cphf_ov, mo_vir)
                                 )

        # Transform the CPHF coefficients to AO:
        uia_ao = np.einsum('mk,xykb,nb->xymn', mo_occ, cphf_ov, mo_vir).reshape((3*natm, nao, nao))

        # create AODensity and Fock matrix objects, contract with ERI
        uia_ao_list = list([uia_ao[x] for x in range(natm * 3)])
        ao_density_uia = AODensityMatrix(uia_ao_list, denmat.rest)

        fock_uia = AOFockMatrix(ao_density_uia)

        fock_flag = fockmat.rgenjk
        for i in range(natm*3):
            fock_uia.set_fock_type(fock_flag, i)

        eri_drv = ElectronRepulsionIntegralsDriver(self.comm)
        screening = eri_drv.compute(get_qq_scheme(self.scf_drv.qq_type),
                                    self.scf_drv.eri_thresh, molecule, ao_basis)
        eri_drv.compute(fock_uia, ao_density_uia, molecule, ao_basis, screening)

        # TODO: can this be done in a different way?
        fock_uia_numpy = np.zeros((natm,3,nao,nao))
        for i in range(natm):
            for x in range(3):
                fock_uia_numpy[i,x] = fock_uia.to_numpy(3*i + x)

        fock_cphf_oo = np.einsum('mi,xymn,nj->xyij', mo_occ, fock_uij, mo_occ)

        fock_cphf_ov = ( np.einsum('mi,xymn,nj->xyij', mo_occ,
                                   fock_uia_numpy, mo_occ)
                        +np.einsum('mj,xymn,ni->xyij', mo_occ,
                                   fock_uia_numpy, mo_occ)
                        )
        # Construct the derivative of the omega multipliers:
        perturbed_omega_ao = - ( orben_perturbed_density
                                + np.einsum('mi,xyij,nj->xymn', mo_occ,
                                            fock_deriv_oo, mo_occ)
                                -0.5*np.einsum('mi,xyij,nj->xymn', mo_occ,
                                                orben_ovlp_deriv_oo, mo_occ)
                                + 2*np.einsum('mi,xyij,nj->xymn', mo_occ,
                                            fock_cphf_oo, mo_occ)
                                + np.einsum('mi,xyij,nj->xymn', mo_occ,
                                            fock_cphf_ov, mo_occ)
                                )

        # First integral derivatives: partial Fock and overlap matrix derivatives
        hessian_first_integral_derivatives = np.zeros((natm, natm, 3, 3))

        for i in range(natm):
            # upper triangular part
            for j in range(i, natm):
                # First derivative of the Fock matrix
                fock_deriv_j = fock_deriv(molecule, ao_basis, density, j)
                # First derivative of overlap matrix
                ovlp_deriv_j = overlap_deriv(molecule, ao_basis, j)
                # Add the contribution of the perturbed density matrix
                hessian_first_integral_derivatives[i,j] += np.einsum('xmn,ymn->xy', 2*perturbed_density[i], fock_deriv_j)
                hessian_first_integral_derivatives[i,j] += np.einsum('xmn,ymn->xy', 2*perturbed_omega_ao[i], ovlp_deriv_j)

            # lower triangular part
            for j in range(i):
                hessian_first_integral_derivatives[i,j] += hessian_first_integral_derivatives[j,i].T

        return hessian_first_integral_derivatives


    def compute_furche(self, molecule, ao_basis, cphf_rhs, cphf_oo, cphf_ov):
        """
        Computes the analytical nuclear Hessian the Furche/Ahlrichs way.
        Chem. Phys. Lett. 362, 511–518 (2002).
        DOI: 10.1016/S0009-2614(02)01084-9

        :param molecule:
            The molecule.
        :param ao_basis:
            The AO basis set.
        :param cphf_rhs:
            The RHS of the CPHF equations in MO basis.
        :param cphf_oo:
            The oo block of the CPHF coefficients.
        :param cphf_ov:
            The ov block of the CPHF coefficients.
        """

        natm = molecule.number_of_atoms()
        nocc = molecule.number_of_alpha_electrons()
        mo = self.scf_drv.scf_tensors['C']
        mo_occ = mo[:, :nocc].copy()
        nao = mo.shape[0]
        density = self.scf_drv.scf_tensors['D_alpha']
        mo_energies = self.scf_drv.scf_tensors['E']
        eocc = mo_energies[:nocc]
        omega_ao = - np.linalg.multi_dot([mo_occ, np.diag(eocc), mo_occ.T])
        eri_drv = ElectronRepulsionIntegralsDriver(self.comm)
        screening = eri_drv.compute(get_qq_scheme(self.scf_drv.qq_type),
                                    self.scf_drv.eri_thresh, molecule, ao_basis)

        # RHS contracted with CPHF coefficients (ov)
        #hessian_cphf_coeff_rhs = np.zeros((natm, natm, 3, 3))
        #for i in range(natm):
        #    for j in range(natm):
        #        hessian_cphf_coeff_rhs[i,j] += 4*np.einsum('xia,yia->xy', cphf_ov[i], cphf_rhs[j])
        hessian_cphf_coeff_rhs = 4 * np.einsum('ixka,jyka->ijxy', cphf_ov, cphf_rhs)

        # First integral derivatives: partial Fock and overlap matrix derivatives
        hessian_first_integral_derivatives = np.zeros((natm, natm, 3, 3))
        for i in range(natm):
            fock_deriv_i = fock_deriv(molecule, ao_basis, density, i)
            ovlp_deriv_i = overlap_deriv(molecule, ao_basis, i)
            # upper triangular part
            for j in range(i, natm):
                ovlp_deriv_j = overlap_deriv(molecule, ao_basis, j)
                fock_deriv_j = fock_deriv(molecule, ao_basis, density, j)
                Fix_Sjy = np.zeros((3,3))
                Fjy_Six = np.zeros((3,3))
                Six_Sjy = np.zeros((3,3))
                for x in range(3):
                    for y in range(3):
                        Fix_Sjy[x,y] = np.trace(np.linalg.multi_dot([density, fock_deriv_i[x], density, ovlp_deriv_j[y]]))
                        Fjy_Six[x,y] = np.trace(np.linalg.multi_dot([density, fock_deriv_j[y], density, ovlp_deriv_i[x]]))
                        Six_Sjy[x,y] = ( 2*np.trace(np.linalg.multi_dot([omega_ao, ovlp_deriv_i[x], density, ovlp_deriv_j[y]])) )
                hessian_first_integral_derivatives[i,j] += -2 * (Fix_Sjy + Fjy_Six + Six_Sjy)
                ##hessian_first_integral_derivatives[i,j] += 2*(-np.einsum('xmn,ykl,mk,nl->xy', fock_deriv_i, ovlp_deriv_j, density, density)
                ##                                                           -np.einsum('ymn,xkl,mk,nl->xy', fock_deriv_j, ovlp_deriv_i, density, density)
                ##                                                           - np.einsum('xmn,ykl,mk,nl->xy', ovlp_deriv_i, ovlp_deriv_j, omega_ao, density)
                ##                                                           - np.einsum('xmn,ykl,nl,mk->xy', ovlp_deriv_i, ovlp_deriv_j, omega_ao, density)
                ##                                                  )
                ##
                ##
                ##

            # lower triangular part
            for j in range(i):
                hessian_first_integral_derivatives[i,j] += hessian_first_integral_derivatives[j,i].T

        # Overlap derivative with ERIs
        hessian_eri_overlap = np.zeros((natm, natm, 3, 3))
        for i in range(natm):
            ovlp_deriv_i = overlap_deriv(molecule, ao_basis, i)
            # upper triangular part
            for j in range(i, natm):
                ovlp_deriv_j = overlap_deriv(molecule, ao_basis, j)

                # Overlap derivative contracted with two density matrices
                ##P_P_Six = np.zeros((3, nao, nao))
                ##P_P_Sjy = np.zeros((3, nao, nao))
                ##for x in range(3):
                ##    P_P_Six[x] = (np.linalg.multi_dot([density, ovlp_deriv_i[x], density]))
                ##    P_P_Sjy[x] = (np.linalg.multi_dot([density, ovlp_deriv_j[x], density]))
                P_P_Six = np.einsum('mn,xnk,kl->xml', density, ovlp_deriv_i, density)
                P_P_Sjy = np.einsum('mn,xnk,kl->xml', density, ovlp_deriv_j, density)

                # Create a list of 2D numpy arrays to create AODensityMatrix objects
                # and calculate auxiliary Fock matrix with ERI driver
                P_P_Six_ao_list = list([P_P_Six[x] for x in range(3)])
                P_P_Six_dm_ao = AODensityMatrix(P_P_Six_ao_list, denmat.rest)
                P_P_Six_fock_ao = AOFockMatrix(P_P_Six_dm_ao)
                eri_drv.compute(P_P_Six_fock_ao, P_P_Six_dm_ao, molecule, ao_basis, screening)

                # Convert the auxiliary Fock matrices to numpy arrays for further use
                np_P_P_Six_fock = np.zeros((3,nao,nao))
                for k in range(3):
                    np_P_P_Six_fock[k] = P_P_Six_fock_ao.to_numpy(k)

                hessian_eri_overlap[i,j] += 2.0 * np.einsum('xmn,ymn->xy', np_P_P_Six_fock, P_P_Sjy)

            # lower triangular part
            for j in range(i):
                hessian_eri_overlap[i,j] += hessian_eri_overlap[j,i].T

        # return the sum of the three contributions
        return (hessian_cphf_coeff_rhs + hessian_first_integral_derivatives
                         + hessian_eri_overlap)


    def compute_cphf(self, molecule, ao_basis):
        """
        Computes the coupled-perturbed Hartree-Fock (CPHF) coefficients.

        :param molecule:
            The molecule.
        :param ao_basis:
            The AO basis set.

        :return:
            A dictionary containing the RHS and solution (ov block)
            of the CPHF equations, the derivative of the AO Fock matrix,
            partial derivative of the overlap matrix (oo block),
            and an auxiliary Fock matrix (oo block of the CPHF coefficients
            contracted with the two-electron integrals).
        """

        # TODO: remove; add profiler
        start = tm.time()
        density = self.scf_drv.scf_tensors['D_alpha']
        #overlap = self.scf_drv.scf_tensors['S']
        natm = molecule.number_of_atoms()
        mo = self.scf_drv.scf_tensors['C_alpha']
        nao = mo.shape[0]
        nocc = molecule.number_of_alpha_electrons()
        nvir = nao - nocc
        mo_occ = mo[:, :nocc]
        mo_vir = mo[:, nocc:]
        mo_energies = self.scf_drv.scf_tensors['E']
        eocc = mo_energies[:nocc]
        eoo = eocc.reshape(-1, 1) + eocc #ei+ej
        omega_ao = - np.linalg.multi_dot([mo_occ, np.diag(eocc), mo_occ.T])
        evir = mo_energies[nocc:]
        eov = eocc.reshape(-1, 1) - evir


        # preparing the CPHF RHS

        ovlp_deriv_ao = np.zeros((natm, 3, nao, nao))
        fock_deriv_ao = np.zeros((natm, 3, nao, nao))

        # import the integral derivatives
        for i in range(natm):
            ovlp_deriv_ao[i] = overlap_deriv(molecule, ao_basis, i)
            fock_deriv_ao[i] = fock_deriv(molecule, ao_basis, density, i)

        # transform integral derivatives to MO basis
        ovlp_deriv_ov = np.einsum('mi,xymn,na->xyia', mo_occ, ovlp_deriv_ao, mo_vir)
        ovlp_deriv_oo = np.einsum('mi,xymn,nj->xyij', mo_occ, ovlp_deriv_ao, mo_occ)
        fock_deriv_ov = np.einsum('mi,xymn,na->xyia', mo_occ, fock_deriv_ao, mo_vir)
        orben_ovlp_deriv_ov = np.einsum('i,xyia->xyia', eocc, ovlp_deriv_ov)

        # the oo part of the CPHF coefficients in AO basis,
        # transforming the oo overlap derivative back to AO basis (not equal to the initial one)
        uij_ao = np.einsum('mi,axij,nj->axmn', mo_occ, -0.5 * ovlp_deriv_oo, mo_occ).reshape((3*natm, nao, nao))
        uij_ao_list = list([uij_ao[x] for x in range(natm * 3)])

        # create AODensity and Fock matrix objects, contract with ERI
        ao_density_uij = AODensityMatrix(uij_ao_list, denmat.rest)
        fock_uij = AOFockMatrix(ao_density_uij)
        #fock_flag = fockmat.rgenjk
        #fock_uij.set_fock_type(fock_flag, 1)
        eri_drv = ElectronRepulsionIntegralsDriver(self.comm)
        screening = eri_drv.compute(get_qq_scheme(self.scf_drv.qq_type),
                                    self.scf_drv.eri_thresh, molecule, ao_basis)
        eri_drv.compute(fock_uij, ao_density_uij, molecule, ao_basis, screening)

        # TODO: how can this be done better?
        fock_uij_numpy = np.zeros((natm,3,nao,nao))
        for i in range(natm):
            for x in range(3):
                fock_uij_numpy[i,x] = fock_uij.to_numpy(3*i + x)

        # transform to MO basis
        fock_uij_mo = np.einsum('mi,axmn,nb->axib', mo_occ, fock_uij_numpy, mo_vir)

        # sum up the terms of the RHS
        cphf_rhs = fock_deriv_ov - orben_ovlp_deriv_ov + 2 * fock_uij_mo

        # Create an OrbitalResponse object
        # TODO: not needed anymore, but copy cphf_dict
        orbrsp_drv = OrbitalResponse(self.comm, self.ostream)
        orbrsp_drv.update_settings(orbrsp_dict=self.cphf_dict, method_dict=self.method_dict)
        dft_dict = orbrsp_drv.init_dft(molecule, self.scf_drv.scf_tensors)

        cphf_ov = np.zeros((natm, 3, nocc, nvir))

        # Solve the CPHF equations
        cphf_ov = self.solve_cphf(molecule, ao_basis, self.scf_drv.scf_tensors,
                                  cphf_rhs.reshape(3*natm, nocc, nvir), # TODO: possibly change the shape
                                  dft_dict, self.profiler)
        ###for i in range(natm):
        ###    for x in range(3):
        ###        # Call compute_lambda for all atoms and coordinates
        ###        cphf_ov[i,x] = orbrsp_drv.compute_lambda(molecule, ao_basis,
        ###                                                 self.scf_drv.scf_tensors,
        ###                                                 cphf_rhs[i,x], dft_dict, self.profiler)


        stop = tm.time()
        print("CPHF took %.2f s." % (stop-start) )
        return {
            'cphf_ov': cphf_ov,
            'cphf_rhs': cphf_rhs,
            'ovlp_deriv_oo': ovlp_deriv_oo,
            'fock_deriv_ao': fock_deriv_ao,
            'fock_uij': fock_uij_numpy,
        }

    def solve_cphf(self, molecule, ao_basis, scf_tensors, cphf_rhs, dft_dict, profiler):
        """
        Solves the CPHF equations for all atomic coordinates to obtain the ov block
        of the CPHF coefficients.

        :param molecule:
            The molecule.
        :param ao_basis:
            The AO basis set.
        :param scf_tensors:
            The tensors from the converged SCF calculation.
        :param cphf_rhs:
            The right-hand side of the CPHF equations for all atomic coordinates.
        :param dft_dict:
            The dictionary of DFT settings.
        :param profiler:
            The profiler.

        :returns:
            The ov block of the CPHF coefficients.
        """

        # count variable for conjugate gradient iterations
        self.iter_count = 0

        nocc = molecule.number_of_alpha_electrons()
        natm = molecule.number_of_atoms()

        if self.rank == mpi_master():
            mo = scf_tensors['C']
            nao = mo.shape[0]
            mo_energies = scf_tensors['E']

            mo_occ = mo[:, :nocc].copy()
            mo_vir = mo[:, nocc:].copy()
            nvir = mo_vir.shape[1]

            eocc = mo_energies[:nocc]
            evir = mo_energies[nocc:]
            eov = eocc.reshape(-1, 1) - evir
        else:
            nvir = None
            eov = None
        nvir = self.comm.bcast(nvir, root=mpi_master())
        eov = self.comm.bcast(eov, root=mpi_master())

        # Calculate the initial guess for the CPHF coefficients given by
        # the RHS divided by orbital-energy differences
        cphf_guess = cphf_rhs / eov

        # TODO: delete this variable
        self.initial_guess = cphf_guess

        print("CPHF guess 0:\n", cphf_guess[0])
        if self.rank == mpi_master():
            # Create AODensityMatrix object from CPHF guess in AO
            cphf_ao = np.einsum('mi,xia,na->xmn', mo_occ, cphf_guess, mo_vir)
            cphf_ao_list = list([cphf_ao[x] for x in range(3*natm)])
            # create AODensityMatrix object
            ao_density_cphf = AODensityMatrix(cphf_ao_list, denmat.rest)

            #TODO: remove this variable:
            self.cphf_ao_list = cphf_ao_list
        else:
            ao_density_cphf = AODensityMatrix()
        ao_density_cphf.broadcast(self.rank, self.comm)

        # ERI driver
        eri_drv = ElectronRepulsionIntegralsDriver(self.comm)
        screening = eri_drv.compute(get_qq_scheme(self.scf_drv.qq_type),
                                    self.scf_drv.eri_thresh, molecule, ao_basis)

        molgrid = dft_dict['molgrid']
        gs_density = dft_dict['gs_density']

        # Create a Fock Matrix Object (initialized with zeros)
        fock_cphf = AOFockMatrix(ao_density_cphf)
        fock_flag = fockmat.rgenjk

        # TODO: can this be done more efficiently, Zilvinas? (C layer)
        for i in range(3*natm):
            if self.dft:
                if self.xcfun.is_hybrid():
                    fock_flag = fockmat.rgenjkx
                    fact_xc = self.xcfun.get_frac_exact_exchange()
                    fock_cphf.set_scale_factor(fact_xc, i)
                else:
                    fock_flag = fockmat.rgenj

            fock_cphf.set_fock_type(fock_flag, i)

        # Matrix-vector product of orbital Hessian with trial vector
        def cphf_matvec(v):
            """
            Function to carry out matrix multiplication of CPHF coefficient
            vector with orbital Hessian matrix.
            """

            profiler.start_timer(self.iter_count, 'CG')

            # Create AODensityMatrix object from lambda in AO
            if self.rank == mpi_master():
                cphf_ao = np.einsum('mi,xia,na->xmn', mo_occ, v.reshape(3*natm, nocc, nvir), mo_vir)
                cphf_ao_list = list([cphf_ao[x] for x in range(3*natm)])
                ao_density_cphf = AODensityMatrix(cphf_ao_list, denmat.rest)
            else:
                ao_density_cphf = AODensityMatrix()
            ao_density_cphf.broadcast(self.rank, self.comm)

            eri_drv.compute(fock_cphf, ao_density_cphf, molecule, ao_basis,
                            screening)
            if self.dft:
                #t0 = tm.time()
                if not self.xcfun.is_hybrid():
                    for i in range(3*natm):
                        fock_cphf.scale(2.0, i)
                xc_drv = XCIntegrator(self.comm)
                molgrid.distribute(self.rank, self.nodes, self.comm)
                xc_drv.integrate(fock_cphf, ao_density_cphf, gs_density,
                                 molecule, ao_basis, molgrid,
                                 self.xcfun.get_func_label())
                #if timing_dict is not None:
                #    timing_dict['DFT'] = tm.time() - t0

            fock_cphf.reduce_sum(self.rank, self.nodes, self.comm)

            # Transform to MO basis (symmetrized w.r.t. occ. and virt.)
            # and add diagonal part
            if self.rank == mpi_master():
                fock_cphf_numpy = np.zeros((3*natm,nao,nao))
                for i in range(3*natm):
                    fock_cphf_numpy[i] = fock_cphf.to_numpy(i)

                    if i == 0 and self.iter_count == 0:
                        print("Diagonal part in MO:\n")
                        print(v.reshape(3*natm, nocc, nvir)[i])
                        print()
                        print("eov:\n")
                        print(eov)

                cphf_mo = (-np.einsum('mi,xmn,na->xia', mo_occ, fock_cphf_numpy, mo_vir)
                          - np.einsum('ma,xmn,ni->xia', mo_vir, fock_cphf_numpy, mo_occ)
                          + v.reshape(3*natm, nocc, nvir) * eov)
            else:
                cphf_mo = None

            cphf_mo = self.comm.bcast(cphf_mo, root=mpi_master())

            profiler.stop_timer(self.iter_count, 'CG')

            profiler.check_memory_usage(
                'CG Iteration {:d}'.format(self.iter_count + 1))

            profiler.print_memory_tracing(self.ostream)

            # increase iteration counter every time this function is called
            self.iter_count += 1

            return cphf_mo.reshape(3 * natm * nocc * nvir)

        print("CPHF guess:\n", cphf_guess[0])
        self.matvec_init_guess = cphf_matvec(cphf_guess.reshape(3*natm * nocc * nvir))

        # Matrix-vector product for preconditioner using the
        # inverse of the diagonal (i.e. eocc - evir)
        def precond_matvec(v):
            """
            Function that defines the matrix-vector product
            required by the pre-conditioner for the conjugate gradient.
            It is an approximation for the inverse of matrix A in Ax = b.
            """
            current_v = v.reshape(3*natm, nocc, nvir)
            M_dot_v = current_v / eov

            return M_dot_v.reshape(3 * natm * nocc * nvir)

        # 5) Define the linear operators and run conjugate gradient
        LinOp = linalg.LinearOperator((3*natm * nocc * nvir, 3*natm * nocc * nvir),
                                      matvec=cphf_matvec)
        PrecondOp = linalg.LinearOperator((3*natm * nocc * nvir, 3*natm * nocc * nvir),
                                          matvec=precond_matvec)

        b = cphf_rhs.reshape(3*natm * nocc * nvir)
        x0 = cphf_guess.reshape(3*natm * nocc * nvir)

        cphf_coefficients_ov, cg_conv = linalg.cg(A=LinOp,
                                                b=b,
                                                x0=x0,
                                                M=PrecondOp,
                                                tol=self.conv_thresh,
                                                atol=0,
                                                maxiter=self.max_iter)

        self.is_converged = (cg_conv == 0)

        return cphf_coefficients_ov.reshape(natm, 3, nocc, nvir)



    def compute_numerical(self, molecule, ao_basis, min_basis=None):
        """
        Performs calculation of numerical Hessian.

        :param molecule:
            The molecule.
        :param ao_basis:
            The AO basis set.
        :param min_basis:
            The minimal AO basis set.
        """

        # settings dictionary for gradient driver
        grad_dict = dict(self.freq_dict)
        if self.numerical_grad:
            grad_dict['numerical'] = 'yes'
            warn_msg = '*** Warning: Numerical Hessian will be calculated based on numerical gradient.'
            self.ostream.print_header(warn_msg.ljust(56))
            warn_msg = '  This takes a long time and has limited accuracy.'
            self.ostream.print_header(warn_msg.ljust(56))
            self.ostream.print_blank()
            self.ostream.flush()
        else:
            grad_dict['numerical'] = 'no'

        scf_ostream_state = self.scf_drv.ostream.state
        self.scf_drv.ostream.state = False

        # atom labels
        labels = molecule.get_labels()

        # number of atoms
        natm = molecule.number_of_atoms()

        # atom coordinates (nx3)
        coords = molecule.get_coordinates()

        # gradient driver
        grad_drv = ScfGradientDriver(self.scf_drv, self.scf_drv.comm, self.scf_drv.ostream)
        grad_drv.update_settings(grad_dict, self.method_dict)

        # number of atomic orbitals
        nao = self.scf_drv.scf_tensors['D_alpha'].shape[0]

        # Hessian
        hessian = np.zeros((natm, 3, natm, 3))

        # Gradient of the MO energies/coefficients and (energy-weighted) density matrix
        self.orben_grad = np.zeros((natm, 3, nao))
        self.mo_grad = np.zeros((natm, 3, nao, nao))
        self.density_grad = np.zeros((natm, 3, nao, nao))
        self.omega_grad = np.zeros((natm, 3, nao, nao))

        # First-order properties for gradient of dipole moment
        prop = FirstOrderProperties(self.comm, self.ostream)
        # numerical gradient (3 dipole components, no. atoms x 3 atom coords)
        #self.dipole_gradient = np.zeros((3, natm, 3))
        self.dipole_gradient = np.zeros((3, 3 * natm))

        # If Raman intensities are calculated, set up LR solver and member variable
        if self.do_raman:
            # linear response driver for polarizability calculation
            lr_drv = LinearResponseSolver(self.comm, self.scf_drv.ostream)
            #lr_ostream_state = lr_drv.ostream.state
            #lr_drv.ostream.state = False
            # polarizability: 3 coordinates x 3 coordinates (ignoring frequencies)
            # polarizability gradient: dictionary goes through 3 coordinates x 3 coordinates
            # each entry having values for no. atoms x 3 coordinates
            self.pol_gradient = np.zeros((3, 3, 3 * natm))
            # dictionary to translate from numbers to operator components 'xyz'
            component_dict = {0: 'x', 1: 'y', 2: 'z'}


        if not self.do_four_point:
            for i in range(natm):
                for d in range(3):
                    coords[i, d] += self.delta_h
                    new_mol = Molecule(labels, coords, units='au')
                    self.scf_drv.compute(new_mol, ao_basis, min_basis)
                    grad_drv.compute(new_mol, ao_basis)
                    grad_plus = grad_drv.get_gradient()

                    #density = 2.0 * self.scf_drv.scf_tensors['D_alpha']
                    orben_plus = self.scf_drv.scf_tensors['E']
                    mo_plus = self.scf_drv.scf_tensors['C_alpha']
                    # *2 for alpha+beta density
                    density_plus = 2.0 * self.scf_drv.scf_tensors['D_alpha']
                    omega_plus = grad_drv.omega_ao
                    #prop.compute(new_mol, ao_basis, density)
                    prop.compute(new_mol, ao_basis, density_plus)
                    mu_plus = prop.get_property('dipole moment')

                    if self.do_raman:
                        lr_drv.is_converged = False
                        lr_results_p = lr_drv.compute(new_mol, ao_basis,
                                                           self.scf_drv.scf_tensors)


                    coords[i, d] -= 2.0 * self.delta_h
                    new_mol = Molecule(labels, coords, units='au')
                    self.scf_drv.compute(new_mol, ao_basis, min_basis)
                    grad_drv.compute(new_mol, ao_basis)
                    grad_minus = grad_drv.get_gradient()

                    #density = 2.0 * self.scf_drv.scf_tensors['D_alpha']
                    orben_minus = self.scf_drv.scf_tensors['E']
                    mo_minus = self.scf_drv.scf_tensors['C_alpha']
                    # *2 for alpha+beta density
                    density_minus = 2.0 * self.scf_drv.scf_tensors['D_alpha']
                    omega_minus = grad_drv.omega_ao
                    #prop.compute(new_mol, ao_basis, density)
                    prop.compute(new_mol, ao_basis, density_minus)
                    mu_minus = prop.get_property('dipole moment')

                    if self.do_raman:
                        lr_drv.is_converged = False
                        lr_results_m = lr_drv.compute(new_mol, ao_basis,
                                                           self.scf_drv.scf_tensors)
                        for aop in range(3):
                            #for bop in lr_drv.b_components:
                            for bop in range(3):
                                self.pol_gradient[aop, bop, 3*i + d] = (
                                    ( lr_results_p['response_functions'][component_dict[aop], component_dict[bop], 0.0]
                                    - lr_results_m['response_functions'][component_dict[aop], component_dict[bop], 0.0] ) /
                                    (2.0 * self.delta_h) )


                    for c in range(3):
                        self.dipole_gradient[c, 3*i + d] = (mu_plus[c] - mu_minus[c]) / (2.0 * self.delta_h)
                    coords[i, d] += self.delta_h
                    hessian[i, d, :, :] = (grad_plus - grad_minus) / (2.0 * self.delta_h)
                    self.orben_grad[i, d, :] = (orben_plus - orben_minus) / (2.0 * self.delta_h)
                    self.mo_grad[i, d, :, :] = (mo_plus - mo_minus) / (2.0 * self.delta_h)
                    self.density_grad[i, d, :, :] = (density_plus - density_minus) / (2.0 * self.delta_h)
                    self.omega_grad[i, d, :, :] = (omega_plus - omega_minus) / (2.0 * self.delta_h)


        else:
            # Four-point numerical derivative approximation
            # for debugging of analytical Hessian:
            # [ f(x - 2h) - 8 f(x - h) + 8 f(x + h) - f(x + 2h) ] / ( 12h )
            for i in range(natm):
                for d in range(3):
                    coords[i, d] += self.delta_h
                    new_mol = Molecule(labels, coords, units='au')
                    self.scf_drv.compute(new_mol, ao_basis, min_basis)
                    grad_drv.compute(new_mol, ao_basis)
                    grad_plus1 = grad_drv.get_gradient()

                    density = 2.0 * self.scf_drv.scf_tensors['D_alpha']
                    prop.compute(new_mol, ao_basis, density)
                    mu_plus1 = prop.get_property('dipole moment')

                    if self.do_raman:
                        lr_drv.is_converged = False
                        lr_results_p1 = lr_drv.compute(new_mol, ao_basis,
                                                           self.scf_drv.scf_tensors)

                    coords[i, d] += self.delta_h
                    new_mol = Molecule(labels, coords, units='au')
                    self.scf_drv.compute(new_mol, ao_basis, min_basis)
                    grad_drv.compute(new_mol, ao_basis)
                    grad_plus2 = grad_drv.get_gradient()

                    density = 2.0 * self.scf_drv.scf_tensors['D_alpha']
                    prop.compute(new_mol, ao_basis, density)
                    mu_plus2 = prop.get_property('dipole moment')

                    if self.do_raman:
                        lr_drv.is_converged = False
                        lr_results_p2 = lr_drv.compute(new_mol, ao_basis,
                                                           self.scf_drv.scf_tensors)

                    coords[i, d] -= 3.0 * self.delta_h
                    new_mol = Molecule(labels, coords, units='au')
                    self.scf_drv.compute(new_mol, ao_basis, min_basis)
                    grad_drv.compute(new_mol, ao_basis)
                    grad_minus1 = grad_drv.get_gradient()

                    density = 2.0 * self.scf_drv.scf_tensors['D_alpha']
                    prop.compute(new_mol, ao_basis, density)
                    mu_minus1 = prop.get_property('dipole moment')

                    if self.do_raman:
                        lr_drv.is_converged = False
                        lr_results_m1 = lr_drv.compute(new_mol, ao_basis,
                                                           self.scf_drv.scf_tensors)

                    coords[i, d] -= self.delta_h
                    new_mol = Molecule(labels, coords, units='au')
                    self.scf_drv.compute(new_mol, ao_basis, min_basis)
                    grad_drv.compute(new_mol, ao_basis)
                    grad_minus2 = grad_drv.get_gradient()

                    density = 2.0 * self.scf_drv.scf_tensors['D_alpha']
                    prop.compute(new_mol, ao_basis, density)
                    mu_minus2 = prop.get_property('dipole moment')

                    if self.do_raman:
                        lr_drv.is_converged = False
                        lr_results_m2 = lr_drv.compute(new_mol, ao_basis,
                                                           self.scf_drv.scf_tensors)

                        for aop in range(3):
                            for bop in range(3):
                                self.pol_gradient[aop, bop, 3*i + d] = (
                                    ( lr_results_m2['response_functions'][component_dict[aop], component_dict[bop], 0.0]
                                    - 8 * lr_results_m1['response_functions'][component_dict[aop], component_dict[bop], 0.0]
                                    + 8 * lr_results_p1['response_functions'][component_dict[aop], component_dict[bop], 0.0]
                                    - lr_results_p2['response_functions'][component_dict[aop], component_dict[bop], 0.0] ) /
                                    (12.0 * self.delta_h) )

                    for c in range(3):
                        self.dipole_gradient[c, 3*i + d] = (mu_minus2[c] - 8.0 * mu_minus1[c]
                                                         + 8.0 * mu_plus1[c] - mu_plus2[c]) / (12.0 * self.delta_h)
                    coords[i, d] += 2.0 * self.delta_h
                    # f'(x) ~ [ f(x - 2h) - 8 f(x - h) + 8 f(x + h) - f(x + 2h) ] / ( 12h )
                    hessian[i, d] = (grad_minus2 - 8.0 * grad_minus1
                                           + 8.0 * grad_plus1 - grad_plus2) / (12.0 * self.delta_h)

        # reshaped Hessian as member variable
        self.hessian = hessian.reshape(3*natm, 3*natm)

        #self.ostream.print_blank()

        self.scf_drv.compute(molecule, ao_basis, min_basis)
        self.scf_drv.ostream.state = scf_ostream_state

    def compute_dipole_gradient(self, molecule, ao_basis, perturbed_density):
        """
        Computes the analytical gradient of the dipole moment.

        :param molecule:
            The molecule.
        :param ao_basis:
            The AO basis set.
        :param perturbed_density:
            The perturbed density matrix.
        """

        # Number of atoms and atomic charges
        natm = molecule.number_of_atoms()
        nuclear_charges = molecule.elem_ids_to_numpy()

        density = self.scf_drv.scf_tensors['D_alpha']

        # Dipole integrals
        dipole_drv = ElectricDipoleIntegralsDriver(self.comm)
        dipole_mats = dipole_drv.compute(molecule, ao_basis)
        dipole_ints = np.array((dipole_mats.x_to_numpy(), dipole_mats.y_to_numpy(),
                       dipole_mats.z_to_numpy()))

        # Initialize a local dipole gradient to zero
        dipole_gradient = np.zeros((3, natm, 3))

        # Put the nuclear contributions to the right place
        natm_zeros = np.zeros((natm))
        dipole_gradient[0] = np.vstack((nuclear_charges, natm_zeros, natm_zeros)).T
        dipole_gradient[1] = np.vstack((natm_zeros, nuclear_charges, natm_zeros)).T
        dipole_gradient[2] = np.vstack((natm_zeros, natm_zeros, nuclear_charges)).T

        # TODO: replace once analytical integral derivatives are available
        dipole_integrals_deriv = self.compute_dipole_integral_derivatives(molecule, ao_basis)

        # Add the electronic contributions
        dipole_gradient += -2 * (np.einsum('mn,caxmn->cax', density, dipole_integrals_deriv)
                           + np.einsum('axmn,cmn->cax', perturbed_density, dipole_ints)
                           )

        self.dipole_gradient = dipole_gradient.reshape(3, 3 * natm)


    def compute_dipole_integral_derivatives(self, molecule, ao_basis):
        """
        Computes numerical derivatives of dipole integrals.

        :param molecule:
            The molecule.
        :param ao_basis:
            The AO basis set.

        :return:
            The dipole integral derivatives.
        """

        # atom labels
        labels = molecule.get_labels()

        # number of atoms
        natm = molecule.number_of_atoms()

        # atom coordinates (nx3)
        coords = molecule.get_coordinates()

        # number of atomic orbitals
        nao = self.scf_drv.scf_tensors['D_alpha'].shape[0]

        # Dipole integrals driver
        dipole_drv = ElectricDipoleIntegralsDriver(self.comm)

        # 3 dipole components x No. atoms x 3 atomic coordinates x No. basis x No. basis
        dipole_integrals_gradient = np.zeros((3, natm, 3, nao, nao))

        # smaller delta_h values can be used here
        local_delta_h = 0.01 * self.delta_h

        for i in range(natm):
            for d in range(3):
                coords[i, d] += local_delta_h
                new_mol = Molecule(labels, coords, units='au')

                dipole_mats_p = dipole_drv.compute(new_mol, ao_basis)
                dipole_ints_p = (dipole_mats_p.x_to_numpy(), dipole_mats_p.y_to_numpy(),
                               dipole_mats_p.z_to_numpy())

                coords[i, d] -= 2.0 * local_delta_h
                new_mol = Molecule(labels, coords, units='au')

                dipole_mats_m = dipole_drv.compute(new_mol, ao_basis)
                dipole_ints_m = (dipole_mats_m.x_to_numpy(), dipole_mats_m.y_to_numpy(),
                               dipole_mats_m.z_to_numpy())

                for c in range(3):
                    dipole_integrals_gradient[c, i, d] = ( dipole_ints_p[c] - dipole_ints_m[c] ) / (2.0 * local_delta_h)

        return dipole_integrals_gradient


