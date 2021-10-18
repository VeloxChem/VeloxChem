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
from .veloxchemlib import ElectronRepulsionIntegralsDriver
from .veloxchemlib import denmat
from .veloxchemlib import fockmat
from .veloxchemlib import AODensityMatrix
from .veloxchemlib import AOFockMatrix
from .veloxchemlib import ElectricDipoleIntegralsDriver

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
        overlap = self.scf_drv.scf_tensors['S']
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
        orbrsp_drv = OrbitalResponse(self.comm, self.ostream)
        orbrsp_drv.update_settings(orbrsp_dict=self.cphf_dict, method_dict=self.method_dict)
        dft_dict = orbrsp_drv.init_dft(molecule, self.scf_drv.scf_tensors)

        cphf_ov = np.zeros((natm, 3, nocc, nvir))

        for i in range(natm):
            for x in range(3):
                # Call compute_lambda for all atoms and coordinates
                cphf_ov[i,x] = orbrsp_drv.compute_lambda(molecule, ao_basis,
                                                         self.scf_drv.scf_tensors,
                                                         cphf_rhs[i,x], dft_dict, self.profiler)


        if self.pople:
            self.compute_pople(molecule, ao_basis, -0.5 * ovlp_deriv_oo, cphf_ov, fock_uij_numpy)
        else:
            hessian_first_order_derivatives = self.compute_furche(molecule, ao_basis,
                                                         cphf_rhs, -0.5 * ovlp_deriv_oo, cphf_ov)

        hessian_2nd_order_derivatives = np.zeros((natm, natm, 3, 3))
        for i in range(natm):
            for j in range(natm):
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

        hessian_nuclear_nuclear = self.hess_nuc_contrib(molecule)

        self.hessian = ( hessian_first_order_derivatives + hessian_2nd_order_derivatives
                       + hessian_nuclear_nuclear ).transpose(0,2,1,3).reshape(3*natm, 3*natm)


    def compute_pople(self, molecule, ao_basis, cphf_oo, cphf_ov, fock_uij):
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
        perturbed_density = ( 2 * np.einsum('mj,xyij,ni->xymn',
                                        mo_occ, cphf_oo, mo_occ)
                              + np.einsum('ma,xyia,ni->xymn',
                                        mo_vir, cphf_ov, mo_occ)
                              + np.einsum('mi,xyia,na->xymn',
                                        mo_occ, cphf_ov, mo_vir)
                            )
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

        fock_uia_numpy = np.zeros((natm,3,nao,nao))
        for i in range(natm):
            for x in range(3):
                fock_uia_numpy[i,x] = fock_uia.to_numpy(3*i + x)

        fock_cphf_oo = np.einsum('mi,xymn,nj->xyij', mo_occ, fock_uij, mo_occ)
        
        fock_cphf_ov = ( np.einsum('mi,xymn,nj->xyij', mo_occ, fock_uia_numpy, mo_occ)
                        +np.einsum('mj,xymn,ni->xyij', mo_occ, fock_uia_numpy, mo_occ)
                        )
        
        # TODO: derivative of "epsilon" is missing. CPHF coefficients need to be contracted with ERI


        # analytical gradient
        ###self.gradient = np.zeros((natm, 3))

        ###for i in range(natm):
        ###    d_ovlp = overlap_deriv(molecule, ao_basis, i)
        ###    d_fock = fock_deriv(molecule, ao_basis, one_pdm_ao, i)
        ###    d_eri = eri_deriv(molecule, ao_basis, i)

        ###    self.gradient[i] += ( 2.0*np.einsum('mn,xmn->x', one_pdm_ao, d_fock)
        ###                    +2.0*np.einsum('mn,xmn->x', epsilon_dm_ao, d_ovlp)
        ###                    -2.0*np.einsum('mt,np,xmtnp->x', one_pdm_ao, one_pdm_ao, d_eri)
        ###                    +1.0*np.einsum('mt,np,xmnpt->x', one_pdm_ao, one_pdm_ao, d_eri)
        ###                    )

        ###self.gradient += self.grad_nuc_contrib(molecule)

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
        hessian_cphf_coeff_rhs = np.zeros((natm, natm, 3, 3))
        for i in range(natm):
            for j in range(natm):
                hessian_cphf_coeff_rhs[i,j] += 4*np.einsum('xia,yia->xy', cphf_ov[i], cphf_rhs[j])

        # First integral derivatives: partial Fock and overlap matrix derivatives
        hessian_first_integral_derivatives = np.zeros((natm, natm, 3, 3))
        for i in range(natm):
            fock_deriv_i = fock_deriv(molecule, ao_basis, density, i)
            ovlp_deriv_i = overlap_deriv(molecule, ao_basis, i)
            for j in range(natm):
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

        # Overlap derivative with ERIs
        hessian_eri_overlap = np.zeros((natm, natm, 3, 3))
        for i in range(natm):
            ovlp_deriv_i = overlap_deriv(molecule, ao_basis, i)
            for j in range(natm):
                ovlp_deriv_j = overlap_deriv(molecule, ao_basis, j)

                # Overlap derivative contracted with two density matrices
                P_P_Six = np.zeros((3, nao, nao))
                P_P_Sjy = np.zeros((3, nao, nao))
                for x in range(3):
                    P_P_Six[x] = (np.linalg.multi_dot([density, ovlp_deriv_i[x], density]))
                    P_P_Sjy[x] = (np.linalg.multi_dot([density, ovlp_deriv_j[x], density]))

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

        # return the sum of the three contributions
        return (hessian_cphf_coeff_rhs + hessian_first_integral_derivatives
                         + hessian_eri_overlap)


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

        for i in range(natm):
            for d in range(3):
                coords[i, d] += self.delta_h
                new_mol = Molecule(labels, coords, units='au')

                dipole_mats_p = dipole_drv.compute(new_mol, ao_basis)
                dipole_ints_p = (dipole_mats_p.x_to_numpy(), dipole_mats_p.y_to_numpy(),
                               dipole_mats_p.z_to_numpy())

                coords[i, d] -= 2.0 * self.delta_h
                new_mol = Molecule(labels, coords, units='au')

                dipole_mats_m = dipole_drv.compute(new_mol, ao_basis)
                dipole_ints_m = (dipole_mats_m.x_to_numpy(), dipole_mats_m.y_to_numpy(),
                               dipole_mats_m.z_to_numpy())

                for c in range(3):
                    dipole_integrals_gradient[c, i, d] = ( dipole_ints_p[c] - dipole_ints_m[c] ) / (2.0 * self.delta_h)

        return dipole_integrals_gradient


