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

import time as tm
import numpy as np

from .molecule import Molecule
from .gradientdriver import GradientDriver
from .outputstream import OutputStream
from .firstorderprop import FirstOrderProperties
from .veloxchemlib import mpi_master

# For PySCF integral derivatives
from .import_from_pyscf import overlap_deriv
from .import_from_pyscf import fock_deriv
from .import_from_pyscf import eri_deriv
from .import_from_pyscf import hcore_deriv
from .import_from_pyscf import vxc_deriv

class ScfGradientDriver(GradientDriver):
    """
    Implements SCF gradient driver.

    :param scf_drv:
        The SCF driver.
    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.

    Instance variables
        - scf_drv: The SCF driver.
        - delta_h: The displacement for finite difference.
    """

    def __init__(self, scf_drv, comm=None, ostream=None):
        """
        Initializes gradient driver.
        """

        super().__init__(scf_drv, comm, ostream)

        self.flag = 'SCF Gradient Driver'
        self.scf_drv = scf_drv
        self.delta_h = 0.001
        self.add_xc_grad = True


    def update_settings(self, grad_dict, method_dict):
        """
        Updates settings in ScfGradientDriver.

        :param grad_dict:
            The input dictionary of gradient settings group.
        :param method_dict:
            The input dicitonary of method settings group.
        """

        # update settings in parent class
        super().update_settings(grad_dict, method_dict)
        if 'add_xc_grad' in grad_dict:
            key = grad_dict['add_xc_grad']
            if key in ['n', 'no']:
                self.add_xc_grad = False


    def compute(self, molecule, ao_basis, min_basis=None):
        """
        Computes the analytical or numerical nuclear gradient.

        :param molecule:
            The molecule.
        :param ao_basis:
            The AO basis set.
        :param min_basis:
            The minimal AO basis set.
        """
        if self.rank == mpi_master():
            self.print_header()

        start_time = tm.time()

        if self.numerical:
            scf_ostream_state = self.scf_drv.ostream.state
            self.scf_drv.ostream.state = False
            self.compute_numerical(molecule, ao_basis, min_basis)
            self.scf_drv.ostream.state = scf_ostream_state

        else:
            self.compute_analytical(molecule, ao_basis)

        if self.rank == mpi_master():
            # print gradient
            self.print_geometry(molecule)
            self.print_gradient(molecule)

            valstr = '*** Time spent in gradient calculation: '
            valstr += '{:.2f} sec ***'.format(tm.time() - start_time)
            self.ostream.print_header(valstr)
            self.ostream.print_blank()
            self.ostream.flush()

    def compute_analytical(self, molecule, ao_basis):
        """
        Computes the analytical nuclear gradient.
        So far only for restricted Hartree-Fock with PySCF integral derivatives...

        :param molecule:
            The molecule.
        :param ao_basis:
            The AO basis set.
        """
        if self.rank == mpi_master():
            natm = molecule.number_of_atoms()
            nocc = molecule.number_of_alpha_electrons()
            mo = self.scf_drv.scf_tensors['C_alpha']
            mo_occ = mo[:, :nocc].copy()
            one_pdm_ao = self.scf_drv.scf_tensors['D_alpha']
            mo_energies = self.scf_drv.scf_tensors['E_alpha']
            eocc = mo_energies[:nocc]
            eo_diag = np.diag(eocc)
            epsilon_dm_ao = - np.linalg.multi_dot([mo_occ, eo_diag, mo_occ.T])

            # analytical gradient
            self.gradient = np.zeros((natm, 3))
            # exc_gradient = np.zeros((natm, 3))

            for i in range(natm):
                d_ovlp = overlap_deriv(molecule, ao_basis, i)
                #d_fock = fock_deriv(molecule, ao_basis, one_pdm_ao, i) #TODO: remove
                d_hcore = hcore_deriv(molecule, ao_basis, i)
                d_eri = eri_deriv(molecule, ao_basis, i)

                self.gradient[i] = ( 2.0 * np.einsum('mn,xmn->x', one_pdm_ao, d_hcore)
                                + 2.0 * np.einsum('mn,xmn->x', epsilon_dm_ao, d_ovlp)
                                )
                if self.dft:
                    d_vxc = vxc_deriv(molecule, ao_basis, one_pdm_ao,
                                      self.xcfun.get_func_label(), i,
                                      self.grid_level)
                    if self.xcfun.is_hybrid():
                        fact_xc = self.xcfun.get_frac_exact_exchange()
                    else:
                        fact_xc = 0
                    self.gradient[i] += (
                                + 2.0 * np.einsum('mt,np,xmtnp->x', one_pdm_ao, one_pdm_ao, d_eri)
                                - fact_xc * np.einsum('mt,np,xmnpt->x', one_pdm_ao, one_pdm_ao, d_eri)
                                )
                    if self.add_xc_grad:
                        self.gradient[i] += 2.0 * np.einsum('mn,xmn->x', one_pdm_ao, d_vxc)
                else:
                    self.gradient[i] += (
                                + 2.0 * np.einsum('mt,np,xmtnp->x', one_pdm_ao, one_pdm_ao, d_eri)
                                - 1.0 * np.einsum('mt,np,xmnpt->x', one_pdm_ao, one_pdm_ao, d_eri)
                                )

            self.gradient += self.grad_nuc_contrib(molecule)
            self.omega_ao = epsilon_dm_ao #TODO remove again

# TODO: delete commented out code; numerical gradient is now included
# in gradientdriver.
#    def compute_numerical(self, molecule, ao_basis, min_basis=None):
#        """
#        Performs calculation of numerical gradient.
#
#        :param molecule:
#            The molecule.
#        :param ao_basis:
#            The AO basis set.
#        :param min_basis:
#            The minimal AO basis set.
#        """
#
#        scf_ostream_state = self.scf_drv.ostream.state
#        self.scf_drv.ostream.state = False
#
#        # atom labels
#        labels = molecule.get_labels()
#
#        # atom coordinates (nx3)
#        coords = molecule.get_coordinates()
#        # the number of atoms
#        natm = molecule.number_of_atoms()
#        # the number of atomic orbitals
#        # nao = self.scf_drv.scf_tensors["D_alpha"].shape[0]
#
#        self.gradient = np.zeros((natm, 3))
#
#        # Gradient of the exchange and correlation energy
#        # TODO: remove vxc numerical gradient
#        if self.dft:
#            self.xc_gradient = np.zeros((natm, 3))
#            # self.vxc_gradient = np.zeros((natm, 3, nao, nao))
#
#        # First-order properties for gradient of dipole moment
#        if self.dipole_deriv:
#            prop = FirstOrderProperties(self.comm, self.ostream)
#            # numerical gradient (3 dipole components x no. atoms x 3 atom coords)
#            self.dipole_gradient = np.zeros((3, natm, 3))
#
#        if not self.do_four_point:
#            for i in range(natm):
#                for d in range(3):
#                    coords[i, d] += self.delta_h
#                    new_mol = Molecule(labels, coords, units='au')
#                    self.scf_drv.compute(new_mol, ao_basis, min_basis)
#                    e_plus = self.scf_drv.get_scf_energy()
#                    # TODO: numerical derivative of the xc energy;
#                    # remove when analytical derivative is working.
#                    if self.dft:
#                        xc_plus = self.scf_drv.scf_tensors['xc_energy']
#                        # vxc_plus = self.scf_drv.scf_tensors['vxc_mat'].get_matrix().to_numpy()
#
#                    if self.dipole_deriv:
#                        prop.compute_scf_prop(new_mol, ao_basis, self.scf_drv.scf_tensors)
#                        mu_plus = prop.get_property('dipole moment')
#
#                    coords[i, d] -= 2.0 * self.delta_h
#                    new_mol = Molecule(labels, coords, units='au')
#                    self.scf_drv.compute(new_mol, ao_basis, min_basis)
#                    e_minus = self.scf_drv.get_scf_energy()
#                    # TODO: numerical derivative of the xc energy;
#                    # remove when analytical derivative is working.
#                    if self.dft:
#                        xc_minus = self.scf_drv.scf_tensors['xc_energy']
#                        # vxc_minus = self.scf_drv.scf_tensors['vxc_mat'].get_matrix().to_numpy()
#
#                    if self.dipole_deriv:
#                        prop.compute_scf_prop(new_mol, ao_basis, self.scf_drv.scf_tensors)
#                        mu_minus = prop.get_property('dipole moment')
#
#                        for c in range(3):
#                            self.dipole_gradient[c, i, d] = (mu_plus[c] - mu_minus[c]) / (2.0 * self.delta_h)
#
#                    coords[i, d] += self.delta_h
#                    self.gradient[i, d] = (e_plus - e_minus) / (2.0 * self.delta_h)
#                    # TODO: numerical derivative of the xc energy;
#                    # remove when analytical derivative is working.
#                    if self.dft:
#                        self.xc_gradient[i, d] = (xc_plus - xc_minus) / (2.0 * self.delta_h)
#                        # self.vxc_gradient[i, d] = (vxc_plus - vxc_minus) / (2.0 * self.delta_h)
#        else:
#            # Four-point numerical derivative approximation
#            # for debugging of analytical gradient:
#            # [ f(x - 2h) - 8 f(x - h) + 8 f(x + h) - f(x + 2h) ] / ( 12h )
#            for i in range(natm):
#                for d in range(3):
#                    coords[i, d] += self.delta_h
#                    new_mol = Molecule(labels, coords, units='au')
#                    self.scf_drv.compute(new_mol, ao_basis, min_basis)
#                    e_plus1 = self.scf_drv.get_scf_energy()
#
#                    if self.dipole_deriv:
#                        prop.compute_scf_prop(new_mol, ao_basis, self.scf_drv.scf_tensors)
#                        mu_plus1 = prop.get_property('dipole moment')
#
#                    coords[i, d] += self.delta_h
#                    new_mol = Molecule(labels, coords, units='au')
#                    self.scf_drv.compute(new_mol, ao_basis, min_basis)
#                    e_plus2 = self.scf_drv.get_scf_energy()
#
#                    if self.dipole_deriv:
#                        prop.compute_scf_prop(new_mol, ao_basis, self.scf_drv.scf_tensors)
#                        mu_plus2 = prop.get_property('dipole moment')
#
#                    coords[i, d] -= 3.0 * self.delta_h
#                    new_mol = Molecule(labels, coords, units='au')
#                    self.scf_drv.compute(new_mol, ao_basis, min_basis)
#                    e_minus1 = self.scf_drv.get_scf_energy()
#
#                    if self.dipole_deriv:
#                        prop.compute_scf_prop(new_mol, ao_basis, self.scf_drv.scf_tensors)
#                        mu_minus1 = prop.get_property('dipole moment')
#
#                    coords[i, d] -= self.delta_h
#                    new_mol = Molecule(labels, coords, units='au')
#                    self.scf_drv.compute(new_mol, ao_basis, min_basis)
#                    e_minus2 = self.scf_drv.get_scf_energy()
#
#                    if self.dipole_deriv:
#                        prop.compute_scf_prop(new_mol, ao_basis, self.scf_drv.scf_tensors)
#                        mu_minus2 = prop.get_property('dipole moment')
#
#                        for c in range(3):
#                            self.dipole_gradient[c, i, d] = (mu_minus2[c] - 8.0 * mu_minus1[c]
#                                                             + 8.0 * mu_plus1[c] - mu_plus2[c]) / (12.0 * self.delta_h)
#
#
#                    coords[i, d] += 2.0 * self.delta_h
#                    # f'(x) ~ [ f(x - 2h) - 8 f(x - h) + 8 f(x + h) - f(x + 2h) ] / ( 12h )
#                    self.gradient[i, d] = (e_minus2 - 8.0 * e_minus1
#                                           + 8.0 * e_plus1 - e_plus2) / (12.0 * self.delta_h)
#
#        self.ostream.print_blank() # TODO: figure out if this is needed here.
#
#        self.scf_drv.compute(molecule, ao_basis, min_basis)
#        self.scf_drv.ostream.state = scf_ostream_state

    def compute_energy(self, molecule, ao_basis, min_basis=None):
        """
        Computes the energy at current geometry.

        :param molecule:
            The molecule.
        :param ao_basis:
            The AO basis set.
        :param min_basis:
            The minimal AO basis set.

        :return:
            The energy.
        """

        self.scf_drv.compute(molecule, ao_basis, min_basis)
        return self.scf_drv.get_scf_energy()
