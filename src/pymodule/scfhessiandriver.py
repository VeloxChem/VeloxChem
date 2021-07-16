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

# For PySCF integral derivatives
from .import_from_pyscf import overlap_deriv
from .import_from_pyscf import fock_deriv
from .import_from_pyscf import eri_deriv


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
    """

    def __init__(self, scf_drv, comm=None, ostream=None):
        """
        Initializes SCF Hessian driver.
        """

        super().__init__(comm, ostream)

        self.flag = 'SCF Hessian Driver'
        self.scf_drv = scf_drv


    def update_settings(self, method_dict, freq_dict=None):
        """
        Updates settings in ScfHessianDriver.

        :param method_dict:
            The input dictionary of method settings group.
        :param freq_dict:
            The input dictionary of Hessian/frequency settings group.
        """

        super().update_settings(method_dict, freq_dict)

        # The electronic energy
        self.elec_energy = self.scf_drv.get_scf_energy()


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

        if True: #self.numerical:
            self.compute_numerical(molecule, ao_basis, min_basis)
        #else:
        #    self.compute_analytical(molecule, ao_basis)

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

###    def compute_analytical(self, molecule, ao_basis):
###        """
###        Computes the analytical nuclear gradient.
###        So far only for restricted Hartree-Fock with PySCF integral derivatives...
###
###        :param molecule:
###            The molecule.
###        :param ao_basis:
###            The AO basis set.
###        """
###        natm = molecule.number_of_atoms()
###        nocc = molecule.number_of_alpha_electrons()
###        mo = self.scf_drv.scf_tensors['C']
###        mo_occ = mo[:, :nocc].copy()
###        one_pdm_ao = self.scf_drv.scf_tensors['D_alpha']
###        mo_energies = self.scf_drv.scf_tensors['E']
###        eocc = mo_energies[:nocc]
###        eo_diag = np.diag(eocc)
###        epsilon_dm_ao = - np.linalg.multi_dot([mo_occ, eo_diag, mo_occ.T])
###
###        # analytical gradient
###        self.gradient = np.zeros((natm, 3))
###
###        for i in range(natm):
###            d_ovlp = overlap_deriv(molecule, ao_basis, i)
###            d_fock = fock_deriv(molecule, ao_basis, one_pdm_ao, i)
###            d_eri = eri_deriv(molecule, ao_basis, i)
###
###            self.gradient[i] += ( 2.0*np.einsum('mn,xmn->x', one_pdm_ao, d_fock)
###                            +2.0*np.einsum('mn,xmn->x', epsilon_dm_ao, d_ovlp)
###                            -2.0*np.einsum('mt,np,xmtnp->x', one_pdm_ao, one_pdm_ao, d_eri)
###                            +1.0*np.einsum('mt,np,xmnpt->x', one_pdm_ao, one_pdm_ao, d_eri)
###                            )
###
###        self.gradient += self.grad_nuc_contrib(molecule)

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

        # Hessian
        hessian = np.zeros((natm, 3, natm, 3))

        # First-order properties for gradient of dipole moment
        prop = FirstOrderProperties(self.comm, self.ostream)
        # numerical gradient (3 dipole components, no. atoms x 3 atom coords)
        #self.dipole_gradient = np.zeros((3, molecule.number_of_atoms(), 3))
        self.dipole_gradient = np.zeros((3, 3 * molecule.number_of_atoms()))

        if not self.do_four_point:
            for i in range(molecule.number_of_atoms()):
                for d in range(3):
                    coords[i, d] += self.delta_h
                    new_mol = Molecule(labels, coords, units='au')
                    self.scf_drv.compute(new_mol, ao_basis, min_basis)
                    grad_drv.compute(new_mol, ao_basis)
                    grad_plus = grad_drv.get_gradient()

                    density = 2.0 * self.scf_drv.scf_tensors['D_alpha']
                    prop.compute(new_mol, ao_basis, density)
                    mu_plus = prop.get_property('dipole moment')


                    coords[i, d] -= 2.0 * self.delta_h
                    new_mol = Molecule(labels, coords, units='au')
                    self.scf_drv.compute(new_mol, ao_basis, min_basis)
                    grad_drv.compute(new_mol, ao_basis)
                    grad_minus = grad_drv.get_gradient()

                    density = 2.0 * self.scf_drv.scf_tensors['D_alpha']
                    prop.compute(new_mol, ao_basis, density)
                    mu_minus = prop.get_property('dipole moment')

                    for c in range(3):
                        self.dipole_gradient[c, 3*i + d] = (mu_plus[c] - mu_minus[c]) / (2.0 * self.delta_h)
                    coords[i, d] += self.delta_h
                    hessian[i, d, :, :] = (grad_plus - grad_minus) / (2.0 * self.delta_h)


        else:
            # Four-point numerical derivative approximation
            # for debugging of analytical Hessian:
            # [ f(x - 2h) - 8 f(x - h) + 8 f(x + h) - f(x + 2h) ] / ( 12h )
            for i in range(molecule.number_of_atoms()):
                for d in range(3):
                    coords[i, d] += self.delta_h
                    new_mol = Molecule(labels, coords, units='au')
                    self.scf_drv.compute(new_mol, ao_basis, min_basis)
                    grad_drv.compute(new_mol, ao_basis)
                    grad_plus1 = grad_drv.get_gradient()

                    density = 2.0 * self.scf_drv.scf_tensors['D_alpha']
                    prop.compute(new_mol, ao_basis, density)
                    mu_plus1 = prop.get_property('dipole moment')

                    coords[i, d] += self.delta_h
                    new_mol = Molecule(labels, coords, units='au')
                    self.scf_drv.compute(new_mol, ao_basis, min_basis)
                    grad_drv.compute(new_mol, ao_basis)
                    grad_plus2 = grad_drv.get_gradient()

                    density = 2.0 * self.scf_drv.scf_tensors['D_alpha']
                    prop.compute(new_mol, ao_basis, density)
                    mu_plus2 = prop.get_property('dipole moment')

                    coords[i, d] -= 3.0 * self.delta_h
                    new_mol = Molecule(labels, coords, units='au')
                    self.scf_drv.compute(new_mol, ao_basis, min_basis)
                    grad_drv.compute(new_mol, ao_basis)
                    grad_minus1 = grad_drv.get_gradient()

                    density = 2.0 * self.scf_drv.scf_tensors['D_alpha']
                    prop.compute(new_mol, ao_basis, density)
                    mu_minus1 = prop.get_property('dipole moment')

                    coords[i, d] -= self.delta_h
                    new_mol = Molecule(labels, coords, units='au')
                    self.scf_drv.compute(new_mol, ao_basis, min_basis)
                    grad_drv.compute(new_mol, ao_basis)
                    grad_minus2 = grad_drv.get_gradient()

                    density = 2.0 * self.scf_drv.scf_tensors['D_alpha']
                    prop.compute(new_mol, ao_basis, density)
                    mu_minus2 = prop.get_property('dipole moment')

                    for c in range(3):
                        self.dipole_gradient[c, i, d] = (mu_minus2[c] - 8.0 * mu_minus1[c]
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

