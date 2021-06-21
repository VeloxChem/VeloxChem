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
from .scfgradientdriver import ScfGradientDriver
from .outputstream import OutputStream
#from .firstorderprop import FirstOrderProperties

# For PySCF integral derivatives
from .import_from_pyscf import overlap_deriv
from .import_from_pyscf import fock_deriv
from .import_from_pyscf import eri_deriv


class ScfHessianDriver:
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
        - delta_h: The displacement for finite difference.
    """

    def __init__(self, scf_drv, comm=None, ostream=None):
        """
        Initializes gradient driver.
        """

        if comm is None:
            self.comm = MPI.COMM_WORLD

        if ostream is None:
            self.ostream = OutputStream(sys.stdout)

        self.flag = 'SCF Hessian Driver'
        self.hessian = None
        self.scf_drv = scf_drv
        self.delta_h = 0.001
        self.numerical = True #False
        # flag for two-point or four-point approximation
        self.do_four_point = False


        # Flag for numerical derivative of dipole moment
        self.dipole_deriv = False
        self.dipole_gradient = None

    def update_settings(self, method_dict):
        """
        Updates settings in ScfHessianDriver.

        TODO: Add new gradient settings group?

        :param method_dict:
            The input dicitonary of method settings group.
        """
        # Analytical DFT gradient is not implemented yet
        if 'xcfun' in method_dict:
            if method_dict['xcfun'] is not None:
                self.numerical = True
        if 'dft' in method_dict:
            key = method_dict['dft'].lower()
            self.numerical = True if key in ['yes', 'y'] else False

        # TODO: possibly put settings in new input group
        if 'numerical_grad' in method_dict:
            key = method_dict['numerical_grad'].lower()
            self.numerical = True if key in ['yes', 'y'] else False


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
        self.print_header()

        start_time = tm.time()

        if True: #self.numerical:
            self.compute_numerical(molecule, ao_basis, min_basis)
        #else:
        #    self.compute_analytical(molecule, ao_basis)

        # print gradient
        self.print_geometry(molecule)
        self.print_gradient(molecule)

        valstr = '*** Time spent in gradient calculation: '
        valstr += '{:.2f} sec ***'.format(tm.time() - start_time)
        self.ostream.print_header(valstr)
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

        # Hessian
        self.hessian = np.zeros((natm, 3, natm, 3))


        if not self.do_four_point:
            for i in range(molecule.number_of_atoms()):
                for d in range(3):
                    coords[i, d] += self.delta_h
                    new_mol = Molecule(labels, coords, units='au')
                    self.scf_drv.compute(new_mol, ao_basis, min_basis)
                    grad_drv.compute(new_mol, ao_basis)
                    #print("Gradient[%d, %d] =" % (i, d), self.grad_drv.print_gradient(molecule))
                    grad_plus = grad_drv.get_gradient()


                    coords[i, d] -= 2.0 * self.delta_h
                    new_mol = Molecule(labels, coords, units='au')
                    self.scf_drv.compute(new_mol, ao_basis, min_basis)
                    grad_drv.compute(new_mol, ao_basis)
                    #print("Gradient[%d, %d] =" % (i, d), self.grad_drv.print_gradient(molecule))
                    grad_minus = grad_drv.get_gradient()

                    coords[i, d] += self.delta_h
                    self.hessian[i, d, :, :] = (grad_plus - grad_minus) / (2.0 * self.delta_h)


        else:
            # Four-point numerical derivative approximation
            # for debugging of analytical Hessian:
            # [ f(x - 2h) - 8 f(x - h) + 8 f(x + h) - f(x + 2h) ] / ( 12h )
            for i in range(molecule.number_of_atoms()):
                for d in range(3):
                    coords[i, d] += self.delta_h
                    new_mol = Molecule(labels, coords, units='au')
                    self.scf_drv.compute(new_mol, ao_basis, min_basis)
                    grad_drv.compute(molecule, ao_basis)
                    grad_plus1 = grad_drv.get_gradient()

                    coords[i, d] += self.delta_h
                    new_mol = Molecule(labels, coords, units='au')
                    self.scf_drv.compute(new_mol, ao_basis, min_basis)
                    grad_drv.compute(molecule, ao_basis)
                    grad_plus2 = grad_drv.get_gradient()

                    coords[i, d] -= 3.0 * self.delta_h
                    new_mol = Molecule(labels, coords, units='au')
                    self.scf_drv.compute(new_mol, ao_basis, min_basis)
                    grad_drv.compute(molecule, ao_basis)
                    grad_minus1 = grad_drv.get_gradient()

                    coords[i, d] -= self.delta_h
                    new_mol = Molecule(labels, coords, units='au')
                    self.scf_drv.compute(new_mol, ao_basis, min_basis)
                    grad_drv.compute(molecule, ao_basis)
                    grad_minus2 = grad_drv.get_gradient()

                    coords[i, d] += 2.0 * self.delta_h
                    # f'(x) ~ [ f(x - 2h) - 8 f(x - h) + 8 f(x + h) - f(x + 2h) ] / ( 12h )
                    self.hessian[i, d] = (grad_minus2 - 8.0 * grad_minus1
                                           + 8.0 * grad_plus1 - grad_plus2) / (12.0 * self.delta_h)

        self.ostream.print_blank()

        self.scf_drv.compute(molecule, ao_basis, min_basis)
        self.scf_drv.ostream.state = scf_ostream_state

    def hess_nuc_contrib(self, molecule):
        """
        Calculates the contribution of the nuclear-nuclear repulsion
        to the analytical nuclear Hessian.

        :param molecule:
            The molecule.

        :return:
            The nuclear contribution to the Hessian.
        """
        # number of atoms
        natm = molecule.number_of_atoms()

        # nuclear repulsion energy contribution to Hessian
        nuc_contrib = np.zeros((natm, 3, natm, 3))

        # atom coordinates (nx3)
        coords = molecule.get_coordinates()

        # atomic charges
        nuclear_charges = molecule.elem_ids_to_numpy()

        # loop over all distinct atom pairs and add energy contribution
        for i in range(natm):
            z_a = nuclear_charges[i]
            r_a = coords[i]
            for j in range(natm):
                if i != j:
                    z_b = nuclear_charges[j]
                    r_b = coords[j]
                    r = np.sqrt(np.dot(r_a - r_b, r_a - r_b))
                    # nondiagonal parts
                    nuc_contrib[i, 0, j, 0] = ( z_a * z_b * ( 1 / r**3 )
                                              - 3*(r_b[0] - r_a[0])**2 / r**5 )
                    nuc_contrib[i, 1, j, 1] = ( z_a * z_b * ( 1 / r**3 )
                                              - 3*(r_b[1] - r_a[1])**2 / r**5 )
                    nuc_contrib[i, 2, j, 2] = ( z_a * z_b * ( 1 / r**3 )
                                              - 3*(r_b[2] - r_a[2])**2 / r**5 )
                    # add the diagonal contribution
                    nuc_contrib[i, 0, i, 0] += ( -z_a * z_b * ( 1 / r**3 )
                                              + 3*(r_b[0] - r_a[0])**2 / r**5 )
                    nuc_contrib[i, 1, i, 1] += ( -z_a * z_b * ( 1 / r**3 )
                                              + 3*(r_b[1] - r_a[1])**2 / r**5 )
                    nuc_contrib[i, 2, i, 2] += ( z_a * z_b * ( 1 / r**3 )
                                              + 3*(r_b[2] - r_a[2])**2 / r**5 )

        return nuc_contrib


    def get_hessian(self):
        """
        Gets the Hessian.

        :return:
            The Hessian.
        """

        return self.hessian

    def print_geometry(self, molecule):
        """
        Prints the geometry.

        :param molecule:
            The molecule.
        """

        self.ostream.print_block(molecule.get_string())

