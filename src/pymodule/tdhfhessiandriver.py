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
from .tdhfgradientdriver import TdhfGradientDriver
from .scfhessiandriver import ScfHessianDriver
from .outputstream import OutputStream
#from .firstorderprop import FirstOrderProperties

# For PySCF integral derivatives
from .import_from_pyscf import overlap_deriv
from .import_from_pyscf import fock_deriv
from .import_from_pyscf import eri_deriv


class TdhfHessianDriver(ScfHessianDriver):
    """
    Implements the Hessian at the TDHF and CIS levels of theory.

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

        self.flag = 'TDHF Hessian Driver'
        self.hessian = None
        self.scf_drv = scf_drv
        self.delta_h = 0.001
        self.numerical = True #False
        # flag for two-point or four-point approximation
        self.do_four_point = False


        # Flag for numerical derivative of dipole moment
        self.dipole_deriv = False
        self.dipole_gradient = None

    def update_settings(self, method_dict, rsp_dict):
        """
        Updates settings in ScfHessianDriver.

        TODO: Add new gradient settings group?

        :param method_dict:
            The input dicitonary of method settings group.
        :param response_dict:
            The input dicitonary of response settings group.
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

        if 'n_state_deriv' in rsp_dict:
            # user gives '1' for first excited state, but internal index is 0
            self.n_state_deriv = int(rsp_dict['n_state_deriv']) - 1

        self.rsp_dict = dict(rsp_dict)
        self.method_dict = dict(method_dict)

    def compute(self, molecule, ao_basis, rsp_drv, min_basis=None):
        """
        Computes the numerical nuclear Hessian at the RPA or TDA level of theory.

        :param molecule:
            The molecule.
        :param ao_basis:
            The AO basis set.
        :param rsp_drv:
            The RPA or TDA driver.
        :param min_basis:
            The minimal AO basis set.
        """
        self.print_header()

        start_time = tm.time()

        if True: #self.numerical:
            self.compute_numerical(molecule, ao_basis, rsp_drv, min_basis)
        #else:
        #    self.compute_analytical(molecule, ao_basis)

        # print Hessian
        #self.print_geometry(molecule)
        #self.print_gradient(molecule)

        valstr = '*** Time spent in Hessian calculation: '
        valstr += '{:.2f} sec ***'.format(tm.time() - start_time)
        self.ostream.print_header(valstr)
        self.ostream.print_blank()
        self.ostream.flush()

###    def compute_analytical(self, molecule, ao_basis):
###        """
###        Computes the analytical nuclear Hessian.
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

    def compute_numerical(self, molecule, ao_basis, rsp_drv, min_basis=None):
        """
        Performs calculation of numerical Hessian.

        :param molecule:
            The molecule.
        :param ao_basis:
            The AO basis set.
        :param rsp_drv:
            The RPA or TDA driver.
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
        grad_drv = TdhfGradientDriver(self.scf_drv, self.scf_drv.comm, self.scf_drv.ostream)
        grad_drv.update_settings(self.rsp_dict, self.method_dict)

        # Hessian
        hessian = np.zeros((natm, 3, natm, 3))


        if not self.do_four_point:
            for i in range(molecule.number_of_atoms()):
                for d in range(3):
                    coords[i, d] += self.delta_h
                    new_mol = Molecule(labels, coords, units='au')
                    self.scf_drv.compute(new_mol, ao_basis, min_basis)
                    scf_tensors = self.scf_drv.scf_tensors
                    rsp_drv.is_converged = False  # only needed for RPA
                    rsp_results = rsp_drv.compute(new_mol, ao_basis,
                                                       scf_tensors)
                    grad_drv.compute(new_mol, ao_basis, rsp_drv, rsp_results)
                    #print("Gradient[%d, %d] =" % (i, d), self.grad_drv.print_gradient(molecule))
                    grad_plus = grad_drv.get_gradient()


                    coords[i, d] -= 2.0 * self.delta_h
                    new_mol = Molecule(labels, coords, units='au')
                    self.scf_drv.compute(new_mol, ao_basis, min_basis)
                    scf_tensors = self.scf_drv.scf_tensors
                    rsp_drv.is_converged = False  # only needed for RPA
                    rsp_results = rsp_drv.compute(new_mol, ao_basis,
                                                       scf_tensors)
                    grad_drv.compute(new_mol, ao_basis, rsp_drv, rsp_results)
                    #print("Gradient[%d, %d] =" % (i, d), self.grad_drv.print_gradient(molecule))
                    grad_minus = grad_drv.get_gradient()

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
                    scf_tensors = self.scf_drv.scf_tensors
                    rsp_drv.is_converged = False  # only needed for RPA
                    rsp_results = rsp_drv.compute(new_mol, ao_basis,
                                                       scf_tensors)
                    grad_drv.compute(new_mol, ao_basis, rsp_drv, rsp_results)
                    grad_plus1 = grad_drv.get_gradient()

                    coords[i, d] += self.delta_h
                    new_mol = Molecule(labels, coords, units='au')
                    self.scf_drv.compute(new_mol, ao_basis, min_basis)
                    scf_tensors = self.scf_drv.scf_tensors
                    rsp_drv.is_converged = False  # only needed for RPA
                    rsp_results = rsp_drv.compute(new_mol, ao_basis,
                                                       scf_tensors)
                    grad_drv.compute(new_mol, ao_basis, rsp_drv, rsp_results)
                    grad_plus2 = grad_drv.get_gradient()

                    coords[i, d] -= 3.0 * self.delta_h
                    new_mol = Molecule(labels, coords, units='au')
                    self.scf_drv.compute(new_mol, ao_basis, min_basis)
                    scf_tensors = self.scf_drv.scf_tensors
                    rsp_drv.is_converged = False  # only needed for RPA
                    rsp_results = rsp_drv.compute(new_mol, ao_basis,
                                                       scf_tensors)
                    grad_drv.compute(new_mol, ao_basis, rsp_drv, rsp_results)
                    grad_minus1 = grad_drv.get_gradient()

                    coords[i, d] -= self.delta_h
                    new_mol = Molecule(labels, coords, units='au')
                    self.scf_drv.compute(new_mol, ao_basis, min_basis)
                    scf_tensors = self.scf_drv.scf_tensors
                    rsp_drv.is_converged = False  # only needed for RPA
                    rsp_results = rsp_drv.compute(new_mol, ao_basis,
                                                       scf_tensors)
                    grad_drv.compute(new_mol, ao_basis, rsp_drv, rsp_results)
                    grad_minus2 = grad_drv.get_gradient()

                    coords[i, d] += 2.0 * self.delta_h
                    # f'(x) ~ [ f(x - 2h) - 8 f(x - h) + 8 f(x + h) - f(x + 2h) ] / ( 12h )
                    hessian[i, d] = (grad_minus2 - 8.0 * grad_minus1
                                           + 8.0 * grad_plus1 - grad_plus2) / (12.0 * self.delta_h)

        # reshaped Hessian as member variable
        self.hessian = hessian.reshape(3*natm, 3*natm)

        self.ostream.print_blank()

        self.scf_drv.compute(molecule, ao_basis, min_basis)
        self.scf_drv.ostream.state = scf_ostream_state

