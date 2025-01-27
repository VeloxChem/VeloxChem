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

from .veloxchemlib import mpi_master
from .molecule import Molecule
from .hessiandriver import HessianDriver
from .errorhandler import assert_msg_critical

class TddftHessianDriver(HessianDriver):
    """
    Implements the TDDFT Hessian driver which is used
    to determine the Hessian for a specific excited state.

    :param scf_drv:
        The SCF driver.
    :param rsp_drv:
        The response driver.
    :param tddft_grad_drv:
        The TDDFT or TDHF gradient driver 

    Instance variables
        - hessian: The Hessian in Hartree per Bohr**2.
    """

    def __init__(self, scf_drv, rsp_drv, tddft_grad_drv):
        """
        Initializes the TDDFT Hessian driver.
        """
        super().__init__(scf_drv.comm, scf_drv.ostream)
        if scf_drv._dft:
            self.flag = "TDDFT Hessian Driver"
        else:
            self.flag = "TDHF Hessian Driver"
        self.scf_drv = scf_drv
        self.rsp_drv = rsp_drv
        self.tddft_grad_drv = tddft_grad_drv
        self.do_print_hessian = False


    def update_settings(self, method_dict, hessian_dict=None):
        """
        Updates settings in TddftHessianDriver.

        :param method_dict:
            The input dictionary of method settings group.
        :param rsp_dict:
            The input dictionary of response settings.
        :param hessian_dict:
            The input dictionary of Hessian settings group.
        :param orbrsp_dict:
            The input dictionary of orbital response settings.
        """
        if hessian_dict is None:
            hessian_dict = {}

        super().update_settings(method_dict, hessian_dict)

    def compute(self, molecule, basis):
        """
        Performs the computation of the TDDFT Hessian.

        :param molecule:
            The molecule.
        :param basis:
            The AO basis set.
        """

        if self.rank == mpi_master():
            self.print_header()

        start_time = tm.time()

        self.compute_numerical_with_analytical_gradient(molecule, basis)

        if self.rank == mpi_master():
            if self.do_print_hessian:
                self.print_geometry(molecule)
                self.ostream.print_blank()
                title = "Numerical TDDFT Hessian based on analytical gradient (Hartree/Bohr**2)"
                self.print_hessian(molecule, title=title)

            valstr = '*** Time spent in Hessian calculation: '
            valstr += '{:.2f} sec ***'.format(tm.time() - start_time)
            self.ostream.print_header(valstr)
            self.ostream.print_blank()
            self.ostream.print_blank()
            self.ostream.flush()

    
    def compute_numerical_with_analytical_gradient(self, molecule, basis):
        """
        Computes the numerical Hessian based on the analytical TDDFT gradient.

        :param molecule:
            The molecule.
        :param basis:
            The AO basis set.
        """
        # atom ids
        labels = molecule.get_labels()

        # atom coordinates (nx3)
        coords = molecule.get_coordinates_in_bohr()
        atom_basis_labels = molecule.get_atom_basis_labels()

        # charge and spin multiplicity
        charge = molecule.get_charge()
        multiplicity = molecule.get_multiplicity()

        # number of atoms
        natm = molecule.number_of_atoms()

        # numerical Hessian
        hessian = np.zeros((natm, 3, natm, 3))

        scf_drv = self.scf_drv
        rsp_drv = self.rsp_drv
        tddft_grad_drv = self.tddft_grad_drv

        # mute ostreams
        scf_drv.ostream.mute()
        rsp_drv.ostream.mute()
        tddft_grad_drv.ostream.mute()

        for i in range(molecule.number_of_atoms()):
            for d in range(3):
                coords[i, d] += self.delta_h
                new_mol = Molecule(labels, coords, 'au', atom_basis_labels)
                new_mol.set_charge(charge)
                new_mol.set_multiplicity(multiplicity)

                grad_plus = self.compute_gradient(new_mol, basis, scf_drv, rsp_drv, tddft_grad_drv)

                coords[i, d] -= 2.0 * self.delta_h
                new_mol = Molecule(labels, coords, 'au', atom_basis_labels)
                new_mol.set_charge(charge)
                new_mol.set_multiplicity(multiplicity)

                grad_minus = self.compute_gradient(new_mol, basis, scf_drv, rsp_drv, tddft_grad_drv)

                if self.do_four_point:
                    coords[i, d] -= self.delta_h
                    new_mol = Molecule(labels, coords, 'au', atom_basis_labels)
                    new_mol.set_charge(charge)
                    new_mol.set_multiplicity(multiplicity)

                    grad_minus2 = self.compute_gradient(new_mol, basis, scf_drv, rsp_drv, tddft_grad_drv)

                    coords[i, d] += 4.0 * self.delta_h
                    new_mol = Molecule(labels, coords, 'au', atom_basis_labels)
                    new_mol.set_charge(charge)
                    new_mol.set_multiplicity(multiplicity)

                    grad_plus2 = self.compute_gradient(new_mol, basis, scf_drv, rsp_drv, tddft_grad_drv)

                    coords[i, d] -= 2.0 * self.delta_h
                    hessian[i, d] = (
                        (grad_minus2 - 8 * grad_minus + 8 * grad_plus - grad_plus2) /
                        (12.0 * self.delta_h))
                else:
                    coords[i, d] += self.delta_h
                    hessian[i, d] = ((grad_plus - grad_minus) /
                                           (2.0 * self.delta_h))

        # restore energy and gradient
        self.compute_gradient(molecule, basis, scf_drv, rsp_drv, tddft_grad_drv)

        # unmute ostreams
        scf_drv.ostream.unmute()
        rsp_drv.ostream.unmute()
        tddft_grad_drv.ostream.unmute()

        # save Hessian in the usual shape
        self.hessian = hessian.reshape((natm*3, natm*3))


    def compute_gradient(self, molecule, basis, scf_drv, rsp_drv, tddft_grad_drv):
        """
        Computes the TDDFT gradient at the current molecular geometry.

        :param molecule:
            The molecule.
        :param basis:
            The AO basis set.
        """
        scf_drv.restart = False
        scf_results = scf_drv.compute(molecule, basis)
        assert_msg_critical(scf_drv.is_converged,
                            'TddftHessianDriver: SCF did not converge')

        rsp_drv.restart = False
        rsp_results = rsp_drv.compute(molecule, basis, scf_results)
        assert_msg_critical(rsp_drv.is_converged,
                            'TddftHessianDriver: response did not converge')

        tddft_grad_drv.compute_analytical(molecule, basis, rsp_results)

        if self.rank == mpi_master():
            # Multiple excited states can be computed simultaneously.
            # For the numerical Hessian, take the first excited state in the list
            gradient = tddft_grad_drv.gradient[0]
        else:
            gradient = None

        return gradient

    def print_header(self):
        """
        Prints Hessian calculation setup details to output stream.
        """

        str_width = 70

        self.ostream.print_blank()
        self.ostream.print_header(self.flag)
        self.ostream.print_header((len(self.flag) + 2) * '=')
        self.ostream.flush()

        cur_str = 'Hessian Type                    : '
        cur_str += 'Numerical using analytical gradient'
        cur_str2 = 'Numerical Method                : '
        if self.do_four_point:
            cur_str2 += 'Five-Point Stencil'
        else:
            cur_str2 += 'Symmetric Difference Quotient'
        cur_str3 = 'Finite Difference Step Size     : '
        cur_str3 += str(self.delta_h) + ' a.u.'
        if self.tddft_grad_drv.state_deriv_index is None:
            s = 1
        else:
            s = self.tddft_grad_drv.state_deriv_index[0]

        cur_str4 = "Exited State                    : %d" % s

        self.ostream.print_blank()
        self.ostream.print_header(cur_str.ljust(str_width))
        self.ostream.print_header(cur_str2.ljust(str_width))
        self.ostream.print_header(cur_str3.ljust(str_width))
        self.ostream.print_header(cur_str4.ljust(str_width))

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
