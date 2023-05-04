#
#                           VELOXCHEM 1.0-RC3
#         ----------------------------------------------------
#                     An Electronic Structure Code
#
#  Copyright © 2018-2022 by VeloxChem developers. All rights reserved.
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
import sys

from .veloxchemlib import GridDriver, XCMolecularGradient
from .veloxchemlib import mpi_master
from .veloxchemlib import parse_xc_func
from .outputstream import OutputStream
from .molecule import Molecule
from .errorhandler import assert_msg_critical


class GradientDriver:
    """
    Implements gradient driver.

    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.

    Instance variables
        - gradient: The gradient.
        - flag: The type of gradient driver.
        - numerical: Perform numerical gradient calculation.
        - do_four_point: Perform numerical gradient calculation using the
          five-point stencil method?
        - dft: The flag for running DFT.
        - grid_level: The accuracy level of DFT grid.
        - xcfun: The XC functional.
    """

    def __init__(self, comm=None, ostream=None):
        """
        Initializes gradient driver.
        """

        if comm is None:
            comm = MPI.COMM_WORLD

        if ostream is None:
            if comm.Get_rank() == mpi_master():
                ostream = OutputStream(sys.stdout)
            else:
                ostream = OutputStream(None)

        self.comm = comm
        self.rank = self.comm.Get_rank()
        self.nodes = self.comm.Get_size()

        self.ostream = ostream

        self.gradient = None
        self.flag = None

        self.numerical = False
        self.do_four_point = False

        self.dft = False
        self.grid_level = 4
        self.xcfun = None

    def update_settings(self, grad_dict, method_dict):
        """
        Updates settings in GradientDriver.

        :param grad_dict:
            The input dictionary of gradient settings group.
        :param method_dict:
            The input dicitonary of method settings group.
        """

        if 'do_four_point' in grad_dict:
            key = grad_dict['do_four_point'].lower()
            self.do_four_point = (key in ['yes', 'y'])
            if self.do_four_point:
                self.numerical = True

        if 'numerical' in grad_dict:
            key = grad_dict['numerical'].lower()
            self.numerical = (key in ['yes', 'y'])

        # TODO: use parse_input and _dft_sanity_check

        if 'dft' in method_dict:
            key = method_dict['dft'].lower()
            self.dft = (key in ['yes', 'y'])
        if 'grid_level' in method_dict:
            self.grid_level = int(method_dict['grid_level'])
        if 'xcfun' in method_dict:
            if 'dft' not in method_dict:
                self.dft = True
            self.xcfun = parse_xc_func(method_dict['xcfun'].upper())
            assert_msg_critical(not self.xcfun.is_undefined(),
                                'Gradient driver: Undefined XC functional')

        if 'delta_h' in grad_dict:
            self.delta_h = float(grad_dict['delta_h'])

    def compute(self, molecule, *args):
        """
        Performs calculation of numerical gradient.

        :param molecule:
            The molecule.
        :param args:
            The arguments.
        """

        return

    def compute_numerical(self, molecule, *args):
        """
        Performs calculation of numerical gradient.

        :param molecule:
            The molecule.
        :param args:
            The same arguments as the "compute" function.
        """

        # atom labels
        labels = molecule.get_labels()

        # atom coordinates (nx3)
        coords = molecule.get_coordinates()

        # numerical gradient
        self.gradient = np.zeros((molecule.number_of_atoms(), 3))

        for i in range(molecule.number_of_atoms()):
            for d in range(3):
                coords[i, d] += self.delta_h
                new_mol = Molecule(labels, coords, units='au')
                e_plus = self.compute_energy(new_mol, *args)

                coords[i, d] -= 2.0 * self.delta_h
                new_mol = Molecule(labels, coords, units='au')
                e_minus = self.compute_energy(new_mol, *args)

                if self.do_four_point:
                    coords[i, d] -= self.delta_h
                    new_mol = Molecule(labels, coords, units='au')
                    e_minus2 = self.compute_energy(new_mol, *args)

                    coords[i, d] += 4.0 * self.delta_h
                    new_mol = Molecule(labels, coords, units='au')
                    e_plus2 = self.compute_energy(new_mol, *args)

                    coords[i, d] -= 2.0 * self.delta_h
                    self.gradient[i, d] = (
                        (e_minus2 - 8 * e_minus + 8 * e_plus - e_plus2) /
                        (12.0 * self.delta_h))
                else:
                    coords[i, d] += self.delta_h
                    self.gradient[i, d] = ((e_plus - e_minus) /
                                           (2.0 * self.delta_h))

        # restore energy driver
        self.compute_energy(molecule, *args)

    def compute_energy(self, molecule, *args):
        """
        Computes the energy at current geometry.

        :param molecule:
            The molecule.
        :param args:
            The same arguments as the "compute" function.

        :return:
            The energy.
        """

        return

    def grad_vxc_contrib(self, molecule, ao_basis, rhow_density, gs_density,
                         xcfun_label):
        """
        Calculates the vxc exchange-correlation contribution to the gradient.

        :param molecule:
            The molecule.
        :param ao_basis:
            The AO basis set.
        :param rhow_density:
            The perturbed density (to be contracted with GTO gradient).
        :param gs_density:
            The ground state density (to generate functional derivative).
        :param xcfun_label:
            The label of the xc functional.

        :return:
            The vxc exchange-correlation contribution to the gradient.
        """

        grid_drv = GridDriver(self.comm)
        grid_drv.set_level(self.grid_level)
        mol_grid = grid_drv.generate(molecule)

        xc_molgrad_drv = XCMolecularGradient(self.comm)
        vxc_contrib = xc_molgrad_drv.integrate_vxc_gradient(
            molecule, ao_basis, rhow_density, gs_density, mol_grid, xcfun_label)
        vxc_contrib = self.comm.reduce(vxc_contrib, root=mpi_master())

        return vxc_contrib

    def grad_vxc2_contrib(self, molecule, ao_basis, rhow_den_1, rhow_den_2,
                          gs_density, xcfun_label):
        """
        Calculates the 2nd-order exchange-correlation contribution to the gradient.

        :param molecule:
            The molecule.
        :param ao_basis:
            The AO basis set.
        :param rhow_den_1:
            The perturbed density.
        :param rhow_den_2:
            The perturbed density (to be contracted with GTO gradient).
        :param gs_density:
            The ground state density (to generate functional derivative).
        :param xcfun_label:
            The label of the xc functional.

        :return:
            The 2nd-order exchange-correlation contribution to the gradient.
        """

        grid_drv = GridDriver(self.comm)
        grid_drv.set_level(self.grid_level)
        mol_grid = grid_drv.generate(molecule)

        xc_molgrad_drv = XCMolecularGradient(self.comm)
        vxc2_contrib = xc_molgrad_drv.integrate_fxc_gradient(
            molecule, ao_basis, rhow_den_1, rhow_den_2, gs_density, mol_grid,
            xcfun_label)
        vxc2_contrib = self.comm.reduce(vxc2_contrib, root=mpi_master())

        return vxc2_contrib

    def grad_vxc3_contrib(self, molecule, ao_basis, rhow_den_1, rhow_den_2,
                          gs_density, xcfun_label):
        """
        Calculates the 3rd-order exchange-correlation contribution to the
        gradient.

        :param molecule:
            The molecule.
        :param ao_basis:
            The AO basis set.
        :param rhow_den_1:
            The perturbed density.
        :param rhow_den_2:
            The perturbed density.
        :param gs_density:
            The ground state density (to generate functional derivative and
            to be contracted with GTO gradient).
        :param xcfun_label:
            The label of the xc functional.

        :return:
            The 3rd-order exchange-correlation contribution to the gradient.
        """

        grid_drv = GridDriver(self.comm)
        grid_drv.set_level(self.grid_level)
        mol_grid = grid_drv.generate(molecule)

        xc_molgrad_drv = XCMolecularGradient(self.comm)
        vxc3_contrib = xc_molgrad_drv.integrate_kxc_gradient(
            molecule, ao_basis, rhow_den_1, rhow_den_2, gs_density, mol_grid,
            xcfun_label)
        vxc3_contrib = self.comm.reduce(vxc3_contrib, root=mpi_master())

        return vxc3_contrib

    def grad_tddft_xc_contrib(self, molecule, ao_basis, rhow_den, xmy_den,
                              gs_density, xcfun_label):
        """
        Calculates exchange-correlation contribution to tddft gradient.

        :param molecule:
            The molecule.
        :param ao_basis:
            The AO basis set.
        :param rhow_den:
            The perturbed density.
        :param xmy_den:
            The X-Y density.
        :param gs_density:
            The ground state density.
        :param xcfun_label:
            The label of the xc functional.

        :return:
            The exchange-correlation contribution to tddft gradient.
        """

        grid_drv = GridDriver(self.comm)
        grid_drv.set_level(self.grid_level)
        mol_grid = grid_drv.generate(molecule)

        xcgrad_drv = XCMolecularGradient(self.comm)
        tddft_xcgrad = xcgrad_drv.integrate_vxc_gradient(
            molecule, ao_basis, rhow_den, gs_density, mol_grid, xcfun_label)
        tddft_xcgrad += xcgrad_drv.integrate_fxc_gradient(
            molecule, ao_basis, rhow_den, gs_density, gs_density, mol_grid,
            xcfun_label)
        tddft_xcgrad += xcgrad_drv.integrate_fxc_gradient(
            molecule, ao_basis, xmy_den, xmy_den, gs_density, mol_grid,
            xcfun_label)
        tddft_xcgrad += xcgrad_drv.integrate_kxc_gradient(
            molecule, ao_basis, xmy_den, xmy_den, gs_density, mol_grid,
            xcfun_label)
        tddft_xcgrad = self.comm.reduce(tddft_xcgrad, root=mpi_master())

        return tddft_xcgrad

    def grad_nuc_contrib(self, molecule):
        """
        Calculates the contribution of the nuclear-nuclear repulsion
        to the analytical nuclear gradient.

        :param molecule:
            The molecule.

        :return:
            The nuclear contribution to the gradient.
        """

        # number of atoms
        natm = molecule.number_of_atoms()

        # nuclear repulsion energy contribution to nuclear gradient
        nuc_contrib = np.zeros((natm, 3))

        # atom coordinates (nx3)
        coords = molecule.get_coordinates()

        # atomic charges
        nuclear_charges = molecule.elem_ids_to_numpy()

        # loop over all distinct atom pairs and add energy contribution
        for i in range(natm):
            z_a = nuclear_charges[i]
            r_a = coords[i]
            for j in range(i + 1, natm):
                z_b = nuclear_charges[j]
                r_b = coords[j]
                r = np.sqrt(np.dot(r_a - r_b, r_a - r_b))
                f_ij = z_a * z_b * (r_b - r_a) / r**3
                nuc_contrib[i] += f_ij
                nuc_contrib[j] -= f_ij

        return nuc_contrib

    def get_gradient(self):
        """
        Gets the gradient.

        :return:
            The gradient.
        """

        return self.gradient

    def print_geometry(self, molecule):
        """
        Prints the geometry.

        :param molecule:
            The molecule.
        """

        self.ostream.print_block(molecule.get_string())

    def print_gradient(self, molecule):
        """
        Prints the gradient.

        :param molecule:
            The molecule.
        """

        labels = molecule.get_labels()

        if self.numerical:
            title = 'Numerical '
        else:
            title = 'Analytical '

        title += 'Gradient (Hartree/Bohr)'
        self.ostream.print_header(title)
        self.ostream.print_header('-' * (len(title) + 2))
        self.ostream.print_blank()

        valstr = '  Atom '
        valstr += '{:>20s}  '.format('Gradient X')
        valstr += '{:>20s}  '.format('Gradient Y')
        valstr += '{:>20s}  '.format('Gradient Z')
        self.ostream.print_header(valstr)
        self.ostream.print_blank()

        for i in range(molecule.number_of_atoms()):
            valstr = '  {:<4s}'.format(labels[i])
            for d in range(3):
                valstr += '{:22.12f}'.format(self.gradient[i, d])
            self.ostream.print_header(valstr)

        self.ostream.print_blank()
        self.ostream.flush()

    def print_header(self, state_deriv_index=None):
        """
        Prints gradient calculation setup details to output stream.
        """

        str_width = 60

        self.ostream.print_blank()
        self.ostream.print_header(self.flag)
        self.ostream.print_header((len(self.flag) + 2) * '=')
        self.ostream.flush()

        cur_str = 'Gradient Type               : '

        if self.numerical:
            cur_str += 'Numerical'
            cur_str2 = 'Numerical Method            : '
            if self.do_four_point:
                cur_str2 += 'Five-Point Stencil'
            else:
                cur_str2 += 'Symmetric Difference Quotient'
            cur_str3 = 'Finite Difference Step Size : '
            cur_str3 += str(self.delta_h) + ' a.u.'
        else:
            cur_str += 'Analytical'

        if self.numerical or state_deriv_index is not None:
            self.ostream.print_blank()
            self.ostream.print_header(cur_str.ljust(str_width))
        if self.numerical:
            self.ostream.print_header(cur_str2.ljust(str_width))
            self.ostream.print_header(cur_str3.ljust(str_width))

        if state_deriv_index is not None:
            cur_str = 'Excited State of Interest   : S' + str(state_deriv_index)
            self.ostream.print_header(cur_str.ljust(str_width))

        self.ostream.print_blank()
        self.ostream.flush()
