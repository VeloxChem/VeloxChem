#
#                                   VELOXCHEM
#              ----------------------------------------------------
#                          An Electronic Structure Code
#
#  SPDX-License-Identifier: BSD-3-Clause
#
#  Copyright 2018-2025 VeloxChem developers
#
#  Redistribution and use in source and binary forms, with or without modification,
#  are permitted provided that the following conditions are met:
#
#  1. Redistributions of source code must retain the above copyright notice, this
#     list of conditions and the following disclaimer.
#  2. Redistributions in binary form must reproduce the above copyright notice,
#     this list of conditions and the following disclaimer in the documentation
#     and/or other materials provided with the distribution.
#  3. Neither the name of the copyright holder nor the names of its contributors
#     may be used to endorse or promote products derived from this software without
#     specific prior written permission.
#
#  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
#  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
#  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
#  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
#  FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
#  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
#  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
#  HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
#  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
#  OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

from mpi4py import MPI
import numpy as np
import sys

from .veloxchemlib import XCMolecularGradient
from .veloxchemlib import mpi_master
from .griddriver import GridDriver
from .outputstream import OutputStream
from .molecule import Molecule
from .dftutils import get_default_grid_level
from .inputparser import parse_input
from .sanitychecks import dft_sanity_check


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
        self.delta_h = 0.001

        self._dft = False
        self.grid_level = None
        self.xcfun = None

        self.potfile = None

        self.checkpoint_file = None

        # verbosity of output (1-3)
        self.print_level = 1

        self._input_keywords = {
            'gradient': {
                'numerical': ('bool', 'do numerical integration'),
                'do_four_point':
                    ('bool', 'do four-point numerical integration'),
                'delta_h': ('float', 'the displacement for finite difference'),
            },
            'method_settings': {
                'xcfun': ('str_upper', 'exchange-correlation functional'),
                'grid_level': ('int', 'accuracy level of DFT grid'),
                'dispersion': ('bool', 'do dft-d4 dispersion correction'),
            }
        }

    def update_settings(self, grad_dict, method_dict):
        """
        Updates settings in GradientDriver.

        :param grad_dict:
            The input dictionary of gradient settings group.
        :param method_dict:
            The input dicitonary of method settings group.
        """

        grad_keywords = {
            key: val[0] for key, val in self._input_keywords['gradient'].items()
        }

        parse_input(self, grad_keywords, grad_dict)

        # sanity check
        if (self.do_four_point and not self.numerical):
            self.numerical = True

        method_keywords = {
            key: val[0]
            for key, val in self._input_keywords['method_settings'].items()
        }

        parse_input(self, method_keywords, method_dict)

        dft_sanity_check(self, 'update_settings')

    def compute(self, molecule, *args):
        """
        Performs calculation of nuclear gradient.

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

        # atom ids
        labels = molecule.get_labels()

        # atom coordinates (nx3)
        coords = molecule.get_coordinates_in_bohr()
        atom_basis_labels = molecule.get_atom_basis_labels()

        # charge and spin multiplicity
        charge = molecule.get_charge()
        multiplicity = molecule.get_multiplicity()

        # numerical gradient
        self.gradient = np.zeros((molecule.number_of_atoms(), 3))

        natoms = molecule.number_of_atoms()

        for i in range(natoms):

            self.ostream.unmute()
            self.ostream.print_info(f'Processing atom {i + 1}/{natoms}...')
            self.ostream.flush()
            self.ostream.mute()

            for d in range(3):
                coords[i, d] += self.delta_h
                new_mol = Molecule(labels, coords, 'au', atom_basis_labels)
                new_mol.set_charge(charge)
                new_mol.set_multiplicity(multiplicity)
                e_plus = self.compute_energy(new_mol, *args)

                coords[i, d] -= 2.0 * self.delta_h
                new_mol = Molecule(labels, coords, 'au', atom_basis_labels)
                new_mol.set_charge(charge)
                new_mol.set_multiplicity(multiplicity)
                e_minus = self.compute_energy(new_mol, *args)

                if self.do_four_point:
                    coords[i, d] -= self.delta_h
                    new_mol = Molecule(labels, coords, 'au', atom_basis_labels)
                    new_mol.set_charge(charge)
                    new_mol.set_multiplicity(multiplicity)
                    e_minus2 = self.compute_energy(new_mol, *args)

                    coords[i, d] += 4.0 * self.delta_h
                    new_mol = Molecule(labels, coords, 'au', atom_basis_labels)
                    new_mol.set_charge(charge)
                    new_mol.set_multiplicity(multiplicity)
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
        Calculates the vxc = (d exc / d rho) exchange-correlation contribution
        to the gradient.

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

        # TODO: take mol_grid from arguments

        grid_drv = GridDriver(self.comm)
        grid_level = (get_default_grid_level(self.xcfun)
                      if self.grid_level is None else self.grid_level)
        grid_drv.set_level(grid_level)
        mol_grid = grid_drv.generate(molecule)

        xc_molgrad_drv = XCMolecularGradient()
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

        # TODO: take mol_grid from arguments

        grid_drv = GridDriver(self.comm)
        grid_level = (get_default_grid_level(self.xcfun)
                      if self.grid_level is None else self.grid_level)
        grid_drv.set_level(grid_level)
        mol_grid = grid_drv.generate(molecule)

        xc_molgrad_drv = XCMolecularGradient()
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

        # TODO: take mol_grid from arguments

        grid_drv = GridDriver(self.comm)
        grid_level = (get_default_grid_level(self.xcfun)
                      if self.grid_level is None else self.grid_level)
        grid_drv.set_level(grid_level)
        mol_grid = grid_drv.generate(molecule)

        xc_molgrad_drv = XCMolecularGradient()
        vxc3_contrib = xc_molgrad_drv.integrate_kxc_gradient(
            molecule, ao_basis, rhow_den_1, rhow_den_2, gs_density, mol_grid,
            xcfun_label)
        vxc3_contrib = self.comm.reduce(vxc3_contrib, root=mpi_master())

        return vxc3_contrib

    def grad_tddft_xc_contrib(self, molecule, ao_basis, rhow_den, xmy_den,
                              gs_density, xcfun_label):
        """
        Calculates exchange-correlation contribution to TDDFT gradient.

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
            The exchange-correlation contribution to TDDFT gradient.
        """

        # TODO: take mol_grid from arguments

        grid_drv = GridDriver(self.comm)
        grid_level = (get_default_grid_level(self.xcfun)
                      if self.grid_level is None else self.grid_level)
        grid_drv.set_level(grid_level)
        mol_grid = grid_drv.generate(molecule)

        xcgrad_drv = XCMolecularGradient()
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
        coords = molecule.get_coordinates_in_bohr()

        # atomic charges
        nuclear_charges = molecule.get_element_ids()

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

        return self.gradient.copy()

    def print_geometry(self, molecule):
        """
        Prints the geometry.

        :param molecule:
            The molecule.
        """

        self.ostream.print_block(molecule.get_string())
        self.ostream.flush()

    def print_gradient(self, molecule, state_deriv_index=None):
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

        if (state_deriv_index is None):
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
        else:
            gradient_shape = self.gradient.shape
            for index, s in enumerate(state_deriv_index):
                self.ostream.print_blank()
                state = 'Excited State %d' % s
                self.ostream.print_header(state)
                self.ostream.print_blank()
                valstr = '  Atom '
                valstr += '{:>20s}  '.format('Gradient X')
                valstr += '{:>20s}  '.format('Gradient Y')
                valstr += '{:>20s}  '.format('Gradient Z')
                self.ostream.print_header(valstr)
                self.ostream.print_blank()
                for i in range(molecule.number_of_atoms()):
                    valstr = '  {:<4s}'.format(labels[i])
                    if len(gradient_shape) == 2:
                        for d in range(3):
                            valstr += '{:22.12f}'.format(self.gradient[i, d])
                    else:
                        for d in range(3):
                            valstr += '{:22.12f}'.format(self.gradient[index, i,
                                                                       d])
                    self.ostream.print_header(valstr)
        self.ostream.print_blank()
        self.ostream.flush()

    def print_header(self, state_deriv_index=None):
        """
        Prints gradient calculation setup details to output stream.
        """

        str_width = 70

        self.ostream.print_blank()
        self.ostream.print_header(self.flag)
        self.ostream.print_header((len(self.flag) + 2) * '=')
        self.ostream.flush()

        cur_str = 'Gradient Type                   : '

        if self.numerical:
            cur_str += 'Numerical'
            cur_str2 = 'Numerical Method                : '
            if self.do_four_point:
                cur_str2 += 'Five-Point Stencil'
            else:
                cur_str2 += 'Symmetric Difference Quotient'
            cur_str3 = 'Finite Difference Step Size     : '
            cur_str3 += str(self.delta_h) + ' a.u.'
        else:
            cur_str += 'Analytical'

        self.ostream.print_blank()
        self.ostream.print_header(cur_str.ljust(str_width))
        if self.numerical:
            self.ostream.print_header(cur_str2.ljust(str_width))
            self.ostream.print_header(cur_str3.ljust(str_width))

        if self._dft:
            cur_str = 'Exchange-Correlation Functional : '
            cur_str += self.xcfun.get_func_label().upper()
            self.ostream.print_header(cur_str.ljust(str_width))
            grid_level = (get_default_grid_level(self.xcfun)
                          if self.grid_level is None else self.grid_level)
            cur_str = 'Molecular Grid Level            : ' + str(grid_level)
            self.ostream.print_header(cur_str.ljust(str_width))

        if state_deriv_index is not None:
            if isinstance(state_deriv_index, int):
                cur_str = ('Excited State of Interest       : ' +
                           str(state_deriv_index + 1))
                self.ostream.print_header(cur_str.ljust(str_width))
            elif isinstance(state_deriv_index, (list, tuple)):
                state_deriv_index_str = [str(s) for s in state_deriv_index]
                states_txt = ", ".join(state_deriv_index_str)
                cur_str = 'Excited States of Interest      : ' + states_txt
                self.ostream.print_header(cur_str.ljust(str_width))

        self.ostream.print_blank()
        self.ostream.flush()
