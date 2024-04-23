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
from contextlib import redirect_stderr
from io import StringIO
import numpy as np
import sys

from .veloxchemlib import GridDriver, XCMolecularHessian
from .veloxchemlib import (mpi_master, bohr_in_angstrom, avogadro_constant,
                           fine_structure_constant, electron_mass_in_amu,
                           amu_in_kg, speed_of_light_in_vacuum_in_SI)
from .veloxchemlib import parse_xc_func
from .outputstream import OutputStream
from .errorhandler import assert_msg_critical
from .dftutils import get_default_grid_level
from .inputparser import parse_input
from .sanitychecks import dft_sanity_check

with redirect_stderr(StringIO()) as fg_err:
    import geometric


class HessianDriver:
    """
    Implements the Hessian driver.

    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.

    Instance variables
        - hessian: The Hessian in Hartree per Bohr**2.
        - dipole_gradient: The gradient of the dipole moment.
        - flag: The type of Hessian driver.
        - numerical: Perform numerical Hessian calculation.
        - delta_h: Nuclear displacement for finite differences.
        - do_four_point: Perform four-point numerical approximation.
        - do_print_hessian: Flag for printing the Hessian.
        - do_dipole_gradient
    """

    def __init__(self, comm=None, ostream=None):
        """
        Initializes Hessian driver.
        """

        if comm is None:
            comm = MPI.COMM_WORLD

        if ostream is None:
            if comm.Get_rank() == mpi_master():
                ostream = OutputStream(sys.stdout)
            else:
                ostream = OutputStream(None)

        # MPI information
        self.comm = comm
        self.rank = self.comm.Get_rank()
        self.nodes = self.comm.Get_size()

        self.ostream = ostream

        self.hessian = None

        self.do_dipole_gradient = False
        self.dipole_gradient = None
        self.flag = None

        self.numerical = False
        self.delta_h = 0.001

        # flag for two-point or four-point approximation
        self.do_four_point = False

        # flag for printing the Hessian
        self.do_print_hessian = False

        self._dft = False
        self.grid_level = None
        self.xcfun = None

        # Timing and profiling
        self.timing = False
        self.profiling = False
        self.memory_profiling = False
        self.memory_tracing = False

        self._input_keywords = {
            'hessian': {
                'numerical': ('bool', 'do numerical hessian'),
                'do_four_point': ('bool', 'do four-point numerical integration'),
                'numerical_grad': ('bool', 'whether the gradient is numerical'),
                'do_print_hessian': ('bool', 'whether to print the Hessian'),
                'do_dipole_gradient': ('bool', 'whether to compute the dipole gradient'),
                'timing': ('bool', 'whether timing is needed'),
                'profiling': ('bool', 'whether profiling is needed'),
                'memory_profiling': ('bool', 'whether to profile memory'),
                'memory_tracing': ('bool', 'whether to trace memory'),
                },
            'method_settings': {
                'xcfun': ('str_upper', 'exchange-correlation functional'),
                'grid_level': ('int', 'accuracy level of DFT grid'),
                }
            }

    def update_settings(self, method_dict=None, hess_dict=None):
        """
        Updates settings in HessianDriver.

        :param method_dict:
            The input dictionary of method settings group.
        :param hess_dict:
            The input dictionary of Hessian settings group.
        """

        if method_dict is None:
            method_dict = {}
        if hess_dict is None:
            hess_dict = {}

        hess_keywords = {
            key: val[0] for key, val in
            self._input_keywords['hessian'].items()
        }

        parse_input(self, hess_keywords, hess_dict)

        method_keywords = {
            key: val[0]
            for key, val in self._input_keywords['method_settings'].items()
        }

        parse_input(self, method_keywords, method_dict)

        dft_sanity_check(self, 'update_settings')

        self.method_dict = dict(method_dict)
        self.hess_dict = dict(hess_dict)

    def compute(self, molecule, ao_basis=None, min_basis=None):
        """
        Performs calculation of molecular Hessian.

        :param molecule:
            The molecule.
        :param ao_basis:
            The AO basis set.
        :param min_basis:
            The minimal AO basis set.
        """

        return

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
        nuc_contrib = np.zeros((natm, natm, 3, 3))

        # atom coordinates (nx3)
        coords = molecule.get_coordinates_in_bohr()

        # atomic charges
        nuclear_charges = molecule.elem_ids_to_numpy()

        # loop over all distinct atom pairs and add energy contribution
        for i in range(natm):
            z_a = nuclear_charges[i]
            r_a = coords[i]
            # upper triangular part
            for j in range(natm):
                if i != j:
                    z_b = nuclear_charges[j]
                    r_b = coords[j]
                    r = np.sqrt(np.dot(r_a - r_b, r_a - r_b))
                    for k in range(3):
                        for l in range(3):

                            # off-diagonal parts
                            nuc_contrib[i, j, k, l] = -3 * z_a * z_b * (
                                r_b[k] - r_a[k]) * (r_b[l] - r_a[l]) / r**5
                            if k == l:
                                nuc_contrib[i, j, k, l] += z_a * z_b / r**3

                            # add the diagonal contribution
                            nuc_contrib[i, i, k, l] += 3 * z_a * z_b * (
                                r_b[k] - r_a[k]) * (r_b[l] - r_a[l]) / r**5
                            if k == l:
                                nuc_contrib[i, i, k, l] -= z_a * z_b / r**3

        return nuc_contrib

    def hess_xc_contrib(self, molecule, basis, gs_density, xcfun_label):
        """
        Calculates the exchange-correlation contribution to the analytical
        nuclear Hessian.

        :param molecule:
            The molecule.
        :param basis:
            The AO basis set.
        :param gs_density:
            The ground state AO density matrix object.
        :param xcfun_label:
            The name of the exchange-correlation functional.

        :return:
            The exchange-correlation contribution to the Hessian.
        """

        xc_mol_hess = XCMolecularHessian(self.comm)

        grid_drv = GridDriver(self.comm)
        grid_level = (get_default_grid_level(self.xcfun)
                      if self.grid_level is None else self.grid_level)
        grid_drv.set_level(grid_level)
        mol_grid = grid_drv.generate(molecule)

        exc_hessian = xc_mol_hess.integrate_vxc_hessian(molecule, basis,
                                                        gs_density, mol_grid,
                                                        xcfun_label)
        exc_hessian += xc_mol_hess.integrate_fxc_hessian(
            molecule, basis, gs_density, mol_grid, xcfun_label)
        exc_hessian = self.comm.reduce(exc_hessian, root=mpi_master())

        return exc_hessian

    def vxc_fock_grad_xc_contrib(self, molecule, basis, gs_density, xcfun_label,
                                 atom_idx):
        """
        Calculates the exchange-correlation contribution to the analytical
        nuclear gradient of Vxc Fock matrix element w.r.t. a given atom.

        :param molecule:
            The molecule.
        :param basis:
            The AO basis set.
        :param gs_density:
            The ground state AO density matrix object.
        :param xcfun_label:
            The name of the exchange-correlation functional.
        :param atom_idx:
            The index (0-based) of the atom.

        :return:
            The exchange-correlation contribution to the nuclear gradient of
            Vxc Fock matrix element w.r.t. a given atom.
        """

        xc_mol_hess = XCMolecularHessian(self.comm)

        grid_drv = GridDriver(self.comm)
        grid_level = (get_default_grid_level(self.xcfun)
                      if self.grid_level is None else self.grid_level)
        grid_drv.set_level(grid_level)
        mol_grid = grid_drv.generate(molecule)

        vxc_grad_atom = xc_mol_hess.integrate_vxc_fock_gradient(
            molecule, basis, gs_density, mol_grid, xcfun_label, atom_idx)
        vxc_grad_atom = self.comm.reduce(vxc_grad_atom, root=mpi_master())

        return vxc_grad_atom

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

    def print_hessian(self, molecule):
        """
        Prints the Hessian.

        :param molecule:
            The molecule.
        """

        # atom labels
        labels = molecule.get_labels()

        if self.numerical:
            title = 'Numerical '
        else:
            title = 'Analytical '

        title += 'Hessian (Hartree/Bohr**2)'
        self.ostream.print_header(title)
        self.ostream.print_header('-' * (len(title) + 2))
        self.ostream.print_blank()

        natm = molecule.number_of_atoms()

        for k in range(0, natm, 2):

            valstr = '{:15s}'.format('  Coord. ')

            coord_dict = {0: '(x)', 1: '(y)', 2: '(z)'}
            end = k + 2
            if k + 2 > natm:
                end = natm
            for i in range(k, end):
                for di in range(3):
                    valstr += '{:16s}'.format('' + str(i + 1) + ' ' +
                                              labels[i] + coord_dict[di] + '')

            self.ostream.print_line(valstr)
            self.ostream.print_blank()

            for i in range(natm):
                for di in range(3):
                    valstr = '  {:7s}'.format(
                        str(i + 1) + ' ' + labels[i] + coord_dict[di])
                    for j in range(k, end):
                        for dj in range(3):
                            valstr += '{:16.8f}'.format(
                                self.hessian[i * 3 + di, j * 3 + dj])
                    self.ostream.print_line(valstr)
            self.ostream.print_blank()
            self.ostream.print_blank()
        self.ostream.flush()

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

        self.ostream.print_blank()
        self.ostream.flush()

