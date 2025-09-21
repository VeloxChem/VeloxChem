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

from .veloxchemlib import mpi_master
from .outputstream import OutputStream
from .dftutils import get_default_grid_level
from .inputparser import parse_input
from .sanitychecks import dft_sanity_check
from .molecule import Molecule

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

        self.potfile = None

        # Timing and profiling
        self.timing = False
        self.profiling = False
        self.memory_profiling = False
        self.memory_tracing = False

        self._input_keywords = {
            'hessian': {
                'numerical': ('bool', 'do numerical hessian'),
                'do_four_point': ('bool', 'do four-point numerical integration'),
                'delta_h': ('float', 'step size for numerical integration'),
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
        nuclear_charges = molecule.get_element_ids()

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
        self.ostream.flush()

    def print_hessian(self, molecule, title=None):
        """
        Prints the Hessian.

        :param molecule:
            The molecule.
        :param title:
            The title.
        """

        # atom labels
        labels = molecule.get_labels()

        # TODO: move to child classes
        if title is None:
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

    def compute_numerical(self, molecule, *args):
        """
        Computes the numerical Hessian based on the analytical gradient.

        :param molecule:
            The molecule.
        :param args:
            The same arguments as the compute function.
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

        # numerical dipole moment gradient
        if self.do_dipole_gradient:
            dipole_gradient = np.zeros((natm, 3, 3))
        else:
            dipole_gradient = None

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

                grad_plus = self.compute_gradient(new_mol, *args)

                if self.do_dipole_gradient:
                    dipmom_plus = (
                         self.compute_electric_dipole_moment(new_mol, *args))

                coords[i, d] -= 2.0 * self.delta_h
                new_mol = Molecule(labels, coords, 'au', atom_basis_labels)
                new_mol.set_charge(charge)
                new_mol.set_multiplicity(multiplicity)

                grad_minus = self.compute_gradient(new_mol, *args)

                if self.do_dipole_gradient:
                    dipmom_minus = (
                         self.compute_electric_dipole_moment(new_mol, *args))


                if self.do_four_point:
                    coords[i, d] -= self.delta_h
                    new_mol = Molecule(labels, coords, 'au', atom_basis_labels)
                    new_mol.set_charge(charge)
                    new_mol.set_multiplicity(multiplicity)

                    grad_minus2 = self.compute_gradient(new_mol, *args)

                    if self.do_dipole_gradient:
                        dipmom_minus2 = (
                           self.compute_electric_dipole_moment(new_mol, *args))

                    coords[i, d] += 4.0 * self.delta_h
                    new_mol = Molecule(labels, coords, 'au', atom_basis_labels)
                    new_mol.set_charge(charge)
                    new_mol.set_multiplicity(multiplicity)

                    grad_plus2 = self.compute_gradient(new_mol, *args)

                    if self.do_dipole_gradient:
                        dipmom_plus2 = (
                           self.compute_electric_dipole_moment(new_mol, *args))

                    coords[i, d] -= 2.0 * self.delta_h
                    hessian[i, d] = ((grad_minus2 - 8 * grad_minus +
                                      8 * grad_plus - grad_plus2) /
                                     (12.0 * self.delta_h))
                    if self.do_dipole_gradient:
                        dipole_gradient[i, d] = (
                            (dipmom_minus2 - 8 * dipmom_minus +
                             8 * dipmom_plus - dipmom_plus2) /
                            (12.0 * self.delta_h))
                else:
                    coords[i, d] += self.delta_h
                    hessian[i, d] = ((grad_plus - grad_minus) /
                                     (2.0 * self.delta_h))
                    if self.do_dipole_gradient:
                        dipole_gradient[i, d] = ((dipmom_plus - dipmom_minus) /
                                                 (2.0 * self.delta_h))

        # save energy for thermodynamics
        # and restore scf_tensors to results for the original geometry.
        self.elec_energy = self.compute_energy(molecule, *args)

        # save Hessian in the usual shape
        self.hessian = hessian.reshape((natm * 3, natm * 3))

        if self.do_dipole_gradient:
            # save the dipole moment gradient in the expected shape
            self.dipole_gradient = dipole_gradient.transpose(2, 0, 1).reshape(
                3, natm * 3)

        # unmute ostream
        # TODO: how should the ostream be handled properly?
        self.ostream.unmute()


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

    def compute_gradient(self, molecule, *args):
        """
        Computes the gradient at current geometry.

        :param molecule:
            The molecule.
        :param args:
            The same arguments as the "compute" function.

        :return:
            The gradient.
        """

        return

    def compute_electric_dipole_moment(self, molecule, *args):
        """
        Computes the electric dipole moment at current geometry.

        :param molecule:
            The molecule.
        :param args:
            The same arguments as the "compute" function.

        :return:
            The electric dipole moment.
        """

        return

