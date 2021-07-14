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

import numpy as np
import sys
from mpi4py import MPI

from .outputstream import OutputStream

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
    """

    def __init__(self, comm=None, ostream=None):
        """
        Initializes gradient driver.
        """

        if comm is None:
            comm = MPI.COMM_WORLD

        if ostream is None:
            ostream = OutputStream(sys.stdout)

        self.comm = comm
        self.ostream = ostream

        self.gradient = None
        self.flag = None
        self.numerical = False
        # flag for two-point or four-point approximation
        self.do_four_point = False

    def update_settings(self, grad_dict, method_dict):
        """
        Updates settings in GradientDriver.

        :param grad_dict:
            The input dictionary of gradient settings group.
        :param method_dict:
            The input dicitonary of method settings group.
        """

        # Numerical gradient?
        if 'numerical' in grad_dict:
            key = grad_dict['numerical'].lower()
            self.numerical = True if key in ['yes', 'y'] else False

        # TODO: Analytical DFT gradient is not implemented yet
        is_dft = False
        if 'xcfun' in method_dict:
            if method_dict['xcfun'] is not None:
                is_dft = True
                #self.numerical = True
        if 'dft' in method_dict:
            key = method_dict['dft'].lower()
            #self.numerical = True if key in ['yes', 'y'] else False
            is_dft = True if key in ['yes', 'y'] else False

        if is_dft and not self.numerical:
            self.numerical = True
            warn_msg = '*** Warning: Analytical DFT gradient is not yet implemented.'
            self.ostream.print_blank()
            self.ostream.print_header(warn_msg.ljust(56))
            warn_msg = '              Gradient will be calculated numerically instead.'
            self.ostream.print_header(warn_msg.ljust(56))
            self.ostream.flush()

        # step size for finite differences
        if 'delta_h' in grad_dict:
            self.delta_h = float(grad_dict['delta_h'])

        # TODO: if this is True, numerical must also be True
        if 'do_four_point' in grad_dict:
            key = grad_dict['do_four_point'].lower()
            self.do_four_point = True if key in ['yes', 'y'] else False
            # if four-point is desired, numerical is also set to True
            if self.do_four_point:
                self.numerical = True

        # Numerical derivative of dipole moment
        if 'dipole_deriv' in grad_dict:
            key = grad_dict['dipole_deriv'].lower()
            self.dipole_deriv = True if key in ['yes', 'y'] else False
            if self.dipole_deriv and not self.numerical:
                self.ostream.print_blank()
                warn_msg = '*** Warning: Dipole moment derivatives requested.'
                self.ostream.print_header(warn_msg.ljust(56))
                warn_msg = '             Gradient will be calculated numerically.'
                self.ostream.print_header(warn_msg.ljust(56))
                self.numerical = True


    def compute(self, molecule, ao_basis=None, min_basis=None):
        """
        Performs calculation of numerical gradient.

        :param molecule:
            The molecule.
        :param ao_basis:
            The AO basis set.
        :param min_basis:
            The minimal AO basis set.
        """

        return

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
            for j in range(i+1, natm):
                z_b = nuclear_charges[j]
                r_b = coords[j]
                r = np.sqrt(np.dot(r_a - r_b, r_a - r_b))
                f_ij = z_a * z_b * (r_b - r_a) / r**3
                nuc_contrib[i] += f_ij
                nuc_contrib[j] -= f_ij #z_a * z_b * (r_a - r_b) / r**3

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

        # atom labels
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


    def print_header(self):
        """
        Prints gradient calculation setup details to output stream.
        """

        self.ostream.print_blank()
        self.ostream.print_header(self.flag)
        self.ostream.print_header((len(self.flag) + 2) * '=')
        self.ostream.print_blank()
        self.ostream.flush()
