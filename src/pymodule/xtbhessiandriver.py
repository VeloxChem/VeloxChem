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

import numpy as np
import time as tm

from .veloxchemlib import mpi_master
from .molecule import Molecule
from .hessiandriver import HessianDriver
from .xtbdriver import XTBDriver


class XTBHessianDriver(HessianDriver):
    """
    Implements XTB Hessian driver.

    :param xtb_drv:
        The XTB driver.
    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.

    Instance variables
        - xtb_drv: The XTB driver.
    """

    def __init__(self, xtb_drv, comm=None, ostream=None):
        """
        Initializes XTB Hessian driver.
        """

        super().__init__(comm, ostream)

        self.flag = 'XTB Hessian Driver'
        self.xtb_drv = xtb_drv

    def update_settings(self, method_dict, freq_dict=None):
        """
        Updates settings in XTBHessianDriver.

        :param method_dict:
            The input dictionary of method settings group.
        :param freq_dict:
            The input dictionary of Hessian/frequency settings group.
        """

        super().update_settings(method_dict, freq_dict)

        if freq_dict is None:
            freq_dict = {}

        # The electronic energy
        self.elec_energy = self.xtb_drv.get_energy()

    def compute(self, molecule):
        """
        Computes the numerical nuclear Hessian.

        :param molecule:
            The molecule.
        """
        self.print_header()

        start_time = tm.time()

        self.compute_numerical(molecule)

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

    def compute_numerical(self, molecule):
        """
        Performs calculation of numerical Hessian.

        :param molecule:
            The molecule.
        """

        xtb_ostream_state = self.ostream.state
        self.ostream.state = False

        # atom labels
        labels = molecule.get_labels()

        # number of atoms
        natm = molecule.number_of_atoms()

        # atom coordinates (nx3)
        coords = molecule.get_coordinates()

        if self.rank == mpi_master():
            # Hessian
            hessian = np.zeros((natm, 3, natm, 3))

            # numerical dipole gradient (3 dipole components,
            # no. atoms x 3 atom coords)
            self.dipole_gradient = np.zeros((3, 3 * natm))

        if not self.do_four_point:
            for i in range(natm):
                for d in range(3):
                    coords[i, d] += self.delta_h
                    new_mol = Molecule(labels, coords, units='au')
                    # create a new XTB driver object;
                    # without this the energy is always zero...;
                    xtb_drv = XTBDriver(self.comm)
                    xtb_drv.compute(new_mol, self.ostream)

                    grad_plus = xtb_drv.get_gradient()

                    mu_plus = xtb_drv.get_dipole()

                    coords[i, d] -= 2.0 * self.delta_h
                    new_mol = Molecule(labels, coords, units='au')
                    xtb_drv = XTBDriver(self.comm)
                    xtb_drv.compute(new_mol, self.ostream)

                    grad_minus = xtb_drv.get_gradient()

                    mu_minus = xtb_drv.get_dipole()

                    coords[i, d] += self.delta_h

                    if self.rank == mpi_master():
                        for c in range(3):
                            self.dipole_gradient[c, 3 * i + d] = (
                                mu_plus[c] - mu_minus[c]) / (2.0 * self.delta_h)

                        hessian[i,
                                d, :, :] = (grad_plus -
                                            grad_minus) / (2.0 * self.delta_h)

        else:
            # Four-point numerical derivative approximation
            # [ f(x - 2h) - 8 f(x - h) + 8 f(x + h) - f(x + 2h) ] / ( 12h )
            for i in range(natm):
                for d in range(3):
                    coords[i, d] += self.delta_h
                    new_mol = Molecule(labels, coords, units='au')
                    xtb_drv = XTBDriver(self.comm)
                    xtb_drv.compute(new_mol, self.ostream)

                    grad_plus1 = xtb_drv.get_gradient()

                    mu_plus1 = xtb_drv.get_dipole()

                    coords[i, d] += self.delta_h
                    new_mol = Molecule(labels, coords, units='au')
                    xtb_drv = XTBDriver(self.comm)
                    xtb_drv.compute(new_mol, self.ostream)

                    grad_plus2 = xtb_drv.get_gradient()

                    mu_plus2 = xtb_drv.get_dipole()

                    coords[i, d] -= 3.0 * self.delta_h
                    new_mol = Molecule(labels, coords, units='au')
                    xtb_drv = XTBDriver(self.comm)
                    xtb_drv.compute(new_mol, self.ostream)

                    grad_minus1 = xtb_drv.get_gradient()

                    mu_minus1 = xtb_drv.get_dipole()

                    coords[i, d] -= self.delta_h
                    new_mol = Molecule(labels, coords, units='au')
                    xtb_drv = XTBDriver(self.comm)
                    xtb_drv.compute(new_mol, self.ostream)

                    grad_minus2 = xtb_drv.get_gradient()

                    mu_minus2 = xtb_drv.get_dipole()

                    coords[i, d] += 2.0 * self.delta_h

                    if self.rank == mpi_master():
                        for c in range(3):
                            self.dipole_gradient[
                                c, 3 * i +
                                d] = (mu_minus2[c] - 8.0 * mu_minus1[c] +
                                      8.0 * mu_plus1[c] -
                                      mu_plus2[c]) / (12.0 * self.delta_h)

                        # f'(x) ~ [ f(x - 2h) - 8 f(x - h)
                        # + 8 f(x + h) - f(x + 2h) ] / ( 12h )
                        hessian[i,
                                d, :, :] = (grad_minus2 - 8.0 * grad_minus1 +
                                            8.0 * grad_plus1 -
                                            grad_plus2) / (12.0 * self.delta_h)

        # reshaped Hessian as member variable
        if self.rank == mpi_master():
            self.hessian = hessian.reshape(3 * natm, 3 * natm)
        else:
            self.hessian = None

        self.ostream.state = xtb_ostream_state
