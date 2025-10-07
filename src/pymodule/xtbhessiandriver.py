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

import numpy as np
import time as tm

from .veloxchemlib import mpi_master
from .molecule import Molecule
from .hessiandriver import HessianDriver
from .xtbdriver import XtbDriver


class XtbHessianDriver(HessianDriver):
    """
    Implements XTB Hessian driver.

    :param xtb_drv:
        The XTB driver.

    Instance variables
        - flag: The driver flag.
    """

    def __init__(self, xtb_drv):
        """
        Initializes XTB Hessian driver.
        """

        super().__init__(xtb_drv.comm, xtb_drv.ostream)
        self.roots_to_follow = 0
        self.xtb_driver = xtb_drv
        self.flag = 'XTB Hessian Driver'
        self.numerical = True

    def update_settings(self, method_dict, hess_dict=None, cphf_dict=None):
        """
        Updates settings in XtbHessianDriver.

        :param method_dict:
            The input dictionary of method settings group.
        :param hess_dict:
            The input dictionary of Hessian settings group.
        :param cphf_dict:
            Dummy input to avoid crash in vibrationalanalysis driver
        """

        super().update_settings(method_dict, hess_dict)

        if hess_dict is None:
            hess_dict = {}
        if cphf_dict is None:
            cphf_dict = {}

    def compute(self, molecule, *args):
        """
        Computes the numerical nuclear Hessian.

        :param molecule:
            The molecule.
        :param args:
            Redundant args from call in vibrationalanalysis
        """

        self.print_header()

        start_time = tm.time()

        self.ostream.mute()
        self.xtb_driver.compute(molecule)
        self.ostream.unmute()

        self.elec_energy = self.xtb_driver.get_energy()

        self.compute_numerical(molecule)

        # print Hessian
        if self.do_print_hessian:
            self.print_geometry(molecule)
            self.ostream.print_blank()
            self.print_hessian(molecule)

        valstr = '*** Time spent in Hessian calculation: '
        valstr += '{:.2f} sec ***'.format(tm.time() - start_time)
        self.ostream.print_header(valstr)
        self.ostream.print_blank()
        self.ostream.flush()

    def compute_numerical(self, molecule):
        """
        Performs calculation of numerical Hessian.

        :param molecule:
            The molecule.
        """

        self.ostream.mute()

        # atom labels
        labels = molecule.get_labels()

        # number of atoms
        natm = molecule.number_of_atoms()

        # atom coordinates (nx3)
        coords = molecule.get_coordinates_in_bohr()

        # charge and spin multiplicity
        charge = molecule.get_charge()
        multiplicity = molecule.get_multiplicity()

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
                    new_mol.set_charge(charge)
                    new_mol.set_multiplicity(multiplicity)
                    # create a new XTB driver object;
                    # without this the energy is always zero...;
                    new_xtb_drv = XtbDriver(self.comm, self.ostream)
                    new_xtb_drv.set_method(self.xtb_driver.get_method())
                    new_xtb_drv.compute(new_mol)

                    grad_plus = new_xtb_drv.get_gradient()

                    mu_plus = new_xtb_drv.get_dipole()

                    coords[i, d] -= 2.0 * self.delta_h
                    new_mol = Molecule(labels, coords, units='au')
                    new_mol.set_charge(charge)
                    new_mol.set_multiplicity(multiplicity)
                    new_xtb_drv = XtbDriver(self.comm, self.ostream)
                    new_xtb_drv.set_method(self.xtb_driver.get_method())
                    new_xtb_drv.compute(new_mol)

                    grad_minus = new_xtb_drv.get_gradient()

                    mu_minus = new_xtb_drv.get_dipole()

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
                    new_mol.set_charge(charge)
                    new_mol.set_multiplicity(multiplicity)
                    new_xtb_drv = XtbDriver(self.comm, self.ostream)
                    new_xtb_drv.set_method(self.xtb_driver.get_method())
                    new_xtb_drv.compute(new_mol)

                    grad_plus1 = new_xtb_drv.get_gradient()

                    mu_plus1 = new_xtb_drv.get_dipole()

                    coords[i, d] += self.delta_h
                    new_mol = Molecule(labels, coords, units='au')
                    new_mol.set_charge(charge)
                    new_mol.set_multiplicity(multiplicity)
                    new_xtb_drv = XtbDriver(self.comm, self.ostream)
                    new_xtb_drv.set_method(self.xtb_driver.get_method())
                    new_xtb_drv.compute(new_mol)

                    grad_plus2 = new_xtb_drv.get_gradient()

                    mu_plus2 = new_xtb_drv.get_dipole()

                    coords[i, d] -= 3.0 * self.delta_h
                    new_mol = Molecule(labels, coords, units='au')
                    new_mol.set_charge(charge)
                    new_mol.set_multiplicity(multiplicity)
                    new_xtb_drv = XtbDriver(self.comm, self.ostream)
                    new_xtb_drv.set_method(self.xtb_driver.get_method())
                    new_xtb_drv.compute(new_mol)

                    grad_minus1 = new_xtb_drv.get_gradient()

                    mu_minus1 = new_xtb_drv.get_dipole()

                    coords[i, d] -= self.delta_h
                    new_mol = Molecule(labels, coords, units='au')
                    new_mol.set_charge(charge)
                    new_mol.set_multiplicity(multiplicity)
                    new_xtb_drv = XtbDriver(self.comm, self.ostream)
                    new_xtb_drv.set_method(self.xtb_driver.get_method())
                    new_xtb_drv.compute(new_mol)

                    grad_minus2 = new_xtb_drv.get_gradient()

                    mu_minus2 = new_xtb_drv.get_dipole()

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

        self.ostream.unmute()
