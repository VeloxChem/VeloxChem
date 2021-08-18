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
from .hessiandriver import HessianDriver
from .xtbgradientdriver import XTBGradientDriver
from .xtbdriver import XTBDriver
from .outputstream import OutputStream
from .firstorderprop import FirstOrderProperties #TODO: remove?
from .lrsolver import LinearResponseSolver #TODO: remove?

# For PySCF integral derivatives
#from .import_from_pyscf import overlap_deriv 
#from .import_from_pyscf import fock_deriv
#from .import_from_pyscf import eri_deriv


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
        - xtb_drv: The SCF driver.
        - do_raman: Additionally calculate Raman intensities
                (at significantly higher computational cost).
    """

    def __init__(self, xtb_drv, comm=None, ostream=None):
        """
        Initializes XTB Hessian driver.
        """

        super().__init__(comm, ostream)

        self.flag = 'XTB Hessian Driver'
        self.xtb_drv = xtb_drv
        self.do_raman = False


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

        # check if Raman intensities are to be calculated (numerically)
        if 'do_raman' in freq_dict:
            key = freq_dict['do_raman'].lower()
            self.do_raman = True if key in ['yes', 'y'] else False

        # The electronic energy 
        self.elec_energy = self.xtb_drv.get_energy()


    def compute(self, molecule, ao_basis, min_basis=None):
        """
        Computes the numerical nuclear Hessian.

        :param molecule:
            The molecule.
        :param ao_basis:
            The AO basis set.
        :param min_basis:
            The minimal AO basis set.
        """
        self.print_header()

        start_time = tm.time()

        self.compute_numerical(molecule, ao_basis, min_basis)

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

        # settings dictionary for gradient driver
        grad_dict = dict(self.freq_dict)

        #TODO -- check if working
        xtb_ostream_state = self.ostream.state
        self.ostream.state = False

        # atom labels
        labels = molecule.get_labels()

        # number of atoms
        natm = molecule.number_of_atoms()

        # atom coordinates (nx3)
        coords = molecule.get_coordinates()

        # Hessian
        hessian = np.zeros((natm, 3, natm, 3))

        # numerical dipole gradient (3 dipole components, no. atoms x 3 atom coords)
        self.dipole_gradient = np.zeros((3, 3 * molecule.number_of_atoms()))

       # # If Raman intensities are calculated, set up LR solver and member variable
       # if self.do_raman:
       #     # linear response driver for polarizability calculation
       #     lr_drv = LinearResponseSolver(self.comm, self.scf_drv.ostream)
       #     #lr_ostream_state = lr_drv.ostream.state
       #     #lr_drv.ostream.state = False
       #     # polarizability: 3 coordinates x 3 coordinates (ignoring frequencies)
       #     # polarizability gradient: dictionary goes through 3 coordinates x 3 coordinates
       #     # each entry having values for no. atoms x 3 coordinates
       #     self.pol_gradient = np.zeros((3, 3, 3 * natm))
       #     # dictionary to translate from numbers to operator components 'xyz'
       #     component_dict = {0: 'x', 1: 'y', 2: 'z'}


        if not self.do_four_point:
            for i in range(natm):
                for d in range(3):
                    coords[i, d] += self.delta_h
                    new_mol = Molecule(labels, coords, units='au')
                    # create a new XTB driver object;
                    # without this the energy is always zero...;
                    xtb_drv = XTBDriver(self.comm)
                    xtb_drv.compute(new_mol, self.ostream)
                    # create a new gradient driver
                    grad_drv = XTBGradientDriver(xtb_drv, self.comm,
                                                 self.ostream)
                    grad_drv.compute(new_mol)
                    grad_plus = grad_drv.get_gradient()

                    mu_plus = xtb_drv.get_dipole()

                    #if self.do_raman:
                    # TODO: add get_polarizability in the C layer
                    # if this is possible to get from XTB

                    coords[i, d] -= 2.0 * self.delta_h
                    new_mol = Molecule(labels, coords, units='au')
                    xtb_drv = XTBDriver(self.comm)
                    xtb_drv.compute(new_mol, self.ostream)
                    # create a new gradient driver
                    grad_drv = XTBGradientDriver(xtb_drv, self.comm,
                                                 self.ostream)
                    grad_drv.compute(new_mol)
                    grad_minus = grad_drv.get_gradient()

                    mu_minus = xtb_drv.get_dipole()

                    #if self.do_raman:
                    # TODO: add get_polarizability in the C layer
                    # if this is possible to get from XTB
                    for c in range(3):
                        self.dipole_gradient[c, 3*i + d] = (mu_plus[c] - mu_minus[c]) / (2.0 * self.delta_h)
                    coords[i, d] += self.delta_h
                    hessian[i, d, :, :] = (grad_plus - grad_minus) / (2.0 * self.delta_h)


        else:
            # Four-point numerical derivative approximation
            # for debugging of analytical Hessian:
            # [ f(x - 2h) - 8 f(x - h) + 8 f(x + h) - f(x + 2h) ] / ( 12h )
            for i in range(natm):
                for d in range(3):
                    coords[i, d] += self.delta_h
                    new_mol = Molecule(labels, coords, units='au')
                    xtb_drv = XTBDriver(self.comm)
                    xtb_drv.compute(new_mol, self.ostream)
                    # create a new gradient driver
                    grad_drv = XTBGradientDriver(xtb_drv, self.comm,
                                                 self.ostream)
                    grad_drv.compute(new_mol)
                    grad_plus1 = grad_drv.get_gradient()

                    mu_plus1 = xtb_drv.get_dipole()

                    #if self.do_raman:
                    # TODO: add get_polarizability in the C layer
                    # if this is possible to get from XTB

                    coords[i, d] += self.delta_h
                    new_mol = Molecule(labels, coords, units='au')
                    xtb_drv = XTBDriver(self.comm)
                    xtb_drv.compute(new_mol, self.ostream)
                    # create a new gradient driver
                    grad_drv = XTBGradientDriver(xtb_drv, self.comm,
                                                 self.ostream)
                    grad_drv.compute(new_mol)
                    grad_plus2 = grad_drv.get_gradient()

                    mu_plus2 = xtb_drv.get_dipole()

                    #if self.do_raman:
                    # TODO: add get_polarizability in the C layer
                    # if this is possible to get from XTB

                    coords[i, d] -= 3.0 * self.delta_h
                    new_mol = Molecule(labels, coords, units='au')
                    xtb_drv = XTBDriver(self.comm)
                    xtb_drv.compute(new_mol, self.ostream)
                    # create a new gradient driver
                    grad_drv = XTBGradientDriver(xtb_drv, self.comm,
                                                 self.ostream)
                    grad_drv.compute(new_mol)
                    grad_minus1 = grad_drv.get_gradient()

                    mu_minus1 = xtb_drv.get_dipole()

                    #if self.do_raman:
                    # TODO: add get_polarizability in the C layer
                    # if this is possible to get from XTB

                    coords[i, d] -= self.delta_h
                    new_mol = Molecule(labels, coords, units='au')
                    xtb_drv = XTBDriver(self.comm)
                    xtb_drv.compute(new_mol, self.ostream)
                    # create a new gradient driver
                    grad_drv = XTBGradientDriver(xtb_drv, self.comm,
                                                 self.ostream)
                    grad_drv.compute(new_mol)
                    grad_minus2 = grad_drv.get_gradient()

                    mu_minus2 = xtb_drv.get_dipole()

                    #if self.do_raman:
                    # TODO: add get_polarizability in the C layer
                    # if this is possible to get from XTB
                    for c in range(3):
                        self.dipole_gradient[c, 3*i + d] = (mu_minus2[c] - 8.0 * mu_minus1[c]
                                                         + 8.0 * mu_plus1[c] - mu_plus2[c]) / (12.0 * self.delta_h)

                    coords[i, d] += 2.0 * self.delta_h
                    # f'(x) ~ [ f(x - 2h) - 8 f(x - h) + 8 f(x + h) - f(x + 2h) ] / ( 12h )
                    hessian[i, d] = (grad_minus2 - 8.0 * grad_minus1
                                           + 8.0 * grad_plus1 - grad_plus2) / (12.0 * self.delta_h)

        # reshaped Hessian as member variable
        self.hessian = hessian.reshape(3*natm, 3*natm)

        #self.ostream.print_blank()

        #self.xtb_drv.compute(molecule, self.ostream)
        self.ostream.state = xtb_ostream_state

