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
import time as tm

from .veloxchemlib import mpi_master
from .molecule import Molecule
from .hessiandriver import HessianDriver
from .firstorderprop import FirstOrderProperties
from .lrsolver import LinearResponseSolver


class ScfHessianDriver(HessianDriver):
    """
    Implements SCF Hessian driver.

    :param scf_drv:
        The SCF driver.
    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.

    Instance variables
        - scf_drv: The SCF driver.
        - do_raman: Additionally calculate Raman intensities
                (at significantly higher computational cost).
    """

    def __init__(self, scf_drv, comm=None, ostream=None):
        """
        Initializes SCF Hessian driver.
        """

        super().__init__(comm, ostream)

        self.flag = 'SCF Hessian Driver'
        self.scf_drv = scf_drv
        self.do_raman = False
        self.pople = False

        # Solver setup
        self.conv_thresh = 1.0e-4
        self.max_iter = 50
        self.iter_count = 0
        self.is_converged = False

    def update_settings(self, method_dict, freq_dict=None, cphf_dict=None):
        """
        Updates settings in ScfHessianDriver.

        :param method_dict:
            The input dictionary of method settings group.
        :param freq_dict:
            The input dictionary of Hessian/frequency settings group.
        :param cphf_dict:
            The input dictionary of CPHF (orbital response) settings.
        """

        super().update_settings(method_dict, freq_dict)

        if freq_dict is None:
            freq_dict = {}

        # Settings for orbital response module
        if cphf_dict is None:
            cphf_dict = {}
        self.cphf_dict = dict(cphf_dict)

        if 'conv_thresh' in cphf_dict:
            self.conv_thresh = float(cphf_dict['conv_thresh'])

        if 'max_iter' in cphf_dict:
            self.max_iter = int(cphf_dict['max_iter'])

        # check if Raman intensities are to be calculated (numerically)
        if 'do_raman' in freq_dict:
            key = freq_dict['do_raman'].lower()
            self.do_raman = (key in ['yes', 'y'])

        # The electronic energy
        self.elec_energy = self.scf_drv.get_scf_energy()

    def compute(self, molecule, ao_basis, min_basis=None):
        """
        Computes the analytical or numerical nuclear Hessian.

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

    def compute_numerical(self, molecule, ao_basis, min_basis):
        """
        Performs the calculation of a numerical Hessian based only
        on the energy.

        :param molecule:
            The molecule.
        :param ao_basis:
            The AO basis set.
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

        # Hessian
        hessian = np.zeros((natm, 3, natm, 3))

        # First-order properties for gradient of dipole moment
        prop = FirstOrderProperties(self.comm, self.ostream)
        # numerical gradient (3 dipole components, no. atoms x 3 atom coords)
        if self.rank == mpi_master():
            self.dipole_gradient = np.zeros((3, 3 * natm))

        # If Raman intensities are calculated,
        # set up LR solver and member variable
        if self.do_raman:
            # linear response driver for polarizability calculation
            lr_drv = LinearResponseSolver(self.comm, self.scf_drv.ostream)
            # polarizability: 3 coordinates x 3 coordinates
            # (ignoring frequencies)
            # polarizability gradient: dictionary goes through
            # 3 coordinates x 3 coordinates
            # each entry having values for no. atoms x 3 coordinates
            if self.rank == mpi_master():
                self.polarizability_gradient = np.zeros((3, 3, 3 * natm))

        self.scf_drv.compute(molecule, ao_basis, min_basis)
        energy_0 = self.scf_drv.get_scf_energy()

        for i in range(natm):
            for x in range(3):
                # Plus x
                coords[i, x] += self.delta_h
                new_mol = Molecule(labels, coords, units='au')
                self.scf_drv.compute(new_mol, ao_basis, min_basis)
                energy_ixp = self.scf_drv.get_scf_energy()

                prop.compute_scf_prop(new_mol, ao_basis,
                                      self.scf_drv.scf_tensors)
                if self.rank == mpi_master():
                    mu_plus = prop.get_property('dipole moment')

                if self.do_raman:
                    lr_drv.is_converged = False
                    lr_results_p = lr_drv.compute(new_mol, ao_basis,
                                                  self.scf_drv.scf_tensors)

                # Minus x
                coords[i, x] -= 2.0 * self.delta_h
                new_mol = Molecule(labels, coords, units='au')
                self.scf_drv.compute(new_mol, ao_basis, min_basis)
                energy_ixm = self.scf_drv.get_scf_energy()

                prop.compute_scf_prop(new_mol, ao_basis,
                                      self.scf_drv.scf_tensors)
                if self.rank == mpi_master():
                    mu_minus = prop.get_property('dipole moment')

                if self.do_raman:
                    lr_drv.is_converged = False
                    lr_results_m = lr_drv.compute(new_mol, ao_basis,
                                                  self.scf_drv.scf_tensors)
                    if self.rank == mpi_master():
                        for ind_aop, aop in enumerate('xyz'):
                            for ind_bop, bop in enumerate('xyz'):
                                # TODO: here the freq. is hard-coded to 0.0!
                                comp_plus = (
                                    lr_results_p['response_functions'][aop, bop,
                                                                       0.0])
                                comp_minus = (
                                    lr_results_m['response_functions'][aop, bop,
                                                                       0.0])
                                self.polarizability_gradient[
                                    ind_aop, ind_bop,
                                    3 * i + x] = ((comp_plus - comp_minus) /
                                                  (2.0 * self.delta_h))

                if self.rank == mpi_master():
                    for c in range(3):
                        self.dipole_gradient[c, 3 * i +
                                             x] = ((mu_plus[c] - mu_minus[c]) /
                                                   (2.0 * self.delta_h))

                hessian[i, x, i,
                        x] = ((energy_ixp - 2 * energy_0 + energy_ixm) /
                              self.delta_h**2)
                coords[i, x] += self.delta_h

                for j in range(i, natm):
                    for y in range(3):
                        if (j == i and x != y) or (j != i):
                            # Plus y
                            coords[j, y] += self.delta_h
                            new_mol = Molecule(labels, coords, units='au')
                            self.scf_drv.compute(new_mol, ao_basis, min_basis)
                            energy_jyp = self.scf_drv.get_scf_energy()

                            # Plus x, plus y
                            coords[i, x] += self.delta_h
                            new_mol = Molecule(labels, coords, units='au')
                            self.scf_drv.compute(new_mol, ao_basis, min_basis)
                            energy_ixp_jyp = self.scf_drv.get_scf_energy()
                            coords[i, x] -= self.delta_h

                            # Minus y
                            coords[j, y] -= 2.0 * self.delta_h
                            new_mol = Molecule(labels, coords, units='au')
                            self.scf_drv.compute(new_mol, ao_basis, min_basis)
                            energy_jym = self.scf_drv.get_scf_energy()

                            # Minus x, minus y:
                            coords[i, x] -= self.delta_h
                            new_mol = Molecule(labels, coords, units='au')
                            self.scf_drv.compute(new_mol, ao_basis, min_basis)
                            energy_ixm_jym = self.scf_drv.get_scf_energy()

                            coords[i, x] += self.delta_h
                            coords[j, y] += self.delta_h

                            hessian[i, x, j, y] = (
                                (energy_ixp_jyp - energy_ixp - energy_jyp +
                                 2 * energy_0 - energy_ixm - energy_jym +
                                 energy_ixm_jym) / (2 * self.delta_h**2))
                            hessian[j, y, i, x] = hessian[i, x, j, y]
        # reshaped Hessian as member variable
        self.hessian = hessian.reshape(3 * natm, 3 * natm)

        # restore scf_drv to initial state
        self.scf_drv.compute(molecule, ao_basis, min_basis)
        self.scf_drv.ostream.state = scf_ostream_state
