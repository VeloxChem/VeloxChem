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
from geometric import normal_modes

class HessianDriver:
    """
    Implements the Hessian driver.

    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.

    Instance variables
        - hessian: The Hessian.
        - flag: The type of Hessian driver.
        - numerical: Perform numerical Hessian calculation.
        - numerical_grad: Perform numerical gradient calculation.
        - delta_h: Nuclear displacement for finite differences.
        - do_four_point: Perform four-point numerical approximation.
        - print_vib_analysis: Print vibrational analysis (frequencies and normal modes)
        - do_print_hessian: Flag for printing the Hessian.
        - elec_energy: The (total) electronic energy.
        - temperature: The temperature (in K) used for thermodynamic analysis.
        - pressure: The pressure (in bar) used for thermodynamic analysis.
    """

    def __init__(self, comm=None, ostream=None):
        """
        Initializes Hessian driver.
        """

        if comm is None:
            comm = MPI.COMM_WORLD

        if ostream is None:
            ostream = OutputStream(sys.stdout)

        self.comm = comm
        self.ostream = ostream

        self.hessian = None
        self.flag = None
        self.numerical = True
        self.numerical_grad = False
        self.delta_h = 0.001
        # flag for two-point or four-point approximation
        self.do_four_point = False
        self.print_vib_analysis = True
        # flag for printing the Hessian
        self.do_print_hessian = False

        # Thermodynamics
        self.elec_energy = 0.0
        self.temperature = 300
        self.pressure = 1.0


    def update_settings(self, method_dict, freq_dict=None):
        """
        Updates settings in HessianDriver.

        :param method_dict:
            The input dictionary of method settings group.
        :param freq_dict:
            The input dictionary of Hessian/frequency settings group.
        """

        if freq_dict is None:
            freq_dict = {}

        # check if Hessianis to be calculated numerically
        if 'numerical' in freq_dict:
            key = freq_dict['numerical'].lower()
            self.numerical = True if key in ['yes', 'y'] else False
            # TODO: analytical Hessian not yet implemented
            if not self.numerical:
                self.numerical = True
                warn_msg = '*** Warning: Analytical Hessian is not yet implemented.'
                self.ostream.print_header(warn_msg.ljust(56))
                warn_msg = '    Hessian will be calculated numerically instead.'
                self.ostream.print_header(warn_msg.ljust(56))
                self.ostream.print_blank()
                self.ostream.flush()

        # check if gradient is to be calculated numerically
        if 'numerical_grad' in freq_dict:
            key = freq_dict['numerical_grad'].lower()
            self.numerical_grad = True if key in ['yes', 'y'] else False

        # if gradient is calculated numerically, so is the Hessian
        if self.numerical_grad:
            self.numerical = True

        # Analytical DFT gradient/Hessian is not implemented yet
        if 'xcfun' in method_dict:
            if method_dict['xcfun'] is not None:
                self.numerical = True
                self.numerical_grad = True
        if 'dft' in method_dict:
            key = method_dict['dft'].lower()
            if key in ['yes', 'y']:
                self.numerical = True
                self.numerical_grad = True

        # print vibrational analysis (frequencies and normal modes)
        if 'print_vib_analysis' in freq_dict:
            key = freq_dict['print_vib_analysis'].lower()
            self.print_vib_analysis = True if key in ['yes', 'y'] else False

        # print the Hessian (not mass-weighted)
        if 'do_print_hessian' in freq_dict:
            key = freq_dict['do_print_hessian'].lower()
            self.do_print_hessian = True if key in ['yes', 'y'] else False

        if 'temperature' in freq_dict:
            self.temperature = float(freq_dict['temperature'])
        if 'pressure' in freq_dict:
            self.pressure = float(freq_dict['pressure'])

        self.method_dict = dict(method_dict)
        self.freq_dict = dict(freq_dict)


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

    def vibrational_analysis(self, molecule, ao_basis, filename=None):
        """
        Performs vibrational analysis (frequencies and normal modes)
        based on the molecular Hessian employing the geomeTRIC module.

        :param molecule:
            The molecule.
        :param ao_basis:
            The AO basis set.
        :param filename:
            Filename where thermodynamic properties are saved by geomeTRIC.
        """

        # number of atoms, elements, and coordinates
        natm = molecule.number_of_atoms()
        elem = molecule.get_labels()
        coords = molecule.get_coordinates().reshape(natm*3)

        normal_modes.frequency_analysis(coords, hessian_drv.hessian, elem,
                        energy=self.elec_energy, temperatur=self.temperature,
                        pressure=self.pressure, outfnm=filename)


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
        nuc_contrib = np.zeros((natm, 3, natm, 3))

        # atom coordinates (nx3)
        coords = molecule.get_coordinates()

        # atomic charges
        nuclear_charges = molecule.elem_ids_to_numpy()

        # loop over all distinct atom pairs and add energy contribution
        for i in range(natm):
            z_a = nuclear_charges[i]
            r_a = coords[i]
            for j in range(natm):
                if i != j:
                    z_b = nuclear_charges[j]
                    r_b = coords[j]
                    r = np.sqrt(np.dot(r_a - r_b, r_a - r_b))
                    for k in range(3):
                        for l in range(3):
                    # off-diagonal parts
                            nuc_contrib[i, k, j, l] = - 3*z_a*z_b*(r_b[k] - r_a[k])*(r_b[l] - r_a[l]) / r**5
                            if k == l:
                                nuc_contrib[i, k, j, l] += z_a * z_b / r**3

                    # add the diagonal contribution
                            nuc_contrib[i, k, i, l] += 3*z_a*z_b*(r_b[k] - r_a[k])*(r_b[l] - r_a[l]) / r**5
                            if k == l:
                                nuc_contrib[i, k, i, l] -= z_a * z_b / r**3


        return nuc_contrib.reshape(3*natm, 3*natm)


    def diagonalize_hessian(self, molecule, basis, rsp_drv=None):
        """
        Diagonalizes the Hessian matrix to obtain force constants as eigenvalues
        (and hence vibrational frequencies) and normal modes as eigenvectors.

        :param molecule:
            The molecule.
        :param basis:
            The AO basis set.
        :param rsp_drv:
            The RPA or TDA driver (for TdhfHessianDriver).

        :return:
            The vibrational frequencies in atomic units.
        """
        # compute the Hessian if not done already
        if self.hessian is None:
            if rsp_drv is None:
                self.compute(molecule, basis)
            else:
                # TdhfHessianDriver inherits this function
                self.compute(molecule, basis, rsp_drv)

        natm = molecule.number_of_atoms()

        # take the square root of the atomic masses, repeat each three times (xyz components)
        # then form the direct product matrix
        masses_sqrt = np.sqrt(molecule.masses_to_numpy())
        masses_repeat = np.repeat(masses_sqrt, 3)
        masses_matrix = masses_repeat.reshape(-1, 1) * masses_repeat

        # reshape the Hessian as 3Nx3N and mass-weight it
        reshaped_hessian = self.hessian.reshape(3*natm, 3*natm)
        self.mass_weighted_hessian = reshaped_hessian / masses_matrix

        # diagonalize the mass-weighted Hessian
        hessian_eigvals, hessian_eigvecs = np.linalg.eigh(self.mass_weighted_hessian)

        # the first 6 elements should be close to zero (translation & rotation)
        return hessian_eigvals[6:]


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

            valstr = '{:15s}'.format('Coord. ')

            coord_dict = {0:'(x)', 1:'(y)', 2:'(z)'}
            end = k + 2
            if k + 2 > natm:
                 end = natm
            for i in range(k,end):
                for di in range(3):
                    valstr += '{:16s}'.format(''+str(i+1)+' '+labels[i]+coord_dict[di]+'')

            self.ostream.print_line(valstr)
            self.ostream.print_blank()


            for i in range(natm):
                for di in range(3):
                    valstr = '  {:7s}'.format(str(i+1)+' '+labels[i]+coord_dict[di])
                    for j in range(k,end):
                        #print("j = %d" % j)
                        for dj in range(3):
                            valstr += '{:16.8f}'.format(self.hessian[i*3+di,j*3+dj])
                    self.ostream.print_line(valstr)
            self.ostream.print_blank()
            self.ostream.print_blank()
        # self.ostream.print_blank()
        self.ostream.flush()


    def print_header(self):
        """
        Prints Hessian calculation setup details to output stream.
        """

        self.ostream.print_blank()
        self.ostream.print_header(self.flag)
        self.ostream.print_header((len(self.flag) + 2) * '=')
        self.ostream.print_blank()
        self.ostream.flush()



