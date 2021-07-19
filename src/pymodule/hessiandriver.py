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
from .veloxchemlib import bohr_in_angstroms
from .veloxchemlib import dipole_in_debye

import geometric


class HessianDriver:
    """
    Implements the Hessian driver.

    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.

    Instance variables
        - hessian: The Hessian.
        - mass_weighted_hessian: The mass-weighted Hessian.
        - reduced_masses: The reduced masses of the normal modes.
        - force_constants: The force constants.
        - frequencies: The vibrational frequencies.
        - normal_modes: The vibrational normal modes in
                        (non-mass-weighted) Cartesian coordinates.
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
        self.mass_weighted_hessian = None
        self.reduced_masses = None
        self.force_constants = None
        self.frequencies = None
        self.normal_modes = None
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
                self.ostream.print_blank()
                warn_msg = '*** Warning: Analytical Hessian is not yet implemented.'
                self.ostream.print_header(warn_msg.ljust(56))
                warn_msg = '    Hessian will be calculated numerically instead.'
                self.ostream.print_header(warn_msg.ljust(56))
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

    def vibrational_analysis(self, molecule, ao_basis, filename=None, rsp_drv=None):
        """
        Performs vibrational analysis (frequencies and normal modes)
        based on the molecular Hessian employing the geomeTRIC module.

        :param molecule:
            The molecule.
        :param ao_basis:
            The AO basis set.
        :param filename:
            Filename where thermodynamic properties are saved by geomeTRIC.
        :param rsp_drv:
            The response driver (for excited state vibrational analysis).
        """

        # number of atoms, elements, and coordinates
        natm = molecule.number_of_atoms()
        elem = molecule.get_labels()
        coords = molecule.get_coordinates().reshape(natm*3)

        # square root of atomic masses needed for computing the reduced mass
        masses_sqrt = np.sqrt(molecule.masses_to_numpy())
        masses_sqrt_repeat = np.repeat(masses_sqrt, 3)

        self.frequencies, self.normal_modes, gibbs_energy = (
                        geometric.normal_modes.frequency_analysis(coords,
                                self.hessian, elem,
                                energy=self.elec_energy,
                                temperature=self.temperature,
                                pressure=self.pressure, outfnm=filename)
                                                    )


        # Diagonalizes Hessian and calculates the reduced masses
        self.diagonalize_hessian(molecule, ao_basis, rsp_drv)

        # Constants and conversion factors
        # TODO: get these from the proper place.
        c = 2.99792458e8
        cm_to_m = 1e-2
        amu_to_kg = 1.6605390666e-27
        N_to_mdyne = 1e8
        m_to_A = 1e10

        self.force_constants = ( 4.0 * np.pi**2
                            * (c * (self.frequencies / cm_to_m) )**2
                            * self.reduced_masses * amu_to_kg
                            )  * ( N_to_mdyne / m_to_A )

        number_of_modes = len(self.frequencies)

        natoms = molecule.number_of_atoms()
        atom_symbol = molecule.get_labels()

        title = 'Vibrational Analysis'
        self.ostream.print_header(title)
        self.ostream.print_header('=' * (len(title) + 2))
        self.ostream.print_blank()

        valstr = 'Frequencies (in cm**-1), reduced masses (in amu),'
        valstr += ' force constants (in mdyne/A), and normal modes.'
        self.ostream.print_header(valstr)
        self.ostream.print_blank()

        for k in range(0, number_of_modes, 3):

            end = k + 3
            if k + 3 > number_of_modes:
                 end = number_of_modes

            # Print indices and frequencies:
            index_string = '{:17s}'.format('  Index: ')
            freq_string =  '{:17s}'.format('  Frequency: ')
            mass_string =  '{:17s}'.format('  Reduced mass: ')
            force_cnst_string =  '{:17s}'.format('  Force constant:')
            normal_mode_string = '{:17s}'.format('  Normal mode: ')
            for i in range(k, end):
                index_string += '{:^31d}'.format(i+1)
                freq_string +=  '{:^31.2f}'.format(self.frequencies[i])
                mass_string += '{:^31.4f}'.format(self.reduced_masses[i])
                force_cnst_string += '{:^31.4f}'.format(self.force_constants[i])
                normal_mode_string += '{:^30s}{:>1s}'.format('X         Y         Z','|')
            self.ostream.print_line(index_string)
            self.ostream.print_line(freq_string)
            self.ostream.print_line(mass_string)
            self.ostream.print_line(force_cnst_string)
            self.ostream.print_line(normal_mode_string)

            # Print normal modes:
            for atom_index in range(natoms):
                valstr =  '{:17s}'.format('  '+str(atom_index+1)+' '+atom_symbol[atom_index])

                for j in range(k, end):
                    valstr += '{:^10.4f}'.format(self.normal_modes[j][3*atom_index]) # X
                    valstr += '{:^10.4f}'.format(self.normal_modes[j][3*atom_index+1]) # Y
                    valstr += '{:^10.4f}{:>1s}'.format(self.normal_modes[j][3*atom_index+2],'|') # Z
                self.ostream.print_line(valstr)

            self.ostream.print_blank()

        self.ostream.print_blank()

        self.ostream.flush()


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
        masses_sqrt_repeat = np.repeat(masses_sqrt, 3)
        masses_matrix = masses_sqrt_repeat.reshape(-1, 1) * masses_sqrt_repeat

        # reshape the Hessian as 3Nx3N and mass-weigh it
        reshaped_hessian = self.hessian.reshape(3*natm, 3*natm)
        self.mass_weighted_hessian = reshaped_hessian / masses_matrix

        # diagonalize the mass-weighted Hessian
        hessian_eigvals, hessian_eigvecs = np.linalg.eigh(self.mass_weighted_hessian)

        cart_normal_modes = np.einsum('k,ki->ki',
                                      1/masses_sqrt_repeat,
                                      hessian_eigvecs)


        reduced_masses = 1.0/(np.einsum('ki->i',cart_normal_modes**2))
        self.cart_normal_modes = cart_normal_modes

        # TODO: linear molecules only for diatomics, needs to be generalized
        if natm == 2:
            self.reduced_masses = reduced_masses[5:]
            return hessian_eigvals[5:]
        else:
            self.reduced_masses = reduced_masses[6:]
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

            valstr = '{:15s}'.format('  Coord. ')

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



