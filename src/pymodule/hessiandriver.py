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

from .veloxchemlib import (mpi_master, bohr_in_angstroms, avogadro_constant,
                           fine_structure_constant, electron_mass_in_amu,
                           amu_in_kg, speed_of_light_in_vacuum_in_SI)
from .outputstream import OutputStream
from .errorhandler import assert_msg_critical

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
        - mass_weighted_hessian: The mass-weighted Hessian
          in Hartree / (amu * Bohr**2).
        - reduced_masses: The reduced masses of the normal modes in amu.
        - force_constants: The force constants in mdyn/Angstrom.
        - frequencies: The vibrational frequencies in cm**-1.
        - normal_modes: The normalized vibrational normal modes in
                        (non-mass-weighted) Cartesian coordinates.
        - cart_normal_modes: The non-normalized vibrational modes in
                        (non-mass-weighted) Cartesian coordinates in 1/sqrt(amu).
        - dipole_gradient: The gradient of the dipole moment.
        - ir_intensities: The IR intensities in km/mol.
        - polarizability_gradient: The gradient of the polarizability.
        - raman_intensities: The Raman intensities (in ??).
        - flag: The type of Hessian driver.
        - numerical: Perform numerical Hessian calculation.
        - delta_h: Nuclear displacement for finite differences.
        - do_four_point: Perform four-point numerical approximation.
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
        self.mass_weighted_hessian = None
        self.reduced_masses = None
        self.force_constants = None
        self.frequencies = None
        self.normal_modes = None  # normalized, not mass-weighted
        self.cart_normal_modes = None  # neither normalized nor mass-weighted
        self.dipole_gradient = None
        self.ir_intensities = None
        self.polarizability_gradient = None
        self.raman_intensities = None
        self.flag = None

        self.numerical = True
        self.delta_h = 0.001

        # flag for two-point or four-point approximation
        self.do_four_point = False

        # flag for printing the Hessian
        self.do_print_hessian = False
        self.print_depolarization_ratio = False

        # Thermodynamics
        self.elec_energy = 0.0
        self.temperature = 298.15
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
            self.numerical = (key in ['yes', 'y'])

        if 'do_four_point' in freq_dict:
            key = freq_dict['do_four_point'].lower()
            self.do_four_point = (key in ['yes', 'y'])

        if self.do_four_point:
            self.numerical = True

        # print the Hessian (not mass-weighted)
        if 'do_print_hessian' in freq_dict:
            key = freq_dict['do_print_hessian'].lower()
            self.do_print_hessian = (key in ['yes', 'y'])

        # print the depolarization ratio, parallel, and perpendicular
        # Raman activities
        if 'print_depolarization_ratio' in freq_dict:
            key = freq_dict['print_depolarization_ratio'].lower()
            self.print_depolarization_ratio = (key in ['yes', 'y'])

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

    def vibrational_analysis(self, molecule, filename=None, rsp_drv=None):
        """
        Performs vibrational analysis (frequencies and normal modes)
        based on the molecular Hessian employing the geomeTRIC module:
        J. Chem, Phys. 144, 214108. DOI: 10.1063/1.4952956

        :param molecule:
            The molecule.
        :param filename:
            Filename where thermodynamic properties are saved by geomeTRIC.
        :param rsp_drv:
            The response driver (for excited state vibrational analysis).
        """

        geometric_repo = 'https://github.com/leeping/geomeTRIC.git'
        err_msg = ('The installed geometric package does not support\n' +
                   '  vibrational analysis. Please install the latest\n' +
                   '  geometric via\n' +
                   f'  python3 -m pip install git+{geometric_repo}\n')
        assert_msg_critical(hasattr(geometric, 'normal_modes'), err_msg)

        # number of atoms, elements, and coordinates
        natm = molecule.number_of_atoms()
        elem = molecule.get_labels()
        coords = molecule.get_coordinates().reshape(natm * 3)

        self.frequencies, self.normal_modes, gibbs_energy = (
            geometric.normal_modes.frequency_analysis(
                coords,
                self.hessian,
                elem,
                energy=self.elec_energy,
                temperature=self.temperature,
                pressure=self.pressure,
                outfnm=filename,
                normalized=False))

        # Diagonalizes Hessian and calculates the reduced masses
        self.reduced_masses = 1.0 / (np.einsum('ki->i', self.normal_modes.T**2))

        # Constants and conversion factors
        c = speed_of_light_in_vacuum_in_SI()
        alpha = fine_structure_constant()
        bohr_in_km = bohr_in_angstroms() * 1e-13
        cm_to_m = 1e-2  # centimeters in meters
        N_to_mdyne = 1e+8  # Newton in milli dyne
        m_to_A = 1e+10  # meters in Angstroms
        raman_conversion_factor = 0.078424

        # Conversion factor of IR intensity to km/mol
        conv_ir_ea0amu2kmmol = (electron_mass_in_amu() * avogadro_constant() *
                                alpha**2 * bohr_in_km * np.pi / 3.0)

        # Calculate force constants
        self.force_constants = (4.0 * np.pi**2 *
                                (c * (self.frequencies / cm_to_m))**2 *
                                self.reduced_masses *
                                amu_in_kg()) * (N_to_mdyne / m_to_A)

        natoms = molecule.number_of_atoms()
        atom_symbol = molecule.get_labels()
        number_of_modes = len(self.frequencies)

        # Calculate IR intensities (for ground state only)
        if self.dipole_gradient is not None:
            ir_trans_dipole = self.dipole_gradient.dot(self.normal_modes.T)
            ir_intensity_au_amu = np.array([
                np.linalg.norm(ir_trans_dipole[:, x])**2
                for x in range(ir_trans_dipole.shape[1])
            ])

            self.ir_intensities = ir_intensity_au_amu * conv_ir_ea0amu2kmmol

        # Calculate Raman intensities, if applicable
        if self.polarizability_gradient is not None:
            raman_transmom = np.einsum('xyi,ik->xyk',
                                       self.polarizability_gradient,
                                       self.normal_modes.T,
                                       optimize=True)
            # Calculate rotational invariants
            alpha_bar = np.zeros((number_of_modes))
            gamma_bar_sq = np.zeros((number_of_modes))
            for i in range(3):
                alpha_bar += raman_transmom[i, i] / 3
                for j in range(i + 1, 3):
                    gamma_bar_sq += 0.5 * (
                        raman_transmom[i, i] -
                        raman_transmom[j, j])**2 + 3 * raman_transmom[i, j]**2

            alpha_bar_sq = alpha_bar**2

            if self.print_depolarization_ratio:
                int_pol = 45 * alpha_bar_sq + 4 * gamma_bar_sq
                int_depol = 3 * gamma_bar_sq
                depol_ratio = int_depol / int_pol

            self.raman_intensities = (
                45 * alpha_bar_sq + 7 * gamma_bar_sq) * raman_conversion_factor

        # Now we can normalize the normal modes -- as done in geomeTRIC
        self.normal_modes /= np.linalg.norm(self.normal_modes,
                                            axis=1)[:, np.newaxis]

        title = 'Vibrational Analysis'
        self.ostream.print_header(title)
        self.ostream.print_header('=' * (len(title) + 2))
        self.ostream.print_blank()

        valstr = 'Harmonic frequencies (in cm**-1), force constants '
        valstr += '(in mdyne/A), reduced masses (in amu),'
        self.ostream.print_header(valstr)
        if self.ir_intensities is not None:
            valstr = ' IR intensities (in km/mol),'
            if self.raman_intensities is not None:
                valstr += ' Raman scattering activities (in A**4/amu),'
            self.ostream.print_header(valstr)
            if (self.raman_intensities is not None and
                    self.print_depolarization_ratio):
                valstr = ' parallel and perpendicular Raman '
                valstr += 'scattering activities,'
                valstr += ' depolarization ratios,'
                self.ostream.print_header(valstr)
        valstr = 'and Cartesian normal mode displacements.'
        self.ostream.print_header(valstr)
        self.ostream.print_blank()
        self.ostream.print_blank()

        for k in range(0, number_of_modes, 3):

            end = k + 3
            if k + 3 > number_of_modes:
                end = number_of_modes

            # Print indices and frequencies:
            index_string = '{:17s}'.format('  Index: ')
            freq_string = '{:17s}'.format('  Frequency: ')
            mass_string = '{:17s}'.format('  Reduced mass: ')
            force_cnst_string = '{:17s}'.format('  Force constant:')
            if self.ir_intensities is not None:
                ir_intens_string = '{:17s}'.format('  IR intensity:')
                if self.raman_intensities is not None:
                    raman_intens_string = '{:17s}'.format('  Raman activ.:')
                    if self.print_depolarization_ratio:
                        raman_parallel_str = '{:17s}'.format(
                            '  Parallel Raman:')
                        raman_perpendicular_str = '{:17s}'.format(
                            '  Perp. Raman:')
                        depolarization_str = '{:17s}'.format('  Depol. ratio:')

            normal_mode_string = '{:17s}'.format('  Normal mode: ')
            for i in range(k, end):
                index_string += '{:^31d}'.format(i + 1)
                freq_string += '{:^31.2f}'.format(self.frequencies[i])
                mass_string += '{:^31.4f}'.format(self.reduced_masses[i])
                force_cnst_string += '{:^31.4f}'.format(self.force_constants[i])
                if self.ir_intensities is not None:
                    ir_intens_string += '{:^31.4f}'.format(
                        self.ir_intensities[i])
                    if self.raman_intensities is not None:
                        raman_intens_string += '{:^31.4f}'.format(
                            self.raman_intensities[i])
                        if self.print_depolarization_ratio:
                            raman_parallel_str += '{:^31.4f}'.format(int_pol[i])
                            raman_perpendicular_str += '{:^31.4f}'.format(
                                int_depol[i])
                            depolarization_str += '{:^31.4f}'.format(
                                depol_ratio[i])

                normal_mode_string += '{:^30s}{:>1s}'.format(
                    'X         Y         Z', '|')
            self.ostream.print_line(index_string)
            self.ostream.print_line(freq_string)
            self.ostream.print_line(force_cnst_string)
            self.ostream.print_line(mass_string)
            if self.ir_intensities is not None:
                self.ostream.print_line(ir_intens_string)
                if self.raman_intensities is not None:
                    self.ostream.print_line(raman_intens_string)
                    if self.print_depolarization_ratio:
                        self.ostream.print_line(raman_parallel_str)
                        self.ostream.print_line(raman_perpendicular_str)
                        self.ostream.print_line(depolarization_str)
            self.ostream.print_line(normal_mode_string)

            # Print normal modes:
            for atom_index in range(natoms):
                valstr = '{:17s}'.format('  ' + str(atom_index + 1) + ' ' +
                                         atom_symbol[atom_index])

                for j in range(k, end):
                    valstr += '{:^10.4f}'.format(
                        self.normal_modes[j][3 * atom_index])  # X
                    valstr += '{:^10.4f}'.format(
                        self.normal_modes[j][3 * atom_index + 1])  # Y
                    valstr += '{:^10.4f}{:>1s}'.format(
                        self.normal_modes[j][3 * atom_index + 2], '|')  # Z
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
        nuc_contrib = np.zeros((natm, natm, 3, 3))

        # atom coordinates (nx3)
        coords = molecule.get_coordinates()

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

        self.ostream.print_blank()
        self.ostream.print_header(self.flag)
        self.ostream.print_header((len(self.flag) + 2) * '=')
        self.ostream.print_blank()
        self.ostream.flush()
