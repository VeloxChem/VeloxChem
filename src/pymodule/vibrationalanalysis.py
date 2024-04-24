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

from .scfgradientdriver import ScfGradientDriver
from .scfhessiandriver import ScfHessianDriver
from .polarizabilitygradient import PolarizabilityGradient
from .lrsolver import LinearResponseSolver
from .cppsolver import ComplexResponse
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


class VibrationalAnalysis:
    """
    Implements the vibrational analysis driver.

    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.
    :param hess_drv:
        The Hessian driver.

    Instance variables
        - hessian: The Hessian in Hartree per Bohr**2.
        - mass_weighted_hessian: The mass-weighted Hessian
          in Hartree / (amu * Bohr**2).
        - reduced_masses: The reduced masses of the normal modes in amu.
        - force_constants: The force constants in mdyn/Angstrom.
        - vib_frequencies: The vibrational frequencies in cm**-1.
        - normal_modes: The normalized vibrational normal modes in
                        (non-mass-weighted) Cartesian coordinates.
        - cart_normal_modes: The non-normalized vibrational modes in
                        (non-mass-weighted) Cartesian coordinates in 1/sqrt(amu).
        - dipole_gradient: The gradient of the dipole moment.
        - ir_intensities: The IR intensities in km/mol.
        - polarizability_gradient: The gradient of the polarizability.
        - raman_intensities: The Raman intensities (in A**4/amu).
        - frequencies: the frequency/ies of external electric field (for resonance Raman)
        - flag: The type of Hessian driver.
        - numerical_hessian: Perform numerical Hessian calculation.
        - numerical_raman: Perform numerical polarizability gradient calculation.
        - do_four_point_hessian: Perform four-point numerical approximation.
        - do_four_point_raman: Perform four-point numerical approximation.
        - do_print_hessian: Flag for printing the Hessian.
        - elec_energy: The (total) electronic energy.
        - temperature: The temperature (in K) used for thermodynamic analysis.
        - pressure: The pressure (in bar) used for thermodynamic analysis.
        - do_raman: Calculate Raman activity
        - do_ir: Calculate IR intensities
    """

    def __init__(self, hess_drv, comm=None, ostream=None):
        """
        Initializes vibrational analysis driver.
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

        self.hessian_driver = hess_drv
        self.scf_driver = hess_drv.scf_driver

        self.hessian = None
        self.mass_weighted_hessian = None # FIXME unused variable
        self.reduced_masses = None
        self.force_constants = None
        self.vib_frequencies = None

        self.normal_modes = None  # normalized, not mass-weighted
        self.cart_normal_modes = None  # neither normalized nor mass-weighted
        self.dipole_gradient = None
        self.ir_intensities = None
        self.polarizability_gradient = None
        self.raman_intensities = None

        self.flag = 'Vibrational Analysis Driver'

        # flag for numerical Hessian and pol. gradient
        self.numerical_hessian = False
        self.numerical_polgrad = False

        # flag for two-point or four-point approximation
        self.do_four_point_hessian = False
        self.do_four_point_polgrad = False

        self.do_ir = True
        self.do_raman = False
        self.frequencies = (0,)

        # flag for printing
        self.do_print_hessian = False
        self.do_print_polgrad = False
        self.print_depolarization_ratio = False

        # Thermodynamics
        self.elec_energy = 0.0
        self.temperature = 298.15
        self.pressure = 1.0

        self._input_keywords = {
            'vibrational': {
                'numerical_hessian': ('bool', 'do numerical hessian'),
                'numerical_polgrad': ('bool', 'do numerical polarizability gradient'),
                'do_four_point_hessian': ('bool', 'do four-point numerical integration'),
                'do_four_point_polgrad': ('bool', 'do four-point numerical integration'),
                'do_raman': ('bool', 'whether to calculate Raman activity'),
                'do_ir': ('bool', 'whether to calculate IR intensities'),
                'do_print_hessian': ('bool', 'whether to print the Hessian'),
                'do_print_polgrad': ('bool', 'whether to print the pol. gradient'),
                'print_depolarization_ratio': ('bool', 'whether to print Raman depolarization ratio'),
                'temperature': ('float', 'the temperature'),
                'pressure': ('float', 'the pressure'),
                'frequencies': ('seq_range', 'frequencies of external electric field'),
                },
            }

    def update_settings(self, method_dict=None, vib_dict=None, hess_dict=None, cphf_dict=None,
                        rsp_dict=None, polgrad_dict=None):
        """
        Updates settings in HessianDriver.

        :param method_dict:
            The input dictionary of method settings group.
        :param hess_dict:
            The input dictionary of Hessian settings group.
        :param cphf_dict:
            The input dictionary of CPHF (orbital response) settings.
        :param rsp_dict:
            The input dictionary for linear response settings
            (needed to compute the polarizability gradient).
        :param polgrad_dict:
            The input dictionary for the polarizability gradient
            (needed to compute Raman activity).
        """

        if method_dict is None:
            method_dict = {}
        if vib_dict is None:
            vib_dict = {}

        vib_keywords = {
            key: val[0] for key, val in
            self._input_keywords['vibrational'].items()
        }

        parse_input(self, vib_keywords, vib_dict)

        # Settings for property modules
        if hess_dict is None:
            hess_dict = {}
        if cphf_dict is None:
            cphf_dict = {}
        if rsp_dict is None:
            rsp_dict = {}
        if polgrad_dict is None:
            polgrad_dict = {}

        self.method_dict = dict(method_dict)
        self.vib_dict = dict(vib_dict)
        self.hess_dict = dict(hess_dict)
        self.cphf_dict = dict(cphf_dict)
        self.rsp_dict = dict(rsp_dict)
        self.polgrad_dict = dict(polgrad_dict)

    def compute(self, molecule, ao_basis=None, min_basis=None):
        """
        Performs vibrational analysis (frequencies and normal modes)
        based on the molecular Hessian employing the geomeTRIC module:
        J. Chem. Phys. 144, 214108, DOI: 10.1063/1.4952956

        :param molecule:
            The molecule.
        :param ao_basis:
            The AO basis set.
        :param min_basis:
            The minimal AO basis set.
        """

        # Compute the Hessian
        self.compute_hessian(molecule, ao_basis)

        # Compute the polarizability gradient for Raman intensities
        if self.do_raman:
            self.compute_polarizability_gradient(molecule, ao_basis)

        if self.rank == mpi_master():
            # Get vibrational frequencies and normal modes
            # TODO what to do about that filename input
            self.frequency_analysis(molecule, filename=None)

            # Calculate force constants
            self.force_constants = self.calculate_force_constant()

            # Calculate the gradient of the dipole moment for IR intensities
            if self.do_ir:
                self.ir_intensities = self.calculate_ir_intensity(self.normal_modes)

            # Calculate the analytical polarizability gradient for Raman intensities
            if self.do_raman:
                #self.compute_polarizability_gradient(molecule, ao_basis)
                self.raman_intensities, depol_ratio = self.calculate_raman_activity(
                        self.normal_modes)

            # Print the vibrational properties
            self.print_vibrational_analysis(molecule)

        return

    def print_vibrational_analysis(self, molecule, filename=None, rsp_drv=None):
        """
        Prints the results from the vibrational analysis.

        :param molecule:
            The molecule.
        :param filename:
            Filename where thermodynamic properties are saved by geomeTRIC.
        :param rsp_drv:
            The response driver (for excited state vibrational analysis).
        """

        # number of atoms, elements, and coordinates
        natm = molecule.number_of_atoms()
        elem = molecule.get_labels()
        number_of_modes = len(self.vib_frequencies)

        # Normalize the normal modes -- as done in geomeTRIC
        self.normal_modes /= np.linalg.norm(self.normal_modes,
                                            axis=1)[:, np.newaxis]

        title = 'Vibrational Analysis'
        self.ostream.print_header(title)
        self.ostream.print_header('=' * (len(title) + 2))
        self.ostream.print_blank()

        width = 52
        for k in range(number_of_modes):

            # Print indices and vibrational frequencies:
            index_string = '{:22s}{:d}'.format('Vibrational Mode', k + 1)
            self.ostream.print_header(index_string.ljust(width))
            self.ostream.print_header('-' * width)

            freq_string = '{:22s}{:20.2f}  {:8s}'.format(
                'Harmonic frequency:', self.vib_frequencies[k], 'cm**-1')
            self.ostream.print_header(freq_string.ljust(width))

            mass_string = '{:22s}{:20.4f}  {:8s}'.format(
                'Reduced mass:', self.reduced_masses[k], 'amu')
            self.ostream.print_header(mass_string.ljust(width))

            force_cnst_string = '{:22s}{:20.4f}  {:8s}'.format(
                'Force constant:', self.force_constants[k], 'mdyne/A')
            self.ostream.print_header(force_cnst_string.ljust(width))

            if self.ir_intensities is not None:
                ir_intens_string = '{:22s}{:20.4f}  {:8s}'.format(
                    'IR intensity:', self.ir_intensities[k], 'km/mol')
                self.ostream.print_header(ir_intens_string.ljust(width))

            if self.raman_intensities is not None:
                freq_unit = ' a.u.'
                freqs = list(self.raman_intensities.keys())
                for freq in freqs:
                    if freq == 0.0:
                        this_freq = 'static'
                    else:
                        this_freq = str(round(freq,4)) + freq_unit
                    raman_intens_string = '{:16s} {:12s} {:12.4f}  {:8s}'.format(
                            'Raman activity:', this_freq, self.raman_intensities[freq][k],
                            'A**4/amu')
                    self.ostream.print_header(raman_intens_string.ljust(width))

                if self.print_depolarization_ratio:
                    raman_parallel_str = '{:22s}{:20.4f}  {:8s}'.format(
                        'Parallel Raman:', int_pol[k], 'A**4/amu')
                    self.ostream.print_header(
                        raman_parallel_str.ljust(width))

                    raman_perpendicular_str = '{:22s}{:20.4f}  {:8s}'.format(
                        'Perpendicular Raman:', int_depol[k], 'A**4/amu')
                    self.ostream.print_header(
                        raman_perpendicular_str.ljust(width))

                    depolarization_str = '{:22s}{:20.4f}'.format(
                        'Depolarization ratio:', depol_ratio[k])
                    self.ostream.print_header(
                        depolarization_str.ljust(width))

            normal_mode_string = '{:22s}'.format('Normal mode:')
            self.ostream.print_header(normal_mode_string.ljust(width))

            normal_mode_string = '{:16s}{:>12s}{:>12s}{:>12s}'.format(
                '', 'X', 'Y', 'Z')
            self.ostream.print_header(normal_mode_string.ljust(width))

            # Print normal modes:
            for atom_index in range(natm):
                valstr = '{:<8d}'.format(atom_index + 1)
                valstr += '{:<8s}'.format(elem[atom_index])
                valstr += '{:12.4f}'.format(
                    self.normal_modes[k][atom_index * 3 + 0])
                valstr += '{:12.4f}'.format(
                    self.normal_modes[k][atom_index * 3 + 1])
                valstr += '{:12.4f}'.format(
                    self.normal_modes[k][atom_index * 3 + 2])
                self.ostream.print_header(valstr.ljust(width))

            self.ostream.print_blank()
            self.ostream.print_blank()

        self.ostream.flush()

    def frequency_analysis(self, molecule, filename=None):
        """
        Runs the frequency analysis from geomeTRIC to obtain vibrational
        frequencies and normal modes.

        :param molecule:
            The molecule.
        """

        err_msg = ('The installed geometric package does not support\n' +
                   '  vibrational analysis. Please install the latest\n' +
                   '  geometric via pip or conda.\n')
        assert_msg_critical(hasattr(geometric, 'normal_modes'), err_msg)

        # number of atoms, elements, and coordinates
        natm = molecule.number_of_atoms()
        elem = molecule.get_labels()
        coords = molecule.get_coordinates_in_bohr().reshape(natm * 3)

        # Save the electronic energy
        self.elec_energy = self.scf_driver.get_scf_energy()

        self.vib_frequencies, self.normal_modes, gibbs_energy = (
            geometric.normal_modes.frequency_analysis(
                coords,
                self.hessian,
                elem,
                energy=self.elec_energy,
                temperature=self.temperature,
                pressure=self.pressure,
                outfnm=filename,
                normalized=False))

    def calculate_raman_activity(self, normal_modes):
        """
        Calculates the Raman activity. Static or dynamic (resonant).

        :param normal_modes:
            The vibrational normal modes.

        :return:
            The raman activity.
        """

        # Get information from polarizability gradient dictionary
        raman_conversion_factor = self.get_conversion_factor('raman')

        # frequency of electric field
        freqs = list(self.polarizability_gradient.keys())
        nfreq = len(freqs)

        number_of_modes = len(self.vib_frequencies)

        # dictionary for Raman intensities
        raman_intensities = {}
        depol_ratio = None

        for freq in freqs:
            # get gradient for current frequency
            current_polarizability_gradient = self.polarizability_gradient[freq]
            size_x = current_polarizability_gradient.shape[0]
            size_y = current_polarizability_gradient.shape[1]
            size_k = self.normal_modes.shape[0]

            # einsum 'xyi,ik->xyk'
            raman_transmom = np.matmul(
                current_polarizability_gradient.reshape(size_x * size_y, -1),
                normal_modes.T).reshape(size_x, size_y, size_k)

            # Calculate rotational invariants
            alpha_bar = np.zeros((number_of_modes))
            gamma_bar_sq = np.zeros((number_of_modes))
            for i in range(3):
                alpha_bar += raman_transmom[i, i] / 3
                for j in range(i + 1, 3):
                    gamma_bar_sq += (0.5 * (raman_transmom[i, i] - raman_transmom[j, j])**2
                                     + 3 * raman_transmom[i, j]**2)

            alpha_bar_sq = alpha_bar**2
            raman_intensities[freq] = (
                45 * alpha_bar_sq + 7 * gamma_bar_sq) * raman_conversion_factor

            if self.print_depolarization_ratio and (freq == 0.0): # TODO dynamic also?
                int_pol = 45 * alpha_bar_sq + 4 * gamma_bar_sq
                int_depol = 3 * gamma_bar_sq
                depol_ratio = int_depol / int_pol

        return raman_intensities, depol_ratio

    def calculate_ir_intensity(self, normal_modes):
        """
        Calculates the IR intensity of the normal modes.

        :param normal_modes:
            The vibrational normal modes.

        :return:
            The IR intensities.
        """

        conv_ir_ea0amu2kmmol = self.get_conversion_factor('ir')

        ir_trans_dipole = self.dipole_gradient.dot(normal_modes.T)
        ir_intensity_au_amu = np.array([
            np.linalg.norm(ir_trans_dipole[:, x])**2
            for x in range(ir_trans_dipole.shape[1])
        ])

        return ir_intensity_au_amu * conv_ir_ea0amu2kmmol

    def calculate_force_constant(self):
        """
        Calculates force constants.

        :return:
            Force constants of vibrational normal modes.
        """

        # Constants and conversion factors
        c = speed_of_light_in_vacuum_in_SI()
        cm_to_m = 1e-2  # centimeters in meters
        N_to_mdyne = 1e+8  # Newton in milli dyne
        m_to_A = 1e+10  # meters in Angstroms

        # Diagonalizes Hessian and calculates the reduced masses
        # einsum 'ki->i'
        # TODO is this a good place for reduced masses?
        self.reduced_masses = 1.0 / np.sum(self.normal_modes.T**2, axis=0)

        force_constants = (4.0 * np.pi**2 *
                                (c * (self.vib_frequencies / cm_to_m))**2 *
                                self.reduced_masses *
                                amu_in_kg()) * (N_to_mdyne / m_to_A)
        return force_constants

    def get_conversion_factor(self, prop):
        """
        Calculates conversion factor for IR or Raman.

        :param prop:
            Which property conversion factor.

        :return:
            Conversion factor.

        """

        # Constants and conversion factors
        alpha = fine_structure_constant()
        bohr_in_km = bohr_in_angstrom() * 1e-13

        # Conversion factor of IR intensity to km/mol
        if (prop == 'ir'):
            conv_factor = conv_ir_ea0amu2kmmol = (electron_mass_in_amu() * avogadro_constant()
                                    * alpha**2 * bohr_in_km * np.pi / 3.0)
        elif (prop == 'raman'):
            conv_factor = 0.078424
            
        return conv_factor

    def compute_hessian(self, molecule, ao_basis):
        """
        Directs the calculation of the nuclear Hessian

        :param molecule:
            The molecule.
        :param ao_basis:
            The AO basis set.
        """

        hessian_drv = self.hessian_driver
        hessian_drv.update_settings(self.method_dict, self.hess_dict, cphf_dict=self.cphf_dict)

        # Transfer settings for vibrational task to Hessian driver
        hessian_drv.numerical = self.numerical_hessian
        hessian_drv.do_four_point = self.do_four_point_hessian
        hessian_drv.do_dipole_gradient = self.do_ir
        hessian_drv.do_print_hessian = self.do_print_hessian

        hessian_drv.compute(molecule, ao_basis)

        self.hessian = hessian_drv.hessian
        self.dipole_gradient = hessian_drv.dipole_gradient

    def compute_polarizability_gradient(self, molecule, ao_basis):
        """
        Directs the calculation of the polarizability gradient
        needed for Raman activity.

        :param molecule:
            The molecule.
        :param ao_basis:
            The AO basis set.
        """

        scf_tensors = self.scf_driver.scf_tensors

        # Perform a linear response calculation
        lr_drv = LinearResponseSolver()
        lr_drv.update_settings(self.rsp_dict, self.method_dict)
        if 'frequencies' not in self.rsp_dict:
            lr_drv.frequencies = self.frequencies
        lr_results = lr_drv.compute(molecule, ao_basis, scf_tensors)

        # Set up the polarizability gradient driver and run
        polgrad_drv = PolarizabilityGradient(self.comm, self.ostream)
        polgrad_drv.update_settings(self.polgrad_dict,
                                    orbrsp_dict = self.cphf_dict,
                                    method_dict = self.method_dict, 
                                    scf_drv = self.scf_driver)

        # Transfer settings for vibrational task to polgrad driver
        polgrad_drv.numerical = self.numerical_polgrad
        polgrad_drv.do_four_point = self.do_four_point_polgrad
        polgrad_drv.do_print_polgrad = self.do_print_polgrad
        # TODO will we do it this way around?
        if 'frequencies' not in self.polgrad_dict:
            polgrad_drv.frequencies = self.frequencies

        polgrad_drv.compute(molecule, ao_basis, scf_tensors, lr_results)

        self.polarizability_gradient = polgrad_drv.polgradient

