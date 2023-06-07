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
import numpy as np
import sys

from .veloxchemlib import Molecule
from .veloxchemlib import (mpi_master, bohr_in_angstroms,
                           hartree_in_wavenumbers, hartree_in_ev,
                           amu_in_electron_masses,
                           boltzmann_in_hartreeperkelvin)
from .outputstream import OutputStream
from .scfrestdriver import ScfRestrictedDriver
from .firstorderprop import FirstOrderProperties
from .rspabsorption import Absorption
from .errorhandler import assert_msg_critical
from .inputparser import parse_input, print_keywords
from .dftutils import get_default_grid_level


class NumerovDriver:
    """
    Implements the calculation of (ro-)vibronic spectra for diatomic
    molecules using the Numerov procedure for numerically solving the
    one-dimensional Schrödinger equation.

    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.

    Instance variables:
        - pec_displacements: The potential energy (PEC) scanning range.
        - pec_potential: The type of potential to fit the scan.
        - pec_data: The external PEC data flag.
        - pec_energies: The PEC energies.
        - pec_properties: The corresponding PEC properties.
        - reduced_mass: The reduced mass.
        - eq_bond_len: The ground state equilibrium bond length.
        - el_transition: The electronic transition flag.
        - final_el_state: The final state of the electronic transition.
        - exc_conv_thresh: The convergence threshold for the excited state calculation.
        - n_vib_states: The number of vibrational states per electronic state.
        - temp: The temperature.
        - n_rot_states: The number of rotational states.
        - conv_thresh: The Numerov convergence threshold.
        - max_iter: The maximum number of iterations per vibronic state.
        - cur_iter: The current Numerov iteration per vibronic state.
        - p_margin: The margin of the numerical grid right of the scanned PEC.
        - n_margin: The margin of the numerical grid left of the scanned PEC.
        - steps_per_au: The number of grid points per Bohr radius.
        - is_converged: The convergence flag for the Numerov solver.
        - comm: The MPI communicator.
        - rank: The rank of MPI process.
        - nodes: The number of MPI processes.
        - ostream: The output stream.
    """

    def __init__(self, comm=None, ostream=None):
        """
        Initializes the Numerov driver.
        """

        if comm is None:
            comm = MPI.COMM_WORLD

        if ostream is None:
            if comm.Get_rank() == mpi_master():
                ostream = OutputStream(sys.stdout)
            else:
                ostream = OutputStream(None)

        # Potential energy curve (PEC) scan parameters
        self.pec_displacements = list(np.arange(-0.7, 2.01, 0.1))
        self.pec_potential = 'morse'
        self.pec_data = False
        self.pec_energies = None
        self.pec_properties = None

        # spectroscopy parameters
        self.reduced_mass = None
        self.eq_bond_len = None
        # electronic
        self.el_transition = False
        self.final_el_state = 1
        self.exc_conv_thresh = 1.0e-4
        # vibrationial
        self.n_vib_states = 5
        # rotational
        self.temp = 298.15
        self.n_rot_states = 15

        # numerov solver setup
        self.conv_thresh = 1.0e-12
        self.max_iter = 1000
        self.cur_iter = 0
        self.p_margin = 1.5
        self.n_margin = 0.5
        self.steps_per_au = 500
        self.is_converged = False

        # mpi information
        self.comm = comm
        self.rank = self.comm.Get_rank()
        self.nodes = self.comm.Get_size()

        # output stream
        self.ostream = ostream

        # input keywords
        self.input_keywords = {
            'numerov': {
                'pec_displacements': ('seq_range', 'PEC scanning range [bohr]'),
                'pec_potential': ('str_lower', 'potential type for fitting'),
                'reduced_mass': ('float', 'reduced mass of the molecule'),
                'el_transition': ('bool', 'include an electronic transition'),
                'final_el_state':
                    ('int', 'final state of the electronic transition'),
                'exc_conv_thresh':
                    ('float', 'excited state calculation threshold'),
                'n_vib_states':
                    ('int', 'number of vibrational states to be resolved'),
                'temp': ('float', 'temperature in Kelvin'),
                'n_rot_states':
                    ('int', 'number of rotational states to be resolved'),
                'conv_thresh':
                    ('float', 'convergence threshold for vibronic energies'),
                'max_iter': ('int', 'max. iteration number per vibronic state'),
                'n_margin': ('float', 'negative margin of the grid [bohr]'),
                'p_margin': ('float', 'positive margin of the grid [bohr]'),
                'steps_per_au':
                    ('int', 'number of grid points per Bohr radius'),
            }
        }

    def print_keywords(self):
        """
        Prints input keywords in Numerov driver.
        """

        print_keywords(self.input_keywords, self.ostream)

    def update_settings(self, numerov_dict, scf_dict=None, method_dict=None):
        """
        Updates settings in Numerov driver.

        :param numerov_dict:
            The dictionary of numerov input.
        :param scf_dict:
            The dictionary of scf settings.
        :param method_dict:
            The dictionary of method settings.
        """

        numerov_keywords = {
            key: val[0] for key, val in self.input_keywords['numerov'].items()
        }

        parse_input(self, numerov_keywords, numerov_dict)

        if scf_dict is not None:
            self.scf_dict = dict(scf_dict)

        if method_dict is not None:
            self.method_dict = dict(method_dict)

    def compute(self, molecule=None, ao_basis=None, min_basis=None):
        """
        Computes vibronic wave functions and energy levels and the IR or UV/Vis
        spectrum for a diatomic molecule.

        :param molecule:
            The molecule.
        :param ao_basis:
            The AO basis set.
        :param min_basis:
            The minimal AO basis set.

        :return:
            The vibronic wave functions, the energy levels, the excitation
            energies and oscillator strengths.
        """

        try:
            from scipy.optimize import curve_fit
            from scipy import interpolate
        except ImportError:
            raise ImportError('Unable to import scipy. Please install scipy ' +
                              'via \'python3 -m pip install scipy\'')

        if not self.pec_data:
            # carry out PEC scan
            self.pec_energies, self.properties = self.generate_pec(
                molecule, ao_basis, min_basis)

        assert_msg_critical(
            self.reduced_mass,
            'NumerovDriver.compute: Reduced mass is not defined')

        # calculate vibronic wave functions for ground state PEC

        # shift data points
        i_pec = np.array(self.pec_energies['i']) - min(self.pec_energies['i'])

        # fit potential
        if self.pec_potential == 'morse':
            i_pec_pot_params, cv = curve_fit(self.morse, self.pec_displacements,
                                             i_pec)
        elif self.pec_potential == 'harmonic':
            i_pec_pot_params = np.polyfit(self.pec_displacements, i_pec, 2)

        # start numerov procedure
        self.print_numerov_header()
        i_grid, i_pot, i_psi, i_vib_levels = self.solve_numerov(
            i_pec_pot_params)

        # check for convergence
        assert_msg_critical(
            self.is_converged,
            'NumerovDriver.compute: Vibronic wave functions for ' +
            'initial electronic state not converged')
        self.print_numerov_convergence('initial')

        if self.el_transition:
            # calculate vibronic wave functions for excited state PEC

            # shift data points
            f_pec = (np.array(self.pec_energies['f']) -
                     min(self.pec_energies['f']))

            # fit potential
            if self.pec_potential == 'morse':
                f_pec_pot_params, cv = curve_fit(self.morse,
                                                 self.pec_displacements, f_pec)
            elif self.pec_potential == 'harmonic':
                f_pec_pot_params = np.polyfit(self.pec_displacements, f_pec, 2)
            else:
                assert_msg_critical(
                    self.pec_potential in ['morse', 'harmonic'],
                    'NumerovDriver.compute: Invalid potential shape {}'.format(
                        self.pec_potential))

            # start numerov procedure
            f_grid, f_pot, f_psi, f_vib_levels = self.solve_numerov(
                f_pec_pot_params)

            # check for convergence
            assert_msg_critical(
                self.is_converged,
                'NumerovDriver.compute: Vibronic wave functions for ' +
                'final electronic state not converged')
            self.print_numerov_convergence('final')

        # correct equilibrium bond length
        self.eq_bond_len += i_grid[np.argmin(i_pot)]

        # calculate spectra

        # discretize property curves
        prop_curves = {}
        if not self.el_transition:
            # use 6th degree polynomial for one elctronic state
            x_prop_coefs = np.polyfit(self.pec_displacements,
                                      np.array(self.properties)[:, 0], 6)
            y_prop_coefs = np.polyfit(self.pec_displacements,
                                      np.array(self.properties)[:, 1], 6)
            z_prop_coefs = np.polyfit(self.pec_displacements,
                                      np.array(self.properties)[:, 2], 6)

            prop_curves['x'] = np.polyval(x_prop_coefs, i_grid)
            prop_curves['y'] = np.polyval(y_prop_coefs, i_grid)
            prop_curves['z'] = np.polyval(z_prop_coefs, i_grid)
        else:
            # use isotropic average to identify outliers
            iso = np.array([np.sum(np.array(i)**2) for i in self.properties])

            # filter outliers potentially caused by conical intersections
            mask = self.mask_outliers(iso)

            displacements = np.ma.masked_array(self.pec_displacements,
                                               mask).compressed()
            x_prop = np.ma.masked_array(np.array(self.properties)[:, 0],
                                        mask).compressed()
            y_prop = np.ma.masked_array(np.array(self.properties)[:, 1],
                                        mask).compressed()
            z_prop = np.ma.masked_array(np.array(self.properties)[:, 2],
                                        mask).compressed()

            x_prop_coefs = interpolate.splrep(displacements, x_prop, s=0)
            prop_curves['x'] = interpolate.splev(i_grid, x_prop_coefs, der=0)
            y_prop_coefs = interpolate.splrep(displacements, y_prop, s=0)
            prop_curves['y'] = interpolate.splev(i_grid, y_prop_coefs, der=0)
            z_prop_coefs = interpolate.splrep(displacements, z_prop, s=0)
            prop_curves['z'] = interpolate.splev(i_grid, z_prop_coefs, der=0)

        # IR for a single electronic state
        if not self.el_transition:
            spectrum = self.get_IR_spectrum(i_psi, i_vib_levels, prop_curves)

            self.print_IR_spectrum(spectrum)

            return {
                'grid': i_grid,
                'potential': i_pot,
                'psi': i_psi,
                'vib_levels': i_vib_levels,
                'excitation_energies': spectrum[0],
                'oscillator_strengths': spectrum[1],
            }

        # UV/Vis for two electronic states
        else:
            spectrum = self.get_UV_spectrum(i_pot, i_psi, i_vib_levels, f_pot,
                                            f_psi, f_vib_levels, prop_curves)

            self.print_UV_spectrum(spectrum)

            return {
                'grid': i_grid,
                'i_potential': i_pot,
                'i_psi': i_psi,
                'i_vib_levels': i_vib_levels,
                'f_potential': f_pot,
                'f_psi': f_psi,
                'f_vib_levels': f_vib_levels,
                'excitation_energies': spectrum[0],
                'oscillator_strengths': spectrum[1],
            }

    def generate_pec(self, molecule, ao_basis, min_basis=None):
        """
        Carries out a potential energy curve scan of a diatomic molecule for the
        electronic ground state and an electronically excited state.

        :param molecule:
            The molecule.
        :param ao_basis:
            The AO basis set.
        :param min_basis:
            The minimal AO basis set.

        :return:
            The total energies and dipoles or transition moments.
        """

        # sanity checks
        assert_msg_critical(
            all([molecule, ao_basis]),
            'NumerovDriver.generate_pec: No molecule or basis set defined')
        assert_msg_critical(
            molecule.number_of_atoms() == 2,
            'NumerovDriver.generate_pec: Only applicable to diatomic molecules')

        # get data for PEC scan
        mol_coords = molecule.get_coordinates_in_bohr()
        self.eq_bond_len = np.linalg.norm(mol_coords[0] - mol_coords[1])
        bond_lengths = self.pec_displacements + self.eq_bond_len

        atom1, atom2 = molecule.get_labels()

        self.calculate_reduced_mass(molecule)

        # calculate PEC(s) for initial (i) and final (f) state
        pec_energies = {'i': [], 'f': []}
        props = []

        f_minimum_found = False

        # initiate unrestricted SCF driver
        scf_drv = ScfRestrictedDriver(self.comm, OutputStream(None))
        scf_drv.update_settings(self.scf_dict, self.method_dict)

        scf_prop = FirstOrderProperties(self.comm, OutputStream(None))

        # PEC scan
        self.print_PEC_header(scf_drv)

        for n, x in enumerate(bond_lengths):
            geometry = Molecule.read_str('{}  0 0 0\n{}  {} 0 0'.format(
                atom1, atom2, x * bohr_in_angstroms()))

            scf_drv.compute(geometry, ao_basis, min_basis)
            scf_prop.compute_scf_prop(geometry, ao_basis, scf_drv.scf_tensors)

            # save energies and dipole moments
            pec_energies['i'].append(scf_drv.scf_energy)
            if not self.el_transition:
                if self.rank == mpi_master():
                    dipole_moment = scf_prop.get_property('dipole moment')
                else:
                    dipole_moment = None
                dipole_moment = self.comm.bcast(dipole_moment,
                                                root=mpi_master())
                props.append(dipole_moment)

            # if an electronic transition is included
            if self.el_transition:
                # calculate the PEC for final state
                excited_state = self.final_el_state
                correct_state_found = False
                prev_energy = None
                if pec_energies['f']:
                    prev_energy = pec_energies['f'][-1]

                # try to find a smooth PEC by jumping up in states
                # after an avoided crossing
                while not correct_state_found:
                    # twice the number of excited states to consider
                    # for unrestricted case
                    rsp_prop = Absorption(
                        {
                            'nstates': 2 * excited_state - 1,
                            'conv_thresh': self.exc_conv_thresh
                        }, self.method_dict)

                    rsp_prop.init_driver(self.comm, OutputStream(None))
                    rsp_prop.compute(geometry, ao_basis, scf_drv.scf_tensors)

                    total_energy = (rsp_prop.rsp_property['eigenvalues'][-1] +
                                    scf_drv.scf_energy)

                    # detect PEC minimum
                    if prev_energy and not f_minimum_found:
                        if prev_energy < total_energy:
                            f_minimum_found = True
                    # total energy expected to rise after the minimum
                    if f_minimum_found:
                        if prev_energy > total_energy:
                            # jumping up one state if not
                            excited_state += 1
                            continue
                    correct_state_found = True

                # save energies and transition moments
                pec_energies['f'].append(total_energy)

                # assume degeneracy
                iso = 2.0 * np.sum(
                    rsp_prop.rsp_property['electric_transition_dipoles'][-1]**2)
                average = np.sqrt(iso / 3.0)
                props.append(np.array([average] * 3))

            self.print_PEC_iteration(n)

        if self.el_transition:
            assert_msg_critical(
                f_minimum_found,
                'NumerovDriver.generate_pec: No minimum was found for ' +
                'the final state PEC')

        self.print_PEC_convergence(bond_lengths, list(pec_energies.values()))

        return pec_energies, props

    def solve_numerov(self, pot_params):
        """
        Calculates the vibronic wave functions and energy levels for a given
        Morse potential using the numerical Numerov procedure.

        :param pot_params:
            The potential parameters.

        :return:
            The grid points, wave functions and energies.
        """

        # reset convergence flag
        self.is_converged = False

        # define grid
        interval = (self.n_margin + -self.pec_displacements[0] +
                    self.pec_displacements[-1] + self.p_margin)
        n_steps = self.steps_per_au * interval

        n_grid = int(n_steps + 1)

        dx = interval / n_steps
        h2 = dx**2

        # discretize potential
        grid = np.zeros(n_grid)
        pot = np.zeros(n_grid)
        for i in range(n_grid):
            grid[i] = (self.pec_displacements[0] - self.n_margin) + i * dx
            if self.pec_potential == 'morse':
                pot[i] = self.morse(grid[i], pot_params[0], pot_params[1],
                                    pot_params[2], pot_params[3])
            elif self.pec_potential == 'harmonic':
                pot[i] = np.polyval(pot_params, grid[i])

        # initialize arrays
        g = np.zeros(n_grid)
        psi_v = np.zeros(n_grid)
        psi = np.zeros([self.n_vib_states, n_grid])
        energies = np.zeros(self.n_vib_states)

        # initialize convergence variables
        e_guess = 1.0e-4
        e_step = 1.0e-4
        n_nodes_prev = 0

        while self.cur_iter < self.max_iter:
            self.cur_iter += 1

            # initialize node counter and last numerically relevant grid point
            n_nodes = 0
            i_last = n_grid

            # define function g for the Numerov procedure
            g = self.reduced_mass * 2.0 * (pot - e_guess)

            # manually set the first two data points (boundary condition)
            psi_v[0] = 0.0
            psi_v[1] = 1.0e-6

            # calculate all other data points with the Numerov procedure
            for i in range(2, n_grid):
                fac1 = 2.0 + h2 * 5.0 / 6.0 * g[i - 1]
                fac2 = 1.0 - h2 * g[i - 2] / 12.0

                numerator = psi_v[i - 1] * fac1 - psi_v[i - 2] * fac2
                denominator = 1.0 - g[i] * h2 / 12.0

                psi_v[i] = numerator / denominator

                # check for node
                if (psi_v[i] / psi_v[i - 1]) < 0.0:
                    n_nodes += 1
                    i_last = i

            # remove numerical noise
            psi_v[i_last + 1:] = 0.0

            if (abs(e_step) < self.conv_thresh) and (n_nodes > n_nodes_prev):
                # convergence of a vibrational state is reached
                self.print_numerov_iteration(n_nodes - 1)

                # normalize the wave function
                psi_v /= np.sqrt(np.dot(psi_v, psi_v))

                # save psi and energy
                psi[n_nodes - 1] = psi_v
                energies[n_nodes - 1] = e_guess

                # reset energy step and iteration
                e_step = 1.0e-4
                n_nodes_prev += 1
                self.cur_iter = 0

                if n_nodes_prev == self.n_vib_states:
                    # all states are converged
                    self.is_converged = True
                    break

            # if energy guess too high/low and moving upwards/downwards
            if (n_nodes > n_nodes_prev and
                    e_step > 0.0) or (n_nodes <= n_nodes_prev and e_step < 0.0):
                # reduce step size and change direction
                e_step /= -10.0

            # update the energy guess
            e_guess += e_step

        return grid, pot, psi, energies

    def morse(self, r, De, a, re, v):
        """
        Calculates a Morse potential.

        :param r:
            The r bond distance.
        :param De:
            The well depth.
        :param a:
            The potential width parameter.
        :param re:
            The equilibrium bond distance.
        :param v:
            The off-set parameter.
        """

        return (De * (np.exp(-2.0 * a * (r - re)) - 2.0 * np.exp(-a *
                                                                 (r - re))) + v)

    def mask_outliers(self, y, m=15.0):
        """
        Masks the outliers from an array.

        :param y:
            The data values.
        :param m:
            The cutoff value.

        :return:
            The inverted mask.
        """

        dev = np.abs(y - np.median(y))
        dev_med = np.median(dev)
        s = dev / dev_med if dev_med else 0.0

        return [s > m]

    def get_IR_spectrum(self, psi, vib_levels, dmc):
        """
        Calculates the rotationally resolved IR spectrum for the transition from
        the vibronic ground state to the first vibronically excited state.

        :param psi:
            The wave functions of the vibronic states.
        :param vib_levels:
            The energy levels of the vibronic states.
        :param dmc:
            The discretized dipole moment curve.

        :return:
            The excitation energies and oscillator strengths.
        """

        exc_energies = {}
        osc_str = {}

        # only transition 0->1 significant
        trans_moms = [np.dot(psi[1], dmc[i] * psi[0]) for i in 'xyz']
        vib_exc_energy = vib_levels[1] - vib_levels[0]
        vib_osc_str = (2.0 / 3.0 * vib_exc_energy *
                       np.sum([tm**2 for tm in trans_moms]))

        # (forbidden) Q branch
        exc_energies['Q'] = vib_exc_energy * hartree_in_wavenumbers()
        osc_str['Q'] = vib_osc_str

        # rotational resolution
        inertia = self.reduced_mass * self.eq_bond_len**2
        B = 1.0 / (2.0 * inertia)
        kB = boltzmann_in_hartreeperkelvin()

        # R branch
        exc_energies['R'] = np.array([])
        osc_str['R'] = np.array([])

        for j in range(0, self.n_rot_states):
            exc_energies['R'] = np.append(exc_energies['R'],
                                          (vib_exc_energy + 2.0 * B *
                                           (j + 1)) * hartree_in_wavenumbers())

            R_j = (2 * j + 1) * np.exp((-B * j * (j + 1)) / (kB * self.temp))
            osc_str['R'] = np.append(osc_str['R'], R_j)

        Z_r = np.sum(osc_str['R'])
        osc_str['R'] *= 0.5 * vib_osc_str / Z_r

        # P branch
        exc_energies['P'] = np.array([])
        osc_str['P'] = np.array([])

        for j in range(1, self.n_rot_states):
            exc_energies['P'] = np.append(exc_energies['P'],
                                          (vib_exc_energy - 2.0 * B * j) *
                                          hartree_in_wavenumbers())

            P_j = (2 * j + 1) * np.exp((-B * j * (j + 1)) / (kB * self.temp))
            osc_str['P'] = np.append(osc_str['P'], P_j)

        Z_p = np.sum(osc_str['P'])
        osc_str['P'] *= 0.5 * vib_osc_str / Z_p

        return (exc_energies, osc_str)

    def get_UV_spectrum(self, i_pot, i_psi, i_vib_levels, f_pot, f_psi,
                        f_vib_levels, tmc):
        """
        Calculates the UV/vis absorption and emission spectra for all calculated
        vibrational levels.

        :param i_pot:
            The Morse potential of the initial electronic state.
        :param i_psi:
            The vibronic wave functions of the initial electronic state.
        :param i_vib_levels:
            The vibronic energy levels of the initial electronic state.
        :param f_pot:
            The Morse potential of the final electronic state.
        :param f_psi:
            The vibronic wave functions of the final electronic state.
        :param f_vib_levels:
            The vibronic energy levels of the final electronic state.
        :param tmc:
            The discretized transition dipole moment curve.

        :return:
            The excitation energies and oscillator strengths.
        """

        exc_energies = {}
        osc_str = {}

        # potential minima difference
        f_pot_min = np.min(f_pot + np.min(self.pec_energies['f']))
        i_pot_min = np.min(i_pot + np.min(self.pec_energies['i']))
        pure_el_exc = f_pot_min - i_pot_min

        a = 'absorption'
        e = 'emission'

        exc_energies[a] = np.array([])
        osc_str[a] = np.array([])

        exc_energies[e] = np.array([])
        osc_str[e] = np.array([])

        for i in range(len(i_vib_levels)):
            # absorption spectrum
            abs_vib_trans_moms = [
                np.dot(i_psi[0], tmc[j] * f_psi[i]) for j in 'xyz'
            ]
            abs_energy = (f_vib_levels[i] + pure_el_exc) - i_vib_levels[0]
            abs_f = 2.0 / 3.0 * abs_energy * np.sum(
                [tm**2 for tm in abs_vib_trans_moms])

            exc_energies[a] = np.append(exc_energies[a],
                                        abs_energy * hartree_in_wavenumbers())
            osc_str[a] = np.append(osc_str[a], abs_f)

            # emission spectrum
            em_vib_trans_moms = [
                np.dot(f_psi[0], tmc[j] * i_psi[i]) for j in 'xyz'
            ]
            em_energy = (f_vib_levels[0] + pure_el_exc) - i_vib_levels[i]
            em_f = 2.0 / 3.0 * em_energy * np.sum(
                [tm**2 for tm in em_vib_trans_moms])

            exc_energies[e] = np.append(exc_energies[e],
                                        em_energy * hartree_in_wavenumbers())
            osc_str[e] = np.append(osc_str[e], em_f)

        return (exc_energies, osc_str)

    def read_pec_data(self, bond_lengths, properties, *pec_energies):
        """
        Reads potential energy curve data manually.

        :param bond_lengths:
            The bond lengths of the PEC scan in bohr radii.
        :param properties:
            The dipoles or transition dipole moments.
        :param pec_energies:
            The total or relative energies of the initial (and final) electronic
            state of the PEC scan.
        """

        # sanity checks
        for i in [bond_lengths, properties, *pec_energies]:
            assert_msg_critical(
                isinstance(i, (list, np.ndarray)),
                'NumerovDriver.read_pec_data: at least one input argument ' +
                'is not a list or array')

        assert_msg_critical(
            all([
                len(i) == len(bond_lengths)
                for i in [properties, *pec_energies]
            ]),
            'NumerovDriver.read_pec_data: input arrays are not of the same length'
        )

        self.properties = properties
        self.pec_energies = {
            key: value for key, value in zip(['i', 'f'], [*pec_energies])
        }
        self.eq_bond_len = (bond_lengths[np.argmin(self.pec_energies['i'])] /
                            bohr_in_angstroms())
        self.pec_displacements = (np.array(bond_lengths) / bohr_in_angstroms() -
                                  self.eq_bond_len)
        self.pec_data = True

    def calculate_reduced_mass(self, molecule):
        """
        Calculates the reduced mass in amu.

        :param molecule:
            The molecule.
        """

        m1, m2 = molecule.masses_to_numpy()
        mu = (m1 * m2) / (m1 + m2)

        self.set_reduced_mass(mu)

    def set_reduced_mass(self, reduced_mass):
        """
        Sets the reduced mass in a.u.

        :param reduced_mass:
            The reduced mass in amu.
        """

        self.reduced_mass = reduced_mass * amu_in_electron_masses()

    def print_PEC_header(self, scf_drv):
        """
        Prints the potential energy curve scan header to output stream.

        :param scf_drv:
            The SCF driver object.
        """

        title = 'Potential Energy Curve Driver Setup'
        self.ostream.print_blank()
        self.ostream.print_header(title)
        self.ostream.print_header('=' * (len(title) + 2))
        self.ostream.print_blank()

        str_width = 62

        # print SCF info
        cur_str = 'Number of Geometries          : '
        cur_str += str(len(self.pec_displacements))
        self.ostream.print_header(cur_str.ljust(str_width))
        cur_str = 'Wave Function Model           : '
        cur_str += scf_drv.get_scf_type_str()
        self.ostream.print_header(cur_str.ljust(str_width))
        cur_str = 'SCF Convergece Threshold      : {:.1e}'.format(
            scf_drv.conv_thresh)
        self.ostream.print_header(cur_str.ljust(str_width))

        if scf_drv._dft:
            cur_str = 'DFT Functional                : '
            cur_str += scf_drv.xcfun.get_func_label().upper()
            self.ostream.print_header(cur_str.ljust(str_width))
            grid_level = (get_default_grid_level(scf_drv.xcfun)
                          if scf_drv.grid_level is None else scf_drv.grid_level)
            cur_str = 'Molecular Grid Level          : ' + str(grid_level)
            self.ostream.print_header(cur_str.ljust(str_width))

        # print excited state info
        if self.el_transition:
            cur_str = 'Targeted Exicted State        : ' + str(
                self.final_el_state)
            self.ostream.print_header(cur_str.ljust(str_width))
            cur_str = 'Excited State Threshold       : {:.1e}'.format(
                self.exc_conv_thresh)
            self.ostream.print_header(cur_str.ljust(str_width))

        self.ostream.print_blank()
        self.ostream.flush()

    def print_PEC_iteration(self, n):
        """
        Prints convergence statement for a single potential energy curve data
        point to output stream.

        :param n:
            The index of the converged data point.
        """

        width = 92

        output_header = '*** Geometry {}/{}:'.format(
            n + 1, str(len(self.pec_displacements)))
        output_header += '  Initial state converged'
        if self.el_transition:
            output_header += ',  Final state converged'
        self.ostream.print_header(output_header.ljust(width))
        self.ostream.flush()

    def print_PEC_convergence(self, bond_lengths, pec_energies):
        """
        Prints potential energy curve data to output stream.

        :param bond_lengths:
            The bond lengths of the PEC scan.
        :param pec_energies:
            The dictionary of total energies of the intitial and final
            electronic state.
        """

        gs_energies = pec_energies[0] - np.min(pec_energies[0])
        if pec_energies[1]:
            es_energies = pec_energies[1] - np.min(pec_energies[0])
        else:
            es_energies = [None] * len(gs_energies)

        width = 92
        title = 'Potential Energy Curve Information'

        self.ostream.print_blank()
        self.ostream.print_header(title.ljust(width))
        self.ostream.print_header(('-' * len(title)).ljust(width))
        cur_str = ' ' * 30
        cur_str += '{:>16s}'.format('rel. GS energy')
        if pec_energies[1]:
            cur_str += '{:>16s}'.format('rel. ES energy')
        self.ostream.print_header(cur_str.ljust(width))
        for bl, gs_e, es_e in zip(bond_lengths, gs_energies, es_energies):
            cur_str = 'Bond distance: {:10.4f} angs '.format(
                bl * bohr_in_angstroms())
            cur_str += '{:12.5f} eV '.format(gs_e * hartree_in_ev())
            if es_e:
                cur_str += '{:12.5f} eV'.format(es_e * hartree_in_ev())
            self.ostream.print_header(cur_str.ljust(width))
        self.ostream.print_blank()
        self.ostream.flush()

    def print_numerov_header(self):
        """
        Prints the Numerov solver header to output stream.
        """

        title = 'Numerov Driver Setup'
        self.ostream.print_blank()
        self.ostream.print_header(title)
        self.ostream.print_header('=' * (len(title) + 2))
        self.ostream.print_blank()

        str_width = 40

        cur_str = 'Number of Vibrational States : ' + str(self.n_vib_states)
        self.ostream.print_header(cur_str.ljust(str_width))
        cur_str = 'Energy Convergence Threshold : {:.1e}'.format(
            self.conv_thresh)
        self.ostream.print_header(cur_str.ljust(str_width))
        cur_str = 'Max. Number of Iterations    : ' + str(self.max_iter)
        self.ostream.print_header(cur_str.ljust(str_width))
        cur_str = 'Grid Points per Bohr Radius  : ' + str(self.steps_per_au)
        self.ostream.print_header(cur_str.ljust(str_width))
        cur_str = 'Margin right of scanned PEC  : {:1.2f}'.format(self.p_margin)
        self.ostream.print_header(cur_str.ljust(str_width))
        cur_str = 'Margin left of scanned PEC   : {:1.2f}'.format(self.n_margin)
        self.ostream.print_header(cur_str.ljust(str_width))

        self.ostream.print_blank()
        self.ostream.flush()

    def print_numerov_iteration(self, v):
        """
        Prints convergence statement for a single vibrational state to output
        stream.

        :param v:
            The converged vibrational state.
        """

        width = 92

        output_header = '*** Vibronic state {} converged '.format(v)
        output_header += 'in {} iterations'.format(self.cur_iter)
        self.ostream.print_header(output_header.ljust(width))
        self.ostream.flush()

    def print_numerov_convergence(self, state):
        """
        Prints convergence statement for the Numerov procedure to output
        stream.

        :param state:
            The converged electronic state.
        """

        width = 92

        self.ostream.print_blank()
        output_header = '*** All vibronic states converged for '
        output_header += 'the {} electronic state'.format(state)
        self.ostream.print_header(output_header.ljust(width))
        self.ostream.print_blank()
        self.ostream.flush()

    def print_IR_spectrum(self, spectrum):
        """
        Prints the IR spectrum to output stream.

        :param spectrum:
            The spectrum results dictionary.
        """

        width = 92

        title = 'Rovibronic IR Spectrum'
        self.ostream.print_header(title.ljust(width))
        self.ostream.print_header(('-' * len(title)).ljust(width))

        cur_str = ' ' * 8
        cur_str += '{:>21s}'.format('Excitation energy')
        cur_str += '{:>23s}'.format('Oscillator strength')
        self.ostream.print_header(cur_str.ljust(width))

        cur_str = '{:>8s}'.format('P Branch')
        self.ostream.print_header(cur_str.ljust(width))
        for exc, osc in zip(spectrum[0]['P'], spectrum[1]['P']):
            cur_str = ' ' * 8
            cur_str += '{:16.2f} cm-1'.format(exc)
            cur_str += ' ' * 12
            cur_str += '{:.5e}'.format(osc)
            self.ostream.print_header(cur_str.ljust(width))

        cur_str = '{:>8s}'.format('Q Branch')
        self.ostream.print_header(cur_str.ljust(width))
        cur_str = ' ' * 8
        cur_str += '{:16.2f} cm-1'.format(spectrum[0]['Q'])
        cur_str += ' ' * 12
        cur_str += '{:.5e}'.format(spectrum[1]['Q'])
        self.ostream.print_header(cur_str.ljust(width))

        cur_str = '{:>8s}'.format('R Branch')
        self.ostream.print_header(cur_str.ljust(width))
        for exc, osc in zip(spectrum[0]['R'], spectrum[1]['R']):
            cur_str = ' ' * 8
            cur_str += '{:16.2f} cm-1'.format(exc)
            cur_str += ' ' * 12
            cur_str += '{:.5e}'.format(osc)
            self.ostream.print_header(cur_str.ljust(width))

        self.ostream.print_blank()
        self.ostream.flush()

    def print_UV_spectrum(self, spectrum):
        """
        Prints the UV/vis spectrum to output stream.

        :param spectrum:
            The spectrum results dictionary.
        """

        width = 92

        a = 'absorption'
        e = 'emission'

        title = 'UV/vis Absorption/Emission Spectrum'
        self.ostream.print_header(title.ljust(width))
        self.ostream.print_header(('-' * len(title)).ljust(width))

        cur_str = ' ' * 10
        cur_str += '{:>21s}'.format('Excitation energy')
        cur_str += '{:>23s}'.format('Oscillator strength')
        self.ostream.print_header(cur_str.ljust(width))

        cur_str = '{:>10s}'.format('Absorption')
        self.ostream.print_header(cur_str.ljust(width))
        for exc, osc in zip(spectrum[0][a], spectrum[1][a]):
            cur_str = ' ' * 10
            cur_str += '{:16.2f} cm-1'.format(exc)
            cur_str += ' ' * 12
            cur_str += '{:.5e}'.format(osc)
            self.ostream.print_header(cur_str.ljust(width))

        cur_str = '{:>8s}'.format('Emission')
        self.ostream.print_header(cur_str.ljust(width))
        for exc, osc in zip(spectrum[0][e], spectrum[1][e]):
            cur_str = ' ' * 10
            cur_str += '{:16.2f} cm-1'.format(exc)
            cur_str += ' ' * 12
            cur_str += '{:.5e}'.format(osc)
            self.ostream.print_header(cur_str.ljust(width))

        self.ostream.print_blank()
        self.ostream.flush()
