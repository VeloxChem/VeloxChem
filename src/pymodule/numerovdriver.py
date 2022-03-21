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
import os
import sys

from scipy.optimize import curve_fit

from .scfunrestdriver import ScfUnrestrictedDriver
from .scffirstorderprop import ScfFirstOrderProperties
from .rspabsorption import Absorption
from .inputparser import parse_input, print_keywords, parse_seq_range
from .errorhandler import assert_msg_critical
from .outputstream import OutputStream

from .veloxchemlib import Molecule
from .veloxchemlib import (bohr_in_angstroms,
                           hartree_in_wavenumbers,
                           hartree_in_ev)
from .veloxchemlib import boltzmann_in_evperkelvin

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
        - work in progress...

    """

    def __init__(self, comm=None, ostream=None):
        """
        Initializes the Numerov driver.
        """

        # Potential energy curve (PEC) scan parameters
        #self.pec_start = 0.7
        #self.pec_end = 2.0
        #self.pec_step = 0.1
        self.pec_displacements = parse_seq_range('-0.7 - 2.0 (0.1)')
        self.pec_energies = None
        self.pec_properties = None

        # spectroscopy parameters
        self.reduced_mass = None
        self.eq_bond_len = None
        # electronic
        self.el_transition = False
        #self.initial_state = 0
        self.final_el_state = 1
        # vibrationial
        self.n_vib_states = 5
        # rotational
        self.rotation = True
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
        if comm is None:
            comm = MPI.COMM_WORLD
        self.comm = comm
        self.rank = self.comm.Get_rank()
        self.nodes = self.comm.Get_size()

        # output streams
        if ostream is None:
            ostream = OutputStream(sys.stdout)
        self.ostream = ostream
        self.silent = OutputStream(os.devnull)

        # input keywords
        self.input_keywords = {
            'numerov': {
                #'pec_start':
                #    ('float', 'start distance for PEC screening [bohr]'),
                #'pec_end':
                #   ('float', 'end distance for PEC screening [bohr]'),
                #'pec_step':
                #    ('float', 'step size for PEC screening [bohr]'),
                'pec_displacements':
                    ('seq_range', 'PEC screening range [bohr]'),
                'n_vib_states':
                    ('int', 'number of vibrational states to be resolved'),
                'reduced_mass':
                    ('float', 'reduced mass of the molecule'),
                'el_transition':
                    ('bool', 'include an electronic transition'),
                #'intial_state':
                #    ('int', 'initial state of the electronic transition'),
                'final_el_state':
                    ('int', 'final state of the electronic transition'),
                'rotation':
                    ('bool', 'enable/disable rotational resolution'),
                'temp':
                    ('float', 'temperature in Kelvin'),
                'conv_thresh':
                    ('float', 'convergence threshold for vibrational energies'),
                'max_iter':
                    ('int', 'maximum number of iteration per vibrational state'),
                'n_margin':
                    ('float', 'negative margin of the numerical grid [bohr]'),
                'p_margin':
                    ('float', 'positive margin of the numerical grid [bohr]'),
                'steps_per_au':
                    ('int', 'number of grid points per bohr radius'),
            }
        }

    def print_keywords(self):
        """
        Prints input keywords in Numerov driver.
        """

        print_keywords(self.input_keywords, self.ostream)

    def update_settings(self, numerov_dict, scf_dict, method_dict=None):
        """
        Updates settings in Numerov driver.

        :param numerov_dict:
            The dictionary of numerov input.
        :param method_dict:
            The dictionary of method settings.
        """

        numerov_keywords = {
            key: val[0] for key, val in self.input_keywords['numerov'].items()
        }

        parse_input(self, numerov_keywords, numerov_dict)

        if method_dict is not None:
            self.method_dict = dict(method_dict)

        if scf_dict is not None:
            self.scf_dict = dict(scf_dict)

    def compute(self, molecule, ao_basis, min_basis=None):
        """
        To be defined...
        """

        # what else is needed if the energy profiles are provided?


        if not self.pec_energies and not self.pec_properties:
            self.pec_energies, self.properties = self.generate_pec(molecule, ao_basis, min_basis)

        #print('ground state energies:', self.pec_energies['i'])
        #print('ground state dipole moments:', dipmoms)

        #print('excited state energies:', self.pec_energies['f'])
        #print('bond lengths:',bond_lengths)
        #print('displacements:',self.pec_displacements)


        # calculate vibronic wave functions using the Numerov method for PEC(s)

        i_pec = np.array(self.pec_energies['i']) - min(self.pec_energies['i'])

        # use morse potential: figure out a nicer way
        tstart = [1.e+3, 1, 3, 0]

        i_pec_morse_params, pcov = curve_fit(self.morse, self.pec_displacements, i_pec, p0 = tstart, maxfev=40000000)

        i_grid, i_psi, i_vib_energies = self.solve_numerov(i_pec_morse_params)

        # check for convergence
        assert_msg_critical(self.is_converged,
           'Vibronic wave functions for initial electronic state not converged')

        #print('energies:', i_vib_energies)

        if self.el_transition:
            f_pec = np.array(pec_energies['f']) - min(pec_energies['f'])

            f_pec_morse_params, pcov = curve_fit(self.morse, self.pec_displacements, f_pec, p0 = tstart, maxfev=40000000)

            f_grid, f_psi, f_vib_energies = self.solve_numerov(f_pec_morse_params)

            # check for convergence
            assert_msg_critical(self.is_converged,
             'Vibronic wave functions for final electronic state not converged')

            #print('final state energies:', f_vib_energies)

        # calculate spectra

        # discretize property curves
        x_prop_coefs = np.polyfit(self.pec_displacements, np.array(self.properties)[:,0], 6)
        y_prop_coefs = np.polyfit(self.pec_displacements, np.array(self.properties)[:,1], 6)
        z_prop_coefs = np.polyfit(self.pec_displacements, np.array(self.properties)[:,2], 6)

        dipmom_curves = {}
        dipmom_curves['x'] = np.polyval(x_prop_coefs, i_grid)
        dipmom_curves['y'] = np.polyval(y_prop_coefs, i_grid)
        dipmom_curves['z'] = np.polyval(z_prop_coefs, i_grid)


        # IR for a single state
        if not self.el_transition:
            spectrum = self.get_IR_spectrum(i_vib_energies, i_psi, dipmom_curves)

            return {
                'grid': i_grid,
                'psi': i_psi,
                'vib_levels': i_vib_energies,
                'excitation_energies': spectrum[0],
                'oscillator_strenths': spectrum[1],
            }

        # UV/Vis for two states
        else:
            UV_energies, UV_intensities = None, None

            return {
                'grid': i_grid,
                'f_grid': f_grid,
                'i_psi': i_psi,
                'i_vib_levels': i_vib_energies,
                'f_spi': f_psi,
                'f_vib_levels': f_vib_energies,
                'UV_energies': UV_energies,
                'UV_intensities': UV_intensities,
            }

    def generate_pec(self, molecule, ao_basis, min_basis=None):

        # sanity check
        assert_msg_critical(molecule.number_of_atoms() == 2,
                            'Only applicable to diatomic molecules')

        # get data for PEC screening
        mol_coords = molecule.get_coordinates()
        self.eq_bond_len = np.linalg.norm(mol_coords[0] - mol_coords[1])
        bond_lengths = self.pec_displacements + self.eq_bond_len

        atom1, atom2 = molecule.get_labels()

        #### extract reduced mass here with the help of the atom labels
        if not self.reduced_mass:
            self.reduced_mass = 0.9796 * 1822.888479031408

        # calculate PEC(s) for initial (i) and final (f) state
        pec_energies = {'i': [], 'f': []}
        dipoles = []

        f_minimum_found = False

        # initiate unrestricted SCF driver
        scf_drv = ScfUnrestrictedDriver(self.comm, self.silent)
        scf_drv.update_settings(self.scf_dict, self.method_dict)

        scf_prop = ScfFirstOrderProperties(self.comm, self.silent)

        # PEC scan
        for x in bond_lengths:
            geometry = Molecule.read_str(
                """{} 0  0 0
                   {} {} 0 0""".format(atom1, atom2, x * bohr_in_angstroms())
            )

            scf_drv.compute(geometry, ao_basis, min_basis)
            scf_prop.compute(geometry, ao_basis, scf_drv.scf_tensors)

            # save energies and dipole moments
            pec_energies['i'].append(scf_drv.iter_data[-1]['energy'])
            dipoles.append(scf_prop.get_property('dipole moment'))

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
                while correct_state_found == False:
                    # twice the number of excited states to consider
                    # for unrestricted case
                    rsp_drv = Absorption({'nstates': 2*excited_state - 1}, self.method_dict)
    
                    rsp_drv.init_driver(self.comm, self.silent)
                    rsp_drv.compute(geometry, ao_basis, scf_drv.scf_tensors)
    
                    total_energy = (rsp_drv.rsp_property['eigenvalues'][-1]
                                    + scf_drv.iter_data[-1]['energy'])
   
                    # detect PEC minimum
                    if prev_energy and f_minimum_found == False:
                        if prev_energy < total_energy:
                            minima_found[state] = True
                    # total energy expected to rise after the minimum
                    if f_minimum_found:
                        if prev_energy > total_energy:
                            # jumping up one state if not
                            excited_state += 1
                            continue
                    correct_state_found = True

                pec_energies['f'].append(total_energy)
                dipoles.append(rsp_drv.rsp_property['electric_transition_dipoles'][-1])

        if self.el_transition:
            assert_msg_critical(f_minimum_found,
                'No minimum was found for the final state PEC')

        return pec_energies, dipoles

    def solve_numerov(self, morse_params):

        # reset convergence flag
        self.is_converged = False

        # define grid
        interval = self.n_margin + -self.pec_displacements[0] + self.pec_displacements[-1] + self.p_margin
        n_steps = self.steps_per_au * interval

        n_grid = int(n_steps + 1)

        dx = interval / n_steps
        h2 = dx**2

        # discretize potential
        grid = np.zeros(n_grid)
        pec = np.zeros(n_grid)
        for i in range(n_grid):
            grid[i] = (self.pec_displacements[0] - self.n_margin) + i * dx
            pec[i] = self.morse(grid[i], morse_params[0], morse_params[1], morse_params[2], morse_params[3])

        # initialize arrays
        g = np.zeros(n_grid)
        psi_v = np.zeros(n_grid)
        psi = np.zeros([self.n_vib_states, n_grid])
        energies = np.zeros(self.n_vib_states)

        # initialize convergence variables
        e_guess = 1.0e-4
        e_step = 1.0e-4

        n_nodes_previous = 0

        while self.cur_iter < self.max_iter:

            self.cur_iter += 1

            # initialize node counter and last numerically relevant grid point
            n_nodes = 0
            i_last = n_grid

            # define function g for the Numerov procedure
            g = self.reduced_mass * 2.0 * (pec - e_guess)

            # manually set the first two data points (boundary condition)
            psi_v[0] = 0.0
            psi_v[1] = 1.0e-6

            # calculate all other data points with the Numerov procedure
            for i in range(2,n_grid):
                fac1 = 2.0 + h2 * 5.0 / 6.0 * g[i-1]
                fac2 = 1.0 - h2 * g[i-2] / 12.0

                numerator = psi_v[i-1] * fac1 - psi_v[i-2] * fac2
                denominator = 1.0 - g[i] * h2 / 12.0

                psi_v[i] = numerator / denominator

                # check for node
                if (psi_v[i] / psi_v[i-1]) < 0.0:
                    n_nodes += 1
                    i_last = i

            # remove numerical noise
            psi_v[i_last+1:] = 0.0

            if (abs(e_step) < self.conv_thresh) and (n_nodes > n_nodes_previous):
                # convergence of a vibrational state is reached
                self.print_iteration(n_nodes-1)

                # normalize the wave function
                psi_v /= np.sqrt(np.dot(psi_v,psi_v))

                # save psi and energy
                psi[n_nodes-1] = psi_v
                energies[n_nodes-1] = e_guess

                # reset energy step and iteration
                e_step = 1.0e-4
                n_nodes_previous += 1
                self.cur_iter = 0

                if n_nodes_previous == self.n_vib_states:
                    # all states are converged
                    self.is_converged = True
                    break

            # if energy guess too high/low and moving upwards/downwards
            if ((n_nodes > n_nodes_previous and e_step > 0.0) or
                (n_nodes <= n_nodes_previous and e_step < 0.0)):
                # reduce step size and change direction
                e_step /= -10.0

            # update the energy guess
            e_guess += e_step

        return grid, psi, energies

    def print_iteration(self, v):

        width = 92

        output_header = '*** Vibronic state {} converged '.format(v)
        output_header += 'in {} iterations'.format(self.cur_iter)
        self.ostream.print_header(output_header.ljust(width))
        self.ostream.flush()


    def morse(self, x, q, m, u, v):
        return (q * (np.exp(-2*m*(x-u))-2*np.exp(-m*(x-u))) + v)


    def get_IR_spectrum(self, vib_energies, psi, dmc):

        exc_energies = {}
        osc_str = {}

        # only transition 0->1 significant
        trans_moms = [np.dot(psi[1], dmc[i] * psi[0]) for i in 'xyz']
        vib_exc_energy = vib_energies[1] - vib_energies[0]
        vib_osc_str = 2.0 / 3.0 * vib_exc_energy * np.sum([tm**2 for tm in trans_moms])

        # (forbidden) Q branch
        exc_energies['Q'] = vib_exc_energy * hartree_in_wavenumbers()
        osc_str['Q'] = vib_osc_str

        # rotational resolution
        inertia = self.reduced_mass * self.eq_bond_len**2
        B = 1.0 / (2.0 * inertia)

        # R branch
        exc_energies['R'] = np.array([])
        osc_str['R'] = np.array([])

        for j in range(0,self.n_rot_states):
            exc_energies['R'] = np.append(exc_energies['R'], (vib_exc_energy + 2.0 * B * (j+1)) * hartree_in_wavenumbers())

            R_j = (2*j + 1) * np.exp((-B * j*(j+1)) / (boltzmann_in_evperkelvin() / hartree_in_ev() * self.temp))
            osc_str['R'] = np.append(osc_str['R'], R_j)

        Z_r = np.sum(osc_str['R'])
        osc_str['R'] *= 0.5 * vib_osc_str / Z_r

        # P branch
        exc_energies['P'] = np.array([])
        osc_str['P'] = np.array([])

        for j in range(1,self.n_rot_states):
            exc_energies['P'] = np.append(exc_energies['P'], (vib_exc_energy - 2.0 * B * j) * hartree_in_wavenumbers())

            P_j = (2*j + 1) * np.exp((-B * j*(j+1)) / (boltzmann_in_evperkelvin() / hartree_in_ev() * self.temp))
            osc_str['P'] = np.append(osc_str['P'], P_j)

        Z_p = np.sum(osc_str['P'])
        osc_str['P'] *= 0.5 * vib_osc_str / Z_p

        return (exc_energies, osc_str)

    def read_pec_data(self, bond_lengths, properties, *pec_energies):

        # sanity checks
        for i in [bond_lengths, properties, *pec_energies]:
            assert_msg_critical(isinstance(i, (list, np.ndarray)),
                'at least one input argument is not a list or array')

        assert_msg_critical(all([len(i) == len(displacements) for i in
                                [properties, *pec_energies]]),
            'input arrays are not of the same length')

        self.properties = properties
        self.pec_energies = {key: value for key, value in zip(['i', 'f'], [*pec_energies])}
        self.eq_bond_len = bond_lengths[np.argmin(self.pec_energies['i'])]
        self.pec_displacements = bond_lengths - self.eq_bond_len

    def set_reduced_mass(self, reduced_mass):

        self.reduced_mass = reduced_mass



