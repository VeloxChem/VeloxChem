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
import os

from scipy.optimize import curve_fit

from .scfunrestdriver import ScfUnrestrictedDriver
from .scffirstorderprop import ScfFirstOrderProperties
from .rspabsorption import Absorption
from .inputparser import parse_input, print_keywords
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
        self.pec_start = 0.7
        self.pec_end = 2.0
        self.pec_step = 0.1

        # spectroscopy parameters
        self.n_vib_states = 5

        self.electronic_excitation = False
        self.initial_state = 0
        self.final_state = 1

        self.rotation = True
        self.temp = 298.15
        self.n_rot_states = 15

        # solver setup
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

        # output streams
        self.ostream = ostream
        self.silent = OutputStream(os.devnull)

        # input keywords
        self.input_keywords = {
            'numerov': {
                'pec_start':
                    ('float', 'start distance for PEC screening [bohr]'),
                'pec_end':
                   ('float', 'end distance for PEC screening [bohr]'),
                'pec_step':
                    ('float', 'step size for PEC screening [bohr]'),
                'n_vib_states':
                    ('int', 'number of vibrational states to be resolved'),
                'electronic_excitation':
                    ('bool', 'include an electronic excitation'),
                'intial_state':
                    ('int', 'initial state of the electronic transition'),
                'final_state':
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

        if 'electronic_exciation' in numerov_keywords:
            if self.electronic_exciation == 'yes':
                self.electronic_excitation = True

        if method_dict is not None:
            self.method_dict = dict(method_dict)

        if scf_dict is not None:
            self.scf_dict = dict(scf_dict)

    def compute(self, molecule, ao_basis, min_basis=None):
        """
        To be defined...
        """

        # sanity check
        assert_msg_critical(molecule.number_of_atoms() == 2,
                            'Only applicable to diatomic molecules')

        # get data for PEC screening
        mol_coords = molecule.get_coordinates()
        eq_bond_len = np.linalg.norm(mol_coords[0] - mol_coords[1])

        bond_lengths = np.arange(eq_bond_len - self.pec_start,
                                 eq_bond_len + self.pec_end + self.pec_step,
                                 self.pec_step)
        displacements = bond_lengths - eq_bond_len

        atom1, atom2 = molecule.get_labels()

        # extract reduced mass here with the help of the atom labels
        self.reduced_mass = 0.9796 * 1822.888479031408

        # calculate PEC(s) for initial (i) and final (f) state
        pec_energies = {'i': [], 'f': []}

        dipmoms = []
        transmoms = []

        minima_found = {'i': False, 'f': False}

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

            # if initial (or only) state is the electronic ground state
            if self.initial_state == 0:
                # save energies and dipole moments
                pec_energies['i'].append(scf_drv.iter_data[-1]['energy'])
                dipmoms.append(scf_prop.get_property('dipole moment')[0])

            # if an electronic excitation is included
            if self.electronic_excitation:
                # calculate the PEC for initial (if not GS) and final state 
                for state, n in {'i': self.initial_state, 'f': self.final_state}.items():
                    if n > 0:
                        excited_state = n
                        correct_state_found = False
                        prev_energy = None
                        if pec_energies[state]:
                            prev_energy = pec_energies[state][-1]

                        # try to find a smooth PEC by jumping up in states
                        # after an ISC
                        while correct_state_found == False:
                            # twice the number of excited states to consider
                            # for unrestricted case
                            rsp_drv = Absorption({'nstates': 2*excited_state - 1}, self.method_dict)
    
                            rsp_drv.init_driver(self.comm, self.silent)
                            rsp_drv.compute(geometry, ao_basis, scf_drv.scf_tensors)
    
                            total_energy = (rsp_drv.rsp_property['eigenvalues'][-1]
                                            + scf_drv.iter_data[-1]['energy'])
                            #print('exc_energies:', rsp_drv.rsp_property['eigenvalues'])
                            #print('elec_trans_dipoles:',rsp_drv.rsp_property['electric_transition_dipoles'])
   
                            # detect PEC minimum
                            if prev_energy and minima_found[state] == False:
                                if prev_energy < total_energy:
                                    minima_found[state] = True
                            # total energy expected to rise right of the minimum
                            if minima_found[state]:
                                if prev_energy > total_energy:
                                    # jumping up one state if not
                                    excited_state += 1
                                    continue
                            correct_state_found = True

                        pec_energies[state].append(total_energy)
                        transmoms.append(rsp_drv.rsp_property['electric_transition_dipoles'])

        if self.initial_state != 0:
            assert_msg_critical(minima_found['i'],
                'No minimum was found for the initial state PEC')
        if self.electronic_excitation:
            assert_msg_critical(minima_found['f'],
                'No minimum was found for the final state PEC')

        #print('ground state energies:', i_pec_energies)
        #print('ground state dipole moments:', dipmoms)

        #print('excited state energies:', f_pec_energies)
        #print('bond lengths:',bond_lengths)
        #print('displacements:',displacements)


        # calculate vibronic wave functions using the Numerov method for PEC(s)

        i_pec = np.array(pec_energies['i']) - min(pec_energies['i'])

        # use morse potential: figure out a nicer way
        tstart = [1.e+3, 1, 3, 0]

        i_pec_morse_params, pcov = curve_fit(self.morse, displacements, i_pec, p0 = tstart, maxfev=40000000)

        i_grid, i_psi, i_vib_energies = self.solve_numerov(i_pec_morse_params)

        # check for convergence
        assert_msg_critical(self.is_converged,
            'Vibronic wave functions for initial electronic state not converged')

        print('energies:', i_vib_energies)

        if self.electronic_excitation:
            f_pec = np.array(pec_energies['f']) - min(pec_energies['f'])

            f_pec_morse_params, pcov = curve_fit(self.morse, displacements, f_pec, p0 = tstart, maxfev=40000000)

            f_grid, f_psi, f_vib_energies = self.solve_numerov(f_pec_morse_params)

            # check for convergence
            assert_msg_critical(self.is_converged,
                'Vibronic wave functions for final electronic state not converged')

            print('final state energies:', f_vib_energies)

        # calculate spectra

        # IR for a single state
        if self.electronic_excitation == False:
            # (so far) only implemented for the electronic ground state
            if self.initial_state == 0:
                i_dip_coefs = np.polyfit(displacements, dipmoms, 6)
                dmc = np.polyval(i_dip_coefs, i_grid)

                spectrum = self.get_spectrum(i_vib_energies, i_vib_energies, i_psi, i_psi, dmc, eq_bond_len)

            if self.rotation:
                return {
                    'grid': i_grid,
                    'psi': i_psi,
                    'vib_levels': i_vib_energies,
                    'IR_energies': spectrum[0],
                    'IR_intensities': spectrum[1],
                    'r_energies': spectrum[2],
                    'r_intensities': spectrum[3],
                    'p_energies': spectrum[4],
                    'p_intensities': spectrum[5],
                }
            else:
                return {
                    'grid': i_grid,
                    'psi': i_psi,
                    'vib_levels': i_vib_energies,
                    'IR_energies': spectrum[0],
                    'IR_intensities': spectrum[1],
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


        # To-Do:
        # Implement rotational resolution

        print('running through')




    def solve_numerov(self, morse_params):

        # reset convergence flag
        self.is_converged = False

        # define grid
        interval = self.n_margin + self.pec_start + self.pec_end + self.p_margin
        n_steps = self.steps_per_au * interval

        n_grid = int(n_steps + 1)

        dx = interval / n_steps
        h2 = dx**2

        # discretize potential
        grid = np.zeros(n_grid)
        pec = np.zeros(n_grid)
        for i in range(n_grid):
            grid[i] = -(self.pec_start + self.n_margin) + i * dx
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
                # write "print_iteration" function
                print('state {} converged'.format(n_nodes))

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


    def morse(self, x, q, m, u, v):
        return (q * (np.exp(-2*m*(x-u))-2*np.exp(-m*(x-u))) + v)


    def get_spectrum(self, i_vib_energies, f_rel_vib_energies, i_psi, f_psi, dmc, eq_bond_len):


        excitation_energies = {}
        oscillator_strengths = {}

        if self.rotation:
            r_energies = {}
            r_intensities = {}
            p_energies = {}
            p_intensities = {}

        for state in range(1,f_psi.shape[0]):

            transition = '0->{}'.format(state)

            trans_mom = np.dot(f_psi[state], dmc * i_psi[0])
            exc_energy = f_rel_vib_energies[state] - i_vib_energies[0]
            osc_str = 2.0 / 3.0 * exc_energy * trans_mom**2

            excitation_energies[transition] = exc_energy * hartree_in_wavenumbers()
            oscillator_strengths[transition] = osc_str

            if self.rotation:
                # rotational resolution
                inertia = self.reduced_mass * eq_bond_len**2
                B = 1.0 / (2.0 * inertia)

                # R branch
                r_energies[transition] = []
                r_intensities[transition] = []

                for j in range(0,self.n_rot_states):
                    r_energies[transition].append((exc_energy + 2.0 * B * (j+1)) * hartree_in_wavenumbers())

                    r_j = (2*j + 1) * np.exp((-B * j*(j+1)) / (boltzmann_in_evperkelvin() / hartree_in_ev() * self.temp))
                    r_intensities[transition].append(r_j)

                Z_r = np.sum(r_intensities[transition])
                r_intensities[transition] /= Z_r * 0.5 * osc_str

                # P branch
                p_energies[transition] = []
                p_intensities[transition] = []

                for j in range(1,self.n_rot_states):
                    p_energies[transition].append((exc_energy - 2.0 * B * j) * hartree_in_wavenumbers())

                    p_j = (2*j + 1) * np.exp((-B * j*(j+1)) / (boltzmann_in_evperkelvin() / hartree_in_ev() * self.temp))
                    p_intensities[transition].append(p_j)    

                Z_p = np.sum(p_intensities[transition])
                p_intensities[transition] /= Z_p * 0.5 * osc_str

        if self.rotation:
            return (excitation_energies, oscillator_strengths, r_energies,
                r_intensities, p_energies, p_intensities)

        return (excitation_energies, oscillator_strengths)


