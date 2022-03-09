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

from .scfunrestdriver import ScfUnrestrictedDriver
from .scffirstorderprop import ScfFirstOrderProperties
from .rspabsorption import Absorption
from .inputparser import parse_input, print_keywords
from .errorhandler import assert_msg_critical

from .veloxchemlib import Molecule
from .veloxchemlib import bohr_in_angstroms

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

    def __init__(self, comm, ostream):
        """
        Initializes the Numerov driver.
        """

        # Potential energy curve (PEC) scan parameters
        self.pec_start = 0.7
        self.pec_end = 2.0
        self.pec_step = 0.1

        # spectroscopy parameters
        self.n_vib_states = 5

        self.vibronic = False
        self.exc_state = 1

        # solver setup
        self.conv_thresh = 1.0e-12
        self.max_iter = 1000
        self.cur_iter = 0
        self.displacement = 2.0
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
                'start_dist':
                    ('float', 'start distance for PEC screening in angstroms'),
                'end_dist':
                   ('float', 'end distance for PEC screening in angstroms'),
                'step_size':
                    ('float', 'step size for PEC screening in angstroms'),
                'n_vib_states':
                    ('int', 'number of vibrational states to be resolved'),
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


        # calculate total energies and dipoles
        gs_energies = []
        exc_energies = []

        gs_dipoles = []

        # initiate SCF driver
        scf_drv = ScfUnrestrictedDriver(self.comm, self.ostream)
        scf_drv.update_settings(self.scf_dict, self.method_dict)

        scf_prop = ScfFirstOrderProperties(self.comm, self.ostream)

        # initiate response driver
        if self.vibronic:
            rsp_drv = Absorption({'nstates': self.exc_state}, self.method_dict)

        # calculate PEC(s)
        for x in bond_lengths:
            geometry = Molecule.read_str(
                """{} 0  0 0
                   {} {} 0 0""".format(atom1, atom2, x * bohr_in_angstroms())
            )

            scf_drv.compute(geometry, ao_basis, min_basis)
            gs_energies.append(scf_drv.iter_data[-1]['energy'])

            scf_prop.compute(geometry, ao_basis, scf_drv.scf_tensors)
            gs_dipoles.append(scf_prop.get_property('dipole moment')[0])

            if self.vibronic:
                rsp_drv.init_driver(self.comm, self.ostream)
                rsp_drv.compute(geometry, ao_basis, scf_drv.scf_tensors)

                print(rsp_drv.rsp_property['eigenvalues'])
                exc_energies.append(rsp_drv.rsp_property['eigenvalues'][0]
                                    + scf_drv.iter_data[-1]['energy'])

        print('ground state energies:', gs_energies)
        print('ground state dipole moments:', gs_dipoles)

        ## after calculation of energy profiles:

        # calculating the vibrational wave functions using Numerov for every potential

        # setting up the input parameters:

        print('bond lengths:',bond_lengths)
        print('displacements:',displacements)

        gs_pec = np.array(gs_energies) - min(gs_energies)
        #print(gs_pec)
        gs_pec_poly_coefs = np.polyfit(displacements, gs_pec, 6)

        gs_dip_poly_coefs = np.polyfit(displacements, gs_dipoles, 6)

        print('pec coefficients:',gs_pec_poly_coefs)
        print('dip coefficients:',gs_dip_poly_coefs)

        if self.vibronic:
            es_pec = np.array(es_energies) - min(es_energies)
            es_pec_poly_coefs = np.polyfit(displacements, es_pec, 6)

        gs_n_grid, gs_psi, gs_psi_sq, gs_vib_energies, gs_tran_moms = self.solve_numerov(gs_pec_poly_coefs, gs_dip_poly_coefs)

        print('energies:',gs_vib_energies)

        print('transition moments:', gs_tran_moms)


        print('first_transition:', gs_vib_energies[1] - gs_vib_energies[0])



        # To-Do:
        # -sanity check if minimum lies in defined interval (after calculation)


        print('running through')


    def solve_numerov(self, pec_poly_coefs, dip_poly_coefs):

        # define grid
        dis = self.displacement
        n_steps = self.steps_per_au * 2.0 * self.displacement

        n_grid = int(n_steps + 1)

        dx = 2.0 * dis / n_steps
        h2 = dx**2

        # discretize potential
        grid = np.zeros(n_grid)
        pec = np.zeros(n_grid)
        dipmom = np.zeros(n_grid)
        for i in range(n_grid):
            grid[i] = -dis + i * dx
            pec[i] = np.polyval(pec_poly_coefs, grid[i])
            dipmom[i] = np.polyval(dip_poly_coefs, grid[i])

        g = np.zeros(n_grid)
        psi_v = np.zeros(n_grid)
        psi = np.zeros([self.n_vib_states, n_grid])
        psi_sq = np.zeros([self.n_vib_states, n_grid])
        energies = np.zeros(self.n_vib_states)
        tran_moms = np.zeros(self.n_vib_states)

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
                print('state {} converged'.format(n_nodes))

                # normalize the wave function
                psi_v /= np.sqrt(np.dot(psi_v,psi_v))

                psi[n_nodes-1] = psi_v
                psi_sq[n_nodes-1] = psi_v * psi_v
                energies[n_nodes-1] = e_guess
                tran_moms[n_nodes-1] = np.dot(psi_v, dipmom * psi_v)

                # reset energy step and iteration
                e_step = 1.0e-4
                n_nodes_previous += 1
                self.cur_iter = 0

                if n_nodes_previous == self.n_vib_states:
                    self.converged = True
                    break

            # if energy guess too high/low and moving upwards/downwards
            if ((n_nodes > n_nodes_previous and e_step > 0.0) or
                (n_nodes <= n_nodes_previous and e_step < 0.0)):
                # reduce step size and change direction
                e_step /= -10.0

            # update the energy guess
            e_guess += e_step


        return n_grid, psi, psi_sq, energies, tran_moms



