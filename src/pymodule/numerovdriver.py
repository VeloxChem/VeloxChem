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

from .scfrestdriver import ScfRestrictedDriver
from .rspabsorption import Absorption
from .inputparser import parse_input, print_keywords
from .errorhandler import assert_msg_critical

from .veloxchemlib import Molecule

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
        self.start_dist = 1.0
        self.end_dist = 1.5
        self.step_size = 0.05

        # spectroscopy parameters
        self.n_vib_states = 5
        print('test')

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

        # get atom labels, PEC screening bond lenghts
        atom1, atom2 = molecule.get_labels()
        bond_lengths = np.arange(self.start_dist, self.end_dist + self.step_size,
                                 self.step_size)
        gs_energies = []
        exc_energies = []

        # initiate SCF driver
        scf_drv = ScfRestrictedDriver(self.comm, self.ostream)
        scf_drv.update_settings(self.scf_dict, self.method_dict)

        # initiate response driver
        rsp_drv = Absorption({'nstates': self.n_vib_states}, self.method_dict)

        # calculate PEC
        for x in bond_lengths:
            geometry = Molecule.read_str(
                """{} 0  0 0
                   {} {} 0 0""".format(atom1, atom2, x)
            )

            scf_drv.compute(geometry, ao_basis, min_basis)

            gs_energies.append(scf_drv.iter_data[-1]['energy'])

            rsp_drv.init_driver(self.comm, self.ostream)
            rsp_drv.compute(geometry, ao_basis, scf_drv.scf_tensors)

            print(rsp_drv.rsp_property['eigenvalues'])
            exc_energies.append(rsp_drv.rsp_property['eigenvalues'][0]
                                + scf_drv.iter_data[-1]['energy'])

        print('ground state energies:', gs_energies)
        print('excited state energies:', exc_energies)

        ## after calculation of energy profiles:

        # calculating the vibrational wave functions using Numerov


        # To-Do:
        # -sanity check if minimum lies in defined interval (after calculation)


        print('running through')






