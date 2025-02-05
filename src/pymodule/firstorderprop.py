#
#                              VELOXCHEM
#         ----------------------------------------------------
#                     An Electronic Structure Code
#
#  Copyright © 2018-2024 by VeloxChem developers. All rights reserved.
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

from .veloxchemlib import mpi_master
from .veloxchemlib import dipole_in_debye
from .outputstream import OutputStream
from .oneeints import compute_electric_dipole_integrals


class FirstOrderProperties:
    """
    Implements first-order properties.

    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.

    Instance variables:
        - properties: The dictionary of properties.
    """

    def __init__(self, comm=None, ostream=None):
        """
        Initializes SCF first-order properties.
        """

        if comm is None:
            comm = MPI.COMM_WORLD

        if ostream is None:
            if comm.Get_rank() == mpi_master():
                ostream = OutputStream(sys.stdout)
            else:
                ostream = OutputStream(None)

        self.comm = comm
        self.rank = comm.Get_rank()

        self.ostream = ostream

        self.properties = {}

    def compute_scf_prop(self, molecule, basis, scf_results):
        """
        Computes first-order properties for SCF.

        :param molecule:
            The molecule
        :param basis:
            The AO basis set.
        :param scf_results:
            The dictionary containing SCF results.
        """

        if self.rank == mpi_master():
            total_density = scf_results['D_alpha'] + scf_results['D_beta']
        else:
            total_density = None

        self.compute(molecule, basis, total_density)

    def compute(self, molecule, basis, total_density):
        """
        Computes first-order properties.

        :param molecule:
            The molecule
        :param basis:
            The AO basis set.
        :param total_density:
            The total electron density.
        """

        if molecule.get_charge() != 0:
            coords = molecule.get_coordinates_in_bohr()
            nuclear_charges = molecule.get_element_ids()
            origin = np.sum(coords.T * nuclear_charges,
                            axis=1) / np.sum(nuclear_charges)
        else:
            origin = np.zeros(3)

        dipole_moment = None
        if self.rank == mpi_master():
            dipole_ints = compute_electric_dipole_integrals(
                molecule, basis, origin)

            # electronic contribution

            # multiple states:
            if len(total_density.shape) > 2:
                dof = total_density.shape[0]
                electronic_dipole = np.zeros((dof, 3))
                for i in range(dof):
                    electronic_dipole[i] = np.array([
                        np.sum(dipole_ints[d] * total_density[i])
                        for d in range(3)
                    ])

            # or single state:
            else:
                electronic_dipole = np.array(
                    [np.sum(dipole_ints[d] * total_density) for d in range(3)])

            # nuclear contribution
            coords = molecule.get_coordinates_in_bohr()
            nuclear_charges = molecule.get_element_ids()
            nuclear_dipole = np.sum((coords - origin).T * nuclear_charges,
                                    axis=1)

            dipole_moment = nuclear_dipole + electronic_dipole
        dipole_moment = self.comm.bcast(dipole_moment, root=mpi_master())
        self.properties['dipole moment'] = dipole_moment
        self.properties['dipole_moment'] = dipole_moment

    def get_property(self, key):
        """
        Gets first-order property.

        :param key:
            The name of the property.

        :return:
            The property.
        """

        return self.properties[key]

    def print_properties(self, molecule, title=None, states=None):
        """
        Prints first-order properties.

        :param molecule:
            The molecule.
        :param title:
            The title of the printout, giving information about the method,
            relaxed/unrelaxed density, which excited state.
        :param states:
            The excited states for which first order properties ar calculated.
        """

        if title is None:
            title = "Ground State Dipole Moment"
        self.ostream.print_blank()

        self.ostream.print_header(title)
        self.ostream.print_header('-' * (len(title) + 2))

        # prints warning if the molecule is charged
        if molecule.get_charge() != 0:
            warn_msg = '*** Warning: Molecule has non-zero charge. Dipole'
            self.ostream.print_header(warn_msg.ljust(56))
            warn_msg = '    moment will be dependent on the choice of origin.'
            self.ostream.print_header(warn_msg.ljust(56))
            warn_msg = '    Center of nuclear charge is chosen as the origin.'
            self.ostream.print_header(warn_msg.ljust(56))

        self.ostream.print_blank()

        dip = self.properties['dipole_moment']
        if states is None:
            dip_au = list(dip) + [np.linalg.norm(dip)]
            dip_debye = [m * dipole_in_debye() for m in dip_au]

            for i, a in enumerate(['  X', '  Y', '  Z', 'Total']):
                valstr = '{:<5s} :'.format(a)
                valstr += '{:17.6f} a.u.'.format(dip_au[i])
                valstr += '{:17.6f} Debye   '.format(dip_debye[i])
                self.ostream.print_header(valstr)
        else:
            for index, s in enumerate(states):
                self.ostream.print_blank()
                state_text = 'Excited State %d' % s
                self.ostream.print_header(state_text)
                self.ostream.print_blank()
                dip_au = list(dip[index]) + [np.linalg.norm(dip[index])]
                dip_debye = [m * dipole_in_debye() for m in dip_au]

                for i, a in enumerate(['  X', '  Y', '  Z', 'Total']):
                    valstr = '{:<5s} :'.format(a)
                    valstr += '{:17.6f} a.u.'.format(dip_au[i])
                    valstr += '{:17.6f} Debye   '.format(dip_debye[i])
                    self.ostream.print_header(valstr)

        self.ostream.print_blank()
        self.ostream.flush()
