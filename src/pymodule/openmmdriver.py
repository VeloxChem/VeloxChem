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
import sys

from .veloxchemlib import mpi_master, bohr_in_angstrom, hartree_in_kjpermol
from .outputstream import OutputStream


class OpenMMDriver:
    """
    Implements OpenMM driver.

    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.

    Instance variables
        - integrator_temperature: The temperature.
        - integrator_friction_coeff: The friction coefficient.
        - integrator_step_size: The step size.
        - energy: The potential energy.
        - gradient: The gradient.
    """

    def __init__(self, comm=None, ostream=None):
        """
        Initializes OpenMM driver.
        """

        from openmm.unit import picosecond, picoseconds, kelvin

        if comm is None:
            comm = MPI.COMM_WORLD

        if ostream is None:
            if comm.Get_rank() == mpi_master():
                ostream = OutputStream(sys.stdout)
            else:
                ostream = OutputStream(None)

        self.comm = comm
        self.rank = comm.Get_rank()
        self.nodes = comm.Get_size()

        self.ostream = ostream

        self.integrator_temperature = 298.15 * kelvin
        self.integrator_friction_coeff = 1 / picosecond
        self.integrator_step_size = 0.004 * picoseconds

        self.energy = None
        self.gradient = None

    def add_topology(self, top_file, include_dir=None):
        """
        Adds topology file.

        :param top_file:
            The name of the topology file.
        :param include_dir:
            The include directory for Gromacs topology files.
        """

        from openmm.app import GromacsTopFile, NoCutoff, Simulation
        from openmm import LangevinIntegrator

        if self.rank == mpi_master():
            top = GromacsTopFile(top_file, includeDir=include_dir)

            system = top.createSystem(NoCutoff)
            integrator = LangevinIntegrator(
                self.integrator_temperature,
                self.integrator_friction_coeff,
                self.integrator_step_size,
            )

            self.simulation = Simulation(top.topology, system, integrator)

    def compute(self, molecule):
        """
        Computes MM potential energy and gradient.

        :param molecule:
            The molecule.
        """

        from openmm.unit import nanometer, md_unit_system

        if self.rank == mpi_master():
            coords_nm = molecule.get_coordinates_in_angstrom() * 0.1

            self.simulation.context.setPositions(coords_nm * nanometer)
            state = self.simulation.context.getState(
                getPositions=True,
                getEnergy=True,
                getForces=True,
            )

            self.energy = state.getPotentialEnergy().value_in_unit_system(
                md_unit_system)
            self.gradient = -1.0 * state.getForces(
                asNumpy=True).value_in_unit_system(md_unit_system)

            # convert to a.u.
            self.energy /= hartree_in_kjpermol()

            # convert to a.u.
            self.gradient /= (hartree_in_kjpermol() * 10.0 / bohr_in_angstrom())

        self.energy = self.comm.bcast(self.energy, root=mpi_master())
        self.gradient = self.comm.bcast(self.gradient, root=mpi_master())

    def get_energy(self):
        """
        Gets MM potential energy.

        :return:
            The potential energy.
        """

        return self.energy

    def get_gradient(self):
        """
        Gets MM gradient.

        :return:
            The gradient.
        """

        return self.gradient
