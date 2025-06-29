#
#                                   VELOXCHEM
#              ----------------------------------------------------
#                          An Electronic Structure Code
#
#  SPDX-License-Identifier: BSD-3-Clause
#
#  Copyright 2018-2025 VeloxChem developers
#
#  Redistribution and use in source and binary forms, with or without modification,
#  are permitted provided that the following conditions are met:
#
#  1. Redistributions of source code must retain the above copyright notice, this
#     list of conditions and the following disclaimer.
#  2. Redistributions in binary form must reproduce the above copyright notice,
#     this list of conditions and the following disclaimer in the documentation
#     and/or other materials provided with the distribution.
#  3. Neither the name of the copyright holder nor the names of its contributors
#     may be used to endorse or promote products derived from this software without
#     specific prior written permission.
#
#  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
#  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
#  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
#  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
#  FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
#  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
#  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
#  HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
#  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
#  OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

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
