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

from importlib.metadata import version

import sys

from .errorhandler import assert_msg_critical
from .evbsystembuilder import EvbForceGroup

try:
    import openmm.app as mmapp
    import openmm as mm
except ImportError:
    pass


class EvbReporter():

    def __init__(
        self,
        file,
        report_interval,
        reactant_ff,
        product_ff,
        topology,
        Lambda,
        outputstream,
        force_file=None,
        append=False,
    ):

        # OpenMM HIP version is slighly older and uses a different format for reporters
        if version('openmm') < '8.2':
            outputstream.print_info(
                'Older version of OpenMM detected. Using tuple format for returning reporter information.')
            outputstream.flush()
            self.use_tuple = True
        else:
            self.use_tuple = False

        self.E_out = open(file, 'a' if append else 'w')
        self.report_interval = report_interval

        self.Lambda = Lambda
        self.simulations = {}

        self.simulations.update({
            'reactant_pes':
            self._get_simulation(
                topology,
                reactant_ff,
                EvbForceGroup.pes_force_groups(),
            ),
            'product_pes':
            self._get_simulation(
                topology,
                product_ff,
                EvbForceGroup.pes_force_groups(),
            ),
            'reactant_integrator':
            self._get_simulation(
                topology,
                reactant_ff,
                EvbForceGroup.integration_force_groups(),
            ),
            'product_integrator':
            self._get_simulation(
                topology,
                product_ff,
                EvbForceGroup.integration_force_groups(),
            ),
            'constraints':
            self._get_simulation(
                topology,
                reactant_ff,
                EvbForceGroup.CONSTRAINT.value,
            ),
        })

        if force_file is not None:
            self.forces = True
            self.F_out = open(force_file, 'a' if append else 'w')
        else:
            self.forces = False

        if not append:
            header = "Lambda, reactant PES, product PES, reactant integration, product integration, E_m, constraints\n"
            self.E_out.write(header)

    def __del__(self):
        self.E_out.close()

    @staticmethod
    def _get_simulation(topology, ff, forcegroups):
        integrator = mm.VerletIntegrator(1)
        integrator.setIntegrationForceGroups(forcegroups)
        return mmapp.Simulation(topology, ff, integrator)

    def describeNextReport(self, simulation):
        steps = self.report_interval - simulation.currentStep % self.report_interval

        if self.use_tuple:
            return (steps, True, False, self.forces, True, True)  #steps, positions, velocities, forces, energy, pbc
        else:
            if self.forces:
                include = ['forces', 'positions', 'energy']
            else:
                include = ['positions', 'energy']

            return {'steps': steps, 'periodic': True, 'include': include}

    def report(self, simulation, sim_state):

        positions = sim_state.getPositions(asNumpy=True)
        Em = sim_state.getPotentialEnergy().value_in_unit(mm.unit.kilojoules_per_mole)

        E_pes_reactant = self._get_energy(
            self.simulations['reactant_pes'],
            positions,
            EvbForceGroup.pes_force_groups(),
        )
        E_pes_product = self._get_energy(
            self.simulations['product_pes'],
            positions,
            EvbForceGroup.pes_force_groups(),
        )
        E_int_reactant = self._get_energy(
            self.simulations['reactant_integrator'],
            positions,
            EvbForceGroup.integration_force_groups(),
        )
        E_int_product = self._get_energy(
            self.simulations['product_integrator'],
            positions,
            EvbForceGroup.integration_force_groups(),
        )

        line = f"{self.Lambda}, {E_pes_reactant}, {E_pes_product}, {E_int_reactant}, {E_int_product}, {Em}"
        line += '\n'

        # if self.forces:
        #     forces = state.getForces(asNumpy=True)
        #     for i in range(forces.shape[0]):
        #         line += f", {forces[i][0]}, {forces[i][1]}, {forces[i][2]}"
        self.E_out.write(line)

    def _get_energy(self, simulation, positions, forcegroups):
        simulation.context.setPositions(positions)
        state = simulation.context.getState(
            getEnergy=True,
            groups=forcegroups,
        )
        return state.getPotentialEnergy().value_in_unit(mm.unit.kilojoules_per_mole)
