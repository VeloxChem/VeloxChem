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
import numpy as np

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
        energy_file,
        report_interval,
        reactant_ff,
        product_ff,
        topology,
        Lambda,
        outputstream,
        force_file=None,
        velocity_file=None,
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

        self.E_out = open(energy_file, 'a' if append else 'w')
        self.report_interval = report_interval

        self.Lambda = Lambda
        self.simulation_dicts = {}

        self.simulation_dicts.update({
            'reactant_pes':
            self._get_simulation_dict(
                topology,
                reactant_ff,
                EvbForceGroup.pes_force_groups(),
            ),
            'product_pes':
            self._get_simulation_dict(
                topology,
                product_ff,
                EvbForceGroup.pes_force_groups(),
            ),
            'reactant_integrator':
            self._get_simulation_dict(
                topology,
                reactant_ff,
                EvbForceGroup.integration_force_groups(),
            ),
            'product_integrator':
            self._get_simulation_dict(
                topology,
                product_ff,
                EvbForceGroup.integration_force_groups(),
            ),
            'constraints':
            self._get_simulation_dict(
                topology,
                reactant_ff,
                EvbForceGroup.CONSTRAINT.value,
            ),
        })

        if not append:
            header = "Lambda, reactant PES, product PES, reactant integration, product integration, E_m, Constraints\n"
            self.E_out.write(header)

        if force_file is None:
            self.forces = False
        else:
            self.forces = True
            self.F_out = open(force_file, 'a' if append else 'w')
            if not append:
                header = "Lambda, "
                for j in range(topology.getNumAtoms()):
                    header += f"F(x, {j}), F(y, {j}), F(z, {j}), norm({j}), "
                header = header[:-2] + '\n'
                self.F_out.write(header)

        if velocity_file is None:
            self.velocities = False
        else:
            self.velocities = True
            self.v_out = open(velocity_file, 'a' if append else 'w')
            if not append:
                header = "Lambda, "
                for j in range(topology.getNumAtoms()):
                    header += f"V(x, {j}), V(y, {j}), V(z, {j}), "
                header = header[:-2] + '\n'
                self.v_out.write(header)

    def __del__(self):
        self.E_out.close()

    @staticmethod
    def _get_simulation_dict(topology, ff, forcegroups):
        integrator = mm.VerletIntegrator(1)
        integrator.setIntegrationForceGroups(forcegroups)
        simulation = mmapp.Simulation(topology, ff, integrator)
        return {"simulation": simulation, "forcegroups": forcegroups}

    def describeNextReport(self, simulation):
        steps = self.report_interval - simulation.currentStep % self.report_interval

        if self.use_tuple:
            return (steps, True, self.velocities, self.forces, True, True
                    )  #steps, positions, velocities, forces, energy, pbc
        else:
            include = ['positions', 'energy']
            if self.velocities:
                include.append('velocities')
            if self.forces:
                include.append('forces')

            return {'steps': steps, 'periodic': True, 'include': include}

    def report(self, simulation, state):

        positions = state.getPositions(asNumpy=True)
        line = f"{self.Lambda}"
        for simulation_dict in self.simulation_dicts.values():
            E = self._get_energy(simulation_dict, positions)
            line += f", {E}"

        Em = state.getPotentialEnergy().value_in_unit(mm.unit.kilojoules_per_mole)
        line += f", {Em}\n"
        self.E_out.write(line)

        if self.forces:
            forces = state.getForces(asNumpy=True)
            norms = np.linalg.norm(forces, axis=1)
            line = f"{self.Lambda}"
            kjpermolenm = mm.unit.kilojoules_per_mole / mm.unit.nanometer
            for i in range(forces.shape[0]):
                line += f", {forces[i][0].value_in_unit(kjpermolenm)}, {forces[i][1].value_in_unit(kjpermolenm)}, {forces[i][2].value_in_unit(kjpermolenm)}, {norms[i]}"
            line += '\n'
            self.F_out.write(line)
        if self.velocities:
            velocities = state.getVelocities(asNumpy=True)
            line = f"{self.Lambda}"
            nmperps = mm.unit.nanometer / mm.unit.picosecond
            for i in range(velocities.shape[0]):
                line += f", {velocities[i][0].value_in_unit(nmperps)}, {velocities[i][1].value_in_unit(nmperps)}, {velocities[i][2].value_in_unit(nmperps)}"
            line += '\n'
            self.v_out.write(line)

    def _get_energy(self, simulation_dict, positions):
        simulation = simulation_dict['simulation']
        forcegroups = simulation_dict['forcegroups']
        simulation.context.setPositions(positions)
        state = simulation.context.getState(
            getEnergy=True,
            groups=forcegroups,
        )
        return state.getPotentialEnergy().value_in_unit(mm.unit.kilojoules_per_mole)
