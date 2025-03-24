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
        forcegroup_file=None,
        force_file=None,
        velocity_file=None,
        append=False,
    ):

        # OpenMM HIP version is slighly older and uses a different format for reporters
        if version('openmm') < '8.2':
            outputstream.print_info(
                'Older version of OpenMM detected. Using tuple format for returning reporter information.'
            )
            outputstream.flush()
            self.use_tuple = True
        else:
            self.use_tuple = False

        self.E_out = open(energy_file, 'a' if append else 'w')
        self.report_interval = report_interval

        self.Lambda = Lambda
        self.simulation_dicts = {}

        self.simulation_dicts.update({
            'reactant_pes': {
                "simulation":
                mmapp.Simulation(topology, reactant_ff, mm.VerletIntegrator(1)),
                "forcegroups":
                EvbForceGroup.pes_force_groups(),
            },
            'product_pes': {
                "simulation":
                mmapp.Simulation(topology, product_ff, mm.VerletIntegrator(1)),
                "forcegroups":
                EvbForceGroup.pes_force_groups(),
            },
            'reactant_integration': {
                "simulation":
                mmapp.Simulation(topology, reactant_ff, mm.VerletIntegrator(1)),
                "forcegroups":
                EvbForceGroup.integration_force_groups(),
            },
            'product_integration': {
                "simulation":
                mmapp.Simulation(topology, product_ff, mm.VerletIntegrator(1)),
                "forcegroups":
                EvbForceGroup.integration_force_groups(),
            },
        })

        if not append:
            header = "Lambda, reactant PES, product PES, reactant integration, product integration, E_m_pes, E_m_int\n"
            self.E_out.write(header)

        if force_file is None:
            self.report_forces = False
        else:
            self.report_forces = True
            self.F_out = open(force_file, 'a' if append else 'w')
            if not append:
                header = "Lambda, "
                for j in range(topology.getNumAtoms()):
                    header += f"F(x, {j}), F(y, {j}), F(z, {j}), norm({j}), "
                header = header[:-2] + '\n'
                self.F_out.write(header)

        if velocity_file is None:
            self.report_velocities = False
        else:
            self.report_velocities = True
            self.v_out = open(velocity_file, 'a' if append else 'w')
            if not append:
                header = "Lambda, "
                for j in range(topology.getNumAtoms()):
                    header += f"V(x, {j}), V(y, {j}), V(z, {j}), "
                header = header[:-2] + '\n'
                self.v_out.write(header)

        if forcegroup_file is None:
            self.report_forcegroups = False
        else:
            self.report_forcegroups = True
            self.FG_out = open(forcegroup_file, 'a' if append else 'w')
            if not append:
                header = ""
                for fg in EvbForceGroup:
                    header += f"{fg.name}, "
                header = header[:-2] + '\n'
                self.FG_out.write(header)

    def __del__(self):
        self.E_out.close()

    def describeNextReport(self, simulation):
        steps = self.report_interval - simulation.currentStep % self.report_interval

        if self.use_tuple:
            return (steps, True, self.report_velocities, self.report_forces,
                    True, True
                    )  #steps, positions, velocities, forces, energy, pbc
        else:
            include = ['positions', 'energy']
            if self.report_velocities:
                include.append('velocities')
            if self.report_forces:
                include.append('forces')

            return {'steps': steps, 'periodic': True, 'include': include}

    def report(self, simulation, state):

        positions = state.getPositions(asNumpy=True)
        boxvectors = state.getPeriodicBoxVectors(asNumpy=True)
        # apparently this is necessary
        # simulation.context.setPositions(positions)
        # simulation.context.setPeriodicBoxVectors(*boxvectors)

        line = f"{self.Lambda}"
        E = []
        for simulation_dict in self.simulation_dicts.values():
            simulation_dict['simulation'].context.setPositions(positions)
            simulation_dict['simulation'].context.setPeriodicBoxVectors(
                *boxvectors)
            _state = simulation_dict['simulation'].context.getState(
                getEnergy=True,
                groups=simulation_dict['forcegroups'],
            )
            E.append(_state.getPotentialEnergy().value_in_unit(
                mm.unit.kilojoules_per_mole))

        Em_pes = simulation.context.getState(
            getEnergy=True,
            groups=EvbForceGroup.pes_force_groups(),
        ).getPotentialEnergy().value_in_unit(mm.unit.kilojoules_per_mole)

        Em_int = simulation.context.getState(
            getEnergy=True,
            groups=EvbForceGroup.integration_force_groups(),
        ).getPotentialEnergy().value_in_unit(mm.unit.kilojoules_per_mole)

        for e in E:
            line += f", {e}"
        line += f", {Em_pes}, {Em_int}\n"

        # verification
        #todo so this keeps failing for condensed systems, why
        pes_recalc = (1 - self.Lambda) * E[0] + self.Lambda * E[1]
        # pes_recalc = E[0]
        self.E_out.write(line)

        if abs(pes_recalc - Em_pes) > 1e1:
            # just go by every forcegroup, only shas to be done in pes

            simulations = [
                simulation,
                self.simulation_dicts['reactant_pes']['simulation'],
                self.simulation_dicts['product_pes']['simulation'],
            ]
            forcegroups = EvbForceGroup.pes_force_groups()
            E_contributions = np.zeros((len(simulations), len(forcegroups)))
            for si, simulation in enumerate(simulations):
                for fg, forcegroup in enumerate(forcegroups):
                    e = simulation.context.getState(
                        getEnergy=True,
                        groups=set([forcegroup]),
                    ).getPotentialEnergy().value_in_unit(
                        mm.unit.kilojoules_per_mole)
                    E_contributions[si, fg] = e
            dif = E_contributions[0] - (1 - self.Lambda) * E_contributions[
                1] + self.Lambda * E_contributions[2]
            dif = E_contributions[0] - E_contributions[1]
            print(
                f"Lambda: {self.Lambda}, pes difference: {Em_pes - pes_recalc}, NB force difference: {dif[3]}"
            )

            context0 = simulations[0].context
            sys0 = context0.getSystem()
            state0 = context0.getState(
                positions=True,
                energy=True,
                groups=EvbForceGroup.pes_force_groups(),
            )
            pos0 = state0.getPositions(asNumpy=True).value_in_unit(
                mm.unit.nanometer)
            E0 = state0.getPotentialEnergy().value_in_unit(
                mm.unit.kilojoules_per_mole)

            context1 = simulations[1].context
            sys1 = context1.getSystem()

            state1 = context1.getState(
                positions=True,
                energy=True,
                groups=EvbForceGroup.pes_force_groups(),
            )

            pos1 = state1.getPositions(asNumpy=True).value_in_unit(
                mm.unit.nanometer)
            E1 = state1.getPotentialEnergy().value_in_unit(
                mm.unit.kilojoules_per_mole)

            sys_string0 = mm.XmlSerializer.serialize(sys0)
            sys_string1 = mm.XmlSerializer.serialize(sys1)
            state_string0 = mm.XmlSerializer.serialize(state0)
            state_string1 = mm.XmlSerializer.serialize(state1)
            print(
                f"Systems equal: {sys_string0==sys_string1}, States equal: {state_string0==state_string1}, Positions equal: {np.all(pos0== pos1)}, Energies difference: {E0 - E1}"
            )
            state_string0 = state_string0.splitlines()
            state_string1 = state_string1.splitlines()
            for i, (line0, line1) in enumerate(zip(state_string0,
                                                   state_string1)):
                if line0 != line1:
                    print(f"line {i}: {line0} vs {line1}")

            # what line are sys0 and sys1 different

            pass

        # if abs(int_recalc - Em_int) > 1e-6:
        #     print(f"int difference: {Em_int - int_recalc}")

        if self.report_forces:
            forces = state.getForces(asNumpy=True)
            norms = np.linalg.norm(forces, axis=1)
            line = f"{self.Lambda}"
            kjpermolenm = mm.unit.kilojoules_per_mole / mm.unit.nanometer
            for i in range(forces.shape[0]):
                line += f", {forces[i][0].value_in_unit(kjpermolenm)}, {forces[i][1].value_in_unit(kjpermolenm)}, {forces[i][2].value_in_unit(kjpermolenm)}, {norms[i]}"
            line += '\n'
            self.F_out.write(line)

        if self.report_velocities:
            velocities = state.getVelocities(asNumpy=True)
            line = f"{self.Lambda}"
            nmperps = mm.unit.nanometer / mm.unit.picosecond
            for i in range(velocities.shape[0]):
                line += f", {velocities[i][0].value_in_unit(nmperps)}, {velocities[i][1].value_in_unit(nmperps)}, {velocities[i][2].value_in_unit(nmperps)}"
            line += '\n'
            self.v_out.write(line)

        if self.report_forcegroups:
            line = ""
            for fg in EvbForceGroup:
                e = simulation.context.getState(
                    getEnergy=True,
                    groups=set([fg.value]),
                ).getPotentialEnergy().value_in_unit(
                    mm.unit.kilojoules_per_mole)
                line += f"{e}, "
            line = line[:-2] + '\n'
            self.FG_out.write(line)
