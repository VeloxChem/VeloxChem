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

from importlib.metadata import version

import sys
import numpy as np
import copy

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
        debug=False,
    ):
        self.ostream = outputstream
        self.debug = debug
        # OpenMM HIP version is slighly older and uses a different format for reporters
        if version('openmm') < '8.2':
            if not append:
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
            no_ext = forcegroup_file.split('.')[0]
            ext = forcegroup_file.split('.')[1]
            rea_fg = no_ext + '_rea.' + ext
            pro_fg = no_ext + '_pro.' + ext

            self.FG_out = open(forcegroup_file, 'a' if append else 'w')
            self.rea_FG_out = open(rea_fg, 'a' if append else 'w')
            self.pro_FG_out = open(pro_fg, 'a' if append else 'w')

            if not append:
                self.FG_out.write(EvbForceGroup.get_header())
                self.rea_FG_out.write(EvbForceGroup.get_header())
                self.pro_FG_out.write(EvbForceGroup.get_header())



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
        line = f"{self.Lambda}"
        E = []
        for simulation_dict in self.simulation_dicts.values():
            simulation_dict['simulation'].context.setState(state)
            e = simulation_dict['simulation'].context.getState(
                getEnergy=True,
                groups=simulation_dict['forcegroups'],
            ).getPotentialEnergy().value_in_unit(mm.unit.kilojoules_per_mole)
            E.append(e)

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
        self.E_out.write(line)

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

        fg_E = []
        if self.report_forcegroups:
            line = ""
            for fg in EvbForceGroup:
                e = simulation.context.getState(
                    getEnergy=True,
                    groups=set([fg.value]),
                ).getPotentialEnergy().value_in_unit(
                    mm.unit.kilojoules_per_mole)
                fg_E.append(e)
                if e> 1e9:
                    raise ValueError(
                        f"Force group {fg.name}({fg.value}) energy is too large: {e}")
                line += f"{e}, "
            line = line[:-2] + '\n'
            self.FG_out.write(line)

            line = ""
            reasim = self.simulation_dicts['reactant_integration']['simulation']
            for fg in EvbForceGroup:
                e = reasim.context.getState(
                    getEnergy=True,
                    groups=set([fg.value]),
                ).getPotentialEnergy().value_in_unit(
                    mm.unit.kilojoules_per_mole)
                fg_E.append(e)
                if e> 1e9:
                    raise ValueError(
                        f"Force group {fg.name}({fg.value}) energy is too large: {e}")
                line += f"{e}, "
            line = line[:-2] + '\n'
            self.rea_FG_out.write(line)

            line = ""
            prosim = self.simulation_dicts['reactant_integration']['simulation']
            for fg in EvbForceGroup:
                e = prosim.context.getState(
                    getEnergy=True,
                    groups=set([fg.value]),
                ).getPotentialEnergy().value_in_unit(
                    mm.unit.kilojoules_per_mole)
                fg_E.append(e)
                if e> 1e9:
                    raise ValueError(
                        f"Force group {fg.name}({fg.value}) energy is too large: {e}")
                line += f"{e}, "
            line = line[:-2] + '\n'
            self.pro_FG_out.write(line)


        if self.debug:

            # verification
            pes_recalc = (1 - self.Lambda) * E[0] + self.Lambda * E[1]

            if abs(pes_recalc - Em_pes) > 1e1:
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
                self.ostream.print_info(
                    f"Lambda: {self.Lambda}, pes difference: {Em_pes - pes_recalc}, NB force difference: {dif[3]}"
                )

            if self.report_forcegroups:
                pes_fg_recalc = 0
                int_fg_recalc = 0
                for i, e in enumerate(fg_E):
                    if i + 1 in EvbForceGroup.pes_force_groups():
                        pes_fg_recalc += e
                    if i + 1 in EvbForceGroup.integration_force_groups():
                        int_fg_recalc += e
                    
                if abs(pes_fg_recalc - Em_pes) > 1e-1:
                    self.ostream.print_info(
                        f"Force group pes difference: {Em_pes - pes_fg_recalc}")
                if abs(int_fg_recalc - Em_int) > 1e-1:
                    self.ostream.print_info(
                        f"Force group int difference: {Em_int - int_fg_recalc}")
                self.ostream.flush()
