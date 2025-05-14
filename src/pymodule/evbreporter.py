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
            no_ext = '.'.join(forcegroup_file.split('.')[:-1])
            ext = forcegroup_file.split('.')[-1]
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
            include = ['energy']
            if self.report_velocities:
                include.append('velocities')
            if self.report_forces:
                include.append('forces')

            return {'steps': steps, 'periodic': True, 'include': include}

    def report(self, simulation, state):
        E = []
        for simulation_dict in self.simulation_dicts.values():
            e = self._get_potential_energy(simulation_dict['simulation'],
                                           simulation_dict['forcegroups'],
                                           state)
            E.append(e)
        E1_pes = E[0]
        E2_pes = E[1]
        E1_int = E[2]
        E2_int = E[3]

        Em_pes = self._get_potential_energy(simulation, 'pes')
        Em_int = self._get_potential_energy(simulation, 'int')

        line = f"{self.Lambda}"
        for e in E:
            line += f", {e:.10e}"
        line += f", {Em_pes:.10e}, {Em_int:.10e}\n"
        self.E_out.write(line)

        Em_fg = []
        E1_fg = []
        E2_fg = []
        if self.report_forcegroups:
            line = ""
            reasim = self.simulation_dicts['reactant_integration']['simulation']
            prosim = self.simulation_dicts['product_integration']['simulation']
            for fg in EvbForceGroup:
                em = self._get_potential_energy(simulation, fg)
                e1 = self._get_potential_energy(reasim, fg)
                e2 = self._get_potential_energy(prosim, fg)

                Em_fg.append(em)
                E1_fg.append(e1)
                E2_fg.append(e2)

                if em > 1e9:
                    raise ValueError(
                        f"Force group {fg.name}({fg.value}) energy is too large: {em}"
                    )

            Em_line = ""
            E1_line = ""
            E2_line = ""

            for em, e1, e2 in zip(Em_fg, E1_fg, E2_fg):
                Em_line += f"{em}, "
                E1_line += f"{e1}, "
                E2_line += f"{e2}, "
            Em_line = Em_line[:-2] + '\n'
            E1_line = E1_line[:-2] + '\n'
            E2_line = E2_line[:-2] + '\n'
            self.FG_out.write(Em_line)
            self.rea_FG_out.write(E1_line)
            self.pro_FG_out.write(E2_line)

        if self.report_forces:
            forces = state.getForces(asNumpy=True)
            norms = np.linalg.norm(forces, axis=1)
            line = f"{self.Lambda}"
            kjpermolenm = mm.unit.kilojoules_per_mole / mm.unit.nanometer
            for i in range(forces.shape[0]):
                line += f", {forces[i][0].value_in_unit(kjpermolenm):.5e}, {forces[i][1].value_in_unit(kjpermolenm):.5e}, {forces[i][2].value_in_unit(kjpermolenm):.5e}, {norms[i]:.5e}"
            line += '\n'
            self.F_out.write(line)

        if self.report_velocities:
            velocities = state.getVelocities(asNumpy=True)
            line = f"{self.Lambda}"
            nmperps = mm.unit.nanometer / mm.unit.picosecond
            for i in range(velocities.shape[0]):
                line += f", {velocities[i][0].value_in_unit(nmperps):.5e}, {velocities[i][1].value_in_unit(nmperps):.5e}, {velocities[i][2].value_in_unit(nmperps):.5e}"
            line += '\n'
            self.v_out.write(line)

        if self.debug:
            # em = (1 - l) * E1 + l * E2
            pes_recalc = (1 - self.Lambda) * E1_pes + self.Lambda * E2_pes
            int_recalc = (1 - self.Lambda) * E1_int + self.Lambda * E2_int
            if abs(pes_recalc - Em_pes) > 1e-1:
                self.ostream.print_info(
                    f"Em pes recalculation is not consistent: {pes_recalc:.3f}(recalc)!={Em_pes:.3f}(em)"
                )

            if abs(int_recalc - Em_int) > 1e-1:
                self.ostream.print_info(
                    f"Em int recalculation is not consistent: {int_recalc:.3f}(recalc)!={Em_int:.3f}(em)"
                )

            if self.report_forcegroups:
                Em_pes_fg = 0
                E1_pes_fg = 0
                E2_pes_fg = 0
                Em_int_fg = 0
                E1_int_fg = 0
                E2_int_fg = 0
                for em, e1, e2, fg in zip(Em_fg, E1_fg, E2_fg, EvbForceGroup):
                    em_rec = (1 - self.Lambda) * e1 + self.Lambda * e2
                    # em = (1 - l) * E1 + l * E2 per forcegroup
                    if abs(em_rec - em) > 1e-2:
                        self.ostream.print_info(
                            f"em recalculation for force group {fg.name}({fg.value}) is not consistent: {em_rec:.3f}(recalc) != {em:.3f}(em)"
                        )

                    if fg.value in EvbForceGroup.pes_force_groups():
                        Em_pes_fg += em
                        E1_pes_fg += e1
                        E2_pes_fg += e2
                    if fg.value in EvbForceGroup.integration_force_groups():
                        Em_int_fg += em
                        E1_int_fg += e1
                        E2_int_fg += e2

                # check if forcegroups add up to E1, E2 and Em
                labels = [
                    "Em_pes", "E1_pes", "E2_pes", "Em_int", "E1_int", "E2_int"
                ]
                fg_recalc = [
                    Em_pes_fg, E1_pes_fg, E2_pes_fg, Em_int_fg, E1_int_fg,
                    E2_int_fg
                ]
                original = [Em_pes, E1_pes, E2_pes, Em_int, E1_int, E2_int]
                for label, recalc, orig in zip(labels, fg_recalc, original):
                    if abs(recalc - orig) > 1e-1:
                        self.ostream.print_info(
                            f"Force group summing {label} energy is not consistent: {recalc:.3f}(recalc) != {orig:.3f}(original)"
                        )
            self.ostream.flush()

    @staticmethod
    def _get_potential_energy(simulation, forcegroups=None, state=None):
        """
        Get the potential energy of the system.
        """
        # return 0
        if state is not None:
            simulation.context.setState(state)

        if forcegroups is None:
            return simulation.context.getState(
                getEnergy=True, ).getPotentialEnergy().value_in_unit(
                    mm.unit.kilojoules_per_mole)
        else:
            if isinstance(forcegroups, str):
                if forcegroups == 'pes':
                    forcegroups = EvbForceGroup.pes_force_groups()
                elif forcegroups == 'int':
                    forcegroups = EvbForceGroup.integration_force_groups()
                else:
                    raise ValueError(
                        f"Unknown force group: {forcegroups}. Use 'pes' or 'int'."
                    )
            if isinstance(forcegroups, EvbForceGroup):
                forcegroups = [forcegroups.value]

            if not isinstance(forcegroups, set):
                forcegroups = set(forcegroups)

            return simulation.context.getState(
                getEnergy=True,
                groups=forcegroups,
            ).getPotentialEnergy().value_in_unit(mm.unit.kilojoules_per_mole)
