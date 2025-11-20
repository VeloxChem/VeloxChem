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
from .evbsystembuilder import EvbForceGroup, EvbSystemBuilder

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
        systems,
        topology,
        lambda_val,
        outputstream,
        forming_bonds=None,
        breaking_bonds=None,
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

        self.out_streams = []

        self.E_out = open(energy_file, 'a' if append else 'w')
        self.out_streams.append(self.E_out)
        self.report_interval = report_interval

        self.lambda_val = lambda_val

        self.simulations = {}
        for name, system in systems.items():
            sim = mmapp.Simulation(topology, system, mm.VerletIntegrator(1))
            self.simulations.update({name: sim})

        if not append:
            header = "Lambda, reactant PES, product PES, reactant integration, product integration, Em \n"
            self.E_out.write(header)

        if force_file is None:
            self.report_forces = False
        else:
            self.report_forces = True
            self.F_out = open(force_file, 'a' if append else 'w')
            self.out_streams.append(self.F_out)
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
            self.out_streams.append(self.v_out)
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
            self.out_streams.append(self.FG_out)
            self.out_streams.append(self.rea_FG_out)
            self.out_streams.append(self.pro_FG_out)
            if not append:
                fg_header = EvbForceGroup.get_header()
                self.FG_out.write(fg_header)
                self.rea_FG_out.write(fg_header)
                self.pro_FG_out.write(fg_header)

        self.decomp_names = [s for s in systems if 'decomp' in str(s)]
        self.report_nb_decomp = False
        if len(self.decomp_names)>0:
            self.report_nb_decomp = True
            dir = '/'.join(energy_file.split('/')[:-1])
            filename = dir + '/NB_decompositions.csv'
            self.decomp_out = open(filename,'a' if append else 'w')
            self.out_streams.append(self.decomp_out)
            if not append:
                header = ", ".join(self.decomp_names)
                header+='\n'
                self.decomp_out.write(header)
        
        self.report_bonded_decomp = False
        if 'reactant_bonded_decomp' in systems.keys():
            if forming_bonds is None or breaking_bonds is None:
                self.ostream.print_warning(
                    "Formed and broken bonds need to be supplied to do bonded decomposition"
                )
                self.ostream.flush()
                return
            self.report_bonded_decomp = True
            active_atoms = []
            for bond in forming_bonds + breaking_bonds:
                if bond[0] not in active_atoms:
                    active_atoms.append(bond[0])
                if bond[1] not in active_atoms:
                    active_atoms.append(bond[1])

            self.reactant_params = self._get_bonded_decomp_params(systems['reactant'],active_atoms)
            self.product_params = self._get_bonded_decomp_params(systems['product'],active_atoms)
            self.measure_params = set()
            for force in list(self.reactant_params.values()) + list(self.product_params.values()):
                for params in force.values():
                    self.measure_params.add(params[0])
            # self.measure_params = sorted(self.measure_params)
            self.measure_params = sorted(self.measure_params, key=lambda x: (len(x), x))
            dir = '/'.join(energy_file.split('/')[:-1])
            filename = dir + '/bonded_E1_decomp.csv'
            self.bonded_E1_decomp_out = open(filename, 'a' if append else 'w')
            filename = dir + '/bonded_E2_decomp.csv'
            self.bonded_E2_decomp_out = open(filename, 'a' if append else 'w')
            filename = dir + '/bonded_params.csv'
            self.bonded_params_out = open(filename, 'a' if append else 'w')
            self.out_streams.append(self.bonded_E1_decomp_out)
            self.out_streams.append(self.bonded_E2_decomp_out)
            self.out_streams.append(self.bonded_params_out)
            if not append:
                rea_header = ""
                for force in self.reactant_params.values():
                    for param in force.values():
                        rea_header += f"{param[0]}, "
                rea_header = rea_header[:-2] + '\n'
                self.bonded_E1_decomp_out.write(rea_header)
                pro_header = ""
                for force in self.product_params.values():
                    for param in force.values():
                        pro_header += f"{param[0]}, "
                pro_header = pro_header[:-2] + '\n'
                self.bonded_E2_decomp_out.write(pro_header)
                params_header = ""
                for param in self.measure_params:
                    params_header += str(param) + ", "
                params_header = params_header[:-2] + '\n'
                self.bonded_params_out.write(params_header)

        for stream in self.out_streams:
            stream.flush()

    @staticmethod
    def _get_bonded_decomp_params(system,active_atoms):
        params = {}
        for force in system.getForces():
            if force.getForceGroup() == EvbForceGroup.REA_MORSE_BOND.value:
                morse_params = {}
                for i in range(force.getNumBonds()):
                    p1,p2,(D,a,r) = force.getBondParameters(i)
                    if p1 in active_atoms or p2 in active_atoms:
                        morse_params[i] = ((p1, p2), (D, a, r))
                params.update({force.getName(): morse_params})

            if force.getForceGroup() == EvbForceGroup.REA_HARM_BOND_DYNAMIC.value or force.getForceGroup() == EvbForceGroup.REA_HARM_BOND_STATIC.value:
                harm_params = {}
                for i in range(force.getNumBonds()):
                    p1,p2,r,k = force.getBondParameters(i)
                    if p1 in active_atoms or p2 in active_atoms:
                        harm_params[i] = ((p1, p2), (r,k))
                params.update({force.getName(): harm_params})

            if force.getForceGroup() == EvbForceGroup.REA_ANGLE.value:
                angle_params = {}
                for i in range(force.getNumAngles()):
                    p1,p2,p3,theta,k = force.getAngleParameters(i)
                    if p1 in active_atoms or p2 in active_atoms or p3 in active_atoms:
                        angle_params[i] = ((p1, p2, p3), (theta, k))
                params.update({force.getName(): angle_params})
            
            if force.getForceGroup() == EvbForceGroup.REA_TORSION.value or force.getForceGroup() == EvbForceGroup.REA_IMP.value:
                torsion_params = {}
                for i in range(force.getNumTorsions()):
                    p1,p2,p3,p4,periodicity,phase,k = force.getTorsionParameters(i)
                    if p1 in active_atoms or p2 in active_atoms or p3 in active_atoms or p4 in active_atoms:
                        torsion_params[i] = ((p1, p2, p3, p4), (periodicity, phase, k))
                params.update({force.getName(): torsion_params})
        return params

    def __del__(self):
        for stream in self.out_streams:
            if not stream.closed:
                stream.close()

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
            if self.report_bonded_decomp:
                include.append('positions')

            return {'steps': steps, 'periodic': True, 'include': include}

    def report(self, simulation, state):
        E = {}
        for name, sim in self.simulations.items():
            if name == 'reactant_bonded_decomp' or name == 'product_bonded_decomp':
                # Skip bonded decomposition systems
                continue
            e = self._get_potential_energy(
                sim,
                state=state,
            )
            E.update({name:e})
        E1_pes = E['reactant']
        E2_pes = E['product']
        E1_int = E[0]
        E2_int = E[1]


        Em = E1_pes * (1 - self.lambda_val) + E2_pes * self.lambda_val
        line = f"{self.lambda_val}, {E1_pes:.10e}, {E2_pes:.10e}, {E1_int:.10e}, {E2_int:.10e}, {Em:.10e} \n"
        self.E_out.write(line)

        # Em_dif = abs(Em_int - E1_int * (1 - self.lambda_val) +
        #              E2_int * self.lambda_val)

        # if Em_dif > 1e-2:
        #     self.ostream.print_info(str(Em_dif))
        #     self.ostream.flush()
        # assert  <1e-2

        Em_fg = []
        E1_fg = []
        E2_fg = []
        if self.report_forcegroups:
            line = ""
            reasim = self.simulations['reactant']
            prosim = self.simulations['product']
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

        if self.report_bonded_decomp:
            reasim = self.simulations['reactant_bonded_decomp']
            E1 = self._get_bonded_decomp_energy(reasim, state,self.reactant_params)
            line = ", ".join([f"{e:.10e}" for e in E1]) + '\n'
            self.bonded_E1_decomp_out.write(line)
            pro_sim = self.simulations['product_bonded_decomp']
            E2 = self._get_bonded_decomp_energy(pro_sim, state,self.product_params)
            line = ", ".join([f"{e:.10e}" for e in E2]) + '\n'
            self.bonded_E2_decomp_out.write(line)
            
            positions = state.getPositions(asNumpy=True)
            line = ""
            for i, param in enumerate(self.measure_params):
                val = ""
                if len(param) == 2:
                    val = EvbSystemBuilder.measure_length(positions[param[0]], positions[param[1]],)
                elif len(param) == 3:
                    val = EvbSystemBuilder.measure_angle(
                        positions[param[0]], positions[param[1]], positions[param[2]],)
                elif len(param) == 4:
                    val = EvbSystemBuilder.measure_dihedral(
                        positions[param[0]], positions[param[1]], positions[param[2]], positions[param[3]],)
                line += f"{val:.10e}, "
            line = line[:-2] + '\n'
            self.bonded_params_out.write(line)
                

        if self.report_forces:
            forces = state.getForces(asNumpy=True)
            norms = np.linalg.norm(forces, axis=1)
            line = f"{self.lambda_val}"
            kjpermolenm = mm.unit.kilojoules_per_mole / mm.unit.nanometer
            for i in range(forces.shape[0]):
                line += f", {forces[i][0].value_in_unit(kjpermolenm):.5e}, {forces[i][1].value_in_unit(kjpermolenm):.5e}, {forces[i][2].value_in_unit(kjpermolenm):.5e}, {norms[i]:.5e}"
            line += '\n'
            self.F_out.write(line)

        if self.report_velocities:
            velocities = state.getVelocities(asNumpy=True)
            line = f"{self.lambda_val}"
            nmperps = mm.unit.nanometer / mm.unit.picosecond
            for i in range(velocities.shape[0]):
                line += f", {velocities[i][0].value_in_unit(nmperps):.5e}, {velocities[i][1].value_in_unit(nmperps):.5e}, {velocities[i][2].value_in_unit(nmperps):.5e}"
            line += '\n'
            self.v_out.write(line)
        
        if self.report_nb_decomp:
            line = ""
            for name in self.decomp_names:
                line += f"{E[name]:.10e}, "
            line= line[:-2]
            line+= '\n'
            self.decomp_out.write(line)
        
        for stream in self.out_streams:
            stream.flush()

    @staticmethod
    def _get_bonded_decomp_energy(sim,state,parameters):
        E = []
        for force, params in zip(sim.system.getForces(),parameters.values()):
            for i, param in params.items():
                if force.getForceGroup() == EvbForceGroup.REA_MORSE_BOND.value:
                    force.setBondParameters(
                        i,
                        param[0][0],
                        param[0][1],
                        param[1],
                    )
                    force.updateParametersInContext(sim.context)
                    e = EvbReporter._get_potential_energy(
                        sim, state=state,)
                    force.setBondParameters(
                        i,
                        param[0][0],
                        param[0][1],
                        (0,0,param[1][2]),
                    )
                elif force.getForceGroup() == EvbForceGroup.REA_HARM_BOND_DYNAMIC.value or force.getForceGroup() == EvbForceGroup.REA_HARM_BOND_STATIC.value:
                    force.setBondParameters(
                        i,
                        param[0][0],
                        param[0][1],
                        param[1][0],
                        param[1][1],
                    )
                    force.updateParametersInContext(sim.context)
                    e = EvbReporter._get_potential_energy(
                        sim, state=state,)
                    force.setBondParameters(
                        i,
                        param[0][0],
                        param[0][1],
                        param[1][0], 0,
                    )
                elif force.getForceGroup() == EvbForceGroup.REA_ANGLE.value:
                    force.setAngleParameters(
                        i,
                        param[0][0],
                        param[0][1],
                        param[0][2],
                        param[1][0],
                        param[1][1],
                    )
                    force.updateParametersInContext(sim.context)
                    e = EvbReporter._get_potential_energy(
                        sim, state=state,)
                    force.setAngleParameters(
                        i,
                        param[0][0],
                        param[0][1],
                        param[0][2],
                        param[1][0], 0,
                    )
                elif force.getForceGroup() == EvbForceGroup.REA_TORSION.value or force.getForceGroup() == EvbForceGroup.REA_IMP.value:
                    force.setTorsionParameters(
                        i,
                        param[0][0],
                        param[0][1],
                        param[0][2],
                        param[0][3],
                        param[1][0],
                        param[1][1],
                        param[1][2],
                    )
                    force.updateParametersInContext(sim.context)
                    e = EvbReporter._get_potential_energy(
                        sim, state=state,)
                    force.setTorsionParameters(
                        i,
                        param[0][0],
                        param[0][1],
                        param[0][2],
                        param[0][3],
                        param[1][0], 0, 0,
                    )
                else:
                    e = 0
                force.updateParametersInContext(sim.context)
                E.append(e)
        return E

    @staticmethod
    def _get_potential_energy(simulation, forcegroups=None, state=None):
        """
        Get the potential energy of the system.
        """
        # return 0
        if state is not None:
            try:
                simulation.context.setState(state)
            except:
                # Decomposition systems which have the barostat removed will throw an error on the above case
                simulation.context.setPositions(state.getPositions())

        if forcegroups is None:
            return simulation.context.getState(
                getEnergy=True, ).getPotentialEnergy().value_in_unit(
                    mm.unit.kilojoules_per_mole)
        else:
            if isinstance(forcegroups, EvbForceGroup):
                forcegroups = set([forcegroups.value] )

            return simulation.context.getState(
                getEnergy=True,
                groups=forcegroups,
            ).getPotentialEnergy().value_in_unit(mm.unit.kilojoules_per_mole)
