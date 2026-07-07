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
import re
import time
import gc

import numpy as np
from pathlib import Path
from mpi4py import MPI

from .errorhandler import assert_msg_critical
from .reactionsystembuilder import EvbForceGroup, ReactionSystemBuilder

try:
    import openmm.app as mmapp
    import openmm as mm
except ImportError:
    pass


class EvbReporter:

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
        replica=0,
        direction=0,
        cpu_threads=None,
        enable_bonded_decomp=True,
        defer_open=False,
        core_outputs=True,
    ):
        # core_outputs: when False, the energy / force / velocity / NB-decomp
        #   files are neither opened nor written (only force-group / bonded
        #   decomposition are). Used by the master-side reporter in async mode
        #   so it does not collide with the worker's Energies.csv.
        # cpu_threads: number of OpenMM CPU-platform threads for the energy
        #   simulations (None = OpenMM default). Used by the async reporter
        #   worker to saturate the CPU on a single node.
        # enable_bonded_decomp: allow the (expensive, GPU-system-dependent)
        #   bonded decomposition. The async worker disables it; the master-side
        #   synchronous reporter keeps it.
        # defer_open: if True, do not open the output files in __init__; the
        #   caller (worker) opens them per window via _open_outputs(append).
        self.ostream = outputstream
        self.debug = debug
        # Replica index and sweep direction (0 = forward l: 0->1, 1 = backward
        # l: 1->0) are stamped onto every reported row so that all replicas and
        # both directions can share a single output file.
        self.replica = replica
        self.direction = direction
        # OpenMM HIP version is slighly older and uses a different format for reporters
        raw_openmm_version = version('openmm')
        match = re.match(r"^(\d+)\.(\d+)", raw_openmm_version)
        if not match:
            raise RuntimeError("Cannot parse required major.minor version from "
                               f"OpenMM: {raw_openmm_version!r}")
        openmm_version = int(match.group(1)), int(match.group(2))
        if openmm_version < (8, 2):
            if not append:
                outputstream.print_info(
                    'Older version of OpenMM detected. Using tuple format for returning reporter information.'
                )
                outputstream.flush()
            self.use_tuple = True
        else:
            self.use_tuple = False

        self.out_streams = []
        self.report_interval = report_interval
        self.lambda_val = lambda_val
        self.topology = topology
        self.num_atoms = topology.getNumAtoms()

        # Store output paths / flags; actual file opening is deferred to
        # _open_outputs so the worker can (re)open per window with an append
        # flag while the CPU simulations below are built only once.
        self.core_outputs = core_outputs
        self._energy_file = energy_file
        self._force_file = force_file
        self._velocity_file = velocity_file
        self._forcegroup_file = forcegroup_file
        self.report_forces = force_file is not None and core_outputs
        self.report_velocities = velocity_file is not None and core_outputs
        self.report_forcegroups = forcegroup_file is not None

        # These auxiliary simulations only evaluate potential energies; run
        # them on the CPU platform so they don't each reserve a separate CUDA
        # context. Building them on the GPU can crash the GPU.
        cpu_platform = mm.Platform.getPlatformByName('CPU')
        if cpu_threads is not None:
            cpu_platform.setPropertyDefaultValue('Threads', str(cpu_threads))
        # Only the reactant/product PES and the two integration endpoints
        # (lambda 0 / 1) feed the energy output; the intermediate lambda-window
        # systems would be evaluated every frame and discarded, so skip building
        # them. Decomposition systems (nb / bonded) are kept because the
        # decomposition reporters below reference them by name.
        self.simulations = {}
        core_names = ['reactant', 'product', 0, 1]
        for name, system in systems.items():
            if name not in core_names and 'decomp' not in str(name):
                continue
            sim = mmapp.Simulation(topology,
                                   system,
                                   mm.LangevinIntegrator(1, 1, 1),
                                   platform=cpu_platform)
            self.simulations.update({name: sim})

        self.decomp_names = [s for s in systems if 'decomp' in str(s)]
        self.report_nb_decomp = len(self.decomp_names) > 0 and core_outputs

        # Bonded-decomposition parameters depend only on the systems, so they
        # are computed once here; the per-window file handling lives in
        # _open_outputs.
        self.report_bonded_decomp = False
        if enable_bonded_decomp and 'reactant_bonded_decomp' in systems.keys():
            if forming_bonds is None or breaking_bonds is None:
                self.ostream.print_warning(
                    "Formed and broken bonds need to be supplied to do bonded decomposition"
                )
                self.ostream.flush()
            else:
                self.report_bonded_decomp = True
                active_atoms = []
                for bond in forming_bonds + breaking_bonds:
                    if bond[0] not in active_atoms:
                        active_atoms.append(bond[0])
                    if bond[1] not in active_atoms:
                        active_atoms.append(bond[1])

                self.reactant_params = self._get_bonded_decomp_params(
                    systems['reactant'], active_atoms)
                self.product_params = self._get_bonded_decomp_params(
                    systems['product'], active_atoms)
                self.measure_params = set()
                for force in list(self.reactant_params.values()) + list(
                        self.product_params.values()):
                    for params in force.values():
                        self.measure_params.add(params[0])
                self.measure_params = sorted(self.measure_params,
                                             key=lambda x: (len(x), x))

        if not defer_open:
            self._open_outputs(append)

    @staticmethod
    def energies_header():
        """Header line for Energies.csv (shared by every recalculation mode)."""
        return ("Lambda, reactant PES, product PES, reactant integration, "
                "product integration, Em, replica, direction \n")

    @staticmethod
    def format_energies_row(lambda_val, E1_pes, E2_pes, E1_int, E2_int, replica,
                            direction):
        """Format one Energies.csv row. Em is the EVB-mixed diabatic energy."""
        Em = E1_pes * (1 - lambda_val) + E2_pes * lambda_val
        return (f"{lambda_val}, {E1_pes:.10e}, {E2_pes:.10e}, {E1_int:.10e}, "
                f"{E2_int:.10e}, {Em:.10e}, {replica}, {direction} \n")

    @staticmethod
    def forces_header(num_atoms):
        """Header line for Forces.csv (shared by every recalculation mode)."""
        header = "Lambda, replica, direction, "
        for j in range(num_atoms):
            header += f"F(x, {j}), F(y, {j}), F(z, {j}), norm({j}), "
        return header[:-2] + '\n'

    @staticmethod
    def format_forces_row(lambda_val, replica, direction, forces):
        """Format one Forces.csv row from an (N,3) force array (kJ/mol/nm)."""
        norms = np.linalg.norm(forces, axis=1)
        line = f"{lambda_val}, {replica}, {direction}"
        for i in range(forces.shape[0]):
            line += (f", {forces[i][0]:.5e}, {forces[i][1]:.5e}, "
                     f"{forces[i][2]:.5e}, {norms[i]:.5e}")
        return line + '\n'

    def _open_outputs(self, append):
        """Open (truncate + write headers, or append) all enabled output files.

        Called once by __init__ in synchronous mode, and once by the async
        worker on the first window (append=False truncates and writes headers;
        the worker keeps the streams open and appends for subsequent windows).
        """
        mode = 'a' if append else 'w'
        self.out_streams = []

        if self.core_outputs:
            self.E_out = open(self._energy_file, mode)
            self.out_streams.append(self.E_out)
            if not append:
                self.E_out.write(self.energies_header())

        if self.report_forces:
            self.F_out = open(self._force_file, mode)
            self.out_streams.append(self.F_out)
            if not append:
                self.F_out.write(self.forces_header(self.num_atoms))

        if self.report_velocities:
            self.v_out = open(self._velocity_file, mode)
            self.out_streams.append(self.v_out)
            if not append:
                header = "Lambda, replica, direction, "
                for j in range(self.num_atoms):
                    header += f"V(x, {j}), V(y, {j}), V(z, {j}), "
                self.v_out.write(header[:-2] + '\n')

        if self.report_forcegroups:
            forcegroup_path = Path(self._forcegroup_file)
            rea_fg = str(
                forcegroup_path.with_name(forcegroup_path.stem + '_rea' +
                                          forcegroup_path.suffix))
            pro_fg = str(
                forcegroup_path.with_name(forcegroup_path.stem + '_pro' +
                                          forcegroup_path.suffix))
            self.FG_out = open(self._forcegroup_file, mode)
            self.rea_FG_out = open(rea_fg, mode)
            self.pro_FG_out = open(pro_fg, mode)
            self.out_streams.append(self.FG_out)
            self.out_streams.append(self.rea_FG_out)
            self.out_streams.append(self.pro_FG_out)
            if not append:
                fg_header = EvbForceGroup.get_header()
                self.FG_out.write(fg_header)
                self.rea_FG_out.write(fg_header)
                self.pro_FG_out.write(fg_header)

        if self.report_nb_decomp:
            output_dir = Path(self._energy_file).parent
            filename = str(output_dir / 'NB_decompositions.csv')
            self.decomp_out = open(filename, mode)
            self.out_streams.append(self.decomp_out)
            if not append:
                self.decomp_out.write(", ".join(self.decomp_names) + '\n')

        if self.report_bonded_decomp:
            output_dir = Path(self._energy_file).parent
            filename = str(output_dir / 'bonded_E1_decomp.csv')
            self.bonded_E1_decomp_out = open(filename, mode)
            self.bonded_E2_decomp_out = open(filename, mode)
            self.bonded_params_out = open(filename, mode)
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
    def _get_bonded_decomp_params(system, active_atoms):
        params = {}
        for force in system.getForces():
            if force.getForceGroup() == EvbForceGroup.REA_MORSE_BOND.value:
                morse_params = {}
                for i in range(force.getNumBonds()):
                    p1, p2, (D, a, r) = force.getBondParameters(i)
                    if p1 in active_atoms or p2 in active_atoms:
                        morse_params[i] = ((p1, p2), (D, a, r))
                params.update({force.getName(): morse_params})

            if force.getForceGroup(
            ) == EvbForceGroup.REA_HARM_BOND_DYNAMIC.value or force.getForceGroup(
            ) == EvbForceGroup.REA_HARM_BOND_STATIC.value:
                harm_params = {}
                for i in range(force.getNumBonds()):
                    p1, p2, r, k = force.getBondParameters(i)
                    if p1 in active_atoms or p2 in active_atoms:
                        harm_params[i] = ((p1, p2), (r, k))
                params.update({force.getName(): harm_params})

            if force.getForceGroup() == EvbForceGroup.REA_ANGLE.value:
                angle_params = {}
                for i in range(force.getNumAngles()):
                    p1, p2, p3, theta, k = force.getAngleParameters(i)
                    if p1 in active_atoms or p2 in active_atoms or p3 in active_atoms:
                        angle_params[i] = ((p1, p2, p3), (theta, k))
                params.update({force.getName(): angle_params})

            if force.getForceGroup(
            ) == EvbForceGroup.REA_TORSION.value or force.getForceGroup(
            ) == EvbForceGroup.REA_IMP.value:
                torsion_params = {}
                for i in range(force.getNumTorsions()):
                    p1, p2, p3, p4, periodicity, phase, k = force.getTorsionParameters(
                        i)
                    if p1 in active_atoms or p2 in active_atoms or p3 in active_atoms or p4 in active_atoms:
                        torsion_params[i] = ((p1, p2, p3, p4), (periodicity,
                                                                phase, k))
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
                    )  # steps, positions, velocities, forces, energy, pbc
        else:
            include = ['energy', 'positions']
            if self.report_velocities:
                include.append('velocities')
            if self.report_forces:
                include.append('forces')

            return {'steps': steps, 'periodic': True, 'include': include}

    @staticmethod
    def _apply_positions(simulation, positions, box):
        """Set positions (nm) and periodic box (nm) on a sim from numpy.

        box may be None for a non-periodic system (e.g. vacuum), in which case
        the context keeps its default box (irrelevant to a NoCutoff energy).
        """
        nm = mm.unit.nanometer
        if box is not None:
            a, b, c = box
            simulation.context.setPeriodicBoxVectors(
                mm.Vec3(*a) * nm,
                mm.Vec3(*b) * nm,
                mm.Vec3(*c) * nm,
            )
        simulation.context.setPositions(positions * nm)

    def _write_core(self, positions, box, velocities=None, forces=None):
        """Write the per-frame core outputs (Energies / Forces / Velocities /
        NB-decomposition) from raw numpy arrays. This is the part that is
        offloaded to the async reporter worker; it needs no live GPU simulation.

        positions: (N,3) nm, box: (3,3) nm, velocities: (N,3) nm/ps,
        forces: (N,3) kJ/mol/nm.
        """
        E = {}
        for name, sim in self.simulations.items():
            self._apply_positions(sim, positions, box)
            E[name] = sim.context.getState(
                getEnergy=True).getPotentialEnergy().value_in_unit(
                    mm.unit.kilojoules_per_mole)
        # When core outputs are disabled (master-side force-group-only reporter)
        # we still apply positions above so the force-group block can read the
        # CPU sims, but skip writing the core CSV rows.
        if not self.core_outputs:
            return

        E1_pes = E['reactant']
        E2_pes = E['product']
        E1_int = E[0]
        E2_int = E[1]

        self.E_out.write(
            self.format_energies_row(self.lambda_val, E1_pes, E2_pes, E1_int,
                                     E2_int, self.replica, self.direction))

        if self.report_forces and forces is not None:
            self.F_out.write(
                self.format_forces_row(self.lambda_val, self.replica,
                                       self.direction, forces))

        if self.report_velocities and velocities is not None:
            line = f"{self.lambda_val}, {self.replica}, {self.direction}"
            for i in range(velocities.shape[0]):
                line += f", {velocities[i][0]:.5e}, {velocities[i][1]:.5e}, {velocities[i][2]:.5e}"
            self.v_out.write(line + '\n')

        if self.report_nb_decomp:
            line = ""
            for name in self.decomp_names:
                line += f"{E[name]:.10e}, "
            self.decomp_out.write(line[:-2] + '\n')

    def report(self, simulation, state):
        positions = state.getPositions(asNumpy=True).value_in_unit(
            mm.unit.nanometer)
        box = state.getPeriodicBoxVectors(asNumpy=True).value_in_unit(
            mm.unit.nanometer)
        velocities = None
        forces = None
        if self.report_velocities:
            velocities = state.getVelocities(asNumpy=True).value_in_unit(
                mm.unit.nanometer / mm.unit.picosecond)
        if self.report_forces:
            forces = state.getForces(asNumpy=True).value_in_unit(
                mm.unit.kilojoules_per_mole / mm.unit.nanometer)
        # Core outputs; this also applies the current positions to every CPU
        # simulation, which the force-group block below relies on.
        self._write_core(positions, box, velocities, forces)

        # Force-group and bonded-decomposition outputs depend on the live GPU
        # simulation / full state and stay synchronous.
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
            E1 = self._get_bonded_decomp_energy(reasim, state,
                                                self.reactant_params)
            line = ", ".join([f"{e:.10e}" for e in E1]) + '\n'
            self.bonded_E1_decomp_out.write(line)
            pro_sim = self.simulations['product_bonded_decomp']
            E2 = self._get_bonded_decomp_energy(pro_sim, state,
                                                self.product_params)
            line = ", ".join([f"{e:.10e}" for e in E2]) + '\n'
            self.bonded_E2_decomp_out.write(line)

            positions = state.getPositions(asNumpy=True)
            line = ""
            for i, param in enumerate(self.measure_params):
                val = ""
                if len(param) == 2:
                    val = ReactionSystemBuilder.measure_length(
                        positions[param[0]],
                        positions[param[1]],
                    )
                elif len(param) == 3:
                    val = ReactionSystemBuilder.measure_angle(
                        positions[param[0]],
                        positions[param[1]],
                        positions[param[2]],
                    )
                elif len(param) == 4:
                    val = ReactionSystemBuilder.measure_dihedral(
                        positions[param[0]],
                        positions[param[1]],
                        positions[param[2]],
                        positions[param[3]],
                    )
                line += f"{val:.10e}, "
            line = line[:-2] + '\n'
            self.bonded_params_out.write(line)

        for stream in self.out_streams:
            stream.flush()

    @staticmethod
    def _get_bonded_decomp_energy(sim, state, parameters):
        E = []
        for force, params in zip(sim.system.getForces(), parameters.values()):
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
                        sim,
                        state=state,
                    )
                    force.setBondParameters(
                        i,
                        param[0][0],
                        param[0][1],
                        (0, 0, param[1][2]),
                    )
                elif force.getForceGroup(
                ) == EvbForceGroup.REA_HARM_BOND_DYNAMIC.value or force.getForceGroup(
                ) == EvbForceGroup.REA_HARM_BOND_STATIC.value:
                    force.setBondParameters(
                        i,
                        param[0][0],
                        param[0][1],
                        param[1][0],
                        param[1][1],
                    )
                    force.updateParametersInContext(sim.context)
                    e = EvbReporter._get_potential_energy(
                        sim,
                        state=state,
                    )
                    force.setBondParameters(
                        i,
                        param[0][0],
                        param[0][1],
                        param[1][0],
                        0,
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
                        sim,
                        state=state,
                    )
                    force.setAngleParameters(
                        i,
                        param[0][0],
                        param[0][1],
                        param[0][2],
                        param[1][0],
                        0,
                    )
                elif force.getForceGroup(
                ) == EvbForceGroup.REA_TORSION.value or force.getForceGroup(
                ) == EvbForceGroup.REA_IMP.value:
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
                        sim,
                        state=state,
                    )
                    force.setTorsionParameters(
                        i,
                        param[0][0],
                        param[0][1],
                        param[0][2],
                        param[0][3],
                        param[1][0],
                        0,
                        0,
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
            except Exception:
                # Decomposition systems which have the barostat removed will throw an error on the above case
                simulation.context.setPositions(state.getPositions())

        if forcegroups is None:
            return simulation.context.getState(
                getEnergy=True, ).getPotentialEnergy().value_in_unit(
                    mm.unit.kilojoules_per_mole)
        else:
            if isinstance(forcegroups, EvbForceGroup):
                forcegroups = set([forcegroups.value])

            return simulation.context.getState(
                getEnergy=True,
                groups=forcegroups,
            ).getPotentialEnergy().value_in_unit(mm.unit.kilojoules_per_mole)


class EvbGpuRecalculator:
    """Batched GPU recalculation of EVB energies for deferred mode.

    Instead of re-scoring each sampled frame inline (blocking the GPU) or on a
    dedicated CPU worker rank, deferred mode samples a whole lambda window first
    and then hands the window's frames here. Because only one OpenMM context
    should live on the GPU at a time, each system is built, evaluated over every
    frame, and torn down before the next system is built - so the GPU holds a
    single context at any moment. The output files (Energies / optional Forces /
    NB decomposition) match those the synchronous ``EvbReporter`` writes.
    """

    def __init__(self,
                 data_folder,
                 systems,
                 topology,
                 simulation_factory,
                 report_forces=False):
        # simulation_factory(system) -> a live mmapp.Simulation on the desired
        # (GPU) platform. The factory owns platform selection; this class only
        # ever builds one context at a time and tears it down before the next.
        self.data_folder = Path(data_folder)
        self.systems = systems
        self.num_atoms = topology.getNumAtoms()
        self._make_sim = simulation_factory
        self.report_forces = report_forces

        # Same system selection as EvbReporter.simulations: the reactant/product
        # PES and the two integration endpoints (0/1), plus any nb/bonded
        # decomposition systems.
        core_names = ['reactant', 'product', 0, 1]
        self.core_names = [n for n in core_names if n in systems]
        self.decomp_names = [s for s in systems if 'decomp' in str(s)]

        self._opened = False
        self.E_out = None
        self.F_out = None
        self.decomp_out = None
        self.out_streams = []

    def _open_outputs(self, append):
        mode = 'a' if append else 'w'
        self.E_out = open(self.data_folder / "Energies.csv", mode)
        self.out_streams.append(self.E_out)
        if not append:
            self.E_out.write(EvbReporter.energies_header())
        if self.report_forces:
            self.F_out = open(self.data_folder / "Forces.csv", mode)
            self.out_streams.append(self.F_out)
            if not append:
                self.F_out.write(EvbReporter.forces_header(self.num_atoms))
        if self.decomp_names:
            self.decomp_out = open(self.data_folder / "NB_decompositions.csv",
                                   mode)
            self.out_streams.append(self.decomp_out)
            if not append:
                self.decomp_out.write(
                    ", ".join([str(n) for n in self.decomp_names]) + '\n')
        self._opened = True

    def _evaluate_system(self, system, frames, want_forces=False):
        """Evaluate one system over all frames using a single GPU context."""
        sim = self._make_sim(system)
        energies = []
        forces = []
        try:
            for pos, box in frames:
                EvbReporter._apply_positions(sim, pos, box)
                state = sim.context.getState(getEnergy=True,
                                             getForces=want_forces)
                energies.append(state.getPotentialEnergy().value_in_unit(
                    mm.unit.kilojoules_per_mole))
                if want_forces:
                    forces.append(
                        state.getForces(asNumpy=True).value_in_unit(
                            mm.unit.kilojoules_per_mole / mm.unit.nanometer))
        finally:
            # Release this GPU context before the next system is built.
            del sim
            gc.collect()
        return energies, forces

    def recalc_batch(self, frame_meta, frames, append):
        """Re-score a batch of frames (a whole replica's worth) and append rows.

        frames: list of (positions_nm (N,3), box_nm (3,3) or None) in nm.
        frame_meta: list of (lambda_val, replica, direction) aligned with
        frames. Each core / decomposition system is evaluated once over the
        entire batch through a single GPU context (built and torn down once),
        so context creation is amortized over all of the replica's frames.
        """
        if not self._opened:
            self._open_outputs(append)

        # The reactant/product PES and the two integration endpoints are the
        # same system regardless of window, so each is evaluated once over the
        # whole batch (one GPU context per system, reused across every frame).
        E = {}
        for name in self.core_names:
            E[name], _ = self._evaluate_system(self.systems[name], frames)
        decomp_E = {}
        for name in self.decomp_names:
            decomp_E[name], _ = self._evaluate_system(self.systems[name],
                                                      frames)

        forces_per_frame = None
        if self.report_forces:
            forces_per_frame = self._evaluate_forces_by_lambda(
                frame_meta, frames)

        for i, (lambda_val, replica, direction) in enumerate(frame_meta):
            self.E_out.write(
                EvbReporter.format_energies_row(lambda_val, E['reactant'][i],
                                                E['product'][i], E[0][i],
                                                E[1][i], replica, direction))
            if self.report_forces:
                self.F_out.write(
                    EvbReporter.format_forces_row(lambda_val, replica,
                                                  direction,
                                                  forces_per_frame[i]))
            if self.decomp_names:
                line = ", ".join(
                    [f"{decomp_E[name][i]:.10e}" for name in self.decomp_names])
                self.decomp_out.write(line + '\n')

        for s in self.out_streams:
            s.flush()

    def _evaluate_forces_by_lambda(self, frame_meta, frames):
        """Forces come from the per-lambda mixed system, which differs per
        window; build each distinct lambda's system once and evaluate the frames
        that belong to it (frames are not necessarily contiguous per lambda, as
        a replica revisits each lambda in both sweep directions)."""
        forces = [None] * len(frames)
        by_lambda = {}
        for i, (lambda_val, _, _) in enumerate(frame_meta):
            by_lambda.setdefault(lambda_val, []).append(i)
        for lambda_val, idxs in by_lambda.items():
            sub = [frames[i] for i in idxs]
            _, sub_forces = self._evaluate_system(self.systems[lambda_val],
                                                  sub,
                                                  want_forces=True)
            for j, i in enumerate(idxs):
                forces[i] = sub_forces[j]
        return forces

    def close(self):
        for s in self.out_streams:
            if not s.closed:
                s.close()


class _EvbReporterMPI:
    """Shared MPI protocol for the async EVB reporter client/server.

    Every message is a single float64 buffer whose word [0] is the message
    type. These are class constants rather than module globals so the protocol
    is self-contained on the two classes that use it.
    """

    _TAG = 7788
    _MSG_DATA = 0.0
    _MSG_BEGIN = 1.0
    _MSG_END = 2.0
    _MSG_TERMINATE = 3.0
    # Slot-free acknowledgement sent by the worker back to the master once it
    # has consumed a shared-memory ring slot (backpressure release).
    _MSG_FREE = 4.0

    @staticmethod
    def _openmm_uses_tuple():
        """Older OpenMM (HIP, < 8.2) uses the tuple describeNextReport format."""
        raw_openmm_version = version('openmm')
        match = re.match(r"^(\d+)\.(\d+)", raw_openmm_version)
        if not match:
            raise RuntimeError("Cannot parse required major.minor version from "
                               f"OpenMM: {raw_openmm_version!r}")
        return (int(match.group(1)), int(match.group(2))) < (8, 2)

    @classmethod
    def send_terminate(cls, comm, dest):
        """Tell the reporter worker to close its files and stop serving."""
        buf = np.array([cls._MSG_TERMINATE], dtype=np.float64)
        comm.Send([buf, MPI.DOUBLE], dest=dest, tag=cls._TAG)


class EvbReporterClient(_EvbReporterMPI):
    """Producer-side OpenMM reporter for asynchronous EVB energy reporting.

    Attached to the GPU sampling simulation on the master rank. Both ranks live
    on the same node, so each frame's positions/box (and velocities/forces if
    requested) are written directly into a shared-memory ring buffer
    (``MPI.Win.Allocate_shared``); only a tiny ``[_MSG_DATA, step, slot]`` control
    message is sent to wake the worker. The GPU integrator keeps stepping and the
    producer blocks only when the ring is full because the worker has fallen
    behind (backpressure), signalled by the worker's ``_MSG_FREE`` acks.
    """

    def __init__(self, comm, dest, report_interval, num_atoms, lambda_val,
                 replica, direction, append, report_forces, report_velocities,
                 shm_ring, shm_win):
        self.comm = comm
        self.dest = dest
        self.report_interval = report_interval
        self.num_atoms = num_atoms
        self.lambda_val = lambda_val
        self.replica = replica
        self.direction = direction
        self.report_forces = report_forces
        self.report_velocities = report_velocities
        self.use_tuple = self._openmm_uses_tuple()

        # Shared-memory ring buffer (queue_depth, max_len) mapped into both the
        # master and the worker; each row holds one packed frame. The slot-free
        # acks give bounded-queue backpressure (a slot may not be overwritten
        # until the worker has consumed its previous contents).
        self._ring = shm_ring
        self._win = shm_win
        self._queue_depth = shm_ring.shape[0]
        self._pending = [False] * self._queue_depth
        self._slot = 0
        self._free_buf = np.empty(3, dtype=np.float64)
        self._ctrl_buf = np.empty(3, dtype=np.float64)

        # Crude profiling: wall time this client spends handling the reporter
        # (buffer packing + MPI sends + any backpressure stalls). Accumulated
        # per window; the master reads it to separate GPU sampling from the
        # master->reporter communication overhead. comm_time is the whole-body
        # total; stall_time / send_time break it down so the derived
        # extract+pack = comm_time - stall_time - send_time is visible:
        #   stall_time  - blocked waiting for the worker to free a ring slot
        #                 (i.e. the reporter is the bottleneck, not transfer),
        #   send_time   - the shared-memory fence + tiny control Send,
        #   stall_count - number of report() calls that had to wait.
        self.comm_time = 0.0
        self.stall_time = 0.0
        self.send_time = 0.0
        self.stall_count = 0

        # Tell the worker to (re)configure its output files for this window.
        self._send_control(self._MSG_BEGIN, [
            lambda_val,
            float(replica),
            float(direction),
            1.0 if append else 0.0,
            1.0 if report_forces else 0.0,
            1.0 if report_velocities else 0.0,
        ])

    def _recv_free(self):
        """Block for one slot-free ack and mark that slot writable again."""
        self.comm.Recv([self._free_buf, MPI.DOUBLE],
                       source=self.dest,
                       tag=self._TAG)
        self._pending[int(self._free_buf[2])] = False

    def _drain_pending(self):
        """Wait for all outstanding slot-free acks (in-order consumption)."""
        for _ in range(sum(self._pending)):
            self._recv_free()

    def _send_control(self, msgtype, payload):
        t0 = time.perf_counter()
        # Control messages must not overtake outstanding data frames, so drain
        # the slot-free acks first (also releases all ring slots).
        self._drain_pending()
        buf = np.empty(1 + len(payload), dtype=np.float64)
        buf[0] = msgtype
        if payload:
            buf[1:] = payload
        self.comm.Send([buf, MPI.DOUBLE], dest=self.dest, tag=self._TAG)
        self.comm_time += time.perf_counter() - t0

    def describeNextReport(self, simulation):
        steps = self.report_interval - simulation.currentStep % self.report_interval
        if self.use_tuple:
            return (steps, True, self.report_velocities, self.report_forces,
                    False, True)
        include = ['positions']
        if self.report_velocities:
            include.append('velocities')
        if self.report_forces:
            include.append('forces')
        return {'steps': steps, 'periodic': True, 'include': include}

    def report(self, simulation, state):
        t0 = time.perf_counter()
        positions = state.getPositions(asNumpy=True).value_in_unit(
            mm.unit.nanometer)
        box = state.getPeriodicBoxVectors(asNumpy=True).value_in_unit(
            mm.unit.nanometer)

        slot = self._slot
        if self._pending[slot]:
            # Backpressure: stall until the worker frees this slot. In-order
            # consumption means the next free ack that unblocks us may name an
            # earlier slot, so keep receiving until this one is released.
            self.stall_count += 1
            t_stall = time.perf_counter()
            while self._pending[slot]:
                self._recv_free()
            self.stall_time += time.perf_counter() - t_stall
        buf = self._ring[slot]

        n = self.num_atoms
        buf[0] = self._MSG_DATA
        buf[1] = float(simulation.currentStep)
        buf[2] = float(n)
        off = 3
        buf[off:off + 9] = np.asarray(box).reshape(9)
        off += 9
        buf[off:off + 3 * n] = np.asarray(positions).reshape(-1)
        off += 3 * n
        if self.report_velocities:
            vel = state.getVelocities(asNumpy=True).value_in_unit(
                mm.unit.nanometer / mm.unit.picosecond)
            buf[off:off + 3 * n] = np.asarray(vel).reshape(-1)
            off += 3 * n
        if self.report_forces:
            frc = state.getForces(asNumpy=True).value_in_unit(
                mm.unit.kilojoules_per_mole / mm.unit.nanometer)
            buf[off:off + 3 * n] = np.asarray(frc).reshape(-1)
            off += 3 * n

        # Release fence: make the writes above visible before the worker is told
        # (via the control message) that this slot is ready to read.
        t_send = time.perf_counter()
        self._win.Sync()
        self._ctrl_buf[0] = self._MSG_DATA
        self._ctrl_buf[1] = float(simulation.currentStep)
        self._ctrl_buf[2] = float(slot)
        self.comm.Send([self._ctrl_buf, MPI.DOUBLE],
                       dest=self.dest,
                       tag=self._TAG)
        self.send_time += time.perf_counter() - t_send
        self._pending[slot] = True
        self._slot = (slot + 1) % self._queue_depth
        self.comm_time += time.perf_counter() - t0

    def finalize(self):
        """Flush this window: drain outstanding frames, then signal END."""
        self._send_control(self._MSG_END, [])


class EvbReporterServer(_EvbReporterMPI):
    """Consumer-side worker for asynchronous EVB energy reporting.

    Builds the CPU energy-evaluation simulations once (multithreaded so it
    saturates the CPU while the GPU rank samples), owns the Energies / Forces /
    Velocities / NB-decomposition output files, and serves snapshots from the
    master in order until it receives TERMINATE. Force-group / bonded
    decomposition outputs are handled synchronously on the master and are not
    produced here.
    """

    def __init__(self,
                 comm,
                 source,
                 energy_file,
                 report_interval,
                 systems,
                 topology,
                 outputstream,
                 shm_ring,
                 shm_win,
                 forming_bonds=None,
                 breaking_bonds=None,
                 force_file=None,
                 velocity_file=None,
                 cpu_threads=None):
        self.comm = comm
        self.source = source
        self.num_atoms = topology.getNumAtoms()
        # Shared-memory ring buffer mapped from the master; frames are read in
        # place (no MPI transfer of the bulk payload).
        self._ring = shm_ring
        self._win = shm_win
        self.engine = EvbReporter(
            energy_file,
            report_interval,
            systems,
            topology,
            lambda_val=0.0,
            outputstream=outputstream,
            forming_bonds=forming_bonds,
            breaking_bonds=breaking_bonds,
            forcegroup_file=None,
            force_file=force_file,
            velocity_file=velocity_file,
            append=False,
            cpu_threads=cpu_threads,
            enable_bonded_decomp=False,
            defer_open=True,
        )
        self._opened = False
        self._has_vel = False
        self._has_forces = False

        # Crude profiling accumulators (wall seconds). recv_time is time blocked
        # in MPI Recv (idle waiting for the master, i.e. reporter head-room;
        # ~0 means the reporter is saturated and gating the pipeline);
        # energy_time is time spent evaluating energies in _write_core;
        # energy_time_max is the slowest single frame; n_frames counts processed
        # data frames; n_sims is the number of CPU energy sims evaluated per
        # frame (explains the per-frame cost).
        self.recv_time = 0.0
        self.energy_time = 0.0
        self.energy_time_max = 0.0
        self.n_frames = 0
        self.n_sims = len(self.engine.simulations)

    def serve(self):
        eng = self.engine
        # Control messages are tiny now (bulk data lives in shared memory); the
        # largest is the BEGIN payload (7 words). DATA carries [_MSG_DATA, step,
        # slot]; the frame itself is read in place from the shared ring.
        ctrlbuf = np.empty(7, dtype=np.float64)
        freebuf = np.empty(3, dtype=np.float64)
        freebuf[0] = self._MSG_FREE
        freebuf[1] = 0.0
        while True:
            t_recv = time.perf_counter()
            self.comm.Recv([ctrlbuf, MPI.DOUBLE],
                           source=self.source,
                           tag=self._TAG)
            self.recv_time += time.perf_counter() - t_recv
            msgtype = ctrlbuf[0]

            if msgtype == self._MSG_TERMINATE:
                for s in eng.out_streams:
                    if not s.closed:
                        s.close()
                # Hand the accumulated reporter-side timing back to the master
                # so it can print a single consolidated profile (this rank's
                # ostream is muted).
                timing = np.array([
                    self.energy_time, self.recv_time,
                    float(self.n_frames),
                    float(self.n_sims), self.energy_time_max
                ],
                                  dtype=np.float64)
                self.comm.Send([timing, MPI.DOUBLE],
                               dest=self.source,
                               tag=self._TAG)
                break

            if msgtype == self._MSG_BEGIN:
                lam, rep, direction, append, rforces, rvel = ctrlbuf[1:7]
                eng.lambda_val = lam
                eng.replica = int(rep)
                eng.direction = int(direction)
                self._has_forces = rforces > 0.5
                self._has_vel = rvel > 0.5
                if not self._opened:
                    eng._open_outputs(append=append > 0.5)
                    self._opened = True
                continue

            if msgtype == self._MSG_END:
                for s in eng.out_streams:
                    s.flush()
                continue

            # _MSG_DATA: read the frame in place from the shared ring slot.
            slot = int(ctrlbuf[2])
            # Acquire fence: ensure the master's writes to this slot are visible.
            self._win.Sync()
            row = self._ring[slot]
            n = int(row[2])
            off = 3
            box = row[off:off + 9].reshape(3, 3)
            off += 9
            positions = row[off:off + 3 * n].reshape(n, 3)
            off += 3 * n
            velocities = None
            forces = None
            if self._has_vel:
                velocities = row[off:off + 3 * n].reshape(n, 3)
                off += 3 * n
            if self._has_forces:
                forces = row[off:off + 3 * n].reshape(n, 3)
                off += 3 * n
            t_energy = time.perf_counter()
            eng._write_core(positions, box, velocities, forces)
            dt_energy = time.perf_counter() - t_energy
            self.energy_time += dt_energy
            if dt_energy > self.energy_time_max:
                self.energy_time_max = dt_energy
            self.n_frames += 1
            for s in eng.out_streams:
                s.flush()
            # Release the slot back to the master (backpressure ack).
            freebuf[2] = float(slot)
            self.comm.Send([freebuf, MPI.DOUBLE],
                           dest=self.source,
                           tag=self._TAG)
