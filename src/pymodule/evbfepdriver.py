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
from pathlib import Path
import numpy as np
import time
import sys
import os
import glob

from .veloxchemlib import mpi_master
from .outputstream import OutputStream
from .evbreporter import EvbReporter
from .errorhandler import assert_msg_critical
from .evbsystembuilder import EvbForceGroup

try:
    import openmm as mm
    import openmm.app as mmapp
    import openmm.unit as mmunit
except ImportError:
    pass


class EvbFepDriver():

    def __init__(self, comm=None, ostream=None):
        '''
        Initialize the EVB driver class.
        '''

        assert_msg_critical('openmm' in sys.modules,
                            'openmm is required for EvbFepDriver.')

        if comm is None:
            comm = MPI.COMM_WORLD

        if ostream is None:
            if comm.Get_rank() == mpi_master():
                ostream = OutputStream(sys.stdout)
            else:
                ostream = OutputStream(None)

        # output stream
        self.ostream = ostream

        # mpi information
        self.comm = comm
        self.rank = self.comm.Get_rank()
        self.nodes = self.comm.Get_size()

        self.systems: dict = None
        self.topology: mmapp.Topology = None
        self.Lambda: list = None
        self.crash_reporting_interval: int = 1

        self.constrain_H: bool = False
        self.report_forces: bool = False
        self.report_velocities: bool = False
        self.report_forcegroups: bool = True
        self.debug: bool = False
        self.save_frames: int = 1000

        self.isothermal: bool = False
        self.isobaric: bool = False

    def run_FEP(
        self,
        equil_NVT_steps,
        equil_NPT_steps,
        total_sample_steps,
        write_step,
        initial_equil_NVT_steps,
        initial_equil_NPT_steps,
        step_size,
        equil_step_size,
        Lambda,
        configuration,
        platform,
    ):
        self.equil_NVT_steps = equil_NVT_steps
        self.equil_NPT_steps = equil_NPT_steps
        self.total_sample_steps = total_sample_steps
        self.write_step = write_step
        self.initial_equil_NVT_steps = initial_equil_NVT_steps
        self.initial_equil_NPT_steps = initial_equil_NPT_steps
        self.step_size = step_size
        self.equil_step_size = equil_step_size
        self.platform = platform

        assert_msg_critical('openmm' in sys.modules,
                            'openmm is required for EvbFepDriver.')

        systems = configuration["systems"]
        topology = configuration["topology"]
        temperature = configuration["temperature"]
        initial_positions = configuration["initial_positions"]

        cwd = Path.cwd()
        self.run_folder = cwd / configuration["run_folder"]
        self.data_folder = cwd / configuration["data_folder"]
        self.Lambda = Lambda
        self.systems = systems
        self.topology = topology

        temperature = configuration.get('temperature', -1)
        pressure = configuration.get('pressure', -1)

        if temperature > 0:
            self.isothermal = True

        if pressure > 0:
            self.isobaric = True

        self.ostream.print_info(
            f"Ensemble info: Isobaric {self.isobaric}, Isothermal {self.isothermal}"
        )
        self.ostream.flush()

        self.integrator_temperature = temperature * mmunit.kelvin  #type: ignore
        self.integrator_friction_coeff = 0.5 / mmunit.picosecond

        assert (total_sample_steps % write_step == 0
                ), "write_step must be a factor of total_sample_steps"
        assert (total_sample_steps >= 2 *
                write_step), "total_sample_steps must be at least 2*write_step"

        self.total_snapshots = total_sample_steps / write_step * len(
            self.Lambda)
        self.ostream.print_info(f"Lambda: {self.Lambda}")
        info = f"Total lambda points: {len(self.Lambda)}, NVT equilibration steps: {equil_NVT_steps}, NPT equiliberation steps: {equil_NPT_steps}, total sample steps: {total_sample_steps}, write step: {write_step}, step size: {step_size}\n"
        info += f"Snapshots per lambda: {total_sample_steps / write_step}, snapshots to be recorded: {self.total_snapshots}\n"
        info += f"System time per snapshot: {step_size * write_step} ps, system time per frame: {step_size * total_sample_steps} ps, total system time: {step_size * total_sample_steps * len(self.Lambda)} ps"
        self.ostream.print_info(info)
        self.ostream.flush()

        timer = Timer(len(self.Lambda))

        self.traj_roporter = mmapp.XTCReporter(
            str(self.data_folder / "trajectory.xtc"),
            write_step,
        )
        timer.start()
        positions = initial_positions

        for i, l in enumerate(self.Lambda):
            if l > 0:
                estimated_time_remaining = timer.calculate_remaining(i)
                time_estimate_str = ", " + timer.get_time_str(
                    estimated_time_remaining)
                self.ostream.print_info(f"lambda = {l}" + time_estimate_str)
            else:
                self.ostream.print_info(f"lambda = {l}")
            system = systems[l]

            if self.constrain_H:
                self.ostream.print_info(
                    "Constraining all bonds involving H atoms")
                system = self._constrain_H_bonds(system)

            equil_state = self._minimize_and_equilibrate(system, l, positions)

            if self.constrain_H:
                info = "Removing constraints involving H atoms"
                self.ostream.print_info(info)
                for i in range(system.getNumConstraints()):
                    system.removeConstraint(0)

            state = self._sample(system, l, equil_state)
            positions = state.getPositions()
        self.ostream.flush()

    def _minimize_and_equilibrate(self, system, l, positions):
        equil_simulation = self._get_simulation(system, self.equil_step_size)
        equil_simulation.context.setPositions(positions)

        if l == 0:
            platformname = equil_simulation.context.getPlatform()
            self.ostream.print_info(
                f"Running FEP on platform: {platformname.getName()}")
            self.ostream.flush()
        self.ostream.print_info("Minimizing energy")
        self.ostream.flush()
        equil_simulation.minimizeEnergy()
        minim_positions = equil_simulation.context.getState(
            getPositions=True, enforcePeriodicBox=True).getPositions()

        mmapp.PDBFile.writeFile(
            self.topology,
            np.array(minim_positions.value_in_unit(mm.unit.angstrom)),
            open(self.run_folder / f"minim_{l:.3f}.pdb", "w"),
        )
        if self.debug:
            self._save_state(equil_simulation, f"minim_{l:.3f}")
            equil_simulation.reporters.append(
                mmapp.PDBReporter(
                    str(self.run_folder / f"traj_equil_{l:.3f}.pdb"),
                    self.write_step,
                    enforcePeriodicBox=True,
                ))
            f_file = str(self.run_folder / f"forces_equil_{l:.3f}.csv")
            v_file = str(self.run_folder / f"velocities_equil_{l:.3f}.csv")
            g_file = str(self.run_folder / f"forcegroups_equil_{l:.3f}.csv")
            equil_simulation.reporters.append(
                EvbReporter(
                    str(self.run_folder / f"energies_equil_{l:.3f}.csv"),
                    self.write_step,
                    self.systems[0],
                    self.systems[1],
                    self.topology,
                    l,
                    self.ostream,
                    force_file=f_file,
                    velocity_file=v_file,
                    forcegroup_file=g_file,
                    append=False,
                ))

        sz = self.equil_step_size * mmunit.picoseconds
        equil_simulation.integrator.setStepSize(sz)
        self.ostream.print_info(f"Equilibration with step size {sz}")

        equil_reporter = mmapp.StateDataReporter(
            str(self.run_folder / f"equil_data_{l:.3f}.csv"),
            1,
            step=True,
            potentialEnergy=True,
            kineticEnergy=True,
            temperature=True,
            volume=True,
            density=True,
            append=False,
        )
        equil_simulation.reporters.append(equil_reporter)
        barostat = [
            force for force in equil_simulation.system.getForces()
            if isinstance(force, mm.MonteCarloBarostat)
        ][0]
        if l == 0:
            if self.isobaric:
                barostat.setFrequency(0)
                self.ostream.print_info(
                    f"Running initial NVT equilibration for {self.initial_equil_NVT_steps} steps"
                )
                self.ostream.flush()

                self._safe_step(equil_simulation, self.initial_equil_NVT_steps)
                self.ostream.print_info(
                    f"Running initial NPT equilibration for {self.initial_equil_NPT_steps} steps"
                )
                self.ostream.flush()
                barostat.setFrequency(25)
                self._safe_step(equil_simulation, self.initial_equil_NPT_steps)
            else:
                barostat.setFrequency(0)
                self.ostream.print_info(
                    f"Running initial equilibration for {self.initial_equil_NVT_steps+self.initial_equil_NPT_steps} steps"
                )
                self.ostream.flush()
                equil_simulation.integrator.setIntegrationForceGroups(
                    EvbForceGroup.integration_force_groups())
                self._safe_step(
                    equil_simulation,
                    self.initial_equil_NVT_steps + self.initial_equil_NPT_steps)

        if self.isobaric:
            self.ostream.print_info(
                f"Running NVT equilibration for {self.equil_NVT_steps} steps")
            self.ostream.flush()
            barostat.setFrequency(0)
            self._safe_step(equil_simulation, self.equil_NVT_steps)

            self.ostream.print_info(
                f"Running NPT equilibration for {self.equil_NPT_steps} steps")
            self.ostream.flush()
            barostat.setFrequency(25)
            self._safe_step(equil_simulation, self.equil_NPT_steps)
        else:
            self.ostream.print_info(
                f"Running equilibration for {self.equil_NVT_steps+self.equil_NPT_steps} steps"
            )
            self.ostream.flush()
            self._safe_step(equil_simulation,
                            self.equil_NVT_steps + self.equil_NPT_steps)

        equil_state = equil_simulation.context.getState(
            getPositions=True,
            getVelocities=True,
            getForces=True,
            getEnergy=True,
            getParameters=True,
            getParameterDerivatives=True,
            getIntegratorParameters=True,
            enforcePeriodicBox=True,
        )

        if self.debug:
            self._save_state(
                equil_simulation,
                f"equil_state_{l:.3f}",
                xml=False,
            )

        return equil_state

    def _sample(self, system, l, initial_state):
        run_simulation = self._get_simulation(system, self.equil_step_size)
        run_simulation.reporters.append(self.traj_roporter)
        sz = self.step_size * mmunit.picoseconds
        run_simulation.integrator.setStepSize(sz)
        run_simulation.context.setState(initial_state)
        if l == 0:
            append = False
        else:
            append = True

        state_reporter = mmapp.StateDataReporter(
            str(self.data_folder / "Data_combined.csv"),
            self.write_step,
            step=True,
            potentialEnergy=True,
            kineticEnergy=True,
            temperature=True,
            volume=True,
            density=True,
            append=append,
        )
        run_simulation.reporters.append(state_reporter)

        if self.report_forces or self.debug:
            f_file = str(self.data_folder / f"Forces.csv")
        else:
            f_file = None
        if self.report_velocities or self.debug:
            v_file = str(self.data_folder / f"Velocities.csv")
        else:
            v_file = None
        if self.report_forcegroups or self.debug:
            g_file = str(self.data_folder / f"ForceGroups.csv")
        else:
            g_file = None

        evb_reporter = EvbReporter(
            str(self.data_folder / "Energies.csv"),
            self.write_step,
            self.systems[0],
            self.systems[1],
            self.topology,
            l,
            self.ostream,
            forcegroup_file=g_file,
            velocity_file=v_file,
            force_file=f_file,
            append=append,
            debug=self.debug,
        )
        run_simulation.reporters.append(evb_reporter)

        self.ostream.print_info(
            f"Running sampling for {self.total_sample_steps} steps with step size {run_simulation.integrator.getStepSize()}"
        )
        self.ostream.flush()
        states = self._safe_step(run_simulation, self.total_sample_steps)
        return states[-1]

    def _get_simulation(self, system, step_size):
        if self.isothermal:
            integrator = mm.LangevinMiddleIntegrator(
                self.integrator_temperature,
                self.integrator_friction_coeff,
                step_size * mmunit.picoseconds,
            )
        else:
            integrator = mm.VerletIntegrator(step_size)
        integrator.setIntegrationForceGroups(
            EvbForceGroup.integration_force_groups())

        if self.platform is not None:
            simulation = mmapp.Simulation(
                self.topology,
                system,
                integrator,
                mm.Platform.getPlatformByName(self.platform),
            )
        else:
            simulation = mmapp.Simulation(
                self.topology,
                system,
                integrator,
            )
        return simulation

    def _constrain_H_bonds(self, system):
        harm_bond_forces = [
            force for force in system.getForces()
            if isinstance(force, mm.HarmonicBondForce)
        ]
        constrained_H_bonds = []
        count = 0
        for harmbond in harm_bond_forces:
            for i in range(harmbond.getNumBonds()):
                particle1, particle2, length, k = harmbond.getBondParameters(i)
                if (system.getParticleMass(particle1).value_in_unit(
                        mmunit.dalton) - 1.007947 < 0.01
                        or system.getParticleMass(particle2).value_in_unit(
                            mmunit.dalton) - 1.007947 < 0.01):
                    H_bond = sorted((particle1, particle2))
                    if H_bond not in constrained_H_bonds:
                        constrained_H_bonds.append(H_bond)
                        system.addConstraint(particle1, particle2, length)
                        count += 1
        self.ostream.print_info(f"Constrained {count} bonds involving H atoms")
        return system

    def _save_state(self, simulation, name, xml=True, chk=True):
        if xml:
            chk_file = str(self.run_folder / f"{name}.chk")
            simulation.saveCheckpoint(chk_file)
        if chk:
            xml_file = str(self.run_folder / f"{name}.xml")
            simulation.saveState(xml_file)

    def _safe_step(self, simulation, steps):
        states = []
        potwarning = False
        for i in range(steps):
            try:
                simulation.step(1)
            except Exception as e:
                self.ostream.print_warning(
                    f"Error during simulation step {i}: {e}")
                self.ostream.flush()
                self._save_states(states, simulation, i)
                raise e

            state = simulation.context.getState(
                getPositions=True,
                getVelocities=True,
                getForces=True,
                getEnergy=True,
                enforcePeriodicBox=True,
            )
            kin = state.getKineticEnergy()
            kin = kin.value_in_unit(mmunit.kilojoule_per_mole)
            pot = state.getPotentialEnergy()
            pot = pot.value_in_unit(mmunit.kilojoule_per_mole)

            states.append(state)
            if len(states) > self.save_frames:
                states.pop(0)
        return states

    def _save_states(self, states, simulation, step):
        self.ostream.print_info(f"Saving last {len(states)} states")
        self.ostream.flush()
        cwd = Path.cwd()
        path = cwd / self.run_folder

        energies = np.zeros((len(states), len(EvbForceGroup) + 4))

        for j, state in enumerate(states):
            step_num = step - len(states) + j
            name = f"state_step_{j}_{step_num}"
            with open(path / f"{name}.xml", "w") as f:
                f.write(mm.XmlSerializer.serialize(state))

            minim_positions = state.getPositions()
            mmapp.PDBFile.writeFile(
                self.topology,
                np.array(minim_positions.value_in_unit(mm.unit.angstrom)),
                open(self.run_folder / f"{name}.pdb", "w"),
            )

            kin = state.getKineticEnergy()
            kin = kin.value_in_unit(mmunit.kilojoule_per_mole)
            pot = state.getPotentialEnergy()
            pot = pot.value_in_unit(mmunit.kilojoule_per_mole)
            vol = state.getPeriodicBoxVolume()
            vol = vol.value_in_unit(mmunit.nanometer**3)

            energies[j, 0] = step_num
            energies[j, 1] = kin
            energies[j, 2] = pot
            energies[j, 3] = vol

            simulation.context.setState(state)
            for k, fg in enumerate(EvbForceGroup):
                fg_state = simulation.context.getState(
                    getEnergy=True,
                    groups=set([fg.value]),
                )
                energy = fg_state.getPotentialEnergy()
                energy = energy.value_in_unit(mmunit.kilojoule_per_mole)
                energies[j, k + 4] = energy

        # Combine all saved PDB files into one and remove the sigle ones
        input_folder = self.run_folder  # replace this
        output_file = "combined_crash.pdb"
        pdb_pattern = "state_step_*.pdb"
        pdb_files = sorted(glob.glob(os.path.join(input_folder, pdb_pattern)))
        with open(self.run_folder / output_file, 'w') as outfile:
            for model_number, pdb_file in enumerate(pdb_files, start=1):
                outfile.write(f"MODEL     {model_number}\n")
                with open(pdb_file, 'r') as infile:
                    for line in infile:
                        if line.startswith(('ATOM', 'HETATM', 'TER',
                                            'END')):  # Skip headers/footers
                            outfile.write(line)
                outfile.write("ENDMDL\n")
                os.remove(pdb_file)

        header = "step, kinetic, potential, volume,"
        header += EvbForceGroup.get_header()
        np.savetxt(
            path / f"crash_energies.csv",
            energies,
            delimiter=",",
            header=header,
            fmt="%.5e",
        )


class Timer:

    def __init__(self, total_iterations):
        self.start_time = time.time()
        self.total_iterations = total_iterations
        self.times = []

    def start(self):
        self.start_time = time.time()

    def calculate_remaining(self, iteration):
        end_time = time.time()
        elapsed_time = end_time - self.start_time

        self.times.append(elapsed_time)

        avg_time_per_iteration = sum(self.times) / len(self.times)
        remaining_iterations = self.total_iterations - iteration
        estimated_time_remaining = avg_time_per_iteration * remaining_iterations

        self.start_time = time.time()  # reset start time for next iteration
        return estimated_time_remaining

    def get_time_str(self, estimated_time_remaining):
        hours, minutes, seconds = self.convert_seconds(estimated_time_remaining)
        time_str = "estimated time remaining: "
        if hours > 0:
            time_str += "{} hour(s), ".format(hours)
        if minutes > 0:
            time_str += "{} minute(s), ".format(minutes)
        time_str += "{} second(s)".format(seconds)
        return time_str

    @staticmethod
    def convert_seconds(seconds):
        hours = seconds // 3600  # calculate hours
        seconds %= 3600  # update seconds remaining
        minutes = seconds // 60  # calculate minutes
        seconds %= 60  # update seconds remaining
        return int(hours), int(minutes), int(seconds)
