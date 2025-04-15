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

    def run_FEP(
        self,
        equilibration_steps,
        total_sample_steps,
        write_step,
        lambda_0_equilibration_steps,
        step_size,
        equil_step_size,
        initial_equil_step_size,
        Lambda,
        configuration,
        platform,
    ):

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

        assert (total_sample_steps % write_step == 0
                ), "write_step must be a factor of total_sample_steps"
        assert (total_sample_steps >= 2 *
                write_step), "total_sample_steps must be at least 2*write_step"

        self.total_snapshots = total_sample_steps / write_step * len(
            self.Lambda)
        simulation_steps = (total_sample_steps + equilibration_steps) * len(
            self.Lambda) + lambda_0_equilibration_steps
        self.ostream.print_info(f"Lambda: {self.Lambda}")
        info = f"Total lambda points: {len(self.Lambda)}, equilibration steps: {equilibration_steps}, total sample steps: {total_sample_steps}, write step: {write_step}, step size: {step_size}\n"
        info += f"Snapshots per lambda: {total_sample_steps / write_step}, snapshots to be recorded: {self.total_snapshots}, total simulation steps: {simulation_steps}\n"
        info += f"System time per snapshot: {step_size * write_step} ps, system time per frame: {step_size * total_sample_steps} ps, total system time: {step_size * total_sample_steps * len(self.Lambda)} ps"
        self.ostream.print_info(info)
        self.ostream.flush()

        integrator_temperature = temperature * mmunit.kelvin  #type: ignore
        integrator_friction_coeff = 1 / mmunit.picosecond

        timer = Timer(len(self.Lambda))

        traj_roporter = mmapp.XTCReporter(
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

            equil_integrator = mm.LangevinMiddleIntegrator(
                integrator_temperature,
                integrator_friction_coeff,
                initial_equil_step_size * mmunit.picoseconds,
            )
            equil_integrator.setIntegrationForceGroups(
                EvbForceGroup.integration_force_groups())
            if platform is not None:
                equil_simulation = mmapp.Simulation(
                    topology,
                    system,
                    equil_integrator,
                    mm.Platform.getPlatformByName(platform),
                )
            else:
                equil_simulation = mmapp.Simulation(
                    topology,
                    system,
                    equil_integrator,
                )
            equil_simulation.context.setPositions(positions)

            if i == 0:
                platformname = equil_simulation.context.getPlatform()
                self.ostream.print_info(
                    f"Running FEP on platform: {platformname.getName()}")
                self.ostream.flush()

            if self.constrain_H:
                self.ostream.print_info(
                    "Constraining all bonds involving H atoms")
                system = self._constrain_H_bonds(system)

            self.ostream.print_info("Minimizing energy")
            self.ostream.flush()
            equil_simulation.minimizeEnergy()
            minim_positions = equil_simulation.context.getState(
                getPositions=True, enforcePeriodicBox=True).getPositions()
            mmapp.PDBFile.writeFile(
                topology,
                np.array(minim_positions.value_in_unit(mm.unit.angstrom)),
                open(self.run_folder / f"minim_{l:.3f}.pdb", "w"),
            )
            if self.debug:
                self._save_state(equil_simulation, f"minim_{l:.3f}")

            sz = equil_step_size * mmunit.picoseconds
            equil_simulation.integrator.setStepSize(sz)
            self.ostream.print_info(
                f"Running equilibration with step size {equil_simulation.integrator.getStepSize()}"
            )
            self.ostream.flush()
            if self.debug:
                equil_simulation.reporters.append(
                    mmapp.PDBReporter(
                        str(self.run_folder / f"traj_equil_{l:.3f}.pdb"),
                        write_step,
                        enforcePeriodicBox=True,
                    ))
                f_file = str(self.run_folder / f"forces_equil_{l:.3f}.csv")
                v_file = str(self.run_folder / f"velocities_equil_{l:.3f}.csv")
                g_file = str(self.run_folder / f"forcegroups_equil_{l:.3f}.csv")
                equil_simulation.reporters.append(
                    EvbReporter(
                        str(self.run_folder / f"energies_equil_{l:.3f}.csv"),
                        write_step,
                        systems[0],
                        systems[1],
                        topology,
                        l,
                        self.ostream,
                        force_file=f_file,
                        velocity_file=v_file,
                        forcegroup_file=g_file,
                        append=False,
                    ))

            states = self._safe_step(equil_simulation, equilibration_steps)

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

            run_integrator = mm.LangevinMiddleIntegrator(
                integrator_temperature,
                integrator_friction_coeff,
                equil_step_size * mmunit.picoseconds,
            )
            run_integrator.setIntegrationForceGroups(
                EvbForceGroup.integration_force_groups())

            if self.constrain_H:
                info = "Removing constraints involving H atoms"
                self.ostream.print_info(info)
                for i in range(system.getNumConstraints()):
                    system.removeConstraint(0)

            if platform is not None:
                run_simulation = mmapp.Simulation(
                    topology,
                    system,
                    run_integrator,
                    mm.Platform.getPlatformByName(platform),
                )
            else:
                run_simulation = mmapp.Simulation(
                    topology,
                    system,
                    run_integrator,
                )

            run_simulation.reporters.append(traj_roporter)
            sz = step_size * mmunit.picoseconds
            run_simulation.integrator.setStepSize(sz)
            run_simulation.context.setState(equil_state)
            if l == 0:
                append = False
            else:
                append = True

            state_reporter = mmapp.StateDataReporter(
                str(self.data_folder / "Data_combined.csv"),
                write_step,
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
                write_step,
                systems[0],
                systems[1],
                topology,
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
                f"Running sampling with step size {run_simulation.integrator.getStepSize()}"
            )
            self.ostream.flush()
            states = self._safe_step(run_simulation, total_sample_steps)

            positions = states[-1].getPositions()
        self.ostream.flush()

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
            )
            kin = state.getKineticEnergy()
            kin = kin.value_in_unit(mmunit.kilojoule_per_mole)
            pot = state.getPotentialEnergy()
            pot = pot.value_in_unit(mmunit.kilojoule_per_mole)

            # self.ostream.print_info(
            #     f"Step {i}, kinetic energy: {kin:.5f} kJ/mol, potential energy: {pot:.5f} kJ/mol"
            # )
            self.ostream.flush()

            states.append(state)
            if len(states) > self.save_frames:
                states.pop(0)

            if pot > 0 and not potwarning:
                
                self.ostream.print_warning(
                    f"Potential energy is positive: {pot:.5f} kJ/mol."
                )
                potwarning = True
                # self._save_states(states, simulation, i)
                # raise RuntimeError(
                #     f"Potential energy is positive: {pot:.5f} kJ/mol. Simulation crashed"
                # )

        return states

    def _save_states(self, states, simulation, step):
        self.ostream.print_info(f"Saving last {len(states)} states")
        self.ostream.flush()
        cwd = Path.cwd()
        path = cwd / self.run_folder

        energies = np.zeros((len(states), len(EvbForceGroup) + 3))

        for j, state in enumerate(states):
            step_num = step - len(states) + j
            with open(path / f"state_step_{step_num}.xml", "w") as f:
                f.write(mm.XmlSerializer.serialize(state))

            minim_positions = state.getPositions()
            mmapp.PDBFile.writeFile(
                self.topology,
                np.array(minim_positions.value_in_unit(mm.unit.angstrom)),
                open(self.run_folder / f"state_step_{step_num}.pdb", "w"),
            )

            kin = state.getKineticEnergy()
            kin = kin.value_in_unit(mmunit.kilojoule_per_mole)
            pot = state.getPotentialEnergy()
            pot = pot.value_in_unit(mmunit.kilojoule_per_mole)

            energies[j, 0] = step_num
            energies[j, 1] = kin
            energies[j, 2] = pot

            simulation.context.setState(state)
            energies[j, 0]
            for k, fg in enumerate(EvbForceGroup):
                fg_state = simulation.context.getState(
                    getEnergy=True,
                    groups=set([fg.value]),
                )
                energy = fg_state.getPotentialEnergy()
                energy = energy.value_in_unit(mmunit.kilojoule_per_mole)
                energies[j, k + 3] = energy
        header = "step,kinetic,potential,"
        header += ",".join([fg.name for fg in EvbForceGroup])
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
