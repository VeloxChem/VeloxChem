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

from .veloxchemlib import mpi_master
from .outputstream import OutputStream
from .evbreporter import EvbReporter
from .errorhandler import assert_msg_critical

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

        assert_msg_critical('openmm' in sys.modules, 'openmm is required for EvbFepDriver.')

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

        self.constrain_H: bool = True

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
        configuration
    ):

        assert_msg_critical('openmm' in sys.modules, 'openmm is required for EvbFepDriver.')

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

        assert (total_sample_steps % write_step == 0), "write_step must be a factor of total_sample_steps"
        assert (total_sample_steps >= 2 * write_step), "total_sample_steps must be at least 2*write_step"

        self.total_snapshots = total_sample_steps / write_step * len(self.Lambda)
        self.ostream.print_info(f"Lambda: {self.Lambda}")
        self.ostream.print_info(f"Total lambda points: {len(self.Lambda)}")
        self.ostream.print_info(f"Snapshots per lambda: {total_sample_steps / write_step}")
        self.ostream.print_info(f"Snapshots to be recorded: {self.total_snapshots}")
        self.ostream.print_info(f"Total simulation steps: {(total_sample_steps + equilibration_steps) * len(self.Lambda) + lambda_0_equilibration_steps}")
        self.ostream.print_info(f"System time per snapshot: {step_size * write_step} ps")
        self.ostream.print_info(f"System time per frame: {step_size * total_sample_steps} ps")
        self.ostream.print_info(f"Total system time: {step_size * total_sample_steps * len(self.Lambda)} ps",)
        self.ostream.flush()
        integrator_temperature = temperature * mmunit.kelvin  #type: ignore
        integrator_friction_coeff = 1 / mmunit.picosecond

        timer = Timer(len(self.Lambda))

        traj_roporter = mmapp.XTCReporter(
            str(self.data_folder / "trajectory.xtc"),
            write_step,
        )
        timer.start()
        for i, l in enumerate(self.Lambda):
            system = systems[l]
            integrator = mm.LangevinMiddleIntegrator(
                integrator_temperature,
                integrator_friction_coeff,
                equil_step_size * mmunit.picoseconds,
            )

            if l>0:
                estimated_time_remaining = timer.calculate_remaining(i)
                time_estimate_str = ", " + timer.get_time_str(estimated_time_remaining)
                self.ostream.print_info(f"lambda = {l}" + time_estimate_str)
            else:
                self.ostream.print_info(f"lambda = {l}")
            
            constrained_H_bonds = []
            if self.constrain_H:
                harm_bond_forces = [
                    force for force in system.getForces() if isinstance(force, mm.HarmonicBondForce)
                ]

                count = 0
                for harmbond in harm_bond_forces:
                    for i in range(harmbond.getNumBonds()):
                        particle1, particle2, length, k = harmbond.getBondParameters(i)
                        if (system.getParticleMass(particle1).value_in_unit(mmunit.dalton) - 1.007947 < 0.01 or
                                system.getParticleMass(particle2).value_in_unit(mmunit.dalton) - 1.007947 < 0.01):
                            H_bond = sorted((particle1, particle2))
                            if H_bond not in constrained_H_bonds:
                                constrained_H_bonds.append(H_bond)
                                system.addConstraint(particle1, particle2, length)
                                count += 1
                self.ostream.print_info(f"Constrained {count} bonds involving H atoms ")

            simulation = mmapp.Simulation(
                topology,
                system,
                integrator,
            )

            simulation.context.setPositions(initial_positions)
            simulation.reporters.append(mmapp.XTCReporter(
                str(self.run_folder / f"minim_{l:.3f}.xtc"),
                write_step,
            ))

            self.ostream.print_info("Minimizing energy")
            self.ostream.flush()
            simulation.minimizeEnergy()

            if l == 0:
                simulation.integrator.setStepSize(initial_equil_step_size * mmunit.picoseconds)
                self.ostream.print_info(f"Running initial equilibration with step size {simulation.integrator.getStepSize()}")
                self.ostream.flush()
                simulation.step(lambda_0_equilibration_steps)
                timer.start()

            # equilibrate
            simulation.integrator.setStepSize(equil_step_size * mmunit.picoseconds)
            self.ostream.print_info(f"Running equilibration with step size {simulation.integrator.getStepSize()}")
            self.ostream.flush()
            simulation.step(equilibration_steps)

            equil_positions = simulation.context.getState(getPositions=True).getPositions()

            if self.constrain_H:
                self.ostream.print_info("Removing constraints")
                for i in range(system.getNumConstraints()):
                    system.removeConstraint(0)

            # Updating the system requires recreating the simulation, and integrators can only be bound to one context
            integrator = mm.LangevinMiddleIntegrator(
                integrator_temperature,
                integrator_friction_coeff,
                step_size * mmunit.picoseconds,
            )
            runsimulation = mmapp.Simulation(
                topology,
                system,
                integrator,
            )
            runsimulation.context.setPositions(equil_positions)
            runsimulation.reporters.append(traj_roporter)

            if l == 0:
                append = False
            else:
                append = True
            runsimulation.reporters.append(EvbReporter(
                str(self.data_folder / "Energies.dat"),
                write_step,
                systems["reactant"],
                systems["product"],
                systems[0],
                systems[1],
                topology,
                l,
                self.ostream,
                append=append,
            ))

            runsimulation.reporters.append(
                mmapp.StateDataReporter(
                    str(self.data_folder / "Data_combined.dat"),
                    write_step,
                    step=True,
                    potentialEnergy=True,
                    kineticEnergy=True,
                    temperature=True,
                    volume=True,
                    density=True,
                    append=append,
                ))

            self.ostream.print_info(f"Running sampling with step size {runsimulation.integrator.getStepSize()}")
            self.ostream.flush()    
            runsimulation.step(total_sample_steps)
            state = runsimulation.context.getState(getPositions=True)
            positions = state.getPositions()
            if np.any(np.array(positions.value_in_unit(mmunit.nanometer)) > 100):
                self.ostream.print_info("Warning: Some positions are larger than 100 nm, system is probably unstable")

        self.ostream.print_info("Merging output files")
        self.ostream.flush()

#?Utility class for predicting approximate runtime of the simulation, do we already have dependencies for this so that this can be scrapped?
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

        # print(f"Iteration: {iteration}")
        # print(f"Elapsed time: {elapsed_time}")
        # print(f"Avg time per iteration: {avg_time_per_iteration}")
        # print(f"Remaining iterations: {remaining_iterations}")
        # print(f"Estimated time remaining: {estimated_time_remaining}")
        # print(f"Times list: {self.times}")
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
