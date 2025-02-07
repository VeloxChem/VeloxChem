import sys
from pathlib import Path
# from EVB.timer import Timer
# from EVB.system_builder import System_builder

import numpy as np
import openmm as mm
import openmm.app as mmapp
import openmm.unit as mmunit

from mpi4py import MPI

from .veloxchemlib import mpi_master
from .outputstream import OutputStream
from .evbreporter import EvbReporter

import time


class FepDriver():

    def __init__(self, comm=None, ostream=None):
        '''
        Initialize the EVB driver class.
        '''

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
        equilliberation_steps,
        total_sample_steps,
        write_step,
        lambda_0_equilliberation_steps,
        step_size,
        equil_step_size,
        initial_equil_step_size,
        Lambda,
        configuration
    ):
        systems = configuration["systems"]
        topology = configuration["topology"]
        temperature = configuration["temperature"]
        initial_positions = configuration["initial_positions"]
        run_folder = configuration["run_folder"]
        data_folder = configuration["data_folder"]

        self.run_folder = run_folder
        self.data_folder = data_folder
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
        self.ostream.print_info(f"Total simulation steps: {(total_sample_steps + equilliberation_steps) * len(self.Lambda) + lambda_0_equilliberation_steps}")
        self.ostream.print_info(f"System time per snapshot: {step_size * write_step} ps")
        self.ostream.print_info(f"System time per frame: {step_size * total_sample_steps} ps")
        self.ostream.print_info(f"Total system time: {step_size * total_sample_steps * len(self.Lambda)} ps",)
        integrator_temperature = temperature * mmunit.kelvin  #type: ignore
        integrator_friction_coeff = 1 / mmunit.picosecond
        # integrator_step_size = 0.001 * unit.picoseconds

        # if self.init_positions == None:
        #     self.init_positions = rea_ff_gen.molecule.get_coordinates_in_angstrom()*0.1

        timer = Timer(len(self.Lambda))
        estimated_time_remaining = None

        for i, l in enumerate(self.Lambda):
            timer.start()
            system = systems[l]
            integrator = mm.LangevinMiddleIntegrator(
                integrator_temperature,
                integrator_friction_coeff,
                equil_step_size * mmunit.picoseconds,
            )

            if estimated_time_remaining:
                time_estimate_str = ", " + timer.get_time_str(estimated_time_remaining)
            else:
                time_estimate_str = ""

            self.ostream.print_info(f"lambda = {l}" + time_estimate_str)
            constrained_H_bonds = []
            if self.constrain_H:
                harm_bond_forces = [
                    force for force in system.getForces() if isinstance(force, mm.HarmonicBondForce)
                ]

                count = 0
                for harmbond in harm_bond_forces:
                    for i in range(harmbond.getNumBonds()):
                        particle1, particle2, length, k = harmbond.getBondParameters(i)

                        if system.getParticleMass(particle1).value_in_unit(
                                mmunit.dalton) == 1.007947 or system.getParticleMass(particle2).value_in_unit(
                                    mmunit.dalton) == 1.007947: #todo more robust hydrogen check
                            H_bond = sorted((particle1, particle2))
                            if H_bond not in constrained_H_bonds:
                                constrained_H_bonds.append(H_bond)
                                system.addConstraint(particle1, particle2, length)
                                count += 1
                self.ostream.print_info(f"Constrained {count} bonds involving H atoms ")
                # with open("constrained_.xml", mode="w", encoding="utf-8") as output:
                #     output.write(mm.XmlSerializer.serialize(system))


            # platform = mm.Platform.getPlatformByName("CUDA")
            # properties = {'Precision': 'single'}
            simulation = mmapp.Simulation(
                topology,
                system,
                integrator,
                # platform,
                # properties,
            )

            simulation.context.setPositions(initial_positions)
            # simulation.context.setParameter("lambda", l)
            simulation.reporters.append(mmapp.XTCReporter(
                f"{self.run_folder}/minim_{l:.3f}.xtc",
                write_step,
            ))

            self.ostream.print_info("Minimizing energy")
            simulation.minimizeEnergy()

            if l == 0:
                simulation.integrator.setStepSize(initial_equil_step_size * mmunit.picoseconds)
                self.ostream.print_info(f"Running initial equilliberation with step size {simulation.integrator.getStepSize()}")
                self.ostream.flush()
                simulation.step(lambda_0_equilliberation_steps)
                timer.start()

            # Equiliberate
            # todo write these reporters on my own
            # todo add lambda value to the reporter
            simulation.integrator.setStepSize(equil_step_size * mmunit.picoseconds)
            self.ostream.print_info(f"Running equilliberation with step size {simulation.integrator.getStepSize()}")
            self.ostream.flush()
            simulation.step(equilliberation_steps)

            equil_positions = simulation.context.getState(getPositions=True).getPositions()

            if self.constrain_H:
                self.ostream.print_info("Removing constraints")
                for i in range(system.getNumConstraints()):
                    #Removing a constraint will change the indexing of the constraints, so always remove the first one until there are none left
                    system.removeConstraint(0)
                # with open("unconstrained_.xml", mode="w", encoding="utf-8") as output:
                #     output.write(mm.XmlSerializer.serialize(system))

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

            if l == 0:
                append = False
            else:
                append = True

            # Append option gives precision mismatch error, investigate or file bug report
            runsimulation.reporters.append(mmapp.XTCReporter(
                f"{self.run_folder}/traj{l:.3f}.xtc",
                write_step,
            ))

            # runsimulation.reporters.append(mmapp.XTCReporter(
            #     f"{self.data_folder}/trajectory.xtc",
            #     write_step,
            #     append=append,
            # ))

            runsimulation.reporters.append(EvbReporter(
                f"{self.data_folder}/Energies.dat",
                write_step,
                systems["reactant"],
                systems["product"],
                systems[0],
                systems[1],
                topology,
                l,
                append=append,
            ))

            runsimulation.reporters.append(
                mmapp.StateDataReporter(
                    f"{self.data_folder}/Data_combined.dat",
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
            estimated_time_remaining = timer.stop_and_calculate(i + 1)

        self.ostream.print_info("Merging output files")
        # self.merge_traj_pdb()
        # self.merge_state()
        self.ostream.flush()

    def merge_traj_pdb(self):
        self.ostream.print_info("merging pdb files")
        output = ""
        traj_path = Path().cwd() / self.data_folder / "traj_combined.pdb"
        with open(traj_path, "w", encoding="utf-8") as file:
            file.write(output)
        frame = 1
        crystline = None
        for l in self.Lambda:
            self.ostream.print_info(f"Lambda = {l}")
            filename = f"{self.run_folder}/traj{l:.3f}.pdb"

            with open(filename, "r", encoding="utf-8") as file:
                file_contents = file.read()

            # self.ostream.print_info(file_contents)
            for line in file_contents.split("\n"):
                parts = line.split()
                if len(parts) == 0:
                    continue
                if parts[0] == "REMARK":
                    continue
                if parts[0] == "MODEL":
                    line = f"MODEL{frame: >9}\n"
                    frame += 1
                if parts[0] == "CONECT":
                    continue
                if parts[0] == "END":
                    continue
                if parts[0] == "CRYST1":
                    crystline = line
                    continue

                output += line + "\n"
            # Write every frame seperately already to the file and empty the output string, otherwise the output string will become too large to handle nicely
            with open(traj_path, "a", encoding="utf-8") as file:
                file.write(output)
            output = ""
        if crystline:
            output += crystline + "\n"
        output += "END"
        with open(traj_path, "a", encoding="utf-8") as file:
            file.write(output)

#?Utility class for predicting approximate runtime of the simulation, do we already have dependencies for this so that this can be scrapped?
class Timer:

    def __init__(self, total_iterations):
        self.start_time = time.time()
        self.total_iterations = total_iterations
        self.times = []

    def start(self):
        self.start_time = time.time()

    def stop_and_calculate(self, iteration):
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

    def print_progress(self, iteration, message):
        progress = iteration / self.total_iterations
        bar_length = 20
        block = int(round(bar_length * progress))
        text = "\rProgress: [{0}] {1}% {2}".format("#" * block + "-" * (bar_length - block), round(progress * 100, 2),
                                                   message)
        sys.stdout.write(text)
        sys.stdout.flush()

    @staticmethod
    def convert_seconds(seconds):
        hours = seconds // 3600  # calculate hours
        seconds %= 3600  # update seconds remaining
        minutes = seconds // 60  # calculate minutes
        seconds %= 60  # update seconds remaining
        return int(hours), int(minutes), int(seconds)
