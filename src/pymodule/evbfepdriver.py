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
from .errorhandler import assert_msg_critical, print_exception_if_debug
from .reactionsystembuilder import EvbForceGroup

try:
    import openmm as mm
    import openmm.app as mmapp
    import openmm.unit as mmunit
except ImportError:
    pass


class EvbFepDriver:

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

        self.isothermal: bool = False
        self.isobaric: bool = False

        self.temperature = -1
        self.pressure = -1

        # Nose hoover options
        # a default of tau = 1000*dt is on the safe side, See discussion on Tdam: https://docs.lammps.org/fix_nh.html
        self.nhc_frequency = 1.0  # 1/ps,
        self.nhc_small_length = 3
        self.nhc_bulk_length = 1

        self.langevin_friction = 1.0  # 1/ps
        # See https://docs.openmm.org/latest/api-python/generated/openmm.openmm.VariableLangevinIntegrator.html
        self.langevin_tolerance = 0.001

        self.NVT_integrator = "Langevin"

        self.sample_steps = 250000
        self.write_step = 1000
        self.equil_NVT_steps = 50000
        self.equil_NPT_steps = 50000
        self.initial_equil_NVT_steps = 150000
        self.initial_equil_NPT_steps = 150000
        self.skip_initial_equil = False
        self.n_replicas = 1
        self.warmup_NVT_steps = -1
        self.warmup_NPT_steps = -1
        self.step_size = 0.001  # ps
        self.equil_step_size = 0.001  # ps
        self.minimize_every_lambda: bool = False

        self.crash_reporting_interval: int = 1
        self.constrain_H: bool = False
        self.report_forces: bool = False
        self.report_velocities: bool = False
        self.report_forcegroups: bool = False

        self.debug: bool = False
        self.save_frames: int = 1000
        self.save_crash_pdb: bool = True
        self.save_crash_xml: bool = True
        self.save_equil_traj: bool = False
        self.save_equil_pdb: bool = True
        self.xml_crash_save_interval: int = 50
        self.pdb_crash_save_interval: int = 1
        self.pdb_equil_start_temp = 10  # kelvin
        self.pdb_equil_temp_step = 50  # kelvin
        self.pdb_temperatures = []

        self.pdb = None
        self.safe_step: bool = False
        self._safe_step_batch = 10

        self.keywords = {
            "langevin_friction": {
                "type": float
            },
            "langevin_tolerance": {
                "type": float
            },
            "nhc_frequency": {
                "type": float
            },
            "nhc_small_length": {
                "type": int
            },
            "nhc_bulk_length": {
                "type": int
            },
            "temperature": {
                "type": float
            },
            "pressure": {
                "type": float
            },
            "equil_NVT_steps": {
                "type": int
            },
            "equil_NPT_steps": {
                "type": int
            },
            "initial_equil_NVT_steps": {
                "type": int
            },
            "initial_equil_NPT_steps": {
                "type": int
            },
            "skip_initial_equil": {
                "type": bool
            },
            "n_replicas": {
                "type": int
            },
            "warmup_NVT_steps": {
                "type": int
            },
            "warmup_NPT_steps": {
                "type": int
            },
            "sample_steps": {
                "type": int
            },
            "write_step": {
                "type": int
            },
            "step_size": {
                "type": float
            },
            "equil_step_size": {
                "type": float
            },
            "minimize_every_lambda": {
                "type": bool
            },
            "crash_reporting_interval": {
                "type": int
            },
            "constrain_H": {
                "type": bool
            },
            "report_forces": {
                "type": bool
            },
            "report_velocities": {
                "type": bool
            },
            "report_forcegroups": {
                "type": bool
            },
            "debug": {
                "type": bool
            },
            "save_frames": {
                "type": int
            },
            "save_crash_pdb": {
                "type": bool
            },
            "save_crash_xml": {
                "type": bool
            },
            "save_equil_traj": {
                "type": bool
            },
            "save_equil_pdb": {
                "type": bool
            },
            "xml_crash_save_interval": {
                "type": int
            },
            "pdb_crash_save_interval": {
                "type": int
            },
            "NVT_integrator": {
                "type": str
            },
            "pdb": {
                "type": str
            },
            "pdb_temperatures": {
                "type": list
            },
            "pdb_equil_temp_step": {
                "type": int
            },
            "pdb_equil_start_temp": {
                "type": int
            },
            "safe_step": {
                "type": bool
            },
        }

    def run_replicas(
        self,
        Lambda,
        configuration,
        platform,
        platform_properties,
    ):
        # todo add this to the configuration keywords

        self.platform = platform
        self.platform_properties = platform_properties
        if self.platform is None and self.platform_properties is not None:
            self.ostream.print_warning(
                "Platform properties are set, but no platform is specified. Platform properties will be ignored."
            )
            self.platform_properties = None

        for keyword, value in self.keywords.items():
            if keyword in configuration:
                if (not isinstance(configuration[keyword], value["type"])
                        and not ((isinstance(configuration[keyword], int)
                                  and value["type"] == float))):
                    raise ValueError(
                        f"Configuration option {keyword} should be of type {value['type']}"
                    )
                else:
                    setattr(self, keyword, configuration[keyword])
                    self.ostream.print_info(
                        f"{keyword}: {getattr(self, keyword)}")
            else:
                self.ostream.print_info(
                    f"{keyword}: {getattr(self, keyword)} (default)")
        self.forming_bonds = configuration.get("forming_bonds")
        self.breaking_bonds = configuration.get("breaking_bonds")

        self.ostream.flush()
        assert_msg_critical('openmm' in sys.modules,
                            'openmm is required for EvbFepDriver.')
        systems = configuration["systems"]
        topology = configuration["topology"]
        initial_positions = configuration["initial_positions"]

        cwd = Path.cwd()
        self.run_folder = cwd / configuration["run_folder"]
        self.data_folder = cwd / configuration["data_folder"]
        self.Lambda = Lambda
        self.systems = systems
        self.topology = topology

        if self.temperature > 0:
            self.isothermal = True

        if self.pressure > 0:
            self.isobaric = True

        if self.pdb is not None:
            if len(self.pdb_temperatures) == 0:
                self.pdb_temperatures = list(
                    np.arange(self.pdb_equil_start_temp, self.temperature,
                              self.pdb_equil_temp_step))

            if self.pdb_temperatures[-1] != self.temperature:
                self.pdb_temperatures.append(self.temperature)
            self.ostream.print_info(
                f"Set PDB equilibration temperatures to {self.pdb_temperatures}"
            )

        self.ostream.flush()

        assert (
            self.sample_steps %
            self.write_step == 0), "write_step must be a factor of sample_steps"
        assert (self.sample_steps >= 2 *
                self.write_step), "sample_steps must be at least 2*write_step"
        assert self.n_replicas >= 1, "n_replicas must be at least 1"

        # Each replica performs a forward (l: 0->1) and a backward (l: 1->0)
        # sweep, hence the factor 2 * n_replicas.
        self.passes_per_replica = 2
        n_passes = self.passes_per_replica * self.n_replicas
        self.total_snapshots = self.sample_steps / self.write_step * len(
            self.Lambda) * n_passes
        self.ostream.print_info(f"Lambda: {np.array(self.Lambda)}")
        info = f"Replicas: {self.n_replicas} (forward + backward each), total lambda sweeps: {n_passes}\n"
        info += f"Total lambda points: {len(self.Lambda)}, NVT equilibration steps: {self.equil_NVT_steps}, NPT equilibration steps: {self.equil_NPT_steps}, total sample steps: {self.sample_steps}, write step: {self.write_step}, step size: {self.step_size}\n"
        info += f"Snapshots per lambda: {self.sample_steps / self.write_step}, snapshots to be recorded: {self.total_snapshots}\n"
        info += f"System time per snapshot: {self.step_size * self.write_step} ps, system time per frame: {self.step_size * self.sample_steps} ps, total system time: {self.step_size * self.sample_steps * len(self.Lambda) * n_passes} ps"
        self.ostream.print_info(
            f"Ensemble info: Isobaric {self.isobaric}, Isothermal {self.isothermal}"
        )
        self.total_steps = 0
        self.run_steps = 0
        # Per-lambda equilibration + sampling, summed over every lambda window in
        # a single sweep, then scaled by the number of sweeps.
        self.pass_steps = (self.equil_NVT_steps + self.sample_steps) * len(
            self.Lambda)
        if self.isobaric:
            self.pass_steps += self.equil_NPT_steps * len(self.Lambda)
        self.total_steps += self.pass_steps * n_passes

        # Initial equilibration / PDB warm-up. This is only executed once, for
        # the l == 0 window and only when skip_initial_equil is not set, so
        # mirror exactly the step counts that _initial_equilibrate will run.
        if not self.skip_initial_equil and 0 in self.Lambda:
            if self.pdb is None:
                self.total_steps += self.initial_equil_NVT_steps
                if self.isobaric:
                    self.total_steps += self.initial_equil_NPT_steps
            else:
                for T in self.pdb_temperatures:
                    if T == self.temperature:
                        self.total_steps += self.initial_equil_NVT_steps
                        if self.isobaric:
                            self.total_steps += self.initial_equil_NPT_steps
                    else:
                        self.total_steps += (self.warmup_NVT_steps
                                             if self.warmup_NVT_steps > 0 else
                                             self.initial_equil_NVT_steps)
                        if self.isobaric:
                            self.total_steps += (self.warmup_NPT_steps if
                                                 self.warmup_NPT_steps > 0 else
                                                 self.initial_equil_NPT_steps)

        self.ostream.print_info(info)
        self.ostream.flush()

        # Global timer spanning all replicas / sweeps for total-progress ETA.
        self.timer = Timer(self.total_steps)

        # Single trajectory file, appended across all sweeps. Frames stay
        # row-aligned with Energies.csv (each row tagged with replica/direction).
        self.traj_roporter = mmapp.XTCReporter(
            str(self.data_folder / "trajectory.xtc"),
            self.write_step,
        )
        # Only the very first written window truncates the output files and
        # writes their headers; everything afterwards appends.
        self._first_write = True
        self.timer.start()

        positions = initial_positions * 0.1
        velocities = None

        # Initial equilibration, performed once on the l == 0 system. The
        # resulting positions / velocities seed the first forward sweep.
        if not self.skip_initial_equil and 0 in self.Lambda:
            system = systems[0]
            if self.constrain_H:
                self.ostream.print_info(
                    "Constraining all bonds involving H atoms")
                system = self._constrain_H_bonds(system)
            state = self._initial_equilibrate(system, positions)
            positions = state.getPositions()
            velocities = state.getVelocities()
            if self.constrain_H:
                self.ostream.print_info(
                    "Removing constraints involving H atoms")
                for _ in range(system.getNumConstraints()):
                    system.removeConstraint(0)
        elif self.skip_initial_equil:
            self.ostream.print_info(
                "Skipping initial equilibration because skip_initial_equil is set to True."
            )

        forward_order = list(self.Lambda)
        backward_order = list(reversed(self.Lambda))
        for replica in range(self.n_replicas):
            # Forward sweep (l: 0 -> 1), direction = 0.
            positions, velocities = self.run_FEP(forward_order, replica, 0,
                                                 positions, velocities)
            # Backward sweep (l: 1 -> 0), direction = 1, seeded by the final
            # state of the forward sweep.
            positions, velocities = self.run_FEP(backward_order, replica, 1,
                                                 positions, velocities)
        self.ostream.flush()

    def run_FEP(self, lambda_order, replica, direction, positions, velocities):
        """Run a single lambda sweep (one replica, one direction).

        Equilibrates then samples each lambda window in ``lambda_order``,
        carrying positions/velocities forward, and returns the final
        positions/velocities so the caller can chain the next sweep.
        """
        dir_label = "forward" if direction == 0 else "backward"
        self.ostream.print_info(
            f"Replica {replica} {dir_label} sweep (lambda {lambda_order[0]} -> {lambda_order[-1]})"
        )
        self.ostream.flush()

        # Per-sweep timer for the ETA within this replica's FEP.
        pass_timer = Timer(self.pass_steps)
        pass_timer.start()
        pass_run_steps = 0
        window_steps = self.equil_NVT_steps + self.sample_steps
        if self.isobaric:
            window_steps += self.equil_NPT_steps

        for l in lambda_order:
            system = self.systems[l]

            if self.constrain_H:
                self.ostream.print_info(
                    "Constraining all bonds involving H atoms")
                system = self._constrain_H_bonds(system)

            equil_state = self._equilibrate(system, l, positions, velocities)

            if self.constrain_H:
                self.ostream.print_info(
                    "Removing constraints involving H atoms")
                for _ in range(system.getNumConstraints()):
                    system.removeConstraint(0)

            state = self._sample(system, l, equil_state, replica, direction)
            positions = state.getPositions()
            velocities = state.getVelocities()

            pass_run_steps += window_steps
            self._print_progress(pass_timer, pass_run_steps, replica, dir_label,
                                 l)

        return positions, velocities

    def _print_progress(self, pass_timer, pass_run_steps, replica, dir_label,
                        l):
        # Per-sweep ETA (this replica/direction) and global ETA (all sweeps).
        pass_remaining, _ = pass_timer.calculate_remaining(pass_run_steps)
        total_remaining, _ = self.timer.calculate_remaining(self.run_steps)
        pass_str = pass_timer.get_time_str(pass_remaining)
        total_str = self.timer.get_time_str(total_remaining)
        msg = (
            f"Replica {replica} {dir_label} l={l}: "
            f"sweep {pass_run_steps}/{self.pass_steps} steps, ETA {pass_str} | "
            f"total {self.run_steps}/{self.total_steps} steps, ETA {total_str}")
        self.ostream.print_info(msg)
        self.ostream.flush()

    def _initial_equilibrate(self, system, positions):
        simulation = self._get_simulation(system, self.equil_step_size)
        simulation.context.setPositions(positions)

        platformname = simulation.context.getPlatform()
        self.ostream.print_info(
            f"Running FEP on platform: {platformname.getName()}")
        self.ostream.flush()

        self._minimize(simulation)
        if self.save_equil_traj:
            equil_traj_reporter = mmapp.XTCReporter(
                str(self.run_folder / "equil_traj_initial.xtc"),
                self.write_step,
                enforcePeriodicBox=True,
            )
            simulation.reporters.append(equil_traj_reporter)
        if self.pdb is None:
            if self.isobaric:
                barostat = self._get_barostat(simulation)
                barostat.setFrequency(0)
                self._safe_step(simulation, self.initial_equil_NVT_steps,
                                "initial NVT equilibration")
                barostat.setFrequency(25)
                self._safe_step(simulation, self.initial_equil_NPT_steps,
                                "initial NPT equilibration")
            else:
                self._safe_step(simulation, self.initial_equil_NVT_steps,
                                "initial equilibration")
        else:
            self.ostream.print_info(
                f"Perfoming PDB warmup with T-vector {np.array(self.pdb_temperatures)}"
            )
            self.ostream.flush()
            for T in self.pdb_temperatures:
                simulation.integrator.setTemperature(T)
                if T == self.temperature:
                    NVT_steps = self.initial_equil_NVT_steps
                    NPT_steps = self.initial_equil_NPT_steps
                else:
                    NVT_steps = self.warmup_NVT_steps if self.warmup_NVT_steps > 0 else self.initial_equil_NVT_steps
                    NPT_steps = self.warmup_NPT_steps if self.warmup_NPT_steps > 0 else self.initial_equil_NPT_steps

                if self.isobaric:
                    barostat = self._get_barostat(simulation)
                    barostat.setFrequency(0)
                    self.ostream.print_info(
                        f"Running PDB warmup NVT equilibration T = {T}")
                    self.ostream.flush()
                    self._safe_step(simulation, NVT_steps,
                                    f"PDB warmup NVT equilibration T = {T}")
                    barostat.setFrequency(25)
                    self.ostream.print_info(
                        f"Running PDB warmup NPT equilibration T = {T}")
                    self.ostream.flush()
                    self._safe_step(simulation, NPT_steps,
                                    f"PDB warmup NPT equilibration T = {T}")
                else:
                    self.ostream.print_info(
                        f"Running PDB warmup NVT equilibration T = {T}")
                    self.ostream.flush()
                    self._safe_step(simulation, NVT_steps,
                                    f"PDB warmup NVT equilibration T = {T}")
                if self.debug:
                    self._save_state(
                        simulation,
                        f"warmup_equil_state_{T:.1f}K",
                        xml=False,
                        chk=False,
                        cif=True,
                    )

        equil_state = simulation.context.getState(
            getPositions=True,
            getVelocities=True,
            getEnergy=True,
            # getParameters=True,
            # getParameterDerivatives=True,
            # getIntegratorParameters=True,
            enforcePeriodicBox=True,
        )
        self._save_state(
            simulation,
            "equil_state_initial",
            xml=False,
            chk=True,
            cif=False,
        )
        return equil_state

    def _equilibrate(self, system, l, positions, velocities=None):
        simulation = self._get_simulation(system, self.equil_step_size)
        simulation.context.setPositions(positions)
        if velocities is not None:
            simulation.context.setVelocities(velocities)
        if self.minimize_every_lambda:
            self._minimize(simulation)

        if self.debug:
            self._save_state(simulation, f"pre_equil_{l:.3f}")
            simulation.reporters.append(
                mmapp.PDBReporter(
                    str(self.run_folder / f"traj_equil_{l:.3f}.pdb"),
                    self.write_step,
                    enforcePeriodicBox=True,
                ))
            f_file = str(self.run_folder / f"forces_equil_{l:.3f}.csv")
            v_file = str(self.run_folder / f"velocities_equil_{l:.3f}.csv")
            g_file = str(self.run_folder / f"forcegroups_equil_{l:.3f}.csv")
            simulation.reporters.append(
                EvbReporter(
                    str(self.run_folder / f"energies_equil_{l:.3f}.csv"),
                    self.write_step,
                    self.systems,
                    self.topology,
                    l,
                    self.ostream,
                    forming_bonds=self.forming_bonds,
                    breaking_bonds=self.breaking_bonds,
                    force_file=f_file,
                    velocity_file=v_file,
                    forcegroup_file=g_file,
                    append=False,
                ))

        sz = self.equil_step_size * mmunit.picoseconds
        simulation.integrator.setStepSize(sz)
        # self.ostream.print_info(f"Equilibration with step size {sz}")

        if self.save_equil_traj:
            equil_traj_reporter = mmapp.XTCReporter(
                str(self.run_folder / f"equil_traj_{l:.3f}.xtc"),
                self.write_step,
                enforcePeriodicBox=True,
            )
            simulation.reporters.append(equil_traj_reporter)

        if self.isobaric:
            barostat = self._get_barostat(simulation)
            barostat.setFrequency(0)
            self._safe_step(simulation, self.equil_NVT_steps,
                            "NVT equilibration")
            barostat.setFrequency(25)
            self._safe_step(simulation, self.equil_NPT_steps,
                            "NPT equilibration")
        else:
            self._safe_step(simulation, self.equil_NVT_steps, "equilibration")

        equil_state = simulation.context.getState(
            getPositions=True,
            getVelocities=True,
            # getForces=True,
            getEnergy=True,
            # getParameters=True,
            # getParameterDerivatives=True,
            # getIntegratorParameters=True,
            enforcePeriodicBox=True,
        )
        if self.save_equil_pdb:

            self._save_state(
                simulation,
                f"equil_state_{l:.3f}",
                chk=False,
                xml=False,
                cif=True,
            )

        return equil_state

    def _minimize(self, simulation, filename=None):
        self.ostream.print_info("Minimizing energy")
        self.ostream.flush()
        simulation.minimizeEnergy()
        positions = simulation.context.getState(
            getPositions=True, enforcePeriodicBox=True).getPositions()
        if filename is not None:
            mmapp.PDBxFile.writeFile(
                self.topology,
                np.array(positions.value_in_unit(mm.unit.angstrom)),
                open(self.run_folder / f"minim_{filename}.cif", "w"),
            )

    @staticmethod
    def _get_barostat(simulation):
        return [
            force for force in simulation.system.getForces()
            if isinstance(force, mm.MonteCarloBarostat)
            or isinstance(force, mm.MonteCarloFlexibleBarostat)
        ][0]

    def _sample(self, system, l, initial_state, replica=0, direction=0):
        simulation = self._get_simulation(system, self.equil_step_size)
        simulation.reporters.append(self.traj_roporter)

        sz = self.step_size * mmunit.picoseconds
        simulation.integrator.setStepSize(sz)
        simulation.context.setState(initial_state)
        # Truncate + write headers only on the very first sampled window across
        # all replicas/sweeps; append for every window afterwards.
        append = not self._first_write
        self._first_write = False
        state_reporter = self._get_data_reporter(
            self.data_folder,
            "Data_combined.csv",
            append,
        )

        evb_reporter = self._get_evb_reporter(
            self.data_folder,
            l,
            append=append,
            replica=replica,
            direction=direction,
        )
        simulation.reporters.append(state_reporter)
        simulation.reporters.append(evb_reporter)
        self.ostream.flush()
        states = self._safe_step(simulation, self.sample_steps, "sampling")
        return states[-1]

    def _get_simulation(self, system, step_size):
        if self.isothermal:
            if self.NVT_integrator == "variable Langevin":
                integrator = mm.VariableLangevinIntegrator(
                    self.temperature * mmunit.kelvin,  #type: ignore
                    self.langevin_friction / mmunit.picosecond,  #type: ignore
                    self.langevin_tolerance,
                )
            elif self.NVT_integrator == "Langevin":
                integrator = mm.LangevinMiddleIntegrator(
                    self.temperature * mmunit.kelvin,  # type: ignore
                    self.langevin_friction / mmunit.picosecond,  # type: ignore
                    step_size * mmunit.picoseconds,
                )
            elif self.NVT_integrator == "Nose-Hoover":
                if system.getNumParticles() > 100:
                    chain_length = self.nhc_bulk_length
                else:
                    chain_length = self.nhc_small_length

                integrator = mm.NoseHooverIntegrator(
                    self.temperature * mmunit.kelvin,
                    self.nhc_frequency / mmunit.picosecond,
                    step_size * mmunit.picoseconds,
                    chain_length,
                )
            else:
                assert False, "NVT-integrator should be either 'variable Langevin', 'Langevin' or 'Nose-Hoover'"
        else:
            integrator = mm.VerletIntegrator(step_size)

        if self.platform is not None:
            simulation = mmapp.Simulation(
                self.topology,
                system,
                integrator,
                platform=mm.Platform.getPlatformByName(self.platform),
                platformProperties=self.platform_properties,
            )
        else:
            simulation = mmapp.Simulation(
                self.topology,
                system,
                integrator,
            )

        return simulation

    def _get_data_reporter(
        self,
        folder,
        name,
        append=False,
        write_step=-1,
    ):
        if write_step == -1:
            write_step = self.write_step
        state_reporter = mmapp.StateDataReporter(
            str(folder / name),
            self.write_step,
            step=True,
            potentialEnergy=True,
            kineticEnergy=True,
            temperature=True,
            volume=True,
            density=True,
            append=append,
        )
        return state_reporter

    def _get_evb_reporter(
        self,
        folder,
        l,
        name_suffix="",
        append=False,
        write_step=-1,
        replica=0,
        direction=0,
    ):
        if self.report_forces or self.debug:
            f_file = str(folder / f"Forces{name_suffix}.csv")
        else:
            f_file = None
        if self.report_velocities or self.debug:
            v_file = str(folder / f"Velocities{name_suffix}.csv")
        else:
            v_file = None
        if self.report_forcegroups or self.debug:
            g_file = str(folder / f"ForceGroups{name_suffix}.csv")
        else:
            g_file = None

        evb_reporter = EvbReporter(
            str(folder / f"Energies{name_suffix}.csv"),
            self.write_step,
            self.systems,
            self.topology,
            l,
            self.ostream,
            forming_bonds=self.forming_bonds,
            breaking_bonds=self.breaking_bonds,
            forcegroup_file=g_file,
            velocity_file=v_file,
            force_file=f_file,
            append=append,
            replica=replica,
            direction=direction,
        )
        return evb_reporter

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

    def _save_state(self,
                    simulation,
                    name,
                    xml=True,
                    chk=True,
                    cif=True,
                    pdb=False):
        if xml:
            chk_file = str(self.run_folder / f"{name}.chk")
            simulation.saveCheckpoint(chk_file)
        if chk:
            xml_file = str(self.run_folder / f"{name}.xml")
            simulation.saveState(xml_file)
        if cif or pdb:
            state = simulation.context.getState(getPositions=True,
                                                enforcePeriodicBox=True)
            positions = np.array(state.getPositions().value_in_unit(
                mm.unit.angstrom))
            if pdb:
                mmapp.PDBFile.writeFile(
                    self.topology,
                    positions,
                    open(self.run_folder / f"{name}.pdb", "w"),
                )
            if cif:
                mmapp.PDBxFile.writeFile(
                    self.topology,
                    positions,
                    open(self.run_folder / f"{name}.cif", "w"),
                )

    def _safe_step(self, simulation, steps, name=""):
        if not self.safe_step:
            self.run_steps += steps
            simulation.step(steps)
            state = simulation.context.getState(
                getPositions=True,
                getVelocities=True,
                getForces=self.report_forces or self.debug,
                getEnergy=True,
                enforcePeriodicBox=True,
            )
            return [state]

        self.ostream.print_info(
            f"Running {name} for {steps} steps for total time: {steps*self.step_size} ps with step size {self.step_size} ps"
        )
        self.ostream.flush()
        states = []
        # potwarning = False
        if steps % self._safe_step_batch != 0:
            self.ostream.print_warning(
                f"Steps {steps} is not a multiple of safe step batch {self._safe_step_batch}, rounding down to {steps - steps % self._safe_step_batch}"
            )
            steps -= steps % self._safe_step_batch
        for i in range(steps // self._safe_step_batch):
            try:
                self.run_steps += self._safe_step_batch
                simulation.step(self._safe_step_batch)
            except Exception as e:
                self.ostream.print_warning(
                    f"Error during simulation step {i}: {e}")
                self.ostream.flush()
                self._save_states(states, simulation, i)
                raise

            #todo only querry velocities and forces if reporting them or debugging
            state = simulation.context.getState(
                getPositions=True,
                getVelocities=True,
                getForces=self.report_forces or self.debug,
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

            if self.save_crash_xml and j % self.xml_crash_save_interval == 0:
                xml_name = f"state_step_{j}_{step_num}"
                with (path / f"{xml_name}.xml").open("w") as f:
                    f.write(mm.XmlSerializer.serialize(state))

            if self.save_crash_pdb and j % self.pdb_crash_save_interval == 0:
                pdb_name = f"state_step_{step_num}"
                positions = np.array(state.getPositions().value_in_unit(
                    mm.unit.angstrom))

                # Make sure that openmm PDB writing doesn't crash when encountering too large values
                positions = np.clip(positions, -9999998, 99999998)
                positions[positions == np.inf] = 99999998
                positions[positions == -np.inf] = -9999998
                positions[positions == np.nan] = 0
                positions = np.nan_to_num(positions,
                                          posinf=99999998,
                                          neginf=-9999998)

                mmapp.PDBFile.writeFile(
                    self.topology,
                    positions,
                    open(self.run_folder / f"{pdb_name}.pdb", "w"),
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
            try:
                simulation.context.setState(state)
                for k, fg in enumerate(EvbForceGroup):
                    fg_state = simulation.context.getState(
                        getEnergy=True,
                        groups=set([fg.value]),
                    )
                    energy = fg_state.getPotentialEnergy()
                    energy = energy.value_in_unit(mmunit.kilojoule_per_mole)
                    energies[j, k + 4] = energy
            except Exception:
                self.ostream.print_warning(
                    "Encountered error while saving forcegroups, continuing without forcegroups"
                )

        # Combine all saved PDB files into one and remove the sigle ones
        if self.save_crash_pdb:

            output_file = "combined_crash.pdb"
            pdb_pattern = "state_step_*.pdb"
            pdb_files = sorted(self.run_folder.glob(pdb_pattern))
            with (self.data_folder / output_file).open("w") as outfile:
                for model_number, pdb_file in enumerate(pdb_files, start=1):
                    outfile.write(f"MODEL     {model_number}\n")
                    with pdb_file.open("r") as infile:
                        for line in infile:
                            if line.startswith(
                                ('ATOM', 'HETATM', 'TER',
                                 'END')) or (model_number == 1
                                             and line.startswith("CRYST1")
                                             ):  # Skip headers/footers
                                outfile.write(line)
                    outfile.write("ENDMDL\n")
                    pdb_file.unlink()

        header = "step, kinetic, potential, volume,"
        header += EvbForceGroup.get_header()
        np.savetxt(
            self.data_folder / "crash_energies.csv",
            energies,
            delimiter=",",
            header=header,
            fmt="%.5e",
        )


class Timer:

    def __init__(self, total_steps):
        self.start_time = time.time()
        self.total_steps = total_steps
        self.times = []
        self.step = 0

    def start(self):
        self.start_time = time.time()

    def print_time_str(self, step, ostream, message=""):
        time_remaining, elapsed_time = self.calculate_remaining(step)
        total_time = sum(self.times)
        remain_str = self.get_time_str(time_remaining)
        elapsed_str = self.get_time_str(elapsed_time)
        total_elapsed_str = self.get_time_str(total_time)
        print_str = f"Estimated remaining time: {remain_str}, elapsed time since last timing: {elapsed_str}, total elapsed time: {total_elapsed_str}"
        if message:
            print_str = f"{message} {print_str}"
        print_str += f" (step {step} of {self.total_steps} total steps)"
        ostream.print_info(print_str)
        ostream.flush()

    def calculate_remaining(self, step):
        end_time = time.time()
        self.step = step
        elapsed_time = end_time - self.start_time

        self.times.append(elapsed_time)

        if step > 0:
            avg_time_per_step = sum(self.times) / step
        else:
            avg_time_per_step = 0.0
        remaining_steps = max(self.total_steps - self.step, 0)
        estimated_time_remaining = avg_time_per_step * remaining_steps

        self.start_time = time.time()  # reset start time for next step
        return estimated_time_remaining, elapsed_time

    def get_time_str(self, time):

        hours, minutes, seconds = self.convert_seconds(time)
        time_str = f"{hours:02}:{minutes:02}:{seconds:02}"
        return time_str

    @staticmethod
    def convert_seconds(seconds):
        hours = seconds // 3600  # calculate hours
        seconds %= 3600  # update seconds remaining
        minutes = seconds // 60  # calculate minutes
        seconds %= 60  # update seconds remaining
        return int(hours), int(minutes), int(seconds)
