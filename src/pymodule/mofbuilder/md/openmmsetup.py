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

try:
    from openmm.app import (
        GromacsGroFile,
        GromacsTopFile,
        PDBFile,
        PDBReporter,
        NoCutoff,
        Simulation,
        PME,
        HBonds,
        StateDataReporter,
    )
    from openmm import LangevinIntegrator, MonteCarloBarostat
    from openmm.unit import (
        atmosphere,
        picosecond,
        picoseconds,
        kelvin,
        femtoseconds,
        nanometer,
    )
except ImportError:
    pass
from mpi4py import MPI
import sys
from typing import Optional, List, Any, Union
from ...outputstream import OutputStream
from ...veloxchemlib import mpi_master
from ...errorhandler import assert_msg_critical


class OpenmmSetup:
    """Set up and execute GROMACS/OPENMM simulations for molecular dynamics.

    Handles system initialization, energy minimization, equilibration in NVT/NPT ensembles,
    and trajectory/log file writing, with parallel output support via MPI.

    Attributes:
        comm (Any): MPI communicator.
        rank (int): MPI process rank.
        nodes (int): Number of MPI processes.
        ostream (OutputStream): Output stream for printing/logging.
        system_pbc (bool): Whether to use periodic boundary conditions.
        gro_file (str): Path to GROMACS .gro coordinates file.
        top_file (str): Path to GROMACS .top topology file.
        temperature_K (float): Simulation temperature in Kelvin.
        timestep_fs (float): Integration timestep in femtoseconds.
        gro: Loaded GROMACS structure file (GromacsGroFile).
        top: Loaded GROMACS topology file (GromacsTopFile).
        positions: Atomic positions.
        barostat_added (bool): Tracks if barostat has been added for NPT.
        simulation: The OpenMM Simulation object (initialized on demand).
        temperature (openmm.unit.Quantity): Temperature (set during pipeline).
        timestep (openmm.unit.Quantity): Timestep (set during pipeline).
        system: OpenMM System object (set during pipeline).
    """

    def __init__(
        self,
        gro_file: str,
        top_file: str,
        temperature_K: float = 300,
        timestep_fs: float = 1,
        comm: Optional[Any] = None,
        ostream: Optional[OutputStream] = None,
    ):
        """Initializes the OpenmmSetup simulation context.

        Args:
            gro_file (str): Path to input .gro coordinate file.
            top_file (str): Path to input .top topology file.
            temperature_K (float, optional): Simulation temperature in Kelvin. Defaults to 300.
            timestep_fs (float, optional): Integration timestep (femtoseconds). Defaults to 1.
            comm (Any, optional): MPI communicator. Defaults to MPI.COMM_WORLD.
            ostream (OutputStream, optional): Output stream for logging.

        """
        assert_msg_critical(
            "openmm" in sys.modules,
            "OpenMM is required for OpenmmSetup.")

        self.comm = comm or MPI.COMM_WORLD
        self.rank = self.comm.Get_rank()
        self.nodes = self.comm.Get_size()
        self.ostream = ostream or OutputStream(
            sys.stdout if self.rank == mpi_master() else None
        )
        self.system_pbc = True  # assume PBC for MD

        self.gro_file = gro_file
        self.top_file = top_file
        self.temperature_K = temperature_K
        self.timestep_fs = timestep_fs
        self.simulation = None  # type: Optional[Simulation]

        # Load GROMACS files
        self.gro = GromacsGroFile(gro_file)
        self.top = GromacsTopFile(
            top_file, periodicBoxVectors=self.gro.getPeriodicBoxVectors(), includeDir="."
        )

        self.positions = self.gro.positions
        self.barostat_added = False

    def _steps_from_ps(self, ps_time: float) -> int:
        """Convert a time in picoseconds to number of simulation steps.

        Args:
            ps_time (float): Duration (picoseconds).

        Returns:
            int: Number of steps for current timestep.
        """
        return int(ps_time / self.timestep.value_in_unit(picoseconds))

    def run_em(self, output_prefix: str = "em", whole_traj_file: Optional[str] = None) -> None:
        """Performs energy minimization using OpenMM and writes output.

        Writes a minimized structure (PDB) to disk and logs potential/temperature.
        Optionally appends coordinates to a trajectory PDB.

        Args:
            output_prefix (str, optional): Prefix for output and log files.
            whole_traj_file (Optional[str], optional): If provided, append structure to this PDB.
        """
        print("Starting energy minimization...")
        if self.simulation is None:
            integrator = LangevinIntegrator(
                self.temperature, 1 / picosecond, self.timestep
            )
            self.simulation = Simulation(self.top.topology, self.system, integrator)
            self.simulation.context.setPositions(self.positions)
        self.simulation.minimizeEnergy()
        self.positions = self.simulation.context.getState(
            getPositions=True
        ).getPositions()

        # Add log
        self.simulation.reporters.append(
            StateDataReporter(
                f"{output_prefix}.log", 1, step=True, potentialEnergy=True, temperature=True
            )
        )
        # Save minimized structure
        em_file = f"{output_prefix}_minimized.pdb"
        PDBFile.writeFile(
            self.simulation.topology, self.positions, open(em_file, "w")
        )
        print("Energy minimization done.")

        # Append to whole trajectory if requested
        if whole_traj_file:
            PDBFile.writeFile(
                self.simulation.topology, self.positions, open(whole_traj_file, "a")
            )

    def _run_md(
        self,
        mode: str = "nvt",
        time_ps: float = 100,
        output_prefix: str = "traj",
        record_interval_ps: float = 1,
        whole_traj_file: Optional[str] = None,
    ) -> None:
        """Runs an MD segment in either NVT or NPT ensemble.

        Args:
            mode (str, optional): 'nvt' or 'npt'.
            time_ps (float, optional): Duration to run (picoseconds).
            output_prefix (str, optional): Output file prefix.
            record_interval_ps (float, optional): Trajectory/log frame interval (picoseconds).
            whole_traj_file (Optional[str], optional): If given, appends frames to this PDB.
        """
        nsteps = self._steps_from_ps(time_ps)
        steps_per_frame = self._steps_from_ps(record_interval_ps)
        if mode.lower() == "npt" and not self.barostat_added:
            barostat = MonteCarloBarostat(1 * atmosphere, self.temperature)
            self.system.addForce(barostat)
            self.barostat_added = True
        if self.simulation is None:
            integrator = LangevinIntegrator(
                self.temperature, 1 / picosecond, self.timestep
            )
            self.simulation = Simulation(self.top.topology, self.system, integrator)
            self.simulation.context.setPositions(self.positions)

        # Set reporters
        if whole_traj_file:
            self.simulation.reporters.append(
                PDBReporter(whole_traj_file, steps_per_frame)
            )
        else:
            self.simulation.reporters.append(
                PDBReporter(f"{output_prefix}.pdb", steps_per_frame)
            )

        self.simulation.reporters.append(
            StateDataReporter(
                f"{output_prefix}.log",
                steps_per_frame,
                step=True,
                potentialEnergy=True,
                temperature=True,
                density=(mode.lower() == "npt"),
            )
        )

        print(f"Running {mode.upper()} for {time_ps} ps ...")
        self.simulation.step(nsteps)
        self.positions = self.simulation.context.getState(
            getPositions=True
        ).getPositions()
        print(f"{mode.upper()} done.")

    def run_nvt(
        self,
        nvt_time: float = 100,
        output_prefix: str = "nvt",
        record_interval: float = 1,
        whole_traj_file: Optional[str] = None,
    ) -> None:
        """Run an NVT (constant volume) MD segment.

        Args:
            nvt_time (float, optional): Simulation time (ps).
            output_prefix (str, optional): Output file prefix.
            record_interval (float, optional): Frame/log interval (ps).
            whole_traj_file (Optional[str], optional): Append all frames to this PDB if provided.
        """
        self._run_md(
            mode="nvt",
            time_ps=nvt_time,
            output_prefix=output_prefix,
            record_interval_ps=record_interval,
            whole_traj_file=whole_traj_file,
        )

    def run_npt(
        self,
        npt_time: float = 100,
        output_prefix: str = "npt",
        record_interval: float = 1,
        whole_traj_file: Optional[str] = None,
    ) -> None:
        """Run an NPT (constant pressure) MD segment.

        Args:
            npt_time (float, optional): Simulation time (ps).
            output_prefix (str, optional): Output file prefix.
            record_interval (float, optional): Frame/log interval (ps).
            whole_traj_file (Optional[str], optional): Append all frames to this PDB if provided.
        """
        self._run_md(
            mode="npt",
            time_ps=npt_time,
            output_prefix=output_prefix,
            record_interval_ps=record_interval,
            whole_traj_file=whole_traj_file,
        )

    def run_pipeline(
        self,
        steps: Optional[List[str]] = None,
        nvt_time: float = 100,
        npt_time: float = 100,
        output_prefix: str = "mdsim",
        record_interval: float = 1,
        whole_traj: bool = False,
    ) -> None:
        """Executes an MD workflow: minimization, NVT, and/or NPT on demand.

        The method prepares the system, then runs the requested steps in order.
        If whole_traj is enabled, all coordinates are appended to a single trajectory PDB.

        Args:
            steps (Optional[List[str]], optional): List of stages to run (e.g., ['em', 'nvt', 'npt']).
            nvt_time (float, optional): Time for NVT segment (ps).
            npt_time (float, optional): Time for NPT segment (ps).
            output_prefix (str, optional): Prefix for all outputs/logs.
            record_interval (float, optional): Frame/log interval (ps).
            whole_traj (bool, optional): If True, write all frames to a single PDB.

        Note:
            The system and OpenMM integrator are initialized in this function.
        """
        if steps is None:
            steps = ["em", "nvt", "npt"]
        self.temperature = self.temperature_K * kelvin
        self.timestep = self.timestep_fs * femtoseconds
        whole_traj_file = f"{output_prefix}_whole.pdb" if whole_traj else None

        self.system = self.top.createSystem(
            nonbondedMethod=PME if self.system_pbc else NoCutoff,
            nonbondedCutoff=1 * nanometer,
            constraints=HBonds,
            rigidWater=True,
        )

        for step in steps:
            if step.lower() == "em":
                self.run_em(output_prefix=f"{output_prefix}_em", whole_traj_file=whole_traj_file)
            elif step.lower() == "nvt":
                self.run_nvt(
                    nvt_time=nvt_time,
                    output_prefix=f"{output_prefix}_nvt",
                    record_interval=record_interval,
                    whole_traj_file=whole_traj_file,
                )
            elif step.lower() == "npt":
                self.run_npt(
                    npt_time=npt_time,
                    output_prefix=f"{output_prefix}_npt",
                    record_interval=record_interval,
                    whole_traj_file=whole_traj_file,
                )
            else:
                print(f"Unknown step: {step}")
