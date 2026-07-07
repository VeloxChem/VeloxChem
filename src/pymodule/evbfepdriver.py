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
import gc

from .veloxchemlib import mpi_master
from .outputstream import OutputStream
from .evbreporter import (EvbReporter, EvbReporterClient, EvbReporterServer,
                          EvbGpuRecalculator)
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

        # Recalculation mode selects how per-frame EVB energies are produced:
        #   "auto"     - offload to a reporter worker if >1 rank, else in-process
        #   "inline"   - always in-process on the CPU (blocks GPU sampling)
        #   "offload"  - always use the 2-rank async CPU reporter worker
        #   "deferred" - sample without recalculation, then re-score each window's
        #                frames in one batched GPU pass read back from the trajectory
        self.recalc_mode: str = "auto"

        self.debug: bool = False
        self.save_frames: int = 1000
        self.save_crash_pdb: bool = True
        self.save_crash_xml: bool = True
        self.save_equil_traj: bool = False
        self.save_equil_pdb: bool = False
        self.xml_crash_save_interval: int = 50
        self.pdb_crash_save_interval: int = 1
        self.pdb_equil_start_temp = 10  # kelvin
        self.pdb_equil_temp_step = 50  # kelvin
        self.pdb_temperatures = []
        self.isobaric: bool = False

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
            "recalc_mode": {
                "type": str
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
            "isobaric": {
                "type": bool
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
        cwd = Path.cwd()
        self.run_folder = cwd / configuration["run_folder"]
        self.data_folder = cwd / configuration["data_folder"]
        self.Lambda = Lambda
        # systems / topology are present on the master (built in-process) and on
        # the async reporter worker (loaded from disk in EvbDriver.run_FEP).
        # Ranks that only synchronise never touch them. initial_positions is
        # read later, on the master path only.
        self.systems = configuration.get("systems")
        self.topology = configuration.get("topology")

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
        n_lambda = len(self.Lambda)
        snapshots_per_lambda = self.sample_steps / self.write_step
        self.total_snapshots = snapshots_per_lambda * n_lambda * n_passes

        # Simulated system time (ps) at various granularities.
        time_per_snapshot = self.step_size * self.write_step
        time_per_lambda = self.step_size * self.sample_steps
        time_per_sweep = time_per_lambda * n_lambda
        total_time = time_per_sweep * n_passes

        self.ostream.print_info(f"Lambda: {np.array(self.Lambda)}")
        info = "\n".join([
            f"Replicas: {self.n_replicas} (forward + backward each)",
            f"Total lambda sweeps: {n_passes}",
            f"Total lambda points: {n_lambda}",
            f"NVT equilibration steps: {self.equil_NVT_steps}",
            f"NPT equilibration steps: {self.equil_NPT_steps}",
            f"Total sample steps: {self.sample_steps}",
            f"Write step: {self.write_step}",
            f"Step size: {self.step_size} ps",
            f"Snapshots per lambda-frame: {snapshots_per_lambda}",
            f"Snapshots to be recorded: {self.total_snapshots}",
            f"System time per snapshot: {time_per_snapshot} ps",
            f"System time per lambda-frame: {time_per_lambda} ps",
            f"System time per lambda sweep: {time_per_sweep} ps",
            f"total system time: {total_time} ps",
        ])

        self.ostream.print_info(info)
        self.ostream.flush()

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

        # Resolve the recalculation mode. "deferred" is orthogonal to rank count:
        # it samples each window without recalculation and then re-scores that
        # window's frames in one batched GPU pass, so it forces the synchronous
        # single-master control flow (no reporter worker).
        assert_msg_critical(
            self.recalc_mode in ("auto", "inline", "offload", "deferred"),
            f"Unknown recalc_mode '{self.recalc_mode}'. Expected one of "
            "'auto', 'inline', 'offload', 'deferred'.")
        self._deferred = self.recalc_mode == "deferred"
        if self.recalc_mode == "offload" and self.nodes <= 1:
            self.ostream.print_warning(
                "recalc_mode='offload' requires >1 MPI rank; falling back to "
                "in-process synchronous reporting.")

        # Asynchronous reporting: with >=2 MPI ranks on the node, rank
        # (master+1) becomes a dedicated CPU reporter worker that owns the
        # Energies/Forces/Velocities files, while the master drives the GPU and
        # offloads each frame to it. With a single rank the reporting stays
        # synchronous and in-process (unchanged behaviour).
        self._async = (self.nodes > 1 and not self._deferred
                       and self.recalc_mode in ("auto", "offload"))

        if self._deferred:
            # Guard combinations that a positions-only trajectory cannot serve.
            assert_msg_critical(
                not (self.report_velocities or self.debug),
                "recalc_mode='deferred' cannot report velocities: they are not "
                "stored in the XTC trajectory. Use recalc_mode 'inline'/'offload'."
            )
            assert_msg_critical(
                not (self.report_forcegroups or self.debug),
                "recalc_mode='deferred' does not support force-group / bonded "
                "decomposition. Use recalc_mode 'inline'/'offload'.")
            if self.rank == mpi_master():
                self.ostream.print_info(
                    "Running in deferred mode: sampling without recalculation, "
                    "then re-scoring each window on the GPU from the trajectory."
                )
                if self.nodes > 1:
                    self.ostream.print_warning(
                        f"{self.nodes} MPI ranks present but recalc_mode="
                        "'deferred' uses only the master; the rest will idle.")
        elif self._async and self.rank == mpi_master():
            self.ostream.print_info(
                f"Running in asynchronous mode with {self.nodes} MPI ranks.")
            if self.nodes > 2:
                self.ostream.print_warning(
                    "More than 2 MPI ranks are present. Only the master and one reporter worker are used; the rest will idle."
                )
        if not self._async and not self._deferred:
            self.ostream.print_info(
                "Running in synchronous mode (single rank or no MPI).")

        self._reporter_worker = mpi_master() + 1
        # Allocate the shared-memory ring buffer (collective over COMM_WORLD via
        # the Split below) before the master and reporter worker diverge. The
        # bulk per-frame payload travels through this window; only tiny control
        # messages stay on the point-to-point channel.
        if self._async:
            self._setup_shared_ring()
        if self._async and self.rank != mpi_master():
            if self.rank == self._reporter_worker:
                self._serve_reporter()
                self._free_shared_ring()
            # The master signals completion with a matching Barrier.
            self.comm.Barrier()
            return
        # Deferred mode runs entirely on the master; any extra ranks idle here
        # and rejoin at the closing Barrier below.
        if self._deferred and self.rank != mpi_master():
            self.comm.Barrier()
            return

        # Master path only: the initial positions are needed to seed the GPU
        # sampling (the reporter worker returned above and never uses them).
        initial_positions = configuration["initial_positions"]

        # Global timer spanning all replicas / sweeps for total-progress ETA.
        self.timer = Timer(self.total_steps)

        # Single trajectory file, appended across all sweeps. Frames stay
        # row-aligned with Energies.csv (each row tagged with replica/direction).
        # In deferred mode the trajectory reporter is recreated per window (so
        # the file is flushed before it is read back), so skip the shared one.
        self.traj_roporter = None
        if not self._deferred:
            self.traj_roporter = mmapp.XTCReporter(
                str(self.data_folder / "trajectory.xtc"),
                self.write_step,
            )
        # Only the very first written window truncates the output files and
        # writes their headers; everything afterwards appends.
        self._first_write = True

        # Deferred mode: a single recalculator owns the Energies/Forces/NB
        # outputs across all windows and re-scores each window on the GPU after
        # its sampling context has been released.
        self._gpu_recalc = None
        if self._deferred:
            self._gpu_recalc = EvbGpuRecalculator(
                self.data_folder,
                self.systems,
                self.topology,
                lambda system: self._get_simulation(system, self.step_size),
                report_forces=self.report_forces,
            )
        self.timer.start()

        positions = initial_positions * 0.1
        velocities = None

        # Initial equilibration, performed once on the l == 0 system. The
        # resulting positions / velocities seed the first forward sweep.
        if not self.skip_initial_equil and 0 in self.Lambda:
            system = self.systems[0]
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

        # Crude profiling accumulators over all sampled lambda windows. Each
        # sampling window's wall time is split into pure GPU sampling and the
        # master->reporter communication that happens inside the sampling loop,
        # so that gpu + communication == window (checked below).
        self._t_window_total = 0.0
        self._t_gpu_total = 0.0
        self._t_comm_total = 0.0
        self._t_stall_total = 0.0
        self._t_send_total = 0.0
        self._n_windows = 0
        # Deferred mode: separate accumulator for the post-sampling GPU
        # recalculation (runs after each window, not overlapped with sampling).
        self._t_recalc_total = 0.0

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

        # Tear down the reporter worker and synchronise all ranks.
        reporter_timing = None
        if self._async:
            EvbReporterClient.send_terminate(self.comm, self._reporter_worker)
            # Collect the reporter-side timing sent back after TERMINATE.
            reporter_timing = np.empty(5, dtype=np.float64)
            self.comm.Recv([reporter_timing, MPI.DOUBLE],
                           source=self._reporter_worker,
                           tag=EvbReporterClient._TAG)
            self._free_shared_ring()
            self.comm.Barrier()
        elif self._deferred and self.nodes > 1:
            # Release the idle ranks that returned early in deferred mode.
            self.comm.Barrier()

        if self._gpu_recalc is not None:
            self._gpu_recalc.close()

        self._print_timing_summary(reporter_timing)
        self.ostream.flush()

    def _setup_shared_ring(self, queue_depth=3):
        """Allocate the master<->reporter shared-memory ring buffer.

        Collective over COMM_WORLD (all ranks call the Split). Only the master
        and the reporter worker actually map the window; any extra idle ranks
        get COMM_NULL and skip it. The two participants are assumed to be on the
        same node (single-node async reporting), so MPI.Win.Allocate_shared maps
        one physical buffer into both processes.
        """
        self._shm_win = None
        self._shm_comm = None
        self._pair_comm = None
        color = (0 if self.rank in (mpi_master(),
                                    self._reporter_worker) else MPI.UNDEFINED)
        pair_comm = self.comm.Split(color, self.rank)
        if pair_comm == MPI.COMM_NULL:
            return
        shm_comm = pair_comm.Split_type(MPI.COMM_TYPE_SHARED)
        # Map each shm member's world rank to its shm-local rank so we can query
        # the master's (the producer's) segment.
        world_ranks = shm_comm.allgather(self.rank)
        assert_msg_critical(
            mpi_master() in world_ranks,
            "EVB async reporter requires the master and reporter worker to share "
            "a node for shared-memory reporting.")
        producer_shm = world_ranks.index(mpi_master())

        num_atoms = self.topology.getNumAtoms()
        # header(3) + box(9) + positions/velocities/forces (9N), fixed layout.
        max_len = 3 + 9 + 9 * num_atoms
        itemsize = MPI.DOUBLE.Get_size()
        nbytes = (queue_depth * max_len *
                  itemsize) if self.rank == mpi_master() else 0
        win = MPI.Win.Allocate_shared(nbytes, itemsize, comm=shm_comm)
        buf, _ = win.Shared_query(producer_shm)
        ring = np.frombuffer(buf,
                             dtype=np.float64).reshape(queue_depth, max_len)
        win.Lock_all()

        self._pair_comm = pair_comm
        self._shm_comm = shm_comm
        self._shm_win = win
        self._shm_ring = ring
        self._queue_depth = queue_depth

    def _free_shared_ring(self):
        """Tear down the shared-memory window (participants only)."""
        if getattr(self, "_shm_win", None) is None:
            return
        self._shm_win.Unlock_all()
        self._shm_win.Free()
        self._shm_comm.Free()
        self._pair_comm.Free()
        self._shm_win = None

    def _serve_reporter(self):
        """Reporter-worker entry point (non-master rank, async mode).

        Builds the CPU energy sims once and serves frames from the master until
        TERMINATE. Mirrors the file selection of _get_evb_reporter so the worker
        writes exactly the files the synchronous path would.
        """
        force_file = (str(self.data_folder / "Forces.csv")
                      if self.report_forces or self.debug else None)
        velocity_file = (str(self.data_folder / "Velocities.csv")
                         if self.report_velocities or self.debug else None)

        server = EvbReporterServer(
            self.comm,
            mpi_master(),
            str(self.data_folder / "Energies.csv"),
            self.write_step,
            self.systems,
            self.topology,
            self.ostream,
            self._shm_ring,
            self._shm_win,
            forming_bonds=self.forming_bonds,
            breaking_bonds=self.breaking_bonds,
            force_file=force_file,
            velocity_file=velocity_file,
        )
        server.serve()

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

            equil_state = self._equilibrate(system,
                                            f"r{replica}_d{direction}_l{l:.3}",
                                            positions, velocities)

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
            barostat = self._get_barostat(simulation)
            if barostat:
                barostat.setFrequency(0)
                self._safe_step(simulation, self.initial_equil_NVT_steps,
                                "initial NVT equilibration")
                barostat.setFrequency(25)
                self._safe_step(simulation, self.initial_equil_NPT_steps,
                                "initial NPT equilibration")
                barostat.setFrequency(0)
            else:
                self._safe_step(simulation, self.initial_equil_NVT_steps,
                                "initial NVT equilibration")

        else:
            self.ostream.print_info(
                f"Perfoming PDB warmup with T-vector {np.array(self.pdb_temperatures)}"
            )
            self.ostream.flush()
            for T in self.pdb_temperatures:
                simulation.integrator.setTemperature(T)
                if T == self.temperature:
                    NVT_steps = self.initial_equil_NVT_steps
                else:
                    if self.warmup_NVT_steps > 0:
                        NVT_steps = self.warmup_NVT_steps
                    else:
                        NVT_steps = self.initial_equil_NVT_steps

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

            barostat = self._get_barostat(simulation)
            if barostat:
                barostat.setFrequency(25)
                self._safe_step(
                    simulation, self.initial_equil_NPT_steps,
                    f"NPT equilibration at target temperature T = {self.temperature}"
                )
                barostat.setFrequency(0)

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

    def _equilibrate(self, system, name, positions, velocities=None):
        simulation = self._get_simulation(system, self.equil_step_size)
        simulation.context.setPositions(positions)
        if velocities is not None:
            simulation.context.setVelocities(velocities)
        if self.minimize_every_lambda:
            self._minimize(simulation)

        if self.debug:
            self._save_state(simulation, f"pre_equil_{name}")
            simulation.reporters.append(
                mmapp.PDBReporter(
                    str(self.run_folder / f"traj_equil_{name}.pdb"),
                    self.write_step,
                    enforcePeriodicBox=True,
                ))
            f_file = str(self.run_folder / f"forces_equil_{name}.csv")
            v_file = str(self.run_folder / f"velocities_equil_{name}.csv")
            g_file = str(self.run_folder / f"forcegroups_equil_{name}.csv")
            simulation.reporters.append(
                EvbReporter(
                    str(self.run_folder / f"energies_equil_{name}.csv"),
                    self.write_step,
                    self.systems,
                    self.topology,
                    name,
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
                str(self.run_folder / f"equil_traj_{name}.xtc"),
                self.write_step,
                enforcePeriodicBox=True,
            )
            simulation.reporters.append(equil_traj_reporter)

        if self.isobaric:
            barostat = self._get_barostat(simulation)

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
                f"equil_state_{name}",
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
        barostats = [
            force for force in simulation.system.getForces()
            if isinstance(force, mm.MonteCarloBarostat)
            or isinstance(force, mm.MonteCarloFlexibleBarostat)
        ]
        if len(barostats) == 0:
            return None
        else:
            return barostats[0]

    def _sample(self, system, l, initial_state, replica=0, direction=0):
        if self._deferred:
            return self._sample_deferred(system, l, initial_state, replica,
                                         direction)
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

        evb_reporters = self._make_evb_reporters(l, append, replica, direction)
        simulation.reporters.append(state_reporter)
        for reporter in evb_reporters:
            simulation.reporters.append(reporter)
        self.ostream.flush()

        # Snapshot the client's communication clock before sampling; the BEGIN
        # control message it already sent should not count against this window's
        # in-loop communication.
        client = next(
            (r for r in evb_reporters if isinstance(r, EvbReporterClient)),
            None)
        comm_before = client.comm_time if client is not None else 0.0
        stall_before = client.stall_time if client is not None else 0.0
        send_before = client.send_time if client is not None else 0.0

        # Time the sampling window. The communication done inside the sampling
        # loop (client.report -> MPI send) is a subset of this wall time, so
        # gpu_time = window - comm_during and gpu + comm == window exactly.
        t_window_start = time.perf_counter()
        states = self._safe_step(simulation, self.sample_steps, "sampling")
        t_window = time.perf_counter() - t_window_start

        # Break communication into: stall (blocked on the reporter freeing a
        # slot), transfer (shared-memory fence + control Send), and the derived
        # remainder extract+pack (pulling the state and filling the ring).
        comm_during = (client.comm_time -
                       comm_before) if client is not None else 0.0
        stall_during = (client.stall_time -
                        stall_before) if client is not None else 0.0
        send_during = (client.send_time -
                       send_before) if client is not None else 0.0

        # Close out the async window (wait for in-flight sends, signal END).
        for reporter in evb_reporters:
            if isinstance(reporter, EvbReporterClient):
                reporter.finalize()

        self._accumulate_window_timing(l, replica, direction, t_window,
                                       comm_during, stall_during, send_during)
        return states[-1]

    def _sample_deferred(self,
                         system,
                         l,
                         initial_state,
                         replica=0,
                         direction=0):
        """Sample one window without recalculation, then re-score its frames on
        the GPU from the trajectory once the sampling context is released.

        Only the trajectory and StateData reporters are attached during
        sampling (no EVB energy evaluation), so the GPU is never stalled. After
        sampling, the sampling context is torn down and the window's frames are
        read back from ``trajectory.xtc`` and re-scored in one batched GPU pass.
        """
        simulation = self._get_simulation(system, self.equil_step_size)

        append = not self._first_write
        self._first_write = False

        # Per-window trajectory reporter appended into the single trajectory.xtc.
        # Recreated (and dropped) each window so the file is flushed to disk
        # before the recalculation reads it back. The first written window
        # truncates the file; the rest append.
        traj_reporter = mmapp.XTCReporter(
            str(self.data_folder / "trajectory.xtc"),
            self.write_step,
            append=append,
        )
        simulation.reporters.append(traj_reporter)

        sz = self.step_size * mmunit.picoseconds
        simulation.integrator.setStepSize(sz)
        simulation.context.setState(initial_state)

        state_reporter = self._get_data_reporter(
            self.data_folder,
            "Data_combined.csv",
            append,
        )
        simulation.reporters.append(state_reporter)
        self.ostream.flush()

        t_window_start = time.perf_counter()
        states = self._safe_step(simulation, self.sample_steps, "sampling")
        t_window = time.perf_counter() - t_window_start

        # Detached snapshot to seed the next window; survives context teardown.
        final_state = states[-1]

        # Release the sampling GPU context (and flush/close the trajectory) so
        # the recalculation can build its own context: only one OpenMM context
        # should live on the GPU at a time.
        simulation.reporters.clear()
        del simulation, states, state_reporter, traj_reporter
        gc.collect()

        t_recalc_start = time.perf_counter()
        self._recalc_window_gpu(l, replica, direction, append)
        t_recalc = time.perf_counter() - t_recalc_start

        self._accumulate_deferred_timing(l, replica, direction, t_window,
                                         t_recalc)
        return final_state

    def _recalc_window_gpu(self, l, replica, direction, append):
        """Read this window's frames back from the trajectory and re-score them
        on the GPU via the batched EvbGpuRecalculator."""
        import MDAnalysis as mda

        snapshots = int(self.sample_steps // self.write_step)
        # Use the in-memory OpenMM topology directly (no dependency on a parsed
        # topology file); the XTC supplies the per-frame coordinates.
        universe = mda.Universe(
            self.topology,
            str(self.data_folder / "trajectory.xtc"),
            topology_format='OPENMMTOPOLOGY',
        )
        # Positions/box arrive from MDAnalysis in angstrom; _apply_positions
        # expects nanometers.
        frames = []
        for ts in universe.trajectory[-snapshots:]:
            pos_nm = np.array(ts.positions, dtype=np.float64) * 0.1
            # None for a non-periodic (vacuum) trajectory; a (3,3) box for a
            # periodic (solvated) one.
            box = ts.triclinic_dimensions
            box_nm = None if box is None else np.array(box,
                                                       dtype=np.float64) * 0.1
            frames.append((pos_nm, box_nm))

        assert_msg_critical(
            len(frames) == snapshots,
            f"Deferred recalculation expected {snapshots} trajectory frames for "
            f"window l={l} but read {len(frames)}.")

        self._gpu_recalc.recalc_window(l, replica, direction, frames, append)

    def _accumulate_deferred_timing(self, l, replica, direction, t_window,
                                    t_recalc):
        """Record and print the sampling vs recalculation split for one window
        in deferred mode."""
        self._t_window_total += t_window
        self._t_gpu_total += t_window
        self._t_recalc_total += t_recalc
        self._n_windows += 1
        self.ostream.print_info(
            f"[timing] window l={l} rep={replica} dir={direction} (deferred): "
            f"sampling(GPU) {t_window:.3f}s + recalculation(GPU) "
            f"{t_recalc:.3f}s = {t_window + t_recalc:.3f}s")
        self.ostream.flush()

    def _accumulate_window_timing(self, l, replica, direction, t_window,
                                  comm_during, stall_during, send_during):
        """Record and print the timing breakdown for one sampling window."""
        gpu_time = t_window - comm_during
        extract_during = comm_during - stall_during - send_during
        self._t_window_total += t_window
        self._t_gpu_total += gpu_time
        self._t_comm_total += comm_during
        self._t_stall_total += stall_during
        self._t_send_total += send_during
        self._n_windows += 1
        mode = "async" if self._async else "sync"
        self.ostream.print_info(
            f"[timing] window l={l} rep={replica} dir={direction} ({mode}): "
            f"total={t_window:.3f}s = sampling(GPU) {gpu_time:.3f}s "
            f"+ communication {comm_during:.3f}s "
            f"(stall {stall_during:.3f}s + transfer {send_during:.3f}s "
            f"+ extract {extract_during:.3f}s) "
            f"(check gpu+comm={gpu_time + comm_during:.3f}s)")
        self.ostream.flush()

    def _print_timing_summary(self, reporter_timing):
        """Print the aggregate profiling summary and consistency checks."""
        n = self._n_windows
        if n == 0:
            return
        if self._deferred:
            recalc = self._t_recalc_total
            sampling = self._t_gpu_total
            total = sampling + recalc
            per_win = recalc / n if n > 0 else 0.0
            self.ostream.print_info(
                f"[timing] Totals over {n} windows (deferred): "
                f"sampling(GPU)={sampling:.3f}s + recalculation(GPU)="
                f"{recalc:.3f}s = {total:.3f}s "
                f"({per_win:.3f}s recalc/window). Sampling ran without any "
                "in-loop recalculation; the GPU re-scored each window in a "
                "single batched pass afterwards.")
            self.ostream.flush()
            return
        summed = self._t_gpu_total + self._t_comm_total
        extract_total = (self._t_comm_total - self._t_stall_total -
                         self._t_send_total)
        self.ostream.print_info(f"[timing] Totals over {n} sampling windows: "
                                f"window={self._t_window_total:.3f}s, "
                                f"sampling(GPU)={self._t_gpu_total:.3f}s, "
                                f"communication={self._t_comm_total:.3f}s")
        self.ostream.print_info(
            f"[timing] Communication breakdown: "
            f"stall(reporter backpressure)={self._t_stall_total:.3f}s + "
            f"transfer(fence+send)={self._t_send_total:.3f}s + "
            f"extract+pack={extract_total:.3f}s")
        self.ostream.print_info(
            f"[timing] Consistency: sampling+communication={summed:.3f}s vs "
            f"summed window total={self._t_window_total:.3f}s "
            f"(diff {summed - self._t_window_total:+.3e}s)")
        if reporter_timing is not None:
            (energy_time, recv_time, n_frames, n_sims,
             energy_time_max) = reporter_timing
            per_frame = energy_time / n_frames if n_frames > 0 else 0.0
            self.ostream.print_info(
                f"[timing] Reporter rank: energy calculation={energy_time:.3f}s "
                f"over {int(n_frames)} frames "
                f"({per_frame * 1e3:.1f} ms/frame avg, "
                f"{energy_time_max * 1e3:.1f} ms/frame max, "
                f"{int(n_sims)} CPU sims/frame), "
                f"idle in Recv={recv_time:.3f}s")
            self.ostream.print_info(
                "[timing] MPI benefit: the "
                f"{energy_time:.3f}s of energy evaluation ran in parallel with "
                f"the {self._t_gpu_total:.3f}s of GPU sampling at a cost of "
                f"{self._t_comm_total:.3f}s communication; synchronous "
                "reporting would have added the energy time to every window.")
            # Verdict: what actually gates the pipeline. Stall dominating the
            # communication total means the reporter (CPU energy eval) cannot
            # keep up with GPU sampling; transfer dominating would mean the data
            # movement itself is the cost.
            if self._t_comm_total > 0:
                if self._t_stall_total >= self._t_send_total + extract_total:
                    verdict = (
                        "reporter-bound: the CPU energy evaluation is the "
                        "bottleneck (backpressure stall dominates); "
                        "communication/transfer is not the problem")
                else:
                    verdict = ("transfer-bound: data movement dominates the "
                               "communication cost")
                self.ostream.print_info(f"[timing] Verdict: {verdict}.")
        if self.report_forcegroups or self.debug:
            self.ostream.print_warning(
                "[timing] A master-side synchronous reporter is active "
                "(report_forcegroups/debug): its CPU energy work runs inside the "
                "sampling step and inflates the reported sampling(GPU) figure.")
        self.ostream.flush()

    def _make_evb_reporters(self, l, append, replica, direction):
        """Build the EVB reporters for one sampling window.

        Synchronous mode returns a single in-process EvbReporter. Async mode
        returns an EvbReporterClient that offloads the core energies to the
        worker, plus (only when force-groups / bonded decomposition are
        requested) a master-side EvbReporter that writes just those files.
        """
        if not self._async:
            return [
                self._get_evb_reporter(self.data_folder,
                                       l,
                                       append=append,
                                       replica=replica,
                                       direction=direction)
            ]

        client = EvbReporterClient(
            self.comm,
            self._reporter_worker,
            self.write_step,
            self.topology.getNumAtoms(),
            l,
            replica,
            direction,
            append,
            self.report_forces or self.debug,
            self.report_velocities or self.debug,
            self._shm_ring,
            self._shm_win,
        )
        reporters = [client]

        if self.report_forcegroups or self.debug:
            fg_reporter = EvbReporter(
                str(self.data_folder / "Energies.csv"),
                self.write_step,
                self.systems,
                self.topology,
                l,
                self.ostream,
                forming_bonds=self.forming_bonds,
                breaking_bonds=self.breaking_bonds,
                forcegroup_file=str(self.data_folder / "ForceGroups.csv"),
                append=append,
                replica=replica,
                direction=direction,
                core_outputs=False,
            )
            reporters.append(fg_reporter)
        return reporters

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
        self.ostream.print_info(
            f"Running {name} for {steps} steps for total time: {steps*self.step_size} ps with step size {self.step_size} ps"
        )

        self.ostream.flush()
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
