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
"""OpenMM MD and solvation FEP for aqueous CDAM systems.

``md_settings`` resolves the MD protocol (production default, or the
``pilot`` and ``fast_mode`` presets) and ``run_md`` minimizes, heats (NVT),
equilibrates the density (NPT) and records production frames in memory.

``run_fep`` decouples the ion from water in two stages (charge, then LJ) with
VeloxChem's ``SolvationFepDriver``. The explicit 12-6-4 C4 term follows the
charge schedule, and the vacuum legs are skipped for a rigid CDAM cluster.
``finite_size_correction`` adds the analytical periodic-box and summation
corrections to the raw FEP result.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from copy import deepcopy
from typing import Any
import collections
import logging
import math
import os
import re
import tempfile
from pathlib import Path

import numpy as np
try:
    import openmm as mm
    import openmm.app as app
    import openmm.unit as unit
except ImportError:
    pass

from ..solvationfepdriver import SolvationFepDriver
from ..veloxchemlib import bohr_in_angstrom, hartree_in_kjpermol
from .candidate import CDAMCandidate
from .generator import C4_FORCE_NAME
from .models import (ANGSTROM_TO_NM, PROGRESS_STATUS, default_ostream,
                     print_settings)

# ==========================================================================
# md
# ==========================================================================

TEMPERATURE_K = 298.15
PRESSURE_ATM = 1.0
FRICTION_PER_PS = 1.0
BAROSTAT_INTERVAL = 25
HEATING_STAGES = 10  # NVT heating ramps 0 K -> T in this many steps

# Production sampling is the default.
DEFAULT_MD_OPTIONS = {
    "timestep_fs": 2.0,
    "minimization_iterations": 1000,
    "nvt_steps": 5000,
    "npt_steps": 5000,
    "production_steps": 250000,
    "report_interval": 250,
    "random_seed": 0,
}
# pilot=True: 0.05 ps production, a workflow test only; values are not meaningful.
PILOT_MD_OPTIONS = {
    "platform": "CPU",
    "platform_properties": {
        "Threads": "2"
    },
    "timestep_fs": 0.25,
    "nvt_steps": 100,
    "npt_steps": 100,
    "production_steps": 200,
    "report_interval": 20,
    "random_seed": 42,
}
# fast_mode=True: short but usable MD (10 ps NVT, 10 ps NPT, 50 ps production)
# for screening; check the final parameters with the production protocol.
FAST_MD_OPTIONS = {
    "production_steps": 25000,
    "random_seed": 42,
}
MD_OPTION_NAMES = frozenset({
    "nvt_steps",
    "npt_steps",
    "production_steps",
    "report_interval",
    "minimization_iterations",
    "timestep_fs",
    "random_seed",
    "platform",
    "platform_properties",
    "pilot",
    "fast_mode",
})


def md_settings(options: dict[str, Any] | None = None) -> dict[str, Any]:
    """Production defaults, then the requested preset, then explicit overrides."""
    options = dict(options or {})
    settings = dict(DEFAULT_MD_OPTIONS)
    if options.get("pilot") and options.get("fast_mode"):
        raise ValueError("choose pilot or fast_mode, not both")
    if options.get("pilot"):
        settings.update(PILOT_MD_OPTIONS)
    if options.get("fast_mode"):
        settings.update(FAST_MD_OPTIONS)
    settings.update(options)
    if "platform_properties" in settings:
        settings["platform_properties"] = dict(settings["platform_properties"])
    return settings


def check_md_settings(settings: dict[str, Any]) -> None:
    """Validate resolved MD settings; raises ``ValueError`` on unknown names or bad step counts.

    Production must be long enough to record at least two frames.
    """
    unknown = set(settings) - MD_OPTION_NAMES
    if unknown:
        raise ValueError(f"unknown MD options: {sorted(unknown)}")
    timestep = float(settings["timestep_fs"])
    if not np.isfinite(timestep) or timestep <= 0:
        raise ValueError("timestep_fs must be positive and finite")
    for name in ("nvt_steps", "npt_steps", "minimization_iterations"):
        if not isinstance(settings[name], int) or settings[name] < 0:
            raise ValueError(f"{name} must be a nonnegative integer")
    for name in ("production_steps", "report_interval"):
        if not isinstance(settings[name], int) or settings[name] < 1:
            raise ValueError(f"{name} must be a positive integer")
    if settings["production_steps"] < 2 * settings["report_interval"]:
        raise ValueError(
            "production_steps must record at least two trajectory frames")


# OpenMM's CPU platform uses every core by default. For a CDAM water box
# (~1300 atoms) PME scales only to ~8 threads, and 128 threads run 2x slower
# than 8. OPENMM_CPU_THREADS overrides this default.
DEFAULT_CPU_THREADS = 8


def cpu_threads() -> str:
    """Number of OpenMM CPU threads: OPENMM_CPU_THREADS, else up to 8."""
    return os.environ.get("OPENMM_CPU_THREADS") or str(
        min(DEFAULT_CPU_THREADS, os.cpu_count() or 1))


@lru_cache(maxsize=None)
def default_platform() -> tuple[mm.Platform, tuple]:
    """Fastest available OpenMM platform and its default properties (as items).

    Also sets the default thread count of the CPU platform, so that contexts
    created elsewhere (the solvation FEP driver) do not use every core.
    """
    mm.Platform.getPlatformByName("CPU").setPropertyDefaultValue(
        "Threads", cpu_threads())
    for name in ("HIP", "CUDA", "OpenCL", "CPU"):
        try:
            platform = mm.Platform.getPlatformByName(name)
        except Exception:
            continue
        if name in ("CUDA", "HIP"):
            return platform, (("Precision", "mixed"),)
        if name == "CPU":
            return platform, (("Threads", cpu_threads()),)
        return platform, ()
    raise RuntimeError("no OpenMM platform is available")


def platform_for(
        settings: dict[str, Any]) -> tuple[mm.Platform, dict[str, str]]:
    """The platform named in ``settings``, else the default one."""
    if "platform" in settings:
        properties = dict(settings.get("platform_properties", {}))
        if settings["platform"] == "CPU":
            properties.setdefault("Threads", cpu_threads())
        return mm.Platform.getPlatformByName(settings["platform"]), properties
    platform, properties = default_platform()
    return platform, dict(settings.get("platform_properties", dict(properties)))


class MemoryReporter:
    """OpenMM reporter that keeps positions and box vectors (nm) in memory."""

    def __init__(self, interval: int):
        """Record a frame every ``interval`` MD steps."""
        self._interval = interval
        self.frames, self.boxes = [], []

    def describeNextReport(self, simulation):
        """Ask OpenMM for positions only, at the next multiple of the interval."""
        steps = self._interval - simulation.currentStep % self._interval
        return (steps, True, False, False, False, None)

    def report(self, simulation, state):
        """Append this frame and its box vectors (nm)."""
        self.frames.append(
            state.getPositions(asNumpy=True).value_in_unit(unit.nanometers))
        self.boxes.append(
            state.getPeriodicBoxVectors(asNumpy=True).value_in_unit(
                unit.nanometers))


@dataclass
class Trajectory:
    """In-memory production MD: final state plus recorded frames and box vectors."""

    final_state: mm.State
    frames: np.ndarray  # (n_frames, n_atoms, 3) nm
    boxes: np.ndarray  # (n_frames, 3, 3) nm, box vectors as rows


def _add_barostat(system: mm.System, simulation: app.Simulation,
                  seed: int) -> None:
    """Add a 1 atm, 298.15 K Monte Carlo barostat and reinitialize the context, keeping state.

    The barostat seed is offset from the integrator seed; 0 lets OpenMM choose.
    """
    barostat = mm.MonteCarloBarostat(PRESSURE_ATM * unit.atmospheres,
                                     TEMPERATURE_K * unit.kelvin,
                                     BAROSTAT_INTERVAL)
    barostat.setRandomNumberSeed(seed + 1 if seed else 0)
    system.addForce(barostat)
    simulation.context.reinitialize(preserveState=True)


def run_md(system: mm.System,
           topology,
           positions,
           options: dict[str, Any] | None = None,
           *,
           ostream=None) -> Trajectory | None:
    """Minimize, heat (NVT), equilibrate (NPT), then record production frames.

    Zero NVT/NPT steps skip that stage; zero minimization iterations minimize
    to OpenMM's tolerance. ``system`` is modified (barostat), so pass a copy.
    Returns None if OpenMM fails (e.g. the system blows up).
    """
    ostream = default_ostream(ostream)
    settings = md_settings(options)
    check_md_settings(settings)
    platform, properties = platform_for(settings)
    seed = settings["random_seed"]
    nvt, npt, production = settings["nvt_steps"], settings[
        "npt_steps"], settings["production_steps"]
    interval = settings["report_interval"]
    timestep = float(settings["timestep_fs"]) * unit.femtoseconds
    step_ps = timestep.value_in_unit(unit.picoseconds)
    npt_production = bool(npt)

    print_settings(ostream, "OpenMM MD Setup", [
        ("Protocol", "pilot (workflow test only)" if settings.get("pilot")
         else "fast mode (short MD)" if settings.get("fast_mode") else
         "production"),
        ("Platform", f"{platform.getName()} {properties}"
         if properties else platform.getName()),
        ("Number of Atoms", topology.getNumAtoms()),
        ("Temperature", f"{TEMPERATURE_K:.2f} K"),
        ("Pressure (NPT)", f"{PRESSURE_ATM:.2f} atm"),
        ("Friction", f"{FRICTION_PER_PS:.2f} 1/ps"),
        ("Time Step", f"{settings['timestep_fs']:g} fs"),
        ("Minimization Iterations", settings["minimization_iterations"] or
         "until converged"),
        ("NVT Heating", f"{nvt} steps ({nvt * step_ps:g} ps)"),
        ("NPT Equilibration", f"{npt} steps ({npt * step_ps:g} ps)"),
        ("Production", f"{production} steps ({production * step_ps:g} ps), "
         f"{'NPT' if npt_production else 'NVT'} ensemble"),
        ("Report Interval",
         f"{interval} steps ({production // interval} frames)"),
        ("Random Seed", seed if seed else "0 (OpenMM chooses)"),
    ])

    def info(message):
        """Print a progress message and flush the stream."""
        ostream.print_info(message)
        ostream.flush()

    try:
        # A caller-supplied template may carry a barostat: NVT runs without it.
        for index in reversed(range(system.getNumForces())):
            if isinstance(system.getForce(index), mm.MonteCarloBarostat):
                system.removeForce(index)
        integrator = mm.LangevinMiddleIntegrator(
            TEMPERATURE_K * unit.kelvin, FRICTION_PER_PS / unit.picoseconds,
            timestep)
        integrator.setRandomNumberSeed(seed)
        simulation = app.Simulation(topology, system, integrator, platform,
                                    properties)
        simulation.context.setPositions(positions)

        info("Minimizing energy...")
        simulation.minimizeEnergy(
            maxIterations=settings["minimization_iterations"])

        info(
            f"Heating (NVT) from 0K to {TEMPERATURE_K} K for {nvt * step_ps:g} ps..."
        )
        simulation.context.setVelocitiesToTemperature(10.0 * unit.kelvin, seed)
        if nvt:
            stages = min(HEATING_STAGES, nvt)
            chunks = np.array_split(np.arange(nvt), stages)
            for stage, steps in enumerate(chunks, start=1):
                temperature = TEMPERATURE_K * stage / stages
                integrator.setTemperature(temperature * unit.kelvin)
                simulation.step(len(steps))
        else:
            simulation.context.setVelocitiesToTemperature(
                TEMPERATURE_K * unit.kelvin, seed)

        info(f"Density Equilibration (NPT) for {npt * step_ps:g} ps...")
        if npt_production:
            _add_barostat(system, simulation, seed)
        if npt:
            simulation.step(npt)

        reporter = MemoryReporter(interval)
        simulation.reporters.append(reporter)
        info(f"Running production MD for {production * step_ps:g} ps...")
        simulation.step(production)
        return Trajectory(simulation.context.getState(getPositions=True),
                          np.array(reporter.frames), np.array(reporter.boxes))
    except mm.OpenMMException as error:
        ostream.print_warning(f"MD simulation crashed: {error}")
        return None


# ==========================================================================
# hfe
# ==========================================================================

# Production FEP protocol (overridable through fep_options).
NUM_STEPS = 50000
NUM_SNAPSHOTS = 500
TIMESTEP_FS = 2.0
STAGE1_LAMBDAS = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,
                  1.0]  # charge
STAGE2_LAMBDAS = [1.0, 0.8, 0.6, 0.4, 0.2, 0.0]  # LJ
DRIVER_OPTIONS = ("num_steps", "num_snapshots", "num_equil_steps",
                  "num_em_steps", "platform")
# The applied FEP protocol (SolvationFepDriver equilibration defaults included).
FEP_DEFAULTS = {
    "num_steps": NUM_STEPS,
    "num_snapshots": NUM_SNAPSHOTS,
    "num_equil_steps": 5000,
    "num_em_steps": 0,
    "timestep_fs": TIMESTEP_FS,
    "lambdas_stage1": STAGE1_LAMBDAS,
    "lambdas_stage2": STAGE2_LAMBDAS,
    "equilibrate_first": True,
}


# quick_dirty=True (pilot runs only): the production lambda windows with 0.5
# ps of sampling each (2 fs steps) and one MBAR solve per stage; tests that the
# FEP runs, the HFE value is not meaningful.
QUICK_DIRTY_FEP_OPTIONS = {
    "num_steps": 250,
    "num_snapshots": 10,
    "num_equil_steps": 100,
}


def fep_settings(options: dict[str, Any] | None = None) -> dict[str, Any]:
    """The production FEP protocol (or the quick_dirty preset), then explicit
    overrides."""
    options = dict(options or {})
    settings = deepcopy(FEP_DEFAULTS)
    if options.get("quick_dirty"):
        settings.update(deepcopy(QUICK_DIRTY_FEP_OPTIONS))
    settings.update(options)
    return settings

# Finite-size correction constants (OpenMM units: nm, e, kJ/mol).
# 1/(4 pi eps0) in kJ/mol nm/e^2, from VeloxChem atomic-unit conversions.
COULOMB_CONSTANT = hartree_in_kjpermol() * bohr_in_angstrom() * ANGSTROM_TO_NM
WATER_DIELECTRIC = 78.39  # as used throughout VeloxChem (e.g. OpenMMDynamics)
WIGNER_LATTICE_SUM = -2.837297
# Quadrupole-moment trace of each water model, e nm^2 (TIP3P if unknown).
WATER_QUADRUPOLE_TRACE = {
    "tip3p": 0.00764,
    "spce": 0.00848,
    "spc/e": 0.00848,
    "tip4pew": 0.00944,
    "tip4p-ew": 0.00944,
    "tip4p": 0.00930,
}


def check_fep_backend() -> None:
    """Raise ``RuntimeError`` if ``SolvationFepDriver`` cannot run FEP from OpenMM files."""
    if not callable(
            getattr(SolvationFepDriver, "_scale_coulomb_global_parameters",
                    None)):
        raise RuntimeError(
            "SolvationFepDriver lacks skip_vacuum/coulomb_global_parameters "
            "support")


def is_rigid_cluster(topology) -> bool:
    """True if every atom is within two bonds of every other (pure 1-2/1-3 topology).

    Such a cluster has no intramolecular nonbonded terms, so its vacuum
    decoupling free energy is zero.
    """
    neighbors = collections.defaultdict(list)
    for bond in topology.bonds():
        neighbors[bond.atom1.index].append(bond.atom2.index)
        neighbors[bond.atom2.index].append(bond.atom1.index)
    n_atoms = topology.getNumAtoms()
    if n_atoms <= 3:
        return True
    for start in range(n_atoms):
        depth = {start: 0}
        queue = collections.deque([start])
        while queue:
            current = queue.popleft()
            if depth[current] > 2:
                return False
            for neighbor in neighbors[current]:
                if neighbor not in depth:
                    depth[neighbor] = depth[current] + 1
                    queue.append(neighbor)
    return len(depth) == n_atoms


def _make_driver(fep_options: dict[str, Any], ostream) -> SolvationFepDriver:
    """``SolvationFepDriver`` with the production protocol, overridden by ``fep_options``.

    ``DRIVER_OPTIONS``, the lambdas and ``timestep_fs`` are applied; other keys are
    ignored here.
    """
    settings = fep_settings(fep_options)
    driver = (_QuickMbarFepDriver(ostream=ostream)
              if settings.get("quick_dirty") else
              SolvationFepDriver(ostream=ostream))
    for name in DRIVER_OPTIONS:
        if name in settings:
            setattr(driver, name, settings[name])
    driver.lambdas_stage1 = list(settings["lambdas_stage1"])
    driver.lambdas_stage2 = list(settings["lambdas_stage2"])
    driver.timestep = float(settings["timestep_fs"]) * unit.femtoseconds
    return driver


def run_fep(systems, final_state: mm.State | None, work_dir: Path,
            fep_options: dict[str, Any], batch, *, ostream) -> dict[str, Any]:
    """Run the two-stage solvation FEP and return the VeloxChem output dictionary.

    ``systems`` is the calculation's ``OpenMMSystems``; the solvated box (at
    ``final_state`` if given) and the candidate are written to a temporary
    directory under ``work_dir`` for VeloxChem's file-based backend. Raises
    ``RuntimeError`` unless the CDAM cluster is rigid (purely 1-2/1-3 bonded),
    because only then are the vacuum legs zero and skipped.
    """
    default_platform()  # CPU thread default for the FEP contexts
    ostream.print_info("Handing CDAM system to SolvationFepDriver...")
    ostream.flush()
    system = systems.solvated_system
    driver = _make_driver(fep_options, ostream)
    driver.ostream = _FepStatus(ostream, driver)
    PROGRESS_STATUS.set("FEP setup")

    c4_active = any(
        force.getName() == C4_FORCE_NAME for force in system.getForces())
    if c4_active:
        # The 12-6-4 C4 term is a separate CustomNonbondedForce,
        # -(1-c4_lambda)*C4_i*C4_j/r^4, unknown to SolvationFepDriver. It is
        # an induced-dipole (charge-driven) attraction, so it is switched off
        # with the charges: c4_lambda follows lambda in stage 1 and stays 1.0
        # (off) in stage 2. Without this the decoupled ion keeps an
        # unopposed -C4/r^4 pull on water, giving a wrong HFE or a crash.
        driver.coulomb_global_parameters = ["c4_lambda"]

    if not is_rigid_cluster(systems.vacuum_topology):
        raise RuntimeError(
            "file-based FEP needs a rigid (1-2/1-3 bonded) CDAM cluster")
    # Every metal-dummy and dummy-dummy pair is a 1-2/1-3 exclusion and C4
    # acts only on water oxygens, so the vacuum legs (stages 3 and 4) have
    # nothing to decouple and are exactly 0.0 kJ/mol; skip them to save time.
    ostream.print_info("CDAM cluster is purely 1-2/1-3 bonded. "
                       "Bypassing vacuum stages (dG = 0.0 kJ/mol).")
    driver.skip_vacuum = True

    step_ps = driver.timestep.value_in_unit(unit.picoseconds)
    print_settings(ostream, "Solvation FEP Setup", [
        ("Backend", "SolvationFepDriver (OpenMM files)"),
        ("Protocol", "quick_dirty (workflow test only)"
         if fep_options.get("quick_dirty") else "production"),
        ("Temperature", f"{getattr(driver, 'temperature', 298.15)} K"),
        ("Time Step",
         f"{driver.timestep.value_in_unit(unit.femtoseconds):g} fs"),
        ("Minimization Steps", getattr(driver, "num_em_steps", None)),
        ("Equilibration per Lambda", f"{driver.num_equil_steps} steps "
         f"({driver.num_equil_steps * step_ps:g} ps)" if hasattr(
             driver, "num_equil_steps") else None),
        ("Sampling per Lambda",
         f"{driver.num_steps} steps ({driver.num_steps * step_ps:g} ps)"),
        ("Snapshots per Lambda", driver.num_snapshots),
        ("Stage 1 Lambdas", driver.lambdas_stage1),
        ("Stage 2 Lambdas", driver.lambdas_stage2),
        ("12-6-4 C4 Term Scaled", c4_active),
        ("Vacuum Stages", "skipped (rigid 1-2/1-3 cluster)"),
        ("Platform", getattr(driver, "platform", None) or
         "SolvationFepDriver default"),
        ("Custom FEP Options", fep_options or None),
    ])

    work_dir.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="autocdam-fep-",
                                     dir=work_dir) as temp_name:
        temp_dir = Path(temp_name)
        pdb_path, xml_path = CDAMCandidate.from_generated_batch(
            batch).write_files(temp_dir, stem="solute")
        system_pdb, system_xml = temp_dir / "system.pdb", temp_dir / "system.xml"
        positions = (final_state.getPositions()
                     if final_state is not None else systems.solvated_positions)
        with system_pdb.open("w", encoding="utf-8") as output:
            app.PDBFile.writeFile(systems.solvated_topology,
                                  positions,
                                  output,
                                  keepIds=True)
        if not any(
                isinstance(force, mm.MonteCarloBarostat)
                for force in system.getForces()):
            system.addForce(
                mm.MonteCarloBarostat(PRESSURE_ATM * unit.atmospheres,
                                      TEMPERATURE_K * unit.kelvin, 25))
        system_xml.write_text(mm.XmlSerializer.serialize(system),
                              encoding="utf-8")
        with _PymbarLog() as mbar_log:
            output = driver.compute_solvation_from_openmm_files(
                solute_pdb=str(pdb_path),
                solute_xml=str(xml_path),
                system_pdb=str(system_pdb),
                system_xml=str(system_xml))
        usable, summary = mbar_log.summary()
        output["mbar"] = {"usable": usable, "summary": summary,
                          "solves": mbar_log.solves,
                          "unconverged_solves": mbar_log.unconverged}
        (ostream.print_info if usable else ostream.print_warning)(summary)
        ostream.flush()
        return output


class _QuickMbarFepDriver(SolvationFepDriver):
    """
    SolvationFepDriver with one MBAR solve per stage (analytical uncertainty)
    instead of 50 bootstrap solves, for the quick_dirty workflow test.
    """

    def _calculate_free_energy(self, u_kn):
        from pymbar import MBAR

        mbar = MBAR(u_kn, N_k=self.N_k, solver_protocol='robust')
        delta_f = mbar.compute_free_energy_differences()
        delta_f['overlap_matrix'] = np.asarray(mbar.compute_overlap()['matrix'],
                                               dtype=float)
        delta_f['effective_sample_sizes'] = np.asarray(
            mbar.compute_effective_sample_number(), dtype=float)
        return delta_f


class _FepStatus:
    """
    Output stream of the FEP driver that also tells the progress bar which
    stage and lambda window is running.
    """

    _PATTERNS = (
        (re.compile(r"Generating systems for stage (\d+)"), "setup"),
        (re.compile(r"Running lambda = (\S+), stage = (\d+)"), "lambda"),
        (re.compile(r"Recalculating energies for stage (\d+)"), "energies"),
        (re.compile(r"free energy with MBAR for stage (\d+)"), "MBAR"),
    )

    def __init__(self, ostream, driver):
        self._ostream = ostream
        self._driver = driver

    def __getattr__(self, name):
        return getattr(self._ostream, name)

    def print_info(self, line):
        for pattern, step in self._PATTERNS:
            match = pattern.search(line)
            if match is None:
                continue
            if step == "lambda":
                stage = int(match.group(2))
                lambdas = getattr(self._driver, f"lambdas_stage{stage}", [])
                index = next((i for i, lam in enumerate(lambdas)
                              if f"{lam}" == match.group(1).rstrip(".,")),
                             None)
                step = (f"lambda {index + 1}/{len(lambdas)}"
                        if index is not None else f"lambda {match.group(1)}")
            else:
                stage = int(match.group(1))
            PROGRESS_STATUS.set(f"FEP stage {stage}/2: {step}")
            break
        self._ostream.print_info(line)


class _PymbarLog(logging.Handler):
    """
    Keeps pymbar's solver messages (and its JAX notice) off the console during
    the FEP and sums them up in one conclusion.
    """

    _GRADIENT = re.compile(r"smallest gradient (\S+) norm")
    _FINAL = re.compile(r"Final gradient norm: (\S+)")
    _DELTA = re.compile(r"max_delta = (\S+),")

    def __init__(self):
        super().__init__(logging.INFO)
        self.solves = 0
        self.unconverged = 0
        self.gradients = []
        self.deltas = []
        self._logger = logging.getLogger("pymbar")

    def emit(self, record):
        message = record.getMessage()
        if self._FINAL.search(message):
            self.solves += 1
        if "No solution found" in message:
            self.unconverged += 1
        for pattern, values in ((self._GRADIENT, self.gradients),
                                (self._DELTA, self.deltas)):
            match = pattern.search(message)
            if match:
                try:
                    values.append(float(match.group(1)))
                except ValueError:
                    pass

    def summary(self):
        """
        Gets the conclusion of the MBAR solves.

        :return:
            (usable, text): whether the free energies can be used, and why.
        """

        solves = max(self.solves, self.unconverged)
        if not self.unconverged:
            return True, f"MBAR converged ({solves} solve(s))"
        change = max(self.deltas, default=0.0)
        gradient = max(self.gradients, default=float("nan"))
        if change <= 1e-10:
            return True, (
                f"MBAR usable: {self.unconverged} of {solves} solve(s) ended "
                f"at gradient {gradient:.1e}, above pymbar's 1e-12 tolerance, "
                "but the free energies had stopped changing")
        return False, (
            f"MBAR not converged: {self.unconverged} of {solves} solve(s), "
            f"free energies still changing by {change:.1e}; check the HFE")

    def __enter__(self):
        self._propagate = self._logger.propagate
        self._level = self._logger.level
        self._logger.addHandler(self)
        self._logger.propagate = False
        self._logger.setLevel(logging.INFO)
        return self

    def __exit__(self, *exc):
        self._logger.removeHandler(self)
        self._logger.propagate = self._propagate
        self._logger.setLevel(self._level)
        return False


def finite_size_correction(formal_charge: float, box_length_nm: float,
                           r_ion_nm: float, solvent_density_nm3: float,
                           water_model: str) -> float:
    """Periodic-box (Eq. 5) plus summation-scheme (Eq. 6) correction, kJ/mol."""
    q, L, r1, rho = float(formal_charge), float(box_length_nm), float(
        r_ion_nm), float(solvent_density_nm3)
    k_e, eps = COULOMB_CONSTANT, WATER_DIELECTRIC
    gamma = WATER_QUADRUPOLE_TRACE.get(water_model.lower(),
                                       WATER_QUADRUPOLE_TRACE["tip3p"])

    dg_pbc = (k_e * q**2 /
              (2.0 * L)) * (-WIGNER_LATTICE_SUM / eps + (1.0 - 1.0 / eps) *
                            ((4.0 * math.pi / 3.0) * (r1 / L)**2 -
                             (16.0 * math.pi**2 / 45.0) * (r1 / L)**5))
    ion_fraction = (4.0 * math.pi * r1**3) / (3.0 * L**3)
    dg_summation = (-q * (1.0 - ion_fraction) *
                    ((2.0 * math.pi / 3.0) * k_e * rho * gamma) -
                    q * ion_fraction * (16.8 - 2.6 / r1))
    return dg_pbc + dg_summation


def fep_uncertainty(fep_output: dict[str, Any]) -> float | None:
    """Stage uncertainties added in quadrature, else the output's own estimate."""
    stages = [
        float(stage["Uncertainty"])
        for key, stage in fep_output.items()
        if key.startswith("Stage ") and "Uncertainty" in stage
    ]
    if stages:
        return math.sqrt(sum(value * value for value in stages))
    if "uncertainty" in fep_output:
        return float(fep_output["uncertainty"])
    return None
