"""Dev tests (WP-E3, defects E3-2 and E3-3): the periodic box and the
thermodynamic ensemble of every MD stage of EvbFepDriver.run_replicas.

E3-2: the box was not carried from one stage to the next. Every window started
in the system's default box with positions from the previous (NPT) box, which
squeezed DhlA by 0.5 % per edge (+3e5 kJ/mol) and crashed 3 of ~540 window
starts ('Particle coordinate is NaN').

E3-3: the window ensemble did not follow `isobaric`. The systems carry a
MonteCarloBarostat at frequency 25 and only the initial equilibration changed
it, on the shared system objects, so every window except lambda 0 ran NPT
while `isobaric` was False, the lambda 0 window ran NVT, and the PDB warm-up
ran NPT at the target temperature's barostat.

A spy on EvbFepDriver._safe_step records, for every stage in the order it runs,
the box before and after it and the barostat frequency in effect.
"""

import os
from pathlib import Path

import numpy as np
import pytest

pytest.importorskip('openmm')

import openmm as mm
import openmm.app as mmapp
import openmm.unit as mmunit

from veloxchem.evbdriver import EvbDriver
from veloxchem.evbfepdriver import EvbFepDriver
from veloxchem.outputstream import OutputStream

from test_evb_helper import EvbTestHelper

LAMBDAS = [0.0, 0.5, 1.0]
N_REPLICAS = 2
BAROSTAT_TYPES = (mm.MonteCarloBarostat, mm.MonteCarloFlexibleBarostat)


def _box(simulation):
    state = simulation.context.getState()
    return np.array(
        state.getPeriodicBoxVectors(asNumpy=True).value_in_unit(
            mmunit.nanometer))


def _frequency(simulation):
    freqs = [
        f.getFrequency()
        for f in simulation.system.getForces()
        if isinstance(f, BAROSTAT_TYPES)
    ]
    return freqs[0] if freqs else None


class StageSpy:
    """Wraps EvbFepDriver._safe_step and records every MD stage."""

    def __init__(self):
        self.stages = []
        self.driver = None
        self._original = EvbFepDriver._safe_step

    def install(self, patch):
        spy = self

        def wrapped(driver, simulation, steps, name=""):
            spy.driver = driver
            temperature = simulation.integrator.getTemperature()
            record = {
                "name": name,
                "steps": steps,
                "box0": _box(simulation),
                "frequency": _frequency(simulation),
                "temperature": temperature.value_in_unit(mmunit.kelvin),
            }
            result = spy._original(driver, simulation, steps, name)
            record["box1"] = _box(simulation)
            spy.stages.append(record)
            return result

        patch.setattr(EvbFepDriver, "_safe_step", wrapped)


def _config(name, **overrides):
    config = {
        "name": name,
        "temperature": 300.0,
        # Step counts are the minimum: the barostat attempts a move every 25
        # steps, so an NPT stage needs 50 for two attempts.
        "sample_steps": 25,
        "write_step": 5,
        "equil_NVT_steps": 25,
        "equil_NPT_steps": 50,
        "initial_equil_NVT_steps": 50,
        "initial_equil_NPT_steps": 250,
        "minimize_every_lambda": False,
        "step_size": 0.001,
        "n_replicas": 1,
        "report_integration_energies": False,
    }
    config.update(overrides)
    return config


def _seeded_driver():
    pair = EvbTestHelper.evb_ff_pair()
    EVB = EvbDriver(ostream=OutputStream(None))
    EVB.water_model = "cspce"  # flexible: no constrained minimisation
    # Identical states: every lambda system is the same Hamiltonian, so the
    # windows need no minimisation (a reacting pair blows up on the 0 -> 0.5
    # jump without minimize_every_lambda, WP-B4 trap) and the test is cheap.
    EvbTestHelper.evb_seed_states(EVB, [pair.reactant, pair.reactant])
    return EVB


def _in_folder(folder, function):
    previous = Path.cwd()
    os.chdir(folder)
    try:
        return function()
    finally:
        os.chdir(previous)


def _run(folder, preset, run_overrides=None, **overrides):
    """Build and run one solvated configuration under the spy."""

    def run():
        EVB = _seeded_driver()
        conf = EVB.default_system_configurations(preset)
        conf.update(_config(preset, padding=1.0, **overrides))
        EVB.build_systems(configurations=[conf], Lambda=LAMBDAS)
        conf = EVB.system_confs[0]
        for c in [conf] + conf.get("steps", []):
            c.update(run_overrides or {})
        spy = StageSpy()
        with pytest.MonkeyPatch.context() as patch:
            spy.install(patch)
            EVB.run_FEP(platform=EvbTestHelper.evb_platform())
        return EVB, spy

    return _in_folder(folder, run)


@pytest.fixture(scope="module")
def nvt_run(tmp_path_factory):
    """water_NVT preset: pressure 1 (barostat in the systems), isobaric off,
    one replica."""
    return _run(tmp_path_factory.mktemp("box_nvt"), "water_NVT")


@pytest.fixture(scope="module")
def npt_run(tmp_path_factory):
    """water_NPT preset (isobaric on), through the PDB warm-up branch.

    `pdb` is only a flag inside EvbFepDriver (it selects the warm-up), so it is
    set after the build on a water box. Two replicas, so the hand-over from
    one replica's backward sweep to the next forward sweep is covered."""
    warmup = {
        "pdb": "flag-only",
        "pdb_temperatures": [150.0, 300.0],
        "warmup_NVT_steps": 50,
    }
    return _run(tmp_path_factory.mktemp("box_npt"),
                "water_NPT",
                n_replicas=N_REPLICAS,
                run_overrides=warmup)


def _is_initial(stage):
    name = stage["name"]
    return name.startswith("initial") or "warmup" in name or (
        "target temperature" in name)


def _window_stages(spy):
    return [s for s in spy.stages if not _is_initial(s)]


def _assert_chain(stages):
    for k in range(1, len(stages)):
        np.testing.assert_allclose(
            stages[k]["box0"],
            stages[k - 1]["box1"],
            rtol=1e-12,
            atol=1e-12,
            err_msg=(f"stage {k} ({stages[k]['name']}) did not start in the "
                     f"box stage {k - 1} ({stages[k - 1]['name']}) ended in"))


@pytest.mark.parametrize("which", ["nvt", "npt"])
def test_box_is_carried_through_every_stage(which, request):
    """E3-2: each stage starts in the box the previous one ended in: initial
    equilibration -> windows -> sweeps -> replicas."""
    _, spy = request.getfixturevalue(f"{which}_run")
    n_windows = len(LAMBDAS) * 2 * spy.driver.n_replicas
    assert len(_window_stages(spy)) == 2 * n_windows
    _assert_chain(spy.stages)


@pytest.mark.parametrize("which", ["nvt", "npt"])
def test_initial_npt_changed_the_box(which, request):
    """Precondition for the chain test: the initial NPT stage moved the box
    away from the built one, so a reset would be seen."""
    _, spy = request.getfixturevalue(f"{which}_run")
    npt = [s for s in spy.stages if _is_initial(s) and "NPT" in s["name"]]
    assert len(npt) == 1
    assert npt[0]["frequency"] == 25
    assert not np.allclose(npt[0]["box1"], npt[0]["box0"], rtol=0, atol=1e-9)


def test_nvt_windows_run_at_constant_volume(nvt_run):
    """E3-3: isobaric False -> barostat off in every window stage (equilibration
    and sampling, every lambda), at the volume the initial NPT ended at."""
    _, spy = nvt_run
    initial = [s for s in spy.stages if _is_initial(s)]
    assert [s["frequency"] for s in initial] == [0, 25]
    windows = _window_stages(spy)
    assert all(s["frequency"] == 0 for s in windows), [
        (s["name"], s["frequency"]) for s in windows
    ]
    for s in windows:
        np.testing.assert_allclose(s["box0"], initial[-1]["box1"], atol=1e-12)
        np.testing.assert_allclose(s["box1"], initial[-1]["box1"], atol=1e-12)


def test_npt_windows_run_at_constant_pressure(npt_run):
    """E3-3: isobaric True -> barostat on in every window stage, equilibration
    (equil_NPT_steps) and sampling alike, and the volume moves."""
    _, spy = npt_run
    windows = _window_stages(spy)
    assert all(s["frequency"] == 25 for s in windows), [
        (s["name"], s["frequency"]) for s in windows
    ]
    equil = [s for s in windows if s["name"] != "sampling"]
    assert all(s["name"] == "NPT equilibration" for s in equil)
    assert all(s["steps"] == 50 for s in equil)
    moved = [
        not np.allclose(s["box1"], s["box0"], rtol=0, atol=1e-9)
        for s in windows
    ]
    assert sum(moved) >= 1


def test_pdb_warmup_is_nvt_then_npt_at_target(npt_run):
    """E3-3: the PDB warm-up stages below the target temperature run NVT
    (barostat off: its temperature is the target's), then one NPT stage at the
    target temperature."""
    _, spy = npt_run
    initial = [s for s in spy.stages if _is_initial(s)]
    assert [s["temperature"] for s in initial] == [150.0, 300.0, 300.0]
    assert [s["frequency"] for s in initial] == [0, 0, 25]


@pytest.mark.parametrize("which", ["nvt", "npt"])
def test_step_accounting_matches_the_steps_run(which, request):
    """The progress/ETA totals count exactly the steps that are run."""
    _, spy = request.getfixturevalue(f"{which}_run")
    driver = spy.driver
    assert driver.run_steps == sum(s["steps"] for s in spy.stages)
    assert driver.total_steps == driver.run_steps


@pytest.mark.parametrize("which", ["nvt", "npt"])
def test_system_objects_keep_their_barostat_frequency(which, request):
    """Stages set the frequency they need; the built systems are left as they
    were built (25), whatever ran before."""
    EVB, _ = request.getfixturevalue(f"{which}_run")
    systems = EVB.system_confs[0]["systems"]
    for key, system in systems.items():
        freqs = [
            f.getFrequency()
            for f in system.getForces()
            if isinstance(f, BAROSTAT_TYPES)
        ]
        assert freqs == [25], (key, freqs)


def test_isobaric_without_positive_pressure_is_refused(tmp_path):
    """isobaric switched on at run time on a solvated build without a pressure
    (barostat at -1 bar; the builder adds one to every solvated system) stops
    with a clear error. Vacuum + isobaric is refused by the builder."""

    def run():
        EVB = _seeded_driver()
        conf = EVB.default_system_configurations("water_NVT")
        conf.update(_config("water", padding=1.0, pressure=-1))
        EVB.build_systems(configurations=[conf], Lambda=LAMBDAS)
        conf = EVB.system_confs[0]
        for c in [conf] + conf.get("steps", []):
            c["isobaric"] = True
        EVB.run_FEP(platform=EvbTestHelper.evb_platform())

    with pytest.raises(Exception, match="isobaric is set"):
        _in_folder(tmp_path, run)


def test_restart_starts_in_the_restart_box(tmp_path, nvt_run):
    """A restart from a cif (load_initialisation restart_pdb) starts in the
    cell of that file, not in the built box."""
    EVB, spy = nvt_run
    data_folder = Path(EVB.system_confs[0]["data_folder"])
    if not data_folder.is_absolute():
        data_folder = spy.driver.data_folder
    box = spy.stages[-1]["box1"] * 1.01
    topology = mmapp.PDBxFile(str(data_folder / "topology.cif")).getTopology()
    positions = mmapp.PDBxFile(str(data_folder /
                                   "topology.cif")).getPositions(asNumpy=True)
    topology.setPeriodicBoxVectors(box * mmunit.nanometer)
    with open(tmp_path / "restart.cif", "w") as f:
        mmapp.PDBxFile.writeFile(topology, positions, f)
    import shutil
    copy = tmp_path / "water"
    shutil.copytree(data_folder, copy)
    shutil.copy(tmp_path / "restart.cif", copy / "restart.cif")

    def run():
        EVB2 = EvbDriver(ostream=OutputStream(None))
        EVB2.load_initialisation(str(copy),
                                 load_systems=True,
                                 restart_pdb="restart.cif",
                                 rename_data=True)
        conf = EVB2.system_confs[0]
        for c in [conf] + conf["steps"]:
            c.update({
                "n_replicas": 1,
                "sample_steps": 4,
                "write_step": 2,
                "equil_NVT_steps": 2,
                # the cif holds the raw solvated build (overlapping waters)
                "minimize_every_lambda": True,
                "report_integration_energies": False,
            })
        spy2 = StageSpy()
        with pytest.MonkeyPatch.context() as patch:
            spy2.install(patch)
            EVB2.run_FEP(platform=EvbTestHelper.evb_platform())
        return spy2

    spy2 = _in_folder(tmp_path, run)
    # the cif cell has 1e-4 Angstrom resolution
    np.testing.assert_allclose(spy2.stages[0]["box0"], box, atol=1e-4)
    _assert_chain(spy2.stages)


def test_saved_cif_carries_the_state_box(tmp_path, nvt_run):
    """_save_state writes the cell of the state it saves, so a restart from
    one of its cif files gets the right box."""
    EVB, spy = nvt_run
    driver = spy.driver
    system = EVB.system_confs[0]["systems"][0.5]
    simulation = driver._get_simulation(system, 0.001)
    box = spy.stages[-1]["box1"] * 1.02
    simulation.context.setPeriodicBoxVectors(*(box * mmunit.nanometer))
    simulation.context.setPositions(
        mmapp.PDBxFile(str(driver.data_folder / "topology.cif")).getPositions())
    old = driver.run_folder
    driver.run_folder = tmp_path
    try:
        driver._save_state(simulation, "probe", xml=False, chk=False, cif=True)
    finally:
        driver.run_folder = old
    read = mmapp.PDBxFile(str(tmp_path / "probe.cif")).getTopology()
    cell = np.array(
        read.getPeriodicBoxVectors().value_in_unit(mmunit.nanometer))
    np.testing.assert_allclose(cell, box, atol=1e-4)
    built = np.array(driver.topology.getPeriodicBoxVectors().value_in_unit(
        mmunit.nanometer))
    assert not np.allclose(built, box, atol=1e-4)
