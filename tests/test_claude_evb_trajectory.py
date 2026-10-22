"""Dev tests (WP-B4, defect B3-2, feature F4): trajectory frames in full
precision for the deferred recalculation.

OpenMM's XTC writer rounds positions to 0.001 nm, so every energy the deferred
recalculation (and compute_force_groups) evaluates on trajectory.xtc carries
noise. trajectory_format="dcd" writes float32 frames instead, read back through
EvbFepDriver._open_trajectory.

One solvated run writes trajectory.dcd plus a shadow xtc of the same states (a
test-only reporter writes both), so the xtc and dcd frames compare directly.
The reference is Ep in Data_combined.csv: the StateDataReporter's potential
energy of the sampling context at the same step, from exact positions.
"""

import os
import shutil
from pathlib import Path

import numpy as np
import pytest

pytest.importorskip('openmm')
pytest.importorskip('MDAnalysis')

import openmm as mm
import openmm.app as mmapp
import openmm.unit as mmunit
from MDAnalysis.lib.formats.libmdaxdr import XTCFile

from veloxchem import evbfepdriver
from veloxchem.evbdriver import EvbDriver
from veloxchem.evbfepdriver import EvbFepDriver
from veloxchem.outputstream import OutputStream

from test_evb_helper import EvbTestHelper

LAMBDAS = [0.0, 0.5, 1.0]
N_REPLICAS = 2
FRAMES_PER_WINDOW = 10


class ShadowDCDReporter(mmapp.DCDReporter):
    """DCDReporter that also writes the same states to <name>_shadow.xtc."""

    def __init__(self, file, reportInterval, **kwargs):
        super().__init__(file, reportInterval, **kwargs)
        self._xtc = mmapp.XTCReporter(file.replace(".dcd", "_shadow.xtc"),
                                      reportInterval)

    def report(self, simulation, state):
        super().report(simulation, state)
        self._xtc.report(simulation, state)


def _config(name, **overrides):
    config = {
        "name": name,
        "temperature": 300.0,
        "sample_steps": 5 * FRAMES_PER_WINDOW,
        "write_step": 5,
        "equil_NVT_steps": 20,
        "equil_NPT_steps": 0,
        "initial_equil_NVT_steps": 50,
        "initial_equil_NPT_steps": 0,
        "minimize_every_lambda": True,
        "step_size": 0.001,
        "n_replicas": N_REPLICAS,
        "trajectory_format": "dcd",
    }
    config.update(overrides)
    return config


def _seeded_driver():
    pair = EvbTestHelper.evb_ff_pair()
    EVB = EvbDriver(ostream=OutputStream(None))
    EVB.water_model = "cspce"
    EvbTestHelper.evb_seed_states(EVB, [pair.reactant, pair.product])
    return EVB


def _in_folder(folder, function):
    previous = Path.cwd()
    os.chdir(folder)
    try:
        return function()
    finally:
        os.chdir(previous)


@pytest.fixture(scope="module")
def water_run(tmp_path_factory):
    """Solvated ethanol (~1300 atoms), 3 lambda, 2 replicas, both directions."""
    folder = tmp_path_factory.mktemp("traj_water")

    def run():
        EVB = _seeded_driver()
        conf = EVB.default_system_configurations("water_NVT")
        conf.update(_config("water", padding=1.0))
        EVB.build_systems(configurations=[conf], Lambda=LAMBDAS)
        with pytest.MonkeyPatch.context() as patch:
            patch.setattr(evbfepdriver.mmapp, "DCDReporter", ShadowDCDReporter)
            EVB.run_FEP(platform=EvbTestHelper.evb_platform())
        return EVB, folder / EVB.system_confs[0]["data_folder"]

    EVB, data_folder = _in_folder(folder, run)
    return folder, EVB, data_folder


@pytest.fixture(scope="module")
def vacuum_run(tmp_path_factory):
    folder = tmp_path_factory.mktemp("traj_vacuum")

    def run():
        EVB = _seeded_driver()
        EVB.build_systems(configurations=[
            _config("vacuum", n_replicas=1, minimize_every_lambda=False)
        ],
                          Lambda=LAMBDAS)
        EVB.run_FEP(platform=EvbTestHelper.evb_platform())
        return folder / EVB.system_confs[0]["data_folder"]

    return _in_folder(folder, run)


def _read_frames(path, periodic):
    with EvbFepDriver._open_trajectory(path, periodic) as (n_frames, read):
        return n_frames, list(read(0, n_frames))


def _read_shadow_xtc(path):
    with XTCFile(str(path), mode='r') as reader:
        return [(np.array(f.x, dtype=np.float64),
                 np.array(f.box, dtype=np.float64)) for f in reader]


def _energy_rows(data_folder):
    data = np.genfromtxt(data_folder / "Energies.csv",
                         delimiter=",",
                         skip_header=1)
    return data.reshape(-1, 8)


def _ep(data_folder):
    # Data_combined.csv: Step, Potential Energy, ...; one row per frame.
    return np.genfromtxt(data_folder / "Data_combined.csv",
                         delimiter=",",
                         skip_header=1)[:, 1]


def _own_lambda_energies(data_folder, frames, lambdas):
    """Potential energy of each frame in its own lambda's integration system,
    on the platform the run used."""
    platform_name = EvbTestHelper.evb_platform()
    contexts = {}
    energies = []
    for (pos, box), lam in zip(frames, lambdas):
        if lam not in contexts:
            xml = (data_folder / "run" / f"{lam:.3f}_sys.xml").read_text()
            system = mm.XmlSerializer.deserialize(xml)
            integrator = mm.VerletIntegrator(0.001)
            if platform_name is None:
                contexts[lam] = mm.Context(system, integrator)
            else:
                contexts[lam] = mm.Context(
                    system, integrator,
                    mm.Platform.getPlatformByName(platform_name))
        context = contexts[lam]
        if box is not None:
            context.setPeriodicBoxVectors(*box)
        context.setPositions(pos)
        energies.append(context.getState(getEnergy=True).getPotentialEnergy(
        ).value_in_unit(mmunit.kilojoule_per_mole))
    return np.array(energies)


@pytest.fixture(scope="module")
def own_lambda(water_run):
    """Own-lambda energies minus Ep, for the xtc and the dcd frames."""
    _, _, data_folder = water_run
    lambdas = _energy_rows(data_folder)[:, 0]
    ep = _ep(data_folder)
    _, dcd = _read_frames(data_folder / "trajectory.dcd", True)
    xtc = _read_shadow_xtc(data_folder / "trajectory_shadow.xtc")
    return {
        "xtc": _own_lambda_energies(data_folder, xtc, lambdas) - ep,
        "dcd": _own_lambda_energies(data_folder, dcd, lambdas) - ep,
        "lag": _own_lambda_energies(data_folder, xtc[:-1], lambdas[:-1])
               - ep[1:],
    }


def test_q1_xtc_frames_are_rounded(own_lambda):
    """B3-2: energies on xtc frames scatter by kJ/mol around Ep. The one-row
    lag is far worse, so the rows are aligned."""
    xtc = own_lambda["xtc"]
    print(f"xtc - Ep: mean {xtc.mean():.3f} std {xtc.std():.3f} "
          f"max |dev| {np.abs(xtc - xtc.mean()).max():.3f} kJ/mol")
    print(f"lag - Ep: std {own_lambda['lag'].std():.1f} kJ/mol")
    assert xtc.std() >= 1.0
    assert own_lambda["lag"].std() > 10 * xtc.std()


def test_q2_dcd_frames_reproduce_ep(own_lambda):
    dcd = own_lambda["dcd"]
    std = dcd.std()
    worst = np.abs(dcd).max()
    print(f"dcd - Ep: mean {dcd.mean():.2e} std {std:.2e} "
          f"max {worst:.2e} kJ/mol")
    assert std <= max(0.05, own_lambda["xtc"].std() / 300)
    assert worst <= 10 * max(std, 1e-6)


def test_q2_recalculated_integration_energy_is_ep(water_run):
    """The code's own Energies.csv: at lambda 0 (1) the reactant (product)
    integration energy is the sampled energy of that frame."""
    _, _, data_folder = water_run
    rows = _energy_rows(data_folder)
    ep = _ep(data_folder)
    for lam, column in ((0.0, 3), (1.0, 4)):
        sel = rows[:, 0] == lam
        diff = rows[sel, column] - ep[sel]
        assert np.abs(diff).max() < 0.05, (lam, np.abs(diff).max())


def test_q3_dcd_round_trip_against_xtc(water_run):
    _, _, data_folder = water_run
    n_frames, dcd = _read_frames(data_folder / "trajectory.dcd", True)
    xtc = _read_shadow_xtc(data_folder / "trajectory_shadow.xtc")
    assert n_frames == len(xtc)
    for (pos_d, box_d), (pos_x, box_x) in zip(dcd, xtc):
        assert pos_d.dtype == np.float64 and pos_d.shape == pos_x.shape
        assert box_d.shape == (3, 3)
        assert np.abs(pos_d - pos_x).max() <= 5e-4 + 1e-6
        assert np.abs(box_d - box_x).max() <= 1e-5


def test_q3_vacuum_box_is_none(vacuum_run):
    n_frames, frames = _read_frames(vacuum_run / "trajectory.dcd", False)
    assert n_frames == len(LAMBDAS) * 2 * FRAMES_PER_WINDOW
    assert all(box is None for _, box in frames)
    assert np.isfinite(_energy_rows(vacuum_run)[:, 1:6]).all()
    assert not (vacuum_run / "trajectory.xtc").exists()


def test_q4_frame_count_matches_energy_rows(water_run):
    _, EVB, data_folder = water_run
    rows = _energy_rows(data_folder)
    n_frames, _ = _read_frames(data_folder / "trajectory.dcd", True)
    assert n_frames == len(rows) == (len(LAMBDAS) * 2 * N_REPLICAS *
                                     FRAMES_PER_WINDOW)
    assert set(rows[:, 6]) == set(range(N_REPLICAS))
    assert set(rows[:, 7]) == {0, 1}
    assert not (data_folder / "trajectory.xtc").exists()
    with open(data_folder / "options.json") as f:
        assert '"trajectory_format": "dcd"' in f.read()


def test_q4_compute_force_groups_on_dcd(water_run):
    folder, EVB, data_folder = water_run
    _in_folder(folder, lambda: EVB.compute_force_groups(
        platform=EvbTestHelper.evb_platform()))
    n_rows = len(_energy_rows(data_folder))
    for name in ("ForceGroups_rea.csv", "ForceGroups_pro.csv"):
        data = np.genfromtxt(data_folder / name, delimiter=",", skip_header=1)
        assert len(data) == n_rows
        assert np.isfinite(data).all()


def test_q4_reloaded_run_appends(water_run, tmp_path):
    """The B3-1 path: reload with restart_pdb from the last sampled frame.
    The old trajectory is moved aside (rename_data), the new one holds exactly
    the new run's frames, starting at the restart structure."""
    _, _, data_folder = water_run
    copy = tmp_path / "water"
    shutil.copytree(data_folder, copy)
    n_old, frames = _read_frames(copy / "trajectory.dcd", True)
    restart_nm, box = frames[-1]
    topology = mmapp.PDBxFile(str(copy / "topology.cif")).getTopology()
    topology.setPeriodicBoxVectors(box)
    with open(copy / "restart.cif", "w") as f:
        mmapp.PDBxFile.writeFile(topology, restart_nm * 10 * mmunit.angstrom,
                                 f)

    def run():
        EVB = EvbDriver(ostream=OutputStream(None))
        EVB.load_initialisation(str(copy),
                                load_systems=True,
                                restart_pdb="restart.cif",
                                rename_data=True)
        conf = EVB.system_confs[0]
        short = {
            "n_replicas": 1,
            "sample_steps": 4,
            "write_step": 2,
            "equil_NVT_steps": 2,
            "minimize_every_lambda": False,
        }
        for c in [conf] + conf["steps"]:
            c.update(short)
        EVB.run_FEP(platform=EvbTestHelper.evb_platform())

    _in_folder(tmp_path, run)
    assert (copy / "OLD_trajectory.dcd").exists()
    rows = _energy_rows(copy)
    n_new, new_frames = _read_frames(copy / "trajectory.dcd", True)
    assert n_new == len(rows) == len(LAMBDAS) * 2 * 2
    assert np.isfinite(rows[:, 1:6]).all()
    # RMS minimum-image displacement of the first new frame (4 fs later) from
    # the restart; a 10x shrunk start would be off by nanometres.
    lengths = np.diag(box)
    delta = new_frames[0][0] - restart_nm
    delta -= lengths * np.round(delta / lengths)
    assert np.sqrt((delta**2).sum(axis=1).mean()) < 0.01
