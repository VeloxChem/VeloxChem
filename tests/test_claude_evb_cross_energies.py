"""Dev tests (WP-B5a, feature F5): cross energies of every frame in every
integration system of the ladder.

With report_integration_energies, the deferred recalculation also scores each
frame in the integration system of every lambda of the step and writes
IntegrationEnergies.csv, row-aligned with Energies.csv. The loader reads it
into configuration_results[name]["U_int"] (n_frames x len(Lambda)), which MBAR
needs to treat the integration systems as the sampled states (B5b).

Checks: the file layout (Q1), the endpoint columns against Energies.csv and the
own-lambda column against Ep (Q2), and the loader, its lambda selection and the
h5 round trip (Q3), for a vacuum, a solvated and a multi-step run, plus the
recalculator alone on the PETase path.
"""

import json
import os
import shutil
from pathlib import Path

import numpy as np
import pytest

pytest.importorskip('openmm')

import openmm as mm
import openmm.app as mmapp
import openmm.unit as mmunit

from veloxchem.evbdriver import EvbDriver
from veloxchem.evbreporter import EvbRecalculator
from veloxchem.outputstream import OutputStream
from veloxchem.reactionsystembuilder import ReactionSystemBuilder

from test_evb_helper import EvbTestHelper
from test_evb_multistep import MultiStepFixtures

LAMBDAS = [0.0, 0.4, 0.7, 1.0]
PATH_LAMBDAS = [0.0, 0.5, 1.0]
FRAMES_PER_WINDOW = 10


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
        "n_replicas": 2,
        "trajectory_format": "dcd",
        "report_integration_energies": True,
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


def _run(folder, configuration, lambdas=LAMBDAS):

    def run():
        EVB = _seeded_driver()
        EVB.build_systems(configurations=[configuration], Lambda=lambdas)
        EVB.run_FEP(platform=EvbTestHelper.evb_platform())
        return EVB, folder / EVB.system_confs[0]["data_folder"]

    return _in_folder(folder, run)


@pytest.fixture(scope="module")
def vacuum_run(tmp_path_factory):
    folder = tmp_path_factory.mktemp("cross_vacuum")
    EVB, data_folder = _run(folder, _config("vacuum"))
    return folder, EVB, data_folder


@pytest.fixture(scope="module")
def water_run(tmp_path_factory):
    """Solvated ethanol, 3 lambda, 1 replica. The water preset with 50 equil
    steps needs minimize_every_lambda and the 0.5 window (B4 trap)."""
    folder = tmp_path_factory.mktemp("cross_water")

    def configuration():
        conf = EvbDriver(ostream=OutputStream(None)
                        ).default_system_configurations("water_NVT")
        conf.update(_config("water", padding=1.0, n_replicas=1))
        return conf

    EVB, data_folder = _run(folder, configuration(), [0.0, 0.5, 1.0])
    return folder, EVB, data_folder


@pytest.fixture(scope="module")
def path_run(tmp_path_factory):
    """The synthetic three-state path ([ethanol, ethanol, ethene+water])."""
    folder = tmp_path_factory.mktemp("cross_path")
    states = MultiStepFixtures._synthetic_three_states()

    def run():
        EVB = MultiStepFixtures._seed_driver(states)
        EVB.build_systems(configurations=[
            _config("path",
                    n_replicas=1,
                    initial_equil_NVT_steps=0,
                    skip_initial_equil=True)
        ],
                          Lambda=PATH_LAMBDAS)
        EVB.run_FEP(platform=EvbTestHelper.evb_platform())
        return EVB, EVB.system_confs[0]

    EVB, conf = _in_folder(folder, run)
    return folder, EVB, conf


def _energies(data_folder):
    """Energies.csv as (header, rows)."""
    path = data_folder / "Energies.csv"
    header = [c.strip() for c in path.read_text().splitlines()[0].split(",")]
    return header, np.loadtxt(path, skiprows=1, delimiter=",", ndmin=2)


def _integration(data_folder):
    path = data_folder / "IntegrationEnergies.csv"
    header = [c.strip() for c in path.read_text().splitlines()[0].split(",")]
    return header, np.loadtxt(path, skiprows=1, delimiter=",", ndmin=2)


def _ep(data_folder):
    return np.loadtxt(data_folder / "Data_combined.csv",
                      skiprows=1,
                      delimiter=",",
                      ndmin=2)[:, 1]


def _own_lambda_column(header, rows):
    labels = header[3:]
    own = [labels.index(f"{lam:.3f}") for lam in rows[:, 0]]
    return rows[np.arange(len(rows)), 3 + np.array(own)]


def _relative(a, b):
    return np.abs(a - b) / np.maximum(np.abs(b), 1.0)


@pytest.mark.parametrize("run_name", ["vacuum_run", "water_run"])
def test_q1_layout(run_name, request):
    _, _, data_folder = request.getfixturevalue(run_name)
    lambdas = json.loads((data_folder / "options.json").read_text())["Lambda"]
    E_header, E_rows = _energies(data_folder)
    U_header, U_rows = _integration(data_folder)

    assert U_header == ["Lambda", "replica", "direction"
                       ] + [f"{lam:.3f}" for lam in lambdas]
    assert U_rows.shape == (len(E_rows), 3 + len(lambdas))
    for name, column in (("Lambda", 0), ("replica", 1), ("direction", 2)):
        assert np.array_equal(U_rows[:, column], E_rows[:, E_header.index(name)])
    assert np.all(np.isfinite(U_rows))


@pytest.mark.parametrize("run_name", ["vacuum_run", "water_run"])
def test_q2_endpoints_equal_energies_columns(run_name, request):
    _, EVB, data_folder = request.getfixturevalue(run_name)
    E_header, E_rows = _energies(data_folder)
    U_header, U_rows = _integration(data_folder)
    for label, name in (("0.000", "reactant integration"),
                        ("1.000", "product integration")):
        U = U_rows[:, U_header.index(label)]
        E = E_rows[:, E_header.index(name)]
        assert _relative(U, E).max() <= 1e-9


def test_q2_own_lambda_equals_ep_vacuum(vacuum_run):
    """float32 frames on a single-precision platform: ~6e-5 kJ/mol std, half
    from the frames, half from the platform (work/WP-B5/gas_ep_precision.py),
    so the bound is absolute."""
    _, _, data_folder = vacuum_run
    U_header, U_rows = _integration(data_folder)
    own = _own_lambda_column(U_header, U_rows)
    ep = _ep(data_folder)
    deviation = own - ep
    print(f"vacuum own - Ep: max |dev| {np.abs(deviation).max():.2e} kJ/mol, "
          f"max rel {_relative(own, ep).max():.2e}")
    assert np.abs(deviation).max() <= 1e-3


def test_q2_own_lambda_equals_ep_water(water_run):
    """dcd frames: B4 measured std 3e-3, max 8e-3 kJ/mol."""
    _, _, data_folder = water_run
    U_header, U_rows = _integration(data_folder)
    deviation = _own_lambda_column(U_header, U_rows) - _ep(data_folder)
    print(f"water own - Ep: std {deviation.std():.2e} max |dev| "
          f"{np.abs(deviation).max():.2e} kJ/mol")
    assert np.abs(deviation).max() <= 0.05


def test_q2_interior_columns_differ(vacuum_run):
    """Each column is its own system: no two lambda columns coincide."""
    _, _, data_folder = vacuum_run
    _, U_rows = _integration(data_folder)
    columns = U_rows[:, 3:]
    for a in range(columns.shape[1]):
        for b in range(a + 1, columns.shape[1]):
            assert np.abs(columns[:, a] - columns[:, b]).max() > 1e-3


def _load(folder, EVB):
    return _in_folder(folder,
                      lambda: EVB._load_output_from_folders(1, False, 1))


@pytest.mark.parametrize("run_name", ["vacuum_run", "water_run"])
def test_q3_loader(run_name, request):
    folder, EVB, data_folder = request.getfixturevalue(run_name)
    results = _load(folder, EVB)
    result = results["configuration_results"][EVB.system_confs[0]["name"]]
    _, U_rows = _integration(data_folder)

    assert result["U_int"].shape == (len(result["Ep"]), len(results["Lambda"]))
    assert np.array_equal(result["U_int"], U_rows[:, 3:])
    assert np.allclose(result["U_int"][:, 0], result["E1_int"], rtol=1e-12)
    assert np.allclose(result["U_int"][:, -1], result["E2_int"], rtol=1e-12)


def test_q3_loader_follows_the_kept_ladder(vacuum_run, tmp_path):
    """options.json listing fewer lambda than the files hold (the B1
    fixtures): frames and U_int columns both follow the kept ladder."""
    folder, EVB, data_folder = vacuum_run
    copy = tmp_path / data_folder.name
    shutil.copytree(data_folder, copy)
    options = json.loads((copy / "options.json").read_text())
    options["Lambda"] = [0.0, 0.7, 1.0]
    (copy / "options.json").write_text(json.dumps(options))

    specific, common = EVB._load_output_files(
        str(copy / "Energies.csv"),
        str(copy / "Data_combined.csv"),
        str(copy / "options.json"),
        integration_file=str(copy / "IntegrationEnergies.csv"),
    )
    U_header, U_rows = _integration(data_folder)
    kept = np.isin(U_rows[:, 0], [0.0, 0.7, 1.0])
    columns = [U_header.index(label) for label in ("0.000", "0.700", "1.000")]

    assert specific["U_int"].shape == (kept.sum(), 3)
    assert np.array_equal(specific["U_int"], U_rows[np.ix_(kept, columns)])
    assert np.array_equal(common["Lambda_frame"], U_rows[kept, 0])


def test_q3_mismatched_rows_are_refused(vacuum_run, tmp_path):
    _, EVB, data_folder = vacuum_run
    copy = tmp_path / data_folder.name
    shutil.copytree(data_folder, copy)
    lines = (copy / "IntegrationEnergies.csv").read_text().splitlines()
    (copy / "IntegrationEnergies.csv").write_text("\n".join(lines[:-1]) + "\n")
    with pytest.raises(Exception, match="rows"):
        EVB._load_output_files(
            str(copy / "Energies.csv"),
            str(copy / "Data_combined.csv"),
            str(copy / "options.json"),
            integration_file=str(copy / "IntegrationEnergies.csv"),
        )


def test_q3_h5_round_trip(vacuum_run, tmp_path):
    folder, EVB, _ = vacuum_run
    results = _load(folder, EVB)
    name = EVB.system_confs[0]["name"]
    _in_folder(tmp_path, lambda: EVB._save_dict_as_h5(results, "cross"))
    loaded = EvbDriver._load_dict_from_h5(tmp_path / "cross.h5")
    U = loaded["configuration_results"][name]["U_int"]
    assert U.dtype == np.float64
    assert np.array_equal(U, results["configuration_results"][name]["U_int"])


def test_q3_without_the_keyword_nothing_changes(tmp_path):
    EVB, data_folder = _run(
        tmp_path,
        _config("plain", n_replicas=1, report_integration_energies=False),
        [0.0, 1.0])
    assert not (data_folder / "IntegrationEnergies.csv").exists()
    results = _load(tmp_path, EVB)
    assert "U_int" not in results["configuration_results"]["plain"]


def test_multistep_per_step_ladder(path_run):
    folder, EVB, conf = path_run
    results = _load(folder, EVB)
    steps = results["configuration_results"][conf["name"]]["steps"]
    assert len(steps) == 2

    for index, (step_conf, step) in enumerate(zip(conf["steps"], steps)):
        data_folder = folder / step_conf["data_folder"]
        assert data_folder.name == f"step_{index:02d}"
        U_header, U_rows = _integration(data_folder)
        assert U_header[3:] == [f"{lam:.3f}" for lam in PATH_LAMBDAS]

        assert step["U_int"].shape == (len(step["Ep"]), len(PATH_LAMBDAS))
        assert np.allclose(step["U_int"][:, 0], step["E1_int"], rtol=1e-12)
        assert np.allclose(step["U_int"][:, -1], step["E2_int"], rtol=1e-12)
        own = _own_lambda_column(U_header, U_rows)
        assert np.abs(own - _ep(data_folder)).max() <= 1e-3


def test_recalculator_on_the_petase_path(tmp_path):
    """Energies only: the recalculator's columns are the step's own lambda
    systems, in ladder order, on perturbed PETase geometries."""
    states = EvbTestHelper.evb_petase_ms_states()
    forming, breaking, reactive, registry = MultiStepFixtures._summarise(states)
    builder = ReactionSystemBuilder(ostream=OutputStream(None))
    builder.water_model = "cspce"
    step_systems, pes_systems, topology, positions = (
        builder.build_multi_step_systems(
            states,
            LAMBDAS,
            {
                "name": "ms",
                "temperature": 300.0
            },
            reactive_bonds=reactive,
            bond_registry=registry,
        ))
    step = 1
    systems = {
        **step_systems[step],
        **{f"state_{k}": system for k, system in pes_systems.items()}
    }
    platform = mm.Platform.getPlatformByName("Reference")

    def simulation(system):
        return mmapp.Simulation(topology, system, mm.VerletIntegrator(0.001),
                                platform)

    recalc = EvbRecalculator(tmp_path,
                             systems,
                             topology,
                             simulation,
                             n_states=len(states),
                             step=step,
                             state_left=step,
                             integration_lambdas=LAMBDAS)
    rng = np.random.default_rng(1)
    frames = [(positions * 0.1 + rng.normal(0, 0.003, positions.shape), None)
              for _ in range(2)]
    recalc.recalc_batch([(0.4, 0, 0), (0.7, 0, 1)], frames, append=False)
    recalc.close()

    U_header, U_rows = _integration(tmp_path)
    assert U_header[3:] == [f"{lam:.3f}" for lam in LAMBDAS]
    for k, lam in enumerate(LAMBDAS):
        for i, (pos, _) in enumerate(frames):
            context = mm.Context(step_systems[step][lam],
                                 mm.VerletIntegrator(0.001), platform)
            context.setPositions(pos)
            expected = context.getState(getEnergy=True).getPotentialEnergy(
            ).value_in_unit(mmunit.kilojoule_per_mole)
            assert U_rows[i, 3 + k] == pytest.approx(expected, rel=1e-9)

    E_header, E_rows = _energies(tmp_path)
    assert np.allclose(U_rows[:, 3], E_rows[:, E_header.index("left integration")],
                       rtol=1e-12)
    assert np.allclose(U_rows[:, -1],
                       E_rows[:, E_header.index("right integration")],
                       rtol=1e-12)
