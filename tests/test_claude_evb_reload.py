"""Dev tests (WP-B3a defect B3-1): reloaded start positions keep the units
EvbFepDriver.run_replicas expects.

build_systems hands run_replicas initial positions in Angstrom (it converts
with * 0.1). load_initialisation(load_top=True) and the restart_pdb path used
to store the cif positions in nm, so a reloaded run started 10x shrunk (a
solvated one fails with 'Particle coordinate is NaN').
"""

import os
from pathlib import Path

import numpy as np
import pytest

from veloxchem.evbdriver import EvbDriver
from veloxchem.outputstream import OutputStream

from test_evb_helper import EvbTestHelper


@pytest.fixture(scope="module")
def built(tmp_path_factory):
    folder = tmp_path_factory.mktemp("reload")
    previous = Path.cwd()
    os.chdir(folder)
    try:
        ff_pair = EvbTestHelper.evb_ff_pair()
        EVB = EvbDriver(ostream=OutputStream(None))
        EVB.water_model = "cspce"
        EVB.states = [ff_pair.reactant, ff_pair.product]
        EVB.forming_bonds = [set()]
        EVB.breaking_bonds = [set()]
        EVB.build_systems(configurations=["vacuum"], Lambda=[0.0, 0.4, 1.0])
        conf = EVB.system_confs[0]
        built_positions = np.array(conf["steps"][0]["initial_positions"])
        data_folder = folder / conf["data_folder"]
    finally:
        os.chdir(previous)
    return data_folder, built_positions


def _reload(data_folder, **kwargs):
    EVB = EvbDriver(ostream=OutputStream(None))
    EVB.load_initialisation(str(data_folder), **kwargs)
    return EVB.system_confs[0]


def test_load_top_positions_are_angstrom(built):
    data_folder, built_positions = built
    conf = _reload(data_folder, load_systems=False, load_top=True)
    for c in [conf] + conf["steps"]:
        loaded = np.asarray(c["initial_positions"])
        # cif stores Angstrom with 4 decimals (some writers 3).
        assert np.allclose(loaded, built_positions, atol=2e-3)


def test_restart_pdb_positions_are_angstrom(built):
    data_folder, built_positions = built
    conf = _reload(data_folder, restart_pdb="topology.cif")
    loaded = np.asarray(conf["steps"][0]["initial_positions"])
    assert np.allclose(loaded, built_positions, atol=2e-3)
    assert conf["skip_initial_equil"] is True
