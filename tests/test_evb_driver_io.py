"""Smoke tests for EvbDriver orchestration and file I/O (no MD, no QM).

Covers default configuration presets, the h5 save/load round-trip, options.json
creation/merging, and the build_systems -> load_initialisation folder round-trip
(seeded from committed force fields so no QM runs).
"""

import json
from pathlib import Path

import numpy as np
import pytest

from veloxchem.evbdriver import EvbDriver

from veloxchem.outputstream import OutputStream

from test_evb_helper import EvbTestHelper

pytestmark = [pytest.mark.timeconsuming]

pytest.importorskip('openmm')


class TestDefaultConfigurations:

    @pytest.mark.parametrize("name", [
        "vacuum",
        "vacuum_NVT",
        "vacuum_NVE",
        "water",
        "water_NVT",
        "water_NPT",
        "E_field",
        "no_reactant",
        "ts_guesser",
        "debug",
    ])
    def test_known_config_has_name(self, name):
        EVB = EvbDriver(ostream=OutputStream(None))
        conf = EVB.default_system_configurations(name)
        assert isinstance(conf, dict)
        assert "name" in conf

    @pytest.mark.parametrize("name", ["vacuum", "water", "water_NPT", "debug"])
    def test_temperature_configs(self, name):
        EVB = EvbDriver(ostream=OutputStream(None))
        conf = EVB.default_system_configurations(name)
        assert "temperature" in conf

    def test_unknown_config_raises(self):
        EVB = EvbDriver(ostream=OutputStream(None))
        with pytest.raises(ValueError):
            EVB.default_system_configurations("definitely_not_a_solvent_xyz")


class TestH5RoundTrip:

    def test_save_load_roundtrip(self):
        EVB = EvbDriver(ostream=OutputStream(None))
        data = {
            "scalar_int": 3,
            "scalar_float": 2.5,
            "flag": True,
            "vector": np.array([1.0, 2.0, 3.0]),
            "nested": {
                "matrix": np.arange(6.0).reshape(2, 3),
                "label": "abc",
            },
        }
        with EvbTestHelper.evb_chdir_tmp():
            EVB._save_dict_as_h5(data, "roundtrip")
            assert Path("roundtrip.h5").exists()

            loaded = EVB._load_dict_from_h5(Path("roundtrip.h5"))

        assert int(loaded["scalar_int"]) == 3
        assert abs(float(loaded["scalar_float"]) - 2.5) < 1e-12
        assert bool(loaded["flag"]) is True
        assert np.allclose(loaded["vector"], [1.0, 2.0, 3.0])
        assert np.allclose(loaded["nested"]["matrix"],
                           np.arange(6.0).reshape(2, 3))


class TestOptionsJson:

    def test_create_and_merge(self):
        EVB = EvbDriver(ostream=OutputStream(None))
        with EvbTestHelper.evb_chdir_tmp():
            folder = Path("cfg_data")
            folder.mkdir()
            conf = {"data_folder": str(folder)}

            EVB.update_options_json({"a": 1, "b": 2}, conf)
            options_path = folder / "options.json"
            assert options_path.exists()

            EVB.update_options_json({"b": 3, "c": 4}, conf)
            with options_path.open() as handle:
                merged = json.load(handle)

        assert merged == {"a": 1, "b": 3, "c": 4}


class TestBuildAndLoad:
    """build_systems writes a data folder that load_initialisation reads back."""

    def _seed_driver(self, ff_pair):
        EVB = EvbDriver(ostream=OutputStream(None))
        EVB.water_model = "cspce"
        EVB.states = [ff_pair.reactant, ff_pair.product]
        EVB.forming_bonds = [set()]
        EVB.breaking_bonds = [set()]
        return EVB

    def test_build_then_load(self):
        with EvbTestHelper.evb_chdir_tmp():
            EVB = self._seed_driver(EvbTestHelper.evb_ff_pair())
            EVB.build_systems(configurations=["vacuum"], Lambda=[0, 0.4, 1])

            conf = EVB.system_confs[0]
            data_folder = Path(conf["data_folder"])
            assert data_folder.is_dir()
            assert (data_folder / "options.json").exists()
            assert (data_folder / "topology.cif").exists()
            assert (data_folder / "reactant_ff_data.json").exists()
            assert (data_folder / "run").is_dir()

            # A fresh driver should be able to load the written folder.
            EVB2 = EvbDriver(ostream=OutputStream(None))
            EVB2.load_initialisation(str(data_folder), "reloaded")

        assert len(EVB2.system_confs) == 1
        assert EVB2.lambda_vec == [0, 0.4, 1]
        assert EVB2.temperature == conf["temperature"]
