"""Cheap-sampling tests: run a tiny but real FEP end to end.

A microscopic vacuum FEP on the 9-atom ethanol -> ethene + water reaction (few
lambda windows, a handful of sampling steps) exercises the FEP driver and the
reporters. Assertions are structural / sanity-band (files exist, expected row
counts, finite energies) rather than golden values, because the Langevin
integrator is not reproducible bit-for-bit across platforms.
"""

import csv
import math
from pathlib import Path

import pytest

from veloxchem.evbdriver import EvbDriver

from veloxchem.outputstream import OutputStream

from test_evb_helper import evb_chdir_tmp, evb_ff_pair

pytestmark = [pytest.mark.timeconsuming]

pytest.importorskip('openmm')


def _tiny_config(name, **overrides):
    config = {
        "name": name,
        "temperature": 300.0,
        "sample_steps": 20,
        "write_step": 10,
        "equil_NVT_steps": 2,
        "equil_NPT_steps": 0,
        "initial_equil_NVT_steps": 0,
        "initial_equil_NPT_steps": 0,
        "skip_initial_equil": True,
        "step_size": 0.001,
        "n_replicas": 1,
    }
    config.update(overrides)
    return config


def _seed_driver(ff_pair):
    EVB = EvbDriver(ostream=OutputStream(None))
    EVB.water_model = "cspce"
    EVB.reactant = ff_pair.reactant
    EVB.product = ff_pair.product
    EVB.forming_bonds = set()
    EVB.breaking_bonds = set()
    return EVB


# NOTE: use float lambda endpoints. The FEP driver formats lambda as f"{l:.3}",
# which raises on integer 0/1 (evbfepdriver.py). Production always passes
# np.linspace floats, so this is a latent bug rather than something under test.
def _run_tiny_fep(ff_pair, config, lambdas=(0.0, 0.5, 1.0)):
    EVB = _seed_driver(ff_pair)
    EVB.build_systems(configurations=[config], Lambda=list(lambdas))
    EVB.run_FEP()
    return EVB, Path(EVB.system_confs[0]["data_folder"])


def _finite_energy_rows(csv_path):
    """Return the number of data rows and assert every numeric cell is finite."""
    with open(csv_path, newline="") as handle:
        reader = csv.reader(handle)
        rows = [row for row in reader if row]
    # header + at least one data row
    n_data = 0
    for row in rows[1:]:
        n_data += 1
        for cell in row:
            try:
                val = float(cell)
            except ValueError:
                continue
            assert math.isfinite(val), f"non-finite value {cell} in {csv_path}"
    return n_data


class TestTinyFep:

    @pytest.mark.parametrize("recalc_mode", ["inline", "deferred"])
    def test_fep_smoke(self, recalc_mode):
        config = _tiny_config(f"tiny_{recalc_mode}", recalc_mode=recalc_mode)

        with evb_chdir_tmp():
            _, data_folder = _run_tiny_fep(evb_ff_pair(), config)

            energies = data_folder / "Energies.csv"
            combined = data_folder / "Data_combined.csv"
            assert energies.exists(), f"missing {energies}"
            assert combined.exists(), f"missing {combined}"

            assert _finite_energy_rows(energies) > 0
            assert _finite_energy_rows(combined) > 0
