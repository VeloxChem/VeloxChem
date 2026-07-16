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

import numpy as np
import pytest

from veloxchem.evbdriver import EvbDriver
from veloxchem.errorhandler import VeloxChemError

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

    @pytest.mark.parametrize("recalc_mode", ["inline", "deferred"])
    def test_compute_force_groups(self, recalc_mode):
        # recalc_mode='deferred' used to be incompatible with in-loop force
        # groups; compute_force_groups() runs identically afterwards
        # regardless of how sampling was done, since it replays the trajectory
        # rather than hooking into the sampling loop.
        config = _tiny_config(f"tiny_fg_{recalc_mode}", recalc_mode=recalc_mode)

        with evb_chdir_tmp():
            EVB, data_folder = _run_tiny_fep(evb_ff_pair(), config)

            EVB.compute_force_groups()

            rea_file = data_folder / "ForceGroups_rea.csv"
            pro_file = data_folder / "ForceGroups_pro.csv"
            assert rea_file.exists(), f"missing {rea_file}"
            assert pro_file.exists(), f"missing {pro_file}"

            energies = data_folder / "Energies.csv"
            n_energy_rows = _finite_energy_rows(energies)
            n_rea_rows = _finite_energy_rows(rea_file)
            n_pro_rows = _finite_energy_rows(pro_file)
            assert n_rea_rows == n_energy_rows
            assert n_pro_rows == n_energy_rows

            from veloxchem.reactionsystembuilder import EvbForceGroup
            expected_header = EvbForceGroup.get_header().strip()
            with open(rea_file) as handle:
                assert handle.readline().strip() == expected_header
            with open(pro_file) as handle:
                assert handle.readline().strip() == expected_header

    def test_compute_force_groups_bonded_decomp_requires_decomp_systems(self):
        # The tiny fixture builds systems with decompose_bonded left at its
        # default (False), so reactant_bonded_decomp/product_bonded_decomp
        # are absent; bonded_decomp=True should fail with a clear error
        # rather than a raw KeyError.
        config = _tiny_config("tiny_fg_bonded_guard")

        with evb_chdir_tmp():
            EVB, _ = _run_tiny_fep(evb_ff_pair(), config)

            with pytest.raises(VeloxChemError, match="bonded_decomp"):
                EVB.compute_force_groups(bonded_decomp=True)
