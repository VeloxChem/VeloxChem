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

    def test_compute_force_groups_bonded_decomp(self):
        # Test if compute_force_groups(decompose_bonded=True) regenerates the bonded
        # decomposition systems and writes the expected files.
        config = _tiny_config("tiny_fg_bonded_regen")

        with evb_chdir_tmp():
            EVB, data_folder = _run_tiny_fep(evb_ff_pair(), config)
            run_folder = Path(EVB.system_confs[0]["run_folder"])

            EVB.compute_force_groups(decompose_bonded=True)

            for fname in ("bonded_E1_decomp.csv", "bonded_E2_decomp.csv",
                          "bonded_params.csv"):
                assert (data_folder / fname).exists(), f"missing {fname}"

            # Regenerated systems are persisted under unique names so a later
            # load_initialisation(load_systems=True) picks them up.
            for name in ("reactant_bonded_decomp", "product_bonded_decomp"):
                assert (run_folder / f"{name}_sys.xml").exists(), \
                    f"missing regenerated {name}_sys.xml"

    def test_reconstruct_decomp_atom_indices_matches_builder(self):
        # The load-bearing guarantee for post-hoc NB decomposition: the
        # reaction/solvent atom partition reconstructed from the (reloaded)
        # topology must reproduce exactly what the builder derived at build
        # time. No FEP/MD needed - just a solvated build.
        from veloxchem.reactionsystembuilder import ReactionSystemBuilder
        from veloxchem.evbfepdriver import EvbFepDriver

        config = {
            "name": "water_recon",
            "temperature": 300.0,
            "solvent": "cspce",
            "pressure": 1,
            "padding": 0.5,
            "ion_count": 0,
            "neutralize": False,
        }

        with evb_chdir_tmp():
            pair = evb_ff_pair()
            builder = ReactionSystemBuilder(ostream=OutputStream(None))
            builder.water_model = "cspce"
            _, topology, _ = builder.build_systems(pair.reactant, pair.product,
                                                   [0.0, 0.5, 1.0], config)

            FEP = EvbFepDriver(ostream=OutputStream(None))
            rec_idx, solv_idx = FEP._reconstruct_decomp_atom_indices(topology)

            assert rec_idx == sorted(builder.reaction_atom_indices)
            assert solv_idx == sorted(builder.solvent_atom_ids)
            # sanity: solute + solvent partition the whole topology
            assert len(rec_idx) + len(solv_idx) == topology.getNumAtoms()
            assert len(solv_idx) > 0

    # TODO: base this on a pre-generated trajectory
    # @pytest.mark.timeconsuming
    # def test_compute_force_groups_nb_decomp(self):
    #     # End-to-end nonbonded decomposition against a solvated configuration:
    #     # decompose_nb reconstructs the reaction/solvent partition from the
    #     # topology, regenerates the decomp_* systems, saves them to run_folder,
    #     # and writes NB_decompositions.csv. The trajectory is fabricated from the
    #     # system's own (stable, packed) initial positions rather than sampled,
    #     # so the test exercises the decomposition plumbing without depending on
    #     # a freshly-solvated box equilibrating without blowing up.
    #     import openmm as mm
    #     from MDAnalysis.lib.formats.libmdaxdr import XTCFile

    #     config = {
    #         "name": "water_nb",
    #         "temperature": 300.0,
    #         "solvent": "cspce",
    #         "pressure": 1,
    #         "padding": 1.0,
    #         "ion_count": 0,
    #         "neutralize": False,
    #         "step_size": 0.001,
    #     }

    #     with evb_chdir_tmp():
    #         EVB = _seed_driver(evb_ff_pair())
    #         EVB.build_systems(configurations=[config], Lambda=[0.0, 0.5, 1.0])
    #         conf = EVB.system_confs[0]
    #         data_folder = Path(conf["data_folder"])
    #         run_folder = Path(conf["run_folder"])

    #         # Fabricate a single-frame trajectory.xtc (nm) from the packed
    #         # initial positions, plus a row-aligned Energies.csv.
    #         positions = conf["initial_positions"]
    #         try:
    #             pos_nm = np.asarray(positions.value_in_unit(mm.unit.nanometer),
    #                                 dtype=np.float32)
    #         except AttributeError:
    #             pos_nm = np.asarray(positions, dtype=np.float32)
    #         box_nm = np.asarray(
    #             conf["topology"].getPeriodicBoxVectors().value_in_unit(
    #                 mm.unit.nanometer),
    #             dtype=np.float32)
    #         with XTCFile(str(data_folder / "trajectory.xtc"), mode='w') as xtc:
    #             xtc.write(pos_nm, box_nm, 1, 0.0, 1000.0)
    #         with open(data_folder / "Energies.csv", "w") as handle:
    #             handle.write("Lambda, reactant PES, product PES, reactant "
    #                          "integration, product integration, Em, replica, "
    #                          "direction\n")
    #             handle.write("0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0, 0\n")

    #         EVB.compute_force_groups(decompose_nb=[[0]])

    #         nb_file = data_folder / "NB_decompositions.csv"
    #         assert nb_file.exists(), f"missing {nb_file}"
    #         assert _finite_energy_rows(nb_file) == 1

    #         # Reactant and product halves must be symmetric so the midpoint
    #         # split in _load_output_files lines columns up.
    #         with open(nb_file) as handle:
    #             names = handle.readline().strip().split(",")
    #         rea_cols = [n for n in names if n.strip().startswith("decomp_rea_")]
    #         pro_cols = [n for n in names if n.strip().startswith("decomp_pro_")]
    #         assert len(rea_cols) == len(pro_cols) > 0
    #         assert names[:len(names) // 2] == rea_cols

    #         # The regenerated nonbonded decomposition systems are persisted for
    #         # both states, so a later reload picks them up.
    #         assert list(run_folder.glob("decomp_rea_*_sys.xml"))
    #         assert list(run_folder.glob("decomp_pro_*_sys.xml"))
