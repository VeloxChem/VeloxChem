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
"""Nonbonded curve mechanics; synthetic recovery is not physical validation."""
from copy import deepcopy

import numpy as np
import pytest
import veloxchem as vlx

from veloxchem.autocdam import Calculator, MetalComplexLJDriver, QMScanDataset
from veloxchem.autocdam.models import read_h5_group, write_h5_groups


@pytest.fixture
def baseline(tmp_path):
    # Synthetic Zn–O environment with a bonded term and fixed nonzero charges.
    molecule = vlx.Molecule.read_xyz_string("2\nsynthetic\nZn 0 0 0\nO 2 0 0\n")
    molecule.set_charge(0)
    molecule.set_multiplicity(1)
    ff = vlx.MMForceFieldGenerator(
        ostream=vlx.OutputStream(str(tmp_path / "baseline.log")))
    ff.molecule = molecule
    ff.atoms = {
        0: {
            "type": "Zn",
            "name": "Zn1",
            "mass": 65.38,
            "charge": 0.5,
            "sigma": 0.22,
            "epsilon": 0.8,
            "comment": "synthetic baseline"
        },
        1: {
            "type": "o",
            "name": "O1",
            "mass": 16.0,
            "charge": -0.5,
            "sigma": 0.31,
            "epsilon": 0.4,
            "comment": "synthetic baseline"
        }
    }
    ff.bonds = {(0, 1): {"equilibrium": 0.2, "force_constant": 1000.0}}
    ff.angles, ff.dihedrals, ff.impropers, ff.pairs = {}, {}, {}, {}
    ff.partial_charges = np.array([0.5, -0.5])
    return ff


def make_driver(baseline, tmp_path, **overrides):
    study = dict(molecule=baseline.molecule,
                 target_atom_index=0,
                 baseline_force_field=baseline)
    settings = dict(probe="He",
                    probe_parameters={
                        "sigma_angstrom": 2.0,
                        "epsilon_kjmol": 0.3,
                        "source": "synthetic test probe"
                    },
                    scan_distances_angstrom=np.linspace(2.0, 5.0, 14),
                    scan_directions=[[0, 0, 1], [0, 1, 0]],
                    clash_cutoff_scale=0,
                    sigma_bounds_angstrom=[1.0, 4.0],
                    epsilon_bounds_kjmol=[0.01, 5.0],
                    energy_tolerance_kjmol=1e-5,
                    filename=str(tmp_path / "run"))
    for key, value in overrides.items():
        (study if key in study else settings)[key] = value
    driver = MetalComplexLJDriver(
        ostream=vlx.OutputStream(str(tmp_path / "driver.log")))
    driver.update_settings(settings)
    driver.setup(**study)
    return driver


def analytical_interaction(driver, coordinates, sigma=2.6, epsilon=1.4):
    """Independent all-host pair sum for a neutral disconnected probe."""
    total = 0.0
    probe = driver.study.probe_parameters
    for i, atom in driver.study.baseline_force_field.atoms.items():
        r = np.linalg.norm(np.array(coordinates[-1]) - coordinates[i])
        host_sigma = sigma if i == driver.study.target_atom_index else atom[
            "sigma"] * 10
        host_epsilon = epsilon if i == driver.study.target_atom_index else atom[
            "epsilon"]
        mixed_sigma = (host_sigma + probe["sigma_angstrom"]) / 2
        mixed_epsilon = np.sqrt(host_epsilon * probe["epsilon_kjmol"])
        total += 4 * mixed_epsilon * ((mixed_sigma / r)**12 -
                                      (mixed_sigma / r)**6)
    return total


def synthesize_qm(driver, fail_one=False):
    calls = []

    class SyntheticQM:

        def single_point(self, molecule, *, settings):
            calls.append(molecule.number_of_atoms())
            if molecule.number_of_atoms() == 2:
                return -100.0
            if molecule.number_of_atoms() == 1:
                return -2.0
            coords = molecule.get_coordinates_in_angstrom().tolist()
            if fail_one and len(calls) == 3:
                raise RuntimeError("SCF did not converge")
            return -102.0 + analytical_interaction(
                driver, coords) / vlx.hartree_in_kjpermol()

    driver.calculator = SyntheticQM()
    driver.compute_qm_dataset()
    return calls


class TestAutoCdamMetalComplexLj:

    def test_recovery_preserves_surroundings_and_exports(
            self, baseline, tmp_path):
        original = deepcopy(baseline.atoms)
        driver = make_driver(baseline, tmp_path)
        calls = synthesize_qm(driver)
        result = driver.fit()
        assert calls.count(1) == calls.count(2) == 1
        assert result.converged
        assert result.parameters["metal_lj_sigma"] == pytest.approx(2.6,
                                                                    abs=1e-5)
        assert result.parameters["metal_lj_epsilon"] == pytest.approx(1.4,
                                                                      abs=1e-5)
        assert result.metrics["rmse_kjmol"] < 1e-6
        assert baseline.atoms == original
        assert driver.baseline_force_field.atoms == original
        fitted = driver.fitted_force_field()
        assert fitted.atoms[1] == original[1]
        assert fitted.atoms[0]["charge"] == original[0]["charge"]
        assert fitted.bonds == baseline.bonds
        assert driver._combined_mm.force_field["bonds"] == baseline.bonds
        assert len(driver._combined_mm.force_field["atoms"]) == 3
        engine = vlx.MMDriver(ostream=driver.ostream)
        engine.load_force_field(fitted)
        assert engine.force_field["atoms"][0]["sigma"] == pytest.approx(0.26)
        path = driver.export_fitted_model(tmp_path / "fitted_complex.xml")
        import xml.etree.ElementTree as ET
        root = ET.parse(path).getroot()
        atoms = root.findall("NonbondedForce/Atom")
        assert float(atoms[0].attrib["sigma"]) == pytest.approx(0.26)
        assert float(atoms[1].attrib["sigma"]) == original[1]["sigma"]
        assert path.with_suffix(".fit.json").exists()
        saved = read_h5_group(tmp_path / "run", "metal_complex_lj/result")
        assert saved["diagnostics"]["scientific_validation"] == "not performed"
        assert len(saved["diagnostics"]["qm_interaction_kjmol"]) == len(
            saved["diagnostics"]["mm_interaction_kjmol"])

    def test_mm_matches_independent_pair_sum_and_qm_stays_fixed(
            self, baseline, tmp_path):
        driver = make_driver(baseline, tmp_path)
        synthesize_qm(driver)
        raw = [row.qm_combined_energy_hartree for row in driver.dataset.records]
        result = driver.evaluate_mm({
            "metal_lj_sigma": 2.6,
            "metal_lj_epsilon": 1.4
        })
        expected = [
            analytical_interaction(driver, row.coordinates_angstrom)
            for row in driver.dataset.records
        ]
        np.testing.assert_allclose(result.values["MM_INTERACTION_ENERGY"],
                                   expected,
                                   atol=1e-10)
        # Cross-term evaluation must equal full native-MM subtraction, even when
        # the host has bonded and charged internal interactions.
        driver._host_mm.compute(driver.study.molecule)
        host_energy = driver._host_mm.get_energy()
        driver._combined_mm.compute(
            driver._combined_molecule(driver.dataset.records[0]))
        subtracted = (driver._combined_mm.get_energy() -
                      host_energy) * vlx.hartree_in_kjpermol()
        assert subtracted == pytest.approx(expected[0], abs=1e-10)
        changed = driver.evaluate_mm({
            "metal_lj_sigma": 2.2,
            "metal_lj_epsilon": 0.8
        })
        assert not np.allclose(changed.values["MM_INTERACTION_ENERGY"],
                               expected)
        assert raw == [
            row.qm_combined_energy_hartree for row in driver.dataset.records
        ]

    def test_portable_dataset_failures_reweighting_and_main_driver(
            self, baseline, tmp_path):
        driver = make_driver(baseline, tmp_path)
        synthesize_qm(driver, fail_one=True)
        assert not driver.dataset.records[0].usable
        assert "SCF did not converge" in driver.dataset.records[0].reject_reason
        dataset_path = driver.save_dataset(tmp_path / "portable")
        saved_json = read_h5_group(dataset_path, "metal_complex_lj/qm_scan")
        assert "fingerprint" not in saved_json
        weights = {driver.dataset.records[1].record_id: 3.0}

        class NoQM:

            def single_point(self, *args, **kwargs):
                raise AssertionError("loaded references must not rerun QM")

        main = make_driver(baseline, tmp_path, point_weights=weights)
        main.calculator = NoQM()
        main.load_dataset(dataset_path)
        result = main.fit()
        assert result.diagnostics["n_rejected"] == 1
        assert result.converged
        assert main.evaluate_mm(
            result.parameters
        ).raw_data["parameters"] == result.parameters.to_dict()
        assert main.evaluate().objective < 1
        assert main.export_fitted_model(tmp_path / "main.xml").exists()
        mismatch = make_driver(baseline,
                               tmp_path,
                               qm_settings={"basis": "STO-3G"})
        assert mismatch.load_dataset(
            dataset_path).records == driver.dataset.records
        # Metadata is advisory; callers own correctness.
        saved_json["metadata"] = {}
        write_h5_groups(dataset_path, {"metal_complex_lj/qm_scan": saved_json})
        loaded = QMScanDataset.load(dataset_path)
        assert loaded.records == driver.dataset.records
        assert loaded.metadata == {}

    def test_auto_and_clash_scans_leave_input_fixed(self, baseline, tmp_path):
        driver = make_driver(baseline,
                             tmp_path,
                             scan_directions="auto",
                             clash_cutoff_scale=0.85)
        coords = baseline.molecule.get_coordinates_in_angstrom().copy()
        rows = driver.generate_scans()
        np.testing.assert_allclose(rows[0].coordinates_angstrom[:-1], coords)
        assert np.isclose(np.linalg.norm(rows[0].direction), 1)
        assert rows[0].direction[0] < 0
        clash = make_driver(baseline,
                            tmp_path,
                            scan_distances_angstrom=[2.0, 3.0],
                            scan_directions=[[1, 0, 0]],
                            clash_cutoff_scale=0.85)
        assert not clash.generate_scans()[0].accepted
        np.testing.assert_allclose(
            baseline.molecule.get_coordinates_in_angstrom(), coords)

    def test_generator_default_calls_and_probe_defaults(self, baseline,
                                                        tmp_path):
        driver = make_driver(baseline,
                             tmp_path,
                             baseline_force_field=None,
                             probe_parameters=None)
        calls = []

        class Generator:

            def create_topology(self, molecule, **kwargs):
                calls.append(kwargs)
                self.molecule = molecule
                if molecule.number_of_atoms() == 2:
                    self.__dict__.update(baseline.__dict__)
                else:
                    self.atoms = {
                        0: {
                            "sigma": 0.2,
                            "epsilon": 0.3,
                            "charge": 0,
                            "comment": "probe default"
                        }
                    }

        driver.force_field_factory = Generator
        driver.prepare_baseline()
        assert calls == [{}, {"resp": False}]
        assert driver.probe_atom["comment"] == "probe default"
        driver.prepare_baseline()
        assert len(calls) == 2

    @pytest.mark.parametrize("overrides", [
        {
            "target_atom_index": -1
        },
        {
            "probe": "Ar"
        },
        {
            "scan_distances_angstrom": [0]
        },
        {
            "scan_distances_angstrom": [2, 2]
        },
        {
            "scan_directions": [[0, 0, 0]]
        },
        {
            "qm_settings": {
                "typo": 1
            }
        },
    ])
    def test_invalid_study_inputs(self, baseline, tmp_path, overrides):
        with pytest.raises(ValueError):
            make_driver(baseline, tmp_path, **overrides)

    def test_zero_weights_missing_data_and_study_mutation(
            self, baseline, tmp_path):
        driver = make_driver(baseline, tmp_path)
        with pytest.raises(RuntimeError, match="no accepted"):
            driver.evaluate_mm()
        synthesize_qm(driver)
        driver.point_weights = {
            row.record_id: 0 for row in driver.dataset.records
        }
        with pytest.raises(ValueError, match="positive weight"):
            driver.fit()
        driver.qm_settings = {"basis": "STO-3G"}
        with pytest.raises(ValueError, match="changed"):
            driver.generate_scans()

    def test_real_single_point_settings_cache(self, tmp_path):
        stream = vlx.OutputStream(str(tmp_path / "scf.log"))
        try:
            calculator = Calculator(ostream=stream)
            mol = vlx.Molecule.read_xyz_string("1\n\nHe 0 0 0\n")
            metadata = {}
            first = calculator.single_point(mol,
                                            settings={"basis": "STO-3G"},
                                            metadata=metadata)
            assert np.isfinite(first) and metadata[
                "scf_converged"] and not metadata["cache_hit"]
            calculator.single_point(mol,
                                    settings={"basis": "STO-3G"},
                                    metadata=metadata)
            assert metadata["cache_hit"]
            calculator.single_point(mol,
                                    settings={
                                        "basis": "STO-3G",
                                        "conv_thresh": 1e-5
                                    },
                                    metadata=metadata)
            assert not metadata["cache_hit"]
            with pytest.raises(ValueError, match="unsupported"):
                calculator.single_point(mol, settings={"typo": 1})
        finally:
            stream.close()

    @pytest.mark.parametrize("probe", ["He", "Ne"])
    def test_real_probe_defaults_and_combined_qm(self, baseline, tmp_path,
                                                 probe):
        # Tiny native-QM geometry catches basis-label assembly and He/Ne support.
        driver = make_driver(baseline,
                             tmp_path,
                             probe=probe,
                             probe_parameters=None,
                             scan_distances_angstrom=[5.0],
                             scan_directions=[[0, 0, 1]],
                             qm_settings={
                                 "basis": "STO-3G",
                                 "max_iter": 150
                             })
        # Use the generated host's exact input geometry/state but skip a metal SCF
        # here: a one-He host suffices to verify construction of the probe basis.
        helium = vlx.Molecule.read_xyz_string("1\n\nHe 0 0 0\n")
        qm_driver = MetalComplexLJDriver(ostream=driver.ostream)
        qm_driver.update_settings({
            "probe": probe,
            "scan_distances_angstrom": [5],
            "scan_directions": [[0, 0, 1]],
            "qm_settings": {"basis": "STO-3G"},
            "filename": str(tmp_path / probe)
        })
        dataset = qm_driver.compute_qm_dataset(helium, 0)
        assert dataset.records[0].usable
        driver.prepare_baseline()
        assert driver.probe_atom["charge"] == 0
        assert driver.probe_atom["sigma"] > 0 and driver.probe_atom[
            "epsilon"] > 0
        assert driver.probe_atom["comment"] in {"UFF", "GAFF", "TM"}

    def test_baseline_metadata_is_advisory_on_reload(self, baseline, tmp_path):
        original = make_driver(baseline, tmp_path)
        synthesize_qm(original)
        original.prepare_baseline()
        saved = original.save_dataset(tmp_path / "with_baseline")
        assert "fingerprint" not in read_h5_group(
            saved, "metal_complex_lj/qm_scan")["metadata"]["baseline"]
        modified = deepcopy(baseline.atoms)
        modified[1]["epsilon"] *= 2
        other_ff = deepcopy(baseline.atoms)
        baseline.atoms = modified
        other = make_driver(baseline, tmp_path)
        other.prepare_baseline()
        assert other.load_dataset(saved).records == original.dataset.records
        other.save_dataset(saved)
        assert "fingerprint" not in read_h5_group(
            saved, "metal_complex_lj/qm_scan")["metadata"]["baseline"]
        baseline.atoms = other_ff

    def test_qm_reference_failure_is_saved(self, baseline, tmp_path):
        driver = make_driver(baseline, tmp_path)

        class BrokenQM:

            def single_point(self, *args, **kwargs):
                raise RuntimeError("reference failed")

        driver.calculator = BrokenQM()
        with pytest.raises(RuntimeError, match="dataset saved"):
            driver.compute_qm_dataset()
        data = QMScanDataset.load(tmp_path / "run")
        assert data.metadata["reference_failure"] == "reference failed"
        assert not any(row.usable for row in data.records)

    def test_de_fit_reports_both_optimizer_stages(self, baseline, tmp_path):
        driver = make_driver(baseline, tmp_path)
        synthesize_qm(driver)
        driver.update_settings({
            "strategy": "differential_evolution+least_squares",
            "de_maxiter": 2,
            "de_popsize": 4
        })
        result = driver.fit()
        assert result.parameters["metal_lj_sigma"] == pytest.approx(2.6,
                                                                    abs=1e-5)
        assert result.diagnostics["de_success"] is False
        assert result.diagnostics["optimizer_success"] is True
        assert result.converged

    @pytest.mark.parametrize(
        "xyz, accepted",
        [
            ("Zn 0 0 0\nO 1.9 0.8 0.3\n",
             True),  # not exact in binary; bohr round trip
            ("Zn 0 0 0\nO 1.9 0.81 0.3\n",
             False),  # a genuinely different geometry
        ])
    def test_supplied_baseline_matches_within_coordinate_tolerance(
            self, baseline, tmp_path, xyz, accepted):
        study_molecule = vlx.Molecule.read_xyz_string(
            "2\nstudy\nZn 0 0 0\nO 1.9 0.8 0.3\n")
        supplied = vlx.Molecule.read_xyz_string(f"2\nsupplied\n{xyz}")
        for molecule in (study_molecule, supplied):
            molecule.set_charge(0)
            molecule.set_multiplicity(1)
        baseline.molecule = supplied
        driver = make_driver(baseline, tmp_path, molecule=study_molecule)
        if accepted:
            driver.prepare_baseline()
            assert driver.baseline_force_field.atoms[0][
                "sigma"] == pytest.approx(0.22)
        else:
            with pytest.raises(ValueError,
                               match="does not match the exact input molecule"):
                driver.prepare_baseline()
