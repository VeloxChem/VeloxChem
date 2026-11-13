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

from veloxchem.autocdam import AutoCDAM, CDAMCandidate


class TestAutoCdamDriverSmoke:

    def test_autocdam_physical_fit_runs_end_to_end(self, tmp_path):
        import math
        import random
        from veloxchem import OutputStream
        from veloxchem.autocdam import FitResult
        from veloxchem.autocdam.models import read_h5_group

        random.seed(42)
        ostream = OutputStream(str(tmp_path / "fit.log"))
        try:
            cdam = AutoCDAM(ostream=ostream)
            cdam.update_settings({
                "search_method": "staged",
                "max_iterations": 1,
                "samples_per_stage": 2,
                "pilot": True,
                "work_dir": str(tmp_path),
                "filename": str(tmp_path / "run"),
            })
            cdam.fit(CDAMCandidate("Zn", 2), {
                "IOD": {
                    "value": 2.09,
                    "tolerance": 0.03
                },
                "CN": {
                    "value": 6.0,
                    "tolerance": 0.25
                }
            })
            result = cdam.fitting_result()
            assert result.iterations == 1
            assert result.diagnostics["n_calculations"] == 2
            assert math.isfinite(result.objective)
            assert set(result.diagnostics["best_properties"]) == {"IOD", "CN"}
            best = next(
                sample for sample in cdam.fitting_driver.table.samples
                if sample.sample_id == result.diagnostics["best_sample_id"])
            assert best.parameters.values == result.parameters.values
            assert best.calculation.values == result.diagnostics[
                "best_properties"]
            assert best.calculation.raw_data["md_analysis"]["n_frames"] >= 2
            assert best.calculation.metadata["IOD"][
                "method"] == "OpenMM MD trajectory analysis"
            assert best.calculation.metadata["IOD"]["md_options"][
                "timestep_fs"] == 0.25
            assert best.calculation.metadata["IOD"]["pilot"] is True
            assert best.calculation.metadata["IOD"]["fast_mode"] is False
            assert read_h5_group(tmp_path / "run.h5",
                                 "autocdam/checkpoint") is not None
            assert FitResult.load(tmp_path / "run").parameters.to_dict(
            ) == result.parameters.to_dict()
        finally:
            ostream.close()

    def test_custom_evaluator_receives_unmodified_calculator_result(
            self) -> None:
        from veloxchem.autocdam import CalculationResult, EvaluationResult, ParameterSet

        calculation = CalculationResult(
            values={"IOD": 2.1},
            metadata={"IOD": {
                "uncertainty": 0.05
            }},
            raw_data={"trajectory": "caller-owned"},
        )

        class SuppliedCalculator:

            def calculate(self, batch, parameters, targets):
                assert isinstance(parameters, ParameterSet)
                return calculation

        class SuppliedEvaluator:

            def evaluate(self, result, targets):
                assert result is calculation
                assert targets[0].name == "IOD"
                return EvaluationResult(objective=7.0,
                                        metrics={"custom_loss": 7.0})

        cdam = AutoCDAM()
        cdam.calculator = SuppliedCalculator()
        cdam.evaluator = SuppliedEvaluator()
        cdam.fit(
            CDAMCandidate("Zn", 2, metal_lj_sigma=2.0, metal_lj_epsilon=0.2),
            {"IOD": {
                "reference": 2.0,
                "tolerance": 0.1
            }},
            fit_parameters={
                "dummy_charge": {
                    "initial": 0.5,
                    "bounds": [0.0, 1.0]
                }
            })
        result = cdam.fitting_result()
        assert result.objective == 7.0
        assert result.metrics == {"custom_loss": 7.0}


class _ToyCalculator:
    """Analytic IOD and CN of a CDAM proposal; no MD."""

    def calculate(self, batch, parameters, targets, **options):
        from veloxchem.autocdam import CalculationResult

        sigma = batch.system.nonbonded_parameters["metal_lj_sigma_angstrom"]
        return CalculationResult(values={"IOD": 1.30 + 0.30 * sigma, "CN": 6.0},
                                 metadata={"IOD": {"unit": "angstrom"},
                                           "CN": {"unit": "count"}})


class TestAutoCdamAttributeInterface:

    def test_fit_from_molecule_target_and_write_openmm_files(self, tmp_path):
        import pytest
        from veloxchem import Molecule, OutputStream
        from veloxchem.errorhandler import VeloxChemError

        mol = Molecule.read_smiles("[Al]")
        mol.set_charge(3)

        cdamfit = AutoCDAM(ostream=OutputStream(None))
        with pytest.raises(VeloxChemError, match="before write_openmm_files"):
            cdamfit.write_openmm_files("x")
        cdamfit.target = mol
        cdamfit.water_model = "spce"
        cdamfit.coordination_type = "octahedral"
        cdamfit.dummy_atom_charge = 0.9
        cdamfit.calculator = _ToyCalculator()
        cdamfit.cache = False
        cdamfit.search_method = "staged"
        cdamfit.samples_per_stage = 4
        cdamfit.max_iterations = 1
        # the trailing comma makes a one-element tuple; it is unwrapped
        cdamfit.properties = {
            "IOD": {"value": 1.90, "tolerance": 0.05, "unit": "angstrom",
                    "weight": 1.0},
            "CN": {"value": 6.0, "tolerance": 0.3, "weight": 1.0}},
        cdamfit.fit()
        result = cdamfit.fitting_result()

        assert cdamfit.candidate.metal == "Al"
        assert cdamfit.candidate.formal_charge == 3
        assert cdamfit.candidate.water_model == "spce"
        assert cdamfit.candidate.geometry == "octahedral"
        assert cdamfit.candidate.dummy_charge == 0.9
        assert "dummy_charge" not in _fitted(cdamfit)
        assert cdamfit._fit_result is result
        assert result.diagnostics["n_calculations"] == 4

        cdamfit.write_openmm_files(str(tmp_path / "al"), mol_name="ALX")
        pdb, xml = tmp_path / "al.pdb", tmp_path / "al.xml"
        assert " ALX A" in pdb.read_text()
        assert '<Residue name="ALX">' in xml.read_text()
        sigma = result.parameters.to_dict()["metal_lj_sigma"]
        assert f'sigma="{sigma / 10:.12g}"' in xml.read_text()

        cdamfit.fitting_report()
        assert not list(tmp_path.glob("*.png"))
        cdamfit.plot_optimization(str(tmp_path / "al"))
        assert (tmp_path / "al_optimization.png").is_file()

    def test_molecule_target_must_be_a_charged_ion(self):
        import pytest
        from veloxchem import Molecule, OutputStream
        from veloxchem.errorhandler import VeloxChemError

        cdamfit = AutoCDAM(ostream=OutputStream(None))
        cdamfit.properties = {"IOD": {"value": 1.9, "tolerance": 0.05}}
        cdamfit.target = Molecule.read_smiles("[Al]")
        with pytest.raises(VeloxChemError, match="set the ion charge"):
            cdamfit.fit()
        cdamfit.target = Molecule.read_smiles("O")
        with pytest.raises(VeloxChemError, match="single metal ion"):
            cdamfit.fit()
        cdamfit.target = None
        with pytest.raises(VeloxChemError, match="both target and properties"):
            cdamfit.fit()
        cdamfit.target = Molecule.read_smiles("[Al+3]")
        cdamfit.coordination_type = "hexagonal"
        with pytest.raises(VeloxChemError, match="coordination_type must be"):
            cdamfit.fit()

    def test_dummy_charge_is_fitted_only_when_a_target_determines_it(self):
        import pytest
        from veloxchem import OutputStream
        from veloxchem.errorhandler import VeloxChemError

        iod = {"IOD": {"value": 1.9, "tolerance": 0.05}}
        hfe = {"HFE": {"value": -4525.0, "tolerance": 50.0, "unit": "kJ/mol"}}
        cdamfit = AutoCDAM(ostream=OutputStream(None))

        cdamfit.setup(CDAMCandidate("Al", 3), iod)
        assert set(_fitted(cdamfit)) == {"metal_lj_sigma", "metal_lj_epsilon"}
        assert cdamfit.candidate.dummy_charge == 0.5

        cdamfit.setup(CDAMCandidate("Al", 3), {**iod, **hfe})
        assert _fitted(cdamfit)["dummy_charge"] == pytest.approx((0.4, 1.0))

        cdamfit.setup(CDAMCandidate("Zn", 2, geometry="tetrahedral"), hfe)
        assert _fitted(cdamfit)["dummy_charge"] == pytest.approx((0.4, 1.0))

        cdamfit.dummy_atom_charge = 0.7
        cdamfit.setup(CDAMCandidate("Al", 3), {**iod, **hfe})
        assert "dummy_charge" not in _fitted(cdamfit)
        assert cdamfit.candidate.dummy_charge == 0.7
        with pytest.raises(VeloxChemError, match="also fit it"):
            cdamfit.setup(CDAMCandidate("Al", 3), iod,
                          parameters={"dummy_charge": (0.4, 1.0)})

    def test_lj12_6_4_fits_c4_by_default(self):
        import pytest
        from veloxchem import OutputStream

        iod = {"IOD": {"value": 1.9, "tolerance": 0.05}}
        cdamfit = AutoCDAM(ostream=OutputStream(None))
        cdamfit.setup(CDAMCandidate("Al", 3), iod)
        assert "c4_metal_water" not in _fitted(cdamfit)

        cdamfit.setup(CDAMCandidate("Al", 3, potential_form="lj12_6_4"), iod)
        spec = next(s for s in cdamfit.parameter_specs
                    if s.name == "c4_metal_water")
        assert spec.optimizable and spec.expandable
        assert tuple(spec.bounds) == pytest.approx((0.0, 200.0))
        assert tuple(spec.hard_bounds) == pytest.approx((0.0, 3000.0))

        cdamfit.setup(CDAMCandidate("Al", 3, potential_form="lj12_6_4"), iod,
                      parameters={"sigma": (1.5, 2.5), "c4": 50.0})
        spec = next(s for s in cdamfit.parameter_specs
                    if s.name == "c4_metal_water")
        assert not spec.optimizable and spec.initial == 50.0

    def test_metal_dummy_distance_is_fixed_or_fitted(self):
        import pytest
        from veloxchem import OutputStream
        from veloxchem.errorhandler import VeloxChemError

        iod = {"IOD": {"value": 1.9, "tolerance": 0.05}}
        cdamfit = AutoCDAM(ostream=OutputStream(None))
        cdamfit.setup(CDAMCandidate("Al", 3), iod)
        assert cdamfit.candidate.metal_dummy_distance == 0.9
        assert "metal_dummy_distance" not in _fitted(cdamfit)

        cdamfit.metal_dummy_distance = 0.95
        cdamfit.setup(CDAMCandidate("Al", 3), iod)
        assert cdamfit.candidate.metal_dummy_distance == 0.95
        assert "metal_dummy_distance" not in _fitted(cdamfit)

        cdamfit.metal_dummy_distance = (0.7, 1.1)
        cdamfit.setup(CDAMCandidate("Al", 3), iod)
        assert _fitted(cdamfit)["metal_dummy_distance"] == (0.7, 1.1)
        with pytest.raises(VeloxChemError, match="disagree"):
            cdamfit.setup(CDAMCandidate("Al", 3), iod,
                          parameters={"sigma": (1.5, 2.5)})
        cdamfit.metal_dummy_distance = (1.1, 0.7)
        with pytest.raises(VeloxChemError, match="lower, upper"):
            cdamfit.setup(CDAMCandidate("Al", 3), iod)

    def test_study_id_is_reproducible_and_kept_in_the_run_file(self, tmp_path):
        from veloxchem import Molecule, OutputStream
        from veloxchem.autocdam.models import read_h5_group

        def run(filename, dummy_charge=0.9):
            cdamfit = AutoCDAM(ostream=OutputStream(None))
            cdamfit.target = Molecule.read_smiles("[Al+3]")
            cdamfit.calculator = _ToyCalculator()
            cdamfit.dummy_atom_charge = dummy_charge
            cdamfit.samples_per_stage = 2
            cdamfit.max_iterations = 1
            cdamfit.filename = str(filename)
            cdamfit.properties = {"IOD": {"value": 1.9, "tolerance": 0.05}}
            cdamfit.fit()
            return cdamfit

        first = run(tmp_path / "Al0.9d0.9e")
        assert (tmp_path / "Al0.9d0.9e.h5").is_file()
        assert first.study_id.startswith("Al3_octahedral_tip3p_lj12_6_")
        assert first.fitting_result().diagnostics["study_id"] == first.study_id
        saved = read_h5_group(tmp_path / "Al0.9d0.9e.h5", "autocdam/checkpoint")
        assert first.study_id.endswith(saved["study_hash"])
        assert run(tmp_path / "again").study_id == first.study_id
        assert run(tmp_path / "other", 0.8).study_id != first.study_id

        text = repr(first.fitting_result())
        assert "metal_lj_sigma" in text and "angstrom" in text
        assert "kcal/mol" in text and first.study_id in text

    def test_pilot_and_fast_mode_md_protocols(self):
        import pytest
        from veloxchem import OutputStream
        from veloxchem.autocdam.md import md_settings
        from veloxchem.errorhandler import VeloxChemError

        pilot = md_settings({"pilot": True})
        assert pilot["production_steps"] * pilot["timestep_fs"] / 1000 == 0.05
        with pytest.raises(ValueError, match="not both"):
            md_settings({"pilot": True, "fast_mode": True})

        md = md_settings({"fast_mode": True})
        assert md["timestep_fs"] == 2.0
        assert md["production_steps"] * md["timestep_fs"] / 1000 == 50.0
        assert md["nvt_steps"] == md_settings(None)["nvt_steps"]

        cdamfit = AutoCDAM(ostream=OutputStream(None))
        cdamfit.fast_mode = True
        cdamfit.setup(CDAMCandidate("Al", 3),
                      {"IOD": {"value": 1.9, "tolerance": 0.05}})
        assert cdamfit.calculator_options["md_options"] == md

        cdamfit.pilot = True
        with pytest.raises(VeloxChemError, match="not both"):
            cdamfit.setup(CDAMCandidate("Al", 3),
                          {"IOD": {"value": 1.9, "tolerance": 0.05}})


def _fitted(cdamfit):
    """Bounds of the fitted parameters by name."""
    return {spec.name: tuple(spec.bounds) for spec in cdamfit.parameter_specs
            if spec.optimizable}
