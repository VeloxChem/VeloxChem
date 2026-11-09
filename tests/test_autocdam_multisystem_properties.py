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
"""Step 5: structural properties in framework systems and multisystem fitting."""
import json
from pathlib import Path

import numpy as np
import pytest
import veloxchem as vlx

from veloxchem.autocdam import AutoCDAM, CalculationResult, CDAMCandidate, FitTarget, FrameworkSystem
from veloxchem.autocdam.calculator import Calculator
from veloxchem.autocdam.multisystem import MultiSystemCalculator
from veloxchem.autocdam.framework import compare_to_baseline

STRUCTURES = Path(__file__).parent / "data"
CHARGES = json.loads(
    (STRUCTURES /
     "autocdam_al_node_fragment_charges.json").read_text())["charges_e"]


def molecule():
    mol = vlx.Molecule.read_xyz_file(
        str(STRUCTURES / "autocdam_al_node_fragment.xyz"))
    mol.set_charge(-1)
    return mol


def node(**kwargs):
    kwargs.setdefault("representation", "cdam")
    return FrameworkSystem(structure=molecule(),
                           charges=CHARGES,
                           target_metal="Al",
                           **kwargs)


def target(name, prop, value=0.0, tolerance=0.1, **extra):
    return FitTarget.from_mapping(name, {
        "property": prop,
        "value": value,
        "tolerance": tolerance,
        **extra
    })


def point(sigma=2.0, epsilon=0.2, q=0.5):
    return CDAMCandidate("Al",
                         3,
                         metal_lj_sigma=sigma,
                         metal_lj_epsilon=epsilon,
                         dummy_charge=q)


STRUCTURAL = [
    target("al_o_carboxylate", "distance", atoms=["Al", "O3"]),
    target("al_o_hydroxide", "distance", atoms=["Al", "O1"]),
    target("o_al_o", "angle", atoms=["O1", "Al", "O2"]),
    target("cn", "coordination_number", cutoff_angstrom=2.5),
    target("shell_rmsd", "coordination_rmsd", selection="first_shell"),
]


class FakeAqueous(Calculator):
    """Aqueous IOD/CN surface depending on the shared CDAM parameters."""

    def __init__(self):
        super().__init__()
        self.seen = []

    def calculate(self, batch, parameters, targets, **kwargs):
        self.seen.append(
            (parameters["metal_lj_sigma"], sorted(t.name for t in targets)))
        sigma = batch.system.nonbonded_parameters["metal_lj_sigma_angstrom"]
        q = batch.system.charges[1]
        values = {"IOD": 0.9 + 0.5 * sigma, "CN": 4.0 + 4.0 * q}
        return CalculationResult(
            values={t.name: values[t.name] for t in targets},
            metadata={t.name: {
                "quality_status": "usable"
            } for t in targets})


class SpyFramework(FrameworkSystem):

    def __post_init__(self):
        super().__post_init__()
        self.sigmas = []

    def calculate(self, candidate, targets, protocol=None):
        self.sigmas.append(candidate.metal_lj_sigma)
        return super().calculate(candidate, targets, protocol=protocol)


class TestAutoCdamMultisystemProperties:

    @pytest.mark.parametrize("representation", ["atom", "cdam"])
    def test_minimized_structural_properties(self, representation):
        system = node(representation=representation)
        result = system.calculate(point(), STRUCTURAL)
        assert not result.failed, result.failure_reason
        values = result.values
        assert 1.5 < values["al_o_carboxylate"] < 3.0 and 1.5 < values[
            "al_o_hydroxide"] < 3.0
        assert 90.0 < values["o_al_o"] <= 180.0
        assert values["cn"] == int(values["cn"]) and 0 <= values["cn"] <= 8
        assert values["shell_rmsd"] >= 0.0
        meta = result.metadata["al_o_carboxylate"]
        assert meta["geometry"] == "minimized" and meta["unit"] == "angstrom"
        assert meta["non_target_preserved"] and meta[
            "representation"] == representation

    def test_larger_metal_sigma_pushes_ligands_out(self):
        system = node()
        probe = [STRUCTURAL[0]]
        short = system.calculate(point(sigma=1.5),
                                 probe).values["al_o_carboxylate"]
        long = system.calculate(point(sigma=2.6),
                                probe).values["al_o_carboxylate"]
        assert long > short

    def test_md_average_is_distinguished_from_minimized_geometry(self):
        system = node(
            protocol={
                "equilibration_steps": 100,
                "production_steps": 400,
                "sample_interval": 100,
                "timestep_fs": 0.5
            })
        result = system.calculate(point(), [
            target("d_md", "distance", atoms=["Al", "O3"], method="md"),
            target("d_min", "distance", atoms=["Al", "O3"]),
        ])
        assert not result.failed, result.failure_reason
        assert result.metadata["d_md"]["geometry"] == "md_average"
        assert result.metadata["d_md"]["n_frames"] == 4
        assert result.metadata["d_min"]["geometry"] == "minimized"

    def test_periodic_distances_use_the_minimum_image(self):
        system = node(box_vectors=np.eye(3) * 25.0)
        baseline = system.prepare_baseline()
        xyz = np.array(baseline.positions_nm)
        measured = system.measure(STRUCTURAL[0], xyz, site=None)
        xyz[5] += baseline.box_vectors_nm[0]  # carboxylate O moved by one cell
        assert system.measure(STRUCTURAL[0], xyz,
                              site=None) == pytest.approx(measured)
        assert measured == pytest.approx(1.90, abs=1e-6)

    def test_coordination_rmsd_against_a_shuffled_rotated_reference(self):
        system = node()
        baseline = system.prepare_baseline()
        atoms = [0, *system.first_shell(6)]
        local = (baseline.positions_nm[atoms] - baseline.positions_nm[0]) * 10.0
        angle = np.radians(37.0)
        rotation = np.array([[np.cos(angle), -np.sin(angle), 0],
                             [np.sin(angle), np.cos(angle), 0], [0, 0, 1]])
        order = [0, 4, 2, 6, 1, 5, 3]  # ligands in a different order
        reference = [["Al" if i == 0 else "O", *(rotation @ local[i] + 3.0)]
                     for i in order]
        shell = target("shell",
                       "coordination_rmsd",
                       reference_structure=reference)
        assert system.measure(shell, baseline.positions_nm,
                              site=None) == pytest.approx(0.0, abs=1e-8)
        wrong = target("bad",
                       "coordination_rmsd",
                       reference_structure=[["Al", 0, 0, 0]] * 7)
        with pytest.raises(ValueError, match="elements"):
            system.measure(wrong, baseline.positions_nm, site=None)

    @pytest.mark.parametrize("bad, message", [
        (target("d", "distance", atoms=["Al"]), "needs 2 atoms"),
        (target("a", "angle", atoms=["O1", "Al", "Q9"]), "unknown atom label"),
        (target("x", "dihedral", atoms=[0, 1, 2, 3]), "supports"),
        (target("m", "distance", atoms=[0, 1], method="anneal"), "method"),
    ])
    def test_invalid_structural_targets_are_reported(self, bad, message):
        with pytest.raises(ValueError, match=message):
            node().check_targets([bad])

    def test_router_merges_systems_and_names_failures(self):
        mof = node()
        systems = {"water": None, "mof": mof}
        router = MultiSystemCalculator(systems,
                                       FakeAqueous(),
                                       default_system="water")
        targets = [
            FitTarget.from_mapping(
                "water_iod", {
                    "system": "water",
                    "property": "IOD",
                    "value": 1.9,
                    "tolerance": 0.03
                }),
            FitTarget.from_mapping("iod_again", {
                "property": "IOD",
                "value": 1.9,
                "tolerance": 0.1
            }),
            FitTarget.from_mapping(
                "al_o", {
                    "system": "mof",
                    "property": "distance",
                    "atoms": ["Al", "O3"],
                    "value": 1.9,
                    "tolerance": 0.05
                })
        ]
        router.check_targets(targets)
        batch = point(sigma=2.0).to_generated_batch()
        from veloxchem.autocdam import ParameterSet
        result = router.calculate(batch, ParameterSet({"metal_lj_sigma": 2.0}),
                                  targets)
        assert not result.failed
        assert result.values["water_iod"] == result.values[
            "iod_again"] == pytest.approx(1.9)
        # one aqueous IOD calculation
        assert router.aqueous.seen == [(2.0, ["IOD"])]
        assert result.metadata["al_o"]["system"] == "mof"
        broken = FitTarget.from_mapping(
            "bad", {
                "system": "mof",
                "property": "distance",
                "atoms": ["Al", "O3"],
                "value": 1.9,
                "tolerance": 0.05
            })
        wrong_metal = CDAMCandidate("Zn",
                                    2,
                                    metal_lj_sigma=2.0,
                                    metal_lj_epsilon=0.2)
        failed = router.calculate(wrong_metal.to_generated_batch(),
                                  ParameterSet(), [broken])
        assert failed.failed and failed.failure_reason.startswith("[mof]")
        assert "does not match target" in failed.failure_reason

    def test_aqueous_targets_sharing_a_property_keep_their_own_units(self):
        router = MultiSystemCalculator({"water": None},
                                       FakeAqueous(),
                                       default_system="water")
        targets = [
            FitTarget.from_mapping("iod_nm", {
                "property": "IOD",
                "value": 0.19,
                "tolerance": 0.003,
                "unit": "nm"
            }),
            FitTarget.from_mapping(
                "iod_a", {
                    "property": "IOD",
                    "value": 1.9,
                    "tolerance": 0.03,
                    "unit": "angstrom"
                })
        ]
        from veloxchem.autocdam import ParameterSet
        result = router.calculate(
            point(sigma=2.0).to_generated_batch(),
            ParameterSet({"metal_lj_sigma": 2.0}), targets)
        assert not result.failed, result.failure_reason
        assert result.values["iod_a"] == pytest.approx(1.9)
        assert result.values["iod_nm"] == pytest.approx(0.19)
        assert result.metadata["iod_nm"]["unit"] == "nm"
        assert result.metadata["iod_a"]["unit"] == "angstrom"
        assert router.aqueous.seen == [(2.0, ["IOD"])]  # still one calculation

    def test_unknown_system_or_aqueous_property_fails_before_fitting(self):
        unknown_system = {
            "d": {
                "system": "mof",
                "property": "distance",
                "atoms": [0, 1],
                "value": 1.9,
                "tolerance": 0.1
            }
        }
        with pytest.raises(ValueError, match="unknown system"):
            AutoCDAM().setup(CDAMCandidate("Al", 3),
                             unknown_system,
                             parameters={"sigma": (1.4, 2.6)},
                             systems={"water": {
                                 "type": "aqueous"
                             }})
        aqueous_distance = {
            "d": {
                "property": "distance",
                "value": 1.9,
                "tolerance": 0.1
            }
        }
        with pytest.raises(ValueError, match="aqueous systems support"):
            AutoCDAM().setup(CDAMCandidate("Al", 3),
                             aqueous_distance,
                             parameters={"sigma": (1.4, 2.6)},
                             systems={"water": {
                                 "type": "aqueous"
                             }})

    def test_atom_representation_flags_unidentifiable_dummy_charge(self):
        parameters = {"sigma": (1.4, 2.6), "dummy_charge": (0.2, 0.7)}
        targets = {
            "d": {
                "property": "distance",
                "atoms": ["Al", "O3"],
                "value": 1.93,
                "tolerance": 0.05
            }
        }
        fit = AutoCDAM()
        fit.setup(CDAMCandidate("Al", 3),
                  targets,
                  parameters=parameters,
                  systems={"mof": node(representation="atom")})
        assert any(
            "dummy_charge" in warning for warning in fit.parameter_warnings)

    def test_one_parameter_set_is_fitted_against_water_and_mof_together(self):
        mof = SpyFramework(structure=molecule(),
                           charges=CHARGES,
                           target_metal="Al",
                           representation="cdam",
                           name="mof")
        aqueous = FakeAqueous()
        systems = {"water": {"type": "aqueous"}, "mof": mof}
        parameters = {
            "sigma": (1.5, 2.5),
            "dummy_charge": (0.3, 0.7),
            "metal_dummy_distance": 0.9
        }
        targets = {
            "IOD": {
                "system": "water",
                "property": "IOD",
                "value": 1.9,
                "tolerance": 0.03,
                "role": "required"
            },
            "CN": {
                "system": "water",
                "property": "CN",
                "value": 6.0,
                "tolerance": 0.2
            },
            "al_o": {
                "system": "mof",
                "property": "distance",
                "atoms": ["Al", "O3"],
                "value": 1.93,
                "tolerance": 0.05,
                "weight": 2.0
            },
            "shell": {
                "system": "mof",
                "property": "coordination_rmsd",
                "value": 0.0,
                "tolerance": 0.1,
                "role": "monitor"
            },
        }
        fit = AutoCDAM()
        fit.calculator = aqueous
        fit.update_settings({
            "search_method": "staged",
            "max_iterations": 1,
            "samples_per_stage": 4
        })
        fit.setup(CDAMCandidate("Al", 3),
                  targets,
                  parameters=parameters,
                  systems=systems)
        baseline_xml = mof.prepare_baseline().system_xml
        result = fit.fit()
        samples = fit.fitting_driver.table.samples
        assert len(samples) == 4 and all(
            not s.calculation.failed for s in samples)
        sigmas = [s.parameters["metal_lj_sigma"] for s in samples]
        assert sigmas == mof.sigmas[:4]
        assert [seen[0] for seen in aqueous.seen] == mof.sigmas[:4]
        assert all(
            set(s.calculation.values) == {"IOD", "CN", "al_o"} for s in samples)
        assert "shell" in result.diagnostics["monitor"]
        report = result.diagnostics["systems"]["mof"]
        assert report["non_target_preserved"] and report[
            "n_checked_calculations"] == 4
        assert mof.baseline_builds == 1 and mof.prepare_baseline(
        ).system_xml == baseline_xml
        system, _, _ = mof.build(
            CDAMCandidate(
                "Al", 3, **{
                    "metal_lj_sigma": result.parameters["metal_lj_sigma"],
                    "dummy_charge": result.parameters["dummy_charge"]
                }))
        assert compare_to_baseline(system, mof.prepare_baseline(),
                                   changed={0}) == []

    @pytest.mark.timeconsuming
    def test_real_aqueous_md_and_framework_in_one_calculation(self, tmp_path):
        systems = {"water": {"type": "aqueous"}, "mof": node()}
        targets = {
            "IOD": {
                "system": "water",
                "value": 1.9,
                "tolerance": 0.05
            },
            "al_o": {
                "system": "mof",
                "property": "distance",
                "atoms": ["Al", "O3"],
                "value": 1.93,
                "tolerance": 0.05
            }
        }
        fit = AutoCDAM()
        fit.update_settings({"pilot": True, "work_dir": str(tmp_path)})
        fit.setup(CDAMCandidate("Al", 3),
                  targets,
                  parameters={"sigma": (1.5, 2.5)},
                  systems=systems)
        calculation = fit.calculate({"metal_lj_sigma": 2.0})
        assert not calculation.failed, calculation.failure_reason
        assert {"IOD", "al_o"} <= set(calculation.values)
