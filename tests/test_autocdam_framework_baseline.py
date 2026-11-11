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
"""Step 3: VeloxChem baseline generation and target-only parameter updates."""
import json
from pathlib import Path

import numpy as np
import pytest
import veloxchem as vlx

from veloxchem.autocdam import CDAMCandidate, FrameworkSystem, ParameterSpec
from veloxchem.autocdam.framework import compare_to_baseline, covalent_groups
from veloxchem.autocdam.geometry import atom_labels, minimum_image
from veloxchem.errorhandler import VeloxChemError

STRUCTURES = Path(__file__).parent / "data"
CHARGES = json.loads(
    (STRUCTURES /
     "autocdam_al_node_fragment_charges.json").read_text())["charges_e"]
# hydroxides (O1, O2) and formates (O3, O5) are -1; the waters are neutral
GROUP_CHARGES = {"O1": -1, "O2": -1, "O3": -1, "O5": -1}


def node_molecule(copies=1, shift=9.0):
    lines = (STRUCTURES /
             "autocdam_al_node_fragment.xyz").read_text().splitlines()[2:]
    rows = []
    for copy in range(copies):
        for line in lines:
            element, *xyz = line.split()
            rows.append(
                f"{element} {float(xyz[0]) + copy * shift} {xyz[1]} {xyz[2]}")
    molecule = vlx.Molecule.read_xyz_string(f"{len(rows)}\n\n" +
                                            "\n".join(rows))
    molecule.set_charge(-copies)
    return molecule


def node_system(**kwargs):
    kwargs.setdefault("charges", CHARGES)
    kwargs.setdefault("target_metal", "Al")
    return FrameworkSystem(structure=node_molecule(), **kwargs)


def nonbonded(system):
    import openmm as mm
    return next(
        f for f in system.getForces() if isinstance(f, mm.NonbondedForce))


def particle(system, index):
    q, sigma, epsilon = nonbonded(system).getParticleParameters(index)
    return q._value, sigma._value, epsilon._value


def point(sigma=2.2, epsilon=0.1):
    return CDAMCandidate("Al",
                         3,
                         metal_lj_sigma=sigma,
                         metal_lj_epsilon=epsilon)


class TestAutoCdamFrameworkBaseline:

    def test_labels_are_veloxchem_atom_names(self):
        molecule = vlx.Molecule.read_str("Al 0 0 0\nO 2 0 0\nH 3 0 0\nO 0 2 0")
        assert atom_labels(molecule) == ["Al", "O1", "H", "O2"]

    def test_veloxchem_baseline_keeps_charges_and_removes_metal_bonded_terms(
            self):
        system = node_system()
        baseline = system.prepare_baseline()
        assert "VeloxChem" in baseline.source and "user-supplied" in baseline.source
        assert baseline.charges == pytest.approx(CHARGES, abs=1e-6)
        assert baseline.total_charge == pytest.approx(-1.0)
        assert baseline.removed_metal_terms["bonds"] == 6  # six Al-O contacts
        assert not any(0 in bond for bond in baseline.bonds)
        assert len(baseline.bonds) == 12  # 2 O-H, 2x3 formate, 2x2 water
        omm = baseline.new_system()
        exceptions = [
            nonbonded(omm).getExceptionParameters(i)[:2]
            for i in range(nonbonded(omm).getNumExceptions())
        ]
        assert exceptions and not any(0 in pair for pair in exceptions)
        # VeloxChem's own LJ for every non-target atom is retained exactly.
        assert particle(omm, 1)[1] > 0 and particle(omm, 1)[2] > 0

    def test_candidate_updates_only_the_target_and_preserve_everything_else(
            self):
        system = node_system()
        baseline = system.prepare_baseline()
        for sigma, epsilon in ((1.6, 0.05), (2.4, 0.4)):
            omm, positions, site = system.build(point(sigma, epsilon))
            q, s, e = particle(omm, 0)
            assert q == pytest.approx(
                baseline.charges[0])  # charge fixed in atom mode
            assert s == pytest.approx(sigma * 0.1) and e == pytest.approx(
                epsilon * 4.184)
            assert compare_to_baseline(omm, baseline, changed={0}) == []
            assert np.allclose(positions, baseline.positions_nm)
            total = sum(
                particle(omm, i)[0] for i in range(omm.getNumParticles()))
            assert total == pytest.approx(-1.0)

    def test_comparison_detects_a_changed_neighbour(self):
        system = node_system()
        omm, _, _ = system.build(point())
        q, s, e = particle(omm, 1)
        nonbonded(omm).setParticleParameters(1, q + 0.1, s, e)
        assert "nonbonded parameters of particle 1" in compare_to_baseline(
            omm, system.prepare_baseline(), changed={0})

    def test_target_label_selects_the_right_metal_among_several(self):
        charges = CHARGES * 2
        system = FrameworkSystem(structure=node_molecule(copies=2),
                                 charges=charges,
                                 target_metal="Al2")
        baseline = system.prepare_baseline()
        assert system.target_index == 19 and baseline.metal_indices == (0, 19)
        omm, _, _ = system.build(point(1.7, 0.3))
        assert particle(omm, 19)[1] == pytest.approx(0.17)
        assert particle(omm,
                        0) == pytest.approx(particle(baseline.new_system(), 0))
        assert compare_to_baseline(omm, baseline, changed={19}) == []

    def test_baseline_is_generated_once_and_never_mutated(self, monkeypatch):
        calls = []
        original = vlx.MMForceFieldGenerator.create_topology

        def counted(self, *args, **kwargs):
            calls.append(1)
            return original(self, *args, **kwargs)

        monkeypatch.setattr(vlx.MMForceFieldGenerator, "create_topology",
                            counted)
        system = node_system()
        baseline = system.prepare_baseline()
        snapshot = baseline.system_xml
        for sigma in (1.5, 1.9, 2.5):
            omm, positions, _ = system.build(point(sigma))
            positions += 1.0  # candidate copies are writable
        assert len(calls) == 1 and system.baseline_builds == 1
        assert system.prepare_baseline(
        ) is baseline and baseline.system_xml == snapshot
        with pytest.raises(ValueError):
            baseline.positions_nm[0, 0] = 5.0  # snapshot is read-only

    def test_supplied_force_field_is_copied_not_modified(self):
        ff = vlx.MMForceFieldGenerator()
        ff.ostream.mute()
        ff.partial_charges = np.array(CHARGES)
        ff.create_topology(node_molecule(), resp=True)
        before = {i: dict(atom) for i, atom in ff.atoms.items()}
        system = FrameworkSystem(structure=node_molecule(),
                                 target_metal=0,
                                 baseline_force_field=ff)
        assert system.prepare_baseline(
        ).source == "supplied VeloxChem force field"
        system.build(point())
        assert {i: dict(atom) for i, atom in ff.atoms.items()} == before

    def test_atom_representation_reports_ineffective_cdam_parameters(self):
        specs = [
            ParameterSpec("metal_lj_sigma", 2.0, (1.0, 3.0)),
            ParameterSpec("dummy_charge", 0.5, (0.2, 0.7)),
            ParameterSpec("metal_dummy_distance", 0.9, optimizable=False)
        ]
        assert node_system().ineffective_parameters(specs) == ["dummy_charge"]

    @pytest.mark.parametrize("kwargs, error, message", [
        ({
            "charges": None
        }, VeloxChemError, "choose baseline charges"),
        ({
            "target_metal": "O1"
        }, VeloxChemError, "not a metal"),
        ({
            "target_metal": "AL9"
        }, ValueError, "unknown atom label"),
        ({
            "charges": [0.0] * 19
        }, ValueError, "sum to"),
        ({
            "charges": [0.1] * 3
        }, ValueError, "19 finite values"),
        ({
            "representation": "united"
        }, VeloxChemError, "representation"),
        ({
            "charges": "groups"
        }, VeloxChemError, "needs the formal metal charges"),
        ({
            "charges": "groups",
            "metal_charges": {"O1": 3}
        }, ValueError, "not a metal"),
        ({
            "charges": "groups",
            "metal_charges": {"Al": 3}
        }, ValueError, "odd number of electrons"),
        ({
            "charges": "groups",
            "metal_charges": {"Al": 2},
            "group_charges": GROUP_CHARGES
        }, ValueError, "sum to"),
        ({
            "charges": "groups",
            "metal_charges": {"Al": 3},
            "group_charges": {"O1": -1, "H1": 0}
        }, ValueError, "two charges"),
        ({
            "water_model": "tip4pew"
        }, VeloxChemError, "water_model"),
    ])
    def test_invalid_framework_inputs_fail_clearly(self, kwargs, error,
                                                   message):
        with pytest.raises(error, match=message):
            node_system(**kwargs).prepare_baseline()

    def test_group_charge_scheme_gives_formal_metal_and_integer_groups(self):
        node = FrameworkSystem(structure=node_molecule(),  # default charges
                               target_metal="Al",
                               metal_charges={"Al": 3},
                               group_charges=GROUP_CHARGES)
        assert node.charges == "groups"
        assert covalent_groups(node.molecule) == [[1, 2], [3, 4],
                                                  [5, 6, 7, 8],
                                                  [9, 10, 11, 12],
                                                  [13, 14, 15],
                                                  [16, 17, 18]]
        charges = np.array(node.prepare_baseline().charges)
        assert charges[0] == pytest.approx(3.0)
        for group, total in (([1, 2], -1), ([3, 4], -1), ([5, 6, 7, 8], -1),
                             ([9, 10, 11, 12], -1), ([13, 14, 15], 0),
                             ([16, 17, 18], 0)):
            assert charges[group].sum() == pytest.approx(total, abs=1e-6)
        assert charges[[13, 14, 15]] == pytest.approx([-0.8476, 0.4238, 0.4238])
        # the test data are this scheme's output
        assert charges == pytest.approx(CHARGES, abs=1e-4)
        same = node_system(charges="groups",
                           metal_charges={0: 3},
                           group_charges=GROUP_CHARGES,
                           water_model="tip3p").prepare_baseline().charges
        assert same[13:16] == pytest.approx((-0.834, 0.417, 0.417))
        assert same[:13] == pytest.approx(tuple(charges[:13]))

    def test_precomputed_group_charges_are_checked_not_recomputed(
            self, monkeypatch):

        def no_resp(*args, **kwargs):
            raise AssertionError("RESP must not run for precomputed charges")

        monkeypatch.setattr(FrameworkSystem, "_isolated_group_resp",
                            staticmethod(no_resp))
        node = node_system(charges=CHARGES,
                           metal_charges={"Al": 3},
                           group_charges=GROUP_CHARGES)
        assert node.prepare_baseline().charges == pytest.approx(tuple(CHARGES))
        # the old illustrative charges (Al +1.7, OH- -0.6) are not this scheme
        illustrative = [1.7, -1.0, 0.4, -1.0, 0.4, -0.7, 0.6, -0.7, 0.05,
                        -0.7, 0.6, -0.7, 0.05, -0.8, 0.4, 0.4, -0.8, 0.4, 0.4]
        with pytest.raises(ValueError, match="does not match the group scheme"):
            node_system(charges=illustrative,
                        metal_charges={"Al": 3},
                        group_charges=GROUP_CHARGES).prepare_baseline()

    def test_periodic_cell_uses_pme_and_minimum_image(self):
        import openmm as mm

        box = np.eye(3) * 25.0
        system = node_system(box_vectors=box)
        baseline = system.prepare_baseline()
        assert baseline.periodic
        force = nonbonded(baseline.new_system())
        assert force.getNonbondedMethod() == mm.NonbondedForce.PME
        assert force.getCutoffDistance()._value == pytest.approx(1.0)
        delta = minimum_image(np.array([[2.3, -2.4, 0.1]]), box * 0.1)
        assert delta == pytest.approx(np.array([[-0.2, 0.1, 0.1]]))

    def test_bond_across_the_cell_boundary_is_refused(self):
        molecule = node_molecule()
        coordinates = molecule.get_coordinates_in_angstrom()
        coordinates[2] = [-5.4, 0.87,
                          0.0]  # hydroxide H wrapped to the far side
        moved = vlx.Molecule(molecule.get_labels(), coordinates, "angstrom")
        moved.set_charge(-1)
        with pytest.raises(ValueError, match="periodic boundary"):
            FrameworkSystem(structure=moved,
                            charges=CHARGES,
                            target_metal="Al",
                            box_vectors=np.eye(3) * 8.0).prepare_baseline()

    @pytest.mark.timeconsuming
    def test_resp_charges_are_generated_by_veloxchem(self):
        baseline = node_system(charges="resp").prepare_baseline()
        assert "RESP" in baseline.source
        assert baseline.total_charge == pytest.approx(-1.0, abs=1e-5)
        assert any(abs(q) > 0.1 for q in baseline.charges)
