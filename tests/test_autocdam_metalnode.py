from pathlib import Path
import json

import numpy as np
import pytest
import veloxchem as vlx

from veloxchem.autocdam import AutoCDAM, CDAMCandidate, MetalNode
from veloxchem.autocdam.multisystem import build_systems, expand_system_objects
from veloxchem.autocdam.models import coerce_targets
from veloxchem.errorhandler import VeloxChemError

DATA = Path(__file__).parent / "data"
XYZ = DATA / "autocdam_al_node_fragment.xyz"
CHARGES = json.loads(
    (DATA / "autocdam_al_node_fragment_charges.json").read_text())["charges_e"]


def al_node(**kwargs):
    return MetalNode.read_xyz_file(XYZ, charge=-1, **kwargs)


def two_al_node():
    """Two copies of the fragment, 9 A apart: two Al of the same element."""
    atoms = XYZ.read_text().splitlines()[2:]
    rows = atoms + [
        f"{e} {float(x) + 9.0} {y} {z}"
        for e, x, y, z in (line.split() for line in atoms)
    ]
    molecule = vlx.Molecule.read_xyz_string(f"{len(rows)}\n\n" +
                                            "\n".join(rows))
    return MetalNode(molecule, charge=-2)


def ammine_node():
    """The fragment with one water replaced by NH3 (not typed automatically)."""
    lines = XYZ.read_text().splitlines()
    atoms = lines[2:]
    atoms[13:16] = [
        "N     0.00000    0.00000    2.00000",
        "H     0.94000    0.00000    2.38000",
        "H    -0.47000    0.81400    2.38000",
        "H    -0.47000   -0.81400    2.38000",
    ]
    text = f"{len(atoms)}\n\n" + "\n".join(atoms)
    molecule = vlx.Molecule.read_xyz_string(text)
    return MetalNode(molecule, charge=-1)


class TestMetalNodeResidues:

    def test_residues_are_typed_with_charges_and_donors(self):
        node = al_node()
        kinds = [(r.kind, r.formula, r.charge) for r in node.residues]
        assert kinds == [("hydroxide", "HO", -1), ("hydroxide", "HO", -1),
                         ("carboxylate", "CHO2", -1),
                         ("carboxylate", "CHO2", -1), ("water", "H2O", 0),
                         ("water", "H2O", 0)]
        donors = [node.labels[d] for r in node.residues for d in r.donors]
        assert donors == ["O1", "O2", "O3", "O5", "O7", "O8"]
        assert all(r.bridging == 0 for r in node.residues)
        assert node._metal_charges() == {0: 3.0}
        node.check()

    def test_residue_table_prints_types_and_caps(self, capsys):
        al_node().print_residues()
        out = capsys.readouterr().out
        assert "carboxylate (formate)" in out and "hydroxide" in out
        assert "Unrecognized" not in out

    def test_unrecognized_residue_stops_until_the_user_sets_it(self):
        node = ammine_node()
        unknown = [r for r in node.residues if r.kind is None]
        assert len(unknown) == 1 and unknown[0].formula == "H3N"
        with pytest.raises(VeloxChemError, match="set_residue_type"):
            node.check()
        with pytest.raises(VeloxChemError, match="formal charge"):
            node.set_residue_type(unknown[0].number, "ammine")
        node.set_residue_type(unknown[0].number, "ammine", charge=0)
        node.check()
        assert unknown[0].label == "ammine"

    def test_charges_that_do_not_add_up_are_refused(self):
        node = al_node(metal_charges={"Al": 2})
        with pytest.raises(VeloxChemError, match="add up"):
            node.check()


class TestMetalNodeTargets:

    def test_grouped_distances_and_cis_trans_angle_rmsds(self):
        targets = al_node().targets("Al in MOF", {"weight": 2.0})
        assert sorted(targets) == [
            "Al in MOF: Al-O(carboxylate)", "Al in MOF: Al-O(hydroxide)",
            "Al in MOF: Al-O(water)", "Al in MOF: cis O-Al-O RMSD",
            "Al in MOF: trans O-Al-O RMSD"
        ]
        hydroxide = targets["Al in MOF: Al-O(hydroxide)"]
        assert hydroxide["atoms"] == [["Al", "O1"], ["Al", "O2"]]
        assert hydroxide["value"] == pytest.approx(1.9)
        assert hydroxide["tolerance"] == 0.05 and hydroxide["weight"] == 2.0
        assert len(targets["Al in MOF: cis O-Al-O RMSD"]["atoms"]) == 12
        assert len(targets["Al in MOF: trans O-Al-O RMSD"]["atoms"]) == 3
        assert targets["Al in MOF: cis O-Al-O RMSD"]["value"] == 0.0

    def test_targets_cover_every_metal_of_the_fitted_element(self):
        node = two_al_node()
        assert node.metal == 0 and node.fitted == (0, 19)
        targets = node.targets("node", {})
        hydroxide = targets["node: Al-O(hydroxide)"]
        assert hydroxide["atoms"] == [["Al1", "O1"], ["Al1", "O2"],
                                      ["Al2", "O9"], ["Al2", "O10"]]
        assert hydroxide["value"] == pytest.approx(1.9)
        cis = targets["node: cis O-Al-O RMSD"]["atoms"]
        assert len(cis) == 24 and {t[1] for t in cis} == {"Al1", "Al2"}
        assert len(targets["node: trans O-Al-O RMSD"]["atoms"]) == 6

    def test_reference_by_residue_type_and_tolerances(self):
        targets = al_node().targets(
            "node", {
                "reference": {"hydroxide": 1.87, "carboxylate": 1.93},
                "tolerance": {"distance": 0.02, "angle": 5.0},
                "role": "required"
            })
        assert targets["node: Al-O(hydroxide)"]["value"] == 1.87
        assert targets["node: Al-O(carboxylate)"]["value"] == 1.93
        assert targets["node: Al-O(water)"]["value"] == pytest.approx(1.9)
        assert targets["node: cis O-Al-O RMSD"]["tolerance"] == 5.0
        assert all(t["role"] == "required" for t in targets.values())

    @pytest.mark.parametrize("values, message", [
        ({"reference": {"oxo": 1.8}}, "no donor of type"),
        ({"atoms": ["Al", "O1"]}, "unsupported keys"),
    ])
    def test_bad_node_settings_are_reported(self, values, message):
        with pytest.raises(VeloxChemError, match=message):
            al_node().targets("node", values)


class TestMetalNodeInStudy:

    def test_expansion_gives_a_framework_system_and_ordinary_targets(self):
        node = al_node(charges=CHARGES)
        properties = {
            "IOD": {"value": 1.9, "tolerance": 0.05, "role": "required"},
            "Al in MOF": {"system": node, "weight": 1.0},
        }
        expanded, systems = expand_system_objects(properties, None, "tip3p")
        assert systems["water"] == {"type": "aqueous"}
        spec = systems["node"]
        assert spec["representation"] == "cdam" and spec["target_metal"] == "Al"
        assert spec["charges"] == CHARGES and "metal_charges" not in spec
        assert "IOD" in expanded and "Al in MOF" not in expanded
        targets = coerce_targets(expanded)
        assert all(t.metadata.get("system") in (None, "node") for t in targets)

    def test_group_scheme_and_charge_file(self, tmp_path):
        spec = al_node().system(water_model="TIP3P")
        assert spec["charges"] == "groups" and spec["water_model"] == "tip3p"
        assert spec["metal_charges"] == {"Al": 3.0}
        assert spec["group_charges"] == {"O1": -1, "O2": -1, "O3": -1, "O5": -1}
        charge_file = tmp_path / "charges.json"
        charge_file.write_text(json.dumps({
            "charges_e": CHARGES, "metal_charges": {"Al": 3},
            "group_charges": {"O1": -1}, "water_model": "spce"}))
        spec = al_node(charges=str(charge_file)).system(water_model="tip3p")
        assert spec["charges"] == CHARGES and spec["water_model"] == "spce"
        assert spec["group_charges"] == {"O1": -1}

    def test_waters_need_a_supported_water_model(self):
        with pytest.raises(VeloxChemError, match="water residues"):
            al_node().system(water_model="tip4pew")

    def test_node_targets_are_measured_on_the_minimized_framework(self):
        node = al_node(charges=CHARGES)
        properties = {"Al in MOF": {"system": node}}
        expanded, systems = expand_system_objects(properties, None, "spce")
        assert "water" not in systems
        targets = coerce_targets(expanded)
        router = build_systems(systems, targets, None, None)
        framework = router.systems["node"]
        result = framework.calculate(
            CDAMCandidate("Al", 3, metal_lj_sigma=2.0, metal_lj_epsilon=0.2,
                          dummy_charge=0.5, bond_force_constant=3000.0),
            targets)
        assert not result.failed, result.failure_reason
        values = result.values
        for name in ("hydroxide", "carboxylate", "water"):
            assert 1.5 < values[f"Al in MOF: Al-O({name})"] < 3.0
        assert 0.0 <= values["Al in MOF: cis O-Al-O RMSD"] < 45.0
        assert result.metadata["Al in MOF: cis O-Al-O RMSD"]["unit"] == "degree"

    def test_driver_setup_registers_the_node(self):
        fit = AutoCDAM()
        fit.update_settings({"cheap_first": False})
        fit.setup(CDAMCandidate("Al", 3, water_model="spce"), {
            "IOD": {"value": 1.9, "tolerance": 0.05, "role": "required"},
            "Al in MOF": {"system": al_node(charges=CHARGES)},
        })
        names = [t.name for t in fit.targets]
        assert "IOD" in names and "Al in MOF: Al-O(hydroxide)" in names
        assert set(fit.system_router.systems) == {"water", "node"}
