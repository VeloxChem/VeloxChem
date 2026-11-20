import json
from pathlib import Path

import numpy as np
import openmm as mm
import pytest

from veloxchem.autocdam import (AutoCDAM, CDAMCandidate, CoordinateGenerator,
                                MetalNode)
from veloxchem.autocdam.md import run_md
from veloxchem.autocdam.models import coerce_targets, minimum_bond_force_constant
from veloxchem.autocdam.multisystem import build_systems, expand_system_objects
from veloxchem.outputstream import OutputStream

DATA = Path(__file__).parent / "data"
XYZ = DATA / "autocdam_al_node_fragment.xyz"
CHARGES = json.loads(
    (DATA / "autocdam_al_node_fragment_charges.json").read_text())["charges_e"]

# Al3+ of Golo, Ahlquist, Su, ACS Omega 10, 3801 (2025), Table 1:
# Rmin/2 0.6158 A, epsilon 29.9565 kcal/mol, dummy 0.73 e at 0.74 A,
# Kb 1950 kcal/mol/A^2, D-M-D 200 kcal/mol/rad^2 (Amber, E = K x^2; AutoCDAM
# uses OpenMM's E = k/2 x^2, so k = 2 K), harmonic M-D bonds
PAPER_AL = {
    "metal_lj_sigma": 2 * 0.6158 / 2**(1 / 6),
    "metal_lj_epsilon": 29.9565,
    "dummy_charge": 0.73,
    "metal_dummy_distance": 0.74,
    "bond_force_constant": 3900.0,
    "angle_force_constant": 400.0,
}


def al(bond, **kwargs):
    return CDAMCandidate("Al", 3, metal_dummy_bond=bond,
                         **{**PAPER_AL, **kwargs})


def site_constraints(system, n_site=7):
    return [i for i in range(system.getNumConstraints())
            if system.getConstraintParameters(i)[0] < n_site and
            system.getConstraintParameters(i)[1] < n_site]


class TestHarmonicMetalDummyBond:

    def test_constrained_candidates_keep_their_batch_data(self):
        bonded = al("constraint").to_generated_batch().system.bonded_parameters
        assert "metal_dummy_bond" not in bonded
        batch = al("harmonic").to_generated_batch()
        assert batch.system.bonded_parameters["metal_dummy_bond"] == "harmonic"
        assert CDAMCandidate.from_generated_batch(
            batch).metal_dummy_bond == "harmonic"
        with pytest.raises(ValueError, match="metal_dummy_bond"):
            al("rigid")

    @pytest.mark.parametrize("bond, n_constraints", [("constraint", 6),
                                                     ("harmonic", 0)])
    def test_water_box_has_md_constraints_only_when_constrained(
            self, bond, n_constraints):
        systems = CoordinateGenerator(
            ostream=OutputStream(None)).create_openmm_template(al(bond))
        system = systems.solvated_system
        assert len(site_constraints(system)) == n_constraints
        bonds = [f for f in system.getForces()
                 if isinstance(f, mm.HarmonicBondForce)][0]
        site = [bonds.getBondParameters(i) for i in range(bonds.getNumBonds())
                if bonds.getBondParameters(i)[0] < 7]
        assert len(site) == 6
        assert site[0][3]._value == pytest.approx(3900.0 * 418.4)

    def test_paper_al_dummies_stay_on_their_bonds_in_water(self):
        systems = CoordinateGenerator(
            ostream=OutputStream(None)).create_openmm_template(al("harmonic"))
        trajectory = run_md(systems.solvated_system, systems.solvated_topology,
                            systems.solvated_positions, {
                                "timestep_fs": 1.0,
                                "nvt_steps": 500,
                                "npt_steps": 0,
                                "production_steps": 1000,
                                "report_interval": 100,
                                "random_seed": 7,
                                "platform": "CPU",
                            }, ostream=OutputStream(None))
        assert trajectory is not None
        frames = trajectory.frames
        d = np.linalg.norm(frames[:, 1:7] - frames[:, :1], axis=-1) * 10.0
        assert np.all(np.abs(d - 0.74) < 0.2), d.round(3)

    def test_harmonic_is_the_default(self):
        from veloxchem.autocdam import FrameworkSystem
        assert AutoCDAM().metal_dummy_bond == "harmonic"
        assert CDAMCandidate("Al", 3).metal_dummy_bond == "harmonic"
        assert CDAMCandidate("Al", 3).bond_force_constant == 2000.0
        assert FrameworkSystem.metal_dummy_bond == "harmonic"

    def test_framework_object_must_match_the_driver(self):
        from veloxchem import Molecule
        from veloxchem.autocdam import FrameworkSystem
        from veloxchem.errorhandler import VeloxChemError
        mol = Molecule.read_xyz_file(str(XYZ))
        mol.set_charge(-1)
        mof = FrameworkSystem(structure=mol, charges=CHARGES,
                              target_metal="Al", representation="cdam",
                              name="mof")
        fit = AutoCDAM(ostream=OutputStream(None))
        fit.metal_dummy_bond = "constraint"
        targets = {"al_o": {"property": "distance", "atoms": ["Al", "O3"],
                            "value": 1.93, "tolerance": 0.05}}
        with pytest.raises(VeloxChemError, match="set the same value"):
            fit.setup(CDAMCandidate("Al", 3), targets, systems={"mof": mof})

    def test_driver_sets_the_bond_time_step_and_frameworks(self):
        fit = AutoCDAM()
        fit.metal_dummy_bond = "harmonic"
        options = fit._calculator_options()
        assert options["md_options"]["timestep_fs"] == 1.0
        assert options["fep_options"]["timestep_fs"] == 1.0
        fit.md_options = {"timestep_fs": 0.5}
        assert fit._calculator_options()["md_options"]["timestep_fs"] == 0.5
        fit.md_options, fit.pilot = None, True
        assert "timestep_fs" not in fit._calculator_options().get(
            "md_options", {})
        fit.metal_dummy_bond = "constraint"
        fit.pilot = False
        assert "md_options" not in fit._calculator_options()

    def test_bond_force_constant_keyword_sets_the_candidate(self):
        fit = AutoCDAM()
        fit.update_settings({"metal_dummy_bond": "harmonic",
                             "bond_force_constant": 5000.0,
                             "cheap_first": False})
        fit.setup(CDAMCandidate("Al", 3), {
            "IOD": {"value": 1.9, "tolerance": 0.05, "role": "required"}})
        assert fit.candidate.bond_force_constant == 5000.0
        assert fit.candidate.metal_dummy_bond == "harmonic"
        assert "bond_force_constant" not in {
            s.name for s in fit.parameter_specs if s.optimizable}

    def test_distance_and_bond_constant_are_searched_when_none(self):
        fit = AutoCDAM()
        fit.update_settings({"metal_dummy_bond": "harmonic",
                             "cheap_first": False})
        fit.setup(CDAMCandidate("Al", 3), {
            "IOD": {"value": 1.9, "tolerance": 0.05, "role": "required"}})
        bounds = {s.name: tuple(s.bounds) for s in fit.parameter_specs
                  if s.optimizable}
        assert bounds["metal_dummy_distance"] == (0.70, 1.00)
        assert bounds["bond_force_constant"] == (1000.0, 3000.0)
        # constrained M-D: the bond constant has no effect and is not searched
        fit.metal_dummy_bond = "constraint"
        fit.setup(CDAMCandidate("Al", 3), {
            "IOD": {"value": 1.9, "tolerance": 0.05, "role": "required"}})
        assert "bond_force_constant" not in {
            s.name for s in fit.parameter_specs if s.optimizable}


class TestBondStabilityCheck:

    def test_minimum_bond_constant(self):
        # dummy 0.94 e at 0.71 A, hydroxide O -1.212 e at 1.80 A: collapsed at
        # k = 1950 in the k1950 fit (Al-O(hydroxide) 1.505 A)
        k = minimum_bond_force_constant(0.94, -1.212, 1.80, 0.71)
        assert k == pytest.approx(1972.0, rel=0.01)
        assert minimum_bond_force_constant(0.94, 0.4, 1.80, 0.71) == 0.0
        assert minimum_bond_force_constant(0.94, -1.0, 0.70, 0.71) == np.inf

    def test_node_fails_weak_bonds_before_any_md(self):
        node = MetalNode.read_xyz_file(XYZ, charge=-1, charges=CHARGES)
        targets = coerce_targets(
            expand_system_objects({"Al in MOF": {"system": node}}, None,
                                  "spce")[0])
        spec = node.system(water_model="spce")
        spec.update(type="framework", metal_dummy_bond="harmonic")
        router = build_systems({"node": spec}, targets, None, None)
        framework = router.systems["node"]
        strong = al("harmonic", dummy_charge=0.4, metal_dummy_distance=0.7,
                    bond_force_constant=1950.0)
        weak = al("harmonic", dummy_charge=0.94, metal_dummy_distance=0.71,
                  bond_force_constant=1950.0)
        assert framework.bond_check(strong) is None
        assert "too weak" in framework.bond_check(weak)
        result = router.calculate(weak.to_generated_batch(), None, targets)
        assert result.failed and "too weak" in result.failure_reason
        # constrained M-D: no check
        framework.metal_dummy_bond = "constraint"
        assert framework.bond_check(weak) is None

    def test_search_proposes_only_points_that_pass(self):
        fit = AutoCDAM()
        fit.update_settings({"metal_dummy_bond": "harmonic",
                             "bond_force_constant": 1950.0,
                             "cheap_first": False, "initial_samples": 9})
        node = MetalNode.read_xyz_file(XYZ, charge=-1, charges=CHARGES)
        fit.setup(CDAMCandidate("Al", 3, water_model="spce"), {
            "IOD": {"value": 1.9, "tolerance": 0.05, "role": "required"},
            "Al in MOF": {"system": node}})
        search = fit._make_fitting_driver({})
        design = search._design()
        assert len(design) == 9
        points = [search._from_unit(u) for u in design]
        assert all(fit._candidate_check(p) is None for p in points)
        # the plain Sobol design of the same size has failing points
        search.feasible = None
        plain = [search._from_unit(u) for u in search._design()]
        assert any(fit._candidate_check(p) for p in plain)
