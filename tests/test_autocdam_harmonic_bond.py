import numpy as np
import openmm as mm
import pytest

from veloxchem.autocdam import AutoCDAM, CDAMCandidate, CoordinateGenerator
from veloxchem.autocdam.md import run_md
from veloxchem.outputstream import OutputStream

# Al3+ of Golo, Ahlquist, Su, ACS Omega 10, 3801 (2025), Table 1:
# Rmin/2 0.6158 A, epsilon 29.9565 kcal/mol, dummy 0.73 e at 0.74 A,
# Kb 1950 kcal/mol/A^2, D-M-D 200 kcal/mol/rad^2, harmonic M-D bonds
PAPER_AL = {
    "metal_lj_sigma": 2 * 0.6158 / 2**(1 / 6),
    "metal_lj_epsilon": 29.9565,
    "dummy_charge": 0.73,
    "metal_dummy_distance": 0.74,
    "bond_force_constant": 1950.0,
    "angle_force_constant": 200.0,
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
        assert site[0][3]._value == pytest.approx(1950.0 * 418.4)

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
