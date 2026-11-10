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
"""Step 4: CDAM internal constraints and nonbonded MOF-node coordination.

The Al(OH)2(HCOO)2(H2O)2(-1) fragment has one metal, coordinating carboxylate
oxygens, two hydroxides, water, and a small linker fragment.
"""
import itertools
import json
from pathlib import Path

import numpy as np
import openmm as mm
import openmm.unit as unit
import pytest
import veloxchem as vlx

from veloxchem.autocdam import CDAMCandidate, FrameworkSystem
from veloxchem.autocdam.framework import compare_to_baseline
from veloxchem.autocdam.geometry import minimum_image

STRUCTURES = Path(__file__).parent / "data"
CHARGES = json.loads(
    (STRUCTURES /
     "autocdam_al_node_fragment_charges.json").read_text())["charges_e"]
HYDROXIDE_O, CARBOXYLATE_O = 1, 5
FIRST_SHELL_O = [1, 3, 5, 9, 13, 16]  # hydroxide, carboxylate, water O on Al


def molecule():
    mol = vlx.Molecule.read_xyz_file(
        str(STRUCTURES / "autocdam_al_node_fragment.xyz"))
    mol.set_charge(-1)
    return mol


@pytest.fixture(scope="module", params=["finite", "periodic"])
def node(request):
    box = np.eye(3) * 25.0 if request.param == "periodic" else None
    return FrameworkSystem(structure=molecule(),
                           charges=CHARGES,
                           target_metal="Al",
                           representation="cdam",
                           box_vectors=box)


def point(**values):
    settings = dict(metal_lj_sigma=2.0, metal_lj_epsilon=0.2, dummy_charge=0.5)
    settings.update(values)
    return CDAMCandidate("Al", 3, **settings)


def nonbonded(system):
    return next(
        f for f in system.getForces() if isinstance(f, mm.NonbondedForce))


def context(system, positions, box=None):
    integrator = mm.LangevinMiddleIntegrator(300 * unit.kelvin,
                                             5 / unit.picosecond,
                                             0.5 * unit.femtosecond)
    integrator.setRandomNumberSeed(3)
    ctx = mm.Context(system, integrator,
                     mm.Platform.getPlatformByName("Reference"))
    if box is not None:
        ctx.setPeriodicBoxVectors(*box)
    ctx.setPositions(positions)
    return ctx


def energy(system, positions):
    return context(
        system, positions).getState(getEnergy=True).getPotentialEnergy()._value


def expected_exceptions(baseline, site):
    """1-2/1-3/1-4 pairs of the covalent graph without metals, plus CDAM internals."""
    neighbours = {i: set() for i in range(len(baseline.elements))}
    for left, right in baseline.bonds:
        neighbours[left].add(right)
        neighbours[right].add(left)
    pairs = set()
    for start in neighbours:
        frontier, seen = {start}, {start}
        for _ in range(3):
            frontier = {n for atom in frontier for n in neighbours[atom]} - seen
            seen |= frontier
            pairs |= {tuple(sorted((start, other))) for other in frontier}
    members = [site["metal"], *site["dummies"]]
    pairs |= {
        tuple(sorted(pair)) for pair in itertools.combinations(members, 2)
    }
    return pairs


class TestAutoCdamFrameworkCdam:

    def test_site_charge_and_system_charge_are_conserved(self, node):
        for q in (0.2, 0.7):
            system, _, site = node.build(point(dummy_charge=q))
            force = nonbonded(system)
            charges = [
                force.getParticleParameters(i)[0]._value
                for i in range(system.getNumParticles())
            ]
            assert sum(charges) == pytest.approx(-1.0)
            assert charges[0] + 6 * q == pytest.approx(
                site["metal_site_charge"])
            assert site["metal_site_charge"] == pytest.approx(CHARGES[0])

    def test_dummies_point_at_the_first_coordination_shell(self, node):
        system, positions, site = node.build(point())
        shell = node.first_shell(6)
        vectors = positions[shell] - positions[0]
        for dummy in site["dummies"]:
            bond = positions[dummy] - positions[0]
            direction = bond / np.linalg.norm(bond)
            cosines = vectors @ direction / np.linalg.norm(vectors, axis=1)
            assert cosines.max() > 0.99
        assert sorted(shell) == [1, 3, 5, 9, 13, 16]

    def test_coordination_is_nonbonded_and_exclusions_are_exact(self, node):
        baseline = node.prepare_baseline()
        system, _, site = node.build(point())
        members = {site["metal"], *site["dummies"]}
        for force in system.getForces():
            if isinstance(force, mm.HarmonicBondForce):
                for i in range(force.getNumBonds()):
                    assert not members & set(force.getBondParameters(i)[:2])
            elif isinstance(force, mm.PeriodicTorsionForce):
                for i in range(force.getNumTorsions()):
                    assert not members & set(force.getTorsionParameters(i)[:4])
            elif isinstance(force, mm.HarmonicAngleForce):
                for i in range(force.getNumAngles()):
                    atoms = set(force.getAngleParameters(i)[:3])
                    assert not atoms & members or atoms <= members  # only CDAM internal angles
        force = nonbonded(system)
        exceptions = {
            tuple(sorted(force.getExceptionParameters(i)[:2]))
            for i in range(force.getNumExceptions())
        }
        assert exceptions == expected_exceptions(baseline, site)
        for ligand in (HYDROXIDE_O, CARBOXYLATE_O):
            assert all(
                tuple(sorted((m, ligand))) not in exceptions for m in members)

    def test_internal_cdam_interactions_are_not_double_counted(self, node):
        system, _, site = node.build(point())
        force = nonbonded(system)
        members = {site["metal"], *site["dummies"]}
        for i in range(force.getNumExceptions()):
            left, right, charge_product, _, epsilon = force.getExceptionParameters(
                i)
            if {left, right} <= members:
                assert charge_product._value == 0.0 and epsilon._value == 0.0
        constrained = {
            tuple(sorted(system.getConstraintParameters(i)[:2]))
            for i in range(system.getNumConstraints())
        }
        assert constrained == {(0, d) for d in site["dummies"]}
        assert not any(
            f.getName() == "CDAM_metal_dummy" for f in system.getForces())

    def test_dummy_oxygen_electrostatics_and_parameters_change_the_energy(
            self, node):
        _, positions, _ = node.build(point())
        energies = {}
        for name, values in {
                "ref": {},
                "q": {
                    "dummy_charge": 0.3
                },
                "sigma": {
                    "metal_lj_sigma": 2.4
                },
                "eps": {
                    "metal_lj_epsilon": 0.4
                }
        }.items():
            system, _, _ = node.build(point(**values))
            energies[name] = energy(system, positions)
        assert all(np.isfinite(list(energies.values())))
        assert len({round(value, 6) for value in energies.values()}) == 4
        # Dummies interact with the hydroxide oxygen: removing its charge changes the energy.
        system, _, site = node.build(point())
        force = nonbonded(system)
        q, sigma, eps = force.getParticleParameters(HYDROXIDE_O)
        before = energy(system, positions)
        force.setParticleParameters(HYDROXIDE_O, 0.0, sigma, eps)
        assert energy(system, positions) != pytest.approx(before)

    def test_non_target_parameters_stay_at_baseline(self, node):
        system, _, site = node.build(
            point(dummy_charge=0.65, metal_lj_sigma=1.7))
        assert compare_to_baseline(system, node.prepare_baseline(),
                                   changed={0}) == []

    def test_md_keeps_distances_geometry_and_covalent_bonds(self, node):
        baseline = node.prepare_baseline()
        system, positions, site = node.build(point(metal_dummy_distance=0.9))
        box = None if baseline.box_vectors_nm is None else baseline.box_vectors_nm * unit.nanometer
        ctx = context(system, positions, box)
        mm.LocalEnergyMinimizer.minimize(ctx, 10.0, 500)
        dummies = site["dummies"]
        reference = np.array(
            [[d1, d2] for d1, d2 in itertools.combinations(range(6), 2)])
        angle_errors, energies = [], []
        for _ in range(20):
            ctx.getIntegrator().step(50)
            state = ctx.getState(getPositions=True,
                                 getEnergy=True,
                                 getForces=True)
            xyz = state.getPositions(asNumpy=True)._value
            energies.append(state.getPotentialEnergy()._value)
            assert np.isfinite(state.getForces(asNumpy=True)._value).all()
            vectors = minimum_image(xyz[dummies] - xyz[0],
                                    baseline.box_vectors_nm)
            assert np.linalg.norm(vectors, axis=1) == pytest.approx(0.09,
                                                                    abs=1e-5)
            units = vectors / np.linalg.norm(vectors, axis=1)[:, None]
            angles = np.degrees(
                np.arccos(
                    np.clip(
                        np.sum(units[reference[:, 0]] * units[reference[:, 1]],
                               axis=1), -1, 1)))
            ideal = np.array([90.0 if a < 135 else 180.0 for a in angles])
            angle_errors.append(np.sqrt(np.mean((angles - ideal)**2)))
            ligand_o = minimum_image(xyz[FIRST_SHELL_O] - xyz[0],
                                     baseline.box_vectors_nm)
            # the nonbonded first shell stays on the metal
            assert np.linalg.norm(ligand_o, axis=1).max() < 0.25
            for left, right in baseline.bonds:  # every covalent bond survives
                length = np.linalg.norm(
                    minimum_image(xyz[right] - xyz[left],
                                  baseline.box_vectors_nm))
                assert 0.08 < length < 0.16
        assert max(angle_errors) < 12.0
        assert np.isfinite(energies).all() and np.ptp(energies) < 500.0

    def test_metal_dummy_distance_change_rebuilds_constraints(self, node):
        for distance in (0.8, 1.0):
            system, positions, site = node.build(
                point(metal_dummy_distance=distance))
            lengths = [
                system.getConstraintParameters(i)[2]._value
                for i in range(system.getNumConstraints())
            ]
            assert lengths == pytest.approx([distance * 0.1] * 6)
            assert np.linalg.norm(positions[site["dummies"][0]] -
                                  positions[0]) == pytest.approx(distance * 0.1)

    def test_harmonic_metal_dummy_bond_is_explicit_and_flexible(self):
        node = FrameworkSystem(structure=molecule(),
                               charges=CHARGES,
                               target_metal="Al",
                               representation="cdam",
                               metal_dummy_bond="harmonic")
        system, _, site = node.build(point())
        assert system.getNumConstraints() == 0
        bonds = next(
            f for f in system.getForces() if f.getName() == "CDAM_metal_dummy")
        assert bonds.getNumBonds() == 6

    def test_c4_acts_only_between_metal_and_water_oxygens(self):
        node = FrameworkSystem(structure=molecule(),
                               charges=CHARGES,
                               target_metal="Al",
                               representation="cdam")
        assert node.c4_oxygens() == {13, 16}
        candidate = point(c4_metal_water=10.0, potential_form="lj12_6_4")
        system, positions, site = node.build(candidate)
        force = next(
            f for f in system.getForces() if f.getName() == "CDAM_1264")
        assert force.getNumInteractionGroups() == 1
        group = force.getInteractionGroupParameters(0)
        assert set(group[0]) == {0} and set(group[1]) == {13, 16}
        no_c4, _, _ = node.build(point())
        assert energy(system,
                      positions) < energy(no_c4,
                                          positions)  # attractive -C4/r^4
