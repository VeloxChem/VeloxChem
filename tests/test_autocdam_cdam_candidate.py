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

import math

import pytest

from veloxchem.autocdam import CDAMCandidate


def _zinc_candidate(**overrides) -> CDAMCandidate:
    values = {
        "metal": "Zn",
        "formal_charge": 2,
        "metal_lj_sigma": 2.0,
        "metal_lj_epsilon": 0.2,
    }
    values.update(overrides)
    return CDAMCandidate(**values)


class TestAutoCdamCdamCandidate:

    def test_octahedral_candidate_derives_the_central_charge_and_geometry(
            self) -> None:
        candidate = _zinc_candidate(dummy_charge=0.5, metal_dummy_distance=0.9)

        assert candidate.coordination_number == 6
        assert candidate.central_charge == -1.0
        assert candidate.coordinates() == [
            [0.0, 0.0, 0.0],
            [0.9, 0.0, 0.0],
            [-0.9, 0.0, 0.0],
            [0.0, 0.9, 0.0],
            [0.0, -0.9, 0.0],
            [0.0, 0.0, 0.9],
            [0.0, 0.0, -0.9],
        ]

        batch = candidate.to_generated_batch()
        assert batch.n_configurations == 1
        assert len(batch.system.atoms) == 7
        assert sum(batch.system.charges) == pytest.approx(2.0)
        assert batch.metadata["candidate_id"] == candidate.candidate_id

    @pytest.mark.parametrize("geometry",
                             ["square_antiprismatic", "dodecahedral"])
    def test_zirconium_eight_coordinate_geometries_have_a_shared_dummy_radius(
            self, geometry) -> None:
        candidate = CDAMCandidate(
            metal="Zr",
            formal_charge=4,
            geometry=geometry,
            dummy_charge=0.55,
            metal_dummy_distance=0.90,
            metal_lj_sigma=2.0,
            metal_lj_epsilon=0.20,
        )

        coordinates = candidate.coordinates()
        batch = candidate.to_generated_batch()

        assert candidate.coordination_number == 8
        assert candidate.central_charge == pytest.approx(-0.40)
        assert len(coordinates) == 9
        assert all(
            math.dist(coordinates[0], coordinate) == pytest.approx(0.90)
            for coordinate in coordinates[1:])
        assert len(batch.system.atoms) == 9
        assert sum(batch.system.charges) == pytest.approx(4.0)

    def test_candidate_writes_pdb_and_openmm_xml_with_unit_conversions(
            self, tmp_path) -> None:
        candidate = _zinc_candidate(
            potential_form="lj12_6_4",
            c4_metal_water=10.0,
            water_model="TIP3P",
        )

        pdb_path, xml_path = candidate.write_files(tmp_path)
        pdb = pdb_path.read_text(encoding="utf-8")
        xml = xml_path.read_text(encoding="utf-8")

        assert pdb_path.name == f"{candidate.candidate_id}.pdb"
        assert pdb.count("HETATM") == 7
        assert "   0.900   0.000   0.000" in pdb
        assert xml.count("<Type ") == 7
        assert xml.count("<Angle ") == 15
        assert 'type="CDAM_Zn" charge="-1" sigma="0.2" epsilon="0.8368"' in xml
        assert '<AutoCDAM potential_form="lj12_6_4" water_model="TIP3P" c4_metal_water="0.004184"' in xml
        assert 'lj_mixing_rule="lorentz_berthelot"' in xml
        assert "sigma_ij=(sigma_i+sigma_j)/2" in xml
        assert '<NonbondedForce coulomb14scale="0.833333" lj14scale="0.5">' in xml

    def test_candidate_rejects_geometry_and_potential_inconsistencies(
            self) -> None:
        with pytest.raises(ValueError, match="coordination_number=6"):
            _zinc_candidate(coordination_number=4)
        with pytest.raises(ValueError, match="nonnegative"):
            _zinc_candidate(dummy_charge=-0.01)
        with pytest.raises(ValueError, match="requires c4"):
            _zinc_candidate(c4_metal_water=1.0)
        with pytest.raises(ValueError, match="requires metal_lj"):
            CDAMCandidate(metal="Zn", formal_charge=2,
                          metal_lj_sigma=None).write_files("unused")

    def test_lj12_6_artifacts_are_loadable_by_openmm_when_available(
            self, tmp_path) -> None:
        app = pytest.importorskip("openmm.app")
        candidate = _zinc_candidate()
        pdb_path, xml_path = candidate.write_files(tmp_path)

        force_field = app.ForceField(str(xml_path))
        pdb = app.PDBFile(str(pdb_path))
        system = force_field.createSystem(pdb.topology)

        assert system.getNumParticles() == 7
        force_types = {
            type(system.getForce(index)).__name__
            for index in range(system.getNumForces())
        }
        assert {"HarmonicBondForce", "HarmonicAngleForce",
                "NonbondedForce"} <= force_types

    def test_large_dummy_charge_preserves_total_ion_charge(self) -> None:
        import xml.etree.ElementTree as ET

        candidate = _zinc_candidate(dummy_charge=1.5)
        batch = candidate.to_generated_batch()
        assert candidate.central_charge == pytest.approx(-7.0)
        assert sum(batch.system.charges) == pytest.approx(2.0)
        atoms = ET.fromstring(
            candidate.to_openmm_xml()).findall("NonbondedForce/Atom")
        assert all(float(atom.attrib["sigma"]) == 0.0 for atom in atoms[1:])
        assert all(float(atom.attrib["epsilon"]) == 0.0 for atom in atoms[1:])

    def test_tip3p_candidate_xml_can_be_combined_with_water_forcefield(
            self, tmp_path) -> None:
        app = pytest.importorskip("openmm.app")
        candidate = _zinc_candidate(water_model="TIP3P")
        _, xml_path = candidate.write_files(tmp_path)

        # Both NonbondedForce sections must agree on 1-4 scales before OpenMM can
        # apply its built-in Lorentz-Berthelot mixing to unlike atom types.
        forcefield = app.ForceField(str(xml_path), "tip3p.xml")
        assert forcefield is not None

    @pytest.mark.parametrize("water_model", ["tip3p", "spce", "tip4pew"])
    def test_exported_xml_uses_openmm_default_mixing_for_metal_water(
            self, tmp_path, water_model):
        """Check actual unlike-pair energies after loading the exported force field.

        Removing charges isolates LJ: only the central metal and water oxygen have
        nonzero epsilon. The expected energy uses Lorentz-Berthelot, without any
        custom mixing expression or AutoCDAM metadata parser.
        """
        import openmm as mm
        import openmm.app as app
        import openmm.unit as unit

        candidate = _zinc_candidate(water_model=water_model)
        pdb_path, xml_path = candidate.write_files(tmp_path)
        pdb = app.PDBFile(str(pdb_path))
        forcefield = app.ForceField(str(xml_path), f"{water_model}.xml")
        modeller = app.Modeller(pdb.topology, pdb.positions)

        water = app.Topology()
        residue = water.addResidue("HOH", water.addChain())
        oxygen = water.addAtom("O", app.element.oxygen, residue)
        h1 = water.addAtom("H1", app.element.hydrogen, residue)
        h2 = water.addAtom("H2", app.element.hydrogen, residue)
        water.addBond(oxygen, h1)
        water.addBond(oxygen, h2)
        water_modeller = app.Modeller(water, [
            mm.Vec3(0.30, 0, 0),
            mm.Vec3(0.36, 0.07, 0),
            mm.Vec3(0.36, -0.07, 0)
        ] * unit.nanometer)
        # Add TIP4P-Ew's virtual charge site before combining with the CDAM sites,
        # which are explicit bonded atoms rather than water virtual sites.
        water_modeller.addExtraParticles(forcefield)
        modeller.add(water_modeller.topology, water_modeller.positions)
        system = forcefield.createSystem(modeller.topology,
                                         nonbondedMethod=app.NoCutoff)
        nonbonded = next(force for force in system.getForces()
                         if isinstance(force, mm.NonbondedForce))
        nonbonded.setForceGroup(1)
        for index in range(nonbonded.getNumParticles()):
            _, sigma, epsilon = nonbonded.getParticleParameters(index)
            nonbonded.setParticleParameters(index, 0, sigma, epsilon)

        oxygen_index = next(atom.index
                            for atom in modeller.topology.atoms()
                            if atom.residue.name == "HOH" and atom.name == "O")
        _, sigma_m, epsilon_m = nonbonded.getParticleParameters(0)
        _, sigma_o, epsilon_o = nonbonded.getParticleParameters(oxygen_index)
        sigma = (sigma_m + sigma_o).value_in_unit(unit.nanometer) / 2
        epsilon = math.sqrt(
            epsilon_m.value_in_unit(unit.kilojoule_per_mole) *
            epsilon_o.value_in_unit(unit.kilojoule_per_mole))
        r = 0.30
        expected = 4 * epsilon * ((sigma / r)**12 - (sigma / r)**6)

        integrator = mm.VerletIntegrator(0.001)
        context = mm.Context(system, integrator,
                             mm.Platform.getPlatformByName("Reference"))
        try:
            context.setPositions(modeller.positions)
            context.computeVirtualSites()
            actual = context.getState(getEnergy=True, groups={
                1
            }).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
            assert actual == pytest.approx(expected, rel=1e-10, abs=1e-10)
        finally:
            del context, integrator
