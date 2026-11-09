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
"""Check actual OpenMM forces/energies and object ownership without running MD."""
import math
import numpy as np
import openmm as mm
import openmm.app as app
import openmm.unit as unit
import pytest

from veloxchem.autocdam import Calculator, CDAMCandidate, CoordinateGenerator, OpenMMSystems, ParameterSet, FitTarget
from veloxchem.autocdam.md import run_md
from veloxchem.autocdam.models import convert_units


@pytest.fixture
def setup(ostream):
    candidate = CDAMCandidate("Al",
                              3,
                              water_model="tip3p",
                              metal_lj_sigma=2,
                              metal_lj_epsilon=0.2,
                              potential_form="lj12_6_4")
    generator = CoordinateGenerator(ostream=ostream)
    # Small deterministic solute + water system: no stochastic solvent packing.
    topology = app.Topology()
    residue = topology.addResidue("CDM", topology.addChain())
    metal = topology.addAtom("Al", app.element.aluminum, residue)
    for i in range(6):
        dummy = topology.addAtom(f"D{i+1}", None, residue)
        topology.addBond(metal, dummy)
    vacuum_topology = topology
    vacuum = mm.System()
    nb = mm.NonbondedForce()
    for charge, mass in zip([0] + [0.5] * 6, [8.9815385] + [3] * 6):
        vacuum.addParticle(mass)
        nb.addParticle(charge, 0.2, 0.0)
    nb.createExceptionsFromBonds([(0, i) for i in range(1, 7)], 0, 0)
    vacuum.addForce(nb)
    bonds = mm.HarmonicBondForce()
    for i in range(1, 7):
        bonds.addBond(0, i, 0.09, 1)
    vacuum.addForce(bonds)
    solvated_topology = app.Topology()
    chain = solvated_topology.addChain()
    solute = solvated_topology.addResidue("CDM", chain)
    for atom in topology.atoms():
        solvated_topology.addAtom(atom.name, atom.element, solute)
    sol_atoms = list(solvated_topology.atoms())
    for i in range(1, 7):
        solvated_topology.addBond(sol_atoms[0], sol_atoms[i])
    water = solvated_topology.addResidue("HOH", chain)
    solvated_topology.addAtom("O", app.element.oxygen, water)
    solvated = mm.XmlSerializer.deserialize(mm.XmlSerializer.serialize(vacuum))
    solvated.addParticle(16)
    next(f for f in solvated.getForces()
         if isinstance(f, mm.NonbondedForce)).addParticle(-0.834, 0.315, 0.636)
    positions = np.array(candidate.coordinates()) * 0.1 * unit.nanometer
    solv_positions = np.vstack((positions.value_in_unit(
        unit.nanometer), [[0.3, 0.01, 0.02]])) * unit.nanometer
    template = OpenMMSystems(vacuum, vacuum_topology, positions, solvated,
                             solvated_topology, solv_positions)
    return candidate, generator, template


def energy(system, positions):
    integrator = mm.VerletIntegrator(0.001)
    context = mm.Context(system, integrator,
                         mm.Platform.getPlatformByName("Reference"))
    context.setPositions(positions)
    result = context.getState(
        getEnergy=True).getPotentialEnergy().value_in_unit(
            unit.kilojoule_per_mole)
    del context, integrator
    return result


class TestAutoCdamOpenmmBatch:

    def test_batch_updates_both_systems_and_preserves_template(
            self, setup, tmp_path, monkeypatch):
        candidate, generator, template = setup
        monkeypatch.chdir(tmp_path)
        before = mm.XmlSerializer.serialize(template.solvated_system)
        proposals = [
            ParameterSet({
                "dummy_charge": 0.4,
                "metal_lj_sigma": 2.2,
                "metal_lj_epsilon": 0.3,
                "metal_dummy_distance": 1.0,
                "c4_metal_water": 2
            }),
            ParameterSet({
                "dummy_charge": 0.5,
                "metal_lj_sigma": 2.8,
                "metal_lj_epsilon": 0.1,
                "metal_dummy_distance": 0.8,
                "c4_metal_water": 5
            })
        ]
        batches = generator.generate_parameter_batch(candidate,
                                                     proposals,
                                                     openmm_template=template)
        energies = []
        for batch, proposal in zip(batches, proposals):
            for system in [
                    batch.openmm_systems.vacuum_system,
                    batch.openmm_systems.solvated_system
            ]:
                force = next(f for f in system.getForces()
                             if isinstance(f, mm.NonbondedForce))
                charge, sigma, epsilon = force.getParticleParameters(0)
                assert charge.value_in_unit(
                    unit.elementary_charge) == pytest.approx(
                        3 - 6 * proposal["dummy_charge"])
                assert sigma.value_in_unit(unit.nanometer) == pytest.approx(
                    proposal["metal_lj_sigma"] * 0.1)
                assert epsilon.value_in_unit(
                    unit.kilojoule_per_mole) == pytest.approx(
                        proposal["metal_lj_epsilon"] * 4.184)
                assert system.getNumConstraints() == 6
                assert system.getConstraintParameters(0)[2].value_in_unit(
                    unit.nanometer) == pytest.approx(
                        proposal["metal_dummy_distance"] * 0.1)
                c4 = next(
                    f for f in system.getForces() if f.getName() == "CDAM_1264")
                assert c4.getParticleParameters(0)[0] == pytest.approx(
                    proposal["c4_metal_water"] * 4.184e-4)
            systems = batch.openmm_systems
            assert np.linalg.norm(systems.vacuum_positions[1].value_in_unit(
                unit.nanometer)) == pytest.approx(
                    proposal["metal_dummy_distance"] * 0.1)
            energies.append(
                energy(systems.solvated_system, systems.solvated_positions))
        assert all(np.isfinite(energies))
        assert energies[0] != pytest.approx(energies[1])
        assert mm.XmlSerializer.serialize(template.solvated_system) == before
        assert list(tmp_path.iterdir()) == []
        # Repeated updates replace constraints and C4 instead of adding duplicates.
        regenerated = generator.update_openmm_systems(batches[0].openmm_systems,
                                                      batches[1])
        assert regenerated.vacuum_system.getNumConstraints() == 6
        assert sum(f.getName() == "CDAM_1264"
                   for f in regenerated.solvated_system.getForces()) == 1

    def test_calculator_consumes_objects_without_loading_files(
            self, setup, monkeypatch, tmp_path):
        candidate, generator, template = setup
        point = ParameterSet({"metal_lj_sigma": 2.5})
        batch = generator.generate_parameter_batch(candidate, [point],
                                                   openmm_template=template)[0]

        def unexpected(*args, **kwargs):
            raise AssertionError("object path must not read files")

        monkeypatch.setattr(
            "veloxchem.autocdam.candidate.CDAMCandidate.write_files",
            unexpected)
        monkeypatch.setattr(Calculator, "_initialize_platform",
                            lambda self: None)

        def fake_md(self, systems, md_targets, options=None):
            system = systems.solvated_system
            nb = next(f for f in system.getForces()
                      if isinstance(f, mm.NonbondedForce))
            assert nb.getParticleParameters(0)[1].value_in_unit(
                unit.nanometer) == pytest.approx(0.25)
            system.addForce(mm.MonteCarloBarostat(1, 298))
            return object(), {"IOD": 2.0}

        monkeypatch.setattr(Calculator, "_execute_md_block", fake_md)
        before = mm.XmlSerializer.serialize(
            batch.openmm_systems.solvated_system)
        calculator = Calculator()
        for repeat in range(2):
            result = calculator.calculate(
                batch,
                point, [FitTarget("IOD", 2, 0.1),
                        FitTarget("QM_ENERGY", 0, 1)],
                work_dir=tmp_path)
            assert result.values["IOD"] == 2
            assert np.isfinite(result.values["QM_ENERGY"])
            assert result.metadata["QM_ENERGY"]["cache_hit"] == bool(repeat)
        assert mm.XmlSerializer.serialize(
            batch.openmm_systems.solvated_system) == before
        assert list(tmp_path.iterdir()) == []
        with pytest.raises(ValueError, match="differ from the generated"):
            Calculator().calculate(batch, ParameterSet({"metal_lj_sigma": 9}),
                                   [FitTarget("IOD", 2, 0.1)])

    @pytest.mark.parametrize("water_model", ["tip3p", "spce", "tip4pew"])
    def test_builds_supported_templates_without_files(self, water_model,
                                                      tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        candidate = CDAMCandidate("Al",
                                  3,
                                  water_model=water_model,
                                  metal_lj_sigma=2,
                                  metal_lj_epsilon=0.2)
        template = CoordinateGenerator().create_openmm_template(candidate)
        assert template.vacuum_system.getNumParticles() == 7
        assert template.solvated_system.getNumParticles() > 7
        assert template.solvated_system.getNumConstraints() > 6
        assert template.vacuum_system.getNumConstraints() == 6
        assert list(tmp_path.iterdir()) == []

    def test_property_units_match_reference_units(self):
        values = {
            "HFE": -418.4,
            "IOD": 2.0,
            "CN": 6.0,
            "SHELL_ANGLE_RMSD": 180.0
        }
        metadata = {
            "HFE": {
                "uncertainty": 4.184
            },
            "IOD": {
                "spread": 0.1
            },
            "CN": {},
            "SHELL_ANGLE_RMSD": {}
        }
        convert_units(values, metadata, [
            FitTarget("HFE", -100, 1, unit="kcal/mol"),
            FitTarget("IOD", 0.2, 0.01, unit="nm"),
            FitTarget("CN", 6, 0.1),
            FitTarget("SHELL_ANGLE_RMSD", 0.0, 0.1, unit="rad"),
        ])
        assert values["HFE"] == pytest.approx(-100)
        assert values["IOD"] == pytest.approx(0.2)
        assert metadata["HFE"]["uncertainty"] == pytest.approx(1)
        assert metadata["IOD"]["spread"] == pytest.approx(0.01)
        assert values["SHELL_ANGLE_RMSD"] == pytest.approx(math.pi)

    def test_md_can_skip_equilibration_and_use_unlimited_minimization(self):
        """Run actual OpenMM to check the zero-stage path and saved production frames."""
        topology = app.Topology()
        residue = topology.addResidue("CDM", topology.addChain())
        topology.addAtom("Al", app.element.aluminum, residue)
        system = mm.System()
        system.addParticle(27)
        system.setDefaultPeriodicBoxVectors(mm.Vec3(2, 0, 0), mm.Vec3(0, 2, 0),
                                            mm.Vec3(0, 0, 2))
        positions = [[1, 1, 1]] * unit.nanometer
        trajectory = run_md(
            system, topology, positions, {
                "platform": "Reference",
                "nvt_steps": 0,
                "npt_steps": 0,
                "minimization_iterations": 0,
                "production_steps": 4,
                "report_interval": 2,
                "random_seed": 42
            })
        frames, boxes = trajectory.frames, trajectory.boxes
        assert trajectory.final_state is not None
        assert frames.shape == (2, 1, 3)
        assert boxes.shape == (2, 3, 3)
        assert np.isfinite(frames).all()
        assert not np.array_equal(frames[0], frames[1])
        assert not any(
            isinstance(force, mm.MonteCarloBarostat)
            for force in system.getForces())

    @pytest.mark.parametrize(
        "option", ["nvt_steps", "npt_steps", "minimization_iterations"])
    def test_md_rejects_negative_stage_lengths(self, option):
        with pytest.raises(ValueError, match="nonnegative integer"):
            run_md(None, None, None, {option: -1})

    def test_file_only_fep_adapter_preserves_proposal_and_requests_c4_scaling(
            self, setup, monkeypatch, tmp_path):
        from pathlib import Path
        import xml.etree.ElementTree as ET
        import veloxchem.autocdam.md as module
        candidate, generator, template = setup
        proposal = ParameterSet({
            "dummy_charge": 0.4,
            "metal_lj_sigma": 2.5,
            "c4_metal_water": 5.0
        })
        batch = generator.generate_parameter_batch(candidate, [proposal],
                                                   openmm_template=template)[0]
        seen = []

        class FileFep:

            def __init__(self, ostream):
                self.skip_vacuum = False
                self.coulomb_global_parameters = []

            def _scale_coulomb_global_parameters(self, system, lambda_val):
                raise AssertionError("not called by the adapter")

            def compute_solvation_from_openmm_files(self, **kwargs):
                seen.append(Path(kwargs["system_xml"]).parent)
                metal = ET.parse(
                    kwargs["solute_xml"]).getroot().find("NonbondedForce/Atom")
                assert float(metal.attrib["sigma"]) == pytest.approx(0.25)
                assert float(metal.attrib["charge"]) == pytest.approx(0.6)
                assert self.skip_vacuum
                assert self.coulomb_global_parameters == ["c4_lambda"]
                system = mm.XmlSerializer.deserialize(
                    Path(kwargs["system_xml"]).read_text())
                force = next(
                    f for f in system.getForces() if f.getName() == "CDAM_1264")
                assert force.getGlobalParameterName(0) == "c4_lambda"
                assert force.getParticleParameters(0)[0] == pytest.approx(
                    5 * 4.184e-4)
                return {"free_energy": -42.0}

        monkeypatch.setattr(module, "SolvationFepDriver", FileFep)
        monkeypatch.setattr(Calculator, "_initialize_platform",
                            lambda self: None)
        result = Calculator().calculate(
            batch,
            proposal, [FitTarget("HFE", -42, 1)],
            work_dir=tmp_path,
            fep_options={"equilibrate_first": False})
        assert result.metadata["HFE"]["quality_status"] == "unreviewed"
        assert seen and not seen[0].exists()
        assert list(tmp_path.iterdir()) == []

    def test_c4_acts_only_between_the_metal_and_water_oxygens(self, setup):
        """Regression: water O-O pairs must not feel the C4 term (solvent collapse)."""
        candidate, generator, template = setup
        for c4 in (0.0, 200.0):
            batch = generator.generate_parameter_batch(
                candidate, [ParameterSet({"c4_metal_water": c4})],
                openmm_template=template)[0]
            systems = batch.openmm_systems
            system = mm.XmlSerializer.deserialize(
                mm.XmlSerializer.serialize(systems.solvated_system))
            force = next(
                f for f in system.getForces() if f.getName() == "CDAM_1264")
            assert force.getNumInteractionGroups() == 1
            for i, f in enumerate(system.getForces()):
                f.setForceGroup(1 if f.getName() == "CDAM_1264" else 0)
            context = mm.Context(system, mm.VerletIntegrator(0.001),
                                 mm.Platform.getPlatformByName("Reference"))
            context.setPositions(systems.solvated_positions)
            energy = context.getState(getEnergy=True, groups={
                1
            }).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
            positions = np.array(
                systems.solvated_positions.value_in_unit(unit.nanometer))
            oxygens = [
                a.index
                for a in systems.solvated_topology.atoms()
                if a.element is not None and a.element.symbol == "O"
            ]
            r = np.linalg.norm(positions[oxygens] - positions[0], axis=1)
            cutoff = force.getCutoffDistance().value_in_unit(unit.nanometer)
            expected = -c4 * 4.184e-4 * np.sum(1.0 / r[r < cutoff]**4)
            assert energy == pytest.approx(expected, abs=1e-6)

    def test_hfe_alone_is_equilibrated_by_md_before_fep(self, setup,
                                                        monkeypatch, tmp_path):
        """Regression: a monitor-only HFE must not start FEP from the minimized template."""
        candidate, generator, template = setup
        batch = generator.generate_parameter_batch(candidate, [ParameterSet()],
                                                   openmm_template=template)[0]
        calls = []

        def md_block(self, systems, md_targets, options=None):
            calls.append(md_targets)
            return None, {}  # stop here: only the order matters

        monkeypatch.setattr(Calculator, "_initialize_platform",
                            lambda self: None)
        monkeypatch.setattr(Calculator, "_execute_md_block", md_block)
        monkeypatch.setattr(
            "veloxchem.autocdam.calculator.run_fep",
            lambda *a, **k: pytest.fail("FEP ran before equilibration"))
        result = Calculator().calculate(batch,
                                        ParameterSet(),
                                        [FitTarget("HFE", -42, 1)],
                                        work_dir=tmp_path)
        assert calls == [[]] and result.failed
