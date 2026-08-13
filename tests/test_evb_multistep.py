"""Multi-step EVB: a reaction path of N states run as N-1 chained FEPs.

The point of the feature is that the potential energy surfaces E_1 ... E_N are
one and the same quantity in every step, so most of what is checked here is
consistency across steps rather than any particular number: the states share one
atom ordering, the integration Hamiltonian agrees where two steps meet, and a
pair that reacts anywhere along the path is described the same way everywhere.

The three-state fixture is the PETase mechanism (int1 -> int2 -> int3); the
cheap end-to-end run uses a synthetic path built from the ethanol fixture.
"""

import copy
import csv
import math
from pathlib import Path

import numpy as np
import pytest

from veloxchem.evbdriver import EvbDriver
from veloxchem.outputstream import OutputStream
from veloxchem.reaffbuilder import ReactionForceFieldBuilder
from veloxchem.reactionsystembuilder import (ReactionSystemBuilder,
                                             EvbForceGroup)

from test_evb_helper import EvbTestHelper

pytest.importorskip('openmm')

import openmm as mm
import openmm.app as mmapp


class MultiStepFixtures:
    """Shared setup for the multi-step tests."""

    @staticmethod
    def _summarise(states):
        """Forming / breaking bonds per step and the path-global reactive set."""
        builder = ReactionForceFieldBuilder(ostream=OutputStream(None))
        forming, breaking = [], []
        for left, right in zip(states[:-1], states[1:]):
            formed, broken = builder._summarise_reaction(left, right)
            forming.append(formed)
            breaking.append(broken)
        reactive = set()
        for formed, broken in zip(forming, breaking):
            reactive |= set(formed) | set(broken)
        registry = builder._build_reactive_bond_registry(states, reactive)
        return forming, breaking, reactive, registry

    @staticmethod
    def _seed_driver(states):
        """An EvbDriver primed as if build_ffs_from_molecules had run.

        The fixtures hold ready-made force fields, so the (QM-free but slow)
        matching step is skipped.
        """
        forming, breaking, reactive, registry = MultiStepFixtures._summarise(states)

        EVB = EvbDriver(ostream=OutputStream(None))
        EVB.water_model = "cspce"
        EVB.states = states
        EVB.forming_bonds = forming
        EVB.breaking_bonds = breaking
        EVB.reactive_bonds = reactive
        EVB.bond_registry = registry
        return EVB

    @staticmethod
    def _synthetic_three_states():
        """[ethanol, ethanol, ethene+water]: an identity step followed by the real
        reaction. Cheap, and every reactive pair is either bonded on both sides of
        the first step or on neither, which is what only a longer path can produce.
        """
        pair = EvbTestHelper.evb_ff_pair()
        return [pair.reactant, copy.deepcopy(pair.reactant), pair.product]

    @staticmethod
    def _energy(topology, system, positions_nm):
        simulation = mmapp.Simulation(topology, system,
                                      mm.VerletIntegrator(0.001),
                                      mm.Platform.getPlatformByName("CPU"))
        simulation.context.setPositions(positions_nm)
        energy = simulation.context.getState(
            getEnergy=True).getPotentialEnergy().value_in_unit(
                mm.unit.kilojoules_per_mole)
        del simulation
        return energy

    @staticmethod
    def _morse_pairs(system):
        """The atom pairs carrying a morse bond in a PES system."""
        pairs = set()
        for force in system.getForces():
            if force.getForceGroup() != EvbForceGroup.REA_MORSE_BOND.value:
                continue
            for i in range(force.getNumBonds()):
                p1, p2, _ = force.getBondParameters(i)
                pairs.add(tuple(sorted((p1, p2))))
        return pairs

class TestForceFieldPath:

    def test_states_share_one_atom_ordering(self):
        # The invariant the whole feature rests on: two independently built
        # steps would carry the same atoms in a different sequence, and nothing
        # downstream would notice.
        states = EvbTestHelper.evb_petase_ms_states()
        ReactionSystemBuilder.assert_consistent_atom_ordering(states, "test")

        reference = states[0].molecule.get_element_ids()
        for state in states[1:]:
            assert list(state.molecule.get_element_ids()) == list(reference)

    def test_permuted_state_is_rejected(self):
        states = EvbTestHelper.evb_petase_ms_states()
        molecule = states[1].molecule
        elements = list(molecule.get_element_ids())
        # Find two atoms of different elements and swap their coordinates, which
        # is what an independently built state amounts to.
        first = 0
        second = next(i for i, e in enumerate(elements) if e != elements[0])
        coords = molecule.get_coordinates_in_angstrom()
        from veloxchem.molecule import Molecule
        from veloxchem.veloxchemlib import Point
        permuted = Molecule()
        order = list(range(len(elements)))
        order[first], order[second] = order[second], order[first]
        for i in order:
            permuted.add_atom(int(elements[i]), Point(coords[i]), 'angstrom')
        states[1].molecule = permuted

        with pytest.raises(Exception, match="Atom ordering mismatch"):
            ReactionSystemBuilder.assert_consistent_atom_ordering(states, "test")

    def test_two_state_wrapper_marshals_arguments(self):
        # build_force_fields is a thin wrapper over build_many_force_fields; the
        # only thing that can go wrong is the packing, so check that the tuple
        # it unpacks to lines up with the dict the general builder returns.
        from veloxchem.molecule import Molecule
        from test_evb_helper import EvbTestHelper

        data_dir = EvbTestHelper.evb_data_dir()
        reactant = Molecule.read_xyz_file(str(data_dir /
                                              "evb_ethanol_input.xyz"))
        product = [
            Molecule.read_xyz_file(str(data_dir / "evb_ethene_input.xyz")),
            Molecule.read_xyz_file(str(data_dir / "evb_water_input.xyz")),
        ]

        def _builder():
            builder = ReactionForceFieldBuilder(ostream=OutputStream(None))
            builder.calculate_resp = False
            builder.optimize_ff = False
            return builder

        rea_ff, pro_ff, forming, breaking, rea_ffs, pro_ffs, mapping = (
            _builder().build_force_fields(reactant=reactant, product=product))
        results = _builder().build_many_force_fields(states=[reactant, product])

        assert len(results['states']) == 2
        assert len(results['forming_bonds']) == 1
        assert len(results['mappings']) == 1
        assert forming == results['forming_bonds'][0]
        assert breaking == results['breaking_bonds'][0]
        assert mapping == results['mappings'][0]
        assert len(rea_ffs) == len(results['fragments'][0])
        assert len(pro_ffs) == len(results['fragments'][1])
        assert set(rea_ff.bonds) == set(results['states'][0].bonds)
        assert set(pro_ff.bonds) == set(results['states'][1].bonds)
        assert results['reactive_bonds'] == (set(forming) | set(breaking))

    def test_registry_averages_force_constants(self):
        states = EvbTestHelper.evb_petase_ms_states()
        _, _, reactive, registry = MultiStepFixtures._summarise(states)

        assert set(registry) == reactive
        for pair, entry in registry.items():
            bonded = [k for k, s in enumerate(states) if pair in s.bonds]
            assert entry['states_bonded'] == bonded
            assert bonded, f"{pair} is bonded in no state"
            expected = np.mean(
                [states[k].bonds[pair]['force_constant'] for k in bonded])
            assert entry['fc_ref'] == pytest.approx(expected)

    def test_registry_covers_every_step(self):
        # The PETase path is the interesting case: one pair forms in step 1 and
        # breaks again in step 2, another is bonded on both sides of step 1 and
        # breaks in step 2, and a third is bonded on neither side of step 1.
        states = EvbTestHelper.evb_petase_ms_states()
        forming, breaking, reactive, _ = MultiStepFixtures._summarise(states)

        assert len(forming) == 2 and len(breaking) == 2
        bonded_in_both = [
            pair for pair in reactive
            if pair in states[0].bonds and pair in states[1].bonds
        ]
        bonded_in_neither = [
            pair for pair in reactive
            if pair not in states[0].bonds and pair not in states[1].bonds
        ]
        assert bonded_in_both, "fixture no longer covers the bonded-in-both case"
        assert bonded_in_neither, (
            "fixture no longer covers the bonded-in-neither case")


class TestOneStepIsJustAShortPath:
    """A reactant -> product reaction goes through the same code as a longer
    path, and must still look exactly like it always did from the outside."""

    def test_reactant_and_product_are_views_on_states(self):
        states = MultiStepFixtures._synthetic_three_states()
        EVB = MultiStepFixtures._seed_driver(states)

        assert EVB.reactant is states[0]
        assert EVB.product is states[-1]

        # Read-only: the path is the single source of truth.
        with pytest.raises(AttributeError):
            EVB.reactant = states[1]

    def test_bond_lists_are_per_step_whatever_the_length(self):
        three = MultiStepFixtures._synthetic_three_states()
        for states in (three[:2], three):
            EVB = MultiStepFixtures._seed_driver(states)
            assert isinstance(EVB.forming_bonds, list)
            assert isinstance(EVB.breaking_bonds, list)
            assert len(EVB.forming_bonds) == len(states) - 1
            assert len(EVB.breaking_bonds) == len(states) - 1
            for bonds in EVB.forming_bonds + EVB.breaking_bonds:
                assert isinstance(bonds, set)

    @pytest.mark.timeconsuming
    def test_one_step_keeps_the_two_state_layout(self):
        pair = EvbTestHelper.evb_ff_pair()
        states = [pair.reactant, pair.product]
        config = {
            "name": "one_step",
            "temperature": 300.0,
            "sample_steps": 20,
            "write_step": 10,
            "equil_NVT_steps": 2,
            "equil_NPT_steps": 0,
            "initial_equil_NVT_steps": 0,
            "initial_equil_NPT_steps": 0,
            "skip_initial_equil": True,
            "step_size": 0.001,
            "n_replicas": 1,
            "minimize_every_lambda": True,
        }

        with EvbTestHelper.evb_chdir_tmp():
            EVB = MultiStepFixtures._seed_driver(states)
            EVB.build_systems(configurations=[config], Lambda=[0.0, 0.5, 1.0])
            EVB.run_FEP()
            conf = EVB.system_confs[0]

            # One step, and it writes into the configuration folder itself
            # rather than a step_00 subfolder.
            assert len(conf['steps']) == 1
            root = Path(conf['data_folder'])
            assert conf['steps'][0]['data_folder'] == str(root)
            assert not (root / "step_00").exists()

            # Unprefixed lambda systems and the reactant / product surfaces.
            run_folder = root / "run"
            assert (run_folder / "0.000_sys.xml").is_file()
            assert (run_folder / "reactant_sys.xml").is_file()
            assert (run_folder / "product_sys.xml").is_file()
            assert not any(run_folder.glob("step0_*"))
            assert not any(run_folder.glob("state_*"))

            # And the Energies.csv layout the analysis code still expects.
            with open(root / "Energies.csv", newline="") as handle:
                header = [cell.strip() for cell in next(csv.reader(handle))]
            assert header == [
                "Lambda", "reactant PES", "product PES", "reactant integration",
                "product integration", "Em", "replica", "direction"
            ]


class TestMultiStepSystems:

    def _build(self, states, configuration=None):
        forming, breaking, reactive, registry = MultiStepFixtures._summarise(states)
        builder = ReactionSystemBuilder(ostream=OutputStream(None))
        builder.water_model = "cspce"
        configuration = configuration or {
            "name": "ms",
            "temperature": 300.0,
        }
        step_systems, pes_systems, topology, positions = (
            builder.build_multi_step_systems(
                states,
                [0.0, 0.5, 1.0],
                configuration,
                reactive_bonds=reactive,
                bond_registry=registry,
            ))
        return step_systems, pes_systems, topology, positions

    def test_junction_hamiltonians_agree(self):
        states = EvbTestHelper.evb_petase_ms_states()
        step_systems, _, topology, positions = self._build(states)
        positions_nm = positions * 0.1

        assert len(step_systems) == len(states) - 1
        for step in range(len(step_systems) - 1):
            left = MultiStepFixtures._energy(topology, step_systems[step][1.0], positions_nm)
            right = MultiStepFixtures._energy(topology, step_systems[step + 1][0.0], positions_nm)
            assert left == pytest.approx(right, abs=1e-6), (
                f"integration Hamiltonian jumps between step {step} (l=1) and "
                f"step {step + 1} (l=0)")

    def test_one_pes_system_per_state(self):
        states = EvbTestHelper.evb_petase_ms_states()
        _, pes_systems, _, _ = self._build(states)
        assert sorted(pes_systems) == list(range(len(states)))

    def test_morse_bonds_follow_the_path_global_set(self):
        # A pair that reacts anywhere along the path is a morse bond in every
        # state that has it bonded - including a step where it does not react,
        # which is exactly what makes E_k independent of the step it comes from.
        states = EvbTestHelper.evb_petase_ms_states()
        _, _, reactive, _ = MultiStepFixtures._summarise(states)
        _, pes_systems, _, _ = self._build(states)

        for index, state in enumerate(states):
            expected = {pair for pair in reactive if pair in state.bonds}
            assert MultiStepFixtures._morse_pairs(pes_systems[index]) == expected, (
                f"state {index} does not describe the reacting pairs it has "
                "bonded as morse bonds")

    def test_identical_states_give_identical_energies(self):
        states = MultiStepFixtures._synthetic_three_states()
        _, pes_systems, topology, positions = self._build(states)
        positions_nm = positions * 0.1

        # States 0 and 1 are the same molecule, so the surfaces must agree.
        first = MultiStepFixtures._energy(topology, pes_systems[0], positions_nm)
        second = MultiStepFixtures._energy(topology, pes_systems[1], positions_nm)
        assert first == pytest.approx(second, abs=1e-9)

    def test_solvation_happens_once(self):
        states = MultiStepFixtures._synthetic_three_states()
        configuration = {
            "name": "ms_water",
            "temperature": 300.0,
            "solvent": "cspce",
            "padding": 0.6,
            "pressure": 1.0,
        }
        with EvbTestHelper.evb_chdir_tmp():
            step_systems, pes_systems, topology, _ = self._build(
                states, configuration)

        counts = {system.getNumParticles() for step in step_systems
                  for system in step.values()}
        counts |= {system.getNumParticles() for system in pes_systems.values()}
        assert len(counts) == 1, (
            f"steps disagree on the number of particles: {counts}")
        assert topology.getNumAtoms() == counts.pop()

        residues = [r for r in topology.residues() if r.name == "SOL"]
        assert residues, "no solvent was added"

    def test_cnt_is_refused(self):
        states = MultiStepFixtures._synthetic_three_states()
        with pytest.raises(Exception, match="not.*supported"):
            self._build(states, {
                "name": "ms_cnt",
                "temperature": 300.0,
                "CNT": True,
            })


@pytest.mark.timeconsuming
class TestMultiStepFep:

    def _tiny_config(self):
        return {
            "name": "ms_tiny",
            "temperature": 300.0,
            "sample_steps": 20,
            "write_step": 10,
            "equil_NVT_steps": 2,
            "equil_NPT_steps": 0,
            "initial_equil_NVT_steps": 0,
            "initial_equil_NPT_steps": 0,
            "skip_initial_equil": True,
            "step_size": 0.001,
            "n_replicas": 1,
            # Two equilibration steps are not enough to keep a Langevin run on
            # an unrelaxed geometry from occasionally blowing up to NaN; the
            # minimisation costs nothing for nine atoms and makes the run
            # reliable without changing what is being tested.
            "minimize_every_lambda": True,
        }

    def _run(self, states):
        EVB = MultiStepFixtures._seed_driver(states)
        EVB.build_systems(configurations=[self._tiny_config()],
                          Lambda=[0.0, 0.5, 1.0])
        EVB.run_FEP()
        return EVB, EVB.system_confs[0]

    def test_path_runs_and_chains(self):
        states = MultiStepFixtures._synthetic_three_states()
        with EvbTestHelper.evb_chdir_tmp():
            EVB, conf = self._run(states)

            assert len(conf['steps']) == len(states) - 1
            root = Path(conf['data_folder'])

            # One topology and one run folder for the whole path.
            assert (root / "topology.cif").is_file()
            run_folder = root / "run"
            for index in range(len(states)):
                assert (run_folder / f"state_{index}_sys.xml").is_file()
            for step in range(len(states) - 1):
                assert (run_folder / f"step{step}_0.000_sys.xml").is_file()

            for index, step_conf in enumerate(conf['steps']):
                folder = Path(step_conf['data_folder'])
                energies = folder / "Energies.csv"
                assert energies.is_file(), f"missing {energies}"

                with open(energies, newline="") as handle:
                    rows = [row for row in csv.reader(handle) if row]
                header = [cell.strip() for cell in rows[0]]
                assert header == [
                    "Lambda", "step", "state_0 PES", "state_1 PES",
                    "state_2 PES", "left integration", "right integration",
                    "Em", "replica", "direction"
                ]
                assert len(rows) > 1
                for row in rows[1:]:
                    values = [float(cell) for cell in row]
                    assert all(math.isfinite(value) for value in values)
                    assert int(values[1]) == index

                # Every step but the last hands its lambda=1 configuration to
                # the next one.
                if index + 1 < len(conf['steps']):
                    assert (folder / "final_state.xml").is_file()
                    assert conf['steps'][index + 1]['initial_state'] == str(
                        folder / "final_state.xml")
                    assert conf['steps'][index + 1]['skip_initial_equil']

    def test_adjacent_reporting_narrows_the_columns(self):
        states = MultiStepFixtures._synthetic_three_states()
        config = self._tiny_config()
        config['pes_states_to_report'] = 'adjacent'
        with EvbTestHelper.evb_chdir_tmp():
            EVB = MultiStepFixtures._seed_driver(states)
            EVB.build_systems(configurations=[config],
                              Lambda=[0.0, 0.5, 1.0])
            EVB.run_FEP()
            conf = EVB.system_confs[0]

            folder = Path(conf['steps'][1]['data_folder'])
            with open(folder / "Energies.csv", newline="") as handle:
                header = [cell.strip() for cell in next(csv.reader(handle))]
            assert header == [
                "Lambda", "step", "state_1 PES", "state_2 PES",
                "left integration", "right integration", "Em", "replica",
                "direction"
            ]

    def test_reload_rebuilds_the_steps(self):
        states = MultiStepFixtures._synthetic_three_states()
        with EvbTestHelper.evb_chdir_tmp():
            EVB, conf = self._run(states)
            data_folder = conf['data_folder']

            reloaded = EvbDriver(ostream=OutputStream(None))
            reloaded.load_initialisation(data_folder,
                                         load_systems=True,
                                         load_top=True)

            assert len(reloaded.system_confs) == 1
            loaded = reloaded.system_confs[0]
            assert loaded['n_states'] == len(states)
            assert len(loaded['steps']) == len(states) - 1
            assert reloaded.lambda_vec == EVB.lambda_vec

            for step, step_conf in enumerate(loaded['steps']):
                systems = step_conf['systems']
                # That step's lambda systems plus every state's PES system.
                for lam in reloaded.lambda_vec:
                    assert lam in systems, f"step {step} is missing lambda {lam}"
                for index in range(len(states)):
                    assert f"state_{index}" in systems

    def test_output_loads_back_per_step(self):
        states = MultiStepFixtures._synthetic_three_states()
        with EvbTestHelper.evb_chdir_tmp():
            EVB, conf = self._run(states)

            results = EVB._load_output_from_folders(1, False, 1)
            loaded = results["configuration_results"][conf["name"]]

            assert loaded["n_states"] == len(states)
            assert len(loaded["steps"]) == len(states) - 1

            for index, step in enumerate(loaded["steps"]):
                # Every state's surface is read back, not just the two
                # flanking this step.
                assert step["E_pes_names"].split(",") == [
                    f"state_{k}" for k in range(len(states))
                ]
                assert step["E_pes"].shape[0] == len(states)
                assert step["step"] == index
                assert np.all(step["Step_frame"] == index)

                # E1 / E2 are this step's own flanking states, picked by name.
                assert np.allclose(step["E1_pes"], step["E_pes"][index])
                assert np.allclose(step["E2_pes"], step["E_pes"][index + 1])

                for key in ("E1_int", "E2_int", "E_m_pes", "Ep", "Ek"):
                    assert np.all(np.isfinite(step[key]))
                assert len(step["E1_pes"]) == len(results["Lambda_frame"])

            # States 0 and 1 are the same molecule here, so the surfaces they
            # contribute must agree wherever both were evaluated.
            for step in loaded["steps"]:
                assert np.allclose(step["E_pes"][0], step["E_pes"][1])

    def test_energy_profiles_refuse_a_path(self):
        states = MultiStepFixtures._synthetic_three_states()
        with EvbTestHelper.evb_chdir_tmp():
            EVB, _ = self._run(states)
            with pytest.raises(Exception, match="cannot yet fit"):
                EVB.compute_energy_profiles(barrier=100.0, free_energy=0.0)

    def test_compute_force_groups_is_refused(self):
        states = MultiStepFixtures._synthetic_three_states()
        with EvbTestHelper.evb_chdir_tmp():
            EVB, _ = self._run(states)
            with pytest.raises(Exception, match="not available"):
                EVB.compute_force_groups()
