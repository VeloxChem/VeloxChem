from mpi4py import MPI
import numpy as np
from pathlib import Path
import pytest

from veloxchem.molecule import Molecule
from veloxchem.mmforcefieldgenerator import MMForceFieldGenerator
from veloxchem.solvationbuilder import SolvationBuilder
from veloxchem.errorhandler import VeloxChemError

pytest.importorskip("rdkit")


class RecordingOutput:

    def __init__(self):
        self.infos = []

    def print_header(self, text):
        pass

    def print_info(self, text):
        self.infos.append(text)

    def print_blank(self):
        pass

    def flush(self):
        pass


def _make_molecule(xyz_string, charge=0):

    molecule = Molecule.read_smiles(xyz_string)
    molecule.set_charge(charge)
    return molecule


def _make_fluoride_solute(charge=-1):

    return _make_molecule("[F]", charge=charge)


def _make_methane_solute(charge=0):

    return _make_molecule("C", charge=charge)


def _make_water():

    return _make_molecule("O")


def _make_carbon_dioxide():

    return _make_molecule("O=C=O")


def _make_ibuprofenate_solute():

    return _make_molecule("CC(C)CC1=CC=C(C=C1)C(C)C(=O)[O-]")


def _make_ammonium_solute(charge=1):

    return _make_molecule("[NH4+]", charge=charge)


def _make_propylene_glycol():

    return _make_molecule("CC(CO)O")


def _target_density_for_count(solute, solvent, desired_count, box):

    builder = SolvationBuilder(ostream=RecordingOutput())
    solute_volume_nm3 = builder._get_volume(solute) * 1.0e-3
    volume_nm3 = box[0] * box[1] * box[2] * 1.0e-3 - solute_volume_nm3
    mass_per_molecule = sum(solvent.get_masses()) * 1.0e-3 / 6.022e23
    density = desired_count * mass_per_molecule / (volume_nm3 * 1.0e-27)

    # Add a tiny margin so int(mols_per_nm3 * volume_nm3) stays at desired_count.
    return density * 1.0001


def _system_density(solute, solvents, counts, box):

    total_mass = sum(solute.get_masses())
    for solvent, count in zip(solvents, counts):
        total_mass += count * sum(solvent.get_masses())

    mass_kg = total_mass * 1.0e-3 / 6.02214076e23
    volume_m3 = box[0] * box[1] * box[2] * 1.0e-30

    return mass_kg / volume_m3


def _make_forcefield(molecule, water_model=None):

    ff_gen = MMForceFieldGenerator()
    ff_gen.ostream.mute()

    if water_model is None:
        ff_gen.create_topology(molecule, resp=False)
    else:
        ff_gen.create_topology(molecule, water_model=water_model)

    return ff_gen


@pytest.mark.skipif(MPI.COMM_WORLD.Get_size() != 1,
                    reason='SolvationBuilder tests require a single MPI rank')
class TestSolvationBuilder:

    def test_solvate_replaces_solvent_slots_with_counterions(self):

        np.random.seed(0)
        solute = _make_fluoride_solute()
        solvent = _make_water()
        box = [25.0, 25.0, 25.0]
        target_density = _target_density_for_count(solute, solvent, 3, box)
        builder = SolvationBuilder(ostream=RecordingOutput())

        builder.solvate(solute,
                        solvent='other',
                        solvent_molecule=solvent,
                        target_density=target_density,
                        neutralize=True,
                        box=box)

        system_labels = [atom[1] for atom in builder.system]

        assert builder.added_solvent_counts == [2]
        assert builder.added_counterions == 1
        assert system_labels.count('O') == 2
        assert system_labels.count('Na') == 1

    def test_solvate_keeps_full_solvent_count_without_neutralization(self):

        np.random.seed(0)
        solute = _make_methane_solute()
        solvent = _make_water()
        box = [10.0, 10.0, 10.0]
        target_density = _target_density_for_count(solute, solvent, 3, box)
        builder = SolvationBuilder(ostream=RecordingOutput())

        builder.solvate(solute,
                        solvent='other',
                        solvent_molecule=solvent,
                        target_density=target_density,
                        neutralize=False,
                        box=box)

        system_labels = [atom[1] for atom in builder.system]

        assert builder.added_solvent_counts == [3]
        assert builder.added_counterions == 0
        assert system_labels.count('O') == 3
        assert 'Na' not in system_labels

    def test_solvate_disables_random_rotation_for_linear_solvent(self):

        np.random.seed(1)
        solute = _make_methane_solute()
        linear_solvent = _make_carbon_dioxide()
        box = [10.0, 10.0, 10.0]
        target_density = _target_density_for_count(solute, linear_solvent, 1,
                                                   box)
        builder = SolvationBuilder(ostream=RecordingOutput())

        builder.solvate(solute,
                        solvent='other',
                        solvent_molecule=linear_solvent,
                        target_density=target_density,
                        neutralize=False,
                        box=box)

        assert builder.random_rotation is False

    def test_solvate_restores_random_rotation_for_non_linear_solvent(self):

        np.random.seed(1)
        solute = _make_methane_solute()
        linear_solvent = _make_carbon_dioxide()
        bent_solvent = _make_water()
        box = [10.0, 10.0, 10.0]
        linear_density = _target_density_for_count(solute, linear_solvent, 1,
                                                   box)
        bent_density = _target_density_for_count(solute, bent_solvent, 1, box)
        builder = SolvationBuilder(ostream=RecordingOutput())

        builder.solvate(solute,
                        solvent='other',
                        solvent_molecule=linear_solvent,
                        target_density=linear_density,
                        neutralize=False,
                        box=box)
        builder.solvate(solute,
                        solvent='other',
                        solvent_molecule=bent_solvent,
                        target_density=bent_density,
                        neutralize=False,
                        box=box)

        assert builder.random_rotation is True

    def test_solvent_properties_returns_known_ethanol_values(self):

        builder = SolvationBuilder(ostream=RecordingOutput())

        mols_per_nm3, density, smiles_code = builder._solvent_properties(
            'ethanol')

        assert mols_per_nm3 == 10.3
        assert density == 789
        assert smiles_code == 'CCO'

    def test_solvate_rejects_unknown_solvent_name(self):

        solute = _make_methane_solute()
        builder = SolvationBuilder(ostream=RecordingOutput())

        with pytest.raises(ValueError, match='Unsupported solvent'):
            builder.solvate(solute,
                            solvent='not-a-solvent',
                            neutralize=False,
                            box=[10.0, 10.0, 10.0])

    def test_solvate_reports_partial_fill_after_repeated_insertion_failures(
            self):

        solute = _make_methane_solute()
        solvent = _make_water()
        box = [10.0, 10.0, 10.0]
        target_density = _target_density_for_count(solute, solvent, 5, box)
        builder = SolvationBuilder(ostream=RecordingOutput())
        builder.failures_factor = 0.4

        call_counter = {'count': 0}

        def fake_insert(molecule, tree):
            call_counter['count'] += 1
            if call_counter['count'] == 1:
                shift = np.array([0.5, 0.5, 0.5])
                coords = molecule.get_coordinates_in_angstrom() + shift
                return [(len(builder.system), label, coord)
                        for label, coord in zip(molecule.get_labels(), coords)]
            return None

        builder._insert_molecule = fake_insert

        builder.solvate(solute,
                        solvent='other',
                        solvent_molecule=solvent,
                        target_density=target_density,
                        neutralize=False,
                        box=box)

        assert builder.added_solvent_counts == [1]
        assert ('Solvated system with 1 solvent molecules out of 5 requested'
                in builder.ostream.infos)

    def test_solvate_tracks_only_successfully_inserted_counterions(self):

        np.random.seed(0)
        solute = _make_fluoride_solute()
        solvent = _make_water()
        box = [10.0, 10.0, 10.0]
        target_density = _target_density_for_count(solute, solvent, 1, box)
        builder = SolvationBuilder(ostream=RecordingOutput())

        call_counter = {'count': 0}

        def fake_insert(molecule, tree):
            call_counter['count'] += 1
            if molecule.get_labels() == ['Na']:
                return None

            shift = np.array([0.5, 0.5, 0.5])
            coords = molecule.get_coordinates_in_angstrom() + shift
            return [(len(builder.system), label, coord)
                    for label, coord in zip(molecule.get_labels(), coords)]

        builder._insert_molecule = fake_insert

        builder.solvate(solute,
                        solvent='other',
                        solvent_molecule=solvent,
                        target_density=target_density,
                        neutralize=True,
                        box=box)

        system_labels = [atom[1] for atom in builder.system]

        assert builder.added_solvent_counts == [0]
        assert builder.added_counterions == 0
        assert 'Na' not in system_labels
        assert ('Solvated system with 0 solvent molecules out of 1 requested'
                in builder.ostream.infos)

    def test_solvate_resets_stale_equilibration_flag(self):

        np.random.seed(0)
        solute = _make_methane_solute()
        solvent = _make_water()
        box = [25.0, 25.0, 25.0]
        target_density = _target_density_for_count(solute, solvent, 1, box)
        builder = SolvationBuilder(ostream=RecordingOutput())
        builder.equilibration_flag = True

        builder.solvate(solute,
                        solvent='other',
                        solvent_molecule=solvent,
                        target_density=target_density,
                        neutralize=False,
                        box=box)

        assert builder.equilibration_flag is False

    def test_solvate_reports_skipped_equilibration_when_files_are_missing(
            self, monkeypatch):

        np.random.seed(0)
        solute = _make_methane_solute()
        solvent = _make_water()
        box = [10.0, 10.0, 10.0]
        target_density = _target_density_for_count(solute, solvent, 1, box)
        builder = SolvationBuilder(ostream=RecordingOutput())

        def fake_perform_equilibration(water_model=None, steps=None):
            raise ValueError('missing forcefield include')

        monkeypatch.setattr(builder, 'perform_equilibration',
                            fake_perform_equilibration)

        builder.solvate(solute,
                        solvent='other',
                        solvent_molecule=solvent,
                        target_density=target_density,
                        neutralize=False,
                        equilibrate=True,
                        equilibration_steps=1,
                        box=box)

        assert builder.equilibration_flag is False
        assert 'Equilibration skipped due to missing files' in builder.ostream.infos
        assert 'Equilibrating the system' not in builder.ostream.infos

    def test_unlink_if_exists_handles_present_and_missing_files(self, tmp_path):

        builder = SolvationBuilder(ostream=RecordingOutput())
        existing_file = tmp_path / 'equilibrated_system.pdb'
        missing_file = tmp_path / 'missing_equilibrated_system.pdb'
        existing_file.write_text('temporary test file')

        builder._unlink_if_exists(existing_file)
        builder._unlink_if_exists(missing_file)

        assert not existing_file.exists()
        assert not missing_file.exists()

    def test_counterion_molecules_rejects_unknown_ion_name(self):

        builder = SolvationBuilder(ostream=RecordingOutput())
        builder.ion_name = 'Ca'

        with pytest.raises(VeloxChemError, match='Unsupported counterion'):
            builder._counterion_molecules()

    def test_solvate_reports_neutral_solute_when_neutralization_is_requested(
            self):

        np.random.seed(0)
        solute = _make_methane_solute()
        solvent = _make_water()
        box = [10.0, 10.0, 10.0]
        target_density = _target_density_for_count(solute, solvent, 1, box)
        builder = SolvationBuilder(ostream=RecordingOutput())

        builder.solvate(solute,
                        solvent='other',
                        solvent_molecule=solvent,
                        target_density=target_density,
                        neutralize=True,
                        box=box)

        assert builder.added_counterions == 0
        assert builder.counterion is None
        assert 'The solute is neutral, no counterions will be added' in builder.ostream.infos

    def test_solvate_positive_solute_uses_negative_counterions(self):

        np.random.seed(0)
        solute = _make_ammonium_solute()
        solvent = _make_water()
        box = [10.0, 10.0, 10.0]
        target_density = _target_density_for_count(solute, solvent, 3, box)
        builder = SolvationBuilder(ostream=RecordingOutput())

        builder.solvate(solute,
                        solvent='other',
                        solvent_molecule=solvent,
                        target_density=target_density,
                        neutralize=True,
                        box=box)

        system_labels = [atom[1] for atom in builder.system]

        assert builder.added_solvent_counts == [2]
        assert builder.added_counterions == 1
        assert builder.ion_name == builder.ncharge
        assert system_labels.count('O') == 2
        assert system_labels.count('Cl') == 1
        assert 'The solute has a charge of 1, adding 1Cl counterions' in builder.ostream.infos

    def test_solvate_rejects_box_smaller_than_solute_volume(self):

        solute = _make_methane_solute()
        solvent = _make_water()
        builder = SolvationBuilder(ostream=RecordingOutput())

        with pytest.raises(VeloxChemError,
                           match='available solvent volume must be positive'):
            builder.solvate(solute,
                            solvent='other',
                            solvent_molecule=solvent,
                            target_density=1000,
                            neutralize=False,
                            box=[1.0, 1.0, 1.0])

    @pytest.mark.parametrize('target_density', [0, -100])
    def test_solvate_rejects_non_positive_target_density(self, target_density):

        solute = _make_methane_solute()
        solvent = _make_water()
        builder = SolvationBuilder(ostream=RecordingOutput())

        with pytest.raises(VeloxChemError,
                           match='target density must be positive'):
            builder.solvate(solute,
                            solvent='other',
                            solvent_molecule=solvent,
                            target_density=target_density,
                            neutralize=False,
                            box=[10.0, 10.0, 10.0])

    @pytest.mark.parametrize(
        ('solvent', 'solvent_molecule', 'target_density', 'message'),
        [
            ('other', None, 1000,
             "The solvent molecule must be provided if the solvent is 'other'"),
            ('other', _make_water(), None,
             "The target density must be provided if the solvent is 'other'"),
            ('itself', None, None,
             "The target density must be provided if the solvent is 'itself'"),
        ])
    def test_solvate_validates_required_typical_solvent_inputs(
            self, solvent, solvent_molecule, target_density, message):

        solute = _make_methane_solute()
        builder = SolvationBuilder(ostream=RecordingOutput())

        with pytest.raises(ValueError, match=message):
            builder.solvate(solute,
                            solvent=solvent,
                            solvent_molecule=solvent_molecule,
                            target_density=target_density,
                            neutralize=False,
                            box=[10.0, 10.0, 10.0])

    @pytest.mark.parametrize('proportion', ([0, 0], [1, 0], [1, -1]))
    def test_custom_solvate_rejects_non_positive_proportions(self, proportion):

        solute = _make_methane_solute()
        solvents = [_make_water(), _make_molecule('CO')]
        builder = SolvationBuilder(ostream=RecordingOutput())

        with pytest.raises(VeloxChemError,
                           match='proportions must be positive'):
            builder.custom_solvate(solute,
                                   solvents=solvents,
                                   proportion=proportion,
                                   box_size=[10.0, 10.0, 10.0])

    def test_custom_solvate_rejects_box_smaller_than_solute_volume(self):

        solute = _make_methane_solute()
        solvents = [_make_water(), _make_molecule('CO')]
        builder = SolvationBuilder(ostream=RecordingOutput())

        with pytest.raises(VeloxChemError,
                           match='available solvent volume must be positive'):
            builder.custom_solvate(solute,
                                   solvents=solvents,
                                   proportion=[1, 1],
                                   box_size=[1.0, 1.0, 1.0])

    def test_custom_solvate_rejects_box_without_solvent_capacity(self):

        solute = _make_methane_solute()
        solvents = [_make_water(), _make_molecule('CO')]
        builder = SolvationBuilder(ostream=RecordingOutput())

        with pytest.raises(VeloxChemError,
                           match='too small to fit any solvent molecules'):
            builder.custom_solvate(solute,
                                   solvents=solvents,
                                   proportion=[1, 1],
                                   box_size=[4.0, 4.0, 4.0])

    def test_custom_solvate_rejects_solvents_without_nonhydrogen_atoms(self):

        solute = _make_methane_solute()
        hydrogen = _make_molecule('[H][H]')
        builder = SolvationBuilder(ostream=RecordingOutput())

        with pytest.raises(VeloxChemError,
                           match='at least one non-hydrogen atom'):
            builder.custom_solvate(solute,
                                   solvents=[hydrogen],
                                   proportion=[1],
                                   box_size=[10.0, 10.0, 10.0])

    def test_custom_solvate_preserves_small_mixture_without_zeroing_components(
            self):

        solute = _make_methane_solute()
        solvents = [_make_water(), _make_molecule('CO')]
        builder = SolvationBuilder(ostream=RecordingOutput())

        insert_counter = {'count': 0}

        def fake_insert(molecule, tree):
            shift = np.array([0.5 + insert_counter['count'], 0.5, 0.5])
            insert_counter['count'] += 1
            coords = molecule.get_coordinates_in_angstrom() + shift
            return [(len(builder.system), label, coord)
                    for label, coord in zip(molecule.get_labels(), coords)]

        builder._insert_molecule = fake_insert

        builder.custom_solvate(solute,
                               solvents=solvents,
                               proportion=[1, 2],
                               box_size=[5.9, 5.9, 5.9])

        assert builder.quantities == [1, 1]
        assert builder.added_solvent_counts == [1, 1]

    def test_custom_solvate_stops_after_repeated_insertion_failures(self):

        solute = _make_methane_solute()
        builder = SolvationBuilder(ostream=RecordingOutput())
        builder.failures_factor = 0.4

        call_counter = {'count': 0}

        def fake_insert(molecule, tree):
            call_counter['count'] += 1
            if call_counter['count'] == 1:
                shift = np.array([0.5, 0.5, 0.5])
                coords = molecule.get_coordinates_in_angstrom() + shift
                return [(len(builder.system), label, coord)
                        for label, coord in zip(molecule.get_labels(), coords)]
            if call_counter['count'] > 3:
                raise RuntimeError(
                    'custom_solvate kept inserting after repeated failures')
            return None

        builder._insert_molecule = fake_insert

        builder.custom_solvate(solute,
                               solvents=[_make_water()],
                               proportion=[1],
                               box_size=[6.0, 6.0, 6.0])

        assert builder.added_solvent_counts == [1]
        assert ('Solvated system with 1 solvent molecules out of 5 requested'
                in builder.ostream.infos)

    def test_pack_solvent_molecules_uses_dynamic_batch_sizing(self):

        builder = SolvationBuilder(ostream=RecordingOutput())
        builder.solvents = [_make_water()]
        builder.quantities = [builder.acceleration_threshold]

        insert_counter = {'count': 0}
        seen_trees = []

        class FakeKDTree:

            def __init__(self, coordinates):
                self.coordinates = coordinates

        def fake_insert(molecule, tree):
            seen_trees.append(tree)
            atom_idx = insert_counter['count']
            insert_counter['count'] += 1
            return [(atom_idx, molecule.get_labels()[0],
                     np.array([float(atom_idx), 0.0, 0.0]))]

        builder._insert_molecule = fake_insert

        tree = builder._pack_solvent_molecules(FakeKDTree,
                                               None,
                                               use_dynamic_batch=True,
                                               max_batch_size=20,
                                               min_batch_size=10)

        assert builder.added_solvent_counts == [builder.acceleration_threshold]
        assert insert_counter['count'] == builder.acceleration_threshold
        assert seen_trees[:20] == [None] * 20
        assert isinstance(seen_trees[20], FakeKDTree)
        assert isinstance(tree, FakeKDTree)
        assert builder._compute_batch_size(added_count=500,
                                           max_batch_size=20,
                                           min_batch_size=10,
                                           total_quantity=1000) == 10

    def test_custom_solvate_keeps_mixture_components_close_for_ibuprofenate(
            self):

        np.random.seed(0)
        solute = _make_ibuprofenate_solute()
        solvents = [_make_water(), _make_propylene_glycol()]
        builder = SolvationBuilder(ostream=RecordingOutput())

        builder.custom_solvate(solute,
                               solvents=solvents,
                               proportion=[50, 50],
                               box_size=[25.0, 25.0, 25.0])

        first_count, second_count = builder.added_solvent_counts
        reference_density = _system_density(solute, solvents,
                                            [first_count, second_count],
                                            [25.0, 25.0, 25.0])
        actual_density = builder._compute_density(25.0 * 25.0 * 25.0 * 1.0e-3)

        assert first_count >= 50
        assert second_count >= 50
        assert abs(first_count - second_count) <= int(
            np.ceil(0.05 * max(first_count, second_count)))
        assert actual_density == pytest.approx(reference_density, rel=0.01)

    def test_write_gromacs_files_writes_system_topology_and_gro(self, tmp_path):

        np.random.seed(0)
        solute = _make_fluoride_solute()
        solvent = _make_water()
        box = [10.0, 10.0, 10.0]
        target_density = _target_density_for_count(solute, solvent, 3, box)
        builder = SolvationBuilder(ostream=RecordingOutput())

        builder.solvate(solute,
                        solvent='other',
                        solvent_molecule=solvent,
                        target_density=target_density,
                        neutralize=True,
                        box=box)

        solute_ff = _make_forcefield(solute)
        solvent_ff = _make_forcefield(solvent, water_model='cspce')

        builder.workdir = str(tmp_path)
        builder.write_gromacs_files(solute_ff=solute_ff,
                                    solvent_ffs=[solvent_ff])

        system_top = (tmp_path / 'system.top').read_text()
        solute_top = (tmp_path / 'solute.top').read_text()
        solute_itp = (tmp_path / 'solute.itp').read_text()
        system_gro = (tmp_path / 'system.gro').read_text()

        assert (tmp_path / 'solute.gro').exists()
        assert (tmp_path / 'solvent_1.itp').exists()
        assert '#include "solute.itp"' in system_top
        assert '#include "solvent_1.itp"' in system_top
        assert f'#include "{builder.parent_forcefield}.ff/forcefield.itp"' in system_top
        assert f'#include "{builder.parent_forcefield}.ff/ions.itp"' in system_top
        assert 'MOL 1' in system_top
        assert 'SOL1 2' in system_top
        assert 'NA 1' in system_top
        assert '[ atomtypes ]' in solute_top
        assert '[ atomtypes ]' in system_top
        assert '[ atomtypes ]' not in solute_itp
        assert 'SOL1' in system_gro
        assert 'NA' in system_gro

    def test_write_gromacs_files_itself_updates_liquid_topology(self, tmp_path):

        np.random.seed(0)
        solute = _make_methane_solute()
        box = [10.0, 10.0, 10.0]
        target_density = _target_density_for_count(solute, solute, 2, box)
        builder = SolvationBuilder(ostream=RecordingOutput())

        builder.solvate(solute,
                        solvent='itself',
                        target_density=target_density,
                        neutralize=False,
                        box=box)

        solute_ff = _make_forcefield(solute)

        builder.workdir = tmp_path
        builder.write_gromacs_files(solute_ff=solute_ff)

        liquid_top = (tmp_path / 'liquid.top').read_text()
        liquid_gro = (tmp_path / 'liquid.gro').read_text()

        assert f'MOL               {builder.added_solvent_counts[0] + 1}' in liquid_top
        assert liquid_gro.splitlines()[1] == str(
            builder.system_molecule.number_of_atoms())
        assert 'MOL' in liquid_gro

    def test_write_openmm_files_writes_system_pdb_and_cleans_stale_file(
            self, tmp_path):

        np.random.seed(1)
        solute = _make_methane_solute()
        solvent = _make_water()
        box = [10.0, 10.0, 10.0]
        target_density = _target_density_for_count(solute, solvent, 1, box)
        builder = SolvationBuilder(ostream=RecordingOutput())

        builder.solvate(solute,
                        solvent='other',
                        solvent_molecule=solvent,
                        target_density=target_density,
                        neutralize=False,
                        box=box)
        builder.equilibration_flag = True

        solute_ff = _make_forcefield(solute)
        solvent_ff = _make_forcefield(solvent, water_model='cspce')

        builder.workdir = tmp_path
        stale_pdb = tmp_path / 'equilibrated_system.pdb'
        stale_pdb.write_text('stale file')

        builder.write_openmm_files(solute_ff=solute_ff, solvent_ffs=[solvent_ff])

        system_pdb = (tmp_path / 'system.pdb').read_text()

        assert not stale_pdb.exists()
        assert (tmp_path / 'solute.xml').exists()
        assert (tmp_path / 'solvent_1.xml').exists()
        assert 'HEADER    Generated by VeloxChem' in system_pdb
        assert 'HETATM' in system_pdb
        assert 'CONECT' in system_pdb
        assert 'S01' in system_pdb

    def test_write_openmm_files_itself_honors_write_pdb_only(self, tmp_path):

        np.random.seed(1)
        solute = _make_methane_solute()
        box = [10.0, 10.0, 10.0]
        target_density = _target_density_for_count(solute, solute, 1, box)
        builder = SolvationBuilder(ostream=RecordingOutput())

        builder.solvate(solute,
                        solvent='itself',
                        target_density=target_density,
                        neutralize=False,
                        box=box)
        builder.write_pdb_only = True

        solute_ff = _make_forcefield(solute)

        builder.workdir = str(tmp_path)
        builder.write_openmm_files(solute_ff=solute_ff)

        liquid_pdb = (tmp_path / 'liquid.pdb').read_text()

        assert (tmp_path / 'liquid.pdb').exists()
        assert not (tmp_path / 'liquid.xml').exists()
        assert 'HEADER    Generated by VeloxChem' in liquid_pdb
        assert 'CONECT' in liquid_pdb

    def test_solvate_equilibrates_with_explicit_small_step_count(self,
                                                                 tmp_path):

        pytest.importorskip("openmm")

        np.random.seed(0)
        solute = _make_methane_solute()
        solvent = _make_water()
        box = [25.0, 25.0, 25.0]
        target_density = _target_density_for_count(solute, solvent, 1, box)
        builder = SolvationBuilder(ostream=RecordingOutput())
        builder.workdir = tmp_path

        builder.solvate(solute,
                        solvent='other',
                        solvent_molecule=solvent,
                        target_density=target_density,
                        neutralize=False,
                        equilibrate=True,
                        equilibration_steps=1,
                        box=box)

        assert builder.equilibration_flag is True
        assert (tmp_path / 'equilibrated_system.pdb').exists()
        assert 'Equilibrating the system' in builder.ostream.infos
        assert 'Duration: 0.001 ps' in builder.ostream.infos

    def test_perform_equilibration_supports_itself_for_non_water_system(
            self, tmp_path):

        pytest.importorskip("openmm")

        np.random.seed(0)
        solute = _make_methane_solute()
        box = [25.0, 25.0, 25.0]
        target_density = _target_density_for_count(solute, solute, 1, box)
        builder = SolvationBuilder(ostream=RecordingOutput())
        builder.workdir = tmp_path

        builder.solvate(solute,
                        solvent='itself',
                        target_density=target_density,
                        neutralize=False,
                        box=box)
        builder.perform_equilibration(steps=1)

        assert (tmp_path / 'equilibrated_system.pdb').exists()
        assert builder.system_molecule.number_of_atoms() > 0
        assert not (tmp_path / 'liquid.gro').exists()
        assert not (tmp_path / 'liquid.top').exists()

    def test_perform_equilibration_supports_named_water_solvent(self,
                                                                tmp_path):

        pytest.importorskip("openmm")

        np.random.seed(0)
        solute = _make_methane_solute()
        box = [20.5, 20.5, 20.5]
        builder = SolvationBuilder(ostream=RecordingOutput())
        builder.workdir = tmp_path

        builder.solvate(solute,
                        solvent='cspce',
                        neutralize=False,
                        box=box)
        builder.perform_equilibration(steps=1)

        assert (tmp_path / 'equilibrated_system.pdb').exists()
        assert builder.system_molecule.number_of_atoms() > 0
        assert not (tmp_path / 'system.gro').exists()
        assert not (tmp_path / 'system.top').exists()

    def test_perform_equilibration_requires_water_model_for_pure_water_itself(
            self):

        np.random.seed(0)
        water = _make_water()
        box = [10.0, 10.0, 10.0]
        target_density = _target_density_for_count(water, water, 2, box)
        builder = SolvationBuilder(ostream=RecordingOutput())

        builder.solvate(water,
                        solvent='itself',
                        target_density=target_density,
                        neutralize=False,
                        box=box)

        with pytest.raises(VeloxChemError, match='water_model must be provided'):
            builder.perform_equilibration()
