from mpi4py import MPI
import numpy as np
import pytest

from veloxchem.molecule import Molecule
from veloxchem.solvationbuilder import SolvationBuilder

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

    return _make_molecule("F", charge=charge)


def _make_methane_solute(charge=0):

    return _make_molecule("C", charge=charge)


def _make_water():

    return _make_molecule("O")


def _make_carbon_dioxide():

    return _make_molecule("O=C=O")


def _target_density_for_count(solute, solvent, desired_count, box):

    builder = SolvationBuilder(ostream=RecordingOutput())
    solute_volume_nm3 = builder._get_volume(solute) * 1.0e-3
    volume_nm3 = box[0] * box[1] * box[2] * 1.0e-3 - solute_volume_nm3
    mass_per_molecule = sum(solvent.get_masses()) * 1.0e-3 / 6.022e23
    density = desired_count * mass_per_molecule / (volume_nm3 * 1.0e-27)

    # Add a tiny margin so int(mols_per_nm3 * volume_nm3) stays at desired_count.
    return density * 1.0001


@pytest.mark.skipif(MPI.COMM_WORLD.Get_size() != 1,
                    reason='SolvationBuilder tests require a single MPI rank')
class TestSolvationBuilder:

    def test_solvate_replaces_solvent_slots_with_counterions(self):

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
        assert 'Failed to pack 4 out of 5 molecules after 2 attempts' in builder.ostream.infos

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
        assert 'Failed to pack 1 out of 1 molecules after 1 attempts' in builder.ostream.infos

    def test_solvate_resets_stale_equilibration_flag(self):

        np.random.seed(0)
        solute = _make_methane_solute()
        solvent = _make_water()
        box = [10.0, 10.0, 10.0]
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

        with pytest.raises(AssertionError, match='Unsupported counterion'):
            builder._counterion_molecules()

    def test_solvate_rejects_box_smaller_than_solute_volume(self):

        solute = _make_methane_solute()
        solvent = _make_water()
        builder = SolvationBuilder(ostream=RecordingOutput())

        with pytest.raises(AssertionError,
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

        with pytest.raises(AssertionError,
                           match='target density must be positive'):
            builder.solvate(solute,
                            solvent='other',
                            solvent_molecule=solvent,
                            target_density=target_density,
                            neutralize=False,
                            box=[10.0, 10.0, 10.0])

    @pytest.mark.parametrize('proportion', ([0, 0], [1, 0], [1, -1]))
    def test_custom_solvate_rejects_non_positive_proportions(self, proportion):

        solute = _make_methane_solute()
        solvents = [_make_water(), _make_molecule('CO')]
        builder = SolvationBuilder(ostream=RecordingOutput())

        with pytest.raises(AssertionError,
                           match='proportions must be positive'):
            builder.custom_solvate(solute,
                                   solvents=solvents,
                                   proportion=proportion,
                                   box_size=[10.0, 10.0, 10.0])

    def test_custom_solvate_rejects_box_smaller_than_solute_volume(self):

        solute = _make_methane_solute()
        solvents = [_make_water(), _make_molecule('CO')]
        builder = SolvationBuilder(ostream=RecordingOutput())

        with pytest.raises(AssertionError,
                           match='available solvent volume must be positive'):
            builder.custom_solvate(solute,
                                   solvents=solvents,
                                   proportion=[1, 1],
                                   box_size=[1.0, 1.0, 1.0])

    def test_custom_solvate_rejects_box_without_solvent_capacity(self):

        solute = _make_methane_solute()
        solvents = [_make_water(), _make_molecule('CO')]
        builder = SolvationBuilder(ostream=RecordingOutput())

        with pytest.raises(AssertionError,
                           match='too small to fit any solvent molecules'):
            builder.custom_solvate(solute,
                                   solvents=solvents,
                                   proportion=[1, 1],
                                   box_size=[4.0, 4.0, 4.0])

    def test_custom_solvate_rejects_solvents_without_nonhydrogen_atoms(self):

        solute = _make_methane_solute()
        hydrogen = _make_molecule('[H][H]')
        builder = SolvationBuilder(ostream=RecordingOutput())

        with pytest.raises(AssertionError,
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
        assert 'Failed to pack 4 out of 5 molecules after 2 attempts' in builder.ostream.infos
