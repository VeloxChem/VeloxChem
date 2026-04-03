import numpy as np
import pytest

from veloxchem.molecule import Molecule
from veloxchem.solvationbuilder import SolvationBuilder


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


def test_solvate_replaces_solvent_slots_with_counterions():

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


def test_solvate_keeps_full_solvent_count_without_neutralization():

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


def test_solvate_disables_random_rotation_for_linear_solvent():

    np.random.seed(1)
    solute = _make_methane_solute()
    linear_solvent = _make_carbon_dioxide()
    box = [10.0, 10.0, 10.0]
    target_density = _target_density_for_count(solute, linear_solvent, 1, box)
    builder = SolvationBuilder(ostream=RecordingOutput())

    builder.solvate(solute,
                    solvent='other',
                    solvent_molecule=linear_solvent,
                    target_density=target_density,
                    neutralize=False,
                    box=box)

    assert builder.random_rotation is False


def test_solvate_restores_random_rotation_for_non_linear_solvent():

    np.random.seed(1)
    solute = _make_methane_solute()
    linear_solvent = _make_carbon_dioxide()
    bent_solvent = _make_water()
    box = [10.0, 10.0, 10.0]
    linear_density = _target_density_for_count(solute, linear_solvent, 1, box)
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


def test_solvent_properties_returns_known_ethanol_values():

    builder = SolvationBuilder(ostream=RecordingOutput())

    mols_per_nm3, density, smiles_code = builder._solvent_properties('ethanol')

    assert mols_per_nm3 == 10.3
    assert density == 789
    assert smiles_code == 'CCO'


def test_solvate_rejects_unknown_solvent_name():

    solute = _make_methane_solute()
    builder = SolvationBuilder(ostream=RecordingOutput())

    with pytest.raises(ValueError, match='Unsupported solvent'):
        builder.solvate(solute,
                        solvent='not-a-solvent',
                        neutralize=False,
                        box=[10.0, 10.0, 10.0])
