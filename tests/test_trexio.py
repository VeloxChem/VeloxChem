from pathlib import Path
import importlib
import importlib.util
import sys

from mpi4py import MPI
import numpy as np
import pytest

pytest.importorskip('trexio')
import trexio as pytrexio

from veloxchem import (Molecule, MolecularBasis, MolecularOrbitals,
                       ScfRestrictedDriver, molorb, mpi_master)
from veloxchem.errorhandler import VeloxChemError


def _load_veloxchem_trexio():

    source_path = Path(__file__).resolve().parents[1]
    module_path = source_path / 'src' / 'pymodule' / 'trexio.py'

    if module_path.is_file():
        spec = importlib.util.spec_from_file_location('veloxchem.trexio',
                                                      module_path)
        module = importlib.util.module_from_spec(spec)
        sys.modules['veloxchem.trexio'] = module
        spec.loader.exec_module(module)
        return module

    return importlib.import_module('veloxchem.trexio')


vlx_trexio = _load_veloxchem_trexio()
read_molecular_orbitals = vlx_trexio.read_molecular_orbitals
read_molecule = vlx_trexio.read_molecule
read_molecule_and_basis = vlx_trexio.read_molecule_and_basis
read_trexio = vlx_trexio.read_trexio
write_trexio = vlx_trexio.write_trexio


def _get_h2_and_basis():

    xyz_string = """2
    h2
    H 0.0 0.0 0.0
    H 0.0 0.0 0.74
    """

    molecule = Molecule.read_xyz_string(xyz_string)
    basis = MolecularBasis.read(molecule, 'sto-3g', ostream=None)

    return molecule, basis


def _get_water_and_basis():

    xyz_string = """3
    water
    O 0.0 0.0 0.0
    H 0.0 0.0 0.9584
    H 0.0 0.9277 -0.2396
    """

    molecule = Molecule.read_xyz_string(xyz_string)
    basis = MolecularBasis.read(molecule, 'sto-3g', ostream=None)

    return molecule, basis


def _get_auh2_and_ecp_basis():

    coord_string = """
    Au   0.000   0.000   0.000
     H   0.200   1.800  -1.100
     H   0.100  -0.900  -1.000
    """

    molecule = Molecule.read_str(coord_string, 'au')
    basis = MolecularBasis.read(molecule, 'def2-svp', ostream=None)

    return molecule, basis


def _basis_shell_data(molecule, basis):

    basis_sets = basis.basis_sets()
    basis_indices = basis.basis_sets_indices()
    shell_data = []

    for atom_index in range(molecule.number_of_atoms()):
        atom_basis = basis_sets[basis_indices[atom_index]]
        for basis_function in atom_basis.get_basis_functions():
            shell_data.append((
                atom_index,
                basis_function.get_angular_momentum(),
                np.array(basis_function.get_exponents()),
                np.array(basis_function.get_normalization_factors()),
            ))

    return shell_data


def _assert_molecules_equal(actual, expected):

    assert actual.number_of_atoms() == expected.number_of_atoms()
    assert actual.get_labels() == expected.get_labels()
    np.testing.assert_allclose(actual.get_coordinates_in_bohr(),
                               expected.get_coordinates_in_bohr())
    assert actual.get_charge() == expected.get_charge()
    assert actual.get_multiplicity() == expected.get_multiplicity()
    assert actual.number_of_alpha_electrons(
    ) == expected.number_of_alpha_electrons()
    assert actual.number_of_beta_electrons(
    ) == expected.number_of_beta_electrons()


def _assert_bases_equal(actual_basis, expected_basis, molecule):

    assert actual_basis.get_dimensions_of_basis(
    ) == expected_basis.get_dimensions_of_basis()
    assert (actual_basis.get_dimensions_of_primitive_basis() ==
            expected_basis.get_dimensions_of_primitive_basis())

    actual_shells = _basis_shell_data(molecule, actual_basis)
    expected_shells = _basis_shell_data(molecule, expected_basis)

    assert len(actual_shells) == len(expected_shells)
    for actual, expected in zip(actual_shells, expected_shells):
        assert actual[0] == expected[0]
        assert actual[1] == expected[1]
        np.testing.assert_allclose(actual[2], expected[2])
        np.testing.assert_allclose(actual[3], expected[3])


def _core_potential_signature(potential):

    return (
        np.array(potential.get_exponents(), dtype=float),
        np.array(potential.get_factors(), dtype=float),
        np.array(potential.get_radial_orders(), dtype=int),
    )


def _assert_ecp_equal(actual_basis, expected_basis):

    actual_sets = actual_basis.basis_sets()
    expected_sets = expected_basis.basis_sets()
    actual_indices = actual_basis.basis_sets_indices()
    expected_indices = expected_basis.basis_sets_indices()

    assert len(actual_indices) == len(expected_indices)

    for actual_basis_index, expected_basis_index in zip(actual_indices,
                                                        expected_indices):
        actual_atom_basis = actual_sets[actual_basis_index]
        expected_atom_basis = expected_sets[expected_basis_index]

        assert actual_atom_basis.has_ecp() == expected_atom_basis.has_ecp()
        if not expected_atom_basis.has_ecp():
            continue

        actual_ecp = actual_atom_basis.get_ecp_potential()
        expected_ecp = expected_atom_basis.get_ecp_potential()

        assert (actual_ecp.number_of_core_electrons() ==
                expected_ecp.number_of_core_electrons())
        assert list(actual_ecp.get_angular_momentums()) == list(
            expected_ecp.get_angular_momentums())

        actual_local = _core_potential_signature(actual_ecp.get_local_potential())
        expected_local = _core_potential_signature(
            expected_ecp.get_local_potential())
        for actual_values, expected_values in zip(actual_local, expected_local):
            np.testing.assert_allclose(actual_values, expected_values)

        actual_projected = actual_ecp.get_projected_potentials()
        expected_projected = expected_ecp.get_projected_potentials()
        assert len(actual_projected) == len(expected_projected)
        for actual_potential, expected_potential in zip(actual_projected,
                                                        expected_projected):
            actual_signature = _core_potential_signature(actual_potential)
            expected_signature = _core_potential_signature(expected_potential)
            for actual_values, expected_values in zip(actual_signature,
                                                      expected_signature):
                np.testing.assert_allclose(actual_values, expected_values)


def test_write_read_trexio_molecule_and_basis(tmp_path):

    if MPI.COMM_WORLD.Get_rank() != mpi_master():
        return

    molecule, basis = _get_h2_and_basis()
    trexio_file = Path(tmp_path) / 'h2.trexio'

    write_trexio(str(trexio_file), molecule, basis)
    actual_molecule, actual_basis = read_molecule_and_basis(str(trexio_file))

    _assert_molecules_equal(actual_molecule, molecule)
    _assert_bases_equal(actual_basis, basis, molecule)


def test_write_read_trexio_text_backend(tmp_path):

    if MPI.COMM_WORLD.Get_rank() != mpi_master():
        return

    molecule, basis = _get_h2_and_basis()
    trexio_file = Path(tmp_path) / 'h2_text.trexio'

    write_trexio(str(trexio_file), molecule, basis, backend='text')
    actual_molecule, actual_basis = read_molecule_and_basis(str(trexio_file),
                                                            backend='text')

    _assert_molecules_equal(actual_molecule, molecule)
    _assert_bases_equal(actual_basis, basis, molecule)


def test_write_trexio_portable_ao_basis_metadata(tmp_path):

    if MPI.COMM_WORLD.Get_rank() != mpi_master():
        return

    molecule, basis = _get_water_and_basis()
    trexio_file = Path(tmp_path) / 'water.trexio'

    write_trexio(str(trexio_file), molecule, basis)

    tf = pytrexio.File(str(trexio_file), 'r', pytrexio.TREXIO_HDF5)
    try:
        ao_shell = np.array(pytrexio.read_ao_shell(tf), dtype=int)
        basis_coefficient = np.array(pytrexio.read_basis_coefficient(tf),
                                     dtype=float)
        basis_prim_factor = np.array(pytrexio.read_basis_prim_factor(tf),
                                     dtype=float)
    finally:
        tf.close()

    np.testing.assert_array_equal(ao_shell, np.array([0, 1, 2, 2, 2, 3, 4]))
    assert not np.allclose(basis_prim_factor, 1.0)

    expected_coefficients = []
    expected_prim_factors = []
    for _, angular_momentum, exponents, norms in _basis_shell_data(
            molecule, basis):
        primitive_factors = vlx_trexio._primitive_normalization_factors(
            exponents, angular_momentum)
        expected_coefficients.extend((norms / primitive_factors).tolist())
        expected_prim_factors.extend(primitive_factors.tolist())

    np.testing.assert_allclose(basis_coefficient, expected_coefficients)
    np.testing.assert_allclose(basis_prim_factor, expected_prim_factors)


def test_write_read_trexio_basis_with_ecp(tmp_path):

    if MPI.COMM_WORLD.Get_rank() != mpi_master():
        return

    molecule, basis = _get_auh2_and_ecp_basis()
    trexio_file = Path(tmp_path) / 'auh2_ecp.trexio'

    write_trexio(str(trexio_file), molecule, basis)
    actual_molecule, actual_basis = read_molecule_and_basis(str(trexio_file))

    assert actual_molecule.number_of_atoms() == molecule.number_of_atoms()
    assert actual_molecule.get_labels() == molecule.get_labels()
    np.testing.assert_allclose(actual_molecule.get_coordinates_in_bohr(),
                               molecule.get_coordinates_in_bohr())
    assert actual_molecule.get_multiplicity() == molecule.get_multiplicity()
    _assert_bases_equal(actual_basis, basis, molecule)
    assert actual_basis.has_ecp()
    assert actual_basis.get_number_of_ecp_core_electrons(
    ) == basis.get_number_of_ecp_core_electrons()
    _assert_ecp_equal(actual_basis, basis)


def test_write_trexio_ecp_metadata(tmp_path):

    if MPI.COMM_WORLD.Get_rank() != mpi_master():
        return

    molecule, basis = _get_auh2_and_ecp_basis()
    trexio_file = Path(tmp_path) / 'auh2_ecp_metadata.trexio'

    write_trexio(str(trexio_file), molecule, basis)

    tf = pytrexio.File(str(trexio_file), 'r', pytrexio.TREXIO_HDF5)
    try:
        actual_z_core = np.array(pytrexio.read_ecp_z_core(tf), dtype=int)
        actual_max_ang = np.array(pytrexio.read_ecp_max_ang_mom_plus_1(tf),
                                  dtype=int)
        actual_nucleus_index = np.array(pytrexio.read_ecp_nucleus_index(tf),
                                        dtype=int)
        actual_ang_mom = np.array(pytrexio.read_ecp_ang_mom(tf), dtype=int)
        actual_exponent = np.array(pytrexio.read_ecp_exponent(tf), dtype=float)
        actual_coefficient = np.array(pytrexio.read_ecp_coefficient(tf),
                                      dtype=float)
        actual_power = np.array(pytrexio.read_ecp_power(tf), dtype=int)
    finally:
        tf.close()

    expected_z_core = []
    expected_max_ang = []
    expected_nucleus_index = []
    expected_ang_mom = []
    expected_exponent = []
    expected_coefficient = []
    expected_power = []

    basis_sets = basis.basis_sets()
    basis_indices = basis.basis_sets_indices()
    for atom_index in range(molecule.number_of_atoms()):
        atom_basis = basis_sets[basis_indices[atom_index]]
        if not atom_basis.has_ecp():
            expected_z_core.append(0)
            expected_max_ang.append(0)
            continue

        atom_ecp = atom_basis.get_ecp_potential()
        projected_angular = list(atom_ecp.get_angular_momentums())
        local_angular = max(projected_angular) + 1 if projected_angular else 0

        expected_z_core.append(int(atom_ecp.number_of_core_electrons()))
        expected_max_ang.append(local_angular)

        local_potential = atom_ecp.get_local_potential()
        for exponent, coefficient, power in zip(local_potential.get_exponents(),
                                                local_potential.get_factors(),
                                                local_potential.get_radial_orders()):
            expected_nucleus_index.append(atom_index)
            expected_ang_mom.append(local_angular)
            expected_exponent.append(exponent)
            expected_coefficient.append(coefficient)
            expected_power.append(power)

        for proj_ang, projected_potential in zip(projected_angular,
                                                 atom_ecp.get_projected_potentials()):
            for exponent, coefficient, power in zip(
                    projected_potential.get_exponents(),
                    projected_potential.get_factors(),
                    projected_potential.get_radial_orders()):
                expected_nucleus_index.append(atom_index)
                expected_ang_mom.append(proj_ang)
                expected_exponent.append(exponent)
                expected_coefficient.append(coefficient)
                expected_power.append(power)

    np.testing.assert_array_equal(actual_z_core, np.array(expected_z_core))
    np.testing.assert_array_equal(actual_max_ang, np.array(expected_max_ang))
    np.testing.assert_array_equal(actual_nucleus_index,
                                  np.array(expected_nucleus_index))
    np.testing.assert_array_equal(actual_ang_mom, np.array(expected_ang_mom))
    np.testing.assert_allclose(actual_exponent, np.array(expected_exponent,
                                                         dtype=float))
    np.testing.assert_allclose(actual_coefficient,
                               np.array(expected_coefficient, dtype=float))
    np.testing.assert_array_equal(actual_power, np.array(expected_power))


def test_write_read_trexio_restricted_scf(tmp_path):

    if MPI.COMM_WORLD.Get_rank() != mpi_master():
        return

    molecule, basis = _get_h2_and_basis()
    scf_drv = ScfRestrictedDriver()
    scf_drv.ostream.mute()
    scf_results = scf_drv.compute(molecule, basis)
    trexio_file = Path(tmp_path) / 'h2_scf.trexio'

    write_trexio(str(trexio_file), molecule, basis, scf_results=scf_results)

    actual_molecule = read_molecule(str(trexio_file))
    actual_orbitals = read_molecular_orbitals(str(trexio_file))
    data = read_trexio(str(trexio_file))

    _assert_molecules_equal(actual_molecule, molecule)
    assert 'basis' in data
    assert 'molecular_orbitals' in data

    assert actual_orbitals.get_orbitals_type() == molorb.rest
    np.testing.assert_allclose(actual_orbitals.alpha_to_numpy(),
                               scf_results['C_alpha'])
    np.testing.assert_allclose(actual_orbitals.ea_to_numpy(),
                               scf_results['E_alpha'])
    np.testing.assert_allclose(actual_orbitals.occa_to_numpy(),
                               scf_results['occ_alpha'])
    np.testing.assert_allclose(actual_orbitals.occb_to_numpy(),
                               scf_results['occ_beta'])


def test_write_read_trexio_unrestricted_orbitals(tmp_path):

    if MPI.COMM_WORLD.Get_rank() != mpi_master():
        return

    molecule, basis = _get_h2_and_basis()
    norb = basis.get_dimensions_of_basis()
    c_alpha = np.eye(norb)
    c_beta = np.fliplr(np.eye(norb))
    e_alpha = np.arange(norb, dtype=float)
    e_beta = e_alpha + 0.5
    occ_alpha = np.array([1.0, 0.0])
    occ_beta = np.array([0.0, 1.0])
    mol_orbs = MolecularOrbitals([c_alpha, c_beta], [e_alpha, e_beta],
                                 [occ_alpha, occ_beta], molorb.unrest)
    trexio_file = Path(tmp_path) / 'h2_unrest.trexio'

    write_trexio(str(trexio_file), molecule, basis, mol_orbs=mol_orbs)
    actual_orbitals = read_molecular_orbitals(str(trexio_file))

    assert actual_orbitals.get_orbitals_type() == molorb.unrest
    np.testing.assert_allclose(actual_orbitals.alpha_to_numpy(), c_alpha)
    np.testing.assert_allclose(actual_orbitals.beta_to_numpy(), c_beta)
    np.testing.assert_allclose(actual_orbitals.ea_to_numpy(), e_alpha)
    np.testing.assert_allclose(actual_orbitals.eb_to_numpy(), e_beta)
    np.testing.assert_allclose(actual_orbitals.occa_to_numpy(), occ_alpha)
    np.testing.assert_allclose(actual_orbitals.occb_to_numpy(), occ_beta)


def test_write_read_trexio_restricted_open_shell_orbitals(tmp_path):

    if MPI.COMM_WORLD.Get_rank() != mpi_master():
        return

    molecule, basis = _get_h2_and_basis()
    norb = basis.get_dimensions_of_basis()
    coefficients = np.eye(norb)
    energies = np.arange(norb, dtype=float)
    occ_alpha = np.array([1.0, 0.0])
    occ_beta = np.array([0.0, 0.0])
    mol_orbs = MolecularOrbitals([coefficients], [energies],
                                 [occ_alpha, occ_beta], molorb.restopen)
    trexio_file = Path(tmp_path) / 'h2_restopen.trexio'

    write_trexio(str(trexio_file), molecule, basis, mol_orbs=mol_orbs)
    actual_orbitals = read_molecular_orbitals(str(trexio_file))

    assert actual_orbitals.get_orbitals_type() == molorb.restopen
    np.testing.assert_allclose(actual_orbitals.alpha_to_numpy(), coefficients)
    np.testing.assert_allclose(actual_orbitals.beta_to_numpy(), coefficients)
    np.testing.assert_allclose(actual_orbitals.ea_to_numpy(), energies)
    np.testing.assert_allclose(actual_orbitals.eb_to_numpy(), energies)
    np.testing.assert_allclose(actual_orbitals.occa_to_numpy(), occ_alpha)
    np.testing.assert_allclose(actual_orbitals.occb_to_numpy(), occ_beta)


def test_primitive_normalization_rejects_unsupported_high_l():

    exponents = np.array([0.5, 1.0], dtype=float)

    with pytest.raises(VeloxChemError,
                       match='does not support primitive normalization'):
        vlx_trexio._primitive_normalization_factors(exponents, 7)
