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


def test_write_read_trexio_pyscf_scf(tmp_path):

    if MPI.COMM_WORLD.Get_rank() != mpi_master():
        return

    pyscf = pytest.importorskip('pyscf')
    from pyscf import gto, scf

    molecule, basis = _get_h2_and_basis()
    py_mol = gto.M(atom='H 0.0 0.0 0.0; H 0.0 0.0 0.74',
                   basis='sto-3g',
                   unit='Angstrom',
                   verbose=0)
    py_scf = scf.RHF(py_mol)
    py_energy = py_scf.kernel()

    scf_results = {
        'C_alpha': np.array(py_scf.mo_coeff, dtype=float),
        'E_alpha': np.array(py_scf.mo_energy, dtype=float),
        'occ_alpha': np.array(py_scf.mo_occ, dtype=float) / 2.0,
        'S': np.array(py_scf.get_ovlp(), dtype=float),
        'scf_energy': float(py_energy),
    }
    trexio_file = Path(tmp_path) / 'h2_pyscf.trexio'

    write_trexio(str(trexio_file), molecule, basis, scf_results=scf_results)

    actual_molecule, actual_basis = read_molecule_and_basis(str(trexio_file))
    actual_orbitals = read_molecular_orbitals(str(trexio_file))
    data = read_trexio(str(trexio_file))

    _assert_molecules_equal(actual_molecule, molecule)
    _assert_bases_equal(actual_basis, basis, molecule)
    assert actual_orbitals.get_orbitals_type() == molorb.rest
    np.testing.assert_allclose(actual_orbitals.alpha_to_numpy(),
                               scf_results['C_alpha'])
    np.testing.assert_allclose(actual_orbitals.ea_to_numpy(),
                               scf_results['E_alpha'])
    np.testing.assert_allclose(actual_orbitals.occa_to_numpy(),
                               scf_results['occ_alpha'])
    np.testing.assert_allclose(data['scf_energy'], py_energy)

    assert pyscf.__version__


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
