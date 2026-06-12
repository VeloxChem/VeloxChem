from pathlib import Path

from mpi4py import MPI
import h5py
import numpy as np
import pytest

from veloxchem import (Molecule, MolecularBasis, OptimizationDriver, MpiTask,
                       OutputStream, ScfGradientDriver, ScfRestrictedDriver,
                       ScfUnrestrictedDriver, mpi_master)
from veloxchem.cppsolver import ComplexResponseSolver
from veloxchem.c6driver import C6Driver
from veloxchem.errorhandler import VeloxChemError
from veloxchem.lrsolver import LinearResponseSolver
from veloxchem.lreigensolver import LinearResponseEigenSolver
from veloxchem.lrsolverunrest import LinearResponseUnrestrictedSolver
from veloxchem.resultsio import (read_results, write_results_to_hdf5,
                                 write_rsp_full_solution_to_hdf5,
                                 write_rsp_results_to_hdf5,
                                 write_scf_results_to_hdf5)
from veloxchem.tdacppsolver import ComplexResponseTdaSolver
from veloxchem.tdaeigensolver import TdaEigenSolver
from veloxchem.vibrationalanalysis import VibrationalAnalysis


def _get_water_and_basis():

    xyz_string = """3
    xyz
    O   -0.1858140  -1.1749469   0.7662596
    H   -0.1285513  -0.8984365   1.6808606
    H   -0.0582782  -0.3702550   0.2638279
    """

    molecule = Molecule.read_xyz_string(xyz_string)
    basis = MolecularBasis.read(molecule, 'sto-3g', ostream=None)

    return molecule, basis


def _run_water_rsp_roundtrip(tmp_path, solver_cls, save_solutions=True):

    molecule, basis = _get_water_and_basis()
    base = str(Path(tmp_path) / solver_cls.__name__.lower())

    scf_drv = ScfRestrictedDriver()
    scf_drv.ostream.mute()
    base = scf_drv.comm.bcast(base, root=mpi_master())
    scf_drv.filename = base
    scf_results = scf_drv.compute(molecule, basis)

    rsp_drv = solver_cls()
    rsp_drv.ostream.mute()
    rsp_drv.filename = base
    rsp_drv.nstates = 2
    rsp_drv.max_subspace_dim = 10
    rsp_drv.save_solutions = save_solutions
    if solver_cls is LinearResponseEigenSolver:
        rsp_drv.esa = True
    rsp_results = rsp_drv.compute(molecule, basis, scf_results)

    return rsp_drv, rsp_results, base + '.h5'


def _run_cpp_rsp_roundtrip(tmp_path, solver_cls):

    xyz_string = """6
    xyz
    O -3.42904  1.55532  0.01546
    C -1.99249  1.74379  0.02665
    H -1.74709  2.74160  0.44749
    H -1.59636  1.67836 -1.00861
    H -1.51398  0.95881  0.64937
    H -3.84726  2.33620 -0.34927
    """

    molecule = Molecule.read_xyz_string(xyz_string)
    basis = MolecularBasis.read(molecule, 'def2-svp', ostream=None)
    base = str(Path(tmp_path) / solver_cls.__name__.lower())

    scf_drv = ScfRestrictedDriver()
    scf_drv.ostream.mute()
    base = scf_drv.comm.bcast(base, root=mpi_master())
    scf_drv.filename = base
    scf_drv.xcfun = 'hf'
    scf_drv.acc_type = 'l2_c2diis'
    scf_results = scf_drv.compute(molecule, basis)

    rsp_drv = solver_cls()
    rsp_drv.ostream.mute()
    rsp_drv.filename = base
    rsp_drv.property = 'absorption'
    rsp_drv.frequencies = [0.39]
    rsp_drv.damping = 0.02
    rsp_drv.non_equilibrium_solv = False
    rsp_drv.ri_auxiliary_basis = 'def2-universal-jfit'
    rsp_drv.conv_thresh = 1.0e-5
    rsp_drv.max_iter = 120
    rsp_results = rsp_drv.compute(molecule, basis, scf_results)

    return rsp_drv, rsp_results, base + '.h5'


def _run_lr_rsp_roundtrip(tmp_path, solver_cls):

    xyz_string = """3
    xyz
    O   -0.1858140  -1.1749469   0.7662596
    H   -0.1285513  -0.8984365   1.6808606
    H   -0.0582782  -0.3702550   0.2638279
    """

    molecule = Molecule.read_xyz_string(xyz_string)
    if solver_cls is LinearResponseUnrestrictedSolver:
        molecule.set_charge(1)
        molecule.set_multiplicity(2)

    basis = MolecularBasis.read(molecule, 'sto-3g', ostream=None)
    base = str(Path(tmp_path) / solver_cls.__name__.lower())

    scf_cls = (ScfUnrestrictedDriver
               if solver_cls is LinearResponseUnrestrictedSolver else
               ScfRestrictedDriver)
    scf_drv = scf_cls()
    scf_drv.ostream.mute()
    base = scf_drv.comm.bcast(base, root=mpi_master())
    scf_drv.filename = base
    scf_results = scf_drv.compute(molecule, basis)

    rsp_drv = solver_cls()
    rsp_drv.ostream.mute()
    rsp_drv.filename = base
    rsp_drv.a_components = 'xz'
    rsp_drv.b_components = 'yz'
    rsp_drv.frequencies = [0.0, 0.05]
    rsp_drv.non_equilibrium_solv = False
    rsp_drv.conv_thresh = 1.0e-5
    rsp_drv.max_iter = 120
    rsp_results = rsp_drv.compute(molecule, basis, scf_results)

    return rsp_drv, rsp_results, base + '.h5'


def _run_c6_rsp_roundtrip(tmp_path):

    molecule, basis = _get_water_and_basis()
    base = str(Path(tmp_path) / 'c6driver')

    scf_drv = ScfRestrictedDriver()
    scf_drv.ostream.mute()
    base = scf_drv.comm.bcast(base, root=mpi_master())
    scf_drv.filename = base
    scf_results = scf_drv.compute(molecule, basis)

    rsp_drv = C6Driver()
    rsp_drv.ostream.mute()
    rsp_drv.filename = base
    rsp_drv.n_points = 3
    rsp_drv.conv_thresh = 1.0e-5
    rsp_drv.max_iter = 120
    rsp_results = rsp_drv.compute(molecule, basis, scf_results)

    return rsp_drv, rsp_results, base + '.h5'


def _assert_roundtrip_equal(expected, actual):

    def assert_value_equal(expected_value, actual_value):
        if expected_value is None:
            assert actual_value is None
        elif isinstance(expected_value, np.ndarray):
            np.testing.assert_allclose(actual_value, expected_value)
        elif isinstance(expected_value, dict):
            assert isinstance(actual_value, dict)
            assert actual_value.keys() == expected_value.keys()
            for key in expected_value:
                assert_value_equal(expected_value[key], actual_value[key])
        elif isinstance(expected_value, list):
            assert isinstance(actual_value, list)
            assert len(actual_value) == len(expected_value)
            for expected_item, actual_item in zip(expected_value, actual_value):
                assert_value_equal(expected_item, actual_item)
        elif isinstance(expected_value, tuple):
            assert isinstance(actual_value, tuple)
            assert len(actual_value) == len(expected_value)
            for expected_item, actual_item in zip(expected_value, actual_value):
                assert_value_equal(expected_item, actual_item)
        else:
            assert actual_value == expected_value

    assert expected.keys() == actual.keys()

    for key, expected_value in expected.items():
        assert_value_equal(expected_value, actual[key])


def _assert_opt_results_roundtrip_equal(expected, actual):

    expected_keys = set(expected.keys()) - {'final_molecule'}
    assert actual.keys() == expected_keys

    for key in expected_keys:
        expected_value = expected[key]
        actual_value = actual[key]

        if isinstance(expected_value, np.ndarray):
            np.testing.assert_allclose(actual_value, expected_value)
        elif isinstance(expected_value, list):
            assert isinstance(actual_value, list)
            assert len(actual_value) == len(expected_value)
            if expected_value and isinstance(expected_value[0], str):
                assert actual_value == expected_value
            else:
                np.testing.assert_allclose(np.array(actual_value),
                                           np.array(expected_value))
        else:
            assert actual_value == expected_value


def test_write_scf_results_to_hdf5_stores_all_entries_with_metadata(tmp_path):

    if MPI.COMM_WORLD.Get_rank() != mpi_master():
        return

    h5file = Path(tmp_path) / 'scf_results.h5'
    with h5py.File(h5file, 'w'):
        pass

    scf_results = {
        'scf_type': 'restricted',
        'scf_energy': -75.123,
        'restart': False,
        'grid_level': 4,
        'xcfun': 'b3lyp',
        'scf_history_energy': np.array([-75.0, -75.123]),
        'scf_history_gradient_norm': np.array([1.0e-3, 1.0e-6]),
        'S': np.eye(2),
        'dipole_moment': np.array([0.1, 0.2, 0.3]),
        'F': (np.eye(2), np.eye(2) * 2.0),
        'potfile': None,
    }

    write_scf_results_to_hdf5(str(h5file), scf_results)

    with h5py.File(h5file, 'r') as h5f:
        assert 'scf' in h5f

        scf_group = h5f['scf']
        assert scf_group.attrs['value_type'] == 'dict'

        assert scf_group['scf_type'].attrs['value_type'] == 'str'
        assert scf_group['scf_type'][()].decode('utf-8') == 'restricted'

        assert scf_group['scf_energy'].attrs['value_type'] == 'float'
        assert scf_group['restart'].attrs['value_type'] == 'bool'
        assert scf_group['grid_level'].attrs['value_type'] == 'int'

        assert scf_group['S'].attrs['value_type'] == 'ndarray'
        np.testing.assert_allclose(np.array(scf_group['S']), np.eye(2))

        assert scf_group['potfile'].attrs['value_type'] == 'none'
        assert scf_group['potfile'].shape == (0,)

        f_group = scf_group['F']
        assert f_group.attrs['value_type'] == 'tuple'
        assert f_group.attrs['length'] == 2
        np.testing.assert_allclose(np.array(f_group['0']), np.eye(2))
        np.testing.assert_allclose(np.array(f_group['1']), np.eye(2) * 2.0)

        scf_history_ene = scf_group['scf_history_energy']
        assert scf_history_ene.attrs['value_type'] == 'ndarray'
        np.testing.assert_allclose(scf_history_ene, np.array([-75.0, -75.123]))

        scf_history_grad = scf_group['scf_history_gradient_norm']
        assert scf_history_grad.attrs['value_type'] == 'ndarray'
        np.testing.assert_allclose(scf_history_grad, np.array([1.0e-3, 1.0e-6]))


def test_write_results_to_hdf5_stores_requested_group(tmp_path):

    if MPI.COMM_WORLD.Get_rank() != mpi_master():
        return

    h5file = Path(tmp_path) / 'generic_results.h5'
    with h5py.File(h5file, 'w') as h5f:
        h5f.create_dataset('basis_set', data=np.bytes_(['def2-svp']))

    vib_results = {
        'frequencies': np.array([100.0, 200.0]),
        'intensities': [0.1, 0.2],
        'projected': True,
    }

    write_results_to_hdf5(str(h5file),
                          'vib',
                          vib_results,
                          value_label='vibrational result')

    with h5py.File(h5file, 'r') as h5f:
        assert 'vib' in h5f
        assert h5f['vib'].attrs['value_type'] == 'dict'
        assert 'basis_set' in h5f

    recovered = read_results(str(h5file), 'vib')
    np.testing.assert_allclose(recovered['frequencies'],
                               vib_results['frequencies'])
    assert recovered['intensities'] == vib_results['intensities']
    assert recovered['projected'] is vib_results['projected']


def test_write_rsp_results_pairs_solution_keys_and_matrix(tmp_path):

    if MPI.COMM_WORLD.Get_rank() != mpi_master():
        return

    h5file = Path(tmp_path) / 'rsp_solutions.h5'
    with h5py.File(h5file, 'w') as h5f:
        h5f.create_group('rsp')

    solution_keys = ['S1', 'S2']
    solutions = [np.arange(3.0), np.arange(3.0, 6.0)]

    for index, solution in enumerate(solutions):
        write_rsp_full_solution_to_hdf5(str(h5file), solution, index,
                                        len(solutions))

    write_rsp_results_to_hdf5(str(h5file), {
        'rsp_type': 'tda',
        'full_solutions_keys': solution_keys,
    })

    recovered = read_results(str(h5file), 'rsp')

    assert set(recovered) == {'rsp_type', 'S1', 'S2'}
    np.testing.assert_allclose(recovered['S1'], solutions[0])
    np.testing.assert_allclose(recovered['S2'], solutions[1])


def test_write_rsp_full_solution_requires_numpy_array(tmp_path):

    if MPI.COMM_WORLD.Get_rank() != mpi_master():
        return

    h5file = Path(tmp_path) / 'rsp_solutions.h5'
    with h5py.File(h5file, 'w') as h5f:
        h5f.create_group('rsp')

    with pytest.raises(VeloxChemError, match='must be a NumPy array'):
        write_rsp_full_solution_to_hdf5(str(h5file), [1.0, 2.0], 0, 1)


def test_tda_rsp_omits_solution_keys_when_saving_is_disabled(tmp_path):

    rsp_drv, rsp_results, h5file = _run_water_rsp_roundtrip(
        tmp_path, TdaEigenSolver, save_solutions=False)

    if rsp_drv.rank == mpi_master():
        assert 'full_solutions_keys' not in rsp_results

        with h5py.File(h5file, 'r') as h5f:
            assert 'full_solutions_keys' not in h5f['rsp']
            assert 'full_solutions_matrix' not in h5f['rsp']


def test_write_results_to_hdf5_roundtrips_dict_with_non_string_keys(tmp_path):

    if MPI.COMM_WORLD.Get_rank() != mpi_master():
        return

    h5file = Path(tmp_path) / 'dict_entries_results.h5'
    with h5py.File(h5file, 'w'):
        pass

    results = {
        'polarizability_gradient': np.array([
            np.arange(18.0).reshape(3, 3, 2),
            np.arange(18.0, 36.0).reshape(3, 3, 2),
        ]),
    }

    write_results_to_hdf5(str(h5file),
                          'vib',
                          results,
                          value_label='vibrational result')

    with h5py.File(h5file, 'r') as h5f:
        polgrad_group = h5f['vib/polarizability_gradient']
        assert polgrad_group.attrs['value_type'] == 'ndarray'

    recovered = read_results(str(h5file), 'vib')
    assert recovered.keys() == results.keys()
    np.testing.assert_allclose(recovered['polarizability_gradient'],
                               results['polarizability_gradient'])


def test_read_results_roundtrips_only_requested_group(tmp_path):

    if MPI.COMM_WORLD.Get_rank() != mpi_master():
        return

    h5file = Path(tmp_path) / 'scf_results.h5'
    with h5py.File(h5file, 'w') as h5f:
        h5f.create_dataset('basis_set', data=np.bytes_(['def2-svp']))
        h5f.create_dataset('nuclear_charges', data=np.array([8, 1, 1]))

    scf_results = {
        'scf_type': 'restricted',
        'scf_energy': -75.123,
        'restart': False,
        'grid_level': 4,
        'xcfun': 'b3lyp',
        'scf_history_energy': np.array([-75.0, -75.123]),
        'scf_history_gradient_norm': np.array([1.0e-3, 1.0e-6]),
        'S': np.eye(2),
        'dipole_moment': np.array([0.1, 0.2, 0.3]),
        'F': (np.eye(2), np.eye(2) * 2.0),
        'potfile': None,
    }

    write_scf_results_to_hdf5(str(h5file), scf_results)

    recovered = read_results(str(h5file), 'scf')

    assert 'basis_set' not in recovered
    assert 'nuclear_charges' not in recovered
    assert recovered['scf_type'] == scf_results['scf_type']
    assert recovered['scf_energy'] == scf_results['scf_energy']
    assert recovered['restart'] == scf_results['restart']
    assert recovered['grid_level'] == scf_results['grid_level']
    assert recovered['xcfun'] == scf_results['xcfun']
    assert recovered['potfile'] is None
    np.testing.assert_allclose(recovered['S'], scf_results['S'])
    np.testing.assert_allclose(recovered['dipole_moment'],
                               scf_results['dipole_moment'])
    assert isinstance(recovered['F'], tuple)
    np.testing.assert_allclose(recovered['F'][0], scf_results['F'][0])
    np.testing.assert_allclose(recovered['F'][1], scf_results['F'][1])


def test_read_results_roundtrips_tda_rsp_and_preserves_legacy_solution_vectors(
        tmp_path):

    rsp_drv, rsp_results, h5file = _run_water_rsp_roundtrip(
        tmp_path, TdaEigenSolver)

    if MPI.COMM_WORLD.Get_rank() == mpi_master():
        assert 'full_solutions_keys' not in rsp_results

        recovered = read_results(h5file, 'rsp')

        assert 'eigenvectors' not in recovered
        assert 'full_solutions_matrix' not in recovered
        assert 'full_solutions_keys' not in recovered

        expected_rsp = {
            key: value
            for key, value in rsp_results.items()
            if key not in ['eigenvectors', 'full_solutions_matrix', 'full_solutions_keys']
        }
        recovered_rsp = {
            key: value
            for key, value in recovered.items() if key not in {'S1', 'S2'}
        }

        _assert_roundtrip_equal(expected_rsp, recovered_rsp)

        np.testing.assert_allclose(recovered['S1'],
                                   rsp_results['eigenvectors'][:, 0])
        np.testing.assert_allclose(recovered['S2'],
                                   rsp_results['eigenvectors'][:, 1])


@pytest.mark.skipif(MPI.COMM_WORLD.Get_size() != 1,
                    reason='runs only on a single MPI rank')
def test_read_results_roundtrips_rpa_rsp_and_preserves_legacy_solution_vectors(
        tmp_path):

    rsp_drv, rsp_results, h5file = _run_water_rsp_roundtrip(
        tmp_path, LinearResponseEigenSolver)

    assert 'full_solutions_keys' not in rsp_results

    recovered = read_results(h5file, 'rsp')

    assert 'eigenvectors_distributed' not in recovered
    assert 'full_solutions_matrix' not in recovered
    assert 'full_solutions_keys' not in recovered

    expected_rsp = {
        key: value
        for key, value in rsp_results.items()
        if key not in ['eigenvectors_distributed', 'full_solutions_matrix', 'full_solutions_keys']
    }
    recovered_rsp = {
        key: value
        for key, value in recovered.items() if key not in {'S1', 'S2'}
    }

    _assert_roundtrip_equal(expected_rsp, recovered_rsp)

    np.testing.assert_allclose(
        recovered['S1'],
        rsp_drv.get_full_solution_vector(
            rsp_results['eigenvectors_distributed'][0]))
    np.testing.assert_allclose(
        recovered['S2'],
        rsp_drv.get_full_solution_vector(
            rsp_results['eigenvectors_distributed'][1]))


@pytest.mark.skipif(MPI.COMM_WORLD.Get_size() != 1,
                    reason='runs only on a single MPI rank')
@pytest.mark.parametrize('solver_cls',
                         [ComplexResponseSolver, ComplexResponseTdaSolver])
def test_read_results_roundtrips_cpp_rsp_and_preserves_legacy_solution_vectors(
        tmp_path, solver_cls):

    rsp_drv, rsp_results, h5file = _run_cpp_rsp_roundtrip(tmp_path, solver_cls)

    assert 'full_solutions_keys' not in rsp_results

    recovered = read_results(h5file, 'rsp')

    expected_rsp_type = (
        'tdacpp' if solver_cls is ComplexResponseTdaSolver else 'cpp')
    assert rsp_results['rsp_type'] == expected_rsp_type
    assert recovered['rsp_type'] == expected_rsp_type

    assert 'solutions' not in recovered
    assert 'full_solutions_matrix' not in recovered
    assert 'full_solutions_keys' not in recovered

    expected_rsp = {
        key: value
        for key, value in rsp_results.items()
        if key not in ['solutions', 'full_solutions_matrix', 'full_solutions_keys']
    }
    spectrum = rsp_drv.get_spectrum(rsp_results, 'au')
    expected_rsp['sigma'] = np.array(spectrum['y_data'])

    recovered_solution_keys = {
        key
        for key in recovered
        if key not in expected_rsp
    }
    expected_solution_keys = {
        f'{bop}_{w:.8f}'
        for (aop, bop, w) in rsp_results['response_functions']
    }

    assert recovered_solution_keys == expected_solution_keys
    _assert_roundtrip_equal(expected_rsp,
                            {key: value for key, value in recovered.items()
                             if key in expected_rsp})

    for key in rsp_results['solutions']:
        full_vec = rsp_drv.get_full_solution_vector(
            rsp_results['solutions'][key])
        flat_key = f'{key[0]}_{key[1]:.8f}'
        np.testing.assert_allclose(recovered[flat_key], full_vec)


@pytest.mark.skipif(MPI.COMM_WORLD.Get_size() != 1,
                    reason='runs only on a single MPI rank')
@pytest.mark.parametrize(
    'solver_cls',
    [LinearResponseSolver, LinearResponseUnrestrictedSolver])
def test_read_results_roundtrips_lr_rsp_and_preserves_solution_vectors(
        tmp_path, solver_cls):

    rsp_drv, rsp_results, h5file = _run_lr_rsp_roundtrip(tmp_path, solver_cls)

    assert 'full_solutions_keys' not in rsp_results

    recovered = read_results(h5file, 'rsp')

    assert rsp_results['rsp_type'] == 'lr'
    assert recovered['rsp_type'] == 'lr'
    assert 'solutions' not in recovered
    assert 'full_solutions_matrix' not in recovered
    assert 'full_solutions_keys' not in recovered

    expected_rsp = {
        key: value
        for key, value in rsp_results.items()
        if key not in ['solutions', 'full_solutions_matrix', 'full_solutions_keys']
    }
    recovered_solution_keys = {
        key
        for key in recovered
        if key not in expected_rsp
    }
    expected_solution_keys = {
        f'{bop}_{w:.8f}' for bop, w in rsp_results['solutions']
    }

    assert recovered_solution_keys == expected_solution_keys
    _assert_roundtrip_equal(expected_rsp,
                            {key: value for key, value in recovered.items()
                             if key in expected_rsp})

    for key, solution in rsp_results['solutions'].items():
        full_vec = rsp_drv.get_full_solution_vector(solution)
        flat_key = f'{key[0]}_{key[1]:.8f}'
        np.testing.assert_allclose(recovered[flat_key], full_vec)


@pytest.mark.skipif(MPI.COMM_WORLD.Get_size() != 1,
                    reason='runs only on a single MPI rank')
def test_read_results_roundtrips_c6_rsp_and_preserves_solution_vectors(
        tmp_path):

    rsp_drv, rsp_results, h5file = _run_c6_rsp_roundtrip(tmp_path)

    recovered = read_results(h5file, 'rsp')

    assert rsp_results['rsp_type'] == 'c6'
    assert recovered['rsp_type'] == 'c6'
    assert set(rsp_results) == {
        'c6',
        'n_points',
        'response_functions',
        'solutions',
        'rsp_type',
        'w0',
    }
    assert rsp_results['n_points'] == 3
    assert rsp_results['w0'] == 0.3
    assert 'solutions' not in recovered
    assert 'full_solutions_matrix' not in recovered
    assert 'full_solutions_keys' not in recovered

    expected_rsp = {
        key: value
        for key, value in rsp_results.items()
        if key not in ['solutions', 'full_solutions_matrix', 'full_solutions_keys']
    }
    recovered_solution_keys = {
        key
        for key in recovered
        if key not in expected_rsp
    }
    expected_solution_keys = {
        f'{bop}_{iw:.8f}' for bop, iw in rsp_results['solutions']
    }

    assert recovered_solution_keys == expected_solution_keys
    _assert_roundtrip_equal(expected_rsp,
                            {key: value for key, value in recovered.items()
                             if key in expected_rsp})

    for key, solution in rsp_results['solutions'].items():
        full_vec = rsp_drv.get_full_solution_vector(solution)
        flat_key = f'{key[0]}_{key[1]:.8f}'
        np.testing.assert_allclose(recovered[flat_key], full_vec)


def test_scf_results_hdf5_roundtrip_with_water_calculation(tmp_path):

    molecule, basis = _get_water_and_basis()
    comm = MPI.COMM_WORLD
    filename = str(tmp_path / 'water_scf_results')
    filename = comm.bcast(filename, root=mpi_master())

    scf_drv = ScfRestrictedDriver()
    scf_drv.ostream.mute()
    scf_drv.filename = filename
    scf_results = scf_drv.compute(molecule, basis)

    if scf_drv.rank == mpi_master():
        h5file = Path(f'{filename}.h5')

        recovered = read_results(str(h5file), 'scf')

        _assert_roundtrip_equal(scf_results, recovered)

        assert 'basis_set' not in recovered
        assert 'nuclear_charges' not in recovered
        assert 'scf_history' not in recovered
        for key in scf_drv.history[0]:
            expected = np.array([step[key] for step in scf_drv.history])
            np.testing.assert_allclose(recovered[f'scf_history_{key}'],
                                       expected)


@pytest.mark.skipif(MPI.COMM_WORLD.Get_size() != 1,
                    reason='runs only on a single MPI rank')
def test_read_results_ignores_legacy_untyped_rsp_nto_group(tmp_path):

    h5file = tmp_path / 'legacy_rsp_nto.h5'

    with h5py.File(h5file, 'w') as h5f:
        rsp_group = h5f.create_group('rsp')
        rsp_group.attrs['value_type'] = 'dict'
        rsp_group.attrs['dict_storage'] = 'named'

        eigenvalues = rsp_group.create_dataset('eigenvalues',
                                               data=np.array([1.2, 1.5]))
        eigenvalues.attrs['value_type'] = 'ndarray'
        eigenvalues.attrs['array_data_type'] = 'float64'

        rsp_group.create_dataset('S1', data=np.array([0.3, 0.4, 0.5]))

        nto_group = rsp_group.create_group('nto')
        nto_group.create_dataset('NTO_S1_alpha_energies',
                                 data=np.array([0.1, 0.2, 0.3]))
        nto_group.create_dataset('NTO_S1_scf_type',
                                 data=np.array([b'restricted']))

    recovered = read_results(str(h5file), 'rsp')

    assert set(recovered) == {'S1', 'eigenvalues'}
    np.testing.assert_allclose(recovered['S1'], np.array([0.3, 0.4, 0.5]))
    np.testing.assert_allclose(recovered['eigenvalues'], np.array([1.2, 1.5]))


def test_opt_results_hdf5_roundtrip_with_nh3_optimization(tmp_path):

    here = Path(__file__).parent
    inpfile = str(here / 'data' / 'nh3.inp')

    task = MpiTask([inpfile, None])
    task.input_dict['scf']['checkpoint_file'] = None
    task.input_dict['method_settings']['basis'] = 'STO-3G'

    filename = str(tmp_path / 'nh3_opt_results')
    filename = task.mpi_comm.bcast(filename, root=mpi_master())
    task.input_dict['filename'] = filename

    if task.mpi_rank == mpi_master():
        task.ao_basis = MolecularBasis.read(task.molecule, 'STO-3G',
                                            ostream=None)
    else:
        task.ao_basis = None
    task.ao_basis = task.mpi_comm.bcast(task.ao_basis, root=mpi_master())

    scf_drv = ScfRestrictedDriver(task.mpi_comm, task.ostream)
    scf_drv.update_settings(task.input_dict['scf'],
                            task.input_dict['method_settings'])
    _ = scf_drv.compute(task.molecule, task.ao_basis)

    grad_drv = ScfGradientDriver(scf_drv)
    opt_drv = OptimizationDriver(grad_drv)
    opt_drv.update_settings({
        'coordsys': 'tric',
        'filename': filename,
    })
    opt_results = opt_drv.compute(task.molecule, task.ao_basis)

    if opt_drv.rank == mpi_master():
        recovered = read_results(f'{filename}.h5', 'opt')

        _assert_opt_results_roundtrip_equal(opt_results, recovered)
        assert 'xyz' not in recovered
        assert 'final_molecule' not in recovered


def test_opt_results_hdf5_roundtrip_with_real_scan(tmp_path):

    molecule = Molecule.read_xyz_string("""2
    H2
    H 0.0 0.0 0.0
    H 0.0 0.0 0.8
    """)
    filename = str(tmp_path / 'opt_scan_real')
    filename = MPI.COMM_WORLD.bcast(filename, root=mpi_master())

    ostream = OutputStream(None)
    scf_drv = ScfRestrictedDriver(MPI.COMM_WORLD, ostream)
    scf_drv.filename = filename
    grad_drv = ScfGradientDriver(scf_drv)
    opt_drv = OptimizationDriver(grad_drv)
    basis = MolecularBasis.read(molecule, 'STO-3G', ostream=None)

    opt_drv.ostream.mute()
    opt_drv.filename = filename
    opt_drv.constraints = ['scan distance 1 2 0.6 1.0 3']
    opt_drv.conv_maxiter = True

    opt_results = opt_drv.compute(molecule, basis)

    if opt_drv.rank == mpi_master():
        recovered = read_results(f'{filename}.h5', 'opt')

        _assert_opt_results_roundtrip_equal(opt_results, recovered)
        assert 'xyz' not in recovered
        assert 'final_molecule' not in recovered

        fpath = Path('scan-final.xyz')
        if fpath.is_file():
            fpath.unlink()


def test_vib_results_hdf5_roundtrip_with_water_calculation(tmp_path):

    here = Path(__file__).parent
    inpfile = str(here / 'data' / 'water_hessian_scf.inp')

    task = MpiTask([inpfile, None])
    filename = str(tmp_path / 'water_vib_results')
    filename = task.mpi_comm.bcast(filename, root=mpi_master())

    scf_drv = ScfRestrictedDriver(task.mpi_comm, task.ostream)
    scf_drv.acc_type = 'l2_c2diis'
    scf_drv.filename = filename
    _ = scf_drv.compute(task.molecule, task.ao_basis)

    vib_drv = VibrationalAnalysis(scf_drv)
    vib_drv.update_settings({}, {
        'do_ir': 'yes',
        'do_raman': 'yes',
        'numerical_hessian': 'no',
        'numerical_raman': 'no',
        'filename': filename,
    })
    vib_drv.ostream.mute()
    vib_results = vib_drv.compute(task.molecule, task.ao_basis)

    if task.mpi_rank == mpi_master():
        recovered = read_results(f'{filename}.h5', 'vib')
        _assert_roundtrip_equal(vib_results, recovered)
        assert 'basis_set' not in recovered
        assert 'nuclear_charges' not in recovered
