from pathlib import Path

from mpi4py import MPI
import h5py
import numpy as np

from veloxchem import (Molecule, MolecularBasis, OptimizationDriver, MpiTask,
                       OutputStream, ScfGradientDriver, ScfRestrictedDriver,
                       mpi_master)
from veloxchem.lreigensolver import LinearResponseEigenSolver
from veloxchem.resultsio import (read_results, write_results_to_hdf5,
                                 write_scf_results_to_hdf5)
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


def _run_water_rsp_roundtrip(tmp_path, solver_cls):

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
    if solver_cls is LinearResponseEigenSolver:
        rsp_drv.esa = True
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
        'scf_history': [{
            'energy': -75.0,
            'gradient_norm': 1.0e-3,
        }, {
            'energy': -75.123,
            'gradient_norm': 1.0e-6,
        }],
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

        history_group = scf_group['scf_history']
        assert history_group.attrs['value_type'] == 'list'
        assert history_group.attrs['length'] == 2
        assert history_group['0'].attrs['value_type'] == 'dict'
        assert history_group['0']['energy'].attrs['value_type'] == 'float'
        assert history_group['0']['gradient_norm'].attrs['value_type'] == 'float'


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


def test_write_results_to_hdf5_roundtrips_dict_with_non_string_keys(tmp_path):

    if MPI.COMM_WORLD.Get_rank() != mpi_master():
        return

    h5file = Path(tmp_path) / 'dict_entries_results.h5'
    with h5py.File(h5file, 'w'):
        pass

    results = {
        'polarizability_gradient': {
            0.0: np.arange(18.0).reshape(2, 3, 3),
            0.2: np.arange(18.0, 36.0).reshape(2, 3, 3),
        }
    }

    write_results_to_hdf5(str(h5file),
                          'vib',
                          results,
                          value_label='vibrational result')

    with h5py.File(h5file, 'r') as h5f:
        polgrad_group = h5f['vib/polarizability_gradient']
        assert polgrad_group.attrs['value_type'] == 'dict'
        assert polgrad_group.attrs['dict_storage'] == 'entries'
        assert polgrad_group.attrs['length'] == 2

    recovered = read_results(str(h5file), 'vib')
    assert recovered.keys() == results.keys()
    assert recovered['polarizability_gradient'].keys() == (
        results['polarizability_gradient'].keys())
    for key in results['polarizability_gradient']:
        np.testing.assert_allclose(recovered['polarizability_gradient'][key],
                                   results['polarizability_gradient'][key])


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
        'scf_history': [{
            'energy': -75.0,
            'gradient_norm': 1.0e-3,
        }, {
            'energy': -75.123,
            'gradient_norm': 1.0e-6,
        }],
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
        recovered = read_results(h5file, 'rsp')

        assert 'eigenvectors' not in recovered

        expected_rsp = {
            key: value
            for key, value in rsp_results.items() if key != 'eigenvectors'
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


def test_read_results_roundtrips_rpa_rsp_and_preserves_legacy_solution_vectors(
        tmp_path):

    rsp_drv, rsp_results, h5file = _run_water_rsp_roundtrip(
        tmp_path, LinearResponseEigenSolver)

    if MPI.COMM_WORLD.Get_rank() == mpi_master():
        recovered = read_results(h5file, 'rsp')

        assert 'eigenvectors_distributed' not in recovered

        expected_rsp = {
            key: value
            for key, value in rsp_results.items()
            if key != 'eigenvectors_distributed'
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
        assert recovered['scf_history'] == scf_drv.history


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
