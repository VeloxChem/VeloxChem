from pathlib import Path

from mpi4py import MPI
import h5py
import numpy as np

from veloxchem import mpi_master
from veloxchem.resultsio import write_scf_results_to_hdf5


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
        'S': np.eye(2),
        'dipole_moment': np.array([0.1, 0.2, 0.3]),
        'F': (np.eye(2), np.eye(2) * 2.0),
        'potfile': None,
    }
    scf_history = [{
        'energy': -75.0,
        'gradient_norm': 1.0e-3,
    }, {
        'energy': -75.123,
        'gradient_norm': 1.0e-6,
    }]

    write_scf_results_to_hdf5(str(h5file), scf_results, scf_history)

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
        assert history_group.attrs['value_type'] == 'dict'
        assert history_group['energy'].attrs['value_type'] == 'ndarray'
        np.testing.assert_allclose(np.array(history_group['energy']),
                                   np.array([-75.0, -75.123]))
        np.testing.assert_allclose(np.array(history_group['gradient_norm']),
                                   np.array([1.0e-3, 1.0e-6]))
