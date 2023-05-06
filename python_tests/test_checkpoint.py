from mpi4py import MPI
from pathlib import Path
from random import choice
import numpy as np

from veloxchem.veloxchemlib import is_mpi_master
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.outputstream import OutputStream
from veloxchem.distributedarray import DistributedArray
from veloxchem.linearsolver import LinearSolver
from veloxchem.checkpoint import (write_rsp_hdf5, read_rsp_hdf5, check_rsp_hdf5,
                                  write_rsp_solution, write_distributed_focks,
                                  read_distributed_focks,
                                  check_distributed_focks)


class TestCheckpoint:

    def get_molecule_and_basis(self):

        mol_str = """
            O   0.0   0.0   0.0
            H   0.0   1.4   1.1
            H   0.0  -1.4   1.1
        """
        mol = Molecule.read_str(mol_str, units='bohr')
        bas = MolecularBasis.read(mol, 'aug-cc-pvdz', ostream=None)

        return mol, bas

    def test_rsp_checkpoint(self):

        dft_dict = {'dft_func_label': 'HF'}
        pe_dict = {'potfile_text': ''}
        ostream = OutputStream(None)

        arrays = [np.random.rand(10, 7), np.random.rand(10, 7)]
        labels = ['b', 'e2b']

        mol, bas = self.get_molecule_and_basis()

        if is_mpi_master():
            here = Path(__file__).parent
            random_string = ''.join([choice('abcdef123456') for i in range(8)])
            fpath = here / 'inputs' / f'vlx_checkpoint_{random_string}.h5'
            fname = str(fpath)

            # test writing
            success = write_rsp_hdf5(fname, arrays, labels, mol, bas, dft_dict,
                                     pe_dict, ostream)
            assert success

            # test reading
            read_arrays = read_rsp_hdf5(fname, labels, mol, bas, dft_dict,
                                        pe_dict, ostream)
            for a, b in zip(read_arrays, arrays):
                assert np.max(np.abs(a - b)) < 1.0e-12

            # test appending
            new_array = np.random.rand(10, 20)
            new_label = 'c'
            write_rsp_solution(fname, new_label, new_array)

            # test hdf5
            arrays.append(new_array)
            labels.append(new_label)
            match = check_rsp_hdf5(fname, labels, mol, bas, dft_dict, pe_dict)
            assert match

            match = check_rsp_hdf5(fname, ['bger', 'e2bger', 'c'], mol, bas,
                                   dft_dict, pe_dict)
            assert not match

            if fpath.is_file():
                fpath.unlink()

    def test_rsp_checkpoint_distributed(self):

        dft_dict = {'dft_func_label': 'HF'}
        pe_dict = {'potfile_text': ''}

        comm = MPI.COMM_WORLD
        ostream = OutputStream(None)
        solver = LinearSolver(comm, ostream)

        labels = ['bger', 'bung', 'e2bger', 'e2bung']
        solver._dist_bger = DistributedArray(np.random.rand(100, 7), comm)
        solver._dist_bung = DistributedArray(np.random.rand(100, 7), comm)
        solver._dist_e2bger = DistributedArray(np.random.rand(100, 7), comm)
        solver._dist_e2bung = DistributedArray(np.random.rand(100, 7), comm)
        solver.nonlinear = False

        backup_data = {
            'bger': solver._dist_bger.data.copy(),
            'bung': solver._dist_bung.data.copy(),
            'e2bger': solver._dist_e2bger.data.copy(),
            'e2bung': solver._dist_e2bung.data.copy(),
        }

        mol, bas = self.get_molecule_and_basis()

        here = Path(__file__).parent
        fpath = here / 'inputs' / 'rsp_with_random_values_for_test.h5'
        fname = str(fpath)

        solver.checkpoint_file = fname
        solver._write_checkpoint(mol, bas, dft_dict, pe_dict, labels)

        solver._read_checkpoint(labels)
        assert np.max(
            np.abs(backup_data['bger'] - solver._dist_bger.data)) < 1.0e-12
        assert np.max(
            np.abs(backup_data['bung'] - solver._dist_bung.data)) < 1.0e-12
        assert np.max(
            np.abs(backup_data['e2bger'] - solver._dist_e2bger.data)) < 1.0e-12
        assert np.max(
            np.abs(backup_data['e2bung'] - solver._dist_e2bung.data)) < 1.0e-12

        if is_mpi_master(comm) and fpath.is_file():
            fpath.unlink()

    def test_fock_checkpoint(self):

        comm = MPI.COMM_WORLD
        ostream = OutputStream(None)

        freqs = [0.05, 0.1]
        keys = ['f_sig_xx', 'f_lamtau_yy']
        key_freq_pairs = [(key, w) for w in freqs for key in keys]

        wrong_key_freq_pairs_1 = [(key, w) for w in freqs for key in ['x', 'y']]
        wrong_key_freq_pairs_2 = [(key, w) for w in [0.1, 0.2] for key in keys]

        dist_focks = DistributedArray(np.random.rand(81, len(key_freq_pairs)),
                                      comm)

        mol, bas = self.get_molecule_and_basis()

        here = Path(__file__).parent
        fpath = here / 'inputs' / 'fock_with_random_values_for_test.h5'
        fname = str(fpath)

        # test writing
        success = write_distributed_focks(fname, dist_focks, key_freq_pairs,
                                          comm, ostream)
        assert success

        # test reading
        read_focks = read_distributed_focks(fname, comm, ostream)

        for fock_index in range(len(key_freq_pairs)):
            a_dist = DistributedArray(dist_focks.data[:, fock_index],
                                      comm,
                                      distribute=False)
            b_dist = DistributedArray(read_focks.data[:, fock_index],
                                      comm,
                                      distribute=False)
            a_full = a_dist.get_full_vector()
            b_full = b_dist.get_full_vector()
            if is_mpi_master(comm):
                assert np.max(np.abs(a_full - b_full)) < 1.0e-12

        # test hdf5
        if is_mpi_master(comm):

            valid_checkpoint = check_distributed_focks(fname, key_freq_pairs)
            assert valid_checkpoint

            valid_checkpoint = check_distributed_focks(fname,
                                                       wrong_key_freq_pairs_1)
            assert not valid_checkpoint

            valid_checkpoint = check_distributed_focks(fname,
                                                       wrong_key_freq_pairs_2)
            assert not valid_checkpoint

            if fpath.is_file():
                fpath.unlink()

    def test_distributed_array(self):

        comm = MPI.COMM_WORLD

        dist_array = DistributedArray(np.random.rand(81, 2), comm)
        label = 'distributed array'

        here = Path(__file__).parent
        fpath = here / 'inputs' / 'array_with_random_values_for_test.h5'
        fname = str(fpath)

        dist_array.append_to_hdf5_file(fname, label)

        read_array = DistributedArray.read_from_hdf5_file(fname, label, comm)

        for col in range(dist_array.shape(1)):
            a = dist_array.get_full_vector(col)
            b = read_array.get_full_vector(col)

            if is_mpi_master(comm):
                assert np.max(np.abs(a - b)) < 1.0e-12

                if fpath.is_file():
                    fpath.unlink()
