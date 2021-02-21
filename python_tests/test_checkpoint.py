from mpi4py import MPI
from pathlib import Path
import numpy as np
import unittest
import tempfile

from veloxchem.veloxchemlib import is_mpi_master
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.outputstream import OutputStream
from veloxchem.distributedarray import DistributedArray
from veloxchem.linearsolver import LinearSolver
from veloxchem.checkpoint import write_rsp_hdf5
from veloxchem.checkpoint import read_rsp_hdf5
from veloxchem.checkpoint import check_rsp_hdf5
from veloxchem.checkpoint import append_rsp_solution_hdf5
from veloxchem.checkpoint import write_distributed_focks
from veloxchem.checkpoint import read_distributed_focks
from veloxchem.checkpoint import check_distributed_focks


class TestCheckpoint(unittest.TestCase):

    def get_molecule_and_basis(self):

        mol_str = """
            O   0.0   0.0   0.0
            H   0.0   1.4   1.1
            H   0.0  -1.4   1.1
        """
        mol = Molecule.read_str(mol_str, units='bohr')
        bas = MolecularBasis.read(mol, 'aug-cc-pvdz')

        return mol, bas

    def test_rsp_checkpoint(self):

        dft_dict = {'dft_func_label': 'HF'}
        pe_dict = {'potfile_text': ''}
        ostream = OutputStream(None)

        arrays = [np.random.rand(10, 7), np.random.rand(10, 7)]
        labels = ['b', 'e2b']

        mol, bas = self.get_molecule_and_basis()

        with tempfile.TemporaryDirectory() as temp_dir:
            if not is_mpi_master():
                return

            fname = str(Path(temp_dir, 'rsp.h5'))

            # test writing
            success = write_rsp_hdf5(fname, arrays, labels, mol, bas, dft_dict,
                                     pe_dict, ostream)
            self.assertTrue(success)

            # test reading
            read_arrays = read_rsp_hdf5(fname, labels, mol, bas, dft_dict,
                                        pe_dict, ostream)
            for a, b in zip(read_arrays, arrays):
                self.assertTrue(np.max(np.abs(a - b)) < 1.0e-12)

            # test appending
            new_array = np.random.rand(10, 20)
            new_label = 'c'
            append_rsp_solution_hdf5(fname, new_label, new_array)

            # test hdf5
            arrays.append(new_array)
            labels.append(new_label)
            match = check_rsp_hdf5(fname, labels, mol, bas, dft_dict, pe_dict)
            self.assertTrue(match)

            match = check_rsp_hdf5(fname, ['bger', 'e2bger', 'c'], mol, bas,
                                   dft_dict, pe_dict)
            self.assertFalse(match)

    def test_rsp_checkpoint_distributed(self):

        dft_dict = {'dft_func_label': 'HF'}
        pe_dict = {'potfile_text': ''}

        comm = MPI.COMM_WORLD
        ostream = OutputStream(None)
        solver = LinearSolver(comm, ostream)

        labels = ['bger', 'bung', 'e2bger', 'e2bung']
        solver.dist_bger = DistributedArray(np.random.rand(10, 7), comm)
        solver.dist_bung = DistributedArray(np.random.rand(10, 7), comm)
        solver.dist_e2bger = DistributedArray(np.random.rand(10, 7), comm)
        solver.dist_e2bung = DistributedArray(np.random.rand(10, 7), comm)
        solver.nonlinear = False

        backup_data = {
            'bger': solver.dist_bger.data.copy(),
            'bung': solver.dist_bung.data.copy(),
            'e2bger': solver.dist_e2bger.data.copy(),
            'e2bung': solver.dist_e2bung.data.copy(),
        }

        mol, bas = self.get_molecule_and_basis()

        with tempfile.TemporaryDirectory() as temp_dir:
            fname = str(Path(temp_dir, 'rsp.h5'))

            solver.checkpoint_file = fname
            solver.write_checkpoint(mol, bas, dft_dict, pe_dict, labels)

            solver.read_vectors(labels)
            self.assertTrue(
                np.max(np.abs(backup_data['bger'] -
                              solver.dist_bger.data)) < 1.0e-12)
            self.assertTrue(
                np.max(np.abs(backup_data['bung'] -
                              solver.dist_bung.data)) < 1.0e-12)
            self.assertTrue(
                np.max(np.abs(backup_data['e2bger'] -
                              solver.dist_e2bger.data)) < 1.0e-12)
            self.assertTrue(
                np.max(np.abs(backup_data['e2bung'] -
                              solver.dist_e2bung.data)) < 1.0e-12)

    def test_fock_checkpoint(self):

        comm = MPI.COMM_WORLD
        ostream = OutputStream(None)

        freqs = [0.05, 0.1]
        keys = ['f_sig_xx', 'f_lamtau_yy']

        focks = {
            key:
            {w: DistributedArray(np.random.rand(81,), comm) for w in freqs}
            for key in keys
        }

        mol, bas = self.get_molecule_and_basis()

        with tempfile.TemporaryDirectory() as temp_dir:
            fname = str(Path(temp_dir, 'fock.h5'))

            # test writing
            success = write_distributed_focks(fname, focks, keys, freqs, comm,
                                              ostream)
            self.assertTrue(success)

            # test reading
            read_focks = read_distributed_focks(fname, keys, freqs, comm,
                                                ostream)
            for key in keys:
                for w in freqs:
                    a = read_focks[key][w].get_full_vector()
                    b = focks[key][w].get_full_vector()
                    if is_mpi_master(comm):
                        self.assertTrue(np.max(np.abs(a - b)) < 1.0e-12)

            # test hdf5
            if is_mpi_master(comm):

                valid_checkpoint = check_distributed_focks(fname, keys, freqs)
                self.assertTrue(valid_checkpoint)

                valid_checkpoint = check_distributed_focks(
                    fname, ['x', 'y'], freqs)
                self.assertFalse(valid_checkpoint)

                valid_checkpoint = check_distributed_focks(
                    fname, keys, [0.1, 0.2])
                self.assertFalse(valid_checkpoint)

    def test_distributed_array(self):

        comm = MPI.COMM_WORLD

        dist_array = DistributedArray(np.random.rand(81, 2), comm)
        label = 'distributed array'

        with tempfile.TemporaryDirectory() as temp_dir:
            fname = str(Path(temp_dir, 'array.h5'))

            dist_array.append_to_hdf5_file(fname, label)

            read_array = DistributedArray.read_from_hdf5_file(
                fname, label, comm)

            for col in range(dist_array.shape(1)):
                a = dist_array.get_full_vector(col)
                b = read_array.get_full_vector(col)
                if is_mpi_master(comm):
                    self.assertTrue(np.max(np.abs(a - b)) < 1.0e-12)


if __name__ == "__main__":
    unittest.main()
