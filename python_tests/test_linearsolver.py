from mpi4py import MPI
import numpy as np
import time as tm
import unittest

from veloxchem.veloxchemlib import denmat
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.outputstream import OutputStream
from veloxchem.aodensitymatrix import AODensityMatrix
from veloxchem.aofockmatrix import AOFockMatrix
from veloxchem.linearsolver import LinearSolver


class TestLinearSolver(unittest.TestCase):

    def get_molecule_and_basis(self):

        mol_str = """
            O   0.0   0.0   0.0
            H   0.0   1.4   1.1
            H   0.0  -1.4   1.1
        """
        mol = Molecule.read_str(mol_str, units='bohr')
        bas = MolecularBasis.read(mol, 'aug-cc-pvdz')

        return mol, bas

    def test_comp_fock_split_comm(self):

        comm = MPI.COMM_WORLD
        ostream = OutputStream(None)

        mol, bas = self.get_molecule_and_basis()
        nao = bas.get_dimensions_of_basis(mol)

        dmat = np.diag(np.ones(nao))
        dens = AODensityMatrix([dmat], denmat.rest)
        fock = AOFockMatrix(dens)

        eri_dict = {'screening': None}
        dft_dict = {'molgrid': None, 'gs_density': None}
        pe_dict = {'V_es': None, 'pe_drv': None}

        solver = LinearSolver(comm, ostream)
        solver.comp_lr_fock_split_comm(fock, dens, mol, bas, eri_dict, dft_dict,
                                       pe_dict)

        self.assertEqual(fock.alpha_to_numpy(0).shape, dmat.shape)

    def test_need_graceful_exit(self):

        comm = MPI.COMM_WORLD
        ostream = OutputStream(None)

        solver = LinearSolver(comm, ostream)
        solver.maximum_hours = 1.0
        solver.program_start_time = tm.time()

        self.assertTrue(solver.need_graceful_exit(1.5))
        self.assertFalse(solver.need_graceful_exit(0.5))


if __name__ == "__main__":
    unittest.main()
