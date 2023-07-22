from mpi4py import MPI

from veloxchem.veloxchemlib import Matrices
from veloxchem.matrix import Matrix
from veloxchem.veloxchemlib import make_matrix
from veloxchem.veloxchemlib import mat_t
from veloxchem.veloxchemlib import MolecularBasis
from veloxchem.veloxchemlib import Molecule
from tester import Tester


class TestMatrices:

    def get_matrix(self):

        h2ostr = """O   0.000   0.000  -1.000
                    H   0.000   1.400  -2.100
                    H   0.000  -1.400  -2.100"""

        mol = Molecule.read_str(h2ostr, 'au')

        bas = MolecularBasis.read(mol, 'DEF2-SVP', 'basis', ostream=None)
        
        return make_matrix(bas, mat_t.symm)
        
    def test_get_matrix_number(self):

        mats = Matrices({0 : self.get_matrix(),
                         2 : self.get_matrix()})
        mats.zero()

        mat_a = self.get_matrix()
        mat_a.zero()
        
        Tester.compare_matrices(mat_a, mats.get_matrix(0))
        Tester.compare_matrices(mat_a, mats.get_matrix(2))
        
    def test_get_matrix_label(self):

        mats = Matrices({0 : self.get_matrix(),
                         1 : self.get_matrix(),
                         2 : self.get_matrix()})
        mats.zero()

        mat_a = self.get_matrix()
        mat_a.zero()
        
        Tester.compare_matrices(mat_a, mats.get_matrix('x'))
        Tester.compare_matrices(mat_a, mats.get_matrix('y'))
        Tester.compare_matrices(mat_a, mats.get_matrix('z'))
        
    def test_get_matrix_atomic_label(self):

        mats = Matrices({0 : self.get_matrix(),
                         2 : self.get_matrix(),
                         5 : self.get_matrix()})
        mats.zero()

        mat_a = self.get_matrix()
        mat_a.zero()
        
        Tester.compare_matrices(mat_a, mats.get_matrix(0, 'x'))
        Tester.compare_matrices(mat_a, mats.get_matrix(0, 'y'))
        Tester.compare_matrices(mat_a, mats.get_matrix(0, 'z'))
        
    def test_get_matrix_pair_of_atomic_labels(self):

        mats = Matrices({0 : self.get_matrix(),
                         2 : self.get_matrix(),
                         5 : self.get_matrix()})
        mats.zero()

        mat_a = self.get_matrix()
        mat_a.zero()
        
        Tester.compare_matrices(mat_a, mats.get_matrix(0, 'x', 0, 'x'))
        Tester.compare_matrices(mat_a, mats.get_matrix(0, 'x', 1, 'x'))
        Tester.compare_matrices(mat_a, mats.get_matrix(0, 'x', 0, 'y'))
        
    def test_get_matrix_triple_of_atomic_labels(self):

        mats = Matrices({0 : self.get_matrix(),
                         2 : self.get_matrix(),
                         5 : self.get_matrix()})
        mats.zero()

        mat_a = self.get_matrix()
        mat_a.zero()
        
        Tester.compare_matrices(mat_a, mats.get_matrix(0, 'x', 0, 'x', 0, 'x'))
        Tester.compare_matrices(mat_a, mats.get_matrix(0, 'x', 0, 'x', 1, 'x'))
        Tester.compare_matrices(mat_a, mats.get_matrix(0, 'x', 0, 'x', 0, 'y'))
        
    def test_get_matrix_triple_of_atomic_labels(self):

        mats = Matrices({0 : self.get_matrix(),
                         2 : self.get_matrix(),
                         5 : self.get_matrix()})
        mats.zero()

        mat_a = self.get_matrix()
        mat_a.zero()
        
        Tester.compare_matrices(mat_a, mats.get_matrix(0, 'x', 0, 'x', 0, 'x', 0, 'x'))
        Tester.compare_matrices(mat_a, mats.get_matrix(0, 'x', 0, 'x', 1, 'x', 0, 'x'))
        Tester.compare_matrices(mat_a, mats.get_matrix(0, 'x', 0, 'x', 0, 'x', 1, 'x'))
