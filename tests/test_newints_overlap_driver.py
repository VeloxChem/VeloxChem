import numpy as np

from veloxchem.veloxchemlib import newints
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis


class TestNewIntsOverlapDriver:

    def water_sto3g(self):

        xyz = """3

        O    0.000000    0.000000    0.000000
        H    0.000000    0.000000    0.950000
        H    0.895670    0.000000   -0.316663
        """
        mol = Molecule.read_xyz_string(xyz)
        bas = MolecularBasis.read(mol, "sto-3g", ostream=None)
        return mol, bas

    def test_compute_returns_symmetric_sparse_matrix(self):

        mol, bas = self.water_sto3g()
        drv = newints.OverlapDriver()
        smat = drv.compute(mol, bas, 1.0e-12)

        assert isinstance(smat, newints.SparseMatrix)
        assert smat.symmetry() == newints.SymmetryType.symmetric

    def test_compute_block_structure(self):

        # the compute loop builds the block structure; kernels are still stubs
        mol, bas = self.water_sto3g()
        smat = newints.OverlapDriver().compute(mol, bas, 1.0e-12)

        # water / STO-3G shells (newints idx): O[0,1,2 = s,s,p], H1[5 = s], H2[6 = s]
        expected_keys = [
            (0, 0), (0, 1), (1, 1), (2, 2),   # O-O, same atom, l == l'
            (0, 5), (1, 5), (2, 5),            # O-H1
            (0, 6), (1, 6), (2, 6),            # O-H2
            (5, 5),                            # H1-H1
            (5, 6),                            # H1-H2
            (6, 6),                            # H2-H2
        ]
        assert smat.number_of_blocks() == len(expected_keys)
        assert smat.keys() == sorted(expected_keys)

    def test_compute_block_layout(self):

        mol, bas = self.water_sto3g()
        smat = newints.OverlapDriver().compute(mol, bas, 1.0e-12)

        # diagonal (i, i) blocks packed lower-triangular; off-diagonal blocks full
        assert smat.block((0, 0)).kind == newints.Kind.lower_triangular
        assert smat.block((2, 2)).kind == newints.Kind.lower_triangular  # p-p diagonal
        assert (smat.block((2, 2)).nrows, smat.block((2, 2)).ncols) == (3, 3)
        assert smat.block((0, 1)).kind == newints.Kind.full              # s-s same atom
        assert smat.block((2, 5)).kind == newints.Kind.full              # p(O)-s(H1)
        assert (smat.block((2, 5)).nrows, smat.block((2, 5)).ncols) == (3, 1)

    def test_compute_values_are_zero_stub(self):

        # kernels are stubs -> all integral values are zero for now
        mol, bas = self.water_sto3g()
        smat = newints.OverlapDriver().compute(mol, bas, 1.0e-12)

        dense = smat.to_dense(bas).to_numpy()
        assert dense.shape == (7, 7)
        assert np.allclose(dense, 0.0)
