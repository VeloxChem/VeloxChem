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

    def test_skeleton_returns_empty_matrix(self):

        # the compute routine is a skeleton: no integrals are filled yet
        mol, bas = self.water_sto3g()
        drv = newints.OverlapDriver()
        smat = drv.compute(mol, bas, 1.0e-12)

        assert smat.number_of_blocks() == 0
