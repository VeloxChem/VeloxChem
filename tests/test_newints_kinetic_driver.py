import numpy as np

from veloxchem.veloxchemlib import newints
from veloxchem.veloxchemlib import KineticEnergyDriver as RefKineticEnergyDriver
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis


class TestNewIntsKineticEnergyDriver:

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
        drv = newints.KineticEnergyDriver()
        smat = drv.compute(mol, bas, 1.0e-12)

        assert isinstance(smat, newints.SparseMatrix)
        assert smat.symmetry() == newints.SymmetryType.symmetric

    def test_compute_block_structure(self):

        mol, bas = self.water_sto3g()
        smat = newints.KineticEnergyDriver().compute(mol, bas, 1.0e-12)

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
        smat = newints.KineticEnergyDriver().compute(mol, bas, 1.0e-12)

        # diagonal (i, i) blocks packed lower-triangular; off-diagonal blocks full
        assert smat.block((0, 0)).kind == newints.Kind.lower_triangular
        assert smat.block((2, 2)).kind == newints.Kind.lower_triangular  # p-p diagonal
        assert (smat.block((2, 2)).nrows, smat.block((2, 2)).ncols) == (3, 3)
        assert smat.block((0, 1)).kind == newints.Kind.full              # s-s same atom
        assert smat.block((2, 5)).kind == newints.Kind.full              # p(O)-s(H1)
        assert (smat.block((2, 5)).nrows, smat.block((2, 5)).ncols) == (3, 1)

    def test_matches_reference_driver(self):

        # same-atom diagonal value + two-center dispatch (incl. the transpose path
        # for l_a < l_b) must reproduce the legacy driver to machine precision.
        mol, bas = self.water_sto3g()
        tnew = newints.KineticEnergyDriver().compute(mol, bas, 1.0e-12).to_dense(bas).to_numpy()
        tref = RefKineticEnergyDriver().compute(mol, bas).to_numpy()

        assert tnew.shape == (7, 7)
        assert np.allclose(tnew, tref)
        # the two-center O-H block is computed, not zero
        assert tnew[0, 2] != 0.0

    def test_matches_reference_high_angular_momentum(self):

        # exercise the d/f/g off-diagonal kernels and the l_a<l_b transpose path
        # via the (l_a, l_b) dispatch (l <= 4 is the supported range)
        xyz = """2

        N    0.000000    0.000000    0.000000
        N    0.000000    0.000000    1.100000
        """
        mol = Molecule.read_xyz_string(xyz)
        bas = MolecularBasis.read(mol, "cc-pvqz", ostream=None)
        tnew = newints.KineticEnergyDriver().compute(mol, bas, 1.0e-14).to_dense(bas).to_numpy()
        tref = RefKineticEnergyDriver().compute(mol, bas).to_numpy()

        assert np.allclose(tnew, tref, atol=1.0e-10)

    def test_screener_filters_offdiagonal_blocks(self):

        # the screener applies only to different-atom pairs; a huge threshold
        # screens all of them out, leaving only the same-atom blocks
        mol, bas = self.water_sto3g()
        smat = newints.KineticEnergyDriver().compute(mol, bas, 1.0e10)

        same_atom = [(0, 0), (0, 1), (1, 1), (2, 2),  # O
                     (5, 5),                           # H1
                     (6, 6)]                           # H2
        assert smat.keys() == sorted(same_atom)

    def test_small_threshold_keeps_all_blocks(self):

        # a tiny threshold screens nothing -> full block structure (13 blocks)
        mol, bas = self.water_sto3g()
        smat = newints.KineticEnergyDriver().compute(mol, bas, 1.0e-12)
        assert smat.number_of_blocks() == 13
