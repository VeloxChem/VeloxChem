import numpy as np

from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.veloxchemlib import OverlapDriver
from veloxchem.tabulalib import TabulaOverlapDriver


class TestTabulaOverlap:
    """Tests for tabula::OverlapDriver — the full late-contraction overlap
    recursion (seed ladder, contraction, single-centre MD recursion,
    Cartesian-to-spherical assembly, scatter). The matrix is validated
    element-wise against VeloxChem's OverlapDriver, which shares the AO
    ordering."""

    WATER = ("O 0.000  0.000  0.000\n"
             "H 0.000  1.430  1.107\n"
             "H 0.000 -1.430  1.107")

    def reference(self, mol, bas):

        return OverlapDriver().compute(mol, bas).full_matrix().to_numpy()

    def check(self, xyz, basis_label):

        mol = Molecule.read_str(xyz, 'au')
        bas = MolecularBasis.read(mol, basis_label, ostream=None)

        tabula_s = TabulaOverlapDriver().compute(mol, bas).to_numpy()
        vlx_s = self.reference(mol, bas)

        assert tabula_s.shape == vlx_s.shape
        assert np.allclose(tabula_s, vlx_s, 0.0, 1.0e-10)

    def test_h2_sto3g(self):

        # s functions only
        self.check("H 0.0 0.0 0.0\nH 0.0 0.0 1.4", 'STO-3G')

    def test_water_sto3g(self):

        # s and p functions
        self.check(self.WATER, 'STO-3G')

    def test_water_def2svp(self):

        # s, p and d functions
        self.check(self.WATER, 'def2-SVP')

    def test_screening_threshold(self):

        # a small screening threshold keeps the result within tolerance of
        # the unscreened overlap
        mol = Molecule.read_str(self.WATER, 'au')
        bas = MolecularBasis.read(mol, 'def2-SVP', ostream=None)

        screened = TabulaOverlapDriver().compute(mol, bas, threshold=1.0e-12)
        assert np.allclose(screened.to_numpy(), self.reference(mol, bas),
                           0.0, 1.0e-8)
