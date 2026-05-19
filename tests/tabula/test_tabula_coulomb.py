import numpy as np

from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.veloxchemlib import TwoCenterElectronRepulsionDriver
from veloxchem.tabulalib import TabulaCoulombDriver


class TestTabulaCoulomb:
    """Tests for tabula::CoulombDriver — the two-center Coulomb integral
    (a|b) = ∫∫ χ_a(r₁)·r₁₂⁻¹·χ_c(r₂). The matrix is validated element-wise
    against VeloxChem's TwoCenterElectronRepulsionDriver, which shares the
    AO ordering."""

    WATER = ("O 0.000  0.000  0.000\n"
             "H 0.000  1.430  1.107\n"
             "H 0.000 -1.430  1.107")

    def reference(self, mol, bas):

        return TwoCenterElectronRepulsionDriver().compute(mol, bas).full_matrix().to_numpy()

    def check(self, xyz, basis_label):

        mol = Molecule.read_str(xyz, 'au')
        bas = MolecularBasis.read(mol, basis_label, ostream=None)

        tabula_j = TabulaCoulombDriver().compute(mol, bas).to_numpy()
        vlx_j = self.reference(mol, bas)

        assert tabula_j.shape == vlx_j.shape
        assert np.allclose(tabula_j, vlx_j, 1.0e-9, 1.0e-9)

    def test_h2_sto3g(self):

        # s functions only
        self.check("H 0.0 0.0 0.0\nH 0.0 0.0 1.4", 'STO-3G')

    def test_water_sto3g(self):

        # s and p functions
        self.check(self.WATER, 'STO-3G')

    def test_water_def2svp(self):

        # s, p and d functions
        self.check(self.WATER, 'def2-SVP')

    def test_water_def2qzvp(self):

        # s, p, d, f and g functions — the full l = 0..4 transform range
        self.check(self.WATER, 'def2-QZVP')

    def test_screening_threshold(self):

        # a small screening threshold keeps the result within tolerance of
        # the unscreened Coulomb matrix
        mol = Molecule.read_str(self.WATER, 'au')
        bas = MolecularBasis.read(mol, 'def2-SVP', ostream=None)

        screened = TabulaCoulombDriver().compute(mol, bas, threshold=1.0e-12)
        unscreened = TabulaCoulombDriver().compute(mol, bas).to_numpy()
        assert np.allclose(screened.to_numpy(), unscreened, 0.0, 1.0e-8)
