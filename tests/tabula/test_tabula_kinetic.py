import numpy as np

from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.veloxchemlib import KineticEnergyDriver
from veloxchem.tabulalib import TabulaKineticDriver


class TestTabulaKinetic:
    """Tests for tabula::KineticDriver — the late-contraction kinetic-energy
    recursion (kinetic seed ladder, contraction, single-centre MD recursion,
    Cartesian-to-spherical assembly, scatter). The matrix is validated
    element-wise against VeloxChem's KineticEnergyDriver, which shares the AO
    ordering."""

    WATER = ("O 0.000  0.000  0.000\n"
             "H 0.000  1.430  1.107\n"
             "H 0.000 -1.430  1.107")

    def reference(self, mol, bas):

        return KineticEnergyDriver().compute(mol, bas).full_matrix().to_numpy()

    def check(self, xyz, basis_label):

        mol = Molecule.read_str(xyz, 'au')
        bas = MolecularBasis.read(mol, basis_label, ostream=None)

        tabula_t = TabulaKineticDriver().compute(mol, bas).to_numpy()
        vlx_t = self.reference(mol, bas)

        assert tabula_t.shape == vlx_t.shape
        assert np.allclose(tabula_t, vlx_t, 1.0e-9, 1.0e-9)

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
        # the unscreened kinetic-energy matrix
        mol = Molecule.read_str(self.WATER, 'au')
        bas = MolecularBasis.read(mol, 'def2-SVP', ostream=None)

        screened = TabulaKineticDriver().compute(mol, bas, threshold=1.0e-12)
        unscreened = TabulaKineticDriver().compute(mol, bas).to_numpy()
        assert np.allclose(screened.to_numpy(), unscreened, 0.0, 1.0e-8)
