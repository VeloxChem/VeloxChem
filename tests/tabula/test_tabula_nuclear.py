import numpy as np

from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem import NuclearPotentialDriver
from veloxchem.tabulalib import TabulaNuclearAttractionDriver


WATER = (
    "O 0.000  0.000 -1.000\n"
    "H 0.000  1.400 -2.100\n"
    "H 0.000 -1.400 -2.100"
)


def _vlx_numpy(matrix):
    return matrix.full_matrix().to_numpy()


def _reference(driver, mol, bas):
    # tolerate either binding argument order
    try:
        return _vlx_numpy(driver.compute(mol, bas))
    except Exception:
        return _vlx_numpy(driver.compute(bas, mol))


def _reference_external(driver, mol, bas, mags, coords):
    try:
        return _vlx_numpy(driver.compute(mol, bas, mags, coords))
    except Exception:
        return _vlx_numpy(driver.compute(mags, coords, bas, mol))


class TestTabulaNuclearAttraction:
    """Validates NuclearAttractionDriver against VeloxChem's
    NuclearPotentialDriver. Tabula matches VeloxChem's positive convention
    (Σ_N Z_N (a|1/|r−N||b), no physical minus). Agreement sits at the Boys
    floor — not bit-exact, as the two Boys implementations differ."""

    def _setup(self, basis_label):
        mol = Molecule.read_str(WATER, "au")
        bas = MolecularBasis.read(mol, basis_label, ostream=None)
        return mol, bas

    def test_water_svp(self):
        # s / p / d
        mol, bas = self._setup("def2-svp")
        tab = TabulaNuclearAttractionDriver().compute(mol, bas).to_numpy()
        ref = _reference(NuclearPotentialDriver(), mol, bas)
        assert np.allclose(tab, ref, 0.0, 1.0e-9)

    def test_water_qzvp(self):
        # exercises the d / f / g (maxL=4) kernels
        mol, bas = self._setup("def2-qzvp")
        tab = TabulaNuclearAttractionDriver().compute(mol, bas).to_numpy()
        ref = _reference(NuclearPotentialDriver(), mol, bas)
        assert np.allclose(tab, ref, 0.0, 1.0e-8)

    def test_external_charges(self):
        # the QM/MM point-charge path
        mol, bas = self._setup("def2-qzvp")
        mags = [-0.834, 0.417, 0.417]
        coords = [[1.5, 1.0, 0.5], [-1.0, 2.0, -0.5], [0.3, -1.4, 1.2]]

        tab = TabulaNuclearAttractionDriver().compute_external(mol, bas, mags, coords).to_numpy()
        ref = _reference_external(NuclearPotentialDriver(), mol, bas, mags, coords)
        assert np.allclose(tab, ref, 0.0, 1.0e-8)

    def test_auto_exact_for_small_molecule(self):
        # a small molecule (< 100 atoms): the auto default (threshold < 0) is
        # exact dense, bit-for-bit identical to an explicit threshold of 0
        mol, bas = self._setup("def2-svp")
        drv = TabulaNuclearAttractionDriver()
        auto = drv.compute(mol, bas).to_numpy()
        exact = drv.compute(mol, bas, 0.0).to_numpy()
        assert np.array_equal(auto, exact)

    def test_auto_screen_for_large_molecule(self):
        # a 100-atom stretched H chain (>= 100 atoms) triggers the auto-screen
        # default; the screened-dense matrix stays within the conservative auto
        # threshold of the exact one, and did drop far blocks
        xyz = "\n".join(f"H 0.0 0.0 {i * 2.6:.3f}" for i in range(100))
        mol = Molecule.read_str(xyz, "au")
        bas = MolecularBasis.read(mol, "def2-svp", ostream=None)
        drv = TabulaNuclearAttractionDriver()
        auto = drv.compute(mol, bas).to_numpy()        # default: auto-screen
        exact = drv.compute(mol, bas, 0.0).to_numpy()  # explicit exact dense
        assert np.allclose(auto, exact, 0.0, 1.0e-12)
        assert not np.array_equal(auto, exact)
