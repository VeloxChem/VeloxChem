import numpy as np

from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.veloxchemlib import TabulaOverlapDriver, TabulaSymmetry


class TestTabulaOverlap:
    """Skeleton tests for tabula::OverlapDriver.

    The per-shell-pair primitive recursion is still a stub, so the overlap
    comes out all zeros. These tests exercise the scaffold — basis-block
    intake, the block-pair loop, the AO scatter, the matrix dimension and
    symmetry — on both s-only and s+p molecules. Once the custom recursion is
    in, the zero checks become comparisons against VeloxChem's overlap.
    """

    def test_h2_sto3g(self):

        mol = Molecule.read_molecule_string(
            "H 0.0 0.0 0.0\nH 0.0 0.0 0.74", units='angstrom')
        bas = MolecularBasis.read(mol, 'STO-3G', ostream=None)

        S = TabulaOverlapDriver().compute(mol, bas)

        # H2 / STO-3G -> 2 s-type AOs
        assert S.rows() == 2 and S.columns() == 2
        assert S.symmetry() == TabulaSymmetry.symmetric
        # the primitive recursion is stubbed -> all zeros
        assert np.allclose(S.to_numpy(), 0.0, 1.0e-13, 1.0e-13)

    def test_water_sto3g(self):

        mol = Molecule.read_molecule_string(
            "O 0.000  0.000  0.000\n"
            "H 0.000  0.757  0.587\n"
            "H 0.000 -0.757  0.587", units='angstrom')
        bas = MolecularBasis.read(mol, 'STO-3G', ostream=None)

        S = TabulaOverlapDriver().compute(mol, bas)

        # H2O / STO-3G -> O(1s,2s,2p) = 5 + 2 x H(1s) = 2 -> 7 AOs;
        # exercises the (2l+1)-component scatter of the p-block
        assert S.rows() == 7 and S.columns() == 7
        assert S.symmetry() == TabulaSymmetry.symmetric
        assert np.allclose(S.to_numpy(), 0.0, 1.0e-13, 1.0e-13)
