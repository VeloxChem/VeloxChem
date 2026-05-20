import numpy as np

from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.tabulalib import TabulaNuclearAttractionDriver


class TestTabulaNuclearAttractionSparse:
    """Tests for NuclearAttractionDriver.compute_sparse — the block-sparse
    nuclear-attraction matrix, validated against the driver's own dense path.
    (Lossless screening of a long-range point-source operator is a later pass;
    here the unscreened block-sparse scatter / reconstruction is checked.)"""

    def _setup(self, xyz, basis_label):
        mol = Molecule.read_str(xyz, "au")
        bas = MolecularBasis.read(mol, basis_label, ostream=None)
        dense = TabulaNuclearAttractionDriver().compute(mol, bas).to_numpy()
        return mol, bas, dense

    def test_unscreened_matches_dense(self):
        # threshold 0 stores every atom-pair block — to_dense reproduces the
        # dense matrix, validating the block-sparse scatter and reconstruction
        mol, bas, dense = self._setup(
            "O 0.000  0.000 -1.000\n"
            "H 0.000  1.400 -2.100\n"
            "H 0.000 -1.400 -2.100", "def2-svp")

        sparse = TabulaNuclearAttractionDriver().compute_sparse(mol, bas, 0.0)
        reconstructed = sparse.to_dense().to_numpy()

        assert reconstructed.shape == dense.shape
        assert np.allclose(reconstructed, dense, 0.0, 1.0e-12)

    def test_screened_within_threshold(self):
        # a stretched hydrogen chain — far-atom blocks screen out (the bra–ket
        # overlap decays); the conservative point-source bound is sound, so the
        # deviation from the unscreened matrix stays within the threshold
        xyz = "\n".join(f"H 0.0 0.0 {i * 2.6:.3f}" for i in range(16))
        mol, bas, dense = self._setup(xyz, "def2-svp")
        dimension = dense.shape[0]

        for threshold in (1.0e-9, 1.0e-6):
            sparse = TabulaNuclearAttractionDriver().compute_sparse(mol, bas, threshold)
            reconstructed = sparse.to_dense().to_numpy()

            deviation = float(np.max(np.abs(reconstructed - dense)))
            assert deviation <= threshold

            # screening dropped far blocks — the footprint is below the dense one
            assert sparse.stored_element_count() < dimension * dimension
