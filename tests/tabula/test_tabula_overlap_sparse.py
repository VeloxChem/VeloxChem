import numpy as np

from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.veloxchemlib import OverlapDriver
from veloxchem.tabulalib import TabulaOverlapDriver


class TestTabulaOverlapSparse:
    """Tests for OverlapDriver.compute_sparse — the block-sparse overlap."""

    def _setup(self, xyz, basis_label):
        mol = Molecule.read_str(xyz, "au")
        bas = MolecularBasis.read(mol, basis_label, ostream=None)
        reference = OverlapDriver().compute(mol, bas).full_matrix().to_numpy()
        return mol, bas, reference

    def test_unscreened_matches_dense(self):

        # threshold 0 stores every atom-pair block — to_dense reproduces the
        # exact overlap, validating the block-sparse scatter and reconstruction
        mol, bas, reference = self._setup(
            "O 0.000  0.000 -1.000\n"
            "H 0.000  1.400 -2.100\n"
            "H 0.000 -1.400 -2.100", "def2-svp")

        sparse = TabulaOverlapDriver().compute_sparse(mol, bas, 0.0)
        dense = sparse.to_dense().to_numpy()

        assert dense.shape == reference.shape
        assert np.allclose(dense, reference, 0.0, 1.0e-12)

    def test_screened_within_threshold(self):

        # a stretched hydrogen chain — far-atom blocks screen out; a sound
        # screener bounds the deviation by the threshold
        xyz = "\n".join(f"H 0.0 0.0 {i * 2.4:.3f}" for i in range(16))
        mol, bas, reference = self._setup(xyz, "def2-svp")
        dimension = reference.shape[0]

        for threshold in (1.0e-12, 1.0e-9):
            sparse = TabulaOverlapDriver().compute_sparse(mol, bas, threshold)
            dense = sparse.to_dense().to_numpy()

            deviation = float(np.max(np.abs(dense - reference)))
            assert deviation <= threshold

            # screening dropped blocks — the footprint is below the dense one
            assert sparse.stored_element_count() < dimension * dimension
