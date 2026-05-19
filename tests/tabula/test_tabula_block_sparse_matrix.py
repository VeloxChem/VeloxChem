import numpy as np

from veloxchem.veloxchemlib import (
    TabulaBlockSparseMatrix,
    TabulaMixedPrecisionBlockSparseMatrix,
)


class TestTabulaBlockSparseMatrix:
    """Tests for tabula::BlockSparseMatrix and
    tabula::MixedPrecisionBlockSparseMatrix — Tabula's block-sparse storage."""

    def make_two_group_matrix(self):
        """A 4x4 matrix over two AO groups of 2 orbitals each. Global AO order
        is interleaved (group 0 -> AOs 0,2; group 1 -> AOs 1,3) to exercise
        the non-contiguous group mapping. All three group pairs are stored."""

        group_global_ao = [[0, 2], [1, 3]]
        group_pairs = [(0, 0), (1, 0), (1, 1)]
        return TabulaBlockSparseMatrix(4, group_global_ao, group_pairs)

    def test_construction(self):

        mat = self.make_two_group_matrix()
        assert mat.dimension() == 4
        assert mat.number_of_groups() == 2
        assert mat.number_of_blocks() == 3
        # 3 blocks of 2x2 -> 12 stored scalars
        assert mat.stored_element_count() == 12

    def test_block_descriptors(self):

        mat = self.make_two_group_matrix()
        b0, b1, b2 = (mat.block(i) for i in range(3))
        assert (b0.group_a, b0.group_b, b0.offset) == (0, 0, 0)
        assert (b1.group_a, b1.group_b, b1.offset) == (1, 0, 4)
        assert (b2.group_a, b2.group_b, b2.offset) == (1, 1, 8)
        assert b1.row_count == 2 and b1.column_count == 2

    def test_block_value_access(self):

        mat = self.make_two_group_matrix()
        mat.set_value(1, 0, 1, 3.5)
        assert mat.value(1, 0, 1) == 3.5
        # an untouched element stays zero
        assert mat.value(1, 1, 0) == 0.0

    def test_to_dense_reconstruction(self):

        mat = self.make_two_group_matrix()
        # diagonal block (0,0): group 0 AOs are 0 and 2
        mat.set_value(0, 0, 0, 1.0)
        mat.set_value(0, 0, 1, 2.0)
        mat.set_value(0, 1, 0, 2.0)
        mat.set_value(0, 1, 1, 3.0)
        # off-diagonal block (1,0): rows = group 1 AOs (1,3), cols = group 0 (0,2)
        mat.set_value(1, 0, 0, 4.0)
        mat.set_value(1, 1, 1, 5.0)

        dense = mat.to_dense().to_numpy()
        assert dense.shape == (4, 4)
        # diagonal block scattered to AOs {0,2}
        assert dense[0, 0] == 1.0 and dense[0, 2] == 2.0
        assert dense[2, 0] == 2.0 and dense[2, 2] == 3.0
        # off-diagonal block (1,0) and its transpose
        assert dense[1, 0] == 4.0 and dense[0, 1] == 4.0
        assert dense[3, 2] == 5.0 and dense[2, 3] == 5.0
        # the reconstructed matrix is symmetric
        assert np.allclose(dense, dense.T, 1.0e-13, 1.0e-13, False)

    def test_mixed_precision_classification(self):

        mat = self.make_two_group_matrix()
        # diagonal blocks 0 (0,0) and 2 (1,1): large values
        mat.set_value(0, 0, 0, 10.0)
        mat.set_value(2, 0, 0, 10.0)
        # off-diagonal block 1 (1,0): all tiny -> below threshold
        mat.set_value(1, 0, 0, 1.0e-9)

        mixed = TabulaMixedPrecisionBlockSparseMatrix(mat, 1.0e-6)
        assert mixed.precision_threshold() == 1.0e-6
        # the two diagonal blocks stay double; the small off-diagonal demotes
        assert mixed.single_block_count() == 1
        assert mixed.double_block_count() == 2
        # the demoted block is the off-diagonal one
        single = [mixed.block(i) for i in range(3) if mixed.block(i).is_single_precision]
        assert len(single) == 1
        assert single[0].group_a != single[0].group_b

    def test_mixed_precision_diagonal_never_demotes(self):

        mat = self.make_two_group_matrix()
        # even a tiny diagonal block must stay double
        mat.set_value(0, 0, 0, 1.0e-12)
        mixed = TabulaMixedPrecisionBlockSparseMatrix(mat, 1.0e-6)
        for i in range(mixed.number_of_blocks()):
            blk = mixed.block(i)
            if blk.group_a == blk.group_b:
                assert not blk.is_single_precision

    def test_mixed_precision_footprint(self):

        mat = self.make_two_group_matrix()
        mat.set_value(1, 0, 0, 1.0e-9)  # demote the single off-diagonal block
        mixed = TabulaMixedPrecisionBlockSparseMatrix(mat, 1.0e-6)
        # 12 scalars total; one 2x2 block (4) is float, two (8) are double
        assert mixed.stored_element_count() == 12
        assert mixed.stored_byte_count() == 8 * 8 + 4 * 4

    def test_mixed_precision_to_dense(self):

        mat = self.make_two_group_matrix()
        mat.set_value(0, 0, 0, 2.0)
        mat.set_value(1, 0, 1, 0.5)
        mat.set_value(2, 1, 1, 7.0)
        mixed = TabulaMixedPrecisionBlockSparseMatrix(mat, 1.0e-6)
        # the dense reconstruction matches the source block-sparse matrix
        assert np.allclose(mixed.to_dense().to_numpy(),
                           mat.to_dense().to_numpy(), 1.0e-7, 1.0e-7, False)
