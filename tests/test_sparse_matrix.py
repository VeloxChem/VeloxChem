import copy

import numpy as np

from veloxchem.veloxchemlib import BasisFunction
from veloxchem.veloxchemlib import ScreenedBasisFunctionPair
from veloxchem.veloxchemlib import SparseMatrix
from veloxchem.veloxchemlib import mat_t


class TestSparseMatrix:

    def make_pair(self, la, lb, npairs):

        # a screened basis function pair whose block layout is
        # (2*la+1)*(2*lb+1) rows and npairs columns
        bra = BasisFunction([1.0], [1.0], la)
        ket = BasisFunction([1.0], [1.0], lb)
        zeros = [0.0] * npairs
        atoms = list(range(npairs))

        return ScreenedBasisFunctionPair(bra, 0, ket, 0, zeros, zeros, zeros,
                                         zeros, zeros, zeros, atoms, atoms)

    def make_pairs(self):

        # (la, lb, npairs) -> block rows x cols
        return [
            self.make_pair(0, 0, 2),  # 1 x 2
            self.make_pair(1, 1, 3),  # 9 x 3
            self.make_pair(0, 1, 4),  # 3 x 4
        ]

    def test_type(self):

        pairs = self.make_pairs()
        assert SparseMatrix(pairs, mat_t.general).type() == mat_t.general
        assert SparseMatrix(pairs, mat_t.symmetric).type() == mat_t.symmetric
        assert SparseMatrix(pairs,
                            mat_t.antisymmetric).type() == mat_t.antisymmetric

    def test_default_constructor(self):

        sm = SparseMatrix()
        assert sm.type() == mat_t.general
        assert sm.number_of_keys() == 0
        assert sm.number_of_blocks() == 0
        assert sm.keys() == []

    def test_layout(self):

        sm = SparseMatrix(self.make_pairs(), mat_t.general)
        assert sm.number_of_keys() == 3
        assert (sm.block_rows(0), sm.block_columns(0)) == (1, 2)
        assert (sm.block_rows(1), sm.block_columns(1)) == (9, 3)
        assert (sm.block_rows(2), sm.block_columns(2)) == (3, 4)

    def test_lazy_blocks(self):

        sm = SparseMatrix(self.make_pairs(), mat_t.general)

        # blocks are created lazily on first access
        assert sm.number_of_blocks() == 0
        assert not sm.has_block(1)

        blk = sm.block(1)
        assert sm.has_block(1)
        assert sm.number_of_blocks() == 1
        assert sm.keys() == [1]
        assert blk.number_of_rows() == 9
        assert blk.number_of_columns() == 3
        assert np.allclose(blk.to_numpy(), np.zeros((9, 3)))

        sm.block(0)
        assert sm.number_of_blocks() == 2
        assert sm.keys() == [0, 1]

    def test_zero(self):

        sm = SparseMatrix(self.make_pairs(), mat_t.general)
        sm.block(2)
        sm.zero()
        assert np.allclose(sm.block(2).to_numpy(), np.zeros((3, 4)))

    def test_equality(self):

        pairs = self.make_pairs()
        a = SparseMatrix(pairs, mat_t.symmetric)
        b = SparseMatrix(pairs, mat_t.symmetric)
        assert a == b

        c = SparseMatrix(pairs, mat_t.general)
        assert a != c

    def test_copy(self):

        a = SparseMatrix(self.make_pairs(), mat_t.symmetric)
        a.block(0)
        b = copy.deepcopy(a)
        assert a == b
        assert b.has_block(0)
