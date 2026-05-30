from veloxchem.veloxchemlib import newints


class TestSparseMatrix:

    def test_default_is_general_and_empty(self):

        m = newints.SparseMatrix()
        assert m.symmetry() == newints.SymmetryType.general
        assert m.number_of_blocks() == 0

    def test_symmetry_ctor_and_setter(self):

        m = newints.SparseMatrix(newints.SymmetryType.symmetric)
        assert m.symmetry() == newints.SymmetryType.symmetric
        m.set_symmetry(newints.SymmetryType.antisymmetric)
        assert m.symmetry() == newints.SymmetryType.antisymmetric

    def test_add_and_block(self):

        m = newints.SparseMatrix()
        # an s-p block: nrows = 2*0+1 = 1, ncols = 2*1+1 = 3
        m.add(0, 2, newints.Block(1, 3, [0.1, 0.2, 0.3]))
        assert m.number_of_blocks() == 1
        assert m.contains((0, 2))
        assert not m.contains((2, 0))
        blk = m.block((0, 2))
        assert blk.nrows == 1
        assert blk.ncols == 3
        assert blk.data == [0.1, 0.2, 0.3]

    def test_add_by_key(self):

        m = newints.SparseMatrix()
        m.add((1, 4), newints.Block(3, 5, [0.0] * 15))
        assert m.contains((1, 4))
        assert m.block((1, 4)).ncols == 5

    def test_absent_block_is_none(self):

        m = newints.SparseMatrix()
        assert m.block((5, 5)) is None

    def test_add_overwrites(self):

        m = newints.SparseMatrix()
        m.add(1, 1, newints.Block(1, 1, [1.0]))
        m.add(1, 1, newints.Block(1, 1, [2.0]))
        assert m.number_of_blocks() == 1
        assert m.block((1, 1)).data == [2.0]

    def test_keys_are_ascending(self):

        m = newints.SparseMatrix()
        m.add(2, 0, newints.Block(5, 1, [0.0] * 5))
        m.add(0, 0, newints.Block(1, 1, [1.0]))
        m.add(0, 1, newints.Block(1, 3, [0.0] * 3))
        assert m.keys() == [(0, 0), (0, 1), (2, 0)]

    def test_zero_keeps_structure(self):

        m = newints.SparseMatrix()
        m.add(0, 0, newints.Block(2, 2, [1.0, 2.0, 3.0, 4.0]))
        m.zero()
        assert m.number_of_blocks() == 1
        assert m.block((0, 0)).data == [0.0, 0.0, 0.0, 0.0]

    def test_block_reference_mutation_persists(self):

        m = newints.SparseMatrix()
        m.add(0, 0, newints.Block(1, 1, [1.0]))
        blk = m.block((0, 0))
        blk.data = [9.0]
        assert m.block((0, 0)).data == [9.0]

    def test_equality(self):

        a = newints.SparseMatrix(newints.SymmetryType.symmetric)
        b = newints.SparseMatrix(newints.SymmetryType.symmetric)
        a.add(0, 0, newints.Block(1, 1, [1.0]))
        b.add(0, 0, newints.Block(1, 1, [1.0]))
        assert a == b

        # differing symmetry breaks equality
        b.set_symmetry(newints.SymmetryType.general)
        assert not (a == b)

        # differing block values break equality
        b.set_symmetry(newints.SymmetryType.symmetric)
        b.add(0, 0, newints.Block(1, 1, [2.0]))
        assert not (a == b)
