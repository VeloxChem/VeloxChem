import numpy as np

from veloxchem.veloxchemlib import newints
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis


def water_sto3g():
    xyz = """3

    O    0.000000    0.000000    0.000000
    H    0.000000    0.000000    0.950000
    H    0.895670    0.000000   -0.316663
    """
    mol = Molecule.read_xyz_string(xyz)
    bas = MolecularBasis.read(mol, "sto-3g", ostream=None)
    return mol, bas


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

    # --- triangular storage (symmetric / antisymmetric) ---

    def test_symmetric_stores_upper_triangle_only(self):

        m = newints.SparseMatrix(newints.SymmetryType.symmetric)
        m.add(0, 2, newints.Block(1, 3, [0.1, 0.2, 0.3]))
        assert m.contains((0, 2))
        assert not m.contains((2, 0))  # partner is implied, not stored
        assert m.number_of_blocks() == 1

    def test_symmetric_add_canonicalizes_lower_to_upper(self):

        m = newints.SparseMatrix(newints.SymmetryType.symmetric)
        # add at (2,0) with a 3x1 block -> transposed and stored at (0,2)
        m.add(2, 0, newints.Block(3, 1, [0.1, 0.2, 0.3]))
        assert m.contains((0, 2))
        assert not m.contains((2, 0))
        blk = m.block((0, 2))
        assert blk.nrows == 1 and blk.ncols == 3
        assert blk.kind == newints.Kind.full
        assert blk.data == [0.1, 0.2, 0.3]

    def test_antisymmetric_add_canonicalizes_with_negation(self):

        m = newints.SparseMatrix(newints.SymmetryType.antisymmetric)
        m.add(2, 0, newints.Block(3, 1, [0.1, 0.2, 0.3]))
        blk = m.block((0, 2))
        assert blk.nrows == 1 and blk.ncols == 3
        assert blk.data == [-0.1, -0.2, -0.3]

    def test_symmetric_diagonal_block_packed(self):

        m = newints.SparseMatrix(newints.SymmetryType.symmetric)
        full = [1.0, 2.0, 3.0,
                2.0, 4.0, 5.0,
                3.0, 5.0, 6.0]
        m.add(0, 0, newints.Block(3, 3, full))
        blk = m.block((0, 0))
        assert blk.kind == newints.Kind.lower_triangular
        assert blk.nrows == 3 and blk.ncols == 3
        # packed lower triangle (r >= c): a00, a10, a11, a20, a21, a22
        assert blk.data == [1.0, 2.0, 4.0, 3.0, 5.0, 6.0]

    def test_antisymmetric_diagonal_block_packed(self):

        m = newints.SparseMatrix(newints.SymmetryType.antisymmetric)
        full = [0.0, 7.0,
                -7.0, 0.0]
        m.add(0, 0, newints.Block(2, 2, full))
        blk = m.block((0, 0))
        assert blk.kind == newints.Kind.lower_triangular
        assert blk.data == [0.0, -7.0, 0.0]

    def test_general_stores_as_given(self):

        m = newints.SparseMatrix()  # general
        m.add(2, 0, newints.Block(1, 1, [5.0]))
        assert m.contains((2, 0))           # no canonicalization
        assert m.block((0, 0)) is None
        assert m.block((2, 0)).kind == newints.Kind.full

    # --- conversion to dense matrix in VeloxChem ordering ---

    def test_to_dense_diagonal_placement(self):

        # newints offsets for water/STO-3G: shells at 0,1,2,5,6 with l = 0,0,1,0,0
        # VeloxChem AO mapping: 0->0, 1->1, 2(p)->[4,5,6], 5(H1 s)->2, 6(H2 s)->3
        _, bas = water_sto3g()
        m = newints.SparseMatrix(newints.SymmetryType.symmetric)
        m.add(0, 0, newints.Block(1, 1, [10.0]))
        m.add(1, 1, newints.Block(1, 1, [11.0]))
        m.add(2, 2, newints.Block(3, 3, [20.0, 0.0, 0.0,
                                         0.0, 21.0, 0.0,
                                         0.0, 0.0, 22.0]))
        m.add(5, 5, newints.Block(1, 1, [13.0]))
        m.add(6, 6, newints.Block(1, 1, [14.0]))

        dense = m.to_dense(bas).to_numpy()
        assert dense.shape == (7, 7)

        expected = np.zeros(7)
        expected[0] = 10.0
        expected[1] = 11.0
        expected[2] = 13.0
        expected[3] = 14.0
        expected[4] = 20.0
        expected[5] = 21.0
        expected[6] = 22.0
        assert np.allclose(np.diag(dense), expected)
        # everything off the diagonal is zero
        assert np.allclose(dense - np.diag(np.diag(dense)), 0.0)

    def test_to_dense_symmetric_offdiagonal_mirror(self):

        _, bas = water_sto3g()
        m = newints.SparseMatrix(newints.SymmetryType.symmetric)
        # s(shell 0) - p(shell 2) block; shell0 -> AO 0, shell2 -> AO [4,5,6]
        m.add(0, 2, newints.Block(1, 3, [1.0, 2.0, 3.0]))
        dense = m.to_dense(bas).to_numpy()

        assert dense[0, 4] == 1.0 and dense[0, 5] == 2.0 and dense[0, 6] == 3.0
        assert dense[4, 0] == 1.0 and dense[5, 0] == 2.0 and dense[6, 0] == 3.0
        assert np.allclose(dense, dense.T)

    def test_to_dense_antisymmetric_offdiagonal_mirror(self):

        _, bas = water_sto3g()
        m = newints.SparseMatrix(newints.SymmetryType.antisymmetric)
        m.add(0, 2, newints.Block(1, 3, [1.0, 2.0, 3.0]))
        dense = m.to_dense(bas).to_numpy()

        assert dense[0, 4] == 1.0 and dense[0, 5] == 2.0 and dense[0, 6] == 3.0
        assert dense[4, 0] == -1.0 and dense[5, 0] == -2.0 and dense[6, 0] == -3.0
        assert np.allclose(dense, -dense.T)
