import numpy as np
import pytest

from veloxchem.veloxchemlib import PackedMatrix, mat_t


class TestPackedMatrixInverse:
    """Tests the inversion of a symmetric matrix in the packed format."""

    def _to_packed(self, dense):
        """Stores a dense symmetric matrix in the packed format."""

        ndim = dense.shape[0]

        matrix = PackedMatrix(ndim, ndim, mat_t.symmetric)

        matrix.from_numpy(np.ascontiguousarray(dense, dtype=np.float64))

        return matrix

    def _random_spd(self, ndim, seed):
        """Makes a random symmetric positive definite matrix."""

        rng = np.random.default_rng(seed)

        amat = rng.standard_normal((ndim, ndim))

        return amat @ amat.T + ndim * np.eye(ndim)

    def test_inverse_is_symmetric_and_exact(self):

        for ndim in [1, 2, 3, 7, 32, 65]:

            dense = self._random_spd(ndim, seed=ndim)

            matrix = self._to_packed(dense)

            inverse = matrix.invert()

            assert inverse.number_of_rows() == ndim
            assert inverse.number_of_columns() == ndim
            assert inverse.get_type() == mat_t.symmetric

            computed = inverse.to_numpy()

            expected = np.linalg.inv(dense)

            assert np.allclose(computed, expected, rtol=1.0e-10, atol=1.0e-12)

            # the inverse must be symmetric to the last bit, as only one
            # triangle of it is stored

            assert np.array_equal(computed, computed.T)

    def test_residual_of_inverse(self):

        ndim = 128

        dense = self._random_spd(ndim, seed=7)

        inverse = self._to_packed(dense).invert().to_numpy()

        residual = dense @ inverse - np.eye(ndim)

        assert np.max(np.abs(residual)) < 1.0e-10

    def test_indefinite_matrix(self):
        """The Cholesky factorization fails and the inversion falls back."""

        ndim = 24

        rng = np.random.default_rng(11)

        amat = rng.standard_normal((ndim, ndim))

        dense = amat + amat.T

        # a symmetric matrix with a random sign pattern of its eigenvalues is
        # invertible but is not positive definite

        assert np.min(np.linalg.eigvalsh(dense)) < 0.0

        inverse = self._to_packed(dense).invert().to_numpy()

        assert np.allclose(inverse, np.linalg.inv(dense), rtol=1.0e-8, atol=1.0e-10)

    def test_cholesky_inverse(self):
        """The inverted Cholesky factor, which is what the fitting needs."""

        for ndim in [1, 2, 3, 7, 32, 65]:

            dense = self._random_spd(ndim, seed=ndim)

            factor = self._to_packed(dense).cholesky_inverse()

            assert factor.get_type() == mat_t.lower_triangular

            # the packed storage is the triangle, not the square

            assert factor.number_of_elements() == ndim * (ndim + 1) // 2

            computed = factor.to_numpy()

            # nothing above the diagonal, rather than the mirrored elements

            assert np.array_equal(computed, np.tril(computed))

            assert np.allclose(computed, np.linalg.inv(np.linalg.cholesky(dense)),
                               rtol=1.0e-10, atol=1.0e-12)

    def test_cholesky_inverse_gives_the_inverse(self):
        """L inverted, transposed, times itself is the inverted matrix. This is the
        property the resolution of the identity rests on."""

        for ndim in [3, 17, 64]:

            dense = self._random_spd(ndim, seed=ndim + 1)

            factor = self._to_packed(dense).cholesky_inverse().to_numpy()

            assert np.allclose(factor.T @ factor, np.linalg.inv(dense),
                               rtol=1.0e-9, atol=1.0e-11)

    def test_inverse_square_root(self):
        """Multiplied by itself it is the inverse, and it is symmetric."""

        for ndim in [1, 2, 3, 7, 32, 65]:

            dense = self._random_spd(ndim, seed=ndim + 2)

            root = self._to_packed(dense).inverse_square_root(1.0e-12)

            assert root.get_type() == mat_t.symmetric

            computed = root.to_numpy()

            # symmetric to the last bit, as one triangle is stored

            assert np.array_equal(computed, computed.T)

            assert np.allclose(computed @ computed, np.linalg.inv(dense),
                               rtol=1.0e-9, atol=1.0e-11)

    def test_inverse_square_root_drops_null_directions(self):
        """A matrix which is not of full rank has no Cholesky factor to invert, and
        this is what the resolution of the identity uses in its place."""

        ndim, defect = 30, 4

        rng = np.random.default_rng(97)

        amat = rng.standard_normal((ndim, ndim))

        eigvals, eigvecs = np.linalg.eigh(amat @ amat.T + ndim * np.eye(ndim))

        # push some directions far below the threshold

        eigvals[:defect] = 1.0e-16

        dense = eigvecs @ np.diag(eigvals) @ eigvecs.T
        dense = 0.5 * (dense + dense.T)

        root = self._to_packed(dense).inverse_square_root(1.0e-12).to_numpy()

        assert np.array_equal(root, root.T)

        # the directions which were dropped are gone from the result

        assert np.linalg.matrix_rank(root) == ndim - defect

        # and it is the inverse on the directions which remain, which is to say
        # that the matrix times it is the projector onto them

        projector = root @ dense @ root

        assert np.max(np.abs(projector @ projector - projector)) < 1.0e-10

        assert np.allclose(root @ root, np.linalg.pinv(dense, rcond=1.0e-12),
                           rtol=1.0e-8, atol=1.0e-10)

    def test_inverse_square_root_keeps_everything_above_the_threshold(self):
        """Nothing is dropped from a well conditioned matrix."""

        ndim = 24

        dense = self._random_spd(ndim, seed=5)

        root = self._to_packed(dense).inverse_square_root(1.0e-12).to_numpy()

        assert np.linalg.matrix_rank(root) == ndim

    def test_cholesky_inverse_refuses_an_indefinite_matrix(self):
        """It has to be an exception. The resolution of the identity catches it and
        inverts the square root instead, which it could not do if the interpreter
        were ended."""

        ndim = 12

        rng = np.random.default_rng(83)

        amat = rng.standard_normal((ndim, ndim))

        dense = amat + amat.T

        assert np.min(np.linalg.eigvalsh(dense)) < 0.0

        with pytest.raises(RuntimeError, match="not positive definite"):
            self._to_packed(dense).cholesky_inverse()
