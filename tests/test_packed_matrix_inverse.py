import numpy as np

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
