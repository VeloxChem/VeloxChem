import numpy as np
import pytest

from veloxchem.diis import Diis


class TestDiis:
    """
    Unit tests for the DIIS accelerator.
    """

    @staticmethod
    def _reference_weights(error_vectors):
        """
        Independently solves the augmented DIIS system for the weights.
        """

        n_vecs = len(error_vectors)

        bmat = np.zeros((n_vecs + 1, n_vecs + 1))
        for i in range(n_vecs):
            for j in range(n_vecs):
                bmat[i, j] = np.vdot(error_vectors[i], error_vectors[j])
        bmat[n_vecs, :n_vecs] = -1.0
        bmat[:n_vecs, n_vecs] = -1.0
        bmat[n_vecs, n_vecs] = 0.0

        bvec = np.zeros(n_vecs + 1)
        bvec[n_vecs] = -1.0

        return np.linalg.solve(bmat, bvec)[:n_vecs]

    def test_get_projected_fock(self):

        fa = np.array([[1.2, 0.3], [0.3, 0.4]])
        fb = np.array([[0.8, 0.1], [0.1, 0.2]])
        da = np.array([[1.0, 0.0], [0.0, 0.0]])
        db = np.array([[0.3, 0.0], [0.0, 0.1]])
        s = np.eye(2)

        projected = Diis.get_projected_fock(fa, fb, da, db, s)

        f0 = 0.5 * (fa + fb)
        inactive = np.matmul(s, db)
        active = np.matmul(s, da - db)
        virtual = np.eye(2) - np.matmul(s, da)
        expected = f0 + np.linalg.multi_dot([inactive, fb - f0, active.T])
        expected += np.linalg.multi_dot([active, fa - f0, virtual.T])
        expected += (expected - f0).T

        assert np.allclose(projected, expected)
        assert np.allclose(projected, projected.T)

    def test_compute_weights_and_effective_fock(self):

        fock_mats = [
            np.array([[1.0, 0.2], [0.2, 0.4]]),
            np.array([[0.8, 0.1], [0.1, 0.3]]),
            np.array([[0.6, 0.05], [0.05, 0.2]]),
        ]
        error_vectors = [
            np.array([[1.0, 0.0], [0.0, 0.0]]),
            np.array([[1.0, 1.0], [0.0, 0.0]]),
            np.array([[0.0, 0.0], [1.0, 1.0]]),
        ]
        den_mat = (np.array([[1.0, 0.0], [0.0, 0.0]]),)
        ovl_mat = np.eye(2)

        diis = Diis(max_err_vecs=5, diis_thresh=1.0, scf_type='restricted')
        for fock_mat, e_vec in zip(fock_mats, error_vectors):
            diis.store_diis_data((fock_mat,), den_mat, ovl_mat, e_vec, 0.0)

        weights = diis.compute_weights()
        ref_weights = self._reference_weights(error_vectors)

        assert np.allclose(weights, ref_weights)
        assert weights.sum() == pytest.approx(1.0)

        expected = sum(w * f for w, f in zip(weights, fock_mats))
        effective = diis.get_effective_fock((fock_mats[-1],))

        assert isinstance(effective, tuple)
        assert len(effective) == 1
        assert np.allclose(effective[0], expected)

    def test_effective_fock_return_types(self):

        fock_a = np.array([[1.0, 0.2], [0.2, 0.4]])
        fock_b = np.array([[0.8, 0.1], [0.1, 0.3]])
        den_a = np.array([[1.0, 0.0], [0.0, 0.0]])
        den_b = np.array([[0.3, 0.0], [0.0, 0.1]])
        ovl = np.eye(2)

        err_1 = np.array([[0.1, 0.0], [0.0, -0.1]])
        err_2 = np.array([[0.0, 0.2], [0.2, 0.0]])

        # Empty subspace: fall back to the current Fock/Kohn-Sham matrices.
        diis = Diis(max_err_vecs=3, diis_thresh=1.0, scf_type='restricted')
        effective = diis.get_effective_fock((fock_a,))
        assert isinstance(effective, tuple)
        assert len(effective) == 1
        assert np.allclose(effective[0], fock_a)

        scf_data = [
            ('restricted', (fock_a,), (den_a,), 1),
            ('unrestricted', (fock_a, fock_b), (den_a, den_b), 2),
            ('restricted_openshell', (fock_a, fock_b), (den_a, den_b), 1),
        ]

        for scf_type, fock_mat, den_mat, n_matrices in scf_data:

            if scf_type == 'unrestricted':
                e_vec_1 = np.vstack((err_1, err_2))
                e_vec_2 = np.vstack((err_2, err_1))
            else:
                e_vec_1 = err_1
                e_vec_2 = err_2

            diis = Diis(max_err_vecs=3, diis_thresh=1.0, scf_type=scf_type)

            # Single stored error vector: return a tuple.
            diis.store_diis_data(fock_mat, den_mat, ovl, e_vec_1, 0.0)
            effective = diis.get_effective_fock(fock_mat)
            assert isinstance(effective, tuple)
            assert len(effective) == n_matrices

            # Multiple stored error vectors: return a tuple as well.
            diis.store_diis_data(fock_mat, den_mat, ovl, e_vec_2, 0.0)
            effective = diis.get_effective_fock(fock_mat)
            assert isinstance(effective, tuple)
            assert len(effective) == n_matrices

    def test_clear_resets_for_reuse(self):

        fock_a = np.array([[1.0, 0.2], [0.2, 0.4]])
        den_a = np.array([[1.0, 0.0], [0.0, 0.0]])
        ovl = np.eye(2)
        err = np.array([[0.1, 0.0], [0.0, -0.1]])

        diis = Diis(max_err_vecs=3, diis_thresh=1.0, scf_type='restricted')
        diis.store_diis_data((fock_a,), (den_a,), ovl, err, 0.0)
        assert len(diis.error_vectors) == 1

        diis.clear()

        # Constructor settings survive the reset and the B matrix is
        # reallocated, so the instance remains usable.
        assert diis.max_err_vecs == 3
        assert diis.diis_thresh == 1.0
        assert diis.scf_type == 'restricted'
        assert diis.b_matrix.shape == (3, 3)
        assert np.allclose(diis.b_matrix, 0.0)
        assert len(diis.error_vectors) == 0
        assert len(diis.fock_matrices) == 0

        diis.store_diis_data((fock_a,), (den_a,), ovl, err, 0.0)
        effective = diis.get_effective_fock((fock_a,))
        assert isinstance(effective, tuple)
        assert len(effective) == 1

    def test_b_matrix_matches_brute_force_after_eviction(self):

        max_err_vecs = 3
        diis = Diis(max_err_vecs, 1.0, 'restricted')
        den_mat = (np.array([[1.0, 0.0], [0.0, 0.0]]),)
        ovl_mat = np.eye(2)

        stored = []
        for idx in range(2 * max_err_vecs):
            e_vec = np.array([[1.0 + idx, 0.1 * idx],
                              [0.2 * idx, -0.5 - 0.5 * idx]])
            fock_mat = np.array([[1.0 + 0.1 * idx, 0.2],
                                 [0.2, 0.4 + 0.1 * idx]])

            if len(stored) == max_err_vecs:
                stored.pop(0)
            stored.append(e_vec)

            diis.store_diis_data((fock_mat,), den_mat, ovl_mat, e_vec, 0.0)

            assert len(diis.error_vectors) == len(stored)
            for i in range(len(stored)):
                for j in range(len(stored)):
                    ref = np.vdot(stored[i], stored[j])
                    assert diis.b_matrix[i, j] == pytest.approx(ref)

        # The oldest vector must have been evicted from the window.
        assert len(stored) == max_err_vecs
        assert np.allclose(diis.error_vectors[0], stored[0])

    def test_degenerate_error_vectors_give_valid_weights(self):

        fock_a = np.array([[1.0, 0.2], [0.2, 0.4]])
        fock_b = np.array([[0.8, 0.1], [0.1, 0.3]])
        den_mat = (np.array([[1.0, 0.0], [0.0, 0.0]]),)
        ovl_mat = np.eye(2)
        e_vec = np.array([[0.1, 0.0], [0.0, -0.1]])

        diis = Diis(max_err_vecs=3, diis_thresh=1.0, scf_type='restricted')
        diis.store_diis_data((fock_a,), den_mat, ovl_mat, e_vec, 0.0)
        diis.store_diis_data((fock_b,), den_mat, ovl_mat, e_vec, 0.0)

        # Duplicate error vectors make the augmented system singular, so
        # compute_weights falls back to the pseudoinverse. The weights must
        # still be finite and normalized.
        weights = diis.compute_weights()

        assert np.all(np.isfinite(weights))
        assert weights.sum() == pytest.approx(1.0)

        expected = weights[0] * fock_a + weights[1] * fock_b
        effective = diis.get_effective_fock((fock_b,))

        assert np.allclose(effective[0], expected)

    def test_single_error_vector_returns_stored_matrices(self):

        fock_a = np.array([[1.0, 0.2], [0.2, 0.4]])
        fock_b = np.array([[0.8, 0.1], [0.1, 0.3]])
        den_a = np.array([[1.0, 0.0], [0.0, 0.0]])
        den_b = np.array([[0.3, 0.0], [0.0, 0.1]])
        ovl_mat = np.eye(2)
        e_vec = np.array([[0.1, 0.0], [0.0, -0.1]])

        # Restricted: the stored Fock matrix is returned, not the current one.
        diis = Diis(max_err_vecs=3, diis_thresh=1.0, scf_type='restricted')
        diis.store_diis_data((fock_a,), (den_a,), ovl_mat, e_vec, 0.0)
        effective = diis.get_effective_fock((fock_b,))

        assert len(effective) == 1
        assert np.allclose(effective[0], fock_a)

        # Unrestricted: both stored spin Fock matrices are returned.
        diis = Diis(max_err_vecs=3, diis_thresh=1.0, scf_type='unrestricted')
        diis.store_diis_data((fock_a, fock_b), (den_a, den_b), ovl_mat,
                             np.vstack((e_vec, 2.0 * e_vec)), 0.0)
        effective = diis.get_effective_fock((fock_b, fock_a))

        assert len(effective) == 2
        assert np.allclose(effective[0], fock_a)
        assert np.allclose(effective[1], fock_b)

        # Restricted open-shell: the stored projected Fock matrix is returned.
        diis = Diis(max_err_vecs=3,
                    diis_thresh=1.0,
                    scf_type='restricted_openshell')
        diis.store_diis_data((fock_a, fock_b), (den_a, den_b), ovl_mat, e_vec,
                             0.0)
        effective = diis.get_effective_fock((fock_b, fock_a))
        expected = Diis.get_projected_fock(fock_a, fock_b, den_a, den_b,
                                           ovl_mat)

        assert len(effective) == 1
        assert np.allclose(effective[0], expected)
