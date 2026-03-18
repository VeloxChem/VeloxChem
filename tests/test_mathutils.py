import math
import numpy as np
import pytest
from mpi4py import MPI

from veloxchem.mathutils import safe_arccos, safe_arcsin
from veloxchem.mathutils import safe_eigh, safe_solve


class TestMathUtils:

    def test_safe_arccos_handles_regular_and_clipped_inputs(self):

        assert safe_arccos(0.5) == pytest.approx(math.acos(0.5))
        assert safe_arccos(1.0 + 5.0e-13) == pytest.approx(0.0)
        assert safe_arccos(-1.0 - 5.0e-13) == pytest.approx(math.pi)

    @pytest.mark.skipif(MPI.COMM_WORLD.Get_size() > 1,
                        reason="pytest.raises only valid in serial")
    def test_safe_arccos_rejects_invalid_input(self):

        with pytest.raises(AssertionError, match="arccos: Invalid cosine value"):
            safe_arccos(1.0 + 1.0e-8)

    def test_safe_arcsin_handles_regular_and_clipped_inputs(self):

        assert safe_arcsin(0.5) == pytest.approx(math.asin(0.5))
        assert safe_arcsin(1.0 + 5.0e-13) == pytest.approx(math.pi / 2.0)
        assert safe_arcsin(-1.0 - 5.0e-13) == pytest.approx(-math.pi / 2.0)

    @pytest.mark.skipif(MPI.COMM_WORLD.Get_size() > 1,
                        reason="pytest.raises only valid in serial")
    def test_safe_arcsin_rejects_invalid_input(self):

        with pytest.raises(AssertionError, match="arcsin: Invalid sine value"):
            safe_arcsin(-1.0 - 1.0e-8)

    def test_safe_solve_matches_numpy_for_nonsingular_system(self):

        mat = np.array([
            [3.0, 1.0],
            [1.0, 2.0],
        ])
        rhs = np.array([7.0, 5.0])

        solution = safe_solve(mat, rhs)

        np.testing.assert_allclose(solution, np.linalg.solve(mat, rhs))

    def test_safe_solve_falls_back_to_pseudoinverse_for_singular_system(self):

        mat = np.array([
            [1.0, 2.0],
            [2.0, 4.0],
        ])
        rhs = np.array([3.0, 6.0])

        solution = safe_solve(mat, rhs)

        np.testing.assert_allclose(solution, np.matmul(np.linalg.pinv(mat), rhs))
        np.testing.assert_allclose(np.matmul(mat, solution), rhs)

    def test_safe_eigh_screens_linear_dependency_and_builds_pseudoinverse(self):

        basis_vectors = np.array([
            [1.0, 0.0, 1.0, 2.0],
            [0.0, 1.0, 1.0, 2.0],
        ])
        overlap = np.matmul(basis_vectors.T, basis_vectors)

        eigvals, eigvecs = safe_eigh(overlap, 1.0e-10)

        assert eigvecs.shape == (4, 2)
        np.testing.assert_allclose(eigvals, np.array([1.0, 11.0]))
        np.testing.assert_allclose(np.matmul(eigvecs.T, eigvecs),
                                   np.eye(2),
                                   atol=1.0e-12)

        overlap_pinv = np.matmul(
            np.matmul(eigvecs, np.diag(1.0 / eigvals)), eigvecs.T)
        ref_pinv = np.linalg.pinv(overlap, hermitian=True)

        np.testing.assert_allclose(overlap_pinv, ref_pinv, atol=1.0e-12)
        np.testing.assert_allclose(np.matmul(np.matmul(overlap, overlap_pinv),
                                             overlap),
                                   overlap,
                                   atol=1.0e-12)

    def test_safe_eigh_keeps_square_eigensystem_when_no_screening_is_needed(
            self):

        overlap = np.array([
            [2.0, 0.3, 0.1],
            [0.3, 1.5, 0.4],
            [0.1, 0.4, 1.2],
        ])

        eigvals, eigvecs = safe_eigh(overlap, 1.0e-10)

        assert eigvecs.shape == (3, 3)
        np.testing.assert_array_less(1.0e-10, eigvals)
        np.testing.assert_allclose(np.matmul(eigvecs.T, eigvecs),
                                   np.eye(3),
                                   atol=1.0e-12)

        overlap_inv = np.matmul(
            np.matmul(eigvecs, np.diag(1.0 / eigvals)), eigvecs.T)
        ref_inv = np.linalg.inv(overlap)

        np.testing.assert_allclose(overlap_inv, ref_inv, atol=1.0e-12)
        np.testing.assert_allclose(np.matmul(overlap, overlap_inv),
                                   np.eye(3),
                                   atol=1.0e-12)
