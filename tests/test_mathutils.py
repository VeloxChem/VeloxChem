import math
import numpy as np
import pytest
from mpi4py import MPI

from veloxchem.mathutils import safe_arccos, safe_arcsin
from veloxchem.mathutils import (safe_solve, screened_eigh,
                                 symmetric_matrix_function,
                                 solve_in_orthogonal_basis)


class TestMathUtils:

    def test_safe_arccos_handles_regular_and_clipped_inputs(self):

        assert safe_arccos(0.5) == pytest.approx(math.acos(0.5))
        assert safe_arccos(1.0 + 5.0e-13) == pytest.approx(0.0)
        assert safe_arccos(-1.0 - 5.0e-13) == pytest.approx(math.pi)

    @pytest.mark.skipif(MPI.COMM_WORLD.Get_size() > 1,
                        reason="pytest.raises only valid in serial")
    def test_safe_arccos_rejects_invalid_input(self):

        with pytest.raises(AssertionError,
                           match="arccos: Invalid cosine value"):
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

        np.testing.assert_allclose(solution, np.matmul(np.linalg.pinv(mat),
                                                       rhs))
        np.testing.assert_allclose(np.matmul(mat, solution), rhs)

    def test_screened_eigh_screens_linear_dependency_and_builds_pseudoinverse(
            self):

        basis_vectors = np.array([
            [1.0, 0.0, 1.0, 2.0],
            [0.0, 1.0, 1.0, 2.0],
        ])
        overlap = np.matmul(basis_vectors.T, basis_vectors)

        eigvals, eigvecs = screened_eigh(overlap, 1.0e-10)

        assert eigvecs.shape == (4, 2)
        np.testing.assert_allclose(eigvals, np.array([1.0, 11.0]))
        np.testing.assert_allclose(np.matmul(eigvecs.T, eigvecs),
                                   np.eye(2),
                                   atol=1.0e-12)

        overlap_pinv = np.linalg.multi_dot(
            [eigvecs, np.diag(1.0 / eigvals), eigvecs.T])
        ref_pinv = np.linalg.pinv(overlap, hermitian=True)

        np.testing.assert_allclose(overlap_pinv, ref_pinv, atol=1.0e-12)
        np.testing.assert_allclose(np.linalg.multi_dot(
            [overlap, overlap_pinv, overlap]),
                                   overlap,
                                   atol=1.0e-12)

    def test_screened_eigh_keeps_square_eigensystem_when_no_screening_is_needed(
            self):

        overlap = np.array([
            [2.0, 0.3, 0.1],
            [0.3, 1.5, 0.4],
            [0.1, 0.4, 1.2],
        ])

        eigvals, eigvecs = screened_eigh(overlap, 1.0e-10)

        assert eigvecs.shape == (3, 3)
        np.testing.assert_array_less(1.0e-10, eigvals)
        np.testing.assert_allclose(np.matmul(eigvecs.T, eigvecs),
                                   np.eye(3),
                                   atol=1.0e-12)

        overlap_inv = np.linalg.multi_dot(
            [eigvecs, np.diag(1.0 / eigvals), eigvecs.T])
        ref_inv = np.linalg.inv(overlap)

        np.testing.assert_allclose(overlap_inv, ref_inv, atol=1.0e-12)
        np.testing.assert_allclose(np.matmul(overlap, overlap_inv),
                                   np.eye(3),
                                   atol=1.0e-12)

    def test_screened_eigh_supports_descending_order(self):

        mat = np.diag([1.0, 3.0, 2.0])

        eigvals, eigvecs = screened_eigh(mat, thresh=None, descending=True)

        np.testing.assert_allclose(eigvals, np.array([3.0, 2.0, 1.0]))
        np.testing.assert_allclose(np.matmul(eigvecs.T, eigvecs),
                                   np.eye(3),
                                   atol=1.0e-12)
        np.testing.assert_allclose(np.linalg.multi_dot(
            [eigvecs, np.diag(eigvals), eigvecs.T]),
                                   mat,
                                   atol=1.0e-12)

    def test_screened_eigh_allows_disabling_screening(self):

        mat = np.diag([0.0, 2.0, 5.0])

        eigvals, eigvecs = screened_eigh(mat, thresh=None)

        np.testing.assert_allclose(eigvals, np.array([0.0, 2.0, 5.0]))
        assert eigvecs.shape == (3, 3)
        np.testing.assert_allclose(np.matmul(eigvecs.T, eigvecs),
                                   np.eye(3),
                                   atol=1.0e-12)

    def test_symmetric_matrix_function_computes_matrix_square_root(self):

        mat = np.array([
            [2.0, 0.3, 0.1],
            [0.3, 1.5, 0.4],
            [0.1, 0.4, 1.2],
        ])

        sqrt_mat = symmetric_matrix_function(mat, np.sqrt)

        np.testing.assert_allclose(sqrt_mat, sqrt_mat.T, atol=1.0e-12)
        np.testing.assert_allclose(np.matmul(sqrt_mat, sqrt_mat),
                                   mat,
                                   atol=1.0e-12)

    def test_symmetric_matrix_function_builds_screened_pseudoinverse(self):

        basis_vectors = np.array([
            [1.0, 0.0, 1.0, 2.0],
            [0.0, 1.0, 1.0, 2.0],
        ])
        overlap = np.matmul(basis_vectors.T, basis_vectors)

        overlap_pinv = symmetric_matrix_function(overlap,
                                                 lambda x: 1.0 / x,
                                                 thresh=1.0e-10)

        ref_pinv = np.linalg.pinv(overlap, hermitian=True)

        np.testing.assert_allclose(overlap_pinv, ref_pinv, atol=1.0e-12)
        np.testing.assert_allclose(np.linalg.multi_dot(
            [overlap, overlap_pinv, overlap]),
                                   overlap,
                                   atol=1.0e-12)

    def test_solve_in_orthogonal_basis_matches_direct_solution(self):

        overlap = np.array([
            [1.8, 0.2, 0.1],
            [0.2, 1.4, 0.3],
            [0.1, 0.3, 1.2],
        ])
        fock = np.array([
            [0.9, 0.1, 0.2],
            [0.1, 1.4, 0.3],
            [0.2, 0.3, 1.8],
        ])

        s_eigvals, s_eigvecs = np.linalg.eigh(overlap)
        transform = s_eigvecs * (1.0 / np.sqrt(s_eigvals))

        eigvals, coeffs = solve_in_orthogonal_basis(fock, transform)

        ortho_fock = np.linalg.multi_dot([transform.T, fock, transform])
        ref_eigvals, ref_eigvecs = np.linalg.eigh(ortho_fock)
        ref_coeffs = np.matmul(transform, ref_eigvecs)

        np.testing.assert_allclose(eigvals, ref_eigvals, atol=1.0e-12)
        np.testing.assert_allclose(np.linalg.multi_dot(
            [coeffs.T, overlap, coeffs]),
                                   np.eye(3),
                                   atol=1.0e-12)
        np.testing.assert_allclose(
            np.linalg.multi_dot([coeffs, np.diag(eigvals), coeffs.T]),
            np.linalg.multi_dot(
                [ref_coeffs, np.diag(ref_eigvals), ref_coeffs.T]),
            atol=1.0e-12)

    def test_solve_in_orthogonal_basis_supports_descending_order(self):

        mat = np.diag([1.0, 3.0, 2.0])
        transform = np.eye(3)

        eigvals, coeffs = solve_in_orthogonal_basis(mat,
                                                    transform,
                                                    descending=True)

        np.testing.assert_allclose(eigvals, np.array([3.0, 2.0, 1.0]))
        np.testing.assert_allclose(np.matmul(coeffs.T, coeffs),
                                   np.eye(3),
                                   atol=1.0e-12)
        np.testing.assert_allclose(np.linalg.multi_dot(
            [coeffs, np.diag(eigvals), coeffs.T]),
                                   mat,
                                   atol=1.0e-12)
