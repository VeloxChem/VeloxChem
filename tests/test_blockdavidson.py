import numpy as np
import pytest

from veloxchem.blockdavidson import BlockDavidsonSolver


@pytest.mark.solvers
class TestBlockDavidsonSolver:

    @staticmethod
    def build_trial_subspace():

        vecs = np.array([
            [1.0, 1.0, 0.0, 0.0],
            [0.0, 1.0, 1.0, 0.0],
            [0.0, 0.0, 1.0, 1.0],
            [1.0, 0.0, 0.0, 1.0],
            [0.0, 1.0, 0.0, 1.0],
            [1.0, 0.0, 1.0, 0.0],
        ])
        qmat, rmat_not_used = np.linalg.qr(vecs)
        return qmat[:, :4]

    def build_diagonal_problem(self):

        amat = np.diag(np.array([1.0, 2.0, 4.0, 5.0, 7.0, 8.0]))
        diag_mat = np.diag(amat)
        trial_mat = self.build_trial_subspace()
        sigma_mat = np.matmul(amat, trial_mat)

        return amat, diag_mat, trial_mat, sigma_mat

    def test_add_iteration_data_stacks_columns(self):

        solver = BlockDavidsonSolver()
        trial_mat = self.build_trial_subspace()
        sigma_mat = np.matmul(np.eye(trial_mat.shape[0]), trial_mat)

        solver.add_iteration_data(sigma_mat[:, :2], trial_mat[:, :2], 2)

        assert solver.neigenpairs == 2
        assert np.allclose(solver.trial_matrices, trial_mat[:, :2])
        assert np.allclose(solver.sigma_matrices, sigma_mat[:, :2])

        solver.add_iteration_data(sigma_mat[:, 2:], trial_mat[:, 2:], 2)

        assert solver.reduced_space_size() == 4
        assert np.allclose(solver.trial_matrices, trial_mat)
        assert np.allclose(solver.sigma_matrices, sigma_mat)

    def test_default_subspace_thresholds_and_no_collapse_below_limit(self):

        solver = BlockDavidsonSolver()

        assert solver._get_max_subspace_dim() is None

        amat_not_used, diag_mat_not_used, trial_mat, sigma_mat = (
            self.build_diagonal_problem())
        solver.add_iteration_data(sigma_mat[:, :2], trial_mat[:, :2], 2)

        assert solver._get_max_subspace_dim() == 40
        assert solver._get_collapse_nvec() == 4
        assert not solver.should_collapse()

        solver.max_subspace_dim = 2
        assert not solver.should_collapse()

        solver.add_iteration_data(sigma_mat[:, 2:3], trial_mat[:, 2:3], 2)
        assert solver.should_collapse()

    def test_residual_and_convergence_accessors(self):

        amat_not_used, diag_mat_not_used, trial_mat, sigma_mat = (
            self.build_diagonal_problem())
        solver = BlockDavidsonSolver(max_subspace_dim=10, collapse_nvec=2)
        solver.add_iteration_data(sigma_mat[:, :2], trial_mat[:, :2], 2)
        solver.add_iteration_data(sigma_mat[:, 2:], trial_mat[:, 2:], 2)

        solver.comp_residual_vectors()

        eigs, norms = solver.get_eigenvalues()
        rmax, rmin = solver.max_min_residual_norms()

        assert eigs.shape == (2,)
        assert norms.shape == (2,)
        assert rmax == pytest.approx(np.amax(norms))
        assert rmin == pytest.approx(np.amin(norms))
        assert solver.check_convergence(rmax + 1.0e-12)
        assert not solver.check_convergence(rmin - 1.0e-12)

    def test_comp_trial_vectors_and_project_trial_vectors(self):

        amat = np.diag(np.array([1.0, 2.0, 3.0, 6.0]))
        diag_mat = np.diag(amat)
        vecs = np.array([
            [1.0, 1.0],
            [1.0, 0.0],
            [0.0, 1.0],
            [0.0, 1.0],
        ])
        trial_mat, rmat_not_used = np.linalg.qr(vecs)
        sigma_mat = np.matmul(amat, trial_mat)

        solver = BlockDavidsonSolver(max_subspace_dim=4, collapse_nvec=2)
        solver.add_iteration_data(sigma_mat, trial_mat, 2)
        solver.comp_residual_vectors()

        trial_updates = solver.comp_trial_vectors(diag_mat)

        assert np.all(np.isfinite(trial_updates))
        assert np.allclose(np.linalg.norm(trial_updates, axis=0),
                           np.ones(2),
                           atol=1.0e-10)

        trial_candidates = np.array([
            [1.0, 0.0],
            [0.0, 1.0],
            [1.0, 1.0],
            [1.0, -1.0],
        ])
        projected_updates = solver.project_trial_vectors(trial_candidates)

        assert projected_updates.shape[0] == 4
        assert projected_updates.shape[1] == 1
        # test orthogonality
        assert np.allclose(np.matmul(solver.trial_matrices.T,
                                     projected_updates),
                           np.zeros((2, projected_updates.shape[1])),
                           atol=1.0e-10)
        # test normalization
        assert np.allclose(np.linalg.norm(projected_updates, axis=0),
                           np.ones(projected_updates.shape[1]),
                           atol=1.0e-10)

    def test_collapse_subspace(self):

        amat, diag_mat_not_used, trial_mat, sigma_mat = (
            self.build_diagonal_problem())

        solver = BlockDavidsonSolver(max_subspace_dim=3, collapse_nvec=2)
        solver.add_iteration_data(sigma_mat[:, :2], trial_mat[:, :2], 2)
        solver.add_iteration_data(sigma_mat[:, 2:], trial_mat[:, 2:], 2)
        solver.comp_residual_vectors()

        expected_eigs = solver.residual_eigs.copy()
        expected_ritz = solver.ritz_vectors.copy()

        solver.collapse_subspace()

        assert solver.collapsed_subspace
        assert solver.collapsed_from_dim == 4
        assert solver.collapsed_to_dim == 2
        assert solver.reduced_space_size() == 2
        assert np.allclose(solver.residual_eigs, expected_eigs)
        assert np.allclose(solver.trial_matrices, expected_ritz)
        assert np.allclose(solver.sigma_matrices,
                           np.matmul(amat, solver.trial_matrices))

    def test_compute_triggers_size_based_collapse(self):

        amat_not_used, diag_mat, trial_mat, sigma_mat = (
            self.build_diagonal_problem())

        solver = BlockDavidsonSolver(max_subspace_dim=3, collapse_nvec=2)
        solver.add_iteration_data(sigma_mat[:, :2], trial_mat[:, :2], 2)
        solver.add_iteration_data(sigma_mat[:, 2:], trial_mat[:, 2:], 2)

        new_trials = solver.compute(diag_mat)

        assert solver.collapsed_subspace
        assert solver.collapsed_from_dim == 4
        assert solver.collapsed_to_dim == 2
        assert solver.reduced_space_size() == 2
        assert new_trials.shape[0] == trial_mat.shape[0]
        assert np.all(np.isfinite(new_trials))
