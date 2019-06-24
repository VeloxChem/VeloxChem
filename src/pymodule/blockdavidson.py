import numpy as np


class BlockDavidsonSolver:
    """Implements block Davidson solver.

    Implements block Davidson solver for symmetric eigenvalues/eigenvectors
    problems i.e. A * X = e X.
    Reference: SIAM J. Sci. Comput. 15 (1994), 62-76.
    """

    def __init__(self):
        """Initializes block Davidson solver.

        Initializes block Davidson sover to default setup.

        Attributes
        ----------
        neigenpairs
            The number of eigenpairs of A matrix determined by solver.
        sigma_matrices
            The sigma vectors A * X in matrix format {A * X_0, A * X_1,...}.
        trial_matrices
            The trial vectors in matrix format {X_0, X_1,...}.
        """

        # excited states information
        self.neigenpairs = None

        # solver storage setup
        self.sigma_matrices = None
        self.trial_matrices = None

        # residual data
        self.residual_matrices = None
        self.residual_norms = None
        self.residual_eigs = None

        # Ritz data
        self.ritz_vectors = None

    def add_iteration_data(self, sig_mat, trial_mat, iteration):
        """Adds sigma and trial vectors to sigma and trial matrices.

        Add sigma and trial vectors to sigma and trial matrices for specific
        iteration.

        Parameters
        ----------
        sig_mat
            The sigma vectors {A * X_i} i = 0..l.
        trial_mat
            The trial vectors {X_i} i = 0..l.
        iteration
            The index of block Davidson iteration.
        """

        if iteration == 0:
            self.neigenpairs = sig_mat.shape[1]
            self.sigma_matrices = sig_mat
            self.trial_matrices = trial_mat
        else:
            self.sigma_matrices = np.hstack((self.sigma_matrices, sig_mat))
            self.trial_matrices = np.hstack((self.trial_matrices, trial_mat))

    def compute(self, diag_mat):
        """Computes new set of trial vectors by solving eigenvalues problem in
        reduced space.

        Computes new set of trial vectors by solving eigenvalues problem by
        diagonalizing interaction matrix in reduced space.

        Parameters
        ----------
        diag_mat
            The approximate diagonal of symmetric A matrix as column matrix.

        Returns
        -------
        numpy.ndarray
            The new set of trial vectors.
        """

        self.comp_residual_vectors()

        tvecs = self.comp_trial_vectors(diag_mat)

        self.project_trial_vectors(tvecs)
        self.orthogonalize_gram_schmidt(tvecs)

        return tvecs

    def check_convergence(self, conv_thresh):
        """Checks if residual norm is bellow given convergence threshold for all
            eigenpairs.

        Checks if residual norm is bellow given convergence threshold for all
        requested eigenpairs of symmetric A matrix.

        Parameters
        ----------
        conv_thresh
            The approximate diagonal of symmetric A matrix as column matrix.

        Returns
        -------
        bool
            The true if residual norms are converged for all eigenpairs,
            false otherwise.
        """

        for rval in self.residual_norms:
            if rval > conv_thresh:
                return False

        return True

    def reduced_space_size(self):
        """Gets number of trial vectors in residual space.

        Gets number of trial vectors in residual space.

        Returns
        -------
        int
            The number of trial vectors in residual space.
        """

        return self.trial_matrices.shape[1]

    def max_min_residual_norms(self):
        """Determines maximum and minumum residual norms.

        Determines maximum and minumum residual norms within set of requested
        eigenpairs.

        Returns
        -------
        tuple
            tuple (max. residual norm, min. residual norm).
        """

        return (np.amax(self.residual_norms), np.amin(self.residual_norms))

    def get_eigenvalues(self):
        """Gets smallest eigenvalues of symmetric A matrix and associated
        residual norms.

        Gets smallest eigenvalues of symmetric A matrix and associated
        residual norms.

        Returns
        -------
        tuple
            tuple (eigenvalues, residual norms).
        """

        return (self.residual_eigs, self.residual_norms)

    def comp_residual_vectors(self):
        """Computes residual vectors and their norms.

        Computes residual vectors and their norms by diagonalizing interaction
        matrix in reduced space.
        """

        rlmat = np.matmul(self.trial_matrices.transpose(), self.sigma_matrices)

        reigs, rvecs = np.linalg.eigh(rlmat)
        yvecs = rvecs[:, :self.neigenpairs]
        self.residual_eigs = reigs[:self.neigenpairs]

        self.ritz_vectors = np.matmul(self.trial_matrices, yvecs)

        self.residual_matrices = self.ritz_vectors.copy()
        for i in range(self.neigenpairs):
            self.residual_matrices[:, i] *= self.residual_eigs[i]
        self.residual_matrices -= np.matmul(self.sigma_matrices, yvecs)

        self.residual_norms = np.linalg.norm(self.residual_matrices, axis=0)

    def comp_trial_vectors(self, diag_mat):
        """Computes new trial vectors.

        Computes new trial vectors by applying Davidson preconditioner to
        residual vectors.

        Parameters
        ----------
        diag_mat
            The approximate diagonal of symmetric A matrix as column matrix.

        Returns
        -------
        numpy.ndarray
            The trial vectors.
        """

        tvecs = self.residual_matrices.copy()
        for i in range(self.neigenpairs):
            pmat = np.full(diag_mat.shape, self.residual_eigs[i]) - diag_mat
            pmat[:, 0] = 1.0 / pmat[:, 0]
            tvecs[:, i] *= pmat[:, 0]

        self.norm_vectors(tvecs)

        return tvecs

    def project_trial_vectors(self, tvecs):
        """Projects out trial vector components already present in reduced
        space expansion.

        Projects out trial vector components already present in reduced space
        expansion and renormalizes resulting trial vectors.

        Parameters
        ----------
        tvecs
            The trial vectors.
        """

        for i in range(tvecs.shape[1]):
            uvec = np.zeros((tvecs.shape[0], 1))
            for j in range(self.trial_matrices.shape[1]):
                f = np.dot(self.trial_matrices[:, j], tvecs[:, i])
                uvec[:, 0] += f * self.trial_matrices[:, j]
            tvecs[:, i] -= uvec[:, 0]

        self.norm_vectors(tvecs)

    def orthogonalize_gram_schmidt(self, tvecs):
        """Applies modified Gram Schmidt orthogonalization to trial vectors.

        Applies modified Gram Schmidt orthogonalization to trial vectors.

        Parameters
        ----------
        tvecs
            The trial vectors.
        """

        if tvecs.shape[1] > 0:

            f = 1.0 / np.linalg.norm(tvecs[:, 0])
            tvecs[:, 0] *= f

            for i in range(1, tvecs.shape[1]):
                for j in range(i):
                    f = np.dot(tvecs[:, i], tvecs[:, j]) / np.dot(
                        tvecs[:, j], tvecs[:, j])
                    tvecs[:, i] -= f * tvecs[:, j]
                f = 1.0 / np.linalg.norm(tvecs[:, i])
                tvecs[:, i] *= f

    def norm_vectors(self, tvecs):
        """Normalizes to trial vectors.

        Normalizes set of tria vectors by inner product.

        Parameters
        ----------
        tvecs
            The trial vectors.
        """

        fnorms = np.linalg.norm(tvecs, axis=0)

        for i in range(tvecs.shape[1]):
            tvecs[:, i] *= 1.0 / fnorms[i]
