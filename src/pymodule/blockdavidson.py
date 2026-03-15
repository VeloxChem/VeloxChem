#
#                                   VELOXCHEM
#              ----------------------------------------------------
#                          An Electronic Structure Code
#
#  SPDX-License-Identifier: BSD-3-Clause
#
#  Copyright 2018-2025 VeloxChem developers
#
#  Redistribution and use in source and binary forms, with or without modification,
#  are permitted provided that the following conditions are met:
#
#  1. Redistributions of source code must retain the above copyright notice, this
#     list of conditions and the following disclaimer.
#  2. Redistributions in binary form must reproduce the above copyright notice,
#     this list of conditions and the following disclaimer in the documentation
#     and/or other materials provided with the distribution.
#  3. Neither the name of the copyright holder nor the names of its contributors
#     may be used to endorse or promote products derived from this software without
#     specific prior written permission.
#
#  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
#  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
#  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
#  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
#  FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
#  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
#  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
#  HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
#  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
#  OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

import numpy as np

from .linearsolver import LinearSolver


class BlockDavidsonSolver:
    """
    Implements block Davidson solver for symmetric eigenvalues/eigenvectors
    problems i.e. A * X = e X. Reference: SIAM J. Sci. Comput. 15 (1994),
    62-76.

    Instance variables
        - neigenpairs: The number of eigenpairs of A matrix determined by
          solver.
        - sigma_matrices: The sigma vectors A * X in matrix format {A * X_0, A
          * X_1,...}.
        - trial_matrices: The trial vectors in matrix format {X_0, X_1,...}.
    """

    def __init__(self, max_subspace_dim=None, collapse_nvec=None):
        """
        Initializes block Davidson sover to default setup.
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
        # subspace restart setup
        self.max_subspace_dim = max_subspace_dim
        self.collapse_nvec = collapse_nvec
        self.collapsed_subspace = False
        self.collapsed_from_dim = None
        self.collapsed_to_dim = None

    def add_iteration_data(self, sig_mat, trial_mat, neigenpairs):
        """
        Add sigma and trial vectors to sigma and trial matrices for specific
        iteration.

        :param neigenpairs:
            The number of eigenpairs.
        :param sig_mat:
            The sigma vectors {A * X_i} i = 0..l.
        :param trial_mat:
            The trial vectors {X_i} i = 0..l.
        """

        self.neigenpairs = neigenpairs

        if self.sigma_matrices is None:
            self.sigma_matrices = sig_mat
        else:
            self.sigma_matrices = np.hstack((self.sigma_matrices, sig_mat))

        if self.trial_matrices is None:
            self.trial_matrices = trial_mat
        else:
            self.trial_matrices = np.hstack((self.trial_matrices, trial_mat))

    def compute(self, diag_mat):
        """
        Computes new set of trial vectors by solving eigenvalues problem by
        diagonalizing interaction matrix in reduced space.

        :param diag_mat:
            The approximate diagonal of symmetric A matrix as column matrix.

        :return:
            The new set of trial vectors.
        """

        self.comp_residual_vectors()
        self.collapsed_subspace = False
        self.collapsed_from_dim = None
        self.collapsed_to_dim = None

        if self.should_collapse():
            self.collapse_subspace()

        tvecs = self.comp_trial_vectors(diag_mat)
        tvecs = self.project_trial_vectors(tvecs)

        return tvecs

    def check_convergence(self, conv_thresh):
        """
        Checks if residual norm is bellow given convergence threshold for all
        requested eigenpairs of symmetric A matrix.

        :param conv_thresh:
            The convergence threshold for the solver.

        :return:
            True if residual norms are converged for all eigenpairs, False
            otherwise.
        """

        for rval in self.residual_norms:
            if rval > conv_thresh:
                return False

        return True

    def reduced_space_size(self):
        """
        Gets number of trial vectors in residual space.

        :return:
            The number of trial vectors in residual space.
        """

        return self.trial_matrices.shape[1]

    def should_collapse(self):
        """
        Checks whether the reduced-space basis should be collapsed.

        :return:
            True if the reduced-space dimension exceeds the collapse threshold.
        """

        max_subspace_dim = self._get_max_subspace_dim()

        return (self.trial_matrices is not None and
                max_subspace_dim is not None and
                self.reduced_space_size() > max_subspace_dim)

    def max_min_residual_norms(self):
        """
        Determines maximum and minumum residual norms within set of requested
        eigenpairs.

        :return:
            tuple (max. residual norm, min. residual norm).
        """

        return (np.amax(self.residual_norms), np.amin(self.residual_norms))

    def get_eigenvalues(self):
        """
        Gets smallest eigenvalues of symmetric A matrix and associated
        residual norms.

        :return:
            tuple (eigenvalues, residual norms).
        """

        return (self.residual_eigs, self.residual_norms)

    def comp_residual_vectors(self):
        """
        Computes residual vectors and their norms by diagonalizing interaction
        matrix in reduced space.
        """

        rlmat = np.matmul(self.trial_matrices.T, self.sigma_matrices)

        reigs, rvecs = np.linalg.eigh(rlmat)
        yvecs = rvecs[:, :self.neigenpairs].copy()
        self.residual_eigs = reigs[:self.neigenpairs]
        self.ritz_vectors = np.matmul(self.trial_matrices, yvecs)

        self.residual_matrices = self.ritz_vectors.copy()
        for i in range(self.neigenpairs):
            self.residual_matrices[:, i] *= self.residual_eigs[i]
        self.residual_matrices -= np.matmul(self.sigma_matrices, yvecs)

        self.residual_norms = np.linalg.norm(self.residual_matrices, axis=0)

    def collapse_subspace(self):
        """
        Collapses the reduced space to the leading Ritz vectors.
        """

        keep_dim = self._get_collapse_nvec()
        keep_dim = min(keep_dim, self.reduced_space_size())

        rlmat = np.matmul(self.trial_matrices.T, self.sigma_matrices)
        reigs, rvecs = np.linalg.eigh(rlmat)
        coeffs = rvecs[:, :keep_dim].copy()

        self.collapsed_subspace = True
        self.collapsed_from_dim = self.reduced_space_size()
        self.collapsed_to_dim = keep_dim

        self.trial_matrices = np.matmul(self.trial_matrices, coeffs)
        self.sigma_matrices = np.matmul(self.sigma_matrices, coeffs)

        # Preserve the current Ritz and residual information for the requested
        # eigenpairs in the collapsed basis.
        self.ritz_vectors = self.trial_matrices[:, :self.neigenpairs].copy()
        self.residual_eigs = reigs[:self.neigenpairs]

    def _get_max_subspace_dim(self):
        """
        Gets the maximum allowed reduced-space dimension.
        """

        if self.max_subspace_dim is not None:
            return self.max_subspace_dim

        if self.neigenpairs is None:
            return None

        return 20 * self.neigenpairs

    def _get_collapse_nvec(self):
        """
        Gets the number of vectors retained when collapsing the subspace.
        """

        if self.collapse_nvec is not None:
            return max(self.neigenpairs, self.collapse_nvec)

        return 2 * self.neigenpairs

    def comp_trial_vectors(self, diag_mat):
        """
        Computes new trial vectors by applying Davidson preconditioner to
        residual vectors.

        :param diag_mat:
            The approximate diagonal of symmetric A matrix as column matrix.

        :return:
            The trial vectors.
        """

        tvecs = self.residual_matrices.copy()
        for i in range(self.neigenpairs):
            tvecs[:, i] /= (self.residual_eigs[i] - diag_mat)

        tvecs = LinearSolver.normalize(tvecs)

        return tvecs

    def project_trial_vectors(self, tvecs):
        """
        Projects out trial vector components already present in reduced
        space expansion and renormalizes resulting trial vectors.

        :param tvecs:
            The trial vectors.
        """

        tvecs = tvecs - np.matmul(self.trial_matrices,
                                  np.matmul(self.trial_matrices.T, tvecs))

        tvecs = LinearSolver.remove_linear_dependence(tvecs, 1.0e-6)
        tvecs = LinearSolver.orthogonalize_gram_schmidt(tvecs)
        tvecs = LinearSolver.normalize(tvecs)

        return tvecs
