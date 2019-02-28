import numpy as np

class BlockDavidsonSolver:
    """Implements block Davidson solver.
        
    Implements block Davidson solver for symmetric eigenvalues/eigenvectors
    problems i.e. A * X = e X.
    """

    def __init__(self, nstates, nmatrices):
        """Initializes block Davidson solver.
            
        Initializes block Davidson sover to default setup.
        """
    
        # excited states information
        self.nstates = nstates
        
        # solver storage setup
        self.nmatrices = nmatrices
        self.sigma_matrices = None
        self.trial_matrices = None
        
        # residue data
        self.resid_matrices = None
        self.resid_norms = None
        self.resid_eigs = None
    
        # Ritz data
        self.ritz_vectors = None
    
    
    def add_iteration_data(self, sig_mat, trial_mat, iter):
    
        if self.nstates > 0:
            
            if iter == 0:
                self.sigma_matrices = sig_mat
                self.trial_matrices = trial_mat
            else:
                self.sigma_matrices = np.hstack((self.sigma_matrices, sig_mat))
                self.trial_matrices = np.hstack((self.trial_matrices, trial_mat))

    def compute(self, diag_mat):
    
        tvecs = []
    
        if self.nstates > 0:
            
            self.comp_residues()
            
            tvecs = self.comp_trial_vectors(diag_mat)
            self.ortho_trial_vectors(tvecs)
            self.ortho_gram_schmidt(tvecs)
            
            return tvecs
        
        return tvecs
    
    def check_convergence(self, conv_thresh):
    
        for rval in self.resid_norms:
            if rval > conv_thresh:
                return False
    
        return True
    
    def red_space_size(self):
        return self.trial_matrices.shape[1]
    
    def max_min_residues(self):
        return (np.amax(self.resid_norms),np.amin(self.resid_norms))
    
    def get_eigenvalues(self):
        return (self.resid_eigs, self.resid_norms)
    
    def comp_residues(self):

        if self.nstates > 0:
            
            rlmat = np.matmul(self.trial_matrices.transpose(),
                              self.sigma_matrices)
            
            reigs, rvecs = np.linalg.eigh(rlmat)
            yvecs = rvecs[:, :self.nstates]
            self.resid_eigs = reigs[:self.nstates]
            
            self.ritz_vectors = np.matmul(self.trial_matrices, yvecs)
            
            self.resid_matrices = self.ritz_vectors.copy()
            for i in range(self.nstates):
                self.resid_matrices[:,i] *= self.resid_eigs[i]
            self.resid_matrices -= np.matmul(self.sigma_matrices, yvecs)
            self.resid_norms = np.linalg.norm(self.resid_matrices, axis=0)

    def comp_trial_vectors(self, diag_mat):

        if self.nstates > 0:
           
            tvecs = self.resid_matrices.copy()
            for i in range(self.nstates):
                pmat = np.full(diag_mat.shape, self.resid_eigs[i]) - diag_mat
                pmat[:,0] = 1.0/pmat[:,0]
                tvecs[:, i] *= pmat[:,0]
        
            fnorms = np.linalg.norm(tvecs, axis=0)
            for i in range(self.nstates):
                tvecs[:,i] *= 1.0/fnorms[i]
    
            return tvecs;

        return None

    def ortho_trial_vectors(self, tvecs):

        if self.nstates > 0:

            for i in range(self.nstates):
                uvec = np.zeros((tvecs.shape[0],1))
                for j in range(self.trial_matrices.shape[1]):
                    f = np.dot(self.trial_matrices[:,j],tvecs[:,i])
                    uvec[:,0] += f * self.trial_matrices[:,j]
                tvecs[:,i] -= uvec[:,0]

            fnorms = np.linalg.norm(tvecs, axis=0)
            for i in range(self.nstates):
                tvecs[:,i] *= 1.0/fnorms[i]

    def ortho_gram_schmidt(self, tvecs):
        
        if tvecs.shape[1] > 0:
            
            f = 1.0 / np.linalg.norm(tvecs[:,0])
            tvecs[:,0] *= f
            
            for i in range(1,tvecs.shape[1]):
                for j in range(i):
                    f = np.dot(tvecs[:,i],tvecs[:,j])/np.dot(tvecs[:,j],tvecs[:,j])
                    tvecs[:,i] -= f * tvecs[:,j]
                f = 1.0 / np.linalg.norm(tvecs[:,i])
                tvecs[:,i] *= f








            








