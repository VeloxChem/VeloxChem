import numpy as np

class CTwoDiis:
    """Implements direct inversion of the iterative subspace.
        
    Implements direct inversion of the iterative subspace in C2 form proposed by
    H. Seller.
        
    Attributes
    ----------
    error_vectors
        The list of error vectors.
    """
    
    def __init__(self):
        """Initializes iterative subspace.
            
        Initializes iterative subspace by setting list of error vectors to
        empty list.
        """
    
        self.error_vectors = []
    
    def compute_error_vectors(self, fock_matrices, density_matrices,
                              overlap_matrix, oao_matrix):
        """Computes error vectors.
        
        Computes error vectors for list of AO Fock matrices using (FDS - SDF)
        in orthogonal AO basis.
        
        Parameters
        ----------
        fock_matrices
            The list of AO Fock matrices.
        density_matrices
            The list of AO density matrices.
        overlap_matrix
            The overlap matrix.
        oao_matrix
            The orthogonalization matrix.
        """

        smat = overlap_matrix.to_numpy()
        tmat = oao_matrix.to_numpy()
    
        self.error_vectors.clear()
        
        for fmat, dmat in zip(fock_matrices, density_matrices):
            
            fa = np.matmul(fmat, np.matmul(dmat, smat))
            fb = np.matmul(smat, np.matmul(dmat, fmat))
            
            self.error_vectors.append(np.matmul(tmat.transpose(),
                                                np.matmul(fa - fb, tmat)))

    def compute_weights(self):
        """Computes C2-DIIS weights.
        
        Computes C2-DIIS weights from error vectors using H. Sellers method
        (Int. J. Quantum Chem., vol. 45, pp. 31-41, 1993.)
        
        Returns
        -------
        numpy.ndarray
            The array of C2-DIIS weights with smallest residual error.
        """
    
        bmat = self.comp_bmatrix()
        
        beigs, bvecs = np.linalg.eigh(bmat)
        
        weights = self.pick_weights(self.norm_bvectors(bvecs))
        
        return weights

    def comp_bmatrix(self):
        """Computes B-matrix.
            
        Computes B-matrix of C2-DIIS method using error vectors.
            
        Returns
        -------
        numpy.ndarray
            The B-matrix.
        """
        
        bdim = len(self.error_vectors)
        
        bmat = np.zeros(shape=(bdim, bdim), dtype=float)
        bidx = np.triu_indices(bdim)
    
        for i, j in zip(bidx[0], bidx[1]):

            fij = np.vdot(self.error_vectors[i], self.error_vectors[j])
                
            bmat[i][j] = fij
            bmat[j][i] = fij
        
        return bmat

    def norm_bvectors(self, bvectors):
        """Normalizes B-matrix eigenvectors.
            
        Normalizes B-matrix eigenvectors by rescaling them to 1.0.
        
        Parameters
        ----------
        bvectors
            The array of B-matrix eigenvectors.
            
        Returns
        -------
        numpy.ndarray
            The normalized B-matrix eigenvectors.
        """
    
        sum_vecs = np.sum(bvectors, axis=0)
    
        norm_vecs = []
        
        for i in range(len(sum_vecs)):
            
            if abs(sum_vecs[i]) > 1.0e-6:
                norm_vecs.append(bvectors[:,i] / sum_vecs[i])
    
        return norm_vecs

    def pick_weights(self, weights):
        """Picks normalize B-matrix eigenvector with smallest residual error.
            
        Picks normalize B-matrix eigenvector with smallest residual error by
        computing residual error for all eigenvectors of B_matrix.
            
        Parameters
        ----------
        bvectors
            The array of B-matrix eigenvectors.
            
        Returns
        -------
        numpy.ndarray
            The normalized B-matrix eigenvector.
        """
    
        fmin = 1.0e+8
        wmin = weights[0]
        
        for w in weights:
            
            evec = np.zeros(self.error_vectors[0].shape, dtype=float)
            
            for f, v in zip(w, self.error_vectors):
                evec = evec + f * v

            fact = np.vdot(evec, evec)

            if fmin > fact:
                fmin = fact
                wmin = w
    
        return wmin
    
