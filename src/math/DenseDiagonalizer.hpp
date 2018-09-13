//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#ifndef DenseDiagonalizer_hpp
#define DenseDiagonalizer_hpp

#include "DenseMatrix.hpp"
#include "MemBlock.hpp"

/**
 Class CDenseDiagonalizer provides methods for diagonalization of dense real
 symmetric matrices and for manipulation with computed eigenvalues and
 eigenvectors.
 
 @author Z. Rinkevicius
 */
class CDenseDiagonalizer
{
    /**
     The state of dense diagonalizer object : true - no errors,
     false - otherwise.
     */
    bool _state;
    
    /**
     The availability of eigenvalues and eigenvectors: true - available,
     false - otherwise.
     */
    bool _isSolved; 
    
    /**
     The temporary dense matrix used by diagonalization routine.
     */
    CDenseMatrix _matrix;
    
    /**
     The eigenvectors of dense matrix.
     */
    CDenseMatrix _eigenVectors;
    
    /**
     The eigenvalues of dense matrix.
     */
    CMemBlock<double> _eigenValues;
    
public:
    
    /**
     Creates an dense diagonalizer object.
     */
    CDenseDiagonalizer();
    
    /**
     Destroys a dense diagonalizer object.
     */
    ~CDenseDiagonalizer();
    
    /**
     Diagonalizes dense matrix and stores eigenvalues/eigenvectors into internal
     data buffers.

     @param matrix the dense matrix.
     */
    void diagonalize(const CDenseMatrix& matrix);
    
    /**
     Gets state of dense diagonalizer object.

     @return true if no errors, false otherwise.
     */
    bool getState() const;
    
    /**
     Gets eigenvectors of dense matrix.

     @return the eigenvectors of matrix.
     */
    CDenseMatrix getEigenVectors() const;
    
    /**
     Gets eigenvalues of dense matrix.

     @return the eigenvalues of matrix.
     */
    CMemBlock<double> getEigenValues() const;
    
    /**
     Computes A^-1/2 matrix using eigenvalues and eigenvectors of A matrix.

     @return the A^-1/2 matrix.
     */
    CDenseMatrix getInvertedSqrtMatrix() const;
    
    /**
     Computes A^-1 matrix using eigenvalues and eigenvectors of A matrix.

     @return the A^-1 matrix.
     */
    CDenseMatrix getInvertedMatrix() const;
};

#endif /* DenseDiagonalizer_hpp */
