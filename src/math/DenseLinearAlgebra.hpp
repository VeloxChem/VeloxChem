//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#ifndef DenseLinearAlgebra_hpp
#define DenseLinearAlgebra_hpp

#include "DenseMatrix.hpp"
#include "MemBlock.hpp"

namespace denblas { // denblas namespace
    
    /**
     Computes matrix multiplication: A * B.

     @param matrixA the matrix A.
     @param matrixB the matrix B
     @return the matrix A * B.
     */
    CDenseMatrix multAB(const CDenseMatrix& matrixA,
                        const CDenseMatrix& matrixB);
    
    /**
     Computes matrix multiplication: A * B^T.
     
     @param matrixA the matrix A.
     @param matrixB the matrix B
     @return the matrix A * B^T.
     */
    CDenseMatrix multABt(const CDenseMatrix& matrixA,
                         const CDenseMatrix& matrixB);
    
    /**
     Computes matrix multiplication: A^T * B.
     
     @param matrixA the matrix A.
     @param matrixB the matrix B
     @return the matrix A^T * B.
     */
    CDenseMatrix multAtB(const CDenseMatrix& matrixA,
                         const CDenseMatrix& matrixB);
    
    /**
     Computes diagonal matrix and matrix multiplication: diag(M) * A.

     @param diagonal the diagonal matrix.
     @param matrix the square matrix.
     @return the matrix diag(M) * A.
     */
    CDenseMatrix multDiagByA(const CMemBlock<double>& diagonal,
                             const CDenseMatrix&      matrix);
    
    /**
     Computes diagonal matrix and matrix multiplication: diag(M) * A^T.
     
     @param diagonal the diagonal matrix.
     @param matrix the square matrix.
     @return the matrix diag(M) * A^T.
     */
    CDenseMatrix multDiagByAt(const CMemBlock<double>& diagonal,
                              const CDenseMatrix&      matrix);

    /**
     Computes matrix substraction: A - B.
     
     @param matrixA the matrix A.
     @param matrixB the matrix B
     @return the matrix A - B.
     */
    CDenseMatrix subAB(const CDenseMatrix& matrixA,
                       const CDenseMatrix& matrixB);
    
    /**
     Computes matrix addition: A + factor * B.
     
     @param matrixA the matrix A.
     @param matrixB the matrix B
     @param factor the scaling factor of matrix B.
     @return the matrix A +  factor * B.
     */
    CDenseMatrix addAB(const CDenseMatrix& matrixA,
                       const CDenseMatrix& matrixB,
                       const double        factor);
    
    /**
     Computes matrix multiplication: C = C + alpha * A * B^T.
     
     @param matrixC the matrix C.
     @param alpha the scaling factor for A * B^T.
     @param matrixA the matrix A.
     @param matrixB the matrix B
     */
    void multABt(      CDenseMatrix& matrixC,
                 const double        alpha,
                 const CDenseMatrix& matrixA,
                 const CDenseMatrix& matrixB);
    
} // denblas namespace

#endif /* DenseLinearAlgebra_hpp */
