//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#ifndef DenseLinearAlgebra_hpp
#define DenseLinearAlgebra_hpp

#include "DenseMatrix.hpp"
#include "MemBlock.hpp"

namespace denblas {  // denblas namespace

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
 Computes matrix multiplication: C = beta C + alpha * A * B^T.

 @param matrixC the matrix C.
 @param alpha the scaling factor for A * B^T.
 @param beta the scaling factor for matrix C.
 @param matrixA the matrix A.
 @param matrixB the matrix B.
 */
void multABt(      CDenseMatrix& matrixC,
             const double        alpha,
             const double        beta,
             const CDenseMatrix& matrixA,
             const CDenseMatrix& matrixB);
    
/**
 Computes matrix multiplication: C = beta C + alpha * A * B^T.

 @param matrixC the pointer to matrix C.
 @param alpha the scaling factor for A * B^T.
 @param beta the scaling factor for matrix C.
 @param matrixA the matrix A.
 @param matrixB the matrix B.
 */
void multABt(      double*       matrixC,
             const double        alpha,
             const double        beta,
             const CDenseMatrix& matrixA,
             const CDenseMatrix& matrixB);
    
/**
Computes matrix multiplication: C = beta C + alpha * A^T * B.
     
@param matrixC the pointer to matrix C.
@param alpha the scaling factor for A^T * B.
@param beta the scaling factor for matrix C.
@param matrixA the matrix A.
@param matrixB the matrix B.
*/
void multAtB(      double*       matrixC,
             const double        alpha,
             const double        beta,
             const CDenseMatrix& matrixA,
             const CDenseMatrix& matrixB);
    
/**
 Computes matrix multiplication: C = beta C + alpha * A * B.

 @param matrixC the pointer to matrix C.
 @param alpha the scaling factor for A * B.
 @param beta the scaling factor for matrix C.
 @param matrixA the matrix A.
 @param matrixB the matrix B
 */
void multAB(      double*       matrixC,
            const double        alpha,
            const double        beta,
            const CDenseMatrix& matrixA,
            const CDenseMatrix& matrixB);
    
/**
 Computes dot product of two vectors.

 @param vectorA the vector A.
 @param vectorB the vector B.
 @return the dot product of A and B.
 */
double dot(const CMemBlock<double>& vectorA,
           const CMemBlock<double>& vectorB);

/**
 Computes dot product of vector and column vector.

 @param vectorA the vector A.
 @param matrixB the column matrix B.
 @return the dot product of vector A and column matrix B.
 */
double dot(const CMemBlock<double>& vectorA,
           const CDenseMatrix&      matrixB);

/**
 Computes trace of matrix.

 @param matrix the matrix.
 @return the trace of matrix.
 */
double trace(const CDenseMatrix& matrix);

/**
 Computes trace of matrix multiplication A * B.

 @param matrixA the matrix A.
 @param matrixB the matrix B.
 @return the trace of A * B.
 */
double trace(const CDenseMatrix& matrixA,
             const CDenseMatrix& matrixB);

}  // namespace denblas

#endif /* DenseLinearAlgebra_hpp */
