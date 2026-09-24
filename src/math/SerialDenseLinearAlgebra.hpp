//
//                                   VELOXCHEM
//              ----------------------------------------------------
//                          An Electronic Structure Code
//
//  SPDX-License-Identifier: BSD-3-Clause
//
//  Copyright 2018-2025 VeloxChem developers
//
//  Redistribution and use in source and binary forms, with or without modification,
//  are permitted provided that the following conditions are met:
//
//  1. Redistributions of source code must retain the above copyright notice, this
//     list of conditions and the following disclaimer.
//  2. Redistributions in binary form must reproduce the above copyright notice,
//     this list of conditions and the following disclaimer in the documentation
//     and/or other materials provided with the distribution.
//  3. Neither the name of the copyright holder nor the names of its contributors
//     may be used to endorse or promote products derived from this software without
//     specific prior written permission.
//
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
//  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
//  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
//  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
//  FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
//  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
//  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
//  HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
//  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
//  OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#ifndef SerialDenseLinearAlgebra_hpp
#define SerialDenseLinearAlgebra_hpp

#include <cstddef>

#include "DenseMatrix.hpp"
#include "SubMatrix.hpp"

/**
 Serial dense linear algebra.

 All routines in this namespace are serial: each call runs on the calling
 thread only and never starts threads of its own. The implementation uses
 Eigen. The routines are intended for the serial inner loops of code that is
 parallelized at a higher level with OpenMP or MPI, so they must not be
 replaced by, or internally route to, a threaded BLAS.
 */
namespace sdenblas {  // sdenblas namespace

/**
 Computes matrix multiplication: A * B.

 @param matrixA the matrix A.
 @param matrixB the matrix B
 @return the matrix A * B.
 */
auto serialMultAB(const CDenseMatrix& matrixA, const CDenseMatrix& matrixB) -> CDenseMatrix;

/**
 Computes matrix multiplication: A * B^T.

 @param matrixA the matrix A.
 @param matrixB the matrix B
 @return the matrix A * B^T.
 */
auto serialMultABt(const CDenseMatrix& matrixA, const CDenseMatrix& matrixB) -> CDenseMatrix;

/**
 Computes matrix multiplication: A^T * B.

 @param matrixA the matrix A.
 @param matrixB the matrix B
 @return the matrix A^T * B.
 */
auto serialMultAtB(const CDenseMatrix& matrixA, const CDenseMatrix& matrixB) -> CDenseMatrix;

/**
 Computes matrix addition: A + B.

 @param matrixA the matrix A.
 @param matrixB the matrix B
 @return the matrix A + B.
 */
auto serialAddAB(const CDenseMatrix& matrixA, const CDenseMatrix& matrixB, const double factor) -> CDenseMatrix;

/**
 Adds matrix B to matrix A.

 @param matrixA the matrix A.
 @param matrixB the matrix B
 */
auto serialInPlaceAddAB(CDenseMatrix& matrixA, const CDenseMatrix& matrixB, const double factor=1.0) -> void;

/**
 Solves Ax=b.

 @param mat the matrix A.
 @param vec the vector b.
 @return the solution vector x.
 */
auto serialSolve(const CDenseMatrix& mat, const std::vector<double>& vec) -> std::vector<double>;

/**
 Computes matrix multiplication: C = alpha * A * B + beta * C, for row major
 matrices with explicit leading dimensions.

 @param nrows the number of rows of A and C.
 @param ncols the number of columns of B and C.
 @param nsums the number of columns of A and rows of B.
 @param alpha the factor of the product.
 @param matrixA the values of A, as a row major array with leading dimension lda.
 @param lda the leading dimension of A.
 @param matrixB the values of B, as a row major array with leading dimension ldb.
 @param ldb the leading dimension of B.
 @param beta the factor of C.
 @param matrixC the values of C, as a row major array with leading dimension ldc.
 @param ldc the leading dimension of C.
 */
auto serialMultAB(const size_t  nrows,
                  const size_t  ncols,
                  const size_t  nsums,
                  const double  alpha,
                  const double *matrixA,
                  const size_t  lda,
                  const double *matrixB,
                  const size_t  ldb,
                  const double  beta,
                  double       *matrixC,
                  const size_t  ldc) -> void;

/**
 Computes matrix multiplication: C = alpha * A * B^T + beta * C, for row major
 matrices with explicit leading dimensions.

 @param nrows the number of rows of A and C.
 @param ncols the number of columns of B and C.
 @param nsums the number of columns of A and B.
 @param alpha the factor of the product.
 @param matrixA the values of A, as a row major array with leading dimension lda.
 @param lda the leading dimension of A.
 @param matrixB the values of B, as a row major array with leading dimension ldb.
 @param ldb the leading dimension of B.
 @param beta the factor of C.
 @param matrixC the values of C, as a row major array with leading dimension ldc.
 @param ldc the leading dimension of C.
 */
auto serialMultABt(const size_t  nrows,
                   const size_t  ncols,
                   const size_t  nsums,
                   const double  alpha,
                   const double *matrixA,
                   const size_t  lda,
                   const double *matrixB,
                   const size_t  ldb,
                   const double  beta,
                   double       *matrixC,
                   const size_t  ldc) -> void;

/**
 Adds a symmetric rank k update: C += alpha * A * A^T, into the lower triangle of
 the row major matrix C.

 @param n the number of rows of A and of C.
 @param k the number of columns of A.
 @param alpha the factor of the update.
 @param matrixA the values of A, as a row major array with leading dimension lda.
 @param lda the leading dimension of A.
 @param matrixC the values of C, as a row major array with leading dimension ldc.
 @param ldc the leading dimension of C.
 */
auto serialRankUpdate(const size_t  n,
                      const size_t  k,
                      const double  alpha,
                      const double *matrixA,
                      const size_t  lda,
                      double       *matrixC,
                      const size_t  ldc) -> void;

/**
 Solves a triangular system with many right hand sides: factor * X = values, or
 factor^T * X = values when the factor is transposed, for a row major lower
 triangular factor and row major right hand sides. The solution is written over
 the right hand sides in place.

 @param nrows the number of rows and columns of the factor, and the rows of values.
 @param ncols the number of columns of values.
 @param factor the values of the lower triangular factor, as a row major array
 with leading dimension ldf.
 @param ldf the leading dimension of the factor.
 @param values the values of the right hand sides, as a row major array with
 leading dimension ldv, overwritten by the solution.
 @param ldv the leading dimension of the values.
 @param transposed whether to solve against the transpose of the factor.
 */
auto serialSolveTriangular(const size_t  nrows,
                           const size_t  ncols,
                           const double *factor,
                           const size_t  ldf,
                           double       *values,
                           const size_t  ldv,
                           const bool    transposed) -> void;

/**
 Computes matrix multiplication: C += A * B.

 @param matrixC the matrix C.
 @param matrixA the matrix A.
 @param matrixB the matrix B
 */
auto serialMultAB(CSubMatrix& matrixC, const CSubMatrix& matrixA, const CSubMatrix& matrixB) -> void;

/**
 Computes matrix multiplication: C += A^t * B.

 @param matrixC the matrix C.
 @param matrixA the matrix A.
 @param matrixB the matrix B
 */
auto serialMultAtB(CSubMatrix& matrixC, const CSubMatrix& matrixA, const CSubMatrix& matrixB) -> void;

}  // namespace sdenblas

#endif /* SerialDenseLinearAlgebra_hpp */
