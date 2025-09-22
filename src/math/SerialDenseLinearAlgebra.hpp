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

#include "DenseMatrix.hpp"
#include "SubMatrix.hpp"

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
 Computes matrix multiplication: C = A * B.

 @param matrixC the matrix C.
 @param matrixA the matrix A.
 @param matrixB the matrix B
 */
auto serialMultAB(CSubMatrix& matrixC, const CSubMatrix& matrixA, const CSubMatrix& matrixB) -> void;

/**
 Computes matrix multiplication: C = A^t * B.

 @param matrixC the matrix C.
 @param matrixA the matrix A.
 @param matrixB the matrix B
 */
auto serialMultAtB(CSubMatrix& matrixC, const CSubMatrix& matrixA, const CSubMatrix& matrixB) -> void;

}  // namespace sdenblas

#endif /* SerialDenseLinearAlgebra_hpp */
