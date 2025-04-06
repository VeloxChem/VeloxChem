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

#ifndef DenseLinearAlgebra_hpp
#define DenseLinearAlgebra_hpp

#include "DenseMatrix.hpp"

namespace denblas {  // denblas namespace

/**
 Computes matrix multiplication: A * B.

 @param matrixA the matrix A.
 @param matrixB the matrix B
 @return the matrix A * B.
 */
auto multAB(const CDenseMatrix& matrixA, const CDenseMatrix& matrixB) -> CDenseMatrix;

/**
 Computes matrix multiplication: A * B^T.

 @param matrixA the matrix A.
 @param matrixB the matrix B
 @return the matrix A * B^T.
 */
auto multABt(const CDenseMatrix& matrixA, const CDenseMatrix& matrixB) -> CDenseMatrix;

/**
 Computes matrix multiplication: A^T * B.

 @param matrixA the matrix A.
 @param matrixB the matrix B
 @return the matrix A^T * B.
 */
auto multAtB(const CDenseMatrix& matrixA, const CDenseMatrix& matrixB) -> CDenseMatrix;

/**
 Computes matrix substraction: A - B.

 @param matrixA the matrix A.
 @param matrixB the matrix B
 @return the matrix A - B.
 */
auto subAB(const CDenseMatrix& matrixA, const CDenseMatrix& matrixB) -> CDenseMatrix;

/**
 Computes matrix addition: A + factor * B.

 @param matrixA the matrix A.
 @param matrixB the matrix B
 @param factor the scaling factor of matrix B.
 @return the matrix A +  factor * B.
 */
auto addAB(const CDenseMatrix& matrixA, const CDenseMatrix& matrixB, const double factor) -> CDenseMatrix;

/**
 Computes dot product of two vectors.

 @param vectorA the vector A.
 @param vectorB the vector B.
 @return the dot product of A and B.
 */
auto dot(const std::vector<double>& vectorA, const std::vector<double>& vectorB) -> double;

/**
 Computes trace of matrix.

 @param matrix the matrix.
 @return the trace of matrix.
 */
auto trace(const CDenseMatrix& matrix) -> double;

}  // namespace denblas

#endif /* DenseLinearAlgebra_hpp */
