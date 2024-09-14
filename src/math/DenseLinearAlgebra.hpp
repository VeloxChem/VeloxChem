//
//                              VELOXCHEM
//         ----------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2023 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact
//
//  SPDX-License-Identifier: LGPL-3.0-or-later
//
//  This file is part of VeloxChem.
//
//  VeloxChem is free software: you can redistribute it and/or modify it under
//  the terms of the GNU Lesser General Public License as published by the Free
//  Software Foundation, either version 3 of the License, or (at your option)
//  any later version.
//
//  VeloxChem is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
//  License for more details.
//
//  You should have received a copy of the GNU Lesser General Public License
//  along with VeloxChem. If not, see <https://www.gnu.org/licenses/>.

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
 Computes diagonal matrix and matrix multiplication: diag(M) * A.

 @param diagonal the diagonal matrix.
 @param matrix the square matrix.
 @return the matrix diag(M) * A.
 */
auto multDiagByA(const std::vector<double>& diagonal, const CDenseMatrix& matrix) -> CDenseMatrix;

/**
 Computes diagonal matrix and matrix multiplication: diag(M) * A^T.

 @param diagonal the diagonal matrix.
 @param matrix the square matrix.
 @return the matrix diag(M) * A^T.
 */
auto multDiagByAt(const std::vector<double>& diagonal, const CDenseMatrix& matrix) -> CDenseMatrix;

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
