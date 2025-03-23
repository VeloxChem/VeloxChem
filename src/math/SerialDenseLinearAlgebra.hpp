//
//                              VELOXCHEM
//         ----------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2024 by VeloxChem developers. All rights reserved.
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

#ifndef SerialDenseLinearAlgebra_hpp
#define SerialDenseLinearAlgebra_hpp

#include "DenseMatrix.hpp"

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

}  // namespace sdenblas

#endif /* SerialDenseLinearAlgebra_hpp */
