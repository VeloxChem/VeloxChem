//
//                           VELOXCHEM 1.0-RC2
//         ----------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2021 by VeloxChem developers. All rights reserved.
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

#ifndef DftSubMatrix_hpp
#define DftSubMatrix_hpp

#include <vector>

#include "AOKohnShamMatrix.hpp"
#include "DenseMatrix.hpp"

namespace dftsubmat {  // dftsubmat namespace

/**
 Gets sub AO density matrices from AO density matrices.

 @param densityMatrix the AO density matrix.
 @param aoIndices the index mapping from submatrix to full matrix.
 @return the sub AO density matrices.
 */
auto getSubDensityMatrix(const CDenseMatrix& densityMatrix, const std::vector<int64_t>& aoIndices) -> CDenseMatrix;

/**
 Distributes partial matrix to full matrix.

 @param matrix the full matrix.
 @param subMatrix the partial matrix.
 @param aoIndices the index mapping from partial matrix to full matrix.
 */
auto distributeSubMatrixToDenseMatrix(CDenseMatrix& matrix, const CDenseMatrix& subMatrix, const std::vector<int64_t>& aoIndices) -> void;

/**
 Distributes partial matrix to AO Kohn-Sham matrix.

 @param matrix the AO Kohn-Sham matrix.
 @param subMatrix the partial matrix.
 @param aoIndices the index mapping from partial matrix to full matrix.
 */
auto distributeSubMatrixToKohnSham(CAOKohnShamMatrix& aoKohnShamMatrix, const CDenseMatrix& subMatrix, const std::vector<int64_t>& aoIndices) -> void;

}  // namespace dftsubmat

#endif /* DftSubMatrix_hpp */
