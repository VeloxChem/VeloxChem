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

#ifndef DftSubMatrix_hpp
#define DftSubMatrix_hpp

#include <vector>

#include "AOKohnShamMatrix.hpp"
#include "Dense4DTensor.hpp"
#include "DenseMatrix.hpp"

namespace dftsubmat {  // dftsubmat namespace

/**
 Gets sub AO density matrices from AO density matrices.

 @param densityPointer the pointer to AO density matrix.
 @param aoIndices the index mapping from submatrix to full matrix.
 @param naos the number of AOs.
 @return the sub AO density matrices.
 */
auto getSubDensityMatrix(const double* densityPointer, const std::vector<int>& aoIndices, const int naos) -> CDenseMatrix;

/**
 Gets sub matrix from an arbitrary matrix by slicing the columns.

 @param denseMatrix the matrix to slice.
 @param aoIndices the index mapping from submatrix to full matrix.
 @param naos the number of indices in full matrix.
 @return the sub matrix.
 */
auto getSubMatrixByColumnSlicing(const CDenseMatrix& denseMatrix, const std::vector<int>& aoIndices, const int naos) -> CDenseMatrix;

/**
 Distributes partial matrix to AO Kohn-Sham matrix.

 @param aoKohnShamMatrix the AO Kohn-Sham matrix.
 @param subMatrix the partial matrix.
 @param aoIndices the index mapping from partial matrix to full matrix.
 */
auto distributeSubMatrixToKohnSham(CAOKohnShamMatrix& aoKohnShamMatrix, const CDenseMatrix& subMatrix, const std::vector<int>& aoIndices) -> void;

/**
 Distributes partial matrices to AO Kohn-Sham matrix.

 @param aoKohnShamMatrix the AO Kohn-Sham matrix.
 @param subMatrices the partial matrices.
 @param aoIndices the index mapping from partial matrix to full matrix.
 */
auto distributeSubMatrixToKohnSham(CAOKohnShamMatrix&               aoKohnShamMatrix,
                                   const std::vector<CDenseMatrix>& subMatrices,
                                   const std::vector<int>&          aoIndices) -> void;

/**
 Distributes partial matrices to AO Fock matrix.

 @param aoFockPointers the pointers to AO Fock matrices.
 @param fockIndex the index of the AO Fock matrix.
 @param subMatrices the partial matrices.
 @param aoIndices the index mapping from partial matrix to full matrix.
 @param naos the number of AOs.
 */
auto distributeSubMatrixToFock(const std::vector<double*>& aoFockPointers,
                               const int                   fockIndex,
                               const CDenseMatrix&         subMatrix,
                               const std::vector<int>&     aoIndices,
                               const int                   naos) -> void;

/**
 Distributes partial Wxc tensor (pair functional) to full Wxc tensor.

 @param fullTensor the full Wxc tensor.
 @param subMatrix the partial Wxc matrix.
 @param aoIndices the index mapping from partial to full.
 @param aoCount the number of indices in partial matrix.
 */
auto distributeSubmatrixTo4DTensor(CDense4DTensor& fullTensor, const CDenseMatrix& subMatrix, const std::vector<int>& aoIndices) -> void;

}  // namespace dftsubmat

#endif /* DftSubMatrix_hpp */
