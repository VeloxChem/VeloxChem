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

#include "AODensityMatrix.hpp"
#include "AOKohnShamMatrix.hpp"
#include "DenseMatrix.hpp"

namespace dftsubmat {  // dftsubmat namespace

/**
 Gets sub AO density matrices from AO density matrices.

 @param densityMatrix the AO density matrix.
 @param densityIndex the index of density matrix.
 @param densitySpin the spin of density matrix.
 @param aoIndices the index mapping from submatrix to full matrix.
 @return the sub AO density matrices.
 */
auto getSubDensityMatrix(const CAODensityMatrix& densityMatrix,
                         const int               densityIndex,
                         const std::string&      densitySpin,
                         const std::vector<int>& aoIndices) -> CDenseMatrix;

/**
 Gets sub AO density matrices from AO density matrices.

 @param densityPointer the pointer to AO density matrix.
 @param aoIndices the index mapping from submatrix to full matrix.
 @param naos the number of AOs.
 @return the sub AO density matrices.
 */
auto getSubDensityMatrix(const double* densityPointer, const std::vector<int>& aoIndices, const int naos) -> CDenseMatrix;

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

}  // namespace dftsubmat

#endif /* DftSubMatrix_hpp */
