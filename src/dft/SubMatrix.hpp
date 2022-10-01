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

#ifndef SubMatrix_hpp
#define SubMatrix_hpp

#include <vector>

#include "AOFockMatrix.hpp"
#include "AOKohnShamMatrix.hpp"
#include "DenseMatrix.hpp"

namespace submat {  // submat namespace

/**
 Gets submatrix from AO density matrix.

 @param densityMatrix the AO density matrix.
 @param densityIndex the density index.
 @param aoIndices the index mapping from submatrix to full matrix.
 @param aoCount the number of indices in submatrix.
 @param nAOs the number of indices in full matrix.
 @return the submatrix.
 */
CDenseMatrix getSubDensityMatrix(const CAODensityMatrix&     densityMatrix,
                                 const int32_t               densityIndex,
                                 const std::vector<int32_t>& aoIndices,
                                 const int32_t               aoCount,
                                 const int32_t               nAOs);

/**
 Gets sub AO density matrices from AO density matrices.

 @param densityMatrix the AO density matrix.
 @param aoIndices the index mapping from submatrix to full matrix.
 @param aoCount the number of indices in submatrix.
 @return the sub AO density matrices.
 */
CAODensityMatrix getSubDensityMatrix(const CAODensityMatrix&     densityMatrix,
                                     const std::vector<int32_t>& aoIndices,
                                     const int32_t               aoCount);

/**
 Distributes partial Vxc matrix to full AO Kohn-Sham matrix.

 @param matVxc the AO Kohn-Sham matrix.
 @param partialMatVxc the partial Vxc matrix.
 @param aoIndices the index mapping from partial matrix to full matrix.
 @param aoCount the number of indices in partial matrix.
 @param nAOs the number of indices in full matrix.
 */
void distributeSubMatrixToKohnSham(CAOKohnShamMatrix&          aoKohnShamMatrix,
                                   const CDenseMatrix&         subMatrix,
                                   const std::vector<int32_t>& aoIndices,
                                   const int32_t               aoCount,
                                   const int32_t               nAOs);

/**
 Distributes partial Fxc or Kxc matrix to full AO Fock matrix.

 @param aoFockMatrix the AO Fock matrix.
 @param fockIndex the index of Fock matrix.
 @param partialMatFxc the partial Fxc or Kxc matrix.
 @param aoIndices the index mapping from partial matrix to full matrix.
 @param aoCount the number of indices in partial matrix.
 @param nAOs the number of indices in full matrix.
 */
void distributeSubMatrixToFock(CAOFockMatrix&              aoFockMatrix,
                               const int32_t               fockIndex,
                               const CDenseMatrix&         subMatrix,
                               const std::vector<int32_t>& aoIndices,
                               const int32_t               aoCount,
                               const int32_t               nAOs);

}  // namespace submat

#endif /* SubMatrix_hpp */
