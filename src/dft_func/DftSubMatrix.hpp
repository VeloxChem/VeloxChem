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

#ifndef DftSubMatrix_hpp
#define DftSubMatrix_hpp

#include <vector>

#include "AODensityMatrix.hpp"
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
 Gets sub AO density matrices from AO density matrices.

 @param densityPointers the pointers to AO density matrix.
 @param aoIndices the index mapping from submatrix to full matrix.
 @return the sub AO density matrices.
 */
auto getSubAODensityMatrix(const std::vector<const double*>& densityPointers, const std::vector<int>& aoIndices, const int naos) -> CAODensityMatrix;

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

auto distributeSubMatrixToKohnSham(CAOKohnShamMatrix&      aoKohnShamMatrix,
                                   const CDenseMatrix&     subMatrix_a,
                                   const CDenseMatrix&     subMatrix_b,
                                   const std::vector<int>& aoIndices) -> void;

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

void distributeSubMatrixToDenseMatrix(CDenseMatrix&           matrix,
                                      const CDenseMatrix&     subMatrix,
                                      const std::vector<int>& aoIndices,
                                      const int               naos);

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
