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

#ifndef SparseMatrixFunc_hpp
#define SparseMatrixFunc_hpp

#include <vector>

#include "DenseMatrix.hpp"
#include "MolecularBasis.hpp"
#include "ScreenedBasisFunctionPair.hpp"
#include "SparseMatrix.hpp"

namespace sparsematfunc {  // sparsematfunc namespace

/// @brief Converts a block sparse matrix to a full sized dense matrix in
/// VeloxChem atomic orbital ordering. Each block of the sparse matrix is
/// scattered into the dense matrix using the screened basis function pair that
/// defines it: the bra and ket angular momenta and local basis function
/// indices, and, per surviving atom pair, the bra and ket atoms. The atomic
/// orbital index of a contracted function of angular momentum l at molecule
/// order position p and spherical component c is
///   dimensions_of_basis(l) + c * number_of_basis_functions(l) + p.
/// Block rows are bra-major, i.e. row = m_a * (2 l_b + 1) + m_b. For symmetric
/// and antisymmetric matrices the transposed element is filled with the same
/// sign and the opposite sign respectively (off-diagonal entries only).
/// @param matrix The block sparse matrix.
/// @param screened_pairs The vector of screened basis function pairs defining
/// the block layout (indexed by the block keys).
/// @param basis The molecular basis defining the atomic orbital ordering.
/// @return The full sized dense matrix.
auto to_dense_matrix(const CSparseMatrix                            &matrix,
                     const std::vector<CScreenedBasisFunctionPair>  &screened_pairs,
                     const CMolecularBasis                          &basis) -> CDenseMatrix;

}  // namespace sparsematfunc

#endif /* SparseMatrixFunc_hpp */
