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

#ifndef ObaraSaikaFunc_hpp
#define ObaraSaikaFunc_hpp

#include "DenseMatrix.hpp"
#include "ScreenedBasisFunctionPair.hpp"

namespace osfunc {  // osfunc namespace

/// @brief Computes the Obara-Saika PA distances for a screened basis function
/// pair: PA = -beta / (alpha + beta) * (A - B), where alpha and beta are the
/// bra and ket primitive exponents and A and B are the bra and ket centers.
/// The result is stored in a dense matrix with 3 * nprim_bra * nprim_ket rows
/// and number_of_pairs columns. Rows are component-major:
///   row = c * (nprim_bra * nprim_ket) + i_bra * nprim_ket + j_ket,
/// with c the Cartesian component (0 = x, 1 = y, 2 = z), i_bra the bra
/// primitive index and j_ket the ket primitive index. Each column is a
/// surviving atom pair.
/// @param pair The screened basis function pair.
/// @return The dense matrix of PA distances.
auto compute_pa(const CScreenedBasisFunctionPair &pair) -> CDenseMatrix;

/// @brief Computes the Obara-Saika PB distances for a screened basis function
/// pair: PB = alpha / (alpha + beta) * (A - B), where alpha and beta are the
/// bra and ket primitive exponents and A and B are the bra and ket centers.
/// The result is stored in a dense matrix with 3 * nprim_bra * nprim_ket rows
/// and number_of_pairs columns, with the same component-major row layout as
/// compute_pa.
/// @param pair The screened basis function pair.
/// @return The dense matrix of PB distances.
auto compute_pb(const CScreenedBasisFunctionPair &pair) -> CDenseMatrix;

}  // namespace osfunc

#endif /* ObaraSaikaFunc_hpp */
