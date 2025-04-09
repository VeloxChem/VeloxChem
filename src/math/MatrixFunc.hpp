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

#ifndef MatrixFunc_hpp
#define MatrixFunc_hpp

#include "Matrix.hpp"
#include "MolecularBasis.hpp"

namespace matfunc {

/// @brief Creates a matrix.
/// @param basis The molecular basis to create matrix.
/// @param mtype  The matrix type of created matrix.
/// @return The matrix.
auto make_matrix(const CMolecularBasis& basis, const mat_t mtype) -> CMatrix;

/// @brief Creates a matrix.
/// @param bra_basis The molecular basis to create matrix (rows data).
/// @param ket_basis The molecular basis to create matrix (columns data).
/// @return The matrix.
auto make_matrix(const CMolecularBasis& bra_basis, const CMolecularBasis& ket_basis) -> CMatrix;

/// @brief Creates matrix.
/// @param label The exchange-correlation functional type label.
/// @param nrows The number of rows in submatrix.
/// @param ncols The number of columns in submatrix.
/// @return The matrix.
auto make_matrix(const std::string& label, const size_t nrows, const size_t ncols) -> CMatrix;

}  // namespace matfunc

#endif /* MatrixFunc_hpp */
