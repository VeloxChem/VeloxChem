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

#ifndef fika_veloxchem_order_hpp
#define fika_veloxchem_order_hpp

#include <cstddef>
#include <vector>

#include "fika_molecular_basis.hpp"
#include "fika_dense_matrix.hpp"

namespace fika {

/// VeloxChem index of each basis function of `basis` (indexed by fika's function index).
/// VeloxChem orders the functions by angular momentum l, then component m = -l..l, then atom,
/// then contracted function (the shells of an atom with equal l in basis order); fika by atom,
/// shell, contracted function, then m. Both use the same real solid harmonics (order and sign of
/// the m components), so the orders differ by this permutation only.
auto veloxchem_order(const MolecularBasis& basis) -> std::vector<std::size_t>;

/// `matrix` (over the functions of `basis` in VeloxChem's order) in fika's order, of the same
/// symmetry: result(i, j) = matrix(v_i, v_j) with v = veloxchem_order(basis). Throws
/// std::invalid_argument unless the matrix is function_count() x function_count().
auto veloxchem_to_fika(const DenseMatrix& matrix, const MolecularBasis& basis) -> DenseMatrix;

/// The inverse of veloxchem_to_fika: `matrix` in fika's order to VeloxChem's, result(v_i, v_j) =
/// matrix(i, j). Throws std::invalid_argument unless the matrix is function_count() x
/// function_count().
auto fika_to_veloxchem(const DenseMatrix& matrix, const MolecularBasis& basis) -> DenseMatrix;

}  // namespace fika

#endif  // fika_veloxchem_order_hpp
