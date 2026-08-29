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


#ifndef SimdBoysFunc_hpp
#define SimdBoysFunc_hpp

#include "SimdVariableMatrix.hpp"

namespace simdfunc {  // simdfunc namespace

/// @brief Computes Boys function values for the arguments of a matrix.
/// @param matrix The matrix with the arguments of Boys function in block zero
/// and the values of Boys function of order zero to N in blocks one to N + 1,
/// where N is the order of Boys function.
/// @note The order of Boys function follows from the number of blocks, which is
/// the order plus two. The arguments in block zero are not modified.
/// @note The values are computed from the Taylor expansion on a grid of
/// arguments below 35.05 and from the asymptotic expansion above it, with the
/// remaining orders following from the downward and upward recursions
/// respectively. The relative accuracy is better than 1.0e-14 for every order up
/// to the highest supported one and for every argument.
auto compute_boys_function(CSimdVariableMatrix &matrix) -> void;

}  // namespace simdfunc

#endif /* SimdBoysFunc_hpp */
