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

#ifndef LocalCorePotentialPrimRecSS_hpp
#define LocalCorePotentialPrimRecSS_hpp

#include "SimdArray.hpp"

namespace ecprec {  // ovlrec namespace

/// @brief Computes primitive [S|U_L|S]  integrals for set of data buffers.
/// @param pbuffer  The primitive integrals buffer.
/// @param idx_lpot_ss The index of integrals in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_r The vector of coordinates R.
/// @param r_a The Cartesian A point coordinates.
/// @param a_exp The primitive basis function exponent on center A.
/// @param a_norm The primitive basis function normalization factor on center A.
/// @param c_exp The local core potential exponent on center C.
/// @param c_fact The local core potential factor on center C.
auto comp_prim_local_core_potential_ss(      CSimdArray<double>& pbuffer,
                                       const size_t              idx_lpot_ss,
                                             CSimdArray<double>& factors,
                                       const size_t              idx_r,
                                       const TPoint<double>&     r_a,
                                       const double              a_exp,
                                       const double              a_norm,
                                       const double              c_exp,
                                       const double              c_fact) -> void;
}  // namespace ecprec

#endif /* LocalCorePotentialPrimRecSS_hpp */
