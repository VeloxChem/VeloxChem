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



#ifndef SimdOverlapPrimRecSLP_hpp
#define SimdOverlapPrimRecSLP_hpp

#include "SimdPrimitives.hpp"

namespace simdovl {  // simdovl namespace

/// @brief Accumulates the prefactor of one pair of primitives of a combination of
/// one basis function of zero angular momentum and one of angular momentum one.
/// @param prim The accumulator of the combination, a row of the primitive buffer.
/// @param ab_2 The squared distances of the atom pairs, a row of the coordinates.
/// @param pair The exponents, the normalization factors and the number of atom
/// pairs of the pair of primitives.
/// @param on_ket True when the angular momentum sits on ket side, false when it
/// sits on bra side.
/// @note The harmonic of the vector between the atoms factors out of the sum over
/// the pairs of primitives, so only its scalar coefficient is accumulated here and
/// the angular components are formed once by the caller. The displacement of the
/// Gaussian product center from the atom carrying the angular momentum is a / p
/// times the vector between the atoms on ket side and -b / p on bra side, raised to
/// the angular momentum, which is the only way the two orders differ.
/// @note Both pointers are rows of a matrix of cache line aligned rows, which the
/// loop relies on to load and store with aligned instructions. A caller which hands
/// over anything else breaks that promise silently.
auto compute_prim_slp_overlap(double *prim, const double *ab_2, const simdfunc::CPrimitivePair &pair, const bool on_ket)
    -> void;

}  // namespace simdovl

#endif /* SimdOverlapPrimRecSLP_hpp */
