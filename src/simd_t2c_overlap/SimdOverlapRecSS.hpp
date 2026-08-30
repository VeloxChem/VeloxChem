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



#ifndef SimdOverlapRecSS_hpp
#define SimdOverlapRecSS_hpp

#include <cstddef>

#include "BasisFunction.hpp"
#include "SimdMatrix.hpp"

namespace simdovl {  // simdovl namespace

/// @brief Computes the overlap integrals of a combination of basis functions of
/// zero angular momentum on bra and ket sides.
/// @param values The values of the combination of basis functions in the values
/// block of the sparsity pattern.
/// @param nvalues The number of values to compute, i.e. the number of atom pairs
/// surviving the screening of the combination of basis functions.
/// @param bra The basis function on bra side.
/// @param ket The basis function on ket side.
/// @param coordinates The coordinates of the atom pairs, as seven rows ordered by
/// ascending interatomic distance, the last of which holds the squared distance
/// of the atom pair. Only the leading nvalues columns are read, as the atom pairs
/// surviving the screening are the closest ones.
/// @param threshold The screening threshold of the integrals.
/// @note The integrals of all pairs of primitives are accumulated in a single
/// row, which spans only the atom pairs reached by the furthest reaching pair of
/// primitives.
/// @note The pairs of primitives are screened with the threshold of the integrals
/// divided by their number, so that the accumulated contributions of the
/// discarded pairs of primitives stay below the threshold of the integrals.
auto compute_ss_overlap(double               *values,
                        const size_t          nvalues,
                        const CBasisFunction &bra,
                        const CBasisFunction &ket,
                        const CSimdMatrix    &coordinates,
                        const double          threshold) -> void;

}  // namespace simdovl

#endif /* SimdOverlapRecSS_hpp */
