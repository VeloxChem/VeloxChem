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


#ifndef SimdPrimitives_hpp
#define SimdPrimitives_hpp

#include <algorithm>
#include <cstddef>
#include <ranges>
#include <string>
#include <vector>

#include "ErrorHandler.hpp"
#include "SimdMatrix.hpp"

namespace simdfunc {  // simdfunc namespace

/// @brief Creates the buffer the pairs of primitives accumulate their
/// contributions in.
/// @param dimensions The number of atom pairs each pair of primitives reaches.
/// @param nrows The number of accumulators of the kernel.
/// @return The zeroed matrix of nrows rows spanning the atom pairs reached by the
/// pair of primitives reaching furthest, empty when no pair of primitives reaches
/// any atom pair.
/// @note The buffer spans the atom pairs reached by the pair of primitives
/// reaching furthest, which is searched for rather than assumed. The primitives
/// are sorted by descending exponent, but the bound of a pair of primitives
/// carries their prefactor as well as their decay, so a tighter pair with a larger
/// prefactor reaches further than a more diffuse pair with a smaller one, and the
/// last pair is not always the furthest reaching.
inline auto
make_primitive_buffer(const std::vector<size_t> &dimensions, const size_t nrows) -> CSimdMatrix
{
    errors::assertMsgCritical(nrows > 0, std::string("SimdPrimitives.make_primitive_buffer: Number of rows must be positive"));

    if (dimensions.empty()) return CSimdMatrix();

    const auto nmax = *std::ranges::max_element(dimensions);

    if (nmax == 0) return CSimdMatrix();

    auto matrix = CSimdMatrix(nrows, nmax);

    matrix.zero();

    return matrix;
}

}  // namespace simdfunc

#endif /* SimdPrimitives_hpp */
