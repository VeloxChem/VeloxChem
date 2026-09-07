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



#include "SimdOverlapRecSS.hpp"

#include <algorithm>
#include <ranges>
#include <string>

#include "ErrorHandler.hpp"
#include "ScreeningFunc.hpp"
#include "SimdDimensions.hpp"
#include "SimdOverlapPrimRecSS.hpp"
#include "SimdPrimitives.hpp"

namespace simdovl {  // simdovl namespace

auto
compute_ss_overlap(double               *values,
                   const size_t          nvalues,
                   const CBasisFunction &bra,
                   const CBasisFunction &ket,
                   const CSimdMatrix    &coordinates,
                   const double          threshold) -> void
{
    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("SimdOverlapFunc.compute_ss_overlap: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto nprims = bra.exponents().size() * ket.exponents().size();

    // NOTE: the pairs of primitives are screened with the threshold of the
    // integrals divided by their number, as their contributions accumulate into
    // a single value and the error of the sum is bounded by the number of terms.

    const auto dimensions = simdfunc::make_column_dimensions(
        bra, ket, nvalues, coordinates, screenfunc::two_center_overlap_primitive_bound, threshold / static_cast<double>(nprims));

    auto buffer = simdfunc::make_primitive_buffer(dimensions, 1);

    if (buffer.number_of_columns() == 0)
    {
        std::fill(values, values + nvalues, 0.0);

        return;
    }

    // NOTE: the integrals of all pairs of primitives are accumulated in a single
    // row, which starts at a cache line boundary and spans only the atom pairs
    // reached by the furthest reaching pair of primitives.

    const auto nmax = buffer.number_of_columns();

    auto *prim = buffer.data(0);

    // NOTE: the squared distances of the atom pairs are carried by the
    // coordinates, so that they are formed once for the whole block instead of
    // once for every combination of basis functions.

    const auto *ab_2 = coordinates.data(9);

    // accumulate the integrals of each pair of primitives

    simdfunc::accumulate_primitives(
        bra, ket, dimensions, [&](const simdfunc::CPrimitivePair &pair) { compute_prim_ss_overlap(prim, ab_2, pair); });

    // NOTE: the atom pairs beyond the reach of every pair of primitives have no
    // contribution and are set to zero.

    std::copy(prim, prim + nmax, values);

    std::fill(values + nmax, values + nvalues, 0.0);
}

}  // namespace simdovl
