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
#include <cmath>
#include <string>

#include "ErrorHandler.hpp"
#include "MathConst.hpp"
#include "ScreeningFunc.hpp"
#include "SimdAlign.hpp"
#include "SimdDimensions.hpp"
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
    if ((bra.get_angular_momentum() != 0) || (ket.get_angular_momentum() != 0))
    {
        errors::assertMsgCritical(
            false, std::string("SimdOverlapFunc.compute_ss_overlap: Basis functions must be of zero angular momentum"));
    }

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

    constexpr auto fpi = mathconst::pi_value();

    // accumulate the integrals of each pair of primitives

    simdfunc::accumulate_primitives(bra, ket, dimensions, [&](const simdfunc::CPrimitivePair &pair) {
        const auto ncols = pair.ncols;

        const auto fexp = pair.aexp + pair.bexp;

        const auto fmu = pair.aexp * pair.bexp / fexp;

        const auto fovl = fpi / fexp;

        const auto ffact = pair.anorm * pair.bnorm * fovl * std::sqrt(fovl);

        // NOTE: the row of the buffer and the row of the coordinates start at a
        // cache line boundary, so the loop is vectorized with aligned loads and
        // stores. A pair of primitives contributes only to the atom pairs it
        // reaches, so the loop shortens as the primitives get tighter.

        // NOTE: the exponential is issued as a call to the vector math library of
        // the platform, so it does not break the vectorization of the loop. It
        // costs about nine tenths of the loop, which is therefore bound by the
        // throughput of the exponential and not by the memory it touches.

#pragma omp simd aligned(prim, ab_2 : simd::cache_line_size())
        for (size_t k = 0; k < ncols; k++)
        {
            prim[k] += ffact * std::exp(-fmu * ab_2[k]);
        }
    });

    // NOTE: the atom pairs beyond the reach of every pair of primitives have no
    // contribution and are set to zero.

    std::copy(prim, prim + nmax, values);

    std::fill(values + nmax, values + nvalues, 0.0);
}

}  // namespace simdovl
