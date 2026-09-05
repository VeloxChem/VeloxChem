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



#include "SimdOverlapRecHS.hpp"

#include <algorithm>
#include <cmath>
#include <ranges>
#include <string>

#include "ErrorHandler.hpp"
#include "MathConst.hpp"
#include "ScreeningFunc.hpp"
#include "SimdAlign.hpp"
#include "SimdDimensions.hpp"
#include "SimdPrimitives.hpp"
#include "SimdStorage.hpp"

namespace simdovl {  // simdovl namespace

auto
compute_hs_overlap(double               *values,
                   const size_t          nvalues,
                   const CBasisFunction &bra,
                   const CBasisFunction &ket,
                   const CSimdMatrix    &coordinates,
                   const double          threshold) -> void
{
    if ((bra.get_angular_momentum() != 5) || (ket.get_angular_momentum() != 0))
    {
        errors::assertMsgCritical(
            false, std::string("SimdOverlapRecHS.compute_hs_overlap: Basis functions must be of angular momenta five and zero"));
    }

    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("SimdOverlapRecHS.compute_hs_overlap: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto nprims = bra.exponents().size() * ket.exponents().size();

    // NOTE: the pairs of primitives are screened with the threshold of the
    // integrals divided by their number, as their contributions accumulate into
    // a single value and the error of the sum is bounded by the number of terms.

    const auto dimensions = simdfunc::make_column_dimensions(
        bra, ket, nvalues, coordinates, screenfunc::two_center_overlap_primitive_bound, threshold / static_cast<double>(nprims));

    // NOTE: the buffer holds the prefactor shared by the angular components in
    // its first row and the integrals of the components in the rows which follow,
    // as the harmonic factors out of the sum over the pairs of primitives and
    // multiplies the accumulated prefactor once.

    auto buffer = simdfunc::make_primitive_buffer(dimensions, 12);

    if (buffer.number_of_columns() == 0)
    {
        simdfunc::store_components(values, nvalues, buffer, 1, 11);

        return;
    }

    const auto nmax = buffer.number_of_columns();

    auto *prim = buffer.data(0);

    auto *out_m5 = buffer.data(1);
    auto *out_m4 = buffer.data(2);
    auto *out_m3 = buffer.data(3);
    auto *out_m2 = buffer.data(4);
    auto *out_m1 = buffer.data(5);
    auto *out_0 = buffer.data(6);
    auto *out_p1 = buffer.data(7);
    auto *out_p2 = buffer.data(8);
    auto *out_p3 = buffer.data(9);
    auto *out_p4 = buffer.data(10);
    auto *out_p5 = buffer.data(11);

    // NOTE: the components of the vector between the atoms and its squared length
    // are carried by the coordinates, so the harmonic below is formed from rows
    // which are already in place.

    const auto *ab_x = coordinates.data(6);
    const auto *ab_y = coordinates.data(7);
    const auto *ab_z = coordinates.data(8);

    const auto *ab_2 = coordinates.data(9);

    constexpr auto fpi = mathconst::pi_value();

    // accumulate the prefactor of each pair of primitives

    simdfunc::accumulate_primitives(bra, ket, dimensions, [&](const simdfunc::CPrimitivePair &pair) {
        const auto ncols = pair.ncols;

        const auto fexp = pair.aexp + pair.bexp;

        const auto fmu = pair.aexp * pair.bexp / fexp;

        const auto fovl = fpi / fexp;

        // NOTE: the harmonic sits on the bra side, so the displacement is -(b / p) times
        // the vector between the atoms and the prefactor carries that ratio, sign
        // included, raised to the power five.

        const auto fr = -pair.bexp / fexp;

        const auto ffact = pair.anorm * pair.bnorm * fovl * std::sqrt(fovl) * fr * fr * fr * fr * fr;

#pragma omp simd aligned(prim, ab_2 : simd::cache_line_size())
        for (size_t k = 0; k < ncols; k++)
        {
            prim[k] += ffact * std::exp(-fmu * ab_2[k]);
        }
    });

    // NOTE: the integrals of the angular components are the accumulated prefactor
    // times the components of the harmonic, formed in one pass over the rows of
    // the buffer and of the coordinates, all of which start at a cache line
    // boundary.

    // NOTE: the components are formed in 3 loops, as the vectorizer runs out
    // of registers with all 11 of them in one. Only the accumulated prefactor and
    // the vector between the atoms are loaded by more than one loop.

#pragma omp simd aligned(out_m5, out_m4, out_m3, out_m2, prim, ab_x, ab_y, ab_z, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];
        const auto r_2 = ab_2[k];

        out_m5[k] = prim[k] * (std::sqrt(12.3046875) * x * x * x * x * y - std::sqrt(49.21875) * x * x * y * y * y + std::sqrt(0.4921875) * y * y * y * y * y);

        out_m4[k] = prim[k] * (std::sqrt(78.75) * x * x * x * y * z - std::sqrt(78.75) * x * y * y * y * z);

        out_m3[k] = prim[k] * (std::sqrt(199.3359375) * x * x * y * z * z - std::sqrt(2.4609375) * x * x * y * r_2 - std::sqrt(22.1484375) * y * y * y * z * z + std::sqrt(0.2734375) * y * y * y * r_2);

        out_m2[k] = prim[k] * (std::sqrt(236.25) * x * y * z * z * z - std::sqrt(26.25) * x * y * z * r_2);
    }

#pragma omp simd aligned(out_m1, out_0, out_p1, out_p2, prim, ab_x, ab_y, ab_z, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];
        const auto r_2 = ab_2[k];

        out_m1[k] = prim[k] * (std::sqrt(103.359375) * y * z * z * z * z - std::sqrt(45.9375) * y * z * z * r_2 + std::sqrt(0.234375) * y * r_2 * r_2);

        out_0[k] = prim[k] * (7.875 * z * z * z * z * z - 8.75 * z * z * z * r_2 + 1.875 * z * r_2 * r_2);

        out_p1[k] = prim[k] * (std::sqrt(103.359375) * x * z * z * z * z - std::sqrt(45.9375) * x * z * z * r_2 + std::sqrt(0.234375) * x * r_2 * r_2);

        out_p2[k] = prim[k] * (std::sqrt(59.0625) * x * x * z * z * z - std::sqrt(6.5625) * x * x * z * r_2 - std::sqrt(59.0625) * y * y * z * z * z + std::sqrt(6.5625) * y * y * z * r_2);
    }

#pragma omp simd aligned(out_p3, out_p4, out_p5, prim, ab_x, ab_y, ab_z, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];
        const auto r_2 = ab_2[k];

        out_p3[k] = prim[k] * (std::sqrt(22.1484375) * x * x * x * z * z - std::sqrt(0.2734375) * x * x * x * r_2 - std::sqrt(199.3359375) * x * y * y * z * z + std::sqrt(2.4609375) * x * y * y * r_2);

        out_p4[k] = prim[k] * (std::sqrt(4.921875) * x * x * x * x * z - std::sqrt(177.1875) * x * x * y * y * z + std::sqrt(4.921875) * y * y * y * y * z);

        out_p5[k] = prim[k] * (std::sqrt(0.4921875) * x * x * x * x * x - std::sqrt(49.21875) * x * x * x * y * y + std::sqrt(12.3046875) * x * y * y * y * y);
    }

    simdfunc::store_components(values, nvalues, buffer, 1, 11);
}

}  // namespace simdovl
