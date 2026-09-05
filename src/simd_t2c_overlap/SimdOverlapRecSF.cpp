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



#include "SimdOverlapRecSF.hpp"

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
compute_sf_overlap(double               *values,
                   const size_t          nvalues,
                   const CBasisFunction &bra,
                   const CBasisFunction &ket,
                   const CSimdMatrix    &coordinates,
                   const double          threshold) -> void
{
    if ((bra.get_angular_momentum() != 0) || (ket.get_angular_momentum() != 3))
    {
        errors::assertMsgCritical(
            false, std::string("SimdOverlapRecSF.compute_sf_overlap: Basis functions must be of angular momenta zero and three"));
    }

    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("SimdOverlapRecSF.compute_sf_overlap: Number of values exceeds number of atom pairs"));
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

    auto buffer = simdfunc::make_primitive_buffer(dimensions, 8);

    if (buffer.number_of_columns() == 0)
    {
        simdfunc::store_components(values, nvalues, buffer, 1, 7);

        return;
    }

    const auto nmax = buffer.number_of_columns();

    auto *prim = buffer.data(0);

    auto *out_m3 = buffer.data(1);
    auto *out_m2 = buffer.data(2);
    auto *out_m1 = buffer.data(3);
    auto *out_0 = buffer.data(4);
    auto *out_p1 = buffer.data(5);
    auto *out_p2 = buffer.data(6);
    auto *out_p3 = buffer.data(7);

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

        // NOTE: the harmonic sits on the ket side, so the Gaussian product center is
        // displaced from it by (a / p) times the vector between the atoms and the
        // prefactor carries that ratio raised to the power three.

        const auto fr = pair.aexp / fexp;

        const auto ffact = pair.anorm * pair.bnorm * fovl * std::sqrt(fovl) * fr * fr * fr;

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

#pragma omp simd aligned(out_m3, out_m2, out_m1, out_0, out_p1, out_p2, out_p3, prim, ab_x, ab_y, ab_z, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];
        const auto r_2 = ab_2[k];

        out_m3[k] = prim[k] * (std::sqrt(5.625) * x * x * y - std::sqrt(0.625) * y * y * y);

        out_m2[k] = prim[k] * (std::sqrt(15.0) * x * y * z);

        out_m1[k] = prim[k] * (std::sqrt(9.375) * y * z * z - std::sqrt(0.375) * y * r_2);

        out_0[k] = prim[k] * (2.5 * z * z * z - 1.5 * z * r_2);

        out_p1[k] = prim[k] * (std::sqrt(9.375) * x * z * z - std::sqrt(0.375) * x * r_2);

        out_p2[k] = prim[k] * (std::sqrt(3.75) * x * x * z - std::sqrt(3.75) * y * y * z);

        out_p3[k] = prim[k] * (std::sqrt(0.625) * x * x * x - std::sqrt(5.625) * x * y * y);
    }

    simdfunc::store_components(values, nvalues, buffer, 1, 7);
}

}  // namespace simdovl
