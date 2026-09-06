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



#include "SimdOverlapRecSLI.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <ranges>
#include <string>

#include "ErrorHandler.hpp"
#include "MathConst.hpp"
#include "ScreeningFunc.hpp"
#include "SimdAlign.hpp"
#include "SimdDimensions.hpp"
#include "SimdPrimitives.hpp"

namespace simdovl {  // simdovl namespace

auto
compute_sli_overlap(double               *values,
                    const size_t          nvalues,
                    const CBasisFunction &bra,
                    const CBasisFunction &ket,
                    const CSimdMatrix    &coordinates,
                    const double          threshold) -> void
{
    const auto lbra = bra.get_angular_momentum();

    const auto lket = ket.get_angular_momentum();

    if (!(((lbra == 0) && (lket == 6)) || ((lbra == 6) && (lket == 0))))
    {
        errors::assertMsgCritical(
            false, std::string("SimdOverlapRecSLI.compute_sli_overlap: Basis functions must be of angular momenta zero and six"));
    }

    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("SimdOverlapRecSLI.compute_sli_overlap: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto nprims = bra.exponents().size() * ket.exponents().size();

    // NOTE: the pairs of primitives are screened with the threshold of the
    // integrals divided by their number, as their contributions accumulate into
    // a single value and the error of the sum is bounded by the number of terms.

    const auto dimensions = simdfunc::make_column_dimensions(
        bra, ket, nvalues, coordinates, screenfunc::two_center_overlap_primitive_bound, threshold / static_cast<double>(nprims));

    // NOTE: the buffer holds the contracted prefactor alone, as the harmonic
    // factors out of the sum over the pairs of primitives and the integrals of the
    // angular components are formed straight into the values.

    auto buffer = simdfunc::make_primitive_buffer(dimensions, 1);

    if (buffer.number_of_columns() == 0)
    {
        std::fill(values, values + 13 * nvalues, 0.0);

        return;
    }

    const auto nmax = buffer.number_of_columns();

    auto *pe_0 = buffer.data(0);

    // NOTE: the components of the vector between the atoms and its squared length
    // are carried by the coordinates, so the angular half below reads rows which
    // are already in place.

    const auto *ab_x = coordinates.data(6);
    const auto *ab_y = coordinates.data(7);
    const auto *ab_z = coordinates.data(8);

    const auto *ab_2 = coordinates.data(9);

    constexpr auto fpi = mathconst::pi_value();

    // NOTE: the harmonic sits on whichever side carries the angular momentum, and
    // the Gaussian product center is displaced from it by (a / p) times the vector
    // between the atoms when that is the ket side and by -(b / p) when it is the bra
    // side. The order is therefore settled once here and not inside any loop.

    const auto on_ket = (lbra == 0);

    // accumulate the prefactor of each pair of primitives

    simdfunc::accumulate_primitives(bra, ket, dimensions, [&](const simdfunc::CPrimitivePair &pair) {
        const auto ncols = pair.ncols;

        const auto fexp = pair.aexp + pair.bexp;

        const auto fmu = pair.aexp * pair.bexp / fexp;

        const auto fovl = fpi / fexp;

        const auto fbase = pair.anorm * pair.bnorm * fovl * std::sqrt(fovl);

        const auto fr = on_ket ? (pair.aexp / fexp) : (-pair.bexp / fexp);

        const auto ffact = fbase * fr * fr * fr * fr * fr * fr;

#pragma omp simd aligned(pe_0, ab_2 : simd::cache_line_size())
        for (size_t k = 0; k < ncols; k++)
        {
            pe_0[k] += ffact * std::exp(-fmu * ab_2[k]);
        }
    });

    // NOTE: the rows of the values are not aligned, as they start at the offset of
    // this combination of basis functions in the values block, so they are kept out
    // of the aligned clauses below.

    auto *pc_0 = values + 0 * nvalues;
    auto *pc_1 = values + 1 * nvalues;
    auto *pc_2 = values + 2 * nvalues;
    auto *pc_3 = values + 3 * nvalues;
    auto *pc_4 = values + 4 * nvalues;
    auto *pc_5 = values + 5 * nvalues;
    auto *pc_6 = values + 6 * nvalues;
    auto *pc_7 = values + 7 * nvalues;
    auto *pc_8 = values + 8 * nvalues;
    auto *pc_9 = values + 9 * nvalues;
    auto *pc_10 = values + 10 * nvalues;
    auto *pc_11 = values + 11 * nvalues;
    auto *pc_12 = values + 12 * nvalues;

    // NOTE: the components are formed in 4 loops, as the vectorizer runs out
    // of registers with all of them in one. Only the prefactor and the vector
    // between the atoms are loaded by more than one loop.

#pragma omp simd aligned(pe_0, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];

        pc_0[k] = e_0 * (std::sqrt(16.2421875) * x * x * x * x * x * y - std::sqrt(180.46875) * x * x * x * y * y * y + std::sqrt(16.2421875) * x * y * y * y * y * y);

        pc_1[k] = e_0 * (std::sqrt(135.3515625) * x * x * x * x * y * z - std::sqrt(541.40625) * x * x * y * y * y * z + std::sqrt(5.4140625) * y * y * y * y * y * z);

        pc_2[k] = e_0 * (-std::sqrt(3.9375) * x * x * x * x * x * y + std::sqrt(393.75) * x * x * x * y * z * z + std::sqrt(3.9375) * x * y * y * y * y * y - std::sqrt(393.75) * x * y * y * y * z * z);

        pc_3[k] = e_0 * (-std::sqrt(66.4453125) * x * x * x * x * y * z - std::sqrt(29.53125) * x * x * y * y * y * z + std::sqrt(472.5) * x * x * y * z * z * z + std::sqrt(7.3828125) * y * y * y * y * y * z - std::sqrt(52.5) * y * y * y * z * z * z);
    }

#pragma omp simd aligned(pe_0, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];

        pc_4[k] = e_0 * (std::sqrt(0.8203125) * x * x * x * x * x * y + std::sqrt(3.28125) * x * x * x * y * y * y - std::sqrt(210.0) * x * x * x * y * z * z + std::sqrt(0.8203125) * x * y * y * y * y * y - std::sqrt(210.0) * x * y * y * y * z * z + std::sqrt(210.0) * x * y * z * z * z * z);

        pc_5[k] = e_0 * (std::sqrt(8.203125) * x * x * x * x * y * z + std::sqrt(32.8125) * x * x * y * y * y * z - std::sqrt(131.25) * x * x * y * z * z * z + std::sqrt(8.203125) * y * y * y * y * y * z - std::sqrt(131.25) * y * y * y * z * z * z + std::sqrt(21.0) * y * z * z * z * z * z);

        pc_6[k] = e_0 * (-0.3125 * x * x * x * x * x * x - 0.9375 * x * x * x * x * y * y + 5.625 * x * x * x * x * z * z - 0.9375 * x * x * y * y * y * y + 11.25 * x * x * y * y * z * z - 7.5 * x * x * z * z * z * z - 0.3125 * y * y * y * y * y * y + 5.625 * y * y * y * y * z * z - 7.5 * y * y * z * z * z * z + z * z * z * z * z * z);

        pc_7[k] = e_0 * (std::sqrt(8.203125) * x * x * x * x * x * z + std::sqrt(32.8125) * x * x * x * y * y * z - std::sqrt(131.25) * x * x * x * z * z * z + std::sqrt(8.203125) * x * y * y * y * y * z - std::sqrt(131.25) * x * y * y * z * z * z + std::sqrt(21.0) * x * z * z * z * z * z);
    }

#pragma omp simd aligned(pe_0, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];

        pc_8[k] = e_0 * (std::sqrt(0.205078125) * x * x * x * x * x * x + std::sqrt(0.205078125) * x * x * x * x * y * y - std::sqrt(52.5) * x * x * x * x * z * z - std::sqrt(0.205078125) * x * x * y * y * y * y + std::sqrt(52.5) * x * x * z * z * z * z - std::sqrt(0.205078125) * y * y * y * y * y * y + std::sqrt(52.5) * y * y * y * y * z * z - std::sqrt(52.5) * y * y * z * z * z * z);

        pc_9[k] = e_0 * (-std::sqrt(7.3828125) * x * x * x * x * x * z + std::sqrt(29.53125) * x * x * x * y * y * z + std::sqrt(52.5) * x * x * x * z * z * z + std::sqrt(66.4453125) * x * y * y * y * y * z - std::sqrt(472.5) * x * y * y * z * z * z);

        pc_10[k] = e_0 * (-std::sqrt(0.24609375) * x * x * x * x * x * x + std::sqrt(6.15234375) * x * x * x * x * y * y + std::sqrt(24.609375) * x * x * x * x * z * z + std::sqrt(6.15234375) * x * x * y * y * y * y - std::sqrt(885.9375) * x * x * y * y * z * z - std::sqrt(0.24609375) * y * y * y * y * y * y + std::sqrt(24.609375) * y * y * y * y * z * z);

        pc_11[k] = e_0 * (std::sqrt(5.4140625) * x * x * x * x * x * z - std::sqrt(541.40625) * x * x * x * y * y * z + std::sqrt(135.3515625) * x * y * y * y * y * z);
    }

#pragma omp simd aligned(pe_0, ab_x, ab_y : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];

        const auto e_0 = pe_0[k];

        pc_12[k] = e_0 * (std::sqrt(0.451171875) * x * x * x * x * x * x - std::sqrt(101.513671875) * x * x * x * x * y * y + std::sqrt(101.513671875) * x * x * y * y * y * y - std::sqrt(0.451171875) * y * y * y * y * y * y);
    }

    // NOTE: the atom pairs beyond the reach of every pair of primitives have no
    // contribution and are set to zero.

    for (size_t m = 0; m < 13; m++)
    {
        auto *pv = values + m * nvalues;

        std::fill(pv + nmax, pv + nvalues, 0.0);
    }
}

}  // namespace simdovl
