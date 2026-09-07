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

#include "SimdOverlapRecSLF.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <ranges>
#include <string>

#include "ErrorHandler.hpp"
#include "ScreeningFunc.hpp"
#include "SimdAlign.hpp"
#include "SimdDimensions.hpp"
#include "SimdOverlapPrimRecSLF.hpp"
#include "SimdPrimitives.hpp"

namespace simdovl {  // simdovl namespace

auto
compute_slf_overlap(double               *values,
                    const size_t          nvalues,
                    const CBasisFunction &bra,
                    const CBasisFunction &ket,
                    const CSimdMatrix    &coordinates,
                    const double          threshold) -> void
{
    // NOTE: the side which carries the angular momentum is read off the bra alone,
    // as the dispatcher reaches this kernel only for the two orders of one basis
    // function of zero angular momentum and one of higher.

    const auto lbra = bra.get_angular_momentum();

    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("SimdOverlapRecSLF.compute_slf_overlap: Number of values exceeds number of atom pairs"));
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
        std::fill(values, values + 7 * nvalues, 0.0);

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

    // NOTE: the harmonic sits on whichever side carries the angular momentum, and
    // the Gaussian product center is displaced from it by (a / p) times the vector
    // between the atoms when that is the ket side and by -(b / p) when it is the bra
    // side. The order is therefore settled once here and not inside any loop.

    const auto on_ket = (lbra == 0);

    // accumulate the prefactor of each pair of primitives

    simdfunc::accumulate_primitives(bra, ket, dimensions, [&](const simdfunc::CPrimitivePair &pair) {
        compute_prim_slf_overlap(pe_0, ab_2, pair, on_ket);
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

    // NOTE: the components are formed in 2 loops, as the vectorizer runs out
    // of registers with all of them in one. Only the prefactor and the vector
    // between the atoms are loaded by more than one loop.

#pragma omp simd aligned(pe_0, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto x2 = x * x;
        const auto y2 = y * y;
        const auto y3 = y2 * y;

        const auto e_0 = pe_0[k];

        pc_0[k] = e_0 * (std::sqrt(5.625) * x2 * y
                       - std::sqrt(0.625) * y3);

        pc_1[k] = e_0 * (std::sqrt(15.0) * x * y * z);

        pc_2[k] = e_0 * (-std::sqrt(0.375) * x2 * y
                       - std::sqrt(0.375) * y3
                       + std::sqrt(6.0) * y * z * z);

        pc_3[k] = e_0 * (-1.5 * x2 * z
                       - 1.5 * y2 * z
                       + z * z * z);
    }

#pragma omp simd aligned(pe_0, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto x2 = x * x;
        const auto x3 = x2 * x;
        const auto y2 = y * y;

        const auto e_0 = pe_0[k];

        pc_4[k] = e_0 * (-std::sqrt(0.375) * x3
                       - std::sqrt(0.375) * x * y2
                       + std::sqrt(6.0) * x * z * z);

        pc_5[k] = e_0 * (std::sqrt(3.75) * x2 * z
                       - std::sqrt(3.75) * y2 * z);

        pc_6[k] = e_0 * (std::sqrt(0.625) * x3
                       - std::sqrt(5.625) * x * y2);
    }

    // NOTE: the atom pairs beyond the reach of every pair of primitives have no
    // contribution and are set to zero.

    for (size_t m = 0; m < 7; m++)
    {
        auto *pv = values + m * nvalues;

        std::fill(pv + nmax, pv + nvalues, 0.0);
    }
}

}  // namespace simdovl
