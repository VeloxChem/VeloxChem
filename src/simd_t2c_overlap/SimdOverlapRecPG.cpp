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



#include "SimdOverlapRecPG.hpp"

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
compute_pg_overlap(double               *values,
                   const size_t          nvalues,
                   const CBasisFunction &bra,
                   const CBasisFunction &ket,
                   const CSimdMatrix    &coordinates,
                   const double          threshold) -> void
{
    if ((bra.get_angular_momentum() != 1) || (ket.get_angular_momentum() != 4))
    {
        errors::assertMsgCritical(
            false, std::string("SimdOverlapRecPG.compute_pg_overlap: Basis functions must be of angular momenta one and four"));
    }

    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("SimdOverlapRecPG.compute_pg_overlap: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto nprims = bra.exponents().size() * ket.exponents().size();

    // NOTE: the pairs of primitives are screened with the threshold of the
    // integrals divided by their number, as their contributions accumulate into
    // a single value and the error of the sum is bounded by the number of terms.

    const auto dimensions = simdfunc::make_column_dimensions(
        bra, ket, nvalues, coordinates, screenfunc::two_center_overlap_primitive_bound, threshold / static_cast<double>(nprims));

    // NOTE: the buffer holds the contracted prefactors of the terms alone, as the
    // integrals of the angular components are formed straight into the values and
    // are not written a second time.

    auto buffer = simdfunc::make_primitive_buffer(dimensions, 2);

    if (buffer.number_of_columns() == 0)
    {
        std::fill(values, values + 27 * nvalues, 0.0);

        return;
    }

    const auto nmax = buffer.number_of_columns();

    auto *pe_0 = buffer.data(0);
    auto *pe_1 = buffer.data(1);

    // NOTE: the components of the vector between the atoms and its squared length
    // are carried by the coordinates, so the angular half below reads rows which
    // are already in place.

    const auto *ab_x = coordinates.data(6);
    const auto *ab_y = coordinates.data(7);
    const auto *ab_z = coordinates.data(8);

    const auto *ab_2 = coordinates.data(9);

    constexpr auto fpi = mathconst::pi_value();

    // accumulate the prefactor of each term over the pairs of primitives

    simdfunc::accumulate_primitives(bra, ket, dimensions, [&](const simdfunc::CPrimitivePair &pair) {
        const auto ncols = pair.ncols;

        const auto fexp = pair.aexp + pair.bexp;

        const auto fmu = pair.aexp * pair.bexp / fexp;

        const auto fovl = fpi / fexp;

        const auto fbase = pair.anorm * pair.bnorm * fovl * std::sqrt(fovl);

        // NOTE: the Gaussian product center is displaced from the atom on bra side
        // by fal times the vector between the atoms and from the atom on ket side by
        // fbe times it, and fh is the second moment the integration over that center
        // leaves behind.

        const auto fal = -pair.bexp / fexp;

        const auto fbe = pair.aexp / fexp;

        const auto fh = 0.5 / fexp;

        const auto f_0 = fbase * fal * fbe * fbe * fbe * fbe;

        const auto f_1 = fbase * fbe * fbe * fbe * fh;

        // NOTE: the exponential depends on the pair of primitives alone, so it is
        // evaluated once and shared by the prefactors of all terms.

#pragma omp simd aligned(pe_0, pe_1, ab_2 : simd::cache_line_size())
        for (size_t k = 0; k < ncols; k++)
        {
            const auto fss = std::exp(-fmu * ab_2[k]);

            pe_0[k] += f_0 * fss;
            pe_1[k] += f_1 * fss;
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
    auto *pc_13 = values + 13 * nvalues;
    auto *pc_14 = values + 14 * nvalues;
    auto *pc_15 = values + 15 * nvalues;
    auto *pc_16 = values + 16 * nvalues;
    auto *pc_17 = values + 17 * nvalues;
    auto *pc_18 = values + 18 * nvalues;
    auto *pc_19 = values + 19 * nvalues;
    auto *pc_20 = values + 20 * nvalues;
    auto *pc_21 = values + 21 * nvalues;
    auto *pc_22 = values + 22 * nvalues;
    auto *pc_23 = values + 23 * nvalues;
    auto *pc_24 = values + 24 * nvalues;
    auto *pc_25 = values + 25 * nvalues;
    auto *pc_26 = values + 26 * nvalues;

    // NOTE: the components are formed in 7 loops, as the vectorizer runs out
    // of registers with all of them in one. Only the prefactors and the vector
    // between the atoms are loaded by more than one loop.

#pragma omp simd aligned(pe_0, pe_1, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];

        pc_0[k] = e_0 * (std::sqrt(8.75) * x * x * x * y * y - std::sqrt(8.75) * x * y * y * y * y) + e_1 * (std::sqrt(8.75) * x * x * x - std::sqrt(78.75) * x * y * y);

        pc_1[k] = e_0 * (std::sqrt(39.375) * x * x * y * y * z - std::sqrt(4.375) * y * y * y * y * z) + e_1 * (std::sqrt(39.375) * x * x * z - std::sqrt(39.375) * y * y * z);

        pc_2[k] = e_0 * (-std::sqrt(1.25) * x * x * x * y * y - std::sqrt(1.25) * x * y * y * y * y + std::sqrt(45.0) * x * y * y * z * z) + e_1 * (-std::sqrt(1.25) * x * x * x - std::sqrt(11.25) * x * y * y + std::sqrt(45.0) * x * z * z);

        pc_3[k] = e_0 * (-std::sqrt(5.625) * x * x * y * y * z - std::sqrt(5.625) * y * y * y * y * z + std::sqrt(10.0) * y * y * z * z * z) + e_1 * (-std::sqrt(5.625) * x * x * z - std::sqrt(50.625) * y * y * z + std::sqrt(10.0) * z * z * z);
    }

#pragma omp simd aligned(pe_0, pe_1, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];

        pc_4[k] = e_0 * (0.375 * x * x * x * x * y + 0.75 * x * x * y * y * y - 3.0 * x * x * y * z * z + 0.375 * y * y * y * y * y - 3.0 * y * y * y * z * z + y * z * z * z * z) + e_1 * (1.5 * x * x * y + 1.5 * y * y * y - 6.0 * y * z * z);

        pc_5[k] = e_0 * (-std::sqrt(5.625) * x * x * x * y * z - std::sqrt(5.625) * x * y * y * y * z + std::sqrt(10.0) * x * y * z * z * z) + e_1 * (-std::sqrt(22.5) * x * y * z);

        pc_6[k] = e_0 * (-std::sqrt(0.3125) * x * x * x * x * y + std::sqrt(11.25) * x * x * y * z * z + std::sqrt(0.3125) * y * y * y * y * y - std::sqrt(11.25) * y * y * y * z * z) + e_1 * (std::sqrt(5.0) * y * y * y - std::sqrt(45.0) * y * z * z);

        pc_7[k] = e_0 * (std::sqrt(4.375) * x * x * x * y * z - std::sqrt(39.375) * x * y * y * y * z) + e_1 * (-std::sqrt(157.5) * x * y * z);
    }

#pragma omp simd aligned(pe_0, pe_1, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];

        pc_8[k] = e_0 * (std::sqrt(0.546875) * x * x * x * x * y - std::sqrt(19.6875) * x * x * y * y * y + std::sqrt(0.546875) * y * y * y * y * y) + e_1 * (-std::sqrt(78.75) * x * x * y + std::sqrt(8.75) * y * y * y);

        pc_9[k] = e_0 * (std::sqrt(8.75) * x * x * x * y * z - std::sqrt(8.75) * x * y * y * y * z);

        pc_10[k] = e_0 * (std::sqrt(39.375) * x * x * y * z * z - std::sqrt(4.375) * y * y * y * z * z) + e_1 * (std::sqrt(39.375) * x * x * y - std::sqrt(4.375) * y * y * y);

        pc_11[k] = e_0 * (-std::sqrt(1.25) * x * x * x * y * z - std::sqrt(1.25) * x * y * y * y * z + std::sqrt(45.0) * x * y * z * z * z) + e_1 * (std::sqrt(180.0) * x * y * z);
    }

#pragma omp simd aligned(pe_0, pe_1, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];

        pc_12[k] = e_0 * (-std::sqrt(5.625) * x * x * y * z * z - std::sqrt(5.625) * y * y * y * z * z + std::sqrt(10.0) * y * z * z * z * z) + e_1 * (-std::sqrt(5.625) * x * x * y - std::sqrt(5.625) * y * y * y + std::sqrt(90.0) * y * z * z);

        pc_13[k] = e_0 * (0.375 * x * x * x * x * z + 0.75 * x * x * y * y * z - 3.0 * x * x * z * z * z + 0.375 * y * y * y * y * z - 3.0 * y * y * z * z * z + z * z * z * z * z) + e_1 * (-6.0 * x * x * z - 6.0 * y * y * z + 4.0 * z * z * z);

        pc_14[k] = e_0 * (-std::sqrt(5.625) * x * x * x * z * z - std::sqrt(5.625) * x * y * y * z * z + std::sqrt(10.0) * x * z * z * z * z) + e_1 * (-std::sqrt(5.625) * x * x * x - std::sqrt(5.625) * x * y * y + std::sqrt(90.0) * x * z * z);

        pc_15[k] = e_0 * (-std::sqrt(0.3125) * x * x * x * x * z + std::sqrt(11.25) * x * x * z * z * z + std::sqrt(0.3125) * y * y * y * y * z - std::sqrt(11.25) * y * y * z * z * z) + e_1 * (std::sqrt(45.0) * x * x * z - std::sqrt(45.0) * y * y * z);
    }

#pragma omp simd aligned(pe_0, pe_1, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];

        pc_16[k] = e_0 * (std::sqrt(4.375) * x * x * x * z * z - std::sqrt(39.375) * x * y * y * z * z) + e_1 * (std::sqrt(4.375) * x * x * x - std::sqrt(39.375) * x * y * y);

        pc_17[k] = e_0 * (std::sqrt(0.546875) * x * x * x * x * z - std::sqrt(19.6875) * x * x * y * y * z + std::sqrt(0.546875) * y * y * y * y * z);

        pc_18[k] = e_0 * (std::sqrt(8.75) * x * x * x * x * y - std::sqrt(8.75) * x * x * y * y * y) + e_1 * (std::sqrt(78.75) * x * x * y - std::sqrt(8.75) * y * y * y);

        pc_19[k] = e_0 * (std::sqrt(39.375) * x * x * x * y * z - std::sqrt(4.375) * x * y * y * y * z) + e_1 * (std::sqrt(157.5) * x * y * z);
    }

#pragma omp simd aligned(pe_0, pe_1, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];

        pc_20[k] = e_0 * (-std::sqrt(1.25) * x * x * x * x * y - std::sqrt(1.25) * x * x * y * y * y + std::sqrt(45.0) * x * x * y * z * z) + e_1 * (-std::sqrt(11.25) * x * x * y - std::sqrt(1.25) * y * y * y + std::sqrt(45.0) * y * z * z);

        pc_21[k] = e_0 * (-std::sqrt(5.625) * x * x * x * y * z - std::sqrt(5.625) * x * y * y * y * z + std::sqrt(10.0) * x * y * z * z * z) + e_1 * (-std::sqrt(22.5) * x * y * z);

        pc_22[k] = e_0 * (0.375 * x * x * x * x * x + 0.75 * x * x * x * y * y - 3.0 * x * x * x * z * z + 0.375 * x * y * y * y * y - 3.0 * x * y * y * z * z + x * z * z * z * z) + e_1 * (1.5 * x * x * x + 1.5 * x * y * y - 6.0 * x * z * z);

        pc_23[k] = e_0 * (-std::sqrt(5.625) * x * x * x * x * z - std::sqrt(5.625) * x * x * y * y * z + std::sqrt(10.0) * x * x * z * z * z) + e_1 * (-std::sqrt(50.625) * x * x * z - std::sqrt(5.625) * y * y * z + std::sqrt(10.0) * z * z * z);
    }

#pragma omp simd aligned(pe_0, pe_1, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];

        pc_24[k] = e_0 * (-std::sqrt(0.3125) * x * x * x * x * x + std::sqrt(11.25) * x * x * x * z * z + std::sqrt(0.3125) * x * y * y * y * y - std::sqrt(11.25) * x * y * y * z * z) + e_1 * (-std::sqrt(5.0) * x * x * x + std::sqrt(45.0) * x * z * z);

        pc_25[k] = e_0 * (std::sqrt(4.375) * x * x * x * x * z - std::sqrt(39.375) * x * x * y * y * z) + e_1 * (std::sqrt(39.375) * x * x * z - std::sqrt(39.375) * y * y * z);

        pc_26[k] = e_0 * (std::sqrt(0.546875) * x * x * x * x * x - std::sqrt(19.6875) * x * x * x * y * y + std::sqrt(0.546875) * x * y * y * y * y) + e_1 * (std::sqrt(8.75) * x * x * x - std::sqrt(78.75) * x * y * y);
    }

    // NOTE: the atom pairs beyond the reach of every pair of primitives have no
    // contribution and are set to zero.

    for (size_t m = 0; m < 27; m++)
    {
        auto *pv = values + m * nvalues;

        std::fill(pv + nmax, pv + nvalues, 0.0);
    }
}

}  // namespace simdovl
