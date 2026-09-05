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



#include "SimdOverlapRecDI.hpp"

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
compute_di_overlap(double               *values,
                   const size_t          nvalues,
                   const CBasisFunction &bra,
                   const CBasisFunction &ket,
                   const CSimdMatrix    &coordinates,
                   const double          threshold) -> void
{
    if ((bra.get_angular_momentum() != 2) || (ket.get_angular_momentum() != 6))
    {
        errors::assertMsgCritical(
            false, std::string("SimdOverlapRecDI.compute_di_overlap: Basis functions must be of angular momenta two and six"));
    }

    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("SimdOverlapRecDI.compute_di_overlap: Number of values exceeds number of atom pairs"));
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

    auto buffer = simdfunc::make_primitive_buffer(dimensions, 3);

    if (buffer.number_of_columns() == 0)
    {
        std::fill(values, values + 65 * nvalues, 0.0);

        return;
    }

    const auto nmax = buffer.number_of_columns();

    auto *pe_0 = buffer.data(0);
    auto *pe_1 = buffer.data(1);
    auto *pe_2 = buffer.data(2);

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

        const auto f_0 = fbase * fal * fal * fbe * fbe * fbe * fbe * fbe * fbe;

        const auto f_1 = fbase * fal * fbe * fbe * fbe * fbe * fbe * fh;

        const auto f_2 = fbase * fbe * fbe * fbe * fbe * fh * fh;

        // NOTE: the exponential depends on the pair of primitives alone, so it is
        // evaluated once and shared by the prefactors of all terms.

#pragma omp simd aligned(pe_0, pe_1, pe_2, ab_2 : simd::cache_line_size())
        for (size_t k = 0; k < ncols; k++)
        {
            const auto fss = std::exp(-fmu * ab_2[k]);

            pe_0[k] += f_0 * fss;
            pe_1[k] += f_1 * fss;
            pe_2[k] += f_2 * fss;
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
    auto *pc_27 = values + 27 * nvalues;
    auto *pc_28 = values + 28 * nvalues;
    auto *pc_29 = values + 29 * nvalues;
    auto *pc_30 = values + 30 * nvalues;
    auto *pc_31 = values + 31 * nvalues;
    auto *pc_32 = values + 32 * nvalues;
    auto *pc_33 = values + 33 * nvalues;
    auto *pc_34 = values + 34 * nvalues;
    auto *pc_35 = values + 35 * nvalues;
    auto *pc_36 = values + 36 * nvalues;
    auto *pc_37 = values + 37 * nvalues;
    auto *pc_38 = values + 38 * nvalues;
    auto *pc_39 = values + 39 * nvalues;
    auto *pc_40 = values + 40 * nvalues;
    auto *pc_41 = values + 41 * nvalues;
    auto *pc_42 = values + 42 * nvalues;
    auto *pc_43 = values + 43 * nvalues;
    auto *pc_44 = values + 44 * nvalues;
    auto *pc_45 = values + 45 * nvalues;
    auto *pc_46 = values + 46 * nvalues;
    auto *pc_47 = values + 47 * nvalues;
    auto *pc_48 = values + 48 * nvalues;
    auto *pc_49 = values + 49 * nvalues;
    auto *pc_50 = values + 50 * nvalues;
    auto *pc_51 = values + 51 * nvalues;
    auto *pc_52 = values + 52 * nvalues;
    auto *pc_53 = values + 53 * nvalues;
    auto *pc_54 = values + 54 * nvalues;
    auto *pc_55 = values + 55 * nvalues;
    auto *pc_56 = values + 56 * nvalues;
    auto *pc_57 = values + 57 * nvalues;
    auto *pc_58 = values + 58 * nvalues;
    auto *pc_59 = values + 59 * nvalues;
    auto *pc_60 = values + 60 * nvalues;
    auto *pc_61 = values + 61 * nvalues;
    auto *pc_62 = values + 62 * nvalues;
    auto *pc_63 = values + 63 * nvalues;
    auto *pc_64 = values + 64 * nvalues;

    // NOTE: the components are formed in 17 loops, as the vectorizer runs out
    // of registers with all of them in one. Only the prefactors and the vector
    // between the atoms are loaded by more than one loop.

#pragma omp simd aligned(pe_0, pe_1, pe_2, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];

        pc_0[k] = e_0 * (std::sqrt(48.7265625) * x * x * x * x * x * x * y * y - std::sqrt(541.40625) * x * x * x * x * y * y * y * y + std::sqrt(48.7265625) * x * x * y * y * y * y * y * y) + e_1 * (std::sqrt(48.7265625) * x * x * x * x * x * x - std::sqrt(1218.1640625) * x * x * x * x * y * y - std::sqrt(1218.1640625) * x * x * y * y * y * y + std::sqrt(48.7265625) * y * y * y * y * y * y) + e_2 * (std::sqrt(1218.1640625) * x * x * x * x - std::sqrt(43853.90625) * x * x * y * y + std::sqrt(1218.1640625) * y * y * y * y);

        pc_1[k] = e_0 * (std::sqrt(406.0546875) * x * x * x * x * x * y * y * z - std::sqrt(1624.21875) * x * x * x * y * y * y * y * z + std::sqrt(16.2421875) * x * y * y * y * y * y * y * z) + e_1 * (std::sqrt(406.0546875) * x * x * x * x * x * z - std::sqrt(1624.21875) * x * x * x * y * y * z - std::sqrt(3654.4921875) * x * y * y * y * y * z) + e_2 * (std::sqrt(6496.875) * x * x * x * z - std::sqrt(58471.875) * x * y * y * z);

        pc_2[k] = e_0 * (-std::sqrt(11.8125) * x * x * x * x * x * x * y * y + std::sqrt(1181.25) * x * x * x * x * y * y * z * z + std::sqrt(11.8125) * x * x * y * y * y * y * y * y - std::sqrt(1181.25) * x * x * y * y * y * y * z * z) + e_1 * (-std::sqrt(11.8125) * x * x * x * x * x * x - std::sqrt(295.3125) * x * x * x * x * y * y + std::sqrt(1181.25) * x * x * x * x * z * z + std::sqrt(295.3125) * x * x * y * y * y * y + std::sqrt(11.8125) * y * y * y * y * y * y - std::sqrt(1181.25) * y * y * y * y * z * z) + e_2 * (-std::sqrt(295.3125) * x * x * x * x + std::sqrt(10631.25) * x * x * z * z + std::sqrt(295.3125) * y * y * y * y - std::sqrt(10631.25) * y * y * z * z);

        pc_3[k] = e_0 * (-std::sqrt(199.3359375) * x * x * x * x * x * y * y * z - std::sqrt(88.59375) * x * x * x * y * y * y * y * z + std::sqrt(1417.5) * x * x * x * y * y * z * z * z + std::sqrt(22.1484375) * x * y * y * y * y * y * y * z - std::sqrt(157.5) * x * y * y * y * y * z * z * z) + e_1 * (-std::sqrt(199.3359375) * x * x * x * x * x * z - std::sqrt(7176.09375) * x * x * x * y * y * z + std::sqrt(1417.5) * x * x * x * z * z * z + std::sqrt(22.1484375) * x * y * y * y * y * z + std::sqrt(1417.5) * x * y * y * z * z * z) + e_2 * (-std::sqrt(3189.375) * x * x * x * z - std::sqrt(3189.375) * x * y * y * z + std::sqrt(5670.0) * x * z * z * z);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];

        pc_4[k] = e_0 * (std::sqrt(2.4609375) * x * x * x * x * x * x * y * y + std::sqrt(9.84375) * x * x * x * x * y * y * y * y - std::sqrt(630.0) * x * x * x * x * y * y * z * z + std::sqrt(2.4609375) * x * x * y * y * y * y * y * y - std::sqrt(630.0) * x * x * y * y * y * y * z * z + std::sqrt(630.0) * x * x * y * y * z * z * z * z) + e_1 * (std::sqrt(2.4609375) * x * x * x * x * x * x + std::sqrt(297.7734375) * x * x * x * x * y * y - std::sqrt(630.0) * x * x * x * x * z * z + std::sqrt(297.7734375) * x * x * y * y * y * y - std::sqrt(22680.0) * x * x * y * y * z * z + std::sqrt(630.0) * x * x * z * z * z * z + std::sqrt(2.4609375) * y * y * y * y * y * y - std::sqrt(630.0) * y * y * y * y * z * z + std::sqrt(630.0) * y * y * z * z * z * z) + e_2 * (std::sqrt(61.5234375) * x * x * x * x + std::sqrt(797.34375) * x * x * y * y - std::sqrt(5670.0) * x * x * z * z + std::sqrt(61.5234375) * y * y * y * y - std::sqrt(5670.0) * y * y * z * z + std::sqrt(630.0) * z * z * z * z);

        pc_5[k] = e_0 * (std::sqrt(24.609375) * x * x * x * x * x * y * y * z + std::sqrt(98.4375) * x * x * x * y * y * y * y * z - std::sqrt(393.75) * x * x * x * y * y * z * z * z + std::sqrt(24.609375) * x * y * y * y * y * y * y * z - std::sqrt(393.75) * x * y * y * y * y * z * z * z + std::sqrt(63.0) * x * y * y * z * z * z * z * z) + e_1 * (std::sqrt(24.609375) * x * x * x * x * x * z + std::sqrt(2460.9375) * x * x * x * y * y * z - std::sqrt(393.75) * x * x * x * z * z * z + std::sqrt(1993.359375) * x * y * y * y * y * z - std::sqrt(9843.75) * x * y * y * z * z * z + std::sqrt(63.0) * x * z * z * z * z * z) + e_2 * (std::sqrt(393.75) * x * x * x * z + std::sqrt(3543.75) * x * y * y * z - std::sqrt(1575.0) * x * z * z * z);

        pc_6[k] = e_0 * (-std::sqrt(0.29296875) * x * x * x * x * x * x * x * y - std::sqrt(2.63671875) * x * x * x * x * x * y * y * y + std::sqrt(94.921875) * x * x * x * x * x * y * z * z - std::sqrt(2.63671875) * x * x * x * y * y * y * y * y + std::sqrt(379.6875) * x * x * x * y * y * y * z * z - std::sqrt(168.75) * x * x * x * y * z * z * z * z - std::sqrt(0.29296875) * x * y * y * y * y * y * y * y + std::sqrt(94.921875) * x * y * y * y * y * y * z * z - std::sqrt(168.75) * x * y * y * y * z * z * z * z + std::sqrt(3.0) * x * y * z * z * z * z * z * z) + e_1 * (-std::sqrt(42.1875) * x * x * x * x * x * y - std::sqrt(168.75) * x * x * x * y * y * y + std::sqrt(6075.0) * x * x * x * y * z * z - std::sqrt(42.1875) * x * y * y * y * y * y + std::sqrt(6075.0) * x * y * y * y * z * z - std::sqrt(2700.0) * x * y * z * z * z * z) + e_2 * (-std::sqrt(168.75) * x * x * x * y - std::sqrt(168.75) * x * y * y * y + std::sqrt(6075.0) * x * y * z * z);

        pc_7[k] = e_0 * (std::sqrt(24.609375) * x * x * x * x * x * x * y * z + std::sqrt(98.4375) * x * x * x * x * y * y * y * z - std::sqrt(393.75) * x * x * x * x * y * z * z * z + std::sqrt(24.609375) * x * x * y * y * y * y * y * z - std::sqrt(393.75) * x * x * y * y * y * z * z * z + std::sqrt(63.0) * x * x * y * z * z * z * z * z) + e_1 * (std::sqrt(1993.359375) * x * x * x * x * y * z + std::sqrt(2460.9375) * x * x * y * y * y * z - std::sqrt(9843.75) * x * x * y * z * z * z + std::sqrt(24.609375) * y * y * y * y * y * z - std::sqrt(393.75) * y * y * y * z * z * z + std::sqrt(63.0) * y * z * z * z * z * z) + e_2 * (std::sqrt(3543.75) * x * x * y * z + std::sqrt(393.75) * y * y * y * z - std::sqrt(1575.0) * y * z * z * z);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];

        pc_8[k] = e_0 * (std::sqrt(0.615234375) * x * x * x * x * x * x * x * y + std::sqrt(0.615234375) * x * x * x * x * x * y * y * y - std::sqrt(157.5) * x * x * x * x * x * y * z * z - std::sqrt(0.615234375) * x * x * x * y * y * y * y * y + std::sqrt(157.5) * x * x * x * y * z * z * z * z - std::sqrt(0.615234375) * x * y * y * y * y * y * y * y + std::sqrt(157.5) * x * y * y * y * y * y * z * z - std::sqrt(157.5) * x * y * y * y * z * z * z * z) + e_1 * (std::sqrt(39.375) * x * x * x * x * x * y - std::sqrt(2520.0) * x * x * x * y * z * z - std::sqrt(39.375) * x * y * y * y * y * y + std::sqrt(2520.0) * x * y * y * y * z * z) + e_2 * (std::sqrt(39.375) * x * x * x * y - std::sqrt(39.375) * x * y * y * y);

        pc_9[k] = e_0 * (-std::sqrt(22.1484375) * x * x * x * x * x * x * y * z + std::sqrt(88.59375) * x * x * x * x * y * y * y * z + std::sqrt(157.5) * x * x * x * x * y * z * z * z + std::sqrt(199.3359375) * x * x * y * y * y * y * y * z - std::sqrt(1417.5) * x * x * y * y * y * z * z * z) + e_1 * (-std::sqrt(22.1484375) * x * x * x * x * y * z + std::sqrt(7176.09375) * x * x * y * y * y * z - std::sqrt(1417.5) * x * x * y * z * z * z + std::sqrt(199.3359375) * y * y * y * y * y * z - std::sqrt(1417.5) * y * y * y * z * z * z) + e_2 * (std::sqrt(3189.375) * x * x * y * z + std::sqrt(3189.375) * y * y * y * z - std::sqrt(5670.0) * y * z * z * z);

        pc_10[k] = e_0 * (-std::sqrt(0.73828125) * x * x * x * x * x * x * x * y + std::sqrt(18.45703125) * x * x * x * x * x * y * y * y + std::sqrt(73.828125) * x * x * x * x * x * y * z * z + std::sqrt(18.45703125) * x * x * x * y * y * y * y * y - std::sqrt(2657.8125) * x * x * x * y * y * y * z * z - std::sqrt(0.73828125) * x * y * y * y * y * y * y * y + std::sqrt(73.828125) * x * y * y * y * y * y * z * z) + e_1 * (std::sqrt(11.8125) * x * x * x * x * x * y + std::sqrt(1181.25) * x * x * x * y * y * y - std::sqrt(4725.0) * x * x * x * y * z * z + std::sqrt(11.8125) * x * y * y * y * y * y - std::sqrt(4725.0) * x * y * y * y * z * z) + e_2 * (std::sqrt(1181.25) * x * x * x * y + std::sqrt(1181.25) * x * y * y * y - std::sqrt(42525.0) * x * y * z * z);

        pc_11[k] = e_0 * (std::sqrt(16.2421875) * x * x * x * x * x * x * y * z - std::sqrt(1624.21875) * x * x * x * x * y * y * y * z + std::sqrt(406.0546875) * x * x * y * y * y * y * y * z) + e_1 * (-std::sqrt(3654.4921875) * x * x * x * x * y * z - std::sqrt(1624.21875) * x * x * y * y * y * z + std::sqrt(406.0546875) * y * y * y * y * y * z) + e_2 * (-std::sqrt(58471.875) * x * x * y * z + std::sqrt(6496.875) * y * y * y * z);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];

        pc_12[k] = e_0 * (std::sqrt(1.353515625) * x * x * x * x * x * x * x * y - std::sqrt(304.541015625) * x * x * x * x * x * y * y * y + std::sqrt(304.541015625) * x * x * x * y * y * y * y * y - std::sqrt(1.353515625) * x * y * y * y * y * y * y * y) + e_1 * (-std::sqrt(779.625) * x * x * x * x * x * y + std::sqrt(779.625) * x * y * y * y * y * y) + e_2 * (-std::sqrt(19490.625) * x * x * x * y + std::sqrt(19490.625) * x * y * y * y);

        pc_13[k] = e_0 * (std::sqrt(48.7265625) * x * x * x * x * x * y * y * z - std::sqrt(541.40625) * x * x * x * y * y * y * y * z + std::sqrt(48.7265625) * x * y * y * y * y * y * y * z) + e_1 * (std::sqrt(48.7265625) * x * x * x * x * x * z - std::sqrt(4872.65625) * x * x * x * y * y * z + std::sqrt(1218.1640625) * x * y * y * y * y * z);

        pc_14[k] = e_0 * (std::sqrt(406.0546875) * x * x * x * x * y * y * z * z - std::sqrt(1624.21875) * x * x * y * y * y * y * z * z + std::sqrt(16.2421875) * y * y * y * y * y * y * z * z) + e_1 * (std::sqrt(406.0546875) * x * x * x * x * y * y + std::sqrt(406.0546875) * x * x * x * x * z * z - std::sqrt(1624.21875) * x * x * y * y * y * y - std::sqrt(14617.96875) * x * x * y * y * z * z + std::sqrt(16.2421875) * y * y * y * y * y * y + std::sqrt(406.0546875) * y * y * y * y * z * z) + e_2 * (std::sqrt(406.0546875) * x * x * x * x - std::sqrt(14617.96875) * x * x * y * y + std::sqrt(406.0546875) * y * y * y * y);

        pc_15[k] = e_0 * (-std::sqrt(11.8125) * x * x * x * x * x * y * y * z + std::sqrt(1181.25) * x * x * x * y * y * z * z * z + std::sqrt(11.8125) * x * y * y * y * y * y * y * z - std::sqrt(1181.25) * x * y * y * y * y * z * z * z) + e_1 * (-std::sqrt(11.8125) * x * x * x * x * x * z + std::sqrt(4725.0) * x * x * x * y * y * z + std::sqrt(1181.25) * x * x * x * z * z * z - std::sqrt(2657.8125) * x * y * y * y * y * z - std::sqrt(10631.25) * x * y * y * z * z * z) + e_2 * (std::sqrt(4725.0) * x * x * x * z - std::sqrt(42525.0) * x * y * y * z);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];

        pc_16[k] = e_0 * (-std::sqrt(199.3359375) * x * x * x * x * y * y * z * z - std::sqrt(88.59375) * x * x * y * y * y * y * z * z + std::sqrt(1417.5) * x * x * y * y * z * z * z * z + std::sqrt(22.1484375) * y * y * y * y * y * y * z * z - std::sqrt(157.5) * y * y * y * y * z * z * z * z) + e_1 * (-std::sqrt(199.3359375) * x * x * x * x * y * y - std::sqrt(199.3359375) * x * x * x * x * z * z - std::sqrt(88.59375) * x * x * y * y * y * y + std::sqrt(7176.09375) * x * x * y * y * z * z + std::sqrt(1417.5) * x * x * z * z * z * z + std::sqrt(22.1484375) * y * y * y * y * y * y - std::sqrt(199.3359375) * y * y * y * y * z * z - std::sqrt(1417.5) * y * y * z * z * z * z) + e_2 * (-std::sqrt(199.3359375) * x * x * x * x - std::sqrt(797.34375) * x * x * y * y + std::sqrt(12757.5) * x * x * z * z + std::sqrt(553.7109375) * y * y * y * y - std::sqrt(12757.5) * y * y * z * z);

        pc_17[k] = e_0 * (std::sqrt(2.4609375) * x * x * x * x * x * y * y * z + std::sqrt(9.84375) * x * x * x * y * y * y * y * z - std::sqrt(630.0) * x * x * x * y * y * z * z * z + std::sqrt(2.4609375) * x * y * y * y * y * y * y * z - std::sqrt(630.0) * x * y * y * y * y * z * z * z + std::sqrt(630.0) * x * y * y * z * z * z * z * z) + e_1 * (std::sqrt(2.4609375) * x * x * x * x * x * z - std::sqrt(1663.59375) * x * x * x * y * y * z - std::sqrt(630.0) * x * x * x * z * z * z - std::sqrt(1794.0234375) * x * y * y * y * y * z + std::sqrt(630.0) * x * y * y * z * z * z + std::sqrt(630.0) * x * z * z * z * z * z) + e_2 * (-std::sqrt(2520.0) * x * x * x * z - std::sqrt(22680.0) * x * y * y * z + std::sqrt(10080.0) * x * z * z * z);

        pc_18[k] = e_0 * (std::sqrt(24.609375) * x * x * x * x * y * y * z * z + std::sqrt(98.4375) * x * x * y * y * y * y * z * z - std::sqrt(393.75) * x * x * y * y * z * z * z * z + std::sqrt(24.609375) * y * y * y * y * y * y * z * z - std::sqrt(393.75) * y * y * y * y * z * z * z * z + std::sqrt(63.0) * y * y * z * z * z * z * z * z) + e_1 * (std::sqrt(24.609375) * x * x * x * x * y * y + std::sqrt(24.609375) * x * x * x * x * z * z + std::sqrt(98.4375) * x * x * y * y * y * y - std::sqrt(885.9375) * x * x * y * y * z * z - std::sqrt(393.75) * x * x * z * z * z * z + std::sqrt(24.609375) * y * y * y * y * y * y - std::sqrt(1205.859375) * y * y * y * y * z * z - std::sqrt(393.75) * y * y * z * z * z * z + std::sqrt(63.0) * z * z * z * z * z * z) + e_2 * (std::sqrt(24.609375) * x * x * x * x + std::sqrt(885.9375) * x * x * y * y - std::sqrt(3543.75) * x * x * z * z + std::sqrt(615.234375) * y * y * y * y - std::sqrt(31893.75) * y * y * z * z + std::sqrt(1575.0) * z * z * z * z);

        pc_19[k] = e_0 * (-std::sqrt(0.29296875) * x * x * x * x * x * x * y * z - std::sqrt(2.63671875) * x * x * x * x * y * y * y * z + std::sqrt(94.921875) * x * x * x * x * y * z * z * z - std::sqrt(2.63671875) * x * x * y * y * y * y * y * z + std::sqrt(379.6875) * x * x * y * y * y * z * z * z - std::sqrt(168.75) * x * x * y * z * z * z * z * z - std::sqrt(0.29296875) * y * y * y * y * y * y * y * z + std::sqrt(94.921875) * y * y * y * y * y * z * z * z - std::sqrt(168.75) * y * y * y * z * z * z * z * z + std::sqrt(3.0) * y * z * z * z * z * z * z * z) + e_1 * (std::sqrt(263.671875) * x * x * x * x * y * z + std::sqrt(1054.6875) * x * x * y * y * y * z - std::sqrt(168.75) * x * x * y * z * z * z + std::sqrt(263.671875) * y * y * y * y * y * z - std::sqrt(168.75) * y * y * y * z * z * z - std::sqrt(243.0) * y * z * z * z * z * z) + e_2 * (std::sqrt(6075.0) * x * x * y * z + std::sqrt(6075.0) * y * y * y * z - std::sqrt(10800.0) * y * z * z * z);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];

        pc_20[k] = e_0 * (std::sqrt(24.609375) * x * x * x * x * x * y * z * z + std::sqrt(98.4375) * x * x * x * y * y * y * z * z - std::sqrt(393.75) * x * x * x * y * z * z * z * z + std::sqrt(24.609375) * x * y * y * y * y * y * z * z - std::sqrt(393.75) * x * y * y * y * z * z * z * z + std::sqrt(63.0) * x * y * z * z * z * z * z * z) + e_1 * (std::sqrt(24.609375) * x * x * x * x * x * y + std::sqrt(98.4375) * x * x * x * y * y * y - std::sqrt(1575.0) * x * x * x * y * z * z + std::sqrt(24.609375) * x * y * y * y * y * y - std::sqrt(1575.0) * x * y * y * y * z * z) + e_2 * (std::sqrt(393.75) * x * x * x * y + std::sqrt(393.75) * x * y * y * y - std::sqrt(14175.0) * x * y * z * z);

        pc_21[k] = e_0 * (std::sqrt(0.615234375) * x * x * x * x * x * x * y * z + std::sqrt(0.615234375) * x * x * x * x * y * y * y * z - std::sqrt(157.5) * x * x * x * x * y * z * z * z - std::sqrt(0.615234375) * x * x * y * y * y * y * y * z + std::sqrt(157.5) * x * x * y * z * z * z * z * z - std::sqrt(0.615234375) * y * y * y * y * y * y * y * z + std::sqrt(157.5) * y * y * y * y * y * z * z * z - std::sqrt(157.5) * y * y * y * z * z * z * z * z) + e_1 * (-std::sqrt(553.7109375) * x * x * x * x * y * z - std::sqrt(9.84375) * x * x * y * y * y * z + std::sqrt(2520.0) * x * x * y * z * z * z + std::sqrt(415.8984375) * y * y * y * y * y * z - std::sqrt(630.0) * y * z * z * z * z * z) + e_2 * (std::sqrt(10080.0) * y * y * y * z - std::sqrt(10080.0) * y * z * z * z);

        pc_22[k] = e_0 * (-std::sqrt(22.1484375) * x * x * x * x * x * y * z * z + std::sqrt(88.59375) * x * x * x * y * y * y * z * z + std::sqrt(157.5) * x * x * x * y * z * z * z * z + std::sqrt(199.3359375) * x * y * y * y * y * y * z * z - std::sqrt(1417.5) * x * y * y * y * z * z * z * z) + e_1 * (-std::sqrt(22.1484375) * x * x * x * x * x * y + std::sqrt(88.59375) * x * x * x * y * y * y + std::sqrt(3189.375) * x * x * x * y * z * z + std::sqrt(199.3359375) * x * y * y * y * y * y - std::sqrt(3189.375) * x * y * y * y * z * z - std::sqrt(5670.0) * x * y * z * z * z * z) + e_2 * (std::sqrt(354.375) * x * x * x * y + std::sqrt(3189.375) * x * y * y * y - std::sqrt(51030.0) * x * y * z * z);

        pc_23[k] = e_0 * (-std::sqrt(0.73828125) * x * x * x * x * x * x * y * z + std::sqrt(18.45703125) * x * x * x * x * y * y * y * z + std::sqrt(73.828125) * x * x * x * x * y * z * z * z + std::sqrt(18.45703125) * x * x * y * y * y * y * y * z - std::sqrt(2657.8125) * x * x * y * y * y * z * z * z - std::sqrt(0.73828125) * y * y * y * y * y * y * y * z + std::sqrt(73.828125) * y * y * y * y * y * z * z * z) + e_1 * (std::sqrt(664.453125) * x * x * x * x * y * z - std::sqrt(7382.8125) * x * x * y * y * y * z - std::sqrt(10631.25) * x * x * y * z * z * z + std::sqrt(144.703125) * y * y * y * y * y * z + std::sqrt(1181.25) * y * y * y * z * z * z) + e_2 * (-std::sqrt(42525.0) * x * x * y * z + std::sqrt(4725.0) * y * y * y * z);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];

        pc_24[k] = e_0 * (std::sqrt(16.2421875) * x * x * x * x * x * y * z * z - std::sqrt(1624.21875) * x * x * x * y * y * y * z * z + std::sqrt(406.0546875) * x * y * y * y * y * y * z * z) + e_1 * (std::sqrt(16.2421875) * x * x * x * x * x * y - std::sqrt(1624.21875) * x * x * x * y * y * y - std::sqrt(6496.875) * x * x * x * y * z * z + std::sqrt(406.0546875) * x * y * y * y * y * y + std::sqrt(6496.875) * x * y * y * y * z * z) + e_2 * (-std::sqrt(6496.875) * x * x * x * y + std::sqrt(6496.875) * x * y * y * y);

        pc_25[k] = e_0 * (std::sqrt(1.353515625) * x * x * x * x * x * x * y * z - std::sqrt(304.541015625) * x * x * x * x * y * y * y * z + std::sqrt(304.541015625) * x * x * y * y * y * y * y * z - std::sqrt(1.353515625) * y * y * y * y * y * y * y * z) + e_1 * (-std::sqrt(1218.1640625) * x * x * x * x * y * z + std::sqrt(4872.65625) * x * x * y * y * y * z - std::sqrt(48.7265625) * y * y * y * y * y * z);

        pc_26[k] = e_0 * (-std::sqrt(4.060546875) * x * x * x * x * x * x * x * y + std::sqrt(22.107421875) * x * x * x * x * x * y * y * y + std::sqrt(16.2421875) * x * x * x * x * x * y * z * z + std::sqrt(22.107421875) * x * x * x * y * y * y * y * y - std::sqrt(180.46875) * x * x * x * y * y * y * z * z - std::sqrt(4.060546875) * x * y * y * y * y * y * y * y + std::sqrt(16.2421875) * x * y * y * y * y * y * z * z) + e_1 * (-std::sqrt(584.71875) * x * x * x * x * x * y + std::sqrt(6496.875) * x * x * x * y * y * y - std::sqrt(584.71875) * x * y * y * y * y * y);

        pc_27[k] = e_0 * (-std::sqrt(33.837890625) * x * x * x * x * x * x * y * z + std::sqrt(33.837890625) * x * x * x * x * y * y * y * z + std::sqrt(135.3515625) * x * x * x * x * y * z * z * z + std::sqrt(109.634765625) * x * x * y * y * y * y * y * z - std::sqrt(541.40625) * x * x * y * y * y * z * z * z - std::sqrt(1.353515625) * y * y * y * y * y * y * y * z + std::sqrt(5.4140625) * y * y * y * y * y * z * z * z) + e_1 * (-std::sqrt(1218.1640625) * x * x * x * x * y * z + std::sqrt(4872.65625) * x * x * y * y * y * z - std::sqrt(48.7265625) * y * y * y * y * y * z);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];

        pc_28[k] = e_0 * (std::sqrt(0.984375) * x * x * x * x * x * x * x * y + std::sqrt(0.984375) * x * x * x * x * x * y * y * y - std::sqrt(141.75) * x * x * x * x * x * y * z * z - std::sqrt(0.984375) * x * x * x * y * y * y * y * y + std::sqrt(393.75) * x * x * x * y * z * z * z * z - std::sqrt(0.984375) * x * y * y * y * y * y * y * y + std::sqrt(141.75) * x * y * y * y * y * y * z * z - std::sqrt(393.75) * x * y * y * y * z * z * z * z) + e_1 * (std::sqrt(141.75) * x * x * x * x * x * y - std::sqrt(141.75) * x * y * y * y * y * y) + e_2 * (std::sqrt(3543.75) * x * x * x * y - std::sqrt(3543.75) * x * y * y * y);

        pc_29[k] = e_0 * (std::sqrt(16.611328125) * x * x * x * x * x * x * y * z + std::sqrt(46.142578125) * x * x * x * x * y * y * y * z - std::sqrt(361.7578125) * x * x * x * x * y * z * z * z + std::sqrt(1.845703125) * x * x * y * y * y * y * y * z - std::sqrt(160.78125) * x * x * y * y * y * z * z * z + std::sqrt(472.5) * x * x * y * z * z * z * z * z - std::sqrt(1.845703125) * y * y * y * y * y * y * y * z + std::sqrt(40.1953125) * y * y * y * y * y * z * z * z - std::sqrt(52.5) * y * y * y * z * z * z * z * z) + e_1 * (std::sqrt(598.0078125) * x * x * x * x * y * z + std::sqrt(265.78125) * x * x * y * y * y * z + std::sqrt(4252.5) * x * x * y * z * z * z - std::sqrt(66.4453125) * y * y * y * y * y * z - std::sqrt(472.5) * y * y * y * z * z * z) + e_2 * (std::sqrt(38272.5) * x * x * y * z - std::sqrt(4252.5) * y * y * y * z);

        pc_30[k] = e_0 * (-std::sqrt(0.205078125) * x * x * x * x * x * x * x * y - std::sqrt(1.845703125) * x * x * x * x * x * y * y * y + std::sqrt(66.4453125) * x * x * x * x * x * y * z * z - std::sqrt(1.845703125) * x * x * x * y * y * y * y * y + std::sqrt(265.78125) * x * x * x * y * y * y * z * z - std::sqrt(472.5) * x * x * x * y * z * z * z * z - std::sqrt(0.205078125) * x * y * y * y * y * y * y * y + std::sqrt(66.4453125) * x * y * y * y * y * y * z * z - std::sqrt(472.5) * x * y * y * y * z * z * z * z + std::sqrt(210.0) * x * y * z * z * z * z * z * z) + e_1 * (-std::sqrt(29.53125) * x * x * x * x * x * y - std::sqrt(118.125) * x * x * x * y * y * y - std::sqrt(29.53125) * x * y * y * y * y * y + std::sqrt(7560.0) * x * y * z * z * z * z) + e_2 * (-std::sqrt(1890.0) * x * x * x * y - std::sqrt(1890.0) * x * y * y * y + std::sqrt(68040.0) * x * y * z * z);

        pc_31[k] = e_0 * (-std::sqrt(2.05078125) * x * x * x * x * x * x * y * z - std::sqrt(18.45703125) * x * x * x * x * y * y * y * z + std::sqrt(73.828125) * x * x * x * x * y * z * z * z - std::sqrt(18.45703125) * x * x * y * y * y * y * y * z + std::sqrt(295.3125) * x * x * y * y * y * z * z * z - std::sqrt(189.0) * x * x * y * z * z * z * z * z - std::sqrt(2.05078125) * y * y * y * y * y * y * y * z + std::sqrt(73.828125) * y * y * y * y * y * z * z * z - std::sqrt(189.0) * y * y * y * z * z * z * z * z + std::sqrt(21.0) * y * z * z * z * z * z * z * z) + e_1 * (-std::sqrt(73.828125) * x * x * x * x * y * z - std::sqrt(295.3125) * x * x * y * y * y * z - std::sqrt(1181.25) * x * x * y * z * z * z - std::sqrt(73.828125) * y * y * y * y * y * z - std::sqrt(1181.25) * y * y * y * z * z * z + std::sqrt(1701.0) * y * z * z * z * z * z) + e_2 * (-std::sqrt(10631.25) * x * x * y * z - std::sqrt(10631.25) * y * y * y * z + std::sqrt(18900.0) * y * z * z * z);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];

        pc_32[k] = e_0 * (0.15625 * x * x * x * x * x * x * x * x + 0.625 * x * x * x * x * x * x * y * y - 3.125 * x * x * x * x * x * x * z * z + 0.9375 * x * x * x * x * y * y * y * y - 9.375 * x * x * x * x * y * y * z * z + 9.375 * x * x * x * x * z * z * z * z + 0.625 * x * x * y * y * y * y * y * y - 9.375 * x * x * y * y * y * y * z * z + 18.75 * x * x * y * y * z * z * z * z - 8.0 * x * x * z * z * z * z * z * z + 0.15625 * y * y * y * y * y * y * y * y - 3.125 * y * y * y * y * y * y * z * z + 9.375 * y * y * y * y * z * z * z * z - 8.0 * y * y * z * z * z * z * z * z + z * z * z * z * z * z * z * z) + e_1 * (1.875 * x * x * x * x * x * x + 5.625 * x * x * x * x * y * y + 5.625 * x * x * y * y * y * y - 45.0 * x * x * z * z * z * z + 1.875 * y * y * y * y * y * y - 45.0 * y * y * z * z * z * z + 12.0 * z * z * z * z * z * z) + e_2 * (16.875 * x * x * x * x + 33.75 * x * x * y * y - 135.0 * x * x * z * z + 16.875 * y * y * y * y - 135.0 * y * y * z * z + 45.0 * z * z * z * z);

        pc_33[k] = e_0 * (-std::sqrt(2.05078125) * x * x * x * x * x * x * x * z - std::sqrt(18.45703125) * x * x * x * x * x * y * y * z + std::sqrt(73.828125) * x * x * x * x * x * z * z * z - std::sqrt(18.45703125) * x * x * x * y * y * y * y * z + std::sqrt(295.3125) * x * x * x * y * y * z * z * z - std::sqrt(189.0) * x * x * x * z * z * z * z * z - std::sqrt(2.05078125) * x * y * y * y * y * y * y * z + std::sqrt(73.828125) * x * y * y * y * y * z * z * z - std::sqrt(189.0) * x * y * y * z * z * z * z * z + std::sqrt(21.0) * x * z * z * z * z * z * z * z) + e_1 * (-std::sqrt(73.828125) * x * x * x * x * x * z - std::sqrt(295.3125) * x * x * x * y * y * z - std::sqrt(1181.25) * x * x * x * z * z * z - std::sqrt(73.828125) * x * y * y * y * y * z - std::sqrt(1181.25) * x * y * y * z * z * z + std::sqrt(1701.0) * x * z * z * z * z * z) + e_2 * (-std::sqrt(10631.25) * x * x * x * z - std::sqrt(10631.25) * x * y * y * z + std::sqrt(18900.0) * x * z * z * z);

        pc_34[k] = e_0 * (-std::sqrt(0.05126953125) * x * x * x * x * x * x * x * x - std::sqrt(0.205078125) * x * x * x * x * x * x * y * y + std::sqrt(16.611328125) * x * x * x * x * x * x * z * z + std::sqrt(16.611328125) * x * x * x * x * y * y * z * z - std::sqrt(118.125) * x * x * x * x * z * z * z * z + std::sqrt(0.205078125) * x * x * y * y * y * y * y * y - std::sqrt(16.611328125) * x * x * y * y * y * y * z * z + std::sqrt(52.5) * x * x * z * z * z * z * z * z + std::sqrt(0.05126953125) * y * y * y * y * y * y * y * y - std::sqrt(16.611328125) * y * y * y * y * y * y * z * z + std::sqrt(118.125) * y * y * y * y * z * z * z * z - std::sqrt(52.5) * y * y * z * z * z * z * z * z) + e_1 * (-std::sqrt(7.3828125) * x * x * x * x * x * x - std::sqrt(7.3828125) * x * x * x * x * y * y + std::sqrt(7.3828125) * x * x * y * y * y * y + std::sqrt(1890.0) * x * x * z * z * z * z + std::sqrt(7.3828125) * y * y * y * y * y * y - std::sqrt(1890.0) * y * y * z * z * z * z) + e_2 * (-std::sqrt(472.5) * x * x * x * x + std::sqrt(17010.0) * x * x * z * z + std::sqrt(472.5) * y * y * y * y - std::sqrt(17010.0) * y * y * z * z);

        pc_35[k] = e_0 * (std::sqrt(1.845703125) * x * x * x * x * x * x * x * z - std::sqrt(1.845703125) * x * x * x * x * x * y * y * z - std::sqrt(40.1953125) * x * x * x * x * x * z * z * z - std::sqrt(46.142578125) * x * x * x * y * y * y * y * z + std::sqrt(160.78125) * x * x * x * y * y * z * z * z + std::sqrt(52.5) * x * x * x * z * z * z * z * z - std::sqrt(16.611328125) * x * y * y * y * y * y * y * z + std::sqrt(361.7578125) * x * y * y * y * y * z * z * z - std::sqrt(472.5) * x * y * y * z * z * z * z * z) + e_1 * (std::sqrt(66.4453125) * x * x * x * x * x * z - std::sqrt(265.78125) * x * x * x * y * y * z + std::sqrt(472.5) * x * x * x * z * z * z - std::sqrt(598.0078125) * x * y * y * y * y * z - std::sqrt(4252.5) * x * y * y * z * z * z) + e_2 * (std::sqrt(4252.5) * x * x * x * z - std::sqrt(38272.5) * x * y * y * z);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];

        pc_36[k] = e_0 * (std::sqrt(0.0615234375) * x * x * x * x * x * x * x * x - std::sqrt(0.984375) * x * x * x * x * x * x * y * y - std::sqrt(8.859375) * x * x * x * x * x * x * z * z - std::sqrt(6.15234375) * x * x * x * x * y * y * y * y + std::sqrt(221.484375) * x * x * x * x * y * y * z * z + std::sqrt(24.609375) * x * x * x * x * z * z * z * z - std::sqrt(0.984375) * x * x * y * y * y * y * y * y + std::sqrt(221.484375) * x * x * y * y * y * y * z * z - std::sqrt(885.9375) * x * x * y * y * z * z * z * z + std::sqrt(0.0615234375) * y * y * y * y * y * y * y * y - std::sqrt(8.859375) * y * y * y * y * y * y * z * z + std::sqrt(24.609375) * y * y * y * y * z * z * z * z) + e_1 * (std::sqrt(8.859375) * x * x * x * x * x * x - std::sqrt(221.484375) * x * x * x * x * y * y - std::sqrt(221.484375) * x * x * y * y * y * y + std::sqrt(8.859375) * y * y * y * y * y * y) + e_2 * (std::sqrt(221.484375) * x * x * x * x - std::sqrt(7973.4375) * x * x * y * y + std::sqrt(221.484375) * y * y * y * y);

        pc_37[k] = e_0 * (-std::sqrt(1.353515625) * x * x * x * x * x * x * x * z + std::sqrt(109.634765625) * x * x * x * x * x * y * y * z + std::sqrt(5.4140625) * x * x * x * x * x * z * z * z + std::sqrt(33.837890625) * x * x * x * y * y * y * y * z - std::sqrt(541.40625) * x * x * x * y * y * z * z * z - std::sqrt(33.837890625) * x * y * y * y * y * y * y * z + std::sqrt(135.3515625) * x * y * y * y * y * z * z * z) + e_1 * (-std::sqrt(48.7265625) * x * x * x * x * x * z + std::sqrt(4872.65625) * x * x * x * y * y * z - std::sqrt(1218.1640625) * x * y * y * y * y * z);

        pc_38[k] = e_0 * (-std::sqrt(0.11279296875) * x * x * x * x * x * x * x * x + std::sqrt(22.107421875) * x * x * x * x * x * x * y * y + std::sqrt(0.451171875) * x * x * x * x * x * x * z * z - std::sqrt(101.513671875) * x * x * x * x * y * y * z * z - std::sqrt(22.107421875) * x * x * y * y * y * y * y * y + std::sqrt(101.513671875) * x * x * y * y * y * y * z * z + std::sqrt(0.11279296875) * y * y * y * y * y * y * y * y - std::sqrt(0.451171875) * y * y * y * y * y * y * z * z) + e_1 * (-std::sqrt(16.2421875) * x * x * x * x * x * x + std::sqrt(3654.4921875) * x * x * x * x * y * y - std::sqrt(3654.4921875) * x * x * y * y * y * y + std::sqrt(16.2421875) * y * y * y * y * y * y);

        pc_39[k] = e_0 * (std::sqrt(48.7265625) * x * x * x * x * x * x * y * z - std::sqrt(541.40625) * x * x * x * x * y * y * y * z + std::sqrt(48.7265625) * x * x * y * y * y * y * y * z) + e_1 * (std::sqrt(1218.1640625) * x * x * x * x * y * z - std::sqrt(4872.65625) * x * x * y * y * y * z + std::sqrt(48.7265625) * y * y * y * y * y * z);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];

        pc_40[k] = e_0 * (std::sqrt(406.0546875) * x * x * x * x * x * y * z * z - std::sqrt(1624.21875) * x * x * x * y * y * y * z * z + std::sqrt(16.2421875) * x * y * y * y * y * y * z * z) + e_1 * (std::sqrt(406.0546875) * x * x * x * x * x * y - std::sqrt(1624.21875) * x * x * x * y * y * y + std::sqrt(6496.875) * x * x * x * y * z * z + std::sqrt(16.2421875) * x * y * y * y * y * y - std::sqrt(6496.875) * x * y * y * y * z * z) + e_2 * (std::sqrt(6496.875) * x * x * x * y - std::sqrt(6496.875) * x * y * y * y);

        pc_41[k] = e_0 * (-std::sqrt(11.8125) * x * x * x * x * x * x * y * z + std::sqrt(1181.25) * x * x * x * x * y * z * z * z + std::sqrt(11.8125) * x * x * y * y * y * y * y * z - std::sqrt(1181.25) * x * x * y * y * y * z * z * z) + e_1 * (std::sqrt(2657.8125) * x * x * x * x * y * z - std::sqrt(4725.0) * x * x * y * y * y * z + std::sqrt(10631.25) * x * x * y * z * z * z + std::sqrt(11.8125) * y * y * y * y * y * z - std::sqrt(1181.25) * y * y * y * z * z * z) + e_2 * (std::sqrt(42525.0) * x * x * y * z - std::sqrt(4725.0) * y * y * y * z);

        pc_42[k] = e_0 * (-std::sqrt(199.3359375) * x * x * x * x * x * y * z * z - std::sqrt(88.59375) * x * x * x * y * y * y * z * z + std::sqrt(1417.5) * x * x * x * y * z * z * z * z + std::sqrt(22.1484375) * x * y * y * y * y * y * z * z - std::sqrt(157.5) * x * y * y * y * z * z * z * z) + e_1 * (-std::sqrt(199.3359375) * x * x * x * x * x * y - std::sqrt(88.59375) * x * x * x * y * y * y + std::sqrt(3189.375) * x * x * x * y * z * z + std::sqrt(22.1484375) * x * y * y * y * y * y - std::sqrt(3189.375) * x * y * y * y * z * z + std::sqrt(5670.0) * x * y * z * z * z * z) + e_2 * (-std::sqrt(3189.375) * x * x * x * y - std::sqrt(354.375) * x * y * y * y + std::sqrt(51030.0) * x * y * z * z);

        pc_43[k] = e_0 * (std::sqrt(2.4609375) * x * x * x * x * x * x * y * z + std::sqrt(9.84375) * x * x * x * x * y * y * y * z - std::sqrt(630.0) * x * x * x * x * y * z * z * z + std::sqrt(2.4609375) * x * x * y * y * y * y * y * z - std::sqrt(630.0) * x * x * y * y * y * z * z * z + std::sqrt(630.0) * x * x * y * z * z * z * z * z) + e_1 * (-std::sqrt(1794.0234375) * x * x * x * x * y * z - std::sqrt(1663.59375) * x * x * y * y * y * z + std::sqrt(630.0) * x * x * y * z * z * z + std::sqrt(2.4609375) * y * y * y * y * y * z - std::sqrt(630.0) * y * y * y * z * z * z + std::sqrt(630.0) * y * z * z * z * z * z) + e_2 * (-std::sqrt(22680.0) * x * x * y * z - std::sqrt(2520.0) * y * y * y * z + std::sqrt(10080.0) * y * z * z * z);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];

        pc_44[k] = e_0 * (std::sqrt(24.609375) * x * x * x * x * x * y * z * z + std::sqrt(98.4375) * x * x * x * y * y * y * z * z - std::sqrt(393.75) * x * x * x * y * z * z * z * z + std::sqrt(24.609375) * x * y * y * y * y * y * z * z - std::sqrt(393.75) * x * y * y * y * z * z * z * z + std::sqrt(63.0) * x * y * z * z * z * z * z * z) + e_1 * (std::sqrt(24.609375) * x * x * x * x * x * y + std::sqrt(98.4375) * x * x * x * y * y * y - std::sqrt(1575.0) * x * x * x * y * z * z + std::sqrt(24.609375) * x * y * y * y * y * y - std::sqrt(1575.0) * x * y * y * y * z * z) + e_2 * (std::sqrt(393.75) * x * x * x * y + std::sqrt(393.75) * x * y * y * y - std::sqrt(14175.0) * x * y * z * z);

        pc_45[k] = e_0 * (-std::sqrt(0.29296875) * x * x * x * x * x * x * x * z - std::sqrt(2.63671875) * x * x * x * x * x * y * y * z + std::sqrt(94.921875) * x * x * x * x * x * z * z * z - std::sqrt(2.63671875) * x * x * x * y * y * y * y * z + std::sqrt(379.6875) * x * x * x * y * y * z * z * z - std::sqrt(168.75) * x * x * x * z * z * z * z * z - std::sqrt(0.29296875) * x * y * y * y * y * y * y * z + std::sqrt(94.921875) * x * y * y * y * y * z * z * z - std::sqrt(168.75) * x * y * y * z * z * z * z * z + std::sqrt(3.0) * x * z * z * z * z * z * z * z) + e_1 * (std::sqrt(263.671875) * x * x * x * x * x * z + std::sqrt(1054.6875) * x * x * x * y * y * z - std::sqrt(168.75) * x * x * x * z * z * z + std::sqrt(263.671875) * x * y * y * y * y * z - std::sqrt(168.75) * x * y * y * z * z * z - std::sqrt(243.0) * x * z * z * z * z * z) + e_2 * (std::sqrt(6075.0) * x * x * x * z + std::sqrt(6075.0) * x * y * y * z - std::sqrt(10800.0) * x * z * z * z);

        pc_46[k] = e_0 * (std::sqrt(24.609375) * x * x * x * x * x * x * z * z + std::sqrt(98.4375) * x * x * x * x * y * y * z * z - std::sqrt(393.75) * x * x * x * x * z * z * z * z + std::sqrt(24.609375) * x * x * y * y * y * y * z * z - std::sqrt(393.75) * x * x * y * y * z * z * z * z + std::sqrt(63.0) * x * x * z * z * z * z * z * z) + e_1 * (std::sqrt(24.609375) * x * x * x * x * x * x + std::sqrt(98.4375) * x * x * x * x * y * y - std::sqrt(1205.859375) * x * x * x * x * z * z + std::sqrt(24.609375) * x * x * y * y * y * y - std::sqrt(885.9375) * x * x * y * y * z * z - std::sqrt(393.75) * x * x * z * z * z * z + std::sqrt(24.609375) * y * y * y * y * z * z - std::sqrt(393.75) * y * y * z * z * z * z + std::sqrt(63.0) * z * z * z * z * z * z) + e_2 * (std::sqrt(615.234375) * x * x * x * x + std::sqrt(885.9375) * x * x * y * y - std::sqrt(31893.75) * x * x * z * z + std::sqrt(24.609375) * y * y * y * y - std::sqrt(3543.75) * y * y * z * z + std::sqrt(1575.0) * z * z * z * z);

        pc_47[k] = e_0 * (std::sqrt(0.615234375) * x * x * x * x * x * x * x * z + std::sqrt(0.615234375) * x * x * x * x * x * y * y * z - std::sqrt(157.5) * x * x * x * x * x * z * z * z - std::sqrt(0.615234375) * x * x * x * y * y * y * y * z + std::sqrt(157.5) * x * x * x * z * z * z * z * z - std::sqrt(0.615234375) * x * y * y * y * y * y * y * z + std::sqrt(157.5) * x * y * y * y * y * z * z * z - std::sqrt(157.5) * x * y * y * z * z * z * z * z) + e_1 * (-std::sqrt(415.8984375) * x * x * x * x * x * z + std::sqrt(9.84375) * x * x * x * y * y * z + std::sqrt(553.7109375) * x * y * y * y * y * z - std::sqrt(2520.0) * x * y * y * z * z * z + std::sqrt(630.0) * x * z * z * z * z * z) + e_2 * (-std::sqrt(10080.0) * x * x * x * z + std::sqrt(10080.0) * x * z * z * z);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];

        pc_48[k] = e_0 * (-std::sqrt(22.1484375) * x * x * x * x * x * x * z * z + std::sqrt(88.59375) * x * x * x * x * y * y * z * z + std::sqrt(157.5) * x * x * x * x * z * z * z * z + std::sqrt(199.3359375) * x * x * y * y * y * y * z * z - std::sqrt(1417.5) * x * x * y * y * z * z * z * z) + e_1 * (-std::sqrt(22.1484375) * x * x * x * x * x * x + std::sqrt(88.59375) * x * x * x * x * y * y + std::sqrt(199.3359375) * x * x * x * x * z * z + std::sqrt(199.3359375) * x * x * y * y * y * y - std::sqrt(7176.09375) * x * x * y * y * z * z + std::sqrt(1417.5) * x * x * z * z * z * z + std::sqrt(199.3359375) * y * y * y * y * z * z - std::sqrt(1417.5) * y * y * z * z * z * z) + e_2 * (-std::sqrt(553.7109375) * x * x * x * x + std::sqrt(797.34375) * x * x * y * y + std::sqrt(12757.5) * x * x * z * z + std::sqrt(199.3359375) * y * y * y * y - std::sqrt(12757.5) * y * y * z * z);

        pc_49[k] = e_0 * (-std::sqrt(0.73828125) * x * x * x * x * x * x * x * z + std::sqrt(18.45703125) * x * x * x * x * x * y * y * z + std::sqrt(73.828125) * x * x * x * x * x * z * z * z + std::sqrt(18.45703125) * x * x * x * y * y * y * y * z - std::sqrt(2657.8125) * x * x * x * y * y * z * z * z - std::sqrt(0.73828125) * x * y * y * y * y * y * y * z + std::sqrt(73.828125) * x * y * y * y * y * z * z * z) + e_1 * (std::sqrt(144.703125) * x * x * x * x * x * z - std::sqrt(7382.8125) * x * x * x * y * y * z + std::sqrt(1181.25) * x * x * x * z * z * z + std::sqrt(664.453125) * x * y * y * y * y * z - std::sqrt(10631.25) * x * y * y * z * z * z) + e_2 * (std::sqrt(4725.0) * x * x * x * z - std::sqrt(42525.0) * x * y * y * z);

        pc_50[k] = e_0 * (std::sqrt(16.2421875) * x * x * x * x * x * x * z * z - std::sqrt(1624.21875) * x * x * x * x * y * y * z * z + std::sqrt(406.0546875) * x * x * y * y * y * y * z * z) + e_1 * (std::sqrt(16.2421875) * x * x * x * x * x * x - std::sqrt(1624.21875) * x * x * x * x * y * y + std::sqrt(406.0546875) * x * x * x * x * z * z + std::sqrt(406.0546875) * x * x * y * y * y * y - std::sqrt(14617.96875) * x * x * y * y * z * z + std::sqrt(406.0546875) * y * y * y * y * z * z) + e_2 * (std::sqrt(406.0546875) * x * x * x * x - std::sqrt(14617.96875) * x * x * y * y + std::sqrt(406.0546875) * y * y * y * y);

        pc_51[k] = e_0 * (std::sqrt(1.353515625) * x * x * x * x * x * x * x * z - std::sqrt(304.541015625) * x * x * x * x * x * y * y * z + std::sqrt(304.541015625) * x * x * x * y * y * y * y * z - std::sqrt(1.353515625) * x * y * y * y * y * y * y * z) + e_1 * (std::sqrt(48.7265625) * x * x * x * x * x * z - std::sqrt(4872.65625) * x * x * x * y * y * z + std::sqrt(1218.1640625) * x * y * y * y * y * z);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];

        pc_52[k] = e_0 * (std::sqrt(12.181640625) * x * x * x * x * x * x * x * y - std::sqrt(228.744140625) * x * x * x * x * x * y * y * y + std::sqrt(228.744140625) * x * x * x * y * y * y * y * y - std::sqrt(12.181640625) * x * y * y * y * y * y * y * y) + e_1 * (std::sqrt(779.625) * x * x * x * x * x * y - std::sqrt(779.625) * x * y * y * y * y * y) + e_2 * (std::sqrt(19490.625) * x * x * x * y - std::sqrt(19490.625) * x * y * y * y);

        pc_53[k] = e_0 * (std::sqrt(101.513671875) * x * x * x * x * x * x * y * z - std::sqrt(913.623046875) * x * x * x * x * y * y * y * z + std::sqrt(491.326171875) * x * x * y * y * y * y * y * z - std::sqrt(4.060546875) * y * y * y * y * y * y * y * z) + e_1 * (std::sqrt(3654.4921875) * x * x * x * x * y * z + std::sqrt(1624.21875) * x * x * y * y * y * z - std::sqrt(406.0546875) * y * y * y * y * y * z) + e_2 * (std::sqrt(58471.875) * x * x * y * z - std::sqrt(6496.875) * y * y * y * z);

        pc_54[k] = e_0 * (-std::sqrt(2.953125) * x * x * x * x * x * x * x * y + std::sqrt(2.953125) * x * x * x * x * x * y * y * y + std::sqrt(295.3125) * x * x * x * x * x * y * z * z + std::sqrt(2.953125) * x * x * x * y * y * y * y * y - std::sqrt(1181.25) * x * x * x * y * y * y * z * z - std::sqrt(2.953125) * x * y * y * y * y * y * y * y + std::sqrt(295.3125) * x * y * y * y * y * y * z * z) + e_1 * (-std::sqrt(189.0) * x * x * x * x * x * y + std::sqrt(4725.0) * x * x * x * y * z * z - std::sqrt(189.0) * x * y * y * y * y * y + std::sqrt(4725.0) * x * y * y * y * z * z) + e_2 * (-std::sqrt(1181.25) * x * x * x * y - std::sqrt(1181.25) * x * y * y * y + std::sqrt(42525.0) * x * y * z * z);

        pc_55[k] = e_0 * (-std::sqrt(49.833984375) * x * x * x * x * x * x * y * z + std::sqrt(5.537109375) * x * x * x * x * y * y * y * z + std::sqrt(354.375) * x * x * x * x * y * z * z * z + std::sqrt(49.833984375) * x * x * y * y * y * y * y * z - std::sqrt(630.0) * x * x * y * y * y * z * z * z - std::sqrt(5.537109375) * y * y * y * y * y * y * y * z + std::sqrt(39.375) * y * y * y * y * y * z * z * z) + e_1 * (-std::sqrt(1794.0234375) * x * x * x * x * y * z + std::sqrt(88.59375) * x * x * y * y * y * z + std::sqrt(1417.5) * x * x * y * z * z * z - std::sqrt(553.7109375) * y * y * y * y * y * z + std::sqrt(1417.5) * y * y * y * z * z * z) + e_2 * (-std::sqrt(3189.375) * x * x * y * z - std::sqrt(3189.375) * y * y * y * z + std::sqrt(5670.0) * y * z * z * z);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];

        pc_56[k] = e_0 * (std::sqrt(0.615234375) * x * x * x * x * x * x * x * y + std::sqrt(0.615234375) * x * x * x * x * x * y * y * y - std::sqrt(157.5) * x * x * x * x * x * y * z * z - std::sqrt(0.615234375) * x * x * x * y * y * y * y * y + std::sqrt(157.5) * x * x * x * y * z * z * z * z - std::sqrt(0.615234375) * x * y * y * y * y * y * y * y + std::sqrt(157.5) * x * y * y * y * y * y * z * z - std::sqrt(157.5) * x * y * y * y * z * z * z * z) + e_1 * (std::sqrt(39.375) * x * x * x * x * x * y - std::sqrt(2520.0) * x * x * x * y * z * z - std::sqrt(39.375) * x * y * y * y * y * y + std::sqrt(2520.0) * x * y * y * y * z * z) + e_2 * (std::sqrt(39.375) * x * x * x * y - std::sqrt(39.375) * x * y * y * y);

        pc_57[k] = e_0 * (std::sqrt(6.15234375) * x * x * x * x * x * x * y * z + std::sqrt(6.15234375) * x * x * x * x * y * y * y * z - std::sqrt(98.4375) * x * x * x * x * y * z * z * z - std::sqrt(6.15234375) * x * x * y * y * y * y * y * z + std::sqrt(15.75) * x * x * y * z * z * z * z * z - std::sqrt(6.15234375) * y * y * y * y * y * y * y * z + std::sqrt(98.4375) * y * y * y * y * y * z * z * z - std::sqrt(15.75) * y * y * y * z * z * z * z * z) + e_1 * (std::sqrt(221.484375) * x * x * x * x * y * z - std::sqrt(98.4375) * x * x * y * y * y * z - std::sqrt(393.75) * x * x * y * z * z * z - std::sqrt(615.234375) * y * y * y * y * y * z + std::sqrt(3543.75) * y * y * y * z * z * z - std::sqrt(63.0) * y * z * z * z * z * z) + e_2 * (-std::sqrt(1575.0) * y * y * y * z + std::sqrt(1575.0) * y * z * z * z);

        pc_58[k] = e_0 * (-std::sqrt(0.0732421875) * x * x * x * x * x * x * x * x - std::sqrt(0.29296875) * x * x * x * x * x * x * y * y + std::sqrt(23.73046875) * x * x * x * x * x * x * z * z + std::sqrt(23.73046875) * x * x * x * x * y * y * z * z - std::sqrt(42.1875) * x * x * x * x * z * z * z * z + std::sqrt(0.29296875) * x * x * y * y * y * y * y * y - std::sqrt(23.73046875) * x * x * y * y * y * y * z * z + std::sqrt(0.75) * x * x * z * z * z * z * z * z + std::sqrt(0.0732421875) * y * y * y * y * y * y * y * y - std::sqrt(23.73046875) * y * y * y * y * y * y * z * z + std::sqrt(42.1875) * y * y * y * y * z * z * z * z - std::sqrt(0.75) * y * y * z * z * z * z * z * z) + e_1 * (-std::sqrt(10.546875) * x * x * x * x * x * x - std::sqrt(10.546875) * x * x * x * x * y * y + std::sqrt(1518.75) * x * x * x * x * z * z + std::sqrt(10.546875) * x * x * y * y * y * y - std::sqrt(675.0) * x * x * z * z * z * z + std::sqrt(10.546875) * y * y * y * y * y * y - std::sqrt(1518.75) * y * y * y * y * z * z + std::sqrt(675.0) * y * y * z * z * z * z) + e_2 * (-std::sqrt(42.1875) * x * x * x * x + std::sqrt(1518.75) * x * x * z * z + std::sqrt(42.1875) * y * y * y * y - std::sqrt(1518.75) * y * y * z * z);

        pc_59[k] = e_0 * (std::sqrt(6.15234375) * x * x * x * x * x * x * x * z + std::sqrt(6.15234375) * x * x * x * x * x * y * y * z - std::sqrt(98.4375) * x * x * x * x * x * z * z * z - std::sqrt(6.15234375) * x * x * x * y * y * y * y * z + std::sqrt(15.75) * x * x * x * z * z * z * z * z - std::sqrt(6.15234375) * x * y * y * y * y * y * y * z + std::sqrt(98.4375) * x * y * y * y * y * z * z * z - std::sqrt(15.75) * x * y * y * z * z * z * z * z) + e_1 * (std::sqrt(615.234375) * x * x * x * x * x * z + std::sqrt(98.4375) * x * x * x * y * y * z - std::sqrt(3543.75) * x * x * x * z * z * z - std::sqrt(221.484375) * x * y * y * y * y * z + std::sqrt(393.75) * x * y * y * z * z * z + std::sqrt(63.0) * x * z * z * z * z * z) + e_2 * (std::sqrt(1575.0) * x * x * x * z - std::sqrt(1575.0) * x * z * z * z);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];

        pc_60[k] = e_0 * (std::sqrt(0.15380859375) * x * x * x * x * x * x * x * x - std::sqrt(39.375) * x * x * x * x * x * x * z * z - std::sqrt(0.615234375) * x * x * x * x * y * y * y * y + std::sqrt(39.375) * x * x * x * x * y * y * z * z + std::sqrt(39.375) * x * x * x * x * z * z * z * z + std::sqrt(39.375) * x * x * y * y * y * y * z * z - std::sqrt(157.5) * x * x * y * y * z * z * z * z + std::sqrt(0.15380859375) * y * y * y * y * y * y * y * y - std::sqrt(39.375) * y * y * y * y * y * y * z * z + std::sqrt(39.375) * y * y * y * y * z * z * z * z) + e_1 * (std::sqrt(22.1484375) * x * x * x * x * x * x + std::sqrt(2.4609375) * x * x * x * x * y * y - std::sqrt(2520.0) * x * x * x * x * z * z + std::sqrt(2.4609375) * x * x * y * y * y * y + std::sqrt(630.0) * x * x * z * z * z * z + std::sqrt(22.1484375) * y * y * y * y * y * y - std::sqrt(2520.0) * y * y * y * y * z * z + std::sqrt(630.0) * y * y * z * z * z * z) + e_2 * (std::sqrt(120.5859375) * x * x * x * x + std::sqrt(88.59375) * x * x * y * y - std::sqrt(5670.0) * x * x * z * z + std::sqrt(120.5859375) * y * y * y * y - std::sqrt(5670.0) * y * y * z * z + std::sqrt(630.0) * z * z * z * z);

        pc_61[k] = e_0 * (-std::sqrt(5.537109375) * x * x * x * x * x * x * x * z + std::sqrt(49.833984375) * x * x * x * x * x * y * y * z + std::sqrt(39.375) * x * x * x * x * x * z * z * z + std::sqrt(5.537109375) * x * x * x * y * y * y * y * z - std::sqrt(630.0) * x * x * x * y * y * z * z * z - std::sqrt(49.833984375) * x * y * y * y * y * y * y * z + std::sqrt(354.375) * x * y * y * y * y * z * z * z) + e_1 * (-std::sqrt(553.7109375) * x * x * x * x * x * z + std::sqrt(88.59375) * x * x * x * y * y * z + std::sqrt(1417.5) * x * x * x * z * z * z - std::sqrt(1794.0234375) * x * y * y * y * y * z + std::sqrt(1417.5) * x * y * y * z * z * z) + e_2 * (-std::sqrt(3189.375) * x * x * x * z - std::sqrt(3189.375) * x * y * y * z + std::sqrt(5670.0) * x * z * z * z);

        pc_62[k] = e_0 * (-std::sqrt(0.1845703125) * x * x * x * x * x * x * x * x + std::sqrt(6.64453125) * x * x * x * x * x * x * y * y + std::sqrt(18.45703125) * x * x * x * x * x * x * z * z - std::sqrt(904.39453125) * x * x * x * x * y * y * z * z - std::sqrt(6.64453125) * x * x * y * y * y * y * y * y + std::sqrt(904.39453125) * x * x * y * y * y * y * z * z + std::sqrt(0.1845703125) * y * y * y * y * y * y * y * y - std::sqrt(18.45703125) * y * y * y * y * y * y * z * z) + e_1 * (-std::sqrt(26.578125) * x * x * x * x * x * x + std::sqrt(73.828125) * x * x * x * x * y * y + std::sqrt(1181.25) * x * x * x * x * z * z - std::sqrt(73.828125) * x * x * y * y * y * y + std::sqrt(26.578125) * y * y * y * y * y * y - std::sqrt(1181.25) * y * y * y * y * z * z) + e_2 * (-std::sqrt(295.3125) * x * x * x * x + std::sqrt(10631.25) * x * x * z * z + std::sqrt(295.3125) * y * y * y * y - std::sqrt(10631.25) * y * y * z * z);

        pc_63[k] = e_0 * (std::sqrt(4.060546875) * x * x * x * x * x * x * x * z - std::sqrt(491.326171875) * x * x * x * x * x * y * y * z + std::sqrt(913.623046875) * x * x * x * y * y * y * y * z - std::sqrt(101.513671875) * x * y * y * y * y * y * y * z) + e_1 * (std::sqrt(406.0546875) * x * x * x * x * x * z - std::sqrt(1624.21875) * x * x * x * y * y * z - std::sqrt(3654.4921875) * x * y * y * y * y * z) + e_2 * (std::sqrt(6496.875) * x * x * x * z - std::sqrt(58471.875) * x * y * y * z);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, ab_x, ab_y : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];

        pc_64[k] = e_0 * (std::sqrt(0.33837890625) * x * x * x * x * x * x * x * x - std::sqrt(86.625) * x * x * x * x * x * x * y * y + std::sqrt(304.541015625) * x * x * x * x * y * y * y * y - std::sqrt(86.625) * x * x * y * y * y * y * y * y + std::sqrt(0.33837890625) * y * y * y * y * y * y * y * y) + e_1 * (std::sqrt(48.7265625) * x * x * x * x * x * x - std::sqrt(1218.1640625) * x * x * x * x * y * y - std::sqrt(1218.1640625) * x * x * y * y * y * y + std::sqrt(48.7265625) * y * y * y * y * y * y) + e_2 * (std::sqrt(1218.1640625) * x * x * x * x - std::sqrt(43853.90625) * x * x * y * y + std::sqrt(1218.1640625) * y * y * y * y);
    }

    // NOTE: the atom pairs beyond the reach of every pair of primitives have no
    // contribution and are set to zero.

    for (size_t m = 0; m < 65; m++)
    {
        auto *pv = values + m * nvalues;

        std::fill(pv + nmax, pv + nvalues, 0.0);
    }
}

}  // namespace simdovl
