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



#include "SimdOverlapRecIF.hpp"

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
compute_if_overlap(double               *values,
                   const size_t          nvalues,
                   const CBasisFunction &bra,
                   const CBasisFunction &ket,
                   const CSimdMatrix    &coordinates,
                   const double          threshold) -> void
{
    if ((bra.get_angular_momentum() != 6) || (ket.get_angular_momentum() != 3))
    {
        errors::assertMsgCritical(
            false, std::string("SimdOverlapRecIF.compute_if_overlap: Basis functions must be of angular momenta six and three"));
    }

    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("SimdOverlapRecIF.compute_if_overlap: Number of values exceeds number of atom pairs"));
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

    auto buffer = simdfunc::make_primitive_buffer(dimensions, 4);

    if (buffer.number_of_columns() == 0)
    {
        std::fill(values, values + 91 * nvalues, 0.0);

        return;
    }

    const auto nmax = buffer.number_of_columns();

    auto *pe_0 = buffer.data(0);
    auto *pe_1 = buffer.data(1);
    auto *pe_2 = buffer.data(2);
    auto *pe_3 = buffer.data(3);

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

        const auto f_0 = fbase * fal * fal * fal * fal * fal * fal * fbe * fbe * fbe;

        const auto f_1 = fbase * fal * fal * fal * fal * fal * fbe * fbe * fh;

        const auto f_2 = fbase * fal * fal * fal * fal * fbe * fh * fh;

        const auto f_3 = fbase * fal * fal * fal * fh * fh * fh;

        // NOTE: the exponential depends on the pair of primitives alone, so it is
        // evaluated once and shared by the prefactors of all terms.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ab_2 : simd::cache_line_size())
        for (size_t k = 0; k < ncols; k++)
        {
            const auto fss = std::exp(-fmu * ab_2[k]);

            pe_0[k] += f_0 * fss;
            pe_1[k] += f_1 * fss;
            pe_2[k] += f_2 * fss;
            pe_3[k] += f_3 * fss;
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
    auto *pc_65 = values + 65 * nvalues;
    auto *pc_66 = values + 66 * nvalues;
    auto *pc_67 = values + 67 * nvalues;
    auto *pc_68 = values + 68 * nvalues;
    auto *pc_69 = values + 69 * nvalues;
    auto *pc_70 = values + 70 * nvalues;
    auto *pc_71 = values + 71 * nvalues;
    auto *pc_72 = values + 72 * nvalues;
    auto *pc_73 = values + 73 * nvalues;
    auto *pc_74 = values + 74 * nvalues;
    auto *pc_75 = values + 75 * nvalues;
    auto *pc_76 = values + 76 * nvalues;
    auto *pc_77 = values + 77 * nvalues;
    auto *pc_78 = values + 78 * nvalues;
    auto *pc_79 = values + 79 * nvalues;
    auto *pc_80 = values + 80 * nvalues;
    auto *pc_81 = values + 81 * nvalues;
    auto *pc_82 = values + 82 * nvalues;
    auto *pc_83 = values + 83 * nvalues;
    auto *pc_84 = values + 84 * nvalues;
    auto *pc_85 = values + 85 * nvalues;
    auto *pc_86 = values + 86 * nvalues;
    auto *pc_87 = values + 87 * nvalues;
    auto *pc_88 = values + 88 * nvalues;
    auto *pc_89 = values + 89 * nvalues;
    auto *pc_90 = values + 90 * nvalues;

    // NOTE: the components are formed in 51 loops, as the vectorizer runs out
    // of registers with all of them in one. Only the prefactors and the vector
    // between the atoms are loaded by more than one loop.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        pc_0[k] = e_0 * (std::sqrt(91.3623046875) * x * x * x * x * x * x * x * y * y - std::sqrt(1228.3154296875) * x * x * x * x * x * y * y * y * y + std::sqrt(407.1826171875) * x * x * x * y * y * y * y * y * y - std::sqrt(10.1513671875) * x * y * y * y * y * y * y * y * y) + e_1 * (std::sqrt(91.3623046875) * x * x * x * x * x * x * x - std::sqrt(91.3623046875) * x * x * x * x * x * y * y - std::sqrt(2284.0576171875) * x * x * x * y * y * y * y - std::sqrt(822.2607421875) * x * y * y * y * y * y * y) + e_2 * (std::sqrt(9136.23046875) * x * x * x * x * x - std::sqrt(36544.921875) * x * x * x * y * y - std::sqrt(82226.07421875) * x * y * y * y * y) + e_3 * (std::sqrt(64968.75) * x * x * x - std::sqrt(584718.75) * x * y * y);

        pc_1[k] = e_0 * (std::sqrt(243.6328125) * x * x * x * x * x * x * y * y * z - std::sqrt(2707.03125) * x * x * x * x * y * y * y * y * z + std::sqrt(243.6328125) * x * x * y * y * y * y * y * y * z) + e_1 * (std::sqrt(243.6328125) * x * x * x * x * x * x * z - std::sqrt(6090.8203125) * x * x * x * x * y * y * z - std::sqrt(6090.8203125) * x * x * y * y * y * y * z + std::sqrt(243.6328125) * y * y * y * y * y * y * z) + e_2 * (std::sqrt(6090.8203125) * x * x * x * x * z - std::sqrt(219269.53125) * x * x * y * y * z + std::sqrt(6090.8203125) * y * y * y * y * z);

        pc_2[k] = e_0 * (-std::sqrt(6.0908203125) * x * x * x * x * x * x * x * y * y + std::sqrt(33.1611328125) * x * x * x * x * x * y * y * y * y + std::sqrt(97.453125) * x * x * x * x * x * y * y * z * z + std::sqrt(33.1611328125) * x * x * x * y * y * y * y * y * y - std::sqrt(1082.8125) * x * x * x * y * y * y * y * z * z - std::sqrt(6.0908203125) * x * y * y * y * y * y * y * y * y + std::sqrt(97.453125) * x * y * y * y * y * y * y * z * z) + e_1 * (-std::sqrt(6.0908203125) * x * x * x * x * x * x * x - std::sqrt(54.8173828125) * x * x * x * x * x * y * y + std::sqrt(97.453125) * x * x * x * x * x * z * z + std::sqrt(12333.9111328125) * x * x * x * y * y * y * y - std::sqrt(9745.3125) * x * x * x * y * y * z * z - std::sqrt(1760.2470703125) * x * y * y * y * y * y * y + std::sqrt(2436.328125) * x * y * y * y * y * z * z) + e_2 * (-std::sqrt(609.08203125) * x * x * x * x * x + std::sqrt(60908.203125) * x * x * x * y * y - std::sqrt(15227.05078125) * x * y * y * y * y);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        pc_3[k] = e_0 * (-std::sqrt(36.544921875) * x * x * x * x * x * x * x * y * z + std::sqrt(198.966796875) * x * x * x * x * x * y * y * y * z + std::sqrt(16.2421875) * x * x * x * x * x * y * z * z * z + std::sqrt(198.966796875) * x * x * x * y * y * y * y * y * z - std::sqrt(180.46875) * x * x * x * y * y * y * z * z * z - std::sqrt(36.544921875) * x * y * y * y * y * y * y * y * z + std::sqrt(16.2421875) * x * y * y * y * y * y * z * z * z) + e_1 * (-std::sqrt(5262.46875) * x * x * x * x * x * y * z + std::sqrt(58471.875) * x * x * x * y * y * y * z - std::sqrt(5262.46875) * x * y * y * y * y * y * z);

        pc_4[k] = e_0 * (-std::sqrt(6.0908203125) * x * x * x * x * x * x * x * x * y + std::sqrt(33.1611328125) * x * x * x * x * x * x * y * y * y + std::sqrt(97.453125) * x * x * x * x * x * x * y * z * z + std::sqrt(33.1611328125) * x * x * x * x * y * y * y * y * y - std::sqrt(1082.8125) * x * x * x * x * y * y * y * z * z - std::sqrt(6.0908203125) * x * x * y * y * y * y * y * y * y + std::sqrt(97.453125) * x * x * y * y * y * y * y * z * z) + e_1 * (-std::sqrt(1760.2470703125) * x * x * x * x * x * x * y + std::sqrt(12333.9111328125) * x * x * x * x * y * y * y + std::sqrt(2436.328125) * x * x * x * x * y * z * z - std::sqrt(54.8173828125) * x * x * y * y * y * y * y - std::sqrt(9745.3125) * x * x * y * y * y * z * z - std::sqrt(6.0908203125) * y * y * y * y * y * y * y + std::sqrt(97.453125) * y * y * y * y * y * z * z) + e_2 * (-std::sqrt(15227.05078125) * x * x * x * x * y + std::sqrt(60908.203125) * x * x * y * y * y - std::sqrt(609.08203125) * y * y * y * y * y);

        pc_5[k] = e_0 * (std::sqrt(60.908203125) * x * x * x * x * x * x * x * y * z - std::sqrt(1143.720703125) * x * x * x * x * x * y * y * y * z + std::sqrt(1143.720703125) * x * x * x * y * y * y * y * y * z - std::sqrt(60.908203125) * x * y * y * y * y * y * y * y * z) + e_1 * (std::sqrt(3898.125) * x * x * x * x * x * y * z - std::sqrt(3898.125) * x * y * y * y * y * y * z) + e_2 * (std::sqrt(97453.125) * x * x * x * y * z - std::sqrt(97453.125) * x * y * y * y * z);

        pc_6[k] = e_0 * (std::sqrt(10.1513671875) * x * x * x * x * x * x * x * x * y - std::sqrt(407.1826171875) * x * x * x * x * x * x * y * y * y + std::sqrt(1228.3154296875) * x * x * x * x * y * y * y * y * y - std::sqrt(91.3623046875) * x * x * y * y * y * y * y * y * y) + e_1 * (std::sqrt(822.2607421875) * x * x * x * x * x * x * y + std::sqrt(2284.0576171875) * x * x * x * x * y * y * y + std::sqrt(91.3623046875) * x * x * y * y * y * y * y - std::sqrt(91.3623046875) * y * y * y * y * y * y * y) + e_2 * (std::sqrt(82226.07421875) * x * x * x * x * y + std::sqrt(36544.921875) * x * x * y * y * y - std::sqrt(9136.23046875) * y * y * y * y * y) + e_3 * (std::sqrt(584718.75) * x * x * y - std::sqrt(64968.75) * y * y * y);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        pc_7[k] = e_0 * (std::sqrt(761.3525390625) * x * x * x * x * x * x * y * y * z - std::sqrt(4145.1416015625) * x * x * x * x * y * y * y * y * z + std::sqrt(571.8603515625) * x * x * y * y * y * y * y * y * z - std::sqrt(3.3837890625) * y * y * y * y * y * y * y * y * z) + e_1 * (std::sqrt(761.3525390625) * x * x * x * x * x * x * z + std::sqrt(761.3525390625) * x * x * x * x * y * y * z - std::sqrt(761.3525390625) * x * x * y * y * y * y * z - std::sqrt(761.3525390625) * y * y * y * y * y * y * z) + e_2 * (std::sqrt(48726.5625) * x * x * x * x * z - std::sqrt(48726.5625) * y * y * y * y * z) + e_3 * (std::sqrt(194906.25) * x * x * z - std::sqrt(194906.25) * y * y * z);

        pc_8[k] = e_0 * (std::sqrt(2030.2734375) * x * x * x * x * x * y * y * z * z - std::sqrt(8121.09375) * x * x * x * y * y * y * y * z * z + std::sqrt(81.2109375) * x * y * y * y * y * y * y * z * z) + e_1 * (std::sqrt(2030.2734375) * x * x * x * x * x * y * y + std::sqrt(2030.2734375) * x * x * x * x * x * z * z - std::sqrt(8121.09375) * x * x * x * y * y * y * y - std::sqrt(8121.09375) * x * x * x * y * y * z * z + std::sqrt(81.2109375) * x * y * y * y * y * y * y - std::sqrt(18272.4609375) * x * y * y * y * y * z * z) + e_2 * (std::sqrt(2030.2734375) * x * x * x * x * x - std::sqrt(8121.09375) * x * x * x * y * y + std::sqrt(32484.375) * x * x * x * z * z - std::sqrt(18272.4609375) * x * y * y * y * y - std::sqrt(292359.375) * x * y * y * z * z) + e_3 * (std::sqrt(32484.375) * x * x * x - std::sqrt(292359.375) * x * y * y);

        pc_9[k] = e_0 * (-std::sqrt(50.7568359375) * x * x * x * x * x * x * y * y * z + std::sqrt(50.7568359375) * x * x * x * x * y * y * y * y * z + std::sqrt(812.109375) * x * x * x * x * y * y * z * z * z + std::sqrt(164.4521484375) * x * x * y * y * y * y * y * y * z - std::sqrt(3248.4375) * x * x * y * y * y * y * z * z * z - std::sqrt(2.0302734375) * y * y * y * y * y * y * y * y * z + std::sqrt(32.484375) * y * y * y * y * y * y * z * z * z) + e_1 * (-std::sqrt(50.7568359375) * x * x * x * x * x * x * z + std::sqrt(456.8115234375) * x * x * x * x * y * y * z + std::sqrt(812.109375) * x * x * x * x * z * z * z + std::sqrt(4111.3037109375) * x * x * y * y * y * y * z - std::sqrt(29235.9375) * x * x * y * y * z * z * z - std::sqrt(99.4833984375) * y * y * y * y * y * y * z + std::sqrt(812.109375) * y * y * y * y * z * z * z);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        pc_10[k] = e_0 * (-std::sqrt(304.541015625) * x * x * x * x * x * x * y * z * z + std::sqrt(304.541015625) * x * x * x * x * y * y * y * z * z + std::sqrt(135.3515625) * x * x * x * x * y * z * z * z * z + std::sqrt(986.712890625) * x * x * y * y * y * y * y * z * z - std::sqrt(541.40625) * x * x * y * y * y * z * z * z * z - std::sqrt(12.181640625) * y * y * y * y * y * y * y * z * z + std::sqrt(5.4140625) * y * y * y * y * y * z * z * z * z) + e_1 * (-std::sqrt(304.541015625) * x * x * x * x * x * x * y + std::sqrt(304.541015625) * x * x * x * x * y * y * y - std::sqrt(19490.625) * x * x * x * x * y * z * z + std::sqrt(986.712890625) * x * x * y * y * y * y * y + std::sqrt(77962.5) * x * x * y * y * y * z * z - std::sqrt(12.181640625) * y * y * y * y * y * y * y - std::sqrt(779.625) * y * y * y * y * y * z * z) + e_2 * (-std::sqrt(30454.1015625) * x * x * x * x * y + std::sqrt(121816.40625) * x * x * y * y * y - std::sqrt(1218.1640625) * y * y * y * y * y);

        pc_11[k] = e_0 * (-std::sqrt(50.7568359375) * x * x * x * x * x * x * x * y * z + std::sqrt(50.7568359375) * x * x * x * x * x * y * y * y * z + std::sqrt(812.109375) * x * x * x * x * x * y * z * z * z + std::sqrt(164.4521484375) * x * x * x * y * y * y * y * y * z - std::sqrt(3248.4375) * x * x * x * y * y * y * z * z * z - std::sqrt(2.0302734375) * x * y * y * y * y * y * y * y * z + std::sqrt(32.484375) * x * y * y * y * y * y * z * z * z) + e_1 * (-std::sqrt(1827.24609375) * x * x * x * x * x * y * z + std::sqrt(812.109375) * x * x * x * y * y * y * z + std::sqrt(12993.75) * x * x * x * y * z * z * z + std::sqrt(657.80859375) * x * y * y * y * y * y * z - std::sqrt(12993.75) * x * y * y * y * z * z * z);

        pc_12[k] = e_0 * (std::sqrt(507.568359375) * x * x * x * x * x * x * y * z * z - std::sqrt(4568.115234375) * x * x * x * x * y * y * y * z * z + std::sqrt(2456.630859375) * x * x * y * y * y * y * y * z * z - std::sqrt(20.302734375) * y * y * y * y * y * y * y * z * z) + e_1 * (std::sqrt(507.568359375) * x * x * x * x * x * x * y - std::sqrt(4568.115234375) * x * x * x * x * y * y * y + std::sqrt(18272.4609375) * x * x * x * x * y * z * z + std::sqrt(2456.630859375) * x * x * y * y * y * y * y + std::sqrt(8121.09375) * x * x * y * y * y * z * z - std::sqrt(20.302734375) * y * y * y * y * y * y * y - std::sqrt(2030.2734375) * y * y * y * y * y * z * z) + e_2 * (std::sqrt(18272.4609375) * x * x * x * x * y + std::sqrt(8121.09375) * x * x * y * y * y + std::sqrt(292359.375) * x * x * y * z * z - std::sqrt(2030.2734375) * y * y * y * y * y - std::sqrt(32484.375) * y * y * y * z * z) + e_3 * (std::sqrt(292359.375) * x * x * y - std::sqrt(32484.375) * y * y * y);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        pc_13[k] = e_0 * (std::sqrt(84.5947265625) * x * x * x * x * x * x * x * y * z - std::sqrt(2114.8681640625) * x * x * x * x * x * y * y * y * z + std::sqrt(3251.8212890625) * x * x * x * y * y * y * y * y * z - std::sqrt(30.4541015625) * x * y * y * y * y * y * y * y * z) + e_1 * (std::sqrt(3045.41015625) * x * x * x * x * x * y * z + std::sqrt(12181.640625) * x * x * x * y * y * y * z + std::sqrt(3045.41015625) * x * y * y * y * y * y * z) + e_2 * (std::sqrt(194906.25) * x * x * x * y * z + std::sqrt(194906.25) * x * y * y * y * z) + e_3 * (std::sqrt(779625.0) * x * y * z);

        pc_14[k] = e_0 * (-std::sqrt(22.1484375) * x * x * x * x * x * x * x * y * y + std::sqrt(2.4609375) * x * x * x * x * x * y * y * y * y + std::sqrt(2214.84375) * x * x * x * x * x * y * y * z * z + std::sqrt(22.1484375) * x * x * x * y * y * y * y * y * y - std::sqrt(3937.5) * x * x * x * y * y * y * y * z * z - std::sqrt(2.4609375) * x * y * y * y * y * y * y * y * y + std::sqrt(246.09375) * x * y * y * y * y * y * y * z * z) + e_1 * (-std::sqrt(22.1484375) * x * x * x * x * x * x * x - std::sqrt(1794.0234375) * x * x * x * x * x * y * y + std::sqrt(2214.84375) * x * x * x * x * x * z * z + std::sqrt(553.7109375) * x * x * x * y * y * y * y + std::sqrt(8859.375) * x * x * x * y * y * z * z - std::sqrt(199.3359375) * x * y * y * y * y * y * y + std::sqrt(2214.84375) * x * y * y * y * y * z * z) + e_2 * (-std::sqrt(2214.84375) * x * x * x * x * x - std::sqrt(8859.375) * x * x * x * y * y + std::sqrt(79734.375) * x * x * x * z * z - std::sqrt(2214.84375) * x * y * y * y * y + std::sqrt(79734.375) * x * y * y * z * z) + e_3 * (-std::sqrt(8859.375) * x * x * x - std::sqrt(8859.375) * x * y * y + std::sqrt(141750.0) * x * z * z);

        pc_15[k] = e_0 * (-std::sqrt(59.0625) * x * x * x * x * x * x * y * y * z + std::sqrt(5906.25) * x * x * x * x * y * y * z * z * z + std::sqrt(59.0625) * x * x * y * y * y * y * y * y * z - std::sqrt(5906.25) * x * x * y * y * y * y * z * z * z) + e_1 * (-std::sqrt(59.0625) * x * x * x * x * x * x * z + std::sqrt(13289.0625) * x * x * x * x * y * y * z + std::sqrt(5906.25) * x * x * x * x * z * z * z - std::sqrt(13289.0625) * x * x * y * y * y * y * z + std::sqrt(59.0625) * y * y * y * y * y * y * z - std::sqrt(5906.25) * y * y * y * y * z * z * z) + e_2 * (std::sqrt(13289.0625) * x * x * x * x * z + std::sqrt(53156.25) * x * x * z * z * z - std::sqrt(13289.0625) * y * y * y * y * z - std::sqrt(53156.25) * y * y * z * z * z) + e_3 * (std::sqrt(212625.0) * x * x * z - std::sqrt(212625.0) * y * y * z);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        pc_16[k] = e_0 * (std::sqrt(1.4765625) * x * x * x * x * x * x * x * y * y + std::sqrt(1.4765625) * x * x * x * x * x * y * y * y * y - std::sqrt(289.40625) * x * x * x * x * x * y * y * z * z - std::sqrt(1.4765625) * x * x * x * y * y * y * y * y * y + std::sqrt(2362.5) * x * x * x * y * y * z * z * z * z - std::sqrt(1.4765625) * x * y * y * y * y * y * y * y * y + std::sqrt(289.40625) * x * y * y * y * y * y * y * z * z - std::sqrt(2362.5) * x * y * y * y * y * z * z * z * z) + e_1 * (std::sqrt(1.4765625) * x * x * x * x * x * x * x + std::sqrt(249.5390625) * x * x * x * x * x * y * y - std::sqrt(289.40625) * x * x * x * x * x * z * z - std::sqrt(36.9140625) * x * x * x * y * y * y * y + std::sqrt(14765.625) * x * x * x * y * y * z * z + std::sqrt(2362.5) * x * x * x * z * z * z * z - std::sqrt(426.7265625) * x * y * y * y * y * y * y - std::sqrt(1328.90625) * x * y * y * y * y * z * z - std::sqrt(21262.5) * x * y * y * z * z * z * z) + e_2 * (std::sqrt(147.65625) * x * x * x * x * x + std::sqrt(14765.625) * x * x * x * y * y + std::sqrt(14765.625) * x * x * x * z * z - std::sqrt(33222.65625) * x * y * y * y * y - std::sqrt(132890.625) * x * y * y * z * z) + e_3 * (std::sqrt(14765.625) * x * x * x - std::sqrt(132890.625) * x * y * y);

        pc_17[k] = e_0 * (std::sqrt(8.859375) * x * x * x * x * x * x * x * y * z + std::sqrt(8.859375) * x * x * x * x * x * y * y * y * z - std::sqrt(1008.0) * x * x * x * x * x * y * z * z * z - std::sqrt(8.859375) * x * x * x * y * y * y * y * y * z + std::sqrt(393.75) * x * x * x * y * z * z * z * z * z - std::sqrt(8.859375) * x * y * y * y * y * y * y * y * z + std::sqrt(1008.0) * x * y * y * y * y * y * z * z * z - std::sqrt(393.75) * x * y * y * y * z * z * z * z * z) + e_1 * (-std::sqrt(567.0) * x * x * x * x * x * y * z - std::sqrt(14175.0) * x * x * x * y * z * z * z + std::sqrt(567.0) * x * y * y * y * y * y * z + std::sqrt(14175.0) * x * y * y * y * z * z * z) + e_2 * (-std::sqrt(88593.75) * x * x * x * y * z + std::sqrt(88593.75) * x * y * y * y * z);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        pc_18[k] = e_0 * (std::sqrt(1.4765625) * x * x * x * x * x * x * x * x * y + std::sqrt(1.4765625) * x * x * x * x * x * x * y * y * y - std::sqrt(289.40625) * x * x * x * x * x * x * y * z * z - std::sqrt(1.4765625) * x * x * x * x * y * y * y * y * y + std::sqrt(2362.5) * x * x * x * x * y * z * z * z * z - std::sqrt(1.4765625) * x * x * y * y * y * y * y * y * y + std::sqrt(289.40625) * x * x * y * y * y * y * y * z * z - std::sqrt(2362.5) * x * x * y * y * y * z * z * z * z) + e_1 * (std::sqrt(426.7265625) * x * x * x * x * x * x * y + std::sqrt(36.9140625) * x * x * x * x * y * y * y + std::sqrt(1328.90625) * x * x * x * x * y * z * z - std::sqrt(249.5390625) * x * x * y * y * y * y * y - std::sqrt(14765.625) * x * x * y * y * y * z * z + std::sqrt(21262.5) * x * x * y * z * z * z * z - std::sqrt(1.4765625) * y * y * y * y * y * y * y + std::sqrt(289.40625) * y * y * y * y * y * z * z - std::sqrt(2362.5) * y * y * y * z * z * z * z) + e_2 * (std::sqrt(33222.65625) * x * x * x * x * y - std::sqrt(14765.625) * x * x * y * y * y + std::sqrt(132890.625) * x * x * y * z * z - std::sqrt(147.65625) * y * y * y * y * y - std::sqrt(14765.625) * y * y * y * z * z) + e_3 * (std::sqrt(132890.625) * x * x * y - std::sqrt(14765.625) * y * y * y);

        pc_19[k] = e_0 * (-std::sqrt(14.765625) * x * x * x * x * x * x * x * y * z + std::sqrt(14.765625) * x * x * x * x * x * y * y * y * z + std::sqrt(1476.5625) * x * x * x * x * x * y * z * z * z + std::sqrt(14.765625) * x * x * x * y * y * y * y * y * z - std::sqrt(5906.25) * x * x * x * y * y * y * z * z * z - std::sqrt(14.765625) * x * y * y * y * y * y * y * y * z + std::sqrt(1476.5625) * x * y * y * y * y * y * z * z * z) + e_1 * (std::sqrt(2126.25) * x * x * x * x * x * y * z - std::sqrt(23625.0) * x * x * x * y * y * y * z + std::sqrt(23625.0) * x * x * x * y * z * z * z + std::sqrt(2126.25) * x * y * y * y * y * y * z + std::sqrt(23625.0) * x * y * y * y * z * z * z) + e_2 * (std::sqrt(53156.25) * x * x * x * y * z + std::sqrt(53156.25) * x * y * y * y * z + std::sqrt(212625.0) * x * y * z * z * z) + e_3 * (std::sqrt(850500.0) * x * y * z);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        pc_20[k] = e_0 * (-std::sqrt(2.4609375) * x * x * x * x * x * x * x * x * y + std::sqrt(22.1484375) * x * x * x * x * x * x * y * y * y + std::sqrt(246.09375) * x * x * x * x * x * x * y * z * z + std::sqrt(2.4609375) * x * x * x * x * y * y * y * y * y - std::sqrt(3937.5) * x * x * x * x * y * y * y * z * z - std::sqrt(22.1484375) * x * x * y * y * y * y * y * y * y + std::sqrt(2214.84375) * x * x * y * y * y * y * y * z * z) + e_1 * (-std::sqrt(199.3359375) * x * x * x * x * x * x * y + std::sqrt(553.7109375) * x * x * x * x * y * y * y + std::sqrt(2214.84375) * x * x * x * x * y * z * z - std::sqrt(1794.0234375) * x * x * y * y * y * y * y + std::sqrt(8859.375) * x * x * y * y * y * z * z - std::sqrt(22.1484375) * y * y * y * y * y * y * y + std::sqrt(2214.84375) * y * y * y * y * y * z * z) + e_2 * (-std::sqrt(2214.84375) * x * x * x * x * y - std::sqrt(8859.375) * x * x * y * y * y + std::sqrt(79734.375) * x * x * y * z * z - std::sqrt(2214.84375) * y * y * y * y * y + std::sqrt(79734.375) * y * y * y * z * z) + e_3 * (-std::sqrt(8859.375) * x * x * y - std::sqrt(8859.375) * y * y * y + std::sqrt(141750.0) * y * z * z);

        pc_21[k] = e_0 * (-std::sqrt(373.7548828125) * x * x * x * x * x * x * y * y * z - std::sqrt(41.5283203125) * x * x * x * x * y * y * y * y * z + std::sqrt(2657.8125) * x * x * x * x * y * y * z * z * z + std::sqrt(115.3564453125) * x * x * y * y * y * y * y * y * z - std::sqrt(1181.25) * x * x * y * y * y * y * z * z * z - std::sqrt(4.6142578125) * y * y * y * y * y * y * y * y * z + std::sqrt(32.8125) * y * y * y * y * y * y * z * z * z) + e_1 * (-std::sqrt(373.7548828125) * x * x * x * x * x * x * z - std::sqrt(30274.1455078125) * x * x * x * x * y * y * z + std::sqrt(2657.8125) * x * x * x * x * z * z * z + std::sqrt(373.7548828125) * x * x * y * y * y * y * z + std::sqrt(10631.25) * x * x * y * y * z * z * z - std::sqrt(1038.2080078125) * y * y * y * y * y * y * z + std::sqrt(2657.8125) * y * y * y * y * z * z * z) + e_2 * (-std::sqrt(23920.3125) * x * x * x * x * z - std::sqrt(95681.25) * x * x * y * y * z + std::sqrt(42525.0) * x * x * z * z * z - std::sqrt(23920.3125) * y * y * y * y * z + std::sqrt(42525.0) * y * y * z * z * z) + e_3 * (-std::sqrt(42525.0) * x * x * z - std::sqrt(42525.0) * y * y * z + std::sqrt(18900.0) * z * z * z);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        pc_22[k] = e_0 * (-std::sqrt(996.6796875) * x * x * x * x * x * y * y * z * z - std::sqrt(442.96875) * x * x * x * y * y * y * y * z * z + std::sqrt(7087.5) * x * x * x * y * y * z * z * z * z + std::sqrt(110.7421875) * x * y * y * y * y * y * y * z * z - std::sqrt(787.5) * x * y * y * y * y * z * z * z * z) + e_1 * (-std::sqrt(996.6796875) * x * x * x * x * x * y * y - std::sqrt(996.6796875) * x * x * x * x * x * z * z - std::sqrt(442.96875) * x * x * x * y * y * y * y + std::sqrt(3986.71875) * x * x * x * y * y * z * z + std::sqrt(7087.5) * x * x * x * z * z * z * z + std::sqrt(110.7421875) * x * y * y * y * y * y * y - std::sqrt(5426.3671875) * x * y * y * y * y * z * z + std::sqrt(7087.5) * x * y * y * z * z * z * z) + e_2 * (-std::sqrt(996.6796875) * x * x * x * x * x - std::sqrt(35880.46875) * x * x * x * y * y + std::sqrt(15946.875) * x * x * x * z * z + std::sqrt(110.7421875) * x * y * y * y * y + std::sqrt(15946.875) * x * y * y * z * z + std::sqrt(28350.0) * x * z * z * z * z) + e_3 * (-std::sqrt(15946.875) * x * x * x - std::sqrt(15946.875) * x * y * y + std::sqrt(255150.0) * x * z * z);

        pc_23[k] = e_0 * (std::sqrt(24.9169921875) * x * x * x * x * x * x * y * y * z + std::sqrt(69.2138671875) * x * x * x * x * y * y * y * y * z - std::sqrt(1107.421875) * x * x * x * x * y * y * z * z * z + std::sqrt(2.7685546875) * x * x * y * y * y * y * y * y * z - std::sqrt(492.1875) * x * x * y * y * y * y * z * z * z + std::sqrt(2835.0) * x * x * y * y * z * z * z * z * z - std::sqrt(2.7685546875) * y * y * y * y * y * y * y * y * z + std::sqrt(123.046875) * y * y * y * y * y * y * z * z * z - std::sqrt(315.0) * y * y * y * y * z * z * z * z * z) + e_1 * (std::sqrt(24.9169921875) * x * x * x * x * x * x * z + std::sqrt(622.9248046875) * x * x * x * x * y * y * z - std::sqrt(1107.421875) * x * x * x * x * z * z * z + std::sqrt(69.2138671875) * x * x * y * y * y * y * z + std::sqrt(39867.1875) * x * x * y * y * z * z * z + std::sqrt(2835.0) * x * x * z * z * z * z * z - std::sqrt(135.6591796875) * y * y * y * y * y * y * z - std::sqrt(1107.421875) * y * y * y * y * z * z * z - std::sqrt(2835.0) * y * y * z * z * z * z * z) + e_2 * (std::sqrt(159468.75) * x * x * y * y * z + std::sqrt(70875.0) * x * x * z * z * z - std::sqrt(17718.75) * y * y * y * y * z - std::sqrt(70875.0) * y * y * z * z * z) + e_3 * (std::sqrt(159468.75) * x * x * z - std::sqrt(159468.75) * y * y * z);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        pc_24[k] = e_0 * (std::sqrt(149.501953125) * x * x * x * x * x * x * y * z * z + std::sqrt(415.283203125) * x * x * x * x * y * y * y * z * z - std::sqrt(1661.1328125) * x * x * x * x * y * z * z * z * z + std::sqrt(16.611328125) * x * x * y * y * y * y * y * z * z - std::sqrt(738.28125) * x * x * y * y * y * z * z * z * z + std::sqrt(472.5) * x * x * y * z * z * z * z * z * z - std::sqrt(16.611328125) * y * y * y * y * y * y * y * z * z + std::sqrt(184.5703125) * y * y * y * y * y * z * z * z * z - std::sqrt(52.5) * y * y * y * z * z * z * z * z * z) + e_1 * (std::sqrt(149.501953125) * x * x * x * x * x * x * y + std::sqrt(415.283203125) * x * x * x * x * y * y * y + std::sqrt(16.611328125) * x * x * y * y * y * y * y - std::sqrt(16.611328125) * y * y * y * y * y * y * y) + e_2 * (std::sqrt(14950.1953125) * x * x * x * x * y + std::sqrt(6644.53125) * x * x * y * y * y - std::sqrt(1661.1328125) * y * y * y * y * y) + e_3 * (std::sqrt(106312.5) * x * x * y - std::sqrt(11812.5) * y * y * y);

        pc_25[k] = e_0 * (std::sqrt(24.9169921875) * x * x * x * x * x * x * x * y * z + std::sqrt(69.2138671875) * x * x * x * x * x * y * y * y * z - std::sqrt(1107.421875) * x * x * x * x * x * y * z * z * z + std::sqrt(2.7685546875) * x * x * x * y * y * y * y * y * z - std::sqrt(492.1875) * x * x * x * y * y * y * z * z * z + std::sqrt(2835.0) * x * x * x * y * z * z * z * z * z - std::sqrt(2.7685546875) * x * y * y * y * y * y * y * y * z + std::sqrt(123.046875) * x * y * y * y * y * y * z * z * z - std::sqrt(315.0) * x * y * y * y * z * z * z * z * z) + e_1 * (std::sqrt(897.01171875) * x * x * x * x * x * y * z + std::sqrt(1107.421875) * x * x * x * y * y * y * z + std::sqrt(17718.75) * x * x * x * y * z * z * z + std::sqrt(11.07421875) * x * y * y * y * y * y * z - std::sqrt(17718.75) * x * y * y * y * z * z * z + std::sqrt(11340.0) * x * y * z * z * z * z * z) + e_2 * (std::sqrt(159468.75) * x * x * x * y * z - std::sqrt(17718.75) * x * y * y * y * z + std::sqrt(283500.0) * x * y * z * z * z) + e_3 * (std::sqrt(637875.0) * x * y * z);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        pc_26[k] = e_0 * (-std::sqrt(249.169921875) * x * x * x * x * x * x * y * z * z + std::sqrt(27.685546875) * x * x * x * x * y * y * y * z * z + std::sqrt(1771.875) * x * x * x * x * y * z * z * z * z + std::sqrt(249.169921875) * x * x * y * y * y * y * y * z * z - std::sqrt(3150.0) * x * x * y * y * y * z * z * z * z - std::sqrt(27.685546875) * y * y * y * y * y * y * y * z * z + std::sqrt(196.875) * y * y * y * y * y * z * z * z * z) + e_1 * (-std::sqrt(249.169921875) * x * x * x * x * x * x * y + std::sqrt(27.685546875) * x * x * x * x * y * y * y + std::sqrt(996.6796875) * x * x * x * x * y * z * z + std::sqrt(249.169921875) * x * x * y * y * y * y * y - std::sqrt(21705.46875) * x * x * y * y * y * z * z + std::sqrt(7087.5) * x * x * y * z * z * z * z - std::sqrt(27.685546875) * y * y * y * y * y * y * y - std::sqrt(110.7421875) * y * y * y * y * y * z * z + std::sqrt(7087.5) * y * y * y * z * z * z * z) + e_2 * (-std::sqrt(8970.1171875) * x * x * x * x * y + std::sqrt(442.96875) * x * x * y * y * y + std::sqrt(15946.875) * x * x * y * z * z - std::sqrt(2768.5546875) * y * y * y * y * y + std::sqrt(15946.875) * y * y * y * z * z + std::sqrt(28350.0) * y * z * z * z * z) + e_3 * (-std::sqrt(15946.875) * x * x * y - std::sqrt(15946.875) * y * y * y + std::sqrt(255150.0) * y * z * z);

        pc_27[k] = e_0 * (-std::sqrt(41.5283203125) * x * x * x * x * x * x * x * y * z + std::sqrt(226.0986328125) * x * x * x * x * x * y * y * y * z + std::sqrt(295.3125) * x * x * x * x * x * y * z * z * z + std::sqrt(226.0986328125) * x * x * x * y * y * y * y * y * z - std::sqrt(3281.25) * x * x * x * y * y * y * z * z * z - std::sqrt(41.5283203125) * x * y * y * y * y * y * y * y * z + std::sqrt(295.3125) * x * y * y * y * y * y * z * z * z) + e_1 * (-std::sqrt(1495.01953125) * x * x * x * x * x * y * z + std::sqrt(16611.328125) * x * x * x * y * y * y * z - std::sqrt(1495.01953125) * x * y * y * y * y * y * z);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        pc_28[k] = e_0 * (std::sqrt(4.6142578125) * x * x * x * x * x * x * x * y * y + std::sqrt(12.8173828125) * x * x * x * x * x * y * y * y * y - std::sqrt(1181.25) * x * x * x * x * x * y * y * z * z + std::sqrt(0.5126953125) * x * x * x * y * y * y * y * y * y - std::sqrt(525.0) * x * x * x * y * y * y * y * z * z + std::sqrt(1181.25) * x * x * x * y * y * z * z * z * z - std::sqrt(0.5126953125) * x * y * y * y * y * y * y * y * y + std::sqrt(131.25) * x * y * y * y * y * y * y * z * z - std::sqrt(131.25) * x * y * y * y * y * z * z * z * z) + e_1 * (std::sqrt(4.6142578125) * x * x * x * x * x * x * x + std::sqrt(1038.2080078125) * x * x * x * x * x * y * y - std::sqrt(1181.25) * x * x * x * x * x * z * z + std::sqrt(558.3251953125) * x * x * x * y * y * y * y - std::sqrt(75600.0) * x * x * x * y * y * z * z + std::sqrt(1181.25) * x * x * x * z * z * z * z - std::sqrt(41.5283203125) * x * y * y * y * y * y * y + std::sqrt(1181.25) * x * y * y * y * y * z * z + std::sqrt(1181.25) * x * y * y * z * z * z * z) + e_2 * (std::sqrt(461.42578125) * x * x * x * x * x + std::sqrt(8933.203125) * x * x * x * y * y - std::sqrt(42525.0) * x * x * x * z * z + std::sqrt(18.45703125) * x * y * y * y * y - std::sqrt(42525.0) * x * y * y * z * z + std::sqrt(4725.0) * x * z * z * z * z) + e_3 * (std::sqrt(1181.25) * x * x * x + std::sqrt(1181.25) * x * y * y - std::sqrt(18900.0) * x * z * z);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        pc_29[k] = e_0 * (std::sqrt(12.3046875) * x * x * x * x * x * x * y * y * z + std::sqrt(49.21875) * x * x * x * x * y * y * y * y * z - std::sqrt(3150.0) * x * x * x * x * y * y * z * z * z + std::sqrt(12.3046875) * x * x * y * y * y * y * y * y * z - std::sqrt(3150.0) * x * x * y * y * y * y * z * z * z + std::sqrt(3150.0) * x * x * y * y * z * z * z * z * z) + e_1 * (std::sqrt(12.3046875) * x * x * x * x * x * x * z - std::sqrt(5426.3671875) * x * x * x * x * y * y * z - std::sqrt(3150.0) * x * x * x * x * z * z * z - std::sqrt(5426.3671875) * x * x * y * y * y * y * z - std::sqrt(12600.0) * x * x * y * y * z * z * z + std::sqrt(3150.0) * x * x * z * z * z * z * z + std::sqrt(12.3046875) * y * y * y * y * y * y * z - std::sqrt(3150.0) * y * y * y * y * z * z * z + std::sqrt(3150.0) * y * y * z * z * z * z * z) + e_2 * (-std::sqrt(8970.1171875) * x * x * x * x * z - std::sqrt(372536.71875) * x * x * y * y * z + std::sqrt(3150.0) * x * x * z * z * z - std::sqrt(8970.1171875) * y * y * y * y * z + std::sqrt(3150.0) * y * y * z * z * z + std::sqrt(3150.0) * z * z * z * z * z) + e_3 * (-std::sqrt(113400.0) * x * x * z - std::sqrt(113400.0) * y * y * z + std::sqrt(50400.0) * z * z * z);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        pc_30[k] = e_0 * (-std::sqrt(0.3076171875) * x * x * x * x * x * x * x * y * y - std::sqrt(2.7685546875) * x * x * x * x * x * y * y * y * y + std::sqrt(123.046875) * x * x * x * x * x * y * y * z * z - std::sqrt(2.7685546875) * x * x * x * y * y * y * y * y * y + std::sqrt(492.1875) * x * x * x * y * y * y * y * z * z - std::sqrt(1968.75) * x * x * x * y * y * z * z * z * z - std::sqrt(0.3076171875) * x * y * y * y * y * y * y * y * y + std::sqrt(123.046875) * x * y * y * y * y * y * y * z * z - std::sqrt(1968.75) * x * y * y * y * y * z * z * z * z + std::sqrt(1260.0) * x * y * y * z * z * z * z * z * z) + e_1 * (-std::sqrt(0.3076171875) * x * x * x * x * x * x * x - std::sqrt(111.0498046875) * x * x * x * x * x * y * y + std::sqrt(123.046875) * x * x * x * x * x * z * z - std::sqrt(376.8310546875) * x * x * x * y * y * y * y - std::sqrt(492.1875) * x * x * x * y * y * z * z - std::sqrt(1968.75) * x * x * x * z * z * z * z - std::sqrt(88.9013671875) * x * y * y * y * y * y * y - std::sqrt(1107.421875) * x * y * y * y * y * z * z + std::sqrt(17718.75) * x * y * y * z * z * z * z + std::sqrt(1260.0) * x * z * z * z * z * z * z) + e_2 * (-std::sqrt(30.76171875) * x * x * x * x * x - std::sqrt(14888.671875) * x * x * x * y * y - std::sqrt(7875.0) * x * x * x * z * z - std::sqrt(13565.91796875) * x * y * y * y * y + std::sqrt(70875.0) * x * y * y * z * z + std::sqrt(70875.0) * x * z * z * z * z) + e_3 * (-std::sqrt(7875.0) * x * x * x - std::sqrt(70875.0) * x * y * y + std::sqrt(283500.0) * x * z * z);

        pc_31[k] = e_0 * (-std::sqrt(1.845703125) * x * x * x * x * x * x * x * y * z - std::sqrt(16.611328125) * x * x * x * x * x * y * y * y * z + std::sqrt(512.6953125) * x * x * x * x * x * y * z * z * z - std::sqrt(16.611328125) * x * x * x * y * y * y * y * y * z + std::sqrt(2050.78125) * x * x * x * y * y * y * z * z * z - std::sqrt(1312.5) * x * x * x * y * z * z * z * z * z - std::sqrt(1.845703125) * x * y * y * y * y * y * y * y * z + std::sqrt(512.6953125) * x * y * y * y * y * y * z * z * z - std::sqrt(1312.5) * x * y * y * y * z * z * z * z * z + std::sqrt(210.0) * x * y * z * z * z * z * z * z * z) + e_1 * (std::sqrt(738.28125) * x * x * x * x * x * y * z + std::sqrt(2953.125) * x * x * x * y * y * y * z + std::sqrt(738.28125) * x * y * y * y * y * y * z + std::sqrt(7560.0) * x * y * z * z * z * z * z) + e_2 * (std::sqrt(47250.0) * x * x * x * y * z + std::sqrt(47250.0) * x * y * y * y * z + std::sqrt(189000.0) * x * y * z * z * z) + e_3 * (std::sqrt(756000.0) * x * y * z);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        pc_32[k] = e_0 * (-std::sqrt(0.3076171875) * x * x * x * x * x * x * x * x * y - std::sqrt(2.7685546875) * x * x * x * x * x * x * y * y * y + std::sqrt(123.046875) * x * x * x * x * x * x * y * z * z - std::sqrt(2.7685546875) * x * x * x * x * y * y * y * y * y + std::sqrt(492.1875) * x * x * x * x * y * y * y * z * z - std::sqrt(1968.75) * x * x * x * x * y * z * z * z * z - std::sqrt(0.3076171875) * x * x * y * y * y * y * y * y * y + std::sqrt(123.046875) * x * x * y * y * y * y * y * z * z - std::sqrt(1968.75) * x * x * y * y * y * z * z * z * z + std::sqrt(1260.0) * x * x * y * z * z * z * z * z * z) + e_1 * (-std::sqrt(88.9013671875) * x * x * x * x * x * x * y - std::sqrt(376.8310546875) * x * x * x * x * y * y * y - std::sqrt(1107.421875) * x * x * x * x * y * z * z - std::sqrt(111.0498046875) * x * x * y * y * y * y * y - std::sqrt(492.1875) * x * x * y * y * y * z * z + std::sqrt(17718.75) * x * x * y * z * z * z * z - std::sqrt(0.3076171875) * y * y * y * y * y * y * y + std::sqrt(123.046875) * y * y * y * y * y * z * z - std::sqrt(1968.75) * y * y * y * z * z * z * z + std::sqrt(1260.0) * y * z * z * z * z * z * z) + e_2 * (-std::sqrt(13565.91796875) * x * x * x * x * y - std::sqrt(14888.671875) * x * x * y * y * y + std::sqrt(70875.0) * x * x * y * z * z - std::sqrt(30.76171875) * y * y * y * y * y - std::sqrt(7875.0) * y * y * y * z * z + std::sqrt(70875.0) * y * z * z * z * z) + e_3 * (-std::sqrt(70875.0) * x * x * y - std::sqrt(7875.0) * y * y * y + std::sqrt(283500.0) * y * z * z);

        pc_33[k] = e_0 * (std::sqrt(3.076171875) * x * x * x * x * x * x * x * y * z + std::sqrt(3.076171875) * x * x * x * x * x * y * y * y * z - std::sqrt(787.5) * x * x * x * x * x * y * z * z * z - std::sqrt(3.076171875) * x * x * x * y * y * y * y * y * z + std::sqrt(787.5) * x * x * x * y * z * z * z * z * z - std::sqrt(3.076171875) * x * y * y * y * y * y * y * y * z + std::sqrt(787.5) * x * y * y * y * y * y * z * z * z - std::sqrt(787.5) * x * y * y * y * z * z * z * z * z) + e_1 * (-std::sqrt(1771.875) * x * x * x * x * x * y * z + std::sqrt(1771.875) * x * y * y * y * y * y * z) + e_2 * (-std::sqrt(44296.875) * x * x * x * y * z + std::sqrt(44296.875) * x * y * y * y * z);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        pc_34[k] = e_0 * (std::sqrt(0.5126953125) * x * x * x * x * x * x * x * x * y - std::sqrt(0.5126953125) * x * x * x * x * x * x * y * y * y - std::sqrt(131.25) * x * x * x * x * x * x * y * z * z - std::sqrt(12.8173828125) * x * x * x * x * y * y * y * y * y + std::sqrt(525.0) * x * x * x * x * y * y * y * z * z + std::sqrt(131.25) * x * x * x * x * y * z * z * z * z - std::sqrt(4.6142578125) * x * x * y * y * y * y * y * y * y + std::sqrt(1181.25) * x * x * y * y * y * y * y * z * z - std::sqrt(1181.25) * x * x * y * y * y * z * z * z * z) + e_1 * (std::sqrt(41.5283203125) * x * x * x * x * x * x * y - std::sqrt(558.3251953125) * x * x * x * x * y * y * y - std::sqrt(1181.25) * x * x * x * x * y * z * z - std::sqrt(1038.2080078125) * x * x * y * y * y * y * y + std::sqrt(75600.0) * x * x * y * y * y * z * z - std::sqrt(1181.25) * x * x * y * z * z * z * z - std::sqrt(4.6142578125) * y * y * y * y * y * y * y + std::sqrt(1181.25) * y * y * y * y * y * z * z - std::sqrt(1181.25) * y * y * y * z * z * z * z) + e_2 * (-std::sqrt(18.45703125) * x * x * x * x * y - std::sqrt(8933.203125) * x * x * y * y * y + std::sqrt(42525.0) * x * x * y * z * z - std::sqrt(461.42578125) * y * y * y * y * y + std::sqrt(42525.0) * y * y * y * z * z - std::sqrt(4725.0) * y * z * z * z * z) + e_3 * (-std::sqrt(1181.25) * x * x * y - std::sqrt(1181.25) * y * y * y + std::sqrt(18900.0) * y * z * z);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        pc_35[k] = e_0 * (std::sqrt(46.142578125) * x * x * x * x * x * x * y * y * z + std::sqrt(128.173828125) * x * x * x * x * y * y * y * y * z - std::sqrt(738.28125) * x * x * x * x * y * y * z * z * z + std::sqrt(5.126953125) * x * x * y * y * y * y * y * y * z - std::sqrt(328.125) * x * x * y * y * y * y * z * z * z + std::sqrt(118.125) * x * x * y * y * z * z * z * z * z - std::sqrt(5.126953125) * y * y * y * y * y * y * y * y * z + std::sqrt(82.03125) * y * y * y * y * y * y * z * z * z - std::sqrt(13.125) * y * y * y * y * z * z * z * z * z) + e_1 * (std::sqrt(46.142578125) * x * x * x * x * x * x * z + std::sqrt(7798.095703125) * x * x * x * x * y * y * z - std::sqrt(738.28125) * x * x * x * x * z * z * z + std::sqrt(2260.986328125) * x * x * y * y * y * y * z - std::sqrt(26578.125) * x * x * y * y * z * z * z + std::sqrt(118.125) * x * x * z * z * z * z * z - std::sqrt(1153.564453125) * y * y * y * y * y * y * z + std::sqrt(6644.53125) * y * y * y * y * z * z * z - std::sqrt(118.125) * y * y * z * z * z * z * z) + e_2 * (std::sqrt(2953.125) * x * x * x * x * z + std::sqrt(26578.125) * x * x * y * y * z - std::sqrt(11812.5) * x * x * z * z * z - std::sqrt(11812.5) * y * y * y * y * z + std::sqrt(11812.5) * y * y * z * z * z) + e_3 * (std::sqrt(2953.125) * x * x * z - std::sqrt(2953.125) * y * y * z);

        pc_36[k] = e_0 * (std::sqrt(123.046875) * x * x * x * x * x * y * y * z * z + std::sqrt(492.1875) * x * x * x * y * y * y * y * z * z - std::sqrt(1968.75) * x * x * x * y * y * z * z * z * z + std::sqrt(123.046875) * x * y * y * y * y * y * y * z * z - std::sqrt(1968.75) * x * y * y * y * y * z * z * z * z + std::sqrt(315.0) * x * y * y * z * z * z * z * z * z) + e_1 * (std::sqrt(123.046875) * x * x * x * x * x * y * y + std::sqrt(123.046875) * x * x * x * x * x * z * z + std::sqrt(492.1875) * x * x * x * y * y * y * y - std::sqrt(492.1875) * x * x * x * y * y * z * z - std::sqrt(1968.75) * x * x * x * z * z * z * z + std::sqrt(123.046875) * x * y * y * y * y * y * y - std::sqrt(1107.421875) * x * y * y * y * y * z * z - std::sqrt(17718.75) * x * y * y * z * z * z * z + std::sqrt(315.0) * x * z * z * z * z * z * z) + e_2 * (std::sqrt(123.046875) * x * x * x * x * x + std::sqrt(12304.6875) * x * x * x * y * y - std::sqrt(7875.0) * x * x * x * z * z + std::sqrt(9966.796875) * x * y * y * y * y - std::sqrt(283500.0) * x * y * y * z * z) + e_3 * (std::sqrt(1968.75) * x * x * x + std::sqrt(17718.75) * x * y * y - std::sqrt(70875.0) * x * z * z);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        pc_37[k] = e_0 * (-std::sqrt(3.076171875) * x * x * x * x * x * x * y * y * z - std::sqrt(27.685546875) * x * x * x * x * y * y * y * y * z + std::sqrt(196.875) * x * x * x * x * y * y * z * z * z - std::sqrt(27.685546875) * x * x * y * y * y * y * y * y * z + std::sqrt(787.5) * x * x * y * y * y * y * z * z * z - std::sqrt(952.875) * x * x * y * y * z * z * z * z * z - std::sqrt(3.076171875) * y * y * y * y * y * y * y * y * z + std::sqrt(196.875) * y * y * y * y * y * y * z * z * z - std::sqrt(952.875) * y * y * y * y * z * z * z * z * z + std::sqrt(126.0) * y * y * z * z * z * z * z * z * z) + e_1 * (-std::sqrt(3.076171875) * x * x * x * x * x * x * z - std::sqrt(249.169921875) * x * x * x * x * y * y * z + std::sqrt(196.875) * x * x * x * x * z * z * z - std::sqrt(692.138671875) * x * x * y * y * y * y * z - std::sqrt(3150.0) * x * x * y * y * z * z * z - std::sqrt(952.875) * x * x * z * z * z * z * z - std::sqrt(150.732421875) * y * y * y * y * y * y * z - std::sqrt(4921.875) * y * y * y * y * z * z * z + std::sqrt(385.875) * y * y * z * z * z * z * z + std::sqrt(126.0) * z * z * z * z * z * z * z) + e_2 * (-std::sqrt(44296.875) * x * x * y * y * z - std::sqrt(19687.5) * x * x * z * z * z - std::sqrt(44296.875) * y * y * y * y * z - std::sqrt(19687.5) * y * y * z * z * z + std::sqrt(12600.0) * z * z * z * z * z) + e_3 * (-std::sqrt(44296.875) * x * x * z - std::sqrt(398671.875) * y * y * z + std::sqrt(78750.0) * z * z * z);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        pc_38[k] = e_0 * (-std::sqrt(18.45703125) * x * x * x * x * x * x * y * z * z - std::sqrt(166.11328125) * x * x * x * x * y * y * y * z * z + std::sqrt(401.953125) * x * x * x * x * y * z * z * z * z - std::sqrt(166.11328125) * x * x * y * y * y * y * y * z * z + std::sqrt(1607.8125) * x * x * y * y * y * z * z * z * z - std::sqrt(336.0) * x * x * y * z * z * z * z * z * z - std::sqrt(18.45703125) * y * y * y * y * y * y * y * z * z + std::sqrt(401.953125) * y * y * y * y * y * z * z * z * z - std::sqrt(336.0) * y * y * y * z * z * z * z * z * z + std::sqrt(21.0) * y * z * z * z * z * z * z * z * z) + e_1 * (-std::sqrt(18.45703125) * x * x * x * x * x * x * y - std::sqrt(166.11328125) * x * x * x * x * y * y * y + std::sqrt(295.3125) * x * x * x * x * y * z * z - std::sqrt(166.11328125) * x * x * y * y * y * y * y + std::sqrt(1181.25) * x * x * y * y * y * z * z - std::sqrt(1181.25) * x * x * y * z * z * z * z - std::sqrt(18.45703125) * y * y * y * y * y * y * y + std::sqrt(295.3125) * y * y * y * y * y * z * z - std::sqrt(1181.25) * y * y * y * z * z * z * z + std::sqrt(3024.0) * y * z * z * z * z * z * z) + e_2 * (-std::sqrt(1845.703125) * x * x * x * x * y - std::sqrt(7382.8125) * x * x * y * y * y - std::sqrt(1845.703125) * y * y * y * y * y + std::sqrt(118125.0) * y * z * z * z * z) + e_3 * (-std::sqrt(29531.25) * x * x * y - std::sqrt(29531.25) * y * y * y + std::sqrt(472500.0) * y * z * z);

        pc_39[k] = e_0 * (-std::sqrt(3.076171875) * x * x * x * x * x * x * x * y * z - std::sqrt(27.685546875) * x * x * x * x * x * y * y * y * z + std::sqrt(196.875) * x * x * x * x * x * y * z * z * z - std::sqrt(27.685546875) * x * x * x * y * y * y * y * y * z + std::sqrt(787.5) * x * x * x * y * y * y * z * z * z - std::sqrt(952.875) * x * x * x * y * z * z * z * z * z - std::sqrt(3.076171875) * x * y * y * y * y * y * y * y * z + std::sqrt(196.875) * x * y * y * y * y * y * z * z * z - std::sqrt(952.875) * x * y * y * y * z * z * z * z * z + std::sqrt(126.0) * x * y * z * z * z * z * z * z * z) + e_1 * (-std::sqrt(110.7421875) * x * x * x * x * x * y * z - std::sqrt(442.96875) * x * x * x * y * y * y * z - std::sqrt(7087.5) * x * x * x * y * z * z * z - std::sqrt(110.7421875) * x * y * y * y * y * y * z - std::sqrt(7087.5) * x * y * y * y * z * z * z + std::sqrt(2551.5) * x * y * z * z * z * z * z) + e_2 * (-std::sqrt(44296.875) * x * x * x * y * z - std::sqrt(44296.875) * x * y * y * y * z) + e_3 * (-std::sqrt(177187.5) * x * y * z);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        pc_40[k] = e_0 * (std::sqrt(30.76171875) * x * x * x * x * x * x * y * z * z + std::sqrt(30.76171875) * x * x * x * x * y * y * y * z * z - std::sqrt(492.1875) * x * x * x * x * y * z * z * z * z - std::sqrt(30.76171875) * x * x * y * y * y * y * y * z * z + std::sqrt(78.75) * x * x * y * z * z * z * z * z * z - std::sqrt(30.76171875) * y * y * y * y * y * y * y * z * z + std::sqrt(492.1875) * y * y * y * y * y * z * z * z * z - std::sqrt(78.75) * y * y * y * z * z * z * z * z * z) + e_1 * (std::sqrt(30.76171875) * x * x * x * x * x * x * y + std::sqrt(30.76171875) * x * x * x * x * y * y * y - std::sqrt(1107.421875) * x * x * x * x * y * z * z - std::sqrt(30.76171875) * x * x * y * y * y * y * y - std::sqrt(492.1875) * x * x * y * y * y * z * z - std::sqrt(30.76171875) * y * y * y * y * y * y * y + std::sqrt(123.046875) * y * y * y * y * y * z * z + std::sqrt(7875.0) * y * y * y * z * z * z * z - std::sqrt(315.0) * y * z * z * z * z * z * z) + e_2 * (std::sqrt(1107.421875) * x * x * x * x * y - std::sqrt(492.1875) * x * x * y * y * y - std::sqrt(17718.75) * x * x * y * z * z - std::sqrt(3076.171875) * y * y * y * y * y + std::sqrt(96468.75) * y * y * y * z * z) + e_3 * (-std::sqrt(7875.0) * y * y * y + std::sqrt(70875.0) * y * z * z);

        pc_41[k] = e_0 * (std::sqrt(5.126953125) * x * x * x * x * x * x * x * y * z - std::sqrt(5.126953125) * x * x * x * x * x * y * y * y * z - std::sqrt(82.03125) * x * x * x * x * x * y * z * z * z - std::sqrt(128.173828125) * x * x * x * y * y * y * y * y * z + std::sqrt(328.125) * x * x * x * y * y * y * z * z * z + std::sqrt(13.125) * x * x * x * y * z * z * z * z * z - std::sqrt(46.142578125) * x * y * y * y * y * y * y * y * z + std::sqrt(738.28125) * x * y * y * y * y * y * z * z * z - std::sqrt(118.125) * x * y * y * y * z * z * z * z * z) + e_1 * (std::sqrt(184.5703125) * x * x * x * x * x * y * z - std::sqrt(6644.53125) * x * x * x * y * y * y * z - std::sqrt(9043.9453125) * x * y * y * y * y * y * z + std::sqrt(47250.0) * x * y * y * y * z * z * z - std::sqrt(472.5) * x * y * z * z * z * z * z) + e_2 * (-std::sqrt(2953.125) * x * x * x * y * z - std::sqrt(73828.125) * x * y * y * y * z + std::sqrt(47250.0) * x * y * z * z * z) + e_3 * (-std::sqrt(11812.5) * x * y * z);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        pc_42[k] = e_0 * (-std::sqrt(0.54931640625) * x * x * x * x * x * x * x * x * y - std::sqrt(3.90625) * x * x * x * x * x * x * y * y * y + std::sqrt(177.978515625) * x * x * x * x * x * x * y * z * z - std::sqrt(2.197265625) * x * x * x * x * y * y * y * y * y + std::sqrt(494.384765625) * x * x * x * x * y * y * y * z * z - std::sqrt(316.40625) * x * x * x * x * y * z * z * z * z + std::sqrt(19.775390625) * x * x * y * y * y * y * y * z * z - std::sqrt(140.625) * x * x * y * y * y * z * z * z * z + std::sqrt(5.625) * x * x * y * z * z * z * z * z * z + std::sqrt(0.06103515625) * y * y * y * y * y * y * y * y * y - std::sqrt(19.775390625) * y * y * y * y * y * y * y * z * z + std::sqrt(35.15625) * y * y * y * y * y * z * z * z * z - std::sqrt(0.625) * y * y * y * z * z * z * z * z * z) + e_1 * (-std::sqrt(177.978515625) * x * x * x * x * x * x * y - std::sqrt(494.384765625) * x * x * x * x * y * y * y + std::sqrt(25628.90625) * x * x * x * x * y * z * z - std::sqrt(19.775390625) * x * x * y * y * y * y * y + std::sqrt(11390.625) * x * x * y * y * y * z * z - std::sqrt(11390.625) * x * x * y * z * z * z * z + std::sqrt(19.775390625) * y * y * y * y * y * y * y - std::sqrt(2847.65625) * y * y * y * y * y * z * z + std::sqrt(1265.625) * y * y * y * z * z * z * z) + e_2 * (-std::sqrt(2847.65625) * x * x * x * x * y - std::sqrt(1265.625) * x * x * y * y * y + std::sqrt(102515.625) * x * x * y * z * z + std::sqrt(316.40625) * y * y * y * y * y - std::sqrt(11390.625) * y * y * y * z * z) + e_3 * (-std::sqrt(1265.625) * x * x * y + std::sqrt(140.625) * y * y * y);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        pc_43[k] = e_0 * (-std::sqrt(1.46484375) * x * x * x * x * x * x * x * y * z - std::sqrt(13.18359375) * x * x * x * x * x * y * y * y * z + std::sqrt(474.609375) * x * x * x * x * x * y * z * z * z - std::sqrt(13.18359375) * x * x * x * y * y * y * y * y * z + std::sqrt(1898.4375) * x * x * x * y * y * y * z * z * z - std::sqrt(843.75) * x * x * x * y * z * z * z * z * z - std::sqrt(1.46484375) * x * y * y * y * y * y * y * y * z + std::sqrt(474.609375) * x * y * y * y * y * y * z * z * z - std::sqrt(843.75) * x * y * y * y * z * z * z * z * z + std::sqrt(15.0) * x * y * z * z * z * z * z * z * z) + e_1 * (std::sqrt(843.75) * x * x * x * x * x * y * z + std::sqrt(3375.0) * x * x * x * y * y * y * z + std::sqrt(3375.0) * x * x * x * y * z * z * z + std::sqrt(843.75) * x * y * y * y * y * y * z + std::sqrt(3375.0) * x * y * y * y * z * z * z - std::sqrt(8640.0) * x * y * z * z * z * z * z) + e_2 * (std::sqrt(102093.75) * x * x * x * y * z + std::sqrt(102093.75) * x * y * y * y * z - std::sqrt(84375.0) * x * y * z * z * z) + e_3 * (std::sqrt(121500.0) * x * y * z);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        pc_44[k] = e_0 * (std::sqrt(0.03662109375) * x * x * x * x * x * x * x * x * y + std::sqrt(0.5859375) * x * x * x * x * x * x * y * y * y - std::sqrt(17.724609375) * x * x * x * x * x * x * y * z * z + std::sqrt(1.318359375) * x * x * x * x * y * y * y * y * y - std::sqrt(159.521484375) * x * x * x * x * y * y * y * z * z + std::sqrt(337.5) * x * x * x * x * y * z * z * z * z + std::sqrt(0.5859375) * x * x * y * y * y * y * y * y * y - std::sqrt(159.521484375) * x * x * y * y * y * y * y * z * z + std::sqrt(1350.0) * x * x * y * y * y * z * z * z * z - std::sqrt(360.375) * x * x * y * z * z * z * z * z * z + std::sqrt(0.03662109375) * y * y * y * y * y * y * y * y * y - std::sqrt(17.724609375) * y * y * y * y * y * y * y * z * z + std::sqrt(337.5) * y * y * y * y * y * z * z * z * z - std::sqrt(360.375) * y * y * y * z * z * z * z * z * z + std::sqrt(6.0) * y * z * z * z * z * z * z * z * z) + e_1 * (std::sqrt(11.865234375) * x * x * x * x * x * x * y + std::sqrt(106.787109375) * x * x * x * x * y * y * y + std::sqrt(84.375) * x * x * x * x * y * z * z + std::sqrt(106.787109375) * x * x * y * y * y * y * y + std::sqrt(337.5) * x * x * y * y * y * z * z - std::sqrt(4134.375) * x * x * y * z * z * z * z + std::sqrt(11.865234375) * y * y * y * y * y * y * y + std::sqrt(84.375) * y * y * y * y * y * z * z - std::sqrt(4134.375) * y * y * y * z * z * z * z - std::sqrt(54.0) * y * z * z * z * z * z * z) + e_2 * (std::sqrt(2109.375) * x * x * x * x * y + std::sqrt(8437.5) * x * x * y * y * y - std::sqrt(18984.375) * x * x * y * z * z + std::sqrt(2109.375) * y * y * y * y * y - std::sqrt(18984.375) * y * y * y * z * z - std::sqrt(33750.0) * y * z * z * z * z) + e_3 * (std::sqrt(18984.375) * x * x * y + std::sqrt(18984.375) * y * y * y - std::sqrt(303750.0) * y * z * z);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        pc_45[k] = e_0 * (0.46875 * x * x * x * x * x * x * x * x * z + 1.875 * x * x * x * x * x * x * y * y * z - 8.75 * x * x * x * x * x * x * z * z * z + 2.8125 * x * x * x * x * y * y * y * y * z - 26.25 * x * x * x * x * y * y * z * z * z + 16.875 * x * x * x * x * z * z * z * z * z + 1.875 * x * x * y * y * y * y * y * y * z - 26.25 * x * x * y * y * y * y * z * z * z + 33.75 * x * x * y * y * z * z * z * z * z - 9.0 * x * x * z * z * z * z * z * z * z + 0.46875 * y * y * y * y * y * y * y * y * z - 8.75 * y * y * y * y * y * y * z * z * z + 16.875 * y * y * y * y * z * z * z * z * z - 9.0 * y * y * z * z * z * z * z * z * z + z * z * z * z * z * z * z * z * z) + e_1 * (-11.25 * x * x * x * x * x * x * z - 33.75 * x * x * x * x * y * y * z + 11.25 * x * x * x * x * z * z * z - 33.75 * x * x * y * y * y * y * z + 22.5 * x * x * y * y * z * z * z - 54.0 * x * x * z * z * z * z * z - 11.25 * y * y * y * y * y * y * z + 11.25 * y * y * y * y * z * z * z - 54.0 * y * y * z * z * z * z * z + 18.0 * z * z * z * z * z * z * z) + e_2 * (-84.375 * x * x * x * x * z - 168.75 * x * x * y * y * z - 225.0 * x * x * z * z * z - 84.375 * y * y * y * y * z - 225.0 * y * y * z * z * z + 135.0 * z * z * z * z * z) + e_3 * (-450.0 * x * x * z - 450.0 * y * y * z + 300.0 * z * z * z);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        pc_46[k] = e_0 * (std::sqrt(0.03662109375) * x * x * x * x * x * x * x * x * x + std::sqrt(0.5859375) * x * x * x * x * x * x * x * y * y - std::sqrt(17.724609375) * x * x * x * x * x * x * x * z * z + std::sqrt(1.318359375) * x * x * x * x * x * y * y * y * y - std::sqrt(159.521484375) * x * x * x * x * x * y * y * z * z + std::sqrt(337.5) * x * x * x * x * x * z * z * z * z + std::sqrt(0.5859375) * x * x * x * y * y * y * y * y * y - std::sqrt(159.521484375) * x * x * x * y * y * y * y * z * z + std::sqrt(1350.0) * x * x * x * y * y * z * z * z * z - std::sqrt(360.375) * x * x * x * z * z * z * z * z * z + std::sqrt(0.03662109375) * x * y * y * y * y * y * y * y * y - std::sqrt(17.724609375) * x * y * y * y * y * y * y * z * z + std::sqrt(337.5) * x * y * y * y * y * z * z * z * z - std::sqrt(360.375) * x * y * y * z * z * z * z * z * z + std::sqrt(6.0) * x * z * z * z * z * z * z * z * z) + e_1 * (std::sqrt(11.865234375) * x * x * x * x * x * x * x + std::sqrt(106.787109375) * x * x * x * x * x * y * y + std::sqrt(84.375) * x * x * x * x * x * z * z + std::sqrt(106.787109375) * x * x * x * y * y * y * y + std::sqrt(337.5) * x * x * x * y * y * z * z - std::sqrt(4134.375) * x * x * x * z * z * z * z + std::sqrt(11.865234375) * x * y * y * y * y * y * y + std::sqrt(84.375) * x * y * y * y * y * z * z - std::sqrt(4134.375) * x * y * y * z * z * z * z - std::sqrt(54.0) * x * z * z * z * z * z * z) + e_2 * (std::sqrt(2109.375) * x * x * x * x * x + std::sqrt(8437.5) * x * x * x * y * y - std::sqrt(18984.375) * x * x * x * z * z + std::sqrt(2109.375) * x * y * y * y * y - std::sqrt(18984.375) * x * y * y * z * z - std::sqrt(33750.0) * x * z * z * z * z) + e_3 * (std::sqrt(18984.375) * x * x * x + std::sqrt(18984.375) * x * y * y - std::sqrt(303750.0) * x * z * z);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        pc_47[k] = e_0 * (-std::sqrt(0.3662109375) * x * x * x * x * x * x * x * x * z - std::sqrt(1.46484375) * x * x * x * x * x * x * y * y * z + std::sqrt(118.65234375) * x * x * x * x * x * x * z * z * z + std::sqrt(118.65234375) * x * x * x * x * y * y * z * z * z - std::sqrt(210.9375) * x * x * x * x * z * z * z * z * z + std::sqrt(1.46484375) * x * x * y * y * y * y * y * y * z - std::sqrt(118.65234375) * x * x * y * y * y * y * z * z * z + std::sqrt(3.75) * x * x * z * z * z * z * z * z * z + std::sqrt(0.3662109375) * y * y * y * y * y * y * y * y * z - std::sqrt(118.65234375) * y * y * y * y * y * y * z * z * z + std::sqrt(210.9375) * y * y * y * y * z * z * z * z * z - std::sqrt(3.75) * y * y * z * z * z * z * z * z * z) + e_1 * (std::sqrt(210.9375) * x * x * x * x * x * x * z + std::sqrt(210.9375) * x * x * x * x * y * y * z + std::sqrt(843.75) * x * x * x * x * z * z * z - std::sqrt(210.9375) * x * x * y * y * y * y * z - std::sqrt(2160.0) * x * x * z * z * z * z * z - std::sqrt(210.9375) * y * y * y * y * y * y * z - std::sqrt(843.75) * y * y * y * y * z * z * z + std::sqrt(2160.0) * y * y * z * z * z * z * z) + e_2 * (std::sqrt(25523.4375) * x * x * x * x * z - std::sqrt(21093.75) * x * x * z * z * z - std::sqrt(25523.4375) * y * y * y * y * z + std::sqrt(21093.75) * y * y * z * z * z) + e_3 * (std::sqrt(30375.0) * x * x * z - std::sqrt(30375.0) * y * y * z);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        pc_48[k] = e_0 * (-std::sqrt(0.06103515625) * x * x * x * x * x * x * x * x * x + std::sqrt(19.775390625) * x * x * x * x * x * x * x * z * z + std::sqrt(2.197265625) * x * x * x * x * x * y * y * y * y - std::sqrt(19.775390625) * x * x * x * x * x * y * y * z * z - std::sqrt(35.15625) * x * x * x * x * x * z * z * z * z + std::sqrt(3.90625) * x * x * x * y * y * y * y * y * y - std::sqrt(494.384765625) * x * x * x * y * y * y * y * z * z + std::sqrt(140.625) * x * x * x * y * y * z * z * z * z + std::sqrt(0.625) * x * x * x * z * z * z * z * z * z + std::sqrt(0.54931640625) * x * y * y * y * y * y * y * y * y - std::sqrt(177.978515625) * x * y * y * y * y * y * y * z * z + std::sqrt(316.40625) * x * y * y * y * y * z * z * z * z - std::sqrt(5.625) * x * y * y * z * z * z * z * z * z) + e_1 * (-std::sqrt(19.775390625) * x * x * x * x * x * x * x + std::sqrt(19.775390625) * x * x * x * x * x * y * y + std::sqrt(2847.65625) * x * x * x * x * x * z * z + std::sqrt(494.384765625) * x * x * x * y * y * y * y - std::sqrt(11390.625) * x * x * x * y * y * z * z - std::sqrt(1265.625) * x * x * x * z * z * z * z + std::sqrt(177.978515625) * x * y * y * y * y * y * y - std::sqrt(25628.90625) * x * y * y * y * y * z * z + std::sqrt(11390.625) * x * y * y * z * z * z * z) + e_2 * (-std::sqrt(316.40625) * x * x * x * x * x + std::sqrt(1265.625) * x * x * x * y * y + std::sqrt(11390.625) * x * x * x * z * z + std::sqrt(2847.65625) * x * y * y * y * y - std::sqrt(102515.625) * x * y * y * z * z) + e_3 * (-std::sqrt(140.625) * x * x * x + std::sqrt(1265.625) * x * y * y);

        pc_49[k] = e_0 * (std::sqrt(46.142578125) * x * x * x * x * x * x * x * y * z + std::sqrt(128.173828125) * x * x * x * x * x * y * y * y * z - std::sqrt(738.28125) * x * x * x * x * x * y * z * z * z + std::sqrt(5.126953125) * x * x * x * y * y * y * y * y * z - std::sqrt(328.125) * x * x * x * y * y * y * z * z * z + std::sqrt(118.125) * x * x * x * y * z * z * z * z * z - std::sqrt(5.126953125) * x * y * y * y * y * y * y * y * z + std::sqrt(82.03125) * x * y * y * y * y * y * z * z * z - std::sqrt(13.125) * x * y * y * y * z * z * z * z * z) + e_1 * (std::sqrt(9043.9453125) * x * x * x * x * x * y * z + std::sqrt(6644.53125) * x * x * x * y * y * y * z - std::sqrt(47250.0) * x * x * x * y * z * z * z - std::sqrt(184.5703125) * x * y * y * y * y * y * z + std::sqrt(472.5) * x * y * z * z * z * z * z) + e_2 * (std::sqrt(73828.125) * x * x * x * y * z + std::sqrt(2953.125) * x * y * y * y * z - std::sqrt(47250.0) * x * y * z * z * z) + e_3 * (std::sqrt(11812.5) * x * y * z);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        pc_50[k] = e_0 * (std::sqrt(123.046875) * x * x * x * x * x * x * y * z * z + std::sqrt(492.1875) * x * x * x * x * y * y * y * z * z - std::sqrt(1968.75) * x * x * x * x * y * z * z * z * z + std::sqrt(123.046875) * x * x * y * y * y * y * y * z * z - std::sqrt(1968.75) * x * x * y * y * y * z * z * z * z + std::sqrt(315.0) * x * x * y * z * z * z * z * z * z) + e_1 * (std::sqrt(123.046875) * x * x * x * x * x * x * y + std::sqrt(492.1875) * x * x * x * x * y * y * y - std::sqrt(1107.421875) * x * x * x * x * y * z * z + std::sqrt(123.046875) * x * x * y * y * y * y * y - std::sqrt(492.1875) * x * x * y * y * y * z * z - std::sqrt(17718.75) * x * x * y * z * z * z * z + std::sqrt(123.046875) * y * y * y * y * y * z * z - std::sqrt(1968.75) * y * y * y * z * z * z * z + std::sqrt(315.0) * y * z * z * z * z * z * z) + e_2 * (std::sqrt(9966.796875) * x * x * x * x * y + std::sqrt(12304.6875) * x * x * y * y * y - std::sqrt(283500.0) * x * x * y * z * z + std::sqrt(123.046875) * y * y * y * y * y - std::sqrt(7875.0) * y * y * y * z * z) + e_3 * (std::sqrt(17718.75) * x * x * y + std::sqrt(1968.75) * y * y * y - std::sqrt(70875.0) * y * z * z);

        pc_51[k] = e_0 * (-std::sqrt(3.076171875) * x * x * x * x * x * x * x * y * z - std::sqrt(27.685546875) * x * x * x * x * x * y * y * y * z + std::sqrt(196.875) * x * x * x * x * x * y * z * z * z - std::sqrt(27.685546875) * x * x * x * y * y * y * y * y * z + std::sqrt(787.5) * x * x * x * y * y * y * z * z * z - std::sqrt(952.875) * x * x * x * y * z * z * z * z * z - std::sqrt(3.076171875) * x * y * y * y * y * y * y * y * z + std::sqrt(196.875) * x * y * y * y * y * y * z * z * z - std::sqrt(952.875) * x * y * y * y * z * z * z * z * z + std::sqrt(126.0) * x * y * z * z * z * z * z * z * z) + e_1 * (-std::sqrt(110.7421875) * x * x * x * x * x * y * z - std::sqrt(442.96875) * x * x * x * y * y * y * z - std::sqrt(7087.5) * x * x * x * y * z * z * z - std::sqrt(110.7421875) * x * y * y * y * y * y * z - std::sqrt(7087.5) * x * y * y * y * z * z * z + std::sqrt(2551.5) * x * y * z * z * z * z * z) + e_2 * (-std::sqrt(44296.875) * x * x * x * y * z - std::sqrt(44296.875) * x * y * y * y * z) + e_3 * (-std::sqrt(177187.5) * x * y * z);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        pc_52[k] = e_0 * (-std::sqrt(18.45703125) * x * x * x * x * x * x * x * z * z - std::sqrt(166.11328125) * x * x * x * x * x * y * y * z * z + std::sqrt(401.953125) * x * x * x * x * x * z * z * z * z - std::sqrt(166.11328125) * x * x * x * y * y * y * y * z * z + std::sqrt(1607.8125) * x * x * x * y * y * z * z * z * z - std::sqrt(336.0) * x * x * x * z * z * z * z * z * z - std::sqrt(18.45703125) * x * y * y * y * y * y * y * z * z + std::sqrt(401.953125) * x * y * y * y * y * z * z * z * z - std::sqrt(336.0) * x * y * y * z * z * z * z * z * z + std::sqrt(21.0) * x * z * z * z * z * z * z * z * z) + e_1 * (-std::sqrt(18.45703125) * x * x * x * x * x * x * x - std::sqrt(166.11328125) * x * x * x * x * x * y * y + std::sqrt(295.3125) * x * x * x * x * x * z * z - std::sqrt(166.11328125) * x * x * x * y * y * y * y + std::sqrt(1181.25) * x * x * x * y * y * z * z - std::sqrt(1181.25) * x * x * x * z * z * z * z - std::sqrt(18.45703125) * x * y * y * y * y * y * y + std::sqrt(295.3125) * x * y * y * y * y * z * z - std::sqrt(1181.25) * x * y * y * z * z * z * z + std::sqrt(3024.0) * x * z * z * z * z * z * z) + e_2 * (-std::sqrt(1845.703125) * x * x * x * x * x - std::sqrt(7382.8125) * x * x * x * y * y - std::sqrt(1845.703125) * x * y * y * y * y + std::sqrt(118125.0) * x * z * z * z * z) + e_3 * (-std::sqrt(29531.25) * x * x * x - std::sqrt(29531.25) * x * y * y + std::sqrt(472500.0) * x * z * z);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        pc_53[k] = e_0 * (-std::sqrt(3.076171875) * x * x * x * x * x * x * x * x * z - std::sqrt(27.685546875) * x * x * x * x * x * x * y * y * z + std::sqrt(196.875) * x * x * x * x * x * x * z * z * z - std::sqrt(27.685546875) * x * x * x * x * y * y * y * y * z + std::sqrt(787.5) * x * x * x * x * y * y * z * z * z - std::sqrt(952.875) * x * x * x * x * z * z * z * z * z - std::sqrt(3.076171875) * x * x * y * y * y * y * y * y * z + std::sqrt(196.875) * x * x * y * y * y * y * z * z * z - std::sqrt(952.875) * x * x * y * y * z * z * z * z * z + std::sqrt(126.0) * x * x * z * z * z * z * z * z * z) + e_1 * (-std::sqrt(150.732421875) * x * x * x * x * x * x * z - std::sqrt(692.138671875) * x * x * x * x * y * y * z - std::sqrt(4921.875) * x * x * x * x * z * z * z - std::sqrt(249.169921875) * x * x * y * y * y * y * z - std::sqrt(3150.0) * x * x * y * y * z * z * z + std::sqrt(385.875) * x * x * z * z * z * z * z - std::sqrt(3.076171875) * y * y * y * y * y * y * z + std::sqrt(196.875) * y * y * y * y * z * z * z - std::sqrt(952.875) * y * y * z * z * z * z * z + std::sqrt(126.0) * z * z * z * z * z * z * z) + e_2 * (-std::sqrt(44296.875) * x * x * x * x * z - std::sqrt(44296.875) * x * x * y * y * z - std::sqrt(19687.5) * x * x * z * z * z - std::sqrt(19687.5) * y * y * z * z * z + std::sqrt(12600.0) * z * z * z * z * z) + e_3 * (-std::sqrt(398671.875) * x * x * z - std::sqrt(44296.875) * y * y * z + std::sqrt(78750.0) * z * z * z);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        pc_54[k] = e_0 * (std::sqrt(30.76171875) * x * x * x * x * x * x * x * z * z + std::sqrt(30.76171875) * x * x * x * x * x * y * y * z * z - std::sqrt(492.1875) * x * x * x * x * x * z * z * z * z - std::sqrt(30.76171875) * x * x * x * y * y * y * y * z * z + std::sqrt(78.75) * x * x * x * z * z * z * z * z * z - std::sqrt(30.76171875) * x * y * y * y * y * y * y * z * z + std::sqrt(492.1875) * x * y * y * y * y * z * z * z * z - std::sqrt(78.75) * x * y * y * z * z * z * z * z * z) + e_1 * (std::sqrt(30.76171875) * x * x * x * x * x * x * x + std::sqrt(30.76171875) * x * x * x * x * x * y * y - std::sqrt(123.046875) * x * x * x * x * x * z * z - std::sqrt(30.76171875) * x * x * x * y * y * y * y + std::sqrt(492.1875) * x * x * x * y * y * z * z - std::sqrt(7875.0) * x * x * x * z * z * z * z - std::sqrt(30.76171875) * x * y * y * y * y * y * y + std::sqrt(1107.421875) * x * y * y * y * y * z * z + std::sqrt(315.0) * x * z * z * z * z * z * z) + e_2 * (std::sqrt(3076.171875) * x * x * x * x * x + std::sqrt(492.1875) * x * x * x * y * y - std::sqrt(96468.75) * x * x * x * z * z - std::sqrt(1107.421875) * x * y * y * y * y + std::sqrt(17718.75) * x * y * y * z * z) + e_3 * (std::sqrt(7875.0) * x * x * x - std::sqrt(70875.0) * x * z * z);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        pc_55[k] = e_0 * (std::sqrt(5.126953125) * x * x * x * x * x * x * x * x * z - std::sqrt(5.126953125) * x * x * x * x * x * x * y * y * z - std::sqrt(82.03125) * x * x * x * x * x * x * z * z * z - std::sqrt(128.173828125) * x * x * x * x * y * y * y * y * z + std::sqrt(328.125) * x * x * x * x * y * y * z * z * z + std::sqrt(13.125) * x * x * x * x * z * z * z * z * z - std::sqrt(46.142578125) * x * x * y * y * y * y * y * y * z + std::sqrt(738.28125) * x * x * y * y * y * y * z * z * z - std::sqrt(118.125) * x * x * y * y * z * z * z * z * z) + e_1 * (std::sqrt(1153.564453125) * x * x * x * x * x * x * z - std::sqrt(2260.986328125) * x * x * x * x * y * y * z - std::sqrt(6644.53125) * x * x * x * x * z * z * z - std::sqrt(7798.095703125) * x * x * y * y * y * y * z + std::sqrt(26578.125) * x * x * y * y * z * z * z + std::sqrt(118.125) * x * x * z * z * z * z * z - std::sqrt(46.142578125) * y * y * y * y * y * y * z + std::sqrt(738.28125) * y * y * y * y * z * z * z - std::sqrt(118.125) * y * y * z * z * z * z * z) + e_2 * (std::sqrt(11812.5) * x * x * x * x * z - std::sqrt(26578.125) * x * x * y * y * z - std::sqrt(11812.5) * x * x * z * z * z - std::sqrt(2953.125) * y * y * y * y * z + std::sqrt(11812.5) * y * y * z * z * z) + e_3 * (std::sqrt(2953.125) * x * x * z - std::sqrt(2953.125) * y * y * z);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        pc_56[k] = e_0 * (std::sqrt(1.153564453125) * x * x * x * x * x * x * x * x * y + std::sqrt(0.5126953125) * x * x * x * x * x * x * y * y * y - std::sqrt(295.3125) * x * x * x * x * x * x * y * z * z - std::sqrt(2.05078125) * x * x * x * x * y * y * y * y * y + std::sqrt(32.8125) * x * x * x * x * y * y * y * z * z + std::sqrt(295.3125) * x * x * x * x * y * z * z * z * z - std::sqrt(0.5126953125) * x * x * y * y * y * y * y * y * y + std::sqrt(295.3125) * x * x * y * y * y * y * y * z * z - std::sqrt(525.0) * x * x * y * y * y * z * z * z * z + std::sqrt(0.128173828125) * y * y * y * y * y * y * y * y * y - std::sqrt(32.8125) * y * y * y * y * y * y * y * z * z + std::sqrt(32.8125) * y * y * y * y * y * z * z * z * z) + e_1 * (std::sqrt(226.0986328125) * x * x * x * x * x * x * y + std::sqrt(4.6142578125) * x * x * x * x * y * y * y - std::sqrt(18900.0) * x * x * x * x * y * z * z - std::sqrt(41.5283203125) * x * x * y * y * y * y * y + std::sqrt(4725.0) * x * x * y * y * y * z * z + std::sqrt(1181.25) * x * x * y * z * z * z * z + std::sqrt(41.5283203125) * y * y * y * y * y * y * y - std::sqrt(4725.0) * y * y * y * y * y * z * z + std::sqrt(1181.25) * y * y * y * z * z * z * z) + e_2 * (std::sqrt(2233.30078125) * x * x * x * x * y + std::sqrt(73.828125) * x * x * y * y * y - std::sqrt(42525.0) * x * x * y * z * z + std::sqrt(904.39453125) * y * y * y * y * y - std::sqrt(42525.0) * y * y * y * z * z + std::sqrt(4725.0) * y * z * z * z * z) + e_3 * (std::sqrt(1181.25) * x * x * y + std::sqrt(1181.25) * y * y * y - std::sqrt(18900.0) * y * z * z);

        pc_57[k] = e_0 * (std::sqrt(3.076171875) * x * x * x * x * x * x * x * y * z + std::sqrt(3.076171875) * x * x * x * x * x * y * y * y * z - std::sqrt(787.5) * x * x * x * x * x * y * z * z * z - std::sqrt(3.076171875) * x * x * x * y * y * y * y * y * z + std::sqrt(787.5) * x * x * x * y * z * z * z * z * z - std::sqrt(3.076171875) * x * y * y * y * y * y * y * y * z + std::sqrt(787.5) * x * y * y * y * y * y * z * z * z - std::sqrt(787.5) * x * y * y * y * z * z * z * z * z) + e_1 * (-std::sqrt(1771.875) * x * x * x * x * x * y * z + std::sqrt(1771.875) * x * y * y * y * y * y * z) + e_2 * (-std::sqrt(44296.875) * x * x * x * y * z + std::sqrt(44296.875) * x * y * y * y * z);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        pc_58[k] = e_0 * (-std::sqrt(0.076904296875) * x * x * x * x * x * x * x * x * y - std::sqrt(0.3076171875) * x * x * x * x * x * x * y * y * y + std::sqrt(30.76171875) * x * x * x * x * x * x * y * z * z + std::sqrt(30.76171875) * x * x * x * x * y * y * y * z * z - std::sqrt(492.1875) * x * x * x * x * y * z * z * z * z + std::sqrt(0.3076171875) * x * x * y * y * y * y * y * y * y - std::sqrt(30.76171875) * x * x * y * y * y * y * y * z * z + std::sqrt(315.0) * x * x * y * z * z * z * z * z * z + std::sqrt(0.076904296875) * y * y * y * y * y * y * y * y * y - std::sqrt(30.76171875) * y * y * y * y * y * y * y * z * z + std::sqrt(492.1875) * y * y * y * y * y * z * z * z * z - std::sqrt(315.0) * y * y * y * z * z * z * z * z * z) + e_1 * (-std::sqrt(15.0732421875) * x * x * x * x * x * x * y - std::sqrt(7.6904296875) * x * x * x * x * y * y * y - std::sqrt(1107.421875) * x * x * x * x * y * z * z + std::sqrt(37.2216796875) * x * x * y * y * y * y * y - std::sqrt(492.1875) * x * x * y * y * y * z * z + std::sqrt(17718.75) * x * x * y * z * z * z * z + std::sqrt(24.9169921875) * y * y * y * y * y * y * y + std::sqrt(123.046875) * y * y * y * y * y * z * z - std::sqrt(1968.75) * y * y * y * z * z * z * z - std::sqrt(1260.0) * y * z * z * z * z * z * z) + e_2 * (-std::sqrt(2491.69921875) * x * x * x * x * y + std::sqrt(123.046875) * x * x * y * y * y + std::sqrt(70875.0) * x * x * y * z * z + std::sqrt(3722.16796875) * y * y * y * y * y - std::sqrt(7875.0) * y * y * y * z * z - std::sqrt(70875.0) * y * z * z * z * z) + e_3 * (std::sqrt(31500.0) * y * y * y - std::sqrt(283500.0) * y * z * z);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        pc_59[k] = e_0 * (-std::sqrt(0.46142578125) * x * x * x * x * x * x * x * x * z - std::sqrt(1.845703125) * x * x * x * x * x * x * y * y * z + std::sqrt(128.173828125) * x * x * x * x * x * x * z * z * z + std::sqrt(128.173828125) * x * x * x * x * y * y * z * z * z - std::sqrt(328.125) * x * x * x * x * z * z * z * z * z + std::sqrt(1.845703125) * x * x * y * y * y * y * y * y * z - std::sqrt(128.173828125) * x * x * y * y * y * y * z * z * z + std::sqrt(52.5) * x * x * z * z * z * z * z * z * z + std::sqrt(0.46142578125) * y * y * y * y * y * y * y * y * z - std::sqrt(128.173828125) * y * y * y * y * y * y * z * z * z + std::sqrt(328.125) * y * y * y * y * z * z * z * z * z - std::sqrt(52.5) * y * y * z * z * z * z * z * z * z) + e_1 * (std::sqrt(184.5703125) * x * x * x * x * x * x * z + std::sqrt(184.5703125) * x * x * x * x * y * y * z - std::sqrt(184.5703125) * x * x * y * y * y * y * z + std::sqrt(1890.0) * x * x * z * z * z * z * z - std::sqrt(184.5703125) * y * y * y * y * y * y * z - std::sqrt(1890.0) * y * y * z * z * z * z * z) + e_2 * (std::sqrt(11812.5) * x * x * x * x * z + std::sqrt(47250.0) * x * x * z * z * z - std::sqrt(11812.5) * y * y * y * y * z - std::sqrt(47250.0) * y * y * z * z * z) + e_3 * (std::sqrt(189000.0) * x * x * z - std::sqrt(189000.0) * y * y * z);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        pc_60[k] = e_0 * (-std::sqrt(0.076904296875) * x * x * x * x * x * x * x * x * x - std::sqrt(0.3076171875) * x * x * x * x * x * x * x * y * y + std::sqrt(30.76171875) * x * x * x * x * x * x * x * z * z + std::sqrt(30.76171875) * x * x * x * x * x * y * y * z * z - std::sqrt(492.1875) * x * x * x * x * x * z * z * z * z + std::sqrt(0.3076171875) * x * x * x * y * y * y * y * y * y - std::sqrt(30.76171875) * x * x * x * y * y * y * y * z * z + std::sqrt(315.0) * x * x * x * z * z * z * z * z * z + std::sqrt(0.076904296875) * x * y * y * y * y * y * y * y * y - std::sqrt(30.76171875) * x * y * y * y * y * y * y * z * z + std::sqrt(492.1875) * x * y * y * y * y * z * z * z * z - std::sqrt(315.0) * x * y * y * z * z * z * z * z * z) + e_1 * (-std::sqrt(24.9169921875) * x * x * x * x * x * x * x - std::sqrt(37.2216796875) * x * x * x * x * x * y * y - std::sqrt(123.046875) * x * x * x * x * x * z * z + std::sqrt(7.6904296875) * x * x * x * y * y * y * y + std::sqrt(492.1875) * x * x * x * y * y * z * z + std::sqrt(1968.75) * x * x * x * z * z * z * z + std::sqrt(15.0732421875) * x * y * y * y * y * y * y + std::sqrt(1107.421875) * x * y * y * y * y * z * z - std::sqrt(17718.75) * x * y * y * z * z * z * z + std::sqrt(1260.0) * x * z * z * z * z * z * z) + e_2 * (-std::sqrt(3722.16796875) * x * x * x * x * x - std::sqrt(123.046875) * x * x * x * y * y + std::sqrt(7875.0) * x * x * x * z * z + std::sqrt(2491.69921875) * x * y * y * y * y - std::sqrt(70875.0) * x * y * y * z * z + std::sqrt(70875.0) * x * z * z * z * z) + e_3 * (-std::sqrt(31500.0) * x * x * x + std::sqrt(283500.0) * x * z * z);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        pc_61[k] = e_0 * (std::sqrt(0.76904296875) * x * x * x * x * x * x * x * x * z - std::sqrt(196.875) * x * x * x * x * x * x * z * z * z - std::sqrt(3.076171875) * x * x * x * x * y * y * y * y * z + std::sqrt(196.875) * x * x * x * x * y * y * z * z * z + std::sqrt(196.875) * x * x * x * x * z * z * z * z * z + std::sqrt(196.875) * x * x * y * y * y * y * z * z * z - std::sqrt(787.5) * x * x * y * y * z * z * z * z * z + std::sqrt(0.76904296875) * y * y * y * y * y * y * y * y * z - std::sqrt(196.875) * y * y * y * y * y * y * z * z * z + std::sqrt(196.875) * y * y * y * y * z * z * z * z * z) + e_1 * (-std::sqrt(307.6171875) * x * x * x * x * x * x * z + std::sqrt(996.6796875) * x * x * x * x * y * y * z - std::sqrt(3150.0) * x * x * x * x * z * z * z + std::sqrt(996.6796875) * x * x * y * y * y * y * z - std::sqrt(12600.0) * x * x * y * y * z * z * z + std::sqrt(3150.0) * x * x * z * z * z * z * z - std::sqrt(307.6171875) * y * y * y * y * y * y * z - std::sqrt(3150.0) * y * y * y * y * z * z * z + std::sqrt(3150.0) * y * y * z * z * z * z * z) + e_2 * (-std::sqrt(39977.9296875) * x * x * x * x * z + std::sqrt(442.96875) * x * x * y * y * z + std::sqrt(3150.0) * x * x * z * z * z - std::sqrt(39977.9296875) * y * y * y * y * z + std::sqrt(3150.0) * y * y * z * z * z + std::sqrt(3150.0) * z * z * z * z * z) + e_3 * (-std::sqrt(113400.0) * x * x * z - std::sqrt(113400.0) * y * y * z + std::sqrt(50400.0) * z * z * z);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        pc_62[k] = e_0 * (std::sqrt(0.128173828125) * x * x * x * x * x * x * x * x * x - std::sqrt(0.5126953125) * x * x * x * x * x * x * x * y * y - std::sqrt(32.8125) * x * x * x * x * x * x * x * z * z - std::sqrt(2.05078125) * x * x * x * x * x * y * y * y * y + std::sqrt(295.3125) * x * x * x * x * x * y * y * z * z + std::sqrt(32.8125) * x * x * x * x * x * z * z * z * z + std::sqrt(0.5126953125) * x * x * x * y * y * y * y * y * y + std::sqrt(32.8125) * x * x * x * y * y * y * y * z * z - std::sqrt(525.0) * x * x * x * y * y * z * z * z * z + std::sqrt(1.153564453125) * x * y * y * y * y * y * y * y * y - std::sqrt(295.3125) * x * y * y * y * y * y * y * z * z + std::sqrt(295.3125) * x * y * y * y * y * z * z * z * z) + e_1 * (std::sqrt(41.5283203125) * x * x * x * x * x * x * x - std::sqrt(41.5283203125) * x * x * x * x * x * y * y - std::sqrt(4725.0) * x * x * x * x * x * z * z + std::sqrt(4.6142578125) * x * x * x * y * y * y * y + std::sqrt(4725.0) * x * x * x * y * y * z * z + std::sqrt(1181.25) * x * x * x * z * z * z * z + std::sqrt(226.0986328125) * x * y * y * y * y * y * y - std::sqrt(18900.0) * x * y * y * y * y * z * z + std::sqrt(1181.25) * x * y * y * z * z * z * z) + e_2 * (std::sqrt(904.39453125) * x * x * x * x * x + std::sqrt(73.828125) * x * x * x * y * y - std::sqrt(42525.0) * x * x * x * z * z + std::sqrt(2233.30078125) * x * y * y * y * y - std::sqrt(42525.0) * x * y * y * z * z + std::sqrt(4725.0) * x * z * z * z * z) + e_3 * (std::sqrt(1181.25) * x * x * x + std::sqrt(1181.25) * x * y * y - std::sqrt(18900.0) * x * z * z);

        pc_63[k] = e_0 * (-std::sqrt(41.5283203125) * x * x * x * x * x * x * x * y * z + std::sqrt(226.0986328125) * x * x * x * x * x * y * y * y * z + std::sqrt(295.3125) * x * x * x * x * x * y * z * z * z + std::sqrt(226.0986328125) * x * x * x * y * y * y * y * y * z - std::sqrt(3281.25) * x * x * x * y * y * y * z * z * z - std::sqrt(41.5283203125) * x * y * y * y * y * y * y * y * z + std::sqrt(295.3125) * x * y * y * y * y * y * z * z * z) + e_1 * (-std::sqrt(1495.01953125) * x * x * x * x * x * y * z + std::sqrt(16611.328125) * x * x * x * y * y * y * z - std::sqrt(1495.01953125) * x * y * y * y * y * y * z);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        pc_64[k] = e_0 * (-std::sqrt(110.7421875) * x * x * x * x * x * x * y * z * z + std::sqrt(442.96875) * x * x * x * x * y * y * y * z * z + std::sqrt(787.5) * x * x * x * x * y * z * z * z * z + std::sqrt(996.6796875) * x * x * y * y * y * y * y * z * z - std::sqrt(7087.5) * x * x * y * y * y * z * z * z * z) + e_1 * (-std::sqrt(110.7421875) * x * x * x * x * x * x * y + std::sqrt(442.96875) * x * x * x * x * y * y * y + std::sqrt(5426.3671875) * x * x * x * x * y * z * z + std::sqrt(996.6796875) * x * x * y * y * y * y * y - std::sqrt(3986.71875) * x * x * y * y * y * z * z - std::sqrt(7087.5) * x * x * y * z * z * z * z + std::sqrt(996.6796875) * y * y * y * y * y * z * z - std::sqrt(7087.5) * y * y * y * z * z * z * z) + e_2 * (-std::sqrt(110.7421875) * x * x * x * x * y + std::sqrt(35880.46875) * x * x * y * y * y - std::sqrt(15946.875) * x * x * y * z * z + std::sqrt(996.6796875) * y * y * y * y * y - std::sqrt(15946.875) * y * y * y * z * z - std::sqrt(28350.0) * y * z * z * z * z) + e_3 * (std::sqrt(15946.875) * x * x * y + std::sqrt(15946.875) * y * y * y - std::sqrt(255150.0) * y * z * z);

        pc_65[k] = e_0 * (std::sqrt(2.7685546875) * x * x * x * x * x * x * x * y * z - std::sqrt(2.7685546875) * x * x * x * x * x * y * y * y * z - std::sqrt(123.046875) * x * x * x * x * x * y * z * z * z - std::sqrt(69.2138671875) * x * x * x * y * y * y * y * y * z + std::sqrt(492.1875) * x * x * x * y * y * y * z * z * z + std::sqrt(315.0) * x * x * x * y * z * z * z * z * z - std::sqrt(24.9169921875) * x * y * y * y * y * y * y * y * z + std::sqrt(1107.421875) * x * y * y * y * y * y * z * z * z - std::sqrt(2835.0) * x * y * y * y * z * z * z * z * z) + e_1 * (-std::sqrt(11.07421875) * x * x * x * x * x * y * z - std::sqrt(1107.421875) * x * x * x * y * y * y * z + std::sqrt(17718.75) * x * x * x * y * z * z * z - std::sqrt(897.01171875) * x * y * y * y * y * y * z - std::sqrt(17718.75) * x * y * y * y * z * z * z - std::sqrt(11340.0) * x * y * z * z * z * z * z) + e_2 * (std::sqrt(17718.75) * x * x * x * y * z - std::sqrt(159468.75) * x * y * y * y * z - std::sqrt(283500.0) * x * y * z * z * z) + e_3 * (-std::sqrt(637875.0) * x * y * z);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        pc_66[k] = e_0 * (std::sqrt(16.611328125) * x * x * x * x * x * x * x * z * z - std::sqrt(16.611328125) * x * x * x * x * x * y * y * z * z - std::sqrt(184.5703125) * x * x * x * x * x * z * z * z * z - std::sqrt(415.283203125) * x * x * x * y * y * y * y * z * z + std::sqrt(738.28125) * x * x * x * y * y * z * z * z * z + std::sqrt(52.5) * x * x * x * z * z * z * z * z * z - std::sqrt(149.501953125) * x * y * y * y * y * y * y * z * z + std::sqrt(1661.1328125) * x * y * y * y * y * z * z * z * z - std::sqrt(472.5) * x * y * y * z * z * z * z * z * z) + e_1 * (std::sqrt(16.611328125) * x * x * x * x * x * x * x - std::sqrt(16.611328125) * x * x * x * x * x * y * y - std::sqrt(415.283203125) * x * x * x * y * y * y * y - std::sqrt(149.501953125) * x * y * y * y * y * y * y) + e_2 * (std::sqrt(1661.1328125) * x * x * x * x * x - std::sqrt(6644.53125) * x * x * x * y * y - std::sqrt(14950.1953125) * x * y * y * y * y) + e_3 * (std::sqrt(11812.5) * x * x * x - std::sqrt(106312.5) * x * y * y);

        pc_67[k] = e_0 * (std::sqrt(2.7685546875) * x * x * x * x * x * x * x * x * z - std::sqrt(2.7685546875) * x * x * x * x * x * x * y * y * z - std::sqrt(123.046875) * x * x * x * x * x * x * z * z * z - std::sqrt(69.2138671875) * x * x * x * x * y * y * y * y * z + std::sqrt(492.1875) * x * x * x * x * y * y * z * z * z + std::sqrt(315.0) * x * x * x * x * z * z * z * z * z - std::sqrt(24.9169921875) * x * x * y * y * y * y * y * y * z + std::sqrt(1107.421875) * x * x * y * y * y * y * z * z * z - std::sqrt(2835.0) * x * x * y * y * z * z * z * z * z) + e_1 * (std::sqrt(135.6591796875) * x * x * x * x * x * x * z - std::sqrt(69.2138671875) * x * x * x * x * y * y * z + std::sqrt(1107.421875) * x * x * x * x * z * z * z - std::sqrt(622.9248046875) * x * x * y * y * y * y * z - std::sqrt(39867.1875) * x * x * y * y * z * z * z + std::sqrt(2835.0) * x * x * z * z * z * z * z - std::sqrt(24.9169921875) * y * y * y * y * y * y * z + std::sqrt(1107.421875) * y * y * y * y * z * z * z - std::sqrt(2835.0) * y * y * z * z * z * z * z) + e_2 * (std::sqrt(17718.75) * x * x * x * x * z - std::sqrt(159468.75) * x * x * y * y * z + std::sqrt(70875.0) * x * x * z * z * z - std::sqrt(70875.0) * y * y * z * z * z) + e_3 * (std::sqrt(159468.75) * x * x * z - std::sqrt(159468.75) * y * y * z);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        pc_68[k] = e_0 * (-std::sqrt(27.685546875) * x * x * x * x * x * x * x * z * z + std::sqrt(249.169921875) * x * x * x * x * x * y * y * z * z + std::sqrt(196.875) * x * x * x * x * x * z * z * z * z + std::sqrt(27.685546875) * x * x * x * y * y * y * y * z * z - std::sqrt(3150.0) * x * x * x * y * y * z * z * z * z - std::sqrt(249.169921875) * x * y * y * y * y * y * y * z * z + std::sqrt(1771.875) * x * y * y * y * y * z * z * z * z) + e_1 * (-std::sqrt(27.685546875) * x * x * x * x * x * x * x + std::sqrt(249.169921875) * x * x * x * x * x * y * y - std::sqrt(110.7421875) * x * x * x * x * x * z * z + std::sqrt(27.685546875) * x * x * x * y * y * y * y - std::sqrt(21705.46875) * x * x * x * y * y * z * z + std::sqrt(7087.5) * x * x * x * z * z * z * z - std::sqrt(249.169921875) * x * y * y * y * y * y * y + std::sqrt(996.6796875) * x * y * y * y * y * z * z + std::sqrt(7087.5) * x * y * y * z * z * z * z) + e_2 * (-std::sqrt(2768.5546875) * x * x * x * x * x + std::sqrt(442.96875) * x * x * x * y * y + std::sqrt(15946.875) * x * x * x * z * z - std::sqrt(8970.1171875) * x * y * y * y * y + std::sqrt(15946.875) * x * y * y * z * z + std::sqrt(28350.0) * x * z * z * z * z) + e_3 * (-std::sqrt(15946.875) * x * x * x - std::sqrt(15946.875) * x * y * y + std::sqrt(255150.0) * x * z * z);

        pc_69[k] = e_0 * (-std::sqrt(4.6142578125) * x * x * x * x * x * x * x * x * z + std::sqrt(115.3564453125) * x * x * x * x * x * x * y * y * z + std::sqrt(32.8125) * x * x * x * x * x * x * z * z * z - std::sqrt(41.5283203125) * x * x * x * x * y * y * y * y * z - std::sqrt(1181.25) * x * x * x * x * y * y * z * z * z - std::sqrt(373.7548828125) * x * x * y * y * y * y * y * y * z + std::sqrt(2657.8125) * x * x * y * y * y * y * z * z * z) + e_1 * (-std::sqrt(1038.2080078125) * x * x * x * x * x * x * z + std::sqrt(373.7548828125) * x * x * x * x * y * y * z + std::sqrt(2657.8125) * x * x * x * x * z * z * z - std::sqrt(30274.1455078125) * x * x * y * y * y * y * z + std::sqrt(10631.25) * x * x * y * y * z * z * z - std::sqrt(373.7548828125) * y * y * y * y * y * y * z + std::sqrt(2657.8125) * y * y * y * y * z * z * z) + e_2 * (-std::sqrt(23920.3125) * x * x * x * x * z - std::sqrt(95681.25) * x * x * y * y * z + std::sqrt(42525.0) * x * x * z * z * z - std::sqrt(23920.3125) * y * y * y * y * z + std::sqrt(42525.0) * y * y * z * z * z) + e_3 * (-std::sqrt(42525.0) * x * x * z - std::sqrt(42525.0) * y * y * z + std::sqrt(18900.0) * z * z * z);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        pc_70[k] = e_0 * (-std::sqrt(1.38427734375) * x * x * x * x * x * x * x * x * y + std::sqrt(39.375) * x * x * x * x * x * x * y * y * y + std::sqrt(138.427734375) * x * x * x * x * x * x * y * z * z + std::sqrt(15.380859375) * x * x * x * x * y * y * y * y * y - std::sqrt(5552.490234375) * x * x * x * x * y * y * y * z * z - std::sqrt(9.84375) * x * x * y * y * y * y * y * y * y + std::sqrt(1245.849609375) * x * x * y * y * y * y * y * z * z + std::sqrt(0.15380859375) * y * y * y * y * y * y * y * y * y - std::sqrt(15.380859375) * y * y * y * y * y * y * y * z * z) + e_1 * (-std::sqrt(5.537109375) * x * x * x * x * x * x * y + std::sqrt(3460.693359375) * x * x * x * x * y * y * y - std::sqrt(2214.84375) * x * x * x * x * y * z * z - std::sqrt(49.833984375) * x * x * y * y * y * y * y - std::sqrt(8859.375) * x * x * y * y * y * z * z + std::sqrt(49.833984375) * y * y * y * y * y * y * y - std::sqrt(2214.84375) * y * y * y * y * y * z * z) + e_2 * (std::sqrt(2214.84375) * x * x * x * x * y + std::sqrt(8859.375) * x * x * y * y * y - std::sqrt(79734.375) * x * x * y * z * z + std::sqrt(2214.84375) * y * y * y * y * y - std::sqrt(79734.375) * y * y * y * z * z) + e_3 * (std::sqrt(8859.375) * x * x * y + std::sqrt(8859.375) * y * y * y - std::sqrt(141750.0) * y * z * z);

        pc_71[k] = e_0 * (-std::sqrt(3.69140625) * x * x * x * x * x * x * x * y * z + std::sqrt(92.28515625) * x * x * x * x * x * y * y * y * z + std::sqrt(369.140625) * x * x * x * x * x * y * z * z * z + std::sqrt(92.28515625) * x * x * x * y * y * y * y * y * z - std::sqrt(13289.0625) * x * x * x * y * y * y * z * z * z - std::sqrt(3.69140625) * x * y * y * y * y * y * y * y * z + std::sqrt(369.140625) * x * y * y * y * y * y * z * z * z) + e_1 * (std::sqrt(2126.25) * x * x * x * x * x * y * z - std::sqrt(23625.0) * x * x * x * y * y * y * z - std::sqrt(23625.0) * x * x * x * y * z * z * z + std::sqrt(2126.25) * x * y * y * y * y * y * z - std::sqrt(23625.0) * x * y * y * y * z * z * z) + e_2 * (-std::sqrt(53156.25) * x * x * x * y * z - std::sqrt(53156.25) * x * y * y * y * z - std::sqrt(212625.0) * x * y * z * z * z) + e_3 * (-std::sqrt(850500.0) * x * y * z);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        pc_72[k] = e_0 * (std::sqrt(0.09228515625) * x * x * x * x * x * x * x * x * y - std::sqrt(1.4765625) * x * x * x * x * x * x * y * y * y - std::sqrt(18.087890625) * x * x * x * x * x * x * y * z * z - std::sqrt(9.228515625) * x * x * x * x * y * y * y * y * y + std::sqrt(452.197265625) * x * x * x * x * y * y * y * z * z + std::sqrt(147.65625) * x * x * x * x * y * z * z * z * z - std::sqrt(1.4765625) * x * x * y * y * y * y * y * y * y + std::sqrt(452.197265625) * x * x * y * y * y * y * y * z * z - std::sqrt(5315.625) * x * x * y * y * y * z * z * z * z + std::sqrt(0.09228515625) * y * y * y * y * y * y * y * y * y - std::sqrt(18.087890625) * y * y * y * y * y * y * y * z * z + std::sqrt(147.65625) * y * y * y * y * y * z * z * z * z) + e_1 * (std::sqrt(0.369140625) * x * x * x * x * x * x * y - std::sqrt(747.509765625) * x * x * x * x * y * y * y + std::sqrt(5315.625) * x * x * x * x * y * z * z - std::sqrt(505.353515625) * x * x * y * y * y * y * y - std::sqrt(9450.0) * x * x * y * y * y * z * z - std::sqrt(21262.5) * x * x * y * z * z * z * z + std::sqrt(29.900390625) * y * y * y * y * y * y * y + std::sqrt(23.625) * y * y * y * y * y * z * z + std::sqrt(2362.5) * y * y * y * z * z * z * z) + e_2 * (-std::sqrt(59062.5) * x * x * y * y * y - std::sqrt(132890.625) * x * x * y * z * z + std::sqrt(2362.5) * y * y * y * y * y + std::sqrt(14765.625) * y * y * y * z * z) + e_3 * (-std::sqrt(132890.625) * x * x * y + std::sqrt(14765.625) * y * y * y);
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

        pc_73[k] = e_0 * (std::sqrt(0.5537109375) * x * x * x * x * x * x * x * x * z - std::sqrt(8.859375) * x * x * x * x * x * x * y * y * z - std::sqrt(63.0) * x * x * x * x * x * x * z * z * z - std::sqrt(55.37109375) * x * x * x * x * y * y * y * y * z + std::sqrt(1575.0) * x * x * x * x * y * y * z * z * z + std::sqrt(24.609375) * x * x * x * x * z * z * z * z * z - std::sqrt(8.859375) * x * x * y * y * y * y * y * y * z + std::sqrt(1575.0) * x * x * y * y * y * y * z * z * z - std::sqrt(885.9375) * x * x * y * y * z * z * z * z * z + std::sqrt(0.5537109375) * y * y * y * y * y * y * y * y * z - std::sqrt(63.0) * y * y * y * y * y * y * z * z * z + std::sqrt(24.609375) * y * y * y * y * z * z * z * z * z) + e_1 * (-std::sqrt(35.4375) * x * x * x * x * x * x * z + std::sqrt(885.9375) * x * x * x * x * y * y * z - std::sqrt(885.9375) * x * x * x * x * z * z * z + std::sqrt(885.9375) * x * x * y * y * y * y * z + std::sqrt(31893.75) * x * x * y * y * z * z * z - std::sqrt(35.4375) * y * y * y * y * y * y * z - std::sqrt(885.9375) * y * y * y * y * z * z * z) + e_2 * (-std::sqrt(5537.109375) * x * x * x * x * z + std::sqrt(199335.9375) * x * x * y * y * z - std::sqrt(5537.109375) * y * y * y * y * z);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        pc_74[k] = e_0 * (std::sqrt(0.09228515625) * x * x * x * x * x * x * x * x * x - std::sqrt(1.4765625) * x * x * x * x * x * x * x * y * y - std::sqrt(18.087890625) * x * x * x * x * x * x * x * z * z - std::sqrt(9.228515625) * x * x * x * x * x * y * y * y * y + std::sqrt(452.197265625) * x * x * x * x * x * y * y * z * z + std::sqrt(147.65625) * x * x * x * x * x * z * z * z * z - std::sqrt(1.4765625) * x * x * x * y * y * y * y * y * y + std::sqrt(452.197265625) * x * x * x * y * y * y * y * z * z - std::sqrt(5315.625) * x * x * x * y * y * z * z * z * z + std::sqrt(0.09228515625) * x * y * y * y * y * y * y * y * y - std::sqrt(18.087890625) * x * y * y * y * y * y * y * z * z + std::sqrt(147.65625) * x * y * y * y * y * z * z * z * z) + e_1 * (std::sqrt(29.900390625) * x * x * x * x * x * x * x - std::sqrt(505.353515625) * x * x * x * x * x * y * y + std::sqrt(23.625) * x * x * x * x * x * z * z - std::sqrt(747.509765625) * x * x * x * y * y * y * y - std::sqrt(9450.0) * x * x * x * y * y * z * z + std::sqrt(2362.5) * x * x * x * z * z * z * z + std::sqrt(0.369140625) * x * y * y * y * y * y * y + std::sqrt(5315.625) * x * y * y * y * y * z * z - std::sqrt(21262.5) * x * y * y * z * z * z * z) + e_2 * (std::sqrt(2362.5) * x * x * x * x * x - std::sqrt(59062.5) * x * x * x * y * y + std::sqrt(14765.625) * x * x * x * z * z - std::sqrt(132890.625) * x * y * y * z * z) + e_3 * (std::sqrt(14765.625) * x * x * x - std::sqrt(132890.625) * x * y * y);

        pc_75[k] = e_0 * (-std::sqrt(0.9228515625) * x * x * x * x * x * x * x * x * z + std::sqrt(33.22265625) * x * x * x * x * x * x * y * y * z + std::sqrt(92.28515625) * x * x * x * x * x * x * z * z * z - std::sqrt(4521.97265625) * x * x * x * x * y * y * z * z * z - std::sqrt(33.22265625) * x * x * y * y * y * y * y * y * z + std::sqrt(4521.97265625) * x * x * y * y * y * y * z * z * z + std::sqrt(0.9228515625) * y * y * y * y * y * y * y * y * z - std::sqrt(92.28515625) * y * y * y * y * y * y * z * z * z) + e_1 * (std::sqrt(59.0625) * x * x * x * x * x * x * z - std::sqrt(13289.0625) * x * x * x * x * y * y * z + std::sqrt(5906.25) * x * x * x * x * z * z * z + std::sqrt(13289.0625) * x * x * y * y * y * y * z - std::sqrt(59.0625) * y * y * y * y * y * y * z - std::sqrt(5906.25) * y * y * y * y * z * z * z) + e_2 * (std::sqrt(13289.0625) * x * x * x * x * z + std::sqrt(53156.25) * x * x * z * z * z - std::sqrt(13289.0625) * y * y * y * y * z - std::sqrt(53156.25) * y * y * z * z * z) + e_3 * (std::sqrt(212625.0) * x * x * z - std::sqrt(212625.0) * y * y * z);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        pc_76[k] = e_0 * (-std::sqrt(0.15380859375) * x * x * x * x * x * x * x * x * x + std::sqrt(9.84375) * x * x * x * x * x * x * x * y * y + std::sqrt(15.380859375) * x * x * x * x * x * x * x * z * z - std::sqrt(15.380859375) * x * x * x * x * x * y * y * y * y - std::sqrt(1245.849609375) * x * x * x * x * x * y * y * z * z - std::sqrt(39.375) * x * x * x * y * y * y * y * y * y + std::sqrt(5552.490234375) * x * x * x * y * y * y * y * z * z + std::sqrt(1.38427734375) * x * y * y * y * y * y * y * y * y - std::sqrt(138.427734375) * x * y * y * y * y * y * y * z * z) + e_1 * (-std::sqrt(49.833984375) * x * x * x * x * x * x * x + std::sqrt(49.833984375) * x * x * x * x * x * y * y + std::sqrt(2214.84375) * x * x * x * x * x * z * z - std::sqrt(3460.693359375) * x * x * x * y * y * y * y + std::sqrt(8859.375) * x * x * x * y * y * z * z + std::sqrt(5.537109375) * x * y * y * y * y * y * y + std::sqrt(2214.84375) * x * y * y * y * y * z * z) + e_2 * (-std::sqrt(2214.84375) * x * x * x * x * x - std::sqrt(8859.375) * x * x * x * y * y + std::sqrt(79734.375) * x * x * x * z * z - std::sqrt(2214.84375) * x * y * y * y * y + std::sqrt(79734.375) * x * y * y * z * z) + e_3 * (-std::sqrt(8859.375) * x * x * x - std::sqrt(8859.375) * x * y * y + std::sqrt(141750.0) * x * z * z);

        pc_77[k] = e_0 * (std::sqrt(30.4541015625) * x * x * x * x * x * x * x * y * z - std::sqrt(3251.8212890625) * x * x * x * x * x * y * y * y * z + std::sqrt(2114.8681640625) * x * x * x * y * y * y * y * y * z - std::sqrt(84.5947265625) * x * y * y * y * y * y * y * y * z) + e_1 * (-std::sqrt(3045.41015625) * x * x * x * x * x * y * z - std::sqrt(12181.640625) * x * x * x * y * y * y * z - std::sqrt(3045.41015625) * x * y * y * y * y * y * z) + e_2 * (-std::sqrt(194906.25) * x * x * x * y * z - std::sqrt(194906.25) * x * y * y * y * z) + e_3 * (-std::sqrt(779625.0) * x * y * z);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        pc_78[k] = e_0 * (std::sqrt(81.2109375) * x * x * x * x * x * x * y * z * z - std::sqrt(8121.09375) * x * x * x * x * y * y * y * z * z + std::sqrt(2030.2734375) * x * x * y * y * y * y * y * z * z) + e_1 * (std::sqrt(81.2109375) * x * x * x * x * x * x * y - std::sqrt(8121.09375) * x * x * x * x * y * y * y - std::sqrt(18272.4609375) * x * x * x * x * y * z * z + std::sqrt(2030.2734375) * x * x * y * y * y * y * y - std::sqrt(8121.09375) * x * x * y * y * y * z * z + std::sqrt(2030.2734375) * y * y * y * y * y * z * z) + e_2 * (-std::sqrt(18272.4609375) * x * x * x * x * y - std::sqrt(8121.09375) * x * x * y * y * y - std::sqrt(292359.375) * x * x * y * z * z + std::sqrt(2030.2734375) * y * y * y * y * y + std::sqrt(32484.375) * y * y * y * z * z) + e_3 * (-std::sqrt(292359.375) * x * x * y + std::sqrt(32484.375) * y * y * y);

        pc_79[k] = e_0 * (-std::sqrt(2.0302734375) * x * x * x * x * x * x * x * y * z + std::sqrt(164.4521484375) * x * x * x * x * x * y * y * y * z + std::sqrt(32.484375) * x * x * x * x * x * y * z * z * z + std::sqrt(50.7568359375) * x * x * x * y * y * y * y * y * z - std::sqrt(3248.4375) * x * x * x * y * y * y * z * z * z - std::sqrt(50.7568359375) * x * y * y * y * y * y * y * y * z + std::sqrt(812.109375) * x * y * y * y * y * y * z * z * z) + e_1 * (std::sqrt(657.80859375) * x * x * x * x * x * y * z + std::sqrt(812.109375) * x * x * x * y * y * y * z - std::sqrt(12993.75) * x * x * x * y * z * z * z - std::sqrt(1827.24609375) * x * y * y * y * y * y * z + std::sqrt(12993.75) * x * y * y * y * z * z * z);

        pc_80[k] = e_0 * (-std::sqrt(12.181640625) * x * x * x * x * x * x * x * z * z + std::sqrt(986.712890625) * x * x * x * x * x * y * y * z * z + std::sqrt(5.4140625) * x * x * x * x * x * z * z * z * z + std::sqrt(304.541015625) * x * x * x * y * y * y * y * z * z - std::sqrt(541.40625) * x * x * x * y * y * z * z * z * z - std::sqrt(304.541015625) * x * y * y * y * y * y * y * z * z + std::sqrt(135.3515625) * x * y * y * y * y * z * z * z * z) + e_1 * (-std::sqrt(12.181640625) * x * x * x * x * x * x * x + std::sqrt(986.712890625) * x * x * x * x * x * y * y - std::sqrt(779.625) * x * x * x * x * x * z * z + std::sqrt(304.541015625) * x * x * x * y * y * y * y + std::sqrt(77962.5) * x * x * x * y * y * z * z - std::sqrt(304.541015625) * x * y * y * y * y * y * y - std::sqrt(19490.625) * x * y * y * y * y * z * z) + e_2 * (-std::sqrt(1218.1640625) * x * x * x * x * x + std::sqrt(121816.40625) * x * x * x * y * y - std::sqrt(30454.1015625) * x * y * y * y * y);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        pc_81[k] = e_0 * (-std::sqrt(2.0302734375) * x * x * x * x * x * x * x * x * z + std::sqrt(164.4521484375) * x * x * x * x * x * x * y * y * z + std::sqrt(32.484375) * x * x * x * x * x * x * z * z * z + std::sqrt(50.7568359375) * x * x * x * x * y * y * y * y * z - std::sqrt(3248.4375) * x * x * x * x * y * y * z * z * z - std::sqrt(50.7568359375) * x * x * y * y * y * y * y * y * z + std::sqrt(812.109375) * x * x * y * y * y * y * z * z * z) + e_1 * (-std::sqrt(99.4833984375) * x * x * x * x * x * x * z + std::sqrt(4111.3037109375) * x * x * x * x * y * y * z + std::sqrt(812.109375) * x * x * x * x * z * z * z + std::sqrt(456.8115234375) * x * x * y * y * y * y * z - std::sqrt(29235.9375) * x * x * y * y * z * z * z - std::sqrt(50.7568359375) * y * y * y * y * y * y * z + std::sqrt(812.109375) * y * y * y * y * z * z * z);

        pc_82[k] = e_0 * (std::sqrt(20.302734375) * x * x * x * x * x * x * x * z * z - std::sqrt(2456.630859375) * x * x * x * x * x * y * y * z * z + std::sqrt(4568.115234375) * x * x * x * y * y * y * y * z * z - std::sqrt(507.568359375) * x * y * y * y * y * y * y * z * z) + e_1 * (std::sqrt(20.302734375) * x * x * x * x * x * x * x - std::sqrt(2456.630859375) * x * x * x * x * x * y * y + std::sqrt(2030.2734375) * x * x * x * x * x * z * z + std::sqrt(4568.115234375) * x * x * x * y * y * y * y - std::sqrt(8121.09375) * x * x * x * y * y * z * z - std::sqrt(507.568359375) * x * y * y * y * y * y * y - std::sqrt(18272.4609375) * x * y * y * y * y * z * z) + e_2 * (std::sqrt(2030.2734375) * x * x * x * x * x - std::sqrt(8121.09375) * x * x * x * y * y + std::sqrt(32484.375) * x * x * x * z * z - std::sqrt(18272.4609375) * x * y * y * y * y - std::sqrt(292359.375) * x * y * y * z * z) + e_3 * (std::sqrt(32484.375) * x * x * x - std::sqrt(292359.375) * x * y * y);

        pc_83[k] = e_0 * (std::sqrt(3.3837890625) * x * x * x * x * x * x * x * x * z - std::sqrt(571.8603515625) * x * x * x * x * x * x * y * y * z + std::sqrt(4145.1416015625) * x * x * x * x * y * y * y * y * z - std::sqrt(761.3525390625) * x * x * y * y * y * y * y * y * z) + e_1 * (std::sqrt(761.3525390625) * x * x * x * x * x * x * z + std::sqrt(761.3525390625) * x * x * x * x * y * y * z - std::sqrt(761.3525390625) * x * x * y * y * y * y * z - std::sqrt(761.3525390625) * y * y * y * y * y * y * z) + e_2 * (std::sqrt(48726.5625) * x * x * x * x * z - std::sqrt(48726.5625) * y * y * y * y * z) + e_3 * (std::sqrt(194906.25) * x * x * z - std::sqrt(194906.25) * y * y * z);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        pc_84[k] = e_0 * (std::sqrt(2.537841796875) * x * x * x * x * x * x * x * x * y - std::sqrt(596.6748046875) * x * x * x * x * x * x * y * y * y + std::sqrt(1015.13671875) * x * x * x * x * y * y * y * y * y - std::sqrt(91.3623046875) * x * x * y * y * y * y * y * y * y + std::sqrt(0.281982421875) * y * y * y * y * y * y * y * y * y) + e_1 * (-std::sqrt(822.2607421875) * x * x * x * x * x * x * y - std::sqrt(2284.0576171875) * x * x * x * x * y * y * y - std::sqrt(91.3623046875) * x * x * y * y * y * y * y + std::sqrt(91.3623046875) * y * y * y * y * y * y * y) + e_2 * (-std::sqrt(82226.07421875) * x * x * x * x * y - std::sqrt(36544.921875) * x * x * y * y * y + std::sqrt(9136.23046875) * y * y * y * y * y) + e_3 * (-std::sqrt(584718.75) * x * x * y + std::sqrt(64968.75) * y * y * y);

        pc_85[k] = e_0 * (std::sqrt(6.767578125) * x * x * x * x * x * x * x * y * z - std::sqrt(1522.705078125) * x * x * x * x * x * y * y * y * z + std::sqrt(1522.705078125) * x * x * x * y * y * y * y * y * z - std::sqrt(6.767578125) * x * y * y * y * y * y * y * y * z) + e_1 * (-std::sqrt(3898.125) * x * x * x * x * x * y * z + std::sqrt(3898.125) * x * y * y * y * y * y * z) + e_2 * (-std::sqrt(97453.125) * x * x * x * y * z + std::sqrt(97453.125) * x * y * y * y * z);

        pc_86[k] = e_0 * (-std::sqrt(0.169189453125) * x * x * x * x * x * x * x * x * y + std::sqrt(33.1611328125) * x * x * x * x * x * x * y * y * y + std::sqrt(2.70703125) * x * x * x * x * x * x * y * z * z - std::sqrt(609.08203125) * x * x * x * x * y * y * y * z * z - std::sqrt(33.1611328125) * x * x * y * y * y * y * y * y * y + std::sqrt(609.08203125) * x * x * y * y * y * y * y * z * z + std::sqrt(0.169189453125) * y * y * y * y * y * y * y * y * y - std::sqrt(2.70703125) * y * y * y * y * y * y * y * z * z) + e_1 * (std::sqrt(54.8173828125) * x * x * x * x * x * x * y + std::sqrt(3806.7626953125) * x * x * x * x * y * y * y - std::sqrt(2436.328125) * x * x * x * x * y * z * z - std::sqrt(9264.1376953125) * x * x * y * y * y * y * y + std::sqrt(9745.3125) * x * x * y * y * y * z * z + std::sqrt(54.8173828125) * y * y * y * y * y * y * y - std::sqrt(97.453125) * y * y * y * y * y * z * z) + e_2 * (std::sqrt(15227.05078125) * x * x * x * x * y - std::sqrt(60908.203125) * x * x * y * y * y + std::sqrt(609.08203125) * y * y * y * y * y);
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

        pc_87[k] = e_0 * (-std::sqrt(1.01513671875) * x * x * x * x * x * x * x * x * z + std::sqrt(198.966796875) * x * x * x * x * x * x * y * y * z + std::sqrt(0.451171875) * x * x * x * x * x * x * z * z * z - std::sqrt(101.513671875) * x * x * x * x * y * y * z * z * z - std::sqrt(198.966796875) * x * x * y * y * y * y * y * y * z + std::sqrt(101.513671875) * x * x * y * y * y * y * z * z * z + std::sqrt(1.01513671875) * y * y * y * y * y * y * y * y * z - std::sqrt(0.451171875) * y * y * y * y * y * y * z * z * z) + e_1 * (-std::sqrt(146.1796875) * x * x * x * x * x * x * z + std::sqrt(32890.4296875) * x * x * x * x * y * y * z - std::sqrt(32890.4296875) * x * x * y * y * y * y * z + std::sqrt(146.1796875) * y * y * y * y * y * y * z);

        pc_88[k] = e_0 * (-std::sqrt(0.169189453125) * x * x * x * x * x * x * x * x * x + std::sqrt(33.1611328125) * x * x * x * x * x * x * x * y * y + std::sqrt(2.70703125) * x * x * x * x * x * x * x * z * z - std::sqrt(609.08203125) * x * x * x * x * x * y * y * z * z - std::sqrt(33.1611328125) * x * x * x * y * y * y * y * y * y + std::sqrt(609.08203125) * x * x * x * y * y * y * y * z * z + std::sqrt(0.169189453125) * x * y * y * y * y * y * y * y * y - std::sqrt(2.70703125) * x * y * y * y * y * y * y * z * z) + e_1 * (-std::sqrt(54.8173828125) * x * x * x * x * x * x * x + std::sqrt(9264.1376953125) * x * x * x * x * x * y * y + std::sqrt(97.453125) * x * x * x * x * x * z * z - std::sqrt(3806.7626953125) * x * x * x * y * y * y * y - std::sqrt(9745.3125) * x * x * x * y * y * z * z - std::sqrt(54.8173828125) * x * y * y * y * y * y * y + std::sqrt(2436.328125) * x * y * y * y * y * z * z) + e_2 * (-std::sqrt(609.08203125) * x * x * x * x * x + std::sqrt(60908.203125) * x * x * x * y * y - std::sqrt(15227.05078125) * x * y * y * y * y);

        pc_89[k] = e_0 * (std::sqrt(1.69189453125) * x * x * x * x * x * x * x * x * z - std::sqrt(433.125) * x * x * x * x * x * x * y * y * z + std::sqrt(1522.705078125) * x * x * x * x * y * y * y * y * z - std::sqrt(433.125) * x * x * y * y * y * y * y * y * z + std::sqrt(1.69189453125) * y * y * y * y * y * y * y * y * z) + e_1 * (std::sqrt(243.6328125) * x * x * x * x * x * x * z - std::sqrt(6090.8203125) * x * x * x * x * y * y * z - std::sqrt(6090.8203125) * x * x * y * y * y * y * z + std::sqrt(243.6328125) * y * y * y * y * y * y * z) + e_2 * (std::sqrt(6090.8203125) * x * x * x * x * z - std::sqrt(219269.53125) * x * x * y * y * z + std::sqrt(6090.8203125) * y * y * y * y * z);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ab_x, ab_y : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        pc_90[k] = e_0 * (std::sqrt(0.281982421875) * x * x * x * x * x * x * x * x * x - std::sqrt(91.3623046875) * x * x * x * x * x * x * x * y * y + std::sqrt(1015.13671875) * x * x * x * x * x * y * y * y * y - std::sqrt(596.6748046875) * x * x * x * y * y * y * y * y * y + std::sqrt(2.537841796875) * x * y * y * y * y * y * y * y * y) + e_1 * (std::sqrt(91.3623046875) * x * x * x * x * x * x * x - std::sqrt(91.3623046875) * x * x * x * x * x * y * y - std::sqrt(2284.0576171875) * x * x * x * y * y * y * y - std::sqrt(822.2607421875) * x * y * y * y * y * y * y) + e_2 * (std::sqrt(9136.23046875) * x * x * x * x * x - std::sqrt(36544.921875) * x * x * x * y * y - std::sqrt(82226.07421875) * x * y * y * y * y) + e_3 * (std::sqrt(64968.75) * x * x * x - std::sqrt(584718.75) * x * y * y);
    }

    // NOTE: the atom pairs beyond the reach of every pair of primitives have no
    // contribution and are set to zero.

    for (size_t m = 0; m < 91; m++)
    {
        auto *pv = values + m * nvalues;

        std::fill(pv + nmax, pv + nvalues, 0.0);
    }
}

}  // namespace simdovl
