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



#include "SimdOverlapRecIH.hpp"

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
compute_ih_overlap(double               *values,
                   const size_t          nvalues,
                   const CBasisFunction &bra,
                   const CBasisFunction &ket,
                   const CSimdMatrix    &coordinates,
                   const double          threshold) -> void
{
    if ((bra.get_angular_momentum() != 6) || (ket.get_angular_momentum() != 5))
    {
        errors::assertMsgCritical(
            false, std::string("SimdOverlapRecIH.compute_ih_overlap: Basis functions must be of angular momenta six and five"));
    }

    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("SimdOverlapRecIH.compute_ih_overlap: Number of values exceeds number of atom pairs"));
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

    auto buffer = simdfunc::make_primitive_buffer(dimensions, 6);

    if (buffer.number_of_columns() == 0)
    {
        std::fill(values, values + 143 * nvalues, 0.0);

        return;
    }

    const auto nmax = buffer.number_of_columns();

    auto *pe_0 = buffer.data(0);
    auto *pe_1 = buffer.data(1);
    auto *pe_2 = buffer.data(2);
    auto *pe_3 = buffer.data(3);
    auto *pe_4 = buffer.data(4);
    auto *pe_5 = buffer.data(5);

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

        const auto f_0 = fbase * fal * fal * fal * fal * fal * fal * fbe * fbe * fbe * fbe * fbe;

        const auto f_1 = fbase * fal * fal * fal * fal * fal * fbe * fbe * fbe * fbe * fh;

        const auto f_2 = fbase * fal * fal * fal * fal * fbe * fbe * fbe * fh * fh;

        const auto f_3 = fbase * fal * fal * fal * fbe * fbe * fh * fh * fh;

        const auto f_4 = fbase * fal * fal * fbe * fh * fh * fh * fh;

        const auto f_5 = fbase * fal * fh * fh * fh * fh * fh;

        // NOTE: the exponential depends on the pair of primitives alone, so it is
        // evaluated once and shared by the prefactors of all terms.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ab_2 : simd::cache_line_size())
        for (size_t k = 0; k < ncols; k++)
        {
            const auto fss = std::exp(-fmu * ab_2[k]);

            pe_0[k] += f_0 * fss;
            pe_1[k] += f_1 * fss;
            pe_2[k] += f_2 * fss;
            pe_3[k] += f_3 * fss;
            pe_4[k] += f_4 * fss;
            pe_5[k] += f_5 * fss;
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
    auto *pc_91 = values + 91 * nvalues;
    auto *pc_92 = values + 92 * nvalues;
    auto *pc_93 = values + 93 * nvalues;
    auto *pc_94 = values + 94 * nvalues;
    auto *pc_95 = values + 95 * nvalues;
    auto *pc_96 = values + 96 * nvalues;
    auto *pc_97 = values + 97 * nvalues;
    auto *pc_98 = values + 98 * nvalues;
    auto *pc_99 = values + 99 * nvalues;
    auto *pc_100 = values + 100 * nvalues;
    auto *pc_101 = values + 101 * nvalues;
    auto *pc_102 = values + 102 * nvalues;
    auto *pc_103 = values + 103 * nvalues;
    auto *pc_104 = values + 104 * nvalues;
    auto *pc_105 = values + 105 * nvalues;
    auto *pc_106 = values + 106 * nvalues;
    auto *pc_107 = values + 107 * nvalues;
    auto *pc_108 = values + 108 * nvalues;
    auto *pc_109 = values + 109 * nvalues;
    auto *pc_110 = values + 110 * nvalues;
    auto *pc_111 = values + 111 * nvalues;
    auto *pc_112 = values + 112 * nvalues;
    auto *pc_113 = values + 113 * nvalues;
    auto *pc_114 = values + 114 * nvalues;
    auto *pc_115 = values + 115 * nvalues;
    auto *pc_116 = values + 116 * nvalues;
    auto *pc_117 = values + 117 * nvalues;
    auto *pc_118 = values + 118 * nvalues;
    auto *pc_119 = values + 119 * nvalues;
    auto *pc_120 = values + 120 * nvalues;
    auto *pc_121 = values + 121 * nvalues;
    auto *pc_122 = values + 122 * nvalues;
    auto *pc_123 = values + 123 * nvalues;
    auto *pc_124 = values + 124 * nvalues;
    auto *pc_125 = values + 125 * nvalues;
    auto *pc_126 = values + 126 * nvalues;
    auto *pc_127 = values + 127 * nvalues;
    auto *pc_128 = values + 128 * nvalues;
    auto *pc_129 = values + 129 * nvalues;
    auto *pc_130 = values + 130 * nvalues;
    auto *pc_131 = values + 131 * nvalues;
    auto *pc_132 = values + 132 * nvalues;
    auto *pc_133 = values + 133 * nvalues;
    auto *pc_134 = values + 134 * nvalues;
    auto *pc_135 = values + 135 * nvalues;
    auto *pc_136 = values + 136 * nvalues;
    auto *pc_137 = values + 137 * nvalues;
    auto *pc_138 = values + 138 * nvalues;
    auto *pc_139 = values + 139 * nvalues;
    auto *pc_140 = values + 140 * nvalues;
    auto *pc_141 = values + 141 * nvalues;
    auto *pc_142 = values + 142 * nvalues;

    // NOTE: the components are formed in 134 loops, as the vectorizer runs out
    // of registers with all of them in one. Only the prefactors and the vector
    // between the atoms are loaded by more than one loop.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        pc_0[k] = e_0 * (std::sqrt(199.85504150390625) * x * x * x * x * x * x * x * x * x * y * y - std::sqrt(5684.765625) * x * x * x * x * x * x * x * y * y * y * y + std::sqrt(12367.918212890625) * x * x * x * x * x * y * y * y * y * y * y - std::sqrt(1421.19140625) * x * x * x * y * y * y * y * y * y * y * y + std::sqrt(7.99420166015625) * x * y * y * y * y * y * y * y * y * y * y) + e_1 * (std::sqrt(199.85504150390625) * x * x * x * x * x * x * x * x * x + std::sqrt(3197.6806640625) * x * x * x * x * x * x * x * y * y + std::sqrt(7194.781494140625) * x * x * x * x * x * y * y * y * y + std::sqrt(3197.6806640625) * x * x * x * y * y * y * y * y * y + std::sqrt(199.85504150390625) * x * y * y * y * y * y * y * y * y) + e_2 * (std::sqrt(79942.0166015625) * x * x * x * x * x * x * x + std::sqrt(719478.1494140625) * x * x * x * x * x * y * y + std::sqrt(719478.1494140625) * x * x * x * y * y * y * y + std::sqrt(79942.0166015625) * x * y * y * y * y * y * y) + e_3 * (std::sqrt(5116289.0625) * x * x * x * x * x + std::sqrt(20465156.25) * x * x * x * y * y + std::sqrt(5116289.0625) * x * y * y * y * y) + e_4 * (std::sqrt(46046601.5625) * x * x * x + std::sqrt(46046601.5625) * x * y * y) + e_5 * (std::sqrt(29469825.0) * x);

        pc_1[k] = e_0 * (std::sqrt(1279.072265625) * x * x * x * x * x * x * x * x * y * y * z - std::sqrt(24018.134765625) * x * x * x * x * x * x * y * y * y * y * z + std::sqrt(24018.134765625) * x * x * x * x * y * y * y * y * y * y * z - std::sqrt(1279.072265625) * x * x * y * y * y * y * y * y * y * y * z) + e_1 * (std::sqrt(1279.072265625) * x * x * x * x * x * x * x * x * z + std::sqrt(5116.2890625) * x * x * x * x * x * x * y * y * z - std::sqrt(5116.2890625) * x * x * y * y * y * y * y * y * z - std::sqrt(1279.072265625) * y * y * y * y * y * y * y * y * z) + e_2 * (std::sqrt(287791.259765625) * x * x * x * x * x * x * z + std::sqrt(287791.259765625) * x * x * x * x * y * y * z - std::sqrt(287791.259765625) * x * x * y * y * y * y * z - std::sqrt(287791.259765625) * y * y * y * y * y * y * z) + e_3 * (std::sqrt(8186062.5) * x * x * x * x * z - std::sqrt(8186062.5) * y * y * y * y * z) + e_4 * (std::sqrt(18418640.625) * x * x * z - std::sqrt(18418640.625) * y * y * z);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        pc_2[k] = e_0 * (-std::sqrt(39.97100830078125) * x * x * x * x * x * x * x * x * x * y * y + std::sqrt(284.23828125) * x * x * x * x * x * x * x * y * y * y * y + std::sqrt(2558.14453125) * x * x * x * x * x * x * x * y * y * z * z + std::sqrt(96.719970703125) * x * x * x * x * x * y * y * y * y * y * y - std::sqrt(34392.83203125) * x * x * x * x * x * y * y * y * y * z * z - std::sqrt(126.328125) * x * x * x * y * y * y * y * y * y * y * y + std::sqrt(11401.11328125) * x * x * x * y * y * y * y * y * y * z * z + std::sqrt(4.44122314453125) * x * y * y * y * y * y * y * y * y * y * y - std::sqrt(284.23828125) * x * y * y * y * y * y * y * y * y * z * z) + e_1 * (-std::sqrt(39.97100830078125) * x * x * x * x * x * x * x * x * x - std::sqrt(5755.8251953125) * x * x * x * x * x * x * x * y * y + std::sqrt(2558.14453125) * x * x * x * x * x * x * x * z * z + std::sqrt(99927.52075195312) * x * x * x * x * x * y * y * y * y - std::sqrt(2558.14453125) * x * x * x * x * x * y * y * z * z - std::sqrt(12009.0673828125) * x * x * x * y * y * y * y * y * y - std::sqrt(63953.61328125) * x * x * x * y * y * y * y * z * z + std::sqrt(1958.5794067382812) * x * y * y * y * y * y * y * y * y - std::sqrt(23023.30078125) * x * y * y * y * y * y * y * z * z) + e_2 * (-std::sqrt(15988.4033203125) * x * x * x * x * x * x * x + std::sqrt(15988.4033203125) * x * x * x * x * x * y * y + std::sqrt(255814.453125) * x * x * x * x * x * z * z + std::sqrt(399710.0830078125) * x * x * x * y * y * y * y - std::sqrt(1023257.8125) * x * x * x * y * y * z * z + std::sqrt(143895.6298828125) * x * y * y * y * y * y * y - std::sqrt(2302330.078125) * x * y * y * y * y * z * z) + e_3 * (-std::sqrt(454781.25) * x * x * x * x * x + std::sqrt(1819125.0) * x * x * x * y * y + std::sqrt(1819125.0) * x * x * x * z * z + std::sqrt(4093031.25) * x * y * y * y * y - std::sqrt(16372125.0) * x * y * y * z * z) + e_4 * (-std::sqrt(1023257.8125) * x * x * x + std::sqrt(9209320.3125) * x * y * y);
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

        pc_3[k] = e_0 * (-std::sqrt(426.357421875) * x * x * x * x * x * x * x * x * y * y * z + std::sqrt(2321.279296875) * x * x * x * x * x * x * y * y * y * y * z + std::sqrt(1705.4296875) * x * x * x * x * x * x * y * y * z * z * z + std::sqrt(2321.279296875) * x * x * x * x * y * y * y * y * y * y * z - std::sqrt(18949.21875) * x * x * x * x * y * y * y * y * z * z * z - std::sqrt(426.357421875) * x * x * y * y * y * y * y * y * y * y * z + std::sqrt(1705.4296875) * x * x * y * y * y * y * y * y * z * z * z) + e_1 * (-std::sqrt(426.357421875) * x * x * x * x * x * x * x * x * z - std::sqrt(27286.875) * x * x * x * x * x * x * y * y * z + std::sqrt(1705.4296875) * x * x * x * x * x * x * z * z * z + std::sqrt(1065893.5546875) * x * x * x * x * y * y * y * y * z - std::sqrt(42635.7421875) * x * x * x * x * y * y * z * z * z - std::sqrt(27286.875) * x * x * y * y * y * y * y * y * z - std::sqrt(42635.7421875) * x * x * y * y * y * y * z * z * z - std::sqrt(426.357421875) * y * y * y * y * y * y * y * y * z + std::sqrt(1705.4296875) * y * y * y * y * y * y * z * z * z) + e_2 * (-std::sqrt(95930.419921875) * x * x * x * x * x * x * z + std::sqrt(2398260.498046875) * x * x * x * x * y * y * z + std::sqrt(42635.7421875) * x * x * x * x * z * z * z + std::sqrt(2398260.498046875) * x * x * y * y * y * y * z - std::sqrt(1534886.71875) * x * x * y * y * z * z * z - std::sqrt(95930.419921875) * y * y * y * y * y * y * z + std::sqrt(42635.7421875) * y * y * y * y * z * z * z) + e_3 * (-std::sqrt(682171.875) * x * x * x * x * z + std::sqrt(24558187.5) * x * x * y * y * z - std::sqrt(682171.875) * y * y * y * y * z);
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

        pc_4[k] = e_0 * (std::sqrt(3.8067626953125) * x * x * x * x * x * x * x * x * x * y * y - std::sqrt(6.767578125) * x * x * x * x * x * x * x * y * y * y * y - std::sqrt(548.173828125) * x * x * x * x * x * x * x * y * y * z * z - std::sqrt(82.90283203125) * x * x * x * x * x * y * y * y * y * y * y + std::sqrt(2984.501953125) * x * x * x * x * x * y * y * y * y * z * z + std::sqrt(243.6328125) * x * x * x * x * x * y * y * z * z * z * z - std::sqrt(6.767578125) * x * x * x * y * y * y * y * y * y * y * y + std::sqrt(2984.501953125) * x * x * x * y * y * y * y * y * y * z * z - std::sqrt(2707.03125) * x * x * x * y * y * y * y * z * z * z * z + std::sqrt(3.8067626953125) * x * y * y * y * y * y * y * y * y * y * y - std::sqrt(548.173828125) * x * y * y * y * y * y * y * y * y * z * z + std::sqrt(243.6328125) * x * y * y * y * y * y * y * z * z * z * z) + e_1 * (std::sqrt(3.8067626953125) * x * x * x * x * x * x * x * x * x + std::sqrt(974.53125) * x * x * x * x * x * x * x * y * y - std::sqrt(548.173828125) * x * x * x * x * x * x * x * z * z - std::sqrt(18653.13720703125) * x * x * x * x * x * y * y * y * y - std::sqrt(4933.564453125) * x * x * x * x * x * y * y * z * z + std::sqrt(243.6328125) * x * x * x * x * x * z * z * z * z - std::sqrt(11938.0078125) * x * x * x * y * y * y * y * y * y + std::sqrt(1110052.001953125) * x * x * x * y * y * y * y * z * z - std::sqrt(24363.28125) * x * x * x * y * y * z * z * z * z + std::sqrt(3201.4874267578125) * x * y * y * y * y * y * y * y * y - std::sqrt(158422.236328125) * x * y * y * y * y * y * y * z * z + std::sqrt(6090.8203125) * x * y * y * y * y * z * z * z * z) + e_2 * (std::sqrt(1522.705078125) * x * x * x * x * x * x * x - std::sqrt(13704.345703125) * x * x * x * x * x * y * y - std::sqrt(54817.3828125) * x * x * x * x * x * z * z - std::sqrt(951690.673828125) * x * x * x * y * y * y * y + std::sqrt(5481738.28125) * x * x * x * y * y * z * z + std::sqrt(184247.314453125) * x * y * y * y * y * y * y - std::sqrt(1370434.5703125) * x * y * y * y * y * z * z) + e_3 * (std::sqrt(24363.28125) * x * x * x * x * x - std::sqrt(2436328.125) * x * x * x * y * y + std::sqrt(609082.03125) * x * y * y * y * y);
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

        pc_5[k] = e_0 * (std::sqrt(57.1014404296875) * x * x * x * x * x * x * x * x * x * y * z - std::sqrt(101.513671875) * x * x * x * x * x * x * x * y * y * y * z - std::sqrt(406.0546875) * x * x * x * x * x * x * x * y * z * z * z - std::sqrt(1243.54248046875) * x * x * x * x * x * y * y * y * y * y * z + std::sqrt(2210.7421875) * x * x * x * x * x * y * y * y * z * z * z + std::sqrt(16.2421875) * x * x * x * x * x * y * z * z * z * z * z - std::sqrt(101.513671875) * x * x * x * y * y * y * y * y * y * y * z + std::sqrt(2210.7421875) * x * x * x * y * y * y * y * y * z * z * z - std::sqrt(180.46875) * x * x * x * y * y * y * z * z * z * z * z + std::sqrt(57.1014404296875) * x * y * y * y * y * y * y * y * y * y * z - std::sqrt(406.0546875) * x * y * y * y * y * y * y * y * z * z * z + std::sqrt(16.2421875) * x * y * y * y * y * y * z * z * z * z * z) + e_1 * (std::sqrt(32890.4296875) * x * x * x * x * x * x * x * y * z - std::sqrt(179070.1171875) * x * x * x * x * x * y * y * y * z - std::sqrt(58471.875) * x * x * x * x * x * y * z * z * z - std::sqrt(179070.1171875) * x * x * x * y * y * y * y * y * z + std::sqrt(649687.5) * x * x * x * y * y * y * z * z * z + std::sqrt(32890.4296875) * x * y * y * y * y * y * y * y * z - std::sqrt(58471.875) * x * y * y * y * y * y * z * z * z) + e_2 * (std::sqrt(822260.7421875) * x * x * x * x * x * y * z - std::sqrt(9136230.46875) * x * x * x * y * y * y * z + std::sqrt(822260.7421875) * x * y * y * y * y * y * z);
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

        pc_6[k] = e_0 * (std::sqrt(3.8067626953125) * x * x * x * x * x * x * x * x * x * x * y - std::sqrt(6.767578125) * x * x * x * x * x * x * x * x * y * y * y - std::sqrt(548.173828125) * x * x * x * x * x * x * x * x * y * z * z - std::sqrt(82.90283203125) * x * x * x * x * x * x * y * y * y * y * y + std::sqrt(2984.501953125) * x * x * x * x * x * x * y * y * y * z * z + std::sqrt(243.6328125) * x * x * x * x * x * x * y * z * z * z * z - std::sqrt(6.767578125) * x * x * x * x * y * y * y * y * y * y * y + std::sqrt(2984.501953125) * x * x * x * x * y * y * y * y * y * z * z - std::sqrt(2707.03125) * x * x * x * x * y * y * y * z * z * z * z + std::sqrt(3.8067626953125) * x * x * y * y * y * y * y * y * y * y * y - std::sqrt(548.173828125) * x * x * y * y * y * y * y * y * y * z * z + std::sqrt(243.6328125) * x * x * y * y * y * y * y * z * z * z * z) + e_1 * (std::sqrt(3201.4874267578125) * x * x * x * x * x * x * x * x * y - std::sqrt(11938.0078125) * x * x * x * x * x * x * y * y * y - std::sqrt(158422.236328125) * x * x * x * x * x * x * y * z * z - std::sqrt(18653.13720703125) * x * x * x * x * y * y * y * y * y + std::sqrt(1110052.001953125) * x * x * x * x * y * y * y * z * z + std::sqrt(6090.8203125) * x * x * x * x * y * z * z * z * z + std::sqrt(974.53125) * x * x * y * y * y * y * y * y * y - std::sqrt(4933.564453125) * x * x * y * y * y * y * y * z * z - std::sqrt(24363.28125) * x * x * y * y * y * z * z * z * z + std::sqrt(3.8067626953125) * y * y * y * y * y * y * y * y * y - std::sqrt(548.173828125) * y * y * y * y * y * y * y * z * z + std::sqrt(243.6328125) * y * y * y * y * y * z * z * z * z) + e_2 * (std::sqrt(184247.314453125) * x * x * x * x * x * x * y - std::sqrt(951690.673828125) * x * x * x * x * y * y * y - std::sqrt(1370434.5703125) * x * x * x * x * y * z * z - std::sqrt(13704.345703125) * x * x * y * y * y * y * y + std::sqrt(5481738.28125) * x * x * y * y * y * z * z + std::sqrt(1522.705078125) * y * y * y * y * y * y * y - std::sqrt(54817.3828125) * y * y * y * y * y * z * z) + e_3 * (std::sqrt(609082.03125) * x * x * x * x * y - std::sqrt(2436328.125) * x * x * y * y * y + std::sqrt(24363.28125) * y * y * y * y * y);
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

        pc_7[k] = e_0 * (-std::sqrt(106.58935546875) * x * x * x * x * x * x * x * x * x * y * z + std::sqrt(1184.326171875) * x * x * x * x * x * x * x * y * y * y * z + std::sqrt(426.357421875) * x * x * x * x * x * x * x * y * z * z * z - std::sqrt(8006.044921875) * x * x * x * x * x * y * y * y * z * z * z - std::sqrt(1184.326171875) * x * x * x * y * y * y * y * y * y * y * z + std::sqrt(8006.044921875) * x * x * x * y * y * y * y * y * z * z * z + std::sqrt(106.58935546875) * x * y * y * y * y * y * y * y * y * y * z - std::sqrt(426.357421875) * x * y * y * y * y * y * y * y * z * z * z) + e_1 * (-std::sqrt(42635.7421875) * x * x * x * x * x * x * x * y * z + std::sqrt(206356.9921875) * x * x * x * x * x * y * y * y * z + std::sqrt(27286.875) * x * x * x * x * x * y * z * z * z - std::sqrt(206356.9921875) * x * x * x * y * y * y * y * y * z + std::sqrt(42635.7421875) * x * y * y * y * y * y * y * y * z - std::sqrt(27286.875) * x * y * y * y * y * y * z * z * z) + e_2 * (-std::sqrt(1534886.71875) * x * x * x * x * x * y * z + std::sqrt(682171.875) * x * x * x * y * z * z * z + std::sqrt(1534886.71875) * x * y * y * y * y * y * z - std::sqrt(682171.875) * x * y * y * y * z * z * z) + e_3 * (-std::sqrt(10914750.0) * x * x * x * y * z + std::sqrt(10914750.0) * x * y * y * y * z);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        pc_8[k] = e_0 * (-std::sqrt(4.44122314453125) * x * x * x * x * x * x * x * x * x * x * y + std::sqrt(126.328125) * x * x * x * x * x * x * x * x * y * y * y + std::sqrt(284.23828125) * x * x * x * x * x * x * x * x * y * z * z - std::sqrt(96.719970703125) * x * x * x * x * x * x * y * y * y * y * y - std::sqrt(11401.11328125) * x * x * x * x * x * x * y * y * y * z * z - std::sqrt(284.23828125) * x * x * x * x * y * y * y * y * y * y * y + std::sqrt(34392.83203125) * x * x * x * x * y * y * y * y * y * z * z + std::sqrt(39.97100830078125) * x * x * y * y * y * y * y * y * y * y * y - std::sqrt(2558.14453125) * x * x * y * y * y * y * y * y * y * z * z) + e_1 * (-std::sqrt(1958.5794067382812) * x * x * x * x * x * x * x * x * y + std::sqrt(12009.0673828125) * x * x * x * x * x * x * y * y * y + std::sqrt(23023.30078125) * x * x * x * x * x * x * y * z * z - std::sqrt(99927.52075195312) * x * x * x * x * y * y * y * y * y + std::sqrt(63953.61328125) * x * x * x * x * y * y * y * z * z + std::sqrt(5755.8251953125) * x * x * y * y * y * y * y * y * y + std::sqrt(2558.14453125) * x * x * y * y * y * y * y * z * z + std::sqrt(39.97100830078125) * y * y * y * y * y * y * y * y * y - std::sqrt(2558.14453125) * y * y * y * y * y * y * y * z * z) + e_2 * (-std::sqrt(143895.6298828125) * x * x * x * x * x * x * y - std::sqrt(399710.0830078125) * x * x * x * x * y * y * y + std::sqrt(2302330.078125) * x * x * x * x * y * z * z - std::sqrt(15988.4033203125) * x * x * y * y * y * y * y + std::sqrt(1023257.8125) * x * x * y * y * y * z * z + std::sqrt(15988.4033203125) * y * y * y * y * y * y * y - std::sqrt(255814.453125) * y * y * y * y * y * z * z) + e_3 * (-std::sqrt(4093031.25) * x * x * x * x * y - std::sqrt(1819125.0) * x * x * y * y * y + std::sqrt(16372125.0) * x * x * y * z * z + std::sqrt(454781.25) * y * y * y * y * y - std::sqrt(1819125.0) * y * y * y * z * z) + e_4 * (-std::sqrt(9209320.3125) * x * x * y + std::sqrt(1023257.8125) * y * y * y);

        pc_9[k] = e_0 * (std::sqrt(79.9420166015625) * x * x * x * x * x * x * x * x * x * y * z - std::sqrt(6963.837890625) * x * x * x * x * x * x * x * y * y * y * z + std::sqrt(38691.93603515625) * x * x * x * x * x * y * y * y * y * y * z - std::sqrt(6963.837890625) * x * x * x * y * y * y * y * y * y * y * z + std::sqrt(79.9420166015625) * x * y * y * y * y * y * y * y * y * y * z) + e_1 * (std::sqrt(5116.2890625) * x * x * x * x * x * x * x * y * z + std::sqrt(46046.6015625) * x * x * x * x * x * y * y * y * z + std::sqrt(46046.6015625) * x * x * x * y * y * y * y * y * z + std::sqrt(5116.2890625) * x * y * y * y * y * y * y * y * z) + e_2 * (std::sqrt(1151165.0390625) * x * x * x * x * x * y * z + std::sqrt(4604660.15625) * x * x * x * y * y * y * z + std::sqrt(1151165.0390625) * x * y * y * y * y * y * z) + e_3 * (std::sqrt(32744250.0) * x * x * x * y * z + std::sqrt(32744250.0) * x * y * y * y * z) + e_4 * (std::sqrt(73674562.5) * x * y * z);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        pc_10[k] = e_0 * (std::sqrt(7.99420166015625) * x * x * x * x * x * x * x * x * x * x * y - std::sqrt(1421.19140625) * x * x * x * x * x * x * x * x * y * y * y + std::sqrt(12367.918212890625) * x * x * x * x * x * x * y * y * y * y * y - std::sqrt(5684.765625) * x * x * x * x * y * y * y * y * y * y * y + std::sqrt(199.85504150390625) * x * x * y * y * y * y * y * y * y * y * y) + e_1 * (std::sqrt(199.85504150390625) * x * x * x * x * x * x * x * x * y + std::sqrt(3197.6806640625) * x * x * x * x * x * x * y * y * y + std::sqrt(7194.781494140625) * x * x * x * x * y * y * y * y * y + std::sqrt(3197.6806640625) * x * x * y * y * y * y * y * y * y + std::sqrt(199.85504150390625) * y * y * y * y * y * y * y * y * y) + e_2 * (std::sqrt(79942.0166015625) * x * x * x * x * x * x * y + std::sqrt(719478.1494140625) * x * x * x * x * y * y * y + std::sqrt(719478.1494140625) * x * x * y * y * y * y * y + std::sqrt(79942.0166015625) * y * y * y * y * y * y * y) + e_3 * (std::sqrt(5116289.0625) * x * x * x * x * y + std::sqrt(20465156.25) * x * x * y * y * y + std::sqrt(5116289.0625) * y * y * y * y * y) + e_4 * (std::sqrt(46046601.5625) * x * x * y + std::sqrt(46046601.5625) * y * y * y) + e_5 * (std::sqrt(29469825.0) * y);

        pc_11[k] = e_0 * (std::sqrt(1665.4586791992188) * x * x * x * x * x * x * x * x * y * y * z - std::sqrt(26647.3388671875) * x * x * x * x * x * x * y * y * y * y * z + std::sqrt(32243.280029296875) * x * x * x * x * y * y * y * y * y * y * z - std::sqrt(1065.8935546875) * x * x * y * y * y * y * y * y * y * y * z + std::sqrt(2.66473388671875) * y * y * y * y * y * y * y * y * y * y * z) + e_1 * (std::sqrt(1665.4586791992188) * x * x * x * x * x * x * x * x * z + std::sqrt(26647.3388671875) * x * x * x * x * x * x * y * y * z + std::sqrt(59956.512451171875) * x * x * x * x * y * y * y * y * z + std::sqrt(26647.3388671875) * x * x * y * y * y * y * y * y * z + std::sqrt(1665.4586791992188) * y * y * y * y * y * y * y * y * z) + e_2 * (std::sqrt(426357.421875) * x * x * x * x * x * x * z + std::sqrt(3837216.796875) * x * x * x * x * y * y * z + std::sqrt(3837216.796875) * x * x * y * y * y * y * z + std::sqrt(426357.421875) * y * y * y * y * y * y * z) + e_3 * (std::sqrt(15348867.1875) * x * x * x * x * z + std::sqrt(61395468.75) * x * x * y * y * z + std::sqrt(15348867.1875) * y * y * y * y * z) + e_4 * (std::sqrt(61395468.75) * x * x * z + std::sqrt(61395468.75) * y * y * z) + e_5 * (std::sqrt(9823275.0) * z);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        pc_12[k] = e_0 * (std::sqrt(10658.935546875) * x * x * x * x * x * x * x * y * y * z * z - std::sqrt(95930.419921875) * x * x * x * x * x * y * y * y * y * z * z + std::sqrt(51589.248046875) * x * x * x * y * y * y * y * y * y * z * z - std::sqrt(426.357421875) * x * y * y * y * y * y * y * y * y * z * z) + e_1 * (std::sqrt(10658.935546875) * x * x * x * x * x * x * x * y * y + std::sqrt(10658.935546875) * x * x * x * x * x * x * x * z * z - std::sqrt(95930.419921875) * x * x * x * x * x * y * y * y * y + std::sqrt(95930.419921875) * x * x * x * x * x * y * y * z * z + std::sqrt(51589.248046875) * x * x * x * y * y * y * y * y * y + std::sqrt(95930.419921875) * x * x * x * y * y * y * y * z * z - std::sqrt(426.357421875) * x * y * y * y * y * y * y * y * y + std::sqrt(10658.935546875) * x * y * y * y * y * y * y * z * z) + e_2 * (std::sqrt(10658.935546875) * x * x * x * x * x * x * x + std::sqrt(95930.419921875) * x * x * x * x * x * y * y + std::sqrt(1534886.71875) * x * x * x * x * x * z * z + std::sqrt(95930.419921875) * x * x * x * y * y * y * y + std::sqrt(6139546.875) * x * x * x * y * y * z * z + std::sqrt(10658.935546875) * x * y * y * y * y * y * y + std::sqrt(1534886.71875) * x * y * y * y * y * z * z) + e_3 * (std::sqrt(1534886.71875) * x * x * x * x * x + std::sqrt(6139546.875) * x * x * x * y * y + std::sqrt(24558187.5) * x * x * x * z * z + std::sqrt(1534886.71875) * x * y * y * y * y + std::sqrt(24558187.5) * x * y * y * z * z) + e_4 * (std::sqrt(24558187.5) * x * x * x + std::sqrt(24558187.5) * x * y * y + std::sqrt(24558187.5) * x * z * z) + e_5 * (std::sqrt(24558187.5) * x);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        pc_13[k] = e_0 * (-std::sqrt(333.09173583984375) * x * x * x * x * x * x * x * x * y * y * z + std::sqrt(592.1630859375) * x * x * x * x * x * x * y * y * y * y * z + std::sqrt(21317.87109375) * x * x * x * x * x * x * y * y * z * z * z + std::sqrt(716.517333984375) * x * x * x * x * y * y * y * y * y * y * z - std::sqrt(116063.96484375) * x * x * x * x * y * y * y * y * z * z * z - std::sqrt(213.1787109375) * x * x * y * y * y * y * y * y * y * y * z + std::sqrt(16012.08984375) * x * x * y * y * y * y * y * y * z * z * z + std::sqrt(1.48040771484375) * y * y * y * y * y * y * y * y * y * y * z - std::sqrt(94.74609375) * y * y * y * y * y * y * y * y * z * z * z) + e_1 * (-std::sqrt(333.09173583984375) * x * x * x * x * x * x * x * x * z + std::sqrt(5329.4677734375) * x * x * x * x * x * x * y * y * z + std::sqrt(21317.87109375) * x * x * x * x * x * x * z * z * z - std::sqrt(65285.980224609375) * x * x * x * x * y * y * y * y * z + std::sqrt(21317.87109375) * x * x * x * x * y * y * z * z * z + std::sqrt(17267.4755859375) * x * x * y * y * y * y * y * y * z - std::sqrt(21317.87109375) * x * x * y * y * y * y * z * z * z + std::sqrt(119.91302490234375) * y * y * y * y * y * y * y * y * z - std::sqrt(21317.87109375) * y * y * y * y * y * y * z * z * z) + e_2 * (std::sqrt(1364343.75) * x * x * x * x * z * z * z - std::sqrt(1364343.75) * y * y * y * y * z * z * z) + e_3 * (std::sqrt(1364343.75) * x * x * x * x * z + std::sqrt(5457375.0) * x * x * z * z * z - std::sqrt(1364343.75) * y * y * y * y * z - std::sqrt(5457375.0) * y * y * z * z * z) + e_4 * (std::sqrt(12279093.75) * x * x * z - std::sqrt(12279093.75) * y * y * z);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        pc_14[k] = e_0 * (-std::sqrt(3552.978515625) * x * x * x * x * x * x * x * y * y * z * z + std::sqrt(3552.978515625) * x * x * x * x * x * y * y * y * y * z * z + std::sqrt(14211.9140625) * x * x * x * x * x * y * y * z * z * z * z + std::sqrt(11511.650390625) * x * x * x * y * y * y * y * y * y * z * z - std::sqrt(56847.65625) * x * x * x * y * y * y * y * z * z * z * z - std::sqrt(142.119140625) * x * y * y * y * y * y * y * y * y * z * z + std::sqrt(568.4765625) * x * y * y * y * y * y * y * z * z * z * z) + e_1 * (-std::sqrt(3552.978515625) * x * x * x * x * x * x * x * y * y - std::sqrt(3552.978515625) * x * x * x * x * x * x * x * z * z + std::sqrt(3552.978515625) * x * x * x * x * x * y * y * y * y - std::sqrt(31976.806640625) * x * x * x * x * x * y * y * z * z + std::sqrt(14211.9140625) * x * x * x * x * x * z * z * z * z + std::sqrt(11511.650390625) * x * x * x * y * y * y * y * y * y + std::sqrt(600453.369140625) * x * x * x * y * y * y * y * z * z - std::sqrt(56847.65625) * x * x * x * y * y * z * z * z * z - std::sqrt(142.119140625) * x * y * y * y * y * y * y * y * y + std::sqrt(17196.416015625) * x * y * y * y * y * y * y * z * z - std::sqrt(127907.2265625) * x * y * y * y * y * z * z * z * z) + e_2 * (-std::sqrt(3552.978515625) * x * x * x * x * x * x * x - std::sqrt(287791.259765625) * x * x * x * x * x * y * y - std::sqrt(127907.2265625) * x * x * x * x * x * z * z + std::sqrt(2220611.572265625) * x * x * x * y * y * y * y + std::sqrt(511628.90625) * x * x * x * y * y * z * z + std::sqrt(227390.625) * x * x * x * z * z * z * z + std::sqrt(3552.978515625) * x * y * y * y * y * y * y + std::sqrt(1151165.0390625) * x * y * y * y * y * z * z - std::sqrt(2046515.625) * x * y * y * z * z * z * z) + e_3 * (-std::sqrt(511628.90625) * x * x * x * x * x + std::sqrt(2046515.625) * x * x * x * y * y + std::sqrt(4604660.15625) * x * y * y * y * y) + e_4 * (-std::sqrt(2046515.625) * x * x * x + std::sqrt(18418640.625) * x * y * y);
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

        pc_15[k] = e_0 * (std::sqrt(31.7230224609375) * x * x * x * x * x * x * x * x * y * y * z - std::sqrt(4568.115234375) * x * x * x * x * x * x * y * y * z * z * z - std::sqrt(248.70849609375) * x * x * x * x * y * y * y * y * y * y * z + std::sqrt(4568.115234375) * x * x * x * x * y * y * y * y * z * z * z + std::sqrt(2030.2734375) * x * x * x * x * y * y * z * z * z * z * z - std::sqrt(81.2109375) * x * x * y * y * y * y * y * y * y * y * z + std::sqrt(14800.693359375) * x * x * y * y * y * y * y * y * z * z * z - std::sqrt(8121.09375) * x * x * y * y * y * y * z * z * z * z * z + std::sqrt(1.2689208984375) * y * y * y * y * y * y * y * y * y * y * z - std::sqrt(182.724609375) * y * y * y * y * y * y * y * y * z * z * z + std::sqrt(81.2109375) * y * y * y * y * y * y * z * z * z * z * z) + e_1 * (std::sqrt(31.7230224609375) * x * x * x * x * x * x * x * x * z - std::sqrt(2030.2734375) * x * x * x * x * x * x * y * y * z - std::sqrt(4568.115234375) * x * x * x * x * x * x * z * z * z - std::sqrt(1142.02880859375) * x * x * x * x * y * y * y * y * z - std::sqrt(24870.849609375) * x * x * x * x * y * y * z * z * z + std::sqrt(2030.2734375) * x * x * x * x * z * z * z * z * z + std::sqrt(324.84375) * x * x * y * y * y * y * y * y * z + std::sqrt(1766845.458984375) * x * x * y * y * y * y * z * z * z - std::sqrt(73089.84375) * x * x * y * y * z * z * z * z * z + std::sqrt(1.2689208984375) * y * y * y * y * y * y * y * y * z - std::sqrt(27794.443359375) * y * y * y * y * y * y * z * z * z + std::sqrt(2030.2734375) * y * y * y * y * z * z * z * z * z) + e_2 * (-std::sqrt(2030.2734375) * x * x * x * x * x * x * z - std::sqrt(456811.5234375) * x * x * x * x * y * y * z - std::sqrt(129937.5) * x * x * x * x * z * z * z + std::sqrt(4111303.7109375) * x * x * y * y * y * y * z + std::sqrt(4677750.0) * x * x * y * y * z * z * z - std::sqrt(50756.8359375) * y * y * y * y * y * y * z - std::sqrt(129937.5) * y * y * y * y * z * z * z) + e_3 * (-std::sqrt(657808.59375) * x * x * x * x * z + std::sqrt(23681109.375) * x * x * y * y * z - std::sqrt(657808.59375) * y * y * y * y * z);
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

        pc_16[k] = e_0 * (std::sqrt(475.8453369140625) * x * x * x * x * x * x * x * x * y * z * z - std::sqrt(3383.7890625) * x * x * x * x * x * x * y * z * z * z * z - std::sqrt(3730.62744140625) * x * x * x * x * y * y * y * y * y * z * z + std::sqrt(3383.7890625) * x * x * x * x * y * y * y * z * z * z * z + std::sqrt(135.3515625) * x * x * x * x * y * z * z * z * z * z * z - std::sqrt(1218.1640625) * x * x * y * y * y * y * y * y * y * z * z + std::sqrt(10963.4765625) * x * x * y * y * y * y * y * z * z * z * z - std::sqrt(541.40625) * x * x * y * y * y * z * z * z * z * z * z + std::sqrt(19.0338134765625) * y * y * y * y * y * y * y * y * y * z * z - std::sqrt(135.3515625) * y * y * y * y * y * y * y * z * z * z * z + std::sqrt(5.4140625) * y * y * y * y * y * z * z * z * z * z * z) + e_1 * (std::sqrt(475.8453369140625) * x * x * x * x * x * x * x * x * y + std::sqrt(68521.728515625) * x * x * x * x * x * x * y * z * z - std::sqrt(3730.62744140625) * x * x * x * x * y * y * y * y * y - std::sqrt(68521.728515625) * x * x * x * x * y * y * y * z * z - std::sqrt(274086.9140625) * x * x * x * x * y * z * z * z * z - std::sqrt(1218.1640625) * x * x * y * y * y * y * y * y * y - std::sqrt(222010.400390625) * x * x * y * y * y * y * y * z * z + std::sqrt(1096347.65625) * x * x * y * y * y * z * z * z * z + std::sqrt(19.0338134765625) * y * y * y * y * y * y * y * y * y + std::sqrt(2740.869140625) * y * y * y * y * y * y * y * z * z - std::sqrt(10963.4765625) * y * y * y * y * y * z * z * z * z) + e_2 * (std::sqrt(190338.134765625) * x * x * x * x * x * x * y - std::sqrt(190338.134765625) * x * x * x * x * y * y * y - std::sqrt(616695.556640625) * x * x * y * y * y * y * y + std::sqrt(7613.525390625) * y * y * y * y * y * y * y) + e_3 * (std::sqrt(3045410.15625) * x * x * x * x * y - std::sqrt(12181640.625) * x * x * y * y * y + std::sqrt(121816.40625) * y * y * y * y * y);
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

        pc_17[k] = e_0 * (std::sqrt(31.7230224609375) * x * x * x * x * x * x * x * x * x * y * z - std::sqrt(4568.115234375) * x * x * x * x * x * x * x * y * z * z * z - std::sqrt(248.70849609375) * x * x * x * x * x * y * y * y * y * y * z + std::sqrt(4568.115234375) * x * x * x * x * x * y * y * y * z * z * z + std::sqrt(2030.2734375) * x * x * x * x * x * y * z * z * z * z * z - std::sqrt(81.2109375) * x * x * x * y * y * y * y * y * y * y * z + std::sqrt(14800.693359375) * x * x * x * y * y * y * y * y * z * z * z - std::sqrt(8121.09375) * x * x * x * y * y * y * z * z * z * z * z + std::sqrt(1.2689208984375) * x * y * y * y * y * y * y * y * y * y * z - std::sqrt(182.724609375) * x * y * y * y * y * y * y * y * z * z * z + std::sqrt(81.2109375) * x * y * y * y * y * y * z * z * z * z * z) + e_1 * (std::sqrt(2030.2734375) * x * x * x * x * x * y * y * y * z - std::sqrt(586749.0234375) * x * x * x * x * x * y * z * z * z + std::sqrt(324.84375) * x * x * x * y * y * y * y * y * z + std::sqrt(982652.34375) * x * x * x * y * y * y * z * z * z + std::sqrt(32484.375) * x * x * x * y * z * z * z * z * z - std::sqrt(730.8984375) * x * y * y * y * y * y * y * y * z + std::sqrt(29317.1484375) * x * y * y * y * y * y * z * z * z - std::sqrt(32484.375) * x * y * y * y * z * z * z * z * z) + e_2 * (-std::sqrt(1169437.5) * x * x * x * x * x * y * z + std::sqrt(3248437.5) * x * x * x * y * y * y * z - std::sqrt(2079000.0) * x * x * x * y * z * z * z + std::sqrt(2079000.0) * x * y * y * y * z * z * z) + e_3 * (-std::sqrt(10524937.5) * x * x * x * y * z + std::sqrt(10524937.5) * x * y * y * y * z);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        pc_18[k] = e_0 * (-std::sqrt(888.24462890625) * x * x * x * x * x * x * x * x * y * z * z + std::sqrt(3552.978515625) * x * x * x * x * x * x * y * y * y * z * z + std::sqrt(3552.978515625) * x * x * x * x * x * x * y * z * z * z * z + std::sqrt(568.4765625) * x * x * x * x * y * y * y * y * y * z * z - std::sqrt(31976.806640625) * x * x * x * x * y * y * y * z * z * z * z - std::sqrt(3552.978515625) * x * x * y * y * y * y * y * y * y * z * z + std::sqrt(17196.416015625) * x * x * y * y * y * y * y * z * z * z * z + std::sqrt(35.52978515625) * y * y * y * y * y * y * y * y * y * z * z - std::sqrt(142.119140625) * y * y * y * y * y * y * y * z * z * z * z) + e_1 * (-std::sqrt(888.24462890625) * x * x * x * x * x * x * x * x * y + std::sqrt(3552.978515625) * x * x * x * x * x * x * y * y * y - std::sqrt(88824.462890625) * x * x * x * x * x * x * y * z * z + std::sqrt(568.4765625) * x * x * x * x * y * y * y * y * y + std::sqrt(3552.978515625) * x * x * x * x * y * y * y * z * z + std::sqrt(127907.2265625) * x * x * x * x * y * z * z * z * z - std::sqrt(3552.978515625) * x * x * y * y * y * y * y * y * y - std::sqrt(103604.853515625) * x * x * y * y * y * y * y * z * z + std::sqrt(56847.65625) * x * x * y * y * y * z * z * z * z + std::sqrt(35.52978515625) * y * y * y * y * y * y * y * y * y + std::sqrt(6963.837890625) * y * y * y * y * y * y * y * z * z - std::sqrt(14211.9140625) * y * y * y * y * y * z * z * z * z) + e_2 * (-std::sqrt(227390.625) * x * x * x * x * x * x * y + std::sqrt(355297.8515625) * x * x * x * x * y * y * y - std::sqrt(1151165.0390625) * x * x * x * x * y * z * z - std::sqrt(511628.90625) * x * x * y * y * y * y * y - std::sqrt(511628.90625) * x * x * y * y * y * z * z + std::sqrt(2046515.625) * x * x * y * z * z * z * z + std::sqrt(14211.9140625) * y * y * y * y * y * y * y + std::sqrt(127907.2265625) * y * y * y * y * y * z * z - std::sqrt(227390.625) * y * y * y * z * z * z * z) + e_3 * (-std::sqrt(4604660.15625) * x * x * x * x * y - std::sqrt(2046515.625) * x * x * y * y * y + std::sqrt(511628.90625) * y * y * y * y * y) + e_4 * (-std::sqrt(18418640.625) * x * x * y + std::sqrt(2046515.625) * y * y * y);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        pc_19[k] = e_0 * (-std::sqrt(37.01019287109375) * x * x * x * x * x * x * x * x * x * y * z + std::sqrt(592.1630859375) * x * x * x * x * x * x * x * y * y * y * z + std::sqrt(2368.65234375) * x * x * x * x * x * x * x * y * z * z * z - std::sqrt(53.294677734375) * x * x * x * x * x * y * y * y * y * y * z - std::sqrt(59216.30859375) * x * x * x * x * x * y * y * y * z * z * z - std::sqrt(1160.6396484375) * x * x * x * y * y * y * y * y * y * y * z + std::sqrt(91050.99609375) * x * x * x * y * y * y * y * y * z * z * z + std::sqrt(13.32366943359375) * x * y * y * y * y * y * y * y * y * y * z - std::sqrt(852.71484375) * x * y * y * y * y * y * y * y * z * z * z) + e_1 * (-std::sqrt(85271.484375) * x * x * x * x * x * y * y * y * z + std::sqrt(85271.484375) * x * x * x * x * x * y * z * z * z + std::sqrt(13643.4375) * x * x * x * y * y * y * y * y * z + std::sqrt(341085.9375) * x * x * x * y * y * y * z * z * z - std::sqrt(3410.859375) * x * y * y * y * y * y * y * y * z + std::sqrt(85271.484375) * x * y * y * y * y * y * z * z * z) + e_2 * (std::sqrt(5457375.0) * x * x * x * y * z * z * z + std::sqrt(5457375.0) * x * y * y * y * z * z * z) + e_3 * (std::sqrt(5457375.0) * x * x * x * y * z + std::sqrt(5457375.0) * x * y * y * y * z + std::sqrt(21829500.0) * x * y * z * z * z) + e_4 * (std::sqrt(49116375.0) * x * y * z);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        pc_20[k] = e_0 * (std::sqrt(666.1834716796875) * x * x * x * x * x * x * x * x * y * z * z - std::sqrt(42635.7421875) * x * x * x * x * x * x * y * y * y * z * z + std::sqrt(116075.80810546875) * x * x * x * x * y * y * y * y * y * z * z - std::sqrt(6821.71875) * x * x * y * y * y * y * y * y * y * z * z + std::sqrt(26.6473388671875) * y * y * y * y * y * y * y * y * y * z * z) + e_1 * (std::sqrt(666.1834716796875) * x * x * x * x * x * x * x * x * y - std::sqrt(42635.7421875) * x * x * x * x * x * x * y * y * y + std::sqrt(10658.935546875) * x * x * x * x * x * x * y * z * z + std::sqrt(116075.80810546875) * x * x * x * x * y * y * y * y * y + std::sqrt(95930.419921875) * x * x * x * x * y * y * y * z * z - std::sqrt(6821.71875) * x * x * y * y * y * y * y * y * y + std::sqrt(95930.419921875) * x * x * y * y * y * y * y * z * z + std::sqrt(26.6473388671875) * y * y * y * y * y * y * y * y * y + std::sqrt(10658.935546875) * y * y * y * y * y * y * y * z * z) + e_2 * (std::sqrt(10658.935546875) * x * x * x * x * x * x * y + std::sqrt(95930.419921875) * x * x * x * x * y * y * y + std::sqrt(1534886.71875) * x * x * x * x * y * z * z + std::sqrt(95930.419921875) * x * x * y * y * y * y * y + std::sqrt(6139546.875) * x * x * y * y * y * z * z + std::sqrt(10658.935546875) * y * y * y * y * y * y * y + std::sqrt(1534886.71875) * y * y * y * y * y * z * z) + e_3 * (std::sqrt(1534886.71875) * x * x * x * x * y + std::sqrt(6139546.875) * x * x * y * y * y + std::sqrt(24558187.5) * x * x * y * z * z + std::sqrt(1534886.71875) * y * y * y * y * y + std::sqrt(24558187.5) * y * y * y * z * z) + e_4 * (std::sqrt(24558187.5) * x * x * y + std::sqrt(24558187.5) * y * y * y + std::sqrt(24558187.5) * y * z * z) + e_5 * (std::sqrt(24558187.5) * y);

        pc_21[k] = e_0 * (std::sqrt(66.61834716796875) * x * x * x * x * x * x * x * x * x * y * z - std::sqrt(9593.0419921875) * x * x * x * x * x * x * x * y * y * y * z + std::sqrt(42305.315185546875) * x * x * x * x * x * y * y * y * y * y * z - std::sqrt(9593.0419921875) * x * x * x * y * y * y * y * y * y * y * z + std::sqrt(66.61834716796875) * x * y * y * y * y * y * y * y * y * y * z);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        pc_22[k] = e_0 * (-std::sqrt(48.44970703125) * x * x * x * x * x * x * x * x * x * y * y + std::sqrt(193.798828125) * x * x * x * x * x * x * x * y * y * y * y + std::sqrt(4844.970703125) * x * x * x * x * x * x * x * y * y * z * z + std::sqrt(31.0078125) * x * x * x * x * x * y * y * y * y * y * y - std::sqrt(43604.736328125) * x * x * x * x * x * y * y * y * y * z * z - std::sqrt(193.798828125) * x * x * x * y * y * y * y * y * y * y * y + std::sqrt(23449.658203125) * x * x * x * y * y * y * y * y * y * z * z + std::sqrt(1.93798828125) * x * y * y * y * y * y * y * y * y * y * y - std::sqrt(193.798828125) * x * y * y * y * y * y * y * y * y * z * z) + e_1 * (-std::sqrt(48.44970703125) * x * x * x * x * x * x * x * x * x - std::sqrt(9496.142578125) * x * x * x * x * x * x * x * y * y + std::sqrt(4844.970703125) * x * x * x * x * x * x * x * z * z + std::sqrt(27907.03125) * x * x * x * x * x * y * y * y * y + std::sqrt(43604.736328125) * x * x * x * x * x * y * y * z * z - std::sqrt(32752.001953125) * x * x * x * y * y * y * y * y * y + std::sqrt(43604.736328125) * x * x * x * y * y * y * y * z * z + std::sqrt(48.44970703125) * x * y * y * y * y * y * y * y * y + std::sqrt(4844.970703125) * x * y * y * y * y * y * y * z * z) + e_2 * (-std::sqrt(19379.8828125) * x * x * x * x * x * x * x - std::sqrt(174418.9453125) * x * x * x * x * x * y * y + std::sqrt(697675.78125) * x * x * x * x * x * z * z - std::sqrt(174418.9453125) * x * x * x * y * y * y * y + std::sqrt(2790703.125) * x * x * x * y * y * z * z - std::sqrt(19379.8828125) * x * y * y * y * y * y * y + std::sqrt(697675.78125) * x * y * y * y * y * z * z) + e_3 * (-std::sqrt(697675.78125) * x * x * x * x * x - std::sqrt(2790703.125) * x * x * x * y * y + std::sqrt(11162812.5) * x * x * x * z * z - std::sqrt(697675.78125) * x * y * y * y * y + std::sqrt(11162812.5) * x * y * y * z * z) + e_4 * (-std::sqrt(2790703.125) * x * x * x - std::sqrt(2790703.125) * x * y * y + std::sqrt(11162812.5) * x * z * z) + e_5 * (-std::sqrt(446512.5) * x);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        pc_23[k] = e_0 * (-std::sqrt(310.078125) * x * x * x * x * x * x * x * x * y * y * z + std::sqrt(310.078125) * x * x * x * x * x * x * y * y * y * y * z + std::sqrt(31007.8125) * x * x * x * x * x * x * y * y * z * z * z + std::sqrt(310.078125) * x * x * x * x * y * y * y * y * y * y * z - std::sqrt(124031.25) * x * x * x * x * y * y * y * y * z * z * z - std::sqrt(310.078125) * x * x * y * y * y * y * y * y * y * y * z + std::sqrt(31007.8125) * x * x * y * y * y * y * y * y * z * z * z) + e_1 * (-std::sqrt(310.078125) * x * x * x * x * x * x * x * x * z + std::sqrt(19845.0) * x * x * x * x * x * x * y * y * z + std::sqrt(31007.8125) * x * x * x * x * x * x * z * z * z - std::sqrt(279070.3125) * x * x * x * x * y * y * y * y * z + std::sqrt(279070.3125) * x * x * x * x * y * y * z * z * z + std::sqrt(19845.0) * x * x * y * y * y * y * y * y * z + std::sqrt(279070.3125) * x * x * y * y * y * y * z * z * z - std::sqrt(310.078125) * y * y * y * y * y * y * y * y * z + std::sqrt(31007.8125) * y * y * y * y * y * y * z * z * z) + e_2 * (std::sqrt(7751.953125) * x * x * x * x * x * x * z + std::sqrt(69767.578125) * x * x * x * x * y * y * z + std::sqrt(2511632.8125) * x * x * x * x * z * z * z + std::sqrt(69767.578125) * x * x * y * y * y * y * z + std::sqrt(10046531.25) * x * x * y * y * z * z * z + std::sqrt(7751.953125) * y * y * y * y * y * y * z + std::sqrt(2511632.8125) * y * y * y * y * z * z * z) + e_3 * (std::sqrt(4465125.0) * x * x * x * x * z + std::sqrt(17860500.0) * x * x * y * y * z + std::sqrt(17860500.0) * x * x * z * z * z + std::sqrt(4465125.0) * y * y * y * y * z + std::sqrt(17860500.0) * y * y * z * z * z) + e_4 * (std::sqrt(54697781.25) * x * x * z + std::sqrt(54697781.25) * y * y * z + std::sqrt(4465125.0) * z * z * z) + e_5 * (std::sqrt(17860500.0) * z);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        pc_24[k] = e_0 * (std::sqrt(9.68994140625) * x * x * x * x * x * x * x * x * x * y * y + std::sqrt(4.306640625) * x * x * x * x * x * x * x * y * y * y * y - std::sqrt(3139.541015625) * x * x * x * x * x * x * x * y * y * z * z - std::sqrt(17.2265625) * x * x * x * x * x * y * y * y * y * y * y + std::sqrt(348.837890625) * x * x * x * x * x * y * y * y * y * z * z + std::sqrt(62015.625) * x * x * x * x * x * y * y * z * z * z * z - std::sqrt(4.306640625) * x * x * x * y * y * y * y * y * y * y * y + std::sqrt(3139.541015625) * x * x * x * y * y * y * y * y * y * z * z - std::sqrt(110250.0) * x * x * x * y * y * y * y * z * z * z * z + std::sqrt(1.07666015625) * x * y * y * y * y * y * y * y * y * y * y - std::sqrt(348.837890625) * x * y * y * y * y * y * y * y * y * z * z + std::sqrt(6890.625) * x * y * y * y * y * y * y * z * z * z * z) + e_1 * (std::sqrt(9.68994140625) * x * x * x * x * x * x * x * x * x + std::sqrt(4689.931640625) * x * x * x * x * x * x * x * y * y - std::sqrt(3139.541015625) * x * x * x * x * x * x * x * z * z + std::sqrt(184535.244140625) * x * x * x * x * x * y * y * z * z + std::sqrt(62015.625) * x * x * x * x * x * z * z * z * z - std::sqrt(1899.228515625) * x * x * x * y * y * y * y * y * y - std::sqrt(931203.369140625) * x * x * x * y * y * y * y * z * z + std::sqrt(248062.5) * x * x * x * y * y * z * z * z * z + std::sqrt(474.80712890625) * x * y * y * y * y * y * y * y * y + std::sqrt(20503.916015625) * x * y * y * y * y * y * y * z * z + std::sqrt(62015.625) * x * y * y * y * y * z * z * z * z) + e_2 * (std::sqrt(3875.9765625) * x * x * x * x * x * x * x + std::sqrt(872094.7265625) * x * x * x * x * x * y * y + std::sqrt(139535.15625) * x * x * x * x * x * z * z - std::sqrt(655040.0390625) * x * x * x * y * y * y * y + std::sqrt(558140.625) * x * x * x * y * y * z * z + std::sqrt(2232562.5) * x * x * x * z * z * z * z + std::sqrt(96899.4140625) * x * y * y * y * y * y * y + std::sqrt(139535.15625) * x * y * y * y * y * z * z + std::sqrt(2232562.5) * x * y * y * z * z * z * z) + e_3 * (std::sqrt(759691.40625) * x * x * x * x * x + std::sqrt(3038765.625) * x * x * x * y * y + std::sqrt(20093062.5) * x * x * x * z * z + std::sqrt(759691.40625) * x * y * y * y * y + std::sqrt(20093062.5) * x * y * y * z * z + std::sqrt(3969000.0) * x * z * z * z * z) + e_4 * (std::sqrt(13953515.625) * x * x * x + std::sqrt(13953515.625) * x * y * y + std::sqrt(55814062.5) * x * z * z) + e_5 * (std::sqrt(20093062.5) * x);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        pc_25[k] = e_0 * (std::sqrt(103.359375) * x * x * x * x * x * x * x * x * y * y * z + std::sqrt(103.359375) * x * x * x * x * x * x * y * y * y * y * z - std::sqrt(14883.75) * x * x * x * x * x * x * y * y * z * z * z - std::sqrt(103.359375) * x * x * x * x * y * y * y * y * y * y * z + std::sqrt(41343.75) * x * x * x * x * y * y * z * z * z * z * z - std::sqrt(103.359375) * x * x * y * y * y * y * y * y * y * y * z + std::sqrt(14883.75) * x * x * y * y * y * y * y * y * z * z * z - std::sqrt(41343.75) * x * x * y * y * y * y * z * z * z * z * z) + e_1 * (std::sqrt(103.359375) * x * x * x * x * x * x * x * x * z - std::sqrt(413.4375) * x * x * x * x * x * x * y * y * z - std::sqrt(14883.75) * x * x * x * x * x * x * z * z * z + std::sqrt(41343.75) * x * x * x * x * y * y * z * z * z + std::sqrt(41343.75) * x * x * x * x * z * z * z * z * z + std::sqrt(413.4375) * x * x * y * y * y * y * y * y * z - std::sqrt(41343.75) * x * x * y * y * y * y * z * z * z - std::sqrt(103.359375) * y * y * y * y * y * y * y * y * z + std::sqrt(14883.75) * y * y * y * y * y * y * z * z * z - std::sqrt(41343.75) * y * y * y * y * z * z * z * z * z) + e_2 * (-std::sqrt(2583.984375) * x * x * x * x * x * x * z + std::sqrt(23255.859375) * x * x * x * x * y * y * z + std::sqrt(41343.75) * x * x * x * x * z * z * z - std::sqrt(23255.859375) * x * x * y * y * y * y * z + std::sqrt(372093.75) * x * x * z * z * z * z * z + std::sqrt(2583.984375) * y * y * y * y * y * y * z - std::sqrt(41343.75) * y * y * y * y * z * z * z - std::sqrt(372093.75) * y * y * z * z * z * z * z) + e_3 * (std::sqrt(5953500.0) * x * x * z * z * z - std::sqrt(5953500.0) * y * y * z * z * z) + e_4 * (std::sqrt(3348843.75) * x * x * z - std::sqrt(3348843.75) * y * y * z);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        pc_26[k] = e_0 * (-std::sqrt(0.9228515625) * x * x * x * x * x * x * x * x * x * y * y - std::sqrt(3.69140625) * x * x * x * x * x * x * x * y * y * y * y + std::sqrt(446.66015625) * x * x * x * x * x * x * x * y * y * z * z + std::sqrt(446.66015625) * x * x * x * x * x * y * y * y * y * z * z - std::sqrt(15120.0) * x * x * x * x * x * y * y * z * z * z * z + std::sqrt(3.69140625) * x * x * x * y * y * y * y * y * y * y * y - std::sqrt(446.66015625) * x * x * x * y * y * y * y * y * y * z * z + std::sqrt(5906.25) * x * x * x * y * y * z * z * z * z * z * z + std::sqrt(0.9228515625) * x * y * y * y * y * y * y * y * y * y * y - std::sqrt(446.66015625) * x * y * y * y * y * y * y * y * y * z * z + std::sqrt(15120.0) * x * y * y * y * y * y * y * z * z * z * z - std::sqrt(5906.25) * x * y * y * y * y * z * z * z * z * z * z) + e_1 * (-std::sqrt(0.9228515625) * x * x * x * x * x * x * x * x * x - std::sqrt(623.84765625) * x * x * x * x * x * x * x * y * y + std::sqrt(446.66015625) * x * x * x * x * x * x * x * z * z - std::sqrt(369.140625) * x * x * x * x * x * y * y * y * y - std::sqrt(27940.25390625) * x * x * x * x * x * y * y * z * z - std::sqrt(15120.0) * x * x * x * x * x * z * z * z * z + std::sqrt(1066.81640625) * x * x * x * y * y * y * y * y * y - std::sqrt(11166.50390625) * x * x * x * y * y * y * y * z * z - std::sqrt(5906.25) * x * x * x * y * y * z * z * z * z + std::sqrt(5906.25) * x * x * x * z * z * z * z * z * z + std::sqrt(776.1181640625) * x * y * y * y * y * y * y * y * y + std::sqrt(6825.41015625) * x * y * y * y * y * y * y * z * z + std::sqrt(478406.25) * x * y * y * y * y * z * z * z * z - std::sqrt(53156.25) * x * y * y * z * z * z * z * z * z) + e_2 * (-std::sqrt(369.140625) * x * x * x * x * x * x * x - std::sqrt(162791.015625) * x * x * x * x * x * y * y - std::sqrt(53156.25) * x * x * x * x * x * z * z + std::sqrt(9228.515625) * x * x * x * y * y * y * y - std::sqrt(1913625.0) * x * x * x * y * y * z * z - std::sqrt(5906.25) * x * x * x * z * z * z * z + std::sqrt(230712.890625) * x * y * y * y * y * y * y + std::sqrt(6431906.25) * x * y * y * y * y * z * z + std::sqrt(53156.25) * x * y * y * z * z * z * z) + e_3 * (-std::sqrt(119601.5625) * x * x * x * x * x - std::sqrt(2604656.25) * x * x * x * y * y - std::sqrt(1913625.0) * x * x * x * z * z + std::sqrt(11176101.5625) * x * y * y * y * y + std::sqrt(17222625.0) * x * y * y * z * z) + e_4 * (-std::sqrt(2604656.25) * x * x * x + std::sqrt(23441906.25) * x * y * y);
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

        pc_27[k] = e_0 * (-std::sqrt(13.8427734375) * x * x * x * x * x * x * x * x * x * y * z - std::sqrt(55.37109375) * x * x * x * x * x * x * x * y * y * y * z + std::sqrt(2220.99609375) * x * x * x * x * x * x * x * y * z * z * z + std::sqrt(2220.99609375) * x * x * x * x * x * y * y * y * z * z * z - std::sqrt(10241.4375) * x * x * x * x * x * y * z * z * z * z * z + std::sqrt(55.37109375) * x * x * x * y * y * y * y * y * y * y * z - std::sqrt(2220.99609375) * x * x * x * y * y * y * y * y * z * z * z + std::sqrt(393.75) * x * x * x * y * z * z * z * z * z * z * z + std::sqrt(13.8427734375) * x * y * y * y * y * y * y * y * y * y * z - std::sqrt(2220.99609375) * x * y * y * y * y * y * y * y * z * z * z + std::sqrt(10241.4375) * x * y * y * y * y * y * z * z * z * z * z - std::sqrt(393.75) * x * y * y * y * z * z * z * z * z * z * z) + e_1 * (-std::sqrt(221.484375) * x * x * x * x * x * x * x * y * z - std::sqrt(221.484375) * x * x * x * x * x * y * y * y * z + std::sqrt(14175.0) * x * x * x * x * x * y * z * z * z + std::sqrt(221.484375) * x * x * x * y * y * y * y * y * z - std::sqrt(354375.0) * x * x * x * y * z * z * z * z * z + std::sqrt(221.484375) * x * y * y * y * y * y * y * y * z - std::sqrt(14175.0) * x * y * y * y * y * y * z * z * z + std::sqrt(354375.0) * x * y * y * y * z * z * z * z * z) + e_2 * (-std::sqrt(5670000.0) * x * x * x * y * z * z * z + std::sqrt(5670000.0) * x * y * y * y * z * z * z) + e_3 * (-std::sqrt(5670000.0) * x * x * x * y * z + std::sqrt(5670000.0) * x * y * y * y * z);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        pc_28[k] = e_0 * (-std::sqrt(0.9228515625) * x * x * x * x * x * x * x * x * x * x * y - std::sqrt(3.69140625) * x * x * x * x * x * x * x * x * y * y * y + std::sqrt(446.66015625) * x * x * x * x * x * x * x * x * y * z * z + std::sqrt(446.66015625) * x * x * x * x * x * x * y * y * y * z * z - std::sqrt(15120.0) * x * x * x * x * x * x * y * z * z * z * z + std::sqrt(3.69140625) * x * x * x * x * y * y * y * y * y * y * y - std::sqrt(446.66015625) * x * x * x * x * y * y * y * y * y * z * z + std::sqrt(5906.25) * x * x * x * x * y * z * z * z * z * z * z + std::sqrt(0.9228515625) * x * x * y * y * y * y * y * y * y * y * y - std::sqrt(446.66015625) * x * x * y * y * y * y * y * y * y * z * z + std::sqrt(15120.0) * x * x * y * y * y * y * y * z * z * z * z - std::sqrt(5906.25) * x * x * y * y * y * z * z * z * z * z * z) + e_1 * (-std::sqrt(776.1181640625) * x * x * x * x * x * x * x * x * y - std::sqrt(1066.81640625) * x * x * x * x * x * x * y * y * y - std::sqrt(6825.41015625) * x * x * x * x * x * x * y * z * z + std::sqrt(369.140625) * x * x * x * x * y * y * y * y * y + std::sqrt(11166.50390625) * x * x * x * x * y * y * y * z * z - std::sqrt(478406.25) * x * x * x * x * y * z * z * z * z + std::sqrt(623.84765625) * x * x * y * y * y * y * y * y * y + std::sqrt(27940.25390625) * x * x * y * y * y * y * y * z * z + std::sqrt(5906.25) * x * x * y * y * y * z * z * z * z + std::sqrt(53156.25) * x * x * y * z * z * z * z * z * z + std::sqrt(0.9228515625) * y * y * y * y * y * y * y * y * y - std::sqrt(446.66015625) * y * y * y * y * y * y * y * z * z + std::sqrt(15120.0) * y * y * y * y * y * z * z * z * z - std::sqrt(5906.25) * y * y * y * z * z * z * z * z * z) + e_2 * (-std::sqrt(230712.890625) * x * x * x * x * x * x * y - std::sqrt(9228.515625) * x * x * x * x * y * y * y - std::sqrt(6431906.25) * x * x * x * x * y * z * z + std::sqrt(162791.015625) * x * x * y * y * y * y * y + std::sqrt(1913625.0) * x * x * y * y * y * z * z - std::sqrt(53156.25) * x * x * y * z * z * z * z + std::sqrt(369.140625) * y * y * y * y * y * y * y + std::sqrt(53156.25) * y * y * y * y * y * z * z + std::sqrt(5906.25) * y * y * y * z * z * z * z) + e_3 * (-std::sqrt(11176101.5625) * x * x * x * x * y + std::sqrt(2604656.25) * x * x * y * y * y - std::sqrt(17222625.0) * x * x * y * z * z + std::sqrt(119601.5625) * y * y * y * y * y + std::sqrt(1913625.0) * y * y * y * z * z) + e_4 * (-std::sqrt(23441906.25) * x * x * y + std::sqrt(2604656.25) * y * y * y);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        pc_29[k] = e_0 * (std::sqrt(25.83984375) * x * x * x * x * x * x * x * x * x * y * z - std::sqrt(3720.9375) * x * x * x * x * x * x * x * y * z * z * z - std::sqrt(103.359375) * x * x * x * x * x * y * y * y * y * y * z + std::sqrt(3720.9375) * x * x * x * x * x * y * y * y * z * z * z + std::sqrt(10335.9375) * x * x * x * x * x * y * z * z * z * z * z + std::sqrt(3720.9375) * x * x * x * y * y * y * y * y * z * z * z - std::sqrt(41343.75) * x * x * x * y * y * y * z * z * z * z * z + std::sqrt(25.83984375) * x * y * y * y * y * y * y * y * y * y * z - std::sqrt(3720.9375) * x * y * y * y * y * y * y * y * z * z * z + std::sqrt(10335.9375) * x * y * y * y * y * y * z * z * z * z * z) + e_1 * (std::sqrt(6615.0) * x * x * x * x * x * y * y * y * z - std::sqrt(6615.0) * x * x * x * x * x * y * z * z * z + std::sqrt(6615.0) * x * x * x * y * y * y * y * y * z - std::sqrt(661500.0) * x * x * x * y * y * y * z * z * z + std::sqrt(165375.0) * x * x * x * y * z * z * z * z * z - std::sqrt(6615.0) * x * y * y * y * y * y * z * z * z + std::sqrt(165375.0) * x * y * y * y * z * z * z * z * z) + e_2 * (-std::sqrt(165375.0) * x * x * x * y * y * y * z + std::sqrt(165375.0) * x * x * x * y * z * z * z + std::sqrt(165375.0) * x * y * y * y * z * z * z + std::sqrt(1488375.0) * x * y * z * z * z * z * z) + e_3 * (std::sqrt(23814000.0) * x * y * z * z * z) + e_4 * (std::sqrt(13395375.0) * x * y * z);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        pc_30[k] = e_0 * (std::sqrt(1.07666015625) * x * x * x * x * x * x * x * x * x * x * y - std::sqrt(4.306640625) * x * x * x * x * x * x * x * x * y * y * y - std::sqrt(348.837890625) * x * x * x * x * x * x * x * x * y * z * z - std::sqrt(17.2265625) * x * x * x * x * x * x * y * y * y * y * y + std::sqrt(3139.541015625) * x * x * x * x * x * x * y * y * y * z * z + std::sqrt(6890.625) * x * x * x * x * x * x * y * z * z * z * z + std::sqrt(4.306640625) * x * x * x * x * y * y * y * y * y * y * y + std::sqrt(348.837890625) * x * x * x * x * y * y * y * y * y * z * z - std::sqrt(110250.0) * x * x * x * x * y * y * y * z * z * z * z + std::sqrt(9.68994140625) * x * x * y * y * y * y * y * y * y * y * y - std::sqrt(3139.541015625) * x * x * y * y * y * y * y * y * y * z * z + std::sqrt(62015.625) * x * x * y * y * y * y * y * z * z * z * z) + e_1 * (std::sqrt(474.80712890625) * x * x * x * x * x * x * x * x * y - std::sqrt(1899.228515625) * x * x * x * x * x * x * y * y * y + std::sqrt(20503.916015625) * x * x * x * x * x * x * y * z * z - std::sqrt(931203.369140625) * x * x * x * x * y * y * y * z * z + std::sqrt(62015.625) * x * x * x * x * y * z * z * z * z + std::sqrt(4689.931640625) * x * x * y * y * y * y * y * y * y + std::sqrt(184535.244140625) * x * x * y * y * y * y * y * z * z + std::sqrt(248062.5) * x * x * y * y * y * z * z * z * z + std::sqrt(9.68994140625) * y * y * y * y * y * y * y * y * y - std::sqrt(3139.541015625) * y * y * y * y * y * y * y * z * z + std::sqrt(62015.625) * y * y * y * y * y * z * z * z * z) + e_2 * (std::sqrt(96899.4140625) * x * x * x * x * x * x * y - std::sqrt(655040.0390625) * x * x * x * x * y * y * y + std::sqrt(139535.15625) * x * x * x * x * y * z * z + std::sqrt(872094.7265625) * x * x * y * y * y * y * y + std::sqrt(558140.625) * x * x * y * y * y * z * z + std::sqrt(2232562.5) * x * x * y * z * z * z * z + std::sqrt(3875.9765625) * y * y * y * y * y * y * y + std::sqrt(139535.15625) * y * y * y * y * y * z * z + std::sqrt(2232562.5) * y * y * y * z * z * z * z) + e_3 * (std::sqrt(759691.40625) * x * x * x * x * y + std::sqrt(3038765.625) * x * x * y * y * y + std::sqrt(20093062.5) * x * x * y * z * z + std::sqrt(759691.40625) * y * y * y * y * y + std::sqrt(20093062.5) * y * y * y * z * z + std::sqrt(3969000.0) * y * z * z * z * z) + e_4 * (std::sqrt(13953515.625) * x * x * y + std::sqrt(13953515.625) * y * y * y + std::sqrt(55814062.5) * y * z * z) + e_5 * (std::sqrt(20093062.5) * y);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        pc_31[k] = e_0 * (-std::sqrt(19.3798828125) * x * x * x * x * x * x * x * x * x * y * z + std::sqrt(697.67578125) * x * x * x * x * x * x * x * y * y * y * z + std::sqrt(1937.98828125) * x * x * x * x * x * x * x * y * z * z * z - std::sqrt(94961.42578125) * x * x * x * x * x * y * y * y * z * z * z - std::sqrt(697.67578125) * x * x * x * y * y * y * y * y * y * y * z + std::sqrt(94961.42578125) * x * x * x * y * y * y * y * y * z * z * z + std::sqrt(19.3798828125) * x * y * y * y * y * y * y * y * y * y * z - std::sqrt(1937.98828125) * x * y * y * y * y * y * y * y * z * z * z) + e_1 * (std::sqrt(2790.703125) * x * x * x * x * x * x * x * y * z - std::sqrt(136744.453125) * x * x * x * x * x * y * y * y * z + std::sqrt(136744.453125) * x * x * x * y * y * y * y * y * z - std::sqrt(2790.703125) * x * y * y * y * y * y * y * y * z);

        pc_32[k] = e_0 * (-std::sqrt(1.93798828125) * x * x * x * x * x * x * x * x * x * x * y + std::sqrt(193.798828125) * x * x * x * x * x * x * x * x * y * y * y + std::sqrt(193.798828125) * x * x * x * x * x * x * x * x * y * z * z - std::sqrt(31.0078125) * x * x * x * x * x * x * y * y * y * y * y - std::sqrt(23449.658203125) * x * x * x * x * x * x * y * y * y * z * z - std::sqrt(193.798828125) * x * x * x * x * y * y * y * y * y * y * y + std::sqrt(43604.736328125) * x * x * x * x * y * y * y * y * y * z * z + std::sqrt(48.44970703125) * x * x * y * y * y * y * y * y * y * y * y - std::sqrt(4844.970703125) * x * x * y * y * y * y * y * y * y * z * z) + e_1 * (-std::sqrt(48.44970703125) * x * x * x * x * x * x * x * x * y + std::sqrt(32752.001953125) * x * x * x * x * x * x * y * y * y - std::sqrt(4844.970703125) * x * x * x * x * x * x * y * z * z - std::sqrt(27907.03125) * x * x * x * x * y * y * y * y * y - std::sqrt(43604.736328125) * x * x * x * x * y * y * y * z * z + std::sqrt(9496.142578125) * x * x * y * y * y * y * y * y * y - std::sqrt(43604.736328125) * x * x * y * y * y * y * y * z * z + std::sqrt(48.44970703125) * y * y * y * y * y * y * y * y * y - std::sqrt(4844.970703125) * y * y * y * y * y * y * y * z * z) + e_2 * (std::sqrt(19379.8828125) * x * x * x * x * x * x * y + std::sqrt(174418.9453125) * x * x * x * x * y * y * y - std::sqrt(697675.78125) * x * x * x * x * y * z * z + std::sqrt(174418.9453125) * x * x * y * y * y * y * y - std::sqrt(2790703.125) * x * x * y * y * y * z * z + std::sqrt(19379.8828125) * y * y * y * y * y * y * y - std::sqrt(697675.78125) * y * y * y * y * y * z * z) + e_3 * (std::sqrt(697675.78125) * x * x * x * x * y + std::sqrt(2790703.125) * x * x * y * y * y - std::sqrt(11162812.5) * x * x * y * z * z + std::sqrt(697675.78125) * y * y * y * y * y - std::sqrt(11162812.5) * y * y * y * z * z) + e_4 * (std::sqrt(2790703.125) * x * x * y + std::sqrt(2790703.125) * y * y * y - std::sqrt(11162812.5) * y * z * z) + e_5 * (std::sqrt(446512.5) * y);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        pc_33[k] = e_0 * (-std::sqrt(817.5888061523438) * x * x * x * x * x * x * x * x * y * y * z + std::sqrt(1453.4912109375) * x * x * x * x * x * x * y * y * y * y * z + std::sqrt(5813.96484375) * x * x * x * x * x * x * y * y * z * z * z + std::sqrt(1758.724365234375) * x * x * x * x * y * y * y * y * y * y * z - std::sqrt(31653.80859375) * x * x * x * x * y * y * y * y * z * z * z - std::sqrt(523.2568359375) * x * x * y * y * y * y * y * y * y * y * z + std::sqrt(4366.93359375) * x * x * y * y * y * y * y * y * z * z * z + std::sqrt(3.63372802734375) * y * y * y * y * y * y * y * y * y * y * z - std::sqrt(25.83984375) * y * y * y * y * y * y * y * y * z * z * z) + e_1 * (-std::sqrt(817.5888061523438) * x * x * x * x * x * x * x * x * z - std::sqrt(117732.7880859375) * x * x * x * x * x * x * y * y * z + std::sqrt(5813.96484375) * x * x * x * x * x * x * z * z * z + std::sqrt(445131.6833496094) * x * x * x * x * y * y * y * y * z + std::sqrt(5813.96484375) * x * x * x * x * y * y * z * z * z - std::sqrt(36337.2802734375) * x * x * y * y * y * y * y * y * z - std::sqrt(5813.96484375) * x * x * y * y * y * y * z * z * z + std::sqrt(2271.0800170898438) * y * y * y * y * y * y * y * y * z - std::sqrt(5813.96484375) * y * y * y * y * y * y * z * z * z) + e_2 * (-std::sqrt(209302.734375) * x * x * x * x * x * x * z - std::sqrt(209302.734375) * x * x * x * x * y * y * z + std::sqrt(372093.75) * x * x * x * x * z * z * z + std::sqrt(209302.734375) * x * x * y * y * y * y * z + std::sqrt(209302.734375) * y * y * y * y * y * y * z - std::sqrt(372093.75) * y * y * y * y * z * z * z) + e_3 * (-std::sqrt(3348843.75) * x * x * x * x * z + std::sqrt(1488375.0) * x * x * z * z * z + std::sqrt(3348843.75) * y * y * y * y * z - std::sqrt(1488375.0) * y * y * z * z * z) + e_4 * (-std::sqrt(3348843.75) * x * x * z + std::sqrt(3348843.75) * y * y * z);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        pc_34[k] = e_0 * (-std::sqrt(5232.568359375) * x * x * x * x * x * x * x * y * y * z * z + std::sqrt(581.396484375) * x * x * x * x * x * y * y * y * y * z * z + std::sqrt(37209.375) * x * x * x * x * x * y * y * z * z * z * z + std::sqrt(5232.568359375) * x * x * x * y * y * y * y * y * y * z * z - std::sqrt(66150.0) * x * x * x * y * y * y * y * z * z * z * z - std::sqrt(581.396484375) * x * y * y * y * y * y * y * y * y * z * z + std::sqrt(4134.375) * x * y * y * y * y * y * y * z * z * z * z) + e_1 * (-std::sqrt(5232.568359375) * x * x * x * x * x * x * x * y * y - std::sqrt(5232.568359375) * x * x * x * x * x * x * x * z * z + std::sqrt(581.396484375) * x * x * x * x * x * y * y * y * y - std::sqrt(47093.115234375) * x * x * x * x * x * y * y * z * z + std::sqrt(37209.375) * x * x * x * x * x * z * z * z * z + std::sqrt(5232.568359375) * x * x * x * y * y * y * y * y * y - std::sqrt(47093.115234375) * x * x * x * y * y * y * y * z * z + std::sqrt(148837.5) * x * x * x * y * y * z * z * z * z - std::sqrt(581.396484375) * x * y * y * y * y * y * y * y * y - std::sqrt(5232.568359375) * x * y * y * y * y * y * y * z * z + std::sqrt(37209.375) * x * y * y * y * y * z * z * z * z) + e_2 * (-std::sqrt(5232.568359375) * x * x * x * x * x * x * x - std::sqrt(633140.771484375) * x * x * x * x * x * y * y - std::sqrt(83721.09375) * x * x * x * x * x * z * z + std::sqrt(307558.740234375) * x * x * x * y * y * y * y - std::sqrt(334884.375) * x * x * x * y * y * z * z + std::sqrt(1339537.5) * x * x * x * z * z * z * z - std::sqrt(70348.974609375) * x * y * y * y * y * y * y - std::sqrt(83721.09375) * x * y * y * y * y * z * z + std::sqrt(1339537.5) * x * y * y * z * z * z * z) + e_3 * (-std::sqrt(753489.84375) * x * x * x * x * x - std::sqrt(3013959.375) * x * x * x * y * y + std::sqrt(1339537.5) * x * x * x * z * z - std::sqrt(753489.84375) * x * y * y * y * y + std::sqrt(1339537.5) * x * y * y * z * z + std::sqrt(2381400.0) * x * z * z * z * z) + e_4 * (-std::sqrt(5358150.0) * x * x * x - std::sqrt(5358150.0) * x * y * y + std::sqrt(12055837.5) * x * z * z) + e_5 * (-std::sqrt(1339537.5) * x);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        pc_35[k] = e_0 * (std::sqrt(163.51776123046875) * x * x * x * x * x * x * x * x * y * y * z + std::sqrt(290.6982421875) * x * x * x * x * x * x * y * y * y * y * z - std::sqrt(18604.6875) * x * x * x * x * x * x * y * y * z * z * z - std::sqrt(8.074951171875) * x * x * x * x * y * y * y * y * y * y * z - std::sqrt(2067.1875) * x * x * x * x * y * y * y * y * z * z * z + std::sqrt(74418.75) * x * x * x * x * y * y * z * z * z * z * z - std::sqrt(32.2998046875) * x * x * y * y * y * y * y * y * y * y * z + std::sqrt(5742.1875) * x * x * y * y * y * y * y * y * z * z * z - std::sqrt(33075.0) * x * x * y * y * y * y * z * z * z * z * z + std::sqrt(2.01873779296875) * y * y * y * y * y * y * y * y * y * y * z - std::sqrt(229.6875) * y * y * y * y * y * y * y * y * z * z * z + std::sqrt(918.75) * y * y * y * y * y * y * z * z * z * z * z) + e_1 * (std::sqrt(163.51776123046875) * x * x * x * x * x * x * x * x * z + std::sqrt(2616.2841796875) * x * x * x * x * x * x * y * y * z - std::sqrt(18604.6875) * x * x * x * x * x * x * z * z * z + std::sqrt(5886.639404296875) * x * x * x * x * y * y * y * y * z + std::sqrt(167442.1875) * x * x * x * x * y * y * z * z * z + std::sqrt(74418.75) * x * x * x * x * z * z * z * z * z + std::sqrt(2616.2841796875) * x * x * y * y * y * y * y * y * z - std::sqrt(911629.6875) * x * x * y * y * y * y * z * z * z + std::sqrt(297675.0) * x * x * y * y * z * z * z * z * z + std::sqrt(163.51776123046875) * y * y * y * y * y * y * y * y * z - std::sqrt(2067.1875) * y * y * y * y * y * y * z * z * z + std::sqrt(74418.75) * y * y * y * y * z * z * z * z * z) + e_2 * (std::sqrt(1506979.6875) * x * x * x * x * y * y * z + std::sqrt(297675.0) * x * x * x * x * z * z * z - std::sqrt(669768.75) * x * x * y * y * y * y * z + std::sqrt(1190700.0) * x * x * y * y * z * z * z + std::sqrt(1190700.0) * x * x * z * z * z * z * z + std::sqrt(18604.6875) * y * y * y * y * y * y * z + std::sqrt(297675.0) * y * y * y * y * z * z * z + std::sqrt(1190700.0) * y * y * z * z * z * z * z) + e_3 * (std::sqrt(911629.6875) * x * x * x * x * z + std::sqrt(3646518.75) * x * x * y * y * z + std::sqrt(25930800.0) * x * x * z * z * z + std::sqrt(911629.6875) * y * y * y * y * z + std::sqrt(25930800.0) * y * y * z * z * z + std::sqrt(529200.0) * z * z * z * z * z) + e_4 * (std::sqrt(32818668.75) * x * x * z + std::sqrt(32818668.75) * y * y * z + std::sqrt(19051200.0) * z * z * z) + e_5 * (std::sqrt(24111675.0) * z);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        pc_36[k] = e_0 * (std::sqrt(1744.189453125) * x * x * x * x * x * x * x * y * y * z * z + std::sqrt(4844.970703125) * x * x * x * x * x * y * y * y * y * z * z - std::sqrt(37984.5703125) * x * x * x * x * x * y * y * z * z * z * z + std::sqrt(193.798828125) * x * x * x * y * y * y * y * y * y * z * z - std::sqrt(16882.03125) * x * x * x * y * y * y * y * z * z * z * z + std::sqrt(49612.5) * x * x * x * y * y * z * z * z * z * z * z - std::sqrt(193.798828125) * x * y * y * y * y * y * y * y * y * z * z + std::sqrt(4220.5078125) * x * y * y * y * y * y * y * z * z * z * z - std::sqrt(5512.5) * x * y * y * y * y * z * z * z * z * z * z) + e_1 * (std::sqrt(1744.189453125) * x * x * x * x * x * x * x * y * y + std::sqrt(1744.189453125) * x * x * x * x * x * x * x * z * z + std::sqrt(4844.970703125) * x * x * x * x * x * y * y * y * y + std::sqrt(15697.705078125) * x * x * x * x * x * y * y * z * z - std::sqrt(37984.5703125) * x * x * x * x * x * z * z * z * z + std::sqrt(193.798828125) * x * x * x * y * y * y * y * y * y + std::sqrt(15697.705078125) * x * x * x * y * y * y * y * z * z + std::sqrt(375194.53125) * x * x * x * y * y * z * z * z * z + std::sqrt(49612.5) * x * x * x * z * z * z * z * z * z - std::sqrt(193.798828125) * x * y * y * y * y * y * y * y * y + std::sqrt(1744.189453125) * x * y * y * y * y * y * y * z * z - std::sqrt(279845.5078125) * x * y * y * y * y * z * z * z * z + std::sqrt(49612.5) * x * y * y * z * z * z * z * z * z) + e_2 * (std::sqrt(1744.189453125) * x * x * x * x * x * x * x + std::sqrt(504070.751953125) * x * x * x * x * x * y * y - std::sqrt(6976.7578125) * x * x * x * x * x * z * z + std::sqrt(265310.595703125) * x * x * x * y * y * y * y + std::sqrt(8065132.03125) * x * x * x * y * y * z * z + std::sqrt(1004653.125) * x * x * x * z * z * z * z - std::sqrt(23449.658203125) * x * y * y * y * y * y * y - std::sqrt(1179072.0703125) * x * y * y * y * y * z * z + std::sqrt(1004653.125) * x * y * y * z * z * z * z + std::sqrt(198450.0) * x * z * z * z * z * z * z) + e_3 * (std::sqrt(251163.28125) * x * x * x * x * x + std::sqrt(18865153.125) * x * x * x * y * y + std::sqrt(7144200.0) * x * x * x * z * z - std::sqrt(375194.53125) * x * y * y * y * y + std::sqrt(7144200.0) * x * y * y * z * z + std::sqrt(12700800.0) * x * z * z * z * z) + e_4 * (std::sqrt(9041878.125) * x * x * x + std::sqrt(9041878.125) * x * y * y + std::sqrt(64297800.0) * x * z * z) + e_5 * (std::sqrt(16074450.0) * x);
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

        pc_37[k] = e_0 * (-std::sqrt(15.5731201171875) * x * x * x * x * x * x * x * x * y * y * z - std::sqrt(110.7421875) * x * x * x * x * x * x * y * y * y * y * z + std::sqrt(3349.951171875) * x * x * x * x * x * x * y * y * z * z * z - std::sqrt(62.29248046875) * x * x * x * x * y * y * y * y * y * y * z + std::sqrt(9305.419921875) * x * x * x * x * y * y * y * y * z * z * z - std::sqrt(24916.9921875) * x * x * x * x * y * y * z * z * z * z * z + std::sqrt(372.216796875) * x * x * y * y * y * y * y * y * z * z * z - std::sqrt(11074.21875) * x * x * y * y * y * y * z * z * z * z * z + std::sqrt(7087.5) * x * x * y * y * z * z * z * z * z * z * z + std::sqrt(1.7303466796875) * y * y * y * y * y * y * y * y * y * y * z - std::sqrt(372.216796875) * y * y * y * y * y * y * y * y * z * z * z + std::sqrt(2768.5546875) * y * y * y * y * y * y * z * z * z * z * z - std::sqrt(787.5) * y * y * y * y * z * z * z * z * z * z * z) + e_1 * (-std::sqrt(15.5731201171875) * x * x * x * x * x * x * x * x * z + std::sqrt(3349.951171875) * x * x * x * x * x * x * z * z * z + std::sqrt(173.03466796875) * x * x * x * x * y * y * y * y * z - std::sqrt(17303.466796875) * x * x * x * x * y * y * z * z * z - std::sqrt(24916.9921875) * x * x * x * x * z * z * z * z * z + std::sqrt(110.7421875) * x * x * y * y * y * y * y * y * z - std::sqrt(33914.794921875) * x * x * y * y * y * y * z * z * z + std::sqrt(35880.46875) * x * x * y * y * z * z * z * z * z + std::sqrt(7087.5) * x * x * z * z * z * z * z * z * z + std::sqrt(1.7303466796875) * y * y * y * y * y * y * y * y * z + std::sqrt(27.685546875) * y * y * y * y * y * y * z * z * z + std::sqrt(8970.1171875) * y * y * y * y * z * z * z * z * z - std::sqrt(7087.5) * y * y * z * z * z * z * z * z * z) + e_2 * (std::sqrt(996.6796875) * x * x * x * x * x * x * z - std::sqrt(24916.9921875) * x * x * x * x * y * y * z - std::sqrt(177187.5) * x * x * x * x * z * z * z - std::sqrt(24916.9921875) * x * x * y * y * y * y * z + std::sqrt(255150.0) * x * x * z * z * z * z * z + std::sqrt(996.6796875) * y * y * y * y * y * y * z + std::sqrt(177187.5) * y * y * y * y * z * z * z - std::sqrt(255150.0) * y * y * z * z * z * z * z) + e_3 * (-std::sqrt(99667.96875) * x * x * x * x * z - std::sqrt(398671.875) * x * x * y * y * z + std::sqrt(708750.0) * x * x * z * z * z + std::sqrt(276855.46875) * y * y * y * y * z - std::sqrt(708750.0) * y * y * z * z * z);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        pc_38[k] = e_0 * (-std::sqrt(233.5968017578125) * x * x * x * x * x * x * x * x * y * z * z - std::sqrt(1661.1328125) * x * x * x * x * x * x * y * y * y * z * z + std::sqrt(6644.53125) * x * x * x * x * x * x * y * z * z * z * z - std::sqrt(934.38720703125) * x * x * x * x * y * y * y * y * y * z * z + std::sqrt(18457.03125) * x * x * x * x * y * y * y * z * z * z * z - std::sqrt(13650.8203125) * x * x * x * x * y * z * z * z * z * z * z + std::sqrt(738.28125) * x * x * y * y * y * y * y * z * z * z * z - std::sqrt(6067.03125) * x * x * y * y * y * z * z * z * z * z * z + std::sqrt(472.5) * x * x * y * z * z * z * z * z * z * z * z + std::sqrt(25.9552001953125) * y * y * y * y * y * y * y * y * y * z * z - std::sqrt(738.28125) * y * y * y * y * y * y * y * z * z * z * z + std::sqrt(1516.7578125) * y * y * y * y * y * z * z * z * z * z * z - std::sqrt(52.5) * y * y * y * z * z * z * z * z * z * z * z) + e_1 * (-std::sqrt(233.5968017578125) * x * x * x * x * x * x * x * x * y - std::sqrt(1661.1328125) * x * x * x * x * x * x * y * y * y - std::sqrt(3737.548828125) * x * x * x * x * x * x * y * z * z - std::sqrt(934.38720703125) * x * x * x * x * y * y * y * y * y - std::sqrt(10382.080078125) * x * x * x * x * y * y * y * z * z - std::sqrt(14950.1953125) * x * x * x * x * y * z * z * z * z - std::sqrt(415.283203125) * x * x * y * y * y * y * y * z * z - std::sqrt(6644.53125) * x * x * y * y * y * z * z * z * z - std::sqrt(106312.5) * x * x * y * z * z * z * z * z * z + std::sqrt(25.9552001953125) * y * y * y * y * y * y * y * y * y + std::sqrt(415.283203125) * y * y * y * y * y * y * y * z * z + std::sqrt(1661.1328125) * y * y * y * y * y * z * z * z * z + std::sqrt(11812.5) * y * y * y * z * z * z * z * z * z) + e_2 * (-std::sqrt(93438.720703125) * x * x * x * x * x * x * y - std::sqrt(259552.001953125) * x * x * x * x * y * y * y - std::sqrt(956812.5) * x * x * x * x * y * z * z - std::sqrt(10382.080078125) * x * x * y * y * y * y * y - std::sqrt(425250.0) * x * x * y * y * y * z * z - std::sqrt(8611312.5) * x * x * y * z * z * z * z + std::sqrt(10382.080078125) * y * y * y * y * y * y * y + std::sqrt(106312.5) * y * y * y * y * y * z * z + std::sqrt(956812.5) * y * y * y * z * z * z * z) + e_3 * (-std::sqrt(5588050.78125) * x * x * x * x * y - std::sqrt(2483578.125) * x * x * y * y * y - std::sqrt(71867250.0) * x * x * y * z * z + std::sqrt(620894.53125) * y * y * y * y * y + std::sqrt(7985250.0) * y * y * y * z * z) + e_4 * (-std::sqrt(46883812.5) * x * x * y + std::sqrt(5209312.5) * y * y * y);
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

        pc_39[k] = e_0 * (-std::sqrt(15.5731201171875) * x * x * x * x * x * x * x * x * x * y * z - std::sqrt(110.7421875) * x * x * x * x * x * x * x * y * y * y * z + std::sqrt(3349.951171875) * x * x * x * x * x * x * x * y * z * z * z - std::sqrt(62.29248046875) * x * x * x * x * x * y * y * y * y * y * z + std::sqrt(9305.419921875) * x * x * x * x * x * y * y * y * z * z * z - std::sqrt(24916.9921875) * x * x * x * x * x * y * z * z * z * z * z + std::sqrt(372.216796875) * x * x * x * y * y * y * y * y * z * z * z - std::sqrt(11074.21875) * x * x * x * y * y * y * z * z * z * z * z + std::sqrt(7087.5) * x * x * x * y * z * z * z * z * z * z * z + std::sqrt(1.7303466796875) * x * y * y * y * y * y * y * y * y * y * z - std::sqrt(372.216796875) * x * y * y * y * y * y * y * y * z * z * z + std::sqrt(2768.5546875) * x * y * y * y * y * y * z * z * z * z * z - std::sqrt(787.5) * x * y * y * y * z * z * z * z * z * z * z) + e_1 * (-std::sqrt(110.7421875) * x * x * x * x * x * y * y * y * z - std::sqrt(5426.3671875) * x * x * x * x * x * y * z * z * z - std::sqrt(442.96875) * x * x * x * y * y * y * y * y * z + std::sqrt(11074.21875) * x * x * x * y * y * y * z * z * z - std::sqrt(15946.875) * x * x * x * y * z * z * z * z * z - std::sqrt(110.7421875) * x * y * y * y * y * y * y * y * z + std::sqrt(32004.4921875) * x * y * y * y * y * y * z * z * z - std::sqrt(143521.875) * x * y * y * y * z * z * z * z * z + std::sqrt(28350.0) * x * y * z * z * z * z * z * z * z) + e_2 * (-std::sqrt(15946.875) * x * x * x * x * x * y * z - std::sqrt(708750.0) * x * x * x * y * z * z * z + std::sqrt(15946.875) * x * y * y * y * y * y * z - std::sqrt(708750.0) * x * y * y * y * z * z * z + std::sqrt(1020600.0) * x * y * z * z * z * z * z) + e_3 * (-std::sqrt(1594687.5) * x * x * x * y * z - std::sqrt(177187.5) * x * y * y * y * z + std::sqrt(2835000.0) * x * y * z * z * z);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        pc_40[k] = e_0 * (std::sqrt(436.04736328125) * x * x * x * x * x * x * x * x * y * z * z + std::sqrt(193.798828125) * x * x * x * x * x * x * y * y * y * z * z - std::sqrt(9496.142578125) * x * x * x * x * x * x * y * z * z * z * z - std::sqrt(775.1953125) * x * x * x * x * y * y * y * y * y * z * z + std::sqrt(1055.126953125) * x * x * x * x * y * y * y * z * z * z * z + std::sqrt(12403.125) * x * x * x * x * y * z * z * z * z * z * z - std::sqrt(193.798828125) * x * x * y * y * y * y * y * y * y * z * z + std::sqrt(9496.142578125) * x * x * y * y * y * y * y * z * z * z * z - std::sqrt(22050.0) * x * x * y * y * y * z * z * z * z * z * z + std::sqrt(48.44970703125) * y * y * y * y * y * y * y * y * y * z * z - std::sqrt(1055.126953125) * y * y * y * y * y * y * y * z * z * z * z + std::sqrt(1378.125) * y * y * y * y * y * z * z * z * z * z * z) + e_1 * (std::sqrt(436.04736328125) * x * x * x * x * x * x * x * x * y + std::sqrt(193.798828125) * x * x * x * x * x * x * y * y * y + std::sqrt(1744.189453125) * x * x * x * x * x * x * y * z * z - std::sqrt(775.1953125) * x * x * x * x * y * y * y * y * y + std::sqrt(15697.705078125) * x * x * x * x * y * y * y * z * z + std::sqrt(93798.6328125) * x * x * x * x * y * z * z * z * z - std::sqrt(193.798828125) * x * x * y * y * y * y * y * y * y + std::sqrt(15697.705078125) * x * x * y * y * y * y * y * z * z - std::sqrt(1119382.03125) * x * x * y * y * y * z * z * z * z + std::sqrt(49612.5) * x * x * y * z * z * z * z * z * z + std::sqrt(48.44970703125) * y * y * y * y * y * y * y * y * y + std::sqrt(1744.189453125) * y * y * y * y * y * y * y * z * z - std::sqrt(775.1953125) * y * y * y * y * y * z * z * z * z + std::sqrt(49612.5) * y * y * y * z * z * z * z * z * z) + e_2 * (std::sqrt(111628.125) * x * x * x * x * x * x * y + std::sqrt(775.1953125) * x * x * x * x * y * y * y + std::sqrt(2016283.0078125) * x * x * x * x * y * z * z - std::sqrt(27907.03125) * x * x * y * y * y * y * y - std::sqrt(4716288.28125) * x * x * y * y * y * z * z + std::sqrt(1004653.125) * x * x * y * z * z * z * z + std::sqrt(19379.8828125) * y * y * y * y * y * y * y + std::sqrt(174418.9453125) * y * y * y * y * y * z * z + std::sqrt(1004653.125) * y * y * y * z * z * z * z + std::sqrt(198450.0) * y * z * z * z * z * z * z) + e_3 * (std::sqrt(4716288.28125) * x * x * x * x * y - std::sqrt(1500778.125) * x * x * y * y * y + std::sqrt(7144200.0) * x * x * y * z * z + std::sqrt(1119382.03125) * y * y * y * y * y + std::sqrt(7144200.0) * y * y * y * z * z + std::sqrt(12700800.0) * y * z * z * z * z) + e_4 * (std::sqrt(9041878.125) * x * x * y + std::sqrt(9041878.125) * y * y * y + std::sqrt(64297800.0) * y * z * z) + e_5 * (std::sqrt(16074450.0) * y);
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

        pc_41[k] = e_0 * (std::sqrt(18.16864013671875) * x * x * x * x * x * x * x * x * x * y * z - std::sqrt(32.2998046875) * x * x * x * x * x * x * x * y * y * y * z - std::sqrt(2067.1875) * x * x * x * x * x * x * x * y * z * z * z - std::sqrt(395.672607421875) * x * x * x * x * x * y * y * y * y * y * z + std::sqrt(11254.6875) * x * x * x * x * x * y * y * y * z * z * z + std::sqrt(8268.75) * x * x * x * x * x * y * z * z * z * z * z - std::sqrt(32.2998046875) * x * x * x * y * y * y * y * y * y * y * z + std::sqrt(11254.6875) * x * x * x * y * y * y * y * y * z * z * z - std::sqrt(91875.0) * x * x * x * y * y * y * z * z * z * z * z + std::sqrt(18.16864013671875) * x * y * y * y * y * y * y * y * y * y * z - std::sqrt(2067.1875) * x * y * y * y * y * y * y * y * z * z * z + std::sqrt(8268.75) * x * y * y * y * y * y * z * z * z * z * z) + e_1 * (std::sqrt(74418.75) * x * x * x * x * x * y * z * z * z - std::sqrt(826875.0) * x * x * x * y * y * y * z * z * z + std::sqrt(74418.75) * x * y * y * y * y * y * z * z * z) + e_2 * (std::sqrt(167442.1875) * x * x * x * x * x * y * z - std::sqrt(1860468.75) * x * x * x * y * y * y * z + std::sqrt(167442.1875) * x * y * y * y * y * y * z);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        pc_42[k] = e_0 * (-std::sqrt(327.0355224609375) * x * x * x * x * x * x * x * x * y * z * z + std::sqrt(9302.34375) * x * x * x * x * x * x * y * y * y * z * z + std::sqrt(2325.5859375) * x * x * x * x * x * x * y * z * z * z * z + std::sqrt(3633.72802734375) * x * x * x * x * y * y * y * y * y * z * z - std::sqrt(93281.8359375) * x * x * x * x * y * y * y * z * z * z * z - std::sqrt(2325.5859375) * x * x * y * y * y * y * y * y * y * z * z + std::sqrt(20930.2734375) * x * x * y * y * y * y * y * z * z * z * z + std::sqrt(36.3372802734375) * y * y * y * y * y * y * y * y * y * z * z - std::sqrt(258.3984375) * y * y * y * y * y * y * y * z * z * z * z) + e_1 * (-std::sqrt(327.0355224609375) * x * x * x * x * x * x * x * x * y + std::sqrt(9302.34375) * x * x * x * x * x * x * y * y * y + std::sqrt(5232.568359375) * x * x * x * x * x * x * y * z * z + std::sqrt(3633.72802734375) * x * x * x * x * y * y * y * y * y + std::sqrt(47093.115234375) * x * x * x * x * y * y * y * z * z - std::sqrt(37209.375) * x * x * x * x * y * z * z * z * z - std::sqrt(2325.5859375) * x * x * y * y * y * y * y * y * y + std::sqrt(47093.115234375) * x * x * y * y * y * y * y * z * z - std::sqrt(148837.5) * x * x * y * y * y * z * z * z * z + std::sqrt(36.3372802734375) * y * y * y * y * y * y * y * y * y + std::sqrt(5232.568359375) * y * y * y * y * y * y * y * z * z - std::sqrt(37209.375) * y * y * y * y * y * z * z * z * z) + e_2 * (-std::sqrt(5232.568359375) * x * x * x * x * x * x * y + std::sqrt(1284304.833984375) * x * x * x * x * y * y * y + std::sqrt(83721.09375) * x * x * x * x * y * z * z - std::sqrt(47093.115234375) * x * x * y * y * y * y * y + std::sqrt(334884.375) * x * x * y * y * y * z * z - std::sqrt(1339537.5) * x * x * y * z * z * z * z + std::sqrt(14534.912109375) * y * y * y * y * y * y * y + std::sqrt(83721.09375) * y * y * y * y * y * z * z - std::sqrt(1339537.5) * y * y * y * z * z * z * z) + e_3 * (std::sqrt(753489.84375) * x * x * x * x * y + std::sqrt(3013959.375) * x * x * y * y * y - std::sqrt(1339537.5) * x * x * y * z * z + std::sqrt(753489.84375) * y * y * y * y * y - std::sqrt(1339537.5) * y * y * y * z * z - std::sqrt(2381400.0) * y * z * z * z * z) + e_4 * (std::sqrt(5358150.0) * x * x * y + std::sqrt(5358150.0) * y * y * y - std::sqrt(12055837.5) * y * z * z) + e_5 * (std::sqrt(1339537.5) * y);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        pc_43[k] = e_0 * (-std::sqrt(32.70355224609375) * x * x * x * x * x * x * x * x * x * y * z + std::sqrt(2848.8427734375) * x * x * x * x * x * x * x * y * y * y * z + std::sqrt(232.55859375) * x * x * x * x * x * x * x * y * z * z * z + std::sqrt(130.814208984375) * x * x * x * x * x * y * y * y * y * y * z - std::sqrt(24832.08984375) * x * x * x * x * x * y * y * y * z * z * z - std::sqrt(1453.4912109375) * x * x * x * y * y * y * y * y * y * y * z + std::sqrt(16149.90234375) * x * x * x * y * y * y * y * y * z * z * z + std::sqrt(90.84320068359375) * x * y * y * y * y * y * y * y * y * y * z - std::sqrt(645.99609375) * x * y * y * y * y * y * y * y * z * z * z) + e_1 * (std::sqrt(581396.484375) * x * x * x * x * x * y * y * y * z - std::sqrt(23255.859375) * x * x * x * x * x * y * z * z * z - std::sqrt(93023.4375) * x * x * x * y * y * y * y * y * z - std::sqrt(93023.4375) * x * x * x * y * y * y * z * z * z + std::sqrt(23255.859375) * x * y * y * y * y * y * y * y * z - std::sqrt(23255.859375) * x * y * y * y * y * y * z * z * z) + e_2 * (std::sqrt(837210.9375) * x * x * x * x * x * y * z + std::sqrt(3348843.75) * x * x * x * y * y * y * z - std::sqrt(1488375.0) * x * x * x * y * z * z * z + std::sqrt(837210.9375) * x * y * y * y * y * y * z - std::sqrt(1488375.0) * x * y * y * y * z * z * z) + e_3 * (std::sqrt(13395375.0) * x * x * x * y * z + std::sqrt(13395375.0) * x * y * y * y * z - std::sqrt(5953500.0) * x * y * z * z * z) + e_4 * (std::sqrt(13395375.0) * x * y * z);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        pc_44[k] = e_0 * (std::sqrt(10.09368896484375) * x * x * x * x * x * x * x * x * x * y * y - std::sqrt(2583.984375) * x * x * x * x * x * x * x * y * y * z * z - std::sqrt(79.134521484375) * x * x * x * x * x * y * y * y * y * y * y + std::sqrt(2583.984375) * x * x * x * x * x * y * y * y * y * z * z + std::sqrt(2583.984375) * x * x * x * x * x * y * y * z * z * z * z - std::sqrt(25.83984375) * x * x * x * y * y * y * y * y * y * y * y + std::sqrt(8372.109375) * x * x * x * y * y * y * y * y * y * z * z - std::sqrt(10335.9375) * x * x * x * y * y * y * y * z * z * z * z + std::sqrt(0.40374755859375) * x * y * y * y * y * y * y * y * y * y * y - std::sqrt(103.359375) * x * y * y * y * y * y * y * y * y * z * z + std::sqrt(103.359375) * x * y * y * y * y * y * y * z * z * z * z) + e_1 * (std::sqrt(10.09368896484375) * x * x * x * x * x * x * x * x * x + std::sqrt(4037.4755859375) * x * x * x * x * x * x * x * y * y - std::sqrt(2583.984375) * x * x * x * x * x * x * x * z * z - std::sqrt(6823.333740234375) * x * x * x * x * x * y * y * y * y - std::sqrt(209302.734375) * x * x * x * x * x * y * y * z * z + std::sqrt(2583.984375) * x * x * x * x * x * z * z * z * z - std::sqrt(19541.3818359375) * x * x * x * y * y * y * y * y * y + std::sqrt(1614990.234375) * x * x * x * y * y * y * y * z * z - std::sqrt(10335.9375) * x * x * x * y * y * z * z * z * z + std::sqrt(10.09368896484375) * x * y * y * y * y * y * y * y * y + std::sqrt(2583.984375) * x * y * y * y * y * y * y * z * z - std::sqrt(23255.859375) * x * y * y * y * y * z * z * z * z) + e_2 * (std::sqrt(4037.4755859375) * x * x * x * x * x * x * x + std::sqrt(36337.2802734375) * x * x * x * x * x * y * y - std::sqrt(372093.75) * x * x * x * x * x * z * z - std::sqrt(682333.3740234375) * x * x * x * y * y * y * y + std::sqrt(1488375.0) * x * x * x * y * y * z * z + std::sqrt(41343.75) * x * x * x * z * z * z * z - std::sqrt(19541.3818359375) * x * y * y * y * y * y * y + std::sqrt(3348843.75) * x * y * y * y * y * z * z - std::sqrt(372093.75) * x * y * y * z * z * z * z) + e_3 * (std::sqrt(93023.4375) * x * x * x * x * x - std::sqrt(372093.75) * x * x * x * y * y - std::sqrt(1488375.0) * x * x * x * z * z - std::sqrt(837210.9375) * x * y * y * y * y + std::sqrt(13395375.0) * x * y * y * z * z) + e_4 * (std::sqrt(93023.4375) * x * x * x - std::sqrt(837210.9375) * x * y * y);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        pc_45[k] = e_0 * (std::sqrt(64.599609375) * x * x * x * x * x * x * x * x * y * y * z + std::sqrt(64.599609375) * x * x * x * x * x * x * y * y * y * y * z - std::sqrt(16537.5) * x * x * x * x * x * x * y * y * z * z * z - std::sqrt(64.599609375) * x * x * x * x * y * y * y * y * y * y * z + std::sqrt(16537.5) * x * x * x * x * y * y * z * z * z * z * z - std::sqrt(64.599609375) * x * x * y * y * y * y * y * y * y * y * z + std::sqrt(16537.5) * x * x * y * y * y * y * y * y * z * z * z - std::sqrt(16537.5) * x * x * y * y * y * y * z * z * z * z * z) + e_1 * (std::sqrt(64.599609375) * x * x * x * x * x * x * x * x * z - std::sqrt(12661.5234375) * x * x * x * x * x * x * y * y * z - std::sqrt(16537.5) * x * x * x * x * x * x * z * z * z - std::sqrt(413437.5) * x * x * x * x * y * y * z * z * z + std::sqrt(16537.5) * x * x * x * x * z * z * z * z * z + std::sqrt(12661.5234375) * x * x * y * y * y * y * y * y * z + std::sqrt(413437.5) * x * x * y * y * y * y * z * z * z - std::sqrt(64.599609375) * y * y * y * y * y * y * y * y * z + std::sqrt(16537.5) * y * y * y * y * y * y * z * z * z - std::sqrt(16537.5) * y * y * y * y * z * z * z * z * z) + e_2 * (-std::sqrt(18669.287109375) * x * x * x * x * x * x * z - std::sqrt(3270355.224609375) * x * x * x * x * y * y * z - std::sqrt(413437.5) * x * x * x * x * z * z * z + std::sqrt(3270355.224609375) * x * x * y * y * y * y * z + std::sqrt(148837.5) * x * x * z * z * z * z * z + std::sqrt(18669.287109375) * y * y * y * y * y * y * z + std::sqrt(413437.5) * y * y * y * y * z * z * z - std::sqrt(148837.5) * y * y * z * z * z * z * z) + e_3 * (-std::sqrt(3720937.5) * x * x * x * x * z + std::sqrt(3720937.5) * y * y * y * y * z) + e_4 * (-std::sqrt(8372109.375) * x * x * z + std::sqrt(8372109.375) * y * y * z);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        pc_46[k] = e_0 * (-std::sqrt(2.01873779296875) * x * x * x * x * x * x * x * x * x * y * y - std::sqrt(14.35546875) * x * x * x * x * x * x * x * y * y * y * y + std::sqrt(1162.79296875) * x * x * x * x * x * x * x * y * y * z * z - std::sqrt(8.074951171875) * x * x * x * x * x * y * y * y * y * y * y + std::sqrt(3229.98046875) * x * x * x * x * x * y * y * y * y * z * z - std::sqrt(41860.546875) * x * x * x * x * x * y * y * z * z * z * z + std::sqrt(129.19921875) * x * x * x * y * y * y * y * y * y * z * z - std::sqrt(18604.6875) * x * x * x * y * y * y * y * z * z * z * z + std::sqrt(33075.0) * x * x * x * y * y * z * z * z * z * z * z + std::sqrt(0.22430419921875) * x * y * y * y * y * y * y * y * y * y * y - std::sqrt(129.19921875) * x * y * y * y * y * y * y * y * y * z * z + std::sqrt(4651.171875) * x * y * y * y * y * y * y * z * z * z * z - std::sqrt(3675.0) * x * y * y * y * y * z * z * z * z * z * z) + e_1 * (-std::sqrt(2.01873779296875) * x * x * x * x * x * x * x * x * x - std::sqrt(1582.6904296875) * x * x * x * x * x * x * x * y * y + std::sqrt(1162.79296875) * x * x * x * x * x * x * x * z * z - std::sqrt(4271.649169921875) * x * x * x * x * x * y * y * y * y - std::sqrt(29069.82421875) * x * x * x * x * x * y * y * z * z - std::sqrt(41860.546875) * x * x * x * x * x * z * z * z * z - std::sqrt(290.6982421875) * x * x * x * y * y * y * y * y * y - std::sqrt(6330.76171875) * x * x * x * y * y * y * y * z * z - std::sqrt(18604.6875) * x * x * x * y * y * z * z * z * z + std::sqrt(33075.0) * x * x * x * z * z * z * z * z * z + std::sqrt(98.91815185546875) * x * y * y * y * y * y * y * y * y + std::sqrt(15633.10546875) * x * y * y * y * y * y * y * z * z - std::sqrt(87338.671875) * x * y * y * y * y * z * z * z * z + std::sqrt(33075.0) * x * y * y * z * z * z * z * z * z) + e_2 * (-std::sqrt(807.4951171875) * x * x * x * x * x * x * x - std::sqrt(488663.7451171875) * x * x * x * x * x * y * y - std::sqrt(116279.296875) * x * x * x * x * x * z * z - std::sqrt(244477.2216796875) * x * x * x * y * y * y * y - std::sqrt(2251167.1875) * x * x * x * y * y * z * z + std::sqrt(74418.75) * x * x * x * z * z * z * z + std::sqrt(31040.1123046875) * x * y * y * y * y * y * y - std::sqrt(4651.171875) * x * y * y * y * y * z * z + std::sqrt(74418.75) * x * y * y * z * z * z * z + std::sqrt(132300.0) * x * z * z * z * z * z * z) + e_3 * (-std::sqrt(297675.0) * x * x * x * x * x - std::sqrt(14586075.0) * x * x * x * y * y - std::sqrt(1190700.0) * x * x * x * z * z + std::sqrt(132300.0) * x * y * y * y * y - std::sqrt(1190700.0) * x * y * y * z * z + std::sqrt(4762800.0) * x * z * z * z * z) + e_4 * (-std::sqrt(6716292.1875) * x * x * x - std::sqrt(6716292.1875) * x * y * y + std::sqrt(4762800.0) * x * z * z) + e_5 * (-std::sqrt(2679075.0) * x);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        pc_47[k] = e_0 * (-std::sqrt(21.533203125) * x * x * x * x * x * x * x * x * y * y * z - std::sqrt(193.798828125) * x * x * x * x * x * x * y * y * y * y * z + std::sqrt(6976.7578125) * x * x * x * x * x * x * y * y * z * z * z - std::sqrt(193.798828125) * x * x * x * x * y * y * y * y * y * y * z + std::sqrt(27907.03125) * x * x * x * x * y * y * y * y * z * z * z - std::sqrt(49612.5) * x * x * x * x * y * y * z * z * z * z * z - std::sqrt(21.533203125) * x * x * y * y * y * y * y * y * y * y * z + std::sqrt(6976.7578125) * x * x * y * y * y * y * y * y * z * z * z - std::sqrt(49612.5) * x * x * y * y * y * y * z * z * z * z * z + std::sqrt(22050.0) * x * x * y * y * z * z * z * z * z * z * z) + e_1 * (-std::sqrt(21.533203125) * x * x * x * x * x * x * x * x * z + std::sqrt(1378.125) * x * x * x * x * x * x * y * y * z + std::sqrt(6976.7578125) * x * x * x * x * x * x * z * z * z + std::sqrt(6976.7578125) * x * x * x * x * y * y * y * y * z + std::sqrt(775.1953125) * x * x * x * x * y * y * z * z * z - std::sqrt(49612.5) * x * x * x * x * z * z * z * z * z + std::sqrt(1378.125) * x * x * y * y * y * y * y * y * z + std::sqrt(775.1953125) * x * x * y * y * y * y * z * z * z + std::sqrt(198450.0) * x * x * y * y * z * z * z * z * z + std::sqrt(22050.0) * x * x * z * z * z * z * z * z * z - std::sqrt(21.533203125) * y * y * y * y * y * y * y * y * z + std::sqrt(6976.7578125) * y * y * y * y * y * y * z * z * z - std::sqrt(49612.5) * y * y * y * y * z * z * z * z * z + std::sqrt(22050.0) * y * y * z * z * z * z * z * z * z) + e_2 * (std::sqrt(6223.095703125) * x * x * x * x * x * x * z + std::sqrt(325775.830078125) * x * x * x * x * y * y * z - std::sqrt(224031.4453125) * x * x * x * x * z * z * z + std::sqrt(325775.830078125) * x * x * y * y * y * y * z + std::sqrt(5733344.53125) * x * x * y * y * z * z * z + std::sqrt(1240312.5) * x * x * z * z * z * z * z + std::sqrt(6223.095703125) * y * y * y * y * y * y * z - std::sqrt(224031.4453125) * y * y * y * y * z * z * z + std::sqrt(1240312.5) * y * y * z * z * z * z * z + std::sqrt(22050.0) * z * z * z * z * z * z * z) + e_3 * (std::sqrt(12403.125) * x * x * x * x * z + std::sqrt(21879112.5) * x * x * y * y * z + std::sqrt(12700800.0) * x * x * z * z * z + std::sqrt(12403.125) * y * y * y * y * z + std::sqrt(12700800.0) * y * y * z * z * z + std::sqrt(3175200.0) * z * z * z * z * z) + e_4 * (std::sqrt(16074450.0) * x * x * z + std::sqrt(16074450.0) * y * y * z + std::sqrt(38896200.0) * z * z * z) + e_5 * (std::sqrt(28576800.0) * z);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        pc_48[k] = e_0 * (std::sqrt(0.1922607421875) * x * x * x * x * x * x * x * x * x * y * y + std::sqrt(3.076171875) * x * x * x * x * x * x * x * y * y * y * y - std::sqrt(150.732421875) * x * x * x * x * x * x * x * y * y * z * z + std::sqrt(6.92138671875) * x * x * x * x * x * y * y * y * y * y * y - std::sqrt(1356.591796875) * x * x * x * x * x * y * y * y * y * z * z + std::sqrt(8970.1171875) * x * x * x * x * x * y * y * z * z * z * z + std::sqrt(3.076171875) * x * x * x * y * y * y * y * y * y * y * y - std::sqrt(1356.591796875) * x * x * x * y * y * y * y * y * y * z * z + std::sqrt(35880.46875) * x * x * x * y * y * y * y * z * z * z * z - std::sqrt(19687.5) * x * x * x * y * y * z * z * z * z * z * z + std::sqrt(0.1922607421875) * x * y * y * y * y * y * y * y * y * y * y - std::sqrt(150.732421875) * x * y * y * y * y * y * y * y * y * z * z + std::sqrt(8970.1171875) * x * y * y * y * y * y * y * z * z * z * z - std::sqrt(19687.5) * x * y * y * y * y * z * z * z * z * z * z + std::sqrt(3150.0) * x * y * y * z * z * z * z * z * z * z * z) + e_1 * (std::sqrt(0.1922607421875) * x * x * x * x * x * x * x * x * x + std::sqrt(196.875) * x * x * x * x * x * x * x * y * y - std::sqrt(150.732421875) * x * x * x * x * x * x * x * z * z + std::sqrt(1557.31201171875) * x * x * x * x * x * y * y * y * y + std::sqrt(8001.123046875) * x * x * x * x * x * y * y * z * z + std::sqrt(8970.1171875) * x * x * x * x * x * z * z * z * z + std::sqrt(1488.8671875) * x * x * x * y * y * y * y * y * y + std::sqrt(46539.404296875) * x * x * x * y * y * y * y * z * z - std::sqrt(442.96875) * x * x * x * y * y * z * z * z * z - std::sqrt(19687.5) * x * x * x * z * z * z * z * z * z + std::sqrt(161.6912841796875) * x * y * y * y * y * y * y * y * y + std::sqrt(12996.826171875) * x * y * y * y * y * y * y * z * z - std::sqrt(13399.8046875) * x * y * y * y * y * z * z * z * z + std::sqrt(95287.5) * x * y * y * z * z * z * z * z * z + std::sqrt(3150.0) * x * z * z * z * z * z * z * z * z) + e_2 * (std::sqrt(76.904296875) * x * x * x * x * x * x * x + std::sqrt(96373.388671875) * x * x * x * x * x * y * y + std::sqrt(39977.9296875) * x * x * x * x * x * z * z + std::sqrt(353516.748046875) * x * x * x * y * y * y * y + std::sqrt(1063567.96875) * x * x * x * y * y * z * z - std::sqrt(347287.5) * x * x * x * z * z * z * z + std::sqrt(85791.357421875) * x * y * y * y * y * y * y + std::sqrt(691141.9921875) * x * y * y * y * y * z * z + std::sqrt(3749287.5) * x * y * y * z * z * z * z + std::sqrt(532350.0) * x * z * z * z * z * z * z) + e_3 * (std::sqrt(53599.21875) * x * x * x * x * x + std::sqrt(6593146.875) * x * x * x * y * y - std::sqrt(28350.0) * x * x * x * z * z + std::sqrt(5457817.96875) * x * y * y * y * y + std::sqrt(43120350.0) * x * y * y * z * z + std::sqrt(13721400.0) * x * z * z * z * z) + e_4 * (std::sqrt(1389150.0) * x * x * x + std::sqrt(50009400.0) * x * y * y + std::sqrt(50009400.0) * x * z * z) + e_5 * (std::sqrt(12502350.0) * x);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        pc_49[k] = e_0 * (std::sqrt(2.8839111328125) * x * x * x * x * x * x * x * x * x * y * z + std::sqrt(46.142578125) * x * x * x * x * x * x * x * y * y * y * z - std::sqrt(1004.8828125) * x * x * x * x * x * x * x * y * z * z * z + std::sqrt(103.82080078125) * x * x * x * x * x * y * y * y * y * y * z - std::sqrt(9043.9453125) * x * x * x * x * x * y * y * y * z * z * z + std::sqrt(10107.0703125) * x * x * x * x * x * y * z * z * z * z * z + std::sqrt(46.142578125) * x * x * x * y * y * y * y * y * y * y * z - std::sqrt(9043.9453125) * x * x * x * y * y * y * y * y * z * z * z + std::sqrt(40428.28125) * x * x * x * y * y * y * z * z * z * z * z - std::sqrt(7560.0) * x * x * x * y * z * z * z * z * z * z * z + std::sqrt(2.8839111328125) * x * y * y * y * y * y * y * y * y * y * z - std::sqrt(1004.8828125) * x * y * y * y * y * y * y * y * z * z * z + std::sqrt(10107.0703125) * x * y * y * y * y * y * z * z * z * z * z - std::sqrt(7560.0) * x * y * y * y * z * z * z * z * z * z * z + std::sqrt(210.0) * x * y * z * z * z * z * z * z * z * z * z) + e_1 * (-std::sqrt(184.5703125) * x * x * x * x * x * x * x * y * z - std::sqrt(1661.1328125) * x * x * x * x * x * y * y * y * z + std::sqrt(2953.125) * x * x * x * x * x * y * z * z * z - std::sqrt(1661.1328125) * x * x * x * y * y * y * y * y * z + std::sqrt(11812.5) * x * x * x * y * y * y * z * z * z - std::sqrt(47250.0) * x * x * x * y * z * z * z * z * z - std::sqrt(184.5703125) * x * y * y * y * y * y * y * y * z + std::sqrt(2953.125) * x * y * y * y * y * y * z * z * z - std::sqrt(47250.0) * x * y * y * y * z * z * z * z * z) + e_2 * (-std::sqrt(14950.1953125) * x * x * x * x * x * y * z - std::sqrt(59800.78125) * x * x * x * y * y * y * z - std::sqrt(425250.0) * x * x * x * y * z * z * z - std::sqrt(14950.1953125) * x * y * y * y * y * y * z - std::sqrt(425250.0) * x * y * y * y * z * z * z - std::sqrt(425250.0) * x * y * z * z * z * z * z) + e_3 * (-std::sqrt(1701000.0) * x * x * x * y * z - std::sqrt(1701000.0) * x * y * y * y * z - std::sqrt(12096000.0) * x * y * z * z * z) + e_4 * (-std::sqrt(20837250.0) * x * y * z);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        pc_50[k] = e_0 * (std::sqrt(0.1922607421875) * x * x * x * x * x * x * x * x * x * x * y + std::sqrt(3.076171875) * x * x * x * x * x * x * x * x * y * y * y - std::sqrt(150.732421875) * x * x * x * x * x * x * x * x * y * z * z + std::sqrt(6.92138671875) * x * x * x * x * x * x * y * y * y * y * y - std::sqrt(1356.591796875) * x * x * x * x * x * x * y * y * y * z * z + std::sqrt(8970.1171875) * x * x * x * x * x * x * y * z * z * z * z + std::sqrt(3.076171875) * x * x * x * x * y * y * y * y * y * y * y - std::sqrt(1356.591796875) * x * x * x * x * y * y * y * y * y * z * z + std::sqrt(35880.46875) * x * x * x * x * y * y * y * z * z * z * z - std::sqrt(19687.5) * x * x * x * x * y * z * z * z * z * z * z + std::sqrt(0.1922607421875) * x * x * y * y * y * y * y * y * y * y * y - std::sqrt(150.732421875) * x * x * y * y * y * y * y * y * y * z * z + std::sqrt(8970.1171875) * x * x * y * y * y * y * y * z * z * z * z - std::sqrt(19687.5) * x * x * y * y * y * z * z * z * z * z * z + std::sqrt(3150.0) * x * x * y * z * z * z * z * z * z * z * z) + e_1 * (std::sqrt(161.6912841796875) * x * x * x * x * x * x * x * x * y + std::sqrt(1488.8671875) * x * x * x * x * x * x * y * y * y + std::sqrt(12996.826171875) * x * x * x * x * x * x * y * z * z + std::sqrt(1557.31201171875) * x * x * x * x * y * y * y * y * y + std::sqrt(46539.404296875) * x * x * x * x * y * y * y * z * z - std::sqrt(13399.8046875) * x * x * x * x * y * z * z * z * z + std::sqrt(196.875) * x * x * y * y * y * y * y * y * y + std::sqrt(8001.123046875) * x * x * y * y * y * y * y * z * z - std::sqrt(442.96875) * x * x * y * y * y * z * z * z * z + std::sqrt(95287.5) * x * x * y * z * z * z * z * z * z + std::sqrt(0.1922607421875) * y * y * y * y * y * y * y * y * y - std::sqrt(150.732421875) * y * y * y * y * y * y * y * z * z + std::sqrt(8970.1171875) * y * y * y * y * y * z * z * z * z - std::sqrt(19687.5) * y * y * y * z * z * z * z * z * z + std::sqrt(3150.0) * y * z * z * z * z * z * z * z * z) + e_2 * (std::sqrt(85791.357421875) * x * x * x * x * x * x * y + std::sqrt(353516.748046875) * x * x * x * x * y * y * y + std::sqrt(691141.9921875) * x * x * x * x * y * z * z + std::sqrt(96373.388671875) * x * x * y * y * y * y * y + std::sqrt(1063567.96875) * x * x * y * y * y * z * z + std::sqrt(3749287.5) * x * x * y * z * z * z * z + std::sqrt(76.904296875) * y * y * y * y * y * y * y + std::sqrt(39977.9296875) * y * y * y * y * y * z * z - std::sqrt(347287.5) * y * y * y * z * z * z * z + std::sqrt(532350.0) * y * z * z * z * z * z * z) + e_3 * (std::sqrt(5457817.96875) * x * x * x * x * y + std::sqrt(6593146.875) * x * x * y * y * y + std::sqrt(43120350.0) * x * x * y * z * z + std::sqrt(53599.21875) * y * y * y * y * y - std::sqrt(28350.0) * y * y * y * z * z + std::sqrt(13721400.0) * y * z * z * z * z) + e_4 * (std::sqrt(50009400.0) * x * x * y + std::sqrt(1389150.0) * y * y * y + std::sqrt(50009400.0) * y * z * z) + e_5 * (std::sqrt(12502350.0) * y);
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

        pc_51[k] = e_0 * (-std::sqrt(5.38330078125) * x * x * x * x * x * x * x * x * x * y * z - std::sqrt(21.533203125) * x * x * x * x * x * x * x * y * y * y * z + std::sqrt(1744.189453125) * x * x * x * x * x * x * x * y * z * z * z + std::sqrt(1744.189453125) * x * x * x * x * x * y * y * y * z * z * z - std::sqrt(12403.125) * x * x * x * x * x * y * z * z * z * z * z + std::sqrt(21.533203125) * x * x * x * y * y * y * y * y * y * y * z - std::sqrt(1744.189453125) * x * x * x * y * y * y * y * y * z * z * z + std::sqrt(5512.5) * x * x * x * y * z * z * z * z * z * z * z + std::sqrt(5.38330078125) * x * y * y * y * y * y * y * y * y * y * z - std::sqrt(1744.189453125) * x * y * y * y * y * y * y * y * z * z * z + std::sqrt(12403.125) * x * y * y * y * y * y * z * z * z * z * z - std::sqrt(5512.5) * x * y * y * y * z * z * z * z * z * z * z) + e_1 * (std::sqrt(775.1953125) * x * x * x * x * x * x * x * y * z + std::sqrt(775.1953125) * x * x * x * x * x * y * y * y * z - std::sqrt(12403.125) * x * x * x * x * x * y * z * z * z - std::sqrt(775.1953125) * x * x * x * y * y * y * y * y * z + std::sqrt(198450.0) * x * x * x * y * z * z * z * z * z - std::sqrt(775.1953125) * x * y * y * y * y * y * y * y * z + std::sqrt(12403.125) * x * y * y * y * y * y * z * z * z - std::sqrt(198450.0) * x * y * y * y * z * z * z * z * z) + e_2 * (std::sqrt(27907.03125) * x * x * x * x * x * y * z + std::sqrt(2790703.125) * x * x * x * y * z * z * z - std::sqrt(27907.03125) * x * y * y * y * y * y * z - std::sqrt(2790703.125) * x * y * y * y * z * z * z) + e_3 * (std::sqrt(4961250.0) * x * x * x * y * z - std::sqrt(4961250.0) * x * y * y * y * z);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        pc_52[k] = e_0 * (-std::sqrt(0.22430419921875) * x * x * x * x * x * x * x * x * x * x * y + std::sqrt(129.19921875) * x * x * x * x * x * x * x * x * y * z * z + std::sqrt(8.074951171875) * x * x * x * x * x * x * y * y * y * y * y - std::sqrt(129.19921875) * x * x * x * x * x * x * y * y * y * z * z - std::sqrt(4651.171875) * x * x * x * x * x * x * y * z * z * z * z + std::sqrt(14.35546875) * x * x * x * x * y * y * y * y * y * y * y - std::sqrt(3229.98046875) * x * x * x * x * y * y * y * y * y * z * z + std::sqrt(18604.6875) * x * x * x * x * y * y * y * z * z * z * z + std::sqrt(3675.0) * x * x * x * x * y * z * z * z * z * z * z + std::sqrt(2.01873779296875) * x * x * y * y * y * y * y * y * y * y * y - std::sqrt(1162.79296875) * x * x * y * y * y * y * y * y * y * z * z + std::sqrt(41860.546875) * x * x * y * y * y * y * y * z * z * z * z - std::sqrt(33075.0) * x * x * y * y * y * z * z * z * z * z * z) + e_1 * (-std::sqrt(98.91815185546875) * x * x * x * x * x * x * x * x * y + std::sqrt(290.6982421875) * x * x * x * x * x * x * y * y * y - std::sqrt(15633.10546875) * x * x * x * x * x * x * y * z * z + std::sqrt(4271.649169921875) * x * x * x * x * y * y * y * y * y + std::sqrt(6330.76171875) * x * x * x * x * y * y * y * z * z + std::sqrt(87338.671875) * x * x * x * x * y * z * z * z * z + std::sqrt(1582.6904296875) * x * x * y * y * y * y * y * y * y + std::sqrt(29069.82421875) * x * x * y * y * y * y * y * z * z + std::sqrt(18604.6875) * x * x * y * y * y * z * z * z * z - std::sqrt(33075.0) * x * x * y * z * z * z * z * z * z + std::sqrt(2.01873779296875) * y * y * y * y * y * y * y * y * y - std::sqrt(1162.79296875) * y * y * y * y * y * y * y * z * z + std::sqrt(41860.546875) * y * y * y * y * y * z * z * z * z - std::sqrt(33075.0) * y * y * y * z * z * z * z * z * z) + e_2 * (-std::sqrt(31040.1123046875) * x * x * x * x * x * x * y + std::sqrt(244477.2216796875) * x * x * x * x * y * y * y + std::sqrt(4651.171875) * x * x * x * x * y * z * z + std::sqrt(488663.7451171875) * x * x * y * y * y * y * y + std::sqrt(2251167.1875) * x * x * y * y * y * z * z - std::sqrt(74418.75) * x * x * y * z * z * z * z + std::sqrt(807.4951171875) * y * y * y * y * y * y * y + std::sqrt(116279.296875) * y * y * y * y * y * z * z - std::sqrt(74418.75) * y * y * y * z * z * z * z - std::sqrt(132300.0) * y * z * z * z * z * z * z) + e_3 * (-std::sqrt(132300.0) * x * x * x * x * y + std::sqrt(14586075.0) * x * x * y * y * y + std::sqrt(1190700.0) * x * x * y * z * z + std::sqrt(297675.0) * y * y * y * y * y + std::sqrt(1190700.0) * y * y * y * z * z - std::sqrt(4762800.0) * y * z * z * z * z) + e_4 * (std::sqrt(6716292.1875) * x * x * y + std::sqrt(6716292.1875) * y * y * y - std::sqrt(4762800.0) * y * z * z) + e_5 * (std::sqrt(2679075.0) * y);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        pc_53[k] = e_0 * (std::sqrt(4.0374755859375) * x * x * x * x * x * x * x * x * x * y * z - std::sqrt(64.599609375) * x * x * x * x * x * x * x * y * y * y * z - std::sqrt(1033.59375) * x * x * x * x * x * x * x * y * z * z * z - std::sqrt(403.74755859375) * x * x * x * x * x * y * y * y * y * y * z + std::sqrt(25839.84375) * x * x * x * x * x * y * y * y * z * z * z + std::sqrt(1033.59375) * x * x * x * x * x * y * z * z * z * z * z - std::sqrt(64.599609375) * x * x * x * y * y * y * y * y * y * y * z + std::sqrt(25839.84375) * x * x * x * y * y * y * y * y * z * z * z - std::sqrt(37209.375) * x * x * x * y * y * y * z * z * z * z * z + std::sqrt(4.0374755859375) * x * y * y * y * y * y * y * y * y * y * z - std::sqrt(1033.59375) * x * y * y * y * y * y * y * y * z * z * z + std::sqrt(1033.59375) * x * y * y * y * y * y * z * z * z * z * z) + e_1 * (-std::sqrt(2325.5859375) * x * x * x * x * x * x * x * y * z + std::sqrt(12661.5234375) * x * x * x * x * x * y * y * y * z + std::sqrt(16537.5) * x * x * x * x * x * y * z * z * z + std::sqrt(12661.5234375) * x * x * x * y * y * y * y * y * z + std::sqrt(1653750.0) * x * x * x * y * y * y * z * z * z - std::sqrt(66150.0) * x * x * x * y * z * z * z * z * z - std::sqrt(2325.5859375) * x * y * y * y * y * y * y * y * z + std::sqrt(16537.5) * x * y * y * y * y * y * z * z * z - std::sqrt(66150.0) * x * y * y * y * z * z * z * z * z) + e_2 * (-std::sqrt(20930.2734375) * x * x * x * x * x * y * z + std::sqrt(9328183.59375) * x * x * x * y * y * y * z + std::sqrt(1653750.0) * x * x * x * y * z * z * z - std::sqrt(20930.2734375) * x * y * y * y * y * y * z + std::sqrt(1653750.0) * x * y * y * y * z * z * z - std::sqrt(595350.0) * x * y * z * z * z * z * z) + e_3 * (std::sqrt(14883750.0) * x * x * x * y * z + std::sqrt(14883750.0) * x * y * y * y * z) + e_4 * (std::sqrt(33488437.5) * x * y * z);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        pc_54[k] = e_0 * (std::sqrt(0.40374755859375) * x * x * x * x * x * x * x * x * x * x * y - std::sqrt(25.83984375) * x * x * x * x * x * x * x * x * y * y * y - std::sqrt(103.359375) * x * x * x * x * x * x * x * x * y * z * z - std::sqrt(79.134521484375) * x * x * x * x * x * x * y * y * y * y * y + std::sqrt(8372.109375) * x * x * x * x * x * x * y * y * y * z * z + std::sqrt(103.359375) * x * x * x * x * x * x * y * z * z * z * z + std::sqrt(2583.984375) * x * x * x * x * y * y * y * y * y * z * z - std::sqrt(10335.9375) * x * x * x * x * y * y * y * z * z * z * z + std::sqrt(10.09368896484375) * x * x * y * y * y * y * y * y * y * y * y - std::sqrt(2583.984375) * x * x * y * y * y * y * y * y * y * z * z + std::sqrt(2583.984375) * x * x * y * y * y * y * y * z * z * z * z) + e_1 * (std::sqrt(10.09368896484375) * x * x * x * x * x * x * x * x * y - std::sqrt(19541.3818359375) * x * x * x * x * x * x * y * y * y + std::sqrt(2583.984375) * x * x * x * x * x * x * y * z * z - std::sqrt(6823.333740234375) * x * x * x * x * y * y * y * y * y + std::sqrt(1614990.234375) * x * x * x * x * y * y * y * z * z - std::sqrt(23255.859375) * x * x * x * x * y * z * z * z * z + std::sqrt(4037.4755859375) * x * x * y * y * y * y * y * y * y - std::sqrt(209302.734375) * x * x * y * y * y * y * y * z * z - std::sqrt(10335.9375) * x * x * y * y * y * z * z * z * z + std::sqrt(10.09368896484375) * y * y * y * y * y * y * y * y * y - std::sqrt(2583.984375) * y * y * y * y * y * y * y * z * z + std::sqrt(2583.984375) * y * y * y * y * y * z * z * z * z) + e_2 * (-std::sqrt(19541.3818359375) * x * x * x * x * x * x * y - std::sqrt(682333.3740234375) * x * x * x * x * y * y * y + std::sqrt(3348843.75) * x * x * x * x * y * z * z + std::sqrt(36337.2802734375) * x * x * y * y * y * y * y + std::sqrt(1488375.0) * x * x * y * y * y * z * z - std::sqrt(372093.75) * x * x * y * z * z * z * z + std::sqrt(4037.4755859375) * y * y * y * y * y * y * y - std::sqrt(372093.75) * y * y * y * y * y * z * z + std::sqrt(41343.75) * y * y * y * z * z * z * z) + e_3 * (-std::sqrt(837210.9375) * x * x * x * x * y - std::sqrt(372093.75) * x * x * y * y * y + std::sqrt(13395375.0) * x * x * y * z * z + std::sqrt(93023.4375) * y * y * y * y * y - std::sqrt(1488375.0) * y * y * y * z * z) + e_4 * (-std::sqrt(837210.9375) * x * x * y + std::sqrt(93023.4375) * y * y * y);
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

        pc_55[k] = e_0 * (std::sqrt(100.9368896484375) * x * x * x * x * x * x * x * x * y * y * z - std::sqrt(1614.990234375) * x * x * x * x * x * x * y * y * z * z * z - std::sqrt(791.34521484375) * x * x * x * x * y * y * y * y * y * y * z + std::sqrt(1614.990234375) * x * x * x * x * y * y * y * y * z * z * z + std::sqrt(258.3984375) * x * x * x * x * y * y * z * z * z * z * z - std::sqrt(258.3984375) * x * x * y * y * y * y * y * y * y * y * z + std::sqrt(5232.568359375) * x * x * y * y * y * y * y * y * z * z * z - std::sqrt(1033.59375) * x * x * y * y * y * y * z * z * z * z * z + std::sqrt(4.0374755859375) * y * y * y * y * y * y * y * y * y * y * z - std::sqrt(64.599609375) * y * y * y * y * y * y * y * y * z * z * z + std::sqrt(10.3359375) * y * y * y * y * y * y * z * z * z * z * z) + e_1 * (std::sqrt(100.9368896484375) * x * x * x * x * x * x * x * x * z + std::sqrt(25839.84375) * x * x * x * x * x * x * y * y * z - std::sqrt(1614.990234375) * x * x * x * x * x * x * z * z * z - std::sqrt(90843.20068359375) * x * x * x * x * y * y * y * y * z - std::sqrt(40374.755859375) * x * x * x * x * y * y * z * z * z + std::sqrt(258.3984375) * x * x * x * x * z * z * z * z * z - std::sqrt(161499.0234375) * x * x * y * y * y * y * y * y * z + std::sqrt(1009368.896484375) * x * x * y * y * y * y * z * z * z - std::sqrt(9302.34375) * x * x * y * y * z * z * z * z * z + std::sqrt(2523.4222412109375) * y * y * y * y * y * y * y * y * z - std::sqrt(14534.912109375) * y * y * y * y * y * y * z * z * z + std::sqrt(258.3984375) * y * y * y * y * z * z * z * z * z) + e_2 * (std::sqrt(25839.84375) * x * x * x * x * x * x * z - std::sqrt(103359.375) * x * x * x * x * z * z * z - std::sqrt(5813964.84375) * x * x * y * y * y * y * z + std::sqrt(3720937.5) * x * x * y * y * z * z * z + std::sqrt(103359.375) * y * y * y * y * y * y * z - std::sqrt(103359.375) * y * y * y * y * z * z * z) + e_3 * (std::sqrt(232558.59375) * x * x * x * x * z - std::sqrt(8372109.375) * x * x * y * y * z + std::sqrt(232558.59375) * y * y * y * y * z);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        pc_56[k] = e_0 * (std::sqrt(645.99609375) * x * x * x * x * x * x * x * y * y * z * z + std::sqrt(645.99609375) * x * x * x * x * x * y * y * y * y * z * z - std::sqrt(10335.9375) * x * x * x * x * x * y * y * z * z * z * z - std::sqrt(645.99609375) * x * x * x * y * y * y * y * y * y * z * z + std::sqrt(1653.75) * x * x * x * y * y * z * z * z * z * z * z - std::sqrt(645.99609375) * x * y * y * y * y * y * y * y * y * z * z + std::sqrt(10335.9375) * x * y * y * y * y * y * y * z * z * z * z - std::sqrt(1653.75) * x * y * y * y * y * z * z * z * z * z * z) + e_1 * (std::sqrt(645.99609375) * x * x * x * x * x * x * x * y * y + std::sqrt(645.99609375) * x * x * x * x * x * x * x * z * z + std::sqrt(645.99609375) * x * x * x * x * x * y * y * y * y + std::sqrt(5813.96484375) * x * x * x * x * x * y * y * z * z - std::sqrt(10335.9375) * x * x * x * x * x * z * z * z * z - std::sqrt(645.99609375) * x * x * x * y * y * y * y * y * y - std::sqrt(16149.90234375) * x * x * x * y * y * y * y * z * z - std::sqrt(165375.0) * x * x * x * y * y * z * z * z * z + std::sqrt(1653.75) * x * x * x * z * z * z * z * z * z - std::sqrt(645.99609375) * x * y * y * y * y * y * y * y * y - std::sqrt(31653.80859375) * x * y * y * y * y * y * y * z * z + std::sqrt(837210.9375) * x * y * y * y * y * z * z * z * z - std::sqrt(14883.75) * x * y * y * z * z * z * z * z * z) + e_2 * (std::sqrt(645.99609375) * x * x * x * x * x * x * x + std::sqrt(145349.12109375) * x * x * x * x * x * y * y - std::sqrt(16149.90234375) * x * x * x * y * y * y * y - std::sqrt(1488375.0) * x * x * x * y * y * z * z - std::sqrt(165375.0) * x * x * x * z * z * z * z - std::sqrt(233204.58984375) * x * y * y * y * y * y * y + std::sqrt(1488375.0) * x * y * y * y * y * z * z + std::sqrt(1488375.0) * x * y * y * z * z * z * z) + e_3 * (std::sqrt(93023.4375) * x * x * x * x * x + std::sqrt(372093.75) * x * x * x * y * y - std::sqrt(1488375.0) * x * x * x * z * z - std::sqrt(4558148.4375) * x * y * y * y * y + std::sqrt(13395375.0) * x * y * y * z * z) + e_4 * (std::sqrt(372093.75) * x * x * x - std::sqrt(3348843.75) * x * y * y);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        pc_57[k] = e_0 * (-std::sqrt(20.1873779296875) * x * x * x * x * x * x * x * x * y * y * z - std::sqrt(143.5546875) * x * x * x * x * x * x * y * y * y * y * z + std::sqrt(2906.982421875) * x * x * x * x * x * x * y * y * z * z * z - std::sqrt(80.74951171875) * x * x * x * x * y * y * y * y * y * y * z + std::sqrt(8074.951171875) * x * x * x * x * y * y * y * y * z * z * z - std::sqrt(22790.7421875) * x * x * x * x * y * y * z * z * z * z * z + std::sqrt(322.998046875) * x * x * y * y * y * y * y * y * z * z * z - std::sqrt(10129.21875) * x * x * y * y * y * y * z * z * z * z * z + std::sqrt(3307.5) * x * x * y * y * z * z * z * z * z * z * z + std::sqrt(2.2430419921875) * y * y * y * y * y * y * y * y * y * y * z - std::sqrt(322.998046875) * y * y * y * y * y * y * y * y * z * z * z + std::sqrt(2532.3046875) * y * y * y * y * y * y * z * z * z * z * z - std::sqrt(367.5) * y * y * y * y * z * z * z * z * z * z * z) + e_1 * (-std::sqrt(20.1873779296875) * x * x * x * x * x * x * x * x * z - std::sqrt(1291.9921875) * x * x * x * x * x * x * y * y * z + std::sqrt(2906.982421875) * x * x * x * x * x * x * z * z * z - std::sqrt(2018.73779296875) * x * x * x * x * y * y * y * y * z - std::sqrt(26162.841796875) * x * x * x * x * y * y * z * z * z - std::sqrt(22790.7421875) * x * x * x * x * z * z * z * z * z - std::sqrt(39082.763671875) * x * x * y * y * y * y * z * z * z - std::sqrt(91162.96875) * x * x * y * y * z * z * z * z * z + std::sqrt(3307.5) * x * x * z * z * z * z * z * z * z + std::sqrt(181.6864013671875) * y * y * y * y * y * y * y * y * z + std::sqrt(322.998046875) * y * y * y * y * y * y * z * z * z + std::sqrt(63307.6171875) * y * y * y * y * z * z * z * z * z - std::sqrt(3307.5) * y * y * z * z * z * z * z * z * z) + e_2 * (-std::sqrt(418605.46875) * x * x * x * x * y * y * z - std::sqrt(186046.875) * x * x * x * x * z * z * z - std::sqrt(186046.875) * x * x * y * y * y * y * z - std::sqrt(6697687.5) * x * x * y * y * z * z * z + std::sqrt(46511.71875) * y * y * y * y * y * y * z + std::sqrt(1674421.875) * y * y * y * y * z * z * z) + e_3 * (-std::sqrt(418605.46875) * x * x * x * x * z - std::sqrt(22511671.875) * x * x * y * y * z - std::sqrt(2976750.0) * x * x * z * z * z + std::sqrt(4966417.96875) * y * y * y * y * z + std::sqrt(2976750.0) * y * y * z * z * z) + e_4 * (-std::sqrt(11907000.0) * x * x * z + std::sqrt(11907000.0) * y * y * z);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        pc_58[k] = e_0 * (-std::sqrt(215.33203125) * x * x * x * x * x * x * x * y * y * z * z - std::sqrt(1937.98828125) * x * x * x * x * x * y * y * y * y * z * z + std::sqrt(7751.953125) * x * x * x * x * x * y * y * z * z * z * z - std::sqrt(1937.98828125) * x * x * x * y * y * y * y * y * y * z * z + std::sqrt(31007.8125) * x * x * x * y * y * y * y * z * z * z * z - std::sqrt(19845.0) * x * x * x * y * y * z * z * z * z * z * z - std::sqrt(215.33203125) * x * y * y * y * y * y * y * y * y * z * z + std::sqrt(7751.953125) * x * y * y * y * y * y * y * z * z * z * z - std::sqrt(19845.0) * x * y * y * y * y * z * z * z * z * z * z + std::sqrt(2205.0) * x * y * y * z * z * z * z * z * z * z * z) + e_1 * (-std::sqrt(215.33203125) * x * x * x * x * x * x * x * y * y - std::sqrt(215.33203125) * x * x * x * x * x * x * x * z * z - std::sqrt(1937.98828125) * x * x * x * x * x * y * y * y * y - std::sqrt(1937.98828125) * x * x * x * x * x * y * y * z * z + std::sqrt(7751.953125) * x * x * x * x * x * z * z * z * z - std::sqrt(1937.98828125) * x * x * x * y * y * y * y * y * y - std::sqrt(1937.98828125) * x * x * x * y * y * y * y * z * z - std::sqrt(31007.8125) * x * x * x * y * y * z * z * z * z - std::sqrt(19845.0) * x * x * x * z * z * z * z * z * z - std::sqrt(215.33203125) * x * y * y * y * y * y * y * y * y - std::sqrt(215.33203125) * x * y * y * y * y * y * y * z * z - std::sqrt(69767.578125) * x * y * y * y * y * z * z * z * z + std::sqrt(2205.0) * x * y * y * z * z * z * z * z * z + std::sqrt(2205.0) * x * z * z * z * z * z * z * z * z) + e_2 * (-std::sqrt(215.33203125) * x * x * x * x * x * x * x - std::sqrt(94961.42578125) * x * x * x * x * x * y * y + std::sqrt(7751.953125) * x * x * x * x * x * z * z - std::sqrt(327520.01953125) * x * x * x * y * y * y * y - std::sqrt(775195.3125) * x * x * x * y * y * z * z - std::sqrt(496125.0) * x * x * x * z * z * z * z - std::sqrt(77734.86328125) * x * y * y * y * y * y * y - std::sqrt(937986.328125) * x * y * y * y * y * z * z - std::sqrt(496125.0) * x * y * y * z * z * z * z + std::sqrt(220500.0) * x * z * z * z * z * z * z) + e_3 * (-std::sqrt(31007.8125) * x * x * x * x * x - std::sqrt(6077531.25) * x * x * x * y * y - std::sqrt(1984500.0) * x * x * x * z * z - std::sqrt(5240320.3125) * x * y * y * y * y - std::sqrt(17860500.0) * x * y * y * z * z + std::sqrt(1984500.0) * x * z * z * z * z) + e_4 * (-std::sqrt(1984500.0) * x * x * x - std::sqrt(40186125.0) * x * y * y) + e_5 * (-std::sqrt(4465125.0) * x);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        pc_59[k] = e_0 * (std::sqrt(1.922607421875) * x * x * x * x * x * x * x * x * y * y * z + std::sqrt(30.76171875) * x * x * x * x * x * x * y * y * y * y * z - std::sqrt(492.1875) * x * x * x * x * x * x * y * y * z * z * z + std::sqrt(69.2138671875) * x * x * x * x * y * y * y * y * y * y * z - std::sqrt(4429.6875) * x * x * x * x * y * y * y * y * z * z * z + std::sqrt(6378.75) * x * x * x * x * y * y * z * z * z * z * z + std::sqrt(30.76171875) * x * x * y * y * y * y * y * y * y * y * z - std::sqrt(4429.6875) * x * x * y * y * y * y * y * y * z * z * z + std::sqrt(25515.0) * x * x * y * y * y * y * z * z * z * z * z - std::sqrt(5040.0) * x * x * y * y * z * z * z * z * z * z * z + std::sqrt(1.922607421875) * y * y * y * y * y * y * y * y * y * y * z - std::sqrt(492.1875) * y * y * y * y * y * y * y * y * z * z * z + std::sqrt(6378.75) * y * y * y * y * y * y * z * z * z * z * z - std::sqrt(5040.0) * y * y * y * y * z * z * z * z * z * z * z + std::sqrt(315.0) * y * y * z * z * z * z * z * z * z * z * z) + e_1 * (std::sqrt(1.922607421875) * x * x * x * x * x * x * x * x * z + std::sqrt(30.76171875) * x * x * x * x * x * x * y * y * z - std::sqrt(492.1875) * x * x * x * x * x * x * z * z * z + std::sqrt(69.2138671875) * x * x * x * x * y * y * y * y * z + std::sqrt(4429.6875) * x * x * x * x * y * y * z * z * z + std::sqrt(6378.75) * x * x * x * x * z * z * z * z * z + std::sqrt(30.76171875) * x * x * y * y * y * y * y * y * z + std::sqrt(39867.1875) * x * x * y * y * y * y * z * z * z - std::sqrt(2835.0) * x * x * y * y * z * z * z * z * z - std::sqrt(5040.0) * x * x * z * z * z * z * z * z * z + std::sqrt(1.922607421875) * y * y * y * y * y * y * y * y * z + std::sqrt(12304.6875) * y * y * y * y * y * y * z * z * z - std::sqrt(17718.75) * y * y * y * y * z * z * z * z * z + std::sqrt(20160.0) * y * y * z * z * z * z * z * z * z + std::sqrt(315.0) * z * z * z * z * z * z * z * z * z) + e_2 * (-std::sqrt(123.046875) * x * x * x * x * x * x * z + std::sqrt(27685.546875) * x * x * x * x * y * y * z + std::sqrt(70875.0) * x * x * x * x * z * z * z + std::sqrt(133998.046875) * x * x * y * y * y * y * z + std::sqrt(283500.0) * x * x * y * y * z * z * z - std::sqrt(283500.0) * x * x * z * z * z * z * z + std::sqrt(35560.546875) * y * y * y * y * y * y * z + std::sqrt(70875.0) * y * y * y * y * z * z * z + std::sqrt(1134000.0) * y * y * z * z * z * z * z + std::sqrt(126000.0) * z * z * z * z * z * z * z) + e_3 * (std::sqrt(70875.0) * x * x * x * x * z + std::sqrt(2551500.0) * x * x * y * y * z - std::sqrt(1134000.0) * x * x * z * z * z + std::sqrt(1771875.0) * y * y * y * y * z + std::sqrt(18144000.0) * y * y * z * z * z + std::sqrt(7087500.0) * z * z * z * z * z) + e_4 * (std::sqrt(31255875.0) * y * y * z + std::sqrt(55566000.0) * z * z * z) + e_5 * (std::sqrt(31255875.0) * z);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        pc_60[k] = e_0 * (std::sqrt(28.839111328125) * x * x * x * x * x * x * x * x * y * z * z + std::sqrt(461.42578125) * x * x * x * x * x * x * y * y * y * z * z - std::sqrt(1281.73828125) * x * x * x * x * x * x * y * z * z * z * z + std::sqrt(1038.2080078125) * x * x * x * x * y * y * y * y * y * z * z - std::sqrt(11535.64453125) * x * x * x * x * y * y * y * z * z * z * z + std::sqrt(4725.0) * x * x * x * x * y * z * z * z * z * z * z + std::sqrt(461.42578125) * x * x * y * y * y * y * y * y * y * z * z - std::sqrt(11535.64453125) * x * x * y * y * y * y * y * z * z * z * z + std::sqrt(18900.0) * x * x * y * y * y * z * z * z * z * z * z - std::sqrt(1181.25) * x * x * y * z * z * z * z * z * z * z * z + std::sqrt(28.839111328125) * y * y * y * y * y * y * y * y * y * z * z - std::sqrt(1281.73828125) * y * y * y * y * y * y * y * z * z * z * z + std::sqrt(4725.0) * y * y * y * y * y * z * z * z * z * z * z - std::sqrt(1181.25) * y * y * y * z * z * z * z * z * z * z * z + std::sqrt(21.0) * y * z * z * z * z * z * z * z * z * z * z) + e_1 * (std::sqrt(28.839111328125) * x * x * x * x * x * x * x * x * y + std::sqrt(461.42578125) * x * x * x * x * x * x * y * y * y + std::sqrt(1038.2080078125) * x * x * x * x * y * y * y * y * y + std::sqrt(29531.25) * x * x * x * x * y * z * z * z * z + std::sqrt(461.42578125) * x * x * y * y * y * y * y * y * y + std::sqrt(118125.0) * x * x * y * y * y * z * z * z * z - std::sqrt(18900.0) * x * x * y * z * z * z * z * z * z + std::sqrt(28.839111328125) * y * y * y * y * y * y * y * y * y + std::sqrt(29531.25) * y * y * y * y * y * z * z * z * z - std::sqrt(18900.0) * y * y * y * z * z * z * z * z * z + std::sqrt(4725.0) * y * z * z * z * z * z * z * z * z) + e_2 * (std::sqrt(11535.64453125) * x * x * x * x * x * x * y + std::sqrt(103820.80078125) * x * x * x * x * y * y * y + std::sqrt(265781.25) * x * x * x * x * y * z * z + std::sqrt(103820.80078125) * x * x * y * y * y * y * y + std::sqrt(1063125.0) * x * x * y * y * y * z * z + std::sqrt(11535.64453125) * y * y * y * y * y * y * y + std::sqrt(265781.25) * y * y * y * y * y * z * z + std::sqrt(472500.0) * y * z * z * z * z * z * z) + e_3 * (std::sqrt(1063125.0) * x * x * x * x * y + std::sqrt(4252500.0) * x * x * y * y * y + std::sqrt(4252500.0) * x * x * y * z * z + std::sqrt(1063125.0) * y * y * y * y * y + std::sqrt(4252500.0) * y * y * y * z * z + std::sqrt(11812500.0) * y * z * z * z * z) + e_4 * (std::sqrt(13023281.25) * x * x * y + std::sqrt(13023281.25) * y * y * y + std::sqrt(52093125.0) * y * z * z) + e_5 * (std::sqrt(18753525.0) * y);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        pc_61[k] = e_0 * (std::sqrt(1.922607421875) * x * x * x * x * x * x * x * x * x * y * z + std::sqrt(30.76171875) * x * x * x * x * x * x * x * y * y * y * z - std::sqrt(492.1875) * x * x * x * x * x * x * x * y * z * z * z + std::sqrt(69.2138671875) * x * x * x * x * x * y * y * y * y * y * z - std::sqrt(4429.6875) * x * x * x * x * x * y * y * y * z * z * z + std::sqrt(6378.75) * x * x * x * x * x * y * z * z * z * z * z + std::sqrt(30.76171875) * x * x * x * y * y * y * y * y * y * y * z - std::sqrt(4429.6875) * x * x * x * y * y * y * y * y * z * z * z + std::sqrt(25515.0) * x * x * x * y * y * y * z * z * z * z * z - std::sqrt(5040.0) * x * x * x * y * z * z * z * z * z * z * z + std::sqrt(1.922607421875) * x * y * y * y * y * y * y * y * y * y * z - std::sqrt(492.1875) * x * y * y * y * y * y * y * y * z * z * z + std::sqrt(6378.75) * x * y * y * y * y * y * z * z * z * z * z - std::sqrt(5040.0) * x * y * y * y * z * z * z * z * z * z * z + std::sqrt(315.0) * x * y * z * z * z * z * z * z * z * z * z) + e_1 * (std::sqrt(17718.75) * x * x * x * x * x * y * z * z * z + std::sqrt(70875.0) * x * x * x * y * y * y * z * z * z - std::sqrt(45360.0) * x * x * x * y * z * z * z * z * z + std::sqrt(17718.75) * x * y * y * y * y * y * z * z * z - std::sqrt(45360.0) * x * y * y * y * z * z * z * z * z + std::sqrt(45360.0) * x * y * z * z * z * z * z * z * z) + e_2 * (std::sqrt(39867.1875) * x * x * x * x * x * y * z + std::sqrt(159468.75) * x * x * x * y * y * y * z + std::sqrt(39867.1875) * x * y * y * y * y * y * z + std::sqrt(2551500.0) * x * y * z * z * z * z * z) + e_3 * (std::sqrt(1134000.0) * x * x * x * y * z + std::sqrt(1134000.0) * x * y * y * y * z + std::sqrt(28350000.0) * x * y * z * z * z) + e_4 * (std::sqrt(31255875.0) * x * y * z);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        pc_62[k] = e_0 * (-std::sqrt(53.8330078125) * x * x * x * x * x * x * x * x * y * z * z - std::sqrt(215.33203125) * x * x * x * x * x * x * y * y * y * z * z + std::sqrt(1937.98828125) * x * x * x * x * x * x * y * z * z * z * z + std::sqrt(1937.98828125) * x * x * x * x * y * y * y * z * z * z * z - std::sqrt(4961.25) * x * x * x * x * y * z * z * z * z * z * z + std::sqrt(215.33203125) * x * x * y * y * y * y * y * y * y * z * z - std::sqrt(1937.98828125) * x * x * y * y * y * y * y * z * z * z * z + std::sqrt(551.25) * x * x * y * z * z * z * z * z * z * z * z + std::sqrt(53.8330078125) * y * y * y * y * y * y * y * y * y * z * z - std::sqrt(1937.98828125) * y * y * y * y * y * y * y * z * z * z * z + std::sqrt(4961.25) * y * y * y * y * y * z * z * z * z * z * z - std::sqrt(551.25) * y * y * y * z * z * z * z * z * z * z * z) + e_1 * (-std::sqrt(53.8330078125) * x * x * x * x * x * x * x * x * y - std::sqrt(215.33203125) * x * x * x * x * x * x * y * y * y + std::sqrt(215.33203125) * x * x * x * x * x * x * y * z * z + std::sqrt(1937.98828125) * x * x * x * x * y * y * y * z * z - std::sqrt(69767.578125) * x * x * x * x * y * z * z * z * z + std::sqrt(215.33203125) * x * x * y * y * y * y * y * y * y + std::sqrt(1937.98828125) * x * x * y * y * y * y * y * z * z - std::sqrt(31007.8125) * x * x * y * y * y * z * z * z * z + std::sqrt(55125.0) * x * x * y * z * z * z * z * z * z + std::sqrt(53.8330078125) * y * y * y * y * y * y * y * y * y + std::sqrt(215.33203125) * y * y * y * y * y * y * y * z * z + std::sqrt(7751.953125) * y * y * y * y * y * z * z * z * z + std::sqrt(2205.0) * y * y * y * z * z * z * z * z * z - std::sqrt(2205.0) * y * z * z * z * z * z * z * z * z) + e_2 * (-std::sqrt(13781.25) * x * x * x * x * x * x * y - std::sqrt(7751.953125) * x * x * x * x * y * y * y - std::sqrt(379845.703125) * x * x * x * x * y * z * z + std::sqrt(31007.8125) * x * x * y * y * y * y * y - std::sqrt(31007.8125) * x * x * y * y * y * z * z + std::sqrt(496125.0) * x * x * y * z * z * z * z + std::sqrt(21533.203125) * y * y * y * y * y * y * y + std::sqrt(193798.828125) * y * y * y * y * y * z * z + std::sqrt(496125.0) * y * y * y * z * z * z * z - std::sqrt(220500.0) * y * z * z * z * z * z * z) + e_3 * (-std::sqrt(775195.3125) * x * x * x * x * y + std::sqrt(124031.25) * x * x * y * y * y + std::sqrt(1519382.8125) * y * y * y * y * y + std::sqrt(7938000.0) * y * y * y * z * z - std::sqrt(1984500.0) * y * z * z * z * z) + e_4 * (-std::sqrt(1116281.25) * x * x * y + std::sqrt(15007781.25) * y * y * y) + e_5 * (std::sqrt(4465125.0) * y);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        pc_63[k] = e_0 * (-std::sqrt(2.2430419921875) * x * x * x * x * x * x * x * x * x * y * z + std::sqrt(322.998046875) * x * x * x * x * x * x * x * y * z * z * z + std::sqrt(80.74951171875) * x * x * x * x * x * y * y * y * y * y * z - std::sqrt(322.998046875) * x * x * x * x * x * y * y * y * z * z * z - std::sqrt(2532.3046875) * x * x * x * x * x * y * z * z * z * z * z + std::sqrt(143.5546875) * x * x * x * y * y * y * y * y * y * y * z - std::sqrt(8074.951171875) * x * x * x * y * y * y * y * y * z * z * z + std::sqrt(10129.21875) * x * x * x * y * y * y * z * z * z * z * z + std::sqrt(367.5) * x * x * x * y * z * z * z * z * z * z * z + std::sqrt(20.1873779296875) * x * y * y * y * y * y * y * y * y * y * z - std::sqrt(2906.982421875) * x * y * y * y * y * y * y * y * z * z * z + std::sqrt(22790.7421875) * x * y * y * y * y * y * z * z * z * z * z - std::sqrt(3307.5) * x * y * y * y * z * z * z * z * z * z * z) + e_1 * (std::sqrt(1291.9921875) * x * x * x * x * x * y * y * y * z - std::sqrt(32299.8046875) * x * x * x * x * x * y * z * z * z + std::sqrt(5167.96875) * x * x * x * y * y * y * y * y * z - std::sqrt(5167.96875) * x * x * x * y * y * y * z * z * z + std::sqrt(40516.875) * x * x * x * y * z * z * z * z * z + std::sqrt(1291.9921875) * x * y * y * y * y * y * y * y * z + std::sqrt(11627.9296875) * x * y * y * y * y * y * z * z * z + std::sqrt(364651.875) * x * y * y * y * z * z * z * z * z - std::sqrt(13230.0) * x * y * z * z * z * z * z * z * z) + e_2 * (-std::sqrt(46511.71875) * x * x * x * x * x * y * z + std::sqrt(186046.875) * x * x * x * y * y * y * z + std::sqrt(418605.46875) * x * y * y * y * y * y * z + std::sqrt(11907000.0) * x * y * y * y * z * z * z) + e_3 * (-std::sqrt(82687.5) * x * x * x * y * z + std::sqrt(36465187.5) * x * y * y * y * z + std::sqrt(11907000.0) * x * y * z * z * z) + e_4 * (std::sqrt(47628000.0) * x * y * z);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        pc_64[k] = e_0 * (std::sqrt(40.374755859375) * x * x * x * x * x * x * x * x * y * z * z - std::sqrt(645.99609375) * x * x * x * x * x * x * y * y * y * z * z - std::sqrt(645.99609375) * x * x * x * x * x * x * y * z * z * z * z - std::sqrt(4037.4755859375) * x * x * x * x * y * y * y * y * y * z * z + std::sqrt(16149.90234375) * x * x * x * x * y * y * y * z * z * z * z + std::sqrt(103.359375) * x * x * x * x * y * z * z * z * z * z * z - std::sqrt(645.99609375) * x * x * y * y * y * y * y * y * y * z * z + std::sqrt(16149.90234375) * x * x * y * y * y * y * y * z * z * z * z - std::sqrt(3720.9375) * x * x * y * y * y * z * z * z * z * z * z + std::sqrt(40.374755859375) * y * y * y * y * y * y * y * y * y * z * z - std::sqrt(645.99609375) * y * y * y * y * y * y * y * z * z * z * z + std::sqrt(103.359375) * y * y * y * y * y * z * z * z * z * z * z) + e_1 * (std::sqrt(40.374755859375) * x * x * x * x * x * x * x * x * y - std::sqrt(645.99609375) * x * x * x * x * x * x * y * y * y - std::sqrt(2583.984375) * x * x * x * x * x * x * y * z * z - std::sqrt(4037.4755859375) * x * x * x * x * y * y * y * y * y - std::sqrt(64599.609375) * x * x * x * x * y * y * y * z * z + std::sqrt(23255.859375) * x * x * x * x * y * z * z * z * z - std::sqrt(645.99609375) * x * x * y * y * y * y * y * y * y - std::sqrt(23255.859375) * x * x * y * y * y * y * y * z * z + std::sqrt(1250648.4375) * x * x * y * y * y * z * z * z * z - std::sqrt(14883.75) * x * x * y * z * z * z * z * z * z + std::sqrt(40.374755859375) * y * y * y * y * y * y * y * y * y + std::sqrt(2583.984375) * y * y * y * y * y * y * y * z * z - std::sqrt(64599.609375) * y * y * y * y * y * z * z * z * z + std::sqrt(1653.75) * y * y * y * z * z * z * z * z * z) + e_2 * (std::sqrt(645.99609375) * x * x * x * x * x * x * y - std::sqrt(403747.55859375) * x * x * x * x * y * y * y - std::sqrt(93023.4375) * x * x * x * x * y * z * z - std::sqrt(284884.27734375) * x * x * y * y * y * y * y + std::sqrt(3348843.75) * x * x * y * y * y * z * z + std::sqrt(1488375.0) * x * x * y * z * z * z * z + std::sqrt(16149.90234375) * y * y * y * y * y * y * y - std::sqrt(93023.4375) * y * y * y * y * y * z * z - std::sqrt(165375.0) * y * y * y * z * z * z * z) + e_3 * (-std::sqrt(372093.75) * x * x * x * x * y - std::sqrt(5953500.0) * x * x * y * y * y + std::sqrt(13395375.0) * x * x * y * z * z + std::sqrt(372093.75) * y * y * y * y * y - std::sqrt(1488375.0) * y * y * y * z * z) + e_4 * (-std::sqrt(3348843.75) * x * x * y + std::sqrt(372093.75) * y * y * y);
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

        pc_65[k] = e_0 * (std::sqrt(4.0374755859375) * x * x * x * x * x * x * x * x * x * y * z - std::sqrt(258.3984375) * x * x * x * x * x * x * x * y * y * y * z - std::sqrt(64.599609375) * x * x * x * x * x * x * x * y * z * z * z - std::sqrt(791.34521484375) * x * x * x * x * x * y * y * y * y * y * z + std::sqrt(5232.568359375) * x * x * x * x * x * y * y * y * z * z * z + std::sqrt(10.3359375) * x * x * x * x * x * y * z * z * z * z * z + std::sqrt(1614.990234375) * x * x * x * y * y * y * y * y * z * z * z - std::sqrt(1033.59375) * x * x * x * y * y * y * z * z * z * z * z + std::sqrt(100.9368896484375) * x * y * y * y * y * y * y * y * y * y * z - std::sqrt(1614.990234375) * x * y * y * y * y * y * y * y * z * z * z + std::sqrt(258.3984375) * x * y * y * y * y * y * z * z * z * z * z) + e_1 * (-std::sqrt(161499.0234375) * x * x * x * x * x * y * y * y * z + std::sqrt(6459.9609375) * x * x * x * x * x * y * z * z * z - std::sqrt(25839.84375) * x * x * x * y * y * y * y * y * z + std::sqrt(645996.09375) * x * x * x * y * y * y * z * z * z - std::sqrt(4134.375) * x * x * x * y * z * z * z * z * z + std::sqrt(58139.6484375) * x * y * y * y * y * y * y * y * z - std::sqrt(316538.0859375) * x * y * y * y * y * y * z * z * z + std::sqrt(4134.375) * x * y * y * y * z * z * z * z * z) + e_2 * (-std::sqrt(232558.59375) * x * x * x * x * x * y * z - std::sqrt(2583984.375) * x * x * x * y * y * y * z + std::sqrt(1653750.0) * x * x * x * y * z * z * z + std::sqrt(2093027.34375) * x * y * y * y * y * y * z - std::sqrt(1653750.0) * x * y * y * y * z * z * z) + e_3 * (-std::sqrt(3720937.5) * x * x * x * y * z + std::sqrt(3720937.5) * x * y * y * y * z);
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

        pc_66[k] = e_0 * (-std::sqrt(1.201629638671875) * x * x * x * x * x * x * x * x * x * x * y - std::sqrt(1.201629638671875) * x * x * x * x * x * x * x * x * y * y * y + std::sqrt(389.3280029296875) * x * x * x * x * x * x * x * x * y * z * z + std::sqrt(9.4207763671875) * x * x * x * x * x * x * y * y * y * y * y - std::sqrt(692.138671875) * x * x * x * x * x * x * y * z * z * z * z + std::sqrt(23.2635498046875) * x * x * x * x * y * y * y * y * y * y * y - std::sqrt(3052.33154296875) * x * x * x * x * y * y * y * y * y * z * z + std::sqrt(692.138671875) * x * x * x * x * y * y * y * z * z * z * z + std::sqrt(12.3046875) * x * x * x * x * y * z * z * z * z * z * z + std::sqrt(2.355194091796875) * x * x * y * y * y * y * y * y * y * y * y - std::sqrt(996.6796875) * x * x * y * y * y * y * y * y * y * z * z + std::sqrt(2242.529296875) * x * x * y * y * y * y * y * z * z * z * z - std::sqrt(49.21875) * x * x * y * y * y * z * z * z * z * z * z - std::sqrt(0.048065185546875) * y * y * y * y * y * y * y * y * y * y * y + std::sqrt(15.5731201171875) * y * y * y * y * y * y * y * y * y * z * z - std::sqrt(27.685546875) * y * y * y * y * y * y * y * z * z * z * z + std::sqrt(0.4921875) * y * y * y * y * y * z * z * z * z * z * z) + e_1 * (-std::sqrt(1081.4666748046875) * x * x * x * x * x * x * x * x * y + std::sqrt(155731.201171875) * x * x * x * x * x * x * y * z * z + std::sqrt(8478.69873046875) * x * x * x * x * y * y * y * y * y - std::sqrt(155731.201171875) * x * x * x * x * y * y * y * z * z - std::sqrt(69213.8671875) * x * x * x * x * y * z * z * z * z + std::sqrt(2768.5546875) * x * x * y * y * y * y * y * y * y - std::sqrt(504569.091796875) * x * x * y * y * y * y * y * z * z + std::sqrt(276855.46875) * x * x * y * y * y * z * z * z * z - std::sqrt(43.2586669921875) * y * y * y * y * y * y * y * y * y + std::sqrt(6229.248046875) * y * y * y * y * y * y * y * z * z - std::sqrt(2768.5546875) * y * y * y * y * y * z * z * z * z) + e_2 * (-std::sqrt(69213.8671875) * x * x * x * x * x * x * y + std::sqrt(69213.8671875) * x * x * x * x * y * y * y + std::sqrt(2491699.21875) * x * x * x * x * y * z * z + std::sqrt(224252.9296875) * x * x * y * y * y * y * y - std::sqrt(9966796.875) * x * x * y * y * y * z * z - std::sqrt(2768.5546875) * y * y * y * y * y * y * y + std::sqrt(99667.96875) * y * y * y * y * y * z * z) + e_3 * (-std::sqrt(276855.46875) * x * x * x * x * y + std::sqrt(1107421.875) * x * x * y * y * y - std::sqrt(11074.21875) * y * y * y * y * y);
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

        pc_67[k] = e_0 * (-std::sqrt(7.6904296875) * x * x * x * x * x * x * x * x * x * y * z - std::sqrt(30.76171875) * x * x * x * x * x * x * x * y * y * y * z + std::sqrt(2491.69921875) * x * x * x * x * x * x * x * y * z * z * z + std::sqrt(2491.69921875) * x * x * x * x * x * y * y * y * z * z * z - std::sqrt(4429.6875) * x * x * x * x * x * y * z * z * z * z * z + std::sqrt(30.76171875) * x * x * x * y * y * y * y * y * y * y * z - std::sqrt(2491.69921875) * x * x * x * y * y * y * y * y * z * z * z + std::sqrt(78.75) * x * x * x * y * z * z * z * z * z * z * z + std::sqrt(7.6904296875) * x * y * y * y * y * y * y * y * y * y * z - std::sqrt(2491.69921875) * x * y * y * y * y * y * y * y * z * z * z + std::sqrt(4429.6875) * x * y * y * y * y * y * z * z * z * z * z - std::sqrt(78.75) * x * y * y * y * z * z * z * z * z * z * z) + e_1 * (std::sqrt(1107.421875) * x * x * x * x * x * x * x * y * z + std::sqrt(1107.421875) * x * x * x * x * x * y * y * y * z + std::sqrt(283500.0) * x * x * x * x * x * y * z * z * z - std::sqrt(1107.421875) * x * x * x * y * y * y * y * y * z - std::sqrt(229635.0) * x * x * x * y * z * z * z * z * z - std::sqrt(1107.421875) * x * y * y * y * y * y * y * y * z - std::sqrt(283500.0) * x * y * y * y * y * y * z * z * z + std::sqrt(229635.0) * x * y * y * y * z * z * z * z * z) + e_2 * (std::sqrt(1435218.75) * x * x * x * x * x * y * z + std::sqrt(70875.0) * x * x * x * y * z * z * z - std::sqrt(1435218.75) * x * y * y * y * y * y * z - std::sqrt(70875.0) * x * y * y * y * z * z * z) + e_3 * (std::sqrt(18144000.0) * x * x * x * y * z - std::sqrt(18144000.0) * x * y * y * y * z);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        pc_68[k] = e_0 * (std::sqrt(0.240325927734375) * x * x * x * x * x * x * x * x * x * x * y + std::sqrt(3.231048583984375) * x * x * x * x * x * x * x * x * y * y * y - std::sqrt(162.4603271484375) * x * x * x * x * x * x * x * x * y * z * z + std::sqrt(5.2337646484375) * x * x * x * x * x * x * y * y * y * y * y - std::sqrt(1155.2734375) * x * x * x * x * x * x * y * y * y * z * z + std::sqrt(6782.958984375) * x * x * x * x * x * x * y * z * z * z * z + std::sqrt(0.9613037109375) * x * x * x * x * y * y * y * y * y * y * y - std::sqrt(649.84130859375) * x * x * x * x * y * y * y * y * y * z * z + std::sqrt(18841.552734375) * x * x * x * x * y * y * y * z * z * z * z - std::sqrt(9157.1484375) * x * x * x * x * y * z * z * z * z * z * z - std::sqrt(0.026702880859375) * x * x * y * y * y * y * y * y * y * y * y + std::sqrt(753.662109375) * x * x * y * y * y * y * y * z * z * z * z - std::sqrt(4069.84375) * x * x * y * y * y * z * z * z * z * z * z + std::sqrt(157.5) * x * x * y * z * z * z * z * z * z * z * z - std::sqrt(0.026702880859375) * y * y * y * y * y * y * y * y * y * y * y + std::sqrt(18.0511474609375) * y * y * y * y * y * y * y * y * y * z * z - std::sqrt(753.662109375) * y * y * y * y * y * y * y * z * z * z * z + std::sqrt(1017.4609375) * y * y * y * y * y * z * z * z * z * z * z - std::sqrt(17.5) * y * y * y * z * z * z * z * z * z * z * z) + e_1 * (std::sqrt(216.2933349609375) * x * x * x * x * x * x * x * x * y + std::sqrt(1538.0859375) * x * x * x * x * x * x * y * y * y + std::sqrt(1245.849609375) * x * x * x * x * x * x * y * z * z + std::sqrt(865.17333984375) * x * x * x * x * y * y * y * y * y + std::sqrt(3460.693359375) * x * x * x * x * y * y * y * z * z + std::sqrt(44850.5859375) * x * x * x * x * y * z * z * z * z + std::sqrt(138.427734375) * x * x * y * y * y * y * y * z * z + std::sqrt(19933.59375) * x * x * y * y * y * z * z * z * z - std::sqrt(171517.5) * x * x * y * z * z * z * z * z * z - std::sqrt(24.0325927734375) * y * y * y * y * y * y * y * y * y - std::sqrt(138.427734375) * y * y * y * y * y * y * y * z * z - std::sqrt(4983.3984375) * y * y * y * y * y * z * z * z * z + std::sqrt(19057.5) * y * y * y * z * z * z * z * z * z) + e_2 * (std::sqrt(79734.375) * x * x * x * x * x * x * y + std::sqrt(221484.375) * x * x * x * x * y * y * y + std::sqrt(976746.09375) * x * x * x * x * y * z * z + std::sqrt(8859.375) * x * x * y * y * y * y * y + std::sqrt(434109.375) * x * x * y * y * y * z * z - std::sqrt(5103000.0) * x * x * y * z * z * z * z - std::sqrt(8859.375) * y * y * y * y * y * y * y - std::sqrt(108527.34375) * y * y * y * y * y * z * z + std::sqrt(567000.0) * y * y * y * z * z * z * z) + e_3 * (std::sqrt(4892589.84375) * x * x * x * x * y + std::sqrt(2174484.375) * x * x * y * y * y - std::sqrt(3543750.0) * x * x * y * z * z - std::sqrt(543621.09375) * y * y * y * y * y + std::sqrt(393750.0) * y * y * y * z * z) + e_4 * (std::sqrt(15627937.5) * x * x * y - std::sqrt(1736437.5) * y * y * y);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        pc_69[k] = e_0 * (std::sqrt(2.5634765625) * x * x * x * x * x * x * x * x * x * y * z + std::sqrt(41.015625) * x * x * x * x * x * x * x * y * y * y * z - std::sqrt(1025.390625) * x * x * x * x * x * x * x * y * z * z * z + std::sqrt(92.28515625) * x * x * x * x * x * y * y * y * y * y * z - std::sqrt(9228.515625) * x * x * x * x * x * y * y * y * z * z * z + std::sqrt(9228.515625) * x * x * x * x * x * y * z * z * z * z * z + std::sqrt(41.015625) * x * x * x * y * y * y * y * y * y * y * z - std::sqrt(9228.515625) * x * x * x * y * y * y * y * y * z * z * z + std::sqrt(36914.0625) * x * x * x * y * y * y * z * z * z * z * z - std::sqrt(6720.0) * x * x * x * y * z * z * z * z * z * z * z + std::sqrt(2.5634765625) * x * y * y * y * y * y * y * y * y * y * z - std::sqrt(1025.390625) * x * y * y * y * y * y * y * y * z * z * z + std::sqrt(9228.515625) * x * y * y * y * y * y * z * z * z * z * z - std::sqrt(6720.0) * x * y * y * y * z * z * z * z * z * z * z + std::sqrt(105.0) * x * y * z * z * z * z * z * z * z * z * z) + e_1 * (-std::sqrt(369.140625) * x * x * x * x * x * x * x * y * z - std::sqrt(3322.265625) * x * x * x * x * x * y * y * y * z - std::sqrt(3322.265625) * x * x * x * y * y * y * y * y * z - std::sqrt(34020.0) * x * x * x * y * z * z * z * z * z - std::sqrt(369.140625) * x * y * y * y * y * y * y * y * z - std::sqrt(34020.0) * x * y * y * y * z * z * z * z * z - std::sqrt(15120.0) * x * y * z * z * z * z * z * z * z) + e_2 * (-std::sqrt(83056.640625) * x * x * x * x * x * y * z - std::sqrt(332226.5625) * x * x * x * y * y * y * z - std::sqrt(850500.0) * x * x * x * y * z * z * z - std::sqrt(83056.640625) * x * y * y * y * y * y * z - std::sqrt(850500.0) * x * y * y * y * z * z * z - std::sqrt(3402000.0) * x * y * z * z * z * z * z) + e_3 * (-std::sqrt(6048000.0) * x * x * x * y * z - std::sqrt(6048000.0) * x * y * y * y * z - std::sqrt(63882000.0) * x * y * z * z * z) + e_4 * (-std::sqrt(93767625.0) * x * y * z);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        pc_70[k] = e_0 * (-std::sqrt(0.02288818359375) * x * x * x * x * x * x * x * x * x * x * y - std::sqrt(0.57220458984375) * x * x * x * x * x * x * x * x * y * y * y + std::sqrt(20.599365234375) * x * x * x * x * x * x * x * x * y * z * z - std::sqrt(2.288818359375) * x * x * x * x * x * x * y * y * y * y * y + std::sqrt(329.58984375) * x * x * x * x * x * x * y * y * y * z * z - std::sqrt(1407.71484375) * x * x * x * x * x * x * y * z * z * z * z - std::sqrt(2.288818359375) * x * x * x * x * y * y * y * y * y * y * y + std::sqrt(741.5771484375) * x * x * x * x * y * y * y * y * y * z * z - std::sqrt(12669.43359375) * x * x * x * x * y * y * y * z * z * z * z + std::sqrt(4335.0) * x * x * x * x * y * z * z * z * z * z * z - std::sqrt(0.57220458984375) * x * x * y * y * y * y * y * y * y * y * y + std::sqrt(329.58984375) * x * x * y * y * y * y * y * y * y * z * z - std::sqrt(12669.43359375) * x * x * y * y * y * y * y * z * z * z * z + std::sqrt(17340.0) * x * x * y * y * y * z * z * z * z * z * z - std::sqrt(1215.0) * x * x * y * z * z * z * z * z * z * z * z - std::sqrt(0.02288818359375) * y * y * y * y * y * y * y * y * y * y * y + std::sqrt(20.599365234375) * y * y * y * y * y * y * y * y * y * z * z - std::sqrt(1407.71484375) * y * y * y * y * y * y * y * z * z * z * z + std::sqrt(4335.0) * y * y * y * y * y * z * z * z * z * z * z - std::sqrt(1215.0) * y * y * y * z * z * z * z * z * z * z * z + std::sqrt(15.0) * y * z * z * z * z * z * z * z * z * z * z) + e_1 * (-std::sqrt(20.599365234375) * x * x * x * x * x * x * x * x * y - std::sqrt(329.58984375) * x * x * x * x * x * x * y * y * y - std::sqrt(1898.4375) * x * x * x * x * x * x * y * z * z - std::sqrt(741.5771484375) * x * x * x * x * y * y * y * y * y - std::sqrt(17085.9375) * x * x * x * x * y * y * y * z * z + std::sqrt(7593.75) * x * x * x * x * y * z * z * z * z - std::sqrt(329.58984375) * x * x * y * y * y * y * y * y * y - std::sqrt(17085.9375) * x * x * y * y * y * y * y * z * z + std::sqrt(30375.0) * x * x * y * y * y * z * z * z * z - std::sqrt(34560.0) * x * x * y * z * z * z * z * z * z - std::sqrt(20.599365234375) * y * y * y * y * y * y * y * y * y - std::sqrt(1898.4375) * y * y * y * y * y * y * y * z * z + std::sqrt(7593.75) * y * y * y * y * y * z * z * z * z - std::sqrt(34560.0) * y * y * y * z * z * z * z * z * z + std::sqrt(1215.0) * y * z * z * z * z * z * z * z * z) + e_2 * (-std::sqrt(12669.43359375) * x * x * x * x * x * x * y - std::sqrt(114024.90234375) * x * x * x * x * y * y * y - std::sqrt(68343.75) * x * x * x * x * y * z * z - std::sqrt(114024.90234375) * x * x * y * y * y * y * y - std::sqrt(273375.0) * x * x * y * y * y * z * z - std::sqrt(759375.0) * x * x * y * z * z * z * z - std::sqrt(12669.43359375) * y * y * y * y * y * y * y - std::sqrt(68343.75) * y * y * y * y * y * z * z - std::sqrt(759375.0) * y * y * y * z * z * z * z + std::sqrt(13500.0) * y * z * z * z * z * z * z) + e_3 * (-std::sqrt(975375.0) * x * x * x * x * y - std::sqrt(3901500.0) * x * x * y * y * y - std::sqrt(7776000.0) * x * x * y * z * z - std::sqrt(975375.0) * y * y * y * y * y - std::sqrt(7776000.0) * y * y * y * z * z - std::sqrt(337500.0) * y * z * z * z * z) + e_4 * (-std::sqrt(13395375.0) * x * x * y - std::sqrt(13395375.0) * y * y * y - std::sqrt(13395375.0) * y * z * z) + e_5 * (-std::sqrt(13395375.0) * y);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        pc_71[k] = e_0 * (-0.5859375 * x * x * x * x * x * x * x * x * x * x * z - 2.9296875 * x * x * x * x * x * x * x * x * y * y * z + 12.109375 * x * x * x * x * x * x * x * x * z * z * z - 5.859375 * x * x * x * x * x * x * y * y * y * y * z + 48.4375 * x * x * x * x * x * x * y * y * z * z * z - 42.5 * x * x * x * x * x * x * z * z * z * z * z - 5.859375 * x * x * x * x * y * y * y * y * y * y * z + 72.65625 * x * x * x * x * y * y * y * y * z * z * z - 127.5 * x * x * x * x * y * y * z * z * z * z * z + 45.0 * x * x * x * x * z * z * z * z * z * z * z - 2.9296875 * x * x * y * y * y * y * y * y * y * y * z + 48.4375 * x * x * y * y * y * y * y * y * z * z * z - 127.5 * x * x * y * y * y * y * z * z * z * z * z + 90.0 * x * x * y * y * z * z * z * z * z * z * z - 12.5 * x * x * z * z * z * z * z * z * z * z * z - 0.5859375 * y * y * y * y * y * y * y * y * y * y * z + 12.109375 * y * y * y * y * y * y * y * y * z * z * z - 42.5 * y * y * y * y * y * y * z * z * z * z * z + 45.0 * y * y * y * y * z * z * z * z * z * z * z - 12.5 * y * y * z * z * z * z * z * z * z * z * z + z * z * z * z * z * z * z * z * z * z * z) + e_1 * (7.03125 * x * x * x * x * x * x * x * x * z + 28.125 * x * x * x * x * x * x * y * y * z - 37.5 * x * x * x * x * x * x * z * z * z + 42.1875 * x * x * x * x * y * y * y * y * z - 112.5 * x * x * x * x * y * y * z * z * z + 180.0 * x * x * x * x * z * z * z * z * z + 28.125 * x * x * y * y * y * y * y * y * z - 112.5 * x * x * y * y * y * y * z * z * z + 360.0 * x * x * y * y * z * z * z * z * z - 90.0 * x * x * z * z * z * z * z * z * z + 7.03125 * y * y * y * y * y * y * y * y * z - 37.5 * y * y * y * y * y * y * z * z * z + 180.0 * y * y * y * y * z * z * z * z * z - 90.0 * y * y * z * z * z * z * z * z * z + 30.0 * z * z * z * z * z * z * z * z * z) + e_2 * (56.25 * x * x * x * x * x * x * z + 168.75 * x * x * x * x * y * y * z + 562.5 * x * x * x * x * z * z * z + 168.75 * x * x * y * y * y * y * z + 1125.0 * x * x * y * y * z * z * z - 225.0 * x * x * z * z * z * z * z + 56.25 * y * y * y * y * y * y * z + 562.5 * y * y * y * y * z * z * z - 225.0 * y * y * z * z * z * z * z + 450.0 * z * z * z * z * z * z * z) + e_3 * (900.0 * x * x * x * x * z + 1800.0 * x * x * y * y * z + 750.0 * x * x * z * z * z + 900.0 * y * y * y * y * z + 750.0 * y * y * z * z * z + 3000.0 * z * z * z * z * z) + e_4 * (2362.5 * x * x * z + 2362.5 * y * y * z + 7875.0 * z * z * z) + e_5 * (5670.0 * z);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        pc_72[k] = e_0 * (-std::sqrt(0.02288818359375) * x * x * x * x * x * x * x * x * x * x * x - std::sqrt(0.57220458984375) * x * x * x * x * x * x * x * x * x * y * y + std::sqrt(20.599365234375) * x * x * x * x * x * x * x * x * x * z * z - std::sqrt(2.288818359375) * x * x * x * x * x * x * x * y * y * y * y + std::sqrt(329.58984375) * x * x * x * x * x * x * x * y * y * z * z - std::sqrt(1407.71484375) * x * x * x * x * x * x * x * z * z * z * z - std::sqrt(2.288818359375) * x * x * x * x * x * y * y * y * y * y * y + std::sqrt(741.5771484375) * x * x * x * x * x * y * y * y * y * z * z - std::sqrt(12669.43359375) * x * x * x * x * x * y * y * z * z * z * z + std::sqrt(4335.0) * x * x * x * x * x * z * z * z * z * z * z - std::sqrt(0.57220458984375) * x * x * x * y * y * y * y * y * y * y * y + std::sqrt(329.58984375) * x * x * x * y * y * y * y * y * y * z * z - std::sqrt(12669.43359375) * x * x * x * y * y * y * y * z * z * z * z + std::sqrt(17340.0) * x * x * x * y * y * z * z * z * z * z * z - std::sqrt(1215.0) * x * x * x * z * z * z * z * z * z * z * z - std::sqrt(0.02288818359375) * x * y * y * y * y * y * y * y * y * y * y + std::sqrt(20.599365234375) * x * y * y * y * y * y * y * y * y * z * z - std::sqrt(1407.71484375) * x * y * y * y * y * y * y * z * z * z * z + std::sqrt(4335.0) * x * y * y * y * y * z * z * z * z * z * z - std::sqrt(1215.0) * x * y * y * z * z * z * z * z * z * z * z + std::sqrt(15.0) * x * z * z * z * z * z * z * z * z * z * z) + e_1 * (-std::sqrt(20.599365234375) * x * x * x * x * x * x * x * x * x - std::sqrt(329.58984375) * x * x * x * x * x * x * x * y * y - std::sqrt(1898.4375) * x * x * x * x * x * x * x * z * z - std::sqrt(741.5771484375) * x * x * x * x * x * y * y * y * y - std::sqrt(17085.9375) * x * x * x * x * x * y * y * z * z + std::sqrt(7593.75) * x * x * x * x * x * z * z * z * z - std::sqrt(329.58984375) * x * x * x * y * y * y * y * y * y - std::sqrt(17085.9375) * x * x * x * y * y * y * y * z * z + std::sqrt(30375.0) * x * x * x * y * y * z * z * z * z - std::sqrt(34560.0) * x * x * x * z * z * z * z * z * z - std::sqrt(20.599365234375) * x * y * y * y * y * y * y * y * y - std::sqrt(1898.4375) * x * y * y * y * y * y * y * z * z + std::sqrt(7593.75) * x * y * y * y * y * z * z * z * z - std::sqrt(34560.0) * x * y * y * z * z * z * z * z * z + std::sqrt(1215.0) * x * z * z * z * z * z * z * z * z) + e_2 * (-std::sqrt(12669.43359375) * x * x * x * x * x * x * x - std::sqrt(114024.90234375) * x * x * x * x * x * y * y - std::sqrt(68343.75) * x * x * x * x * x * z * z - std::sqrt(114024.90234375) * x * x * x * y * y * y * y - std::sqrt(273375.0) * x * x * x * y * y * z * z - std::sqrt(759375.0) * x * x * x * z * z * z * z - std::sqrt(12669.43359375) * x * y * y * y * y * y * y - std::sqrt(68343.75) * x * y * y * y * y * z * z - std::sqrt(759375.0) * x * y * y * z * z * z * z + std::sqrt(13500.0) * x * z * z * z * z * z * z) + e_3 * (-std::sqrt(975375.0) * x * x * x * x * x - std::sqrt(3901500.0) * x * x * x * y * y - std::sqrt(7776000.0) * x * x * x * z * z - std::sqrt(975375.0) * x * y * y * y * y - std::sqrt(7776000.0) * x * y * y * z * z - std::sqrt(337500.0) * x * z * z * z * z) + e_4 * (-std::sqrt(13395375.0) * x * x * x - std::sqrt(13395375.0) * x * y * y - std::sqrt(13395375.0) * x * z * z) + e_5 * (-std::sqrt(13395375.0) * x);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        pc_73[k] = e_0 * (std::sqrt(0.640869140625) * x * x * x * x * x * x * x * x * x * x * z + std::sqrt(5.767822265625) * x * x * x * x * x * x * x * x * y * y * z - std::sqrt(256.34765625) * x * x * x * x * x * x * x * x * z * z * z + std::sqrt(2.5634765625) * x * x * x * x * x * x * y * y * y * y * z - std::sqrt(1025.390625) * x * x * x * x * x * x * y * y * z * z * z + std::sqrt(2307.12890625) * x * x * x * x * x * x * z * z * z * z * z - std::sqrt(2.5634765625) * x * x * x * x * y * y * y * y * y * y * z + std::sqrt(2307.12890625) * x * x * x * x * y * y * z * z * z * z * z - std::sqrt(1680.0) * x * x * x * x * z * z * z * z * z * z * z - std::sqrt(5.767822265625) * x * x * y * y * y * y * y * y * y * y * z + std::sqrt(1025.390625) * x * x * y * y * y * y * y * y * z * z * z - std::sqrt(2307.12890625) * x * x * y * y * y * y * z * z * z * z * z + std::sqrt(26.25) * x * x * z * z * z * z * z * z * z * z * z - std::sqrt(0.640869140625) * y * y * y * y * y * y * y * y * y * y * z + std::sqrt(256.34765625) * y * y * y * y * y * y * y * y * z * z * z - std::sqrt(2307.12890625) * y * y * y * y * y * y * z * z * z * z * z + std::sqrt(1680.0) * y * y * y * y * z * z * z * z * z * z * z - std::sqrt(26.25) * y * y * z * z * z * z * z * z * z * z * z) + e_1 * (-std::sqrt(92.28515625) * x * x * x * x * x * x * x * x * z - std::sqrt(369.140625) * x * x * x * x * x * x * y * y * z - std::sqrt(8505.0) * x * x * x * x * z * z * z * z * z + std::sqrt(369.140625) * x * x * y * y * y * y * y * y * z - std::sqrt(3780.0) * x * x * z * z * z * z * z * z * z + std::sqrt(92.28515625) * y * y * y * y * y * y * y * y * z + std::sqrt(8505.0) * y * y * y * y * z * z * z * z * z + std::sqrt(3780.0) * y * y * z * z * z * z * z * z * z) + e_2 * (-std::sqrt(20764.16015625) * x * x * x * x * x * x * z - std::sqrt(20764.16015625) * x * x * x * x * y * y * z - std::sqrt(212625.0) * x * x * x * x * z * z * z + std::sqrt(20764.16015625) * x * x * y * y * y * y * z - std::sqrt(850500.0) * x * x * z * z * z * z * z + std::sqrt(20764.16015625) * y * y * y * y * y * y * z + std::sqrt(212625.0) * y * y * y * y * z * z * z + std::sqrt(850500.0) * y * y * z * z * z * z * z) + e_3 * (-std::sqrt(1512000.0) * x * x * x * x * z - std::sqrt(15970500.0) * x * x * z * z * z + std::sqrt(1512000.0) * y * y * y * y * z + std::sqrt(15970500.0) * y * y * z * z * z) + e_4 * (-std::sqrt(23441906.25) * x * x * z + std::sqrt(23441906.25) * y * y * z);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        pc_74[k] = e_0 * (std::sqrt(0.026702880859375) * x * x * x * x * x * x * x * x * x * x * x + std::sqrt(0.026702880859375) * x * x * x * x * x * x * x * x * x * y * y - std::sqrt(18.0511474609375) * x * x * x * x * x * x * x * x * x * z * z - std::sqrt(0.9613037109375) * x * x * x * x * x * x * x * y * y * y * y + std::sqrt(753.662109375) * x * x * x * x * x * x * x * z * z * z * z - std::sqrt(5.2337646484375) * x * x * x * x * x * y * y * y * y * y * y + std::sqrt(649.84130859375) * x * x * x * x * x * y * y * y * y * z * z - std::sqrt(753.662109375) * x * x * x * x * x * y * y * z * z * z * z - std::sqrt(1017.4609375) * x * x * x * x * x * z * z * z * z * z * z - std::sqrt(3.231048583984375) * x * x * x * y * y * y * y * y * y * y * y + std::sqrt(1155.2734375) * x * x * x * y * y * y * y * y * y * z * z - std::sqrt(18841.552734375) * x * x * x * y * y * y * y * z * z * z * z + std::sqrt(4069.84375) * x * x * x * y * y * z * z * z * z * z * z + std::sqrt(17.5) * x * x * x * z * z * z * z * z * z * z * z - std::sqrt(0.240325927734375) * x * y * y * y * y * y * y * y * y * y * y + std::sqrt(162.4603271484375) * x * y * y * y * y * y * y * y * y * z * z - std::sqrt(6782.958984375) * x * y * y * y * y * y * y * z * z * z * z + std::sqrt(9157.1484375) * x * y * y * y * y * z * z * z * z * z * z - std::sqrt(157.5) * x * y * y * z * z * z * z * z * z * z * z) + e_1 * (std::sqrt(24.0325927734375) * x * x * x * x * x * x * x * x * x + std::sqrt(138.427734375) * x * x * x * x * x * x * x * z * z - std::sqrt(865.17333984375) * x * x * x * x * x * y * y * y * y - std::sqrt(138.427734375) * x * x * x * x * x * y * y * z * z + std::sqrt(4983.3984375) * x * x * x * x * x * z * z * z * z - std::sqrt(1538.0859375) * x * x * x * y * y * y * y * y * y - std::sqrt(3460.693359375) * x * x * x * y * y * y * y * z * z - std::sqrt(19933.59375) * x * x * x * y * y * z * z * z * z - std::sqrt(19057.5) * x * x * x * z * z * z * z * z * z - std::sqrt(216.2933349609375) * x * y * y * y * y * y * y * y * y - std::sqrt(1245.849609375) * x * y * y * y * y * y * y * z * z - std::sqrt(44850.5859375) * x * y * y * y * y * z * z * z * z + std::sqrt(171517.5) * x * y * y * z * z * z * z * z * z) + e_2 * (std::sqrt(8859.375) * x * x * x * x * x * x * x - std::sqrt(8859.375) * x * x * x * x * x * y * y + std::sqrt(108527.34375) * x * x * x * x * x * z * z - std::sqrt(221484.375) * x * x * x * y * y * y * y - std::sqrt(434109.375) * x * x * x * y * y * z * z - std::sqrt(567000.0) * x * x * x * z * z * z * z - std::sqrt(79734.375) * x * y * y * y * y * y * y - std::sqrt(976746.09375) * x * y * y * y * y * z * z + std::sqrt(5103000.0) * x * y * y * z * z * z * z) + e_3 * (std::sqrt(543621.09375) * x * x * x * x * x - std::sqrt(2174484.375) * x * x * x * y * y - std::sqrt(393750.0) * x * x * x * z * z - std::sqrt(4892589.84375) * x * y * y * y * y + std::sqrt(3543750.0) * x * y * y * z * z) + e_4 * (std::sqrt(1736437.5) * x * x * x - std::sqrt(15627937.5) * x * y * y);
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

        pc_75[k] = e_0 * (-std::sqrt(0.48065185546875) * x * x * x * x * x * x * x * x * x * x * z + std::sqrt(4.32586669921875) * x * x * x * x * x * x * x * x * y * y * z + std::sqrt(155.731201171875) * x * x * x * x * x * x * x * x * z * z * z + std::sqrt(94.207763671875) * x * x * x * x * x * x * y * y * y * y * z - std::sqrt(2491.69921875) * x * x * x * x * x * x * y * y * z * z * z - std::sqrt(276.85546875) * x * x * x * x * x * x * z * z * z * z * z + std::sqrt(94.207763671875) * x * x * x * x * y * y * y * y * y * y * z - std::sqrt(15573.1201171875) * x * x * x * x * y * y * y * y * z * z * z + std::sqrt(6921.38671875) * x * x * x * x * y * y * z * z * z * z * z + std::sqrt(4.921875) * x * x * x * x * z * z * z * z * z * z * z + std::sqrt(4.32586669921875) * x * x * y * y * y * y * y * y * y * y * z - std::sqrt(2491.69921875) * x * x * y * y * y * y * y * y * z * z * z + std::sqrt(6921.38671875) * x * x * y * y * y * y * z * z * z * z * z - std::sqrt(177.1875) * x * x * y * y * z * z * z * z * z * z * z - std::sqrt(0.48065185546875) * y * y * y * y * y * y * y * y * y * y * z + std::sqrt(155.731201171875) * y * y * y * y * y * y * y * y * z * z * z - std::sqrt(276.85546875) * y * y * y * y * y * y * z * z * z * z * z + std::sqrt(4.921875) * y * y * y * y * z * z * z * z * z * z * z) + e_1 * (std::sqrt(69.2138671875) * x * x * x * x * x * x * x * x * z - std::sqrt(1107.421875) * x * x * x * x * x * x * y * y * z + std::sqrt(17718.75) * x * x * x * x * x * x * z * z * z - std::sqrt(6921.38671875) * x * x * x * x * y * y * y * y * z - std::sqrt(442968.75) * x * x * x * x * y * y * z * z * z - std::sqrt(14352.1875) * x * x * x * x * z * z * z * z * z - std::sqrt(1107.421875) * x * x * y * y * y * y * y * y * z - std::sqrt(442968.75) * x * x * y * y * y * y * z * z * z + std::sqrt(516678.75) * x * x * y * y * z * z * z * z * z + std::sqrt(69.2138671875) * y * y * y * y * y * y * y * y * z + std::sqrt(17718.75) * y * y * y * y * y * y * z * z * z - std::sqrt(14352.1875) * y * y * y * y * z * z * z * z * z) + e_2 * (std::sqrt(89701.171875) * x * x * x * x * x * x * z - std::sqrt(2242529.296875) * x * x * x * x * y * y * z + std::sqrt(4429.6875) * x * x * x * x * z * z * z - std::sqrt(2242529.296875) * x * x * y * y * y * y * z - std::sqrt(159468.75) * x * x * y * y * z * z * z + std::sqrt(89701.171875) * y * y * y * y * y * y * z + std::sqrt(4429.6875) * y * y * y * y * z * z * z) + e_3 * (std::sqrt(1134000.0) * x * x * x * x * z - std::sqrt(40824000.0) * x * x * y * y * z + std::sqrt(1134000.0) * y * y * y * y * z);
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

        pc_76[k] = e_0 * (-std::sqrt(0.048065185546875) * x * x * x * x * x * x * x * x * x * x * x + std::sqrt(2.355194091796875) * x * x * x * x * x * x * x * x * x * y * y + std::sqrt(15.5731201171875) * x * x * x * x * x * x * x * x * x * z * z + std::sqrt(23.2635498046875) * x * x * x * x * x * x * x * y * y * y * y - std::sqrt(996.6796875) * x * x * x * x * x * x * x * y * y * z * z - std::sqrt(27.685546875) * x * x * x * x * x * x * x * z * z * z * z + std::sqrt(9.4207763671875) * x * x * x * x * x * y * y * y * y * y * y - std::sqrt(3052.33154296875) * x * x * x * x * x * y * y * y * y * z * z + std::sqrt(2242.529296875) * x * x * x * x * x * y * y * z * z * z * z + std::sqrt(0.4921875) * x * x * x * x * x * z * z * z * z * z * z - std::sqrt(1.201629638671875) * x * x * x * y * y * y * y * y * y * y * y + std::sqrt(692.138671875) * x * x * x * y * y * y * y * z * z * z * z - std::sqrt(49.21875) * x * x * x * y * y * z * z * z * z * z * z - std::sqrt(1.201629638671875) * x * y * y * y * y * y * y * y * y * y * y + std::sqrt(389.3280029296875) * x * y * y * y * y * y * y * y * y * z * z - std::sqrt(692.138671875) * x * y * y * y * y * y * y * z * z * z * z + std::sqrt(12.3046875) * x * y * y * y * y * z * z * z * z * z * z) + e_1 * (-std::sqrt(43.2586669921875) * x * x * x * x * x * x * x * x * x + std::sqrt(2768.5546875) * x * x * x * x * x * x * x * y * y + std::sqrt(6229.248046875) * x * x * x * x * x * x * x * z * z + std::sqrt(8478.69873046875) * x * x * x * x * x * y * y * y * y - std::sqrt(504569.091796875) * x * x * x * x * x * y * y * z * z - std::sqrt(2768.5546875) * x * x * x * x * x * z * z * z * z - std::sqrt(155731.201171875) * x * x * x * y * y * y * y * z * z + std::sqrt(276855.46875) * x * x * x * y * y * z * z * z * z - std::sqrt(1081.4666748046875) * x * y * y * y * y * y * y * y * y + std::sqrt(155731.201171875) * x * y * y * y * y * y * y * z * z - std::sqrt(69213.8671875) * x * y * y * y * y * z * z * z * z) + e_2 * (-std::sqrt(2768.5546875) * x * x * x * x * x * x * x + std::sqrt(224252.9296875) * x * x * x * x * x * y * y + std::sqrt(99667.96875) * x * x * x * x * x * z * z + std::sqrt(69213.8671875) * x * x * x * y * y * y * y - std::sqrt(9966796.875) * x * x * x * y * y * z * z - std::sqrt(69213.8671875) * x * y * y * y * y * y * y + std::sqrt(2491699.21875) * x * y * y * y * y * z * z) + e_3 * (-std::sqrt(11074.21875) * x * x * x * x * x + std::sqrt(1107421.875) * x * x * x * y * y - std::sqrt(276855.46875) * x * y * y * y * y);
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

        pc_77[k] = e_0 * (std::sqrt(100.9368896484375) * x * x * x * x * x * x * x * x * x * y * z - std::sqrt(1614.990234375) * x * x * x * x * x * x * x * y * z * z * z - std::sqrt(791.34521484375) * x * x * x * x * x * y * y * y * y * y * z + std::sqrt(1614.990234375) * x * x * x * x * x * y * y * y * z * z * z + std::sqrt(258.3984375) * x * x * x * x * x * y * z * z * z * z * z - std::sqrt(258.3984375) * x * x * x * y * y * y * y * y * y * y * z + std::sqrt(5232.568359375) * x * x * x * y * y * y * y * y * z * z * z - std::sqrt(1033.59375) * x * x * x * y * y * y * z * z * z * z * z + std::sqrt(4.0374755859375) * x * y * y * y * y * y * y * y * y * y * z - std::sqrt(64.599609375) * x * y * y * y * y * y * y * y * z * z * z + std::sqrt(10.3359375) * x * y * y * y * y * y * z * z * z * z * z) + e_1 * (std::sqrt(58139.6484375) * x * x * x * x * x * x * x * y * z - std::sqrt(25839.84375) * x * x * x * x * x * y * y * y * z - std::sqrt(316538.0859375) * x * x * x * x * x * y * z * z * z - std::sqrt(161499.0234375) * x * x * x * y * y * y * y * y * z + std::sqrt(645996.09375) * x * x * x * y * y * y * z * z * z + std::sqrt(4134.375) * x * x * x * y * z * z * z * z * z + std::sqrt(6459.9609375) * x * y * y * y * y * y * z * z * z - std::sqrt(4134.375) * x * y * y * y * z * z * z * z * z) + e_2 * (std::sqrt(2093027.34375) * x * x * x * x * x * y * z - std::sqrt(2583984.375) * x * x * x * y * y * y * z - std::sqrt(1653750.0) * x * x * x * y * z * z * z - std::sqrt(232558.59375) * x * y * y * y * y * y * z + std::sqrt(1653750.0) * x * y * y * y * z * z * z) + e_3 * (std::sqrt(3720937.5) * x * x * x * y * z - std::sqrt(3720937.5) * x * y * y * y * z);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        pc_78[k] = e_0 * (std::sqrt(645.99609375) * x * x * x * x * x * x * x * x * y * z * z + std::sqrt(645.99609375) * x * x * x * x * x * x * y * y * y * z * z - std::sqrt(10335.9375) * x * x * x * x * x * x * y * z * z * z * z - std::sqrt(645.99609375) * x * x * x * x * y * y * y * y * y * z * z + std::sqrt(1653.75) * x * x * x * x * y * z * z * z * z * z * z - std::sqrt(645.99609375) * x * x * y * y * y * y * y * y * y * z * z + std::sqrt(10335.9375) * x * x * y * y * y * y * y * z * z * z * z - std::sqrt(1653.75) * x * x * y * y * y * z * z * z * z * z * z) + e_1 * (std::sqrt(645.99609375) * x * x * x * x * x * x * x * x * y + std::sqrt(645.99609375) * x * x * x * x * x * x * y * y * y + std::sqrt(31653.80859375) * x * x * x * x * x * x * y * z * z - std::sqrt(645.99609375) * x * x * x * x * y * y * y * y * y + std::sqrt(16149.90234375) * x * x * x * x * y * y * y * z * z - std::sqrt(837210.9375) * x * x * x * x * y * z * z * z * z - std::sqrt(645.99609375) * x * x * y * y * y * y * y * y * y - std::sqrt(5813.96484375) * x * x * y * y * y * y * y * z * z + std::sqrt(165375.0) * x * x * y * y * y * z * z * z * z + std::sqrt(14883.75) * x * x * y * z * z * z * z * z * z - std::sqrt(645.99609375) * y * y * y * y * y * y * y * z * z + std::sqrt(10335.9375) * y * y * y * y * y * z * z * z * z - std::sqrt(1653.75) * y * y * y * z * z * z * z * z * z) + e_2 * (std::sqrt(233204.58984375) * x * x * x * x * x * x * y + std::sqrt(16149.90234375) * x * x * x * x * y * y * y - std::sqrt(1488375.0) * x * x * x * x * y * z * z - std::sqrt(145349.12109375) * x * x * y * y * y * y * y + std::sqrt(1488375.0) * x * x * y * y * y * z * z - std::sqrt(1488375.0) * x * x * y * z * z * z * z - std::sqrt(645.99609375) * y * y * y * y * y * y * y + std::sqrt(165375.0) * y * y * y * z * z * z * z) + e_3 * (std::sqrt(4558148.4375) * x * x * x * x * y - std::sqrt(372093.75) * x * x * y * y * y - std::sqrt(13395375.0) * x * x * y * z * z - std::sqrt(93023.4375) * y * y * y * y * y + std::sqrt(1488375.0) * y * y * y * z * z) + e_4 * (std::sqrt(3348843.75) * x * x * y - std::sqrt(372093.75) * y * y * y);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        pc_79[k] = e_0 * (-std::sqrt(20.1873779296875) * x * x * x * x * x * x * x * x * x * y * z - std::sqrt(143.5546875) * x * x * x * x * x * x * x * y * y * y * z + std::sqrt(2906.982421875) * x * x * x * x * x * x * x * y * z * z * z - std::sqrt(80.74951171875) * x * x * x * x * x * y * y * y * y * y * z + std::sqrt(8074.951171875) * x * x * x * x * x * y * y * y * z * z * z - std::sqrt(22790.7421875) * x * x * x * x * x * y * z * z * z * z * z + std::sqrt(322.998046875) * x * x * x * y * y * y * y * y * z * z * z - std::sqrt(10129.21875) * x * x * x * y * y * y * z * z * z * z * z + std::sqrt(3307.5) * x * x * x * y * z * z * z * z * z * z * z + std::sqrt(2.2430419921875) * x * y * y * y * y * y * y * y * y * y * z - std::sqrt(322.998046875) * x * y * y * y * y * y * y * y * z * z * z + std::sqrt(2532.3046875) * x * y * y * y * y * y * z * z * z * z * z - std::sqrt(367.5) * x * y * y * y * z * z * z * z * z * z * z) + e_1 * (-std::sqrt(1291.9921875) * x * x * x * x * x * x * x * y * z - std::sqrt(5167.96875) * x * x * x * x * x * y * y * y * z - std::sqrt(11627.9296875) * x * x * x * x * x * y * z * z * z - std::sqrt(1291.9921875) * x * x * x * y * y * y * y * y * z + std::sqrt(5167.96875) * x * x * x * y * y * y * z * z * z - std::sqrt(364651.875) * x * x * x * y * z * z * z * z * z + std::sqrt(32299.8046875) * x * y * y * y * y * y * z * z * z - std::sqrt(40516.875) * x * y * y * y * z * z * z * z * z + std::sqrt(13230.0) * x * y * z * z * z * z * z * z * z) + e_2 * (-std::sqrt(418605.46875) * x * x * x * x * x * y * z - std::sqrt(186046.875) * x * x * x * y * y * y * z - std::sqrt(11907000.0) * x * x * x * y * z * z * z + std::sqrt(46511.71875) * x * y * y * y * y * y * z) + e_3 * (-std::sqrt(36465187.5) * x * x * x * y * z + std::sqrt(82687.5) * x * y * y * y * z - std::sqrt(11907000.0) * x * y * z * z * z) + e_4 * (-std::sqrt(47628000.0) * x * y * z);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        pc_80[k] = e_0 * (-std::sqrt(215.33203125) * x * x * x * x * x * x * x * x * y * z * z - std::sqrt(1937.98828125) * x * x * x * x * x * x * y * y * y * z * z + std::sqrt(7751.953125) * x * x * x * x * x * x * y * z * z * z * z - std::sqrt(1937.98828125) * x * x * x * x * y * y * y * y * y * z * z + std::sqrt(31007.8125) * x * x * x * x * y * y * y * z * z * z * z - std::sqrt(19845.0) * x * x * x * x * y * z * z * z * z * z * z - std::sqrt(215.33203125) * x * x * y * y * y * y * y * y * y * z * z + std::sqrt(7751.953125) * x * x * y * y * y * y * y * z * z * z * z - std::sqrt(19845.0) * x * x * y * y * y * z * z * z * z * z * z + std::sqrt(2205.0) * x * x * y * z * z * z * z * z * z * z * z) + e_1 * (-std::sqrt(215.33203125) * x * x * x * x * x * x * x * x * y - std::sqrt(1937.98828125) * x * x * x * x * x * x * y * y * y - std::sqrt(215.33203125) * x * x * x * x * x * x * y * z * z - std::sqrt(1937.98828125) * x * x * x * x * y * y * y * y * y - std::sqrt(1937.98828125) * x * x * x * x * y * y * y * z * z - std::sqrt(69767.578125) * x * x * x * x * y * z * z * z * z - std::sqrt(215.33203125) * x * x * y * y * y * y * y * y * y - std::sqrt(1937.98828125) * x * x * y * y * y * y * y * z * z - std::sqrt(31007.8125) * x * x * y * y * y * z * z * z * z + std::sqrt(2205.0) * x * x * y * z * z * z * z * z * z - std::sqrt(215.33203125) * y * y * y * y * y * y * y * z * z + std::sqrt(7751.953125) * y * y * y * y * y * z * z * z * z - std::sqrt(19845.0) * y * y * y * z * z * z * z * z * z + std::sqrt(2205.0) * y * z * z * z * z * z * z * z * z) + e_2 * (-std::sqrt(77734.86328125) * x * x * x * x * x * x * y - std::sqrt(327520.01953125) * x * x * x * x * y * y * y - std::sqrt(937986.328125) * x * x * x * x * y * z * z - std::sqrt(94961.42578125) * x * x * y * y * y * y * y - std::sqrt(775195.3125) * x * x * y * y * y * z * z - std::sqrt(496125.0) * x * x * y * z * z * z * z - std::sqrt(215.33203125) * y * y * y * y * y * y * y + std::sqrt(7751.953125) * y * y * y * y * y * z * z - std::sqrt(496125.0) * y * y * y * z * z * z * z + std::sqrt(220500.0) * y * z * z * z * z * z * z) + e_3 * (-std::sqrt(5240320.3125) * x * x * x * x * y - std::sqrt(6077531.25) * x * x * y * y * y - std::sqrt(17860500.0) * x * x * y * z * z - std::sqrt(31007.8125) * y * y * y * y * y - std::sqrt(1984500.0) * y * y * y * z * z + std::sqrt(1984500.0) * y * z * z * z * z) + e_4 * (-std::sqrt(40186125.0) * x * x * y - std::sqrt(1984500.0) * y * y * y) + e_5 * (-std::sqrt(4465125.0) * y);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        pc_81[k] = e_0 * (std::sqrt(1.922607421875) * x * x * x * x * x * x * x * x * x * y * z + std::sqrt(30.76171875) * x * x * x * x * x * x * x * y * y * y * z - std::sqrt(492.1875) * x * x * x * x * x * x * x * y * z * z * z + std::sqrt(69.2138671875) * x * x * x * x * x * y * y * y * y * y * z - std::sqrt(4429.6875) * x * x * x * x * x * y * y * y * z * z * z + std::sqrt(6378.75) * x * x * x * x * x * y * z * z * z * z * z + std::sqrt(30.76171875) * x * x * x * y * y * y * y * y * y * y * z - std::sqrt(4429.6875) * x * x * x * y * y * y * y * y * z * z * z + std::sqrt(25515.0) * x * x * x * y * y * y * z * z * z * z * z - std::sqrt(5040.0) * x * x * x * y * z * z * z * z * z * z * z + std::sqrt(1.922607421875) * x * y * y * y * y * y * y * y * y * y * z - std::sqrt(492.1875) * x * y * y * y * y * y * y * y * z * z * z + std::sqrt(6378.75) * x * y * y * y * y * y * z * z * z * z * z - std::sqrt(5040.0) * x * y * y * y * z * z * z * z * z * z * z + std::sqrt(315.0) * x * y * z * z * z * z * z * z * z * z * z) + e_1 * (std::sqrt(17718.75) * x * x * x * x * x * y * z * z * z + std::sqrt(70875.0) * x * x * x * y * y * y * z * z * z - std::sqrt(45360.0) * x * x * x * y * z * z * z * z * z + std::sqrt(17718.75) * x * y * y * y * y * y * z * z * z - std::sqrt(45360.0) * x * y * y * y * z * z * z * z * z + std::sqrt(45360.0) * x * y * z * z * z * z * z * z * z) + e_2 * (std::sqrt(39867.1875) * x * x * x * x * x * y * z + std::sqrt(159468.75) * x * x * x * y * y * y * z + std::sqrt(39867.1875) * x * y * y * y * y * y * z + std::sqrt(2551500.0) * x * y * z * z * z * z * z) + e_3 * (std::sqrt(1134000.0) * x * x * x * y * z + std::sqrt(1134000.0) * x * y * y * y * z + std::sqrt(28350000.0) * x * y * z * z * z) + e_4 * (std::sqrt(31255875.0) * x * y * z);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        pc_82[k] = e_0 * (std::sqrt(28.839111328125) * x * x * x * x * x * x * x * x * x * z * z + std::sqrt(461.42578125) * x * x * x * x * x * x * x * y * y * z * z - std::sqrt(1281.73828125) * x * x * x * x * x * x * x * z * z * z * z + std::sqrt(1038.2080078125) * x * x * x * x * x * y * y * y * y * z * z - std::sqrt(11535.64453125) * x * x * x * x * x * y * y * z * z * z * z + std::sqrt(4725.0) * x * x * x * x * x * z * z * z * z * z * z + std::sqrt(461.42578125) * x * x * x * y * y * y * y * y * y * z * z - std::sqrt(11535.64453125) * x * x * x * y * y * y * y * z * z * z * z + std::sqrt(18900.0) * x * x * x * y * y * z * z * z * z * z * z - std::sqrt(1181.25) * x * x * x * z * z * z * z * z * z * z * z + std::sqrt(28.839111328125) * x * y * y * y * y * y * y * y * y * z * z - std::sqrt(1281.73828125) * x * y * y * y * y * y * y * z * z * z * z + std::sqrt(4725.0) * x * y * y * y * y * z * z * z * z * z * z - std::sqrt(1181.25) * x * y * y * z * z * z * z * z * z * z * z + std::sqrt(21.0) * x * z * z * z * z * z * z * z * z * z * z) + e_1 * (std::sqrt(28.839111328125) * x * x * x * x * x * x * x * x * x + std::sqrt(461.42578125) * x * x * x * x * x * x * x * y * y + std::sqrt(1038.2080078125) * x * x * x * x * x * y * y * y * y + std::sqrt(29531.25) * x * x * x * x * x * z * z * z * z + std::sqrt(461.42578125) * x * x * x * y * y * y * y * y * y + std::sqrt(118125.0) * x * x * x * y * y * z * z * z * z - std::sqrt(18900.0) * x * x * x * z * z * z * z * z * z + std::sqrt(28.839111328125) * x * y * y * y * y * y * y * y * y + std::sqrt(29531.25) * x * y * y * y * y * z * z * z * z - std::sqrt(18900.0) * x * y * y * z * z * z * z * z * z + std::sqrt(4725.0) * x * z * z * z * z * z * z * z * z) + e_2 * (std::sqrt(11535.64453125) * x * x * x * x * x * x * x + std::sqrt(103820.80078125) * x * x * x * x * x * y * y + std::sqrt(265781.25) * x * x * x * x * x * z * z + std::sqrt(103820.80078125) * x * x * x * y * y * y * y + std::sqrt(1063125.0) * x * x * x * y * y * z * z + std::sqrt(11535.64453125) * x * y * y * y * y * y * y + std::sqrt(265781.25) * x * y * y * y * y * z * z + std::sqrt(472500.0) * x * z * z * z * z * z * z) + e_3 * (std::sqrt(1063125.0) * x * x * x * x * x + std::sqrt(4252500.0) * x * x * x * y * y + std::sqrt(4252500.0) * x * x * x * z * z + std::sqrt(1063125.0) * x * y * y * y * y + std::sqrt(4252500.0) * x * y * y * z * z + std::sqrt(11812500.0) * x * z * z * z * z) + e_4 * (std::sqrt(13023281.25) * x * x * x + std::sqrt(13023281.25) * x * y * y + std::sqrt(52093125.0) * x * z * z) + e_5 * (std::sqrt(18753525.0) * x);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        pc_83[k] = e_0 * (std::sqrt(1.922607421875) * x * x * x * x * x * x * x * x * x * x * z + std::sqrt(30.76171875) * x * x * x * x * x * x * x * x * y * y * z - std::sqrt(492.1875) * x * x * x * x * x * x * x * x * z * z * z + std::sqrt(69.2138671875) * x * x * x * x * x * x * y * y * y * y * z - std::sqrt(4429.6875) * x * x * x * x * x * x * y * y * z * z * z + std::sqrt(6378.75) * x * x * x * x * x * x * z * z * z * z * z + std::sqrt(30.76171875) * x * x * x * x * y * y * y * y * y * y * z - std::sqrt(4429.6875) * x * x * x * x * y * y * y * y * z * z * z + std::sqrt(25515.0) * x * x * x * x * y * y * z * z * z * z * z - std::sqrt(5040.0) * x * x * x * x * z * z * z * z * z * z * z + std::sqrt(1.922607421875) * x * x * y * y * y * y * y * y * y * y * z - std::sqrt(492.1875) * x * x * y * y * y * y * y * y * z * z * z + std::sqrt(6378.75) * x * x * y * y * y * y * z * z * z * z * z - std::sqrt(5040.0) * x * x * y * y * z * z * z * z * z * z * z + std::sqrt(315.0) * x * x * z * z * z * z * z * z * z * z * z) + e_1 * (std::sqrt(1.922607421875) * x * x * x * x * x * x * x * x * z + std::sqrt(30.76171875) * x * x * x * x * x * x * y * y * z + std::sqrt(12304.6875) * x * x * x * x * x * x * z * z * z + std::sqrt(69.2138671875) * x * x * x * x * y * y * y * y * z + std::sqrt(39867.1875) * x * x * x * x * y * y * z * z * z - std::sqrt(17718.75) * x * x * x * x * z * z * z * z * z + std::sqrt(30.76171875) * x * x * y * y * y * y * y * y * z + std::sqrt(4429.6875) * x * x * y * y * y * y * z * z * z - std::sqrt(2835.0) * x * x * y * y * z * z * z * z * z + std::sqrt(20160.0) * x * x * z * z * z * z * z * z * z + std::sqrt(1.922607421875) * y * y * y * y * y * y * y * y * z - std::sqrt(492.1875) * y * y * y * y * y * y * z * z * z + std::sqrt(6378.75) * y * y * y * y * z * z * z * z * z - std::sqrt(5040.0) * y * y * z * z * z * z * z * z * z + std::sqrt(315.0) * z * z * z * z * z * z * z * z * z) + e_2 * (std::sqrt(35560.546875) * x * x * x * x * x * x * z + std::sqrt(133998.046875) * x * x * x * x * y * y * z + std::sqrt(70875.0) * x * x * x * x * z * z * z + std::sqrt(27685.546875) * x * x * y * y * y * y * z + std::sqrt(283500.0) * x * x * y * y * z * z * z + std::sqrt(1134000.0) * x * x * z * z * z * z * z - std::sqrt(123.046875) * y * y * y * y * y * y * z + std::sqrt(70875.0) * y * y * y * y * z * z * z - std::sqrt(283500.0) * y * y * z * z * z * z * z + std::sqrt(126000.0) * z * z * z * z * z * z * z) + e_3 * (std::sqrt(1771875.0) * x * x * x * x * z + std::sqrt(2551500.0) * x * x * y * y * z + std::sqrt(18144000.0) * x * x * z * z * z + std::sqrt(70875.0) * y * y * y * y * z - std::sqrt(1134000.0) * y * y * z * z * z + std::sqrt(7087500.0) * z * z * z * z * z) + e_4 * (std::sqrt(31255875.0) * x * x * z + std::sqrt(55566000.0) * z * z * z) + e_5 * (std::sqrt(31255875.0) * z);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        pc_84[k] = e_0 * (-std::sqrt(53.8330078125) * x * x * x * x * x * x * x * x * x * z * z - std::sqrt(215.33203125) * x * x * x * x * x * x * x * y * y * z * z + std::sqrt(1937.98828125) * x * x * x * x * x * x * x * z * z * z * z + std::sqrt(1937.98828125) * x * x * x * x * x * y * y * z * z * z * z - std::sqrt(4961.25) * x * x * x * x * x * z * z * z * z * z * z + std::sqrt(215.33203125) * x * x * x * y * y * y * y * y * y * z * z - std::sqrt(1937.98828125) * x * x * x * y * y * y * y * z * z * z * z + std::sqrt(551.25) * x * x * x * z * z * z * z * z * z * z * z + std::sqrt(53.8330078125) * x * y * y * y * y * y * y * y * y * z * z - std::sqrt(1937.98828125) * x * y * y * y * y * y * y * z * z * z * z + std::sqrt(4961.25) * x * y * y * y * y * z * z * z * z * z * z - std::sqrt(551.25) * x * y * y * z * z * z * z * z * z * z * z) + e_1 * (-std::sqrt(53.8330078125) * x * x * x * x * x * x * x * x * x - std::sqrt(215.33203125) * x * x * x * x * x * x * x * y * y - std::sqrt(215.33203125) * x * x * x * x * x * x * x * z * z - std::sqrt(1937.98828125) * x * x * x * x * x * y * y * z * z - std::sqrt(7751.953125) * x * x * x * x * x * z * z * z * z + std::sqrt(215.33203125) * x * x * x * y * y * y * y * y * y - std::sqrt(1937.98828125) * x * x * x * y * y * y * y * z * z + std::sqrt(31007.8125) * x * x * x * y * y * z * z * z * z - std::sqrt(2205.0) * x * x * x * z * z * z * z * z * z + std::sqrt(53.8330078125) * x * y * y * y * y * y * y * y * y - std::sqrt(215.33203125) * x * y * y * y * y * y * y * z * z + std::sqrt(69767.578125) * x * y * y * y * y * z * z * z * z - std::sqrt(55125.0) * x * y * y * z * z * z * z * z * z + std::sqrt(2205.0) * x * z * z * z * z * z * z * z * z) + e_2 * (-std::sqrt(21533.203125) * x * x * x * x * x * x * x - std::sqrt(31007.8125) * x * x * x * x * x * y * y - std::sqrt(193798.828125) * x * x * x * x * x * z * z + std::sqrt(7751.953125) * x * x * x * y * y * y * y + std::sqrt(31007.8125) * x * x * x * y * y * z * z - std::sqrt(496125.0) * x * x * x * z * z * z * z + std::sqrt(13781.25) * x * y * y * y * y * y * y + std::sqrt(379845.703125) * x * y * y * y * y * z * z - std::sqrt(496125.0) * x * y * y * z * z * z * z + std::sqrt(220500.0) * x * z * z * z * z * z * z) + e_3 * (-std::sqrt(1519382.8125) * x * x * x * x * x - std::sqrt(124031.25) * x * x * x * y * y - std::sqrt(7938000.0) * x * x * x * z * z + std::sqrt(775195.3125) * x * y * y * y * y + std::sqrt(1984500.0) * x * z * z * z * z) + e_4 * (-std::sqrt(15007781.25) * x * x * x + std::sqrt(1116281.25) * x * y * y) + e_5 * (-std::sqrt(4465125.0) * x);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        pc_85[k] = e_0 * (-std::sqrt(2.2430419921875) * x * x * x * x * x * x * x * x * x * x * z + std::sqrt(322.998046875) * x * x * x * x * x * x * x * x * z * z * z + std::sqrt(80.74951171875) * x * x * x * x * x * x * y * y * y * y * z - std::sqrt(322.998046875) * x * x * x * x * x * x * y * y * z * z * z - std::sqrt(2532.3046875) * x * x * x * x * x * x * z * z * z * z * z + std::sqrt(143.5546875) * x * x * x * x * y * y * y * y * y * y * z - std::sqrt(8074.951171875) * x * x * x * x * y * y * y * y * z * z * z + std::sqrt(10129.21875) * x * x * x * x * y * y * z * z * z * z * z + std::sqrt(367.5) * x * x * x * x * z * z * z * z * z * z * z + std::sqrt(20.1873779296875) * x * x * y * y * y * y * y * y * y * y * z - std::sqrt(2906.982421875) * x * x * y * y * y * y * y * y * z * z * z + std::sqrt(22790.7421875) * x * x * y * y * y * y * z * z * z * z * z - std::sqrt(3307.5) * x * x * y * y * z * z * z * z * z * z * z) + e_1 * (-std::sqrt(181.6864013671875) * x * x * x * x * x * x * x * x * z - std::sqrt(322.998046875) * x * x * x * x * x * x * z * z * z + std::sqrt(2018.73779296875) * x * x * x * x * y * y * y * y * z + std::sqrt(39082.763671875) * x * x * x * x * y * y * z * z * z - std::sqrt(63307.6171875) * x * x * x * x * z * z * z * z * z + std::sqrt(1291.9921875) * x * x * y * y * y * y * y * y * z + std::sqrt(26162.841796875) * x * x * y * y * y * y * z * z * z + std::sqrt(91162.96875) * x * x * y * y * z * z * z * z * z + std::sqrt(3307.5) * x * x * z * z * z * z * z * z * z + std::sqrt(20.1873779296875) * y * y * y * y * y * y * y * y * z - std::sqrt(2906.982421875) * y * y * y * y * y * y * z * z * z + std::sqrt(22790.7421875) * y * y * y * y * z * z * z * z * z - std::sqrt(3307.5) * y * y * z * z * z * z * z * z * z) + e_2 * (-std::sqrt(46511.71875) * x * x * x * x * x * x * z + std::sqrt(186046.875) * x * x * x * x * y * y * z - std::sqrt(1674421.875) * x * x * x * x * z * z * z + std::sqrt(418605.46875) * x * x * y * y * y * y * z + std::sqrt(6697687.5) * x * x * y * y * z * z * z + std::sqrt(186046.875) * y * y * y * y * z * z * z) + e_3 * (-std::sqrt(4966417.96875) * x * x * x * x * z + std::sqrt(22511671.875) * x * x * y * y * z - std::sqrt(2976750.0) * x * x * z * z * z + std::sqrt(418605.46875) * y * y * y * y * z + std::sqrt(2976750.0) * y * y * z * z * z) + e_4 * (-std::sqrt(11907000.0) * x * x * z + std::sqrt(11907000.0) * y * y * z);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        pc_86[k] = e_0 * (std::sqrt(40.374755859375) * x * x * x * x * x * x * x * x * x * z * z - std::sqrt(645.99609375) * x * x * x * x * x * x * x * y * y * z * z - std::sqrt(645.99609375) * x * x * x * x * x * x * x * z * z * z * z - std::sqrt(4037.4755859375) * x * x * x * x * x * y * y * y * y * z * z + std::sqrt(16149.90234375) * x * x * x * x * x * y * y * z * z * z * z + std::sqrt(103.359375) * x * x * x * x * x * z * z * z * z * z * z - std::sqrt(645.99609375) * x * x * x * y * y * y * y * y * y * z * z + std::sqrt(16149.90234375) * x * x * x * y * y * y * y * z * z * z * z - std::sqrt(3720.9375) * x * x * x * y * y * z * z * z * z * z * z + std::sqrt(40.374755859375) * x * y * y * y * y * y * y * y * y * z * z - std::sqrt(645.99609375) * x * y * y * y * y * y * y * z * z * z * z + std::sqrt(103.359375) * x * y * y * y * y * z * z * z * z * z * z) + e_1 * (std::sqrt(40.374755859375) * x * x * x * x * x * x * x * x * x - std::sqrt(645.99609375) * x * x * x * x * x * x * x * y * y + std::sqrt(2583.984375) * x * x * x * x * x * x * x * z * z - std::sqrt(4037.4755859375) * x * x * x * x * x * y * y * y * y - std::sqrt(23255.859375) * x * x * x * x * x * y * y * z * z - std::sqrt(64599.609375) * x * x * x * x * x * z * z * z * z - std::sqrt(645.99609375) * x * x * x * y * y * y * y * y * y - std::sqrt(64599.609375) * x * x * x * y * y * y * y * z * z + std::sqrt(1250648.4375) * x * x * x * y * y * z * z * z * z + std::sqrt(1653.75) * x * x * x * z * z * z * z * z * z + std::sqrt(40.374755859375) * x * y * y * y * y * y * y * y * y - std::sqrt(2583.984375) * x * y * y * y * y * y * y * z * z + std::sqrt(23255.859375) * x * y * y * y * y * z * z * z * z - std::sqrt(14883.75) * x * y * y * z * z * z * z * z * z) + e_2 * (std::sqrt(16149.90234375) * x * x * x * x * x * x * x - std::sqrt(284884.27734375) * x * x * x * x * x * y * y - std::sqrt(93023.4375) * x * x * x * x * x * z * z - std::sqrt(403747.55859375) * x * x * x * y * y * y * y + std::sqrt(3348843.75) * x * x * x * y * y * z * z - std::sqrt(165375.0) * x * x * x * z * z * z * z + std::sqrt(645.99609375) * x * y * y * y * y * y * y - std::sqrt(93023.4375) * x * y * y * y * y * z * z + std::sqrt(1488375.0) * x * y * y * z * z * z * z) + e_3 * (std::sqrt(372093.75) * x * x * x * x * x - std::sqrt(5953500.0) * x * x * x * y * y - std::sqrt(1488375.0) * x * x * x * z * z - std::sqrt(372093.75) * x * y * y * y * y + std::sqrt(13395375.0) * x * y * y * z * z) + e_4 * (std::sqrt(372093.75) * x * x * x - std::sqrt(3348843.75) * x * y * y);
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

        pc_87[k] = e_0 * (std::sqrt(4.0374755859375) * x * x * x * x * x * x * x * x * x * x * z - std::sqrt(258.3984375) * x * x * x * x * x * x * x * x * y * y * z - std::sqrt(64.599609375) * x * x * x * x * x * x * x * x * z * z * z - std::sqrt(791.34521484375) * x * x * x * x * x * x * y * y * y * y * z + std::sqrt(5232.568359375) * x * x * x * x * x * x * y * y * z * z * z + std::sqrt(10.3359375) * x * x * x * x * x * x * z * z * z * z * z + std::sqrt(1614.990234375) * x * x * x * x * y * y * y * y * z * z * z - std::sqrt(1033.59375) * x * x * x * x * y * y * z * z * z * z * z + std::sqrt(100.9368896484375) * x * x * y * y * y * y * y * y * y * y * z - std::sqrt(1614.990234375) * x * x * y * y * y * y * y * y * z * z * z + std::sqrt(258.3984375) * x * x * y * y * y * y * z * z * z * z * z) + e_1 * (std::sqrt(2523.4222412109375) * x * x * x * x * x * x * x * x * z - std::sqrt(161499.0234375) * x * x * x * x * x * x * y * y * z - std::sqrt(14534.912109375) * x * x * x * x * x * x * z * z * z - std::sqrt(90843.20068359375) * x * x * x * x * y * y * y * y * z + std::sqrt(1009368.896484375) * x * x * x * x * y * y * z * z * z + std::sqrt(258.3984375) * x * x * x * x * z * z * z * z * z + std::sqrt(25839.84375) * x * x * y * y * y * y * y * y * z - std::sqrt(40374.755859375) * x * x * y * y * y * y * z * z * z - std::sqrt(9302.34375) * x * x * y * y * z * z * z * z * z + std::sqrt(100.9368896484375) * y * y * y * y * y * y * y * y * z - std::sqrt(1614.990234375) * y * y * y * y * y * y * z * z * z + std::sqrt(258.3984375) * y * y * y * y * z * z * z * z * z) + e_2 * (std::sqrt(103359.375) * x * x * x * x * x * x * z - std::sqrt(5813964.84375) * x * x * x * x * y * y * z - std::sqrt(103359.375) * x * x * x * x * z * z * z + std::sqrt(3720937.5) * x * x * y * y * z * z * z + std::sqrt(25839.84375) * y * y * y * y * y * y * z - std::sqrt(103359.375) * y * y * y * y * z * z * z) + e_3 * (std::sqrt(232558.59375) * x * x * x * x * z - std::sqrt(8372109.375) * x * x * y * y * z + std::sqrt(232558.59375) * y * y * y * y * z);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        pc_88[k] = e_0 * (std::sqrt(2.5234222412109375) * x * x * x * x * x * x * x * x * x * x * y - std::sqrt(2.5234222412109375) * x * x * x * x * x * x * x * x * y * y * y - std::sqrt(645.99609375) * x * x * x * x * x * x * x * x * y * z * z - std::sqrt(19.78363037109375) * x * x * x * x * x * x * y * y * y * y * y + std::sqrt(2583.984375) * x * x * x * x * x * x * y * y * y * z * z + std::sqrt(645.99609375) * x * x * x * x * x * x * y * z * z * z * z + std::sqrt(3.63372802734375) * x * x * x * x * y * y * y * y * y * y * y + std::sqrt(413.4375) * x * x * x * x * y * y * y * y * y * z * z - std::sqrt(5813.96484375) * x * x * x * x * y * y * y * z * z * z * z + std::sqrt(8.175888061523438) * x * x * y * y * y * y * y * y * y * y * y - std::sqrt(2583.984375) * x * x * y * y * y * y * y * y * y * z * z + std::sqrt(3126.62109375) * x * x * y * y * y * y * y * z * z * z * z - std::sqrt(0.1009368896484375) * y * y * y * y * y * y * y * y * y * y * y + std::sqrt(25.83984375) * y * y * y * y * y * y * y * y * y * z * z - std::sqrt(25.83984375) * y * y * y * y * y * y * y * z * z * z * z) + e_1 * (std::sqrt(1705.8334350585938) * x * x * x * x * x * x * x * x * y - std::sqrt(1453.4912109375) * x * x * x * x * x * x * y * y * y - std::sqrt(165375.0) * x * x * x * x * x * x * y * z * z - std::sqrt(40.374755859375) * x * x * x * x * y * y * y * y * y + std::sqrt(258398.4375) * x * x * x * x * y * y * y * z * z + std::sqrt(23255.859375) * x * x * x * x * y * z * z * z * z + std::sqrt(4037.4755859375) * x * x * y * y * y * y * y * y * y - std::sqrt(372093.75) * x * x * y * y * y * y * y * z * z + std::sqrt(10335.9375) * x * x * y * y * y * z * z * z * z - std::sqrt(90.84320068359375) * y * y * y * y * y * y * y * y * y + std::sqrt(10335.9375) * y * y * y * y * y * y * y * z * z - std::sqrt(2583.984375) * y * y * y * y * y * z * z * z * z) + e_2 * (std::sqrt(100936.8896484375) * x * x * x * x * x * x * y - std::sqrt(4037.4755859375) * x * x * x * x * y * y * y - std::sqrt(3348843.75) * x * x * x * x * y * z * z + std::sqrt(117732.7880859375) * x * x * y * y * y * y * y - std::sqrt(1488375.0) * x * x * y * y * y * z * z + std::sqrt(372093.75) * x * x * y * z * z * z * z - std::sqrt(7913.4521484375) * y * y * y * y * y * y * y + std::sqrt(372093.75) * y * y * y * y * y * z * z - std::sqrt(41343.75) * y * y * y * z * z * z * z) + e_3 * (std::sqrt(837210.9375) * x * x * x * x * y + std::sqrt(372093.75) * x * x * y * y * y - std::sqrt(13395375.0) * x * x * y * z * z - std::sqrt(93023.4375) * y * y * y * y * y + std::sqrt(1488375.0) * y * y * y * z * z) + e_4 * (std::sqrt(837210.9375) * x * x * y - std::sqrt(93023.4375) * y * y * y);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        pc_89[k] = e_0 * (std::sqrt(16.14990234375) * x * x * x * x * x * x * x * x * x * y * z - std::sqrt(4134.375) * x * x * x * x * x * x * x * y * z * z * z - std::sqrt(64.599609375) * x * x * x * x * x * y * y * y * y * y * z + std::sqrt(4134.375) * x * x * x * x * x * y * y * y * z * z * z + std::sqrt(4134.375) * x * x * x * x * x * y * z * z * z * z * z + std::sqrt(4134.375) * x * x * x * y * y * y * y * y * z * z * z - std::sqrt(16537.5) * x * x * x * y * y * y * z * z * z * z * z + std::sqrt(16.14990234375) * x * y * y * y * y * y * y * y * y * y * z - std::sqrt(4134.375) * x * y * y * y * y * y * y * y * z * z * z + std::sqrt(4134.375) * x * y * y * y * y * y * z * z * z * z * z) + e_1 * (-std::sqrt(2325.5859375) * x * x * x * x * x * x * x * y * z + std::sqrt(12661.5234375) * x * x * x * x * x * y * y * y * z - std::sqrt(264600.0) * x * x * x * x * x * y * z * z * z + std::sqrt(12661.5234375) * x * x * x * y * y * y * y * y * z + std::sqrt(66150.0) * x * x * x * y * z * z * z * z * z - std::sqrt(2325.5859375) * x * y * y * y * y * y * y * y * z - std::sqrt(264600.0) * x * y * y * y * y * y * z * z * z + std::sqrt(66150.0) * x * y * y * y * z * z * z * z * z) + e_2 * (-std::sqrt(1230234.9609375) * x * x * x * x * x * y * z + std::sqrt(1266152.34375) * x * x * x * y * y * y * z - std::sqrt(1653750.0) * x * x * x * y * z * z * z - std::sqrt(1230234.9609375) * x * y * y * y * y * y * z - std::sqrt(1653750.0) * x * y * y * y * z * z * z + std::sqrt(595350.0) * x * y * z * z * z * z * z) + e_3 * (-std::sqrt(14883750.0) * x * x * x * y * z - std::sqrt(14883750.0) * x * y * y * y * z) + e_4 * (-std::sqrt(33488437.5) * x * y * z);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        pc_90[k] = e_0 * (-std::sqrt(0.5046844482421875) * x * x * x * x * x * x * x * x * x * x * y - std::sqrt(1.4019012451171875) * x * x * x * x * x * x * x * x * y * y * y + std::sqrt(290.6982421875) * x * x * x * x * x * x * x * x * y * z * z + std::sqrt(0.22430419921875) * x * x * x * x * x * x * y * y * y * y * y + std::sqrt(129.19921875) * x * x * x * x * x * x * y * y * y * z * z - std::sqrt(10465.13671875) * x * x * x * x * x * x * y * z * z * z * z + std::sqrt(2.01873779296875) * x * x * x * x * y * y * y * y * y * y * y - std::sqrt(516.796875) * x * x * x * x * y * y * y * y * y * z * z + std::sqrt(1162.79296875) * x * x * x * x * y * y * y * z * z * z * z + std::sqrt(8268.75) * x * x * x * x * y * z * z * z * z * z * z + std::sqrt(0.0560760498046875) * x * x * y * y * y * y * y * y * y * y * y - std::sqrt(129.19921875) * x * x * y * y * y * y * y * y * y * z * z + std::sqrt(10465.13671875) * x * x * y * y * y * y * y * z * z * z * z - std::sqrt(14700.0) * x * x * y * y * y * z * z * z * z * z * z - std::sqrt(0.0560760498046875) * y * y * y * y * y * y * y * y * y * y * y + std::sqrt(32.2998046875) * y * y * y * y * y * y * y * y * y * z * z - std::sqrt(1162.79296875) * y * y * y * y * y * y * y * z * z * z * z + std::sqrt(918.75) * y * y * y * y * y * z * z * z * z * z * z) + e_1 * (-std::sqrt(341.16668701171875) * x * x * x * x * x * x * x * x * y - std::sqrt(290.6982421875) * x * x * x * x * x * x * y * y * y - std::sqrt(10465.13671875) * x * x * x * x * x * x * y * z * z + std::sqrt(201.873779296875) * x * x * x * x * y * y * y * y * y + std::sqrt(21834.66796875) * x * x * x * x * y * y * y * z * z - std::sqrt(4651.171875) * x * x * x * x * y * z * z * z * z + std::sqrt(32.2998046875) * x * x * y * y * y * y * y * y * y + std::sqrt(56976.85546875) * x * x * y * y * y * y * y * z * z - std::sqrt(349354.6875) * x * x * y * y * y * z * z * z * z + std::sqrt(33075.0) * x * x * y * z * z * z * z * z * z - std::sqrt(50.46844482421875) * y * y * y * y * y * y * y * y * y - std::sqrt(129.19921875) * y * y * y * y * y * y * y * z * z - std::sqrt(25323.046875) * y * y * y * y * y * z * z * z * z + std::sqrt(33075.0) * y * y * y * z * z * z * z * z * z) + e_2 * (-std::sqrt(112435.6201171875) * x * x * x * x * x * x * y + std::sqrt(290.6982421875) * x * x * x * x * y * y * y - std::sqrt(562791.796875) * x * x * x * x * y * z * z + std::sqrt(49128.0029296875) * x * x * y * y * y * y * y - std::sqrt(18604.6875) * x * x * y * y * y * z * z + std::sqrt(74418.75) * x * x * y * z * z * z * z - std::sqrt(17086.5966796875) * y * y * y * y * y * y * y - std::sqrt(227907.421875) * y * y * y * y * y * z * z + std::sqrt(74418.75) * y * y * y * z * z * z * z + std::sqrt(132300.0) * y * z * z * z * z * z * z) + e_3 * (-std::sqrt(3646518.75) * x * x * x * x * y + std::sqrt(529200.0) * x * x * y * y * y - std::sqrt(1190700.0) * x * x * y * z * z - std::sqrt(1000518.75) * y * y * y * y * y - std::sqrt(1190700.0) * y * y * y * z * z + std::sqrt(4762800.0) * y * z * z * z * z) + e_4 * (-std::sqrt(6716292.1875) * x * x * y - std::sqrt(6716292.1875) * y * y * y + std::sqrt(4762800.0) * y * z * z) + e_5 * (-std::sqrt(2679075.0) * y);
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

        pc_91[k] = e_0 * (-std::sqrt(5.38330078125) * x * x * x * x * x * x * x * x * x * y * z - std::sqrt(21.533203125) * x * x * x * x * x * x * x * y * y * y * z + std::sqrt(1744.189453125) * x * x * x * x * x * x * x * y * z * z * z + std::sqrt(1744.189453125) * x * x * x * x * x * y * y * y * z * z * z - std::sqrt(12403.125) * x * x * x * x * x * y * z * z * z * z * z + std::sqrt(21.533203125) * x * x * x * y * y * y * y * y * y * y * z - std::sqrt(1744.189453125) * x * x * x * y * y * y * y * y * z * z * z + std::sqrt(5512.5) * x * x * x * y * z * z * z * z * z * z * z + std::sqrt(5.38330078125) * x * y * y * y * y * y * y * y * y * y * z - std::sqrt(1744.189453125) * x * y * y * y * y * y * y * y * z * z * z + std::sqrt(12403.125) * x * y * y * y * y * y * z * z * z * z * z - std::sqrt(5512.5) * x * y * y * y * z * z * z * z * z * z * z) + e_1 * (std::sqrt(775.1953125) * x * x * x * x * x * x * x * y * z + std::sqrt(775.1953125) * x * x * x * x * x * y * y * y * z - std::sqrt(12403.125) * x * x * x * x * x * y * z * z * z - std::sqrt(775.1953125) * x * x * x * y * y * y * y * y * z + std::sqrt(198450.0) * x * x * x * y * z * z * z * z * z - std::sqrt(775.1953125) * x * y * y * y * y * y * y * y * z + std::sqrt(12403.125) * x * y * y * y * y * y * z * z * z - std::sqrt(198450.0) * x * y * y * y * z * z * z * z * z) + e_2 * (std::sqrt(27907.03125) * x * x * x * x * x * y * z + std::sqrt(2790703.125) * x * x * x * y * z * z * z - std::sqrt(27907.03125) * x * y * y * y * y * y * z - std::sqrt(2790703.125) * x * y * y * y * z * z * z) + e_3 * (std::sqrt(4961250.0) * x * x * x * y * z - std::sqrt(4961250.0) * x * y * y * y * z);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        pc_92[k] = e_0 * (std::sqrt(0.048065185546875) * x * x * x * x * x * x * x * x * x * x * y + std::sqrt(0.432586669921875) * x * x * x * x * x * x * x * x * y * y * y - std::sqrt(37.68310546875) * x * x * x * x * x * x * x * x * y * z * z + std::sqrt(0.1922607421875) * x * x * x * x * x * x * y * y * y * y * y - std::sqrt(150.732421875) * x * x * x * x * x * x * y * y * y * z * z + std::sqrt(2242.529296875) * x * x * x * x * x * x * y * z * z * z * z - std::sqrt(0.1922607421875) * x * x * x * x * y * y * y * y * y * y * y + std::sqrt(2242.529296875) * x * x * x * x * y * y * y * z * z * z * z - std::sqrt(4921.875) * x * x * x * x * y * z * z * z * z * z * z - std::sqrt(0.432586669921875) * x * x * y * y * y * y * y * y * y * y * y + std::sqrt(150.732421875) * x * x * y * y * y * y * y * y * y * z * z - std::sqrt(2242.529296875) * x * x * y * y * y * y * y * z * z * z * z + std::sqrt(787.5) * x * x * y * z * z * z * z * z * z * z * z - std::sqrt(0.048065185546875) * y * y * y * y * y * y * y * y * y * y * y + std::sqrt(37.68310546875) * y * y * y * y * y * y * y * y * y * z * z - std::sqrt(2242.529296875) * y * y * y * y * y * y * y * z * z * z * z + std::sqrt(4921.875) * y * y * y * y * y * z * z * z * z * z * z - std::sqrt(787.5) * y * y * y * z * z * z * z * z * z * z * z) + e_1 * (std::sqrt(32.4920654296875) * x * x * x * x * x * x * x * x * y + std::sqrt(110.7421875) * x * x * x * x * x * x * y * y * y + std::sqrt(5687.841796875) * x * x * x * x * x * x * y * z * z - std::sqrt(6.92138671875) * x * x * x * x * y * y * y * y * y + std::sqrt(9994.482421875) * x * x * x * x * y * y * y * z * z - std::sqrt(39977.9296875) * x * x * x * x * y * z * z * z * z - std::sqrt(196.875) * x * x * y * y * y * y * y * y * y - std::sqrt(692.138671875) * x * x * y * y * y * y * y * z * z - std::sqrt(35880.46875) * x * x * y * y * y * z * z * z * z + std::sqrt(133087.5) * x * x * y * z * z * z * z * z * z - std::sqrt(43.2586669921875) * y * y * y * y * y * y * y * y * y - std::sqrt(2587.060546875) * y * y * y * y * y * y * y * z * z + std::sqrt(110.7421875) * y * y * y * y * y * z * z * z * z - std::sqrt(7087.5) * y * y * y * z * z * z * z * z * z - std::sqrt(3150.0) * y * z * z * z * z * z * z * z * z) + e_2 * (std::sqrt(17767.96875) * x * x * x * x * x * x * y + std::sqrt(13399.8046875) * x * x * x * x * y * y * y + std::sqrt(13399.8046875) * x * x * x * x * y * z * z - std::sqrt(28350.0) * x * x * y * y * y * y * y - std::sqrt(159911.71875) * x * x * y * y * y * z * z + std::sqrt(3430350.0) * x * x * y * z * z * z * z - std::sqrt(22751.3671875) * y * y * y * y * y * y * y - std::sqrt(265891.9921875) * y * y * y * y * y * z * z - std::sqrt(453600.0) * y * y * y * z * z * z * z - std::sqrt(532350.0) * y * z * z * z * z * z * z) + e_3 * (std::sqrt(673755.46875) * x * x * x * x * y - std::sqrt(214396.875) * x * x * y * y * y + std::sqrt(12502350.0) * x * x * y * z * z - std::sqrt(1648286.71875) * y * y * y * y * y - std::sqrt(10234350.0) * y * y * y * z * z - std::sqrt(13721400.0) * y * z * z * z * z) + e_4 * (std::sqrt(3125587.5) * x * x * y - std::sqrt(17017087.5) * y * y * y - std::sqrt(50009400.0) * y * z * z) + e_5 * (-std::sqrt(12502350.0) * y);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        pc_93[k] = e_0 * (std::sqrt(0.720977783203125) * x * x * x * x * x * x * x * x * x * x * z + std::sqrt(6.488800048828125) * x * x * x * x * x * x * x * x * y * y * z - std::sqrt(251.220703125) * x * x * x * x * x * x * x * x * z * z * z + std::sqrt(2.8839111328125) * x * x * x * x * x * x * y * y * y * y * z - std::sqrt(1004.8828125) * x * x * x * x * x * x * y * y * z * z * z + std::sqrt(2526.767578125) * x * x * x * x * x * x * z * z * z * z * z - std::sqrt(2.8839111328125) * x * x * x * x * y * y * y * y * y * y * z + std::sqrt(2526.767578125) * x * x * x * x * y * y * z * z * z * z * z - std::sqrt(1890.0) * x * x * x * x * z * z * z * z * z * z * z - std::sqrt(6.488800048828125) * x * x * y * y * y * y * y * y * y * y * z + std::sqrt(1004.8828125) * x * x * y * y * y * y * y * y * z * z * z - std::sqrt(2526.767578125) * x * x * y * y * y * y * z * z * z * z * z + std::sqrt(52.5) * x * x * z * z * z * z * z * z * z * z * z - std::sqrt(0.720977783203125) * y * y * y * y * y * y * y * y * y * y * z + std::sqrt(251.220703125) * y * y * y * y * y * y * y * y * z * z * z - std::sqrt(2526.767578125) * y * y * y * y * y * y * z * z * z * z * z + std::sqrt(1890.0) * y * y * y * y * z * z * z * z * z * z * z - std::sqrt(52.5) * y * y * z * z * z * z * z * z * z * z * z) + e_1 * (-std::sqrt(46.142578125) * x * x * x * x * x * x * x * x * z - std::sqrt(184.5703125) * x * x * x * x * x * x * y * y * z + std::sqrt(738.28125) * x * x * x * x * x * x * z * z * z + std::sqrt(738.28125) * x * x * x * x * y * y * z * z * z - std::sqrt(11812.5) * x * x * x * x * z * z * z * z * z + std::sqrt(184.5703125) * x * x * y * y * y * y * y * y * z - std::sqrt(738.28125) * x * x * y * y * y * y * z * z * z + std::sqrt(46.142578125) * y * y * y * y * y * y * y * y * z - std::sqrt(738.28125) * y * y * y * y * y * y * z * z * z + std::sqrt(11812.5) * y * y * y * y * z * z * z * z * z) + e_2 * (-std::sqrt(3737.548828125) * x * x * x * x * x * x * z - std::sqrt(3737.548828125) * x * x * x * x * y * y * z - std::sqrt(106312.5) * x * x * x * x * z * z * z + std::sqrt(3737.548828125) * x * x * y * y * y * y * z - std::sqrt(106312.5) * x * x * z * z * z * z * z + std::sqrt(3737.548828125) * y * y * y * y * y * y * z + std::sqrt(106312.5) * y * y * y * y * z * z * z + std::sqrt(106312.5) * y * y * z * z * z * z * z) + e_3 * (-std::sqrt(425250.0) * x * x * x * x * z - std::sqrt(3024000.0) * x * x * z * z * z + std::sqrt(425250.0) * y * y * y * y * z + std::sqrt(3024000.0) * y * y * z * z * z) + e_4 * (-std::sqrt(5209312.5) * x * x * z + std::sqrt(5209312.5) * y * y * z);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        pc_94[k] = e_0 * (std::sqrt(0.048065185546875) * x * x * x * x * x * x * x * x * x * x * x + std::sqrt(0.432586669921875) * x * x * x * x * x * x * x * x * x * y * y - std::sqrt(37.68310546875) * x * x * x * x * x * x * x * x * x * z * z + std::sqrt(0.1922607421875) * x * x * x * x * x * x * x * y * y * y * y - std::sqrt(150.732421875) * x * x * x * x * x * x * x * y * y * z * z + std::sqrt(2242.529296875) * x * x * x * x * x * x * x * z * z * z * z - std::sqrt(0.1922607421875) * x * x * x * x * x * y * y * y * y * y * y + std::sqrt(2242.529296875) * x * x * x * x * x * y * y * z * z * z * z - std::sqrt(4921.875) * x * x * x * x * x * z * z * z * z * z * z - std::sqrt(0.432586669921875) * x * x * x * y * y * y * y * y * y * y * y + std::sqrt(150.732421875) * x * x * x * y * y * y * y * y * y * z * z - std::sqrt(2242.529296875) * x * x * x * y * y * y * y * z * z * z * z + std::sqrt(787.5) * x * x * x * z * z * z * z * z * z * z * z - std::sqrt(0.048065185546875) * x * y * y * y * y * y * y * y * y * y * y + std::sqrt(37.68310546875) * x * y * y * y * y * y * y * y * y * z * z - std::sqrt(2242.529296875) * x * y * y * y * y * y * y * z * z * z * z + std::sqrt(4921.875) * x * y * y * y * y * z * z * z * z * z * z - std::sqrt(787.5) * x * y * y * z * z * z * z * z * z * z * z) + e_1 * (std::sqrt(43.2586669921875) * x * x * x * x * x * x * x * x * x + std::sqrt(196.875) * x * x * x * x * x * x * x * y * y + std::sqrt(2587.060546875) * x * x * x * x * x * x * x * z * z + std::sqrt(6.92138671875) * x * x * x * x * x * y * y * y * y + std::sqrt(692.138671875) * x * x * x * x * x * y * y * z * z - std::sqrt(110.7421875) * x * x * x * x * x * z * z * z * z - std::sqrt(110.7421875) * x * x * x * y * y * y * y * y * y - std::sqrt(9994.482421875) * x * x * x * y * y * y * y * z * z + std::sqrt(35880.46875) * x * x * x * y * y * z * z * z * z + std::sqrt(7087.5) * x * x * x * z * z * z * z * z * z - std::sqrt(32.4920654296875) * x * y * y * y * y * y * y * y * y - std::sqrt(5687.841796875) * x * y * y * y * y * y * y * z * z + std::sqrt(39977.9296875) * x * y * y * y * y * z * z * z * z - std::sqrt(133087.5) * x * y * y * z * z * z * z * z * z + std::sqrt(3150.0) * x * z * z * z * z * z * z * z * z) + e_2 * (std::sqrt(22751.3671875) * x * x * x * x * x * x * x + std::sqrt(28350.0) * x * x * x * x * x * y * y + std::sqrt(265891.9921875) * x * x * x * x * x * z * z - std::sqrt(13399.8046875) * x * x * x * y * y * y * y + std::sqrt(159911.71875) * x * x * x * y * y * z * z + std::sqrt(453600.0) * x * x * x * z * z * z * z - std::sqrt(17767.96875) * x * y * y * y * y * y * y - std::sqrt(13399.8046875) * x * y * y * y * y * z * z - std::sqrt(3430350.0) * x * y * y * z * z * z * z + std::sqrt(532350.0) * x * z * z * z * z * z * z) + e_3 * (std::sqrt(1648286.71875) * x * x * x * x * x + std::sqrt(214396.875) * x * x * x * y * y + std::sqrt(10234350.0) * x * x * x * z * z - std::sqrt(673755.46875) * x * y * y * y * y - std::sqrt(12502350.0) * x * y * y * z * z + std::sqrt(13721400.0) * x * z * z * z * z) + e_4 * (std::sqrt(17017087.5) * x * x * x - std::sqrt(3125587.5) * x * y * y + std::sqrt(50009400.0) * x * z * z) + e_5 * (std::sqrt(12502350.0) * x);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        pc_95[k] = e_0 * (-std::sqrt(1.3458251953125) * x * x * x * x * x * x * x * x * x * x * z - std::sqrt(1.3458251953125) * x * x * x * x * x * x * x * x * y * y * z + std::sqrt(436.04736328125) * x * x * x * x * x * x * x * x * z * z * z + std::sqrt(5.38330078125) * x * x * x * x * x * x * y * y * y * y * z - std::sqrt(3100.78125) * x * x * x * x * x * x * z * z * z * z * z + std::sqrt(5.38330078125) * x * x * x * x * y * y * y * y * y * y * z - std::sqrt(1744.189453125) * x * x * x * x * y * y * y * y * z * z * z + std::sqrt(3100.78125) * x * x * x * x * y * y * z * z * z * z * z + std::sqrt(1378.125) * x * x * x * x * z * z * z * z * z * z * z - std::sqrt(1.3458251953125) * x * x * y * y * y * y * y * y * y * y * z + std::sqrt(3100.78125) * x * x * y * y * y * y * z * z * z * z * z - std::sqrt(5512.5) * x * x * y * y * z * z * z * z * z * z * z - std::sqrt(1.3458251953125) * y * y * y * y * y * y * y * y * y * y * z + std::sqrt(436.04736328125) * y * y * y * y * y * y * y * y * z * z * z - std::sqrt(3100.78125) * y * y * y * y * y * y * z * z * z * z * z + std::sqrt(1378.125) * y * y * y * y * z * z * z * z * z * z * z) + e_1 * (std::sqrt(86.1328125) * x * x * x * x * x * x * x * x * z - std::sqrt(344.53125) * x * x * x * x * x * x * y * y * z + std::sqrt(775.1953125) * x * x * x * x * x * x * z * z * z - std::sqrt(3100.78125) * x * x * x * x * y * y * y * y * z + std::sqrt(93798.6328125) * x * x * x * x * y * y * z * z * z - std::sqrt(344.53125) * x * x * y * y * y * y * y * y * z + std::sqrt(93798.6328125) * x * x * y * y * y * y * z * z * z - std::sqrt(793800.0) * x * x * y * y * z * z * z * z * z + std::sqrt(22050.0) * x * x * z * z * z * z * z * z * z + std::sqrt(86.1328125) * y * y * y * y * y * y * y * y * z + std::sqrt(775.1953125) * y * y * y * y * y * y * z * z * z + std::sqrt(22050.0) * y * y * z * z * z * z * z * z * z) + e_2 * (std::sqrt(26378.173828125) * x * x * x * x * x * x * z + std::sqrt(23449.658203125) * x * x * x * x * y * y * z + std::sqrt(131008.0078125) * x * x * x * x * z * z * z + std::sqrt(23449.658203125) * x * x * y * y * y * y * z - std::sqrt(6849625.78125) * x * x * y * y * z * z * z + std::sqrt(1240312.5) * x * x * z * z * z * z * z + std::sqrt(26378.173828125) * y * y * y * y * y * y * z + std::sqrt(131008.0078125) * y * y * y * y * z * z * z + std::sqrt(1240312.5) * y * y * z * z * z * z * z + std::sqrt(22050.0) * z * z * z * z * z * z * z) + e_3 * (std::sqrt(1500778.125) * x * x * x * x * z - std::sqrt(4018612.5) * x * x * y * y * z + std::sqrt(12700800.0) * x * x * z * z * z + std::sqrt(1500778.125) * y * y * y * y * z + std::sqrt(12700800.0) * y * y * z * z * z + std::sqrt(3175200.0) * z * z * z * z * z) + e_4 * (std::sqrt(16074450.0) * x * x * z + std::sqrt(16074450.0) * y * y * z + std::sqrt(38896200.0) * z * z * z) + e_5 * (std::sqrt(28576800.0) * z);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        pc_96[k] = e_0 * (-std::sqrt(0.0560760498046875) * x * x * x * x * x * x * x * x * x * x * x + std::sqrt(0.0560760498046875) * x * x * x * x * x * x * x * x * x * y * y + std::sqrt(32.2998046875) * x * x * x * x * x * x * x * x * x * z * z + std::sqrt(2.01873779296875) * x * x * x * x * x * x * x * y * y * y * y - std::sqrt(129.19921875) * x * x * x * x * x * x * x * y * y * z * z - std::sqrt(1162.79296875) * x * x * x * x * x * x * x * z * z * z * z + std::sqrt(0.22430419921875) * x * x * x * x * x * y * y * y * y * y * y - std::sqrt(516.796875) * x * x * x * x * x * y * y * y * y * z * z + std::sqrt(10465.13671875) * x * x * x * x * x * y * y * z * z * z * z + std::sqrt(918.75) * x * x * x * x * x * z * z * z * z * z * z - std::sqrt(1.4019012451171875) * x * x * x * y * y * y * y * y * y * y * y + std::sqrt(129.19921875) * x * x * x * y * y * y * y * y * y * z * z + std::sqrt(1162.79296875) * x * x * x * y * y * y * y * z * z * z * z - std::sqrt(14700.0) * x * x * x * y * y * z * z * z * z * z * z - std::sqrt(0.5046844482421875) * x * y * y * y * y * y * y * y * y * y * y + std::sqrt(290.6982421875) * x * y * y * y * y * y * y * y * y * z * z - std::sqrt(10465.13671875) * x * y * y * y * y * y * y * z * z * z * z + std::sqrt(8268.75) * x * y * y * y * y * z * z * z * z * z * z) + e_1 * (-std::sqrt(50.46844482421875) * x * x * x * x * x * x * x * x * x + std::sqrt(32.2998046875) * x * x * x * x * x * x * x * y * y - std::sqrt(129.19921875) * x * x * x * x * x * x * x * z * z + std::sqrt(201.873779296875) * x * x * x * x * x * y * y * y * y + std::sqrt(56976.85546875) * x * x * x * x * x * y * y * z * z - std::sqrt(25323.046875) * x * x * x * x * x * z * z * z * z - std::sqrt(290.6982421875) * x * x * x * y * y * y * y * y * y + std::sqrt(21834.66796875) * x * x * x * y * y * y * y * z * z - std::sqrt(349354.6875) * x * x * x * y * y * z * z * z * z + std::sqrt(33075.0) * x * x * x * z * z * z * z * z * z - std::sqrt(341.16668701171875) * x * y * y * y * y * y * y * y * y - std::sqrt(10465.13671875) * x * y * y * y * y * y * y * z * z - std::sqrt(4651.171875) * x * y * y * y * y * z * z * z * z + std::sqrt(33075.0) * x * y * y * z * z * z * z * z * z) + e_2 * (-std::sqrt(17086.5966796875) * x * x * x * x * x * x * x + std::sqrt(49128.0029296875) * x * x * x * x * x * y * y - std::sqrt(227907.421875) * x * x * x * x * x * z * z + std::sqrt(290.6982421875) * x * x * x * y * y * y * y - std::sqrt(18604.6875) * x * x * x * y * y * z * z + std::sqrt(74418.75) * x * x * x * z * z * z * z - std::sqrt(112435.6201171875) * x * y * y * y * y * y * y - std::sqrt(562791.796875) * x * y * y * y * y * z * z + std::sqrt(74418.75) * x * y * y * z * z * z * z + std::sqrt(132300.0) * x * z * z * z * z * z * z) + e_3 * (-std::sqrt(1000518.75) * x * x * x * x * x + std::sqrt(529200.0) * x * x * x * y * y - std::sqrt(1190700.0) * x * x * x * z * z - std::sqrt(3646518.75) * x * y * y * y * y - std::sqrt(1190700.0) * x * y * y * z * z + std::sqrt(4762800.0) * x * z * z * z * z) + e_4 * (-std::sqrt(6716292.1875) * x * x * x - std::sqrt(6716292.1875) * x * y * y + std::sqrt(4762800.0) * x * z * z) + e_5 * (-std::sqrt(2679075.0) * x);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        pc_97[k] = e_0 * (std::sqrt(1.009368896484375) * x * x * x * x * x * x * x * x * x * x * z - std::sqrt(25.234222412109375) * x * x * x * x * x * x * x * x * y * y * z - std::sqrt(258.3984375) * x * x * x * x * x * x * x * x * z * z * z - std::sqrt(36.3372802734375) * x * x * x * x * x * x * y * y * y * y * z + std::sqrt(9302.34375) * x * x * x * x * x * x * y * y * z * z * z + std::sqrt(258.3984375) * x * x * x * x * x * x * z * z * z * z * z + std::sqrt(36.3372802734375) * x * x * x * x * y * y * y * y * y * y * z - std::sqrt(12661.5234375) * x * x * x * x * y * y * z * z * z * z * z + std::sqrt(25.234222412109375) * x * x * y * y * y * y * y * y * y * y * z - std::sqrt(9302.34375) * x * x * y * y * y * y * y * y * z * z * z + std::sqrt(12661.5234375) * x * x * y * y * y * y * z * z * z * z * z - std::sqrt(1.009368896484375) * y * y * y * y * y * y * y * y * y * y * z + std::sqrt(258.3984375) * y * y * y * y * y * y * y * y * z * z * z - std::sqrt(258.3984375) * y * y * y * y * y * y * z * z * z * z * z) + e_1 * (-std::sqrt(64.599609375) * x * x * x * x * x * x * x * x * z + std::sqrt(12661.5234375) * x * x * x * x * x * x * y * y * z - std::sqrt(37209.375) * x * x * x * x * x * x * z * z * z + std::sqrt(103359.375) * x * x * x * x * y * y * z * z * z + std::sqrt(16537.5) * x * x * x * x * z * z * z * z * z - std::sqrt(12661.5234375) * x * x * y * y * y * y * y * y * z - std::sqrt(103359.375) * x * x * y * y * y * y * z * z * z + std::sqrt(64.599609375) * y * y * y * y * y * y * y * y * z + std::sqrt(37209.375) * y * y * y * y * y * y * z * z * z - std::sqrt(16537.5) * y * y * y * y * z * z * z * z * z) + e_2 * (-std::sqrt(119444.677734375) * x * x * x * x * x * x * z + std::sqrt(1758724.365234375) * x * x * x * x * y * y * z - std::sqrt(413437.5) * x * x * x * x * z * z * z - std::sqrt(1758724.365234375) * x * x * y * y * y * y * z + std::sqrt(148837.5) * x * x * z * z * z * z * z + std::sqrt(119444.677734375) * y * y * y * y * y * y * z + std::sqrt(413437.5) * y * y * y * y * z * z * z - std::sqrt(148837.5) * y * y * z * z * z * z * z) + e_3 * (-std::sqrt(3720937.5) * x * x * x * x * z + std::sqrt(3720937.5) * y * y * y * y * z) + e_4 * (-std::sqrt(8372109.375) * x * x * z + std::sqrt(8372109.375) * y * y * z);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        pc_98[k] = e_0 * (std::sqrt(0.1009368896484375) * x * x * x * x * x * x * x * x * x * x * x - std::sqrt(8.175888061523438) * x * x * x * x * x * x * x * x * x * y * y - std::sqrt(25.83984375) * x * x * x * x * x * x * x * x * x * z * z - std::sqrt(3.63372802734375) * x * x * x * x * x * x * x * y * y * y * y + std::sqrt(2583.984375) * x * x * x * x * x * x * x * y * y * z * z + std::sqrt(25.83984375) * x * x * x * x * x * x * x * z * z * z * z + std::sqrt(19.78363037109375) * x * x * x * x * x * y * y * y * y * y * y - std::sqrt(413.4375) * x * x * x * x * x * y * y * y * y * z * z - std::sqrt(3126.62109375) * x * x * x * x * x * y * y * z * z * z * z + std::sqrt(2.5234222412109375) * x * x * x * y * y * y * y * y * y * y * y - std::sqrt(2583.984375) * x * x * x * y * y * y * y * y * y * z * z + std::sqrt(5813.96484375) * x * x * x * y * y * y * y * z * z * z * z - std::sqrt(2.5234222412109375) * x * y * y * y * y * y * y * y * y * y * y + std::sqrt(645.99609375) * x * y * y * y * y * y * y * y * y * z * z - std::sqrt(645.99609375) * x * y * y * y * y * y * y * z * z * z * z) + e_1 * (std::sqrt(90.84320068359375) * x * x * x * x * x * x * x * x * x - std::sqrt(4037.4755859375) * x * x * x * x * x * x * x * y * y - std::sqrt(10335.9375) * x * x * x * x * x * x * x * z * z + std::sqrt(40.374755859375) * x * x * x * x * x * y * y * y * y + std::sqrt(372093.75) * x * x * x * x * x * y * y * z * z + std::sqrt(2583.984375) * x * x * x * x * x * z * z * z * z + std::sqrt(1453.4912109375) * x * x * x * y * y * y * y * y * y - std::sqrt(258398.4375) * x * x * x * y * y * y * y * z * z - std::sqrt(10335.9375) * x * x * x * y * y * z * z * z * z - std::sqrt(1705.8334350585938) * x * y * y * y * y * y * y * y * y + std::sqrt(165375.0) * x * y * y * y * y * y * y * z * z - std::sqrt(23255.859375) * x * y * y * y * y * z * z * z * z) + e_2 * (std::sqrt(7913.4521484375) * x * x * x * x * x * x * x - std::sqrt(117732.7880859375) * x * x * x * x * x * y * y - std::sqrt(372093.75) * x * x * x * x * x * z * z + std::sqrt(4037.4755859375) * x * x * x * y * y * y * y + std::sqrt(1488375.0) * x * x * x * y * y * z * z + std::sqrt(41343.75) * x * x * x * z * z * z * z - std::sqrt(100936.8896484375) * x * y * y * y * y * y * y + std::sqrt(3348843.75) * x * y * y * y * y * z * z - std::sqrt(372093.75) * x * y * y * z * z * z * z) + e_3 * (std::sqrt(93023.4375) * x * x * x * x * x - std::sqrt(372093.75) * x * x * x * y * y - std::sqrt(1488375.0) * x * x * x * z * z - std::sqrt(837210.9375) * x * y * y * y * y + std::sqrt(13395375.0) * x * y * y * z * z) + e_4 * (std::sqrt(93023.4375) * x * x * x - std::sqrt(837210.9375) * x * y * y);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        pc_99[k] = e_0 * (-std::sqrt(90.84320068359375) * x * x * x * x * x * x * x * x * x * y * z + std::sqrt(1453.4912109375) * x * x * x * x * x * x * x * y * y * y * z + std::sqrt(645.99609375) * x * x * x * x * x * x * x * y * z * z * z - std::sqrt(130.814208984375) * x * x * x * x * x * y * y * y * y * y * z - std::sqrt(16149.90234375) * x * x * x * x * x * y * y * y * z * z * z - std::sqrt(2848.8427734375) * x * x * x * y * y * y * y * y * y * y * z + std::sqrt(24832.08984375) * x * x * x * y * y * y * y * y * z * z * z + std::sqrt(32.70355224609375) * x * y * y * y * y * y * y * y * y * y * z - std::sqrt(232.55859375) * x * y * y * y * y * y * y * y * z * z * z) + e_1 * (-std::sqrt(23255.859375) * x * x * x * x * x * x * x * y * z + std::sqrt(93023.4375) * x * x * x * x * x * y * y * y * z + std::sqrt(23255.859375) * x * x * x * x * x * y * z * z * z - std::sqrt(581396.484375) * x * x * x * y * y * y * y * y * z + std::sqrt(93023.4375) * x * x * x * y * y * y * z * z * z + std::sqrt(23255.859375) * x * y * y * y * y * y * z * z * z) + e_2 * (-std::sqrt(837210.9375) * x * x * x * x * x * y * z - std::sqrt(3348843.75) * x * x * x * y * y * y * z + std::sqrt(1488375.0) * x * x * x * y * z * z * z - std::sqrt(837210.9375) * x * y * y * y * y * y * z + std::sqrt(1488375.0) * x * y * y * y * z * z * z) + e_3 * (-std::sqrt(13395375.0) * x * x * x * y * z - std::sqrt(13395375.0) * x * y * y * y * z + std::sqrt(5953500.0) * x * y * z * z * z) + e_4 * (-std::sqrt(13395375.0) * x * y * z);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        pc_100[k] = e_0 * (-std::sqrt(581.396484375) * x * x * x * x * x * x * x * x * y * z * z + std::sqrt(5232.568359375) * x * x * x * x * x * x * y * y * y * z * z + std::sqrt(4134.375) * x * x * x * x * x * x * y * z * z * z * z + std::sqrt(581.396484375) * x * x * x * x * y * y * y * y * y * z * z - std::sqrt(66150.0) * x * x * x * x * y * y * y * z * z * z * z - std::sqrt(5232.568359375) * x * x * y * y * y * y * y * y * y * z * z + std::sqrt(37209.375) * x * x * y * y * y * y * y * z * z * z * z) + e_1 * (-std::sqrt(581.396484375) * x * x * x * x * x * x * x * x * y + std::sqrt(5232.568359375) * x * x * x * x * x * x * y * y * y - std::sqrt(5232.568359375) * x * x * x * x * x * x * y * z * z + std::sqrt(581.396484375) * x * x * x * x * y * y * y * y * y - std::sqrt(47093.115234375) * x * x * x * x * y * y * y * z * z + std::sqrt(37209.375) * x * x * x * x * y * z * z * z * z - std::sqrt(5232.568359375) * x * x * y * y * y * y * y * y * y - std::sqrt(47093.115234375) * x * x * y * y * y * y * y * z * z + std::sqrt(148837.5) * x * x * y * y * y * z * z * z * z - std::sqrt(5232.568359375) * y * y * y * y * y * y * y * z * z + std::sqrt(37209.375) * y * y * y * y * y * z * z * z * z) + e_2 * (-std::sqrt(70348.974609375) * x * x * x * x * x * x * y + std::sqrt(307558.740234375) * x * x * x * x * y * y * y - std::sqrt(83721.09375) * x * x * x * x * y * z * z - std::sqrt(633140.771484375) * x * x * y * y * y * y * y - std::sqrt(334884.375) * x * x * y * y * y * z * z + std::sqrt(1339537.5) * x * x * y * z * z * z * z - std::sqrt(5232.568359375) * y * y * y * y * y * y * y - std::sqrt(83721.09375) * y * y * y * y * y * z * z + std::sqrt(1339537.5) * y * y * y * z * z * z * z) + e_3 * (-std::sqrt(753489.84375) * x * x * x * x * y - std::sqrt(3013959.375) * x * x * y * y * y + std::sqrt(1339537.5) * x * x * y * z * z - std::sqrt(753489.84375) * y * y * y * y * y + std::sqrt(1339537.5) * y * y * y * z * z + std::sqrt(2381400.0) * y * z * z * z * z) + e_4 * (-std::sqrt(5358150.0) * x * x * y - std::sqrt(5358150.0) * y * y * y + std::sqrt(12055837.5) * y * z * z) + e_5 * (-std::sqrt(1339537.5) * y);
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

        pc_101[k] = e_0 * (std::sqrt(18.16864013671875) * x * x * x * x * x * x * x * x * x * y * z - std::sqrt(32.2998046875) * x * x * x * x * x * x * x * y * y * y * z - std::sqrt(2067.1875) * x * x * x * x * x * x * x * y * z * z * z - std::sqrt(395.672607421875) * x * x * x * x * x * y * y * y * y * y * z + std::sqrt(11254.6875) * x * x * x * x * x * y * y * y * z * z * z + std::sqrt(8268.75) * x * x * x * x * x * y * z * z * z * z * z - std::sqrt(32.2998046875) * x * x * x * y * y * y * y * y * y * y * z + std::sqrt(11254.6875) * x * x * x * y * y * y * y * y * z * z * z - std::sqrt(91875.0) * x * x * x * y * y * y * z * z * z * z * z + std::sqrt(18.16864013671875) * x * y * y * y * y * y * y * y * y * y * z - std::sqrt(2067.1875) * x * y * y * y * y * y * y * y * z * z * z + std::sqrt(8268.75) * x * y * y * y * y * y * z * z * z * z * z) + e_1 * (std::sqrt(74418.75) * x * x * x * x * x * y * z * z * z - std::sqrt(826875.0) * x * x * x * y * y * y * z * z * z + std::sqrt(74418.75) * x * y * y * y * y * y * z * z * z) + e_2 * (std::sqrt(167442.1875) * x * x * x * x * x * y * z - std::sqrt(1860468.75) * x * x * x * y * y * y * z + std::sqrt(167442.1875) * x * y * y * y * y * y * z);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        pc_102[k] = e_0 * (std::sqrt(193.798828125) * x * x * x * x * x * x * x * x * y * z * z - std::sqrt(193.798828125) * x * x * x * x * x * x * y * y * y * z * z - std::sqrt(4220.5078125) * x * x * x * x * x * x * y * z * z * z * z - std::sqrt(4844.970703125) * x * x * x * x * y * y * y * y * y * z * z + std::sqrt(16882.03125) * x * x * x * x * y * y * y * z * z * z * z + std::sqrt(5512.5) * x * x * x * x * y * z * z * z * z * z * z - std::sqrt(1744.189453125) * x * x * y * y * y * y * y * y * y * z * z + std::sqrt(37984.5703125) * x * x * y * y * y * y * y * z * z * z * z - std::sqrt(49612.5) * x * x * y * y * y * z * z * z * z * z * z) + e_1 * (std::sqrt(193.798828125) * x * x * x * x * x * x * x * x * y - std::sqrt(193.798828125) * x * x * x * x * x * x * y * y * y - std::sqrt(1744.189453125) * x * x * x * x * x * x * y * z * z - std::sqrt(4844.970703125) * x * x * x * x * y * y * y * y * y - std::sqrt(15697.705078125) * x * x * x * x * y * y * y * z * z + std::sqrt(279845.5078125) * x * x * x * x * y * z * z * z * z - std::sqrt(1744.189453125) * x * x * y * y * y * y * y * y * y - std::sqrt(15697.705078125) * x * x * y * y * y * y * y * z * z - std::sqrt(375194.53125) * x * x * y * y * y * z * z * z * z - std::sqrt(49612.5) * x * x * y * z * z * z * z * z * z - std::sqrt(1744.189453125) * y * y * y * y * y * y * y * z * z + std::sqrt(37984.5703125) * y * y * y * y * y * z * z * z * z - std::sqrt(49612.5) * y * y * y * z * z * z * z * z * z) + e_2 * (std::sqrt(23449.658203125) * x * x * x * x * x * x * y - std::sqrt(265310.595703125) * x * x * x * x * y * y * y + std::sqrt(1179072.0703125) * x * x * x * x * y * z * z - std::sqrt(504070.751953125) * x * x * y * y * y * y * y - std::sqrt(8065132.03125) * x * x * y * y * y * z * z - std::sqrt(1004653.125) * x * x * y * z * z * z * z - std::sqrt(1744.189453125) * y * y * y * y * y * y * y + std::sqrt(6976.7578125) * y * y * y * y * y * z * z - std::sqrt(1004653.125) * y * y * y * z * z * z * z - std::sqrt(198450.0) * y * z * z * z * z * z * z) + e_3 * (std::sqrt(375194.53125) * x * x * x * x * y - std::sqrt(18865153.125) * x * x * y * y * y - std::sqrt(7144200.0) * x * x * y * z * z - std::sqrt(251163.28125) * y * y * y * y * y - std::sqrt(7144200.0) * y * y * y * z * z - std::sqrt(12700800.0) * y * z * z * z * z) + e_4 * (-std::sqrt(9041878.125) * x * x * y - std::sqrt(9041878.125) * y * y * y - std::sqrt(64297800.0) * y * z * z) + e_5 * (-std::sqrt(16074450.0) * y);
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

        pc_103[k] = e_0 * (-std::sqrt(1.7303466796875) * x * x * x * x * x * x * x * x * x * y * z + std::sqrt(372.216796875) * x * x * x * x * x * x * x * y * z * z * z + std::sqrt(62.29248046875) * x * x * x * x * x * y * y * y * y * y * z - std::sqrt(372.216796875) * x * x * x * x * x * y * y * y * z * z * z - std::sqrt(2768.5546875) * x * x * x * x * x * y * z * z * z * z * z + std::sqrt(110.7421875) * x * x * x * y * y * y * y * y * y * y * z - std::sqrt(9305.419921875) * x * x * x * y * y * y * y * y * z * z * z + std::sqrt(11074.21875) * x * x * x * y * y * y * z * z * z * z * z + std::sqrt(787.5) * x * x * x * y * z * z * z * z * z * z * z + std::sqrt(15.5731201171875) * x * y * y * y * y * y * y * y * y * y * z - std::sqrt(3349.951171875) * x * y * y * y * y * y * y * y * z * z * z + std::sqrt(24916.9921875) * x * y * y * y * y * y * z * z * z * z * z - std::sqrt(7087.5) * x * y * y * y * z * z * z * z * z * z * z) + e_1 * (std::sqrt(110.7421875) * x * x * x * x * x * x * x * y * z + std::sqrt(442.96875) * x * x * x * x * x * y * y * y * z - std::sqrt(32004.4921875) * x * x * x * x * x * y * z * z * z + std::sqrt(110.7421875) * x * x * x * y * y * y * y * y * z - std::sqrt(11074.21875) * x * x * x * y * y * y * z * z * z + std::sqrt(143521.875) * x * x * x * y * z * z * z * z * z + std::sqrt(5426.3671875) * x * y * y * y * y * y * z * z * z + std::sqrt(15946.875) * x * y * y * y * z * z * z * z * z - std::sqrt(28350.0) * x * y * z * z * z * z * z * z * z) + e_2 * (-std::sqrt(15946.875) * x * x * x * x * x * y * z + std::sqrt(708750.0) * x * x * x * y * z * z * z + std::sqrt(15946.875) * x * y * y * y * y * y * z + std::sqrt(708750.0) * x * y * y * y * z * z * z - std::sqrt(1020600.0) * x * y * z * z * z * z * z) + e_3 * (std::sqrt(177187.5) * x * x * x * y * z + std::sqrt(1594687.5) * x * y * y * y * z - std::sqrt(2835000.0) * x * y * z * z * z);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        pc_104[k] = e_0 * (-std::sqrt(25.9552001953125) * x * x * x * x * x * x * x * x * x * z * z + std::sqrt(738.28125) * x * x * x * x * x * x * x * z * z * z * z + std::sqrt(934.38720703125) * x * x * x * x * x * y * y * y * y * z * z - std::sqrt(738.28125) * x * x * x * x * x * y * y * z * z * z * z - std::sqrt(1516.7578125) * x * x * x * x * x * z * z * z * z * z * z + std::sqrt(1661.1328125) * x * x * x * y * y * y * y * y * y * z * z - std::sqrt(18457.03125) * x * x * x * y * y * y * y * z * z * z * z + std::sqrt(6067.03125) * x * x * x * y * y * z * z * z * z * z * z + std::sqrt(52.5) * x * x * x * z * z * z * z * z * z * z * z + std::sqrt(233.5968017578125) * x * y * y * y * y * y * y * y * y * z * z - std::sqrt(6644.53125) * x * y * y * y * y * y * y * z * z * z * z + std::sqrt(13650.8203125) * x * y * y * y * y * z * z * z * z * z * z - std::sqrt(472.5) * x * y * y * z * z * z * z * z * z * z * z) + e_1 * (-std::sqrt(25.9552001953125) * x * x * x * x * x * x * x * x * x - std::sqrt(415.283203125) * x * x * x * x * x * x * x * z * z + std::sqrt(934.38720703125) * x * x * x * x * x * y * y * y * y + std::sqrt(415.283203125) * x * x * x * x * x * y * y * z * z - std::sqrt(1661.1328125) * x * x * x * x * x * z * z * z * z + std::sqrt(1661.1328125) * x * x * x * y * y * y * y * y * y + std::sqrt(10382.080078125) * x * x * x * y * y * y * y * z * z + std::sqrt(6644.53125) * x * x * x * y * y * z * z * z * z - std::sqrt(11812.5) * x * x * x * z * z * z * z * z * z + std::sqrt(233.5968017578125) * x * y * y * y * y * y * y * y * y + std::sqrt(3737.548828125) * x * y * y * y * y * y * y * z * z + std::sqrt(14950.1953125) * x * y * y * y * y * z * z * z * z + std::sqrt(106312.5) * x * y * y * z * z * z * z * z * z) + e_2 * (-std::sqrt(10382.080078125) * x * x * x * x * x * x * x + std::sqrt(10382.080078125) * x * x * x * x * x * y * y - std::sqrt(106312.5) * x * x * x * x * x * z * z + std::sqrt(259552.001953125) * x * x * x * y * y * y * y + std::sqrt(425250.0) * x * x * x * y * y * z * z - std::sqrt(956812.5) * x * x * x * z * z * z * z + std::sqrt(93438.720703125) * x * y * y * y * y * y * y + std::sqrt(956812.5) * x * y * y * y * y * z * z + std::sqrt(8611312.5) * x * y * y * z * z * z * z) + e_3 * (-std::sqrt(620894.53125) * x * x * x * x * x + std::sqrt(2483578.125) * x * x * x * y * y - std::sqrt(7985250.0) * x * x * x * z * z + std::sqrt(5588050.78125) * x * y * y * y * y + std::sqrt(71867250.0) * x * y * y * z * z) + e_4 * (-std::sqrt(5209312.5) * x * x * x + std::sqrt(46883812.5) * x * y * y);
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

        pc_105[k] = e_0 * (-std::sqrt(1.7303466796875) * x * x * x * x * x * x * x * x * x * x * z + std::sqrt(372.216796875) * x * x * x * x * x * x * x * x * z * z * z + std::sqrt(62.29248046875) * x * x * x * x * x * x * y * y * y * y * z - std::sqrt(372.216796875) * x * x * x * x * x * x * y * y * z * z * z - std::sqrt(2768.5546875) * x * x * x * x * x * x * z * z * z * z * z + std::sqrt(110.7421875) * x * x * x * x * y * y * y * y * y * y * z - std::sqrt(9305.419921875) * x * x * x * x * y * y * y * y * z * z * z + std::sqrt(11074.21875) * x * x * x * x * y * y * z * z * z * z * z + std::sqrt(787.5) * x * x * x * x * z * z * z * z * z * z * z + std::sqrt(15.5731201171875) * x * x * y * y * y * y * y * y * y * y * z - std::sqrt(3349.951171875) * x * x * y * y * y * y * y * y * z * z * z + std::sqrt(24916.9921875) * x * x * y * y * y * y * z * z * z * z * z - std::sqrt(7087.5) * x * x * y * y * z * z * z * z * z * z * z) + e_1 * (-std::sqrt(1.7303466796875) * x * x * x * x * x * x * x * x * z - std::sqrt(110.7421875) * x * x * x * x * x * x * y * y * z - std::sqrt(27.685546875) * x * x * x * x * x * x * z * z * z - std::sqrt(173.03466796875) * x * x * x * x * y * y * y * y * z + std::sqrt(33914.794921875) * x * x * x * x * y * y * z * z * z - std::sqrt(8970.1171875) * x * x * x * x * z * z * z * z * z + std::sqrt(17303.466796875) * x * x * y * y * y * y * z * z * z - std::sqrt(35880.46875) * x * x * y * y * z * z * z * z * z + std::sqrt(7087.5) * x * x * z * z * z * z * z * z * z + std::sqrt(15.5731201171875) * y * y * y * y * y * y * y * y * z - std::sqrt(3349.951171875) * y * y * y * y * y * y * z * z * z + std::sqrt(24916.9921875) * y * y * y * y * z * z * z * z * z - std::sqrt(7087.5) * y * y * z * z * z * z * z * z * z) + e_2 * (-std::sqrt(996.6796875) * x * x * x * x * x * x * z + std::sqrt(24916.9921875) * x * x * x * x * y * y * z - std::sqrt(177187.5) * x * x * x * x * z * z * z + std::sqrt(24916.9921875) * x * x * y * y * y * y * z + std::sqrt(255150.0) * x * x * z * z * z * z * z - std::sqrt(996.6796875) * y * y * y * y * y * y * z + std::sqrt(177187.5) * y * y * y * y * z * z * z - std::sqrt(255150.0) * y * y * z * z * z * z * z) + e_3 * (-std::sqrt(276855.46875) * x * x * x * x * z + std::sqrt(398671.875) * x * x * y * y * z + std::sqrt(708750.0) * x * x * z * z * z + std::sqrt(99667.96875) * y * y * y * y * z - std::sqrt(708750.0) * y * y * z * z * z);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        pc_106[k] = e_0 * (std::sqrt(48.44970703125) * x * x * x * x * x * x * x * x * x * z * z - std::sqrt(193.798828125) * x * x * x * x * x * x * x * y * y * z * z - std::sqrt(1055.126953125) * x * x * x * x * x * x * x * z * z * z * z - std::sqrt(775.1953125) * x * x * x * x * x * y * y * y * y * z * z + std::sqrt(9496.142578125) * x * x * x * x * x * y * y * z * z * z * z + std::sqrt(1378.125) * x * x * x * x * x * z * z * z * z * z * z + std::sqrt(193.798828125) * x * x * x * y * y * y * y * y * y * z * z + std::sqrt(1055.126953125) * x * x * x * y * y * y * y * z * z * z * z - std::sqrt(22050.0) * x * x * x * y * y * z * z * z * z * z * z + std::sqrt(436.04736328125) * x * y * y * y * y * y * y * y * y * z * z - std::sqrt(9496.142578125) * x * y * y * y * y * y * y * z * z * z * z + std::sqrt(12403.125) * x * y * y * y * y * z * z * z * z * z * z) + e_1 * (std::sqrt(48.44970703125) * x * x * x * x * x * x * x * x * x - std::sqrt(193.798828125) * x * x * x * x * x * x * x * y * y + std::sqrt(1744.189453125) * x * x * x * x * x * x * x * z * z - std::sqrt(775.1953125) * x * x * x * x * x * y * y * y * y + std::sqrt(15697.705078125) * x * x * x * x * x * y * y * z * z - std::sqrt(775.1953125) * x * x * x * x * x * z * z * z * z + std::sqrt(193.798828125) * x * x * x * y * y * y * y * y * y + std::sqrt(15697.705078125) * x * x * x * y * y * y * y * z * z - std::sqrt(1119382.03125) * x * x * x * y * y * z * z * z * z + std::sqrt(49612.5) * x * x * x * z * z * z * z * z * z + std::sqrt(436.04736328125) * x * y * y * y * y * y * y * y * y + std::sqrt(1744.189453125) * x * y * y * y * y * y * y * z * z + std::sqrt(93798.6328125) * x * y * y * y * y * z * z * z * z + std::sqrt(49612.5) * x * y * y * z * z * z * z * z * z) + e_2 * (std::sqrt(19379.8828125) * x * x * x * x * x * x * x - std::sqrt(27907.03125) * x * x * x * x * x * y * y + std::sqrt(174418.9453125) * x * x * x * x * x * z * z + std::sqrt(775.1953125) * x * x * x * y * y * y * y - std::sqrt(4716288.28125) * x * x * x * y * y * z * z + std::sqrt(1004653.125) * x * x * x * z * z * z * z + std::sqrt(111628.125) * x * y * y * y * y * y * y + std::sqrt(2016283.0078125) * x * y * y * y * y * z * z + std::sqrt(1004653.125) * x * y * y * z * z * z * z + std::sqrt(198450.0) * x * z * z * z * z * z * z) + e_3 * (std::sqrt(1119382.03125) * x * x * x * x * x - std::sqrt(1500778.125) * x * x * x * y * y + std::sqrt(7144200.0) * x * x * x * z * z + std::sqrt(4716288.28125) * x * y * y * y * y + std::sqrt(7144200.0) * x * y * y * z * z + std::sqrt(12700800.0) * x * z * z * z * z) + e_4 * (std::sqrt(9041878.125) * x * x * x + std::sqrt(9041878.125) * x * y * y + std::sqrt(64297800.0) * x * z * z) + e_5 * (std::sqrt(16074450.0) * x);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        pc_107[k] = e_0 * (std::sqrt(2.01873779296875) * x * x * x * x * x * x * x * x * x * x * z - std::sqrt(32.2998046875) * x * x * x * x * x * x * x * x * y * y * z - std::sqrt(229.6875) * x * x * x * x * x * x * x * x * z * z * z - std::sqrt(8.074951171875) * x * x * x * x * x * x * y * y * y * y * z + std::sqrt(5742.1875) * x * x * x * x * x * x * y * y * z * z * z + std::sqrt(918.75) * x * x * x * x * x * x * z * z * z * z * z + std::sqrt(290.6982421875) * x * x * x * x * y * y * y * y * y * y * z - std::sqrt(2067.1875) * x * x * x * x * y * y * y * y * z * z * z - std::sqrt(33075.0) * x * x * x * x * y * y * z * z * z * z * z + std::sqrt(163.51776123046875) * x * x * y * y * y * y * y * y * y * y * z - std::sqrt(18604.6875) * x * x * y * y * y * y * y * y * z * z * z + std::sqrt(74418.75) * x * x * y * y * y * y * z * z * z * z * z) + e_1 * (std::sqrt(163.51776123046875) * x * x * x * x * x * x * x * x * z + std::sqrt(2616.2841796875) * x * x * x * x * x * x * y * y * z - std::sqrt(2067.1875) * x * x * x * x * x * x * z * z * z + std::sqrt(5886.639404296875) * x * x * x * x * y * y * y * y * z - std::sqrt(911629.6875) * x * x * x * x * y * y * z * z * z + std::sqrt(74418.75) * x * x * x * x * z * z * z * z * z + std::sqrt(2616.2841796875) * x * x * y * y * y * y * y * y * z + std::sqrt(167442.1875) * x * x * y * y * y * y * z * z * z + std::sqrt(297675.0) * x * x * y * y * z * z * z * z * z + std::sqrt(163.51776123046875) * y * y * y * y * y * y * y * y * z - std::sqrt(18604.6875) * y * y * y * y * y * y * z * z * z + std::sqrt(74418.75) * y * y * y * y * z * z * z * z * z) + e_2 * (std::sqrt(18604.6875) * x * x * x * x * x * x * z - std::sqrt(669768.75) * x * x * x * x * y * y * z + std::sqrt(297675.0) * x * x * x * x * z * z * z + std::sqrt(1506979.6875) * x * x * y * y * y * y * z + std::sqrt(1190700.0) * x * x * y * y * z * z * z + std::sqrt(1190700.0) * x * x * z * z * z * z * z + std::sqrt(297675.0) * y * y * y * y * z * z * z + std::sqrt(1190700.0) * y * y * z * z * z * z * z) + e_3 * (std::sqrt(911629.6875) * x * x * x * x * z + std::sqrt(3646518.75) * x * x * y * y * z + std::sqrt(25930800.0) * x * x * z * z * z + std::sqrt(911629.6875) * y * y * y * y * z + std::sqrt(25930800.0) * y * y * z * z * z + std::sqrt(529200.0) * z * z * z * z * z) + e_4 * (std::sqrt(32818668.75) * x * x * z + std::sqrt(32818668.75) * y * y * z + std::sqrt(19051200.0) * z * z * z) + e_5 * (std::sqrt(24111675.0) * z);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        pc_108[k] = e_0 * (-std::sqrt(36.3372802734375) * x * x * x * x * x * x * x * x * x * z * z + std::sqrt(2325.5859375) * x * x * x * x * x * x * x * y * y * z * z + std::sqrt(258.3984375) * x * x * x * x * x * x * x * z * z * z * z - std::sqrt(3633.72802734375) * x * x * x * x * x * y * y * y * y * z * z - std::sqrt(20930.2734375) * x * x * x * x * x * y * y * z * z * z * z - std::sqrt(9302.34375) * x * x * x * y * y * y * y * y * y * z * z + std::sqrt(93281.8359375) * x * x * x * y * y * y * y * z * z * z * z + std::sqrt(327.0355224609375) * x * y * y * y * y * y * y * y * y * z * z - std::sqrt(2325.5859375) * x * y * y * y * y * y * y * z * z * z * z) + e_1 * (-std::sqrt(36.3372802734375) * x * x * x * x * x * x * x * x * x + std::sqrt(2325.5859375) * x * x * x * x * x * x * x * y * y - std::sqrt(5232.568359375) * x * x * x * x * x * x * x * z * z - std::sqrt(3633.72802734375) * x * x * x * x * x * y * y * y * y - std::sqrt(47093.115234375) * x * x * x * x * x * y * y * z * z + std::sqrt(37209.375) * x * x * x * x * x * z * z * z * z - std::sqrt(9302.34375) * x * x * x * y * y * y * y * y * y - std::sqrt(47093.115234375) * x * x * x * y * y * y * y * z * z + std::sqrt(148837.5) * x * x * x * y * y * z * z * z * z + std::sqrt(327.0355224609375) * x * y * y * y * y * y * y * y * y - std::sqrt(5232.568359375) * x * y * y * y * y * y * y * z * z + std::sqrt(37209.375) * x * y * y * y * y * z * z * z * z) + e_2 * (-std::sqrt(14534.912109375) * x * x * x * x * x * x * x + std::sqrt(47093.115234375) * x * x * x * x * x * y * y - std::sqrt(83721.09375) * x * x * x * x * x * z * z - std::sqrt(1284304.833984375) * x * x * x * y * y * y * y - std::sqrt(334884.375) * x * x * x * y * y * z * z + std::sqrt(1339537.5) * x * x * x * z * z * z * z + std::sqrt(5232.568359375) * x * y * y * y * y * y * y - std::sqrt(83721.09375) * x * y * y * y * y * z * z + std::sqrt(1339537.5) * x * y * y * z * z * z * z) + e_3 * (-std::sqrt(753489.84375) * x * x * x * x * x - std::sqrt(3013959.375) * x * x * x * y * y + std::sqrt(1339537.5) * x * x * x * z * z - std::sqrt(753489.84375) * x * y * y * y * y + std::sqrt(1339537.5) * x * y * y * z * z + std::sqrt(2381400.0) * x * z * z * z * z) + e_4 * (-std::sqrt(5358150.0) * x * x * x - std::sqrt(5358150.0) * x * y * y + std::sqrt(12055837.5) * x * z * z) + e_5 * (-std::sqrt(1339537.5) * x);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        pc_109[k] = e_0 * (-std::sqrt(3.63372802734375) * x * x * x * x * x * x * x * x * x * x * z + std::sqrt(523.2568359375) * x * x * x * x * x * x * x * x * y * y * z + std::sqrt(25.83984375) * x * x * x * x * x * x * x * x * z * z * z - std::sqrt(1758.724365234375) * x * x * x * x * x * x * y * y * y * y * z - std::sqrt(4366.93359375) * x * x * x * x * x * x * y * y * z * z * z - std::sqrt(1453.4912109375) * x * x * x * x * y * y * y * y * y * y * z + std::sqrt(31653.80859375) * x * x * x * x * y * y * y * y * z * z * z + std::sqrt(817.5888061523438) * x * x * y * y * y * y * y * y * y * y * z - std::sqrt(5813.96484375) * x * x * y * y * y * y * y * y * z * z * z) + e_1 * (-std::sqrt(2271.0800170898438) * x * x * x * x * x * x * x * x * z + std::sqrt(36337.2802734375) * x * x * x * x * x * x * y * y * z + std::sqrt(5813.96484375) * x * x * x * x * x * x * z * z * z - std::sqrt(445131.6833496094) * x * x * x * x * y * y * y * y * z + std::sqrt(5813.96484375) * x * x * x * x * y * y * z * z * z + std::sqrt(117732.7880859375) * x * x * y * y * y * y * y * y * z - std::sqrt(5813.96484375) * x * x * y * y * y * y * z * z * z + std::sqrt(817.5888061523438) * y * y * y * y * y * y * y * y * z - std::sqrt(5813.96484375) * y * y * y * y * y * y * z * z * z) + e_2 * (-std::sqrt(209302.734375) * x * x * x * x * x * x * z - std::sqrt(209302.734375) * x * x * x * x * y * y * z + std::sqrt(372093.75) * x * x * x * x * z * z * z + std::sqrt(209302.734375) * x * x * y * y * y * y * z + std::sqrt(209302.734375) * y * y * y * y * y * y * z - std::sqrt(372093.75) * y * y * y * y * z * z * z) + e_3 * (-std::sqrt(3348843.75) * x * x * x * x * z + std::sqrt(1488375.0) * x * x * z * z * z + std::sqrt(3348843.75) * y * y * y * y * z - std::sqrt(1488375.0) * y * y * z * z * z) + e_4 * (-std::sqrt(3348843.75) * x * x * z + std::sqrt(3348843.75) * y * y * z);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        pc_110[k] = e_0 * (-std::sqrt(3.028106689453125) * x * x * x * x * x * x * x * x * x * x * y + std::sqrt(148.37722778320312) * x * x * x * x * x * x * x * x * y * y * y + std::sqrt(302.8106689453125) * x * x * x * x * x * x * x * x * y * z * z - std::sqrt(81.8800048828125) * x * x * x * x * x * x * y * y * y * y * y - std::sqrt(19379.8828125) * x * x * x * x * x * x * y * y * y * z * z - std::sqrt(302.8106689453125) * x * x * x * x * y * y * y * y * y * y * y + std::sqrt(52761.73095703125) * x * x * x * x * y * y * y * y * y * z * z + std::sqrt(27.252960205078125) * x * x * y * y * y * y * y * y * y * y * y - std::sqrt(3100.78125) * x * x * y * y * y * y * y * y * y * z * z - std::sqrt(0.121124267578125) * y * y * y * y * y * y * y * y * y * y * y + std::sqrt(12.1124267578125) * y * y * y * y * y * y * y * y * y * z * z) + e_1 * (-std::sqrt(593.5089111328125) * x * x * x * x * x * x * x * x * y + std::sqrt(12403.125) * x * x * x * x * x * x * y * y * y + std::sqrt(4844.970703125) * x * x * x * x * x * x * y * z * z - std::sqrt(73692.00439453125) * x * x * x * x * y * y * y * y * y + std::sqrt(43604.736328125) * x * x * x * x * y * y * y * z * z + std::sqrt(775.1953125) * x * x * y * y * y * y * y * y * y + std::sqrt(43604.736328125) * x * x * y * y * y * y * y * z * z - std::sqrt(109.0118408203125) * y * y * y * y * y * y * y * y * y + std::sqrt(4844.970703125) * y * y * y * y * y * y * y * z * z) + e_2 * (-std::sqrt(19379.8828125) * x * x * x * x * x * x * y - std::sqrt(174418.9453125) * x * x * x * x * y * y * y + std::sqrt(697675.78125) * x * x * x * x * y * z * z - std::sqrt(174418.9453125) * x * x * y * y * y * y * y + std::sqrt(2790703.125) * x * x * y * y * y * z * z - std::sqrt(19379.8828125) * y * y * y * y * y * y * y + std::sqrt(697675.78125) * y * y * y * y * y * z * z) + e_3 * (-std::sqrt(697675.78125) * x * x * x * x * y - std::sqrt(2790703.125) * x * x * y * y * y + std::sqrt(11162812.5) * x * x * y * z * z - std::sqrt(697675.78125) * y * y * y * y * y + std::sqrt(11162812.5) * y * y * y * z * z) + e_4 * (-std::sqrt(2790703.125) * x * x * y - std::sqrt(2790703.125) * y * y * y + std::sqrt(11162812.5) * y * z * z) + e_5 * (-std::sqrt(446512.5) * y);

        pc_111[k] = e_0 * (-std::sqrt(19.3798828125) * x * x * x * x * x * x * x * x * x * y * z + std::sqrt(697.67578125) * x * x * x * x * x * x * x * y * y * y * z + std::sqrt(1937.98828125) * x * x * x * x * x * x * x * y * z * z * z - std::sqrt(94961.42578125) * x * x * x * x * x * y * y * y * z * z * z - std::sqrt(697.67578125) * x * x * x * y * y * y * y * y * y * y * z + std::sqrt(94961.42578125) * x * x * x * y * y * y * y * y * z * z * z + std::sqrt(19.3798828125) * x * y * y * y * y * y * y * y * y * y * z - std::sqrt(1937.98828125) * x * y * y * y * y * y * y * y * z * z * z) + e_1 * (std::sqrt(2790.703125) * x * x * x * x * x * x * x * y * z - std::sqrt(136744.453125) * x * x * x * x * x * y * y * y * z + std::sqrt(136744.453125) * x * x * x * y * y * y * y * y * z - std::sqrt(2790.703125) * x * y * y * y * y * y * y * y * z);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        pc_112[k] = e_0 * (std::sqrt(0.605621337890625) * x * x * x * x * x * x * x * x * x * x * y - std::sqrt(11.372222900390625) * x * x * x * x * x * x * x * x * y * y * y - std::sqrt(196.2213134765625) * x * x * x * x * x * x * x * x * y * z * z - std::sqrt(45.4888916015625) * x * x * x * x * x * x * y * y * y * y * y + std::sqrt(5581.40625) * x * x * x * x * x * x * y * y * y * z * z + std::sqrt(3875.9765625) * x * x * x * x * x * x * y * z * z * z * z - std::sqrt(0.2691650390625) * x * x * x * x * y * y * y * y * y * y * y + std::sqrt(2180.23681640625) * x * x * x * x * y * y * y * y * y * z * z - std::sqrt(155469.7265625) * x * x * x * x * y * y * y * z * z * z * z + std::sqrt(3.297271728515625) * x * x * y * y * y * y * y * y * y * y * y - std::sqrt(1395.3515625) * x * x * y * y * y * y * y * y * y * z * z + std::sqrt(34883.7890625) * x * x * y * y * y * y * y * z * z * z * z - std::sqrt(0.067291259765625) * y * y * y * y * y * y * y * y * y * y * y + std::sqrt(21.8023681640625) * y * y * y * y * y * y * y * y * y * z * z - std::sqrt(430.6640625) * y * y * y * y * y * y * y * z * z * z * z) + e_1 * (std::sqrt(118.7017822265625) * x * x * x * x * x * x * x * x * y - std::sqrt(7596.9140625) * x * x * x * x * x * x * y * y * y + std::sqrt(42209.384765625) * x * x * x * x * x * x * y * z * z - std::sqrt(4273.26416015625) * x * x * x * x * y * y * y * y * y - std::sqrt(605621.337890625) * x * x * x * x * y * y * y * z * z - std::sqrt(62015.625) * x * x * x * x * y * z * z * z * z + std::sqrt(620.15625) * x * x * y * y * y * y * y * y * y + std::sqrt(379884.462890625) * x * x * y * y * y * y * y * z * z - std::sqrt(248062.5) * x * x * y * y * y * z * z * z * z - std::sqrt(60.5621337890625) * y * y * y * y * y * y * y * y * y + std::sqrt(38.759765625) * y * y * y * y * y * y * y * z * z - std::sqrt(62015.625) * y * y * y * y * y * z * z * z * z) + e_2 * (std::sqrt(15503.90625) * x * x * x * x * x * x * y - std::sqrt(1875972.65625) * x * x * x * x * y * y * y - std::sqrt(139535.15625) * x * x * x * x * y * z * z + std::sqrt(139535.15625) * x * x * y * y * y * y * y - std::sqrt(558140.625) * x * x * y * y * y * z * z - std::sqrt(2232562.5) * x * x * y * z * z * z * z - std::sqrt(15503.90625) * y * y * y * y * y * y * y - std::sqrt(139535.15625) * y * y * y * y * y * z * z - std::sqrt(2232562.5) * y * y * y * z * z * z * z) + e_3 * (-std::sqrt(759691.40625) * x * x * x * x * y - std::sqrt(3038765.625) * x * x * y * y * y - std::sqrt(20093062.5) * x * x * y * z * z - std::sqrt(759691.40625) * y * y * y * y * y - std::sqrt(20093062.5) * y * y * y * z * z - std::sqrt(3969000.0) * y * z * z * z * z) + e_4 * (-std::sqrt(13953515.625) * x * x * y - std::sqrt(13953515.625) * y * y * y - std::sqrt(55814062.5) * y * z * z) + e_5 * (-std::sqrt(20093062.5) * y);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        pc_113[k] = e_0 * (std::sqrt(6.4599609375) * x * x * x * x * x * x * x * x * x * y * z - std::sqrt(103.359375) * x * x * x * x * x * x * x * y * y * y * z - std::sqrt(930.234375) * x * x * x * x * x * x * x * y * z * z * z - std::sqrt(645.99609375) * x * x * x * x * x * y * y * y * y * y * z + std::sqrt(23255.859375) * x * x * x * x * x * y * y * y * z * z * z + std::sqrt(2583.984375) * x * x * x * x * x * y * z * z * z * z * z - std::sqrt(103.359375) * x * x * x * y * y * y * y * y * y * y * z + std::sqrt(23255.859375) * x * x * x * y * y * y * y * y * z * z * z - std::sqrt(93023.4375) * x * x * x * y * y * y * z * z * z * z * z + std::sqrt(6.4599609375) * x * y * y * y * y * y * y * y * y * y * z - std::sqrt(930.234375) * x * y * y * y * y * y * y * y * z * z * z + std::sqrt(2583.984375) * x * y * y * y * y * y * z * z * z * z * z) + e_1 * (-std::sqrt(930.234375) * x * x * x * x * x * x * x * y * z - std::sqrt(103.359375) * x * x * x * x * x * y * y * y * z + std::sqrt(105840.0) * x * x * x * x * x * y * z * z * z - std::sqrt(103.359375) * x * x * x * y * y * y * y * y * z - std::sqrt(165375.0) * x * x * x * y * z * z * z * z * z - std::sqrt(930.234375) * x * y * y * y * y * y * y * y * z + std::sqrt(105840.0) * x * y * y * y * y * y * z * z * z - std::sqrt(165375.0) * x * y * y * y * z * z * z * z * z) + e_2 * (std::sqrt(23255.859375) * x * x * x * x * x * y * z - std::sqrt(10335.9375) * x * x * x * y * y * y * z - std::sqrt(165375.0) * x * x * x * y * z * z * z + std::sqrt(23255.859375) * x * y * y * y * y * y * z - std::sqrt(165375.0) * x * y * y * y * z * z * z - std::sqrt(1488375.0) * x * y * z * z * z * z * z) + e_3 * (-std::sqrt(23814000.0) * x * y * z * z * z) + e_4 * (-std::sqrt(13395375.0) * x * y * z);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        pc_114[k] = e_0 * (-std::sqrt(0.05767822265625) * x * x * x * x * x * x * x * x * x * x * y + std::sqrt(0.51910400390625) * x * x * x * x * x * x * x * x * y * y * y + std::sqrt(27.916259765625) * x * x * x * x * x * x * x * x * y * z * z + std::sqrt(11.304931640625) * x * x * x * x * x * x * y * y * y * y * y - std::sqrt(446.66015625) * x * x * x * x * x * x * y * y * y * z * z - std::sqrt(945.0) * x * x * x * x * x * x * y * z * z * z * z + std::sqrt(11.304931640625) * x * x * x * x * y * y * y * y * y * y * y - std::sqrt(2791.6259765625) * x * x * x * x * y * y * y * y * y * z * z + std::sqrt(23625.0) * x * x * x * x * y * y * y * z * z * z * z + std::sqrt(369.140625) * x * x * x * x * y * z * z * z * z * z * z + std::sqrt(0.51910400390625) * x * x * y * y * y * y * y * y * y * y * y - std::sqrt(446.66015625) * x * x * y * y * y * y * y * y * y * z * z + std::sqrt(23625.0) * x * x * y * y * y * y * y * z * z * z * z - std::sqrt(13289.0625) * x * x * y * y * y * z * z * z * z * z * z - std::sqrt(0.05767822265625) * y * y * y * y * y * y * y * y * y * y * y + std::sqrt(27.916259765625) * y * y * y * y * y * y * y * y * y * z * z - std::sqrt(945.0) * y * y * y * y * y * y * y * z * z * z * z + std::sqrt(369.140625) * y * y * y * y * y * z * z * z * z * z * z) + e_1 * (-std::sqrt(11.304931640625) * x * x * x * x * x * x * x * x * y + std::sqrt(1066.81640625) * x * x * x * x * x * x * y * y * y - std::sqrt(9981.5625) * x * x * x * x * x * x * y * z * z + std::sqrt(4652.0947265625) * x * x * x * x * y * y * y * y * y + std::sqrt(5906.25) * x * x * x * x * y * y * y * z * z + std::sqrt(83056.640625) * x * x * x * x * y * z * z * z * z + std::sqrt(623.84765625) * x * x * y * y * y * y * y * y * y + std::sqrt(26046.5625) * x * x * y * y * y * y * y * z * z + std::sqrt(533039.0625) * x * x * y * y * y * z * z * z * z - std::sqrt(53156.25) * x * x * y * z * z * z * z * z * z - std::sqrt(51.910400390625) * y * y * y * y * y * y * y * y * y - std::sqrt(236.25) * y * y * y * y * y * y * y * z * z - std::sqrt(41476.640625) * y * y * y * y * y * z * z * z * z + std::sqrt(5906.25) * y * y * y * z * z * z * z * z * z) + e_2 * (-std::sqrt(2307.12890625) * x * x * x * x * x * x * y + std::sqrt(389904.78515625) * x * x * x * x * y * y * y + std::sqrt(53156.25) * x * x * x * x * y * z * z + std::sqrt(299834.47265625) * x * x * y * y * y * y * y + std::sqrt(10418625.0) * x * x * y * y * y * z * z + std::sqrt(53156.25) * x * x * y * z * z * z * z - std::sqrt(15596.19140625) * y * y * y * y * y * y * y - std::sqrt(478406.25) * y * y * y * y * y * z * z - std::sqrt(5906.25) * y * y * y * z * z * z * z) + e_3 * (std::sqrt(212625.0) * x * x * x * x * y + std::sqrt(17222625.0) * x * x * y * y * y + std::sqrt(17222625.0) * x * x * y * z * z - std::sqrt(850500.0) * y * y * y * y * y - std::sqrt(1913625.0) * y * y * y * z * z) + e_4 * (std::sqrt(23441906.25) * x * x * y - std::sqrt(2604656.25) * y * y * y);
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

        pc_115[k] = e_0 * (-std::sqrt(0.86517333984375) * x * x * x * x * x * x * x * x * x * x * z + std::sqrt(7.78656005859375) * x * x * x * x * x * x * x * x * y * y * z + std::sqrt(138.812255859375) * x * x * x * x * x * x * x * x * z * z * z + std::sqrt(169.573974609375) * x * x * x * x * x * x * y * y * y * y * z - std::sqrt(2220.99609375) * x * x * x * x * x * x * y * y * z * z * z - std::sqrt(640.08984375) * x * x * x * x * x * x * z * z * z * z * z + std::sqrt(169.573974609375) * x * x * x * x * y * y * y * y * y * y * z - std::sqrt(13881.2255859375) * x * x * x * x * y * y * y * y * z * z * z + std::sqrt(16002.24609375) * x * x * x * x * y * y * z * z * z * z * z + std::sqrt(24.609375) * x * x * x * x * z * z * z * z * z * z * z + std::sqrt(7.78656005859375) * x * x * y * y * y * y * y * y * y * y * z - std::sqrt(2220.99609375) * x * x * y * y * y * y * y * y * z * z * z + std::sqrt(16002.24609375) * x * x * y * y * y * y * z * z * z * z * z - std::sqrt(885.9375) * x * x * y * y * z * z * z * z * z * z * z - std::sqrt(0.86517333984375) * y * y * y * y * y * y * y * y * y * y * z + std::sqrt(138.812255859375) * y * y * y * y * y * y * y * y * z * z * z - std::sqrt(640.08984375) * y * y * y * y * y * y * z * z * z * z * z + std::sqrt(24.609375) * y * y * y * y * z * z * z * z * z * z * z) + e_1 * (-std::sqrt(13.8427734375) * x * x * x * x * x * x * x * x * z + std::sqrt(221.484375) * x * x * x * x * x * x * y * y * z + std::sqrt(885.9375) * x * x * x * x * x * x * z * z * z + std::sqrt(1384.27734375) * x * x * x * x * y * y * y * y * z - std::sqrt(22148.4375) * x * x * x * x * y * y * z * z * z - std::sqrt(22148.4375) * x * x * x * x * z * z * z * z * z + std::sqrt(221.484375) * x * x * y * y * y * y * y * y * z - std::sqrt(22148.4375) * x * x * y * y * y * y * z * z * z + std::sqrt(797343.75) * x * x * y * y * z * z * z * z * z - std::sqrt(13.8427734375) * y * y * y * y * y * y * y * y * z + std::sqrt(885.9375) * y * y * y * y * y * y * z * z * z - std::sqrt(22148.4375) * y * y * y * y * z * z * z * z * z) + e_2 * (-std::sqrt(354375.0) * x * x * x * x * z * z * z + std::sqrt(12757500.0) * x * x * y * y * z * z * z - std::sqrt(354375.0) * y * y * y * y * z * z * z) + e_3 * (-std::sqrt(354375.0) * x * x * x * x * z + std::sqrt(12757500.0) * x * x * y * y * z - std::sqrt(354375.0) * y * y * y * y * z);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        pc_116[k] = e_0 * (-std::sqrt(0.05767822265625) * x * x * x * x * x * x * x * x * x * x * x + std::sqrt(0.51910400390625) * x * x * x * x * x * x * x * x * x * y * y + std::sqrt(27.916259765625) * x * x * x * x * x * x * x * x * x * z * z + std::sqrt(11.304931640625) * x * x * x * x * x * x * x * y * y * y * y - std::sqrt(446.66015625) * x * x * x * x * x * x * x * y * y * z * z - std::sqrt(945.0) * x * x * x * x * x * x * x * z * z * z * z + std::sqrt(11.304931640625) * x * x * x * x * x * y * y * y * y * y * y - std::sqrt(2791.6259765625) * x * x * x * x * x * y * y * y * y * z * z + std::sqrt(23625.0) * x * x * x * x * x * y * y * z * z * z * z + std::sqrt(369.140625) * x * x * x * x * x * z * z * z * z * z * z + std::sqrt(0.51910400390625) * x * x * x * y * y * y * y * y * y * y * y - std::sqrt(446.66015625) * x * x * x * y * y * y * y * y * y * z * z + std::sqrt(23625.0) * x * x * x * y * y * y * y * z * z * z * z - std::sqrt(13289.0625) * x * x * x * y * y * z * z * z * z * z * z - std::sqrt(0.05767822265625) * x * y * y * y * y * y * y * y * y * y * y + std::sqrt(27.916259765625) * x * y * y * y * y * y * y * y * y * z * z - std::sqrt(945.0) * x * y * y * y * y * y * y * z * z * z * z + std::sqrt(369.140625) * x * y * y * y * y * z * z * z * z * z * z) + e_1 * (-std::sqrt(51.910400390625) * x * x * x * x * x * x * x * x * x + std::sqrt(623.84765625) * x * x * x * x * x * x * x * y * y - std::sqrt(236.25) * x * x * x * x * x * x * x * z * z + std::sqrt(4652.0947265625) * x * x * x * x * x * y * y * y * y + std::sqrt(26046.5625) * x * x * x * x * x * y * y * z * z - std::sqrt(41476.640625) * x * x * x * x * x * z * z * z * z + std::sqrt(1066.81640625) * x * x * x * y * y * y * y * y * y + std::sqrt(5906.25) * x * x * x * y * y * y * y * z * z + std::sqrt(533039.0625) * x * x * x * y * y * z * z * z * z + std::sqrt(5906.25) * x * x * x * z * z * z * z * z * z - std::sqrt(11.304931640625) * x * y * y * y * y * y * y * y * y - std::sqrt(9981.5625) * x * y * y * y * y * y * y * z * z + std::sqrt(83056.640625) * x * y * y * y * y * z * z * z * z - std::sqrt(53156.25) * x * y * y * z * z * z * z * z * z) + e_2 * (-std::sqrt(15596.19140625) * x * x * x * x * x * x * x + std::sqrt(299834.47265625) * x * x * x * x * x * y * y - std::sqrt(478406.25) * x * x * x * x * x * z * z + std::sqrt(389904.78515625) * x * x * x * y * y * y * y + std::sqrt(10418625.0) * x * x * x * y * y * z * z - std::sqrt(5906.25) * x * x * x * z * z * z * z - std::sqrt(2307.12890625) * x * y * y * y * y * y * y + std::sqrt(53156.25) * x * y * y * y * y * z * z + std::sqrt(53156.25) * x * y * y * z * z * z * z) + e_3 * (-std::sqrt(850500.0) * x * x * x * x * x + std::sqrt(17222625.0) * x * x * x * y * y - std::sqrt(1913625.0) * x * x * x * z * z + std::sqrt(212625.0) * x * y * y * y * y + std::sqrt(17222625.0) * x * y * y * z * z) + e_4 * (-std::sqrt(2604656.25) * x * x * x + std::sqrt(23441906.25) * x * y * y);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        pc_117[k] = e_0 * (std::sqrt(1.614990234375) * x * x * x * x * x * x * x * x * x * x * z - std::sqrt(40.374755859375) * x * x * x * x * x * x * x * x * y * y * z - std::sqrt(232.55859375) * x * x * x * x * x * x * x * x * z * z * z - std::sqrt(58.1396484375) * x * x * x * x * x * x * y * y * y * y * z + std::sqrt(8372.109375) * x * x * x * x * x * x * y * y * z * z * z + std::sqrt(645.99609375) * x * x * x * x * x * x * z * z * z * z * z + std::sqrt(58.1396484375) * x * x * x * x * y * y * y * y * y * y * z - std::sqrt(31653.80859375) * x * x * x * x * y * y * z * z * z * z * z + std::sqrt(40.374755859375) * x * x * y * y * y * y * y * y * y * y * z - std::sqrt(8372.109375) * x * x * y * y * y * y * y * y * z * z * z + std::sqrt(31653.80859375) * x * x * y * y * y * y * z * z * z * z * z - std::sqrt(1.614990234375) * y * y * y * y * y * y * y * y * y * y * z + std::sqrt(232.55859375) * y * y * y * y * y * y * y * y * z * z * z - std::sqrt(645.99609375) * y * y * y * y * y * y * z * z * z * z * z) + e_1 * (std::sqrt(25.83984375) * x * x * x * x * x * x * x * x * z + std::sqrt(2583.984375) * x * x * x * x * x * x * y * y * z - std::sqrt(6615.0) * x * x * x * x * x * x * z * z * z - std::sqrt(165375.0) * x * x * x * x * y * y * z * z * z + std::sqrt(41343.75) * x * x * x * x * z * z * z * z * z - std::sqrt(2583.984375) * x * x * y * y * y * y * y * y * z + std::sqrt(165375.0) * x * x * y * y * y * y * z * z * z - std::sqrt(25.83984375) * y * y * y * y * y * y * y * y * z + std::sqrt(6615.0) * y * y * y * y * y * y * z * z * z - std::sqrt(41343.75) * y * y * y * y * z * z * z * z * z) + e_2 * (-std::sqrt(645.99609375) * x * x * x * x * x * x * z - std::sqrt(52325.68359375) * x * x * x * x * y * y * z + std::sqrt(41343.75) * x * x * x * x * z * z * z + std::sqrt(52325.68359375) * x * x * y * y * y * y * z + std::sqrt(372093.75) * x * x * z * z * z * z * z + std::sqrt(645.99609375) * y * y * y * y * y * y * z - std::sqrt(41343.75) * y * y * y * y * z * z * z - std::sqrt(372093.75) * y * y * z * z * z * z * z) + e_3 * (std::sqrt(5953500.0) * x * x * z * z * z - std::sqrt(5953500.0) * y * y * z * z * z) + e_4 * (std::sqrt(3348843.75) * x * x * z - std::sqrt(3348843.75) * y * y * z);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        pc_118[k] = e_0 * (std::sqrt(0.067291259765625) * x * x * x * x * x * x * x * x * x * x * x - std::sqrt(3.297271728515625) * x * x * x * x * x * x * x * x * x * y * y - std::sqrt(21.8023681640625) * x * x * x * x * x * x * x * x * x * z * z + std::sqrt(0.2691650390625) * x * x * x * x * x * x * x * y * y * y * y + std::sqrt(1395.3515625) * x * x * x * x * x * x * x * y * y * z * z + std::sqrt(430.6640625) * x * x * x * x * x * x * x * z * z * z * z + std::sqrt(45.4888916015625) * x * x * x * x * x * y * y * y * y * y * y - std::sqrt(2180.23681640625) * x * x * x * x * x * y * y * y * y * z * z - std::sqrt(34883.7890625) * x * x * x * x * x * y * y * z * z * z * z + std::sqrt(11.372222900390625) * x * x * x * y * y * y * y * y * y * y * y - std::sqrt(5581.40625) * x * x * x * y * y * y * y * y * y * z * z + std::sqrt(155469.7265625) * x * x * x * y * y * y * y * z * z * z * z - std::sqrt(0.605621337890625) * x * y * y * y * y * y * y * y * y * y * y + std::sqrt(196.2213134765625) * x * y * y * y * y * y * y * y * y * z * z - std::sqrt(3875.9765625) * x * y * y * y * y * y * y * z * z * z * z) + e_1 * (std::sqrt(60.5621337890625) * x * x * x * x * x * x * x * x * x - std::sqrt(620.15625) * x * x * x * x * x * x * x * y * y - std::sqrt(38.759765625) * x * x * x * x * x * x * x * z * z + std::sqrt(4273.26416015625) * x * x * x * x * x * y * y * y * y - std::sqrt(379884.462890625) * x * x * x * x * x * y * y * z * z + std::sqrt(62015.625) * x * x * x * x * x * z * z * z * z + std::sqrt(7596.9140625) * x * x * x * y * y * y * y * y * y + std::sqrt(605621.337890625) * x * x * x * y * y * y * y * z * z + std::sqrt(248062.5) * x * x * x * y * y * z * z * z * z - std::sqrt(118.7017822265625) * x * y * y * y * y * y * y * y * y - std::sqrt(42209.384765625) * x * y * y * y * y * y * y * z * z + std::sqrt(62015.625) * x * y * y * y * y * z * z * z * z) + e_2 * (std::sqrt(15503.90625) * x * x * x * x * x * x * x - std::sqrt(139535.15625) * x * x * x * x * x * y * y + std::sqrt(139535.15625) * x * x * x * x * x * z * z + std::sqrt(1875972.65625) * x * x * x * y * y * y * y + std::sqrt(558140.625) * x * x * x * y * y * z * z + std::sqrt(2232562.5) * x * x * x * z * z * z * z - std::sqrt(15503.90625) * x * y * y * y * y * y * y + std::sqrt(139535.15625) * x * y * y * y * y * z * z + std::sqrt(2232562.5) * x * y * y * z * z * z * z) + e_3 * (std::sqrt(759691.40625) * x * x * x * x * x + std::sqrt(3038765.625) * x * x * x * y * y + std::sqrt(20093062.5) * x * x * x * z * z + std::sqrt(759691.40625) * x * y * y * y * y + std::sqrt(20093062.5) * x * y * y * z * z + std::sqrt(3969000.0) * x * z * z * z * z) + e_4 * (std::sqrt(13953515.625) * x * x * x + std::sqrt(13953515.625) * x * y * y + std::sqrt(55814062.5) * x * z * z) + e_5 * (std::sqrt(20093062.5) * x);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        pc_119[k] = e_0 * (-std::sqrt(1.21124267578125) * x * x * x * x * x * x * x * x * x * x * z + std::sqrt(146.56036376953125) * x * x * x * x * x * x * x * x * y * y * z + std::sqrt(121.124267578125) * x * x * x * x * x * x * x * x * z * z * z - std::sqrt(818.800048828125) * x * x * x * x * x * x * y * y * y * y * z - std::sqrt(17441.89453125) * x * x * x * x * x * x * y * y * z * z * z - std::sqrt(818.800048828125) * x * x * x * x * y * y * y * y * y * y * z + std::sqrt(174903.4423828125) * x * x * x * x * y * y * y * y * z * z * z + std::sqrt(146.56036376953125) * x * x * y * y * y * y * y * y * y * y * z - std::sqrt(17441.89453125) * x * x * y * y * y * y * y * y * z * z * z - std::sqrt(1.21124267578125) * y * y * y * y * y * y * y * y * y * y * z + std::sqrt(121.124267578125) * y * y * y * y * y * y * y * y * z * z * z) + e_1 * (-std::sqrt(19.3798828125) * x * x * x * x * x * x * x * x * z - std::sqrt(52403.203125) * x * x * x * x * x * x * y * y * z + std::sqrt(31007.8125) * x * x * x * x * x * x * z * z * z + std::sqrt(156977.05078125) * x * x * x * x * y * y * y * y * z + std::sqrt(279070.3125) * x * x * x * x * y * y * z * z * z - std::sqrt(52403.203125) * x * x * y * y * y * y * y * y * z + std::sqrt(279070.3125) * x * x * y * y * y * y * z * z * z - std::sqrt(19.3798828125) * y * y * y * y * y * y * y * y * z + std::sqrt(31007.8125) * y * y * y * y * y * y * z * z * z) + e_2 * (std::sqrt(7751.953125) * x * x * x * x * x * x * z + std::sqrt(69767.578125) * x * x * x * x * y * y * z + std::sqrt(2511632.8125) * x * x * x * x * z * z * z + std::sqrt(69767.578125) * x * x * y * y * y * y * z + std::sqrt(10046531.25) * x * x * y * y * z * z * z + std::sqrt(7751.953125) * y * y * y * y * y * y * z + std::sqrt(2511632.8125) * y * y * y * y * z * z * z) + e_3 * (std::sqrt(4465125.0) * x * x * x * x * z + std::sqrt(17860500.0) * x * x * y * y * z + std::sqrt(17860500.0) * x * x * z * z * z + std::sqrt(4465125.0) * y * y * y * y * z + std::sqrt(17860500.0) * y * y * z * z * z) + e_4 * (std::sqrt(54697781.25) * x * x * z + std::sqrt(54697781.25) * y * y * z + std::sqrt(4465125.0) * z * z * z) + e_5 * (std::sqrt(17860500.0) * z);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        pc_120[k] = e_0 * (-std::sqrt(0.121124267578125) * x * x * x * x * x * x * x * x * x * x * x + std::sqrt(27.252960205078125) * x * x * x * x * x * x * x * x * x * y * y + std::sqrt(12.1124267578125) * x * x * x * x * x * x * x * x * x * z * z - std::sqrt(302.8106689453125) * x * x * x * x * x * x * x * y * y * y * y - std::sqrt(3100.78125) * x * x * x * x * x * x * x * y * y * z * z - std::sqrt(81.8800048828125) * x * x * x * x * x * y * y * y * y * y * y + std::sqrt(52761.73095703125) * x * x * x * x * x * y * y * y * y * z * z + std::sqrt(148.37722778320312) * x * x * x * y * y * y * y * y * y * y * y - std::sqrt(19379.8828125) * x * x * x * y * y * y * y * y * y * z * z - std::sqrt(3.028106689453125) * x * y * y * y * y * y * y * y * y * y * y + std::sqrt(302.8106689453125) * x * y * y * y * y * y * y * y * y * z * z) + e_1 * (-std::sqrt(109.0118408203125) * x * x * x * x * x * x * x * x * x + std::sqrt(775.1953125) * x * x * x * x * x * x * x * y * y + std::sqrt(4844.970703125) * x * x * x * x * x * x * x * z * z - std::sqrt(73692.00439453125) * x * x * x * x * x * y * y * y * y + std::sqrt(43604.736328125) * x * x * x * x * x * y * y * z * z + std::sqrt(12403.125) * x * x * x * y * y * y * y * y * y + std::sqrt(43604.736328125) * x * x * x * y * y * y * y * z * z - std::sqrt(593.5089111328125) * x * y * y * y * y * y * y * y * y + std::sqrt(4844.970703125) * x * y * y * y * y * y * y * z * z) + e_2 * (-std::sqrt(19379.8828125) * x * x * x * x * x * x * x - std::sqrt(174418.9453125) * x * x * x * x * x * y * y + std::sqrt(697675.78125) * x * x * x * x * x * z * z - std::sqrt(174418.9453125) * x * x * x * y * y * y * y + std::sqrt(2790703.125) * x * x * x * y * y * z * z - std::sqrt(19379.8828125) * x * y * y * y * y * y * y + std::sqrt(697675.78125) * x * y * y * y * y * z * z) + e_3 * (-std::sqrt(697675.78125) * x * x * x * x * x - std::sqrt(2790703.125) * x * x * x * y * y + std::sqrt(11162812.5) * x * x * x * z * z - std::sqrt(697675.78125) * x * y * y * y * y + std::sqrt(11162812.5) * x * y * y * z * z) + e_4 * (-std::sqrt(2790703.125) * x * x * x - std::sqrt(2790703.125) * x * y * y + std::sqrt(11162812.5) * x * z * z) + e_5 * (-std::sqrt(446512.5) * x);

        pc_121[k] = e_0 * (std::sqrt(66.61834716796875) * x * x * x * x * x * x * x * x * x * y * z - std::sqrt(9593.0419921875) * x * x * x * x * x * x * x * y * y * y * z + std::sqrt(42305.315185546875) * x * x * x * x * x * y * y * y * y * y * z - std::sqrt(9593.0419921875) * x * x * x * y * y * y * y * y * y * y * z + std::sqrt(66.61834716796875) * x * y * y * y * y * y * y * y * y * y * z);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        pc_122[k] = e_0 * (std::sqrt(426.357421875) * x * x * x * x * x * x * x * x * y * z * z - std::sqrt(51589.248046875) * x * x * x * x * x * x * y * y * y * z * z + std::sqrt(95930.419921875) * x * x * x * x * y * y * y * y * y * z * z - std::sqrt(10658.935546875) * x * x * y * y * y * y * y * y * y * z * z) + e_1 * (std::sqrt(426.357421875) * x * x * x * x * x * x * x * x * y - std::sqrt(51589.248046875) * x * x * x * x * x * x * y * y * y - std::sqrt(10658.935546875) * x * x * x * x * x * x * y * z * z + std::sqrt(95930.419921875) * x * x * x * x * y * y * y * y * y - std::sqrt(95930.419921875) * x * x * x * x * y * y * y * z * z - std::sqrt(10658.935546875) * x * x * y * y * y * y * y * y * y - std::sqrt(95930.419921875) * x * x * y * y * y * y * y * z * z - std::sqrt(10658.935546875) * y * y * y * y * y * y * y * z * z) + e_2 * (-std::sqrt(10658.935546875) * x * x * x * x * x * x * y - std::sqrt(95930.419921875) * x * x * x * x * y * y * y - std::sqrt(1534886.71875) * x * x * x * x * y * z * z - std::sqrt(95930.419921875) * x * x * y * y * y * y * y - std::sqrt(6139546.875) * x * x * y * y * y * z * z - std::sqrt(10658.935546875) * y * y * y * y * y * y * y - std::sqrt(1534886.71875) * y * y * y * y * y * z * z) + e_3 * (-std::sqrt(1534886.71875) * x * x * x * x * y - std::sqrt(6139546.875) * x * x * y * y * y - std::sqrt(24558187.5) * x * x * y * z * z - std::sqrt(1534886.71875) * y * y * y * y * y - std::sqrt(24558187.5) * y * y * y * z * z) + e_4 * (-std::sqrt(24558187.5) * x * x * y - std::sqrt(24558187.5) * y * y * y - std::sqrt(24558187.5) * y * z * z) + e_5 * (-std::sqrt(24558187.5) * y);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        pc_123[k] = e_0 * (-std::sqrt(13.32366943359375) * x * x * x * x * x * x * x * x * x * y * z + std::sqrt(1160.6396484375) * x * x * x * x * x * x * x * y * y * y * z + std::sqrt(852.71484375) * x * x * x * x * x * x * x * y * z * z * z + std::sqrt(53.294677734375) * x * x * x * x * x * y * y * y * y * y * z - std::sqrt(91050.99609375) * x * x * x * x * x * y * y * y * z * z * z - std::sqrt(592.1630859375) * x * x * x * y * y * y * y * y * y * y * z + std::sqrt(59216.30859375) * x * x * x * y * y * y * y * y * z * z * z + std::sqrt(37.01019287109375) * x * y * y * y * y * y * y * y * y * y * z - std::sqrt(2368.65234375) * x * y * y * y * y * y * y * y * z * z * z) + e_1 * (std::sqrt(3410.859375) * x * x * x * x * x * x * x * y * z - std::sqrt(13643.4375) * x * x * x * x * x * y * y * y * z - std::sqrt(85271.484375) * x * x * x * x * x * y * z * z * z + std::sqrt(85271.484375) * x * x * x * y * y * y * y * y * z - std::sqrt(341085.9375) * x * x * x * y * y * y * z * z * z - std::sqrt(85271.484375) * x * y * y * y * y * y * z * z * z) + e_2 * (-std::sqrt(5457375.0) * x * x * x * y * z * z * z - std::sqrt(5457375.0) * x * y * y * y * z * z * z) + e_3 * (-std::sqrt(5457375.0) * x * x * x * y * z - std::sqrt(5457375.0) * x * y * y * y * z - std::sqrt(21829500.0) * x * y * z * z * z) + e_4 * (-std::sqrt(49116375.0) * x * y * z);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        pc_124[k] = e_0 * (-std::sqrt(142.119140625) * x * x * x * x * x * x * x * x * y * z * z + std::sqrt(11511.650390625) * x * x * x * x * x * x * y * y * y * z * z + std::sqrt(568.4765625) * x * x * x * x * x * x * y * z * z * z * z + std::sqrt(3552.978515625) * x * x * x * x * y * y * y * y * y * z * z - std::sqrt(56847.65625) * x * x * x * x * y * y * y * z * z * z * z - std::sqrt(3552.978515625) * x * x * y * y * y * y * y * y * y * z * z + std::sqrt(14211.9140625) * x * x * y * y * y * y * y * z * z * z * z) + e_1 * (-std::sqrt(142.119140625) * x * x * x * x * x * x * x * x * y + std::sqrt(11511.650390625) * x * x * x * x * x * x * y * y * y + std::sqrt(17196.416015625) * x * x * x * x * x * x * y * z * z + std::sqrt(3552.978515625) * x * x * x * x * y * y * y * y * y + std::sqrt(600453.369140625) * x * x * x * x * y * y * y * z * z - std::sqrt(127907.2265625) * x * x * x * x * y * z * z * z * z - std::sqrt(3552.978515625) * x * x * y * y * y * y * y * y * y - std::sqrt(31976.806640625) * x * x * y * y * y * y * y * z * z - std::sqrt(56847.65625) * x * x * y * y * y * z * z * z * z - std::sqrt(3552.978515625) * y * y * y * y * y * y * y * z * z + std::sqrt(14211.9140625) * y * y * y * y * y * z * z * z * z) + e_2 * (std::sqrt(3552.978515625) * x * x * x * x * x * x * y + std::sqrt(2220611.572265625) * x * x * x * x * y * y * y + std::sqrt(1151165.0390625) * x * x * x * x * y * z * z - std::sqrt(287791.259765625) * x * x * y * y * y * y * y + std::sqrt(511628.90625) * x * x * y * y * y * z * z - std::sqrt(2046515.625) * x * x * y * z * z * z * z - std::sqrt(3552.978515625) * y * y * y * y * y * y * y - std::sqrt(127907.2265625) * y * y * y * y * y * z * z + std::sqrt(227390.625) * y * y * y * z * z * z * z) + e_3 * (std::sqrt(4604660.15625) * x * x * x * x * y + std::sqrt(2046515.625) * x * x * y * y * y - std::sqrt(511628.90625) * y * y * y * y * y) + e_4 * (std::sqrt(18418640.625) * x * x * y - std::sqrt(2046515.625) * y * y * y);
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

        pc_125[k] = e_0 * (std::sqrt(1.2689208984375) * x * x * x * x * x * x * x * x * x * y * z - std::sqrt(81.2109375) * x * x * x * x * x * x * x * y * y * y * z - std::sqrt(182.724609375) * x * x * x * x * x * x * x * y * z * z * z - std::sqrt(248.70849609375) * x * x * x * x * x * y * y * y * y * y * z + std::sqrt(14800.693359375) * x * x * x * x * x * y * y * y * z * z * z + std::sqrt(81.2109375) * x * x * x * x * x * y * z * z * z * z * z + std::sqrt(4568.115234375) * x * x * x * y * y * y * y * y * z * z * z - std::sqrt(8121.09375) * x * x * x * y * y * y * z * z * z * z * z + std::sqrt(31.7230224609375) * x * y * y * y * y * y * y * y * y * y * z - std::sqrt(4568.115234375) * x * y * y * y * y * y * y * y * z * z * z + std::sqrt(2030.2734375) * x * y * y * y * y * y * z * z * z * z * z) + e_1 * (-std::sqrt(730.8984375) * x * x * x * x * x * x * x * y * z + std::sqrt(324.84375) * x * x * x * x * x * y * y * y * z + std::sqrt(29317.1484375) * x * x * x * x * x * y * z * z * z + std::sqrt(2030.2734375) * x * x * x * y * y * y * y * y * z + std::sqrt(982652.34375) * x * x * x * y * y * y * z * z * z - std::sqrt(32484.375) * x * x * x * y * z * z * z * z * z - std::sqrt(586749.0234375) * x * y * y * y * y * y * z * z * z + std::sqrt(32484.375) * x * y * y * y * z * z * z * z * z) + e_2 * (std::sqrt(3248437.5) * x * x * x * y * y * y * z + std::sqrt(2079000.0) * x * x * x * y * z * z * z - std::sqrt(1169437.5) * x * y * y * y * y * y * z - std::sqrt(2079000.0) * x * y * y * y * z * z * z) + e_3 * (std::sqrt(10524937.5) * x * x * x * y * z - std::sqrt(10524937.5) * x * y * y * y * z);
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

        pc_126[k] = e_0 * (std::sqrt(19.0338134765625) * x * x * x * x * x * x * x * x * x * z * z - std::sqrt(1218.1640625) * x * x * x * x * x * x * x * y * y * z * z - std::sqrt(135.3515625) * x * x * x * x * x * x * x * z * z * z * z - std::sqrt(3730.62744140625) * x * x * x * x * x * y * y * y * y * z * z + std::sqrt(10963.4765625) * x * x * x * x * x * y * y * z * z * z * z + std::sqrt(5.4140625) * x * x * x * x * x * z * z * z * z * z * z + std::sqrt(3383.7890625) * x * x * x * y * y * y * y * z * z * z * z - std::sqrt(541.40625) * x * x * x * y * y * z * z * z * z * z * z + std::sqrt(475.8453369140625) * x * y * y * y * y * y * y * y * y * z * z - std::sqrt(3383.7890625) * x * y * y * y * y * y * y * z * z * z * z + std::sqrt(135.3515625) * x * y * y * y * y * z * z * z * z * z * z) + e_1 * (std::sqrt(19.0338134765625) * x * x * x * x * x * x * x * x * x - std::sqrt(1218.1640625) * x * x * x * x * x * x * x * y * y + std::sqrt(2740.869140625) * x * x * x * x * x * x * x * z * z - std::sqrt(3730.62744140625) * x * x * x * x * x * y * y * y * y - std::sqrt(222010.400390625) * x * x * x * x * x * y * y * z * z - std::sqrt(10963.4765625) * x * x * x * x * x * z * z * z * z - std::sqrt(68521.728515625) * x * x * x * y * y * y * y * z * z + std::sqrt(1096347.65625) * x * x * x * y * y * z * z * z * z + std::sqrt(475.8453369140625) * x * y * y * y * y * y * y * y * y + std::sqrt(68521.728515625) * x * y * y * y * y * y * y * z * z - std::sqrt(274086.9140625) * x * y * y * y * y * z * z * z * z) + e_2 * (std::sqrt(7613.525390625) * x * x * x * x * x * x * x - std::sqrt(616695.556640625) * x * x * x * x * x * y * y - std::sqrt(190338.134765625) * x * x * x * y * y * y * y + std::sqrt(190338.134765625) * x * y * y * y * y * y * y) + e_3 * (std::sqrt(121816.40625) * x * x * x * x * x - std::sqrt(12181640.625) * x * x * x * y * y + std::sqrt(3045410.15625) * x * y * y * y * y);
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

        pc_127[k] = e_0 * (std::sqrt(1.2689208984375) * x * x * x * x * x * x * x * x * x * x * z - std::sqrt(81.2109375) * x * x * x * x * x * x * x * x * y * y * z - std::sqrt(182.724609375) * x * x * x * x * x * x * x * x * z * z * z - std::sqrt(248.70849609375) * x * x * x * x * x * x * y * y * y * y * z + std::sqrt(14800.693359375) * x * x * x * x * x * x * y * y * z * z * z + std::sqrt(81.2109375) * x * x * x * x * x * x * z * z * z * z * z + std::sqrt(4568.115234375) * x * x * x * x * y * y * y * y * z * z * z - std::sqrt(8121.09375) * x * x * x * x * y * y * z * z * z * z * z + std::sqrt(31.7230224609375) * x * x * y * y * y * y * y * y * y * y * z - std::sqrt(4568.115234375) * x * x * y * y * y * y * y * y * z * z * z + std::sqrt(2030.2734375) * x * x * y * y * y * y * z * z * z * z * z) + e_1 * (std::sqrt(1.2689208984375) * x * x * x * x * x * x * x * x * z + std::sqrt(324.84375) * x * x * x * x * x * x * y * y * z - std::sqrt(27794.443359375) * x * x * x * x * x * x * z * z * z - std::sqrt(1142.02880859375) * x * x * x * x * y * y * y * y * z + std::sqrt(1766845.458984375) * x * x * x * x * y * y * z * z * z + std::sqrt(2030.2734375) * x * x * x * x * z * z * z * z * z - std::sqrt(2030.2734375) * x * x * y * y * y * y * y * y * z - std::sqrt(24870.849609375) * x * x * y * y * y * y * z * z * z - std::sqrt(73089.84375) * x * x * y * y * z * z * z * z * z + std::sqrt(31.7230224609375) * y * y * y * y * y * y * y * y * z - std::sqrt(4568.115234375) * y * y * y * y * y * y * z * z * z + std::sqrt(2030.2734375) * y * y * y * y * z * z * z * z * z) + e_2 * (-std::sqrt(50756.8359375) * x * x * x * x * x * x * z + std::sqrt(4111303.7109375) * x * x * x * x * y * y * z - std::sqrt(129937.5) * x * x * x * x * z * z * z - std::sqrt(456811.5234375) * x * x * y * y * y * y * z + std::sqrt(4677750.0) * x * x * y * y * z * z * z - std::sqrt(2030.2734375) * y * y * y * y * y * y * z - std::sqrt(129937.5) * y * y * y * y * z * z * z) + e_3 * (-std::sqrt(657808.59375) * x * x * x * x * z + std::sqrt(23681109.375) * x * x * y * y * z - std::sqrt(657808.59375) * y * y * y * y * z);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        pc_128[k] = e_0 * (-std::sqrt(35.52978515625) * x * x * x * x * x * x * x * x * x * z * z + std::sqrt(3552.978515625) * x * x * x * x * x * x * x * y * y * z * z + std::sqrt(142.119140625) * x * x * x * x * x * x * x * z * z * z * z - std::sqrt(568.4765625) * x * x * x * x * x * y * y * y * y * z * z - std::sqrt(17196.416015625) * x * x * x * x * x * y * y * z * z * z * z - std::sqrt(3552.978515625) * x * x * x * y * y * y * y * y * y * z * z + std::sqrt(31976.806640625) * x * x * x * y * y * y * y * z * z * z * z + std::sqrt(888.24462890625) * x * y * y * y * y * y * y * y * y * z * z - std::sqrt(3552.978515625) * x * y * y * y * y * y * y * z * z * z * z) + e_1 * (-std::sqrt(35.52978515625) * x * x * x * x * x * x * x * x * x + std::sqrt(3552.978515625) * x * x * x * x * x * x * x * y * y - std::sqrt(6963.837890625) * x * x * x * x * x * x * x * z * z - std::sqrt(568.4765625) * x * x * x * x * x * y * y * y * y + std::sqrt(103604.853515625) * x * x * x * x * x * y * y * z * z + std::sqrt(14211.9140625) * x * x * x * x * x * z * z * z * z - std::sqrt(3552.978515625) * x * x * x * y * y * y * y * y * y - std::sqrt(3552.978515625) * x * x * x * y * y * y * y * z * z - std::sqrt(56847.65625) * x * x * x * y * y * z * z * z * z + std::sqrt(888.24462890625) * x * y * y * y * y * y * y * y * y + std::sqrt(88824.462890625) * x * y * y * y * y * y * y * z * z - std::sqrt(127907.2265625) * x * y * y * y * y * z * z * z * z) + e_2 * (-std::sqrt(14211.9140625) * x * x * x * x * x * x * x + std::sqrt(511628.90625) * x * x * x * x * x * y * y - std::sqrt(127907.2265625) * x * x * x * x * x * z * z - std::sqrt(355297.8515625) * x * x * x * y * y * y * y + std::sqrt(511628.90625) * x * x * x * y * y * z * z + std::sqrt(227390.625) * x * x * x * z * z * z * z + std::sqrt(227390.625) * x * y * y * y * y * y * y + std::sqrt(1151165.0390625) * x * y * y * y * y * z * z - std::sqrt(2046515.625) * x * y * y * z * z * z * z) + e_3 * (-std::sqrt(511628.90625) * x * x * x * x * x + std::sqrt(2046515.625) * x * x * x * y * y + std::sqrt(4604660.15625) * x * y * y * y * y) + e_4 * (-std::sqrt(2046515.625) * x * x * x + std::sqrt(18418640.625) * x * y * y);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        pc_129[k] = e_0 * (-std::sqrt(1.48040771484375) * x * x * x * x * x * x * x * x * x * x * z + std::sqrt(213.1787109375) * x * x * x * x * x * x * x * x * y * y * z + std::sqrt(94.74609375) * x * x * x * x * x * x * x * x * z * z * z - std::sqrt(716.517333984375) * x * x * x * x * x * x * y * y * y * y * z - std::sqrt(16012.08984375) * x * x * x * x * x * x * y * y * z * z * z - std::sqrt(592.1630859375) * x * x * x * x * y * y * y * y * y * y * z + std::sqrt(116063.96484375) * x * x * x * x * y * y * y * y * z * z * z + std::sqrt(333.09173583984375) * x * x * y * y * y * y * y * y * y * y * z - std::sqrt(21317.87109375) * x * x * y * y * y * y * y * y * z * z * z) + e_1 * (-std::sqrt(119.91302490234375) * x * x * x * x * x * x * x * x * z - std::sqrt(17267.4755859375) * x * x * x * x * x * x * y * y * z + std::sqrt(21317.87109375) * x * x * x * x * x * x * z * z * z + std::sqrt(65285.980224609375) * x * x * x * x * y * y * y * y * z + std::sqrt(21317.87109375) * x * x * x * x * y * y * z * z * z - std::sqrt(5329.4677734375) * x * x * y * y * y * y * y * y * z - std::sqrt(21317.87109375) * x * x * y * y * y * y * z * z * z + std::sqrt(333.09173583984375) * y * y * y * y * y * y * y * y * z - std::sqrt(21317.87109375) * y * y * y * y * y * y * z * z * z) + e_2 * (std::sqrt(1364343.75) * x * x * x * x * z * z * z - std::sqrt(1364343.75) * y * y * y * y * z * z * z) + e_3 * (std::sqrt(1364343.75) * x * x * x * x * z + std::sqrt(5457375.0) * x * x * z * z * z - std::sqrt(1364343.75) * y * y * y * y * z - std::sqrt(5457375.0) * y * y * z * z * z) + e_4 * (std::sqrt(12279093.75) * x * x * z - std::sqrt(12279093.75) * y * y * z);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        pc_130[k] = e_0 * (std::sqrt(26.6473388671875) * x * x * x * x * x * x * x * x * x * z * z - std::sqrt(6821.71875) * x * x * x * x * x * x * x * y * y * z * z + std::sqrt(116075.80810546875) * x * x * x * x * x * y * y * y * y * z * z - std::sqrt(42635.7421875) * x * x * x * y * y * y * y * y * y * z * z + std::sqrt(666.1834716796875) * x * y * y * y * y * y * y * y * y * z * z) + e_1 * (std::sqrt(26.6473388671875) * x * x * x * x * x * x * x * x * x - std::sqrt(6821.71875) * x * x * x * x * x * x * x * y * y + std::sqrt(10658.935546875) * x * x * x * x * x * x * x * z * z + std::sqrt(116075.80810546875) * x * x * x * x * x * y * y * y * y + std::sqrt(95930.419921875) * x * x * x * x * x * y * y * z * z - std::sqrt(42635.7421875) * x * x * x * y * y * y * y * y * y + std::sqrt(95930.419921875) * x * x * x * y * y * y * y * z * z + std::sqrt(666.1834716796875) * x * y * y * y * y * y * y * y * y + std::sqrt(10658.935546875) * x * y * y * y * y * y * y * z * z) + e_2 * (std::sqrt(10658.935546875) * x * x * x * x * x * x * x + std::sqrt(95930.419921875) * x * x * x * x * x * y * y + std::sqrt(1534886.71875) * x * x * x * x * x * z * z + std::sqrt(95930.419921875) * x * x * x * y * y * y * y + std::sqrt(6139546.875) * x * x * x * y * y * z * z + std::sqrt(10658.935546875) * x * y * y * y * y * y * y + std::sqrt(1534886.71875) * x * y * y * y * y * z * z) + e_3 * (std::sqrt(1534886.71875) * x * x * x * x * x + std::sqrt(6139546.875) * x * x * x * y * y + std::sqrt(24558187.5) * x * x * x * z * z + std::sqrt(1534886.71875) * x * y * y * y * y + std::sqrt(24558187.5) * x * y * y * z * z) + e_4 * (std::sqrt(24558187.5) * x * x * x + std::sqrt(24558187.5) * x * y * y + std::sqrt(24558187.5) * x * z * z) + e_5 * (std::sqrt(24558187.5) * x);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        pc_131[k] = e_0 * (std::sqrt(2.66473388671875) * x * x * x * x * x * x * x * x * x * x * z - std::sqrt(1065.8935546875) * x * x * x * x * x * x * x * x * y * y * z + std::sqrt(32243.280029296875) * x * x * x * x * x * x * y * y * y * y * z - std::sqrt(26647.3388671875) * x * x * x * x * y * y * y * y * y * y * z + std::sqrt(1665.4586791992188) * x * x * y * y * y * y * y * y * y * y * z) + e_1 * (std::sqrt(1665.4586791992188) * x * x * x * x * x * x * x * x * z + std::sqrt(26647.3388671875) * x * x * x * x * x * x * y * y * z + std::sqrt(59956.512451171875) * x * x * x * x * y * y * y * y * z + std::sqrt(26647.3388671875) * x * x * y * y * y * y * y * y * z + std::sqrt(1665.4586791992188) * y * y * y * y * y * y * y * y * z) + e_2 * (std::sqrt(426357.421875) * x * x * x * x * x * x * z + std::sqrt(3837216.796875) * x * x * x * x * y * y * z + std::sqrt(3837216.796875) * x * x * y * y * y * y * z + std::sqrt(426357.421875) * y * y * y * y * y * y * z) + e_3 * (std::sqrt(15348867.1875) * x * x * x * x * z + std::sqrt(61395468.75) * x * x * y * y * z + std::sqrt(15348867.1875) * y * y * y * y * z) + e_4 * (std::sqrt(61395468.75) * x * x * z + std::sqrt(61395468.75) * y * y * z) + e_5 * (std::sqrt(9823275.0) * z);

        pc_132[k] = e_0 * (std::sqrt(5.5515289306640625) * x * x * x * x * x * x * x * x * x * x * y - std::sqrt(1604.391860961914) * x * x * x * x * x * x * x * x * y * y * y + std::sqrt(11341.995666503906) * x * x * x * x * x * x * y * y * y * y * y - std::sqrt(6417.567443847656) * x * x * x * x * y * y * y * y * y * y * y + std::sqrt(138.78822326660156) * x * x * y * y * y * y * y * y * y * y * y - std::sqrt(0.2220611572265625) * y * y * y * y * y * y * y * y * y * y * y) + e_1 * (-std::sqrt(199.85504150390625) * x * x * x * x * x * x * x * x * y - std::sqrt(3197.6806640625) * x * x * x * x * x * x * y * y * y - std::sqrt(7194.781494140625) * x * x * x * x * y * y * y * y * y - std::sqrt(3197.6806640625) * x * x * y * y * y * y * y * y * y - std::sqrt(199.85504150390625) * y * y * y * y * y * y * y * y * y) + e_2 * (-std::sqrt(79942.0166015625) * x * x * x * x * x * x * y - std::sqrt(719478.1494140625) * x * x * x * x * y * y * y - std::sqrt(719478.1494140625) * x * x * y * y * y * y * y - std::sqrt(79942.0166015625) * y * y * y * y * y * y * y) + e_3 * (-std::sqrt(5116289.0625) * x * x * x * x * y - std::sqrt(20465156.25) * x * x * y * y * y - std::sqrt(5116289.0625) * y * y * y * y * y) + e_4 * (-std::sqrt(46046601.5625) * x * x * y - std::sqrt(46046601.5625) * y * y * y) + e_5 * (-std::sqrt(29469825.0) * y);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        pc_133[k] = e_0 * (std::sqrt(35.52978515625) * x * x * x * x * x * x * x * x * x * y * z - std::sqrt(9095.625) * x * x * x * x * x * x * x * y * y * y * z + std::sqrt(31976.806640625) * x * x * x * x * x * y * y * y * y * y * z - std::sqrt(9095.625) * x * x * x * y * y * y * y * y * y * y * z + std::sqrt(35.52978515625) * x * y * y * y * y * y * y * y * y * y * z) + e_1 * (-std::sqrt(5116.2890625) * x * x * x * x * x * x * x * y * z - std::sqrt(46046.6015625) * x * x * x * x * x * y * y * y * z - std::sqrt(46046.6015625) * x * x * x * y * y * y * y * y * z - std::sqrt(5116.2890625) * x * y * y * y * y * y * y * y * z) + e_2 * (-std::sqrt(1151165.0390625) * x * x * x * x * x * y * z - std::sqrt(4604660.15625) * x * x * x * y * y * y * z - std::sqrt(1151165.0390625) * x * y * y * y * y * y * z) + e_3 * (-std::sqrt(32744250.0) * x * x * x * y * z - std::sqrt(32744250.0) * x * y * y * y * z) + e_4 * (-std::sqrt(73674562.5) * x * y * z);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        pc_134[k] = e_0 * (-std::sqrt(1.1103057861328125) * x * x * x * x * x * x * x * x * x * x * y + std::sqrt(228.1061553955078) * x * x * x * x * x * x * x * x * y * y * y + std::sqrt(71.0595703125) * x * x * x * x * x * x * x * x * y * z * z - std::sqrt(24.17999267578125) * x * x * x * x * x * x * y * y * y * y * y - std::sqrt(16706.89453125) * x * x * x * x * x * x * y * y * y * z * z - std::sqrt(217.61993408203125) * x * x * x * x * y * y * y * y * y * y * y + std::sqrt(28423.828125) * x * x * x * x * y * y * y * y * y * z * z + std::sqrt(35.65315246582031) * x * x * y * y * y * y * y * y * y * y * y - std::sqrt(2558.14453125) * x * x * y * y * y * y * y * y * y * z * z - std::sqrt(0.1233673095703125) * y * y * y * y * y * y * y * y * y * y * y + std::sqrt(7.8955078125) * y * y * y * y * y * y * y * y * y * z * z) + e_1 * (std::sqrt(39.97100830078125) * x * x * x * x * x * x * x * x * y + std::sqrt(59761.0986328125) * x * x * x * x * x * x * y * y * y - std::sqrt(23023.30078125) * x * x * x * x * x * x * y * z * z - std::sqrt(46206.485595703125) * x * x * x * x * y * y * y * y * y - std::sqrt(63953.61328125) * x * x * x * x * y * y * y * z * z + std::sqrt(5755.8251953125) * x * x * y * y * y * y * y * y * y - std::sqrt(2558.14453125) * x * x * y * y * y * y * y * z * z - std::sqrt(111.03057861328125) * y * y * y * y * y * y * y * y * y + std::sqrt(2558.14453125) * y * y * y * y * y * y * y * z * z) + e_2 * (std::sqrt(143895.6298828125) * x * x * x * x * x * x * y + std::sqrt(399710.0830078125) * x * x * x * x * y * y * y - std::sqrt(2302330.078125) * x * x * x * x * y * z * z + std::sqrt(15988.4033203125) * x * x * y * y * y * y * y - std::sqrt(1023257.8125) * x * x * y * y * y * z * z - std::sqrt(15988.4033203125) * y * y * y * y * y * y * y + std::sqrt(255814.453125) * y * y * y * y * y * z * z) + e_3 * (std::sqrt(4093031.25) * x * x * x * x * y + std::sqrt(1819125.0) * x * x * y * y * y - std::sqrt(16372125.0) * x * x * y * z * z - std::sqrt(454781.25) * y * y * y * y * y + std::sqrt(1819125.0) * y * y * y * z * z) + e_4 * (std::sqrt(9209320.3125) * x * x * y - std::sqrt(1023257.8125) * y * y * y);
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

        pc_135[k] = e_0 * (-std::sqrt(11.84326171875) * x * x * x * x * x * x * x * x * x * y * z + std::sqrt(2321.279296875) * x * x * x * x * x * x * x * y * y * y * z + std::sqrt(47.373046875) * x * x * x * x * x * x * x * y * z * z * z - std::sqrt(10658.935546875) * x * x * x * x * x * y * y * y * z * z * z - std::sqrt(2321.279296875) * x * x * x * y * y * y * y * y * y * y * z + std::sqrt(10658.935546875) * x * x * x * y * y * y * y * y * z * z * z + std::sqrt(11.84326171875) * x * y * y * y * y * y * y * y * y * y * z - std::sqrt(47.373046875) * x * y * y * y * y * y * y * y * z * z * z) + e_1 * (std::sqrt(1705.4296875) * x * x * x * x * x * x * x * y * z + std::sqrt(492869.1796875) * x * x * x * x * x * y * y * y * z - std::sqrt(27286.875) * x * x * x * x * x * y * z * z * z - std::sqrt(492869.1796875) * x * x * x * y * y * y * y * y * z - std::sqrt(1705.4296875) * x * y * y * y * y * y * y * y * z + std::sqrt(27286.875) * x * y * y * y * y * y * z * z * z) + e_2 * (std::sqrt(1534886.71875) * x * x * x * x * x * y * z - std::sqrt(682171.875) * x * x * x * y * z * z * z - std::sqrt(1534886.71875) * x * y * y * y * y * y * z + std::sqrt(682171.875) * x * y * y * y * z * z * z) + e_3 * (std::sqrt(10914750.0) * x * x * x * y * z - std::sqrt(10914750.0) * x * y * y * y * z);
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

        pc_136[k] = e_0 * (std::sqrt(0.105743408203125) * x * x * x * x * x * x * x * x * x * x * y - std::sqrt(17.870635986328125) * x * x * x * x * x * x * x * x * y * y * y - std::sqrt(15.22705078125) * x * x * x * x * x * x * x * x * y * z * z - std::sqrt(20.7257080078125) * x * x * x * x * x * x * y * y * y * y * y + std::sqrt(2984.501953125) * x * x * x * x * x * x * y * y * y * z * z + std::sqrt(6.767578125) * x * x * x * x * x * x * y * z * z * z * z + std::sqrt(20.7257080078125) * x * x * x * x * y * y * y * y * y * y * y - std::sqrt(1522.705078125) * x * x * x * x * y * y * y * z * z * z * z + std::sqrt(17.870635986328125) * x * x * y * y * y * y * y * y * y * y * y - std::sqrt(2984.501953125) * x * x * y * y * y * y * y * y * y * z * z + std::sqrt(1522.705078125) * x * x * y * y * y * y * y * z * z * z * z - std::sqrt(0.105743408203125) * y * y * y * y * y * y * y * y * y * y * y + std::sqrt(15.22705078125) * y * y * y * y * y * y * y * y * y * z * z - std::sqrt(6.767578125) * y * y * y * y * y * y * y * z * z * z * z) + e_1 * (-std::sqrt(3.8067626953125) * x * x * x * x * x * x * x * x * y - std::sqrt(11938.0078125) * x * x * x * x * x * x * y * y * y + std::sqrt(4933.564453125) * x * x * x * x * x * x * y * z * z + std::sqrt(746.12548828125) * x * x * x * x * y * y * y * y * y + std::sqrt(342608.642578125) * x * x * x * x * y * y * y * z * z - std::sqrt(6090.8203125) * x * x * x * x * y * z * z * z * z + std::sqrt(15592.5) * x * x * y * y * y * y * y * y * y - std::sqrt(833772.392578125) * x * x * y * y * y * y * y * z * z + std::sqrt(24363.28125) * x * x * y * y * y * z * z * z * z - std::sqrt(95.1690673828125) * y * y * y * y * y * y * y * y * y + std::sqrt(4933.564453125) * y * y * y * y * y * y * y * z * z - std::sqrt(243.6328125) * y * y * y * y * y * z * z * z * z) + e_2 * (-std::sqrt(24363.28125) * x * x * x * x * x * x * y - std::sqrt(152270.5078125) * x * x * x * x * y * y * y + std::sqrt(1370434.5703125) * x * x * x * x * y * z * z + std::sqrt(877078.125) * x * x * y * y * y * y * y - std::sqrt(5481738.28125) * x * x * y * y * y * z * z - std::sqrt(6090.8203125) * y * y * y * y * y * y * y + std::sqrt(54817.3828125) * y * y * y * y * y * z * z) + e_3 * (-std::sqrt(609082.03125) * x * x * x * x * y + std::sqrt(2436328.125) * x * x * y * y * y - std::sqrt(24363.28125) * y * y * y * y * y);
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

        pc_137[k] = e_0 * (std::sqrt(1.586151123046875) * x * x * x * x * x * x * x * x * x * x * z - std::sqrt(268.0595397949219) * x * x * x * x * x * x * x * x * y * y * z - std::sqrt(11.279296875) * x * x * x * x * x * x * x * x * z * z * z - std::sqrt(310.8856201171875) * x * x * x * x * x * x * y * y * y * y * z + std::sqrt(2210.7421875) * x * x * x * x * x * x * y * y * z * z * z + std::sqrt(0.451171875) * x * x * x * x * x * x * z * z * z * z * z + std::sqrt(310.8856201171875) * x * x * x * x * y * y * y * y * y * y * z - std::sqrt(101.513671875) * x * x * x * x * y * y * z * z * z * z * z + std::sqrt(268.0595397949219) * x * x * y * y * y * y * y * y * y * y * z - std::sqrt(2210.7421875) * x * x * y * y * y * y * y * y * z * z * z + std::sqrt(101.513671875) * x * x * y * y * y * y * z * z * z * z * z - std::sqrt(1.586151123046875) * y * y * y * y * y * y * y * y * y * y * z + std::sqrt(11.279296875) * y * y * y * y * y * y * y * y * z * z * z - std::sqrt(0.451171875) * y * y * y * y * y * y * z * z * z * z * z) + e_1 * (std::sqrt(913.623046875) * x * x * x * x * x * x * x * x * z - std::sqrt(179070.1171875) * x * x * x * x * x * x * y * y * z - std::sqrt(1624.21875) * x * x * x * x * x * x * z * z * z + std::sqrt(365449.21875) * x * x * x * x * y * y * z * z * z + std::sqrt(179070.1171875) * x * x * y * y * y * y * y * y * z - std::sqrt(365449.21875) * x * x * y * y * y * y * z * z * z - std::sqrt(913.623046875) * y * y * y * y * y * y * y * y * z + std::sqrt(1624.21875) * y * y * y * y * y * y * z * z * z) + e_2 * (std::sqrt(22840.576171875) * x * x * x * x * x * x * z - std::sqrt(5139129.638671875) * x * x * x * x * y * y * z + std::sqrt(5139129.638671875) * x * x * y * y * y * y * z - std::sqrt(22840.576171875) * y * y * y * y * y * y * z);
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

        pc_138[k] = e_0 * (std::sqrt(0.105743408203125) * x * x * x * x * x * x * x * x * x * x * x - std::sqrt(17.870635986328125) * x * x * x * x * x * x * x * x * x * y * y - std::sqrt(15.22705078125) * x * x * x * x * x * x * x * x * x * z * z - std::sqrt(20.7257080078125) * x * x * x * x * x * x * x * y * y * y * y + std::sqrt(2984.501953125) * x * x * x * x * x * x * x * y * y * z * z + std::sqrt(6.767578125) * x * x * x * x * x * x * x * z * z * z * z + std::sqrt(20.7257080078125) * x * x * x * x * x * y * y * y * y * y * y - std::sqrt(1522.705078125) * x * x * x * x * x * y * y * z * z * z * z + std::sqrt(17.870635986328125) * x * x * x * y * y * y * y * y * y * y * y - std::sqrt(2984.501953125) * x * x * x * y * y * y * y * y * y * z * z + std::sqrt(1522.705078125) * x * x * x * y * y * y * y * z * z * z * z - std::sqrt(0.105743408203125) * x * y * y * y * y * y * y * y * y * y * y + std::sqrt(15.22705078125) * x * y * y * y * y * y * y * y * y * z * z - std::sqrt(6.767578125) * x * y * y * y * y * y * y * z * z * z * z) + e_1 * (std::sqrt(95.1690673828125) * x * x * x * x * x * x * x * x * x - std::sqrt(15592.5) * x * x * x * x * x * x * x * y * y - std::sqrt(4933.564453125) * x * x * x * x * x * x * x * z * z - std::sqrt(746.12548828125) * x * x * x * x * x * y * y * y * y + std::sqrt(833772.392578125) * x * x * x * x * x * y * y * z * z + std::sqrt(243.6328125) * x * x * x * x * x * z * z * z * z + std::sqrt(11938.0078125) * x * x * x * y * y * y * y * y * y - std::sqrt(342608.642578125) * x * x * x * y * y * y * y * z * z - std::sqrt(24363.28125) * x * x * x * y * y * z * z * z * z + std::sqrt(3.8067626953125) * x * y * y * y * y * y * y * y * y - std::sqrt(4933.564453125) * x * y * y * y * y * y * y * z * z + std::sqrt(6090.8203125) * x * y * y * y * y * z * z * z * z) + e_2 * (std::sqrt(6090.8203125) * x * x * x * x * x * x * x - std::sqrt(877078.125) * x * x * x * x * x * y * y - std::sqrt(54817.3828125) * x * x * x * x * x * z * z + std::sqrt(152270.5078125) * x * x * x * y * y * y * y + std::sqrt(5481738.28125) * x * x * x * y * y * z * z + std::sqrt(24363.28125) * x * y * y * y * y * y * y - std::sqrt(1370434.5703125) * x * y * y * y * y * z * z) + e_3 * (std::sqrt(24363.28125) * x * x * x * x * x - std::sqrt(2436328.125) * x * x * x * y * y + std::sqrt(609082.03125) * x * y * y * y * y);
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

        pc_139[k] = e_0 * (-std::sqrt(2.9608154296875) * x * x * x * x * x * x * x * x * x * x * z + std::sqrt(666.1834716796875) * x * x * x * x * x * x * x * x * y * y * z + std::sqrt(11.84326171875) * x * x * x * x * x * x * x * x * z * z * z - std::sqrt(580.31982421875) * x * x * x * x * x * x * y * y * y * y * z - std::sqrt(3031.875) * x * x * x * x * x * x * y * y * z * z * z - std::sqrt(580.31982421875) * x * x * x * x * y * y * y * y * y * y * z + std::sqrt(10658.935546875) * x * x * x * x * y * y * y * y * z * z * z + std::sqrt(666.1834716796875) * x * x * y * y * y * y * y * y * y * y * z - std::sqrt(3031.875) * x * x * y * y * y * y * y * y * z * z * z - std::sqrt(2.9608154296875) * y * y * y * y * y * y * y * y * y * y * z + std::sqrt(11.84326171875) * y * y * y * y * y * y * y * y * z * z * z) + e_1 * (-std::sqrt(1705.4296875) * x * x * x * x * x * x * x * x * z + std::sqrt(170542.96875) * x * x * x * x * x * x * y * y * z + std::sqrt(1705.4296875) * x * x * x * x * x * x * z * z * z - std::sqrt(170542.96875) * x * x * x * x * y * y * y * y * z - std::sqrt(42635.7421875) * x * x * x * x * y * y * z * z * z + std::sqrt(170542.96875) * x * x * y * y * y * y * y * y * z - std::sqrt(42635.7421875) * x * x * y * y * y * y * z * z * z - std::sqrt(1705.4296875) * y * y * y * y * y * y * y * y * z + std::sqrt(1705.4296875) * y * y * y * y * y * y * z * z * z) + e_2 * (-std::sqrt(95930.419921875) * x * x * x * x * x * x * z + std::sqrt(2398260.498046875) * x * x * x * x * y * y * z + std::sqrt(42635.7421875) * x * x * x * x * z * z * z + std::sqrt(2398260.498046875) * x * x * y * y * y * y * z - std::sqrt(1534886.71875) * x * x * y * y * z * z * z - std::sqrt(95930.419921875) * y * y * y * y * y * y * z + std::sqrt(42635.7421875) * y * y * y * y * z * z * z) + e_3 * (-std::sqrt(682171.875) * x * x * x * x * z + std::sqrt(24558187.5) * x * x * y * y * z - std::sqrt(682171.875) * y * y * y * y * z);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        pc_140[k] = e_0 * (-std::sqrt(0.1233673095703125) * x * x * x * x * x * x * x * x * x * x * x + std::sqrt(35.65315246582031) * x * x * x * x * x * x * x * x * x * y * y + std::sqrt(7.8955078125) * x * x * x * x * x * x * x * x * x * z * z - std::sqrt(217.61993408203125) * x * x * x * x * x * x * x * y * y * y * y - std::sqrt(2558.14453125) * x * x * x * x * x * x * x * y * y * z * z - std::sqrt(24.17999267578125) * x * x * x * x * x * y * y * y * y * y * y + std::sqrt(28423.828125) * x * x * x * x * x * y * y * y * y * z * z + std::sqrt(228.1061553955078) * x * x * x * y * y * y * y * y * y * y * y - std::sqrt(16706.89453125) * x * x * x * y * y * y * y * y * y * z * z - std::sqrt(1.1103057861328125) * x * y * y * y * y * y * y * y * y * y * y + std::sqrt(71.0595703125) * x * y * y * y * y * y * y * y * y * z * z) + e_1 * (-std::sqrt(111.03057861328125) * x * x * x * x * x * x * x * x * x + std::sqrt(5755.8251953125) * x * x * x * x * x * x * x * y * y + std::sqrt(2558.14453125) * x * x * x * x * x * x * x * z * z - std::sqrt(46206.485595703125) * x * x * x * x * x * y * y * y * y - std::sqrt(2558.14453125) * x * x * x * x * x * y * y * z * z + std::sqrt(59761.0986328125) * x * x * x * y * y * y * y * y * y - std::sqrt(63953.61328125) * x * x * x * y * y * y * y * z * z + std::sqrt(39.97100830078125) * x * y * y * y * y * y * y * y * y - std::sqrt(23023.30078125) * x * y * y * y * y * y * y * z * z) + e_2 * (-std::sqrt(15988.4033203125) * x * x * x * x * x * x * x + std::sqrt(15988.4033203125) * x * x * x * x * x * y * y + std::sqrt(255814.453125) * x * x * x * x * x * z * z + std::sqrt(399710.0830078125) * x * x * x * y * y * y * y - std::sqrt(1023257.8125) * x * x * x * y * y * z * z + std::sqrt(143895.6298828125) * x * y * y * y * y * y * y - std::sqrt(2302330.078125) * x * y * y * y * y * z * z) + e_3 * (-std::sqrt(454781.25) * x * x * x * x * x + std::sqrt(1819125.0) * x * x * x * y * y + std::sqrt(1819125.0) * x * x * x * z * z + std::sqrt(4093031.25) * x * y * y * y * y - std::sqrt(16372125.0) * x * y * y * z * z) + e_4 * (-std::sqrt(1023257.8125) * x * x * x + std::sqrt(9209320.3125) * x * y * y);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ab_x, ab_y, ab_z : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];
        const auto z = ab_z[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        pc_141[k] = e_0 * (std::sqrt(2.220611572265625) * x * x * x * x * x * x * x * x * x * x * z - std::sqrt(979.2897033691406) * x * x * x * x * x * x * x * x * y * y * z + std::sqrt(24950.791625976562) * x * x * x * x * x * x * y * y * y * y * z - std::sqrt(24950.791625976562) * x * x * x * x * y * y * y * y * y * y * z + std::sqrt(979.2897033691406) * x * x * y * y * y * y * y * y * y * y * z - std::sqrt(2.220611572265625) * y * y * y * y * y * y * y * y * y * y * z) + e_1 * (std::sqrt(1279.072265625) * x * x * x * x * x * x * x * x * z + std::sqrt(5116.2890625) * x * x * x * x * x * x * y * y * z - std::sqrt(5116.2890625) * x * x * y * y * y * y * y * y * z - std::sqrt(1279.072265625) * y * y * y * y * y * y * y * y * z) + e_2 * (std::sqrt(287791.259765625) * x * x * x * x * x * x * z + std::sqrt(287791.259765625) * x * x * x * x * y * y * z - std::sqrt(287791.259765625) * x * x * y * y * y * y * z - std::sqrt(287791.259765625) * y * y * y * y * y * y * z) + e_3 * (std::sqrt(8186062.5) * x * x * x * x * z - std::sqrt(8186062.5) * y * y * y * y * z) + e_4 * (std::sqrt(18418640.625) * x * x * z - std::sqrt(18418640.625) * y * y * z);

        pc_142[k] = e_0 * (std::sqrt(0.2220611572265625) * x * x * x * x * x * x * x * x * x * x * x - std::sqrt(138.78822326660156) * x * x * x * x * x * x * x * x * x * y * y + std::sqrt(6417.567443847656) * x * x * x * x * x * x * x * y * y * y * y - std::sqrt(11341.995666503906) * x * x * x * x * x * y * y * y * y * y * y + std::sqrt(1604.391860961914) * x * x * x * y * y * y * y * y * y * y * y - std::sqrt(5.5515289306640625) * x * y * y * y * y * y * y * y * y * y * y) + e_1 * (std::sqrt(199.85504150390625) * x * x * x * x * x * x * x * x * x + std::sqrt(3197.6806640625) * x * x * x * x * x * x * x * y * y + std::sqrt(7194.781494140625) * x * x * x * x * x * y * y * y * y + std::sqrt(3197.6806640625) * x * x * x * y * y * y * y * y * y + std::sqrt(199.85504150390625) * x * y * y * y * y * y * y * y * y) + e_2 * (std::sqrt(79942.0166015625) * x * x * x * x * x * x * x + std::sqrt(719478.1494140625) * x * x * x * x * x * y * y + std::sqrt(719478.1494140625) * x * x * x * y * y * y * y + std::sqrt(79942.0166015625) * x * y * y * y * y * y * y) + e_3 * (std::sqrt(5116289.0625) * x * x * x * x * x + std::sqrt(20465156.25) * x * x * x * y * y + std::sqrt(5116289.0625) * x * y * y * y * y) + e_4 * (std::sqrt(46046601.5625) * x * x * x + std::sqrt(46046601.5625) * x * y * y) + e_5 * (std::sqrt(29469825.0) * x);
    }

    // NOTE: the atom pairs beyond the reach of every pair of primitives have no
    // contribution and are set to zero.

    for (size_t m = 0; m < 143; m++)
    {
        auto *pv = values + m * nvalues;

        std::fill(pv + nmax, pv + nvalues, 0.0);
    }
}

}  // namespace simdovl
