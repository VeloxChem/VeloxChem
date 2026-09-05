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



#include "SimdOverlapRecGI.hpp"

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
compute_gi_overlap(double               *values,
                   const size_t          nvalues,
                   const CBasisFunction &bra,
                   const CBasisFunction &ket,
                   const CSimdMatrix    &coordinates,
                   const double          threshold) -> void
{
    if ((bra.get_angular_momentum() != 4) || (ket.get_angular_momentum() != 6))
    {
        errors::assertMsgCritical(
            false, std::string("SimdOverlapRecGI.compute_gi_overlap: Basis functions must be of angular momenta four and six"));
    }

    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("SimdOverlapRecGI.compute_gi_overlap: Number of values exceeds number of atom pairs"));
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

    auto buffer = simdfunc::make_primitive_buffer(dimensions, 5);

    if (buffer.number_of_columns() == 0)
    {
        std::fill(values, values + 117 * nvalues, 0.0);

        return;
    }

    const auto nmax = buffer.number_of_columns();

    auto *pe_0 = buffer.data(0);
    auto *pe_1 = buffer.data(1);
    auto *pe_2 = buffer.data(2);
    auto *pe_3 = buffer.data(3);
    auto *pe_4 = buffer.data(4);

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

        const auto f_0 = fbase * fal * fal * fal * fal * fbe * fbe * fbe * fbe * fbe * fbe;

        const auto f_1 = fbase * fal * fal * fal * fbe * fbe * fbe * fbe * fbe * fh;

        const auto f_2 = fbase * fal * fal * fbe * fbe * fbe * fbe * fh * fh;

        const auto f_3 = fbase * fal * fbe * fbe * fbe * fh * fh * fh;

        const auto f_4 = fbase * fbe * fbe * fh * fh * fh * fh;

        // NOTE: the exponential depends on the pair of primitives alone, so it is
        // evaluated once and shared by the prefactors of all terms.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ab_2 : simd::cache_line_size())
        for (size_t k = 0; k < ncols; k++)
        {
            const auto fss = std::exp(-fmu * ab_2[k]);

            pe_0[k] += f_0 * fss;
            pe_1[k] += f_1 * fss;
            pe_2[k] += f_2 * fss;
            pe_3[k] += f_3 * fss;
            pe_4[k] += f_4 * fss;
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

    // NOTE: the components are formed in 94 loops, as the vectorizer runs out
    // of registers with all of them in one. Only the prefactors and the vector
    // between the atoms are loaded by more than one loop.

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

        pc_0[k] = e_0 * (std::sqrt(142.119140625) * x * x * x * x * x * x * x * x * y * y - std::sqrt(2668.681640625) * x * x * x * x * x * x * y * y * y * y + std::sqrt(2668.681640625) * x * x * x * x * y * y * y * y * y * y - std::sqrt(142.119140625) * x * x * y * y * y * y * y * y * y * y) + e_1 * (std::sqrt(142.119140625) * x * x * x * x * x * x * x * x + std::sqrt(568.4765625) * x * x * x * x * x * x * y * y - std::sqrt(568.4765625) * x * x * y * y * y * y * y * y - std::sqrt(142.119140625) * y * y * y * y * y * y * y * y) + e_2 * (std::sqrt(31976.806640625) * x * x * x * x * x * x + std::sqrt(31976.806640625) * x * x * x * x * y * y - std::sqrt(31976.806640625) * x * x * y * y * y * y - std::sqrt(31976.806640625) * y * y * y * y * y * y) + e_3 * (std::sqrt(909562.5) * x * x * x * x - std::sqrt(909562.5) * y * y * y * y) + e_4 * (std::sqrt(2046515.625) * x * x - std::sqrt(2046515.625) * y * y);

        pc_1[k] = e_0 * (std::sqrt(1184.326171875) * x * x * x * x * x * x * x * y * y * z - std::sqrt(10658.935546875) * x * x * x * x * x * y * y * y * y * z + std::sqrt(5732.138671875) * x * x * x * y * y * y * y * y * y * z - std::sqrt(47.373046875) * x * y * y * y * y * y * y * y * y * z) + e_1 * (std::sqrt(1184.326171875) * x * x * x * x * x * x * x * z + std::sqrt(10658.935546875) * x * x * x * x * x * y * y * z + std::sqrt(10658.935546875) * x * x * x * y * y * y * y * z + std::sqrt(1184.326171875) * x * y * y * y * y * y * y * z) + e_2 * (std::sqrt(170542.96875) * x * x * x * x * x * z + std::sqrt(682171.875) * x * x * x * y * y * z + std::sqrt(170542.96875) * x * y * y * y * y * z) + e_3 * (std::sqrt(2728687.5) * x * x * x * z + std::sqrt(2728687.5) * x * y * y * z) + e_4 * (std::sqrt(2728687.5) * x * z);
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

        pc_2[k] = e_0 * (-std::sqrt(34.453125) * x * x * x * x * x * x * x * x * y * y + std::sqrt(34.453125) * x * x * x * x * x * x * y * y * y * y + std::sqrt(3445.3125) * x * x * x * x * x * x * y * y * z * z + std::sqrt(34.453125) * x * x * x * x * y * y * y * y * y * y - std::sqrt(13781.25) * x * x * x * x * y * y * y * y * z * z - std::sqrt(34.453125) * x * x * y * y * y * y * y * y * y * y + std::sqrt(3445.3125) * x * x * y * y * y * y * y * y * z * z) + e_1 * (-std::sqrt(34.453125) * x * x * x * x * x * x * x * x - std::sqrt(4961.25) * x * x * x * x * x * x * y * y + std::sqrt(3445.3125) * x * x * x * x * x * x * z * z + std::sqrt(3445.3125) * x * x * x * x * y * y * y * y + std::sqrt(31007.8125) * x * x * x * x * y * y * z * z - std::sqrt(4961.25) * x * x * y * y * y * y * y * y + std::sqrt(31007.8125) * x * x * y * y * y * y * z * z - std::sqrt(34.453125) * y * y * y * y * y * y * y * y + std::sqrt(3445.3125) * y * y * y * y * y * y * z * z) + e_2 * (-std::sqrt(7751.953125) * x * x * x * x * x * x - std::sqrt(69767.578125) * x * x * x * x * y * y + std::sqrt(279070.3125) * x * x * x * x * z * z - std::sqrt(69767.578125) * x * x * y * y * y * y + std::sqrt(1116281.25) * x * x * y * y * z * z - std::sqrt(7751.953125) * y * y * y * y * y * y + std::sqrt(279070.3125) * y * y * y * y * z * z) + e_3 * (-std::sqrt(124031.25) * x * x * x * x - std::sqrt(496125.0) * x * x * y * y + std::sqrt(1984500.0) * x * x * z * z - std::sqrt(124031.25) * y * y * y * y + std::sqrt(1984500.0) * y * y * z * z) + e_4 * (-std::sqrt(124031.25) * x * x - std::sqrt(124031.25) * y * y + std::sqrt(496125.0) * z * z);
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

        pc_3[k] = e_0 * (-std::sqrt(581.396484375) * x * x * x * x * x * x * x * y * y * z + std::sqrt(64.599609375) * x * x * x * x * x * y * y * y * y * z + std::sqrt(4134.375) * x * x * x * x * x * y * y * z * z * z + std::sqrt(581.396484375) * x * x * x * y * y * y * y * y * y * z - std::sqrt(7350.0) * x * x * x * y * y * y * y * z * z * z - std::sqrt(64.599609375) * x * y * y * y * y * y * y * y * y * z + std::sqrt(459.375) * x * y * y * y * y * y * y * z * z * z) + e_1 * (-std::sqrt(581.396484375) * x * x * x * x * x * x * x * z - std::sqrt(70348.974609375) * x * x * x * x * x * y * y * z + std::sqrt(4134.375) * x * x * x * x * x * z * z * z + std::sqrt(34173.193359375) * x * x * x * y * y * y * y * z + std::sqrt(16537.5) * x * x * x * y * y * z * z * z - std::sqrt(7816.552734375) * x * y * y * y * y * y * y * z + std::sqrt(4134.375) * x * y * y * y * y * z * z * z) + e_2 * (-std::sqrt(83721.09375) * x * x * x * x * x * z - std::sqrt(334884.375) * x * x * x * y * y * z + std::sqrt(148837.5) * x * x * x * z * z * z - std::sqrt(83721.09375) * x * y * y * y * y * z + std::sqrt(148837.5) * x * y * y * z * z * z) + e_3 * (-std::sqrt(595350.0) * x * x * x * z - std::sqrt(595350.0) * x * y * y * z + std::sqrt(264600.0) * x * z * z * z) + e_4 * (-std::sqrt(148837.5) * x * z);
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

        pc_4[k] = e_0 * (std::sqrt(7.177734375) * x * x * x * x * x * x * x * x * y * y + std::sqrt(7.177734375) * x * x * x * x * x * x * y * y * y * y - std::sqrt(1837.5) * x * x * x * x * x * x * y * y * z * z - std::sqrt(7.177734375) * x * x * x * x * y * y * y * y * y * y + std::sqrt(1837.5) * x * x * x * x * y * y * z * z * z * z - std::sqrt(7.177734375) * x * x * y * y * y * y * y * y * y * y + std::sqrt(1837.5) * x * x * y * y * y * y * y * y * z * z - std::sqrt(1837.5) * x * x * y * y * y * y * z * z * z * z) + e_1 * (std::sqrt(7.177734375) * x * x * x * x * x * x * x * x + std::sqrt(2325.5859375) * x * x * x * x * x * x * y * y - std::sqrt(1837.5) * x * x * x * x * x * x * z * z - std::sqrt(148837.5) * x * x * x * x * y * y * z * z + std::sqrt(1837.5) * x * x * x * x * z * z * z * z - std::sqrt(2325.5859375) * x * x * y * y * y * y * y * y + std::sqrt(148837.5) * x * x * y * y * y * y * z * z - std::sqrt(7.177734375) * y * y * y * y * y * y * y * y + std::sqrt(1837.5) * y * y * y * y * y * y * z * z - std::sqrt(1837.5) * y * y * y * y * z * z * z * z) + e_2 * (std::sqrt(1614.990234375) * x * x * x * x * x * x + std::sqrt(28488.427734375) * x * x * x * x * y * y - std::sqrt(148837.5) * x * x * x * x * z * z - std::sqrt(28488.427734375) * x * x * y * y * y * y + std::sqrt(16537.5) * x * x * z * z * z * z - std::sqrt(1614.990234375) * y * y * y * y * y * y + std::sqrt(148837.5) * y * y * y * y * z * z - std::sqrt(16537.5) * y * y * z * z * z * z) + e_3 * (std::sqrt(16537.5) * x * x * x * x - std::sqrt(264600.0) * x * x * z * z - std::sqrt(16537.5) * y * y * y * y + std::sqrt(264600.0) * y * y * z * z) + e_4 * (std::sqrt(4134.375) * x * x - std::sqrt(4134.375) * y * y);
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

        pc_5[k] = e_0 * (std::sqrt(71.77734375) * x * x * x * x * x * x * x * y * y * z + std::sqrt(71.77734375) * x * x * x * x * x * y * y * y * y * z - std::sqrt(1148.4375) * x * x * x * x * x * y * y * z * z * z - std::sqrt(71.77734375) * x * x * x * y * y * y * y * y * y * z + std::sqrt(183.75) * x * x * x * y * y * z * z * z * z * z - std::sqrt(71.77734375) * x * y * y * y * y * y * y * y * y * z + std::sqrt(1148.4375) * x * y * y * y * y * y * y * z * z * z - std::sqrt(183.75) * x * y * y * y * y * z * z * z * z * z) + e_1 * (std::sqrt(71.77734375) * x * x * x * x * x * x * x * z + std::sqrt(16149.90234375) * x * x * x * x * x * y * y * z - std::sqrt(1148.4375) * x * x * x * x * x * z * z * z - std::sqrt(1794.43359375) * x * x * x * y * y * y * y * z - std::sqrt(41343.75) * x * x * x * y * y * z * z * z + std::sqrt(183.75) * x * x * x * z * z * z * z * z - std::sqrt(25911.62109375) * x * y * y * y * y * y * y * z + std::sqrt(138960.9375) * x * y * y * y * y * z * z * z - std::sqrt(1653.75) * x * y * y * z * z * z * z * z) + e_2 * (std::sqrt(10335.9375) * x * x * x * x * x * z + std::sqrt(41343.75) * x * x * x * y * y * z - std::sqrt(41343.75) * x * x * x * z * z * z - std::sqrt(506460.9375) * x * y * y * y * y * z + std::sqrt(372093.75) * x * y * y * z * z * z) + e_3 * (std::sqrt(41343.75) * x * x * x * z - std::sqrt(372093.75) * x * y * y * z);
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

        pc_6[k] = e_0 * (-std::sqrt(0.8544921875) * x * x * x * x * x * x * x * x * x * y - std::sqrt(3.41796875) * x * x * x * x * x * x * x * y * y * y + std::sqrt(276.85546875) * x * x * x * x * x * x * x * y * z * z + std::sqrt(276.85546875) * x * x * x * x * x * y * y * y * z * z - std::sqrt(492.1875) * x * x * x * x * x * y * z * z * z * z + std::sqrt(3.41796875) * x * x * x * y * y * y * y * y * y * y - std::sqrt(276.85546875) * x * x * x * y * y * y * y * y * z * z + std::sqrt(8.75) * x * x * x * y * z * z * z * z * z * z + std::sqrt(0.8544921875) * x * y * y * y * y * y * y * y * y * y - std::sqrt(276.85546875) * x * y * y * y * y * y * y * y * z * z + std::sqrt(492.1875) * x * y * y * y * y * y * z * z * z * z - std::sqrt(8.75) * x * y * y * y * z * z * z * z * z * z) + e_1 * (-std::sqrt(492.1875) * x * x * x * x * x * x * x * y - std::sqrt(492.1875) * x * x * x * x * x * y * y * y + std::sqrt(70875.0) * x * x * x * x * x * y * z * z + std::sqrt(492.1875) * x * x * x * y * y * y * y * y - std::sqrt(31500.0) * x * x * x * y * z * z * z * z + std::sqrt(492.1875) * x * y * y * y * y * y * y * y - std::sqrt(70875.0) * x * y * y * y * y * y * z * z + std::sqrt(31500.0) * x * y * y * y * z * z * z * z) + e_2 * (-std::sqrt(17718.75) * x * x * x * x * x * y + std::sqrt(637875.0) * x * x * x * y * z * z + std::sqrt(17718.75) * x * y * y * y * y * y - std::sqrt(637875.0) * x * y * y * y * z * z) + e_3 * (-std::sqrt(31500.0) * x * x * x * y + std::sqrt(31500.0) * x * y * y * y);
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

        pc_7[k] = e_0 * (std::sqrt(71.77734375) * x * x * x * x * x * x * x * x * y * z + std::sqrt(71.77734375) * x * x * x * x * x * x * y * y * y * z - std::sqrt(1148.4375) * x * x * x * x * x * x * y * z * z * z - std::sqrt(71.77734375) * x * x * x * x * y * y * y * y * y * z + std::sqrt(183.75) * x * x * x * x * y * z * z * z * z * z - std::sqrt(71.77734375) * x * x * y * y * y * y * y * y * y * z + std::sqrt(1148.4375) * x * x * y * y * y * y * y * z * z * z - std::sqrt(183.75) * x * x * y * y * y * z * z * z * z * z) + e_1 * (std::sqrt(25911.62109375) * x * x * x * x * x * x * y * z + std::sqrt(1794.43359375) * x * x * x * x * y * y * y * z - std::sqrt(138960.9375) * x * x * x * x * y * z * z * z - std::sqrt(16149.90234375) * x * x * y * y * y * y * y * z + std::sqrt(41343.75) * x * x * y * y * y * z * z * z + std::sqrt(1653.75) * x * x * y * z * z * z * z * z - std::sqrt(71.77734375) * y * y * y * y * y * y * y * z + std::sqrt(1148.4375) * y * y * y * y * y * z * z * z - std::sqrt(183.75) * y * y * y * z * z * z * z * z) + e_2 * (std::sqrt(506460.9375) * x * x * x * x * y * z - std::sqrt(41343.75) * x * x * y * y * y * z - std::sqrt(372093.75) * x * x * y * z * z * z - std::sqrt(10335.9375) * y * y * y * y * y * z + std::sqrt(41343.75) * y * y * y * z * z * z) + e_3 * (std::sqrt(372093.75) * x * x * y * z - std::sqrt(41343.75) * y * y * y * z);
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

        pc_8[k] = e_0 * (std::sqrt(1.79443359375) * x * x * x * x * x * x * x * x * x * y - std::sqrt(459.375) * x * x * x * x * x * x * x * y * z * z - std::sqrt(7.177734375) * x * x * x * x * x * y * y * y * y * y + std::sqrt(459.375) * x * x * x * x * x * y * y * y * z * z + std::sqrt(459.375) * x * x * x * x * x * y * z * z * z * z + std::sqrt(459.375) * x * x * x * y * y * y * y * y * z * z - std::sqrt(1837.5) * x * x * x * y * y * y * z * z * z * z + std::sqrt(1.79443359375) * x * y * y * y * y * y * y * y * y * y - std::sqrt(459.375) * x * y * y * y * y * y * y * y * z * z + std::sqrt(459.375) * x * y * y * y * y * y * z * z * z * z) + e_1 * (std::sqrt(717.7734375) * x * x * x * x * x * x * x * y - std::sqrt(28.7109375) * x * x * x * x * x * y * y * y - std::sqrt(66150.0) * x * x * x * x * x * y * z * z - std::sqrt(28.7109375) * x * x * x * y * y * y * y * y + std::sqrt(29400.0) * x * x * x * y * y * y * z * z + std::sqrt(7350.0) * x * x * x * y * z * z * z * z + std::sqrt(717.7734375) * x * y * y * y * y * y * y * y - std::sqrt(66150.0) * x * y * y * y * y * y * z * z + std::sqrt(7350.0) * x * y * y * y * z * z * z * z) + e_2 * (std::sqrt(20930.2734375) * x * x * x * x * x * y + std::sqrt(1033.59375) * x * x * x * y * y * y - std::sqrt(595350.0) * x * x * x * y * z * z + std::sqrt(20930.2734375) * x * y * y * y * y * y - std::sqrt(595350.0) * x * y * y * y * z * z + std::sqrt(66150.0) * x * y * z * z * z * z) + e_3 * (std::sqrt(66150.0) * x * x * x * y + std::sqrt(66150.0) * x * y * y * y - std::sqrt(1058400.0) * x * y * z * z) + e_4 * (std::sqrt(16537.5) * x * y);
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

        pc_9[k] = e_0 * (-std::sqrt(64.599609375) * x * x * x * x * x * x * x * x * y * z + std::sqrt(581.396484375) * x * x * x * x * x * x * y * y * y * z + std::sqrt(459.375) * x * x * x * x * x * x * y * z * z * z + std::sqrt(64.599609375) * x * x * x * x * y * y * y * y * y * z - std::sqrt(7350.0) * x * x * x * x * y * y * y * z * z * z - std::sqrt(581.396484375) * x * x * y * y * y * y * y * y * y * z + std::sqrt(4134.375) * x * x * y * y * y * y * y * z * z * z) + e_1 * (-std::sqrt(7816.552734375) * x * x * x * x * x * x * y * z + std::sqrt(34173.193359375) * x * x * x * x * y * y * y * z + std::sqrt(4134.375) * x * x * x * x * y * z * z * z - std::sqrt(70348.974609375) * x * x * y * y * y * y * y * z + std::sqrt(16537.5) * x * x * y * y * y * z * z * z - std::sqrt(581.396484375) * y * y * y * y * y * y * y * z + std::sqrt(4134.375) * y * y * y * y * y * z * z * z) + e_2 * (-std::sqrt(83721.09375) * x * x * x * x * y * z - std::sqrt(334884.375) * x * x * y * y * y * z + std::sqrt(148837.5) * x * x * y * z * z * z - std::sqrt(83721.09375) * y * y * y * y * y * z + std::sqrt(148837.5) * y * y * y * z * z * z) + e_3 * (-std::sqrt(595350.0) * x * x * y * z - std::sqrt(595350.0) * y * y * y * z + std::sqrt(264600.0) * y * z * z * z) + e_4 * (-std::sqrt(148837.5) * y * z);

        pc_10[k] = e_0 * (-std::sqrt(2.1533203125) * x * x * x * x * x * x * x * x * x * y + std::sqrt(77.51953125) * x * x * x * x * x * x * x * y * y * y + std::sqrt(215.33203125) * x * x * x * x * x * x * x * y * z * z - std::sqrt(10551.26953125) * x * x * x * x * x * y * y * y * z * z - std::sqrt(77.51953125) * x * x * x * y * y * y * y * y * y * y + std::sqrt(10551.26953125) * x * x * x * y * y * y * y * y * z * z + std::sqrt(2.1533203125) * x * y * y * y * y * y * y * y * y * y - std::sqrt(215.33203125) * x * y * y * y * y * y * y * y * z * z) + e_1 * (-std::sqrt(137.8125) * x * x * x * x * x * x * x * y + std::sqrt(6752.8125) * x * x * x * x * x * y * y * y - std::sqrt(6752.8125) * x * x * x * y * y * y * y * y + std::sqrt(137.8125) * x * y * y * y * y * y * y * y);
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

        pc_11[k] = e_0 * (std::sqrt(47.373046875) * x * x * x * x * x * x * x * x * y * z - std::sqrt(5732.138671875) * x * x * x * x * x * x * y * y * y * z + std::sqrt(10658.935546875) * x * x * x * x * y * y * y * y * y * z - std::sqrt(1184.326171875) * x * x * y * y * y * y * y * y * y * z) + e_1 * (-std::sqrt(1184.326171875) * x * x * x * x * x * x * y * z - std::sqrt(10658.935546875) * x * x * x * x * y * y * y * z - std::sqrt(10658.935546875) * x * x * y * y * y * y * y * z - std::sqrt(1184.326171875) * y * y * y * y * y * y * y * z) + e_2 * (-std::sqrt(170542.96875) * x * x * x * x * y * z - std::sqrt(682171.875) * x * x * y * y * y * z - std::sqrt(170542.96875) * y * y * y * y * y * z) + e_3 * (-std::sqrt(2728687.5) * x * x * y * z - std::sqrt(2728687.5) * y * y * y * z) + e_4 * (-std::sqrt(2728687.5) * y * z);

        pc_12[k] = e_0 * (std::sqrt(3.94775390625) * x * x * x * x * x * x * x * x * x * y - std::sqrt(1010.625) * x * x * x * x * x * x * x * y * y * y + std::sqrt(3552.978515625) * x * x * x * x * x * y * y * y * y * y - std::sqrt(1010.625) * x * x * x * y * y * y * y * y * y * y + std::sqrt(3.94775390625) * x * y * y * y * y * y * y * y * y * y) + e_1 * (-std::sqrt(568.4765625) * x * x * x * x * x * x * x * y - std::sqrt(5116.2890625) * x * x * x * x * x * y * y * y - std::sqrt(5116.2890625) * x * x * x * y * y * y * y * y - std::sqrt(568.4765625) * x * y * y * y * y * y * y * y) + e_2 * (-std::sqrt(127907.2265625) * x * x * x * x * x * y - std::sqrt(511628.90625) * x * x * x * y * y * y - std::sqrt(127907.2265625) * x * y * y * y * y * y) + e_3 * (-std::sqrt(3638250.0) * x * x * x * y - std::sqrt(3638250.0) * x * y * y * y) + e_4 * (-std::sqrt(8186062.5) * x * y);

        pc_13[k] = e_0 * (std::sqrt(639.5361328125) * x * x * x * x * x * x * x * y * y * z - std::sqrt(8598.2080078125) * x * x * x * x * x * y * y * y * y * z + std::sqrt(2850.2783203125) * x * x * x * y * y * y * y * y * y * z - std::sqrt(71.0595703125) * x * y * y * y * y * y * y * y * y * z) + e_1 * (std::sqrt(639.5361328125) * x * x * x * x * x * x * x * z - std::sqrt(639.5361328125) * x * x * x * x * x * y * y * z - std::sqrt(15988.4033203125) * x * x * x * y * y * y * y * z - std::sqrt(5755.8251953125) * x * y * y * y * y * y * y * z) + e_2 * (std::sqrt(63953.61328125) * x * x * x * x * x * z - std::sqrt(255814.453125) * x * x * x * y * y * z - std::sqrt(575582.51953125) * x * y * y * y * y * z) + e_3 * (std::sqrt(454781.25) * x * x * x * z - std::sqrt(4093031.25) * x * y * y * z);
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

        pc_14[k] = e_0 * (std::sqrt(5329.4677734375) * x * x * x * x * x * x * y * y * z * z - std::sqrt(29015.9912109375) * x * x * x * x * y * y * y * y * z * z + std::sqrt(4003.0224609375) * x * x * y * y * y * y * y * y * z * z - std::sqrt(23.6865234375) * y * y * y * y * y * y * y * y * z * z) + e_1 * (std::sqrt(5329.4677734375) * x * x * x * x * x * x * y * y + std::sqrt(5329.4677734375) * x * x * x * x * x * x * z * z - std::sqrt(29015.9912109375) * x * x * x * x * y * y * y * y + std::sqrt(5329.4677734375) * x * x * x * x * y * y * z * z + std::sqrt(4003.0224609375) * x * x * y * y * y * y * y * y - std::sqrt(5329.4677734375) * x * x * y * y * y * y * z * z - std::sqrt(23.6865234375) * y * y * y * y * y * y * y * y - std::sqrt(5329.4677734375) * y * y * y * y * y * y * z * z) + e_2 * (std::sqrt(5329.4677734375) * x * x * x * x * x * x + std::sqrt(5329.4677734375) * x * x * x * x * y * y + std::sqrt(341085.9375) * x * x * x * x * z * z - std::sqrt(5329.4677734375) * x * x * y * y * y * y - std::sqrt(5329.4677734375) * y * y * y * y * y * y - std::sqrt(341085.9375) * y * y * y * y * z * z) + e_3 * (std::sqrt(341085.9375) * x * x * x * x + std::sqrt(1364343.75) * x * x * z * z - std::sqrt(341085.9375) * y * y * y * y - std::sqrt(1364343.75) * y * y * z * z) + e_4 * (std::sqrt(1364343.75) * x * x - std::sqrt(1364343.75) * y * y);

        pc_15[k] = e_0 * (-std::sqrt(155.0390625) * x * x * x * x * x * x * x * y * y * z + std::sqrt(17.2265625) * x * x * x * x * x * y * y * y * y * z + std::sqrt(15503.90625) * x * x * x * x * x * y * y * z * z * z + std::sqrt(155.0390625) * x * x * x * y * y * y * y * y * y * z - std::sqrt(27562.5) * x * x * x * y * y * y * y * z * z * z - std::sqrt(17.2265625) * x * y * y * y * y * y * y * y * y * z + std::sqrt(1722.65625) * x * y * y * y * y * y * y * z * z * z) + e_1 * (-std::sqrt(155.0390625) * x * x * x * x * x * x * x * z + std::sqrt(18759.7265625) * x * x * x * x * x * y * y * z + std::sqrt(15503.90625) * x * x * x * x * x * z * z * z - std::sqrt(72782.2265625) * x * x * x * y * y * y * y * z + std::sqrt(62015.625) * x * x * x * y * y * z * z * z + std::sqrt(2084.4140625) * x * y * y * y * y * y * y * z + std::sqrt(15503.90625) * x * y * y * y * y * z * z * z) + e_2 * (std::sqrt(15503.90625) * x * x * x * x * x * z + std::sqrt(62015.625) * x * x * x * y * y * z + std::sqrt(558140.625) * x * x * x * z * z * z + std::sqrt(15503.90625) * x * y * y * y * y * z + std::sqrt(558140.625) * x * y * y * z * z * z) + e_3 * (std::sqrt(1550390.625) * x * x * x * z + std::sqrt(1550390.625) * x * y * y * z + std::sqrt(992250.0) * x * z * z * z) + e_4 * (std::sqrt(3969000.0) * x * z);
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

        pc_16[k] = e_0 * (-std::sqrt(2616.2841796875) * x * x * x * x * x * x * y * y * z * z - std::sqrt(290.6982421875) * x * x * x * x * y * y * y * y * z * z + std::sqrt(18604.6875) * x * x * x * x * y * y * z * z * z * z + std::sqrt(807.4951171875) * x * x * y * y * y * y * y * y * z * z - std::sqrt(8268.75) * x * x * y * y * y * y * z * z * z * z - std::sqrt(32.2998046875) * y * y * y * y * y * y * y * y * z * z + std::sqrt(229.6875) * y * y * y * y * y * y * z * z * z * z) + e_1 * (-std::sqrt(2616.2841796875) * x * x * x * x * x * x * y * y - std::sqrt(2616.2841796875) * x * x * x * x * x * x * z * z - std::sqrt(290.6982421875) * x * x * x * x * y * y * y * y - std::sqrt(2616.2841796875) * x * x * x * x * y * y * z * z + std::sqrt(18604.6875) * x * x * x * x * z * z * z * z + std::sqrt(807.4951171875) * x * x * y * y * y * y * y * y - std::sqrt(49128.0029296875) * x * x * y * y * y * y * z * z + std::sqrt(74418.75) * x * x * y * y * z * z * z * z - std::sqrt(32.2998046875) * y * y * y * y * y * y * y * y - std::sqrt(1582.6904296875) * y * y * y * y * y * y * z * z + std::sqrt(18604.6875) * y * y * y * y * z * z * z * z) + e_2 * (-std::sqrt(2616.2841796875) * x * x * x * x * x * x - std::sqrt(211919.0185546875) * x * x * x * x * y * y + std::sqrt(2616.2841796875) * x * x * y * y * y * y + std::sqrt(297675.0) * x * x * z * z * z * z - std::sqrt(7267.4560546875) * y * y * y * y * y * y + std::sqrt(297675.0) * y * y * z * z * z * z) + e_3 * (-std::sqrt(167442.1875) * x * x * x * x - std::sqrt(669768.75) * x * x * y * y + std::sqrt(1190700.0) * x * x * z * z - std::sqrt(167442.1875) * y * y * y * y + std::sqrt(1190700.0) * y * y * z * z + std::sqrt(132300.0) * z * z * z * z) + e_4 * (-std::sqrt(297675.0) * x * x - std::sqrt(297675.0) * y * y + std::sqrt(1190700.0) * z * z);
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

        pc_17[k] = e_0 * (std::sqrt(32.2998046875) * x * x * x * x * x * x * x * y * y * z + std::sqrt(89.7216796875) * x * x * x * x * x * y * y * y * y * z - std::sqrt(8268.75) * x * x * x * x * x * y * y * z * z * z + std::sqrt(3.5888671875) * x * x * x * y * y * y * y * y * y * z - std::sqrt(3675.0) * x * x * x * y * y * y * y * z * z * z + std::sqrt(8268.75) * x * x * x * y * y * z * z * z * z * z - std::sqrt(3.5888671875) * x * y * y * y * y * y * y * y * y * z + std::sqrt(918.75) * x * y * y * y * y * y * y * z * z * z - std::sqrt(918.75) * x * y * y * y * y * z * z * z * z * z) + e_1 * (std::sqrt(32.2998046875) * x * x * x * x * x * x * x * z - std::sqrt(9334.6435546875) * x * x * x * x * x * y * y * z - std::sqrt(8268.75) * x * x * x * x * x * z * z * z - std::sqrt(3448.9013671875) * x * x * x * y * y * y * y * z - std::sqrt(132300.0) * x * x * x * y * y * z * z * z + std::sqrt(8268.75) * x * x * x * z * z * z * z * z + std::sqrt(1898.5107421875) * x * y * y * y * y * y * y * z - std::sqrt(918.75) * x * y * y * y * y * z * z * z + std::sqrt(8268.75) * x * y * y * z * z * z * z * z) + e_2 * (-std::sqrt(15633.10546875) * x * x * x * x * x * z - std::sqrt(1451682.421875) * x * x * x * y * y * z - std::sqrt(33075.0) * x * x * x * z * z * z + std::sqrt(37338.57421875) * x * y * y * y * y * z - std::sqrt(33075.0) * x * y * y * z * z * z + std::sqrt(33075.0) * x * z * z * z * z * z) + e_3 * (-std::sqrt(1000518.75) * x * x * x * z - std::sqrt(1000518.75) * x * y * y * z + std::sqrt(132300.0) * x * z * z * z) + e_4 * (-std::sqrt(529200.0) * x * z);
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

        pc_18[k] = e_0 * (std::sqrt(322.998046875) * x * x * x * x * x * x * y * y * z * z + std::sqrt(897.216796875) * x * x * x * x * y * y * y * y * z * z - std::sqrt(5167.96875) * x * x * x * x * y * y * z * z * z * z + std::sqrt(35.888671875) * x * x * y * y * y * y * y * y * z * z - std::sqrt(2296.875) * x * x * y * y * y * y * z * z * z * z + std::sqrt(826.875) * x * x * y * y * z * z * z * z * z * z - std::sqrt(35.888671875) * y * y * y * y * y * y * y * y * z * z + std::sqrt(574.21875) * y * y * y * y * y * y * z * z * z * z - std::sqrt(91.875) * y * y * y * y * z * z * z * z * z * z) + e_1 * (std::sqrt(322.998046875) * x * x * x * x * x * x * y * y + std::sqrt(322.998046875) * x * x * x * x * x * x * z * z + std::sqrt(897.216796875) * x * x * x * x * y * y * y * y + std::sqrt(322.998046875) * x * x * x * x * y * y * z * z - std::sqrt(5167.96875) * x * x * x * x * z * z * z * z + std::sqrt(35.888671875) * x * x * y * y * y * y * y * y - std::sqrt(322.998046875) * x * x * y * y * y * y * z * z - std::sqrt(82687.5) * x * x * y * y * z * z * z * z + std::sqrt(826.875) * x * x * z * z * z * z * z * z - std::sqrt(35.888671875) * y * y * y * y * y * y * y * y - std::sqrt(322.998046875) * y * y * y * y * y * y * z * z + std::sqrt(28136.71875) * y * y * y * y * z * z * z * z - std::sqrt(826.875) * y * y * z * z * z * z * z * z) + e_2 * (std::sqrt(322.998046875) * x * x * x * x * x * x + std::sqrt(54586.669921875) * x * x * x * x * y * y - std::sqrt(5167.96875) * x * x * x * x * z * z + std::sqrt(15826.904296875) * x * x * y * y * y * y - std::sqrt(744187.5) * x * x * y * y * z * z - std::sqrt(20671.875) * x * x * z * z * z * z - std::sqrt(8074.951171875) * y * y * y * y * y * y + std::sqrt(129199.21875) * y * y * y * y * z * z + std::sqrt(20671.875) * y * y * z * z * z * z) + e_3 * (std::sqrt(20671.875) * x * x * x * x + std::sqrt(186046.875) * x * x * y * y - std::sqrt(516796.875) * x * x * z * z - std::sqrt(82687.5) * y * y * y * y + std::sqrt(516796.875) * y * y * z * z) + e_4 * (std::sqrt(20671.875) * x * x - std::sqrt(20671.875) * y * y);
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

        pc_19[k] = e_0 * (-std::sqrt(3.84521484375) * x * x * x * x * x * x * x * x * y * z - std::sqrt(27.34375) * x * x * x * x * x * x * y * y * y * z + std::sqrt(1245.849609375) * x * x * x * x * x * x * y * z * z * z - std::sqrt(15.380859375) * x * x * x * x * y * y * y * y * y * z + std::sqrt(3460.693359375) * x * x * x * x * y * y * y * z * z * z - std::sqrt(2214.84375) * x * x * x * x * y * z * z * z * z * z + std::sqrt(138.427734375) * x * x * y * y * y * y * y * z * z * z - std::sqrt(984.375) * x * x * y * y * y * z * z * z * z * z + std::sqrt(39.375) * x * x * y * z * z * z * z * z * z * z + std::sqrt(0.42724609375) * y * y * y * y * y * y * y * y * y * z - std::sqrt(138.427734375) * y * y * y * y * y * y * y * z * z * z + std::sqrt(246.09375) * y * y * y * y * y * z * z * z * z * z - std::sqrt(4.375) * y * y * y * z * z * z * z * z * z * z) + e_1 * (std::sqrt(1245.849609375) * x * x * x * x * x * x * y * z + std::sqrt(3460.693359375) * x * x * x * x * y * y * y * z + std::sqrt(55371.09375) * x * x * x * x * y * z * z * z + std::sqrt(138.427734375) * x * x * y * y * y * y * y * z + std::sqrt(24609.375) * x * x * y * y * y * z * z * z - std::sqrt(59889.375) * x * x * y * z * z * z * z * z - std::sqrt(138.427734375) * y * y * y * y * y * y * y * z - std::sqrt(6152.34375) * y * y * y * y * y * z * z * z + std::sqrt(6654.375) * y * y * y * z * z * z * z * z) + e_2 * (std::sqrt(498339.84375) * x * x * x * x * y * z + std::sqrt(221484.375) * x * x * y * y * y * z - std::sqrt(79734.375) * x * x * y * z * z * z - std::sqrt(55371.09375) * y * y * y * y * y * z + std::sqrt(8859.375) * y * y * y * z * z * z) + e_3 * (std::sqrt(2560359.375) * x * x * y * z - std::sqrt(284484.375) * y * y * y * z);
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

        pc_20[k] = e_0 * (std::sqrt(322.998046875) * x * x * x * x * x * x * x * y * z * z + std::sqrt(897.216796875) * x * x * x * x * x * y * y * y * z * z - std::sqrt(5167.96875) * x * x * x * x * x * y * z * z * z * z + std::sqrt(35.888671875) * x * x * x * y * y * y * y * y * z * z - std::sqrt(2296.875) * x * x * x * y * y * y * z * z * z * z + std::sqrt(826.875) * x * x * x * y * z * z * z * z * z * z - std::sqrt(35.888671875) * x * y * y * y * y * y * y * y * z * z + std::sqrt(574.21875) * x * y * y * y * y * y * z * z * z * z - std::sqrt(91.875) * x * y * y * y * z * z * z * z * z * z) + e_1 * (std::sqrt(322.998046875) * x * x * x * x * x * x * x * y + std::sqrt(897.216796875) * x * x * x * x * x * y * y * y + std::sqrt(1291.9921875) * x * x * x * x * x * y * z * z + std::sqrt(35.888671875) * x * x * x * y * y * y * y * y + std::sqrt(5167.96875) * x * x * x * y * y * y * z * z - std::sqrt(186046.875) * x * x * x * y * z * z * z * z - std::sqrt(35.888671875) * x * y * y * y * y * y * y * y + std::sqrt(1291.9921875) * x * y * y * y * y * y * z * z - std::sqrt(2296.875) * x * y * y * y * z * z * z * z + std::sqrt(3307.5) * x * y * z * z * z * z * z * z) + e_2 * (std::sqrt(63307.6171875) * x * x * x * x * x * y + std::sqrt(46511.71875) * x * x * x * y * y * y - std::sqrt(1012921.875) * x * x * x * y * z * z - std::sqrt(1291.9921875) * x * y * y * y * y * y + std::sqrt(20671.875) * x * y * y * y * z * z - std::sqrt(82687.5) * x * y * z * z * z * z) + e_3 * (std::sqrt(516796.875) * x * x * x * y + std::sqrt(20671.875) * x * y * y * y - std::sqrt(2067187.5) * x * y * z * z) + e_4 * (std::sqrt(82687.5) * x * y);
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

        pc_21[k] = e_0 * (std::sqrt(8.074951171875) * x * x * x * x * x * x * x * x * y * z + std::sqrt(3.5888671875) * x * x * x * x * x * x * y * y * y * z - std::sqrt(2067.1875) * x * x * x * x * x * x * y * z * z * z - std::sqrt(14.35546875) * x * x * x * x * y * y * y * y * y * z + std::sqrt(229.6875) * x * x * x * x * y * y * y * z * z * z + std::sqrt(2067.1875) * x * x * x * x * y * z * z * z * z * z - std::sqrt(3.5888671875) * x * x * y * y * y * y * y * y * y * z + std::sqrt(2067.1875) * x * x * y * y * y * y * y * z * z * z - std::sqrt(3675.0) * x * x * y * y * y * z * z * z * z * z + std::sqrt(0.897216796875) * y * y * y * y * y * y * y * y * y * z - std::sqrt(229.6875) * y * y * y * y * y * y * y * z * z * z + std::sqrt(229.6875) * y * y * y * y * y * z * z * z * z * z) + e_1 * (-std::sqrt(2616.2841796875) * x * x * x * x * x * x * y * z + std::sqrt(1295.5810546875) * x * x * x * x * y * y * y * z - std::sqrt(33075.0) * x * x * x * x * y * z * z * z + std::sqrt(5458.6669921875) * x * x * y * y * y * y * y * z - std::sqrt(3675.0) * x * x * y * y * y * z * z * z + std::sqrt(8268.75) * x * x * y * z * z * z * z * z - std::sqrt(175.8544921875) * y * y * y * y * y * y * y * z - std::sqrt(14700.0) * y * y * y * y * y * z * z * z + std::sqrt(8268.75) * y * y * y * z * z * z * z * z) + e_2 * (-std::sqrt(362920.60546875) * x * x * x * x * y * z + std::sqrt(149354.296875) * x * x * y * y * y * z - std::sqrt(33075.0) * x * x * y * z * z * z - std::sqrt(80749.51171875) * y * y * y * y * y * z - std::sqrt(33075.0) * y * y * y * z * z * z + std::sqrt(33075.0) * y * z * z * z * z * z) + e_3 * (-std::sqrt(1000518.75) * x * x * y * z - std::sqrt(1000518.75) * y * y * y * z + std::sqrt(132300.0) * y * z * z * z) + e_4 * (-std::sqrt(529200.0) * y * z);

        pc_22[k] = e_0 * (-std::sqrt(290.6982421875) * x * x * x * x * x * x * x * y * z * z + std::sqrt(1582.6904296875) * x * x * x * x * x * y * y * y * z * z + std::sqrt(2067.1875) * x * x * x * x * x * y * z * z * z * z + std::sqrt(1582.6904296875) * x * x * x * y * y * y * y * y * z * z - std::sqrt(22968.75) * x * x * x * y * y * y * z * z * z * z - std::sqrt(290.6982421875) * x * y * y * y * y * y * y * y * z * z + std::sqrt(2067.1875) * x * y * y * y * y * y * z * z * z * z) + e_1 * (-std::sqrt(290.6982421875) * x * x * x * x * x * x * x * y + std::sqrt(1582.6904296875) * x * x * x * x * x * y * y * y + std::sqrt(1162.79296875) * x * x * x * x * x * y * z * z + std::sqrt(1582.6904296875) * x * x * x * y * y * y * y * y - std::sqrt(12919.921875) * x * x * x * y * y * y * z * z - std::sqrt(290.6982421875) * x * y * y * y * y * y * y * y + std::sqrt(1162.79296875) * x * y * y * y * y * y * z * z) + e_2 * (-std::sqrt(10465.13671875) * x * x * x * x * x * y + std::sqrt(116279.296875) * x * x * x * y * y * y - std::sqrt(10465.13671875) * x * y * y * y * y * y);
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

        pc_23[k] = e_0 * (-std::sqrt(9.68994140625) * x * x * x * x * x * x * x * x * y * z + std::sqrt(275.625) * x * x * x * x * x * x * y * y * y * z + std::sqrt(968.994140625) * x * x * x * x * x * x * y * z * z * z + std::sqrt(107.666015625) * x * x * x * x * y * y * y * y * y * z - std::sqrt(38867.431640625) * x * x * x * x * y * y * y * z * z * z - std::sqrt(68.90625) * x * x * y * y * y * y * y * y * y * z + std::sqrt(8720.947265625) * x * x * y * y * y * y * y * z * z * z + std::sqrt(1.07666015625) * y * y * y * y * y * y * y * y * y * z - std::sqrt(107.666015625) * y * y * y * y * y * y * y * z * z * z) + e_1 * (std::sqrt(3139.541015625) * x * x * x * x * x * x * y * z - std::sqrt(56955.322265625) * x * x * x * x * y * y * y * z - std::sqrt(15503.90625) * x * x * x * x * y * z * z * z + std::sqrt(28255.869140625) * x * x * y * y * y * y * y * z - std::sqrt(62015.625) * x * x * y * y * y * z * z * z - std::sqrt(4.306640625) * y * y * y * y * y * y * y * z - std::sqrt(15503.90625) * y * y * y * y * y * z * z * z) + e_2 * (-std::sqrt(15503.90625) * x * x * x * x * y * z - std::sqrt(62015.625) * x * x * y * y * y * z - std::sqrt(558140.625) * x * x * y * z * z * z - std::sqrt(15503.90625) * y * y * y * y * y * z - std::sqrt(558140.625) * y * y * y * z * z * z) + e_3 * (-std::sqrt(1550390.625) * x * x * y * z - std::sqrt(1550390.625) * y * y * y * z - std::sqrt(992250.0) * y * z * z * z) + e_4 * (-std::sqrt(3969000.0) * y * z);

        pc_24[k] = e_0 * (std::sqrt(213.1787109375) * x * x * x * x * x * x * x * y * z * z - std::sqrt(22762.7490234375) * x * x * x * x * x * y * y * y * z * z + std::sqrt(14804.0771484375) * x * x * x * y * y * y * y * y * z * z - std::sqrt(592.1630859375) * x * y * y * y * y * y * y * y * z * z) + e_1 * (std::sqrt(213.1787109375) * x * x * x * x * x * x * x * y - std::sqrt(22762.7490234375) * x * x * x * x * x * y * y * y - std::sqrt(21317.87109375) * x * x * x * x * x * y * z * z + std::sqrt(14804.0771484375) * x * x * x * y * y * y * y * y - std::sqrt(85271.484375) * x * x * x * y * y * y * z * z - std::sqrt(592.1630859375) * x * y * y * y * y * y * y * y - std::sqrt(21317.87109375) * x * y * y * y * y * y * z * z) + e_2 * (-std::sqrt(21317.87109375) * x * x * x * x * x * y - std::sqrt(85271.484375) * x * x * x * y * y * y - std::sqrt(1364343.75) * x * x * x * y * z * z - std::sqrt(21317.87109375) * x * y * y * y * y * y - std::sqrt(1364343.75) * x * y * y * y * z * z) + e_3 * (-std::sqrt(1364343.75) * x * x * x * y - std::sqrt(1364343.75) * x * y * y * y - std::sqrt(5457375.0) * x * y * z * z) + e_4 * (-std::sqrt(5457375.0) * x * y);
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

        pc_25[k] = e_0 * (std::sqrt(17.764892578125) * x * x * x * x * x * x * x * x * y * z - std::sqrt(4176.7236328125) * x * x * x * x * x * x * y * y * y * z + std::sqrt(7105.95703125) * x * x * x * x * y * y * y * y * y * z - std::sqrt(639.5361328125) * x * x * y * y * y * y * y * y * y * z + std::sqrt(1.973876953125) * y * y * y * y * y * y * y * y * y * z) + e_1 * (-std::sqrt(5755.8251953125) * x * x * x * x * x * x * y * z - std::sqrt(15988.4033203125) * x * x * x * x * y * y * y * z - std::sqrt(639.5361328125) * x * x * y * y * y * y * y * z + std::sqrt(639.5361328125) * y * y * y * y * y * y * y * z) + e_2 * (-std::sqrt(575582.51953125) * x * x * x * x * y * z - std::sqrt(255814.453125) * x * x * y * y * y * z + std::sqrt(63953.61328125) * y * y * y * y * y * z) + e_3 * (-std::sqrt(4093031.25) * x * x * y * z + std::sqrt(454781.25) * y * y * y * z);

        pc_26[k] = e_0 * (-std::sqrt(20.302734375) * x * x * x * x * x * x * x * x * y * y + std::sqrt(110.537109375) * x * x * x * x * x * x * y * y * y * y + std::sqrt(730.8984375) * x * x * x * x * x * x * y * y * z * z + std::sqrt(110.537109375) * x * x * x * x * y * y * y * y * y * y - std::sqrt(8121.09375) * x * x * x * x * y * y * y * y * z * z - std::sqrt(20.302734375) * x * x * y * y * y * y * y * y * y * y + std::sqrt(730.8984375) * x * x * y * y * y * y * y * y * z * z) + e_1 * (-std::sqrt(20.302734375) * x * x * x * x * x * x * x * x - std::sqrt(1299.375) * x * x * x * x * x * x * y * y + std::sqrt(730.8984375) * x * x * x * x * x * x * z * z + std::sqrt(50756.8359375) * x * x * x * x * y * y * y * y - std::sqrt(18272.4609375) * x * x * x * x * y * y * z * z - std::sqrt(1299.375) * x * x * y * y * y * y * y * y - std::sqrt(18272.4609375) * x * x * y * y * y * y * z * z - std::sqrt(20.302734375) * y * y * y * y * y * y * y * y + std::sqrt(730.8984375) * y * y * y * y * y * y * z * z) + e_2 * (-std::sqrt(4568.115234375) * x * x * x * x * x * x + std::sqrt(114202.880859375) * x * x * x * x * y * y + std::sqrt(18272.4609375) * x * x * x * x * z * z + std::sqrt(114202.880859375) * x * x * y * y * y * y - std::sqrt(657808.59375) * x * x * y * y * z * z - std::sqrt(4568.115234375) * y * y * y * y * y * y + std::sqrt(18272.4609375) * y * y * y * y * z * z) + e_3 * (-std::sqrt(32484.375) * x * x * x * x + std::sqrt(1169437.5) * x * x * y * y - std::sqrt(32484.375) * y * y * y * y);
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

        pc_27[k] = e_0 * (-std::sqrt(169.189453125) * x * x * x * x * x * x * x * y * y * z + std::sqrt(169.189453125) * x * x * x * x * x * y * y * y * y * z + std::sqrt(6090.8203125) * x * x * x * x * x * y * y * z * z * z + std::sqrt(548.173828125) * x * x * x * y * y * y * y * y * y * z - std::sqrt(24363.28125) * x * x * x * y * y * y * y * z * z * z - std::sqrt(6.767578125) * x * y * y * y * y * y * y * y * y * z + std::sqrt(243.6328125) * x * y * y * y * y * y * y * z * z * z) + e_1 * (-std::sqrt(169.189453125) * x * x * x * x * x * x * x * z + std::sqrt(1522.705078125) * x * x * x * x * x * y * y * z + std::sqrt(6090.8203125) * x * x * x * x * x * z * z * z + std::sqrt(169.189453125) * x * x * x * y * y * y * y * z - std::sqrt(24363.28125) * x * x * x * y * y * z * z * z + std::sqrt(1955.830078125) * x * y * y * y * y * y * y * z - std::sqrt(54817.3828125) * x * y * y * y * y * z * z * z) + e_2 * (std::sqrt(97453.125) * x * x * x * z * z * z - std::sqrt(877078.125) * x * y * y * z * z * z) + e_3 * (std::sqrt(97453.125) * x * x * x * z - std::sqrt(877078.125) * x * y * y * z);
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

        pc_28[k] = e_0 * (std::sqrt(4.921875) * x * x * x * x * x * x * x * x * y * y + std::sqrt(4.921875) * x * x * x * x * x * x * y * y * y * y - std::sqrt(1260.0) * x * x * x * x * x * x * y * y * z * z - std::sqrt(4.921875) * x * x * x * x * y * y * y * y * y * y + std::sqrt(17718.75) * x * x * x * x * y * y * z * z * z * z - std::sqrt(4.921875) * x * x * y * y * y * y * y * y * y * y + std::sqrt(1260.0) * x * x * y * y * y * y * y * y * z * z - std::sqrt(17718.75) * x * x * y * y * y * y * z * z * z * z) + e_1 * (std::sqrt(4.921875) * x * x * x * x * x * x * x * x + std::sqrt(1594.6875) * x * x * x * x * x * x * y * y - std::sqrt(1260.0) * x * x * x * x * x * x * z * z + std::sqrt(70875.0) * x * x * x * x * y * y * z * z + std::sqrt(17718.75) * x * x * x * x * z * z * z * z - std::sqrt(1594.6875) * x * x * y * y * y * y * y * y - std::sqrt(70875.0) * x * x * y * y * y * y * z * z - std::sqrt(4.921875) * y * y * y * y * y * y * y * y + std::sqrt(1260.0) * y * y * y * y * y * y * z * z - std::sqrt(17718.75) * y * y * y * y * z * z * z * z) + e_2 * (std::sqrt(1107.421875) * x * x * x * x * x * x + std::sqrt(187154.296875) * x * x * x * x * y * y + std::sqrt(70875.0) * x * x * x * x * z * z - std::sqrt(187154.296875) * x * x * y * y * y * y + std::sqrt(159468.75) * x * x * z * z * z * z - std::sqrt(1107.421875) * y * y * y * y * y * y - std::sqrt(70875.0) * y * y * y * y * z * z - std::sqrt(159468.75) * y * y * z * z * z * z) + e_3 * (std::sqrt(159468.75) * x * x * x * x + std::sqrt(1771875.0) * x * x * z * z - std::sqrt(159468.75) * y * y * y * y - std::sqrt(1771875.0) * y * y * z * z) + e_4 * (std::sqrt(868218.75) * x * x - std::sqrt(868218.75) * y * y);
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

        pc_29[k] = e_0 * (std::sqrt(83.056640625) * x * x * x * x * x * x * x * y * y * z + std::sqrt(230.712890625) * x * x * x * x * x * y * y * y * y * z - std::sqrt(6238.4765625) * x * x * x * x * x * y * y * z * z * z + std::sqrt(9.228515625) * x * x * x * y * y * y * y * y * y * z - std::sqrt(2772.65625) * x * x * x * y * y * y * y * z * z * z + std::sqrt(21262.5) * x * x * x * y * y * z * z * z * z * z - std::sqrt(9.228515625) * x * y * y * y * y * y * y * y * y * z + std::sqrt(693.1640625) * x * y * y * y * y * y * y * z * z * z - std::sqrt(2362.5) * x * y * y * y * y * z * z * z * z * z) + e_1 * (std::sqrt(83.056640625) * x * x * x * x * x * x * x * z + std::sqrt(2076.416015625) * x * x * x * x * x * y * y * z - std::sqrt(6238.4765625) * x * x * x * x * x * z * z * z + std::sqrt(1559.619140625) * x * x * x * y * y * y * y * z + std::sqrt(124178.90625) * x * x * x * y * y * z * z * z + std::sqrt(21262.5) * x * x * x * z * z * z * z * z + std::sqrt(9.228515625) * x * y * y * y * y * y * y * z - std::sqrt(62052.5390625) * x * y * y * y * y * z * z * z + std::sqrt(21262.5) * x * y * y * z * z * z * z * z) + e_2 * (std::sqrt(765450.0) * x * x * x * y * y * z + std::sqrt(260465.625) * x * x * x * z * z * z - std::sqrt(85050.0) * x * y * y * y * y * z + std::sqrt(260465.625) * x * y * y * z * z * z + std::sqrt(85050.0) * x * z * z * z * z * z) + e_3 * (std::sqrt(643190.625) * x * x * x * z + std::sqrt(643190.625) * x * y * y * z + std::sqrt(2731050.0) * x * z * z * z) + e_4 * (std::sqrt(4167450.0) * x * z);
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

        pc_30[k] = e_0 * (-std::sqrt(1.025390625) * x * x * x * x * x * x * x * x * y * y - std::sqrt(9.228515625) * x * x * x * x * x * x * y * y * y * y + std::sqrt(496.2890625) * x * x * x * x * x * x * y * y * z * z - std::sqrt(9.228515625) * x * x * x * x * y * y * y * y * y * y + std::sqrt(1985.15625) * x * x * x * x * y * y * y * y * z * z - std::sqrt(12862.5) * x * x * x * x * y * y * z * z * z * z - std::sqrt(1.025390625) * x * x * y * y * y * y * y * y * y * y + std::sqrt(496.2890625) * x * x * y * y * y * y * y * y * z * z - std::sqrt(12862.5) * x * x * y * y * y * y * z * z * z * z + std::sqrt(9450.0) * x * x * y * y * z * z * z * z * z * z) + e_1 * (-std::sqrt(1.025390625) * x * x * x * x * x * x * x * x - std::sqrt(590.625) * x * x * x * x * x * x * y * y + std::sqrt(496.2890625) * x * x * x * x * x * x * z * z - std::sqrt(2169.7265625) * x * x * x * x * y * y * y * y - std::sqrt(6238.4765625) * x * x * x * x * y * y * z * z - std::sqrt(12862.5) * x * x * x * x * z * z * z * z - std::sqrt(590.625) * x * x * y * y * y * y * y * y - std::sqrt(6238.4765625) * x * x * y * y * y * y * z * z + std::sqrt(9450.0) * x * x * y * y * z * z * z * z + std::sqrt(9450.0) * x * x * z * z * z * z * z * z - std::sqrt(1.025390625) * y * y * y * y * y * y * y * y + std::sqrt(496.2890625) * y * y * y * y * y * y * z * z - std::sqrt(12862.5) * y * y * y * y * z * z * z * z + std::sqrt(9450.0) * y * y * z * z * z * z * z * z) + e_2 * (-std::sqrt(230.712890625) * x * x * x * x * x * x - std::sqrt(130685.009765625) * x * x * x * x * y * y - std::sqrt(45219.7265625) * x * x * x * x * z * z - std::sqrt(130685.009765625) * x * x * y * y * y * y - std::sqrt(33222.65625) * x * x * y * y * z * z + std::sqrt(191362.5) * x * x * z * z * z * z - std::sqrt(230.712890625) * y * y * y * y * y * y - std::sqrt(45219.7265625) * y * y * y * y * z * z + std::sqrt(191362.5) * y * y * z * z * z * z + std::sqrt(9450.0) * z * z * z * z * z * z) + e_3 * (-std::sqrt(71465.625) * x * x * x * x - std::sqrt(2270362.5) * x * x * y * y + std::sqrt(151200.0) * x * x * z * z - std::sqrt(71465.625) * y * y * y * y + std::sqrt(151200.0) * y * y * z * z + std::sqrt(604800.0) * z * z * z * z) + e_4 * (-std::sqrt(463050.0) * x * x - std::sqrt(463050.0) * y * y + std::sqrt(1852200.0) * z * z);
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

        pc_31[k] = e_0 * (-std::sqrt(10.25390625) * x * x * x * x * x * x * x * y * y * z - std::sqrt(92.28515625) * x * x * x * x * x * y * y * y * y * z + std::sqrt(1025.390625) * x * x * x * x * x * y * y * z * z * z - std::sqrt(92.28515625) * x * x * x * y * y * y * y * y * y * z + std::sqrt(4101.5625) * x * x * x * y * y * y * y * z * z * z - std::sqrt(6720.0) * x * x * x * y * y * z * z * z * z * z - std::sqrt(10.25390625) * x * y * y * y * y * y * y * y * y * z + std::sqrt(1025.390625) * x * y * y * y * y * y * y * z * z * z - std::sqrt(6720.0) * x * y * y * y * y * z * z * z * z * z + std::sqrt(945.0) * x * y * y * z * z * z * z * z * z * z) + e_1 * (-std::sqrt(10.25390625) * x * x * x * x * x * x * x * z - std::sqrt(830.56640625) * x * x * x * x * x * y * y * z + std::sqrt(1025.390625) * x * x * x * x * x * z * z * z - std::sqrt(2307.12890625) * x * x * x * y * y * y * y * z - std::sqrt(13289.0625) * x * x * x * y * y * z * z * z - std::sqrt(6720.0) * x * x * x * z * z * z * z * z - std::sqrt(502.44140625) * x * y * y * y * y * y * y * z - std::sqrt(21697.265625) * x * y * y * y * y * z * z * z - std::sqrt(8505.0) * x * y * y * z * z * z * z * z + std::sqrt(945.0) * x * z * z * z * z * z * z * z) + e_2 * (-std::sqrt(212625.0) * x * x * x * y * y * z - std::sqrt(94500.0) * x * x * x * z * z * z - std::sqrt(212625.0) * x * y * y * y * y * z - std::sqrt(1157625.0) * x * y * y * z * z * z + std::sqrt(23625.0) * x * z * z * z * z * z) + e_3 * (-std::sqrt(212625.0) * x * x * x * z - std::sqrt(6048000.0) * x * y * y * z - std::sqrt(23625.0) * x * z * z * z) + e_4 * (-std::sqrt(1157625.0) * x * z);
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

        pc_32[k] = e_0 * (std::sqrt(0.1220703125) * x * x * x * x * x * x * x * x * x * y + std::sqrt(1.953125) * x * x * x * x * x * x * x * y * y * y - std::sqrt(70.3125) * x * x * x * x * x * x * x * y * z * z + std::sqrt(4.39453125) * x * x * x * x * x * y * y * y * y * y - std::sqrt(632.8125) * x * x * x * x * x * y * y * y * z * z + std::sqrt(2126.953125) * x * x * x * x * x * y * z * z * z * z + std::sqrt(1.953125) * x * x * x * y * y * y * y * y * y * y - std::sqrt(632.8125) * x * x * x * y * y * y * y * y * z * z + std::sqrt(8507.8125) * x * x * x * y * y * y * z * z * z * z - std::sqrt(2645.0) * x * x * x * y * z * z * z * z * z * z + std::sqrt(0.1220703125) * x * y * y * y * y * y * y * y * y * y - std::sqrt(70.3125) * x * y * y * y * y * y * y * y * z * z + std::sqrt(2126.953125) * x * y * y * y * y * y * z * z * z * z - std::sqrt(2645.0) * x * y * y * y * z * z * z * z * z * z + std::sqrt(45.0) * x * y * z * z * z * z * z * z * z * z) + e_1 * (std::sqrt(70.3125) * x * x * x * x * x * x * x * y + std::sqrt(632.8125) * x * x * x * x * x * y * y * y + std::sqrt(632.8125) * x * x * x * x * x * y * z * z + std::sqrt(632.8125) * x * x * x * y * y * y * y * y + std::sqrt(2531.25) * x * x * x * y * y * y * z * z - std::sqrt(1125.0) * x * x * x * y * z * z * z * z + std::sqrt(70.3125) * x * y * y * y * y * y * y * y + std::sqrt(632.8125) * x * y * y * y * y * y * z * z - std::sqrt(1125.0) * x * y * y * y * z * z * z * z - std::sqrt(14580.0) * x * y * z * z * z * z * z * z) + e_2 * (std::sqrt(19142.578125) * x * x * x * x * x * y + std::sqrt(76570.3125) * x * x * x * y * y * y + std::sqrt(10125.0) * x * x * x * y * z * z + std::sqrt(19142.578125) * x * y * y * y * y * y + std::sqrt(10125.0) * x * y * y * y * z * z - std::sqrt(1012500.0) * x * y * z * z * z * z) + e_3 * (std::sqrt(595125.0) * x * x * x * y + std::sqrt(595125.0) * x * y * y * y - std::sqrt(3280500.0) * x * y * z * z) + e_4 * (std::sqrt(496125.0) * x * y);
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

        pc_33[k] = e_0 * (-std::sqrt(10.25390625) * x * x * x * x * x * x * x * x * y * z - std::sqrt(92.28515625) * x * x * x * x * x * x * y * y * y * z + std::sqrt(1025.390625) * x * x * x * x * x * x * y * z * z * z - std::sqrt(92.28515625) * x * x * x * x * y * y * y * y * y * z + std::sqrt(4101.5625) * x * x * x * x * y * y * y * z * z * z - std::sqrt(6720.0) * x * x * x * x * y * z * z * z * z * z - std::sqrt(10.25390625) * x * x * y * y * y * y * y * y * y * z + std::sqrt(1025.390625) * x * x * y * y * y * y * y * z * z * z - std::sqrt(6720.0) * x * x * y * y * y * z * z * z * z * z + std::sqrt(945.0) * x * x * y * z * z * z * z * z * z * z) + e_1 * (-std::sqrt(502.44140625) * x * x * x * x * x * x * y * z - std::sqrt(2307.12890625) * x * x * x * x * y * y * y * z - std::sqrt(21697.265625) * x * x * x * x * y * z * z * z - std::sqrt(830.56640625) * x * x * y * y * y * y * y * z - std::sqrt(13289.0625) * x * x * y * y * y * z * z * z - std::sqrt(8505.0) * x * x * y * z * z * z * z * z - std::sqrt(10.25390625) * y * y * y * y * y * y * y * z + std::sqrt(1025.390625) * y * y * y * y * y * z * z * z - std::sqrt(6720.0) * y * y * y * z * z * z * z * z + std::sqrt(945.0) * y * z * z * z * z * z * z * z) + e_2 * (-std::sqrt(212625.0) * x * x * x * x * y * z - std::sqrt(212625.0) * x * x * y * y * y * z - std::sqrt(1157625.0) * x * x * y * z * z * z - std::sqrt(94500.0) * y * y * y * z * z * z + std::sqrt(23625.0) * y * z * z * z * z * z) + e_3 * (-std::sqrt(6048000.0) * x * x * y * z - std::sqrt(212625.0) * y * y * y * z - std::sqrt(23625.0) * y * z * z * z) + e_4 * (-std::sqrt(1157625.0) * y * z);
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

        pc_34[k] = e_0 * (-std::sqrt(0.25634765625) * x * x * x * x * x * x * x * x * x * y - std::sqrt(1.025390625) * x * x * x * x * x * x * x * y * y * y + std::sqrt(124.072265625) * x * x * x * x * x * x * x * y * z * z + std::sqrt(124.072265625) * x * x * x * x * x * y * y * y * z * z - std::sqrt(3215.625) * x * x * x * x * x * y * z * z * z * z + std::sqrt(1.025390625) * x * x * x * y * y * y * y * y * y * y - std::sqrt(124.072265625) * x * x * x * y * y * y * y * y * z * z + std::sqrt(2362.5) * x * x * x * y * z * z * z * z * z * z + std::sqrt(0.25634765625) * x * y * y * y * y * y * y * y * y * y - std::sqrt(124.072265625) * x * y * y * y * y * y * y * y * z * z + std::sqrt(3215.625) * x * y * y * y * y * y * z * z * z * z - std::sqrt(2362.5) * x * y * y * y * z * z * z * z * z * z) + e_1 * (-std::sqrt(102.5390625) * x * x * x * x * x * x * x * y - std::sqrt(102.5390625) * x * x * x * x * x * y * y * y - std::sqrt(5315.625) * x * x * x * x * x * y * z * z + std::sqrt(102.5390625) * x * x * x * y * y * y * y * y + std::sqrt(26250.0) * x * x * x * y * z * z * z * z + std::sqrt(102.5390625) * x * y * y * y * y * y * y * y + std::sqrt(5315.625) * x * y * y * y * y * y * z * z - std::sqrt(26250.0) * x * y * y * y * z * z * z * z) + e_2 * (-std::sqrt(24953.90625) * x * x * x * x * x * y + std::sqrt(14765.625) * x * x * x * y * z * z + std::sqrt(24953.90625) * x * y * y * y * y * y - std::sqrt(14765.625) * x * y * y * y * z * z) + e_3 * (-std::sqrt(236250.0) * x * x * x * y + std::sqrt(236250.0) * x * y * y * y);
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

        pc_35[k] = e_0 * (std::sqrt(9.228515625) * x * x * x * x * x * x * x * x * y * z - std::sqrt(9.228515625) * x * x * x * x * x * x * y * y * y * z - std::sqrt(693.1640625) * x * x * x * x * x * x * y * z * z * z - std::sqrt(230.712890625) * x * x * x * x * y * y * y * y * y * z + std::sqrt(2772.65625) * x * x * x * x * y * y * y * z * z * z + std::sqrt(2362.5) * x * x * x * x * y * z * z * z * z * z - std::sqrt(83.056640625) * x * x * y * y * y * y * y * y * y * z + std::sqrt(6238.4765625) * x * x * y * y * y * y * y * z * z * z - std::sqrt(21262.5) * x * x * y * y * y * z * z * z * z * z) + e_1 * (-std::sqrt(9.228515625) * x * x * x * x * x * x * y * z - std::sqrt(1559.619140625) * x * x * x * x * y * y * y * z + std::sqrt(62052.5390625) * x * x * x * x * y * z * z * z - std::sqrt(2076.416015625) * x * x * y * y * y * y * y * z - std::sqrt(124178.90625) * x * x * y * y * y * z * z * z - std::sqrt(21262.5) * x * x * y * z * z * z * z * z - std::sqrt(83.056640625) * y * y * y * y * y * y * y * z + std::sqrt(6238.4765625) * y * y * y * y * y * z * z * z - std::sqrt(21262.5) * y * y * y * z * z * z * z * z) + e_2 * (std::sqrt(85050.0) * x * x * x * x * y * z - std::sqrt(765450.0) * x * x * y * y * y * z - std::sqrt(260465.625) * x * x * y * z * z * z - std::sqrt(260465.625) * y * y * y * z * z * z - std::sqrt(85050.0) * y * z * z * z * z * z) + e_3 * (-std::sqrt(643190.625) * x * x * y * z - std::sqrt(643190.625) * y * y * y * z - std::sqrt(2731050.0) * y * z * z * z) + e_4 * (-std::sqrt(4167450.0) * y * z);
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

        pc_36[k] = e_0 * (std::sqrt(0.3076171875) * x * x * x * x * x * x * x * x * x * y - std::sqrt(4.921875) * x * x * x * x * x * x * x * y * y * y - std::sqrt(78.75) * x * x * x * x * x * x * x * y * z * z - std::sqrt(30.76171875) * x * x * x * x * x * y * y * y * y * y + std::sqrt(1968.75) * x * x * x * x * x * y * y * y * z * z + std::sqrt(1107.421875) * x * x * x * x * x * y * z * z * z * z - std::sqrt(4.921875) * x * x * x * y * y * y * y * y * y * y + std::sqrt(1968.75) * x * x * x * y * y * y * y * y * z * z - std::sqrt(39867.1875) * x * x * x * y * y * y * z * z * z * z + std::sqrt(0.3076171875) * x * y * y * y * y * y * y * y * y * y - std::sqrt(78.75) * x * y * y * y * y * y * y * y * z * z + std::sqrt(1107.421875) * x * y * y * y * y * y * z * z * z * z) + e_1 * (std::sqrt(19.6875) * x * x * x * x * x * x * x * y - std::sqrt(3327.1875) * x * x * x * x * x * y * y * y + std::sqrt(21439.6875) * x * x * x * x * x * y * z * z - std::sqrt(3327.1875) * x * x * x * y * y * y * y * y - std::sqrt(96468.75) * x * x * x * y * y * y * z * z - std::sqrt(70875.0) * x * x * x * y * z * z * z * z + std::sqrt(19.6875) * x * y * y * y * y * y * y * y + std::sqrt(21439.6875) * x * y * y * y * y * y * z * z - std::sqrt(70875.0) * x * y * y * y * z * z * z * z) + e_2 * (std::sqrt(1107.421875) * x * x * x * x * x * y - std::sqrt(535992.1875) * x * x * x * y * y * y - std::sqrt(283500.0) * x * x * x * y * z * z + std::sqrt(1107.421875) * x * y * y * y * y * y - std::sqrt(283500.0) * x * y * y * y * z * z - std::sqrt(637875.0) * x * y * z * z * z * z) + e_3 * (-std::sqrt(637875.0) * x * x * x * y - std::sqrt(637875.0) * x * y * y * y - std::sqrt(7087500.0) * x * y * z * z) + e_4 * (-std::sqrt(3472875.0) * x * y);
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

        pc_37[k] = e_0 * (-std::sqrt(6.767578125) * x * x * x * x * x * x * x * x * y * z + std::sqrt(548.173828125) * x * x * x * x * x * x * y * y * y * z + std::sqrt(243.6328125) * x * x * x * x * x * x * y * z * z * z + std::sqrt(169.189453125) * x * x * x * x * y * y * y * y * y * z - std::sqrt(24363.28125) * x * x * x * x * y * y * y * z * z * z - std::sqrt(169.189453125) * x * x * y * y * y * y * y * y * y * z + std::sqrt(6090.8203125) * x * x * y * y * y * y * y * z * z * z) + e_1 * (std::sqrt(1955.830078125) * x * x * x * x * x * x * y * z + std::sqrt(169.189453125) * x * x * x * x * y * y * y * z - std::sqrt(54817.3828125) * x * x * x * x * y * z * z * z + std::sqrt(1522.705078125) * x * x * y * y * y * y * y * z - std::sqrt(24363.28125) * x * x * y * y * y * z * z * z - std::sqrt(169.189453125) * y * y * y * y * y * y * y * z + std::sqrt(6090.8203125) * y * y * y * y * y * z * z * z) + e_2 * (-std::sqrt(877078.125) * x * x * y * z * z * z + std::sqrt(97453.125) * y * y * y * z * z * z) + e_3 * (-std::sqrt(877078.125) * x * x * y * z + std::sqrt(97453.125) * y * y * y * z);

        pc_38[k] = e_0 * (-std::sqrt(0.56396484375) * x * x * x * x * x * x * x * x * x * y + std::sqrt(110.537109375) * x * x * x * x * x * x * x * y * y * y + std::sqrt(20.302734375) * x * x * x * x * x * x * x * y * z * z - std::sqrt(4568.115234375) * x * x * x * x * x * y * y * y * z * z - std::sqrt(110.537109375) * x * x * x * y * y * y * y * y * y * y + std::sqrt(4568.115234375) * x * x * x * y * y * y * y * y * z * z + std::sqrt(0.56396484375) * x * y * y * y * y * y * y * y * y * y - std::sqrt(20.302734375) * x * y * y * y * y * y * y * y * z * z) + e_1 * (std::sqrt(81.2109375) * x * x * x * x * x * x * x * y + std::sqrt(23469.9609375) * x * x * x * x * x * y * y * y - std::sqrt(11694.375) * x * x * x * x * x * y * z * z - std::sqrt(23469.9609375) * x * x * x * y * y * y * y * y - std::sqrt(81.2109375) * x * y * y * y * y * y * y * y + std::sqrt(11694.375) * x * y * y * y * y * y * z * z) + e_2 * (std::sqrt(73089.84375) * x * x * x * x * x * y - std::sqrt(292359.375) * x * x * x * y * z * z - std::sqrt(73089.84375) * x * y * y * y * y * y + std::sqrt(292359.375) * x * y * y * y * z * z) + e_3 * (std::sqrt(519750.0) * x * x * x * y - std::sqrt(519750.0) * x * y * y * y);
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

        pc_39[k] = e_0 * (-std::sqrt(91.3623046875) * x * x * x * x * x * x * x * y * y * z + std::sqrt(497.4169921875) * x * x * x * x * x * y * y * y * y * z + std::sqrt(162.421875) * x * x * x * x * x * y * y * z * z * z + std::sqrt(497.4169921875) * x * x * x * y * y * y * y * y * y * z - std::sqrt(1804.6875) * x * x * x * y * y * y * y * z * z * z - std::sqrt(91.3623046875) * x * y * y * y * y * y * y * y * y * z + std::sqrt(162.421875) * x * y * y * y * y * y * y * z * z * z) + e_1 * (-std::sqrt(91.3623046875) * x * x * x * x * x * x * x * z - std::sqrt(822.2607421875) * x * x * x * x * x * y * y * z + std::sqrt(162.421875) * x * x * x * x * x * z * z * z + std::sqrt(185008.6669921875) * x * x * x * y * y * y * y * z - std::sqrt(16242.1875) * x * x * x * y * y * z * z * z - std::sqrt(26403.7060546875) * x * y * y * y * y * y * y * z + std::sqrt(4060.546875) * x * y * y * y * y * z * z * z) + e_2 * (-std::sqrt(9136.23046875) * x * x * x * x * x * z + std::sqrt(913623.046875) * x * x * x * y * y * z - std::sqrt(228405.76171875) * x * y * y * y * y * z);

        pc_40[k] = e_0 * (-std::sqrt(761.3525390625) * x * x * x * x * x * x * y * y * z * z + std::sqrt(761.3525390625) * x * x * x * x * y * y * y * y * z * z + std::sqrt(1353.515625) * x * x * x * x * y * y * z * z * z * z + std::sqrt(2466.7822265625) * x * x * y * y * y * y * y * y * z * z - std::sqrt(5414.0625) * x * x * y * y * y * y * z * z * z * z - std::sqrt(30.4541015625) * y * y * y * y * y * y * y * y * z * z + std::sqrt(54.140625) * y * y * y * y * y * y * z * z * z * z) + e_1 * (-std::sqrt(761.3525390625) * x * x * x * x * x * x * y * y - std::sqrt(761.3525390625) * x * x * x * x * x * x * z * z + std::sqrt(761.3525390625) * x * x * x * x * y * y * y * y - std::sqrt(761.3525390625) * x * x * x * x * y * y * z * z + std::sqrt(1353.515625) * x * x * x * x * z * z * z * z + std::sqrt(2466.7822265625) * x * x * y * y * y * y * y * y + std::sqrt(220030.8837890625) * x * x * y * y * y * y * z * z - std::sqrt(48726.5625) * x * x * y * y * z * z * z * z - std::sqrt(30.4541015625) * y * y * y * y * y * y * y * y - std::sqrt(3684.9462890625) * y * y * y * y * y * y * z * z + std::sqrt(1353.515625) * y * y * y * y * z * z * z * z) + e_2 * (-std::sqrt(761.3525390625) * x * x * x * x * x * x - std::sqrt(19033.8134765625) * x * x * x * x * y * y - std::sqrt(12181.640625) * x * x * x * x * z * z + std::sqrt(475845.3369140625) * x * x * y * y * y * y + std::sqrt(438539.0625) * x * x * y * y * z * z - std::sqrt(6852.1728515625) * y * y * y * y * y * y - std::sqrt(12181.640625) * y * y * y * y * z * z) + e_3 * (-std::sqrt(48726.5625) * x * x * x * x + std::sqrt(1754156.25) * x * x * y * y - std::sqrt(48726.5625) * y * y * y * y);
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

        pc_41[k] = e_0 * (std::sqrt(22.1484375) * x * x * x * x * x * x * x * y * y * z + std::sqrt(22.1484375) * x * x * x * x * x * y * y * y * y * z - std::sqrt(2844.84375) * x * x * x * x * x * y * y * z * z * z - std::sqrt(22.1484375) * x * x * x * y * y * y * y * y * y * z + std::sqrt(3937.5) * x * x * x * y * y * z * z * z * z * z - std::sqrt(22.1484375) * x * y * y * y * y * y * y * y * y * z + std::sqrt(2844.84375) * x * y * y * y * y * y * y * z * z * z - std::sqrt(3937.5) * x * y * y * y * y * z * z * z * z * z) + e_1 * (std::sqrt(22.1484375) * x * x * x * x * x * x * x * z - std::sqrt(1085.2734375) * x * x * x * x * x * y * y * z - std::sqrt(2844.84375) * x * x * x * x * x * z * z * z - std::sqrt(553.7109375) * x * x * x * y * y * y * y * z + std::sqrt(8859.375) * x * x * x * y * y * z * z * z + std::sqrt(3937.5) * x * x * x * z * z * z * z * z + std::sqrt(199.3359375) * x * y * y * y * y * y * y * z + std::sqrt(29777.34375) * x * y * y * y * y * z * z * z - std::sqrt(35437.5) * x * y * y * z * z * z * z * z) + e_2 * (-std::sqrt(2214.84375) * x * x * x * x * x * z - std::sqrt(8859.375) * x * x * x * y * y * z + std::sqrt(8859.375) * x * x * x * z * z * z + std::sqrt(108527.34375) * x * y * y * y * y * z - std::sqrt(79734.375) * x * y * y * z * z * z) + e_3 * (-std::sqrt(8859.375) * x * x * x * z + std::sqrt(79734.375) * x * y * y * z);
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

        pc_42[k] = e_0 * (std::sqrt(373.7548828125) * x * x * x * x * x * x * y * y * z * z + std::sqrt(1038.2080078125) * x * x * x * x * y * y * y * y * z * z - std::sqrt(5980.078125) * x * x * x * x * y * y * z * z * z * z + std::sqrt(41.5283203125) * x * x * y * y * y * y * y * y * z * z - std::sqrt(2657.8125) * x * x * y * y * y * y * z * z * z * z + std::sqrt(4725.0) * x * x * y * y * z * z * z * z * z * z - std::sqrt(41.5283203125) * y * y * y * y * y * y * y * y * z * z + std::sqrt(664.453125) * y * y * y * y * y * y * z * z * z * z - std::sqrt(525.0) * y * y * y * y * z * z * z * z * z * z) + e_1 * (std::sqrt(373.7548828125) * x * x * x * x * x * x * y * y + std::sqrt(373.7548828125) * x * x * x * x * x * x * z * z + std::sqrt(1038.2080078125) * x * x * x * x * y * y * y * y + std::sqrt(373.7548828125) * x * x * x * x * y * y * z * z - std::sqrt(5980.078125) * x * x * x * x * z * z * z * z + std::sqrt(41.5283203125) * x * x * y * y * y * y * y * y - std::sqrt(373.7548828125) * x * x * y * y * y * y * z * z + std::sqrt(66445.3125) * x * x * y * y * z * z * z * z + std::sqrt(4725.0) * x * x * z * z * z * z * z * z - std::sqrt(41.5283203125) * y * y * y * y * y * y * y * y - std::sqrt(373.7548828125) * y * y * y * y * y * y * z * z - std::sqrt(73.828125) * y * y * y * y * z * z * z * z - std::sqrt(4725.0) * y * y * z * z * z * z * z * z) + e_2 * (std::sqrt(373.7548828125) * x * x * x * x * x * x + std::sqrt(63164.5751953125) * x * x * x * x * y * y - std::sqrt(5980.078125) * x * x * x * x * z * z + std::sqrt(18313.9892578125) * x * x * y * y * y * y + std::sqrt(598007.8125) * x * x * y * y * z * z + std::sqrt(170100.0) * x * x * z * z * z * z - std::sqrt(9343.8720703125) * y * y * y * y * y * y - std::sqrt(32558.203125) * y * y * y * y * z * z - std::sqrt(170100.0) * y * y * z * z * z * z) + e_3 * (std::sqrt(23920.3125) * x * x * x * x + std::sqrt(1063125.0) * x * x * y * y + std::sqrt(861131.25) * x * x * z * z - std::sqrt(248357.8125) * y * y * y * y - std::sqrt(861131.25) * y * y * z * z) + e_4 * (std::sqrt(520931.25) * x * x - std::sqrt(520931.25) * y * y);
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

        pc_43[k] = e_0 * (-std::sqrt(4.6142578125) * x * x * x * x * x * x * x * y * y * z - std::sqrt(41.5283203125) * x * x * x * x * x * y * y * y * y * z + std::sqrt(1386.328125) * x * x * x * x * x * y * y * z * z * z - std::sqrt(41.5283203125) * x * x * x * y * y * y * y * y * y * z + std::sqrt(5545.3125) * x * x * x * y * y * y * y * z * z * z - std::sqrt(6431.25) * x * x * x * y * y * z * z * z * z * z - std::sqrt(4.6142578125) * x * y * y * y * y * y * y * y * y * z + std::sqrt(1386.328125) * x * y * y * y * y * y * y * z * z * z - std::sqrt(6431.25) * x * y * y * y * y * z * z * z * z * z + std::sqrt(2100.0) * x * y * y * z * z * z * z * z * z * z) + e_1 * (-std::sqrt(4.6142578125) * x * x * x * x * x * x * x * z + std::sqrt(779.8095703125) * x * x * x * x * x * y * y * z + std::sqrt(1386.328125) * x * x * x * x * x * z * z * z + std::sqrt(3880.5908203125) * x * x * x * y * y * y * y * z + std::sqrt(295.3125) * x * x * x * y * y * z * z * z - std::sqrt(6431.25) * x * x * x * z * z * z * z * z + std::sqrt(1038.2080078125) * x * y * y * y * y * y * y * z - std::sqrt(401.953125) * x * y * y * y * y * z * z * z + std::sqrt(57881.25) * x * y * y * z * z * z * z * z + std::sqrt(2100.0) * x * z * z * z * z * z * z * z) + e_2 * (std::sqrt(2233.30078125) * x * x * x * x * x * z + std::sqrt(124105.078125) * x * x * x * y * y * z - std::sqrt(42525.0) * x * x * x * z * z * z + std::sqrt(93041.89453125) * x * y * y * y * y * z + std::sqrt(1365525.0) * x * y * y * z * z * z + std::sqrt(231525.0) * x * z * z * z * z * z) + e_3 * (std::sqrt(4725.0) * x * x * x * z + std::sqrt(4540725.0) * x * y * y * z + std::sqrt(3194100.0) * x * z * z * z) + e_4 * (std::sqrt(3704400.0) * x * z);
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

        pc_44[k] = e_0 * (-std::sqrt(46.142578125) * x * x * x * x * x * x * y * y * z * z - std::sqrt(415.283203125) * x * x * x * x * y * y * y * y * z * z + std::sqrt(1312.5) * x * x * x * x * y * y * z * z * z * z - std::sqrt(415.283203125) * x * x * y * y * y * y * y * y * z * z + std::sqrt(5250.0) * x * x * y * y * y * y * z * z * z * z - std::sqrt(2218.125) * x * x * y * y * z * z * z * z * z * z - std::sqrt(46.142578125) * y * y * y * y * y * y * y * y * z * z + std::sqrt(1312.5) * y * y * y * y * y * y * z * z * z * z - std::sqrt(2218.125) * y * y * y * y * z * z * z * z * z * z + std::sqrt(210.0) * y * y * z * z * z * z * z * z * z * z) + e_1 * (-std::sqrt(46.142578125) * x * x * x * x * x * x * y * y - std::sqrt(46.142578125) * x * x * x * x * x * x * z * z - std::sqrt(415.283203125) * x * x * x * x * y * y * y * y - std::sqrt(46.142578125) * x * x * x * x * y * y * z * z + std::sqrt(1312.5) * x * x * x * x * z * z * z * z - std::sqrt(415.283203125) * x * x * y * y * y * y * y * y + std::sqrt(46.142578125) * x * x * y * y * y * y * z * z - std::sqrt(2953.125) * x * x * y * y * z * z * z * z - std::sqrt(2218.125) * x * x * z * z * z * z * z * z - std::sqrt(46.142578125) * y * y * y * y * y * y * y * y + std::sqrt(46.142578125) * y * y * y * y * y * y * z * z - std::sqrt(8203.125) * y * y * y * y * z * z * z * z + std::sqrt(5788.125) * y * y * z * z * z * z * z * z + std::sqrt(210.0) * z * z * z * z * z * z * z * z) + e_2 * (-std::sqrt(46.142578125) * x * x * x * x * x * x - std::sqrt(13335.205078125) * x * x * x * x * y * y + std::sqrt(2953.125) * x * x * x * x * z * z - std::sqrt(44343.017578125) * x * x * y * y * y * y - std::sqrt(26578.125) * x * x * y * y * z * z - std::sqrt(73828.125) * x * x * z * z * z * z - std::sqrt(10382.080078125) * y * y * y * y * y * y - std::sqrt(47250.0) * y * y * y * y * z * z + std::sqrt(73828.125) * y * y * z * z * z * z + std::sqrt(47250.0) * z * z * z * z * z * z) + e_3 * (-std::sqrt(2953.125) * x * x * x * x - std::sqrt(499078.125) * x * x * y * y - std::sqrt(239203.125) * x * x * z * z - std::sqrt(425250.0) * y * y * y * y + std::sqrt(2953.125) * y * y * z * z + std::sqrt(1181250.0) * z * z * z * z) + e_4 * (-std::sqrt(144703.125) * x * x - std::sqrt(1302328.125) * y * y + std::sqrt(2315250.0) * z * z);
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

        pc_45[k] = e_0 * (std::sqrt(0.54931640625) * x * x * x * x * x * x * x * x * y * z + std::sqrt(8.7890625) * x * x * x * x * x * x * y * y * y * z - std::sqrt(205.322265625) * x * x * x * x * x * x * y * z * z * z + std::sqrt(19.775390625) * x * x * x * x * y * y * y * y * y * z - std::sqrt(1847.900390625) * x * x * x * x * y * y * y * z * z * z + std::sqrt(1265.625) * x * x * x * x * y * z * z * z * z * z + std::sqrt(8.7890625) * x * x * y * y * y * y * y * y * y * z - std::sqrt(1847.900390625) * x * x * y * y * y * y * y * z * z * z + std::sqrt(5062.5) * x * x * y * y * y * z * z * z * z * z - std::sqrt(680.625) * x * x * y * z * z * z * z * z * z * z + std::sqrt(0.54931640625) * y * y * y * y * y * y * y * y * y * z - std::sqrt(205.322265625) * y * y * y * y * y * y * y * z * z * z + std::sqrt(1265.625) * y * y * y * y * y * z * z * z * z * z - std::sqrt(680.625) * y * y * y * z * z * z * z * z * z * z + std::sqrt(10.0) * y * z * z * z * z * z * z * z * z * z) + e_1 * (-std::sqrt(177.978515625) * x * x * x * x * x * x * y * z - std::sqrt(1601.806640625) * x * x * x * x * y * y * y * z + std::sqrt(140.625) * x * x * x * x * y * z * z * z - std::sqrt(1601.806640625) * x * x * y * y * y * y * y * z + std::sqrt(562.5) * x * x * y * y * y * z * z * z - std::sqrt(14630.625) * x * x * y * z * z * z * z * z - std::sqrt(177.978515625) * y * y * y * y * y * y * y * z + std::sqrt(140.625) * y * y * y * y * y * z * z * z - std::sqrt(14630.625) * y * y * y * z * z * z * z * z + std::sqrt(90.0) * y * z * z * z * z * z * z * z) + e_2 * (-std::sqrt(20250.0) * x * x * x * x * y * z - std::sqrt(81000.0) * x * x * y * y * y * z - std::sqrt(284765.625) * x * x * y * z * z * z - std::sqrt(20250.0) * y * y * y * y * y * z - std::sqrt(284765.625) * y * y * y * z * z * z - std::sqrt(20250.0) * y * z * z * z * z * z) + e_3 * (-std::sqrt(1216265.625) * x * x * y * z - std::sqrt(1216265.625) * y * y * y * z - std::sqrt(1406250.0) * y * z * z * z) + e_4 * (-std::sqrt(3969000.0) * y * z);
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

        pc_46[k] = e_0 * (-std::sqrt(46.142578125) * x * x * x * x * x * x * x * y * z * z - std::sqrt(415.283203125) * x * x * x * x * x * y * y * y * z * z + std::sqrt(1312.5) * x * x * x * x * x * y * z * z * z * z - std::sqrt(415.283203125) * x * x * x * y * y * y * y * y * z * z + std::sqrt(5250.0) * x * x * x * y * y * y * z * z * z * z - std::sqrt(2218.125) * x * x * x * y * z * z * z * z * z * z - std::sqrt(46.142578125) * x * y * y * y * y * y * y * y * z * z + std::sqrt(1312.5) * x * y * y * y * y * y * z * z * z * z - std::sqrt(2218.125) * x * y * y * y * z * z * z * z * z * z + std::sqrt(210.0) * x * y * z * z * z * z * z * z * z * z) + e_1 * (-std::sqrt(46.142578125) * x * x * x * x * x * x * x * y - std::sqrt(415.283203125) * x * x * x * x * x * y * y * y + std::sqrt(184.5703125) * x * x * x * x * x * y * z * z - std::sqrt(415.283203125) * x * x * x * y * y * y * y * y + std::sqrt(738.28125) * x * x * x * y * y * y * z * z - std::sqrt(16078.125) * x * x * x * y * z * z * z * z - std::sqrt(46.142578125) * x * y * y * y * y * y * y * y + std::sqrt(184.5703125) * x * y * y * y * y * y * z * z - std::sqrt(16078.125) * x * y * y * y * z * z * z * z + std::sqrt(15172.5) * x * y * z * z * z * z * z * z) + e_2 * (-std::sqrt(9043.9453125) * x * x * x * x * x * y - std::sqrt(36175.78125) * x * x * x * y * y * y - std::sqrt(73828.125) * x * x * x * y * z * z - std::sqrt(9043.9453125) * x * y * y * y * y * y - std::sqrt(73828.125) * x * y * y * y * z * z + std::sqrt(295312.5) * x * y * z * z * z * z) + e_3 * (-std::sqrt(357328.125) * x * x * x * y - std::sqrt(357328.125) * x * y * y * y + std::sqrt(295312.5) * x * y * z * z) + e_4 * (-std::sqrt(578812.5) * x * y);
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

        pc_47[k] = e_0 * (-std::sqrt(1.153564453125) * x * x * x * x * x * x * x * x * y * z - std::sqrt(4.6142578125) * x * x * x * x * x * x * y * y * y * z + std::sqrt(346.58203125) * x * x * x * x * x * x * y * z * z * z + std::sqrt(346.58203125) * x * x * x * x * y * y * y * z * z * z - std::sqrt(1607.8125) * x * x * x * x * y * z * z * z * z * z + std::sqrt(4.6142578125) * x * x * y * y * y * y * y * y * y * z - std::sqrt(346.58203125) * x * x * y * y * y * y * y * z * z * z + std::sqrt(525.0) * x * x * y * z * z * z * z * z * z * z + std::sqrt(1.153564453125) * y * y * y * y * y * y * y * y * y * z - std::sqrt(346.58203125) * y * y * y * y * y * y * y * z * z * z + std::sqrt(1607.8125) * y * y * y * y * y * z * z * z * z * z - std::sqrt(525.0) * y * y * y * z * z * z * z * z * z * z) + e_1 * (std::sqrt(373.7548828125) * x * x * x * x * x * x * y * z + std::sqrt(558.3251953125) * x * x * x * x * y * y * y * z - std::sqrt(4339.453125) * x * x * x * x * y * z * z * z - std::sqrt(115.3564453125) * x * x * y * y * y * y * y * z - std::sqrt(5545.3125) * x * x * y * y * y * z * z * z + std::sqrt(57881.25) * x * x * y * z * z * z * z * z - std::sqrt(226.0986328125) * y * y * y * y * y * y * y * z - std::sqrt(73.828125) * y * y * y * y * y * z * z * z - std::sqrt(6431.25) * y * y * y * z * z * z * z * z - std::sqrt(2100.0) * y * z * z * z * z * z * z * z) + e_2 * (std::sqrt(6662.98828125) * x * x * x * x * y * z - std::sqrt(8933.203125) * x * x * y * y * y * z + std::sqrt(798525.0) * x * x * y * z * z * z - std::sqrt(31026.26953125) * y * y * y * y * y * z - std::sqrt(231525.0) * y * y * y * z * z * z - std::sqrt(231525.0) * y * z * z * z * z * z) + e_3 * (std::sqrt(926100.0) * x * x * y * z - std::sqrt(1209600.0) * y * y * y * z - std::sqrt(3194100.0) * y * z * z * z) + e_4 * (-std::sqrt(3704400.0) * y * z);
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

        pc_48[k] = e_0 * (std::sqrt(41.5283203125) * x * x * x * x * x * x * x * y * z * z - std::sqrt(41.5283203125) * x * x * x * x * x * y * y * y * z * z - std::sqrt(664.453125) * x * x * x * x * x * y * z * z * z * z - std::sqrt(1038.2080078125) * x * x * x * y * y * y * y * y * z * z + std::sqrt(2657.8125) * x * x * x * y * y * y * z * z * z * z + std::sqrt(525.0) * x * x * x * y * z * z * z * z * z * z - std::sqrt(373.7548828125) * x * y * y * y * y * y * y * y * z * z + std::sqrt(5980.078125) * x * y * y * y * y * y * z * z * z * z - std::sqrt(4725.0) * x * y * y * y * z * z * z * z * z * z) + e_1 * (std::sqrt(41.5283203125) * x * x * x * x * x * x * x * y - std::sqrt(41.5283203125) * x * x * x * x * x * y * y * y - std::sqrt(1495.01953125) * x * x * x * x * x * y * z * z - std::sqrt(1038.2080078125) * x * x * x * y * y * y * y * y - std::sqrt(5980.078125) * x * x * x * y * y * y * z * z + std::sqrt(57881.25) * x * x * x * y * z * z * z * z - std::sqrt(373.7548828125) * x * y * y * y * y * y * y * y - std::sqrt(1495.01953125) * x * y * y * y * y * y * z * z - std::sqrt(10631.25) * x * y * y * y * z * z * z * z - std::sqrt(18900.0) * x * y * z * z * z * z * z * z) + e_2 * (std::sqrt(1495.01953125) * x * x * x * x * x * y - std::sqrt(53820.703125) * x * x * x * y * y * y + std::sqrt(170100.0) * x * x * x * y * z * z - std::sqrt(73255.95703125) * x * y * y * y * y * y - std::sqrt(382725.0) * x * y * y * y * z * z - std::sqrt(680400.0) * x * y * z * z * z * z) + e_3 * (std::sqrt(1181.25) * x * x * x * y - std::sqrt(1796681.25) * x * y * y * y - std::sqrt(3444525.0) * x * y * z * z) + e_4 * (-std::sqrt(2083725.0) * x * y);
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

        pc_49[k] = e_0 * (std::sqrt(1.38427734375) * x * x * x * x * x * x * x * x * y * z - std::sqrt(22.1484375) * x * x * x * x * x * x * y * y * y * z - std::sqrt(177.802734375) * x * x * x * x * x * x * y * z * z * z - std::sqrt(138.427734375) * x * x * x * x * y * y * y * y * y * z + std::sqrt(4445.068359375) * x * x * x * x * y * y * y * z * z * z + std::sqrt(246.09375) * x * x * x * x * y * z * z * z * z * z - std::sqrt(22.1484375) * x * x * y * y * y * y * y * y * y * z + std::sqrt(4445.068359375) * x * x * y * y * y * y * y * z * z * z - std::sqrt(8859.375) * x * x * y * y * y * z * z * z * z * z + std::sqrt(1.38427734375) * y * y * y * y * y * y * y * y * y * z - std::sqrt(177.802734375) * y * y * y * y * y * y * y * z * z * z + std::sqrt(246.09375) * y * y * y * y * y * z * z * z * z * z) + e_1 * (-std::sqrt(448.505859375) * x * x * x * x * x * x * y * z + std::sqrt(138.427734375) * x * x * x * x * y * y * y * z + std::sqrt(24609.375) * x * x * x * x * y * z * z * z + std::sqrt(935.771484375) * x * x * y * y * y * y * y * z + std::sqrt(15750.0) * x * x * y * y * y * z * z * z - std::sqrt(35437.5) * x * x * y * z * z * z * z * z - std::sqrt(5.537109375) * y * y * y * y * y * y * y * z - std::sqrt(3189.375) * y * y * y * y * y * z * z * z + std::sqrt(3937.5) * y * y * y * z * z * z * z * z) + e_2 * (std::sqrt(8859.375) * x * x * x * x * y * z + std::sqrt(141750.0) * x * x * y * y * y * z - std::sqrt(79734.375) * x * x * y * z * z * z - std::sqrt(8859.375) * y * y * y * y * y * z + std::sqrt(8859.375) * y * y * y * z * z * z) + e_3 * (std::sqrt(79734.375) * x * x * y * z - std::sqrt(8859.375) * y * y * y * z);
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

        pc_50[k] = e_0 * (-std::sqrt(30.4541015625) * x * x * x * x * x * x * x * y * z * z + std::sqrt(2466.7822265625) * x * x * x * x * x * y * y * y * z * z + std::sqrt(54.140625) * x * x * x * x * x * y * z * z * z * z + std::sqrt(761.3525390625) * x * x * x * y * y * y * y * y * z * z - std::sqrt(5414.0625) * x * x * x * y * y * y * z * z * z * z - std::sqrt(761.3525390625) * x * y * y * y * y * y * y * y * z * z + std::sqrt(1353.515625) * x * y * y * y * y * y * z * z * z * z) + e_1 * (-std::sqrt(30.4541015625) * x * x * x * x * x * x * x * y + std::sqrt(2466.7822265625) * x * x * x * x * x * y * y * y + std::sqrt(5969.00390625) * x * x * x * x * x * y * z * z + std::sqrt(761.3525390625) * x * x * x * y * y * y * y * y + std::sqrt(109634.765625) * x * x * x * y * y * y * z * z - std::sqrt(21656.25) * x * x * x * y * z * z * z * z - std::sqrt(761.3525390625) * x * y * y * y * y * y * y * y - std::sqrt(76135.25390625) * x * y * y * y * y * y * z * z + std::sqrt(21656.25) * x * y * y * y * z * z * z * z) + e_2 * (std::sqrt(3045.41015625) * x * x * x * x * x * y + std::sqrt(304541.015625) * x * x * x * y * y * y + std::sqrt(194906.25) * x * x * x * y * z * z - std::sqrt(149225.09765625) * x * y * y * y * y * y - std::sqrt(194906.25) * x * y * y * y * z * z) + e_3 * (std::sqrt(779625.0) * x * x * x * y - std::sqrt(779625.0) * x * y * y * y);

        pc_51[k] = e_0 * (-std::sqrt(2.537841796875) * x * x * x * x * x * x * x * x * y * z + std::sqrt(497.4169921875) * x * x * x * x * x * x * y * y * y * z + std::sqrt(4.51171875) * x * x * x * x * x * x * y * z * z * z - std::sqrt(1015.13671875) * x * x * x * x * y * y * y * z * z * z - std::sqrt(497.4169921875) * x * x * y * y * y * y * y * y * y * z + std::sqrt(1015.13671875) * x * x * y * y * y * y * y * z * z * z + std::sqrt(2.537841796875) * y * y * y * y * y * y * y * y * y * z - std::sqrt(4.51171875) * y * y * y * y * y * y * y * z * z * z) + e_1 * (std::sqrt(822.2607421875) * x * x * x * x * x * x * y * z + std::sqrt(57101.4404296875) * x * x * x * x * y * y * y * z - std::sqrt(4060.546875) * x * x * x * x * y * z * z * z - std::sqrt(138962.0654296875) * x * x * y * y * y * y * y * z + std::sqrt(16242.1875) * x * x * y * y * y * z * z * z + std::sqrt(822.2607421875) * y * y * y * y * y * y * y * z - std::sqrt(162.421875) * y * y * y * y * y * z * z * z) + e_2 * (std::sqrt(228405.76171875) * x * x * x * x * y * z - std::sqrt(913623.046875) * x * x * y * y * y * z + std::sqrt(9136.23046875) * y * y * y * y * y * z);
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

        pc_52[k] = e_0 * (std::sqrt(2.2840576171875) * x * x * x * x * x * x * x * x * x * y - std::sqrt(4.060546875) * x * x * x * x * x * x * x * y * y * y - std::sqrt(146.1796875) * x * x * x * x * x * x * x * y * z * z - std::sqrt(49.74169921875) * x * x * x * x * x * y * y * y * y * y + std::sqrt(795.8671875) * x * x * x * x * x * y * y * y * z * z + std::sqrt(16.2421875) * x * x * x * x * x * y * z * z * z * z - std::sqrt(4.060546875) * x * x * x * y * y * y * y * y * y * y + std::sqrt(795.8671875) * x * x * x * y * y * y * y * y * z * z - std::sqrt(180.46875) * x * x * x * y * y * y * z * z * z * z + std::sqrt(2.2840576171875) * x * y * y * y * y * y * y * y * y * y - std::sqrt(146.1796875) * x * y * y * y * y * y * y * y * z * z + std::sqrt(16.2421875) * x * y * y * y * y * y * z * z * z * z) + e_1 * (std::sqrt(1315.6171875) * x * x * x * x * x * x * x * y - std::sqrt(7162.8046875) * x * x * x * x * x * y * y * y - std::sqrt(21049.875) * x * x * x * x * x * y * z * z - std::sqrt(7162.8046875) * x * x * x * y * y * y * y * y + std::sqrt(233887.5) * x * x * x * y * y * y * z * z + std::sqrt(1315.6171875) * x * y * y * y * y * y * y * y - std::sqrt(21049.875) * x * y * y * y * y * y * z * z) + e_2 * (std::sqrt(32890.4296875) * x * x * x * x * x * y - std::sqrt(365449.21875) * x * x * x * y * y * y + std::sqrt(32890.4296875) * x * y * y * y * y * y);

        pc_53[k] = e_0 * (std::sqrt(19.0338134765625) * x * x * x * x * x * x * x * x * y * z - std::sqrt(1218.1640625) * x * x * x * x * x * x * y * z * z * z - std::sqrt(149.22509765625) * x * x * x * x * y * y * y * y * y * z + std::sqrt(1218.1640625) * x * x * x * x * y * y * y * z * z * z + std::sqrt(135.3515625) * x * x * x * x * y * z * z * z * z * z - std::sqrt(48.7265625) * x * x * y * y * y * y * y * y * y * z + std::sqrt(3946.8515625) * x * x * y * y * y * y * y * z * z * z - std::sqrt(541.40625) * x * x * y * y * y * z * z * z * z * z + std::sqrt(0.7613525390625) * y * y * y * y * y * y * y * y * y * z - std::sqrt(48.7265625) * y * y * y * y * y * y * y * z * z * z + std::sqrt(5.4140625) * y * y * y * y * y * z * z * z * z * z) + e_1 * (std::sqrt(304.541015625) * x * x * x * x * x * x * y * z - std::sqrt(304.541015625) * x * x * x * x * y * y * y * z - std::sqrt(91497.65625) * x * x * x * x * y * z * z * z - std::sqrt(986.712890625) * x * x * y * y * y * y * y * z + std::sqrt(365990.625) * x * x * y * y * y * z * z * z + std::sqrt(12.181640625) * y * y * y * y * y * y * y * z - std::sqrt(3659.90625) * y * y * y * y * y * z * z * z) + e_2 * (-std::sqrt(121816.40625) * x * x * x * x * y * z + std::sqrt(487265.625) * x * x * y * y * y * z - std::sqrt(4872.65625) * y * y * y * y * y * z);
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

        pc_54[k] = e_0 * (-std::sqrt(0.5537109375) * x * x * x * x * x * x * x * x * x * y - std::sqrt(2.21484375) * x * x * x * x * x * x * x * y * y * y + std::sqrt(179.40234375) * x * x * x * x * x * x * x * y * z * z + std::sqrt(179.40234375) * x * x * x * x * x * y * y * y * z * z - std::sqrt(3783.9375) * x * x * x * x * x * y * z * z * z * z + std::sqrt(2.21484375) * x * x * x * y * y * y * y * y * y * y - std::sqrt(179.40234375) * x * x * x * y * y * y * y * y * z * z + std::sqrt(393.75) * x * x * x * y * z * z * z * z * z * z + std::sqrt(0.5537109375) * x * y * y * y * y * y * y * y * y * y - std::sqrt(179.40234375) * x * y * y * y * y * y * y * y * z * z + std::sqrt(3783.9375) * x * y * y * y * y * y * z * z * z * z - std::sqrt(393.75) * x * y * y * y * z * z * z * z * z * z) + e_1 * (-std::sqrt(318.9375) * x * x * x * x * x * x * x * y - std::sqrt(318.9375) * x * x * x * x * x * y * y * y - std::sqrt(2268.0) * x * x * x * x * x * y * z * z + std::sqrt(318.9375) * x * x * x * y * y * y * y * y - std::sqrt(100800.0) * x * x * x * y * z * z * z * z + std::sqrt(318.9375) * x * y * y * y * y * y * y * y + std::sqrt(2268.0) * x * y * y * y * y * y * z * z + std::sqrt(100800.0) * x * y * y * y * z * z * z * z) + e_2 * (-std::sqrt(56700.0) * x * x * x * x * x * y - std::sqrt(1417500.0) * x * x * x * y * z * z + std::sqrt(56700.0) * x * y * y * y * y * y + std::sqrt(1417500.0) * x * y * y * y * z * z) + e_3 * (-std::sqrt(1417500.0) * x * x * x * y + std::sqrt(1417500.0) * x * y * y * y);
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

        pc_55[k] = e_0 * (-std::sqrt(9.3438720703125) * x * x * x * x * x * x * x * x * y * z - std::sqrt(66.4453125) * x * x * x * x * x * x * y * y * y * z + std::sqrt(1063.125) * x * x * x * x * x * x * y * z * z * z - std::sqrt(37.37548828125) * x * x * x * x * y * y * y * y * y * z + std::sqrt(2953.125) * x * x * x * x * y * y * y * z * z * z - std::sqrt(5382.0703125) * x * x * x * x * y * z * z * z * z * z + std::sqrt(118.125) * x * x * y * y * y * y * y * z * z * z - std::sqrt(2392.03125) * x * x * y * y * y * z * z * z * z * z + std::sqrt(472.5) * x * x * y * z * z * z * z * z * z * z + std::sqrt(1.0382080078125) * y * y * y * y * y * y * y * y * y * z - std::sqrt(118.125) * y * y * y * y * y * y * y * z * z * z + std::sqrt(598.0078125) * y * y * y * y * y * z * z * z * z * z - std::sqrt(52.5) * y * y * y * z * z * z * z * z * z * z) + e_1 * (-std::sqrt(149.501953125) * x * x * x * x * x * x * y * z - std::sqrt(415.283203125) * x * x * x * x * y * y * y * z - std::sqrt(6644.53125) * x * x * x * x * y * z * z * z - std::sqrt(16.611328125) * x * x * y * y * y * y * y * z - std::sqrt(2953.125) * x * x * y * y * y * z * z * z - std::sqrt(17010.0) * x * x * y * z * z * z * z * z + std::sqrt(16.611328125) * y * y * y * y * y * y * y * z + std::sqrt(738.28125) * y * y * y * y * y * z * z * z + std::sqrt(1890.0) * y * y * y * z * z * z * z * z) + e_2 * (-std::sqrt(59800.78125) * x * x * x * x * y * z - std::sqrt(26578.125) * x * x * y * y * y * z - std::sqrt(956812.5) * x * x * y * z * z * z + std::sqrt(6644.53125) * y * y * y * y * y * z + std::sqrt(106312.5) * y * y * y * z * z * z) + e_3 * (-std::sqrt(2657812.5) * x * x * y * z + std::sqrt(295312.5) * y * y * y * z);
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

        pc_56[k] = e_0 * (std::sqrt(0.1153564453125) * x * x * x * x * x * x * x * x * x * y + std::sqrt(1.845703125) * x * x * x * x * x * x * x * y * y * y - std::sqrt(66.4453125) * x * x * x * x * x * x * x * y * z * z + std::sqrt(4.15283203125) * x * x * x * x * x * y * y * y * y * y - std::sqrt(598.0078125) * x * x * x * x * x * y * y * y * z * z + std::sqrt(2481.4453125) * x * x * x * x * x * y * z * z * z * z + std::sqrt(1.845703125) * x * x * x * y * y * y * y * y * y * y - std::sqrt(598.0078125) * x * x * x * y * y * y * y * y * z * z + std::sqrt(9925.78125) * x * x * x * y * y * y * z * z * z * z - std::sqrt(3360.0) * x * x * x * y * z * z * z * z * z * z + std::sqrt(0.1153564453125) * x * y * y * y * y * y * y * y * y * y - std::sqrt(66.4453125) * x * y * y * y * y * y * y * y * z * z + std::sqrt(2481.4453125) * x * y * y * y * y * y * z * z * z * z - std::sqrt(3360.0) * x * y * y * y * z * z * z * z * z * z + std::sqrt(210.0) * x * y * z * z * z * z * z * z * z * z) + e_1 * (std::sqrt(66.4453125) * x * x * x * x * x * x * x * y + std::sqrt(598.0078125) * x * x * x * x * x * y * y * y + std::sqrt(2953.125) * x * x * x * x * x * y * z * z + std::sqrt(598.0078125) * x * x * x * y * y * y * y * y + std::sqrt(11812.5) * x * x * x * y * y * y * z * z - std::sqrt(5250.0) * x * x * x * y * z * z * z * z + std::sqrt(66.4453125) * x * y * y * y * y * y * y * y + std::sqrt(2953.125) * x * y * y * y * y * y * z * z - std::sqrt(5250.0) * x * y * y * y * z * z * z * z + std::sqrt(3360.0) * x * y * z * z * z * z * z * z) + e_2 * (std::sqrt(22333.0078125) * x * x * x * x * x * y + std::sqrt(89332.03125) * x * x * x * y * y * y + std::sqrt(47250.0) * x * x * x * y * z * z + std::sqrt(22333.0078125) * x * y * y * y * y * y + std::sqrt(47250.0) * x * y * y * y * z * z + std::sqrt(47250.0) * x * y * z * z * z * z) + e_3 * (std::sqrt(756000.0) * x * x * x * y + std::sqrt(756000.0) * x * y * y * y + std::sqrt(756000.0) * x * y * z * z) + e_4 * (std::sqrt(2315250.0) * x * y);
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

        pc_57[k] = e_0 * (std::sqrt(1.153564453125) * x * x * x * x * x * x * x * x * y * z + std::sqrt(18.45703125) * x * x * x * x * x * x * y * y * y * z - std::sqrt(166.11328125) * x * x * x * x * x * x * y * z * z * z + std::sqrt(41.5283203125) * x * x * x * x * y * y * y * y * y * z - std::sqrt(1495.01953125) * x * x * x * x * y * y * y * z * z * z + std::sqrt(1517.25) * x * x * x * x * y * z * z * z * z * z + std::sqrt(18.45703125) * x * x * y * y * y * y * y * y * y * z - std::sqrt(1495.01953125) * x * x * y * y * y * y * y * z * z * z + std::sqrt(6069.0) * x * x * y * y * y * z * z * z * z * z - std::sqrt(635.25) * x * x * y * z * z * z * z * z * z * z + std::sqrt(1.153564453125) * y * y * y * y * y * y * y * y * y * z - std::sqrt(166.11328125) * y * y * y * y * y * y * y * z * z * z + std::sqrt(1517.25) * y * y * y * y * y * z * z * z * z * z - std::sqrt(635.25) * y * y * y * z * z * z * z * z * z * z + std::sqrt(21.0) * y * z * z * z * z * z * z * z * z * z) + e_1 * (std::sqrt(18.45703125) * x * x * x * x * x * x * y * z + std::sqrt(166.11328125) * x * x * x * x * y * y * y * z + std::sqrt(6431.25) * x * x * x * x * y * z * z * z + std::sqrt(166.11328125) * x * x * y * y * y * y * y * z + std::sqrt(25725.0) * x * x * y * y * y * z * z * z - std::sqrt(3827.25) * x * x * y * z * z * z * z * z + std::sqrt(18.45703125) * y * y * y * y * y * y * y * z + std::sqrt(6431.25) * y * y * y * y * y * z * z * z - std::sqrt(3827.25) * y * y * y * z * z * z * z * z + std::sqrt(4116.0) * y * z * z * z * z * z * z * z) + e_2 * (std::sqrt(29531.25) * x * x * x * x * y * z + std::sqrt(118125.0) * x * x * y * y * y * z + std::sqrt(29531.25) * x * x * y * z * z * z + std::sqrt(29531.25) * y * y * y * y * y * z + std::sqrt(29531.25) * y * y * y * z * z * z + std::sqrt(302400.0) * y * z * z * z * z * z) + e_3 * (std::sqrt(738281.25) * x * x * y * z + std::sqrt(738281.25) * y * y * y * z + std::sqrt(4252500.0) * y * z * z * z) + e_4 * (std::sqrt(5788125.0) * y * z);
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

        pc_58[k] = e_0 * (-0.1171875 * x * x * x * x * x * x * x * x * x * x - 0.5859375 * x * x * x * x * x * x * x * x * y * y + 3.046875 * x * x * x * x * x * x * x * x * z * z - 1.171875 * x * x * x * x * x * x * y * y * y * y + 12.1875 * x * x * x * x * x * x * y * y * z * z - 20.0 * x * x * x * x * x * x * z * z * z * z - 1.171875 * x * x * x * x * y * y * y * y * y * y + 18.28125 * x * x * x * x * y * y * y * y * z * z - 60.0 * x * x * x * x * y * y * z * z * z * z + 28.5 * x * x * x * x * z * z * z * z * z * z - 0.5859375 * x * x * y * y * y * y * y * y * y * y + 12.1875 * x * x * y * y * y * y * y * y * z * z - 60.0 * x * x * y * y * y * y * z * z * z * z + 57.0 * x * x * y * y * z * z * z * z * z * z - 10.5 * x * x * z * z * z * z * z * z * z * z - 0.1171875 * y * y * y * y * y * y * y * y * y * y + 3.046875 * y * y * y * y * y * y * y * y * z * z - 20.0 * y * y * y * y * y * y * z * z * z * z + 28.5 * y * y * y * y * z * z * z * z * z * z - 10.5 * y * y * z * z * z * z * z * z * z * z + z * z * z * z * z * z * z * z * z * z) + e_1 * (-2.8125 * x * x * x * x * x * x * x * x - 11.25 * x * x * x * x * x * x * y * y - 22.5 * x * x * x * x * x * x * z * z - 16.875 * x * x * x * x * y * y * y * y - 67.5 * x * x * x * x * y * y * z * z + 67.5 * x * x * x * x * z * z * z * z - 11.25 * x * x * y * y * y * y * y * y - 67.5 * x * x * y * y * y * y * z * z + 135.0 * x * x * y * y * z * z * z * z - 66.0 * x * x * z * z * z * z * z * z - 2.8125 * y * y * y * y * y * y * y * y - 22.5 * y * y * y * y * y * y * z * z + 67.5 * y * y * y * y * z * z * z * z - 66.0 * y * y * z * z * z * z * z * z + 24.0 * z * z * z * z * z * z * z * z) + e_2 * (-56.25 * x * x * x * x * x * x - 168.75 * x * x * x * x * y * y - 168.75 * x * x * y * y * y * y - 225.0 * x * x * z * z * z * z - 56.25 * y * y * y * y * y * y - 225.0 * y * y * z * z * z * z + 270.0 * z * z * z * z * z * z) + e_3 * (-337.5 * x * x * x * x - 675.0 * x * x * y * y - 450.0 * x * x * z * z - 337.5 * y * y * y * y - 450.0 * y * y * z * z + 1200.0 * z * z * z * z) + e_4 * (-787.5 * x * x - 787.5 * y * y + 1575.0 * z * z);
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

        pc_59[k] = e_0 * (std::sqrt(1.153564453125) * x * x * x * x * x * x * x * x * x * z + std::sqrt(18.45703125) * x * x * x * x * x * x * x * y * y * z - std::sqrt(166.11328125) * x * x * x * x * x * x * x * z * z * z + std::sqrt(41.5283203125) * x * x * x * x * x * y * y * y * y * z - std::sqrt(1495.01953125) * x * x * x * x * x * y * y * z * z * z + std::sqrt(1517.25) * x * x * x * x * x * z * z * z * z * z + std::sqrt(18.45703125) * x * x * x * y * y * y * y * y * y * z - std::sqrt(1495.01953125) * x * x * x * y * y * y * y * z * z * z + std::sqrt(6069.0) * x * x * x * y * y * z * z * z * z * z - std::sqrt(635.25) * x * x * x * z * z * z * z * z * z * z + std::sqrt(1.153564453125) * x * y * y * y * y * y * y * y * y * z - std::sqrt(166.11328125) * x * y * y * y * y * y * y * z * z * z + std::sqrt(1517.25) * x * y * y * y * y * z * z * z * z * z - std::sqrt(635.25) * x * y * y * z * z * z * z * z * z * z + std::sqrt(21.0) * x * z * z * z * z * z * z * z * z * z) + e_1 * (std::sqrt(18.45703125) * x * x * x * x * x * x * x * z + std::sqrt(166.11328125) * x * x * x * x * x * y * y * z + std::sqrt(6431.25) * x * x * x * x * x * z * z * z + std::sqrt(166.11328125) * x * x * x * y * y * y * y * z + std::sqrt(25725.0) * x * x * x * y * y * z * z * z - std::sqrt(3827.25) * x * x * x * z * z * z * z * z + std::sqrt(18.45703125) * x * y * y * y * y * y * y * z + std::sqrt(6431.25) * x * y * y * y * y * z * z * z - std::sqrt(3827.25) * x * y * y * z * z * z * z * z + std::sqrt(4116.0) * x * z * z * z * z * z * z * z) + e_2 * (std::sqrt(29531.25) * x * x * x * x * x * z + std::sqrt(118125.0) * x * x * x * y * y * z + std::sqrt(29531.25) * x * x * x * z * z * z + std::sqrt(29531.25) * x * y * y * y * y * z + std::sqrt(29531.25) * x * y * y * z * z * z + std::sqrt(302400.0) * x * z * z * z * z * z) + e_3 * (std::sqrt(738281.25) * x * x * x * z + std::sqrt(738281.25) * x * y * y * z + std::sqrt(4252500.0) * x * z * z * z) + e_4 * (std::sqrt(5788125.0) * x * z);
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

        pc_60[k] = e_0 * (std::sqrt(0.028839111328125) * x * x * x * x * x * x * x * x * x * x + std::sqrt(0.259552001953125) * x * x * x * x * x * x * x * x * y * y - std::sqrt(16.611328125) * x * x * x * x * x * x * x * x * z * z + std::sqrt(0.1153564453125) * x * x * x * x * x * x * y * y * y * y - std::sqrt(66.4453125) * x * x * x * x * x * x * y * y * z * z + std::sqrt(620.361328125) * x * x * x * x * x * x * z * z * z * z - std::sqrt(0.1153564453125) * x * x * x * x * y * y * y * y * y * y + std::sqrt(620.361328125) * x * x * x * x * y * y * z * z * z * z - std::sqrt(840.0) * x * x * x * x * z * z * z * z * z * z - std::sqrt(0.259552001953125) * x * x * y * y * y * y * y * y * y * y + std::sqrt(66.4453125) * x * x * y * y * y * y * y * y * z * z - std::sqrt(620.361328125) * x * x * y * y * y * y * z * z * z * z + std::sqrt(52.5) * x * x * z * z * z * z * z * z * z * z - std::sqrt(0.028839111328125) * y * y * y * y * y * y * y * y * y * y + std::sqrt(16.611328125) * y * y * y * y * y * y * y * y * z * z - std::sqrt(620.361328125) * y * y * y * y * y * y * z * z * z * z + std::sqrt(840.0) * y * y * y * y * z * z * z * z * z * z - std::sqrt(52.5) * y * y * z * z * z * z * z * z * z * z) + e_1 * (std::sqrt(16.611328125) * x * x * x * x * x * x * x * x + std::sqrt(66.4453125) * x * x * x * x * x * x * y * y + std::sqrt(738.28125) * x * x * x * x * x * x * z * z + std::sqrt(738.28125) * x * x * x * x * y * y * z * z - std::sqrt(1312.5) * x * x * x * x * z * z * z * z - std::sqrt(66.4453125) * x * x * y * y * y * y * y * y - std::sqrt(738.28125) * x * x * y * y * y * y * z * z + std::sqrt(840.0) * x * x * z * z * z * z * z * z - std::sqrt(16.611328125) * y * y * y * y * y * y * y * y - std::sqrt(738.28125) * y * y * y * y * y * y * z * z + std::sqrt(1312.5) * y * y * y * y * z * z * z * z - std::sqrt(840.0) * y * y * z * z * z * z * z * z) + e_2 * (std::sqrt(5583.251953125) * x * x * x * x * x * x + std::sqrt(5583.251953125) * x * x * x * x * y * y + std::sqrt(11812.5) * x * x * x * x * z * z - std::sqrt(5583.251953125) * x * x * y * y * y * y + std::sqrt(11812.5) * x * x * z * z * z * z - std::sqrt(5583.251953125) * y * y * y * y * y * y - std::sqrt(11812.5) * y * y * y * y * z * z - std::sqrt(11812.5) * y * y * z * z * z * z) + e_3 * (std::sqrt(189000.0) * x * x * x * x + std::sqrt(189000.0) * x * x * z * z - std::sqrt(189000.0) * y * y * y * y - std::sqrt(189000.0) * y * y * z * z) + e_4 * (std::sqrt(578812.5) * x * x - std::sqrt(578812.5) * y * y);
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

        pc_61[k] = e_0 * (-std::sqrt(1.0382080078125) * x * x * x * x * x * x * x * x * x * z + std::sqrt(118.125) * x * x * x * x * x * x * x * z * z * z + std::sqrt(37.37548828125) * x * x * x * x * x * y * y * y * y * z - std::sqrt(118.125) * x * x * x * x * x * y * y * z * z * z - std::sqrt(598.0078125) * x * x * x * x * x * z * z * z * z * z + std::sqrt(66.4453125) * x * x * x * y * y * y * y * y * y * z - std::sqrt(2953.125) * x * x * x * y * y * y * y * z * z * z + std::sqrt(2392.03125) * x * x * x * y * y * z * z * z * z * z + std::sqrt(52.5) * x * x * x * z * z * z * z * z * z * z + std::sqrt(9.3438720703125) * x * y * y * y * y * y * y * y * y * z - std::sqrt(1063.125) * x * y * y * y * y * y * y * z * z * z + std::sqrt(5382.0703125) * x * y * y * y * y * z * z * z * z * z - std::sqrt(472.5) * x * y * y * z * z * z * z * z * z * z) + e_1 * (-std::sqrt(16.611328125) * x * x * x * x * x * x * x * z + std::sqrt(16.611328125) * x * x * x * x * x * y * y * z - std::sqrt(738.28125) * x * x * x * x * x * z * z * z + std::sqrt(415.283203125) * x * x * x * y * y * y * y * z + std::sqrt(2953.125) * x * x * x * y * y * z * z * z - std::sqrt(1890.0) * x * x * x * z * z * z * z * z + std::sqrt(149.501953125) * x * y * y * y * y * y * y * z + std::sqrt(6644.53125) * x * y * y * y * y * z * z * z + std::sqrt(17010.0) * x * y * y * z * z * z * z * z) + e_2 * (-std::sqrt(6644.53125) * x * x * x * x * x * z + std::sqrt(26578.125) * x * x * x * y * y * z - std::sqrt(106312.5) * x * x * x * z * z * z + std::sqrt(59800.78125) * x * y * y * y * y * z + std::sqrt(956812.5) * x * y * y * z * z * z) + e_3 * (-std::sqrt(295312.5) * x * x * x * z + std::sqrt(2657812.5) * x * y * y * z);
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

        pc_62[k] = e_0 * (-std::sqrt(0.03460693359375) * x * x * x * x * x * x * x * x * x * x + std::sqrt(0.31146240234375) * x * x * x * x * x * x * x * x * y * y + std::sqrt(11.212646484375) * x * x * x * x * x * x * x * x * z * z + std::sqrt(6.782958984375) * x * x * x * x * x * x * y * y * y * y - std::sqrt(179.40234375) * x * x * x * x * x * x * y * y * z * z - std::sqrt(236.49609375) * x * x * x * x * x * x * z * z * z * z + std::sqrt(6.782958984375) * x * x * x * x * y * y * y * y * y * y - std::sqrt(1121.2646484375) * x * x * x * x * y * y * y * y * z * z + std::sqrt(5912.40234375) * x * x * x * x * y * y * z * z * z * z + std::sqrt(24.609375) * x * x * x * x * z * z * z * z * z * z + std::sqrt(0.31146240234375) * x * x * y * y * y * y * y * y * y * y - std::sqrt(179.40234375) * x * x * y * y * y * y * y * y * z * z + std::sqrt(5912.40234375) * x * x * y * y * y * y * z * z * z * z - std::sqrt(885.9375) * x * x * y * y * z * z * z * z * z * z - std::sqrt(0.03460693359375) * y * y * y * y * y * y * y * y * y * y + std::sqrt(11.212646484375) * y * y * y * y * y * y * y * y * z * z - std::sqrt(236.49609375) * y * y * y * y * y * y * z * z * z * z + std::sqrt(24.609375) * y * y * y * y * z * z * z * z * z * z) + e_1 * (-std::sqrt(19.93359375) * x * x * x * x * x * x * x * x + std::sqrt(318.9375) * x * x * x * x * x * x * y * y - std::sqrt(141.75) * x * x * x * x * x * x * z * z + std::sqrt(1993.359375) * x * x * x * x * y * y * y * y + std::sqrt(3543.75) * x * x * x * x * y * y * z * z - std::sqrt(6300.0) * x * x * x * x * z * z * z * z + std::sqrt(318.9375) * x * x * y * y * y * y * y * y + std::sqrt(3543.75) * x * x * y * y * y * y * z * z + std::sqrt(226800.0) * x * x * y * y * z * z * z * z - std::sqrt(19.93359375) * y * y * y * y * y * y * y * y - std::sqrt(141.75) * y * y * y * y * y * y * z * z - std::sqrt(6300.0) * y * y * y * y * z * z * z * z) + e_2 * (-std::sqrt(3543.75) * x * x * x * x * x * x + std::sqrt(88593.75) * x * x * x * x * y * y - std::sqrt(88593.75) * x * x * x * x * z * z + std::sqrt(88593.75) * x * x * y * y * y * y + std::sqrt(3189375.0) * x * x * y * y * z * z - std::sqrt(3543.75) * y * y * y * y * y * y - std::sqrt(88593.75) * y * y * y * y * z * z) + e_3 * (-std::sqrt(88593.75) * x * x * x * x + std::sqrt(3189375.0) * x * x * y * y - std::sqrt(88593.75) * y * y * y * y);
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

        pc_63[k] = e_0 * (std::sqrt(0.7613525390625) * x * x * x * x * x * x * x * x * x * z - std::sqrt(48.7265625) * x * x * x * x * x * x * x * y * y * z - std::sqrt(48.7265625) * x * x * x * x * x * x * x * z * z * z - std::sqrt(149.22509765625) * x * x * x * x * x * y * y * y * y * z + std::sqrt(3946.8515625) * x * x * x * x * x * y * y * z * z * z + std::sqrt(5.4140625) * x * x * x * x * x * z * z * z * z * z + std::sqrt(1218.1640625) * x * x * x * y * y * y * y * z * z * z - std::sqrt(541.40625) * x * x * x * y * y * z * z * z * z * z + std::sqrt(19.0338134765625) * x * y * y * y * y * y * y * y * y * z - std::sqrt(1218.1640625) * x * y * y * y * y * y * y * z * z * z + std::sqrt(135.3515625) * x * y * y * y * y * z * z * z * z * z) + e_1 * (std::sqrt(12.181640625) * x * x * x * x * x * x * x * z - std::sqrt(986.712890625) * x * x * x * x * x * y * y * z - std::sqrt(3659.90625) * x * x * x * x * x * z * z * z - std::sqrt(304.541015625) * x * x * x * y * y * y * y * z + std::sqrt(365990.625) * x * x * x * y * y * z * z * z + std::sqrt(304.541015625) * x * y * y * y * y * y * y * z - std::sqrt(91497.65625) * x * y * y * y * y * z * z * z) + e_2 * (-std::sqrt(4872.65625) * x * x * x * x * x * z + std::sqrt(487265.625) * x * x * x * y * y * z - std::sqrt(121816.40625) * x * y * y * y * y * z);

        pc_64[k] = e_0 * (std::sqrt(0.063446044921875) * x * x * x * x * x * x * x * x * x * x - std::sqrt(10.722381591796875) * x * x * x * x * x * x * x * x * y * y - std::sqrt(4.060546875) * x * x * x * x * x * x * x * x * z * z - std::sqrt(12.4354248046875) * x * x * x * x * x * x * y * y * y * y + std::sqrt(795.8671875) * x * x * x * x * x * x * y * y * z * z + std::sqrt(0.451171875) * x * x * x * x * x * x * z * z * z * z + std::sqrt(12.4354248046875) * x * x * x * x * y * y * y * y * y * y - std::sqrt(101.513671875) * x * x * x * x * y * y * z * z * z * z + std::sqrt(10.722381591796875) * x * x * y * y * y * y * y * y * y * y - std::sqrt(795.8671875) * x * x * y * y * y * y * y * y * z * z + std::sqrt(101.513671875) * x * x * y * y * y * y * z * z * z * z - std::sqrt(0.063446044921875) * y * y * y * y * y * y * y * y * y * y + std::sqrt(4.060546875) * y * y * y * y * y * y * y * y * z * z - std::sqrt(0.451171875) * y * y * y * y * y * y * z * z * z * z) + e_1 * (std::sqrt(36.544921875) * x * x * x * x * x * x * x * x - std::sqrt(7162.8046875) * x * x * x * x * x * x * y * y - std::sqrt(584.71875) * x * x * x * x * x * x * z * z + std::sqrt(131561.71875) * x * x * x * x * y * y * z * z + std::sqrt(7162.8046875) * x * x * y * y * y * y * y * y - std::sqrt(131561.71875) * x * x * y * y * y * y * z * z - std::sqrt(36.544921875) * y * y * y * y * y * y * y * y + std::sqrt(584.71875) * y * y * y * y * y * y * z * z) + e_2 * (std::sqrt(913.623046875) * x * x * x * x * x * x - std::sqrt(205565.185546875) * x * x * x * x * y * y + std::sqrt(205565.185546875) * x * x * y * y * y * y - std::sqrt(913.623046875) * y * y * y * y * y * y);
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

        pc_65[k] = e_0 * (-std::sqrt(91.3623046875) * x * x * x * x * x * x * x * x * y * z + std::sqrt(497.4169921875) * x * x * x * x * x * x * y * y * y * z + std::sqrt(162.421875) * x * x * x * x * x * x * y * z * z * z + std::sqrt(497.4169921875) * x * x * x * x * y * y * y * y * y * z - std::sqrt(1804.6875) * x * x * x * x * y * y * y * z * z * z - std::sqrt(91.3623046875) * x * x * y * y * y * y * y * y * y * z + std::sqrt(162.421875) * x * x * y * y * y * y * y * z * z * z) + e_1 * (-std::sqrt(26403.7060546875) * x * x * x * x * x * x * y * z + std::sqrt(185008.6669921875) * x * x * x * x * y * y * y * z + std::sqrt(4060.546875) * x * x * x * x * y * z * z * z - std::sqrt(822.2607421875) * x * x * y * y * y * y * y * z - std::sqrt(16242.1875) * x * x * y * y * y * z * z * z - std::sqrt(91.3623046875) * y * y * y * y * y * y * y * z + std::sqrt(162.421875) * y * y * y * y * y * z * z * z) + e_2 * (-std::sqrt(228405.76171875) * x * x * x * x * y * z + std::sqrt(913623.046875) * x * x * y * y * y * z - std::sqrt(9136.23046875) * y * y * y * y * y * z);

        pc_66[k] = e_0 * (-std::sqrt(761.3525390625) * x * x * x * x * x * x * x * y * z * z + std::sqrt(761.3525390625) * x * x * x * x * x * y * y * y * z * z + std::sqrt(1353.515625) * x * x * x * x * x * y * z * z * z * z + std::sqrt(2466.7822265625) * x * x * x * y * y * y * y * y * z * z - std::sqrt(5414.0625) * x * x * x * y * y * y * z * z * z * z - std::sqrt(30.4541015625) * x * y * y * y * y * y * y * y * z * z + std::sqrt(54.140625) * x * y * y * y * y * y * z * z * z * z) + e_1 * (-std::sqrt(761.3525390625) * x * x * x * x * x * x * x * y + std::sqrt(761.3525390625) * x * x * x * x * x * y * y * y - std::sqrt(76135.25390625) * x * x * x * x * x * y * z * z + std::sqrt(2466.7822265625) * x * x * x * y * y * y * y * y + std::sqrt(109634.765625) * x * x * x * y * y * y * z * z + std::sqrt(21656.25) * x * x * x * y * z * z * z * z - std::sqrt(30.4541015625) * x * y * y * y * y * y * y * y + std::sqrt(5969.00390625) * x * y * y * y * y * y * z * z - std::sqrt(21656.25) * x * y * y * y * z * z * z * z) + e_2 * (-std::sqrt(149225.09765625) * x * x * x * x * x * y + std::sqrt(304541.015625) * x * x * x * y * y * y - std::sqrt(194906.25) * x * x * x * y * z * z + std::sqrt(3045.41015625) * x * y * y * y * y * y + std::sqrt(194906.25) * x * y * y * y * z * z) + e_3 * (-std::sqrt(779625.0) * x * x * x * y + std::sqrt(779625.0) * x * y * y * y);
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

        pc_67[k] = e_0 * (std::sqrt(22.1484375) * x * x * x * x * x * x * x * x * y * z + std::sqrt(22.1484375) * x * x * x * x * x * x * y * y * y * z - std::sqrt(2844.84375) * x * x * x * x * x * x * y * z * z * z - std::sqrt(22.1484375) * x * x * x * x * y * y * y * y * y * z + std::sqrt(3937.5) * x * x * x * x * y * z * z * z * z * z - std::sqrt(22.1484375) * x * x * y * y * y * y * y * y * y * z + std::sqrt(2844.84375) * x * x * y * y * y * y * y * z * z * z - std::sqrt(3937.5) * x * x * y * y * y * z * z * z * z * z) + e_1 * (-std::sqrt(199.3359375) * x * x * x * x * x * x * y * z + std::sqrt(553.7109375) * x * x * x * x * y * y * y * z - std::sqrt(29777.34375) * x * x * x * x * y * z * z * z + std::sqrt(1085.2734375) * x * x * y * y * y * y * y * z - std::sqrt(8859.375) * x * x * y * y * y * z * z * z + std::sqrt(35437.5) * x * x * y * z * z * z * z * z - std::sqrt(22.1484375) * y * y * y * y * y * y * y * z + std::sqrt(2844.84375) * y * y * y * y * y * z * z * z - std::sqrt(3937.5) * y * y * y * z * z * z * z * z) + e_2 * (-std::sqrt(108527.34375) * x * x * x * x * y * z + std::sqrt(8859.375) * x * x * y * y * y * z + std::sqrt(79734.375) * x * x * y * z * z * z + std::sqrt(2214.84375) * y * y * y * y * y * z - std::sqrt(8859.375) * y * y * y * z * z * z) + e_3 * (-std::sqrt(79734.375) * x * x * y * z + std::sqrt(8859.375) * y * y * y * z);
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

        pc_68[k] = e_0 * (std::sqrt(373.7548828125) * x * x * x * x * x * x * x * y * z * z + std::sqrt(1038.2080078125) * x * x * x * x * x * y * y * y * z * z - std::sqrt(5980.078125) * x * x * x * x * x * y * z * z * z * z + std::sqrt(41.5283203125) * x * x * x * y * y * y * y * y * z * z - std::sqrt(2657.8125) * x * x * x * y * y * y * z * z * z * z + std::sqrt(4725.0) * x * x * x * y * z * z * z * z * z * z - std::sqrt(41.5283203125) * x * y * y * y * y * y * y * y * z * z + std::sqrt(664.453125) * x * y * y * y * y * y * z * z * z * z - std::sqrt(525.0) * x * y * y * y * z * z * z * z * z * z) + e_1 * (std::sqrt(373.7548828125) * x * x * x * x * x * x * x * y + std::sqrt(1038.2080078125) * x * x * x * x * x * y * y * y + std::sqrt(1495.01953125) * x * x * x * x * x * y * z * z + std::sqrt(41.5283203125) * x * x * x * y * y * y * y * y + std::sqrt(5980.078125) * x * x * x * y * y * y * z * z + std::sqrt(10631.25) * x * x * x * y * z * z * z * z - std::sqrt(41.5283203125) * x * y * y * y * y * y * y * y + std::sqrt(1495.01953125) * x * y * y * y * y * y * z * z - std::sqrt(57881.25) * x * y * y * y * z * z * z * z + std::sqrt(18900.0) * x * y * z * z * z * z * z * z) + e_2 * (std::sqrt(73255.95703125) * x * x * x * x * x * y + std::sqrt(53820.703125) * x * x * x * y * y * y + std::sqrt(382725.0) * x * x * x * y * z * z - std::sqrt(1495.01953125) * x * y * y * y * y * y - std::sqrt(170100.0) * x * y * y * y * z * z + std::sqrt(680400.0) * x * y * z * z * z * z) + e_3 * (std::sqrt(1796681.25) * x * x * x * y - std::sqrt(1181.25) * x * y * y * y + std::sqrt(3444525.0) * x * y * z * z) + e_4 * (std::sqrt(2083725.0) * x * y);
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

        pc_69[k] = e_0 * (-std::sqrt(4.6142578125) * x * x * x * x * x * x * x * x * y * z - std::sqrt(41.5283203125) * x * x * x * x * x * x * y * y * y * z + std::sqrt(1386.328125) * x * x * x * x * x * x * y * z * z * z - std::sqrt(41.5283203125) * x * x * x * x * y * y * y * y * y * z + std::sqrt(5545.3125) * x * x * x * x * y * y * y * z * z * z - std::sqrt(6431.25) * x * x * x * x * y * z * z * z * z * z - std::sqrt(4.6142578125) * x * x * y * y * y * y * y * y * y * z + std::sqrt(1386.328125) * x * x * y * y * y * y * y * z * z * z - std::sqrt(6431.25) * x * x * y * y * y * z * z * z * z * z + std::sqrt(2100.0) * x * x * y * z * z * z * z * z * z * z) + e_1 * (std::sqrt(1038.2080078125) * x * x * x * x * x * x * y * z + std::sqrt(3880.5908203125) * x * x * x * x * y * y * y * z - std::sqrt(401.953125) * x * x * x * x * y * z * z * z + std::sqrt(779.8095703125) * x * x * y * y * y * y * y * z + std::sqrt(295.3125) * x * x * y * y * y * z * z * z + std::sqrt(57881.25) * x * x * y * z * z * z * z * z - std::sqrt(4.6142578125) * y * y * y * y * y * y * y * z + std::sqrt(1386.328125) * y * y * y * y * y * z * z * z - std::sqrt(6431.25) * y * y * y * z * z * z * z * z + std::sqrt(2100.0) * y * z * z * z * z * z * z * z) + e_2 * (std::sqrt(93041.89453125) * x * x * x * x * y * z + std::sqrt(124105.078125) * x * x * y * y * y * z + std::sqrt(1365525.0) * x * x * y * z * z * z + std::sqrt(2233.30078125) * y * y * y * y * y * z - std::sqrt(42525.0) * y * y * y * z * z * z + std::sqrt(231525.0) * y * z * z * z * z * z) + e_3 * (std::sqrt(4540725.0) * x * x * y * z + std::sqrt(4725.0) * y * y * y * z + std::sqrt(3194100.0) * y * z * z * z) + e_4 * (std::sqrt(3704400.0) * y * z);
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

        pc_70[k] = e_0 * (-std::sqrt(46.142578125) * x * x * x * x * x * x * x * y * z * z - std::sqrt(415.283203125) * x * x * x * x * x * y * y * y * z * z + std::sqrt(1312.5) * x * x * x * x * x * y * z * z * z * z - std::sqrt(415.283203125) * x * x * x * y * y * y * y * y * z * z + std::sqrt(5250.0) * x * x * x * y * y * y * z * z * z * z - std::sqrt(2218.125) * x * x * x * y * z * z * z * z * z * z - std::sqrt(46.142578125) * x * y * y * y * y * y * y * y * z * z + std::sqrt(1312.5) * x * y * y * y * y * y * z * z * z * z - std::sqrt(2218.125) * x * y * y * y * z * z * z * z * z * z + std::sqrt(210.0) * x * y * z * z * z * z * z * z * z * z) + e_1 * (-std::sqrt(46.142578125) * x * x * x * x * x * x * x * y - std::sqrt(415.283203125) * x * x * x * x * x * y * y * y + std::sqrt(184.5703125) * x * x * x * x * x * y * z * z - std::sqrt(415.283203125) * x * x * x * y * y * y * y * y + std::sqrt(738.28125) * x * x * x * y * y * y * z * z - std::sqrt(16078.125) * x * x * x * y * z * z * z * z - std::sqrt(46.142578125) * x * y * y * y * y * y * y * y + std::sqrt(184.5703125) * x * y * y * y * y * y * z * z - std::sqrt(16078.125) * x * y * y * y * z * z * z * z + std::sqrt(15172.5) * x * y * z * z * z * z * z * z) + e_2 * (-std::sqrt(9043.9453125) * x * x * x * x * x * y - std::sqrt(36175.78125) * x * x * x * y * y * y - std::sqrt(73828.125) * x * x * x * y * z * z - std::sqrt(9043.9453125) * x * y * y * y * y * y - std::sqrt(73828.125) * x * y * y * y * z * z + std::sqrt(295312.5) * x * y * z * z * z * z) + e_3 * (-std::sqrt(357328.125) * x * x * x * y - std::sqrt(357328.125) * x * y * y * y + std::sqrt(295312.5) * x * y * z * z) + e_4 * (-std::sqrt(578812.5) * x * y);
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

        pc_71[k] = e_0 * (std::sqrt(0.54931640625) * x * x * x * x * x * x * x * x * x * z + std::sqrt(8.7890625) * x * x * x * x * x * x * x * y * y * z - std::sqrt(205.322265625) * x * x * x * x * x * x * x * z * z * z + std::sqrt(19.775390625) * x * x * x * x * x * y * y * y * y * z - std::sqrt(1847.900390625) * x * x * x * x * x * y * y * z * z * z + std::sqrt(1265.625) * x * x * x * x * x * z * z * z * z * z + std::sqrt(8.7890625) * x * x * x * y * y * y * y * y * y * z - std::sqrt(1847.900390625) * x * x * x * y * y * y * y * z * z * z + std::sqrt(5062.5) * x * x * x * y * y * z * z * z * z * z - std::sqrt(680.625) * x * x * x * z * z * z * z * z * z * z + std::sqrt(0.54931640625) * x * y * y * y * y * y * y * y * y * z - std::sqrt(205.322265625) * x * y * y * y * y * y * y * z * z * z + std::sqrt(1265.625) * x * y * y * y * y * z * z * z * z * z - std::sqrt(680.625) * x * y * y * z * z * z * z * z * z * z + std::sqrt(10.0) * x * z * z * z * z * z * z * z * z * z) + e_1 * (-std::sqrt(177.978515625) * x * x * x * x * x * x * x * z - std::sqrt(1601.806640625) * x * x * x * x * x * y * y * z + std::sqrt(140.625) * x * x * x * x * x * z * z * z - std::sqrt(1601.806640625) * x * x * x * y * y * y * y * z + std::sqrt(562.5) * x * x * x * y * y * z * z * z - std::sqrt(14630.625) * x * x * x * z * z * z * z * z - std::sqrt(177.978515625) * x * y * y * y * y * y * y * z + std::sqrt(140.625) * x * y * y * y * y * z * z * z - std::sqrt(14630.625) * x * y * y * z * z * z * z * z + std::sqrt(90.0) * x * z * z * z * z * z * z * z) + e_2 * (-std::sqrt(20250.0) * x * x * x * x * x * z - std::sqrt(81000.0) * x * x * x * y * y * z - std::sqrt(284765.625) * x * x * x * z * z * z - std::sqrt(20250.0) * x * y * y * y * y * z - std::sqrt(284765.625) * x * y * y * z * z * z - std::sqrt(20250.0) * x * z * z * z * z * z) + e_3 * (-std::sqrt(1216265.625) * x * x * x * z - std::sqrt(1216265.625) * x * y * y * z - std::sqrt(1406250.0) * x * z * z * z) + e_4 * (-std::sqrt(3969000.0) * x * z);
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

        pc_72[k] = e_0 * (-std::sqrt(46.142578125) * x * x * x * x * x * x * x * x * z * z - std::sqrt(415.283203125) * x * x * x * x * x * x * y * y * z * z + std::sqrt(1312.5) * x * x * x * x * x * x * z * z * z * z - std::sqrt(415.283203125) * x * x * x * x * y * y * y * y * z * z + std::sqrt(5250.0) * x * x * x * x * y * y * z * z * z * z - std::sqrt(2218.125) * x * x * x * x * z * z * z * z * z * z - std::sqrt(46.142578125) * x * x * y * y * y * y * y * y * z * z + std::sqrt(1312.5) * x * x * y * y * y * y * z * z * z * z - std::sqrt(2218.125) * x * x * y * y * z * z * z * z * z * z + std::sqrt(210.0) * x * x * z * z * z * z * z * z * z * z) + e_1 * (-std::sqrt(46.142578125) * x * x * x * x * x * x * x * x - std::sqrt(415.283203125) * x * x * x * x * x * x * y * y + std::sqrt(46.142578125) * x * x * x * x * x * x * z * z - std::sqrt(415.283203125) * x * x * x * x * y * y * y * y + std::sqrt(46.142578125) * x * x * x * x * y * y * z * z - std::sqrt(8203.125) * x * x * x * x * z * z * z * z - std::sqrt(46.142578125) * x * x * y * y * y * y * y * y - std::sqrt(46.142578125) * x * x * y * y * y * y * z * z - std::sqrt(2953.125) * x * x * y * y * z * z * z * z + std::sqrt(5788.125) * x * x * z * z * z * z * z * z - std::sqrt(46.142578125) * y * y * y * y * y * y * z * z + std::sqrt(1312.5) * y * y * y * y * z * z * z * z - std::sqrt(2218.125) * y * y * z * z * z * z * z * z + std::sqrt(210.0) * z * z * z * z * z * z * z * z) + e_2 * (-std::sqrt(10382.080078125) * x * x * x * x * x * x - std::sqrt(44343.017578125) * x * x * x * x * y * y - std::sqrt(47250.0) * x * x * x * x * z * z - std::sqrt(13335.205078125) * x * x * y * y * y * y - std::sqrt(26578.125) * x * x * y * y * z * z + std::sqrt(73828.125) * x * x * z * z * z * z - std::sqrt(46.142578125) * y * y * y * y * y * y + std::sqrt(2953.125) * y * y * y * y * z * z - std::sqrt(73828.125) * y * y * z * z * z * z + std::sqrt(47250.0) * z * z * z * z * z * z) + e_3 * (-std::sqrt(425250.0) * x * x * x * x - std::sqrt(499078.125) * x * x * y * y + std::sqrt(2953.125) * x * x * z * z - std::sqrt(2953.125) * y * y * y * y - std::sqrt(239203.125) * y * y * z * z + std::sqrt(1181250.0) * z * z * z * z) + e_4 * (-std::sqrt(1302328.125) * x * x - std::sqrt(144703.125) * y * y + std::sqrt(2315250.0) * z * z);
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

        pc_73[k] = e_0 * (-std::sqrt(1.153564453125) * x * x * x * x * x * x * x * x * x * z - std::sqrt(4.6142578125) * x * x * x * x * x * x * x * y * y * z + std::sqrt(346.58203125) * x * x * x * x * x * x * x * z * z * z + std::sqrt(346.58203125) * x * x * x * x * x * y * y * z * z * z - std::sqrt(1607.8125) * x * x * x * x * x * z * z * z * z * z + std::sqrt(4.6142578125) * x * x * x * y * y * y * y * y * y * z - std::sqrt(346.58203125) * x * x * x * y * y * y * y * z * z * z + std::sqrt(525.0) * x * x * x * z * z * z * z * z * z * z + std::sqrt(1.153564453125) * x * y * y * y * y * y * y * y * y * z - std::sqrt(346.58203125) * x * y * y * y * y * y * y * z * z * z + std::sqrt(1607.8125) * x * y * y * y * y * z * z * z * z * z - std::sqrt(525.0) * x * y * y * z * z * z * z * z * z * z) + e_1 * (std::sqrt(226.0986328125) * x * x * x * x * x * x * x * z + std::sqrt(115.3564453125) * x * x * x * x * x * y * y * z + std::sqrt(73.828125) * x * x * x * x * x * z * z * z - std::sqrt(558.3251953125) * x * x * x * y * y * y * y * z + std::sqrt(5545.3125) * x * x * x * y * y * z * z * z + std::sqrt(6431.25) * x * x * x * z * z * z * z * z - std::sqrt(373.7548828125) * x * y * y * y * y * y * y * z + std::sqrt(4339.453125) * x * y * y * y * y * z * z * z - std::sqrt(57881.25) * x * y * y * z * z * z * z * z + std::sqrt(2100.0) * x * z * z * z * z * z * z * z) + e_2 * (std::sqrt(31026.26953125) * x * x * x * x * x * z + std::sqrt(8933.203125) * x * x * x * y * y * z + std::sqrt(231525.0) * x * x * x * z * z * z - std::sqrt(6662.98828125) * x * y * y * y * y * z - std::sqrt(798525.0) * x * y * y * z * z * z + std::sqrt(231525.0) * x * z * z * z * z * z) + e_3 * (std::sqrt(1209600.0) * x * x * x * z - std::sqrt(926100.0) * x * y * y * z + std::sqrt(3194100.0) * x * z * z * z) + e_4 * (std::sqrt(3704400.0) * x * z);
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

        pc_74[k] = e_0 * (std::sqrt(41.5283203125) * x * x * x * x * x * x * x * x * z * z - std::sqrt(41.5283203125) * x * x * x * x * x * x * y * y * z * z - std::sqrt(664.453125) * x * x * x * x * x * x * z * z * z * z - std::sqrt(1038.2080078125) * x * x * x * x * y * y * y * y * z * z + std::sqrt(2657.8125) * x * x * x * x * y * y * z * z * z * z + std::sqrt(525.0) * x * x * x * x * z * z * z * z * z * z - std::sqrt(373.7548828125) * x * x * y * y * y * y * y * y * z * z + std::sqrt(5980.078125) * x * x * y * y * y * y * z * z * z * z - std::sqrt(4725.0) * x * x * y * y * z * z * z * z * z * z) + e_1 * (std::sqrt(41.5283203125) * x * x * x * x * x * x * x * x - std::sqrt(41.5283203125) * x * x * x * x * x * x * y * y + std::sqrt(373.7548828125) * x * x * x * x * x * x * z * z - std::sqrt(1038.2080078125) * x * x * x * x * y * y * y * y + std::sqrt(373.7548828125) * x * x * x * x * y * y * z * z + std::sqrt(73.828125) * x * x * x * x * z * z * z * z - std::sqrt(373.7548828125) * x * x * y * y * y * y * y * y - std::sqrt(373.7548828125) * x * x * y * y * y * y * z * z - std::sqrt(66445.3125) * x * x * y * y * z * z * z * z + std::sqrt(4725.0) * x * x * z * z * z * z * z * z - std::sqrt(373.7548828125) * y * y * y * y * y * y * z * z + std::sqrt(5980.078125) * y * y * y * y * z * z * z * z - std::sqrt(4725.0) * y * y * z * z * z * z * z * z) + e_2 * (std::sqrt(9343.8720703125) * x * x * x * x * x * x - std::sqrt(18313.9892578125) * x * x * x * x * y * y + std::sqrt(32558.203125) * x * x * x * x * z * z - std::sqrt(63164.5751953125) * x * x * y * y * y * y - std::sqrt(598007.8125) * x * x * y * y * z * z + std::sqrt(170100.0) * x * x * z * z * z * z - std::sqrt(373.7548828125) * y * y * y * y * y * y + std::sqrt(5980.078125) * y * y * y * y * z * z - std::sqrt(170100.0) * y * y * z * z * z * z) + e_3 * (std::sqrt(248357.8125) * x * x * x * x - std::sqrt(1063125.0) * x * x * y * y + std::sqrt(861131.25) * x * x * z * z - std::sqrt(23920.3125) * y * y * y * y - std::sqrt(861131.25) * y * y * z * z) + e_4 * (std::sqrt(520931.25) * x * x - std::sqrt(520931.25) * y * y);
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

        pc_75[k] = e_0 * (std::sqrt(1.38427734375) * x * x * x * x * x * x * x * x * x * z - std::sqrt(22.1484375) * x * x * x * x * x * x * x * y * y * z - std::sqrt(177.802734375) * x * x * x * x * x * x * x * z * z * z - std::sqrt(138.427734375) * x * x * x * x * x * y * y * y * y * z + std::sqrt(4445.068359375) * x * x * x * x * x * y * y * z * z * z + std::sqrt(246.09375) * x * x * x * x * x * z * z * z * z * z - std::sqrt(22.1484375) * x * x * x * y * y * y * y * y * y * z + std::sqrt(4445.068359375) * x * x * x * y * y * y * y * z * z * z - std::sqrt(8859.375) * x * x * x * y * y * z * z * z * z * z + std::sqrt(1.38427734375) * x * y * y * y * y * y * y * y * y * z - std::sqrt(177.802734375) * x * y * y * y * y * y * y * z * z * z + std::sqrt(246.09375) * x * y * y * y * y * z * z * z * z * z) + e_1 * (-std::sqrt(5.537109375) * x * x * x * x * x * x * x * z + std::sqrt(935.771484375) * x * x * x * x * x * y * y * z - std::sqrt(3189.375) * x * x * x * x * x * z * z * z + std::sqrt(138.427734375) * x * x * x * y * y * y * y * z + std::sqrt(15750.0) * x * x * x * y * y * z * z * z + std::sqrt(3937.5) * x * x * x * z * z * z * z * z - std::sqrt(448.505859375) * x * y * y * y * y * y * y * z + std::sqrt(24609.375) * x * y * y * y * y * z * z * z - std::sqrt(35437.5) * x * y * y * z * z * z * z * z) + e_2 * (-std::sqrt(8859.375) * x * x * x * x * x * z + std::sqrt(141750.0) * x * x * x * y * y * z + std::sqrt(8859.375) * x * x * x * z * z * z + std::sqrt(8859.375) * x * y * y * y * y * z - std::sqrt(79734.375) * x * y * y * z * z * z) + e_3 * (-std::sqrt(8859.375) * x * x * x * z + std::sqrt(79734.375) * x * y * y * z);
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

        pc_76[k] = e_0 * (-std::sqrt(30.4541015625) * x * x * x * x * x * x * x * x * z * z + std::sqrt(2466.7822265625) * x * x * x * x * x * x * y * y * z * z + std::sqrt(54.140625) * x * x * x * x * x * x * z * z * z * z + std::sqrt(761.3525390625) * x * x * x * x * y * y * y * y * z * z - std::sqrt(5414.0625) * x * x * x * x * y * y * z * z * z * z - std::sqrt(761.3525390625) * x * x * y * y * y * y * y * y * z * z + std::sqrt(1353.515625) * x * x * y * y * y * y * z * z * z * z) + e_1 * (-std::sqrt(30.4541015625) * x * x * x * x * x * x * x * x + std::sqrt(2466.7822265625) * x * x * x * x * x * x * y * y - std::sqrt(3684.9462890625) * x * x * x * x * x * x * z * z + std::sqrt(761.3525390625) * x * x * x * x * y * y * y * y + std::sqrt(220030.8837890625) * x * x * x * x * y * y * z * z + std::sqrt(1353.515625) * x * x * x * x * z * z * z * z - std::sqrt(761.3525390625) * x * x * y * y * y * y * y * y - std::sqrt(761.3525390625) * x * x * y * y * y * y * z * z - std::sqrt(48726.5625) * x * x * y * y * z * z * z * z - std::sqrt(761.3525390625) * y * y * y * y * y * y * z * z + std::sqrt(1353.515625) * y * y * y * y * z * z * z * z) + e_2 * (-std::sqrt(6852.1728515625) * x * x * x * x * x * x + std::sqrt(475845.3369140625) * x * x * x * x * y * y - std::sqrt(12181.640625) * x * x * x * x * z * z - std::sqrt(19033.8134765625) * x * x * y * y * y * y + std::sqrt(438539.0625) * x * x * y * y * z * z - std::sqrt(761.3525390625) * y * y * y * y * y * y - std::sqrt(12181.640625) * y * y * y * y * z * z) + e_3 * (-std::sqrt(48726.5625) * x * x * x * x + std::sqrt(1754156.25) * x * x * y * y - std::sqrt(48726.5625) * y * y * y * y);

        pc_77[k] = e_0 * (-std::sqrt(2.537841796875) * x * x * x * x * x * x * x * x * x * z + std::sqrt(497.4169921875) * x * x * x * x * x * x * x * y * y * z + std::sqrt(4.51171875) * x * x * x * x * x * x * x * z * z * z - std::sqrt(1015.13671875) * x * x * x * x * x * y * y * z * z * z - std::sqrt(497.4169921875) * x * x * x * y * y * y * y * y * y * z + std::sqrt(1015.13671875) * x * x * x * y * y * y * y * z * z * z + std::sqrt(2.537841796875) * x * y * y * y * y * y * y * y * y * z - std::sqrt(4.51171875) * x * y * y * y * y * y * y * z * z * z) + e_1 * (-std::sqrt(822.2607421875) * x * x * x * x * x * x * x * z + std::sqrt(138962.0654296875) * x * x * x * x * x * y * y * z + std::sqrt(162.421875) * x * x * x * x * x * z * z * z - std::sqrt(57101.4404296875) * x * x * x * y * y * y * y * z - std::sqrt(16242.1875) * x * x * x * y * y * z * z * z - std::sqrt(822.2607421875) * x * y * y * y * y * y * y * z + std::sqrt(4060.546875) * x * y * y * y * y * z * z * z) + e_2 * (-std::sqrt(9136.23046875) * x * x * x * x * x * z + std::sqrt(913623.046875) * x * x * x * y * y * z - std::sqrt(228405.76171875) * x * y * y * y * y * z);
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

        pc_78[k] = e_0 * (-std::sqrt(5.07568359375) * x * x * x * x * x * x * x * x * x * y + std::sqrt(56.396484375) * x * x * x * x * x * x * x * y * y * y + std::sqrt(182.724609375) * x * x * x * x * x * x * x * y * z * z - std::sqrt(3431.162109375) * x * x * x * x * x * y * y * y * z * z - std::sqrt(56.396484375) * x * x * x * y * y * y * y * y * y * y + std::sqrt(3431.162109375) * x * x * x * y * y * y * y * y * z * z + std::sqrt(5.07568359375) * x * y * y * y * y * y * y * y * y * y - std::sqrt(182.724609375) * x * y * y * y * y * y * y * y * z * z) + e_1 * (-std::sqrt(2030.2734375) * x * x * x * x * x * x * x * y + std::sqrt(9826.5234375) * x * x * x * x * x * y * y * y + std::sqrt(11694.375) * x * x * x * x * x * y * z * z - std::sqrt(9826.5234375) * x * x * x * y * y * y * y * y + std::sqrt(2030.2734375) * x * y * y * y * y * y * y * y - std::sqrt(11694.375) * x * y * y * y * y * y * z * z) + e_2 * (-std::sqrt(73089.84375) * x * x * x * x * x * y + std::sqrt(292359.375) * x * x * x * y * z * z + std::sqrt(73089.84375) * x * y * y * y * y * y - std::sqrt(292359.375) * x * y * y * y * z * z) + e_3 * (-std::sqrt(519750.0) * x * x * x * y + std::sqrt(519750.0) * x * y * y * y);

        pc_79[k] = e_0 * (-std::sqrt(42.29736328125) * x * x * x * x * x * x * x * x * y * z + std::sqrt(169.189453125) * x * x * x * x * x * x * y * y * y * z + std::sqrt(1522.705078125) * x * x * x * x * x * x * y * z * z * z + std::sqrt(27.0703125) * x * x * x * x * y * y * y * y * y * z - std::sqrt(13704.345703125) * x * x * x * x * y * y * y * z * z * z - std::sqrt(169.189453125) * x * x * y * y * y * y * y * y * y * z + std::sqrt(7369.892578125) * x * x * y * y * y * y * y * z * z * z + std::sqrt(1.69189453125) * y * y * y * y * y * y * y * y * y * z - std::sqrt(60.908203125) * y * y * y * y * y * y * y * z * z * z) + e_1 * (-std::sqrt(676.7578125) * x * x * x * x * x * x * y * z - std::sqrt(10828.125) * x * x * x * x * y * y * y * z + std::sqrt(54817.3828125) * x * x * x * x * y * z * z * z + std::sqrt(243.6328125) * x * x * y * y * y * y * y * z + std::sqrt(24363.28125) * x * x * y * y * y * z * z * z + std::sqrt(108.28125) * y * y * y * y * y * y * y * z - std::sqrt(6090.8203125) * y * y * y * y * y * z * z * z) + e_2 * (std::sqrt(877078.125) * x * x * y * z * z * z - std::sqrt(97453.125) * y * y * y * z * z * z) + e_3 * (std::sqrt(877078.125) * x * x * y * z - std::sqrt(97453.125) * y * y * y * z);
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

        pc_80[k] = e_0 * (std::sqrt(1.23046875) * x * x * x * x * x * x * x * x * x * y - std::sqrt(315.0) * x * x * x * x * x * x * x * y * z * z - std::sqrt(4.921875) * x * x * x * x * x * y * y * y * y * y + std::sqrt(315.0) * x * x * x * x * x * y * y * y * z * z + std::sqrt(4429.6875) * x * x * x * x * x * y * z * z * z * z + std::sqrt(315.0) * x * x * x * y * y * y * y * y * z * z - std::sqrt(17718.75) * x * x * x * y * y * y * z * z * z * z + std::sqrt(1.23046875) * x * y * y * y * y * y * y * y * y * y - std::sqrt(315.0) * x * y * y * y * y * y * y * y * z * z + std::sqrt(4429.6875) * x * y * y * y * y * y * z * z * z * z) + e_1 * (std::sqrt(492.1875) * x * x * x * x * x * x * x * y - std::sqrt(19.6875) * x * x * x * x * x * y * y * y + std::sqrt(6378.75) * x * x * x * x * x * y * z * z - std::sqrt(19.6875) * x * x * x * y * y * y * y * y - std::sqrt(196875.0) * x * x * x * y * y * y * z * z + std::sqrt(70875.0) * x * x * x * y * z * z * z * z + std::sqrt(492.1875) * x * y * y * y * y * y * y * y + std::sqrt(6378.75) * x * y * y * y * y * y * z * z + std::sqrt(70875.0) * x * y * y * y * z * z * z * z) + e_2 * (std::sqrt(70875.0) * x * x * x * x * x * y - std::sqrt(70875.0) * x * x * x * y * y * y + std::sqrt(283500.0) * x * x * x * y * z * z + std::sqrt(70875.0) * x * y * y * y * y * y + std::sqrt(283500.0) * x * y * y * y * z * z + std::sqrt(637875.0) * x * y * z * z * z * z) + e_3 * (std::sqrt(637875.0) * x * x * x * y + std::sqrt(637875.0) * x * y * y * y + std::sqrt(7087500.0) * x * y * z * z) + e_4 * (std::sqrt(3472875.0) * x * y);
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

        pc_81[k] = e_0 * (std::sqrt(20.76416015625) * x * x * x * x * x * x * x * x * y * z + std::sqrt(9.228515625) * x * x * x * x * x * x * y * y * y * z - std::sqrt(1559.619140625) * x * x * x * x * x * x * y * z * z * z - std::sqrt(36.9140625) * x * x * x * x * y * y * y * y * y * z + std::sqrt(173.291015625) * x * x * x * x * y * y * y * z * z * z + std::sqrt(5315.625) * x * x * x * x * y * z * z * z * z * z - std::sqrt(9.228515625) * x * x * y * y * y * y * y * y * y * z + std::sqrt(1559.619140625) * x * x * y * y * y * y * y * z * z * z - std::sqrt(9450.0) * x * x * y * y * y * z * z * z * z * z + std::sqrt(2.30712890625) * y * y * y * y * y * y * y * y * y * z - std::sqrt(173.291015625) * y * y * y * y * y * y * y * z * z * z + std::sqrt(590.625) * y * y * y * y * y * z * z * z * z * z) + e_1 * (std::sqrt(332.2265625) * x * x * x * x * x * x * y * z + std::sqrt(590.625) * x * x * x * x * y * y * y * z + std::sqrt(31044.7265625) * x * x * x * x * y * z * z * z + std::sqrt(332.2265625) * x * x * y * y * y * y * y * z - std::sqrt(248210.15625) * x * x * y * y * y * z * z * z + std::sqrt(21262.5) * x * x * y * z * z * z * z * z + std::sqrt(147.65625) * y * y * y * y * y * y * y * z + std::sqrt(36.9140625) * y * y * y * y * y * z * z * z + std::sqrt(21262.5) * y * y * y * z * z * z * z * z) + e_2 * (std::sqrt(191362.5) * x * x * x * x * y * z - std::sqrt(340200.0) * x * x * y * y * y * z + std::sqrt(260465.625) * x * x * y * z * z * z + std::sqrt(21262.5) * y * y * y * y * y * z + std::sqrt(260465.625) * y * y * y * z * z * z + std::sqrt(85050.0) * y * z * z * z * z * z) + e_3 * (std::sqrt(643190.625) * x * x * y * z + std::sqrt(643190.625) * y * y * y * z + std::sqrt(2731050.0) * y * z * z * z) + e_4 * (std::sqrt(4167450.0) * y * z);
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

        pc_82[k] = e_0 * (-std::sqrt(0.25634765625) * x * x * x * x * x * x * x * x * x * y - std::sqrt(1.025390625) * x * x * x * x * x * x * x * y * y * y + std::sqrt(124.072265625) * x * x * x * x * x * x * x * y * z * z + std::sqrt(124.072265625) * x * x * x * x * x * y * y * y * z * z - std::sqrt(3215.625) * x * x * x * x * x * y * z * z * z * z + std::sqrt(1.025390625) * x * x * x * y * y * y * y * y * y * y - std::sqrt(124.072265625) * x * x * x * y * y * y * y * y * z * z + std::sqrt(2362.5) * x * x * x * y * z * z * z * z * z * z + std::sqrt(0.25634765625) * x * y * y * y * y * y * y * y * y * y - std::sqrt(124.072265625) * x * y * y * y * y * y * y * y * z * z + std::sqrt(3215.625) * x * y * y * y * y * y * z * z * z * z - std::sqrt(2362.5) * x * y * y * y * z * z * z * z * z * z) + e_1 * (-std::sqrt(102.5390625) * x * x * x * x * x * x * x * y - std::sqrt(102.5390625) * x * x * x * x * x * y * y * y - std::sqrt(5315.625) * x * x * x * x * x * y * z * z + std::sqrt(102.5390625) * x * x * x * y * y * y * y * y + std::sqrt(26250.0) * x * x * x * y * z * z * z * z + std::sqrt(102.5390625) * x * y * y * y * y * y * y * y + std::sqrt(5315.625) * x * y * y * y * y * y * z * z - std::sqrt(26250.0) * x * y * y * y * z * z * z * z) + e_2 * (-std::sqrt(24953.90625) * x * x * x * x * x * y + std::sqrt(14765.625) * x * x * x * y * z * z + std::sqrt(24953.90625) * x * y * y * y * y * y - std::sqrt(14765.625) * x * y * y * y * z * z) + e_3 * (-std::sqrt(236250.0) * x * x * x * y + std::sqrt(236250.0) * x * y * y * y);
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

        pc_83[k] = e_0 * (-std::sqrt(2.5634765625) * x * x * x * x * x * x * x * x * y * z - std::sqrt(10.25390625) * x * x * x * x * x * x * y * y * y * z + std::sqrt(256.34765625) * x * x * x * x * x * x * y * z * z * z + std::sqrt(256.34765625) * x * x * x * x * y * y * y * z * z * z - std::sqrt(1680.0) * x * x * x * x * y * z * z * z * z * z + std::sqrt(10.25390625) * x * x * y * y * y * y * y * y * y * z - std::sqrt(256.34765625) * x * x * y * y * y * y * y * z * z * z + std::sqrt(236.25) * x * x * y * z * z * z * z * z * z * z + std::sqrt(2.5634765625) * y * y * y * y * y * y * y * y * y * z - std::sqrt(256.34765625) * y * y * y * y * y * y * y * z * z * z + std::sqrt(1680.0) * y * y * y * y * y * z * z * z * z * z - std::sqrt(236.25) * y * y * y * z * z * z * z * z * z * z) + e_1 * (-std::sqrt(41.015625) * x * x * x * x * x * x * y * z - std::sqrt(14806.640625) * x * x * x * x * y * z * z * z + std::sqrt(369.140625) * x * x * y * y * y * y * y * z - std::sqrt(4101.5625) * x * x * y * y * y * z * z * z + std::sqrt(5906.25) * x * x * y * z * z * z * z * z + std::sqrt(164.0625) * y * y * y * y * y * y * y * z + std::sqrt(3322.265625) * y * y * y * y * y * z * z * z + std::sqrt(7586.25) * y * y * y * z * z * z * z * z - std::sqrt(945.0) * y * z * z * z * z * z * z * z) + e_2 * (-std::sqrt(53156.25) * x * x * x * x * y * z - std::sqrt(5906.25) * x * x * y * z * z * z + std::sqrt(53156.25) * y * y * y * y * y * z + std::sqrt(478406.25) * y * y * y * z * z * z - std::sqrt(23625.0) * y * z * z * z * z * z) + e_3 * (-std::sqrt(289406.25) * x * x * y * z + std::sqrt(2132156.25) * y * y * y * z + std::sqrt(23625.0) * y * z * z * z) + e_4 * (std::sqrt(1157625.0) * y * z);
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

        pc_84[k] = e_0 * (std::sqrt(0.030517578125) * x * x * x * x * x * x * x * x * x * x + std::sqrt(0.274658203125) * x * x * x * x * x * x * x * x * y * y - std::sqrt(17.578125) * x * x * x * x * x * x * x * x * z * z + std::sqrt(0.1220703125) * x * x * x * x * x * x * y * y * y * y - std::sqrt(70.3125) * x * x * x * x * x * x * y * y * z * z + std::sqrt(531.73828125) * x * x * x * x * x * x * z * z * z * z - std::sqrt(0.1220703125) * x * x * x * x * y * y * y * y * y * y + std::sqrt(531.73828125) * x * x * x * x * y * y * z * z * z * z - std::sqrt(661.25) * x * x * x * x * z * z * z * z * z * z - std::sqrt(0.274658203125) * x * x * y * y * y * y * y * y * y * y + std::sqrt(70.3125) * x * x * y * y * y * y * y * y * z * z - std::sqrt(531.73828125) * x * x * y * y * y * y * z * z * z * z + std::sqrt(11.25) * x * x * z * z * z * z * z * z * z * z - std::sqrt(0.030517578125) * y * y * y * y * y * y * y * y * y * y + std::sqrt(17.578125) * y * y * y * y * y * y * y * y * z * z - std::sqrt(531.73828125) * y * y * y * y * y * y * z * z * z * z + std::sqrt(661.25) * y * y * y * y * z * z * z * z * z * z - std::sqrt(11.25) * y * y * z * z * z * z * z * z * z * z) + e_1 * (std::sqrt(17.578125) * x * x * x * x * x * x * x * x + std::sqrt(70.3125) * x * x * x * x * x * x * y * y + std::sqrt(158.203125) * x * x * x * x * x * x * z * z + std::sqrt(158.203125) * x * x * x * x * y * y * z * z - std::sqrt(281.25) * x * x * x * x * z * z * z * z - std::sqrt(70.3125) * x * x * y * y * y * y * y * y - std::sqrt(158.203125) * x * x * y * y * y * y * z * z - std::sqrt(3645.0) * x * x * z * z * z * z * z * z - std::sqrt(17.578125) * y * y * y * y * y * y * y * y - std::sqrt(158.203125) * y * y * y * y * y * y * z * z + std::sqrt(281.25) * y * y * y * y * z * z * z * z + std::sqrt(3645.0) * y * y * z * z * z * z * z * z) + e_2 * (std::sqrt(4785.64453125) * x * x * x * x * x * x + std::sqrt(4785.64453125) * x * x * x * x * y * y + std::sqrt(2531.25) * x * x * x * x * z * z - std::sqrt(4785.64453125) * x * x * y * y * y * y - std::sqrt(253125.0) * x * x * z * z * z * z - std::sqrt(4785.64453125) * y * y * y * y * y * y - std::sqrt(2531.25) * y * y * y * y * z * z + std::sqrt(253125.0) * y * y * z * z * z * z) + e_3 * (std::sqrt(148781.25) * x * x * x * x - std::sqrt(820125.0) * x * x * z * z - std::sqrt(148781.25) * y * y * y * y + std::sqrt(820125.0) * y * y * z * z) + e_4 * (std::sqrt(124031.25) * x * x - std::sqrt(124031.25) * y * y);
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

        pc_85[k] = e_0 * (-std::sqrt(2.5634765625) * x * x * x * x * x * x * x * x * x * z - std::sqrt(10.25390625) * x * x * x * x * x * x * x * y * y * z + std::sqrt(256.34765625) * x * x * x * x * x * x * x * z * z * z + std::sqrt(256.34765625) * x * x * x * x * x * y * y * z * z * z - std::sqrt(1680.0) * x * x * x * x * x * z * z * z * z * z + std::sqrt(10.25390625) * x * x * x * y * y * y * y * y * y * z - std::sqrt(256.34765625) * x * x * x * y * y * y * y * z * z * z + std::sqrt(236.25) * x * x * x * z * z * z * z * z * z * z + std::sqrt(2.5634765625) * x * y * y * y * y * y * y * y * y * z - std::sqrt(256.34765625) * x * y * y * y * y * y * y * z * z * z + std::sqrt(1680.0) * x * y * y * y * y * z * z * z * z * z - std::sqrt(236.25) * x * y * y * z * z * z * z * z * z * z) + e_1 * (-std::sqrt(164.0625) * x * x * x * x * x * x * x * z - std::sqrt(369.140625) * x * x * x * x * x * y * y * z - std::sqrt(3322.265625) * x * x * x * x * x * z * z * z + std::sqrt(4101.5625) * x * x * x * y * y * z * z * z - std::sqrt(7586.25) * x * x * x * z * z * z * z * z + std::sqrt(41.015625) * x * y * y * y * y * y * y * z + std::sqrt(14806.640625) * x * y * y * y * y * z * z * z - std::sqrt(5906.25) * x * y * y * z * z * z * z * z + std::sqrt(945.0) * x * z * z * z * z * z * z * z) + e_2 * (-std::sqrt(53156.25) * x * x * x * x * x * z - std::sqrt(478406.25) * x * x * x * z * z * z + std::sqrt(53156.25) * x * y * y * y * y * z + std::sqrt(5906.25) * x * y * y * z * z * z + std::sqrt(23625.0) * x * z * z * z * z * z) + e_3 * (-std::sqrt(2132156.25) * x * x * x * z + std::sqrt(289406.25) * x * y * y * z - std::sqrt(23625.0) * x * z * z * z) + e_4 * (-std::sqrt(1157625.0) * x * z);
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

        pc_86[k] = e_0 * (-std::sqrt(0.0640869140625) * x * x * x * x * x * x * x * x * x * x - std::sqrt(0.0640869140625) * x * x * x * x * x * x * x * x * y * y + std::sqrt(31.01806640625) * x * x * x * x * x * x * x * x * z * z + std::sqrt(0.25634765625) * x * x * x * x * x * x * y * y * y * y - std::sqrt(803.90625) * x * x * x * x * x * x * z * z * z * z + std::sqrt(0.25634765625) * x * x * x * x * y * y * y * y * y * y - std::sqrt(124.072265625) * x * x * x * x * y * y * y * y * z * z + std::sqrt(803.90625) * x * x * x * x * y * y * z * z * z * z + std::sqrt(590.625) * x * x * x * x * z * z * z * z * z * z - std::sqrt(0.0640869140625) * x * x * y * y * y * y * y * y * y * y + std::sqrt(803.90625) * x * x * y * y * y * y * z * z * z * z - std::sqrt(2362.5) * x * x * y * y * z * z * z * z * z * z - std::sqrt(0.0640869140625) * y * y * y * y * y * y * y * y * y * y + std::sqrt(31.01806640625) * y * y * y * y * y * y * y * y * z * z - std::sqrt(803.90625) * y * y * y * y * y * y * z * z * z * z + std::sqrt(590.625) * y * y * y * y * z * z * z * z * z * z) + e_1 * (-std::sqrt(36.9140625) * x * x * x * x * x * x * x * x - std::sqrt(16.40625) * x * x * x * x * x * x * y * y - std::sqrt(200.9765625) * x * x * x * x * x * x * z * z + std::sqrt(16.40625) * x * x * x * x * y * y * y * y + std::sqrt(10668.1640625) * x * x * x * x * y * y * z * z - std::sqrt(1050.0) * x * x * x * x * z * z * z * z - std::sqrt(16.40625) * x * x * y * y * y * y * y * y + std::sqrt(10668.1640625) * x * x * y * y * y * y * z * z - std::sqrt(151200.0) * x * x * y * y * z * z * z * z + std::sqrt(9450.0) * x * x * z * z * z * z * z * z - std::sqrt(36.9140625) * y * y * y * y * y * y * y * y - std::sqrt(200.9765625) * y * y * y * y * y * y * z * z - std::sqrt(1050.0) * y * y * y * y * z * z * z * z + std::sqrt(9450.0) * y * y * z * z * z * z * z * z) + e_2 * (-std::sqrt(8868.603515625) * x * x * x * x * x * x + std::sqrt(1116.650390625) * x * x * x * x * y * y - std::sqrt(23071.2890625) * x * x * x * x * z * z + std::sqrt(1116.650390625) * x * x * y * y * y * y - std::sqrt(299003.90625) * x * x * y * y * z * z + std::sqrt(191362.5) * x * x * z * z * z * z - std::sqrt(8868.603515625) * y * y * y * y * y * y - std::sqrt(23071.2890625) * y * y * y * y * z * z + std::sqrt(191362.5) * y * y * z * z * z * z + std::sqrt(9450.0) * z * z * z * z * z * z) + e_3 * (-std::sqrt(260465.625) * x * x * x * x - std::sqrt(2362.5) * x * x * y * y + std::sqrt(151200.0) * x * x * z * z - std::sqrt(260465.625) * y * y * y * y + std::sqrt(151200.0) * y * y * z * z + std::sqrt(604800.0) * z * z * z * z) + e_4 * (-std::sqrt(463050.0) * x * x - std::sqrt(463050.0) * y * y + std::sqrt(1852200.0) * z * z);
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

        pc_87[k] = e_0 * (std::sqrt(2.30712890625) * x * x * x * x * x * x * x * x * x * z - std::sqrt(9.228515625) * x * x * x * x * x * x * x * y * y * z - std::sqrt(173.291015625) * x * x * x * x * x * x * x * z * z * z - std::sqrt(36.9140625) * x * x * x * x * x * y * y * y * y * z + std::sqrt(1559.619140625) * x * x * x * x * x * y * y * z * z * z + std::sqrt(590.625) * x * x * x * x * x * z * z * z * z * z + std::sqrt(9.228515625) * x * x * x * y * y * y * y * y * y * z + std::sqrt(173.291015625) * x * x * x * y * y * y * y * z * z * z - std::sqrt(9450.0) * x * x * x * y * y * z * z * z * z * z + std::sqrt(20.76416015625) * x * y * y * y * y * y * y * y * y * z - std::sqrt(1559.619140625) * x * y * y * y * y * y * y * z * z * z + std::sqrt(5315.625) * x * y * y * y * y * z * z * z * z * z) + e_1 * (std::sqrt(147.65625) * x * x * x * x * x * x * x * z + std::sqrt(332.2265625) * x * x * x * x * x * y * y * z + std::sqrt(36.9140625) * x * x * x * x * x * z * z * z + std::sqrt(590.625) * x * x * x * y * y * y * y * z - std::sqrt(248210.15625) * x * x * x * y * y * z * z * z + std::sqrt(21262.5) * x * x * x * z * z * z * z * z + std::sqrt(332.2265625) * x * y * y * y * y * y * y * z + std::sqrt(31044.7265625) * x * y * y * y * y * z * z * z + std::sqrt(21262.5) * x * y * y * z * z * z * z * z) + e_2 * (std::sqrt(21262.5) * x * x * x * x * x * z - std::sqrt(340200.0) * x * x * x * y * y * z + std::sqrt(260465.625) * x * x * x * z * z * z + std::sqrt(191362.5) * x * y * y * y * y * z + std::sqrt(260465.625) * x * y * y * z * z * z + std::sqrt(85050.0) * x * z * z * z * z * z) + e_3 * (std::sqrt(643190.625) * x * x * x * z + std::sqrt(643190.625) * x * y * y * z + std::sqrt(2731050.0) * x * z * z * z) + e_4 * (std::sqrt(4167450.0) * x * z);
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

        pc_88[k] = e_0 * (std::sqrt(0.076904296875) * x * x * x * x * x * x * x * x * x * x - std::sqrt(1.922607421875) * x * x * x * x * x * x * x * x * y * y - std::sqrt(19.6875) * x * x * x * x * x * x * x * x * z * z - std::sqrt(2.7685546875) * x * x * x * x * x * x * y * y * y * y + std::sqrt(708.75) * x * x * x * x * x * x * y * y * z * z + std::sqrt(276.85546875) * x * x * x * x * x * x * z * z * z * z + std::sqrt(2.7685546875) * x * x * x * x * y * y * y * y * y * y - std::sqrt(13565.91796875) * x * x * x * x * y * y * z * z * z * z + std::sqrt(1.922607421875) * x * x * y * y * y * y * y * y * y * y - std::sqrt(708.75) * x * x * y * y * y * y * y * y * z * z + std::sqrt(13565.91796875) * x * x * y * y * y * y * z * z * z * z - std::sqrt(0.076904296875) * y * y * y * y * y * y * y * y * y * y + std::sqrt(19.6875) * y * y * y * y * y * y * y * y * z * z - std::sqrt(276.85546875) * y * y * y * y * y * y * z * z * z * z) + e_1 * (std::sqrt(44.296875) * x * x * x * x * x * x * x * x - std::sqrt(492.1875) * x * x * x * x * x * x * y * y + std::sqrt(4.921875) * x * x * x * x * x * x * z * z - std::sqrt(89701.171875) * x * x * x * x * y * y * z * z + std::sqrt(17718.75) * x * x * x * x * z * z * z * z + std::sqrt(492.1875) * x * x * y * y * y * y * y * y + std::sqrt(89701.171875) * x * x * y * y * y * y * z * z - std::sqrt(44.296875) * y * y * y * y * y * y * y * y - std::sqrt(4.921875) * y * y * y * y * y * y * z * z - std::sqrt(17718.75) * y * y * y * y * z * z * z * z) + e_2 * (std::sqrt(6921.38671875) * x * x * x * x * x * x - std::sqrt(99944.82421875) * x * x * x * x * y * y + std::sqrt(70875.0) * x * x * x * x * z * z + std::sqrt(99944.82421875) * x * x * y * y * y * y + std::sqrt(159468.75) * x * x * z * z * z * z - std::sqrt(6921.38671875) * y * y * y * y * y * y - std::sqrt(70875.0) * y * y * y * y * z * z - std::sqrt(159468.75) * y * y * z * z * z * z) + e_3 * (std::sqrt(159468.75) * x * x * x * x + std::sqrt(1771875.0) * x * x * z * z - std::sqrt(159468.75) * y * y * y * y - std::sqrt(1771875.0) * y * y * z * z) + e_4 * (std::sqrt(868218.75) * x * x - std::sqrt(868218.75) * y * y);
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

        pc_89[k] = e_0 * (-std::sqrt(1.69189453125) * x * x * x * x * x * x * x * x * x * z + std::sqrt(169.189453125) * x * x * x * x * x * x * x * y * y * z + std::sqrt(60.908203125) * x * x * x * x * x * x * x * z * z * z - std::sqrt(27.0703125) * x * x * x * x * x * y * y * y * y * z - std::sqrt(7369.892578125) * x * x * x * x * x * y * y * z * z * z - std::sqrt(169.189453125) * x * x * x * y * y * y * y * y * y * z + std::sqrt(13704.345703125) * x * x * x * y * y * y * y * z * z * z + std::sqrt(42.29736328125) * x * y * y * y * y * y * y * y * y * z - std::sqrt(1522.705078125) * x * y * y * y * y * y * y * z * z * z) + e_1 * (-std::sqrt(108.28125) * x * x * x * x * x * x * x * z - std::sqrt(243.6328125) * x * x * x * x * x * y * y * z + std::sqrt(6090.8203125) * x * x * x * x * x * z * z * z + std::sqrt(10828.125) * x * x * x * y * y * y * y * z - std::sqrt(24363.28125) * x * x * x * y * y * z * z * z + std::sqrt(676.7578125) * x * y * y * y * y * y * y * z - std::sqrt(54817.3828125) * x * y * y * y * y * z * z * z) + e_2 * (std::sqrt(97453.125) * x * x * x * z * z * z - std::sqrt(877078.125) * x * y * y * z * z * z) + e_3 * (std::sqrt(97453.125) * x * x * x * z - std::sqrt(877078.125) * x * y * y * z);
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

        pc_90[k] = e_0 * (-std::sqrt(0.1409912109375) * x * x * x * x * x * x * x * x * x * x + std::sqrt(31.7230224609375) * x * x * x * x * x * x * x * x * y * y + std::sqrt(5.07568359375) * x * x * x * x * x * x * x * x * z * z - std::sqrt(27.63427734375) * x * x * x * x * x * x * y * y * y * y - std::sqrt(1299.375) * x * x * x * x * x * x * y * y * z * z - std::sqrt(27.63427734375) * x * x * x * x * y * y * y * y * y * y + std::sqrt(4568.115234375) * x * x * x * x * y * y * y * y * z * z + std::sqrt(31.7230224609375) * x * x * y * y * y * y * y * y * y * y - std::sqrt(1299.375) * x * x * y * y * y * y * y * y * z * z - std::sqrt(0.1409912109375) * y * y * y * y * y * y * y * y * y * y + std::sqrt(5.07568359375) * y * y * y * y * y * y * y * y * z * z) + e_1 * (-std::sqrt(81.2109375) * x * x * x * x * x * x * x * x + std::sqrt(8121.09375) * x * x * x * x * x * x * y * y + std::sqrt(730.8984375) * x * x * x * x * x * x * z * z - std::sqrt(8121.09375) * x * x * x * x * y * y * y * y - std::sqrt(18272.4609375) * x * x * x * x * y * y * z * z + std::sqrt(8121.09375) * x * x * y * y * y * y * y * y - std::sqrt(18272.4609375) * x * x * y * y * y * y * z * z - std::sqrt(81.2109375) * y * y * y * y * y * y * y * y + std::sqrt(730.8984375) * y * y * y * y * y * y * z * z) + e_2 * (-std::sqrt(4568.115234375) * x * x * x * x * x * x + std::sqrt(114202.880859375) * x * x * x * x * y * y + std::sqrt(18272.4609375) * x * x * x * x * z * z + std::sqrt(114202.880859375) * x * x * y * y * y * y - std::sqrt(657808.59375) * x * x * y * y * z * z - std::sqrt(4568.115234375) * y * y * y * y * y * y + std::sqrt(18272.4609375) * y * y * y * y * z * z) + e_3 * (-std::sqrt(32484.375) * x * x * x * x + std::sqrt(1169437.5) * x * x * y * y - std::sqrt(32484.375) * y * y * y * y);

        pc_91[k] = e_0 * (std::sqrt(71.0595703125) * x * x * x * x * x * x * x * x * y * z - std::sqrt(2850.2783203125) * x * x * x * x * x * x * y * y * y * z + std::sqrt(8598.2080078125) * x * x * x * x * y * y * y * y * y * z - std::sqrt(639.5361328125) * x * x * y * y * y * y * y * y * y * z) + e_1 * (std::sqrt(5755.8251953125) * x * x * x * x * x * x * y * z + std::sqrt(15988.4033203125) * x * x * x * x * y * y * y * z + std::sqrt(639.5361328125) * x * x * y * y * y * y * y * z - std::sqrt(639.5361328125) * y * y * y * y * y * y * y * z) + e_2 * (std::sqrt(575582.51953125) * x * x * x * x * y * z + std::sqrt(255814.453125) * x * x * y * y * y * z - std::sqrt(63953.61328125) * y * y * y * y * y * z) + e_3 * (std::sqrt(4093031.25) * x * x * y * z - std::sqrt(454781.25) * y * y * y * z);
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

        pc_92[k] = e_0 * (std::sqrt(592.1630859375) * x * x * x * x * x * x * x * y * z * z - std::sqrt(14804.0771484375) * x * x * x * x * x * y * y * y * z * z + std::sqrt(22762.7490234375) * x * x * x * y * y * y * y * y * z * z - std::sqrt(213.1787109375) * x * y * y * y * y * y * y * y * z * z) + e_1 * (std::sqrt(592.1630859375) * x * x * x * x * x * x * x * y - std::sqrt(14804.0771484375) * x * x * x * x * x * y * y * y + std::sqrt(21317.87109375) * x * x * x * x * x * y * z * z + std::sqrt(22762.7490234375) * x * x * x * y * y * y * y * y + std::sqrt(85271.484375) * x * x * x * y * y * y * z * z - std::sqrt(213.1787109375) * x * y * y * y * y * y * y * y + std::sqrt(21317.87109375) * x * y * y * y * y * y * z * z) + e_2 * (std::sqrt(21317.87109375) * x * x * x * x * x * y + std::sqrt(85271.484375) * x * x * x * y * y * y + std::sqrt(1364343.75) * x * x * x * y * z * z + std::sqrt(21317.87109375) * x * y * y * y * y * y + std::sqrt(1364343.75) * x * y * y * y * z * z) + e_3 * (std::sqrt(1364343.75) * x * x * x * y + std::sqrt(1364343.75) * x * y * y * y + std::sqrt(5457375.0) * x * y * z * z) + e_4 * (std::sqrt(5457375.0) * x * y);

        pc_93[k] = e_0 * (-std::sqrt(17.2265625) * x * x * x * x * x * x * x * x * y * z + std::sqrt(155.0390625) * x * x * x * x * x * x * y * y * y * z + std::sqrt(1722.65625) * x * x * x * x * x * x * y * z * z * z + std::sqrt(17.2265625) * x * x * x * x * y * y * y * y * y * z - std::sqrt(27562.5) * x * x * x * x * y * y * y * z * z * z - std::sqrt(155.0390625) * x * x * y * y * y * y * y * y * y * z + std::sqrt(15503.90625) * x * x * y * y * y * y * y * z * z * z) + e_1 * (std::sqrt(2084.4140625) * x * x * x * x * x * x * y * z - std::sqrt(72782.2265625) * x * x * x * x * y * y * y * z + std::sqrt(15503.90625) * x * x * x * x * y * z * z * z + std::sqrt(18759.7265625) * x * x * y * y * y * y * y * z + std::sqrt(62015.625) * x * x * y * y * y * z * z * z - std::sqrt(155.0390625) * y * y * y * y * y * y * y * z + std::sqrt(15503.90625) * y * y * y * y * y * z * z * z) + e_2 * (std::sqrt(15503.90625) * x * x * x * x * y * z + std::sqrt(62015.625) * x * x * y * y * y * z + std::sqrt(558140.625) * x * x * y * z * z * z + std::sqrt(15503.90625) * y * y * y * y * y * z + std::sqrt(558140.625) * y * y * y * z * z * z) + e_3 * (std::sqrt(1550390.625) * x * x * y * z + std::sqrt(1550390.625) * y * y * y * z + std::sqrt(992250.0) * y * z * z * z) + e_4 * (std::sqrt(3969000.0) * y * z);
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

        pc_94[k] = e_0 * (-std::sqrt(290.6982421875) * x * x * x * x * x * x * x * y * z * z + std::sqrt(1582.6904296875) * x * x * x * x * x * y * y * y * z * z + std::sqrt(2067.1875) * x * x * x * x * x * y * z * z * z * z + std::sqrt(1582.6904296875) * x * x * x * y * y * y * y * y * z * z - std::sqrt(22968.75) * x * x * x * y * y * y * z * z * z * z - std::sqrt(290.6982421875) * x * y * y * y * y * y * y * y * z * z + std::sqrt(2067.1875) * x * y * y * y * y * y * z * z * z * z) + e_1 * (-std::sqrt(290.6982421875) * x * x * x * x * x * x * x * y + std::sqrt(1582.6904296875) * x * x * x * x * x * y * y * y + std::sqrt(1162.79296875) * x * x * x * x * x * y * z * z + std::sqrt(1582.6904296875) * x * x * x * y * y * y * y * y - std::sqrt(12919.921875) * x * x * x * y * y * y * z * z - std::sqrt(290.6982421875) * x * y * y * y * y * y * y * y + std::sqrt(1162.79296875) * x * y * y * y * y * y * z * z) + e_2 * (-std::sqrt(10465.13671875) * x * x * x * x * x * y + std::sqrt(116279.296875) * x * x * x * y * y * y - std::sqrt(10465.13671875) * x * y * y * y * y * y);

        pc_95[k] = e_0 * (std::sqrt(3.5888671875) * x * x * x * x * x * x * x * x * y * z - std::sqrt(3.5888671875) * x * x * x * x * x * x * y * y * y * z - std::sqrt(918.75) * x * x * x * x * x * x * y * z * z * z - std::sqrt(89.7216796875) * x * x * x * x * y * y * y * y * y * z + std::sqrt(3675.0) * x * x * x * x * y * y * y * z * z * z + std::sqrt(918.75) * x * x * x * x * y * z * z * z * z * z - std::sqrt(32.2998046875) * x * x * y * y * y * y * y * y * y * z + std::sqrt(8268.75) * x * x * y * y * y * y * y * z * z * z - std::sqrt(8268.75) * x * x * y * y * y * z * z * z * z * z) + e_1 * (-std::sqrt(1898.5107421875) * x * x * x * x * x * x * y * z + std::sqrt(3448.9013671875) * x * x * x * x * y * y * y * z + std::sqrt(918.75) * x * x * x * x * y * z * z * z + std::sqrt(9334.6435546875) * x * x * y * y * y * y * y * z + std::sqrt(132300.0) * x * x * y * y * y * z * z * z - std::sqrt(8268.75) * x * x * y * z * z * z * z * z - std::sqrt(32.2998046875) * y * y * y * y * y * y * y * z + std::sqrt(8268.75) * y * y * y * y * y * z * z * z - std::sqrt(8268.75) * y * y * y * z * z * z * z * z) + e_2 * (-std::sqrt(37338.57421875) * x * x * x * x * y * z + std::sqrt(1451682.421875) * x * x * y * y * y * z + std::sqrt(33075.0) * x * x * y * z * z * z + std::sqrt(15633.10546875) * y * y * y * y * y * z + std::sqrt(33075.0) * y * y * y * z * z * z - std::sqrt(33075.0) * y * z * z * z * z * z) + e_3 * (std::sqrt(1000518.75) * x * x * y * z + std::sqrt(1000518.75) * y * y * y * z - std::sqrt(132300.0) * y * z * z * z) + e_4 * (std::sqrt(529200.0) * y * z);
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

        pc_96[k] = e_0 * (std::sqrt(35.888671875) * x * x * x * x * x * x * x * y * z * z - std::sqrt(35.888671875) * x * x * x * x * x * y * y * y * z * z - std::sqrt(574.21875) * x * x * x * x * x * y * z * z * z * z - std::sqrt(897.216796875) * x * x * x * y * y * y * y * y * z * z + std::sqrt(2296.875) * x * x * x * y * y * y * z * z * z * z + std::sqrt(91.875) * x * x * x * y * z * z * z * z * z * z - std::sqrt(322.998046875) * x * y * y * y * y * y * y * y * z * z + std::sqrt(5167.96875) * x * y * y * y * y * y * z * z * z * z - std::sqrt(826.875) * x * y * y * y * z * z * z * z * z * z) + e_1 * (std::sqrt(35.888671875) * x * x * x * x * x * x * x * y - std::sqrt(35.888671875) * x * x * x * x * x * y * y * y - std::sqrt(1291.9921875) * x * x * x * x * x * y * z * z - std::sqrt(897.216796875) * x * x * x * y * y * y * y * y - std::sqrt(5167.96875) * x * x * x * y * y * y * z * z + std::sqrt(2296.875) * x * x * x * y * z * z * z * z - std::sqrt(322.998046875) * x * y * y * y * y * y * y * y - std::sqrt(1291.9921875) * x * y * y * y * y * y * z * z + std::sqrt(186046.875) * x * y * y * y * z * z * z * z - std::sqrt(3307.5) * x * y * z * z * z * z * z * z) + e_2 * (std::sqrt(1291.9921875) * x * x * x * x * x * y - std::sqrt(46511.71875) * x * x * x * y * y * y - std::sqrt(20671.875) * x * x * x * y * z * z - std::sqrt(63307.6171875) * x * y * y * y * y * y + std::sqrt(1012921.875) * x * y * y * y * z * z + std::sqrt(82687.5) * x * y * z * z * z * z) + e_3 * (-std::sqrt(20671.875) * x * x * x * y - std::sqrt(516796.875) * x * y * y * y + std::sqrt(2067187.5) * x * y * z * z) + e_4 * (-std::sqrt(82687.5) * x * y);
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

        pc_97[k] = e_0 * (-std::sqrt(0.42724609375) * x * x * x * x * x * x * x * x * x * z + std::sqrt(138.427734375) * x * x * x * x * x * x * x * z * z * z + std::sqrt(15.380859375) * x * x * x * x * x * y * y * y * y * z - std::sqrt(138.427734375) * x * x * x * x * x * y * y * z * z * z - std::sqrt(246.09375) * x * x * x * x * x * z * z * z * z * z + std::sqrt(27.34375) * x * x * x * y * y * y * y * y * y * z - std::sqrt(3460.693359375) * x * x * x * y * y * y * y * z * z * z + std::sqrt(984.375) * x * x * x * y * y * z * z * z * z * z + std::sqrt(4.375) * x * x * x * z * z * z * z * z * z * z + std::sqrt(3.84521484375) * x * y * y * y * y * y * y * y * y * z - std::sqrt(1245.849609375) * x * y * y * y * y * y * y * z * z * z + std::sqrt(2214.84375) * x * y * y * y * y * z * z * z * z * z - std::sqrt(39.375) * x * y * y * z * z * z * z * z * z * z) + e_1 * (std::sqrt(138.427734375) * x * x * x * x * x * x * x * z - std::sqrt(138.427734375) * x * x * x * x * x * y * y * z + std::sqrt(6152.34375) * x * x * x * x * x * z * z * z - std::sqrt(3460.693359375) * x * x * x * y * y * y * y * z - std::sqrt(24609.375) * x * x * x * y * y * z * z * z - std::sqrt(6654.375) * x * x * x * z * z * z * z * z - std::sqrt(1245.849609375) * x * y * y * y * y * y * y * z - std::sqrt(55371.09375) * x * y * y * y * y * z * z * z + std::sqrt(59889.375) * x * y * y * z * z * z * z * z) + e_2 * (std::sqrt(55371.09375) * x * x * x * x * x * z - std::sqrt(221484.375) * x * x * x * y * y * z - std::sqrt(8859.375) * x * x * x * z * z * z - std::sqrt(498339.84375) * x * y * y * y * y * z + std::sqrt(79734.375) * x * y * y * z * z * z) + e_3 * (std::sqrt(284484.375) * x * x * x * z - std::sqrt(2560359.375) * x * y * y * z);
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

        pc_98[k] = e_0 * (std::sqrt(35.888671875) * x * x * x * x * x * x * x * x * z * z - std::sqrt(35.888671875) * x * x * x * x * x * x * y * y * z * z - std::sqrt(574.21875) * x * x * x * x * x * x * z * z * z * z - std::sqrt(897.216796875) * x * x * x * x * y * y * y * y * z * z + std::sqrt(2296.875) * x * x * x * x * y * y * z * z * z * z + std::sqrt(91.875) * x * x * x * x * z * z * z * z * z * z - std::sqrt(322.998046875) * x * x * y * y * y * y * y * y * z * z + std::sqrt(5167.96875) * x * x * y * y * y * y * z * z * z * z - std::sqrt(826.875) * x * x * y * y * z * z * z * z * z * z) + e_1 * (std::sqrt(35.888671875) * x * x * x * x * x * x * x * x - std::sqrt(35.888671875) * x * x * x * x * x * x * y * y + std::sqrt(322.998046875) * x * x * x * x * x * x * z * z - std::sqrt(897.216796875) * x * x * x * x * y * y * y * y + std::sqrt(322.998046875) * x * x * x * x * y * y * z * z - std::sqrt(28136.71875) * x * x * x * x * z * z * z * z - std::sqrt(322.998046875) * x * x * y * y * y * y * y * y - std::sqrt(322.998046875) * x * x * y * y * y * y * z * z + std::sqrt(82687.5) * x * x * y * y * z * z * z * z + std::sqrt(826.875) * x * x * z * z * z * z * z * z - std::sqrt(322.998046875) * y * y * y * y * y * y * z * z + std::sqrt(5167.96875) * y * y * y * y * z * z * z * z - std::sqrt(826.875) * y * y * z * z * z * z * z * z) + e_2 * (std::sqrt(8074.951171875) * x * x * x * x * x * x - std::sqrt(15826.904296875) * x * x * x * x * y * y - std::sqrt(129199.21875) * x * x * x * x * z * z - std::sqrt(54586.669921875) * x * x * y * y * y * y + std::sqrt(744187.5) * x * x * y * y * z * z - std::sqrt(20671.875) * x * x * z * z * z * z - std::sqrt(322.998046875) * y * y * y * y * y * y + std::sqrt(5167.96875) * y * y * y * y * z * z + std::sqrt(20671.875) * y * y * z * z * z * z) + e_3 * (std::sqrt(82687.5) * x * x * x * x - std::sqrt(186046.875) * x * x * y * y - std::sqrt(516796.875) * x * x * z * z - std::sqrt(20671.875) * y * y * y * y + std::sqrt(516796.875) * y * y * z * z) + e_4 * (std::sqrt(20671.875) * x * x - std::sqrt(20671.875) * y * y);
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

        pc_99[k] = e_0 * (std::sqrt(0.897216796875) * x * x * x * x * x * x * x * x * x * z - std::sqrt(3.5888671875) * x * x * x * x * x * x * x * y * y * z - std::sqrt(229.6875) * x * x * x * x * x * x * x * z * z * z - std::sqrt(14.35546875) * x * x * x * x * x * y * y * y * y * z + std::sqrt(2067.1875) * x * x * x * x * x * y * y * z * z * z + std::sqrt(229.6875) * x * x * x * x * x * z * z * z * z * z + std::sqrt(3.5888671875) * x * x * x * y * y * y * y * y * y * z + std::sqrt(229.6875) * x * x * x * y * y * y * y * z * z * z - std::sqrt(3675.0) * x * x * x * y * y * z * z * z * z * z + std::sqrt(8.074951171875) * x * y * y * y * y * y * y * y * y * z - std::sqrt(2067.1875) * x * y * y * y * y * y * y * z * z * z + std::sqrt(2067.1875) * x * y * y * y * y * z * z * z * z * z) + e_1 * (-std::sqrt(175.8544921875) * x * x * x * x * x * x * x * z + std::sqrt(5458.6669921875) * x * x * x * x * x * y * y * z - std::sqrt(14700.0) * x * x * x * x * x * z * z * z + std::sqrt(1295.5810546875) * x * x * x * y * y * y * y * z - std::sqrt(3675.0) * x * x * x * y * y * z * z * z + std::sqrt(8268.75) * x * x * x * z * z * z * z * z - std::sqrt(2616.2841796875) * x * y * y * y * y * y * y * z - std::sqrt(33075.0) * x * y * y * y * y * z * z * z + std::sqrt(8268.75) * x * y * y * z * z * z * z * z) + e_2 * (-std::sqrt(80749.51171875) * x * x * x * x * x * z + std::sqrt(149354.296875) * x * x * x * y * y * z - std::sqrt(33075.0) * x * x * x * z * z * z - std::sqrt(362920.60546875) * x * y * y * y * y * z - std::sqrt(33075.0) * x * y * y * z * z * z + std::sqrt(33075.0) * x * z * z * z * z * z) + e_3 * (-std::sqrt(1000518.75) * x * x * x * z - std::sqrt(1000518.75) * x * y * y * z + std::sqrt(132300.0) * x * z * z * z) + e_4 * (-std::sqrt(529200.0) * x * z);
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

        pc_100[k] = e_0 * (-std::sqrt(32.2998046875) * x * x * x * x * x * x * x * x * z * z + std::sqrt(807.4951171875) * x * x * x * x * x * x * y * y * z * z + std::sqrt(229.6875) * x * x * x * x * x * x * z * z * z * z - std::sqrt(290.6982421875) * x * x * x * x * y * y * y * y * z * z - std::sqrt(8268.75) * x * x * x * x * y * y * z * z * z * z - std::sqrt(2616.2841796875) * x * x * y * y * y * y * y * y * z * z + std::sqrt(18604.6875) * x * x * y * y * y * y * z * z * z * z) + e_1 * (-std::sqrt(32.2998046875) * x * x * x * x * x * x * x * x + std::sqrt(807.4951171875) * x * x * x * x * x * x * y * y - std::sqrt(1582.6904296875) * x * x * x * x * x * x * z * z - std::sqrt(290.6982421875) * x * x * x * x * y * y * y * y - std::sqrt(49128.0029296875) * x * x * x * x * y * y * z * z + std::sqrt(18604.6875) * x * x * x * x * z * z * z * z - std::sqrt(2616.2841796875) * x * x * y * y * y * y * y * y - std::sqrt(2616.2841796875) * x * x * y * y * y * y * z * z + std::sqrt(74418.75) * x * x * y * y * z * z * z * z - std::sqrt(2616.2841796875) * y * y * y * y * y * y * z * z + std::sqrt(18604.6875) * y * y * y * y * z * z * z * z) + e_2 * (-std::sqrt(7267.4560546875) * x * x * x * x * x * x + std::sqrt(2616.2841796875) * x * x * x * x * y * y - std::sqrt(211919.0185546875) * x * x * y * y * y * y + std::sqrt(297675.0) * x * x * z * z * z * z - std::sqrt(2616.2841796875) * y * y * y * y * y * y + std::sqrt(297675.0) * y * y * z * z * z * z) + e_3 * (-std::sqrt(167442.1875) * x * x * x * x - std::sqrt(669768.75) * x * x * y * y + std::sqrt(1190700.0) * x * x * z * z - std::sqrt(167442.1875) * y * y * y * y + std::sqrt(1190700.0) * y * y * z * z + std::sqrt(132300.0) * z * z * z * z) + e_4 * (-std::sqrt(297675.0) * x * x - std::sqrt(297675.0) * y * y + std::sqrt(1190700.0) * z * z);
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

        pc_101[k] = e_0 * (-std::sqrt(1.07666015625) * x * x * x * x * x * x * x * x * x * z + std::sqrt(68.90625) * x * x * x * x * x * x * x * y * y * z + std::sqrt(107.666015625) * x * x * x * x * x * x * x * z * z * z - std::sqrt(107.666015625) * x * x * x * x * x * y * y * y * y * z - std::sqrt(8720.947265625) * x * x * x * x * x * y * y * z * z * z - std::sqrt(275.625) * x * x * x * y * y * y * y * y * y * z + std::sqrt(38867.431640625) * x * x * x * y * y * y * y * z * z * z + std::sqrt(9.68994140625) * x * y * y * y * y * y * y * y * y * z - std::sqrt(968.994140625) * x * y * y * y * y * y * y * z * z * z) + e_1 * (std::sqrt(4.306640625) * x * x * x * x * x * x * x * z - std::sqrt(28255.869140625) * x * x * x * x * x * y * y * z + std::sqrt(15503.90625) * x * x * x * x * x * z * z * z + std::sqrt(56955.322265625) * x * x * x * y * y * y * y * z + std::sqrt(62015.625) * x * x * x * y * y * z * z * z - std::sqrt(3139.541015625) * x * y * y * y * y * y * y * z + std::sqrt(15503.90625) * x * y * y * y * y * z * z * z) + e_2 * (std::sqrt(15503.90625) * x * x * x * x * x * z + std::sqrt(62015.625) * x * x * x * y * y * z + std::sqrt(558140.625) * x * x * x * z * z * z + std::sqrt(15503.90625) * x * y * y * y * y * z + std::sqrt(558140.625) * x * y * y * z * z * z) + e_3 * (std::sqrt(1550390.625) * x * x * x * z + std::sqrt(1550390.625) * x * y * y * z + std::sqrt(992250.0) * x * z * z * z) + e_4 * (std::sqrt(3969000.0) * x * z);
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

        pc_102[k] = e_0 * (std::sqrt(23.6865234375) * x * x * x * x * x * x * x * x * z * z - std::sqrt(4003.0224609375) * x * x * x * x * x * x * y * y * z * z + std::sqrt(29015.9912109375) * x * x * x * x * y * y * y * y * z * z - std::sqrt(5329.4677734375) * x * x * y * y * y * y * y * y * z * z) + e_1 * (std::sqrt(23.6865234375) * x * x * x * x * x * x * x * x - std::sqrt(4003.0224609375) * x * x * x * x * x * x * y * y + std::sqrt(5329.4677734375) * x * x * x * x * x * x * z * z + std::sqrt(29015.9912109375) * x * x * x * x * y * y * y * y + std::sqrt(5329.4677734375) * x * x * x * x * y * y * z * z - std::sqrt(5329.4677734375) * x * x * y * y * y * y * y * y - std::sqrt(5329.4677734375) * x * x * y * y * y * y * z * z - std::sqrt(5329.4677734375) * y * y * y * y * y * y * z * z) + e_2 * (std::sqrt(5329.4677734375) * x * x * x * x * x * x + std::sqrt(5329.4677734375) * x * x * x * x * y * y + std::sqrt(341085.9375) * x * x * x * x * z * z - std::sqrt(5329.4677734375) * x * x * y * y * y * y - std::sqrt(5329.4677734375) * y * y * y * y * y * y - std::sqrt(341085.9375) * y * y * y * y * z * z) + e_3 * (std::sqrt(341085.9375) * x * x * x * x + std::sqrt(1364343.75) * x * x * z * z - std::sqrt(341085.9375) * y * y * y * y - std::sqrt(1364343.75) * y * y * z * z) + e_4 * (std::sqrt(1364343.75) * x * x - std::sqrt(1364343.75) * y * y);

        pc_103[k] = e_0 * (std::sqrt(1.973876953125) * x * x * x * x * x * x * x * x * x * z - std::sqrt(639.5361328125) * x * x * x * x * x * x * x * y * y * z + std::sqrt(7105.95703125) * x * x * x * x * x * y * y * y * y * z - std::sqrt(4176.7236328125) * x * x * x * y * y * y * y * y * y * z + std::sqrt(17.764892578125) * x * y * y * y * y * y * y * y * y * z) + e_1 * (std::sqrt(639.5361328125) * x * x * x * x * x * x * x * z - std::sqrt(639.5361328125) * x * x * x * x * x * y * y * z - std::sqrt(15988.4033203125) * x * x * x * y * y * y * y * z - std::sqrt(5755.8251953125) * x * y * y * y * y * y * y * z) + e_2 * (std::sqrt(63953.61328125) * x * x * x * x * x * z - std::sqrt(255814.453125) * x * x * x * y * y * z - std::sqrt(575582.51953125) * x * y * y * y * y * z) + e_3 * (std::sqrt(454781.25) * x * x * x * z - std::sqrt(4093031.25) * x * y * y * z);
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

        pc_104[k] = e_0 * (std::sqrt(8.8824462890625) * x * x * x * x * x * x * x * x * x * y - std::sqrt(773.759765625) * x * x * x * x * x * x * x * y * y * y + std::sqrt(4299.10400390625) * x * x * x * x * x * y * y * y * y * y - std::sqrt(773.759765625) * x * x * x * y * y * y * y * y * y * y + std::sqrt(8.8824462890625) * x * y * y * y * y * y * y * y * y * y) + e_1 * (std::sqrt(568.4765625) * x * x * x * x * x * x * x * y + std::sqrt(5116.2890625) * x * x * x * x * x * y * y * y + std::sqrt(5116.2890625) * x * x * x * y * y * y * y * y + std::sqrt(568.4765625) * x * y * y * y * y * y * y * y) + e_2 * (std::sqrt(127907.2265625) * x * x * x * x * x * y + std::sqrt(511628.90625) * x * x * x * y * y * y + std::sqrt(127907.2265625) * x * y * y * y * y * y) + e_3 * (std::sqrt(3638250.0) * x * x * x * y + std::sqrt(3638250.0) * x * y * y * y) + e_4 * (std::sqrt(8186062.5) * x * y);

        pc_105[k] = e_0 * (std::sqrt(74.0203857421875) * x * x * x * x * x * x * x * x * y * z - std::sqrt(4737.3046875) * x * x * x * x * x * x * y * y * y * z + std::sqrt(12897.31201171875) * x * x * x * x * y * y * y * y * y * z - std::sqrt(757.96875) * x * x * y * y * y * y * y * y * y * z + std::sqrt(2.9608154296875) * y * y * y * y * y * y * y * y * y * z) + e_1 * (std::sqrt(1184.326171875) * x * x * x * x * x * x * y * z + std::sqrt(10658.935546875) * x * x * x * x * y * y * y * z + std::sqrt(10658.935546875) * x * x * y * y * y * y * y * z + std::sqrt(1184.326171875) * y * y * y * y * y * y * y * z) + e_2 * (std::sqrt(170542.96875) * x * x * x * x * y * z + std::sqrt(682171.875) * x * x * y * y * y * z + std::sqrt(170542.96875) * y * y * y * y * y * z) + e_3 * (std::sqrt(2728687.5) * x * x * y * z + std::sqrt(2728687.5) * y * y * y * z) + e_4 * (std::sqrt(2728687.5) * y * z);

        pc_106[k] = e_0 * (-std::sqrt(2.1533203125) * x * x * x * x * x * x * x * x * x * y + std::sqrt(77.51953125) * x * x * x * x * x * x * x * y * y * y + std::sqrt(215.33203125) * x * x * x * x * x * x * x * y * z * z - std::sqrt(10551.26953125) * x * x * x * x * x * y * y * y * z * z - std::sqrt(77.51953125) * x * x * x * y * y * y * y * y * y * y + std::sqrt(10551.26953125) * x * x * x * y * y * y * y * y * z * z + std::sqrt(2.1533203125) * x * y * y * y * y * y * y * y * y * y - std::sqrt(215.33203125) * x * y * y * y * y * y * y * y * z * z) + e_1 * (-std::sqrt(137.8125) * x * x * x * x * x * x * x * y + std::sqrt(6752.8125) * x * x * x * x * x * y * y * y - std::sqrt(6752.8125) * x * x * x * y * y * y * y * y + std::sqrt(137.8125) * x * y * y * y * y * y * y * y);
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

        pc_107[k] = e_0 * (-std::sqrt(36.3372802734375) * x * x * x * x * x * x * x * x * y * z + std::sqrt(1033.59375) * x * x * x * x * x * x * y * y * y * z + std::sqrt(258.3984375) * x * x * x * x * x * x * y * z * z * z + std::sqrt(403.74755859375) * x * x * x * x * y * y * y * y * y * z - std::sqrt(10364.6484375) * x * x * x * x * y * y * y * z * z * z - std::sqrt(258.3984375) * x * x * y * y * y * y * y * y * y * z + std::sqrt(2325.5859375) * x * x * y * y * y * y * y * z * z * z + std::sqrt(4.0374755859375) * y * y * y * y * y * y * y * y * y * z - std::sqrt(28.7109375) * y * y * y * y * y * y * y * z * z * z) + e_1 * (-std::sqrt(581.396484375) * x * x * x * x * x * x * y * z + std::sqrt(142700.537109375) * x * x * x * x * y * y * y * z - std::sqrt(4134.375) * x * x * x * x * y * z * z * z - std::sqrt(5232.568359375) * x * x * y * y * y * y * y * z - std::sqrt(16537.5) * x * x * y * y * y * z * z * z + std::sqrt(1614.990234375) * y * y * y * y * y * y * y * z - std::sqrt(4134.375) * y * y * y * y * y * z * z * z) + e_2 * (std::sqrt(83721.09375) * x * x * x * x * y * z + std::sqrt(334884.375) * x * x * y * y * y * z - std::sqrt(148837.5) * x * x * y * z * z * z + std::sqrt(83721.09375) * y * y * y * y * y * z - std::sqrt(148837.5) * y * y * y * z * z * z) + e_3 * (std::sqrt(595350.0) * x * x * y * z + std::sqrt(595350.0) * y * y * y * z - std::sqrt(264600.0) * y * z * z * z) + e_4 * (std::sqrt(148837.5) * y * z);
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

        pc_108[k] = e_0 * (std::sqrt(0.4486083984375) * x * x * x * x * x * x * x * x * x * y - std::sqrt(7.177734375) * x * x * x * x * x * x * x * y * y * y - std::sqrt(114.84375) * x * x * x * x * x * x * x * y * z * z - std::sqrt(44.86083984375) * x * x * x * x * x * y * y * y * y * y + std::sqrt(2871.09375) * x * x * x * x * x * y * y * y * z * z + std::sqrt(114.84375) * x * x * x * x * x * y * z * z * z * z - std::sqrt(7.177734375) * x * x * x * y * y * y * y * y * y * y + std::sqrt(2871.09375) * x * x * x * y * y * y * y * y * z * z - std::sqrt(4134.375) * x * x * x * y * y * y * z * z * z * z + std::sqrt(0.4486083984375) * x * y * y * y * y * y * y * y * y * y - std::sqrt(114.84375) * x * y * y * y * y * y * y * y * z * z + std::sqrt(114.84375) * x * y * y * y * y * y * z * z * z * z) + e_1 * (std::sqrt(28.7109375) * x * x * x * x * x * x * x * y - std::sqrt(4852.1484375) * x * x * x * x * x * y * y * y - std::sqrt(4852.1484375) * x * x * x * y * y * y * y * y + std::sqrt(470400.0) * x * x * x * y * y * y * z * z - std::sqrt(7350.0) * x * x * x * y * z * z * z * z + std::sqrt(28.7109375) * x * y * y * y * y * y * y * y - std::sqrt(7350.0) * x * y * y * y * z * z * z * z) + e_2 * (-std::sqrt(2325.5859375) * x * x * x * x * x * y - std::sqrt(125064.84375) * x * x * x * y * y * y + std::sqrt(595350.0) * x * x * x * y * z * z - std::sqrt(2325.5859375) * x * y * y * y * y * y + std::sqrt(595350.0) * x * y * y * y * z * z - std::sqrt(66150.0) * x * y * z * z * z * z) + e_3 * (-std::sqrt(66150.0) * x * x * x * y - std::sqrt(66150.0) * x * y * y * y + std::sqrt(1058400.0) * x * y * z * z) + e_4 * (-std::sqrt(16537.5) * x * y);
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

        pc_109[k] = e_0 * (std::sqrt(4.486083984375) * x * x * x * x * x * x * x * x * y * z - std::sqrt(71.77734375) * x * x * x * x * x * x * y * y * y * z - std::sqrt(71.77734375) * x * x * x * x * x * x * y * z * z * z - std::sqrt(448.6083984375) * x * x * x * x * y * y * y * y * y * z + std::sqrt(1794.43359375) * x * x * x * x * y * y * y * z * z * z + std::sqrt(11.484375) * x * x * x * x * y * z * z * z * z * z - std::sqrt(71.77734375) * x * x * y * y * y * y * y * y * y * z + std::sqrt(1794.43359375) * x * x * y * y * y * y * y * z * z * z - std::sqrt(413.4375) * x * x * y * y * y * z * z * z * z * z + std::sqrt(4.486083984375) * y * y * y * y * y * y * y * y * y * z - std::sqrt(71.77734375) * y * y * y * y * y * y * y * z * z * z + std::sqrt(11.484375) * y * y * y * y * y * z * z * z * z * z) + e_1 * (std::sqrt(71.77734375) * x * x * x * x * x * x * y * z - std::sqrt(44860.83984375) * x * x * x * x * y * y * y * z + std::sqrt(1148.4375) * x * x * x * x * y * z * z * z - std::sqrt(31653.80859375) * x * x * y * y * y * y * y * z + std::sqrt(225093.75) * x * x * y * y * y * z * z * z - std::sqrt(1653.75) * x * x * y * z * z * z * z * z + std::sqrt(1794.43359375) * y * y * y * y * y * y * y * z - std::sqrt(10335.9375) * y * y * y * y * y * z * z * z + std::sqrt(183.75) * y * y * y * z * z * z * z * z) + e_2 * (-std::sqrt(41343.75) * x * x * x * x * y * z - std::sqrt(661500.0) * x * x * y * y * y * z + std::sqrt(372093.75) * x * x * y * z * z * z + std::sqrt(41343.75) * y * y * y * y * y * z - std::sqrt(41343.75) * y * y * y * z * z * z) + e_3 * (-std::sqrt(372093.75) * x * x * y * z + std::sqrt(41343.75) * y * y * y * z);
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

        pc_110[k] = e_0 * (-std::sqrt(0.05340576171875) * x * x * x * x * x * x * x * x * x * x + std::sqrt(0.48065185546875) * x * x * x * x * x * x * x * x * y * y + std::sqrt(17.303466796875) * x * x * x * x * x * x * x * x * z * z + std::sqrt(10.467529296875) * x * x * x * x * x * x * y * y * y * y - std::sqrt(276.85546875) * x * x * x * x * x * x * y * y * z * z - std::sqrt(30.76171875) * x * x * x * x * x * x * z * z * z * z + std::sqrt(10.467529296875) * x * x * x * x * y * y * y * y * y * y - std::sqrt(1730.3466796875) * x * x * x * x * y * y * y * y * z * z + std::sqrt(769.04296875) * x * x * x * x * y * y * z * z * z * z + std::sqrt(0.546875) * x * x * x * x * z * z * z * z * z * z + std::sqrt(0.48065185546875) * x * x * y * y * y * y * y * y * y * y - std::sqrt(276.85546875) * x * x * y * y * y * y * y * y * z * z + std::sqrt(769.04296875) * x * x * y * y * y * y * z * z * z * z - std::sqrt(19.6875) * x * x * y * y * z * z * z * z * z * z - std::sqrt(0.05340576171875) * y * y * y * y * y * y * y * y * y * y + std::sqrt(17.303466796875) * y * y * y * y * y * y * y * y * z * z - std::sqrt(30.76171875) * y * y * y * y * y * y * z * z * z * z + std::sqrt(0.546875) * y * y * y * y * z * z * z * z * z * z) + e_1 * (-std::sqrt(30.76171875) * x * x * x * x * x * x * x * x + std::sqrt(492.1875) * x * x * x * x * x * x * y * y + std::sqrt(4429.6875) * x * x * x * x * x * x * z * z + std::sqrt(3076.171875) * x * x * x * x * y * y * y * y - std::sqrt(110742.1875) * x * x * x * x * y * y * z * z - std::sqrt(1968.75) * x * x * x * x * z * z * z * z + std::sqrt(492.1875) * x * x * y * y * y * y * y * y - std::sqrt(110742.1875) * x * x * y * y * y * y * z * z + std::sqrt(70875.0) * x * x * y * y * z * z * z * z - std::sqrt(30.76171875) * y * y * y * y * y * y * y * y + std::sqrt(4429.6875) * y * y * y * y * y * y * z * z - std::sqrt(1968.75) * y * y * y * y * z * z * z * z) + e_2 * (-std::sqrt(1107.421875) * x * x * x * x * x * x + std::sqrt(27685.546875) * x * x * x * x * y * y + std::sqrt(39867.1875) * x * x * x * x * z * z + std::sqrt(27685.546875) * x * x * y * y * y * y - std::sqrt(1435218.75) * x * x * y * y * z * z - std::sqrt(1107.421875) * y * y * y * y * y * y + std::sqrt(39867.1875) * y * y * y * y * z * z) + e_3 * (-std::sqrt(1968.75) * x * x * x * x + std::sqrt(70875.0) * x * x * y * y - std::sqrt(1968.75) * y * y * y * y);
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

        pc_111[k] = e_0 * (std::sqrt(4.486083984375) * x * x * x * x * x * x * x * x * x * z - std::sqrt(71.77734375) * x * x * x * x * x * x * x * y * y * z - std::sqrt(71.77734375) * x * x * x * x * x * x * x * z * z * z - std::sqrt(448.6083984375) * x * x * x * x * x * y * y * y * y * z + std::sqrt(1794.43359375) * x * x * x * x * x * y * y * z * z * z + std::sqrt(11.484375) * x * x * x * x * x * z * z * z * z * z - std::sqrt(71.77734375) * x * x * x * y * y * y * y * y * y * z + std::sqrt(1794.43359375) * x * x * x * y * y * y * y * z * z * z - std::sqrt(413.4375) * x * x * x * y * y * z * z * z * z * z + std::sqrt(4.486083984375) * x * y * y * y * y * y * y * y * y * z - std::sqrt(71.77734375) * x * y * y * y * y * y * y * z * z * z + std::sqrt(11.484375) * x * y * y * y * y * z * z * z * z * z) + e_1 * (std::sqrt(1794.43359375) * x * x * x * x * x * x * x * z - std::sqrt(31653.80859375) * x * x * x * x * x * y * y * z - std::sqrt(10335.9375) * x * x * x * x * x * z * z * z - std::sqrt(44860.83984375) * x * x * x * y * y * y * y * z + std::sqrt(225093.75) * x * x * x * y * y * z * z * z + std::sqrt(183.75) * x * x * x * z * z * z * z * z + std::sqrt(71.77734375) * x * y * y * y * y * y * y * z + std::sqrt(1148.4375) * x * y * y * y * y * z * z * z - std::sqrt(1653.75) * x * y * y * z * z * z * z * z) + e_2 * (std::sqrt(41343.75) * x * x * x * x * x * z - std::sqrt(661500.0) * x * x * x * y * y * z - std::sqrt(41343.75) * x * x * x * z * z * z - std::sqrt(41343.75) * x * y * y * y * y * z + std::sqrt(372093.75) * x * y * y * z * z * z) + e_3 * (std::sqrt(41343.75) * x * x * x * z - std::sqrt(372093.75) * x * y * y * z);
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

        pc_112[k] = e_0 * (std::sqrt(0.112152099609375) * x * x * x * x * x * x * x * x * x * x - std::sqrt(2.803802490234375) * x * x * x * x * x * x * x * x * y * y - std::sqrt(28.7109375) * x * x * x * x * x * x * x * x * z * z - std::sqrt(4.0374755859375) * x * x * x * x * x * x * y * y * y * y + std::sqrt(1033.59375) * x * x * x * x * x * x * y * y * z * z + std::sqrt(28.7109375) * x * x * x * x * x * x * z * z * z * z + std::sqrt(4.0374755859375) * x * x * x * x * y * y * y * y * y * y - std::sqrt(1406.8359375) * x * x * x * x * y * y * z * z * z * z + std::sqrt(2.803802490234375) * x * x * y * y * y * y * y * y * y * y - std::sqrt(1033.59375) * x * x * y * y * y * y * y * y * z * z + std::sqrt(1406.8359375) * x * x * y * y * y * y * z * z * z * z - std::sqrt(0.112152099609375) * y * y * y * y * y * y * y * y * y * y + std::sqrt(28.7109375) * y * y * y * y * y * y * y * y * z * z - std::sqrt(28.7109375) * y * y * y * y * y * y * z * z * z * z) + e_1 * (std::sqrt(64.599609375) * x * x * x * x * x * x * x * x - std::sqrt(717.7734375) * x * x * x * x * x * x * y * y - std::sqrt(7350.0) * x * x * x * x * x * x * z * z + std::sqrt(66150.0) * x * x * x * x * y * y * z * z + std::sqrt(1837.5) * x * x * x * x * z * z * z * z + std::sqrt(717.7734375) * x * x * y * y * y * y * y * y - std::sqrt(66150.0) * x * x * y * y * y * y * z * z - std::sqrt(64.599609375) * y * y * y * y * y * y * y * y + std::sqrt(7350.0) * y * y * y * y * y * y * z * z - std::sqrt(1837.5) * y * y * y * y * z * z * z * z) + e_2 * (std::sqrt(3165.380859375) * x * x * x * x * x * x - std::sqrt(5232.568359375) * x * x * x * x * y * y - std::sqrt(148837.5) * x * x * x * x * z * z + std::sqrt(5232.568359375) * x * x * y * y * y * y + std::sqrt(16537.5) * x * x * z * z * z * z - std::sqrt(3165.380859375) * y * y * y * y * y * y + std::sqrt(148837.5) * y * y * y * y * z * z - std::sqrt(16537.5) * y * y * z * z * z * z) + e_3 * (std::sqrt(16537.5) * x * x * x * x - std::sqrt(264600.0) * x * x * z * z - std::sqrt(16537.5) * y * y * y * y + std::sqrt(264600.0) * y * y * z * z) + e_4 * (std::sqrt(4134.375) * x * x - std::sqrt(4134.375) * y * y);
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

        pc_113[k] = e_0 * (-std::sqrt(4.0374755859375) * x * x * x * x * x * x * x * x * x * z + std::sqrt(258.3984375) * x * x * x * x * x * x * x * y * y * z + std::sqrt(28.7109375) * x * x * x * x * x * x * x * z * z * z - std::sqrt(403.74755859375) * x * x * x * x * x * y * y * y * y * z - std::sqrt(2325.5859375) * x * x * x * x * x * y * y * z * z * z - std::sqrt(1033.59375) * x * x * x * y * y * y * y * y * y * z + std::sqrt(10364.6484375) * x * x * x * y * y * y * y * z * z * z + std::sqrt(36.3372802734375) * x * y * y * y * y * y * y * y * y * z - std::sqrt(258.3984375) * x * y * y * y * y * y * y * z * z * z) + e_1 * (-std::sqrt(1614.990234375) * x * x * x * x * x * x * x * z + std::sqrt(5232.568359375) * x * x * x * x * x * y * y * z + std::sqrt(4134.375) * x * x * x * x * x * z * z * z - std::sqrt(142700.537109375) * x * x * x * y * y * y * y * z + std::sqrt(16537.5) * x * x * x * y * y * z * z * z + std::sqrt(581.396484375) * x * y * y * y * y * y * y * z + std::sqrt(4134.375) * x * y * y * y * y * z * z * z) + e_2 * (-std::sqrt(83721.09375) * x * x * x * x * x * z - std::sqrt(334884.375) * x * x * x * y * y * z + std::sqrt(148837.5) * x * x * x * z * z * z - std::sqrt(83721.09375) * x * y * y * y * y * z + std::sqrt(148837.5) * x * y * y * z * z * z) + e_3 * (-std::sqrt(595350.0) * x * x * x * z - std::sqrt(595350.0) * x * y * y * z + std::sqrt(264600.0) * x * z * z * z) + e_4 * (-std::sqrt(148837.5) * x * z);
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

        pc_114[k] = e_0 * (-std::sqrt(0.13458251953125) * x * x * x * x * x * x * x * x * x * x + std::sqrt(16.28448486328125) * x * x * x * x * x * x * x * x * y * y + std::sqrt(13.458251953125) * x * x * x * x * x * x * x * x * z * z - std::sqrt(90.977783203125) * x * x * x * x * x * x * y * y * y * y - std::sqrt(1937.98828125) * x * x * x * x * x * x * y * y * z * z - std::sqrt(90.977783203125) * x * x * x * x * y * y * y * y * y * y + std::sqrt(19433.7158203125) * x * x * x * x * y * y * y * y * z * z + std::sqrt(16.28448486328125) * x * x * y * y * y * y * y * y * y * y - std::sqrt(1937.98828125) * x * x * y * y * y * y * y * y * z * z - std::sqrt(0.13458251953125) * y * y * y * y * y * y * y * y * y * y + std::sqrt(13.458251953125) * y * y * y * y * y * y * y * y * z * z) + e_1 * (-std::sqrt(77.51953125) * x * x * x * x * x * x * x * x + std::sqrt(137.8125) * x * x * x * x * x * x * y * y + std::sqrt(3445.3125) * x * x * x * x * x * x * z * z - std::sqrt(21533.203125) * x * x * x * x * y * y * y * y + std::sqrt(31007.8125) * x * x * x * x * y * y * z * z + std::sqrt(137.8125) * x * x * y * y * y * y * y * y + std::sqrt(31007.8125) * x * x * y * y * y * y * z * z - std::sqrt(77.51953125) * y * y * y * y * y * y * y * y + std::sqrt(3445.3125) * y * y * y * y * y * y * z * z) + e_2 * (-std::sqrt(7751.953125) * x * x * x * x * x * x - std::sqrt(69767.578125) * x * x * x * x * y * y + std::sqrt(279070.3125) * x * x * x * x * z * z - std::sqrt(69767.578125) * x * x * y * y * y * y + std::sqrt(1116281.25) * x * x * y * y * z * z - std::sqrt(7751.953125) * y * y * y * y * y * y + std::sqrt(279070.3125) * y * y * y * y * z * z) + e_3 * (-std::sqrt(124031.25) * x * x * x * x - std::sqrt(496125.0) * x * x * y * y + std::sqrt(1984500.0) * x * x * z * z - std::sqrt(124031.25) * y * y * y * y + std::sqrt(1984500.0) * y * y * z * z) + e_4 * (-std::sqrt(124031.25) * x * x - std::sqrt(124031.25) * y * y + std::sqrt(496125.0) * z * z);
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

        pc_115[k] = e_0 * (std::sqrt(2.9608154296875) * x * x * x * x * x * x * x * x * x * z - std::sqrt(757.96875) * x * x * x * x * x * x * x * y * y * z + std::sqrt(12897.31201171875) * x * x * x * x * x * y * y * y * y * z - std::sqrt(4737.3046875) * x * x * x * y * y * y * y * y * y * z + std::sqrt(74.0203857421875) * x * y * y * y * y * y * y * y * y * z) + e_1 * (std::sqrt(1184.326171875) * x * x * x * x * x * x * x * z + std::sqrt(10658.935546875) * x * x * x * x * x * y * y * z + std::sqrt(10658.935546875) * x * x * x * y * y * y * y * z + std::sqrt(1184.326171875) * x * y * y * y * y * y * y * z) + e_2 * (std::sqrt(170542.96875) * x * x * x * x * x * z + std::sqrt(682171.875) * x * x * x * y * y * z + std::sqrt(170542.96875) * x * y * y * y * y * z) + e_3 * (std::sqrt(2728687.5) * x * x * x * z + std::sqrt(2728687.5) * x * y * y * z) + e_4 * (std::sqrt(2728687.5) * x * z);

        pc_116[k] = e_0 * (std::sqrt(0.246734619140625) * x * x * x * x * x * x * x * x * x * x - std::sqrt(108.80996704101562) * x * x * x * x * x * x * x * x * y * y + std::sqrt(2772.3101806640625) * x * x * x * x * x * x * y * y * y * y - std::sqrt(2772.3101806640625) * x * x * x * x * y * y * y * y * y * y + std::sqrt(108.80996704101562) * x * x * y * y * y * y * y * y * y * y - std::sqrt(0.246734619140625) * y * y * y * y * y * y * y * y * y * y) + e_1 * (std::sqrt(142.119140625) * x * x * x * x * x * x * x * x + std::sqrt(568.4765625) * x * x * x * x * x * x * y * y - std::sqrt(568.4765625) * x * x * y * y * y * y * y * y - std::sqrt(142.119140625) * y * y * y * y * y * y * y * y) + e_2 * (std::sqrt(31976.806640625) * x * x * x * x * x * x + std::sqrt(31976.806640625) * x * x * x * x * y * y - std::sqrt(31976.806640625) * x * x * y * y * y * y - std::sqrt(31976.806640625) * y * y * y * y * y * y) + e_3 * (std::sqrt(909562.5) * x * x * x * x - std::sqrt(909562.5) * y * y * y * y) + e_4 * (std::sqrt(2046515.625) * x * x - std::sqrt(2046515.625) * y * y);
    }

    // NOTE: the atom pairs beyond the reach of every pair of primitives have no
    // contribution and are set to zero.

    for (size_t m = 0; m < 117; m++)
    {
        auto *pv = values + m * nvalues;

        std::fill(pv + nmax, pv + nvalues, 0.0);
    }
}

}  // namespace simdovl
