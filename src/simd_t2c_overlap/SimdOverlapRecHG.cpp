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



#include "SimdOverlapRecHG.hpp"

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
compute_hg_overlap(double               *values,
                   const size_t          nvalues,
                   const CBasisFunction &bra,
                   const CBasisFunction &ket,
                   const CSimdMatrix    &coordinates,
                   const double          threshold) -> void
{
    if ((bra.get_angular_momentum() != 5) || (ket.get_angular_momentum() != 4))
    {
        errors::assertMsgCritical(
            false, std::string("SimdOverlapRecHG.compute_hg_overlap: Basis functions must be of angular momenta five and four"));
    }

    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("SimdOverlapRecHG.compute_hg_overlap: Number of values exceeds number of atom pairs"));
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
        std::fill(values, values + 99 * nvalues, 0.0);

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

        const auto f_0 = fbase * fal * fal * fal * fal * fal * fbe * fbe * fbe * fbe;

        const auto f_1 = fbase * fal * fal * fal * fal * fbe * fbe * fbe * fh;

        const auto f_2 = fbase * fal * fal * fal * fbe * fbe * fh * fh;

        const auto f_3 = fbase * fal * fal * fbe * fh * fh * fh;

        const auto f_4 = fbase * fal * fh * fh * fh * fh;

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

    // NOTE: the components are formed in 25 loops, as the vectorizer runs out
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

        pc_0[k] = e_0 * (std::sqrt(107.666015625) * x * x * x * x * x * x * x * y * y - std::sqrt(968.994140625) * x * x * x * x * x * y * y * y * y + std::sqrt(521.103515625) * x * x * x * y * y * y * y * y * y - std::sqrt(4.306640625) * x * y * y * y * y * y * y * y * y) + e_1 * (std::sqrt(107.666015625) * x * x * x * x * x * x * x + std::sqrt(968.994140625) * x * x * x * x * x * y * y + std::sqrt(968.994140625) * x * x * x * y * y * y * y + std::sqrt(107.666015625) * x * y * y * y * y * y * y) + e_2 * (std::sqrt(15503.90625) * x * x * x * x * x + std::sqrt(62015.625) * x * x * x * y * y + std::sqrt(15503.90625) * x * y * y * y * y) + e_3 * (std::sqrt(248062.5) * x * x * x + std::sqrt(248062.5) * x * y * y) + e_4 * (std::sqrt(248062.5) * x);

        pc_1[k] = e_0 * (std::sqrt(484.4970703125) * x * x * x * x * x * x * y * y * z - std::sqrt(2637.8173828125) * x * x * x * x * y * y * y * y * z + std::sqrt(363.9111328125) * x * x * y * y * y * y * y * y * z - std::sqrt(2.1533203125) * y * y * y * y * y * y * y * y * z) + e_1 * (std::sqrt(484.4970703125) * x * x * x * x * x * x * z + std::sqrt(484.4970703125) * x * x * x * x * y * y * z - std::sqrt(484.4970703125) * x * x * y * y * y * y * z - std::sqrt(484.4970703125) * y * y * y * y * y * y * z) + e_2 * (std::sqrt(31007.8125) * x * x * x * x * z - std::sqrt(31007.8125) * y * y * y * y * z) + e_3 * (std::sqrt(124031.25) * x * x * z - std::sqrt(124031.25) * y * y * z);

        pc_2[k] = e_0 * (-std::sqrt(15.380859375) * x * x * x * x * x * x * x * y * y + std::sqrt(15.380859375) * x * x * x * x * x * y * y * y * y + std::sqrt(553.7109375) * x * x * x * x * x * y * y * z * z + std::sqrt(49.833984375) * x * x * x * y * y * y * y * y * y - std::sqrt(2214.84375) * x * x * x * y * y * y * y * z * z - std::sqrt(0.615234375) * x * y * y * y * y * y * y * y * y + std::sqrt(22.1484375) * x * y * y * y * y * y * y * z * z) + e_1 * (-std::sqrt(15.380859375) * x * x * x * x * x * x * x - std::sqrt(1245.849609375) * x * x * x * x * x * y * y + std::sqrt(553.7109375) * x * x * x * x * x * z * z + std::sqrt(9613.037109375) * x * x * x * y * y * y * y - std::sqrt(2214.84375) * x * x * x * y * y * z * z + std::sqrt(15.380859375) * x * y * y * y * y * y * y - std::sqrt(4983.3984375) * x * y * y * y * y * z * z) + e_2 * (-std::sqrt(2214.84375) * x * x * x * x * x + std::sqrt(8859.375) * x * x * x * y * y + std::sqrt(8859.375) * x * x * x * z * z + std::sqrt(19933.59375) * x * y * y * y * y - std::sqrt(79734.375) * x * y * y * z * z) + e_3 * (-std::sqrt(8859.375) * x * x * x + std::sqrt(79734.375) * x * y * y);

        pc_3[k] = e_0 * (-std::sqrt(69.2138671875) * x * x * x * x * x * x * y * y * z + std::sqrt(69.2138671875) * x * x * x * x * y * y * y * y * z + std::sqrt(123.046875) * x * x * x * x * y * y * z * z * z + std::sqrt(224.2529296875) * x * x * y * y * y * y * y * y * z - std::sqrt(492.1875) * x * x * y * y * y * y * z * z * z - std::sqrt(2.7685546875) * y * y * y * y * y * y * y * y * z + std::sqrt(4.921875) * y * y * y * y * y * y * z * z * z) + e_1 * (-std::sqrt(69.2138671875) * x * x * x * x * x * x * z - std::sqrt(1730.3466796875) * x * x * x * x * y * y * z + std::sqrt(123.046875) * x * x * x * x * z * z * z + std::sqrt(43258.6669921875) * x * x * y * y * y * y * z - std::sqrt(4429.6875) * x * x * y * y * z * z * z - std::sqrt(622.9248046875) * y * y * y * y * y * y * z + std::sqrt(123.046875) * y * y * y * y * z * z * z) + e_2 * (-std::sqrt(4429.6875) * x * x * x * x * z + std::sqrt(159468.75) * x * x * y * y * z - std::sqrt(4429.6875) * y * y * y * y * z);
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

        pc_4[k] = e_0 * (std::sqrt(1.7303466796875) * x * x * x * x * x * x * x * x * y - std::sqrt(110.7421875) * x * x * x * x * x * x * y * z * z - std::sqrt(13.56591796875) * x * x * x * x * y * y * y * y * y + std::sqrt(110.7421875) * x * x * x * x * y * y * y * z * z + std::sqrt(12.3046875) * x * x * x * x * y * z * z * z * z - std::sqrt(4.4296875) * x * x * y * y * y * y * y * y * y + std::sqrt(358.8046875) * x * x * y * y * y * y * y * z * z - std::sqrt(49.21875) * x * x * y * y * y * z * z * z * z + std::sqrt(0.0692138671875) * y * y * y * y * y * y * y * y * y - std::sqrt(4.4296875) * y * y * y * y * y * y * y * z * z + std::sqrt(0.4921875) * y * y * y * y * y * z * z * z * z) + e_1 * (std::sqrt(692.138671875) * x * x * x * x * x * x * y - std::sqrt(692.138671875) * x * x * x * x * y * y * y - std::sqrt(11074.21875) * x * x * x * x * y * z * z - std::sqrt(2242.529296875) * x * x * y * y * y * y * y + std::sqrt(44296.875) * x * x * y * y * y * z * z + std::sqrt(27.685546875) * y * y * y * y * y * y * y - std::sqrt(442.96875) * y * y * y * y * y * z * z) + e_2 * (std::sqrt(11074.21875) * x * x * x * x * y - std::sqrt(44296.875) * x * x * y * y * y + std::sqrt(442.96875) * y * y * y * y * y);

        pc_5[k] = e_0 * (-std::sqrt(69.2138671875) * x * x * x * x * x * x * x * y * z + std::sqrt(69.2138671875) * x * x * x * x * x * y * y * y * z + std::sqrt(123.046875) * x * x * x * x * x * y * z * z * z + std::sqrt(224.2529296875) * x * x * x * y * y * y * y * y * z - std::sqrt(492.1875) * x * x * x * y * y * y * z * z * z - std::sqrt(2.7685546875) * x * y * y * y * y * y * y * y * z + std::sqrt(4.921875) * x * y * y * y * y * y * z * z * z) + e_1 * (-std::sqrt(13565.91796875) * x * x * x * x * x * y * z + std::sqrt(27685.546875) * x * x * x * y * y * y * z + std::sqrt(1968.75) * x * x * x * y * z * z * z + std::sqrt(276.85546875) * x * y * y * y * y * y * z - std::sqrt(1968.75) * x * y * y * y * z * z * z) + e_2 * (-std::sqrt(70875.0) * x * x * x * y * z + std::sqrt(70875.0) * x * y * y * y * z);

        pc_6[k] = e_0 * (-std::sqrt(3.84521484375) * x * x * x * x * x * x * x * x * y + std::sqrt(15.380859375) * x * x * x * x * x * x * y * y * y + std::sqrt(138.427734375) * x * x * x * x * x * x * y * z * z + std::sqrt(2.4609375) * x * x * x * x * y * y * y * y * y - std::sqrt(1245.849609375) * x * x * x * x * y * y * y * z * z - std::sqrt(15.380859375) * x * x * y * y * y * y * y * y * y + std::sqrt(669.990234375) * x * x * y * y * y * y * y * z * z + std::sqrt(0.15380859375) * y * y * y * y * y * y * y * y * y - std::sqrt(5.537109375) * y * y * y * y * y * y * y * z * z) + e_1 * (-std::sqrt(984.375) * x * x * x * x * x * x * y + std::sqrt(1538.0859375) * x * x * x * x * y * y * y + std::sqrt(4983.3984375) * x * x * x * x * y * z * z - std::sqrt(2214.84375) * x * x * y * y * y * y * y + std::sqrt(2214.84375) * x * x * y * y * y * z * z + std::sqrt(61.5234375) * y * y * y * y * y * y * y - std::sqrt(553.7109375) * y * y * y * y * y * z * z) + e_2 * (-std::sqrt(19933.59375) * x * x * x * x * y - std::sqrt(8859.375) * x * x * y * y * y + std::sqrt(79734.375) * x * x * y * z * z + std::sqrt(2214.84375) * y * y * y * y * y - std::sqrt(8859.375) * y * y * y * z * z) + e_3 * (-std::sqrt(79734.375) * x * x * y + std::sqrt(8859.375) * y * y * y);

        pc_7[k] = e_0 * (std::sqrt(53.8330078125) * x * x * x * x * x * x * x * y * z - std::sqrt(1345.8251953125) * x * x * x * x * x * y * y * y * z + std::sqrt(2069.3408203125) * x * x * x * y * y * y * y * y * z - std::sqrt(19.3798828125) * x * y * y * y * y * y * y * y * z) + e_1 * (std::sqrt(1937.98828125) * x * x * x * x * x * y * z + std::sqrt(7751.953125) * x * x * x * y * y * y * z + std::sqrt(1937.98828125) * x * y * y * y * y * y * z) + e_2 * (std::sqrt(124031.25) * x * x * x * y * z + std::sqrt(124031.25) * x * y * y * y * z) + e_3 * (std::sqrt(496125.0) * x * y * z);
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

        pc_8[k] = e_0 * (std::sqrt(6.7291259765625) * x * x * x * x * x * x * x * x * y - std::sqrt(430.6640625) * x * x * x * x * x * x * y * y * y + std::sqrt(1172.48291015625) * x * x * x * x * y * y * y * y * y - std::sqrt(68.90625) * x * x * y * y * y * y * y * y * y + std::sqrt(0.2691650390625) * y * y * y * y * y * y * y * y * y) + e_1 * (std::sqrt(107.666015625) * x * x * x * x * x * x * y + std::sqrt(968.994140625) * x * x * x * x * y * y * y + std::sqrt(968.994140625) * x * x * y * y * y * y * y + std::sqrt(107.666015625) * y * y * y * y * y * y * y) + e_2 * (std::sqrt(15503.90625) * x * x * x * x * y + std::sqrt(62015.625) * x * x * y * y * y + std::sqrt(15503.90625) * y * y * y * y * y) + e_3 * (std::sqrt(248062.5) * x * x * y + std::sqrt(248062.5) * y * y * y) + e_4 * (std::sqrt(248062.5) * y);

        pc_9[k] = e_0 * (26.25 * x * x * x * x * x * x * y * y * z - 52.5 * x * x * x * x * y * y * y * y * z + 26.25 * x * x * y * y * y * y * y * y * z) + e_1 * (26.25 * x * x * x * x * x * x * z + 78.75 * x * x * x * x * y * y * z + 78.75 * x * x * y * y * y * y * z + 26.25 * y * y * y * y * y * y * z) + e_2 * (236.25 * x * x * x * x * z + 472.5 * x * x * y * y * z + 236.25 * y * y * y * y * z) + e_3 * (630.0 * x * x * z + 630.0 * y * y * z) + e_4 * (315.0 * z);

        pc_10[k] = e_0 * (std::sqrt(3100.78125) * x * x * x * x * x * y * y * z * z - std::sqrt(5512.5) * x * x * x * y * y * y * y * z * z + std::sqrt(344.53125) * x * y * y * y * y * y * y * z * z) + e_1 * (std::sqrt(3100.78125) * x * x * x * x * x * y * y + std::sqrt(3100.78125) * x * x * x * x * x * z * z - std::sqrt(5512.5) * x * x * x * y * y * y * y + std::sqrt(12403.125) * x * x * x * y * y * z * z + std::sqrt(344.53125) * x * y * y * y * y * y * y + std::sqrt(3100.78125) * x * y * y * y * y * z * z) + e_2 * (std::sqrt(3100.78125) * x * x * x * x * x + std::sqrt(12403.125) * x * x * x * y * y + std::sqrt(111628.125) * x * x * x * z * z + std::sqrt(3100.78125) * x * y * y * y * y + std::sqrt(111628.125) * x * y * y * z * z) + e_3 * (std::sqrt(111628.125) * x * x * x + std::sqrt(111628.125) * x * y * y + std::sqrt(198450.0) * x * z * z) + e_4 * (std::sqrt(198450.0) * x);

        pc_11[k] = e_0 * (-std::sqrt(98.4375) * x * x * x * x * x * x * y * y * z + std::sqrt(3543.75) * x * x * x * x * y * y * z * z * z + std::sqrt(98.4375) * x * x * y * y * y * y * y * y * z - std::sqrt(3543.75) * x * x * y * y * y * y * z * z * z) + e_1 * (-std::sqrt(98.4375) * x * x * x * x * x * x * z + std::sqrt(885.9375) * x * x * x * x * y * y * z + std::sqrt(3543.75) * x * x * x * x * z * z * z - std::sqrt(885.9375) * x * x * y * y * y * y * z + std::sqrt(98.4375) * y * y * y * y * y * y * z - std::sqrt(3543.75) * y * y * y * y * z * z * z) + e_2 * (std::sqrt(885.9375) * x * x * x * x * z + std::sqrt(31893.75) * x * x * z * z * z - std::sqrt(885.9375) * y * y * y * y * z - std::sqrt(31893.75) * y * y * z * z * z) + e_3 * (std::sqrt(56700.0) * x * x * z - std::sqrt(56700.0) * y * y * z);
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

        pc_12[k] = e_0 * (-std::sqrt(442.96875) * x * x * x * x * x * y * y * z * z + std::sqrt(787.5) * x * x * x * y * y * z * z * z * z + std::sqrt(442.96875) * x * y * y * y * y * y * y * z * z - std::sqrt(787.5) * x * y * y * y * y * z * z * z * z) + e_1 * (-std::sqrt(442.96875) * x * x * x * x * x * y * y - std::sqrt(442.96875) * x * x * x * x * x * z * z - std::sqrt(1771.875) * x * x * x * y * y * z * z + std::sqrt(787.5) * x * x * x * z * z * z * z + std::sqrt(442.96875) * x * y * y * y * y * y * y + std::sqrt(21705.46875) * x * y * y * y * y * z * z - std::sqrt(7087.5) * x * y * y * z * z * z * z) + e_2 * (-std::sqrt(442.96875) * x * x * x * x * x - std::sqrt(15946.875) * x * x * x * y * y - std::sqrt(1771.875) * x * x * x * z * z + std::sqrt(53599.21875) * x * y * y * y * y + std::sqrt(15946.875) * x * y * y * z * z) + e_3 * (-std::sqrt(15946.875) * x * x * x + std::sqrt(143521.875) * x * y * y);

        pc_13[k] = e_0 * (std::sqrt(11.07421875) * x * x * x * x * x * x * x * y * z + std::sqrt(11.07421875) * x * x * x * x * x * y * y * y * z - std::sqrt(708.75) * x * x * x * x * x * y * z * z * z - std::sqrt(11.07421875) * x * x * x * y * y * y * y * y * z + std::sqrt(78.75) * x * x * x * y * z * z * z * z * z - std::sqrt(11.07421875) * x * y * y * y * y * y * y * y * z + std::sqrt(708.75) * x * y * y * y * y * y * z * z * z - std::sqrt(78.75) * x * y * y * y * z * z * z * z * z) + e_1 * (-std::sqrt(31500.0) * x * x * x * y * z * z * z + std::sqrt(31500.0) * x * y * y * y * z * z * z) + e_2 * (-std::sqrt(70875.0) * x * x * x * y * z + std::sqrt(70875.0) * x * y * y * y * z);

        pc_14[k] = e_0 * (-std::sqrt(442.96875) * x * x * x * x * x * x * y * z * z + std::sqrt(787.5) * x * x * x * x * y * z * z * z * z + std::sqrt(442.96875) * x * x * y * y * y * y * y * z * z - std::sqrt(787.5) * x * x * y * y * y * z * z * z * z) + e_1 * (-std::sqrt(442.96875) * x * x * x * x * x * x * y - std::sqrt(21705.46875) * x * x * x * x * y * z * z + std::sqrt(442.96875) * x * x * y * y * y * y * y + std::sqrt(1771.875) * x * x * y * y * y * z * z + std::sqrt(7087.5) * x * x * y * z * z * z * z + std::sqrt(442.96875) * y * y * y * y * y * z * z - std::sqrt(787.5) * y * y * y * z * z * z * z) + e_2 * (-std::sqrt(53599.21875) * x * x * x * x * y + std::sqrt(15946.875) * x * x * y * y * y - std::sqrt(15946.875) * x * x * y * z * z + std::sqrt(442.96875) * y * y * y * y * y + std::sqrt(1771.875) * y * y * y * z * z) + e_3 * (-std::sqrt(143521.875) * x * x * y + std::sqrt(15946.875) * y * y * y);

        pc_15[k] = e_0 * (-std::sqrt(24.609375) * x * x * x * x * x * x * x * y * z + std::sqrt(24.609375) * x * x * x * x * x * y * y * y * z + std::sqrt(885.9375) * x * x * x * x * x * y * z * z * z + std::sqrt(24.609375) * x * x * x * y * y * y * y * y * z - std::sqrt(3543.75) * x * x * x * y * y * y * z * z * z - std::sqrt(24.609375) * x * y * y * y * y * y * y * y * z + std::sqrt(885.9375) * x * y * y * y * y * y * z * z * z) + e_1 * (-std::sqrt(6300.0) * x * x * x * y * y * y * z + std::sqrt(14175.0) * x * x * x * y * z * z * z + std::sqrt(14175.0) * x * y * y * y * z * z * z) + e_2 * (std::sqrt(3543.75) * x * x * x * y * z + std::sqrt(3543.75) * x * y * y * y * z + std::sqrt(127575.0) * x * y * z * z * z) + e_3 * (std::sqrt(226800.0) * x * y * z);
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

        pc_16[k] = e_0 * (std::sqrt(344.53125) * x * x * x * x * x * x * y * z * z - std::sqrt(5512.5) * x * x * x * x * y * y * y * z * z + std::sqrt(3100.78125) * x * x * y * y * y * y * y * z * z) + e_1 * (std::sqrt(344.53125) * x * x * x * x * x * x * y - std::sqrt(5512.5) * x * x * x * x * y * y * y + std::sqrt(3100.78125) * x * x * x * x * y * z * z + std::sqrt(3100.78125) * x * x * y * y * y * y * y + std::sqrt(12403.125) * x * x * y * y * y * z * z + std::sqrt(3100.78125) * y * y * y * y * y * z * z) + e_2 * (std::sqrt(3100.78125) * x * x * x * x * y + std::sqrt(12403.125) * x * x * y * y * y + std::sqrt(111628.125) * x * x * y * z * z + std::sqrt(3100.78125) * y * y * y * y * y + std::sqrt(111628.125) * y * y * y * z * z) + e_3 * (std::sqrt(111628.125) * x * x * y + std::sqrt(111628.125) * y * y * y + std::sqrt(198450.0) * y * z * z) + e_4 * (std::sqrt(198450.0) * y);

        pc_17[k] = e_0 * (6.5625 * x * x * x * x * x * x * x * y * z - 45.9375 * x * x * x * x * x * y * y * y * z + 45.9375 * x * x * x * y * y * y * y * y * z - 6.5625 * x * y * y * y * y * y * y * y * z);

        pc_18[k] = e_0 * (-std::sqrt(21.533203125) * x * x * x * x * x * x * x * y * y + std::sqrt(2.392578125) * x * x * x * x * x * y * y * y * y + std::sqrt(1378.125) * x * x * x * x * x * y * y * z * z + std::sqrt(21.533203125) * x * x * x * y * y * y * y * y * y - std::sqrt(2450.0) * x * x * x * y * y * y * y * z * z - std::sqrt(2.392578125) * x * y * y * y * y * y * y * y * y + std::sqrt(153.125) * x * y * y * y * y * y * y * z * z) + e_1 * (-std::sqrt(21.533203125) * x * x * x * x * x * x * x - std::sqrt(2605.517578125) * x * x * x * x * x * y * y + std::sqrt(1378.125) * x * x * x * x * x * z * z + std::sqrt(1265.673828125) * x * x * x * y * y * y * y + std::sqrt(5512.5) * x * x * x * y * y * z * z - std::sqrt(289.501953125) * x * y * y * y * y * y * y + std::sqrt(1378.125) * x * y * y * y * y * z * z) + e_2 * (-std::sqrt(3100.78125) * x * x * x * x * x - std::sqrt(12403.125) * x * x * x * y * y + std::sqrt(49612.5) * x * x * x * z * z - std::sqrt(3100.78125) * x * y * y * y * y + std::sqrt(49612.5) * x * y * y * z * z) + e_3 * (-std::sqrt(22050.0) * x * x * x - std::sqrt(22050.0) * x * y * y + std::sqrt(88200.0) * x * z * z) + e_4 * (-std::sqrt(5512.5) * x);

        pc_19[k] = e_0 * (-9.84375 * x * x * x * x * x * x * y * y * z - 3.28125 * x * x * x * x * y * y * y * y * z + 78.75 * x * x * x * x * y * y * z * z * z + 5.46875 * x * x * y * y * y * y * y * y * z - 52.5 * x * x * y * y * y * y * z * z * z - 1.09375 * y * y * y * y * y * y * y * y * z + 8.75 * y * y * y * y * y * y * z * z * z) + e_1 * (-9.84375 * x * x * x * x * x * x * z + 68.90625 * x * x * x * x * y * y * z + 78.75 * x * x * x * x * z * z * z - 95.15625 * x * x * y * y * y * y * z + 157.5 * x * x * y * y * z * z * z + 1.09375 * y * y * y * y * y * y * z + 78.75 * y * y * y * y * z * z * z) + e_2 * (78.75 * x * x * x * x * z + 157.5 * x * x * y * y * z + 315.0 * x * x * z * z * z + 78.75 * y * y * y * y * z + 315.0 * y * y * z * z * z) + e_3 * (525.0 * x * x * z + 525.0 * y * y * z + 210.0 * z * z * z) + e_4 * (420.0 * z);
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

        pc_20[k] = e_0 * (std::sqrt(3.076171875) * x * x * x * x * x * x * x * y * y + std::sqrt(8.544921875) * x * x * x * x * x * y * y * y * y - std::sqrt(602.9296875) * x * x * x * x * x * y * y * z * z + std::sqrt(0.341796875) * x * x * x * y * y * y * y * y * y - std::sqrt(267.96875) * x * x * x * y * y * y * y * z * z + std::sqrt(7087.5) * x * x * x * y * y * z * z * z * z - std::sqrt(0.341796875) * x * y * y * y * y * y * y * y * y + std::sqrt(66.9921875) * x * y * y * y * y * y * y * z * z - std::sqrt(787.5) * x * y * y * y * y * z * z * z * z) + e_1 * (std::sqrt(3.076171875) * x * x * x * x * x * x * x + std::sqrt(889.013671875) * x * x * x * x * x * y * y - std::sqrt(602.9296875) * x * x * x * x * x * z * z + std::sqrt(467.919921875) * x * x * x * y * y * y * y + std::sqrt(26036.71875) * x * x * x * y * y * z * z + std::sqrt(7087.5) * x * x * x * z * z * z * z - std::sqrt(41.357421875) * x * y * y * y * y * y * y - std::sqrt(8970.1171875) * x * y * y * y * y * z * z + std::sqrt(7087.5) * x * y * y * z * z * z * z) + e_2 * (std::sqrt(442.96875) * x * x * x * x * x + std::sqrt(86821.875) * x * x * x * y * y + std::sqrt(44296.875) * x * x * x * z * z - std::sqrt(3986.71875) * x * y * y * y * y + std::sqrt(44296.875) * x * y * y * z * z + std::sqrt(28350.0) * x * z * z * z * z) + e_3 * (std::sqrt(56896.875) * x * x * x + std::sqrt(56896.875) * x * y * y + std::sqrt(381150.0) * x * z * z) + e_4 * (std::sqrt(154350.0) * x);

        pc_21[k] = e_0 * (std::sqrt(13.8427734375) * x * x * x * x * x * x * y * y * z + std::sqrt(38.4521484375) * x * x * x * x * y * y * y * y * z - std::sqrt(1205.859375) * x * x * x * x * y * y * z * z * z + std::sqrt(1.5380859375) * x * x * y * y * y * y * y * y * z - std::sqrt(535.9375) * x * x * y * y * y * y * z * z * z + std::sqrt(1575.0) * x * x * y * y * z * z * z * z * z - std::sqrt(1.5380859375) * y * y * y * y * y * y * y * y * z + std::sqrt(133.984375) * y * y * y * y * y * y * z * z * z - std::sqrt(175.0) * y * y * y * y * z * z * z * z * z) + e_1 * (std::sqrt(13.8427734375) * x * x * x * x * x * x * z - std::sqrt(124.5849609375) * x * x * x * x * y * y * z - std::sqrt(1205.859375) * x * x * x * x * z * z * z - std::sqrt(186.1083984375) * x * x * y * y * y * y * z + std::sqrt(2460.9375) * x * x * y * y * z * z * z + std::sqrt(1575.0) * x * x * z * z * z * z * z + std::sqrt(1.5380859375) * y * y * y * y * y * y * z + std::sqrt(330.859375) * y * y * y * y * z * z * z - std::sqrt(1575.0) * y * y * z * z * z * z * z) + e_2 * (-std::sqrt(885.9375) * x * x * x * x * z + std::sqrt(14175.0) * x * x * z * z * z + std::sqrt(885.9375) * y * y * y * y * z - std::sqrt(14175.0) * y * y * z * z * z) + e_3 * (std::sqrt(3543.75) * x * x * z - std::sqrt(3543.75) * y * y * z);

        pc_22[k] = e_0 * (-std::sqrt(0.3460693359375) * x * x * x * x * x * x * x * x * y - std::sqrt(2.4609375) * x * x * x * x * x * x * y * y * y + std::sqrt(88.59375) * x * x * x * x * x * x * y * z * z - std::sqrt(1.38427734375) * x * x * x * x * y * y * y * y * y + std::sqrt(246.09375) * x * x * x * x * y * y * y * z * z - std::sqrt(1538.0859375) * x * x * x * x * y * z * z * z * z + std::sqrt(9.84375) * x * x * y * y * y * y * y * z * z - std::sqrt(683.59375) * x * x * y * y * y * z * z * z * z + std::sqrt(157.5) * x * x * y * z * z * z * z * z * z + std::sqrt(0.0384521484375) * y * y * y * y * y * y * y * y * y - std::sqrt(9.84375) * y * y * y * y * y * y * y * z * z + std::sqrt(170.8984375) * y * y * y * y * y * z * z * z * z - std::sqrt(17.5) * y * y * y * z * z * z * z * z * z) + e_1 * (-std::sqrt(138.427734375) * x * x * x * x * x * x * y - std::sqrt(384.521484375) * x * x * x * x * y * y * y - std::sqrt(2214.84375) * x * x * x * x * y * z * z - std::sqrt(15.380859375) * x * x * y * y * y * y * y - std::sqrt(984.375) * x * x * y * y * y * z * z - std::sqrt(15750.0) * x * x * y * z * z * z * z + std::sqrt(15.380859375) * y * y * y * y * y * y * y + std::sqrt(246.09375) * y * y * y * y * y * z * z + std::sqrt(1750.0) * y * y * y * z * z * z * z) + e_2 * (-std::sqrt(19933.59375) * x * x * x * x * y - std::sqrt(8859.375) * x * x * y * y * y - std::sqrt(318937.5) * x * x * y * z * z + std::sqrt(2214.84375) * y * y * y * y * y + std::sqrt(35437.5) * y * y * y * z * z) + e_3 * (-std::sqrt(318937.5) * x * x * y + std::sqrt(35437.5) * y * y * y);

        pc_23[k] = e_0 * (std::sqrt(13.8427734375) * x * x * x * x * x * x * x * y * z + std::sqrt(38.4521484375) * x * x * x * x * x * y * y * y * z - std::sqrt(1205.859375) * x * x * x * x * x * y * z * z * z + std::sqrt(1.5380859375) * x * x * x * y * y * y * y * y * z - std::sqrt(535.9375) * x * x * x * y * y * y * z * z * z + std::sqrt(1575.0) * x * x * x * y * z * z * z * z * z - std::sqrt(1.5380859375) * x * y * y * y * y * y * y * y * z + std::sqrt(133.984375) * x * y * y * y * y * y * z * z * z - std::sqrt(175.0) * x * y * y * y * z * z * z * z * z) + e_1 * (-std::sqrt(55.37109375) * x * x * x * x * x * y * z + std::sqrt(24.609375) * x * x * x * y * y * y * z - std::sqrt(393.75) * x * x * x * y * z * z * z + std::sqrt(153.80859375) * x * y * y * y * y * y * z - std::sqrt(7393.75) * x * y * y * y * z * z * z + std::sqrt(6300.0) * x * y * z * z * z * z * z) + e_2 * (-std::sqrt(3543.75) * x * x * x * y * z - std::sqrt(3543.75) * x * y * y * y * z + std::sqrt(56700.0) * x * y * z * z * z) + e_3 * (std::sqrt(14175.0) * x * y * z);
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

        pc_24[k] = e_0 * (std::sqrt(0.76904296875) * x * x * x * x * x * x * x * x * y + std::sqrt(0.341796875) * x * x * x * x * x * x * y * y * y - std::sqrt(150.732421875) * x * x * x * x * x * x * y * z * z - std::sqrt(1.3671875) * x * x * x * x * y * y * y * y * y + std::sqrt(16.748046875) * x * x * x * x * y * y * y * z * z + std::sqrt(1771.875) * x * x * x * x * y * z * z * z * z - std::sqrt(0.341796875) * x * x * y * y * y * y * y * y * y + std::sqrt(150.732421875) * x * x * y * y * y * y * y * z * z - std::sqrt(3150.0) * x * x * y * y * y * z * z * z * z + std::sqrt(0.08544921875) * y * y * y * y * y * y * y * y * y - std::sqrt(16.748046875) * y * y * y * y * y * y * y * z * z + std::sqrt(196.875) * y * y * y * y * y * z * z * z * z) + e_1 * (std::sqrt(196.875) * x * x * x * x * x * x * y + std::sqrt(1.3671875) * x * x * x * x * y * y * y + std::sqrt(6509.1796875) * x * x * x * x * y * z * z - std::sqrt(49.21875) * x * x * y * y * y * y * y - std::sqrt(35880.46875) * x * x * y * y * y * z * z + std::sqrt(7087.5) * x * x * y * z * z * z * z + std::sqrt(34.1796875) * y * y * y * y * y * y * y + std::sqrt(110.7421875) * y * y * y * y * y * z * z + std::sqrt(7087.5) * y * y * y * z * z * z * z) + e_2 * (std::sqrt(21705.46875) * x * x * x * x * y - std::sqrt(15946.875) * x * x * y * y * y + std::sqrt(44296.875) * x * x * y * z * z + std::sqrt(3986.71875) * y * y * y * y * y + std::sqrt(44296.875) * y * y * y * z * z + std::sqrt(28350.0) * y * z * z * z * z) + e_3 * (std::sqrt(56896.875) * x * x * y + std::sqrt(56896.875) * y * y * y + std::sqrt(381150.0) * y * z * z) + e_4 * (std::sqrt(154350.0) * y);

        pc_25[k] = e_0 * (-3.28125 * x * x * x * x * x * x * x * y * z + 7.65625 * x * x * x * x * x * y * y * y * z + 26.25 * x * x * x * x * x * y * z * z * z + 7.65625 * x * x * x * y * y * y * y * y * z - 87.5 * x * x * x * y * y * y * z * z * z - 3.28125 * x * y * y * y * y * y * y * y * z + 26.25 * x * y * y * y * y * y * z * z * z) + e_1 * (32.8125 * x * x * x * x * x * y * z - 109.375 * x * x * x * y * y * y * z + 32.8125 * x * y * y * y * y * y * z);

        pc_26[k] = e_0 * (-std::sqrt(1.3458251953125) * x * x * x * x * x * x * x * x * y + std::sqrt(38.28125) * x * x * x * x * x * x * y * y * y + std::sqrt(86.1328125) * x * x * x * x * x * x * y * z * z + std::sqrt(14.95361328125) * x * x * x * x * y * y * y * y * y - std::sqrt(3454.8828125) * x * x * x * x * y * y * y * z * z - std::sqrt(9.5703125) * x * x * y * y * y * y * y * y * y + std::sqrt(775.1953125) * x * x * y * y * y * y * y * z * z + std::sqrt(0.1495361328125) * y * y * y * y * y * y * y * y * y - std::sqrt(9.5703125) * y * y * y * y * y * y * y * z * z) + e_1 * (-std::sqrt(21.533203125) * x * x * x * x * x * x * y + std::sqrt(5285.205078125) * x * x * x * x * y * y * y - std::sqrt(1378.125) * x * x * x * x * y * z * z - std::sqrt(193.798828125) * x * x * y * y * y * y * y - std::sqrt(5512.5) * x * x * y * y * y * z * z + std::sqrt(59.814453125) * y * y * y * y * y * y * y - std::sqrt(1378.125) * y * y * y * y * y * z * z) + e_2 * (std::sqrt(3100.78125) * x * x * x * x * y + std::sqrt(12403.125) * x * x * y * y * y - std::sqrt(49612.5) * x * x * y * z * z + std::sqrt(3100.78125) * y * y * y * y * y - std::sqrt(49612.5) * y * y * y * z * z) + e_3 * (std::sqrt(22050.0) * x * x * y + std::sqrt(22050.0) * y * y * y - std::sqrt(88200.0) * y * z * z) + e_4 * (std::sqrt(5512.5) * y);

        pc_27[k] = e_0 * (-std::sqrt(229.6875) * x * x * x * x * x * x * y * y * z + std::sqrt(918.75) * x * x * x * x * y * y * z * z * z + std::sqrt(229.6875) * x * x * y * y * y * y * y * y * z - std::sqrt(918.75) * x * x * y * y * y * y * z * z * z) + e_1 * (-std::sqrt(229.6875) * x * x * x * x * x * x * z - std::sqrt(18604.6875) * x * x * x * x * y * y * z + std::sqrt(918.75) * x * x * x * x * z * z * z + std::sqrt(18604.6875) * x * x * y * y * y * y * z + std::sqrt(229.6875) * y * y * y * y * y * y * z - std::sqrt(918.75) * y * y * y * y * z * z * z) + e_2 * (-std::sqrt(18604.6875) * x * x * x * x * z + std::sqrt(8268.75) * x * x * z * z * z + std::sqrt(18604.6875) * y * y * y * y * z - std::sqrt(8268.75) * y * y * z * z * z) + e_3 * (-std::sqrt(33075.0) * x * x * z + std::sqrt(33075.0) * y * y * z);
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

        pc_28[k] = e_0 * (-std::sqrt(1033.59375) * x * x * x * x * x * y * y * z * z - std::sqrt(459.375) * x * x * x * y * y * y * y * z * z + std::sqrt(4134.375) * x * x * x * y * y * z * z * z * z + std::sqrt(114.84375) * x * y * y * y * y * y * y * z * z - std::sqrt(459.375) * x * y * y * y * y * z * z * z * z) + e_1 * (-std::sqrt(1033.59375) * x * x * x * x * x * y * y - std::sqrt(1033.59375) * x * x * x * x * x * z * z - std::sqrt(459.375) * x * x * x * y * y * y * y - std::sqrt(4134.375) * x * x * x * y * y * z * z + std::sqrt(4134.375) * x * x * x * z * z * z * z + std::sqrt(114.84375) * x * y * y * y * y * y * y - std::sqrt(1033.59375) * x * y * y * y * y * z * z + std::sqrt(4134.375) * x * y * y * z * z * z * z) + e_2 * (-std::sqrt(1033.59375) * x * x * x * x * x - std::sqrt(66150.0) * x * x * x * y * y + std::sqrt(1033.59375) * x * y * y * y * y + std::sqrt(16537.5) * x * z * z * z * z) + e_3 * (-std::sqrt(37209.375) * x * x * x - std::sqrt(37209.375) * x * y * y + std::sqrt(66150.0) * x * z * z) + e_4 * (-std::sqrt(16537.5) * x);

        pc_29[k] = e_0 * (std::sqrt(32.8125) * x * x * x * x * x * x * y * y * z + std::sqrt(131.25) * x * x * x * x * y * y * y * y * z - std::sqrt(2100.0) * x * x * x * x * y * y * z * z * z + std::sqrt(32.8125) * x * x * y * y * y * y * y * y * z - std::sqrt(2100.0) * x * x * y * y * y * y * z * z * z + std::sqrt(4725.0) * x * x * y * y * z * z * z * z * z) + e_1 * (std::sqrt(32.8125) * x * x * x * x * x * x * z + std::sqrt(295.3125) * x * x * x * x * y * y * z - std::sqrt(2100.0) * x * x * x * x * z * z * z + std::sqrt(295.3125) * x * x * y * y * y * y * z + std::sqrt(18900.0) * x * x * y * y * z * z * z + std::sqrt(4725.0) * x * x * z * z * z * z * z + std::sqrt(32.8125) * y * y * y * y * y * y * z - std::sqrt(2100.0) * y * y * y * y * z * z * z + std::sqrt(4725.0) * y * y * z * z * z * z * z) + e_2 * (-std::sqrt(295.3125) * x * x * x * x * z + std::sqrt(95681.25) * x * x * y * y * z + std::sqrt(75600.0) * x * x * z * z * z - std::sqrt(295.3125) * y * y * y * y * z + std::sqrt(75600.0) * y * y * z * z * z + std::sqrt(4725.0) * z * z * z * z * z) + e_3 * (std::sqrt(118125.0) * x * x * z + std::sqrt(118125.0) * y * y * z + std::sqrt(170100.0) * z * z * z) + e_4 * (std::sqrt(231525.0) * z);

        pc_30[k] = e_0 * (std::sqrt(147.65625) * x * x * x * x * x * y * y * z * z + std::sqrt(590.625) * x * x * x * y * y * y * y * z * z - std::sqrt(1640.625) * x * x * x * y * y * z * z * z * z + std::sqrt(147.65625) * x * y * y * y * y * y * y * z * z - std::sqrt(1640.625) * x * y * y * y * y * z * z * z * z + std::sqrt(1050.0) * x * y * y * z * z * z * z * z * z) + e_1 * (std::sqrt(147.65625) * x * x * x * x * x * y * y + std::sqrt(147.65625) * x * x * x * x * x * z * z + std::sqrt(590.625) * x * x * x * y * y * y * y + std::sqrt(590.625) * x * x * x * y * y * z * z - std::sqrt(1640.625) * x * x * x * z * z * z * z + std::sqrt(147.65625) * x * y * y * y * y * y * y + std::sqrt(147.65625) * x * y * y * y * y * z * z + std::sqrt(14765.625) * x * y * y * z * z * z * z + std::sqrt(1050.0) * x * z * z * z * z * z * z) + e_2 * (std::sqrt(147.65625) * x * x * x * x * x + std::sqrt(21262.5) * x * x * x * y * y - std::sqrt(2362.5) * x * x * x * z * z + std::sqrt(17866.40625) * x * y * y * y * y + std::sqrt(191362.5) * x * y * y * z * z + std::sqrt(59062.5) * x * z * z * z * z) + e_3 * (std::sqrt(5315.625) * x * x * x + std::sqrt(312440.625) * x * y * y + std::sqrt(340200.0) * x * z * z) + e_4 * (std::sqrt(115762.5) * x);

        pc_31[k] = e_0 * (-std::sqrt(3.69140625) * x * x * x * x * x * x * x * y * z - std::sqrt(33.22265625) * x * x * x * x * x * y * y * y * z + std::sqrt(369.140625) * x * x * x * x * x * y * z * z * z - std::sqrt(33.22265625) * x * x * x * y * y * y * y * y * z + std::sqrt(1476.5625) * x * x * x * y * y * y * z * z * z - std::sqrt(1286.25) * x * x * x * y * z * z * z * z * z - std::sqrt(3.69140625) * x * y * y * y * y * y * y * y * z + std::sqrt(369.140625) * x * y * y * y * y * y * z * z * z - std::sqrt(1286.25) * x * y * y * y * z * z * z * z * z + std::sqrt(105.0) * x * y * z * z * z * z * z * z * z) + e_1 * (-std::sqrt(2625.0) * x * x * x * y * z * z * z - std::sqrt(2625.0) * x * y * y * y * z * z * z) + e_2 * (-std::sqrt(5906.25) * x * x * x * y * z - std::sqrt(5906.25) * x * y * y * y * z - std::sqrt(23625.0) * x * y * z * z * z) + e_3 * (-std::sqrt(94500.0) * x * y * z);
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

        pc_32[k] = e_0 * (std::sqrt(147.65625) * x * x * x * x * x * x * y * z * z + std::sqrt(590.625) * x * x * x * x * y * y * y * z * z - std::sqrt(1640.625) * x * x * x * x * y * z * z * z * z + std::sqrt(147.65625) * x * x * y * y * y * y * y * z * z - std::sqrt(1640.625) * x * x * y * y * y * z * z * z * z + std::sqrt(1050.0) * x * x * y * z * z * z * z * z * z) + e_1 * (std::sqrt(147.65625) * x * x * x * x * x * x * y + std::sqrt(590.625) * x * x * x * x * y * y * y + std::sqrt(147.65625) * x * x * x * x * y * z * z + std::sqrt(147.65625) * x * x * y * y * y * y * y + std::sqrt(590.625) * x * x * y * y * y * z * z + std::sqrt(14765.625) * x * x * y * z * z * z * z + std::sqrt(147.65625) * y * y * y * y * y * z * z - std::sqrt(1640.625) * y * y * y * z * z * z * z + std::sqrt(1050.0) * y * z * z * z * z * z * z) + e_2 * (std::sqrt(17866.40625) * x * x * x * x * y + std::sqrt(21262.5) * x * x * y * y * y + std::sqrt(191362.5) * x * x * y * z * z + std::sqrt(147.65625) * y * y * y * y * y - std::sqrt(2362.5) * y * y * y * z * z + std::sqrt(59062.5) * y * z * z * z * z) + e_3 * (std::sqrt(312440.625) * x * x * y + std::sqrt(5315.625) * y * y * y + std::sqrt(340200.0) * y * z * z) + e_4 * (std::sqrt(115762.5) * y);

        pc_33[k] = e_0 * (std::sqrt(8.203125) * x * x * x * x * x * x * x * y * z + std::sqrt(8.203125) * x * x * x * x * x * y * y * y * z - std::sqrt(525.0) * x * x * x * x * x * y * z * z * z - std::sqrt(8.203125) * x * x * x * y * y * y * y * y * z + std::sqrt(1181.25) * x * x * x * y * z * z * z * z * z - std::sqrt(8.203125) * x * y * y * y * y * y * y * y * z + std::sqrt(525.0) * x * y * y * y * y * y * z * z * z - std::sqrt(1181.25) * x * y * y * y * z * z * z * z * z) + e_1 * (std::sqrt(13125.0) * x * x * x * y * z * z * z - std::sqrt(13125.0) * x * y * y * y * z * z * z) + e_2 * (std::sqrt(29531.25) * x * x * x * y * z - std::sqrt(29531.25) * x * y * y * y * z);

        pc_34[k] = e_0 * (-std::sqrt(114.84375) * x * x * x * x * x * x * y * z * z + std::sqrt(459.375) * x * x * x * x * y * y * y * z * z + std::sqrt(459.375) * x * x * x * x * y * z * z * z * z + std::sqrt(1033.59375) * x * x * y * y * y * y * y * z * z - std::sqrt(4134.375) * x * x * y * y * y * z * z * z * z) + e_1 * (-std::sqrt(114.84375) * x * x * x * x * x * x * y + std::sqrt(459.375) * x * x * x * x * y * y * y + std::sqrt(1033.59375) * x * x * x * x * y * z * z + std::sqrt(1033.59375) * x * x * y * y * y * y * y + std::sqrt(4134.375) * x * x * y * y * y * z * z - std::sqrt(4134.375) * x * x * y * z * z * z * z + std::sqrt(1033.59375) * y * y * y * y * y * z * z - std::sqrt(4134.375) * y * y * y * z * z * z * z) + e_2 * (-std::sqrt(1033.59375) * x * x * x * x * y + std::sqrt(66150.0) * x * x * y * y * y + std::sqrt(1033.59375) * y * y * y * y * y - std::sqrt(16537.5) * y * z * z * z * z) + e_3 * (std::sqrt(37209.375) * x * x * y + std::sqrt(37209.375) * y * y * y - std::sqrt(66150.0) * y * z * z) + e_4 * (std::sqrt(16537.5) * y);

        pc_35[k] = e_0 * (-std::sqrt(14.35546875) * x * x * x * x * x * x * x * y * z + std::sqrt(358.88671875) * x * x * x * x * x * y * y * y * z + std::sqrt(57.421875) * x * x * x * x * x * y * z * z * z + std::sqrt(358.88671875) * x * x * x * y * y * y * y * y * z - std::sqrt(2067.1875) * x * x * x * y * y * y * z * z * z - std::sqrt(14.35546875) * x * y * y * y * y * y * y * y * z + std::sqrt(57.421875) * x * y * y * y * y * y * z * z * z) + e_1 * (std::sqrt(58800.0) * x * x * x * y * y * y * z - std::sqrt(3675.0) * x * x * x * y * z * z * z - std::sqrt(3675.0) * x * y * y * y * z * z * z) + e_2 * (std::sqrt(74418.75) * x * x * x * y * z + std::sqrt(74418.75) * x * y * y * y * z - std::sqrt(33075.0) * x * y * z * z * z) + e_3 * (std::sqrt(132300.0) * x * y * z);
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

        pc_36[k] = e_0 * (std::sqrt(2.05078125) * x * x * x * x * x * x * x * y * y + std::sqrt(2.05078125) * x * x * x * x * x * y * y * y * y - std::sqrt(295.3125) * x * x * x * x * x * y * y * z * z - std::sqrt(2.05078125) * x * x * x * y * y * y * y * y * y + std::sqrt(131.25) * x * x * x * y * y * z * z * z * z - std::sqrt(2.05078125) * x * y * y * y * y * y * y * y * y + std::sqrt(295.3125) * x * y * y * y * y * y * y * z * z - std::sqrt(131.25) * x * y * y * y * y * z * z * z * z) + e_1 * (std::sqrt(2.05078125) * x * x * x * x * x * x * x + std::sqrt(461.42578125) * x * x * x * x * x * y * y - std::sqrt(295.3125) * x * x * x * x * x * z * z - std::sqrt(51.26953125) * x * x * x * y * y * y * y - std::sqrt(10631.25) * x * x * x * y * y * z * z + std::sqrt(131.25) * x * x * x * z * z * z * z - std::sqrt(740.33203125) * x * y * y * y * y * y * y + std::sqrt(35732.8125) * x * y * y * y * y * z * z - std::sqrt(1181.25) * x * y * y * z * z * z * z) + e_2 * (std::sqrt(295.3125) * x * x * x * x * x + std::sqrt(1181.25) * x * x * x * y * y - std::sqrt(10631.25) * x * x * x * z * z - std::sqrt(14470.3125) * x * y * y * y * y + std::sqrt(95681.25) * x * y * y * z * z) + e_3 * (std::sqrt(1181.25) * x * x * x - std::sqrt(10631.25) * x * y * y);

        pc_37[k] = e_0 * (std::sqrt(9.228515625) * x * x * x * x * x * x * y * y * z + std::sqrt(25.634765625) * x * x * x * x * y * y * y * y * z - std::sqrt(1328.90625) * x * x * x * x * y * y * z * z * z + std::sqrt(1.025390625) * x * x * y * y * y * y * y * y * z - std::sqrt(590.625) * x * x * y * y * y * y * z * z * z + std::sqrt(590.625) * x * x * y * y * z * z * z * z * z - std::sqrt(1.025390625) * y * y * y * y * y * y * y * y * z + std::sqrt(147.65625) * y * y * y * y * y * y * z * z * z - std::sqrt(65.625) * y * y * y * y * z * z * z * z * z) + e_1 * (std::sqrt(9.228515625) * x * x * x * x * x * x * z - std::sqrt(1116.650390625) * x * x * x * x * y * y * z - std::sqrt(1328.90625) * x * x * x * x * z * z * z - std::sqrt(747.509765625) * x * x * y * y * y * y * z - std::sqrt(14765.625) * x * x * y * y * z * z * z + std::sqrt(590.625) * x * x * z * z * z * z * z + std::sqrt(83.056640625) * y * y * y * y * y * y * z + std::sqrt(5922.65625) * y * y * y * y * z * z * z - std::sqrt(590.625) * y * y * z * z * z * z * z) + e_2 * (-std::sqrt(2362.5) * x * x * x * x * z - std::sqrt(132890.625) * x * x * y * y * z - std::sqrt(2362.5) * x * x * z * z * z + std::sqrt(28940.625) * y * y * y * y * z + std::sqrt(2362.5) * y * y * z * z * z) + e_3 * (-std::sqrt(71465.625) * x * x * z + std::sqrt(71465.625) * y * y * z);

        pc_38[k] = e_0 * (-std::sqrt(0.29296875) * x * x * x * x * x * x * x * y * y - std::sqrt(2.63671875) * x * x * x * x * x * y * y * y * y + std::sqrt(94.921875) * x * x * x * x * x * y * y * z * z - std::sqrt(2.63671875) * x * x * x * y * y * y * y * y * y + std::sqrt(379.6875) * x * x * x * y * y * y * y * z * z - std::sqrt(1875.0) * x * x * x * y * y * z * z * z * z - std::sqrt(0.29296875) * x * y * y * y * y * y * y * y * y + std::sqrt(94.921875) * x * y * y * y * y * y * y * z * z - std::sqrt(1875.0) * x * y * y * y * y * z * z * z * z + std::sqrt(675.0) * x * y * y * z * z * z * z * z * z) + e_1 * (-std::sqrt(0.29296875) * x * x * x * x * x * x * x - std::sqrt(129.19921875) * x * x * x * x * x * y * y + std::sqrt(94.921875) * x * x * x * x * x * z * z - std::sqrt(445.60546875) * x * x * x * y * y * y * y - std::sqrt(2067.1875) * x * x * x * y * y * z * z - std::sqrt(1875.0) * x * x * x * z * z * z * z - std::sqrt(105.76171875) * x * y * y * y * y * y * y - std::sqrt(3048.046875) * x * y * y * y * y * z * z + std::sqrt(675.0) * x * z * z * z * z * z * z) + e_2 * (-std::sqrt(42.1875) * x * x * x * x * x - std::sqrt(20418.75) * x * x * x * y * y - std::sqrt(10800.0) * x * x * x * z * z - std::sqrt(18604.6875) * x * y * y * y * y - std::sqrt(54675.0) * x * y * y * z * z + std::sqrt(16875.0) * x * z * z * z * z) + e_3 * (-std::sqrt(10800.0) * x * x * x - std::sqrt(243675.0) * x * y * y + std::sqrt(6075.0) * x * z * z) + e_4 * (-std::sqrt(33075.0) * x);

        pc_39[k] = e_0 * (-std::sqrt(1.318359375) * x * x * x * x * x * x * y * y * z - std::sqrt(11.865234375) * x * x * x * x * y * y * y * y * z + std::sqrt(234.375) * x * x * x * x * y * y * z * z * z - std::sqrt(11.865234375) * x * x * y * y * y * y * y * y * z + std::sqrt(937.5) * x * x * y * y * y * y * z * z * z - std::sqrt(759.375) * x * x * y * y * z * z * z * z * z - std::sqrt(1.318359375) * y * y * y * y * y * y * y * y * z + std::sqrt(234.375) * y * y * y * y * y * y * z * z * z - std::sqrt(759.375) * y * y * y * y * z * z * z * z * z + std::sqrt(150.0) * y * y * z * z * z * z * z * z * z) + e_1 * (-std::sqrt(1.318359375) * x * x * x * x * x * x * z + std::sqrt(64.599609375) * x * x * x * x * y * y * z + std::sqrt(234.375) * x * x * x * x * z * z * z + std::sqrt(381.005859375) * x * x * y * y * y * y * z - std::sqrt(759.375) * x * x * z * z * z * z * z + std::sqrt(106.787109375) * y * y * y * y * y * y * z - std::sqrt(234.375) * y * y * y * y * z * z * z + std::sqrt(4134.375) * y * y * z * z * z * z * z + std::sqrt(150.0) * z * z * z * z * z * z * z) + e_2 * (std::sqrt(337.5) * x * x * x * x * z + std::sqrt(6834.375) * x * x * y * y * z - std::sqrt(8437.5) * x * x * z * z * z + std::sqrt(4134.375) * y * y * y * y * z + std::sqrt(75937.5) * y * y * z * z * z + std::sqrt(21600.0) * z * z * z * z * z) + e_3 * (-std::sqrt(759.375) * x * x * z + std::sqrt(186384.375) * y * y * z + std::sqrt(303750.0) * z * z * z) + e_4 * (std::sqrt(264600.0) * z);
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

        pc_40[k] = e_0 * (std::sqrt(0.032958984375) * x * x * x * x * x * x * x * x * y + std::sqrt(0.52734375) * x * x * x * x * x * x * y * y * y - std::sqrt(13.18359375) * x * x * x * x * x * x * y * z * z + std::sqrt(1.1865234375) * x * x * x * x * y * y * y * y * y - std::sqrt(118.65234375) * x * x * x * x * y * y * y * z * z + std::sqrt(375.0) * x * x * x * x * y * z * z * z * z + std::sqrt(0.52734375) * x * x * y * y * y * y * y * y * y - std::sqrt(118.65234375) * x * x * y * y * y * y * y * z * z + std::sqrt(1500.0) * x * x * y * y * y * z * z * z * z - std::sqrt(303.75) * x * x * y * z * z * z * z * z * z + std::sqrt(0.032958984375) * y * y * y * y * y * y * y * y * y - std::sqrt(13.18359375) * y * y * y * y * y * y * y * z * z + std::sqrt(375.0) * y * y * y * y * y * z * z * z * z - std::sqrt(303.75) * y * y * y * z * z * z * z * z * z + std::sqrt(15.0) * y * z * z * z * z * z * z * z * z) + e_1 * (std::sqrt(13.18359375) * x * x * x * x * x * x * y + std::sqrt(118.65234375) * x * x * x * x * y * y * y + std::sqrt(843.75) * x * x * x * x * y * z * z + std::sqrt(118.65234375) * x * x * y * y * y * y * y + std::sqrt(3375.0) * x * x * y * y * y * z * z - std::sqrt(843.75) * x * x * y * z * z * z * z + std::sqrt(13.18359375) * y * y * y * y * y * y * y + std::sqrt(843.75) * y * y * y * y * y * z * z - std::sqrt(843.75) * y * y * y * z * z * z * z + std::sqrt(1500.0) * y * z * z * z * z * z * z) + e_2 * (std::sqrt(3375.0) * x * x * x * x * y + std::sqrt(13500.0) * x * x * y * y * y + std::sqrt(7593.75) * x * x * y * z * z + std::sqrt(3375.0) * y * y * y * y * y + std::sqrt(7593.75) * y * y * y * z * z + std::sqrt(54000.0) * y * z * z * z * z) + e_3 * (std::sqrt(68343.75) * x * x * y + std::sqrt(68343.75) * y * y * y + std::sqrt(337500.0) * y * z * z) + e_4 * (std::sqrt(165375.0) * y);

        pc_41[k] = e_0 * (-std::sqrt(1.318359375) * x * x * x * x * x * x * x * y * z - std::sqrt(11.865234375) * x * x * x * x * x * y * y * y * z + std::sqrt(234.375) * x * x * x * x * x * y * z * z * z - std::sqrt(11.865234375) * x * x * x * y * y * y * y * y * z + std::sqrt(937.5) * x * x * x * y * y * y * z * z * z - std::sqrt(759.375) * x * x * x * y * z * z * z * z * z - std::sqrt(1.318359375) * x * y * y * y * y * y * y * y * z + std::sqrt(234.375) * x * y * y * y * y * y * z * z * z - std::sqrt(759.375) * x * y * y * y * z * z * z * z * z + std::sqrt(150.0) * x * y * z * z * z * z * z * z * z) + e_1 * (std::sqrt(131.8359375) * x * x * x * x * x * y * z + std::sqrt(527.34375) * x * x * x * y * y * y * z - std::sqrt(937.5) * x * x * x * y * z * z * z + std::sqrt(131.8359375) * x * y * y * y * y * y * z - std::sqrt(937.5) * x * y * y * y * z * z * z + std::sqrt(8437.5) * x * y * z * z * z * z * z) + e_2 * (std::sqrt(2109.375) * x * x * x * y * z + std::sqrt(2109.375) * x * y * y * y * z + std::sqrt(135000.0) * x * y * z * z * z) + e_3 * (std::sqrt(210937.5) * x * y * z);

        pc_42[k] = e_0 * (-std::sqrt(0.0732421875) * x * x * x * x * x * x * x * x * y - std::sqrt(0.29296875) * x * x * x * x * x * x * y * y * y + std::sqrt(23.73046875) * x * x * x * x * x * x * y * z * z + std::sqrt(23.73046875) * x * x * x * x * y * y * y * z * z - std::sqrt(468.75) * x * x * x * x * y * z * z * z * z + std::sqrt(0.29296875) * x * x * y * y * y * y * y * y * y - std::sqrt(23.73046875) * x * x * y * y * y * y * y * z * z + std::sqrt(168.75) * x * x * y * z * z * z * z * z * z + std::sqrt(0.0732421875) * y * y * y * y * y * y * y * y * y - std::sqrt(23.73046875) * y * y * y * y * y * y * y * z * z + std::sqrt(468.75) * y * y * y * y * y * z * z * z * z - std::sqrt(168.75) * y * y * y * z * z * z * z * z * z) + e_1 * (-std::sqrt(18.75) * x * x * x * x * x * x * y - std::sqrt(10.546875) * x * x * x * x * y * y * y - std::sqrt(1782.421875) * x * x * x * x * y * z * z + std::sqrt(42.1875) * x * x * y * y * y * y * y - std::sqrt(379.6875) * x * x * y * y * y * z * z + std::sqrt(4218.75) * x * x * y * z * z * z * z + std::sqrt(29.296875) * y * y * y * y * y * y * y + std::sqrt(516.796875) * y * y * y * y * y * z * z + std::sqrt(468.75) * y * y * y * z * z * z * z - std::sqrt(675.0) * y * z * z * z * z * z * z) + e_2 * (-std::sqrt(3417.1875) * x * x * x * x * y + std::sqrt(168.75) * x * x * y * y * y + std::sqrt(1518.75) * x * x * y * z * z + std::sqrt(5104.6875) * y * y * y * y * y + std::sqrt(28518.75) * y * y * y * z * z - std::sqrt(16875.0) * y * z * z * z * z) + e_3 * (-std::sqrt(8268.75) * x * x * y + std::sqrt(89268.75) * y * y * y - std::sqrt(6075.0) * y * z * z) + e_4 * (std::sqrt(33075.0) * y);

        pc_43[k] = e_0 * (std::sqrt(1.025390625) * x * x * x * x * x * x * x * y * z - std::sqrt(1.025390625) * x * x * x * x * x * y * y * y * z - std::sqrt(147.65625) * x * x * x * x * x * y * z * z * z - std::sqrt(25.634765625) * x * x * x * y * y * y * y * y * z + std::sqrt(590.625) * x * x * x * y * y * y * z * z * z + std::sqrt(65.625) * x * x * x * y * z * z * z * z * z - std::sqrt(9.228515625) * x * y * y * y * y * y * y * y * z + std::sqrt(1328.90625) * x * y * y * y * y * y * z * z * z - std::sqrt(590.625) * x * y * y * y * z * z * z * z * z) + e_1 * (-std::sqrt(332.2265625) * x * x * x * x * x * y * z + std::sqrt(147.65625) * x * x * x * y * y * y * z + std::sqrt(1050.0) * x * x * x * y * z * z * z + std::sqrt(922.8515625) * x * y * y * y * y * y * z + std::sqrt(37800.0) * x * y * y * y * z * z * z - std::sqrt(2362.5) * x * y * z * z * z * z * z) + e_2 * (-std::sqrt(590.625) * x * x * x * y * z + std::sqrt(213215.625) * x * y * y * y * z + std::sqrt(9450.0) * x * y * z * z * z) + e_3 * (std::sqrt(285862.5) * x * y * z);
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

        pc_44[k] = e_0 * (std::sqrt(0.128173828125) * x * x * x * x * x * x * x * x * y - std::sqrt(2.05078125) * x * x * x * x * x * x * y * y * y - std::sqrt(18.45703125) * x * x * x * x * x * x * y * z * z - std::sqrt(12.8173828125) * x * x * x * x * y * y * y * y * y + std::sqrt(461.42578125) * x * x * x * x * y * y * y * z * z + std::sqrt(8.203125) * x * x * x * x * y * z * z * z * z - std::sqrt(2.05078125) * x * x * y * y * y * y * y * y * y + std::sqrt(461.42578125) * x * x * y * y * y * y * y * z * z - std::sqrt(295.3125) * x * x * y * y * y * z * z * z * z + std::sqrt(0.128173828125) * y * y * y * y * y * y * y * y * y - std::sqrt(18.45703125) * y * y * y * y * y * y * y * z * z + std::sqrt(8.203125) * y * y * y * y * y * z * z * z * z) + e_1 * (std::sqrt(2.05078125) * x * x * x * x * x * x * y - std::sqrt(1281.73828125) * x * x * x * x * y * y * y + std::sqrt(295.3125) * x * x * x * x * y * z * z - std::sqrt(904.39453125) * x * x * y * y * y * y * y + std::sqrt(57881.25) * x * x * y * y * y * z * z - std::sqrt(1181.25) * x * x * y * z * z * z * z + std::sqrt(51.26953125) * y * y * y * y * y * y * y - std::sqrt(2657.8125) * y * y * y * y * y * z * z + std::sqrt(131.25) * y * y * y * z * z * z * z) + e_2 * (-std::sqrt(1181.25) * x * x * x * x * y - std::sqrt(18900.0) * x * x * y * y * y + std::sqrt(95681.25) * x * x * y * z * z + std::sqrt(1181.25) * y * y * y * y * y - std::sqrt(10631.25) * y * y * y * z * z) + e_3 * (-std::sqrt(10631.25) * x * x * y + std::sqrt(1181.25) * y * y * y);

        pc_45[k] = e_0 * (std::sqrt(30.76171875) * x * x * x * x * x * x * x * y * z + std::sqrt(30.76171875) * x * x * x * x * x * y * y * y * z - std::sqrt(218.75) * x * x * x * x * x * y * z * z * z - std::sqrt(30.76171875) * x * x * x * y * y * y * y * y * z + std::sqrt(8.75) * x * x * x * y * z * z * z * z * z - std::sqrt(30.76171875) * x * y * y * y * y * y * y * y * z + std::sqrt(218.75) * x * y * y * y * y * y * z * z * z - std::sqrt(8.75) * x * y * y * y * z * z * z * z * z) + e_1 * (std::sqrt(7875.0) * x * x * x * x * x * y * z - std::sqrt(14000.0) * x * x * x * y * z * z * z - std::sqrt(7875.0) * x * y * y * y * y * y * z + std::sqrt(14000.0) * x * y * y * y * z * z * z) + e_2 * (std::sqrt(70875.0) * x * x * x * y * z - std::sqrt(70875.0) * x * y * y * y * z);

        pc_46[k] = e_0 * (std::sqrt(138.427734375) * x * x * x * x * x * x * y * z * z + std::sqrt(384.521484375) * x * x * x * x * y * y * y * z * z - std::sqrt(984.375) * x * x * x * x * y * z * z * z * z + std::sqrt(15.380859375) * x * x * y * y * y * y * y * z * z - std::sqrt(437.5) * x * x * y * y * y * z * z * z * z + std::sqrt(39.375) * x * x * y * z * z * z * z * z * z - std::sqrt(15.380859375) * y * y * y * y * y * y * y * z * z + std::sqrt(109.375) * y * y * y * y * y * z * z * z * z - std::sqrt(4.375) * y * y * y * z * z * z * z * z * z) + e_1 * (std::sqrt(138.427734375) * x * x * x * x * x * x * y + std::sqrt(384.521484375) * x * x * x * x * y * y * y + std::sqrt(2214.84375) * x * x * x * x * y * z * z + std::sqrt(15.380859375) * x * x * y * y * y * y * y + std::sqrt(984.375) * x * x * y * y * y * z * z - std::sqrt(24609.375) * x * x * y * z * z * z * z - std::sqrt(15.380859375) * y * y * y * y * y * y * y - std::sqrt(246.09375) * y * y * y * y * y * z * z + std::sqrt(2734.375) * y * y * y * z * z * z * z) + e_2 * (std::sqrt(19933.59375) * x * x * x * x * y + std::sqrt(8859.375) * x * x * y * y * y - std::sqrt(79734.375) * x * x * y * z * z - std::sqrt(2214.84375) * y * y * y * y * y + std::sqrt(8859.375) * y * y * y * z * z) + e_3 * (std::sqrt(79734.375) * x * x * y - std::sqrt(8859.375) * y * y * y);

        pc_47[k] = e_0 * (-std::sqrt(4.39453125) * x * x * x * x * x * x * x * y * z - std::sqrt(39.55078125) * x * x * x * x * x * y * y * y * z + std::sqrt(330.078125) * x * x * x * x * x * y * z * z * z - std::sqrt(39.55078125) * x * x * x * y * y * y * y * y * z + std::sqrt(1320.3125) * x * x * x * y * y * y * z * z * z - std::sqrt(1201.25) * x * x * x * y * z * z * z * z * z - std::sqrt(4.39453125) * x * y * y * y * y * y * y * y * z + std::sqrt(330.078125) * x * y * y * y * y * y * z * z * z - std::sqrt(1201.25) * x * y * y * y * z * z * z * z * z + std::sqrt(45.0) * x * y * z * z * z * z * z * z * z) + e_1 * (-std::sqrt(70.3125) * x * x * x * x * x * y * z - std::sqrt(281.25) * x * x * x * y * y * y * z - std::sqrt(3125.0) * x * x * x * y * z * z * z - std::sqrt(70.3125) * x * y * y * y * y * y * z - std::sqrt(3125.0) * x * y * y * y * z * z * z - std::sqrt(4500.0) * x * y * z * z * z * z * z) + e_2 * (-std::sqrt(22781.25) * x * x * x * y * z - std::sqrt(22781.25) * x * y * y * y * z - std::sqrt(253125.0) * x * y * z * z * z) + e_3 * (-std::sqrt(648000.0) * x * y * z);
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

        pc_48[k] = e_0 * (-std::sqrt(19.775390625) * x * x * x * x * x * x * y * z * z - std::sqrt(177.978515625) * x * x * x * x * y * y * y * z * z + std::sqrt(316.40625) * x * x * x * x * y * z * z * z * z - std::sqrt(177.978515625) * x * x * y * y * y * y * y * z * z + std::sqrt(1265.625) * x * x * y * y * y * z * z * z * z - std::sqrt(330.625) * x * x * y * z * z * z * z * z * z - std::sqrt(19.775390625) * y * y * y * y * y * y * y * z * z + std::sqrt(316.40625) * y * y * y * y * y * z * z * z * z - std::sqrt(330.625) * y * y * y * z * z * z * z * z * z + std::sqrt(10.0) * y * z * z * z * z * z * z * z * z) + e_1 * (-std::sqrt(19.775390625) * x * x * x * x * x * x * y - std::sqrt(177.978515625) * x * x * x * x * y * y * y - std::sqrt(177.978515625) * x * x * y * y * y * y * y - std::sqrt(3515.625) * x * x * y * z * z * z * z - std::sqrt(19.775390625) * y * y * y * y * y * y * y - std::sqrt(3515.625) * y * y * y * z * z * z * z + std::sqrt(250.0) * y * z * z * z * z * z * z) + e_2 * (-std::sqrt(2847.65625) * x * x * x * x * y - std::sqrt(11390.625) * x * x * y * y * y - std::sqrt(31640.625) * x * x * y * z * z - std::sqrt(2847.65625) * y * y * y * y * y - std::sqrt(31640.625) * y * y * y * z * z) + e_3 * (-std::sqrt(74390.625) * x * x * y - std::sqrt(74390.625) * y * y * y - std::sqrt(56250.0) * y * z * z) + e_4 * (-std::sqrt(110250.0) * y);

        pc_49[k] = e_0 * (0.703125 * x * x * x * x * x * x * x * x * z + 2.8125 * x * x * x * x * x * x * y * y * z - 7.5 * x * x * x * x * x * x * z * z * z + 4.21875 * x * x * x * x * y * y * y * y * z - 22.5 * x * x * x * x * y * y * z * z * z + 17.25 * x * x * x * x * z * z * z * z * z + 2.8125 * x * x * y * y * y * y * y * y * z - 22.5 * x * x * y * y * y * y * z * z * z + 34.5 * x * x * y * y * z * z * z * z * z - 8.0 * x * x * z * z * z * z * z * z * z + 0.703125 * y * y * y * y * y * y * y * y * z - 7.5 * y * y * y * y * y * y * z * z * z + 17.25 * y * y * y * y * z * z * z * z * z - 8.0 * y * y * z * z * z * z * z * z * z + z * z * z * z * z * z * z * z * z) + e_1 * (37.5 * x * x * x * x * z * z * z + 75.0 * x * x * y * y * z * z * z - 30.0 * x * x * z * z * z * z * z + 37.5 * y * y * y * y * z * z * z - 30.0 * y * y * z * z * z * z * z + 20.0 * z * z * z * z * z * z * z) + e_2 * (56.25 * x * x * x * x * z + 112.5 * x * x * y * y * z + 56.25 * y * y * y * y * z + 180.0 * z * z * z * z * z) + e_3 * (150.0 * x * x * z + 150.0 * y * y * z + 600.0 * z * z * z) + e_4 * (525.0 * z);

        pc_50[k] = e_0 * (-std::sqrt(19.775390625) * x * x * x * x * x * x * x * z * z - std::sqrt(177.978515625) * x * x * x * x * x * y * y * z * z + std::sqrt(316.40625) * x * x * x * x * x * z * z * z * z - std::sqrt(177.978515625) * x * x * x * y * y * y * y * z * z + std::sqrt(1265.625) * x * x * x * y * y * z * z * z * z - std::sqrt(330.625) * x * x * x * z * z * z * z * z * z - std::sqrt(19.775390625) * x * y * y * y * y * y * y * z * z + std::sqrt(316.40625) * x * y * y * y * y * z * z * z * z - std::sqrt(330.625) * x * y * y * z * z * z * z * z * z + std::sqrt(10.0) * x * z * z * z * z * z * z * z * z) + e_1 * (-std::sqrt(19.775390625) * x * x * x * x * x * x * x - std::sqrt(177.978515625) * x * x * x * x * x * y * y - std::sqrt(177.978515625) * x * x * x * y * y * y * y - std::sqrt(3515.625) * x * x * x * z * z * z * z - std::sqrt(19.775390625) * x * y * y * y * y * y * y - std::sqrt(3515.625) * x * y * y * z * z * z * z + std::sqrt(250.0) * x * z * z * z * z * z * z) + e_2 * (-std::sqrt(2847.65625) * x * x * x * x * x - std::sqrt(11390.625) * x * x * x * y * y - std::sqrt(31640.625) * x * x * x * z * z - std::sqrt(2847.65625) * x * y * y * y * y - std::sqrt(31640.625) * x * y * y * z * z) + e_3 * (-std::sqrt(74390.625) * x * x * x - std::sqrt(74390.625) * x * y * y - std::sqrt(56250.0) * x * z * z) + e_4 * (-std::sqrt(110250.0) * x);

        pc_51[k] = e_0 * (-std::sqrt(1.0986328125) * x * x * x * x * x * x * x * x * z - std::sqrt(4.39453125) * x * x * x * x * x * x * y * y * z + std::sqrt(82.51953125) * x * x * x * x * x * x * z * z * z + std::sqrt(82.51953125) * x * x * x * x * y * y * z * z * z - std::sqrt(300.3125) * x * x * x * x * z * z * z * z * z + std::sqrt(4.39453125) * x * x * y * y * y * y * y * y * z - std::sqrt(82.51953125) * x * x * y * y * y * y * z * z * z + std::sqrt(11.25) * x * x * z * z * z * z * z * z * z + std::sqrt(1.0986328125) * y * y * y * y * y * y * y * y * z - std::sqrt(82.51953125) * y * y * y * y * y * y * z * z * z + std::sqrt(300.3125) * y * y * y * y * z * z * z * z * z - std::sqrt(11.25) * y * y * z * z * z * z * z * z * z) + e_1 * (-std::sqrt(17.578125) * x * x * x * x * x * x * z - std::sqrt(17.578125) * x * x * x * x * y * y * z - std::sqrt(781.25) * x * x * x * x * z * z * z + std::sqrt(17.578125) * x * x * y * y * y * y * z - std::sqrt(1125.0) * x * x * z * z * z * z * z + std::sqrt(17.578125) * y * y * y * y * y * y * z + std::sqrt(781.25) * y * y * y * y * z * z * z + std::sqrt(1125.0) * y * y * z * z * z * z * z) + e_2 * (-std::sqrt(5695.3125) * x * x * x * x * z - std::sqrt(63281.25) * x * x * z * z * z + std::sqrt(5695.3125) * y * y * y * y * z + std::sqrt(63281.25) * y * y * z * z * z) + e_3 * (-std::sqrt(162000.0) * x * x * z + std::sqrt(162000.0) * y * y * z);
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

        pc_52[k] = e_0 * (std::sqrt(15.380859375) * x * x * x * x * x * x * x * z * z - std::sqrt(15.380859375) * x * x * x * x * x * y * y * z * z - std::sqrt(109.375) * x * x * x * x * x * z * z * z * z - std::sqrt(384.521484375) * x * x * x * y * y * y * y * z * z + std::sqrt(437.5) * x * x * x * y * y * z * z * z * z + std::sqrt(4.375) * x * x * x * z * z * z * z * z * z - std::sqrt(138.427734375) * x * y * y * y * y * y * y * z * z + std::sqrt(984.375) * x * y * y * y * y * z * z * z * z - std::sqrt(39.375) * x * y * y * z * z * z * z * z * z) + e_1 * (std::sqrt(15.380859375) * x * x * x * x * x * x * x - std::sqrt(15.380859375) * x * x * x * x * x * y * y + std::sqrt(246.09375) * x * x * x * x * x * z * z - std::sqrt(384.521484375) * x * x * x * y * y * y * y - std::sqrt(984.375) * x * x * x * y * y * z * z - std::sqrt(2734.375) * x * x * x * z * z * z * z - std::sqrt(138.427734375) * x * y * y * y * y * y * y - std::sqrt(2214.84375) * x * y * y * y * y * z * z + std::sqrt(24609.375) * x * y * y * z * z * z * z) + e_2 * (std::sqrt(2214.84375) * x * x * x * x * x - std::sqrt(8859.375) * x * x * x * y * y - std::sqrt(8859.375) * x * x * x * z * z - std::sqrt(19933.59375) * x * y * y * y * y + std::sqrt(79734.375) * x * y * y * z * z) + e_3 * (std::sqrt(8859.375) * x * x * x - std::sqrt(79734.375) * x * y * y);

        pc_53[k] = e_0 * (std::sqrt(1.922607421875) * x * x * x * x * x * x * x * x * z - std::sqrt(30.76171875) * x * x * x * x * x * x * y * y * z - std::sqrt(13.671875) * x * x * x * x * x * x * z * z * z - std::sqrt(192.2607421875) * x * x * x * x * y * y * y * y * z + std::sqrt(341.796875) * x * x * x * x * y * y * z * z * z + std::sqrt(0.546875) * x * x * x * x * z * z * z * z * z - std::sqrt(30.76171875) * x * x * y * y * y * y * y * y * z + std::sqrt(341.796875) * x * x * y * y * y * y * z * z * z - std::sqrt(19.6875) * x * x * y * y * z * z * z * z * z + std::sqrt(1.922607421875) * y * y * y * y * y * y * y * y * z - std::sqrt(13.671875) * y * y * y * y * y * y * z * z * z + std::sqrt(0.546875) * y * y * y * y * z * z * z * z * z) + e_1 * (std::sqrt(492.1875) * x * x * x * x * x * x * z - std::sqrt(12304.6875) * x * x * x * x * y * y * z - std::sqrt(875.0) * x * x * x * x * z * z * z - std::sqrt(12304.6875) * x * x * y * y * y * y * z + std::sqrt(31500.0) * x * x * y * y * z * z * z + std::sqrt(492.1875) * y * y * y * y * y * y * z - std::sqrt(875.0) * y * y * y * y * z * z * z) + e_2 * (std::sqrt(4429.6875) * x * x * x * x * z - std::sqrt(159468.75) * x * x * y * y * z + std::sqrt(4429.6875) * y * y * y * y * z);

        pc_54[k] = e_0 * (std::sqrt(2.05078125) * x * x * x * x * x * x * x * x * y + std::sqrt(2.05078125) * x * x * x * x * x * x * y * y * y - std::sqrt(295.3125) * x * x * x * x * x * x * y * z * z - std::sqrt(2.05078125) * x * x * x * x * y * y * y * y * y + std::sqrt(131.25) * x * x * x * x * y * z * z * z * z - std::sqrt(2.05078125) * x * x * y * y * y * y * y * y * y + std::sqrt(295.3125) * x * x * y * y * y * y * y * z * z - std::sqrt(131.25) * x * x * y * y * y * z * z * z * z) + e_1 * (std::sqrt(740.33203125) * x * x * x * x * x * x * y + std::sqrt(51.26953125) * x * x * x * x * y * y * y - std::sqrt(35732.8125) * x * x * x * x * y * z * z - std::sqrt(461.42578125) * x * x * y * y * y * y * y + std::sqrt(10631.25) * x * x * y * y * y * z * z + std::sqrt(1181.25) * x * x * y * z * z * z * z - std::sqrt(2.05078125) * y * y * y * y * y * y * y + std::sqrt(295.3125) * y * y * y * y * y * z * z - std::sqrt(131.25) * y * y * y * z * z * z * z) + e_2 * (std::sqrt(14470.3125) * x * x * x * x * y - std::sqrt(1181.25) * x * x * y * y * y - std::sqrt(95681.25) * x * x * y * z * z - std::sqrt(295.3125) * y * y * y * y * y + std::sqrt(10631.25) * y * y * y * z * z) + e_3 * (std::sqrt(10631.25) * x * x * y - std::sqrt(1181.25) * y * y * y);

        pc_55[k] = e_0 * (std::sqrt(9.228515625) * x * x * x * x * x * x * x * y * z + std::sqrt(25.634765625) * x * x * x * x * x * y * y * y * z - std::sqrt(1328.90625) * x * x * x * x * x * y * z * z * z + std::sqrt(1.025390625) * x * x * x * y * y * y * y * y * z - std::sqrt(590.625) * x * x * x * y * y * y * z * z * z + std::sqrt(590.625) * x * x * x * y * z * z * z * z * z - std::sqrt(1.025390625) * x * y * y * y * y * y * y * y * z + std::sqrt(147.65625) * x * y * y * y * y * y * z * z * z - std::sqrt(65.625) * x * y * y * y * z * z * z * z * z) + e_1 * (-std::sqrt(922.8515625) * x * x * x * x * x * y * z - std::sqrt(147.65625) * x * x * x * y * y * y * z - std::sqrt(37800.0) * x * x * x * y * z * z * z + std::sqrt(332.2265625) * x * y * y * y * y * y * z - std::sqrt(1050.0) * x * y * y * y * z * z * z + std::sqrt(2362.5) * x * y * z * z * z * z * z) + e_2 * (-std::sqrt(213215.625) * x * x * x * y * z + std::sqrt(590.625) * x * y * y * y * z - std::sqrt(9450.0) * x * y * z * z * z) + e_3 * (-std::sqrt(285862.5) * x * y * z);
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

        pc_56[k] = e_0 * (-std::sqrt(0.29296875) * x * x * x * x * x * x * x * x * y - std::sqrt(2.63671875) * x * x * x * x * x * x * y * y * y + std::sqrt(94.921875) * x * x * x * x * x * x * y * z * z - std::sqrt(2.63671875) * x * x * x * x * y * y * y * y * y + std::sqrt(379.6875) * x * x * x * x * y * y * y * z * z - std::sqrt(1875.0) * x * x * x * x * y * z * z * z * z - std::sqrt(0.29296875) * x * x * y * y * y * y * y * y * y + std::sqrt(94.921875) * x * x * y * y * y * y * y * z * z - std::sqrt(1875.0) * x * x * y * y * y * z * z * z * z + std::sqrt(675.0) * x * x * y * z * z * z * z * z * z) + e_1 * (-std::sqrt(105.76171875) * x * x * x * x * x * x * y - std::sqrt(445.60546875) * x * x * x * x * y * y * y - std::sqrt(3048.046875) * x * x * x * x * y * z * z - std::sqrt(129.19921875) * x * x * y * y * y * y * y - std::sqrt(2067.1875) * x * x * y * y * y * z * z - std::sqrt(0.29296875) * y * y * y * y * y * y * y + std::sqrt(94.921875) * y * y * y * y * y * z * z - std::sqrt(1875.0) * y * y * y * z * z * z * z + std::sqrt(675.0) * y * z * z * z * z * z * z) + e_2 * (-std::sqrt(18604.6875) * x * x * x * x * y - std::sqrt(20418.75) * x * x * y * y * y - std::sqrt(54675.0) * x * x * y * z * z - std::sqrt(42.1875) * y * y * y * y * y - std::sqrt(10800.0) * y * y * y * z * z + std::sqrt(16875.0) * y * z * z * z * z) + e_3 * (-std::sqrt(243675.0) * x * x * y - std::sqrt(10800.0) * y * y * y + std::sqrt(6075.0) * y * z * z) + e_4 * (-std::sqrt(33075.0) * y);

        pc_57[k] = e_0 * (-std::sqrt(1.318359375) * x * x * x * x * x * x * x * y * z - std::sqrt(11.865234375) * x * x * x * x * x * y * y * y * z + std::sqrt(234.375) * x * x * x * x * x * y * z * z * z - std::sqrt(11.865234375) * x * x * x * y * y * y * y * y * z + std::sqrt(937.5) * x * x * x * y * y * y * z * z * z - std::sqrt(759.375) * x * x * x * y * z * z * z * z * z - std::sqrt(1.318359375) * x * y * y * y * y * y * y * y * z + std::sqrt(234.375) * x * y * y * y * y * y * z * z * z - std::sqrt(759.375) * x * y * y * y * z * z * z * z * z + std::sqrt(150.0) * x * y * z * z * z * z * z * z * z) + e_1 * (std::sqrt(131.8359375) * x * x * x * x * x * y * z + std::sqrt(527.34375) * x * x * x * y * y * y * z - std::sqrt(937.5) * x * x * x * y * z * z * z + std::sqrt(131.8359375) * x * y * y * y * y * y * z - std::sqrt(937.5) * x * y * y * y * z * z * z + std::sqrt(8437.5) * x * y * z * z * z * z * z) + e_2 * (std::sqrt(2109.375) * x * x * x * y * z + std::sqrt(2109.375) * x * y * y * y * z + std::sqrt(135000.0) * x * y * z * z * z) + e_3 * (std::sqrt(210937.5) * x * y * z);

        pc_58[k] = e_0 * (std::sqrt(0.032958984375) * x * x * x * x * x * x * x * x * x + std::sqrt(0.52734375) * x * x * x * x * x * x * x * y * y - std::sqrt(13.18359375) * x * x * x * x * x * x * x * z * z + std::sqrt(1.1865234375) * x * x * x * x * x * y * y * y * y - std::sqrt(118.65234375) * x * x * x * x * x * y * y * z * z + std::sqrt(375.0) * x * x * x * x * x * z * z * z * z + std::sqrt(0.52734375) * x * x * x * y * y * y * y * y * y - std::sqrt(118.65234375) * x * x * x * y * y * y * y * z * z + std::sqrt(1500.0) * x * x * x * y * y * z * z * z * z - std::sqrt(303.75) * x * x * x * z * z * z * z * z * z + std::sqrt(0.032958984375) * x * y * y * y * y * y * y * y * y - std::sqrt(13.18359375) * x * y * y * y * y * y * y * z * z + std::sqrt(375.0) * x * y * y * y * y * z * z * z * z - std::sqrt(303.75) * x * y * y * z * z * z * z * z * z + std::sqrt(15.0) * x * z * z * z * z * z * z * z * z) + e_1 * (std::sqrt(13.18359375) * x * x * x * x * x * x * x + std::sqrt(118.65234375) * x * x * x * x * x * y * y + std::sqrt(843.75) * x * x * x * x * x * z * z + std::sqrt(118.65234375) * x * x * x * y * y * y * y + std::sqrt(3375.0) * x * x * x * y * y * z * z - std::sqrt(843.75) * x * x * x * z * z * z * z + std::sqrt(13.18359375) * x * y * y * y * y * y * y + std::sqrt(843.75) * x * y * y * y * y * z * z - std::sqrt(843.75) * x * y * y * z * z * z * z + std::sqrt(1500.0) * x * z * z * z * z * z * z) + e_2 * (std::sqrt(3375.0) * x * x * x * x * x + std::sqrt(13500.0) * x * x * x * y * y + std::sqrt(7593.75) * x * x * x * z * z + std::sqrt(3375.0) * x * y * y * y * y + std::sqrt(7593.75) * x * y * y * z * z + std::sqrt(54000.0) * x * z * z * z * z) + e_3 * (std::sqrt(68343.75) * x * x * x + std::sqrt(68343.75) * x * y * y + std::sqrt(337500.0) * x * z * z) + e_4 * (std::sqrt(165375.0) * x);

        pc_59[k] = e_0 * (-std::sqrt(1.318359375) * x * x * x * x * x * x * x * x * z - std::sqrt(11.865234375) * x * x * x * x * x * x * y * y * z + std::sqrt(234.375) * x * x * x * x * x * x * z * z * z - std::sqrt(11.865234375) * x * x * x * x * y * y * y * y * z + std::sqrt(937.5) * x * x * x * x * y * y * z * z * z - std::sqrt(759.375) * x * x * x * x * z * z * z * z * z - std::sqrt(1.318359375) * x * x * y * y * y * y * y * y * z + std::sqrt(234.375) * x * x * y * y * y * y * z * z * z - std::sqrt(759.375) * x * x * y * y * z * z * z * z * z + std::sqrt(150.0) * x * x * z * z * z * z * z * z * z) + e_1 * (std::sqrt(106.787109375) * x * x * x * x * x * x * z + std::sqrt(381.005859375) * x * x * x * x * y * y * z - std::sqrt(234.375) * x * x * x * x * z * z * z + std::sqrt(64.599609375) * x * x * y * y * y * y * z + std::sqrt(4134.375) * x * x * z * z * z * z * z - std::sqrt(1.318359375) * y * y * y * y * y * y * z + std::sqrt(234.375) * y * y * y * y * z * z * z - std::sqrt(759.375) * y * y * z * z * z * z * z + std::sqrt(150.0) * z * z * z * z * z * z * z) + e_2 * (std::sqrt(4134.375) * x * x * x * x * z + std::sqrt(6834.375) * x * x * y * y * z + std::sqrt(75937.5) * x * x * z * z * z + std::sqrt(337.5) * y * y * y * y * z - std::sqrt(8437.5) * y * y * z * z * z + std::sqrt(21600.0) * z * z * z * z * z) + e_3 * (std::sqrt(186384.375) * x * x * z - std::sqrt(759.375) * y * y * z + std::sqrt(303750.0) * z * z * z) + e_4 * (std::sqrt(264600.0) * z);
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

        pc_60[k] = e_0 * (-std::sqrt(0.0732421875) * x * x * x * x * x * x * x * x * x - std::sqrt(0.29296875) * x * x * x * x * x * x * x * y * y + std::sqrt(23.73046875) * x * x * x * x * x * x * x * z * z + std::sqrt(23.73046875) * x * x * x * x * x * y * y * z * z - std::sqrt(468.75) * x * x * x * x * x * z * z * z * z + std::sqrt(0.29296875) * x * x * x * y * y * y * y * y * y - std::sqrt(23.73046875) * x * x * x * y * y * y * y * z * z + std::sqrt(168.75) * x * x * x * z * z * z * z * z * z + std::sqrt(0.0732421875) * x * y * y * y * y * y * y * y * y - std::sqrt(23.73046875) * x * y * y * y * y * y * y * z * z + std::sqrt(468.75) * x * y * y * y * y * z * z * z * z - std::sqrt(168.75) * x * y * y * z * z * z * z * z * z) + e_1 * (-std::sqrt(29.296875) * x * x * x * x * x * x * x - std::sqrt(42.1875) * x * x * x * x * x * y * y - std::sqrt(516.796875) * x * x * x * x * x * z * z + std::sqrt(10.546875) * x * x * x * y * y * y * y + std::sqrt(379.6875) * x * x * x * y * y * z * z - std::sqrt(468.75) * x * x * x * z * z * z * z + std::sqrt(18.75) * x * y * y * y * y * y * y + std::sqrt(1782.421875) * x * y * y * y * y * z * z - std::sqrt(4218.75) * x * y * y * z * z * z * z + std::sqrt(675.0) * x * z * z * z * z * z * z) + e_2 * (-std::sqrt(5104.6875) * x * x * x * x * x - std::sqrt(168.75) * x * x * x * y * y - std::sqrt(28518.75) * x * x * x * z * z + std::sqrt(3417.1875) * x * y * y * y * y - std::sqrt(1518.75) * x * y * y * z * z + std::sqrt(16875.0) * x * z * z * z * z) + e_3 * (-std::sqrt(89268.75) * x * x * x + std::sqrt(8268.75) * x * y * y + std::sqrt(6075.0) * x * z * z) + e_4 * (-std::sqrt(33075.0) * x);

        pc_61[k] = e_0 * (std::sqrt(1.025390625) * x * x * x * x * x * x * x * x * z - std::sqrt(1.025390625) * x * x * x * x * x * x * y * y * z - std::sqrt(147.65625) * x * x * x * x * x * x * z * z * z - std::sqrt(25.634765625) * x * x * x * x * y * y * y * y * z + std::sqrt(590.625) * x * x * x * x * y * y * z * z * z + std::sqrt(65.625) * x * x * x * x * z * z * z * z * z - std::sqrt(9.228515625) * x * x * y * y * y * y * y * y * z + std::sqrt(1328.90625) * x * x * y * y * y * y * z * z * z - std::sqrt(590.625) * x * x * y * y * z * z * z * z * z) + e_1 * (-std::sqrt(83.056640625) * x * x * x * x * x * x * z + std::sqrt(747.509765625) * x * x * x * x * y * y * z - std::sqrt(5922.65625) * x * x * x * x * z * z * z + std::sqrt(1116.650390625) * x * x * y * y * y * y * z + std::sqrt(14765.625) * x * x * y * y * z * z * z + std::sqrt(590.625) * x * x * z * z * z * z * z - std::sqrt(9.228515625) * y * y * y * y * y * y * z + std::sqrt(1328.90625) * y * y * y * y * z * z * z - std::sqrt(590.625) * y * y * z * z * z * z * z) + e_2 * (-std::sqrt(28940.625) * x * x * x * x * z + std::sqrt(132890.625) * x * x * y * y * z - std::sqrt(2362.5) * x * x * z * z * z + std::sqrt(2362.5) * y * y * y * y * z + std::sqrt(2362.5) * y * y * z * z * z) + e_3 * (-std::sqrt(71465.625) * x * x * z + std::sqrt(71465.625) * y * y * z);

        pc_62[k] = e_0 * (std::sqrt(0.128173828125) * x * x * x * x * x * x * x * x * x - std::sqrt(2.05078125) * x * x * x * x * x * x * x * y * y - std::sqrt(18.45703125) * x * x * x * x * x * x * x * z * z - std::sqrt(12.8173828125) * x * x * x * x * x * y * y * y * y + std::sqrt(461.42578125) * x * x * x * x * x * y * y * z * z + std::sqrt(8.203125) * x * x * x * x * x * z * z * z * z - std::sqrt(2.05078125) * x * x * x * y * y * y * y * y * y + std::sqrt(461.42578125) * x * x * x * y * y * y * y * z * z - std::sqrt(295.3125) * x * x * x * y * y * z * z * z * z + std::sqrt(0.128173828125) * x * y * y * y * y * y * y * y * y - std::sqrt(18.45703125) * x * y * y * y * y * y * y * z * z + std::sqrt(8.203125) * x * y * y * y * y * z * z * z * z) + e_1 * (std::sqrt(51.26953125) * x * x * x * x * x * x * x - std::sqrt(904.39453125) * x * x * x * x * x * y * y - std::sqrt(2657.8125) * x * x * x * x * x * z * z - std::sqrt(1281.73828125) * x * x * x * y * y * y * y + std::sqrt(57881.25) * x * x * x * y * y * z * z + std::sqrt(131.25) * x * x * x * z * z * z * z + std::sqrt(2.05078125) * x * y * y * y * y * y * y + std::sqrt(295.3125) * x * y * y * y * y * z * z - std::sqrt(1181.25) * x * y * y * z * z * z * z) + e_2 * (std::sqrt(1181.25) * x * x * x * x * x - std::sqrt(18900.0) * x * x * x * y * y - std::sqrt(10631.25) * x * x * x * z * z - std::sqrt(1181.25) * x * y * y * y * y + std::sqrt(95681.25) * x * y * y * z * z) + e_3 * (std::sqrt(1181.25) * x * x * x - std::sqrt(10631.25) * x * y * y);

        pc_63[k] = e_0 * (-std::sqrt(57.421875) * x * x * x * x * x * x * x * y * z + std::sqrt(57.421875) * x * x * x * x * x * y * y * y * z + std::sqrt(229.6875) * x * x * x * x * x * y * z * z * z + std::sqrt(57.421875) * x * x * x * y * y * y * y * y * z - std::sqrt(918.75) * x * x * x * y * y * y * z * z * z - std::sqrt(57.421875) * x * y * y * y * y * y * y * y * z + std::sqrt(229.6875) * x * y * y * y * y * y * z * z * z) + e_1 * (-std::sqrt(8268.75) * x * x * x * x * x * y * z + std::sqrt(3675.0) * x * x * x * y * y * y * z + std::sqrt(3675.0) * x * x * x * y * z * z * z - std::sqrt(8268.75) * x * y * y * y * y * y * z + std::sqrt(3675.0) * x * y * y * y * z * z * z) + e_2 * (-std::sqrt(74418.75) * x * x * x * y * z - std::sqrt(74418.75) * x * y * y * y * z + std::sqrt(33075.0) * x * y * z * z * z) + e_3 * (-std::sqrt(132300.0) * x * y * z);
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

        pc_64[k] = e_0 * (-std::sqrt(258.3984375) * x * x * x * x * x * x * y * z * z + std::sqrt(28.7109375) * x * x * x * x * y * y * y * z * z + std::sqrt(1033.59375) * x * x * x * x * y * z * z * z * z + std::sqrt(258.3984375) * x * x * y * y * y * y * y * z * z - std::sqrt(1837.5) * x * x * y * y * y * z * z * z * z - std::sqrt(28.7109375) * y * y * y * y * y * y * y * z * z + std::sqrt(114.84375) * y * y * y * y * y * z * z * z * z) + e_1 * (-std::sqrt(258.3984375) * x * x * x * x * x * x * y + std::sqrt(28.7109375) * x * x * x * x * y * y * y - std::sqrt(1033.59375) * x * x * x * x * y * z * z + std::sqrt(258.3984375) * x * x * y * y * y * y * y - std::sqrt(4134.375) * x * x * y * y * y * z * z + std::sqrt(4134.375) * x * x * y * z * z * z * z - std::sqrt(28.7109375) * y * y * y * y * y * y * y - std::sqrt(1033.59375) * y * y * y * y * y * z * z + std::sqrt(4134.375) * y * y * y * z * z * z * z) + e_2 * (-std::sqrt(16537.5) * x * x * x * x * y + std::sqrt(4134.375) * x * x * y * y * y - std::sqrt(4134.375) * y * y * y * y * y + std::sqrt(16537.5) * y * z * z * z * z) + e_3 * (-std::sqrt(37209.375) * x * x * y - std::sqrt(37209.375) * y * y * y + std::sqrt(66150.0) * y * z * z) + e_4 * (-std::sqrt(16537.5) * y);

        pc_65[k] = e_0 * (std::sqrt(8.203125) * x * x * x * x * x * x * x * y * z + std::sqrt(8.203125) * x * x * x * x * x * y * y * y * z - std::sqrt(525.0) * x * x * x * x * x * y * z * z * z - std::sqrt(8.203125) * x * x * x * y * y * y * y * y * z + std::sqrt(1181.25) * x * x * x * y * z * z * z * z * z - std::sqrt(8.203125) * x * y * y * y * y * y * y * y * z + std::sqrt(525.0) * x * y * y * y * y * y * z * z * z - std::sqrt(1181.25) * x * y * y * y * z * z * z * z * z) + e_1 * (std::sqrt(13125.0) * x * x * x * y * z * z * z - std::sqrt(13125.0) * x * y * y * y * z * z * z) + e_2 * (std::sqrt(29531.25) * x * x * x * y * z - std::sqrt(29531.25) * x * y * y * y * z);

        pc_66[k] = e_0 * (std::sqrt(36.9140625) * x * x * x * x * x * x * y * z * z + std::sqrt(36.9140625) * x * x * x * x * y * y * y * z * z - std::sqrt(410.15625) * x * x * x * x * y * z * z * z * z - std::sqrt(36.9140625) * x * x * y * y * y * y * y * z * z + std::sqrt(262.5) * x * x * y * z * z * z * z * z * z - std::sqrt(36.9140625) * y * y * y * y * y * y * y * z * z + std::sqrt(410.15625) * y * y * y * y * y * z * z * z * z - std::sqrt(262.5) * y * y * y * z * z * z * z * z * z) + e_1 * (std::sqrt(36.9140625) * x * x * x * x * x * x * y + std::sqrt(36.9140625) * x * x * x * x * y * y * y - std::sqrt(147.65625) * x * x * x * x * y * z * z - std::sqrt(36.9140625) * x * x * y * y * y * y * y - std::sqrt(590.625) * x * x * y * y * y * z * z + std::sqrt(14765.625) * x * x * y * z * z * z * z - std::sqrt(36.9140625) * y * y * y * y * y * y * y - std::sqrt(147.65625) * y * y * y * y * y * z * z - std::sqrt(1640.625) * y * y * y * z * z * z * z - std::sqrt(1050.0) * y * z * z * z * z * z * z) + e_2 * (std::sqrt(2362.5) * x * x * x * x * y - std::sqrt(590.625) * x * x * y * y * y + std::sqrt(85050.0) * x * x * y * z * z - std::sqrt(5315.625) * y * y * y * y * y - std::sqrt(37800.0) * y * y * y * z * z - std::sqrt(59062.5) * y * z * z * z * z) + e_3 * (std::sqrt(28940.625) * x * x * y - std::sqrt(99815.625) * y * y * y - std::sqrt(340200.0) * y * z * z) + e_4 * (-std::sqrt(115762.5) * y);

        pc_67[k] = e_0 * (-std::sqrt(0.9228515625) * x * x * x * x * x * x * x * x * z - std::sqrt(3.69140625) * x * x * x * x * x * x * y * y * z + std::sqrt(92.28515625) * x * x * x * x * x * x * z * z * z + std::sqrt(92.28515625) * x * x * x * x * y * y * z * z * z - std::sqrt(321.5625) * x * x * x * x * z * z * z * z * z + std::sqrt(3.69140625) * x * x * y * y * y * y * y * y * z - std::sqrt(92.28515625) * x * x * y * y * y * y * z * z * z + std::sqrt(26.25) * x * x * z * z * z * z * z * z * z + std::sqrt(0.9228515625) * y * y * y * y * y * y * y * y * z - std::sqrt(92.28515625) * y * y * y * y * y * y * z * z * z + std::sqrt(321.5625) * y * y * y * y * z * z * z * z * z - std::sqrt(26.25) * y * y * z * z * z * z * z * z * z) + e_1 * (-std::sqrt(656.25) * x * x * x * x * z * z * z + std::sqrt(656.25) * y * y * y * y * z * z * z) + e_2 * (-std::sqrt(1476.5625) * x * x * x * x * z - std::sqrt(5906.25) * x * x * z * z * z + std::sqrt(1476.5625) * y * y * y * y * z + std::sqrt(5906.25) * y * y * z * z * z) + e_3 * (-std::sqrt(23625.0) * x * x * z + std::sqrt(23625.0) * y * y * z);
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

        pc_68[k] = e_0 * (std::sqrt(36.9140625) * x * x * x * x * x * x * x * z * z + std::sqrt(36.9140625) * x * x * x * x * x * y * y * z * z - std::sqrt(410.15625) * x * x * x * x * x * z * z * z * z - std::sqrt(36.9140625) * x * x * x * y * y * y * y * z * z + std::sqrt(262.5) * x * x * x * z * z * z * z * z * z - std::sqrt(36.9140625) * x * y * y * y * y * y * y * z * z + std::sqrt(410.15625) * x * y * y * y * y * z * z * z * z - std::sqrt(262.5) * x * y * y * z * z * z * z * z * z) + e_1 * (std::sqrt(36.9140625) * x * x * x * x * x * x * x + std::sqrt(36.9140625) * x * x * x * x * x * y * y + std::sqrt(147.65625) * x * x * x * x * x * z * z - std::sqrt(36.9140625) * x * x * x * y * y * y * y + std::sqrt(590.625) * x * x * x * y * y * z * z + std::sqrt(1640.625) * x * x * x * z * z * z * z - std::sqrt(36.9140625) * x * y * y * y * y * y * y + std::sqrt(147.65625) * x * y * y * y * y * z * z - std::sqrt(14765.625) * x * y * y * z * z * z * z + std::sqrt(1050.0) * x * z * z * z * z * z * z) + e_2 * (std::sqrt(5315.625) * x * x * x * x * x + std::sqrt(590.625) * x * x * x * y * y + std::sqrt(37800.0) * x * x * x * z * z - std::sqrt(2362.5) * x * y * y * y * y - std::sqrt(85050.0) * x * y * y * z * z + std::sqrt(59062.5) * x * z * z * z * z) + e_3 * (std::sqrt(99815.625) * x * x * x - std::sqrt(28940.625) * x * y * y + std::sqrt(340200.0) * x * z * z) + e_4 * (std::sqrt(115762.5) * x);

        pc_69[k] = e_0 * (std::sqrt(2.05078125) * x * x * x * x * x * x * x * x * z - std::sqrt(131.25) * x * x * x * x * x * x * z * z * z - std::sqrt(8.203125) * x * x * x * x * y * y * y * y * z + std::sqrt(131.25) * x * x * x * x * y * y * z * z * z + std::sqrt(295.3125) * x * x * x * x * z * z * z * z * z + std::sqrt(131.25) * x * x * y * y * y * y * z * z * z - std::sqrt(1181.25) * x * x * y * y * z * z * z * z * z + std::sqrt(2.05078125) * y * y * y * y * y * y * y * y * z - std::sqrt(131.25) * y * y * y * y * y * y * z * z * z + std::sqrt(295.3125) * y * y * y * y * z * z * z * z * z) + e_1 * (std::sqrt(32.8125) * x * x * x * x * x * x * z + std::sqrt(295.3125) * x * x * x * x * y * y * z + std::sqrt(131.25) * x * x * x * x * z * z * z + std::sqrt(295.3125) * x * x * y * y * y * y * z - std::sqrt(42525.0) * x * x * y * y * z * z * z + std::sqrt(4725.0) * x * x * z * z * z * z * z + std::sqrt(32.8125) * y * y * y * y * y * y * z + std::sqrt(131.25) * y * y * y * y * z * z * z + std::sqrt(4725.0) * y * y * z * z * z * z * z) + e_2 * (std::sqrt(4725.0) * x * x * x * x * z - std::sqrt(42525.0) * x * x * y * y * z + std::sqrt(75600.0) * x * x * z * z * z + std::sqrt(4725.0) * y * y * y * y * z + std::sqrt(75600.0) * y * y * z * z * z + std::sqrt(4725.0) * z * z * z * z * z) + e_3 * (std::sqrt(118125.0) * x * x * z + std::sqrt(118125.0) * y * y * z + std::sqrt(170100.0) * z * z * z) + e_4 * (std::sqrt(231525.0) * z);

        pc_70[k] = e_0 * (-std::sqrt(28.7109375) * x * x * x * x * x * x * x * z * z + std::sqrt(258.3984375) * x * x * x * x * x * y * y * z * z + std::sqrt(114.84375) * x * x * x * x * x * z * z * z * z + std::sqrt(28.7109375) * x * x * x * y * y * y * y * z * z - std::sqrt(1837.5) * x * x * x * y * y * z * z * z * z - std::sqrt(258.3984375) * x * y * y * y * y * y * y * z * z + std::sqrt(1033.59375) * x * y * y * y * y * z * z * z * z) + e_1 * (-std::sqrt(28.7109375) * x * x * x * x * x * x * x + std::sqrt(258.3984375) * x * x * x * x * x * y * y - std::sqrt(1033.59375) * x * x * x * x * x * z * z + std::sqrt(28.7109375) * x * x * x * y * y * y * y - std::sqrt(4134.375) * x * x * x * y * y * z * z + std::sqrt(4134.375) * x * x * x * z * z * z * z - std::sqrt(258.3984375) * x * y * y * y * y * y * y - std::sqrt(1033.59375) * x * y * y * y * y * z * z + std::sqrt(4134.375) * x * y * y * z * z * z * z) + e_2 * (-std::sqrt(4134.375) * x * x * x * x * x + std::sqrt(4134.375) * x * x * x * y * y - std::sqrt(16537.5) * x * y * y * y * y + std::sqrt(16537.5) * x * z * z * z * z) + e_3 * (-std::sqrt(37209.375) * x * x * x - std::sqrt(37209.375) * x * y * y + std::sqrt(66150.0) * x * z * z) + e_4 * (-std::sqrt(16537.5) * x);

        pc_71[k] = e_0 * (-std::sqrt(3.5888671875) * x * x * x * x * x * x * x * x * z + std::sqrt(129.19921875) * x * x * x * x * x * x * y * y * z + std::sqrt(14.35546875) * x * x * x * x * x * x * z * z * z - std::sqrt(703.41796875) * x * x * x * x * y * y * z * z * z - std::sqrt(129.19921875) * x * x * y * y * y * y * y * y * z + std::sqrt(703.41796875) * x * x * y * y * y * y * z * z * z + std::sqrt(3.5888671875) * y * y * y * y * y * y * y * y * z - std::sqrt(14.35546875) * y * y * y * y * y * y * z * z * z) + e_1 * (-std::sqrt(918.75) * x * x * x * x * x * x * z + std::sqrt(8268.75) * x * x * x * x * y * y * z + std::sqrt(918.75) * x * x * x * x * z * z * z - std::sqrt(8268.75) * x * x * y * y * y * y * z + std::sqrt(918.75) * y * y * y * y * y * y * z - std::sqrt(918.75) * y * y * y * y * z * z * z) + e_2 * (-std::sqrt(18604.6875) * x * x * x * x * z + std::sqrt(8268.75) * x * x * z * z * z + std::sqrt(18604.6875) * y * y * y * y * z - std::sqrt(8268.75) * y * y * z * z * z) + e_3 * (-std::sqrt(33075.0) * x * x * z + std::sqrt(33075.0) * y * y * z);
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

        pc_72[k] = e_0 * (-std::sqrt(2.392578125) * x * x * x * x * x * x * x * x * y + std::sqrt(21.533203125) * x * x * x * x * x * x * y * y * y + std::sqrt(153.125) * x * x * x * x * x * x * y * z * z + std::sqrt(2.392578125) * x * x * x * x * y * y * y * y * y - std::sqrt(2450.0) * x * x * x * x * y * y * y * z * z - std::sqrt(21.533203125) * x * x * y * y * y * y * y * y * y + std::sqrt(1378.125) * x * x * y * y * y * y * y * z * z) + e_1 * (-std::sqrt(289.501953125) * x * x * x * x * x * x * y + std::sqrt(1265.673828125) * x * x * x * x * y * y * y + std::sqrt(1378.125) * x * x * x * x * y * z * z - std::sqrt(2605.517578125) * x * x * y * y * y * y * y + std::sqrt(5512.5) * x * x * y * y * y * z * z - std::sqrt(21.533203125) * y * y * y * y * y * y * y + std::sqrt(1378.125) * y * y * y * y * y * z * z) + e_2 * (-std::sqrt(3100.78125) * x * x * x * x * y - std::sqrt(12403.125) * x * x * y * y * y + std::sqrt(49612.5) * x * x * y * z * z - std::sqrt(3100.78125) * y * y * y * y * y + std::sqrt(49612.5) * y * y * y * z * z) + e_3 * (-std::sqrt(22050.0) * x * x * y - std::sqrt(22050.0) * y * y * y + std::sqrt(88200.0) * y * z * z) + e_4 * (-std::sqrt(5512.5) * y);

        pc_73[k] = e_0 * (-3.28125 * x * x * x * x * x * x * x * y * z + 7.65625 * x * x * x * x * x * y * y * y * z + 26.25 * x * x * x * x * x * y * z * z * z + 7.65625 * x * x * x * y * y * y * y * y * z - 87.5 * x * x * x * y * y * y * z * z * z - 3.28125 * x * y * y * y * y * y * y * y * z + 26.25 * x * y * y * y * y * y * z * z * z) + e_1 * (32.8125 * x * x * x * x * x * y * z - 109.375 * x * x * x * y * y * y * z + 32.8125 * x * y * y * y * y * y * z);

        pc_74[k] = e_0 * (std::sqrt(0.341796875) * x * x * x * x * x * x * x * x * y - std::sqrt(0.341796875) * x * x * x * x * x * x * y * y * y - std::sqrt(66.9921875) * x * x * x * x * x * x * y * z * z - std::sqrt(8.544921875) * x * x * x * x * y * y * y * y * y + std::sqrt(267.96875) * x * x * x * x * y * y * y * z * z + std::sqrt(787.5) * x * x * x * x * y * z * z * z * z - std::sqrt(3.076171875) * x * x * y * y * y * y * y * y * y + std::sqrt(602.9296875) * x * x * y * y * y * y * y * z * z - std::sqrt(7087.5) * x * x * y * y * y * z * z * z * z) + e_1 * (std::sqrt(41.357421875) * x * x * x * x * x * x * y - std::sqrt(467.919921875) * x * x * x * x * y * y * y + std::sqrt(8970.1171875) * x * x * x * x * y * z * z - std::sqrt(889.013671875) * x * x * y * y * y * y * y - std::sqrt(26036.71875) * x * x * y * y * y * z * z - std::sqrt(7087.5) * x * x * y * z * z * z * z - std::sqrt(3.076171875) * y * y * y * y * y * y * y + std::sqrt(602.9296875) * y * y * y * y * y * z * z - std::sqrt(7087.5) * y * y * y * z * z * z * z) + e_2 * (std::sqrt(3986.71875) * x * x * x * x * y - std::sqrt(86821.875) * x * x * y * y * y - std::sqrt(44296.875) * x * x * y * z * z - std::sqrt(442.96875) * y * y * y * y * y - std::sqrt(44296.875) * y * y * y * z * z - std::sqrt(28350.0) * y * z * z * z * z) + e_3 * (-std::sqrt(56896.875) * x * x * y - std::sqrt(56896.875) * y * y * y - std::sqrt(381150.0) * y * z * z) + e_4 * (-std::sqrt(154350.0) * y);

        pc_75[k] = e_0 * (std::sqrt(1.5380859375) * x * x * x * x * x * x * x * y * z - std::sqrt(1.5380859375) * x * x * x * x * x * y * y * y * z - std::sqrt(133.984375) * x * x * x * x * x * y * z * z * z - std::sqrt(38.4521484375) * x * x * x * y * y * y * y * y * z + std::sqrt(535.9375) * x * x * x * y * y * y * z * z * z + std::sqrt(175.0) * x * x * x * y * z * z * z * z * z - std::sqrt(13.8427734375) * x * y * y * y * y * y * y * y * z + std::sqrt(1205.859375) * x * y * y * y * y * y * z * z * z - std::sqrt(1575.0) * x * y * y * y * z * z * z * z * z) + e_1 * (-std::sqrt(153.80859375) * x * x * x * x * x * y * z - std::sqrt(24.609375) * x * x * x * y * y * y * z + std::sqrt(7393.75) * x * x * x * y * z * z * z + std::sqrt(55.37109375) * x * y * y * y * y * y * z + std::sqrt(393.75) * x * y * y * y * z * z * z - std::sqrt(6300.0) * x * y * z * z * z * z * z) + e_2 * (std::sqrt(3543.75) * x * x * x * y * z + std::sqrt(3543.75) * x * y * y * y * z - std::sqrt(56700.0) * x * y * z * z * z) + e_3 * (-std::sqrt(14175.0) * x * y * z);
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

        pc_76[k] = e_0 * (-std::sqrt(0.0384521484375) * x * x * x * x * x * x * x * x * x + std::sqrt(9.84375) * x * x * x * x * x * x * x * z * z + std::sqrt(1.38427734375) * x * x * x * x * x * y * y * y * y - std::sqrt(9.84375) * x * x * x * x * x * y * y * z * z - std::sqrt(170.8984375) * x * x * x * x * x * z * z * z * z + std::sqrt(2.4609375) * x * x * x * y * y * y * y * y * y - std::sqrt(246.09375) * x * x * x * y * y * y * y * z * z + std::sqrt(683.59375) * x * x * x * y * y * z * z * z * z + std::sqrt(17.5) * x * x * x * z * z * z * z * z * z + std::sqrt(0.3460693359375) * x * y * y * y * y * y * y * y * y - std::sqrt(88.59375) * x * y * y * y * y * y * y * z * z + std::sqrt(1538.0859375) * x * y * y * y * y * z * z * z * z - std::sqrt(157.5) * x * y * y * z * z * z * z * z * z) + e_1 * (-std::sqrt(15.380859375) * x * x * x * x * x * x * x + std::sqrt(15.380859375) * x * x * x * x * x * y * y - std::sqrt(246.09375) * x * x * x * x * x * z * z + std::sqrt(384.521484375) * x * x * x * y * y * y * y + std::sqrt(984.375) * x * x * x * y * y * z * z - std::sqrt(1750.0) * x * x * x * z * z * z * z + std::sqrt(138.427734375) * x * y * y * y * y * y * y + std::sqrt(2214.84375) * x * y * y * y * y * z * z + std::sqrt(15750.0) * x * y * y * z * z * z * z) + e_2 * (-std::sqrt(2214.84375) * x * x * x * x * x + std::sqrt(8859.375) * x * x * x * y * y - std::sqrt(35437.5) * x * x * x * z * z + std::sqrt(19933.59375) * x * y * y * y * y + std::sqrt(318937.5) * x * y * y * z * z) + e_3 * (-std::sqrt(35437.5) * x * x * x + std::sqrt(318937.5) * x * y * y);

        pc_77[k] = e_0 * (std::sqrt(1.5380859375) * x * x * x * x * x * x * x * x * z - std::sqrt(1.5380859375) * x * x * x * x * x * x * y * y * z - std::sqrt(133.984375) * x * x * x * x * x * x * z * z * z - std::sqrt(38.4521484375) * x * x * x * x * y * y * y * y * z + std::sqrt(535.9375) * x * x * x * x * y * y * z * z * z + std::sqrt(175.0) * x * x * x * x * z * z * z * z * z - std::sqrt(13.8427734375) * x * x * y * y * y * y * y * y * z + std::sqrt(1205.859375) * x * x * y * y * y * y * z * z * z - std::sqrt(1575.0) * x * x * y * y * z * z * z * z * z) + e_1 * (-std::sqrt(1.5380859375) * x * x * x * x * x * x * z + std::sqrt(186.1083984375) * x * x * x * x * y * y * z - std::sqrt(330.859375) * x * x * x * x * z * z * z + std::sqrt(124.5849609375) * x * x * y * y * y * y * z - std::sqrt(2460.9375) * x * x * y * y * z * z * z + std::sqrt(1575.0) * x * x * z * z * z * z * z - std::sqrt(13.8427734375) * y * y * y * y * y * y * z + std::sqrt(1205.859375) * y * y * y * y * z * z * z - std::sqrt(1575.0) * y * y * z * z * z * z * z) + e_2 * (-std::sqrt(885.9375) * x * x * x * x * z + std::sqrt(14175.0) * x * x * z * z * z + std::sqrt(885.9375) * y * y * y * y * z - std::sqrt(14175.0) * y * y * z * z * z) + e_3 * (std::sqrt(3543.75) * x * x * z - std::sqrt(3543.75) * y * y * z);

        pc_78[k] = e_0 * (std::sqrt(0.08544921875) * x * x * x * x * x * x * x * x * x - std::sqrt(0.341796875) * x * x * x * x * x * x * x * y * y - std::sqrt(16.748046875) * x * x * x * x * x * x * x * z * z - std::sqrt(1.3671875) * x * x * x * x * x * y * y * y * y + std::sqrt(150.732421875) * x * x * x * x * x * y * y * z * z + std::sqrt(196.875) * x * x * x * x * x * z * z * z * z + std::sqrt(0.341796875) * x * x * x * y * y * y * y * y * y + std::sqrt(16.748046875) * x * x * x * y * y * y * y * z * z - std::sqrt(3150.0) * x * x * x * y * y * z * z * z * z + std::sqrt(0.76904296875) * x * y * y * y * y * y * y * y * y - std::sqrt(150.732421875) * x * y * y * y * y * y * y * z * z + std::sqrt(1771.875) * x * y * y * y * y * z * z * z * z) + e_1 * (std::sqrt(34.1796875) * x * x * x * x * x * x * x - std::sqrt(49.21875) * x * x * x * x * x * y * y + std::sqrt(110.7421875) * x * x * x * x * x * z * z + std::sqrt(1.3671875) * x * x * x * y * y * y * y - std::sqrt(35880.46875) * x * x * x * y * y * z * z + std::sqrt(7087.5) * x * x * x * z * z * z * z + std::sqrt(196.875) * x * y * y * y * y * y * y + std::sqrt(6509.1796875) * x * y * y * y * y * z * z + std::sqrt(7087.5) * x * y * y * z * z * z * z) + e_2 * (std::sqrt(3986.71875) * x * x * x * x * x - std::sqrt(15946.875) * x * x * x * y * y + std::sqrt(44296.875) * x * x * x * z * z + std::sqrt(21705.46875) * x * y * y * y * y + std::sqrt(44296.875) * x * y * y * z * z + std::sqrt(28350.0) * x * z * z * z * z) + e_3 * (std::sqrt(56896.875) * x * x * x + std::sqrt(56896.875) * x * y * y + std::sqrt(381150.0) * x * z * z) + e_4 * (std::sqrt(154350.0) * x);

        pc_79[k] = e_0 * (-1.09375 * x * x * x * x * x * x * x * x * z + 5.46875 * x * x * x * x * x * x * y * y * z + 8.75 * x * x * x * x * x * x * z * z * z - 3.28125 * x * x * x * x * y * y * y * y * z - 52.5 * x * x * x * x * y * y * z * z * z - 9.84375 * x * x * y * y * y * y * y * y * z + 78.75 * x * x * y * y * y * y * z * z * z) + e_1 * (1.09375 * x * x * x * x * x * x * z - 95.15625 * x * x * x * x * y * y * z + 78.75 * x * x * x * x * z * z * z + 68.90625 * x * x * y * y * y * y * z + 157.5 * x * x * y * y * z * z * z - 9.84375 * y * y * y * y * y * y * z + 78.75 * y * y * y * y * z * z * z) + e_2 * (78.75 * x * x * x * x * z + 157.5 * x * x * y * y * z + 315.0 * x * x * z * z * z + 78.75 * y * y * y * y * z + 315.0 * y * y * z * z * z) + e_3 * (525.0 * x * x * z + 525.0 * y * y * z + 210.0 * z * z * z) + e_4 * (420.0 * z);
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

        pc_80[k] = e_0 * (-std::sqrt(0.1495361328125) * x * x * x * x * x * x * x * x * x + std::sqrt(9.5703125) * x * x * x * x * x * x * x * y * y + std::sqrt(9.5703125) * x * x * x * x * x * x * x * z * z - std::sqrt(14.95361328125) * x * x * x * x * x * y * y * y * y - std::sqrt(775.1953125) * x * x * x * x * x * y * y * z * z - std::sqrt(38.28125) * x * x * x * y * y * y * y * y * y + std::sqrt(3454.8828125) * x * x * x * y * y * y * y * z * z + std::sqrt(1.3458251953125) * x * y * y * y * y * y * y * y * y - std::sqrt(86.1328125) * x * y * y * y * y * y * y * z * z) + e_1 * (-std::sqrt(59.814453125) * x * x * x * x * x * x * x + std::sqrt(193.798828125) * x * x * x * x * x * y * y + std::sqrt(1378.125) * x * x * x * x * x * z * z - std::sqrt(5285.205078125) * x * x * x * y * y * y * y + std::sqrt(5512.5) * x * x * x * y * y * z * z + std::sqrt(21.533203125) * x * y * y * y * y * y * y + std::sqrt(1378.125) * x * y * y * y * y * z * z) + e_2 * (-std::sqrt(3100.78125) * x * x * x * x * x - std::sqrt(12403.125) * x * x * x * y * y + std::sqrt(49612.5) * x * x * x * z * z - std::sqrt(3100.78125) * x * y * y * y * y + std::sqrt(49612.5) * x * y * y * z * z) + e_3 * (-std::sqrt(22050.0) * x * x * x - std::sqrt(22050.0) * x * y * y + std::sqrt(88200.0) * x * z * z) + e_4 * (-std::sqrt(5512.5) * x);

        pc_81[k] = e_0 * (6.5625 * x * x * x * x * x * x * x * y * z - 45.9375 * x * x * x * x * x * y * y * y * z + 45.9375 * x * x * x * y * y * y * y * y * z - 6.5625 * x * y * y * y * y * y * y * y * z);

        pc_82[k] = e_0 * (std::sqrt(193.798828125) * x * x * x * x * x * x * y * z * z - std::sqrt(7773.486328125) * x * x * x * x * y * y * y * z * z + std::sqrt(1744.189453125) * x * x * y * y * y * y * y * z * z - std::sqrt(21.533203125) * y * y * y * y * y * y * y * z * z) + e_1 * (std::sqrt(193.798828125) * x * x * x * x * x * x * y - std::sqrt(7773.486328125) * x * x * x * x * y * y * y - std::sqrt(3100.78125) * x * x * x * x * y * z * z + std::sqrt(1744.189453125) * x * x * y * y * y * y * y - std::sqrt(12403.125) * x * x * y * y * y * z * z - std::sqrt(21.533203125) * y * y * y * y * y * y * y - std::sqrt(3100.78125) * y * y * y * y * y * z * z) + e_2 * (-std::sqrt(3100.78125) * x * x * x * x * y - std::sqrt(12403.125) * x * x * y * y * y - std::sqrt(111628.125) * x * x * y * z * z - std::sqrt(3100.78125) * y * y * y * y * y - std::sqrt(111628.125) * y * y * y * z * z) + e_3 * (-std::sqrt(111628.125) * x * x * y - std::sqrt(111628.125) * y * y * y - std::sqrt(198450.0) * y * z * z) + e_4 * (-std::sqrt(198450.0) * y);

        pc_83[k] = e_0 * (-std::sqrt(6.15234375) * x * x * x * x * x * x * x * y * z + std::sqrt(153.80859375) * x * x * x * x * x * y * y * y * z + std::sqrt(221.484375) * x * x * x * x * x * y * z * z * z + std::sqrt(153.80859375) * x * x * x * y * y * y * y * y * z - std::sqrt(7973.4375) * x * x * x * y * y * y * z * z * z - std::sqrt(6.15234375) * x * y * y * y * y * y * y * y * z + std::sqrt(221.484375) * x * y * y * y * y * y * z * z * z) + e_1 * (std::sqrt(885.9375) * x * x * x * x * x * y * z - std::sqrt(393.75) * x * x * x * y * y * y * z - std::sqrt(14175.0) * x * x * x * y * z * z * z + std::sqrt(885.9375) * x * y * y * y * y * y * z - std::sqrt(14175.0) * x * y * y * y * z * z * z) + e_2 * (-std::sqrt(3543.75) * x * x * x * y * z - std::sqrt(3543.75) * x * y * y * y * z - std::sqrt(127575.0) * x * y * z * z * z) + e_3 * (-std::sqrt(226800.0) * x * y * z);
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

        pc_84[k] = e_0 * (-std::sqrt(27.685546875) * x * x * x * x * x * x * y * z * z + std::sqrt(692.138671875) * x * x * x * x * y * y * y * z * z + std::sqrt(49.21875) * x * x * x * x * y * z * z * z * z + std::sqrt(692.138671875) * x * x * y * y * y * y * y * z * z - std::sqrt(1771.875) * x * x * y * y * y * z * z * z * z - std::sqrt(27.685546875) * y * y * y * y * y * y * y * z * z + std::sqrt(49.21875) * y * y * y * y * y * z * z * z * z) + e_1 * (-std::sqrt(27.685546875) * x * x * x * x * x * x * y + std::sqrt(692.138671875) * x * x * x * x * y * y * y + std::sqrt(1771.875) * x * x * x * x * y * z * z + std::sqrt(692.138671875) * x * x * y * y * y * y * y + std::sqrt(28350.0) * x * x * y * y * y * z * z - std::sqrt(7087.5) * x * x * y * z * z * z * z - std::sqrt(27.685546875) * y * y * y * y * y * y * y - std::sqrt(1771.875) * y * y * y * y * y * z * z + std::sqrt(787.5) * y * y * y * z * z * z * z) + e_2 * (std::sqrt(442.96875) * x * x * x * x * y + std::sqrt(86821.875) * x * x * y * y * y + std::sqrt(15946.875) * x * x * y * z * z - std::sqrt(3986.71875) * y * y * y * y * y - std::sqrt(1771.875) * y * y * y * z * z) + e_3 * (std::sqrt(143521.875) * x * x * y - std::sqrt(15946.875) * y * y * y);

        pc_85[k] = e_0 * (std::sqrt(0.692138671875) * x * x * x * x * x * x * x * x * z - std::sqrt(11.07421875) * x * x * x * x * x * x * y * y * z - std::sqrt(44.296875) * x * x * x * x * x * x * z * z * z - std::sqrt(69.2138671875) * x * x * x * x * y * y * y * y * z + std::sqrt(1107.421875) * x * x * x * x * y * y * z * z * z + std::sqrt(4.921875) * x * x * x * x * z * z * z * z * z - std::sqrt(11.07421875) * x * x * y * y * y * y * y * y * z + std::sqrt(1107.421875) * x * x * y * y * y * y * z * z * z - std::sqrt(177.1875) * x * x * y * y * z * z * z * z * z + std::sqrt(0.692138671875) * y * y * y * y * y * y * y * y * z - std::sqrt(44.296875) * y * y * y * y * y * y * z * z * z + std::sqrt(4.921875) * y * y * y * y * z * z * z * z * z) + e_1 * (-std::sqrt(1968.75) * x * x * x * x * z * z * z + std::sqrt(70875.0) * x * x * y * y * z * z * z - std::sqrt(1968.75) * y * y * y * y * z * z * z) + e_2 * (-std::sqrt(4429.6875) * x * x * x * x * z + std::sqrt(159468.75) * x * x * y * y * z - std::sqrt(4429.6875) * y * y * y * y * z);

        pc_86[k] = e_0 * (-std::sqrt(27.685546875) * x * x * x * x * x * x * x * z * z + std::sqrt(692.138671875) * x * x * x * x * x * y * y * z * z + std::sqrt(49.21875) * x * x * x * x * x * z * z * z * z + std::sqrt(692.138671875) * x * x * x * y * y * y * y * z * z - std::sqrt(1771.875) * x * x * x * y * y * z * z * z * z - std::sqrt(27.685546875) * x * y * y * y * y * y * y * z * z + std::sqrt(49.21875) * x * y * y * y * y * z * z * z * z) + e_1 * (-std::sqrt(27.685546875) * x * x * x * x * x * x * x + std::sqrt(692.138671875) * x * x * x * x * x * y * y - std::sqrt(1771.875) * x * x * x * x * x * z * z + std::sqrt(692.138671875) * x * x * x * y * y * y * y + std::sqrt(28350.0) * x * x * x * y * y * z * z + std::sqrt(787.5) * x * x * x * z * z * z * z - std::sqrt(27.685546875) * x * y * y * y * y * y * y + std::sqrt(1771.875) * x * y * y * y * y * z * z - std::sqrt(7087.5) * x * y * y * z * z * z * z) + e_2 * (-std::sqrt(3986.71875) * x * x * x * x * x + std::sqrt(86821.875) * x * x * x * y * y - std::sqrt(1771.875) * x * x * x * z * z + std::sqrt(442.96875) * x * y * y * y * y + std::sqrt(15946.875) * x * y * y * z * z) + e_3 * (-std::sqrt(15946.875) * x * x * x + std::sqrt(143521.875) * x * y * y);

        pc_87[k] = e_0 * (-std::sqrt(1.5380859375) * x * x * x * x * x * x * x * x * z + std::sqrt(55.37109375) * x * x * x * x * x * x * y * y * z + std::sqrt(55.37109375) * x * x * x * x * x * x * z * z * z - std::sqrt(2713.18359375) * x * x * x * x * y * y * z * z * z - std::sqrt(55.37109375) * x * x * y * y * y * y * y * y * z + std::sqrt(2713.18359375) * x * x * y * y * y * y * z * z * z + std::sqrt(1.5380859375) * y * y * y * y * y * y * y * y * z - std::sqrt(55.37109375) * y * y * y * y * y * y * z * z * z) + e_1 * (-std::sqrt(24.609375) * x * x * x * x * x * x * z - std::sqrt(1993.359375) * x * x * x * x * y * y * z + std::sqrt(3543.75) * x * x * x * x * z * z * z + std::sqrt(1993.359375) * x * x * y * y * y * y * z + std::sqrt(24.609375) * y * y * y * y * y * y * z - std::sqrt(3543.75) * y * y * y * y * z * z * z) + e_2 * (std::sqrt(885.9375) * x * x * x * x * z + std::sqrt(31893.75) * x * x * z * z * z - std::sqrt(885.9375) * y * y * y * y * z - std::sqrt(31893.75) * y * y * z * z * z) + e_3 * (std::sqrt(56700.0) * x * x * z - std::sqrt(56700.0) * y * y * z);
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

        pc_88[k] = e_0 * (std::sqrt(21.533203125) * x * x * x * x * x * x * x * z * z - std::sqrt(1744.189453125) * x * x * x * x * x * y * y * z * z + std::sqrt(7773.486328125) * x * x * x * y * y * y * y * z * z - std::sqrt(193.798828125) * x * y * y * y * y * y * y * z * z) + e_1 * (std::sqrt(21.533203125) * x * x * x * x * x * x * x - std::sqrt(1744.189453125) * x * x * x * x * x * y * y + std::sqrt(3100.78125) * x * x * x * x * x * z * z + std::sqrt(7773.486328125) * x * x * x * y * y * y * y + std::sqrt(12403.125) * x * x * x * y * y * z * z - std::sqrt(193.798828125) * x * y * y * y * y * y * y + std::sqrt(3100.78125) * x * y * y * y * y * z * z) + e_2 * (std::sqrt(3100.78125) * x * x * x * x * x + std::sqrt(12403.125) * x * x * x * y * y + std::sqrt(111628.125) * x * x * x * z * z + std::sqrt(3100.78125) * x * y * y * y * y + std::sqrt(111628.125) * x * y * y * z * z) + e_3 * (std::sqrt(111628.125) * x * x * x + std::sqrt(111628.125) * x * y * y + std::sqrt(198450.0) * x * z * z) + e_4 * (std::sqrt(198450.0) * x);

        pc_89[k] = e_0 * (1.640625 * x * x * x * x * x * x * x * x * z - 19.6875 * x * x * x * x * x * x * y * y * z + 62.34375 * x * x * x * x * y * y * y * y * z - 19.6875 * x * x * y * y * y * y * y * y * z + 1.640625 * y * y * y * y * y * y * y * y * z) + e_1 * (26.25 * x * x * x * x * x * x * z + 78.75 * x * x * x * x * y * y * z + 78.75 * x * x * y * y * y * y * z + 26.25 * y * y * y * y * y * y * z) + e_2 * (236.25 * x * x * x * x * z + 472.5 * x * x * y * y * z + 236.25 * y * y * y * y * z) + e_3 * (630.0 * x * x * z + 630.0 * y * y * z) + e_4 * (315.0 * z);

        pc_90[k] = e_0 * (std::sqrt(4.306640625) * x * x * x * x * x * x * x * x * y - std::sqrt(521.103515625) * x * x * x * x * x * x * y * y * y + std::sqrt(968.994140625) * x * x * x * x * y * y * y * y * y - std::sqrt(107.666015625) * x * x * y * y * y * y * y * y * y) + e_1 * (-std::sqrt(107.666015625) * x * x * x * x * x * x * y - std::sqrt(968.994140625) * x * x * x * x * y * y * y - std::sqrt(968.994140625) * x * x * y * y * y * y * y - std::sqrt(107.666015625) * y * y * y * y * y * y * y) + e_2 * (-std::sqrt(15503.90625) * x * x * x * x * y - std::sqrt(62015.625) * x * x * y * y * y - std::sqrt(15503.90625) * y * y * y * y * y) + e_3 * (-std::sqrt(248062.5) * x * x * y - std::sqrt(248062.5) * y * y * y) + e_4 * (-std::sqrt(248062.5) * y);

        pc_91[k] = e_0 * (std::sqrt(19.3798828125) * x * x * x * x * x * x * x * y * z - std::sqrt(2069.3408203125) * x * x * x * x * x * y * y * y * z + std::sqrt(1345.8251953125) * x * x * x * y * y * y * y * y * z - std::sqrt(53.8330078125) * x * y * y * y * y * y * y * y * z) + e_1 * (-std::sqrt(1937.98828125) * x * x * x * x * x * y * z - std::sqrt(7751.953125) * x * x * x * y * y * y * z - std::sqrt(1937.98828125) * x * y * y * y * y * y * z) + e_2 * (-std::sqrt(124031.25) * x * x * x * y * z - std::sqrt(124031.25) * x * y * y * y * z) + e_3 * (-std::sqrt(496125.0) * x * y * z);
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

        pc_92[k] = e_0 * (-std::sqrt(0.615234375) * x * x * x * x * x * x * x * x * y + std::sqrt(49.833984375) * x * x * x * x * x * x * y * y * y + std::sqrt(22.1484375) * x * x * x * x * x * x * y * z * z + std::sqrt(15.380859375) * x * x * x * x * y * y * y * y * y - std::sqrt(2214.84375) * x * x * x * x * y * y * y * z * z - std::sqrt(15.380859375) * x * x * y * y * y * y * y * y * y + std::sqrt(553.7109375) * x * x * y * y * y * y * y * z * z) + e_1 * (std::sqrt(15.380859375) * x * x * x * x * x * x * y + std::sqrt(9613.037109375) * x * x * x * x * y * y * y - std::sqrt(4983.3984375) * x * x * x * x * y * z * z - std::sqrt(1245.849609375) * x * x * y * y * y * y * y - std::sqrt(2214.84375) * x * x * y * y * y * z * z - std::sqrt(15.380859375) * y * y * y * y * y * y * y + std::sqrt(553.7109375) * y * y * y * y * y * z * z) + e_2 * (std::sqrt(19933.59375) * x * x * x * x * y + std::sqrt(8859.375) * x * x * y * y * y - std::sqrt(79734.375) * x * x * y * z * z - std::sqrt(2214.84375) * y * y * y * y * y + std::sqrt(8859.375) * y * y * y * z * z) + e_3 * (std::sqrt(79734.375) * x * x * y - std::sqrt(8859.375) * y * y * y);

        pc_93[k] = e_0 * (-std::sqrt(2.7685546875) * x * x * x * x * x * x * x * y * z + std::sqrt(224.2529296875) * x * x * x * x * x * y * y * y * z + std::sqrt(4.921875) * x * x * x * x * x * y * z * z * z + std::sqrt(69.2138671875) * x * x * x * y * y * y * y * y * z - std::sqrt(492.1875) * x * x * x * y * y * y * z * z * z - std::sqrt(69.2138671875) * x * y * y * y * y * y * y * y * z + std::sqrt(123.046875) * x * y * y * y * y * y * z * z * z) + e_1 * (std::sqrt(276.85546875) * x * x * x * x * x * y * z + std::sqrt(27685.546875) * x * x * x * y * y * y * z - std::sqrt(1968.75) * x * x * x * y * z * z * z - std::sqrt(13565.91796875) * x * y * y * y * y * y * z + std::sqrt(1968.75) * x * y * y * y * z * z * z) + e_2 * (std::sqrt(70875.0) * x * x * x * y * z - std::sqrt(70875.0) * x * y * y * y * z);

        pc_94[k] = e_0 * (std::sqrt(0.0692138671875) * x * x * x * x * x * x * x * x * x - std::sqrt(4.4296875) * x * x * x * x * x * x * x * y * y - std::sqrt(4.4296875) * x * x * x * x * x * x * x * z * z - std::sqrt(13.56591796875) * x * x * x * x * x * y * y * y * y + std::sqrt(358.8046875) * x * x * x * x * x * y * y * z * z + std::sqrt(0.4921875) * x * x * x * x * x * z * z * z * z + std::sqrt(110.7421875) * x * x * x * y * y * y * y * z * z - std::sqrt(49.21875) * x * x * x * y * y * z * z * z * z + std::sqrt(1.7303466796875) * x * y * y * y * y * y * y * y * y - std::sqrt(110.7421875) * x * y * y * y * y * y * y * z * z + std::sqrt(12.3046875) * x * y * y * y * y * z * z * z * z) + e_1 * (std::sqrt(27.685546875) * x * x * x * x * x * x * x - std::sqrt(2242.529296875) * x * x * x * x * x * y * y - std::sqrt(442.96875) * x * x * x * x * x * z * z - std::sqrt(692.138671875) * x * x * x * y * y * y * y + std::sqrt(44296.875) * x * x * x * y * y * z * z + std::sqrt(692.138671875) * x * y * y * y * y * y * y - std::sqrt(11074.21875) * x * y * y * y * y * z * z) + e_2 * (std::sqrt(442.96875) * x * x * x * x * x - std::sqrt(44296.875) * x * x * x * y * y + std::sqrt(11074.21875) * x * y * y * y * y);

        pc_95[k] = e_0 * (-std::sqrt(2.7685546875) * x * x * x * x * x * x * x * x * z + std::sqrt(224.2529296875) * x * x * x * x * x * x * y * y * z + std::sqrt(4.921875) * x * x * x * x * x * x * z * z * z + std::sqrt(69.2138671875) * x * x * x * x * y * y * y * y * z - std::sqrt(492.1875) * x * x * x * x * y * y * z * z * z - std::sqrt(69.2138671875) * x * x * y * y * y * y * y * y * z + std::sqrt(123.046875) * x * x * y * y * y * y * z * z * z) + e_1 * (-std::sqrt(622.9248046875) * x * x * x * x * x * x * z + std::sqrt(43258.6669921875) * x * x * x * x * y * y * z + std::sqrt(123.046875) * x * x * x * x * z * z * z - std::sqrt(1730.3466796875) * x * x * y * y * y * y * z - std::sqrt(4429.6875) * x * x * y * y * z * z * z - std::sqrt(69.2138671875) * y * y * y * y * y * y * z + std::sqrt(123.046875) * y * y * y * y * z * z * z) + e_2 * (-std::sqrt(4429.6875) * x * x * x * x * z + std::sqrt(159468.75) * x * x * y * y * z - std::sqrt(4429.6875) * y * y * y * y * z);
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

        pc_96[k] = e_0 * (-std::sqrt(0.15380859375) * x * x * x * x * x * x * x * x * x + std::sqrt(15.380859375) * x * x * x * x * x * x * x * y * y + std::sqrt(5.537109375) * x * x * x * x * x * x * x * z * z - std::sqrt(2.4609375) * x * x * x * x * x * y * y * y * y - std::sqrt(669.990234375) * x * x * x * x * x * y * y * z * z - std::sqrt(15.380859375) * x * x * x * y * y * y * y * y * y + std::sqrt(1245.849609375) * x * x * x * y * y * y * y * z * z + std::sqrt(3.84521484375) * x * y * y * y * y * y * y * y * y - std::sqrt(138.427734375) * x * y * y * y * y * y * y * z * z) + e_1 * (-std::sqrt(61.5234375) * x * x * x * x * x * x * x + std::sqrt(2214.84375) * x * x * x * x * x * y * y + std::sqrt(553.7109375) * x * x * x * x * x * z * z - std::sqrt(1538.0859375) * x * x * x * y * y * y * y - std::sqrt(2214.84375) * x * x * x * y * y * z * z + std::sqrt(984.375) * x * y * y * y * y * y * y - std::sqrt(4983.3984375) * x * y * y * y * y * z * z) + e_2 * (-std::sqrt(2214.84375) * x * x * x * x * x + std::sqrt(8859.375) * x * x * x * y * y + std::sqrt(8859.375) * x * x * x * z * z + std::sqrt(19933.59375) * x * y * y * y * y - std::sqrt(79734.375) * x * y * y * z * z) + e_3 * (-std::sqrt(8859.375) * x * x * x + std::sqrt(79734.375) * x * y * y);

        pc_97[k] = e_0 * (std::sqrt(2.1533203125) * x * x * x * x * x * x * x * x * z - std::sqrt(363.9111328125) * x * x * x * x * x * x * y * y * z + std::sqrt(2637.8173828125) * x * x * x * x * y * y * y * y * z - std::sqrt(484.4970703125) * x * x * y * y * y * y * y * y * z) + e_1 * (std::sqrt(484.4970703125) * x * x * x * x * x * x * z + std::sqrt(484.4970703125) * x * x * x * x * y * y * z - std::sqrt(484.4970703125) * x * x * y * y * y * y * z - std::sqrt(484.4970703125) * y * y * y * y * y * y * z) + e_2 * (std::sqrt(31007.8125) * x * x * x * x * z - std::sqrt(31007.8125) * y * y * y * y * z) + e_3 * (std::sqrt(124031.25) * x * x * z - std::sqrt(124031.25) * y * y * z);

        pc_98[k] = e_0 * (std::sqrt(0.2691650390625) * x * x * x * x * x * x * x * x * x - std::sqrt(68.90625) * x * x * x * x * x * x * x * y * y + std::sqrt(1172.48291015625) * x * x * x * x * x * y * y * y * y - std::sqrt(430.6640625) * x * x * x * y * y * y * y * y * y + std::sqrt(6.7291259765625) * x * y * y * y * y * y * y * y * y) + e_1 * (std::sqrt(107.666015625) * x * x * x * x * x * x * x + std::sqrt(968.994140625) * x * x * x * x * x * y * y + std::sqrt(968.994140625) * x * x * x * y * y * y * y + std::sqrt(107.666015625) * x * y * y * y * y * y * y) + e_2 * (std::sqrt(15503.90625) * x * x * x * x * x + std::sqrt(62015.625) * x * x * x * y * y + std::sqrt(15503.90625) * x * y * y * y * y) + e_3 * (std::sqrt(248062.5) * x * x * x + std::sqrt(248062.5) * x * y * y) + e_4 * (std::sqrt(248062.5) * x);
    }

    // NOTE: the atom pairs beyond the reach of every pair of primitives have no
    // contribution and are set to zero.

    for (size_t m = 0; m < 99; m++)
    {
        auto *pv = values + m * nvalues;

        std::fill(pv + nmax, pv + nvalues, 0.0);
    }
}

}  // namespace simdovl
