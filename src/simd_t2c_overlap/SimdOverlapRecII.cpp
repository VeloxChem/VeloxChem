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



#include "SimdOverlapRecII.hpp"

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
compute_ii_overlap(double               *values,
                   const size_t          nvalues,
                   const CBasisFunction &bra,
                   const CBasisFunction &ket,
                   const CSimdMatrix    &coordinates,
                   const double          threshold) -> void
{
    if ((bra.get_angular_momentum() != 6) || (ket.get_angular_momentum() != 6))
    {
        errors::assertMsgCritical(
            false, std::string("SimdOverlapRecII.compute_ii_overlap: Basis functions must be of angular momenta six and six"));
    }

    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("SimdOverlapRecII.compute_ii_overlap: Number of values exceeds number of atom pairs"));
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

    auto buffer = simdfunc::make_primitive_buffer(dimensions, 7);

    if (buffer.number_of_columns() == 0)
    {
        std::fill(values, values + 169 * nvalues, 0.0);

        return;
    }

    const auto nmax = buffer.number_of_columns();

    auto *pe_0 = buffer.data(0);
    auto *pe_1 = buffer.data(1);
    auto *pe_2 = buffer.data(2);
    auto *pe_3 = buffer.data(3);
    auto *pe_4 = buffer.data(4);
    auto *pe_5 = buffer.data(5);
    auto *pe_6 = buffer.data(6);

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

        const auto f_0 = fbase * fal * fal * fal * fal * fal * fal * fbe * fbe * fbe * fbe * fbe * fbe;

        const auto f_1 = fbase * fal * fal * fal * fal * fal * fbe * fbe * fbe * fbe * fbe * fh;

        const auto f_2 = fbase * fal * fal * fal * fal * fbe * fbe * fbe * fbe * fh * fh;

        const auto f_3 = fbase * fal * fal * fal * fbe * fbe * fbe * fh * fh * fh;

        const auto f_4 = fbase * fal * fal * fbe * fbe * fh * fh * fh * fh;

        const auto f_5 = fbase * fal * fbe * fh * fh * fh * fh * fh;

        const auto f_6 = fbase * fh * fh * fh * fh * fh * fh;

        // NOTE: the exponential depends on the pair of primitives alone, so it is
        // evaluated once and shared by the prefactors of all terms.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ab_2 : simd::cache_line_size())
        for (size_t k = 0; k < ncols; k++)
        {
            const auto fss = std::exp(-fmu * ab_2[k]);

            pe_0[k] += f_0 * fss;
            pe_1[k] += f_1 * fss;
            pe_2[k] += f_2 * fss;
            pe_3[k] += f_3 * fss;
            pe_4[k] += f_4 * fss;
            pe_5[k] += f_5 * fss;
            pe_6[k] += f_6 * fss;
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
    auto *pc_13 = values + 14 * nvalues;
    auto *pc_14 = values + 15 * nvalues;
    auto *pc_15 = values + 16 * nvalues;
    auto *pc_16 = values + 17 * nvalues;
    auto *pc_17 = values + 18 * nvalues;
    auto *pc_18 = values + 19 * nvalues;
    auto *pc_19 = values + 20 * nvalues;
    auto *pc_20 = values + 21 * nvalues;
    auto *pc_21 = values + 22 * nvalues;
    auto *pc_22 = values + 23 * nvalues;
    auto *pc_23 = values + 24 * nvalues;
    auto *pc_24 = values + 25 * nvalues;
    auto *pc_25 = values + 28 * nvalues;
    auto *pc_26 = values + 29 * nvalues;
    auto *pc_27 = values + 30 * nvalues;
    auto *pc_28 = values + 31 * nvalues;
    auto *pc_29 = values + 32 * nvalues;
    auto *pc_30 = values + 33 * nvalues;
    auto *pc_31 = values + 34 * nvalues;
    auto *pc_32 = values + 35 * nvalues;
    auto *pc_33 = values + 36 * nvalues;
    auto *pc_34 = values + 37 * nvalues;
    auto *pc_35 = values + 38 * nvalues;
    auto *pc_36 = values + 42 * nvalues;
    auto *pc_37 = values + 43 * nvalues;
    auto *pc_38 = values + 44 * nvalues;
    auto *pc_39 = values + 45 * nvalues;
    auto *pc_40 = values + 46 * nvalues;
    auto *pc_41 = values + 47 * nvalues;
    auto *pc_42 = values + 48 * nvalues;
    auto *pc_43 = values + 49 * nvalues;
    auto *pc_44 = values + 50 * nvalues;
    auto *pc_45 = values + 51 * nvalues;
    auto *pc_46 = values + 56 * nvalues;
    auto *pc_47 = values + 57 * nvalues;
    auto *pc_48 = values + 58 * nvalues;
    auto *pc_49 = values + 59 * nvalues;
    auto *pc_50 = values + 60 * nvalues;
    auto *pc_51 = values + 61 * nvalues;
    auto *pc_52 = values + 62 * nvalues;
    auto *pc_53 = values + 63 * nvalues;
    auto *pc_54 = values + 64 * nvalues;
    auto *pc_55 = values + 70 * nvalues;
    auto *pc_56 = values + 71 * nvalues;
    auto *pc_57 = values + 72 * nvalues;
    auto *pc_58 = values + 73 * nvalues;
    auto *pc_59 = values + 74 * nvalues;
    auto *pc_60 = values + 75 * nvalues;
    auto *pc_61 = values + 76 * nvalues;
    auto *pc_62 = values + 77 * nvalues;
    auto *pc_63 = values + 84 * nvalues;
    auto *pc_64 = values + 85 * nvalues;
    auto *pc_65 = values + 86 * nvalues;
    auto *pc_66 = values + 87 * nvalues;
    auto *pc_67 = values + 88 * nvalues;
    auto *pc_68 = values + 89 * nvalues;
    auto *pc_69 = values + 90 * nvalues;
    auto *pc_70 = values + 98 * nvalues;
    auto *pc_71 = values + 99 * nvalues;
    auto *pc_72 = values + 100 * nvalues;
    auto *pc_73 = values + 101 * nvalues;
    auto *pc_74 = values + 102 * nvalues;
    auto *pc_75 = values + 103 * nvalues;
    auto *pc_76 = values + 112 * nvalues;
    auto *pc_77 = values + 113 * nvalues;
    auto *pc_78 = values + 114 * nvalues;
    auto *pc_79 = values + 115 * nvalues;
    auto *pc_80 = values + 116 * nvalues;
    auto *pc_81 = values + 126 * nvalues;
    auto *pc_82 = values + 127 * nvalues;
    auto *pc_83 = values + 128 * nvalues;
    auto *pc_84 = values + 129 * nvalues;
    auto *pc_85 = values + 140 * nvalues;
    auto *pc_86 = values + 141 * nvalues;
    auto *pc_87 = values + 142 * nvalues;
    auto *pc_88 = values + 154 * nvalues;
    auto *pc_89 = values + 155 * nvalues;
    auto *pc_90 = values + 168 * nvalues;

    // NOTE: the components are formed in 88 loops, as the vectorizer runs out
    // of registers with all of them in one. Only the prefactors and the vector
    // between the atoms are loaded by more than one loop.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ab_x, ab_y, ab_z : simd::cache_line_size())
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
        const auto e_6 = pe_6[k];

        pc_0[k] = e_0 * (16.2421875 * x * x * x * x * x * x * x * x * x * x * y * y - 108.28125 * x * x * x * x * x * x * x * x * y * y * y * y + 212.953125 * x * x * x * x * x * x * y * y * y * y * y * y - 108.28125 * x * x * x * x * y * y * y * y * y * y * y * y + 16.2421875 * x * x * y * y * y * y * y * y * y * y * y * y) + e_1 * (16.2421875 * x * x * x * x * x * x * x * x * x * x + 81.2109375 * x * x * x * x * x * x * x * x * y * y + 162.421875 * x * x * x * x * x * x * y * y * y * y + 162.421875 * x * x * x * x * y * y * y * y * y * y + 81.2109375 * x * x * y * y * y * y * y * y * y * y + 16.2421875 * y * y * y * y * y * y * y * y * y * y) + e_2 * (406.0546875 * x * x * x * x * x * x * x * x + 1624.21875 * x * x * x * x * x * x * y * y + 2436.328125 * x * x * x * x * y * y * y * y + 1624.21875 * x * x * y * y * y * y * y * y + 406.0546875 * y * y * y * y * y * y * y * y) + e_3 * (4331.25 * x * x * x * x * x * x + 12993.75 * x * x * x * x * y * y + 12993.75 * x * x * y * y * y * y + 4331.25 * y * y * y * y * y * y) + e_4 * (19490.625 * x * x * x * x + 38981.25 * x * x * y * y + 19490.625 * y * y * y * y) + e_5 * (31185.0 * x * x + 31185.0 * y * y) + e_6 * (10395.0);

        pc_1[k] = e_0 * (std::sqrt(2198.4054565429688) * x * x * x * x * x * x * x * x * x * y * y * z - std::sqrt(62532.421875) * x * x * x * x * x * x * x * y * y * y * y * z + std::sqrt(136047.10034179688) * x * x * x * x * x * y * y * y * y * y * y * z - std::sqrt(15633.10546875) * x * x * x * y * y * y * y * y * y * y * y * z + std::sqrt(87.93621826171875) * x * y * y * y * y * y * y * y * y * y * y * z) + e_1 * (std::sqrt(2198.4054565429688) * x * x * x * x * x * x * x * x * x * z + std::sqrt(35174.4873046875) * x * x * x * x * x * x * x * y * y * z + std::sqrt(79142.59643554688) * x * x * x * x * x * y * y * y * y * z + std::sqrt(35174.4873046875) * x * x * x * y * y * y * y * y * y * z + std::sqrt(2198.4054565429688) * x * y * y * y * y * y * y * y * y * z) + e_2 * (std::sqrt(879362.1826171875) * x * x * x * x * x * x * x * z + std::sqrt(7914259.6435546875) * x * x * x * x * x * y * y * z + std::sqrt(7914259.6435546875) * x * x * x * y * y * y * y * z + std::sqrt(879362.1826171875) * x * y * y * y * y * y * y * z) + e_3 * (std::sqrt(56279179.6875) * x * x * x * x * x * z + std::sqrt(225116718.75) * x * x * x * y * y * z + std::sqrt(56279179.6875) * x * y * y * y * y * z) + e_4 * (std::sqrt(506512617.1875) * x * x * x * z + std::sqrt(506512617.1875) * x * y * y * z) + e_5 * (std::sqrt(324168075.0) * x * z);
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

        pc_2[k] = e_0 * (-std::sqrt(63.95361328125) * x * x * x * x * x * x * x * x * x * x * y * y + std::sqrt(710.595703125) * x * x * x * x * x * x * x * x * y * y * y * y + std::sqrt(6395.361328125) * x * x * x * x * x * x * x * x * y * y * z * z - std::sqrt(120090.673828125) * x * x * x * x * x * x * y * y * y * y * z * z - std::sqrt(710.595703125) * x * x * x * x * y * y * y * y * y * y * y * y + std::sqrt(120090.673828125) * x * x * x * x * y * y * y * y * y * y * z * z + std::sqrt(63.95361328125) * x * x * y * y * y * y * y * y * y * y * y * y - std::sqrt(6395.361328125) * x * x * y * y * y * y * y * y * y * y * z * z) + e_1 * (-std::sqrt(63.95361328125) * x * x * x * x * x * x * x * x * x * x - std::sqrt(14389.56298828125) * x * x * x * x * x * x * x * x * y * y + std::sqrt(6395.361328125) * x * x * x * x * x * x * x * x * z * z + std::sqrt(159884.033203125) * x * x * x * x * x * x * y * y * y * y + std::sqrt(25581.4453125) * x * x * x * x * x * x * y * y * z * z - std::sqrt(159884.033203125) * x * x * x * x * y * y * y * y * y * y + std::sqrt(14389.56298828125) * x * x * y * y * y * y * y * y * y * y - std::sqrt(25581.4453125) * x * x * y * y * y * y * y * y * z * z + std::sqrt(63.95361328125) * y * y * y * y * y * y * y * y * y * y - std::sqrt(6395.361328125) * y * y * y * y * y * y * y * y * z * z) + e_2 * (-std::sqrt(39971.00830078125) * x * x * x * x * x * x * x * x - std::sqrt(159884.033203125) * x * x * x * x * x * x * y * y + std::sqrt(1438956.298828125) * x * x * x * x * x * x * z * z + std::sqrt(1438956.298828125) * x * x * x * x * y * y * z * z + std::sqrt(159884.033203125) * x * x * y * y * y * y * y * y - std::sqrt(1438956.298828125) * x * x * y * y * y * y * z * z + std::sqrt(39971.00830078125) * y * y * y * y * y * y * y * y - std::sqrt(1438956.298828125) * y * y * y * y * y * y * z * z) + e_3 * (-std::sqrt(2558144.53125) * x * x * x * x * x * x - std::sqrt(2558144.53125) * x * x * x * x * y * y + std::sqrt(40930312.5) * x * x * x * x * z * z + std::sqrt(2558144.53125) * x * x * y * y * y * y + std::sqrt(2558144.53125) * y * y * y * y * y * y - std::sqrt(40930312.5) * y * y * y * y * z * z) + e_4 * (-std::sqrt(23023300.78125) * x * x * x * x + std::sqrt(92093203.125) * x * x * z * z + std::sqrt(23023300.78125) * y * y * y * y - std::sqrt(92093203.125) * y * y * z * z) + e_5 * (-std::sqrt(14734912.5) * x * x + std::sqrt(14734912.5) * y * y);
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

        pc_3[k] = e_0 * (-std::sqrt(1079.2172241210938) * x * x * x * x * x * x * x * x * x * y * y * z + std::sqrt(7674.43359375) * x * x * x * x * x * x * x * y * y * y * y * z + std::sqrt(7674.43359375) * x * x * x * x * x * x * x * y * y * z * z * z + std::sqrt(2611.439208984375) * x * x * x * x * x * y * y * y * y * y * y * z - std::sqrt(103178.49609375) * x * x * x * x * x * y * y * y * y * z * z * z - std::sqrt(3410.859375) * x * x * x * y * y * y * y * y * y * y * y * z + std::sqrt(34203.33984375) * x * x * x * y * y * y * y * y * y * z * z * z + std::sqrt(119.91302490234375) * x * y * y * y * y * y * y * y * y * y * y * z - std::sqrt(852.71484375) * x * y * y * y * y * y * y * y * y * z * z * z) + e_1 * (-std::sqrt(1079.2172241210938) * x * x * x * x * x * x * x * x * x * z - std::sqrt(155407.2802734375) * x * x * x * x * x * x * x * y * y * z + std::sqrt(7674.43359375) * x * x * x * x * x * x * x * z * z * z + std::sqrt(2698043.0603027344) * x * x * x * x * x * y * y * y * y * z - std::sqrt(7674.43359375) * x * x * x * x * x * y * y * z * z * z - std::sqrt(324244.8193359375) * x * x * x * y * y * y * y * y * y * z - std::sqrt(191860.83984375) * x * x * x * y * y * y * y * z * z * z + std::sqrt(52881.643981933594) * x * y * y * y * y * y * y * y * y * z - std::sqrt(69069.90234375) * x * y * y * y * y * y * y * z * z * z) + e_2 * (-std::sqrt(431686.8896484375) * x * x * x * x * x * x * x * z + std::sqrt(431686.8896484375) * x * x * x * x * x * y * y * z + std::sqrt(767443.359375) * x * x * x * x * x * z * z * z + std::sqrt(10792172.241210938) * x * x * x * y * y * y * y * z - std::sqrt(3069773.4375) * x * x * x * y * y * z * z * z + std::sqrt(3885182.0068359375) * x * y * y * y * y * y * y * z - std::sqrt(6906990.234375) * x * y * y * y * y * z * z * z) + e_3 * (-std::sqrt(12279093.75) * x * x * x * x * x * z + std::sqrt(49116375.0) * x * x * x * y * y * z + std::sqrt(5457375.0) * x * x * x * z * z * z + std::sqrt(110511843.75) * x * y * y * y * y * z - std::sqrt(49116375.0) * x * y * y * z * z * z) + e_4 * (-std::sqrt(27627960.9375) * x * x * x * z + std::sqrt(248651648.4375) * x * y * y * z);
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

        pc_4[k] = e_0 * (std::sqrt(13.32366943359375) * x * x * x * x * x * x * x * x * x * x * y * y - std::sqrt(23.6865234375) * x * x * x * x * x * x * x * x * y * y * y * y - std::sqrt(3410.859375) * x * x * x * x * x * x * x * x * y * y * z * z - std::sqrt(290.159912109375) * x * x * x * x * x * x * y * y * y * y * y * y + std::sqrt(18570.234375) * x * x * x * x * x * x * y * y * y * y * z * z + std::sqrt(3410.859375) * x * x * x * x * x * x * y * y * z * z * z * z - std::sqrt(23.6865234375) * x * x * x * x * y * y * y * y * y * y * y * y + std::sqrt(18570.234375) * x * x * x * x * y * y * y * y * y * y * z * z - std::sqrt(37898.4375) * x * x * x * x * y * y * y * y * z * z * z * z + std::sqrt(13.32366943359375) * x * x * y * y * y * y * y * y * y * y * y * y - std::sqrt(3410.859375) * x * x * y * y * y * y * y * y * y * y * z * z + std::sqrt(3410.859375) * x * x * y * y * y * y * y * y * z * z * z * z) + e_1 * (std::sqrt(13.32366943359375) * x * x * x * x * x * x * x * x * x * x + std::sqrt(5875.738220214844) * x * x * x * x * x * x * x * x * y * y - std::sqrt(3410.859375) * x * x * x * x * x * x * x * x * z * z - std::sqrt(65285.980224609375) * x * x * x * x * x * x * y * y * y * y - std::sqrt(218295.0) * x * x * x * x * x * x * y * y * z * z + std::sqrt(3410.859375) * x * x * x * x * x * x * z * z * z * z - std::sqrt(65285.980224609375) * x * x * x * x * y * y * y * y * y * y + std::sqrt(8527148.4375) * x * x * x * x * y * y * y * y * z * z - std::sqrt(85271.484375) * x * x * x * x * y * y * z * z * z * z + std::sqrt(5875.738220214844) * x * x * y * y * y * y * y * y * y * y - std::sqrt(218295.0) * x * x * y * y * y * y * y * y * z * z - std::sqrt(85271.484375) * x * x * y * y * y * y * z * z * z * z + std::sqrt(13.32366943359375) * y * y * y * y * y * y * y * y * y * y - std::sqrt(3410.859375) * y * y * y * y * y * y * y * y * z * z + std::sqrt(3410.859375) * y * y * y * y * y * y * z * z * z * z) + e_2 * (std::sqrt(8327.293395996094) * x * x * x * x * x * x * x * x + std::sqrt(5329.4677734375) * x * x * x * x * x * x * y * y - std::sqrt(767443.359375) * x * x * x * x * x * x * z * z - std::sqrt(5629250.335693359) * x * x * x * x * y * y * y * y + std::sqrt(19186083.984375) * x * x * x * x * y * y * z * z + std::sqrt(85271.484375) * x * x * x * x * z * z * z * z + std::sqrt(5329.4677734375) * x * x * y * y * y * y * y * y + std::sqrt(19186083.984375) * x * x * y * y * y * y * z * z - std::sqrt(3069773.4375) * x * x * y * y * z * z * z * z + std::sqrt(8327.293395996094) * y * y * y * y * y * y * y * y - std::sqrt(767443.359375) * y * y * y * y * y * y * z * z + std::sqrt(85271.484375) * y * y * y * y * z * z * z * z) + e_3 * (std::sqrt(341085.9375) * x * x * x * x * x * x - std::sqrt(8527148.4375) * x * x * x * x * y * y - std::sqrt(5457375.0) * x * x * x * x * z * z - std::sqrt(8527148.4375) * x * x * y * y * y * y + std::sqrt(196465500.0) * x * x * y * y * z * z + std::sqrt(341085.9375) * y * y * y * y * y * y - std::sqrt(5457375.0) * y * y * y * y * z * z) + e_4 * (std::sqrt(767443.359375) * x * x * x * x - std::sqrt(27627960.9375) * x * x * y * y + std::sqrt(767443.359375) * y * y * y * y);
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

        pc_5[k] = e_0 * (std::sqrt(133.2366943359375) * x * x * x * x * x * x * x * x * x * y * y * z - std::sqrt(236.865234375) * x * x * x * x * x * x * x * y * y * y * y * z - std::sqrt(2131.787109375) * x * x * x * x * x * x * x * y * y * z * z * z - std::sqrt(2901.59912109375) * x * x * x * x * x * y * y * y * y * y * y * z + std::sqrt(11606.396484375) * x * x * x * x * x * y * y * y * y * z * z * z + std::sqrt(341.0859375) * x * x * x * x * x * y * y * z * z * z * z * z - std::sqrt(236.865234375) * x * x * x * y * y * y * y * y * y * y * y * z + std::sqrt(11606.396484375) * x * x * x * y * y * y * y * y * y * z * z * z - std::sqrt(3789.84375) * x * x * x * y * y * y * y * z * z * z * z * z + std::sqrt(133.2366943359375) * x * y * y * y * y * y * y * y * y * y * y * z - std::sqrt(2131.787109375) * x * y * y * y * y * y * y * y * y * z * z * z + std::sqrt(341.0859375) * x * y * y * y * y * y * y * z * z * z * z * z) + e_1 * (std::sqrt(133.2366943359375) * x * x * x * x * x * x * x * x * x * z + std::sqrt(34108.59375) * x * x * x * x * x * x * x * y * y * z - std::sqrt(2131.787109375) * x * x * x * x * x * x * x * z * z * z - std::sqrt(652859.8022460938) * x * x * x * x * x * y * y * y * y * z - std::sqrt(19186.083984375) * x * x * x * x * x * y * y * z * z * z + std::sqrt(341.0859375) * x * x * x * x * x * z * z * z * z * z - std::sqrt(417830.2734375) * x * x * x * y * y * y * y * y * y * z + std::sqrt(4316868.896484375) * x * x * x * y * y * y * y * z * z * z - std::sqrt(34108.59375) * x * x * x * y * y * z * z * z * z * z + std::sqrt(112052.05993652344) * x * y * y * y * y * y * y * y * y * z - std::sqrt(616086.474609375) * x * y * y * y * y * y * y * z * z * z + std::sqrt(8527.1484375) * x * y * y * y * y * z * z * z * z * z) + e_2 * (std::sqrt(53294.677734375) * x * x * x * x * x * x * x * z - std::sqrt(479652.099609375) * x * x * x * x * x * y * y * z - std::sqrt(213178.7109375) * x * x * x * x * x * z * z * z - std::sqrt(33309173.583984375) * x * x * x * y * y * y * y * z + std::sqrt(21317871.09375) * x * x * x * y * y * z * z * z + std::sqrt(6448656.005859375) * x * y * y * y * y * y * y * z - std::sqrt(5329467.7734375) * x * y * y * y * y * z * z * z) + e_3 * (std::sqrt(852714.84375) * x * x * x * x * x * z - std::sqrt(85271484.375) * x * x * x * y * y * z + std::sqrt(21317871.09375) * x * y * y * y * y * z);
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

        pc_6[k] = e_0 * (-std::sqrt(1.586151123046875) * x * x * x * x * x * x * x * x * x * x * x * y + std::sqrt(0.176239013671875) * x * x * x * x * x * x * x * x * x * y * y * y + std::sqrt(513.9129638671875) * x * x * x * x * x * x * x * x * x * y * z * z + std::sqrt(57.1014404296875) * x * x * x * x * x * x * x * y * y * y * y * y - std::sqrt(913.623046875) * x * x * x * x * x * x * x * y * y * y * z * z - std::sqrt(913.623046875) * x * x * x * x * x * x * x * y * z * z * z * z + std::sqrt(57.1014404296875) * x * x * x * x * x * y * y * y * y * y * y * y - std::sqrt(11191.88232421875) * x * x * x * x * x * y * y * y * y * y * z * z + std::sqrt(4974.169921875) * x * x * x * x * x * y * y * y * z * z * z * z + std::sqrt(16.2421875) * x * x * x * x * x * y * z * z * z * z * z * z + std::sqrt(0.176239013671875) * x * x * x * y * y * y * y * y * y * y * y * y - std::sqrt(913.623046875) * x * x * x * y * y * y * y * y * y * y * z * z + std::sqrt(4974.169921875) * x * x * x * y * y * y * y * y * z * z * z * z - std::sqrt(180.46875) * x * x * x * y * y * y * z * z * z * z * z * z - std::sqrt(1.586151123046875) * x * y * y * y * y * y * y * y * y * y * y * y + std::sqrt(513.9129638671875) * x * y * y * y * y * y * y * y * y * y * z * z - std::sqrt(913.623046875) * x * y * y * y * y * y * y * y * z * z * z * z + std::sqrt(16.2421875) * x * y * y * y * y * y * z * z * z * z * z * z) + e_1 * (-std::sqrt(2055.65185546875) * x * x * x * x * x * x * x * x * x * y + std::sqrt(3654.4921875) * x * x * x * x * x * x * x * y * y * y + std::sqrt(296013.8671875) * x * x * x * x * x * x * x * y * z * z + std::sqrt(44767.529296875) * x * x * x * x * x * y * y * y * y * y - std::sqrt(1611631.0546875) * x * x * x * x * x * y * y * y * z * z - std::sqrt(131561.71875) * x * x * x * x * x * y * z * z * z * z + std::sqrt(3654.4921875) * x * x * x * y * y * y * y * y * y * y - std::sqrt(1611631.0546875) * x * x * x * y * y * y * y * y * z * z + std::sqrt(1461796.875) * x * x * x * y * y * y * z * z * z * z - std::sqrt(2055.65185546875) * x * y * y * y * y * y * y * y * y * y + std::sqrt(296013.8671875) * x * y * y * y * y * y * y * y * z * z - std::sqrt(131561.71875) * x * y * y * y * y * y * z * z * z * z) + e_2 * (-std::sqrt(205565.185546875) * x * x * x * x * x * x * x * y + std::sqrt(1119188.232421875) * x * x * x * x * x * y * y * y + std::sqrt(7400346.6796875) * x * x * x * x * x * y * z * z + std::sqrt(1119188.232421875) * x * x * x * y * y * y * y * y - std::sqrt(82226074.21875) * x * x * x * y * y * y * z * z - std::sqrt(205565.185546875) * x * y * y * y * y * y * y * y + std::sqrt(7400346.6796875) * x * y * y * y * y * y * z * z) + e_3 * (-std::sqrt(1461796.875) * x * x * x * x * x * y + std::sqrt(16242187.5) * x * x * x * y * y * y - std::sqrt(1461796.875) * x * y * y * y * y * y);
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

        pc_7[k] = e_0 * (std::sqrt(133.2366943359375) * x * x * x * x * x * x * x * x * x * x * y * z - std::sqrt(236.865234375) * x * x * x * x * x * x * x * x * y * y * y * z - std::sqrt(2131.787109375) * x * x * x * x * x * x * x * x * y * z * z * z - std::sqrt(2901.59912109375) * x * x * x * x * x * x * y * y * y * y * y * z + std::sqrt(11606.396484375) * x * x * x * x * x * x * y * y * y * z * z * z + std::sqrt(341.0859375) * x * x * x * x * x * x * y * z * z * z * z * z - std::sqrt(236.865234375) * x * x * x * x * y * y * y * y * y * y * y * z + std::sqrt(11606.396484375) * x * x * x * x * y * y * y * y * y * z * z * z - std::sqrt(3789.84375) * x * x * x * x * y * y * y * z * z * z * z * z + std::sqrt(133.2366943359375) * x * x * y * y * y * y * y * y * y * y * y * z - std::sqrt(2131.787109375) * x * x * y * y * y * y * y * y * y * z * z * z + std::sqrt(341.0859375) * x * x * y * y * y * y * y * z * z * z * z * z) + e_1 * (std::sqrt(112052.05993652344) * x * x * x * x * x * x * x * x * y * z - std::sqrt(417830.2734375) * x * x * x * x * x * x * y * y * y * z - std::sqrt(616086.474609375) * x * x * x * x * x * x * y * z * z * z - std::sqrt(652859.8022460938) * x * x * x * x * y * y * y * y * y * z + std::sqrt(4316868.896484375) * x * x * x * x * y * y * y * z * z * z + std::sqrt(8527.1484375) * x * x * x * x * y * z * z * z * z * z + std::sqrt(34108.59375) * x * x * y * y * y * y * y * y * y * z - std::sqrt(19186.083984375) * x * x * y * y * y * y * y * z * z * z - std::sqrt(34108.59375) * x * x * y * y * y * z * z * z * z * z + std::sqrt(133.2366943359375) * y * y * y * y * y * y * y * y * y * z - std::sqrt(2131.787109375) * y * y * y * y * y * y * y * z * z * z + std::sqrt(341.0859375) * y * y * y * y * y * z * z * z * z * z) + e_2 * (std::sqrt(6448656.005859375) * x * x * x * x * x * x * y * z - std::sqrt(33309173.583984375) * x * x * x * x * y * y * y * z - std::sqrt(5329467.7734375) * x * x * x * x * y * z * z * z - std::sqrt(479652.099609375) * x * x * y * y * y * y * y * z + std::sqrt(21317871.09375) * x * x * y * y * y * z * z * z + std::sqrt(53294.677734375) * y * y * y * y * y * y * y * z - std::sqrt(213178.7109375) * y * y * y * y * y * z * z * z) + e_3 * (std::sqrt(21317871.09375) * x * x * x * x * y * z - std::sqrt(85271484.375) * x * x * y * y * y * z + std::sqrt(852714.84375) * y * y * y * y * y * z);
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

        pc_8[k] = e_0 * (std::sqrt(3.3309173583984375) * x * x * x * x * x * x * x * x * x * x * x * y - std::sqrt(18.134994506835938) * x * x * x * x * x * x * x * x * x * y * y * y - std::sqrt(852.71484375) * x * x * x * x * x * x * x * x * x * y * z * z - std::sqrt(37.01019287109375) * x * x * x * x * x * x * x * y * y * y * y * y + std::sqrt(9474.609375) * x * x * x * x * x * x * x * y * y * y * z * z + std::sqrt(852.71484375) * x * x * x * x * x * x * x * y * z * z * z * z + std::sqrt(37.01019287109375) * x * x * x * x * x * y * y * y * y * y * y * y - std::sqrt(16012.08984375) * x * x * x * x * x * y * y * y * z * z * z * z + std::sqrt(18.134994506835938) * x * x * x * y * y * y * y * y * y * y * y * y - std::sqrt(9474.609375) * x * x * x * y * y * y * y * y * y * y * z * z + std::sqrt(16012.08984375) * x * x * x * y * y * y * y * y * z * z * z * z - std::sqrt(3.3309173583984375) * x * y * y * y * y * y * y * y * y * y * y * y + std::sqrt(852.71484375) * x * y * y * y * y * y * y * y * y * y * z * z - std::sqrt(852.71484375) * x * y * y * y * y * y * y * y * z * z * z * z) + e_1 * (std::sqrt(3410.859375) * x * x * x * x * x * x * x * x * x * y - std::sqrt(13643.4375) * x * x * x * x * x * x * x * y * y * y - std::sqrt(341085.9375) * x * x * x * x * x * x * x * y * z * z + std::sqrt(1650855.9375) * x * x * x * x * x * y * y * y * z * z + std::sqrt(54573.75) * x * x * x * x * x * y * z * z * z * z + std::sqrt(13643.4375) * x * x * x * y * y * y * y * y * y * y - std::sqrt(1650855.9375) * x * x * x * y * y * y * y * y * z * z - std::sqrt(3410.859375) * x * y * y * y * y * y * y * y * y * y + std::sqrt(341085.9375) * x * y * y * y * y * y * y * y * z * z - std::sqrt(54573.75) * x * y * y * y * y * y * z * z * z * z) + e_2 * (std::sqrt(341085.9375) * x * x * x * x * x * x * x * y - std::sqrt(341085.9375) * x * x * x * x * x * y * y * y - std::sqrt(12279093.75) * x * x * x * x * x * y * z * z + std::sqrt(341085.9375) * x * x * x * y * y * y * y * y + std::sqrt(1364343.75) * x * x * x * y * z * z * z * z - std::sqrt(341085.9375) * x * y * y * y * y * y * y * y + std::sqrt(12279093.75) * x * y * y * y * y * y * z * z - std::sqrt(1364343.75) * x * y * y * y * z * z * z * z) + e_3 * (std::sqrt(5457375.0) * x * x * x * x * x * y - std::sqrt(87318000.0) * x * x * x * y * z * z - std::sqrt(5457375.0) * x * y * y * y * y * y + std::sqrt(87318000.0) * x * y * y * y * z * z) + e_4 * (std::sqrt(12279093.75) * x * x * x * y - std::sqrt(12279093.75) * x * y * y * y);
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

        pc_9[k] = e_0 * (-std::sqrt(119.91302490234375) * x * x * x * x * x * x * x * x * x * x * y * z + std::sqrt(3410.859375) * x * x * x * x * x * x * x * x * y * y * y * z + std::sqrt(852.71484375) * x * x * x * x * x * x * x * x * y * z * z * z - std::sqrt(2611.439208984375) * x * x * x * x * x * x * y * y * y * y * y * z - std::sqrt(34203.33984375) * x * x * x * x * x * x * y * y * y * z * z * z - std::sqrt(7674.43359375) * x * x * x * x * y * y * y * y * y * y * y * z + std::sqrt(103178.49609375) * x * x * x * x * y * y * y * y * y * z * z * z + std::sqrt(1079.2172241210938) * x * x * y * y * y * y * y * y * y * y * y * z - std::sqrt(7674.43359375) * x * x * y * y * y * y * y * y * y * z * z * z) + e_1 * (-std::sqrt(52881.643981933594) * x * x * x * x * x * x * x * x * y * z + std::sqrt(324244.8193359375) * x * x * x * x * x * x * y * y * y * z + std::sqrt(69069.90234375) * x * x * x * x * x * x * y * z * z * z - std::sqrt(2698043.0603027344) * x * x * x * x * y * y * y * y * y * z + std::sqrt(191860.83984375) * x * x * x * x * y * y * y * z * z * z + std::sqrt(155407.2802734375) * x * x * y * y * y * y * y * y * y * z + std::sqrt(7674.43359375) * x * x * y * y * y * y * y * z * z * z + std::sqrt(1079.2172241210938) * y * y * y * y * y * y * y * y * y * z - std::sqrt(7674.43359375) * y * y * y * y * y * y * y * z * z * z) + e_2 * (-std::sqrt(3885182.0068359375) * x * x * x * x * x * x * y * z - std::sqrt(10792172.241210938) * x * x * x * x * y * y * y * z + std::sqrt(6906990.234375) * x * x * x * x * y * z * z * z - std::sqrt(431686.8896484375) * x * x * y * y * y * y * y * z + std::sqrt(3069773.4375) * x * x * y * y * y * z * z * z + std::sqrt(431686.8896484375) * y * y * y * y * y * y * y * z - std::sqrt(767443.359375) * y * y * y * y * y * z * z * z) + e_3 * (-std::sqrt(110511843.75) * x * x * x * x * y * z - std::sqrt(49116375.0) * x * x * y * y * y * z + std::sqrt(49116375.0) * x * x * y * z * z * z + std::sqrt(12279093.75) * y * y * y * y * y * z - std::sqrt(5457375.0) * y * y * y * z * z * z) + e_4 * (-std::sqrt(248651648.4375) * x * x * y * z + std::sqrt(27627960.9375) * y * y * y * z);
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

        pc_10[k] = e_0 * (-std::sqrt(3.997100830078125) * x * x * x * x * x * x * x * x * x * x * x * y + std::sqrt(277.5764465332031) * x * x * x * x * x * x * x * x * x * y * y * y + std::sqrt(399.7100830078125) * x * x * x * x * x * x * x * x * x * y * z * z - std::sqrt(641.3126220703125) * x * x * x * x * x * x * x * y * y * y * y * y - std::sqrt(34819.189453125) * x * x * x * x * x * x * x * y * y * y * z * z - std::sqrt(641.3126220703125) * x * x * x * x * x * y * y * y * y * y * y * y + std::sqrt(193459.68017578125) * x * x * x * x * x * y * y * y * y * y * z * z + std::sqrt(277.5764465332031) * x * x * x * y * y * y * y * y * y * y * y * y - std::sqrt(34819.189453125) * x * x * x * y * y * y * y * y * y * y * z * z - std::sqrt(3.997100830078125) * x * y * y * y * y * y * y * y * y * y * y * y + std::sqrt(399.7100830078125) * x * y * y * y * y * y * y * y * y * y * z * z) + e_1 * (-std::sqrt(1598.84033203125) * x * x * x * x * x * x * x * x * x * y + std::sqrt(25581.4453125) * x * x * x * x * x * x * x * y * y * y + std::sqrt(25581.4453125) * x * x * x * x * x * x * x * y * z * z - std::sqrt(389093.783203125) * x * x * x * x * x * y * y * y * y * y + std::sqrt(230233.0078125) * x * x * x * x * x * y * y * y * z * z + std::sqrt(25581.4453125) * x * x * x * y * y * y * y * y * y * y + std::sqrt(230233.0078125) * x * x * x * y * y * y * y * y * z * z - std::sqrt(1598.84033203125) * x * y * y * y * y * y * y * y * y * y + std::sqrt(25581.4453125) * x * y * y * y * y * y * y * y * z * z) + e_2 * (-std::sqrt(159884.033203125) * x * x * x * x * x * x * x * y - std::sqrt(1438956.298828125) * x * x * x * x * x * y * y * y + std::sqrt(5755825.1953125) * x * x * x * x * x * y * z * z - std::sqrt(1438956.298828125) * x * x * x * y * y * y * y * y + std::sqrt(23023300.78125) * x * x * x * y * y * y * z * z - std::sqrt(159884.033203125) * x * y * y * y * y * y * y * y + std::sqrt(5755825.1953125) * x * y * y * y * y * y * z * z) + e_3 * (-std::sqrt(10232578.125) * x * x * x * x * x * y - std::sqrt(40930312.5) * x * x * x * y * y * y + std::sqrt(163721250.0) * x * x * x * y * z * z - std::sqrt(10232578.125) * x * y * y * y * y * y + std::sqrt(163721250.0) * x * y * y * y * z * z) + e_4 * (-std::sqrt(92093203.125) * x * x * x * y - std::sqrt(92093203.125) * x * y * y * y + std::sqrt(368372812.5) * x * y * z * z) + e_5 * (-std::sqrt(58939650.0) * x * y);
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

        pc_11[k] = e_0 * (std::sqrt(87.93621826171875) * x * x * x * x * x * x * x * x * x * x * y * z - std::sqrt(15633.10546875) * x * x * x * x * x * x * x * x * y * y * y * z + std::sqrt(136047.10034179688) * x * x * x * x * x * x * y * y * y * y * y * z - std::sqrt(62532.421875) * x * x * x * x * y * y * y * y * y * y * y * z + std::sqrt(2198.4054565429688) * x * x * y * y * y * y * y * y * y * y * y * z) + e_1 * (std::sqrt(2198.4054565429688) * x * x * x * x * x * x * x * x * y * z + std::sqrt(35174.4873046875) * x * x * x * x * x * x * y * y * y * z + std::sqrt(79142.59643554688) * x * x * x * x * y * y * y * y * y * z + std::sqrt(35174.4873046875) * x * x * y * y * y * y * y * y * y * z + std::sqrt(2198.4054565429688) * y * y * y * y * y * y * y * y * y * z) + e_2 * (std::sqrt(879362.1826171875) * x * x * x * x * x * x * y * z + std::sqrt(7914259.6435546875) * x * x * x * x * y * y * y * z + std::sqrt(7914259.6435546875) * x * x * y * y * y * y * y * z + std::sqrt(879362.1826171875) * y * y * y * y * y * y * y * z) + e_3 * (std::sqrt(56279179.6875) * x * x * x * x * y * z + std::sqrt(225116718.75) * x * x * y * y * y * z + std::sqrt(56279179.6875) * y * y * y * y * y * z) + e_4 * (std::sqrt(506512617.1875) * x * x * y * z + std::sqrt(506512617.1875) * y * y * y * z) + e_5 * (std::sqrt(324168075.0) * y * z);

        pc_12[k] = e_0 * (2.70703125 * x * x * x * x * x * x * x * x * x * x * x * y - 49.62890625 * x * x * x * x * x * x * x * x * x * y * y * y + 178.6640625 * x * x * x * x * x * x * x * y * y * y * y * y - 178.6640625 * x * x * x * x * x * y * y * y * y * y * y * y + 49.62890625 * x * x * x * y * y * y * y * y * y * y * y * y - 2.70703125 * x * y * y * y * y * y * y * y * y * y * y * y);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ab_x, ab_y, ab_z : simd::cache_line_size())
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
        const auto e_6 = pe_6[k];

        pc_13[k] = e_0 * (135.3515625 * x * x * x * x * x * x * x * x * y * y * z * z - 541.40625 * x * x * x * x * x * x * y * y * y * y * z * z + 595.546875 * x * x * x * x * y * y * y * y * y * y * z * z - 108.28125 * x * x * y * y * y * y * y * y * y * y * z * z + 5.4140625 * y * y * y * y * y * y * y * y * y * y * z * z) + e_1 * (135.3515625 * x * x * x * x * x * x * x * x * y * y + 135.3515625 * x * x * x * x * x * x * x * x * z * z - 541.40625 * x * x * x * x * x * x * y * y * y * y + 541.40625 * x * x * x * x * x * x * y * y * z * z + 595.546875 * x * x * x * x * y * y * y * y * y * y + 812.109375 * x * x * x * x * y * y * y * y * z * z - 108.28125 * x * x * y * y * y * y * y * y * y * y + 541.40625 * x * x * y * y * y * y * y * y * z * z + 5.4140625 * y * y * y * y * y * y * y * y * y * y + 135.3515625 * y * y * y * y * y * y * y * y * z * z) + e_2 * (135.3515625 * x * x * x * x * x * x * x * x + 541.40625 * x * x * x * x * x * x * y * y + 2165.625 * x * x * x * x * x * x * z * z + 812.109375 * x * x * x * x * y * y * y * y + 6496.875 * x * x * x * x * y * y * z * z + 541.40625 * x * x * y * y * y * y * y * y + 6496.875 * x * x * y * y * y * y * z * z + 135.3515625 * y * y * y * y * y * y * y * y + 2165.625 * y * y * y * y * y * y * z * z) + e_3 * (2165.625 * x * x * x * x * x * x + 6496.875 * x * x * x * x * y * y + 12993.75 * x * x * x * x * z * z + 6496.875 * x * x * y * y * y * y + 25987.5 * x * x * y * y * z * z + 2165.625 * y * y * y * y * y * y + 12993.75 * y * y * y * y * z * z) + e_4 * (12993.75 * x * x * x * x + 25987.5 * x * x * y * y + 25987.5 * x * x * z * z + 12993.75 * y * y * y * y + 25987.5 * y * y * z * z) + e_5 * (25987.5 * x * x + 25987.5 * y * y + 10395.0 * z * z) + e_6 * (10395.0);
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

        pc_14[k] = e_0 * (-std::sqrt(532.94677734375) * x * x * x * x * x * x * x * x * x * y * y * z + std::sqrt(2131.787109375) * x * x * x * x * x * x * x * y * y * y * y * z + std::sqrt(53294.677734375) * x * x * x * x * x * x * x * y * y * z * z * z + std::sqrt(341.0859375) * x * x * x * x * x * y * y * y * y * y * y * z - std::sqrt(479652.099609375) * x * x * x * x * x * y * y * y * y * z * z * z - std::sqrt(2131.787109375) * x * x * x * y * y * y * y * y * y * y * y * z + std::sqrt(257946.240234375) * x * x * x * y * y * y * y * y * y * z * z * z + std::sqrt(21.31787109375) * x * y * y * y * y * y * y * y * y * y * y * z - std::sqrt(2131.787109375) * x * y * y * y * y * y * y * y * y * z * z * z) + e_1 * (-std::sqrt(532.94677734375) * x * x * x * x * x * x * x * x * x * z + std::sqrt(19186.083984375) * x * x * x * x * x * x * x * y * y * z + std::sqrt(53294.677734375) * x * x * x * x * x * x * x * z * z * z - std::sqrt(690699.0234375) * x * x * x * x * x * y * y * y * y * z + std::sqrt(479652.099609375) * x * x * x * x * x * y * y * z * z * z + std::sqrt(172674.755859375) * x * x * x * y * y * y * y * y * y * z + std::sqrt(479652.099609375) * x * x * x * y * y * y * y * z * z * z - std::sqrt(4796.52099609375) * x * y * y * y * y * y * y * y * y * z + std::sqrt(53294.677734375) * x * y * y * y * y * y * y * z * z * z) + e_2 * (std::sqrt(7674433.59375) * x * x * x * x * x * z * z * z + std::sqrt(30697734.375) * x * x * x * y * y * z * z * z + std::sqrt(7674433.59375) * x * y * y * y * y * z * z * z) + e_3 * (std::sqrt(7674433.59375) * x * x * x * x * x * z + std::sqrt(30697734.375) * x * x * x * y * y * z + std::sqrt(122790937.5) * x * x * x * z * z * z + std::sqrt(7674433.59375) * x * y * y * y * y * z + std::sqrt(122790937.5) * x * y * y * z * z * z) + e_4 * (std::sqrt(276279609.375) * x * x * x * z + std::sqrt(276279609.375) * x * y * y * z + std::sqrt(122790937.5) * x * z * z * z) + e_5 * (std::sqrt(397842637.5) * x * z);
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

        pc_15[k] = e_0 * (-std::sqrt(8993.476867675781) * x * x * x * x * x * x * x * x * y * y * z * z + std::sqrt(15988.4033203125) * x * x * x * x * x * x * y * y * y * y * z * z + std::sqrt(63953.61328125) * x * x * x * x * x * x * y * y * z * z * z * z + std::sqrt(19345.968017578125) * x * x * x * x * y * y * y * y * y * y * z * z - std::sqrt(348191.89453125) * x * x * x * x * y * y * y * y * z * z * z * z - std::sqrt(5755.8251953125) * x * x * y * y * y * y * y * y * y * y * z * z + std::sqrt(48036.26953125) * x * x * y * y * y * y * y * y * z * z * z * z + std::sqrt(39.97100830078125) * y * y * y * y * y * y * y * y * y * y * z * z - std::sqrt(284.23828125) * y * y * y * y * y * y * y * y * z * z * z * z) + e_1 * (-std::sqrt(8993.476867675781) * x * x * x * x * x * x * x * x * y * y - std::sqrt(8993.476867675781) * x * x * x * x * x * x * x * x * z * z + std::sqrt(15988.4033203125) * x * x * x * x * x * x * y * y * y * y - std::sqrt(143895.6298828125) * x * x * x * x * x * x * y * y * z * z + std::sqrt(63953.61328125) * x * x * x * x * x * x * z * z * z * z + std::sqrt(19345.968017578125) * x * x * x * x * y * y * y * y * y * y + std::sqrt(195857.94067382812) * x * x * x * x * y * y * y * y * z * z + std::sqrt(63953.61328125) * x * x * x * x * y * y * z * z * z * z - std::sqrt(5755.8251953125) * x * x * y * y * y * y * y * y * y * y + std::sqrt(639.5361328125) * x * x * y * y * y * y * y * y * z * z - std::sqrt(63953.61328125) * x * x * y * y * y * y * z * z * z * z + std::sqrt(39.97100830078125) * y * y * y * y * y * y * y * y * y * y + std::sqrt(11551.621398925781) * y * y * y * y * y * y * y * y * z * z - std::sqrt(63953.61328125) * y * y * y * y * y * y * z * z * z * z) + e_2 * (-std::sqrt(8993.476867675781) * x * x * x * x * x * x * x * x - std::sqrt(1295060.6689453125) * x * x * x * x * x * x * y * y - std::sqrt(575582.51953125) * x * x * x * x * x * x * z * z + std::sqrt(4896448.516845703) * x * x * x * x * y * y * y * y - std::sqrt(575582.51953125) * x * x * x * x * y * y * z * z + std::sqrt(4093031.25) * x * x * x * x * z * z * z * z - std::sqrt(399710.0830078125) * x * x * y * y * y * y * y * y + std::sqrt(575582.51953125) * x * x * y * y * y * y * z * z + std::sqrt(24981.88018798828) * y * y * y * y * y * y * y * y + std::sqrt(575582.51953125) * y * y * y * y * y * y * z * z - std::sqrt(4093031.25) * y * y * y * y * z * z * z * z) + e_3 * (-std::sqrt(2302330.078125) * x * x * x * x * x * x - std::sqrt(2302330.078125) * x * x * x * x * y * y + std::sqrt(2302330.078125) * x * x * y * y * y * y + std::sqrt(16372125.0) * x * x * z * z * z * z + std::sqrt(2302330.078125) * y * y * y * y * y * y - std::sqrt(16372125.0) * y * y * z * z * z * z) + e_4 * (-std::sqrt(36837281.25) * x * x * x * x + std::sqrt(36837281.25) * x * x * z * z + std::sqrt(36837281.25) * y * y * y * y - std::sqrt(36837281.25) * y * y * z * z) + e_5 * (-std::sqrt(36837281.25) * x * x + std::sqrt(36837281.25) * y * y);
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

        pc_16[k] = e_0 * (std::sqrt(111.03057861328125) * x * x * x * x * x * x * x * x * x * y * y * z - std::sqrt(28423.828125) * x * x * x * x * x * x * x * y * y * z * z * z - std::sqrt(870.479736328125) * x * x * x * x * x * y * y * y * y * y * y * z + std::sqrt(28423.828125) * x * x * x * x * x * y * y * y * y * z * z * z + std::sqrt(28423.828125) * x * x * x * x * x * y * y * z * z * z * z * z - std::sqrt(284.23828125) * x * x * x * y * y * y * y * y * y * y * y * z + std::sqrt(92093.203125) * x * x * x * y * y * y * y * y * y * z * z * z - std::sqrt(113695.3125) * x * x * x * y * y * y * y * z * z * z * z * z + std::sqrt(4.44122314453125) * x * y * y * y * y * y * y * y * y * y * y * z - std::sqrt(1136.953125) * x * y * y * y * y * y * y * y * y * z * z * z + std::sqrt(1136.953125) * x * y * y * y * y * y * y * z * z * z * z * z) + e_1 * (std::sqrt(111.03057861328125) * x * x * x * x * x * x * x * x * x * z - std::sqrt(15988.4033203125) * x * x * x * x * x * x * x * y * y * z - std::sqrt(28423.828125) * x * x * x * x * x * x * x * z * z * z + std::sqrt(3997.100830078125) * x * x * x * x * x * y * y * y * y * z - std::sqrt(710595.703125) * x * x * x * x * x * y * y * z * z * z + std::sqrt(28423.828125) * x * x * x * x * x * z * z * z * z * z + std::sqrt(20536.2158203125) * x * x * x * y * y * y * y * y * y * z + std::sqrt(8214486.328125) * x * x * x * y * y * y * y * z * z * z - std::sqrt(113695.3125) * x * x * x * y * y * z * z * z * z * z - std::sqrt(3237.6516723632812) * x * y * y * y * y * y * y * y * y * z + std::sqrt(92093.203125) * x * y * y * y * y * y * y * z * z * z - std::sqrt(255814.453125) * x * y * y * y * y * z * z * z * z * z) + e_2 * (-std::sqrt(15988.4033203125) * x * x * x * x * x * x * x * z - std::sqrt(5771813.5986328125) * x * x * x * x * x * y * y * z - std::sqrt(1819125.0) * x * x * x * x * x * z * z * z + std::sqrt(32376516.723632812) * x * x * x * y * y * y * y * z + std::sqrt(7276500.0) * x * x * x * y * y * z * z * z + std::sqrt(454781.25) * x * x * x * z * z * z * z * z - std::sqrt(15988.4033203125) * x * y * y * y * y * y * y * z + std::sqrt(16372125.0) * x * y * y * y * y * z * z * z - std::sqrt(4093031.25) * x * y * y * z * z * z * z * z) + e_3 * (-std::sqrt(9209320.3125) * x * x * x * x * x * z + std::sqrt(36837281.25) * x * x * x * y * y * z - std::sqrt(1819125.0) * x * x * x * z * z * z + std::sqrt(82883882.8125) * x * y * y * y * y * z + std::sqrt(16372125.0) * x * y * y * z * z * z) + e_4 * (-std::sqrt(50139632.8125) * x * x * x * z + std::sqrt(451256695.3125) * x * y * y * z);
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

        pc_17[k] = e_0 * (std::sqrt(1110.3057861328125) * x * x * x * x * x * x * x * x * y * y * z * z - std::sqrt(17764.892578125) * x * x * x * x * x * x * y * y * z * z * z * z - std::sqrt(8704.79736328125) * x * x * x * x * y * y * y * y * y * y * z * z + std::sqrt(17764.892578125) * x * x * x * x * y * y * y * y * z * z * z * z + std::sqrt(2842.3828125) * x * x * x * x * y * y * z * z * z * z * z * z - std::sqrt(2842.3828125) * x * x * y * y * y * y * y * y * y * y * z * z + std::sqrt(57558.251953125) * x * x * y * y * y * y * y * y * z * z * z * z - std::sqrt(11369.53125) * x * x * y * y * y * y * z * z * z * z * z * z + std::sqrt(44.4122314453125) * y * y * y * y * y * y * y * y * y * y * z * z - std::sqrt(710.595703125) * y * y * y * y * y * y * y * y * z * z * z * z + std::sqrt(113.6953125) * y * y * y * y * y * y * z * z * z * z * z * z) + e_1 * (std::sqrt(1110.3057861328125) * x * x * x * x * x * x * x * x * y * y + std::sqrt(1110.3057861328125) * x * x * x * x * x * x * x * x * z * z + std::sqrt(17764.892578125) * x * x * x * x * x * x * y * y * z * z - std::sqrt(17764.892578125) * x * x * x * x * x * x * z * z * z * z - std::sqrt(8704.79736328125) * x * x * x * x * y * y * y * y * y * y - std::sqrt(359739.07470703125) * x * x * x * x * y * y * y * y * z * z - std::sqrt(159884.033203125) * x * x * x * x * y * y * z * z * z * z + std::sqrt(2842.3828125) * x * x * x * x * z * z * z * z * z * z - std::sqrt(2842.3828125) * x * x * y * y * y * y * y * y * y * y - std::sqrt(375905.126953125) * x * x * y * y * y * y * y * y * z * z + std::sqrt(7834317.626953125) * x * x * y * y * y * y * z * z * z * z - std::sqrt(102325.78125) * x * x * y * y * z * z * z * z * z * z + std::sqrt(44.4122314453125) * y * y * y * y * y * y * y * y * y * y + std::sqrt(7505.6671142578125) * y * y * y * y * y * y * y * y * z * z - std::sqrt(120090.673828125) * y * y * y * y * y * y * z * z * z * z + std::sqrt(2842.3828125) * y * y * y * y * z * z * z * z * z * z) + e_2 * (std::sqrt(1110.3057861328125) * x * x * x * x * x * x * x * x + std::sqrt(284238.28125) * x * x * x * x * x * x * y * y + std::sqrt(17764.892578125) * x * x * x * x * x * x * z * z - std::sqrt(999275.2075195312) * x * x * x * x * y * y * y * y - std::sqrt(3997100.830078125) * x * x * x * x * y * y * z * z - std::sqrt(639536.1328125) * x * x * x * x * z * z * z * z - std::sqrt(1776489.2578125) * x * x * y * y * y * y * y * y + std::sqrt(3997100.830078125) * x * x * y * y * y * y * z * z + std::sqrt(23023300.78125) * x * x * y * y * z * z * z * z + std::sqrt(27757.644653320312) * y * y * y * y * y * y * y * y - std::sqrt(17764.892578125) * y * y * y * y * y * y * z * z - std::sqrt(639536.1328125) * y * y * y * y * z * z * z * z) + e_3 * (std::sqrt(284238.28125) * x * x * x * x * x * x - std::sqrt(2558144.53125) * x * x * x * x * z * z - std::sqrt(63953613.28125) * x * x * y * y * y * y + std::sqrt(92093203.125) * x * x * y * y * z * z + std::sqrt(1136953.125) * y * y * y * y * y * y - std::sqrt(2558144.53125) * y * y * y * y * z * z) + e_4 * (std::sqrt(2558144.53125) * x * x * x * x - std::sqrt(92093203.125) * x * x * y * y + std::sqrt(2558144.53125) * y * y * y * y);
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

        pc_18[k] = e_0 * (-std::sqrt(13.217926025390625) * x * x * x * x * x * x * x * x * x * x * y * z - std::sqrt(13.217926025390625) * x * x * x * x * x * x * x * x * y * y * y * z + std::sqrt(4282.6080322265625) * x * x * x * x * x * x * x * x * y * z * z * z + std::sqrt(103.6285400390625) * x * x * x * x * x * x * y * y * y * y * y * z - std::sqrt(7613.525390625) * x * x * x * x * x * x * y * z * z * z * z * z + std::sqrt(255.8990478515625) * x * x * x * x * y * y * y * y * y * y * y * z - std::sqrt(33575.64697265625) * x * x * x * x * y * y * y * y * y * z * z * z + std::sqrt(7613.525390625) * x * x * x * x * y * y * y * z * z * z * z * z + std::sqrt(135.3515625) * x * x * x * x * y * z * z * z * z * z * z * z + std::sqrt(25.907135009765625) * x * x * y * y * y * y * y * y * y * y * y * z - std::sqrt(10963.4765625) * x * x * y * y * y * y * y * y * y * z * z * z + std::sqrt(24667.822265625) * x * x * y * y * y * y * y * z * z * z * z * z - std::sqrt(541.40625) * x * x * y * y * y * z * z * z * z * z * z * z - std::sqrt(0.528717041015625) * y * y * y * y * y * y * y * y * y * y * y * z + std::sqrt(171.3043212890625) * y * y * y * y * y * y * y * y * y * z * z * z - std::sqrt(304.541015625) * y * y * y * y * y * y * y * z * z * z * z * z + std::sqrt(5.4140625) * y * y * y * y * y * z * z * z * z * z * z * z) + e_1 * (std::sqrt(475.8453369140625) * x * x * x * x * x * x * x * x * y * z + std::sqrt(921236.572265625) * x * x * x * x * x * x * y * z * z * z - std::sqrt(3730.62744140625) * x * x * x * x * y * y * y * y * y * z - std::sqrt(921236.572265625) * x * x * x * x * y * y * y * z * z * z - std::sqrt(644408.7890625) * x * x * x * x * y * z * z * z * z * z - std::sqrt(1218.1640625) * x * x * y * y * y * y * y * y * y * z - std::sqrt(2984806.494140625) * x * x * y * y * y * y * y * z * z * z + std::sqrt(2577635.15625) * x * x * y * y * y * z * z * z * z * z + std::sqrt(19.0338134765625) * y * y * y * y * y * y * y * y * y * z + std::sqrt(36849.462890625) * y * y * y * y * y * y * y * z * z * z - std::sqrt(25776.3515625) * y * y * y * y * y * z * z * z * z * z) + e_2 * (std::sqrt(3045410.15625) * x * x * x * x * x * x * y * z - std::sqrt(3045410.15625) * x * x * x * x * y * y * y * z + std::sqrt(3045410.15625) * x * x * x * x * y * z * z * z - std::sqrt(9867128.90625) * x * x * y * y * y * y * y * z - std::sqrt(12181640.625) * x * x * y * y * y * z * z * z + std::sqrt(121816.40625) * y * y * y * y * y * y * y * z + std::sqrt(121816.40625) * y * y * y * y * y * z * z * z) + e_3 * (std::sqrt(76135253.90625) * x * x * x * x * y * z - std::sqrt(304541015.625) * x * x * y * y * y * z + std::sqrt(3045410.15625) * y * y * y * y * y * z);
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

        pc_19[k] = e_0 * (std::sqrt(1110.3057861328125) * x * x * x * x * x * x * x * x * x * y * z * z - std::sqrt(17764.892578125) * x * x * x * x * x * x * x * y * z * z * z * z - std::sqrt(8704.79736328125) * x * x * x * x * x * y * y * y * y * y * z * z + std::sqrt(17764.892578125) * x * x * x * x * x * y * y * y * z * z * z * z + std::sqrt(2842.3828125) * x * x * x * x * x * y * z * z * z * z * z * z - std::sqrt(2842.3828125) * x * x * x * y * y * y * y * y * y * y * z * z + std::sqrt(57558.251953125) * x * x * x * y * y * y * y * y * z * z * z * z - std::sqrt(11369.53125) * x * x * x * y * y * y * z * z * z * z * z * z + std::sqrt(44.4122314453125) * x * y * y * y * y * y * y * y * y * y * z * z - std::sqrt(710.595703125) * x * y * y * y * y * y * y * y * z * z * z * z + std::sqrt(113.6953125) * x * y * y * y * y * y * z * z * z * z * z * z) + e_1 * (std::sqrt(1110.3057861328125) * x * x * x * x * x * x * x * x * x * y + std::sqrt(159884.033203125) * x * x * x * x * x * x * x * y * z * z - std::sqrt(8704.79736328125) * x * x * x * x * x * y * y * y * y * y - std::sqrt(17764.892578125) * x * x * x * x * x * y * y * y * z * z - std::sqrt(2558144.53125) * x * x * x * x * x * y * z * z * z * z - std::sqrt(2842.3828125) * x * x * x * y * y * y * y * y * y * y - std::sqrt(375905.126953125) * x * x * x * y * y * y * y * y * z * z + std::sqrt(4547812.5) * x * x * x * y * y * y * z * z * z * z + std::sqrt(45478.125) * x * x * x * y * z * z * z * z * z * z + std::sqrt(44.4122314453125) * x * y * y * y * y * y * y * y * y * y - std::sqrt(6395.361328125) * x * y * y * y * y * y * y * y * z * z + std::sqrt(102325.78125) * x * y * y * y * y * y * z * z * z * z - std::sqrt(45478.125) * x * y * y * y * z * z * z * z * z * z) + e_2 * (std::sqrt(639536.1328125) * x * x * x * x * x * x * x * y - std::sqrt(284238.28125) * x * x * x * x * x * y * y * y - std::sqrt(639536.1328125) * x * x * x * x * x * y * z * z - std::sqrt(1776489.2578125) * x * x * x * y * y * y * y * y + std::sqrt(7105957.03125) * x * x * x * y * y * y * z * z - std::sqrt(10232578.125) * x * x * x * y * z * z * z * z - std::sqrt(639536.1328125) * x * y * y * y * y * y * z * z + std::sqrt(10232578.125) * x * y * y * y * z * z * z * z) + e_3 * (std::sqrt(23023300.78125) * x * x * x * x * x * y - std::sqrt(28423828.125) * x * x * x * y * y * y - std::sqrt(40930312.5) * x * x * x * y * z * z - std::sqrt(2558144.53125) * x * y * y * y * y * y + std::sqrt(40930312.5) * x * y * y * y * z * z) + e_4 * (std::sqrt(40930312.5) * x * x * x * y - std::sqrt(40930312.5) * x * y * y * y);
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

        pc_20[k] = e_0 * (std::sqrt(27.757644653320312) * x * x * x * x * x * x * x * x * x * x * y * z - std::sqrt(27.757644653320312) * x * x * x * x * x * x * x * x * y * y * y * z - std::sqrt(7105.95703125) * x * x * x * x * x * x * x * x * y * z * z * z - std::sqrt(217.61993408203125) * x * x * x * x * x * x * y * y * y * y * y * z + std::sqrt(28423.828125) * x * x * x * x * x * x * y * y * y * z * z * z + std::sqrt(7105.95703125) * x * x * x * x * x * x * y * z * z * z * z * z + std::sqrt(39.97100830078125) * x * x * x * x * y * y * y * y * y * y * y * z + std::sqrt(4547.8125) * x * x * x * x * y * y * y * y * y * z * z * z - std::sqrt(63953.61328125) * x * x * x * x * y * y * y * z * z * z * z * z + std::sqrt(89.93476867675781) * x * x * y * y * y * y * y * y * y * y * y * z - std::sqrt(28423.828125) * x * x * y * y * y * y * y * y * y * z * z * z + std::sqrt(34392.83203125) * x * x * y * y * y * y * y * z * z * z * z * z - std::sqrt(1.1103057861328125) * y * y * y * y * y * y * y * y * y * y * y * z + std::sqrt(284.23828125) * y * y * y * y * y * y * y * y * y * z * z * z - std::sqrt(284.23828125) * y * y * y * y * y * y * y * z * z * z * z * z) + e_1 * (-std::sqrt(999.2752075195312) * x * x * x * x * x * x * x * x * y * z + std::sqrt(44412.2314453125) * x * x * x * x * x * x * y * y * y * z - std::sqrt(1023257.8125) * x * x * x * x * x * x * y * z * z * z + std::sqrt(12950.606689453125) * x * x * x * x * y * y * y * y * y * z + std::sqrt(454781.25) * x * x * x * x * y * y * y * z * z * z + std::sqrt(255814.453125) * x * x * x * x * y * z * z * z * z * z - std::sqrt(15988.4033203125) * x * x * y * y * y * y * y * y * y * z - std::sqrt(1641760.3125) * x * x * y * y * y * y * y * z * z * z + std::sqrt(113695.3125) * x * x * y * y * y * z * z * z * z * z + std::sqrt(4.44122314453125) * y * y * y * y * y * y * y * y * y * z + std::sqrt(72765.0) * y * y * y * y * y * y * y * z * z * z - std::sqrt(28423.828125) * y * y * y * y * y * z * z * z * z * z) + e_2 * (-std::sqrt(2702040.1611328125) * x * x * x * x * x * x * y * z + std::sqrt(9992752.075195312) * x * x * x * x * y * y * y * z - std::sqrt(16372125.0) * x * x * x * x * y * z * z * z - std::sqrt(8457865.356445312) * x * x * y * y * y * y * y * z - std::sqrt(7276500.0) * x * x * y * y * y * z * z * z + std::sqrt(4093031.25) * x * x * y * z * z * z * z * z + std::sqrt(143895.6298828125) * y * y * y * y * y * y * y * z + std::sqrt(1819125.0) * y * y * y * y * y * z * z * z - std::sqrt(454781.25) * y * y * y * z * z * z * z * z) + e_3 * (-std::sqrt(82883882.8125) * x * x * x * x * y * z - std::sqrt(36837281.25) * x * x * y * y * y * z - std::sqrt(16372125.0) * x * x * y * z * z * z + std::sqrt(9209320.3125) * y * y * y * y * y * z + std::sqrt(1819125.0) * y * y * y * z * z * z) + e_4 * (-std::sqrt(451256695.3125) * x * x * y * z + std::sqrt(50139632.8125) * y * y * y * z);
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

        pc_21[k] = e_0 * (-std::sqrt(999.2752075195312) * x * x * x * x * x * x * x * x * x * y * z * z + std::sqrt(15988.4033203125) * x * x * x * x * x * x * x * y * y * y * z * z + std::sqrt(7105.95703125) * x * x * x * x * x * x * x * y * z * z * z * z - std::sqrt(1438.956298828125) * x * x * x * x * x * y * y * y * y * y * z * z - std::sqrt(177648.92578125) * x * x * x * x * x * y * y * y * z * z * z * z - std::sqrt(31337.2705078125) * x * x * x * y * y * y * y * y * y * y * z * z + std::sqrt(273152.98828125) * x * x * x * y * y * y * y * y * z * z * z * z + std::sqrt(359.73907470703125) * x * y * y * y * y * y * y * y * y * y * z * z - std::sqrt(2558.14453125) * x * y * y * y * y * y * y * y * z * z * z * z) + e_1 * (-std::sqrt(999.2752075195312) * x * x * x * x * x * x * x * x * x * y + std::sqrt(15988.4033203125) * x * x * x * x * x * x * x * y * y * y - std::sqrt(63953.61328125) * x * x * x * x * x * x * x * y * z * z - std::sqrt(1438.956298828125) * x * x * x * x * x * y * y * y * y * y - std::sqrt(63953.61328125) * x * x * x * x * x * y * y * y * z * z + std::sqrt(255814.453125) * x * x * x * x * x * y * z * z * z * z - std::sqrt(31337.2705078125) * x * x * x * y * y * y * y * y * y * y - std::sqrt(923490.17578125) * x * x * x * y * y * y * y * y * z * z + std::sqrt(1023257.8125) * x * x * x * y * y * y * z * z * z * z + std::sqrt(359.73907470703125) * x * y * y * y * y * y * y * y * y * y - std::sqrt(23023.30078125) * x * y * y * y * y * y * y * y * z * z + std::sqrt(255814.453125) * x * y * y * y * y * y * z * z * z * z) + e_2 * (-std::sqrt(255814.453125) * x * x * x * x * x * x * x * y + std::sqrt(1023257.8125) * x * x * x * x * x * y * y * y - std::sqrt(2302330.078125) * x * x * x * x * x * y * z * z - std::sqrt(6395361.328125) * x * x * x * y * y * y * y * y - std::sqrt(9209320.3125) * x * x * x * y * y * y * z * z + std::sqrt(16372125.0) * x * x * x * y * z * z * z * z - std::sqrt(2302330.078125) * x * y * y * y * y * y * z * z + std::sqrt(16372125.0) * x * y * y * y * z * z * z * z) + e_3 * (-std::sqrt(9209320.3125) * x * x * x * x * x * y - std::sqrt(36837281.25) * x * x * x * y * y * y - std::sqrt(9209320.3125) * x * y * y * y * y * y + std::sqrt(65488500.0) * x * y * z * z * z * z) + e_4 * (-std::sqrt(147349125.0) * x * x * x * y - std::sqrt(147349125.0) * x * y * y * y + std::sqrt(147349125.0) * x * y * z * z) + e_5 * (-std::sqrt(147349125.0) * x * y);
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

        pc_22[k] = e_0 * (-std::sqrt(33.309173583984375) * x * x * x * x * x * x * x * x * x * x * y * z + std::sqrt(1632.1495056152344) * x * x * x * x * x * x * x * x * y * y * y * z + std::sqrt(3330.9173583984375) * x * x * x * x * x * x * x * x * y * z * z * z - std::sqrt(900.6800537109375) * x * x * x * x * x * x * y * y * y * y * y * z - std::sqrt(213178.7109375) * x * x * x * x * x * x * y * y * y * z * z * z - std::sqrt(3330.9173583984375) * x * x * x * x * y * y * y * y * y * y * y * z + std::sqrt(580379.0405273438) * x * x * x * x * y * y * y * y * y * z * z * z + std::sqrt(299.7825622558594) * x * x * y * y * y * y * y * y * y * y * y * z - std::sqrt(34108.59375) * x * x * y * y * y * y * y * y * y * z * z * z - std::sqrt(1.332366943359375) * y * y * y * y * y * y * y * y * y * y * y * z + std::sqrt(133.2366943359375) * y * y * y * y * y * y * y * y * y * z * z * z) + e_1 * (std::sqrt(1199.1302490234375) * x * x * x * x * x * x * x * x * y * z - std::sqrt(306977.34375) * x * x * x * x * x * x * y * y * y * z + std::sqrt(53294.677734375) * x * x * x * x * x * x * y * z * z * z + std::sqrt(388518.20068359375) * x * x * x * x * y * y * y * y * y * z + std::sqrt(479652.099609375) * x * x * x * x * y * y * y * z * z * z - std::sqrt(76744.3359375) * x * x * y * y * y * y * y * y * y * z + std::sqrt(479652.099609375) * x * x * y * y * y * y * y * z * z * z - std::sqrt(133.2366943359375) * y * y * y * y * y * y * y * y * y * z + std::sqrt(53294.677734375) * y * y * y * y * y * y * y * z * z * z) + e_2 * (std::sqrt(7674433.59375) * x * x * x * x * y * z * z * z + std::sqrt(30697734.375) * x * x * y * y * y * z * z * z + std::sqrt(7674433.59375) * y * y * y * y * y * z * z * z) + e_3 * (std::sqrt(7674433.59375) * x * x * x * x * y * z + std::sqrt(30697734.375) * x * x * y * y * y * z + std::sqrt(122790937.5) * x * x * y * z * z * z + std::sqrt(7674433.59375) * y * y * y * y * y * z + std::sqrt(122790937.5) * y * y * y * z * z * z) + e_4 * (std::sqrt(276279609.375) * x * x * y * z + std::sqrt(276279609.375) * y * y * y * z + std::sqrt(122790937.5) * y * z * z * z) + e_5 * (std::sqrt(397842637.5) * y * z);

        pc_23[k] = e_0 * (27.0703125 * x * x * x * x * x * x * x * x * x * y * z * z - 324.84375 * x * x * x * x * x * x * x * y * y * y * z * z + 682.171875 * x * x * x * x * x * y * y * y * y * y * z * z - 324.84375 * x * x * x * y * y * y * y * y * y * y * z * z + 27.0703125 * x * y * y * y * y * y * y * y * y * y * z * z) + e_1 * (27.0703125 * x * x * x * x * x * x * x * x * x * y - 324.84375 * x * x * x * x * x * x * x * y * y * y + 682.171875 * x * x * x * x * x * y * y * y * y * y - 324.84375 * x * x * x * y * y * y * y * y * y * y + 27.0703125 * x * y * y * y * y * y * y * y * y * y);
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

        pc_24[k] = e_0 * (std::sqrt(61.06681823730469) * x * x * x * x * x * x * x * x * x * x * y * z - std::sqrt(17648.310470581055) * x * x * x * x * x * x * x * x * y * y * y * z + std::sqrt(124761.95233154297) * x * x * x * x * x * x * y * y * y * y * y * z - std::sqrt(70593.24188232422) * x * x * x * x * y * y * y * y * y * y * y * z + std::sqrt(1526.6704559326172) * x * x * y * y * y * y * y * y * y * y * y * z - std::sqrt(2.4426727294921875) * y * y * y * y * y * y * y * y * y * y * y * z) + e_1 * (-std::sqrt(2198.4054565429688) * x * x * x * x * x * x * x * x * y * z - std::sqrt(35174.4873046875) * x * x * x * x * x * x * y * y * y * z - std::sqrt(79142.59643554688) * x * x * x * x * y * y * y * y * y * z - std::sqrt(35174.4873046875) * x * x * y * y * y * y * y * y * y * z - std::sqrt(2198.4054565429688) * y * y * y * y * y * y * y * y * y * z) + e_2 * (-std::sqrt(879362.1826171875) * x * x * x * x * x * x * y * z - std::sqrt(7914259.6435546875) * x * x * x * x * y * y * y * z - std::sqrt(7914259.6435546875) * x * x * y * y * y * y * y * z - std::sqrt(879362.1826171875) * y * y * y * y * y * y * y * z) + e_3 * (-std::sqrt(56279179.6875) * x * x * x * x * y * z - std::sqrt(225116718.75) * x * x * y * y * y * z - std::sqrt(56279179.6875) * y * y * y * y * y * z) + e_4 * (-std::sqrt(506512617.1875) * x * x * y * z - std::sqrt(506512617.1875) * y * y * y * z) + e_5 * (-std::sqrt(324168075.0) * y * z);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ab_x, ab_y, ab_z : simd::cache_line_size())
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
        const auto e_6 = pe_6[k];

        pc_25[k] = e_0 * (3.9375 * x * x * x * x * x * x * x * x * x * x * y * y - 78.75 * x * x * x * x * x * x * x * x * y * y * z * z - 7.875 * x * x * x * x * x * x * y * y * y * y * y * y + 78.75 * x * x * x * x * x * x * y * y * y * y * z * z + 393.75 * x * x * x * x * x * x * y * y * z * z * z * z + 78.75 * x * x * x * x * y * y * y * y * y * y * z * z - 787.5 * x * x * x * x * y * y * y * y * z * z * z * z + 3.9375 * x * x * y * y * y * y * y * y * y * y * y * y - 78.75 * x * x * y * y * y * y * y * y * y * y * z * z + 393.75 * x * x * y * y * y * y * y * y * z * z * z * z) + e_1 * (3.9375 * x * x * x * x * x * x * x * x * x * x + 98.4375 * x * x * x * x * x * x * x * x * y * y - 78.75 * x * x * x * x * x * x * x * x * z * z - 39.375 * x * x * x * x * x * x * y * y * y * y + 630.0 * x * x * x * x * x * x * y * y * z * z + 393.75 * x * x * x * x * x * x * z * z * z * z - 39.375 * x * x * x * x * y * y * y * y * y * y - 2362.5 * x * x * x * x * y * y * y * y * z * z + 1181.25 * x * x * x * x * y * y * z * z * z * z + 98.4375 * x * x * y * y * y * y * y * y * y * y + 630.0 * x * x * y * y * y * y * y * y * z * z + 1181.25 * x * x * y * y * y * y * z * z * z * z + 3.9375 * y * y * y * y * y * y * y * y * y * y - 78.75 * y * y * y * y * y * y * y * y * z * z + 393.75 * y * y * y * y * y * y * z * z * z * z) + e_2 * (98.4375 * x * x * x * x * x * x * x * x + 1575.0 * x * x * x * x * x * x * y * y + 393.75 * x * x * x * x * x * x * z * z - 1771.875 * x * x * x * x * y * y * y * y + 1181.25 * x * x * x * x * y * y * z * z + 3543.75 * x * x * x * x * z * z * z * z + 1575.0 * x * x * y * y * y * y * y * y + 1181.25 * x * x * y * y * y * y * z * z + 7087.5 * x * x * y * y * z * z * z * z + 98.4375 * y * y * y * y * y * y * y * y + 393.75 * y * y * y * y * y * y * z * z + 3543.75 * y * y * y * y * z * z * z * z) + e_3 * (1575.0 * x * x * x * x * x * x + 4725.0 * x * x * x * x * y * y + 9450.0 * x * x * x * x * z * z + 4725.0 * x * x * y * y * y * y + 18900.0 * x * x * y * y * z * z + 9450.0 * x * x * z * z * z * z + 1575.0 * y * y * y * y * y * y + 9450.0 * y * y * y * y * z * z + 9450.0 * y * y * z * z * z * z) + e_4 * (9450.0 * x * x * x * x + 18900.0 * x * x * y * y + 33075.0 * x * x * z * z + 9450.0 * y * y * y * y + 33075.0 * y * y * z * z + 4725.0 * z * z * z * z) + e_5 * (21735.0 * x * x + 21735.0 * y * y + 18900.0 * z * z) + e_6 * (10395.0);
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

        pc_26[k] = e_0 * (std::sqrt(261.62841796875) * x * x * x * x * x * x * x * x * x * y * y * z + std::sqrt(116.279296875) * x * x * x * x * x * x * x * y * y * y * y * z - std::sqrt(41976.826171875) * x * x * x * x * x * x * x * y * y * z * z * z - std::sqrt(465.1171875) * x * x * x * x * x * y * y * y * y * y * y * z + std::sqrt(4664.091796875) * x * x * x * x * x * y * y * y * y * z * z * z + std::sqrt(186046.875) * x * x * x * x * x * y * y * z * z * z * z * z - std::sqrt(116.279296875) * x * x * x * y * y * y * y * y * y * y * y * z + std::sqrt(41976.826171875) * x * x * x * y * y * y * y * y * y * z * z * z - std::sqrt(330750.0) * x * x * x * y * y * y * y * z * z * z * z * z + std::sqrt(29.06982421875) * x * y * y * y * y * y * y * y * y * y * y * z - std::sqrt(4664.091796875) * x * y * y * y * y * y * y * y * y * z * z * z + std::sqrt(20671.875) * x * y * y * y * y * y * y * z * z * z * z * z) + e_1 * (std::sqrt(261.62841796875) * x * x * x * x * x * x * x * x * x * z + std::sqrt(1046.513671875) * x * x * x * x * x * x * x * y * y * z - std::sqrt(41976.826171875) * x * x * x * x * x * x * x * z * z * z + std::sqrt(11627.9296875) * x * x * x * x * x * y * y * y * y * z + std::sqrt(176860.810546875) * x * x * x * x * x * y * y * z * z * z + std::sqrt(186046.875) * x * x * x * x * x * z * z * z * z * z + std::sqrt(9418.623046875) * x * x * x * y * y * y * y * y * y * z - std::sqrt(3979658.935546875) * x * x * x * y * y * y * y * z * z * z + std::sqrt(744187.5) * x * x * x * y * y * z * z * z * z * z + std::sqrt(29.06982421875) * x * y * y * y * y * y * y * y * y * z + std::sqrt(19651.201171875) * x * y * y * y * y * y * y * z * z * z + std::sqrt(186046.875) * x * y * y * y * y * z * z * z * z * z) + e_2 * (std::sqrt(1674421.875) * x * x * x * x * x * y * y * z + std::sqrt(46511.71875) * x * x * x * x * x * z * z * z - std::sqrt(2976750.0) * x * x * x * y * y * y * y * z + std::sqrt(186046.875) * x * x * x * y * y * z * z * z + std::sqrt(6697687.5) * x * x * x * z * z * z * z * z + std::sqrt(186046.875) * x * y * y * y * y * y * y * z + std::sqrt(46511.71875) * x * y * y * y * y * z * z * z + std::sqrt(6697687.5) * x * y * y * z * z * z * z * z) + e_3 * (std::sqrt(418605.46875) * x * x * x * x * x * z + std::sqrt(1674421.875) * x * x * x * y * y * z + std::sqrt(90046687.5) * x * x * x * z * z * z + std::sqrt(418605.46875) * x * y * y * y * y * z + std::sqrt(90046687.5) * x * y * y * z * z * z + std::sqrt(11907000.0) * x * z * z * z * z * z) + e_4 * (std::sqrt(82046671.875) * x * x * x * z + std::sqrt(82046671.875) * x * y * y * z + std::sqrt(328186687.5) * x * z * z * z) + e_5 * (std::sqrt(328186687.5) * x * z);
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

        pc_27[k] = e_0 * (-std::sqrt(3.22998046875) * x * x * x * x * x * x * x * x * x * x * y * y - std::sqrt(12.919921875) * x * x * x * x * x * x * x * x * y * y * y * y + std::sqrt(2183.466796875) * x * x * x * x * x * x * x * x * y * y * z * z + std::sqrt(2183.466796875) * x * x * x * x * x * x * y * y * y * y * z * z - std::sqrt(100051.875) * x * x * x * x * x * x * y * y * z * z * z * z + std::sqrt(12.919921875) * x * x * x * x * y * y * y * y * y * y * y * y - std::sqrt(2183.466796875) * x * x * x * x * y * y * y * y * y * y * z * z + std::sqrt(82687.5) * x * x * x * x * y * y * z * z * z * z * z * z + std::sqrt(3.22998046875) * x * x * y * y * y * y * y * y * y * y * y * y - std::sqrt(2183.466796875) * x * x * y * y * y * y * y * y * y * y * z * z + std::sqrt(100051.875) * x * x * y * y * y * y * y * y * z * z * z * z - std::sqrt(82687.5) * x * x * y * y * y * y * z * z * z * z * z * z) + e_1 * (-std::sqrt(3.22998046875) * x * x * x * x * x * x * x * x * x * x - std::sqrt(3104.01123046875) * x * x * x * x * x * x * x * x * y * y + std::sqrt(2183.466796875) * x * x * x * x * x * x * x * x * z * z - std::sqrt(2906.982421875) * x * x * x * x * x * x * y * y * y * y - std::sqrt(95555.7421875) * x * x * x * x * x * x * y * y * z * z - std::sqrt(100051.875) * x * x * x * x * x * x * z * z * z * z + std::sqrt(2906.982421875) * x * x * x * x * y * y * y * y * y * y - std::sqrt(186046.875) * x * x * x * x * y * y * z * z * z * z + std::sqrt(82687.5) * x * x * x * x * z * z * z * z * z * z + std::sqrt(3104.01123046875) * x * x * y * y * y * y * y * y * y * y + std::sqrt(95555.7421875) * x * x * y * y * y * y * y * y * z * z + std::sqrt(186046.875) * x * x * y * y * y * y * z * z * z * z + std::sqrt(3.22998046875) * y * y * y * y * y * y * y * y * y * y - std::sqrt(2183.466796875) * y * y * y * y * y * y * y * y * z * z + std::sqrt(100051.875) * y * y * y * y * y * y * z * z * z * z - std::sqrt(82687.5) * y * y * y * y * z * z * z * z * z * z) + e_2 * (-std::sqrt(2018.73779296875) * x * x * x * x * x * x * x * x - std::sqrt(1201875.732421875) * x * x * x * x * x * x * y * y - std::sqrt(201873.779296875) * x * x * x * x * x * x * z * z - std::sqrt(13049444.091796875) * x * x * x * x * y * y * z * z - std::sqrt(186046.875) * x * x * x * x * z * z * z * z + std::sqrt(1201875.732421875) * x * x * y * y * y * y * y * y + std::sqrt(13049444.091796875) * x * x * y * y * y * y * z * z + std::sqrt(744187.5) * x * x * z * z * z * z * z * z + std::sqrt(2018.73779296875) * y * y * y * y * y * y * y * y + std::sqrt(201873.779296875) * y * y * y * y * y * y * z * z + std::sqrt(186046.875) * y * y * y * y * z * z * z * z - std::sqrt(744187.5) * y * y * z * z * z * z * z * z) + e_3 * (-std::sqrt(873386.71875) * x * x * x * x * x * x - std::sqrt(44697761.71875) * x * x * x * x * y * y - std::sqrt(18604687.5) * x * x * x * x * z * z + std::sqrt(44697761.71875) * x * x * y * y * y * y + std::sqrt(11907000.0) * x * x * z * z * z * z + std::sqrt(873386.71875) * y * y * y * y * y * y + std::sqrt(18604687.5) * y * y * y * y * z * z - std::sqrt(11907000.0) * y * y * z * z * z * z) + e_4 * (-std::sqrt(39116355.46875) * x * x * x * x - std::sqrt(1674421.875) * x * x * z * z + std::sqrt(39116355.46875) * y * y * y * y + std::sqrt(1674421.875) * y * y * z * z) + e_5 * (-std::sqrt(60279187.5) * x * x + std::sqrt(60279187.5) * y * y);
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

        pc_28[k] = e_0 * (-std::sqrt(32.2998046875) * x * x * x * x * x * x * x * x * x * y * y * z - std::sqrt(129.19921875) * x * x * x * x * x * x * x * y * y * y * y * z + std::sqrt(6330.76171875) * x * x * x * x * x * x * x * y * y * z * z * z + std::sqrt(6330.76171875) * x * x * x * x * x * y * y * y * y * z * z * z - std::sqrt(55896.75) * x * x * x * x * x * y * y * z * z * z * z * z + std::sqrt(129.19921875) * x * x * x * y * y * y * y * y * y * y * y * z - std::sqrt(6330.76171875) * x * x * x * y * y * y * y * y * y * z * z * z + std::sqrt(8268.75) * x * x * x * y * y * z * z * z * z * z * z * z + std::sqrt(32.2998046875) * x * y * y * y * y * y * y * y * y * y * y * z - std::sqrt(6330.76171875) * x * y * y * y * y * y * y * y * y * z * z * z + std::sqrt(55896.75) * x * y * y * y * y * y * y * z * z * z * z * z - std::sqrt(8268.75) * x * y * y * y * y * z * z * z * z * z * z * z) + e_1 * (-std::sqrt(32.2998046875) * x * x * x * x * x * x * x * x * x * z - std::sqrt(1162.79296875) * x * x * x * x * x * x * x * y * y * z + std::sqrt(6330.76171875) * x * x * x * x * x * x * x * z * z * z - std::sqrt(46640.91796875) * x * x * x * x * x * y * y * z * z * z - std::sqrt(55896.75) * x * x * x * x * x * z * z * z * z * z + std::sqrt(6330.76171875) * x * x * x * y * y * y * y * y * y * z - std::sqrt(158269.04296875) * x * x * x * y * y * y * y * z * z * z - std::sqrt(206718.75) * x * x * x * y * y * z * z * z * z * z + std::sqrt(8268.75) * x * x * x * z * z * z * z * z * z * z + std::sqrt(2616.2841796875) * x * y * y * y * y * y * y * y * y * z - std::sqrt(10465.13671875) * x * y * y * y * y * y * y * z * z * z + std::sqrt(2679075.0) * x * y * y * y * y * z * z * z * z * z - std::sqrt(74418.75) * x * y * y * z * z * z * z * z * z * z) + e_2 * (-std::sqrt(465117.1875) * x * x * x * x * x * y * y * z - std::sqrt(206718.75) * x * x * x * x * x * z * z * z - std::sqrt(20671875.0) * x * x * x * y * y * z * z * z - std::sqrt(206718.75) * x * x * x * z * z * z * z * z + std::sqrt(465117.1875) * x * y * y * y * y * y * y * z + std::sqrt(46511718.75) * x * y * y * y * y * z * z * z + std::sqrt(1860468.75) * x * y * y * z * z * z * z * z) + e_3 * (-std::sqrt(465117.1875) * x * x * x * x * x * z - std::sqrt(46511718.75) * x * x * x * y * y * z - std::sqrt(20671875.0) * x * x * x * z * z * z + std::sqrt(104651367.1875) * x * y * y * y * y * z + std::sqrt(186046875.0) * x * y * y * z * z * z) + e_4 * (-std::sqrt(46511718.75) * x * x * x * z + std::sqrt(418605468.75) * x * y * y * z);
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

        pc_29[k] = e_0 * (std::sqrt(0.384521484375) * x * x * x * x * x * x * x * x * x * x * x * y + std::sqrt(3.460693359375) * x * x * x * x * x * x * x * x * x * y * y * y - std::sqrt(301.46484375) * x * x * x * x * x * x * x * x * x * y * z * z + std::sqrt(1.5380859375) * x * x * x * x * x * x * x * y * y * y * y * y - std::sqrt(1205.859375) * x * x * x * x * x * x * x * y * y * y * z * z + std::sqrt(16002.24609375) * x * x * x * x * x * x * x * y * z * z * z * z - std::sqrt(1.5380859375) * x * x * x * x * x * y * y * y * y * y * y * y + std::sqrt(16002.24609375) * x * x * x * x * x * y * y * y * z * z * z * z - std::sqrt(22743.0) * x * x * x * x * x * y * z * z * z * z * z * z - std::sqrt(3.460693359375) * x * x * x * y * y * y * y * y * y * y * y * y + std::sqrt(1205.859375) * x * x * x * y * y * y * y * y * y * y * z * z - std::sqrt(16002.24609375) * x * x * x * y * y * y * y * y * z * z * z * z + std::sqrt(393.75) * x * x * x * y * z * z * z * z * z * z * z * z - std::sqrt(0.384521484375) * x * y * y * y * y * y * y * y * y * y * y * y + std::sqrt(301.46484375) * x * y * y * y * y * y * y * y * y * y * z * z - std::sqrt(16002.24609375) * x * y * y * y * y * y * y * y * z * z * z * z + std::sqrt(22743.0) * x * y * y * y * y * y * z * z * z * z * z * z - std::sqrt(393.75) * x * y * y * y * z * z * z * z * z * z * z * z) + e_1 * (std::sqrt(498.33984375) * x * x * x * x * x * x * x * x * x * y + std::sqrt(1993.359375) * x * x * x * x * x * x * x * y * y * y + std::sqrt(885.9375) * x * x * x * x * x * x * x * y * z * z + std::sqrt(885.9375) * x * x * x * x * x * y * y * y * z * z + std::sqrt(598893.75) * x * x * x * x * x * y * z * z * z * z - std::sqrt(1993.359375) * x * x * x * y * y * y * y * y * y * y - std::sqrt(885.9375) * x * x * x * y * y * y * y * y * z * z - std::sqrt(907200.0) * x * x * x * y * z * z * z * z * z * z - std::sqrt(498.33984375) * x * y * y * y * y * y * y * y * y * y - std::sqrt(885.9375) * x * y * y * y * y * y * y * y * z * z - std::sqrt(598893.75) * x * y * y * y * y * y * z * z * z * z + std::sqrt(907200.0) * x * y * y * y * z * z * z * z * z * z) + e_2 * (std::sqrt(233942.87109375) * x * x * x * x * x * x * x * y + std::sqrt(233942.87109375) * x * x * x * x * x * y * y * y + std::sqrt(7176093.75) * x * x * x * x * x * y * z * z - std::sqrt(233942.87109375) * x * x * x * y * y * y * y * y - std::sqrt(10719843.75) * x * x * x * y * z * z * z * z - std::sqrt(233942.87109375) * x * y * y * y * y * y * y * y - std::sqrt(7176093.75) * x * y * y * y * y * y * z * z + std::sqrt(10719843.75) * x * y * y * y * z * z * z * z) + e_3 * (std::sqrt(22680000.0) * x * x * x * x * x * y + std::sqrt(5670000.0) * x * x * x * y * z * z - std::sqrt(22680000.0) * x * y * y * y * y * y - std::sqrt(5670000.0) * x * y * y * y * z * z) + e_4 * (std::sqrt(156279375.0) * x * x * x * y - std::sqrt(156279375.0) * x * y * y * y);
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

        pc_30[k] = e_0 * (-std::sqrt(32.2998046875) * x * x * x * x * x * x * x * x * x * x * y * z - std::sqrt(129.19921875) * x * x * x * x * x * x * x * x * y * y * y * z + std::sqrt(6330.76171875) * x * x * x * x * x * x * x * x * y * z * z * z + std::sqrt(6330.76171875) * x * x * x * x * x * x * y * y * y * z * z * z - std::sqrt(55896.75) * x * x * x * x * x * x * y * z * z * z * z * z + std::sqrt(129.19921875) * x * x * x * x * y * y * y * y * y * y * y * z - std::sqrt(6330.76171875) * x * x * x * x * y * y * y * y * y * z * z * z + std::sqrt(8268.75) * x * x * x * x * y * z * z * z * z * z * z * z + std::sqrt(32.2998046875) * x * x * y * y * y * y * y * y * y * y * y * z - std::sqrt(6330.76171875) * x * x * y * y * y * y * y * y * y * z * z * z + std::sqrt(55896.75) * x * x * y * y * y * y * y * z * z * z * z * z - std::sqrt(8268.75) * x * x * y * y * y * z * z * z * z * z * z * z) + e_1 * (-std::sqrt(2616.2841796875) * x * x * x * x * x * x * x * x * y * z - std::sqrt(6330.76171875) * x * x * x * x * x * x * y * y * y * z + std::sqrt(10465.13671875) * x * x * x * x * x * x * y * z * z * z + std::sqrt(158269.04296875) * x * x * x * x * y * y * y * z * z * z - std::sqrt(2679075.0) * x * x * x * x * y * z * z * z * z * z + std::sqrt(1162.79296875) * x * x * y * y * y * y * y * y * y * z + std::sqrt(46640.91796875) * x * x * y * y * y * y * y * z * z * z + std::sqrt(206718.75) * x * x * y * y * y * z * z * z * z * z + std::sqrt(74418.75) * x * x * y * z * z * z * z * z * z * z + std::sqrt(32.2998046875) * y * y * y * y * y * y * y * y * y * z - std::sqrt(6330.76171875) * y * y * y * y * y * y * y * z * z * z + std::sqrt(55896.75) * y * y * y * y * y * z * z * z * z * z - std::sqrt(8268.75) * y * y * y * z * z * z * z * z * z * z) + e_2 * (-std::sqrt(465117.1875) * x * x * x * x * x * x * y * z - std::sqrt(46511718.75) * x * x * x * x * y * z * z * z + std::sqrt(465117.1875) * x * x * y * y * y * y * y * z + std::sqrt(20671875.0) * x * x * y * y * y * z * z * z - std::sqrt(1860468.75) * x * x * y * z * z * z * z * z + std::sqrt(206718.75) * y * y * y * y * y * z * z * z + std::sqrt(206718.75) * y * y * y * z * z * z * z * z) + e_3 * (-std::sqrt(104651367.1875) * x * x * x * x * y * z + std::sqrt(46511718.75) * x * x * y * y * y * z - std::sqrt(186046875.0) * x * x * y * z * z * z + std::sqrt(465117.1875) * y * y * y * y * y * z + std::sqrt(20671875.0) * y * y * y * z * z * z) + e_4 * (-std::sqrt(418605468.75) * x * x * y * z + std::sqrt(46511718.75) * y * y * y * z);
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

        pc_31[k] = e_0 * (-std::sqrt(0.8074951171875) * x * x * x * x * x * x * x * x * x * x * x * y - std::sqrt(0.8074951171875) * x * x * x * x * x * x * x * x * x * y * y * y + std::sqrt(545.86669921875) * x * x * x * x * x * x * x * x * x * y * z * z + std::sqrt(3.22998046875) * x * x * x * x * x * x * x * y * y * y * y * y - std::sqrt(25012.96875) * x * x * x * x * x * x * x * y * z * z * z * z + std::sqrt(3.22998046875) * x * x * x * x * x * y * y * y * y * y * y * y - std::sqrt(2183.466796875) * x * x * x * x * x * y * y * y * y * y * z * z + std::sqrt(25012.96875) * x * x * x * x * x * y * y * y * z * z * z * z + std::sqrt(20671.875) * x * x * x * x * x * y * z * z * z * z * z * z - std::sqrt(0.8074951171875) * x * x * x * y * y * y * y * y * y * y * y * y + std::sqrt(25012.96875) * x * x * x * y * y * y * y * y * z * z * z * z - std::sqrt(82687.5) * x * x * x * y * y * y * z * z * z * z * z * z - std::sqrt(0.8074951171875) * x * y * y * y * y * y * y * y * y * y * y * y + std::sqrt(545.86669921875) * x * y * y * y * y * y * y * y * y * y * z * z - std::sqrt(25012.96875) * x * y * y * y * y * y * y * y * z * z * z * z + std::sqrt(20671.875) * x * y * y * y * y * y * z * z * z * z * z * z) + e_1 * (-std::sqrt(826.875) * x * x * x * x * x * x * x * x * x * y - std::sqrt(206.71875) * x * x * x * x * x * x * x * y * y * y - std::sqrt(11627.9296875) * x * x * x * x * x * x * x * y * z * z + std::sqrt(826.875) * x * x * x * x * x * y * y * y * y * y + std::sqrt(231990.1171875) * x * x * x * x * x * y * y * y * z * z - std::sqrt(476280.0) * x * x * x * x * x * y * z * z * z * z - std::sqrt(206.71875) * x * x * x * y * y * y * y * y * y * y + std::sqrt(231990.1171875) * x * x * x * y * y * y * y * y * z * z - std::sqrt(1323000.0) * x * x * x * y * y * y * z * z * z * z + std::sqrt(330750.0) * x * x * x * y * z * z * z * z * z * z - std::sqrt(826.875) * x * y * y * y * y * y * y * y * y * y - std::sqrt(11627.9296875) * x * y * y * y * y * y * y * y * z * z - std::sqrt(476280.0) * x * y * y * y * y * y * z * z * z * z + std::sqrt(330750.0) * x * y * y * y * z * z * z * z * z * z) + e_2 * (-std::sqrt(351744.873046875) * x * x * x * x * x * x * x * y + std::sqrt(54586.669921875) * x * x * x * x * x * y * y * y - std::sqrt(6151174.8046875) * x * x * x * x * x * y * z * z + std::sqrt(54586.669921875) * x * x * x * y * y * y * y * y + std::sqrt(1865636.71875) * x * x * x * y * y * y * z * z - std::sqrt(744187.5) * x * x * x * y * z * z * z * z - std::sqrt(351744.873046875) * x * y * y * y * y * y * y * y - std::sqrt(6151174.8046875) * x * y * y * y * y * y * z * z - std::sqrt(744187.5) * x * y * y * y * z * z * z * z + std::sqrt(2976750.0) * x * y * z * z * z * z * z * z) + e_3 * (-std::sqrt(22511671.875) * x * x * x * x * x * y + std::sqrt(4051687.5) * x * x * x * y * y * y - std::sqrt(74418750.0) * x * x * x * y * z * z - std::sqrt(22511671.875) * x * y * y * y * y * y - std::sqrt(74418750.0) * x * y * y * y * z * z + std::sqrt(47628000.0) * x * y * z * z * z * z) + e_4 * (-std::sqrt(156465421.875) * x * x * x * y - std::sqrt(156465421.875) * x * y * y * y - std::sqrt(6697687.5) * x * y * z * z) + e_5 * (-std::sqrt(241116750.0) * x * y);
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

        pc_32[k] = e_0 * (std::sqrt(29.06982421875) * x * x * x * x * x * x * x * x * x * x * y * z - std::sqrt(116.279296875) * x * x * x * x * x * x * x * x * y * y * y * z - std::sqrt(4664.091796875) * x * x * x * x * x * x * x * x * y * z * z * z - std::sqrt(465.1171875) * x * x * x * x * x * x * y * y * y * y * y * z + std::sqrt(41976.826171875) * x * x * x * x * x * x * y * y * y * z * z * z + std::sqrt(20671.875) * x * x * x * x * x * x * y * z * z * z * z * z + std::sqrt(116.279296875) * x * x * x * x * y * y * y * y * y * y * y * z + std::sqrt(4664.091796875) * x * x * x * x * y * y * y * y * y * z * z * z - std::sqrt(330750.0) * x * x * x * x * y * y * y * z * z * z * z * z + std::sqrt(261.62841796875) * x * x * y * y * y * y * y * y * y * y * y * z - std::sqrt(41976.826171875) * x * x * y * y * y * y * y * y * y * z * z * z + std::sqrt(186046.875) * x * x * y * y * y * y * y * z * z * z * z * z) + e_1 * (std::sqrt(29.06982421875) * x * x * x * x * x * x * x * x * y * z + std::sqrt(9418.623046875) * x * x * x * x * x * x * y * y * y * z + std::sqrt(19651.201171875) * x * x * x * x * x * x * y * z * z * z + std::sqrt(11627.9296875) * x * x * x * x * y * y * y * y * y * z - std::sqrt(3979658.935546875) * x * x * x * x * y * y * y * z * z * z + std::sqrt(186046.875) * x * x * x * x * y * z * z * z * z * z + std::sqrt(1046.513671875) * x * x * y * y * y * y * y * y * y * z + std::sqrt(176860.810546875) * x * x * y * y * y * y * y * z * z * z + std::sqrt(744187.5) * x * x * y * y * y * z * z * z * z * z + std::sqrt(261.62841796875) * y * y * y * y * y * y * y * y * y * z - std::sqrt(41976.826171875) * y * y * y * y * y * y * y * z * z * z + std::sqrt(186046.875) * y * y * y * y * y * z * z * z * z * z) + e_2 * (std::sqrt(186046.875) * x * x * x * x * x * x * y * z - std::sqrt(2976750.0) * x * x * x * x * y * y * y * z + std::sqrt(46511.71875) * x * x * x * x * y * z * z * z + std::sqrt(1674421.875) * x * x * y * y * y * y * y * z + std::sqrt(186046.875) * x * x * y * y * y * z * z * z + std::sqrt(6697687.5) * x * x * y * z * z * z * z * z + std::sqrt(46511.71875) * y * y * y * y * y * z * z * z + std::sqrt(6697687.5) * y * y * y * z * z * z * z * z) + e_3 * (std::sqrt(418605.46875) * x * x * x * x * y * z + std::sqrt(1674421.875) * x * x * y * y * y * z + std::sqrt(90046687.5) * x * x * y * z * z * z + std::sqrt(418605.46875) * y * y * y * y * y * z + std::sqrt(90046687.5) * y * y * y * z * z * z + std::sqrt(11907000.0) * y * z * z * z * z * z) + e_4 * (std::sqrt(82046671.875) * x * x * y * z + std::sqrt(82046671.875) * y * y * y * z + std::sqrt(328186687.5) * y * z * z * z) + e_5 * (std::sqrt(328186687.5) * y * z);
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

        pc_33[k] = e_0 * (0.984375 * x * x * x * x * x * x * x * x * x * x * x * y - 4.921875 * x * x * x * x * x * x * x * x * x * y * y * y - 19.6875 * x * x * x * x * x * x * x * x * x * y * z * z - 5.90625 * x * x * x * x * x * x * x * y * y * y * y * y + 118.125 * x * x * x * x * x * x * x * y * y * y * z * z + 98.4375 * x * x * x * x * x * x * x * y * z * z * z * z + 5.90625 * x * x * x * x * x * y * y * y * y * y * y * y - 689.0625 * x * x * x * x * x * y * y * y * z * z * z * z + 4.921875 * x * x * x * y * y * y * y * y * y * y * y * y - 118.125 * x * x * x * y * y * y * y * y * y * y * z * z + 689.0625 * x * x * x * y * y * y * y * y * z * z * z * z - 0.984375 * x * y * y * y * y * y * y * y * y * y * y * y + 19.6875 * x * y * y * y * y * y * y * y * y * y * z * z - 98.4375 * x * y * y * y * y * y * y * y * z * z * z * z) + e_1 * (19.6875 * x * x * x * x * x * x * x * x * x * y - 118.125 * x * x * x * x * x * x * x * y * y * y + 236.25 * x * x * x * x * x * x * x * y * z * z - 1653.75 * x * x * x * x * x * y * y * y * z * z + 118.125 * x * x * x * y * y * y * y * y * y * y + 1653.75 * x * x * x * y * y * y * y * y * z * z - 19.6875 * x * y * y * y * y * y * y * y * y * y - 236.25 * x * y * y * y * y * y * y * y * z * z) + e_2 * (295.3125 * x * x * x * x * x * x * x * y - 2067.1875 * x * x * x * x * x * y * y * y + 2067.1875 * x * x * x * y * y * y * y * y - 295.3125 * x * y * y * y * y * y * y * y);
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

        pc_34[k] = e_0 * (-std::sqrt(21.31787109375) * x * x * x * x * x * x * x * x * x * x * y * z + std::sqrt(2131.787109375) * x * x * x * x * x * x * x * x * y * y * y * z + std::sqrt(2131.787109375) * x * x * x * x * x * x * x * x * y * z * z * z - std::sqrt(341.0859375) * x * x * x * x * x * x * y * y * y * y * y * z - std::sqrt(257946.240234375) * x * x * x * x * x * x * y * y * y * z * z * z - std::sqrt(2131.787109375) * x * x * x * x * y * y * y * y * y * y * y * z + std::sqrt(479652.099609375) * x * x * x * x * y * y * y * y * y * z * z * z + std::sqrt(532.94677734375) * x * x * y * y * y * y * y * y * y * y * y * z - std::sqrt(53294.677734375) * x * x * y * y * y * y * y * y * y * z * z * z) + e_1 * (std::sqrt(4796.52099609375) * x * x * x * x * x * x * x * x * y * z - std::sqrt(172674.755859375) * x * x * x * x * x * x * y * y * y * z - std::sqrt(53294.677734375) * x * x * x * x * x * x * y * z * z * z + std::sqrt(690699.0234375) * x * x * x * x * y * y * y * y * y * z - std::sqrt(479652.099609375) * x * x * x * x * y * y * y * z * z * z - std::sqrt(19186.083984375) * x * x * y * y * y * y * y * y * y * z - std::sqrt(479652.099609375) * x * x * y * y * y * y * y * z * z * z + std::sqrt(532.94677734375) * y * y * y * y * y * y * y * y * y * z - std::sqrt(53294.677734375) * y * y * y * y * y * y * y * z * z * z) + e_2 * (-std::sqrt(7674433.59375) * x * x * x * x * y * z * z * z - std::sqrt(30697734.375) * x * x * y * y * y * z * z * z - std::sqrt(7674433.59375) * y * y * y * y * y * z * z * z) + e_3 * (-std::sqrt(7674433.59375) * x * x * x * x * y * z - std::sqrt(30697734.375) * x * x * y * y * y * z - std::sqrt(122790937.5) * x * x * y * z * z * z - std::sqrt(7674433.59375) * y * y * y * y * y * z - std::sqrt(122790937.5) * y * y * y * z * z * z) + e_4 * (-std::sqrt(276279609.375) * x * x * y * z - std::sqrt(276279609.375) * y * y * y * z - std::sqrt(122790937.5) * y * z * z * z) + e_5 * (-std::sqrt(397842637.5) * y * z);
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

        pc_35[k] = e_0 * (-std::sqrt(1.7764892578125) * x * x * x * x * x * x * x * x * x * x * x * y + std::sqrt(399.7100830078125) * x * x * x * x * x * x * x * x * x * y * y * y + std::sqrt(177.64892578125) * x * x * x * x * x * x * x * x * x * y * z * z - std::sqrt(348.19189453125) * x * x * x * x * x * x * x * y * y * y * y * y - std::sqrt(45478.125) * x * x * x * x * x * x * x * y * y * y * z * z - std::sqrt(348.19189453125) * x * x * x * x * x * y * y * y * y * y * y * y + std::sqrt(159884.033203125) * x * x * x * x * x * y * y * y * y * y * z * z + std::sqrt(399.7100830078125) * x * x * x * y * y * y * y * y * y * y * y * y - std::sqrt(45478.125) * x * x * x * y * y * y * y * y * y * y * z * z - std::sqrt(1.7764892578125) * x * y * y * y * y * y * y * y * y * y * y * y + std::sqrt(177.64892578125) * x * y * y * y * y * y * y * y * y * y * z * z) + e_1 * (std::sqrt(102325.78125) * x * x * x * x * x * x * x * y * y * y - std::sqrt(25581.4453125) * x * x * x * x * x * x * x * y * z * z - std::sqrt(147349.125) * x * x * x * x * x * y * y * y * y * y - std::sqrt(230233.0078125) * x * x * x * x * x * y * y * y * z * z + std::sqrt(102325.78125) * x * x * x * y * y * y * y * y * y * y - std::sqrt(230233.0078125) * x * x * x * y * y * y * y * y * z * z - std::sqrt(25581.4453125) * x * y * y * y * y * y * y * y * z * z) + e_2 * (std::sqrt(159884.033203125) * x * x * x * x * x * x * x * y + std::sqrt(1438956.298828125) * x * x * x * x * x * y * y * y - std::sqrt(5755825.1953125) * x * x * x * x * x * y * z * z + std::sqrt(1438956.298828125) * x * x * x * y * y * y * y * y - std::sqrt(23023300.78125) * x * x * x * y * y * y * z * z + std::sqrt(159884.033203125) * x * y * y * y * y * y * y * y - std::sqrt(5755825.1953125) * x * y * y * y * y * y * z * z) + e_3 * (std::sqrt(10232578.125) * x * x * x * x * x * y + std::sqrt(40930312.5) * x * x * x * y * y * y - std::sqrt(163721250.0) * x * x * x * y * z * z + std::sqrt(10232578.125) * x * y * y * y * y * y - std::sqrt(163721250.0) * x * y * y * y * z * z) + e_4 * (std::sqrt(92093203.125) * x * x * x * y + std::sqrt(92093203.125) * x * y * y * y - std::sqrt(368372812.5) * x * y * z * z) + e_5 * (std::sqrt(58939650.0) * x * y);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ab_x, ab_y, ab_z : simd::cache_line_size())
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
        const auto e_6 = pe_6[k];

        pc_36[k] = e_0 * (66.4453125 * x * x * x * x * x * x * x * x * y * y * z * z + 88.59375 * x * x * x * x * x * x * y * y * y * y * z * z - 354.375 * x * x * x * x * x * x * y * y * z * z * z * z - 14.765625 * x * x * x * x * y * y * y * y * y * y * z * z - 118.125 * x * x * x * x * y * y * y * y * z * z * z * z + 472.5 * x * x * x * x * y * y * z * z * z * z * z * z - 29.53125 * x * x * y * y * y * y * y * y * y * y * z * z + 196.875 * x * x * y * y * y * y * y * y * z * z * z * z - 315.0 * x * x * y * y * y * y * z * z * z * z * z * z + 7.3828125 * y * y * y * y * y * y * y * y * y * y * z * z - 39.375 * y * y * y * y * y * y * y * y * z * z * z * z + 52.5 * y * y * y * y * y * y * z * z * z * z * z * z) + e_1 * (66.4453125 * x * x * x * x * x * x * x * x * y * y + 66.4453125 * x * x * x * x * x * x * x * x * z * z + 88.59375 * x * x * x * x * x * x * y * y * y * y + 265.78125 * x * x * x * x * x * x * y * y * z * z - 354.375 * x * x * x * x * x * x * z * z * z * z - 14.765625 * x * x * x * x * y * y * y * y * y * y + 398.671875 * x * x * x * x * y * y * y * y * z * z + 1063.125 * x * x * x * x * y * y * z * z * z * z + 472.5 * x * x * x * x * z * z * z * z * z * z - 29.53125 * x * x * y * y * y * y * y * y * y * y + 265.78125 * x * x * y * y * y * y * y * y * z * z - 2480.625 * x * x * y * y * y * y * z * z * z * z + 945.0 * x * x * y * y * z * z * z * z * z * z + 7.3828125 * y * y * y * y * y * y * y * y * y * y + 66.4453125 * y * y * y * y * y * y * y * y * z * z - 118.125 * y * y * y * y * y * y * z * z * z * z + 472.5 * y * y * y * y * z * z * z * z * z * z) + e_2 * (66.4453125 * x * x * x * x * x * x * x * x + 1328.90625 * x * x * x * x * x * x * y * y + 753.046875 * x * x * x * x * y * y * y * y + 6378.75 * x * x * x * x * y * y * z * z + 1417.5 * x * x * x * x * z * z * z * z - 324.84375 * x * x * y * y * y * y * y * y - 4252.5 * x * x * y * y * y * y * z * z + 2835.0 * x * x * y * y * z * z * z * z + 1890.0 * x * x * z * z * z * z * z * z + 184.5703125 * y * y * y * y * y * y * y * y + 708.75 * y * y * y * y * y * y * z * z + 1417.5 * y * y * y * y * z * z * z * z + 1890.0 * y * y * z * z * z * z * z * z) + e_3 * (1063.125 * x * x * x * x * x * x + 10276.875 * x * x * x * x * y * y + 4961.25 * x * x * x * x * z * z - 1535.625 * x * x * y * y * y * y + 9922.5 * x * x * y * y * z * z + 13230.0 * x * x * z * z * z * z + 1850.625 * y * y * y * y * y * y + 4961.25 * y * y * y * y * z * z + 13230.0 * y * y * z * z * z * z + 1260.0 * z * z * z * z * z * z) + e_4 * (7796.25 * x * x * x * x + 15592.5 * x * x * y * y + 29767.5 * x * x * z * z + 7796.25 * y * y * y * y + 29767.5 * y * y * z * z + 11340.0 * z * z * z * z) + e_5 * (18427.5 * x * x + 18427.5 * y * y + 25515.0 * z * z) + e_6 * (10395.0);
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

        pc_37[k] = e_0 * (-7.3828125 * x * x * x * x * x * x * x * x * x * y * y * z - 19.6875 * x * x * x * x * x * x * x * y * y * y * y * z + 137.8125 * x * x * x * x * x * x * x * y * y * z * z * z - 14.765625 * x * x * x * x * x * y * y * y * y * y * y * z + 229.6875 * x * x * x * x * x * y * y * y * y * z * z * z - 433.125 * x * x * x * x * x * y * y * z * z * z * z * z + 45.9375 * x * x * x * y * y * y * y * y * y * z * z * z - 288.75 * x * x * x * y * y * y * y * z * z * z * z * z + 315.0 * x * x * x * y * y * z * z * z * z * z * z * z + 2.4609375 * x * y * y * y * y * y * y * y * y * y * y * z - 45.9375 * x * y * y * y * y * y * y * y * y * z * z * z + 144.375 * x * y * y * y * y * y * y * z * z * z * z * z - 105.0 * x * y * y * y * y * z * z * z * z * z * z * z) + e_1 * (-7.3828125 * x * x * x * x * x * x * x * x * x * z + 29.53125 * x * x * x * x * x * x * x * y * y * z + 137.8125 * x * x * x * x * x * x * x * z * z * z + 54.140625 * x * x * x * x * x * y * y * y * y * z - 59.0625 * x * x * x * x * x * y * y * z * z * z - 433.125 * x * x * x * x * x * z * z * z * z * z - 9.84375 * x * x * x * y * y * y * y * y * y * z + 98.4375 * x * x * x * y * y * y * y * z * z * z + 551.25 * x * x * x * y * y * z * z * z * z * z + 315.0 * x * x * x * z * z * z * z * z * z * z - 27.0703125 * x * y * y * y * y * y * y * y * y * z + 295.3125 * x * y * y * y * y * y * y * z * z * z - 905.625 * x * y * y * y * y * z * z * z * z * z + 315.0 * x * y * y * z * z * z * z * z * z * z) + e_2 * (88.59375 * x * x * x * x * x * x * x * z + 383.90625 * x * x * x * x * x * y * y * z - 748.125 * x * x * x * x * x * z * z * z + 344.53125 * x * x * x * y * y * y * y * z + 2756.25 * x * x * x * y * y * z * z * z + 1417.5 * x * x * x * z * z * z * z * z + 49.21875 * x * y * y * y * y * y * y * z - 2165.625 * x * y * y * y * y * z * z * z + 1417.5 * x * y * y * z * z * z * z * z + 630.0 * x * z * z * z * z * z * z * z) + e_3 * (4725.0 * x * x * x * y * y * z + 3150.0 * x * x * x * z * z * z - 1575.0 * x * y * y * y * y * z + 3150.0 * x * y * y * z * z * z + 6300.0 * x * z * z * z * z * z) + e_4 * (3543.75 * x * x * x * z + 3543.75 * x * y * y * z + 18900.0 * x * z * z * z) + e_5 * (14175.0 * x * z);
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

        pc_38[k] = e_0 * (-std::sqrt(545.0592041015625) * x * x * x * x * x * x * x * x * y * y * z * z - std::sqrt(3875.9765625) * x * x * x * x * x * x * y * y * y * y * z * z + std::sqrt(24224.853515625) * x * x * x * x * x * x * y * y * z * z * z * z - std::sqrt(2180.23681640625) * x * x * x * x * y * y * y * y * y * y * z * z + std::sqrt(67291.259765625) * x * x * x * x * y * y * y * y * z * z * z * z - std::sqrt(82015.6640625) * x * x * x * x * y * y * z * z * z * z * z * z + std::sqrt(2691.650390625) * x * x * y * y * y * y * y * y * z * z * z * z - std::sqrt(36451.40625) * x * x * y * y * y * y * z * z * z * z * z * z + std::sqrt(9922.5) * x * x * y * y * z * z * z * z * z * z * z * z + std::sqrt(60.5621337890625) * y * y * y * y * y * y * y * y * y * y * z * z - std::sqrt(2691.650390625) * y * y * y * y * y * y * y * y * z * z * z * z + std::sqrt(9112.8515625) * y * y * y * y * y * y * z * z * z * z * z * z - std::sqrt(1102.5) * y * y * y * y * z * z * z * z * z * z * z * z) + e_1 * (-std::sqrt(545.0592041015625) * x * x * x * x * x * x * x * x * y * y - std::sqrt(545.0592041015625) * x * x * x * x * x * x * x * x * z * z - std::sqrt(3875.9765625) * x * x * x * x * x * x * y * y * y * y - std::sqrt(8720.947265625) * x * x * x * x * x * x * y * y * z * z + std::sqrt(24224.853515625) * x * x * x * x * x * x * z * z * z * z - std::sqrt(2180.23681640625) * x * x * x * x * y * y * y * y * y * y - std::sqrt(6056.21337890625) * x * x * x * x * y * y * y * y * z * z - std::sqrt(163760.009765625) * x * x * x * x * y * y * z * z * z * z - std::sqrt(82015.6640625) * x * x * x * x * z * z * z * z * z * z + std::sqrt(968.994140625) * x * x * y * y * y * y * y * y * z * z - std::sqrt(280039.306640625) * x * x * y * y * y * y * z * z * z * z - std::sqrt(5581.40625) * x * x * y * y * z * z * z * z * z * z + std::sqrt(9922.5) * x * x * z * z * z * z * z * z * z * z + std::sqrt(60.5621337890625) * y * y * y * y * y * y * y * y * y * y + std::sqrt(1514.0533447265625) * y * y * y * y * y * y * y * y * z * z + std::sqrt(968.994140625) * y * y * y * y * y * y * z * z * z * z + std::sqrt(96899.4140625) * y * y * y * y * z * z * z * z * z * z - std::sqrt(9922.5) * y * y * z * z * z * z * z * z * z * z) + e_2 * (-std::sqrt(545.0592041015625) * x * x * x * x * x * x * x * x - std::sqrt(313954.1015625) * x * x * x * x * x * x * y * y + std::sqrt(8720.947265625) * x * x * x * x * x * x * z * z - std::sqrt(732801.8188476562) * x * x * x * x * y * y * y * y - std::sqrt(4613381.103515625) * x * x * x * x * y * y * z * z - std::sqrt(1399227.5390625) * x * x * x * x * z * z * z * z - std::sqrt(15503.90625) * x * x * y * y * y * y * y * y - std::sqrt(2520353.759765625) * x * x * y * y * y * y * z * z - std::sqrt(11302347.65625) * x * x * y * y * z * z * z * z + std::sqrt(248062.5) * x * x * z * z * z * z * z * z + std::sqrt(37851.33361816406) * y * y * y * y * y * y * y * y + std::sqrt(427326.416015625) * y * y * y * y * y * y * z * z + std::sqrt(5306211.9140625) * y * y * y * y * z * z * z * z - std::sqrt(248062.5) * y * y * z * z * z * z * z * z) + e_3 * (-std::sqrt(139535.15625) * x * x * x * x * x * x - std::sqrt(27348890.625) * x * x * x * x * y * y - std::sqrt(6837222.65625) * x * x * x * x * z * z - std::sqrt(8201566.40625) * x * x * y * y * y * y - std::sqrt(201488765.625) * x * x * y * y * z * z - std::sqrt(992250.0) * x * x * z * z * z * z + std::sqrt(3969000.0) * y * y * y * y * y * y + std::sqrt(53969097.65625) * y * y * y * y * z * z + std::sqrt(992250.0) * y * y * z * z * z * z) + e_4 * (-std::sqrt(11302347.65625) * x * x * x * x - std::sqrt(246140015.625) * x * x * y * y - std::sqrt(80372250.0) * x * x * z * z + std::sqrt(73814097.65625) * y * y * y * y + std::sqrt(80372250.0) * y * y * z * z) + e_5 * (-std::sqrt(80372250.0) * x * x + std::sqrt(80372250.0) * y * y);
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

        pc_39[k] = e_0 * (std::sqrt(6.488800048828125) * x * x * x * x * x * x * x * x * x * x * y * z + std::sqrt(87.23831176757812) * x * x * x * x * x * x * x * x * y * y * y * z - std::sqrt(2771.4385986328125) * x * x * x * x * x * x * x * x * y * z * z * z + std::sqrt(141.3116455078125) * x * x * x * x * x * x * y * y * y * y * y * z - std::sqrt(19708.0078125) * x * x * x * x * x * x * y * y * y * z * z * z + std::sqrt(33637.939453125) * x * x * x * x * x * x * y * z * z * z * z * z + std::sqrt(25.9552001953125) * x * x * x * x * y * y * y * y * y * y * y * z - std::sqrt(11085.75439453125) * x * x * x * x * y * y * y * y * y * z * z * z + std::sqrt(93438.720703125) * x * x * x * x * y * y * y * z * z * z * z * z - std::sqrt(29302.3828125) * x * x * x * x * y * z * z * z * z * z * z * z - std::sqrt(0.720977783203125) * x * x * y * y * y * y * y * y * y * y * y * z + std::sqrt(3737.548828125) * x * x * y * y * y * y * y * z * z * z * z * z - std::sqrt(13023.28125) * x * x * y * y * y * z * z * z * z * z * z * z + std::sqrt(472.5) * x * x * y * z * z * z * z * z * z * z * z * z - std::sqrt(0.720977783203125) * y * y * y * y * y * y * y * y * y * y * y * z + std::sqrt(307.9376220703125) * y * y * y * y * y * y * y * y * y * z * z * z - std::sqrt(3737.548828125) * y * y * y * y * y * y * y * z * z * z * z * z + std::sqrt(3255.8203125) * y * y * y * y * y * z * z * z * z * z * z * z - std::sqrt(52.5) * y * y * y * z * z * z * z * z * z * z * z * z) + e_1 * (-std::sqrt(233.5968017578125) * x * x * x * x * x * x * x * x * y * z - std::sqrt(1661.1328125) * x * x * x * x * x * x * y * y * y * z - std::sqrt(3737.548828125) * x * x * x * x * x * x * y * z * z * z - std::sqrt(934.38720703125) * x * x * x * x * y * y * y * y * y * z - std::sqrt(10382.080078125) * x * x * x * x * y * y * y * z * z * z + std::sqrt(5382.0703125) * x * x * x * x * y * z * z * z * z * z - std::sqrt(415.283203125) * x * x * y * y * y * y * y * z * z * z + std::sqrt(2392.03125) * x * x * y * y * y * z * z * z * z * z - std::sqrt(344452.5) * x * x * y * z * z * z * z * z * z * z + std::sqrt(25.9552001953125) * y * y * y * y * y * y * y * y * y * z + std::sqrt(415.283203125) * y * y * y * y * y * y * y * z * z * z - std::sqrt(598.0078125) * y * y * y * y * y * z * z * z * z * z + std::sqrt(38272.5) * y * y * y * z * z * z * z * z * z * z) + e_2 * (-std::sqrt(134551.7578125) * x * x * x * x * x * x * y * z - std::sqrt(373754.8828125) * x * x * x * x * y * y * y * z - std::sqrt(59800.78125) * x * x * x * x * y * z * z * z - std::sqrt(14950.1953125) * x * x * y * y * y * y * y * z - std::sqrt(26578.125) * x * x * y * y * y * z * z * z - std::sqrt(34445250.0) * x * x * y * z * z * z * z * z + std::sqrt(14950.1953125) * y * y * y * y * y * y * y * z + std::sqrt(6644.53125) * y * y * y * y * y * z * z * z + std::sqrt(3827250.0) * y * y * y * z * z * z * z * z) + e_3 * (-std::sqrt(7235894.53125) * x * x * x * x * y * z - std::sqrt(3215953.125) * x * x * y * y * y * z - std::sqrt(408665250.0) * x * x * y * z * z * z + std::sqrt(803988.28125) * y * y * y * y * y * z + std::sqrt(45407250.0) * y * y * y * z * z * z) + e_4 * (-std::sqrt(421954312.5) * x * x * y * z + std::sqrt(46883812.5) * y * y * y * z);
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

        pc_40[k] = e_0 * (-std::sqrt(545.0592041015625) * x * x * x * x * x * x * x * x * x * y * z * z - std::sqrt(3875.9765625) * x * x * x * x * x * x * x * y * y * y * z * z + std::sqrt(24224.853515625) * x * x * x * x * x * x * x * y * z * z * z * z - std::sqrt(2180.23681640625) * x * x * x * x * x * y * y * y * y * y * z * z + std::sqrt(67291.259765625) * x * x * x * x * x * y * y * y * z * z * z * z - std::sqrt(82015.6640625) * x * x * x * x * x * y * z * z * z * z * z * z + std::sqrt(2691.650390625) * x * x * x * y * y * y * y * y * z * z * z * z - std::sqrt(36451.40625) * x * x * x * y * y * y * z * z * z * z * z * z + std::sqrt(9922.5) * x * x * x * y * z * z * z * z * z * z * z * z + std::sqrt(60.5621337890625) * x * y * y * y * y * y * y * y * y * y * z * z - std::sqrt(2691.650390625) * x * y * y * y * y * y * y * y * z * z * z * z + std::sqrt(9112.8515625) * x * y * y * y * y * y * z * z * z * z * z * z - std::sqrt(1102.5) * x * y * y * y * z * z * z * z * z * z * z * z) + e_1 * (-std::sqrt(545.0592041015625) * x * x * x * x * x * x * x * x * x * y - std::sqrt(3875.9765625) * x * x * x * x * x * x * x * y * y * y - std::sqrt(8720.947265625) * x * x * x * x * x * x * x * y * z * z - std::sqrt(2180.23681640625) * x * x * x * x * x * y * y * y * y * y - std::sqrt(47480.712890625) * x * x * x * x * x * y * y * y * z * z - std::sqrt(62015.625) * x * x * x * x * x * y * z * z * z * z - std::sqrt(24224.853515625) * x * x * x * y * y * y * y * y * z * z + std::sqrt(62015.625) * x * x * x * y * y * y * z * z * z * z - std::sqrt(419225.625) * x * x * x * y * z * z * z * z * z * z + std::sqrt(60.5621337890625) * x * y * y * y * y * y * y * y * y * y - std::sqrt(968.994140625) * x * y * y * y * y * y * y * y * z * z + std::sqrt(248062.5) * x * y * y * y * y * y * z * z * z * z - std::sqrt(300155.625) * x * y * y * y * z * z * z * z * z * z + std::sqrt(39690.0) * x * y * z * z * z * z * z * z * z * z) + e_2 * (-std::sqrt(313954.1015625) * x * x * x * x * x * x * x * y - std::sqrt(992250.0) * x * x * x * x * x * y * y * y - std::sqrt(4220938.4765625) * x * x * x * x * x * y * z * z - std::sqrt(96899.4140625) * x * x * x * y * y * y * y * y - std::sqrt(1255816.40625) * x * x * x * y * y * y * z * z - std::sqrt(32806265.625) * x * x * x * y * z * z * z * z + std::sqrt(15503.90625) * x * y * y * y * y * y * y * y + std::sqrt(872094.7265625) * x * y * y * y * y * y * z * z - std::sqrt(1550390.625) * x * y * y * y * z * z * z * z + std::sqrt(992250.0) * x * y * z * z * z * z * z * z) + e_3 * (-std::sqrt(31395410.15625) * x * x * x * x * x * y - std::sqrt(22387640.625) * x * x * x * y * y * y - std::sqrt(377303062.5) * x * x * x * y * z * z + std::sqrt(759691.40625) * x * y * y * y * y * y - std::sqrt(248062.5) * x * y * y * y * z * z - std::sqrt(3969000.0) * x * y * z * z * z * z) + e_4 * (-std::sqrt(502326562.5) * x * x * x * y - std::sqrt(2232562.5) * x * y * y * y - std::sqrt(321489000.0) * x * y * z * z) + e_5 * (-std::sqrt(321489000.0) * x * y);
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

        pc_41[k] = e_0 * (-3.69140625 * x * x * x * x * x * x * x * x * x * x * y * z - 6.15234375 * x * x * x * x * x * x * x * x * y * y * y * z + 68.90625 * x * x * x * x * x * x * x * x * y * z * z * z + 2.4609375 * x * x * x * x * x * x * y * y * y * y * y * z + 45.9375 * x * x * x * x * x * x * y * y * y * z * z * z - 216.5625 * x * x * x * x * x * x * y * z * z * z * z * z + 7.3828125 * x * x * x * x * y * y * y * y * y * y * y * z - 91.875 * x * x * x * x * y * y * y * y * y * z * z * z + 72.1875 * x * x * x * x * y * y * y * z * z * z * z * z + 157.5 * x * x * x * x * y * z * z * z * z * z * z * z + 1.23046875 * x * x * y * y * y * y * y * y * y * y * y * z - 45.9375 * x * x * y * y * y * y * y * y * y * z * z * z + 216.5625 * x * x * y * y * y * y * y * z * z * z * z * z - 210.0 * x * x * y * y * y * z * z * z * z * z * z * z - 1.23046875 * y * y * y * y * y * y * y * y * y * y * y * z + 22.96875 * y * y * y * y * y * y * y * y * y * z * z * z - 72.1875 * y * y * y * y * y * y * y * z * z * z * z * z + 52.5 * y * y * y * y * y * z * z * z * z * z * z * z) + e_1 * (22.1484375 * x * x * x * x * x * x * x * x * y * z - 9.84375 * x * x * x * x * x * x * y * y * y * z - 98.4375 * x * x * x * x * x * x * y * z * z * z - 83.671875 * x * x * x * x * y * y * y * y * y * z + 492.1875 * x * x * x * x * y * y * y * z * z * z + 275.625 * x * x * x * x * y * z * z * z * z * z - 49.21875 * x * x * y * y * y * y * y * y * y * z + 649.6875 * x * x * y * y * y * y * y * z * z * z - 1811.25 * x * x * y * y * y * z * z * z * z * z + 315.0 * x * x * y * z * z * z * z * z * z * z + 2.4609375 * y * y * y * y * y * y * y * y * y * z + 59.0625 * y * y * y * y * y * y * y * z * z * z - 196.875 * y * y * y * y * y * z * z * z * z * z + 315.0 * y * y * y * z * z * z * z * z * z * z) + e_2 * (147.65625 * x * x * x * x * x * x * y * z + 246.09375 * x * x * x * x * y * y * y * z + 1378.125 * x * x * x * x * y * z * z * z + 206.71875 * x * x * y * y * y * y * y * z - 4331.25 * x * x * y * y * y * z * z * z + 1417.5 * x * x * y * z * z * z * z * z + 108.28125 * y * y * y * y * y * y * y * z - 39.375 * y * y * y * y * y * z * z * z + 1417.5 * y * y * y * z * z * z * z * z + 630.0 * y * z * z * z * z * z * z * z) + e_3 * (2362.5 * x * x * x * x * y * z - 3150.0 * x * x * y * y * y * z + 3150.0 * x * x * y * z * z * z + 787.5 * y * y * y * y * y * z + 3150.0 * y * y * y * z * z * z + 6300.0 * y * z * z * z * z * z) + e_4 * (3543.75 * x * x * y * z + 3543.75 * y * y * y * z + 18900.0 * y * z * z * z) + e_5 * (14175.0 * y * z);
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

        pc_42[k] = e_0 * (22.1484375 * x * x * x * x * x * x * x * x * x * y * z * z - 29.53125 * x * x * x * x * x * x * x * y * y * y * z * z - 118.125 * x * x * x * x * x * x * x * y * z * z * z * z - 103.359375 * x * x * x * x * x * y * y * y * y * y * z * z + 275.625 * x * x * x * x * x * y * y * y * z * z * z * z + 157.5 * x * x * x * x * x * y * z * z * z * z * z * z - 29.53125 * x * x * x * y * y * y * y * y * y * y * z * z + 275.625 * x * x * x * y * y * y * y * y * z * z * z * z - 525.0 * x * x * x * y * y * y * z * z * z * z * z * z + 22.1484375 * x * y * y * y * y * y * y * y * y * y * z * z - 118.125 * x * y * y * y * y * y * y * y * z * z * z * z + 157.5 * x * y * y * y * y * y * z * z * z * z * z * z) + e_1 * (22.1484375 * x * x * x * x * x * x * x * x * x * y - 29.53125 * x * x * x * x * x * x * x * y * y * y - 103.359375 * x * x * x * x * x * y * y * y * y * y + 708.75 * x * x * x * x * x * y * z * z * z * z - 29.53125 * x * x * x * y * y * y * y * y * y * y - 2362.5 * x * x * x * y * y * y * z * z * z * z + 22.1484375 * x * y * y * y * y * y * y * y * y * y + 708.75 * x * y * y * y * y * y * z * z * z * z) + e_2 * (354.375 * x * x * x * x * x * x * x * y - 826.875 * x * x * x * x * x * y * y * y + 2126.25 * x * x * x * x * x * y * z * z - 826.875 * x * x * x * y * y * y * y * y - 7087.5 * x * x * x * y * y * y * z * z + 354.375 * x * y * y * y * y * y * y * y + 2126.25 * x * y * y * y * y * y * z * z) + e_3 * (2362.5 * x * x * x * x * x * y - 7875.0 * x * x * x * y * y * y + 2362.5 * x * y * y * y * y * y);
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

        pc_43[k] = e_0 * (std::sqrt(16.351776123046875) * x * x * x * x * x * x * x * x * x * x * y * z - std::sqrt(307.0500183105469) * x * x * x * x * x * x * x * x * y * y * y * z - std::sqrt(2623.5516357421875) * x * x * x * x * x * x * x * x * y * z * z * z - std::sqrt(1228.2000732421875) * x * x * x * x * x * x * y * y * y * y * y * z + std::sqrt(74625.46875) * x * x * x * x * x * x * y * y * y * z * z * z + std::sqrt(11627.9296875) * x * x * x * x * x * x * y * z * z * z * z * z - std::sqrt(7.2674560546875) * x * x * x * x * y * y * y * y * y * y * y * z + std::sqrt(29150.57373046875) * x * x * x * x * y * y * y * y * y * z * z * z - std::sqrt(466409.1796875) * x * x * x * x * y * y * y * z * z * z * z * z + std::sqrt(89.02633666992188) * x * x * y * y * y * y * y * y * y * y * y * z - std::sqrt(18656.3671875) * x * x * y * y * y * y * y * y * y * z * z * z + std::sqrt(104651.3671875) * x * x * y * y * y * y * y * z * z * z * z * z - std::sqrt(1.816864013671875) * y * y * y * y * y * y * y * y * y * y * y * z + std::sqrt(291.5057373046875) * y * y * y * y * y * y * y * y * y * z * z * z - std::sqrt(1291.9921875) * y * y * y * y * y * y * y * z * z * z * z * z) + e_1 * (-std::sqrt(588.6639404296875) * x * x * x * x * x * x * x * x * y * z - std::sqrt(465.1171875) * x * x * x * x * x * x * y * y * y * z + std::sqrt(215000.419921875) * x * x * x * x * x * x * y * z * z * z - std::sqrt(4912.80029296875) * x * x * x * x * y * y * y * y * y * z - std::sqrt(1049420.654296875) * x * x * x * x * y * y * y * z * z * z - std::sqrt(186046.875) * x * x * x * x * y * z * z * z * z * z - std::sqrt(7441.875) * x * x * y * y * y * y * y * y * y * z + std::sqrt(1935003.779296875) * x * x * y * y * y * y * y * z * z * z - std::sqrt(744187.5) * x * x * y * y * y * z * z * z * z * z - std::sqrt(181.6864013671875) * y * y * y * y * y * y * y * y * y * z + std::sqrt(14069.794921875) * y * y * y * y * y * y * y * z * z * z - std::sqrt(186046.875) * y * y * y * y * y * z * z * z * z * z) + e_2 * (std::sqrt(104651.3671875) * x * x * x * x * x * x * y * z - std::sqrt(4197682.6171875) * x * x * x * x * y * y * y * z - std::sqrt(46511.71875) * x * x * x * x * y * z * z * z + std::sqrt(941862.3046875) * x * x * y * y * y * y * y * z - std::sqrt(186046.875) * x * x * y * y * y * z * z * z - std::sqrt(6697687.5) * x * x * y * z * z * z * z * z - std::sqrt(11627.9296875) * y * y * y * y * y * y * y * z - std::sqrt(46511.71875) * y * y * y * y * y * z * z * z - std::sqrt(6697687.5) * y * y * y * z * z * z * z * z) + e_3 * (-std::sqrt(418605.46875) * x * x * x * x * y * z - std::sqrt(1674421.875) * x * x * y * y * y * z - std::sqrt(90046687.5) * x * x * y * z * z * z - std::sqrt(418605.46875) * y * y * y * y * y * z - std::sqrt(90046687.5) * y * y * y * z * z * z - std::sqrt(11907000.0) * y * z * z * z * z * z) + e_4 * (-std::sqrt(82046671.875) * x * x * y * z - std::sqrt(82046671.875) * y * y * y * z - std::sqrt(328186687.5) * y * z * z * z) + e_5 * (-std::sqrt(328186687.5) * y * z);
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

        pc_44[k] = e_0 * (-std::sqrt(359.73907470703125) * x * x * x * x * x * x * x * x * x * y * z * z + std::sqrt(31337.2705078125) * x * x * x * x * x * x * x * y * y * y * z * z + std::sqrt(2558.14453125) * x * x * x * x * x * x * x * y * z * z * z * z + std::sqrt(1438.956298828125) * x * x * x * x * x * y * y * y * y * y * z * z - std::sqrt(273152.98828125) * x * x * x * x * x * y * y * y * z * z * z * z - std::sqrt(15988.4033203125) * x * x * x * y * y * y * y * y * y * y * z * z + std::sqrt(177648.92578125) * x * x * x * y * y * y * y * y * z * z * z * z + std::sqrt(999.2752075195312) * x * y * y * y * y * y * y * y * y * y * z * z - std::sqrt(7105.95703125) * x * y * y * y * y * y * y * y * z * z * z * z) + e_1 * (-std::sqrt(359.73907470703125) * x * x * x * x * x * x * x * x * x * y + std::sqrt(31337.2705078125) * x * x * x * x * x * x * x * y * y * y + std::sqrt(23023.30078125) * x * x * x * x * x * x * x * y * z * z + std::sqrt(1438.956298828125) * x * x * x * x * x * y * y * y * y * y + std::sqrt(923490.17578125) * x * x * x * x * x * y * y * y * z * z - std::sqrt(255814.453125) * x * x * x * x * x * y * z * z * z * z - std::sqrt(15988.4033203125) * x * x * x * y * y * y * y * y * y * y + std::sqrt(63953.61328125) * x * x * x * y * y * y * y * y * z * z - std::sqrt(1023257.8125) * x * x * x * y * y * y * z * z * z * z + std::sqrt(999.2752075195312) * x * y * y * y * y * y * y * y * y * y + std::sqrt(63953.61328125) * x * y * y * y * y * y * y * y * z * z - std::sqrt(255814.453125) * x * y * y * y * y * y * z * z * z * z) + e_2 * (std::sqrt(6395361.328125) * x * x * x * x * x * y * y * y + std::sqrt(2302330.078125) * x * x * x * x * x * y * z * z - std::sqrt(1023257.8125) * x * x * x * y * y * y * y * y + std::sqrt(9209320.3125) * x * x * x * y * y * y * z * z - std::sqrt(16372125.0) * x * x * x * y * z * z * z * z + std::sqrt(255814.453125) * x * y * y * y * y * y * y * y + std::sqrt(2302330.078125) * x * y * y * y * y * y * z * z - std::sqrt(16372125.0) * x * y * y * y * z * z * z * z) + e_3 * (std::sqrt(9209320.3125) * x * x * x * x * x * y + std::sqrt(36837281.25) * x * x * x * y * y * y + std::sqrt(9209320.3125) * x * y * y * y * y * y - std::sqrt(65488500.0) * x * y * z * z * z * z) + e_4 * (std::sqrt(147349125.0) * x * x * x * y + std::sqrt(147349125.0) * x * y * y * y - std::sqrt(147349125.0) * x * y * z * z) + e_5 * (std::sqrt(147349125.0) * x * y);
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

        pc_45[k] = e_0 * (-std::sqrt(29.978256225585938) * x * x * x * x * x * x * x * x * x * x * y * z + std::sqrt(6158.866195678711) * x * x * x * x * x * x * x * x * y * y * y * z + std::sqrt(213.1787109375) * x * x * x * x * x * x * x * x * y * z * z * z - std::sqrt(652.8598022460938) * x * x * x * x * x * x * y * y * y * y * y * z - std::sqrt(50120.68359375) * x * x * x * x * x * x * y * y * y * z * z * z - std::sqrt(5875.738220214844) * x * x * x * x * y * y * y * y * y * y * y * z + std::sqrt(85271.484375) * x * x * x * x * y * y * y * y * y * z * z * z + std::sqrt(962.6351165771484) * x * x * y * y * y * y * y * y * y * y * y * z - std::sqrt(7674.43359375) * x * x * y * y * y * y * y * y * y * z * z * z - std::sqrt(3.3309173583984375) * y * y * y * y * y * y * y * y * y * y * y * z + std::sqrt(23.6865234375) * y * y * y * y * y * y * y * y * y * z * z * z) + e_1 * (std::sqrt(1079.2172241210938) * x * x * x * x * x * x * x * x * y * z + std::sqrt(1613549.6630859375) * x * x * x * x * x * x * y * y * y * z - std::sqrt(69069.90234375) * x * x * x * x * x * x * y * z * z * z - std::sqrt(1247575.1110839844) * x * x * x * x * y * y * y * y * y * z - std::sqrt(191860.83984375) * x * x * x * x * y * y * y * z * z * z + std::sqrt(155407.2802734375) * x * x * y * y * y * y * y * y * y * z - std::sqrt(7674.43359375) * x * x * y * y * y * y * y * z * z * z - std::sqrt(2997.8256225585938) * y * y * y * y * y * y * y * y * y * z + std::sqrt(7674.43359375) * y * y * y * y * y * y * y * z * z * z) + e_2 * (std::sqrt(3885182.0068359375) * x * x * x * x * x * x * y * z + std::sqrt(10792172.241210938) * x * x * x * x * y * y * y * z - std::sqrt(6906990.234375) * x * x * x * x * y * z * z * z + std::sqrt(431686.8896484375) * x * x * y * y * y * y * y * z - std::sqrt(3069773.4375) * x * x * y * y * y * z * z * z - std::sqrt(431686.8896484375) * y * y * y * y * y * y * y * z + std::sqrt(767443.359375) * y * y * y * y * y * z * z * z) + e_3 * (std::sqrt(110511843.75) * x * x * x * x * y * z + std::sqrt(49116375.0) * x * x * y * y * y * z - std::sqrt(49116375.0) * x * x * y * z * z * z - std::sqrt(12279093.75) * y * y * y * y * y * z + std::sqrt(5457375.0) * y * y * y * z * z * z) + e_4 * (std::sqrt(248651648.4375) * x * x * y * z - std::sqrt(27627960.9375) * y * y * y * z);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ab_x, ab_y, ab_z : simd::cache_line_size())
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
        const auto e_6 = pe_6[k];

        pc_46[k] = e_0 * (0.8203125 * x * x * x * x * x * x * x * x * x * x * y * y + 3.28125 * x * x * x * x * x * x * x * x * y * y * y * y - 26.25 * x * x * x * x * x * x * x * x * y * y * z * z + 4.921875 * x * x * x * x * x * x * y * y * y * y * y * y - 78.75 * x * x * x * x * x * x * y * y * y * y * z * z + 236.25 * x * x * x * x * x * x * y * y * z * z * z * z + 3.28125 * x * x * x * x * y * y * y * y * y * y * y * y - 78.75 * x * x * x * x * y * y * y * y * y * y * z * z + 472.5 * x * x * x * x * y * y * y * y * z * z * z * z - 420.0 * x * x * x * x * y * y * z * z * z * z * z * z + 0.8203125 * x * x * y * y * y * y * y * y * y * y * y * y - 26.25 * x * x * y * y * y * y * y * y * y * y * z * z + 236.25 * x * x * y * y * y * y * y * y * z * z * z * z - 420.0 * x * x * y * y * y * y * z * z * z * z * z * z + 210.0 * x * x * y * y * z * z * z * z * z * z * z * z) + e_1 * (0.8203125 * x * x * x * x * x * x * x * x * x * x + 30.3515625 * x * x * x * x * x * x * x * x * y * y - 26.25 * x * x * x * x * x * x * x * x * z * z + 86.953125 * x * x * x * x * x * x * y * y * y * y + 210.0 * x * x * x * x * x * x * y * y * z * z + 236.25 * x * x * x * x * x * x * z * z * z * z + 86.953125 * x * x * x * x * y * y * y * y * y * y + 472.5 * x * x * x * x * y * y * y * y * z * z + 78.75 * x * x * x * x * y * y * z * z * z * z - 420.0 * x * x * x * x * z * z * z * z * z * z + 30.3515625 * x * x * y * y * y * y * y * y * y * y + 210.0 * x * x * y * y * y * y * y * y * z * z + 78.75 * x * x * y * y * y * y * z * z * z * z + 840.0 * x * x * y * y * z * z * z * z * z * z + 210.0 * x * x * z * z * z * z * z * z * z * z + 0.8203125 * y * y * y * y * y * y * y * y * y * y - 26.25 * y * y * y * y * y * y * y * y * z * z + 236.25 * y * y * y * y * y * y * z * z * z * z - 420.0 * y * y * y * y * z * z * z * z * z * z + 210.0 * y * y * z * z * z * z * z * z * z * z) + e_2 * (20.5078125 * x * x * x * x * x * x * x * x + 790.78125 * x * x * x * x * x * x * y * y + 446.25 * x * x * x * x * x * x * z * z + 1540.546875 * x * x * x * x * y * y * y * y + 3228.75 * x * x * x * x * y * y * z * z - 1338.75 * x * x * x * x * z * z * z * z + 790.78125 * x * x * y * y * y * y * y * y + 3228.75 * x * x * y * y * y * y * z * z + 6772.5 * x * x * y * y * z * z * z * z + 2100.0 * x * x * z * z * z * z * z * z + 20.5078125 * y * y * y * y * y * y * y * y + 446.25 * y * y * y * y * y * y * z * z - 1338.75 * y * y * y * y * z * z * z * z + 2100.0 * y * y * z * z * z * z * z * z + 210.0 * z * z * z * z * z * z * z * z) + e_3 * (603.75 * x * x * x * x * x * x + 8111.25 * x * x * x * x * y * y + 630.0 * x * x * x * x * z * z + 8111.25 * x * x * y * y * y * y + 26460.0 * x * x * y * y * z * z + 10080.0 * x * x * z * z * z * z + 603.75 * y * y * y * y * y * y + 630.0 * y * y * y * y * z * z + 10080.0 * y * y * z * z * z * z + 3360.0 * z * z * z * z * z * z) + e_4 * (4449.375 * x * x * x * x + 30948.75 * x * x * y * y + 22680.0 * x * x * z * z + 4449.375 * y * y * y * y + 22680.0 * y * y * z * z + 17640.0 * z * z * z * z) + e_5 * (16065.0 * x * x + 16065.0 * y * y + 30240.0 * z * z) + e_6 * (10395.0);
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

        pc_47[k] = e_0 * (std::sqrt(6.7291259765625) * x * x * x * x * x * x * x * x * x * y * y * z + std::sqrt(107.666015625) * x * x * x * x * x * x * x * y * y * y * y * z - std::sqrt(2691.650390625) * x * x * x * x * x * x * x * y * y * z * z * z + std::sqrt(242.24853515625) * x * x * x * x * x * y * y * y * y * y * y * z - std::sqrt(24224.853515625) * x * x * x * x * x * y * y * y * y * z * z * z + std::sqrt(44806.2890625) * x * x * x * x * x * y * y * z * z * z * z * z + std::sqrt(107.666015625) * x * x * x * y * y * y * y * y * y * y * y * z - std::sqrt(24224.853515625) * x * x * x * y * y * y * y * y * y * z * z * z + std::sqrt(179225.15625) * x * x * x * y * y * y * y * z * z * z * z * z - std::sqrt(54022.5) * x * x * x * y * y * z * z * z * z * z * z * z + std::sqrt(6.7291259765625) * x * y * y * y * y * y * y * y * y * y * y * z - std::sqrt(2691.650390625) * x * y * y * y * y * y * y * y * y * z * z * z + std::sqrt(44806.2890625) * x * y * y * y * y * y * y * z * z * z * z * z - std::sqrt(54022.5) * x * y * y * y * y * z * z * z * z * z * z * z + std::sqrt(4410.0) * x * y * y * z * z * z * z * z * z * z * z * z) + e_1 * (std::sqrt(6.7291259765625) * x * x * x * x * x * x * x * x * x * z - std::sqrt(2691.650390625) * x * x * x * x * x * x * x * z * z * z - std::sqrt(242.24853515625) * x * x * x * x * x * y * y * y * y * z + std::sqrt(8720.947265625) * x * x * x * x * x * y * y * z * z * z + std::sqrt(44806.2890625) * x * x * x * x * x * z * z * z * z * z - std::sqrt(430.6640625) * x * x * x * y * y * y * y * y * y * z + std::sqrt(117248.291015625) * x * x * x * y * y * y * y * z * z * z - std::sqrt(50232.65625) * x * x * x * y * y * z * z * z * z * z - std::sqrt(54022.5) * x * x * x * z * z * z * z * z * z * z - std::sqrt(60.5621337890625) * x * y * y * y * y * y * y * y * y * z + std::sqrt(38867.431640625) * x * y * y * y * y * y * y * z * z * z - std::sqrt(189922.8515625) * x * y * y * y * y * z * z * z * z * z + std::sqrt(89302.5) * x * y * y * z * z * z * z * z * z * z + std::sqrt(4410.0) * x * z * z * z * z * z * z * z * z * z) + e_2 * (-std::sqrt(968.994140625) * x * x * x * x * x * x * x * z + std::sqrt(8720.947265625) * x * x * x * x * x * y * y * z + std::sqrt(313954.1015625) * x * x * x * x * x * z * z * z + std::sqrt(78488.525390625) * x * x * x * y * y * y * y * z + std::sqrt(139535.15625) * x * x * x * y * y * z * z * z - std::sqrt(2232562.5) * x * x * x * z * z * z * z * z + std::sqrt(24224.853515625) * x * y * y * y * y * y * y * z - std::sqrt(34883.7890625) * x * y * y * y * y * z * z * z + std::sqrt(2232562.5) * x * y * y * z * z * z * z * z + std::sqrt(992250.0) * x * z * z * z * z * z * z * z) + e_3 * (std::sqrt(139535.15625) * x * x * x * x * x * z + std::sqrt(1550390.625) * x * x * x * y * y * z - std::sqrt(8930250.0) * x * x * x * z * z * z + std::sqrt(759691.40625) * x * y * y * y * y * z + std::sqrt(24806250.0) * x * y * y * z * z * z + std::sqrt(35721000.0) * x * z * z * z * z * z) + e_4 * (-std::sqrt(992250.0) * x * x * x * z + std::sqrt(35721000.0) * x * y * y * z + std::sqrt(194481000.0) * x * z * z * z) + e_5 * (std::sqrt(80372250.0) * x * z);
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

        pc_48[k] = e_0 * (-std::sqrt(0.080108642578125) * x * x * x * x * x * x * x * x * x * x * x * y - std::sqrt(2.002716064453125) * x * x * x * x * x * x * x * x * x * y * y * y + std::sqrt(92.6055908203125) * x * x * x * x * x * x * x * x * x * y * z * z - std::sqrt(8.0108642578125) * x * x * x * x * x * x * x * y * y * y * y * y + std::sqrt(1481.689453125) * x * x * x * x * x * x * x * y * y * y * z * z - std::sqrt(8618.408203125) * x * x * x * x * x * x * x * y * z * z * z * z - std::sqrt(8.0108642578125) * x * x * x * x * x * y * y * y * y * y * y * y + std::sqrt(3333.80126953125) * x * x * x * x * x * y * y * y * y * y * z * z - std::sqrt(77565.673828125) * x * x * x * x * x * y * y * y * z * z * z * z + std::sqrt(36521.1328125) * x * x * x * x * x * y * z * z * z * z * z * z - std::sqrt(2.002716064453125) * x * x * x * y * y * y * y * y * y * y * y * y + std::sqrt(1481.689453125) * x * x * x * y * y * y * y * y * y * y * z * z - std::sqrt(77565.673828125) * x * x * x * y * y * y * y * y * z * z * z * z + std::sqrt(146084.53125) * x * x * x * y * y * y * z * z * z * z * z * z - std::sqrt(15172.5) * x * x * x * y * z * z * z * z * z * z * z * z - std::sqrt(0.080108642578125) * x * y * y * y * y * y * y * y * y * y * y * y + std::sqrt(92.6055908203125) * x * y * y * y * y * y * y * y * y * y * z * z - std::sqrt(8618.408203125) * x * y * y * y * y * y * y * y * z * z * z * z + std::sqrt(36521.1328125) * x * y * y * y * y * y * z * z * z * z * z * z - std::sqrt(15172.5) * x * y * y * y * z * z * z * z * z * z * z * z + std::sqrt(210.0) * x * y * z * z * z * z * z * z * z * z * z * z) + e_1 * (-std::sqrt(103.82080078125) * x * x * x * x * x * x * x * x * x * y - std::sqrt(1661.1328125) * x * x * x * x * x * x * x * y * y * y - std::sqrt(9043.9453125) * x * x * x * x * x * x * x * y * z * z - std::sqrt(3737.548828125) * x * x * x * x * x * y * y * y * y * y - std::sqrt(81395.5078125) * x * x * x * x * x * y * y * y * z * z + std::sqrt(6644.53125) * x * x * x * x * x * y * z * z * z * z - std::sqrt(1661.1328125) * x * x * x * y * y * y * y * y * y * y - std::sqrt(81395.5078125) * x * x * x * y * y * y * y * y * z * z + std::sqrt(26578.125) * x * x * x * y * y * y * z * z * z * z - std::sqrt(153090.0) * x * x * x * y * z * z * z * z * z * z - std::sqrt(103.82080078125) * x * y * y * y * y * y * y * y * y * y - std::sqrt(9043.9453125) * x * y * y * y * y * y * y * y * z * z + std::sqrt(6644.53125) * x * y * y * y * y * y * z * z * z * z - std::sqrt(153090.0) * x * y * y * y * z * z * z * z * z * z - std::sqrt(7560.0) * x * y * z * z * z * z * z * z * z * z) + e_2 * (-std::sqrt(85317.626953125) * x * x * x * x * x * x * x * y - std::sqrt(767858.642578125) * x * x * x * x * x * y * y * y - std::sqrt(1397012.6953125) * x * x * x * x * x * y * z * z - std::sqrt(767858.642578125) * x * x * x * y * y * y * y * y - std::sqrt(5588050.78125) * x * x * x * y * y * y * z * z - std::sqrt(5209312.5) * x * x * x * y * z * z * z * z - std::sqrt(85317.626953125) * x * y * y * y * y * y * y * y - std::sqrt(1397012.6953125) * x * y * y * y * y * y * z * z - std::sqrt(5209312.5) * x * y * y * y * z * z * z * z - std::sqrt(5717250.0) * x * y * z * z * z * z * z * z) + e_3 * (-std::sqrt(10988578.125) * x * x * x * x * x * y - std::sqrt(43954312.5) * x * x * x * y * y * y - std::sqrt(118125000.0) * x * x * x * y * z * z - std::sqrt(10988578.125) * x * y * y * y * y * y - std::sqrt(118125000.0) * x * y * y * y * z * z - std::sqrt(272916000.0) * x * y * z * z * z * z) + e_4 * (-std::sqrt(255256312.5) * x * x * x * y - std::sqrt(255256312.5) * x * y * y * y - std::sqrt(1687817250.0) * x * y * z * z) + e_5 * (-std::sqrt(750141000.0) * x * y);
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

        pc_49[k] = e_0 * (std::sqrt(6.7291259765625) * x * x * x * x * x * x * x * x * x * x * y * z + std::sqrt(107.666015625) * x * x * x * x * x * x * x * x * y * y * y * z - std::sqrt(2691.650390625) * x * x * x * x * x * x * x * x * y * z * z * z + std::sqrt(242.24853515625) * x * x * x * x * x * x * y * y * y * y * y * z - std::sqrt(24224.853515625) * x * x * x * x * x * x * y * y * y * z * z * z + std::sqrt(44806.2890625) * x * x * x * x * x * x * y * z * z * z * z * z + std::sqrt(107.666015625) * x * x * x * x * y * y * y * y * y * y * y * z - std::sqrt(24224.853515625) * x * x * x * x * y * y * y * y * y * z * z * z + std::sqrt(179225.15625) * x * x * x * x * y * y * y * z * z * z * z * z - std::sqrt(54022.5) * x * x * x * x * y * z * z * z * z * z * z * z + std::sqrt(6.7291259765625) * x * x * y * y * y * y * y * y * y * y * y * z - std::sqrt(2691.650390625) * x * x * y * y * y * y * y * y * y * z * z * z + std::sqrt(44806.2890625) * x * x * y * y * y * y * y * z * z * z * z * z - std::sqrt(54022.5) * x * x * y * y * y * z * z * z * z * z * z * z + std::sqrt(4410.0) * x * x * y * z * z * z * z * z * z * z * z * z) + e_1 * (-std::sqrt(60.5621337890625) * x * x * x * x * x * x * x * x * y * z - std::sqrt(430.6640625) * x * x * x * x * x * x * y * y * y * z + std::sqrt(38867.431640625) * x * x * x * x * x * x * y * z * z * z - std::sqrt(242.24853515625) * x * x * x * x * y * y * y * y * y * z + std::sqrt(117248.291015625) * x * x * x * x * y * y * y * z * z * z - std::sqrt(189922.8515625) * x * x * x * x * y * z * z * z * z * z + std::sqrt(8720.947265625) * x * x * y * y * y * y * y * z * z * z - std::sqrt(50232.65625) * x * x * y * y * y * z * z * z * z * z + std::sqrt(89302.5) * x * x * y * z * z * z * z * z * z * z + std::sqrt(6.7291259765625) * y * y * y * y * y * y * y * y * y * z - std::sqrt(2691.650390625) * y * y * y * y * y * y * y * z * z * z + std::sqrt(44806.2890625) * y * y * y * y * y * z * z * z * z * z - std::sqrt(54022.5) * y * y * y * z * z * z * z * z * z * z + std::sqrt(4410.0) * y * z * z * z * z * z * z * z * z * z) + e_2 * (std::sqrt(24224.853515625) * x * x * x * x * x * x * y * z + std::sqrt(78488.525390625) * x * x * x * x * y * y * y * z - std::sqrt(34883.7890625) * x * x * x * x * y * z * z * z + std::sqrt(8720.947265625) * x * x * y * y * y * y * y * z + std::sqrt(139535.15625) * x * x * y * y * y * z * z * z + std::sqrt(2232562.5) * x * x * y * z * z * z * z * z - std::sqrt(968.994140625) * y * y * y * y * y * y * y * z + std::sqrt(313954.1015625) * y * y * y * y * y * z * z * z - std::sqrt(2232562.5) * y * y * y * z * z * z * z * z + std::sqrt(992250.0) * y * z * z * z * z * z * z * z) + e_3 * (std::sqrt(759691.40625) * x * x * x * x * y * z + std::sqrt(1550390.625) * x * x * y * y * y * z + std::sqrt(24806250.0) * x * x * y * z * z * z + std::sqrt(139535.15625) * y * y * y * y * y * z - std::sqrt(8930250.0) * y * y * y * z * z * z + std::sqrt(35721000.0) * y * z * z * z * z * z) + e_4 * (std::sqrt(35721000.0) * x * x * y * z - std::sqrt(992250.0) * y * y * y * z + std::sqrt(194481000.0) * y * z * z * z) + e_5 * (std::sqrt(80372250.0) * y * z);
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

        pc_50[k] = e_0 * (0.41015625 * x * x * x * x * x * x * x * x * x * x * x * y + 1.23046875 * x * x * x * x * x * x * x * x * x * y * y * y - 13.125 * x * x * x * x * x * x * x * x * x * y * z * z + 0.8203125 * x * x * x * x * x * x * x * y * y * y * y * y - 26.25 * x * x * x * x * x * x * x * y * y * y * z * z + 118.125 * x * x * x * x * x * x * x * y * z * z * z * z - 0.8203125 * x * x * x * x * x * y * y * y * y * y * y * y + 118.125 * x * x * x * x * x * y * y * y * z * z * z * z - 210.0 * x * x * x * x * x * y * z * z * z * z * z * z - 1.23046875 * x * x * x * y * y * y * y * y * y * y * y * y + 26.25 * x * x * x * y * y * y * y * y * y * y * z * z - 118.125 * x * x * x * y * y * y * y * y * z * z * z * z + 105.0 * x * x * x * y * z * z * z * z * z * z * z * z - 0.41015625 * x * y * y * y * y * y * y * y * y * y * y * y + 13.125 * x * y * y * y * y * y * y * y * y * y * z * z - 118.125 * x * y * y * y * y * y * y * y * z * z * z * z + 210.0 * x * y * y * y * y * y * z * z * z * z * z * z - 105.0 * x * y * y * y * z * z * z * z * z * z * z * z) + e_1 * (13.125 * x * x * x * x * x * x * x * x * x * y + 26.25 * x * x * x * x * x * x * x * y * y * y + 157.5 * x * x * x * x * x * x * x * y * z * z + 157.5 * x * x * x * x * x * y * y * y * z * z - 315.0 * x * x * x * x * x * y * z * z * z * z - 26.25 * x * x * x * y * y * y * y * y * y * y - 157.5 * x * x * x * y * y * y * y * y * z * z + 840.0 * x * x * x * y * z * z * z * z * z * z - 13.125 * x * y * y * y * y * y * y * y * y * y - 157.5 * x * y * y * y * y * y * y * y * z * z + 315.0 * x * y * y * y * y * y * z * z * z * z - 840.0 * x * y * y * y * z * z * z * z * z * z) + e_2 * (354.375 * x * x * x * x * x * x * x * y + 354.375 * x * x * x * x * x * y * y * y + 945.0 * x * x * x * x * x * y * z * z - 354.375 * x * x * x * y * y * y * y * y + 4725.0 * x * x * x * y * z * z * z * z - 354.375 * x * y * y * y * y * y * y * y - 945.0 * x * y * y * y * y * y * z * z - 4725.0 * x * y * y * y * z * z * z * z) + e_3 * (3150.0 * x * x * x * x * x * y + 12600.0 * x * x * x * y * z * z - 3150.0 * x * y * y * y * y * y - 12600.0 * x * y * y * y * z * z) + e_4 * (11025.0 * x * x * x * y - 11025.0 * x * y * y * y);
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

        pc_51[k] = e_0 * (-2.4609375 * x * x * x * x * x * x * x * x * x * x * y * z + 45.9375 * x * x * x * x * x * x * x * x * y * z * z * z + 14.765625 * x * x * x * x * x * x * y * y * y * y * y * z - 45.9375 * x * x * x * x * x * x * y * y * y * z * z * z - 144.375 * x * x * x * x * x * x * y * z * z * z * z * z + 19.6875 * x * x * x * x * y * y * y * y * y * y * y * z - 229.6875 * x * x * x * x * y * y * y * y * y * z * z * z + 288.75 * x * x * x * x * y * y * y * z * z * z * z * z + 105.0 * x * x * x * x * y * z * z * z * z * z * z * z + 7.3828125 * x * x * y * y * y * y * y * y * y * y * y * z - 137.8125 * x * x * y * y * y * y * y * y * y * z * z * z + 433.125 * x * x * y * y * y * y * y * z * z * z * z * z - 315.0 * x * x * y * y * y * z * z * z * z * z * z * z) + e_1 * (27.0703125 * x * x * x * x * x * x * x * x * y * z + 9.84375 * x * x * x * x * x * x * y * y * y * z - 295.3125 * x * x * x * x * x * x * y * z * z * z - 54.140625 * x * x * x * x * y * y * y * y * y * z - 98.4375 * x * x * x * x * y * y * y * z * z * z + 905.625 * x * x * x * x * y * z * z * z * z * z - 29.53125 * x * x * y * y * y * y * y * y * y * z + 59.0625 * x * x * y * y * y * y * y * z * z * z - 551.25 * x * x * y * y * y * z * z * z * z * z - 315.0 * x * x * y * z * z * z * z * z * z * z + 7.3828125 * y * y * y * y * y * y * y * y * y * z - 137.8125 * y * y * y * y * y * y * y * z * z * z + 433.125 * y * y * y * y * y * z * z * z * z * z - 315.0 * y * y * y * z * z * z * z * z * z * z) + e_2 * (-49.21875 * x * x * x * x * x * x * y * z - 344.53125 * x * x * x * x * y * y * y * z + 2165.625 * x * x * x * x * y * z * z * z - 383.90625 * x * x * y * y * y * y * y * z - 2756.25 * x * x * y * y * y * z * z * z - 1417.5 * x * x * y * z * z * z * z * z - 88.59375 * y * y * y * y * y * y * y * z + 748.125 * y * y * y * y * y * z * z * z - 1417.5 * y * y * y * z * z * z * z * z - 630.0 * y * z * z * z * z * z * z * z) + e_3 * (1575.0 * x * x * x * x * y * z - 4725.0 * x * x * y * y * y * z - 3150.0 * x * x * y * z * z * z - 3150.0 * y * y * y * z * z * z - 6300.0 * y * z * z * z * z * z) + e_4 * (-3543.75 * x * x * y * z - 3543.75 * y * y * y * z - 18900.0 * y * z * z * z) + e_5 * (-14175.0 * y * z);
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

        pc_52[k] = e_0 * (-std::sqrt(0.201873779296875) * x * x * x * x * x * x * x * x * x * x * x * y + std::sqrt(1.816864013671875) * x * x * x * x * x * x * x * x * x * y * y * y + std::sqrt(136.4666748046875) * x * x * x * x * x * x * x * x * x * y * z * z + std::sqrt(39.5672607421875) * x * x * x * x * x * x * x * y * y * y * y * y - std::sqrt(2183.466796875) * x * x * x * x * x * x * x * y * y * y * z * z - std::sqrt(6253.2421875) * x * x * x * x * x * x * x * y * z * z * z * z + std::sqrt(39.5672607421875) * x * x * x * x * x * y * y * y * y * y * y * y - std::sqrt(13646.66748046875) * x * x * x * x * x * y * y * y * y * y * z * z + std::sqrt(156331.0546875) * x * x * x * x * x * y * y * y * z * z * z * z + std::sqrt(5167.96875) * x * x * x * x * x * y * z * z * z * z * z * z + std::sqrt(1.816864013671875) * x * x * x * y * y * y * y * y * y * y * y * y - std::sqrt(2183.466796875) * x * x * x * y * y * y * y * y * y * y * z * z + std::sqrt(156331.0546875) * x * x * x * y * y * y * y * y * z * z * z * z - std::sqrt(186046.875) * x * x * x * y * y * y * z * z * z * z * z * z - std::sqrt(0.201873779296875) * x * y * y * y * y * y * y * y * y * y * y * y + std::sqrt(136.4666748046875) * x * y * y * y * y * y * y * y * y * y * z * z - std::sqrt(6253.2421875) * x * y * y * y * y * y * y * y * z * z * z * z + std::sqrt(5167.96875) * x * y * y * y * y * y * z * z * z * z * z * z) + e_1 * (-std::sqrt(80.74951171875) * x * x * x * x * x * x * x * x * x * y + std::sqrt(4186.0546875) * x * x * x * x * x * x * x * y * y * y - std::sqrt(37674.4921875) * x * x * x * x * x * x * x * y * z * z + std::sqrt(21718.388671875) * x * x * x * x * x * y * y * y * y * y + std::sqrt(49664.1796875) * x * x * x * x * x * y * y * y * z * z + std::sqrt(364651.875) * x * x * x * x * x * y * z * z * z * z + std::sqrt(4186.0546875) * x * x * x * y * y * y * y * y * y * y + std::sqrt(49664.1796875) * x * x * x * y * y * y * y * y * z * z + std::sqrt(2067187.5) * x * x * x * y * y * y * z * z * z * z - std::sqrt(330750.0) * x * x * x * y * z * z * z * z * z * z - std::sqrt(80.74951171875) * x * y * y * y * y * y * y * y * y * y - std::sqrt(37674.4921875) * x * y * y * y * y * y * y * y * z * z + std::sqrt(364651.875) * x * y * y * y * y * y * z * z * z * z - std::sqrt(330750.0) * x * y * y * y * z * z * z * z * z * z) + e_2 * (-std::sqrt(26162.841796875) * x * x * x * x * x * x * x * y + std::sqrt(2333660.888671875) * x * x * x * x * x * y * y * y + std::sqrt(11627.9296875) * x * x * x * x * x * y * z * z + std::sqrt(2333660.888671875) * x * x * x * y * y * y * y * y + std::sqrt(42795949.21875) * x * x * x * y * y * y * z * z + std::sqrt(744187.5) * x * x * x * y * z * z * z * z - std::sqrt(26162.841796875) * x * y * y * y * y * y * y * y + std::sqrt(11627.9296875) * x * y * y * y * y * y * z * z + std::sqrt(744187.5) * x * y * y * y * z * z * z * z - std::sqrt(2976750.0) * x * y * z * z * z * z * z * z) + e_3 * (std::sqrt(186046.875) * x * x * x * x * x * y + std::sqrt(152889187.5) * x * x * x * y * y * y + std::sqrt(74418750.0) * x * x * x * y * z * z + std::sqrt(186046.875) * x * y * y * y * y * y + std::sqrt(74418750.0) * x * y * y * y * z * z - std::sqrt(47628000.0) * x * y * z * z * z * z) + e_4 * (std::sqrt(156465421.875) * x * x * x * y + std::sqrt(156465421.875) * x * y * y * y + std::sqrt(6697687.5) * x * y * z * z) + e_5 * (std::sqrt(241116750.0) * x * y);
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

        pc_53[k] = e_0 * (std::sqrt(4.44122314453125) * x * x * x * x * x * x * x * x * x * x * y * z - std::sqrt(284.23828125) * x * x * x * x * x * x * x * x * y * y * y * z - std::sqrt(1136.953125) * x * x * x * x * x * x * x * x * y * z * z * z - std::sqrt(870.479736328125) * x * x * x * x * x * x * y * y * y * y * y * z + std::sqrt(92093.203125) * x * x * x * x * x * x * y * y * y * z * z * z + std::sqrt(1136.953125) * x * x * x * x * x * x * y * z * z * z * z * z + std::sqrt(28423.828125) * x * x * x * x * y * y * y * y * y * z * z * z - std::sqrt(113695.3125) * x * x * x * x * y * y * y * z * z * z * z * z + std::sqrt(111.03057861328125) * x * x * y * y * y * y * y * y * y * y * y * z - std::sqrt(28423.828125) * x * x * y * y * y * y * y * y * y * z * z * z + std::sqrt(28423.828125) * x * x * y * y * y * y * y * z * z * z * z * z) + e_1 * (-std::sqrt(3237.6516723632812) * x * x * x * x * x * x * x * x * y * z + std::sqrt(20536.2158203125) * x * x * x * x * x * x * y * y * y * z + std::sqrt(92093.203125) * x * x * x * x * x * x * y * z * z * z + std::sqrt(3997.100830078125) * x * x * x * x * y * y * y * y * y * z + std::sqrt(8214486.328125) * x * x * x * x * y * y * y * z * z * z - std::sqrt(255814.453125) * x * x * x * x * y * z * z * z * z * z - std::sqrt(15988.4033203125) * x * x * y * y * y * y * y * y * y * z - std::sqrt(710595.703125) * x * x * y * y * y * y * y * z * z * z - std::sqrt(113695.3125) * x * x * y * y * y * z * z * z * z * z + std::sqrt(111.03057861328125) * y * y * y * y * y * y * y * y * y * z - std::sqrt(28423.828125) * y * y * y * y * y * y * y * z * z * z + std::sqrt(28423.828125) * y * y * y * y * y * z * z * z * z * z) + e_2 * (-std::sqrt(15988.4033203125) * x * x * x * x * x * x * y * z + std::sqrt(32376516.723632812) * x * x * x * x * y * y * y * z + std::sqrt(16372125.0) * x * x * x * x * y * z * z * z - std::sqrt(5771813.5986328125) * x * x * y * y * y * y * y * z + std::sqrt(7276500.0) * x * x * y * y * y * z * z * z - std::sqrt(4093031.25) * x * x * y * z * z * z * z * z - std::sqrt(15988.4033203125) * y * y * y * y * y * y * y * z - std::sqrt(1819125.0) * y * y * y * y * y * z * z * z + std::sqrt(454781.25) * y * y * y * z * z * z * z * z) + e_3 * (std::sqrt(82883882.8125) * x * x * x * x * y * z + std::sqrt(36837281.25) * x * x * y * y * y * z + std::sqrt(16372125.0) * x * x * y * z * z * z - std::sqrt(9209320.3125) * y * y * y * y * y * z - std::sqrt(1819125.0) * y * y * y * z * z * z) + e_4 * (std::sqrt(451256695.3125) * x * x * y * z - std::sqrt(50139632.8125) * y * y * y * z);
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

        pc_54[k] = e_0 * (std::sqrt(0.3701019287109375) * x * x * x * x * x * x * x * x * x * x * x * y - std::sqrt(62.54722595214844) * x * x * x * x * x * x * x * x * x * y * y * y - std::sqrt(94.74609375) * x * x * x * x * x * x * x * x * x * y * z * z - std::sqrt(72.53997802734375) * x * x * x * x * x * x * x * y * y * y * y * y + std::sqrt(18570.234375) * x * x * x * x * x * x * x * y * y * y * z * z + std::sqrt(94.74609375) * x * x * x * x * x * x * x * y * z * z * z * z + std::sqrt(72.53997802734375) * x * x * x * x * x * y * y * y * y * y * y * y - std::sqrt(21317.87109375) * x * x * x * x * x * y * y * y * z * z * z * z + std::sqrt(62.54722595214844) * x * x * x * y * y * y * y * y * y * y * y * y - std::sqrt(18570.234375) * x * x * x * y * y * y * y * y * y * y * z * z + std::sqrt(21317.87109375) * x * x * x * y * y * y * y * y * z * z * z * z - std::sqrt(0.3701019287109375) * x * y * y * y * y * y * y * y * y * y * y * y + std::sqrt(94.74609375) * x * y * y * y * y * y * y * y * y * y * z * z - std::sqrt(94.74609375) * x * y * y * y * y * y * y * y * z * z * z * z) + e_1 * (-std::sqrt(54573.75) * x * x * x * x * x * x * x * y * y * y + std::sqrt(13643.4375) * x * x * x * x * x * x * x * y * z * z + std::sqrt(3942953.4375) * x * x * x * x * x * y * y * y * z * z - std::sqrt(54573.75) * x * x * x * x * x * y * z * z * z * z + std::sqrt(54573.75) * x * x * x * y * y * y * y * y * y * y - std::sqrt(3942953.4375) * x * x * x * y * y * y * y * y * z * z - std::sqrt(13643.4375) * x * y * y * y * y * y * y * y * z * z + std::sqrt(54573.75) * x * y * y * y * y * y * z * z * z * z) + e_2 * (-std::sqrt(85271.484375) * x * x * x * x * x * x * x * y - std::sqrt(2131787.109375) * x * x * x * x * x * y * y * y + std::sqrt(12279093.75) * x * x * x * x * x * y * z * z + std::sqrt(2131787.109375) * x * x * x * y * y * y * y * y - std::sqrt(1364343.75) * x * x * x * y * z * z * z * z + std::sqrt(85271.484375) * x * y * y * y * y * y * y * y - std::sqrt(12279093.75) * x * y * y * y * y * y * z * z + std::sqrt(1364343.75) * x * y * y * y * z * z * z * z) + e_3 * (-std::sqrt(5457375.0) * x * x * x * x * x * y + std::sqrt(87318000.0) * x * x * x * y * z * z + std::sqrt(5457375.0) * x * y * y * y * y * y - std::sqrt(87318000.0) * x * y * y * y * z * z) + e_4 * (-std::sqrt(12279093.75) * x * x * x * y + std::sqrt(12279093.75) * x * y * y * y);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ab_x, ab_y, ab_z : simd::cache_line_size())
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
        const auto e_6 = pe_6[k];

        pc_55[k] = e_0 * (8.203125 * x * x * x * x * x * x * x * x * y * y * z * z + 32.8125 * x * x * x * x * x * x * y * y * y * y * z * z - 65.625 * x * x * x * x * x * x * y * y * z * z * z * z + 49.21875 * x * x * x * x * y * y * y * y * y * y * z * z - 196.875 * x * x * x * x * y * y * y * y * z * z * z * z + 157.5 * x * x * x * x * y * y * z * z * z * z * z * z + 32.8125 * x * x * y * y * y * y * y * y * y * y * z * z - 196.875 * x * x * y * y * y * y * y * y * z * z * z * z + 315.0 * x * x * y * y * y * y * z * z * z * z * z * z - 105.0 * x * x * y * y * z * z * z * z * z * z * z * z + 8.203125 * y * y * y * y * y * y * y * y * y * y * z * z - 65.625 * y * y * y * y * y * y * y * y * z * z * z * z + 157.5 * y * y * y * y * y * y * z * z * z * z * z * z - 105.0 * y * y * y * y * z * z * z * z * z * z * z * z + 21.0 * y * y * z * z * z * z * z * z * z * z * z * z) + e_1 * (8.203125 * x * x * x * x * x * x * x * x * y * y + 8.203125 * x * x * x * x * x * x * x * x * z * z + 32.8125 * x * x * x * x * x * x * y * y * y * y + 32.8125 * x * x * x * x * x * x * y * y * z * z - 65.625 * x * x * x * x * x * x * z * z * z * z + 49.21875 * x * x * x * x * y * y * y * y * y * y + 49.21875 * x * x * x * x * y * y * y * y * z * z + 196.875 * x * x * x * x * y * y * z * z * z * z + 157.5 * x * x * x * x * z * z * z * z * z * z + 32.8125 * x * x * y * y * y * y * y * y * y * y + 32.8125 * x * x * y * y * y * y * y * y * z * z + 590.625 * x * x * y * y * y * y * z * z * z * z - 105.0 * x * x * y * y * z * z * z * z * z * z - 105.0 * x * x * z * z * z * z * z * z * z * z + 8.203125 * y * y * y * y * y * y * y * y * y * y + 8.203125 * y * y * y * y * y * y * y * y * z * z + 328.125 * y * y * y * y * y * y * z * z * z * z - 262.5 * y * y * y * y * z * z * z * z * z * z + 210.0 * y * y * z * z * z * z * z * z * z * z + 21.0 * z * z * z * z * z * z * z * z * z * z) + e_2 * (8.203125 * x * x * x * x * x * x * x * x + 229.6875 * x * x * x * x * x * x * y * y - 65.625 * x * x * x * x * x * x * z * z + 639.84375 * x * x * x * x * y * y * y * y + 984.375 * x * x * x * x * y * y * z * z + 787.5 * x * x * x * x * z * z * z * z + 623.4375 * x * x * y * y * y * y * y * y + 2165.625 * x * x * y * y * y * y * z * z + 1575.0 * x * x * y * y * z * z * z * z - 1050.0 * x * x * z * z * z * z * z * z + 205.078125 * y * y * y * y * y * y * y * y + 1115.625 * y * y * y * y * y * y * z * z + 787.5 * y * y * y * y * z * z * z * z + 2100.0 * y * y * z * z * z * z * z * z + 525.0 * z * z * z * z * z * z * z * z) + e_3 * (131.25 * x * x * x * x * x * x + 2756.25 * x * x * x * x * y * y + 1575.0 * x * x * x * x * z * z + 5118.75 * x * x * y * y * y * y + 9450.0 * x * x * y * y * z * z - 3150.0 * x * x * z * z * z * z + 2493.75 * y * y * y * y * y * y + 7875.0 * y * y * y * y * z * z + 12600.0 * y * y * z * z * z * z + 5250.0 * z * z * z * z * z * z) + e_4 * (1575.0 * x * x * x * x + 14175.0 * x * x * y * y + 12600.0 * y * y * y * y + 33075.0 * y * y * z * z + 22050.0 * z * z * z * z) + e_5 * (4725.0 * x * x + 24570.0 * y * y + 33075.0 * z * z) + e_6 * (10395.0);
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

        pc_56[k] = e_0 * (-std::sqrt(0.80108642578125) * x * x * x * x * x * x * x * x * x * x * y * z - std::sqrt(20.02716064453125) * x * x * x * x * x * x * x * x * y * y * y * z + std::sqrt(387.725830078125) * x * x * x * x * x * x * x * x * y * z * z * z - std::sqrt(80.108642578125) * x * x * x * x * x * x * y * y * y * y * y * z + std::sqrt(6203.61328125) * x * x * x * x * x * x * y * y * y * z * z * z - std::sqrt(7630.95703125) * x * x * x * x * x * x * y * z * z * z * z * z - std::sqrt(80.108642578125) * x * x * x * x * y * y * y * y * y * y * y * z + std::sqrt(13958.1298828125) * x * x * x * x * y * y * y * y * y * z * z * z - std::sqrt(68678.61328125) * x * x * x * x * y * y * y * z * z * z * z * z + std::sqrt(13125.0) * x * x * x * x * y * z * z * z * z * z * z * z - std::sqrt(20.02716064453125) * x * x * y * y * y * y * y * y * y * y * y * z + std::sqrt(6203.61328125) * x * x * y * y * y * y * y * y * y * z * z * z - std::sqrt(68678.61328125) * x * x * y * y * y * y * y * z * z * z * z * z + std::sqrt(52500.0) * x * x * y * y * y * z * z * z * z * z * z * z - std::sqrt(2100.0) * x * x * y * z * z * z * z * z * z * z * z * z - std::sqrt(0.80108642578125) * y * y * y * y * y * y * y * y * y * y * y * z + std::sqrt(387.725830078125) * y * y * y * y * y * y * y * y * y * z * z * z - std::sqrt(7630.95703125) * y * y * y * y * y * y * y * z * z * z * z * z + std::sqrt(13125.0) * y * y * y * y * y * z * z * z * z * z * z * z - std::sqrt(2100.0) * y * y * y * z * z * z * z * z * z * z * z * z + std::sqrt(21.0) * y * z * z * z * z * z * z * z * z * z * z * z) + e_1 * (std::sqrt(28.839111328125) * x * x * x * x * x * x * x * x * y * z + std::sqrt(461.42578125) * x * x * x * x * x * x * y * y * y * z - std::sqrt(7382.8125) * x * x * x * x * x * x * y * z * z * z + std::sqrt(1038.2080078125) * x * x * x * x * y * y * y * y * y * z - std::sqrt(66445.3125) * x * x * x * x * y * y * y * z * z * z + std::sqrt(95681.25) * x * x * x * x * y * z * z * z * z * z + std::sqrt(461.42578125) * x * x * y * y * y * y * y * y * y * z - std::sqrt(66445.3125) * x * x * y * y * y * y * y * z * z * z + std::sqrt(382725.0) * x * x * y * y * y * z * z * z * z * z - std::sqrt(75600.0) * x * x * y * z * z * z * z * z * z * z + std::sqrt(28.839111328125) * y * y * y * y * y * y * y * y * y * z - std::sqrt(7382.8125) * y * y * y * y * y * y * y * z * z * z + std::sqrt(95681.25) * y * y * y * y * y * z * z * z * z * z - std::sqrt(75600.0) * y * y * y * z * z * z * z * z * z * z + std::sqrt(4725.0) * y * z * z * z * z * z * z * z * z * z) + e_2 * (-std::sqrt(461.42578125) * x * x * x * x * x * x * y * z - std::sqrt(4152.83203125) * x * x * x * x * y * y * y * z + std::sqrt(265781.25) * x * x * x * x * y * z * z * z - std::sqrt(4152.83203125) * x * x * y * y * y * y * y * z + std::sqrt(1063125.0) * x * x * y * y * y * z * z * z - std::sqrt(1063125.0) * x * x * y * z * z * z * z * z - std::sqrt(461.42578125) * y * y * y * y * y * y * y * z + std::sqrt(265781.25) * y * y * y * y * y * z * z * z - std::sqrt(1063125.0) * y * y * y * z * z * z * z * z + std::sqrt(472500.0) * y * z * z * z * z * z * z * z) + e_3 * (std::sqrt(118125.0) * x * x * x * x * y * z + std::sqrt(472500.0) * x * x * y * y * y * z - std::sqrt(1890000.0) * x * x * y * z * z * z + std::sqrt(118125.0) * y * y * y * y * y * z - std::sqrt(1890000.0) * y * y * y * z * z * z + std::sqrt(11812500.0) * y * z * z * z * z * z) + e_4 * (std::sqrt(52093125.0) * y * z * z * z) + e_5 * (std::sqrt(18753525.0) * y * z);
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

        pc_57[k] = e_0 * (8.203125 * x * x * x * x * x * x * x * x * x * y * z * z + 32.8125 * x * x * x * x * x * x * x * y * y * y * z * z - 65.625 * x * x * x * x * x * x * x * y * z * z * z * z + 49.21875 * x * x * x * x * x * y * y * y * y * y * z * z - 196.875 * x * x * x * x * x * y * y * y * z * z * z * z + 157.5 * x * x * x * x * x * y * z * z * z * z * z * z + 32.8125 * x * x * x * y * y * y * y * y * y * y * z * z - 196.875 * x * x * x * y * y * y * y * y * z * z * z * z + 315.0 * x * x * x * y * y * y * z * z * z * z * z * z - 105.0 * x * x * x * y * z * z * z * z * z * z * z * z + 8.203125 * x * y * y * y * y * y * y * y * y * y * z * z - 65.625 * x * y * y * y * y * y * y * y * z * z * z * z + 157.5 * x * y * y * y * y * y * z * z * z * z * z * z - 105.0 * x * y * y * y * z * z * z * z * z * z * z * z + 21.0 * x * y * z * z * z * z * z * z * z * z * z * z) + e_1 * (8.203125 * x * x * x * x * x * x * x * x * x * y + 32.8125 * x * x * x * x * x * x * x * y * y * y + 49.21875 * x * x * x * x * x * y * y * y * y * y + 393.75 * x * x * x * x * x * y * z * z * z * z + 32.8125 * x * x * x * y * y * y * y * y * y * y + 787.5 * x * x * x * y * y * y * z * z * z * z - 420.0 * x * x * x * y * z * z * z * z * z * z + 8.203125 * x * y * y * y * y * y * y * y * y * y + 393.75 * x * y * y * y * y * y * z * z * z * z - 420.0 * x * y * y * y * z * z * z * z * z * z + 315.0 * x * y * z * z * z * z * z * z * z * z) + e_2 * (196.875 * x * x * x * x * x * x * x * y + 590.625 * x * x * x * x * x * y * y * y + 1181.25 * x * x * x * x * x * y * z * z + 590.625 * x * x * x * y * y * y * y * y + 2362.5 * x * x * x * y * y * y * z * z + 196.875 * x * y * y * y * y * y * y * y + 1181.25 * x * y * y * y * y * y * z * z + 3150.0 * x * y * z * z * z * z * z * z) + e_3 * (2362.5 * x * x * x * x * x * y + 4725.0 * x * x * x * y * y * y + 6300.0 * x * x * x * y * z * z + 2362.5 * x * y * y * y * y * y + 6300.0 * x * y * y * y * z * z + 15750.0 * x * y * z * z * z * z) + e_4 * (11025.0 * x * x * x * y + 11025.0 * x * y * y * y + 33075.0 * x * y * z * z) + e_5 * (19845.0 * x * y);
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

        pc_58[k] = e_0 * (std::sqrt(1.682281494140625) * x * x * x * x * x * x * x * x * x * x * y * z + std::sqrt(15.140533447265625) * x * x * x * x * x * x * x * x * y * y * y * z - std::sqrt(672.91259765625) * x * x * x * x * x * x * x * x * y * z * z * z + std::sqrt(6.7291259765625) * x * x * x * x * x * x * y * y * y * y * y * z - std::sqrt(2691.650390625) * x * x * x * x * x * x * y * y * y * z * z * z + std::sqrt(11201.572265625) * x * x * x * x * x * x * y * z * z * z * z * z - std::sqrt(6.7291259765625) * x * x * x * x * y * y * y * y * y * y * y * z + std::sqrt(11201.572265625) * x * x * x * x * y * y * y * z * z * z * z * z - std::sqrt(13505.625) * x * x * x * x * y * z * z * z * z * z * z * z - std::sqrt(15.140533447265625) * x * x * y * y * y * y * y * y * y * y * y * z + std::sqrt(2691.650390625) * x * x * y * y * y * y * y * y * y * z * z * z - std::sqrt(11201.572265625) * x * x * y * y * y * y * y * z * z * z * z * z + std::sqrt(1102.5) * x * x * y * z * z * z * z * z * z * z * z * z - std::sqrt(1.682281494140625) * y * y * y * y * y * y * y * y * y * y * y * z + std::sqrt(672.91259765625) * y * y * y * y * y * y * y * y * y * z * z * z - std::sqrt(11201.572265625) * y * y * y * y * y * y * y * z * z * z * z * z + std::sqrt(13505.625) * y * y * y * y * y * z * z * z * z * z * z * z - std::sqrt(1102.5) * y * y * y * z * z * z * z * z * z * z * z * z) + e_1 * (-std::sqrt(60.5621337890625) * x * x * x * x * x * x * x * x * y * z - std::sqrt(430.6640625) * x * x * x * x * x * x * y * y * y * z + std::sqrt(31115.478515625) * x * x * x * x * x * x * y * z * z * z - std::sqrt(242.24853515625) * x * x * x * x * y * y * y * y * y * z + std::sqrt(78488.525390625) * x * x * x * x * y * y * y * z * z * z - std::sqrt(286667.2265625) * x * x * x * x * y * z * z * z * z * z + std::sqrt(968.994140625) * x * x * y * y * y * y * y * z * z * z - std::sqrt(179225.15625) * x * x * y * y * y * z * z * z * z * z + std::sqrt(248062.5) * x * x * y * z * z * z * z * z * z * z + std::sqrt(6.7291259765625) * y * y * y * y * y * y * y * y * y * z - std::sqrt(5275.634765625) * y * y * y * y * y * y * y * z * z * z + std::sqrt(12558.1640625) * y * y * y * y * y * z * z * z * z * z - std::sqrt(1102.5) * y * y * y * z * z * z * z * z * z * z - std::sqrt(4410.0) * y * z * z * z * z * z * z * z * z * z) + e_2 * (std::sqrt(15503.90625) * x * x * x * x * x * x * y * z + std::sqrt(34883.7890625) * x * x * x * x * y * y * y * z - std::sqrt(872094.7265625) * x * x * x * x * y * z * z * z - std::sqrt(1255816.40625) * x * x * y * y * y * z * z * z + std::sqrt(8930250.0) * x * x * y * z * z * z * z * z - std::sqrt(3875.9765625) * y * y * y * y * y * y * y * z - std::sqrt(34883.7890625) * y * y * y * y * y * z * z * z - std::sqrt(992250.0) * y * z * z * z * z * z * z * z) + e_3 * (-std::sqrt(15503.90625) * x * x * x * x * y * z - std::sqrt(558140.625) * x * x * y * y * y * z + std::sqrt(48620250.0) * x * x * y * z * z * z - std::sqrt(387597.65625) * y * y * y * y * y * z - std::sqrt(992250.0) * y * y * y * z * z * z - std::sqrt(35721000.0) * y * z * z * z * z * z) + e_4 * (std::sqrt(20093062.5) * x * x * y * z - std::sqrt(6201562.5) * y * y * y * z - std::sqrt(194481000.0) * y * z * z * z) + e_5 * (-std::sqrt(80372250.0) * y * z);
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

        pc_59[k] = e_0 * (-std::sqrt(60.5621337890625) * x * x * x * x * x * x * x * x * x * y * z * z + std::sqrt(2691.650390625) * x * x * x * x * x * x * x * y * z * z * z * z + std::sqrt(2180.23681640625) * x * x * x * x * x * y * y * y * y * y * z * z - std::sqrt(2691.650390625) * x * x * x * x * x * y * y * y * z * z * z * z - std::sqrt(9112.8515625) * x * x * x * x * x * y * z * z * z * z * z * z + std::sqrt(3875.9765625) * x * x * x * y * y * y * y * y * y * y * z * z - std::sqrt(67291.259765625) * x * x * x * y * y * y * y * y * z * z * z * z + std::sqrt(36451.40625) * x * x * x * y * y * y * z * z * z * z * z * z + std::sqrt(1102.5) * x * x * x * y * z * z * z * z * z * z * z * z + std::sqrt(545.0592041015625) * x * y * y * y * y * y * y * y * y * y * z * z - std::sqrt(24224.853515625) * x * y * y * y * y * y * y * y * z * z * z * z + std::sqrt(82015.6640625) * x * y * y * y * y * y * z * z * z * z * z * z - std::sqrt(9922.5) * x * y * y * y * z * z * z * z * z * z * z * z) + e_1 * (-std::sqrt(60.5621337890625) * x * x * x * x * x * x * x * x * x * y + std::sqrt(968.994140625) * x * x * x * x * x * x * x * y * z * z + std::sqrt(2180.23681640625) * x * x * x * x * x * y * y * y * y * y + std::sqrt(24224.853515625) * x * x * x * x * x * y * y * y * z * z - std::sqrt(248062.5) * x * x * x * x * x * y * z * z * z * z + std::sqrt(3875.9765625) * x * x * x * y * y * y * y * y * y * y + std::sqrt(47480.712890625) * x * x * x * y * y * y * y * y * z * z - std::sqrt(62015.625) * x * x * x * y * y * y * z * z * z * z + std::sqrt(300155.625) * x * x * x * y * z * z * z * z * z * z + std::sqrt(545.0592041015625) * x * y * y * y * y * y * y * y * y * y + std::sqrt(8720.947265625) * x * y * y * y * y * y * y * y * z * z + std::sqrt(62015.625) * x * y * y * y * y * y * z * z * z * z + std::sqrt(419225.625) * x * y * y * y * z * z * z * z * z * z - std::sqrt(39690.0) * x * y * z * z * z * z * z * z * z * z) + e_2 * (-std::sqrt(15503.90625) * x * x * x * x * x * x * x * y + std::sqrt(96899.4140625) * x * x * x * x * x * y * y * y - std::sqrt(872094.7265625) * x * x * x * x * x * y * z * z + std::sqrt(992250.0) * x * x * x * y * y * y * y * y + std::sqrt(1255816.40625) * x * x * x * y * y * y * z * z + std::sqrt(1550390.625) * x * x * x * y * z * z * z * z + std::sqrt(313954.1015625) * x * y * y * y * y * y * y * y + std::sqrt(4220938.4765625) * x * y * y * y * y * y * z * z + std::sqrt(32806265.625) * x * y * y * y * z * z * z * z - std::sqrt(992250.0) * x * y * z * z * z * z * z * z) + e_3 * (-std::sqrt(759691.40625) * x * x * x * x * x * y + std::sqrt(22387640.625) * x * x * x * y * y * y + std::sqrt(248062.5) * x * x * x * y * z * z + std::sqrt(31395410.15625) * x * y * y * y * y * y + std::sqrt(377303062.5) * x * y * y * y * z * z + std::sqrt(3969000.0) * x * y * z * z * z * z) + e_4 * (std::sqrt(2232562.5) * x * x * x * y + std::sqrt(502326562.5) * x * y * y * y + std::sqrt(321489000.0) * x * y * z * z) + e_5 * (std::sqrt(321489000.0) * x * y);
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

        pc_60[k] = e_0 * (-std::sqrt(2.01873779296875) * x * x * x * x * x * x * x * x * x * x * y * z + std::sqrt(18.16864013671875) * x * x * x * x * x * x * x * x * y * y * y * z + std::sqrt(395.672607421875) * x * x * x * x * x * x * x * x * y * z * z * z + std::sqrt(395.672607421875) * x * x * x * x * x * x * y * y * y * y * y * z - std::sqrt(6330.76171875) * x * x * x * x * x * x * y * y * y * z * z * z - std::sqrt(3493.546875) * x * x * x * x * x * x * y * z * z * z * z * z + std::sqrt(395.672607421875) * x * x * x * x * y * y * y * y * y * y * y * z - std::sqrt(39567.2607421875) * x * x * x * x * y * y * y * y * y * z * z * z + std::sqrt(87338.671875) * x * x * x * x * y * y * y * z * z * z * z * z + std::sqrt(516.796875) * x * x * x * x * y * z * z * z * z * z * z * z + std::sqrt(18.16864013671875) * x * x * y * y * y * y * y * y * y * y * y * z - std::sqrt(6330.76171875) * x * x * y * y * y * y * y * y * y * z * z * z + std::sqrt(87338.671875) * x * x * y * y * y * y * y * z * z * z * z * z - std::sqrt(18604.6875) * x * x * y * y * y * z * z * z * z * z * z * z - std::sqrt(2.01873779296875) * y * y * y * y * y * y * y * y * y * y * y * z + std::sqrt(395.672607421875) * y * y * y * y * y * y * y * y * y * z * z * z - std::sqrt(3493.546875) * y * y * y * y * y * y * y * z * z * z * z * z + std::sqrt(516.796875) * y * y * y * y * y * z * z * z * z * z * z * z) + e_1 * (std::sqrt(72.674560546875) * x * x * x * x * x * x * x * x * y * z + std::sqrt(6330.76171875) * x * x * x * x * x * x * y * y * y * z - std::sqrt(74418.75) * x * x * x * x * x * x * y * z * z * z + std::sqrt(14244.2138671875) * x * x * x * x * y * y * y * y * y * z - std::sqrt(51679.6875) * x * x * x * x * y * y * y * z * z * z + std::sqrt(227907.421875) * x * x * x * x * y * z * z * z * z * z + std::sqrt(1162.79296875) * x * x * y * y * y * y * y * y * y * z + std::sqrt(8268.75) * x * x * y * y * y * y * y * z * z * z + std::sqrt(3474942.1875) * x * x * y * y * y * z * z * z * z * z - std::sqrt(74418.75) * x * x * y * z * z * z * z * z * z * z - std::sqrt(201.873779296875) * y * y * y * y * y * y * y * y * y * z + std::sqrt(2067.1875) * y * y * y * y * y * y * y * z * z * z - std::sqrt(219307.921875) * y * y * y * y * y * z * z * z * z * z + std::sqrt(8268.75) * y * y * y * z * z * z * z * z * z * z) + e_2 * (-std::sqrt(29069.82421875) * x * x * x * x * x * x * y * z + std::sqrt(726745.60546875) * x * x * x * x * y * y * y * z + std::sqrt(726745.60546875) * x * x * y * y * y * y * y * z + std::sqrt(82687500.0) * x * x * y * y * y * z * z * z + std::sqrt(1860468.75) * x * x * y * z * z * z * z * z - std::sqrt(29069.82421875) * y * y * y * y * y * y * y * z - std::sqrt(3307500.0) * y * y * y * y * y * z * z * z - std::sqrt(206718.75) * y * y * y * z * z * z * z * z) + e_3 * (std::sqrt(186046875.0) * x * x * y * y * y * z + std::sqrt(186046875.0) * x * x * y * z * z * z - std::sqrt(7441875.0) * y * y * y * y * y * z - std::sqrt(20671875.0) * y * y * y * z * z * z) + e_4 * (std::sqrt(418605468.75) * x * x * y * z - std::sqrt(46511718.75) * y * y * y * z);
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

        pc_61[k] = e_0 * (std::sqrt(44.4122314453125) * x * x * x * x * x * x * x * x * x * y * z * z - std::sqrt(2842.3828125) * x * x * x * x * x * x * x * y * y * y * z * z - std::sqrt(710.595703125) * x * x * x * x * x * x * x * y * z * z * z * z - std::sqrt(8704.79736328125) * x * x * x * x * x * y * y * y * y * y * z * z + std::sqrt(57558.251953125) * x * x * x * x * x * y * y * y * z * z * z * z + std::sqrt(113.6953125) * x * x * x * x * x * y * z * z * z * z * z * z + std::sqrt(17764.892578125) * x * x * x * y * y * y * y * y * z * z * z * z - std::sqrt(11369.53125) * x * x * x * y * y * y * z * z * z * z * z * z + std::sqrt(1110.3057861328125) * x * y * y * y * y * y * y * y * y * y * z * z - std::sqrt(17764.892578125) * x * y * y * y * y * y * y * y * z * z * z * z + std::sqrt(2842.3828125) * x * y * y * y * y * y * z * z * z * z * z * z) + e_1 * (std::sqrt(44.4122314453125) * x * x * x * x * x * x * x * x * x * y - std::sqrt(2842.3828125) * x * x * x * x * x * x * x * y * y * y - std::sqrt(6395.361328125) * x * x * x * x * x * x * x * y * z * z - std::sqrt(8704.79736328125) * x * x * x * x * x * y * y * y * y * y - std::sqrt(375905.126953125) * x * x * x * x * x * y * y * y * z * z + std::sqrt(102325.78125) * x * x * x * x * x * y * z * z * z * z - std::sqrt(17764.892578125) * x * x * x * y * y * y * y * y * z * z + std::sqrt(4547812.5) * x * x * x * y * y * y * z * z * z * z - std::sqrt(45478.125) * x * x * x * y * z * z * z * z * z * z + std::sqrt(1110.3057861328125) * x * y * y * y * y * y * y * y * y * y + std::sqrt(159884.033203125) * x * y * y * y * y * y * y * y * z * z - std::sqrt(2558144.53125) * x * y * y * y * y * y * z * z * z * z + std::sqrt(45478.125) * x * y * y * y * z * z * z * z * z * z) + e_2 * (-std::sqrt(1776489.2578125) * x * x * x * x * x * y * y * y - std::sqrt(639536.1328125) * x * x * x * x * x * y * z * z - std::sqrt(284238.28125) * x * x * x * y * y * y * y * y + std::sqrt(7105957.03125) * x * x * x * y * y * y * z * z + std::sqrt(10232578.125) * x * x * x * y * z * z * z * z + std::sqrt(639536.1328125) * x * y * y * y * y * y * y * y - std::sqrt(639536.1328125) * x * y * y * y * y * y * z * z - std::sqrt(10232578.125) * x * y * y * y * z * z * z * z) + e_3 * (-std::sqrt(2558144.53125) * x * x * x * x * x * y - std::sqrt(28423828.125) * x * x * x * y * y * y + std::sqrt(40930312.5) * x * x * x * y * z * z + std::sqrt(23023300.78125) * x * y * y * y * y * y - std::sqrt(40930312.5) * x * y * y * y * z * z) + e_4 * (-std::sqrt(40930312.5) * x * x * x * y + std::sqrt(40930312.5) * x * y * y * y);
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

        pc_62[k] = e_0 * (std::sqrt(3.701019287109375) * x * x * x * x * x * x * x * x * x * x * y * z - std::sqrt(625.4722595214844) * x * x * x * x * x * x * x * x * y * y * y * z - std::sqrt(59.21630859375) * x * x * x * x * x * x * x * x * y * z * z * z - std::sqrt(725.3997802734375) * x * x * x * x * x * x * y * y * y * y * y * z + std::sqrt(11606.396484375) * x * x * x * x * x * x * y * y * y * z * z * z + std::sqrt(9.474609375) * x * x * x * x * x * x * y * z * z * z * z * z + std::sqrt(725.3997802734375) * x * x * x * x * y * y * y * y * y * y * y * z - std::sqrt(2131.787109375) * x * x * x * x * y * y * y * z * z * z * z * z + std::sqrt(625.4722595214844) * x * x * y * y * y * y * y * y * y * y * y * z - std::sqrt(11606.396484375) * x * x * y * y * y * y * y * y * y * z * z * z + std::sqrt(2131.787109375) * x * x * y * y * y * y * y * z * z * z * z * z - std::sqrt(3.701019287109375) * y * y * y * y * y * y * y * y * y * y * y * z + std::sqrt(59.21630859375) * y * y * y * y * y * y * y * y * y * z * z * z - std::sqrt(9.474609375) * y * y * y * y * y * y * y * z * z * z * z * z) + e_1 * (-std::sqrt(133.2366943359375) * x * x * x * x * x * x * x * x * y * z - std::sqrt(417830.2734375) * x * x * x * x * x * x * y * y * y * z + std::sqrt(19186.083984375) * x * x * x * x * x * x * y * z * z * z + std::sqrt(26114.39208984375) * x * x * x * x * y * y * y * y * y * z + std::sqrt(1332366.943359375) * x * x * x * x * y * y * y * z * z * z - std::sqrt(8527.1484375) * x * x * x * x * y * z * z * z * z * z + std::sqrt(545737.5) * x * x * y * y * y * y * y * y * y * z - std::sqrt(3242448.193359375) * x * x * y * y * y * y * y * z * z * z + std::sqrt(34108.59375) * x * x * y * y * y * z * z * z * z * z - std::sqrt(3330.9173583984375) * y * y * y * y * y * y * y * y * y * z + std::sqrt(19186.083984375) * y * y * y * y * y * y * y * z * z * z - std::sqrt(341.0859375) * y * y * y * y * y * z * z * z * z * z) + e_2 * (-std::sqrt(852714.84375) * x * x * x * x * x * x * y * z - std::sqrt(5329467.7734375) * x * x * x * x * y * y * y * z + std::sqrt(5329467.7734375) * x * x * x * x * y * z * z * z + std::sqrt(30697734.375) * x * x * y * y * y * y * y * z - std::sqrt(21317871.09375) * x * x * y * y * y * z * z * z - std::sqrt(213178.7109375) * y * y * y * y * y * y * y * z + std::sqrt(213178.7109375) * y * y * y * y * y * z * z * z) + e_3 * (-std::sqrt(21317871.09375) * x * x * x * x * y * z + std::sqrt(85271484.375) * x * x * y * y * y * z - std::sqrt(852714.84375) * y * y * y * y * y * z);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ab_x, ab_y, ab_z : simd::cache_line_size())
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
        const auto e_6 = pe_6[k];

        pc_63[k] = e_0 * (0.09765625 * x * x * x * x * x * x * x * x * x * x * x * x + 0.5859375 * x * x * x * x * x * x * x * x * x * x * y * y - 3.515625 * x * x * x * x * x * x * x * x * x * x * z * z + 1.46484375 * x * x * x * x * x * x * x * x * y * y * y * y - 17.578125 * x * x * x * x * x * x * x * x * y * y * z * z + 36.328125 * x * x * x * x * x * x * x * x * z * z * z * z + 1.953125 * x * x * x * x * x * x * y * y * y * y * y * y - 35.15625 * x * x * x * x * x * x * y * y * y * y * z * z + 145.3125 * x * x * x * x * x * x * y * y * z * z * z * z - 85.0 * x * x * x * x * x * x * z * z * z * z * z * z + 1.46484375 * x * x * x * x * y * y * y * y * y * y * y * y - 35.15625 * x * x * x * x * y * y * y * y * y * y * z * z + 217.96875 * x * x * x * x * y * y * y * y * z * z * z * z - 255.0 * x * x * x * x * y * y * z * z * z * z * z * z + 67.5 * x * x * x * x * z * z * z * z * z * z * z * z + 0.5859375 * x * x * y * y * y * y * y * y * y * y * y * y - 17.578125 * x * x * y * y * y * y * y * y * y * y * z * z + 145.3125 * x * x * y * y * y * y * y * y * z * z * z * z - 255.0 * x * x * y * y * y * y * z * z * z * z * z * z + 135.0 * x * x * y * y * z * z * z * z * z * z * z * z - 15.0 * x * x * z * z * z * z * z * z * z * z * z * z + 0.09765625 * y * y * y * y * y * y * y * y * y * y * y * y - 3.515625 * y * y * y * y * y * y * y * y * y * y * z * z + 36.328125 * y * y * y * y * y * y * y * y * z * z * z * z - 85.0 * y * y * y * y * y * y * z * z * z * z * z * z + 67.5 * y * y * y * y * z * z * z * z * z * z * z * z - 15.0 * y * y * z * z * z * z * z * z * z * z * z * z + z * z * z * z * z * z * z * z * z * z * z * z) + e_1 * (3.515625 * x * x * x * x * x * x * x * x * x * x + 17.578125 * x * x * x * x * x * x * x * x * y * y + 42.1875 * x * x * x * x * x * x * x * x * z * z + 35.15625 * x * x * x * x * x * x * y * y * y * y + 168.75 * x * x * x * x * x * x * y * y * z * z - 112.5 * x * x * x * x * x * x * z * z * z * z + 35.15625 * x * x * x * x * y * y * y * y * y * y + 253.125 * x * x * x * x * y * y * y * y * z * z - 337.5 * x * x * x * x * y * y * z * z * z * z + 360.0 * x * x * x * x * z * z * z * z * z * z + 17.578125 * x * x * y * y * y * y * y * y * y * y + 168.75 * x * x * y * y * y * y * y * y * z * z - 337.5 * x * x * y * y * y * y * z * z * z * z + 720.0 * x * x * y * y * z * z * z * z * z * z - 135.0 * x * x * z * z * z * z * z * z * z * z + 3.515625 * y * y * y * y * y * y * y * y * y * y + 42.1875 * y * y * y * y * y * y * y * y * z * z - 112.5 * y * y * y * y * y * y * z * z * z * z + 360.0 * y * y * y * y * z * z * z * z * z * z - 135.0 * y * y * z * z * z * z * z * z * z * z + 36.0 * z * z * z * z * z * z * z * z * z * z) + e_2 * (108.984375 * x * x * x * x * x * x * x * x + 435.9375 * x * x * x * x * x * x * y * y + 337.5 * x * x * x * x * x * x * z * z + 653.90625 * x * x * x * x * y * y * y * y + 1012.5 * x * x * x * x * y * y * z * z + 1687.5 * x * x * x * x * z * z * z * z + 435.9375 * x * x * y * y * y * y * y * y + 1012.5 * x * x * y * y * y * y * z * z + 3375.0 * x * x * y * y * z * z * z * z - 450.0 * x * x * z * z * z * z * z * z + 108.984375 * y * y * y * y * y * y * y * y + 337.5 * y * y * y * y * y * y * z * z + 1687.5 * y * y * y * y * z * z * z * z - 450.0 * y * y * z * z * z * z * z * z + 675.0 * z * z * z * z * z * z * z * z) + e_3 * (1275.0 * x * x * x * x * x * x + 3825.0 * x * x * x * x * y * y + 5400.0 * x * x * x * x * z * z + 3825.0 * x * x * y * y * y * y + 10800.0 * x * x * y * y * z * z + 2250.0 * x * x * z * z * z * z + 1275.0 * y * y * y * y * y * y + 5400.0 * y * y * y * y * z * z + 2250.0 * y * y * z * z * z * z + 6000.0 * z * z * z * z * z * z) + e_4 * (7087.5 * x * x * x * x + 14175.0 * x * x * y * y + 14175.0 * x * x * z * z + 7087.5 * y * y * y * y + 14175.0 * y * y * z * z + 23625.0 * z * z * z * z) + e_5 * (14175.0 * x * x + 14175.0 * y * y + 34020.0 * z * z) + e_6 * (10395.0);
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

        pc_64[k] = e_0 * (-std::sqrt(0.80108642578125) * x * x * x * x * x * x * x * x * x * x * x * z - std::sqrt(20.02716064453125) * x * x * x * x * x * x * x * x * x * y * y * z + std::sqrt(387.725830078125) * x * x * x * x * x * x * x * x * x * z * z * z - std::sqrt(80.108642578125) * x * x * x * x * x * x * x * y * y * y * y * z + std::sqrt(6203.61328125) * x * x * x * x * x * x * x * y * y * z * z * z - std::sqrt(7630.95703125) * x * x * x * x * x * x * x * z * z * z * z * z - std::sqrt(80.108642578125) * x * x * x * x * x * y * y * y * y * y * y * z + std::sqrt(13958.1298828125) * x * x * x * x * x * y * y * y * y * z * z * z - std::sqrt(68678.61328125) * x * x * x * x * x * y * y * z * z * z * z * z + std::sqrt(13125.0) * x * x * x * x * x * z * z * z * z * z * z * z - std::sqrt(20.02716064453125) * x * x * x * y * y * y * y * y * y * y * y * z + std::sqrt(6203.61328125) * x * x * x * y * y * y * y * y * y * z * z * z - std::sqrt(68678.61328125) * x * x * x * y * y * y * y * z * z * z * z * z + std::sqrt(52500.0) * x * x * x * y * y * z * z * z * z * z * z * z - std::sqrt(2100.0) * x * x * x * z * z * z * z * z * z * z * z * z - std::sqrt(0.80108642578125) * x * y * y * y * y * y * y * y * y * y * y * z + std::sqrt(387.725830078125) * x * y * y * y * y * y * y * y * y * z * z * z - std::sqrt(7630.95703125) * x * y * y * y * y * y * y * z * z * z * z * z + std::sqrt(13125.0) * x * y * y * y * y * z * z * z * z * z * z * z - std::sqrt(2100.0) * x * y * y * z * z * z * z * z * z * z * z * z + std::sqrt(21.0) * x * z * z * z * z * z * z * z * z * z * z * z) + e_1 * (std::sqrt(28.839111328125) * x * x * x * x * x * x * x * x * x * z + std::sqrt(461.42578125) * x * x * x * x * x * x * x * y * y * z - std::sqrt(7382.8125) * x * x * x * x * x * x * x * z * z * z + std::sqrt(1038.2080078125) * x * x * x * x * x * y * y * y * y * z - std::sqrt(66445.3125) * x * x * x * x * x * y * y * z * z * z + std::sqrt(95681.25) * x * x * x * x * x * z * z * z * z * z + std::sqrt(461.42578125) * x * x * x * y * y * y * y * y * y * z - std::sqrt(66445.3125) * x * x * x * y * y * y * y * z * z * z + std::sqrt(382725.0) * x * x * x * y * y * z * z * z * z * z - std::sqrt(75600.0) * x * x * x * z * z * z * z * z * z * z + std::sqrt(28.839111328125) * x * y * y * y * y * y * y * y * y * z - std::sqrt(7382.8125) * x * y * y * y * y * y * y * z * z * z + std::sqrt(95681.25) * x * y * y * y * y * z * z * z * z * z - std::sqrt(75600.0) * x * y * y * z * z * z * z * z * z * z + std::sqrt(4725.0) * x * z * z * z * z * z * z * z * z * z) + e_2 * (-std::sqrt(461.42578125) * x * x * x * x * x * x * x * z - std::sqrt(4152.83203125) * x * x * x * x * x * y * y * z + std::sqrt(265781.25) * x * x * x * x * x * z * z * z - std::sqrt(4152.83203125) * x * x * x * y * y * y * y * z + std::sqrt(1063125.0) * x * x * x * y * y * z * z * z - std::sqrt(1063125.0) * x * x * x * z * z * z * z * z - std::sqrt(461.42578125) * x * y * y * y * y * y * y * z + std::sqrt(265781.25) * x * y * y * y * y * z * z * z - std::sqrt(1063125.0) * x * y * y * z * z * z * z * z + std::sqrt(472500.0) * x * z * z * z * z * z * z * z) + e_3 * (std::sqrt(118125.0) * x * x * x * x * x * z + std::sqrt(472500.0) * x * x * x * y * y * z - std::sqrt(1890000.0) * x * x * x * z * z * z + std::sqrt(118125.0) * x * y * y * y * y * z - std::sqrt(1890000.0) * x * y * y * z * z * z + std::sqrt(11812500.0) * x * z * z * z * z * z) + e_4 * (std::sqrt(52093125.0) * x * z * z * z) + e_5 * (std::sqrt(18753525.0) * x * z);
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

        pc_65[k] = e_0 * (-std::sqrt(0.02002716064453125) * x * x * x * x * x * x * x * x * x * x * x * x - std::sqrt(0.3204345703125) * x * x * x * x * x * x * x * x * x * x * y * y + std::sqrt(23.151397705078125) * x * x * x * x * x * x * x * x * x * x * z * z - std::sqrt(0.5006790161132812) * x * x * x * x * x * x * x * x * y * y * y * y + std::sqrt(208.36257934570312) * x * x * x * x * x * x * x * x * y * y * z * z - std::sqrt(2154.60205078125) * x * x * x * x * x * x * x * x * z * z * z * z + std::sqrt(92.6055908203125) * x * x * x * x * x * x * y * y * y * y * z * z - std::sqrt(8618.408203125) * x * x * x * x * x * x * y * y * z * z * z * z + std::sqrt(9130.283203125) * x * x * x * x * x * x * z * z * z * z * z * z + std::sqrt(0.5006790161132812) * x * x * x * x * y * y * y * y * y * y * y * y - std::sqrt(92.6055908203125) * x * x * x * x * y * y * y * y * y * y * z * z + std::sqrt(9130.283203125) * x * x * x * x * y * y * z * z * z * z * z * z - std::sqrt(3793.125) * x * x * x * x * z * z * z * z * z * z * z * z + std::sqrt(0.3204345703125) * x * x * y * y * y * y * y * y * y * y * y * y - std::sqrt(208.36257934570312) * x * x * y * y * y * y * y * y * y * y * z * z + std::sqrt(8618.408203125) * x * x * y * y * y * y * y * y * z * z * z * z - std::sqrt(9130.283203125) * x * x * y * y * y * y * z * z * z * z * z * z + std::sqrt(52.5) * x * x * z * z * z * z * z * z * z * z * z * z + std::sqrt(0.02002716064453125) * y * y * y * y * y * y * y * y * y * y * y * y - std::sqrt(23.151397705078125) * y * y * y * y * y * y * y * y * y * y * z * z + std::sqrt(2154.60205078125) * y * y * y * y * y * y * y * y * z * z * z * z - std::sqrt(9130.283203125) * y * y * y * y * y * y * z * z * z * z * z * z + std::sqrt(3793.125) * y * y * y * y * z * z * z * z * z * z * z * z - std::sqrt(52.5) * y * y * z * z * z * z * z * z * z * z * z * z) + e_1 * (-std::sqrt(25.9552001953125) * x * x * x * x * x * x * x * x * x * x - std::sqrt(233.5968017578125) * x * x * x * x * x * x * x * x * y * y - std::sqrt(2260.986328125) * x * x * x * x * x * x * x * x * z * z - std::sqrt(103.82080078125) * x * x * x * x * x * x * y * y * y * y - std::sqrt(9043.9453125) * x * x * x * x * x * x * y * y * z * z + std::sqrt(1661.1328125) * x * x * x * x * x * x * z * z * z * z + std::sqrt(103.82080078125) * x * x * x * x * y * y * y * y * y * y + std::sqrt(1661.1328125) * x * x * x * x * y * y * z * z * z * z - std::sqrt(38272.5) * x * x * x * x * z * z * z * z * z * z + std::sqrt(233.5968017578125) * x * x * y * y * y * y * y * y * y * y + std::sqrt(9043.9453125) * x * x * y * y * y * y * y * y * z * z - std::sqrt(1661.1328125) * x * x * y * y * y * y * z * z * z * z - std::sqrt(1890.0) * x * x * z * z * z * z * z * z * z * z + std::sqrt(25.9552001953125) * y * y * y * y * y * y * y * y * y * y + std::sqrt(2260.986328125) * y * y * y * y * y * y * y * y * z * z - std::sqrt(1661.1328125) * y * y * y * y * y * y * z * z * z * z + std::sqrt(38272.5) * y * y * y * y * z * z * z * z * z * z + std::sqrt(1890.0) * y * y * z * z * z * z * z * z * z * z) + e_2 * (-std::sqrt(21329.40673828125) * x * x * x * x * x * x * x * x - std::sqrt(85317.626953125) * x * x * x * x * x * x * y * y - std::sqrt(349253.173828125) * x * x * x * x * x * x * z * z - std::sqrt(349253.173828125) * x * x * x * x * y * y * z * z - std::sqrt(1302328.125) * x * x * x * x * z * z * z * z + std::sqrt(85317.626953125) * x * x * y * y * y * y * y * y + std::sqrt(349253.173828125) * x * x * y * y * y * y * z * z - std::sqrt(1429312.5) * x * x * z * z * z * z * z * z + std::sqrt(21329.40673828125) * y * y * y * y * y * y * y * y + std::sqrt(349253.173828125) * y * y * y * y * y * y * z * z + std::sqrt(1302328.125) * y * y * y * y * z * z * z * z + std::sqrt(1429312.5) * y * y * z * z * z * z * z * z) + e_3 * (-std::sqrt(2747144.53125) * x * x * x * x * x * x - std::sqrt(2747144.53125) * x * x * x * x * y * y - std::sqrt(29531250.0) * x * x * x * x * z * z + std::sqrt(2747144.53125) * x * x * y * y * y * y - std::sqrt(68229000.0) * x * x * z * z * z * z + std::sqrt(2747144.53125) * y * y * y * y * y * y + std::sqrt(29531250.0) * y * y * y * y * z * z + std::sqrt(68229000.0) * y * y * z * z * z * z) + e_4 * (-std::sqrt(63814078.125) * x * x * x * x - std::sqrt(421954312.5) * x * x * z * z + std::sqrt(63814078.125) * y * y * y * y + std::sqrt(421954312.5) * y * y * z * z) + e_5 * (-std::sqrt(187535250.0) * x * x + std::sqrt(187535250.0) * y * y);
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

        pc_66[k] = e_0 * (std::sqrt(0.720977783203125) * x * x * x * x * x * x * x * x * x * x * x * z + std::sqrt(0.720977783203125) * x * x * x * x * x * x * x * x * x * y * y * z - std::sqrt(307.9376220703125) * x * x * x * x * x * x * x * x * x * z * z * z - std::sqrt(25.9552001953125) * x * x * x * x * x * x * x * y * y * y * y * z + std::sqrt(3737.548828125) * x * x * x * x * x * x * x * z * z * z * z * z - std::sqrt(141.3116455078125) * x * x * x * x * x * y * y * y * y * y * y * z + std::sqrt(11085.75439453125) * x * x * x * x * x * y * y * y * y * z * z * z - std::sqrt(3737.548828125) * x * x * x * x * x * y * y * z * z * z * z * z - std::sqrt(3255.8203125) * x * x * x * x * x * z * z * z * z * z * z * z - std::sqrt(87.23831176757812) * x * x * x * y * y * y * y * y * y * y * y * z + std::sqrt(19708.0078125) * x * x * x * y * y * y * y * y * y * z * z * z - std::sqrt(93438.720703125) * x * x * x * y * y * y * y * z * z * z * z * z + std::sqrt(13023.28125) * x * x * x * y * y * z * z * z * z * z * z * z + std::sqrt(52.5) * x * x * x * z * z * z * z * z * z * z * z * z - std::sqrt(6.488800048828125) * x * y * y * y * y * y * y * y * y * y * y * z + std::sqrt(2771.4385986328125) * x * y * y * y * y * y * y * y * y * z * z * z - std::sqrt(33637.939453125) * x * y * y * y * y * y * y * z * z * z * z * z + std::sqrt(29302.3828125) * x * y * y * y * y * z * z * z * z * z * z * z - std::sqrt(472.5) * x * y * y * z * z * z * z * z * z * z * z * z) + e_1 * (-std::sqrt(25.9552001953125) * x * x * x * x * x * x * x * x * x * z - std::sqrt(415.283203125) * x * x * x * x * x * x * x * z * z * z + std::sqrt(934.38720703125) * x * x * x * x * x * y * y * y * y * z + std::sqrt(415.283203125) * x * x * x * x * x * y * y * z * z * z + std::sqrt(598.0078125) * x * x * x * x * x * z * z * z * z * z + std::sqrt(1661.1328125) * x * x * x * y * y * y * y * y * y * z + std::sqrt(10382.080078125) * x * x * x * y * y * y * y * z * z * z - std::sqrt(2392.03125) * x * x * x * y * y * z * z * z * z * z - std::sqrt(38272.5) * x * x * x * z * z * z * z * z * z * z + std::sqrt(233.5968017578125) * x * y * y * y * y * y * y * y * y * z + std::sqrt(3737.548828125) * x * y * y * y * y * y * y * z * z * z - std::sqrt(5382.0703125) * x * y * y * y * y * z * z * z * z * z + std::sqrt(344452.5) * x * y * y * z * z * z * z * z * z * z) + e_2 * (-std::sqrt(14950.1953125) * x * x * x * x * x * x * x * z + std::sqrt(14950.1953125) * x * x * x * x * x * y * y * z - std::sqrt(6644.53125) * x * x * x * x * x * z * z * z + std::sqrt(373754.8828125) * x * x * x * y * y * y * y * z + std::sqrt(26578.125) * x * x * x * y * y * z * z * z - std::sqrt(3827250.0) * x * x * x * z * z * z * z * z + std::sqrt(134551.7578125) * x * y * y * y * y * y * y * z + std::sqrt(59800.78125) * x * y * y * y * y * z * z * z + std::sqrt(34445250.0) * x * y * y * z * z * z * z * z) + e_3 * (-std::sqrt(803988.28125) * x * x * x * x * x * z + std::sqrt(3215953.125) * x * x * x * y * y * z - std::sqrt(45407250.0) * x * x * x * z * z * z + std::sqrt(7235894.53125) * x * y * y * y * y * z + std::sqrt(408665250.0) * x * y * y * z * z * z) + e_4 * (-std::sqrt(46883812.5) * x * x * x * z + std::sqrt(421954312.5) * x * y * y * z);
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

        pc_67[k] = e_0 * (std::sqrt(0.0240325927734375) * x * x * x * x * x * x * x * x * x * x * x * x - std::sqrt(0.09613037109375) * x * x * x * x * x * x * x * x * x * x * y * y - std::sqrt(18.841552734375) * x * x * x * x * x * x * x * x * x * x * z * z - std::sqrt(6.9454193115234375) * x * x * x * x * x * x * x * x * y * y * y * y + std::sqrt(169.573974609375) * x * x * x * x * x * x * x * x * y * y * z * z + std::sqrt(1000.140380859375) * x * x * x * x * x * x * x * x * z * z * z * z - std::sqrt(18.841552734375) * x * x * x * x * x * x * y * y * y * y * y * y + std::sqrt(3692.9443359375) * x * x * x * x * x * x * y * y * y * y * z * z - std::sqrt(16002.24609375) * x * x * x * x * x * x * y * y * z * z * z * z - std::sqrt(1421.4375) * x * x * x * x * x * x * z * z * z * z * z * z - std::sqrt(6.9454193115234375) * x * x * x * x * y * y * y * y * y * y * y * y + std::sqrt(3692.9443359375) * x * x * x * x * y * y * y * y * y * y * z * z - std::sqrt(100014.0380859375) * x * x * x * x * y * y * y * y * z * z * z * z + std::sqrt(35535.9375) * x * x * x * x * y * y * z * z * z * z * z * z + std::sqrt(24.609375) * x * x * x * x * z * z * z * z * z * z * z * z - std::sqrt(0.09613037109375) * x * x * y * y * y * y * y * y * y * y * y * y + std::sqrt(169.573974609375) * x * x * y * y * y * y * y * y * y * y * z * z - std::sqrt(16002.24609375) * x * x * y * y * y * y * y * y * z * z * z * z + std::sqrt(35535.9375) * x * x * y * y * y * y * z * z * z * z * z * z - std::sqrt(885.9375) * x * x * y * y * z * z * z * z * z * z * z * z + std::sqrt(0.0240325927734375) * y * y * y * y * y * y * y * y * y * y * y * y - std::sqrt(18.841552734375) * y * y * y * y * y * y * y * y * y * y * z * z + std::sqrt(1000.140380859375) * y * y * y * y * y * y * y * y * z * z * z * z - std::sqrt(1421.4375) * y * y * y * y * y * y * z * z * z * z * z * z + std::sqrt(24.609375) * y * y * y * y * z * z * z * z * z * z * z * z) + e_1 * (std::sqrt(31.146240234375) * x * x * x * x * x * x * x * x * x * x - std::sqrt(280.316162109375) * x * x * x * x * x * x * x * x * y * y + std::sqrt(55.37109375) * x * x * x * x * x * x * x * x * z * z - std::sqrt(6104.6630859375) * x * x * x * x * x * x * y * y * y * y - std::sqrt(885.9375) * x * x * x * x * x * x * y * y * z * z + std::sqrt(37430.859375) * x * x * x * x * x * x * z * z * z * z - std::sqrt(6104.6630859375) * x * x * x * x * y * y * y * y * y * y - std::sqrt(5537.109375) * x * x * x * x * y * y * y * y * z * z - std::sqrt(935771.484375) * x * x * x * x * y * y * z * z * z * z - std::sqrt(56700.0) * x * x * x * x * z * z * z * z * z * z - std::sqrt(280.316162109375) * x * x * y * y * y * y * y * y * y * y - std::sqrt(885.9375) * x * x * y * y * y * y * y * y * z * z - std::sqrt(935771.484375) * x * x * y * y * y * y * z * z * z * z + std::sqrt(2041200.0) * x * x * y * y * z * z * z * z * z * z + std::sqrt(31.146240234375) * y * y * y * y * y * y * y * y * y * y + std::sqrt(55.37109375) * y * y * y * y * y * y * y * y * z * z + std::sqrt(37430.859375) * y * y * y * y * y * y * z * z * z * z - std::sqrt(56700.0) * y * y * y * y * z * z * z * z * z * z) + e_2 * (std::sqrt(14621.429443359375) * x * x * x * x * x * x * x * x - std::sqrt(233942.87109375) * x * x * x * x * x * x * y * y + std::sqrt(448505.859375) * x * x * x * x * x * x * z * z - std::sqrt(1462142.9443359375) * x * x * x * x * y * y * y * y - std::sqrt(11212646.484375) * x * x * x * x * y * y * z * z - std::sqrt(669990.234375) * x * x * x * x * z * z * z * z - std::sqrt(233942.87109375) * x * x * y * y * y * y * y * y - std::sqrt(11212646.484375) * x * x * y * y * y * y * z * z + std::sqrt(24119648.4375) * x * x * y * y * z * z * z * z + std::sqrt(14621.429443359375) * y * y * y * y * y * y * y * y + std::sqrt(448505.859375) * y * y * y * y * y * y * z * z - std::sqrt(669990.234375) * y * y * y * y * z * z * z * z) + e_3 * (std::sqrt(1417500.0) * x * x * x * x * x * x - std::sqrt(35437500.0) * x * x * x * x * y * y + std::sqrt(354375.0) * x * x * x * x * z * z - std::sqrt(35437500.0) * x * x * y * y * y * y - std::sqrt(12757500.0) * x * x * y * y * z * z + std::sqrt(1417500.0) * y * y * y * y * y * y + std::sqrt(354375.0) * y * y * y * y * z * z) + e_4 * (std::sqrt(9767460.9375) * x * x * x * x - std::sqrt(351628593.75) * x * x * y * y + std::sqrt(9767460.9375) * y * y * y * y);
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

        pc_68[k] = e_0 * (-std::sqrt(0.528717041015625) * x * x * x * x * x * x * x * x * x * x * x * z + std::sqrt(25.907135009765625) * x * x * x * x * x * x * x * x * x * y * y * z + std::sqrt(171.3043212890625) * x * x * x * x * x * x * x * x * x * z * z * z + std::sqrt(255.8990478515625) * x * x * x * x * x * x * x * y * y * y * y * z - std::sqrt(10963.4765625) * x * x * x * x * x * x * x * y * y * z * z * z - std::sqrt(304.541015625) * x * x * x * x * x * x * x * z * z * z * z * z + std::sqrt(103.6285400390625) * x * x * x * x * x * y * y * y * y * y * y * z - std::sqrt(33575.64697265625) * x * x * x * x * x * y * y * y * y * z * z * z + std::sqrt(24667.822265625) * x * x * x * x * x * y * y * z * z * z * z * z + std::sqrt(5.4140625) * x * x * x * x * x * z * z * z * z * z * z * z - std::sqrt(13.217926025390625) * x * x * x * y * y * y * y * y * y * y * y * z + std::sqrt(7613.525390625) * x * x * x * y * y * y * y * z * z * z * z * z - std::sqrt(541.40625) * x * x * x * y * y * z * z * z * z * z * z * z - std::sqrt(13.217926025390625) * x * y * y * y * y * y * y * y * y * y * y * z + std::sqrt(4282.6080322265625) * x * y * y * y * y * y * y * y * y * z * z * z - std::sqrt(7613.525390625) * x * y * y * y * y * y * y * z * z * z * z * z + std::sqrt(135.3515625) * x * y * y * y * y * z * z * z * z * z * z * z) + e_1 * (std::sqrt(19.0338134765625) * x * x * x * x * x * x * x * x * x * z - std::sqrt(1218.1640625) * x * x * x * x * x * x * x * y * y * z + std::sqrt(36849.462890625) * x * x * x * x * x * x * x * z * z * z - std::sqrt(3730.62744140625) * x * x * x * x * x * y * y * y * y * z - std::sqrt(2984806.494140625) * x * x * x * x * x * y * y * z * z * z - std::sqrt(25776.3515625) * x * x * x * x * x * z * z * z * z * z - std::sqrt(921236.572265625) * x * x * x * y * y * y * y * z * z * z + std::sqrt(2577635.15625) * x * x * x * y * y * z * z * z * z * z + std::sqrt(475.8453369140625) * x * y * y * y * y * y * y * y * y * z + std::sqrt(921236.572265625) * x * y * y * y * y * y * y * z * z * z - std::sqrt(644408.7890625) * x * y * y * y * y * z * z * z * z * z) + e_2 * (std::sqrt(121816.40625) * x * x * x * x * x * x * x * z - std::sqrt(9867128.90625) * x * x * x * x * x * y * y * z + std::sqrt(121816.40625) * x * x * x * x * x * z * z * z - std::sqrt(3045410.15625) * x * x * x * y * y * y * y * z - std::sqrt(12181640.625) * x * x * x * y * y * z * z * z + std::sqrt(3045410.15625) * x * y * y * y * y * y * y * z + std::sqrt(3045410.15625) * x * y * y * y * y * z * z * z) + e_3 * (std::sqrt(3045410.15625) * x * x * x * x * x * z - std::sqrt(304541015.625) * x * x * x * y * y * z + std::sqrt(76135253.90625) * x * y * y * y * y * z);
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

        pc_69[k] = e_0 * (-std::sqrt(0.04405975341796875) * x * x * x * x * x * x * x * x * x * x * x * x + std::sqrt(6.3446044921875) * x * x * x * x * x * x * x * x * x * x * y * y + std::sqrt(14.275360107421875) * x * x * x * x * x * x * x * x * x * x * z * z + std::sqrt(32.11956024169922) * x * x * x * x * x * x * x * x * y * y * y * y - std::sqrt(2412.535858154297) * x * x * x * x * x * x * x * x * y * y * z * z - std::sqrt(25.37841796875) * x * x * x * x * x * x * x * x * z * z * z * z - std::sqrt(2797.9705810546875) * x * x * x * x * x * x * y * y * y * y * z * z + std::sqrt(4974.169921875) * x * x * x * x * x * x * y * y * z * z * z * z + std::sqrt(0.451171875) * x * x * x * x * x * x * z * z * z * z * z * z - std::sqrt(32.11956024169922) * x * x * x * x * y * y * y * y * y * y * y * y + std::sqrt(2797.9705810546875) * x * x * x * x * y * y * y * y * y * y * z * z - std::sqrt(101.513671875) * x * x * x * x * y * y * z * z * z * z * z * z - std::sqrt(6.3446044921875) * x * x * y * y * y * y * y * y * y * y * y * y + std::sqrt(2412.535858154297) * x * x * y * y * y * y * y * y * y * y * z * z - std::sqrt(4974.169921875) * x * x * y * y * y * y * y * y * z * z * z * z + std::sqrt(101.513671875) * x * x * y * y * y * y * z * z * z * z * z * z + std::sqrt(0.04405975341796875) * y * y * y * y * y * y * y * y * y * y * y * y - std::sqrt(14.275360107421875) * y * y * y * y * y * y * y * y * y * y * z * z + std::sqrt(25.37841796875) * y * y * y * y * y * y * y * y * z * z * z * z - std::sqrt(0.451171875) * y * y * y * y * y * y * z * z * z * z * z * z) + e_1 * (-std::sqrt(57.1014404296875) * x * x * x * x * x * x * x * x * x * x + std::sqrt(9650.143432617188) * x * x * x * x * x * x * x * x * y * y + std::sqrt(8222.607421875) * x * x * x * x * x * x * x * x * z * z + std::sqrt(11191.88232421875) * x * x * x * x * x * x * y * y * y * y - std::sqrt(1611631.0546875) * x * x * x * x * x * x * y * y * z * z - std::sqrt(3654.4921875) * x * x * x * x * x * x * z * z * z * z - std::sqrt(11191.88232421875) * x * x * x * x * y * y * y * y * y * y + std::sqrt(822260.7421875) * x * x * x * x * y * y * z * z * z * z - std::sqrt(9650.143432617188) * x * x * y * y * y * y * y * y * y * y + std::sqrt(1611631.0546875) * x * x * y * y * y * y * y * y * z * z - std::sqrt(822260.7421875) * x * x * y * y * y * y * z * z * z * z + std::sqrt(57.1014404296875) * y * y * y * y * y * y * y * y * y * y - std::sqrt(8222.607421875) * y * y * y * y * y * y * y * y * z * z + std::sqrt(3654.4921875) * y * y * y * y * y * y * z * z * z * z) + e_2 * (-std::sqrt(5710.14404296875) * x * x * x * x * x * x * x * x + std::sqrt(1119188.232421875) * x * x * x * x * x * x * y * y + std::sqrt(205565.185546875) * x * x * x * x * x * x * z * z - std::sqrt(46252166.748046875) * x * x * x * x * y * y * z * z - std::sqrt(1119188.232421875) * x * x * y * y * y * y * y * y + std::sqrt(46252166.748046875) * x * x * y * y * y * y * z * z + std::sqrt(5710.14404296875) * y * y * y * y * y * y * y * y - std::sqrt(205565.185546875) * y * y * y * y * y * y * z * z) + e_3 * (-std::sqrt(40605.46875) * x * x * x * x * x * x + std::sqrt(9136230.46875) * x * x * x * x * y * y - std::sqrt(9136230.46875) * x * x * y * y * y * y + std::sqrt(40605.46875) * y * y * y * y * y * y);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ab_x, ab_y, ab_z : simd::cache_line_size())
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
        const auto e_6 = pe_6[k];

        pc_70[k] = e_0 * (8.203125 * x * x * x * x * x * x * x * x * x * x * z * z + 32.8125 * x * x * x * x * x * x * x * x * y * y * z * z - 65.625 * x * x * x * x * x * x * x * x * z * z * z * z + 49.21875 * x * x * x * x * x * x * y * y * y * y * z * z - 196.875 * x * x * x * x * x * x * y * y * z * z * z * z + 157.5 * x * x * x * x * x * x * z * z * z * z * z * z + 32.8125 * x * x * x * x * y * y * y * y * y * y * z * z - 196.875 * x * x * x * x * y * y * y * y * z * z * z * z + 315.0 * x * x * x * x * y * y * z * z * z * z * z * z - 105.0 * x * x * x * x * z * z * z * z * z * z * z * z + 8.203125 * x * x * y * y * y * y * y * y * y * y * z * z - 65.625 * x * x * y * y * y * y * y * y * z * z * z * z + 157.5 * x * x * y * y * y * y * z * z * z * z * z * z - 105.0 * x * x * y * y * z * z * z * z * z * z * z * z + 21.0 * x * x * z * z * z * z * z * z * z * z * z * z) + e_1 * (8.203125 * x * x * x * x * x * x * x * x * x * x + 32.8125 * x * x * x * x * x * x * x * x * y * y + 8.203125 * x * x * x * x * x * x * x * x * z * z + 49.21875 * x * x * x * x * x * x * y * y * y * y + 32.8125 * x * x * x * x * x * x * y * y * z * z + 328.125 * x * x * x * x * x * x * z * z * z * z + 32.8125 * x * x * x * x * y * y * y * y * y * y + 49.21875 * x * x * x * x * y * y * y * y * z * z + 590.625 * x * x * x * x * y * y * z * z * z * z - 262.5 * x * x * x * x * z * z * z * z * z * z + 8.203125 * x * x * y * y * y * y * y * y * y * y + 32.8125 * x * x * y * y * y * y * y * y * z * z + 196.875 * x * x * y * y * y * y * z * z * z * z - 105.0 * x * x * y * y * z * z * z * z * z * z + 210.0 * x * x * z * z * z * z * z * z * z * z + 8.203125 * y * y * y * y * y * y * y * y * z * z - 65.625 * y * y * y * y * y * y * z * z * z * z + 157.5 * y * y * y * y * z * z * z * z * z * z - 105.0 * y * y * z * z * z * z * z * z * z * z + 21.0 * z * z * z * z * z * z * z * z * z * z) + e_2 * (205.078125 * x * x * x * x * x * x * x * x + 623.4375 * x * x * x * x * x * x * y * y + 1115.625 * x * x * x * x * x * x * z * z + 639.84375 * x * x * x * x * y * y * y * y + 2165.625 * x * x * x * x * y * y * z * z + 787.5 * x * x * x * x * z * z * z * z + 229.6875 * x * x * y * y * y * y * y * y + 984.375 * x * x * y * y * y * y * z * z + 1575.0 * x * x * y * y * z * z * z * z + 2100.0 * x * x * z * z * z * z * z * z + 8.203125 * y * y * y * y * y * y * y * y - 65.625 * y * y * y * y * y * y * z * z + 787.5 * y * y * y * y * z * z * z * z - 1050.0 * y * y * z * z * z * z * z * z + 525.0 * z * z * z * z * z * z * z * z) + e_3 * (2493.75 * x * x * x * x * x * x + 5118.75 * x * x * x * x * y * y + 7875.0 * x * x * x * x * z * z + 2756.25 * x * x * y * y * y * y + 9450.0 * x * x * y * y * z * z + 12600.0 * x * x * z * z * z * z + 131.25 * y * y * y * y * y * y + 1575.0 * y * y * y * y * z * z - 3150.0 * y * y * z * z * z * z + 5250.0 * z * z * z * z * z * z) + e_4 * (12600.0 * x * x * x * x + 14175.0 * x * x * y * y + 33075.0 * x * x * z * z + 1575.0 * y * y * y * y + 22050.0 * z * z * z * z) + e_5 * (24570.0 * x * x + 4725.0 * y * y + 33075.0 * z * z) + e_6 * (10395.0);
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

        pc_71[k] = e_0 * (std::sqrt(1.682281494140625) * x * x * x * x * x * x * x * x * x * x * x * z + std::sqrt(15.140533447265625) * x * x * x * x * x * x * x * x * x * y * y * z - std::sqrt(672.91259765625) * x * x * x * x * x * x * x * x * x * z * z * z + std::sqrt(6.7291259765625) * x * x * x * x * x * x * x * y * y * y * y * z - std::sqrt(2691.650390625) * x * x * x * x * x * x * x * y * y * z * z * z + std::sqrt(11201.572265625) * x * x * x * x * x * x * x * z * z * z * z * z - std::sqrt(6.7291259765625) * x * x * x * x * x * y * y * y * y * y * y * z + std::sqrt(11201.572265625) * x * x * x * x * x * y * y * z * z * z * z * z - std::sqrt(13505.625) * x * x * x * x * x * z * z * z * z * z * z * z - std::sqrt(15.140533447265625) * x * x * x * y * y * y * y * y * y * y * y * z + std::sqrt(2691.650390625) * x * x * x * y * y * y * y * y * y * z * z * z - std::sqrt(11201.572265625) * x * x * x * y * y * y * y * z * z * z * z * z + std::sqrt(1102.5) * x * x * x * z * z * z * z * z * z * z * z * z - std::sqrt(1.682281494140625) * x * y * y * y * y * y * y * y * y * y * y * z + std::sqrt(672.91259765625) * x * y * y * y * y * y * y * y * y * z * z * z - std::sqrt(11201.572265625) * x * y * y * y * y * y * y * z * z * z * z * z + std::sqrt(13505.625) * x * y * y * y * y * z * z * z * z * z * z * z - std::sqrt(1102.5) * x * y * y * z * z * z * z * z * z * z * z * z) + e_1 * (-std::sqrt(6.7291259765625) * x * x * x * x * x * x * x * x * x * z + std::sqrt(5275.634765625) * x * x * x * x * x * x * x * z * z * z + std::sqrt(242.24853515625) * x * x * x * x * x * y * y * y * y * z - std::sqrt(968.994140625) * x * x * x * x * x * y * y * z * z * z - std::sqrt(12558.1640625) * x * x * x * x * x * z * z * z * z * z + std::sqrt(430.6640625) * x * x * x * y * y * y * y * y * y * z - std::sqrt(78488.525390625) * x * x * x * y * y * y * y * z * z * z + std::sqrt(179225.15625) * x * x * x * y * y * z * z * z * z * z + std::sqrt(1102.5) * x * x * x * z * z * z * z * z * z * z + std::sqrt(60.5621337890625) * x * y * y * y * y * y * y * y * y * z - std::sqrt(31115.478515625) * x * y * y * y * y * y * y * z * z * z + std::sqrt(286667.2265625) * x * y * y * y * y * z * z * z * z * z - std::sqrt(248062.5) * x * y * y * z * z * z * z * z * z * z + std::sqrt(4410.0) * x * z * z * z * z * z * z * z * z * z) + e_2 * (std::sqrt(3875.9765625) * x * x * x * x * x * x * x * z + std::sqrt(34883.7890625) * x * x * x * x * x * z * z * z - std::sqrt(34883.7890625) * x * x * x * y * y * y * y * z + std::sqrt(1255816.40625) * x * x * x * y * y * z * z * z - std::sqrt(15503.90625) * x * y * y * y * y * y * y * z + std::sqrt(872094.7265625) * x * y * y * y * y * z * z * z - std::sqrt(8930250.0) * x * y * y * z * z * z * z * z + std::sqrt(992250.0) * x * z * z * z * z * z * z * z) + e_3 * (std::sqrt(387597.65625) * x * x * x * x * x * z + std::sqrt(558140.625) * x * x * x * y * y * z + std::sqrt(992250.0) * x * x * x * z * z * z + std::sqrt(15503.90625) * x * y * y * y * y * z - std::sqrt(48620250.0) * x * y * y * z * z * z + std::sqrt(35721000.0) * x * z * z * z * z * z) + e_4 * (std::sqrt(6201562.5) * x * x * x * z - std::sqrt(20093062.5) * x * y * y * z + std::sqrt(194481000.0) * x * z * z * z) + e_5 * (std::sqrt(80372250.0) * x * z);
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

        pc_72[k] = e_0 * (-std::sqrt(60.5621337890625) * x * x * x * x * x * x * x * x * x * x * z * z + std::sqrt(2691.650390625) * x * x * x * x * x * x * x * x * z * z * z * z + std::sqrt(2180.23681640625) * x * x * x * x * x * x * y * y * y * y * z * z - std::sqrt(2691.650390625) * x * x * x * x * x * x * y * y * z * z * z * z - std::sqrt(9112.8515625) * x * x * x * x * x * x * z * z * z * z * z * z + std::sqrt(3875.9765625) * x * x * x * x * y * y * y * y * y * y * z * z - std::sqrt(67291.259765625) * x * x * x * x * y * y * y * y * z * z * z * z + std::sqrt(36451.40625) * x * x * x * x * y * y * z * z * z * z * z * z + std::sqrt(1102.5) * x * x * x * x * z * z * z * z * z * z * z * z + std::sqrt(545.0592041015625) * x * x * y * y * y * y * y * y * y * y * z * z - std::sqrt(24224.853515625) * x * x * y * y * y * y * y * y * z * z * z * z + std::sqrt(82015.6640625) * x * x * y * y * y * y * z * z * z * z * z * z - std::sqrt(9922.5) * x * x * y * y * z * z * z * z * z * z * z * z) + e_1 * (-std::sqrt(60.5621337890625) * x * x * x * x * x * x * x * x * x * x - std::sqrt(1514.0533447265625) * x * x * x * x * x * x * x * x * z * z + std::sqrt(2180.23681640625) * x * x * x * x * x * x * y * y * y * y - std::sqrt(968.994140625) * x * x * x * x * x * x * y * y * z * z - std::sqrt(968.994140625) * x * x * x * x * x * x * z * z * z * z + std::sqrt(3875.9765625) * x * x * x * x * y * y * y * y * y * y + std::sqrt(6056.21337890625) * x * x * x * x * y * y * y * y * z * z + std::sqrt(280039.306640625) * x * x * x * x * y * y * z * z * z * z - std::sqrt(96899.4140625) * x * x * x * x * z * z * z * z * z * z + std::sqrt(545.0592041015625) * x * x * y * y * y * y * y * y * y * y + std::sqrt(8720.947265625) * x * x * y * y * y * y * y * y * z * z + std::sqrt(163760.009765625) * x * x * y * y * y * y * z * z * z * z + std::sqrt(5581.40625) * x * x * y * y * z * z * z * z * z * z + std::sqrt(9922.5) * x * x * z * z * z * z * z * z * z * z + std::sqrt(545.0592041015625) * y * y * y * y * y * y * y * y * z * z - std::sqrt(24224.853515625) * y * y * y * y * y * y * z * z * z * z + std::sqrt(82015.6640625) * y * y * y * y * z * z * z * z * z * z - std::sqrt(9922.5) * y * y * z * z * z * z * z * z * z * z) + e_2 * (-std::sqrt(37851.33361816406) * x * x * x * x * x * x * x * x + std::sqrt(15503.90625) * x * x * x * x * x * x * y * y - std::sqrt(427326.416015625) * x * x * x * x * x * x * z * z + std::sqrt(732801.8188476562) * x * x * x * x * y * y * y * y + std::sqrt(2520353.759765625) * x * x * x * x * y * y * z * z - std::sqrt(5306211.9140625) * x * x * x * x * z * z * z * z + std::sqrt(313954.1015625) * x * x * y * y * y * y * y * y + std::sqrt(4613381.103515625) * x * x * y * y * y * y * z * z + std::sqrt(11302347.65625) * x * x * y * y * z * z * z * z + std::sqrt(248062.5) * x * x * z * z * z * z * z * z + std::sqrt(545.0592041015625) * y * y * y * y * y * y * y * y - std::sqrt(8720.947265625) * y * y * y * y * y * y * z * z + std::sqrt(1399227.5390625) * y * y * y * y * z * z * z * z - std::sqrt(248062.5) * y * y * z * z * z * z * z * z) + e_3 * (-std::sqrt(3969000.0) * x * x * x * x * x * x + std::sqrt(8201566.40625) * x * x * x * x * y * y - std::sqrt(53969097.65625) * x * x * x * x * z * z + std::sqrt(27348890.625) * x * x * y * y * y * y + std::sqrt(201488765.625) * x * x * y * y * z * z - std::sqrt(992250.0) * x * x * z * z * z * z + std::sqrt(139535.15625) * y * y * y * y * y * y + std::sqrt(6837222.65625) * y * y * y * y * z * z + std::sqrt(992250.0) * y * y * z * z * z * z) + e_4 * (-std::sqrt(73814097.65625) * x * x * x * x + std::sqrt(246140015.625) * x * x * y * y - std::sqrt(80372250.0) * x * x * z * z + std::sqrt(11302347.65625) * y * y * y * y + std::sqrt(80372250.0) * y * y * z * z) + e_5 * (-std::sqrt(80372250.0) * x * x + std::sqrt(80372250.0) * y * y);
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

        pc_73[k] = e_0 * (-std::sqrt(2.01873779296875) * x * x * x * x * x * x * x * x * x * x * x * z + std::sqrt(18.16864013671875) * x * x * x * x * x * x * x * x * x * y * y * z + std::sqrt(395.672607421875) * x * x * x * x * x * x * x * x * x * z * z * z + std::sqrt(395.672607421875) * x * x * x * x * x * x * x * y * y * y * y * z - std::sqrt(6330.76171875) * x * x * x * x * x * x * x * y * y * z * z * z - std::sqrt(3493.546875) * x * x * x * x * x * x * x * z * z * z * z * z + std::sqrt(395.672607421875) * x * x * x * x * x * y * y * y * y * y * y * z - std::sqrt(39567.2607421875) * x * x * x * x * x * y * y * y * y * z * z * z + std::sqrt(87338.671875) * x * x * x * x * x * y * y * z * z * z * z * z + std::sqrt(516.796875) * x * x * x * x * x * z * z * z * z * z * z * z + std::sqrt(18.16864013671875) * x * x * x * y * y * y * y * y * y * y * y * z - std::sqrt(6330.76171875) * x * x * x * y * y * y * y * y * y * z * z * z + std::sqrt(87338.671875) * x * x * x * y * y * y * y * z * z * z * z * z - std::sqrt(18604.6875) * x * x * x * y * y * z * z * z * z * z * z * z - std::sqrt(2.01873779296875) * x * y * y * y * y * y * y * y * y * y * y * z + std::sqrt(395.672607421875) * x * y * y * y * y * y * y * y * y * z * z * z - std::sqrt(3493.546875) * x * y * y * y * y * y * y * z * z * z * z * z + std::sqrt(516.796875) * x * y * y * y * y * z * z * z * z * z * z * z) + e_1 * (-std::sqrt(201.873779296875) * x * x * x * x * x * x * x * x * x * z + std::sqrt(1162.79296875) * x * x * x * x * x * x * x * y * y * z + std::sqrt(2067.1875) * x * x * x * x * x * x * x * z * z * z + std::sqrt(14244.2138671875) * x * x * x * x * x * y * y * y * y * z + std::sqrt(8268.75) * x * x * x * x * x * y * y * z * z * z - std::sqrt(219307.921875) * x * x * x * x * x * z * z * z * z * z + std::sqrt(6330.76171875) * x * x * x * y * y * y * y * y * y * z - std::sqrt(51679.6875) * x * x * x * y * y * y * y * z * z * z + std::sqrt(3474942.1875) * x * x * x * y * y * z * z * z * z * z + std::sqrt(8268.75) * x * x * x * z * z * z * z * z * z * z + std::sqrt(72.674560546875) * x * y * y * y * y * y * y * y * y * z - std::sqrt(74418.75) * x * y * y * y * y * y * y * z * z * z + std::sqrt(227907.421875) * x * y * y * y * y * z * z * z * z * z - std::sqrt(74418.75) * x * y * y * z * z * z * z * z * z * z) + e_2 * (-std::sqrt(29069.82421875) * x * x * x * x * x * x * x * z + std::sqrt(726745.60546875) * x * x * x * x * x * y * y * z - std::sqrt(3307500.0) * x * x * x * x * x * z * z * z + std::sqrt(726745.60546875) * x * x * x * y * y * y * y * z + std::sqrt(82687500.0) * x * x * x * y * y * z * z * z - std::sqrt(206718.75) * x * x * x * z * z * z * z * z - std::sqrt(29069.82421875) * x * y * y * y * y * y * y * z + std::sqrt(1860468.75) * x * y * y * z * z * z * z * z) + e_3 * (-std::sqrt(7441875.0) * x * x * x * x * x * z + std::sqrt(186046875.0) * x * x * x * y * y * z - std::sqrt(20671875.0) * x * x * x * z * z * z + std::sqrt(186046875.0) * x * y * y * z * z * z) + e_4 * (-std::sqrt(46511718.75) * x * x * x * z + std::sqrt(418605468.75) * x * y * y * z);
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

        pc_74[k] = e_0 * (std::sqrt(44.4122314453125) * x * x * x * x * x * x * x * x * x * x * z * z - std::sqrt(2842.3828125) * x * x * x * x * x * x * x * x * y * y * z * z - std::sqrt(710.595703125) * x * x * x * x * x * x * x * x * z * z * z * z - std::sqrt(8704.79736328125) * x * x * x * x * x * x * y * y * y * y * z * z + std::sqrt(57558.251953125) * x * x * x * x * x * x * y * y * z * z * z * z + std::sqrt(113.6953125) * x * x * x * x * x * x * z * z * z * z * z * z + std::sqrt(17764.892578125) * x * x * x * x * y * y * y * y * z * z * z * z - std::sqrt(11369.53125) * x * x * x * x * y * y * z * z * z * z * z * z + std::sqrt(1110.3057861328125) * x * x * y * y * y * y * y * y * y * y * z * z - std::sqrt(17764.892578125) * x * x * y * y * y * y * y * y * z * z * z * z + std::sqrt(2842.3828125) * x * x * y * y * y * y * z * z * z * z * z * z) + e_1 * (std::sqrt(44.4122314453125) * x * x * x * x * x * x * x * x * x * x - std::sqrt(2842.3828125) * x * x * x * x * x * x * x * x * y * y + std::sqrt(7505.6671142578125) * x * x * x * x * x * x * x * x * z * z - std::sqrt(8704.79736328125) * x * x * x * x * x * x * y * y * y * y - std::sqrt(375905.126953125) * x * x * x * x * x * x * y * y * z * z - std::sqrt(120090.673828125) * x * x * x * x * x * x * z * z * z * z - std::sqrt(359739.07470703125) * x * x * x * x * y * y * y * y * z * z + std::sqrt(7834317.626953125) * x * x * x * x * y * y * z * z * z * z + std::sqrt(2842.3828125) * x * x * x * x * z * z * z * z * z * z + std::sqrt(1110.3057861328125) * x * x * y * y * y * y * y * y * y * y + std::sqrt(17764.892578125) * x * x * y * y * y * y * y * y * z * z - std::sqrt(159884.033203125) * x * x * y * y * y * y * z * z * z * z - std::sqrt(102325.78125) * x * x * y * y * z * z * z * z * z * z + std::sqrt(1110.3057861328125) * y * y * y * y * y * y * y * y * z * z - std::sqrt(17764.892578125) * y * y * y * y * y * y * z * z * z * z + std::sqrt(2842.3828125) * y * y * y * y * z * z * z * z * z * z) + e_2 * (std::sqrt(27757.644653320312) * x * x * x * x * x * x * x * x - std::sqrt(1776489.2578125) * x * x * x * x * x * x * y * y - std::sqrt(17764.892578125) * x * x * x * x * x * x * z * z - std::sqrt(999275.2075195312) * x * x * x * x * y * y * y * y + std::sqrt(3997100.830078125) * x * x * x * x * y * y * z * z - std::sqrt(639536.1328125) * x * x * x * x * z * z * z * z + std::sqrt(284238.28125) * x * x * y * y * y * y * y * y - std::sqrt(3997100.830078125) * x * x * y * y * y * y * z * z + std::sqrt(23023300.78125) * x * x * y * y * z * z * z * z + std::sqrt(1110.3057861328125) * y * y * y * y * y * y * y * y + std::sqrt(17764.892578125) * y * y * y * y * y * y * z * z - std::sqrt(639536.1328125) * y * y * y * y * z * z * z * z) + e_3 * (std::sqrt(1136953.125) * x * x * x * x * x * x - std::sqrt(63953613.28125) * x * x * x * x * y * y - std::sqrt(2558144.53125) * x * x * x * x * z * z + std::sqrt(92093203.125) * x * x * y * y * z * z + std::sqrt(284238.28125) * y * y * y * y * y * y - std::sqrt(2558144.53125) * y * y * y * y * z * z) + e_4 * (std::sqrt(2558144.53125) * x * x * x * x - std::sqrt(92093203.125) * x * x * y * y + std::sqrt(2558144.53125) * y * y * y * y);
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

        pc_75[k] = e_0 * (std::sqrt(3.701019287109375) * x * x * x * x * x * x * x * x * x * x * x * z - std::sqrt(625.4722595214844) * x * x * x * x * x * x * x * x * x * y * y * z - std::sqrt(59.21630859375) * x * x * x * x * x * x * x * x * x * z * z * z - std::sqrt(725.3997802734375) * x * x * x * x * x * x * x * y * y * y * y * z + std::sqrt(11606.396484375) * x * x * x * x * x * x * x * y * y * z * z * z + std::sqrt(9.474609375) * x * x * x * x * x * x * x * z * z * z * z * z + std::sqrt(725.3997802734375) * x * x * x * x * x * y * y * y * y * y * y * z - std::sqrt(2131.787109375) * x * x * x * x * x * y * y * z * z * z * z * z + std::sqrt(625.4722595214844) * x * x * x * y * y * y * y * y * y * y * y * z - std::sqrt(11606.396484375) * x * x * x * y * y * y * y * y * y * z * z * z + std::sqrt(2131.787109375) * x * x * x * y * y * y * y * z * z * z * z * z - std::sqrt(3.701019287109375) * x * y * y * y * y * y * y * y * y * y * y * z + std::sqrt(59.21630859375) * x * y * y * y * y * y * y * y * y * z * z * z - std::sqrt(9.474609375) * x * y * y * y * y * y * y * z * z * z * z * z) + e_1 * (std::sqrt(3330.9173583984375) * x * x * x * x * x * x * x * x * x * z - std::sqrt(545737.5) * x * x * x * x * x * x * x * y * y * z - std::sqrt(19186.083984375) * x * x * x * x * x * x * x * z * z * z - std::sqrt(26114.39208984375) * x * x * x * x * x * y * y * y * y * z + std::sqrt(3242448.193359375) * x * x * x * x * x * y * y * z * z * z + std::sqrt(341.0859375) * x * x * x * x * x * z * z * z * z * z + std::sqrt(417830.2734375) * x * x * x * y * y * y * y * y * y * z - std::sqrt(1332366.943359375) * x * x * x * y * y * y * y * z * z * z - std::sqrt(34108.59375) * x * x * x * y * y * z * z * z * z * z + std::sqrt(133.2366943359375) * x * y * y * y * y * y * y * y * y * z - std::sqrt(19186.083984375) * x * y * y * y * y * y * y * z * z * z + std::sqrt(8527.1484375) * x * y * y * y * y * z * z * z * z * z) + e_2 * (std::sqrt(213178.7109375) * x * x * x * x * x * x * x * z - std::sqrt(30697734.375) * x * x * x * x * x * y * y * z - std::sqrt(213178.7109375) * x * x * x * x * x * z * z * z + std::sqrt(5329467.7734375) * x * x * x * y * y * y * y * z + std::sqrt(21317871.09375) * x * x * x * y * y * z * z * z + std::sqrt(852714.84375) * x * y * y * y * y * y * y * z - std::sqrt(5329467.7734375) * x * y * y * y * y * z * z * z) + e_3 * (std::sqrt(852714.84375) * x * x * x * x * x * z - std::sqrt(85271484.375) * x * x * x * y * y * z + std::sqrt(21317871.09375) * x * y * y * y * y * z);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ab_x, ab_y, ab_z : simd::cache_line_size())
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
        const auto e_6 = pe_6[k];

        pc_76[k] = e_0 * (0.205078125 * x * x * x * x * x * x * x * x * x * x * x * x + 0.41015625 * x * x * x * x * x * x * x * x * x * x * y * y - 6.5625 * x * x * x * x * x * x * x * x * x * x * z * z - 0.205078125 * x * x * x * x * x * x * x * x * y * y * y * y - 6.5625 * x * x * x * x * x * x * x * x * y * y * z * z + 59.0625 * x * x * x * x * x * x * x * x * z * z * z * z - 0.8203125 * x * x * x * x * x * x * y * y * y * y * y * y + 13.125 * x * x * x * x * x * x * y * y * y * y * z * z - 105.0 * x * x * x * x * x * x * z * z * z * z * z * z - 0.205078125 * x * x * x * x * y * y * y * y * y * y * y * y + 13.125 * x * x * x * x * y * y * y * y * y * y * z * z - 118.125 * x * x * x * x * y * y * y * y * z * z * z * z + 105.0 * x * x * x * x * y * y * z * z * z * z * z * z + 52.5 * x * x * x * x * z * z * z * z * z * z * z * z + 0.41015625 * x * x * y * y * y * y * y * y * y * y * y * y - 6.5625 * x * x * y * y * y * y * y * y * y * y * z * z + 105.0 * x * x * y * y * y * y * z * z * z * z * z * z - 105.0 * x * x * y * y * z * z * z * z * z * z * z * z + 0.205078125 * y * y * y * y * y * y * y * y * y * y * y * y - 6.5625 * y * y * y * y * y * y * y * y * y * y * z * z + 59.0625 * y * y * y * y * y * y * y * y * z * z * z * z - 105.0 * y * y * y * y * y * y * z * z * z * z * z * z + 52.5 * y * y * y * y * z * z * z * z * z * z * z * z) + e_1 * (7.3828125 * x * x * x * x * x * x * x * x * x * x + 10.6640625 * x * x * x * x * x * x * x * x * y * y + 52.5 * x * x * x * x * x * x * x * x * z * z - 4.921875 * x * x * x * x * x * x * y * y * y * y - 105.0 * x * x * x * x * x * x * y * y * z * z + 78.75 * x * x * x * x * x * x * z * z * z * z - 4.921875 * x * x * x * x * y * y * y * y * y * y - 315.0 * x * x * x * x * y * y * y * y * z * z + 866.25 * x * x * x * x * y * y * z * z * z * z + 10.6640625 * x * x * y * y * y * y * y * y * y * y - 105.0 * x * x * y * y * y * y * y * y * z * z + 866.25 * x * x * y * y * y * y * z * z * z * z - 1680.0 * x * x * y * y * z * z * z * z * z * z + 210.0 * x * x * z * z * z * z * z * z * z * z + 7.3828125 * y * y * y * y * y * y * y * y * y * y + 52.5 * y * y * y * y * y * y * y * y * z * z + 78.75 * y * y * y * y * y * y * z * z * z * z + 210.0 * y * y * z * z * z * z * z * z * z * z) + e_2 * (197.6953125 * x * x * x * x * x * x * x * x + 82.03125 * x * x * x * x * x * x * y * y + 918.75 * x * x * x * x * x * x * z * z - 231.328125 * x * x * x * x * y * y * y * y + 866.25 * x * x * x * x * y * y * z * z + 1023.75 * x * x * x * x * z * z * z * z + 82.03125 * x * x * y * y * y * y * y * y + 866.25 * x * x * y * y * y * y * z * z - 7402.5 * x * x * y * y * z * z * z * z + 2100.0 * x * x * z * z * z * z * z * z + 197.6953125 * y * y * y * y * y * y * y * y + 918.75 * y * y * y * y * y * y * z * z + 1023.75 * y * y * y * y * z * z * z * z + 2100.0 * y * y * z * z * z * z * z * z + 210.0 * z * z * z * z * z * z * z * z) + e_3 * (2178.75 * x * x * x * x * x * x + 236.25 * x * x * x * x * y * y + 6930.0 * x * x * x * x * z * z + 236.25 * x * x * y * y * y * y - 11340.0 * x * x * y * y * z * z + 10080.0 * x * x * z * z * z * z + 2178.75 * y * y * y * y * y * y + 6930.0 * y * y * y * y * z * z + 10080.0 * y * y * z * z * z * z + 3360.0 * z * z * z * z * z * z) + e_4 * (9961.875 * x * x * x * x - 2126.25 * x * x * y * y + 22680.0 * x * x * z * z + 9961.875 * y * y * y * y + 22680.0 * y * y * z * z + 17640.0 * z * z * z * z) + e_5 * (16065.0 * x * x + 16065.0 * y * y + 30240.0 * z * z) + e_6 * (10395.0);
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

        pc_77[k] = e_0 * (-1.23046875 * x * x * x * x * x * x * x * x * x * x * x * z + 1.23046875 * x * x * x * x * x * x * x * x * x * y * y * z + 22.96875 * x * x * x * x * x * x * x * x * x * z * z * z + 7.3828125 * x * x * x * x * x * x * x * y * y * y * y * z - 45.9375 * x * x * x * x * x * x * x * y * y * z * z * z - 72.1875 * x * x * x * x * x * x * x * z * z * z * z * z + 2.4609375 * x * x * x * x * x * y * y * y * y * y * y * z - 91.875 * x * x * x * x * x * y * y * y * y * z * z * z + 216.5625 * x * x * x * x * x * y * y * z * z * z * z * z + 52.5 * x * x * x * x * x * z * z * z * z * z * z * z - 6.15234375 * x * x * x * y * y * y * y * y * y * y * y * z + 45.9375 * x * x * x * y * y * y * y * y * y * z * z * z + 72.1875 * x * x * x * y * y * y * y * z * z * z * z * z - 210.0 * x * x * x * y * y * z * z * z * z * z * z * z - 3.69140625 * x * y * y * y * y * y * y * y * y * y * y * z + 68.90625 * x * y * y * y * y * y * y * y * y * z * z * z - 216.5625 * x * y * y * y * y * y * y * z * z * z * z * z + 157.5 * x * y * y * y * y * z * z * z * z * z * z * z) + e_1 * (2.4609375 * x * x * x * x * x * x * x * x * x * z - 49.21875 * x * x * x * x * x * x * x * y * y * z + 59.0625 * x * x * x * x * x * x * x * z * z * z - 83.671875 * x * x * x * x * x * y * y * y * y * z + 649.6875 * x * x * x * x * x * y * y * z * z * z - 196.875 * x * x * x * x * x * z * z * z * z * z - 9.84375 * x * x * x * y * y * y * y * y * y * z + 492.1875 * x * x * x * y * y * y * y * z * z * z - 1811.25 * x * x * x * y * y * z * z * z * z * z + 315.0 * x * x * x * z * z * z * z * z * z * z + 22.1484375 * x * y * y * y * y * y * y * y * y * z - 98.4375 * x * y * y * y * y * y * y * z * z * z + 275.625 * x * y * y * y * y * z * z * z * z * z + 315.0 * x * y * y * z * z * z * z * z * z * z) + e_2 * (108.28125 * x * x * x * x * x * x * x * z + 206.71875 * x * x * x * x * x * y * y * z - 39.375 * x * x * x * x * x * z * z * z + 246.09375 * x * x * x * y * y * y * y * z - 4331.25 * x * x * x * y * y * z * z * z + 1417.5 * x * x * x * z * z * z * z * z + 147.65625 * x * y * y * y * y * y * y * z + 1378.125 * x * y * y * y * y * z * z * z + 1417.5 * x * y * y * z * z * z * z * z + 630.0 * x * z * z * z * z * z * z * z) + e_3 * (787.5 * x * x * x * x * x * z - 3150.0 * x * x * x * y * y * z + 3150.0 * x * x * x * z * z * z + 2362.5 * x * y * y * y * y * z + 3150.0 * x * y * y * z * z * z + 6300.0 * x * z * z * z * z * z) + e_4 * (3543.75 * x * x * x * z + 3543.75 * x * y * y * z + 18900.0 * x * z * z * z) + e_5 * (14175.0 * x * z);
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

        pc_78[k] = e_0 * (-std::sqrt(0.05046844482421875) * x * x * x * x * x * x * x * x * x * x * x * x + std::sqrt(0.8074951171875) * x * x * x * x * x * x * x * x * x * x * y * y + std::sqrt(34.116668701171875) * x * x * x * x * x * x * x * x * x * x * z * z + std::sqrt(6.106681823730469) * x * x * x * x * x * x * x * x * y * y * y * y - std::sqrt(852.9167175292969) * x * x * x * x * x * x * x * x * y * y * z * z - std::sqrt(1563.310546875) * x * x * x * x * x * x * x * x * z * z * z * z - std::sqrt(1228.2000732421875) * x * x * x * x * x * x * y * y * y * y * z * z + std::sqrt(56279.1796875) * x * x * x * x * x * x * y * y * z * z * z * z + std::sqrt(1291.9921875) * x * x * x * x * x * x * z * z * z * z * z * z - std::sqrt(6.106681823730469) * x * x * x * x * y * y * y * y * y * y * y * y + std::sqrt(1228.2000732421875) * x * x * x * x * y * y * y * y * y * y * z * z - std::sqrt(63307.6171875) * x * x * x * x * y * y * z * z * z * z * z * z - std::sqrt(0.8074951171875) * x * x * y * y * y * y * y * y * y * y * y * y + std::sqrt(852.9167175292969) * x * x * y * y * y * y * y * y * y * y * z * z - std::sqrt(56279.1796875) * x * x * y * y * y * y * y * y * z * z * z * z + std::sqrt(63307.6171875) * x * x * y * y * y * y * z * z * z * z * z * z + std::sqrt(0.05046844482421875) * y * y * y * y * y * y * y * y * y * y * y * y - std::sqrt(34.116668701171875) * y * y * y * y * y * y * y * y * y * y * z * z + std::sqrt(1563.310546875) * y * y * y * y * y * y * y * y * z * z * z * z - std::sqrt(1291.9921875) * y * y * y * y * y * y * z * z * z * z * z * z) + e_1 * (-std::sqrt(65.4071044921875) * x * x * x * x * x * x * x * x * x * x + std::sqrt(679.1033935546875) * x * x * x * x * x * x * x * x * y * y - std::sqrt(12.919921875) * x * x * x * x * x * x * x * x * z * z + std::sqrt(1166.02294921875) * x * x * x * x * x * x * y * y * y * y + std::sqrt(156331.0546875) * x * x * x * x * x * x * y * y * z * z - std::sqrt(109354.21875) * x * x * x * x * x * x * z * z * z * z - std::sqrt(1166.02294921875) * x * x * x * x * y * y * y * y * y * y - std::sqrt(46511.71875) * x * x * x * x * y * y * z * z * z * z + std::sqrt(82687.5) * x * x * x * x * z * z * z * z * z * z - std::sqrt(679.1033935546875) * x * x * y * y * y * y * y * y * y * y - std::sqrt(156331.0546875) * x * x * y * y * y * y * y * y * z * z + std::sqrt(46511.71875) * x * x * y * y * y * y * z * z * z * z + std::sqrt(65.4071044921875) * y * y * y * y * y * y * y * y * y * y + std::sqrt(12.919921875) * y * y * y * y * y * y * y * y * z * z + std::sqrt(109354.21875) * y * y * y * y * y * y * z * z * z * z - std::sqrt(82687.5) * y * y * y * y * z * z * z * z * z * z) + e_2 * (-std::sqrt(29150.57373046875) * x * x * x * x * x * x * x * x + std::sqrt(442184.326171875) * x * x * x * x * x * x * y * y - std::sqrt(713502.685546875) * x * x * x * x * x * x * z * z + std::sqrt(5375010.498046875) * x * x * x * x * y * y * z * z - std::sqrt(186046.875) * x * x * x * x * z * z * z * z - std::sqrt(442184.326171875) * x * x * y * y * y * y * y * y - std::sqrt(5375010.498046875) * x * x * y * y * y * y * z * z + std::sqrt(744187.5) * x * x * z * z * z * z * z * z + std::sqrt(29150.57373046875) * y * y * y * y * y * y * y * y + std::sqrt(713502.685546875) * y * y * y * y * y * y * z * z + std::sqrt(186046.875) * y * y * y * y * z * z * z * z - std::sqrt(744187.5) * y * y * z * z * z * z * z * z) + e_3 * (-std::sqrt(2733855.46875) * x * x * x * x * x * x + std::sqrt(16790730.46875) * x * x * x * x * y * y - std::sqrt(18604687.5) * x * x * x * x * z * z - std::sqrt(16790730.46875) * x * x * y * y * y * y + std::sqrt(11907000.0) * x * x * z * z * z * z + std::sqrt(2733855.46875) * y * y * y * y * y * y + std::sqrt(18604687.5) * y * y * y * y * z * z - std::sqrt(11907000.0) * y * y * z * z * z * z) + e_4 * (-std::sqrt(39116355.46875) * x * x * x * x - std::sqrt(1674421.875) * x * x * z * z + std::sqrt(39116355.46875) * y * y * y * y + std::sqrt(1674421.875) * y * y * z * z) + e_5 * (-std::sqrt(60279187.5) * x * x + std::sqrt(60279187.5) * y * y);
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

        pc_79[k] = e_0 * (std::sqrt(1.1103057861328125) * x * x * x * x * x * x * x * x * x * x * x * z - std::sqrt(89.93476867675781) * x * x * x * x * x * x * x * x * x * y * y * z - std::sqrt(284.23828125) * x * x * x * x * x * x * x * x * x * z * z * z - std::sqrt(39.97100830078125) * x * x * x * x * x * x * x * y * y * y * y * z + std::sqrt(28423.828125) * x * x * x * x * x * x * x * y * y * z * z * z + std::sqrt(284.23828125) * x * x * x * x * x * x * x * z * z * z * z * z + std::sqrt(217.61993408203125) * x * x * x * x * x * y * y * y * y * y * y * z - std::sqrt(4547.8125) * x * x * x * x * x * y * y * y * y * z * z * z - std::sqrt(34392.83203125) * x * x * x * x * x * y * y * z * z * z * z * z + std::sqrt(27.757644653320312) * x * x * x * y * y * y * y * y * y * y * y * z - std::sqrt(28423.828125) * x * x * x * y * y * y * y * y * y * z * z * z + std::sqrt(63953.61328125) * x * x * x * y * y * y * y * z * z * z * z * z - std::sqrt(27.757644653320312) * x * y * y * y * y * y * y * y * y * y * y * z + std::sqrt(7105.95703125) * x * y * y * y * y * y * y * y * y * z * z * z - std::sqrt(7105.95703125) * x * y * y * y * y * y * y * z * z * z * z * z) + e_1 * (-std::sqrt(4.44122314453125) * x * x * x * x * x * x * x * x * x * z + std::sqrt(15988.4033203125) * x * x * x * x * x * x * x * y * y * z - std::sqrt(72765.0) * x * x * x * x * x * x * x * z * z * z - std::sqrt(12950.606689453125) * x * x * x * x * x * y * y * y * y * z + std::sqrt(1641760.3125) * x * x * x * x * x * y * y * z * z * z + std::sqrt(28423.828125) * x * x * x * x * x * z * z * z * z * z - std::sqrt(44412.2314453125) * x * x * x * y * y * y * y * y * y * z - std::sqrt(454781.25) * x * x * x * y * y * y * y * z * z * z - std::sqrt(113695.3125) * x * x * x * y * y * z * z * z * z * z + std::sqrt(999.2752075195312) * x * y * y * y * y * y * y * y * y * z + std::sqrt(1023257.8125) * x * y * y * y * y * y * y * z * z * z - std::sqrt(255814.453125) * x * y * y * y * y * z * z * z * z * z) + e_2 * (-std::sqrt(143895.6298828125) * x * x * x * x * x * x * x * z + std::sqrt(8457865.356445312) * x * x * x * x * x * y * y * z - std::sqrt(1819125.0) * x * x * x * x * x * z * z * z - std::sqrt(9992752.075195312) * x * x * x * y * y * y * y * z + std::sqrt(7276500.0) * x * x * x * y * y * z * z * z + std::sqrt(454781.25) * x * x * x * z * z * z * z * z + std::sqrt(2702040.1611328125) * x * y * y * y * y * y * y * z + std::sqrt(16372125.0) * x * y * y * y * y * z * z * z - std::sqrt(4093031.25) * x * y * y * z * z * z * z * z) + e_3 * (-std::sqrt(9209320.3125) * x * x * x * x * x * z + std::sqrt(36837281.25) * x * x * x * y * y * z - std::sqrt(1819125.0) * x * x * x * z * z * z + std::sqrt(82883882.8125) * x * y * y * y * y * z + std::sqrt(16372125.0) * x * y * y * z * z * z) + e_4 * (-std::sqrt(50139632.8125) * x * x * x * z + std::sqrt(451256695.3125) * x * y * y * z);
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

        pc_80[k] = e_0 * (std::sqrt(0.09252548217773438) * x * x * x * x * x * x * x * x * x * x * x * x - std::sqrt(18.134994506835938) * x * x * x * x * x * x * x * x * x * x * y * y - std::sqrt(23.6865234375) * x * x * x * x * x * x * x * x * x * x * z * z - std::sqrt(0.09252548217773438) * x * x * x * x * x * x * x * x * y * y * y * y + std::sqrt(5329.4677734375) * x * x * x * x * x * x * x * x * y * y * z * z + std::sqrt(23.6865234375) * x * x * x * x * x * x * x * x * z * z * z * z + std::sqrt(72.53997802734375) * x * x * x * x * x * x * y * y * y * y * y * y - std::sqrt(4642.55859375) * x * x * x * x * x * x * y * y * y * y * z * z - std::sqrt(6063.75) * x * x * x * x * x * x * y * y * z * z * z * z - std::sqrt(0.09252548217773438) * x * x * x * x * y * y * y * y * y * y * y * y - std::sqrt(4642.55859375) * x * x * x * x * y * y * y * y * y * y * z * z + std::sqrt(21317.87109375) * x * x * x * x * y * y * y * y * z * z * z * z - std::sqrt(18.134994506835938) * x * x * y * y * y * y * y * y * y * y * y * y + std::sqrt(5329.4677734375) * x * x * y * y * y * y * y * y * y * y * z * z - std::sqrt(6063.75) * x * x * y * y * y * y * y * y * z * z * z * z + std::sqrt(0.09252548217773438) * y * y * y * y * y * y * y * y * y * y * y * y - std::sqrt(23.6865234375) * y * y * y * y * y * y * y * y * y * y * z * z + std::sqrt(23.6865234375) * y * y * y * y * y * y * y * y * z * z * z * z) + e_1 * (std::sqrt(119.91302490234375) * x * x * x * x * x * x * x * x * x * x - std::sqrt(14509.476013183594) * x * x * x * x * x * x * x * x * y * y - std::sqrt(13643.4375) * x * x * x * x * x * x * x * x * z * z + std::sqrt(2611.439208984375) * x * x * x * x * x * x * y * y * y * y + std::sqrt(1364343.75) * x * x * x * x * x * x * y * y * z * z + std::sqrt(3410.859375) * x * x * x * x * x * x * z * z * z * z + std::sqrt(2611.439208984375) * x * x * x * x * y * y * y * y * y * y - std::sqrt(1364343.75) * x * x * x * x * y * y * y * y * z * z - std::sqrt(85271.484375) * x * x * x * x * y * y * z * z * z * z - std::sqrt(14509.476013183594) * x * x * y * y * y * y * y * y * y * y + std::sqrt(1364343.75) * x * x * y * y * y * y * y * y * z * z - std::sqrt(85271.484375) * x * x * y * y * y * y * z * z * z * z + std::sqrt(119.91302490234375) * y * y * y * y * y * y * y * y * y * y - std::sqrt(13643.4375) * y * y * y * y * y * y * y * y * z * z + std::sqrt(3410.859375) * y * y * y * y * y * y * z * z * z * z) + e_2 * (std::sqrt(16321.495056152344) * x * x * x * x * x * x * x * x - std::sqrt(900680.0537109375) * x * x * x * x * x * x * y * y - std::sqrt(767443.359375) * x * x * x * x * x * x * z * z + std::sqrt(33309.173583984375) * x * x * x * x * y * y * y * y + std::sqrt(19186083.984375) * x * x * x * x * y * y * z * z + std::sqrt(85271.484375) * x * x * x * x * z * z * z * z - std::sqrt(900680.0537109375) * x * x * y * y * y * y * y * y + std::sqrt(19186083.984375) * x * x * y * y * y * y * z * z - std::sqrt(3069773.4375) * x * x * y * y * z * z * z * z + std::sqrt(16321.495056152344) * y * y * y * y * y * y * y * y - std::sqrt(767443.359375) * y * y * y * y * y * y * z * z + std::sqrt(85271.484375) * y * y * y * y * z * z * z * z) + e_3 * (std::sqrt(341085.9375) * x * x * x * x * x * x - std::sqrt(8527148.4375) * x * x * x * x * y * y - std::sqrt(5457375.0) * x * x * x * x * z * z - std::sqrt(8527148.4375) * x * x * y * y * y * y + std::sqrt(196465500.0) * x * x * y * y * z * z + std::sqrt(341085.9375) * y * y * y * y * y * y - std::sqrt(5457375.0) * y * y * y * y * z * z) + e_4 * (std::sqrt(767443.359375) * x * x * x * x - std::sqrt(27627960.9375) * x * x * y * y + std::sqrt(767443.359375) * y * y * y * y);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ab_x, ab_y, ab_z : simd::cache_line_size())
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
        const auto e_6 = pe_6[k];

        pc_81[k] = e_0 * (7.3828125 * x * x * x * x * x * x * x * x * x * x * z * z - 29.53125 * x * x * x * x * x * x * x * x * y * y * z * z - 39.375 * x * x * x * x * x * x * x * x * z * z * z * z - 14.765625 * x * x * x * x * x * x * y * y * y * y * z * z + 196.875 * x * x * x * x * x * x * y * y * z * z * z * z + 52.5 * x * x * x * x * x * x * z * z * z * z * z * z + 88.59375 * x * x * x * x * y * y * y * y * y * y * z * z - 118.125 * x * x * x * x * y * y * y * y * z * z * z * z - 315.0 * x * x * x * x * y * y * z * z * z * z * z * z + 66.4453125 * x * x * y * y * y * y * y * y * y * y * z * z - 354.375 * x * x * y * y * y * y * y * y * z * z * z * z + 472.5 * x * x * y * y * y * y * z * z * z * z * z * z) + e_1 * (7.3828125 * x * x * x * x * x * x * x * x * x * x - 29.53125 * x * x * x * x * x * x * x * x * y * y + 66.4453125 * x * x * x * x * x * x * x * x * z * z - 14.765625 * x * x * x * x * x * x * y * y * y * y + 265.78125 * x * x * x * x * x * x * y * y * z * z - 118.125 * x * x * x * x * x * x * z * z * z * z + 88.59375 * x * x * x * x * y * y * y * y * y * y + 398.671875 * x * x * x * x * y * y * y * y * z * z - 2480.625 * x * x * x * x * y * y * z * z * z * z + 472.5 * x * x * x * x * z * z * z * z * z * z + 66.4453125 * x * x * y * y * y * y * y * y * y * y + 265.78125 * x * x * y * y * y * y * y * y * z * z + 1063.125 * x * x * y * y * y * y * z * z * z * z + 945.0 * x * x * y * y * z * z * z * z * z * z + 66.4453125 * y * y * y * y * y * y * y * y * z * z - 354.375 * y * y * y * y * y * y * z * z * z * z + 472.5 * y * y * y * y * z * z * z * z * z * z) + e_2 * (184.5703125 * x * x * x * x * x * x * x * x - 324.84375 * x * x * x * x * x * x * y * y + 708.75 * x * x * x * x * x * x * z * z + 753.046875 * x * x * x * x * y * y * y * y - 4252.5 * x * x * x * x * y * y * z * z + 1417.5 * x * x * x * x * z * z * z * z + 1328.90625 * x * x * y * y * y * y * y * y + 6378.75 * x * x * y * y * y * y * z * z + 2835.0 * x * x * y * y * z * z * z * z + 1890.0 * x * x * z * z * z * z * z * z + 66.4453125 * y * y * y * y * y * y * y * y + 1417.5 * y * y * y * y * z * z * z * z + 1890.0 * y * y * z * z * z * z * z * z) + e_3 * (1850.625 * x * x * x * x * x * x - 1535.625 * x * x * x * x * y * y + 4961.25 * x * x * x * x * z * z + 10276.875 * x * x * y * y * y * y + 9922.5 * x * x * y * y * z * z + 13230.0 * x * x * z * z * z * z + 1063.125 * y * y * y * y * y * y + 4961.25 * y * y * y * y * z * z + 13230.0 * y * y * z * z * z * z + 1260.0 * z * z * z * z * z * z) + e_4 * (7796.25 * x * x * x * x + 15592.5 * x * x * y * y + 29767.5 * x * x * z * z + 7796.25 * y * y * y * y + 29767.5 * y * y * z * z + 11340.0 * z * z * z * z) + e_5 * (18427.5 * x * x + 18427.5 * y * y + 25515.0 * z * z) + e_6 * (10395.0);
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

        pc_82[k] = e_0 * (std::sqrt(1.816864013671875) * x * x * x * x * x * x * x * x * x * x * x * z - std::sqrt(89.02633666992188) * x * x * x * x * x * x * x * x * x * y * y * z - std::sqrt(291.5057373046875) * x * x * x * x * x * x * x * x * x * z * z * z + std::sqrt(7.2674560546875) * x * x * x * x * x * x * x * y * y * y * y * z + std::sqrt(18656.3671875) * x * x * x * x * x * x * x * y * y * z * z * z + std::sqrt(1291.9921875) * x * x * x * x * x * x * x * z * z * z * z * z + std::sqrt(1228.2000732421875) * x * x * x * x * x * y * y * y * y * y * y * z - std::sqrt(29150.57373046875) * x * x * x * x * x * y * y * y * y * z * z * z - std::sqrt(104651.3671875) * x * x * x * x * x * y * y * z * z * z * z * z + std::sqrt(307.0500183105469) * x * x * x * y * y * y * y * y * y * y * y * z - std::sqrt(74625.46875) * x * x * x * y * y * y * y * y * y * z * z * z + std::sqrt(466409.1796875) * x * x * x * y * y * y * y * z * z * z * z * z - std::sqrt(16.351776123046875) * x * y * y * y * y * y * y * y * y * y * y * z + std::sqrt(2623.5516357421875) * x * y * y * y * y * y * y * y * y * z * z * z - std::sqrt(11627.9296875) * x * y * y * y * y * y * y * z * z * z * z * z) + e_1 * (std::sqrt(181.6864013671875) * x * x * x * x * x * x * x * x * x * z + std::sqrt(7441.875) * x * x * x * x * x * x * x * y * y * z - std::sqrt(14069.794921875) * x * x * x * x * x * x * x * z * z * z + std::sqrt(4912.80029296875) * x * x * x * x * x * y * y * y * y * z - std::sqrt(1935003.779296875) * x * x * x * x * x * y * y * z * z * z + std::sqrt(186046.875) * x * x * x * x * x * z * z * z * z * z + std::sqrt(465.1171875) * x * x * x * y * y * y * y * y * y * z + std::sqrt(1049420.654296875) * x * x * x * y * y * y * y * z * z * z + std::sqrt(744187.5) * x * x * x * y * y * z * z * z * z * z + std::sqrt(588.6639404296875) * x * y * y * y * y * y * y * y * y * z - std::sqrt(215000.419921875) * x * y * y * y * y * y * y * z * z * z + std::sqrt(186046.875) * x * y * y * y * y * z * z * z * z * z) + e_2 * (std::sqrt(11627.9296875) * x * x * x * x * x * x * x * z - std::sqrt(941862.3046875) * x * x * x * x * x * y * y * z + std::sqrt(46511.71875) * x * x * x * x * x * z * z * z + std::sqrt(4197682.6171875) * x * x * x * y * y * y * y * z + std::sqrt(186046.875) * x * x * x * y * y * z * z * z + std::sqrt(6697687.5) * x * x * x * z * z * z * z * z - std::sqrt(104651.3671875) * x * y * y * y * y * y * y * z + std::sqrt(46511.71875) * x * y * y * y * y * z * z * z + std::sqrt(6697687.5) * x * y * y * z * z * z * z * z) + e_3 * (std::sqrt(418605.46875) * x * x * x * x * x * z + std::sqrt(1674421.875) * x * x * x * y * y * z + std::sqrt(90046687.5) * x * x * x * z * z * z + std::sqrt(418605.46875) * x * y * y * y * y * z + std::sqrt(90046687.5) * x * y * y * z * z * z + std::sqrt(11907000.0) * x * z * z * z * z * z) + e_4 * (std::sqrt(82046671.875) * x * x * x * z + std::sqrt(82046671.875) * x * y * y * z + std::sqrt(328186687.5) * x * z * z * z) + e_5 * (std::sqrt(328186687.5) * x * z);
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

        pc_83[k] = e_0 * (-std::sqrt(39.97100830078125) * x * x * x * x * x * x * x * x * x * x * z * z + std::sqrt(5755.8251953125) * x * x * x * x * x * x * x * x * y * y * z * z + std::sqrt(284.23828125) * x * x * x * x * x * x * x * x * z * z * z * z - std::sqrt(19345.968017578125) * x * x * x * x * x * x * y * y * y * y * z * z - std::sqrt(48036.26953125) * x * x * x * x * x * x * y * y * z * z * z * z - std::sqrt(15988.4033203125) * x * x * x * x * y * y * y * y * y * y * z * z + std::sqrt(348191.89453125) * x * x * x * x * y * y * y * y * z * z * z * z + std::sqrt(8993.476867675781) * x * x * y * y * y * y * y * y * y * y * z * z - std::sqrt(63953.61328125) * x * x * y * y * y * y * y * y * z * z * z * z) + e_1 * (-std::sqrt(39.97100830078125) * x * x * x * x * x * x * x * x * x * x + std::sqrt(5755.8251953125) * x * x * x * x * x * x * x * x * y * y - std::sqrt(11551.621398925781) * x * x * x * x * x * x * x * x * z * z - std::sqrt(19345.968017578125) * x * x * x * x * x * x * y * y * y * y - std::sqrt(639.5361328125) * x * x * x * x * x * x * y * y * z * z + std::sqrt(63953.61328125) * x * x * x * x * x * x * z * z * z * z - std::sqrt(15988.4033203125) * x * x * x * x * y * y * y * y * y * y - std::sqrt(195857.94067382812) * x * x * x * x * y * y * y * y * z * z + std::sqrt(63953.61328125) * x * x * x * x * y * y * z * z * z * z + std::sqrt(8993.476867675781) * x * x * y * y * y * y * y * y * y * y + std::sqrt(143895.6298828125) * x * x * y * y * y * y * y * y * z * z - std::sqrt(63953.61328125) * x * x * y * y * y * y * z * z * z * z + std::sqrt(8993.476867675781) * y * y * y * y * y * y * y * y * z * z - std::sqrt(63953.61328125) * y * y * y * y * y * y * z * z * z * z) + e_2 * (-std::sqrt(24981.88018798828) * x * x * x * x * x * x * x * x + std::sqrt(399710.0830078125) * x * x * x * x * x * x * y * y - std::sqrt(575582.51953125) * x * x * x * x * x * x * z * z - std::sqrt(4896448.516845703) * x * x * x * x * y * y * y * y - std::sqrt(575582.51953125) * x * x * x * x * y * y * z * z + std::sqrt(4093031.25) * x * x * x * x * z * z * z * z + std::sqrt(1295060.6689453125) * x * x * y * y * y * y * y * y + std::sqrt(575582.51953125) * x * x * y * y * y * y * z * z + std::sqrt(8993.476867675781) * y * y * y * y * y * y * y * y + std::sqrt(575582.51953125) * y * y * y * y * y * y * z * z - std::sqrt(4093031.25) * y * y * y * y * z * z * z * z) + e_3 * (-std::sqrt(2302330.078125) * x * x * x * x * x * x - std::sqrt(2302330.078125) * x * x * x * x * y * y + std::sqrt(2302330.078125) * x * x * y * y * y * y + std::sqrt(16372125.0) * x * x * z * z * z * z + std::sqrt(2302330.078125) * y * y * y * y * y * y - std::sqrt(16372125.0) * y * y * z * z * z * z) + e_4 * (-std::sqrt(36837281.25) * x * x * x * x + std::sqrt(36837281.25) * x * x * z * z + std::sqrt(36837281.25) * y * y * y * y - std::sqrt(36837281.25) * y * y * z * z) + e_5 * (-std::sqrt(36837281.25) * x * x + std::sqrt(36837281.25) * y * y);
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

        pc_84[k] = e_0 * (-std::sqrt(3.3309173583984375) * x * x * x * x * x * x * x * x * x * x * x * z + std::sqrt(962.6351165771484) * x * x * x * x * x * x * x * x * x * y * y * z + std::sqrt(23.6865234375) * x * x * x * x * x * x * x * x * x * z * z * z - std::sqrt(5875.738220214844) * x * x * x * x * x * x * x * y * y * y * y * z - std::sqrt(7674.43359375) * x * x * x * x * x * x * x * y * y * z * z * z - std::sqrt(652.8598022460938) * x * x * x * x * x * y * y * y * y * y * y * z + std::sqrt(85271.484375) * x * x * x * x * x * y * y * y * y * z * z * z + std::sqrt(6158.866195678711) * x * x * x * y * y * y * y * y * y * y * y * z - std::sqrt(50120.68359375) * x * x * x * y * y * y * y * y * y * z * z * z - std::sqrt(29.978256225585938) * x * y * y * y * y * y * y * y * y * y * y * z + std::sqrt(213.1787109375) * x * y * y * y * y * y * y * y * y * z * z * z) + e_1 * (-std::sqrt(2997.8256225585938) * x * x * x * x * x * x * x * x * x * z + std::sqrt(155407.2802734375) * x * x * x * x * x * x * x * y * y * z + std::sqrt(7674.43359375) * x * x * x * x * x * x * x * z * z * z - std::sqrt(1247575.1110839844) * x * x * x * x * x * y * y * y * y * z - std::sqrt(7674.43359375) * x * x * x * x * x * y * y * z * z * z + std::sqrt(1613549.6630859375) * x * x * x * y * y * y * y * y * y * z - std::sqrt(191860.83984375) * x * x * x * y * y * y * y * z * z * z + std::sqrt(1079.2172241210938) * x * y * y * y * y * y * y * y * y * z - std::sqrt(69069.90234375) * x * y * y * y * y * y * y * z * z * z) + e_2 * (-std::sqrt(431686.8896484375) * x * x * x * x * x * x * x * z + std::sqrt(431686.8896484375) * x * x * x * x * x * y * y * z + std::sqrt(767443.359375) * x * x * x * x * x * z * z * z + std::sqrt(10792172.241210938) * x * x * x * y * y * y * y * z - std::sqrt(3069773.4375) * x * x * x * y * y * z * z * z + std::sqrt(3885182.0068359375) * x * y * y * y * y * y * y * z - std::sqrt(6906990.234375) * x * y * y * y * y * z * z * z) + e_3 * (-std::sqrt(12279093.75) * x * x * x * x * x * z + std::sqrt(49116375.0) * x * x * x * y * y * z + std::sqrt(5457375.0) * x * x * x * z * z * z + std::sqrt(110511843.75) * x * y * y * y * y * z - std::sqrt(49116375.0) * x * y * y * z * z * z) + e_4 * (-std::sqrt(27627960.9375) * x * x * x * z + std::sqrt(248651648.4375) * x * y * y * z);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ab_x, ab_y, ab_z : simd::cache_line_size())
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
        const auto e_6 = pe_6[k];

        pc_85[k] = e_0 * (0.24609375 * x * x * x * x * x * x * x * x * x * x * x * x - 2.4609375 * x * x * x * x * x * x * x * x * x * x * y * y - 4.921875 * x * x * x * x * x * x * x * x * x * x * z * z + 3.69140625 * x * x * x * x * x * x * x * x * y * y * y * y + 54.140625 * x * x * x * x * x * x * x * x * y * y * z * z + 24.609375 * x * x * x * x * x * x * x * x * z * z * z * z + 12.796875 * x * x * x * x * x * x * y * y * y * y * y * y - 127.96875 * x * x * x * x * x * x * y * y * y * y * z * z - 295.3125 * x * x * x * x * x * x * y * y * z * z * z * z + 3.69140625 * x * x * x * x * y * y * y * y * y * y * y * y - 127.96875 * x * x * x * x * y * y * y * y * y * y * z * z + 935.15625 * x * x * x * x * y * y * y * y * z * z * z * z - 2.4609375 * x * x * y * y * y * y * y * y * y * y * y * y + 54.140625 * x * x * y * y * y * y * y * y * y * y * z * z - 295.3125 * x * x * y * y * y * y * y * y * z * z * z * z + 0.24609375 * y * y * y * y * y * y * y * y * y * y * y * y - 4.921875 * y * y * y * y * y * y * y * y * y * y * z * z + 24.609375 * y * y * y * y * y * y * y * y * z * z * z * z) + e_1 * (8.859375 * x * x * x * x * x * x * x * x * x * x - 34.453125 * x * x * x * x * x * x * x * x * y * y - 19.6875 * x * x * x * x * x * x * x * x * z * z + 167.34375 * x * x * x * x * x * x * y * y * y * y - 1023.75 * x * x * x * x * x * x * y * y * z * z + 393.75 * x * x * x * x * x * x * z * z * z * z + 167.34375 * x * x * x * x * y * y * y * y * y * y + 1771.875 * x * x * x * x * y * y * y * y * z * z + 1181.25 * x * x * x * x * y * y * z * z * z * z - 34.453125 * x * x * y * y * y * y * y * y * y * y - 1023.75 * x * x * y * y * y * y * y * y * z * z + 1181.25 * x * x * y * y * y * y * z * z * z * z + 8.859375 * y * y * y * y * y * y * y * y * y * y - 19.6875 * y * y * y * y * y * y * y * y * z * z + 393.75 * y * y * y * y * y * y * z * z * z * z) + e_2 * (172.265625 * x * x * x * x * x * x * x * x - 492.1875 * x * x * x * x * x * x * y * y + 393.75 * x * x * x * x * x * x * z * z + 3396.09375 * x * x * x * x * y * y * y * y + 1181.25 * x * x * x * x * y * y * z * z + 3543.75 * x * x * x * x * z * z * z * z - 492.1875 * x * x * y * y * y * y * y * y + 1181.25 * x * x * y * y * y * y * z * z + 7087.5 * x * x * y * y * z * z * z * z + 172.265625 * y * y * y * y * y * y * y * y + 393.75 * y * y * y * y * y * y * z * z + 3543.75 * y * y * y * y * z * z * z * z) + e_3 * (1575.0 * x * x * x * x * x * x + 4725.0 * x * x * x * x * y * y + 9450.0 * x * x * x * x * z * z + 4725.0 * x * x * y * y * y * y + 18900.0 * x * x * y * y * z * z + 9450.0 * x * x * z * z * z * z + 1575.0 * y * y * y * y * y * y + 9450.0 * y * y * y * y * z * z + 9450.0 * y * y * z * z * z * z) + e_4 * (9450.0 * x * x * x * x + 18900.0 * x * x * y * y + 33075.0 * x * x * z * z + 9450.0 * y * y * y * y + 33075.0 * y * y * z * z + 4725.0 * z * z * z * z) + e_5 * (21735.0 * x * x + 21735.0 * y * y + 18900.0 * z * z) + e_6 * (10395.0);
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

        pc_86[k] = e_0 * (-std::sqrt(1.332366943359375) * x * x * x * x * x * x * x * x * x * x * x * z + std::sqrt(299.7825622558594) * x * x * x * x * x * x * x * x * x * y * y * z + std::sqrt(133.2366943359375) * x * x * x * x * x * x * x * x * x * z * z * z - std::sqrt(3330.9173583984375) * x * x * x * x * x * x * x * y * y * y * y * z - std::sqrt(34108.59375) * x * x * x * x * x * x * x * y * y * z * z * z - std::sqrt(900.6800537109375) * x * x * x * x * x * y * y * y * y * y * y * z + std::sqrt(580379.0405273438) * x * x * x * x * x * y * y * y * y * z * z * z + std::sqrt(1632.1495056152344) * x * x * x * y * y * y * y * y * y * y * y * z - std::sqrt(213178.7109375) * x * x * x * y * y * y * y * y * y * z * z * z - std::sqrt(33.309173583984375) * x * y * y * y * y * y * y * y * y * y * y * z + std::sqrt(3330.9173583984375) * x * y * y * y * y * y * y * y * y * z * z * z) + e_1 * (-std::sqrt(133.2366943359375) * x * x * x * x * x * x * x * x * x * z - std::sqrt(76744.3359375) * x * x * x * x * x * x * x * y * y * z + std::sqrt(53294.677734375) * x * x * x * x * x * x * x * z * z * z + std::sqrt(388518.20068359375) * x * x * x * x * x * y * y * y * y * z + std::sqrt(479652.099609375) * x * x * x * x * x * y * y * z * z * z - std::sqrt(306977.34375) * x * x * x * y * y * y * y * y * y * z + std::sqrt(479652.099609375) * x * x * x * y * y * y * y * z * z * z + std::sqrt(1199.1302490234375) * x * y * y * y * y * y * y * y * y * z + std::sqrt(53294.677734375) * x * y * y * y * y * y * y * z * z * z) + e_2 * (std::sqrt(7674433.59375) * x * x * x * x * x * z * z * z + std::sqrt(30697734.375) * x * x * x * y * y * z * z * z + std::sqrt(7674433.59375) * x * y * y * y * y * z * z * z) + e_3 * (std::sqrt(7674433.59375) * x * x * x * x * x * z + std::sqrt(30697734.375) * x * x * x * y * y * z + std::sqrt(122790937.5) * x * x * x * z * z * z + std::sqrt(7674433.59375) * x * y * y * y * y * z + std::sqrt(122790937.5) * x * y * y * z * z * z) + e_4 * (std::sqrt(276279609.375) * x * x * x * z + std::sqrt(276279609.375) * x * y * y * z + std::sqrt(122790937.5) * x * z * z * z) + e_5 * (std::sqrt(397842637.5) * x * z);
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

        pc_87[k] = e_0 * (-std::sqrt(0.11103057861328125) * x * x * x * x * x * x * x * x * x * x * x * x + std::sqrt(44.4122314453125) * x * x * x * x * x * x * x * x * x * x * y * y + std::sqrt(11.103057861328125) * x * x * x * x * x * x * x * x * x * x * z * z - std::sqrt(802.195930480957) * x * x * x * x * x * x * x * x * y * y * y * y - std::sqrt(4896.448516845703) * x * x * x * x * x * x * x * x * y * y * z * z + std::sqrt(124753.95812988281) * x * x * x * x * x * x * y * y * y * y * z * z + std::sqrt(802.195930480957) * x * x * x * x * y * y * y * y * y * y * y * y - std::sqrt(124753.95812988281) * x * x * x * x * y * y * y * y * y * y * z * z - std::sqrt(44.4122314453125) * x * x * y * y * y * y * y * y * y * y * y * y + std::sqrt(4896.448516845703) * x * x * y * y * y * y * y * y * y * y * z * z + std::sqrt(0.11103057861328125) * y * y * y * y * y * y * y * y * y * y * y * y - std::sqrt(11.103057861328125) * y * y * y * y * y * y * y * y * y * y * z * z) + e_1 * (-std::sqrt(143.8956298828125) * x * x * x * x * x * x * x * x * x * x + std::sqrt(3597.3907470703125) * x * x * x * x * x * x * x * x * y * y + std::sqrt(6395.361328125) * x * x * x * x * x * x * x * x * z * z - std::sqrt(193459.68017578125) * x * x * x * x * x * x * y * y * y * y + std::sqrt(25581.4453125) * x * x * x * x * x * x * y * y * z * z + std::sqrt(193459.68017578125) * x * x * x * x * y * y * y * y * y * y - std::sqrt(3597.3907470703125) * x * x * y * y * y * y * y * y * y * y - std::sqrt(25581.4453125) * x * x * y * y * y * y * y * y * z * z + std::sqrt(143.8956298828125) * y * y * y * y * y * y * y * y * y * y - std::sqrt(6395.361328125) * y * y * y * y * y * y * y * y * z * z) + e_2 * (-std::sqrt(39971.00830078125) * x * x * x * x * x * x * x * x - std::sqrt(159884.033203125) * x * x * x * x * x * x * y * y + std::sqrt(1438956.298828125) * x * x * x * x * x * x * z * z + std::sqrt(1438956.298828125) * x * x * x * x * y * y * z * z + std::sqrt(159884.033203125) * x * x * y * y * y * y * y * y - std::sqrt(1438956.298828125) * x * x * y * y * y * y * z * z + std::sqrt(39971.00830078125) * y * y * y * y * y * y * y * y - std::sqrt(1438956.298828125) * y * y * y * y * y * y * z * z) + e_3 * (-std::sqrt(2558144.53125) * x * x * x * x * x * x - std::sqrt(2558144.53125) * x * x * x * x * y * y + std::sqrt(40930312.5) * x * x * x * x * z * z + std::sqrt(2558144.53125) * x * x * y * y * y * y + std::sqrt(2558144.53125) * y * y * y * y * y * y - std::sqrt(40930312.5) * y * y * y * y * z * z) + e_4 * (-std::sqrt(23023300.78125) * x * x * x * x + std::sqrt(92093203.125) * x * x * z * z + std::sqrt(23023300.78125) * y * y * y * y - std::sqrt(92093203.125) * y * y * z * z) + e_5 * (-std::sqrt(14734912.5) * x * x + std::sqrt(14734912.5) * y * y);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ab_x, ab_y, ab_z : simd::cache_line_size())
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
        const auto e_6 = pe_6[k];

        pc_88[k] = e_0 * (5.4140625 * x * x * x * x * x * x * x * x * x * x * z * z - 108.28125 * x * x * x * x * x * x * x * x * y * y * z * z + 595.546875 * x * x * x * x * x * x * y * y * y * y * z * z - 541.40625 * x * x * x * x * y * y * y * y * y * y * z * z + 135.3515625 * x * x * y * y * y * y * y * y * y * y * z * z) + e_1 * (5.4140625 * x * x * x * x * x * x * x * x * x * x - 108.28125 * x * x * x * x * x * x * x * x * y * y + 135.3515625 * x * x * x * x * x * x * x * x * z * z + 595.546875 * x * x * x * x * x * x * y * y * y * y + 541.40625 * x * x * x * x * x * x * y * y * z * z - 541.40625 * x * x * x * x * y * y * y * y * y * y + 812.109375 * x * x * x * x * y * y * y * y * z * z + 135.3515625 * x * x * y * y * y * y * y * y * y * y + 541.40625 * x * x * y * y * y * y * y * y * z * z + 135.3515625 * y * y * y * y * y * y * y * y * z * z) + e_2 * (135.3515625 * x * x * x * x * x * x * x * x + 541.40625 * x * x * x * x * x * x * y * y + 2165.625 * x * x * x * x * x * x * z * z + 812.109375 * x * x * x * x * y * y * y * y + 6496.875 * x * x * x * x * y * y * z * z + 541.40625 * x * x * y * y * y * y * y * y + 6496.875 * x * x * y * y * y * y * z * z + 135.3515625 * y * y * y * y * y * y * y * y + 2165.625 * y * y * y * y * y * y * z * z) + e_3 * (2165.625 * x * x * x * x * x * x + 6496.875 * x * x * x * x * y * y + 12993.75 * x * x * x * x * z * z + 6496.875 * x * x * y * y * y * y + 25987.5 * x * x * y * y * z * z + 2165.625 * y * y * y * y * y * y + 12993.75 * y * y * y * y * z * z) + e_4 * (12993.75 * x * x * x * x + 25987.5 * x * x * y * y + 25987.5 * x * x * z * z + 12993.75 * y * y * y * y + 25987.5 * y * y * z * z) + e_5 * (25987.5 * x * x + 25987.5 * y * y + 10395.0 * z * z) + e_6 * (10395.0);
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

        pc_89[k] = e_0 * (std::sqrt(2.4426727294921875) * x * x * x * x * x * x * x * x * x * x * x * z - std::sqrt(1526.6704559326172) * x * x * x * x * x * x * x * x * x * y * y * z + std::sqrt(70593.24188232422) * x * x * x * x * x * x * x * y * y * y * y * z - std::sqrt(124761.95233154297) * x * x * x * x * x * y * y * y * y * y * y * z + std::sqrt(17648.310470581055) * x * x * x * y * y * y * y * y * y * y * y * z - std::sqrt(61.06681823730469) * x * y * y * y * y * y * y * y * y * y * y * z) + e_1 * (std::sqrt(2198.4054565429688) * x * x * x * x * x * x * x * x * x * z + std::sqrt(35174.4873046875) * x * x * x * x * x * x * x * y * y * z + std::sqrt(79142.59643554688) * x * x * x * x * x * y * y * y * y * z + std::sqrt(35174.4873046875) * x * x * x * y * y * y * y * y * y * z + std::sqrt(2198.4054565429688) * x * y * y * y * y * y * y * y * y * z) + e_2 * (std::sqrt(879362.1826171875) * x * x * x * x * x * x * x * z + std::sqrt(7914259.6435546875) * x * x * x * x * x * y * y * z + std::sqrt(7914259.6435546875) * x * x * x * y * y * y * y * z + std::sqrt(879362.1826171875) * x * y * y * y * y * y * y * z) + e_3 * (std::sqrt(56279179.6875) * x * x * x * x * x * z + std::sqrt(225116718.75) * x * x * x * y * y * z + std::sqrt(56279179.6875) * x * y * y * y * y * z) + e_4 * (std::sqrt(506512617.1875) * x * x * x * z + std::sqrt(506512617.1875) * x * y * y * z) + e_5 * (std::sqrt(324168075.0) * x * z);
    }

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ab_x, ab_y : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto x = ab_x[k];
        const auto y = ab_y[k];

        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];

        pc_90[k] = e_0 * (0.451171875 * x * x * x * x * x * x * x * x * x * x * x * x - 13.53515625 * x * x * x * x * x * x * x * x * x * x * y * y + 115.048828125 * x * x * x * x * x * x * x * x * y * y * y * y - 203.9296875 * x * x * x * x * x * x * y * y * y * y * y * y + 115.048828125 * x * x * x * x * y * y * y * y * y * y * y * y - 13.53515625 * x * x * y * y * y * y * y * y * y * y * y * y + 0.451171875 * y * y * y * y * y * y * y * y * y * y * y * y) + e_1 * (16.2421875 * x * x * x * x * x * x * x * x * x * x + 81.2109375 * x * x * x * x * x * x * x * x * y * y + 162.421875 * x * x * x * x * x * x * y * y * y * y + 162.421875 * x * x * x * x * y * y * y * y * y * y + 81.2109375 * x * x * y * y * y * y * y * y * y * y + 16.2421875 * y * y * y * y * y * y * y * y * y * y) + e_2 * (406.0546875 * x * x * x * x * x * x * x * x + 1624.21875 * x * x * x * x * x * x * y * y + 2436.328125 * x * x * x * x * y * y * y * y + 1624.21875 * x * x * y * y * y * y * y * y + 406.0546875 * y * y * y * y * y * y * y * y) + e_3 * (4331.25 * x * x * x * x * x * x + 12993.75 * x * x * x * x * y * y + 12993.75 * x * x * y * y * y * y + 4331.25 * y * y * y * y * y * y) + e_4 * (19490.625 * x * x * x * x + 38981.25 * x * x * y * y + 19490.625 * y * y * y * y) + e_5 * (31185.0 * x * x + 31185.0 * y * y) + e_6 * (10395.0);
    }

    // NOTE: the values of a combination of angular components are stored as one
    // row of nvalues columns, with the component on bra side running slowest. The
    // rows which the symmetry relates to an already formed one are copied from it,
    // and the atom pairs beyond the reach of every pair of primitives are set to
    // zero.

    const size_t sources[169] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 1, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 2, 15, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 3, 16, 29, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 4, 17, 30, 43, 56, 57, 58, 59, 60, 61, 62, 63, 64, 5, 18, 31, 44, 57, 70, 71, 72, 73, 74, 75, 76, 77, 6, 19, 32, 45, 58, 71, 84, 85, 86, 87, 88, 89, 90, 7, 20, 33, 46, 59, 72, 85, 98, 99, 100, 101, 102, 103, 8, 21, 34, 47, 60, 73, 86, 99, 112, 113, 114, 115, 116, 9, 22, 35, 48, 61, 74, 87, 100, 113, 126, 127, 128, 129, 10, 23, 36, 49, 62, 75, 88, 101, 114, 127, 140, 141, 142, 11, 24, 37, 50, 63, 76, 89, 102, 115, 128, 141, 154, 155, 12, 25, 38, 51, 64, 77, 90, 103, 116, 129, 142, 155, 168};

    for (size_t m = 0; m < 169; m++)
    {
        auto *pv = values + m * nvalues;

        const auto *pc = values + sources[m] * nvalues;

        if (pv != pc) std::copy(pc, pc + nmax, pv);

        std::fill(pv + nmax, pv + nvalues, 0.0);
    }
}

}  // namespace simdovl
