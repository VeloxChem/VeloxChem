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



#include "SimdOverlapRecGH.hpp"

#include <algorithm>
#include <cmath>
#include <string>

#include "ErrorHandler.hpp"
#include "MathConst.hpp"
#include "ScreeningFunc.hpp"
#include "SimdAlign.hpp"
#include "SimdDimensions.hpp"

namespace simdovl {  // simdovl namespace

auto
compute_gh_overlap(double                         *values,
                   const size_t                    nvalues,
                   const CBasisFunction           &bra,
                   const CBasisFunction           &ket,
                   const std::vector<CSimdMatrix> &harmonics,
                   const CSimdMatrix              &coordinates,
                   const double                    threshold) -> void
{
    if ((bra.get_angular_momentum() != 4) || (ket.get_angular_momentum() != 5))
    {
        errors::assertMsgCritical(
            false, std::string("SimdOverlapRecGH.compute_gh_overlap: Basis functions must be of angular momenta four and five"));
    }

    if (harmonics.size() < 9)
    {
        errors::assertMsgCritical(
            false, std::string("SimdOverlapRecGH.compute_gh_overlap: Harmonics must reach angular momentum nine"));
    }

    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("SimdOverlapRecGH.compute_gh_overlap: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nprims = nprim_a * nprim_b;

    // NOTE: the pairs of primitives are screened with the threshold of the
    // integrals divided by their number, as their contributions accumulate into
    // a single value and the error of the sum is bounded by the number of terms.

    const auto dimensions = simdfunc::make_column_dimensions(
        bra, ket, nvalues, coordinates, screenfunc::two_center_overlap_primitive_bound, threshold / static_cast<double>(nprims));

    const auto nmax = dimensions.back();

    if (nmax == 0)
    {
        std::fill(values, values + 99 * nvalues, 0.0);

        return;
    }

    // NOTE: the buffer holds the contracted prefactors of the terms alone, as the
    // integrals of the angular components are formed straight into the values and
    // are not written a second time.

    auto buffer = CSimdMatrix(5, nmax);

    auto *pe_0 = buffer.data(0);
    auto *pe_1 = buffer.data(1);
    auto *pe_2 = buffer.data(2);
    auto *pe_3 = buffer.data(3);
    auto *pe_4 = buffer.data(4);

    std::fill(pe_0, pe_0 + nmax, 0.0);
    std::fill(pe_1, pe_1 + nmax, 0.0);
    std::fill(pe_2, pe_2 + nmax, 0.0);
    std::fill(pe_3, pe_3 + nmax, 0.0);
    std::fill(pe_4, pe_4 + nmax, 0.0);

    const auto *ab_2 = coordinates.data(6);

    constexpr auto fpi = mathconst::pi_value();

    // accumulate the prefactor of each term over the pairs of primitives

    for (size_t i = 0; i < nprim_a; i++)
    {
        const auto aexp = a_exps[i];

        const auto anorm = a_norms[i];

        for (size_t j = 0; j < nprim_b; j++)
        {
            const auto ncols = dimensions[i * nprim_b + j];

            if (ncols == 0) continue;

            const auto bexp = b_exps[j];

            const auto fexp = aexp + bexp;

            const auto fmu = aexp * bexp / fexp;

            const auto fovl = fpi / fexp;

            const auto fbase = anorm * b_norms[j] * fovl * std::sqrt(fovl);

            const auto f_0 = fbase * aexp * fmu / fexp / fexp / fexp / fexp / fexp;

            const auto f_1 = fbase * aexp * fmu * fmu / fexp / fexp / fexp / fexp / fexp;

            const auto f_2 = fbase * aexp * fmu * fmu * fmu / fexp / fexp / fexp / fexp / fexp;

            const auto f_3 = fbase * aexp * fmu * fmu * fmu * fmu / fexp / fexp / fexp / fexp / fexp;

            const auto f_4 = fbase * aexp / fexp / fexp / fexp / fexp / fexp;

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
        }
    }

    // NOTE: the geometry of a term is a solid harmonic of the vector between the
    // atoms times a power of their squared distance.

    const auto *ph1_m1 = harmonics[0].data(0);
    const auto *ph1_0 = harmonics[0].data(1);
    const auto *ph1_p1 = harmonics[0].data(2);
    const auto *ph3_m3 = harmonics[2].data(0);
    const auto *ph3_m2 = harmonics[2].data(1);
    const auto *ph3_m1 = harmonics[2].data(2);
    const auto *ph3_0 = harmonics[2].data(3);
    const auto *ph3_p1 = harmonics[2].data(4);
    const auto *ph3_p2 = harmonics[2].data(5);
    const auto *ph3_p3 = harmonics[2].data(6);
    const auto *ph5_m5 = harmonics[4].data(0);
    const auto *ph5_m4 = harmonics[4].data(1);
    const auto *ph5_m3 = harmonics[4].data(2);
    const auto *ph5_m2 = harmonics[4].data(3);
    const auto *ph5_m1 = harmonics[4].data(4);
    const auto *ph5_0 = harmonics[4].data(5);
    const auto *ph5_p1 = harmonics[4].data(6);
    const auto *ph5_p2 = harmonics[4].data(7);
    const auto *ph5_p3 = harmonics[4].data(8);
    const auto *ph5_p4 = harmonics[4].data(9);
    const auto *ph5_p5 = harmonics[4].data(10);
    const auto *ph7_m7 = harmonics[6].data(0);
    const auto *ph7_m6 = harmonics[6].data(1);
    const auto *ph7_m5 = harmonics[6].data(2);
    const auto *ph7_m4 = harmonics[6].data(3);
    const auto *ph7_m3 = harmonics[6].data(4);
    const auto *ph7_m2 = harmonics[6].data(5);
    const auto *ph7_m1 = harmonics[6].data(6);
    const auto *ph7_0 = harmonics[6].data(7);
    const auto *ph7_p1 = harmonics[6].data(8);
    const auto *ph7_p2 = harmonics[6].data(9);
    const auto *ph7_p3 = harmonics[6].data(10);
    const auto *ph7_p4 = harmonics[6].data(11);
    const auto *ph7_p5 = harmonics[6].data(12);
    const auto *ph7_p6 = harmonics[6].data(13);
    const auto *ph7_p7 = harmonics[6].data(14);
    const auto *ph9_m9 = harmonics[8].data(0);
    const auto *ph9_m8 = harmonics[8].data(1);
    const auto *ph9_m7 = harmonics[8].data(2);
    const auto *ph9_m6 = harmonics[8].data(3);
    const auto *ph9_m5 = harmonics[8].data(4);
    const auto *ph9_m4 = harmonics[8].data(5);
    const auto *ph9_m3 = harmonics[8].data(6);
    const auto *ph9_m2 = harmonics[8].data(7);
    const auto *ph9_m1 = harmonics[8].data(8);
    const auto *ph9_0 = harmonics[8].data(9);
    const auto *ph9_p1 = harmonics[8].data(10);
    const auto *ph9_p2 = harmonics[8].data(11);
    const auto *ph9_p3 = harmonics[8].data(12);
    const auto *ph9_p4 = harmonics[8].data(13);
    const auto *ph9_p5 = harmonics[8].data(14);
    const auto *ph9_p6 = harmonics[8].data(15);
    const auto *ph9_p7 = harmonics[8].data(16);
    const auto *ph9_p8 = harmonics[8].data(17);
    const auto *ph9_p9 = harmonics[8].data(18);

    // NOTE: the rows of the values are not aligned, as they start at the offset
    // of this combination of basis functions in the values block, so they are kept
    // out of the aligned clauses below.

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

    // NOTE: the factors of the terms depend on the angular momenta alone, so they
    // are formed once for the whole matrix instead of once for every atom pair.

    const auto fs_6615_16 = std::sqrt(413.4375);
    const auto fs_19845_8 = std::sqrt(2480.625);
    const auto fs_75_8 = std::sqrt(9.375);
    const auto fs_735_4 = std::sqrt(183.75);
    const auto fs_3645_8 = std::sqrt(455.625);
    const auto fs_875_40898 = std::sqrt(875.0 / 40898.0);
    const auto fs_150_169 = std::sqrt(150.0 / 169.0);
    const auto fs_735_121 = std::sqrt(735.0 / 121.0);
    const auto fs_10 = std::sqrt(10.0);
    const auto fs_49_11819522 = std::sqrt(49.0 / 11819522.0);
    const auto fs_441_2431 = std::sqrt(441.0 / 2431.0);
    const auto fs_1750_5909761 = std::sqrt(1750.0 / 5909761.0);
    const auto fs_2_507 = std::sqrt(2.0 / 507.0);
    const auto fs_980_61347 = std::sqrt(980.0 / 61347.0);
    const auto fs_5_242 = std::sqrt(5.0 / 242.0);
    const auto fs_496125_512 = std::sqrt(968.994140625);
    const auto f_63_2 = 31.5;
    const auto f_15_2 = 7.5;
    const auto f_21 = 21.0;
    const auto f_27_2 = 13.5;
    const auto f_70_143 = 70.0 / 143.0;
    const auto f_30_13 = 30.0 / 13.0;
    const auto f_42_11 = 42.0 / 11.0;
    const auto f_2 = 2.0;
    const auto f_21_2431 = 21.0 / 2431.0;
    const auto fs_245_2431 = std::sqrt(245.0 / 2431.0);
    const auto f_140_2431 = 140.0 / 2431.0;
    const auto f_2_13 = 2.0 / 13.0;
    const auto f_28_143 = 28.0 / 143.0;
    const auto f_1_11 = 1.0 / 11.0;
    const auto f_315_16 = 19.6875;
    const auto fs_1323_4 = std::sqrt(330.75);
    const auto fs_441_8 = std::sqrt(55.125);
    const auto fs_375_8 = std::sqrt(46.875);
    const auto fs_147 = std::sqrt(147.0);
    const auto fs_81_8 = std::sqrt(10.125);
    const auto fs_68600_184041 = std::sqrt(68600.0 / 184041.0);
    const auto fs_1225_858 = std::sqrt(1225.0 / 858.0);
    const auto fs_750_169 = std::sqrt(750.0 / 169.0);
    const auto fs_588_121 = std::sqrt(588.0 / 121.0);
    const auto fs_2_9 = std::sqrt(2.0 / 9.0);
    const auto fs_2205_11819522 = std::sqrt(2205.0 / 11819522.0);
    const auto fs_2205_41327 = std::sqrt(2205.0 / 41327.0);
    const auto fs_274400_53187849 = std::sqrt(274400.0 / 53187849.0);
    const auto fs_2450_123981 = std::sqrt(2450.0 / 123981.0);
    const auto fs_10_507 = std::sqrt(10.0 / 507.0);
    const auto fs_784_61347 = std::sqrt(784.0 / 61347.0);
    const auto fs_1_2178 = std::sqrt(1.0 / 2178.0);
    const auto fs_11025_512 = std::sqrt(21.533203125);
    const auto fs_2205_16 = std::sqrt(137.8125);
    const auto fs_875_16 = std::sqrt(54.6875);
    const auto fs_245_4 = std::sqrt(61.25);
    const auto fs_17150_20449 = std::sqrt(17150.0 / 20449.0);
    const auto fs_350_143 = std::sqrt(350.0 / 143.0);
    const auto fs_875_169 = std::sqrt(875.0 / 169.0);
    const auto fs_245_121 = std::sqrt(245.0 / 121.0);
    const auto fs_735_1074502 = std::sqrt(735.0 / 1074502.0);
    const auto fs_2205_82654 = std::sqrt(2205.0 / 82654.0);
    const auto fs_68600_5909761 = std::sqrt(68600.0 / 5909761.0);
    const auto fs_1400_41327 = std::sqrt(1400.0 / 41327.0);
    const auto fs_35_1521 = std::sqrt(35.0 / 1521.0);
    const auto fs_980_184041 = std::sqrt(980.0 / 184041.0);
    const auto fs_945_32 = std::sqrt(29.53125);
    const auto fs_105_8 = std::sqrt(13.125);
    const auto f_175_143 = 175.0 / 143.0;
    const auto fs_4900_1859 = std::sqrt(4900.0 / 1859.0);
    const auto fs_105_242 = std::sqrt(105.0 / 242.0);
    const auto fs_2205_1074502 = std::sqrt(2205.0 / 1074502.0);
    const auto fs_1029_82654 = std::sqrt(1029.0 / 82654.0);
    const auto f_350_2431 = 350.0 / 2431.0;
    const auto fs_19600_537251 = std::sqrt(19600.0 / 537251.0);
    const auto fs_70_61347 = std::sqrt(70.0 / 61347.0);
    const auto fs_24500_5577 = std::sqrt(24500.0 / 5577.0);
    const auto fs_441_41327 = std::sqrt(441.0 / 41327.0);
    const auto fs_98000_1611753 = std::sqrt(98000.0 / 1611753.0);
    const auto fs_33075_64 = std::sqrt(516.796875);
    const auto fs_525_16 = std::sqrt(32.8125);
    const auto fs_3675_16 = std::sqrt(229.6875);
    const auto fs_23625_163592 = std::sqrt(23625.0 / 163592.0);
    const auto fs_525_169 = std::sqrt(525.0 / 169.0);
    const auto fs_3675_484 = std::sqrt(3675.0 / 484.0);
    const auto fs_49_1074502 = std::sqrt(49.0 / 1074502.0);
    const auto fs_196_2431 = std::sqrt(196.0 / 2431.0);
    const auto fs_23625_11819522 = std::sqrt(23625.0 / 11819522.0);
    const auto fs_7_507 = std::sqrt(7.0 / 507.0);
    const auto fs_1225_61347 = std::sqrt(1225.0 / 61347.0);
    const auto fs_1323_64 = std::sqrt(20.671875);
    const auto fs_3969_2 = std::sqrt(1984.5);
    const auto fs_147_16 = std::sqrt(9.1875);
    const auto fs_729_2 = std::sqrt(364.5);
    const auto fs_92575_163592 = std::sqrt(92575.0 / 163592.0);
    const auto fs_3675_1144 = std::sqrt(3675.0 / 1144.0);
    const auto fs_147_484 = std::sqrt(147.0 / 484.0);
    const auto fs_8 = std::sqrt(8.0);
    const auto fs_1960_5909761 = std::sqrt(1960.0 / 5909761.0);
    const auto fs_3920_41327 = std::sqrt(3920.0 / 41327.0);
    const auto fs_92575_11819522 = std::sqrt(92575.0 / 11819522.0);
    const auto fs_3675_82654 = std::sqrt(3675.0 / 82654.0);
    const auto fs_49_61347 = std::sqrt(49.0 / 61347.0);
    const auto fs_2_121 = std::sqrt(2.0 / 121.0);
    const auto fs_99225_128 = std::sqrt(775.1953125);
    const auto f_63_4 = 15.75;
    const auto f_42 = 42.0;
    const auto f_21_2 = 10.5;
    const auto f_18 = 18.0;
    const auto f_665_429 = 665.0 / 429.0;
    const auto fs_4375_3432 = std::sqrt(4375.0 / 3432.0);
    const auto f_21_11 = 21.0 / 11.0;
    const auto f_8_3 = 8.0 / 3.0;
    const auto f_126_2431 = 126.0 / 2431.0;
    const auto fs_6615_82654 = std::sqrt(6615.0 / 82654.0);
    const auto f_1330_7293 = 1330.0 / 7293.0;
    const auto fs_4375_247962 = std::sqrt(4375.0 / 247962.0);
    const auto f_14_143 = 14.0 / 143.0;
    const auto f_4_33 = 4.0 / 33.0;
    const auto f_105_4 = 26.25;
    const auto f_147_8 = 18.375;
    const auto fs_1323_8 = std::sqrt(165.375);
    const auto fs_125_32 = std::sqrt(3.90625);
    const auto f_49_4 = 12.25;
    const auto fs_243_8 = std::sqrt(30.375);
    const auto fs_214375_122694 = std::sqrt(214375.0 / 122694.0);
    const auto fs_175_3718 = std::sqrt(175.0 / 3718.0);
    const auto fs_125_338 = std::sqrt(125.0 / 338.0);
    const auto f_49_22 = 49.0 / 22.0;
    const auto fs_2_3 = std::sqrt(2.0 / 3.0);
    const auto fs_23520_5909761 = std::sqrt(23520.0 / 5909761.0);
    const auto fs_2352_41327 = std::sqrt(2352.0 / 41327.0);
    const auto fs_428750_17729283 = std::sqrt(428750.0 / 17729283.0);
    const auto fs_350_537251 = std::sqrt(350.0 / 537251.0);
    const auto fs_5_3042 = std::sqrt(5.0 / 3042.0);
    const auto f_49_429 = 49.0 / 429.0;
    const auto fs_1_726 = std::sqrt(1.0 / 726.0);
    const auto fs_33075_512 = std::sqrt(64.599609375);
    const auto fs_38115_128 = std::sqrt(297.7734375);
    const auto fs_4235_32 = std::sqrt(132.34375);
    const auto f_35_26 = 35.0 / 26.0;
    const auto fs_1225_3718 = std::sqrt(1225.0 / 3718.0);
    const auto fs_35_8 = std::sqrt(4.375);
    const auto fs_5145_537251 = std::sqrt(5145.0 / 537251.0);
    const auto fs_1470_41327 = std::sqrt(1470.0 / 41327.0);
    const auto f_35_221 = 35.0 / 221.0;
    const auto fs_2450_537251 = std::sqrt(2450.0 / 537251.0);
    const auto fs_35_3042 = std::sqrt(35.0 / 3042.0);
    const auto fs_14175_64 = std::sqrt(221.484375);
    const auto fs_1575_16 = std::sqrt(98.4375);
    const auto fs_300125_122694 = std::sqrt(300125.0 / 122694.0);
    const auto fs_1575_484 = std::sqrt(1575.0 / 484.0);
    const auto fs_21168_537251 = std::sqrt(21168.0 / 537251.0);
    const auto fs_600250_17729283 = std::sqrt(600250.0 / 17729283.0);
    const auto fs_175_20449 = std::sqrt(175.0 / 20449.0);
    const auto fs_84375_163592 = std::sqrt(84375.0 / 163592.0);
    const auto fs_2625_1144 = std::sqrt(2625.0 / 1144.0);
    const auto fs_147_537251 = std::sqrt(147.0 / 537251.0);
    const auto fs_1372_41327 = std::sqrt(1372.0 / 41327.0);
    const auto fs_84375_11819522 = std::sqrt(84375.0 / 11819522.0);
    const auto fs_2625_82654 = std::sqrt(2625.0 / 82654.0);
    const auto fs_945_4 = std::sqrt(236.25);
    const auto fs_375_16 = std::sqrt(23.4375);
    const auto fs_105 = std::sqrt(105.0);
    const auto fs_54675_40898 = std::sqrt(54675.0 / 40898.0);
    const auto fs_75_286 = std::sqrt(75.0 / 286.0);
    const auto fs_375_169 = std::sqrt(375.0 / 169.0);
    const auto fs_420_121 = std::sqrt(420.0 / 121.0);
    const auto fs_1715_1074502 = std::sqrt(1715.0 / 1074502.0);
    const auto fs_5145_82654 = std::sqrt(5145.0 / 82654.0);
    const auto fs_109350_5909761 = std::sqrt(109350.0 / 5909761.0);
    const auto fs_150_41327 = std::sqrt(150.0 / 41327.0);
    const auto fs_5_507 = std::sqrt(5.0 / 507.0);
    const auto fs_560_61347 = std::sqrt(560.0 / 61347.0);
    const auto fs_15309_64 = std::sqrt(239.203125);
    const auto fs_3087_2 = std::sqrt(1543.5);
    const auto fs_1701_16 = std::sqrt(106.3125);
    const auto fs_567_2 = std::sqrt(283.5);
    const auto fs_15625_8712 = std::sqrt(15625.0 / 8712.0);
    const auto fs_70225_44616 = std::sqrt(70225.0 / 44616.0);
    const auto fs_1701_484 = std::sqrt(1701.0 / 484.0);
    const auto fs_56_9 = std::sqrt(56.0 / 9.0);
    const auto fs_30870_5909761 = std::sqrt(30870.0 / 5909761.0);
    const auto fs_3087_41327 = std::sqrt(3087.0 / 41327.0);
    const auto fs_15625_629442 = std::sqrt(15625.0 / 629442.0);
    const auto fs_70225_3223506 = std::sqrt(70225.0 / 3223506.0);
    const auto fs_189_20449 = std::sqrt(189.0 / 20449.0);
    const auto fs_14_1089 = std::sqrt(14.0 / 1089.0);
    const auto fs_77175_128 = std::sqrt(602.9296875);
    const auto fs_189_16 = std::sqrt(11.8125);
    const auto fs_9261_4 = std::sqrt(2315.25);
    const auto fs_21_4 = std::sqrt(5.25);
    const auto fs_1701_4 = std::sqrt(425.25);
    const auto fs_179200_61347 = std::sqrt(179200.0 / 61347.0);
    const auto fs_2500_1859 = std::sqrt(2500.0 / 1859.0);
    const auto fs_21_121 = std::sqrt(21.0 / 121.0);
    const auto fs_28_3 = std::sqrt(28.0 / 3.0);
    const auto fs_148176_5909761 = std::sqrt(148176.0 / 5909761.0);
    const auto fs_2940_41327 = std::sqrt(2940.0 / 41327.0);
    const auto fs_716800_17729283 = std::sqrt(716800.0 / 17729283.0);
    const auto fs_10000_537251 = std::sqrt(10000.0 / 537251.0);
    const auto fs_28_61347 = std::sqrt(28.0 / 61347.0);
    const auto fs_7_363 = std::sqrt(7.0 / 363.0);
    const auto fs_231525_256 = std::sqrt(904.39453125);
    const auto fs_16641_128 = std::sqrt(130.0078125);
    const auto fs_30375_128 = std::sqrt(237.3046875);
    const auto fs_125_4 = std::sqrt(31.25);
    const auto fs_1849_32 = std::sqrt(57.78125);
    const auto fs_3375_32 = std::sqrt(105.46875);
    const auto fs_243_4 = std::sqrt(60.75);
    const auto fs_147175_245388 = std::sqrt(147175.0 / 245388.0);
    const auto fs_175_484 = std::sqrt(175.0 / 484.0);
    const auto fs_500_169 = std::sqrt(500.0 / 169.0);
    const auto fs_1849_968 = std::sqrt(1849.0 / 968.0);
    const auto fs_3375_968 = std::sqrt(3375.0 / 968.0);
    const auto fs_4_3 = std::sqrt(4.0 / 3.0);
    const auto fs_144060_5909761 = std::sqrt(144060.0 / 5909761.0);
    const auto fs_30870_537251 = std::sqrt(30870.0 / 537251.0);
    const auto fs_147175_17729283 = std::sqrt(147175.0 / 17729283.0);
    const auto fs_175_34969 = std::sqrt(175.0 / 34969.0);
    const auto fs_20_1521 = std::sqrt(20.0 / 1521.0);
    const auto fs_1849_368082 = std::sqrt(1849.0 / 368082.0);
    const auto fs_375_40898 = std::sqrt(375.0 / 40898.0);
    const auto fs_1_363 = std::sqrt(1.0 / 363.0);
    const auto fs_33075_256 = std::sqrt(129.19921875);
    const auto fs_675 = std::sqrt(675.0);
    const auto fs_300 = std::sqrt(300.0);
    const auto fs_1750_61347 = std::sqrt(1750.0 / 61347.0);
    const auto fs_1200_121 = std::sqrt(1200.0 / 121.0);
    const auto fs_43218_537251 = std::sqrt(43218.0 / 537251.0);
    const auto fs_7000_17729283 = std::sqrt(7000.0 / 17729283.0);
    const auto fs_1600_61347 = std::sqrt(1600.0 / 61347.0);
    const auto fs_9375_7436 = std::sqrt(9375.0 / 7436.0);
    const auto fs_3375_1144 = std::sqrt(3375.0 / 1144.0);
    const auto fs_49_41327 = std::sqrt(49.0 / 41327.0);
    const auto fs_9375_537251 = std::sqrt(9375.0 / 537251.0);
    const auto fs_3375_82654 = std::sqrt(3375.0 / 82654.0);
    const auto fs_25515_64 = std::sqrt(398.671875);
    const auto fs_2835_16 = std::sqrt(177.1875);
    const auto fs_1875_968 = std::sqrt(1875.0 / 968.0);
    const auto fs_6075_14872 = std::sqrt(6075.0 / 14872.0);
    const auto fs_2835_484 = std::sqrt(2835.0 / 484.0);
    const auto fs_2940_537251 = std::sqrt(2940.0 / 537251.0);
    const auto fs_1875_69938 = std::sqrt(1875.0 / 69938.0);
    const auto fs_6075_1074502 = std::sqrt(6075.0 / 1074502.0);
    const auto fs_315_20449 = std::sqrt(315.0 / 20449.0);
    const auto fs_945_64 = std::sqrt(14.765625);
    const auto fs_105_16 = std::sqrt(6.5625);
    const auto fs_664225_490776 = std::sqrt(664225.0 / 490776.0);
    const auto fs_25_132 = std::sqrt(25.0 / 132.0);
    const auto fs_105_484 = std::sqrt(105.0 / 484.0);
    const auto fs_15435_1074502 = std::sqrt(15435.0 / 1074502.0);
    const auto fs_664225_35458566 = std::sqrt(664225.0 / 35458566.0);
    const auto fs_25_9537 = std::sqrt(25.0 / 9537.0);
    const auto fs_35_61347 = std::sqrt(35.0 / 61347.0);
    const auto fs_252 = std::sqrt(252.0);
    const auto fs_23625_64 = std::sqrt(369.140625);
    const auto fs_9261_8 = std::sqrt(1157.625);
    const auto fs_875_32 = std::sqrt(27.34375);
    const auto fs_112 = std::sqrt(112.0);
    const auto fs_2625_16 = std::sqrt(164.0625);
    const auto fs_1701_8 = std::sqrt(212.625);
    const auto fs_42025_122694 = std::sqrt(42025.0 / 122694.0);
    const auto fs_46225_40898 = std::sqrt(46225.0 / 40898.0);
    const auto fs_875_338 = std::sqrt(875.0 / 338.0);
    const auto fs_448_121 = std::sqrt(448.0 / 121.0);
    const auto fs_2625_484 = std::sqrt(2625.0 / 484.0);
    const auto fs_14_3 = std::sqrt(14.0 / 3.0);
    const auto fs_164640_5909761 = std::sqrt(164640.0 / 5909761.0);
    const auto fs_35280_537251 = std::sqrt(35280.0 / 537251.0);
    const auto fs_84050_17729283 = std::sqrt(84050.0 / 17729283.0);
    const auto fs_92450_5909761 = std::sqrt(92450.0 / 5909761.0);
    const auto fs_1792_184041 = std::sqrt(1792.0 / 184041.0);
    const auto fs_875_61347 = std::sqrt(875.0 / 61347.0);
    const auto fs_7_726 = std::sqrt(7.0 / 726.0);
    const auto fs_231525_512 = std::sqrt(452.197265625);
    const auto fs_9747_32 = std::sqrt(304.59375);
    const auto fs_28125_128 = std::sqrt(219.7265625);
    const auto fs_2646 = std::sqrt(2646.0);
    const auto fs_1083_8 = std::sqrt(135.375);
    const auto fs_3125_32 = std::sqrt(97.65625);
    const auto fs_486 = std::sqrt(486.0);
    const auto fs_2450_61347 = std::sqrt(2450.0 / 61347.0);
    const auto fs_109375_81796 = std::sqrt(109375.0 / 81796.0);
    const auto fs_1083_242 = std::sqrt(1083.0 / 242.0);
    const auto fs_3125_968 = std::sqrt(3125.0 / 968.0);
    const auto fs_32_3 = std::sqrt(32.0 / 3.0);
    const auto fs_518616_5909761 = std::sqrt(518616.0 / 5909761.0);
    const auto fs_36015_537251 = std::sqrt(36015.0 / 537251.0);
    const auto fs_9800_17729283 = std::sqrt(9800.0 / 17729283.0);
    const auto fs_109375_5909761 = std::sqrt(109375.0 / 5909761.0);
    const auto fs_722_61347 = std::sqrt(722.0 / 61347.0);
    const auto fs_3125_368082 = std::sqrt(3125.0 / 368082.0);
    const auto fs_8_363 = std::sqrt(8.0 / 363.0);
    const auto fs_33075_32 = std::sqrt(1033.59375);
    const auto fs_135_64 = std::sqrt(2.109375);
    const auto fs_2205_2 = std::sqrt(1102.5);
    const auto fs_15_16 = std::sqrt(0.9375);
    const auto fs_405_2 = std::sqrt(202.5);
    const auto fs_462875_368082 = std::sqrt(462875.0 / 368082.0);
    const auto fs_15_484 = std::sqrt(15.0 / 484.0);
    const auto fs_40_9 = std::sqrt(40.0 / 9.0);
    const auto fs_691488_5909761 = std::sqrt(691488.0 / 5909761.0);
    const auto fs_925750_53187849 = std::sqrt(925750.0 / 53187849.0);
    const auto fs_5_61347 = std::sqrt(5.0 / 61347.0);
    const auto fs_10_1089 = std::sqrt(10.0 / 1089.0);
    const auto fs_55125_128 = std::sqrt(430.6640625);
    const auto fs_16875_3718 = std::sqrt(16875.0 / 3718.0);
    const auto fs_343_41327 = std::sqrt(343.0 / 41327.0);
    const auto fs_33750_537251 = std::sqrt(33750.0 / 537251.0);
    const auto fs_6000_1859 = std::sqrt(6000.0 / 1859.0);
    const auto fs_1225_41327 = std::sqrt(1225.0 / 41327.0);
    const auto fs_24000_537251 = std::sqrt(24000.0 / 537251.0);
    const auto fs_14175_16 = std::sqrt(885.9375);
    const auto fs_1575_4 = std::sqrt(393.75);
    const auto fs_66125_122694 = std::sqrt(66125.0 / 122694.0);
    const auto fs_1575_121 = std::sqrt(1575.0 / 121.0);
    const auto fs_33075_537251 = std::sqrt(33075.0 / 537251.0);
    const auto fs_132250_17729283 = std::sqrt(132250.0 / 17729283.0);
    const auto fs_700_20449 = std::sqrt(700.0 / 20449.0);
    const auto f_5_4 = 1.25;
    const auto fs_175_4 = std::sqrt(43.75);
    const auto fs_4000_20449 = std::sqrt(4000.0 / 20449.0);
    const auto f_5_13 = 5.0 / 13.0;
    const auto fs_175_121 = std::sqrt(175.0 / 121.0);
    const auto fs_51450_537251 = std::sqrt(51450.0 / 537251.0);
    const auto fs_16000_5909761 = std::sqrt(16000.0 / 5909761.0);
    const auto f_1_39 = 1.0 / 39.0;
    const auto fs_700_184041 = std::sqrt(700.0 / 184041.0);
    const auto fs_5445_32 = std::sqrt(170.15625);
    const auto fs_6615_4 = std::sqrt(1653.75);
    const auto f_5 = 5.0;
    const auto fs_605_8 = std::sqrt(75.625);
    const auto fs_1215_4 = std::sqrt(303.75);
    const auto fs_875_507 = std::sqrt(875.0 / 507.0);
    const auto f_20_13 = 20.0 / 13.0;
    const auto fs_5_2 = std::sqrt(2.5);
    const auto fs_20_3 = std::sqrt(20.0 / 3.0);
    const auto fs_720300_5909761 = std::sqrt(720300.0 / 5909761.0);
    const auto fs_3500_146523 = std::sqrt(3500.0 / 146523.0);
    const auto f_4_39 = 4.0 / 39.0;
    const auto fs_10_1521 = std::sqrt(10.0 / 1521.0);
    const auto fs_5_363 = std::sqrt(5.0 / 363.0);
    const auto fs_165375_256 = std::sqrt(645.99609375);
    const auto f_45_2 = 22.5;
    const auto f_105_2 = 52.5;
    const auto f_15 = 15.0;
    const auto f_700_429 = 700.0 / 429.0;
    const auto f_30_11 = 30.0 / 11.0;
    const auto f_10_3 = 10.0 / 3.0;
    const auto f_882_2431 = 882.0 / 2431.0;
    const auto f_1400_7293 = 1400.0 / 7293.0;
    const auto f_20_143 = 20.0 / 143.0;
    const auto f_5_33 = 5.0 / 33.0;
    const auto f_525_16 = 32.8125;

    // NOTE: the rows are formed in 59 loops, as the vectorizer runs out of
    // registers with all 99 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph1_0, ph1_p1, ph3_0, ph3_p1, ph5_0, ph5_p1, ph7_0, ph7_p1, ph9_0, ph9_p1, ph9_p8, ph9_p9, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h1_0 = ph1_0[k];
        const auto h1_p1 = ph1_p1[k];
        const auto h3_0 = ph3_0[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h5_0 = ph5_0[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h7_0 = ph7_0[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h9_0 = ph9_0[k];
        const auto h9_p1 = ph9_p1[k];
        const auto h9_p8 = ph9_p8[k];
        const auto h9_p9 = ph9_p9[k];

        pc_0[k] = e_0 * (fs_6615_16 * h3_p1 - fs_19845_8 * r_2 * h1_p1) + e_1 * (fs_75_8 * h5_p1 - fs_735_4 * r_2 * h3_p1 + fs_3645_8 * r_4 * h1_p1) + e_2 * (fs_875_40898 * h7_p1 - fs_150_169 * r_2 * h5_p1 + fs_735_121 * r_4 * h3_p1 - fs_10 * r_6 * h1_p1) + e_3 * (fs_49_11819522 * h9_p1 - fs_441_2431 * h9_p9 - fs_1750_5909761 * r_2 * h7_p1 + fs_2_507 * r_4 * h5_p1 - fs_980_61347 * r_6 * h3_p1 + fs_5_242 * r_8 * h1_p1) + fs_496125_512 * e_4 * h1_p1;

        pc_1[k] = e_0 * (f_63_2 * h3_0 - f_63_2 * r_2 * h1_0) + e_1 * (f_15_2 * h5_0 - f_21 * r_2 * h3_0 + f_27_2 * r_4 * h1_0) + e_2 * (f_70_143 * h7_0 - f_30_13 * r_2 * h5_0 + f_42_11 * r_4 * h3_0 - f_2 * r_6 * h1_0) + e_3 * (f_21_2431 * h9_0 - fs_245_2431 * h9_p8 - f_140_2431 * r_2 * h7_0 + f_2_13 * r_4 * h5_0 - f_28_143 * r_6 * h3_0 + f_1_11 * r_8 * h1_0) + f_315_16 * e_4 * h1_0;
    }

    // NOTE: the rows are formed in 59 loops, as the vectorizer runs out of
    // registers with all 99 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph1_p1, ph3_p1, ph5_p1, ph7_p1, ph7_p7, ph9_p1, ph9_p7, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h1_p1 = ph1_p1[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h7_p7 = ph7_p7[k];
        const auto h9_p1 = ph9_p1[k];
        const auto h9_p7 = ph9_p7[k];

        pc_2[k] = e_0 * (-fs_1323_4 * h3_p1 + fs_441_8 * r_2 * h1_p1) + e_1 * (-fs_375_8 * h5_p1 + fs_147 * r_2 * h3_p1 - fs_81_8 * r_4 * h1_p1) + e_2 * (-fs_68600_184041 * h7_p1 - fs_1225_858 * h7_p7 + fs_750_169 * r_2 * h5_p1 - fs_588_121 * r_4 * h3_p1 + fs_2_9 * r_6 * h1_p1) + e_3 * (-fs_2205_11819522 * h9_p1 - fs_2205_41327 * h9_p7 + fs_274400_53187849 * r_2 * h7_p1 + fs_2450_123981 * r_2 * h7_p7 - fs_10_507 * r_4 * h5_p1 + fs_784_61347 * r_6 * h3_p1 - fs_1_2178 * r_8 * h1_p1) - fs_11025_512 * e_4 * h1_p1;
    }

    // NOTE: the rows are formed in 59 loops, as the vectorizer runs out of
    // registers with all 99 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph3_p2, ph3_p3, ph5_p2, ph5_p3, ph5_p5, ph7_p2, ph7_p3, ph7_p5, ph7_p6, ph9_p2, ph9_p3, ph9_p5, ph9_p6, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h3_p2 = ph3_p2[k];
        const auto h3_p3 = ph3_p3[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h5_p5 = ph5_p5[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h7_p3 = ph7_p3[k];
        const auto h7_p5 = ph7_p5[k];
        const auto h7_p6 = ph7_p6[k];
        const auto h9_p2 = ph9_p2[k];
        const auto h9_p3 = ph9_p3[k];
        const auto h9_p5 = ph9_p5[k];
        const auto h9_p6 = ph9_p6[k];

        pc_3[k] = fs_2205_16 * e_0 * h3_p2 + e_1 * (fs_875_16 * h5_p2 - fs_245_4 * r_2 * h3_p2) + e_2 * (fs_17150_20449 * h7_p2 - fs_350_143 * h7_p6 - fs_875_169 * r_2 * h5_p2 + fs_245_121 * r_4 * h3_p2) + e_3 * (fs_735_1074502 * h9_p2 - fs_2205_82654 * h9_p6 - fs_68600_5909761 * r_2 * h7_p2 + fs_1400_41327 * r_2 * h7_p6 + fs_35_1521 * r_4 * h5_p2 - fs_980_184041 * r_6 * h3_p2);

        pc_4[k] = -fs_945_32 * e_0 * h3_p3 + e_1 * (-fs_375_8 * h5_p3 - fs_75_8 * h5_p5 + fs_105_8 * r_2 * h3_p3) + e_2 * (-f_175_143 * h7_p3 - fs_4900_1859 * h7_p5 + fs_750_169 * r_2 * h5_p3 + fs_150_169 * r_2 * h5_p5 - fs_105_242 * r_4 * h3_p3) + e_3 * (-fs_2205_1074502 * h9_p3 - fs_1029_82654 * h9_p5 + f_350_2431 * r_2 * h7_p3 + fs_19600_537251 * r_2 * h7_p5 - fs_10_507 * r_4 * h5_p3 - fs_2_507 * r_4 * h5_p5 + fs_70_61347 * r_6 * h3_p3);
    }

    // NOTE: the rows are formed in 59 loops, as the vectorizer runs out of
    // registers with all 99 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph3_m3, ph5_m5, ph5_m4, ph5_m3, ph7_m5, ph7_m4, ph7_m3, ph9_m5, ph9_m4, ph9_m3, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h3_m3 = ph3_m3[k];
        const auto h5_m5 = ph5_m5[k];
        const auto h5_m4 = ph5_m4[k];
        const auto h5_m3 = ph5_m3[k];
        const auto h7_m5 = ph7_m5[k];
        const auto h7_m4 = ph7_m4[k];
        const auto h7_m3 = ph7_m3[k];
        const auto h9_m5 = ph9_m5[k];
        const auto h9_m4 = ph9_m4[k];
        const auto h9_m3 = ph9_m3[k];

        pc_5[k] = f_15_2 * e_1 * h5_m4 + e_2 * (fs_24500_5577 * h7_m4 - f_30_13 * r_2 * h5_m4) + e_3 * (fs_441_41327 * h9_m4 - fs_98000_1611753 * r_2 * h7_m4 + f_2_13 * r_4 * h5_m4);

        pc_6[k] = -fs_945_32 * e_0 * h3_m3 + e_1 * (fs_75_8 * h5_m5 - fs_375_8 * h5_m3 + fs_105_8 * r_2 * h3_m3) + e_2 * (fs_4900_1859 * h7_m5 - f_175_143 * h7_m3 - fs_150_169 * r_2 * h5_m5 + fs_750_169 * r_2 * h5_m3 - fs_105_242 * r_4 * h3_m3) + e_3 * (fs_1029_82654 * h9_m5 - fs_2205_1074502 * h9_m3 - fs_19600_537251 * r_2 * h7_m5 + f_350_2431 * r_2 * h7_m3 + fs_2_507 * r_4 * h5_m5 - fs_10_507 * r_4 * h5_m3 + fs_70_61347 * r_6 * h3_m3);
    }

    // NOTE: the rows are formed in 59 loops, as the vectorizer runs out of
    // registers with all 99 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph3_m2, ph5_m2, ph7_m6, ph7_m2, ph9_m6, ph9_m2, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h3_m2 = ph3_m2[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h7_m6 = ph7_m6[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h9_m6 = ph9_m6[k];
        const auto h9_m2 = ph9_m2[k];

        pc_7[k] = fs_2205_16 * e_0 * h3_m2 + e_1 * (fs_875_16 * h5_m2 - fs_245_4 * r_2 * h3_m2) + e_2 * (fs_350_143 * h7_m6 + fs_17150_20449 * h7_m2 - fs_875_169 * r_2 * h5_m2 + fs_245_121 * r_4 * h3_m2) + e_3 * (fs_2205_82654 * h9_m6 + fs_735_1074502 * h9_m2 - fs_1400_41327 * r_2 * h7_m6 - fs_68600_5909761 * r_2 * h7_m2 + fs_35_1521 * r_4 * h5_m2 - fs_980_184041 * r_6 * h3_m2);
    }

    // NOTE: the rows are formed in 59 loops, as the vectorizer runs out of
    // registers with all 99 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph1_m1, ph3_m1, ph5_m1, ph7_m7, ph7_m1, ph9_m9, ph9_m8, ph9_m7, ph9_m1, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h1_m1 = ph1_m1[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h7_m7 = ph7_m7[k];
        const auto h7_m1 = ph7_m1[k];
        const auto h9_m9 = ph9_m9[k];
        const auto h9_m8 = ph9_m8[k];
        const auto h9_m7 = ph9_m7[k];
        const auto h9_m1 = ph9_m1[k];

        pc_8[k] = e_0 * (-fs_1323_4 * h3_m1 + fs_441_8 * r_2 * h1_m1) + e_1 * (-fs_375_8 * h5_m1 + fs_147 * r_2 * h3_m1 - fs_81_8 * r_4 * h1_m1) + e_2 * (fs_1225_858 * h7_m7 - fs_68600_184041 * h7_m1 + fs_750_169 * r_2 * h5_m1 - fs_588_121 * r_4 * h3_m1 + fs_2_9 * r_6 * h1_m1) + e_3 * (fs_2205_41327 * h9_m7 - fs_2205_11819522 * h9_m1 - fs_2450_123981 * r_2 * h7_m7 + fs_274400_53187849 * r_2 * h7_m1 - fs_10_507 * r_4 * h5_m1 + fs_784_61347 * r_6 * h3_m1 - fs_1_2178 * r_8 * h1_m1) - fs_11025_512 * e_4 * h1_m1;

        pc_9[k] = fs_245_2431 * e_3 * h9_m8;

        pc_10[k] = e_0 * (-fs_6615_16 * h3_m1 + fs_19845_8 * r_2 * h1_m1) + e_1 * (-fs_75_8 * h5_m1 + fs_735_4 * r_2 * h3_m1 - fs_3645_8 * r_4 * h1_m1) + e_2 * (-fs_875_40898 * h7_m1 + fs_150_169 * r_2 * h5_m1 - fs_735_121 * r_4 * h3_m1 + fs_10 * r_6 * h1_m1) + e_3 * (fs_441_2431 * h9_m9 - fs_49_11819522 * h9_m1 + fs_1750_5909761 * r_2 * h7_m1 - fs_2_507 * r_4 * h5_m1 + fs_980_61347 * r_6 * h3_m1 - fs_5_242 * r_8 * h1_m1) - fs_496125_512 * e_4 * h1_m1;
    }

    // NOTE: the rows are formed in 59 loops, as the vectorizer runs out of
    // registers with all 99 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph1_p1, ph3_p1, ph3_p2, ph5_p1, ph5_p2, ph7_p1, ph7_p2, ph7_p7, ph9_p1, ph9_p2, ph9_p7, ph9_p8, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h1_p1 = ph1_p1[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h3_p2 = ph3_p2[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h7_p7 = ph7_p7[k];
        const auto h9_p1 = ph9_p1[k];
        const auto h9_p2 = ph9_p2[k];
        const auto h9_p7 = ph9_p7[k];
        const auto h9_p8 = ph9_p8[k];

        pc_11[k] = -fs_33075_64 * e_0 * h3_p2 + e_1 * (-fs_525_16 * h5_p2 + fs_3675_16 * r_2 * h3_p2) + e_2 * (-fs_23625_163592 * h7_p2 + fs_525_169 * r_2 * h5_p2 - fs_3675_484 * r_4 * h3_p2) + e_3 * (-fs_49_1074502 * h9_p2 - fs_196_2431 * h9_p8 + fs_23625_11819522 * r_2 * h7_p2 - fs_7_507 * r_4 * h5_p2 + fs_1225_61347 * r_6 * h3_p2);

        pc_12[k] = e_0 * (-fs_1323_64 * h3_p1 - fs_3969_2 * r_2 * h1_p1) + e_1 * (-fs_375_8 * h5_p1 + fs_147_16 * r_2 * h3_p1 + fs_729_2 * r_4 * h1_p1) + e_2 * (-fs_92575_163592 * h7_p1 + fs_3675_1144 * h7_p7 + fs_750_169 * r_2 * h5_p1 - fs_147_484 * r_4 * h3_p1 - fs_8 * r_6 * h1_p1) + e_3 * (-fs_1960_5909761 * h9_p1 - fs_3920_41327 * h9_p7 + fs_92575_11819522 * r_2 * h7_p1 - fs_3675_82654 * r_2 * h7_p7 - fs_10_507 * r_4 * h5_p1 + fs_49_61347 * r_6 * h3_p1 + fs_2_121 * r_8 * h1_p1) + fs_99225_128 * e_4 * h1_p1;
    }

    // NOTE: the rows are formed in 59 loops, as the vectorizer runs out of
    // registers with all 99 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph1_0, ph3_0, ph5_0, ph7_0, ph7_p6, ph9_0, ph9_p6, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h1_0 = ph1_0[k];
        const auto h3_0 = ph3_0[k];
        const auto h5_0 = ph5_0[k];
        const auto h7_0 = ph7_0[k];
        const auto h7_p6 = ph7_p6[k];
        const auto h9_0 = ph9_0[k];
        const auto h9_p6 = ph9_p6[k];

        pc_13[k] = e_0 * (f_63_4 * h3_0 - f_42 * r_2 * h1_0) + e_1 * (-f_15_2 * h5_0 - f_21_2 * r_2 * h3_0 + f_18 * r_4 * h1_0) + e_2 * (-f_665_429 * h7_0 + fs_4375_3432 * h7_p6 + f_30_13 * r_2 * h5_0 + f_21_11 * r_4 * h3_0 - f_8_3 * r_6 * h1_0) + e_3 * (-f_126_2431 * h9_0 - fs_6615_82654 * h9_p6 + f_1330_7293 * r_2 * h7_0 - fs_4375_247962 * r_2 * h7_p6 - f_2_13 * r_4 * h5_0 - f_14_143 * r_6 * h3_0 + f_4_33 * r_8 * h1_0) + f_105_4 * e_4 * h1_0;
    }

    // NOTE: the rows are formed in 59 loops, as the vectorizer runs out of
    // registers with all 99 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph1_p1, ph3_p1, ph5_p1, ph5_p5, ph7_p1, ph7_p5, ph9_p1, ph9_p5, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h1_p1 = ph1_p1[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h5_p5 = ph5_p5[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h7_p5 = ph7_p5[k];
        const auto h9_p1 = ph9_p1[k];
        const auto h9_p5 = ph9_p5[k];

        pc_14[k] = e_0 * (-f_147_8 * h3_p1 + fs_1323_8 * r_2 * h1_p1) + e_1 * (fs_125_32 * h5_p1 + fs_525_16 * h5_p5 + f_49_4 * r_2 * h3_p1 - fs_243_8 * r_4 * h1_p1) + e_2 * (fs_214375_122694 * h7_p1 + fs_175_3718 * h7_p5 - fs_125_338 * r_2 * h5_p1 - fs_525_169 * r_2 * h5_p5 - f_49_22 * r_4 * h3_p1 + fs_2_3 * r_6 * h1_p1) + e_3 * (fs_23520_5909761 * h9_p1 - fs_2352_41327 * h9_p5 - fs_428750_17729283 * r_2 * h7_p1 - fs_350_537251 * r_2 * h7_p5 + fs_5_3042 * r_4 * h5_p1 + fs_7_507 * r_4 * h5_p5 + f_49_429 * r_6 * h3_p1 - fs_1_726 * r_8 * h1_p1) - fs_33075_512 * e_4 * h1_p1;
    }

    // NOTE: the rows are formed in 59 loops, as the vectorizer runs out of
    // registers with all 99 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph3_m3, ph3_p2, ph5_m3, ph5_p2, ph5_p4, ph7_m3, ph7_p2, ph7_p4, ph9_m3, ph9_p2, ph9_p4, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h3_m3 = ph3_m3[k];
        const auto h3_p2 = ph3_p2[k];
        const auto h5_m3 = ph5_m3[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h5_p4 = ph5_p4[k];
        const auto h7_m3 = ph7_m3[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h7_p4 = ph7_p4[k];
        const auto h9_m3 = ph9_m3[k];
        const auto h9_p2 = ph9_p2[k];
        const auto h9_p4 = ph9_p4[k];

        pc_15[k] = fs_38115_128 * e_0 * h3_p2 + e_1 * (fs_125_32 * h5_p2 + fs_375_8 * h5_p4 - fs_4235_32 * r_2 * h3_p2) + e_2 * (-f_35_26 * h7_p2 - fs_1225_3718 * h7_p4 - fs_125_338 * r_2 * h5_p2 - fs_750_169 * r_2 * h5_p4 + fs_35_8 * r_4 * h3_p2) + e_3 * (-fs_5145_537251 * h9_p2 - fs_1470_41327 * h9_p4 + f_35_221 * r_2 * h7_p2 + fs_2450_537251 * r_2 * h7_p4 + fs_5_3042 * r_4 * h5_p2 + fs_10_507 * r_4 * h5_p4 - fs_35_3042 * r_6 * h3_p2);

        pc_16[k] = -fs_14175_64 * e_0 * h3_m3 + e_1 * (-f_15_2 * h5_m3 + fs_1575_16 * r_2 * h3_m3) + e_2 * (fs_300125_122694 * h7_m3 + f_30_13 * r_2 * h5_m3 - fs_1575_484 * r_4 * h3_m3) + e_3 * (fs_21168_537251 * h9_m3 - fs_600250_17729283 * r_2 * h7_m3 - f_2_13 * r_4 * h5_m3 + fs_175_20449 * r_6 * h3_m3);
    }

    // NOTE: the rows are formed in 59 loops, as the vectorizer runs out of
    // registers with all 99 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph3_m2, ph5_m4, ph5_m2, ph7_m4, ph7_m2, ph9_m4, ph9_m2, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h3_m2 = ph3_m2[k];
        const auto h5_m4 = ph5_m4[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h7_m4 = ph7_m4[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h9_m4 = ph9_m4[k];
        const auto h9_m2 = ph9_m2[k];

        pc_17[k] = fs_38115_128 * e_0 * h3_m2 + e_1 * (-fs_375_8 * h5_m4 + fs_125_32 * h5_m2 - fs_4235_32 * r_2 * h3_m2) + e_2 * (fs_1225_3718 * h7_m4 - f_35_26 * h7_m2 + fs_750_169 * r_2 * h5_m4 - fs_125_338 * r_2 * h5_m2 + fs_35_8 * r_4 * h3_m2) + e_3 * (fs_1470_41327 * h9_m4 - fs_5145_537251 * h9_m2 - fs_2450_537251 * r_2 * h7_m4 + f_35_221 * r_2 * h7_m2 - fs_10_507 * r_4 * h5_m4 + fs_5_3042 * r_4 * h5_m2 - fs_35_3042 * r_6 * h3_m2);
    }

    // NOTE: the rows are formed in 59 loops, as the vectorizer runs out of
    // registers with all 99 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph1_m1, ph3_m1, ph5_m5, ph5_m1, ph7_m7, ph7_m6, ph7_m5, ph7_m1, ph9_m7, ph9_m6, ph9_m5, ph9_m1, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h1_m1 = ph1_m1[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h5_m5 = ph5_m5[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h7_m7 = ph7_m7[k];
        const auto h7_m6 = ph7_m6[k];
        const auto h7_m5 = ph7_m5[k];
        const auto h7_m1 = ph7_m1[k];
        const auto h9_m7 = ph9_m7[k];
        const auto h9_m6 = ph9_m6[k];
        const auto h9_m5 = ph9_m5[k];
        const auto h9_m1 = ph9_m1[k];

        pc_18[k] = e_0 * (-f_147_8 * h3_m1 + fs_1323_8 * r_2 * h1_m1) + e_1 * (-fs_525_16 * h5_m5 + fs_125_32 * h5_m1 + f_49_4 * r_2 * h3_m1 - fs_243_8 * r_4 * h1_m1) + e_2 * (-fs_175_3718 * h7_m5 + fs_214375_122694 * h7_m1 + fs_525_169 * r_2 * h5_m5 - fs_125_338 * r_2 * h5_m1 - f_49_22 * r_4 * h3_m1 + fs_2_3 * r_6 * h1_m1) + e_3 * (fs_2352_41327 * h9_m5 + fs_23520_5909761 * h9_m1 + fs_350_537251 * r_2 * h7_m5 - fs_428750_17729283 * r_2 * h7_m1 - fs_7_507 * r_4 * h5_m5 + fs_5_3042 * r_4 * h5_m1 + f_49_429 * r_6 * h3_m1 - fs_1_726 * r_8 * h1_m1) - fs_33075_512 * e_4 * h1_m1;

        pc_19[k] = -fs_4375_3432 * e_2 * h7_m6 + e_3 * (fs_6615_82654 * h9_m6 + fs_4375_247962 * r_2 * h7_m6);

        pc_20[k] = e_0 * (fs_1323_64 * h3_m1 + fs_3969_2 * r_2 * h1_m1) + e_1 * (fs_375_8 * h5_m1 - fs_147_16 * r_2 * h3_m1 - fs_729_2 * r_4 * h1_m1) + e_2 * (-fs_3675_1144 * h7_m7 + fs_92575_163592 * h7_m1 - fs_750_169 * r_2 * h5_m1 + fs_147_484 * r_4 * h3_m1 + fs_8 * r_6 * h1_m1) + e_3 * (fs_3920_41327 * h9_m7 + fs_1960_5909761 * h9_m1 + fs_3675_82654 * r_2 * h7_m7 - fs_92575_11819522 * r_2 * h7_m1 + fs_10_507 * r_4 * h5_m1 - fs_49_61347 * r_6 * h3_m1 - fs_2_121 * r_8 * h1_m1) - fs_99225_128 * e_4 * h1_m1;
    }

    // NOTE: the rows are formed in 59 loops, as the vectorizer runs out of
    // registers with all 99 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph3_m2, ph3_p3, ph5_m2, ph5_p3, ph7_m2, ph7_p3, ph7_p7, ph9_m8, ph9_m2, ph9_p3, ph9_p7, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h3_m2 = ph3_m2[k];
        const auto h3_p3 = ph3_p3[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h7_p3 = ph7_p3[k];
        const auto h7_p7 = ph7_p7[k];
        const auto h9_m8 = ph9_m8[k];
        const auto h9_m2 = ph9_m2[k];
        const auto h9_p3 = ph9_p3[k];
        const auto h9_p7 = ph9_p7[k];

        pc_21[k] = fs_33075_64 * e_0 * h3_m2 + e_1 * (fs_525_16 * h5_m2 - fs_3675_16 * r_2 * h3_m2) + e_2 * (fs_23625_163592 * h7_m2 - fs_525_169 * r_2 * h5_m2 + fs_3675_484 * r_4 * h3_m2) + e_3 * (fs_196_2431 * h9_m8 + fs_49_1074502 * h9_m2 - fs_23625_11819522 * r_2 * h7_m2 + fs_7_507 * r_4 * h5_m2 - fs_1225_61347 * r_6 * h3_m2);

        pc_22[k] = fs_14175_64 * e_0 * h3_p3 + e_1 * (f_15_2 * h5_p3 - fs_1575_16 * r_2 * h3_p3) + e_2 * (fs_84375_163592 * h7_p3 - fs_2625_1144 * h7_p7 - f_30_13 * r_2 * h5_p3 + fs_1575_484 * r_4 * h3_p3) + e_3 * (fs_147_537251 * h9_p3 - fs_1372_41327 * h9_p7 - fs_84375_11819522 * r_2 * h7_p3 + fs_2625_82654 * r_2 * h7_p7 + f_2_13 * r_4 * h5_p3 - fs_175_20449 * r_6 * h3_p3);
    }

    // NOTE: the rows are formed in 59 loops, as the vectorizer runs out of
    // registers with all 99 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph3_p2, ph5_p2, ph7_p2, ph7_p6, ph9_p2, ph9_p6, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h3_p2 = ph3_p2[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h7_p6 = ph7_p6[k];
        const auto h9_p2 = ph9_p2[k];
        const auto h9_p6 = ph9_p6[k];

        pc_23[k] = -fs_945_4 * e_0 * h3_p2 + e_1 * (fs_375_16 * h5_p2 + fs_105 * r_2 * h3_p2) + e_2 * (fs_54675_40898 * h7_p2 + fs_75_286 * h7_p6 - fs_375_169 * r_2 * h5_p2 - fs_420_121 * r_4 * h3_p2) + e_3 * (fs_1715_1074502 * h9_p2 - fs_5145_82654 * h9_p6 - fs_109350_5909761 * r_2 * h7_p2 - fs_150_41327 * r_2 * h7_p6 + fs_5_507 * r_4 * h5_p2 + fs_560_61347 * r_6 * h3_p2);
    }

    // NOTE: the rows are formed in 59 loops, as the vectorizer runs out of
    // registers with all 99 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph1_p1, ph3_p1, ph5_p5, ph7_p1, ph7_p5, ph9_p1, ph9_p5, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h1_p1 = ph1_p1[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h5_p5 = ph5_p5[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h7_p5 = ph7_p5[k];
        const auto h9_p1 = ph9_p1[k];
        const auto h9_p5 = ph9_p5[k];

        pc_24[k] = e_0 * (-fs_15309_64 * h3_p1 - fs_3087_2 * r_2 * h1_p1) + e_1 * (-f_15_2 * h5_p5 + fs_1701_16 * r_2 * h3_p1 + fs_567_2 * r_4 * h1_p1) + e_2 * (fs_15625_8712 * h7_p1 + fs_70225_44616 * h7_p5 + f_30_13 * r_2 * h5_p5 - fs_1701_484 * r_4 * h3_p1 - fs_56_9 * r_6 * h1_p1) + e_3 * (fs_30870_5909761 * h9_p1 - fs_3087_41327 * h9_p5 - fs_15625_629442 * r_2 * h7_p1 - fs_70225_3223506 * r_2 * h7_p5 - f_2_13 * r_4 * h5_p5 + fs_189_20449 * r_6 * h3_p1 + fs_14_1089 * r_8 * h1_p1) + fs_77175_128 * e_4 * h1_p1;
    }

    // NOTE: the rows are formed in 59 loops, as the vectorizer runs out of
    // registers with all 99 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph1_0, ph3_0, ph5_0, ph5_p4, ph7_0, ph7_p4, ph9_0, ph9_p4, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h1_0 = ph1_0[k];
        const auto h3_0 = ph3_0[k];
        const auto h5_0 = ph5_0[k];
        const auto h5_p4 = ph5_p4[k];
        const auto h7_0 = ph7_0[k];
        const auto h7_p4 = ph7_p4[k];
        const auto h9_0 = ph9_0[k];
        const auto h9_p4 = ph9_p4[k];

        pc_25[k] = e_0 * (-fs_189_16 * h3_0 - fs_9261_4 * r_2 * h1_0) + e_1 * (-fs_525_16 * h5_0 - fs_375_16 * h5_p4 + fs_21_4 * r_2 * h3_0 + fs_1701_4 * r_4 * h1_0) + e_2 * (fs_179200_61347 * h7_0 + fs_2500_1859 * h7_p4 + fs_525_169 * r_2 * h5_0 + fs_375_169 * r_2 * h5_p4 - fs_21_121 * r_4 * h3_0 - fs_28_3 * r_6 * h1_0) + e_3 * (fs_148176_5909761 * h9_0 - fs_2940_41327 * h9_p4 - fs_716800_17729283 * r_2 * h7_0 - fs_10000_537251 * r_2 * h7_p4 - fs_7_507 * r_4 * h5_0 - fs_5_507 * r_4 * h5_p4 + fs_28_61347 * r_6 * h3_0 + fs_7_363 * r_8 * h1_0) + fs_231525_256 * e_4 * h1_0;
    }

    // NOTE: the rows are formed in 59 loops, as the vectorizer runs out of
    // registers with all 99 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph1_p1, ph3_m2, ph3_p1, ph3_p3, ph5_m2, ph5_p1, ph7_m2, ph7_p1, ph7_p3, ph9_m2, ph9_p1, ph9_p3, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h1_p1 = ph1_p1[k];
        const auto h3_m2 = ph3_m2[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h3_p3 = ph3_p3[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h7_p3 = ph7_p3[k];
        const auto h9_m2 = ph9_m2[k];
        const auto h9_p1 = ph9_p1[k];
        const auto h9_p3 = ph9_p3[k];

        pc_26[k] = e_0 * (-fs_16641_128 * h3_p1 - fs_30375_128 * h3_p3 + fs_1323_4 * r_2 * h1_p1) + e_1 * (fs_125_4 * h5_p1 + fs_1849_32 * r_2 * h3_p1 + fs_3375_32 * r_2 * h3_p3 - fs_243_4 * r_4 * h1_p1) + e_2 * (-fs_147175_245388 * h7_p1 + fs_175_484 * h7_p3 - fs_500_169 * r_2 * h5_p1 - fs_1849_968 * r_4 * h3_p1 - fs_3375_968 * r_4 * h3_p3 + fs_4_3 * r_6 * h1_p1) + e_3 * (-fs_144060_5909761 * h9_p1 - fs_30870_537251 * h9_p3 + fs_147175_17729283 * r_2 * h7_p1 - fs_175_34969 * r_2 * h7_p3 + fs_20_1521 * r_4 * h5_p1 + fs_1849_368082 * r_6 * h3_p1 + fs_375_40898 * r_6 * h3_p3 - fs_1_363 * r_8 * h1_p1) - fs_33075_256 * e_4 * h1_p1;

        pc_27[k] = fs_675 * e_0 * h3_m2 + e_1 * (-fs_525_16 * h5_m2 - fs_300 * r_2 * h3_m2) + e_2 * (fs_1750_61347 * h7_m2 + fs_525_169 * r_2 * h5_m2 + fs_1200_121 * r_4 * h3_m2) + e_3 * (fs_43218_537251 * h9_m2 - fs_7000_17729283 * r_2 * h7_m2 - fs_7_507 * r_4 * h5_m2 - fs_1600_61347 * r_6 * h3_m2);
    }

    // NOTE: the rows are formed in 59 loops, as the vectorizer runs out of
    // registers with all 99 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph1_m1, ph3_m3, ph3_m1, ph5_m4, ph5_m1, ph7_m4, ph7_m3, ph7_m1, ph9_m4, ph9_m3, ph9_m1, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h1_m1 = ph1_m1[k];
        const auto h3_m3 = ph3_m3[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h5_m4 = ph5_m4[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h7_m4 = ph7_m4[k];
        const auto h7_m3 = ph7_m3[k];
        const auto h7_m1 = ph7_m1[k];
        const auto h9_m4 = ph9_m4[k];
        const auto h9_m3 = ph9_m3[k];
        const auto h9_m1 = ph9_m1[k];

        pc_28[k] = e_0 * (fs_30375_128 * h3_m3 - fs_16641_128 * h3_m1 + fs_1323_4 * r_2 * h1_m1) + e_1 * (fs_125_4 * h5_m1 - fs_3375_32 * r_2 * h3_m3 + fs_1849_32 * r_2 * h3_m1 - fs_243_4 * r_4 * h1_m1) + e_2 * (-fs_175_484 * h7_m3 - fs_147175_245388 * h7_m1 - fs_500_169 * r_2 * h5_m1 + fs_3375_968 * r_4 * h3_m3 - fs_1849_968 * r_4 * h3_m1 + fs_4_3 * r_6 * h1_m1) + e_3 * (fs_30870_537251 * h9_m3 - fs_144060_5909761 * h9_m1 + fs_175_34969 * r_2 * h7_m3 + fs_147175_17729283 * r_2 * h7_m1 + fs_20_1521 * r_4 * h5_m1 - fs_375_40898 * r_6 * h3_m3 + fs_1849_368082 * r_6 * h3_m1 - fs_1_363 * r_8 * h1_m1) - fs_33075_256 * e_4 * h1_m1;

        pc_29[k] = fs_375_16 * e_1 * h5_m4 + e_2 * (-fs_2500_1859 * h7_m4 - fs_375_169 * r_2 * h5_m4) + e_3 * (fs_2940_41327 * h9_m4 + fs_10000_537251 * r_2 * h7_m4 + fs_5_507 * r_4 * h5_m4);
    }

    // NOTE: the rows are formed in 59 loops, as the vectorizer runs out of
    // registers with all 99 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph1_m1, ph3_m1, ph5_m5, ph7_m5, ph7_m1, ph9_m5, ph9_m1, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h1_m1 = ph1_m1[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h5_m5 = ph5_m5[k];
        const auto h7_m5 = ph7_m5[k];
        const auto h7_m1 = ph7_m1[k];
        const auto h9_m5 = ph9_m5[k];
        const auto h9_m1 = ph9_m1[k];

        pc_30[k] = e_0 * (fs_15309_64 * h3_m1 + fs_3087_2 * r_2 * h1_m1) + e_1 * (f_15_2 * h5_m5 - fs_1701_16 * r_2 * h3_m1 - fs_567_2 * r_4 * h1_m1) + e_2 * (-fs_70225_44616 * h7_m5 - fs_15625_8712 * h7_m1 - f_30_13 * r_2 * h5_m5 + fs_1701_484 * r_4 * h3_m1 + fs_56_9 * r_6 * h1_m1) + e_3 * (fs_3087_41327 * h9_m5 - fs_30870_5909761 * h9_m1 + fs_70225_3223506 * r_2 * h7_m5 + fs_15625_629442 * r_2 * h7_m1 + f_2_13 * r_4 * h5_m5 - fs_189_20449 * r_6 * h3_m1 - fs_14_1089 * r_8 * h1_m1) - fs_77175_128 * e_4 * h1_m1;
    }

    // NOTE: the rows are formed in 59 loops, as the vectorizer runs out of
    // registers with all 99 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph3_m3, ph3_m2, ph5_m3, ph5_m2, ph7_m7, ph7_m6, ph7_m3, ph7_m2, ph9_m7, ph9_m6, ph9_m3, ph9_m2, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h3_m3 = ph3_m3[k];
        const auto h3_m2 = ph3_m2[k];
        const auto h5_m3 = ph5_m3[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h7_m7 = ph7_m7[k];
        const auto h7_m6 = ph7_m6[k];
        const auto h7_m3 = ph7_m3[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h9_m7 = ph9_m7[k];
        const auto h9_m6 = ph9_m6[k];
        const auto h9_m3 = ph9_m3[k];
        const auto h9_m2 = ph9_m2[k];

        pc_31[k] = fs_945_4 * e_0 * h3_m2 + e_1 * (-fs_375_16 * h5_m2 - fs_105 * r_2 * h3_m2) + e_2 * (-fs_75_286 * h7_m6 - fs_54675_40898 * h7_m2 + fs_375_169 * r_2 * h5_m2 + fs_420_121 * r_4 * h3_m2) + e_3 * (fs_5145_82654 * h9_m6 - fs_1715_1074502 * h9_m2 + fs_150_41327 * r_2 * h7_m6 + fs_109350_5909761 * r_2 * h7_m2 - fs_5_507 * r_4 * h5_m2 - fs_560_61347 * r_6 * h3_m2);

        pc_32[k] = -fs_14175_64 * e_0 * h3_m3 + e_1 * (-f_15_2 * h5_m3 + fs_1575_16 * r_2 * h3_m3) + e_2 * (fs_2625_1144 * h7_m7 - fs_84375_163592 * h7_m3 + f_30_13 * r_2 * h5_m3 - fs_1575_484 * r_4 * h3_m3) + e_3 * (fs_1372_41327 * h9_m7 - fs_147_537251 * h9_m3 - fs_2625_82654 * r_2 * h7_m7 + fs_84375_11819522 * r_2 * h7_m3 - f_2_13 * r_4 * h5_m3 + fs_175_20449 * r_6 * h3_m3);
    }

    // NOTE: the rows are formed in 59 loops, as the vectorizer runs out of
    // registers with all 99 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph3_p3, ph5_p4, ph5_p5, ph7_p3, ph7_p4, ph7_p5, ph7_p6, ph9_p3, ph9_p4, ph9_p5, ph9_p6, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h3_p3 = ph3_p3[k];
        const auto h5_p4 = ph5_p4[k];
        const auto h5_p5 = ph5_p5[k];
        const auto h7_p3 = ph7_p3[k];
        const auto h7_p4 = ph7_p4[k];
        const auto h7_p5 = ph7_p5[k];
        const auto h7_p6 = ph7_p6[k];
        const auto h9_p3 = ph9_p3[k];
        const auto h9_p4 = ph9_p4[k];
        const auto h9_p5 = ph9_p5[k];
        const auto h9_p6 = ph9_p6[k];

        pc_33[k] = -f_15_2 * e_1 * h5_p4 + e_2 * (-fs_9375_7436 * h7_p4 - fs_3375_1144 * h7_p6 + f_30_13 * r_2 * h5_p4) + e_3 * (-fs_49_41327 * h9_p4 - fs_1029_82654 * h9_p6 + fs_9375_537251 * r_2 * h7_p4 + fs_3375_82654 * r_2 * h7_p6 - f_2_13 * r_4 * h5_p4);

        pc_34[k] = fs_25515_64 * e_0 * h3_p3 + e_1 * (f_15_2 * h5_p5 - fs_2835_16 * r_2 * h3_p3) + e_2 * (-fs_1875_968 * h7_p3 - fs_6075_14872 * h7_p5 - f_30_13 * r_2 * h5_p5 + fs_2835_484 * r_4 * h3_p3) + e_3 * (-fs_2940_537251 * h9_p3 - fs_1372_41327 * h9_p5 + fs_1875_69938 * r_2 * h7_p3 + fs_6075_1074502 * r_2 * h7_p5 + f_2_13 * r_4 * h5_p5 - fs_315_20449 * r_6 * h3_p3);
    }

    // NOTE: the rows are formed in 59 loops, as the vectorizer runs out of
    // registers with all 99 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph3_p2, ph5_p2, ph7_p2, ph7_p4, ph9_p2, ph9_p4, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h3_p2 = ph3_p2[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h7_p4 = ph7_p4[k];
        const auto h9_p2 = ph9_p2[k];
        const auto h9_p4 = ph9_p4[k];

        pc_35[k] = -fs_945_64 * e_0 * h3_p2 + e_1 * (fs_375_16 * h5_p2 + fs_105_16 * r_2 * h3_p2) + e_2 * (-fs_664225_490776 * h7_p2 + fs_25_132 * h7_p4 - fs_375_169 * r_2 * h5_p2 - fs_105_484 * r_4 * h3_p2) + e_3 * (-fs_15435_1074502 * h9_p2 - fs_2205_41327 * h9_p4 + fs_664225_35458566 * r_2 * h7_p2 - fs_25_9537 * r_2 * h7_p4 + fs_5_507 * r_4 * h5_p2 + fs_35_61347 * r_6 * h3_p2);
    }

    // NOTE: the rows are formed in 59 loops, as the vectorizer runs out of
    // registers with all 99 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph1_p1, ph3_p1, ph3_p3, ph5_p1, ph5_p3, ph7_p1, ph7_p3, ph9_p1, ph9_p3, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h1_p1 = ph1_p1[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h3_p3 = ph3_p3[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h7_p3 = ph7_p3[k];
        const auto h9_p1 = ph9_p1[k];
        const auto h9_p3 = ph9_p3[k];

        pc_36[k] = e_0 * (-fs_252 * h3_p1 + fs_23625_64 * h3_p3 - fs_9261_8 * r_2 * h1_p1) + e_1 * (fs_875_32 * h5_p1 - fs_375_16 * h5_p3 + fs_112 * r_2 * h3_p1 - fs_2625_16 * r_2 * h3_p3 + fs_1701_8 * r_4 * h1_p1) + e_2 * (-fs_42025_122694 * h7_p1 + fs_46225_40898 * h7_p3 - fs_875_338 * r_2 * h5_p1 + fs_375_169 * r_2 * h5_p3 - fs_448_121 * r_4 * h3_p1 + fs_2625_484 * r_4 * h3_p3 - fs_14_3 * r_6 * h1_p1) + e_3 * (-fs_164640_5909761 * h9_p1 - fs_35280_537251 * h9_p3 + fs_84050_17729283 * r_2 * h7_p1 - fs_92450_5909761 * r_2 * h7_p3 + fs_35_3042 * r_4 * h5_p1 - fs_5_507 * r_4 * h5_p3 + fs_1792_184041 * r_6 * h3_p1 - fs_875_61347 * r_6 * h3_p3 + fs_7_726 * r_8 * h1_p1) + fs_231525_512 * e_4 * h1_p1;
    }

    // NOTE: the rows are formed in 59 loops, as the vectorizer runs out of
    // registers with all 99 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph1_0, ph3_0, ph3_p2, ph5_0, ph5_p2, ph7_0, ph7_p2, ph9_0, ph9_p2, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h1_0 = ph1_0[k];
        const auto h3_0 = ph3_0[k];
        const auto h3_p2 = ph3_p2[k];
        const auto h5_0 = ph5_0[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h7_0 = ph7_0[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h9_0 = ph9_0[k];
        const auto h9_p2 = ph9_p2[k];

        pc_37[k] = e_0 * (-fs_9747_32 * h3_0 + fs_28125_128 * h3_p2 - fs_2646 * r_2 * h1_0) + e_1 * (fs_75_8 * h5_0 - fs_875_32 * h5_p2 + fs_1083_8 * r_2 * h3_0 - fs_3125_32 * r_2 * h3_p2 + fs_486 * r_4 * h1_0) + e_2 * (fs_2450_61347 * h7_0 + fs_109375_81796 * h7_p2 - fs_150_169 * r_2 * h5_0 + fs_875_338 * r_2 * h5_p2 - fs_1083_242 * r_4 * h3_0 + fs_3125_968 * r_4 * h3_p2 - fs_32_3 * r_6 * h1_0) + e_3 * (-fs_518616_5909761 * h9_0 - fs_36015_537251 * h9_p2 - fs_9800_17729283 * r_2 * h7_0 - fs_109375_5909761 * r_2 * h7_p2 + fs_2_507 * r_4 * h5_0 - fs_35_3042 * r_4 * h5_p2 + fs_722_61347 * r_6 * h3_0 - fs_3125_368082 * r_6 * h3_p2 + fs_8_363 * r_8 * h1_0) + fs_33075_32 * e_4 * h1_0;
    }

    // NOTE: the rows are formed in 59 loops, as the vectorizer runs out of
    // registers with all 99 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph1_m1, ph3_m2, ph3_m1, ph5_m2, ph5_m1, ph7_m2, ph7_m1, ph9_m2, ph9_m1, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h1_m1 = ph1_m1[k];
        const auto h3_m2 = ph3_m2[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h7_m1 = ph7_m1[k];
        const auto h9_m2 = ph9_m2[k];
        const auto h9_m1 = ph9_m1[k];

        pc_38[k] = e_0 * (-fs_135_64 * h3_m1 + fs_2205_2 * r_2 * h1_m1) + e_1 * (fs_75_8 * h5_m1 + fs_15_16 * r_2 * h3_m1 - fs_405_2 * r_4 * h1_m1) + e_2 * (-fs_462875_368082 * h7_m1 - fs_150_169 * r_2 * h5_m1 - fs_15_484 * r_4 * h3_m1 + fs_40_9 * r_6 * h1_m1) + e_3 * (fs_691488_5909761 * h9_m1 + fs_925750_53187849 * r_2 * h7_m1 + fs_2_507 * r_4 * h5_m1 + fs_5_61347 * r_6 * h3_m1 - fs_10_1089 * r_8 * h1_m1) - fs_55125_128 * e_4 * h1_m1;

        pc_39[k] = -fs_28125_128 * e_0 * h3_m2 + e_1 * (fs_875_32 * h5_m2 + fs_3125_32 * r_2 * h3_m2) + e_2 * (-fs_109375_81796 * h7_m2 - fs_875_338 * r_2 * h5_m2 - fs_3125_968 * r_4 * h3_m2) + e_3 * (fs_36015_537251 * h9_m2 + fs_109375_5909761 * r_2 * h7_m2 + fs_35_3042 * r_4 * h5_m2 + fs_3125_368082 * r_6 * h3_m2);
    }

    // NOTE: the rows are formed in 59 loops, as the vectorizer runs out of
    // registers with all 99 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph1_m1, ph3_m3, ph3_m1, ph5_m3, ph5_m1, ph7_m3, ph7_m1, ph9_m3, ph9_m1, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h1_m1 = ph1_m1[k];
        const auto h3_m3 = ph3_m3[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h5_m3 = ph5_m3[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h7_m3 = ph7_m3[k];
        const auto h7_m1 = ph7_m1[k];
        const auto h9_m3 = ph9_m3[k];
        const auto h9_m1 = ph9_m1[k];

        pc_40[k] = e_0 * (-fs_23625_64 * h3_m3 + fs_252 * h3_m1 + fs_9261_8 * r_2 * h1_m1) + e_1 * (fs_375_16 * h5_m3 - fs_875_32 * h5_m1 + fs_2625_16 * r_2 * h3_m3 - fs_112 * r_2 * h3_m1 - fs_1701_8 * r_4 * h1_m1) + e_2 * (-fs_46225_40898 * h7_m3 + fs_42025_122694 * h7_m1 - fs_375_169 * r_2 * h5_m3 + fs_875_338 * r_2 * h5_m1 - fs_2625_484 * r_4 * h3_m3 + fs_448_121 * r_4 * h3_m1 + fs_14_3 * r_6 * h1_m1) + e_3 * (fs_35280_537251 * h9_m3 + fs_164640_5909761 * h9_m1 + fs_92450_5909761 * r_2 * h7_m3 - fs_84050_17729283 * r_2 * h7_m1 + fs_5_507 * r_4 * h5_m3 - fs_35_3042 * r_4 * h5_m1 + fs_875_61347 * r_6 * h3_m3 - fs_1792_184041 * r_6 * h3_m1 - fs_7_726 * r_8 * h1_m1) - fs_231525_512 * e_4 * h1_m1;
    }

    // NOTE: the rows are formed in 59 loops, as the vectorizer runs out of
    // registers with all 99 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph3_m3, ph3_m2, ph5_m5, ph5_m2, ph7_m5, ph7_m4, ph7_m3, ph7_m2, ph9_m5, ph9_m4, ph9_m3, ph9_m2, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h3_m3 = ph3_m3[k];
        const auto h3_m2 = ph3_m2[k];
        const auto h5_m5 = ph5_m5[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h7_m5 = ph7_m5[k];
        const auto h7_m4 = ph7_m4[k];
        const auto h7_m3 = ph7_m3[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h9_m5 = ph9_m5[k];
        const auto h9_m4 = ph9_m4[k];
        const auto h9_m3 = ph9_m3[k];
        const auto h9_m2 = ph9_m2[k];

        pc_41[k] = fs_945_64 * e_0 * h3_m2 + e_1 * (-fs_375_16 * h5_m2 - fs_105_16 * r_2 * h3_m2) + e_2 * (-fs_25_132 * h7_m4 + fs_664225_490776 * h7_m2 + fs_375_169 * r_2 * h5_m2 + fs_105_484 * r_4 * h3_m2) + e_3 * (fs_2205_41327 * h9_m4 + fs_15435_1074502 * h9_m2 + fs_25_9537 * r_2 * h7_m4 - fs_664225_35458566 * r_2 * h7_m2 - fs_5_507 * r_4 * h5_m2 - fs_35_61347 * r_6 * h3_m2);

        pc_42[k] = -fs_25515_64 * e_0 * h3_m3 + e_1 * (-f_15_2 * h5_m5 + fs_2835_16 * r_2 * h3_m3) + e_2 * (fs_6075_14872 * h7_m5 + fs_1875_968 * h7_m3 + f_30_13 * r_2 * h5_m5 - fs_2835_484 * r_4 * h3_m3) + e_3 * (fs_1372_41327 * h9_m5 + fs_2940_537251 * h9_m3 - fs_6075_1074502 * r_2 * h7_m5 - fs_1875_69938 * r_2 * h7_m3 - f_2_13 * r_4 * h5_m5 + fs_315_20449 * r_6 * h3_m3);
    }

    // NOTE: the rows are formed in 59 loops, as the vectorizer runs out of
    // registers with all 99 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph3_m3, ph5_m5, ph5_m4, ph5_m3, ph7_m6, ph7_m5, ph7_m4, ph7_m3, ph9_m6, ph9_m5, ph9_m4, ph9_m3, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h3_m3 = ph3_m3[k];
        const auto h5_m5 = ph5_m5[k];
        const auto h5_m4 = ph5_m4[k];
        const auto h5_m3 = ph5_m3[k];
        const auto h7_m6 = ph7_m6[k];
        const auto h7_m5 = ph7_m5[k];
        const auto h7_m4 = ph7_m4[k];
        const auto h7_m3 = ph7_m3[k];
        const auto h9_m6 = ph9_m6[k];
        const auto h9_m5 = ph9_m5[k];
        const auto h9_m4 = ph9_m4[k];
        const auto h9_m3 = ph9_m3[k];

        pc_43[k] = f_15_2 * e_1 * h5_m4 + e_2 * (fs_3375_1144 * h7_m6 + fs_9375_7436 * h7_m4 - f_30_13 * r_2 * h5_m4) + e_3 * (fs_1029_82654 * h9_m6 + fs_49_41327 * h9_m4 - fs_3375_82654 * r_2 * h7_m6 - fs_9375_537251 * r_2 * h7_m4 + f_2_13 * r_4 * h5_m4);

        pc_44[k] = f_15_2 * e_1 * h5_m5 + e_2 * (fs_16875_3718 * h7_m5 - f_30_13 * r_2 * h5_m5) + e_3 * (fs_343_41327 * h9_m5 - fs_33750_537251 * r_2 * h7_m5 + f_2_13 * r_4 * h5_m5);

        pc_45[k] = -f_15_2 * e_1 * h5_m4 + e_2 * (fs_6000_1859 * h7_m4 + f_30_13 * r_2 * h5_m4) + e_3 * (fs_1225_41327 * h9_m4 - fs_24000_537251 * r_2 * h7_m4 - f_2_13 * r_4 * h5_m4);

        pc_46[k] = fs_14175_16 * e_0 * h3_m3 + e_1 * (-f_15_2 * h5_m3 - fs_1575_4 * r_2 * h3_m3) + e_2 * (fs_66125_122694 * h7_m3 + f_30_13 * r_2 * h5_m3 + fs_1575_121 * r_4 * h3_m3) + e_3 * (fs_33075_537251 * h9_m3 - fs_132250_17729283 * r_2 * h7_m3 - f_2_13 * r_4 * h5_m3 - fs_700_20449 * r_6 * h3_m3);
    }

    // NOTE: the rows are formed in 59 loops, as the vectorizer runs out of
    // registers with all 99 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph1_m1, ph3_m2, ph3_m1, ph5_m2, ph5_m1, ph7_m2, ph7_m1, ph9_m2, ph9_m1, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h1_m1 = ph1_m1[k];
        const auto h3_m2 = ph3_m2[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h7_m1 = ph7_m1[k];
        const auto h9_m2 = ph9_m2[k];
        const auto h9_m1 = ph9_m1[k];

        pc_47[k] = fs_1575_16 * e_0 * h3_m2 + e_1 * (-f_5_4 * h5_m2 - fs_175_4 * r_2 * h3_m2) + e_2 * (-fs_4000_20449 * h7_m2 + f_5_13 * r_2 * h5_m2 + fs_175_121 * r_4 * h3_m2) + e_3 * (fs_51450_537251 * h9_m2 + fs_16000_5909761 * r_2 * h7_m2 - f_1_39 * r_4 * h5_m2 - fs_700_184041 * r_6 * h3_m2);

        pc_48[k] = e_0 * (-fs_5445_32 * h3_m1 - fs_6615_4 * r_2 * h1_m1) + e_1 * (f_5 * h5_m1 + fs_605_8 * r_2 * h3_m1 + fs_1215_4 * r_4 * h1_m1) + e_2 * (-fs_875_507 * h7_m1 - f_20_13 * r_2 * h5_m1 - fs_5_2 * r_4 * h3_m1 - fs_20_3 * r_6 * h1_m1) + e_3 * (fs_720300_5909761 * h9_m1 + fs_3500_146523 * r_2 * h7_m1 + f_4_39 * r_4 * h5_m1 + fs_10_1521 * r_6 * h3_m1 + fs_5_363 * r_8 * h1_m1) + fs_165375_256 * e_4 * h1_m1;
    }

    // NOTE: the rows are formed in 59 loops, as the vectorizer runs out of
    // registers with all 99 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph1_0, ph1_p1, ph3_0, ph3_p1, ph5_0, ph5_p1, ph7_0, ph7_p1, ph9_0, ph9_p1, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h1_0 = ph1_0[k];
        const auto h1_p1 = ph1_p1[k];
        const auto h3_0 = ph3_0[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h5_0 = ph5_0[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h7_0 = ph7_0[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h9_0 = ph9_0[k];
        const auto h9_p1 = ph9_p1[k];

        pc_49[k] = e_0 * (-f_45_2 * h3_0 - f_105_2 * r_2 * h1_0) + e_1 * (f_15_2 * h5_0 + f_15 * r_2 * h3_0 + f_45_2 * r_4 * h1_0) + e_2 * (-f_700_429 * h7_0 - f_30_13 * r_2 * h5_0 - f_30_11 * r_4 * h3_0 - f_10_3 * r_6 * h1_0) + e_3 * (f_882_2431 * h9_0 + f_1400_7293 * r_2 * h7_0 + f_2_13 * r_4 * h5_0 + f_20_143 * r_6 * h3_0 + f_5_33 * r_8 * h1_0) + f_525_16 * e_4 * h1_0;

        pc_50[k] = e_0 * (-fs_5445_32 * h3_p1 - fs_6615_4 * r_2 * h1_p1) + e_1 * (f_5 * h5_p1 + fs_605_8 * r_2 * h3_p1 + fs_1215_4 * r_4 * h1_p1) + e_2 * (-fs_875_507 * h7_p1 - f_20_13 * r_2 * h5_p1 - fs_5_2 * r_4 * h3_p1 - fs_20_3 * r_6 * h1_p1) + e_3 * (fs_720300_5909761 * h9_p1 + fs_3500_146523 * r_2 * h7_p1 + f_4_39 * r_4 * h5_p1 + fs_10_1521 * r_6 * h3_p1 + fs_5_363 * r_8 * h1_p1) + fs_165375_256 * e_4 * h1_p1;
    }

    // NOTE: the rows are formed in 59 loops, as the vectorizer runs out of
    // registers with all 99 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph3_p2, ph3_p3, ph5_p2, ph5_p3, ph5_p4, ph7_p2, ph7_p3, ph7_p4, ph9_p2, ph9_p3, ph9_p4, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h3_p2 = ph3_p2[k];
        const auto h3_p3 = ph3_p3[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h5_p4 = ph5_p4[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h7_p3 = ph7_p3[k];
        const auto h7_p4 = ph7_p4[k];
        const auto h9_p2 = ph9_p2[k];
        const auto h9_p3 = ph9_p3[k];
        const auto h9_p4 = ph9_p4[k];

        pc_51[k] = fs_1575_16 * e_0 * h3_p2 + e_1 * (-f_5_4 * h5_p2 - fs_175_4 * r_2 * h3_p2) + e_2 * (-fs_4000_20449 * h7_p2 + f_5_13 * r_2 * h5_p2 + fs_175_121 * r_4 * h3_p2) + e_3 * (fs_51450_537251 * h9_p2 + fs_16000_5909761 * r_2 * h7_p2 - f_1_39 * r_4 * h5_p2 - fs_700_184041 * r_6 * h3_p2);

        pc_52[k] = fs_14175_16 * e_0 * h3_p3 + e_1 * (-f_15_2 * h5_p3 - fs_1575_4 * r_2 * h3_p3) + e_2 * (fs_66125_122694 * h7_p3 + f_30_13 * r_2 * h5_p3 + fs_1575_121 * r_4 * h3_p3) + e_3 * (fs_33075_537251 * h9_p3 - fs_132250_17729283 * r_2 * h7_p3 - f_2_13 * r_4 * h5_p3 - fs_700_20449 * r_6 * h3_p3);

        pc_53[k] = -f_15_2 * e_1 * h5_p4 + e_2 * (fs_6000_1859 * h7_p4 + f_30_13 * r_2 * h5_p4) + e_3 * (fs_1225_41327 * h9_p4 - fs_24000_537251 * r_2 * h7_p4 - f_2_13 * r_4 * h5_p4);
    }

    // NOTE: the rows are formed in 59 loops, as the vectorizer runs out of
    // registers with all 99 of them in one.

#pragma omp simd aligned(pe_1, pe_2, pe_3, ph5_m4, ph5_p5, ph7_m6, ph7_m4, ph7_p5, ph9_m6, ph9_m4, ph9_p5, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

        const auto h5_m4 = ph5_m4[k];
        const auto h5_p5 = ph5_p5[k];
        const auto h7_m6 = ph7_m6[k];
        const auto h7_m4 = ph7_m4[k];
        const auto h7_p5 = ph7_p5[k];
        const auto h9_m6 = ph9_m6[k];
        const auto h9_m4 = ph9_m4[k];
        const auto h9_p5 = ph9_p5[k];

        pc_54[k] = f_15_2 * e_1 * h5_p5 + e_2 * (fs_16875_3718 * h7_p5 - f_30_13 * r_2 * h5_p5) + e_3 * (fs_343_41327 * h9_p5 - fs_33750_537251 * r_2 * h7_p5 + f_2_13 * r_4 * h5_p5);

        pc_55[k] = -f_15_2 * e_1 * h5_m4 + e_2 * (fs_3375_1144 * h7_m6 - fs_9375_7436 * h7_m4 + f_30_13 * r_2 * h5_m4) + e_3 * (fs_1029_82654 * h9_m6 - fs_49_41327 * h9_m4 - fs_3375_82654 * r_2 * h7_m6 + fs_9375_537251 * r_2 * h7_m4 - f_2_13 * r_4 * h5_m4);
    }

    // NOTE: the rows are formed in 59 loops, as the vectorizer runs out of
    // registers with all 99 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph3_m3, ph3_m2, ph5_m5, ph5_m2, ph7_m5, ph7_m4, ph7_m3, ph7_m2, ph9_m5, ph9_m4, ph9_m3, ph9_m2, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h3_m3 = ph3_m3[k];
        const auto h3_m2 = ph3_m2[k];
        const auto h5_m5 = ph5_m5[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h7_m5 = ph7_m5[k];
        const auto h7_m4 = ph7_m4[k];
        const auto h7_m3 = ph7_m3[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h9_m5 = ph9_m5[k];
        const auto h9_m4 = ph9_m4[k];
        const auto h9_m3 = ph9_m3[k];
        const auto h9_m2 = ph9_m2[k];

        pc_56[k] = fs_25515_64 * e_0 * h3_m3 + e_1 * (-f_15_2 * h5_m5 - fs_2835_16 * r_2 * h3_m3) + e_2 * (fs_6075_14872 * h7_m5 - fs_1875_968 * h7_m3 + f_30_13 * r_2 * h5_m5 + fs_2835_484 * r_4 * h3_m3) + e_3 * (fs_1372_41327 * h9_m5 - fs_2940_537251 * h9_m3 - fs_6075_1074502 * r_2 * h7_m5 + fs_1875_69938 * r_2 * h7_m3 - f_2_13 * r_4 * h5_m5 - fs_315_20449 * r_6 * h3_m3);

        pc_57[k] = -fs_945_64 * e_0 * h3_m2 + e_1 * (fs_375_16 * h5_m2 + fs_105_16 * r_2 * h3_m2) + e_2 * (-fs_25_132 * h7_m4 - fs_664225_490776 * h7_m2 - fs_375_169 * r_2 * h5_m2 - fs_105_484 * r_4 * h3_m2) + e_3 * (fs_2205_41327 * h9_m4 - fs_15435_1074502 * h9_m2 + fs_25_9537 * r_2 * h7_m4 + fs_664225_35458566 * r_2 * h7_m2 + fs_5_507 * r_4 * h5_m2 + fs_35_61347 * r_6 * h3_m2);
    }

    // NOTE: the rows are formed in 59 loops, as the vectorizer runs out of
    // registers with all 99 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph1_m1, ph3_m3, ph3_m1, ph5_m3, ph5_m1, ph7_m3, ph7_m1, ph9_m3, ph9_m1, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h1_m1 = ph1_m1[k];
        const auto h3_m3 = ph3_m3[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h5_m3 = ph5_m3[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h7_m3 = ph7_m3[k];
        const auto h7_m1 = ph7_m1[k];
        const auto h9_m3 = ph9_m3[k];
        const auto h9_m1 = ph9_m1[k];

        pc_58[k] = e_0 * (-fs_23625_64 * h3_m3 - fs_252 * h3_m1 - fs_9261_8 * r_2 * h1_m1) + e_1 * (fs_375_16 * h5_m3 + fs_875_32 * h5_m1 + fs_2625_16 * r_2 * h3_m3 + fs_112 * r_2 * h3_m1 + fs_1701_8 * r_4 * h1_m1) + e_2 * (-fs_46225_40898 * h7_m3 - fs_42025_122694 * h7_m1 - fs_375_169 * r_2 * h5_m3 - fs_875_338 * r_2 * h5_m1 - fs_2625_484 * r_4 * h3_m3 - fs_448_121 * r_4 * h3_m1 - fs_14_3 * r_6 * h1_m1) + e_3 * (fs_35280_537251 * h9_m3 - fs_164640_5909761 * h9_m1 + fs_92450_5909761 * r_2 * h7_m3 + fs_84050_17729283 * r_2 * h7_m1 + fs_5_507 * r_4 * h5_m3 + fs_35_3042 * r_4 * h5_m1 + fs_875_61347 * r_6 * h3_m3 + fs_1792_184041 * r_6 * h3_m1 + fs_7_726 * r_8 * h1_m1) + fs_231525_512 * e_4 * h1_m1;
    }

    // NOTE: the rows are formed in 59 loops, as the vectorizer runs out of
    // registers with all 99 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph1_p1, ph3_m2, ph3_p1, ph5_m2, ph5_p1, ph7_m2, ph7_p1, ph9_m2, ph9_p1, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h1_p1 = ph1_p1[k];
        const auto h3_m2 = ph3_m2[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h9_m2 = ph9_m2[k];
        const auto h9_p1 = ph9_p1[k];

        pc_59[k] = -fs_28125_128 * e_0 * h3_m2 + e_1 * (fs_875_32 * h5_m2 + fs_3125_32 * r_2 * h3_m2) + e_2 * (-fs_109375_81796 * h7_m2 - fs_875_338 * r_2 * h5_m2 - fs_3125_968 * r_4 * h3_m2) + e_3 * (fs_36015_537251 * h9_m2 + fs_109375_5909761 * r_2 * h7_m2 + fs_35_3042 * r_4 * h5_m2 + fs_3125_368082 * r_6 * h3_m2);

        pc_60[k] = e_0 * (-fs_135_64 * h3_p1 + fs_2205_2 * r_2 * h1_p1) + e_1 * (fs_75_8 * h5_p1 + fs_15_16 * r_2 * h3_p1 - fs_405_2 * r_4 * h1_p1) + e_2 * (-fs_462875_368082 * h7_p1 - fs_150_169 * r_2 * h5_p1 - fs_15_484 * r_4 * h3_p1 + fs_40_9 * r_6 * h1_p1) + e_3 * (fs_691488_5909761 * h9_p1 + fs_925750_53187849 * r_2 * h7_p1 + fs_2_507 * r_4 * h5_p1 + fs_5_61347 * r_6 * h3_p1 - fs_10_1089 * r_8 * h1_p1) - fs_55125_128 * e_4 * h1_p1;
    }

    // NOTE: the rows are formed in 59 loops, as the vectorizer runs out of
    // registers with all 99 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph1_0, ph3_0, ph3_p2, ph5_0, ph5_p2, ph7_0, ph7_p2, ph9_0, ph9_p2, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h1_0 = ph1_0[k];
        const auto h3_0 = ph3_0[k];
        const auto h3_p2 = ph3_p2[k];
        const auto h5_0 = ph5_0[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h7_0 = ph7_0[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h9_0 = ph9_0[k];
        const auto h9_p2 = ph9_p2[k];

        pc_61[k] = e_0 * (-fs_9747_32 * h3_0 - fs_28125_128 * h3_p2 - fs_2646 * r_2 * h1_0) + e_1 * (fs_75_8 * h5_0 + fs_875_32 * h5_p2 + fs_1083_8 * r_2 * h3_0 + fs_3125_32 * r_2 * h3_p2 + fs_486 * r_4 * h1_0) + e_2 * (fs_2450_61347 * h7_0 - fs_109375_81796 * h7_p2 - fs_150_169 * r_2 * h5_0 - fs_875_338 * r_2 * h5_p2 - fs_1083_242 * r_4 * h3_0 - fs_3125_968 * r_4 * h3_p2 - fs_32_3 * r_6 * h1_0) + e_3 * (-fs_518616_5909761 * h9_0 + fs_36015_537251 * h9_p2 - fs_9800_17729283 * r_2 * h7_0 + fs_109375_5909761 * r_2 * h7_p2 + fs_2_507 * r_4 * h5_0 + fs_35_3042 * r_4 * h5_p2 + fs_722_61347 * r_6 * h3_0 + fs_3125_368082 * r_6 * h3_p2 + fs_8_363 * r_8 * h1_0) + fs_33075_32 * e_4 * h1_0;
    }

    // NOTE: the rows are formed in 59 loops, as the vectorizer runs out of
    // registers with all 99 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph1_p1, ph3_p1, ph3_p3, ph5_p1, ph5_p3, ph7_p1, ph7_p3, ph9_p1, ph9_p3, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h1_p1 = ph1_p1[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h3_p3 = ph3_p3[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h7_p3 = ph7_p3[k];
        const auto h9_p1 = ph9_p1[k];
        const auto h9_p3 = ph9_p3[k];

        pc_62[k] = e_0 * (-fs_252 * h3_p1 - fs_23625_64 * h3_p3 - fs_9261_8 * r_2 * h1_p1) + e_1 * (fs_875_32 * h5_p1 + fs_375_16 * h5_p3 + fs_112 * r_2 * h3_p1 + fs_2625_16 * r_2 * h3_p3 + fs_1701_8 * r_4 * h1_p1) + e_2 * (-fs_42025_122694 * h7_p1 - fs_46225_40898 * h7_p3 - fs_875_338 * r_2 * h5_p1 - fs_375_169 * r_2 * h5_p3 - fs_448_121 * r_4 * h3_p1 - fs_2625_484 * r_4 * h3_p3 - fs_14_3 * r_6 * h1_p1) + e_3 * (-fs_164640_5909761 * h9_p1 + fs_35280_537251 * h9_p3 + fs_84050_17729283 * r_2 * h7_p1 + fs_92450_5909761 * r_2 * h7_p3 + fs_35_3042 * r_4 * h5_p1 + fs_5_507 * r_4 * h5_p3 + fs_1792_184041 * r_6 * h3_p1 + fs_875_61347 * r_6 * h3_p3 + fs_7_726 * r_8 * h1_p1) + fs_231525_512 * e_4 * h1_p1;
    }

    // NOTE: the rows are formed in 59 loops, as the vectorizer runs out of
    // registers with all 99 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph3_p2, ph3_p3, ph5_p2, ph5_p5, ph7_p2, ph7_p3, ph7_p4, ph7_p5, ph9_p2, ph9_p3, ph9_p4, ph9_p5, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h3_p2 = ph3_p2[k];
        const auto h3_p3 = ph3_p3[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h5_p5 = ph5_p5[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h7_p3 = ph7_p3[k];
        const auto h7_p4 = ph7_p4[k];
        const auto h7_p5 = ph7_p5[k];
        const auto h9_p2 = ph9_p2[k];
        const auto h9_p3 = ph9_p3[k];
        const auto h9_p4 = ph9_p4[k];
        const auto h9_p5 = ph9_p5[k];

        pc_63[k] = -fs_945_64 * e_0 * h3_p2 + e_1 * (fs_375_16 * h5_p2 + fs_105_16 * r_2 * h3_p2) + e_2 * (-fs_664225_490776 * h7_p2 - fs_25_132 * h7_p4 - fs_375_169 * r_2 * h5_p2 - fs_105_484 * r_4 * h3_p2) + e_3 * (-fs_15435_1074502 * h9_p2 + fs_2205_41327 * h9_p4 + fs_664225_35458566 * r_2 * h7_p2 + fs_25_9537 * r_2 * h7_p4 + fs_5_507 * r_4 * h5_p2 + fs_35_61347 * r_6 * h3_p2);

        pc_64[k] = fs_25515_64 * e_0 * h3_p3 + e_1 * (-f_15_2 * h5_p5 - fs_2835_16 * r_2 * h3_p3) + e_2 * (-fs_1875_968 * h7_p3 + fs_6075_14872 * h7_p5 + f_30_13 * r_2 * h5_p5 + fs_2835_484 * r_4 * h3_p3) + e_3 * (-fs_2940_537251 * h9_p3 + fs_1372_41327 * h9_p5 + fs_1875_69938 * r_2 * h7_p3 - fs_6075_1074502 * r_2 * h7_p5 - f_2_13 * r_4 * h5_p5 - fs_315_20449 * r_6 * h3_p3);
    }

    // NOTE: the rows are formed in 59 loops, as the vectorizer runs out of
    // registers with all 99 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph3_m3, ph5_m3, ph5_p4, ph7_m7, ph7_m3, ph7_p4, ph7_p6, ph9_m7, ph9_m3, ph9_p4, ph9_p6, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h3_m3 = ph3_m3[k];
        const auto h5_m3 = ph5_m3[k];
        const auto h5_p4 = ph5_p4[k];
        const auto h7_m7 = ph7_m7[k];
        const auto h7_m3 = ph7_m3[k];
        const auto h7_p4 = ph7_p4[k];
        const auto h7_p6 = ph7_p6[k];
        const auto h9_m7 = ph9_m7[k];
        const auto h9_m3 = ph9_m3[k];
        const auto h9_p4 = ph9_p4[k];
        const auto h9_p6 = ph9_p6[k];

        pc_65[k] = -f_15_2 * e_1 * h5_p4 + e_2 * (-fs_9375_7436 * h7_p4 + fs_3375_1144 * h7_p6 + f_30_13 * r_2 * h5_p4) + e_3 * (-fs_49_41327 * h9_p4 + fs_1029_82654 * h9_p6 + fs_9375_537251 * r_2 * h7_p4 - fs_3375_82654 * r_2 * h7_p6 - f_2_13 * r_4 * h5_p4);

        pc_66[k] = fs_14175_64 * e_0 * h3_m3 + e_1 * (f_15_2 * h5_m3 - fs_1575_16 * r_2 * h3_m3) + e_2 * (fs_2625_1144 * h7_m7 + fs_84375_163592 * h7_m3 - f_30_13 * r_2 * h5_m3 + fs_1575_484 * r_4 * h3_m3) + e_3 * (fs_1372_41327 * h9_m7 + fs_147_537251 * h9_m3 - fs_2625_82654 * r_2 * h7_m7 - fs_84375_11819522 * r_2 * h7_m3 + f_2_13 * r_4 * h5_m3 - fs_175_20449 * r_6 * h3_m3);
    }

    // NOTE: the rows are formed in 59 loops, as the vectorizer runs out of
    // registers with all 99 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph3_m2, ph5_m2, ph7_m6, ph7_m2, ph9_m6, ph9_m2, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h3_m2 = ph3_m2[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h7_m6 = ph7_m6[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h9_m6 = ph9_m6[k];
        const auto h9_m2 = ph9_m2[k];

        pc_67[k] = -fs_945_4 * e_0 * h3_m2 + e_1 * (fs_375_16 * h5_m2 + fs_105 * r_2 * h3_m2) + e_2 * (-fs_75_286 * h7_m6 + fs_54675_40898 * h7_m2 - fs_375_169 * r_2 * h5_m2 - fs_420_121 * r_4 * h3_m2) + e_3 * (fs_5145_82654 * h9_m6 + fs_1715_1074502 * h9_m2 + fs_150_41327 * r_2 * h7_m6 - fs_109350_5909761 * r_2 * h7_m2 + fs_5_507 * r_4 * h5_m2 + fs_560_61347 * r_6 * h3_m2);
    }

    // NOTE: the rows are formed in 59 loops, as the vectorizer runs out of
    // registers with all 99 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph1_m1, ph3_m1, ph5_m5, ph5_m4, ph7_m5, ph7_m4, ph7_m1, ph9_m5, ph9_m4, ph9_m1, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h1_m1 = ph1_m1[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h5_m5 = ph5_m5[k];
        const auto h5_m4 = ph5_m4[k];
        const auto h7_m5 = ph7_m5[k];
        const auto h7_m4 = ph7_m4[k];
        const auto h7_m1 = ph7_m1[k];
        const auto h9_m5 = ph9_m5[k];
        const auto h9_m4 = ph9_m4[k];
        const auto h9_m1 = ph9_m1[k];

        pc_68[k] = e_0 * (-fs_15309_64 * h3_m1 - fs_3087_2 * r_2 * h1_m1) + e_1 * (f_15_2 * h5_m5 + fs_1701_16 * r_2 * h3_m1 + fs_567_2 * r_4 * h1_m1) + e_2 * (-fs_70225_44616 * h7_m5 + fs_15625_8712 * h7_m1 - f_30_13 * r_2 * h5_m5 - fs_1701_484 * r_4 * h3_m1 - fs_56_9 * r_6 * h1_m1) + e_3 * (fs_3087_41327 * h9_m5 + fs_30870_5909761 * h9_m1 + fs_70225_3223506 * r_2 * h7_m5 - fs_15625_629442 * r_2 * h7_m1 + f_2_13 * r_4 * h5_m5 + fs_189_20449 * r_6 * h3_m1 + fs_14_1089 * r_8 * h1_m1) + fs_77175_128 * e_4 * h1_m1;

        pc_69[k] = fs_375_16 * e_1 * h5_m4 + e_2 * (-fs_2500_1859 * h7_m4 - fs_375_169 * r_2 * h5_m4) + e_3 * (fs_2940_41327 * h9_m4 + fs_10000_537251 * r_2 * h7_m4 + fs_5_507 * r_4 * h5_m4);
    }

    // NOTE: the rows are formed in 59 loops, as the vectorizer runs out of
    // registers with all 99 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph1_m1, ph3_m3, ph3_m1, ph3_p2, ph5_m1, ph5_p2, ph7_m3, ph7_m1, ph7_p2, ph9_m3, ph9_m1, ph9_p2, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h1_m1 = ph1_m1[k];
        const auto h3_m3 = ph3_m3[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h3_p2 = ph3_p2[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h7_m3 = ph7_m3[k];
        const auto h7_m1 = ph7_m1[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h9_m3 = ph9_m3[k];
        const auto h9_m1 = ph9_m1[k];
        const auto h9_p2 = ph9_p2[k];

        pc_70[k] = e_0 * (fs_30375_128 * h3_m3 + fs_16641_128 * h3_m1 - fs_1323_4 * r_2 * h1_m1) + e_1 * (-fs_125_4 * h5_m1 - fs_3375_32 * r_2 * h3_m3 - fs_1849_32 * r_2 * h3_m1 + fs_243_4 * r_4 * h1_m1) + e_2 * (-fs_175_484 * h7_m3 + fs_147175_245388 * h7_m1 + fs_500_169 * r_2 * h5_m1 + fs_3375_968 * r_4 * h3_m3 + fs_1849_968 * r_4 * h3_m1 - fs_4_3 * r_6 * h1_m1) + e_3 * (fs_30870_537251 * h9_m3 + fs_144060_5909761 * h9_m1 + fs_175_34969 * r_2 * h7_m3 - fs_147175_17729283 * r_2 * h7_m1 - fs_20_1521 * r_4 * h5_m1 - fs_375_40898 * r_6 * h3_m3 - fs_1849_368082 * r_6 * h3_m1 + fs_1_363 * r_8 * h1_m1) + fs_33075_256 * e_4 * h1_m1;

        pc_71[k] = fs_675 * e_0 * h3_p2 + e_1 * (-fs_525_16 * h5_p2 - fs_300 * r_2 * h3_p2) + e_2 * (fs_1750_61347 * h7_p2 + fs_525_169 * r_2 * h5_p2 + fs_1200_121 * r_4 * h3_p2) + e_3 * (fs_43218_537251 * h9_p2 - fs_7000_17729283 * r_2 * h7_p2 - fs_7_507 * r_4 * h5_p2 - fs_1600_61347 * r_6 * h3_p2);
    }

    // NOTE: the rows are formed in 59 loops, as the vectorizer runs out of
    // registers with all 99 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph1_p1, ph3_p1, ph3_p3, ph5_p1, ph7_p1, ph7_p3, ph9_p1, ph9_p3, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h1_p1 = ph1_p1[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h3_p3 = ph3_p3[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h7_p3 = ph7_p3[k];
        const auto h9_p1 = ph9_p1[k];
        const auto h9_p3 = ph9_p3[k];

        pc_72[k] = e_0 * (-fs_16641_128 * h3_p1 + fs_30375_128 * h3_p3 + fs_1323_4 * r_2 * h1_p1) + e_1 * (fs_125_4 * h5_p1 + fs_1849_32 * r_2 * h3_p1 - fs_3375_32 * r_2 * h3_p3 - fs_243_4 * r_4 * h1_p1) + e_2 * (-fs_147175_245388 * h7_p1 - fs_175_484 * h7_p3 - fs_500_169 * r_2 * h5_p1 - fs_1849_968 * r_4 * h3_p1 + fs_3375_968 * r_4 * h3_p3 + fs_4_3 * r_6 * h1_p1) + e_3 * (-fs_144060_5909761 * h9_p1 + fs_30870_537251 * h9_p3 + fs_147175_17729283 * r_2 * h7_p1 + fs_175_34969 * r_2 * h7_p3 + fs_20_1521 * r_4 * h5_p1 + fs_1849_368082 * r_6 * h3_p1 - fs_375_40898 * r_6 * h3_p3 - fs_1_363 * r_8 * h1_p1) - fs_33075_256 * e_4 * h1_p1;
    }

    // NOTE: the rows are formed in 59 loops, as the vectorizer runs out of
    // registers with all 99 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph1_0, ph3_0, ph5_0, ph5_p4, ph7_0, ph7_p4, ph9_0, ph9_p4, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h1_0 = ph1_0[k];
        const auto h3_0 = ph3_0[k];
        const auto h5_0 = ph5_0[k];
        const auto h5_p4 = ph5_p4[k];
        const auto h7_0 = ph7_0[k];
        const auto h7_p4 = ph7_p4[k];
        const auto h9_0 = ph9_0[k];
        const auto h9_p4 = ph9_p4[k];

        pc_73[k] = e_0 * (-fs_189_16 * h3_0 - fs_9261_4 * r_2 * h1_0) + e_1 * (-fs_525_16 * h5_0 + fs_375_16 * h5_p4 + fs_21_4 * r_2 * h3_0 + fs_1701_4 * r_4 * h1_0) + e_2 * (fs_179200_61347 * h7_0 - fs_2500_1859 * h7_p4 + fs_525_169 * r_2 * h5_0 - fs_375_169 * r_2 * h5_p4 - fs_21_121 * r_4 * h3_0 - fs_28_3 * r_6 * h1_0) + e_3 * (fs_148176_5909761 * h9_0 + fs_2940_41327 * h9_p4 - fs_716800_17729283 * r_2 * h7_0 + fs_10000_537251 * r_2 * h7_p4 - fs_7_507 * r_4 * h5_0 + fs_5_507 * r_4 * h5_p4 + fs_28_61347 * r_6 * h3_0 + fs_7_363 * r_8 * h1_0) + fs_231525_256 * e_4 * h1_0;
    }

    // NOTE: the rows are formed in 59 loops, as the vectorizer runs out of
    // registers with all 99 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph1_p1, ph3_p1, ph5_p5, ph7_p1, ph7_p5, ph9_p1, ph9_p5, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h1_p1 = ph1_p1[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h5_p5 = ph5_p5[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h7_p5 = ph7_p5[k];
        const auto h9_p1 = ph9_p1[k];
        const auto h9_p5 = ph9_p5[k];

        pc_74[k] = e_0 * (-fs_15309_64 * h3_p1 - fs_3087_2 * r_2 * h1_p1) + e_1 * (f_15_2 * h5_p5 + fs_1701_16 * r_2 * h3_p1 + fs_567_2 * r_4 * h1_p1) + e_2 * (fs_15625_8712 * h7_p1 - fs_70225_44616 * h7_p5 - f_30_13 * r_2 * h5_p5 - fs_1701_484 * r_4 * h3_p1 - fs_56_9 * r_6 * h1_p1) + e_3 * (fs_30870_5909761 * h9_p1 + fs_3087_41327 * h9_p5 - fs_15625_629442 * r_2 * h7_p1 + fs_70225_3223506 * r_2 * h7_p5 + f_2_13 * r_4 * h5_p5 + fs_189_20449 * r_6 * h3_p1 + fs_14_1089 * r_8 * h1_p1) + fs_77175_128 * e_4 * h1_p1;
    }

    // NOTE: the rows are formed in 59 loops, as the vectorizer runs out of
    // registers with all 99 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph3_p2, ph3_p3, ph5_p2, ph5_p3, ph7_p2, ph7_p3, ph7_p6, ph7_p7, ph9_p2, ph9_p3, ph9_p6, ph9_p7, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h3_p2 = ph3_p2[k];
        const auto h3_p3 = ph3_p3[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h7_p3 = ph7_p3[k];
        const auto h7_p6 = ph7_p6[k];
        const auto h7_p7 = ph7_p7[k];
        const auto h9_p2 = ph9_p2[k];
        const auto h9_p3 = ph9_p3[k];
        const auto h9_p6 = ph9_p6[k];
        const auto h9_p7 = ph9_p7[k];

        pc_75[k] = -fs_945_4 * e_0 * h3_p2 + e_1 * (fs_375_16 * h5_p2 + fs_105 * r_2 * h3_p2) + e_2 * (fs_54675_40898 * h7_p2 - fs_75_286 * h7_p6 - fs_375_169 * r_2 * h5_p2 - fs_420_121 * r_4 * h3_p2) + e_3 * (fs_1715_1074502 * h9_p2 + fs_5145_82654 * h9_p6 - fs_109350_5909761 * r_2 * h7_p2 + fs_150_41327 * r_2 * h7_p6 + fs_5_507 * r_4 * h5_p2 + fs_560_61347 * r_6 * h3_p2);

        pc_76[k] = fs_14175_64 * e_0 * h3_p3 + e_1 * (f_15_2 * h5_p3 - fs_1575_16 * r_2 * h3_p3) + e_2 * (fs_84375_163592 * h7_p3 + fs_2625_1144 * h7_p7 - f_30_13 * r_2 * h5_p3 + fs_1575_484 * r_4 * h3_p3) + e_3 * (fs_147_537251 * h9_p3 + fs_1372_41327 * h9_p7 - fs_84375_11819522 * r_2 * h7_p3 - fs_2625_82654 * r_2 * h7_p7 + f_2_13 * r_4 * h5_p3 - fs_175_20449 * r_6 * h3_p3);
    }

    // NOTE: the rows are formed in 59 loops, as the vectorizer runs out of
    // registers with all 99 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph1_m1, ph3_m2, ph3_m1, ph5_m2, ph5_m1, ph7_m7, ph7_m2, ph7_m1, ph9_m8, ph9_m7, ph9_m2, ph9_m1, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h1_m1 = ph1_m1[k];
        const auto h3_m2 = ph3_m2[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h7_m7 = ph7_m7[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h7_m1 = ph7_m1[k];
        const auto h9_m8 = ph9_m8[k];
        const auto h9_m7 = ph9_m7[k];
        const auto h9_m2 = ph9_m2[k];
        const auto h9_m1 = ph9_m1[k];

        pc_77[k] = -fs_33075_64 * e_0 * h3_m2 + e_1 * (-fs_525_16 * h5_m2 + fs_3675_16 * r_2 * h3_m2) + e_2 * (-fs_23625_163592 * h7_m2 + fs_525_169 * r_2 * h5_m2 - fs_3675_484 * r_4 * h3_m2) + e_3 * (fs_196_2431 * h9_m8 - fs_49_1074502 * h9_m2 + fs_23625_11819522 * r_2 * h7_m2 - fs_7_507 * r_4 * h5_m2 + fs_1225_61347 * r_6 * h3_m2);

        pc_78[k] = e_0 * (-fs_1323_64 * h3_m1 - fs_3969_2 * r_2 * h1_m1) + e_1 * (-fs_375_8 * h5_m1 + fs_147_16 * r_2 * h3_m1 + fs_729_2 * r_4 * h1_m1) + e_2 * (-fs_3675_1144 * h7_m7 - fs_92575_163592 * h7_m1 + fs_750_169 * r_2 * h5_m1 - fs_147_484 * r_4 * h3_m1 - fs_8 * r_6 * h1_m1) + e_3 * (fs_3920_41327 * h9_m7 - fs_1960_5909761 * h9_m1 + fs_3675_82654 * r_2 * h7_m7 + fs_92575_11819522 * r_2 * h7_m1 - fs_10_507 * r_4 * h5_m1 + fs_49_61347 * r_6 * h3_m1 + fs_2_121 * r_8 * h1_m1) + fs_99225_128 * e_4 * h1_m1;
    }

    // NOTE: the rows are formed in 59 loops, as the vectorizer runs out of
    // registers with all 99 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph1_m1, ph3_m1, ph5_m5, ph5_m1, ph7_m6, ph7_m5, ph7_m1, ph9_m6, ph9_m5, ph9_m1, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h1_m1 = ph1_m1[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h5_m5 = ph5_m5[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h7_m6 = ph7_m6[k];
        const auto h7_m5 = ph7_m5[k];
        const auto h7_m1 = ph7_m1[k];
        const auto h9_m6 = ph9_m6[k];
        const auto h9_m5 = ph9_m5[k];
        const auto h9_m1 = ph9_m1[k];

        pc_79[k] = -fs_4375_3432 * e_2 * h7_m6 + e_3 * (fs_6615_82654 * h9_m6 + fs_4375_247962 * r_2 * h7_m6);

        pc_80[k] = e_0 * (f_147_8 * h3_m1 - fs_1323_8 * r_2 * h1_m1) + e_1 * (-fs_525_16 * h5_m5 - fs_125_32 * h5_m1 - f_49_4 * r_2 * h3_m1 + fs_243_8 * r_4 * h1_m1) + e_2 * (-fs_175_3718 * h7_m5 - fs_214375_122694 * h7_m1 + fs_525_169 * r_2 * h5_m5 + fs_125_338 * r_2 * h5_m1 + f_49_22 * r_4 * h3_m1 - fs_2_3 * r_6 * h1_m1) + e_3 * (fs_2352_41327 * h9_m5 - fs_23520_5909761 * h9_m1 + fs_350_537251 * r_2 * h7_m5 + fs_428750_17729283 * r_2 * h7_m1 - fs_7_507 * r_4 * h5_m5 - fs_5_3042 * r_4 * h5_m1 - f_49_429 * r_6 * h3_m1 + fs_1_726 * r_8 * h1_m1) + fs_33075_512 * e_4 * h1_m1;
    }

    // NOTE: the rows are formed in 59 loops, as the vectorizer runs out of
    // registers with all 99 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph3_m2, ph3_p3, ph5_m4, ph5_m2, ph5_p3, ph7_m4, ph7_m2, ph7_p3, ph9_m4, ph9_m2, ph9_p3, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h3_m2 = ph3_m2[k];
        const auto h3_p3 = ph3_p3[k];
        const auto h5_m4 = ph5_m4[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h7_m4 = ph7_m4[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h7_p3 = ph7_p3[k];
        const auto h9_m4 = ph9_m4[k];
        const auto h9_m2 = ph9_m2[k];
        const auto h9_p3 = ph9_p3[k];

        pc_81[k] = -fs_38115_128 * e_0 * h3_m2 + e_1 * (-fs_375_8 * h5_m4 - fs_125_32 * h5_m2 + fs_4235_32 * r_2 * h3_m2) + e_2 * (fs_1225_3718 * h7_m4 + f_35_26 * h7_m2 + fs_750_169 * r_2 * h5_m4 + fs_125_338 * r_2 * h5_m2 - fs_35_8 * r_4 * h3_m2) + e_3 * (fs_1470_41327 * h9_m4 + fs_5145_537251 * h9_m2 - fs_2450_537251 * r_2 * h7_m4 - f_35_221 * r_2 * h7_m2 - fs_10_507 * r_4 * h5_m4 - fs_5_3042 * r_4 * h5_m2 + fs_35_3042 * r_6 * h3_m2);

        pc_82[k] = -fs_14175_64 * e_0 * h3_p3 + e_1 * (-f_15_2 * h5_p3 + fs_1575_16 * r_2 * h3_p3) + e_2 * (fs_300125_122694 * h7_p3 + f_30_13 * r_2 * h5_p3 - fs_1575_484 * r_4 * h3_p3) + e_3 * (fs_21168_537251 * h9_p3 - fs_600250_17729283 * r_2 * h7_p3 - f_2_13 * r_4 * h5_p3 + fs_175_20449 * r_6 * h3_p3);
    }

    // NOTE: the rows are formed in 59 loops, as the vectorizer runs out of
    // registers with all 99 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph3_p2, ph5_p2, ph5_p4, ph7_p2, ph7_p4, ph9_p2, ph9_p4, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h3_p2 = ph3_p2[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h5_p4 = ph5_p4[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h7_p4 = ph7_p4[k];
        const auto h9_p2 = ph9_p2[k];
        const auto h9_p4 = ph9_p4[k];

        pc_83[k] = fs_38115_128 * e_0 * h3_p2 + e_1 * (fs_125_32 * h5_p2 - fs_375_8 * h5_p4 - fs_4235_32 * r_2 * h3_p2) + e_2 * (-f_35_26 * h7_p2 + fs_1225_3718 * h7_p4 - fs_125_338 * r_2 * h5_p2 + fs_750_169 * r_2 * h5_p4 + fs_35_8 * r_4 * h3_p2) + e_3 * (-fs_5145_537251 * h9_p2 + fs_1470_41327 * h9_p4 + f_35_221 * r_2 * h7_p2 - fs_2450_537251 * r_2 * h7_p4 + fs_5_3042 * r_4 * h5_p2 - fs_10_507 * r_4 * h5_p4 - fs_35_3042 * r_6 * h3_p2);
    }

    // NOTE: the rows are formed in 59 loops, as the vectorizer runs out of
    // registers with all 99 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph1_p1, ph3_p1, ph5_p1, ph5_p5, ph7_p1, ph7_p5, ph9_p1, ph9_p5, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h1_p1 = ph1_p1[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h5_p5 = ph5_p5[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h7_p5 = ph7_p5[k];
        const auto h9_p1 = ph9_p1[k];
        const auto h9_p5 = ph9_p5[k];

        pc_84[k] = e_0 * (-f_147_8 * h3_p1 + fs_1323_8 * r_2 * h1_p1) + e_1 * (fs_125_32 * h5_p1 - fs_525_16 * h5_p5 + f_49_4 * r_2 * h3_p1 - fs_243_8 * r_4 * h1_p1) + e_2 * (fs_214375_122694 * h7_p1 - fs_175_3718 * h7_p5 - fs_125_338 * r_2 * h5_p1 + fs_525_169 * r_2 * h5_p5 - f_49_22 * r_4 * h3_p1 + fs_2_3 * r_6 * h1_p1) + e_3 * (fs_23520_5909761 * h9_p1 + fs_2352_41327 * h9_p5 - fs_428750_17729283 * r_2 * h7_p1 + fs_350_537251 * r_2 * h7_p5 + fs_5_3042 * r_4 * h5_p1 - fs_7_507 * r_4 * h5_p5 + f_49_429 * r_6 * h3_p1 - fs_1_726 * r_8 * h1_p1) - fs_33075_512 * e_4 * h1_p1;
    }

    // NOTE: the rows are formed in 59 loops, as the vectorizer runs out of
    // registers with all 99 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph1_0, ph3_0, ph5_0, ph7_0, ph7_p6, ph9_0, ph9_p6, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h1_0 = ph1_0[k];
        const auto h3_0 = ph3_0[k];
        const auto h5_0 = ph5_0[k];
        const auto h7_0 = ph7_0[k];
        const auto h7_p6 = ph7_p6[k];
        const auto h9_0 = ph9_0[k];
        const auto h9_p6 = ph9_p6[k];

        pc_85[k] = e_0 * (f_63_4 * h3_0 - f_42 * r_2 * h1_0) + e_1 * (-f_15_2 * h5_0 - f_21_2 * r_2 * h3_0 + f_18 * r_4 * h1_0) + e_2 * (-f_665_429 * h7_0 - fs_4375_3432 * h7_p6 + f_30_13 * r_2 * h5_0 + f_21_11 * r_4 * h3_0 - f_8_3 * r_6 * h1_0) + e_3 * (-f_126_2431 * h9_0 + fs_6615_82654 * h9_p6 + f_1330_7293 * r_2 * h7_0 + fs_4375_247962 * r_2 * h7_p6 - f_2_13 * r_4 * h5_0 - f_14_143 * r_6 * h3_0 + f_4_33 * r_8 * h1_0) + f_105_4 * e_4 * h1_0;
    }

    // NOTE: the rows are formed in 59 loops, as the vectorizer runs out of
    // registers with all 99 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph1_p1, ph3_p1, ph3_p2, ph5_p1, ph5_p2, ph7_p1, ph7_p2, ph7_p7, ph9_p1, ph9_p2, ph9_p7, ph9_p8, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h1_p1 = ph1_p1[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h3_p2 = ph3_p2[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h7_p7 = ph7_p7[k];
        const auto h9_p1 = ph9_p1[k];
        const auto h9_p2 = ph9_p2[k];
        const auto h9_p7 = ph9_p7[k];
        const auto h9_p8 = ph9_p8[k];

        pc_86[k] = e_0 * (-fs_1323_64 * h3_p1 - fs_3969_2 * r_2 * h1_p1) + e_1 * (-fs_375_8 * h5_p1 + fs_147_16 * r_2 * h3_p1 + fs_729_2 * r_4 * h1_p1) + e_2 * (-fs_92575_163592 * h7_p1 - fs_3675_1144 * h7_p7 + fs_750_169 * r_2 * h5_p1 - fs_147_484 * r_4 * h3_p1 - fs_8 * r_6 * h1_p1) + e_3 * (-fs_1960_5909761 * h9_p1 + fs_3920_41327 * h9_p7 + fs_92575_11819522 * r_2 * h7_p1 + fs_3675_82654 * r_2 * h7_p7 - fs_10_507 * r_4 * h5_p1 + fs_49_61347 * r_6 * h3_p1 + fs_2_121 * r_8 * h1_p1) + fs_99225_128 * e_4 * h1_p1;

        pc_87[k] = -fs_33075_64 * e_0 * h3_p2 + e_1 * (-fs_525_16 * h5_p2 + fs_3675_16 * r_2 * h3_p2) + e_2 * (-fs_23625_163592 * h7_p2 + fs_525_169 * r_2 * h5_p2 - fs_3675_484 * r_4 * h3_p2) + e_3 * (-fs_49_1074502 * h9_p2 + fs_196_2431 * h9_p8 + fs_23625_11819522 * r_2 * h7_p2 - fs_7_507 * r_4 * h5_p2 + fs_1225_61347 * r_6 * h3_p2);
    }

    // NOTE: the rows are formed in 59 loops, as the vectorizer runs out of
    // registers with all 99 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph1_m1, ph3_m1, ph5_m1, ph7_m7, ph7_m1, ph9_m9, ph9_m8, ph9_m7, ph9_m1, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h1_m1 = ph1_m1[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h7_m7 = ph7_m7[k];
        const auto h7_m1 = ph7_m1[k];
        const auto h9_m9 = ph9_m9[k];
        const auto h9_m8 = ph9_m8[k];
        const auto h9_m7 = ph9_m7[k];
        const auto h9_m1 = ph9_m1[k];

        pc_88[k] = e_0 * (fs_6615_16 * h3_m1 - fs_19845_8 * r_2 * h1_m1) + e_1 * (fs_75_8 * h5_m1 - fs_735_4 * r_2 * h3_m1 + fs_3645_8 * r_4 * h1_m1) + e_2 * (fs_875_40898 * h7_m1 - fs_150_169 * r_2 * h5_m1 + fs_735_121 * r_4 * h3_m1 - fs_10 * r_6 * h1_m1) + e_3 * (fs_441_2431 * h9_m9 + fs_49_11819522 * h9_m1 - fs_1750_5909761 * r_2 * h7_m1 + fs_2_507 * r_4 * h5_m1 - fs_980_61347 * r_6 * h3_m1 + fs_5_242 * r_8 * h1_m1) + fs_496125_512 * e_4 * h1_m1;

        pc_89[k] = fs_245_2431 * e_3 * h9_m8;

        pc_90[k] = e_0 * (fs_1323_4 * h3_m1 - fs_441_8 * r_2 * h1_m1) + e_1 * (fs_375_8 * h5_m1 - fs_147 * r_2 * h3_m1 + fs_81_8 * r_4 * h1_m1) + e_2 * (fs_1225_858 * h7_m7 + fs_68600_184041 * h7_m1 - fs_750_169 * r_2 * h5_m1 + fs_588_121 * r_4 * h3_m1 - fs_2_9 * r_6 * h1_m1) + e_3 * (fs_2205_41327 * h9_m7 + fs_2205_11819522 * h9_m1 - fs_2450_123981 * r_2 * h7_m7 - fs_274400_53187849 * r_2 * h7_m1 + fs_10_507 * r_4 * h5_m1 - fs_784_61347 * r_6 * h3_m1 + fs_1_2178 * r_8 * h1_m1) + fs_11025_512 * e_4 * h1_m1;
    }

    // NOTE: the rows are formed in 59 loops, as the vectorizer runs out of
    // registers with all 99 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph3_m3, ph3_m2, ph5_m5, ph5_m3, ph5_m2, ph7_m6, ph7_m5, ph7_m3, ph7_m2, ph9_m6, ph9_m5, ph9_m3, ph9_m2, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h3_m3 = ph3_m3[k];
        const auto h3_m2 = ph3_m2[k];
        const auto h5_m5 = ph5_m5[k];
        const auto h5_m3 = ph5_m3[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h7_m6 = ph7_m6[k];
        const auto h7_m5 = ph7_m5[k];
        const auto h7_m3 = ph7_m3[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h9_m6 = ph9_m6[k];
        const auto h9_m5 = ph9_m5[k];
        const auto h9_m3 = ph9_m3[k];
        const auto h9_m2 = ph9_m2[k];

        pc_91[k] = -fs_2205_16 * e_0 * h3_m2 + e_1 * (-fs_875_16 * h5_m2 + fs_245_4 * r_2 * h3_m2) + e_2 * (fs_350_143 * h7_m6 - fs_17150_20449 * h7_m2 + fs_875_169 * r_2 * h5_m2 - fs_245_121 * r_4 * h3_m2) + e_3 * (fs_2205_82654 * h9_m6 - fs_735_1074502 * h9_m2 - fs_1400_41327 * r_2 * h7_m6 + fs_68600_5909761 * r_2 * h7_m2 - fs_35_1521 * r_4 * h5_m2 + fs_980_184041 * r_6 * h3_m2);

        pc_92[k] = fs_945_32 * e_0 * h3_m3 + e_1 * (fs_75_8 * h5_m5 + fs_375_8 * h5_m3 - fs_105_8 * r_2 * h3_m3) + e_2 * (fs_4900_1859 * h7_m5 + f_175_143 * h7_m3 - fs_150_169 * r_2 * h5_m5 - fs_750_169 * r_2 * h5_m3 + fs_105_242 * r_4 * h3_m3) + e_3 * (fs_1029_82654 * h9_m5 + fs_2205_1074502 * h9_m3 - fs_19600_537251 * r_2 * h7_m5 - f_350_2431 * r_2 * h7_m3 + fs_2_507 * r_4 * h5_m5 + fs_10_507 * r_4 * h5_m3 - fs_70_61347 * r_6 * h3_m3);
    }

    // NOTE: the rows are formed in 59 loops, as the vectorizer runs out of
    // registers with all 99 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph3_p3, ph5_p3, ph5_p4, ph5_p5, ph7_p3, ph7_p4, ph7_p5, ph9_p3, ph9_p4, ph9_p5, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h3_p3 = ph3_p3[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h5_p4 = ph5_p4[k];
        const auto h5_p5 = ph5_p5[k];
        const auto h7_p3 = ph7_p3[k];
        const auto h7_p4 = ph7_p4[k];
        const auto h7_p5 = ph7_p5[k];
        const auto h9_p3 = ph9_p3[k];
        const auto h9_p4 = ph9_p4[k];
        const auto h9_p5 = ph9_p5[k];

        pc_93[k] = f_15_2 * e_1 * h5_p4 + e_2 * (fs_24500_5577 * h7_p4 - f_30_13 * r_2 * h5_p4) + e_3 * (fs_441_41327 * h9_p4 - fs_98000_1611753 * r_2 * h7_p4 + f_2_13 * r_4 * h5_p4);

        pc_94[k] = -fs_945_32 * e_0 * h3_p3 + e_1 * (-fs_375_8 * h5_p3 + fs_75_8 * h5_p5 + fs_105_8 * r_2 * h3_p3) + e_2 * (-f_175_143 * h7_p3 + fs_4900_1859 * h7_p5 + fs_750_169 * r_2 * h5_p3 - fs_150_169 * r_2 * h5_p5 - fs_105_242 * r_4 * h3_p3) + e_3 * (-fs_2205_1074502 * h9_p3 + fs_1029_82654 * h9_p5 + f_350_2431 * r_2 * h7_p3 - fs_19600_537251 * r_2 * h7_p5 - fs_10_507 * r_4 * h5_p3 + fs_2_507 * r_4 * h5_p5 + fs_70_61347 * r_6 * h3_p3);
    }

    // NOTE: the rows are formed in 59 loops, as the vectorizer runs out of
    // registers with all 99 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph3_p2, ph5_p2, ph7_p2, ph7_p6, ph9_p2, ph9_p6, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h3_p2 = ph3_p2[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h7_p6 = ph7_p6[k];
        const auto h9_p2 = ph9_p2[k];
        const auto h9_p6 = ph9_p6[k];

        pc_95[k] = fs_2205_16 * e_0 * h3_p2 + e_1 * (fs_875_16 * h5_p2 - fs_245_4 * r_2 * h3_p2) + e_2 * (fs_17150_20449 * h7_p2 + fs_350_143 * h7_p6 - fs_875_169 * r_2 * h5_p2 + fs_245_121 * r_4 * h3_p2) + e_3 * (fs_735_1074502 * h9_p2 + fs_2205_82654 * h9_p6 - fs_68600_5909761 * r_2 * h7_p2 - fs_1400_41327 * r_2 * h7_p6 + fs_35_1521 * r_4 * h5_p2 - fs_980_184041 * r_6 * h3_p2);
    }

    // NOTE: the rows are formed in 59 loops, as the vectorizer runs out of
    // registers with all 99 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph1_p1, ph3_p1, ph5_p1, ph7_p1, ph7_p7, ph9_p1, ph9_p7, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h1_p1 = ph1_p1[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h7_p7 = ph7_p7[k];
        const auto h9_p1 = ph9_p1[k];
        const auto h9_p7 = ph9_p7[k];

        pc_96[k] = e_0 * (-fs_1323_4 * h3_p1 + fs_441_8 * r_2 * h1_p1) + e_1 * (-fs_375_8 * h5_p1 + fs_147 * r_2 * h3_p1 - fs_81_8 * r_4 * h1_p1) + e_2 * (-fs_68600_184041 * h7_p1 + fs_1225_858 * h7_p7 + fs_750_169 * r_2 * h5_p1 - fs_588_121 * r_4 * h3_p1 + fs_2_9 * r_6 * h1_p1) + e_3 * (-fs_2205_11819522 * h9_p1 + fs_2205_41327 * h9_p7 + fs_274400_53187849 * r_2 * h7_p1 - fs_2450_123981 * r_2 * h7_p7 - fs_10_507 * r_4 * h5_p1 + fs_784_61347 * r_6 * h3_p1 - fs_1_2178 * r_8 * h1_p1) - fs_11025_512 * e_4 * h1_p1;
    }

    // NOTE: the rows are formed in 59 loops, as the vectorizer runs out of
    // registers with all 99 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph1_0, ph1_p1, ph3_0, ph3_p1, ph5_0, ph5_p1, ph7_0, ph7_p1, ph9_0, ph9_p1, ph9_p8, ph9_p9, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h1_0 = ph1_0[k];
        const auto h1_p1 = ph1_p1[k];
        const auto h3_0 = ph3_0[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h5_0 = ph5_0[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h7_0 = ph7_0[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h9_0 = ph9_0[k];
        const auto h9_p1 = ph9_p1[k];
        const auto h9_p8 = ph9_p8[k];
        const auto h9_p9 = ph9_p9[k];

        pc_97[k] = e_0 * (f_63_2 * h3_0 - f_63_2 * r_2 * h1_0) + e_1 * (f_15_2 * h5_0 - f_21 * r_2 * h3_0 + f_27_2 * r_4 * h1_0) + e_2 * (f_70_143 * h7_0 - f_30_13 * r_2 * h5_0 + f_42_11 * r_4 * h3_0 - f_2 * r_6 * h1_0) + e_3 * (f_21_2431 * h9_0 + fs_245_2431 * h9_p8 - f_140_2431 * r_2 * h7_0 + f_2_13 * r_4 * h5_0 - f_28_143 * r_6 * h3_0 + f_1_11 * r_8 * h1_0) + f_315_16 * e_4 * h1_0;

        pc_98[k] = e_0 * (fs_6615_16 * h3_p1 - fs_19845_8 * r_2 * h1_p1) + e_1 * (fs_75_8 * h5_p1 - fs_735_4 * r_2 * h3_p1 + fs_3645_8 * r_4 * h1_p1) + e_2 * (fs_875_40898 * h7_p1 - fs_150_169 * r_2 * h5_p1 + fs_735_121 * r_4 * h3_p1 - fs_10 * r_6 * h1_p1) + e_3 * (fs_49_11819522 * h9_p1 + fs_441_2431 * h9_p9 - fs_1750_5909761 * r_2 * h7_p1 + fs_2_507 * r_4 * h5_p1 - fs_980_61347 * r_6 * h3_p1 + fs_5_242 * r_8 * h1_p1) + fs_496125_512 * e_4 * h1_p1;
    }

    // NOTE: the values of a combination of angular components are stored as one
    // row of nvalues columns, with the component on bra side running slowest. The
    // rows which the symmetry relates to an already formed one are copied from it,
    // and the atom pairs beyond the reach of every pair of primitives are set to
    // zero.

    const size_t sources[99] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98};

    for (size_t m = 0; m < 99; m++)
    {
        auto *pv = values + m * nvalues;

        const auto *pc = values + sources[m] * nvalues;

        if (pv != pc) std::copy(pc, pc + nmax, pv);

        std::fill(pv + nmax, pv + nvalues, 0.0);
    }
}

}  // namespace simdovl
