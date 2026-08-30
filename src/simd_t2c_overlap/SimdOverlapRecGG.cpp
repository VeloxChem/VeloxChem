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



#include "SimdOverlapRecGG.hpp"

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
compute_gg_overlap(double                         *values,
                   const size_t                    nvalues,
                   const CBasisFunction           &bra,
                   const CBasisFunction           &ket,
                   const std::vector<CSimdMatrix> &harmonics,
                   const CSimdMatrix              &coordinates,
                   const double                    threshold) -> void
{
    if ((bra.get_angular_momentum() != 4) || (ket.get_angular_momentum() != 4))
    {
        errors::assertMsgCritical(
            false, std::string("SimdOverlapRecGG.compute_gg_overlap: Basis functions must be of angular momenta four and four"));
    }

    if (harmonics.size() < 8)
    {
        errors::assertMsgCritical(
            false, std::string("SimdOverlapRecGG.compute_gg_overlap: Harmonics must reach angular momentum eight"));
    }

    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("SimdOverlapRecGG.compute_gg_overlap: Number of values exceeds number of atom pairs"));
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
        std::fill(values, values + 81 * nvalues, 0.0);

        return;
    }

    // NOTE: the first five rows accumulate the contracted prefactors of the terms,
    // and the remaining 45 rows hold the integrals of the combinations of angular
    // components which are not related by symmetry.

    auto buffer = CSimdMatrix(50, nmax);

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

            const auto f_0 = fbase * fmu / fexp / fexp / fexp / fexp;

            const auto f_1 = fbase * fmu * fmu / fexp / fexp / fexp / fexp;

            const auto f_2 = fbase * fmu * fmu * fmu / fexp / fexp / fexp / fexp;

            const auto f_3 = fbase * fmu * fmu * fmu * fmu / fexp / fexp / fexp / fexp;

            const auto f_4 = fbase / fexp / fexp / fexp / fexp;

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

    const auto *ph2_m2 = harmonics[1].data(0);
    const auto *ph2_m1 = harmonics[1].data(1);
    const auto *ph2_0 = harmonics[1].data(2);
    const auto *ph2_p1 = harmonics[1].data(3);
    const auto *ph2_p2 = harmonics[1].data(4);
    const auto *ph4_m4 = harmonics[3].data(0);
    const auto *ph4_m3 = harmonics[3].data(1);
    const auto *ph4_m2 = harmonics[3].data(2);
    const auto *ph4_m1 = harmonics[3].data(3);
    const auto *ph4_0 = harmonics[3].data(4);
    const auto *ph4_p1 = harmonics[3].data(5);
    const auto *ph4_p2 = harmonics[3].data(6);
    const auto *ph4_p3 = harmonics[3].data(7);
    const auto *ph4_p4 = harmonics[3].data(8);
    const auto *ph6_m6 = harmonics[5].data(0);
    const auto *ph6_m5 = harmonics[5].data(1);
    const auto *ph6_m4 = harmonics[5].data(2);
    const auto *ph6_m3 = harmonics[5].data(3);
    const auto *ph6_m2 = harmonics[5].data(4);
    const auto *ph6_m1 = harmonics[5].data(5);
    const auto *ph6_0 = harmonics[5].data(6);
    const auto *ph6_p1 = harmonics[5].data(7);
    const auto *ph6_p2 = harmonics[5].data(8);
    const auto *ph6_p3 = harmonics[5].data(9);
    const auto *ph6_p4 = harmonics[5].data(10);
    const auto *ph6_p5 = harmonics[5].data(11);
    const auto *ph6_p6 = harmonics[5].data(12);
    const auto *ph8_m8 = harmonics[7].data(0);
    const auto *ph8_m7 = harmonics[7].data(1);
    const auto *ph8_m6 = harmonics[7].data(2);
    const auto *ph8_m5 = harmonics[7].data(3);
    const auto *ph8_m4 = harmonics[7].data(4);
    const auto *ph8_m3 = harmonics[7].data(5);
    const auto *ph8_m2 = harmonics[7].data(6);
    const auto *ph8_m1 = harmonics[7].data(7);
    const auto *ph8_0 = harmonics[7].data(8);
    const auto *ph8_p1 = harmonics[7].data(9);
    const auto *ph8_p2 = harmonics[7].data(10);
    const auto *ph8_p3 = harmonics[7].data(11);
    const auto *ph8_p4 = harmonics[7].data(12);
    const auto *ph8_p5 = harmonics[7].data(13);
    const auto *ph8_p6 = harmonics[7].data(14);
    const auto *ph8_p7 = harmonics[7].data(15);
    const auto *ph8_p8 = harmonics[7].data(16);

    auto *pc_0 = buffer.data(5);
    auto *pc_1 = buffer.data(6);
    auto *pc_2 = buffer.data(7);
    auto *pc_3 = buffer.data(8);
    auto *pc_4 = buffer.data(9);
    auto *pc_5 = buffer.data(10);
    auto *pc_6 = buffer.data(11);
    auto *pc_7 = buffer.data(12);
    auto *pc_8 = buffer.data(13);
    auto *pc_9 = buffer.data(14);
    auto *pc_10 = buffer.data(15);
    auto *pc_11 = buffer.data(16);
    auto *pc_12 = buffer.data(17);
    auto *pc_13 = buffer.data(18);
    auto *pc_14 = buffer.data(19);
    auto *pc_15 = buffer.data(20);
    auto *pc_16 = buffer.data(21);
    auto *pc_17 = buffer.data(22);
    auto *pc_18 = buffer.data(23);
    auto *pc_19 = buffer.data(24);
    auto *pc_20 = buffer.data(25);
    auto *pc_21 = buffer.data(26);
    auto *pc_22 = buffer.data(27);
    auto *pc_23 = buffer.data(28);
    auto *pc_24 = buffer.data(29);
    auto *pc_25 = buffer.data(30);
    auto *pc_26 = buffer.data(31);
    auto *pc_27 = buffer.data(32);
    auto *pc_28 = buffer.data(33);
    auto *pc_29 = buffer.data(34);
    auto *pc_30 = buffer.data(35);
    auto *pc_31 = buffer.data(36);
    auto *pc_32 = buffer.data(37);
    auto *pc_33 = buffer.data(38);
    auto *pc_34 = buffer.data(39);
    auto *pc_35 = buffer.data(40);
    auto *pc_36 = buffer.data(41);
    auto *pc_37 = buffer.data(42);
    auto *pc_38 = buffer.data(43);
    auto *pc_39 = buffer.data(44);
    auto *pc_40 = buffer.data(45);
    auto *pc_41 = buffer.data(46);
    auto *pc_42 = buffer.data(47);
    auto *pc_43 = buffer.data(48);
    auto *pc_44 = buffer.data(49);

    // NOTE: the factors of the terms depend on the angular momenta alone, so they
    // are formed once for the whole matrix instead of once for every atom pair.

    const auto f_35_2 = 17.5;
    const auto f_9_2 = 4.5;
    const auto f_15 = 15.0;
    const auto f_21_2 = 10.5;
    const auto f_10_33 = 10.0 / 33.0;
    const auto f_18_11 = 18.0 / 11.0;
    const auto f_10_3 = 10.0 / 3.0;
    const auto f_2 = 2.0;
    const auto f_7_1287 = 7.0 / 1287.0;
    const auto fs_245_1287 = std::sqrt(245.0 / 1287.0);
    const auto f_4_99 = 4.0 / 99.0;
    const auto f_18_143 = 18.0 / 143.0;
    const auto f_20_99 = 20.0 / 99.0;
    const auto f_1_9 = 1.0 / 9.0;
    const auto f_105_16 = 6.5625;
    const auto fs_3675_32 = std::sqrt(114.84375);
    const auto fs_405_16 = std::sqrt(25.3125);
    const auto fs_675_8 = std::sqrt(84.375);
    const auto fs_175_726 = std::sqrt(175.0 / 726.0);
    const auto fs_405_121 = std::sqrt(405.0 / 121.0);
    const auto fs_25_6 = std::sqrt(25.0 / 6.0);
    const auto fs_49_368082 = std::sqrt(49.0 / 368082.0);
    const auto fs_245_2574 = std::sqrt(245.0 / 2574.0);
    const auto fs_14_3267 = std::sqrt(14.0 / 3267.0);
    const auto fs_405_20449 = std::sqrt(405.0 / 20449.0);
    const auto fs_50_3267 = std::sqrt(50.0 / 3267.0);
    const auto fs_525_16 = std::sqrt(32.8125);
    const auto fs_3645_112 = std::sqrt(3645.0 / 112.0);
    const auto fs_675_28 = std::sqrt(675.0 / 28.0);
    const auto fs_250_363 = std::sqrt(250.0 / 363.0);
    const auto fs_50_33 = std::sqrt(50.0 / 33.0);
    const auto fs_3645_847 = std::sqrt(3645.0 / 847.0);
    const auto fs_25_21 = std::sqrt(25.0 / 21.0);
    const auto fs_245_368082 = std::sqrt(245.0 / 368082.0);
    const auto fs_343_7722 = std::sqrt(343.0 / 7722.0);
    const auto fs_40_3267 = std::sqrt(40.0 / 3267.0);
    const auto fs_8_297 = std::sqrt(8.0 / 297.0);
    const auto fs_3645_143143 = std::sqrt(3645.0 / 143143.0);
    const auto fs_100_22869 = std::sqrt(100.0 / 22869.0);
    const auto fs_500_363 = std::sqrt(500.0 / 363.0);
    const auto fs_25_11 = std::sqrt(25.0 / 11.0);
    const auto fs_245_100386 = std::sqrt(245.0 / 100386.0);
    const auto fs_49_2574 = std::sqrt(49.0 / 2574.0);
    const auto fs_80_3267 = std::sqrt(80.0 / 3267.0);
    const auto fs_4_99 = std::sqrt(4.0 / 99.0);
    const auto fs_500_121 = std::sqrt(500.0 / 121.0);
    const auto fs_245_16731 = std::sqrt(245.0 / 16731.0);
    const auto fs_80_1089 = std::sqrt(80.0 / 1089.0);
    const auto f_35_8 = 4.375;
    const auto f_27_4 = 6.75;
    const auto f_15_4 = 3.75;
    const auto f_85_66 = 85.0 / 66.0;
    const auto fs_175_66 = std::sqrt(175.0 / 66.0);
    const auto f_27_11 = 27.0 / 11.0;
    const auto f_5_6 = 5.0 / 6.0;
    const auto f_56_1287 = 56.0 / 1287.0;
    const auto fs_392_3861 = std::sqrt(392.0 / 3861.0);
    const auto f_17_99 = 17.0 / 99.0;
    const auto fs_14_297 = std::sqrt(14.0 / 297.0);
    const auto f_27_143 = 27.0 / 143.0;
    const auto f_5_99 = 5.0 / 99.0;
    const auto fs_13125_128 = std::sqrt(102.5390625);
    const auto fs_405_112 = std::sqrt(405.0 / 112.0);
    const auto fs_16875_224 = std::sqrt(16875.0 / 224.0);
    const auto fs_4225_2904 = std::sqrt(4225.0 / 2904.0);
    const auto fs_25_44 = std::sqrt(25.0 / 44.0);
    const auto fs_405_847 = std::sqrt(405.0 / 847.0);
    const auto fs_625_168 = std::sqrt(625.0 / 168.0);
    const auto fs_686_184041 = std::sqrt(686.0 / 184041.0);
    const auto fs_98_1287 = std::sqrt(98.0 / 1287.0);
    const auto fs_169_6534 = std::sqrt(169.0 / 6534.0);
    const auto fs_1_99 = std::sqrt(1.0 / 99.0);
    const auto fs_405_143143 = std::sqrt(405.0 / 143143.0);
    const auto fs_625_45738 = std::sqrt(625.0 / 45738.0);
    const auto fs_4725_64 = std::sqrt(73.828125);
    const auto fs_6075_112 = std::sqrt(6075.0 / 112.0);
    const auto fs_375_242 = std::sqrt(375.0 / 242.0);
    const auto f_5_22 = 5.0 / 22.0;
    const auto fs_75_28 = std::sqrt(75.0 / 28.0);
    const auto fs_1960_184041 = std::sqrt(1960.0 / 184041.0);
    const auto fs_784_16731 = std::sqrt(784.0 / 16731.0);
    const auto fs_10_363 = std::sqrt(10.0 / 363.0);
    const auto f_1_33 = 1.0 / 33.0;
    const auto fs_25_2541 = std::sqrt(25.0 / 2541.0);
    const auto fs_625_363 = std::sqrt(625.0 / 363.0);
    const auto fs_2450_50193 = std::sqrt(2450.0 / 50193.0);
    const auto fs_100_3267 = std::sqrt(100.0 / 3267.0);
    const auto f_5 = 5.0;
    const auto f_99_28 = 99.0 / 28.0;
    const auto f_30_7 = 30.0 / 7.0;
    const auto f_5_3 = 5.0 / 3.0;
    const auto fs_175_121 = std::sqrt(175.0 / 121.0);
    const auto f_9_7 = 9.0 / 7.0;
    const auto f_20_21 = 20.0 / 21.0;
    const auto f_196_1287 = 196.0 / 1287.0;
    const auto fs_1372_16731 = std::sqrt(1372.0 / 16731.0);
    const auto f_2_9 = 2.0 / 9.0;
    const auto fs_28_1089 = std::sqrt(28.0 / 1089.0);
    const auto f_9_91 = 9.0 / 91.0;
    const auto f_40_693 = 40.0 / 693.0;
    const auto fs_6075_128 = std::sqrt(47.4609375);
    const auto fs_3645_196 = std::sqrt(3645.0 / 196.0);
    const auto fs_54675_1568 = std::sqrt(54675.0 / 1568.0);
    const auto fs_525_968 = std::sqrt(525.0 / 968.0);
    const auto fs_875_1452 = std::sqrt(875.0 / 1452.0);
    const auto fs_14580_5929 = std::sqrt(14580.0 / 5929.0);
    const auto fs_675_392 = std::sqrt(675.0 / 392.0);
    const auto fs_4802_184041 = std::sqrt(4802.0 / 184041.0);
    const auto fs_3430_50193 = std::sqrt(3430.0 / 50193.0);
    const auto fs_7_726 = std::sqrt(7.0 / 726.0);
    const auto fs_35_3267 = std::sqrt(35.0 / 3267.0);
    const auto fs_14580_1002001 = std::sqrt(14580.0 / 1002001.0);
    const auto fs_75_11858 = std::sqrt(75.0 / 11858.0);
    const auto fs_3375_16 = std::sqrt(210.9375);
    const auto fs_30375_196 = std::sqrt(30375.0 / 196.0);
    const auto fs_375_49 = std::sqrt(375.0 / 49.0);
    const auto fs_17150_184041 = std::sqrt(17150.0 / 184041.0);
    const auto fs_500_17787 = std::sqrt(500.0 / 17787.0);
    const auto f_85_8 = 10.625;
    const auto fs_1875_16 = std::sqrt(117.1875);
    const auto f_81_28 = 81.0 / 28.0;
    const auto f_255_28 = 255.0 / 28.0;
    const auto fs_16875_196 = std::sqrt(16875.0 / 196.0);
    const auto f_5_66 = 5.0 / 66.0;
    const auto fs_875_726 = std::sqrt(875.0 / 726.0);
    const auto f_81_77 = 81.0 / 77.0;
    const auto f_85_42 = 85.0 / 42.0;
    const auto fs_625_147 = std::sqrt(625.0 / 147.0);
    const auto f_392_1287 = 392.0 / 1287.0;
    const auto fs_13720_184041 = std::sqrt(13720.0 / 184041.0);
    const auto f_1_99 = 1.0 / 99.0;
    const auto fs_70_3267 = std::sqrt(70.0 / 3267.0);
    const auto f_81_1001 = 81.0 / 1001.0;
    const auto f_85_693 = 85.0 / 693.0;
    const auto fs_2500_160083 = std::sqrt(2500.0 / 160083.0);
    const auto fs_375_32 = std::sqrt(11.71875);
    const auto fs_3375_392 = std::sqrt(3375.0 / 392.0);
    const auto fs_125_294 = std::sqrt(125.0 / 294.0);
    const auto fs_24010_184041 = std::sqrt(24010.0 / 184041.0);
    const auto fs_250_160083 = std::sqrt(250.0 / 160083.0);
    const auto f_25_2 = 12.5;
    const auto f_81_14 = 81.0 / 14.0;
    const auto f_75_7 = 75.0 / 7.0;
    const auto f_50_33 = 50.0 / 33.0;
    const auto f_162_77 = 162.0 / 77.0;
    const auto f_50_21 = 50.0 / 21.0;
    const auto f_490_1287 = 490.0 / 1287.0;
    const auto f_162_1001 = 162.0 / 1001.0;
    const auto f_100_693 = 100.0 / 693.0;

    // NOTE: the rows are formed in 23 loops, as the vectorizer runs out of
    // registers with all 45 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_0, ph2_p1, ph4_0, ph4_p1, ph6_0, ph6_p1, ph8_0, ph8_p1, ph8_p7, ph8_p8, ab_2, pc_0, pc_1 : simd::cache_line_size())
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

        const auto h2_0 = ph2_0[k];
        const auto h2_p1 = ph2_p1[k];
        const auto h4_0 = ph4_0[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h6_0 = ph6_0[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h8_0 = ph8_0[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h8_p7 = ph8_p7[k];
        const auto h8_p8 = ph8_p8[k];

        pc_0[k] = e_0 * (f_35_2 * h2_0 - f_35_2 * r_2) + e_1 * (f_9_2 * h4_0 - f_15 * r_2 * h2_0 + f_21_2 * r_4) + e_2 * (f_10_33 * h6_0 - f_18_11 * r_2 * h4_0 + f_10_3 * r_4 * h2_0 - f_2 * r_6) + e_3 * (f_7_1287 * h8_0 - fs_245_1287 * h8_p8 - f_4_99 * r_2 * h6_0 + f_18_143 * r_4 * h4_0 - f_20_99 * r_6 * h2_0 + f_1_9 * r_8) + f_105_16 * e_4;

        pc_1[k] = -fs_3675_32 * e_0 * h2_p1 + e_1 * (-fs_405_16 * h4_p1 + fs_675_8 * r_2 * h2_p1) + e_2 * (-fs_175_726 * h6_p1 + fs_405_121 * r_2 * h4_p1 - fs_25_6 * r_4 * h2_p1) + e_3 * (-fs_49_368082 * h8_p1 - fs_245_2574 * h8_p7 + fs_14_3267 * r_2 * h6_p1 - fs_405_20449 * r_4 * h4_p1 + fs_50_3267 * r_6 * h2_p1);
    }

    // NOTE: the rows are formed in 23 loops, as the vectorizer runs out of
    // registers with all 45 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph2_p2, ph4_p2, ph4_p3, ph6_p2, ph6_p3, ph6_p5, ph6_p6, ph8_p2, ph8_p3, ph8_p5, ph8_p6, ab_2, pc_2, pc_3 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h2_p2 = ph2_p2[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h4_p3 = ph4_p3[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h6_p5 = ph6_p5[k];
        const auto h6_p6 = ph6_p6[k];
        const auto h8_p2 = ph8_p2[k];
        const auto h8_p3 = ph8_p3[k];
        const auto h8_p5 = ph8_p5[k];
        const auto h8_p6 = ph8_p6[k];

        pc_2[k] = fs_525_16 * e_0 * h2_p2 + e_1 * (fs_3645_112 * h4_p2 - fs_675_28 * r_2 * h2_p2) + e_2 * (fs_250_363 * h6_p2 - fs_50_33 * h6_p6 - fs_3645_847 * r_2 * h4_p2 + fs_25_21 * r_4 * h2_p2) + e_3 * (fs_245_368082 * h8_p2 - fs_343_7722 * h8_p6 - fs_40_3267 * r_2 * h6_p2 + fs_8_297 * r_2 * h6_p6 + fs_3645_143143 * r_4 * h4_p2 - fs_100_22869 * r_6 * h2_p2);

        pc_3[k] = -fs_405_16 * e_1 * h4_p3 + e_2 * (-fs_500_363 * h6_p3 - fs_25_11 * h6_p5 + fs_405_121 * r_2 * h4_p3) + e_3 * (-fs_245_100386 * h8_p3 - fs_49_2574 * h8_p5 + fs_80_3267 * r_2 * h6_p3 + fs_4_99 * r_2 * h6_p5 - fs_405_20449 * r_4 * h4_p3);
    }

    // NOTE: the rows are formed in 23 loops, as the vectorizer runs out of
    // registers with all 45 of them in one.

#pragma omp simd aligned(pe_1, pe_2, pe_3, ph4_m4, ph4_m3, ph6_m5, ph6_m4, ph6_m3, ph8_m5, ph8_m4, ph8_m3, ab_2, pc_4, pc_5 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

        const auto h4_m4 = ph4_m4[k];
        const auto h4_m3 = ph4_m3[k];
        const auto h6_m5 = ph6_m5[k];
        const auto h6_m4 = ph6_m4[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h8_m5 = ph8_m5[k];
        const auto h8_m4 = ph8_m4[k];
        const auto h8_m3 = ph8_m3[k];

        pc_4[k] = f_9_2 * e_1 * h4_m4 + e_2 * (fs_500_121 * h6_m4 - f_18_11 * r_2 * h4_m4) + e_3 * (fs_245_16731 * h8_m4 - fs_80_1089 * r_2 * h6_m4 + f_18_143 * r_4 * h4_m4);

        pc_5[k] = -fs_405_16 * e_1 * h4_m3 + e_2 * (fs_25_11 * h6_m5 - fs_500_363 * h6_m3 + fs_405_121 * r_2 * h4_m3) + e_3 * (fs_49_2574 * h8_m5 - fs_245_100386 * h8_m3 - fs_4_99 * r_2 * h6_m5 + fs_80_3267 * r_2 * h6_m3 - fs_405_20449 * r_4 * h4_m3);
    }

    // NOTE: the rows are formed in 23 loops, as the vectorizer runs out of
    // registers with all 45 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph2_m2, ph2_m1, ph4_m2, ph4_m1, ph6_m6, ph6_m2, ph6_m1, ph8_m8, ph8_m7, ph8_m6, ph8_m2, ph8_m1, ab_2, pc_6, pc_7, pc_8 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h2_m2 = ph2_m2[k];
        const auto h2_m1 = ph2_m1[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h6_m6 = ph6_m6[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h8_m8 = ph8_m8[k];
        const auto h8_m7 = ph8_m7[k];
        const auto h8_m6 = ph8_m6[k];
        const auto h8_m2 = ph8_m2[k];
        const auto h8_m1 = ph8_m1[k];

        pc_6[k] = fs_525_16 * e_0 * h2_m2 + e_1 * (fs_3645_112 * h4_m2 - fs_675_28 * r_2 * h2_m2) + e_2 * (fs_50_33 * h6_m6 + fs_250_363 * h6_m2 - fs_3645_847 * r_2 * h4_m2 + fs_25_21 * r_4 * h2_m2) + e_3 * (fs_343_7722 * h8_m6 + fs_245_368082 * h8_m2 - fs_8_297 * r_2 * h6_m6 - fs_40_3267 * r_2 * h6_m2 + fs_3645_143143 * r_4 * h4_m2 - fs_100_22869 * r_6 * h2_m2);

        pc_7[k] = -fs_3675_32 * e_0 * h2_m1 + e_1 * (-fs_405_16 * h4_m1 + fs_675_8 * r_2 * h2_m1) + e_2 * (-fs_175_726 * h6_m1 + fs_405_121 * r_2 * h4_m1 - fs_25_6 * r_4 * h2_m1) + e_3 * (fs_245_2574 * h8_m7 - fs_49_368082 * h8_m1 + fs_14_3267 * r_2 * h6_m1 - fs_405_20449 * r_4 * h4_m1 + fs_50_3267 * r_6 * h2_m1);

        pc_8[k] = fs_245_1287 * e_3 * h8_m8;
    }

    // NOTE: the rows are formed in 23 loops, as the vectorizer runs out of
    // registers with all 45 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_0, ph2_p1, ph4_0, ph4_p1, ph6_0, ph6_p1, ph6_p5, ph6_p6, ph8_0, ph8_p1, ph8_p5, ph8_p6, ab_2, pc_9, pc_10 : simd::cache_line_size())
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

        const auto h2_0 = ph2_0[k];
        const auto h2_p1 = ph2_p1[k];
        const auto h4_0 = ph4_0[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h6_0 = ph6_0[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h6_p5 = ph6_p5[k];
        const auto h6_p6 = ph6_p6[k];
        const auto h8_0 = ph8_0[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h8_p5 = ph8_p5[k];
        const auto h8_p6 = ph8_p6[k];

        pc_9[k] = e_0 * (f_35_8 * h2_0 - f_35_2 * r_2) + e_1 * (-f_27_4 * h4_0 - f_15_4 * r_2 * h2_0 + f_21_2 * r_4) + e_2 * (-f_85_66 * h6_0 + fs_175_66 * h6_p6 + f_27_11 * r_2 * h4_0 + f_5_6 * r_4 * h2_0 - f_2 * r_6) + e_3 * (-f_56_1287 * h8_0 - fs_392_3861 * h8_p6 + f_17_99 * r_2 * h6_0 - fs_14_297 * r_2 * h6_p6 - f_27_143 * r_4 * h4_0 - f_5_99 * r_6 * h2_0 + f_1_9 * r_8) + f_105_16 * e_4;

        pc_10[k] = -fs_13125_128 * e_0 * h2_p1 + e_1 * (fs_405_112 * h4_p1 + fs_16875_224 * r_2 * h2_p1) + e_2 * (fs_4225_2904 * h6_p1 + fs_25_44 * h6_p5 - fs_405_847 * r_2 * h4_p1 - fs_625_168 * r_4 * h2_p1) + e_3 * (fs_686_184041 * h8_p1 - fs_98_1287 * h8_p5 - fs_169_6534 * r_2 * h6_p1 - fs_1_99 * r_2 * h6_p5 + fs_405_143143 * r_4 * h4_p1 + fs_625_45738 * r_6 * h2_p1);
    }

    // NOTE: the rows are formed in 23 loops, as the vectorizer runs out of
    // registers with all 45 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph2_p2, ph4_m3, ph4_p2, ph4_p4, ph6_m3, ph6_p2, ph6_p4, ph8_m3, ph8_p2, ph8_p4, ab_2, pc_11, pc_12 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h2_p2 = ph2_p2[k];
        const auto h4_m3 = ph4_m3[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h4_p4 = ph4_p4[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h6_p4 = ph6_p4[k];
        const auto h8_m3 = ph8_m3[k];
        const auto h8_p2 = ph8_p2[k];
        const auto h8_p4 = ph8_p4[k];

        pc_11[k] = fs_4725_64 * e_0 * h2_p2 + e_1 * (fs_405_112 * h4_p2 + fs_405_16 * h4_p4 - fs_6075_112 * r_2 * h2_p2) + e_2 * (-fs_375_242 * h6_p2 - f_5_22 * h6_p4 - fs_405_847 * r_2 * h4_p2 - fs_405_121 * r_2 * h4_p4 + fs_75_28 * r_4 * h2_p2) + e_3 * (-fs_1960_184041 * h8_p2 - fs_784_16731 * h8_p4 + fs_10_363 * r_2 * h6_p2 + f_1_33 * r_2 * h6_p4 + fs_405_143143 * r_4 * h4_p2 + fs_405_20449 * r_4 * h4_p4 - fs_25_2541 * r_6 * h2_p2);

        pc_12[k] = -f_27_4 * e_1 * h4_m3 + e_2 * (fs_625_363 * h6_m3 + f_27_11 * r_2 * h4_m3) + e_3 * (fs_2450_50193 * h8_m3 - fs_100_3267 * r_2 * h6_m3 - f_27_143 * r_4 * h4_m3);
    }

    // NOTE: the rows are formed in 23 loops, as the vectorizer runs out of
    // registers with all 45 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph2_m2, ph2_m1, ph4_m4, ph4_m2, ph4_m1, ph6_m5, ph6_m4, ph6_m2, ph6_m1, ph8_m5, ph8_m4, ph8_m2, ph8_m1, ab_2, pc_13, pc_14 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h2_m2 = ph2_m2[k];
        const auto h2_m1 = ph2_m1[k];
        const auto h4_m4 = ph4_m4[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h6_m5 = ph6_m5[k];
        const auto h6_m4 = ph6_m4[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h8_m5 = ph8_m5[k];
        const auto h8_m4 = ph8_m4[k];
        const auto h8_m2 = ph8_m2[k];
        const auto h8_m1 = ph8_m1[k];

        pc_13[k] = fs_4725_64 * e_0 * h2_m2 + e_1 * (-fs_405_16 * h4_m4 + fs_405_112 * h4_m2 - fs_6075_112 * r_2 * h2_m2) + e_2 * (f_5_22 * h6_m4 - fs_375_242 * h6_m2 + fs_405_121 * r_2 * h4_m4 - fs_405_847 * r_2 * h4_m2 + fs_75_28 * r_4 * h2_m2) + e_3 * (fs_784_16731 * h8_m4 - fs_1960_184041 * h8_m2 - f_1_33 * r_2 * h6_m4 + fs_10_363 * r_2 * h6_m2 - fs_405_20449 * r_4 * h4_m4 + fs_405_143143 * r_4 * h4_m2 - fs_25_2541 * r_6 * h2_m2);

        pc_14[k] = -fs_13125_128 * e_0 * h2_m1 + e_1 * (fs_405_112 * h4_m1 + fs_16875_224 * r_2 * h2_m1) + e_2 * (-fs_25_44 * h6_m5 + fs_4225_2904 * h6_m1 - fs_405_847 * r_2 * h4_m1 - fs_625_168 * r_4 * h2_m1) + e_3 * (fs_98_1287 * h8_m5 + fs_686_184041 * h8_m1 + fs_1_99 * r_2 * h6_m5 - fs_169_6534 * r_2 * h6_m1 + fs_405_143143 * r_4 * h4_m1 + fs_625_45738 * r_6 * h2_m1);
    }

    // NOTE: the rows are formed in 23 loops, as the vectorizer runs out of
    // registers with all 45 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph2_m1, ph4_m1, ph6_m6, ph6_m1, ph8_m7, ph8_m6, ph8_m1, ab_2, pc_15, pc_16 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h2_m1 = ph2_m1[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h6_m6 = ph6_m6[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h8_m7 = ph8_m7[k];
        const auto h8_m6 = ph8_m6[k];
        const auto h8_m1 = ph8_m1[k];

        pc_15[k] = -fs_175_66 * e_2 * h6_m6 + e_3 * (fs_392_3861 * h8_m6 + fs_14_297 * r_2 * h6_m6);

        pc_16[k] = fs_3675_32 * e_0 * h2_m1 + e_1 * (fs_405_16 * h4_m1 - fs_675_8 * r_2 * h2_m1) + e_2 * (fs_175_726 * h6_m1 - fs_405_121 * r_2 * h4_m1 + fs_25_6 * r_4 * h2_m1) + e_3 * (fs_245_2574 * h8_m7 + fs_49_368082 * h8_m1 - fs_14_3267 * r_2 * h6_m1 + fs_405_20449 * r_4 * h4_m1 - fs_50_3267 * r_6 * h2_m1);
    }

    // NOTE: the rows are formed in 23 loops, as the vectorizer runs out of
    // registers with all 45 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_0, ph4_0, ph4_p4, ph6_0, ph6_p4, ph8_0, ph8_p4, ab_2, pc_17 : simd::cache_line_size())
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

        const auto h2_0 = ph2_0[k];
        const auto h4_0 = ph4_0[k];
        const auto h4_p4 = ph4_p4[k];
        const auto h6_0 = ph6_0[k];
        const auto h6_p4 = ph6_p4[k];
        const auto h8_0 = ph8_0[k];
        const auto h8_p4 = ph8_p4[k];

        pc_17[k] = e_0 * (-f_5 * h2_0 - f_35_2 * r_2) + e_1 * (-f_99_28 * h4_0 - fs_3645_112 * h4_p4 + f_30_7 * r_2 * h2_0 + f_21_2 * r_4) + e_2 * (f_5_3 * h6_0 + fs_175_121 * h6_p4 + f_9_7 * r_2 * h4_0 + fs_3645_847 * r_2 * h4_p4 - f_20_21 * r_4 * h2_0 - f_2 * r_6) + e_3 * (f_196_1287 * h8_0 - fs_1372_16731 * h8_p4 - f_2_9 * r_2 * h6_0 - fs_28_1089 * r_2 * h6_p4 - f_9_91 * r_4 * h4_0 - fs_3645_143143 * r_4 * h4_p4 + f_40_693 * r_6 * h2_0 + f_1_9 * r_8) + f_105_16 * e_4;
    }

    // NOTE: the rows are formed in 23 loops, as the vectorizer runs out of
    // registers with all 45 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph2_m2, ph2_p1, ph4_m2, ph4_p1, ph4_p3, ph6_p1, ph6_p3, ph8_m2, ph8_p1, ph8_p3, ab_2, pc_18, pc_19 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h2_m2 = ph2_m2[k];
        const auto h2_p1 = ph2_p1[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h4_p3 = ph4_p3[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h8_m2 = ph8_m2[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h8_p3 = ph8_p3[k];

        pc_18[k] = -fs_6075_128 * e_0 * h2_p1 + e_1 * (fs_3645_196 * h4_p1 - fs_405_112 * h4_p3 + fs_54675_1568 * r_2 * h2_p1) + e_2 * (-fs_525_968 * h6_p1 + fs_875_1452 * h6_p3 - fs_14580_5929 * r_2 * h4_p1 + fs_405_847 * r_2 * h4_p3 - fs_675_392 * r_4 * h2_p1) + e_3 * (-fs_4802_184041 * h8_p1 - fs_3430_50193 * h8_p3 + fs_7_726 * r_2 * h6_p1 - fs_35_3267 * r_2 * h6_p3 + fs_14580_1002001 * r_4 * h4_p1 - fs_405_143143 * r_4 * h4_p3 + fs_75_11858 * r_6 * h2_p1);

        pc_19[k] = fs_3375_16 * e_0 * h2_m2 + e_1 * (-f_99_28 * h4_m2 - fs_30375_196 * r_2 * h2_m2) + e_2 * (f_9_7 * r_2 * h4_m2 + fs_375_49 * r_4 * h2_m2) + e_3 * (fs_17150_184041 * h8_m2 - f_9_91 * r_4 * h4_m2 - fs_500_17787 * r_6 * h2_m2);
    }

    // NOTE: the rows are formed in 23 loops, as the vectorizer runs out of
    // registers with all 45 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph2_m1, ph4_m4, ph4_m3, ph4_m1, ph6_m5, ph6_m4, ph6_m3, ph6_m1, ph8_m5, ph8_m4, ph8_m3, ph8_m1, ab_2, pc_20, pc_21, pc_22 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h2_m1 = ph2_m1[k];
        const auto h4_m4 = ph4_m4[k];
        const auto h4_m3 = ph4_m3[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h6_m5 = ph6_m5[k];
        const auto h6_m4 = ph6_m4[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h8_m5 = ph8_m5[k];
        const auto h8_m4 = ph8_m4[k];
        const auto h8_m3 = ph8_m3[k];
        const auto h8_m1 = ph8_m1[k];

        pc_20[k] = -fs_6075_128 * e_0 * h2_m1 + e_1 * (fs_405_112 * h4_m3 + fs_3645_196 * h4_m1 + fs_54675_1568 * r_2 * h2_m1) + e_2 * (-fs_875_1452 * h6_m3 - fs_525_968 * h6_m1 - fs_405_847 * r_2 * h4_m3 - fs_14580_5929 * r_2 * h4_m1 - fs_675_392 * r_4 * h2_m1) + e_3 * (fs_3430_50193 * h8_m3 - fs_4802_184041 * h8_m1 + fs_35_3267 * r_2 * h6_m3 + fs_7_726 * r_2 * h6_m1 + fs_405_143143 * r_4 * h4_m3 + fs_14580_1002001 * r_4 * h4_m1 + fs_75_11858 * r_6 * h2_m1);

        pc_21[k] = fs_3645_112 * e_1 * h4_m4 + e_2 * (-fs_175_121 * h6_m4 - fs_3645_847 * r_2 * h4_m4) + e_3 * (fs_1372_16731 * h8_m4 + fs_28_1089 * r_2 * h6_m4 + fs_3645_143143 * r_4 * h4_m4);

        pc_22[k] = fs_13125_128 * e_0 * h2_m1 + e_1 * (-fs_405_112 * h4_m1 - fs_16875_224 * r_2 * h2_m1) + e_2 * (-fs_25_44 * h6_m5 - fs_4225_2904 * h6_m1 + fs_405_847 * r_2 * h4_m1 + fs_625_168 * r_4 * h2_m1) + e_3 * (fs_98_1287 * h8_m5 - fs_686_184041 * h8_m1 + fs_1_99 * r_2 * h6_m5 + fs_169_6534 * r_2 * h6_m1 - fs_405_143143 * r_4 * h4_m1 - fs_625_45738 * r_6 * h2_m1);
    }

    // NOTE: the rows are formed in 23 loops, as the vectorizer runs out of
    // registers with all 45 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph2_m2, ph4_m2, ph6_m6, ph6_m2, ph8_m6, ph8_m2, ab_2, pc_23 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h2_m2 = ph2_m2[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h6_m6 = ph6_m6[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h8_m6 = ph8_m6[k];
        const auto h8_m2 = ph8_m2[k];

        pc_23[k] = -fs_525_16 * e_0 * h2_m2 + e_1 * (-fs_3645_112 * h4_m2 + fs_675_28 * r_2 * h2_m2) + e_2 * (fs_50_33 * h6_m6 - fs_250_363 * h6_m2 + fs_3645_847 * r_2 * h4_m2 - fs_25_21 * r_4 * h2_m2) + e_3 * (fs_343_7722 * h8_m6 - fs_245_368082 * h8_m2 - fs_8_297 * r_2 * h6_m6 + fs_40_3267 * r_2 * h6_m2 - fs_3645_143143 * r_4 * h4_m2 + fs_100_22869 * r_6 * h2_m2);
    }

    // NOTE: the rows are formed in 23 loops, as the vectorizer runs out of
    // registers with all 45 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_m1, ph2_0, ph2_p2, ph4_m1, ph4_0, ph4_p2, ph6_m1, ph6_0, ph6_p2, ph8_m1, ph8_0, ph8_p2, ab_2, pc_24, pc_25 : simd::cache_line_size())
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

        const auto h2_m1 = ph2_m1[k];
        const auto h2_0 = ph2_0[k];
        const auto h2_p2 = ph2_p2[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h4_0 = ph4_0[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h6_0 = ph6_0[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h8_m1 = ph8_m1[k];
        const auto h8_0 = ph8_0[k];
        const auto h8_p2 = ph8_p2[k];

        pc_24[k] = e_0 * (-f_85_8 * h2_0 + fs_1875_16 * h2_p2 - f_35_2 * r_2) + e_1 * (f_81_28 * h4_0 - fs_3645_196 * h4_p2 + f_255_28 * r_2 * h2_0 - fs_16875_196 * r_2 * h2_p2 + f_21_2 * r_4) + e_2 * (f_5_66 * h6_0 + fs_875_726 * h6_p2 - f_81_77 * r_2 * h4_0 + fs_14580_5929 * r_2 * h4_p2 - f_85_42 * r_4 * h2_0 + fs_625_147 * r_4 * h2_p2 - f_2 * r_6) + e_3 * (-f_392_1287 * h8_0 - fs_13720_184041 * h8_p2 - f_1_99 * r_2 * h6_0 - fs_70_3267 * r_2 * h6_p2 + f_81_1001 * r_4 * h4_0 - fs_14580_1002001 * r_4 * h4_p2 + f_85_693 * r_6 * h2_0 - fs_2500_160083 * r_6 * h2_p2 + f_1_9 * r_8) + f_105_16 * e_4;

        pc_25[k] = -fs_375_32 * e_0 * h2_m1 + e_1 * (f_81_28 * h4_m1 + fs_3375_392 * r_2 * h2_m1) + e_2 * (-fs_875_726 * h6_m1 - f_81_77 * r_2 * h4_m1 - fs_125_294 * r_4 * h2_m1) + e_3 * (fs_24010_184041 * h8_m1 + fs_70_3267 * r_2 * h6_m1 + f_81_1001 * r_4 * h4_m1 + fs_250_160083 * r_6 * h2_m1);
    }

    // NOTE: the rows are formed in 23 loops, as the vectorizer runs out of
    // registers with all 45 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph2_m2, ph2_m1, ph4_m3, ph4_m2, ph4_m1, ph6_m3, ph6_m2, ph6_m1, ph8_m3, ph8_m2, ph8_m1, ab_2, pc_26, pc_27 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h2_m2 = ph2_m2[k];
        const auto h2_m1 = ph2_m1[k];
        const auto h4_m3 = ph4_m3[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h8_m3 = ph8_m3[k];
        const auto h8_m2 = ph8_m2[k];
        const auto h8_m1 = ph8_m1[k];

        pc_26[k] = -fs_1875_16 * e_0 * h2_m2 + e_1 * (fs_3645_196 * h4_m2 + fs_16875_196 * r_2 * h2_m2) + e_2 * (-fs_875_726 * h6_m2 - fs_14580_5929 * r_2 * h4_m2 - fs_625_147 * r_4 * h2_m2) + e_3 * (fs_13720_184041 * h8_m2 + fs_70_3267 * r_2 * h6_m2 + fs_14580_1002001 * r_4 * h4_m2 + fs_2500_160083 * r_6 * h2_m2);

        pc_27[k] = fs_6075_128 * e_0 * h2_m1 + e_1 * (fs_405_112 * h4_m3 - fs_3645_196 * h4_m1 - fs_54675_1568 * r_2 * h2_m1) + e_2 * (-fs_875_1452 * h6_m3 + fs_525_968 * h6_m1 - fs_405_847 * r_2 * h4_m3 + fs_14580_5929 * r_2 * h4_m1 + fs_675_392 * r_4 * h2_m1) + e_3 * (fs_3430_50193 * h8_m3 + fs_4802_184041 * h8_m1 + fs_35_3267 * r_2 * h6_m3 - fs_7_726 * r_2 * h6_m1 + fs_405_143143 * r_4 * h4_m3 - fs_14580_1002001 * r_4 * h4_m1 - fs_75_11858 * r_6 * h2_m1);
    }

    // NOTE: the rows are formed in 23 loops, as the vectorizer runs out of
    // registers with all 45 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph2_m2, ph4_m4, ph4_m3, ph4_m2, ph6_m5, ph6_m4, ph6_m3, ph6_m2, ph8_m5, ph8_m4, ph8_m3, ph8_m2, ab_2, pc_28, pc_29 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h2_m2 = ph2_m2[k];
        const auto h4_m4 = ph4_m4[k];
        const auto h4_m3 = ph4_m3[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h6_m5 = ph6_m5[k];
        const auto h6_m4 = ph6_m4[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h8_m5 = ph8_m5[k];
        const auto h8_m4 = ph8_m4[k];
        const auto h8_m3 = ph8_m3[k];
        const auto h8_m2 = ph8_m2[k];

        pc_28[k] = -fs_4725_64 * e_0 * h2_m2 + e_1 * (-fs_405_16 * h4_m4 - fs_405_112 * h4_m2 + fs_6075_112 * r_2 * h2_m2) + e_2 * (f_5_22 * h6_m4 + fs_375_242 * h6_m2 + fs_405_121 * r_2 * h4_m4 + fs_405_847 * r_2 * h4_m2 - fs_75_28 * r_4 * h2_m2) + e_3 * (fs_784_16731 * h8_m4 + fs_1960_184041 * h8_m2 - f_1_33 * r_2 * h6_m4 - fs_10_363 * r_2 * h6_m2 - fs_405_20449 * r_4 * h4_m4 - fs_405_143143 * r_4 * h4_m2 + fs_25_2541 * r_6 * h2_m2);

        pc_29[k] = fs_405_16 * e_1 * h4_m3 + e_2 * (fs_25_11 * h6_m5 + fs_500_363 * h6_m3 - fs_405_121 * r_2 * h4_m3) + e_3 * (fs_49_2574 * h8_m5 + fs_245_100386 * h8_m3 - fs_4_99 * r_2 * h6_m5 - fs_80_3267 * r_2 * h6_m3 + fs_405_20449 * r_4 * h4_m3);
    }

    // NOTE: the rows are formed in 23 loops, as the vectorizer runs out of
    // registers with all 45 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_0, ph2_p1, ph2_p2, ph4_0, ph4_p1, ph4_p2, ph6_0, ph6_p1, ph8_0, ph8_p1, ph8_p2, ab_2, pc_30, pc_31, pc_32 : simd::cache_line_size())
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

        const auto h2_0 = ph2_0[k];
        const auto h2_p1 = ph2_p1[k];
        const auto h2_p2 = ph2_p2[k];
        const auto h4_0 = ph4_0[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h6_0 = ph6_0[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h8_0 = ph8_0[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h8_p2 = ph8_p2[k];

        pc_30[k] = e_0 * (-f_25_2 * h2_0 - f_35_2 * r_2) + e_1 * (f_81_14 * h4_0 + f_75_7 * r_2 * h2_0 + f_21_2 * r_4) + e_2 * (-f_50_33 * h6_0 - f_162_77 * r_2 * h4_0 - f_50_21 * r_4 * h2_0 - f_2 * r_6) + e_3 * (f_490_1287 * h8_0 + f_20_99 * r_2 * h6_0 + f_162_1001 * r_4 * h4_0 + f_100_693 * r_6 * h2_0 + f_1_9 * r_8) + f_105_16 * e_4;

        pc_31[k] = -fs_375_32 * e_0 * h2_p1 + e_1 * (f_81_28 * h4_p1 + fs_3375_392 * r_2 * h2_p1) + e_2 * (-fs_875_726 * h6_p1 - f_81_77 * r_2 * h4_p1 - fs_125_294 * r_4 * h2_p1) + e_3 * (fs_24010_184041 * h8_p1 + fs_70_3267 * r_2 * h6_p1 + f_81_1001 * r_4 * h4_p1 + fs_250_160083 * r_6 * h2_p1);

        pc_32[k] = fs_3375_16 * e_0 * h2_p2 + e_1 * (-f_99_28 * h4_p2 - fs_30375_196 * r_2 * h2_p2) + e_2 * (f_9_7 * r_2 * h4_p2 + fs_375_49 * r_4 * h2_p2) + e_3 * (fs_17150_184041 * h8_p2 - f_9_91 * r_4 * h4_p2 - fs_500_17787 * r_6 * h2_p2);
    }

    // NOTE: the rows are formed in 23 loops, as the vectorizer runs out of
    // registers with all 45 of them in one.

#pragma omp simd aligned(pe_1, pe_2, pe_3, ph4_p3, ph4_p4, ph6_p3, ph6_p4, ph8_p3, ph8_p4, ab_2, pc_33, pc_34 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

        const auto h4_p3 = ph4_p3[k];
        const auto h4_p4 = ph4_p4[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h6_p4 = ph6_p4[k];
        const auto h8_p3 = ph8_p3[k];
        const auto h8_p4 = ph8_p4[k];

        pc_33[k] = -f_27_4 * e_1 * h4_p3 + e_2 * (fs_625_363 * h6_p3 + f_27_11 * r_2 * h4_p3) + e_3 * (fs_2450_50193 * h8_p3 - fs_100_3267 * r_2 * h6_p3 - f_27_143 * r_4 * h4_p3);

        pc_34[k] = f_9_2 * e_1 * h4_p4 + e_2 * (fs_500_121 * h6_p4 - f_18_11 * r_2 * h4_p4) + e_3 * (fs_245_16731 * h8_p4 - fs_80_1089 * r_2 * h6_p4 + f_18_143 * r_4 * h4_p4);
    }

    // NOTE: the rows are formed in 23 loops, as the vectorizer runs out of
    // registers with all 45 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_0, ph2_p2, ph4_0, ph4_p2, ph6_0, ph6_p2, ph8_0, ph8_p2, ab_2, pc_35 : simd::cache_line_size())
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

        const auto h2_0 = ph2_0[k];
        const auto h2_p2 = ph2_p2[k];
        const auto h4_0 = ph4_0[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h6_0 = ph6_0[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h8_0 = ph8_0[k];
        const auto h8_p2 = ph8_p2[k];

        pc_35[k] = e_0 * (-f_85_8 * h2_0 - fs_1875_16 * h2_p2 - f_35_2 * r_2) + e_1 * (f_81_28 * h4_0 + fs_3645_196 * h4_p2 + f_255_28 * r_2 * h2_0 + fs_16875_196 * r_2 * h2_p2 + f_21_2 * r_4) + e_2 * (f_5_66 * h6_0 - fs_875_726 * h6_p2 - f_81_77 * r_2 * h4_0 - fs_14580_5929 * r_2 * h4_p2 - f_85_42 * r_4 * h2_0 - fs_625_147 * r_4 * h2_p2 - f_2 * r_6) + e_3 * (-f_392_1287 * h8_0 + fs_13720_184041 * h8_p2 - f_1_99 * r_2 * h6_0 + fs_70_3267 * r_2 * h6_p2 + f_81_1001 * r_4 * h4_0 + fs_14580_1002001 * r_4 * h4_p2 + f_85_693 * r_6 * h2_0 + fs_2500_160083 * r_6 * h2_p2 + f_1_9 * r_8) + f_105_16 * e_4;
    }

    // NOTE: the rows are formed in 23 loops, as the vectorizer runs out of
    // registers with all 45 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph2_p1, ph4_p1, ph4_p3, ph6_p1, ph6_p3, ph8_p1, ph8_p3, ab_2, pc_36 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h2_p1 = ph2_p1[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h4_p3 = ph4_p3[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h8_p3 = ph8_p3[k];

        pc_36[k] = -fs_6075_128 * e_0 * h2_p1 + e_1 * (fs_3645_196 * h4_p1 + fs_405_112 * h4_p3 + fs_54675_1568 * r_2 * h2_p1) + e_2 * (-fs_525_968 * h6_p1 - fs_875_1452 * h6_p3 - fs_14580_5929 * r_2 * h4_p1 - fs_405_847 * r_2 * h4_p3 - fs_675_392 * r_4 * h2_p1) + e_3 * (-fs_4802_184041 * h8_p1 + fs_3430_50193 * h8_p3 + fs_7_726 * r_2 * h6_p1 + fs_35_3267 * r_2 * h6_p3 + fs_14580_1002001 * r_4 * h4_p1 + fs_405_143143 * r_4 * h4_p3 + fs_75_11858 * r_6 * h2_p1);
    }

    // NOTE: the rows are formed in 23 loops, as the vectorizer runs out of
    // registers with all 45 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph2_p2, ph4_p2, ph4_p3, ph4_p4, ph6_p2, ph6_p3, ph6_p4, ph6_p5, ph8_p2, ph8_p3, ph8_p4, ph8_p5, ab_2, pc_37, pc_38 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h2_p2 = ph2_p2[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h4_p3 = ph4_p3[k];
        const auto h4_p4 = ph4_p4[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h6_p4 = ph6_p4[k];
        const auto h6_p5 = ph6_p5[k];
        const auto h8_p2 = ph8_p2[k];
        const auto h8_p3 = ph8_p3[k];
        const auto h8_p4 = ph8_p4[k];
        const auto h8_p5 = ph8_p5[k];

        pc_37[k] = fs_4725_64 * e_0 * h2_p2 + e_1 * (fs_405_112 * h4_p2 - fs_405_16 * h4_p4 - fs_6075_112 * r_2 * h2_p2) + e_2 * (-fs_375_242 * h6_p2 + f_5_22 * h6_p4 - fs_405_847 * r_2 * h4_p2 + fs_405_121 * r_2 * h4_p4 + fs_75_28 * r_4 * h2_p2) + e_3 * (-fs_1960_184041 * h8_p2 + fs_784_16731 * h8_p4 + fs_10_363 * r_2 * h6_p2 - f_1_33 * r_2 * h6_p4 + fs_405_143143 * r_4 * h4_p2 - fs_405_20449 * r_4 * h4_p4 - fs_25_2541 * r_6 * h2_p2);

        pc_38[k] = -fs_405_16 * e_1 * h4_p3 + e_2 * (-fs_500_363 * h6_p3 + fs_25_11 * h6_p5 + fs_405_121 * r_2 * h4_p3) + e_3 * (-fs_245_100386 * h8_p3 + fs_49_2574 * h8_p5 + fs_80_3267 * r_2 * h6_p3 - fs_4_99 * r_2 * h6_p5 - fs_405_20449 * r_4 * h4_p3);
    }

    // NOTE: the rows are formed in 23 loops, as the vectorizer runs out of
    // registers with all 45 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_0, ph4_0, ph4_p4, ph6_0, ph6_p4, ph8_0, ph8_p4, ab_2, pc_39 : simd::cache_line_size())
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

        const auto h2_0 = ph2_0[k];
        const auto h4_0 = ph4_0[k];
        const auto h4_p4 = ph4_p4[k];
        const auto h6_0 = ph6_0[k];
        const auto h6_p4 = ph6_p4[k];
        const auto h8_0 = ph8_0[k];
        const auto h8_p4 = ph8_p4[k];

        pc_39[k] = e_0 * (-f_5 * h2_0 - f_35_2 * r_2) + e_1 * (-f_99_28 * h4_0 + fs_3645_112 * h4_p4 + f_30_7 * r_2 * h2_0 + f_21_2 * r_4) + e_2 * (f_5_3 * h6_0 - fs_175_121 * h6_p4 + f_9_7 * r_2 * h4_0 - fs_3645_847 * r_2 * h4_p4 - f_20_21 * r_4 * h2_0 - f_2 * r_6) + e_3 * (f_196_1287 * h8_0 + fs_1372_16731 * h8_p4 - f_2_9 * r_2 * h6_0 + fs_28_1089 * r_2 * h6_p4 - f_9_91 * r_4 * h4_0 + fs_3645_143143 * r_4 * h4_p4 + f_40_693 * r_6 * h2_0 + f_1_9 * r_8) + f_105_16 * e_4;
    }

    // NOTE: the rows are formed in 23 loops, as the vectorizer runs out of
    // registers with all 45 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph2_p1, ph2_p2, ph4_p1, ph4_p2, ph6_p1, ph6_p2, ph6_p5, ph6_p6, ph8_p1, ph8_p2, ph8_p5, ph8_p6, ab_2, pc_40, pc_41 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h2_p1 = ph2_p1[k];
        const auto h2_p2 = ph2_p2[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h6_p5 = ph6_p5[k];
        const auto h6_p6 = ph6_p6[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h8_p2 = ph8_p2[k];
        const auto h8_p5 = ph8_p5[k];
        const auto h8_p6 = ph8_p6[k];

        pc_40[k] = -fs_13125_128 * e_0 * h2_p1 + e_1 * (fs_405_112 * h4_p1 + fs_16875_224 * r_2 * h2_p1) + e_2 * (fs_4225_2904 * h6_p1 - fs_25_44 * h6_p5 - fs_405_847 * r_2 * h4_p1 - fs_625_168 * r_4 * h2_p1) + e_3 * (fs_686_184041 * h8_p1 + fs_98_1287 * h8_p5 - fs_169_6534 * r_2 * h6_p1 + fs_1_99 * r_2 * h6_p5 + fs_405_143143 * r_4 * h4_p1 + fs_625_45738 * r_6 * h2_p1);

        pc_41[k] = fs_525_16 * e_0 * h2_p2 + e_1 * (fs_3645_112 * h4_p2 - fs_675_28 * r_2 * h2_p2) + e_2 * (fs_250_363 * h6_p2 + fs_50_33 * h6_p6 - fs_3645_847 * r_2 * h4_p2 + fs_25_21 * r_4 * h2_p2) + e_3 * (fs_245_368082 * h8_p2 + fs_343_7722 * h8_p6 - fs_40_3267 * r_2 * h6_p2 - fs_8_297 * r_2 * h6_p6 + fs_3645_143143 * r_4 * h4_p2 - fs_100_22869 * r_6 * h2_p2);
    }

    // NOTE: the rows are formed in 23 loops, as the vectorizer runs out of
    // registers with all 45 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_0, ph2_p1, ph4_0, ph4_p1, ph6_0, ph6_p1, ph6_p6, ph8_0, ph8_p1, ph8_p6, ph8_p7, ph8_p8, ab_2, pc_42, pc_43, pc_44 : simd::cache_line_size())
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

        const auto h2_0 = ph2_0[k];
        const auto h2_p1 = ph2_p1[k];
        const auto h4_0 = ph4_0[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h6_0 = ph6_0[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h6_p6 = ph6_p6[k];
        const auto h8_0 = ph8_0[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h8_p6 = ph8_p6[k];
        const auto h8_p7 = ph8_p7[k];
        const auto h8_p8 = ph8_p8[k];

        pc_42[k] = e_0 * (f_35_8 * h2_0 - f_35_2 * r_2) + e_1 * (-f_27_4 * h4_0 - f_15_4 * r_2 * h2_0 + f_21_2 * r_4) + e_2 * (-f_85_66 * h6_0 - fs_175_66 * h6_p6 + f_27_11 * r_2 * h4_0 + f_5_6 * r_4 * h2_0 - f_2 * r_6) + e_3 * (-f_56_1287 * h8_0 + fs_392_3861 * h8_p6 + f_17_99 * r_2 * h6_0 + fs_14_297 * r_2 * h6_p6 - f_27_143 * r_4 * h4_0 - f_5_99 * r_6 * h2_0 + f_1_9 * r_8) + f_105_16 * e_4;

        pc_43[k] = -fs_3675_32 * e_0 * h2_p1 + e_1 * (-fs_405_16 * h4_p1 + fs_675_8 * r_2 * h2_p1) + e_2 * (-fs_175_726 * h6_p1 + fs_405_121 * r_2 * h4_p1 - fs_25_6 * r_4 * h2_p1) + e_3 * (-fs_49_368082 * h8_p1 + fs_245_2574 * h8_p7 + fs_14_3267 * r_2 * h6_p1 - fs_405_20449 * r_4 * h4_p1 + fs_50_3267 * r_6 * h2_p1);

        pc_44[k] = e_0 * (f_35_2 * h2_0 - f_35_2 * r_2) + e_1 * (f_9_2 * h4_0 - f_15 * r_2 * h2_0 + f_21_2 * r_4) + e_2 * (f_10_33 * h6_0 - f_18_11 * r_2 * h4_0 + f_10_3 * r_4 * h2_0 - f_2 * r_6) + e_3 * (f_7_1287 * h8_0 + fs_245_1287 * h8_p8 - f_4_99 * r_2 * h6_0 + f_18_143 * r_4 * h4_0 - f_20_99 * r_6 * h2_0 + f_1_9 * r_8) + f_105_16 * e_4;
    }

    // NOTE: the values of a combination of angular components are stored as one
    // row of nvalues columns, with the component on bra side running slowest, and
    // the atom pairs beyond the reach of every pair of primitives are set to zero.

    const size_t sources[81] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 1, 9, 10, 11, 12, 13, 14, 15, 16, 2, 10, 17, 18, 19, 20, 21, 22, 23, 3, 11, 18, 24, 25, 26, 27, 28, 29, 4, 12, 19, 25, 30, 31, 32, 33, 34, 5, 13, 20, 26, 31, 35, 36, 37, 38, 6, 14, 21, 27, 32, 36, 39, 40, 41, 7, 15, 22, 28, 33, 37, 40, 42, 43, 8, 16, 23, 29, 34, 38, 41, 43, 44};

    for (size_t m = 0; m < 81; m++)
    {
        const auto *pc = buffer.data(5 + sources[m]);

        std::copy(pc, pc + nmax, values + m * nvalues);

        std::fill(values + m * nvalues + nmax, values + (m + 1) * nvalues, 0.0);
    }
}

}  // namespace simdovl
