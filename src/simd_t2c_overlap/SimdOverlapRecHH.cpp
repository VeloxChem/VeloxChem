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



#include "SimdOverlapRecHH.hpp"

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
compute_hh_overlap(double                         *values,
                   const size_t                    nvalues,
                   const CBasisFunction           &bra,
                   const CBasisFunction           &ket,
                   const std::vector<CSimdMatrix> &harmonics,
                   const CSimdMatrix              &coordinates,
                   const double                    threshold) -> void
{
    if ((bra.get_angular_momentum() != 5) || (ket.get_angular_momentum() != 5))
    {
        errors::assertMsgCritical(
            false, std::string("SimdOverlapRecHH.compute_hh_overlap: Basis functions must be of angular momenta five and five"));
    }

    if (harmonics.size() < 10)
    {
        errors::assertMsgCritical(
            false, std::string("SimdOverlapRecHH.compute_hh_overlap: Harmonics must reach angular momentum ten"));
    }

    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("SimdOverlapRecHH.compute_hh_overlap: Number of values exceeds number of atom pairs"));
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
        std::fill(values, values + 121 * nvalues, 0.0);

        return;
    }

    // NOTE: the first six rows accumulate the contracted prefactors of the terms,
    // and the remaining 66 rows hold the integrals of the combinations of angular
    // components which are not related by symmetry.

    auto buffer = CSimdMatrix(72, nmax);

    auto *pe_0 = buffer.data(0);
    auto *pe_1 = buffer.data(1);
    auto *pe_2 = buffer.data(2);
    auto *pe_3 = buffer.data(3);
    auto *pe_4 = buffer.data(4);
    auto *pe_5 = buffer.data(5);

    std::fill(pe_0, pe_0 + nmax, 0.0);
    std::fill(pe_1, pe_1 + nmax, 0.0);
    std::fill(pe_2, pe_2 + nmax, 0.0);
    std::fill(pe_3, pe_3 + nmax, 0.0);
    std::fill(pe_4, pe_4 + nmax, 0.0);
    std::fill(pe_5, pe_5 + nmax, 0.0);

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

            const auto f_0 = fbase * fmu / fexp / fexp / fexp / fexp / fexp;

            const auto f_1 = fbase * fmu * fmu / fexp / fexp / fexp / fexp / fexp;

            const auto f_2 = fbase * fmu * fmu * fmu / fexp / fexp / fexp / fexp / fexp;

            const auto f_3 = fbase * fmu * fmu * fmu * fmu / fexp / fexp / fexp / fexp / fexp;

            const auto f_4 = fbase * fmu * fmu * fmu * fmu * fmu / fexp / fexp / fexp / fexp / fexp;

            const auto f_5 = fbase / fexp / fexp / fexp / fexp / fexp;

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
    const auto *ph10_m10 = harmonics[9].data(0);
    const auto *ph10_m9 = harmonics[9].data(1);
    const auto *ph10_m8 = harmonics[9].data(2);
    const auto *ph10_m7 = harmonics[9].data(3);
    const auto *ph10_m6 = harmonics[9].data(4);
    const auto *ph10_m5 = harmonics[9].data(5);
    const auto *ph10_m4 = harmonics[9].data(6);
    const auto *ph10_m3 = harmonics[9].data(7);
    const auto *ph10_m2 = harmonics[9].data(8);
    const auto *ph10_m1 = harmonics[9].data(9);
    const auto *ph10_0 = harmonics[9].data(10);
    const auto *ph10_p1 = harmonics[9].data(11);
    const auto *ph10_p2 = harmonics[9].data(12);
    const auto *ph10_p3 = harmonics[9].data(13);
    const auto *ph10_p4 = harmonics[9].data(14);
    const auto *ph10_p5 = harmonics[9].data(15);
    const auto *ph10_p6 = harmonics[9].data(16);
    const auto *ph10_p7 = harmonics[9].data(17);
    const auto *ph10_p8 = harmonics[9].data(18);
    const auto *ph10_p9 = harmonics[9].data(19);
    const auto *ph10_p10 = harmonics[9].data(20);

    auto *pc_0 = buffer.data(6);
    auto *pc_1 = buffer.data(7);
    auto *pc_2 = buffer.data(8);
    auto *pc_3 = buffer.data(9);
    auto *pc_4 = buffer.data(10);
    auto *pc_5 = buffer.data(11);
    auto *pc_6 = buffer.data(12);
    auto *pc_7 = buffer.data(13);
    auto *pc_8 = buffer.data(14);
    auto *pc_9 = buffer.data(15);
    auto *pc_10 = buffer.data(16);
    auto *pc_11 = buffer.data(17);
    auto *pc_12 = buffer.data(18);
    auto *pc_13 = buffer.data(19);
    auto *pc_14 = buffer.data(20);
    auto *pc_15 = buffer.data(21);
    auto *pc_16 = buffer.data(22);
    auto *pc_17 = buffer.data(23);
    auto *pc_18 = buffer.data(24);
    auto *pc_19 = buffer.data(25);
    auto *pc_20 = buffer.data(26);
    auto *pc_21 = buffer.data(27);
    auto *pc_22 = buffer.data(28);
    auto *pc_23 = buffer.data(29);
    auto *pc_24 = buffer.data(30);
    auto *pc_25 = buffer.data(31);
    auto *pc_26 = buffer.data(32);
    auto *pc_27 = buffer.data(33);
    auto *pc_28 = buffer.data(34);
    auto *pc_29 = buffer.data(35);
    auto *pc_30 = buffer.data(36);
    auto *pc_31 = buffer.data(37);
    auto *pc_32 = buffer.data(38);
    auto *pc_33 = buffer.data(39);
    auto *pc_34 = buffer.data(40);
    auto *pc_35 = buffer.data(41);
    auto *pc_36 = buffer.data(42);
    auto *pc_37 = buffer.data(43);
    auto *pc_38 = buffer.data(44);
    auto *pc_39 = buffer.data(45);
    auto *pc_40 = buffer.data(46);
    auto *pc_41 = buffer.data(47);
    auto *pc_42 = buffer.data(48);
    auto *pc_43 = buffer.data(49);
    auto *pc_44 = buffer.data(50);
    auto *pc_45 = buffer.data(51);
    auto *pc_46 = buffer.data(52);
    auto *pc_47 = buffer.data(53);
    auto *pc_48 = buffer.data(54);
    auto *pc_49 = buffer.data(55);
    auto *pc_50 = buffer.data(56);
    auto *pc_51 = buffer.data(57);
    auto *pc_52 = buffer.data(58);
    auto *pc_53 = buffer.data(59);
    auto *pc_54 = buffer.data(60);
    auto *pc_55 = buffer.data(61);
    auto *pc_56 = buffer.data(62);
    auto *pc_57 = buffer.data(63);
    auto *pc_58 = buffer.data(64);
    auto *pc_59 = buffer.data(65);
    auto *pc_60 = buffer.data(66);
    auto *pc_61 = buffer.data(67);
    auto *pc_62 = buffer.data(68);
    auto *pc_63 = buffer.data(69);
    auto *pc_64 = buffer.data(70);
    auto *pc_65 = buffer.data(71);

    // NOTE: the factors of the terms depend on the angular momenta alone, so they
    // are formed once for the whole matrix instead of once for every atom pair.

    const auto f_1575_16 = 98.4375;
    const auto f_135_4 = 33.75;
    const auto f_225_2 = 112.5;
    const auto f_315_4 = 78.75;
    const auto f_75_22 = 75.0 / 22.0;
    const auto f_405_22 = 405.0 / 22.0;
    const auto f_75_2 = 37.5;
    const auto f_45_2 = 22.5;
    const auto f_35_286 = 35.0 / 286.0;
    const auto f_10_11 = 10.0 / 11.0;
    const auto f_405_143 = 405.0 / 143.0;
    const auto f_50_11 = 50.0 / 11.0;
    const auto f_5_2 = 2.5;
    const auto f_63_46189 = 63.0 / 46189.0;
    const auto fs_7938_46189 = std::sqrt(7938.0 / 46189.0);
    const auto f_35_2717 = 35.0 / 2717.0;
    const auto f_10_187 = 10.0 / 187.0;
    const auto f_18_143 = 18.0 / 143.0;
    const auto f_25_143 = 25.0 / 143.0;
    const auto f_1_11 = 1.0 / 11.0;
    const auto f_945_32 = 29.53125;
    const auto fs_1488375_512 = std::sqrt(2906.982421875);
    const auto fs_30375_8 = std::sqrt(3796.875);
    const auto fs_23625_968 = std::sqrt(23625.0 / 968.0);
    const auto fs_3375_8 = std::sqrt(421.875);
    const auto fs_2205_40898 = std::sqrt(2205.0 / 40898.0);
    const auto fs_210_121 = std::sqrt(210.0 / 121.0);
    const auto fs_750_121 = std::sqrt(750.0 / 121.0);
    const auto fs_3969_387895222 = std::sqrt(3969.0 / 387895222.0);
    const auto fs_3969_46189 = std::sqrt(3969.0 / 46189.0);
    const auto fs_4410_7382089 = std::sqrt(4410.0 / 7382089.0);
    const auto fs_210_34969 = std::sqrt(210.0 / 34969.0);
    const auto fs_375_40898 = std::sqrt(375.0 / 40898.0);
    const auto fs_165375_256 = std::sqrt(645.99609375);
    const auto fs_3375_4 = std::sqrt(843.75);
    const auto fs_13125_242 = std::sqrt(13125.0 / 242.0);
    const auto fs_375_4 = std::sqrt(93.75);
    const auto fs_8575_40898 = std::sqrt(8575.0 / 40898.0);
    const auto fs_1225_572 = std::sqrt(1225.0 / 572.0);
    const auto fs_1400_363 = std::sqrt(1400.0 / 363.0);
    const auto fs_500_363 = std::sqrt(500.0 / 363.0);
    const auto fs_11907_193947611 = std::sqrt(11907.0 / 193947611.0);
    const auto fs_35721_877591 = std::sqrt(35721.0 / 877591.0);
    const auto fs_17150_7382089 = std::sqrt(17150.0 / 7382089.0);
    const auto fs_1225_51623 = std::sqrt(1225.0 / 51623.0);
    const auto fs_1400_104907 = std::sqrt(1400.0 / 104907.0);
    const auto fs_125_61347 = std::sqrt(125.0 / 61347.0);
    const auto fs_42525_64 = std::sqrt(664.453125);
    const auto fs_39375_484 = std::sqrt(39375.0 / 484.0);
    const auto fs_382725_1936 = std::sqrt(382725.0 / 1936.0);
    const auto fs_8575_14872 = std::sqrt(8575.0 / 14872.0);
    const auto fs_3675_1144 = std::sqrt(3675.0 / 1144.0);
    const auto fs_700_121 = std::sqrt(700.0 / 121.0);
    const auto fs_382725_81796 = std::sqrt(382725.0 / 81796.0);
    const auto fs_3969_14919047 = std::sqrt(3969.0 / 14919047.0);
    const auto fs_15876_877591 = std::sqrt(15876.0 / 877591.0);
    const auto fs_8575_1342198 = std::sqrt(8575.0 / 1342198.0);
    const auto fs_3675_103246 = std::sqrt(3675.0 / 103246.0);
    const auto fs_700_34969 = std::sqrt(700.0 / 34969.0);
    const auto fs_189_20449 = std::sqrt(189.0 / 20449.0);
    const auto fs_6075_32 = std::sqrt(189.84375);
    const auto fs_84375_968 = std::sqrt(84375.0 / 968.0);
    const auto fs_1125_44 = std::sqrt(1125.0 / 44.0);
    const auto fs_54675_968 = std::sqrt(54675.0 / 968.0);
    const auto fs_18375_14872 = std::sqrt(18375.0 / 14872.0);
    const auto fs_1715_572 = std::sqrt(1715.0 / 572.0);
    const auto fs_20_11 = std::sqrt(20.0 / 11.0);
    const auto fs_54675_40898 = std::sqrt(54675.0 / 40898.0);
    const auto fs_27783_29838094 = std::sqrt(27783.0 / 29838094.0);
    const auto fs_111132_14919047 = std::sqrt(111132.0 / 14919047.0);
    const auto fs_18375_1342198 = std::sqrt(18375.0 / 1342198.0);
    const auto fs_1715_51623 = std::sqrt(1715.0 / 51623.0);
    const auto fs_750_34969 = std::sqrt(750.0 / 34969.0);
    const auto fs_20_3179 = std::sqrt(20.0 / 3179.0);
    const auto fs_54_20449 = std::sqrt(54.0 / 20449.0);
    const auto fs_5625_44 = std::sqrt(5625.0 / 44.0);
    const auto fs_1225_286 = std::sqrt(1225.0 / 286.0);
    const auto fs_100_11 = std::sqrt(100.0 / 11.0);
    const auto fs_83349_14919047 = std::sqrt(83349.0 / 14919047.0);
    const auto fs_2450_51623 = std::sqrt(2450.0 / 51623.0);
    const auto fs_100_3179 = std::sqrt(100.0 / 3179.0);
    const auto f_315_8 = 39.375;
    const auto f_45 = 45.0;
    const auto f_120_11 = 120.0 / 11.0;
    const auto f_15 = 15.0;
    const auto f_217_286 = 217.0 / 286.0;
    const auto fs_2205_572 = std::sqrt(2205.0 / 572.0);
    const auto f_32_11 = 32.0 / 11.0;
    const auto f_20_11 = 20.0 / 11.0;
    const auto f_630_46189 = 630.0 / 46189.0;
    const auto fs_79380_877591 = std::sqrt(79380.0 / 877591.0);
    const auto f_217_2717 = 217.0 / 2717.0;
    const auto fs_2205_51623 = std::sqrt(2205.0 / 51623.0);
    const auto f_32_187 = 32.0 / 187.0;
    const auto f_10_143 = 10.0 / 143.0;
    const auto fs_1620675_512 = std::sqrt(3165.380859375);
    const auto fs_33075_8 = std::sqrt(4134.375);
    const auto fs_525_8 = std::sqrt(65.625);
    const auto fs_3675_8 = std::sqrt(459.375);
    const auto fs_98_121 = std::sqrt(98.0 / 121.0);
    const auto fs_245_286 = std::sqrt(245.0 / 286.0);
    const auto fs_14_3 = std::sqrt(14.0 / 3.0);
    const auto fs_2450_363 = std::sqrt(2450.0 / 363.0);
    const auto fs_178605_387895222 = std::sqrt(178605.0 / 387895222.0);
    const auto fs_59535_877591 = std::sqrt(59535.0 / 877591.0);
    const auto fs_392_43681 = std::sqrt(392.0 / 43681.0);
    const auto fs_490_51623 = std::sqrt(490.0 / 51623.0);
    const auto fs_14_867 = std::sqrt(14.0 / 867.0);
    const auto fs_1225_122694 = std::sqrt(1225.0 / 122694.0);
    const auto fs_30375_64 = std::sqrt(474.609375);
    const auto fs_7875_242 = std::sqrt(7875.0 / 242.0);
    const auto fs_1575_22 = std::sqrt(1575.0 / 22.0);
    const auto fs_273375_1936 = std::sqrt(273375.0 / 1936.0);
    const auto fs_252105_163592 = std::sqrt(252105.0 / 163592.0);
    const auto fs_49_1144 = std::sqrt(49.0 / 1144.0);
    const auto fs_280_121 = std::sqrt(280.0 / 121.0);
    const auto fs_56_11 = std::sqrt(56.0 / 11.0);
    const auto fs_273375_81796 = std::sqrt(273375.0 / 81796.0);
    const auto fs_317520_193947611 = std::sqrt(317520.0 / 193947611.0);
    const auto fs_635040_14919047 = std::sqrt(635040.0 / 14919047.0);
    const auto fs_252105_14764178 = std::sqrt(252105.0 / 14764178.0);
    const auto fs_49_103246 = std::sqrt(49.0 / 103246.0);
    const auto fs_280_34969 = std::sqrt(280.0 / 34969.0);
    const auto fs_56_3179 = std::sqrt(56.0 / 3179.0);
    const auto fs_135_20449 = std::sqrt(135.0 / 20449.0);
    const auto fs_30375_32 = std::sqrt(949.21875);
    const auto fs_1125_968 = std::sqrt(1125.0 / 968.0);
    const auto fs_6075_88 = std::sqrt(6075.0 / 88.0);
    const auto fs_273375_968 = std::sqrt(273375.0 / 968.0);
    const auto fs_3920_1859 = std::sqrt(3920.0 / 1859.0);
    const auto fs_147_143 = std::sqrt(147.0 / 143.0);
    const auto fs_10_121 = std::sqrt(10.0 / 121.0);
    const auto fs_54_11 = std::sqrt(54.0 / 11.0);
    const auto fs_273375_40898 = std::sqrt(273375.0 / 40898.0);
    const auto fs_138915_29838094 = std::sqrt(138915.0 / 29838094.0);
    const auto fs_694575_29838094 = std::sqrt(694575.0 / 29838094.0);
    const auto fs_15680_671099 = std::sqrt(15680.0 / 671099.0);
    const auto fs_588_51623 = std::sqrt(588.0 / 51623.0);
    const auto fs_10_34969 = std::sqrt(10.0 / 34969.0);
    const auto fs_54_3179 = std::sqrt(54.0 / 3179.0);
    const auto fs_270_20449 = std::sqrt(270.0 / 20449.0);
    const auto fs_4500_121 = std::sqrt(4500.0 / 121.0);
    const auto fs_2695_676 = std::sqrt(2695.0 / 676.0);
    const auto fs_320_121 = std::sqrt(320.0 / 121.0);
    const auto fs_333396_14919047 = std::sqrt(333396.0 / 14919047.0);
    const auto fs_2695_61009 = std::sqrt(2695.0 / 61009.0);
    const auto fs_320_34969 = std::sqrt(320.0 / 34969.0);
    const auto f_105_16 = 6.5625;
    const auto f_15_2 = 7.5;
    const auto f_145_22 = 145.0 / 22.0;
    const auto fs_1050_11 = std::sqrt(1050.0 / 11.0);
    const auto f_511_286 = 511.0 / 286.0;
    const auto fs_294_143 = std::sqrt(294.0 / 143.0);
    const auto f_58_33 = 58.0 / 33.0;
    const auto fs_224_33 = std::sqrt(224.0 / 33.0);
    const auto f_10_33 = 10.0 / 33.0;
    const auto f_2835_46189 = 2835.0 / 46189.0;
    const auto fs_1071630_14919047 = std::sqrt(1071630.0 / 14919047.0);
    const auto f_511_2717 = 511.0 / 2717.0;
    const auto fs_1176_51623 = std::sqrt(1176.0 / 51623.0);
    const auto f_58_561 = 58.0 / 561.0;
    const auto fs_224_9537 = std::sqrt(224.0 / 9537.0);
    const auto f_5_429 = 5.0 / 429.0;
    const auto fs_275625_128 = std::sqrt(2153.3203125);
    const auto fs_5625_2 = std::sqrt(2812.5);
    const auto fs_175_242 = std::sqrt(175.0 / 242.0);
    const auto fs_525_44 = std::sqrt(525.0 / 44.0);
    const auto fs_625_2 = std::sqrt(312.5);
    const auto fs_324723_163592 = std::sqrt(324723.0 / 163592.0);
    const auto fs_1029_1144 = std::sqrt(1029.0 / 1144.0);
    const auto fs_56_1089 = std::sqrt(56.0 / 1089.0);
    const auto fs_28_33 = std::sqrt(28.0 / 33.0);
    const auto fs_5000_1089 = std::sqrt(5000.0 / 1089.0);
    const auto fs_1071630_193947611 = std::sqrt(1071630.0 / 193947611.0);
    const auto fs_893025_14919047 = std::sqrt(893025.0 / 14919047.0);
    const auto fs_324723_14764178 = std::sqrt(324723.0 / 14764178.0);
    const auto fs_1029_103246 = std::sqrt(1029.0 / 103246.0);
    const auto fs_56_314721 = std::sqrt(56.0 / 314721.0);
    const auto fs_28_9537 = std::sqrt(28.0 / 9537.0);
    const auto fs_1250_184041 = std::sqrt(1250.0 / 184041.0);
    const auto fs_77175_32 = std::sqrt(2411.71875);
    const auto fs_3150 = std::sqrt(3150.0);
    const auto fs_125_4 = std::sqrt(31.25);
    const auto fs_12675_968 = std::sqrt(12675.0 / 968.0);
    const auto fs_350 = std::sqrt(350.0);
    const auto fs_735_484 = std::sqrt(735.0 / 484.0);
    const auto fs_147_14872 = std::sqrt(147.0 / 14872.0);
    const auto fs_20_9 = std::sqrt(20.0 / 9.0);
    const auto fs_338_363 = std::sqrt(338.0 / 363.0);
    const auto fs_5600_1089 = std::sqrt(5600.0 / 1089.0);
    const auto fs_2500470_193947611 = std::sqrt(2500470.0 / 193947611.0);
    const auto fs_1250235_29838094 = std::sqrt(1250235.0 / 29838094.0);
    const auto fs_735_43681 = std::sqrt(735.0 / 43681.0);
    const auto fs_147_1342198 = std::sqrt(147.0 / 1342198.0);
    const auto fs_20_2601 = std::sqrt(20.0 / 2601.0);
    const auto fs_338_104907 = std::sqrt(338.0 / 104907.0);
    const auto fs_1400_184041 = std::sqrt(1400.0 / 184041.0);
    const auto fs_46875_484 = std::sqrt(46875.0 / 484.0);
    const auto fs_3675_3718 = std::sqrt(3675.0 / 3718.0);
    const auto fs_2500_363 = std::sqrt(2500.0 / 363.0);
    const auto fs_750141_14919047 = std::sqrt(750141.0 / 14919047.0);
    const auto fs_7350_671099 = std::sqrt(7350.0 / 671099.0);
    const auto fs_2500_104907 = std::sqrt(2500.0 / 104907.0);
    const auto f_45_8 = 5.625;
    const auto fs_70875_64 = std::sqrt(1107.421875);
    const auto f_90_11 = 90.0 / 11.0;
    const auto fs_6300_121 = std::sqrt(6300.0 / 121.0);
    const auto f_135_44 = 135.0 / 44.0;
    const auto fs_637875_1936 = std::sqrt(637875.0 / 1936.0);
    const auto f_238_143 = 238.0 / 143.0;
    const auto fs_3087_1859 = std::sqrt(3087.0 / 1859.0);
    const auto f_24_11 = 24.0 / 11.0;
    const auto fs_448_121 = std::sqrt(448.0 / 121.0);
    const auto f_135_286 = 135.0 / 286.0;
    const auto fs_637875_81796 = std::sqrt(637875.0 / 81796.0);
    const auto f_7560_46189 = 7560.0 / 46189.0;
    const auto fs_952560_14919047 = std::sqrt(952560.0 / 14919047.0);
    const auto f_476_2717 = 476.0 / 2717.0;
    const auto fs_12348_671099 = std::sqrt(12348.0 / 671099.0);
    const auto f_24_187 = 24.0 / 187.0;
    const auto fs_448_34969 = std::sqrt(448.0 / 34969.0);
    const auto f_3_143 = 3.0 / 143.0;
    const auto fs_315_20449 = std::sqrt(315.0 / 20449.0);
    const auto fs_231525_256 = std::sqrt(904.39453125);
    const auto fs_70875_128 = std::sqrt(553.7109375);
    const auto fs_10125_128 = std::sqrt(79.1015625);
    const auto fs_4725_4 = std::sqrt(1181.25);
    const auto fs_4800_121 = std::sqrt(4800.0 / 121.0);
    const auto fs_3375_242 = std::sqrt(3375.0 / 242.0);
    const auto fs_637875_3872 = std::sqrt(637875.0 / 3872.0);
    const auto fs_91125_3872 = std::sqrt(91125.0 / 3872.0);
    const auto fs_525_4 = std::sqrt(131.25);
    const auto fs_27783_81796 = std::sqrt(27783.0 / 81796.0);
    const auto fs_6615_7436 = std::sqrt(6615.0 / 7436.0);
    const auto fs_1024_363 = std::sqrt(1024.0 / 363.0);
    const auto fs_120_121 = std::sqrt(120.0 / 121.0);
    const auto fs_637875_163592 = std::sqrt(637875.0 / 163592.0);
    const auto fs_91125_163592 = std::sqrt(91125.0 / 163592.0);
    const auto fs_700_363 = std::sqrt(700.0 / 363.0);
    const auto fs_5000940_193947611 = std::sqrt(5000940.0 / 193947611.0);
    const auto fs_833490_14919047 = std::sqrt(833490.0 / 14919047.0);
    const auto fs_27783_7382089 = std::sqrt(27783.0 / 7382089.0);
    const auto fs_6615_671099 = std::sqrt(6615.0 / 671099.0);
    const auto fs_1024_104907 = std::sqrt(1024.0 / 104907.0);
    const auto fs_120_34969 = std::sqrt(120.0 / 34969.0);
    const auto fs_315_40898 = std::sqrt(315.0 / 40898.0);
    const auto fs_45_40898 = std::sqrt(45.0 / 40898.0);
    const auto fs_175_61347 = std::sqrt(175.0 / 61347.0);
    const auto fs_385875_64 = std::sqrt(6029.296875);
    const auto fs_7875 = std::sqrt(7875.0);
    const auto fs_1250_121 = std::sqrt(1250.0 / 121.0);
    const auto fs_875 = std::sqrt(875.0);
    const auto fs_3675_40898 = std::sqrt(3675.0 / 40898.0);
    const auto fs_800_1089 = std::sqrt(800.0 / 1089.0);
    const auto fs_14000_1089 = std::sqrt(14000.0 / 1089.0);
    const auto fs_16003008_193947611 = std::sqrt(16003008.0 / 193947611.0);
    const auto fs_7350_7382089 = std::sqrt(7350.0 / 7382089.0);
    const auto fs_800_314721 = std::sqrt(800.0 / 314721.0);
    const auto fs_3500_184041 = std::sqrt(3500.0 / 184041.0);
    const auto f_945_16 = 59.0625;
    const auto fs_826875_256 = std::sqrt(3229.98046875);
    const auto fs_10125_16 = std::sqrt(632.8125);
    const auto f_135_2 = 67.5;
    const auto fs_16875_4 = std::sqrt(4218.75);
    const auto f_30_11 = 30.0 / 11.0;
    const auto fs_5250_121 = std::sqrt(5250.0 / 121.0);
    const auto f_135_11 = 135.0 / 11.0;
    const auto fs_91125_484 = std::sqrt(91125.0 / 484.0);
    const auto fs_1875_4 = std::sqrt(468.75);
    const auto f_49_143 = 49.0 / 143.0;
    const auto fs_30870_20449 = std::sqrt(30870.0 / 20449.0);
    const auto f_8_11 = 8.0 / 11.0;
    const auto fs_1120_363 = std::sqrt(1120.0 / 363.0);
    const auto f_270_143 = 270.0 / 143.0;
    const auto fs_91125_20449 = std::sqrt(91125.0 / 20449.0);
    const auto f_13230_46189 = 13230.0 / 46189.0;
    const auto fs_11668860_193947611 = std::sqrt(11668860.0 / 193947611.0);
    const auto f_98_2717 = 98.0 / 2717.0;
    const auto fs_123480_7382089 = std::sqrt(123480.0 / 7382089.0);
    const auto f_8_187 = 8.0 / 187.0;
    const auto fs_1120_104907 = std::sqrt(1120.0 / 104907.0);
    const auto f_12_143 = 12.0 / 143.0;
    const auto fs_180_20449 = std::sqrt(180.0 / 20449.0);
    const auto f_15_143 = 15.0 / 143.0;
    const auto fs_625_61347 = std::sqrt(625.0 / 61347.0);
    const auto fs_55125_256 = std::sqrt(215.33203125);
    const auto fs_1125_4 = std::sqrt(281.25);
    const auto fs_3500_121 = std::sqrt(3500.0 / 121.0);
    const auto fs_36015_20449 = std::sqrt(36015.0 / 20449.0);
    const auto fs_2240_1089 = std::sqrt(2240.0 / 1089.0);
    const auto fs_500_1089 = std::sqrt(500.0 / 1089.0);
    const auto fs_21003948_193947611 = std::sqrt(21003948.0 / 193947611.0);
    const auto fs_144060_7382089 = std::sqrt(144060.0 / 7382089.0);
    const auto fs_2240_314721 = std::sqrt(2240.0 / 314721.0);
    const auto fs_125_184041 = std::sqrt(125.0 / 184041.0);
    const auto f_525_8 = 65.625;
    const auto f_75 = 75.0;
    const auto f_100_11 = 100.0 / 11.0;
    const auto f_25 = 25.0;
    const auto f_245_143 = 245.0 / 143.0;
    const auto f_80_33 = 80.0 / 33.0;
    const auto f_100_33 = 100.0 / 33.0;
    const auto f_15876_46189 = 15876.0 / 46189.0;
    const auto f_490_2717 = 490.0 / 2717.0;
    const auto f_80_561 = 80.0 / 561.0;
    const auto f_50_429 = 50.0 / 429.0;

    // NOTE: the rows are formed in 53 loops, as the vectorizer runs out of
    // registers with all 66 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph2_0, ph4_0, ph6_0, ph8_0, ph10_0, ph10_p10, ab_2, pc_0 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;
        const auto r_10 = r_8 * r_2;

        const auto h2_0 = ph2_0[k];
        const auto h4_0 = ph4_0[k];
        const auto h6_0 = ph6_0[k];
        const auto h8_0 = ph8_0[k];
        const auto h10_0 = ph10_0[k];
        const auto h10_p10 = ph10_p10[k];

        pc_0[k] = e_0 * (f_1575_16 * h2_0 - f_1575_16 * r_2) + e_1 * (f_135_4 * h4_0 - f_225_2 * r_2 * h2_0 + f_315_4 * r_4) + e_2 * (f_75_22 * h6_0 - f_405_22 * r_2 * h4_0 + f_75_2 * r_4 * h2_0 - f_45_2 * r_6) + e_3 * (f_35_286 * h8_0 - f_10_11 * r_2 * h6_0 + f_405_143 * r_4 * h4_0 - f_50_11 * r_6 * h2_0 + f_5_2 * r_8) + e_4 * (f_63_46189 * h10_0 + fs_7938_46189 * h10_p10 - f_35_2717 * r_2 * h8_0 + f_10_187 * r_4 * h6_0 - f_18_143 * r_6 * h4_0 + f_25_143 * r_8 * h2_0 - f_1_11 * r_10) + f_945_32 * e_5;
    }

    // NOTE: the rows are formed in 53 loops, as the vectorizer runs out of
    // registers with all 66 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_p1, ph4_p1, ph6_p1, ph8_p1, ph10_p1, ph10_p9, ab_2, pc_1 : simd::cache_line_size())
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

        const auto h2_p1 = ph2_p1[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h10_p1 = ph10_p1[k];
        const auto h10_p9 = ph10_p9[k];

        pc_1[k] = -fs_1488375_512 * e_0 * h2_p1 + e_1 * (-f_135_4 * h4_p1 + fs_30375_8 * r_2 * h2_p1) + e_2 * (-fs_23625_968 * h6_p1 + f_405_22 * r_2 * h4_p1 - fs_3375_8 * r_4 * h2_p1) + e_3 * (-fs_2205_40898 * h8_p1 + fs_210_121 * r_2 * h6_p1 - f_405_143 * r_4 * h4_p1 + fs_750_121 * r_6 * h2_p1) + e_4 * (-fs_3969_387895222 * h10_p1 + fs_3969_46189 * h10_p9 + fs_4410_7382089 * r_2 * h8_p1 - fs_210_34969 * r_4 * h6_p1 + f_18_143 * r_6 * h4_p1 - fs_375_40898 * r_8 * h2_p1);
    }

    // NOTE: the rows are formed in 53 loops, as the vectorizer runs out of
    // registers with all 66 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_p2, ph4_p2, ph6_p2, ph8_p2, ph8_p8, ph10_p2, ph10_p8, ab_2, pc_2 : simd::cache_line_size())
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

        const auto h2_p2 = ph2_p2[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h8_p2 = ph8_p2[k];
        const auto h8_p8 = ph8_p8[k];
        const auto h10_p2 = ph10_p2[k];
        const auto h10_p8 = ph10_p8[k];

        pc_2[k] = fs_165375_256 * e_0 * h2_p2 + e_1 * (f_135_4 * h4_p2 - fs_3375_4 * r_2 * h2_p2) + e_2 * (fs_13125_242 * h6_p2 - f_405_22 * r_2 * h4_p2 + fs_375_4 * r_4 * h2_p2) + e_3 * (fs_8575_40898 * h8_p2 + fs_1225_572 * h8_p8 - fs_1400_363 * r_2 * h6_p2 + f_405_143 * r_4 * h4_p2 - fs_500_363 * r_6 * h2_p2) + e_4 * (fs_11907_193947611 * h10_p2 + fs_35721_877591 * h10_p8 - fs_17150_7382089 * r_2 * h8_p2 - fs_1225_51623 * r_2 * h8_p8 + fs_1400_104907 * r_4 * h6_p2 - f_18_143 * r_6 * h4_p2 + fs_125_61347 * r_8 * h2_p2);
    }

    // NOTE: the rows are formed in 53 loops, as the vectorizer runs out of
    // registers with all 66 of them in one.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, ph4_p3, ph4_p4, ph6_p3, ph6_p4, ph6_p6, ph8_p3, ph8_p4, ph8_p6, ph8_p7, ph10_p3, ph10_p4, ph10_p6, ph10_p7, ab_2, pc_3, pc_4 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h4_p3 = ph4_p3[k];
        const auto h4_p4 = ph4_p4[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h6_p4 = ph6_p4[k];
        const auto h6_p6 = ph6_p6[k];
        const auto h8_p3 = ph8_p3[k];
        const auto h8_p4 = ph8_p4[k];
        const auto h8_p6 = ph8_p6[k];
        const auto h8_p7 = ph8_p7[k];
        const auto h10_p3 = ph10_p3[k];
        const auto h10_p4 = ph10_p4[k];
        const auto h10_p6 = ph10_p6[k];
        const auto h10_p7 = ph10_p7[k];

        pc_3[k] = -fs_42525_64 * e_1 * h4_p3 + e_2 * (-fs_39375_484 * h6_p3 + fs_382725_1936 * r_2 * h4_p3) + e_3 * (-fs_8575_14872 * h8_p3 + fs_3675_1144 * h8_p7 + fs_700_121 * r_2 * h6_p3 - fs_382725_81796 * r_4 * h4_p3) + e_4 * (-fs_3969_14919047 * h10_p3 + fs_15876_877591 * h10_p7 + fs_8575_1342198 * r_2 * h8_p3 - fs_3675_103246 * r_2 * h8_p7 - fs_700_34969 * r_4 * h6_p3 + fs_189_20449 * r_6 * h4_p3);

        pc_4[k] = fs_6075_32 * e_1 * h4_p4 + e_2 * (fs_84375_968 * h6_p4 + fs_1125_44 * h6_p6 - fs_54675_968 * r_2 * h4_p4) + e_3 * (fs_18375_14872 * h8_p4 + fs_1715_572 * h8_p6 - fs_750_121 * r_2 * h6_p4 - fs_20_11 * r_2 * h6_p6 + fs_54675_40898 * r_4 * h4_p4) + e_4 * (fs_27783_29838094 * h10_p4 + fs_111132_14919047 * h10_p6 - fs_18375_1342198 * r_2 * h8_p4 - fs_1715_51623 * r_2 * h8_p6 + fs_750_34969 * r_4 * h6_p4 + fs_20_3179 * r_4 * h6_p6 - fs_54_20449 * r_6 * h4_p4);
    }

    // NOTE: the rows are formed in 53 loops, as the vectorizer runs out of
    // registers with all 66 of them in one.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, ph4_m4, ph6_m6, ph6_m5, ph6_m4, ph8_m6, ph8_m5, ph8_m4, ph10_m6, ph10_m5, ph10_m4, ab_2, pc_5, pc_6 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h4_m4 = ph4_m4[k];
        const auto h6_m6 = ph6_m6[k];
        const auto h6_m5 = ph6_m5[k];
        const auto h6_m4 = ph6_m4[k];
        const auto h8_m6 = ph8_m6[k];
        const auto h8_m5 = ph8_m5[k];
        const auto h8_m4 = ph8_m4[k];
        const auto h10_m6 = ph10_m6[k];
        const auto h10_m5 = ph10_m5[k];
        const auto h10_m4 = ph10_m4[k];

        pc_5[k] = -fs_5625_44 * e_2 * h6_m5 + e_3 * (-fs_1225_286 * h8_m5 + fs_100_11 * r_2 * h6_m5) + e_4 * (-fs_83349_14919047 * h10_m5 + fs_2450_51623 * r_2 * h8_m5 - fs_100_3179 * r_4 * h6_m5);

        pc_6[k] = fs_6075_32 * e_1 * h4_m4 + e_2 * (-fs_1125_44 * h6_m6 + fs_84375_968 * h6_m4 - fs_54675_968 * r_2 * h4_m4) + e_3 * (-fs_1715_572 * h8_m6 + fs_18375_14872 * h8_m4 + fs_20_11 * r_2 * h6_m6 - fs_750_121 * r_2 * h6_m4 + fs_54675_40898 * r_4 * h4_m4) + e_4 * (-fs_111132_14919047 * h10_m6 + fs_27783_29838094 * h10_m4 + fs_1715_51623 * r_2 * h8_m6 - fs_18375_1342198 * r_2 * h8_m4 - fs_20_3179 * r_4 * h6_m6 + fs_750_34969 * r_4 * h6_m4 - fs_54_20449 * r_6 * h4_m4);
    }

    // NOTE: the rows are formed in 53 loops, as the vectorizer runs out of
    // registers with all 66 of them in one.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, ph4_m3, ph6_m3, ph8_m7, ph8_m3, ph10_m7, ph10_m3, ab_2, pc_7 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h4_m3 = ph4_m3[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h8_m7 = ph8_m7[k];
        const auto h8_m3 = ph8_m3[k];
        const auto h10_m7 = ph10_m7[k];
        const auto h10_m3 = ph10_m3[k];

        pc_7[k] = -fs_42525_64 * e_1 * h4_m3 + e_2 * (-fs_39375_484 * h6_m3 + fs_382725_1936 * r_2 * h4_m3) + e_3 * (-fs_3675_1144 * h8_m7 - fs_8575_14872 * h8_m3 + fs_700_121 * r_2 * h6_m3 - fs_382725_81796 * r_4 * h4_m3) + e_4 * (-fs_15876_877591 * h10_m7 - fs_3969_14919047 * h10_m3 + fs_3675_103246 * r_2 * h8_m7 + fs_8575_1342198 * r_2 * h8_m3 - fs_700_34969 * r_4 * h6_m3 + fs_189_20449 * r_6 * h4_m3);
    }

    // NOTE: the rows are formed in 53 loops, as the vectorizer runs out of
    // registers with all 66 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_m2, ph4_m2, ph6_m2, ph8_m8, ph8_m2, ph10_m8, ph10_m2, ab_2, pc_8 : simd::cache_line_size())
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

        const auto h2_m2 = ph2_m2[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h8_m8 = ph8_m8[k];
        const auto h8_m2 = ph8_m2[k];
        const auto h10_m8 = ph10_m8[k];
        const auto h10_m2 = ph10_m2[k];

        pc_8[k] = fs_165375_256 * e_0 * h2_m2 + e_1 * (f_135_4 * h4_m2 - fs_3375_4 * r_2 * h2_m2) + e_2 * (fs_13125_242 * h6_m2 - f_405_22 * r_2 * h4_m2 + fs_375_4 * r_4 * h2_m2) + e_3 * (-fs_1225_572 * h8_m8 + fs_8575_40898 * h8_m2 - fs_1400_363 * r_2 * h6_m2 + f_405_143 * r_4 * h4_m2 - fs_500_363 * r_6 * h2_m2) + e_4 * (-fs_35721_877591 * h10_m8 + fs_11907_193947611 * h10_m2 + fs_1225_51623 * r_2 * h8_m8 - fs_17150_7382089 * r_2 * h8_m2 + fs_1400_104907 * r_4 * h6_m2 - f_18_143 * r_6 * h4_m2 + fs_125_61347 * r_8 * h2_m2);
    }

    // NOTE: the rows are formed in 53 loops, as the vectorizer runs out of
    // registers with all 66 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_m1, ph4_m1, ph6_m1, ph8_m1, ph10_m10, ph10_m9, ph10_m1, ab_2, pc_9, pc_10 : simd::cache_line_size())
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
        const auto h4_m1 = ph4_m1[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h8_m1 = ph8_m1[k];
        const auto h10_m10 = ph10_m10[k];
        const auto h10_m9 = ph10_m9[k];
        const auto h10_m1 = ph10_m1[k];

        pc_9[k] = -fs_1488375_512 * e_0 * h2_m1 + e_1 * (-f_135_4 * h4_m1 + fs_30375_8 * r_2 * h2_m1) + e_2 * (-fs_23625_968 * h6_m1 + f_405_22 * r_2 * h4_m1 - fs_3375_8 * r_4 * h2_m1) + e_3 * (-fs_2205_40898 * h8_m1 + fs_210_121 * r_2 * h6_m1 - f_405_143 * r_4 * h4_m1 + fs_750_121 * r_6 * h2_m1) + e_4 * (-fs_3969_46189 * h10_m9 - fs_3969_387895222 * h10_m1 + fs_4410_7382089 * r_2 * h8_m1 - fs_210_34969 * r_4 * h6_m1 + f_18_143 * r_6 * h4_m1 - fs_375_40898 * r_8 * h2_m1);

        pc_10[k] = -fs_7938_46189 * e_4 * h10_m10;
    }

    // NOTE: the rows are formed in 53 loops, as the vectorizer runs out of
    // registers with all 66 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph2_0, ph4_0, ph6_0, ph8_0, ph8_p8, ph10_0, ph10_p8, ab_2, pc_11 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;
        const auto r_10 = r_8 * r_2;

        const auto h2_0 = ph2_0[k];
        const auto h4_0 = ph4_0[k];
        const auto h6_0 = ph6_0[k];
        const auto h8_0 = ph8_0[k];
        const auto h8_p8 = ph8_p8[k];
        const auto h10_0 = ph10_0[k];
        const auto h10_p8 = ph10_p8[k];

        pc_11[k] = e_0 * (f_315_8 * h2_0 - f_1575_16 * r_2) + e_1 * (-f_135_4 * h4_0 - f_45 * r_2 * h2_0 + f_315_4 * r_4) + e_2 * (-f_120_11 * h6_0 + f_405_22 * r_2 * h4_0 + f_15 * r_4 * h2_0 - f_45_2 * r_6) + e_3 * (-f_217_286 * h8_0 - fs_2205_572 * h8_p8 + f_32_11 * r_2 * h6_0 - f_405_143 * r_4 * h4_0 - f_20_11 * r_6 * h2_0 + f_5_2 * r_8) + e_4 * (-f_630_46189 * h10_0 + fs_79380_877591 * h10_p8 + f_217_2717 * r_2 * h8_0 + fs_2205_51623 * r_2 * h8_p8 - f_32_187 * r_4 * h6_0 + f_18_143 * r_6 * h4_0 + f_10_143 * r_8 * h2_0 - f_1_11 * r_10) + f_945_32 * e_5;
    }

    // NOTE: the rows are formed in 53 loops, as the vectorizer runs out of
    // registers with all 66 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_p1, ph6_p1, ph8_p1, ph8_p7, ph10_p1, ph10_p7, ab_2, pc_12 : simd::cache_line_size())
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

        const auto h2_p1 = ph2_p1[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h8_p7 = ph8_p7[k];
        const auto h10_p1 = ph10_p1[k];
        const auto h10_p7 = ph10_p7[k];

        pc_12[k] = -fs_1620675_512 * e_0 * h2_p1 + fs_33075_8 * e_1 * r_2 * h2_p1 + e_2 * (fs_525_8 * h6_p1 - fs_3675_8 * r_4 * h2_p1) + e_3 * (fs_98_121 * h8_p1 - fs_245_286 * h8_p7 - fs_14_3 * r_2 * h6_p1 + fs_2450_363 * r_6 * h2_p1) + e_4 * (fs_178605_387895222 * h10_p1 + fs_59535_877591 * h10_p7 - fs_392_43681 * r_2 * h8_p1 + fs_490_51623 * r_2 * h8_p7 + fs_14_867 * r_4 * h6_p1 - fs_1225_122694 * r_8 * h2_p1);
    }

    // NOTE: the rows are formed in 53 loops, as the vectorizer runs out of
    // registers with all 66 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_p2, ph4_p2, ph6_p2, ph6_p6, ph8_p2, ph8_p6, ph10_p2, ph10_p6, ab_2, pc_13 : simd::cache_line_size())
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

        const auto h2_p2 = ph2_p2[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h6_p6 = ph6_p6[k];
        const auto h8_p2 = ph8_p2[k];
        const auto h8_p6 = ph8_p6[k];
        const auto h10_p2 = ph10_p2[k];
        const auto h10_p6 = ph10_p6[k];

        pc_13[k] = f_315_8 * e_0 * h2_p2 + e_1 * (fs_30375_64 * h4_p2 - f_45 * r_2 * h2_p2) + e_2 * (-fs_7875_242 * h6_p2 - fs_1575_22 * h6_p6 - fs_273375_1936 * r_2 * h4_p2 + f_15 * r_4 * h2_p2) + e_3 * (-fs_252105_163592 * h8_p2 + fs_49_1144 * h8_p6 + fs_280_121 * r_2 * h6_p2 + fs_56_11 * r_2 * h6_p6 + fs_273375_81796 * r_4 * h4_p2 - f_20_11 * r_6 * h2_p2) + e_4 * (-fs_317520_193947611 * h10_p2 + fs_635040_14919047 * h10_p6 + fs_252105_14764178 * r_2 * h8_p2 - fs_49_103246 * r_2 * h8_p6 - fs_280_34969 * r_4 * h6_p2 - fs_56_3179 * r_4 * h6_p6 - fs_135_20449 * r_6 * h4_p2 + f_10_143 * r_8 * h2_p2);
    }

    // NOTE: the rows are formed in 53 loops, as the vectorizer runs out of
    // registers with all 66 of them in one.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, ph4_m4, ph4_p3, ph6_m4, ph6_p3, ph6_p5, ph8_m4, ph8_p3, ph8_p5, ph10_m4, ph10_p3, ph10_p5, ab_2, pc_14, pc_15 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h4_m4 = ph4_m4[k];
        const auto h4_p3 = ph4_p3[k];
        const auto h6_m4 = ph6_m4[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h6_p5 = ph6_p5[k];
        const auto h8_m4 = ph8_m4[k];
        const auto h8_p3 = ph8_p3[k];
        const auto h8_p5 = ph8_p5[k];
        const auto h10_m4 = ph10_m4[k];
        const auto h10_p3 = ph10_p3[k];
        const auto h10_p5 = ph10_p5[k];

        pc_14[k] = -fs_30375_32 * e_1 * h4_p3 + e_2 * (fs_1125_968 * h6_p3 - fs_6075_88 * h6_p5 + fs_273375_968 * r_2 * h4_p3) + e_3 * (fs_3920_1859 * h8_p3 + fs_147_143 * h8_p5 - fs_10_121 * r_2 * h6_p3 + fs_54_11 * r_2 * h6_p5 - fs_273375_40898 * r_4 * h4_p3) + e_4 * (fs_138915_29838094 * h10_p3 + fs_694575_29838094 * h10_p5 - fs_15680_671099 * r_2 * h8_p3 - fs_588_51623 * r_2 * h8_p5 + fs_10_34969 * r_4 * h6_p3 - fs_54_3179 * r_4 * h6_p5 + fs_270_20449 * r_6 * h4_p3);

        pc_15[k] = f_135_4 * e_1 * h4_m4 + e_2 * (fs_4500_121 * h6_m4 - f_405_22 * r_2 * h4_m4) + e_3 * (-fs_2695_676 * h8_m4 - fs_320_121 * r_2 * h6_m4 + f_405_143 * r_4 * h4_m4) + e_4 * (-fs_333396_14919047 * h10_m4 + fs_2695_61009 * r_2 * h8_m4 + fs_320_34969 * r_4 * h6_m4 - f_18_143 * r_6 * h4_m4);
    }

    // NOTE: the rows are formed in 53 loops, as the vectorizer runs out of
    // registers with all 66 of them in one.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, ph4_m3, ph6_m5, ph6_m3, ph8_m5, ph8_m3, ph10_m5, ph10_m3, ab_2, pc_16 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h4_m3 = ph4_m3[k];
        const auto h6_m5 = ph6_m5[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h8_m5 = ph8_m5[k];
        const auto h8_m3 = ph8_m3[k];
        const auto h10_m5 = ph10_m5[k];
        const auto h10_m3 = ph10_m3[k];

        pc_16[k] = -fs_30375_32 * e_1 * h4_m3 + e_2 * (fs_6075_88 * h6_m5 + fs_1125_968 * h6_m3 + fs_273375_968 * r_2 * h4_m3) + e_3 * (-fs_147_143 * h8_m5 + fs_3920_1859 * h8_m3 - fs_54_11 * r_2 * h6_m5 - fs_10_121 * r_2 * h6_m3 - fs_273375_40898 * r_4 * h4_m3) + e_4 * (-fs_694575_29838094 * h10_m5 + fs_138915_29838094 * h10_m3 + fs_588_51623 * r_2 * h8_m5 - fs_15680_671099 * r_2 * h8_m3 + fs_54_3179 * r_4 * h6_m5 + fs_10_34969 * r_4 * h6_m3 + fs_270_20449 * r_6 * h4_m3);
    }

    // NOTE: the rows are formed in 53 loops, as the vectorizer runs out of
    // registers with all 66 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_m2, ph4_m2, ph6_m6, ph6_m2, ph8_m6, ph8_m2, ph10_m6, ph10_m2, ab_2, pc_17 : simd::cache_line_size())
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

        const auto h2_m2 = ph2_m2[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h6_m6 = ph6_m6[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h8_m6 = ph8_m6[k];
        const auto h8_m2 = ph8_m2[k];
        const auto h10_m6 = ph10_m6[k];
        const auto h10_m2 = ph10_m2[k];

        pc_17[k] = f_315_8 * e_0 * h2_m2 + e_1 * (fs_30375_64 * h4_m2 - f_45 * r_2 * h2_m2) + e_2 * (fs_1575_22 * h6_m6 - fs_7875_242 * h6_m2 - fs_273375_1936 * r_2 * h4_m2 + f_15 * r_4 * h2_m2) + e_3 * (-fs_49_1144 * h8_m6 - fs_252105_163592 * h8_m2 - fs_56_11 * r_2 * h6_m6 + fs_280_121 * r_2 * h6_m2 + fs_273375_81796 * r_4 * h4_m2 - f_20_11 * r_6 * h2_m2) + e_4 * (-fs_635040_14919047 * h10_m6 - fs_317520_193947611 * h10_m2 + fs_49_103246 * r_2 * h8_m6 + fs_252105_14764178 * r_2 * h8_m2 + fs_56_3179 * r_4 * h6_m6 - fs_280_34969 * r_4 * h6_m2 - fs_135_20449 * r_6 * h4_m2 + f_10_143 * r_8 * h2_m2);
    }

    // NOTE: the rows are formed in 53 loops, as the vectorizer runs out of
    // registers with all 66 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_m1, ph4_m1, ph6_m1, ph8_m8, ph8_m7, ph8_m1, ph10_m9, ph10_m8, ph10_m7, ph10_m1, ab_2, pc_18, pc_19, pc_20 : simd::cache_line_size())
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
        const auto h4_m1 = ph4_m1[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h8_m8 = ph8_m8[k];
        const auto h8_m7 = ph8_m7[k];
        const auto h8_m1 = ph8_m1[k];
        const auto h10_m9 = ph10_m9[k];
        const auto h10_m8 = ph10_m8[k];
        const auto h10_m7 = ph10_m7[k];
        const auto h10_m1 = ph10_m1[k];

        pc_18[k] = -fs_1620675_512 * e_0 * h2_m1 + fs_33075_8 * e_1 * r_2 * h2_m1 + e_2 * (fs_525_8 * h6_m1 - fs_3675_8 * r_4 * h2_m1) + e_3 * (fs_245_286 * h8_m7 + fs_98_121 * h8_m1 - fs_14_3 * r_2 * h6_m1 + fs_2450_363 * r_6 * h2_m1) + e_4 * (-fs_59535_877591 * h10_m7 + fs_178605_387895222 * h10_m1 - fs_490_51623 * r_2 * h8_m7 - fs_392_43681 * r_2 * h8_m1 + fs_14_867 * r_4 * h6_m1 - fs_1225_122694 * r_8 * h2_m1);

        pc_19[k] = fs_2205_572 * e_3 * h8_m8 + e_4 * (-fs_79380_877591 * h10_m8 - fs_2205_51623 * r_2 * h8_m8);

        pc_20[k] = fs_1488375_512 * e_0 * h2_m1 + e_1 * (f_135_4 * h4_m1 - fs_30375_8 * r_2 * h2_m1) + e_2 * (fs_23625_968 * h6_m1 - f_405_22 * r_2 * h4_m1 + fs_3375_8 * r_4 * h2_m1) + e_3 * (fs_2205_40898 * h8_m1 - fs_210_121 * r_2 * h6_m1 + f_405_143 * r_4 * h4_m1 - fs_750_121 * r_6 * h2_m1) + e_4 * (-fs_3969_46189 * h10_m9 + fs_3969_387895222 * h10_m1 - fs_4410_7382089 * r_2 * h8_m1 + fs_210_34969 * r_4 * h6_m1 - f_18_143 * r_6 * h4_m1 + fs_375_40898 * r_8 * h2_m1);
    }

    // NOTE: the rows are formed in 53 loops, as the vectorizer runs out of
    // registers with all 66 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph2_0, ph4_0, ph6_0, ph6_p6, ph8_0, ph8_p6, ph10_0, ph10_p6, ab_2, pc_21 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;
        const auto r_10 = r_8 * r_2;

        const auto h2_0 = ph2_0[k];
        const auto h4_0 = ph4_0[k];
        const auto h6_0 = ph6_0[k];
        const auto h6_p6 = ph6_p6[k];
        const auto h8_0 = ph8_0[k];
        const auto h8_p6 = ph8_p6[k];
        const auto h10_0 = ph10_0[k];
        const auto h10_p6 = ph10_p6[k];

        pc_21[k] = e_0 * (-f_105_16 * h2_0 - f_1575_16 * r_2) + e_1 * (-f_135_4 * h4_0 + f_15_2 * r_2 * h2_0 + f_315_4 * r_4) + e_2 * (f_145_22 * h6_0 + fs_1050_11 * h6_p6 + f_405_22 * r_2 * h4_0 - f_5_2 * r_4 * h2_0 - f_45_2 * r_6) + e_3 * (f_511_286 * h8_0 - fs_294_143 * h8_p6 - f_58_33 * r_2 * h6_0 - fs_224_33 * r_2 * h6_p6 - f_405_143 * r_4 * h4_0 + f_10_33 * r_6 * h2_0 + f_5_2 * r_8) + e_4 * (f_2835_46189 * h10_0 + fs_1071630_14919047 * h10_p6 - f_511_2717 * r_2 * h8_0 + fs_1176_51623 * r_2 * h8_p6 + f_58_561 * r_4 * h6_0 + fs_224_9537 * r_4 * h6_p6 + f_18_143 * r_6 * h4_0 - f_5_429 * r_8 * h2_0 - f_1_11 * r_10) + f_945_32 * e_5;
    }

    // NOTE: the rows are formed in 53 loops, as the vectorizer runs out of
    // registers with all 66 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_p1, ph4_p1, ph6_p1, ph6_p5, ph8_p1, ph8_p5, ph10_p1, ph10_p5, ab_2, pc_22 : simd::cache_line_size())
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

        const auto h2_p1 = ph2_p1[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h6_p5 = ph6_p5[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h8_p5 = ph8_p5[k];
        const auto h10_p1 = ph10_p1[k];
        const auto h10_p5 = ph10_p5[k];

        pc_22[k] = -fs_275625_128 * e_0 * h2_p1 + e_1 * (fs_30375_64 * h4_p1 + fs_5625_2 * r_2 * h2_p1) + e_2 * (fs_175_242 * h6_p1 + fs_525_44 * h6_p5 - fs_273375_1936 * r_2 * h4_p1 - fs_625_2 * r_4 * h2_p1) + e_3 * (-fs_324723_163592 * h8_p1 - fs_1029_1144 * h8_p5 - fs_56_1089 * r_2 * h6_p1 - fs_28_33 * r_2 * h6_p5 + fs_273375_81796 * r_4 * h4_p1 + fs_5000_1089 * r_6 * h2_p1) + e_4 * (-fs_1071630_193947611 * h10_p1 + fs_893025_14919047 * h10_p5 + fs_324723_14764178 * r_2 * h8_p1 + fs_1029_103246 * r_2 * h8_p5 + fs_56_314721 * r_4 * h6_p1 + fs_28_9537 * r_4 * h6_p5 - fs_135_20449 * r_6 * h4_p1 - fs_1250_184041 * r_8 * h2_p1);
    }

    // NOTE: the rows are formed in 53 loops, as the vectorizer runs out of
    // registers with all 66 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_p2, ph4_m3, ph4_p4, ph6_m3, ph6_p2, ph6_p4, ph8_m3, ph8_p2, ph8_p4, ph10_m3, ph10_p2, ph10_p4, ab_2, pc_23, pc_24 : simd::cache_line_size())
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

        const auto h2_p2 = ph2_p2[k];
        const auto h4_m3 = ph4_m3[k];
        const auto h4_p4 = ph4_p4[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h6_p4 = ph6_p4[k];
        const auto h8_m3 = ph8_m3[k];
        const auto h8_p2 = ph8_p2[k];
        const auto h8_p4 = ph8_p4[k];
        const auto h10_m3 = ph10_m3[k];
        const auto h10_p2 = ph10_p2[k];
        const auto h10_p4 = ph10_p4[k];

        pc_23[k] = fs_77175_32 * e_0 * h2_p2 + e_1 * (fs_30375_32 * h4_p4 - fs_3150 * r_2 * h2_p2) + e_2 * (-fs_125_4 * h6_p2 - fs_12675_968 * h6_p4 - fs_273375_968 * r_2 * h4_p4 + fs_350 * r_4 * h2_p2) + e_3 * (fs_735_484 * h8_p2 - fs_147_14872 * h8_p4 + fs_20_9 * r_2 * h6_p2 + fs_338_363 * r_2 * h6_p4 + fs_273375_40898 * r_4 * h4_p4 - fs_5600_1089 * r_6 * h2_p2) + e_4 * (fs_2500470_193947611 * h10_p2 + fs_1250235_29838094 * h10_p4 - fs_735_43681 * r_2 * h8_p2 + fs_147_1342198 * r_2 * h8_p4 - fs_20_2601 * r_4 * h6_p2 - fs_338_104907 * r_4 * h6_p4 - fs_270_20449 * r_6 * h4_p4 + fs_1400_184041 * r_8 * h2_p2);

        pc_24[k] = -f_135_4 * e_1 * h4_m3 + e_2 * (fs_46875_484 * h6_m3 + f_405_22 * r_2 * h4_m3) + e_3 * (-fs_3675_3718 * h8_m3 - fs_2500_363 * r_2 * h6_m3 - f_405_143 * r_4 * h4_m3) + e_4 * (-fs_750141_14919047 * h10_m3 + fs_7350_671099 * r_2 * h8_m3 + fs_2500_104907 * r_4 * h6_m3 + f_18_143 * r_6 * h4_m3);
    }

    // NOTE: the rows are formed in 53 loops, as the vectorizer runs out of
    // registers with all 66 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_m2, ph4_m4, ph6_m4, ph6_m2, ph8_m4, ph8_m2, ph10_m4, ph10_m2, ab_2, pc_25 : simd::cache_line_size())
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

        const auto h2_m2 = ph2_m2[k];
        const auto h4_m4 = ph4_m4[k];
        const auto h6_m4 = ph6_m4[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h8_m4 = ph8_m4[k];
        const auto h8_m2 = ph8_m2[k];
        const auto h10_m4 = ph10_m4[k];
        const auto h10_m2 = ph10_m2[k];

        pc_25[k] = fs_77175_32 * e_0 * h2_m2 + e_1 * (-fs_30375_32 * h4_m4 - fs_3150 * r_2 * h2_m2) + e_2 * (fs_12675_968 * h6_m4 - fs_125_4 * h6_m2 + fs_273375_968 * r_2 * h4_m4 + fs_350 * r_4 * h2_m2) + e_3 * (fs_147_14872 * h8_m4 + fs_735_484 * h8_m2 - fs_338_363 * r_2 * h6_m4 + fs_20_9 * r_2 * h6_m2 - fs_273375_40898 * r_4 * h4_m4 - fs_5600_1089 * r_6 * h2_m2) + e_4 * (-fs_1250235_29838094 * h10_m4 + fs_2500470_193947611 * h10_m2 - fs_147_1342198 * r_2 * h8_m4 - fs_735_43681 * r_2 * h8_m2 + fs_338_104907 * r_4 * h6_m4 - fs_20_2601 * r_4 * h6_m2 + fs_270_20449 * r_6 * h4_m4 + fs_1400_184041 * r_8 * h2_m2);
    }

    // NOTE: the rows are formed in 53 loops, as the vectorizer runs out of
    // registers with all 66 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_m1, ph4_m1, ph6_m6, ph6_m5, ph6_m1, ph8_m6, ph8_m5, ph8_m1, ph10_m6, ph10_m5, ph10_m1, ab_2, pc_26, pc_27 : simd::cache_line_size())
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
        const auto h4_m1 = ph4_m1[k];
        const auto h6_m6 = ph6_m6[k];
        const auto h6_m5 = ph6_m5[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h8_m6 = ph8_m6[k];
        const auto h8_m5 = ph8_m5[k];
        const auto h8_m1 = ph8_m1[k];
        const auto h10_m6 = ph10_m6[k];
        const auto h10_m5 = ph10_m5[k];
        const auto h10_m1 = ph10_m1[k];

        pc_26[k] = -fs_275625_128 * e_0 * h2_m1 + e_1 * (fs_30375_64 * h4_m1 + fs_5625_2 * r_2 * h2_m1) + e_2 * (-fs_525_44 * h6_m5 + fs_175_242 * h6_m1 - fs_273375_1936 * r_2 * h4_m1 - fs_625_2 * r_4 * h2_m1) + e_3 * (fs_1029_1144 * h8_m5 - fs_324723_163592 * h8_m1 + fs_28_33 * r_2 * h6_m5 - fs_56_1089 * r_2 * h6_m1 + fs_273375_81796 * r_4 * h4_m1 + fs_5000_1089 * r_6 * h2_m1) + e_4 * (-fs_893025_14919047 * h10_m5 - fs_1071630_193947611 * h10_m1 - fs_1029_103246 * r_2 * h8_m5 + fs_324723_14764178 * r_2 * h8_m1 - fs_28_9537 * r_4 * h6_m5 + fs_56_314721 * r_4 * h6_m1 - fs_135_20449 * r_6 * h4_m1 - fs_1250_184041 * r_8 * h2_m1);

        pc_27[k] = -fs_1050_11 * e_2 * h6_m6 + e_3 * (fs_294_143 * h8_m6 + fs_224_33 * r_2 * h6_m6) + e_4 * (-fs_1071630_14919047 * h10_m6 - fs_1176_51623 * r_2 * h8_m6 - fs_224_9537 * r_4 * h6_m6);
    }

    // NOTE: the rows are formed in 53 loops, as the vectorizer runs out of
    // registers with all 66 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_m1, ph6_m1, ph8_m7, ph8_m1, ph10_m7, ph10_m1, ab_2, pc_28 : simd::cache_line_size())
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
        const auto h6_m1 = ph6_m1[k];
        const auto h8_m7 = ph8_m7[k];
        const auto h8_m1 = ph8_m1[k];
        const auto h10_m7 = ph10_m7[k];
        const auto h10_m1 = ph10_m1[k];

        pc_28[k] = fs_1620675_512 * e_0 * h2_m1 - fs_33075_8 * e_1 * r_2 * h2_m1 + e_2 * (-fs_525_8 * h6_m1 + fs_3675_8 * r_4 * h2_m1) + e_3 * (fs_245_286 * h8_m7 - fs_98_121 * h8_m1 + fs_14_3 * r_2 * h6_m1 - fs_2450_363 * r_6 * h2_m1) + e_4 * (-fs_59535_877591 * h10_m7 - fs_178605_387895222 * h10_m1 - fs_490_51623 * r_2 * h8_m7 + fs_392_43681 * r_2 * h8_m1 - fs_14_867 * r_4 * h6_m1 + fs_1225_122694 * r_8 * h2_m1);
    }

    // NOTE: the rows are formed in 53 loops, as the vectorizer runs out of
    // registers with all 66 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_m2, ph4_m2, ph6_m2, ph8_m8, ph8_m2, ph10_m8, ph10_m2, ab_2, pc_29 : simd::cache_line_size())
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

        const auto h2_m2 = ph2_m2[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h8_m8 = ph8_m8[k];
        const auto h8_m2 = ph8_m2[k];
        const auto h10_m8 = ph10_m8[k];
        const auto h10_m2 = ph10_m2[k];

        pc_29[k] = -fs_165375_256 * e_0 * h2_m2 + e_1 * (-f_135_4 * h4_m2 + fs_3375_4 * r_2 * h2_m2) + e_2 * (-fs_13125_242 * h6_m2 + f_405_22 * r_2 * h4_m2 - fs_375_4 * r_4 * h2_m2) + e_3 * (-fs_1225_572 * h8_m8 - fs_8575_40898 * h8_m2 + fs_1400_363 * r_2 * h6_m2 - f_405_143 * r_4 * h4_m2 + fs_500_363 * r_6 * h2_m2) + e_4 * (-fs_35721_877591 * h10_m8 - fs_11907_193947611 * h10_m2 + fs_1225_51623 * r_2 * h8_m8 + fs_17150_7382089 * r_2 * h8_m2 - fs_1400_104907 * r_4 * h6_m2 + f_18_143 * r_6 * h4_m2 - fs_125_61347 * r_8 * h2_m2);
    }

    // NOTE: the rows are formed in 53 loops, as the vectorizer runs out of
    // registers with all 66 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph2_0, ph4_0, ph4_p4, ph6_0, ph6_p4, ph8_0, ph8_p4, ph10_0, ph10_p4, ab_2, pc_30 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;
        const auto r_10 = r_8 * r_2;

        const auto h2_0 = ph2_0[k];
        const auto h4_0 = ph4_0[k];
        const auto h4_p4 = ph4_p4[k];
        const auto h6_0 = ph6_0[k];
        const auto h6_p4 = ph6_p4[k];
        const auto h8_0 = ph8_0[k];
        const auto h8_p4 = ph8_p4[k];
        const auto h10_0 = ph10_0[k];
        const auto h10_p4 = ph10_p4[k];

        pc_30[k] = e_0 * (-f_315_8 * h2_0 - f_1575_16 * r_2) + e_1 * (-f_45_8 * h4_0 - fs_70875_64 * h4_p4 + f_45 * r_2 * h2_0 + f_315_4 * r_4) + e_2 * (f_90_11 * h6_0 + fs_6300_121 * h6_p4 + f_135_44 * r_2 * h4_0 + fs_637875_1936 * r_2 * h4_p4 - f_15 * r_4 * h2_0 - f_45_2 * r_6) + e_3 * (-f_238_143 * h8_0 - fs_3087_1859 * h8_p4 - f_24_11 * r_2 * h6_0 - fs_448_121 * r_2 * h6_p4 - f_135_286 * r_4 * h4_0 - fs_637875_81796 * r_4 * h4_p4 + f_20_11 * r_6 * h2_0 + f_5_2 * r_8) + e_4 * (-f_7560_46189 * h10_0 + fs_952560_14919047 * h10_p4 + f_476_2717 * r_2 * h8_0 + fs_12348_671099 * r_2 * h8_p4 + f_24_187 * r_4 * h6_0 + fs_448_34969 * r_4 * h6_p4 + f_3_143 * r_6 * h4_0 + fs_315_20449 * r_6 * h4_p4 - f_10_143 * r_8 * h2_0 - f_1_11 * r_10) + f_945_32 * e_5;
    }

    // NOTE: the rows are formed in 53 loops, as the vectorizer runs out of
    // registers with all 66 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_p1, ph4_p1, ph4_p3, ph6_p1, ph6_p3, ph8_p1, ph8_p3, ph10_p1, ph10_p3, ab_2, pc_31 : simd::cache_line_size())
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

        const auto h2_p1 = ph2_p1[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h4_p3 = ph4_p3[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h8_p3 = ph8_p3[k];
        const auto h10_p1 = ph10_p1[k];
        const auto h10_p3 = ph10_p3[k];

        pc_31[k] = -fs_231525_256 * e_0 * h2_p1 + e_1 * (fs_70875_128 * h4_p1 - fs_10125_128 * h4_p3 + fs_4725_4 * r_2 * h2_p1) + e_2 * (-fs_4800_121 * h6_p1 + fs_3375_242 * h6_p3 - fs_637875_3872 * r_2 * h4_p1 + fs_91125_3872 * r_2 * h4_p3 - fs_525_4 * r_4 * h2_p1) + e_3 * (fs_27783_81796 * h8_p1 - fs_6615_7436 * h8_p3 + fs_1024_363 * r_2 * h6_p1 - fs_120_121 * r_2 * h6_p3 + fs_637875_163592 * r_4 * h4_p1 - fs_91125_163592 * r_4 * h4_p3 + fs_700_363 * r_6 * h2_p1) + e_4 * (fs_5000940_193947611 * h10_p1 + fs_833490_14919047 * h10_p3 - fs_27783_7382089 * r_2 * h8_p1 + fs_6615_671099 * r_2 * h8_p3 - fs_1024_104907 * r_4 * h6_p1 + fs_120_34969 * r_4 * h6_p3 - fs_315_40898 * r_6 * h4_p1 + fs_45_40898 * r_6 * h4_p3 - fs_175_61347 * r_8 * h2_p1);
    }

    // NOTE: the rows are formed in 53 loops, as the vectorizer runs out of
    // registers with all 66 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_m2, ph4_m2, ph6_m2, ph8_m2, ph10_m2, ab_2, pc_32 : simd::cache_line_size())
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

        const auto h2_m2 = ph2_m2[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h8_m2 = ph8_m2[k];
        const auto h10_m2 = ph10_m2[k];

        pc_32[k] = fs_385875_64 * e_0 * h2_m2 + e_1 * (-fs_42525_64 * h4_m2 - fs_7875 * r_2 * h2_m2) + e_2 * (fs_1250_121 * h6_m2 + fs_382725_1936 * r_2 * h4_m2 + fs_875 * r_4 * h2_m2) + e_3 * (fs_3675_40898 * h8_m2 - fs_800_1089 * r_2 * h6_m2 - fs_382725_81796 * r_4 * h4_m2 - fs_14000_1089 * r_6 * h2_m2) + e_4 * (-fs_16003008_193947611 * h10_m2 - fs_7350_7382089 * r_2 * h8_m2 + fs_800_314721 * r_4 * h6_m2 + fs_189_20449 * r_6 * h4_m2 + fs_3500_184041 * r_8 * h2_m2);
    }

    // NOTE: the rows are formed in 53 loops, as the vectorizer runs out of
    // registers with all 66 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_m1, ph4_m3, ph4_m1, ph6_m3, ph6_m1, ph8_m3, ph8_m1, ph10_m3, ph10_m1, ab_2, pc_33 : simd::cache_line_size())
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
        const auto h4_m3 = ph4_m3[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h8_m3 = ph8_m3[k];
        const auto h8_m1 = ph8_m1[k];
        const auto h10_m3 = ph10_m3[k];
        const auto h10_m1 = ph10_m1[k];

        pc_33[k] = -fs_231525_256 * e_0 * h2_m1 + e_1 * (fs_10125_128 * h4_m3 + fs_70875_128 * h4_m1 + fs_4725_4 * r_2 * h2_m1) + e_2 * (-fs_3375_242 * h6_m3 - fs_4800_121 * h6_m1 - fs_91125_3872 * r_2 * h4_m3 - fs_637875_3872 * r_2 * h4_m1 - fs_525_4 * r_4 * h2_m1) + e_3 * (fs_6615_7436 * h8_m3 + fs_27783_81796 * h8_m1 + fs_120_121 * r_2 * h6_m3 + fs_1024_363 * r_2 * h6_m1 + fs_91125_163592 * r_4 * h4_m3 + fs_637875_163592 * r_4 * h4_m1 + fs_700_363 * r_6 * h2_m1) + e_4 * (-fs_833490_14919047 * h10_m3 + fs_5000940_193947611 * h10_m1 - fs_6615_671099 * r_2 * h8_m3 - fs_27783_7382089 * r_2 * h8_m1 - fs_120_34969 * r_4 * h6_m3 - fs_1024_104907 * r_4 * h6_m1 - fs_45_40898 * r_6 * h4_m3 - fs_315_40898 * r_6 * h4_m1 - fs_175_61347 * r_8 * h2_m1);
    }

    // NOTE: the rows are formed in 53 loops, as the vectorizer runs out of
    // registers with all 66 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_m1, ph4_m4, ph4_m1, ph6_m5, ph6_m4, ph6_m1, ph8_m5, ph8_m4, ph8_m1, ph10_m5, ph10_m4, ph10_m1, ab_2, pc_34, pc_35 : simd::cache_line_size())
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
        const auto h4_m4 = ph4_m4[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h6_m5 = ph6_m5[k];
        const auto h6_m4 = ph6_m4[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h8_m5 = ph8_m5[k];
        const auto h8_m4 = ph8_m4[k];
        const auto h8_m1 = ph8_m1[k];
        const auto h10_m5 = ph10_m5[k];
        const auto h10_m4 = ph10_m4[k];
        const auto h10_m1 = ph10_m1[k];

        pc_34[k] = fs_70875_64 * e_1 * h4_m4 + e_2 * (-fs_6300_121 * h6_m4 - fs_637875_1936 * r_2 * h4_m4) + e_3 * (fs_3087_1859 * h8_m4 + fs_448_121 * r_2 * h6_m4 + fs_637875_81796 * r_4 * h4_m4) + e_4 * (-fs_952560_14919047 * h10_m4 - fs_12348_671099 * r_2 * h8_m4 - fs_448_34969 * r_4 * h6_m4 - fs_315_20449 * r_6 * h4_m4);

        pc_35[k] = fs_275625_128 * e_0 * h2_m1 + e_1 * (-fs_30375_64 * h4_m1 - fs_5625_2 * r_2 * h2_m1) + e_2 * (-fs_525_44 * h6_m5 - fs_175_242 * h6_m1 + fs_273375_1936 * r_2 * h4_m1 + fs_625_2 * r_4 * h2_m1) + e_3 * (fs_1029_1144 * h8_m5 + fs_324723_163592 * h8_m1 + fs_28_33 * r_2 * h6_m5 + fs_56_1089 * r_2 * h6_m1 - fs_273375_81796 * r_4 * h4_m1 - fs_5000_1089 * r_6 * h2_m1) + e_4 * (-fs_893025_14919047 * h10_m5 + fs_1071630_193947611 * h10_m1 - fs_1029_103246 * r_2 * h8_m5 - fs_324723_14764178 * r_2 * h8_m1 - fs_28_9537 * r_4 * h6_m5 - fs_56_314721 * r_4 * h6_m1 + fs_135_20449 * r_6 * h4_m1 + fs_1250_184041 * r_8 * h2_m1);
    }

    // NOTE: the rows are formed in 53 loops, as the vectorizer runs out of
    // registers with all 66 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_m2, ph4_m2, ph6_m6, ph6_m2, ph8_m6, ph8_m2, ph10_m6, ph10_m2, ab_2, pc_36 : simd::cache_line_size())
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

        const auto h2_m2 = ph2_m2[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h6_m6 = ph6_m6[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h8_m6 = ph8_m6[k];
        const auto h8_m2 = ph8_m2[k];
        const auto h10_m6 = ph10_m6[k];
        const auto h10_m2 = ph10_m2[k];

        pc_36[k] = -f_315_8 * e_0 * h2_m2 + e_1 * (-fs_30375_64 * h4_m2 + f_45 * r_2 * h2_m2) + e_2 * (fs_1575_22 * h6_m6 + fs_7875_242 * h6_m2 + fs_273375_1936 * r_2 * h4_m2 - f_15 * r_4 * h2_m2) + e_3 * (-fs_49_1144 * h8_m6 + fs_252105_163592 * h8_m2 - fs_56_11 * r_2 * h6_m6 - fs_280_121 * r_2 * h6_m2 - fs_273375_81796 * r_4 * h4_m2 + f_20_11 * r_6 * h2_m2) + e_4 * (-fs_635040_14919047 * h10_m6 + fs_317520_193947611 * h10_m2 + fs_49_103246 * r_2 * h8_m6 - fs_252105_14764178 * r_2 * h8_m2 + fs_56_3179 * r_4 * h6_m6 + fs_280_34969 * r_4 * h6_m2 + fs_135_20449 * r_6 * h4_m2 - f_10_143 * r_8 * h2_m2);
    }

    // NOTE: the rows are formed in 53 loops, as the vectorizer runs out of
    // registers with all 66 of them in one.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, ph4_m3, ph6_m3, ph8_m7, ph8_m3, ph10_m7, ph10_m3, ab_2, pc_37 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h4_m3 = ph4_m3[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h8_m7 = ph8_m7[k];
        const auto h8_m3 = ph8_m3[k];
        const auto h10_m7 = ph10_m7[k];
        const auto h10_m3 = ph10_m3[k];

        pc_37[k] = fs_42525_64 * e_1 * h4_m3 + e_2 * (fs_39375_484 * h6_m3 - fs_382725_1936 * r_2 * h4_m3) + e_3 * (-fs_3675_1144 * h8_m7 + fs_8575_14872 * h8_m3 - fs_700_121 * r_2 * h6_m3 + fs_382725_81796 * r_4 * h4_m3) + e_4 * (-fs_15876_877591 * h10_m7 + fs_3969_14919047 * h10_m3 + fs_3675_103246 * r_2 * h8_m7 - fs_8575_1342198 * r_2 * h8_m3 + fs_700_34969 * r_4 * h6_m3 - fs_189_20449 * r_6 * h4_m3);
    }

    // NOTE: the rows are formed in 53 loops, as the vectorizer runs out of
    // registers with all 66 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph2_0, ph2_p2, ph4_0, ph4_p2, ph6_0, ph6_p2, ph8_0, ph8_p2, ph10_0, ph10_p2, ab_2, pc_38 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;
        const auto r_10 = r_8 * r_2;

        const auto h2_0 = ph2_0[k];
        const auto h2_p2 = ph2_p2[k];
        const auto h4_0 = ph4_0[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h6_0 = ph6_0[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h8_0 = ph8_0[k];
        const auto h8_p2 = ph8_p2[k];
        const auto h10_0 = ph10_0[k];
        const auto h10_p2 = ph10_p2[k];

        pc_38[k] = e_0 * (-f_945_16 * h2_0 + fs_826875_256 * h2_p2 - f_1575_16 * r_2) + e_1 * (f_45_2 * h4_0 - fs_10125_16 * h4_p2 + f_135_2 * r_2 * h2_0 - fs_16875_4 * r_2 * h2_p2 + f_315_4 * r_4) + e_2 * (-f_30_11 * h6_0 + fs_5250_121 * h6_p2 - f_135_11 * r_2 * h4_0 + fs_91125_484 * r_2 * h4_p2 - f_45_2 * r_4 * h2_0 + fs_1875_4 * r_4 * h2_p2 - f_45_2 * r_6) + e_3 * (-f_49_143 * h8_0 - fs_30870_20449 * h8_p2 + f_8_11 * r_2 * h6_0 - fs_1120_363 * r_2 * h6_p2 + f_270_143 * r_4 * h4_0 - fs_91125_20449 * r_4 * h4_p2 + f_30_11 * r_6 * h2_0 - fs_2500_363 * r_6 * h2_p2 + f_5_2 * r_8) + e_4 * (f_13230_46189 * h10_0 + fs_11668860_193947611 * h10_p2 + f_98_2717 * r_2 * h8_0 + fs_123480_7382089 * r_2 * h8_p2 - f_8_187 * r_4 * h6_0 + fs_1120_104907 * r_4 * h6_p2 - f_12_143 * r_6 * h4_0 + fs_180_20449 * r_6 * h4_p2 - f_15_143 * r_8 * h2_0 + fs_625_61347 * r_8 * h2_p2 - f_1_11 * r_10) + f_945_32 * e_5;
    }

    // NOTE: the rows are formed in 53 loops, as the vectorizer runs out of
    // registers with all 66 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_m2, ph2_m1, ph4_m2, ph4_m1, ph6_m2, ph6_m1, ph8_m2, ph8_m1, ph10_m2, ph10_m1, ab_2, pc_39, pc_40 : simd::cache_line_size())
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

        const auto h2_m2 = ph2_m2[k];
        const auto h2_m1 = ph2_m1[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h8_m2 = ph8_m2[k];
        const auto h8_m1 = ph8_m1[k];
        const auto h10_m2 = ph10_m2[k];
        const auto h10_m1 = ph10_m1[k];

        pc_39[k] = -fs_55125_256 * e_0 * h2_m1 + e_1 * (fs_6075_32 * h4_m1 + fs_1125_4 * r_2 * h2_m1) + e_2 * (-fs_3500_121 * h6_m1 - fs_54675_968 * r_2 * h4_m1 - fs_125_4 * r_4 * h2_m1) + e_3 * (fs_36015_20449 * h8_m1 + fs_2240_1089 * r_2 * h6_m1 + fs_54675_40898 * r_4 * h4_m1 + fs_500_1089 * r_6 * h2_m1) + e_4 * (-fs_21003948_193947611 * h10_m1 - fs_144060_7382089 * r_2 * h8_m1 - fs_2240_314721 * r_4 * h6_m1 - fs_54_20449 * r_6 * h4_m1 - fs_125_184041 * r_8 * h2_m1);

        pc_40[k] = -fs_826875_256 * e_0 * h2_m2 + e_1 * (fs_10125_16 * h4_m2 + fs_16875_4 * r_2 * h2_m2) + e_2 * (-fs_5250_121 * h6_m2 - fs_91125_484 * r_2 * h4_m2 - fs_1875_4 * r_4 * h2_m2) + e_3 * (fs_30870_20449 * h8_m2 + fs_1120_363 * r_2 * h6_m2 + fs_91125_20449 * r_4 * h4_m2 + fs_2500_363 * r_6 * h2_m2) + e_4 * (-fs_11668860_193947611 * h10_m2 - fs_123480_7382089 * r_2 * h8_m2 - fs_1120_104907 * r_4 * h6_m2 - fs_180_20449 * r_6 * h4_m2 - fs_625_61347 * r_8 * h2_m2);
    }

    // NOTE: the rows are formed in 53 loops, as the vectorizer runs out of
    // registers with all 66 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_m1, ph4_m3, ph4_m1, ph6_m3, ph6_m1, ph8_m3, ph8_m1, ph10_m3, ph10_m1, ab_2, pc_41 : simd::cache_line_size())
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
        const auto h4_m3 = ph4_m3[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h8_m3 = ph8_m3[k];
        const auto h8_m1 = ph8_m1[k];
        const auto h10_m3 = ph10_m3[k];
        const auto h10_m1 = ph10_m1[k];

        pc_41[k] = fs_231525_256 * e_0 * h2_m1 + e_1 * (fs_10125_128 * h4_m3 - fs_70875_128 * h4_m1 - fs_4725_4 * r_2 * h2_m1) + e_2 * (-fs_3375_242 * h6_m3 + fs_4800_121 * h6_m1 - fs_91125_3872 * r_2 * h4_m3 + fs_637875_3872 * r_2 * h4_m1 + fs_525_4 * r_4 * h2_m1) + e_3 * (fs_6615_7436 * h8_m3 - fs_27783_81796 * h8_m1 + fs_120_121 * r_2 * h6_m3 - fs_1024_363 * r_2 * h6_m1 + fs_91125_163592 * r_4 * h4_m3 - fs_637875_163592 * r_4 * h4_m1 - fs_700_363 * r_6 * h2_m1) + e_4 * (-fs_833490_14919047 * h10_m3 - fs_5000940_193947611 * h10_m1 - fs_6615_671099 * r_2 * h8_m3 + fs_27783_7382089 * r_2 * h8_m1 - fs_120_34969 * r_4 * h6_m3 + fs_1024_104907 * r_4 * h6_m1 - fs_45_40898 * r_6 * h4_m3 + fs_315_40898 * r_6 * h4_m1 + fs_175_61347 * r_8 * h2_m1);
    }

    // NOTE: the rows are formed in 53 loops, as the vectorizer runs out of
    // registers with all 66 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_m2, ph4_m4, ph6_m4, ph6_m2, ph8_m4, ph8_m2, ph10_m4, ph10_m2, ab_2, pc_42 : simd::cache_line_size())
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

        const auto h2_m2 = ph2_m2[k];
        const auto h4_m4 = ph4_m4[k];
        const auto h6_m4 = ph6_m4[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h8_m4 = ph8_m4[k];
        const auto h8_m2 = ph8_m2[k];
        const auto h10_m4 = ph10_m4[k];
        const auto h10_m2 = ph10_m2[k];

        pc_42[k] = -fs_77175_32 * e_0 * h2_m2 + e_1 * (-fs_30375_32 * h4_m4 + fs_3150 * r_2 * h2_m2) + e_2 * (fs_12675_968 * h6_m4 + fs_125_4 * h6_m2 + fs_273375_968 * r_2 * h4_m4 - fs_350 * r_4 * h2_m2) + e_3 * (fs_147_14872 * h8_m4 - fs_735_484 * h8_m2 - fs_338_363 * r_2 * h6_m4 - fs_20_9 * r_2 * h6_m2 - fs_273375_40898 * r_4 * h4_m4 + fs_5600_1089 * r_6 * h2_m2) + e_4 * (-fs_1250235_29838094 * h10_m4 - fs_2500470_193947611 * h10_m2 - fs_147_1342198 * r_2 * h8_m4 + fs_735_43681 * r_2 * h8_m2 + fs_338_104907 * r_4 * h6_m4 + fs_20_2601 * r_4 * h6_m2 + fs_270_20449 * r_6 * h4_m4 - fs_1400_184041 * r_8 * h2_m2);
    }

    // NOTE: the rows are formed in 53 loops, as the vectorizer runs out of
    // registers with all 66 of them in one.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, ph4_m3, ph6_m5, ph6_m3, ph8_m5, ph8_m3, ph10_m5, ph10_m3, ab_2, pc_43 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h4_m3 = ph4_m3[k];
        const auto h6_m5 = ph6_m5[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h8_m5 = ph8_m5[k];
        const auto h8_m3 = ph8_m3[k];
        const auto h10_m5 = ph10_m5[k];
        const auto h10_m3 = ph10_m3[k];

        pc_43[k] = fs_30375_32 * e_1 * h4_m3 + e_2 * (fs_6075_88 * h6_m5 - fs_1125_968 * h6_m3 - fs_273375_968 * r_2 * h4_m3) + e_3 * (-fs_147_143 * h8_m5 - fs_3920_1859 * h8_m3 - fs_54_11 * r_2 * h6_m5 + fs_10_121 * r_2 * h6_m3 + fs_273375_40898 * r_4 * h4_m3) + e_4 * (-fs_694575_29838094 * h10_m5 - fs_138915_29838094 * h10_m3 + fs_588_51623 * r_2 * h8_m5 + fs_15680_671099 * r_2 * h8_m3 + fs_54_3179 * r_4 * h6_m5 - fs_10_34969 * r_4 * h6_m3 - fs_270_20449 * r_6 * h4_m3);
    }

    // NOTE: the rows are formed in 53 loops, as the vectorizer runs out of
    // registers with all 66 of them in one.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, ph4_m4, ph6_m6, ph6_m4, ph8_m6, ph8_m4, ph10_m6, ph10_m4, ab_2, pc_44 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h4_m4 = ph4_m4[k];
        const auto h6_m6 = ph6_m6[k];
        const auto h6_m4 = ph6_m4[k];
        const auto h8_m6 = ph8_m6[k];
        const auto h8_m4 = ph8_m4[k];
        const auto h10_m6 = ph10_m6[k];
        const auto h10_m4 = ph10_m4[k];

        pc_44[k] = -fs_6075_32 * e_1 * h4_m4 + e_2 * (-fs_1125_44 * h6_m6 - fs_84375_968 * h6_m4 + fs_54675_968 * r_2 * h4_m4) + e_3 * (-fs_1715_572 * h8_m6 - fs_18375_14872 * h8_m4 + fs_20_11 * r_2 * h6_m6 + fs_750_121 * r_2 * h6_m4 - fs_54675_40898 * r_4 * h4_m4) + e_4 * (-fs_111132_14919047 * h10_m6 - fs_27783_29838094 * h10_m4 + fs_1715_51623 * r_2 * h8_m6 + fs_18375_1342198 * r_2 * h8_m4 - fs_20_3179 * r_4 * h6_m6 - fs_750_34969 * r_4 * h6_m4 + fs_54_20449 * r_6 * h4_m4);
    }

    // NOTE: the rows are formed in 53 loops, as the vectorizer runs out of
    // registers with all 66 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph2_0, ph2_p1, ph4_0, ph4_p1, ph6_0, ph6_p1, ph8_0, ph8_p1, ph10_0, ph10_p1, ab_2, pc_45, pc_46 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;
        const auto r_10 = r_8 * r_2;

        const auto h2_0 = ph2_0[k];
        const auto h2_p1 = ph2_p1[k];
        const auto h4_0 = ph4_0[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h6_0 = ph6_0[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h8_0 = ph8_0[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h10_0 = ph10_0[k];
        const auto h10_p1 = ph10_p1[k];

        pc_45[k] = e_0 * (-f_525_8 * h2_0 - f_1575_16 * r_2) + e_1 * (f_135_4 * h4_0 + f_75 * r_2 * h2_0 + f_315_4 * r_4) + e_2 * (-f_100_11 * h6_0 - f_405_22 * r_2 * h4_0 - f_25 * r_4 * h2_0 - f_45_2 * r_6) + e_3 * (f_245_143 * h8_0 + f_80_33 * r_2 * h6_0 + f_405_143 * r_4 * h4_0 + f_100_33 * r_6 * h2_0 + f_5_2 * r_8) + e_4 * (-f_15876_46189 * h10_0 - f_490_2717 * r_2 * h8_0 - f_80_561 * r_4 * h6_0 - f_18_143 * r_6 * h4_0 - f_50_429 * r_8 * h2_0 - f_1_11 * r_10) + f_945_32 * e_5;

        pc_46[k] = -fs_55125_256 * e_0 * h2_p1 + e_1 * (fs_6075_32 * h4_p1 + fs_1125_4 * r_2 * h2_p1) + e_2 * (-fs_3500_121 * h6_p1 - fs_54675_968 * r_2 * h4_p1 - fs_125_4 * r_4 * h2_p1) + e_3 * (fs_36015_20449 * h8_p1 + fs_2240_1089 * r_2 * h6_p1 + fs_54675_40898 * r_4 * h4_p1 + fs_500_1089 * r_6 * h2_p1) + e_4 * (-fs_21003948_193947611 * h10_p1 - fs_144060_7382089 * r_2 * h8_p1 - fs_2240_314721 * r_4 * h6_p1 - fs_54_20449 * r_6 * h4_p1 - fs_125_184041 * r_8 * h2_p1);
    }

    // NOTE: the rows are formed in 53 loops, as the vectorizer runs out of
    // registers with all 66 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_p2, ph4_p2, ph4_p3, ph6_p2, ph6_p3, ph8_p2, ph8_p3, ph10_p2, ph10_p3, ab_2, pc_47, pc_48 : simd::cache_line_size())
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

        const auto h2_p2 = ph2_p2[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h4_p3 = ph4_p3[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h8_p2 = ph8_p2[k];
        const auto h8_p3 = ph8_p3[k];
        const auto h10_p2 = ph10_p2[k];
        const auto h10_p3 = ph10_p3[k];

        pc_47[k] = fs_385875_64 * e_0 * h2_p2 + e_1 * (-fs_42525_64 * h4_p2 - fs_7875 * r_2 * h2_p2) + e_2 * (fs_1250_121 * h6_p2 + fs_382725_1936 * r_2 * h4_p2 + fs_875 * r_4 * h2_p2) + e_3 * (fs_3675_40898 * h8_p2 - fs_800_1089 * r_2 * h6_p2 - fs_382725_81796 * r_4 * h4_p2 - fs_14000_1089 * r_6 * h2_p2) + e_4 * (-fs_16003008_193947611 * h10_p2 - fs_7350_7382089 * r_2 * h8_p2 + fs_800_314721 * r_4 * h6_p2 + fs_189_20449 * r_6 * h4_p2 + fs_3500_184041 * r_8 * h2_p2);

        pc_48[k] = -f_135_4 * e_1 * h4_p3 + e_2 * (fs_46875_484 * h6_p3 + f_405_22 * r_2 * h4_p3) + e_3 * (-fs_3675_3718 * h8_p3 - fs_2500_363 * r_2 * h6_p3 - f_405_143 * r_4 * h4_p3) + e_4 * (-fs_750141_14919047 * h10_p3 + fs_7350_671099 * r_2 * h8_p3 + fs_2500_104907 * r_4 * h6_p3 + f_18_143 * r_6 * h4_p3);
    }

    // NOTE: the rows are formed in 53 loops, as the vectorizer runs out of
    // registers with all 66 of them in one.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, ph4_p4, ph6_p4, ph6_p5, ph8_p4, ph8_p5, ph10_p4, ph10_p5, ab_2, pc_49, pc_50 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h4_p4 = ph4_p4[k];
        const auto h6_p4 = ph6_p4[k];
        const auto h6_p5 = ph6_p5[k];
        const auto h8_p4 = ph8_p4[k];
        const auto h8_p5 = ph8_p5[k];
        const auto h10_p4 = ph10_p4[k];
        const auto h10_p5 = ph10_p5[k];

        pc_49[k] = f_135_4 * e_1 * h4_p4 + e_2 * (fs_4500_121 * h6_p4 - f_405_22 * r_2 * h4_p4) + e_3 * (-fs_2695_676 * h8_p4 - fs_320_121 * r_2 * h6_p4 + f_405_143 * r_4 * h4_p4) + e_4 * (-fs_333396_14919047 * h10_p4 + fs_2695_61009 * r_2 * h8_p4 + fs_320_34969 * r_4 * h6_p4 - f_18_143 * r_6 * h4_p4);

        pc_50[k] = -fs_5625_44 * e_2 * h6_p5 + e_3 * (-fs_1225_286 * h8_p5 + fs_100_11 * r_2 * h6_p5) + e_4 * (-fs_83349_14919047 * h10_p5 + fs_2450_51623 * r_2 * h8_p5 - fs_100_3179 * r_4 * h6_p5);
    }

    // NOTE: the rows are formed in 53 loops, as the vectorizer runs out of
    // registers with all 66 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph2_0, ph2_p2, ph4_0, ph4_p2, ph6_0, ph6_p2, ph8_0, ph8_p2, ph10_0, ph10_p2, ab_2, pc_51 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;
        const auto r_10 = r_8 * r_2;

        const auto h2_0 = ph2_0[k];
        const auto h2_p2 = ph2_p2[k];
        const auto h4_0 = ph4_0[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h6_0 = ph6_0[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h8_0 = ph8_0[k];
        const auto h8_p2 = ph8_p2[k];
        const auto h10_0 = ph10_0[k];
        const auto h10_p2 = ph10_p2[k];

        pc_51[k] = e_0 * (-f_945_16 * h2_0 - fs_826875_256 * h2_p2 - f_1575_16 * r_2) + e_1 * (f_45_2 * h4_0 + fs_10125_16 * h4_p2 + f_135_2 * r_2 * h2_0 + fs_16875_4 * r_2 * h2_p2 + f_315_4 * r_4) + e_2 * (-f_30_11 * h6_0 - fs_5250_121 * h6_p2 - f_135_11 * r_2 * h4_0 - fs_91125_484 * r_2 * h4_p2 - f_45_2 * r_4 * h2_0 - fs_1875_4 * r_4 * h2_p2 - f_45_2 * r_6) + e_3 * (-f_49_143 * h8_0 + fs_30870_20449 * h8_p2 + f_8_11 * r_2 * h6_0 + fs_1120_363 * r_2 * h6_p2 + f_270_143 * r_4 * h4_0 + fs_91125_20449 * r_4 * h4_p2 + f_30_11 * r_6 * h2_0 + fs_2500_363 * r_6 * h2_p2 + f_5_2 * r_8) + e_4 * (f_13230_46189 * h10_0 - fs_11668860_193947611 * h10_p2 + f_98_2717 * r_2 * h8_0 - fs_123480_7382089 * r_2 * h8_p2 - f_8_187 * r_4 * h6_0 - fs_1120_104907 * r_4 * h6_p2 - f_12_143 * r_6 * h4_0 - fs_180_20449 * r_6 * h4_p2 - f_15_143 * r_8 * h2_0 - fs_625_61347 * r_8 * h2_p2 - f_1_11 * r_10) + f_945_32 * e_5;
    }

    // NOTE: the rows are formed in 53 loops, as the vectorizer runs out of
    // registers with all 66 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_p1, ph4_p1, ph4_p3, ph6_p1, ph6_p3, ph8_p1, ph8_p3, ph10_p1, ph10_p3, ab_2, pc_52 : simd::cache_line_size())
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

        const auto h2_p1 = ph2_p1[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h4_p3 = ph4_p3[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h8_p3 = ph8_p3[k];
        const auto h10_p1 = ph10_p1[k];
        const auto h10_p3 = ph10_p3[k];

        pc_52[k] = -fs_231525_256 * e_0 * h2_p1 + e_1 * (fs_70875_128 * h4_p1 + fs_10125_128 * h4_p3 + fs_4725_4 * r_2 * h2_p1) + e_2 * (-fs_4800_121 * h6_p1 - fs_3375_242 * h6_p3 - fs_637875_3872 * r_2 * h4_p1 - fs_91125_3872 * r_2 * h4_p3 - fs_525_4 * r_4 * h2_p1) + e_3 * (fs_27783_81796 * h8_p1 + fs_6615_7436 * h8_p3 + fs_1024_363 * r_2 * h6_p1 + fs_120_121 * r_2 * h6_p3 + fs_637875_163592 * r_4 * h4_p1 + fs_91125_163592 * r_4 * h4_p3 + fs_700_363 * r_6 * h2_p1) + e_4 * (fs_5000940_193947611 * h10_p1 - fs_833490_14919047 * h10_p3 - fs_27783_7382089 * r_2 * h8_p1 - fs_6615_671099 * r_2 * h8_p3 - fs_1024_104907 * r_4 * h6_p1 - fs_120_34969 * r_4 * h6_p3 - fs_315_40898 * r_6 * h4_p1 - fs_45_40898 * r_6 * h4_p3 - fs_175_61347 * r_8 * h2_p1);
    }

    // NOTE: the rows are formed in 53 loops, as the vectorizer runs out of
    // registers with all 66 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_p2, ph4_p4, ph6_p2, ph6_p4, ph8_p2, ph8_p4, ph10_p2, ph10_p4, ab_2, pc_53 : simd::cache_line_size())
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

        const auto h2_p2 = ph2_p2[k];
        const auto h4_p4 = ph4_p4[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h6_p4 = ph6_p4[k];
        const auto h8_p2 = ph8_p2[k];
        const auto h8_p4 = ph8_p4[k];
        const auto h10_p2 = ph10_p2[k];
        const auto h10_p4 = ph10_p4[k];

        pc_53[k] = fs_77175_32 * e_0 * h2_p2 + e_1 * (-fs_30375_32 * h4_p4 - fs_3150 * r_2 * h2_p2) + e_2 * (-fs_125_4 * h6_p2 + fs_12675_968 * h6_p4 + fs_273375_968 * r_2 * h4_p4 + fs_350 * r_4 * h2_p2) + e_3 * (fs_735_484 * h8_p2 + fs_147_14872 * h8_p4 + fs_20_9 * r_2 * h6_p2 - fs_338_363 * r_2 * h6_p4 - fs_273375_40898 * r_4 * h4_p4 - fs_5600_1089 * r_6 * h2_p2) + e_4 * (fs_2500470_193947611 * h10_p2 - fs_1250235_29838094 * h10_p4 - fs_735_43681 * r_2 * h8_p2 - fs_147_1342198 * r_2 * h8_p4 - fs_20_2601 * r_4 * h6_p2 + fs_338_104907 * r_4 * h6_p4 + fs_270_20449 * r_6 * h4_p4 + fs_1400_184041 * r_8 * h2_p2);
    }

    // NOTE: the rows are formed in 53 loops, as the vectorizer runs out of
    // registers with all 66 of them in one.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, ph4_p3, ph6_p3, ph6_p5, ph8_p3, ph8_p5, ph10_p3, ph10_p5, ab_2, pc_54 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h4_p3 = ph4_p3[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h6_p5 = ph6_p5[k];
        const auto h8_p3 = ph8_p3[k];
        const auto h8_p5 = ph8_p5[k];
        const auto h10_p3 = ph10_p3[k];
        const auto h10_p5 = ph10_p5[k];

        pc_54[k] = -fs_30375_32 * e_1 * h4_p3 + e_2 * (fs_1125_968 * h6_p3 + fs_6075_88 * h6_p5 + fs_273375_968 * r_2 * h4_p3) + e_3 * (fs_3920_1859 * h8_p3 - fs_147_143 * h8_p5 - fs_10_121 * r_2 * h6_p3 - fs_54_11 * r_2 * h6_p5 - fs_273375_40898 * r_4 * h4_p3) + e_4 * (fs_138915_29838094 * h10_p3 - fs_694575_29838094 * h10_p5 - fs_15680_671099 * r_2 * h8_p3 + fs_588_51623 * r_2 * h8_p5 + fs_10_34969 * r_4 * h6_p3 + fs_54_3179 * r_4 * h6_p5 + fs_270_20449 * r_6 * h4_p3);
    }

    // NOTE: the rows are formed in 53 loops, as the vectorizer runs out of
    // registers with all 66 of them in one.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, ph4_p4, ph6_p4, ph6_p6, ph8_p4, ph8_p6, ph10_p4, ph10_p6, ab_2, pc_55 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h4_p4 = ph4_p4[k];
        const auto h6_p4 = ph6_p4[k];
        const auto h6_p6 = ph6_p6[k];
        const auto h8_p4 = ph8_p4[k];
        const auto h8_p6 = ph8_p6[k];
        const auto h10_p4 = ph10_p4[k];
        const auto h10_p6 = ph10_p6[k];

        pc_55[k] = fs_6075_32 * e_1 * h4_p4 + e_2 * (fs_84375_968 * h6_p4 - fs_1125_44 * h6_p6 - fs_54675_968 * r_2 * h4_p4) + e_3 * (fs_18375_14872 * h8_p4 - fs_1715_572 * h8_p6 - fs_750_121 * r_2 * h6_p4 + fs_20_11 * r_2 * h6_p6 + fs_54675_40898 * r_4 * h4_p4) + e_4 * (fs_27783_29838094 * h10_p4 - fs_111132_14919047 * h10_p6 - fs_18375_1342198 * r_2 * h8_p4 + fs_1715_51623 * r_2 * h8_p6 + fs_750_34969 * r_4 * h6_p4 - fs_20_3179 * r_4 * h6_p6 - fs_54_20449 * r_6 * h4_p4);
    }

    // NOTE: the rows are formed in 53 loops, as the vectorizer runs out of
    // registers with all 66 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph2_0, ph4_0, ph4_p4, ph6_0, ph6_p4, ph8_0, ph8_p4, ph10_0, ph10_p4, ab_2, pc_56 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;
        const auto r_10 = r_8 * r_2;

        const auto h2_0 = ph2_0[k];
        const auto h4_0 = ph4_0[k];
        const auto h4_p4 = ph4_p4[k];
        const auto h6_0 = ph6_0[k];
        const auto h6_p4 = ph6_p4[k];
        const auto h8_0 = ph8_0[k];
        const auto h8_p4 = ph8_p4[k];
        const auto h10_0 = ph10_0[k];
        const auto h10_p4 = ph10_p4[k];

        pc_56[k] = e_0 * (-f_315_8 * h2_0 - f_1575_16 * r_2) + e_1 * (-f_45_8 * h4_0 + fs_70875_64 * h4_p4 + f_45 * r_2 * h2_0 + f_315_4 * r_4) + e_2 * (f_90_11 * h6_0 - fs_6300_121 * h6_p4 + f_135_44 * r_2 * h4_0 - fs_637875_1936 * r_2 * h4_p4 - f_15 * r_4 * h2_0 - f_45_2 * r_6) + e_3 * (-f_238_143 * h8_0 + fs_3087_1859 * h8_p4 - f_24_11 * r_2 * h6_0 + fs_448_121 * r_2 * h6_p4 - f_135_286 * r_4 * h4_0 + fs_637875_81796 * r_4 * h4_p4 + f_20_11 * r_6 * h2_0 + f_5_2 * r_8) + e_4 * (-f_7560_46189 * h10_0 - fs_952560_14919047 * h10_p4 + f_476_2717 * r_2 * h8_0 - fs_12348_671099 * r_2 * h8_p4 + f_24_187 * r_4 * h6_0 - fs_448_34969 * r_4 * h6_p4 + f_3_143 * r_6 * h4_0 - fs_315_20449 * r_6 * h4_p4 - f_10_143 * r_8 * h2_0 - f_1_11 * r_10) + f_945_32 * e_5;
    }

    // NOTE: the rows are formed in 53 loops, as the vectorizer runs out of
    // registers with all 66 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_p1, ph4_p1, ph6_p1, ph6_p5, ph8_p1, ph8_p5, ph10_p1, ph10_p5, ab_2, pc_57 : simd::cache_line_size())
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

        const auto h2_p1 = ph2_p1[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h6_p5 = ph6_p5[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h8_p5 = ph8_p5[k];
        const auto h10_p1 = ph10_p1[k];
        const auto h10_p5 = ph10_p5[k];

        pc_57[k] = -fs_275625_128 * e_0 * h2_p1 + e_1 * (fs_30375_64 * h4_p1 + fs_5625_2 * r_2 * h2_p1) + e_2 * (fs_175_242 * h6_p1 - fs_525_44 * h6_p5 - fs_273375_1936 * r_2 * h4_p1 - fs_625_2 * r_4 * h2_p1) + e_3 * (-fs_324723_163592 * h8_p1 + fs_1029_1144 * h8_p5 - fs_56_1089 * r_2 * h6_p1 + fs_28_33 * r_2 * h6_p5 + fs_273375_81796 * r_4 * h4_p1 + fs_5000_1089 * r_6 * h2_p1) + e_4 * (-fs_1071630_193947611 * h10_p1 - fs_893025_14919047 * h10_p5 + fs_324723_14764178 * r_2 * h8_p1 - fs_1029_103246 * r_2 * h8_p5 + fs_56_314721 * r_4 * h6_p1 - fs_28_9537 * r_4 * h6_p5 - fs_135_20449 * r_6 * h4_p1 - fs_1250_184041 * r_8 * h2_p1);
    }

    // NOTE: the rows are formed in 53 loops, as the vectorizer runs out of
    // registers with all 66 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_p2, ph4_p2, ph6_p2, ph6_p6, ph8_p2, ph8_p6, ph10_p2, ph10_p6, ab_2, pc_58 : simd::cache_line_size())
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

        const auto h2_p2 = ph2_p2[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h6_p6 = ph6_p6[k];
        const auto h8_p2 = ph8_p2[k];
        const auto h8_p6 = ph8_p6[k];
        const auto h10_p2 = ph10_p2[k];
        const auto h10_p6 = ph10_p6[k];

        pc_58[k] = f_315_8 * e_0 * h2_p2 + e_1 * (fs_30375_64 * h4_p2 - f_45 * r_2 * h2_p2) + e_2 * (-fs_7875_242 * h6_p2 + fs_1575_22 * h6_p6 - fs_273375_1936 * r_2 * h4_p2 + f_15 * r_4 * h2_p2) + e_3 * (-fs_252105_163592 * h8_p2 - fs_49_1144 * h8_p6 + fs_280_121 * r_2 * h6_p2 - fs_56_11 * r_2 * h6_p6 + fs_273375_81796 * r_4 * h4_p2 - f_20_11 * r_6 * h2_p2) + e_4 * (-fs_317520_193947611 * h10_p2 - fs_635040_14919047 * h10_p6 + fs_252105_14764178 * r_2 * h8_p2 + fs_49_103246 * r_2 * h8_p6 - fs_280_34969 * r_4 * h6_p2 + fs_56_3179 * r_4 * h6_p6 - fs_135_20449 * r_6 * h4_p2 + f_10_143 * r_8 * h2_p2);
    }

    // NOTE: the rows are formed in 53 loops, as the vectorizer runs out of
    // registers with all 66 of them in one.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, ph4_p3, ph6_p3, ph8_p3, ph8_p7, ph10_p3, ph10_p7, ab_2, pc_59 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h4_p3 = ph4_p3[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h8_p3 = ph8_p3[k];
        const auto h8_p7 = ph8_p7[k];
        const auto h10_p3 = ph10_p3[k];
        const auto h10_p7 = ph10_p7[k];

        pc_59[k] = -fs_42525_64 * e_1 * h4_p3 + e_2 * (-fs_39375_484 * h6_p3 + fs_382725_1936 * r_2 * h4_p3) + e_3 * (-fs_8575_14872 * h8_p3 - fs_3675_1144 * h8_p7 + fs_700_121 * r_2 * h6_p3 - fs_382725_81796 * r_4 * h4_p3) + e_4 * (-fs_3969_14919047 * h10_p3 - fs_15876_877591 * h10_p7 + fs_8575_1342198 * r_2 * h8_p3 + fs_3675_103246 * r_2 * h8_p7 - fs_700_34969 * r_4 * h6_p3 + fs_189_20449 * r_6 * h4_p3);
    }

    // NOTE: the rows are formed in 53 loops, as the vectorizer runs out of
    // registers with all 66 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph2_0, ph4_0, ph6_0, ph6_p6, ph8_0, ph8_p6, ph10_0, ph10_p6, ab_2, pc_60 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;
        const auto r_10 = r_8 * r_2;

        const auto h2_0 = ph2_0[k];
        const auto h4_0 = ph4_0[k];
        const auto h6_0 = ph6_0[k];
        const auto h6_p6 = ph6_p6[k];
        const auto h8_0 = ph8_0[k];
        const auto h8_p6 = ph8_p6[k];
        const auto h10_0 = ph10_0[k];
        const auto h10_p6 = ph10_p6[k];

        pc_60[k] = e_0 * (-f_105_16 * h2_0 - f_1575_16 * r_2) + e_1 * (-f_135_4 * h4_0 + f_15_2 * r_2 * h2_0 + f_315_4 * r_4) + e_2 * (f_145_22 * h6_0 - fs_1050_11 * h6_p6 + f_405_22 * r_2 * h4_0 - f_5_2 * r_4 * h2_0 - f_45_2 * r_6) + e_3 * (f_511_286 * h8_0 + fs_294_143 * h8_p6 - f_58_33 * r_2 * h6_0 + fs_224_33 * r_2 * h6_p6 - f_405_143 * r_4 * h4_0 + f_10_33 * r_6 * h2_0 + f_5_2 * r_8) + e_4 * (f_2835_46189 * h10_0 - fs_1071630_14919047 * h10_p6 - f_511_2717 * r_2 * h8_0 - fs_1176_51623 * r_2 * h8_p6 + f_58_561 * r_4 * h6_0 - fs_224_9537 * r_4 * h6_p6 + f_18_143 * r_6 * h4_0 - f_5_429 * r_8 * h2_0 - f_1_11 * r_10) + f_945_32 * e_5;
    }

    // NOTE: the rows are formed in 53 loops, as the vectorizer runs out of
    // registers with all 66 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_p1, ph6_p1, ph8_p1, ph8_p7, ph10_p1, ph10_p7, ab_2, pc_61 : simd::cache_line_size())
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

        const auto h2_p1 = ph2_p1[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h8_p7 = ph8_p7[k];
        const auto h10_p1 = ph10_p1[k];
        const auto h10_p7 = ph10_p7[k];

        pc_61[k] = -fs_1620675_512 * e_0 * h2_p1 + fs_33075_8 * e_1 * r_2 * h2_p1 + e_2 * (fs_525_8 * h6_p1 - fs_3675_8 * r_4 * h2_p1) + e_3 * (fs_98_121 * h8_p1 + fs_245_286 * h8_p7 - fs_14_3 * r_2 * h6_p1 + fs_2450_363 * r_6 * h2_p1) + e_4 * (fs_178605_387895222 * h10_p1 - fs_59535_877591 * h10_p7 - fs_392_43681 * r_2 * h8_p1 - fs_490_51623 * r_2 * h8_p7 + fs_14_867 * r_4 * h6_p1 - fs_1225_122694 * r_8 * h2_p1);
    }

    // NOTE: the rows are formed in 53 loops, as the vectorizer runs out of
    // registers with all 66 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_p2, ph4_p2, ph6_p2, ph8_p2, ph8_p8, ph10_p2, ph10_p8, ab_2, pc_62 : simd::cache_line_size())
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

        const auto h2_p2 = ph2_p2[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h8_p2 = ph8_p2[k];
        const auto h8_p8 = ph8_p8[k];
        const auto h10_p2 = ph10_p2[k];
        const auto h10_p8 = ph10_p8[k];

        pc_62[k] = fs_165375_256 * e_0 * h2_p2 + e_1 * (f_135_4 * h4_p2 - fs_3375_4 * r_2 * h2_p2) + e_2 * (fs_13125_242 * h6_p2 - f_405_22 * r_2 * h4_p2 + fs_375_4 * r_4 * h2_p2) + e_3 * (fs_8575_40898 * h8_p2 - fs_1225_572 * h8_p8 - fs_1400_363 * r_2 * h6_p2 + f_405_143 * r_4 * h4_p2 - fs_500_363 * r_6 * h2_p2) + e_4 * (fs_11907_193947611 * h10_p2 - fs_35721_877591 * h10_p8 - fs_17150_7382089 * r_2 * h8_p2 + fs_1225_51623 * r_2 * h8_p8 + fs_1400_104907 * r_4 * h6_p2 - f_18_143 * r_6 * h4_p2 + fs_125_61347 * r_8 * h2_p2);
    }

    // NOTE: the rows are formed in 53 loops, as the vectorizer runs out of
    // registers with all 66 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph2_0, ph4_0, ph6_0, ph8_0, ph8_p8, ph10_0, ph10_p8, ab_2, pc_63 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;
        const auto r_10 = r_8 * r_2;

        const auto h2_0 = ph2_0[k];
        const auto h4_0 = ph4_0[k];
        const auto h6_0 = ph6_0[k];
        const auto h8_0 = ph8_0[k];
        const auto h8_p8 = ph8_p8[k];
        const auto h10_0 = ph10_0[k];
        const auto h10_p8 = ph10_p8[k];

        pc_63[k] = e_0 * (f_315_8 * h2_0 - f_1575_16 * r_2) + e_1 * (-f_135_4 * h4_0 - f_45 * r_2 * h2_0 + f_315_4 * r_4) + e_2 * (-f_120_11 * h6_0 + f_405_22 * r_2 * h4_0 + f_15 * r_4 * h2_0 - f_45_2 * r_6) + e_3 * (-f_217_286 * h8_0 + fs_2205_572 * h8_p8 + f_32_11 * r_2 * h6_0 - f_405_143 * r_4 * h4_0 - f_20_11 * r_6 * h2_0 + f_5_2 * r_8) + e_4 * (-f_630_46189 * h10_0 - fs_79380_877591 * h10_p8 + f_217_2717 * r_2 * h8_0 - fs_2205_51623 * r_2 * h8_p8 - f_32_187 * r_4 * h6_0 + f_18_143 * r_6 * h4_0 + f_10_143 * r_8 * h2_0 - f_1_11 * r_10) + f_945_32 * e_5;
    }

    // NOTE: the rows are formed in 53 loops, as the vectorizer runs out of
    // registers with all 66 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_p1, ph4_p1, ph6_p1, ph8_p1, ph10_p1, ph10_p9, ab_2, pc_64 : simd::cache_line_size())
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

        const auto h2_p1 = ph2_p1[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h10_p1 = ph10_p1[k];
        const auto h10_p9 = ph10_p9[k];

        pc_64[k] = -fs_1488375_512 * e_0 * h2_p1 + e_1 * (-f_135_4 * h4_p1 + fs_30375_8 * r_2 * h2_p1) + e_2 * (-fs_23625_968 * h6_p1 + f_405_22 * r_2 * h4_p1 - fs_3375_8 * r_4 * h2_p1) + e_3 * (-fs_2205_40898 * h8_p1 + fs_210_121 * r_2 * h6_p1 - f_405_143 * r_4 * h4_p1 + fs_750_121 * r_6 * h2_p1) + e_4 * (-fs_3969_387895222 * h10_p1 - fs_3969_46189 * h10_p9 + fs_4410_7382089 * r_2 * h8_p1 - fs_210_34969 * r_4 * h6_p1 + f_18_143 * r_6 * h4_p1 - fs_375_40898 * r_8 * h2_p1);
    }

    // NOTE: the rows are formed in 53 loops, as the vectorizer runs out of
    // registers with all 66 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph2_0, ph4_0, ph6_0, ph8_0, ph10_0, ph10_p10, ab_2, pc_65 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;
        const auto r_10 = r_8 * r_2;

        const auto h2_0 = ph2_0[k];
        const auto h4_0 = ph4_0[k];
        const auto h6_0 = ph6_0[k];
        const auto h8_0 = ph8_0[k];
        const auto h10_0 = ph10_0[k];
        const auto h10_p10 = ph10_p10[k];

        pc_65[k] = e_0 * (f_1575_16 * h2_0 - f_1575_16 * r_2) + e_1 * (f_135_4 * h4_0 - f_225_2 * r_2 * h2_0 + f_315_4 * r_4) + e_2 * (f_75_22 * h6_0 - f_405_22 * r_2 * h4_0 + f_75_2 * r_4 * h2_0 - f_45_2 * r_6) + e_3 * (f_35_286 * h8_0 - f_10_11 * r_2 * h6_0 + f_405_143 * r_4 * h4_0 - f_50_11 * r_6 * h2_0 + f_5_2 * r_8) + e_4 * (f_63_46189 * h10_0 - fs_7938_46189 * h10_p10 - f_35_2717 * r_2 * h8_0 + f_10_187 * r_4 * h6_0 - f_18_143 * r_6 * h4_0 + f_25_143 * r_8 * h2_0 - f_1_11 * r_10) + f_945_32 * e_5;
    }

    // NOTE: the values of a combination of angular components are stored as one
    // row of nvalues columns, with the component on bra side running slowest, and
    // the atom pairs beyond the reach of every pair of primitives are set to zero.

    const size_t sources[121] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 1, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 2, 12, 21, 22, 23, 24, 25, 26, 27, 28, 29, 3, 13, 22, 30, 31, 32, 33, 34, 35, 36, 37, 4, 14, 23, 31, 38, 39, 40, 41, 42, 43, 44, 5, 15, 24, 32, 39, 45, 46, 47, 48, 49, 50, 6, 16, 25, 33, 40, 46, 51, 52, 53, 54, 55, 7, 17, 26, 34, 41, 47, 52, 56, 57, 58, 59, 8, 18, 27, 35, 42, 48, 53, 57, 60, 61, 62, 9, 19, 28, 36, 43, 49, 54, 58, 61, 63, 64, 10, 20, 29, 37, 44, 50, 55, 59, 62, 64, 65};

    for (size_t m = 0; m < 121; m++)
    {
        const auto *pc = buffer.data(6 + sources[m]);

        std::copy(pc, pc + nmax, values + m * nvalues);

        std::fill(values + m * nvalues + nmax, values + (m + 1) * nvalues, 0.0);
    }
}

}  // namespace simdovl
