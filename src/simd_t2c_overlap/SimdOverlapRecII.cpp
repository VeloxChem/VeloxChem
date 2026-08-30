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
#include <string>

#include "ErrorHandler.hpp"
#include "MathConst.hpp"
#include "ScreeningFunc.hpp"
#include "SimdAlign.hpp"
#include "SimdDimensions.hpp"

namespace simdovl {  // simdovl namespace

auto
compute_ii_overlap(double                         *values,
                   const size_t                    nvalues,
                   const CBasisFunction           &bra,
                   const CBasisFunction           &ket,
                   const std::vector<CSimdMatrix> &harmonics,
                   const CSimdMatrix              &coordinates,
                   const double                    threshold) -> void
{
    if ((bra.get_angular_momentum() != 6) || (ket.get_angular_momentum() != 6))
    {
        errors::assertMsgCritical(
            false, std::string("SimdOverlapRecII.compute_ii_overlap: Basis functions must be of angular momenta six and six"));
    }

    if (harmonics.size() < 12)
    {
        errors::assertMsgCritical(
            false, std::string("SimdOverlapRecII.compute_ii_overlap: Harmonics must reach angular momentum twelve"));
    }

    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("SimdOverlapRecII.compute_ii_overlap: Number of values exceeds number of atom pairs"));
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
        std::fill(values, values + 169 * nvalues, 0.0);

        return;
    }

    // NOTE: the first seven rows accumulate the contracted prefactors of the terms,
    // and the remaining 91 rows hold the integrals of the combinations of angular
    // components which are not related by symmetry.

    auto buffer = CSimdMatrix(98, nmax);

    auto *pe_0 = buffer.data(0);
    auto *pe_1 = buffer.data(1);
    auto *pe_2 = buffer.data(2);
    auto *pe_3 = buffer.data(3);
    auto *pe_4 = buffer.data(4);
    auto *pe_5 = buffer.data(5);
    auto *pe_6 = buffer.data(6);

    std::fill(pe_0, pe_0 + nmax, 0.0);
    std::fill(pe_1, pe_1 + nmax, 0.0);
    std::fill(pe_2, pe_2 + nmax, 0.0);
    std::fill(pe_3, pe_3 + nmax, 0.0);
    std::fill(pe_4, pe_4 + nmax, 0.0);
    std::fill(pe_5, pe_5 + nmax, 0.0);
    std::fill(pe_6, pe_6 + nmax, 0.0);

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

            const auto f_0 = fbase * fmu / fexp / fexp / fexp / fexp / fexp / fexp;

            const auto f_1 = fbase * fmu * fmu / fexp / fexp / fexp / fexp / fexp / fexp;

            const auto f_2 = fbase * fmu * fmu * fmu / fexp / fexp / fexp / fexp / fexp / fexp;

            const auto f_3 = fbase * fmu * fmu * fmu * fmu / fexp / fexp / fexp / fexp / fexp / fexp;

            const auto f_4 = fbase * fmu * fmu * fmu * fmu * fmu / fexp / fexp / fexp / fexp / fexp / fexp;

            const auto f_5 = fbase * fmu * fmu * fmu * fmu * fmu * fmu / fexp / fexp / fexp / fexp / fexp / fexp;

            const auto f_6 = fbase / fexp / fexp / fexp / fexp / fexp / fexp;

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
    const auto *ph12_m12 = harmonics[11].data(0);
    const auto *ph12_m11 = harmonics[11].data(1);
    const auto *ph12_m10 = harmonics[11].data(2);
    const auto *ph12_m9 = harmonics[11].data(3);
    const auto *ph12_m8 = harmonics[11].data(4);
    const auto *ph12_m7 = harmonics[11].data(5);
    const auto *ph12_m6 = harmonics[11].data(6);
    const auto *ph12_m5 = harmonics[11].data(7);
    const auto *ph12_m4 = harmonics[11].data(8);
    const auto *ph12_m3 = harmonics[11].data(9);
    const auto *ph12_m2 = harmonics[11].data(10);
    const auto *ph12_m1 = harmonics[11].data(11);
    const auto *ph12_0 = harmonics[11].data(12);
    const auto *ph12_p1 = harmonics[11].data(13);
    const auto *ph12_p2 = harmonics[11].data(14);
    const auto *ph12_p3 = harmonics[11].data(15);
    const auto *ph12_p4 = harmonics[11].data(16);
    const auto *ph12_p5 = harmonics[11].data(17);
    const auto *ph12_p6 = harmonics[11].data(18);
    const auto *ph12_p7 = harmonics[11].data(19);
    const auto *ph12_p8 = harmonics[11].data(20);
    const auto *ph12_p9 = harmonics[11].data(21);
    const auto *ph12_p10 = harmonics[11].data(22);
    const auto *ph12_p11 = harmonics[11].data(23);
    const auto *ph12_p12 = harmonics[11].data(24);

    auto *pc_0 = buffer.data(7);
    auto *pc_1 = buffer.data(8);
    auto *pc_2 = buffer.data(9);
    auto *pc_3 = buffer.data(10);
    auto *pc_4 = buffer.data(11);
    auto *pc_5 = buffer.data(12);
    auto *pc_6 = buffer.data(13);
    auto *pc_7 = buffer.data(14);
    auto *pc_8 = buffer.data(15);
    auto *pc_9 = buffer.data(16);
    auto *pc_10 = buffer.data(17);
    auto *pc_11 = buffer.data(18);
    auto *pc_12 = buffer.data(19);
    auto *pc_13 = buffer.data(20);
    auto *pc_14 = buffer.data(21);
    auto *pc_15 = buffer.data(22);
    auto *pc_16 = buffer.data(23);
    auto *pc_17 = buffer.data(24);
    auto *pc_18 = buffer.data(25);
    auto *pc_19 = buffer.data(26);
    auto *pc_20 = buffer.data(27);
    auto *pc_21 = buffer.data(28);
    auto *pc_22 = buffer.data(29);
    auto *pc_23 = buffer.data(30);
    auto *pc_24 = buffer.data(31);
    auto *pc_25 = buffer.data(32);
    auto *pc_26 = buffer.data(33);
    auto *pc_27 = buffer.data(34);
    auto *pc_28 = buffer.data(35);
    auto *pc_29 = buffer.data(36);
    auto *pc_30 = buffer.data(37);
    auto *pc_31 = buffer.data(38);
    auto *pc_32 = buffer.data(39);
    auto *pc_33 = buffer.data(40);
    auto *pc_34 = buffer.data(41);
    auto *pc_35 = buffer.data(42);
    auto *pc_36 = buffer.data(43);
    auto *pc_37 = buffer.data(44);
    auto *pc_38 = buffer.data(45);
    auto *pc_39 = buffer.data(46);
    auto *pc_40 = buffer.data(47);
    auto *pc_41 = buffer.data(48);
    auto *pc_42 = buffer.data(49);
    auto *pc_43 = buffer.data(50);
    auto *pc_44 = buffer.data(51);
    auto *pc_45 = buffer.data(52);
    auto *pc_46 = buffer.data(53);
    auto *pc_47 = buffer.data(54);
    auto *pc_48 = buffer.data(55);
    auto *pc_49 = buffer.data(56);
    auto *pc_50 = buffer.data(57);
    auto *pc_51 = buffer.data(58);
    auto *pc_52 = buffer.data(59);
    auto *pc_53 = buffer.data(60);
    auto *pc_54 = buffer.data(61);
    auto *pc_55 = buffer.data(62);
    auto *pc_56 = buffer.data(63);
    auto *pc_57 = buffer.data(64);
    auto *pc_58 = buffer.data(65);
    auto *pc_59 = buffer.data(66);
    auto *pc_60 = buffer.data(67);
    auto *pc_61 = buffer.data(68);
    auto *pc_62 = buffer.data(69);
    auto *pc_63 = buffer.data(70);
    auto *pc_64 = buffer.data(71);
    auto *pc_65 = buffer.data(72);
    auto *pc_66 = buffer.data(73);
    auto *pc_67 = buffer.data(74);
    auto *pc_68 = buffer.data(75);
    auto *pc_69 = buffer.data(76);
    auto *pc_70 = buffer.data(77);
    auto *pc_71 = buffer.data(78);
    auto *pc_72 = buffer.data(79);
    auto *pc_73 = buffer.data(80);
    auto *pc_74 = buffer.data(81);
    auto *pc_75 = buffer.data(82);
    auto *pc_76 = buffer.data(83);
    auto *pc_77 = buffer.data(84);
    auto *pc_78 = buffer.data(85);
    auto *pc_79 = buffer.data(86);
    auto *pc_80 = buffer.data(87);
    auto *pc_81 = buffer.data(88);
    auto *pc_82 = buffer.data(89);
    auto *pc_83 = buffer.data(90);
    auto *pc_84 = buffer.data(91);
    auto *pc_85 = buffer.data(92);
    auto *pc_86 = buffer.data(93);
    auto *pc_87 = buffer.data(94);
    auto *pc_88 = buffer.data(95);
    auto *pc_89 = buffer.data(96);
    auto *pc_90 = buffer.data(97);

    // NOTE: the factors of the terms depend on the angular momenta alone, so they
    // are formed once for the whole matrix instead of once for every atom pair.

    const auto f_10395_16 = 649.6875;
    const auto f_4455_16 = 278.4375;
    const auto f_7425_8 = 928.125;
    const auto f_75_2 = 37.5;
    const auto f_405_2 = 202.5;
    const auto f_825_2 = 412.5;
    const auto f_495_2 = 247.5;
    const auto f_105_52 = 105.0 / 52.0;
    const auto f_15 = 15.0;
    const auto f_1215_26 = 1215.0 / 26.0;
    const auto f_75 = 75.0;
    const auto f_165_4 = 41.25;
    const auto f_189_4199 = 189.0 / 4199.0;
    const auto f_105_247 = 105.0 / 247.0;
    const auto f_30_17 = 30.0 / 17.0;
    const auto f_54_13 = 54.0 / 13.0;
    const auto f_75_13 = 75.0 / 13.0;
    const auto f_3 = 3.0;
    const auto f_33_96577 = 33.0 / 96577.0;
    const auto fs_15246_96577 = std::sqrt(15246.0 / 96577.0);
    const auto f_378_96577 = 378.0 / 96577.0;
    const auto f_5_247 = 5.0 / 247.0;
    const auto f_20_323 = 20.0 / 323.0;
    const auto f_27_221 = 27.0 / 221.0;
    const auto f_2_13 = 2.0 / 13.0;
    const auto f_1_13 = 1.0 / 13.0;
    const auto f_10395_64 = 162.421875;
    const auto f_10395_32 = 324.84375;
    const auto fs_33078375_512 = std::sqrt(64606.201171875);
    const auto f_7425_16 = 464.0625;
    const auto fs_39375_16 = std::sqrt(2460.9375);
    const auto fs_273375_8 = std::sqrt(34171.875);
    const auto f_825_4 = 206.25;
    const auto fs_33075_2704 = std::sqrt(33075.0 / 2704.0);
    const auto fs_1575_4 = std::sqrt(393.75);
    const auto fs_2460375_1352 = std::sqrt(2460375.0 / 1352.0);
    const auto fs_654885_70526404 = std::sqrt(654885.0 / 70526404.0);
    const auto fs_33075_61009 = std::sqrt(33075.0 / 61009.0);
    const auto fs_1575_289 = std::sqrt(1575.0 / 289.0);
    const auto fs_2430_169 = std::sqrt(2430.0 / 169.0);
    const auto f_75_26 = 75.0 / 26.0;
    const auto fs_1089_1434941066 = std::sqrt(1089.0 / 1434941066.0);
    const auto fs_7623_96577 = std::sqrt(7623.0 / 96577.0);
    const auto fs_654885_9327116929 = std::sqrt(654885.0 / 9327116929.0);
    const auto fs_75_61009 = std::sqrt(75.0 / 61009.0);
    const auto fs_700_104329 = std::sqrt(700.0 / 104329.0);
    const auto fs_1215_97682 = std::sqrt(1215.0 / 97682.0);
    const auto fs_9823275_512 = std::sqrt(19186.083984375);
    const auto fs_27064125_512 = std::sqrt(52859.619140625);
    const auto fs_5011875_128 = std::sqrt(39155.2734375);
    const auto fs_196875_44 = std::sqrt(196875.0 / 44.0);
    const auto fs_2460375_88 = std::sqrt(2460375.0 / 88.0);
    const auto fs_61875_8 = std::sqrt(7734.375);
    const auto fs_1157625_29744 = std::sqrt(1157625.0 / 29744.0);
    const auto fs_7875_11 = std::sqrt(7875.0 / 11.0);
    const auto fs_22143375_14872 = std::sqrt(22143375.0 / 14872.0);
    const auto fs_5625_22 = std::sqrt(5625.0 / 22.0);
    const auto fs_1607445_35263202 = std::sqrt(1607445.0 / 35263202.0);
    const auto fs_11907_4199 = std::sqrt(11907.0 / 4199.0);
    const auto fs_1157625_671099 = std::sqrt(1157625.0 / 671099.0);
    const auto fs_31500_3179 = std::sqrt(31500.0 / 3179.0);
    const auto fs_21870_1859 = std::sqrt(21870.0 / 1859.0);
    const auto fs_5625_3718 = std::sqrt(5625.0 / 3718.0);
    const auto fs_7623_1434941066 = std::sqrt(7623.0 / 1434941066.0);
    const auto fs_83853_2221271 = std::sqrt(83853.0 / 2221271.0);
    const auto fs_3214890_9327116929 = std::sqrt(3214890.0 / 9327116929.0);
    const auto fs_47628_2221271 = std::sqrt(47628.0 / 2221271.0);
    const auto fs_2625_671099 = std::sqrt(2625.0 / 671099.0);
    const auto fs_14000_1147619 = std::sqrt(14000.0 / 1147619.0);
    const auto fs_10935_1074502 = std::sqrt(10935.0 / 1074502.0);
    const auto fs_2_1859 = std::sqrt(2.0 / 1859.0);
    const auto fs_12629925_512 = std::sqrt(24667.822265625);
    const auto fs_118125_22 = std::sqrt(118125.0 / 22.0);
    const auto fs_1148175_88 = std::sqrt(1148175.0 / 88.0);
    const auto fs_231525_2704 = std::sqrt(231525.0 / 2704.0);
    const auto fs_9450_11 = std::sqrt(9450.0 / 11.0);
    const auto fs_10333575_14872 = std::sqrt(10333575.0 / 14872.0);
    const auto fs_214326_1356277 = std::sqrt(214326.0 / 1356277.0);
    const auto fs_35721_8398 = std::sqrt(35721.0 / 8398.0);
    const auto fs_231525_61009 = std::sqrt(231525.0 / 61009.0);
    const auto fs_37800_3179 = std::sqrt(37800.0 / 3179.0);
    const auto fs_10206_1859 = std::sqrt(10206.0 / 1859.0);
    const auto fs_38115_1434941066 = std::sqrt(38115.0 / 1434941066.0);
    const auto fs_38115_2221271 = std::sqrt(38115.0 / 2221271.0);
    const auto fs_857304_717470533 = std::sqrt(857304.0 / 717470533.0);
    const auto fs_71442_2221271 = std::sqrt(71442.0 / 2221271.0);
    const auto fs_525_61009 = std::sqrt(525.0 / 61009.0);
    const auto fs_16800_1147619 = std::sqrt(16800.0 / 1147619.0);
    const auto fs_5103_1074502 = std::sqrt(5103.0 / 1074502.0);
    const auto fs_1403325_256 = std::sqrt(5481.73828125);
    const auto fs_127575_44 = std::sqrt(127575.0 / 44.0);
    const auto fs_385875_2704 = std::sqrt(385875.0 / 2704.0);
    const auto fs_11025_208 = std::sqrt(11025.0 / 208.0);
    const auto fs_1148175_7436 = std::sqrt(1148175.0 / 7436.0);
    const auto fs_583443_1356277 = std::sqrt(583443.0 / 1356277.0);
    const auto fs_321489_79781 = std::sqrt(321489.0 / 79781.0);
    const auto fs_385875_61009 = std::sqrt(385875.0 / 61009.0);
    const auto fs_11025_4693 = std::sqrt(11025.0 / 4693.0);
    const auto fs_2268_1859 = std::sqrt(2268.0 / 1859.0);
    const auto fs_76230_717470533 = std::sqrt(76230.0 / 717470533.0);
    const auto fs_16335_2221271 = std::sqrt(16335.0 / 2221271.0);
    const auto fs_2333772_717470533 = std::sqrt(2333772.0 / 717470533.0);
    const auto fs_1285956_42204149 = std::sqrt(1285956.0 / 42204149.0);
    const auto fs_875_61009 = std::sqrt(875.0 / 61009.0);
    const auto fs_25_4693 = std::sqrt(25.0 / 4693.0);
    const auto fs_567_537251 = std::sqrt(567.0 / 537251.0);
    const auto fs_77175_416 = std::sqrt(77175.0 / 416.0);
    const auto fs_55125_416 = std::sqrt(55125.0 / 416.0);
    const auto fs_5250987_5425108 = std::sqrt(5250987.0 / 5425108.0);
    const auto fs_238140_79781 = std::sqrt(238140.0 / 79781.0);
    const auto fs_77175_9386 = std::sqrt(77175.0 / 9386.0);
    const auto fs_55125_9386 = std::sqrt(55125.0 / 9386.0);
    const auto fs_15246_42204149 = std::sqrt(15246.0 / 42204149.0);
    const auto fs_6534_2221271 = std::sqrt(6534.0 / 2221271.0);
    const auto fs_5250987_717470533 = std::sqrt(5250987.0 / 717470533.0);
    const auto fs_952560_42204149 = std::sqrt(952560.0 / 42204149.0);
    const auto fs_175_9386 = std::sqrt(175.0 / 9386.0);
    const auto fs_125_9386 = std::sqrt(125.0 / 9386.0);
    const auto fs_77175_208 = std::sqrt(77175.0 / 208.0);
    const auto fs_5000940_1356277 = std::sqrt(5000940.0 / 1356277.0);
    const auto fs_77175_4693 = std::sqrt(77175.0 / 4693.0);
    const auto fs_91476_42204149 = std::sqrt(91476.0 / 42204149.0);
    const auto fs_20003760_717470533 = std::sqrt(20003760.0 / 717470533.0);
    const auto fs_175_4693 = std::sqrt(175.0 / 4693.0);
    const auto f_1485_8 = 185.625;
    const auto f_375_4 = 93.75;
    const auto f_135 = 135.0;
    const auto f_525_52 = 525.0 / 52.0;
    const auto f_405_13 = 405.0 / 13.0;
    const auto f_3087_8398 = 3087.0 / 8398.0;
    const auto fs_43659_8398 = std::sqrt(43659.0 / 8398.0);
    const auto f_525_247 = 525.0 / 247.0;
    const auto f_75_17 = 75.0 / 17.0;
    const auto f_36_13 = 36.0 / 13.0;
    const auto f_396_96577 = 396.0 / 96577.0;
    const auto fs_182952_2221271 = std::sqrt(182952.0 / 2221271.0);
    const auto f_3087_96577 = 3087.0 / 96577.0;
    const auto fs_87318_2221271 = std::sqrt(87318.0 / 2221271.0);
    const auto f_25_247 = 25.0 / 247.0;
    const auto f_50_323 = 50.0 / 323.0;
    const auto f_18_221 = 18.0 / 221.0;
    const auto fs_265228425_2048 = std::sqrt(129506.06689453125);
    const auto fs_1002375_256 = std::sqrt(3915.52734375);
    const auto fs_135320625_512 = std::sqrt(264298.095703125);
    const auto fs_1063125_352 = std::sqrt(1063125.0 / 352.0);
    const auto fs_91125_44 = std::sqrt(91125.0 / 44.0);
    const auto fs_1670625_32 = std::sqrt(52207.03125);
    const auto fs_198450_1859 = std::sqrt(198450.0 / 1859.0);
    const auto fs_42525_88 = std::sqrt(42525.0 / 88.0);
    const auto fs_820125_7436 = std::sqrt(820125.0 / 7436.0);
    const auto fs_151875_88 = std::sqrt(151875.0 / 88.0);
    const auto fs_36693405_141052808 = std::sqrt(36693405.0 / 141052808.0);
    const auto fs_19845_16796 = std::sqrt(19845.0 / 16796.0);
    const auto fs_3175200_671099 = std::sqrt(3175200.0 / 671099.0);
    const auto fs_42525_6358 = std::sqrt(42525.0 / 6358.0);
    const auto fs_1620_1859 = std::sqrt(1620.0 / 1859.0);
    const auto fs_151875_14872 = std::sqrt(151875.0 / 14872.0);
    const auto fs_35937_717470533 = std::sqrt(35937.0 / 717470533.0);
    const auto fs_137214_2221271 = std::sqrt(137214.0 / 2221271.0);
    const auto fs_36693405_18654233858 = std::sqrt(36693405.0 / 18654233858.0);
    const auto fs_19845_2221271 = std::sqrt(19845.0 / 2221271.0);
    const auto fs_7200_671099 = std::sqrt(7200.0 / 671099.0);
    const auto fs_9450_1147619 = std::sqrt(9450.0 / 1147619.0);
    const auto fs_405_537251 = std::sqrt(405.0 / 537251.0);
    const auto fs_27_3718 = std::sqrt(27.0 / 3718.0);
    const auto fs_49116375_1024 = std::sqrt(47965.2099609375);
    const auto fs_601425_16 = std::sqrt(37589.0625);
    const auto fs_25059375_256 = std::sqrt(97888.18359375);
    const auto fs_39375_88 = std::sqrt(39375.0 / 88.0);
    const auto fs_218700_11 = std::sqrt(218700.0 / 11.0);
    const auto fs_309375_16 = std::sqrt(19335.9375);
    const auto fs_2083725_14872 = std::sqrt(2083725.0 / 14872.0);
    const auto fs_33075_208 = std::sqrt(33075.0 / 208.0);
    const auto fs_1575_22 = std::sqrt(1575.0 / 22.0);
    const auto fs_1968300_1859 = std::sqrt(1968300.0 / 1859.0);
    const auto fs_28125_44 = std::sqrt(28125.0 / 44.0);
    const auto f_6993_8398 = 6993.0 / 8398.0;
    const auto fs_11907_319124 = std::sqrt(11907.0 / 319124.0);
    const auto fs_4167450_671099 = std::sqrt(4167450.0 / 671099.0);
    const auto fs_33075_4693 = std::sqrt(33075.0 / 4693.0);
    const auto fs_3150_3179 = std::sqrt(3150.0 / 3179.0);
    const auto fs_15552_1859 = std::sqrt(15552.0 / 1859.0);
    const auto fs_28125_7436 = std::sqrt(28125.0 / 7436.0);
    const auto fs_152460_717470533 = std::sqrt(152460.0 / 717470533.0);
    const auto fs_87120_2221271 = std::sqrt(87120.0 / 2221271.0);
    const auto f_6993_96577 = 6993.0 / 96577.0;
    const auto fs_11907_42204149 = std::sqrt(11907.0 / 42204149.0);
    const auto fs_9450_671099 = std::sqrt(9450.0 / 671099.0);
    const auto fs_75_4693 = std::sqrt(75.0 / 4693.0);
    const auto fs_1400_1147619 = std::sqrt(1400.0 / 1147619.0);
    const auto fs_3888_537251 = std::sqrt(3888.0 / 537251.0);
    const auto fs_5_1859 = std::sqrt(5.0 / 1859.0);
    const auto fs_22920975_512 = std::sqrt(44767.529296875);
    const auto fs_2083725_88 = std::sqrt(2083725.0 / 88.0);
    const auto fs_77175_676 = std::sqrt(77175.0 / 676.0);
    const auto fs_18753525_14872 = std::sqrt(18753525.0 / 14872.0);
    const auto fs_3814209_2712554 = std::sqrt(3814209.0 / 2712554.0);
    const auto fs_194481_159562 = std::sqrt(194481.0 / 159562.0);
    const auto fs_308700_61009 = std::sqrt(308700.0 / 61009.0);
    const auto fs_18522_1859 = std::sqrt(18522.0 / 1859.0);
    const auto fs_1029105_1434941066 = std::sqrt(1029105.0 / 1434941066.0);
    const auto fs_49005_2221271 = std::sqrt(49005.0 / 2221271.0);
    const auto fs_7628418_717470533 = std::sqrt(7628418.0 / 717470533.0);
    const auto fs_388962_42204149 = std::sqrt(388962.0 / 42204149.0);
    const auto fs_700_61009 = std::sqrt(700.0 / 61009.0);
    const auto fs_9261_1074502 = std::sqrt(9261.0 / 1074502.0);
    const auto fs_2338875_128 = std::sqrt(18272.4609375);
    const auto fs_212625_22 = std::sqrt(212625.0 / 22.0);
    const auto fs_231525_5408 = std::sqrt(231525.0 / 5408.0);
    const auto fs_1913625_3718 = std::sqrt(1913625.0 / 3718.0);
    const auto fs_24310125_10850216 = std::sqrt(24310125.0 / 10850216.0);
    const auto fs_257985_104329 = std::sqrt(257985.0 / 104329.0);
    const auto fs_231525_122018 = std::sqrt(231525.0 / 122018.0);
    const auto fs_7560_1859 = std::sqrt(7560.0 / 1859.0);
    const auto fs_1463616_717470533 = std::sqrt(1463616.0 / 717470533.0);
    const auto fs_470448_42204149 = std::sqrt(470448.0 / 42204149.0);
    const auto fs_24310125_1434941066 = std::sqrt(24310125.0 / 1434941066.0);
    const auto fs_1031940_55190041 = std::sqrt(1031940.0 / 55190041.0);
    const auto fs_525_122018 = std::sqrt(525.0 / 122018.0);
    const auto fs_1890_537251 = std::sqrt(1890.0 / 537251.0);
    const auto fs_83349_15028 = std::sqrt(83349.0 / 15028.0);
    const auto fs_426888_42204149 = std::sqrt(426888.0 / 42204149.0);
    const auto fs_83349_1987453 = std::sqrt(83349.0 / 1987453.0);
    const auto f_945_16 = 59.0625;
    const auto f_270 = 270.0;
    const auto f_675_8 = 84.375;
    const auto f_150_11 = 150.0 / 11.0;
    const auto f_2160_11 = 2160.0 / 11.0;
    const auto f_9345_572 = 9345.0 / 572.0;
    const auto fs_496125_2288 = std::sqrt(496125.0 / 2288.0);
    const auto f_60_11 = 60.0 / 11.0;
    const auto f_6480_143 = 6480.0 / 143.0;
    const auto f_75_11 = 75.0 / 11.0;
    const auto f_5229_4199 = 5229.0 / 4199.0;
    const auto fs_218295_79781 = std::sqrt(218295.0 / 79781.0);
    const auto f_9345_2717 = 9345.0 / 2717.0;
    const auto fs_496125_51623 = std::sqrt(496125.0 / 51623.0);
    const auto f_120_187 = 120.0 / 187.0;
    const auto f_576_143 = 576.0 / 143.0;
    const auto f_75_143 = 75.0 / 143.0;
    const auto f_2178_96577 = 2178.0 / 96577.0;
    const auto fs_143748_2221271 = std::sqrt(143748.0 / 2221271.0);
    const auto f_10458_96577 = 10458.0 / 96577.0;
    const auto fs_873180_42204149 = std::sqrt(873180.0 / 42204149.0);
    const auto f_445_2717 = 445.0 / 2717.0;
    const auto fs_1125_51623 = std::sqrt(1125.0 / 51623.0);
    const auto f_80_3553 = 80.0 / 3553.0;
    const auto f_288_2431 = 288.0 / 2431.0;
    const auto f_2_143 = 2.0 / 143.0;
    const auto fs_218791125_2048 = std::sqrt(106831.60400390625);
    const auto fs_2679075_256 = std::sqrt(10465.13671875);
    const auto fs_111628125_512 = std::sqrt(218023.681640625);
    const auto fs_4921875_3872 = std::sqrt(4921875.0 / 3872.0);
    const auto fs_2679075_484 = std::sqrt(2679075.0 / 484.0);
    const auto fs_1378125_32 = std::sqrt(43066.40625);
    const auto fs_13395375_163592 = std::sqrt(13395375.0 / 163592.0);
    const auto fs_33075_1144 = std::sqrt(33075.0 / 1144.0);
    const auto fs_196875_968 = std::sqrt(196875.0 / 968.0);
    const auto fs_24111675_81796 = std::sqrt(24111675.0 / 81796.0);
    const auto fs_1378125_968 = std::sqrt(1378125.0 / 968.0);
    const auto fs_220172337_141052808 = std::sqrt(220172337.0 / 141052808.0);
    const auto fs_392931_319124 = std::sqrt(392931.0 / 319124.0);
    const auto fs_26790750_7382089 = std::sqrt(26790750.0 / 7382089.0);
    const auto fs_66150_51623 = std::sqrt(66150.0 / 51623.0);
    const auto fs_196875_69938 = std::sqrt(196875.0 / 69938.0);
    const auto fs_47628_20449 = std::sqrt(47628.0 / 20449.0);
    const auto fs_1378125_163592 = std::sqrt(1378125.0 / 163592.0);
    const auto fs_658845_717470533 = std::sqrt(658845.0 / 717470533.0);
    const auto fs_119790_2221271 = std::sqrt(119790.0 / 2221271.0);
    const auto fs_220172337_18654233858 = std::sqrt(220172337.0 / 18654233858.0);
    const auto fs_392931_42204149 = std::sqrt(392931.0 / 42204149.0);
    const auto fs_60750_7382089 = std::sqrt(60750.0 / 7382089.0);
    const auto fs_150_51623 = std::sqrt(150.0 / 51623.0);
    const auto fs_43750_12623809 = std::sqrt(43750.0 / 12623809.0);
    const auto fs_11907_5909761 = std::sqrt(11907.0 / 5909761.0);
    const auto fs_245_40898 = std::sqrt(245.0 / 40898.0);
    const auto fs_40186125_512 = std::sqrt(78488.525390625);
    const auto fs_3213675_512 = std::sqrt(6276.708984375);
    const auto fs_20503125_128 = std::sqrt(160180.6640625);
    const auto fs_354375_121 = std::sqrt(354375.0 / 121.0);
    const auto fs_3213675_968 = std::sqrt(3213675.0 / 968.0);
    const auto fs_253125_8 = std::sqrt(31640.625);
    const auto fs_231525_20449 = std::sqrt(231525.0 / 20449.0);
    const auto fs_55125_2288 = std::sqrt(55125.0 / 2288.0);
    const auto fs_56700_121 = std::sqrt(56700.0 / 121.0);
    const auto fs_28923075_163592 = std::sqrt(28923075.0 / 163592.0);
    const auto fs_253125_242 = std::sqrt(253125.0 / 242.0);
    const auto fs_80725491_35263202 = std::sqrt(80725491.0 / 35263202.0);
    const auto fs_43659_1356277 = std::sqrt(43659.0 / 1356277.0);
    const auto fs_3704400_7382089 = std::sqrt(3704400.0 / 7382089.0);
    const auto fs_55125_51623 = std::sqrt(55125.0 / 51623.0);
    const auto fs_226800_34969 = std::sqrt(226800.0 / 34969.0);
    const auto fs_28566_20449 = std::sqrt(28566.0 / 20449.0);
    const auto fs_253125_40898 = std::sqrt(253125.0 / 40898.0);
    const auto fs_3773385_1434941066 = std::sqrt(3773385.0 / 1434941066.0);
    const auto fs_1617165_42204149 = std::sqrt(1617165.0 / 42204149.0);
    const auto fs_161450982_9327116929 = std::sqrt(161450982.0 / 9327116929.0);
    const auto fs_174636_717470533 = std::sqrt(174636.0 / 717470533.0);
    const auto fs_8400_7382089 = std::sqrt(8400.0 / 7382089.0);
    const auto fs_125_51623 = std::sqrt(125.0 / 51623.0);
    const auto fs_100800_12623809 = std::sqrt(100800.0 / 12623809.0);
    const auto fs_14283_11819522 = std::sqrt(14283.0 / 11819522.0);
    const auto fs_90_20449 = std::sqrt(90.0 / 20449.0);
    const auto fs_5315625_128 = std::sqrt(41528.3203125);
    const auto fs_5315625_242 = std::sqrt(5315625.0 / 242.0);
    const auto fs_385875_29744 = std::sqrt(385875.0 / 29744.0);
    const auto fs_231525_2288 = std::sqrt(231525.0 / 2288.0);
    const auto fs_47840625_40898 = std::sqrt(47840625.0 / 40898.0);
    const auto fs_26413695_10850216 = std::sqrt(26413695.0 / 10850216.0);
    const auto fs_5282739_10850216 = std::sqrt(5282739.0 / 10850216.0);
    const auto fs_385875_671099 = std::sqrt(385875.0 / 671099.0);
    const auto fs_231525_51623 = std::sqrt(231525.0 / 51623.0);
    const auto fs_189000_20449 = std::sqrt(189000.0 / 20449.0);
    const auto fs_4528062_717470533 = std::sqrt(4528062.0 / 717470533.0);
    const auto fs_1006236_42204149 = std::sqrt(1006236.0 / 42204149.0);
    const auto fs_26413695_1434941066 = std::sqrt(26413695.0 / 1434941066.0);
    const auto fs_5282739_1434941066 = std::sqrt(5282739.0 / 1434941066.0);
    const auto fs_875_671099 = std::sqrt(875.0 / 671099.0);
    const auto fs_525_51623 = std::sqrt(525.0 / 51623.0);
    const auto fs_47250_5909761 = std::sqrt(47250.0 / 5909761.0);
    const auto fs_4465125_64 = std::sqrt(69767.578125);
    const auto fs_4465125_121 = std::sqrt(4465125.0 / 121.0);
    const auto fs_4862025_29744 = std::sqrt(4862025.0 / 29744.0);
    const auto fs_40186125_20449 = std::sqrt(40186125.0 / 20449.0);
    const auto fs_4584195_1356277 = std::sqrt(4584195.0 / 1356277.0);
    const auto fs_4862025_671099 = std::sqrt(4862025.0 / 671099.0);
    const auto fs_317520_20449 = std::sqrt(317520.0 / 20449.0);
    const auto fs_18783072_717470533 = std::sqrt(18783072.0 / 717470533.0);
    const auto fs_18336780_717470533 = std::sqrt(18336780.0 / 717470533.0);
    const auto fs_11025_671099 = std::sqrt(11025.0 / 671099.0);
    const auto fs_79380_5909761 = std::sqrt(79380.0 / 5909761.0);
    const auto f_4725_32 = 147.65625;
    const auto f_1215_8 = 151.875;
    const auto f_3375_16 = 210.9375;
    const auto f_3225_44 = 3225.0 / 44.0;
    const auto f_1215_11 = 1215.0 / 11.0;
    const auto f_1995_572 = 1995.0 / 572.0;
    const auto fs_33075_286 = std::sqrt(33075.0 / 286.0);
    const auto f_645_22 = 645.0 / 22.0;
    const auto f_3645_143 = 3645.0 / 143.0;
    const auto f_375_22 = 375.0 / 22.0;
    const auto f_945_442 = 945.0 / 442.0;
    const auto fs_5893965_2712554 = std::sqrt(5893965.0 / 2712554.0);
    const auto f_105_143 = 105.0 / 143.0;
    const auto fs_264600_51623 = std::sqrt(264600.0 / 51623.0);
    const auto f_645_187 = 645.0 / 187.0;
    const auto f_324_143 = 324.0 / 143.0;
    const auto f_375_286 = 375.0 / 286.0;
    const auto f_7260_96577 = 7260.0 / 96577.0;
    const auto fs_2395800_42204149 = std::sqrt(2395800.0 / 42204149.0);
    const auto f_945_5083 = 945.0 / 5083.0;
    const auto fs_11787930_717470533 = std::sqrt(11787930.0 / 717470533.0);
    const auto f_5_143 = 5.0 / 143.0;
    const auto fs_600_51623 = std::sqrt(600.0 / 51623.0);
    const auto f_430_3553 = 430.0 / 3553.0;
    const auto f_162_2431 = 162.0 / 2431.0;
    const auto fs_66976875_1024 = std::sqrt(65407.1044921875);
    const auto fs_15400125_512 = std::sqrt(30078.369140625);
    const auto fs_34171875_256 = std::sqrt(133483.88671875);
    const auto fs_1063125_1936 = std::sqrt(1063125.0 / 1936.0);
    const auto fs_15400125_968 = std::sqrt(15400125.0 / 968.0);
    const auto fs_421875_16 = std::sqrt(26367.1875);
    const auto f_105_22 = 105.0 / 22.0;
    const auto fs_77175_2288 = std::sqrt(77175.0 / 2288.0);
    const auto fs_42525_484 = std::sqrt(42525.0 / 484.0);
    const auto fs_820125_968 = std::sqrt(820125.0 / 968.0);
    const auto fs_421875_484 = std::sqrt(421875.0 / 484.0);
    const auto fs_159137055_70526404 = std::sqrt(159137055.0 / 70526404.0);
    const auto fs_3274425_2712554 = std::sqrt(3274425.0 / 2712554.0);
    const auto f_210_209 = 210.0 / 209.0;
    const auto fs_77175_51623 = std::sqrt(77175.0 / 51623.0);
    const auto fs_42525_34969 = std::sqrt(42525.0 / 34969.0);
    const auto fs_810_121 = std::sqrt(810.0 / 121.0);
    const auto fs_421875_81796 = std::sqrt(421875.0 / 81796.0);
    const auto fs_9882675_1434941066 = std::sqrt(9882675.0 / 1434941066.0);
    const auto fs_2096325_42204149 = std::sqrt(2096325.0 / 42204149.0);
    const auto fs_159137055_9327116929 = std::sqrt(159137055.0 / 9327116929.0);
    const auto fs_6548850_717470533 = std::sqrt(6548850.0 / 717470533.0);
    const auto f_10_209 = 10.0 / 209.0;
    const auto fs_175_51623 = std::sqrt(175.0 / 51623.0);
    const auto fs_18900_12623809 = std::sqrt(18900.0 / 12623809.0);
    const auto fs_405_69938 = std::sqrt(405.0 / 69938.0);
    const auto fs_75_20449 = std::sqrt(75.0 / 20449.0);
    const auto fs_13395375_128 = std::sqrt(104651.3671875);
    const auto fs_18225_8 = std::sqrt(2278.125);
    const auto fs_6251175_128 = std::sqrt(48837.3046875);
    const auto fs_6834375_32 = std::sqrt(213574.21875);
    const auto fs_145800_121 = std::sqrt(145800.0 / 121.0);
    const auto fs_6251175_242 = std::sqrt(6251175.0 / 242.0);
    const auto fs_84375_2 = std::sqrt(42187.5);
    const auto fs_27860175_327184 = std::sqrt(27860175.0 / 327184.0);
    const auto fs_385875_59488 = std::sqrt(385875.0 / 59488.0);
    const auto fs_1312200_20449 = std::sqrt(1312200.0 / 20449.0);
    const auto fs_56260575_40898 = std::sqrt(56260575.0 / 40898.0);
    const auto fs_168750_121 = std::sqrt(168750.0 / 121.0);
    const auto fs_130977_97682 = std::sqrt(130977.0 / 97682.0);
    const auto fs_1178793_10850216 = std::sqrt(1178793.0 / 10850216.0);
    const auto fs_77175_20449 = std::sqrt(77175.0 / 20449.0);
    const auto fs_385875_1342198 = std::sqrt(385875.0 / 1342198.0);
    const auto fs_10368_20449 = std::sqrt(10368.0 / 20449.0);
    const auto fs_222264_20449 = std::sqrt(222264.0 / 20449.0);
    const auto fs_168750_20449 = std::sqrt(168750.0 / 20449.0);
    const auto fs_10062360_717470533 = std::sqrt(10062360.0 / 717470533.0);
    const auto fs_26832960_717470533 = std::sqrt(26832960.0 / 717470533.0);
    const auto fs_261954_25836889 = std::sqrt(261954.0 / 25836889.0);
    const auto fs_1178793_1434941066 = std::sqrt(1178793.0 / 1434941066.0);
    const auto fs_175_20449 = std::sqrt(175.0 / 20449.0);
    const auto fs_875_1342198 = std::sqrt(875.0 / 1342198.0);
    const auto fs_2592_5909761 = std::sqrt(2592.0 / 5909761.0);
    const auto fs_55566_5909761 = std::sqrt(55566.0 / 5909761.0);
    const auto fs_120_20449 = std::sqrt(120.0 / 20449.0);
    const auto fs_2679075_64 = std::sqrt(41860.546875);
    const auto fs_2679075_121 = std::sqrt(2679075.0 / 121.0);
    const auto fs_540225_3718 = std::sqrt(540225.0 / 3718.0);
    const auto fs_24111675_20449 = std::sqrt(24111675.0 / 20449.0);
    const auto fs_2750517_5425108 = std::sqrt(2750517.0 / 5425108.0);
    const auto fs_4321800_671099 = std::sqrt(4321800.0 / 671099.0);
    const auto fs_190512_20449 = std::sqrt(190512.0 / 20449.0);
    const auto fs_35218260_717470533 = std::sqrt(35218260.0 / 717470533.0);
    const auto fs_2750517_717470533 = std::sqrt(2750517.0 / 717470533.0);
    const auto fs_9800_671099 = std::sqrt(9800.0 / 671099.0);
    const auto fs_47628_5909761 = std::sqrt(47628.0 / 5909761.0);
    const auto f_4725_16 = 295.3125;
    const auto f_495_16 = 30.9375;
    const auto fs_3472875_64 = std::sqrt(54263.671875);
    const auto f_3375_8 = 421.875;
    const auto f_45_2 = 22.5;
    const auto fs_3472875_121 = std::sqrt(3472875.0 / 121.0);
    const auto f_375_2 = 187.5;
    const auto f_7455_572 = 7455.0 / 572.0;
    const auto fs_694575_7436 = std::sqrt(694575.0 / 7436.0);
    const auto f_135_26 = 135.0 / 26.0;
    const auto fs_31255875_20449 = std::sqrt(31255875.0 / 20449.0);
    const auto f_375_11 = 375.0 / 11.0;
    const auto f_6615_4199 = 6615.0 / 4199.0;
    const auto fs_2619540_1356277 = std::sqrt(2619540.0 / 1356277.0);
    const auto f_7455_2717 = 7455.0 / 2717.0;
    const auto fs_2778300_671099 = std::sqrt(2778300.0 / 671099.0);
    const auto f_6_13 = 6.0 / 13.0;
    const auto fs_246960_20449 = std::sqrt(246960.0 / 20449.0);
    const auto f_375_143 = 375.0 / 143.0;
    const auto f_16335_96577 = 16335.0 / 96577.0;
    const auto fs_37733850_717470533 = std::sqrt(37733850.0 / 717470533.0);
    const auto f_13230_96577 = 13230.0 / 96577.0;
    const auto fs_10478160_717470533 = std::sqrt(10478160.0 / 717470533.0);
    const auto f_355_2717 = 355.0 / 2717.0;
    const auto fs_6300_671099 = std::sqrt(6300.0 / 671099.0);
    const auto f_3_221 = 3.0 / 221.0;
    const auto fs_61740_5909761 = std::sqrt(61740.0 / 5909761.0);
    const auto f_10_143 = 10.0 / 143.0;
    const auto fs_13395375_512 = std::sqrt(26162.841796875);
    const auto f_2385_16 = 149.0625;
    const auto fs_694575_256 = std::sqrt(2713.18359375);
    const auto fs_6834375_128 = std::sqrt(53393.5546875);
    const auto fs_590625_242 = std::sqrt(590625.0 / 242.0);
    const auto f_2385_22 = 2385.0 / 22.0;
    const auto fs_694575_484 = std::sqrt(694575.0 / 484.0);
    const auto fs_84375_8 = std::sqrt(10546.875);
    const auto fs_40186125_654368 = std::sqrt(40186125.0 / 654368.0);
    const auto fs_2083725_59488 = std::sqrt(2083725.0 / 59488.0);
    const auto fs_47250_121 = std::sqrt(47250.0 / 121.0);
    const auto f_7155_286 = 7155.0 / 286.0;
    const auto fs_6251175_81796 = std::sqrt(6251175.0 / 81796.0);
    const auto fs_84375_242 = std::sqrt(84375.0 / 242.0);
    const auto fs_3143448_17631601 = std::sqrt(3143448.0 / 17631601.0);
    const auto fs_6417873_5425108 = std::sqrt(6417873.0 / 5425108.0);
    const auto fs_40186125_14764178 = std::sqrt(40186125.0 / 14764178.0);
    const auto fs_2083725_1342198 = std::sqrt(2083725.0 / 1342198.0);
    const auto fs_189000_34969 = std::sqrt(189000.0 / 34969.0);
    const auto f_318_143 = 318.0 / 143.0;
    const auto fs_12348_20449 = std::sqrt(12348.0 / 20449.0);
    const auto fs_84375_40898 = std::sqrt(84375.0 / 40898.0);
    const auto fs_17788815_717470533 = std::sqrt(17788815.0 / 717470533.0);
    const auto fs_33960465_717470533 = std::sqrt(33960465.0 / 717470533.0);
    const auto fs_12573792_9327116929 = std::sqrt(12573792.0 / 9327116929.0);
    const auto fs_6417873_717470533 = std::sqrt(6417873.0 / 717470533.0);
    const auto fs_91125_14764178 = std::sqrt(91125.0 / 14764178.0);
    const auto fs_4725_1342198 = std::sqrt(4725.0 / 1342198.0);
    const auto fs_84000_12623809 = std::sqrt(84000.0 / 12623809.0);
    const auto f_159_2431 = 159.0 / 2431.0;
    const auto fs_3087_5909761 = std::sqrt(3087.0 / 5909761.0);
    const auto fs_30_20449 = std::sqrt(30.0 / 20449.0);
    const auto fs_31255875_128 = std::sqrt(244186.5234375);
    const auto fs_5145525_128 = std::sqrt(40199.4140625);
    const auto fs_15946875_32 = std::sqrt(498339.84375);
    const auto fs_42525_2 = std::sqrt(21262.5);
    const auto fs_196875_2 = std::sqrt(98437.5);
    const auto fs_1620675_327184 = std::sqrt(1620675.0 / 327184.0);
    const auto fs_382725_338 = std::sqrt(382725.0 / 338.0);
    const auto fs_393750_121 = std::sqrt(393750.0 / 121.0);
    const auto fs_5501034_17631601 = std::sqrt(5501034.0 / 17631601.0);
    const auto fs_1620675_7382089 = std::sqrt(1620675.0 / 7382089.0);
    const auto fs_1512_169 = std::sqrt(1512.0 / 169.0);
    const auto fs_393750_20449 = std::sqrt(393750.0 / 20449.0);
    const auto fs_52827390_717470533 = std::sqrt(52827390.0 / 717470533.0);
    const auto fs_22004136_9327116929 = std::sqrt(22004136.0 / 9327116929.0);
    const auto fs_3675_7382089 = std::sqrt(3675.0 / 7382089.0);
    const auto fs_378_48841 = std::sqrt(378.0 / 48841.0);
    const auto fs_280_20449 = std::sqrt(280.0 / 20449.0);
    const auto f_12285_32 = 383.90625;
    const auto fs_131274675_1024 = std::sqrt(128197.9248046875);
    const auto f_180 = 180.0;
    const auto fs_496125_16 = std::sqrt(31007.8125);
    const auto f_8775_16 = 548.4375;
    const auto fs_66976875_256 = std::sqrt(261628.41796875);
    const auto f_1440_11 = 1440.0 / 11.0;
    const auto fs_1984500_121 = std::sqrt(1984500.0 / 121.0);
    const auto f_975_4 = 243.75;
    const auto fs_826875_16 = std::sqrt(51679.6875);
    const auto f_525_286 = 525.0 / 286.0;
    const auto fs_3472875_40898 = std::sqrt(3472875.0 / 40898.0);
    const auto f_4320_143 = 4320.0 / 143.0;
    const auto fs_17860500_20449 = std::sqrt(17860500.0 / 20449.0);
    const auto f_975_22 = 975.0 / 22.0;
    const auto fs_826875_484 = std::sqrt(826875.0 / 484.0);
    const auto f_189_323 = 189.0 / 323.0;
    const auto fs_32089365_17631601 = std::sqrt(32089365.0 / 17631601.0);
    const auto f_1050_2717 = 1050.0 / 2717.0;
    const auto fs_27783000_7382089 = std::sqrt(27783000.0 / 7382089.0);
    const auto f_300_187 = 300.0 / 187.0;
    const auto f_384_143 = 384.0 / 143.0;
    const auto fs_141120_20449 = std::sqrt(141120.0 / 20449.0);
    const auto f_75_22 = 75.0 / 22.0;
    const auto fs_826875_81796 = std::sqrt(826875.0 / 81796.0);
    const auto f_26136_96577 = 26136.0 / 96577.0;
    const auto fs_36224496_717470533 = std::sqrt(36224496.0 / 717470533.0);
    const auto f_378_7429 = 378.0 / 7429.0;
    const auto fs_128357460_9327116929 = std::sqrt(128357460.0 / 9327116929.0);
    const auto f_50_2717 = 50.0 / 2717.0;
    const auto fs_63000_7382089 = std::sqrt(63000.0 / 7382089.0);
    const auto f_200_3553 = 200.0 / 3553.0;
    const auto f_192_2431 = 192.0 / 2431.0;
    const auto fs_35280_5909761 = std::sqrt(35280.0 / 5909761.0);
    const auto f_1_11 = 1.0 / 11.0;
    const auto fs_147_20449 = std::sqrt(147.0 / 20449.0);
    const auto fs_6251175_1024 = std::sqrt(6104.6630859375);
    const auto fs_212625_32 = std::sqrt(6644.53125);
    const auto fs_3189375_256 = std::sqrt(12458.49609375);
    const auto fs_425250_121 = std::sqrt(425250.0 / 121.0);
    const auto fs_5788125_81796 = std::sqrt(5788125.0 / 81796.0);
    const auto fs_3827250_20449 = std::sqrt(3827250.0 / 20449.0);
    const auto fs_39375_484 = std::sqrt(39375.0 / 484.0);
    const auto fs_41257755_17631601 = std::sqrt(41257755.0 / 17631601.0);
    const auto fs_23152500_7382089 = std::sqrt(23152500.0 / 7382089.0);
    const auto fs_30240_20449 = std::sqrt(30240.0 / 20449.0);
    const auto fs_39375_81796 = std::sqrt(39375.0 / 81796.0);
    const auto fs_66411576_717470533 = std::sqrt(66411576.0 / 717470533.0);
    const auto fs_165031020_9327116929 = std::sqrt(165031020.0 / 9327116929.0);
    const auto fs_52500_7382089 = std::sqrt(52500.0 / 7382089.0);
    const auto fs_7560_5909761 = std::sqrt(7560.0 / 5909761.0);
    const auto fs_7_20449 = std::sqrt(7.0 / 20449.0);
    const auto f_6615_16 = 413.4375;
    const auto f_945_4 = 236.25;
    const auto f_4725_8 = 590.625;
    const auto f_750_11 = 750.0 / 11.0;
    const auto f_1890_11 = 1890.0 / 11.0;
    const auto f_525_2 = 262.5;
    const auto f_3675_286 = 3675.0 / 286.0;
    const auto f_300_11 = 300.0 / 11.0;
    const auto f_5670_143 = 5670.0 / 143.0;
    const auto f_525_11 = 525.0 / 11.0;
    const auto f_7938_4199 = 7938.0 / 4199.0;
    const auto f_7350_2717 = 7350.0 / 2717.0;
    const auto f_600_187 = 600.0 / 187.0;
    const auto f_504_143 = 504.0 / 143.0;
    const auto f_525_143 = 525.0 / 143.0;
    const auto f_30492_96577 = 30492.0 / 96577.0;
    const auto f_15876_96577 = 15876.0 / 96577.0;
    const auto f_350_2717 = 350.0 / 2717.0;
    const auto f_400_3553 = 400.0 / 3553.0;
    const auto f_252_2431 = 252.0 / 2431.0;
    const auto f_14_143 = 14.0 / 143.0;

    // NOTE: the rows are formed in 82 loops, as the vectorizer runs out of
    // registers with all 91 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph2_0, ph4_0, ph6_0, ph8_0, ph10_0, ph12_0, ph12_p12, ab_2, pc_0 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;
        const auto r_10 = r_8 * r_2;
        const auto r_12 = r_10 * r_2;

        const auto h2_0 = ph2_0[k];
        const auto h4_0 = ph4_0[k];
        const auto h6_0 = ph6_0[k];
        const auto h8_0 = ph8_0[k];
        const auto h10_0 = ph10_0[k];
        const auto h12_0 = ph12_0[k];
        const auto h12_p12 = ph12_p12[k];

        pc_0[k] = e_0 * (f_10395_16 * h2_0 - f_10395_16 * r_2) + e_1 * (f_4455_16 * h4_0 - f_7425_8 * r_2 * h2_0 + f_10395_16 * r_4) + e_2 * (f_75_2 * h6_0 - f_405_2 * r_2 * h4_0 + f_825_2 * r_4 * h2_0 - f_495_2 * r_6) + e_3 * (f_105_52 * h8_0 - f_15 * r_2 * h6_0 + f_1215_26 * r_4 * h4_0 - f_75 * r_6 * h2_0 + f_165_4 * r_8) + e_4 * (f_189_4199 * h10_0 - f_105_247 * r_2 * h8_0 + f_30_17 * r_4 * h6_0 - f_54_13 * r_6 * h4_0 + f_75_13 * r_8 * h2_0 - f_3 * r_10) + e_5 * (f_33_96577 * h12_0 - fs_15246_96577 * h12_p12 - f_378_96577 * r_2 * h10_0 + f_5_247 * r_4 * h8_0 - f_20_323 * r_6 * h6_0 + f_27_221 * r_8 * h4_0 - f_2_13 * r_10 * h2_0 + f_1_13 * r_12) + f_10395_64 * e_6;
    }

    // NOTE: the rows are formed in 82 loops, as the vectorizer runs out of
    // registers with all 91 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph2_p1, ph4_p1, ph6_p1, ph8_p1, ph10_p1, ph12_p1, ph12_p11, ab_2, pc_1 : simd::cache_line_size())
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

        const auto h2_p1 = ph2_p1[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h10_p1 = ph10_p1[k];
        const auto h12_p1 = ph12_p1[k];
        const auto h12_p11 = ph12_p11[k];

        pc_1[k] = -f_10395_32 * e_0 * h2_p1 + e_1 * (-fs_33078375_512 * h4_p1 + f_7425_16 * r_2 * h2_p1) + e_2 * (-fs_39375_16 * h6_p1 + fs_273375_8 * r_2 * h4_p1 - f_825_4 * r_4 * h2_p1) + e_3 * (-fs_33075_2704 * h8_p1 + fs_1575_4 * r_2 * h6_p1 - fs_2460375_1352 * r_4 * h4_p1 + f_75_2 * r_6 * h2_p1) + e_4 * (-fs_654885_70526404 * h10_p1 + fs_33075_61009 * r_2 * h8_p1 - fs_1575_289 * r_4 * h6_p1 + fs_2430_169 * r_6 * h4_p1 - f_75_26 * r_8 * h2_p1) + e_5 * (-fs_1089_1434941066 * h12_p1 - fs_7623_96577 * h12_p11 + fs_654885_9327116929 * r_2 * h10_p1 - fs_75_61009 * r_4 * h8_p1 + fs_700_104329 * r_6 * h6_p1 - fs_1215_97682 * r_8 * h4_p1 + f_1_13 * r_10 * h2_p1);
    }

    // NOTE: the rows are formed in 82 loops, as the vectorizer runs out of
    // registers with all 91 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph2_p2, ph4_p2, ph6_p2, ph8_p2, ph10_p2, ph10_p10, ph12_p2, ph12_p10, ab_2, pc_2 : simd::cache_line_size())
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

        const auto h2_p2 = ph2_p2[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h8_p2 = ph8_p2[k];
        const auto h10_p2 = ph10_p2[k];
        const auto h10_p10 = ph10_p10[k];
        const auto h12_p2 = ph12_p2[k];
        const auto h12_p10 = ph12_p10[k];

        pc_2[k] = fs_9823275_512 * e_0 * h2_p2 + e_1 * (fs_27064125_512 * h4_p2 - fs_5011875_128 * r_2 * h2_p2) + e_2 * (fs_196875_44 * h6_p2 - fs_2460375_88 * r_2 * h4_p2 + fs_61875_8 * r_4 * h2_p2) + e_3 * (fs_1157625_29744 * h8_p2 - fs_7875_11 * r_2 * h6_p2 + fs_22143375_14872 * r_4 * h4_p2 - fs_5625_22 * r_6 * h2_p2) + e_4 * (fs_1607445_35263202 * h10_p2 - fs_11907_4199 * h10_p10 - fs_1157625_671099 * r_2 * h8_p2 + fs_31500_3179 * r_4 * h6_p2 - fs_21870_1859 * r_6 * h4_p2 + fs_5625_3718 * r_8 * h2_p2) + e_5 * (fs_7623_1434941066 * h12_p2 - fs_83853_2221271 * h12_p10 - fs_3214890_9327116929 * r_2 * h10_p2 + fs_47628_2221271 * r_2 * h10_p10 + fs_2625_671099 * r_4 * h8_p2 - fs_14000_1147619 * r_6 * h6_p2 + fs_10935_1074502 * r_8 * h4_p2 - fs_2_1859 * r_10 * h2_p2);
    }

    // NOTE: the rows are formed in 82 loops, as the vectorizer runs out of
    // registers with all 91 of them in one.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph4_p3, ph6_p3, ph8_p3, ph10_p3, ph10_p9, ph12_p3, ph12_p9, ab_2, pc_3 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h4_p3 = ph4_p3[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h8_p3 = ph8_p3[k];
        const auto h10_p3 = ph10_p3[k];
        const auto h10_p9 = ph10_p9[k];
        const auto h12_p3 = ph12_p3[k];
        const auto h12_p9 = ph12_p9[k];

        pc_3[k] = -fs_12629925_512 * e_1 * h4_p3 + e_2 * (-fs_118125_22 * h6_p3 + fs_1148175_88 * r_2 * h4_p3) + e_3 * (-fs_231525_2704 * h8_p3 + fs_9450_11 * r_2 * h6_p3 - fs_10333575_14872 * r_4 * h4_p3) + e_4 * (-fs_214326_1356277 * h10_p3 - fs_35721_8398 * h10_p9 + fs_231525_61009 * r_2 * h8_p3 - fs_37800_3179 * r_4 * h6_p3 + fs_10206_1859 * r_6 * h4_p3) + e_5 * (-fs_38115_1434941066 * h12_p3 - fs_38115_2221271 * h12_p9 + fs_857304_717470533 * r_2 * h10_p3 + fs_71442_2221271 * r_2 * h10_p9 - fs_525_61009 * r_4 * h8_p3 + fs_16800_1147619 * r_6 * h6_p3 - fs_5103_1074502 * r_8 * h4_p3);
    }

    // NOTE: the rows are formed in 82 loops, as the vectorizer runs out of
    // registers with all 91 of them in one.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph4_p4, ph6_p4, ph8_p4, ph8_p8, ph10_p4, ph10_p8, ph12_p4, ph12_p8, ab_2, pc_4 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h4_p4 = ph4_p4[k];
        const auto h6_p4 = ph6_p4[k];
        const auto h8_p4 = ph8_p4[k];
        const auto h8_p8 = ph8_p8[k];
        const auto h10_p4 = ph10_p4[k];
        const auto h10_p8 = ph10_p8[k];
        const auto h12_p4 = ph12_p4[k];
        const auto h12_p8 = ph12_p8[k];

        pc_4[k] = fs_1403325_256 * e_1 * h4_p4 + e_2 * (fs_196875_44 * h6_p4 - fs_127575_44 * r_2 * h4_p4) + e_3 * (fs_385875_2704 * h8_p4 - fs_11025_208 * h8_p8 - fs_7875_11 * r_2 * h6_p4 + fs_1148175_7436 * r_4 * h4_p4) + e_4 * (fs_583443_1356277 * h10_p4 - fs_321489_79781 * h10_p8 - fs_385875_61009 * r_2 * h8_p4 + fs_11025_4693 * r_2 * h8_p8 + fs_31500_3179 * r_4 * h6_p4 - fs_2268_1859 * r_6 * h4_p4) + e_5 * (fs_76230_717470533 * h12_p4 - fs_16335_2221271 * h12_p8 - fs_2333772_717470533 * r_2 * h10_p4 + fs_1285956_42204149 * r_2 * h10_p8 + fs_875_61009 * r_4 * h8_p4 - fs_25_4693 * r_4 * h8_p8 - fs_14000_1147619 * r_6 * h6_p4 + fs_567_537251 * r_8 * h4_p4);
    }

    // NOTE: the rows are formed in 82 loops, as the vectorizer runs out of
    // registers with all 91 of them in one.

#pragma omp simd aligned(pe_2, pe_3, pe_4, pe_5, ph6_m6, ph6_p5, ph8_m6, ph8_p5, ph8_p7, ph10_m6, ph10_p5, ph10_p7, ph12_m6, ph12_p5, ph12_p7, ab_2, pc_5, pc_6 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h6_m6 = ph6_m6[k];
        const auto h6_p5 = ph6_p5[k];
        const auto h8_m6 = ph8_m6[k];
        const auto h8_p5 = ph8_p5[k];
        const auto h8_p7 = ph8_p7[k];
        const auto h10_m6 = ph10_m6[k];
        const auto h10_p5 = ph10_p5[k];
        const auto h10_p7 = ph10_p7[k];
        const auto h12_m6 = ph12_m6[k];
        const auto h12_p5 = ph12_p5[k];
        const auto h12_p7 = ph12_p7[k];

        pc_5[k] = -fs_39375_16 * e_2 * h6_p5 + e_3 * (-fs_77175_416 * h8_p5 - fs_55125_416 * h8_p7 + fs_1575_4 * r_2 * h6_p5) + e_4 * (-fs_5250987_5425108 * h10_p5 - fs_238140_79781 * h10_p7 + fs_77175_9386 * r_2 * h8_p5 + fs_55125_9386 * r_2 * h8_p7 - fs_1575_289 * r_4 * h6_p5) + e_5 * (-fs_15246_42204149 * h12_p5 - fs_6534_2221271 * h12_p7 + fs_5250987_717470533 * r_2 * h10_p5 + fs_952560_42204149 * r_2 * h10_p7 - fs_175_9386 * r_4 * h8_p5 - fs_125_9386 * r_4 * h8_p7 + fs_700_104329 * r_6 * h6_p5);

        pc_6[k] = f_75_2 * e_2 * h6_m6 + e_3 * (fs_77175_208 * h8_m6 - f_15 * r_2 * h6_m6) + e_4 * (fs_5000940_1356277 * h10_m6 - fs_77175_4693 * r_2 * h8_m6 + f_30_17 * r_4 * h6_m6) + e_5 * (fs_91476_42204149 * h12_m6 - fs_20003760_717470533 * r_2 * h10_m6 + fs_175_4693 * r_4 * h8_m6 - f_20_323 * r_6 * h6_m6);
    }

    // NOTE: the rows are formed in 82 loops, as the vectorizer runs out of
    // registers with all 91 of them in one.

#pragma omp simd aligned(pe_2, pe_3, pe_4, pe_5, ph6_m5, ph8_m7, ph8_m5, ph10_m7, ph10_m5, ph12_m7, ph12_m5, ab_2, pc_7 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h6_m5 = ph6_m5[k];
        const auto h8_m7 = ph8_m7[k];
        const auto h8_m5 = ph8_m5[k];
        const auto h10_m7 = ph10_m7[k];
        const auto h10_m5 = ph10_m5[k];
        const auto h12_m7 = ph12_m7[k];
        const auto h12_m5 = ph12_m5[k];

        pc_7[k] = -fs_39375_16 * e_2 * h6_m5 + e_3 * (fs_55125_416 * h8_m7 - fs_77175_416 * h8_m5 + fs_1575_4 * r_2 * h6_m5) + e_4 * (fs_238140_79781 * h10_m7 - fs_5250987_5425108 * h10_m5 - fs_55125_9386 * r_2 * h8_m7 + fs_77175_9386 * r_2 * h8_m5 - fs_1575_289 * r_4 * h6_m5) + e_5 * (fs_6534_2221271 * h12_m7 - fs_15246_42204149 * h12_m5 - fs_952560_42204149 * r_2 * h10_m7 + fs_5250987_717470533 * r_2 * h10_m5 + fs_125_9386 * r_4 * h8_m7 - fs_175_9386 * r_4 * h8_m5 + fs_700_104329 * r_6 * h6_m5);
    }

    // NOTE: the rows are formed in 82 loops, as the vectorizer runs out of
    // registers with all 91 of them in one.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph4_m4, ph6_m4, ph8_m8, ph8_m4, ph10_m8, ph10_m4, ph12_m8, ph12_m4, ab_2, pc_8 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h4_m4 = ph4_m4[k];
        const auto h6_m4 = ph6_m4[k];
        const auto h8_m8 = ph8_m8[k];
        const auto h8_m4 = ph8_m4[k];
        const auto h10_m8 = ph10_m8[k];
        const auto h10_m4 = ph10_m4[k];
        const auto h12_m8 = ph12_m8[k];
        const auto h12_m4 = ph12_m4[k];

        pc_8[k] = fs_1403325_256 * e_1 * h4_m4 + e_2 * (fs_196875_44 * h6_m4 - fs_127575_44 * r_2 * h4_m4) + e_3 * (fs_11025_208 * h8_m8 + fs_385875_2704 * h8_m4 - fs_7875_11 * r_2 * h6_m4 + fs_1148175_7436 * r_4 * h4_m4) + e_4 * (fs_321489_79781 * h10_m8 + fs_583443_1356277 * h10_m4 - fs_11025_4693 * r_2 * h8_m8 - fs_385875_61009 * r_2 * h8_m4 + fs_31500_3179 * r_4 * h6_m4 - fs_2268_1859 * r_6 * h4_m4) + e_5 * (fs_16335_2221271 * h12_m8 + fs_76230_717470533 * h12_m4 - fs_1285956_42204149 * r_2 * h10_m8 - fs_2333772_717470533 * r_2 * h10_m4 + fs_25_4693 * r_4 * h8_m8 + fs_875_61009 * r_4 * h8_m4 - fs_14000_1147619 * r_6 * h6_m4 + fs_567_537251 * r_8 * h4_m4);
    }

    // NOTE: the rows are formed in 82 loops, as the vectorizer runs out of
    // registers with all 91 of them in one.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph4_m3, ph6_m3, ph8_m3, ph10_m9, ph10_m3, ph12_m9, ph12_m3, ab_2, pc_9 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h4_m3 = ph4_m3[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h8_m3 = ph8_m3[k];
        const auto h10_m9 = ph10_m9[k];
        const auto h10_m3 = ph10_m3[k];
        const auto h12_m9 = ph12_m9[k];
        const auto h12_m3 = ph12_m3[k];

        pc_9[k] = -fs_12629925_512 * e_1 * h4_m3 + e_2 * (-fs_118125_22 * h6_m3 + fs_1148175_88 * r_2 * h4_m3) + e_3 * (-fs_231525_2704 * h8_m3 + fs_9450_11 * r_2 * h6_m3 - fs_10333575_14872 * r_4 * h4_m3) + e_4 * (fs_35721_8398 * h10_m9 - fs_214326_1356277 * h10_m3 + fs_231525_61009 * r_2 * h8_m3 - fs_37800_3179 * r_4 * h6_m3 + fs_10206_1859 * r_6 * h4_m3) + e_5 * (fs_38115_2221271 * h12_m9 - fs_38115_1434941066 * h12_m3 - fs_71442_2221271 * r_2 * h10_m9 + fs_857304_717470533 * r_2 * h10_m3 - fs_525_61009 * r_4 * h8_m3 + fs_16800_1147619 * r_6 * h6_m3 - fs_5103_1074502 * r_8 * h4_m3);
    }

    // NOTE: the rows are formed in 82 loops, as the vectorizer runs out of
    // registers with all 91 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph2_m2, ph4_m2, ph6_m2, ph8_m2, ph10_m10, ph10_m2, ph12_m10, ph12_m2, ab_2, pc_10 : simd::cache_line_size())
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

        const auto h2_m2 = ph2_m2[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h8_m2 = ph8_m2[k];
        const auto h10_m10 = ph10_m10[k];
        const auto h10_m2 = ph10_m2[k];
        const auto h12_m10 = ph12_m10[k];
        const auto h12_m2 = ph12_m2[k];

        pc_10[k] = fs_9823275_512 * e_0 * h2_m2 + e_1 * (fs_27064125_512 * h4_m2 - fs_5011875_128 * r_2 * h2_m2) + e_2 * (fs_196875_44 * h6_m2 - fs_2460375_88 * r_2 * h4_m2 + fs_61875_8 * r_4 * h2_m2) + e_3 * (fs_1157625_29744 * h8_m2 - fs_7875_11 * r_2 * h6_m2 + fs_22143375_14872 * r_4 * h4_m2 - fs_5625_22 * r_6 * h2_m2) + e_4 * (fs_11907_4199 * h10_m10 + fs_1607445_35263202 * h10_m2 - fs_1157625_671099 * r_2 * h8_m2 + fs_31500_3179 * r_4 * h6_m2 - fs_21870_1859 * r_6 * h4_m2 + fs_5625_3718 * r_8 * h2_m2) + e_5 * (fs_83853_2221271 * h12_m10 + fs_7623_1434941066 * h12_m2 - fs_47628_2221271 * r_2 * h10_m10 - fs_3214890_9327116929 * r_2 * h10_m2 + fs_2625_671099 * r_4 * h8_m2 - fs_14000_1147619 * r_6 * h6_m2 + fs_10935_1074502 * r_8 * h4_m2 - fs_2_1859 * r_10 * h2_m2);
    }

    // NOTE: the rows are formed in 82 loops, as the vectorizer runs out of
    // registers with all 91 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph2_m1, ph4_m1, ph6_m1, ph8_m1, ph10_m1, ph12_m12, ph12_m11, ph12_m1, ab_2, pc_11, pc_12 : simd::cache_line_size())
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

        const auto h2_m1 = ph2_m1[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h8_m1 = ph8_m1[k];
        const auto h10_m1 = ph10_m1[k];
        const auto h12_m12 = ph12_m12[k];
        const auto h12_m11 = ph12_m11[k];
        const auto h12_m1 = ph12_m1[k];

        pc_11[k] = -f_10395_32 * e_0 * h2_m1 + e_1 * (-fs_33078375_512 * h4_m1 + f_7425_16 * r_2 * h2_m1) + e_2 * (-fs_39375_16 * h6_m1 + fs_273375_8 * r_2 * h4_m1 - f_825_4 * r_4 * h2_m1) + e_3 * (-fs_33075_2704 * h8_m1 + fs_1575_4 * r_2 * h6_m1 - fs_2460375_1352 * r_4 * h4_m1 + f_75_2 * r_6 * h2_m1) + e_4 * (-fs_654885_70526404 * h10_m1 + fs_33075_61009 * r_2 * h8_m1 - fs_1575_289 * r_4 * h6_m1 + fs_2430_169 * r_6 * h4_m1 - f_75_26 * r_8 * h2_m1) + e_5 * (fs_7623_96577 * h12_m11 - fs_1089_1434941066 * h12_m1 + fs_654885_9327116929 * r_2 * h10_m1 - fs_75_61009 * r_4 * h8_m1 + fs_700_104329 * r_6 * h6_m1 - fs_1215_97682 * r_8 * h4_m1 + f_1_13 * r_10 * h2_m1);

        pc_12[k] = fs_15246_96577 * e_5 * h12_m12;
    }

    // NOTE: the rows are formed in 82 loops, as the vectorizer runs out of
    // registers with all 91 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph2_0, ph4_0, ph6_0, ph8_0, ph10_0, ph10_p10, ph12_0, ph12_p10, ab_2, pc_13 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;
        const auto r_10 = r_8 * r_2;
        const auto r_12 = r_10 * r_2;

        const auto h2_0 = ph2_0[k];
        const auto h4_0 = ph4_0[k];
        const auto h6_0 = ph6_0[k];
        const auto h8_0 = ph8_0[k];
        const auto h10_0 = ph10_0[k];
        const auto h10_p10 = ph10_p10[k];
        const auto h12_0 = ph12_0[k];
        const auto h12_p10 = ph12_p10[k];

        pc_13[k] = e_0 * (f_10395_32 * h2_0 - f_10395_16 * r_2) + e_1 * (-f_1485_8 * h4_0 - f_7425_16 * r_2 * h2_0 + f_10395_16 * r_4) + e_2 * (-f_375_4 * h6_0 + f_135 * r_2 * h4_0 + f_825_4 * r_4 * h2_0 - f_495_2 * r_6) + e_3 * (-f_525_52 * h8_0 + f_75_2 * r_2 * h6_0 - f_405_13 * r_4 * h4_0 - f_75_2 * r_6 * h2_0 + f_165_4 * r_8) + e_4 * (-f_3087_8398 * h10_0 + fs_43659_8398 * h10_p10 + f_525_247 * r_2 * h8_0 - f_75_17 * r_4 * h6_0 + f_36_13 * r_6 * h4_0 + f_75_26 * r_8 * h2_0 - f_3 * r_10) + e_5 * (-f_396_96577 * h12_0 - fs_182952_2221271 * h12_p10 + f_3087_96577 * r_2 * h10_0 - fs_87318_2221271 * r_2 * h10_p10 - f_25_247 * r_4 * h8_0 + f_50_323 * r_6 * h6_0 - f_18_221 * r_8 * h4_0 - f_1_13 * r_10 * h2_0 + f_1_13 * r_12) + f_10395_64 * e_6;
    }

    // NOTE: the rows are formed in 82 loops, as the vectorizer runs out of
    // registers with all 91 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph2_p1, ph4_p1, ph6_p1, ph8_p1, ph10_p1, ph10_p9, ph12_p1, ph12_p9, ab_2, pc_14 : simd::cache_line_size())
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

        const auto h2_p1 = ph2_p1[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h10_p1 = ph10_p1[k];
        const auto h10_p9 = ph10_p9[k];
        const auto h12_p1 = ph12_p1[k];
        const auto h12_p9 = ph12_p9[k];

        pc_14[k] = -fs_265228425_2048 * e_0 * h2_p1 + e_1 * (-fs_1002375_256 * h4_p1 + fs_135320625_512 * r_2 * h2_p1) + e_2 * (fs_1063125_352 * h6_p1 + fs_91125_44 * r_2 * h4_p1 - fs_1670625_32 * r_4 * h2_p1) + e_3 * (fs_198450_1859 * h8_p1 - fs_42525_88 * r_2 * h6_p1 - fs_820125_7436 * r_4 * h4_p1 + fs_151875_88 * r_6 * h2_p1) + e_4 * (fs_36693405_141052808 * h10_p1 + fs_19845_16796 * h10_p9 - fs_3175200_671099 * r_2 * h8_p1 + fs_42525_6358 * r_4 * h6_p1 + fs_1620_1859 * r_6 * h4_p1 - fs_151875_14872 * r_8 * h2_p1) + e_5 * (fs_35937_717470533 * h12_p1 - fs_137214_2221271 * h12_p9 - fs_36693405_18654233858 * r_2 * h10_p1 - fs_19845_2221271 * r_2 * h10_p9 + fs_7200_671099 * r_4 * h8_p1 - fs_9450_1147619 * r_6 * h6_p1 - fs_405_537251 * r_8 * h4_p1 + fs_27_3718 * r_10 * h2_p1);
    }

    // NOTE: the rows are formed in 82 loops, as the vectorizer runs out of
    // registers with all 91 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph2_p2, ph4_p2, ph6_p2, ph8_p2, ph8_p8, ph10_p2, ph10_p8, ph12_p2, ph12_p8, ab_2, pc_15 : simd::cache_line_size())
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

        const auto h2_p2 = ph2_p2[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h8_p2 = ph8_p2[k];
        const auto h8_p8 = ph8_p8[k];
        const auto h10_p2 = ph10_p2[k];
        const auto h10_p8 = ph10_p8[k];
        const auto h12_p2 = ph12_p2[k];
        const auto h12_p8 = ph12_p8[k];

        pc_15[k] = fs_49116375_1024 * e_0 * h2_p2 + e_1 * (fs_601425_16 * h4_p2 - fs_25059375_256 * r_2 * h2_p2) + e_2 * (-fs_39375_88 * h6_p2 - fs_218700_11 * r_2 * h4_p2 + fs_309375_16 * r_4 * h2_p2) + e_3 * (-fs_2083725_14872 * h8_p2 + fs_33075_208 * h8_p8 + fs_1575_22 * r_2 * h6_p2 + fs_1968300_1859 * r_4 * h4_p2 - fs_28125_44 * r_6 * h2_p2) + e_4 * (-f_6993_8398 * h10_p2 - fs_11907_319124 * h10_p8 + fs_4167450_671099 * r_2 * h8_p2 - fs_33075_4693 * r_2 * h8_p8 - fs_3150_3179 * r_4 * h6_p2 - fs_15552_1859 * r_6 * h4_p2 + fs_28125_7436 * r_8 * h2_p2) + e_5 * (-fs_152460_717470533 * h12_p2 - fs_87120_2221271 * h12_p8 + f_6993_96577 * r_2 * h10_p2 + fs_11907_42204149 * r_2 * h10_p8 - fs_9450_671099 * r_4 * h8_p2 + fs_75_4693 * r_4 * h8_p8 + fs_1400_1147619 * r_6 * h6_p2 + fs_3888_537251 * r_8 * h4_p2 - fs_5_1859 * r_10 * h2_p2);
    }

    // NOTE: the rows are formed in 82 loops, as the vectorizer runs out of
    // registers with all 91 of them in one.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph4_p3, ph6_p3, ph8_p3, ph8_p7, ph10_p3, ph10_p7, ph12_p3, ph12_p7, ab_2, pc_16 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h4_p3 = ph4_p3[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h8_p3 = ph8_p3[k];
        const auto h8_p7 = ph8_p7[k];
        const auto h10_p3 = ph10_p3[k];
        const auto h10_p7 = ph10_p7[k];
        const auto h12_p3 = ph12_p3[k];
        const auto h12_p7 = ph12_p7[k];

        pc_16[k] = -fs_22920975_512 * e_1 * h4_p3 + e_2 * (-fs_39375_88 * h6_p3 + fs_2083725_88 * r_2 * h4_p3) + e_3 * (fs_77175_676 * h8_p3 + fs_33075_208 * h8_p7 + fs_1575_22 * r_2 * h6_p3 - fs_18753525_14872 * r_4 * h4_p3) + e_4 * (fs_3814209_2712554 * h10_p3 - fs_194481_159562 * h10_p7 - fs_308700_61009 * r_2 * h8_p3 - fs_33075_4693 * r_2 * h8_p7 - fs_3150_3179 * r_4 * h6_p3 + fs_18522_1859 * r_6 * h4_p3) + e_5 * (fs_1029105_1434941066 * h12_p3 - fs_49005_2221271 * h12_p7 - fs_7628418_717470533 * r_2 * h10_p3 + fs_388962_42204149 * r_2 * h10_p7 + fs_700_61009 * r_4 * h8_p3 + fs_75_4693 * r_4 * h8_p7 + fs_1400_1147619 * r_6 * h6_p3 - fs_9261_1074502 * r_8 * h4_p3);
    }

    // NOTE: the rows are formed in 82 loops, as the vectorizer runs out of
    // registers with all 91 of them in one.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph4_p4, ph6_m5, ph6_p4, ph6_p6, ph8_p4, ph8_p6, ph10_m5, ph10_p4, ph10_p6, ph12_m5, ph12_p4, ph12_p6, ab_2, pc_17, pc_18 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h4_p4 = ph4_p4[k];
        const auto h6_m5 = ph6_m5[k];
        const auto h6_p4 = ph6_p4[k];
        const auto h6_p6 = ph6_p6[k];
        const auto h8_p4 = ph8_p4[k];
        const auto h8_p6 = ph8_p6[k];
        const auto h10_m5 = ph10_m5[k];
        const auto h10_p4 = ph10_p4[k];
        const auto h10_p6 = ph10_p6[k];
        const auto h12_m5 = ph12_m5[k];
        const auto h12_p4 = ph12_p4[k];
        const auto h12_p6 = ph12_p6[k];

        pc_17[k] = fs_2338875_128 * e_1 * h4_p4 + e_2 * (fs_1063125_352 * h6_p4 + fs_39375_16 * h6_p6 - fs_212625_22 * r_2 * h4_p4) + e_3 * (-fs_231525_5408 * h8_p4 + fs_11025_208 * h8_p6 - fs_42525_88 * r_2 * h6_p4 - fs_1575_4 * r_2 * h6_p6 + fs_1913625_3718 * r_4 * h4_p4) + e_4 * (-fs_24310125_10850216 * h10_p4 - fs_257985_104329 * h10_p6 + fs_231525_122018 * r_2 * h8_p4 - fs_11025_4693 * r_2 * h8_p6 + fs_42525_6358 * r_4 * h6_p4 + fs_1575_289 * r_4 * h6_p6 - fs_7560_1859 * r_6 * h4_p4) + e_5 * (-fs_1463616_717470533 * h12_p4 - fs_470448_42204149 * h12_p6 + fs_24310125_1434941066 * r_2 * h10_p4 + fs_1031940_55190041 * r_2 * h10_p6 - fs_525_122018 * r_4 * h8_p4 + fs_25_4693 * r_4 * h8_p6 - fs_9450_1147619 * r_6 * h6_p4 - fs_700_104329 * r_6 * h6_p6 + fs_1890_537251 * r_8 * h4_p4);

        pc_18[k] = -f_375_4 * e_2 * h6_m5 + f_75_2 * e_3 * r_2 * h6_m5 + e_4 * (fs_83349_15028 * h10_m5 - f_75_17 * r_4 * h6_m5) + e_5 * (fs_426888_42204149 * h12_m5 - fs_83349_1987453 * r_2 * h10_m5 + f_50_323 * r_6 * h6_m5);
    }

    // NOTE: the rows are formed in 82 loops, as the vectorizer runs out of
    // registers with all 91 of them in one.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph4_m4, ph6_m6, ph6_m4, ph8_m6, ph8_m4, ph10_m6, ph10_m4, ph12_m6, ph12_m4, ab_2, pc_19 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h4_m4 = ph4_m4[k];
        const auto h6_m6 = ph6_m6[k];
        const auto h6_m4 = ph6_m4[k];
        const auto h8_m6 = ph8_m6[k];
        const auto h8_m4 = ph8_m4[k];
        const auto h10_m6 = ph10_m6[k];
        const auto h10_m4 = ph10_m4[k];
        const auto h12_m6 = ph12_m6[k];
        const auto h12_m4 = ph12_m4[k];

        pc_19[k] = fs_2338875_128 * e_1 * h4_m4 + e_2 * (-fs_39375_16 * h6_m6 + fs_1063125_352 * h6_m4 - fs_212625_22 * r_2 * h4_m4) + e_3 * (-fs_11025_208 * h8_m6 - fs_231525_5408 * h8_m4 + fs_1575_4 * r_2 * h6_m6 - fs_42525_88 * r_2 * h6_m4 + fs_1913625_3718 * r_4 * h4_m4) + e_4 * (fs_257985_104329 * h10_m6 - fs_24310125_10850216 * h10_m4 + fs_11025_4693 * r_2 * h8_m6 + fs_231525_122018 * r_2 * h8_m4 - fs_1575_289 * r_4 * h6_m6 + fs_42525_6358 * r_4 * h6_m4 - fs_7560_1859 * r_6 * h4_m4) + e_5 * (fs_470448_42204149 * h12_m6 - fs_1463616_717470533 * h12_m4 - fs_1031940_55190041 * r_2 * h10_m6 + fs_24310125_1434941066 * r_2 * h10_m4 - fs_25_4693 * r_4 * h8_m6 - fs_525_122018 * r_4 * h8_m4 + fs_700_104329 * r_6 * h6_m6 - fs_9450_1147619 * r_6 * h6_m4 + fs_1890_537251 * r_8 * h4_m4);
    }

    // NOTE: the rows are formed in 82 loops, as the vectorizer runs out of
    // registers with all 91 of them in one.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph4_m3, ph6_m3, ph8_m7, ph8_m3, ph10_m7, ph10_m3, ph12_m7, ph12_m3, ab_2, pc_20 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h4_m3 = ph4_m3[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h8_m7 = ph8_m7[k];
        const auto h8_m3 = ph8_m3[k];
        const auto h10_m7 = ph10_m7[k];
        const auto h10_m3 = ph10_m3[k];
        const auto h12_m7 = ph12_m7[k];
        const auto h12_m3 = ph12_m3[k];

        pc_20[k] = -fs_22920975_512 * e_1 * h4_m3 + e_2 * (-fs_39375_88 * h6_m3 + fs_2083725_88 * r_2 * h4_m3) + e_3 * (-fs_33075_208 * h8_m7 + fs_77175_676 * h8_m3 + fs_1575_22 * r_2 * h6_m3 - fs_18753525_14872 * r_4 * h4_m3) + e_4 * (fs_194481_159562 * h10_m7 + fs_3814209_2712554 * h10_m3 + fs_33075_4693 * r_2 * h8_m7 - fs_308700_61009 * r_2 * h8_m3 - fs_3150_3179 * r_4 * h6_m3 + fs_18522_1859 * r_6 * h4_m3) + e_5 * (fs_49005_2221271 * h12_m7 + fs_1029105_1434941066 * h12_m3 - fs_388962_42204149 * r_2 * h10_m7 - fs_7628418_717470533 * r_2 * h10_m3 - fs_75_4693 * r_4 * h8_m7 + fs_700_61009 * r_4 * h8_m3 + fs_1400_1147619 * r_6 * h6_m3 - fs_9261_1074502 * r_8 * h4_m3);
    }

    // NOTE: the rows are formed in 82 loops, as the vectorizer runs out of
    // registers with all 91 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph2_m2, ph4_m2, ph6_m2, ph8_m8, ph8_m2, ph10_m8, ph10_m2, ph12_m8, ph12_m2, ab_2, pc_21 : simd::cache_line_size())
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

        const auto h2_m2 = ph2_m2[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h8_m8 = ph8_m8[k];
        const auto h8_m2 = ph8_m2[k];
        const auto h10_m8 = ph10_m8[k];
        const auto h10_m2 = ph10_m2[k];
        const auto h12_m8 = ph12_m8[k];
        const auto h12_m2 = ph12_m2[k];

        pc_21[k] = fs_49116375_1024 * e_0 * h2_m2 + e_1 * (fs_601425_16 * h4_m2 - fs_25059375_256 * r_2 * h2_m2) + e_2 * (-fs_39375_88 * h6_m2 - fs_218700_11 * r_2 * h4_m2 + fs_309375_16 * r_4 * h2_m2) + e_3 * (-fs_33075_208 * h8_m8 - fs_2083725_14872 * h8_m2 + fs_1575_22 * r_2 * h6_m2 + fs_1968300_1859 * r_4 * h4_m2 - fs_28125_44 * r_6 * h2_m2) + e_4 * (fs_11907_319124 * h10_m8 - f_6993_8398 * h10_m2 + fs_33075_4693 * r_2 * h8_m8 + fs_4167450_671099 * r_2 * h8_m2 - fs_3150_3179 * r_4 * h6_m2 - fs_15552_1859 * r_6 * h4_m2 + fs_28125_7436 * r_8 * h2_m2) + e_5 * (fs_87120_2221271 * h12_m8 - fs_152460_717470533 * h12_m2 - fs_11907_42204149 * r_2 * h10_m8 + f_6993_96577 * r_2 * h10_m2 - fs_75_4693 * r_4 * h8_m8 - fs_9450_671099 * r_4 * h8_m2 + fs_1400_1147619 * r_6 * h6_m2 + fs_3888_537251 * r_8 * h4_m2 - fs_5_1859 * r_10 * h2_m2);
    }

    // NOTE: the rows are formed in 82 loops, as the vectorizer runs out of
    // registers with all 91 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph2_m1, ph4_m1, ph6_m1, ph8_m1, ph10_m10, ph10_m9, ph10_m1, ph12_m11, ph12_m10, ph12_m9, ph12_m1, ab_2, pc_22, pc_23, pc_24 : simd::cache_line_size())
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

        const auto h2_m1 = ph2_m1[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h8_m1 = ph8_m1[k];
        const auto h10_m10 = ph10_m10[k];
        const auto h10_m9 = ph10_m9[k];
        const auto h10_m1 = ph10_m1[k];
        const auto h12_m11 = ph12_m11[k];
        const auto h12_m10 = ph12_m10[k];
        const auto h12_m9 = ph12_m9[k];
        const auto h12_m1 = ph12_m1[k];

        pc_22[k] = -fs_265228425_2048 * e_0 * h2_m1 + e_1 * (-fs_1002375_256 * h4_m1 + fs_135320625_512 * r_2 * h2_m1) + e_2 * (fs_1063125_352 * h6_m1 + fs_91125_44 * r_2 * h4_m1 - fs_1670625_32 * r_4 * h2_m1) + e_3 * (fs_198450_1859 * h8_m1 - fs_42525_88 * r_2 * h6_m1 - fs_820125_7436 * r_4 * h4_m1 + fs_151875_88 * r_6 * h2_m1) + e_4 * (-fs_19845_16796 * h10_m9 + fs_36693405_141052808 * h10_m1 - fs_3175200_671099 * r_2 * h8_m1 + fs_42525_6358 * r_4 * h6_m1 + fs_1620_1859 * r_6 * h4_m1 - fs_151875_14872 * r_8 * h2_m1) + e_5 * (fs_137214_2221271 * h12_m9 + fs_35937_717470533 * h12_m1 + fs_19845_2221271 * r_2 * h10_m9 - fs_36693405_18654233858 * r_2 * h10_m1 + fs_7200_671099 * r_4 * h8_m1 - fs_9450_1147619 * r_6 * h6_m1 - fs_405_537251 * r_8 * h4_m1 + fs_27_3718 * r_10 * h2_m1);

        pc_23[k] = -fs_43659_8398 * e_4 * h10_m10 + e_5 * (fs_182952_2221271 * h12_m10 + fs_87318_2221271 * r_2 * h10_m10);

        pc_24[k] = f_10395_32 * e_0 * h2_m1 + e_1 * (fs_33078375_512 * h4_m1 - f_7425_16 * r_2 * h2_m1) + e_2 * (fs_39375_16 * h6_m1 - fs_273375_8 * r_2 * h4_m1 + f_825_4 * r_4 * h2_m1) + e_3 * (fs_33075_2704 * h8_m1 - fs_1575_4 * r_2 * h6_m1 + fs_2460375_1352 * r_4 * h4_m1 - f_75_2 * r_6 * h2_m1) + e_4 * (fs_654885_70526404 * h10_m1 - fs_33075_61009 * r_2 * h8_m1 + fs_1575_289 * r_4 * h6_m1 - fs_2430_169 * r_6 * h4_m1 + f_75_26 * r_8 * h2_m1) + e_5 * (fs_7623_96577 * h12_m11 + fs_1089_1434941066 * h12_m1 - fs_654885_9327116929 * r_2 * h10_m1 + fs_75_61009 * r_4 * h8_m1 - fs_700_104329 * r_6 * h6_m1 + fs_1215_97682 * r_8 * h4_m1 - f_1_13 * r_10 * h2_m1);
    }

    // NOTE: the rows are formed in 82 loops, as the vectorizer runs out of
    // registers with all 91 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph2_0, ph4_0, ph6_0, ph8_0, ph8_p8, ph10_0, ph10_p8, ph12_0, ph12_p8, ab_2, pc_25 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;
        const auto r_10 = r_8 * r_2;
        const auto r_12 = r_10 * r_2;

        const auto h2_0 = ph2_0[k];
        const auto h4_0 = ph4_0[k];
        const auto h6_0 = ph6_0[k];
        const auto h8_0 = ph8_0[k];
        const auto h8_p8 = ph8_p8[k];
        const auto h10_0 = ph10_0[k];
        const auto h10_p8 = ph10_p8[k];
        const auto h12_0 = ph12_0[k];
        const auto h12_p8 = ph12_p8[k];

        pc_25[k] = e_0 * (f_945_16 * h2_0 - f_10395_16 * r_2) + e_1 * (-f_270 * h4_0 - f_675_8 * r_2 * h2_0 + f_10395_16 * r_4) + e_2 * (f_150_11 * h6_0 + f_2160_11 * r_2 * h4_0 + f_75_2 * r_4 * h2_0 - f_495_2 * r_6) + e_3 * (f_9345_572 * h8_0 - fs_496125_2288 * h8_p8 - f_60_11 * r_2 * h6_0 - f_6480_143 * r_4 * h4_0 - f_75_11 * r_6 * h2_0 + f_165_4 * r_8) + e_4 * (f_5229_4199 * h10_0 + fs_218295_79781 * h10_p8 - f_9345_2717 * r_2 * h8_0 + fs_496125_51623 * r_2 * h8_p8 + f_120_187 * r_4 * h6_0 + f_576_143 * r_6 * h4_0 + f_75_143 * r_8 * h2_0 - f_3 * r_10) + e_5 * (f_2178_96577 * h12_0 - fs_143748_2221271 * h12_p8 - f_10458_96577 * r_2 * h10_0 - fs_873180_42204149 * r_2 * h10_p8 + f_445_2717 * r_4 * h8_0 - fs_1125_51623 * r_4 * h8_p8 - f_80_3553 * r_6 * h6_0 - f_288_2431 * r_8 * h4_0 - f_2_143 * r_10 * h2_0 + f_1_13 * r_12) + f_10395_64 * e_6;
    }

    // NOTE: the rows are formed in 82 loops, as the vectorizer runs out of
    // registers with all 91 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph2_p1, ph4_p1, ph6_p1, ph8_p1, ph8_p7, ph10_p1, ph10_p7, ph12_p1, ph12_p7, ab_2, pc_26 : simd::cache_line_size())
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

        const auto h2_p1 = ph2_p1[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h8_p7 = ph8_p7[k];
        const auto h10_p1 = ph10_p1[k];
        const auto h10_p7 = ph10_p7[k];
        const auto h12_p1 = ph12_p1[k];
        const auto h12_p7 = ph12_p7[k];

        pc_26[k] = -fs_218791125_2048 * e_0 * h2_p1 + e_1 * (fs_2679075_256 * h4_p1 + fs_111628125_512 * r_2 * h2_p1) + e_2 * (fs_4921875_3872 * h6_p1 - fs_2679075_484 * r_2 * h4_p1 - fs_1378125_32 * r_4 * h2_p1) + e_3 * (-fs_13395375_163592 * h8_p1 - fs_33075_1144 * h8_p7 - fs_196875_968 * r_2 * h6_p1 + fs_24111675_81796 * r_4 * h4_p1 + fs_1378125_968 * r_6 * h2_p1) + e_4 * (-fs_220172337_141052808 * h10_p1 + fs_392931_319124 * h10_p7 + fs_26790750_7382089 * r_2 * h8_p1 + fs_66150_51623 * r_2 * h8_p7 + fs_196875_69938 * r_4 * h6_p1 - fs_47628_20449 * r_6 * h4_p1 - fs_1378125_163592 * r_8 * h2_p1) + e_5 * (-fs_658845_717470533 * h12_p1 - fs_119790_2221271 * h12_p7 + fs_220172337_18654233858 * r_2 * h10_p1 - fs_392931_42204149 * r_2 * h10_p7 - fs_60750_7382089 * r_4 * h8_p1 - fs_150_51623 * r_4 * h8_p7 - fs_43750_12623809 * r_6 * h6_p1 + fs_11907_5909761 * r_8 * h4_p1 + fs_245_40898 * r_10 * h2_p1);
    }

    // NOTE: the rows are formed in 82 loops, as the vectorizer runs out of
    // registers with all 91 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph2_p2, ph4_p2, ph6_p2, ph6_p6, ph8_p2, ph8_p6, ph10_p2, ph10_p6, ph12_p2, ph12_p6, ab_2, pc_27 : simd::cache_line_size())
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

        const auto h2_p2 = ph2_p2[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h6_p6 = ph6_p6[k];
        const auto h8_p2 = ph8_p2[k];
        const auto h8_p6 = ph8_p6[k];
        const auto h10_p2 = ph10_p2[k];
        const auto h10_p6 = ph10_p6[k];
        const auto h12_p2 = ph12_p2[k];
        const auto h12_p6 = ph12_p6[k];

        pc_27[k] = fs_40186125_512 * e_0 * h2_p2 + e_1 * (fs_3213675_512 * h4_p2 - fs_20503125_128 * r_2 * h2_p2) + e_2 * (-fs_354375_121 * h6_p2 - fs_196875_44 * h6_p6 - fs_3213675_968 * r_2 * h4_p2 + fs_253125_8 * r_4 * h2_p2) + e_3 * (fs_231525_20449 * h8_p2 + fs_55125_2288 * h8_p6 + fs_56700_121 * r_2 * h6_p2 + fs_7875_11 * r_2 * h6_p6 + fs_28923075_163592 * r_4 * h4_p2 - fs_253125_242 * r_6 * h2_p2) + e_4 * (fs_80725491_35263202 * h10_p2 + fs_43659_1356277 * h10_p6 - fs_3704400_7382089 * r_2 * h8_p2 - fs_55125_51623 * r_2 * h8_p6 - fs_226800_34969 * r_4 * h6_p2 - fs_31500_3179 * r_4 * h6_p6 - fs_28566_20449 * r_6 * h4_p2 + fs_253125_40898 * r_8 * h2_p2) + e_5 * (fs_3773385_1434941066 * h12_p2 - fs_1617165_42204149 * h12_p6 - fs_161450982_9327116929 * r_2 * h10_p2 - fs_174636_717470533 * r_2 * h10_p6 + fs_8400_7382089 * r_4 * h8_p2 + fs_125_51623 * r_4 * h8_p6 + fs_100800_12623809 * r_6 * h6_p2 + fs_14000_1147619 * r_6 * h6_p6 + fs_14283_11819522 * r_8 * h4_p2 - fs_90_20449 * r_10 * h2_p2);
    }

    // NOTE: the rows are formed in 82 loops, as the vectorizer runs out of
    // registers with all 91 of them in one.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph4_p3, ph6_p3, ph6_p5, ph8_p3, ph8_p5, ph10_p3, ph10_p5, ph12_p3, ph12_p5, ab_2, pc_28 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h4_p3 = ph4_p3[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h6_p5 = ph6_p5[k];
        const auto h8_p3 = ph8_p3[k];
        const auto h8_p5 = ph8_p5[k];
        const auto h10_p3 = ph10_p3[k];
        const auto h10_p5 = ph10_p5[k];
        const auto h12_p3 = ph12_p3[k];
        const auto h12_p5 = ph12_p5[k];

        pc_28[k] = -fs_5315625_128 * e_1 * h4_p3 + e_2 * (fs_4921875_3872 * h6_p3 - fs_1063125_352 * h6_p5 + fs_5315625_242 * r_2 * h4_p3) + e_3 * (fs_385875_29744 * h8_p3 + fs_231525_2288 * h8_p5 - fs_196875_968 * r_2 * h6_p3 + fs_42525_88 * r_2 * h6_p5 - fs_47840625_40898 * r_4 * h4_p3) + e_4 * (-fs_26413695_10850216 * h10_p3 - fs_5282739_10850216 * h10_p5 - fs_385875_671099 * r_2 * h8_p3 - fs_231525_51623 * r_2 * h8_p5 + fs_196875_69938 * r_4 * h6_p3 - fs_42525_6358 * r_4 * h6_p5 + fs_189000_20449 * r_6 * h4_p3) + e_5 * (-fs_4528062_717470533 * h12_p3 - fs_1006236_42204149 * h12_p5 + fs_26413695_1434941066 * r_2 * h10_p3 + fs_5282739_1434941066 * r_2 * h10_p5 + fs_875_671099 * r_4 * h8_p3 + fs_525_51623 * r_4 * h8_p5 - fs_43750_12623809 * r_6 * h6_p3 + fs_9450_1147619 * r_6 * h6_p5 - fs_47250_5909761 * r_8 * h4_p3);
    }

    // NOTE: the rows are formed in 82 loops, as the vectorizer runs out of
    // registers with all 91 of them in one.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph4_m4, ph6_m4, ph8_m4, ph10_m4, ph12_m4, ab_2, pc_29 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h4_m4 = ph4_m4[k];
        const auto h6_m4 = ph6_m4[k];
        const auto h8_m4 = ph8_m4[k];
        const auto h10_m4 = ph10_m4[k];
        const auto h12_m4 = ph12_m4[k];

        pc_29[k] = fs_4465125_64 * e_1 * h4_m4 + e_2 * (f_150_11 * h6_m4 - fs_4465125_121 * r_2 * h4_m4) + e_3 * (-fs_4862025_29744 * h8_m4 - f_60_11 * r_2 * h6_m4 + fs_40186125_20449 * r_4 * h4_m4) + e_4 * (fs_4584195_1356277 * h10_m4 + fs_4862025_671099 * r_2 * h8_m4 + f_120_187 * r_4 * h6_m4 - fs_317520_20449 * r_6 * h4_m4) + e_5 * (fs_18783072_717470533 * h12_m4 - fs_18336780_717470533 * r_2 * h10_m4 - fs_11025_671099 * r_4 * h8_m4 - f_80_3553 * r_6 * h6_m4 + fs_79380_5909761 * r_8 * h4_m4);
    }

    // NOTE: the rows are formed in 82 loops, as the vectorizer runs out of
    // registers with all 91 of them in one.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph4_m3, ph6_m5, ph6_m3, ph8_m5, ph8_m3, ph10_m5, ph10_m3, ph12_m5, ph12_m3, ab_2, pc_30 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h4_m3 = ph4_m3[k];
        const auto h6_m5 = ph6_m5[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h8_m5 = ph8_m5[k];
        const auto h8_m3 = ph8_m3[k];
        const auto h10_m5 = ph10_m5[k];
        const auto h10_m3 = ph10_m3[k];
        const auto h12_m5 = ph12_m5[k];
        const auto h12_m3 = ph12_m3[k];

        pc_30[k] = -fs_5315625_128 * e_1 * h4_m3 + e_2 * (fs_1063125_352 * h6_m5 + fs_4921875_3872 * h6_m3 + fs_5315625_242 * r_2 * h4_m3) + e_3 * (-fs_231525_2288 * h8_m5 + fs_385875_29744 * h8_m3 - fs_42525_88 * r_2 * h6_m5 - fs_196875_968 * r_2 * h6_m3 - fs_47840625_40898 * r_4 * h4_m3) + e_4 * (fs_5282739_10850216 * h10_m5 - fs_26413695_10850216 * h10_m3 + fs_231525_51623 * r_2 * h8_m5 - fs_385875_671099 * r_2 * h8_m3 + fs_42525_6358 * r_4 * h6_m5 + fs_196875_69938 * r_4 * h6_m3 + fs_189000_20449 * r_6 * h4_m3) + e_5 * (fs_1006236_42204149 * h12_m5 - fs_4528062_717470533 * h12_m3 - fs_5282739_1434941066 * r_2 * h10_m5 + fs_26413695_1434941066 * r_2 * h10_m3 - fs_525_51623 * r_4 * h8_m5 + fs_875_671099 * r_4 * h8_m3 - fs_9450_1147619 * r_6 * h6_m5 - fs_43750_12623809 * r_6 * h6_m3 - fs_47250_5909761 * r_8 * h4_m3);
    }

    // NOTE: the rows are formed in 82 loops, as the vectorizer runs out of
    // registers with all 91 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph2_m2, ph4_m2, ph6_m6, ph6_m2, ph8_m6, ph8_m2, ph10_m6, ph10_m2, ph12_m6, ph12_m2, ab_2, pc_31 : simd::cache_line_size())
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

        const auto h2_m2 = ph2_m2[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h6_m6 = ph6_m6[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h8_m6 = ph8_m6[k];
        const auto h8_m2 = ph8_m2[k];
        const auto h10_m6 = ph10_m6[k];
        const auto h10_m2 = ph10_m2[k];
        const auto h12_m6 = ph12_m6[k];
        const auto h12_m2 = ph12_m2[k];

        pc_31[k] = fs_40186125_512 * e_0 * h2_m2 + e_1 * (fs_3213675_512 * h4_m2 - fs_20503125_128 * r_2 * h2_m2) + e_2 * (fs_196875_44 * h6_m6 - fs_354375_121 * h6_m2 - fs_3213675_968 * r_2 * h4_m2 + fs_253125_8 * r_4 * h2_m2) + e_3 * (-fs_55125_2288 * h8_m6 + fs_231525_20449 * h8_m2 - fs_7875_11 * r_2 * h6_m6 + fs_56700_121 * r_2 * h6_m2 + fs_28923075_163592 * r_4 * h4_m2 - fs_253125_242 * r_6 * h2_m2) + e_4 * (-fs_43659_1356277 * h10_m6 + fs_80725491_35263202 * h10_m2 + fs_55125_51623 * r_2 * h8_m6 - fs_3704400_7382089 * r_2 * h8_m2 + fs_31500_3179 * r_4 * h6_m6 - fs_226800_34969 * r_4 * h6_m2 - fs_28566_20449 * r_6 * h4_m2 + fs_253125_40898 * r_8 * h2_m2) + e_5 * (fs_1617165_42204149 * h12_m6 + fs_3773385_1434941066 * h12_m2 + fs_174636_717470533 * r_2 * h10_m6 - fs_161450982_9327116929 * r_2 * h10_m2 - fs_125_51623 * r_4 * h8_m6 + fs_8400_7382089 * r_4 * h8_m2 - fs_14000_1147619 * r_6 * h6_m6 + fs_100800_12623809 * r_6 * h6_m2 + fs_14283_11819522 * r_8 * h4_m2 - fs_90_20449 * r_10 * h2_m2);
    }

    // NOTE: the rows are formed in 82 loops, as the vectorizer runs out of
    // registers with all 91 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph2_m1, ph4_m1, ph6_m1, ph8_m7, ph8_m1, ph10_m7, ph10_m1, ph12_m7, ph12_m1, ab_2, pc_32 : simd::cache_line_size())
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

        const auto h2_m1 = ph2_m1[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h8_m7 = ph8_m7[k];
        const auto h8_m1 = ph8_m1[k];
        const auto h10_m7 = ph10_m7[k];
        const auto h10_m1 = ph10_m1[k];
        const auto h12_m7 = ph12_m7[k];
        const auto h12_m1 = ph12_m1[k];

        pc_32[k] = -fs_218791125_2048 * e_0 * h2_m1 + e_1 * (fs_2679075_256 * h4_m1 + fs_111628125_512 * r_2 * h2_m1) + e_2 * (fs_4921875_3872 * h6_m1 - fs_2679075_484 * r_2 * h4_m1 - fs_1378125_32 * r_4 * h2_m1) + e_3 * (fs_33075_1144 * h8_m7 - fs_13395375_163592 * h8_m1 - fs_196875_968 * r_2 * h6_m1 + fs_24111675_81796 * r_4 * h4_m1 + fs_1378125_968 * r_6 * h2_m1) + e_4 * (-fs_392931_319124 * h10_m7 - fs_220172337_141052808 * h10_m1 - fs_66150_51623 * r_2 * h8_m7 + fs_26790750_7382089 * r_2 * h8_m1 + fs_196875_69938 * r_4 * h6_m1 - fs_47628_20449 * r_6 * h4_m1 - fs_1378125_163592 * r_8 * h2_m1) + e_5 * (fs_119790_2221271 * h12_m7 - fs_658845_717470533 * h12_m1 + fs_392931_42204149 * r_2 * h10_m7 + fs_220172337_18654233858 * r_2 * h10_m1 + fs_150_51623 * r_4 * h8_m7 - fs_60750_7382089 * r_4 * h8_m1 - fs_43750_12623809 * r_6 * h6_m1 + fs_11907_5909761 * r_8 * h4_m1 + fs_245_40898 * r_10 * h2_m1);
    }

    // NOTE: the rows are formed in 82 loops, as the vectorizer runs out of
    // registers with all 91 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph2_m1, ph4_m1, ph6_m1, ph8_m8, ph8_m1, ph10_m9, ph10_m8, ph10_m1, ph12_m9, ph12_m8, ph12_m1, ab_2, pc_33, pc_34 : simd::cache_line_size())
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

        const auto h2_m1 = ph2_m1[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h8_m8 = ph8_m8[k];
        const auto h8_m1 = ph8_m1[k];
        const auto h10_m9 = ph10_m9[k];
        const auto h10_m8 = ph10_m8[k];
        const auto h10_m1 = ph10_m1[k];
        const auto h12_m9 = ph12_m9[k];
        const auto h12_m8 = ph12_m8[k];
        const auto h12_m1 = ph12_m1[k];

        pc_33[k] = fs_496125_2288 * e_3 * h8_m8 + e_4 * (-fs_218295_79781 * h10_m8 - fs_496125_51623 * r_2 * h8_m8) + e_5 * (fs_143748_2221271 * h12_m8 + fs_873180_42204149 * r_2 * h10_m8 + fs_1125_51623 * r_4 * h8_m8);

        pc_34[k] = fs_265228425_2048 * e_0 * h2_m1 + e_1 * (fs_1002375_256 * h4_m1 - fs_135320625_512 * r_2 * h2_m1) + e_2 * (-fs_1063125_352 * h6_m1 - fs_91125_44 * r_2 * h4_m1 + fs_1670625_32 * r_4 * h2_m1) + e_3 * (-fs_198450_1859 * h8_m1 + fs_42525_88 * r_2 * h6_m1 + fs_820125_7436 * r_4 * h4_m1 - fs_151875_88 * r_6 * h2_m1) + e_4 * (-fs_19845_16796 * h10_m9 - fs_36693405_141052808 * h10_m1 + fs_3175200_671099 * r_2 * h8_m1 - fs_42525_6358 * r_4 * h6_m1 - fs_1620_1859 * r_6 * h4_m1 + fs_151875_14872 * r_8 * h2_m1) + e_5 * (fs_137214_2221271 * h12_m9 - fs_35937_717470533 * h12_m1 + fs_19845_2221271 * r_2 * h10_m9 + fs_36693405_18654233858 * r_2 * h10_m1 - fs_7200_671099 * r_4 * h8_m1 + fs_9450_1147619 * r_6 * h6_m1 + fs_405_537251 * r_8 * h4_m1 - fs_27_3718 * r_10 * h2_m1);
    }

    // NOTE: the rows are formed in 82 loops, as the vectorizer runs out of
    // registers with all 91 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph2_m2, ph4_m2, ph6_m2, ph8_m2, ph10_m10, ph10_m2, ph12_m10, ph12_m2, ab_2, pc_35 : simd::cache_line_size())
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

        const auto h2_m2 = ph2_m2[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h8_m2 = ph8_m2[k];
        const auto h10_m10 = ph10_m10[k];
        const auto h10_m2 = ph10_m2[k];
        const auto h12_m10 = ph12_m10[k];
        const auto h12_m2 = ph12_m2[k];

        pc_35[k] = -fs_9823275_512 * e_0 * h2_m2 + e_1 * (-fs_27064125_512 * h4_m2 + fs_5011875_128 * r_2 * h2_m2) + e_2 * (-fs_196875_44 * h6_m2 + fs_2460375_88 * r_2 * h4_m2 - fs_61875_8 * r_4 * h2_m2) + e_3 * (-fs_1157625_29744 * h8_m2 + fs_7875_11 * r_2 * h6_m2 - fs_22143375_14872 * r_4 * h4_m2 + fs_5625_22 * r_6 * h2_m2) + e_4 * (fs_11907_4199 * h10_m10 - fs_1607445_35263202 * h10_m2 + fs_1157625_671099 * r_2 * h8_m2 - fs_31500_3179 * r_4 * h6_m2 + fs_21870_1859 * r_6 * h4_m2 - fs_5625_3718 * r_8 * h2_m2) + e_5 * (fs_83853_2221271 * h12_m10 - fs_7623_1434941066 * h12_m2 - fs_47628_2221271 * r_2 * h10_m10 + fs_3214890_9327116929 * r_2 * h10_m2 - fs_2625_671099 * r_4 * h8_m2 + fs_14000_1147619 * r_6 * h6_m2 - fs_10935_1074502 * r_8 * h4_m2 + fs_2_1859 * r_10 * h2_m2);
    }

    // NOTE: the rows are formed in 82 loops, as the vectorizer runs out of
    // registers with all 91 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph2_0, ph4_0, ph6_0, ph6_p6, ph8_0, ph8_p6, ph10_0, ph10_p6, ph12_0, ph12_p6, ab_2, pc_36 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;
        const auto r_10 = r_8 * r_2;
        const auto r_12 = r_10 * r_2;

        const auto h2_0 = ph2_0[k];
        const auto h4_0 = ph4_0[k];
        const auto h6_0 = ph6_0[k];
        const auto h6_p6 = ph6_p6[k];
        const auto h8_0 = ph8_0[k];
        const auto h8_p6 = ph8_p6[k];
        const auto h10_0 = ph10_0[k];
        const auto h10_p6 = ph10_p6[k];
        const auto h12_0 = ph12_0[k];
        const auto h12_p6 = ph12_p6[k];

        pc_36[k] = e_0 * (-f_4725_32 * h2_0 - f_10395_16 * r_2) + e_1 * (-f_1215_8 * h4_0 + f_3375_16 * r_2 * h2_0 + f_10395_16 * r_4) + e_2 * (f_3225_44 * h6_0 + fs_118125_22 * h6_p6 + f_1215_11 * r_2 * h4_0 - f_375_4 * r_4 * h2_0 - f_495_2 * r_6) + e_3 * (-f_1995_572 * h8_0 - fs_33075_286 * h8_p6 - f_645_22 * r_2 * h6_0 - fs_9450_11 * r_2 * h6_p6 - f_3645_143 * r_4 * h4_0 + f_375_22 * r_6 * h2_0 + f_165_4 * r_8) + e_4 * (-f_945_442 * h10_0 + fs_5893965_2712554 * h10_p6 + f_105_143 * r_2 * h8_0 + fs_264600_51623 * r_2 * h8_p6 + f_645_187 * r_4 * h6_0 + fs_37800_3179 * r_4 * h6_p6 + f_324_143 * r_6 * h4_0 - f_375_286 * r_8 * h2_0 - f_3 * r_10) + e_5 * (-f_7260_96577 * h12_0 - fs_2395800_42204149 * h12_p6 + f_945_5083 * r_2 * h10_0 - fs_11787930_717470533 * r_2 * h10_p6 - f_5_143 * r_4 * h8_0 - fs_600_51623 * r_4 * h8_p6 - f_430_3553 * r_6 * h6_0 - fs_16800_1147619 * r_6 * h6_p6 - f_162_2431 * r_8 * h4_0 + f_5_143 * r_10 * h2_0 + f_1_13 * r_12) + f_10395_64 * e_6;
    }

    // NOTE: the rows are formed in 82 loops, as the vectorizer runs out of
    // registers with all 91 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph2_p1, ph4_p1, ph6_p1, ph6_p5, ph8_p1, ph8_p5, ph10_p1, ph10_p5, ph12_p1, ph12_p5, ab_2, pc_37 : simd::cache_line_size())
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

        const auto h2_p1 = ph2_p1[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h6_p5 = ph6_p5[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h8_p5 = ph8_p5[k];
        const auto h10_p1 = ph10_p1[k];
        const auto h10_p5 = ph10_p5[k];
        const auto h12_p1 = ph12_p1[k];
        const auto h12_p5 = ph12_p5[k];

        pc_37[k] = -fs_66976875_1024 * e_0 * h2_p1 + e_1 * (fs_15400125_512 * h4_p1 + fs_34171875_256 * r_2 * h2_p1) + e_2 * (-fs_1063125_1936 * h6_p1 + fs_39375_88 * h6_p5 - fs_15400125_968 * r_2 * h4_p1 - fs_421875_16 * r_4 * h2_p1) + e_3 * (-f_105_22 * h8_p1 - fs_77175_2288 * h8_p5 + fs_42525_484 * r_2 * h6_p1 - fs_1575_22 * r_2 * h6_p5 + fs_820125_968 * r_4 * h4_p1 + fs_421875_484 * r_6 * h2_p1) + e_4 * (fs_159137055_70526404 * h10_p1 + fs_3274425_2712554 * h10_p5 + f_210_209 * r_2 * h8_p1 + fs_77175_51623 * r_2 * h8_p5 - fs_42525_34969 * r_4 * h6_p1 + fs_3150_3179 * r_4 * h6_p5 - fs_810_121 * r_6 * h4_p1 - fs_421875_81796 * r_8 * h2_p1) + e_5 * (fs_9882675_1434941066 * h12_p1 - fs_2096325_42204149 * h12_p5 - fs_159137055_9327116929 * r_2 * h10_p1 - fs_6548850_717470533 * r_2 * h10_p5 - f_10_209 * r_4 * h8_p1 - fs_175_51623 * r_4 * h8_p5 + fs_18900_12623809 * r_6 * h6_p1 - fs_1400_1147619 * r_6 * h6_p5 + fs_405_69938 * r_8 * h4_p1 + fs_75_20449 * r_10 * h2_p1);
    }

    // NOTE: the rows are formed in 82 loops, as the vectorizer runs out of
    // registers with all 91 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph2_p2, ph4_p2, ph4_p4, ph6_p2, ph6_p4, ph8_p2, ph8_p4, ph10_p2, ph10_p4, ph12_p2, ph12_p4, ab_2, pc_38 : simd::cache_line_size())
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

        const auto h2_p2 = ph2_p2[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h4_p4 = ph4_p4[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h6_p4 = ph6_p4[k];
        const auto h8_p2 = ph8_p2[k];
        const auto h8_p4 = ph8_p4[k];
        const auto h10_p2 = ph10_p2[k];
        const auto h10_p4 = ph10_p4[k];
        const auto h12_p2 = ph12_p2[k];
        const auto h12_p4 = ph12_p4[k];

        pc_38[k] = fs_13395375_128 * e_0 * h2_p2 + e_1 * (-fs_18225_8 * h4_p2 + fs_6251175_128 * h4_p4 - fs_6834375_32 * r_2 * h2_p2) + e_2 * (-fs_1063125_1936 * h6_p2 - fs_4921875_3872 * h6_p4 + fs_145800_121 * r_2 * h4_p2 - fs_6251175_242 * r_2 * h4_p4 + fs_84375_2 * r_4 * h2_p2) + e_3 * (fs_27860175_327184 * h8_p2 + fs_385875_59488 * h8_p4 + fs_42525_484 * r_2 * h6_p2 + fs_196875_968 * r_2 * h6_p4 - fs_1312200_20449 * r_4 * h4_p2 + fs_56260575_40898 * r_4 * h4_p4 - fs_168750_121 * r_6 * h2_p2) + e_4 * (-fs_130977_97682 * h10_p2 + fs_1178793_10850216 * h10_p4 - fs_77175_20449 * r_2 * h8_p2 - fs_385875_1342198 * r_2 * h8_p4 - fs_42525_34969 * r_4 * h6_p2 - fs_196875_69938 * r_4 * h6_p4 + fs_10368_20449 * r_6 * h4_p2 - fs_222264_20449 * r_6 * h4_p4 + fs_168750_20449 * r_8 * h2_p2) + e_5 * (-fs_10062360_717470533 * h12_p2 - fs_26832960_717470533 * h12_p4 + fs_261954_25836889 * r_2 * h10_p2 - fs_1178793_1434941066 * r_2 * h10_p4 + fs_175_20449 * r_4 * h8_p2 + fs_875_1342198 * r_4 * h8_p4 + fs_18900_12623809 * r_6 * h6_p2 + fs_43750_12623809 * r_6 * h6_p4 - fs_2592_5909761 * r_8 * h4_p2 + fs_55566_5909761 * r_8 * h4_p4 - fs_120_20449 * r_10 * h2_p2);
    }

    // NOTE: the rows are formed in 82 loops, as the vectorizer runs out of
    // registers with all 91 of them in one.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph4_m3, ph6_m3, ph8_m3, ph10_m3, ph12_m3, ab_2, pc_39 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h4_m3 = ph4_m3[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h8_m3 = ph8_m3[k];
        const auto h10_m3 = ph10_m3[k];
        const auto h12_m3 = ph12_m3[k];

        pc_39[k] = -fs_2679075_64 * e_1 * h4_m3 + e_2 * (f_3225_44 * h6_m3 + fs_2679075_121 * r_2 * h4_m3) + e_3 * (-fs_540225_3718 * h8_m3 - f_645_22 * r_2 * h6_m3 - fs_24111675_20449 * r_4 * h4_m3) + e_4 * (fs_2750517_5425108 * h10_m3 + fs_4321800_671099 * r_2 * h8_m3 + f_645_187 * r_4 * h6_m3 + fs_190512_20449 * r_6 * h4_m3) + e_5 * (fs_35218260_717470533 * h12_m3 - fs_2750517_717470533 * r_2 * h10_m3 - fs_9800_671099 * r_4 * h8_m3 - f_430_3553 * r_6 * h6_m3 - fs_47628_5909761 * r_8 * h4_m3);
    }

    // NOTE: the rows are formed in 82 loops, as the vectorizer runs out of
    // registers with all 91 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph2_m2, ph4_m4, ph4_m2, ph6_m4, ph6_m2, ph8_m4, ph8_m2, ph10_m4, ph10_m2, ph12_m4, ph12_m2, ab_2, pc_40 : simd::cache_line_size())
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

        const auto h2_m2 = ph2_m2[k];
        const auto h4_m4 = ph4_m4[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h6_m4 = ph6_m4[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h8_m4 = ph8_m4[k];
        const auto h8_m2 = ph8_m2[k];
        const auto h10_m4 = ph10_m4[k];
        const auto h10_m2 = ph10_m2[k];
        const auto h12_m4 = ph12_m4[k];
        const auto h12_m2 = ph12_m2[k];

        pc_40[k] = fs_13395375_128 * e_0 * h2_m2 + e_1 * (-fs_6251175_128 * h4_m4 - fs_18225_8 * h4_m2 - fs_6834375_32 * r_2 * h2_m2) + e_2 * (fs_4921875_3872 * h6_m4 - fs_1063125_1936 * h6_m2 + fs_6251175_242 * r_2 * h4_m4 + fs_145800_121 * r_2 * h4_m2 + fs_84375_2 * r_4 * h2_m2) + e_3 * (-fs_385875_59488 * h8_m4 + fs_27860175_327184 * h8_m2 - fs_196875_968 * r_2 * h6_m4 + fs_42525_484 * r_2 * h6_m2 - fs_56260575_40898 * r_4 * h4_m4 - fs_1312200_20449 * r_4 * h4_m2 - fs_168750_121 * r_6 * h2_m2) + e_4 * (-fs_1178793_10850216 * h10_m4 - fs_130977_97682 * h10_m2 + fs_385875_1342198 * r_2 * h8_m4 - fs_77175_20449 * r_2 * h8_m2 + fs_196875_69938 * r_4 * h6_m4 - fs_42525_34969 * r_4 * h6_m2 + fs_222264_20449 * r_6 * h4_m4 + fs_10368_20449 * r_6 * h4_m2 + fs_168750_20449 * r_8 * h2_m2) + e_5 * (fs_26832960_717470533 * h12_m4 - fs_10062360_717470533 * h12_m2 + fs_1178793_1434941066 * r_2 * h10_m4 + fs_261954_25836889 * r_2 * h10_m2 - fs_875_1342198 * r_4 * h8_m4 + fs_175_20449 * r_4 * h8_m2 - fs_43750_12623809 * r_6 * h6_m4 + fs_18900_12623809 * r_6 * h6_m2 - fs_55566_5909761 * r_8 * h4_m4 - fs_2592_5909761 * r_8 * h4_m2 - fs_120_20449 * r_10 * h2_m2);
    }

    // NOTE: the rows are formed in 82 loops, as the vectorizer runs out of
    // registers with all 91 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph2_m1, ph4_m1, ph6_m5, ph6_m1, ph8_m5, ph8_m1, ph10_m5, ph10_m1, ph12_m5, ph12_m1, ab_2, pc_41 : simd::cache_line_size())
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

        const auto h2_m1 = ph2_m1[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h6_m5 = ph6_m5[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h8_m5 = ph8_m5[k];
        const auto h8_m1 = ph8_m1[k];
        const auto h10_m5 = ph10_m5[k];
        const auto h10_m1 = ph10_m1[k];
        const auto h12_m5 = ph12_m5[k];
        const auto h12_m1 = ph12_m1[k];

        pc_41[k] = -fs_66976875_1024 * e_0 * h2_m1 + e_1 * (fs_15400125_512 * h4_m1 + fs_34171875_256 * r_2 * h2_m1) + e_2 * (-fs_39375_88 * h6_m5 - fs_1063125_1936 * h6_m1 - fs_15400125_968 * r_2 * h4_m1 - fs_421875_16 * r_4 * h2_m1) + e_3 * (fs_77175_2288 * h8_m5 - f_105_22 * h8_m1 + fs_1575_22 * r_2 * h6_m5 + fs_42525_484 * r_2 * h6_m1 + fs_820125_968 * r_4 * h4_m1 + fs_421875_484 * r_6 * h2_m1) + e_4 * (-fs_3274425_2712554 * h10_m5 + fs_159137055_70526404 * h10_m1 - fs_77175_51623 * r_2 * h8_m5 + f_210_209 * r_2 * h8_m1 - fs_3150_3179 * r_4 * h6_m5 - fs_42525_34969 * r_4 * h6_m1 - fs_810_121 * r_6 * h4_m1 - fs_421875_81796 * r_8 * h2_m1) + e_5 * (fs_2096325_42204149 * h12_m5 + fs_9882675_1434941066 * h12_m1 + fs_6548850_717470533 * r_2 * h10_m5 - fs_159137055_9327116929 * r_2 * h10_m1 + fs_175_51623 * r_4 * h8_m5 - f_10_209 * r_4 * h8_m1 + fs_1400_1147619 * r_6 * h6_m5 + fs_18900_12623809 * r_6 * h6_m1 + fs_405_69938 * r_8 * h4_m1 + fs_75_20449 * r_10 * h2_m1);
    }

    // NOTE: the rows are formed in 82 loops, as the vectorizer runs out of
    // registers with all 91 of them in one.

#pragma omp simd aligned(pe_2, pe_3, pe_4, pe_5, ph6_m6, ph8_m6, ph10_m6, ph12_m6, ab_2, pc_42 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h6_m6 = ph6_m6[k];
        const auto h8_m6 = ph8_m6[k];
        const auto h10_m6 = ph10_m6[k];
        const auto h12_m6 = ph12_m6[k];

        pc_42[k] = -fs_118125_22 * e_2 * h6_m6 + e_3 * (fs_33075_286 * h8_m6 + fs_9450_11 * r_2 * h6_m6) + e_4 * (-fs_5893965_2712554 * h10_m6 - fs_264600_51623 * r_2 * h8_m6 - fs_37800_3179 * r_4 * h6_m6) + e_5 * (fs_2395800_42204149 * h12_m6 + fs_11787930_717470533 * r_2 * h10_m6 + fs_600_51623 * r_4 * h8_m6 + fs_16800_1147619 * r_6 * h6_m6);
    }

    // NOTE: the rows are formed in 82 loops, as the vectorizer runs out of
    // registers with all 91 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph2_m1, ph4_m1, ph6_m1, ph8_m7, ph8_m1, ph10_m7, ph10_m1, ph12_m7, ph12_m1, ab_2, pc_43 : simd::cache_line_size())
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

        const auto h2_m1 = ph2_m1[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h8_m7 = ph8_m7[k];
        const auto h8_m1 = ph8_m1[k];
        const auto h10_m7 = ph10_m7[k];
        const auto h10_m1 = ph10_m1[k];
        const auto h12_m7 = ph12_m7[k];
        const auto h12_m1 = ph12_m1[k];

        pc_43[k] = fs_218791125_2048 * e_0 * h2_m1 + e_1 * (-fs_2679075_256 * h4_m1 - fs_111628125_512 * r_2 * h2_m1) + e_2 * (-fs_4921875_3872 * h6_m1 + fs_2679075_484 * r_2 * h4_m1 + fs_1378125_32 * r_4 * h2_m1) + e_3 * (fs_33075_1144 * h8_m7 + fs_13395375_163592 * h8_m1 + fs_196875_968 * r_2 * h6_m1 - fs_24111675_81796 * r_4 * h4_m1 - fs_1378125_968 * r_6 * h2_m1) + e_4 * (-fs_392931_319124 * h10_m7 + fs_220172337_141052808 * h10_m1 - fs_66150_51623 * r_2 * h8_m7 - fs_26790750_7382089 * r_2 * h8_m1 - fs_196875_69938 * r_4 * h6_m1 + fs_47628_20449 * r_6 * h4_m1 + fs_1378125_163592 * r_8 * h2_m1) + e_5 * (fs_119790_2221271 * h12_m7 + fs_658845_717470533 * h12_m1 + fs_392931_42204149 * r_2 * h10_m7 - fs_220172337_18654233858 * r_2 * h10_m1 + fs_150_51623 * r_4 * h8_m7 + fs_60750_7382089 * r_4 * h8_m1 + fs_43750_12623809 * r_6 * h6_m1 - fs_11907_5909761 * r_8 * h4_m1 - fs_245_40898 * r_10 * h2_m1);
    }

    // NOTE: the rows are formed in 82 loops, as the vectorizer runs out of
    // registers with all 91 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph2_m2, ph4_m2, ph6_m2, ph8_m8, ph8_m2, ph10_m8, ph10_m2, ph12_m8, ph12_m2, ab_2, pc_44 : simd::cache_line_size())
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

        const auto h2_m2 = ph2_m2[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h8_m8 = ph8_m8[k];
        const auto h8_m2 = ph8_m2[k];
        const auto h10_m8 = ph10_m8[k];
        const auto h10_m2 = ph10_m2[k];
        const auto h12_m8 = ph12_m8[k];
        const auto h12_m2 = ph12_m2[k];

        pc_44[k] = -fs_49116375_1024 * e_0 * h2_m2 + e_1 * (-fs_601425_16 * h4_m2 + fs_25059375_256 * r_2 * h2_m2) + e_2 * (fs_39375_88 * h6_m2 + fs_218700_11 * r_2 * h4_m2 - fs_309375_16 * r_4 * h2_m2) + e_3 * (-fs_33075_208 * h8_m8 + fs_2083725_14872 * h8_m2 - fs_1575_22 * r_2 * h6_m2 - fs_1968300_1859 * r_4 * h4_m2 + fs_28125_44 * r_6 * h2_m2) + e_4 * (fs_11907_319124 * h10_m8 + f_6993_8398 * h10_m2 + fs_33075_4693 * r_2 * h8_m8 - fs_4167450_671099 * r_2 * h8_m2 + fs_3150_3179 * r_4 * h6_m2 + fs_15552_1859 * r_6 * h4_m2 - fs_28125_7436 * r_8 * h2_m2) + e_5 * (fs_87120_2221271 * h12_m8 + fs_152460_717470533 * h12_m2 - fs_11907_42204149 * r_2 * h10_m8 - f_6993_96577 * r_2 * h10_m2 - fs_75_4693 * r_4 * h8_m8 + fs_9450_671099 * r_4 * h8_m2 - fs_1400_1147619 * r_6 * h6_m2 - fs_3888_537251 * r_8 * h4_m2 + fs_5_1859 * r_10 * h2_m2);
    }

    // NOTE: the rows are formed in 82 loops, as the vectorizer runs out of
    // registers with all 91 of them in one.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph4_m3, ph6_m3, ph8_m3, ph10_m9, ph10_m3, ph12_m9, ph12_m3, ab_2, pc_45 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h4_m3 = ph4_m3[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h8_m3 = ph8_m3[k];
        const auto h10_m9 = ph10_m9[k];
        const auto h10_m3 = ph10_m3[k];
        const auto h12_m9 = ph12_m9[k];
        const auto h12_m3 = ph12_m3[k];

        pc_45[k] = fs_12629925_512 * e_1 * h4_m3 + e_2 * (fs_118125_22 * h6_m3 - fs_1148175_88 * r_2 * h4_m3) + e_3 * (fs_231525_2704 * h8_m3 - fs_9450_11 * r_2 * h6_m3 + fs_10333575_14872 * r_4 * h4_m3) + e_4 * (fs_35721_8398 * h10_m9 + fs_214326_1356277 * h10_m3 - fs_231525_61009 * r_2 * h8_m3 + fs_37800_3179 * r_4 * h6_m3 - fs_10206_1859 * r_6 * h4_m3) + e_5 * (fs_38115_2221271 * h12_m9 + fs_38115_1434941066 * h12_m3 - fs_71442_2221271 * r_2 * h10_m9 - fs_857304_717470533 * r_2 * h10_m3 + fs_525_61009 * r_4 * h8_m3 - fs_16800_1147619 * r_6 * h6_m3 + fs_5103_1074502 * r_8 * h4_m3);
    }

    // NOTE: the rows are formed in 82 loops, as the vectorizer runs out of
    // registers with all 91 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph2_0, ph4_0, ph4_p4, ph6_0, ph6_p4, ph8_0, ph8_p4, ph10_0, ph10_p4, ph12_0, ph12_p4, ab_2, pc_46 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;
        const auto r_10 = r_8 * r_2;
        const auto r_12 = r_10 * r_2;

        const auto h2_0 = ph2_0[k];
        const auto h4_0 = ph4_0[k];
        const auto h4_p4 = ph4_p4[k];
        const auto h6_0 = ph6_0[k];
        const auto h6_p4 = ph6_p4[k];
        const auto h8_0 = ph8_0[k];
        const auto h8_p4 = ph8_p4[k];
        const auto h10_0 = ph10_0[k];
        const auto h10_p4 = ph10_p4[k];
        const auto h12_0 = ph12_0[k];
        const auto h12_p4 = ph12_p4[k];

        pc_46[k] = e_0 * (-f_4725_16 * h2_0 - f_10395_16 * r_2) + e_1 * (f_495_16 * h4_0 - fs_3472875_64 * h4_p4 + f_3375_8 * r_2 * h2_0 + f_10395_16 * r_4) + e_2 * (f_75_2 * h6_0 + fs_354375_121 * h6_p4 - f_45_2 * r_2 * h4_0 + fs_3472875_121 * r_2 * h4_p4 - f_375_2 * r_4 * h2_0 - f_495_2 * r_6) + e_3 * (-f_7455_572 * h8_0 - fs_694575_7436 * h8_p4 - f_15 * r_2 * h6_0 - fs_56700_121 * r_2 * h6_p4 + f_135_26 * r_4 * h4_0 - fs_31255875_20449 * r_4 * h4_p4 + f_375_11 * r_6 * h2_0 + f_165_4 * r_8) + e_4 * (f_6615_4199 * h10_0 + fs_2619540_1356277 * h10_p4 + f_7455_2717 * r_2 * h8_0 + fs_2778300_671099 * r_2 * h8_p4 + f_30_17 * r_4 * h6_0 + fs_226800_34969 * r_4 * h6_p4 - f_6_13 * r_6 * h4_0 + fs_246960_20449 * r_6 * h4_p4 - f_375_143 * r_8 * h2_0 - f_3 * r_10) + e_5 * (f_16335_96577 * h12_0 - fs_37733850_717470533 * h12_p4 - f_13230_96577 * r_2 * h10_0 - fs_10478160_717470533 * r_2 * h10_p4 - f_355_2717 * r_4 * h8_0 - fs_6300_671099 * r_4 * h8_p4 - f_20_323 * r_6 * h6_0 - fs_100800_12623809 * r_6 * h6_p4 + f_3_221 * r_8 * h4_0 - fs_61740_5909761 * r_8 * h4_p4 + f_10_143 * r_10 * h2_0 + f_1_13 * r_12) + f_10395_64 * e_6;
    }

    // NOTE: the rows are formed in 82 loops, as the vectorizer runs out of
    // registers with all 91 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph2_p1, ph4_p1, ph4_p3, ph6_p1, ph6_p3, ph8_p1, ph8_p3, ph10_p1, ph10_p3, ph12_p1, ph12_p3, ab_2, pc_47 : simd::cache_line_size())
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

        const auto h2_p1 = ph2_p1[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h4_p3 = ph4_p3[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h8_p3 = ph8_p3[k];
        const auto h10_p1 = ph10_p1[k];
        const auto h10_p3 = ph10_p3[k];
        const auto h12_p1 = ph12_p1[k];
        const auto h12_p3 = ph12_p3[k];

        pc_47[k] = -fs_13395375_512 * e_0 * h2_p1 + e_1 * (f_2385_16 * h4_p1 - fs_694575_256 * h4_p3 + fs_6834375_128 * r_2 * h2_p1) + e_2 * (-fs_590625_242 * h6_p1 + fs_1063125_1936 * h6_p3 - f_2385_22 * r_2 * h4_p1 + fs_694575_484 * r_2 * h4_p3 - fs_84375_8 * r_4 * h2_p1) + e_3 * (fs_40186125_654368 * h8_p1 - fs_2083725_59488 * h8_p3 + fs_47250_121 * r_2 * h6_p1 - fs_42525_484 * r_2 * h6_p3 + f_7155_286 * r_4 * h4_p1 - fs_6251175_81796 * r_4 * h4_p3 + fs_84375_242 * r_6 * h2_p1) + e_4 * (-fs_3143448_17631601 * h10_p1 + fs_6417873_5425108 * h10_p3 - fs_40186125_14764178 * r_2 * h8_p1 + fs_2083725_1342198 * r_2 * h8_p3 - fs_189000_34969 * r_4 * h6_p1 + fs_42525_34969 * r_4 * h6_p3 - f_318_143 * r_6 * h4_p1 + fs_12348_20449 * r_6 * h4_p3 - fs_84375_40898 * r_8 * h2_p1) + e_5 * (-fs_17788815_717470533 * h12_p1 - fs_33960465_717470533 * h12_p3 + fs_12573792_9327116929 * r_2 * h10_p1 - fs_6417873_717470533 * r_2 * h10_p3 + fs_91125_14764178 * r_4 * h8_p1 - fs_4725_1342198 * r_4 * h8_p3 + fs_84000_12623809 * r_6 * h6_p1 - fs_18900_12623809 * r_6 * h6_p3 + f_159_2431 * r_8 * h4_p1 - fs_3087_5909761 * r_8 * h4_p3 + fs_30_20449 * r_10 * h2_p1);
    }

    // NOTE: the rows are formed in 82 loops, as the vectorizer runs out of
    // registers with all 91 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph2_m2, ph4_m2, ph6_m2, ph8_m2, ph10_m2, ph12_m2, ab_2, pc_48 : simd::cache_line_size())
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

        const auto h2_m2 = ph2_m2[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h8_m2 = ph8_m2[k];
        const auto h10_m2 = ph10_m2[k];
        const auto h12_m2 = ph12_m2[k];

        pc_48[k] = fs_31255875_128 * e_0 * h2_m2 + e_1 * (-fs_5145525_128 * h4_m2 - fs_15946875_32 * r_2 * h2_m2) + e_2 * (f_75_2 * h6_m2 + fs_42525_2 * r_2 * h4_m2 + fs_196875_2 * r_4 * h2_m2) + e_3 * (-fs_1620675_327184 * h8_m2 - f_15 * r_2 * h6_m2 - fs_382725_338 * r_4 * h4_m2 - fs_393750_121 * r_6 * h2_m2) + e_4 * (-fs_5501034_17631601 * h10_m2 + fs_1620675_7382089 * r_2 * h8_m2 + f_30_17 * r_4 * h6_m2 + fs_1512_169 * r_6 * h4_m2 + fs_393750_20449 * r_8 * h2_m2) + e_5 * (fs_52827390_717470533 * h12_m2 + fs_22004136_9327116929 * r_2 * h10_m2 - fs_3675_7382089 * r_4 * h8_m2 - f_20_323 * r_6 * h6_m2 - fs_378_48841 * r_8 * h4_m2 - fs_280_20449 * r_10 * h2_m2);
    }

    // NOTE: the rows are formed in 82 loops, as the vectorizer runs out of
    // registers with all 91 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph2_m1, ph4_m3, ph4_m1, ph6_m3, ph6_m1, ph8_m3, ph8_m1, ph10_m3, ph10_m1, ph12_m3, ph12_m1, ab_2, pc_49 : simd::cache_line_size())
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

        const auto h2_m1 = ph2_m1[k];
        const auto h4_m3 = ph4_m3[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h8_m3 = ph8_m3[k];
        const auto h8_m1 = ph8_m1[k];
        const auto h10_m3 = ph10_m3[k];
        const auto h10_m1 = ph10_m1[k];
        const auto h12_m3 = ph12_m3[k];
        const auto h12_m1 = ph12_m1[k];

        pc_49[k] = -fs_13395375_512 * e_0 * h2_m1 + e_1 * (fs_694575_256 * h4_m3 + f_2385_16 * h4_m1 + fs_6834375_128 * r_2 * h2_m1) + e_2 * (-fs_1063125_1936 * h6_m3 - fs_590625_242 * h6_m1 - fs_694575_484 * r_2 * h4_m3 - f_2385_22 * r_2 * h4_m1 - fs_84375_8 * r_4 * h2_m1) + e_3 * (fs_2083725_59488 * h8_m3 + fs_40186125_654368 * h8_m1 + fs_42525_484 * r_2 * h6_m3 + fs_47250_121 * r_2 * h6_m1 + fs_6251175_81796 * r_4 * h4_m3 + f_7155_286 * r_4 * h4_m1 + fs_84375_242 * r_6 * h2_m1) + e_4 * (-fs_6417873_5425108 * h10_m3 - fs_3143448_17631601 * h10_m1 - fs_2083725_1342198 * r_2 * h8_m3 - fs_40186125_14764178 * r_2 * h8_m1 - fs_42525_34969 * r_4 * h6_m3 - fs_189000_34969 * r_4 * h6_m1 - fs_12348_20449 * r_6 * h4_m3 - f_318_143 * r_6 * h4_m1 - fs_84375_40898 * r_8 * h2_m1) + e_5 * (fs_33960465_717470533 * h12_m3 - fs_17788815_717470533 * h12_m1 + fs_6417873_717470533 * r_2 * h10_m3 + fs_12573792_9327116929 * r_2 * h10_m1 + fs_4725_1342198 * r_4 * h8_m3 + fs_91125_14764178 * r_4 * h8_m1 + fs_18900_12623809 * r_6 * h6_m3 + fs_84000_12623809 * r_6 * h6_m1 + fs_3087_5909761 * r_8 * h4_m3 + f_159_2431 * r_8 * h4_m1 + fs_30_20449 * r_10 * h2_m1);
    }

    // NOTE: the rows are formed in 82 loops, as the vectorizer runs out of
    // registers with all 91 of them in one.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph4_m4, ph6_m4, ph8_m4, ph10_m4, ph12_m4, ab_2, pc_50 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h4_m4 = ph4_m4[k];
        const auto h6_m4 = ph6_m4[k];
        const auto h8_m4 = ph8_m4[k];
        const auto h10_m4 = ph10_m4[k];
        const auto h12_m4 = ph12_m4[k];

        pc_50[k] = fs_3472875_64 * e_1 * h4_m4 + e_2 * (-fs_354375_121 * h6_m4 - fs_3472875_121 * r_2 * h4_m4) + e_3 * (fs_694575_7436 * h8_m4 + fs_56700_121 * r_2 * h6_m4 + fs_31255875_20449 * r_4 * h4_m4) + e_4 * (-fs_2619540_1356277 * h10_m4 - fs_2778300_671099 * r_2 * h8_m4 - fs_226800_34969 * r_4 * h6_m4 - fs_246960_20449 * r_6 * h4_m4) + e_5 * (fs_37733850_717470533 * h12_m4 + fs_10478160_717470533 * r_2 * h10_m4 + fs_6300_671099 * r_4 * h8_m4 + fs_100800_12623809 * r_6 * h6_m4 + fs_61740_5909761 * r_8 * h4_m4);
    }

    // NOTE: the rows are formed in 82 loops, as the vectorizer runs out of
    // registers with all 91 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph2_m1, ph4_m1, ph6_m5, ph6_m1, ph8_m5, ph8_m1, ph10_m5, ph10_m1, ph12_m5, ph12_m1, ab_2, pc_51 : simd::cache_line_size())
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

        const auto h2_m1 = ph2_m1[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h6_m5 = ph6_m5[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h8_m5 = ph8_m5[k];
        const auto h8_m1 = ph8_m1[k];
        const auto h10_m5 = ph10_m5[k];
        const auto h10_m1 = ph10_m1[k];
        const auto h12_m5 = ph12_m5[k];
        const auto h12_m1 = ph12_m1[k];

        pc_51[k] = fs_66976875_1024 * e_0 * h2_m1 + e_1 * (-fs_15400125_512 * h4_m1 - fs_34171875_256 * r_2 * h2_m1) + e_2 * (-fs_39375_88 * h6_m5 + fs_1063125_1936 * h6_m1 + fs_15400125_968 * r_2 * h4_m1 + fs_421875_16 * r_4 * h2_m1) + e_3 * (fs_77175_2288 * h8_m5 + f_105_22 * h8_m1 + fs_1575_22 * r_2 * h6_m5 - fs_42525_484 * r_2 * h6_m1 - fs_820125_968 * r_4 * h4_m1 - fs_421875_484 * r_6 * h2_m1) + e_4 * (-fs_3274425_2712554 * h10_m5 - fs_159137055_70526404 * h10_m1 - fs_77175_51623 * r_2 * h8_m5 - f_210_209 * r_2 * h8_m1 - fs_3150_3179 * r_4 * h6_m5 + fs_42525_34969 * r_4 * h6_m1 + fs_810_121 * r_6 * h4_m1 + fs_421875_81796 * r_8 * h2_m1) + e_5 * (fs_2096325_42204149 * h12_m5 - fs_9882675_1434941066 * h12_m1 + fs_6548850_717470533 * r_2 * h10_m5 + fs_159137055_9327116929 * r_2 * h10_m1 + fs_175_51623 * r_4 * h8_m5 + f_10_209 * r_4 * h8_m1 + fs_1400_1147619 * r_6 * h6_m5 - fs_18900_12623809 * r_6 * h6_m1 - fs_405_69938 * r_8 * h4_m1 - fs_75_20449 * r_10 * h2_m1);
    }

    // NOTE: the rows are formed in 82 loops, as the vectorizer runs out of
    // registers with all 91 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph2_m2, ph4_m2, ph6_m6, ph6_m2, ph8_m6, ph8_m2, ph10_m6, ph10_m2, ph12_m6, ph12_m2, ab_2, pc_52 : simd::cache_line_size())
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

        const auto h2_m2 = ph2_m2[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h6_m6 = ph6_m6[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h8_m6 = ph8_m6[k];
        const auto h8_m2 = ph8_m2[k];
        const auto h10_m6 = ph10_m6[k];
        const auto h10_m2 = ph10_m2[k];
        const auto h12_m6 = ph12_m6[k];
        const auto h12_m2 = ph12_m2[k];

        pc_52[k] = -fs_40186125_512 * e_0 * h2_m2 + e_1 * (-fs_3213675_512 * h4_m2 + fs_20503125_128 * r_2 * h2_m2) + e_2 * (fs_196875_44 * h6_m6 + fs_354375_121 * h6_m2 + fs_3213675_968 * r_2 * h4_m2 - fs_253125_8 * r_4 * h2_m2) + e_3 * (-fs_55125_2288 * h8_m6 - fs_231525_20449 * h8_m2 - fs_7875_11 * r_2 * h6_m6 - fs_56700_121 * r_2 * h6_m2 - fs_28923075_163592 * r_4 * h4_m2 + fs_253125_242 * r_6 * h2_m2) + e_4 * (-fs_43659_1356277 * h10_m6 - fs_80725491_35263202 * h10_m2 + fs_55125_51623 * r_2 * h8_m6 + fs_3704400_7382089 * r_2 * h8_m2 + fs_31500_3179 * r_4 * h6_m6 + fs_226800_34969 * r_4 * h6_m2 + fs_28566_20449 * r_6 * h4_m2 - fs_253125_40898 * r_8 * h2_m2) + e_5 * (fs_1617165_42204149 * h12_m6 - fs_3773385_1434941066 * h12_m2 + fs_174636_717470533 * r_2 * h10_m6 + fs_161450982_9327116929 * r_2 * h10_m2 - fs_125_51623 * r_4 * h8_m6 - fs_8400_7382089 * r_4 * h8_m2 - fs_14000_1147619 * r_6 * h6_m6 - fs_100800_12623809 * r_6 * h6_m2 - fs_14283_11819522 * r_8 * h4_m2 + fs_90_20449 * r_10 * h2_m2);
    }

    // NOTE: the rows are formed in 82 loops, as the vectorizer runs out of
    // registers with all 91 of them in one.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph4_m3, ph6_m3, ph8_m7, ph8_m3, ph10_m7, ph10_m3, ph12_m7, ph12_m3, ab_2, pc_53 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h4_m3 = ph4_m3[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h8_m7 = ph8_m7[k];
        const auto h8_m3 = ph8_m3[k];
        const auto h10_m7 = ph10_m7[k];
        const auto h10_m3 = ph10_m3[k];
        const auto h12_m7 = ph12_m7[k];
        const auto h12_m3 = ph12_m3[k];

        pc_53[k] = fs_22920975_512 * e_1 * h4_m3 + e_2 * (fs_39375_88 * h6_m3 - fs_2083725_88 * r_2 * h4_m3) + e_3 * (-fs_33075_208 * h8_m7 - fs_77175_676 * h8_m3 - fs_1575_22 * r_2 * h6_m3 + fs_18753525_14872 * r_4 * h4_m3) + e_4 * (fs_194481_159562 * h10_m7 - fs_3814209_2712554 * h10_m3 + fs_33075_4693 * r_2 * h8_m7 + fs_308700_61009 * r_2 * h8_m3 + fs_3150_3179 * r_4 * h6_m3 - fs_18522_1859 * r_6 * h4_m3) + e_5 * (fs_49005_2221271 * h12_m7 - fs_1029105_1434941066 * h12_m3 - fs_388962_42204149 * r_2 * h10_m7 + fs_7628418_717470533 * r_2 * h10_m3 - fs_75_4693 * r_4 * h8_m7 - fs_700_61009 * r_4 * h8_m3 - fs_1400_1147619 * r_6 * h6_m3 + fs_9261_1074502 * r_8 * h4_m3);
    }

    // NOTE: the rows are formed in 82 loops, as the vectorizer runs out of
    // registers with all 91 of them in one.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph4_m4, ph6_m4, ph8_m8, ph8_m4, ph10_m8, ph10_m4, ph12_m8, ph12_m4, ab_2, pc_54 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h4_m4 = ph4_m4[k];
        const auto h6_m4 = ph6_m4[k];
        const auto h8_m8 = ph8_m8[k];
        const auto h8_m4 = ph8_m4[k];
        const auto h10_m8 = ph10_m8[k];
        const auto h10_m4 = ph10_m4[k];
        const auto h12_m8 = ph12_m8[k];
        const auto h12_m4 = ph12_m4[k];

        pc_54[k] = -fs_1403325_256 * e_1 * h4_m4 + e_2 * (-fs_196875_44 * h6_m4 + fs_127575_44 * r_2 * h4_m4) + e_3 * (fs_11025_208 * h8_m8 - fs_385875_2704 * h8_m4 + fs_7875_11 * r_2 * h6_m4 - fs_1148175_7436 * r_4 * h4_m4) + e_4 * (fs_321489_79781 * h10_m8 - fs_583443_1356277 * h10_m4 - fs_11025_4693 * r_2 * h8_m8 + fs_385875_61009 * r_2 * h8_m4 - fs_31500_3179 * r_4 * h6_m4 + fs_2268_1859 * r_6 * h4_m4) + e_5 * (fs_16335_2221271 * h12_m8 - fs_76230_717470533 * h12_m4 - fs_1285956_42204149 * r_2 * h10_m8 + fs_2333772_717470533 * r_2 * h10_m4 + fs_25_4693 * r_4 * h8_m8 - fs_875_61009 * r_4 * h8_m4 + fs_14000_1147619 * r_6 * h6_m4 - fs_567_537251 * r_8 * h4_m4);
    }

    // NOTE: the rows are formed in 82 loops, as the vectorizer runs out of
    // registers with all 91 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph2_0, ph2_p2, ph4_0, ph4_p2, ph6_0, ph6_p2, ph8_0, ph8_p2, ph10_0, ph10_p2, ph12_0, ph12_p2, ab_2, pc_55 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;
        const auto r_10 = r_8 * r_2;
        const auto r_12 = r_10 * r_2;

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
        const auto h12_0 = ph12_0[k];
        const auto h12_p2 = ph12_p2[k];

        pc_55[k] = e_0 * (-f_12285_32 * h2_0 + fs_131274675_1024 * h2_p2 - f_10395_16 * r_2) + e_1 * (f_180 * h4_0 - fs_496125_16 * h4_p2 + f_8775_16 * r_2 * h2_0 - fs_66976875_256 * r_2 * h2_p2 + f_10395_16 * r_4) + e_2 * (-f_375_11 * h6_0 + fs_590625_242 * h6_p2 - f_1440_11 * r_2 * h4_0 + fs_1984500_121 * r_2 * h4_p2 - f_975_4 * r_4 * h2_0 + fs_826875_16 * r_4 * h2_p2 - f_495_2 * r_6) + e_3 * (f_525_286 * h8_0 - fs_3472875_40898 * h8_p2 + f_150_11 * r_2 * h6_0 - fs_47250_121 * r_2 * h6_p2 + f_4320_143 * r_4 * h4_0 - fs_17860500_20449 * r_4 * h4_p2 + f_975_22 * r_6 * h2_0 - fs_826875_484 * r_6 * h2_p2 + f_165_4 * r_8) + e_4 * (f_189_323 * h10_0 + fs_32089365_17631601 * h10_p2 - f_1050_2717 * r_2 * h8_0 + fs_27783000_7382089 * r_2 * h8_p2 - f_300_187 * r_4 * h6_0 + fs_189000_34969 * r_4 * h6_p2 - f_384_143 * r_6 * h4_0 + fs_141120_20449 * r_6 * h4_p2 - f_75_22 * r_8 * h2_0 + fs_826875_81796 * r_8 * h2_p2 - f_3 * r_10) + e_5 * (-f_26136_96577 * h12_0 - fs_36224496_717470533 * h12_p2 - f_378_7429 * r_2 * h10_0 - fs_128357460_9327116929 * r_2 * h10_p2 + f_50_2717 * r_4 * h8_0 - fs_63000_7382089 * r_4 * h8_p2 + f_200_3553 * r_6 * h6_0 - fs_84000_12623809 * r_6 * h6_p2 + f_192_2431 * r_8 * h4_0 - fs_35280_5909761 * r_8 * h4_p2 + f_1_11 * r_10 * h2_0 - fs_147_20449 * r_10 * h2_p2 + f_1_13 * r_12) + f_10395_64 * e_6;
    }

    // NOTE: the rows are formed in 82 loops, as the vectorizer runs out of
    // registers with all 91 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph2_m1, ph4_m1, ph6_m1, ph8_m1, ph10_m1, ph12_m1, ab_2, pc_56 : simd::cache_line_size())
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

        const auto h2_m1 = ph2_m1[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h8_m1 = ph8_m1[k];
        const auto h10_m1 = ph10_m1[k];
        const auto h12_m1 = ph12_m1[k];

        pc_56[k] = -fs_6251175_1024 * e_0 * h2_m1 + e_1 * (fs_212625_32 * h4_m1 + fs_3189375_256 * r_2 * h2_m1) + e_2 * (-f_375_11 * h6_m1 - fs_425250_121 * r_2 * h4_m1 - fs_39375_16 * r_4 * h2_m1) + e_3 * (fs_5788125_81796 * h8_m1 + f_150_11 * r_2 * h6_m1 + fs_3827250_20449 * r_4 * h4_m1 + fs_39375_484 * r_6 * h2_m1) + e_4 * (-fs_41257755_17631601 * h10_m1 - fs_23152500_7382089 * r_2 * h8_m1 - f_300_187 * r_4 * h6_m1 - fs_30240_20449 * r_6 * h4_m1 - fs_39375_81796 * r_8 * h2_m1) + e_5 * (fs_66411576_717470533 * h12_m1 + fs_165031020_9327116929 * r_2 * h10_m1 + fs_52500_7382089 * r_4 * h8_m1 + f_200_3553 * r_6 * h6_m1 + fs_7560_5909761 * r_8 * h4_m1 + fs_7_20449 * r_10 * h2_m1);
    }

    // NOTE: the rows are formed in 82 loops, as the vectorizer runs out of
    // registers with all 91 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph2_m2, ph4_m2, ph6_m2, ph8_m2, ph10_m2, ph12_m2, ab_2, pc_57 : simd::cache_line_size())
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

        const auto h2_m2 = ph2_m2[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h8_m2 = ph8_m2[k];
        const auto h10_m2 = ph10_m2[k];
        const auto h12_m2 = ph12_m2[k];

        pc_57[k] = -fs_131274675_1024 * e_0 * h2_m2 + e_1 * (fs_496125_16 * h4_m2 + fs_66976875_256 * r_2 * h2_m2) + e_2 * (-fs_590625_242 * h6_m2 - fs_1984500_121 * r_2 * h4_m2 - fs_826875_16 * r_4 * h2_m2) + e_3 * (fs_3472875_40898 * h8_m2 + fs_47250_121 * r_2 * h6_m2 + fs_17860500_20449 * r_4 * h4_m2 + fs_826875_484 * r_6 * h2_m2) + e_4 * (-fs_32089365_17631601 * h10_m2 - fs_27783000_7382089 * r_2 * h8_m2 - fs_189000_34969 * r_4 * h6_m2 - fs_141120_20449 * r_6 * h4_m2 - fs_826875_81796 * r_8 * h2_m2) + e_5 * (fs_36224496_717470533 * h12_m2 + fs_128357460_9327116929 * r_2 * h10_m2 + fs_63000_7382089 * r_4 * h8_m2 + fs_84000_12623809 * r_6 * h6_m2 + fs_35280_5909761 * r_8 * h4_m2 + fs_147_20449 * r_10 * h2_m2);
    }

    // NOTE: the rows are formed in 82 loops, as the vectorizer runs out of
    // registers with all 91 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph2_m1, ph4_m3, ph4_m1, ph6_m3, ph6_m1, ph8_m3, ph8_m1, ph10_m3, ph10_m1, ph12_m3, ph12_m1, ab_2, pc_58 : simd::cache_line_size())
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

        const auto h2_m1 = ph2_m1[k];
        const auto h4_m3 = ph4_m3[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h8_m3 = ph8_m3[k];
        const auto h8_m1 = ph8_m1[k];
        const auto h10_m3 = ph10_m3[k];
        const auto h10_m1 = ph10_m1[k];
        const auto h12_m3 = ph12_m3[k];
        const auto h12_m1 = ph12_m1[k];

        pc_58[k] = fs_13395375_512 * e_0 * h2_m1 + e_1 * (fs_694575_256 * h4_m3 - f_2385_16 * h4_m1 - fs_6834375_128 * r_2 * h2_m1) + e_2 * (-fs_1063125_1936 * h6_m3 + fs_590625_242 * h6_m1 - fs_694575_484 * r_2 * h4_m3 + f_2385_22 * r_2 * h4_m1 + fs_84375_8 * r_4 * h2_m1) + e_3 * (fs_2083725_59488 * h8_m3 - fs_40186125_654368 * h8_m1 + fs_42525_484 * r_2 * h6_m3 - fs_47250_121 * r_2 * h6_m1 + fs_6251175_81796 * r_4 * h4_m3 - f_7155_286 * r_4 * h4_m1 - fs_84375_242 * r_6 * h2_m1) + e_4 * (-fs_6417873_5425108 * h10_m3 + fs_3143448_17631601 * h10_m1 - fs_2083725_1342198 * r_2 * h8_m3 + fs_40186125_14764178 * r_2 * h8_m1 - fs_42525_34969 * r_4 * h6_m3 + fs_189000_34969 * r_4 * h6_m1 - fs_12348_20449 * r_6 * h4_m3 + f_318_143 * r_6 * h4_m1 + fs_84375_40898 * r_8 * h2_m1) + e_5 * (fs_33960465_717470533 * h12_m3 + fs_17788815_717470533 * h12_m1 + fs_6417873_717470533 * r_2 * h10_m3 - fs_12573792_9327116929 * r_2 * h10_m1 + fs_4725_1342198 * r_4 * h8_m3 - fs_91125_14764178 * r_4 * h8_m1 + fs_18900_12623809 * r_6 * h6_m3 - fs_84000_12623809 * r_6 * h6_m1 + fs_3087_5909761 * r_8 * h4_m3 - f_159_2431 * r_8 * h4_m1 - fs_30_20449 * r_10 * h2_m1);
    }

    // NOTE: the rows are formed in 82 loops, as the vectorizer runs out of
    // registers with all 91 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph2_m2, ph4_m4, ph4_m2, ph6_m4, ph6_m2, ph8_m4, ph8_m2, ph10_m4, ph10_m2, ph12_m4, ph12_m2, ab_2, pc_59 : simd::cache_line_size())
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

        const auto h2_m2 = ph2_m2[k];
        const auto h4_m4 = ph4_m4[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h6_m4 = ph6_m4[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h8_m4 = ph8_m4[k];
        const auto h8_m2 = ph8_m2[k];
        const auto h10_m4 = ph10_m4[k];
        const auto h10_m2 = ph10_m2[k];
        const auto h12_m4 = ph12_m4[k];
        const auto h12_m2 = ph12_m2[k];

        pc_59[k] = -fs_13395375_128 * e_0 * h2_m2 + e_1 * (-fs_6251175_128 * h4_m4 + fs_18225_8 * h4_m2 + fs_6834375_32 * r_2 * h2_m2) + e_2 * (fs_4921875_3872 * h6_m4 + fs_1063125_1936 * h6_m2 + fs_6251175_242 * r_2 * h4_m4 - fs_145800_121 * r_2 * h4_m2 - fs_84375_2 * r_4 * h2_m2) + e_3 * (-fs_385875_59488 * h8_m4 - fs_27860175_327184 * h8_m2 - fs_196875_968 * r_2 * h6_m4 - fs_42525_484 * r_2 * h6_m2 - fs_56260575_40898 * r_4 * h4_m4 + fs_1312200_20449 * r_4 * h4_m2 + fs_168750_121 * r_6 * h2_m2) + e_4 * (-fs_1178793_10850216 * h10_m4 + fs_130977_97682 * h10_m2 + fs_385875_1342198 * r_2 * h8_m4 + fs_77175_20449 * r_2 * h8_m2 + fs_196875_69938 * r_4 * h6_m4 + fs_42525_34969 * r_4 * h6_m2 + fs_222264_20449 * r_6 * h4_m4 - fs_10368_20449 * r_6 * h4_m2 - fs_168750_20449 * r_8 * h2_m2) + e_5 * (fs_26832960_717470533 * h12_m4 + fs_10062360_717470533 * h12_m2 + fs_1178793_1434941066 * r_2 * h10_m4 - fs_261954_25836889 * r_2 * h10_m2 - fs_875_1342198 * r_4 * h8_m4 - fs_175_20449 * r_4 * h8_m2 - fs_43750_12623809 * r_6 * h6_m4 - fs_18900_12623809 * r_6 * h6_m2 - fs_55566_5909761 * r_8 * h4_m4 + fs_2592_5909761 * r_8 * h4_m2 + fs_120_20449 * r_10 * h2_m2);
    }

    // NOTE: the rows are formed in 82 loops, as the vectorizer runs out of
    // registers with all 91 of them in one.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph4_m3, ph6_m5, ph6_m3, ph8_m5, ph8_m3, ph10_m5, ph10_m3, ph12_m5, ph12_m3, ab_2, pc_60 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h4_m3 = ph4_m3[k];
        const auto h6_m5 = ph6_m5[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h8_m5 = ph8_m5[k];
        const auto h8_m3 = ph8_m3[k];
        const auto h10_m5 = ph10_m5[k];
        const auto h10_m3 = ph10_m3[k];
        const auto h12_m5 = ph12_m5[k];
        const auto h12_m3 = ph12_m3[k];

        pc_60[k] = fs_5315625_128 * e_1 * h4_m3 + e_2 * (fs_1063125_352 * h6_m5 - fs_4921875_3872 * h6_m3 - fs_5315625_242 * r_2 * h4_m3) + e_3 * (-fs_231525_2288 * h8_m5 - fs_385875_29744 * h8_m3 - fs_42525_88 * r_2 * h6_m5 + fs_196875_968 * r_2 * h6_m3 + fs_47840625_40898 * r_4 * h4_m3) + e_4 * (fs_5282739_10850216 * h10_m5 + fs_26413695_10850216 * h10_m3 + fs_231525_51623 * r_2 * h8_m5 + fs_385875_671099 * r_2 * h8_m3 + fs_42525_6358 * r_4 * h6_m5 - fs_196875_69938 * r_4 * h6_m3 - fs_189000_20449 * r_6 * h4_m3) + e_5 * (fs_1006236_42204149 * h12_m5 + fs_4528062_717470533 * h12_m3 - fs_5282739_1434941066 * r_2 * h10_m5 - fs_26413695_1434941066 * r_2 * h10_m3 - fs_525_51623 * r_4 * h8_m5 - fs_875_671099 * r_4 * h8_m3 - fs_9450_1147619 * r_6 * h6_m5 + fs_43750_12623809 * r_6 * h6_m3 + fs_47250_5909761 * r_8 * h4_m3);
    }

    // NOTE: the rows are formed in 82 loops, as the vectorizer runs out of
    // registers with all 91 of them in one.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph4_m4, ph6_m6, ph6_m4, ph8_m6, ph8_m4, ph10_m6, ph10_m4, ph12_m6, ph12_m4, ab_2, pc_61 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h4_m4 = ph4_m4[k];
        const auto h6_m6 = ph6_m6[k];
        const auto h6_m4 = ph6_m4[k];
        const auto h8_m6 = ph8_m6[k];
        const auto h8_m4 = ph8_m4[k];
        const auto h10_m6 = ph10_m6[k];
        const auto h10_m4 = ph10_m4[k];
        const auto h12_m6 = ph12_m6[k];
        const auto h12_m4 = ph12_m4[k];

        pc_61[k] = -fs_2338875_128 * e_1 * h4_m4 + e_2 * (-fs_39375_16 * h6_m6 - fs_1063125_352 * h6_m4 + fs_212625_22 * r_2 * h4_m4) + e_3 * (-fs_11025_208 * h8_m6 + fs_231525_5408 * h8_m4 + fs_1575_4 * r_2 * h6_m6 + fs_42525_88 * r_2 * h6_m4 - fs_1913625_3718 * r_4 * h4_m4) + e_4 * (fs_257985_104329 * h10_m6 + fs_24310125_10850216 * h10_m4 + fs_11025_4693 * r_2 * h8_m6 - fs_231525_122018 * r_2 * h8_m4 - fs_1575_289 * r_4 * h6_m6 - fs_42525_6358 * r_4 * h6_m4 + fs_7560_1859 * r_6 * h4_m4) + e_5 * (fs_470448_42204149 * h12_m6 + fs_1463616_717470533 * h12_m4 - fs_1031940_55190041 * r_2 * h10_m6 - fs_24310125_1434941066 * r_2 * h10_m4 - fs_25_4693 * r_4 * h8_m6 + fs_525_122018 * r_4 * h8_m4 + fs_700_104329 * r_6 * h6_m6 + fs_9450_1147619 * r_6 * h6_m4 - fs_1890_537251 * r_8 * h4_m4);
    }

    // NOTE: the rows are formed in 82 loops, as the vectorizer runs out of
    // registers with all 91 of them in one.

#pragma omp simd aligned(pe_2, pe_3, pe_4, pe_5, ph6_m5, ph8_m7, ph8_m5, ph10_m7, ph10_m5, ph12_m7, ph12_m5, ab_2, pc_62 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h6_m5 = ph6_m5[k];
        const auto h8_m7 = ph8_m7[k];
        const auto h8_m5 = ph8_m5[k];
        const auto h10_m7 = ph10_m7[k];
        const auto h10_m5 = ph10_m5[k];
        const auto h12_m7 = ph12_m7[k];
        const auto h12_m5 = ph12_m5[k];

        pc_62[k] = fs_39375_16 * e_2 * h6_m5 + e_3 * (fs_55125_416 * h8_m7 + fs_77175_416 * h8_m5 - fs_1575_4 * r_2 * h6_m5) + e_4 * (fs_238140_79781 * h10_m7 + fs_5250987_5425108 * h10_m5 - fs_55125_9386 * r_2 * h8_m7 - fs_77175_9386 * r_2 * h8_m5 + fs_1575_289 * r_4 * h6_m5) + e_5 * (fs_6534_2221271 * h12_m7 + fs_15246_42204149 * h12_m5 - fs_952560_42204149 * r_2 * h10_m7 - fs_5250987_717470533 * r_2 * h10_m5 + fs_125_9386 * r_4 * h8_m7 + fs_175_9386 * r_4 * h8_m5 - fs_700_104329 * r_6 * h6_m5);
    }

    // NOTE: the rows are formed in 82 loops, as the vectorizer runs out of
    // registers with all 91 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph2_0, ph4_0, ph6_0, ph8_0, ph10_0, ph12_0, ab_2, pc_63 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;
        const auto r_10 = r_8 * r_2;
        const auto r_12 = r_10 * r_2;

        const auto h2_0 = ph2_0[k];
        const auto h4_0 = ph4_0[k];
        const auto h6_0 = ph6_0[k];
        const auto h8_0 = ph8_0[k];
        const auto h10_0 = ph10_0[k];
        const auto h12_0 = ph12_0[k];

        pc_63[k] = e_0 * (-f_6615_16 * h2_0 - f_10395_16 * r_2) + e_1 * (f_945_4 * h4_0 + f_4725_8 * r_2 * h2_0 + f_10395_16 * r_4) + e_2 * (-f_750_11 * h6_0 - f_1890_11 * r_2 * h4_0 - f_525_2 * r_4 * h2_0 - f_495_2 * r_6) + e_3 * (f_3675_286 * h8_0 + f_300_11 * r_2 * h6_0 + f_5670_143 * r_4 * h4_0 + f_525_11 * r_6 * h2_0 + f_165_4 * r_8) + e_4 * (-f_7938_4199 * h10_0 - f_7350_2717 * r_2 * h8_0 - f_600_187 * r_4 * h6_0 - f_504_143 * r_6 * h4_0 - f_525_143 * r_8 * h2_0 - f_3 * r_10) + e_5 * (f_30492_96577 * h12_0 + f_15876_96577 * r_2 * h10_0 + f_350_2717 * r_4 * h8_0 + f_400_3553 * r_6 * h6_0 + f_252_2431 * r_8 * h4_0 + f_14_143 * r_10 * h2_0 + f_1_13 * r_12) + f_10395_64 * e_6;
    }

    // NOTE: the rows are formed in 82 loops, as the vectorizer runs out of
    // registers with all 91 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph2_p1, ph4_p1, ph6_p1, ph8_p1, ph10_p1, ph12_p1, ab_2, pc_64 : simd::cache_line_size())
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

        const auto h2_p1 = ph2_p1[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h10_p1 = ph10_p1[k];
        const auto h12_p1 = ph12_p1[k];

        pc_64[k] = -fs_6251175_1024 * e_0 * h2_p1 + e_1 * (fs_212625_32 * h4_p1 + fs_3189375_256 * r_2 * h2_p1) + e_2 * (-f_375_11 * h6_p1 - fs_425250_121 * r_2 * h4_p1 - fs_39375_16 * r_4 * h2_p1) + e_3 * (fs_5788125_81796 * h8_p1 + f_150_11 * r_2 * h6_p1 + fs_3827250_20449 * r_4 * h4_p1 + fs_39375_484 * r_6 * h2_p1) + e_4 * (-fs_41257755_17631601 * h10_p1 - fs_23152500_7382089 * r_2 * h8_p1 - f_300_187 * r_4 * h6_p1 - fs_30240_20449 * r_6 * h4_p1 - fs_39375_81796 * r_8 * h2_p1) + e_5 * (fs_66411576_717470533 * h12_p1 + fs_165031020_9327116929 * r_2 * h10_p1 + fs_52500_7382089 * r_4 * h8_p1 + f_200_3553 * r_6 * h6_p1 + fs_7560_5909761 * r_8 * h4_p1 + fs_7_20449 * r_10 * h2_p1);
    }

    // NOTE: the rows are formed in 82 loops, as the vectorizer runs out of
    // registers with all 91 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph2_p2, ph4_p2, ph4_p3, ph6_p2, ph6_p3, ph8_p2, ph8_p3, ph10_p2, ph10_p3, ph12_p2, ph12_p3, ab_2, pc_65, pc_66 : simd::cache_line_size())
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

        const auto h2_p2 = ph2_p2[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h4_p3 = ph4_p3[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h8_p2 = ph8_p2[k];
        const auto h8_p3 = ph8_p3[k];
        const auto h10_p2 = ph10_p2[k];
        const auto h10_p3 = ph10_p3[k];
        const auto h12_p2 = ph12_p2[k];
        const auto h12_p3 = ph12_p3[k];

        pc_65[k] = fs_31255875_128 * e_0 * h2_p2 + e_1 * (-fs_5145525_128 * h4_p2 - fs_15946875_32 * r_2 * h2_p2) + e_2 * (f_75_2 * h6_p2 + fs_42525_2 * r_2 * h4_p2 + fs_196875_2 * r_4 * h2_p2) + e_3 * (-fs_1620675_327184 * h8_p2 - f_15 * r_2 * h6_p2 - fs_382725_338 * r_4 * h4_p2 - fs_393750_121 * r_6 * h2_p2) + e_4 * (-fs_5501034_17631601 * h10_p2 + fs_1620675_7382089 * r_2 * h8_p2 + f_30_17 * r_4 * h6_p2 + fs_1512_169 * r_6 * h4_p2 + fs_393750_20449 * r_8 * h2_p2) + e_5 * (fs_52827390_717470533 * h12_p2 + fs_22004136_9327116929 * r_2 * h10_p2 - fs_3675_7382089 * r_4 * h8_p2 - f_20_323 * r_6 * h6_p2 - fs_378_48841 * r_8 * h4_p2 - fs_280_20449 * r_10 * h2_p2);

        pc_66[k] = -fs_2679075_64 * e_1 * h4_p3 + e_2 * (f_3225_44 * h6_p3 + fs_2679075_121 * r_2 * h4_p3) + e_3 * (-fs_540225_3718 * h8_p3 - f_645_22 * r_2 * h6_p3 - fs_24111675_20449 * r_4 * h4_p3) + e_4 * (fs_2750517_5425108 * h10_p3 + fs_4321800_671099 * r_2 * h8_p3 + f_645_187 * r_4 * h6_p3 + fs_190512_20449 * r_6 * h4_p3) + e_5 * (fs_35218260_717470533 * h12_p3 - fs_2750517_717470533 * r_2 * h10_p3 - fs_9800_671099 * r_4 * h8_p3 - f_430_3553 * r_6 * h6_p3 - fs_47628_5909761 * r_8 * h4_p3);
    }

    // NOTE: the rows are formed in 82 loops, as the vectorizer runs out of
    // registers with all 91 of them in one.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph4_p4, ph6_p4, ph6_p5, ph6_p6, ph8_p4, ph8_p6, ph10_p4, ph10_p5, ph10_p6, ph12_p4, ph12_p5, ph12_p6, ab_2, pc_67, pc_68, pc_69 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h4_p4 = ph4_p4[k];
        const auto h6_p4 = ph6_p4[k];
        const auto h6_p5 = ph6_p5[k];
        const auto h6_p6 = ph6_p6[k];
        const auto h8_p4 = ph8_p4[k];
        const auto h8_p6 = ph8_p6[k];
        const auto h10_p4 = ph10_p4[k];
        const auto h10_p5 = ph10_p5[k];
        const auto h10_p6 = ph10_p6[k];
        const auto h12_p4 = ph12_p4[k];
        const auto h12_p5 = ph12_p5[k];
        const auto h12_p6 = ph12_p6[k];

        pc_67[k] = fs_4465125_64 * e_1 * h4_p4 + e_2 * (f_150_11 * h6_p4 - fs_4465125_121 * r_2 * h4_p4) + e_3 * (-fs_4862025_29744 * h8_p4 - f_60_11 * r_2 * h6_p4 + fs_40186125_20449 * r_4 * h4_p4) + e_4 * (fs_4584195_1356277 * h10_p4 + fs_4862025_671099 * r_2 * h8_p4 + f_120_187 * r_4 * h6_p4 - fs_317520_20449 * r_6 * h4_p4) + e_5 * (fs_18783072_717470533 * h12_p4 - fs_18336780_717470533 * r_2 * h10_p4 - fs_11025_671099 * r_4 * h8_p4 - f_80_3553 * r_6 * h6_p4 + fs_79380_5909761 * r_8 * h4_p4);

        pc_68[k] = -f_375_4 * e_2 * h6_p5 + f_75_2 * e_3 * r_2 * h6_p5 + e_4 * (fs_83349_15028 * h10_p5 - f_75_17 * r_4 * h6_p5) + e_5 * (fs_426888_42204149 * h12_p5 - fs_83349_1987453 * r_2 * h10_p5 + f_50_323 * r_6 * h6_p5);

        pc_69[k] = f_75_2 * e_2 * h6_p6 + e_3 * (fs_77175_208 * h8_p6 - f_15 * r_2 * h6_p6) + e_4 * (fs_5000940_1356277 * h10_p6 - fs_77175_4693 * r_2 * h8_p6 + f_30_17 * r_4 * h6_p6) + e_5 * (fs_91476_42204149 * h12_p6 - fs_20003760_717470533 * r_2 * h10_p6 + fs_175_4693 * r_4 * h8_p6 - f_20_323 * r_6 * h6_p6);
    }

    // NOTE: the rows are formed in 82 loops, as the vectorizer runs out of
    // registers with all 91 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph2_0, ph2_p2, ph4_0, ph4_p2, ph6_0, ph6_p2, ph8_0, ph8_p2, ph10_0, ph10_p2, ph12_0, ph12_p2, ab_2, pc_70 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;
        const auto r_10 = r_8 * r_2;
        const auto r_12 = r_10 * r_2;

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
        const auto h12_0 = ph12_0[k];
        const auto h12_p2 = ph12_p2[k];

        pc_70[k] = e_0 * (-f_12285_32 * h2_0 - fs_131274675_1024 * h2_p2 - f_10395_16 * r_2) + e_1 * (f_180 * h4_0 + fs_496125_16 * h4_p2 + f_8775_16 * r_2 * h2_0 + fs_66976875_256 * r_2 * h2_p2 + f_10395_16 * r_4) + e_2 * (-f_375_11 * h6_0 - fs_590625_242 * h6_p2 - f_1440_11 * r_2 * h4_0 - fs_1984500_121 * r_2 * h4_p2 - f_975_4 * r_4 * h2_0 - fs_826875_16 * r_4 * h2_p2 - f_495_2 * r_6) + e_3 * (f_525_286 * h8_0 + fs_3472875_40898 * h8_p2 + f_150_11 * r_2 * h6_0 + fs_47250_121 * r_2 * h6_p2 + f_4320_143 * r_4 * h4_0 + fs_17860500_20449 * r_4 * h4_p2 + f_975_22 * r_6 * h2_0 + fs_826875_484 * r_6 * h2_p2 + f_165_4 * r_8) + e_4 * (f_189_323 * h10_0 - fs_32089365_17631601 * h10_p2 - f_1050_2717 * r_2 * h8_0 - fs_27783000_7382089 * r_2 * h8_p2 - f_300_187 * r_4 * h6_0 - fs_189000_34969 * r_4 * h6_p2 - f_384_143 * r_6 * h4_0 - fs_141120_20449 * r_6 * h4_p2 - f_75_22 * r_8 * h2_0 - fs_826875_81796 * r_8 * h2_p2 - f_3 * r_10) + e_5 * (-f_26136_96577 * h12_0 + fs_36224496_717470533 * h12_p2 - f_378_7429 * r_2 * h10_0 + fs_128357460_9327116929 * r_2 * h10_p2 + f_50_2717 * r_4 * h8_0 + fs_63000_7382089 * r_4 * h8_p2 + f_200_3553 * r_6 * h6_0 + fs_84000_12623809 * r_6 * h6_p2 + f_192_2431 * r_8 * h4_0 + fs_35280_5909761 * r_8 * h4_p2 + f_1_11 * r_10 * h2_0 + fs_147_20449 * r_10 * h2_p2 + f_1_13 * r_12) + f_10395_64 * e_6;
    }

    // NOTE: the rows are formed in 82 loops, as the vectorizer runs out of
    // registers with all 91 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph2_p1, ph4_p1, ph4_p3, ph6_p1, ph6_p3, ph8_p1, ph8_p3, ph10_p1, ph10_p3, ph12_p1, ph12_p3, ab_2, pc_71 : simd::cache_line_size())
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

        const auto h2_p1 = ph2_p1[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h4_p3 = ph4_p3[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h8_p3 = ph8_p3[k];
        const auto h10_p1 = ph10_p1[k];
        const auto h10_p3 = ph10_p3[k];
        const auto h12_p1 = ph12_p1[k];
        const auto h12_p3 = ph12_p3[k];

        pc_71[k] = -fs_13395375_512 * e_0 * h2_p1 + e_1 * (f_2385_16 * h4_p1 + fs_694575_256 * h4_p3 + fs_6834375_128 * r_2 * h2_p1) + e_2 * (-fs_590625_242 * h6_p1 - fs_1063125_1936 * h6_p3 - f_2385_22 * r_2 * h4_p1 - fs_694575_484 * r_2 * h4_p3 - fs_84375_8 * r_4 * h2_p1) + e_3 * (fs_40186125_654368 * h8_p1 + fs_2083725_59488 * h8_p3 + fs_47250_121 * r_2 * h6_p1 + fs_42525_484 * r_2 * h6_p3 + f_7155_286 * r_4 * h4_p1 + fs_6251175_81796 * r_4 * h4_p3 + fs_84375_242 * r_6 * h2_p1) + e_4 * (-fs_3143448_17631601 * h10_p1 - fs_6417873_5425108 * h10_p3 - fs_40186125_14764178 * r_2 * h8_p1 - fs_2083725_1342198 * r_2 * h8_p3 - fs_189000_34969 * r_4 * h6_p1 - fs_42525_34969 * r_4 * h6_p3 - f_318_143 * r_6 * h4_p1 - fs_12348_20449 * r_6 * h4_p3 - fs_84375_40898 * r_8 * h2_p1) + e_5 * (-fs_17788815_717470533 * h12_p1 + fs_33960465_717470533 * h12_p3 + fs_12573792_9327116929 * r_2 * h10_p1 + fs_6417873_717470533 * r_2 * h10_p3 + fs_91125_14764178 * r_4 * h8_p1 + fs_4725_1342198 * r_4 * h8_p3 + fs_84000_12623809 * r_6 * h6_p1 + fs_18900_12623809 * r_6 * h6_p3 + f_159_2431 * r_8 * h4_p1 + fs_3087_5909761 * r_8 * h4_p3 + fs_30_20449 * r_10 * h2_p1);
    }

    // NOTE: the rows are formed in 82 loops, as the vectorizer runs out of
    // registers with all 91 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph2_p2, ph4_p2, ph4_p4, ph6_p2, ph6_p4, ph8_p2, ph8_p4, ph10_p2, ph10_p4, ph12_p2, ph12_p4, ab_2, pc_72 : simd::cache_line_size())
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

        const auto h2_p2 = ph2_p2[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h4_p4 = ph4_p4[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h6_p4 = ph6_p4[k];
        const auto h8_p2 = ph8_p2[k];
        const auto h8_p4 = ph8_p4[k];
        const auto h10_p2 = ph10_p2[k];
        const auto h10_p4 = ph10_p4[k];
        const auto h12_p2 = ph12_p2[k];
        const auto h12_p4 = ph12_p4[k];

        pc_72[k] = fs_13395375_128 * e_0 * h2_p2 + e_1 * (-fs_18225_8 * h4_p2 - fs_6251175_128 * h4_p4 - fs_6834375_32 * r_2 * h2_p2) + e_2 * (-fs_1063125_1936 * h6_p2 + fs_4921875_3872 * h6_p4 + fs_145800_121 * r_2 * h4_p2 + fs_6251175_242 * r_2 * h4_p4 + fs_84375_2 * r_4 * h2_p2) + e_3 * (fs_27860175_327184 * h8_p2 - fs_385875_59488 * h8_p4 + fs_42525_484 * r_2 * h6_p2 - fs_196875_968 * r_2 * h6_p4 - fs_1312200_20449 * r_4 * h4_p2 - fs_56260575_40898 * r_4 * h4_p4 - fs_168750_121 * r_6 * h2_p2) + e_4 * (-fs_130977_97682 * h10_p2 - fs_1178793_10850216 * h10_p4 - fs_77175_20449 * r_2 * h8_p2 + fs_385875_1342198 * r_2 * h8_p4 - fs_42525_34969 * r_4 * h6_p2 + fs_196875_69938 * r_4 * h6_p4 + fs_10368_20449 * r_6 * h4_p2 + fs_222264_20449 * r_6 * h4_p4 + fs_168750_20449 * r_8 * h2_p2) + e_5 * (-fs_10062360_717470533 * h12_p2 + fs_26832960_717470533 * h12_p4 + fs_261954_25836889 * r_2 * h10_p2 + fs_1178793_1434941066 * r_2 * h10_p4 + fs_175_20449 * r_4 * h8_p2 - fs_875_1342198 * r_4 * h8_p4 + fs_18900_12623809 * r_6 * h6_p2 - fs_43750_12623809 * r_6 * h6_p4 - fs_2592_5909761 * r_8 * h4_p2 - fs_55566_5909761 * r_8 * h4_p4 - fs_120_20449 * r_10 * h2_p2);
    }

    // NOTE: the rows are formed in 82 loops, as the vectorizer runs out of
    // registers with all 91 of them in one.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph4_p3, ph6_p3, ph6_p5, ph8_p3, ph8_p5, ph10_p3, ph10_p5, ph12_p3, ph12_p5, ab_2, pc_73 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h4_p3 = ph4_p3[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h6_p5 = ph6_p5[k];
        const auto h8_p3 = ph8_p3[k];
        const auto h8_p5 = ph8_p5[k];
        const auto h10_p3 = ph10_p3[k];
        const auto h10_p5 = ph10_p5[k];
        const auto h12_p3 = ph12_p3[k];
        const auto h12_p5 = ph12_p5[k];

        pc_73[k] = -fs_5315625_128 * e_1 * h4_p3 + e_2 * (fs_4921875_3872 * h6_p3 + fs_1063125_352 * h6_p5 + fs_5315625_242 * r_2 * h4_p3) + e_3 * (fs_385875_29744 * h8_p3 - fs_231525_2288 * h8_p5 - fs_196875_968 * r_2 * h6_p3 - fs_42525_88 * r_2 * h6_p5 - fs_47840625_40898 * r_4 * h4_p3) + e_4 * (-fs_26413695_10850216 * h10_p3 + fs_5282739_10850216 * h10_p5 - fs_385875_671099 * r_2 * h8_p3 + fs_231525_51623 * r_2 * h8_p5 + fs_196875_69938 * r_4 * h6_p3 + fs_42525_6358 * r_4 * h6_p5 + fs_189000_20449 * r_6 * h4_p3) + e_5 * (-fs_4528062_717470533 * h12_p3 + fs_1006236_42204149 * h12_p5 + fs_26413695_1434941066 * r_2 * h10_p3 - fs_5282739_1434941066 * r_2 * h10_p5 + fs_875_671099 * r_4 * h8_p3 - fs_525_51623 * r_4 * h8_p5 - fs_43750_12623809 * r_6 * h6_p3 - fs_9450_1147619 * r_6 * h6_p5 - fs_47250_5909761 * r_8 * h4_p3);
    }

    // NOTE: the rows are formed in 82 loops, as the vectorizer runs out of
    // registers with all 91 of them in one.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph4_p4, ph6_p4, ph6_p6, ph8_p4, ph8_p6, ph10_p4, ph10_p6, ph12_p4, ph12_p6, ab_2, pc_74 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h4_p4 = ph4_p4[k];
        const auto h6_p4 = ph6_p4[k];
        const auto h6_p6 = ph6_p6[k];
        const auto h8_p4 = ph8_p4[k];
        const auto h8_p6 = ph8_p6[k];
        const auto h10_p4 = ph10_p4[k];
        const auto h10_p6 = ph10_p6[k];
        const auto h12_p4 = ph12_p4[k];
        const auto h12_p6 = ph12_p6[k];

        pc_74[k] = fs_2338875_128 * e_1 * h4_p4 + e_2 * (fs_1063125_352 * h6_p4 - fs_39375_16 * h6_p6 - fs_212625_22 * r_2 * h4_p4) + e_3 * (-fs_231525_5408 * h8_p4 - fs_11025_208 * h8_p6 - fs_42525_88 * r_2 * h6_p4 + fs_1575_4 * r_2 * h6_p6 + fs_1913625_3718 * r_4 * h4_p4) + e_4 * (-fs_24310125_10850216 * h10_p4 + fs_257985_104329 * h10_p6 + fs_231525_122018 * r_2 * h8_p4 + fs_11025_4693 * r_2 * h8_p6 + fs_42525_6358 * r_4 * h6_p4 - fs_1575_289 * r_4 * h6_p6 - fs_7560_1859 * r_6 * h4_p4) + e_5 * (-fs_1463616_717470533 * h12_p4 + fs_470448_42204149 * h12_p6 + fs_24310125_1434941066 * r_2 * h10_p4 - fs_1031940_55190041 * r_2 * h10_p6 - fs_525_122018 * r_4 * h8_p4 - fs_25_4693 * r_4 * h8_p6 - fs_9450_1147619 * r_6 * h6_p4 + fs_700_104329 * r_6 * h6_p6 + fs_1890_537251 * r_8 * h4_p4);
    }

    // NOTE: the rows are formed in 82 loops, as the vectorizer runs out of
    // registers with all 91 of them in one.

#pragma omp simd aligned(pe_2, pe_3, pe_4, pe_5, ph6_p5, ph8_p5, ph8_p7, ph10_p5, ph10_p7, ph12_p5, ph12_p7, ab_2, pc_75 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;

        const auto h6_p5 = ph6_p5[k];
        const auto h8_p5 = ph8_p5[k];
        const auto h8_p7 = ph8_p7[k];
        const auto h10_p5 = ph10_p5[k];
        const auto h10_p7 = ph10_p7[k];
        const auto h12_p5 = ph12_p5[k];
        const auto h12_p7 = ph12_p7[k];

        pc_75[k] = -fs_39375_16 * e_2 * h6_p5 + e_3 * (-fs_77175_416 * h8_p5 + fs_55125_416 * h8_p7 + fs_1575_4 * r_2 * h6_p5) + e_4 * (-fs_5250987_5425108 * h10_p5 + fs_238140_79781 * h10_p7 + fs_77175_9386 * r_2 * h8_p5 - fs_55125_9386 * r_2 * h8_p7 - fs_1575_289 * r_4 * h6_p5) + e_5 * (-fs_15246_42204149 * h12_p5 + fs_6534_2221271 * h12_p7 + fs_5250987_717470533 * r_2 * h10_p5 - fs_952560_42204149 * r_2 * h10_p7 - fs_175_9386 * r_4 * h8_p5 + fs_125_9386 * r_4 * h8_p7 + fs_700_104329 * r_6 * h6_p5);
    }

    // NOTE: the rows are formed in 82 loops, as the vectorizer runs out of
    // registers with all 91 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph2_0, ph4_0, ph4_p4, ph6_0, ph6_p4, ph8_0, ph8_p4, ph10_0, ph10_p4, ph12_0, ph12_p4, ab_2, pc_76 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;
        const auto r_10 = r_8 * r_2;
        const auto r_12 = r_10 * r_2;

        const auto h2_0 = ph2_0[k];
        const auto h4_0 = ph4_0[k];
        const auto h4_p4 = ph4_p4[k];
        const auto h6_0 = ph6_0[k];
        const auto h6_p4 = ph6_p4[k];
        const auto h8_0 = ph8_0[k];
        const auto h8_p4 = ph8_p4[k];
        const auto h10_0 = ph10_0[k];
        const auto h10_p4 = ph10_p4[k];
        const auto h12_0 = ph12_0[k];
        const auto h12_p4 = ph12_p4[k];

        pc_76[k] = e_0 * (-f_4725_16 * h2_0 - f_10395_16 * r_2) + e_1 * (f_495_16 * h4_0 + fs_3472875_64 * h4_p4 + f_3375_8 * r_2 * h2_0 + f_10395_16 * r_4) + e_2 * (f_75_2 * h6_0 - fs_354375_121 * h6_p4 - f_45_2 * r_2 * h4_0 - fs_3472875_121 * r_2 * h4_p4 - f_375_2 * r_4 * h2_0 - f_495_2 * r_6) + e_3 * (-f_7455_572 * h8_0 + fs_694575_7436 * h8_p4 - f_15 * r_2 * h6_0 + fs_56700_121 * r_2 * h6_p4 + f_135_26 * r_4 * h4_0 + fs_31255875_20449 * r_4 * h4_p4 + f_375_11 * r_6 * h2_0 + f_165_4 * r_8) + e_4 * (f_6615_4199 * h10_0 - fs_2619540_1356277 * h10_p4 + f_7455_2717 * r_2 * h8_0 - fs_2778300_671099 * r_2 * h8_p4 + f_30_17 * r_4 * h6_0 - fs_226800_34969 * r_4 * h6_p4 - f_6_13 * r_6 * h4_0 - fs_246960_20449 * r_6 * h4_p4 - f_375_143 * r_8 * h2_0 - f_3 * r_10) + e_5 * (f_16335_96577 * h12_0 + fs_37733850_717470533 * h12_p4 - f_13230_96577 * r_2 * h10_0 + fs_10478160_717470533 * r_2 * h10_p4 - f_355_2717 * r_4 * h8_0 + fs_6300_671099 * r_4 * h8_p4 - f_20_323 * r_6 * h6_0 + fs_100800_12623809 * r_6 * h6_p4 + f_3_221 * r_8 * h4_0 + fs_61740_5909761 * r_8 * h4_p4 + f_10_143 * r_10 * h2_0 + f_1_13 * r_12) + f_10395_64 * e_6;
    }

    // NOTE: the rows are formed in 82 loops, as the vectorizer runs out of
    // registers with all 91 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph2_p1, ph4_p1, ph6_p1, ph6_p5, ph8_p1, ph8_p5, ph10_p1, ph10_p5, ph12_p1, ph12_p5, ab_2, pc_77 : simd::cache_line_size())
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

        const auto h2_p1 = ph2_p1[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h6_p5 = ph6_p5[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h8_p5 = ph8_p5[k];
        const auto h10_p1 = ph10_p1[k];
        const auto h10_p5 = ph10_p5[k];
        const auto h12_p1 = ph12_p1[k];
        const auto h12_p5 = ph12_p5[k];

        pc_77[k] = -fs_66976875_1024 * e_0 * h2_p1 + e_1 * (fs_15400125_512 * h4_p1 + fs_34171875_256 * r_2 * h2_p1) + e_2 * (-fs_1063125_1936 * h6_p1 - fs_39375_88 * h6_p5 - fs_15400125_968 * r_2 * h4_p1 - fs_421875_16 * r_4 * h2_p1) + e_3 * (-f_105_22 * h8_p1 + fs_77175_2288 * h8_p5 + fs_42525_484 * r_2 * h6_p1 + fs_1575_22 * r_2 * h6_p5 + fs_820125_968 * r_4 * h4_p1 + fs_421875_484 * r_6 * h2_p1) + e_4 * (fs_159137055_70526404 * h10_p1 - fs_3274425_2712554 * h10_p5 + f_210_209 * r_2 * h8_p1 - fs_77175_51623 * r_2 * h8_p5 - fs_42525_34969 * r_4 * h6_p1 - fs_3150_3179 * r_4 * h6_p5 - fs_810_121 * r_6 * h4_p1 - fs_421875_81796 * r_8 * h2_p1) + e_5 * (fs_9882675_1434941066 * h12_p1 + fs_2096325_42204149 * h12_p5 - fs_159137055_9327116929 * r_2 * h10_p1 + fs_6548850_717470533 * r_2 * h10_p5 - f_10_209 * r_4 * h8_p1 + fs_175_51623 * r_4 * h8_p5 + fs_18900_12623809 * r_6 * h6_p1 + fs_1400_1147619 * r_6 * h6_p5 + fs_405_69938 * r_8 * h4_p1 + fs_75_20449 * r_10 * h2_p1);
    }

    // NOTE: the rows are formed in 82 loops, as the vectorizer runs out of
    // registers with all 91 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph2_p2, ph4_p2, ph6_p2, ph6_p6, ph8_p2, ph8_p6, ph10_p2, ph10_p6, ph12_p2, ph12_p6, ab_2, pc_78 : simd::cache_line_size())
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

        const auto h2_p2 = ph2_p2[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h6_p6 = ph6_p6[k];
        const auto h8_p2 = ph8_p2[k];
        const auto h8_p6 = ph8_p6[k];
        const auto h10_p2 = ph10_p2[k];
        const auto h10_p6 = ph10_p6[k];
        const auto h12_p2 = ph12_p2[k];
        const auto h12_p6 = ph12_p6[k];

        pc_78[k] = fs_40186125_512 * e_0 * h2_p2 + e_1 * (fs_3213675_512 * h4_p2 - fs_20503125_128 * r_2 * h2_p2) + e_2 * (-fs_354375_121 * h6_p2 + fs_196875_44 * h6_p6 - fs_3213675_968 * r_2 * h4_p2 + fs_253125_8 * r_4 * h2_p2) + e_3 * (fs_231525_20449 * h8_p2 - fs_55125_2288 * h8_p6 + fs_56700_121 * r_2 * h6_p2 - fs_7875_11 * r_2 * h6_p6 + fs_28923075_163592 * r_4 * h4_p2 - fs_253125_242 * r_6 * h2_p2) + e_4 * (fs_80725491_35263202 * h10_p2 - fs_43659_1356277 * h10_p6 - fs_3704400_7382089 * r_2 * h8_p2 + fs_55125_51623 * r_2 * h8_p6 - fs_226800_34969 * r_4 * h6_p2 + fs_31500_3179 * r_4 * h6_p6 - fs_28566_20449 * r_6 * h4_p2 + fs_253125_40898 * r_8 * h2_p2) + e_5 * (fs_3773385_1434941066 * h12_p2 + fs_1617165_42204149 * h12_p6 - fs_161450982_9327116929 * r_2 * h10_p2 + fs_174636_717470533 * r_2 * h10_p6 + fs_8400_7382089 * r_4 * h8_p2 - fs_125_51623 * r_4 * h8_p6 + fs_100800_12623809 * r_6 * h6_p2 - fs_14000_1147619 * r_6 * h6_p6 + fs_14283_11819522 * r_8 * h4_p2 - fs_90_20449 * r_10 * h2_p2);
    }

    // NOTE: the rows are formed in 82 loops, as the vectorizer runs out of
    // registers with all 91 of them in one.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph4_p3, ph6_p3, ph8_p3, ph8_p7, ph10_p3, ph10_p7, ph12_p3, ph12_p7, ab_2, pc_79 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h4_p3 = ph4_p3[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h8_p3 = ph8_p3[k];
        const auto h8_p7 = ph8_p7[k];
        const auto h10_p3 = ph10_p3[k];
        const auto h10_p7 = ph10_p7[k];
        const auto h12_p3 = ph12_p3[k];
        const auto h12_p7 = ph12_p7[k];

        pc_79[k] = -fs_22920975_512 * e_1 * h4_p3 + e_2 * (-fs_39375_88 * h6_p3 + fs_2083725_88 * r_2 * h4_p3) + e_3 * (fs_77175_676 * h8_p3 - fs_33075_208 * h8_p7 + fs_1575_22 * r_2 * h6_p3 - fs_18753525_14872 * r_4 * h4_p3) + e_4 * (fs_3814209_2712554 * h10_p3 + fs_194481_159562 * h10_p7 - fs_308700_61009 * r_2 * h8_p3 + fs_33075_4693 * r_2 * h8_p7 - fs_3150_3179 * r_4 * h6_p3 + fs_18522_1859 * r_6 * h4_p3) + e_5 * (fs_1029105_1434941066 * h12_p3 + fs_49005_2221271 * h12_p7 - fs_7628418_717470533 * r_2 * h10_p3 - fs_388962_42204149 * r_2 * h10_p7 + fs_700_61009 * r_4 * h8_p3 - fs_75_4693 * r_4 * h8_p7 + fs_1400_1147619 * r_6 * h6_p3 - fs_9261_1074502 * r_8 * h4_p3);
    }

    // NOTE: the rows are formed in 82 loops, as the vectorizer runs out of
    // registers with all 91 of them in one.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph4_p4, ph6_p4, ph8_p4, ph8_p8, ph10_p4, ph10_p8, ph12_p4, ph12_p8, ab_2, pc_80 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h4_p4 = ph4_p4[k];
        const auto h6_p4 = ph6_p4[k];
        const auto h8_p4 = ph8_p4[k];
        const auto h8_p8 = ph8_p8[k];
        const auto h10_p4 = ph10_p4[k];
        const auto h10_p8 = ph10_p8[k];
        const auto h12_p4 = ph12_p4[k];
        const auto h12_p8 = ph12_p8[k];

        pc_80[k] = fs_1403325_256 * e_1 * h4_p4 + e_2 * (fs_196875_44 * h6_p4 - fs_127575_44 * r_2 * h4_p4) + e_3 * (fs_385875_2704 * h8_p4 + fs_11025_208 * h8_p8 - fs_7875_11 * r_2 * h6_p4 + fs_1148175_7436 * r_4 * h4_p4) + e_4 * (fs_583443_1356277 * h10_p4 + fs_321489_79781 * h10_p8 - fs_385875_61009 * r_2 * h8_p4 - fs_11025_4693 * r_2 * h8_p8 + fs_31500_3179 * r_4 * h6_p4 - fs_2268_1859 * r_6 * h4_p4) + e_5 * (fs_76230_717470533 * h12_p4 + fs_16335_2221271 * h12_p8 - fs_2333772_717470533 * r_2 * h10_p4 - fs_1285956_42204149 * r_2 * h10_p8 + fs_875_61009 * r_4 * h8_p4 + fs_25_4693 * r_4 * h8_p8 - fs_14000_1147619 * r_6 * h6_p4 + fs_567_537251 * r_8 * h4_p4);
    }

    // NOTE: the rows are formed in 82 loops, as the vectorizer runs out of
    // registers with all 91 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph2_0, ph4_0, ph6_0, ph6_p6, ph8_0, ph8_p6, ph10_0, ph10_p6, ph12_0, ph12_p6, ab_2, pc_81 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;
        const auto r_10 = r_8 * r_2;
        const auto r_12 = r_10 * r_2;

        const auto h2_0 = ph2_0[k];
        const auto h4_0 = ph4_0[k];
        const auto h6_0 = ph6_0[k];
        const auto h6_p6 = ph6_p6[k];
        const auto h8_0 = ph8_0[k];
        const auto h8_p6 = ph8_p6[k];
        const auto h10_0 = ph10_0[k];
        const auto h10_p6 = ph10_p6[k];
        const auto h12_0 = ph12_0[k];
        const auto h12_p6 = ph12_p6[k];

        pc_81[k] = e_0 * (-f_4725_32 * h2_0 - f_10395_16 * r_2) + e_1 * (-f_1215_8 * h4_0 + f_3375_16 * r_2 * h2_0 + f_10395_16 * r_4) + e_2 * (f_3225_44 * h6_0 - fs_118125_22 * h6_p6 + f_1215_11 * r_2 * h4_0 - f_375_4 * r_4 * h2_0 - f_495_2 * r_6) + e_3 * (-f_1995_572 * h8_0 + fs_33075_286 * h8_p6 - f_645_22 * r_2 * h6_0 + fs_9450_11 * r_2 * h6_p6 - f_3645_143 * r_4 * h4_0 + f_375_22 * r_6 * h2_0 + f_165_4 * r_8) + e_4 * (-f_945_442 * h10_0 - fs_5893965_2712554 * h10_p6 + f_105_143 * r_2 * h8_0 - fs_264600_51623 * r_2 * h8_p6 + f_645_187 * r_4 * h6_0 - fs_37800_3179 * r_4 * h6_p6 + f_324_143 * r_6 * h4_0 - f_375_286 * r_8 * h2_0 - f_3 * r_10) + e_5 * (-f_7260_96577 * h12_0 + fs_2395800_42204149 * h12_p6 + f_945_5083 * r_2 * h10_0 + fs_11787930_717470533 * r_2 * h10_p6 - f_5_143 * r_4 * h8_0 + fs_600_51623 * r_4 * h8_p6 - f_430_3553 * r_6 * h6_0 + fs_16800_1147619 * r_6 * h6_p6 - f_162_2431 * r_8 * h4_0 + f_5_143 * r_10 * h2_0 + f_1_13 * r_12) + f_10395_64 * e_6;
    }

    // NOTE: the rows are formed in 82 loops, as the vectorizer runs out of
    // registers with all 91 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph2_p1, ph4_p1, ph6_p1, ph8_p1, ph8_p7, ph10_p1, ph10_p7, ph12_p1, ph12_p7, ab_2, pc_82 : simd::cache_line_size())
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

        const auto h2_p1 = ph2_p1[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h8_p7 = ph8_p7[k];
        const auto h10_p1 = ph10_p1[k];
        const auto h10_p7 = ph10_p7[k];
        const auto h12_p1 = ph12_p1[k];
        const auto h12_p7 = ph12_p7[k];

        pc_82[k] = -fs_218791125_2048 * e_0 * h2_p1 + e_1 * (fs_2679075_256 * h4_p1 + fs_111628125_512 * r_2 * h2_p1) + e_2 * (fs_4921875_3872 * h6_p1 - fs_2679075_484 * r_2 * h4_p1 - fs_1378125_32 * r_4 * h2_p1) + e_3 * (-fs_13395375_163592 * h8_p1 + fs_33075_1144 * h8_p7 - fs_196875_968 * r_2 * h6_p1 + fs_24111675_81796 * r_4 * h4_p1 + fs_1378125_968 * r_6 * h2_p1) + e_4 * (-fs_220172337_141052808 * h10_p1 - fs_392931_319124 * h10_p7 + fs_26790750_7382089 * r_2 * h8_p1 - fs_66150_51623 * r_2 * h8_p7 + fs_196875_69938 * r_4 * h6_p1 - fs_47628_20449 * r_6 * h4_p1 - fs_1378125_163592 * r_8 * h2_p1) + e_5 * (-fs_658845_717470533 * h12_p1 + fs_119790_2221271 * h12_p7 + fs_220172337_18654233858 * r_2 * h10_p1 + fs_392931_42204149 * r_2 * h10_p7 - fs_60750_7382089 * r_4 * h8_p1 + fs_150_51623 * r_4 * h8_p7 - fs_43750_12623809 * r_6 * h6_p1 + fs_11907_5909761 * r_8 * h4_p1 + fs_245_40898 * r_10 * h2_p1);
    }

    // NOTE: the rows are formed in 82 loops, as the vectorizer runs out of
    // registers with all 91 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph2_p2, ph4_p2, ph6_p2, ph8_p2, ph8_p8, ph10_p2, ph10_p8, ph12_p2, ph12_p8, ab_2, pc_83 : simd::cache_line_size())
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

        const auto h2_p2 = ph2_p2[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h8_p2 = ph8_p2[k];
        const auto h8_p8 = ph8_p8[k];
        const auto h10_p2 = ph10_p2[k];
        const auto h10_p8 = ph10_p8[k];
        const auto h12_p2 = ph12_p2[k];
        const auto h12_p8 = ph12_p8[k];

        pc_83[k] = fs_49116375_1024 * e_0 * h2_p2 + e_1 * (fs_601425_16 * h4_p2 - fs_25059375_256 * r_2 * h2_p2) + e_2 * (-fs_39375_88 * h6_p2 - fs_218700_11 * r_2 * h4_p2 + fs_309375_16 * r_4 * h2_p2) + e_3 * (-fs_2083725_14872 * h8_p2 - fs_33075_208 * h8_p8 + fs_1575_22 * r_2 * h6_p2 + fs_1968300_1859 * r_4 * h4_p2 - fs_28125_44 * r_6 * h2_p2) + e_4 * (-f_6993_8398 * h10_p2 + fs_11907_319124 * h10_p8 + fs_4167450_671099 * r_2 * h8_p2 + fs_33075_4693 * r_2 * h8_p8 - fs_3150_3179 * r_4 * h6_p2 - fs_15552_1859 * r_6 * h4_p2 + fs_28125_7436 * r_8 * h2_p2) + e_5 * (-fs_152460_717470533 * h12_p2 + fs_87120_2221271 * h12_p8 + f_6993_96577 * r_2 * h10_p2 - fs_11907_42204149 * r_2 * h10_p8 - fs_9450_671099 * r_4 * h8_p2 - fs_75_4693 * r_4 * h8_p8 + fs_1400_1147619 * r_6 * h6_p2 + fs_3888_537251 * r_8 * h4_p2 - fs_5_1859 * r_10 * h2_p2);
    }

    // NOTE: the rows are formed in 82 loops, as the vectorizer runs out of
    // registers with all 91 of them in one.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph4_p3, ph6_p3, ph8_p3, ph10_p3, ph10_p9, ph12_p3, ph12_p9, ab_2, pc_84 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;

        const auto h4_p3 = ph4_p3[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h8_p3 = ph8_p3[k];
        const auto h10_p3 = ph10_p3[k];
        const auto h10_p9 = ph10_p9[k];
        const auto h12_p3 = ph12_p3[k];
        const auto h12_p9 = ph12_p9[k];

        pc_84[k] = -fs_12629925_512 * e_1 * h4_p3 + e_2 * (-fs_118125_22 * h6_p3 + fs_1148175_88 * r_2 * h4_p3) + e_3 * (-fs_231525_2704 * h8_p3 + fs_9450_11 * r_2 * h6_p3 - fs_10333575_14872 * r_4 * h4_p3) + e_4 * (-fs_214326_1356277 * h10_p3 + fs_35721_8398 * h10_p9 + fs_231525_61009 * r_2 * h8_p3 - fs_37800_3179 * r_4 * h6_p3 + fs_10206_1859 * r_6 * h4_p3) + e_5 * (-fs_38115_1434941066 * h12_p3 + fs_38115_2221271 * h12_p9 + fs_857304_717470533 * r_2 * h10_p3 - fs_71442_2221271 * r_2 * h10_p9 - fs_525_61009 * r_4 * h8_p3 + fs_16800_1147619 * r_6 * h6_p3 - fs_5103_1074502 * r_8 * h4_p3);
    }

    // NOTE: the rows are formed in 82 loops, as the vectorizer runs out of
    // registers with all 91 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph2_0, ph4_0, ph6_0, ph8_0, ph8_p8, ph10_0, ph10_p8, ph12_0, ph12_p8, ab_2, pc_85 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;
        const auto r_10 = r_8 * r_2;
        const auto r_12 = r_10 * r_2;

        const auto h2_0 = ph2_0[k];
        const auto h4_0 = ph4_0[k];
        const auto h6_0 = ph6_0[k];
        const auto h8_0 = ph8_0[k];
        const auto h8_p8 = ph8_p8[k];
        const auto h10_0 = ph10_0[k];
        const auto h10_p8 = ph10_p8[k];
        const auto h12_0 = ph12_0[k];
        const auto h12_p8 = ph12_p8[k];

        pc_85[k] = e_0 * (f_945_16 * h2_0 - f_10395_16 * r_2) + e_1 * (-f_270 * h4_0 - f_675_8 * r_2 * h2_0 + f_10395_16 * r_4) + e_2 * (f_150_11 * h6_0 + f_2160_11 * r_2 * h4_0 + f_75_2 * r_4 * h2_0 - f_495_2 * r_6) + e_3 * (f_9345_572 * h8_0 + fs_496125_2288 * h8_p8 - f_60_11 * r_2 * h6_0 - f_6480_143 * r_4 * h4_0 - f_75_11 * r_6 * h2_0 + f_165_4 * r_8) + e_4 * (f_5229_4199 * h10_0 - fs_218295_79781 * h10_p8 - f_9345_2717 * r_2 * h8_0 - fs_496125_51623 * r_2 * h8_p8 + f_120_187 * r_4 * h6_0 + f_576_143 * r_6 * h4_0 + f_75_143 * r_8 * h2_0 - f_3 * r_10) + e_5 * (f_2178_96577 * h12_0 + fs_143748_2221271 * h12_p8 - f_10458_96577 * r_2 * h10_0 + fs_873180_42204149 * r_2 * h10_p8 + f_445_2717 * r_4 * h8_0 + fs_1125_51623 * r_4 * h8_p8 - f_80_3553 * r_6 * h6_0 - f_288_2431 * r_8 * h4_0 - f_2_143 * r_10 * h2_0 + f_1_13 * r_12) + f_10395_64 * e_6;
    }

    // NOTE: the rows are formed in 82 loops, as the vectorizer runs out of
    // registers with all 91 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph2_p1, ph4_p1, ph6_p1, ph8_p1, ph10_p1, ph10_p9, ph12_p1, ph12_p9, ab_2, pc_86 : simd::cache_line_size())
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

        const auto h2_p1 = ph2_p1[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h10_p1 = ph10_p1[k];
        const auto h10_p9 = ph10_p9[k];
        const auto h12_p1 = ph12_p1[k];
        const auto h12_p9 = ph12_p9[k];

        pc_86[k] = -fs_265228425_2048 * e_0 * h2_p1 + e_1 * (-fs_1002375_256 * h4_p1 + fs_135320625_512 * r_2 * h2_p1) + e_2 * (fs_1063125_352 * h6_p1 + fs_91125_44 * r_2 * h4_p1 - fs_1670625_32 * r_4 * h2_p1) + e_3 * (fs_198450_1859 * h8_p1 - fs_42525_88 * r_2 * h6_p1 - fs_820125_7436 * r_4 * h4_p1 + fs_151875_88 * r_6 * h2_p1) + e_4 * (fs_36693405_141052808 * h10_p1 - fs_19845_16796 * h10_p9 - fs_3175200_671099 * r_2 * h8_p1 + fs_42525_6358 * r_4 * h6_p1 + fs_1620_1859 * r_6 * h4_p1 - fs_151875_14872 * r_8 * h2_p1) + e_5 * (fs_35937_717470533 * h12_p1 + fs_137214_2221271 * h12_p9 - fs_36693405_18654233858 * r_2 * h10_p1 + fs_19845_2221271 * r_2 * h10_p9 + fs_7200_671099 * r_4 * h8_p1 - fs_9450_1147619 * r_6 * h6_p1 - fs_405_537251 * r_8 * h4_p1 + fs_27_3718 * r_10 * h2_p1);
    }

    // NOTE: the rows are formed in 82 loops, as the vectorizer runs out of
    // registers with all 91 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph2_p2, ph4_p2, ph6_p2, ph8_p2, ph10_p2, ph10_p10, ph12_p2, ph12_p10, ab_2, pc_87 : simd::cache_line_size())
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

        const auto h2_p2 = ph2_p2[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h8_p2 = ph8_p2[k];
        const auto h10_p2 = ph10_p2[k];
        const auto h10_p10 = ph10_p10[k];
        const auto h12_p2 = ph12_p2[k];
        const auto h12_p10 = ph12_p10[k];

        pc_87[k] = fs_9823275_512 * e_0 * h2_p2 + e_1 * (fs_27064125_512 * h4_p2 - fs_5011875_128 * r_2 * h2_p2) + e_2 * (fs_196875_44 * h6_p2 - fs_2460375_88 * r_2 * h4_p2 + fs_61875_8 * r_4 * h2_p2) + e_3 * (fs_1157625_29744 * h8_p2 - fs_7875_11 * r_2 * h6_p2 + fs_22143375_14872 * r_4 * h4_p2 - fs_5625_22 * r_6 * h2_p2) + e_4 * (fs_1607445_35263202 * h10_p2 + fs_11907_4199 * h10_p10 - fs_1157625_671099 * r_2 * h8_p2 + fs_31500_3179 * r_4 * h6_p2 - fs_21870_1859 * r_6 * h4_p2 + fs_5625_3718 * r_8 * h2_p2) + e_5 * (fs_7623_1434941066 * h12_p2 + fs_83853_2221271 * h12_p10 - fs_3214890_9327116929 * r_2 * h10_p2 - fs_47628_2221271 * r_2 * h10_p10 + fs_2625_671099 * r_4 * h8_p2 - fs_14000_1147619 * r_6 * h6_p2 + fs_10935_1074502 * r_8 * h4_p2 - fs_2_1859 * r_10 * h2_p2);
    }

    // NOTE: the rows are formed in 82 loops, as the vectorizer runs out of
    // registers with all 91 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph2_0, ph4_0, ph6_0, ph8_0, ph10_0, ph10_p10, ph12_0, ph12_p10, ab_2, pc_88 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;
        const auto r_10 = r_8 * r_2;
        const auto r_12 = r_10 * r_2;

        const auto h2_0 = ph2_0[k];
        const auto h4_0 = ph4_0[k];
        const auto h6_0 = ph6_0[k];
        const auto h8_0 = ph8_0[k];
        const auto h10_0 = ph10_0[k];
        const auto h10_p10 = ph10_p10[k];
        const auto h12_0 = ph12_0[k];
        const auto h12_p10 = ph12_p10[k];

        pc_88[k] = e_0 * (f_10395_32 * h2_0 - f_10395_16 * r_2) + e_1 * (-f_1485_8 * h4_0 - f_7425_16 * r_2 * h2_0 + f_10395_16 * r_4) + e_2 * (-f_375_4 * h6_0 + f_135 * r_2 * h4_0 + f_825_4 * r_4 * h2_0 - f_495_2 * r_6) + e_3 * (-f_525_52 * h8_0 + f_75_2 * r_2 * h6_0 - f_405_13 * r_4 * h4_0 - f_75_2 * r_6 * h2_0 + f_165_4 * r_8) + e_4 * (-f_3087_8398 * h10_0 - fs_43659_8398 * h10_p10 + f_525_247 * r_2 * h8_0 - f_75_17 * r_4 * h6_0 + f_36_13 * r_6 * h4_0 + f_75_26 * r_8 * h2_0 - f_3 * r_10) + e_5 * (-f_396_96577 * h12_0 + fs_182952_2221271 * h12_p10 + f_3087_96577 * r_2 * h10_0 + fs_87318_2221271 * r_2 * h10_p10 - f_25_247 * r_4 * h8_0 + f_50_323 * r_6 * h6_0 - f_18_221 * r_8 * h4_0 - f_1_13 * r_10 * h2_0 + f_1_13 * r_12) + f_10395_64 * e_6;
    }

    // NOTE: the rows are formed in 82 loops, as the vectorizer runs out of
    // registers with all 91 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph2_p1, ph4_p1, ph6_p1, ph8_p1, ph10_p1, ph12_p1, ph12_p11, ab_2, pc_89 : simd::cache_line_size())
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

        const auto h2_p1 = ph2_p1[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h10_p1 = ph10_p1[k];
        const auto h12_p1 = ph12_p1[k];
        const auto h12_p11 = ph12_p11[k];

        pc_89[k] = -f_10395_32 * e_0 * h2_p1 + e_1 * (-fs_33078375_512 * h4_p1 + f_7425_16 * r_2 * h2_p1) + e_2 * (-fs_39375_16 * h6_p1 + fs_273375_8 * r_2 * h4_p1 - f_825_4 * r_4 * h2_p1) + e_3 * (-fs_33075_2704 * h8_p1 + fs_1575_4 * r_2 * h6_p1 - fs_2460375_1352 * r_4 * h4_p1 + f_75_2 * r_6 * h2_p1) + e_4 * (-fs_654885_70526404 * h10_p1 + fs_33075_61009 * r_2 * h8_p1 - fs_1575_289 * r_4 * h6_p1 + fs_2430_169 * r_6 * h4_p1 - f_75_26 * r_8 * h2_p1) + e_5 * (-fs_1089_1434941066 * h12_p1 + fs_7623_96577 * h12_p11 + fs_654885_9327116929 * r_2 * h10_p1 - fs_75_61009 * r_4 * h8_p1 + fs_700_104329 * r_6 * h6_p1 - fs_1215_97682 * r_8 * h4_p1 + f_1_13 * r_10 * h2_p1);
    }

    // NOTE: the rows are formed in 82 loops, as the vectorizer runs out of
    // registers with all 91 of them in one.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph2_0, ph4_0, ph6_0, ph8_0, ph10_0, ph12_0, ph12_p12, ab_2, pc_90 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_4 * r_2;
        const auto r_8 = r_6 * r_2;
        const auto r_10 = r_8 * r_2;
        const auto r_12 = r_10 * r_2;

        const auto h2_0 = ph2_0[k];
        const auto h4_0 = ph4_0[k];
        const auto h6_0 = ph6_0[k];
        const auto h8_0 = ph8_0[k];
        const auto h10_0 = ph10_0[k];
        const auto h12_0 = ph12_0[k];
        const auto h12_p12 = ph12_p12[k];

        pc_90[k] = e_0 * (f_10395_16 * h2_0 - f_10395_16 * r_2) + e_1 * (f_4455_16 * h4_0 - f_7425_8 * r_2 * h2_0 + f_10395_16 * r_4) + e_2 * (f_75_2 * h6_0 - f_405_2 * r_2 * h4_0 + f_825_2 * r_4 * h2_0 - f_495_2 * r_6) + e_3 * (f_105_52 * h8_0 - f_15 * r_2 * h6_0 + f_1215_26 * r_4 * h4_0 - f_75 * r_6 * h2_0 + f_165_4 * r_8) + e_4 * (f_189_4199 * h10_0 - f_105_247 * r_2 * h8_0 + f_30_17 * r_4 * h6_0 - f_54_13 * r_6 * h4_0 + f_75_13 * r_8 * h2_0 - f_3 * r_10) + e_5 * (f_33_96577 * h12_0 + fs_15246_96577 * h12_p12 - f_378_96577 * r_2 * h10_0 + f_5_247 * r_4 * h8_0 - f_20_323 * r_6 * h6_0 + f_27_221 * r_8 * h4_0 - f_2_13 * r_10 * h2_0 + f_1_13 * r_12) + f_10395_64 * e_6;
    }

    // NOTE: the values of a combination of angular components are stored as one
    // row of nvalues columns, with the component on bra side running slowest, and
    // the atom pairs beyond the reach of every pair of primitives are set to zero.

    const size_t sources[169] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 1, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 2, 14, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 3, 15, 26, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 4, 16, 27, 37, 46, 47, 48, 49, 50, 51, 52, 53, 54, 5, 17, 28, 38, 47, 55, 56, 57, 58, 59, 60, 61, 62, 6, 18, 29, 39, 48, 56, 63, 64, 65, 66, 67, 68, 69, 7, 19, 30, 40, 49, 57, 64, 70, 71, 72, 73, 74, 75, 8, 20, 31, 41, 50, 58, 65, 71, 76, 77, 78, 79, 80, 9, 21, 32, 42, 51, 59, 66, 72, 77, 81, 82, 83, 84, 10, 22, 33, 43, 52, 60, 67, 73, 78, 82, 85, 86, 87, 11, 23, 34, 44, 53, 61, 68, 74, 79, 83, 86, 88, 89, 12, 24, 35, 45, 54, 62, 69, 75, 80, 84, 87, 89, 90};

    for (size_t m = 0; m < 169; m++)
    {
        const auto *pc = buffer.data(7 + sources[m]);

        std::copy(pc, pc + nmax, values + m * nvalues);

        std::fill(values + m * nvalues + nmax, values + (m + 1) * nvalues, 0.0);
    }
}

}  // namespace simdovl
