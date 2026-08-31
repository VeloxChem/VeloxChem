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



#include "SimdKineticEnergyRecIG.hpp"

#include <algorithm>
#include <ranges>
#include <cmath>
#include <string>

#include "ErrorHandler.hpp"
#include "MathConst.hpp"
#include "ScreeningFunc.hpp"
#include "SimdAlign.hpp"
#include "SimdDimensions.hpp"

namespace simdkin {  // simdkin namespace

auto
compute_ig_kinetic_energy(double                         *values,
                           const size_t                    nvalues,
                           const CBasisFunction           &bra,
                           const CBasisFunction           &ket,
                           const std::vector<CSimdMatrix> &harmonics,
                           const CSimdMatrix              &coordinates,
                           const double                    threshold) -> void
{
    if ((bra.get_angular_momentum() != 6) || (ket.get_angular_momentum() != 4))
    {
        errors::assertMsgCritical(
            false, std::string("SimdKineticEnergyFunc.compute_ig_kinetic_energy: Basis functions must be of angular momenta six and four"));
    }

    if (harmonics.size() < 10)
    {
        errors::assertMsgCritical(
            false, std::string("SimdKineticEnergyFunc.compute_ig_kinetic_energy: Harmonics must reach angular momentum 10"));
    }

    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("SimdKineticEnergyFunc.compute_ig_kinetic_energy: Number of values exceeds number of atom pairs"));
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
        bra, ket, nvalues, coordinates, screenfunc::two_center_kinetic_energy_primitive_bound,
        threshold / static_cast<double>(nprims));

    // NOTE: the buffer spans the atom pairs reached by the pair of primitives
    // reaching furthest, which is searched for rather than assumed. The
    // primitives are sorted by descending exponent, but the bound of a pair of
    // primitives carries their prefactor as well as their decay, so a tighter
    // pair with a larger prefactor reaches further than a more diffuse pair with
    // a smaller one, and the last pair is not always the furthest reaching.

    const auto nmax = *std::ranges::max_element(dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 117 * nvalues, 0.0);

        return;
    }

    // NOTE: the buffer holds the contracted prefactor of each exponent factor of
    // the terms, as the integrals of the angular components are formed straight
    // into the values and are not written a second time.

    auto buffer = CSimdMatrix(6, nmax);

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

    // accumulate the prefactor of each exponent factor over the pairs of primitives

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

            const auto ff_0 = fbase * bexp * bexp * fmu / (fexp * fexp * fexp * fexp * fexp * fexp);

            const auto ff_1 = fbase * bexp * bexp * fmu * fmu / (fexp * fexp * fexp * fexp * fexp * fexp);

            const auto ff_2 = fbase * bexp * bexp * fmu * fmu * fmu / (fexp * fexp * fexp * fexp * fexp * fexp);

            const auto ff_3 = fbase * bexp * bexp * fmu * fmu * fmu * fmu / (fexp * fexp * fexp * fexp * fexp * fexp);

            const auto ff_4 = fbase * bexp * bexp * fmu * fmu * fmu * fmu * fmu / (fexp * fexp * fexp * fexp * fexp * fexp);

            const auto ff_5 = fbase * bexp * bexp * fmu * fmu * fmu * fmu * fmu * fmu / (fexp * fexp * fexp * fexp * fexp * fexp);

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ab_2 : simd::cache_line_size())
            for (size_t k = 0; k < ncols; k++)
            {
                const auto fterm = std::exp(-fmu * ab_2[k]);

                pe_0[k] += ff_0 * fterm;
                pe_1[k] += ff_1 * fterm;
                pe_2[k] += ff_2 * fterm;
                pe_3[k] += ff_3 * fterm;
                pe_4[k] += ff_4 * fterm;
                pe_5[k] += ff_5 * fterm;
            }
        }
    }

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

    // NOTE: the factors of the terms depend on the angular momenta alone, so they
    // are formed once for the whole matrix instead of once for every atom pair.

    const auto f_1007_187 = 1007.0 / 187.0;
    const auto f_1007_22 = 1007.0 / 22.0;
    const auto f_1008_2717 = 1008.0 / 2717.0;
    const auto f_1020_11 = 1020.0 / 11.0;
    const auto f_1026_187 = 1026.0 / 187.0;
    const auto f_106_561 = 106.0 / 561.0;
    const auto f_1125_2 = 1125.0 / 2.0;
    const auto f_114_17 = 114.0 / 17.0;
    const auto f_1156_143 = 1156.0 / 143.0;
    const auto f_1216_187 = 1216.0 / 187.0;
    const auto f_1275_2 = 1275.0 / 2.0;
    const auto f_128_561 = 128.0 / 561.0;
    const auto f_13005_143 = 13005.0 / 143.0;
    const auto f_1360_143 = 1360.0 / 143.0;
    const auto f_136_13 = 136.0 / 13.0;
    const auto f_1520_11 = 1520.0 / 11.0;
    const auto f_15300_143 = 15300.0 / 143.0;
    const auto f_1530_13 = 1530.0 / 13.0;
    const auto f_1575_4 = 1575.0 / 4.0;
    const auto f_1596_187 = 1596.0 / 187.0;
    const auto f_16875_8 = 16875.0 / 8.0;
    const auto f_171_17 = 171.0 / 17.0;
    const auto f_171_2 = 171.0 / 2.0;
    const auto f_1824_187 = 1824.0 / 187.0;
    const auto f_1875_11 = 1875.0 / 11.0;
    const auto f_1875_143 = 1875.0 / 143.0;
    const auto f_1875_2 = 1875.0 / 2.0;
    const auto f_1995_11 = 1995.0 / 11.0;
    const auto f_19_17 = 19.0 / 17.0;
    const auto f_19_2 = 19.0 / 2.0;
    const auto f_204_143 = 204.0 / 143.0;
    const auto f_21168_2717 = 21168.0 / 2717.0;
    const auto f_2280_11 = 2280.0 / 11.0;
    const auto f_2295_143 = 2295.0 / 143.0;
    const auto f_23625_16 = 23625.0 / 16.0;
    const auto f_250_1 = 250.0;
    const auto f_255_11 = 255.0 / 11.0;
    const auto f_255_2 = 255.0 / 2.0;
    const auto f_255_8 = 255.0 / 8.0;
    const auto f_2565_22 = 2565.0 / 22.0;
    const auto f_272_143 = 272.0 / 143.0;
    const auto f_2805_4 = 2805.0 / 4.0;
    const auto f_285_2 = 285.0 / 2.0;
    const auto f_2940_4199 = 2940.0 / 4199.0;
    const auto f_2_143 = 2.0 / 143.0;
    const auto f_2_51 = 2.0 / 51.0;
    const auto f_3060_143 = 3060.0 / 143.0;
    const auto f_33810_4199 = 33810.0 / 4199.0;
    const auto f_34_143 = 34.0 / 143.0;
    const auto f_36_187 = 36.0 / 187.0;
    const auto f_40_143 = 40.0 / 143.0;
    const auto f_40_429 = 40.0 / 429.0;
    const auto f_4335_11 = 4335.0 / 11.0;
    const auto f_4335_8 = 4335.0 / 8.0;
    const auto f_4_13 = 4.0 / 13.0;
    const auto f_4_17 = 4.0 / 17.0;
    const auto f_500_11 = 500.0 / 11.0;
    const auto f_500_143 = 500.0 / 143.0;
    const auto f_5035_44 = 5035.0 / 44.0;
    const auto f_50_143 = 50.0 / 143.0;
    const auto f_5100_11 = 5100.0 / 11.0;
    const auto f_510_1 = 510.0;
    const auto f_513_11 = 513.0 / 11.0;
    const auto f_5292_143 = 5292.0 / 143.0;
    const auto f_56_187 = 56.0 / 187.0;
    const auto f_57_1 = 57.0;
    const auto f_608_11 = 608.0 / 11.0;
    const auto f_64_187 = 64.0 / 187.0;
    const auto f_68_143 = 68.0 / 143.0;
    const auto f_6_143 = 6.0 / 143.0;
    const auto f_6_17 = 6.0 / 17.0;
    const auto f_765_11 = 765.0 / 11.0;
    const auto f_765_143 = 765.0 / 143.0;
    const auto f_765_8 = 765.0 / 8.0;
    const auto f_798_11 = 798.0 / 11.0;
    const auto f_855_4 = 855.0 / 4.0;
    const auto f_8_143 = 8.0 / 143.0;
    const auto f_912_11 = 912.0 / 11.0;
    const auto f_95_4 = 95.0 / 4.0;
    const auto fs_1008_2717_14 = std::sqrt(14224896.0 / 7382089.0);
    const auto fs_101_2717_66 = std::sqrt(61206.0 / 671099.0);
    const auto fs_1020_11_15 = std::sqrt(15606000.0 / 121.0);
    const auto fs_1029_2717_286 = std::sqrt(2117682.0 / 51623.0);
    const auto fs_1029_572_286 = std::sqrt(1058841.0 / 1144.0);
    const auto fs_102_143_42 = std::sqrt(436968.0 / 20449.0);
    const auto fs_105_11_3 = std::sqrt(33075.0 / 121.0);
    const auto fs_105_286_2310 = std::sqrt(1157625.0 / 3718.0);
    const auto fs_105_4199_286 = std::sqrt(242550.0 / 1356277.0);
    const auto fs_10_143_14 = std::sqrt(1400.0 / 20449.0);
    const auto fs_10_143_5 = std::sqrt(500.0 / 20449.0);
    const auto fs_10_143_7 = std::sqrt(700.0 / 20449.0);
    const auto fs_10_2717_2310 = std::sqrt(21000.0 / 671099.0);
    const auto fs_10_429_105 = std::sqrt(3500.0 / 61347.0);
    const auto fs_10_429_15 = std::sqrt(500.0 / 61347.0);
    const auto fs_10_429_165 = std::sqrt(500.0 / 5577.0);
    const auto fs_10_429_210 = std::sqrt(7000.0 / 61347.0);
    const auto fs_10_429_35 = std::sqrt(3500.0 / 184041.0);
    const auto fs_10_429_70 = std::sqrt(7000.0 / 184041.0);
    const auto fs_10_561_210 = std::sqrt(7000.0 / 104907.0);
    const auto fs_111_2717_30 = std::sqrt(369630.0 / 7382089.0);
    const auto fs_1125_16_10 = std::sqrt(6328125.0 / 128.0);
    const auto fs_1125_16_2 = std::sqrt(1265625.0 / 128.0);
    const auto fs_1125_16_330 = std::sqrt(208828125.0 / 128.0);
    const auto fs_1125_16_70 = std::sqrt(44296875.0 / 128.0);
    const auto fs_1125_2_7 = std::sqrt(8859375.0 / 4.0);
    const auto fs_1125_4_30 = std::sqrt(18984375.0 / 8.0);
    const auto fs_1125_4_42 = std::sqrt(26578125.0 / 8.0);
    const auto fs_1125_8_105 = std::sqrt(132890625.0 / 64.0);
    const auto fs_1125_8_15 = std::sqrt(18984375.0 / 64.0);
    const auto fs_1125_8_165 = std::sqrt(208828125.0 / 64.0);
    const auto fs_1125_8_210 = std::sqrt(132890625.0 / 32.0);
    const auto fs_1125_8_35 = std::sqrt(44296875.0 / 64.0);
    const auto fs_1125_8_70 = std::sqrt(44296875.0 / 32.0);
    const auto fs_114_11_2 = std::sqrt(25992.0 / 121.0);
    const auto fs_114_11_33 = std::sqrt(38988.0 / 11.0);
    const auto fs_11676_2717_2 = std::sqrt(272657952.0 / 7382089.0);
    const auto fs_1176_2717_105 = std::sqrt(145212480.0 / 7382089.0);
    const auto fs_1176_2717_143 = std::sqrt(1382976.0 / 51623.0);
    const auto fs_1218_143_11 = std::sqrt(1483524.0 / 1859.0);
    const auto fs_1218_2717_21 = std::sqrt(31154004.0 / 7382089.0);
    const auto fs_12348_2717_3 = std::sqrt(457419312.0 / 7382089.0);
    const auto fs_125_11_105 = std::sqrt(1640625.0 / 121.0);
    const auto fs_125_11_15 = std::sqrt(234375.0 / 121.0);
    const auto fs_125_11_165 = std::sqrt(234375.0 / 11.0);
    const auto fs_125_11_210 = std::sqrt(3281250.0 / 121.0);
    const auto fs_125_11_35 = std::sqrt(546875.0 / 121.0);
    const auto fs_125_11_70 = std::sqrt(1093750.0 / 121.0);
    const auto fs_125_143_105 = std::sqrt(1640625.0 / 20449.0);
    const auto fs_125_143_15 = std::sqrt(234375.0 / 20449.0);
    const auto fs_125_143_165 = std::sqrt(234375.0 / 1859.0);
    const auto fs_125_143_210 = std::sqrt(3281250.0 / 20449.0);
    const auto fs_125_143_35 = std::sqrt(546875.0 / 20449.0);
    const auto fs_125_143_70 = std::sqrt(1093750.0 / 20449.0);
    const auto fs_125_1_30 = std::sqrt(468750.0);
    const auto fs_125_1_42 = std::sqrt(656250.0);
    const auto fs_125_22_10 = std::sqrt(78125.0 / 242.0);
    const auto fs_125_22_2 = std::sqrt(15625.0 / 242.0);
    const auto fs_125_22_330 = std::sqrt(234375.0 / 22.0);
    const auto fs_125_22_70 = std::sqrt(546875.0 / 242.0);
    const auto fs_125_286_10 = std::sqrt(78125.0 / 40898.0);
    const auto fs_125_286_2 = std::sqrt(15625.0 / 40898.0);
    const auto fs_125_286_330 = std::sqrt(234375.0 / 3718.0);
    const auto fs_125_286_70 = std::sqrt(546875.0 / 40898.0);
    const auto fs_125_2_105 = std::sqrt(1640625.0 / 4.0);
    const auto fs_125_2_15 = std::sqrt(234375.0 / 4.0);
    const auto fs_125_2_165 = std::sqrt(2578125.0 / 4.0);
    const auto fs_125_2_210 = std::sqrt(1640625.0 / 2.0);
    const auto fs_125_2_35 = std::sqrt(546875.0 / 4.0);
    const auto fs_125_2_70 = std::sqrt(546875.0 / 2.0);
    const auto fs_125_4_10 = std::sqrt(78125.0 / 8.0);
    const auto fs_125_4_2 = std::sqrt(15625.0 / 8.0);
    const auto fs_125_4_330 = std::sqrt(2578125.0 / 8.0);
    const auto fs_125_4_70 = std::sqrt(546875.0 / 8.0);
    const auto fs_126_209_11 = std::sqrt(15876.0 / 3971.0);
    const auto fs_126_247_273 = std::sqrt(333396.0 / 4693.0);
    const auto fs_126_247_30 = std::sqrt(476280.0 / 61009.0);
    const auto fs_126_247_7 = std::sqrt(111132.0 / 61009.0);
    const auto fs_126_2717_66 = std::sqrt(95256.0 / 671099.0);
    const auto fs_1275_4_3 = std::sqrt(4876875.0 / 16.0);
    const auto fs_127_2717_42 = std::sqrt(677418.0 / 7382089.0);
    const auto fs_12_143_11 = std::sqrt(144.0 / 1859.0);
    const auto fs_12_143_7 = std::sqrt(1008.0 / 20449.0);
    const auto fs_12_209_42 = std::sqrt(6048.0 / 43681.0);
    const auto fs_12_247_91 = std::sqrt(1008.0 / 4693.0);
    const auto fs_12_2717_210 = std::sqrt(30240.0 / 7382089.0);
    const auto fs_1323_286_66 = std::sqrt(5250987.0 / 3718.0);
    const auto fs_133_11_35 = std::sqrt(619115.0 / 121.0);
    const auto fs_133_187_7 = std::sqrt(123823.0 / 34969.0);
    const auto fs_133_22_7 = std::sqrt(123823.0 / 484.0);
    const auto fs_133_374_462 = std::sqrt(371469.0 / 6358.0);
    const auto fs_133_44_462 = std::sqrt(371469.0 / 88.0);
    const auto fs_136_143_30 = std::sqrt(554880.0 / 20449.0);
    const auto fs_136_143_66 = std::sqrt(110976.0 / 1859.0);
    const auto fs_13_187_10 = std::sqrt(1690.0 / 34969.0);
    const auto fs_147_13_13 = std::sqrt(21609.0 / 13.0);
    const auto fs_147_22_21 = std::sqrt(453789.0 / 484.0);
    const auto fs_147_286_110 = std::sqrt(108045.0 / 3718.0);
    const auto fs_147_286_4290 = std::sqrt(324135.0 / 286.0);
    const auto fs_14_187_14 = std::sqrt(2744.0 / 34969.0);
    const auto fs_14_187_3 = std::sqrt(588.0 / 34969.0);
    const auto fs_14_209_21 = std::sqrt(4116.0 / 43681.0);
    const auto fs_14_2717_110 = std::sqrt(1960.0 / 671099.0);
    const auto fs_14_2717_4290 = std::sqrt(5880.0 / 51623.0);
    const auto fs_14_4199_12155 = std::sqrt(10780.0 / 79781.0);
    const auto fs_14_4199_15015 = std::sqrt(226380.0 / 1356277.0);
    const auto fs_14_4199_2310 = std::sqrt(452760.0 / 17631601.0);
    const auto fs_14_4199_25194 = std::sqrt(1176.0 / 4199.0);
    const auto fs_14_4199_273 = std::sqrt(4116.0 / 1356277.0);
    const auto fs_14_4199_3003 = std::sqrt(45276.0 / 1356277.0);
    const auto fs_14_4199_3094 = std::sqrt(2744.0 / 79781.0);
    const auto fs_14_4199_62985 = std::sqrt(2940.0 / 4199.0);
    const auto fs_14_4199_9282 = std::sqrt(8232.0 / 79781.0);
    const auto fs_14_561_7 = std::sqrt(1372.0 / 314721.0);
    const auto fs_1530_11_11 = std::sqrt(2340900.0 / 11.0);
    const auto fs_1530_11_7 = std::sqrt(16386300.0 / 121.0);
    const auto fs_1530_143_30 = std::sqrt(70227000.0 / 20449.0);
    const auto fs_1530_143_66 = std::sqrt(14045400.0 / 1859.0);
    const auto fs_1554_2717_66 = std::sqrt(14489496.0 / 671099.0);
    const auto fs_1575_16_105 = std::sqrt(260465625.0 / 256.0);
    const auto fs_1575_16_15 = std::sqrt(37209375.0 / 256.0);
    const auto fs_1575_16_165 = std::sqrt(409303125.0 / 256.0);
    const auto fs_1575_16_210 = std::sqrt(260465625.0 / 128.0);
    const auto fs_1575_16_35 = std::sqrt(86821875.0 / 256.0);
    const auto fs_1575_16_70 = std::sqrt(86821875.0 / 128.0);
    const auto fs_1575_32_10 = std::sqrt(12403125.0 / 512.0);
    const auto fs_1575_32_2 = std::sqrt(2480625.0 / 512.0);
    const auto fs_1575_32_330 = std::sqrt(409303125.0 / 512.0);
    const auto fs_1575_32_70 = std::sqrt(86821875.0 / 512.0);
    const auto fs_1575_4_7 = std::sqrt(17364375.0 / 16.0);
    const auto fs_1575_8_30 = std::sqrt(37209375.0 / 32.0);
    const auto fs_1575_8_42 = std::sqrt(52093125.0 / 32.0);
    const auto fs_161_4199_12155 = std::sqrt(1425655.0 / 79781.0);
    const auto fs_161_4199_15015 = std::sqrt(29938755.0 / 1356277.0);
    const auto fs_161_4199_2310 = std::sqrt(59877510.0 / 17631601.0);
    const auto fs_161_4199_25194 = std::sqrt(155526.0 / 4199.0);
    const auto fs_161_4199_273 = std::sqrt(544341.0 / 1356277.0);
    const auto fs_161_4199_3003 = std::sqrt(5987751.0 / 1356277.0);
    const auto fs_161_4199_3094 = std::sqrt(362894.0 / 79781.0);
    const auto fs_161_4199_62985 = std::sqrt(388815.0 / 4199.0);
    const auto fs_161_4199_9282 = std::sqrt(1088682.0 / 79781.0);
    const auto fs_161_8398_182 = std::sqrt(181447.0 / 2712554.0);
    const auto fs_161_8398_2 = std::sqrt(25921.0 / 35263202.0);
    const auto fs_161_8398_26 = std::sqrt(25921.0 / 2712554.0);
    const auto fs_161_8398_330 = std::sqrt(4276965.0 / 35263202.0);
    const auto fs_161_8398_910 = std::sqrt(907235.0 / 2712554.0);
    const auto fs_168_2717_15015 = std::sqrt(2963520.0 / 51623.0);
    const auto fs_168_2717_5 = std::sqrt(141120.0 / 7382089.0);
    const auto fs_168_2717_55 = std::sqrt(141120.0 / 671099.0);
    const auto fs_168_2717_6 = std::sqrt(169344.0 / 7382089.0);
    const auto fs_168_4199_154 = std::sqrt(4346496.0 / 17631601.0);
    const auto fs_168_4199_210 = std::sqrt(5927040.0 / 17631601.0);
    const auto fs_168_4199_22 = std::sqrt(620928.0 / 17631601.0);
    const auto fs_168_4199_221 = std::sqrt(28224.0 / 79781.0);
    const auto fs_171_374_154 = std::sqrt(204687.0 / 6358.0);
    const auto fs_171_374_330 = std::sqrt(438615.0 / 6358.0);
    const auto fs_171_44_154 = std::sqrt(204687.0 / 88.0);
    const auto fs_171_44_330 = std::sqrt(438615.0 / 88.0);
    const auto fs_175_2717_10 = std::sqrt(306250.0 / 7382089.0);
    const auto fs_1764_2717_143 = std::sqrt(3111696.0 / 51623.0);
    const auto fs_1764_2717_77 = std::sqrt(21781872.0 / 671099.0);
    const auto fs_1806_2717_33 = std::sqrt(9784908.0 / 671099.0);
    const auto fs_188_2717_21 = std::sqrt(742224.0 / 7382089.0);
    const auto fs_1932_4199_154 = std::sqrt(574824096.0 / 17631601.0);
    const auto fs_1932_4199_210 = std::sqrt(783851040.0 / 17631601.0);
    const auto fs_1932_4199_22 = std::sqrt(82117728.0 / 17631601.0);
    const auto fs_1932_4199_221 = std::sqrt(3732624.0 / 79781.0);
    const auto fs_1995_22_3 = std::sqrt(11940075.0 / 484.0);
    const auto fs_1995_22_5 = std::sqrt(19900125.0 / 484.0);
    const auto fs_1995_44_14 = std::sqrt(27860175.0 / 968.0);
    const auto fs_1995_44_3 = std::sqrt(11940075.0 / 1936.0);
    const auto fs_19_11_210 = std::sqrt(75810.0 / 121.0);
    const auto fs_19_17_42 = std::sqrt(15162.0 / 289.0);
    const auto fs_19_187_2310 = std::sqrt(75810.0 / 3179.0);
    const auto fs_19_22_2310 = std::sqrt(37905.0 / 22.0);
    const auto fs_19_2_42 = std::sqrt(7581.0 / 2.0);
    const auto fs_1_11_6 = std::sqrt(6.0 / 121.0);
    const auto fs_1_17_30 = std::sqrt(30.0 / 289.0);
    const auto fs_1_247_182 = std::sqrt(14.0 / 4693.0);
    const auto fs_1_247_2730 = std::sqrt(210.0 / 4693.0);
    const auto fs_1_2717_1430 = std::sqrt(10.0 / 51623.0);
    const auto fs_2016_2717_55 = std::sqrt(20321280.0 / 671099.0);
    const auto fs_204_143_11 = std::sqrt(41616.0 / 1859.0);
    const auto fs_204_143_30 = std::sqrt(1248480.0 / 20449.0);
    const auto fs_20_143_3 = std::sqrt(1200.0 / 20449.0);
    const auto fs_20_209_3 = std::sqrt(1200.0 / 43681.0);
    const auto fs_20_429_30 = std::sqrt(4000.0 / 61347.0);
    const auto fs_20_429_42 = std::sqrt(5600.0 / 61347.0);
    const auto fs_210_2717_10 = std::sqrt(441000.0 / 7382089.0);
    const auto fs_210_2717_2310 = std::sqrt(9261000.0 / 671099.0);
    const auto fs_210_4199_154 = std::sqrt(6791400.0 / 17631601.0);
    const auto fs_210_4199_42 = std::sqrt(1852200.0 / 17631601.0);
    const auto fs_2121_2717_66 = std::sqrt(26991846.0 / 671099.0);
    const auto fs_2121_572_66 = std::sqrt(13495923.0 / 14872.0);
    const auto fs_2185_88_6 = std::sqrt(14322675.0 / 3872.0);
    const auto fs_21_13_195 = std::sqrt(6615.0 / 13.0);
    const auto fs_21_13_546 = std::sqrt(18522.0 / 13.0);
    const auto fs_21_143_231 = std::sqrt(9261.0 / 1859.0);
    const auto fs_21_143_30030 = std::sqrt(92610.0 / 143.0);
    const auto fs_21_247_182 = std::sqrt(6174.0 / 4693.0);
    const auto fs_21_247_2730 = std::sqrt(92610.0 / 4693.0);
    const auto fs_21_2717_1430 = std::sqrt(4410.0 / 51623.0);
    const auto fs_21_2717_286 = std::sqrt(882.0 / 51623.0);
    const auto fs_21_286_3003 = std::sqrt(9261.0 / 572.0);
    const auto fs_21_4199_10010 = std::sqrt(339570.0 / 1356277.0);
    const auto fs_21_4199_110 = std::sqrt(48510.0 / 17631601.0);
    const auto fs_21_4199_1430 = std::sqrt(48510.0 / 1356277.0);
    const auto fs_21_4199_2 = std::sqrt(882.0 / 17631601.0);
    const auto fs_21_4199_286 = std::sqrt(9702.0 / 1356277.0);
    const auto fs_21_52_182 = std::sqrt(3087.0 / 104.0);
    const auto fs_21_52_2730 = std::sqrt(46305.0 / 104.0);
    const auto fs_21_572_1430 = std::sqrt(2205.0 / 1144.0);
    const auto fs_2205_286_10 = std::sqrt(24310125.0 / 40898.0);
    const auto fs_228_187_2 = std::sqrt(103968.0 / 34969.0);
    const auto fs_228_187_33 = std::sqrt(155952.0 / 3179.0);
    const auto fs_2295_143_11 = std::sqrt(5267025.0 / 1859.0);
    const auto fs_2295_143_30 = std::sqrt(158010750.0 / 20449.0);
    const auto fs_2295_286_42 = std::sqrt(110607525.0 / 40898.0);
    const auto fs_22_247_26 = std::sqrt(968.0 / 4693.0);
    const auto fs_231_26_26 = std::sqrt(53361.0 / 26.0);
    const auto fs_232_2717_11 = std::sqrt(53824.0 / 671099.0);
    const auto fs_2331_2717_30 = std::sqrt(163006830.0 / 7382089.0);
    const auto fs_2331_572_30 = std::sqrt(81503415.0 / 163592.0);
    const auto fs_23_247_6 = std::sqrt(3174.0 / 61009.0);
    const auto fs_23_561_6 = std::sqrt(1058.0 / 104907.0);
    const auto fs_2415_4199_154 = std::sqrt(898162650.0 / 17631601.0);
    const auto fs_2415_4199_42 = std::sqrt(244953450.0 / 17631601.0);
    const auto fs_2415_8398_286 = std::sqrt(64154475.0 / 2712554.0);
    const auto fs_2499_143_3 = std::sqrt(18735003.0 / 20449.0);
    const auto fs_24_143_5 = std::sqrt(2880.0 / 20449.0);
    const auto fs_250_11_30 = std::sqrt(1875000.0 / 121.0);
    const auto fs_250_11_42 = std::sqrt(2625000.0 / 121.0);
    const auto fs_250_143_30 = std::sqrt(1875000.0 / 20449.0);
    const auto fs_250_143_42 = std::sqrt(2625000.0 / 20449.0);
    const auto fs_250_1_7 = std::sqrt(437500.0);
    const auto fs_252_143_14 = std::sqrt(889056.0 / 20449.0);
    const auto fs_252_209_42 = std::sqrt(2667168.0 / 43681.0);
    const auto fs_252_247_91 = std::sqrt(444528.0 / 4693.0);
    const auto fs_252_2717_210 = std::sqrt(13335840.0 / 7382089.0);
    const auto fs_2550_11_3 = std::sqrt(19507500.0 / 121.0);
    const auto fs_255_11_105 = std::sqrt(6827625.0 / 121.0);
    const auto fs_255_11_210 = std::sqrt(13655250.0 / 121.0);
    const auto fs_255_11_42 = std::sqrt(2731050.0 / 121.0);
    const auto fs_255_11_462 = std::sqrt(2731050.0 / 11.0);
    const auto fs_255_2_15 = std::sqrt(975375.0 / 4.0);
    const auto fs_255_4_30 = std::sqrt(975375.0 / 8.0);
    const auto fs_255_4_66 = std::sqrt(2145825.0 / 8.0);
    const auto fs_255_8_105 = std::sqrt(6827625.0 / 64.0);
    const auto fs_255_8_210 = std::sqrt(6827625.0 / 32.0);
    const auto fs_255_8_42 = std::sqrt(1365525.0 / 32.0);
    const auto fs_255_8_462 = std::sqrt(15020775.0 / 32.0);
    const auto fs_2646_2717_66 = std::sqrt(42007896.0 / 671099.0);
    const auto fs_2667_2717_42 = std::sqrt(298741338.0 / 7382089.0);
    const auto fs_2667_572_42 = std::sqrt(149370669.0 / 163592.0);
    const auto fs_266_11_5 = std::sqrt(353780.0 / 121.0);
    const auto fs_266_187_35 = std::sqrt(2476460.0 / 34969.0);
    const auto fs_272_143_15 = std::sqrt(1109760.0 / 20449.0);
    const auto fs_27_2717_858 = std::sqrt(4374.0 / 51623.0);
    const auto fs_280_4199_3 = std::sqrt(235200.0 / 17631601.0);
    const auto fs_285_11_2 = std::sqrt(162450.0 / 121.0);
    const auto fs_285_11_33 = std::sqrt(243675.0 / 11.0);
    const auto fs_285_44_55 = std::sqrt(406125.0 / 176.0);
    const auto fs_285_44_77 = std::sqrt(568575.0 / 176.0);
    const auto fs_285_8_30 = std::sqrt(1218375.0 / 32.0);
    const auto fs_28_187_3 = std::sqrt(2352.0 / 34969.0);
    const auto fs_28_187_5 = std::sqrt(3920.0 / 34969.0);
    const auto fs_28_247_13 = std::sqrt(784.0 / 4693.0);
    const auto fs_28_4199_231 = std::sqrt(181104.0 / 17631601.0);
    const auto fs_28_4199_273 = std::sqrt(16464.0 / 1356277.0);
    const auto fs_28_4199_3003 = std::sqrt(181104.0 / 1356277.0);
    const auto fs_28_4199_455 = std::sqrt(27440.0 / 1356277.0);
    const auto fs_28_4199_4641 = std::sqrt(16464.0 / 79781.0);
    const auto fs_28_4199_5005 = std::sqrt(301840.0 / 1356277.0);
    const auto fs_28_4199_6006 = std::sqrt(362208.0 / 1356277.0);
    const auto fs_28_4199_7293 = std::sqrt(25872.0 / 79781.0);
    const auto fs_28_561_35 = std::sqrt(27440.0 / 314721.0);
    const auto fs_2919_143_2 = std::sqrt(17041122.0 / 20449.0);
    const auto fs_294_143_105 = std::sqrt(9075780.0 / 20449.0);
    const auto fs_294_143_143 = std::sqrt(86436.0 / 143.0);
    const auto fs_294_209_21 = std::sqrt(1815156.0 / 43681.0);
    const auto fs_294_2717_110 = std::sqrt(864360.0 / 671099.0);
    const auto fs_294_2717_4290 = std::sqrt(2593080.0 / 51623.0);
    const auto fs_29_2717_154 = std::sqrt(11774.0 / 671099.0);
    const auto fs_2_143_105 = std::sqrt(420.0 / 20449.0);
    const auto fs_2_143_210 = std::sqrt(840.0 / 20449.0);
    const auto fs_2_143_42 = std::sqrt(168.0 / 20449.0);
    const auto fs_2_143_462 = std::sqrt(168.0 / 1859.0);
    const auto fs_2_187_55 = std::sqrt(20.0 / 3179.0);
    const auto fs_2_187_77 = std::sqrt(28.0 / 3179.0);
    const auto fs_2_2717_3003 = std::sqrt(84.0 / 51623.0);
    const auto fs_2_51_42 = std::sqrt(56.0 / 867.0);
    const auto fs_2_561_2310 = std::sqrt(280.0 / 9537.0);
    const auto fs_3060_11_5 = std::sqrt(46818000.0 / 121.0);
    const auto fs_3060_143_15 = std::sqrt(140454000.0 / 20449.0);
    const auto fs_3087_143_3 = std::sqrt(28588707.0 / 20449.0);
    const auto fs_30_2717_22 = std::sqrt(1800.0 / 671099.0);
    const auto fs_315_286_22 = std::sqrt(99225.0 / 3718.0);
    const auto fs_3220_4199_3 = std::sqrt(31105200.0 / 17631601.0);
    const auto fs_322_4199_231 = std::sqrt(23951004.0 / 17631601.0);
    const auto fs_322_4199_273 = std::sqrt(2177364.0 / 1356277.0);
    const auto fs_322_4199_3003 = std::sqrt(23951004.0 / 1356277.0);
    const auto fs_322_4199_455 = std::sqrt(3628940.0 / 1356277.0);
    const auto fs_322_4199_4641 = std::sqrt(2177364.0 / 79781.0);
    const auto fs_322_4199_5005 = std::sqrt(39918340.0 / 1356277.0);
    const auto fs_322_4199_6006 = std::sqrt(47902008.0 / 1356277.0);
    const auto fs_322_4199_7293 = std::sqrt(3421572.0 / 79781.0);
    const auto fs_3315_16_6 = std::sqrt(32967675.0 / 128.0);
    const auto fs_3315_22_6 = std::sqrt(32967675.0 / 242.0);
    const auto fs_3375_16_110 = std::sqrt(626484375.0 / 128.0);
    const auto fs_3375_16_2 = std::sqrt(11390625.0 / 128.0);
    const auto fs_3375_4_3 = std::sqrt(34171875.0 / 16.0);
    const auto fs_3375_8_14 = std::sqrt(79734375.0 / 32.0);
    const auto fs_3375_8_5 = std::sqrt(56953125.0 / 64.0);
    const auto fs_3375_8_7 = std::sqrt(79734375.0 / 64.0);
    const auto fs_34_11_6 = std::sqrt(6936.0 / 121.0);
    const auto fs_34_247_5 = std::sqrt(5780.0 / 61009.0);
    const auto fs_3528_2717_5 = std::sqrt(62233920.0 / 7382089.0);
    const auto fs_3528_2717_55 = std::sqrt(62233920.0 / 671099.0);
    const auto fs_3528_2717_6 = std::sqrt(74680704.0 / 7382089.0);
    const auto fs_357_26_5 = std::sqrt(637245.0 / 676.0);
    const auto fs_3675_2717_10 = std::sqrt(135056250.0 / 7382089.0);
    const auto fs_3675_572_10 = std::sqrt(67528125.0 / 163592.0);
    const auto fs_3705_88_10 = std::sqrt(68635125.0 / 3872.0);
    const auto fs_375_11_14 = std::sqrt(1968750.0 / 121.0);
    const auto fs_375_11_5 = std::sqrt(703125.0 / 121.0);
    const auto fs_375_11_7 = std::sqrt(984375.0 / 121.0);
    const auto fs_375_143_14 = std::sqrt(1968750.0 / 20449.0);
    const auto fs_375_143_5 = std::sqrt(703125.0 / 20449.0);
    const auto fs_375_143_7 = std::sqrt(984375.0 / 20449.0);
    const auto fs_375_1_3 = std::sqrt(421875.0);
    const auto fs_375_22_110 = std::sqrt(703125.0 / 22.0);
    const auto fs_375_22_2 = std::sqrt(140625.0 / 242.0);
    const auto fs_375_286_110 = std::sqrt(703125.0 / 3718.0);
    const auto fs_375_286_2 = std::sqrt(140625.0 / 40898.0);
    const auto fs_375_2_14 = std::sqrt(984375.0 / 2.0);
    const auto fs_375_2_5 = std::sqrt(703125.0 / 4.0);
    const auto fs_375_2_7 = std::sqrt(984375.0 / 4.0);
    const auto fs_375_4_110 = std::sqrt(7734375.0 / 8.0);
    const auto fs_375_4_2 = std::sqrt(140625.0 / 8.0);
    const auto fs_38_187_210 = std::sqrt(303240.0 / 34969.0);
    const auto fs_392_4199_33 = std::sqrt(5070912.0 / 17631601.0);
    const auto fs_3948_2717_21 = std::sqrt(327320784.0 / 7382089.0);
    const auto fs_398_2717_3 = std::sqrt(475212.0 / 7382089.0);
    const auto fs_399_11_3 = std::sqrt(477603.0 / 121.0);
    const auto fs_399_11_5 = std::sqrt(796005.0 / 121.0);
    const auto fs_399_187_14 = std::sqrt(2228814.0 / 34969.0);
    const auto fs_399_187_3 = std::sqrt(477603.0 / 34969.0);
    const auto fs_399_22_14 = std::sqrt(1114407.0 / 242.0);
    const auto fs_399_22_3 = std::sqrt(477603.0 / 484.0);
    const auto fs_3_143_42 = std::sqrt(378.0 / 20449.0);
    const auto fs_3_187_154 = std::sqrt(126.0 / 3179.0);
    const auto fs_3_187_330 = std::sqrt(270.0 / 3179.0);
    const auto fs_408_143_11 = std::sqrt(166464.0 / 1859.0);
    const auto fs_408_143_7 = std::sqrt(1165248.0 / 20449.0);
    const auto fs_40_429_7 = std::sqrt(11200.0 / 184041.0);
    const auto fs_4179_286_3 = std::sqrt(52392123.0 / 81796.0);
    const auto fs_420_209_3 = std::sqrt(529200.0 / 43681.0);
    const auto fs_42_143_15015 = std::sqrt(185220.0 / 143.0);
    const auto fs_42_2717_3003 = std::sqrt(37044.0 / 51623.0);
    const auto fs_42_4199_1155 = std::sqrt(2037420.0 / 17631601.0);
    const auto fs_42_4199_165 = std::sqrt(291060.0 / 17631601.0);
    const auto fs_42_4199_2002 = std::sqrt(271656.0 / 1356277.0);
    const auto fs_42_4199_22 = std::sqrt(38808.0 / 17631601.0);
    const auto fs_42_4199_2431 = std::sqrt(19404.0 / 79781.0);
    const auto fs_42_4199_4199 = std::sqrt(1764.0 / 4199.0);
    const auto fs_42_4199_5 = std::sqrt(8820.0 / 17631601.0);
    const auto fs_42_4199_715 = std::sqrt(97020.0 / 1356277.0);
    const auto fs_437_374_6 = std::sqrt(572907.0 / 69938.0);
    const auto fs_437_44_6 = std::sqrt(572907.0 / 968.0);
    const auto fs_4410_2717_10 = std::sqrt(194481000.0 / 7382089.0);
    const auto fs_441_143_143 = std::sqrt(194481.0 / 143.0);
    const auto fs_441_143_77 = std::sqrt(1361367.0 / 1859.0);
    const auto fs_441_2717_286 = std::sqrt(388962.0 / 51623.0);
    const auto fs_441_572_286 = std::sqrt(194481.0 / 1144.0);
    const auto fs_4508_4199_33 = std::sqrt(670628112.0 / 17631601.0);
    const auto fs_4590_143_11 = std::sqrt(21068100.0 / 1859.0);
    const auto fs_4590_143_7 = std::sqrt(147476700.0 / 20449.0);
    const auto fs_462_247_26 = std::sqrt(426888.0 / 4693.0);
    const auto fs_46_2717_105 = std::sqrt(222180.0 / 7382089.0);
    const auto fs_4725_16_14 = std::sqrt(156279375.0 / 128.0);
    const auto fs_4725_16_5 = std::sqrt(111628125.0 / 256.0);
    const auto fs_4725_16_7 = std::sqrt(156279375.0 / 256.0);
    const auto fs_4725_32_110 = std::sqrt(1227909375.0 / 512.0);
    const auto fs_4725_32_2 = std::sqrt(22325625.0 / 512.0);
    const auto fs_4725_8_3 = std::sqrt(66976875.0 / 64.0);
    const auto fs_475_44_210 = std::sqrt(23690625.0 / 968.0);
    const auto fs_476_2717_3 = std::sqrt(679728.0 / 7382089.0);
    const auto fs_483_247_6 = std::sqrt(1399734.0 / 61009.0);
    const auto fs_483_286_105 = std::sqrt(24495345.0 / 81796.0);
    const auto fs_483_4199_1155 = std::sqrt(269448795.0 / 17631601.0);
    const auto fs_483_4199_165 = std::sqrt(38492685.0 / 17631601.0);
    const auto fs_483_4199_2002 = std::sqrt(35926506.0 / 1356277.0);
    const auto fs_483_4199_22 = std::sqrt(5132358.0 / 17631601.0);
    const auto fs_483_4199_2431 = std::sqrt(2566179.0 / 79781.0);
    const auto fs_483_4199_4199 = std::sqrt(233289.0 / 4199.0);
    const auto fs_483_4199_5 = std::sqrt(1166445.0 / 17631601.0);
    const auto fs_483_4199_715 = std::sqrt(12830895.0 / 1356277.0);
    const auto fs_483_52_6 = std::sqrt(699867.0 / 1352.0);
    const auto fs_483_8398_10010 = std::sqrt(89816265.0 / 2712554.0);
    const auto fs_483_8398_110 = std::sqrt(12830895.0 / 35263202.0);
    const auto fs_483_8398_1430 = std::sqrt(12830895.0 / 2712554.0);
    const auto fs_483_8398_2 = std::sqrt(233289.0 / 35263202.0);
    const auto fs_483_8398_286 = std::sqrt(2566179.0 / 2712554.0);
    const auto fs_4872_2717_11 = std::sqrt(23736384.0 / 671099.0);
    const auto fs_48_2717_14 = std::sqrt(32256.0 / 7382089.0);
    const auto fs_49_2717_286 = std::sqrt(4802.0 / 51623.0);
    const auto fs_4_143_30 = std::sqrt(480.0 / 20449.0);
    const auto fs_4_143_66 = std::sqrt(96.0 / 1859.0);
    const auto fs_4_247_195 = std::sqrt(240.0 / 4693.0);
    const auto fs_4_247_546 = std::sqrt(672.0 / 4693.0);
    const auto fs_4_2717_231 = std::sqrt(336.0 / 671099.0);
    const auto fs_4_2717_30030 = std::sqrt(3360.0 / 51623.0);
    const auto fs_4_561_210 = std::sqrt(1120.0 / 104907.0);
    const auto fs_500_11_7 = std::sqrt(1750000.0 / 121.0);
    const auto fs_500_143_7 = std::sqrt(1750000.0 / 20449.0);
    const auto fs_504_143_55 = std::sqrt(1270080.0 / 1859.0);
    const auto fs_50_429_7 = std::sqrt(17500.0 / 184041.0);
    const auto fs_510_11_30 = std::sqrt(7803000.0 / 121.0);
    const auto fs_510_11_66 = std::sqrt(1560600.0 / 11.0);
    const auto fs_532_187_5 = std::sqrt(1415120.0 / 34969.0);
    const auto fs_556_2717_2 = std::sqrt(618272.0 / 7382089.0);
    const auto fs_5625_8_7 = std::sqrt(221484375.0 / 64.0);
    const auto fs_567_2717_858 = std::sqrt(1928934.0 / 51623.0);
    const auto fs_567_572_858 = std::sqrt(964467.0 / 1144.0);
    const auto fs_56_2717_105 = std::sqrt(329280.0 / 7382089.0);
    const auto fs_56_2717_143 = std::sqrt(3136.0 / 51623.0);
    const auto fs_56_4199_1430 = std::sqrt(344960.0 / 1356277.0);
    const auto fs_56_4199_3 = std::sqrt(9408.0 / 17631601.0);
    const auto fs_56_4199_385 = std::sqrt(1207360.0 / 17631601.0);
    const auto fs_56_4199_546 = std::sqrt(131712.0 / 1356277.0);
    const auto fs_56_561_5 = std::sqrt(15680.0 / 314721.0);
    const auto fs_57_187_55 = std::sqrt(16245.0 / 3179.0);
    const auto fs_57_187_77 = std::sqrt(22743.0 / 3179.0);
    const auto fs_57_22_55 = std::sqrt(16245.0 / 44.0);
    const auto fs_57_22_77 = std::sqrt(22743.0 / 44.0);
    const auto fs_57_34_30 = std::sqrt(48735.0 / 578.0);
    const auto fs_57_4_30 = std::sqrt(48735.0 / 8.0);
    const auto fs_588_247_13 = std::sqrt(345744.0 / 4693.0);
    const auto fs_588_2717_3 = std::sqrt(1037232.0 / 7382089.0);
    const auto fs_588_4199_22 = std::sqrt(7606368.0 / 17631601.0);
    const auto fs_58_2717_21 = std::sqrt(70644.0 / 7382089.0);
    const auto fs_5_143_110 = std::sqrt(250.0 / 1859.0);
    const auto fs_5_143_2 = std::sqrt(50.0 / 20449.0);
    const auto fs_5_429_10 = std::sqrt(250.0 / 184041.0);
    const auto fs_5_429_2 = std::sqrt(50.0 / 184041.0);
    const auto fs_5_429_330 = std::sqrt(250.0 / 5577.0);
    const auto fs_5_429_70 = std::sqrt(1750.0 / 184041.0);
    const auto fs_609_2717_154 = std::sqrt(5192334.0 / 671099.0);
    const auto fs_609_286_21 = std::sqrt(7788501.0 / 81796.0);
    const auto fs_609_572_154 = std::sqrt(2596167.0 / 14872.0);
    const auto fs_625_11_7 = std::sqrt(2734375.0 / 121.0);
    const auto fs_625_143_7 = std::sqrt(2734375.0 / 20449.0);
    const auto fs_625_2_7 = std::sqrt(2734375.0 / 4.0);
    const auto fs_630_2717_22 = std::sqrt(793800.0 / 671099.0);
    const auto fs_63_11_42 = std::sqrt(166698.0 / 121.0);
    const auto fs_63_13_91 = std::sqrt(27783.0 / 13.0);
    const auto fs_63_143_210 = std::sqrt(833490.0 / 20449.0);
    const auto fs_63_22_11 = std::sqrt(3969.0 / 44.0);
    const auto fs_63_26_273 = std::sqrt(83349.0 / 52.0);
    const auto fs_63_26_30 = std::sqrt(59535.0 / 338.0);
    const auto fs_63_26_7 = std::sqrt(27783.0 / 676.0);
    const auto fs_644_4199_1430 = std::sqrt(45620960.0 / 1356277.0);
    const auto fs_644_4199_3 = std::sqrt(1244208.0 / 17631601.0);
    const auto fs_644_4199_385 = std::sqrt(159673360.0 / 17631601.0);
    const auto fs_644_4199_546 = std::sqrt(17418912.0 / 1356277.0);
    const auto fs_665_11_5 = std::sqrt(2211125.0 / 121.0);
    const auto fs_665_22_35 = std::sqrt(15477875.0 / 484.0);
    const auto fs_665_44_7 = std::sqrt(3095575.0 / 1936.0);
    const auto fs_665_88_462 = std::sqrt(9286725.0 / 352.0);
    const auto fs_6762_4199_22 = std::sqrt(1005942168.0 / 17631601.0);
    const auto fs_680_143_3 = std::sqrt(1387200.0 / 20449.0);
    const auto fs_68_143_105 = std::sqrt(485520.0 / 20449.0);
    const auto fs_68_143_210 = std::sqrt(971040.0 / 20449.0);
    const auto fs_68_143_42 = std::sqrt(194208.0 / 20449.0);
    const auto fs_68_143_462 = std::sqrt(194208.0 / 1859.0);
    const auto fs_6_143_11 = std::sqrt(36.0 / 1859.0);
    const auto fs_6_143_30 = std::sqrt(1080.0 / 20449.0);
    const auto fs_6_209_11 = std::sqrt(36.0 / 3971.0);
    const auto fs_6_247_273 = std::sqrt(756.0 / 4693.0);
    const auto fs_6_247_30 = std::sqrt(1080.0 / 61009.0);
    const auto fs_6_247_7 = std::sqrt(252.0 / 61009.0);
    const auto fs_70_4199_1001 = std::sqrt(377300.0 / 1356277.0);
    const auto fs_70_4199_273 = std::sqrt(102900.0 / 1356277.0);
    const auto fs_714_247_5 = std::sqrt(2548980.0 / 61009.0);
    const auto fs_741_374_10 = std::sqrt(2745405.0 / 69938.0);
    const auto fs_741_44_10 = std::sqrt(2745405.0 / 968.0);
    const auto fs_74_2717_66 = std::sqrt(32856.0 / 671099.0);
    const auto fs_750_11_3 = std::sqrt(1687500.0 / 121.0);
    const auto fs_750_143_3 = std::sqrt(1687500.0 / 20449.0);
    const auto fs_7650_143_3 = std::sqrt(175567500.0 / 20449.0);
    const auto fs_765_11_11 = std::sqrt(585225.0 / 11.0);
    const auto fs_765_11_30 = std::sqrt(17556750.0 / 121.0);
    const auto fs_765_143_105 = std::sqrt(61448625.0 / 20449.0);
    const auto fs_765_143_210 = std::sqrt(122897250.0 / 20449.0);
    const auto fs_765_143_42 = std::sqrt(24579450.0 / 20449.0);
    const auto fs_765_143_462 = std::sqrt(24579450.0 / 1859.0);
    const auto fs_765_16_42 = std::sqrt(12289725.0 / 128.0);
    const auto fs_765_22_42 = std::sqrt(12289725.0 / 242.0);
    const auto fs_765_22_6 = std::sqrt(1755675.0 / 242.0);
    const auto fs_765_2_5 = std::sqrt(2926125.0 / 4.0);
    const auto fs_765_4_11 = std::sqrt(6437475.0 / 16.0);
    const auto fs_765_4_7 = std::sqrt(4096575.0 / 16.0);
    const auto fs_765_8_11 = std::sqrt(6437475.0 / 64.0);
    const auto fs_765_8_30 = std::sqrt(8778375.0 / 32.0);
    const auto fs_777_286_66 = std::sqrt(1811187.0 / 3718.0);
    const auto fs_7875_16_7 = std::sqrt(434109375.0 / 256.0);
    const auto fs_798_187_3 = std::sqrt(1910412.0 / 34969.0);
    const auto fs_798_187_5 = std::sqrt(3184020.0 / 34969.0);
    const auto fs_7_4199_182 = std::sqrt(686.0 / 1356277.0);
    const auto fs_7_4199_2 = std::sqrt(98.0 / 17631601.0);
    const auto fs_7_4199_26 = std::sqrt(98.0 / 1356277.0);
    const auto fs_7_4199_330 = std::sqrt(16170.0 / 17631601.0);
    const auto fs_7_4199_910 = std::sqrt(3430.0 / 1356277.0);
    const auto fs_7_561_462 = std::sqrt(686.0 / 9537.0);
    const auto fs_805_4199_1001 = std::sqrt(49897925.0 / 1356277.0);
    const auto fs_805_4199_273 = std::sqrt(13608525.0 / 1356277.0);
    const auto fs_816_143_5 = std::sqrt(3329280.0 / 20449.0);
    const auto fs_8358_2717_3 = std::sqrt(209568492.0 / 7382089.0);
    const auto fs_84_247_195 = std::sqrt(105840.0 / 4693.0);
    const auto fs_84_247_546 = std::sqrt(296352.0 / 4693.0);
    const auto fs_84_2717_143 = std::sqrt(7056.0 / 51623.0);
    const auto fs_84_2717_231 = std::sqrt(148176.0 / 671099.0);
    const auto fs_84_2717_30030 = std::sqrt(1481760.0 / 51623.0);
    const auto fs_84_2717_77 = std::sqrt(49392.0 / 671099.0);
    const auto fs_84_4199_1155 = std::sqrt(8149680.0 / 17631601.0);
    const auto fs_84_4199_231 = std::sqrt(1629936.0 / 17631601.0);
    const auto fs_84_4199_286 = std::sqrt(155232.0 / 1356277.0);
    const auto fs_855_88_154 = std::sqrt(5117175.0 / 352.0);
    const auto fs_855_88_330 = std::sqrt(10965375.0 / 352.0);
    const auto fs_86_2717_33 = std::sqrt(22188.0 / 671099.0);
    const auto fs_882_143_5 = std::sqrt(3889620.0 / 20449.0);
    const auto fs_882_143_55 = std::sqrt(3889620.0 / 1859.0);
    const auto fs_882_143_6 = std::sqrt(4667544.0 / 20449.0);
    const auto fs_8_143_15 = std::sqrt(960.0 / 20449.0);
    const auto fs_8_187_2 = std::sqrt(128.0 / 34969.0);
    const auto fs_8_187_33 = std::sqrt(192.0 / 3179.0);
    const auto fs_8_2717_15015 = std::sqrt(6720.0 / 51623.0);
    const auto fs_903_286_33 = std::sqrt(2446227.0 / 7436.0);
    const auto fs_9180_143_5 = std::sqrt(421362000.0 / 20449.0);
    const auto fs_95_187_210 = std::sqrt(1895250.0 / 34969.0);
    const auto fs_95_22_210 = std::sqrt(947625.0 / 242.0);
    const auto fs_95_44_2310 = std::sqrt(947625.0 / 88.0);
    const auto fs_95_4_42 = std::sqrt(189525.0 / 8.0);
    const auto fs_966_2717_105 = std::sqrt(97981380.0 / 7382089.0);
    const auto fs_966_4199_1155 = std::sqrt(1077795180.0 / 17631601.0);
    const auto fs_966_4199_231 = std::sqrt(215559036.0 / 17631601.0);
    const auto fs_966_4199_286 = std::sqrt(20529432.0 / 1356277.0);
    const auto fs_96_2717_55 = std::sqrt(46080.0 / 671099.0);
    const auto fs_987_143_21 = std::sqrt(20457549.0 / 20449.0);
    const auto fs_9996_2717_3 = std::sqrt(299760048.0 / 7382089.0);

    // NOTE: the angular components are formed in 93 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph2_p2, ph4_p2, ph4_p3, ph6_p2, ph6_p3, ph8_p2, ph8_p3, ph10_p2, ph10_p3, ph10_p9, ph10_p10, ab_2, pc_0, pc_1 : simd::cache_line_size())
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
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;

        const auto h2_p2 = ph2_p2[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h4_p3 = ph4_p3[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h8_p2 = ph8_p2[k];
        const auto h8_p3 = ph8_p3[k];
        const auto h10_p2 = ph10_p2[k];
        const auto h10_p3 = ph10_p3[k];
        const auto h10_p9 = ph10_p9[k];
        const auto h10_p10 = ph10_p10[k];

        pc_0[k] = e_0 * fs_4725_32_110 * h2_p2 + e_1 * fs_255_4_66 * h4_p2 - e_1 * fs_3375_16_110 * r_2 * h2_p2 + e_2 * fs_285_44_77 * h6_p2 - e_2 * fs_510_11_66 * r_2 * h4_p2 + e_2 * fs_375_4_110 * r_4 * h2_p2 + e_3 * fs_21_143_231 * h8_p2 - e_3 * fs_57_22_77 * r_2 * h6_p2 + e_3 * fs_1530_143_66 * r_4 * h4_p2 - e_3 * fs_375_22_110 * r_6 * h2_p2 + e_4 * fs_161_8398_2 * h10_p2 - e_4 * fs_161_4199_62985 * h10_p10 - e_4 * fs_84_2717_231 * r_2 * h8_p2 + e_4 * fs_57_187_77 * r_4 * h6_p2 - e_4 * fs_136_143_66 * r_6 * h4_p2 + e_4 * fs_375_286_110 * r_8 * h2_p2 - e_5 * fs_7_4199_2 * r_2 * h10_p2 + e_5 * fs_14_4199_62985 * r_2 * h10_p10 + e_5 * fs_4_2717_231 * r_4 * h8_p2 - e_5 * fs_2_187_77 * r_6 * h6_p2 + e_5 * fs_4_143_66 * r_8 * h4_p2 - e_5 * fs_5_143_110 * r_10 * h2_p2;

        pc_1[k] = - e_1 * fs_255_8_462 * h4_p3 - e_2 * fs_855_88_154 * h6_p3 + e_2 * fs_255_11_462 * r_2 * h4_p3 - e_3 * fs_63_26_7 * h8_p3 + e_3 * fs_171_44_154 * r_2 * h6_p3 - e_3 * fs_765_143_462 * r_4 * h4_p3 - e_4 * fs_161_8398_26 * h10_p3 - e_4 * fs_161_4199_25194 * h10_p9 + e_4 * fs_126_247_7 * r_2 * h8_p3 - e_4 * fs_171_374_154 * r_4 * h6_p3 + e_4 * fs_68_143_462 * r_6 * h4_p3 + e_5 * fs_7_4199_26 * r_2 * h10_p3 + e_5 * fs_14_4199_25194 * r_2 * h10_p9 - e_5 * fs_6_247_7 * r_4 * h8_p3 + e_5 * fs_3_187_154 * r_6 * h6_p3 - e_5 * fs_2_143_462 * r_8 * h4_p3;
    }

    // NOTE: the angular components are formed in 93 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph4_p4, ph6_p4, ph6_p5, ph8_p4, ph8_p5, ph8_p7, ph8_p8, ph10_p4, ph10_p5, ph10_p7, ph10_p8, ab_2, pc_2, pc_3 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

        const auto h4_p4 = ph4_p4[k];
        const auto h6_p4 = ph6_p4[k];
        const auto h6_p5 = ph6_p5[k];
        const auto h8_p4 = ph8_p4[k];
        const auto h8_p5 = ph8_p5[k];
        const auto h8_p7 = ph8_p7[k];
        const auto h8_p8 = ph8_p8[k];
        const auto h10_p4 = ph10_p4[k];
        const auto h10_p5 = ph10_p5[k];
        const auto h10_p7 = ph10_p7[k];
        const auto h10_p8 = ph10_p8[k];

        pc_2[k] = e_1 * fs_255_4_66 * h4_p4 + e_2 * fs_855_88_330 * h6_p4 - e_2 * fs_510_11_66 * r_2 * h4_p4 + e_3 * fs_63_26_30 * h8_p4 - e_3 * fs_21_13_546 * h8_p8 - e_3 * fs_171_44_330 * r_2 * h6_p4 + e_3 * fs_1530_143_66 * r_4 * h4_p4 + e_4 * fs_161_8398_182 * h10_p4 - e_4 * fs_161_4199_9282 * h10_p8 - e_4 * fs_126_247_30 * r_2 * h8_p4 + e_4 * fs_84_247_546 * r_2 * h8_p8 + e_4 * fs_171_374_330 * r_4 * h6_p4 - e_4 * fs_136_143_66 * r_6 * h4_p4 - e_5 * fs_7_4199_182 * r_2 * h10_p4 + e_5 * fs_14_4199_9282 * r_2 * h10_p8 + e_5 * fs_6_247_30 * r_4 * h8_p4 - e_5 * fs_4_247_546 * r_4 * h8_p8 - e_5 * fs_3_187_330 * r_6 * h6_p4 + e_5 * fs_4_143_66 * r_8 * h4_p4;

        pc_3[k] = - e_2 * fs_285_8_30 * h6_p5 - e_3 * fs_21_13_195 * h8_p5 - e_3 * fs_63_26_273 * h8_p7 + e_3 * fs_57_4_30 * r_2 * h6_p5 - e_4 * fs_161_8398_910 * h10_p5 - e_4 * fs_161_4199_3094 * h10_p7 + e_4 * fs_84_247_195 * r_2 * h8_p5 + e_4 * fs_126_247_273 * r_2 * h8_p7 - e_4 * fs_57_34_30 * r_4 * h6_p5 + e_5 * fs_7_4199_910 * r_2 * h10_p5 + e_5 * fs_14_4199_3094 * r_2 * h10_p7 - e_5 * fs_4_247_195 * r_4 * h8_p5 - e_5 * fs_6_247_273 * r_4 * h8_p7 + e_5 * fs_1_17_30 * r_6 * h6_p5;
    }

    // NOTE: the angular components are formed in 93 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_2, pe_3, pe_4, pe_5, ph6_m6, ph6_m5, ph8_m7, ph8_m6, ph8_m5, ph10_m7, ph10_m6, ph10_m5, ab_2, pc_4, pc_5 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h6_m6 = ph6_m6[k];
        const auto h6_m5 = ph6_m5[k];
        const auto h8_m7 = ph8_m7[k];
        const auto h8_m6 = ph8_m6[k];
        const auto h8_m5 = ph8_m5[k];
        const auto h10_m7 = ph10_m7[k];
        const auto h10_m6 = ph10_m6[k];
        const auto h10_m5 = ph10_m5[k];

        pc_4[k] = e_2 * f_855_4 * h6_m6 + e_3 * fs_63_13_91 * h8_m6 - e_3 * f_171_2 * r_2 * h6_m6 + e_4 * fs_322_4199_455 * h10_m6 - e_4 * fs_252_247_91 * r_2 * h8_m6 + e_4 * f_171_17 * r_4 * h6_m6 - e_5 * fs_28_4199_455 * r_2 * h10_m6 + e_5 * fs_12_247_91 * r_4 * h8_m6 - e_5 * f_6_17 * r_6 * h6_m6;

        pc_5[k] = - e_2 * fs_285_8_30 * h6_m5 + e_3 * fs_63_26_273 * h8_m7 - e_3 * fs_21_13_195 * h8_m5 + e_3 * fs_57_4_30 * r_2 * h6_m5 + e_4 * fs_161_4199_3094 * h10_m7 - e_4 * fs_161_8398_910 * h10_m5 - e_4 * fs_126_247_273 * r_2 * h8_m7 + e_4 * fs_84_247_195 * r_2 * h8_m5 - e_4 * fs_57_34_30 * r_4 * h6_m5 - e_5 * fs_14_4199_3094 * r_2 * h10_m7 + e_5 * fs_7_4199_910 * r_2 * h10_m5 + e_5 * fs_6_247_273 * r_4 * h8_m7 - e_5 * fs_4_247_195 * r_4 * h8_m5 + e_5 * fs_1_17_30 * r_6 * h6_m5;
    }

    // NOTE: the angular components are formed in 93 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph4_m4, ph4_m3, ph6_m4, ph6_m3, ph8_m8, ph8_m4, ph8_m3, ph10_m9, ph10_m8, ph10_m4, ph10_m3, ab_2, pc_6, pc_7 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

        const auto h4_m4 = ph4_m4[k];
        const auto h4_m3 = ph4_m3[k];
        const auto h6_m4 = ph6_m4[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h8_m8 = ph8_m8[k];
        const auto h8_m4 = ph8_m4[k];
        const auto h8_m3 = ph8_m3[k];
        const auto h10_m9 = ph10_m9[k];
        const auto h10_m8 = ph10_m8[k];
        const auto h10_m4 = ph10_m4[k];
        const auto h10_m3 = ph10_m3[k];

        pc_6[k] = e_1 * fs_255_4_66 * h4_m4 + e_2 * fs_855_88_330 * h6_m4 - e_2 * fs_510_11_66 * r_2 * h4_m4 + e_3 * fs_21_13_546 * h8_m8 + e_3 * fs_63_26_30 * h8_m4 - e_3 * fs_171_44_330 * r_2 * h6_m4 + e_3 * fs_1530_143_66 * r_4 * h4_m4 + e_4 * fs_161_4199_9282 * h10_m8 + e_4 * fs_161_8398_182 * h10_m4 - e_4 * fs_84_247_546 * r_2 * h8_m8 - e_4 * fs_126_247_30 * r_2 * h8_m4 + e_4 * fs_171_374_330 * r_4 * h6_m4 - e_4 * fs_136_143_66 * r_6 * h4_m4 - e_5 * fs_14_4199_9282 * r_2 * h10_m8 - e_5 * fs_7_4199_182 * r_2 * h10_m4 + e_5 * fs_4_247_546 * r_4 * h8_m8 + e_5 * fs_6_247_30 * r_4 * h8_m4 - e_5 * fs_3_187_330 * r_6 * h6_m4 + e_5 * fs_4_143_66 * r_8 * h4_m4;

        pc_7[k] = - e_1 * fs_255_8_462 * h4_m3 - e_2 * fs_855_88_154 * h6_m3 + e_2 * fs_255_11_462 * r_2 * h4_m3 - e_3 * fs_63_26_7 * h8_m3 + e_3 * fs_171_44_154 * r_2 * h6_m3 - e_3 * fs_765_143_462 * r_4 * h4_m3 + e_4 * fs_161_4199_25194 * h10_m9 - e_4 * fs_161_8398_26 * h10_m3 + e_4 * fs_126_247_7 * r_2 * h8_m3 - e_4 * fs_171_374_154 * r_4 * h6_m3 + e_4 * fs_68_143_462 * r_6 * h4_m3 - e_5 * fs_14_4199_25194 * r_2 * h10_m9 + e_5 * fs_7_4199_26 * r_2 * h10_m3 - e_5 * fs_6_247_7 * r_4 * h8_m3 + e_5 * fs_3_187_154 * r_6 * h6_m3 - e_5 * fs_2_143_462 * r_8 * h4_m3;
    }

    // NOTE: the angular components are formed in 93 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph2_m2, ph4_m2, ph6_m2, ph8_m2, ph10_m10, ph10_m2, ab_2, pc_8 : simd::cache_line_size())
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
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;

        const auto h2_m2 = ph2_m2[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h8_m2 = ph8_m2[k];
        const auto h10_m10 = ph10_m10[k];
        const auto h10_m2 = ph10_m2[k];

        pc_8[k] = e_0 * fs_4725_32_110 * h2_m2 + e_1 * fs_255_4_66 * h4_m2 - e_1 * fs_3375_16_110 * r_2 * h2_m2 + e_2 * fs_285_44_77 * h6_m2 - e_2 * fs_510_11_66 * r_2 * h4_m2 + e_2 * fs_375_4_110 * r_4 * h2_m2 + e_3 * fs_21_143_231 * h8_m2 - e_3 * fs_57_22_77 * r_2 * h6_m2 + e_3 * fs_1530_143_66 * r_4 * h4_m2 - e_3 * fs_375_22_110 * r_6 * h2_m2 + e_4 * fs_161_4199_62985 * h10_m10 + e_4 * fs_161_8398_2 * h10_m2 - e_4 * fs_84_2717_231 * r_2 * h8_m2 + e_4 * fs_57_187_77 * r_4 * h6_m2 - e_4 * fs_136_143_66 * r_6 * h4_m2 + e_4 * fs_375_286_110 * r_8 * h2_m2 - e_5 * fs_14_4199_62985 * r_2 * h10_m10 - e_5 * fs_7_4199_2 * r_2 * h10_m2 + e_5 * fs_4_2717_231 * r_4 * h8_m2 - e_5 * fs_2_187_77 * r_6 * h6_m2 + e_5 * fs_4_143_66 * r_8 * h4_m2 - e_5 * fs_5_143_110 * r_10 * h2_m2;
    }

    // NOTE: the angular components are formed in 93 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph2_p1, ph4_p1, ph6_p1, ph8_p1, ph10_p1, ph10_p9, ab_2, pc_9 : simd::cache_line_size())
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
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;

        const auto h2_p1 = ph2_p1[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h10_p1 = ph10_p1[k];
        const auto h10_p9 = ph10_p9[k];

        pc_9[k] = e_0 * fs_1575_32_330 * h2_p1 + e_1 * fs_765_4_11 * h4_p1 - e_1 * fs_1125_16_330 * r_2 * h2_p1 + e_2 * fs_95_44_2310 * h6_p1 - e_2 * fs_1530_11_11 * r_2 * h4_p1 + e_2 * fs_125_4_330 * r_4 * h2_p1 + e_3 * fs_147_286_110 * h8_p1 - e_3 * fs_19_22_2310 * r_2 * h6_p1 + e_3 * fs_4590_143_11 * r_4 * h4_p1 - e_3 * fs_125_22_330 * r_6 * h2_p1 + e_4 * fs_483_8398_2 * h10_p1 - e_4 * fs_483_4199_4199 * h10_p9 - e_4 * fs_294_2717_110 * r_2 * h8_p1 + e_4 * fs_19_187_2310 * r_4 * h6_p1 - e_4 * fs_408_143_11 * r_6 * h4_p1 + e_4 * fs_125_286_330 * r_8 * h2_p1 - e_5 * fs_21_4199_2 * r_2 * h10_p1 + e_5 * fs_42_4199_4199 * r_2 * h10_p9 + e_5 * fs_14_2717_110 * r_4 * h8_p1 - e_5 * fs_2_561_2310 * r_6 * h6_p1 + e_5 * fs_12_143_11 * r_8 * h4_p1 - e_5 * fs_5_429_330 * r_10 * h2_p1;
    }

    // NOTE: the angular components are formed in 93 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph2_p2, ph4_p2, ph6_p2, ph8_p2, ph8_p8, ph10_p2, ph10_p8, ab_2, pc_10 : simd::cache_line_size())
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
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;

        const auto h2_p2 = ph2_p2[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h8_p2 = ph8_p2[k];
        const auto h8_p8 = ph8_p8[k];
        const auto h10_p2 = ph10_p2[k];
        const auto h10_p8 = ph10_p8[k];

        pc_10[k] = e_0 * fs_1575_16_165 * h2_p2 - e_1 * fs_765_8_11 * h4_p2 - e_1 * fs_1125_8_165 * r_2 * h2_p2 - e_2 * fs_665_88_462 * h6_p2 + e_2 * fs_765_11_11 * r_2 * h4_p2 + e_2 * fs_125_2_165 * r_4 * h2_p2 - e_3 * fs_609_572_154 * h8_p2 + e_3 * fs_147_13_13 * h8_p8 + e_3 * fs_133_44_462 * r_2 * h6_p2 - e_3 * fs_2295_143_11 * r_4 * h4_p2 - e_3 * fs_125_11_165 * r_6 * h2_p2 - e_4 * fs_644_4199_3 * h10_p2 - e_4 * fs_1932_4199_221 * h10_p8 + e_4 * fs_609_2717_154 * r_2 * h8_p2 - e_4 * fs_588_247_13 * r_2 * h8_p8 - e_4 * fs_133_374_462 * r_4 * h6_p2 + e_4 * fs_204_143_11 * r_6 * h4_p2 + e_4 * fs_125_143_165 * r_8 * h2_p2 + e_5 * fs_56_4199_3 * r_2 * h10_p2 + e_5 * fs_168_4199_221 * r_2 * h10_p8 - e_5 * fs_29_2717_154 * r_4 * h8_p2 + e_5 * fs_28_247_13 * r_4 * h8_p8 + e_5 * fs_7_561_462 * r_6 * h6_p2 - e_5 * fs_6_143_11 * r_8 * h4_p2 - e_5 * fs_10_429_165 * r_10 * h2_p2;
    }

    // NOTE: the angular components are formed in 93 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph4_p3, ph6_p3, ph8_p3, ph8_p7, ph10_p3, ph10_p7, ab_2, pc_11 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

        const auto h4_p3 = ph4_p3[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h8_p3 = ph8_p3[k];
        const auto h8_p7 = ph8_p7[k];
        const auto h10_p3 = ph10_p3[k];
        const auto h10_p7 = ph10_p7[k];

        pc_11[k] = - e_1 * fs_765_8_11 * h4_p3 + e_2 * fs_285_11_33 * h6_p3 + e_2 * fs_765_11_11 * r_2 * h4_p3 + e_3 * fs_483_52_6 * h8_p3 + e_3 * fs_21_52_182 * h8_p7 - e_3 * fs_114_11_33 * r_2 * h6_p3 - e_3 * fs_2295_143_11 * r_4 * h4_p3 + e_4 * fs_161_4199_273 * h10_p3 - e_4 * fs_322_4199_4641 * h10_p7 - e_4 * fs_483_247_6 * r_2 * h8_p3 - e_4 * fs_21_247_182 * r_2 * h8_p7 + e_4 * fs_228_187_33 * r_4 * h6_p3 + e_4 * fs_204_143_11 * r_6 * h4_p3 - e_5 * fs_14_4199_273 * r_2 * h10_p3 + e_5 * fs_28_4199_4641 * r_2 * h10_p7 + e_5 * fs_23_247_6 * r_4 * h8_p3 + e_5 * fs_1_247_182 * r_4 * h8_p7 - e_5 * fs_8_187_33 * r_6 * h6_p3 - e_5 * fs_6_143_11 * r_8 * h4_p3;
    }

    // NOTE: the angular components are formed in 93 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph4_p4, ph6_m5, ph6_p4, ph6_p6, ph8_m5, ph8_p4, ph8_p6, ph10_m5, ph10_p4, ph10_p6, ab_2, pc_12, pc_13 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

        const auto h4_p4 = ph4_p4[k];
        const auto h6_m5 = ph6_m5[k];
        const auto h6_p4 = ph6_p4[k];
        const auto h6_p6 = ph6_p6[k];
        const auto h8_m5 = ph8_m5[k];
        const auto h8_p4 = ph8_p4[k];
        const auto h8_p6 = ph8_p6[k];
        const auto h10_m5 = ph10_m5[k];
        const auto h10_p4 = ph10_p4[k];
        const auto h10_p6 = ph10_p6[k];

        pc_12[k] = e_1 * fs_765_4_11 * h4_p4 - e_2 * fs_285_44_55 * h6_p4 + e_2 * fs_285_8_30 * h6_p6 - e_2 * fs_1530_11_11 * r_2 * h4_p4 - e_3 * fs_357_26_5 * h8_p4 - e_3 * fs_21_52_2730 * h8_p6 + e_3 * fs_57_22_55 * r_2 * h6_p4 - e_3 * fs_57_4_30 * r_2 * h6_p6 + e_3 * fs_4590_143_11 * r_4 * h4_p4 - e_4 * fs_322_4199_273 * h10_p4 - e_4 * fs_644_4199_546 * h10_p6 + e_4 * fs_714_247_5 * r_2 * h8_p4 + e_4 * fs_21_247_2730 * r_2 * h8_p6 - e_4 * fs_57_187_55 * r_4 * h6_p4 + e_4 * fs_57_34_30 * r_4 * h6_p6 - e_4 * fs_408_143_11 * r_6 * h4_p4 + e_5 * fs_28_4199_273 * r_2 * h10_p4 + e_5 * fs_56_4199_546 * r_2 * h10_p6 - e_5 * fs_34_247_5 * r_4 * h8_p4 - e_5 * fs_1_247_2730 * r_4 * h8_p6 + e_5 * fs_2_187_55 * r_6 * h6_p4 - e_5 * fs_1_17_30 * r_6 * h6_p6 + e_5 * fs_12_143_11 * r_8 * h4_p4;

        pc_13[k] = - e_2 * f_285_2 * h6_m5 + e_3 * fs_231_26_26 * h8_m5 + e_3 * f_57_1 * r_2 * h6_m5 + e_4 * fs_805_4199_273 * h10_m5 - e_4 * fs_462_247_26 * r_2 * h8_m5 - e_4 * f_114_17 * r_4 * h6_m5 - e_5 * fs_70_4199_273 * r_2 * h10_m5 + e_5 * fs_22_247_26 * r_4 * h8_m5 + e_5 * f_4_17 * r_6 * h6_m5;
    }

    // NOTE: the angular components are formed in 93 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph4_m4, ph6_m6, ph6_m4, ph8_m6, ph8_m4, ph10_m6, ph10_m4, ab_2, pc_14 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

        const auto h4_m4 = ph4_m4[k];
        const auto h6_m6 = ph6_m6[k];
        const auto h6_m4 = ph6_m4[k];
        const auto h8_m6 = ph8_m6[k];
        const auto h8_m4 = ph8_m4[k];
        const auto h10_m6 = ph10_m6[k];
        const auto h10_m4 = ph10_m4[k];

        pc_14[k] = e_1 * fs_765_4_11 * h4_m4 - e_2 * fs_285_8_30 * h6_m6 - e_2 * fs_285_44_55 * h6_m4 - e_2 * fs_1530_11_11 * r_2 * h4_m4 + e_3 * fs_21_52_2730 * h8_m6 - e_3 * fs_357_26_5 * h8_m4 + e_3 * fs_57_4_30 * r_2 * h6_m6 + e_3 * fs_57_22_55 * r_2 * h6_m4 + e_3 * fs_4590_143_11 * r_4 * h4_m4 + e_4 * fs_644_4199_546 * h10_m6 - e_4 * fs_322_4199_273 * h10_m4 - e_4 * fs_21_247_2730 * r_2 * h8_m6 + e_4 * fs_714_247_5 * r_2 * h8_m4 - e_4 * fs_57_34_30 * r_4 * h6_m6 - e_4 * fs_57_187_55 * r_4 * h6_m4 - e_4 * fs_408_143_11 * r_6 * h4_m4 - e_5 * fs_56_4199_546 * r_2 * h10_m6 + e_5 * fs_28_4199_273 * r_2 * h10_m4 + e_5 * fs_1_247_2730 * r_4 * h8_m6 - e_5 * fs_34_247_5 * r_4 * h8_m4 + e_5 * fs_1_17_30 * r_6 * h6_m6 + e_5 * fs_2_187_55 * r_6 * h6_m4 + e_5 * fs_12_143_11 * r_8 * h4_m4;
    }

    // NOTE: the angular components are formed in 93 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph4_m3, ph6_m3, ph8_m7, ph8_m3, ph10_m7, ph10_m3, ab_2, pc_15 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

        const auto h4_m3 = ph4_m3[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h8_m7 = ph8_m7[k];
        const auto h8_m3 = ph8_m3[k];
        const auto h10_m7 = ph10_m7[k];
        const auto h10_m3 = ph10_m3[k];

        pc_15[k] = - e_1 * fs_765_8_11 * h4_m3 + e_2 * fs_285_11_33 * h6_m3 + e_2 * fs_765_11_11 * r_2 * h4_m3 - e_3 * fs_21_52_182 * h8_m7 + e_3 * fs_483_52_6 * h8_m3 - e_3 * fs_114_11_33 * r_2 * h6_m3 - e_3 * fs_2295_143_11 * r_4 * h4_m3 + e_4 * fs_322_4199_4641 * h10_m7 + e_4 * fs_161_4199_273 * h10_m3 + e_4 * fs_21_247_182 * r_2 * h8_m7 - e_4 * fs_483_247_6 * r_2 * h8_m3 + e_4 * fs_228_187_33 * r_4 * h6_m3 + e_4 * fs_204_143_11 * r_6 * h4_m3 - e_5 * fs_28_4199_4641 * r_2 * h10_m7 - e_5 * fs_14_4199_273 * r_2 * h10_m3 - e_5 * fs_1_247_182 * r_4 * h8_m7 + e_5 * fs_23_247_6 * r_4 * h8_m3 - e_5 * fs_8_187_33 * r_6 * h6_m3 - e_5 * fs_6_143_11 * r_8 * h4_m3;
    }

    // NOTE: the angular components are formed in 93 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph2_m2, ph4_m2, ph6_m2, ph8_m8, ph8_m2, ph10_m8, ph10_m2, ab_2, pc_16 : simd::cache_line_size())
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
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;

        const auto h2_m2 = ph2_m2[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h8_m8 = ph8_m8[k];
        const auto h8_m2 = ph8_m2[k];
        const auto h10_m8 = ph10_m8[k];
        const auto h10_m2 = ph10_m2[k];

        pc_16[k] = e_0 * fs_1575_16_165 * h2_m2 - e_1 * fs_765_8_11 * h4_m2 - e_1 * fs_1125_8_165 * r_2 * h2_m2 - e_2 * fs_665_88_462 * h6_m2 + e_2 * fs_765_11_11 * r_2 * h4_m2 + e_2 * fs_125_2_165 * r_4 * h2_m2 - e_3 * fs_147_13_13 * h8_m8 - e_3 * fs_609_572_154 * h8_m2 + e_3 * fs_133_44_462 * r_2 * h6_m2 - e_3 * fs_2295_143_11 * r_4 * h4_m2 - e_3 * fs_125_11_165 * r_6 * h2_m2 + e_4 * fs_1932_4199_221 * h10_m8 - e_4 * fs_644_4199_3 * h10_m2 + e_4 * fs_588_247_13 * r_2 * h8_m8 + e_4 * fs_609_2717_154 * r_2 * h8_m2 - e_4 * fs_133_374_462 * r_4 * h6_m2 + e_4 * fs_204_143_11 * r_6 * h4_m2 + e_4 * fs_125_143_165 * r_8 * h2_m2 - e_5 * fs_168_4199_221 * r_2 * h10_m8 + e_5 * fs_56_4199_3 * r_2 * h10_m2 - e_5 * fs_28_247_13 * r_4 * h8_m8 - e_5 * fs_29_2717_154 * r_4 * h8_m2 + e_5 * fs_7_561_462 * r_6 * h6_m2 - e_5 * fs_6_143_11 * r_8 * h4_m2 - e_5 * fs_10_429_165 * r_10 * h2_m2;
    }

    // NOTE: the angular components are formed in 93 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph2_m1, ph4_m1, ph6_m1, ph8_m1, ph10_m9, ph10_m1, ab_2, pc_17 : simd::cache_line_size())
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
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;

        const auto h2_m1 = ph2_m1[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h8_m1 = ph8_m1[k];
        const auto h10_m9 = ph10_m9[k];
        const auto h10_m1 = ph10_m1[k];

        pc_17[k] = e_0 * fs_1575_32_330 * h2_m1 + e_1 * fs_765_4_11 * h4_m1 - e_1 * fs_1125_16_330 * r_2 * h2_m1 + e_2 * fs_95_44_2310 * h6_m1 - e_2 * fs_1530_11_11 * r_2 * h4_m1 + e_2 * fs_125_4_330 * r_4 * h2_m1 + e_3 * fs_147_286_110 * h8_m1 - e_3 * fs_19_22_2310 * r_2 * h6_m1 + e_3 * fs_4590_143_11 * r_4 * h4_m1 - e_3 * fs_125_22_330 * r_6 * h2_m1 + e_4 * fs_483_4199_4199 * h10_m9 + e_4 * fs_483_8398_2 * h10_m1 - e_4 * fs_294_2717_110 * r_2 * h8_m1 + e_4 * fs_19_187_2310 * r_4 * h6_m1 - e_4 * fs_408_143_11 * r_6 * h4_m1 + e_4 * fs_125_286_330 * r_8 * h2_m1 - e_5 * fs_42_4199_4199 * r_2 * h10_m9 - e_5 * fs_21_4199_2 * r_2 * h10_m1 + e_5 * fs_14_2717_110 * r_4 * h8_m1 - e_5 * fs_2_561_2310 * r_6 * h6_m1 + e_5 * fs_12_143_11 * r_8 * h4_m1 - e_5 * fs_5_429_330 * r_10 * h2_m1;
    }

    // NOTE: the angular components are formed in 93 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph2_0, ph4_0, ph6_0, ph8_0, ph8_p8, ph10_0, ph10_p8, ab_2, pc_18 : simd::cache_line_size())
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
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;

        const auto h2_0 = ph2_0[k];
        const auto h4_0 = ph4_0[k];
        const auto h6_0 = ph6_0[k];
        const auto h8_0 = ph8_0[k];
        const auto h8_p8 = ph8_p8[k];
        const auto h10_0 = ph10_0[k];
        const auto h10_p8 = ph10_p8[k];

        pc_18[k] = e_0 * fs_4725_16_5 * h2_0 + e_1 * fs_765_2_5 * h4_0 - e_1 * fs_3375_8_5 * r_2 * h2_0 + e_2 * fs_1995_22_5 * h6_0 - e_2 * fs_3060_11_5 * r_2 * h4_0 + e_2 * fs_375_2_5 * r_4 * h2_0 + e_3 * fs_882_143_5 * h8_0 - e_3 * fs_294_143_143 * h8_p8 - e_3 * fs_399_11_5 * r_2 * h6_0 + e_3 * fs_9180_143_5 * r_4 * h4_0 - e_3 * fs_375_11_5 * r_6 * h2_0 + e_4 * fs_483_4199_5 * h10_0 - e_4 * fs_483_4199_2431 * h10_p8 - e_4 * fs_3528_2717_5 * r_2 * h8_0 + e_4 * fs_1176_2717_143 * r_2 * h8_p8 + e_4 * fs_798_187_5 * r_4 * h6_0 - e_4 * fs_816_143_5 * r_6 * h4_0 + e_4 * fs_375_143_5 * r_8 * h2_0 - e_5 * fs_42_4199_5 * r_2 * h10_0 + e_5 * fs_42_4199_2431 * r_2 * h10_p8 + e_5 * fs_168_2717_5 * r_4 * h8_0 - e_5 * fs_56_2717_143 * r_4 * h8_p8 - e_5 * fs_28_187_5 * r_6 * h6_0 + e_5 * fs_24_143_5 * r_8 * h4_0 - e_5 * fs_10_143_5 * r_10 * h2_0;
    }

    // NOTE: the angular components are formed in 93 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph2_p1, ph4_p1, ph6_p1, ph8_p1, ph8_p7, ph10_p1, ph10_p7, ab_2, pc_19 : simd::cache_line_size())
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
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;

        const auto h2_p1 = ph2_p1[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h8_p7 = ph8_p7[k];
        const auto h10_p1 = ph10_p1[k];
        const auto h10_p7 = ph10_p7[k];

        pc_19[k] = e_0 * fs_1575_8_30 * h2_p1 + e_1 * f_765_8 * h4_p1 - e_1 * fs_1125_4_30 * r_2 * h2_p1 - e_2 * fs_475_44_210 * h6_p1 - e_2 * f_765_11 * r_2 * h4_p1 + e_2 * fs_125_1_30 * r_4 * h2_p1 - e_3 * fs_3675_572_10 * h8_p1 + e_3 * fs_1029_572_286 * h8_p7 + e_3 * fs_95_22_210 * r_2 * h6_p1 + e_3 * f_2295_143 * r_4 * h4_p1 - e_3 * fs_250_11_30 * r_6 * h2_p1 - e_4 * fs_483_4199_22 * h10_p1 - e_4 * fs_322_4199_7293 * h10_p7 + e_4 * fs_3675_2717_10 * r_2 * h8_p1 - e_4 * fs_1029_2717_286 * r_2 * h8_p7 - e_4 * fs_95_187_210 * r_4 * h6_p1 - e_4 * f_204_143 * r_6 * h4_p1 + e_4 * fs_250_143_30 * r_8 * h2_p1 + e_5 * fs_42_4199_22 * r_2 * h10_p1 + e_5 * fs_28_4199_7293 * r_2 * h10_p7 - e_5 * fs_175_2717_10 * r_4 * h8_p1 + e_5 * fs_49_2717_286 * r_4 * h8_p7 + e_5 * fs_10_561_210 * r_6 * h6_p1 + e_5 * f_6_143 * r_8 * h4_p1 - e_5 * fs_20_429_30 * r_10 * h2_p1;
    }

    // NOTE: the angular components are formed in 93 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph2_p2, ph4_p2, ph6_p2, ph6_p6, ph8_p2, ph8_p6, ph10_p2, ph10_p6, ab_2, pc_20 : simd::cache_line_size())
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
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;

        const auto h2_p2 = ph2_p2[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h6_p6 = ph6_p6[k];
        const auto h8_p2 = ph8_p2[k];
        const auto h8_p6 = ph8_p6[k];
        const auto h10_p2 = ph10_p2[k];
        const auto h10_p6 = ph10_p6[k];

        pc_20[k] = e_0 * fs_1575_16_105 * h2_p2 - e_1 * fs_765_4_7 * h4_p2 - e_1 * fs_1125_8_105 * r_2 * h2_p2 + e_2 * fs_2185_88_6 * h6_p2 - e_2 * fs_855_88_330 * h6_p6 + e_2 * fs_1530_11_7 * r_2 * h4_p2 + e_2 * fs_125_2_105 * r_4 * h2_p2 + e_3 * fs_2919_143_2 * h8_p2 + e_3 * fs_21_143_30030 * h8_p6 - e_3 * fs_437_44_6 * r_2 * h6_p2 + e_3 * fs_171_44_330 * r_2 * h6_p6 - e_3 * fs_4590_143_7 * r_4 * h4_p2 - e_3 * fs_125_11_105 * r_6 * h2_p2 + e_4 * fs_322_4199_231 * h10_p2 - e_4 * fs_322_4199_6006 * h10_p6 - e_4 * fs_11676_2717_2 * r_2 * h8_p2 - e_4 * fs_84_2717_30030 * r_2 * h8_p6 + e_4 * fs_437_374_6 * r_4 * h6_p2 - e_4 * fs_171_374_330 * r_4 * h6_p6 + e_4 * fs_408_143_7 * r_6 * h4_p2 + e_4 * fs_125_143_105 * r_8 * h2_p2 - e_5 * fs_28_4199_231 * r_2 * h10_p2 + e_5 * fs_28_4199_6006 * r_2 * h10_p6 + e_5 * fs_556_2717_2 * r_4 * h8_p2 + e_5 * fs_4_2717_30030 * r_4 * h8_p6 - e_5 * fs_23_561_6 * r_6 * h6_p2 + e_5 * fs_3_187_330 * r_6 * h6_p6 - e_5 * fs_12_143_7 * r_8 * h4_p2 - e_5 * fs_10_429_105 * r_10 * h2_p2;
    }

    // NOTE: the angular components are formed in 93 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph4_m4, ph4_p3, ph6_m4, ph6_p3, ph6_p5, ph8_m4, ph8_p3, ph8_p5, ph10_m4, ph10_p3, ph10_p5, ab_2, pc_21, pc_22 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

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

        pc_21[k] = e_1 * f_765_8 * h4_p3 + e_2 * fs_1995_44_3 * h6_p3 + e_2 * fs_285_44_55 * h6_p5 - e_2 * f_765_11 * r_2 * h4_p3 - e_3 * fs_2121_572_66 * h8_p3 + e_3 * fs_21_572_1430 * h8_p5 - e_3 * fs_399_22_3 * r_2 * h6_p3 - e_3 * fs_57_22_55 * r_2 * h6_p5 + e_3 * f_2295_143 * r_4 * h4_p3 - e_4 * fs_161_4199_3003 * h10_p3 - e_4 * fs_161_4199_15015 * h10_p5 + e_4 * fs_2121_2717_66 * r_2 * h8_p3 - e_4 * fs_21_2717_1430 * r_2 * h8_p5 + e_4 * fs_399_187_3 * r_4 * h6_p3 + e_4 * fs_57_187_55 * r_4 * h6_p5 - e_4 * f_204_143 * r_6 * h4_p3 + e_5 * fs_14_4199_3003 * r_2 * h10_p3 + e_5 * fs_14_4199_15015 * r_2 * h10_p5 - e_5 * fs_101_2717_66 * r_4 * h8_p3 + e_5 * fs_1_2717_1430 * r_4 * h8_p5 - e_5 * fs_14_187_3 * r_6 * h6_p3 - e_5 * fs_2_187_55 * r_6 * h6_p5 + e_5 * f_6_143 * r_8 * h4_p3;

        pc_22[k] = e_1 * fs_765_2_5 * h4_m4 - e_2 * f_2280_11 * h6_m4 - e_2 * fs_3060_11_5 * r_2 * h4_m4 + e_3 * fs_1218_143_11 * h8_m4 + e_3 * f_912_11 * r_2 * h6_m4 + e_3 * fs_9180_143_5 * r_4 * h4_m4 + e_4 * fs_161_4199_15015 * h10_m4 - e_4 * fs_4872_2717_11 * r_2 * h8_m4 - e_4 * f_1824_187 * r_4 * h6_m4 - e_4 * fs_816_143_5 * r_6 * h4_m4 - e_5 * fs_14_4199_15015 * r_2 * h10_m4 + e_5 * fs_232_2717_11 * r_4 * h8_m4 + e_5 * f_64_187 * r_6 * h6_m4 + e_5 * fs_24_143_5 * r_8 * h4_m4;
    }

    // NOTE: the angular components are formed in 93 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph4_m3, ph6_m5, ph6_m3, ph8_m5, ph8_m3, ph10_m5, ph10_m3, ab_2, pc_23 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

        const auto h4_m3 = ph4_m3[k];
        const auto h6_m5 = ph6_m5[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h8_m5 = ph8_m5[k];
        const auto h8_m3 = ph8_m3[k];
        const auto h10_m5 = ph10_m5[k];
        const auto h10_m3 = ph10_m3[k];

        pc_23[k] = e_1 * f_765_8 * h4_m3 - e_2 * fs_285_44_55 * h6_m5 + e_2 * fs_1995_44_3 * h6_m3 - e_2 * f_765_11 * r_2 * h4_m3 - e_3 * fs_21_572_1430 * h8_m5 - e_3 * fs_2121_572_66 * h8_m3 + e_3 * fs_57_22_55 * r_2 * h6_m5 - e_3 * fs_399_22_3 * r_2 * h6_m3 + e_3 * f_2295_143 * r_4 * h4_m3 + e_4 * fs_161_4199_15015 * h10_m5 - e_4 * fs_161_4199_3003 * h10_m3 + e_4 * fs_21_2717_1430 * r_2 * h8_m5 + e_4 * fs_2121_2717_66 * r_2 * h8_m3 - e_4 * fs_57_187_55 * r_4 * h6_m5 + e_4 * fs_399_187_3 * r_4 * h6_m3 - e_4 * f_204_143 * r_6 * h4_m3 - e_5 * fs_14_4199_15015 * r_2 * h10_m5 + e_5 * fs_14_4199_3003 * r_2 * h10_m3 - e_5 * fs_1_2717_1430 * r_4 * h8_m5 - e_5 * fs_101_2717_66 * r_4 * h8_m3 + e_5 * fs_2_187_55 * r_6 * h6_m5 - e_5 * fs_14_187_3 * r_6 * h6_m3 + e_5 * f_6_143 * r_8 * h4_m3;
    }

    // NOTE: the angular components are formed in 93 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph2_m2, ph4_m2, ph6_m6, ph6_m2, ph8_m6, ph8_m2, ph10_m6, ph10_m2, ab_2, pc_24 : simd::cache_line_size())
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
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;

        const auto h2_m2 = ph2_m2[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h6_m6 = ph6_m6[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h8_m6 = ph8_m6[k];
        const auto h8_m2 = ph8_m2[k];
        const auto h10_m6 = ph10_m6[k];
        const auto h10_m2 = ph10_m2[k];

        pc_24[k] = e_0 * fs_1575_16_105 * h2_m2 - e_1 * fs_765_4_7 * h4_m2 - e_1 * fs_1125_8_105 * r_2 * h2_m2 + e_2 * fs_855_88_330 * h6_m6 + e_2 * fs_2185_88_6 * h6_m2 + e_2 * fs_1530_11_7 * r_2 * h4_m2 + e_2 * fs_125_2_105 * r_4 * h2_m2 - e_3 * fs_21_143_30030 * h8_m6 + e_3 * fs_2919_143_2 * h8_m2 - e_3 * fs_171_44_330 * r_2 * h6_m6 - e_3 * fs_437_44_6 * r_2 * h6_m2 - e_3 * fs_4590_143_7 * r_4 * h4_m2 - e_3 * fs_125_11_105 * r_6 * h2_m2 + e_4 * fs_322_4199_6006 * h10_m6 + e_4 * fs_322_4199_231 * h10_m2 + e_4 * fs_84_2717_30030 * r_2 * h8_m6 - e_4 * fs_11676_2717_2 * r_2 * h8_m2 + e_4 * fs_171_374_330 * r_4 * h6_m6 + e_4 * fs_437_374_6 * r_4 * h6_m2 + e_4 * fs_408_143_7 * r_6 * h4_m2 + e_4 * fs_125_143_105 * r_8 * h2_m2 - e_5 * fs_28_4199_6006 * r_2 * h10_m6 - e_5 * fs_28_4199_231 * r_2 * h10_m2 - e_5 * fs_4_2717_30030 * r_4 * h8_m6 + e_5 * fs_556_2717_2 * r_4 * h8_m2 - e_5 * fs_3_187_330 * r_6 * h6_m6 - e_5 * fs_23_561_6 * r_6 * h6_m2 - e_5 * fs_12_143_7 * r_8 * h4_m2 - e_5 * fs_10_429_105 * r_10 * h2_m2;
    }

    // NOTE: the angular components are formed in 93 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph2_m1, ph4_m1, ph6_m1, ph8_m8, ph8_m7, ph8_m1, ph10_m8, ph10_m7, ph10_m1, ab_2, pc_25, pc_26 : simd::cache_line_size())
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
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;

        const auto h2_m1 = ph2_m1[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h8_m8 = ph8_m8[k];
        const auto h8_m7 = ph8_m7[k];
        const auto h8_m1 = ph8_m1[k];
        const auto h10_m8 = ph10_m8[k];
        const auto h10_m7 = ph10_m7[k];
        const auto h10_m1 = ph10_m1[k];

        pc_25[k] = e_0 * fs_1575_8_30 * h2_m1 + e_1 * f_765_8 * h4_m1 - e_1 * fs_1125_4_30 * r_2 * h2_m1 - e_2 * fs_475_44_210 * h6_m1 - e_2 * f_765_11 * r_2 * h4_m1 + e_2 * fs_125_1_30 * r_4 * h2_m1 - e_3 * fs_1029_572_286 * h8_m7 - e_3 * fs_3675_572_10 * h8_m1 + e_3 * fs_95_22_210 * r_2 * h6_m1 + e_3 * f_2295_143 * r_4 * h4_m1 - e_3 * fs_250_11_30 * r_6 * h2_m1 + e_4 * fs_322_4199_7293 * h10_m7 - e_4 * fs_483_4199_22 * h10_m1 + e_4 * fs_1029_2717_286 * r_2 * h8_m7 + e_4 * fs_3675_2717_10 * r_2 * h8_m1 - e_4 * fs_95_187_210 * r_4 * h6_m1 - e_4 * f_204_143 * r_6 * h4_m1 + e_4 * fs_250_143_30 * r_8 * h2_m1 - e_5 * fs_28_4199_7293 * r_2 * h10_m7 + e_5 * fs_42_4199_22 * r_2 * h10_m1 - e_5 * fs_49_2717_286 * r_4 * h8_m7 - e_5 * fs_175_2717_10 * r_4 * h8_m1 + e_5 * fs_10_561_210 * r_6 * h6_m1 + e_5 * f_6_143 * r_8 * h4_m1 - e_5 * fs_20_429_30 * r_10 * h2_m1;

        pc_26[k] = e_3 * fs_294_143_143 * h8_m8 + e_4 * fs_483_4199_2431 * h10_m8 - e_4 * fs_1176_2717_143 * r_2 * h8_m8 - e_5 * fs_42_4199_2431 * r_2 * h10_m8 + e_5 * fs_56_2717_143 * r_4 * h8_m8;
    }

    // NOTE: the angular components are formed in 93 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph2_p1, ph4_p1, ph6_p1, ph8_p1, ph8_p7, ph10_p1, ph10_p7, ab_2, pc_27 : simd::cache_line_size())
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
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;

        const auto h2_p1 = ph2_p1[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h8_p7 = ph8_p7[k];
        const auto h10_p1 = ph10_p1[k];
        const auto h10_p7 = ph10_p7[k];

        pc_27[k] = - e_0 * fs_4725_32_2 * h2_p1 - e_1 * fs_255_2_15 * h4_p1 + e_1 * fs_3375_16_2 * r_2 * h2_p1 - e_2 * fs_1995_44_14 * h6_p1 + e_2 * fs_1020_11_15 * r_2 * h4_p1 - e_2 * fs_375_4_2 * r_4 * h2_p1 - e_3 * fs_882_143_6 * h8_p1 - e_3 * fs_147_286_4290 * h8_p7 + e_3 * fs_399_22_14 * r_2 * h6_p1 - e_3 * fs_3060_143_15 * r_4 * h4_p1 + e_3 * fs_375_22_2 * r_6 * h2_p1 - e_4 * fs_161_8398_330 * h10_p1 - e_4 * fs_161_4199_12155 * h10_p7 + e_4 * fs_3528_2717_6 * r_2 * h8_p1 + e_4 * fs_294_2717_4290 * r_2 * h8_p7 - e_4 * fs_399_187_14 * r_4 * h6_p1 + e_4 * fs_272_143_15 * r_6 * h4_p1 - e_4 * fs_375_286_2 * r_8 * h2_p1 + e_5 * fs_7_4199_330 * r_2 * h10_p1 + e_5 * fs_14_4199_12155 * r_2 * h10_p7 - e_5 * fs_168_2717_6 * r_4 * h8_p1 - e_5 * fs_14_2717_4290 * r_4 * h8_p7 + e_5 * fs_14_187_14 * r_6 * h6_p1 - e_5 * fs_8_143_15 * r_8 * h4_p1 + e_5 * fs_5_143_2 * r_10 * h2_p1;
    }

    // NOTE: the angular components are formed in 93 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph2_0, ph4_0, ph6_0, ph6_p6, ph8_0, ph8_p6, ph10_0, ph10_p6, ab_2, pc_28 : simd::cache_line_size())
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
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;

        const auto h2_0 = ph2_0[k];
        const auto h4_0 = ph4_0[k];
        const auto h6_0 = ph6_0[k];
        const auto h6_p6 = ph6_p6[k];
        const auto h8_0 = ph8_0[k];
        const auto h8_p6 = ph8_p6[k];
        const auto h10_0 = ph10_0[k];
        const auto h10_p6 = ph10_p6[k];

        pc_28[k] = e_0 * fs_4725_8_3 * h2_0 + e_1 * fs_1275_4_3 * h4_0 - e_1 * fs_3375_4_3 * r_2 * h2_0 - e_2 * fs_1995_22_3 * h6_0 + e_2 * fs_855_88_154 * h6_p6 - e_2 * fs_2550_11_3 * r_2 * h4_0 + e_2 * fs_375_1_3 * r_4 * h2_0 - e_3 * fs_3087_143_3 * h8_0 + e_3 * fs_441_572_286 * h8_p6 + e_3 * fs_399_11_3 * r_2 * h6_0 - e_3 * fs_171_44_154 * r_2 * h6_p6 + e_3 * fs_7650_143_3 * r_4 * h4_0 - e_3 * fs_750_11_3 * r_6 * h2_0 - e_4 * fs_3220_4199_3 * h10_0 - e_4 * fs_644_4199_1430 * h10_p6 + e_4 * fs_12348_2717_3 * r_2 * h8_0 - e_4 * fs_441_2717_286 * r_2 * h8_p6 - e_4 * fs_798_187_3 * r_4 * h6_0 + e_4 * fs_171_374_154 * r_4 * h6_p6 - e_4 * fs_680_143_3 * r_6 * h4_0 + e_4 * fs_750_143_3 * r_8 * h2_0 + e_5 * fs_280_4199_3 * r_2 * h10_0 + e_5 * fs_56_4199_1430 * r_2 * h10_p6 - e_5 * fs_588_2717_3 * r_4 * h8_0 + e_5 * fs_21_2717_286 * r_4 * h8_p6 + e_5 * fs_28_187_3 * r_6 * h6_0 - e_5 * fs_3_187_154 * r_6 * h6_p6 + e_5 * fs_20_143_3 * r_8 * h4_0 - e_5 * fs_20_143_3 * r_10 * h2_0;
    }

    // NOTE: the angular components are formed in 93 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph2_p1, ph4_p1, ph6_p1, ph6_p5, ph8_p1, ph8_p5, ph10_p1, ph10_p5, ab_2, pc_29 : simd::cache_line_size())
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
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;

        const auto h2_p1 = ph2_p1[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h6_p5 = ph6_p5[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h8_p5 = ph8_p5[k];
        const auto h10_p1 = ph10_p1[k];
        const auto h10_p5 = ph10_p5[k];

        pc_29[k] = e_0 * fs_4725_16_14 * h2_p1 - e_1 * fs_255_8_105 * h4_p1 - e_1 * fs_3375_8_14 * r_2 * h2_p1 - e_2 * fs_285_11_2 * h6_p1 - e_2 * fs_285_11_33 * h6_p5 + e_2 * fs_255_11_105 * r_2 * h4_p1 + e_2 * fs_375_2_14 * r_4 * h2_p1 + e_3 * fs_2667_572_42 * h8_p1 + e_3 * fs_567_572_858 * h8_p5 + e_3 * fs_114_11_2 * r_2 * h6_p1 + e_3 * fs_114_11_33 * r_2 * h6_p5 - e_3 * fs_765_143_105 * r_4 * h4_p1 - e_3 * fs_375_11_14 * r_6 * h2_p1 + e_4 * fs_161_4199_2310 * h10_p1 - e_4 * fs_805_4199_1001 * h10_p5 - e_4 * fs_2667_2717_42 * r_2 * h8_p1 - e_4 * fs_567_2717_858 * r_2 * h8_p5 - e_4 * fs_228_187_2 * r_4 * h6_p1 - e_4 * fs_228_187_33 * r_4 * h6_p5 + e_4 * fs_68_143_105 * r_6 * h4_p1 + e_4 * fs_375_143_14 * r_8 * h2_p1 - e_5 * fs_14_4199_2310 * r_2 * h10_p1 + e_5 * fs_70_4199_1001 * r_2 * h10_p5 + e_5 * fs_127_2717_42 * r_4 * h8_p1 + e_5 * fs_27_2717_858 * r_4 * h8_p5 + e_5 * fs_8_187_2 * r_6 * h6_p1 + e_5 * fs_8_187_33 * r_6 * h6_p5 - e_5 * fs_2_143_105 * r_8 * h4_p1 - e_5 * fs_10_143_14 * r_10 * h2_p1;
    }

    // NOTE: the angular components are formed in 93 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph2_p2, ph4_p2, ph4_p4, ph6_p2, ph6_p4, ph8_p2, ph8_p4, ph10_p2, ph10_p4, ab_2, pc_30 : simd::cache_line_size())
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
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;

        const auto h2_p2 = ph2_p2[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h4_p4 = ph4_p4[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h6_p4 = ph6_p4[k];
        const auto h8_p2 = ph8_p2[k];
        const auto h8_p4 = ph8_p4[k];
        const auto h10_p2 = ph10_p2[k];
        const auto h10_p4 = ph10_p4[k];

        pc_30[k] = e_0 * fs_4725_16_7 * h2_p2 - e_1 * fs_255_8_105 * h4_p2 + e_1 * fs_255_2_15 * h4_p4 - e_1 * fs_3375_8_7 * r_2 * h2_p2 + e_2 * fs_3705_88_10 * h6_p2 - e_2 * fs_1995_44_3 * h6_p4 + e_2 * fs_255_11_105 * r_2 * h4_p2 - e_2 * fs_1020_11_15 * r_2 * h4_p4 + e_2 * fs_375_2_7 * r_4 * h2_p2 - e_3 * fs_2331_572_30 * h8_p2 + e_3 * fs_903_286_33 * h8_p4 - e_3 * fs_741_44_10 * r_2 * h6_p2 + e_3 * fs_399_22_3 * r_2 * h6_p4 - e_3 * fs_765_143_105 * r_4 * h4_p2 + e_3 * fs_3060_143_15 * r_4 * h4_p4 - e_3 * fs_375_11_7 * r_6 * h2_p2 - e_4 * fs_644_4199_385 * h10_p2 - e_4 * fs_322_4199_5005 * h10_p4 + e_4 * fs_2331_2717_30 * r_2 * h8_p2 - e_4 * fs_1806_2717_33 * r_2 * h8_p4 + e_4 * fs_741_374_10 * r_4 * h6_p2 - e_4 * fs_399_187_3 * r_4 * h6_p4 + e_4 * fs_68_143_105 * r_6 * h4_p2 - e_4 * fs_272_143_15 * r_6 * h4_p4 + e_4 * fs_375_143_7 * r_8 * h2_p2 + e_5 * fs_56_4199_385 * r_2 * h10_p2 + e_5 * fs_28_4199_5005 * r_2 * h10_p4 - e_5 * fs_111_2717_30 * r_4 * h8_p2 + e_5 * fs_86_2717_33 * r_4 * h8_p4 - e_5 * fs_13_187_10 * r_6 * h6_p2 + e_5 * fs_14_187_3 * r_6 * h6_p4 - e_5 * fs_2_143_105 * r_8 * h4_p2 + e_5 * fs_8_143_15 * r_8 * h4_p4 - e_5 * fs_10_143_7 * r_10 * h2_p2;
    }

    // NOTE: the angular components are formed in 93 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph4_m3, ph6_m3, ph8_m3, ph10_m3, ab_2, pc_31 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

        const auto h4_m3 = ph4_m3[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h8_m3 = ph8_m3[k];
        const auto h10_m3 = ph10_m3[k];

        pc_31[k] = e_1 * fs_1275_4_3 * h4_m3 - e_2 * f_2565_22 * h6_m3 - e_2 * fs_2550_11_3 * r_2 * h4_m3 + e_3 * fs_315_286_22 * h8_m3 + e_3 * f_513_11 * r_2 * h6_m3 + e_3 * fs_7650_143_3 * r_4 * h4_m3 + e_4 * fs_805_4199_1001 * h10_m3 - e_4 * fs_630_2717_22 * r_2 * h8_m3 - e_4 * f_1026_187 * r_4 * h6_m3 - e_4 * fs_680_143_3 * r_6 * h4_m3 - e_5 * fs_70_4199_1001 * r_2 * h10_m3 + e_5 * fs_30_2717_22 * r_4 * h8_m3 + e_5 * f_36_187 * r_6 * h6_m3 + e_5 * fs_20_143_3 * r_8 * h4_m3;
    }

    // NOTE: the angular components are formed in 93 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph2_m2, ph4_m4, ph4_m2, ph6_m4, ph6_m2, ph8_m4, ph8_m2, ph10_m4, ph10_m2, ab_2, pc_32 : simd::cache_line_size())
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
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;

        const auto h2_m2 = ph2_m2[k];
        const auto h4_m4 = ph4_m4[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h6_m4 = ph6_m4[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h8_m4 = ph8_m4[k];
        const auto h8_m2 = ph8_m2[k];
        const auto h10_m4 = ph10_m4[k];
        const auto h10_m2 = ph10_m2[k];

        pc_32[k] = e_0 * fs_4725_16_7 * h2_m2 - e_1 * fs_255_2_15 * h4_m4 - e_1 * fs_255_8_105 * h4_m2 - e_1 * fs_3375_8_7 * r_2 * h2_m2 + e_2 * fs_1995_44_3 * h6_m4 + e_2 * fs_3705_88_10 * h6_m2 + e_2 * fs_1020_11_15 * r_2 * h4_m4 + e_2 * fs_255_11_105 * r_2 * h4_m2 + e_2 * fs_375_2_7 * r_4 * h2_m2 - e_3 * fs_903_286_33 * h8_m4 - e_3 * fs_2331_572_30 * h8_m2 - e_3 * fs_399_22_3 * r_2 * h6_m4 - e_3 * fs_741_44_10 * r_2 * h6_m2 - e_3 * fs_3060_143_15 * r_4 * h4_m4 - e_3 * fs_765_143_105 * r_4 * h4_m2 - e_3 * fs_375_11_7 * r_6 * h2_m2 + e_4 * fs_322_4199_5005 * h10_m4 - e_4 * fs_644_4199_385 * h10_m2 + e_4 * fs_1806_2717_33 * r_2 * h8_m4 + e_4 * fs_2331_2717_30 * r_2 * h8_m2 + e_4 * fs_399_187_3 * r_4 * h6_m4 + e_4 * fs_741_374_10 * r_4 * h6_m2 + e_4 * fs_272_143_15 * r_6 * h4_m4 + e_4 * fs_68_143_105 * r_6 * h4_m2 + e_4 * fs_375_143_7 * r_8 * h2_m2 - e_5 * fs_28_4199_5005 * r_2 * h10_m4 + e_5 * fs_56_4199_385 * r_2 * h10_m2 - e_5 * fs_86_2717_33 * r_4 * h8_m4 - e_5 * fs_111_2717_30 * r_4 * h8_m2 - e_5 * fs_14_187_3 * r_6 * h6_m4 - e_5 * fs_13_187_10 * r_6 * h6_m2 - e_5 * fs_8_143_15 * r_8 * h4_m4 - e_5 * fs_2_143_105 * r_8 * h4_m2 - e_5 * fs_10_143_7 * r_10 * h2_m2;
    }

    // NOTE: the angular components are formed in 93 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph2_m1, ph4_m1, ph6_m6, ph6_m5, ph6_m1, ph8_m6, ph8_m5, ph8_m1, ph10_m6, ph10_m5, ph10_m1, ab_2, pc_33, pc_34 : simd::cache_line_size())
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
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;

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

        pc_33[k] = e_0 * fs_4725_16_14 * h2_m1 - e_1 * fs_255_8_105 * h4_m1 - e_1 * fs_3375_8_14 * r_2 * h2_m1 + e_2 * fs_285_11_33 * h6_m5 - e_2 * fs_285_11_2 * h6_m1 + e_2 * fs_255_11_105 * r_2 * h4_m1 + e_2 * fs_375_2_14 * r_4 * h2_m1 - e_3 * fs_567_572_858 * h8_m5 + e_3 * fs_2667_572_42 * h8_m1 - e_3 * fs_114_11_33 * r_2 * h6_m5 + e_3 * fs_114_11_2 * r_2 * h6_m1 - e_3 * fs_765_143_105 * r_4 * h4_m1 - e_3 * fs_375_11_14 * r_6 * h2_m1 + e_4 * fs_805_4199_1001 * h10_m5 + e_4 * fs_161_4199_2310 * h10_m1 + e_4 * fs_567_2717_858 * r_2 * h8_m5 - e_4 * fs_2667_2717_42 * r_2 * h8_m1 + e_4 * fs_228_187_33 * r_4 * h6_m5 - e_4 * fs_228_187_2 * r_4 * h6_m1 + e_4 * fs_68_143_105 * r_6 * h4_m1 + e_4 * fs_375_143_14 * r_8 * h2_m1 - e_5 * fs_70_4199_1001 * r_2 * h10_m5 - e_5 * fs_14_4199_2310 * r_2 * h10_m1 - e_5 * fs_27_2717_858 * r_4 * h8_m5 + e_5 * fs_127_2717_42 * r_4 * h8_m1 - e_5 * fs_8_187_33 * r_6 * h6_m5 + e_5 * fs_8_187_2 * r_6 * h6_m1 - e_5 * fs_2_143_105 * r_8 * h4_m1 - e_5 * fs_10_143_14 * r_10 * h2_m1;

        pc_34[k] = - e_2 * fs_855_88_154 * h6_m6 - e_3 * fs_441_572_286 * h8_m6 + e_3 * fs_171_44_154 * r_2 * h6_m6 + e_4 * fs_644_4199_1430 * h10_m6 + e_4 * fs_441_2717_286 * r_2 * h8_m6 - e_4 * fs_171_374_154 * r_4 * h6_m6 - e_5 * fs_56_4199_1430 * r_2 * h10_m6 - e_5 * fs_21_2717_286 * r_4 * h8_m6 + e_5 * fs_3_187_154 * r_6 * h6_m6;
    }

    // NOTE: the angular components are formed in 93 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph2_m1, ph4_m1, ph6_m1, ph8_m7, ph8_m1, ph10_m7, ph10_m1, ab_2, pc_35 : simd::cache_line_size())
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
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;

        const auto h2_m1 = ph2_m1[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h8_m7 = ph8_m7[k];
        const auto h8_m1 = ph8_m1[k];
        const auto h10_m7 = ph10_m7[k];
        const auto h10_m1 = ph10_m1[k];

        pc_35[k] = e_0 * fs_4725_32_2 * h2_m1 + e_1 * fs_255_2_15 * h4_m1 - e_1 * fs_3375_16_2 * r_2 * h2_m1 + e_2 * fs_1995_44_14 * h6_m1 - e_2 * fs_1020_11_15 * r_2 * h4_m1 + e_2 * fs_375_4_2 * r_4 * h2_m1 + e_3 * fs_147_286_4290 * h8_m7 + e_3 * fs_882_143_6 * h8_m1 - e_3 * fs_399_22_14 * r_2 * h6_m1 + e_3 * fs_3060_143_15 * r_4 * h4_m1 - e_3 * fs_375_22_2 * r_6 * h2_m1 + e_4 * fs_161_4199_12155 * h10_m7 + e_4 * fs_161_8398_330 * h10_m1 - e_4 * fs_294_2717_4290 * r_2 * h8_m7 - e_4 * fs_3528_2717_6 * r_2 * h8_m1 + e_4 * fs_399_187_14 * r_4 * h6_m1 - e_4 * fs_272_143_15 * r_6 * h4_m1 + e_4 * fs_375_286_2 * r_8 * h2_m1 - e_5 * fs_14_4199_12155 * r_2 * h10_m7 - e_5 * fs_7_4199_330 * r_2 * h10_m1 + e_5 * fs_14_2717_4290 * r_4 * h8_m7 + e_5 * fs_168_2717_6 * r_4 * h8_m1 - e_5 * fs_14_187_14 * r_6 * h6_m1 + e_5 * fs_8_143_15 * r_8 * h4_m1 - e_5 * fs_5_143_2 * r_10 * h2_m1;
    }

    // NOTE: the angular components are formed in 93 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph2_p2, ph4_p2, ph6_p2, ph6_p6, ph8_p2, ph8_p6, ph10_p2, ph10_p6, ab_2, pc_36 : simd::cache_line_size())
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
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;

        const auto h2_p2 = ph2_p2[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h6_p6 = ph6_p6[k];
        const auto h8_p2 = ph8_p2[k];
        const auto h8_p6 = ph8_p6[k];
        const auto h10_p2 = ph10_p2[k];
        const auto h10_p6 = ph10_p6[k];

        pc_36[k] = e_0 * fs_1575_32_2 * h2_p2 + e_1 * fs_255_4_30 * h4_p2 - e_1 * fs_1125_16_2 * r_2 * h2_p2 + e_2 * fs_665_22_35 * h6_p2 - e_2 * fs_285_44_77 * h6_p6 - e_2 * fs_510_11_30 * r_2 * h4_p2 + e_2 * fs_125_4_2 * r_4 * h2_p2 + e_3 * fs_294_143_105 * h8_p2 - e_3 * fs_441_143_143 * h8_p6 - e_3 * fs_133_11_35 * r_2 * h6_p2 + e_3 * fs_57_22_77 * r_2 * h6_p6 + e_3 * fs_1530_143_30 * r_4 * h4_p2 - e_3 * fs_125_22_2 * r_6 * h2_p2 + e_4 * fs_483_8398_110 * h10_p2 - e_4 * fs_483_4199_715 * h10_p6 - e_4 * fs_1176_2717_105 * r_2 * h8_p2 + e_4 * fs_1764_2717_143 * r_2 * h8_p6 + e_4 * fs_266_187_35 * r_4 * h6_p2 - e_4 * fs_57_187_77 * r_4 * h6_p6 - e_4 * fs_136_143_30 * r_6 * h4_p2 + e_4 * fs_125_286_2 * r_8 * h2_p2 - e_5 * fs_21_4199_110 * r_2 * h10_p2 + e_5 * fs_42_4199_715 * r_2 * h10_p6 + e_5 * fs_56_2717_105 * r_4 * h8_p2 - e_5 * fs_84_2717_143 * r_4 * h8_p6 - e_5 * fs_28_561_35 * r_6 * h6_p2 + e_5 * fs_2_187_77 * r_6 * h6_p6 + e_5 * fs_4_143_30 * r_8 * h4_p2 - e_5 * fs_5_429_2 * r_10 * h2_p2;
    }

    // NOTE: the angular components are formed in 93 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph2_p1, ph4_p1, ph6_p1, ph6_p5, ph8_p1, ph8_p5, ph10_p1, ph10_p5, ab_2, pc_37 : simd::cache_line_size())
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
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;

        const auto h2_p1 = ph2_p1[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h6_p5 = ph6_p5[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h8_p5 = ph8_p5[k];
        const auto h10_p1 = ph10_p1[k];
        const auto h10_p5 = ph10_p5[k];

        pc_37[k] = - e_0 * f_1575_4 * h2_p1 - e_1 * fs_765_8_30 * h4_p1 + e_1 * f_1125_2 * r_2 * h2_p1 + e_2 * fs_665_44_7 * h6_p1 + e_2 * fs_665_88_462 * h6_p5 + e_2 * fs_765_11_30 * r_2 * h4_p1 - e_2 * f_250_1 * r_4 * h2_p1 + e_3 * fs_2499_143_3 * h8_p1 - e_3 * fs_21_286_3003 * h8_p5 - e_3 * fs_133_22_7 * r_2 * h6_p1 - e_3 * fs_133_44_462 * r_2 * h6_p5 - e_3 * fs_2295_143_30 * r_4 * h4_p1 + e_3 * f_500_11 * r_6 * h2_p1 + e_4 * fs_483_4199_165 * h10_p1 - e_4 * fs_2415_8398_286 * h10_p5 - e_4 * fs_9996_2717_3 * r_2 * h8_p1 + e_4 * fs_42_2717_3003 * r_2 * h8_p5 + e_4 * fs_133_187_7 * r_4 * h6_p1 + e_4 * fs_133_374_462 * r_4 * h6_p5 + e_4 * fs_204_143_30 * r_6 * h4_p1 - e_4 * f_500_143 * r_8 * h2_p1 - e_5 * fs_42_4199_165 * r_2 * h10_p1 + e_5 * fs_105_4199_286 * r_2 * h10_p5 + e_5 * fs_476_2717_3 * r_4 * h8_p1 - e_5 * fs_2_2717_3003 * r_4 * h8_p5 - e_5 * fs_14_561_7 * r_6 * h6_p1 - e_5 * fs_7_561_462 * r_6 * h6_p5 - e_5 * fs_6_143_30 * r_8 * h4_p1 + e_5 * f_40_429 * r_10 * h2_p1;
    }

    // NOTE: the angular components are formed in 93 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph2_0, ph4_p4, ph6_0, ph6_p4, ph8_0, ph8_p4, ph10_0, ph10_p4, ab_2, pc_38 : simd::cache_line_size())
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
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;

        const auto h2_0 = ph2_0[k];
        const auto h4_p4 = ph4_p4[k];
        const auto h6_0 = ph6_0[k];
        const auto h6_p4 = ph6_p4[k];
        const auto h8_0 = ph8_0[k];
        const auto h8_p4 = ph8_p4[k];
        const auto h10_0 = ph10_0[k];
        const auto h10_p4 = ph10_p4[k];

        pc_38[k] = e_0 * fs_1575_8_42 * h2_0 - e_1 * fs_255_4_30 * h4_p4 - e_1 * fs_1125_4_42 * r_2 * h2_0 - e_2 * fs_95_4_42 * h6_0 - e_2 * fs_2185_88_6 * h6_p4 + e_2 * fs_510_11_30 * r_2 * h4_p4 + e_2 * fs_125_1_42 * r_4 * h2_0 + e_3 * fs_63_11_42 * h8_0 + e_3 * fs_777_286_66 * h8_p4 + e_3 * fs_19_2_42 * r_2 * h6_0 + e_3 * fs_437_44_6 * r_2 * h6_p4 - e_3 * fs_1530_143_30 * r_4 * h4_p4 - e_3 * fs_250_11_42 * r_6 * h2_0 + e_4 * fs_2415_4199_42 * h10_0 - e_4 * fs_483_8398_10010 * h10_p4 - e_4 * fs_252_209_42 * r_2 * h8_0 - e_4 * fs_1554_2717_66 * r_2 * h8_p4 - e_4 * fs_19_17_42 * r_4 * h6_0 - e_4 * fs_437_374_6 * r_4 * h6_p4 + e_4 * fs_136_143_30 * r_6 * h4_p4 + e_4 * fs_250_143_42 * r_8 * h2_0 - e_5 * fs_210_4199_42 * r_2 * h10_0 + e_5 * fs_21_4199_10010 * r_2 * h10_p4 + e_5 * fs_12_209_42 * r_4 * h8_0 + e_5 * fs_74_2717_66 * r_4 * h8_p4 + e_5 * fs_2_51_42 * r_6 * h6_0 + e_5 * fs_23_561_6 * r_6 * h6_p4 - e_5 * fs_4_143_30 * r_8 * h4_p4 - e_5 * fs_20_429_42 * r_10 * h2_0;
    }

    // NOTE: the angular components are formed in 93 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph2_p1, ph4_p1, ph4_p3, ph6_p1, ph6_p3, ph8_p1, ph8_p3, ph10_p1, ph10_p3, ab_2, pc_39 : simd::cache_line_size())
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
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;

        const auto h2_p1 = ph2_p1[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h4_p3 = ph4_p3[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h8_p3 = ph8_p3[k];
        const auto h10_p1 = ph10_p1[k];
        const auto h10_p3 = ph10_p3[k];

        pc_39[k] = e_0 * fs_1575_4_7 * h2_p1 - e_1 * fs_255_8_210 * h4_p1 + e_1 * fs_765_8_30 * h4_p3 - e_1 * fs_1125_2_7 * r_2 * h2_p1 + e_2 * f_5035_44 * h6_p1 - e_2 * fs_3705_88_10 * h6_p3 + e_2 * fs_255_11_210 * r_2 * h4_p1 - e_2 * fs_765_11_30 * r_2 * h4_p3 + e_2 * fs_250_1_7 * r_4 * h2_p1 - e_3 * fs_609_286_21 * h8_p1 + e_3 * fs_504_143_55 * h8_p3 - e_3 * f_1007_22 * r_2 * h6_p1 + e_3 * fs_741_44_10 * r_2 * h6_p3 - e_3 * fs_765_143_210 * r_4 * h4_p1 + e_3 * fs_2295_143_30 * r_4 * h4_p3 - e_3 * fs_500_11_7 * r_6 * h2_p1 - e_4 * fs_483_4199_1155 * h10_p1 - e_4 * fs_483_8398_10010 * h10_p3 + e_4 * fs_1218_2717_21 * r_2 * h8_p1 - e_4 * fs_2016_2717_55 * r_2 * h8_p3 + e_4 * f_1007_187 * r_4 * h6_p1 - e_4 * fs_741_374_10 * r_4 * h6_p3 + e_4 * fs_68_143_210 * r_6 * h4_p1 - e_4 * fs_204_143_30 * r_6 * h4_p3 + e_4 * fs_500_143_7 * r_8 * h2_p1 + e_5 * fs_42_4199_1155 * r_2 * h10_p1 + e_5 * fs_21_4199_10010 * r_2 * h10_p3 - e_5 * fs_58_2717_21 * r_4 * h8_p1 + e_5 * fs_96_2717_55 * r_4 * h8_p3 - e_5 * f_106_561 * r_6 * h6_p1 + e_5 * fs_13_187_10 * r_6 * h6_p3 - e_5 * fs_2_143_210 * r_8 * h4_p1 + e_5 * fs_6_143_30 * r_8 * h4_p3 - e_5 * fs_40_429_7 * r_10 * h2_p1;
    }

    // NOTE: the angular components are formed in 93 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph2_m2, ph6_m2, ph8_m2, ph10_m2, ab_2, pc_40 : simd::cache_line_size())
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
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;

        const auto h2_m2 = ph2_m2[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h8_m2 = ph8_m2[k];
        const auto h10_m2 = ph10_m2[k];

        pc_40[k] = e_0 * fs_1575_16_70 * h2_m2 - e_1 * fs_1125_8_70 * r_2 * h2_m2 + e_2 * f_95_4 * h6_m2 + e_2 * fs_125_2_70 * r_4 * h2_m2 - e_3 * fs_105_11_3 * h8_m2 - e_3 * f_19_2 * r_2 * h6_m2 - e_3 * fs_125_11_70 * r_6 * h2_m2 + e_4 * fs_2415_4199_154 * h10_m2 + e_4 * fs_420_209_3 * r_2 * h8_m2 + e_4 * f_19_17 * r_4 * h6_m2 + e_4 * fs_125_143_70 * r_8 * h2_m2 - e_5 * fs_210_4199_154 * r_2 * h10_m2 - e_5 * fs_20_209_3 * r_4 * h8_m2 - e_5 * f_2_51 * r_6 * h6_m2 - e_5 * fs_10_429_70 * r_10 * h2_m2;
    }

    // NOTE: the angular components are formed in 93 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph2_m1, ph4_m3, ph4_m1, ph6_m3, ph6_m1, ph8_m3, ph8_m1, ph10_m3, ph10_m1, ab_2, pc_41 : simd::cache_line_size())
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
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;

        const auto h2_m1 = ph2_m1[k];
        const auto h4_m3 = ph4_m3[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h8_m3 = ph8_m3[k];
        const auto h8_m1 = ph8_m1[k];
        const auto h10_m3 = ph10_m3[k];
        const auto h10_m1 = ph10_m1[k];

        pc_41[k] = e_0 * fs_1575_4_7 * h2_m1 - e_1 * fs_765_8_30 * h4_m3 - e_1 * fs_255_8_210 * h4_m1 - e_1 * fs_1125_2_7 * r_2 * h2_m1 + e_2 * fs_3705_88_10 * h6_m3 + e_2 * f_5035_44 * h6_m1 + e_2 * fs_765_11_30 * r_2 * h4_m3 + e_2 * fs_255_11_210 * r_2 * h4_m1 + e_2 * fs_250_1_7 * r_4 * h2_m1 - e_3 * fs_504_143_55 * h8_m3 - e_3 * fs_609_286_21 * h8_m1 - e_3 * fs_741_44_10 * r_2 * h6_m3 - e_3 * f_1007_22 * r_2 * h6_m1 - e_3 * fs_2295_143_30 * r_4 * h4_m3 - e_3 * fs_765_143_210 * r_4 * h4_m1 - e_3 * fs_500_11_7 * r_6 * h2_m1 + e_4 * fs_483_8398_10010 * h10_m3 - e_4 * fs_483_4199_1155 * h10_m1 + e_4 * fs_2016_2717_55 * r_2 * h8_m3 + e_4 * fs_1218_2717_21 * r_2 * h8_m1 + e_4 * fs_741_374_10 * r_4 * h6_m3 + e_4 * f_1007_187 * r_4 * h6_m1 + e_4 * fs_204_143_30 * r_6 * h4_m3 + e_4 * fs_68_143_210 * r_6 * h4_m1 + e_4 * fs_500_143_7 * r_8 * h2_m1 - e_5 * fs_21_4199_10010 * r_2 * h10_m3 + e_5 * fs_42_4199_1155 * r_2 * h10_m1 - e_5 * fs_96_2717_55 * r_4 * h8_m3 - e_5 * fs_58_2717_21 * r_4 * h8_m1 - e_5 * fs_13_187_10 * r_6 * h6_m3 - e_5 * f_106_561 * r_6 * h6_m1 - e_5 * fs_6_143_30 * r_8 * h4_m3 - e_5 * fs_2_143_210 * r_8 * h4_m1 - e_5 * fs_40_429_7 * r_10 * h2_m1;
    }

    // NOTE: the angular components are formed in 93 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph4_m4, ph6_m4, ph8_m4, ph10_m4, ab_2, pc_42 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

        const auto h4_m4 = ph4_m4[k];
        const auto h6_m4 = ph6_m4[k];
        const auto h8_m4 = ph8_m4[k];
        const auto h10_m4 = ph10_m4[k];

        pc_42[k] = e_1 * fs_255_4_30 * h4_m4 + e_2 * fs_2185_88_6 * h6_m4 - e_2 * fs_510_11_30 * r_2 * h4_m4 - e_3 * fs_777_286_66 * h8_m4 - e_3 * fs_437_44_6 * r_2 * h6_m4 + e_3 * fs_1530_143_30 * r_4 * h4_m4 + e_4 * fs_483_8398_10010 * h10_m4 + e_4 * fs_1554_2717_66 * r_2 * h8_m4 + e_4 * fs_437_374_6 * r_4 * h6_m4 - e_4 * fs_136_143_30 * r_6 * h4_m4 - e_5 * fs_21_4199_10010 * r_2 * h10_m4 - e_5 * fs_74_2717_66 * r_4 * h8_m4 - e_5 * fs_23_561_6 * r_6 * h6_m4 + e_5 * fs_4_143_30 * r_8 * h4_m4;
    }

    // NOTE: the angular components are formed in 93 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph2_m1, ph4_m1, ph6_m5, ph6_m1, ph8_m5, ph8_m1, ph10_m5, ph10_m1, ab_2, pc_43 : simd::cache_line_size())
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
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;

        const auto h2_m1 = ph2_m1[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h6_m5 = ph6_m5[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h8_m5 = ph8_m5[k];
        const auto h8_m1 = ph8_m1[k];
        const auto h10_m5 = ph10_m5[k];
        const auto h10_m1 = ph10_m1[k];

        pc_43[k] = e_0 * f_1575_4 * h2_m1 + e_1 * fs_765_8_30 * h4_m1 - e_1 * f_1125_2 * r_2 * h2_m1 - e_2 * fs_665_88_462 * h6_m5 - e_2 * fs_665_44_7 * h6_m1 - e_2 * fs_765_11_30 * r_2 * h4_m1 + e_2 * f_250_1 * r_4 * h2_m1 + e_3 * fs_21_286_3003 * h8_m5 - e_3 * fs_2499_143_3 * h8_m1 + e_3 * fs_133_44_462 * r_2 * h6_m5 + e_3 * fs_133_22_7 * r_2 * h6_m1 + e_3 * fs_2295_143_30 * r_4 * h4_m1 - e_3 * f_500_11 * r_6 * h2_m1 + e_4 * fs_2415_8398_286 * h10_m5 - e_4 * fs_483_4199_165 * h10_m1 - e_4 * fs_42_2717_3003 * r_2 * h8_m5 + e_4 * fs_9996_2717_3 * r_2 * h8_m1 - e_4 * fs_133_374_462 * r_4 * h6_m5 - e_4 * fs_133_187_7 * r_4 * h6_m1 - e_4 * fs_204_143_30 * r_6 * h4_m1 + e_4 * f_500_143 * r_8 * h2_m1 - e_5 * fs_105_4199_286 * r_2 * h10_m5 + e_5 * fs_42_4199_165 * r_2 * h10_m1 + e_5 * fs_2_2717_3003 * r_4 * h8_m5 - e_5 * fs_476_2717_3 * r_4 * h8_m1 + e_5 * fs_7_561_462 * r_6 * h6_m5 + e_5 * fs_14_561_7 * r_6 * h6_m1 + e_5 * fs_6_143_30 * r_8 * h4_m1 - e_5 * f_40_429 * r_10 * h2_m1;
    }

    // NOTE: the angular components are formed in 93 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph2_m2, ph4_m2, ph6_m6, ph6_m2, ph8_m6, ph8_m2, ph10_m6, ph10_m2, ab_2, pc_44 : simd::cache_line_size())
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
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;

        const auto h2_m2 = ph2_m2[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h6_m6 = ph6_m6[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h8_m6 = ph8_m6[k];
        const auto h8_m2 = ph8_m2[k];
        const auto h10_m6 = ph10_m6[k];
        const auto h10_m2 = ph10_m2[k];

        pc_44[k] = - e_0 * fs_1575_32_2 * h2_m2 - e_1 * fs_255_4_30 * h4_m2 + e_1 * fs_1125_16_2 * r_2 * h2_m2 + e_2 * fs_285_44_77 * h6_m6 - e_2 * fs_665_22_35 * h6_m2 + e_2 * fs_510_11_30 * r_2 * h4_m2 - e_2 * fs_125_4_2 * r_4 * h2_m2 + e_3 * fs_441_143_143 * h8_m6 - e_3 * fs_294_143_105 * h8_m2 - e_3 * fs_57_22_77 * r_2 * h6_m6 + e_3 * fs_133_11_35 * r_2 * h6_m2 - e_3 * fs_1530_143_30 * r_4 * h4_m2 + e_3 * fs_125_22_2 * r_6 * h2_m2 + e_4 * fs_483_4199_715 * h10_m6 - e_4 * fs_483_8398_110 * h10_m2 - e_4 * fs_1764_2717_143 * r_2 * h8_m6 + e_4 * fs_1176_2717_105 * r_2 * h8_m2 + e_4 * fs_57_187_77 * r_4 * h6_m6 - e_4 * fs_266_187_35 * r_4 * h6_m2 + e_4 * fs_136_143_30 * r_6 * h4_m2 - e_4 * fs_125_286_2 * r_8 * h2_m2 - e_5 * fs_42_4199_715 * r_2 * h10_m6 + e_5 * fs_21_4199_110 * r_2 * h10_m2 + e_5 * fs_84_2717_143 * r_4 * h8_m6 - e_5 * fs_56_2717_105 * r_4 * h8_m2 - e_5 * fs_2_187_77 * r_6 * h6_m6 + e_5 * fs_28_561_35 * r_6 * h6_m2 - e_5 * fs_4_143_30 * r_8 * h4_m2 + e_5 * fs_5_429_2 * r_10 * h2_m2;
    }

    // NOTE: the angular components are formed in 93 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph4_p3, ph6_p3, ph6_p5, ph8_p3, ph8_p5, ph10_p3, ph10_p5, ab_2, pc_45 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

        const auto h4_p3 = ph4_p3[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h6_p5 = ph6_p5[k];
        const auto h8_p3 = ph8_p3[k];
        const auto h8_p5 = ph8_p5[k];
        const auto h10_p3 = ph10_p3[k];
        const auto h10_p5 = ph10_p5[k];

        pc_45[k] = - e_1 * fs_255_8_42 * h4_p3 - e_2 * fs_1995_44_14 * h6_p3 - e_2 * fs_95_44_2310 * h6_p5 + e_2 * fs_255_11_42 * r_2 * h4_p3 - e_3 * fs_441_143_77 * h8_p3 - e_3 * fs_42_143_15015 * h8_p5 + e_3 * fs_399_22_14 * r_2 * h6_p3 + e_3 * fs_19_22_2310 * r_2 * h6_p5 - e_3 * fs_765_143_42 * r_4 * h4_p3 - e_4 * fs_483_8398_286 * h10_p3 - e_4 * fs_483_8398_1430 * h10_p5 + e_4 * fs_1764_2717_77 * r_2 * h8_p3 + e_4 * fs_168_2717_15015 * r_2 * h8_p5 - e_4 * fs_399_187_14 * r_4 * h6_p3 - e_4 * fs_19_187_2310 * r_4 * h6_p5 + e_4 * fs_68_143_42 * r_6 * h4_p3 + e_5 * fs_21_4199_286 * r_2 * h10_p3 + e_5 * fs_21_4199_1430 * r_2 * h10_p5 - e_5 * fs_84_2717_77 * r_4 * h8_p3 - e_5 * fs_8_2717_15015 * r_4 * h8_p5 + e_5 * fs_14_187_14 * r_6 * h6_p3 + e_5 * fs_2_561_2310 * r_6 * h6_p5 - e_5 * fs_2_143_42 * r_8 * h4_p3;
    }

    // NOTE: the angular components are formed in 93 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph2_p2, ph4_p2, ph4_p4, ph6_p2, ph6_p4, ph8_p2, ph8_p4, ph10_p2, ph10_p4, ab_2, pc_46 : simd::cache_line_size())
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
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;

        const auto h2_p2 = ph2_p2[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h4_p4 = ph4_p4[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h6_p4 = ph6_p4[k];
        const auto h8_p2 = ph8_p2[k];
        const auto h8_p4 = ph8_p4[k];
        const auto h10_p2 = ph10_p2[k];
        const auto h10_p4 = ph10_p4[k];

        pc_46[k] = e_0 * fs_1575_32_10 * h2_p2 + e_1 * fs_3315_16_6 * h4_p2 + e_1 * fs_255_8_42 * h4_p4 - e_1 * fs_1125_16_10 * r_2 * h2_p2 + e_2 * fs_665_44_7 * h6_p2 + e_2 * fs_475_44_210 * h6_p4 - e_2 * fs_3315_22_6 * r_2 * h4_p2 - e_2 * fs_255_11_42 * r_2 * h4_p4 + e_2 * fs_125_4_10 * r_4 * h2_p2 - e_3 * fs_147_22_21 * h8_p2 - e_3 * fs_105_286_2310 * h8_p4 - e_3 * fs_133_22_7 * r_2 * h6_p2 - e_3 * fs_95_22_210 * r_2 * h6_p4 + e_3 * fs_765_22_6 * r_4 * h4_p2 + e_3 * fs_765_143_42 * r_4 * h4_p4 - e_3 * fs_125_22_10 * r_6 * h2_p2 - e_4 * fs_1932_4199_22 * h10_p2 - e_4 * fs_966_4199_286 * h10_p4 + e_4 * fs_294_209_21 * r_2 * h8_p2 + e_4 * fs_210_2717_2310 * r_2 * h8_p4 + e_4 * fs_133_187_7 * r_4 * h6_p2 + e_4 * fs_95_187_210 * r_4 * h6_p4 - e_4 * fs_34_11_6 * r_6 * h4_p2 - e_4 * fs_68_143_42 * r_6 * h4_p4 + e_4 * fs_125_286_10 * r_8 * h2_p2 + e_5 * fs_168_4199_22 * r_2 * h10_p2 + e_5 * fs_84_4199_286 * r_2 * h10_p4 - e_5 * fs_14_209_21 * r_4 * h8_p2 - e_5 * fs_10_2717_2310 * r_4 * h8_p4 - e_5 * fs_14_561_7 * r_6 * h6_p2 - e_5 * fs_10_561_210 * r_6 * h6_p4 + e_5 * fs_1_11_6 * r_8 * h4_p2 + e_5 * fs_2_143_42 * r_8 * h4_p4 - e_5 * fs_5_429_10 * r_10 * h2_p2;
    }

    // NOTE: the angular components are formed in 93 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph2_p1, ph4_p1, ph4_p3, ph6_p1, ph6_p3, ph8_p1, ph8_p3, ph10_p1, ph10_p3, ab_2, pc_47 : simd::cache_line_size())
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
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;

        const auto h2_p1 = ph2_p1[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h4_p3 = ph4_p3[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h8_p3 = ph8_p3[k];
        const auto h10_p1 = ph10_p1[k];
        const auto h10_p3 = ph10_p3[k];

        pc_47[k] = - e_0 * fs_1575_16_35 * h2_p1 - e_1 * fs_765_16_42 * h4_p1 - e_1 * fs_3315_16_6 * h4_p3 + e_1 * fs_1125_8_35 * r_2 * h2_p1 + e_2 * fs_665_11_5 * h6_p1 + e_2 * fs_285_11_2 * h6_p3 + e_2 * fs_765_22_42 * r_2 * h4_p1 + e_2 * fs_3315_22_6 * r_2 * h4_p3 - e_2 * fs_125_2_35 * r_4 * h2_p1 - e_3 * fs_483_286_105 * h8_p1 + e_3 * fs_63_22_11 * h8_p3 - e_3 * fs_266_11_5 * r_2 * h6_p1 - e_3 * fs_114_11_2 * r_2 * h6_p3 - e_3 * fs_2295_286_42 * r_4 * h4_p1 - e_3 * fs_765_22_6 * r_4 * h4_p3 + e_3 * fs_125_11_35 * r_6 * h2_p1 - e_4 * fs_966_4199_231 * h10_p1 - e_4 * fs_483_4199_2002 * h10_p3 + e_4 * fs_966_2717_105 * r_2 * h8_p1 - e_4 * fs_126_209_11 * r_2 * h8_p3 + e_4 * fs_532_187_5 * r_4 * h6_p1 + e_4 * fs_228_187_2 * r_4 * h6_p3 + e_4 * fs_102_143_42 * r_6 * h4_p1 + e_4 * fs_34_11_6 * r_6 * h4_p3 - e_4 * fs_125_143_35 * r_8 * h2_p1 + e_5 * fs_84_4199_231 * r_2 * h10_p1 + e_5 * fs_42_4199_2002 * r_2 * h10_p3 - e_5 * fs_46_2717_105 * r_4 * h8_p1 + e_5 * fs_6_209_11 * r_4 * h8_p3 - e_5 * fs_56_561_5 * r_6 * h6_p1 - e_5 * fs_8_187_2 * r_6 * h6_p3 - e_5 * fs_3_143_42 * r_8 * h4_p1 - e_5 * fs_1_11_6 * r_8 * h4_p3 + e_5 * fs_10_429_35 * r_10 * h2_p1;
    }

    // NOTE: the angular components are formed in 93 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph2_0, ph2_p2, ph4_0, ph4_p2, ph6_0, ph6_p2, ph8_0, ph8_p2, ph10_0, ph10_p2, ab_2, pc_48 : simd::cache_line_size())
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
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;

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

        pc_48[k] = e_0 * fs_1575_16_210 * h2_0 + e_0 * fs_1575_32_70 * h2_p2 - e_1 * fs_255_8_210 * h4_0 + e_1 * fs_765_16_42 * h4_p2 - e_1 * fs_1125_8_210 * r_2 * h2_0 - e_1 * fs_1125_16_70 * r_2 * h2_p2 + e_2 * fs_95_22_210 * h6_0 - e_2 * f_5035_44 * h6_p2 + e_2 * fs_255_11_210 * r_2 * h4_0 - e_2 * fs_765_22_42 * r_2 * h4_p2 + e_2 * fs_125_2_210 * r_4 * h2_0 + e_2 * fs_125_4_70 * r_4 * h2_p2 + e_3 * fs_63_143_210 * h8_0 + e_3 * fs_4179_286_3 * h8_p2 - e_3 * fs_19_11_210 * r_2 * h6_0 + e_3 * f_1007_22 * r_2 * h6_p2 - e_3 * fs_765_143_210 * r_4 * h4_0 + e_3 * fs_2295_286_42 * r_4 * h4_p2 - e_3 * fs_125_11_210 * r_6 * h2_0 - e_3 * fs_125_22_70 * r_6 * h2_p2 - e_4 * fs_1932_4199_210 * h10_0 - e_4 * fs_1932_4199_154 * h10_p2 - e_4 * fs_252_2717_210 * r_2 * h8_0 - e_4 * fs_8358_2717_3 * r_2 * h8_p2 + e_4 * fs_38_187_210 * r_4 * h6_0 - e_4 * f_1007_187 * r_4 * h6_p2 + e_4 * fs_68_143_210 * r_6 * h4_0 - e_4 * fs_102_143_42 * r_6 * h4_p2 + e_4 * fs_125_143_210 * r_8 * h2_0 + e_4 * fs_125_286_70 * r_8 * h2_p2 + e_5 * fs_168_4199_210 * r_2 * h10_0 + e_5 * fs_168_4199_154 * r_2 * h10_p2 + e_5 * fs_12_2717_210 * r_4 * h8_0 + e_5 * fs_398_2717_3 * r_4 * h8_p2 - e_5 * fs_4_561_210 * r_6 * h6_0 + e_5 * f_106_561 * r_6 * h6_p2 - e_5 * fs_2_143_210 * r_8 * h4_0 + e_5 * fs_3_143_42 * r_8 * h4_p2 - e_5 * fs_10_429_210 * r_10 * h2_0 - e_5 * fs_5_429_70 * r_10 * h2_p2;
    }

    // NOTE: the angular components are formed in 93 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph2_m2, ph2_m1, ph4_m2, ph4_m1, ph6_m2, ph6_m1, ph8_m2, ph8_m1, ph10_m2, ph10_m1, ab_2, pc_49, pc_50 : simd::cache_line_size())
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
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;

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

        pc_49[k] = e_0 * fs_7875_16_7 * h2_m1 - e_1 * fs_255_8_210 * h4_m1 - e_1 * fs_5625_8_7 * r_2 * h2_m1 + e_2 * f_1520_11 * h6_m1 + e_2 * fs_255_11_210 * r_2 * h4_m1 + e_2 * fs_625_2_7 * r_4 * h2_m1 - e_3 * fs_987_143_21 * h8_m1 - e_3 * f_608_11 * r_2 * h6_m1 - e_3 * fs_765_143_210 * r_4 * h4_m1 - e_3 * fs_625_11_7 * r_6 * h2_m1 + e_4 * fs_966_4199_1155 * h10_m1 + e_4 * fs_3948_2717_21 * r_2 * h8_m1 + e_4 * f_1216_187 * r_4 * h6_m1 + e_4 * fs_68_143_210 * r_6 * h4_m1 + e_4 * fs_625_143_7 * r_8 * h2_m1 - e_5 * fs_84_4199_1155 * r_2 * h10_m1 - e_5 * fs_188_2717_21 * r_4 * h8_m1 - e_5 * f_128_561 * r_6 * h6_m1 - e_5 * fs_2_143_210 * r_8 * h4_m1 - e_5 * fs_50_429_7 * r_10 * h2_m1;

        pc_50[k] = - e_0 * fs_1575_32_70 * h2_m2 - e_1 * fs_765_16_42 * h4_m2 + e_1 * fs_1125_16_70 * r_2 * h2_m2 + e_2 * f_5035_44 * h6_m2 + e_2 * fs_765_22_42 * r_2 * h4_m2 - e_2 * fs_125_4_70 * r_4 * h2_m2 - e_3 * fs_4179_286_3 * h8_m2 - e_3 * f_1007_22 * r_2 * h6_m2 - e_3 * fs_2295_286_42 * r_4 * h4_m2 + e_3 * fs_125_22_70 * r_6 * h2_m2 + e_4 * fs_1932_4199_154 * h10_m2 + e_4 * fs_8358_2717_3 * r_2 * h8_m2 + e_4 * f_1007_187 * r_4 * h6_m2 + e_4 * fs_102_143_42 * r_6 * h4_m2 - e_4 * fs_125_286_70 * r_8 * h2_m2 - e_5 * fs_168_4199_154 * r_2 * h10_m2 - e_5 * fs_398_2717_3 * r_4 * h8_m2 - e_5 * f_106_561 * r_6 * h6_m2 - e_5 * fs_3_143_42 * r_8 * h4_m2 + e_5 * fs_5_429_70 * r_10 * h2_m2;
    }

    // NOTE: the angular components are formed in 93 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph2_m1, ph4_m3, ph4_m1, ph6_m3, ph6_m1, ph8_m3, ph8_m1, ph10_m3, ph10_m1, ab_2, pc_51 : simd::cache_line_size())
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
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;

        const auto h2_m1 = ph2_m1[k];
        const auto h4_m3 = ph4_m3[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h8_m3 = ph8_m3[k];
        const auto h8_m1 = ph8_m1[k];
        const auto h10_m3 = ph10_m3[k];
        const auto h10_m1 = ph10_m1[k];

        pc_51[k] = e_0 * fs_1575_16_35 * h2_m1 + e_1 * fs_3315_16_6 * h4_m3 + e_1 * fs_765_16_42 * h4_m1 - e_1 * fs_1125_8_35 * r_2 * h2_m1 - e_2 * fs_285_11_2 * h6_m3 - e_2 * fs_665_11_5 * h6_m1 - e_2 * fs_3315_22_6 * r_2 * h4_m3 - e_2 * fs_765_22_42 * r_2 * h4_m1 + e_2 * fs_125_2_35 * r_4 * h2_m1 - e_3 * fs_63_22_11 * h8_m3 + e_3 * fs_483_286_105 * h8_m1 + e_3 * fs_114_11_2 * r_2 * h6_m3 + e_3 * fs_266_11_5 * r_2 * h6_m1 + e_3 * fs_765_22_6 * r_4 * h4_m3 + e_3 * fs_2295_286_42 * r_4 * h4_m1 - e_3 * fs_125_11_35 * r_6 * h2_m1 + e_4 * fs_483_4199_2002 * h10_m3 + e_4 * fs_966_4199_231 * h10_m1 + e_4 * fs_126_209_11 * r_2 * h8_m3 - e_4 * fs_966_2717_105 * r_2 * h8_m1 - e_4 * fs_228_187_2 * r_4 * h6_m3 - e_4 * fs_532_187_5 * r_4 * h6_m1 - e_4 * fs_34_11_6 * r_6 * h4_m3 - e_4 * fs_102_143_42 * r_6 * h4_m1 + e_4 * fs_125_143_35 * r_8 * h2_m1 - e_5 * fs_42_4199_2002 * r_2 * h10_m3 - e_5 * fs_84_4199_231 * r_2 * h10_m1 - e_5 * fs_6_209_11 * r_4 * h8_m3 + e_5 * fs_46_2717_105 * r_4 * h8_m1 + e_5 * fs_8_187_2 * r_6 * h6_m3 + e_5 * fs_56_561_5 * r_6 * h6_m1 + e_5 * fs_1_11_6 * r_8 * h4_m3 + e_5 * fs_3_143_42 * r_8 * h4_m1 - e_5 * fs_10_429_35 * r_10 * h2_m1;
    }

    // NOTE: the angular components are formed in 93 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph2_m2, ph4_m4, ph4_m2, ph6_m4, ph6_m2, ph8_m4, ph8_m2, ph10_m4, ph10_m2, ab_2, pc_52 : simd::cache_line_size())
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
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;

        const auto h2_m2 = ph2_m2[k];
        const auto h4_m4 = ph4_m4[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h6_m4 = ph6_m4[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h8_m4 = ph8_m4[k];
        const auto h8_m2 = ph8_m2[k];
        const auto h10_m4 = ph10_m4[k];
        const auto h10_m2 = ph10_m2[k];

        pc_52[k] = - e_0 * fs_1575_32_10 * h2_m2 - e_1 * fs_255_8_42 * h4_m4 - e_1 * fs_3315_16_6 * h4_m2 + e_1 * fs_1125_16_10 * r_2 * h2_m2 - e_2 * fs_475_44_210 * h6_m4 - e_2 * fs_665_44_7 * h6_m2 + e_2 * fs_255_11_42 * r_2 * h4_m4 + e_2 * fs_3315_22_6 * r_2 * h4_m2 - e_2 * fs_125_4_10 * r_4 * h2_m2 + e_3 * fs_105_286_2310 * h8_m4 + e_3 * fs_147_22_21 * h8_m2 + e_3 * fs_95_22_210 * r_2 * h6_m4 + e_3 * fs_133_22_7 * r_2 * h6_m2 - e_3 * fs_765_143_42 * r_4 * h4_m4 - e_3 * fs_765_22_6 * r_4 * h4_m2 + e_3 * fs_125_22_10 * r_6 * h2_m2 + e_4 * fs_966_4199_286 * h10_m4 + e_4 * fs_1932_4199_22 * h10_m2 - e_4 * fs_210_2717_2310 * r_2 * h8_m4 - e_4 * fs_294_209_21 * r_2 * h8_m2 - e_4 * fs_95_187_210 * r_4 * h6_m4 - e_4 * fs_133_187_7 * r_4 * h6_m2 + e_4 * fs_68_143_42 * r_6 * h4_m4 + e_4 * fs_34_11_6 * r_6 * h4_m2 - e_4 * fs_125_286_10 * r_8 * h2_m2 - e_5 * fs_84_4199_286 * r_2 * h10_m4 - e_5 * fs_168_4199_22 * r_2 * h10_m2 + e_5 * fs_10_2717_2310 * r_4 * h8_m4 + e_5 * fs_14_209_21 * r_4 * h8_m2 + e_5 * fs_10_561_210 * r_6 * h6_m4 + e_5 * fs_14_561_7 * r_6 * h6_m2 - e_5 * fs_2_143_42 * r_8 * h4_m4 - e_5 * fs_1_11_6 * r_8 * h4_m2 + e_5 * fs_5_429_10 * r_10 * h2_m2;
    }

    // NOTE: the angular components are formed in 93 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph4_m4, ph4_m3, ph6_m5, ph6_m4, ph6_m3, ph8_m5, ph8_m4, ph8_m3, ph10_m5, ph10_m4, ph10_m3, ab_2, pc_53, pc_54, pc_55 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

        const auto h4_m4 = ph4_m4[k];
        const auto h4_m3 = ph4_m3[k];
        const auto h6_m5 = ph6_m5[k];
        const auto h6_m4 = ph6_m4[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h8_m5 = ph8_m5[k];
        const auto h8_m4 = ph8_m4[k];
        const auto h8_m3 = ph8_m3[k];
        const auto h10_m5 = ph10_m5[k];
        const auto h10_m4 = ph10_m4[k];
        const auto h10_m3 = ph10_m3[k];

        pc_53[k] = e_1 * fs_255_8_42 * h4_m3 + e_2 * fs_95_44_2310 * h6_m5 + e_2 * fs_1995_44_14 * h6_m3 - e_2 * fs_255_11_42 * r_2 * h4_m3 + e_3 * fs_42_143_15015 * h8_m5 + e_3 * fs_441_143_77 * h8_m3 - e_3 * fs_19_22_2310 * r_2 * h6_m5 - e_3 * fs_399_22_14 * r_2 * h6_m3 + e_3 * fs_765_143_42 * r_4 * h4_m3 + e_4 * fs_483_8398_1430 * h10_m5 + e_4 * fs_483_8398_286 * h10_m3 - e_4 * fs_168_2717_15015 * r_2 * h8_m5 - e_4 * fs_1764_2717_77 * r_2 * h8_m3 + e_4 * fs_19_187_2310 * r_4 * h6_m5 + e_4 * fs_399_187_14 * r_4 * h6_m3 - e_4 * fs_68_143_42 * r_6 * h4_m3 - e_5 * fs_21_4199_1430 * r_2 * h10_m5 - e_5 * fs_21_4199_286 * r_2 * h10_m3 + e_5 * fs_8_2717_15015 * r_4 * h8_m5 + e_5 * fs_84_2717_77 * r_4 * h8_m3 - e_5 * fs_2_561_2310 * r_6 * h6_m5 - e_5 * fs_14_187_14 * r_6 * h6_m3 + e_5 * fs_2_143_42 * r_8 * h4_m3;

        pc_54[k] = e_1 * f_255_2 * h4_m4 + e_2 * fs_1995_22_5 * h6_m4 - e_2 * f_1020_11 * r_2 * h4_m4 + e_3 * fs_882_143_55 * h8_m4 - e_3 * fs_399_11_5 * r_2 * h6_m4 + e_3 * f_3060_143 * r_4 * h4_m4 + e_4 * fs_161_4199_3003 * h10_m4 - e_4 * fs_3528_2717_55 * r_2 * h8_m4 + e_4 * fs_798_187_5 * r_4 * h6_m4 - e_4 * f_272_143 * r_6 * h4_m4 - e_5 * fs_14_4199_3003 * r_2 * h10_m4 + e_5 * fs_168_2717_55 * r_4 * h8_m4 - e_5 * fs_28_187_5 * r_6 * h6_m4 + e_5 * f_8_143 * r_8 * h4_m4;

        pc_55[k] = - e_1 * f_4335_8 * h4_m3 - e_2 * fs_1995_22_3 * h6_m3 + e_2 * f_4335_11 * r_2 * h4_m3 + e_3 * fs_1323_286_66 * h8_m3 + e_3 * fs_399_11_3 * r_2 * h6_m3 - e_3 * f_13005_143 * r_4 * h4_m3 + e_4 * fs_322_4199_3003 * h10_m3 - e_4 * fs_2646_2717_66 * r_2 * h8_m3 - e_4 * fs_798_187_3 * r_4 * h6_m3 + e_4 * f_1156_143 * r_6 * h4_m3 - e_5 * fs_28_4199_3003 * r_2 * h10_m3 + e_5 * fs_126_2717_66 * r_4 * h8_m3 + e_5 * fs_28_187_3 * r_6 * h6_m3 - e_5 * f_34_143 * r_8 * h4_m3;
    }

    // NOTE: the angular components are formed in 93 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph2_m2, ph2_m1, ph4_m2, ph4_m1, ph6_m2, ph6_m1, ph8_m2, ph8_m1, ph10_m2, ph10_m1, ab_2, pc_56, pc_57 : simd::cache_line_size())
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
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;

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

        pc_56[k] = e_0 * fs_1575_16_15 * h2_m2 + e_1 * f_2805_4 * h4_m2 - e_1 * fs_1125_8_15 * r_2 * h2_m2 - e_2 * fs_95_4_42 * h6_m2 - e_2 * f_510_1 * r_2 * h4_m2 + e_2 * fs_125_2_15 * r_4 * h2_m2 + e_3 * fs_252_143_14 * h8_m2 + e_3 * fs_19_2_42 * r_2 * h6_m2 + e_3 * f_1530_13 * r_4 * h4_m2 - e_3 * fs_125_11_15 * r_6 * h2_m2 + e_4 * fs_4508_4199_33 * h10_m2 - e_4 * fs_1008_2717_14 * r_2 * h8_m2 - e_4 * fs_19_17_42 * r_4 * h6_m2 - e_4 * f_136_13 * r_6 * h4_m2 + e_4 * fs_125_143_15 * r_8 * h2_m2 - e_5 * fs_392_4199_33 * r_2 * h10_m2 + e_5 * fs_48_2717_14 * r_4 * h8_m2 + e_5 * fs_2_51_42 * r_6 * h6_m2 + e_5 * f_4_13 * r_8 * h4_m2 - e_5 * fs_10_429_15 * r_10 * h2_m2;

        pc_57[k] = - e_0 * fs_1575_8_30 * h2_m1 + e_1 * f_255_8 * h4_m1 + e_1 * fs_1125_4_30 * r_2 * h2_m1 + e_2 * fs_95_22_210 * h6_m1 - e_2 * f_255_11 * r_2 * h4_m1 - e_2 * fs_125_1_30 * r_4 * h2_m1 - e_3 * fs_2205_286_10 * h8_m1 - e_3 * fs_19_11_210 * r_2 * h6_m1 + e_3 * f_765_143 * r_4 * h4_m1 + e_3 * fs_250_11_30 * r_6 * h2_m1 + e_4 * fs_6762_4199_22 * h10_m1 + e_4 * fs_4410_2717_10 * r_2 * h8_m1 + e_4 * fs_38_187_210 * r_4 * h6_m1 - e_4 * f_68_143 * r_6 * h4_m1 - e_4 * fs_250_143_30 * r_8 * h2_m1 - e_5 * fs_588_4199_22 * r_2 * h10_m1 - e_5 * fs_210_2717_10 * r_4 * h8_m1 - e_5 * fs_4_561_210 * r_6 * h6_m1 + e_5 * f_2_143 * r_8 * h4_m1 + e_5 * fs_20_429_30 * r_10 * h2_m1;
    }

    // NOTE: the angular components are formed in 93 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph2_0, ph2_p1, ph4_0, ph4_p1, ph6_0, ph6_p1, ph8_0, ph8_p1, ph10_0, ph10_p1, ab_2, pc_58, pc_59 : simd::cache_line_size())
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
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;

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

        pc_58[k] = e_0 * f_23625_16 * h2_0 - e_1 * f_1275_2 * h4_0 - e_1 * f_16875_8 * r_2 * h2_0 + e_2 * f_1995_11 * h6_0 + e_2 * f_5100_11 * r_2 * h4_0 + e_2 * f_1875_2 * r_4 * h2_0 - e_3 * f_5292_143 * h8_0 - e_3 * f_798_11 * r_2 * h6_0 - e_3 * f_15300_143 * r_4 * h4_0 - e_3 * f_1875_11 * r_6 * h2_0 + e_4 * f_33810_4199 * h10_0 + e_4 * f_21168_2717 * r_2 * h8_0 + e_4 * f_1596_187 * r_4 * h6_0 + e_4 * f_1360_143 * r_6 * h4_0 + e_4 * f_1875_143 * r_8 * h2_0 - e_5 * f_2940_4199 * r_2 * h10_0 - e_5 * f_1008_2717 * r_4 * h8_0 - e_5 * f_56_187 * r_6 * h6_0 - e_5 * f_40_143 * r_8 * h4_0 - e_5 * f_50_143 * r_10 * h2_0;

        pc_59[k] = - e_0 * fs_1575_8_30 * h2_p1 + e_1 * f_255_8 * h4_p1 + e_1 * fs_1125_4_30 * r_2 * h2_p1 + e_2 * fs_95_22_210 * h6_p1 - e_2 * f_255_11 * r_2 * h4_p1 - e_2 * fs_125_1_30 * r_4 * h2_p1 - e_3 * fs_2205_286_10 * h8_p1 - e_3 * fs_19_11_210 * r_2 * h6_p1 + e_3 * f_765_143 * r_4 * h4_p1 + e_3 * fs_250_11_30 * r_6 * h2_p1 + e_4 * fs_6762_4199_22 * h10_p1 + e_4 * fs_4410_2717_10 * r_2 * h8_p1 + e_4 * fs_38_187_210 * r_4 * h6_p1 - e_4 * f_68_143 * r_6 * h4_p1 - e_4 * fs_250_143_30 * r_8 * h2_p1 - e_5 * fs_588_4199_22 * r_2 * h10_p1 - e_5 * fs_210_2717_10 * r_4 * h8_p1 - e_5 * fs_4_561_210 * r_6 * h6_p1 + e_5 * f_2_143 * r_8 * h4_p1 + e_5 * fs_20_429_30 * r_10 * h2_p1;
    }

    // NOTE: the angular components are formed in 93 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph2_p2, ph4_p2, ph4_p3, ph6_p2, ph6_p3, ph8_p2, ph8_p3, ph10_p2, ph10_p3, ab_2, pc_60, pc_61 : simd::cache_line_size())
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
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;

        const auto h2_p2 = ph2_p2[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h4_p3 = ph4_p3[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h8_p2 = ph8_p2[k];
        const auto h8_p3 = ph8_p3[k];
        const auto h10_p2 = ph10_p2[k];
        const auto h10_p3 = ph10_p3[k];

        pc_60[k] = e_0 * fs_1575_16_15 * h2_p2 + e_1 * f_2805_4 * h4_p2 - e_1 * fs_1125_8_15 * r_2 * h2_p2 - e_2 * fs_95_4_42 * h6_p2 - e_2 * f_510_1 * r_2 * h4_p2 + e_2 * fs_125_2_15 * r_4 * h2_p2 + e_3 * fs_252_143_14 * h8_p2 + e_3 * fs_19_2_42 * r_2 * h6_p2 + e_3 * f_1530_13 * r_4 * h4_p2 - e_3 * fs_125_11_15 * r_6 * h2_p2 + e_4 * fs_4508_4199_33 * h10_p2 - e_4 * fs_1008_2717_14 * r_2 * h8_p2 - e_4 * fs_19_17_42 * r_4 * h6_p2 - e_4 * f_136_13 * r_6 * h4_p2 + e_4 * fs_125_143_15 * r_8 * h2_p2 - e_5 * fs_392_4199_33 * r_2 * h10_p2 + e_5 * fs_48_2717_14 * r_4 * h8_p2 + e_5 * fs_2_51_42 * r_6 * h6_p2 + e_5 * f_4_13 * r_8 * h4_p2 - e_5 * fs_10_429_15 * r_10 * h2_p2;

        pc_61[k] = - e_1 * f_4335_8 * h4_p3 - e_2 * fs_1995_22_3 * h6_p3 + e_2 * f_4335_11 * r_2 * h4_p3 + e_3 * fs_1323_286_66 * h8_p3 + e_3 * fs_399_11_3 * r_2 * h6_p3 - e_3 * f_13005_143 * r_4 * h4_p3 + e_4 * fs_322_4199_3003 * h10_p3 - e_4 * fs_2646_2717_66 * r_2 * h8_p3 - e_4 * fs_798_187_3 * r_4 * h6_p3 + e_4 * f_1156_143 * r_6 * h4_p3 - e_5 * fs_28_4199_3003 * r_2 * h10_p3 + e_5 * fs_126_2717_66 * r_4 * h8_p3 + e_5 * fs_28_187_3 * r_6 * h6_p3 - e_5 * f_34_143 * r_8 * h4_p3;
    }

    // NOTE: the angular components are formed in 93 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph4_m3, ph4_p4, ph6_m5, ph6_m3, ph6_p4, ph8_m5, ph8_m3, ph8_p4, ph10_m5, ph10_m3, ph10_p4, ab_2, pc_62, pc_63 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

        const auto h4_m3 = ph4_m3[k];
        const auto h4_p4 = ph4_p4[k];
        const auto h6_m5 = ph6_m5[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h6_p4 = ph6_p4[k];
        const auto h8_m5 = ph8_m5[k];
        const auto h8_m3 = ph8_m3[k];
        const auto h8_p4 = ph8_p4[k];
        const auto h10_m5 = ph10_m5[k];
        const auto h10_m3 = ph10_m3[k];
        const auto h10_p4 = ph10_p4[k];

        pc_62[k] = e_1 * f_255_2 * h4_p4 + e_2 * fs_1995_22_5 * h6_p4 - e_2 * f_1020_11 * r_2 * h4_p4 + e_3 * fs_882_143_55 * h8_p4 - e_3 * fs_399_11_5 * r_2 * h6_p4 + e_3 * f_3060_143 * r_4 * h4_p4 + e_4 * fs_161_4199_3003 * h10_p4 - e_4 * fs_3528_2717_55 * r_2 * h8_p4 + e_4 * fs_798_187_5 * r_4 * h6_p4 - e_4 * f_272_143 * r_6 * h4_p4 - e_5 * fs_14_4199_3003 * r_2 * h10_p4 + e_5 * fs_168_2717_55 * r_4 * h8_p4 - e_5 * fs_28_187_5 * r_6 * h6_p4 + e_5 * f_8_143 * r_8 * h4_p4;

        pc_63[k] = - e_1 * fs_255_8_42 * h4_m3 + e_2 * fs_95_44_2310 * h6_m5 - e_2 * fs_1995_44_14 * h6_m3 + e_2 * fs_255_11_42 * r_2 * h4_m3 + e_3 * fs_42_143_15015 * h8_m5 - e_3 * fs_441_143_77 * h8_m3 - e_3 * fs_19_22_2310 * r_2 * h6_m5 + e_3 * fs_399_22_14 * r_2 * h6_m3 - e_3 * fs_765_143_42 * r_4 * h4_m3 + e_4 * fs_483_8398_1430 * h10_m5 - e_4 * fs_483_8398_286 * h10_m3 - e_4 * fs_168_2717_15015 * r_2 * h8_m5 + e_4 * fs_1764_2717_77 * r_2 * h8_m3 + e_4 * fs_19_187_2310 * r_4 * h6_m5 - e_4 * fs_399_187_14 * r_4 * h6_m3 + e_4 * fs_68_143_42 * r_6 * h4_m3 - e_5 * fs_21_4199_1430 * r_2 * h10_m5 + e_5 * fs_21_4199_286 * r_2 * h10_m3 + e_5 * fs_8_2717_15015 * r_4 * h8_m5 - e_5 * fs_84_2717_77 * r_4 * h8_m3 - e_5 * fs_2_561_2310 * r_6 * h6_m5 + e_5 * fs_14_187_14 * r_6 * h6_m3 - e_5 * fs_2_143_42 * r_8 * h4_m3;
    }

    // NOTE: the angular components are formed in 93 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph2_m2, ph4_m4, ph4_m2, ph6_m4, ph6_m2, ph8_m4, ph8_m2, ph10_m4, ph10_m2, ab_2, pc_64 : simd::cache_line_size())
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
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;

        const auto h2_m2 = ph2_m2[k];
        const auto h4_m4 = ph4_m4[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h6_m4 = ph6_m4[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h8_m4 = ph8_m4[k];
        const auto h8_m2 = ph8_m2[k];
        const auto h10_m4 = ph10_m4[k];
        const auto h10_m2 = ph10_m2[k];

        pc_64[k] = e_0 * fs_1575_32_10 * h2_m2 - e_1 * fs_255_8_42 * h4_m4 + e_1 * fs_3315_16_6 * h4_m2 - e_1 * fs_1125_16_10 * r_2 * h2_m2 - e_2 * fs_475_44_210 * h6_m4 + e_2 * fs_665_44_7 * h6_m2 + e_2 * fs_255_11_42 * r_2 * h4_m4 - e_2 * fs_3315_22_6 * r_2 * h4_m2 + e_2 * fs_125_4_10 * r_4 * h2_m2 + e_3 * fs_105_286_2310 * h8_m4 - e_3 * fs_147_22_21 * h8_m2 + e_3 * fs_95_22_210 * r_2 * h6_m4 - e_3 * fs_133_22_7 * r_2 * h6_m2 - e_3 * fs_765_143_42 * r_4 * h4_m4 + e_3 * fs_765_22_6 * r_4 * h4_m2 - e_3 * fs_125_22_10 * r_6 * h2_m2 + e_4 * fs_966_4199_286 * h10_m4 - e_4 * fs_1932_4199_22 * h10_m2 - e_4 * fs_210_2717_2310 * r_2 * h8_m4 + e_4 * fs_294_209_21 * r_2 * h8_m2 - e_4 * fs_95_187_210 * r_4 * h6_m4 + e_4 * fs_133_187_7 * r_4 * h6_m2 + e_4 * fs_68_143_42 * r_6 * h4_m4 - e_4 * fs_34_11_6 * r_6 * h4_m2 + e_4 * fs_125_286_10 * r_8 * h2_m2 - e_5 * fs_84_4199_286 * r_2 * h10_m4 + e_5 * fs_168_4199_22 * r_2 * h10_m2 + e_5 * fs_10_2717_2310 * r_4 * h8_m4 - e_5 * fs_14_209_21 * r_4 * h8_m2 + e_5 * fs_10_561_210 * r_6 * h6_m4 - e_5 * fs_14_561_7 * r_6 * h6_m2 - e_5 * fs_2_143_42 * r_8 * h4_m4 + e_5 * fs_1_11_6 * r_8 * h4_m2 - e_5 * fs_5_429_10 * r_10 * h2_m2;
    }

    // NOTE: the angular components are formed in 93 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph2_m1, ph4_m3, ph4_m1, ph6_m3, ph6_m1, ph8_m3, ph8_m1, ph10_m3, ph10_m1, ab_2, pc_65 : simd::cache_line_size())
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
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;

        const auto h2_m1 = ph2_m1[k];
        const auto h4_m3 = ph4_m3[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h8_m3 = ph8_m3[k];
        const auto h8_m1 = ph8_m1[k];
        const auto h10_m3 = ph10_m3[k];
        const auto h10_m1 = ph10_m1[k];

        pc_65[k] = - e_0 * fs_1575_16_35 * h2_m1 + e_1 * fs_3315_16_6 * h4_m3 - e_1 * fs_765_16_42 * h4_m1 + e_1 * fs_1125_8_35 * r_2 * h2_m1 - e_2 * fs_285_11_2 * h6_m3 + e_2 * fs_665_11_5 * h6_m1 - e_2 * fs_3315_22_6 * r_2 * h4_m3 + e_2 * fs_765_22_42 * r_2 * h4_m1 - e_2 * fs_125_2_35 * r_4 * h2_m1 - e_3 * fs_63_22_11 * h8_m3 - e_3 * fs_483_286_105 * h8_m1 + e_3 * fs_114_11_2 * r_2 * h6_m3 - e_3 * fs_266_11_5 * r_2 * h6_m1 + e_3 * fs_765_22_6 * r_4 * h4_m3 - e_3 * fs_2295_286_42 * r_4 * h4_m1 + e_3 * fs_125_11_35 * r_6 * h2_m1 + e_4 * fs_483_4199_2002 * h10_m3 - e_4 * fs_966_4199_231 * h10_m1 + e_4 * fs_126_209_11 * r_2 * h8_m3 + e_4 * fs_966_2717_105 * r_2 * h8_m1 - e_4 * fs_228_187_2 * r_4 * h6_m3 + e_4 * fs_532_187_5 * r_4 * h6_m1 - e_4 * fs_34_11_6 * r_6 * h4_m3 + e_4 * fs_102_143_42 * r_6 * h4_m1 - e_4 * fs_125_143_35 * r_8 * h2_m1 - e_5 * fs_42_4199_2002 * r_2 * h10_m3 + e_5 * fs_84_4199_231 * r_2 * h10_m1 - e_5 * fs_6_209_11 * r_4 * h8_m3 - e_5 * fs_46_2717_105 * r_4 * h8_m1 + e_5 * fs_8_187_2 * r_6 * h6_m3 - e_5 * fs_56_561_5 * r_6 * h6_m1 + e_5 * fs_1_11_6 * r_8 * h4_m3 - e_5 * fs_3_143_42 * r_8 * h4_m1 + e_5 * fs_10_429_35 * r_10 * h2_m1;
    }

    // NOTE: the angular components are formed in 93 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph2_m2, ph2_p1, ph4_m2, ph4_p1, ph6_m2, ph6_p1, ph8_m2, ph8_p1, ph10_m2, ph10_p1, ab_2, pc_66, pc_67 : simd::cache_line_size())
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
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;

        const auto h2_m2 = ph2_m2[k];
        const auto h2_p1 = ph2_p1[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h8_m2 = ph8_m2[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h10_m2 = ph10_m2[k];
        const auto h10_p1 = ph10_p1[k];

        pc_66[k] = - e_0 * fs_1575_32_70 * h2_m2 - e_1 * fs_765_16_42 * h4_m2 + e_1 * fs_1125_16_70 * r_2 * h2_m2 + e_2 * f_5035_44 * h6_m2 + e_2 * fs_765_22_42 * r_2 * h4_m2 - e_2 * fs_125_4_70 * r_4 * h2_m2 - e_3 * fs_4179_286_3 * h8_m2 - e_3 * f_1007_22 * r_2 * h6_m2 - e_3 * fs_2295_286_42 * r_4 * h4_m2 + e_3 * fs_125_22_70 * r_6 * h2_m2 + e_4 * fs_1932_4199_154 * h10_m2 + e_4 * fs_8358_2717_3 * r_2 * h8_m2 + e_4 * f_1007_187 * r_4 * h6_m2 + e_4 * fs_102_143_42 * r_6 * h4_m2 - e_4 * fs_125_286_70 * r_8 * h2_m2 - e_5 * fs_168_4199_154 * r_2 * h10_m2 - e_5 * fs_398_2717_3 * r_4 * h8_m2 - e_5 * f_106_561 * r_6 * h6_m2 - e_5 * fs_3_143_42 * r_8 * h4_m2 + e_5 * fs_5_429_70 * r_10 * h2_m2;

        pc_67[k] = e_0 * fs_7875_16_7 * h2_p1 - e_1 * fs_255_8_210 * h4_p1 - e_1 * fs_5625_8_7 * r_2 * h2_p1 + e_2 * f_1520_11 * h6_p1 + e_2 * fs_255_11_210 * r_2 * h4_p1 + e_2 * fs_625_2_7 * r_4 * h2_p1 - e_3 * fs_987_143_21 * h8_p1 - e_3 * f_608_11 * r_2 * h6_p1 - e_3 * fs_765_143_210 * r_4 * h4_p1 - e_3 * fs_625_11_7 * r_6 * h2_p1 + e_4 * fs_966_4199_1155 * h10_p1 + e_4 * fs_3948_2717_21 * r_2 * h8_p1 + e_4 * f_1216_187 * r_4 * h6_p1 + e_4 * fs_68_143_210 * r_6 * h4_p1 + e_4 * fs_625_143_7 * r_8 * h2_p1 - e_5 * fs_84_4199_1155 * r_2 * h10_p1 - e_5 * fs_188_2717_21 * r_4 * h8_p1 - e_5 * f_128_561 * r_6 * h6_p1 - e_5 * fs_2_143_210 * r_8 * h4_p1 - e_5 * fs_50_429_7 * r_10 * h2_p1;
    }

    // NOTE: the angular components are formed in 93 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph2_0, ph2_p2, ph4_0, ph4_p2, ph6_0, ph6_p2, ph8_0, ph8_p2, ph10_0, ph10_p2, ab_2, pc_68 : simd::cache_line_size())
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
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;

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

        pc_68[k] = e_0 * fs_1575_16_210 * h2_0 - e_0 * fs_1575_32_70 * h2_p2 - e_1 * fs_255_8_210 * h4_0 - e_1 * fs_765_16_42 * h4_p2 - e_1 * fs_1125_8_210 * r_2 * h2_0 + e_1 * fs_1125_16_70 * r_2 * h2_p2 + e_2 * fs_95_22_210 * h6_0 + e_2 * f_5035_44 * h6_p2 + e_2 * fs_255_11_210 * r_2 * h4_0 + e_2 * fs_765_22_42 * r_2 * h4_p2 + e_2 * fs_125_2_210 * r_4 * h2_0 - e_2 * fs_125_4_70 * r_4 * h2_p2 + e_3 * fs_63_143_210 * h8_0 - e_3 * fs_4179_286_3 * h8_p2 - e_3 * fs_19_11_210 * r_2 * h6_0 - e_3 * f_1007_22 * r_2 * h6_p2 - e_3 * fs_765_143_210 * r_4 * h4_0 - e_3 * fs_2295_286_42 * r_4 * h4_p2 - e_3 * fs_125_11_210 * r_6 * h2_0 + e_3 * fs_125_22_70 * r_6 * h2_p2 - e_4 * fs_1932_4199_210 * h10_0 + e_4 * fs_1932_4199_154 * h10_p2 - e_4 * fs_252_2717_210 * r_2 * h8_0 + e_4 * fs_8358_2717_3 * r_2 * h8_p2 + e_4 * fs_38_187_210 * r_4 * h6_0 + e_4 * f_1007_187 * r_4 * h6_p2 + e_4 * fs_68_143_210 * r_6 * h4_0 + e_4 * fs_102_143_42 * r_6 * h4_p2 + e_4 * fs_125_143_210 * r_8 * h2_0 - e_4 * fs_125_286_70 * r_8 * h2_p2 + e_5 * fs_168_4199_210 * r_2 * h10_0 - e_5 * fs_168_4199_154 * r_2 * h10_p2 + e_5 * fs_12_2717_210 * r_4 * h8_0 - e_5 * fs_398_2717_3 * r_4 * h8_p2 - e_5 * fs_4_561_210 * r_6 * h6_0 - e_5 * f_106_561 * r_6 * h6_p2 - e_5 * fs_2_143_210 * r_8 * h4_0 - e_5 * fs_3_143_42 * r_8 * h4_p2 - e_5 * fs_10_429_210 * r_10 * h2_0 + e_5 * fs_5_429_70 * r_10 * h2_p2;
    }

    // NOTE: the angular components are formed in 93 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph2_p1, ph4_p1, ph4_p3, ph6_p1, ph6_p3, ph8_p1, ph8_p3, ph10_p1, ph10_p3, ab_2, pc_69 : simd::cache_line_size())
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
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;

        const auto h2_p1 = ph2_p1[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h4_p3 = ph4_p3[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h8_p3 = ph8_p3[k];
        const auto h10_p1 = ph10_p1[k];
        const auto h10_p3 = ph10_p3[k];

        pc_69[k] = - e_0 * fs_1575_16_35 * h2_p1 - e_1 * fs_765_16_42 * h4_p1 + e_1 * fs_3315_16_6 * h4_p3 + e_1 * fs_1125_8_35 * r_2 * h2_p1 + e_2 * fs_665_11_5 * h6_p1 - e_2 * fs_285_11_2 * h6_p3 + e_2 * fs_765_22_42 * r_2 * h4_p1 - e_2 * fs_3315_22_6 * r_2 * h4_p3 - e_2 * fs_125_2_35 * r_4 * h2_p1 - e_3 * fs_483_286_105 * h8_p1 - e_3 * fs_63_22_11 * h8_p3 - e_3 * fs_266_11_5 * r_2 * h6_p1 + e_3 * fs_114_11_2 * r_2 * h6_p3 - e_3 * fs_2295_286_42 * r_4 * h4_p1 + e_3 * fs_765_22_6 * r_4 * h4_p3 + e_3 * fs_125_11_35 * r_6 * h2_p1 - e_4 * fs_966_4199_231 * h10_p1 + e_4 * fs_483_4199_2002 * h10_p3 + e_4 * fs_966_2717_105 * r_2 * h8_p1 + e_4 * fs_126_209_11 * r_2 * h8_p3 + e_4 * fs_532_187_5 * r_4 * h6_p1 - e_4 * fs_228_187_2 * r_4 * h6_p3 + e_4 * fs_102_143_42 * r_6 * h4_p1 - e_4 * fs_34_11_6 * r_6 * h4_p3 - e_4 * fs_125_143_35 * r_8 * h2_p1 + e_5 * fs_84_4199_231 * r_2 * h10_p1 - e_5 * fs_42_4199_2002 * r_2 * h10_p3 - e_5 * fs_46_2717_105 * r_4 * h8_p1 - e_5 * fs_6_209_11 * r_4 * h8_p3 - e_5 * fs_56_561_5 * r_6 * h6_p1 + e_5 * fs_8_187_2 * r_6 * h6_p3 - e_5 * fs_3_143_42 * r_8 * h4_p1 + e_5 * fs_1_11_6 * r_8 * h4_p3 + e_5 * fs_10_429_35 * r_10 * h2_p1;
    }

    // NOTE: the angular components are formed in 93 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph2_p2, ph4_p2, ph4_p4, ph6_p2, ph6_p4, ph8_p2, ph8_p4, ph10_p2, ph10_p4, ab_2, pc_70 : simd::cache_line_size())
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
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;

        const auto h2_p2 = ph2_p2[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h4_p4 = ph4_p4[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h6_p4 = ph6_p4[k];
        const auto h8_p2 = ph8_p2[k];
        const auto h8_p4 = ph8_p4[k];
        const auto h10_p2 = ph10_p2[k];
        const auto h10_p4 = ph10_p4[k];

        pc_70[k] = e_0 * fs_1575_32_10 * h2_p2 + e_1 * fs_3315_16_6 * h4_p2 - e_1 * fs_255_8_42 * h4_p4 - e_1 * fs_1125_16_10 * r_2 * h2_p2 + e_2 * fs_665_44_7 * h6_p2 - e_2 * fs_475_44_210 * h6_p4 - e_2 * fs_3315_22_6 * r_2 * h4_p2 + e_2 * fs_255_11_42 * r_2 * h4_p4 + e_2 * fs_125_4_10 * r_4 * h2_p2 - e_3 * fs_147_22_21 * h8_p2 + e_3 * fs_105_286_2310 * h8_p4 - e_3 * fs_133_22_7 * r_2 * h6_p2 + e_3 * fs_95_22_210 * r_2 * h6_p4 + e_3 * fs_765_22_6 * r_4 * h4_p2 - e_3 * fs_765_143_42 * r_4 * h4_p4 - e_3 * fs_125_22_10 * r_6 * h2_p2 - e_4 * fs_1932_4199_22 * h10_p2 + e_4 * fs_966_4199_286 * h10_p4 + e_4 * fs_294_209_21 * r_2 * h8_p2 - e_4 * fs_210_2717_2310 * r_2 * h8_p4 + e_4 * fs_133_187_7 * r_4 * h6_p2 - e_4 * fs_95_187_210 * r_4 * h6_p4 - e_4 * fs_34_11_6 * r_6 * h4_p2 + e_4 * fs_68_143_42 * r_6 * h4_p4 + e_4 * fs_125_286_10 * r_8 * h2_p2 + e_5 * fs_168_4199_22 * r_2 * h10_p2 - e_5 * fs_84_4199_286 * r_2 * h10_p4 - e_5 * fs_14_209_21 * r_4 * h8_p2 + e_5 * fs_10_2717_2310 * r_4 * h8_p4 - e_5 * fs_14_561_7 * r_6 * h6_p2 + e_5 * fs_10_561_210 * r_6 * h6_p4 + e_5 * fs_1_11_6 * r_8 * h4_p2 - e_5 * fs_2_143_42 * r_8 * h4_p4 - e_5 * fs_5_429_10 * r_10 * h2_p2;
    }

    // NOTE: the angular components are formed in 93 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph4_p3, ph6_p3, ph6_p5, ph8_p3, ph8_p5, ph10_p3, ph10_p5, ab_2, pc_71 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

        const auto h4_p3 = ph4_p3[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h6_p5 = ph6_p5[k];
        const auto h8_p3 = ph8_p3[k];
        const auto h8_p5 = ph8_p5[k];
        const auto h10_p3 = ph10_p3[k];
        const auto h10_p5 = ph10_p5[k];

        pc_71[k] = - e_1 * fs_255_8_42 * h4_p3 - e_2 * fs_1995_44_14 * h6_p3 + e_2 * fs_95_44_2310 * h6_p5 + e_2 * fs_255_11_42 * r_2 * h4_p3 - e_3 * fs_441_143_77 * h8_p3 + e_3 * fs_42_143_15015 * h8_p5 + e_3 * fs_399_22_14 * r_2 * h6_p3 - e_3 * fs_19_22_2310 * r_2 * h6_p5 - e_3 * fs_765_143_42 * r_4 * h4_p3 - e_4 * fs_483_8398_286 * h10_p3 + e_4 * fs_483_8398_1430 * h10_p5 + e_4 * fs_1764_2717_77 * r_2 * h8_p3 - e_4 * fs_168_2717_15015 * r_2 * h8_p5 - e_4 * fs_399_187_14 * r_4 * h6_p3 + e_4 * fs_19_187_2310 * r_4 * h6_p5 + e_4 * fs_68_143_42 * r_6 * h4_p3 + e_5 * fs_21_4199_286 * r_2 * h10_p3 - e_5 * fs_21_4199_1430 * r_2 * h10_p5 - e_5 * fs_84_2717_77 * r_4 * h8_p3 + e_5 * fs_8_2717_15015 * r_4 * h8_p5 + e_5 * fs_14_187_14 * r_6 * h6_p3 - e_5 * fs_2_561_2310 * r_6 * h6_p5 - e_5 * fs_2_143_42 * r_8 * h4_p3;
    }

    // NOTE: the angular components are formed in 93 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph2_m2, ph4_m2, ph6_m6, ph6_m2, ph8_m6, ph8_m2, ph10_m6, ph10_m2, ab_2, pc_72 : simd::cache_line_size())
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
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;

        const auto h2_m2 = ph2_m2[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h6_m6 = ph6_m6[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h8_m6 = ph8_m6[k];
        const auto h8_m2 = ph8_m2[k];
        const auto h10_m6 = ph10_m6[k];
        const auto h10_m2 = ph10_m2[k];

        pc_72[k] = e_0 * fs_1575_32_2 * h2_m2 + e_1 * fs_255_4_30 * h4_m2 - e_1 * fs_1125_16_2 * r_2 * h2_m2 + e_2 * fs_285_44_77 * h6_m6 + e_2 * fs_665_22_35 * h6_m2 - e_2 * fs_510_11_30 * r_2 * h4_m2 + e_2 * fs_125_4_2 * r_4 * h2_m2 + e_3 * fs_441_143_143 * h8_m6 + e_3 * fs_294_143_105 * h8_m2 - e_3 * fs_57_22_77 * r_2 * h6_m6 - e_3 * fs_133_11_35 * r_2 * h6_m2 + e_3 * fs_1530_143_30 * r_4 * h4_m2 - e_3 * fs_125_22_2 * r_6 * h2_m2 + e_4 * fs_483_4199_715 * h10_m6 + e_4 * fs_483_8398_110 * h10_m2 - e_4 * fs_1764_2717_143 * r_2 * h8_m6 - e_4 * fs_1176_2717_105 * r_2 * h8_m2 + e_4 * fs_57_187_77 * r_4 * h6_m6 + e_4 * fs_266_187_35 * r_4 * h6_m2 - e_4 * fs_136_143_30 * r_6 * h4_m2 + e_4 * fs_125_286_2 * r_8 * h2_m2 - e_5 * fs_42_4199_715 * r_2 * h10_m6 - e_5 * fs_21_4199_110 * r_2 * h10_m2 + e_5 * fs_84_2717_143 * r_4 * h8_m6 + e_5 * fs_56_2717_105 * r_4 * h8_m2 - e_5 * fs_2_187_77 * r_6 * h6_m6 - e_5 * fs_28_561_35 * r_6 * h6_m2 + e_5 * fs_4_143_30 * r_8 * h4_m2 - e_5 * fs_5_429_2 * r_10 * h2_m2;
    }

    // NOTE: the angular components are formed in 93 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph2_m1, ph4_m1, ph6_m5, ph6_m1, ph8_m5, ph8_m1, ph10_m5, ph10_m1, ab_2, pc_73 : simd::cache_line_size())
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
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;

        const auto h2_m1 = ph2_m1[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h6_m5 = ph6_m5[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h8_m5 = ph8_m5[k];
        const auto h8_m1 = ph8_m1[k];
        const auto h10_m5 = ph10_m5[k];
        const auto h10_m1 = ph10_m1[k];

        pc_73[k] = - e_0 * f_1575_4 * h2_m1 - e_1 * fs_765_8_30 * h4_m1 + e_1 * f_1125_2 * r_2 * h2_m1 - e_2 * fs_665_88_462 * h6_m5 + e_2 * fs_665_44_7 * h6_m1 + e_2 * fs_765_11_30 * r_2 * h4_m1 - e_2 * f_250_1 * r_4 * h2_m1 + e_3 * fs_21_286_3003 * h8_m5 + e_3 * fs_2499_143_3 * h8_m1 + e_3 * fs_133_44_462 * r_2 * h6_m5 - e_3 * fs_133_22_7 * r_2 * h6_m1 - e_3 * fs_2295_143_30 * r_4 * h4_m1 + e_3 * f_500_11 * r_6 * h2_m1 + e_4 * fs_2415_8398_286 * h10_m5 + e_4 * fs_483_4199_165 * h10_m1 - e_4 * fs_42_2717_3003 * r_2 * h8_m5 - e_4 * fs_9996_2717_3 * r_2 * h8_m1 - e_4 * fs_133_374_462 * r_4 * h6_m5 + e_4 * fs_133_187_7 * r_4 * h6_m1 + e_4 * fs_204_143_30 * r_6 * h4_m1 - e_4 * f_500_143 * r_8 * h2_m1 - e_5 * fs_105_4199_286 * r_2 * h10_m5 - e_5 * fs_42_4199_165 * r_2 * h10_m1 + e_5 * fs_2_2717_3003 * r_4 * h8_m5 + e_5 * fs_476_2717_3 * r_4 * h8_m1 + e_5 * fs_7_561_462 * r_6 * h6_m5 - e_5 * fs_14_561_7 * r_6 * h6_m1 - e_5 * fs_6_143_30 * r_8 * h4_m1 + e_5 * f_40_429 * r_10 * h2_m1;
    }

    // NOTE: the angular components are formed in 93 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph4_m4, ph6_m4, ph8_m4, ph10_m4, ab_2, pc_74 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

        const auto h4_m4 = ph4_m4[k];
        const auto h6_m4 = ph6_m4[k];
        const auto h8_m4 = ph8_m4[k];
        const auto h10_m4 = ph10_m4[k];

        pc_74[k] = e_1 * fs_255_4_30 * h4_m4 + e_2 * fs_2185_88_6 * h6_m4 - e_2 * fs_510_11_30 * r_2 * h4_m4 - e_3 * fs_777_286_66 * h8_m4 - e_3 * fs_437_44_6 * r_2 * h6_m4 + e_3 * fs_1530_143_30 * r_4 * h4_m4 + e_4 * fs_483_8398_10010 * h10_m4 + e_4 * fs_1554_2717_66 * r_2 * h8_m4 + e_4 * fs_437_374_6 * r_4 * h6_m4 - e_4 * fs_136_143_30 * r_6 * h4_m4 - e_5 * fs_21_4199_10010 * r_2 * h10_m4 - e_5 * fs_74_2717_66 * r_4 * h8_m4 - e_5 * fs_23_561_6 * r_6 * h6_m4 + e_5 * fs_4_143_30 * r_8 * h4_m4;
    }

    // NOTE: the angular components are formed in 93 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph2_m1, ph4_m3, ph4_m1, ph6_m3, ph6_m1, ph8_m3, ph8_m1, ph10_m3, ph10_m1, ab_2, pc_75 : simd::cache_line_size())
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
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;

        const auto h2_m1 = ph2_m1[k];
        const auto h4_m3 = ph4_m3[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h8_m3 = ph8_m3[k];
        const auto h8_m1 = ph8_m1[k];
        const auto h10_m3 = ph10_m3[k];
        const auto h10_m1 = ph10_m1[k];

        pc_75[k] = - e_0 * fs_1575_4_7 * h2_m1 - e_1 * fs_765_8_30 * h4_m3 + e_1 * fs_255_8_210 * h4_m1 + e_1 * fs_1125_2_7 * r_2 * h2_m1 + e_2 * fs_3705_88_10 * h6_m3 - e_2 * f_5035_44 * h6_m1 + e_2 * fs_765_11_30 * r_2 * h4_m3 - e_2 * fs_255_11_210 * r_2 * h4_m1 - e_2 * fs_250_1_7 * r_4 * h2_m1 - e_3 * fs_504_143_55 * h8_m3 + e_3 * fs_609_286_21 * h8_m1 - e_3 * fs_741_44_10 * r_2 * h6_m3 + e_3 * f_1007_22 * r_2 * h6_m1 - e_3 * fs_2295_143_30 * r_4 * h4_m3 + e_3 * fs_765_143_210 * r_4 * h4_m1 + e_3 * fs_500_11_7 * r_6 * h2_m1 + e_4 * fs_483_8398_10010 * h10_m3 + e_4 * fs_483_4199_1155 * h10_m1 + e_4 * fs_2016_2717_55 * r_2 * h8_m3 - e_4 * fs_1218_2717_21 * r_2 * h8_m1 + e_4 * fs_741_374_10 * r_4 * h6_m3 - e_4 * f_1007_187 * r_4 * h6_m1 + e_4 * fs_204_143_30 * r_6 * h4_m3 - e_4 * fs_68_143_210 * r_6 * h4_m1 - e_4 * fs_500_143_7 * r_8 * h2_m1 - e_5 * fs_21_4199_10010 * r_2 * h10_m3 - e_5 * fs_42_4199_1155 * r_2 * h10_m1 - e_5 * fs_96_2717_55 * r_4 * h8_m3 + e_5 * fs_58_2717_21 * r_4 * h8_m1 - e_5 * fs_13_187_10 * r_6 * h6_m3 + e_5 * f_106_561 * r_6 * h6_m1 - e_5 * fs_6_143_30 * r_8 * h4_m3 + e_5 * fs_2_143_210 * r_8 * h4_m1 + e_5 * fs_40_429_7 * r_10 * h2_m1;
    }

    // NOTE: the angular components are formed in 93 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph2_p2, ph6_p2, ph8_p2, ph10_p2, ab_2, pc_76 : simd::cache_line_size())
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
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;

        const auto h2_p2 = ph2_p2[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h8_p2 = ph8_p2[k];
        const auto h10_p2 = ph10_p2[k];

        pc_76[k] = e_0 * fs_1575_16_70 * h2_p2 - e_1 * fs_1125_8_70 * r_2 * h2_p2 + e_2 * f_95_4 * h6_p2 + e_2 * fs_125_2_70 * r_4 * h2_p2 - e_3 * fs_105_11_3 * h8_p2 - e_3 * f_19_2 * r_2 * h6_p2 - e_3 * fs_125_11_70 * r_6 * h2_p2 + e_4 * fs_2415_4199_154 * h10_p2 + e_4 * fs_420_209_3 * r_2 * h8_p2 + e_4 * f_19_17 * r_4 * h6_p2 + e_4 * fs_125_143_70 * r_8 * h2_p2 - e_5 * fs_210_4199_154 * r_2 * h10_p2 - e_5 * fs_20_209_3 * r_4 * h8_p2 - e_5 * f_2_51 * r_6 * h6_p2 - e_5 * fs_10_429_70 * r_10 * h2_p2;
    }

    // NOTE: the angular components are formed in 93 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph2_p1, ph4_p1, ph4_p3, ph6_p1, ph6_p3, ph8_p1, ph8_p3, ph10_p1, ph10_p3, ab_2, pc_77 : simd::cache_line_size())
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
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;

        const auto h2_p1 = ph2_p1[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h4_p3 = ph4_p3[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h8_p3 = ph8_p3[k];
        const auto h10_p1 = ph10_p1[k];
        const auto h10_p3 = ph10_p3[k];

        pc_77[k] = e_0 * fs_1575_4_7 * h2_p1 - e_1 * fs_255_8_210 * h4_p1 - e_1 * fs_765_8_30 * h4_p3 - e_1 * fs_1125_2_7 * r_2 * h2_p1 + e_2 * f_5035_44 * h6_p1 + e_2 * fs_3705_88_10 * h6_p3 + e_2 * fs_255_11_210 * r_2 * h4_p1 + e_2 * fs_765_11_30 * r_2 * h4_p3 + e_2 * fs_250_1_7 * r_4 * h2_p1 - e_3 * fs_609_286_21 * h8_p1 - e_3 * fs_504_143_55 * h8_p3 - e_3 * f_1007_22 * r_2 * h6_p1 - e_3 * fs_741_44_10 * r_2 * h6_p3 - e_3 * fs_765_143_210 * r_4 * h4_p1 - e_3 * fs_2295_143_30 * r_4 * h4_p3 - e_3 * fs_500_11_7 * r_6 * h2_p1 - e_4 * fs_483_4199_1155 * h10_p1 + e_4 * fs_483_8398_10010 * h10_p3 + e_4 * fs_1218_2717_21 * r_2 * h8_p1 + e_4 * fs_2016_2717_55 * r_2 * h8_p3 + e_4 * f_1007_187 * r_4 * h6_p1 + e_4 * fs_741_374_10 * r_4 * h6_p3 + e_4 * fs_68_143_210 * r_6 * h4_p1 + e_4 * fs_204_143_30 * r_6 * h4_p3 + e_4 * fs_500_143_7 * r_8 * h2_p1 + e_5 * fs_42_4199_1155 * r_2 * h10_p1 - e_5 * fs_21_4199_10010 * r_2 * h10_p3 - e_5 * fs_58_2717_21 * r_4 * h8_p1 - e_5 * fs_96_2717_55 * r_4 * h8_p3 - e_5 * f_106_561 * r_6 * h6_p1 - e_5 * fs_13_187_10 * r_6 * h6_p3 - e_5 * fs_2_143_210 * r_8 * h4_p1 - e_5 * fs_6_143_30 * r_8 * h4_p3 - e_5 * fs_40_429_7 * r_10 * h2_p1;
    }

    // NOTE: the angular components are formed in 93 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph2_0, ph4_p4, ph6_0, ph6_p4, ph8_0, ph8_p4, ph10_0, ph10_p4, ab_2, pc_78 : simd::cache_line_size())
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
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;

        const auto h2_0 = ph2_0[k];
        const auto h4_p4 = ph4_p4[k];
        const auto h6_0 = ph6_0[k];
        const auto h6_p4 = ph6_p4[k];
        const auto h8_0 = ph8_0[k];
        const auto h8_p4 = ph8_p4[k];
        const auto h10_0 = ph10_0[k];
        const auto h10_p4 = ph10_p4[k];

        pc_78[k] = e_0 * fs_1575_8_42 * h2_0 + e_1 * fs_255_4_30 * h4_p4 - e_1 * fs_1125_4_42 * r_2 * h2_0 - e_2 * fs_95_4_42 * h6_0 + e_2 * fs_2185_88_6 * h6_p4 - e_2 * fs_510_11_30 * r_2 * h4_p4 + e_2 * fs_125_1_42 * r_4 * h2_0 + e_3 * fs_63_11_42 * h8_0 - e_3 * fs_777_286_66 * h8_p4 + e_3 * fs_19_2_42 * r_2 * h6_0 - e_3 * fs_437_44_6 * r_2 * h6_p4 + e_3 * fs_1530_143_30 * r_4 * h4_p4 - e_3 * fs_250_11_42 * r_6 * h2_0 + e_4 * fs_2415_4199_42 * h10_0 + e_4 * fs_483_8398_10010 * h10_p4 - e_4 * fs_252_209_42 * r_2 * h8_0 + e_4 * fs_1554_2717_66 * r_2 * h8_p4 - e_4 * fs_19_17_42 * r_4 * h6_0 + e_4 * fs_437_374_6 * r_4 * h6_p4 - e_4 * fs_136_143_30 * r_6 * h4_p4 + e_4 * fs_250_143_42 * r_8 * h2_0 - e_5 * fs_210_4199_42 * r_2 * h10_0 - e_5 * fs_21_4199_10010 * r_2 * h10_p4 + e_5 * fs_12_209_42 * r_4 * h8_0 - e_5 * fs_74_2717_66 * r_4 * h8_p4 + e_5 * fs_2_51_42 * r_6 * h6_0 - e_5 * fs_23_561_6 * r_6 * h6_p4 + e_5 * fs_4_143_30 * r_8 * h4_p4 - e_5 * fs_20_429_42 * r_10 * h2_0;
    }

    // NOTE: the angular components are formed in 93 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph2_p1, ph4_p1, ph6_p1, ph6_p5, ph8_p1, ph8_p5, ph10_p1, ph10_p5, ab_2, pc_79 : simd::cache_line_size())
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
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;

        const auto h2_p1 = ph2_p1[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h6_p5 = ph6_p5[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h8_p5 = ph8_p5[k];
        const auto h10_p1 = ph10_p1[k];
        const auto h10_p5 = ph10_p5[k];

        pc_79[k] = - e_0 * f_1575_4 * h2_p1 - e_1 * fs_765_8_30 * h4_p1 + e_1 * f_1125_2 * r_2 * h2_p1 + e_2 * fs_665_44_7 * h6_p1 - e_2 * fs_665_88_462 * h6_p5 + e_2 * fs_765_11_30 * r_2 * h4_p1 - e_2 * f_250_1 * r_4 * h2_p1 + e_3 * fs_2499_143_3 * h8_p1 + e_3 * fs_21_286_3003 * h8_p5 - e_3 * fs_133_22_7 * r_2 * h6_p1 + e_3 * fs_133_44_462 * r_2 * h6_p5 - e_3 * fs_2295_143_30 * r_4 * h4_p1 + e_3 * f_500_11 * r_6 * h2_p1 + e_4 * fs_483_4199_165 * h10_p1 + e_4 * fs_2415_8398_286 * h10_p5 - e_4 * fs_9996_2717_3 * r_2 * h8_p1 - e_4 * fs_42_2717_3003 * r_2 * h8_p5 + e_4 * fs_133_187_7 * r_4 * h6_p1 - e_4 * fs_133_374_462 * r_4 * h6_p5 + e_4 * fs_204_143_30 * r_6 * h4_p1 - e_4 * f_500_143 * r_8 * h2_p1 - e_5 * fs_42_4199_165 * r_2 * h10_p1 - e_5 * fs_105_4199_286 * r_2 * h10_p5 + e_5 * fs_476_2717_3 * r_4 * h8_p1 + e_5 * fs_2_2717_3003 * r_4 * h8_p5 - e_5 * fs_14_561_7 * r_6 * h6_p1 + e_5 * fs_7_561_462 * r_6 * h6_p5 - e_5 * fs_6_143_30 * r_8 * h4_p1 + e_5 * f_40_429 * r_10 * h2_p1;
    }

    // NOTE: the angular components are formed in 93 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph2_p2, ph4_p2, ph6_p2, ph6_p6, ph8_p2, ph8_p6, ph10_p2, ph10_p6, ab_2, pc_80 : simd::cache_line_size())
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
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;

        const auto h2_p2 = ph2_p2[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h6_p6 = ph6_p6[k];
        const auto h8_p2 = ph8_p2[k];
        const auto h8_p6 = ph8_p6[k];
        const auto h10_p2 = ph10_p2[k];
        const auto h10_p6 = ph10_p6[k];

        pc_80[k] = e_0 * fs_1575_32_2 * h2_p2 + e_1 * fs_255_4_30 * h4_p2 - e_1 * fs_1125_16_2 * r_2 * h2_p2 + e_2 * fs_665_22_35 * h6_p2 + e_2 * fs_285_44_77 * h6_p6 - e_2 * fs_510_11_30 * r_2 * h4_p2 + e_2 * fs_125_4_2 * r_4 * h2_p2 + e_3 * fs_294_143_105 * h8_p2 + e_3 * fs_441_143_143 * h8_p6 - e_3 * fs_133_11_35 * r_2 * h6_p2 - e_3 * fs_57_22_77 * r_2 * h6_p6 + e_3 * fs_1530_143_30 * r_4 * h4_p2 - e_3 * fs_125_22_2 * r_6 * h2_p2 + e_4 * fs_483_8398_110 * h10_p2 + e_4 * fs_483_4199_715 * h10_p6 - e_4 * fs_1176_2717_105 * r_2 * h8_p2 - e_4 * fs_1764_2717_143 * r_2 * h8_p6 + e_4 * fs_266_187_35 * r_4 * h6_p2 + e_4 * fs_57_187_77 * r_4 * h6_p6 - e_4 * fs_136_143_30 * r_6 * h4_p2 + e_4 * fs_125_286_2 * r_8 * h2_p2 - e_5 * fs_21_4199_110 * r_2 * h10_p2 - e_5 * fs_42_4199_715 * r_2 * h10_p6 + e_5 * fs_56_2717_105 * r_4 * h8_p2 + e_5 * fs_84_2717_143 * r_4 * h8_p6 - e_5 * fs_28_561_35 * r_6 * h6_p2 - e_5 * fs_2_187_77 * r_6 * h6_p6 + e_5 * fs_4_143_30 * r_8 * h4_p2 - e_5 * fs_5_429_2 * r_10 * h2_p2;
    }

    // NOTE: the angular components are formed in 93 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph2_m1, ph4_m1, ph6_m6, ph6_m1, ph8_m7, ph8_m6, ph8_m1, ph10_m7, ph10_m6, ph10_m1, ab_2, pc_81, pc_82 : simd::cache_line_size())
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
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;

        const auto h2_m1 = ph2_m1[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h6_m6 = ph6_m6[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h8_m7 = ph8_m7[k];
        const auto h8_m6 = ph8_m6[k];
        const auto h8_m1 = ph8_m1[k];
        const auto h10_m7 = ph10_m7[k];
        const auto h10_m6 = ph10_m6[k];
        const auto h10_m1 = ph10_m1[k];

        pc_81[k] = - e_0 * fs_4725_32_2 * h2_m1 - e_1 * fs_255_2_15 * h4_m1 + e_1 * fs_3375_16_2 * r_2 * h2_m1 - e_2 * fs_1995_44_14 * h6_m1 + e_2 * fs_1020_11_15 * r_2 * h4_m1 - e_2 * fs_375_4_2 * r_4 * h2_m1 + e_3 * fs_147_286_4290 * h8_m7 - e_3 * fs_882_143_6 * h8_m1 + e_3 * fs_399_22_14 * r_2 * h6_m1 - e_3 * fs_3060_143_15 * r_4 * h4_m1 + e_3 * fs_375_22_2 * r_6 * h2_m1 + e_4 * fs_161_4199_12155 * h10_m7 - e_4 * fs_161_8398_330 * h10_m1 - e_4 * fs_294_2717_4290 * r_2 * h8_m7 + e_4 * fs_3528_2717_6 * r_2 * h8_m1 - e_4 * fs_399_187_14 * r_4 * h6_m1 + e_4 * fs_272_143_15 * r_6 * h4_m1 - e_4 * fs_375_286_2 * r_8 * h2_m1 - e_5 * fs_14_4199_12155 * r_2 * h10_m7 + e_5 * fs_7_4199_330 * r_2 * h10_m1 + e_5 * fs_14_2717_4290 * r_4 * h8_m7 - e_5 * fs_168_2717_6 * r_4 * h8_m1 + e_5 * fs_14_187_14 * r_6 * h6_m1 - e_5 * fs_8_143_15 * r_8 * h4_m1 + e_5 * fs_5_143_2 * r_10 * h2_m1;

        pc_82[k] = - e_2 * fs_855_88_154 * h6_m6 - e_3 * fs_441_572_286 * h8_m6 + e_3 * fs_171_44_154 * r_2 * h6_m6 + e_4 * fs_644_4199_1430 * h10_m6 + e_4 * fs_441_2717_286 * r_2 * h8_m6 - e_4 * fs_171_374_154 * r_4 * h6_m6 - e_5 * fs_56_4199_1430 * r_2 * h10_m6 - e_5 * fs_21_2717_286 * r_4 * h8_m6 + e_5 * fs_3_187_154 * r_6 * h6_m6;
    }

    // NOTE: the angular components are formed in 93 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph2_m1, ph4_m1, ph6_m5, ph6_m1, ph8_m5, ph8_m1, ph10_m5, ph10_m1, ab_2, pc_83 : simd::cache_line_size())
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
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;

        const auto h2_m1 = ph2_m1[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h6_m5 = ph6_m5[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h8_m5 = ph8_m5[k];
        const auto h8_m1 = ph8_m1[k];
        const auto h10_m5 = ph10_m5[k];
        const auto h10_m1 = ph10_m1[k];

        pc_83[k] = - e_0 * fs_4725_16_14 * h2_m1 + e_1 * fs_255_8_105 * h4_m1 + e_1 * fs_3375_8_14 * r_2 * h2_m1 + e_2 * fs_285_11_33 * h6_m5 + e_2 * fs_285_11_2 * h6_m1 - e_2 * fs_255_11_105 * r_2 * h4_m1 - e_2 * fs_375_2_14 * r_4 * h2_m1 - e_3 * fs_567_572_858 * h8_m5 - e_3 * fs_2667_572_42 * h8_m1 - e_3 * fs_114_11_33 * r_2 * h6_m5 - e_3 * fs_114_11_2 * r_2 * h6_m1 + e_3 * fs_765_143_105 * r_4 * h4_m1 + e_3 * fs_375_11_14 * r_6 * h2_m1 + e_4 * fs_805_4199_1001 * h10_m5 - e_4 * fs_161_4199_2310 * h10_m1 + e_4 * fs_567_2717_858 * r_2 * h8_m5 + e_4 * fs_2667_2717_42 * r_2 * h8_m1 + e_4 * fs_228_187_33 * r_4 * h6_m5 + e_4 * fs_228_187_2 * r_4 * h6_m1 - e_4 * fs_68_143_105 * r_6 * h4_m1 - e_4 * fs_375_143_14 * r_8 * h2_m1 - e_5 * fs_70_4199_1001 * r_2 * h10_m5 + e_5 * fs_14_4199_2310 * r_2 * h10_m1 - e_5 * fs_27_2717_858 * r_4 * h8_m5 - e_5 * fs_127_2717_42 * r_4 * h8_m1 - e_5 * fs_8_187_33 * r_6 * h6_m5 - e_5 * fs_8_187_2 * r_6 * h6_m1 + e_5 * fs_2_143_105 * r_8 * h4_m1 + e_5 * fs_10_143_14 * r_10 * h2_m1;
    }

    // NOTE: the angular components are formed in 93 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph2_m2, ph4_m4, ph4_m2, ph6_m4, ph6_m2, ph8_m4, ph8_m2, ph10_m4, ph10_m2, ab_2, pc_84 : simd::cache_line_size())
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
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;

        const auto h2_m2 = ph2_m2[k];
        const auto h4_m4 = ph4_m4[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h6_m4 = ph6_m4[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h8_m4 = ph8_m4[k];
        const auto h8_m2 = ph8_m2[k];
        const auto h10_m4 = ph10_m4[k];
        const auto h10_m2 = ph10_m2[k];

        pc_84[k] = - e_0 * fs_4725_16_7 * h2_m2 - e_1 * fs_255_2_15 * h4_m4 + e_1 * fs_255_8_105 * h4_m2 + e_1 * fs_3375_8_7 * r_2 * h2_m2 + e_2 * fs_1995_44_3 * h6_m4 - e_2 * fs_3705_88_10 * h6_m2 + e_2 * fs_1020_11_15 * r_2 * h4_m4 - e_2 * fs_255_11_105 * r_2 * h4_m2 - e_2 * fs_375_2_7 * r_4 * h2_m2 - e_3 * fs_903_286_33 * h8_m4 + e_3 * fs_2331_572_30 * h8_m2 - e_3 * fs_399_22_3 * r_2 * h6_m4 + e_3 * fs_741_44_10 * r_2 * h6_m2 - e_3 * fs_3060_143_15 * r_4 * h4_m4 + e_3 * fs_765_143_105 * r_4 * h4_m2 + e_3 * fs_375_11_7 * r_6 * h2_m2 + e_4 * fs_322_4199_5005 * h10_m4 + e_4 * fs_644_4199_385 * h10_m2 + e_4 * fs_1806_2717_33 * r_2 * h8_m4 - e_4 * fs_2331_2717_30 * r_2 * h8_m2 + e_4 * fs_399_187_3 * r_4 * h6_m4 - e_4 * fs_741_374_10 * r_4 * h6_m2 + e_4 * fs_272_143_15 * r_6 * h4_m4 - e_4 * fs_68_143_105 * r_6 * h4_m2 - e_4 * fs_375_143_7 * r_8 * h2_m2 - e_5 * fs_28_4199_5005 * r_2 * h10_m4 - e_5 * fs_56_4199_385 * r_2 * h10_m2 - e_5 * fs_86_2717_33 * r_4 * h8_m4 + e_5 * fs_111_2717_30 * r_4 * h8_m2 - e_5 * fs_14_187_3 * r_6 * h6_m4 + e_5 * fs_13_187_10 * r_6 * h6_m2 - e_5 * fs_8_143_15 * r_8 * h4_m4 + e_5 * fs_2_143_105 * r_8 * h4_m2 + e_5 * fs_10_143_7 * r_10 * h2_m2;
    }

    // NOTE: the angular components are formed in 93 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph4_p3, ph6_p3, ph8_p3, ph10_p3, ab_2, pc_85 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

        const auto h4_p3 = ph4_p3[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h8_p3 = ph8_p3[k];
        const auto h10_p3 = ph10_p3[k];

        pc_85[k] = e_1 * fs_1275_4_3 * h4_p3 - e_2 * f_2565_22 * h6_p3 - e_2 * fs_2550_11_3 * r_2 * h4_p3 + e_3 * fs_315_286_22 * h8_p3 + e_3 * f_513_11 * r_2 * h6_p3 + e_3 * fs_7650_143_3 * r_4 * h4_p3 + e_4 * fs_805_4199_1001 * h10_p3 - e_4 * fs_630_2717_22 * r_2 * h8_p3 - e_4 * f_1026_187 * r_4 * h6_p3 - e_4 * fs_680_143_3 * r_6 * h4_p3 - e_5 * fs_70_4199_1001 * r_2 * h10_p3 + e_5 * fs_30_2717_22 * r_4 * h8_p3 + e_5 * f_36_187 * r_6 * h6_p3 + e_5 * fs_20_143_3 * r_8 * h4_p3;
    }

    // NOTE: the angular components are formed in 93 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph2_p2, ph4_p2, ph4_p4, ph6_p2, ph6_p4, ph8_p2, ph8_p4, ph10_p2, ph10_p4, ab_2, pc_86 : simd::cache_line_size())
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
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;

        const auto h2_p2 = ph2_p2[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h4_p4 = ph4_p4[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h6_p4 = ph6_p4[k];
        const auto h8_p2 = ph8_p2[k];
        const auto h8_p4 = ph8_p4[k];
        const auto h10_p2 = ph10_p2[k];
        const auto h10_p4 = ph10_p4[k];

        pc_86[k] = e_0 * fs_4725_16_7 * h2_p2 - e_1 * fs_255_8_105 * h4_p2 - e_1 * fs_255_2_15 * h4_p4 - e_1 * fs_3375_8_7 * r_2 * h2_p2 + e_2 * fs_3705_88_10 * h6_p2 + e_2 * fs_1995_44_3 * h6_p4 + e_2 * fs_255_11_105 * r_2 * h4_p2 + e_2 * fs_1020_11_15 * r_2 * h4_p4 + e_2 * fs_375_2_7 * r_4 * h2_p2 - e_3 * fs_2331_572_30 * h8_p2 - e_3 * fs_903_286_33 * h8_p4 - e_3 * fs_741_44_10 * r_2 * h6_p2 - e_3 * fs_399_22_3 * r_2 * h6_p4 - e_3 * fs_765_143_105 * r_4 * h4_p2 - e_3 * fs_3060_143_15 * r_4 * h4_p4 - e_3 * fs_375_11_7 * r_6 * h2_p2 - e_4 * fs_644_4199_385 * h10_p2 + e_4 * fs_322_4199_5005 * h10_p4 + e_4 * fs_2331_2717_30 * r_2 * h8_p2 + e_4 * fs_1806_2717_33 * r_2 * h8_p4 + e_4 * fs_741_374_10 * r_4 * h6_p2 + e_4 * fs_399_187_3 * r_4 * h6_p4 + e_4 * fs_68_143_105 * r_6 * h4_p2 + e_4 * fs_272_143_15 * r_6 * h4_p4 + e_4 * fs_375_143_7 * r_8 * h2_p2 + e_5 * fs_56_4199_385 * r_2 * h10_p2 - e_5 * fs_28_4199_5005 * r_2 * h10_p4 - e_5 * fs_111_2717_30 * r_4 * h8_p2 - e_5 * fs_86_2717_33 * r_4 * h8_p4 - e_5 * fs_13_187_10 * r_6 * h6_p2 - e_5 * fs_14_187_3 * r_6 * h6_p4 - e_5 * fs_2_143_105 * r_8 * h4_p2 - e_5 * fs_8_143_15 * r_8 * h4_p4 - e_5 * fs_10_143_7 * r_10 * h2_p2;
    }

    // NOTE: the angular components are formed in 93 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph2_p1, ph4_p1, ph6_p1, ph6_p5, ph8_p1, ph8_p5, ph10_p1, ph10_p5, ab_2, pc_87 : simd::cache_line_size())
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
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;

        const auto h2_p1 = ph2_p1[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h6_p5 = ph6_p5[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h8_p5 = ph8_p5[k];
        const auto h10_p1 = ph10_p1[k];
        const auto h10_p5 = ph10_p5[k];

        pc_87[k] = e_0 * fs_4725_16_14 * h2_p1 - e_1 * fs_255_8_105 * h4_p1 - e_1 * fs_3375_8_14 * r_2 * h2_p1 - e_2 * fs_285_11_2 * h6_p1 + e_2 * fs_285_11_33 * h6_p5 + e_2 * fs_255_11_105 * r_2 * h4_p1 + e_2 * fs_375_2_14 * r_4 * h2_p1 + e_3 * fs_2667_572_42 * h8_p1 - e_3 * fs_567_572_858 * h8_p5 + e_3 * fs_114_11_2 * r_2 * h6_p1 - e_3 * fs_114_11_33 * r_2 * h6_p5 - e_3 * fs_765_143_105 * r_4 * h4_p1 - e_3 * fs_375_11_14 * r_6 * h2_p1 + e_4 * fs_161_4199_2310 * h10_p1 + e_4 * fs_805_4199_1001 * h10_p5 - e_4 * fs_2667_2717_42 * r_2 * h8_p1 + e_4 * fs_567_2717_858 * r_2 * h8_p5 - e_4 * fs_228_187_2 * r_4 * h6_p1 + e_4 * fs_228_187_33 * r_4 * h6_p5 + e_4 * fs_68_143_105 * r_6 * h4_p1 + e_4 * fs_375_143_14 * r_8 * h2_p1 - e_5 * fs_14_4199_2310 * r_2 * h10_p1 - e_5 * fs_70_4199_1001 * r_2 * h10_p5 + e_5 * fs_127_2717_42 * r_4 * h8_p1 - e_5 * fs_27_2717_858 * r_4 * h8_p5 + e_5 * fs_8_187_2 * r_6 * h6_p1 - e_5 * fs_8_187_33 * r_6 * h6_p5 - e_5 * fs_2_143_105 * r_8 * h4_p1 - e_5 * fs_10_143_14 * r_10 * h2_p1;
    }

    // NOTE: the angular components are formed in 93 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph2_0, ph4_0, ph6_0, ph6_p6, ph8_0, ph8_p6, ph10_0, ph10_p6, ab_2, pc_88 : simd::cache_line_size())
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
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;

        const auto h2_0 = ph2_0[k];
        const auto h4_0 = ph4_0[k];
        const auto h6_0 = ph6_0[k];
        const auto h6_p6 = ph6_p6[k];
        const auto h8_0 = ph8_0[k];
        const auto h8_p6 = ph8_p6[k];
        const auto h10_0 = ph10_0[k];
        const auto h10_p6 = ph10_p6[k];

        pc_88[k] = e_0 * fs_4725_8_3 * h2_0 + e_1 * fs_1275_4_3 * h4_0 - e_1 * fs_3375_4_3 * r_2 * h2_0 - e_2 * fs_1995_22_3 * h6_0 - e_2 * fs_855_88_154 * h6_p6 - e_2 * fs_2550_11_3 * r_2 * h4_0 + e_2 * fs_375_1_3 * r_4 * h2_0 - e_3 * fs_3087_143_3 * h8_0 - e_3 * fs_441_572_286 * h8_p6 + e_3 * fs_399_11_3 * r_2 * h6_0 + e_3 * fs_171_44_154 * r_2 * h6_p6 + e_3 * fs_7650_143_3 * r_4 * h4_0 - e_3 * fs_750_11_3 * r_6 * h2_0 - e_4 * fs_3220_4199_3 * h10_0 + e_4 * fs_644_4199_1430 * h10_p6 + e_4 * fs_12348_2717_3 * r_2 * h8_0 + e_4 * fs_441_2717_286 * r_2 * h8_p6 - e_4 * fs_798_187_3 * r_4 * h6_0 - e_4 * fs_171_374_154 * r_4 * h6_p6 - e_4 * fs_680_143_3 * r_6 * h4_0 + e_4 * fs_750_143_3 * r_8 * h2_0 + e_5 * fs_280_4199_3 * r_2 * h10_0 - e_5 * fs_56_4199_1430 * r_2 * h10_p6 - e_5 * fs_588_2717_3 * r_4 * h8_0 - e_5 * fs_21_2717_286 * r_4 * h8_p6 + e_5 * fs_28_187_3 * r_6 * h6_0 + e_5 * fs_3_187_154 * r_6 * h6_p6 + e_5 * fs_20_143_3 * r_8 * h4_0 - e_5 * fs_20_143_3 * r_10 * h2_0;
    }

    // NOTE: the angular components are formed in 93 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph2_p1, ph4_p1, ph6_p1, ph8_m8, ph8_p1, ph8_p7, ph10_m8, ph10_p1, ph10_p7, ab_2, pc_89, pc_90 : simd::cache_line_size())
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
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;

        const auto h2_p1 = ph2_p1[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h8_m8 = ph8_m8[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h8_p7 = ph8_p7[k];
        const auto h10_m8 = ph10_m8[k];
        const auto h10_p1 = ph10_p1[k];
        const auto h10_p7 = ph10_p7[k];

        pc_89[k] = - e_0 * fs_4725_32_2 * h2_p1 - e_1 * fs_255_2_15 * h4_p1 + e_1 * fs_3375_16_2 * r_2 * h2_p1 - e_2 * fs_1995_44_14 * h6_p1 + e_2 * fs_1020_11_15 * r_2 * h4_p1 - e_2 * fs_375_4_2 * r_4 * h2_p1 - e_3 * fs_882_143_6 * h8_p1 + e_3 * fs_147_286_4290 * h8_p7 + e_3 * fs_399_22_14 * r_2 * h6_p1 - e_3 * fs_3060_143_15 * r_4 * h4_p1 + e_3 * fs_375_22_2 * r_6 * h2_p1 - e_4 * fs_161_8398_330 * h10_p1 + e_4 * fs_161_4199_12155 * h10_p7 + e_4 * fs_3528_2717_6 * r_2 * h8_p1 - e_4 * fs_294_2717_4290 * r_2 * h8_p7 - e_4 * fs_399_187_14 * r_4 * h6_p1 + e_4 * fs_272_143_15 * r_6 * h4_p1 - e_4 * fs_375_286_2 * r_8 * h2_p1 + e_5 * fs_7_4199_330 * r_2 * h10_p1 - e_5 * fs_14_4199_12155 * r_2 * h10_p7 - e_5 * fs_168_2717_6 * r_4 * h8_p1 + e_5 * fs_14_2717_4290 * r_4 * h8_p7 + e_5 * fs_14_187_14 * r_6 * h6_p1 - e_5 * fs_8_143_15 * r_8 * h4_p1 + e_5 * fs_5_143_2 * r_10 * h2_p1;

        pc_90[k] = e_3 * fs_294_143_143 * h8_m8 + e_4 * fs_483_4199_2431 * h10_m8 - e_4 * fs_1176_2717_143 * r_2 * h8_m8 - e_5 * fs_42_4199_2431 * r_2 * h10_m8 + e_5 * fs_56_2717_143 * r_4 * h8_m8;
    }

    // NOTE: the angular components are formed in 93 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph2_m1, ph4_m1, ph6_m1, ph8_m7, ph8_m1, ph10_m7, ph10_m1, ab_2, pc_91 : simd::cache_line_size())
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
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;

        const auto h2_m1 = ph2_m1[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h8_m7 = ph8_m7[k];
        const auto h8_m1 = ph8_m1[k];
        const auto h10_m7 = ph10_m7[k];
        const auto h10_m1 = ph10_m1[k];

        pc_91[k] = - e_0 * fs_1575_8_30 * h2_m1 - e_1 * f_765_8 * h4_m1 + e_1 * fs_1125_4_30 * r_2 * h2_m1 + e_2 * fs_475_44_210 * h6_m1 + e_2 * f_765_11 * r_2 * h4_m1 - e_2 * fs_125_1_30 * r_4 * h2_m1 - e_3 * fs_1029_572_286 * h8_m7 + e_3 * fs_3675_572_10 * h8_m1 - e_3 * fs_95_22_210 * r_2 * h6_m1 - e_3 * f_2295_143 * r_4 * h4_m1 + e_3 * fs_250_11_30 * r_6 * h2_m1 + e_4 * fs_322_4199_7293 * h10_m7 + e_4 * fs_483_4199_22 * h10_m1 + e_4 * fs_1029_2717_286 * r_2 * h8_m7 - e_4 * fs_3675_2717_10 * r_2 * h8_m1 + e_4 * fs_95_187_210 * r_4 * h6_m1 + e_4 * f_204_143 * r_6 * h4_m1 - e_4 * fs_250_143_30 * r_8 * h2_m1 - e_5 * fs_28_4199_7293 * r_2 * h10_m7 - e_5 * fs_42_4199_22 * r_2 * h10_m1 - e_5 * fs_49_2717_286 * r_4 * h8_m7 + e_5 * fs_175_2717_10 * r_4 * h8_m1 - e_5 * fs_10_561_210 * r_6 * h6_m1 - e_5 * f_6_143 * r_8 * h4_m1 + e_5 * fs_20_429_30 * r_10 * h2_m1;
    }

    // NOTE: the angular components are formed in 93 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph2_m2, ph4_m2, ph6_m6, ph6_m2, ph8_m6, ph8_m2, ph10_m6, ph10_m2, ab_2, pc_92 : simd::cache_line_size())
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
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;

        const auto h2_m2 = ph2_m2[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h6_m6 = ph6_m6[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h8_m6 = ph8_m6[k];
        const auto h8_m2 = ph8_m2[k];
        const auto h10_m6 = ph10_m6[k];
        const auto h10_m2 = ph10_m2[k];

        pc_92[k] = - e_0 * fs_1575_16_105 * h2_m2 + e_1 * fs_765_4_7 * h4_m2 + e_1 * fs_1125_8_105 * r_2 * h2_m2 + e_2 * fs_855_88_330 * h6_m6 - e_2 * fs_2185_88_6 * h6_m2 - e_2 * fs_1530_11_7 * r_2 * h4_m2 - e_2 * fs_125_2_105 * r_4 * h2_m2 - e_3 * fs_21_143_30030 * h8_m6 - e_3 * fs_2919_143_2 * h8_m2 - e_3 * fs_171_44_330 * r_2 * h6_m6 + e_3 * fs_437_44_6 * r_2 * h6_m2 + e_3 * fs_4590_143_7 * r_4 * h4_m2 + e_3 * fs_125_11_105 * r_6 * h2_m2 + e_4 * fs_322_4199_6006 * h10_m6 - e_4 * fs_322_4199_231 * h10_m2 + e_4 * fs_84_2717_30030 * r_2 * h8_m6 + e_4 * fs_11676_2717_2 * r_2 * h8_m2 + e_4 * fs_171_374_330 * r_4 * h6_m6 - e_4 * fs_437_374_6 * r_4 * h6_m2 - e_4 * fs_408_143_7 * r_6 * h4_m2 - e_4 * fs_125_143_105 * r_8 * h2_m2 - e_5 * fs_28_4199_6006 * r_2 * h10_m6 + e_5 * fs_28_4199_231 * r_2 * h10_m2 - e_5 * fs_4_2717_30030 * r_4 * h8_m6 - e_5 * fs_556_2717_2 * r_4 * h8_m2 - e_5 * fs_3_187_330 * r_6 * h6_m6 + e_5 * fs_23_561_6 * r_6 * h6_m2 + e_5 * fs_12_143_7 * r_8 * h4_m2 + e_5 * fs_10_429_105 * r_10 * h2_m2;
    }

    // NOTE: the angular components are formed in 93 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph4_m3, ph4_p4, ph6_m5, ph6_m3, ph6_p4, ph8_m5, ph8_m3, ph8_p4, ph10_m5, ph10_m3, ph10_p4, ab_2, pc_93, pc_94 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

        const auto h4_m3 = ph4_m3[k];
        const auto h4_p4 = ph4_p4[k];
        const auto h6_m5 = ph6_m5[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h6_p4 = ph6_p4[k];
        const auto h8_m5 = ph8_m5[k];
        const auto h8_m3 = ph8_m3[k];
        const auto h8_p4 = ph8_p4[k];
        const auto h10_m5 = ph10_m5[k];
        const auto h10_m3 = ph10_m3[k];
        const auto h10_p4 = ph10_p4[k];

        pc_93[k] = - e_1 * f_765_8 * h4_m3 - e_2 * fs_285_44_55 * h6_m5 - e_2 * fs_1995_44_3 * h6_m3 + e_2 * f_765_11 * r_2 * h4_m3 - e_3 * fs_21_572_1430 * h8_m5 + e_3 * fs_2121_572_66 * h8_m3 + e_3 * fs_57_22_55 * r_2 * h6_m5 + e_3 * fs_399_22_3 * r_2 * h6_m3 - e_3 * f_2295_143 * r_4 * h4_m3 + e_4 * fs_161_4199_15015 * h10_m5 + e_4 * fs_161_4199_3003 * h10_m3 + e_4 * fs_21_2717_1430 * r_2 * h8_m5 - e_4 * fs_2121_2717_66 * r_2 * h8_m3 - e_4 * fs_57_187_55 * r_4 * h6_m5 - e_4 * fs_399_187_3 * r_4 * h6_m3 + e_4 * f_204_143 * r_6 * h4_m3 - e_5 * fs_14_4199_15015 * r_2 * h10_m5 - e_5 * fs_14_4199_3003 * r_2 * h10_m3 - e_5 * fs_1_2717_1430 * r_4 * h8_m5 + e_5 * fs_101_2717_66 * r_4 * h8_m3 + e_5 * fs_2_187_55 * r_6 * h6_m5 + e_5 * fs_14_187_3 * r_6 * h6_m3 - e_5 * f_6_143 * r_8 * h4_m3;

        pc_94[k] = e_1 * fs_765_2_5 * h4_p4 - e_2 * f_2280_11 * h6_p4 - e_2 * fs_3060_11_5 * r_2 * h4_p4 + e_3 * fs_1218_143_11 * h8_p4 + e_3 * f_912_11 * r_2 * h6_p4 + e_3 * fs_9180_143_5 * r_4 * h4_p4 + e_4 * fs_161_4199_15015 * h10_p4 - e_4 * fs_4872_2717_11 * r_2 * h8_p4 - e_4 * f_1824_187 * r_4 * h6_p4 - e_4 * fs_816_143_5 * r_6 * h4_p4 - e_5 * fs_14_4199_15015 * r_2 * h10_p4 + e_5 * fs_232_2717_11 * r_4 * h8_p4 + e_5 * f_64_187 * r_6 * h6_p4 + e_5 * fs_24_143_5 * r_8 * h4_p4;
    }

    // NOTE: the angular components are formed in 93 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph4_p3, ph6_p3, ph6_p5, ph8_p3, ph8_p5, ph10_p3, ph10_p5, ab_2, pc_95 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

        const auto h4_p3 = ph4_p3[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h6_p5 = ph6_p5[k];
        const auto h8_p3 = ph8_p3[k];
        const auto h8_p5 = ph8_p5[k];
        const auto h10_p3 = ph10_p3[k];
        const auto h10_p5 = ph10_p5[k];

        pc_95[k] = e_1 * f_765_8 * h4_p3 + e_2 * fs_1995_44_3 * h6_p3 - e_2 * fs_285_44_55 * h6_p5 - e_2 * f_765_11 * r_2 * h4_p3 - e_3 * fs_2121_572_66 * h8_p3 - e_3 * fs_21_572_1430 * h8_p5 - e_3 * fs_399_22_3 * r_2 * h6_p3 + e_3 * fs_57_22_55 * r_2 * h6_p5 + e_3 * f_2295_143 * r_4 * h4_p3 - e_4 * fs_161_4199_3003 * h10_p3 + e_4 * fs_161_4199_15015 * h10_p5 + e_4 * fs_2121_2717_66 * r_2 * h8_p3 + e_4 * fs_21_2717_1430 * r_2 * h8_p5 + e_4 * fs_399_187_3 * r_4 * h6_p3 - e_4 * fs_57_187_55 * r_4 * h6_p5 - e_4 * f_204_143 * r_6 * h4_p3 + e_5 * fs_14_4199_3003 * r_2 * h10_p3 - e_5 * fs_14_4199_15015 * r_2 * h10_p5 - e_5 * fs_101_2717_66 * r_4 * h8_p3 - e_5 * fs_1_2717_1430 * r_4 * h8_p5 - e_5 * fs_14_187_3 * r_6 * h6_p3 + e_5 * fs_2_187_55 * r_6 * h6_p5 + e_5 * f_6_143 * r_8 * h4_p3;
    }

    // NOTE: the angular components are formed in 93 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph2_p2, ph4_p2, ph6_p2, ph6_p6, ph8_p2, ph8_p6, ph10_p2, ph10_p6, ab_2, pc_96 : simd::cache_line_size())
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
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;

        const auto h2_p2 = ph2_p2[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h6_p6 = ph6_p6[k];
        const auto h8_p2 = ph8_p2[k];
        const auto h8_p6 = ph8_p6[k];
        const auto h10_p2 = ph10_p2[k];
        const auto h10_p6 = ph10_p6[k];

        pc_96[k] = e_0 * fs_1575_16_105 * h2_p2 - e_1 * fs_765_4_7 * h4_p2 - e_1 * fs_1125_8_105 * r_2 * h2_p2 + e_2 * fs_2185_88_6 * h6_p2 + e_2 * fs_855_88_330 * h6_p6 + e_2 * fs_1530_11_7 * r_2 * h4_p2 + e_2 * fs_125_2_105 * r_4 * h2_p2 + e_3 * fs_2919_143_2 * h8_p2 - e_3 * fs_21_143_30030 * h8_p6 - e_3 * fs_437_44_6 * r_2 * h6_p2 - e_3 * fs_171_44_330 * r_2 * h6_p6 - e_3 * fs_4590_143_7 * r_4 * h4_p2 - e_3 * fs_125_11_105 * r_6 * h2_p2 + e_4 * fs_322_4199_231 * h10_p2 + e_4 * fs_322_4199_6006 * h10_p6 - e_4 * fs_11676_2717_2 * r_2 * h8_p2 + e_4 * fs_84_2717_30030 * r_2 * h8_p6 + e_4 * fs_437_374_6 * r_4 * h6_p2 + e_4 * fs_171_374_330 * r_4 * h6_p6 + e_4 * fs_408_143_7 * r_6 * h4_p2 + e_4 * fs_125_143_105 * r_8 * h2_p2 - e_5 * fs_28_4199_231 * r_2 * h10_p2 - e_5 * fs_28_4199_6006 * r_2 * h10_p6 + e_5 * fs_556_2717_2 * r_4 * h8_p2 - e_5 * fs_4_2717_30030 * r_4 * h8_p6 - e_5 * fs_23_561_6 * r_6 * h6_p2 - e_5 * fs_3_187_330 * r_6 * h6_p6 - e_5 * fs_12_143_7 * r_8 * h4_p2 - e_5 * fs_10_429_105 * r_10 * h2_p2;
    }

    // NOTE: the angular components are formed in 93 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph2_p1, ph4_p1, ph6_p1, ph8_p1, ph8_p7, ph10_p1, ph10_p7, ab_2, pc_97 : simd::cache_line_size())
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
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;

        const auto h2_p1 = ph2_p1[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h8_p7 = ph8_p7[k];
        const auto h10_p1 = ph10_p1[k];
        const auto h10_p7 = ph10_p7[k];

        pc_97[k] = e_0 * fs_1575_8_30 * h2_p1 + e_1 * f_765_8 * h4_p1 - e_1 * fs_1125_4_30 * r_2 * h2_p1 - e_2 * fs_475_44_210 * h6_p1 - e_2 * f_765_11 * r_2 * h4_p1 + e_2 * fs_125_1_30 * r_4 * h2_p1 - e_3 * fs_3675_572_10 * h8_p1 - e_3 * fs_1029_572_286 * h8_p7 + e_3 * fs_95_22_210 * r_2 * h6_p1 + e_3 * f_2295_143 * r_4 * h4_p1 - e_3 * fs_250_11_30 * r_6 * h2_p1 - e_4 * fs_483_4199_22 * h10_p1 + e_4 * fs_322_4199_7293 * h10_p7 + e_4 * fs_3675_2717_10 * r_2 * h8_p1 + e_4 * fs_1029_2717_286 * r_2 * h8_p7 - e_4 * fs_95_187_210 * r_4 * h6_p1 - e_4 * f_204_143 * r_6 * h4_p1 + e_4 * fs_250_143_30 * r_8 * h2_p1 + e_5 * fs_42_4199_22 * r_2 * h10_p1 - e_5 * fs_28_4199_7293 * r_2 * h10_p7 - e_5 * fs_175_2717_10 * r_4 * h8_p1 - e_5 * fs_49_2717_286 * r_4 * h8_p7 + e_5 * fs_10_561_210 * r_6 * h6_p1 + e_5 * f_6_143 * r_8 * h4_p1 - e_5 * fs_20_429_30 * r_10 * h2_p1;
    }

    // NOTE: the angular components are formed in 93 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph2_0, ph4_0, ph6_0, ph8_0, ph8_p8, ph10_0, ph10_p8, ab_2, pc_98 : simd::cache_line_size())
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
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;

        const auto h2_0 = ph2_0[k];
        const auto h4_0 = ph4_0[k];
        const auto h6_0 = ph6_0[k];
        const auto h8_0 = ph8_0[k];
        const auto h8_p8 = ph8_p8[k];
        const auto h10_0 = ph10_0[k];
        const auto h10_p8 = ph10_p8[k];

        pc_98[k] = e_0 * fs_4725_16_5 * h2_0 + e_1 * fs_765_2_5 * h4_0 - e_1 * fs_3375_8_5 * r_2 * h2_0 + e_2 * fs_1995_22_5 * h6_0 - e_2 * fs_3060_11_5 * r_2 * h4_0 + e_2 * fs_375_2_5 * r_4 * h2_0 + e_3 * fs_882_143_5 * h8_0 + e_3 * fs_294_143_143 * h8_p8 - e_3 * fs_399_11_5 * r_2 * h6_0 + e_3 * fs_9180_143_5 * r_4 * h4_0 - e_3 * fs_375_11_5 * r_6 * h2_0 + e_4 * fs_483_4199_5 * h10_0 + e_4 * fs_483_4199_2431 * h10_p8 - e_4 * fs_3528_2717_5 * r_2 * h8_0 - e_4 * fs_1176_2717_143 * r_2 * h8_p8 + e_4 * fs_798_187_5 * r_4 * h6_0 - e_4 * fs_816_143_5 * r_6 * h4_0 + e_4 * fs_375_143_5 * r_8 * h2_0 - e_5 * fs_42_4199_5 * r_2 * h10_0 - e_5 * fs_42_4199_2431 * r_2 * h10_p8 + e_5 * fs_168_2717_5 * r_4 * h8_0 + e_5 * fs_56_2717_143 * r_4 * h8_p8 - e_5 * fs_28_187_5 * r_6 * h6_0 + e_5 * fs_24_143_5 * r_8 * h4_0 - e_5 * fs_10_143_5 * r_10 * h2_0;
    }

    // NOTE: the angular components are formed in 93 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph2_m1, ph4_m1, ph6_m1, ph8_m1, ph10_m9, ph10_m1, ab_2, pc_99 : simd::cache_line_size())
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
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;

        const auto h2_m1 = ph2_m1[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h8_m1 = ph8_m1[k];
        const auto h10_m9 = ph10_m9[k];
        const auto h10_m1 = ph10_m1[k];

        pc_99[k] = - e_0 * fs_1575_32_330 * h2_m1 - e_1 * fs_765_4_11 * h4_m1 + e_1 * fs_1125_16_330 * r_2 * h2_m1 - e_2 * fs_95_44_2310 * h6_m1 + e_2 * fs_1530_11_11 * r_2 * h4_m1 - e_2 * fs_125_4_330 * r_4 * h2_m1 - e_3 * fs_147_286_110 * h8_m1 + e_3 * fs_19_22_2310 * r_2 * h6_m1 - e_3 * fs_4590_143_11 * r_4 * h4_m1 + e_3 * fs_125_22_330 * r_6 * h2_m1 + e_4 * fs_483_4199_4199 * h10_m9 - e_4 * fs_483_8398_2 * h10_m1 + e_4 * fs_294_2717_110 * r_2 * h8_m1 - e_4 * fs_19_187_2310 * r_4 * h6_m1 + e_4 * fs_408_143_11 * r_6 * h4_m1 - e_4 * fs_125_286_330 * r_8 * h2_m1 - e_5 * fs_42_4199_4199 * r_2 * h10_m9 + e_5 * fs_21_4199_2 * r_2 * h10_m1 - e_5 * fs_14_2717_110 * r_4 * h8_m1 + e_5 * fs_2_561_2310 * r_6 * h6_m1 - e_5 * fs_12_143_11 * r_8 * h4_m1 + e_5 * fs_5_429_330 * r_10 * h2_m1;
    }

    // NOTE: the angular components are formed in 93 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph2_m2, ph4_m2, ph6_m2, ph8_m8, ph8_m2, ph10_m8, ph10_m2, ab_2, pc_100 : simd::cache_line_size())
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
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;

        const auto h2_m2 = ph2_m2[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h8_m8 = ph8_m8[k];
        const auto h8_m2 = ph8_m2[k];
        const auto h10_m8 = ph10_m8[k];
        const auto h10_m2 = ph10_m2[k];

        pc_100[k] = - e_0 * fs_1575_16_165 * h2_m2 + e_1 * fs_765_8_11 * h4_m2 + e_1 * fs_1125_8_165 * r_2 * h2_m2 + e_2 * fs_665_88_462 * h6_m2 - e_2 * fs_765_11_11 * r_2 * h4_m2 - e_2 * fs_125_2_165 * r_4 * h2_m2 - e_3 * fs_147_13_13 * h8_m8 + e_3 * fs_609_572_154 * h8_m2 - e_3 * fs_133_44_462 * r_2 * h6_m2 + e_3 * fs_2295_143_11 * r_4 * h4_m2 + e_3 * fs_125_11_165 * r_6 * h2_m2 + e_4 * fs_1932_4199_221 * h10_m8 + e_4 * fs_644_4199_3 * h10_m2 + e_4 * fs_588_247_13 * r_2 * h8_m8 - e_4 * fs_609_2717_154 * r_2 * h8_m2 + e_4 * fs_133_374_462 * r_4 * h6_m2 - e_4 * fs_204_143_11 * r_6 * h4_m2 - e_4 * fs_125_143_165 * r_8 * h2_m2 - e_5 * fs_168_4199_221 * r_2 * h10_m8 - e_5 * fs_56_4199_3 * r_2 * h10_m2 - e_5 * fs_28_247_13 * r_4 * h8_m8 + e_5 * fs_29_2717_154 * r_4 * h8_m2 - e_5 * fs_7_561_462 * r_6 * h6_m2 + e_5 * fs_6_143_11 * r_8 * h4_m2 + e_5 * fs_10_429_165 * r_10 * h2_m2;
    }

    // NOTE: the angular components are formed in 93 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph4_m3, ph6_m3, ph8_m7, ph8_m3, ph10_m7, ph10_m3, ab_2, pc_101 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

        const auto h4_m3 = ph4_m3[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h8_m7 = ph8_m7[k];
        const auto h8_m3 = ph8_m3[k];
        const auto h10_m7 = ph10_m7[k];
        const auto h10_m3 = ph10_m3[k];

        pc_101[k] = e_1 * fs_765_8_11 * h4_m3 - e_2 * fs_285_11_33 * h6_m3 - e_2 * fs_765_11_11 * r_2 * h4_m3 - e_3 * fs_21_52_182 * h8_m7 - e_3 * fs_483_52_6 * h8_m3 + e_3 * fs_114_11_33 * r_2 * h6_m3 + e_3 * fs_2295_143_11 * r_4 * h4_m3 + e_4 * fs_322_4199_4641 * h10_m7 - e_4 * fs_161_4199_273 * h10_m3 + e_4 * fs_21_247_182 * r_2 * h8_m7 + e_4 * fs_483_247_6 * r_2 * h8_m3 - e_4 * fs_228_187_33 * r_4 * h6_m3 - e_4 * fs_204_143_11 * r_6 * h4_m3 - e_5 * fs_28_4199_4641 * r_2 * h10_m7 + e_5 * fs_14_4199_273 * r_2 * h10_m3 - e_5 * fs_1_247_182 * r_4 * h8_m7 - e_5 * fs_23_247_6 * r_4 * h8_m3 + e_5 * fs_8_187_33 * r_6 * h6_m3 + e_5 * fs_6_143_11 * r_8 * h4_m3;
    }

    // NOTE: the angular components are formed in 93 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph4_m4, ph6_m6, ph6_m4, ph6_p5, ph8_m6, ph8_m4, ph8_p5, ph10_m6, ph10_m4, ph10_p5, ab_2, pc_102, pc_103 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

        const auto h4_m4 = ph4_m4[k];
        const auto h6_m6 = ph6_m6[k];
        const auto h6_m4 = ph6_m4[k];
        const auto h6_p5 = ph6_p5[k];
        const auto h8_m6 = ph8_m6[k];
        const auto h8_m4 = ph8_m4[k];
        const auto h8_p5 = ph8_p5[k];
        const auto h10_m6 = ph10_m6[k];
        const auto h10_m4 = ph10_m4[k];
        const auto h10_p5 = ph10_p5[k];

        pc_102[k] = - e_1 * fs_765_4_11 * h4_m4 - e_2 * fs_285_8_30 * h6_m6 + e_2 * fs_285_44_55 * h6_m4 + e_2 * fs_1530_11_11 * r_2 * h4_m4 + e_3 * fs_21_52_2730 * h8_m6 + e_3 * fs_357_26_5 * h8_m4 + e_3 * fs_57_4_30 * r_2 * h6_m6 - e_3 * fs_57_22_55 * r_2 * h6_m4 - e_3 * fs_4590_143_11 * r_4 * h4_m4 + e_4 * fs_644_4199_546 * h10_m6 + e_4 * fs_322_4199_273 * h10_m4 - e_4 * fs_21_247_2730 * r_2 * h8_m6 - e_4 * fs_714_247_5 * r_2 * h8_m4 - e_4 * fs_57_34_30 * r_4 * h6_m6 + e_4 * fs_57_187_55 * r_4 * h6_m4 + e_4 * fs_408_143_11 * r_6 * h4_m4 - e_5 * fs_56_4199_546 * r_2 * h10_m6 - e_5 * fs_28_4199_273 * r_2 * h10_m4 + e_5 * fs_1_247_2730 * r_4 * h8_m6 + e_5 * fs_34_247_5 * r_4 * h8_m4 + e_5 * fs_1_17_30 * r_6 * h6_m6 - e_5 * fs_2_187_55 * r_6 * h6_m4 - e_5 * fs_12_143_11 * r_8 * h4_m4;

        pc_103[k] = - e_2 * f_285_2 * h6_p5 + e_3 * fs_231_26_26 * h8_p5 + e_3 * f_57_1 * r_2 * h6_p5 + e_4 * fs_805_4199_273 * h10_p5 - e_4 * fs_462_247_26 * r_2 * h8_p5 - e_4 * f_114_17 * r_4 * h6_p5 - e_5 * fs_70_4199_273 * r_2 * h10_p5 + e_5 * fs_22_247_26 * r_4 * h8_p5 + e_5 * f_4_17 * r_6 * h6_p5;
    }

    // NOTE: the angular components are formed in 93 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph4_p4, ph6_p4, ph6_p6, ph8_p4, ph8_p6, ph10_p4, ph10_p6, ab_2, pc_104 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

        const auto h4_p4 = ph4_p4[k];
        const auto h6_p4 = ph6_p4[k];
        const auto h6_p6 = ph6_p6[k];
        const auto h8_p4 = ph8_p4[k];
        const auto h8_p6 = ph8_p6[k];
        const auto h10_p4 = ph10_p4[k];
        const auto h10_p6 = ph10_p6[k];

        pc_104[k] = e_1 * fs_765_4_11 * h4_p4 - e_2 * fs_285_44_55 * h6_p4 - e_2 * fs_285_8_30 * h6_p6 - e_2 * fs_1530_11_11 * r_2 * h4_p4 - e_3 * fs_357_26_5 * h8_p4 + e_3 * fs_21_52_2730 * h8_p6 + e_3 * fs_57_22_55 * r_2 * h6_p4 + e_3 * fs_57_4_30 * r_2 * h6_p6 + e_3 * fs_4590_143_11 * r_4 * h4_p4 - e_4 * fs_322_4199_273 * h10_p4 + e_4 * fs_644_4199_546 * h10_p6 + e_4 * fs_714_247_5 * r_2 * h8_p4 - e_4 * fs_21_247_2730 * r_2 * h8_p6 - e_4 * fs_57_187_55 * r_4 * h6_p4 - e_4 * fs_57_34_30 * r_4 * h6_p6 - e_4 * fs_408_143_11 * r_6 * h4_p4 + e_5 * fs_28_4199_273 * r_2 * h10_p4 - e_5 * fs_56_4199_546 * r_2 * h10_p6 - e_5 * fs_34_247_5 * r_4 * h8_p4 + e_5 * fs_1_247_2730 * r_4 * h8_p6 + e_5 * fs_2_187_55 * r_6 * h6_p4 + e_5 * fs_1_17_30 * r_6 * h6_p6 + e_5 * fs_12_143_11 * r_8 * h4_p4;
    }

    // NOTE: the angular components are formed in 93 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph4_p3, ph6_p3, ph8_p3, ph8_p7, ph10_p3, ph10_p7, ab_2, pc_105 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

        const auto h4_p3 = ph4_p3[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h8_p3 = ph8_p3[k];
        const auto h8_p7 = ph8_p7[k];
        const auto h10_p3 = ph10_p3[k];
        const auto h10_p7 = ph10_p7[k];

        pc_105[k] = - e_1 * fs_765_8_11 * h4_p3 + e_2 * fs_285_11_33 * h6_p3 + e_2 * fs_765_11_11 * r_2 * h4_p3 + e_3 * fs_483_52_6 * h8_p3 - e_3 * fs_21_52_182 * h8_p7 - e_3 * fs_114_11_33 * r_2 * h6_p3 - e_3 * fs_2295_143_11 * r_4 * h4_p3 + e_4 * fs_161_4199_273 * h10_p3 + e_4 * fs_322_4199_4641 * h10_p7 - e_4 * fs_483_247_6 * r_2 * h8_p3 + e_4 * fs_21_247_182 * r_2 * h8_p7 + e_4 * fs_228_187_33 * r_4 * h6_p3 + e_4 * fs_204_143_11 * r_6 * h4_p3 - e_5 * fs_14_4199_273 * r_2 * h10_p3 - e_5 * fs_28_4199_4641 * r_2 * h10_p7 + e_5 * fs_23_247_6 * r_4 * h8_p3 - e_5 * fs_1_247_182 * r_4 * h8_p7 - e_5 * fs_8_187_33 * r_6 * h6_p3 - e_5 * fs_6_143_11 * r_8 * h4_p3;
    }

    // NOTE: the angular components are formed in 93 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph2_p2, ph4_p2, ph6_p2, ph8_p2, ph8_p8, ph10_p2, ph10_p8, ab_2, pc_106 : simd::cache_line_size())
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
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;

        const auto h2_p2 = ph2_p2[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h8_p2 = ph8_p2[k];
        const auto h8_p8 = ph8_p8[k];
        const auto h10_p2 = ph10_p2[k];
        const auto h10_p8 = ph10_p8[k];

        pc_106[k] = e_0 * fs_1575_16_165 * h2_p2 - e_1 * fs_765_8_11 * h4_p2 - e_1 * fs_1125_8_165 * r_2 * h2_p2 - e_2 * fs_665_88_462 * h6_p2 + e_2 * fs_765_11_11 * r_2 * h4_p2 + e_2 * fs_125_2_165 * r_4 * h2_p2 - e_3 * fs_609_572_154 * h8_p2 - e_3 * fs_147_13_13 * h8_p8 + e_3 * fs_133_44_462 * r_2 * h6_p2 - e_3 * fs_2295_143_11 * r_4 * h4_p2 - e_3 * fs_125_11_165 * r_6 * h2_p2 - e_4 * fs_644_4199_3 * h10_p2 + e_4 * fs_1932_4199_221 * h10_p8 + e_4 * fs_609_2717_154 * r_2 * h8_p2 + e_4 * fs_588_247_13 * r_2 * h8_p8 - e_4 * fs_133_374_462 * r_4 * h6_p2 + e_4 * fs_204_143_11 * r_6 * h4_p2 + e_4 * fs_125_143_165 * r_8 * h2_p2 + e_5 * fs_56_4199_3 * r_2 * h10_p2 - e_5 * fs_168_4199_221 * r_2 * h10_p8 - e_5 * fs_29_2717_154 * r_4 * h8_p2 - e_5 * fs_28_247_13 * r_4 * h8_p8 + e_5 * fs_7_561_462 * r_6 * h6_p2 - e_5 * fs_6_143_11 * r_8 * h4_p2 - e_5 * fs_10_429_165 * r_10 * h2_p2;
    }

    // NOTE: the angular components are formed in 93 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph2_p1, ph4_p1, ph6_p1, ph8_p1, ph10_p1, ph10_p9, ab_2, pc_107 : simd::cache_line_size())
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
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;

        const auto h2_p1 = ph2_p1[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h10_p1 = ph10_p1[k];
        const auto h10_p9 = ph10_p9[k];

        pc_107[k] = e_0 * fs_1575_32_330 * h2_p1 + e_1 * fs_765_4_11 * h4_p1 - e_1 * fs_1125_16_330 * r_2 * h2_p1 + e_2 * fs_95_44_2310 * h6_p1 - e_2 * fs_1530_11_11 * r_2 * h4_p1 + e_2 * fs_125_4_330 * r_4 * h2_p1 + e_3 * fs_147_286_110 * h8_p1 - e_3 * fs_19_22_2310 * r_2 * h6_p1 + e_3 * fs_4590_143_11 * r_4 * h4_p1 - e_3 * fs_125_22_330 * r_6 * h2_p1 + e_4 * fs_483_8398_2 * h10_p1 + e_4 * fs_483_4199_4199 * h10_p9 - e_4 * fs_294_2717_110 * r_2 * h8_p1 + e_4 * fs_19_187_2310 * r_4 * h6_p1 - e_4 * fs_408_143_11 * r_6 * h4_p1 + e_4 * fs_125_286_330 * r_8 * h2_p1 - e_5 * fs_21_4199_2 * r_2 * h10_p1 - e_5 * fs_42_4199_4199 * r_2 * h10_p9 + e_5 * fs_14_2717_110 * r_4 * h8_p1 - e_5 * fs_2_561_2310 * r_6 * h6_p1 + e_5 * fs_12_143_11 * r_8 * h4_p1 - e_5 * fs_5_429_330 * r_10 * h2_p1;
    }

    // NOTE: the angular components are formed in 93 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph2_m2, ph4_m3, ph4_m2, ph6_m3, ph6_m2, ph8_m3, ph8_m2, ph10_m10, ph10_m9, ph10_m3, ph10_m2, ab_2, pc_108, pc_109 : simd::cache_line_size())
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
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;

        const auto h2_m2 = ph2_m2[k];
        const auto h4_m3 = ph4_m3[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h8_m3 = ph8_m3[k];
        const auto h8_m2 = ph8_m2[k];
        const auto h10_m10 = ph10_m10[k];
        const auto h10_m9 = ph10_m9[k];
        const auto h10_m3 = ph10_m3[k];
        const auto h10_m2 = ph10_m2[k];

        pc_108[k] = - e_0 * fs_4725_32_110 * h2_m2 - e_1 * fs_255_4_66 * h4_m2 + e_1 * fs_3375_16_110 * r_2 * h2_m2 - e_2 * fs_285_44_77 * h6_m2 + e_2 * fs_510_11_66 * r_2 * h4_m2 - e_2 * fs_375_4_110 * r_4 * h2_m2 - e_3 * fs_21_143_231 * h8_m2 + e_3 * fs_57_22_77 * r_2 * h6_m2 - e_3 * fs_1530_143_66 * r_4 * h4_m2 + e_3 * fs_375_22_110 * r_6 * h2_m2 + e_4 * fs_161_4199_62985 * h10_m10 - e_4 * fs_161_8398_2 * h10_m2 + e_4 * fs_84_2717_231 * r_2 * h8_m2 - e_4 * fs_57_187_77 * r_4 * h6_m2 + e_4 * fs_136_143_66 * r_6 * h4_m2 - e_4 * fs_375_286_110 * r_8 * h2_m2 - e_5 * fs_14_4199_62985 * r_2 * h10_m10 + e_5 * fs_7_4199_2 * r_2 * h10_m2 - e_5 * fs_4_2717_231 * r_4 * h8_m2 + e_5 * fs_2_187_77 * r_6 * h6_m2 - e_5 * fs_4_143_66 * r_8 * h4_m2 + e_5 * fs_5_143_110 * r_10 * h2_m2;

        pc_109[k] = e_1 * fs_255_8_462 * h4_m3 + e_2 * fs_855_88_154 * h6_m3 - e_2 * fs_255_11_462 * r_2 * h4_m3 + e_3 * fs_63_26_7 * h8_m3 - e_3 * fs_171_44_154 * r_2 * h6_m3 + e_3 * fs_765_143_462 * r_4 * h4_m3 + e_4 * fs_161_4199_25194 * h10_m9 + e_4 * fs_161_8398_26 * h10_m3 - e_4 * fs_126_247_7 * r_2 * h8_m3 + e_4 * fs_171_374_154 * r_4 * h6_m3 - e_4 * fs_68_143_462 * r_6 * h4_m3 - e_5 * fs_14_4199_25194 * r_2 * h10_m9 - e_5 * fs_7_4199_26 * r_2 * h10_m3 + e_5 * fs_6_247_7 * r_4 * h8_m3 - e_5 * fs_3_187_154 * r_6 * h6_m3 + e_5 * fs_2_143_462 * r_8 * h4_m3;
    }

    // NOTE: the angular components are formed in 93 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph4_m4, ph6_m5, ph6_m4, ph8_m8, ph8_m7, ph8_m5, ph8_m4, ph10_m8, ph10_m7, ph10_m5, ph10_m4, ab_2, pc_110, pc_111 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

        const auto h4_m4 = ph4_m4[k];
        const auto h6_m5 = ph6_m5[k];
        const auto h6_m4 = ph6_m4[k];
        const auto h8_m8 = ph8_m8[k];
        const auto h8_m7 = ph8_m7[k];
        const auto h8_m5 = ph8_m5[k];
        const auto h8_m4 = ph8_m4[k];
        const auto h10_m8 = ph10_m8[k];
        const auto h10_m7 = ph10_m7[k];
        const auto h10_m5 = ph10_m5[k];
        const auto h10_m4 = ph10_m4[k];

        pc_110[k] = - e_1 * fs_255_4_66 * h4_m4 - e_2 * fs_855_88_330 * h6_m4 + e_2 * fs_510_11_66 * r_2 * h4_m4 + e_3 * fs_21_13_546 * h8_m8 - e_3 * fs_63_26_30 * h8_m4 + e_3 * fs_171_44_330 * r_2 * h6_m4 - e_3 * fs_1530_143_66 * r_4 * h4_m4 + e_4 * fs_161_4199_9282 * h10_m8 - e_4 * fs_161_8398_182 * h10_m4 - e_4 * fs_84_247_546 * r_2 * h8_m8 + e_4 * fs_126_247_30 * r_2 * h8_m4 - e_4 * fs_171_374_330 * r_4 * h6_m4 + e_4 * fs_136_143_66 * r_6 * h4_m4 - e_5 * fs_14_4199_9282 * r_2 * h10_m8 + e_5 * fs_7_4199_182 * r_2 * h10_m4 + e_5 * fs_4_247_546 * r_4 * h8_m8 - e_5 * fs_6_247_30 * r_4 * h8_m4 + e_5 * fs_3_187_330 * r_6 * h6_m4 - e_5 * fs_4_143_66 * r_8 * h4_m4;

        pc_111[k] = e_2 * fs_285_8_30 * h6_m5 + e_3 * fs_63_26_273 * h8_m7 + e_3 * fs_21_13_195 * h8_m5 - e_3 * fs_57_4_30 * r_2 * h6_m5 + e_4 * fs_161_4199_3094 * h10_m7 + e_4 * fs_161_8398_910 * h10_m5 - e_4 * fs_126_247_273 * r_2 * h8_m7 - e_4 * fs_84_247_195 * r_2 * h8_m5 + e_4 * fs_57_34_30 * r_4 * h6_m5 - e_5 * fs_14_4199_3094 * r_2 * h10_m7 - e_5 * fs_7_4199_910 * r_2 * h10_m5 + e_5 * fs_6_247_273 * r_4 * h8_m7 + e_5 * fs_4_247_195 * r_4 * h8_m5 - e_5 * fs_1_17_30 * r_6 * h6_m5;
    }

    // NOTE: the angular components are formed in 93 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_2, pe_3, pe_4, pe_5, ph6_p5, ph6_p6, ph8_p5, ph8_p6, ph8_p7, ph10_p5, ph10_p6, ph10_p7, ab_2, pc_112, pc_113 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h6_p5 = ph6_p5[k];
        const auto h6_p6 = ph6_p6[k];
        const auto h8_p5 = ph8_p5[k];
        const auto h8_p6 = ph8_p6[k];
        const auto h8_p7 = ph8_p7[k];
        const auto h10_p5 = ph10_p5[k];
        const auto h10_p6 = ph10_p6[k];
        const auto h10_p7 = ph10_p7[k];

        pc_112[k] = e_2 * f_855_4 * h6_p6 + e_3 * fs_63_13_91 * h8_p6 - e_3 * f_171_2 * r_2 * h6_p6 + e_4 * fs_322_4199_455 * h10_p6 - e_4 * fs_252_247_91 * r_2 * h8_p6 + e_4 * f_171_17 * r_4 * h6_p6 - e_5 * fs_28_4199_455 * r_2 * h10_p6 + e_5 * fs_12_247_91 * r_4 * h8_p6 - e_5 * f_6_17 * r_6 * h6_p6;

        pc_113[k] = - e_2 * fs_285_8_30 * h6_p5 - e_3 * fs_21_13_195 * h8_p5 + e_3 * fs_63_26_273 * h8_p7 + e_3 * fs_57_4_30 * r_2 * h6_p5 - e_4 * fs_161_8398_910 * h10_p5 + e_4 * fs_161_4199_3094 * h10_p7 + e_4 * fs_84_247_195 * r_2 * h8_p5 - e_4 * fs_126_247_273 * r_2 * h8_p7 - e_4 * fs_57_34_30 * r_4 * h6_p5 + e_5 * fs_7_4199_910 * r_2 * h10_p5 - e_5 * fs_14_4199_3094 * r_2 * h10_p7 - e_5 * fs_4_247_195 * r_4 * h8_p5 + e_5 * fs_6_247_273 * r_4 * h8_p7 + e_5 * fs_1_17_30 * r_6 * h6_p5;
    }

    // NOTE: the angular components are formed in 93 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph4_p3, ph4_p4, ph6_p3, ph6_p4, ph8_p3, ph8_p4, ph8_p8, ph10_p3, ph10_p4, ph10_p8, ph10_p9, ab_2, pc_114, pc_115 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

        const auto h4_p3 = ph4_p3[k];
        const auto h4_p4 = ph4_p4[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h6_p4 = ph6_p4[k];
        const auto h8_p3 = ph8_p3[k];
        const auto h8_p4 = ph8_p4[k];
        const auto h8_p8 = ph8_p8[k];
        const auto h10_p3 = ph10_p3[k];
        const auto h10_p4 = ph10_p4[k];
        const auto h10_p8 = ph10_p8[k];
        const auto h10_p9 = ph10_p9[k];

        pc_114[k] = e_1 * fs_255_4_66 * h4_p4 + e_2 * fs_855_88_330 * h6_p4 - e_2 * fs_510_11_66 * r_2 * h4_p4 + e_3 * fs_63_26_30 * h8_p4 + e_3 * fs_21_13_546 * h8_p8 - e_3 * fs_171_44_330 * r_2 * h6_p4 + e_3 * fs_1530_143_66 * r_4 * h4_p4 + e_4 * fs_161_8398_182 * h10_p4 + e_4 * fs_161_4199_9282 * h10_p8 - e_4 * fs_126_247_30 * r_2 * h8_p4 - e_4 * fs_84_247_546 * r_2 * h8_p8 + e_4 * fs_171_374_330 * r_4 * h6_p4 - e_4 * fs_136_143_66 * r_6 * h4_p4 - e_5 * fs_7_4199_182 * r_2 * h10_p4 - e_5 * fs_14_4199_9282 * r_2 * h10_p8 + e_5 * fs_6_247_30 * r_4 * h8_p4 + e_5 * fs_4_247_546 * r_4 * h8_p8 - e_5 * fs_3_187_330 * r_6 * h6_p4 + e_5 * fs_4_143_66 * r_8 * h4_p4;

        pc_115[k] = - e_1 * fs_255_8_462 * h4_p3 - e_2 * fs_855_88_154 * h6_p3 + e_2 * fs_255_11_462 * r_2 * h4_p3 - e_3 * fs_63_26_7 * h8_p3 + e_3 * fs_171_44_154 * r_2 * h6_p3 - e_3 * fs_765_143_462 * r_4 * h4_p3 - e_4 * fs_161_8398_26 * h10_p3 + e_4 * fs_161_4199_25194 * h10_p9 + e_4 * fs_126_247_7 * r_2 * h8_p3 - e_4 * fs_171_374_154 * r_4 * h6_p3 + e_4 * fs_68_143_462 * r_6 * h4_p3 + e_5 * fs_7_4199_26 * r_2 * h10_p3 - e_5 * fs_14_4199_25194 * r_2 * h10_p9 - e_5 * fs_6_247_7 * r_4 * h8_p3 + e_5 * fs_3_187_154 * r_6 * h6_p3 - e_5 * fs_2_143_462 * r_8 * h4_p3;
    }

    // NOTE: the angular components are formed in 93 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph2_p2, ph4_p2, ph6_p2, ph8_p2, ph10_p2, ph10_p10, ab_2, pc_116 : simd::cache_line_size())
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
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;

        const auto h2_p2 = ph2_p2[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h8_p2 = ph8_p2[k];
        const auto h10_p2 = ph10_p2[k];
        const auto h10_p10 = ph10_p10[k];

        pc_116[k] = e_0 * fs_4725_32_110 * h2_p2 + e_1 * fs_255_4_66 * h4_p2 - e_1 * fs_3375_16_110 * r_2 * h2_p2 + e_2 * fs_285_44_77 * h6_p2 - e_2 * fs_510_11_66 * r_2 * h4_p2 + e_2 * fs_375_4_110 * r_4 * h2_p2 + e_3 * fs_21_143_231 * h8_p2 - e_3 * fs_57_22_77 * r_2 * h6_p2 + e_3 * fs_1530_143_66 * r_4 * h4_p2 - e_3 * fs_375_22_110 * r_6 * h2_p2 + e_4 * fs_161_8398_2 * h10_p2 + e_4 * fs_161_4199_62985 * h10_p10 - e_4 * fs_84_2717_231 * r_2 * h8_p2 + e_4 * fs_57_187_77 * r_4 * h6_p2 - e_4 * fs_136_143_66 * r_6 * h4_p2 + e_4 * fs_375_286_110 * r_8 * h2_p2 - e_5 * fs_7_4199_2 * r_2 * h10_p2 - e_5 * fs_14_4199_62985 * r_2 * h10_p10 + e_5 * fs_4_2717_231 * r_4 * h8_p2 - e_5 * fs_2_187_77 * r_6 * h6_p2 + e_5 * fs_4_143_66 * r_8 * h4_p2 - e_5 * fs_5_143_110 * r_10 * h2_p2;
    }

    // NOTE: the values of a combination of angular components are stored as one
    // row of nvalues columns, with the component on bra side running slowest. The
    // rows which the symmetry relates to an already formed one are copied from it,
    // and the atom pairs beyond the reach of every pair of primitives are set to
    // zero.

    const size_t sources[117] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116};

    for (size_t n = 0; n < 117; n++)
    {
        if (sources[n] != n) std::copy(values + sources[n] * nvalues, values + sources[n] * nvalues + nmax, values + n * nvalues);

        std::fill(values + n * nvalues + nmax, values + (n + 1) * nvalues, 0.0);
    }
}

}  // namespace simdkin
