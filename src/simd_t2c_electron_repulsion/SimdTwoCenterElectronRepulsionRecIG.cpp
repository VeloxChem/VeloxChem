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



#include "SimdTwoCenterElectronRepulsionRecIG.hpp"

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

#include "ErrorHandler.hpp"
#include "MathConst.hpp"
#include "SimdAlign.hpp"
#include "SimdBoysFunc.hpp"
#include "SimdVariableMatrix.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_ig_electron_repulsion(double                         *values,
                              const size_t                    nvalues,
                              const CBasisFunction           &bra,
                              const CBasisFunction           &ket,
                              const std::vector<CSimdMatrix> &harmonics,
                              const CSimdMatrix              &coordinates) -> void
{
    if ((bra.get_angular_momentum() != 6) || (ket.get_angular_momentum() != 4))
    {
        errors::assertMsgCritical(
            false, std::string("SimdTwoCenterElectronRepulsionFunc.compute_ig_electron_repulsion: Basis functions must be of angular momenta six and four"));
    }

    if (harmonics.size() < 10)
    {
        errors::assertMsgCritical(
            false, std::string("SimdTwoCenterElectronRepulsionFunc.compute_ig_electron_repulsion: Harmonics must reach angular momentum 10"));
    }

    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("SimdTwoCenterElectronRepulsionFunc.compute_ig_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nprims = nprim_a * nprim_b;

    const auto *ab_2 = coordinates.data(6);

    // NOTE: the pairs of primitives are not screened and every one of them
    // reaches every atom pair, as the Coulomb operator decays as the inverse of
    // the interatomic distance, so every row spans the atom pairs of the block.

    auto boys = CSimdVariableMatrix(std::vector<size_t>(nprims, nvalues), 12);

    // set up the arguments of the Boys function of each pair of primitives

    for (size_t i = 0; i < nprim_a; i++)
    {
        const auto aexp = a_exps[i];

        for (size_t j = 0; j < nprim_b; j++)
        {
            const auto fmu = aexp * b_exps[j] / (aexp + b_exps[j]);

            auto *bargs = boys.data(0, i * nprim_b + j);

#pragma omp simd aligned(bargs, ab_2 : simd::cache_line_size())
            for (size_t k = 0; k < nvalues; k++)
            {
                bargs[k] = fmu * ab_2[k];
            }
        }
    }

    // NOTE: the Boys function of every pair of primitives is computed by one
    // call, which fills the orders 0 to 10 of every row. The terms read the
    // orders 6 to 10 alone, and the orders below them are formed on the
    // way to them by the recursion the Boys function is evaluated with.

    simdfunc::compute_boys_function(boys);

    // NOTE: the buffer holds the contracted prefactor of each exponent factor of
    // the terms, as the integrals of the angular components are formed straight
    // into the values and are not written a second time. Every exponent factor is
    // used with one order of the Boys function alone, so the order follows from
    // the factor and one accumulator per factor suffices.

    auto buffer = CSimdMatrix(5, nvalues);

    auto *pe_0 = buffer.data(0);
    auto *pe_1 = buffer.data(1);
    auto *pe_2 = buffer.data(2);
    auto *pe_3 = buffer.data(3);
    auto *pe_4 = buffer.data(4);

    std::fill(pe_0, pe_0 + nvalues, 0.0);
    std::fill(pe_1, pe_1 + nvalues, 0.0);
    std::fill(pe_2, pe_2 + nvalues, 0.0);
    std::fill(pe_3, pe_3 + nvalues, 0.0);
    std::fill(pe_4, pe_4 + nvalues, 0.0);

    constexpr auto fpi = mathconst::pi_value();

    // NOTE: the two-center repulsion of a pair of s type primitives is two pi to
    // the five halves over the exponents times the square root of their sum,
    // times the Boys function of the order the terms ask for.

    const auto fcoul = 2.0 * fpi * fpi * std::sqrt(fpi);

    // accumulate the prefactor of each exponent factor over the pairs of primitives

    for (size_t i = 0; i < nprim_a; i++)
    {
        const auto aexp = a_exps[i];

        const auto anorm = a_norms[i];

        for (size_t j = 0; j < nprim_b; j++)
        {
            const auto bexp = b_exps[j];

            const auto fexp = aexp + bexp;

            const auto fmu = aexp * bexp / fexp;

            const auto fbase = fcoul * anorm * b_norms[j] / (aexp * bexp * std::sqrt(fexp));

            const auto ff_0 = fbase * bexp * bexp / (fexp * fexp * fexp * fexp * fexp * fexp);

            const auto ff_1 = fbase * bexp * bexp * fmu / (fexp * fexp * fexp * fexp * fexp * fexp);

            const auto ff_2 = fbase * bexp * bexp * fmu * fmu / (fexp * fexp * fexp * fexp * fexp * fexp);

            const auto ff_3 = fbase * bexp * bexp * fmu * fmu * fmu / (fexp * fexp * fexp * fexp * fexp * fexp);

            const auto ff_4 = fbase * bexp * bexp * fmu * fmu * fmu * fmu / (fexp * fexp * fexp * fexp * fexp * fexp);

            const auto *bv_0 = boys.data(7, i * nprim_b + j);

            const auto *bv_1 = boys.data(8, i * nprim_b + j);

            const auto *bv_2 = boys.data(9, i * nprim_b + j);

            const auto *bv_3 = boys.data(10, i * nprim_b + j);

            const auto *bv_4 = boys.data(11, i * nprim_b + j);

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, bv_0, bv_1, bv_2, bv_3, bv_4 : simd::cache_line_size())
            for (size_t k = 0; k < nvalues; k++)
            {
                pe_0[k] += ff_0 * bv_0[k];
                pe_1[k] += ff_1 * bv_1[k];
                pe_2[k] += ff_2 * bv_2[k];
                pe_3[k] += ff_3 * bv_3[k];
                pe_4[k] += ff_4 * bv_4[k];
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

    const auto f_105_11 = 105.0 / 11.0;
    const auto f_105_4 = 105.0 / 4.0;
    const auto f_10_1 = 10.0;
    const auto f_120_11 = 120.0 / 11.0;
    const auto f_135_22 = 135.0 / 22.0;
    const auto f_135_286 = 135.0 / 286.0;
    const auto f_135_44 = 135.0 / 44.0;
    const auto f_1470_4199 = 1470.0 / 4199.0;
    const auto f_1575_16 = 1575.0 / 16.0;
    const auto f_15_2 = 15.0 / 2.0;
    const auto f_15_8 = 15.0 / 8.0;
    const auto f_165_4 = 165.0 / 4.0;
    const auto f_17_143 = 17.0 / 143.0;
    const auto f_18_11 = 18.0 / 11.0;
    const auto f_18_187 = 18.0 / 187.0;
    const auto f_1_143 = 1.0 / 143.0;
    const auto f_1_3 = 1.0 / 3.0;
    const auto f_1_51 = 1.0 / 51.0;
    const auto f_20_143 = 20.0 / 143.0;
    const auto f_20_429 = 20.0 / 429.0;
    const auto f_225_11 = 225.0 / 11.0;
    const auto f_225_2 = 225.0 / 2.0;
    const auto f_252_143 = 252.0 / 143.0;
    const auto f_255_8 = 255.0 / 8.0;
    const auto f_25_143 = 25.0 / 143.0;
    const auto f_265_44 = 265.0 / 44.0;
    const auto f_28_11 = 28.0 / 11.0;
    const auto f_28_187 = 28.0 / 187.0;
    const auto f_2_1 = 2.0;
    const auto f_2_13 = 2.0 / 13.0;
    const auto f_2_17 = 2.0 / 17.0;
    const auto f_30_1 = 30.0;
    const auto f_32_11 = 32.0 / 11.0;
    const auto f_32_187 = 32.0 / 187.0;
    const auto f_3_1 = 3.0;
    const auto f_3_143 = 3.0 / 143.0;
    const auto f_3_17 = 3.0 / 17.0;
    const auto f_40_33 = 40.0 / 33.0;
    const auto f_450_143 = 450.0 / 143.0;
    const auto f_45_11 = 45.0 / 11.0;
    const auto f_45_13 = 45.0 / 13.0;
    const auto f_45_2 = 45.0 / 2.0;
    const auto f_45_286 = 45.0 / 286.0;
    const auto f_45_4 = 45.0 / 4.0;
    const auto f_45_44 = 45.0 / 44.0;
    const auto f_45_8 = 45.0 / 8.0;
    const auto f_4_143 = 4.0 / 143.0;
    const auto f_504_2717 = 504.0 / 2717.0;
    const auto f_50_11 = 50.0 / 11.0;
    const auto f_53_33 = 53.0 / 33.0;
    const auto f_53_561 = 53.0 / 561.0;
    const auto f_5_4 = 5.0 / 4.0;
    const auto f_64_33 = 64.0 / 33.0;
    const auto f_64_561 = 64.0 / 561.0;
    const auto f_75_2 = 75.0 / 2.0;
    const auto f_765_286 = 765.0 / 286.0;
    const auto f_765_44 = 765.0 / 44.0;
    const auto f_80_11 = 80.0 / 11.0;
    const auto f_90_143 = 90.0 / 143.0;
    const auto fs_101_5434_66 = std::sqrt(30603.0 / 1342198.0);
    const auto fs_101_572_66 = std::sqrt(30603.0 / 14872.0);
    const auto fs_105_16_105 = std::sqrt(1157625.0 / 256.0);
    const auto fs_105_16_15 = std::sqrt(165375.0 / 256.0);
    const auto fs_105_16_165 = std::sqrt(1819125.0 / 256.0);
    const auto fs_105_16_210 = std::sqrt(1157625.0 / 128.0);
    const auto fs_105_16_35 = std::sqrt(385875.0 / 256.0);
    const auto fs_105_16_70 = std::sqrt(385875.0 / 128.0);
    const auto fs_105_22_3 = std::sqrt(33075.0 / 484.0);
    const auto fs_105_22_5 = std::sqrt(55125.0 / 484.0);
    const auto fs_105_2717_10 = std::sqrt(110250.0 / 7382089.0);
    const auto fs_105_286_10 = std::sqrt(55125.0 / 40898.0);
    const auto fs_105_32_10 = std::sqrt(55125.0 / 512.0);
    const auto fs_105_32_2 = std::sqrt(11025.0 / 512.0);
    const auto fs_105_32_330 = std::sqrt(1819125.0 / 512.0);
    const auto fs_105_32_70 = std::sqrt(385875.0 / 512.0);
    const auto fs_105_4199_154 = std::sqrt(1697850.0 / 17631601.0);
    const auto fs_105_4199_42 = std::sqrt(463050.0 / 17631601.0);
    const auto fs_105_44_14 = std::sqrt(77175.0 / 968.0);
    const auto fs_105_44_3 = std::sqrt(33075.0 / 1936.0);
    const auto fs_105_4_7 = std::sqrt(77175.0 / 16.0);
    const auto fs_105_8398_286 = std::sqrt(121275.0 / 2712554.0);
    const auto fs_105_8_30 = std::sqrt(165375.0 / 32.0);
    const auto fs_105_8_42 = std::sqrt(231525.0 / 32.0);
    const auto fs_10_11_14 = std::sqrt(1400.0 / 121.0);
    const auto fs_10_11_5 = std::sqrt(500.0 / 121.0);
    const auto fs_10_11_7 = std::sqrt(700.0 / 121.0);
    const auto fs_10_143_3 = std::sqrt(300.0 / 20449.0);
    const auto fs_10_1_7 = std::sqrt(700.0);
    const auto fs_10_209_3 = std::sqrt(300.0 / 43681.0);
    const auto fs_10_33_105 = std::sqrt(3500.0 / 363.0);
    const auto fs_10_33_15 = std::sqrt(500.0 / 363.0);
    const auto fs_10_33_165 = std::sqrt(500.0 / 33.0);
    const auto fs_10_33_210 = std::sqrt(7000.0 / 363.0);
    const auto fs_10_33_35 = std::sqrt(3500.0 / 1089.0);
    const auto fs_10_33_70 = std::sqrt(7000.0 / 1089.0);
    const auto fs_10_429_30 = std::sqrt(1000.0 / 61347.0);
    const auto fs_10_429_42 = std::sqrt(1400.0 / 61347.0);
    const auto fs_111_5434_30 = std::sqrt(184815.0 / 14764178.0);
    const auto fs_111_572_30 = std::sqrt(184815.0 / 163592.0);
    const auto fs_115_88_6 = std::sqrt(39675.0 / 3872.0);
    const auto fs_116_2717_11 = std::sqrt(13456.0 / 671099.0);
    const auto fs_119_143_3 = std::sqrt(42483.0 / 20449.0);
    const auto fs_11_247_26 = std::sqrt(242.0 / 4693.0);
    const auto fs_11_26_26 = std::sqrt(121.0 / 26.0);
    const auto fs_127_5434_42 = std::sqrt(338709.0 / 14764178.0);
    const auto fs_127_572_42 = std::sqrt(338709.0 / 163592.0);
    const auto fs_12_143_14 = std::sqrt(2016.0 / 20449.0);
    const auto fs_12_143_5 = std::sqrt(720.0 / 20449.0);
    const auto fs_135_11_5 = std::sqrt(91125.0 / 121.0);
    const auto fs_135_143_11 = std::sqrt(18225.0 / 1859.0);
    const auto fs_135_143_7 = std::sqrt(127575.0 / 20449.0);
    const auto fs_135_22_11 = std::sqrt(18225.0 / 44.0);
    const auto fs_135_22_7 = std::sqrt(127575.0 / 484.0);
    const auto fs_135_286_11 = std::sqrt(18225.0 / 7436.0);
    const auto fs_135_286_30 = std::sqrt(273375.0 / 40898.0);
    const auto fs_135_44_11 = std::sqrt(18225.0 / 176.0);
    const auto fs_135_44_30 = std::sqrt(273375.0 / 968.0);
    const auto fs_135_572_42 = std::sqrt(382725.0 / 163592.0);
    const auto fs_135_88_42 = std::sqrt(382725.0 / 3872.0);
    const auto fs_139_143_2 = std::sqrt(38642.0 / 20449.0);
    const auto fs_13_22_10 = std::sqrt(845.0 / 242.0);
    const auto fs_13_374_10 = std::sqrt(845.0 / 69938.0);
    const auto fs_140_4199_3 = std::sqrt(58800.0 / 17631601.0);
    const auto fs_147_143_3 = std::sqrt(64827.0 / 20449.0);
    const auto fs_14_11_3 = std::sqrt(588.0 / 121.0);
    const auto fs_14_11_5 = std::sqrt(980.0 / 121.0);
    const auto fs_14_143_105 = std::sqrt(20580.0 / 20449.0);
    const auto fs_14_143_143 = std::sqrt(196.0 / 143.0);
    const auto fs_14_187_3 = std::sqrt(588.0 / 34969.0);
    const auto fs_14_187_5 = std::sqrt(980.0 / 34969.0);
    const auto fs_14_247_13 = std::sqrt(196.0 / 4693.0);
    const auto fs_14_33_35 = std::sqrt(6860.0 / 1089.0);
    const auto fs_14_4199_231 = std::sqrt(45276.0 / 17631601.0);
    const auto fs_14_4199_273 = std::sqrt(4116.0 / 1356277.0);
    const auto fs_14_4199_3003 = std::sqrt(45276.0 / 1356277.0);
    const auto fs_14_4199_455 = std::sqrt(6860.0 / 1356277.0);
    const auto fs_14_4199_4641 = std::sqrt(4116.0 / 79781.0);
    const auto fs_14_4199_5005 = std::sqrt(75460.0 / 1356277.0);
    const auto fs_14_4199_6006 = std::sqrt(90552.0 / 1356277.0);
    const auto fs_14_4199_7293 = std::sqrt(6468.0 / 79781.0);
    const auto fs_14_561_35 = std::sqrt(6860.0 / 314721.0);
    const auto fs_15_11_2 = std::sqrt(450.0 / 121.0);
    const auto fs_15_11_33 = std::sqrt(675.0 / 11.0);
    const auto fs_15_1_3 = std::sqrt(675.0);
    const auto fs_15_1_30 = std::sqrt(6750.0);
    const auto fs_15_1_42 = std::sqrt(9450.0);
    const auto fs_15_2717_22 = std::sqrt(450.0 / 671099.0);
    const auto fs_15_286_22 = std::sqrt(225.0 / 3718.0);
    const auto fs_15_2_105 = std::sqrt(23625.0 / 4.0);
    const auto fs_15_2_14 = std::sqrt(1575.0 / 2.0);
    const auto fs_15_2_15 = std::sqrt(3375.0 / 4.0);
    const auto fs_15_2_165 = std::sqrt(37125.0 / 4.0);
    const auto fs_15_2_210 = std::sqrt(23625.0 / 2.0);
    const auto fs_15_2_35 = std::sqrt(7875.0 / 4.0);
    const auto fs_15_2_5 = std::sqrt(1125.0 / 4.0);
    const auto fs_15_2_7 = std::sqrt(1575.0 / 4.0);
    const auto fs_15_2_70 = std::sqrt(7875.0 / 2.0);
    const auto fs_15_44_55 = std::sqrt(1125.0 / 176.0);
    const auto fs_15_44_77 = std::sqrt(1575.0 / 176.0);
    const auto fs_15_4_10 = std::sqrt(1125.0 / 8.0);
    const auto fs_15_4_110 = std::sqrt(12375.0 / 8.0);
    const auto fs_15_4_2 = std::sqrt(225.0 / 8.0);
    const auto fs_15_4_30 = std::sqrt(3375.0 / 8.0);
    const auto fs_15_4_330 = std::sqrt(37125.0 / 8.0);
    const auto fs_15_4_66 = std::sqrt(7425.0 / 8.0);
    const auto fs_15_4_70 = std::sqrt(7875.0 / 8.0);
    const auto fs_15_8_105 = std::sqrt(23625.0 / 64.0);
    const auto fs_15_8_210 = std::sqrt(23625.0 / 32.0);
    const auto fs_15_8_30 = std::sqrt(3375.0 / 32.0);
    const auto fs_15_8_42 = std::sqrt(4725.0 / 32.0);
    const auto fs_15_8_462 = std::sqrt(51975.0 / 32.0);
    const auto fs_175_5434_10 = std::sqrt(153125.0 / 14764178.0);
    const auto fs_175_572_10 = std::sqrt(153125.0 / 163592.0);
    const auto fs_17_247_5 = std::sqrt(1445.0 / 61009.0);
    const auto fs_17_26_5 = std::sqrt(1445.0 / 676.0);
    const auto fs_195_16_6 = std::sqrt(114075.0 / 128.0);
    const auto fs_195_88_10 = std::sqrt(190125.0 / 3872.0);
    const auto fs_196_4199_33 = std::sqrt(1267728.0 / 17631601.0);
    const auto fs_199_2717_3 = std::sqrt(118803.0 / 7382089.0);
    const auto fs_199_286_3 = std::sqrt(118803.0 / 81796.0);
    const auto fs_1_11_55 = std::sqrt(5.0 / 11.0);
    const auto fs_1_11_77 = std::sqrt(7.0 / 11.0);
    const auto fs_1_13_195 = std::sqrt(15.0 / 13.0);
    const auto fs_1_13_546 = std::sqrt(42.0 / 13.0);
    const auto fs_1_143_105 = std::sqrt(105.0 / 20449.0);
    const auto fs_1_143_210 = std::sqrt(210.0 / 20449.0);
    const auto fs_1_143_231 = std::sqrt(21.0 / 1859.0);
    const auto fs_1_143_30030 = std::sqrt(210.0 / 143.0);
    const auto fs_1_143_42 = std::sqrt(42.0 / 20449.0);
    const auto fs_1_143_462 = std::sqrt(42.0 / 1859.0);
    const auto fs_1_187_55 = std::sqrt(5.0 / 3179.0);
    const auto fs_1_187_77 = std::sqrt(7.0 / 3179.0);
    const auto fs_1_22_6 = std::sqrt(3.0 / 242.0);
    const auto fs_1_2717_3003 = std::sqrt(21.0 / 51623.0);
    const auto fs_1_286_3003 = std::sqrt(21.0 / 572.0);
    const auto fs_1_2_30 = std::sqrt(15.0 / 2.0);
    const auto fs_1_33_2310 = std::sqrt(70.0 / 33.0);
    const auto fs_1_34_30 = std::sqrt(15.0 / 578.0);
    const auto fs_1_3_42 = std::sqrt(14.0 / 3.0);
    const auto fs_1_494_182 = std::sqrt(7.0 / 9386.0);
    const auto fs_1_494_2730 = std::sqrt(105.0 / 9386.0);
    const auto fs_1_51_42 = std::sqrt(14.0 / 867.0);
    const auto fs_1_52_182 = std::sqrt(7.0 / 104.0);
    const auto fs_1_52_2730 = std::sqrt(105.0 / 104.0);
    const auto fs_1_5434_1430 = std::sqrt(5.0 / 103246.0);
    const auto fs_1_561_2310 = std::sqrt(70.0 / 9537.0);
    const auto fs_1_572_1430 = std::sqrt(5.0 / 1144.0);
    const auto fs_20_11_3 = std::sqrt(1200.0 / 121.0);
    const auto fs_20_33_30 = std::sqrt(4000.0 / 363.0);
    const auto fs_20_33_42 = std::sqrt(5600.0 / 363.0);
    const auto fs_20_429_7 = std::sqrt(2800.0 / 184041.0);
    const auto fs_21_143_143 = std::sqrt(441.0 / 143.0);
    const auto fs_21_143_77 = std::sqrt(3087.0 / 1859.0);
    const auto fs_21_4199_1155 = std::sqrt(509355.0 / 17631601.0);
    const auto fs_21_4199_165 = std::sqrt(72765.0 / 17631601.0);
    const auto fs_21_4199_2002 = std::sqrt(67914.0 / 1356277.0);
    const auto fs_21_4199_22 = std::sqrt(9702.0 / 17631601.0);
    const auto fs_21_4199_2431 = std::sqrt(4851.0 / 79781.0);
    const auto fs_21_4199_4199 = std::sqrt(441.0 / 4199.0);
    const auto fs_21_4199_5 = std::sqrt(2205.0 / 17631601.0);
    const auto fs_21_4199_715 = std::sqrt(24255.0 / 1356277.0);
    const auto fs_21_5434_286 = std::sqrt(441.0 / 103246.0);
    const auto fs_21_572_286 = std::sqrt(441.0 / 1144.0);
    const auto fs_21_8398_10010 = std::sqrt(169785.0 / 2712554.0);
    const auto fs_21_8398_110 = std::sqrt(24255.0 / 35263202.0);
    const auto fs_21_8398_1430 = std::sqrt(24255.0 / 2712554.0);
    const auto fs_21_8398_2 = std::sqrt(441.0 / 35263202.0);
    const auto fs_21_8398_286 = std::sqrt(4851.0 / 2712554.0);
    const auto fs_225_143_3 = std::sqrt(151875.0 / 20449.0);
    const auto fs_225_22_3 = std::sqrt(151875.0 / 484.0);
    const auto fs_238_2717_3 = std::sqrt(169932.0 / 7382089.0);
    const auto fs_23_1122_6 = std::sqrt(529.0 / 209814.0);
    const auto fs_23_2717_105 = std::sqrt(55545.0 / 7382089.0);
    const auto fs_23_286_105 = std::sqrt(55545.0 / 81796.0);
    const auto fs_23_494_6 = std::sqrt(1587.0 / 122018.0);
    const auto fs_23_52_6 = std::sqrt(1587.0 / 1352.0);
    const auto fs_23_66_6 = std::sqrt(529.0 / 726.0);
    const auto fs_24_143_55 = std::sqrt(2880.0 / 1859.0);
    const auto fs_24_2717_14 = std::sqrt(8064.0 / 7382089.0);
    const auto fs_25_2_7 = std::sqrt(4375.0 / 4.0);
    const auto fs_25_429_7 = std::sqrt(4375.0 / 184041.0);
    const auto fs_25_44_210 = std::sqrt(65625.0 / 968.0);
    const auto fs_270_143_5 = std::sqrt(364500.0 / 20449.0);
    const auto fs_278_2717_2 = std::sqrt(154568.0 / 7382089.0);
    const auto fs_27_5434_858 = std::sqrt(2187.0 / 103246.0);
    const auto fs_27_572_858 = std::sqrt(2187.0 / 1144.0);
    const auto fs_28_2717_105 = std::sqrt(82320.0 / 7382089.0);
    const auto fs_28_2717_143 = std::sqrt(784.0 / 51623.0);
    const auto fs_28_33_5 = std::sqrt(3920.0 / 1089.0);
    const auto fs_28_4199_1430 = std::sqrt(86240.0 / 1356277.0);
    const auto fs_28_4199_3 = std::sqrt(2352.0 / 17631601.0);
    const auto fs_28_4199_385 = std::sqrt(301840.0 / 17631601.0);
    const auto fs_28_4199_546 = std::sqrt(32928.0 / 1356277.0);
    const auto fs_28_561_5 = std::sqrt(3920.0 / 314721.0);
    const auto fs_294_2717_3 = std::sqrt(259308.0 / 7382089.0);
    const auto fs_294_4199_22 = std::sqrt(1901592.0 / 17631601.0);
    const auto fs_29_2717_21 = std::sqrt(17661.0 / 7382089.0);
    const auto fs_29_286_21 = std::sqrt(17661.0 / 81796.0);
    const auto fs_29_5434_154 = std::sqrt(5887.0 / 1342198.0);
    const auto fs_29_572_154 = std::sqrt(5887.0 / 14872.0);
    const auto fs_2_143_15015 = std::sqrt(420.0 / 143.0);
    const auto fs_2_143_30 = std::sqrt(120.0 / 20449.0);
    const auto fs_2_143_66 = std::sqrt(24.0 / 1859.0);
    const auto fs_2_247_195 = std::sqrt(60.0 / 4693.0);
    const auto fs_2_247_546 = std::sqrt(168.0 / 4693.0);
    const auto fs_2_2717_231 = std::sqrt(84.0 / 671099.0);
    const auto fs_2_2717_30030 = std::sqrt(840.0 / 51623.0);
    const auto fs_2_33_210 = std::sqrt(280.0 / 363.0);
    const auto fs_2_561_210 = std::sqrt(280.0 / 104907.0);
    const auto fs_30_1_7 = std::sqrt(6300.0);
    const auto fs_315_16_14 = std::sqrt(694575.0 / 128.0);
    const auto fs_315_16_5 = std::sqrt(496125.0 / 256.0);
    const auto fs_315_16_7 = std::sqrt(694575.0 / 256.0);
    const auto fs_315_32_110 = std::sqrt(5457375.0 / 512.0);
    const auto fs_315_32_2 = std::sqrt(99225.0 / 512.0);
    const auto fs_315_8_3 = std::sqrt(297675.0 / 64.0);
    const auto fs_35_11_5 = std::sqrt(6125.0 / 121.0);
    const auto fs_35_22_35 = std::sqrt(42875.0 / 484.0);
    const auto fs_35_4199_1001 = std::sqrt(94325.0 / 1356277.0);
    const auto fs_35_4199_273 = std::sqrt(25725.0 / 1356277.0);
    const auto fs_35_44_7 = std::sqrt(8575.0 / 1936.0);
    const auto fs_35_88_462 = std::sqrt(25725.0 / 352.0);
    const auto fs_37_2717_66 = std::sqrt(8214.0 / 671099.0);
    const auto fs_37_286_66 = std::sqrt(4107.0 / 3718.0);
    const auto fs_3_11_42 = std::sqrt(378.0 / 121.0);
    const auto fs_3_13_91 = std::sqrt(63.0 / 13.0);
    const auto fs_3_143_11 = std::sqrt(9.0 / 1859.0);
    const auto fs_3_143_210 = std::sqrt(1890.0 / 20449.0);
    const auto fs_3_143_30 = std::sqrt(270.0 / 20449.0);
    const auto fs_3_209_11 = std::sqrt(9.0 / 3971.0);
    const auto fs_3_22_11 = std::sqrt(9.0 / 44.0);
    const auto fs_3_22_154 = std::sqrt(63.0 / 22.0);
    const auto fs_3_22_330 = std::sqrt(135.0 / 22.0);
    const auto fs_3_247_273 = std::sqrt(189.0 / 4693.0);
    const auto fs_3_247_30 = std::sqrt(270.0 / 61009.0);
    const auto fs_3_247_7 = std::sqrt(63.0 / 61009.0);
    const auto fs_3_26_273 = std::sqrt(189.0 / 52.0);
    const auto fs_3_26_30 = std::sqrt(135.0 / 338.0);
    const auto fs_3_26_7 = std::sqrt(63.0 / 676.0);
    const auto fs_3_286_42 = std::sqrt(189.0 / 40898.0);
    const auto fs_3_374_154 = std::sqrt(63.0 / 6358.0);
    const auto fs_3_374_330 = std::sqrt(135.0 / 6358.0);
    const auto fs_40_33_7 = std::sqrt(11200.0 / 1089.0);
    const auto fs_42_143_5 = std::sqrt(8820.0 / 20449.0);
    const auto fs_42_143_55 = std::sqrt(8820.0 / 1859.0);
    const auto fs_42_143_6 = std::sqrt(10584.0 / 20449.0);
    const auto fs_42_2717_143 = std::sqrt(1764.0 / 51623.0);
    const auto fs_42_2717_77 = std::sqrt(12348.0 / 671099.0);
    const auto fs_42_4199_1155 = std::sqrt(2037420.0 / 17631601.0);
    const auto fs_42_4199_231 = std::sqrt(407484.0 / 17631601.0);
    const auto fs_42_4199_286 = std::sqrt(38808.0 / 1356277.0);
    const auto fs_43_2717_33 = std::sqrt(5547.0 / 671099.0);
    const auto fs_43_286_33 = std::sqrt(5547.0 / 7436.0);
    const auto fs_45_11_15 = std::sqrt(30375.0 / 121.0);
    const auto fs_45_143_30 = std::sqrt(60750.0 / 20449.0);
    const auto fs_45_143_66 = std::sqrt(12150.0 / 1859.0);
    const auto fs_45_16_42 = std::sqrt(42525.0 / 128.0);
    const auto fs_45_1_3 = std::sqrt(6075.0);
    const auto fs_45_22_30 = std::sqrt(30375.0 / 242.0);
    const auto fs_45_22_66 = std::sqrt(6075.0 / 22.0);
    const auto fs_45_286_105 = std::sqrt(212625.0 / 81796.0);
    const auto fs_45_286_210 = std::sqrt(212625.0 / 40898.0);
    const auto fs_45_286_42 = std::sqrt(42525.0 / 40898.0);
    const auto fs_45_286_462 = std::sqrt(42525.0 / 3718.0);
    const auto fs_45_2_14 = std::sqrt(14175.0 / 2.0);
    const auto fs_45_2_5 = std::sqrt(10125.0 / 4.0);
    const auto fs_45_2_7 = std::sqrt(14175.0 / 4.0);
    const auto fs_45_44_105 = std::sqrt(212625.0 / 1936.0);
    const auto fs_45_44_210 = std::sqrt(212625.0 / 968.0);
    const auto fs_45_44_42 = std::sqrt(42525.0 / 968.0);
    const auto fs_45_44_462 = std::sqrt(42525.0 / 88.0);
    const auto fs_45_44_6 = std::sqrt(6075.0 / 968.0);
    const auto fs_45_4_11 = std::sqrt(22275.0 / 16.0);
    const auto fs_45_4_110 = std::sqrt(111375.0 / 8.0);
    const auto fs_45_4_2 = std::sqrt(2025.0 / 8.0);
    const auto fs_45_4_7 = std::sqrt(14175.0 / 16.0);
    const auto fs_45_88_154 = std::sqrt(14175.0 / 352.0);
    const auto fs_45_88_330 = std::sqrt(30375.0 / 352.0);
    const auto fs_45_8_11 = std::sqrt(22275.0 / 64.0);
    const auto fs_45_8_30 = std::sqrt(30375.0 / 32.0);
    const auto fs_47_143_21 = std::sqrt(46389.0 / 20449.0);
    const auto fs_48_2717_55 = std::sqrt(11520.0 / 671099.0);
    const auto fs_49_5434_286 = std::sqrt(2401.0 / 103246.0);
    const auto fs_49_572_286 = std::sqrt(2401.0 / 1144.0);
    const auto fs_4_11_2 = std::sqrt(32.0 / 121.0);
    const auto fs_4_11_33 = std::sqrt(48.0 / 11.0);
    const auto fs_4_143_15 = std::sqrt(240.0 / 20449.0);
    const auto fs_4_187_2 = std::sqrt(32.0 / 34969.0);
    const auto fs_4_187_33 = std::sqrt(48.0 / 3179.0);
    const auto fs_4_2717_15015 = std::sqrt(1680.0 / 51623.0);
    const auto fs_50_33_7 = std::sqrt(17500.0 / 1089.0);
    const auto fs_525_16_7 = std::sqrt(1929375.0 / 256.0);
    const auto fs_585_88_6 = std::sqrt(1026675.0 / 3872.0);
    const auto fs_58_143_11 = std::sqrt(3364.0 / 1859.0);
    const auto fs_5_11_110 = std::sqrt(250.0 / 11.0);
    const auto fs_5_11_2 = std::sqrt(50.0 / 121.0);
    const auto fs_5_11_3 = std::sqrt(75.0 / 121.0);
    const auto fs_5_143_14 = std::sqrt(350.0 / 20449.0);
    const auto fs_5_143_5 = std::sqrt(125.0 / 20449.0);
    const auto fs_5_143_7 = std::sqrt(175.0 / 20449.0);
    const auto fs_5_1_30 = std::sqrt(750.0);
    const auto fs_5_1_42 = std::sqrt(1050.0);
    const auto fs_5_22_210 = std::sqrt(2625.0 / 242.0);
    const auto fs_5_2717_2310 = std::sqrt(5250.0 / 671099.0);
    const auto fs_5_286_110 = std::sqrt(125.0 / 3718.0);
    const auto fs_5_286_2 = std::sqrt(25.0 / 40898.0);
    const auto fs_5_286_2310 = std::sqrt(2625.0 / 3718.0);
    const auto fs_5_2_105 = std::sqrt(2625.0 / 4.0);
    const auto fs_5_2_15 = std::sqrt(375.0 / 4.0);
    const auto fs_5_2_165 = std::sqrt(4125.0 / 4.0);
    const auto fs_5_2_210 = std::sqrt(2625.0 / 2.0);
    const auto fs_5_2_35 = std::sqrt(875.0 / 4.0);
    const auto fs_5_2_70 = std::sqrt(875.0 / 2.0);
    const auto fs_5_33_10 = std::sqrt(250.0 / 1089.0);
    const auto fs_5_33_2 = std::sqrt(50.0 / 1089.0);
    const auto fs_5_33_210 = std::sqrt(1750.0 / 363.0);
    const auto fs_5_33_330 = std::sqrt(250.0 / 33.0);
    const auto fs_5_33_70 = std::sqrt(1750.0 / 1089.0);
    const auto fs_5_429_105 = std::sqrt(875.0 / 61347.0);
    const auto fs_5_429_15 = std::sqrt(125.0 / 61347.0);
    const auto fs_5_429_165 = std::sqrt(125.0 / 5577.0);
    const auto fs_5_429_210 = std::sqrt(1750.0 / 61347.0);
    const auto fs_5_429_35 = std::sqrt(875.0 / 184041.0);
    const auto fs_5_429_70 = std::sqrt(1750.0 / 184041.0);
    const auto fs_5_44_2310 = std::sqrt(2625.0 / 88.0);
    const auto fs_5_4_10 = std::sqrt(125.0 / 8.0);
    const auto fs_5_4_2 = std::sqrt(25.0 / 8.0);
    const auto fs_5_4_330 = std::sqrt(4125.0 / 8.0);
    const auto fs_5_4_42 = std::sqrt(525.0 / 8.0);
    const auto fs_5_4_70 = std::sqrt(875.0 / 8.0);
    const auto fs_5_561_210 = std::sqrt(1750.0 / 104907.0);
    const auto fs_5_858_10 = std::sqrt(125.0 / 368082.0);
    const auto fs_5_858_2 = std::sqrt(25.0 / 368082.0);
    const auto fs_5_858_330 = std::sqrt(125.0 / 11154.0);
    const auto fs_5_858_70 = std::sqrt(875.0 / 368082.0);
    const auto fs_63_2717_66 = std::sqrt(23814.0 / 671099.0);
    const auto fs_63_286_66 = std::sqrt(11907.0 / 3718.0);
    const auto fs_6_143_11 = std::sqrt(36.0 / 1859.0);
    const auto fs_6_143_7 = std::sqrt(252.0 / 20449.0);
    const auto fs_6_209_42 = std::sqrt(1512.0 / 43681.0);
    const auto fs_6_247_91 = std::sqrt(252.0 / 4693.0);
    const auto fs_6_2717_210 = std::sqrt(7560.0 / 7382089.0);
    const auto fs_75_2_7 = std::sqrt(39375.0 / 4.0);
    const auto fs_75_4_3 = std::sqrt(16875.0 / 16.0);
    const auto fs_7_1122_462 = std::sqrt(343.0 / 19074.0);
    const auto fs_7_11_14 = std::sqrt(686.0 / 121.0);
    const auto fs_7_11_3 = std::sqrt(147.0 / 121.0);
    const auto fs_7_13_13 = std::sqrt(49.0 / 13.0);
    const auto fs_7_187_14 = std::sqrt(686.0 / 34969.0);
    const auto fs_7_187_3 = std::sqrt(147.0 / 34969.0);
    const auto fs_7_209_21 = std::sqrt(1029.0 / 43681.0);
    const auto fs_7_22_21 = std::sqrt(1029.0 / 484.0);
    const auto fs_7_2717_110 = std::sqrt(490.0 / 671099.0);
    const auto fs_7_2717_4290 = std::sqrt(1470.0 / 51623.0);
    const auto fs_7_286_110 = std::sqrt(245.0 / 3718.0);
    const auto fs_7_286_4290 = std::sqrt(735.0 / 286.0);
    const auto fs_7_33_7 = std::sqrt(343.0 / 1089.0);
    const auto fs_7_4199_12155 = std::sqrt(2695.0 / 79781.0);
    const auto fs_7_4199_15015 = std::sqrt(56595.0 / 1356277.0);
    const auto fs_7_4199_2310 = std::sqrt(113190.0 / 17631601.0);
    const auto fs_7_4199_25194 = std::sqrt(294.0 / 4199.0);
    const auto fs_7_4199_273 = std::sqrt(1029.0 / 1356277.0);
    const auto fs_7_4199_3003 = std::sqrt(11319.0 / 1356277.0);
    const auto fs_7_4199_3094 = std::sqrt(686.0 / 79781.0);
    const auto fs_7_4199_62985 = std::sqrt(735.0 / 4199.0);
    const auto fs_7_4199_9282 = std::sqrt(2058.0 / 79781.0);
    const auto fs_7_561_7 = std::sqrt(343.0 / 314721.0);
    const auto fs_7_66_462 = std::sqrt(343.0 / 66.0);
    const auto fs_7_8398_182 = std::sqrt(343.0 / 2712554.0);
    const auto fs_7_8398_2 = std::sqrt(49.0 / 35263202.0);
    const auto fs_7_8398_26 = std::sqrt(49.0 / 2712554.0);
    const auto fs_7_8398_330 = std::sqrt(8085.0 / 35263202.0);
    const auto fs_7_8398_910 = std::sqrt(1715.0 / 2712554.0);
    const auto fs_84_2717_5 = std::sqrt(35280.0 / 7382089.0);
    const auto fs_84_2717_55 = std::sqrt(35280.0 / 671099.0);
    const auto fs_84_2717_6 = std::sqrt(42336.0 / 7382089.0);
    const auto fs_84_4199_154 = std::sqrt(1086624.0 / 17631601.0);
    const auto fs_84_4199_210 = std::sqrt(1481760.0 / 17631601.0);
    const auto fs_84_4199_22 = std::sqrt(155232.0 / 17631601.0);
    const auto fs_84_4199_221 = std::sqrt(7056.0 / 79781.0);
    const auto fs_90_143_15 = std::sqrt(121500.0 / 20449.0);
    const auto fs_94_2717_21 = std::sqrt(185556.0 / 7382089.0);

    // NOTE: the angular components are formed in 87 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_p2, ph4_p2, ph4_p3, ph6_p2, ph6_p3, ph8_p2, ph8_p3, ph10_p2, ph10_p3, ph10_p9, ph10_p10, ab_2, pc_0, pc_1 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

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

        pc_0[k] = e_0 * fs_315_32_110 * h2_p2 + e_1 * fs_15_4_66 * h4_p2 - e_1 * fs_45_4_110 * r_2 * h2_p2 + e_2 * fs_15_44_77 * h6_p2 - e_2 * fs_45_22_66 * r_2 * h4_p2 + e_2 * fs_15_4_110 * r_4 * h2_p2 + e_3 * fs_1_143_231 * h8_p2 - e_3 * fs_1_11_77 * r_2 * h6_p2 + e_3 * fs_45_143_66 * r_4 * h4_p2 - e_3 * fs_5_11_110 * r_6 * h2_p2 + e_4 * fs_7_8398_2 * h10_p2 - e_4 * fs_7_4199_62985 * h10_p10 - e_4 * fs_2_2717_231 * r_2 * h8_p2 + e_4 * fs_1_187_77 * r_4 * h6_p2 - e_4 * fs_2_143_66 * r_6 * h4_p2 + e_4 * fs_5_286_110 * r_8 * h2_p2;

        pc_1[k] = - e_1 * fs_15_8_462 * h4_p3 - e_2 * fs_45_88_154 * h6_p3 + e_2 * fs_45_44_462 * r_2 * h4_p3 - e_3 * fs_3_26_7 * h8_p3 + e_3 * fs_3_22_154 * r_2 * h6_p3 - e_3 * fs_45_286_462 * r_4 * h4_p3 - e_4 * fs_7_8398_26 * h10_p3 - e_4 * fs_7_4199_25194 * h10_p9 + e_4 * fs_3_247_7 * r_2 * h8_p3 - e_4 * fs_3_374_154 * r_4 * h6_p3 + e_4 * fs_1_143_462 * r_6 * h4_p3;
    }

    // NOTE: the angular components are formed in 87 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, ph4_p4, ph6_p4, ph6_p5, ph8_p4, ph8_p5, ph8_p7, ph8_p8, ph10_p4, ph10_p5, ph10_p7, ph10_p8, ab_2, pc_2, pc_3 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

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

        pc_2[k] = e_1 * fs_15_4_66 * h4_p4 + e_2 * fs_45_88_330 * h6_p4 - e_2 * fs_45_22_66 * r_2 * h4_p4 + e_3 * fs_3_26_30 * h8_p4 - e_3 * fs_1_13_546 * h8_p8 - e_3 * fs_3_22_330 * r_2 * h6_p4 + e_3 * fs_45_143_66 * r_4 * h4_p4 + e_4 * fs_7_8398_182 * h10_p4 - e_4 * fs_7_4199_9282 * h10_p8 - e_4 * fs_3_247_30 * r_2 * h8_p4 + e_4 * fs_2_247_546 * r_2 * h8_p8 + e_4 * fs_3_374_330 * r_4 * h6_p4 - e_4 * fs_2_143_66 * r_6 * h4_p4;

        pc_3[k] = - e_2 * fs_15_8_30 * h6_p5 - e_3 * fs_1_13_195 * h8_p5 - e_3 * fs_3_26_273 * h8_p7 + e_3 * fs_1_2_30 * r_2 * h6_p5 - e_4 * fs_7_8398_910 * h10_p5 - e_4 * fs_7_4199_3094 * h10_p7 + e_4 * fs_2_247_195 * r_2 * h8_p5 + e_4 * fs_3_247_273 * r_2 * h8_p7 - e_4 * fs_1_34_30 * r_4 * h6_p5;
    }

    // NOTE: the angular components are formed in 87 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_2, pe_3, pe_4, ph6_m6, ph6_m5, ph8_m7, ph8_m6, ph8_m5, ph10_m7, ph10_m6, ph10_m5, ab_2, pc_4, pc_5 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

        const auto h6_m6 = ph6_m6[k];
        const auto h6_m5 = ph6_m5[k];
        const auto h8_m7 = ph8_m7[k];
        const auto h8_m6 = ph8_m6[k];
        const auto h8_m5 = ph8_m5[k];
        const auto h10_m7 = ph10_m7[k];
        const auto h10_m6 = ph10_m6[k];
        const auto h10_m5 = ph10_m5[k];

        pc_4[k] = e_2 * f_45_4 * h6_m6 + e_3 * fs_3_13_91 * h8_m6 - e_3 * f_3_1 * r_2 * h6_m6 + e_4 * fs_14_4199_455 * h10_m6 - e_4 * fs_6_247_91 * r_2 * h8_m6 + e_4 * f_3_17 * r_4 * h6_m6;

        pc_5[k] = - e_2 * fs_15_8_30 * h6_m5 + e_3 * fs_3_26_273 * h8_m7 - e_3 * fs_1_13_195 * h8_m5 + e_3 * fs_1_2_30 * r_2 * h6_m5 + e_4 * fs_7_4199_3094 * h10_m7 - e_4 * fs_7_8398_910 * h10_m5 - e_4 * fs_3_247_273 * r_2 * h8_m7 + e_4 * fs_2_247_195 * r_2 * h8_m5 - e_4 * fs_1_34_30 * r_4 * h6_m5;
    }

    // NOTE: the angular components are formed in 87 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, ph4_m4, ph4_m3, ph6_m4, ph6_m3, ph8_m8, ph8_m4, ph8_m3, ph10_m9, ph10_m8, ph10_m4, ph10_m3, ab_2, pc_6, pc_7 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

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

        pc_6[k] = e_1 * fs_15_4_66 * h4_m4 + e_2 * fs_45_88_330 * h6_m4 - e_2 * fs_45_22_66 * r_2 * h4_m4 + e_3 * fs_1_13_546 * h8_m8 + e_3 * fs_3_26_30 * h8_m4 - e_3 * fs_3_22_330 * r_2 * h6_m4 + e_3 * fs_45_143_66 * r_4 * h4_m4 + e_4 * fs_7_4199_9282 * h10_m8 + e_4 * fs_7_8398_182 * h10_m4 - e_4 * fs_2_247_546 * r_2 * h8_m8 - e_4 * fs_3_247_30 * r_2 * h8_m4 + e_4 * fs_3_374_330 * r_4 * h6_m4 - e_4 * fs_2_143_66 * r_6 * h4_m4;

        pc_7[k] = - e_1 * fs_15_8_462 * h4_m3 - e_2 * fs_45_88_154 * h6_m3 + e_2 * fs_45_44_462 * r_2 * h4_m3 - e_3 * fs_3_26_7 * h8_m3 + e_3 * fs_3_22_154 * r_2 * h6_m3 - e_3 * fs_45_286_462 * r_4 * h4_m3 + e_4 * fs_7_4199_25194 * h10_m9 - e_4 * fs_7_8398_26 * h10_m3 + e_4 * fs_3_247_7 * r_2 * h8_m3 - e_4 * fs_3_374_154 * r_4 * h6_m3 + e_4 * fs_1_143_462 * r_6 * h4_m3;
    }

    // NOTE: the angular components are formed in 87 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_m2, ph2_p1, ph4_m2, ph4_p1, ph6_m2, ph6_p1, ph8_m2, ph8_p1, ph10_m10, ph10_m2, ph10_p1, ph10_p9, ab_2, pc_8, pc_9 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

        const auto h2_m2 = ph2_m2[k];
        const auto h2_p1 = ph2_p1[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h8_m2 = ph8_m2[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h10_m10 = ph10_m10[k];
        const auto h10_m2 = ph10_m2[k];
        const auto h10_p1 = ph10_p1[k];
        const auto h10_p9 = ph10_p9[k];

        pc_8[k] = e_0 * fs_315_32_110 * h2_m2 + e_1 * fs_15_4_66 * h4_m2 - e_1 * fs_45_4_110 * r_2 * h2_m2 + e_2 * fs_15_44_77 * h6_m2 - e_2 * fs_45_22_66 * r_2 * h4_m2 + e_2 * fs_15_4_110 * r_4 * h2_m2 + e_3 * fs_1_143_231 * h8_m2 - e_3 * fs_1_11_77 * r_2 * h6_m2 + e_3 * fs_45_143_66 * r_4 * h4_m2 - e_3 * fs_5_11_110 * r_6 * h2_m2 + e_4 * fs_7_4199_62985 * h10_m10 + e_4 * fs_7_8398_2 * h10_m2 - e_4 * fs_2_2717_231 * r_2 * h8_m2 + e_4 * fs_1_187_77 * r_4 * h6_m2 - e_4 * fs_2_143_66 * r_6 * h4_m2 + e_4 * fs_5_286_110 * r_8 * h2_m2;

        pc_9[k] = e_0 * fs_105_32_330 * h2_p1 + e_1 * fs_45_4_11 * h4_p1 - e_1 * fs_15_4_330 * r_2 * h2_p1 + e_2 * fs_5_44_2310 * h6_p1 - e_2 * fs_135_22_11 * r_2 * h4_p1 + e_2 * fs_5_4_330 * r_4 * h2_p1 + e_3 * fs_7_286_110 * h8_p1 - e_3 * fs_1_33_2310 * r_2 * h6_p1 + e_3 * fs_135_143_11 * r_4 * h4_p1 - e_3 * fs_5_33_330 * r_6 * h2_p1 + e_4 * fs_21_8398_2 * h10_p1 - e_4 * fs_21_4199_4199 * h10_p9 - e_4 * fs_7_2717_110 * r_2 * h8_p1 + e_4 * fs_1_561_2310 * r_4 * h6_p1 - e_4 * fs_6_143_11 * r_6 * h4_p1 + e_4 * fs_5_858_330 * r_8 * h2_p1;
    }

    // NOTE: the angular components are formed in 87 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_p2, ph4_p2, ph6_p2, ph8_p2, ph8_p8, ph10_p2, ph10_p8, ab_2, pc_10 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

        const auto h2_p2 = ph2_p2[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h8_p2 = ph8_p2[k];
        const auto h8_p8 = ph8_p8[k];
        const auto h10_p2 = ph10_p2[k];
        const auto h10_p8 = ph10_p8[k];

        pc_10[k] = e_0 * fs_105_16_165 * h2_p2 - e_1 * fs_45_8_11 * h4_p2 - e_1 * fs_15_2_165 * r_2 * h2_p2 - e_2 * fs_35_88_462 * h6_p2 + e_2 * fs_135_44_11 * r_2 * h4_p2 + e_2 * fs_5_2_165 * r_4 * h2_p2 - e_3 * fs_29_572_154 * h8_p2 + e_3 * fs_7_13_13 * h8_p8 + e_3 * fs_7_66_462 * r_2 * h6_p2 - e_3 * fs_135_286_11 * r_4 * h4_p2 - e_3 * fs_10_33_165 * r_6 * h2_p2 - e_4 * fs_28_4199_3 * h10_p2 - e_4 * fs_84_4199_221 * h10_p8 + e_4 * fs_29_5434_154 * r_2 * h8_p2 - e_4 * fs_14_247_13 * r_2 * h8_p8 - e_4 * fs_7_1122_462 * r_4 * h6_p2 + e_4 * fs_3_143_11 * r_6 * h4_p2 + e_4 * fs_5_429_165 * r_8 * h2_p2;
    }

    // NOTE: the angular components are formed in 87 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, ph4_p3, ph4_p4, ph6_p3, ph6_p4, ph6_p6, ph8_p3, ph8_p4, ph8_p6, ph8_p7, ph10_p3, ph10_p4, ph10_p6, ph10_p7, ab_2, pc_11, pc_12 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

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

        pc_11[k] = - e_1 * fs_45_8_11 * h4_p3 + e_2 * fs_15_11_33 * h6_p3 + e_2 * fs_135_44_11 * r_2 * h4_p3 + e_3 * fs_23_52_6 * h8_p3 + e_3 * fs_1_52_182 * h8_p7 - e_3 * fs_4_11_33 * r_2 * h6_p3 - e_3 * fs_135_286_11 * r_4 * h4_p3 + e_4 * fs_7_4199_273 * h10_p3 - e_4 * fs_14_4199_4641 * h10_p7 - e_4 * fs_23_494_6 * r_2 * h8_p3 - e_4 * fs_1_494_182 * r_2 * h8_p7 + e_4 * fs_4_187_33 * r_4 * h6_p3 + e_4 * fs_3_143_11 * r_6 * h4_p3;

        pc_12[k] = e_1 * fs_45_4_11 * h4_p4 - e_2 * fs_15_44_55 * h6_p4 + e_2 * fs_15_8_30 * h6_p6 - e_2 * fs_135_22_11 * r_2 * h4_p4 - e_3 * fs_17_26_5 * h8_p4 - e_3 * fs_1_52_2730 * h8_p6 + e_3 * fs_1_11_55 * r_2 * h6_p4 - e_3 * fs_1_2_30 * r_2 * h6_p6 + e_3 * fs_135_143_11 * r_4 * h4_p4 - e_4 * fs_14_4199_273 * h10_p4 - e_4 * fs_28_4199_546 * h10_p6 + e_4 * fs_17_247_5 * r_2 * h8_p4 + e_4 * fs_1_494_2730 * r_2 * h8_p6 - e_4 * fs_1_187_55 * r_4 * h6_p4 + e_4 * fs_1_34_30 * r_4 * h6_p6 - e_4 * fs_6_143_11 * r_6 * h4_p4;
    }

    // NOTE: the angular components are formed in 87 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, ph4_m4, ph6_m6, ph6_m5, ph6_m4, ph8_m6, ph8_m5, ph8_m4, ph10_m6, ph10_m5, ph10_m4, ab_2, pc_13, pc_14 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

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

        pc_13[k] = - e_2 * f_15_2 * h6_m5 + e_3 * fs_11_26_26 * h8_m5 + e_3 * f_2_1 * r_2 * h6_m5 + e_4 * fs_35_4199_273 * h10_m5 - e_4 * fs_11_247_26 * r_2 * h8_m5 - e_4 * f_2_17 * r_4 * h6_m5;

        pc_14[k] = e_1 * fs_45_4_11 * h4_m4 - e_2 * fs_15_8_30 * h6_m6 - e_2 * fs_15_44_55 * h6_m4 - e_2 * fs_135_22_11 * r_2 * h4_m4 + e_3 * fs_1_52_2730 * h8_m6 - e_3 * fs_17_26_5 * h8_m4 + e_3 * fs_1_2_30 * r_2 * h6_m6 + e_3 * fs_1_11_55 * r_2 * h6_m4 + e_3 * fs_135_143_11 * r_4 * h4_m4 + e_4 * fs_28_4199_546 * h10_m6 - e_4 * fs_14_4199_273 * h10_m4 - e_4 * fs_1_494_2730 * r_2 * h8_m6 + e_4 * fs_17_247_5 * r_2 * h8_m4 - e_4 * fs_1_34_30 * r_4 * h6_m6 - e_4 * fs_1_187_55 * r_4 * h6_m4 - e_4 * fs_6_143_11 * r_6 * h4_m4;
    }

    // NOTE: the angular components are formed in 87 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, ph4_m3, ph6_m3, ph8_m7, ph8_m3, ph10_m7, ph10_m3, ab_2, pc_15 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h4_m3 = ph4_m3[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h8_m7 = ph8_m7[k];
        const auto h8_m3 = ph8_m3[k];
        const auto h10_m7 = ph10_m7[k];
        const auto h10_m3 = ph10_m3[k];

        pc_15[k] = - e_1 * fs_45_8_11 * h4_m3 + e_2 * fs_15_11_33 * h6_m3 + e_2 * fs_135_44_11 * r_2 * h4_m3 - e_3 * fs_1_52_182 * h8_m7 + e_3 * fs_23_52_6 * h8_m3 - e_3 * fs_4_11_33 * r_2 * h6_m3 - e_3 * fs_135_286_11 * r_4 * h4_m3 + e_4 * fs_14_4199_4641 * h10_m7 + e_4 * fs_7_4199_273 * h10_m3 + e_4 * fs_1_494_182 * r_2 * h8_m7 - e_4 * fs_23_494_6 * r_2 * h8_m3 + e_4 * fs_4_187_33 * r_4 * h6_m3 + e_4 * fs_3_143_11 * r_6 * h4_m3;
    }

    // NOTE: the angular components are formed in 87 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_m2, ph4_m2, ph6_m2, ph8_m8, ph8_m2, ph10_m8, ph10_m2, ab_2, pc_16 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

        const auto h2_m2 = ph2_m2[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h8_m8 = ph8_m8[k];
        const auto h8_m2 = ph8_m2[k];
        const auto h10_m8 = ph10_m8[k];
        const auto h10_m2 = ph10_m2[k];

        pc_16[k] = e_0 * fs_105_16_165 * h2_m2 - e_1 * fs_45_8_11 * h4_m2 - e_1 * fs_15_2_165 * r_2 * h2_m2 - e_2 * fs_35_88_462 * h6_m2 + e_2 * fs_135_44_11 * r_2 * h4_m2 + e_2 * fs_5_2_165 * r_4 * h2_m2 - e_3 * fs_7_13_13 * h8_m8 - e_3 * fs_29_572_154 * h8_m2 + e_3 * fs_7_66_462 * r_2 * h6_m2 - e_3 * fs_135_286_11 * r_4 * h4_m2 - e_3 * fs_10_33_165 * r_6 * h2_m2 + e_4 * fs_84_4199_221 * h10_m8 - e_4 * fs_28_4199_3 * h10_m2 + e_4 * fs_14_247_13 * r_2 * h8_m8 + e_4 * fs_29_5434_154 * r_2 * h8_m2 - e_4 * fs_7_1122_462 * r_4 * h6_m2 + e_4 * fs_3_143_11 * r_6 * h4_m2 + e_4 * fs_5_429_165 * r_8 * h2_m2;
    }

    // NOTE: the angular components are formed in 87 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_m1, ph4_m1, ph6_m1, ph8_m1, ph10_m9, ph10_m1, ab_2, pc_17 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

        const auto h2_m1 = ph2_m1[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h8_m1 = ph8_m1[k];
        const auto h10_m9 = ph10_m9[k];
        const auto h10_m1 = ph10_m1[k];

        pc_17[k] = e_0 * fs_105_32_330 * h2_m1 + e_1 * fs_45_4_11 * h4_m1 - e_1 * fs_15_4_330 * r_2 * h2_m1 + e_2 * fs_5_44_2310 * h6_m1 - e_2 * fs_135_22_11 * r_2 * h4_m1 + e_2 * fs_5_4_330 * r_4 * h2_m1 + e_3 * fs_7_286_110 * h8_m1 - e_3 * fs_1_33_2310 * r_2 * h6_m1 + e_3 * fs_135_143_11 * r_4 * h4_m1 - e_3 * fs_5_33_330 * r_6 * h2_m1 + e_4 * fs_21_4199_4199 * h10_m9 + e_4 * fs_21_8398_2 * h10_m1 - e_4 * fs_7_2717_110 * r_2 * h8_m1 + e_4 * fs_1_561_2310 * r_4 * h6_m1 - e_4 * fs_6_143_11 * r_6 * h4_m1 + e_4 * fs_5_858_330 * r_8 * h2_m1;
    }

    // NOTE: the angular components are formed in 87 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_0, ph4_0, ph6_0, ph8_0, ph8_p8, ph10_0, ph10_p8, ab_2, pc_18 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

        const auto h2_0 = ph2_0[k];
        const auto h4_0 = ph4_0[k];
        const auto h6_0 = ph6_0[k];
        const auto h8_0 = ph8_0[k];
        const auto h8_p8 = ph8_p8[k];
        const auto h10_0 = ph10_0[k];
        const auto h10_p8 = ph10_p8[k];

        pc_18[k] = e_0 * fs_315_16_5 * h2_0 + e_1 * fs_45_2_5 * h4_0 - e_1 * fs_45_2_5 * r_2 * h2_0 + e_2 * fs_105_22_5 * h6_0 - e_2 * fs_135_11_5 * r_2 * h4_0 + e_2 * fs_15_2_5 * r_4 * h2_0 + e_3 * fs_42_143_5 * h8_0 - e_3 * fs_14_143_143 * h8_p8 - e_3 * fs_14_11_5 * r_2 * h6_0 + e_3 * fs_270_143_5 * r_4 * h4_0 - e_3 * fs_10_11_5 * r_6 * h2_0 + e_4 * fs_21_4199_5 * h10_0 - e_4 * fs_21_4199_2431 * h10_p8 - e_4 * fs_84_2717_5 * r_2 * h8_0 + e_4 * fs_28_2717_143 * r_2 * h8_p8 + e_4 * fs_14_187_5 * r_4 * h6_0 - e_4 * fs_12_143_5 * r_6 * h4_0 + e_4 * fs_5_143_5 * r_8 * h2_0;
    }

    // NOTE: the angular components are formed in 87 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_p1, ph4_p1, ph6_p1, ph8_p1, ph8_p7, ph10_p1, ph10_p7, ab_2, pc_19 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

        const auto h2_p1 = ph2_p1[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h8_p7 = ph8_p7[k];
        const auto h10_p1 = ph10_p1[k];
        const auto h10_p7 = ph10_p7[k];

        pc_19[k] = e_0 * fs_105_8_30 * h2_p1 + e_1 * f_45_8 * h4_p1 - e_1 * fs_15_1_30 * r_2 * h2_p1 - e_2 * fs_25_44_210 * h6_p1 - e_2 * f_135_44 * r_2 * h4_p1 + e_2 * fs_5_1_30 * r_4 * h2_p1 - e_3 * fs_175_572_10 * h8_p1 + e_3 * fs_49_572_286 * h8_p7 + e_3 * fs_5_33_210 * r_2 * h6_p1 + e_3 * f_135_286 * r_4 * h4_p1 - e_3 * fs_20_33_30 * r_6 * h2_p1 - e_4 * fs_21_4199_22 * h10_p1 - e_4 * fs_14_4199_7293 * h10_p7 + e_4 * fs_175_5434_10 * r_2 * h8_p1 - e_4 * fs_49_5434_286 * r_2 * h8_p7 - e_4 * fs_5_561_210 * r_4 * h6_p1 - e_4 * f_3_143 * r_6 * h4_p1 + e_4 * fs_10_429_30 * r_8 * h2_p1;
    }

    // NOTE: the angular components are formed in 87 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_p2, ph4_p2, ph6_p2, ph6_p6, ph8_p2, ph8_p6, ph10_p2, ph10_p6, ab_2, pc_20 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

        const auto h2_p2 = ph2_p2[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h6_p6 = ph6_p6[k];
        const auto h8_p2 = ph8_p2[k];
        const auto h8_p6 = ph8_p6[k];
        const auto h10_p2 = ph10_p2[k];
        const auto h10_p6 = ph10_p6[k];

        pc_20[k] = e_0 * fs_105_16_105 * h2_p2 - e_1 * fs_45_4_7 * h4_p2 - e_1 * fs_15_2_105 * r_2 * h2_p2 + e_2 * fs_115_88_6 * h6_p2 - e_2 * fs_45_88_330 * h6_p6 + e_2 * fs_135_22_7 * r_2 * h4_p2 + e_2 * fs_5_2_105 * r_4 * h2_p2 + e_3 * fs_139_143_2 * h8_p2 + e_3 * fs_1_143_30030 * h8_p6 - e_3 * fs_23_66_6 * r_2 * h6_p2 + e_3 * fs_3_22_330 * r_2 * h6_p6 - e_3 * fs_135_143_7 * r_4 * h4_p2 - e_3 * fs_10_33_105 * r_6 * h2_p2 + e_4 * fs_14_4199_231 * h10_p2 - e_4 * fs_14_4199_6006 * h10_p6 - e_4 * fs_278_2717_2 * r_2 * h8_p2 - e_4 * fs_2_2717_30030 * r_2 * h8_p6 + e_4 * fs_23_1122_6 * r_4 * h6_p2 - e_4 * fs_3_374_330 * r_4 * h6_p6 + e_4 * fs_6_143_7 * r_6 * h4_p2 + e_4 * fs_5_429_105 * r_8 * h2_p2;
    }

    // NOTE: the angular components are formed in 87 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, ph4_m4, ph4_p3, ph6_m4, ph6_p3, ph6_p5, ph8_m4, ph8_p3, ph8_p5, ph10_m4, ph10_p3, ph10_p5, ab_2, pc_21, pc_22 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

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

        pc_21[k] = e_1 * f_45_8 * h4_p3 + e_2 * fs_105_44_3 * h6_p3 + e_2 * fs_15_44_55 * h6_p5 - e_2 * f_135_44 * r_2 * h4_p3 - e_3 * fs_101_572_66 * h8_p3 + e_3 * fs_1_572_1430 * h8_p5 - e_3 * fs_7_11_3 * r_2 * h6_p3 - e_3 * fs_1_11_55 * r_2 * h6_p5 + e_3 * f_135_286 * r_4 * h4_p3 - e_4 * fs_7_4199_3003 * h10_p3 - e_4 * fs_7_4199_15015 * h10_p5 + e_4 * fs_101_5434_66 * r_2 * h8_p3 - e_4 * fs_1_5434_1430 * r_2 * h8_p5 + e_4 * fs_7_187_3 * r_4 * h6_p3 + e_4 * fs_1_187_55 * r_4 * h6_p5 - e_4 * f_3_143 * r_6 * h4_p3;

        pc_22[k] = e_1 * fs_45_2_5 * h4_m4 - e_2 * f_120_11 * h6_m4 - e_2 * fs_135_11_5 * r_2 * h4_m4 + e_3 * fs_58_143_11 * h8_m4 + e_3 * f_32_11 * r_2 * h6_m4 + e_3 * fs_270_143_5 * r_4 * h4_m4 + e_4 * fs_7_4199_15015 * h10_m4 - e_4 * fs_116_2717_11 * r_2 * h8_m4 - e_4 * f_32_187 * r_4 * h6_m4 - e_4 * fs_12_143_5 * r_6 * h4_m4;
    }

    // NOTE: the angular components are formed in 87 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, ph4_m3, ph6_m5, ph6_m3, ph8_m5, ph8_m3, ph10_m5, ph10_m3, ab_2, pc_23 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h4_m3 = ph4_m3[k];
        const auto h6_m5 = ph6_m5[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h8_m5 = ph8_m5[k];
        const auto h8_m3 = ph8_m3[k];
        const auto h10_m5 = ph10_m5[k];
        const auto h10_m3 = ph10_m3[k];

        pc_23[k] = e_1 * f_45_8 * h4_m3 - e_2 * fs_15_44_55 * h6_m5 + e_2 * fs_105_44_3 * h6_m3 - e_2 * f_135_44 * r_2 * h4_m3 - e_3 * fs_1_572_1430 * h8_m5 - e_3 * fs_101_572_66 * h8_m3 + e_3 * fs_1_11_55 * r_2 * h6_m5 - e_3 * fs_7_11_3 * r_2 * h6_m3 + e_3 * f_135_286 * r_4 * h4_m3 + e_4 * fs_7_4199_15015 * h10_m5 - e_4 * fs_7_4199_3003 * h10_m3 + e_4 * fs_1_5434_1430 * r_2 * h8_m5 + e_4 * fs_101_5434_66 * r_2 * h8_m3 - e_4 * fs_1_187_55 * r_4 * h6_m5 + e_4 * fs_7_187_3 * r_4 * h6_m3 - e_4 * f_3_143 * r_6 * h4_m3;
    }

    // NOTE: the angular components are formed in 87 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_m2, ph4_m2, ph6_m6, ph6_m2, ph8_m6, ph8_m2, ph10_m6, ph10_m2, ab_2, pc_24 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

        const auto h2_m2 = ph2_m2[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h6_m6 = ph6_m6[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h8_m6 = ph8_m6[k];
        const auto h8_m2 = ph8_m2[k];
        const auto h10_m6 = ph10_m6[k];
        const auto h10_m2 = ph10_m2[k];

        pc_24[k] = e_0 * fs_105_16_105 * h2_m2 - e_1 * fs_45_4_7 * h4_m2 - e_1 * fs_15_2_105 * r_2 * h2_m2 + e_2 * fs_45_88_330 * h6_m6 + e_2 * fs_115_88_6 * h6_m2 + e_2 * fs_135_22_7 * r_2 * h4_m2 + e_2 * fs_5_2_105 * r_4 * h2_m2 - e_3 * fs_1_143_30030 * h8_m6 + e_3 * fs_139_143_2 * h8_m2 - e_3 * fs_3_22_330 * r_2 * h6_m6 - e_3 * fs_23_66_6 * r_2 * h6_m2 - e_3 * fs_135_143_7 * r_4 * h4_m2 - e_3 * fs_10_33_105 * r_6 * h2_m2 + e_4 * fs_14_4199_6006 * h10_m6 + e_4 * fs_14_4199_231 * h10_m2 + e_4 * fs_2_2717_30030 * r_2 * h8_m6 - e_4 * fs_278_2717_2 * r_2 * h8_m2 + e_4 * fs_3_374_330 * r_4 * h6_m6 + e_4 * fs_23_1122_6 * r_4 * h6_m2 + e_4 * fs_6_143_7 * r_6 * h4_m2 + e_4 * fs_5_429_105 * r_8 * h2_m2;
    }

    // NOTE: the angular components are formed in 87 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_m1, ph4_m1, ph6_m1, ph8_m8, ph8_m7, ph8_m1, ph10_m8, ph10_m7, ph10_m1, ab_2, pc_25, pc_26 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

        const auto h2_m1 = ph2_m1[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h8_m8 = ph8_m8[k];
        const auto h8_m7 = ph8_m7[k];
        const auto h8_m1 = ph8_m1[k];
        const auto h10_m8 = ph10_m8[k];
        const auto h10_m7 = ph10_m7[k];
        const auto h10_m1 = ph10_m1[k];

        pc_25[k] = e_0 * fs_105_8_30 * h2_m1 + e_1 * f_45_8 * h4_m1 - e_1 * fs_15_1_30 * r_2 * h2_m1 - e_2 * fs_25_44_210 * h6_m1 - e_2 * f_135_44 * r_2 * h4_m1 + e_2 * fs_5_1_30 * r_4 * h2_m1 - e_3 * fs_49_572_286 * h8_m7 - e_3 * fs_175_572_10 * h8_m1 + e_3 * fs_5_33_210 * r_2 * h6_m1 + e_3 * f_135_286 * r_4 * h4_m1 - e_3 * fs_20_33_30 * r_6 * h2_m1 + e_4 * fs_14_4199_7293 * h10_m7 - e_4 * fs_21_4199_22 * h10_m1 + e_4 * fs_49_5434_286 * r_2 * h8_m7 + e_4 * fs_175_5434_10 * r_2 * h8_m1 - e_4 * fs_5_561_210 * r_4 * h6_m1 - e_4 * f_3_143 * r_6 * h4_m1 + e_4 * fs_10_429_30 * r_8 * h2_m1;

        pc_26[k] = e_3 * fs_14_143_143 * h8_m8 + e_4 * fs_21_4199_2431 * h10_m8 - e_4 * fs_28_2717_143 * r_2 * h8_m8;
    }

    // NOTE: the angular components are formed in 87 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_p1, ph4_p1, ph6_p1, ph8_p1, ph8_p7, ph10_p1, ph10_p7, ab_2, pc_27 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

        const auto h2_p1 = ph2_p1[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h8_p7 = ph8_p7[k];
        const auto h10_p1 = ph10_p1[k];
        const auto h10_p7 = ph10_p7[k];

        pc_27[k] = - e_0 * fs_315_32_2 * h2_p1 - e_1 * fs_15_2_15 * h4_p1 + e_1 * fs_45_4_2 * r_2 * h2_p1 - e_2 * fs_105_44_14 * h6_p1 + e_2 * fs_45_11_15 * r_2 * h4_p1 - e_2 * fs_15_4_2 * r_4 * h2_p1 - e_3 * fs_42_143_6 * h8_p1 - e_3 * fs_7_286_4290 * h8_p7 + e_3 * fs_7_11_14 * r_2 * h6_p1 - e_3 * fs_90_143_15 * r_4 * h4_p1 + e_3 * fs_5_11_2 * r_6 * h2_p1 - e_4 * fs_7_8398_330 * h10_p1 - e_4 * fs_7_4199_12155 * h10_p7 + e_4 * fs_84_2717_6 * r_2 * h8_p1 + e_4 * fs_7_2717_4290 * r_2 * h8_p7 - e_4 * fs_7_187_14 * r_4 * h6_p1 + e_4 * fs_4_143_15 * r_6 * h4_p1 - e_4 * fs_5_286_2 * r_8 * h2_p1;
    }

    // NOTE: the angular components are formed in 87 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_0, ph4_0, ph6_0, ph6_p6, ph8_0, ph8_p6, ph10_0, ph10_p6, ab_2, pc_28 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

        const auto h2_0 = ph2_0[k];
        const auto h4_0 = ph4_0[k];
        const auto h6_0 = ph6_0[k];
        const auto h6_p6 = ph6_p6[k];
        const auto h8_0 = ph8_0[k];
        const auto h8_p6 = ph8_p6[k];
        const auto h10_0 = ph10_0[k];
        const auto h10_p6 = ph10_p6[k];

        pc_28[k] = e_0 * fs_315_8_3 * h2_0 + e_1 * fs_75_4_3 * h4_0 - e_1 * fs_45_1_3 * r_2 * h2_0 - e_2 * fs_105_22_3 * h6_0 + e_2 * fs_45_88_154 * h6_p6 - e_2 * fs_225_22_3 * r_2 * h4_0 + e_2 * fs_15_1_3 * r_4 * h2_0 - e_3 * fs_147_143_3 * h8_0 + e_3 * fs_21_572_286 * h8_p6 + e_3 * fs_14_11_3 * r_2 * h6_0 - e_3 * fs_3_22_154 * r_2 * h6_p6 + e_3 * fs_225_143_3 * r_4 * h4_0 - e_3 * fs_20_11_3 * r_6 * h2_0 - e_4 * fs_140_4199_3 * h10_0 - e_4 * fs_28_4199_1430 * h10_p6 + e_4 * fs_294_2717_3 * r_2 * h8_0 - e_4 * fs_21_5434_286 * r_2 * h8_p6 - e_4 * fs_14_187_3 * r_4 * h6_0 + e_4 * fs_3_374_154 * r_4 * h6_p6 - e_4 * fs_10_143_3 * r_6 * h4_0 + e_4 * fs_10_143_3 * r_8 * h2_0;
    }

    // NOTE: the angular components are formed in 87 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_p1, ph4_p1, ph6_p1, ph6_p5, ph8_p1, ph8_p5, ph10_p1, ph10_p5, ab_2, pc_29 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

        const auto h2_p1 = ph2_p1[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h6_p5 = ph6_p5[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h8_p5 = ph8_p5[k];
        const auto h10_p1 = ph10_p1[k];
        const auto h10_p5 = ph10_p5[k];

        pc_29[k] = e_0 * fs_315_16_14 * h2_p1 - e_1 * fs_15_8_105 * h4_p1 - e_1 * fs_45_2_14 * r_2 * h2_p1 - e_2 * fs_15_11_2 * h6_p1 - e_2 * fs_15_11_33 * h6_p5 + e_2 * fs_45_44_105 * r_2 * h4_p1 + e_2 * fs_15_2_14 * r_4 * h2_p1 + e_3 * fs_127_572_42 * h8_p1 + e_3 * fs_27_572_858 * h8_p5 + e_3 * fs_4_11_2 * r_2 * h6_p1 + e_3 * fs_4_11_33 * r_2 * h6_p5 - e_3 * fs_45_286_105 * r_4 * h4_p1 - e_3 * fs_10_11_14 * r_6 * h2_p1 + e_4 * fs_7_4199_2310 * h10_p1 - e_4 * fs_35_4199_1001 * h10_p5 - e_4 * fs_127_5434_42 * r_2 * h8_p1 - e_4 * fs_27_5434_858 * r_2 * h8_p5 - e_4 * fs_4_187_2 * r_4 * h6_p1 - e_4 * fs_4_187_33 * r_4 * h6_p5 + e_4 * fs_1_143_105 * r_6 * h4_p1 + e_4 * fs_5_143_14 * r_8 * h2_p1;
    }

    // NOTE: the angular components are formed in 87 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_p2, ph4_p2, ph4_p4, ph6_p2, ph6_p4, ph8_p2, ph8_p4, ph10_p2, ph10_p4, ab_2, pc_30 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

        const auto h2_p2 = ph2_p2[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h4_p4 = ph4_p4[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h6_p4 = ph6_p4[k];
        const auto h8_p2 = ph8_p2[k];
        const auto h8_p4 = ph8_p4[k];
        const auto h10_p2 = ph10_p2[k];
        const auto h10_p4 = ph10_p4[k];

        pc_30[k] = e_0 * fs_315_16_7 * h2_p2 - e_1 * fs_15_8_105 * h4_p2 + e_1 * fs_15_2_15 * h4_p4 - e_1 * fs_45_2_7 * r_2 * h2_p2 + e_2 * fs_195_88_10 * h6_p2 - e_2 * fs_105_44_3 * h6_p4 + e_2 * fs_45_44_105 * r_2 * h4_p2 - e_2 * fs_45_11_15 * r_2 * h4_p4 + e_2 * fs_15_2_7 * r_4 * h2_p2 - e_3 * fs_111_572_30 * h8_p2 + e_3 * fs_43_286_33 * h8_p4 - e_3 * fs_13_22_10 * r_2 * h6_p2 + e_3 * fs_7_11_3 * r_2 * h6_p4 - e_3 * fs_45_286_105 * r_4 * h4_p2 + e_3 * fs_90_143_15 * r_4 * h4_p4 - e_3 * fs_10_11_7 * r_6 * h2_p2 - e_4 * fs_28_4199_385 * h10_p2 - e_4 * fs_14_4199_5005 * h10_p4 + e_4 * fs_111_5434_30 * r_2 * h8_p2 - e_4 * fs_43_2717_33 * r_2 * h8_p4 + e_4 * fs_13_374_10 * r_4 * h6_p2 - e_4 * fs_7_187_3 * r_4 * h6_p4 + e_4 * fs_1_143_105 * r_6 * h4_p2 - e_4 * fs_4_143_15 * r_6 * h4_p4 + e_4 * fs_5_143_7 * r_8 * h2_p2;
    }

    // NOTE: the angular components are formed in 87 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, ph4_m3, ph6_m3, ph8_m3, ph10_m3, ab_2, pc_31 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h4_m3 = ph4_m3[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h8_m3 = ph8_m3[k];
        const auto h10_m3 = ph10_m3[k];

        pc_31[k] = e_1 * fs_75_4_3 * h4_m3 - e_2 * f_135_22 * h6_m3 - e_2 * fs_225_22_3 * r_2 * h4_m3 + e_3 * fs_15_286_22 * h8_m3 + e_3 * f_18_11 * r_2 * h6_m3 + e_3 * fs_225_143_3 * r_4 * h4_m3 + e_4 * fs_35_4199_1001 * h10_m3 - e_4 * fs_15_2717_22 * r_2 * h8_m3 - e_4 * f_18_187 * r_4 * h6_m3 - e_4 * fs_10_143_3 * r_6 * h4_m3;
    }

    // NOTE: the angular components are formed in 87 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_m2, ph4_m4, ph4_m2, ph6_m4, ph6_m2, ph8_m4, ph8_m2, ph10_m4, ph10_m2, ab_2, pc_32 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

        const auto h2_m2 = ph2_m2[k];
        const auto h4_m4 = ph4_m4[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h6_m4 = ph6_m4[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h8_m4 = ph8_m4[k];
        const auto h8_m2 = ph8_m2[k];
        const auto h10_m4 = ph10_m4[k];
        const auto h10_m2 = ph10_m2[k];

        pc_32[k] = e_0 * fs_315_16_7 * h2_m2 - e_1 * fs_15_2_15 * h4_m4 - e_1 * fs_15_8_105 * h4_m2 - e_1 * fs_45_2_7 * r_2 * h2_m2 + e_2 * fs_105_44_3 * h6_m4 + e_2 * fs_195_88_10 * h6_m2 + e_2 * fs_45_11_15 * r_2 * h4_m4 + e_2 * fs_45_44_105 * r_2 * h4_m2 + e_2 * fs_15_2_7 * r_4 * h2_m2 - e_3 * fs_43_286_33 * h8_m4 - e_3 * fs_111_572_30 * h8_m2 - e_3 * fs_7_11_3 * r_2 * h6_m4 - e_3 * fs_13_22_10 * r_2 * h6_m2 - e_3 * fs_90_143_15 * r_4 * h4_m4 - e_3 * fs_45_286_105 * r_4 * h4_m2 - e_3 * fs_10_11_7 * r_6 * h2_m2 + e_4 * fs_14_4199_5005 * h10_m4 - e_4 * fs_28_4199_385 * h10_m2 + e_4 * fs_43_2717_33 * r_2 * h8_m4 + e_4 * fs_111_5434_30 * r_2 * h8_m2 + e_4 * fs_7_187_3 * r_4 * h6_m4 + e_4 * fs_13_374_10 * r_4 * h6_m2 + e_4 * fs_4_143_15 * r_6 * h4_m4 + e_4 * fs_1_143_105 * r_6 * h4_m2 + e_4 * fs_5_143_7 * r_8 * h2_m2;
    }

    // NOTE: the angular components are formed in 87 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_m1, ph4_m1, ph6_m6, ph6_m5, ph6_m1, ph8_m6, ph8_m5, ph8_m1, ph10_m6, ph10_m5, ph10_m1, ab_2, pc_33, pc_34 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

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

        pc_33[k] = e_0 * fs_315_16_14 * h2_m1 - e_1 * fs_15_8_105 * h4_m1 - e_1 * fs_45_2_14 * r_2 * h2_m1 + e_2 * fs_15_11_33 * h6_m5 - e_2 * fs_15_11_2 * h6_m1 + e_2 * fs_45_44_105 * r_2 * h4_m1 + e_2 * fs_15_2_14 * r_4 * h2_m1 - e_3 * fs_27_572_858 * h8_m5 + e_3 * fs_127_572_42 * h8_m1 - e_3 * fs_4_11_33 * r_2 * h6_m5 + e_3 * fs_4_11_2 * r_2 * h6_m1 - e_3 * fs_45_286_105 * r_4 * h4_m1 - e_3 * fs_10_11_14 * r_6 * h2_m1 + e_4 * fs_35_4199_1001 * h10_m5 + e_4 * fs_7_4199_2310 * h10_m1 + e_4 * fs_27_5434_858 * r_2 * h8_m5 - e_4 * fs_127_5434_42 * r_2 * h8_m1 + e_4 * fs_4_187_33 * r_4 * h6_m5 - e_4 * fs_4_187_2 * r_4 * h6_m1 + e_4 * fs_1_143_105 * r_6 * h4_m1 + e_4 * fs_5_143_14 * r_8 * h2_m1;

        pc_34[k] = - e_2 * fs_45_88_154 * h6_m6 - e_3 * fs_21_572_286 * h8_m6 + e_3 * fs_3_22_154 * r_2 * h6_m6 + e_4 * fs_28_4199_1430 * h10_m6 + e_4 * fs_21_5434_286 * r_2 * h8_m6 - e_4 * fs_3_374_154 * r_4 * h6_m6;
    }

    // NOTE: the angular components are formed in 87 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_m1, ph4_m1, ph6_m1, ph8_m7, ph8_m1, ph10_m7, ph10_m1, ab_2, pc_35 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

        const auto h2_m1 = ph2_m1[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h8_m7 = ph8_m7[k];
        const auto h8_m1 = ph8_m1[k];
        const auto h10_m7 = ph10_m7[k];
        const auto h10_m1 = ph10_m1[k];

        pc_35[k] = e_0 * fs_315_32_2 * h2_m1 + e_1 * fs_15_2_15 * h4_m1 - e_1 * fs_45_4_2 * r_2 * h2_m1 + e_2 * fs_105_44_14 * h6_m1 - e_2 * fs_45_11_15 * r_2 * h4_m1 + e_2 * fs_15_4_2 * r_4 * h2_m1 + e_3 * fs_7_286_4290 * h8_m7 + e_3 * fs_42_143_6 * h8_m1 - e_3 * fs_7_11_14 * r_2 * h6_m1 + e_3 * fs_90_143_15 * r_4 * h4_m1 - e_3 * fs_5_11_2 * r_6 * h2_m1 + e_4 * fs_7_4199_12155 * h10_m7 + e_4 * fs_7_8398_330 * h10_m1 - e_4 * fs_7_2717_4290 * r_2 * h8_m7 - e_4 * fs_84_2717_6 * r_2 * h8_m1 + e_4 * fs_7_187_14 * r_4 * h6_m1 - e_4 * fs_4_143_15 * r_6 * h4_m1 + e_4 * fs_5_286_2 * r_8 * h2_m1;
    }

    // NOTE: the angular components are formed in 87 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_p2, ph4_p2, ph6_p2, ph6_p6, ph8_p2, ph8_p6, ph10_p2, ph10_p6, ab_2, pc_36 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

        const auto h2_p2 = ph2_p2[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h6_p6 = ph6_p6[k];
        const auto h8_p2 = ph8_p2[k];
        const auto h8_p6 = ph8_p6[k];
        const auto h10_p2 = ph10_p2[k];
        const auto h10_p6 = ph10_p6[k];

        pc_36[k] = e_0 * fs_105_32_2 * h2_p2 + e_1 * fs_15_4_30 * h4_p2 - e_1 * fs_15_4_2 * r_2 * h2_p2 + e_2 * fs_35_22_35 * h6_p2 - e_2 * fs_15_44_77 * h6_p6 - e_2 * fs_45_22_30 * r_2 * h4_p2 + e_2 * fs_5_4_2 * r_4 * h2_p2 + e_3 * fs_14_143_105 * h8_p2 - e_3 * fs_21_143_143 * h8_p6 - e_3 * fs_14_33_35 * r_2 * h6_p2 + e_3 * fs_1_11_77 * r_2 * h6_p6 + e_3 * fs_45_143_30 * r_4 * h4_p2 - e_3 * fs_5_33_2 * r_6 * h2_p2 + e_4 * fs_21_8398_110 * h10_p2 - e_4 * fs_21_4199_715 * h10_p6 - e_4 * fs_28_2717_105 * r_2 * h8_p2 + e_4 * fs_42_2717_143 * r_2 * h8_p6 + e_4 * fs_14_561_35 * r_4 * h6_p2 - e_4 * fs_1_187_77 * r_4 * h6_p6 - e_4 * fs_2_143_30 * r_6 * h4_p2 + e_4 * fs_5_858_2 * r_8 * h2_p2;
    }

    // NOTE: the angular components are formed in 87 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_p1, ph4_p1, ph6_p1, ph6_p5, ph8_p1, ph8_p5, ph10_p1, ph10_p5, ab_2, pc_37 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

        const auto h2_p1 = ph2_p1[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h6_p5 = ph6_p5[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h8_p5 = ph8_p5[k];
        const auto h10_p1 = ph10_p1[k];
        const auto h10_p5 = ph10_p5[k];

        pc_37[k] = - e_0 * f_105_4 * h2_p1 - e_1 * fs_45_8_30 * h4_p1 + e_1 * f_30_1 * r_2 * h2_p1 + e_2 * fs_35_44_7 * h6_p1 + e_2 * fs_35_88_462 * h6_p5 + e_2 * fs_135_44_30 * r_2 * h4_p1 - e_2 * f_10_1 * r_4 * h2_p1 + e_3 * fs_119_143_3 * h8_p1 - e_3 * fs_1_286_3003 * h8_p5 - e_3 * fs_7_33_7 * r_2 * h6_p1 - e_3 * fs_7_66_462 * r_2 * h6_p5 - e_3 * fs_135_286_30 * r_4 * h4_p1 + e_3 * f_40_33 * r_6 * h2_p1 + e_4 * fs_21_4199_165 * h10_p1 - e_4 * fs_105_8398_286 * h10_p5 - e_4 * fs_238_2717_3 * r_2 * h8_p1 + e_4 * fs_1_2717_3003 * r_2 * h8_p5 + e_4 * fs_7_561_7 * r_4 * h6_p1 + e_4 * fs_7_1122_462 * r_4 * h6_p5 + e_4 * fs_3_143_30 * r_6 * h4_p1 - e_4 * f_20_429 * r_8 * h2_p1;
    }

    // NOTE: the angular components are formed in 87 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_0, ph4_p4, ph6_0, ph6_p4, ph8_0, ph8_p4, ph10_0, ph10_p4, ab_2, pc_38 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

        const auto h2_0 = ph2_0[k];
        const auto h4_p4 = ph4_p4[k];
        const auto h6_0 = ph6_0[k];
        const auto h6_p4 = ph6_p4[k];
        const auto h8_0 = ph8_0[k];
        const auto h8_p4 = ph8_p4[k];
        const auto h10_0 = ph10_0[k];
        const auto h10_p4 = ph10_p4[k];

        pc_38[k] = e_0 * fs_105_8_42 * h2_0 - e_1 * fs_15_4_30 * h4_p4 - e_1 * fs_15_1_42 * r_2 * h2_0 - e_2 * fs_5_4_42 * h6_0 - e_2 * fs_115_88_6 * h6_p4 + e_2 * fs_45_22_30 * r_2 * h4_p4 + e_2 * fs_5_1_42 * r_4 * h2_0 + e_3 * fs_3_11_42 * h8_0 + e_3 * fs_37_286_66 * h8_p4 + e_3 * fs_1_3_42 * r_2 * h6_0 + e_3 * fs_23_66_6 * r_2 * h6_p4 - e_3 * fs_45_143_30 * r_4 * h4_p4 - e_3 * fs_20_33_42 * r_6 * h2_0 + e_4 * fs_105_4199_42 * h10_0 - e_4 * fs_21_8398_10010 * h10_p4 - e_4 * fs_6_209_42 * r_2 * h8_0 - e_4 * fs_37_2717_66 * r_2 * h8_p4 - e_4 * fs_1_51_42 * r_4 * h6_0 - e_4 * fs_23_1122_6 * r_4 * h6_p4 + e_4 * fs_2_143_30 * r_6 * h4_p4 + e_4 * fs_10_429_42 * r_8 * h2_0;
    }

    // NOTE: the angular components are formed in 87 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_p1, ph4_p1, ph4_p3, ph6_p1, ph6_p3, ph8_p1, ph8_p3, ph10_p1, ph10_p3, ab_2, pc_39 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

        const auto h2_p1 = ph2_p1[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h4_p3 = ph4_p3[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h8_p3 = ph8_p3[k];
        const auto h10_p1 = ph10_p1[k];
        const auto h10_p3 = ph10_p3[k];

        pc_39[k] = e_0 * fs_105_4_7 * h2_p1 - e_1 * fs_15_8_210 * h4_p1 + e_1 * fs_45_8_30 * h4_p3 - e_1 * fs_30_1_7 * r_2 * h2_p1 + e_2 * f_265_44 * h6_p1 - e_2 * fs_195_88_10 * h6_p3 + e_2 * fs_45_44_210 * r_2 * h4_p1 - e_2 * fs_135_44_30 * r_2 * h4_p3 + e_2 * fs_10_1_7 * r_4 * h2_p1 - e_3 * fs_29_286_21 * h8_p1 + e_3 * fs_24_143_55 * h8_p3 - e_3 * f_53_33 * r_2 * h6_p1 + e_3 * fs_13_22_10 * r_2 * h6_p3 - e_3 * fs_45_286_210 * r_4 * h4_p1 + e_3 * fs_135_286_30 * r_4 * h4_p3 - e_3 * fs_40_33_7 * r_6 * h2_p1 - e_4 * fs_21_4199_1155 * h10_p1 - e_4 * fs_21_8398_10010 * h10_p3 + e_4 * fs_29_2717_21 * r_2 * h8_p1 - e_4 * fs_48_2717_55 * r_2 * h8_p3 + e_4 * f_53_561 * r_4 * h6_p1 - e_4 * fs_13_374_10 * r_4 * h6_p3 + e_4 * fs_1_143_210 * r_6 * h4_p1 - e_4 * fs_3_143_30 * r_6 * h4_p3 + e_4 * fs_20_429_7 * r_8 * h2_p1;
    }

    // NOTE: the angular components are formed in 87 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_m2, ph6_m2, ph8_m2, ph10_m2, ab_2, pc_40 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

        const auto h2_m2 = ph2_m2[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h8_m2 = ph8_m2[k];
        const auto h10_m2 = ph10_m2[k];

        pc_40[k] = e_0 * fs_105_16_70 * h2_m2 - e_1 * fs_15_2_70 * r_2 * h2_m2 + e_2 * f_5_4 * h6_m2 + e_2 * fs_5_2_70 * r_4 * h2_m2 - e_3 * fs_5_11_3 * h8_m2 - e_3 * f_1_3 * r_2 * h6_m2 - e_3 * fs_10_33_70 * r_6 * h2_m2 + e_4 * fs_105_4199_154 * h10_m2 + e_4 * fs_10_209_3 * r_2 * h8_m2 + e_4 * f_1_51 * r_4 * h6_m2 + e_4 * fs_5_429_70 * r_8 * h2_m2;
    }

    // NOTE: the angular components are formed in 87 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_m1, ph4_m3, ph4_m1, ph6_m3, ph6_m1, ph8_m3, ph8_m1, ph10_m3, ph10_m1, ab_2, pc_41 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

        const auto h2_m1 = ph2_m1[k];
        const auto h4_m3 = ph4_m3[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h8_m3 = ph8_m3[k];
        const auto h8_m1 = ph8_m1[k];
        const auto h10_m3 = ph10_m3[k];
        const auto h10_m1 = ph10_m1[k];

        pc_41[k] = e_0 * fs_105_4_7 * h2_m1 - e_1 * fs_45_8_30 * h4_m3 - e_1 * fs_15_8_210 * h4_m1 - e_1 * fs_30_1_7 * r_2 * h2_m1 + e_2 * fs_195_88_10 * h6_m3 + e_2 * f_265_44 * h6_m1 + e_2 * fs_135_44_30 * r_2 * h4_m3 + e_2 * fs_45_44_210 * r_2 * h4_m1 + e_2 * fs_10_1_7 * r_4 * h2_m1 - e_3 * fs_24_143_55 * h8_m3 - e_3 * fs_29_286_21 * h8_m1 - e_3 * fs_13_22_10 * r_2 * h6_m3 - e_3 * f_53_33 * r_2 * h6_m1 - e_3 * fs_135_286_30 * r_4 * h4_m3 - e_3 * fs_45_286_210 * r_4 * h4_m1 - e_3 * fs_40_33_7 * r_6 * h2_m1 + e_4 * fs_21_8398_10010 * h10_m3 - e_4 * fs_21_4199_1155 * h10_m1 + e_4 * fs_48_2717_55 * r_2 * h8_m3 + e_4 * fs_29_2717_21 * r_2 * h8_m1 + e_4 * fs_13_374_10 * r_4 * h6_m3 + e_4 * f_53_561 * r_4 * h6_m1 + e_4 * fs_3_143_30 * r_6 * h4_m3 + e_4 * fs_1_143_210 * r_6 * h4_m1 + e_4 * fs_20_429_7 * r_8 * h2_m1;
    }

    // NOTE: the angular components are formed in 87 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_m1, ph4_m4, ph4_m1, ph6_m5, ph6_m4, ph6_m1, ph8_m5, ph8_m4, ph8_m1, ph10_m5, ph10_m4, ph10_m1, ab_2, pc_42, pc_43 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

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

        pc_42[k] = e_1 * fs_15_4_30 * h4_m4 + e_2 * fs_115_88_6 * h6_m4 - e_2 * fs_45_22_30 * r_2 * h4_m4 - e_3 * fs_37_286_66 * h8_m4 - e_3 * fs_23_66_6 * r_2 * h6_m4 + e_3 * fs_45_143_30 * r_4 * h4_m4 + e_4 * fs_21_8398_10010 * h10_m4 + e_4 * fs_37_2717_66 * r_2 * h8_m4 + e_4 * fs_23_1122_6 * r_4 * h6_m4 - e_4 * fs_2_143_30 * r_6 * h4_m4;

        pc_43[k] = e_0 * f_105_4 * h2_m1 + e_1 * fs_45_8_30 * h4_m1 - e_1 * f_30_1 * r_2 * h2_m1 - e_2 * fs_35_88_462 * h6_m5 - e_2 * fs_35_44_7 * h6_m1 - e_2 * fs_135_44_30 * r_2 * h4_m1 + e_2 * f_10_1 * r_4 * h2_m1 + e_3 * fs_1_286_3003 * h8_m5 - e_3 * fs_119_143_3 * h8_m1 + e_3 * fs_7_66_462 * r_2 * h6_m5 + e_3 * fs_7_33_7 * r_2 * h6_m1 + e_3 * fs_135_286_30 * r_4 * h4_m1 - e_3 * f_40_33 * r_6 * h2_m1 + e_4 * fs_105_8398_286 * h10_m5 - e_4 * fs_21_4199_165 * h10_m1 - e_4 * fs_1_2717_3003 * r_2 * h8_m5 + e_4 * fs_238_2717_3 * r_2 * h8_m1 - e_4 * fs_7_1122_462 * r_4 * h6_m5 - e_4 * fs_7_561_7 * r_4 * h6_m1 - e_4 * fs_3_143_30 * r_6 * h4_m1 + e_4 * f_20_429 * r_8 * h2_m1;
    }

    // NOTE: the angular components are formed in 87 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_m2, ph4_m2, ph6_m6, ph6_m2, ph8_m6, ph8_m2, ph10_m6, ph10_m2, ab_2, pc_44 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

        const auto h2_m2 = ph2_m2[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h6_m6 = ph6_m6[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h8_m6 = ph8_m6[k];
        const auto h8_m2 = ph8_m2[k];
        const auto h10_m6 = ph10_m6[k];
        const auto h10_m2 = ph10_m2[k];

        pc_44[k] = - e_0 * fs_105_32_2 * h2_m2 - e_1 * fs_15_4_30 * h4_m2 + e_1 * fs_15_4_2 * r_2 * h2_m2 + e_2 * fs_15_44_77 * h6_m6 - e_2 * fs_35_22_35 * h6_m2 + e_2 * fs_45_22_30 * r_2 * h4_m2 - e_2 * fs_5_4_2 * r_4 * h2_m2 + e_3 * fs_21_143_143 * h8_m6 - e_3 * fs_14_143_105 * h8_m2 - e_3 * fs_1_11_77 * r_2 * h6_m6 + e_3 * fs_14_33_35 * r_2 * h6_m2 - e_3 * fs_45_143_30 * r_4 * h4_m2 + e_3 * fs_5_33_2 * r_6 * h2_m2 + e_4 * fs_21_4199_715 * h10_m6 - e_4 * fs_21_8398_110 * h10_m2 - e_4 * fs_42_2717_143 * r_2 * h8_m6 + e_4 * fs_28_2717_105 * r_2 * h8_m2 + e_4 * fs_1_187_77 * r_4 * h6_m6 - e_4 * fs_14_561_35 * r_4 * h6_m2 + e_4 * fs_2_143_30 * r_6 * h4_m2 - e_4 * fs_5_858_2 * r_8 * h2_m2;
    }

    // NOTE: the angular components are formed in 87 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, ph4_p3, ph6_p3, ph6_p5, ph8_p3, ph8_p5, ph10_p3, ph10_p5, ab_2, pc_45 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h4_p3 = ph4_p3[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h6_p5 = ph6_p5[k];
        const auto h8_p3 = ph8_p3[k];
        const auto h8_p5 = ph8_p5[k];
        const auto h10_p3 = ph10_p3[k];
        const auto h10_p5 = ph10_p5[k];

        pc_45[k] = - e_1 * fs_15_8_42 * h4_p3 - e_2 * fs_105_44_14 * h6_p3 - e_2 * fs_5_44_2310 * h6_p5 + e_2 * fs_45_44_42 * r_2 * h4_p3 - e_3 * fs_21_143_77 * h8_p3 - e_3 * fs_2_143_15015 * h8_p5 + e_3 * fs_7_11_14 * r_2 * h6_p3 + e_3 * fs_1_33_2310 * r_2 * h6_p5 - e_3 * fs_45_286_42 * r_4 * h4_p3 - e_4 * fs_21_8398_286 * h10_p3 - e_4 * fs_21_8398_1430 * h10_p5 + e_4 * fs_42_2717_77 * r_2 * h8_p3 + e_4 * fs_4_2717_15015 * r_2 * h8_p5 - e_4 * fs_7_187_14 * r_4 * h6_p3 - e_4 * fs_1_561_2310 * r_4 * h6_p5 + e_4 * fs_1_143_42 * r_6 * h4_p3;
    }

    // NOTE: the angular components are formed in 87 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_p2, ph4_p2, ph4_p4, ph6_p2, ph6_p4, ph8_p2, ph8_p4, ph10_p2, ph10_p4, ab_2, pc_46 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

        const auto h2_p2 = ph2_p2[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h4_p4 = ph4_p4[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h6_p4 = ph6_p4[k];
        const auto h8_p2 = ph8_p2[k];
        const auto h8_p4 = ph8_p4[k];
        const auto h10_p2 = ph10_p2[k];
        const auto h10_p4 = ph10_p4[k];

        pc_46[k] = e_0 * fs_105_32_10 * h2_p2 + e_1 * fs_195_16_6 * h4_p2 + e_1 * fs_15_8_42 * h4_p4 - e_1 * fs_15_4_10 * r_2 * h2_p2 + e_2 * fs_35_44_7 * h6_p2 + e_2 * fs_25_44_210 * h6_p4 - e_2 * fs_585_88_6 * r_2 * h4_p2 - e_2 * fs_45_44_42 * r_2 * h4_p4 + e_2 * fs_5_4_10 * r_4 * h2_p2 - e_3 * fs_7_22_21 * h8_p2 - e_3 * fs_5_286_2310 * h8_p4 - e_3 * fs_7_33_7 * r_2 * h6_p2 - e_3 * fs_5_33_210 * r_2 * h6_p4 + e_3 * fs_45_44_6 * r_4 * h4_p2 + e_3 * fs_45_286_42 * r_4 * h4_p4 - e_3 * fs_5_33_10 * r_6 * h2_p2 - e_4 * fs_84_4199_22 * h10_p2 - e_4 * fs_42_4199_286 * h10_p4 + e_4 * fs_7_209_21 * r_2 * h8_p2 + e_4 * fs_5_2717_2310 * r_2 * h8_p4 + e_4 * fs_7_561_7 * r_4 * h6_p2 + e_4 * fs_5_561_210 * r_4 * h6_p4 - e_4 * fs_1_22_6 * r_6 * h4_p2 - e_4 * fs_1_143_42 * r_6 * h4_p4 + e_4 * fs_5_858_10 * r_8 * h2_p2;
    }

    // NOTE: the angular components are formed in 87 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_p1, ph4_p1, ph4_p3, ph6_p1, ph6_p3, ph8_p1, ph8_p3, ph10_p1, ph10_p3, ab_2, pc_47 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

        const auto h2_p1 = ph2_p1[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h4_p3 = ph4_p3[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h8_p3 = ph8_p3[k];
        const auto h10_p1 = ph10_p1[k];
        const auto h10_p3 = ph10_p3[k];

        pc_47[k] = - e_0 * fs_105_16_35 * h2_p1 - e_1 * fs_45_16_42 * h4_p1 - e_1 * fs_195_16_6 * h4_p3 + e_1 * fs_15_2_35 * r_2 * h2_p1 + e_2 * fs_35_11_5 * h6_p1 + e_2 * fs_15_11_2 * h6_p3 + e_2 * fs_135_88_42 * r_2 * h4_p1 + e_2 * fs_585_88_6 * r_2 * h4_p3 - e_2 * fs_5_2_35 * r_4 * h2_p1 - e_3 * fs_23_286_105 * h8_p1 + e_3 * fs_3_22_11 * h8_p3 - e_3 * fs_28_33_5 * r_2 * h6_p1 - e_3 * fs_4_11_2 * r_2 * h6_p3 - e_3 * fs_135_572_42 * r_4 * h4_p1 - e_3 * fs_45_44_6 * r_4 * h4_p3 + e_3 * fs_10_33_35 * r_6 * h2_p1 - e_4 * fs_42_4199_231 * h10_p1 - e_4 * fs_21_4199_2002 * h10_p3 + e_4 * fs_23_2717_105 * r_2 * h8_p1 - e_4 * fs_3_209_11 * r_2 * h8_p3 + e_4 * fs_28_561_5 * r_4 * h6_p1 + e_4 * fs_4_187_2 * r_4 * h6_p3 + e_4 * fs_3_286_42 * r_6 * h4_p1 + e_4 * fs_1_22_6 * r_6 * h4_p3 - e_4 * fs_5_429_35 * r_8 * h2_p1;
    }

    // NOTE: the angular components are formed in 87 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_0, ph2_p2, ph4_0, ph4_p2, ph6_0, ph6_p2, ph8_0, ph8_p2, ph10_0, ph10_p2, ab_2, pc_48 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

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

        pc_48[k] = e_0 * fs_105_16_210 * h2_0 + e_0 * fs_105_32_70 * h2_p2 - e_1 * fs_15_8_210 * h4_0 + e_1 * fs_45_16_42 * h4_p2 - e_1 * fs_15_2_210 * r_2 * h2_0 - e_1 * fs_15_4_70 * r_2 * h2_p2 + e_2 * fs_5_22_210 * h6_0 - e_2 * f_265_44 * h6_p2 + e_2 * fs_45_44_210 * r_2 * h4_0 - e_2 * fs_135_88_42 * r_2 * h4_p2 + e_2 * fs_5_2_210 * r_4 * h2_0 + e_2 * fs_5_4_70 * r_4 * h2_p2 + e_3 * fs_3_143_210 * h8_0 + e_3 * fs_199_286_3 * h8_p2 - e_3 * fs_2_33_210 * r_2 * h6_0 + e_3 * f_53_33 * r_2 * h6_p2 - e_3 * fs_45_286_210 * r_4 * h4_0 + e_3 * fs_135_572_42 * r_4 * h4_p2 - e_3 * fs_10_33_210 * r_6 * h2_0 - e_3 * fs_5_33_70 * r_6 * h2_p2 - e_4 * fs_84_4199_210 * h10_0 - e_4 * fs_84_4199_154 * h10_p2 - e_4 * fs_6_2717_210 * r_2 * h8_0 - e_4 * fs_199_2717_3 * r_2 * h8_p2 + e_4 * fs_2_561_210 * r_4 * h6_0 - e_4 * f_53_561 * r_4 * h6_p2 + e_4 * fs_1_143_210 * r_6 * h4_0 - e_4 * fs_3_286_42 * r_6 * h4_p2 + e_4 * fs_5_429_210 * r_8 * h2_0 + e_4 * fs_5_858_70 * r_8 * h2_p2;
    }

    // NOTE: the angular components are formed in 87 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_m2, ph2_m1, ph4_m2, ph4_m1, ph6_m2, ph6_m1, ph8_m2, ph8_m1, ph10_m2, ph10_m1, ab_2, pc_49, pc_50 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

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

        pc_49[k] = e_0 * fs_525_16_7 * h2_m1 - e_1 * fs_15_8_210 * h4_m1 - e_1 * fs_75_2_7 * r_2 * h2_m1 + e_2 * f_80_11 * h6_m1 + e_2 * fs_45_44_210 * r_2 * h4_m1 + e_2 * fs_25_2_7 * r_4 * h2_m1 - e_3 * fs_47_143_21 * h8_m1 - e_3 * f_64_33 * r_2 * h6_m1 - e_3 * fs_45_286_210 * r_4 * h4_m1 - e_3 * fs_50_33_7 * r_6 * h2_m1 + e_4 * fs_42_4199_1155 * h10_m1 + e_4 * fs_94_2717_21 * r_2 * h8_m1 + e_4 * f_64_561 * r_4 * h6_m1 + e_4 * fs_1_143_210 * r_6 * h4_m1 + e_4 * fs_25_429_7 * r_8 * h2_m1;

        pc_50[k] = - e_0 * fs_105_32_70 * h2_m2 - e_1 * fs_45_16_42 * h4_m2 + e_1 * fs_15_4_70 * r_2 * h2_m2 + e_2 * f_265_44 * h6_m2 + e_2 * fs_135_88_42 * r_2 * h4_m2 - e_2 * fs_5_4_70 * r_4 * h2_m2 - e_3 * fs_199_286_3 * h8_m2 - e_3 * f_53_33 * r_2 * h6_m2 - e_3 * fs_135_572_42 * r_4 * h4_m2 + e_3 * fs_5_33_70 * r_6 * h2_m2 + e_4 * fs_84_4199_154 * h10_m2 + e_4 * fs_199_2717_3 * r_2 * h8_m2 + e_4 * f_53_561 * r_4 * h6_m2 + e_4 * fs_3_286_42 * r_6 * h4_m2 - e_4 * fs_5_858_70 * r_8 * h2_m2;
    }

    // NOTE: the angular components are formed in 87 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_m1, ph4_m3, ph4_m1, ph6_m3, ph6_m1, ph8_m3, ph8_m1, ph10_m3, ph10_m1, ab_2, pc_51 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

        const auto h2_m1 = ph2_m1[k];
        const auto h4_m3 = ph4_m3[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h8_m3 = ph8_m3[k];
        const auto h8_m1 = ph8_m1[k];
        const auto h10_m3 = ph10_m3[k];
        const auto h10_m1 = ph10_m1[k];

        pc_51[k] = e_0 * fs_105_16_35 * h2_m1 + e_1 * fs_195_16_6 * h4_m3 + e_1 * fs_45_16_42 * h4_m1 - e_1 * fs_15_2_35 * r_2 * h2_m1 - e_2 * fs_15_11_2 * h6_m3 - e_2 * fs_35_11_5 * h6_m1 - e_2 * fs_585_88_6 * r_2 * h4_m3 - e_2 * fs_135_88_42 * r_2 * h4_m1 + e_2 * fs_5_2_35 * r_4 * h2_m1 - e_3 * fs_3_22_11 * h8_m3 + e_3 * fs_23_286_105 * h8_m1 + e_3 * fs_4_11_2 * r_2 * h6_m3 + e_3 * fs_28_33_5 * r_2 * h6_m1 + e_3 * fs_45_44_6 * r_4 * h4_m3 + e_3 * fs_135_572_42 * r_4 * h4_m1 - e_3 * fs_10_33_35 * r_6 * h2_m1 + e_4 * fs_21_4199_2002 * h10_m3 + e_4 * fs_42_4199_231 * h10_m1 + e_4 * fs_3_209_11 * r_2 * h8_m3 - e_4 * fs_23_2717_105 * r_2 * h8_m1 - e_4 * fs_4_187_2 * r_4 * h6_m3 - e_4 * fs_28_561_5 * r_4 * h6_m1 - e_4 * fs_1_22_6 * r_6 * h4_m3 - e_4 * fs_3_286_42 * r_6 * h4_m1 + e_4 * fs_5_429_35 * r_8 * h2_m1;
    }

    // NOTE: the angular components are formed in 87 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_m2, ph4_m4, ph4_m2, ph6_m4, ph6_m2, ph8_m4, ph8_m2, ph10_m4, ph10_m2, ab_2, pc_52 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

        const auto h2_m2 = ph2_m2[k];
        const auto h4_m4 = ph4_m4[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h6_m4 = ph6_m4[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h8_m4 = ph8_m4[k];
        const auto h8_m2 = ph8_m2[k];
        const auto h10_m4 = ph10_m4[k];
        const auto h10_m2 = ph10_m2[k];

        pc_52[k] = - e_0 * fs_105_32_10 * h2_m2 - e_1 * fs_15_8_42 * h4_m4 - e_1 * fs_195_16_6 * h4_m2 + e_1 * fs_15_4_10 * r_2 * h2_m2 - e_2 * fs_25_44_210 * h6_m4 - e_2 * fs_35_44_7 * h6_m2 + e_2 * fs_45_44_42 * r_2 * h4_m4 + e_2 * fs_585_88_6 * r_2 * h4_m2 - e_2 * fs_5_4_10 * r_4 * h2_m2 + e_3 * fs_5_286_2310 * h8_m4 + e_3 * fs_7_22_21 * h8_m2 + e_3 * fs_5_33_210 * r_2 * h6_m4 + e_3 * fs_7_33_7 * r_2 * h6_m2 - e_3 * fs_45_286_42 * r_4 * h4_m4 - e_3 * fs_45_44_6 * r_4 * h4_m2 + e_3 * fs_5_33_10 * r_6 * h2_m2 + e_4 * fs_42_4199_286 * h10_m4 + e_4 * fs_84_4199_22 * h10_m2 - e_4 * fs_5_2717_2310 * r_2 * h8_m4 - e_4 * fs_7_209_21 * r_2 * h8_m2 - e_4 * fs_5_561_210 * r_4 * h6_m4 - e_4 * fs_7_561_7 * r_4 * h6_m2 + e_4 * fs_1_143_42 * r_6 * h4_m4 + e_4 * fs_1_22_6 * r_6 * h4_m2 - e_4 * fs_5_858_10 * r_8 * h2_m2;
    }

    // NOTE: the angular components are formed in 87 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, ph4_m4, ph4_m3, ph6_m5, ph6_m4, ph6_m3, ph8_m5, ph8_m4, ph8_m3, ph10_m5, ph10_m4, ph10_m3, ab_2, pc_53, pc_54, pc_55 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

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

        pc_53[k] = e_1 * fs_15_8_42 * h4_m3 + e_2 * fs_5_44_2310 * h6_m5 + e_2 * fs_105_44_14 * h6_m3 - e_2 * fs_45_44_42 * r_2 * h4_m3 + e_3 * fs_2_143_15015 * h8_m5 + e_3 * fs_21_143_77 * h8_m3 - e_3 * fs_1_33_2310 * r_2 * h6_m5 - e_3 * fs_7_11_14 * r_2 * h6_m3 + e_3 * fs_45_286_42 * r_4 * h4_m3 + e_4 * fs_21_8398_1430 * h10_m5 + e_4 * fs_21_8398_286 * h10_m3 - e_4 * fs_4_2717_15015 * r_2 * h8_m5 - e_4 * fs_42_2717_77 * r_2 * h8_m3 + e_4 * fs_1_561_2310 * r_4 * h6_m5 + e_4 * fs_7_187_14 * r_4 * h6_m3 - e_4 * fs_1_143_42 * r_6 * h4_m3;

        pc_54[k] = e_1 * f_15_2 * h4_m4 + e_2 * fs_105_22_5 * h6_m4 - e_2 * f_45_11 * r_2 * h4_m4 + e_3 * fs_42_143_55 * h8_m4 - e_3 * fs_14_11_5 * r_2 * h6_m4 + e_3 * f_90_143 * r_4 * h4_m4 + e_4 * fs_7_4199_3003 * h10_m4 - e_4 * fs_84_2717_55 * r_2 * h8_m4 + e_4 * fs_14_187_5 * r_4 * h6_m4 - e_4 * f_4_143 * r_6 * h4_m4;

        pc_55[k] = - e_1 * f_255_8 * h4_m3 - e_2 * fs_105_22_3 * h6_m3 + e_2 * f_765_44 * r_2 * h4_m3 + e_3 * fs_63_286_66 * h8_m3 + e_3 * fs_14_11_3 * r_2 * h6_m3 - e_3 * f_765_286 * r_4 * h4_m3 + e_4 * fs_14_4199_3003 * h10_m3 - e_4 * fs_63_2717_66 * r_2 * h8_m3 - e_4 * fs_14_187_3 * r_4 * h6_m3 + e_4 * f_17_143 * r_6 * h4_m3;
    }

    // NOTE: the angular components are formed in 87 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_m2, ph2_m1, ph4_m2, ph4_m1, ph6_m2, ph6_m1, ph8_m2, ph8_m1, ph10_m2, ph10_m1, ab_2, pc_56, pc_57 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

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

        pc_56[k] = e_0 * fs_105_16_15 * h2_m2 + e_1 * f_165_4 * h4_m2 - e_1 * fs_15_2_15 * r_2 * h2_m2 - e_2 * fs_5_4_42 * h6_m2 - e_2 * f_45_2 * r_2 * h4_m2 + e_2 * fs_5_2_15 * r_4 * h2_m2 + e_3 * fs_12_143_14 * h8_m2 + e_3 * fs_1_3_42 * r_2 * h6_m2 + e_3 * f_45_13 * r_4 * h4_m2 - e_3 * fs_10_33_15 * r_6 * h2_m2 + e_4 * fs_196_4199_33 * h10_m2 - e_4 * fs_24_2717_14 * r_2 * h8_m2 - e_4 * fs_1_51_42 * r_4 * h6_m2 - e_4 * f_2_13 * r_6 * h4_m2 + e_4 * fs_5_429_15 * r_8 * h2_m2;

        pc_57[k] = - e_0 * fs_105_8_30 * h2_m1 + e_1 * f_15_8 * h4_m1 + e_1 * fs_15_1_30 * r_2 * h2_m1 + e_2 * fs_5_22_210 * h6_m1 - e_2 * f_45_44 * r_2 * h4_m1 - e_2 * fs_5_1_30 * r_4 * h2_m1 - e_3 * fs_105_286_10 * h8_m1 - e_3 * fs_2_33_210 * r_2 * h6_m1 + e_3 * f_45_286 * r_4 * h4_m1 + e_3 * fs_20_33_30 * r_6 * h2_m1 + e_4 * fs_294_4199_22 * h10_m1 + e_4 * fs_105_2717_10 * r_2 * h8_m1 + e_4 * fs_2_561_210 * r_4 * h6_m1 - e_4 * f_1_143 * r_6 * h4_m1 - e_4 * fs_10_429_30 * r_8 * h2_m1;
    }

    // NOTE: the angular components are formed in 87 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_0, ph2_p1, ph4_0, ph4_p1, ph6_0, ph6_p1, ph8_0, ph8_p1, ph10_0, ph10_p1, ab_2, pc_58, pc_59 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

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

        pc_58[k] = e_0 * f_1575_16 * h2_0 - e_1 * f_75_2 * h4_0 - e_1 * f_225_2 * r_2 * h2_0 + e_2 * f_105_11 * h6_0 + e_2 * f_225_11 * r_2 * h4_0 + e_2 * f_75_2 * r_4 * h2_0 - e_3 * f_252_143 * h8_0 - e_3 * f_28_11 * r_2 * h6_0 - e_3 * f_450_143 * r_4 * h4_0 - e_3 * f_50_11 * r_6 * h2_0 + e_4 * f_1470_4199 * h10_0 + e_4 * f_504_2717 * r_2 * h8_0 + e_4 * f_28_187 * r_4 * h6_0 + e_4 * f_20_143 * r_6 * h4_0 + e_4 * f_25_143 * r_8 * h2_0;

        pc_59[k] = - e_0 * fs_105_8_30 * h2_p1 + e_1 * f_15_8 * h4_p1 + e_1 * fs_15_1_30 * r_2 * h2_p1 + e_2 * fs_5_22_210 * h6_p1 - e_2 * f_45_44 * r_2 * h4_p1 - e_2 * fs_5_1_30 * r_4 * h2_p1 - e_3 * fs_105_286_10 * h8_p1 - e_3 * fs_2_33_210 * r_2 * h6_p1 + e_3 * f_45_286 * r_4 * h4_p1 + e_3 * fs_20_33_30 * r_6 * h2_p1 + e_4 * fs_294_4199_22 * h10_p1 + e_4 * fs_105_2717_10 * r_2 * h8_p1 + e_4 * fs_2_561_210 * r_4 * h6_p1 - e_4 * f_1_143 * r_6 * h4_p1 - e_4 * fs_10_429_30 * r_8 * h2_p1;
    }

    // NOTE: the angular components are formed in 87 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_p2, ph4_p2, ph4_p3, ph6_p2, ph6_p3, ph8_p2, ph8_p3, ph10_p2, ph10_p3, ab_2, pc_60, pc_61 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

        const auto h2_p2 = ph2_p2[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h4_p3 = ph4_p3[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h8_p2 = ph8_p2[k];
        const auto h8_p3 = ph8_p3[k];
        const auto h10_p2 = ph10_p2[k];
        const auto h10_p3 = ph10_p3[k];

        pc_60[k] = e_0 * fs_105_16_15 * h2_p2 + e_1 * f_165_4 * h4_p2 - e_1 * fs_15_2_15 * r_2 * h2_p2 - e_2 * fs_5_4_42 * h6_p2 - e_2 * f_45_2 * r_2 * h4_p2 + e_2 * fs_5_2_15 * r_4 * h2_p2 + e_3 * fs_12_143_14 * h8_p2 + e_3 * fs_1_3_42 * r_2 * h6_p2 + e_3 * f_45_13 * r_4 * h4_p2 - e_3 * fs_10_33_15 * r_6 * h2_p2 + e_4 * fs_196_4199_33 * h10_p2 - e_4 * fs_24_2717_14 * r_2 * h8_p2 - e_4 * fs_1_51_42 * r_4 * h6_p2 - e_4 * f_2_13 * r_6 * h4_p2 + e_4 * fs_5_429_15 * r_8 * h2_p2;

        pc_61[k] = - e_1 * f_255_8 * h4_p3 - e_2 * fs_105_22_3 * h6_p3 + e_2 * f_765_44 * r_2 * h4_p3 + e_3 * fs_63_286_66 * h8_p3 + e_3 * fs_14_11_3 * r_2 * h6_p3 - e_3 * f_765_286 * r_4 * h4_p3 + e_4 * fs_14_4199_3003 * h10_p3 - e_4 * fs_63_2717_66 * r_2 * h8_p3 - e_4 * fs_14_187_3 * r_4 * h6_p3 + e_4 * f_17_143 * r_6 * h4_p3;
    }

    // NOTE: the angular components are formed in 87 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, ph4_m3, ph4_p4, ph6_m5, ph6_m3, ph6_p4, ph8_m5, ph8_m3, ph8_p4, ph10_m5, ph10_m3, ph10_p4, ab_2, pc_62, pc_63 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

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

        pc_62[k] = e_1 * f_15_2 * h4_p4 + e_2 * fs_105_22_5 * h6_p4 - e_2 * f_45_11 * r_2 * h4_p4 + e_3 * fs_42_143_55 * h8_p4 - e_3 * fs_14_11_5 * r_2 * h6_p4 + e_3 * f_90_143 * r_4 * h4_p4 + e_4 * fs_7_4199_3003 * h10_p4 - e_4 * fs_84_2717_55 * r_2 * h8_p4 + e_4 * fs_14_187_5 * r_4 * h6_p4 - e_4 * f_4_143 * r_6 * h4_p4;

        pc_63[k] = - e_1 * fs_15_8_42 * h4_m3 + e_2 * fs_5_44_2310 * h6_m5 - e_2 * fs_105_44_14 * h6_m3 + e_2 * fs_45_44_42 * r_2 * h4_m3 + e_3 * fs_2_143_15015 * h8_m5 - e_3 * fs_21_143_77 * h8_m3 - e_3 * fs_1_33_2310 * r_2 * h6_m5 + e_3 * fs_7_11_14 * r_2 * h6_m3 - e_3 * fs_45_286_42 * r_4 * h4_m3 + e_4 * fs_21_8398_1430 * h10_m5 - e_4 * fs_21_8398_286 * h10_m3 - e_4 * fs_4_2717_15015 * r_2 * h8_m5 + e_4 * fs_42_2717_77 * r_2 * h8_m3 + e_4 * fs_1_561_2310 * r_4 * h6_m5 - e_4 * fs_7_187_14 * r_4 * h6_m3 + e_4 * fs_1_143_42 * r_6 * h4_m3;
    }

    // NOTE: the angular components are formed in 87 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_m2, ph4_m4, ph4_m2, ph6_m4, ph6_m2, ph8_m4, ph8_m2, ph10_m4, ph10_m2, ab_2, pc_64 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

        const auto h2_m2 = ph2_m2[k];
        const auto h4_m4 = ph4_m4[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h6_m4 = ph6_m4[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h8_m4 = ph8_m4[k];
        const auto h8_m2 = ph8_m2[k];
        const auto h10_m4 = ph10_m4[k];
        const auto h10_m2 = ph10_m2[k];

        pc_64[k] = e_0 * fs_105_32_10 * h2_m2 - e_1 * fs_15_8_42 * h4_m4 + e_1 * fs_195_16_6 * h4_m2 - e_1 * fs_15_4_10 * r_2 * h2_m2 - e_2 * fs_25_44_210 * h6_m4 + e_2 * fs_35_44_7 * h6_m2 + e_2 * fs_45_44_42 * r_2 * h4_m4 - e_2 * fs_585_88_6 * r_2 * h4_m2 + e_2 * fs_5_4_10 * r_4 * h2_m2 + e_3 * fs_5_286_2310 * h8_m4 - e_3 * fs_7_22_21 * h8_m2 + e_3 * fs_5_33_210 * r_2 * h6_m4 - e_3 * fs_7_33_7 * r_2 * h6_m2 - e_3 * fs_45_286_42 * r_4 * h4_m4 + e_3 * fs_45_44_6 * r_4 * h4_m2 - e_3 * fs_5_33_10 * r_6 * h2_m2 + e_4 * fs_42_4199_286 * h10_m4 - e_4 * fs_84_4199_22 * h10_m2 - e_4 * fs_5_2717_2310 * r_2 * h8_m4 + e_4 * fs_7_209_21 * r_2 * h8_m2 - e_4 * fs_5_561_210 * r_4 * h6_m4 + e_4 * fs_7_561_7 * r_4 * h6_m2 + e_4 * fs_1_143_42 * r_6 * h4_m4 - e_4 * fs_1_22_6 * r_6 * h4_m2 + e_4 * fs_5_858_10 * r_8 * h2_m2;
    }

    // NOTE: the angular components are formed in 87 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_m1, ph4_m3, ph4_m1, ph6_m3, ph6_m1, ph8_m3, ph8_m1, ph10_m3, ph10_m1, ab_2, pc_65 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

        const auto h2_m1 = ph2_m1[k];
        const auto h4_m3 = ph4_m3[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h8_m3 = ph8_m3[k];
        const auto h8_m1 = ph8_m1[k];
        const auto h10_m3 = ph10_m3[k];
        const auto h10_m1 = ph10_m1[k];

        pc_65[k] = - e_0 * fs_105_16_35 * h2_m1 + e_1 * fs_195_16_6 * h4_m3 - e_1 * fs_45_16_42 * h4_m1 + e_1 * fs_15_2_35 * r_2 * h2_m1 - e_2 * fs_15_11_2 * h6_m3 + e_2 * fs_35_11_5 * h6_m1 - e_2 * fs_585_88_6 * r_2 * h4_m3 + e_2 * fs_135_88_42 * r_2 * h4_m1 - e_2 * fs_5_2_35 * r_4 * h2_m1 - e_3 * fs_3_22_11 * h8_m3 - e_3 * fs_23_286_105 * h8_m1 + e_3 * fs_4_11_2 * r_2 * h6_m3 - e_3 * fs_28_33_5 * r_2 * h6_m1 + e_3 * fs_45_44_6 * r_4 * h4_m3 - e_3 * fs_135_572_42 * r_4 * h4_m1 + e_3 * fs_10_33_35 * r_6 * h2_m1 + e_4 * fs_21_4199_2002 * h10_m3 - e_4 * fs_42_4199_231 * h10_m1 + e_4 * fs_3_209_11 * r_2 * h8_m3 + e_4 * fs_23_2717_105 * r_2 * h8_m1 - e_4 * fs_4_187_2 * r_4 * h6_m3 + e_4 * fs_28_561_5 * r_4 * h6_m1 - e_4 * fs_1_22_6 * r_6 * h4_m3 + e_4 * fs_3_286_42 * r_6 * h4_m1 - e_4 * fs_5_429_35 * r_8 * h2_m1;
    }

    // NOTE: the angular components are formed in 87 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_m2, ph2_p1, ph4_m2, ph4_p1, ph6_m2, ph6_p1, ph8_m2, ph8_p1, ph10_m2, ph10_p1, ab_2, pc_66, pc_67 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

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

        pc_66[k] = - e_0 * fs_105_32_70 * h2_m2 - e_1 * fs_45_16_42 * h4_m2 + e_1 * fs_15_4_70 * r_2 * h2_m2 + e_2 * f_265_44 * h6_m2 + e_2 * fs_135_88_42 * r_2 * h4_m2 - e_2 * fs_5_4_70 * r_4 * h2_m2 - e_3 * fs_199_286_3 * h8_m2 - e_3 * f_53_33 * r_2 * h6_m2 - e_3 * fs_135_572_42 * r_4 * h4_m2 + e_3 * fs_5_33_70 * r_6 * h2_m2 + e_4 * fs_84_4199_154 * h10_m2 + e_4 * fs_199_2717_3 * r_2 * h8_m2 + e_4 * f_53_561 * r_4 * h6_m2 + e_4 * fs_3_286_42 * r_6 * h4_m2 - e_4 * fs_5_858_70 * r_8 * h2_m2;

        pc_67[k] = e_0 * fs_525_16_7 * h2_p1 - e_1 * fs_15_8_210 * h4_p1 - e_1 * fs_75_2_7 * r_2 * h2_p1 + e_2 * f_80_11 * h6_p1 + e_2 * fs_45_44_210 * r_2 * h4_p1 + e_2 * fs_25_2_7 * r_4 * h2_p1 - e_3 * fs_47_143_21 * h8_p1 - e_3 * f_64_33 * r_2 * h6_p1 - e_3 * fs_45_286_210 * r_4 * h4_p1 - e_3 * fs_50_33_7 * r_6 * h2_p1 + e_4 * fs_42_4199_1155 * h10_p1 + e_4 * fs_94_2717_21 * r_2 * h8_p1 + e_4 * f_64_561 * r_4 * h6_p1 + e_4 * fs_1_143_210 * r_6 * h4_p1 + e_4 * fs_25_429_7 * r_8 * h2_p1;
    }

    // NOTE: the angular components are formed in 87 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_0, ph2_p2, ph4_0, ph4_p2, ph6_0, ph6_p2, ph8_0, ph8_p2, ph10_0, ph10_p2, ab_2, pc_68 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

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

        pc_68[k] = e_0 * fs_105_16_210 * h2_0 - e_0 * fs_105_32_70 * h2_p2 - e_1 * fs_15_8_210 * h4_0 - e_1 * fs_45_16_42 * h4_p2 - e_1 * fs_15_2_210 * r_2 * h2_0 + e_1 * fs_15_4_70 * r_2 * h2_p2 + e_2 * fs_5_22_210 * h6_0 + e_2 * f_265_44 * h6_p2 + e_2 * fs_45_44_210 * r_2 * h4_0 + e_2 * fs_135_88_42 * r_2 * h4_p2 + e_2 * fs_5_2_210 * r_4 * h2_0 - e_2 * fs_5_4_70 * r_4 * h2_p2 + e_3 * fs_3_143_210 * h8_0 - e_3 * fs_199_286_3 * h8_p2 - e_3 * fs_2_33_210 * r_2 * h6_0 - e_3 * f_53_33 * r_2 * h6_p2 - e_3 * fs_45_286_210 * r_4 * h4_0 - e_3 * fs_135_572_42 * r_4 * h4_p2 - e_3 * fs_10_33_210 * r_6 * h2_0 + e_3 * fs_5_33_70 * r_6 * h2_p2 - e_4 * fs_84_4199_210 * h10_0 + e_4 * fs_84_4199_154 * h10_p2 - e_4 * fs_6_2717_210 * r_2 * h8_0 + e_4 * fs_199_2717_3 * r_2 * h8_p2 + e_4 * fs_2_561_210 * r_4 * h6_0 + e_4 * f_53_561 * r_4 * h6_p2 + e_4 * fs_1_143_210 * r_6 * h4_0 + e_4 * fs_3_286_42 * r_6 * h4_p2 + e_4 * fs_5_429_210 * r_8 * h2_0 - e_4 * fs_5_858_70 * r_8 * h2_p2;
    }

    // NOTE: the angular components are formed in 87 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_p1, ph4_p1, ph4_p3, ph6_p1, ph6_p3, ph8_p1, ph8_p3, ph10_p1, ph10_p3, ab_2, pc_69 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

        const auto h2_p1 = ph2_p1[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h4_p3 = ph4_p3[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h8_p3 = ph8_p3[k];
        const auto h10_p1 = ph10_p1[k];
        const auto h10_p3 = ph10_p3[k];

        pc_69[k] = - e_0 * fs_105_16_35 * h2_p1 - e_1 * fs_45_16_42 * h4_p1 + e_1 * fs_195_16_6 * h4_p3 + e_1 * fs_15_2_35 * r_2 * h2_p1 + e_2 * fs_35_11_5 * h6_p1 - e_2 * fs_15_11_2 * h6_p3 + e_2 * fs_135_88_42 * r_2 * h4_p1 - e_2 * fs_585_88_6 * r_2 * h4_p3 - e_2 * fs_5_2_35 * r_4 * h2_p1 - e_3 * fs_23_286_105 * h8_p1 - e_3 * fs_3_22_11 * h8_p3 - e_3 * fs_28_33_5 * r_2 * h6_p1 + e_3 * fs_4_11_2 * r_2 * h6_p3 - e_3 * fs_135_572_42 * r_4 * h4_p1 + e_3 * fs_45_44_6 * r_4 * h4_p3 + e_3 * fs_10_33_35 * r_6 * h2_p1 - e_4 * fs_42_4199_231 * h10_p1 + e_4 * fs_21_4199_2002 * h10_p3 + e_4 * fs_23_2717_105 * r_2 * h8_p1 + e_4 * fs_3_209_11 * r_2 * h8_p3 + e_4 * fs_28_561_5 * r_4 * h6_p1 - e_4 * fs_4_187_2 * r_4 * h6_p3 + e_4 * fs_3_286_42 * r_6 * h4_p1 - e_4 * fs_1_22_6 * r_6 * h4_p3 - e_4 * fs_5_429_35 * r_8 * h2_p1;
    }

    // NOTE: the angular components are formed in 87 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_p2, ph4_p2, ph4_p4, ph6_p2, ph6_p4, ph8_p2, ph8_p4, ph10_p2, ph10_p4, ab_2, pc_70 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

        const auto h2_p2 = ph2_p2[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h4_p4 = ph4_p4[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h6_p4 = ph6_p4[k];
        const auto h8_p2 = ph8_p2[k];
        const auto h8_p4 = ph8_p4[k];
        const auto h10_p2 = ph10_p2[k];
        const auto h10_p4 = ph10_p4[k];

        pc_70[k] = e_0 * fs_105_32_10 * h2_p2 + e_1 * fs_195_16_6 * h4_p2 - e_1 * fs_15_8_42 * h4_p4 - e_1 * fs_15_4_10 * r_2 * h2_p2 + e_2 * fs_35_44_7 * h6_p2 - e_2 * fs_25_44_210 * h6_p4 - e_2 * fs_585_88_6 * r_2 * h4_p2 + e_2 * fs_45_44_42 * r_2 * h4_p4 + e_2 * fs_5_4_10 * r_4 * h2_p2 - e_3 * fs_7_22_21 * h8_p2 + e_3 * fs_5_286_2310 * h8_p4 - e_3 * fs_7_33_7 * r_2 * h6_p2 + e_3 * fs_5_33_210 * r_2 * h6_p4 + e_3 * fs_45_44_6 * r_4 * h4_p2 - e_3 * fs_45_286_42 * r_4 * h4_p4 - e_3 * fs_5_33_10 * r_6 * h2_p2 - e_4 * fs_84_4199_22 * h10_p2 + e_4 * fs_42_4199_286 * h10_p4 + e_4 * fs_7_209_21 * r_2 * h8_p2 - e_4 * fs_5_2717_2310 * r_2 * h8_p4 + e_4 * fs_7_561_7 * r_4 * h6_p2 - e_4 * fs_5_561_210 * r_4 * h6_p4 - e_4 * fs_1_22_6 * r_6 * h4_p2 + e_4 * fs_1_143_42 * r_6 * h4_p4 + e_4 * fs_5_858_10 * r_8 * h2_p2;
    }

    // NOTE: the angular components are formed in 87 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, ph4_p3, ph6_p3, ph6_p5, ph8_p3, ph8_p5, ph10_p3, ph10_p5, ab_2, pc_71 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h4_p3 = ph4_p3[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h6_p5 = ph6_p5[k];
        const auto h8_p3 = ph8_p3[k];
        const auto h8_p5 = ph8_p5[k];
        const auto h10_p3 = ph10_p3[k];
        const auto h10_p5 = ph10_p5[k];

        pc_71[k] = - e_1 * fs_15_8_42 * h4_p3 - e_2 * fs_105_44_14 * h6_p3 + e_2 * fs_5_44_2310 * h6_p5 + e_2 * fs_45_44_42 * r_2 * h4_p3 - e_3 * fs_21_143_77 * h8_p3 + e_3 * fs_2_143_15015 * h8_p5 + e_3 * fs_7_11_14 * r_2 * h6_p3 - e_3 * fs_1_33_2310 * r_2 * h6_p5 - e_3 * fs_45_286_42 * r_4 * h4_p3 - e_4 * fs_21_8398_286 * h10_p3 + e_4 * fs_21_8398_1430 * h10_p5 + e_4 * fs_42_2717_77 * r_2 * h8_p3 - e_4 * fs_4_2717_15015 * r_2 * h8_p5 - e_4 * fs_7_187_14 * r_4 * h6_p3 + e_4 * fs_1_561_2310 * r_4 * h6_p5 + e_4 * fs_1_143_42 * r_6 * h4_p3;
    }

    // NOTE: the angular components are formed in 87 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_m2, ph4_m2, ph6_m6, ph6_m2, ph8_m6, ph8_m2, ph10_m6, ph10_m2, ab_2, pc_72 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

        const auto h2_m2 = ph2_m2[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h6_m6 = ph6_m6[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h8_m6 = ph8_m6[k];
        const auto h8_m2 = ph8_m2[k];
        const auto h10_m6 = ph10_m6[k];
        const auto h10_m2 = ph10_m2[k];

        pc_72[k] = e_0 * fs_105_32_2 * h2_m2 + e_1 * fs_15_4_30 * h4_m2 - e_1 * fs_15_4_2 * r_2 * h2_m2 + e_2 * fs_15_44_77 * h6_m6 + e_2 * fs_35_22_35 * h6_m2 - e_2 * fs_45_22_30 * r_2 * h4_m2 + e_2 * fs_5_4_2 * r_4 * h2_m2 + e_3 * fs_21_143_143 * h8_m6 + e_3 * fs_14_143_105 * h8_m2 - e_3 * fs_1_11_77 * r_2 * h6_m6 - e_3 * fs_14_33_35 * r_2 * h6_m2 + e_3 * fs_45_143_30 * r_4 * h4_m2 - e_3 * fs_5_33_2 * r_6 * h2_m2 + e_4 * fs_21_4199_715 * h10_m6 + e_4 * fs_21_8398_110 * h10_m2 - e_4 * fs_42_2717_143 * r_2 * h8_m6 - e_4 * fs_28_2717_105 * r_2 * h8_m2 + e_4 * fs_1_187_77 * r_4 * h6_m6 + e_4 * fs_14_561_35 * r_4 * h6_m2 - e_4 * fs_2_143_30 * r_6 * h4_m2 + e_4 * fs_5_858_2 * r_8 * h2_m2;
    }

    // NOTE: the angular components are formed in 87 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_m1, ph4_m4, ph4_m1, ph6_m5, ph6_m4, ph6_m1, ph8_m5, ph8_m4, ph8_m1, ph10_m5, ph10_m4, ph10_m1, ab_2, pc_73, pc_74 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

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

        pc_73[k] = - e_0 * f_105_4 * h2_m1 - e_1 * fs_45_8_30 * h4_m1 + e_1 * f_30_1 * r_2 * h2_m1 - e_2 * fs_35_88_462 * h6_m5 + e_2 * fs_35_44_7 * h6_m1 + e_2 * fs_135_44_30 * r_2 * h4_m1 - e_2 * f_10_1 * r_4 * h2_m1 + e_3 * fs_1_286_3003 * h8_m5 + e_3 * fs_119_143_3 * h8_m1 + e_3 * fs_7_66_462 * r_2 * h6_m5 - e_3 * fs_7_33_7 * r_2 * h6_m1 - e_3 * fs_135_286_30 * r_4 * h4_m1 + e_3 * f_40_33 * r_6 * h2_m1 + e_4 * fs_105_8398_286 * h10_m5 + e_4 * fs_21_4199_165 * h10_m1 - e_4 * fs_1_2717_3003 * r_2 * h8_m5 - e_4 * fs_238_2717_3 * r_2 * h8_m1 - e_4 * fs_7_1122_462 * r_4 * h6_m5 + e_4 * fs_7_561_7 * r_4 * h6_m1 + e_4 * fs_3_143_30 * r_6 * h4_m1 - e_4 * f_20_429 * r_8 * h2_m1;

        pc_74[k] = e_1 * fs_15_4_30 * h4_m4 + e_2 * fs_115_88_6 * h6_m4 - e_2 * fs_45_22_30 * r_2 * h4_m4 - e_3 * fs_37_286_66 * h8_m4 - e_3 * fs_23_66_6 * r_2 * h6_m4 + e_3 * fs_45_143_30 * r_4 * h4_m4 + e_4 * fs_21_8398_10010 * h10_m4 + e_4 * fs_37_2717_66 * r_2 * h8_m4 + e_4 * fs_23_1122_6 * r_4 * h6_m4 - e_4 * fs_2_143_30 * r_6 * h4_m4;
    }

    // NOTE: the angular components are formed in 87 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_m1, ph4_m3, ph4_m1, ph6_m3, ph6_m1, ph8_m3, ph8_m1, ph10_m3, ph10_m1, ab_2, pc_75 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

        const auto h2_m1 = ph2_m1[k];
        const auto h4_m3 = ph4_m3[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h8_m3 = ph8_m3[k];
        const auto h8_m1 = ph8_m1[k];
        const auto h10_m3 = ph10_m3[k];
        const auto h10_m1 = ph10_m1[k];

        pc_75[k] = - e_0 * fs_105_4_7 * h2_m1 - e_1 * fs_45_8_30 * h4_m3 + e_1 * fs_15_8_210 * h4_m1 + e_1 * fs_30_1_7 * r_2 * h2_m1 + e_2 * fs_195_88_10 * h6_m3 - e_2 * f_265_44 * h6_m1 + e_2 * fs_135_44_30 * r_2 * h4_m3 - e_2 * fs_45_44_210 * r_2 * h4_m1 - e_2 * fs_10_1_7 * r_4 * h2_m1 - e_3 * fs_24_143_55 * h8_m3 + e_3 * fs_29_286_21 * h8_m1 - e_3 * fs_13_22_10 * r_2 * h6_m3 + e_3 * f_53_33 * r_2 * h6_m1 - e_3 * fs_135_286_30 * r_4 * h4_m3 + e_3 * fs_45_286_210 * r_4 * h4_m1 + e_3 * fs_40_33_7 * r_6 * h2_m1 + e_4 * fs_21_8398_10010 * h10_m3 + e_4 * fs_21_4199_1155 * h10_m1 + e_4 * fs_48_2717_55 * r_2 * h8_m3 - e_4 * fs_29_2717_21 * r_2 * h8_m1 + e_4 * fs_13_374_10 * r_4 * h6_m3 - e_4 * f_53_561 * r_4 * h6_m1 + e_4 * fs_3_143_30 * r_6 * h4_m3 - e_4 * fs_1_143_210 * r_6 * h4_m1 - e_4 * fs_20_429_7 * r_8 * h2_m1;
    }

    // NOTE: the angular components are formed in 87 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_p2, ph6_p2, ph8_p2, ph10_p2, ab_2, pc_76 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

        const auto h2_p2 = ph2_p2[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h8_p2 = ph8_p2[k];
        const auto h10_p2 = ph10_p2[k];

        pc_76[k] = e_0 * fs_105_16_70 * h2_p2 - e_1 * fs_15_2_70 * r_2 * h2_p2 + e_2 * f_5_4 * h6_p2 + e_2 * fs_5_2_70 * r_4 * h2_p2 - e_3 * fs_5_11_3 * h8_p2 - e_3 * f_1_3 * r_2 * h6_p2 - e_3 * fs_10_33_70 * r_6 * h2_p2 + e_4 * fs_105_4199_154 * h10_p2 + e_4 * fs_10_209_3 * r_2 * h8_p2 + e_4 * f_1_51 * r_4 * h6_p2 + e_4 * fs_5_429_70 * r_8 * h2_p2;
    }

    // NOTE: the angular components are formed in 87 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_p1, ph4_p1, ph4_p3, ph6_p1, ph6_p3, ph8_p1, ph8_p3, ph10_p1, ph10_p3, ab_2, pc_77 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

        const auto h2_p1 = ph2_p1[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h4_p3 = ph4_p3[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h8_p3 = ph8_p3[k];
        const auto h10_p1 = ph10_p1[k];
        const auto h10_p3 = ph10_p3[k];

        pc_77[k] = e_0 * fs_105_4_7 * h2_p1 - e_1 * fs_15_8_210 * h4_p1 - e_1 * fs_45_8_30 * h4_p3 - e_1 * fs_30_1_7 * r_2 * h2_p1 + e_2 * f_265_44 * h6_p1 + e_2 * fs_195_88_10 * h6_p3 + e_2 * fs_45_44_210 * r_2 * h4_p1 + e_2 * fs_135_44_30 * r_2 * h4_p3 + e_2 * fs_10_1_7 * r_4 * h2_p1 - e_3 * fs_29_286_21 * h8_p1 - e_3 * fs_24_143_55 * h8_p3 - e_3 * f_53_33 * r_2 * h6_p1 - e_3 * fs_13_22_10 * r_2 * h6_p3 - e_3 * fs_45_286_210 * r_4 * h4_p1 - e_3 * fs_135_286_30 * r_4 * h4_p3 - e_3 * fs_40_33_7 * r_6 * h2_p1 - e_4 * fs_21_4199_1155 * h10_p1 + e_4 * fs_21_8398_10010 * h10_p3 + e_4 * fs_29_2717_21 * r_2 * h8_p1 + e_4 * fs_48_2717_55 * r_2 * h8_p3 + e_4 * f_53_561 * r_4 * h6_p1 + e_4 * fs_13_374_10 * r_4 * h6_p3 + e_4 * fs_1_143_210 * r_6 * h4_p1 + e_4 * fs_3_143_30 * r_6 * h4_p3 + e_4 * fs_20_429_7 * r_8 * h2_p1;
    }

    // NOTE: the angular components are formed in 87 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_0, ph4_p4, ph6_0, ph6_p4, ph8_0, ph8_p4, ph10_0, ph10_p4, ab_2, pc_78 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

        const auto h2_0 = ph2_0[k];
        const auto h4_p4 = ph4_p4[k];
        const auto h6_0 = ph6_0[k];
        const auto h6_p4 = ph6_p4[k];
        const auto h8_0 = ph8_0[k];
        const auto h8_p4 = ph8_p4[k];
        const auto h10_0 = ph10_0[k];
        const auto h10_p4 = ph10_p4[k];

        pc_78[k] = e_0 * fs_105_8_42 * h2_0 + e_1 * fs_15_4_30 * h4_p4 - e_1 * fs_15_1_42 * r_2 * h2_0 - e_2 * fs_5_4_42 * h6_0 + e_2 * fs_115_88_6 * h6_p4 - e_2 * fs_45_22_30 * r_2 * h4_p4 + e_2 * fs_5_1_42 * r_4 * h2_0 + e_3 * fs_3_11_42 * h8_0 - e_3 * fs_37_286_66 * h8_p4 + e_3 * fs_1_3_42 * r_2 * h6_0 - e_3 * fs_23_66_6 * r_2 * h6_p4 + e_3 * fs_45_143_30 * r_4 * h4_p4 - e_3 * fs_20_33_42 * r_6 * h2_0 + e_4 * fs_105_4199_42 * h10_0 + e_4 * fs_21_8398_10010 * h10_p4 - e_4 * fs_6_209_42 * r_2 * h8_0 + e_4 * fs_37_2717_66 * r_2 * h8_p4 - e_4 * fs_1_51_42 * r_4 * h6_0 + e_4 * fs_23_1122_6 * r_4 * h6_p4 - e_4 * fs_2_143_30 * r_6 * h4_p4 + e_4 * fs_10_429_42 * r_8 * h2_0;
    }

    // NOTE: the angular components are formed in 87 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_p1, ph4_p1, ph6_p1, ph6_p5, ph8_p1, ph8_p5, ph10_p1, ph10_p5, ab_2, pc_79 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

        const auto h2_p1 = ph2_p1[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h6_p5 = ph6_p5[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h8_p5 = ph8_p5[k];
        const auto h10_p1 = ph10_p1[k];
        const auto h10_p5 = ph10_p5[k];

        pc_79[k] = - e_0 * f_105_4 * h2_p1 - e_1 * fs_45_8_30 * h4_p1 + e_1 * f_30_1 * r_2 * h2_p1 + e_2 * fs_35_44_7 * h6_p1 - e_2 * fs_35_88_462 * h6_p5 + e_2 * fs_135_44_30 * r_2 * h4_p1 - e_2 * f_10_1 * r_4 * h2_p1 + e_3 * fs_119_143_3 * h8_p1 + e_3 * fs_1_286_3003 * h8_p5 - e_3 * fs_7_33_7 * r_2 * h6_p1 + e_3 * fs_7_66_462 * r_2 * h6_p5 - e_3 * fs_135_286_30 * r_4 * h4_p1 + e_3 * f_40_33 * r_6 * h2_p1 + e_4 * fs_21_4199_165 * h10_p1 + e_4 * fs_105_8398_286 * h10_p5 - e_4 * fs_238_2717_3 * r_2 * h8_p1 - e_4 * fs_1_2717_3003 * r_2 * h8_p5 + e_4 * fs_7_561_7 * r_4 * h6_p1 - e_4 * fs_7_1122_462 * r_4 * h6_p5 + e_4 * fs_3_143_30 * r_6 * h4_p1 - e_4 * f_20_429 * r_8 * h2_p1;
    }

    // NOTE: the angular components are formed in 87 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_p2, ph4_p2, ph6_p2, ph6_p6, ph8_p2, ph8_p6, ph10_p2, ph10_p6, ab_2, pc_80 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

        const auto h2_p2 = ph2_p2[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h6_p6 = ph6_p6[k];
        const auto h8_p2 = ph8_p2[k];
        const auto h8_p6 = ph8_p6[k];
        const auto h10_p2 = ph10_p2[k];
        const auto h10_p6 = ph10_p6[k];

        pc_80[k] = e_0 * fs_105_32_2 * h2_p2 + e_1 * fs_15_4_30 * h4_p2 - e_1 * fs_15_4_2 * r_2 * h2_p2 + e_2 * fs_35_22_35 * h6_p2 + e_2 * fs_15_44_77 * h6_p6 - e_2 * fs_45_22_30 * r_2 * h4_p2 + e_2 * fs_5_4_2 * r_4 * h2_p2 + e_3 * fs_14_143_105 * h8_p2 + e_3 * fs_21_143_143 * h8_p6 - e_3 * fs_14_33_35 * r_2 * h6_p2 - e_3 * fs_1_11_77 * r_2 * h6_p6 + e_3 * fs_45_143_30 * r_4 * h4_p2 - e_3 * fs_5_33_2 * r_6 * h2_p2 + e_4 * fs_21_8398_110 * h10_p2 + e_4 * fs_21_4199_715 * h10_p6 - e_4 * fs_28_2717_105 * r_2 * h8_p2 - e_4 * fs_42_2717_143 * r_2 * h8_p6 + e_4 * fs_14_561_35 * r_4 * h6_p2 + e_4 * fs_1_187_77 * r_4 * h6_p6 - e_4 * fs_2_143_30 * r_6 * h4_p2 + e_4 * fs_5_858_2 * r_8 * h2_p2;
    }

    // NOTE: the angular components are formed in 87 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_m1, ph4_m1, ph6_m6, ph6_m1, ph8_m7, ph8_m6, ph8_m1, ph10_m7, ph10_m6, ph10_m1, ab_2, pc_81, pc_82 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

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

        pc_81[k] = - e_0 * fs_315_32_2 * h2_m1 - e_1 * fs_15_2_15 * h4_m1 + e_1 * fs_45_4_2 * r_2 * h2_m1 - e_2 * fs_105_44_14 * h6_m1 + e_2 * fs_45_11_15 * r_2 * h4_m1 - e_2 * fs_15_4_2 * r_4 * h2_m1 + e_3 * fs_7_286_4290 * h8_m7 - e_3 * fs_42_143_6 * h8_m1 + e_3 * fs_7_11_14 * r_2 * h6_m1 - e_3 * fs_90_143_15 * r_4 * h4_m1 + e_3 * fs_5_11_2 * r_6 * h2_m1 + e_4 * fs_7_4199_12155 * h10_m7 - e_4 * fs_7_8398_330 * h10_m1 - e_4 * fs_7_2717_4290 * r_2 * h8_m7 + e_4 * fs_84_2717_6 * r_2 * h8_m1 - e_4 * fs_7_187_14 * r_4 * h6_m1 + e_4 * fs_4_143_15 * r_6 * h4_m1 - e_4 * fs_5_286_2 * r_8 * h2_m1;

        pc_82[k] = - e_2 * fs_45_88_154 * h6_m6 - e_3 * fs_21_572_286 * h8_m6 + e_3 * fs_3_22_154 * r_2 * h6_m6 + e_4 * fs_28_4199_1430 * h10_m6 + e_4 * fs_21_5434_286 * r_2 * h8_m6 - e_4 * fs_3_374_154 * r_4 * h6_m6;
    }

    // NOTE: the angular components are formed in 87 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_m1, ph4_m1, ph6_m5, ph6_m1, ph8_m5, ph8_m1, ph10_m5, ph10_m1, ab_2, pc_83 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

        const auto h2_m1 = ph2_m1[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h6_m5 = ph6_m5[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h8_m5 = ph8_m5[k];
        const auto h8_m1 = ph8_m1[k];
        const auto h10_m5 = ph10_m5[k];
        const auto h10_m1 = ph10_m1[k];

        pc_83[k] = - e_0 * fs_315_16_14 * h2_m1 + e_1 * fs_15_8_105 * h4_m1 + e_1 * fs_45_2_14 * r_2 * h2_m1 + e_2 * fs_15_11_33 * h6_m5 + e_2 * fs_15_11_2 * h6_m1 - e_2 * fs_45_44_105 * r_2 * h4_m1 - e_2 * fs_15_2_14 * r_4 * h2_m1 - e_3 * fs_27_572_858 * h8_m5 - e_3 * fs_127_572_42 * h8_m1 - e_3 * fs_4_11_33 * r_2 * h6_m5 - e_3 * fs_4_11_2 * r_2 * h6_m1 + e_3 * fs_45_286_105 * r_4 * h4_m1 + e_3 * fs_10_11_14 * r_6 * h2_m1 + e_4 * fs_35_4199_1001 * h10_m5 - e_4 * fs_7_4199_2310 * h10_m1 + e_4 * fs_27_5434_858 * r_2 * h8_m5 + e_4 * fs_127_5434_42 * r_2 * h8_m1 + e_4 * fs_4_187_33 * r_4 * h6_m5 + e_4 * fs_4_187_2 * r_4 * h6_m1 - e_4 * fs_1_143_105 * r_6 * h4_m1 - e_4 * fs_5_143_14 * r_8 * h2_m1;
    }

    // NOTE: the angular components are formed in 87 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_m2, ph4_m4, ph4_m2, ph6_m4, ph6_m2, ph8_m4, ph8_m2, ph10_m4, ph10_m2, ab_2, pc_84 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

        const auto h2_m2 = ph2_m2[k];
        const auto h4_m4 = ph4_m4[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h6_m4 = ph6_m4[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h8_m4 = ph8_m4[k];
        const auto h8_m2 = ph8_m2[k];
        const auto h10_m4 = ph10_m4[k];
        const auto h10_m2 = ph10_m2[k];

        pc_84[k] = - e_0 * fs_315_16_7 * h2_m2 - e_1 * fs_15_2_15 * h4_m4 + e_1 * fs_15_8_105 * h4_m2 + e_1 * fs_45_2_7 * r_2 * h2_m2 + e_2 * fs_105_44_3 * h6_m4 - e_2 * fs_195_88_10 * h6_m2 + e_2 * fs_45_11_15 * r_2 * h4_m4 - e_2 * fs_45_44_105 * r_2 * h4_m2 - e_2 * fs_15_2_7 * r_4 * h2_m2 - e_3 * fs_43_286_33 * h8_m4 + e_3 * fs_111_572_30 * h8_m2 - e_3 * fs_7_11_3 * r_2 * h6_m4 + e_3 * fs_13_22_10 * r_2 * h6_m2 - e_3 * fs_90_143_15 * r_4 * h4_m4 + e_3 * fs_45_286_105 * r_4 * h4_m2 + e_3 * fs_10_11_7 * r_6 * h2_m2 + e_4 * fs_14_4199_5005 * h10_m4 + e_4 * fs_28_4199_385 * h10_m2 + e_4 * fs_43_2717_33 * r_2 * h8_m4 - e_4 * fs_111_5434_30 * r_2 * h8_m2 + e_4 * fs_7_187_3 * r_4 * h6_m4 - e_4 * fs_13_374_10 * r_4 * h6_m2 + e_4 * fs_4_143_15 * r_6 * h4_m4 - e_4 * fs_1_143_105 * r_6 * h4_m2 - e_4 * fs_5_143_7 * r_8 * h2_m2;
    }

    // NOTE: the angular components are formed in 87 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, ph4_p3, ph6_p3, ph8_p3, ph10_p3, ab_2, pc_85 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h4_p3 = ph4_p3[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h8_p3 = ph8_p3[k];
        const auto h10_p3 = ph10_p3[k];

        pc_85[k] = e_1 * fs_75_4_3 * h4_p3 - e_2 * f_135_22 * h6_p3 - e_2 * fs_225_22_3 * r_2 * h4_p3 + e_3 * fs_15_286_22 * h8_p3 + e_3 * f_18_11 * r_2 * h6_p3 + e_3 * fs_225_143_3 * r_4 * h4_p3 + e_4 * fs_35_4199_1001 * h10_p3 - e_4 * fs_15_2717_22 * r_2 * h8_p3 - e_4 * f_18_187 * r_4 * h6_p3 - e_4 * fs_10_143_3 * r_6 * h4_p3;
    }

    // NOTE: the angular components are formed in 87 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_p2, ph4_p2, ph4_p4, ph6_p2, ph6_p4, ph8_p2, ph8_p4, ph10_p2, ph10_p4, ab_2, pc_86 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

        const auto h2_p2 = ph2_p2[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h4_p4 = ph4_p4[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h6_p4 = ph6_p4[k];
        const auto h8_p2 = ph8_p2[k];
        const auto h8_p4 = ph8_p4[k];
        const auto h10_p2 = ph10_p2[k];
        const auto h10_p4 = ph10_p4[k];

        pc_86[k] = e_0 * fs_315_16_7 * h2_p2 - e_1 * fs_15_8_105 * h4_p2 - e_1 * fs_15_2_15 * h4_p4 - e_1 * fs_45_2_7 * r_2 * h2_p2 + e_2 * fs_195_88_10 * h6_p2 + e_2 * fs_105_44_3 * h6_p4 + e_2 * fs_45_44_105 * r_2 * h4_p2 + e_2 * fs_45_11_15 * r_2 * h4_p4 + e_2 * fs_15_2_7 * r_4 * h2_p2 - e_3 * fs_111_572_30 * h8_p2 - e_3 * fs_43_286_33 * h8_p4 - e_3 * fs_13_22_10 * r_2 * h6_p2 - e_3 * fs_7_11_3 * r_2 * h6_p4 - e_3 * fs_45_286_105 * r_4 * h4_p2 - e_3 * fs_90_143_15 * r_4 * h4_p4 - e_3 * fs_10_11_7 * r_6 * h2_p2 - e_4 * fs_28_4199_385 * h10_p2 + e_4 * fs_14_4199_5005 * h10_p4 + e_4 * fs_111_5434_30 * r_2 * h8_p2 + e_4 * fs_43_2717_33 * r_2 * h8_p4 + e_4 * fs_13_374_10 * r_4 * h6_p2 + e_4 * fs_7_187_3 * r_4 * h6_p4 + e_4 * fs_1_143_105 * r_6 * h4_p2 + e_4 * fs_4_143_15 * r_6 * h4_p4 + e_4 * fs_5_143_7 * r_8 * h2_p2;
    }

    // NOTE: the angular components are formed in 87 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_p1, ph4_p1, ph6_p1, ph6_p5, ph8_p1, ph8_p5, ph10_p1, ph10_p5, ab_2, pc_87 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

        const auto h2_p1 = ph2_p1[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h6_p5 = ph6_p5[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h8_p5 = ph8_p5[k];
        const auto h10_p1 = ph10_p1[k];
        const auto h10_p5 = ph10_p5[k];

        pc_87[k] = e_0 * fs_315_16_14 * h2_p1 - e_1 * fs_15_8_105 * h4_p1 - e_1 * fs_45_2_14 * r_2 * h2_p1 - e_2 * fs_15_11_2 * h6_p1 + e_2 * fs_15_11_33 * h6_p5 + e_2 * fs_45_44_105 * r_2 * h4_p1 + e_2 * fs_15_2_14 * r_4 * h2_p1 + e_3 * fs_127_572_42 * h8_p1 - e_3 * fs_27_572_858 * h8_p5 + e_3 * fs_4_11_2 * r_2 * h6_p1 - e_3 * fs_4_11_33 * r_2 * h6_p5 - e_3 * fs_45_286_105 * r_4 * h4_p1 - e_3 * fs_10_11_14 * r_6 * h2_p1 + e_4 * fs_7_4199_2310 * h10_p1 + e_4 * fs_35_4199_1001 * h10_p5 - e_4 * fs_127_5434_42 * r_2 * h8_p1 + e_4 * fs_27_5434_858 * r_2 * h8_p5 - e_4 * fs_4_187_2 * r_4 * h6_p1 + e_4 * fs_4_187_33 * r_4 * h6_p5 + e_4 * fs_1_143_105 * r_6 * h4_p1 + e_4 * fs_5_143_14 * r_8 * h2_p1;
    }

    // NOTE: the angular components are formed in 87 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_0, ph4_0, ph6_0, ph6_p6, ph8_0, ph8_p6, ph10_0, ph10_p6, ab_2, pc_88 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

        const auto h2_0 = ph2_0[k];
        const auto h4_0 = ph4_0[k];
        const auto h6_0 = ph6_0[k];
        const auto h6_p6 = ph6_p6[k];
        const auto h8_0 = ph8_0[k];
        const auto h8_p6 = ph8_p6[k];
        const auto h10_0 = ph10_0[k];
        const auto h10_p6 = ph10_p6[k];

        pc_88[k] = e_0 * fs_315_8_3 * h2_0 + e_1 * fs_75_4_3 * h4_0 - e_1 * fs_45_1_3 * r_2 * h2_0 - e_2 * fs_105_22_3 * h6_0 - e_2 * fs_45_88_154 * h6_p6 - e_2 * fs_225_22_3 * r_2 * h4_0 + e_2 * fs_15_1_3 * r_4 * h2_0 - e_3 * fs_147_143_3 * h8_0 - e_3 * fs_21_572_286 * h8_p6 + e_3 * fs_14_11_3 * r_2 * h6_0 + e_3 * fs_3_22_154 * r_2 * h6_p6 + e_3 * fs_225_143_3 * r_4 * h4_0 - e_3 * fs_20_11_3 * r_6 * h2_0 - e_4 * fs_140_4199_3 * h10_0 + e_4 * fs_28_4199_1430 * h10_p6 + e_4 * fs_294_2717_3 * r_2 * h8_0 + e_4 * fs_21_5434_286 * r_2 * h8_p6 - e_4 * fs_14_187_3 * r_4 * h6_0 - e_4 * fs_3_374_154 * r_4 * h6_p6 - e_4 * fs_10_143_3 * r_6 * h4_0 + e_4 * fs_10_143_3 * r_8 * h2_0;
    }

    // NOTE: the angular components are formed in 87 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_p1, ph4_p1, ph6_p1, ph8_m8, ph8_p1, ph8_p7, ph10_m8, ph10_p1, ph10_p7, ab_2, pc_89, pc_90 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

        const auto h2_p1 = ph2_p1[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h8_m8 = ph8_m8[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h8_p7 = ph8_p7[k];
        const auto h10_m8 = ph10_m8[k];
        const auto h10_p1 = ph10_p1[k];
        const auto h10_p7 = ph10_p7[k];

        pc_89[k] = - e_0 * fs_315_32_2 * h2_p1 - e_1 * fs_15_2_15 * h4_p1 + e_1 * fs_45_4_2 * r_2 * h2_p1 - e_2 * fs_105_44_14 * h6_p1 + e_2 * fs_45_11_15 * r_2 * h4_p1 - e_2 * fs_15_4_2 * r_4 * h2_p1 - e_3 * fs_42_143_6 * h8_p1 + e_3 * fs_7_286_4290 * h8_p7 + e_3 * fs_7_11_14 * r_2 * h6_p1 - e_3 * fs_90_143_15 * r_4 * h4_p1 + e_3 * fs_5_11_2 * r_6 * h2_p1 - e_4 * fs_7_8398_330 * h10_p1 + e_4 * fs_7_4199_12155 * h10_p7 + e_4 * fs_84_2717_6 * r_2 * h8_p1 - e_4 * fs_7_2717_4290 * r_2 * h8_p7 - e_4 * fs_7_187_14 * r_4 * h6_p1 + e_4 * fs_4_143_15 * r_6 * h4_p1 - e_4 * fs_5_286_2 * r_8 * h2_p1;

        pc_90[k] = e_3 * fs_14_143_143 * h8_m8 + e_4 * fs_21_4199_2431 * h10_m8 - e_4 * fs_28_2717_143 * r_2 * h8_m8;
    }

    // NOTE: the angular components are formed in 87 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_m1, ph4_m1, ph6_m1, ph8_m7, ph8_m1, ph10_m7, ph10_m1, ab_2, pc_91 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

        const auto h2_m1 = ph2_m1[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h8_m7 = ph8_m7[k];
        const auto h8_m1 = ph8_m1[k];
        const auto h10_m7 = ph10_m7[k];
        const auto h10_m1 = ph10_m1[k];

        pc_91[k] = - e_0 * fs_105_8_30 * h2_m1 - e_1 * f_45_8 * h4_m1 + e_1 * fs_15_1_30 * r_2 * h2_m1 + e_2 * fs_25_44_210 * h6_m1 + e_2 * f_135_44 * r_2 * h4_m1 - e_2 * fs_5_1_30 * r_4 * h2_m1 - e_3 * fs_49_572_286 * h8_m7 + e_3 * fs_175_572_10 * h8_m1 - e_3 * fs_5_33_210 * r_2 * h6_m1 - e_3 * f_135_286 * r_4 * h4_m1 + e_3 * fs_20_33_30 * r_6 * h2_m1 + e_4 * fs_14_4199_7293 * h10_m7 + e_4 * fs_21_4199_22 * h10_m1 + e_4 * fs_49_5434_286 * r_2 * h8_m7 - e_4 * fs_175_5434_10 * r_2 * h8_m1 + e_4 * fs_5_561_210 * r_4 * h6_m1 + e_4 * f_3_143 * r_6 * h4_m1 - e_4 * fs_10_429_30 * r_8 * h2_m1;
    }

    // NOTE: the angular components are formed in 87 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_m2, ph4_m2, ph6_m6, ph6_m2, ph8_m6, ph8_m2, ph10_m6, ph10_m2, ab_2, pc_92 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

        const auto h2_m2 = ph2_m2[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h6_m6 = ph6_m6[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h8_m6 = ph8_m6[k];
        const auto h8_m2 = ph8_m2[k];
        const auto h10_m6 = ph10_m6[k];
        const auto h10_m2 = ph10_m2[k];

        pc_92[k] = - e_0 * fs_105_16_105 * h2_m2 + e_1 * fs_45_4_7 * h4_m2 + e_1 * fs_15_2_105 * r_2 * h2_m2 + e_2 * fs_45_88_330 * h6_m6 - e_2 * fs_115_88_6 * h6_m2 - e_2 * fs_135_22_7 * r_2 * h4_m2 - e_2 * fs_5_2_105 * r_4 * h2_m2 - e_3 * fs_1_143_30030 * h8_m6 - e_3 * fs_139_143_2 * h8_m2 - e_3 * fs_3_22_330 * r_2 * h6_m6 + e_3 * fs_23_66_6 * r_2 * h6_m2 + e_3 * fs_135_143_7 * r_4 * h4_m2 + e_3 * fs_10_33_105 * r_6 * h2_m2 + e_4 * fs_14_4199_6006 * h10_m6 - e_4 * fs_14_4199_231 * h10_m2 + e_4 * fs_2_2717_30030 * r_2 * h8_m6 + e_4 * fs_278_2717_2 * r_2 * h8_m2 + e_4 * fs_3_374_330 * r_4 * h6_m6 - e_4 * fs_23_1122_6 * r_4 * h6_m2 - e_4 * fs_6_143_7 * r_6 * h4_m2 - e_4 * fs_5_429_105 * r_8 * h2_m2;
    }

    // NOTE: the angular components are formed in 87 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, ph4_m3, ph4_p4, ph6_m5, ph6_m3, ph6_p4, ph8_m5, ph8_m3, ph8_p4, ph10_m5, ph10_m3, ph10_p4, ab_2, pc_93, pc_94 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

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

        pc_93[k] = - e_1 * f_45_8 * h4_m3 - e_2 * fs_15_44_55 * h6_m5 - e_2 * fs_105_44_3 * h6_m3 + e_2 * f_135_44 * r_2 * h4_m3 - e_3 * fs_1_572_1430 * h8_m5 + e_3 * fs_101_572_66 * h8_m3 + e_3 * fs_1_11_55 * r_2 * h6_m5 + e_3 * fs_7_11_3 * r_2 * h6_m3 - e_3 * f_135_286 * r_4 * h4_m3 + e_4 * fs_7_4199_15015 * h10_m5 + e_4 * fs_7_4199_3003 * h10_m3 + e_4 * fs_1_5434_1430 * r_2 * h8_m5 - e_4 * fs_101_5434_66 * r_2 * h8_m3 - e_4 * fs_1_187_55 * r_4 * h6_m5 - e_4 * fs_7_187_3 * r_4 * h6_m3 + e_4 * f_3_143 * r_6 * h4_m3;

        pc_94[k] = e_1 * fs_45_2_5 * h4_p4 - e_2 * f_120_11 * h6_p4 - e_2 * fs_135_11_5 * r_2 * h4_p4 + e_3 * fs_58_143_11 * h8_p4 + e_3 * f_32_11 * r_2 * h6_p4 + e_3 * fs_270_143_5 * r_4 * h4_p4 + e_4 * fs_7_4199_15015 * h10_p4 - e_4 * fs_116_2717_11 * r_2 * h8_p4 - e_4 * f_32_187 * r_4 * h6_p4 - e_4 * fs_12_143_5 * r_6 * h4_p4;
    }

    // NOTE: the angular components are formed in 87 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, ph4_p3, ph6_p3, ph6_p5, ph8_p3, ph8_p5, ph10_p3, ph10_p5, ab_2, pc_95 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h4_p3 = ph4_p3[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h6_p5 = ph6_p5[k];
        const auto h8_p3 = ph8_p3[k];
        const auto h8_p5 = ph8_p5[k];
        const auto h10_p3 = ph10_p3[k];
        const auto h10_p5 = ph10_p5[k];

        pc_95[k] = e_1 * f_45_8 * h4_p3 + e_2 * fs_105_44_3 * h6_p3 - e_2 * fs_15_44_55 * h6_p5 - e_2 * f_135_44 * r_2 * h4_p3 - e_3 * fs_101_572_66 * h8_p3 - e_3 * fs_1_572_1430 * h8_p5 - e_3 * fs_7_11_3 * r_2 * h6_p3 + e_3 * fs_1_11_55 * r_2 * h6_p5 + e_3 * f_135_286 * r_4 * h4_p3 - e_4 * fs_7_4199_3003 * h10_p3 + e_4 * fs_7_4199_15015 * h10_p5 + e_4 * fs_101_5434_66 * r_2 * h8_p3 + e_4 * fs_1_5434_1430 * r_2 * h8_p5 + e_4 * fs_7_187_3 * r_4 * h6_p3 - e_4 * fs_1_187_55 * r_4 * h6_p5 - e_4 * f_3_143 * r_6 * h4_p3;
    }

    // NOTE: the angular components are formed in 87 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_p2, ph4_p2, ph6_p2, ph6_p6, ph8_p2, ph8_p6, ph10_p2, ph10_p6, ab_2, pc_96 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

        const auto h2_p2 = ph2_p2[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h6_p6 = ph6_p6[k];
        const auto h8_p2 = ph8_p2[k];
        const auto h8_p6 = ph8_p6[k];
        const auto h10_p2 = ph10_p2[k];
        const auto h10_p6 = ph10_p6[k];

        pc_96[k] = e_0 * fs_105_16_105 * h2_p2 - e_1 * fs_45_4_7 * h4_p2 - e_1 * fs_15_2_105 * r_2 * h2_p2 + e_2 * fs_115_88_6 * h6_p2 + e_2 * fs_45_88_330 * h6_p6 + e_2 * fs_135_22_7 * r_2 * h4_p2 + e_2 * fs_5_2_105 * r_4 * h2_p2 + e_3 * fs_139_143_2 * h8_p2 - e_3 * fs_1_143_30030 * h8_p6 - e_3 * fs_23_66_6 * r_2 * h6_p2 - e_3 * fs_3_22_330 * r_2 * h6_p6 - e_3 * fs_135_143_7 * r_4 * h4_p2 - e_3 * fs_10_33_105 * r_6 * h2_p2 + e_4 * fs_14_4199_231 * h10_p2 + e_4 * fs_14_4199_6006 * h10_p6 - e_4 * fs_278_2717_2 * r_2 * h8_p2 + e_4 * fs_2_2717_30030 * r_2 * h8_p6 + e_4 * fs_23_1122_6 * r_4 * h6_p2 + e_4 * fs_3_374_330 * r_4 * h6_p6 + e_4 * fs_6_143_7 * r_6 * h4_p2 + e_4 * fs_5_429_105 * r_8 * h2_p2;
    }

    // NOTE: the angular components are formed in 87 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_p1, ph4_p1, ph6_p1, ph8_p1, ph8_p7, ph10_p1, ph10_p7, ab_2, pc_97 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

        const auto h2_p1 = ph2_p1[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h8_p7 = ph8_p7[k];
        const auto h10_p1 = ph10_p1[k];
        const auto h10_p7 = ph10_p7[k];

        pc_97[k] = e_0 * fs_105_8_30 * h2_p1 + e_1 * f_45_8 * h4_p1 - e_1 * fs_15_1_30 * r_2 * h2_p1 - e_2 * fs_25_44_210 * h6_p1 - e_2 * f_135_44 * r_2 * h4_p1 + e_2 * fs_5_1_30 * r_4 * h2_p1 - e_3 * fs_175_572_10 * h8_p1 - e_3 * fs_49_572_286 * h8_p7 + e_3 * fs_5_33_210 * r_2 * h6_p1 + e_3 * f_135_286 * r_4 * h4_p1 - e_3 * fs_20_33_30 * r_6 * h2_p1 - e_4 * fs_21_4199_22 * h10_p1 + e_4 * fs_14_4199_7293 * h10_p7 + e_4 * fs_175_5434_10 * r_2 * h8_p1 + e_4 * fs_49_5434_286 * r_2 * h8_p7 - e_4 * fs_5_561_210 * r_4 * h6_p1 - e_4 * f_3_143 * r_6 * h4_p1 + e_4 * fs_10_429_30 * r_8 * h2_p1;
    }

    // NOTE: the angular components are formed in 87 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_0, ph4_0, ph6_0, ph8_0, ph8_p8, ph10_0, ph10_p8, ab_2, pc_98 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

        const auto h2_0 = ph2_0[k];
        const auto h4_0 = ph4_0[k];
        const auto h6_0 = ph6_0[k];
        const auto h8_0 = ph8_0[k];
        const auto h8_p8 = ph8_p8[k];
        const auto h10_0 = ph10_0[k];
        const auto h10_p8 = ph10_p8[k];

        pc_98[k] = e_0 * fs_315_16_5 * h2_0 + e_1 * fs_45_2_5 * h4_0 - e_1 * fs_45_2_5 * r_2 * h2_0 + e_2 * fs_105_22_5 * h6_0 - e_2 * fs_135_11_5 * r_2 * h4_0 + e_2 * fs_15_2_5 * r_4 * h2_0 + e_3 * fs_42_143_5 * h8_0 + e_3 * fs_14_143_143 * h8_p8 - e_3 * fs_14_11_5 * r_2 * h6_0 + e_3 * fs_270_143_5 * r_4 * h4_0 - e_3 * fs_10_11_5 * r_6 * h2_0 + e_4 * fs_21_4199_5 * h10_0 + e_4 * fs_21_4199_2431 * h10_p8 - e_4 * fs_84_2717_5 * r_2 * h8_0 - e_4 * fs_28_2717_143 * r_2 * h8_p8 + e_4 * fs_14_187_5 * r_4 * h6_0 - e_4 * fs_12_143_5 * r_6 * h4_0 + e_4 * fs_5_143_5 * r_8 * h2_0;
    }

    // NOTE: the angular components are formed in 87 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_m1, ph4_m1, ph6_m1, ph8_m1, ph10_m9, ph10_m1, ab_2, pc_99 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

        const auto h2_m1 = ph2_m1[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h8_m1 = ph8_m1[k];
        const auto h10_m9 = ph10_m9[k];
        const auto h10_m1 = ph10_m1[k];

        pc_99[k] = - e_0 * fs_105_32_330 * h2_m1 - e_1 * fs_45_4_11 * h4_m1 + e_1 * fs_15_4_330 * r_2 * h2_m1 - e_2 * fs_5_44_2310 * h6_m1 + e_2 * fs_135_22_11 * r_2 * h4_m1 - e_2 * fs_5_4_330 * r_4 * h2_m1 - e_3 * fs_7_286_110 * h8_m1 + e_3 * fs_1_33_2310 * r_2 * h6_m1 - e_3 * fs_135_143_11 * r_4 * h4_m1 + e_3 * fs_5_33_330 * r_6 * h2_m1 + e_4 * fs_21_4199_4199 * h10_m9 - e_4 * fs_21_8398_2 * h10_m1 + e_4 * fs_7_2717_110 * r_2 * h8_m1 - e_4 * fs_1_561_2310 * r_4 * h6_m1 + e_4 * fs_6_143_11 * r_6 * h4_m1 - e_4 * fs_5_858_330 * r_8 * h2_m1;
    }

    // NOTE: the angular components are formed in 87 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_m2, ph4_m2, ph6_m2, ph8_m8, ph8_m2, ph10_m8, ph10_m2, ab_2, pc_100 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

        const auto h2_m2 = ph2_m2[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h8_m8 = ph8_m8[k];
        const auto h8_m2 = ph8_m2[k];
        const auto h10_m8 = ph10_m8[k];
        const auto h10_m2 = ph10_m2[k];

        pc_100[k] = - e_0 * fs_105_16_165 * h2_m2 + e_1 * fs_45_8_11 * h4_m2 + e_1 * fs_15_2_165 * r_2 * h2_m2 + e_2 * fs_35_88_462 * h6_m2 - e_2 * fs_135_44_11 * r_2 * h4_m2 - e_2 * fs_5_2_165 * r_4 * h2_m2 - e_3 * fs_7_13_13 * h8_m8 + e_3 * fs_29_572_154 * h8_m2 - e_3 * fs_7_66_462 * r_2 * h6_m2 + e_3 * fs_135_286_11 * r_4 * h4_m2 + e_3 * fs_10_33_165 * r_6 * h2_m2 + e_4 * fs_84_4199_221 * h10_m8 + e_4 * fs_28_4199_3 * h10_m2 + e_4 * fs_14_247_13 * r_2 * h8_m8 - e_4 * fs_29_5434_154 * r_2 * h8_m2 + e_4 * fs_7_1122_462 * r_4 * h6_m2 - e_4 * fs_3_143_11 * r_6 * h4_m2 - e_4 * fs_5_429_165 * r_8 * h2_m2;
    }

    // NOTE: the angular components are formed in 87 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, ph4_m4, ph4_m3, ph6_m6, ph6_m4, ph6_m3, ph8_m7, ph8_m6, ph8_m4, ph8_m3, ph10_m7, ph10_m6, ph10_m4, ph10_m3, ab_2, pc_101, pc_102 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h4_m4 = ph4_m4[k];
        const auto h4_m3 = ph4_m3[k];
        const auto h6_m6 = ph6_m6[k];
        const auto h6_m4 = ph6_m4[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h8_m7 = ph8_m7[k];
        const auto h8_m6 = ph8_m6[k];
        const auto h8_m4 = ph8_m4[k];
        const auto h8_m3 = ph8_m3[k];
        const auto h10_m7 = ph10_m7[k];
        const auto h10_m6 = ph10_m6[k];
        const auto h10_m4 = ph10_m4[k];
        const auto h10_m3 = ph10_m3[k];

        pc_101[k] = e_1 * fs_45_8_11 * h4_m3 - e_2 * fs_15_11_33 * h6_m3 - e_2 * fs_135_44_11 * r_2 * h4_m3 - e_3 * fs_1_52_182 * h8_m7 - e_3 * fs_23_52_6 * h8_m3 + e_3 * fs_4_11_33 * r_2 * h6_m3 + e_3 * fs_135_286_11 * r_4 * h4_m3 + e_4 * fs_14_4199_4641 * h10_m7 - e_4 * fs_7_4199_273 * h10_m3 + e_4 * fs_1_494_182 * r_2 * h8_m7 + e_4 * fs_23_494_6 * r_2 * h8_m3 - e_4 * fs_4_187_33 * r_4 * h6_m3 - e_4 * fs_3_143_11 * r_6 * h4_m3;

        pc_102[k] = - e_1 * fs_45_4_11 * h4_m4 - e_2 * fs_15_8_30 * h6_m6 + e_2 * fs_15_44_55 * h6_m4 + e_2 * fs_135_22_11 * r_2 * h4_m4 + e_3 * fs_1_52_2730 * h8_m6 + e_3 * fs_17_26_5 * h8_m4 + e_3 * fs_1_2_30 * r_2 * h6_m6 - e_3 * fs_1_11_55 * r_2 * h6_m4 - e_3 * fs_135_143_11 * r_4 * h4_m4 + e_4 * fs_28_4199_546 * h10_m6 + e_4 * fs_14_4199_273 * h10_m4 - e_4 * fs_1_494_2730 * r_2 * h8_m6 - e_4 * fs_17_247_5 * r_2 * h8_m4 - e_4 * fs_1_34_30 * r_4 * h6_m6 + e_4 * fs_1_187_55 * r_4 * h6_m4 + e_4 * fs_6_143_11 * r_6 * h4_m4;
    }

    // NOTE: the angular components are formed in 87 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, ph4_p4, ph6_p4, ph6_p5, ph6_p6, ph8_p4, ph8_p5, ph8_p6, ph10_p4, ph10_p5, ph10_p6, ab_2, pc_103, pc_104 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h4_p4 = ph4_p4[k];
        const auto h6_p4 = ph6_p4[k];
        const auto h6_p5 = ph6_p5[k];
        const auto h6_p6 = ph6_p6[k];
        const auto h8_p4 = ph8_p4[k];
        const auto h8_p5 = ph8_p5[k];
        const auto h8_p6 = ph8_p6[k];
        const auto h10_p4 = ph10_p4[k];
        const auto h10_p5 = ph10_p5[k];
        const auto h10_p6 = ph10_p6[k];

        pc_103[k] = - e_2 * f_15_2 * h6_p5 + e_3 * fs_11_26_26 * h8_p5 + e_3 * f_2_1 * r_2 * h6_p5 + e_4 * fs_35_4199_273 * h10_p5 - e_4 * fs_11_247_26 * r_2 * h8_p5 - e_4 * f_2_17 * r_4 * h6_p5;

        pc_104[k] = e_1 * fs_45_4_11 * h4_p4 - e_2 * fs_15_44_55 * h6_p4 - e_2 * fs_15_8_30 * h6_p6 - e_2 * fs_135_22_11 * r_2 * h4_p4 - e_3 * fs_17_26_5 * h8_p4 + e_3 * fs_1_52_2730 * h8_p6 + e_3 * fs_1_11_55 * r_2 * h6_p4 + e_3 * fs_1_2_30 * r_2 * h6_p6 + e_3 * fs_135_143_11 * r_4 * h4_p4 - e_4 * fs_14_4199_273 * h10_p4 + e_4 * fs_28_4199_546 * h10_p6 + e_4 * fs_17_247_5 * r_2 * h8_p4 - e_4 * fs_1_494_2730 * r_2 * h8_p6 - e_4 * fs_1_187_55 * r_4 * h6_p4 - e_4 * fs_1_34_30 * r_4 * h6_p6 - e_4 * fs_6_143_11 * r_6 * h4_p4;
    }

    // NOTE: the angular components are formed in 87 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, ph4_p3, ph6_p3, ph8_p3, ph8_p7, ph10_p3, ph10_p7, ab_2, pc_105 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h4_p3 = ph4_p3[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h8_p3 = ph8_p3[k];
        const auto h8_p7 = ph8_p7[k];
        const auto h10_p3 = ph10_p3[k];
        const auto h10_p7 = ph10_p7[k];

        pc_105[k] = - e_1 * fs_45_8_11 * h4_p3 + e_2 * fs_15_11_33 * h6_p3 + e_2 * fs_135_44_11 * r_2 * h4_p3 + e_3 * fs_23_52_6 * h8_p3 - e_3 * fs_1_52_182 * h8_p7 - e_3 * fs_4_11_33 * r_2 * h6_p3 - e_3 * fs_135_286_11 * r_4 * h4_p3 + e_4 * fs_7_4199_273 * h10_p3 + e_4 * fs_14_4199_4641 * h10_p7 - e_4 * fs_23_494_6 * r_2 * h8_p3 + e_4 * fs_1_494_182 * r_2 * h8_p7 + e_4 * fs_4_187_33 * r_4 * h6_p3 + e_4 * fs_3_143_11 * r_6 * h4_p3;
    }

    // NOTE: the angular components are formed in 87 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_p2, ph4_p2, ph6_p2, ph8_p2, ph8_p8, ph10_p2, ph10_p8, ab_2, pc_106 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

        const auto h2_p2 = ph2_p2[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h8_p2 = ph8_p2[k];
        const auto h8_p8 = ph8_p8[k];
        const auto h10_p2 = ph10_p2[k];
        const auto h10_p8 = ph10_p8[k];

        pc_106[k] = e_0 * fs_105_16_165 * h2_p2 - e_1 * fs_45_8_11 * h4_p2 - e_1 * fs_15_2_165 * r_2 * h2_p2 - e_2 * fs_35_88_462 * h6_p2 + e_2 * fs_135_44_11 * r_2 * h4_p2 + e_2 * fs_5_2_165 * r_4 * h2_p2 - e_3 * fs_29_572_154 * h8_p2 - e_3 * fs_7_13_13 * h8_p8 + e_3 * fs_7_66_462 * r_2 * h6_p2 - e_3 * fs_135_286_11 * r_4 * h4_p2 - e_3 * fs_10_33_165 * r_6 * h2_p2 - e_4 * fs_28_4199_3 * h10_p2 + e_4 * fs_84_4199_221 * h10_p8 + e_4 * fs_29_5434_154 * r_2 * h8_p2 + e_4 * fs_14_247_13 * r_2 * h8_p8 - e_4 * fs_7_1122_462 * r_4 * h6_p2 + e_4 * fs_3_143_11 * r_6 * h4_p2 + e_4 * fs_5_429_165 * r_8 * h2_p2;
    }

    // NOTE: the angular components are formed in 87 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_m2, ph2_p1, ph4_m2, ph4_p1, ph6_m2, ph6_p1, ph8_m2, ph8_p1, ph10_m10, ph10_m2, ph10_p1, ph10_p9, ab_2, pc_107, pc_108 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

        const auto h2_m2 = ph2_m2[k];
        const auto h2_p1 = ph2_p1[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h8_m2 = ph8_m2[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h10_m10 = ph10_m10[k];
        const auto h10_m2 = ph10_m2[k];
        const auto h10_p1 = ph10_p1[k];
        const auto h10_p9 = ph10_p9[k];

        pc_107[k] = e_0 * fs_105_32_330 * h2_p1 + e_1 * fs_45_4_11 * h4_p1 - e_1 * fs_15_4_330 * r_2 * h2_p1 + e_2 * fs_5_44_2310 * h6_p1 - e_2 * fs_135_22_11 * r_2 * h4_p1 + e_2 * fs_5_4_330 * r_4 * h2_p1 + e_3 * fs_7_286_110 * h8_p1 - e_3 * fs_1_33_2310 * r_2 * h6_p1 + e_3 * fs_135_143_11 * r_4 * h4_p1 - e_3 * fs_5_33_330 * r_6 * h2_p1 + e_4 * fs_21_8398_2 * h10_p1 + e_4 * fs_21_4199_4199 * h10_p9 - e_4 * fs_7_2717_110 * r_2 * h8_p1 + e_4 * fs_1_561_2310 * r_4 * h6_p1 - e_4 * fs_6_143_11 * r_6 * h4_p1 + e_4 * fs_5_858_330 * r_8 * h2_p1;

        pc_108[k] = - e_0 * fs_315_32_110 * h2_m2 - e_1 * fs_15_4_66 * h4_m2 + e_1 * fs_45_4_110 * r_2 * h2_m2 - e_2 * fs_15_44_77 * h6_m2 + e_2 * fs_45_22_66 * r_2 * h4_m2 - e_2 * fs_15_4_110 * r_4 * h2_m2 - e_3 * fs_1_143_231 * h8_m2 + e_3 * fs_1_11_77 * r_2 * h6_m2 - e_3 * fs_45_143_66 * r_4 * h4_m2 + e_3 * fs_5_11_110 * r_6 * h2_m2 + e_4 * fs_7_4199_62985 * h10_m10 - e_4 * fs_7_8398_2 * h10_m2 + e_4 * fs_2_2717_231 * r_2 * h8_m2 - e_4 * fs_1_187_77 * r_4 * h6_m2 + e_4 * fs_2_143_66 * r_6 * h4_m2 - e_4 * fs_5_286_110 * r_8 * h2_m2;
    }

    // NOTE: the angular components are formed in 87 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, ph4_m4, ph4_m3, ph6_m4, ph6_m3, ph8_m8, ph8_m4, ph8_m3, ph10_m9, ph10_m8, ph10_m4, ph10_m3, ab_2, pc_109, pc_110 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

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

        pc_109[k] = e_1 * fs_15_8_462 * h4_m3 + e_2 * fs_45_88_154 * h6_m3 - e_2 * fs_45_44_462 * r_2 * h4_m3 + e_3 * fs_3_26_7 * h8_m3 - e_3 * fs_3_22_154 * r_2 * h6_m3 + e_3 * fs_45_286_462 * r_4 * h4_m3 + e_4 * fs_7_4199_25194 * h10_m9 + e_4 * fs_7_8398_26 * h10_m3 - e_4 * fs_3_247_7 * r_2 * h8_m3 + e_4 * fs_3_374_154 * r_4 * h6_m3 - e_4 * fs_1_143_462 * r_6 * h4_m3;

        pc_110[k] = - e_1 * fs_15_4_66 * h4_m4 - e_2 * fs_45_88_330 * h6_m4 + e_2 * fs_45_22_66 * r_2 * h4_m4 + e_3 * fs_1_13_546 * h8_m8 - e_3 * fs_3_26_30 * h8_m4 + e_3 * fs_3_22_330 * r_2 * h6_m4 - e_3 * fs_45_143_66 * r_4 * h4_m4 + e_4 * fs_7_4199_9282 * h10_m8 - e_4 * fs_7_8398_182 * h10_m4 - e_4 * fs_2_247_546 * r_2 * h8_m8 + e_4 * fs_3_247_30 * r_2 * h8_m4 - e_4 * fs_3_374_330 * r_4 * h6_m4 + e_4 * fs_2_143_66 * r_6 * h4_m4;
    }

    // NOTE: the angular components are formed in 87 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_2, pe_3, pe_4, ph6_m5, ph6_p5, ph6_p6, ph8_m7, ph8_m5, ph8_p5, ph8_p6, ph8_p7, ph10_m7, ph10_m5, ph10_p5, ph10_p6, ph10_p7, ab_2, pc_111, pc_112, pc_113 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

        const auto h6_m5 = ph6_m5[k];
        const auto h6_p5 = ph6_p5[k];
        const auto h6_p6 = ph6_p6[k];
        const auto h8_m7 = ph8_m7[k];
        const auto h8_m5 = ph8_m5[k];
        const auto h8_p5 = ph8_p5[k];
        const auto h8_p6 = ph8_p6[k];
        const auto h8_p7 = ph8_p7[k];
        const auto h10_m7 = ph10_m7[k];
        const auto h10_m5 = ph10_m5[k];
        const auto h10_p5 = ph10_p5[k];
        const auto h10_p6 = ph10_p6[k];
        const auto h10_p7 = ph10_p7[k];

        pc_111[k] = e_2 * fs_15_8_30 * h6_m5 + e_3 * fs_3_26_273 * h8_m7 + e_3 * fs_1_13_195 * h8_m5 - e_3 * fs_1_2_30 * r_2 * h6_m5 + e_4 * fs_7_4199_3094 * h10_m7 + e_4 * fs_7_8398_910 * h10_m5 - e_4 * fs_3_247_273 * r_2 * h8_m7 - e_4 * fs_2_247_195 * r_2 * h8_m5 + e_4 * fs_1_34_30 * r_4 * h6_m5;

        pc_112[k] = e_2 * f_45_4 * h6_p6 + e_3 * fs_3_13_91 * h8_p6 - e_3 * f_3_1 * r_2 * h6_p6 + e_4 * fs_14_4199_455 * h10_p6 - e_4 * fs_6_247_91 * r_2 * h8_p6 + e_4 * f_3_17 * r_4 * h6_p6;

        pc_113[k] = - e_2 * fs_15_8_30 * h6_p5 - e_3 * fs_1_13_195 * h8_p5 + e_3 * fs_3_26_273 * h8_p7 + e_3 * fs_1_2_30 * r_2 * h6_p5 - e_4 * fs_7_8398_910 * h10_p5 + e_4 * fs_7_4199_3094 * h10_p7 + e_4 * fs_2_247_195 * r_2 * h8_p5 - e_4 * fs_3_247_273 * r_2 * h8_p7 - e_4 * fs_1_34_30 * r_4 * h6_p5;
    }

    // NOTE: the angular components are formed in 87 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, ph4_p3, ph4_p4, ph6_p3, ph6_p4, ph8_p3, ph8_p4, ph8_p8, ph10_p3, ph10_p4, ph10_p8, ph10_p9, ab_2, pc_114, pc_115 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

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

        pc_114[k] = e_1 * fs_15_4_66 * h4_p4 + e_2 * fs_45_88_330 * h6_p4 - e_2 * fs_45_22_66 * r_2 * h4_p4 + e_3 * fs_3_26_30 * h8_p4 + e_3 * fs_1_13_546 * h8_p8 - e_3 * fs_3_22_330 * r_2 * h6_p4 + e_3 * fs_45_143_66 * r_4 * h4_p4 + e_4 * fs_7_8398_182 * h10_p4 + e_4 * fs_7_4199_9282 * h10_p8 - e_4 * fs_3_247_30 * r_2 * h8_p4 - e_4 * fs_2_247_546 * r_2 * h8_p8 + e_4 * fs_3_374_330 * r_4 * h6_p4 - e_4 * fs_2_143_66 * r_6 * h4_p4;

        pc_115[k] = - e_1 * fs_15_8_462 * h4_p3 - e_2 * fs_45_88_154 * h6_p3 + e_2 * fs_45_44_462 * r_2 * h4_p3 - e_3 * fs_3_26_7 * h8_p3 + e_3 * fs_3_22_154 * r_2 * h6_p3 - e_3 * fs_45_286_462 * r_4 * h4_p3 - e_4 * fs_7_8398_26 * h10_p3 + e_4 * fs_7_4199_25194 * h10_p9 + e_4 * fs_3_247_7 * r_2 * h8_p3 - e_4 * fs_3_374_154 * r_4 * h6_p3 + e_4 * fs_1_143_462 * r_6 * h4_p3;
    }

    // NOTE: the angular components are formed in 87 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_p2, ph4_p2, ph6_p2, ph8_p2, ph10_p2, ph10_p10, ab_2, pc_116 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

        const auto h2_p2 = ph2_p2[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h8_p2 = ph8_p2[k];
        const auto h10_p2 = ph10_p2[k];
        const auto h10_p10 = ph10_p10[k];

        pc_116[k] = e_0 * fs_315_32_110 * h2_p2 + e_1 * fs_15_4_66 * h4_p2 - e_1 * fs_45_4_110 * r_2 * h2_p2 + e_2 * fs_15_44_77 * h6_p2 - e_2 * fs_45_22_66 * r_2 * h4_p2 + e_2 * fs_15_4_110 * r_4 * h2_p2 + e_3 * fs_1_143_231 * h8_p2 - e_3 * fs_1_11_77 * r_2 * h6_p2 + e_3 * fs_45_143_66 * r_4 * h4_p2 - e_3 * fs_5_11_110 * r_6 * h2_p2 + e_4 * fs_7_8398_2 * h10_p2 + e_4 * fs_7_4199_62985 * h10_p10 - e_4 * fs_2_2717_231 * r_2 * h8_p2 + e_4 * fs_1_187_77 * r_4 * h6_p2 - e_4 * fs_2_143_66 * r_6 * h4_p2 + e_4 * fs_5_286_110 * r_8 * h2_p2;
    }

    // NOTE: the values of a combination of angular components are stored as one
    // row of nvalues columns, with the component on bra side running slowest. The
    // rows which the symmetry relates to an already formed one are copied from it.

    const size_t sources[117] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116};

    for (size_t n = 0; n < 117; n++)
    {
        if (sources[n] != n) std::copy(values + sources[n] * nvalues, values + (sources[n] + 1) * nvalues, values + n * nvalues);
    }
}

}  // namespace simdt2ceri
