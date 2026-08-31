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



#include "SimdTwoCenterElectronRepulsionRecHG.hpp"

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
compute_hg_electron_repulsion(double                         *values,
                              const size_t                    nvalues,
                              const CBasisFunction           &bra,
                              const CBasisFunction           &ket,
                              const std::vector<CSimdMatrix> &harmonics,
                              const CSimdMatrix              &coordinates) -> void
{
    if ((bra.get_angular_momentum() != 5) || (ket.get_angular_momentum() != 4))
    {
        errors::assertMsgCritical(
            false, std::string("SimdTwoCenterElectronRepulsionFunc.compute_hg_electron_repulsion: Basis functions must be of angular momenta five and four"));
    }

    if (harmonics.size() < 9)
    {
        errors::assertMsgCritical(
            false, std::string("SimdTwoCenterElectronRepulsionFunc.compute_hg_electron_repulsion: Harmonics must reach angular momentum 9"));
    }

    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("SimdTwoCenterElectronRepulsionFunc.compute_hg_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    auto boys = CSimdVariableMatrix(std::vector<size_t>(nprims, nvalues), 11);

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
    // call, which fills the orders 0 to 9 of every row. The terms read the
    // orders 5 to 9 alone, and the orders below them are formed on the
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

            const auto ff_0 = fbase * bexp / (fexp * fexp * fexp * fexp * fexp);

            const auto ff_1 = fbase * bexp * fmu / (fexp * fexp * fexp * fexp * fexp);

            const auto ff_2 = fbase * bexp * fmu * fmu / (fexp * fexp * fexp * fexp * fexp);

            const auto ff_3 = fbase * bexp * fmu * fmu * fmu / (fexp * fexp * fexp * fexp * fexp);

            const auto ff_4 = fbase * bexp * fmu * fmu * fmu * fmu / (fexp * fexp * fexp * fexp * fexp);

            const auto *bv_0 = boys.data(6, i * nprim_b + j);

            const auto *bv_1 = boys.data(7, i * nprim_b + j);

            const auto *bv_2 = boys.data(8, i * nprim_b + j);

            const auto *bv_3 = boys.data(9, i * nprim_b + j);

            const auto *bv_4 = boys.data(10, i * nprim_b + j);

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

    const auto f_105_2 = 105.0 / 2.0;
    const auto f_105_4 = 105.0 / 4.0;
    const auto f_10_3 = 10.0 / 3.0;
    const auto f_126_2431 = 126.0 / 2431.0;
    const auto f_1330_7293 = 1330.0 / 7293.0;
    const auto f_1400_7293 = 1400.0 / 7293.0;
    const auto f_140_2431 = 140.0 / 2431.0;
    const auto f_147_8 = 147.0 / 8.0;
    const auto f_14_143 = 14.0 / 143.0;
    const auto f_15_1 = 15.0;
    const auto f_15_2 = 15.0 / 2.0;
    const auto f_175_143 = 175.0 / 143.0;
    const auto f_18_1 = 18.0;
    const auto f_1_11 = 1.0 / 11.0;
    const auto f_1_39 = 1.0 / 39.0;
    const auto f_20_13 = 20.0 / 13.0;
    const auto f_20_143 = 20.0 / 143.0;
    const auto f_21_1 = 21.0;
    const auto f_21_11 = 21.0 / 11.0;
    const auto f_21_2 = 21.0 / 2.0;
    const auto f_21_2431 = 21.0 / 2431.0;
    const auto f_27_2 = 27.0 / 2.0;
    const auto f_28_143 = 28.0 / 143.0;
    const auto f_2_1 = 2.0;
    const auto f_2_13 = 2.0 / 13.0;
    const auto f_30_11 = 30.0 / 11.0;
    const auto f_30_13 = 30.0 / 13.0;
    const auto f_315_16 = 315.0 / 16.0;
    const auto f_350_2431 = 350.0 / 2431.0;
    const auto f_35_221 = 35.0 / 221.0;
    const auto f_35_26 = 35.0 / 26.0;
    const auto f_42_1 = 42.0;
    const auto f_42_11 = 42.0 / 11.0;
    const auto f_45_2 = 45.0 / 2.0;
    const auto f_49_22 = 49.0 / 22.0;
    const auto f_49_4 = 49.0 / 4.0;
    const auto f_49_429 = 49.0 / 429.0;
    const auto f_4_33 = 4.0 / 33.0;
    const auto f_4_39 = 4.0 / 39.0;
    const auto f_525_16 = 525.0 / 16.0;
    const auto f_5_1 = 5.0;
    const auto f_5_13 = 5.0 / 13.0;
    const auto f_5_33 = 5.0 / 33.0;
    const auto f_5_4 = 5.0 / 4.0;
    const auto f_63_2 = 63.0 / 2.0;
    const auto f_63_4 = 63.0 / 4.0;
    const auto f_665_429 = 665.0 / 429.0;
    const auto f_700_429 = 700.0 / 429.0;
    const auto f_70_143 = 70.0 / 143.0;
    const auto f_882_2431 = 882.0 / 2431.0;
    const auto f_8_3 = 8.0 / 3.0;
    const auto fs_100_2431_11 = std::sqrt(10000.0 / 537251.0);
    const auto fs_105_16_10 = std::sqrt(55125.0 / 128.0);
    const auto fs_105_16_14 = std::sqrt(77175.0 / 128.0);
    const auto fs_105_16_15 = std::sqrt(165375.0 / 256.0);
    const auto fs_105_16_21 = std::sqrt(231525.0 / 256.0);
    const auto fs_105_16_3 = std::sqrt(33075.0 / 256.0);
    const auto fs_105_2431_33 = std::sqrt(33075.0 / 537251.0);
    const auto fs_105_32_2 = std::sqrt(11025.0 / 512.0);
    const auto fs_105_32_42 = std::sqrt(231525.0 / 512.0);
    const auto fs_105_32_6 = std::sqrt(33075.0 / 512.0);
    const auto fs_105_8_3 = std::sqrt(33075.0 / 64.0);
    const auto fs_105_8_6 = std::sqrt(33075.0 / 32.0);
    const auto fs_10_13_5 = std::sqrt(500.0 / 169.0);
    const auto fs_10_143_7 = std::sqrt(700.0 / 20449.0);
    const auto fs_10_1_3 = std::sqrt(300.0);
    const auto fs_10_2431_2002 = std::sqrt(1400.0 / 41327.0);
    const auto fs_10_429_7 = std::sqrt(700.0 / 184041.0);
    const auto fs_10_663_105 = std::sqrt(3500.0 / 146523.0);
    const auto fs_10_7293_210 = std::sqrt(7000.0 / 17729283.0);
    const auto fs_115_4862_14 = std::sqrt(92575.0 / 11819522.0);
    const auto fs_115_572_14 = std::sqrt(92575.0 / 163592.0);
    const auto fs_115_7293_30 = std::sqrt(132250.0 / 17729283.0);
    const auto fs_115_7293_70 = std::sqrt(925750.0 / 53187849.0);
    const auto fs_115_858_30 = std::sqrt(66125.0 / 122694.0);
    const auto fs_115_858_70 = std::sqrt(462875.0 / 368082.0);
    const auto fs_11_4_10 = std::sqrt(605.0 / 8.0);
    const auto fs_11_8_70 = std::sqrt(4235.0 / 32.0);
    const auto fs_125_1122_2 = std::sqrt(15625.0 / 629442.0);
    const auto fs_125_132_2 = std::sqrt(15625.0 / 8712.0);
    const auto fs_125_2431_7 = std::sqrt(109375.0 / 5909761.0);
    const auto fs_125_286_7 = std::sqrt(109375.0 / 81796.0);
    const auto fs_129_16_2 = std::sqrt(16641.0 / 128.0);
    const auto fs_135_2431_6 = std::sqrt(109350.0 / 5909761.0);
    const auto fs_135_286_6 = std::sqrt(54675.0 / 40898.0);
    const auto fs_140_2431_11 = std::sqrt(19600.0 / 537251.0);
    const auto fs_140_7293_14 = std::sqrt(274400.0 / 53187849.0);
    const auto fs_140_7293_165 = std::sqrt(98000.0 / 1611753.0);
    const auto fs_145_7293_21 = std::sqrt(147175.0 / 17729283.0);
    const auto fs_145_858_21 = std::sqrt(147175.0 / 245388.0);
    const auto fs_147_2431_22 = std::sqrt(43218.0 / 537251.0);
    const auto fs_14_11_3 = std::sqrt(588.0 / 121.0);
    const auto fs_14_2431_10 = std::sqrt(1960.0 / 5909761.0);
    const auto fs_14_2431_1001 = std::sqrt(1372.0 / 41327.0);
    const auto fs_14_2431_165 = std::sqrt(2940.0 / 537251.0);
    const auto fs_14_2431_2145 = std::sqrt(2940.0 / 41327.0);
    const auto fs_14_2431_2431 = std::sqrt(196.0 / 2431.0);
    const auto fs_14_429_15 = std::sqrt(980.0 / 61347.0);
    const auto fs_14_429_5 = std::sqrt(980.0 / 184041.0);
    const auto fs_15_11_7 = std::sqrt(1575.0 / 121.0);
    const auto fs_15_1_3 = std::sqrt(675.0);
    const auto fs_15_22_7 = std::sqrt(1575.0 / 484.0);
    const auto fs_15_2_7 = std::sqrt(1575.0 / 4.0);
    const auto fs_15_44_30 = std::sqrt(3375.0 / 968.0);
    const auto fs_15_4862_210 = std::sqrt(23625.0 / 11819522.0);
    const auto fs_15_4862_4290 = std::sqrt(3375.0 / 82654.0);
    const auto fs_15_4_7 = std::sqrt(1575.0 / 16.0);
    const auto fs_15_572_210 = std::sqrt(23625.0 / 163592.0);
    const auto fs_15_572_4290 = std::sqrt(3375.0 / 1144.0);
    const auto fs_15_8_105 = std::sqrt(23625.0 / 64.0);
    const auto fs_15_8_30 = std::sqrt(3375.0 / 32.0);
    const auto fs_160_429_21 = std::sqrt(179200.0 / 61347.0);
    const auto fs_16_429_7 = std::sqrt(1792.0 / 184041.0);
    const auto fs_175_7293_42 = std::sqrt(428750.0 / 17729283.0);
    const auto fs_175_858_42 = std::sqrt(214375.0 / 122694.0);
    const auto fs_19_22_6 = std::sqrt(1083.0 / 242.0);
    const auto fs_19_429_6 = std::sqrt(722.0 / 61347.0);
    const auto fs_19_4_6 = std::sqrt(1083.0 / 8.0);
    const auto fs_1_11_2 = std::sqrt(2.0 / 121.0);
    const auto fs_1_11_21 = std::sqrt(21.0 / 121.0);
    const auto fs_1_1_10 = std::sqrt(10.0);
    const auto fs_1_1_105 = std::sqrt(105.0);
    const auto fs_1_22_10 = std::sqrt(5.0 / 242.0);
    const auto fs_1_22_105 = std::sqrt(105.0 / 484.0);
    const auto fs_1_22_15 = std::sqrt(15.0 / 484.0);
    const auto fs_1_22_210 = std::sqrt(105.0 / 242.0);
    const auto fs_1_2_10 = std::sqrt(5.0 / 2.0);
    const auto fs_1_2_21 = std::sqrt(21.0 / 4.0);
    const auto fs_1_33_10 = std::sqrt(10.0 / 1089.0);
    const auto fs_1_33_14 = std::sqrt(14.0 / 1089.0);
    const auto fs_1_33_15 = std::sqrt(5.0 / 363.0);
    const auto fs_1_33_21 = std::sqrt(7.0 / 363.0);
    const auto fs_1_33_3 = std::sqrt(1.0 / 363.0);
    const auto fs_1_39_10 = std::sqrt(10.0 / 1521.0);
    const auto fs_1_39_15 = std::sqrt(5.0 / 507.0);
    const auto fs_1_39_21 = std::sqrt(7.0 / 507.0);
    const auto fs_1_39_30 = std::sqrt(10.0 / 507.0);
    const auto fs_1_39_35 = std::sqrt(35.0 / 1521.0);
    const auto fs_1_39_6 = std::sqrt(2.0 / 507.0);
    const auto fs_1_3_2 = std::sqrt(2.0 / 9.0);
    const auto fs_1_3_42 = std::sqrt(14.0 / 3.0);
    const auto fs_1_3_6 = std::sqrt(2.0 / 3.0);
    const auto fs_1_429_105 = std::sqrt(35.0 / 61347.0);
    const auto fs_1_429_15 = std::sqrt(5.0 / 61347.0);
    const auto fs_1_429_210 = std::sqrt(70.0 / 61347.0);
    const auto fs_1_4_105 = std::sqrt(105.0 / 16.0);
    const auto fs_1_4_15 = std::sqrt(15.0 / 16.0);
    const auto fs_1_4_210 = std::sqrt(105.0 / 8.0);
    const auto fs_1_4_70 = std::sqrt(35.0 / 8.0);
    const auto fs_1_66_2 = std::sqrt(1.0 / 2178.0);
    const auto fs_1_66_42 = std::sqrt(7.0 / 726.0);
    const auto fs_1_66_6 = std::sqrt(1.0 / 726.0);
    const auto fs_1_78_10 = std::sqrt(5.0 / 3042.0);
    const auto fs_1_78_70 = std::sqrt(35.0 / 3042.0);
    const auto fs_205_7293_6 = std::sqrt(84050.0 / 17729283.0);
    const auto fs_205_858_6 = std::sqrt(42025.0 / 122694.0);
    const auto fs_20_11_3 = std::sqrt(1200.0 / 121.0);
    const auto fs_20_143_10 = std::sqrt(4000.0 / 20449.0);
    const auto fs_20_143_165 = std::sqrt(6000.0 / 1859.0);
    const auto fs_215_2431_2 = std::sqrt(92450.0 / 5909761.0);
    const auto fs_215_286_2 = std::sqrt(46225.0 / 40898.0);
    const auto fs_21_1_6 = std::sqrt(2646.0);
    const auto fs_21_2431_1001 = std::sqrt(3087.0 / 41327.0);
    const auto fs_21_2431_143 = std::sqrt(441.0 / 41327.0);
    const auto fs_21_2431_2431 = std::sqrt(441.0 / 2431.0);
    const auto fs_21_2431_70 = std::sqrt(30870.0 / 5909761.0);
    const auto fs_21_2431_715 = std::sqrt(2205.0 / 41327.0);
    const auto fs_21_2431_770 = std::sqrt(30870.0 / 537251.0);
    const auto fs_21_2_10 = std::sqrt(2205.0 / 2.0);
    const auto fs_21_2_14 = std::sqrt(3087.0 / 2.0);
    const auto fs_21_2_15 = std::sqrt(6615.0 / 4.0);
    const auto fs_21_2_21 = std::sqrt(9261.0 / 4.0);
    const auto fs_21_2_3 = std::sqrt(1323.0 / 4.0);
    const auto fs_21_4862_10 = std::sqrt(2205.0 / 11819522.0);
    const auto fs_21_4862_110 = std::sqrt(2205.0 / 1074502.0);
    const auto fs_21_4862_1430 = std::sqrt(2205.0 / 82654.0);
    const auto fs_21_4862_4290 = std::sqrt(6615.0 / 82654.0);
    const auto fs_21_4862_770 = std::sqrt(15435.0 / 1074502.0);
    const auto fs_21_4_15 = std::sqrt(6615.0 / 16.0);
    const auto fs_21_4_2 = std::sqrt(441.0 / 8.0);
    const auto fs_21_4_42 = std::sqrt(9261.0 / 8.0);
    const auto fs_21_4_5 = std::sqrt(2205.0 / 16.0);
    const auto fs_21_4_6 = std::sqrt(1323.0 / 8.0);
    const auto fs_21_8_3 = std::sqrt(1323.0 / 64.0);
    const auto fs_245_7293_30 = std::sqrt(600250.0 / 17729283.0);
    const auto fs_245_858_30 = std::sqrt(300125.0 / 122694.0);
    const auto fs_25_14586_6006 = std::sqrt(4375.0 / 247962.0);
    const auto fs_25_1716_6006 = std::sqrt(4375.0 / 3432.0);
    const auto fs_25_2431_165 = std::sqrt(9375.0 / 537251.0);
    const auto fs_25_286_165 = std::sqrt(9375.0 / 7436.0);
    const auto fs_25_374_6 = std::sqrt(1875.0 / 69938.0);
    const auto fs_25_44_10 = std::sqrt(3125.0 / 968.0);
    const auto fs_25_44_6 = std::sqrt(1875.0 / 968.0);
    const auto fs_25_858_10 = std::sqrt(3125.0 / 368082.0);
    const auto fs_25_8_10 = std::sqrt(3125.0 / 32.0);
    const auto fs_265_14586_66 = std::sqrt(70225.0 / 3223506.0);
    const auto fs_265_1716_66 = std::sqrt(70225.0 / 44616.0);
    const auto fs_27_2_2 = std::sqrt(729.0 / 2.0);
    const auto fs_27_4_10 = std::sqrt(3645.0 / 8.0);
    const auto fs_27_8_21 = std::sqrt(15309.0 / 64.0);
    const auto fs_27_8_35 = std::sqrt(25515.0 / 64.0);
    const auto fs_28_2431_210 = std::sqrt(164640.0 / 5909761.0);
    const auto fs_28_2431_30 = std::sqrt(23520.0 / 5909761.0);
    const auto fs_28_2431_429 = std::sqrt(2352.0 / 41327.0);
    const auto fs_28_2431_715 = std::sqrt(3920.0 / 41327.0);
    const auto fs_28_429_3 = std::sqrt(784.0 / 61347.0);
    const auto fs_294_2431_6 = std::sqrt(518616.0 / 5909761.0);
    const auto fs_2_11_105 = std::sqrt(420.0 / 121.0);
    const auto fs_2_1_2 = std::sqrt(8.0);
    const auto fs_2_33_6 = std::sqrt(8.0 / 363.0);
    const auto fs_2_39_5 = std::sqrt(20.0 / 1521.0);
    const auto fs_2_3_10 = std::sqrt(40.0 / 9.0);
    const auto fs_2_3_14 = std::sqrt(56.0 / 9.0);
    const auto fs_2_3_15 = std::sqrt(20.0 / 3.0);
    const auto fs_2_3_21 = std::sqrt(28.0 / 3.0);
    const auto fs_2_3_3 = std::sqrt(4.0 / 3.0);
    const auto fs_2_429_21 = std::sqrt(28.0 / 61347.0);
    const auto fs_315_16_2 = std::sqrt(99225.0 / 128.0);
    const auto fs_315_32_10 = std::sqrt(496125.0 / 512.0);
    const auto fs_320_7293_21 = std::sqrt(716800.0 / 17729283.0);
    const auto fs_33_16_70 = std::sqrt(38115.0 / 128.0);
    const auto fs_33_8_10 = std::sqrt(5445.0 / 32.0);
    const auto fs_35_143_14 = std::sqrt(17150.0 / 20449.0);
    const auto fs_35_22_3 = std::sqrt(3675.0 / 484.0);
    const auto fs_35_2431_143 = std::sqrt(1225.0 / 41327.0);
    const auto fs_35_2431_22 = std::sqrt(2450.0 / 537251.0);
    const auto fs_35_2431_462 = std::sqrt(51450.0 / 537251.0);
    const auto fs_35_286_22 = std::sqrt(1225.0 / 3718.0);
    const auto fs_35_429_3 = std::sqrt(1225.0 / 61347.0);
    const auto fs_35_429_6 = std::sqrt(2450.0 / 61347.0);
    const auto fs_35_4862_858 = std::sqrt(3675.0 / 82654.0);
    const auto fs_35_4_3 = std::sqrt(3675.0 / 16.0);
    const auto fs_35_572_858 = std::sqrt(3675.0 / 1144.0);
    const auto fs_35_7293_858 = std::sqrt(2450.0 / 123981.0);
    const auto fs_35_858_858 = std::sqrt(1225.0 / 858.0);
    const auto fs_3_143_21 = std::sqrt(189.0 / 20449.0);
    const auto fs_3_143_35 = std::sqrt(315.0 / 20449.0);
    const auto fs_3_2_105 = std::sqrt(945.0 / 4.0);
    const auto fs_3_4_21 = std::sqrt(189.0 / 16.0);
    const auto fs_3_8_105 = std::sqrt(945.0 / 64.0);
    const auto fs_3_8_15 = std::sqrt(135.0 / 64.0);
    const auto fs_3_8_210 = std::sqrt(945.0 / 32.0);
    const auto fs_40_2431_10 = std::sqrt(16000.0 / 5909761.0);
    const auto fs_40_2431_165 = std::sqrt(24000.0 / 537251.0);
    const auto fs_40_429_3 = std::sqrt(1600.0 / 61347.0);
    const auto fs_43_44_2 = std::sqrt(1849.0 / 968.0);
    const auto fs_43_858_2 = std::sqrt(1849.0 / 368082.0);
    const auto fs_43_8_2 = std::sqrt(1849.0 / 32.0);
    const auto fs_45_16_30 = std::sqrt(30375.0 / 128.0);
    const auto fs_45_4862_66 = std::sqrt(6075.0 / 1074502.0);
    const auto fs_45_4_7 = std::sqrt(14175.0 / 16.0);
    const auto fs_45_572_66 = std::sqrt(6075.0 / 14872.0);
    const auto fs_45_8_7 = std::sqrt(14175.0 / 64.0);
    const auto fs_490_2431_3 = std::sqrt(720300.0 / 5909761.0);
    const auto fs_49_2431_165 = std::sqrt(36015.0 / 537251.0);
    const auto fs_4_1_7 = std::sqrt(112.0);
    const auto fs_4_3_6 = std::sqrt(32.0 / 3.0);
    const auto fs_4_429_105 = std::sqrt(560.0 / 61347.0);
    const auto fs_50_143_11 = std::sqrt(2500.0 / 1859.0);
    const auto fs_57_8_6 = std::sqrt(9747.0 / 32.0);
    const auto fs_588_2431_2 = std::sqrt(691488.0 / 5909761.0);
    const auto fs_5_11_7 = std::sqrt(175.0 / 121.0);
    const auto fs_5_13_15 = std::sqrt(375.0 / 169.0);
    const auto fs_5_13_21 = std::sqrt(525.0 / 169.0);
    const auto fs_5_13_30 = std::sqrt(750.0 / 169.0);
    const auto fs_5_13_35 = std::sqrt(875.0 / 169.0);
    const auto fs_5_13_6 = std::sqrt(150.0 / 169.0);
    const auto fs_5_143_2002 = std::sqrt(350.0 / 143.0);
    const auto fs_5_143_7 = std::sqrt(175.0 / 20449.0);
    const auto fs_5_187_7 = std::sqrt(175.0 / 34969.0);
    const auto fs_5_22_105 = std::sqrt(2625.0 / 484.0);
    const auto fs_5_22_7 = std::sqrt(175.0 / 484.0);
    const auto fs_5_2431_154 = std::sqrt(350.0 / 537251.0);
    const auto fs_5_2431_70 = std::sqrt(1750.0 / 5909761.0);
    const auto fs_5_2431_858 = std::sqrt(150.0 / 41327.0);
    const auto fs_5_26_10 = std::sqrt(125.0 / 338.0);
    const auto fs_5_26_70 = std::sqrt(875.0 / 338.0);
    const auto fs_5_286_154 = std::sqrt(175.0 / 3718.0);
    const auto fs_5_286_30 = std::sqrt(375.0 / 40898.0);
    const auto fs_5_286_70 = std::sqrt(875.0 / 40898.0);
    const auto fs_5_286_858 = std::sqrt(75.0 / 286.0);
    const auto fs_5_2_5 = std::sqrt(125.0 / 4.0);
    const auto fs_5_2_7 = std::sqrt(175.0 / 4.0);
    const auto fs_5_39_105 = std::sqrt(875.0 / 507.0);
    const auto fs_5_429_105 = std::sqrt(875.0 / 61347.0);
    const auto fs_5_429_210 = std::sqrt(1750.0 / 61347.0);
    const auto fs_5_4862_30030 = std::sqrt(2625.0 / 82654.0);
    const auto fs_5_4_105 = std::sqrt(2625.0 / 16.0);
    const auto fs_5_4_15 = std::sqrt(375.0 / 16.0);
    const auto fs_5_4_21 = std::sqrt(525.0 / 16.0);
    const auto fs_5_4_30 = std::sqrt(375.0 / 8.0);
    const auto fs_5_4_35 = std::sqrt(875.0 / 16.0);
    const auto fs_5_4_6 = std::sqrt(75.0 / 8.0);
    const auto fs_5_561_33 = std::sqrt(25.0 / 9537.0);
    const auto fs_5_572_30030 = std::sqrt(2625.0 / 1144.0);
    const auto fs_5_66_33 = std::sqrt(25.0 / 132.0);
    const auto fs_5_8_10 = std::sqrt(125.0 / 32.0);
    const auto fs_5_8_70 = std::sqrt(875.0 / 32.0);
    const auto fs_63_2_2 = std::sqrt(3969.0 / 2.0);
    const auto fs_63_4_10 = std::sqrt(19845.0 / 8.0);
    const auto fs_6_1_7 = std::sqrt(252.0);
    const auto fs_70_143_11 = std::sqrt(4900.0 / 1859.0);
    const auto fs_70_2431_14 = std::sqrt(68600.0 / 5909761.0);
    const auto fs_70_429_14 = std::sqrt(68600.0 / 184041.0);
    const auto fs_70_429_165 = std::sqrt(24500.0 / 5577.0);
    const auto fs_70_7293_6 = std::sqrt(9800.0 / 17729283.0);
    const auto fs_75_16_10 = std::sqrt(28125.0 / 128.0);
    const auto fs_75_2431_66 = std::sqrt(33750.0 / 537251.0);
    const auto fs_75_286_66 = std::sqrt(16875.0 / 3718.0);
    const auto fs_75_4862_30 = std::sqrt(84375.0 / 11819522.0);
    const auto fs_75_572_30 = std::sqrt(84375.0 / 163592.0);
    const auto fs_7_11_15 = std::sqrt(735.0 / 121.0);
    const auto fs_7_11_5 = std::sqrt(245.0 / 121.0);
    const auto fs_7_1_3 = std::sqrt(147.0);
    const auto fs_7_22_3 = std::sqrt(147.0 / 484.0);
    const auto fs_7_2431_1001 = std::sqrt(343.0 / 41327.0);
    const auto fs_7_2431_1155 = std::sqrt(5145.0 / 537251.0);
    const auto fs_7_2431_12155 = std::sqrt(245.0 / 2431.0);
    const auto fs_7_2431_143 = std::sqrt(49.0 / 41327.0);
    const auto fs_7_2431_33 = std::sqrt(147.0 / 537251.0);
    const auto fs_7_2431_4290 = std::sqrt(1470.0 / 41327.0);
    const auto fs_7_2_15 = std::sqrt(735.0 / 4.0);
    const auto fs_7_2_5 = std::sqrt(245.0 / 4.0);
    const auto fs_7_429_3 = std::sqrt(49.0 / 61347.0);
    const auto fs_7_4862_2 = std::sqrt(49.0 / 11819522.0);
    const auto fs_7_4862_22 = std::sqrt(49.0 / 1074502.0);
    const auto fs_7_4862_30030 = std::sqrt(5145.0 / 82654.0);
    const auto fs_7_4862_330 = std::sqrt(735.0 / 1074502.0);
    const auto fs_7_4862_6006 = std::sqrt(1029.0 / 82654.0);
    const auto fs_7_4862_770 = std::sqrt(1715.0 / 1074502.0);
    const auto fs_7_4_3 = std::sqrt(147.0 / 16.0);
    const auto fs_815_14586_6 = std::sqrt(664225.0 / 35458566.0);
    const auto fs_815_1716_6 = std::sqrt(664225.0 / 490776.0);
    const auto fs_84_2431_21 = std::sqrt(148176.0 / 5909761.0);
    const auto fs_84_2431_33 = std::sqrt(21168.0 / 537251.0);
    const auto fs_84_2431_55 = std::sqrt(35280.0 / 537251.0);
    const auto fs_8_11_7 = std::sqrt(448.0 / 121.0);
    const auto fs_98_2431_15 = std::sqrt(144060.0 / 5909761.0);
    const auto fs_9_1_6 = std::sqrt(486.0);
    const auto fs_9_22_21 = std::sqrt(1701.0 / 484.0);
    const auto fs_9_22_35 = std::sqrt(2835.0 / 484.0);
    const auto fs_9_2_10 = std::sqrt(405.0 / 2.0);
    const auto fs_9_2_14 = std::sqrt(567.0 / 2.0);
    const auto fs_9_2_15 = std::sqrt(1215.0 / 4.0);
    const auto fs_9_2_21 = std::sqrt(1701.0 / 4.0);
    const auto fs_9_2_3 = std::sqrt(243.0 / 4.0);
    const auto fs_9_4_2 = std::sqrt(81.0 / 8.0);
    const auto fs_9_4_21 = std::sqrt(1701.0 / 16.0);
    const auto fs_9_4_35 = std::sqrt(2835.0 / 16.0);
    const auto fs_9_4_42 = std::sqrt(1701.0 / 8.0);
    const auto fs_9_4_6 = std::sqrt(243.0 / 8.0);

    // NOTE: the angular components are formed in 67 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph1_p1, ph3_p1, ph3_p2, ph5_p1, ph5_p2, ph7_p1, ph7_p2, ph9_p1, ph9_p2, ph9_p8, ph9_p9, ab_2, pc_0, pc_1 : simd::cache_line_size())
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

        const auto h1_p1 = ph1_p1[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h3_p2 = ph3_p2[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h9_p1 = ph9_p1[k];
        const auto h9_p2 = ph9_p2[k];
        const auto h9_p8 = ph9_p8[k];
        const auto h9_p9 = ph9_p9[k];

        pc_0[k] = - e_0 * fs_315_32_10 * h1_p1 - e_1 * fs_21_4_15 * h3_p1 + e_1 * fs_63_4_10 * r_2 * h1_p1 - e_2 * fs_5_4_6 * h5_p1 + e_2 * fs_7_2_15 * r_2 * h3_p1 - e_2 * fs_27_4_10 * r_4 * h1_p1 - e_3 * fs_5_286_70 * h7_p1 + e_3 * fs_5_13_6 * r_2 * h5_p1 - e_3 * fs_7_11_15 * r_4 * h3_p1 + e_3 * fs_1_1_10 * r_6 * h1_p1 - e_4 * fs_7_4862_2 * h9_p1 + e_4 * fs_21_2431_2431 * h9_p9 + e_4 * fs_5_2431_70 * r_2 * h7_p1 - e_4 * fs_1_39_6 * r_4 * h5_p1 + e_4 * fs_14_429_15 * r_6 * h3_p1 - e_4 * fs_1_22_10 * r_8 * h1_p1;

        pc_1[k] = e_1 * fs_105_8_3 * h3_p2 + e_2 * fs_5_4_21 * h5_p2 - e_2 * fs_35_4_3 * r_2 * h3_p2 + e_3 * fs_15_572_210 * h7_p2 - e_3 * fs_5_13_21 * r_2 * h5_p2 + e_3 * fs_35_22_3 * r_4 * h3_p2 + e_4 * fs_7_4862_22 * h9_p2 + e_4 * fs_14_2431_2431 * h9_p8 - e_4 * fs_15_4862_210 * r_2 * h7_p2 + e_4 * fs_1_39_21 * r_4 * h5_p2 - e_4 * fs_35_429_3 * r_6 * h3_p2;
    }

    // NOTE: the angular components are formed in 67 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, ph3_p3, ph5_p3, ph5_p4, ph7_p3, ph7_p4, ph7_p6, ph7_p7, ph9_p3, ph9_p4, ph9_p6, ph9_p7, ab_2, pc_2, pc_3 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h3_p3 = ph3_p3[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h5_p4 = ph5_p4[k];
        const auto h7_p3 = ph7_p3[k];
        const auto h7_p4 = ph7_p4[k];
        const auto h7_p6 = ph7_p6[k];
        const auto h7_p7 = ph7_p7[k];
        const auto h9_p3 = ph9_p3[k];
        const auto h9_p4 = ph9_p4[k];
        const auto h9_p6 = ph9_p6[k];
        const auto h9_p7 = ph9_p7[k];

        pc_2[k] = - e_1 * fs_45_8_7 * h3_p3 - e_2 * f_15_2 * h5_p3 + e_2 * fs_15_4_7 * r_2 * h3_p3 - e_3 * fs_75_572_30 * h7_p3 + e_3 * fs_5_572_30030 * h7_p7 + e_3 * f_30_13 * r_2 * h5_p3 - e_3 * fs_15_22_7 * r_4 * h3_p3 - e_4 * fs_7_2431_33 * h9_p3 + e_4 * fs_14_2431_1001 * h9_p7 + e_4 * fs_75_4862_30 * r_2 * h7_p3 - e_4 * fs_5_4862_30030 * r_2 * h7_p7 - e_4 * f_2_13 * r_4 * h5_p3 + e_4 * fs_5_143_7 * r_6 * h3_p3;

        pc_3[k] = e_2 * f_15_2 * h5_p4 + e_3 * fs_25_286_165 * h7_p4 + e_3 * fs_15_572_4290 * h7_p6 - e_3 * f_30_13 * r_2 * h5_p4 + e_4 * fs_7_2431_143 * h9_p4 + e_4 * fs_7_4862_6006 * h9_p6 - e_4 * fs_25_2431_165 * r_2 * h7_p4 - e_4 * fs_15_4862_4290 * r_2 * h7_p6 + e_4 * f_2_13 * r_4 * h5_p4;
    }

    // NOTE: the angular components are formed in 67 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_2, pe_3, pe_4, ph5_m5, ph5_m4, ph7_m6, ph7_m5, ph7_m4, ph9_m6, ph9_m5, ph9_m4, ab_2, pc_4, pc_5 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

        const auto h5_m5 = ph5_m5[k];
        const auto h5_m4 = ph5_m4[k];
        const auto h7_m6 = ph7_m6[k];
        const auto h7_m5 = ph7_m5[k];
        const auto h7_m4 = ph7_m4[k];
        const auto h9_m6 = ph9_m6[k];
        const auto h9_m5 = ph9_m5[k];
        const auto h9_m4 = ph9_m4[k];

        pc_4[k] = - e_2 * f_15_2 * h5_m5 - e_3 * fs_75_286_66 * h7_m5 + e_3 * f_30_13 * r_2 * h5_m5 - e_4 * fs_7_2431_1001 * h9_m5 + e_4 * fs_75_2431_66 * r_2 * h7_m5 - e_4 * f_2_13 * r_4 * h5_m5;

        pc_5[k] = e_2 * f_15_2 * h5_m4 - e_3 * fs_15_572_4290 * h7_m6 + e_3 * fs_25_286_165 * h7_m4 - e_3 * f_30_13 * r_2 * h5_m4 - e_4 * fs_7_4862_6006 * h9_m6 + e_4 * fs_7_2431_143 * h9_m4 + e_4 * fs_15_4862_4290 * r_2 * h7_m6 - e_4 * fs_25_2431_165 * r_2 * h7_m4 + e_4 * f_2_13 * r_4 * h5_m4;
    }

    // NOTE: the angular components are formed in 67 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, ph3_m3, ph3_m2, ph5_m3, ph5_m2, ph7_m7, ph7_m3, ph7_m2, ph9_m8, ph9_m7, ph9_m3, ph9_m2, ab_2, pc_6, pc_7 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h3_m3 = ph3_m3[k];
        const auto h3_m2 = ph3_m2[k];
        const auto h5_m3 = ph5_m3[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h7_m7 = ph7_m7[k];
        const auto h7_m3 = ph7_m3[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h9_m8 = ph9_m8[k];
        const auto h9_m7 = ph9_m7[k];
        const auto h9_m3 = ph9_m3[k];
        const auto h9_m2 = ph9_m2[k];

        pc_6[k] = - e_1 * fs_45_8_7 * h3_m3 - e_2 * f_15_2 * h5_m3 + e_2 * fs_15_4_7 * r_2 * h3_m3 - e_3 * fs_5_572_30030 * h7_m7 - e_3 * fs_75_572_30 * h7_m3 + e_3 * f_30_13 * r_2 * h5_m3 - e_3 * fs_15_22_7 * r_4 * h3_m3 - e_4 * fs_14_2431_1001 * h9_m7 - e_4 * fs_7_2431_33 * h9_m3 + e_4 * fs_5_4862_30030 * r_2 * h7_m7 + e_4 * fs_75_4862_30 * r_2 * h7_m3 - e_4 * f_2_13 * r_4 * h5_m3 + e_4 * fs_5_143_7 * r_6 * h3_m3;

        pc_7[k] = e_1 * fs_105_8_3 * h3_m2 + e_2 * fs_5_4_21 * h5_m2 - e_2 * fs_35_4_3 * r_2 * h3_m2 + e_3 * fs_15_572_210 * h7_m2 - e_3 * fs_5_13_21 * r_2 * h5_m2 + e_3 * fs_35_22_3 * r_4 * h3_m2 - e_4 * fs_14_2431_2431 * h9_m8 + e_4 * fs_7_4862_22 * h9_m2 - e_4 * fs_15_4862_210 * r_2 * h7_m2 + e_4 * fs_1_39_21 * r_4 * h5_m2 - e_4 * fs_35_429_3 * r_6 * h3_m2;
    }

    // NOTE: the angular components are formed in 67 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph1_m1, ph1_0, ph3_m1, ph3_0, ph5_m1, ph5_0, ph7_m1, ph7_0, ph9_m9, ph9_m1, ph9_0, ph9_p8, ab_2, pc_8, pc_9 : simd::cache_line_size())
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

        const auto h1_m1 = ph1_m1[k];
        const auto h1_0 = ph1_0[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h3_0 = ph3_0[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h5_0 = ph5_0[k];
        const auto h7_m1 = ph7_m1[k];
        const auto h7_0 = ph7_0[k];
        const auto h9_m9 = ph9_m9[k];
        const auto h9_m1 = ph9_m1[k];
        const auto h9_0 = ph9_0[k];
        const auto h9_p8 = ph9_p8[k];

        pc_8[k] = - e_0 * fs_315_32_10 * h1_m1 - e_1 * fs_21_4_15 * h3_m1 + e_1 * fs_63_4_10 * r_2 * h1_m1 - e_2 * fs_5_4_6 * h5_m1 + e_2 * fs_7_2_15 * r_2 * h3_m1 - e_2 * fs_27_4_10 * r_4 * h1_m1 - e_3 * fs_5_286_70 * h7_m1 + e_3 * fs_5_13_6 * r_2 * h5_m1 - e_3 * fs_7_11_15 * r_4 * h3_m1 + e_3 * fs_1_1_10 * r_6 * h1_m1 - e_4 * fs_21_2431_2431 * h9_m9 - e_4 * fs_7_4862_2 * h9_m1 + e_4 * fs_5_2431_70 * r_2 * h7_m1 - e_4 * fs_1_39_6 * r_4 * h5_m1 + e_4 * fs_14_429_15 * r_6 * h3_m1 - e_4 * fs_1_22_10 * r_8 * h1_m1;

        pc_9[k] = - e_0 * f_315_16 * h1_0 - e_1 * f_63_2 * h3_0 + e_1 * f_63_2 * r_2 * h1_0 - e_2 * f_15_2 * h5_0 + e_2 * f_21_1 * r_2 * h3_0 - e_2 * f_27_2 * r_4 * h1_0 - e_3 * f_70_143 * h7_0 + e_3 * f_30_13 * r_2 * h5_0 - e_3 * f_42_11 * r_4 * h3_0 + e_3 * f_2_1 * r_6 * h1_0 - e_4 * f_21_2431 * h9_0 + e_4 * fs_7_2431_12155 * h9_p8 + e_4 * f_140_2431 * r_2 * h7_0 - e_4 * f_2_13 * r_4 * h5_0 + e_4 * f_28_143 * r_6 * h3_0 - e_4 * f_1_11 * r_8 * h1_0;
    }

    // NOTE: the angular components are formed in 67 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph1_p1, ph3_p1, ph5_p1, ph7_p1, ph7_p7, ph9_p1, ph9_p7, ab_2, pc_10 : simd::cache_line_size())
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

        const auto h1_p1 = ph1_p1[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h7_p7 = ph7_p7[k];
        const auto h9_p1 = ph9_p1[k];
        const auto h9_p7 = ph9_p7[k];

        pc_10[k] = - e_0 * fs_315_16_2 * h1_p1 + e_1 * fs_21_8_3 * h3_p1 + e_1 * fs_63_2_2 * r_2 * h1_p1 + e_2 * fs_5_4_30 * h5_p1 - e_2 * fs_7_4_3 * r_2 * h3_p1 - e_2 * fs_27_2_2 * r_4 * h1_p1 + e_3 * fs_115_572_14 * h7_p1 - e_3 * fs_35_572_858 * h7_p7 - e_3 * fs_5_13_30 * r_2 * h5_p1 + e_3 * fs_7_22_3 * r_4 * h3_p1 + e_3 * fs_2_1_2 * r_6 * h1_p1 + e_4 * fs_14_2431_10 * h9_p1 + e_4 * fs_28_2431_715 * h9_p7 - e_4 * fs_115_4862_14 * r_2 * h7_p1 + e_4 * fs_35_4862_858 * r_2 * h7_p7 + e_4 * fs_1_39_30 * r_4 * h5_p1 - e_4 * fs_7_429_3 * r_6 * h3_p1 - e_4 * fs_1_11_2 * r_8 * h1_p1;
    }

    // NOTE: the angular components are formed in 67 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, ph3_p2, ph3_p3, ph5_p2, ph5_p5, ph7_p2, ph7_p3, ph7_p5, ph7_p6, ph9_p2, ph9_p3, ph9_p5, ph9_p6, ab_2, pc_11, pc_12 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h3_p2 = ph3_p2[k];
        const auto h3_p3 = ph3_p3[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h5_p5 = ph5_p5[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h7_p3 = ph7_p3[k];
        const auto h7_p5 = ph7_p5[k];
        const auto h7_p6 = ph7_p6[k];
        const auto h9_p2 = ph9_p2[k];
        const auto h9_p3 = ph9_p3[k];
        const auto h9_p5 = ph9_p5[k];
        const auto h9_p6 = ph9_p6[k];

        pc_11[k] = e_1 * fs_3_2_105 * h3_p2 - e_2 * fs_5_4_15 * h5_p2 - e_2 * fs_1_1_105 * r_2 * h3_p2 - e_3 * fs_135_286_6 * h7_p2 - e_3 * fs_5_286_858 * h7_p6 + e_3 * fs_5_13_15 * r_2 * h5_p2 + e_3 * fs_2_11_105 * r_4 * h3_p2 - e_4 * fs_7_4862_770 * h9_p2 + e_4 * fs_7_4862_30030 * h9_p6 + e_4 * fs_135_2431_6 * r_2 * h7_p2 + e_4 * fs_5_2431_858 * r_2 * h7_p6 - e_4 * fs_1_39_15 * r_4 * h5_p2 - e_4 * fs_4_429_105 * r_6 * h3_p2;

        pc_12[k] = - e_1 * fs_27_8_35 * h3_p3 - e_2 * f_15_2 * h5_p5 + e_2 * fs_9_4_35 * r_2 * h3_p3 + e_3 * fs_25_44_6 * h7_p3 + e_3 * fs_45_572_66 * h7_p5 + e_3 * f_30_13 * r_2 * h5_p5 - e_3 * fs_9_22_35 * r_4 * h3_p3 + e_4 * fs_14_2431_165 * h9_p3 + e_4 * fs_14_2431_1001 * h9_p5 - e_4 * fs_25_374_6 * r_2 * h7_p3 - e_4 * fs_45_4862_66 * r_2 * h7_p5 - e_4 * f_2_13 * r_4 * h5_p5 + e_4 * fs_3_143_35 * r_6 * h3_p3;
    }

    // NOTE: the angular components are formed in 67 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, ph3_m3, ph5_m5, ph5_m4, ph7_m5, ph7_m4, ph7_m3, ph9_m5, ph9_m4, ph9_m3, ab_2, pc_13, pc_14 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h3_m3 = ph3_m3[k];
        const auto h5_m5 = ph5_m5[k];
        const auto h5_m4 = ph5_m4[k];
        const auto h7_m5 = ph7_m5[k];
        const auto h7_m4 = ph7_m4[k];
        const auto h7_m3 = ph7_m3[k];
        const auto h9_m5 = ph9_m5[k];
        const auto h9_m4 = ph9_m4[k];
        const auto h9_m3 = ph9_m3[k];

        pc_13[k] = e_2 * f_15_2 * h5_m4 - e_3 * fs_20_143_165 * h7_m4 - e_3 * f_30_13 * r_2 * h5_m4 - e_4 * fs_35_2431_143 * h9_m4 + e_4 * fs_40_2431_165 * r_2 * h7_m4 + e_4 * f_2_13 * r_4 * h5_m4;

        pc_14[k] = - e_1 * fs_27_8_35 * h3_m3 + e_2 * f_15_2 * h5_m5 + e_2 * fs_9_4_35 * r_2 * h3_m3 - e_3 * fs_45_572_66 * h7_m5 + e_3 * fs_25_44_6 * h7_m3 - e_3 * f_30_13 * r_2 * h5_m5 - e_3 * fs_9_22_35 * r_4 * h3_m3 - e_4 * fs_14_2431_1001 * h9_m5 + e_4 * fs_14_2431_165 * h9_m3 + e_4 * fs_45_4862_66 * r_2 * h7_m5 - e_4 * fs_25_374_6 * r_2 * h7_m3 + e_4 * f_2_13 * r_4 * h5_m5 + e_4 * fs_3_143_35 * r_6 * h3_m3;
    }

    // NOTE: the angular components are formed in 67 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, ph3_m2, ph5_m2, ph7_m6, ph7_m2, ph9_m6, ph9_m2, ab_2, pc_15 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h3_m2 = ph3_m2[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h7_m6 = ph7_m6[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h9_m6 = ph9_m6[k];
        const auto h9_m2 = ph9_m2[k];

        pc_15[k] = e_1 * fs_3_2_105 * h3_m2 - e_2 * fs_5_4_15 * h5_m2 - e_2 * fs_1_1_105 * r_2 * h3_m2 + e_3 * fs_5_286_858 * h7_m6 - e_3 * fs_135_286_6 * h7_m2 + e_3 * fs_5_13_15 * r_2 * h5_m2 + e_3 * fs_2_11_105 * r_4 * h3_m2 - e_4 * fs_7_4862_30030 * h9_m6 - e_4 * fs_7_4862_770 * h9_m2 - e_4 * fs_5_2431_858 * r_2 * h7_m6 + e_4 * fs_135_2431_6 * r_2 * h7_m2 - e_4 * fs_1_39_15 * r_4 * h5_m2 - e_4 * fs_4_429_105 * r_6 * h3_m2;
    }

    // NOTE: the angular components are formed in 67 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph1_m1, ph3_m1, ph5_m1, ph7_m7, ph7_m1, ph9_m8, ph9_m7, ph9_m1, ab_2, pc_16, pc_17 : simd::cache_line_size())
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

        const auto h1_m1 = ph1_m1[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h7_m7 = ph7_m7[k];
        const auto h7_m1 = ph7_m1[k];
        const auto h9_m8 = ph9_m8[k];
        const auto h9_m7 = ph9_m7[k];
        const auto h9_m1 = ph9_m1[k];

        pc_16[k] = - e_0 * fs_315_16_2 * h1_m1 + e_1 * fs_21_8_3 * h3_m1 + e_1 * fs_63_2_2 * r_2 * h1_m1 + e_2 * fs_5_4_30 * h5_m1 - e_2 * fs_7_4_3 * r_2 * h3_m1 - e_2 * fs_27_2_2 * r_4 * h1_m1 + e_3 * fs_35_572_858 * h7_m7 + e_3 * fs_115_572_14 * h7_m1 - e_3 * fs_5_13_30 * r_2 * h5_m1 + e_3 * fs_7_22_3 * r_4 * h3_m1 + e_3 * fs_2_1_2 * r_6 * h1_m1 - e_4 * fs_28_2431_715 * h9_m7 + e_4 * fs_14_2431_10 * h9_m1 - e_4 * fs_35_4862_858 * r_2 * h7_m7 - e_4 * fs_115_4862_14 * r_2 * h7_m1 + e_4 * fs_1_39_30 * r_4 * h5_m1 - e_4 * fs_7_429_3 * r_6 * h3_m1 - e_4 * fs_1_11_2 * r_8 * h1_m1;

        pc_17[k] = - e_4 * fs_7_2431_12155 * h9_m8;
    }

    // NOTE: the angular components are formed in 67 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph1_p1, ph3_p1, ph5_p1, ph7_p1, ph7_p7, ph9_p1, ph9_p7, ab_2, pc_18 : simd::cache_line_size())
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

        const auto h1_p1 = ph1_p1[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h7_p7 = ph7_p7[k];
        const auto h9_p1 = ph9_p1[k];
        const auto h9_p7 = ph9_p7[k];

        pc_18[k] = e_0 * fs_105_32_2 * h1_p1 + e_1 * fs_21_2_3 * h3_p1 - e_1 * fs_21_4_2 * r_2 * h1_p1 + e_2 * fs_5_4_30 * h5_p1 - e_2 * fs_7_1_3 * r_2 * h3_p1 + e_2 * fs_9_4_2 * r_4 * h1_p1 + e_3 * fs_70_429_14 * h7_p1 + e_3 * fs_35_858_858 * h7_p7 - e_3 * fs_5_13_30 * r_2 * h5_p1 + e_3 * fs_14_11_3 * r_4 * h3_p1 - e_3 * fs_1_3_2 * r_6 * h1_p1 + e_4 * fs_21_4862_10 * h9_p1 + e_4 * fs_21_2431_715 * h9_p7 - e_4 * fs_140_7293_14 * r_2 * h7_p1 - e_4 * fs_35_7293_858 * r_2 * h7_p7 + e_4 * fs_1_39_30 * r_4 * h5_p1 - e_4 * fs_28_429_3 * r_6 * h3_p1 + e_4 * fs_1_66_2 * r_8 * h1_p1;
    }

    // NOTE: the angular components are formed in 67 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph1_0, ph3_0, ph5_0, ph7_0, ph7_p6, ph9_0, ph9_p6, ab_2, pc_19 : simd::cache_line_size())
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

        const auto h1_0 = ph1_0[k];
        const auto h3_0 = ph3_0[k];
        const auto h5_0 = ph5_0[k];
        const auto h7_0 = ph7_0[k];
        const auto h7_p6 = ph7_p6[k];
        const auto h9_0 = ph9_0[k];
        const auto h9_p6 = ph9_p6[k];

        pc_19[k] = - e_0 * f_105_4 * h1_0 - e_1 * f_63_4 * h3_0 + e_1 * f_42_1 * r_2 * h1_0 + e_2 * f_15_2 * h5_0 + e_2 * f_21_2 * r_2 * h3_0 - e_2 * f_18_1 * r_4 * h1_0 + e_3 * f_665_429 * h7_0 - e_3 * fs_25_1716_6006 * h7_p6 - e_3 * f_30_13 * r_2 * h5_0 - e_3 * f_21_11 * r_4 * h3_0 + e_3 * f_8_3 * r_6 * h1_0 + e_4 * f_126_2431 * h9_0 + e_4 * fs_21_4862_4290 * h9_p6 - e_4 * f_1330_7293 * r_2 * h7_0 + e_4 * fs_25_14586_6006 * r_2 * h7_p6 + e_4 * f_2_13 * r_4 * h5_0 + e_4 * f_14_143 * r_6 * h3_0 - e_4 * f_4_33 * r_8 * h1_0;
    }

    // NOTE: the angular components are formed in 67 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph1_p1, ph3_p1, ph5_p5, ph7_p1, ph7_p5, ph9_p1, ph9_p5, ab_2, pc_20 : simd::cache_line_size())
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

        const auto h1_p1 = ph1_p1[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h5_p5 = ph5_p5[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h7_p5 = ph7_p5[k];
        const auto h9_p1 = ph9_p1[k];
        const auto h9_p5 = ph9_p5[k];

        pc_20[k] = - e_0 * fs_105_16_14 * h1_p1 + e_1 * fs_27_8_21 * h3_p1 + e_1 * fs_21_2_14 * r_2 * h1_p1 + e_2 * f_15_2 * h5_p5 - e_2 * fs_9_4_21 * r_2 * h3_p1 - e_2 * fs_9_2_14 * r_4 * h1_p1 - e_3 * fs_125_132_2 * h7_p1 - e_3 * fs_265_1716_66 * h7_p5 - e_3 * f_30_13 * r_2 * h5_p5 + e_3 * fs_9_22_21 * r_4 * h3_p1 + e_3 * fs_2_3_14 * r_6 * h1_p1 - e_4 * fs_21_2431_70 * h9_p1 + e_4 * fs_21_2431_1001 * h9_p5 + e_4 * fs_125_1122_2 * r_2 * h7_p1 + e_4 * fs_265_14586_66 * r_2 * h7_p5 + e_4 * f_2_13 * r_4 * h5_p5 - e_4 * fs_3_143_21 * r_6 * h3_p1 - e_4 * fs_1_33_14 * r_8 * h1_p1;
    }

    // NOTE: the angular components are formed in 67 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, ph3_m3, ph3_p2, ph5_m3, ph5_p2, ph7_m3, ph7_p2, ph7_p4, ph9_m3, ph9_p2, ph9_p4, ab_2, pc_21, pc_22 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h3_m3 = ph3_m3[k];
        const auto h3_p2 = ph3_p2[k];
        const auto h5_m3 = ph5_m3[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h7_m3 = ph7_m3[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h7_p4 = ph7_p4[k];
        const auto h9_m3 = ph9_m3[k];
        const auto h9_p2 = ph9_p2[k];
        const auto h9_p4 = ph9_p4[k];

        pc_21[k] = e_1 * fs_3_8_105 * h3_p2 - e_2 * fs_5_4_15 * h5_p2 - e_2 * fs_1_4_105 * r_2 * h3_p2 + e_3 * fs_815_1716_6 * h7_p2 - e_3 * fs_5_66_33 * h7_p4 + e_3 * fs_5_13_15 * r_2 * h5_p2 + e_3 * fs_1_22_105 * r_4 * h3_p2 + e_4 * fs_21_4862_770 * h9_p2 + e_4 * fs_21_2431_715 * h9_p4 - e_4 * fs_815_14586_6 * r_2 * h7_p2 + e_4 * fs_5_561_33 * r_2 * h7_p4 - e_4 * fs_1_39_15 * r_4 * h5_p2 - e_4 * fs_1_429_105 * r_6 * h3_p2;

        pc_22[k] = - e_1 * fs_45_4_7 * h3_m3 + e_2 * f_15_2 * h5_m3 + e_2 * fs_15_2_7 * r_2 * h3_m3 - e_3 * fs_115_858_30 * h7_m3 - e_3 * f_30_13 * r_2 * h5_m3 - e_3 * fs_15_11_7 * r_4 * h3_m3 - e_4 * fs_105_2431_33 * h9_m3 + e_4 * fs_115_7293_30 * r_2 * h7_m3 + e_4 * f_2_13 * r_4 * h5_m3 + e_4 * fs_10_143_7 * r_6 * h3_m3;
    }

    // NOTE: the angular components are formed in 67 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, ph3_m2, ph5_m2, ph7_m4, ph7_m2, ph9_m4, ph9_m2, ab_2, pc_23 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h3_m2 = ph3_m2[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h7_m4 = ph7_m4[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h9_m4 = ph9_m4[k];
        const auto h9_m2 = ph9_m2[k];

        pc_23[k] = e_1 * fs_3_8_105 * h3_m2 - e_2 * fs_5_4_15 * h5_m2 - e_2 * fs_1_4_105 * r_2 * h3_m2 + e_3 * fs_5_66_33 * h7_m4 + e_3 * fs_815_1716_6 * h7_m2 + e_3 * fs_5_13_15 * r_2 * h5_m2 + e_3 * fs_1_22_105 * r_4 * h3_m2 - e_4 * fs_21_2431_715 * h9_m4 + e_4 * fs_21_4862_770 * h9_m2 - e_4 * fs_5_561_33 * r_2 * h7_m4 - e_4 * fs_815_14586_6 * r_2 * h7_m2 - e_4 * fs_1_39_15 * r_4 * h5_m2 - e_4 * fs_1_429_105 * r_6 * h3_m2;
    }

    // NOTE: the angular components are formed in 67 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph1_m1, ph3_m1, ph5_m5, ph7_m6, ph7_m5, ph7_m1, ph9_m6, ph9_m5, ph9_m1, ab_2, pc_24, pc_25 : simd::cache_line_size())
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

        const auto h1_m1 = ph1_m1[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h5_m5 = ph5_m5[k];
        const auto h7_m6 = ph7_m6[k];
        const auto h7_m5 = ph7_m5[k];
        const auto h7_m1 = ph7_m1[k];
        const auto h9_m6 = ph9_m6[k];
        const auto h9_m5 = ph9_m5[k];
        const auto h9_m1 = ph9_m1[k];

        pc_24[k] = - e_0 * fs_105_16_14 * h1_m1 + e_1 * fs_27_8_21 * h3_m1 + e_1 * fs_21_2_14 * r_2 * h1_m1 - e_2 * f_15_2 * h5_m5 - e_2 * fs_9_4_21 * r_2 * h3_m1 - e_2 * fs_9_2_14 * r_4 * h1_m1 + e_3 * fs_265_1716_66 * h7_m5 - e_3 * fs_125_132_2 * h7_m1 + e_3 * f_30_13 * r_2 * h5_m5 + e_3 * fs_9_22_21 * r_4 * h3_m1 + e_3 * fs_2_3_14 * r_6 * h1_m1 - e_4 * fs_21_2431_1001 * h9_m5 - e_4 * fs_21_2431_70 * h9_m1 - e_4 * fs_265_14586_66 * r_2 * h7_m5 + e_4 * fs_125_1122_2 * r_2 * h7_m1 - e_4 * f_2_13 * r_4 * h5_m5 - e_4 * fs_3_143_21 * r_6 * h3_m1 - e_4 * fs_1_33_14 * r_8 * h1_m1;

        pc_25[k] = e_3 * fs_25_1716_6006 * h7_m6 - e_4 * fs_21_4862_4290 * h9_m6 - e_4 * fs_25_14586_6006 * r_2 * h7_m6;
    }

    // NOTE: the angular components are formed in 67 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph1_m1, ph3_m1, ph5_m1, ph7_m7, ph7_m1, ph9_m7, ph9_m1, ab_2, pc_26 : simd::cache_line_size())
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

        const auto h1_m1 = ph1_m1[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h7_m7 = ph7_m7[k];
        const auto h7_m1 = ph7_m1[k];
        const auto h9_m7 = ph9_m7[k];
        const auto h9_m1 = ph9_m1[k];

        pc_26[k] = - e_0 * fs_105_32_2 * h1_m1 - e_1 * fs_21_2_3 * h3_m1 + e_1 * fs_21_4_2 * r_2 * h1_m1 - e_2 * fs_5_4_30 * h5_m1 + e_2 * fs_7_1_3 * r_2 * h3_m1 - e_2 * fs_9_4_2 * r_4 * h1_m1 - e_3 * fs_35_858_858 * h7_m7 - e_3 * fs_70_429_14 * h7_m1 + e_3 * fs_5_13_30 * r_2 * h5_m1 - e_3 * fs_14_11_3 * r_4 * h3_m1 + e_3 * fs_1_3_2 * r_6 * h1_m1 - e_4 * fs_21_2431_715 * h9_m7 - e_4 * fs_21_4862_10 * h9_m1 + e_4 * fs_35_7293_858 * r_2 * h7_m7 + e_4 * fs_140_7293_14 * r_2 * h7_m1 - e_4 * fs_1_39_30 * r_4 * h5_m1 + e_4 * fs_28_429_3 * r_6 * h3_m1 - e_4 * fs_1_66_2 * r_8 * h1_m1;
    }

    // NOTE: the angular components are formed in 67 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, ph3_p2, ph5_p2, ph7_p2, ph7_p6, ph9_p2, ph9_p6, ab_2, pc_27 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h3_p2 = ph3_p2[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h7_p6 = ph7_p6[k];
        const auto h9_p2 = ph9_p2[k];
        const auto h9_p6 = ph9_p6[k];

        pc_27[k] = - e_1 * fs_21_4_5 * h3_p2 - e_2 * fs_5_4_35 * h5_p2 + e_2 * fs_7_2_5 * r_2 * h3_p2 - e_3 * fs_35_143_14 * h7_p2 + e_3 * fs_5_143_2002 * h7_p6 + e_3 * fs_5_13_35 * r_2 * h5_p2 - e_3 * fs_7_11_5 * r_4 * h3_p2 - e_4 * fs_7_4862_330 * h9_p2 + e_4 * fs_21_4862_1430 * h9_p6 + e_4 * fs_70_2431_14 * r_2 * h7_p2 - e_4 * fs_10_2431_2002 * r_2 * h7_p6 - e_4 * fs_1_39_35 * r_4 * h5_p2 + e_4 * fs_14_429_5 * r_6 * h3_p2;
    }

    // NOTE: the angular components are formed in 67 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph1_p1, ph3_p1, ph5_p1, ph5_p5, ph7_p1, ph7_p5, ph9_p1, ph9_p5, ab_2, pc_28 : simd::cache_line_size())
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

        const auto h1_p1 = ph1_p1[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h5_p5 = ph5_p5[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h7_p5 = ph7_p5[k];
        const auto h9_p1 = ph9_p1[k];
        const auto h9_p5 = ph9_p5[k];

        pc_28[k] = e_0 * fs_105_32_6 * h1_p1 + e_1 * f_147_8 * h3_p1 - e_1 * fs_21_4_6 * r_2 * h1_p1 - e_2 * fs_5_8_10 * h5_p1 - e_2 * fs_5_4_21 * h5_p5 - e_2 * f_49_4 * r_2 * h3_p1 + e_2 * fs_9_4_6 * r_4 * h1_p1 - e_3 * fs_175_858_42 * h7_p1 - e_3 * fs_5_286_154 * h7_p5 + e_3 * fs_5_26_10 * r_2 * h5_p1 + e_3 * fs_5_13_21 * r_2 * h5_p5 + e_3 * f_49_22 * r_4 * h3_p1 - e_3 * fs_1_3_6 * r_6 * h1_p1 - e_4 * fs_28_2431_30 * h9_p1 + e_4 * fs_28_2431_429 * h9_p5 + e_4 * fs_175_7293_42 * r_2 * h7_p1 + e_4 * fs_5_2431_154 * r_2 * h7_p5 - e_4 * fs_1_78_10 * r_4 * h5_p1 - e_4 * fs_1_39_21 * r_4 * h5_p5 - e_4 * f_49_429 * r_6 * h3_p1 + e_4 * fs_1_66_6 * r_8 * h1_p1;
    }

    // NOTE: the angular components are formed in 67 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph1_0, ph3_0, ph5_0, ph5_p4, ph7_0, ph7_p4, ph9_0, ph9_p4, ab_2, pc_29 : simd::cache_line_size())
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

        const auto h1_0 = ph1_0[k];
        const auto h3_0 = ph3_0[k];
        const auto h5_0 = ph5_0[k];
        const auto h5_p4 = ph5_p4[k];
        const auto h7_0 = ph7_0[k];
        const auto h7_p4 = ph7_p4[k];
        const auto h9_0 = ph9_0[k];
        const auto h9_p4 = ph9_p4[k];

        pc_29[k] = - e_0 * fs_105_16_21 * h1_0 + e_1 * fs_3_4_21 * h3_0 + e_1 * fs_21_2_21 * r_2 * h1_0 + e_2 * fs_5_4_21 * h5_0 + e_2 * fs_5_4_15 * h5_p4 - e_2 * fs_1_2_21 * r_2 * h3_0 - e_2 * fs_9_2_21 * r_4 * h1_0 - e_3 * fs_160_429_21 * h7_0 - e_3 * fs_50_143_11 * h7_p4 - e_3 * fs_5_13_21 * r_2 * h5_0 - e_3 * fs_5_13_15 * r_2 * h5_p4 + e_3 * fs_1_11_21 * r_4 * h3_0 + e_3 * fs_2_3_21 * r_6 * h1_0 - e_4 * fs_84_2431_21 * h9_0 + e_4 * fs_14_2431_2145 * h9_p4 + e_4 * fs_320_7293_21 * r_2 * h7_0 + e_4 * fs_100_2431_11 * r_2 * h7_p4 + e_4 * fs_1_39_21 * r_4 * h5_0 + e_4 * fs_1_39_15 * r_4 * h5_p4 - e_4 * fs_2_429_21 * r_6 * h3_0 - e_4 * fs_1_33_21 * r_8 * h1_0;
    }

    // NOTE: the angular components are formed in 67 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph1_p1, ph3_p1, ph3_p3, ph5_p1, ph5_p3, ph7_p1, ph7_p3, ph9_p1, ph9_p3, ab_2, pc_30 : simd::cache_line_size())
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

        const auto h1_p1 = ph1_p1[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h3_p3 = ph3_p3[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h7_p3 = ph7_p3[k];
        const auto h9_p1 = ph9_p1[k];
        const auto h9_p3 = ph9_p3[k];

        pc_30[k] = - e_0 * fs_105_32_42 * h1_p1 + e_1 * fs_6_1_7 * h3_p1 - e_1 * fs_15_8_105 * h3_p3 + e_1 * fs_21_4_42 * r_2 * h1_p1 - e_2 * fs_5_8_70 * h5_p1 + e_2 * fs_5_4_15 * h5_p3 - e_2 * fs_4_1_7 * r_2 * h3_p1 + e_2 * fs_5_4_105 * r_2 * h3_p3 - e_2 * fs_9_4_42 * r_4 * h1_p1 + e_3 * fs_205_858_6 * h7_p1 - e_3 * fs_215_286_2 * h7_p3 + e_3 * fs_5_26_70 * r_2 * h5_p1 - e_3 * fs_5_13_15 * r_2 * h5_p3 + e_3 * fs_8_11_7 * r_4 * h3_p1 - e_3 * fs_5_22_105 * r_4 * h3_p3 + e_3 * fs_1_3_42 * r_6 * h1_p1 + e_4 * fs_28_2431_210 * h9_p1 + e_4 * fs_84_2431_55 * h9_p3 - e_4 * fs_205_7293_6 * r_2 * h7_p1 + e_4 * fs_215_2431_2 * r_2 * h7_p3 - e_4 * fs_1_78_70 * r_4 * h5_p1 + e_4 * fs_1_39_15 * r_4 * h5_p3 - e_4 * fs_16_429_7 * r_6 * h3_p1 + e_4 * fs_5_429_105 * r_6 * h3_p3 - e_4 * fs_1_66_42 * r_8 * h1_p1;
    }

    // NOTE: the angular components are formed in 67 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, ph3_m2, ph5_m2, ph7_m2, ph9_m2, ab_2, pc_31 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h3_m2 = ph3_m2[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h9_m2 = ph9_m2[k];

        pc_31[k] = - e_1 * fs_15_4_7 * h3_m2 + e_2 * f_5_4 * h5_m2 + e_2 * fs_5_2_7 * r_2 * h3_m2 + e_3 * fs_20_143_10 * h7_m2 - e_3 * f_5_13 * r_2 * h5_m2 - e_3 * fs_5_11_7 * r_4 * h3_m2 - e_4 * fs_35_2431_462 * h9_m2 - e_4 * fs_40_2431_10 * r_2 * h7_m2 + e_4 * f_1_39 * r_4 * h5_m2 + e_4 * fs_10_429_7 * r_6 * h3_m2;
    }

    // NOTE: the angular components are formed in 67 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph1_m1, ph3_m3, ph3_m1, ph5_m4, ph5_m3, ph5_m1, ph7_m4, ph7_m3, ph7_m1, ph9_m4, ph9_m3, ph9_m1, ab_2, pc_32, pc_33 : simd::cache_line_size())
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

        const auto h1_m1 = ph1_m1[k];
        const auto h3_m3 = ph3_m3[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h5_m4 = ph5_m4[k];
        const auto h5_m3 = ph5_m3[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h7_m4 = ph7_m4[k];
        const auto h7_m3 = ph7_m3[k];
        const auto h7_m1 = ph7_m1[k];
        const auto h9_m4 = ph9_m4[k];
        const auto h9_m3 = ph9_m3[k];
        const auto h9_m1 = ph9_m1[k];

        pc_32[k] = - e_0 * fs_105_32_42 * h1_m1 + e_1 * fs_15_8_105 * h3_m3 + e_1 * fs_6_1_7 * h3_m1 + e_1 * fs_21_4_42 * r_2 * h1_m1 - e_2 * fs_5_4_15 * h5_m3 - e_2 * fs_5_8_70 * h5_m1 - e_2 * fs_5_4_105 * r_2 * h3_m3 - e_2 * fs_4_1_7 * r_2 * h3_m1 - e_2 * fs_9_4_42 * r_4 * h1_m1 + e_3 * fs_215_286_2 * h7_m3 + e_3 * fs_205_858_6 * h7_m1 + e_3 * fs_5_13_15 * r_2 * h5_m3 + e_3 * fs_5_26_70 * r_2 * h5_m1 + e_3 * fs_5_22_105 * r_4 * h3_m3 + e_3 * fs_8_11_7 * r_4 * h3_m1 + e_3 * fs_1_3_42 * r_6 * h1_m1 - e_4 * fs_84_2431_55 * h9_m3 + e_4 * fs_28_2431_210 * h9_m1 - e_4 * fs_215_2431_2 * r_2 * h7_m3 - e_4 * fs_205_7293_6 * r_2 * h7_m1 - e_4 * fs_1_39_15 * r_4 * h5_m3 - e_4 * fs_1_78_70 * r_4 * h5_m1 - e_4 * fs_5_429_105 * r_6 * h3_m3 - e_4 * fs_16_429_7 * r_6 * h3_m1 - e_4 * fs_1_66_42 * r_8 * h1_m1;

        pc_33[k] = - e_2 * fs_5_4_15 * h5_m4 + e_3 * fs_50_143_11 * h7_m4 + e_3 * fs_5_13_15 * r_2 * h5_m4 - e_4 * fs_14_2431_2145 * h9_m4 - e_4 * fs_100_2431_11 * r_2 * h7_m4 - e_4 * fs_1_39_15 * r_4 * h5_m4;
    }

    // NOTE: the angular components are formed in 67 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph1_m1, ph3_m1, ph5_m5, ph5_m1, ph7_m5, ph7_m1, ph9_m5, ph9_m1, ab_2, pc_34 : simd::cache_line_size())
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

        const auto h1_m1 = ph1_m1[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h5_m5 = ph5_m5[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h7_m5 = ph7_m5[k];
        const auto h7_m1 = ph7_m1[k];
        const auto h9_m5 = ph9_m5[k];
        const auto h9_m1 = ph9_m1[k];

        pc_34[k] = - e_0 * fs_105_32_6 * h1_m1 - e_1 * f_147_8 * h3_m1 + e_1 * fs_21_4_6 * r_2 * h1_m1 + e_2 * fs_5_4_21 * h5_m5 + e_2 * fs_5_8_10 * h5_m1 + e_2 * f_49_4 * r_2 * h3_m1 - e_2 * fs_9_4_6 * r_4 * h1_m1 + e_3 * fs_5_286_154 * h7_m5 + e_3 * fs_175_858_42 * h7_m1 - e_3 * fs_5_13_21 * r_2 * h5_m5 - e_3 * fs_5_26_10 * r_2 * h5_m1 - e_3 * f_49_22 * r_4 * h3_m1 + e_3 * fs_1_3_6 * r_6 * h1_m1 - e_4 * fs_28_2431_429 * h9_m5 + e_4 * fs_28_2431_30 * h9_m1 - e_4 * fs_5_2431_154 * r_2 * h7_m5 - e_4 * fs_175_7293_42 * r_2 * h7_m1 + e_4 * fs_1_39_21 * r_4 * h5_m5 + e_4 * fs_1_78_10 * r_4 * h5_m1 + e_4 * f_49_429 * r_6 * h3_m1 - e_4 * fs_1_66_6 * r_8 * h1_m1;
    }

    // NOTE: the angular components are formed in 67 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, ph3_m2, ph3_p3, ph5_m2, ph5_p3, ph5_p5, ph7_m6, ph7_m2, ph7_p3, ph7_p5, ph9_m6, ph9_m2, ph9_p3, ph9_p5, ab_2, pc_35, pc_36 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h3_m2 = ph3_m2[k];
        const auto h3_p3 = ph3_p3[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h5_p5 = ph5_p5[k];
        const auto h7_m6 = ph7_m6[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h7_p3 = ph7_p3[k];
        const auto h7_p5 = ph7_p5[k];
        const auto h9_m6 = ph9_m6[k];
        const auto h9_m2 = ph9_m2[k];
        const auto h9_p3 = ph9_p3[k];
        const auto h9_p5 = ph9_p5[k];

        pc_35[k] = e_1 * fs_21_4_5 * h3_m2 + e_2 * fs_5_4_35 * h5_m2 - e_2 * fs_7_2_5 * r_2 * h3_m2 - e_3 * fs_5_143_2002 * h7_m6 + e_3 * fs_35_143_14 * h7_m2 - e_3 * fs_5_13_35 * r_2 * h5_m2 + e_3 * fs_7_11_5 * r_4 * h3_m2 - e_4 * fs_21_4862_1430 * h9_m6 + e_4 * fs_7_4862_330 * h9_m2 + e_4 * fs_10_2431_2002 * r_2 * h7_m6 - e_4 * fs_70_2431_14 * r_2 * h7_m2 + e_4 * fs_1_39_35 * r_4 * h5_m2 - e_4 * fs_14_429_5 * r_6 * h3_m2;

        pc_36[k] = e_1 * fs_3_8_210 * h3_p3 + e_2 * fs_5_4_30 * h5_p3 + e_2 * fs_5_4_6 * h5_p5 - e_2 * fs_1_4_210 * r_2 * h3_p3 + e_3 * f_175_143 * h7_p3 + e_3 * fs_70_143_11 * h7_p5 - e_3 * fs_5_13_30 * r_2 * h5_p3 - e_3 * fs_5_13_6 * r_2 * h5_p5 + e_3 * fs_1_22_210 * r_4 * h3_p3 + e_4 * fs_21_4862_110 * h9_p3 + e_4 * fs_7_4862_6006 * h9_p5 - e_4 * f_350_2431 * r_2 * h7_p3 - e_4 * fs_140_2431_11 * r_2 * h7_p5 + e_4 * fs_1_39_30 * r_4 * h5_p3 + e_4 * fs_1_39_6 * r_4 * h5_p5 - e_4 * fs_1_429_210 * r_6 * h3_p3;
    }

    // NOTE: the angular components are formed in 67 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, ph3_p2, ph5_p2, ph5_p4, ph7_p2, ph7_p4, ph9_p2, ph9_p4, ab_2, pc_37 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h3_p2 = ph3_p2[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h5_p4 = ph5_p4[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h7_p4 = ph7_p4[k];
        const auto h9_p2 = ph9_p2[k];
        const auto h9_p4 = ph9_p4[k];

        pc_37[k] = - e_1 * fs_33_16_70 * h3_p2 - e_2 * fs_5_8_10 * h5_p2 - e_2 * fs_5_4_30 * h5_p4 + e_2 * fs_11_8_70 * r_2 * h3_p2 + e_3 * f_35_26 * h7_p2 + e_3 * fs_35_286_22 * h7_p4 + e_3 * fs_5_26_10 * r_2 * h5_p2 + e_3 * fs_5_13_30 * r_2 * h5_p4 - e_3 * fs_1_4_70 * r_4 * h3_p2 + e_4 * fs_7_2431_1155 * h9_p2 + e_4 * fs_7_2431_4290 * h9_p4 - e_4 * f_35_221 * r_2 * h7_p2 - e_4 * fs_35_2431_22 * r_2 * h7_p4 - e_4 * fs_1_78_10 * r_4 * h5_p2 - e_4 * fs_1_39_30 * r_4 * h5_p4 + e_4 * fs_1_78_70 * r_6 * h3_p2;
    }

    // NOTE: the angular components are formed in 67 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph1_p1, ph3_p1, ph3_p3, ph5_p1, ph7_p1, ph7_p3, ph9_p1, ph9_p3, ab_2, pc_38 : simd::cache_line_size())
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

        const auto h1_p1 = ph1_p1[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h3_p3 = ph3_p3[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h7_p3 = ph7_p3[k];
        const auto h9_p1 = ph9_p1[k];
        const auto h9_p3 = ph9_p3[k];

        pc_38[k] = e_0 * fs_105_16_3 * h1_p1 + e_1 * fs_129_16_2 * h3_p1 + e_1 * fs_45_16_30 * h3_p3 - e_1 * fs_21_2_3 * r_2 * h1_p1 - e_2 * fs_5_2_5 * h5_p1 - e_2 * fs_43_8_2 * r_2 * h3_p1 - e_2 * fs_15_8_30 * r_2 * h3_p3 + e_2 * fs_9_2_3 * r_4 * h1_p1 + e_3 * fs_145_858_21 * h7_p1 - e_3 * fs_5_22_7 * h7_p3 + e_3 * fs_10_13_5 * r_2 * h5_p1 + e_3 * fs_43_44_2 * r_4 * h3_p1 + e_3 * fs_15_44_30 * r_4 * h3_p3 - e_3 * fs_2_3_3 * r_6 * h1_p1 + e_4 * fs_98_2431_15 * h9_p1 + e_4 * fs_21_2431_770 * h9_p3 - e_4 * fs_145_7293_21 * r_2 * h7_p1 + e_4 * fs_5_187_7 * r_2 * h7_p3 - e_4 * fs_2_39_5 * r_4 * h5_p1 - e_4 * fs_43_858_2 * r_6 * h3_p1 - e_4 * fs_5_286_30 * r_6 * h3_p3 + e_4 * fs_1_33_3 * r_8 * h1_p1;
    }

    // NOTE: the angular components are formed in 67 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph1_0, ph3_0, ph3_p2, ph5_0, ph5_p2, ph7_0, ph7_p2, ph9_0, ph9_p2, ab_2, pc_39 : simd::cache_line_size())
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

        const auto h1_0 = ph1_0[k];
        const auto h3_0 = ph3_0[k];
        const auto h3_p2 = ph3_p2[k];
        const auto h5_0 = ph5_0[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h7_0 = ph7_0[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h9_0 = ph9_0[k];
        const auto h9_p2 = ph9_p2[k];

        pc_39[k] = - e_0 * fs_105_8_6 * h1_0 + e_1 * fs_57_8_6 * h3_0 - e_1 * fs_75_16_10 * h3_p2 + e_1 * fs_21_1_6 * r_2 * h1_0 - e_2 * fs_5_4_6 * h5_0 + e_2 * fs_5_8_70 * h5_p2 - e_2 * fs_19_4_6 * r_2 * h3_0 + e_2 * fs_25_8_10 * r_2 * h3_p2 - e_2 * fs_9_1_6 * r_4 * h1_0 - e_3 * fs_35_429_6 * h7_0 - e_3 * fs_125_286_7 * h7_p2 + e_3 * fs_5_13_6 * r_2 * h5_0 - e_3 * fs_5_26_70 * r_2 * h5_p2 + e_3 * fs_19_22_6 * r_4 * h3_0 - e_3 * fs_25_44_10 * r_4 * h3_p2 + e_3 * fs_4_3_6 * r_6 * h1_0 + e_4 * fs_294_2431_6 * h9_0 + e_4 * fs_49_2431_165 * h9_p2 + e_4 * fs_70_7293_6 * r_2 * h7_0 + e_4 * fs_125_2431_7 * r_2 * h7_p2 - e_4 * fs_1_39_6 * r_4 * h5_0 + e_4 * fs_1_78_70 * r_4 * h5_p2 - e_4 * fs_19_429_6 * r_6 * h3_0 + e_4 * fs_25_858_10 * r_6 * h3_p2 - e_4 * fs_2_33_6 * r_8 * h1_0;
    }

    // NOTE: the angular components are formed in 67 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph1_m1, ph3_m2, ph3_m1, ph5_m2, ph5_m1, ph7_m2, ph7_m1, ph9_m2, ph9_m1, ab_2, pc_40, pc_41 : simd::cache_line_size())
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

        const auto h1_m1 = ph1_m1[k];
        const auto h3_m2 = ph3_m2[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h7_m1 = ph7_m1[k];
        const auto h9_m2 = ph9_m2[k];
        const auto h9_m1 = ph9_m1[k];

        pc_40[k] = - e_0 * fs_105_16_15 * h1_m1 + e_1 * fs_33_8_10 * h3_m1 + e_1 * fs_21_2_15 * r_2 * h1_m1 - e_2 * f_5_1 * h5_m1 - e_2 * fs_11_4_10 * r_2 * h3_m1 - e_2 * fs_9_2_15 * r_4 * h1_m1 + e_3 * fs_5_39_105 * h7_m1 + e_3 * f_20_13 * r_2 * h5_m1 + e_3 * fs_1_2_10 * r_4 * h3_m1 + e_3 * fs_2_3_15 * r_6 * h1_m1 - e_4 * fs_490_2431_3 * h9_m1 - e_4 * fs_10_663_105 * r_2 * h7_m1 - e_4 * f_4_39 * r_4 * h5_m1 - e_4 * fs_1_39_10 * r_6 * h3_m1 - e_4 * fs_1_33_15 * r_8 * h1_m1;

        pc_41[k] = e_1 * fs_75_16_10 * h3_m2 - e_2 * fs_5_8_70 * h5_m2 - e_2 * fs_25_8_10 * r_2 * h3_m2 + e_3 * fs_125_286_7 * h7_m2 + e_3 * fs_5_26_70 * r_2 * h5_m2 + e_3 * fs_25_44_10 * r_4 * h3_m2 - e_4 * fs_49_2431_165 * h9_m2 - e_4 * fs_125_2431_7 * r_2 * h7_m2 - e_4 * fs_1_78_70 * r_4 * h5_m2 - e_4 * fs_25_858_10 * r_6 * h3_m2;
    }

    // NOTE: the angular components are formed in 67 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph1_m1, ph3_m3, ph3_m1, ph5_m1, ph7_m3, ph7_m1, ph9_m3, ph9_m1, ab_2, pc_42 : simd::cache_line_size())
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

        const auto h1_m1 = ph1_m1[k];
        const auto h3_m3 = ph3_m3[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h7_m3 = ph7_m3[k];
        const auto h7_m1 = ph7_m1[k];
        const auto h9_m3 = ph9_m3[k];
        const auto h9_m1 = ph9_m1[k];

        pc_42[k] = - e_0 * fs_105_16_3 * h1_m1 - e_1 * fs_45_16_30 * h3_m3 - e_1 * fs_129_16_2 * h3_m1 + e_1 * fs_21_2_3 * r_2 * h1_m1 + e_2 * fs_5_2_5 * h5_m1 + e_2 * fs_15_8_30 * r_2 * h3_m3 + e_2 * fs_43_8_2 * r_2 * h3_m1 - e_2 * fs_9_2_3 * r_4 * h1_m1 + e_3 * fs_5_22_7 * h7_m3 - e_3 * fs_145_858_21 * h7_m1 - e_3 * fs_10_13_5 * r_2 * h5_m1 - e_3 * fs_15_44_30 * r_4 * h3_m3 - e_3 * fs_43_44_2 * r_4 * h3_m1 + e_3 * fs_2_3_3 * r_6 * h1_m1 - e_4 * fs_21_2431_770 * h9_m3 - e_4 * fs_98_2431_15 * h9_m1 - e_4 * fs_5_187_7 * r_2 * h7_m3 + e_4 * fs_145_7293_21 * r_2 * h7_m1 + e_4 * fs_2_39_5 * r_4 * h5_m1 + e_4 * fs_5_286_30 * r_6 * h3_m3 + e_4 * fs_43_858_2 * r_6 * h3_m1 - e_4 * fs_1_33_3 * r_8 * h1_m1;
    }

    // NOTE: the angular components are formed in 67 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, ph3_m2, ph5_m4, ph5_m2, ph7_m4, ph7_m2, ph9_m4, ph9_m2, ab_2, pc_43 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h3_m2 = ph3_m2[k];
        const auto h5_m4 = ph5_m4[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h7_m4 = ph7_m4[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h9_m4 = ph9_m4[k];
        const auto h9_m2 = ph9_m2[k];

        pc_43[k] = e_1 * fs_33_16_70 * h3_m2 + e_2 * fs_5_4_30 * h5_m4 + e_2 * fs_5_8_10 * h5_m2 - e_2 * fs_11_8_70 * r_2 * h3_m2 - e_3 * fs_35_286_22 * h7_m4 - e_3 * f_35_26 * h7_m2 - e_3 * fs_5_13_30 * r_2 * h5_m4 - e_3 * fs_5_26_10 * r_2 * h5_m2 + e_3 * fs_1_4_70 * r_4 * h3_m2 - e_4 * fs_7_2431_4290 * h9_m4 - e_4 * fs_7_2431_1155 * h9_m2 + e_4 * fs_35_2431_22 * r_2 * h7_m4 + e_4 * f_35_221 * r_2 * h7_m2 + e_4 * fs_1_39_30 * r_4 * h5_m4 + e_4 * fs_1_78_10 * r_4 * h5_m2 - e_4 * fs_1_78_70 * r_6 * h3_m2;
    }

    // NOTE: the angular components are formed in 67 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, ph3_m3, ph5_m5, ph5_m4, ph5_m3, ph7_m5, ph7_m4, ph7_m3, ph9_m5, ph9_m4, ph9_m3, ab_2, pc_44, pc_45, pc_46 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

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

        pc_44[k] = - e_1 * fs_3_8_210 * h3_m3 - e_2 * fs_5_4_6 * h5_m5 - e_2 * fs_5_4_30 * h5_m3 + e_2 * fs_1_4_210 * r_2 * h3_m3 - e_3 * fs_70_143_11 * h7_m5 - e_3 * f_175_143 * h7_m3 + e_3 * fs_5_13_6 * r_2 * h5_m5 + e_3 * fs_5_13_30 * r_2 * h5_m3 - e_3 * fs_1_22_210 * r_4 * h3_m3 - e_4 * fs_7_4862_6006 * h9_m5 - e_4 * fs_21_4862_110 * h9_m3 + e_4 * fs_140_2431_11 * r_2 * h7_m5 + e_4 * f_350_2431 * r_2 * h7_m3 - e_4 * fs_1_39_6 * r_4 * h5_m5 - e_4 * fs_1_39_30 * r_4 * h5_m3 + e_4 * fs_1_429_210 * r_6 * h3_m3;

        pc_45[k] = - e_2 * f_15_2 * h5_m4 - e_3 * fs_70_429_165 * h7_m4 + e_3 * f_30_13 * r_2 * h5_m4 - e_4 * fs_21_2431_143 * h9_m4 + e_4 * fs_140_7293_165 * r_2 * h7_m4 - e_4 * f_2_13 * r_4 * h5_m4;

        pc_46[k] = e_1 * fs_45_8_7 * h3_m3 + e_2 * f_15_2 * h5_m3 - e_2 * fs_15_4_7 * r_2 * h3_m3 - e_3 * fs_245_858_30 * h7_m3 - e_3 * f_30_13 * r_2 * h5_m3 + e_3 * fs_15_22_7 * r_4 * h3_m3 - e_4 * fs_84_2431_33 * h9_m3 + e_4 * fs_245_7293_30 * r_2 * h7_m3 + e_4 * f_2_13 * r_4 * h5_m3 - e_4 * fs_5_143_7 * r_6 * h3_m3;
    }

    // NOTE: the angular components are formed in 67 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph1_m1, ph3_m2, ph3_m1, ph5_m2, ph5_m1, ph7_m2, ph7_m1, ph9_m2, ph9_m1, ab_2, pc_47, pc_48 : simd::cache_line_size())
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

        const auto h1_m1 = ph1_m1[k];
        const auto h3_m2 = ph3_m2[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h7_m1 = ph7_m1[k];
        const auto h9_m2 = ph9_m2[k];
        const auto h9_m1 = ph9_m1[k];

        pc_47[k] = - e_1 * fs_15_1_3 * h3_m2 + e_2 * fs_5_4_21 * h5_m2 + e_2 * fs_10_1_3 * r_2 * h3_m2 - e_3 * fs_5_429_210 * h7_m2 - e_3 * fs_5_13_21 * r_2 * h5_m2 - e_3 * fs_20_11_3 * r_4 * h3_m2 - e_4 * fs_147_2431_22 * h9_m2 + e_4 * fs_10_7293_210 * r_2 * h7_m2 + e_4 * fs_1_39_21 * r_4 * h5_m2 + e_4 * fs_40_429_3 * r_6 * h3_m2;

        pc_48[k] = e_0 * fs_105_16_10 * h1_m1 + e_1 * fs_3_8_15 * h3_m1 - e_1 * fs_21_2_10 * r_2 * h1_m1 - e_2 * fs_5_4_6 * h5_m1 - e_2 * fs_1_4_15 * r_2 * h3_m1 + e_2 * fs_9_2_10 * r_4 * h1_m1 + e_3 * fs_115_858_70 * h7_m1 + e_3 * fs_5_13_6 * r_2 * h5_m1 + e_3 * fs_1_22_15 * r_4 * h3_m1 - e_3 * fs_2_3_10 * r_6 * h1_m1 - e_4 * fs_588_2431_2 * h9_m1 - e_4 * fs_115_7293_70 * r_2 * h7_m1 - e_4 * fs_1_39_6 * r_4 * h5_m1 - e_4 * fs_1_429_15 * r_6 * h3_m1 + e_4 * fs_1_33_10 * r_8 * h1_m1;
    }

    // NOTE: the angular components are formed in 67 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph1_0, ph1_p1, ph3_0, ph3_p1, ph5_0, ph5_p1, ph7_0, ph7_p1, ph9_0, ph9_p1, ab_2, pc_49, pc_50 : simd::cache_line_size())
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

        pc_49[k] = - e_0 * f_525_16 * h1_0 + e_1 * f_45_2 * h3_0 + e_1 * f_105_2 * r_2 * h1_0 - e_2 * f_15_2 * h5_0 - e_2 * f_15_1 * r_2 * h3_0 - e_2 * f_45_2 * r_4 * h1_0 + e_3 * f_700_429 * h7_0 + e_3 * f_30_13 * r_2 * h5_0 + e_3 * f_30_11 * r_4 * h3_0 + e_3 * f_10_3 * r_6 * h1_0 - e_4 * f_882_2431 * h9_0 - e_4 * f_1400_7293 * r_2 * h7_0 - e_4 * f_2_13 * r_4 * h5_0 - e_4 * f_20_143 * r_6 * h3_0 - e_4 * f_5_33 * r_8 * h1_0;

        pc_50[k] = e_0 * fs_105_16_10 * h1_p1 + e_1 * fs_3_8_15 * h3_p1 - e_1 * fs_21_2_10 * r_2 * h1_p1 - e_2 * fs_5_4_6 * h5_p1 - e_2 * fs_1_4_15 * r_2 * h3_p1 + e_2 * fs_9_2_10 * r_4 * h1_p1 + e_3 * fs_115_858_70 * h7_p1 + e_3 * fs_5_13_6 * r_2 * h5_p1 + e_3 * fs_1_22_15 * r_4 * h3_p1 - e_3 * fs_2_3_10 * r_6 * h1_p1 - e_4 * fs_588_2431_2 * h9_p1 - e_4 * fs_115_7293_70 * r_2 * h7_p1 - e_4 * fs_1_39_6 * r_4 * h5_p1 - e_4 * fs_1_429_15 * r_6 * h3_p1 + e_4 * fs_1_33_10 * r_8 * h1_p1;
    }

    // NOTE: the angular components are formed in 67 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, ph3_p2, ph3_p3, ph5_p2, ph5_p3, ph5_p4, ph7_p2, ph7_p3, ph7_p4, ph9_p2, ph9_p3, ph9_p4, ab_2, pc_51, pc_52, pc_53 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

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

        pc_51[k] = - e_1 * fs_15_1_3 * h3_p2 + e_2 * fs_5_4_21 * h5_p2 + e_2 * fs_10_1_3 * r_2 * h3_p2 - e_3 * fs_5_429_210 * h7_p2 - e_3 * fs_5_13_21 * r_2 * h5_p2 - e_3 * fs_20_11_3 * r_4 * h3_p2 - e_4 * fs_147_2431_22 * h9_p2 + e_4 * fs_10_7293_210 * r_2 * h7_p2 + e_4 * fs_1_39_21 * r_4 * h5_p2 + e_4 * fs_40_429_3 * r_6 * h3_p2;

        pc_52[k] = e_1 * fs_45_8_7 * h3_p3 + e_2 * f_15_2 * h5_p3 - e_2 * fs_15_4_7 * r_2 * h3_p3 - e_3 * fs_245_858_30 * h7_p3 - e_3 * f_30_13 * r_2 * h5_p3 + e_3 * fs_15_22_7 * r_4 * h3_p3 - e_4 * fs_84_2431_33 * h9_p3 + e_4 * fs_245_7293_30 * r_2 * h7_p3 + e_4 * f_2_13 * r_4 * h5_p3 - e_4 * fs_5_143_7 * r_6 * h3_p3;

        pc_53[k] = - e_2 * f_15_2 * h5_p4 - e_3 * fs_70_429_165 * h7_p4 + e_3 * f_30_13 * r_2 * h5_p4 - e_4 * fs_21_2431_143 * h9_p4 + e_4 * fs_140_7293_165 * r_2 * h7_p4 - e_4 * f_2_13 * r_4 * h5_p4;
    }

    // NOTE: the angular components are formed in 67 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, ph3_m3, ph5_m5, ph5_m3, ph7_m5, ph7_m3, ph9_m5, ph9_m3, ab_2, pc_54 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h3_m3 = ph3_m3[k];
        const auto h5_m5 = ph5_m5[k];
        const auto h5_m3 = ph5_m3[k];
        const auto h7_m5 = ph7_m5[k];
        const auto h7_m3 = ph7_m3[k];
        const auto h9_m5 = ph9_m5[k];
        const auto h9_m3 = ph9_m3[k];

        pc_54[k] = e_1 * fs_3_8_210 * h3_m3 - e_2 * fs_5_4_6 * h5_m5 + e_2 * fs_5_4_30 * h5_m3 - e_2 * fs_1_4_210 * r_2 * h3_m3 - e_3 * fs_70_143_11 * h7_m5 + e_3 * f_175_143 * h7_m3 + e_3 * fs_5_13_6 * r_2 * h5_m5 - e_3 * fs_5_13_30 * r_2 * h5_m3 + e_3 * fs_1_22_210 * r_4 * h3_m3 - e_4 * fs_7_4862_6006 * h9_m5 + e_4 * fs_21_4862_110 * h9_m3 + e_4 * fs_140_2431_11 * r_2 * h7_m5 - e_4 * f_350_2431 * r_2 * h7_m3 - e_4 * fs_1_39_6 * r_4 * h5_m5 + e_4 * fs_1_39_30 * r_4 * h5_m3 - e_4 * fs_1_429_210 * r_6 * h3_m3;
    }

    // NOTE: the angular components are formed in 67 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, ph3_m2, ph5_m4, ph5_m2, ph7_m4, ph7_m2, ph9_m4, ph9_m2, ab_2, pc_55 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h3_m2 = ph3_m2[k];
        const auto h5_m4 = ph5_m4[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h7_m4 = ph7_m4[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h9_m4 = ph9_m4[k];
        const auto h9_m2 = ph9_m2[k];

        pc_55[k] = - e_1 * fs_33_16_70 * h3_m2 + e_2 * fs_5_4_30 * h5_m4 - e_2 * fs_5_8_10 * h5_m2 + e_2 * fs_11_8_70 * r_2 * h3_m2 - e_3 * fs_35_286_22 * h7_m4 + e_3 * f_35_26 * h7_m2 - e_3 * fs_5_13_30 * r_2 * h5_m4 + e_3 * fs_5_26_10 * r_2 * h5_m2 - e_3 * fs_1_4_70 * r_4 * h3_m2 - e_4 * fs_7_2431_4290 * h9_m4 + e_4 * fs_7_2431_1155 * h9_m2 + e_4 * fs_35_2431_22 * r_2 * h7_m4 - e_4 * f_35_221 * r_2 * h7_m2 + e_4 * fs_1_39_30 * r_4 * h5_m4 - e_4 * fs_1_78_10 * r_4 * h5_m2 + e_4 * fs_1_78_70 * r_6 * h3_m2;
    }

    // NOTE: the angular components are formed in 67 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph1_m1, ph3_m3, ph3_m2, ph3_m1, ph5_m2, ph5_m1, ph7_m3, ph7_m2, ph7_m1, ph9_m3, ph9_m2, ph9_m1, ab_2, pc_56, pc_57 : simd::cache_line_size())
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

        const auto h1_m1 = ph1_m1[k];
        const auto h3_m3 = ph3_m3[k];
        const auto h3_m2 = ph3_m2[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h7_m3 = ph7_m3[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h7_m1 = ph7_m1[k];
        const auto h9_m3 = ph9_m3[k];
        const auto h9_m2 = ph9_m2[k];
        const auto h9_m1 = ph9_m1[k];

        pc_56[k] = e_0 * fs_105_16_3 * h1_m1 - e_1 * fs_45_16_30 * h3_m3 + e_1 * fs_129_16_2 * h3_m1 - e_1 * fs_21_2_3 * r_2 * h1_m1 - e_2 * fs_5_2_5 * h5_m1 + e_2 * fs_15_8_30 * r_2 * h3_m3 - e_2 * fs_43_8_2 * r_2 * h3_m1 + e_2 * fs_9_2_3 * r_4 * h1_m1 + e_3 * fs_5_22_7 * h7_m3 + e_3 * fs_145_858_21 * h7_m1 + e_3 * fs_10_13_5 * r_2 * h5_m1 - e_3 * fs_15_44_30 * r_4 * h3_m3 + e_3 * fs_43_44_2 * r_4 * h3_m1 - e_3 * fs_2_3_3 * r_6 * h1_m1 - e_4 * fs_21_2431_770 * h9_m3 + e_4 * fs_98_2431_15 * h9_m1 - e_4 * fs_5_187_7 * r_2 * h7_m3 - e_4 * fs_145_7293_21 * r_2 * h7_m1 - e_4 * fs_2_39_5 * r_4 * h5_m1 + e_4 * fs_5_286_30 * r_6 * h3_m3 - e_4 * fs_43_858_2 * r_6 * h3_m1 + e_4 * fs_1_33_3 * r_8 * h1_m1;

        pc_57[k] = e_1 * fs_75_16_10 * h3_m2 - e_2 * fs_5_8_70 * h5_m2 - e_2 * fs_25_8_10 * r_2 * h3_m2 + e_3 * fs_125_286_7 * h7_m2 + e_3 * fs_5_26_70 * r_2 * h5_m2 + e_3 * fs_25_44_10 * r_4 * h3_m2 - e_4 * fs_49_2431_165 * h9_m2 - e_4 * fs_125_2431_7 * r_2 * h7_m2 - e_4 * fs_1_78_70 * r_4 * h5_m2 - e_4 * fs_25_858_10 * r_6 * h3_m2;
    }

    // NOTE: the angular components are formed in 67 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph1_p1, ph3_p1, ph5_p1, ph7_p1, ph9_p1, ab_2, pc_58 : simd::cache_line_size())
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

        const auto h1_p1 = ph1_p1[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h9_p1 = ph9_p1[k];

        pc_58[k] = - e_0 * fs_105_16_15 * h1_p1 + e_1 * fs_33_8_10 * h3_p1 + e_1 * fs_21_2_15 * r_2 * h1_p1 - e_2 * f_5_1 * h5_p1 - e_2 * fs_11_4_10 * r_2 * h3_p1 - e_2 * fs_9_2_15 * r_4 * h1_p1 + e_3 * fs_5_39_105 * h7_p1 + e_3 * f_20_13 * r_2 * h5_p1 + e_3 * fs_1_2_10 * r_4 * h3_p1 + e_3 * fs_2_3_15 * r_6 * h1_p1 - e_4 * fs_490_2431_3 * h9_p1 - e_4 * fs_10_663_105 * r_2 * h7_p1 - e_4 * f_4_39 * r_4 * h5_p1 - e_4 * fs_1_39_10 * r_6 * h3_p1 - e_4 * fs_1_33_15 * r_8 * h1_p1;
    }

    // NOTE: the angular components are formed in 67 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph1_0, ph3_0, ph3_p2, ph5_0, ph5_p2, ph7_0, ph7_p2, ph9_0, ph9_p2, ab_2, pc_59 : simd::cache_line_size())
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

        const auto h1_0 = ph1_0[k];
        const auto h3_0 = ph3_0[k];
        const auto h3_p2 = ph3_p2[k];
        const auto h5_0 = ph5_0[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h7_0 = ph7_0[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h9_0 = ph9_0[k];
        const auto h9_p2 = ph9_p2[k];

        pc_59[k] = - e_0 * fs_105_8_6 * h1_0 + e_1 * fs_57_8_6 * h3_0 + e_1 * fs_75_16_10 * h3_p2 + e_1 * fs_21_1_6 * r_2 * h1_0 - e_2 * fs_5_4_6 * h5_0 - e_2 * fs_5_8_70 * h5_p2 - e_2 * fs_19_4_6 * r_2 * h3_0 - e_2 * fs_25_8_10 * r_2 * h3_p2 - e_2 * fs_9_1_6 * r_4 * h1_0 - e_3 * fs_35_429_6 * h7_0 + e_3 * fs_125_286_7 * h7_p2 + e_3 * fs_5_13_6 * r_2 * h5_0 + e_3 * fs_5_26_70 * r_2 * h5_p2 + e_3 * fs_19_22_6 * r_4 * h3_0 + e_3 * fs_25_44_10 * r_4 * h3_p2 + e_3 * fs_4_3_6 * r_6 * h1_0 + e_4 * fs_294_2431_6 * h9_0 - e_4 * fs_49_2431_165 * h9_p2 + e_4 * fs_70_7293_6 * r_2 * h7_0 - e_4 * fs_125_2431_7 * r_2 * h7_p2 - e_4 * fs_1_39_6 * r_4 * h5_0 - e_4 * fs_1_78_70 * r_4 * h5_p2 - e_4 * fs_19_429_6 * r_6 * h3_0 - e_4 * fs_25_858_10 * r_6 * h3_p2 - e_4 * fs_2_33_6 * r_8 * h1_0;
    }

    // NOTE: the angular components are formed in 67 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph1_p1, ph3_p1, ph3_p3, ph5_p1, ph7_p1, ph7_p3, ph9_p1, ph9_p3, ab_2, pc_60 : simd::cache_line_size())
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

        const auto h1_p1 = ph1_p1[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h3_p3 = ph3_p3[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h7_p3 = ph7_p3[k];
        const auto h9_p1 = ph9_p1[k];
        const auto h9_p3 = ph9_p3[k];

        pc_60[k] = e_0 * fs_105_16_3 * h1_p1 + e_1 * fs_129_16_2 * h3_p1 - e_1 * fs_45_16_30 * h3_p3 - e_1 * fs_21_2_3 * r_2 * h1_p1 - e_2 * fs_5_2_5 * h5_p1 - e_2 * fs_43_8_2 * r_2 * h3_p1 + e_2 * fs_15_8_30 * r_2 * h3_p3 + e_2 * fs_9_2_3 * r_4 * h1_p1 + e_3 * fs_145_858_21 * h7_p1 + e_3 * fs_5_22_7 * h7_p3 + e_3 * fs_10_13_5 * r_2 * h5_p1 + e_3 * fs_43_44_2 * r_4 * h3_p1 - e_3 * fs_15_44_30 * r_4 * h3_p3 - e_3 * fs_2_3_3 * r_6 * h1_p1 + e_4 * fs_98_2431_15 * h9_p1 - e_4 * fs_21_2431_770 * h9_p3 - e_4 * fs_145_7293_21 * r_2 * h7_p1 - e_4 * fs_5_187_7 * r_2 * h7_p3 - e_4 * fs_2_39_5 * r_4 * h5_p1 - e_4 * fs_43_858_2 * r_6 * h3_p1 + e_4 * fs_5_286_30 * r_6 * h3_p3 + e_4 * fs_1_33_3 * r_8 * h1_p1;
    }

    // NOTE: the angular components are formed in 67 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, ph3_p2, ph5_p2, ph5_p4, ph7_p2, ph7_p4, ph9_p2, ph9_p4, ab_2, pc_61 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h3_p2 = ph3_p2[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h5_p4 = ph5_p4[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h7_p4 = ph7_p4[k];
        const auto h9_p2 = ph9_p2[k];
        const auto h9_p4 = ph9_p4[k];

        pc_61[k] = - e_1 * fs_33_16_70 * h3_p2 - e_2 * fs_5_8_10 * h5_p2 + e_2 * fs_5_4_30 * h5_p4 + e_2 * fs_11_8_70 * r_2 * h3_p2 + e_3 * f_35_26 * h7_p2 - e_3 * fs_35_286_22 * h7_p4 + e_3 * fs_5_26_10 * r_2 * h5_p2 - e_3 * fs_5_13_30 * r_2 * h5_p4 - e_3 * fs_1_4_70 * r_4 * h3_p2 + e_4 * fs_7_2431_1155 * h9_p2 - e_4 * fs_7_2431_4290 * h9_p4 - e_4 * f_35_221 * r_2 * h7_p2 + e_4 * fs_35_2431_22 * r_2 * h7_p4 - e_4 * fs_1_78_10 * r_4 * h5_p2 + e_4 * fs_1_39_30 * r_4 * h5_p4 + e_4 * fs_1_78_70 * r_6 * h3_p2;
    }

    // NOTE: the angular components are formed in 67 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, ph3_m2, ph3_p3, ph5_m2, ph5_p3, ph5_p5, ph7_m6, ph7_m2, ph7_p3, ph7_p5, ph9_m6, ph9_m2, ph9_p3, ph9_p5, ab_2, pc_62, pc_63 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h3_m2 = ph3_m2[k];
        const auto h3_p3 = ph3_p3[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h5_p5 = ph5_p5[k];
        const auto h7_m6 = ph7_m6[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h7_p3 = ph7_p3[k];
        const auto h7_p5 = ph7_p5[k];
        const auto h9_m6 = ph9_m6[k];
        const auto h9_m2 = ph9_m2[k];
        const auto h9_p3 = ph9_p3[k];
        const auto h9_p5 = ph9_p5[k];

        pc_62[k] = e_1 * fs_3_8_210 * h3_p3 + e_2 * fs_5_4_30 * h5_p3 - e_2 * fs_5_4_6 * h5_p5 - e_2 * fs_1_4_210 * r_2 * h3_p3 + e_3 * f_175_143 * h7_p3 - e_3 * fs_70_143_11 * h7_p5 - e_3 * fs_5_13_30 * r_2 * h5_p3 + e_3 * fs_5_13_6 * r_2 * h5_p5 + e_3 * fs_1_22_210 * r_4 * h3_p3 + e_4 * fs_21_4862_110 * h9_p3 - e_4 * fs_7_4862_6006 * h9_p5 - e_4 * f_350_2431 * r_2 * h7_p3 + e_4 * fs_140_2431_11 * r_2 * h7_p5 + e_4 * fs_1_39_30 * r_4 * h5_p3 - e_4 * fs_1_39_6 * r_4 * h5_p5 - e_4 * fs_1_429_210 * r_6 * h3_p3;

        pc_63[k] = - e_1 * fs_21_4_5 * h3_m2 - e_2 * fs_5_4_35 * h5_m2 + e_2 * fs_7_2_5 * r_2 * h3_m2 - e_3 * fs_5_143_2002 * h7_m6 - e_3 * fs_35_143_14 * h7_m2 + e_3 * fs_5_13_35 * r_2 * h5_m2 - e_3 * fs_7_11_5 * r_4 * h3_m2 - e_4 * fs_21_4862_1430 * h9_m6 - e_4 * fs_7_4862_330 * h9_m2 + e_4 * fs_10_2431_2002 * r_2 * h7_m6 + e_4 * fs_70_2431_14 * r_2 * h7_m2 - e_4 * fs_1_39_35 * r_4 * h5_m2 + e_4 * fs_14_429_5 * r_6 * h3_m2;
    }

    // NOTE: the angular components are formed in 67 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph1_m1, ph3_m1, ph5_m5, ph5_m4, ph5_m1, ph7_m5, ph7_m4, ph7_m1, ph9_m5, ph9_m4, ph9_m1, ab_2, pc_64, pc_65 : simd::cache_line_size())
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

        const auto h1_m1 = ph1_m1[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h5_m5 = ph5_m5[k];
        const auto h5_m4 = ph5_m4[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h7_m5 = ph7_m5[k];
        const auto h7_m4 = ph7_m4[k];
        const auto h7_m1 = ph7_m1[k];
        const auto h9_m5 = ph9_m5[k];
        const auto h9_m4 = ph9_m4[k];
        const auto h9_m1 = ph9_m1[k];

        pc_64[k] = e_0 * fs_105_32_6 * h1_m1 + e_1 * f_147_8 * h3_m1 - e_1 * fs_21_4_6 * r_2 * h1_m1 + e_2 * fs_5_4_21 * h5_m5 - e_2 * fs_5_8_10 * h5_m1 - e_2 * f_49_4 * r_2 * h3_m1 + e_2 * fs_9_4_6 * r_4 * h1_m1 + e_3 * fs_5_286_154 * h7_m5 - e_3 * fs_175_858_42 * h7_m1 - e_3 * fs_5_13_21 * r_2 * h5_m5 + e_3 * fs_5_26_10 * r_2 * h5_m1 + e_3 * f_49_22 * r_4 * h3_m1 - e_3 * fs_1_3_6 * r_6 * h1_m1 - e_4 * fs_28_2431_429 * h9_m5 - e_4 * fs_28_2431_30 * h9_m1 - e_4 * fs_5_2431_154 * r_2 * h7_m5 + e_4 * fs_175_7293_42 * r_2 * h7_m1 + e_4 * fs_1_39_21 * r_4 * h5_m5 - e_4 * fs_1_78_10 * r_4 * h5_m1 - e_4 * f_49_429 * r_6 * h3_m1 + e_4 * fs_1_66_6 * r_8 * h1_m1;

        pc_65[k] = - e_2 * fs_5_4_15 * h5_m4 + e_3 * fs_50_143_11 * h7_m4 + e_3 * fs_5_13_15 * r_2 * h5_m4 - e_4 * fs_14_2431_2145 * h9_m4 - e_4 * fs_100_2431_11 * r_2 * h7_m4 - e_4 * fs_1_39_15 * r_4 * h5_m4;
    }

    // NOTE: the angular components are formed in 67 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph1_m1, ph3_m3, ph3_m1, ph5_m3, ph5_m1, ph7_m3, ph7_m1, ph9_m3, ph9_m1, ab_2, pc_66 : simd::cache_line_size())
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

        const auto h1_m1 = ph1_m1[k];
        const auto h3_m3 = ph3_m3[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h5_m3 = ph5_m3[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h7_m3 = ph7_m3[k];
        const auto h7_m1 = ph7_m1[k];
        const auto h9_m3 = ph9_m3[k];
        const auto h9_m1 = ph9_m1[k];

        pc_66[k] = e_0 * fs_105_32_42 * h1_m1 + e_1 * fs_15_8_105 * h3_m3 - e_1 * fs_6_1_7 * h3_m1 - e_1 * fs_21_4_42 * r_2 * h1_m1 - e_2 * fs_5_4_15 * h5_m3 + e_2 * fs_5_8_70 * h5_m1 - e_2 * fs_5_4_105 * r_2 * h3_m3 + e_2 * fs_4_1_7 * r_2 * h3_m1 + e_2 * fs_9_4_42 * r_4 * h1_m1 + e_3 * fs_215_286_2 * h7_m3 - e_3 * fs_205_858_6 * h7_m1 + e_3 * fs_5_13_15 * r_2 * h5_m3 - e_3 * fs_5_26_70 * r_2 * h5_m1 + e_3 * fs_5_22_105 * r_4 * h3_m3 - e_3 * fs_8_11_7 * r_4 * h3_m1 - e_3 * fs_1_3_42 * r_6 * h1_m1 - e_4 * fs_84_2431_55 * h9_m3 - e_4 * fs_28_2431_210 * h9_m1 - e_4 * fs_215_2431_2 * r_2 * h7_m3 + e_4 * fs_205_7293_6 * r_2 * h7_m1 - e_4 * fs_1_39_15 * r_4 * h5_m3 + e_4 * fs_1_78_70 * r_4 * h5_m1 - e_4 * fs_5_429_105 * r_6 * h3_m3 + e_4 * fs_16_429_7 * r_6 * h3_m1 + e_4 * fs_1_66_42 * r_8 * h1_m1;
    }

    // NOTE: the angular components are formed in 67 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, ph3_p2, ph5_p2, ph7_p2, ph9_p2, ab_2, pc_67 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h3_p2 = ph3_p2[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h9_p2 = ph9_p2[k];

        pc_67[k] = - e_1 * fs_15_4_7 * h3_p2 + e_2 * f_5_4 * h5_p2 + e_2 * fs_5_2_7 * r_2 * h3_p2 + e_3 * fs_20_143_10 * h7_p2 - e_3 * f_5_13 * r_2 * h5_p2 - e_3 * fs_5_11_7 * r_4 * h3_p2 - e_4 * fs_35_2431_462 * h9_p2 - e_4 * fs_40_2431_10 * r_2 * h7_p2 + e_4 * f_1_39 * r_4 * h5_p2 + e_4 * fs_10_429_7 * r_6 * h3_p2;
    }

    // NOTE: the angular components are formed in 67 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph1_p1, ph3_p1, ph3_p3, ph5_p1, ph5_p3, ph7_p1, ph7_p3, ph9_p1, ph9_p3, ab_2, pc_68 : simd::cache_line_size())
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

        const auto h1_p1 = ph1_p1[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h3_p3 = ph3_p3[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h7_p3 = ph7_p3[k];
        const auto h9_p1 = ph9_p1[k];
        const auto h9_p3 = ph9_p3[k];

        pc_68[k] = - e_0 * fs_105_32_42 * h1_p1 + e_1 * fs_6_1_7 * h3_p1 + e_1 * fs_15_8_105 * h3_p3 + e_1 * fs_21_4_42 * r_2 * h1_p1 - e_2 * fs_5_8_70 * h5_p1 - e_2 * fs_5_4_15 * h5_p3 - e_2 * fs_4_1_7 * r_2 * h3_p1 - e_2 * fs_5_4_105 * r_2 * h3_p3 - e_2 * fs_9_4_42 * r_4 * h1_p1 + e_3 * fs_205_858_6 * h7_p1 + e_3 * fs_215_286_2 * h7_p3 + e_3 * fs_5_26_70 * r_2 * h5_p1 + e_3 * fs_5_13_15 * r_2 * h5_p3 + e_3 * fs_8_11_7 * r_4 * h3_p1 + e_3 * fs_5_22_105 * r_4 * h3_p3 + e_3 * fs_1_3_42 * r_6 * h1_p1 + e_4 * fs_28_2431_210 * h9_p1 - e_4 * fs_84_2431_55 * h9_p3 - e_4 * fs_205_7293_6 * r_2 * h7_p1 - e_4 * fs_215_2431_2 * r_2 * h7_p3 - e_4 * fs_1_78_70 * r_4 * h5_p1 - e_4 * fs_1_39_15 * r_4 * h5_p3 - e_4 * fs_16_429_7 * r_6 * h3_p1 - e_4 * fs_5_429_105 * r_6 * h3_p3 - e_4 * fs_1_66_42 * r_8 * h1_p1;
    }

    // NOTE: the angular components are formed in 67 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph1_0, ph3_0, ph5_0, ph5_p4, ph7_0, ph7_p4, ph9_0, ph9_p4, ab_2, pc_69 : simd::cache_line_size())
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

        const auto h1_0 = ph1_0[k];
        const auto h3_0 = ph3_0[k];
        const auto h5_0 = ph5_0[k];
        const auto h5_p4 = ph5_p4[k];
        const auto h7_0 = ph7_0[k];
        const auto h7_p4 = ph7_p4[k];
        const auto h9_0 = ph9_0[k];
        const auto h9_p4 = ph9_p4[k];

        pc_69[k] = - e_0 * fs_105_16_21 * h1_0 + e_1 * fs_3_4_21 * h3_0 + e_1 * fs_21_2_21 * r_2 * h1_0 + e_2 * fs_5_4_21 * h5_0 - e_2 * fs_5_4_15 * h5_p4 - e_2 * fs_1_2_21 * r_2 * h3_0 - e_2 * fs_9_2_21 * r_4 * h1_0 - e_3 * fs_160_429_21 * h7_0 + e_3 * fs_50_143_11 * h7_p4 - e_3 * fs_5_13_21 * r_2 * h5_0 + e_3 * fs_5_13_15 * r_2 * h5_p4 + e_3 * fs_1_11_21 * r_4 * h3_0 + e_3 * fs_2_3_21 * r_6 * h1_0 - e_4 * fs_84_2431_21 * h9_0 - e_4 * fs_14_2431_2145 * h9_p4 + e_4 * fs_320_7293_21 * r_2 * h7_0 - e_4 * fs_100_2431_11 * r_2 * h7_p4 + e_4 * fs_1_39_21 * r_4 * h5_0 - e_4 * fs_1_39_15 * r_4 * h5_p4 - e_4 * fs_2_429_21 * r_6 * h3_0 - e_4 * fs_1_33_21 * r_8 * h1_0;
    }

    // NOTE: the angular components are formed in 67 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph1_p1, ph3_p1, ph5_p1, ph5_p5, ph7_p1, ph7_p5, ph9_p1, ph9_p5, ab_2, pc_70 : simd::cache_line_size())
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

        const auto h1_p1 = ph1_p1[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h5_p5 = ph5_p5[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h7_p5 = ph7_p5[k];
        const auto h9_p1 = ph9_p1[k];
        const auto h9_p5 = ph9_p5[k];

        pc_70[k] = e_0 * fs_105_32_6 * h1_p1 + e_1 * f_147_8 * h3_p1 - e_1 * fs_21_4_6 * r_2 * h1_p1 - e_2 * fs_5_8_10 * h5_p1 + e_2 * fs_5_4_21 * h5_p5 - e_2 * f_49_4 * r_2 * h3_p1 + e_2 * fs_9_4_6 * r_4 * h1_p1 - e_3 * fs_175_858_42 * h7_p1 + e_3 * fs_5_286_154 * h7_p5 + e_3 * fs_5_26_10 * r_2 * h5_p1 - e_3 * fs_5_13_21 * r_2 * h5_p5 + e_3 * f_49_22 * r_4 * h3_p1 - e_3 * fs_1_3_6 * r_6 * h1_p1 - e_4 * fs_28_2431_30 * h9_p1 - e_4 * fs_28_2431_429 * h9_p5 + e_4 * fs_175_7293_42 * r_2 * h7_p1 - e_4 * fs_5_2431_154 * r_2 * h7_p5 - e_4 * fs_1_78_10 * r_4 * h5_p1 + e_4 * fs_1_39_21 * r_4 * h5_p5 - e_4 * f_49_429 * r_6 * h3_p1 + e_4 * fs_1_66_6 * r_8 * h1_p1;
    }

    // NOTE: the angular components are formed in 67 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, ph3_p2, ph5_p2, ph7_p2, ph7_p6, ph9_p2, ph9_p6, ab_2, pc_71 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h3_p2 = ph3_p2[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h7_p6 = ph7_p6[k];
        const auto h9_p2 = ph9_p2[k];
        const auto h9_p6 = ph9_p6[k];

        pc_71[k] = - e_1 * fs_21_4_5 * h3_p2 - e_2 * fs_5_4_35 * h5_p2 + e_2 * fs_7_2_5 * r_2 * h3_p2 - e_3 * fs_35_143_14 * h7_p2 - e_3 * fs_5_143_2002 * h7_p6 + e_3 * fs_5_13_35 * r_2 * h5_p2 - e_3 * fs_7_11_5 * r_4 * h3_p2 - e_4 * fs_7_4862_330 * h9_p2 - e_4 * fs_21_4862_1430 * h9_p6 + e_4 * fs_70_2431_14 * r_2 * h7_p2 + e_4 * fs_10_2431_2002 * r_2 * h7_p6 - e_4 * fs_1_39_35 * r_4 * h5_p2 + e_4 * fs_14_429_5 * r_6 * h3_p2;
    }

    // NOTE: the angular components are formed in 67 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph1_m1, ph3_m1, ph5_m1, ph7_m7, ph7_m6, ph7_m1, ph9_m7, ph9_m6, ph9_m1, ab_2, pc_72, pc_73 : simd::cache_line_size())
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

        const auto h1_m1 = ph1_m1[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h7_m7 = ph7_m7[k];
        const auto h7_m6 = ph7_m6[k];
        const auto h7_m1 = ph7_m1[k];
        const auto h9_m7 = ph9_m7[k];
        const auto h9_m6 = ph9_m6[k];
        const auto h9_m1 = ph9_m1[k];

        pc_72[k] = e_0 * fs_105_32_2 * h1_m1 + e_1 * fs_21_2_3 * h3_m1 - e_1 * fs_21_4_2 * r_2 * h1_m1 + e_2 * fs_5_4_30 * h5_m1 - e_2 * fs_7_1_3 * r_2 * h3_m1 + e_2 * fs_9_4_2 * r_4 * h1_m1 - e_3 * fs_35_858_858 * h7_m7 + e_3 * fs_70_429_14 * h7_m1 - e_3 * fs_5_13_30 * r_2 * h5_m1 + e_3 * fs_14_11_3 * r_4 * h3_m1 - e_3 * fs_1_3_2 * r_6 * h1_m1 - e_4 * fs_21_2431_715 * h9_m7 + e_4 * fs_21_4862_10 * h9_m1 + e_4 * fs_35_7293_858 * r_2 * h7_m7 - e_4 * fs_140_7293_14 * r_2 * h7_m1 + e_4 * fs_1_39_30 * r_4 * h5_m1 - e_4 * fs_28_429_3 * r_6 * h3_m1 + e_4 * fs_1_66_2 * r_8 * h1_m1;

        pc_73[k] = e_3 * fs_25_1716_6006 * h7_m6 - e_4 * fs_21_4862_4290 * h9_m6 - e_4 * fs_25_14586_6006 * r_2 * h7_m6;
    }

    // NOTE: the angular components are formed in 67 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph1_m1, ph3_m1, ph5_m5, ph7_m5, ph7_m1, ph9_m5, ph9_m1, ab_2, pc_74 : simd::cache_line_size())
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

        const auto h1_m1 = ph1_m1[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h5_m5 = ph5_m5[k];
        const auto h7_m5 = ph7_m5[k];
        const auto h7_m1 = ph7_m1[k];
        const auto h9_m5 = ph9_m5[k];
        const auto h9_m1 = ph9_m1[k];

        pc_74[k] = e_0 * fs_105_16_14 * h1_m1 - e_1 * fs_27_8_21 * h3_m1 - e_1 * fs_21_2_14 * r_2 * h1_m1 - e_2 * f_15_2 * h5_m5 + e_2 * fs_9_4_21 * r_2 * h3_m1 + e_2 * fs_9_2_14 * r_4 * h1_m1 + e_3 * fs_265_1716_66 * h7_m5 + e_3 * fs_125_132_2 * h7_m1 + e_3 * f_30_13 * r_2 * h5_m5 - e_3 * fs_9_22_21 * r_4 * h3_m1 - e_3 * fs_2_3_14 * r_6 * h1_m1 - e_4 * fs_21_2431_1001 * h9_m5 + e_4 * fs_21_2431_70 * h9_m1 - e_4 * fs_265_14586_66 * r_2 * h7_m5 - e_4 * fs_125_1122_2 * r_2 * h7_m1 - e_4 * f_2_13 * r_4 * h5_m5 + e_4 * fs_3_143_21 * r_6 * h3_m1 + e_4 * fs_1_33_14 * r_8 * h1_m1;
    }

    // NOTE: the angular components are formed in 67 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, ph3_m2, ph3_p3, ph5_m2, ph5_p3, ph7_m4, ph7_m2, ph7_p3, ph9_m4, ph9_m2, ph9_p3, ab_2, pc_75, pc_76 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h3_m2 = ph3_m2[k];
        const auto h3_p3 = ph3_p3[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h7_m4 = ph7_m4[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h7_p3 = ph7_p3[k];
        const auto h9_m4 = ph9_m4[k];
        const auto h9_m2 = ph9_m2[k];
        const auto h9_p3 = ph9_p3[k];

        pc_75[k] = - e_1 * fs_3_8_105 * h3_m2 + e_2 * fs_5_4_15 * h5_m2 + e_2 * fs_1_4_105 * r_2 * h3_m2 + e_3 * fs_5_66_33 * h7_m4 - e_3 * fs_815_1716_6 * h7_m2 - e_3 * fs_5_13_15 * r_2 * h5_m2 - e_3 * fs_1_22_105 * r_4 * h3_m2 - e_4 * fs_21_2431_715 * h9_m4 - e_4 * fs_21_4862_770 * h9_m2 - e_4 * fs_5_561_33 * r_2 * h7_m4 + e_4 * fs_815_14586_6 * r_2 * h7_m2 + e_4 * fs_1_39_15 * r_4 * h5_m2 + e_4 * fs_1_429_105 * r_6 * h3_m2;

        pc_76[k] = - e_1 * fs_45_4_7 * h3_p3 + e_2 * f_15_2 * h5_p3 + e_2 * fs_15_2_7 * r_2 * h3_p3 - e_3 * fs_115_858_30 * h7_p3 - e_3 * f_30_13 * r_2 * h5_p3 - e_3 * fs_15_11_7 * r_4 * h3_p3 - e_4 * fs_105_2431_33 * h9_p3 + e_4 * fs_115_7293_30 * r_2 * h7_p3 + e_4 * f_2_13 * r_4 * h5_p3 + e_4 * fs_10_143_7 * r_6 * h3_p3;
    }

    // NOTE: the angular components are formed in 67 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, ph3_p2, ph5_p2, ph7_p2, ph7_p4, ph9_p2, ph9_p4, ab_2, pc_77 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h3_p2 = ph3_p2[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h7_p4 = ph7_p4[k];
        const auto h9_p2 = ph9_p2[k];
        const auto h9_p4 = ph9_p4[k];

        pc_77[k] = e_1 * fs_3_8_105 * h3_p2 - e_2 * fs_5_4_15 * h5_p2 - e_2 * fs_1_4_105 * r_2 * h3_p2 + e_3 * fs_815_1716_6 * h7_p2 + e_3 * fs_5_66_33 * h7_p4 + e_3 * fs_5_13_15 * r_2 * h5_p2 + e_3 * fs_1_22_105 * r_4 * h3_p2 + e_4 * fs_21_4862_770 * h9_p2 - e_4 * fs_21_2431_715 * h9_p4 - e_4 * fs_815_14586_6 * r_2 * h7_p2 - e_4 * fs_5_561_33 * r_2 * h7_p4 - e_4 * fs_1_39_15 * r_4 * h5_p2 - e_4 * fs_1_429_105 * r_6 * h3_p2;
    }

    // NOTE: the angular components are formed in 67 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph1_p1, ph3_p1, ph5_p5, ph7_p1, ph7_p5, ph9_p1, ph9_p5, ab_2, pc_78 : simd::cache_line_size())
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

        const auto h1_p1 = ph1_p1[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h5_p5 = ph5_p5[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h7_p5 = ph7_p5[k];
        const auto h9_p1 = ph9_p1[k];
        const auto h9_p5 = ph9_p5[k];

        pc_78[k] = - e_0 * fs_105_16_14 * h1_p1 + e_1 * fs_27_8_21 * h3_p1 + e_1 * fs_21_2_14 * r_2 * h1_p1 - e_2 * f_15_2 * h5_p5 - e_2 * fs_9_4_21 * r_2 * h3_p1 - e_2 * fs_9_2_14 * r_4 * h1_p1 - e_3 * fs_125_132_2 * h7_p1 + e_3 * fs_265_1716_66 * h7_p5 + e_3 * f_30_13 * r_2 * h5_p5 + e_3 * fs_9_22_21 * r_4 * h3_p1 + e_3 * fs_2_3_14 * r_6 * h1_p1 - e_4 * fs_21_2431_70 * h9_p1 - e_4 * fs_21_2431_1001 * h9_p5 + e_4 * fs_125_1122_2 * r_2 * h7_p1 - e_4 * fs_265_14586_66 * r_2 * h7_p5 - e_4 * f_2_13 * r_4 * h5_p5 - e_4 * fs_3_143_21 * r_6 * h3_p1 - e_4 * fs_1_33_14 * r_8 * h1_p1;
    }

    // NOTE: the angular components are formed in 67 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph1_0, ph3_0, ph5_0, ph7_0, ph7_p6, ph9_0, ph9_p6, ab_2, pc_79 : simd::cache_line_size())
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

        const auto h1_0 = ph1_0[k];
        const auto h3_0 = ph3_0[k];
        const auto h5_0 = ph5_0[k];
        const auto h7_0 = ph7_0[k];
        const auto h7_p6 = ph7_p6[k];
        const auto h9_0 = ph9_0[k];
        const auto h9_p6 = ph9_p6[k];

        pc_79[k] = - e_0 * f_105_4 * h1_0 - e_1 * f_63_4 * h3_0 + e_1 * f_42_1 * r_2 * h1_0 + e_2 * f_15_2 * h5_0 + e_2 * f_21_2 * r_2 * h3_0 - e_2 * f_18_1 * r_4 * h1_0 + e_3 * f_665_429 * h7_0 + e_3 * fs_25_1716_6006 * h7_p6 - e_3 * f_30_13 * r_2 * h5_0 - e_3 * f_21_11 * r_4 * h3_0 + e_3 * f_8_3 * r_6 * h1_0 + e_4 * f_126_2431 * h9_0 - e_4 * fs_21_4862_4290 * h9_p6 - e_4 * f_1330_7293 * r_2 * h7_0 - e_4 * fs_25_14586_6006 * r_2 * h7_p6 + e_4 * f_2_13 * r_4 * h5_0 + e_4 * f_14_143 * r_6 * h3_0 - e_4 * f_4_33 * r_8 * h1_0;
    }

    // NOTE: the angular components are formed in 67 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph1_p1, ph3_p1, ph5_p1, ph7_p1, ph7_p7, ph9_m8, ph9_p1, ph9_p7, ab_2, pc_80, pc_81 : simd::cache_line_size())
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

        const auto h1_p1 = ph1_p1[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h7_p7 = ph7_p7[k];
        const auto h9_m8 = ph9_m8[k];
        const auto h9_p1 = ph9_p1[k];
        const auto h9_p7 = ph9_p7[k];

        pc_80[k] = e_0 * fs_105_32_2 * h1_p1 + e_1 * fs_21_2_3 * h3_p1 - e_1 * fs_21_4_2 * r_2 * h1_p1 + e_2 * fs_5_4_30 * h5_p1 - e_2 * fs_7_1_3 * r_2 * h3_p1 + e_2 * fs_9_4_2 * r_4 * h1_p1 + e_3 * fs_70_429_14 * h7_p1 - e_3 * fs_35_858_858 * h7_p7 - e_3 * fs_5_13_30 * r_2 * h5_p1 + e_3 * fs_14_11_3 * r_4 * h3_p1 - e_3 * fs_1_3_2 * r_6 * h1_p1 + e_4 * fs_21_4862_10 * h9_p1 - e_4 * fs_21_2431_715 * h9_p7 - e_4 * fs_140_7293_14 * r_2 * h7_p1 + e_4 * fs_35_7293_858 * r_2 * h7_p7 + e_4 * fs_1_39_30 * r_4 * h5_p1 - e_4 * fs_28_429_3 * r_6 * h3_p1 + e_4 * fs_1_66_2 * r_8 * h1_p1;

        pc_81[k] = - e_4 * fs_7_2431_12155 * h9_m8;
    }

    // NOTE: the angular components are formed in 67 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph1_m1, ph3_m1, ph5_m1, ph7_m7, ph7_m1, ph9_m7, ph9_m1, ab_2, pc_82 : simd::cache_line_size())
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

        const auto h1_m1 = ph1_m1[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h7_m7 = ph7_m7[k];
        const auto h7_m1 = ph7_m1[k];
        const auto h9_m7 = ph9_m7[k];
        const auto h9_m1 = ph9_m1[k];

        pc_82[k] = e_0 * fs_315_16_2 * h1_m1 - e_1 * fs_21_8_3 * h3_m1 - e_1 * fs_63_2_2 * r_2 * h1_m1 - e_2 * fs_5_4_30 * h5_m1 + e_2 * fs_7_4_3 * r_2 * h3_m1 + e_2 * fs_27_2_2 * r_4 * h1_m1 + e_3 * fs_35_572_858 * h7_m7 - e_3 * fs_115_572_14 * h7_m1 + e_3 * fs_5_13_30 * r_2 * h5_m1 - e_3 * fs_7_22_3 * r_4 * h3_m1 - e_3 * fs_2_1_2 * r_6 * h1_m1 - e_4 * fs_28_2431_715 * h9_m7 - e_4 * fs_14_2431_10 * h9_m1 - e_4 * fs_35_4862_858 * r_2 * h7_m7 + e_4 * fs_115_4862_14 * r_2 * h7_m1 - e_4 * fs_1_39_30 * r_4 * h5_m1 + e_4 * fs_7_429_3 * r_6 * h3_m1 + e_4 * fs_1_11_2 * r_8 * h1_m1;
    }

    // NOTE: the angular components are formed in 67 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, ph3_m3, ph3_m2, ph5_m5, ph5_m2, ph7_m6, ph7_m5, ph7_m3, ph7_m2, ph9_m6, ph9_m5, ph9_m3, ph9_m2, ab_2, pc_83, pc_84 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h3_m3 = ph3_m3[k];
        const auto h3_m2 = ph3_m2[k];
        const auto h5_m5 = ph5_m5[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h7_m6 = ph7_m6[k];
        const auto h7_m5 = ph7_m5[k];
        const auto h7_m3 = ph7_m3[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h9_m6 = ph9_m6[k];
        const auto h9_m5 = ph9_m5[k];
        const auto h9_m3 = ph9_m3[k];
        const auto h9_m2 = ph9_m2[k];

        pc_83[k] = - e_1 * fs_3_2_105 * h3_m2 + e_2 * fs_5_4_15 * h5_m2 + e_2 * fs_1_1_105 * r_2 * h3_m2 + e_3 * fs_5_286_858 * h7_m6 + e_3 * fs_135_286_6 * h7_m2 - e_3 * fs_5_13_15 * r_2 * h5_m2 - e_3 * fs_2_11_105 * r_4 * h3_m2 - e_4 * fs_7_4862_30030 * h9_m6 + e_4 * fs_7_4862_770 * h9_m2 - e_4 * fs_5_2431_858 * r_2 * h7_m6 - e_4 * fs_135_2431_6 * r_2 * h7_m2 + e_4 * fs_1_39_15 * r_4 * h5_m2 + e_4 * fs_4_429_105 * r_6 * h3_m2;

        pc_84[k] = e_1 * fs_27_8_35 * h3_m3 + e_2 * f_15_2 * h5_m5 - e_2 * fs_9_4_35 * r_2 * h3_m3 - e_3 * fs_45_572_66 * h7_m5 - e_3 * fs_25_44_6 * h7_m3 - e_3 * f_30_13 * r_2 * h5_m5 + e_3 * fs_9_22_35 * r_4 * h3_m3 - e_4 * fs_14_2431_1001 * h9_m5 - e_4 * fs_14_2431_165 * h9_m3 + e_4 * fs_45_4862_66 * r_2 * h7_m5 + e_4 * fs_25_374_6 * r_2 * h7_m3 + e_4 * f_2_13 * r_4 * h5_m5 - e_4 * fs_3_143_35 * r_6 * h3_m3;
    }

    // NOTE: the angular components are formed in 67 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, ph3_p3, ph5_p4, ph5_p5, ph7_p3, ph7_p4, ph7_p5, ph9_p3, ph9_p4, ph9_p5, ab_2, pc_85, pc_86 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h3_p3 = ph3_p3[k];
        const auto h5_p4 = ph5_p4[k];
        const auto h5_p5 = ph5_p5[k];
        const auto h7_p3 = ph7_p3[k];
        const auto h7_p4 = ph7_p4[k];
        const auto h7_p5 = ph7_p5[k];
        const auto h9_p3 = ph9_p3[k];
        const auto h9_p4 = ph9_p4[k];
        const auto h9_p5 = ph9_p5[k];

        pc_85[k] = e_2 * f_15_2 * h5_p4 - e_3 * fs_20_143_165 * h7_p4 - e_3 * f_30_13 * r_2 * h5_p4 - e_4 * fs_35_2431_143 * h9_p4 + e_4 * fs_40_2431_165 * r_2 * h7_p4 + e_4 * f_2_13 * r_4 * h5_p4;

        pc_86[k] = - e_1 * fs_27_8_35 * h3_p3 + e_2 * f_15_2 * h5_p5 + e_2 * fs_9_4_35 * r_2 * h3_p3 + e_3 * fs_25_44_6 * h7_p3 - e_3 * fs_45_572_66 * h7_p5 - e_3 * f_30_13 * r_2 * h5_p5 - e_3 * fs_9_22_35 * r_4 * h3_p3 + e_4 * fs_14_2431_165 * h9_p3 - e_4 * fs_14_2431_1001 * h9_p5 - e_4 * fs_25_374_6 * r_2 * h7_p3 + e_4 * fs_45_4862_66 * r_2 * h7_p5 + e_4 * f_2_13 * r_4 * h5_p5 + e_4 * fs_3_143_35 * r_6 * h3_p3;
    }

    // NOTE: the angular components are formed in 67 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, ph3_p2, ph5_p2, ph7_p2, ph7_p6, ph9_p2, ph9_p6, ab_2, pc_87 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h3_p2 = ph3_p2[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h7_p6 = ph7_p6[k];
        const auto h9_p2 = ph9_p2[k];
        const auto h9_p6 = ph9_p6[k];

        pc_87[k] = e_1 * fs_3_2_105 * h3_p2 - e_2 * fs_5_4_15 * h5_p2 - e_2 * fs_1_1_105 * r_2 * h3_p2 - e_3 * fs_135_286_6 * h7_p2 + e_3 * fs_5_286_858 * h7_p6 + e_3 * fs_5_13_15 * r_2 * h5_p2 + e_3 * fs_2_11_105 * r_4 * h3_p2 - e_4 * fs_7_4862_770 * h9_p2 - e_4 * fs_7_4862_30030 * h9_p6 + e_4 * fs_135_2431_6 * r_2 * h7_p2 - e_4 * fs_5_2431_858 * r_2 * h7_p6 - e_4 * fs_1_39_15 * r_4 * h5_p2 - e_4 * fs_4_429_105 * r_6 * h3_p2;
    }

    // NOTE: the angular components are formed in 67 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph1_p1, ph3_p1, ph5_p1, ph7_p1, ph7_p7, ph9_p1, ph9_p7, ab_2, pc_88 : simd::cache_line_size())
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

        const auto h1_p1 = ph1_p1[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h7_p7 = ph7_p7[k];
        const auto h9_p1 = ph9_p1[k];
        const auto h9_p7 = ph9_p7[k];

        pc_88[k] = - e_0 * fs_315_16_2 * h1_p1 + e_1 * fs_21_8_3 * h3_p1 + e_1 * fs_63_2_2 * r_2 * h1_p1 + e_2 * fs_5_4_30 * h5_p1 - e_2 * fs_7_4_3 * r_2 * h3_p1 - e_2 * fs_27_2_2 * r_4 * h1_p1 + e_3 * fs_115_572_14 * h7_p1 + e_3 * fs_35_572_858 * h7_p7 - e_3 * fs_5_13_30 * r_2 * h5_p1 + e_3 * fs_7_22_3 * r_4 * h3_p1 + e_3 * fs_2_1_2 * r_6 * h1_p1 + e_4 * fs_14_2431_10 * h9_p1 - e_4 * fs_28_2431_715 * h9_p7 - e_4 * fs_115_4862_14 * r_2 * h7_p1 - e_4 * fs_35_4862_858 * r_2 * h7_p7 + e_4 * fs_1_39_30 * r_4 * h5_p1 - e_4 * fs_7_429_3 * r_6 * h3_p1 - e_4 * fs_1_11_2 * r_8 * h1_p1;
    }

    // NOTE: the angular components are formed in 67 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph1_m1, ph1_0, ph3_m1, ph3_0, ph5_m1, ph5_0, ph7_m1, ph7_0, ph9_m9, ph9_m1, ph9_0, ph9_p8, ab_2, pc_89, pc_90 : simd::cache_line_size())
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

        const auto h1_m1 = ph1_m1[k];
        const auto h1_0 = ph1_0[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h3_0 = ph3_0[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h5_0 = ph5_0[k];
        const auto h7_m1 = ph7_m1[k];
        const auto h7_0 = ph7_0[k];
        const auto h9_m9 = ph9_m9[k];
        const auto h9_m1 = ph9_m1[k];
        const auto h9_0 = ph9_0[k];
        const auto h9_p8 = ph9_p8[k];

        pc_89[k] = - e_0 * f_315_16 * h1_0 - e_1 * f_63_2 * h3_0 + e_1 * f_63_2 * r_2 * h1_0 - e_2 * f_15_2 * h5_0 + e_2 * f_21_1 * r_2 * h3_0 - e_2 * f_27_2 * r_4 * h1_0 - e_3 * f_70_143 * h7_0 + e_3 * f_30_13 * r_2 * h5_0 - e_3 * f_42_11 * r_4 * h3_0 + e_3 * f_2_1 * r_6 * h1_0 - e_4 * f_21_2431 * h9_0 - e_4 * fs_7_2431_12155 * h9_p8 + e_4 * f_140_2431 * r_2 * h7_0 - e_4 * f_2_13 * r_4 * h5_0 + e_4 * f_28_143 * r_6 * h3_0 - e_4 * f_1_11 * r_8 * h1_0;

        pc_90[k] = e_0 * fs_315_32_10 * h1_m1 + e_1 * fs_21_4_15 * h3_m1 - e_1 * fs_63_4_10 * r_2 * h1_m1 + e_2 * fs_5_4_6 * h5_m1 - e_2 * fs_7_2_15 * r_2 * h3_m1 + e_2 * fs_27_4_10 * r_4 * h1_m1 + e_3 * fs_5_286_70 * h7_m1 - e_3 * fs_5_13_6 * r_2 * h5_m1 + e_3 * fs_7_11_15 * r_4 * h3_m1 - e_3 * fs_1_1_10 * r_6 * h1_m1 - e_4 * fs_21_2431_2431 * h9_m9 + e_4 * fs_7_4862_2 * h9_m1 - e_4 * fs_5_2431_70 * r_2 * h7_m1 + e_4 * fs_1_39_6 * r_4 * h5_m1 - e_4 * fs_14_429_15 * r_6 * h3_m1 + e_4 * fs_1_22_10 * r_8 * h1_m1;
    }

    // NOTE: the angular components are formed in 67 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, ph3_m3, ph3_m2, ph5_m3, ph5_m2, ph7_m7, ph7_m3, ph7_m2, ph9_m8, ph9_m7, ph9_m3, ph9_m2, ab_2, pc_91, pc_92 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h3_m3 = ph3_m3[k];
        const auto h3_m2 = ph3_m2[k];
        const auto h5_m3 = ph5_m3[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h7_m7 = ph7_m7[k];
        const auto h7_m3 = ph7_m3[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h9_m8 = ph9_m8[k];
        const auto h9_m7 = ph9_m7[k];
        const auto h9_m3 = ph9_m3[k];
        const auto h9_m2 = ph9_m2[k];

        pc_91[k] = - e_1 * fs_105_8_3 * h3_m2 - e_2 * fs_5_4_21 * h5_m2 + e_2 * fs_35_4_3 * r_2 * h3_m2 - e_3 * fs_15_572_210 * h7_m2 + e_3 * fs_5_13_21 * r_2 * h5_m2 - e_3 * fs_35_22_3 * r_4 * h3_m2 - e_4 * fs_14_2431_2431 * h9_m8 - e_4 * fs_7_4862_22 * h9_m2 + e_4 * fs_15_4862_210 * r_2 * h7_m2 - e_4 * fs_1_39_21 * r_4 * h5_m2 + e_4 * fs_35_429_3 * r_6 * h3_m2;

        pc_92[k] = e_1 * fs_45_8_7 * h3_m3 + e_2 * f_15_2 * h5_m3 - e_2 * fs_15_4_7 * r_2 * h3_m3 - e_3 * fs_5_572_30030 * h7_m7 + e_3 * fs_75_572_30 * h7_m3 - e_3 * f_30_13 * r_2 * h5_m3 + e_3 * fs_15_22_7 * r_4 * h3_m3 - e_4 * fs_14_2431_1001 * h9_m7 + e_4 * fs_7_2431_33 * h9_m3 + e_4 * fs_5_4862_30030 * r_2 * h7_m7 - e_4 * fs_75_4862_30 * r_2 * h7_m3 + e_4 * f_2_13 * r_4 * h5_m3 - e_4 * fs_5_143_7 * r_6 * h3_m3;
    }

    // NOTE: the angular components are formed in 67 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_2, pe_3, pe_4, ph5_m4, ph5_p4, ph5_p5, ph7_m6, ph7_m4, ph7_p4, ph7_p5, ph7_p6, ph9_m6, ph9_m4, ph9_p4, ph9_p5, ph9_p6, ab_2, pc_93, pc_94, pc_95 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

        const auto h5_m4 = ph5_m4[k];
        const auto h5_p4 = ph5_p4[k];
        const auto h5_p5 = ph5_p5[k];
        const auto h7_m6 = ph7_m6[k];
        const auto h7_m4 = ph7_m4[k];
        const auto h7_p4 = ph7_p4[k];
        const auto h7_p5 = ph7_p5[k];
        const auto h7_p6 = ph7_p6[k];
        const auto h9_m6 = ph9_m6[k];
        const auto h9_m4 = ph9_m4[k];
        const auto h9_p4 = ph9_p4[k];
        const auto h9_p5 = ph9_p5[k];
        const auto h9_p6 = ph9_p6[k];

        pc_93[k] = - e_2 * f_15_2 * h5_m4 - e_3 * fs_15_572_4290 * h7_m6 - e_3 * fs_25_286_165 * h7_m4 + e_3 * f_30_13 * r_2 * h5_m4 - e_4 * fs_7_4862_6006 * h9_m6 - e_4 * fs_7_2431_143 * h9_m4 + e_4 * fs_15_4862_4290 * r_2 * h7_m6 + e_4 * fs_25_2431_165 * r_2 * h7_m4 - e_4 * f_2_13 * r_4 * h5_m4;

        pc_94[k] = - e_2 * f_15_2 * h5_p5 - e_3 * fs_75_286_66 * h7_p5 + e_3 * f_30_13 * r_2 * h5_p5 - e_4 * fs_7_2431_1001 * h9_p5 + e_4 * fs_75_2431_66 * r_2 * h7_p5 - e_4 * f_2_13 * r_4 * h5_p5;

        pc_95[k] = e_2 * f_15_2 * h5_p4 + e_3 * fs_25_286_165 * h7_p4 - e_3 * fs_15_572_4290 * h7_p6 - e_3 * f_30_13 * r_2 * h5_p4 + e_4 * fs_7_2431_143 * h9_p4 - e_4 * fs_7_4862_6006 * h9_p6 - e_4 * fs_25_2431_165 * r_2 * h7_p4 + e_4 * fs_15_4862_4290 * r_2 * h7_p6 + e_4 * f_2_13 * r_4 * h5_p4;
    }

    // NOTE: the angular components are formed in 67 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, ph3_p2, ph3_p3, ph5_p2, ph5_p3, ph7_p2, ph7_p3, ph7_p7, ph9_p2, ph9_p3, ph9_p7, ph9_p8, ab_2, pc_96, pc_97 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h3_p2 = ph3_p2[k];
        const auto h3_p3 = ph3_p3[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h7_p3 = ph7_p3[k];
        const auto h7_p7 = ph7_p7[k];
        const auto h9_p2 = ph9_p2[k];
        const auto h9_p3 = ph9_p3[k];
        const auto h9_p7 = ph9_p7[k];
        const auto h9_p8 = ph9_p8[k];

        pc_96[k] = - e_1 * fs_45_8_7 * h3_p3 - e_2 * f_15_2 * h5_p3 + e_2 * fs_15_4_7 * r_2 * h3_p3 - e_3 * fs_75_572_30 * h7_p3 - e_3 * fs_5_572_30030 * h7_p7 + e_3 * f_30_13 * r_2 * h5_p3 - e_3 * fs_15_22_7 * r_4 * h3_p3 - e_4 * fs_7_2431_33 * h9_p3 - e_4 * fs_14_2431_1001 * h9_p7 + e_4 * fs_75_4862_30 * r_2 * h7_p3 + e_4 * fs_5_4862_30030 * r_2 * h7_p7 - e_4 * f_2_13 * r_4 * h5_p3 + e_4 * fs_5_143_7 * r_6 * h3_p3;

        pc_97[k] = e_1 * fs_105_8_3 * h3_p2 + e_2 * fs_5_4_21 * h5_p2 - e_2 * fs_35_4_3 * r_2 * h3_p2 + e_3 * fs_15_572_210 * h7_p2 - e_3 * fs_5_13_21 * r_2 * h5_p2 + e_3 * fs_35_22_3 * r_4 * h3_p2 + e_4 * fs_7_4862_22 * h9_p2 - e_4 * fs_14_2431_2431 * h9_p8 - e_4 * fs_15_4862_210 * r_2 * h7_p2 + e_4 * fs_1_39_21 * r_4 * h5_p2 - e_4 * fs_35_429_3 * r_6 * h3_p2;
    }

    // NOTE: the angular components are formed in 67 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph1_p1, ph3_p1, ph5_p1, ph7_p1, ph9_p1, ph9_p9, ab_2, pc_98 : simd::cache_line_size())
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

        const auto h1_p1 = ph1_p1[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h9_p1 = ph9_p1[k];
        const auto h9_p9 = ph9_p9[k];

        pc_98[k] = - e_0 * fs_315_32_10 * h1_p1 - e_1 * fs_21_4_15 * h3_p1 + e_1 * fs_63_4_10 * r_2 * h1_p1 - e_2 * fs_5_4_6 * h5_p1 + e_2 * fs_7_2_15 * r_2 * h3_p1 - e_2 * fs_27_4_10 * r_4 * h1_p1 - e_3 * fs_5_286_70 * h7_p1 + e_3 * fs_5_13_6 * r_2 * h5_p1 - e_3 * fs_7_11_15 * r_4 * h3_p1 + e_3 * fs_1_1_10 * r_6 * h1_p1 - e_4 * fs_7_4862_2 * h9_p1 - e_4 * fs_21_2431_2431 * h9_p9 + e_4 * fs_5_2431_70 * r_2 * h7_p1 - e_4 * fs_1_39_6 * r_4 * h5_p1 + e_4 * fs_14_429_15 * r_6 * h3_p1 - e_4 * fs_1_22_10 * r_8 * h1_p1;
    }

    // NOTE: the values of a combination of angular components are stored as one
    // row of nvalues columns, with the component on bra side running slowest. The
    // rows which the symmetry relates to an already formed one are copied from it.

    const size_t sources[99] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98};

    for (size_t n = 0; n < 99; n++)
    {
        if (sources[n] != n) std::copy(values + sources[n] * nvalues, values + (sources[n] + 1) * nvalues, values + n * nvalues);
    }
}

}  // namespace simdt2ceri
