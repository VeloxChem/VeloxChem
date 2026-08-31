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



#include "SimdTwoCenterElectronRepulsionRecFI.hpp"

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
compute_fi_electron_repulsion(double                         *values,
                              const size_t                    nvalues,
                              const CBasisFunction           &bra,
                              const CBasisFunction           &ket,
                              const std::vector<CSimdMatrix> &harmonics,
                              const CSimdMatrix              &coordinates) -> void
{
    if ((bra.get_angular_momentum() != 3) || (ket.get_angular_momentum() != 6))
    {
        errors::assertMsgCritical(
            false, std::string("SimdTwoCenterElectronRepulsionFunc.compute_fi_electron_repulsion: Basis functions must be of angular momenta three and six"));
    }

    if (harmonics.size() < 9)
    {
        errors::assertMsgCritical(
            false, std::string("SimdTwoCenterElectronRepulsionFunc.compute_fi_electron_repulsion: Harmonics must reach angular momentum 9"));
    }

    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("SimdTwoCenterElectronRepulsionFunc.compute_fi_electron_repulsion: Number of values exceeds number of atom pairs"));
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
    // orders 6 to 9 alone, and the orders below them are formed on the
    // way to them by the recursion the Boys function is evaluated with.

    simdfunc::compute_boys_function(boys);

    // NOTE: the buffer holds the contracted prefactor of each exponent factor of
    // the terms, as the integrals of the angular components are formed straight
    // into the values and are not written a second time. Every exponent factor is
    // used with one order of the Boys function alone, so the order follows from
    // the factor and one accumulator per factor suffices.

    auto buffer = CSimdMatrix(4, nvalues);

    auto *pe_0 = buffer.data(0);
    auto *pe_1 = buffer.data(1);
    auto *pe_2 = buffer.data(2);
    auto *pe_3 = buffer.data(3);

    std::fill(pe_0, pe_0 + nvalues, 0.0);
    std::fill(pe_1, pe_1 + nvalues, 0.0);
    std::fill(pe_2, pe_2 + nvalues, 0.0);
    std::fill(pe_3, pe_3 + nvalues, 0.0);

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

            const auto ff_0 = fbase * aexp * aexp * aexp / (fexp * fexp * fexp * fexp * fexp * fexp);

            const auto ff_1 = fbase * aexp * aexp * aexp * fmu / (fexp * fexp * fexp * fexp * fexp * fexp);

            const auto ff_2 = fbase * aexp * aexp * aexp * fmu * fmu / (fexp * fexp * fexp * fexp * fexp * fexp);

            const auto ff_3 = fbase * aexp * aexp * aexp * fmu * fmu * fmu / (fexp * fexp * fexp * fexp * fexp * fexp);

            const auto *bv_0 = boys.data(7, i * nprim_b + j);

            const auto *bv_1 = boys.data(8, i * nprim_b + j);

            const auto *bv_2 = boys.data(9, i * nprim_b + j);

            const auto *bv_3 = boys.data(10, i * nprim_b + j);

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, bv_0, bv_1, bv_2, bv_3 : simd::cache_line_size())
            for (size_t k = 0; k < nvalues; k++)
            {
                pe_0[k] += ff_0 * bv_0[k];
                pe_1[k] += ff_1 * bv_1[k];
                pe_2[k] += ff_2 * bv_2[k];
                pe_3[k] += ff_3 * bv_3[k];
            }
        }
    }

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

    // NOTE: the factors of the terms depend on the angular momenta alone, so they
    // are formed once for the whole matrix instead of once for every atom pair.

    const auto f_100_429 = 100.0 / 429.0;
    const auto f_10_143 = 10.0 / 143.0;
    const auto f_15_11 = 15.0 / 11.0;
    const auto f_15_13 = 15.0 / 13.0;
    const auto f_15_2 = 15.0 / 2.0;
    const auto f_15_4 = 15.0 / 4.0;
    const auto f_15_8 = 15.0 / 8.0;
    const auto f_1_13 = 1.0 / 13.0;
    const auto f_225_8 = 225.0 / 8.0;
    const auto f_252_143 = 252.0 / 143.0;
    const auto f_25_1 = 25.0;
    const auto f_25_143 = 25.0 / 143.0;
    const auto f_2_13 = 2.0 / 13.0;
    const auto f_30_13 = 30.0 / 13.0;
    const auto f_35_13 = 35.0 / 13.0;
    const auto f_35_4 = 35.0 / 4.0;
    const auto f_3_17 = 3.0 / 17.0;
    const auto f_3_2 = 3.0 / 2.0;
    const auto f_45_4 = 45.0 / 4.0;
    const auto f_504_2431 = 504.0 / 2431.0;
    const auto f_50_11 = 50.0 / 11.0;
    const auto f_5_22 = 5.0 / 22.0;
    const auto f_5_4 = 5.0 / 4.0;
    const auto f_5_429 = 5.0 / 429.0;
    const auto f_6_221 = 6.0 / 221.0;
    const auto f_75_2 = 75.0 / 2.0;
    const auto f_75_22 = 75.0 / 22.0;
    const auto f_75_4 = 75.0 / 4.0;
    const auto f_7_39 = 7.0 / 39.0;
    const auto f_84_221 = 84.0 / 221.0;
    const auto fs_105_143_3 = std::sqrt(33075.0 / 20449.0);
    const auto fs_10_11_14 = std::sqrt(1400.0 / 121.0);
    const auto fs_10_13_2 = std::sqrt(200.0 / 169.0);
    const auto fs_10_13_21 = std::sqrt(2100.0 / 169.0);
    const auto fs_10_13_7 = std::sqrt(700.0 / 169.0);
    const auto fs_10_143_7 = std::sqrt(700.0 / 20449.0);
    const auto fs_10_429_21 = std::sqrt(700.0 / 61347.0);
    const auto fs_111_4862_30 = std::sqrt(184815.0 / 11819522.0);
    const auto fs_111_572_30 = std::sqrt(184815.0 / 163592.0);
    const auto fs_120_2431_10 = std::sqrt(144000.0 / 5909761.0);
    const auto fs_126_2431_5 = std::sqrt(79380.0 / 5909761.0);
    const auto fs_129_143_3 = std::sqrt(49923.0 / 20449.0);
    const auto fs_12_221_7 = std::sqrt(1008.0 / 48841.0);
    const auto fs_12_2431_1430 = std::sqrt(1440.0 / 41327.0);
    const auto fs_135_2431_11 = std::sqrt(18225.0 / 537251.0);
    const auto fs_135_286_11 = std::sqrt(18225.0 / 7436.0);
    const auto fs_135_4862_10 = std::sqrt(91125.0 / 11819522.0);
    const auto fs_135_572_10 = std::sqrt(91125.0 / 163592.0);
    const auto fs_147_2431_5 = std::sqrt(108045.0 / 5909761.0);
    const auto fs_147_286_5 = std::sqrt(108045.0 / 81796.0);
    const auto fs_14_221_30 = std::sqrt(5880.0 / 48841.0);
    const auto fs_150_2431_14 = std::sqrt(315000.0 / 5909761.0);
    const auto fs_15_11_7 = std::sqrt(1575.0 / 121.0);
    const auto fs_15_13_11 = std::sqrt(2475.0 / 169.0);
    const auto fs_15_13_3 = std::sqrt(675.0 / 169.0);
    const auto fs_15_13_5 = std::sqrt(1125.0 / 169.0);
    const auto fs_15_13_7 = std::sqrt(1575.0 / 169.0);
    const auto fs_15_16_14 = std::sqrt(1575.0 / 128.0);
    const auto fs_15_16_210 = std::sqrt(23625.0 / 128.0);
    const auto fs_15_221_7 = std::sqrt(1575.0 / 48841.0);
    const auto fs_15_22_21 = std::sqrt(4725.0 / 484.0);
    const auto fs_15_2431_15 = std::sqrt(3375.0 / 5909761.0);
    const auto fs_15_2431_66 = std::sqrt(1350.0 / 537251.0);
    const auto fs_15_26_10 = std::sqrt(1125.0 / 338.0);
    const auto fs_15_26_22 = std::sqrt(2475.0 / 338.0);
    const auto fs_15_286_15 = std::sqrt(3375.0 / 81796.0);
    const auto fs_15_286_66 = std::sqrt(675.0 / 3718.0);
    const auto fs_15_2_14 = std::sqrt(1575.0 / 2.0);
    const auto fs_15_2_7 = std::sqrt(1575.0 / 4.0);
    const auto fs_15_442_66 = std::sqrt(7425.0 / 97682.0);
    const auto fs_15_4862_330 = std::sqrt(3375.0 / 1074502.0);
    const auto fs_15_4_11 = std::sqrt(2475.0 / 16.0);
    const auto fs_15_4_21 = std::sqrt(4725.0 / 16.0);
    const auto fs_15_4_3 = std::sqrt(675.0 / 16.0);
    const auto fs_15_4_5 = std::sqrt(1125.0 / 16.0);
    const auto fs_15_4_7 = std::sqrt(1575.0 / 16.0);
    const auto fs_15_572_330 = std::sqrt(3375.0 / 14872.0);
    const auto fs_15_8_10 = std::sqrt(1125.0 / 32.0);
    const auto fs_15_8_105 = std::sqrt(23625.0 / 64.0);
    const auto fs_15_8_14 = std::sqrt(1575.0 / 32.0);
    const auto fs_15_8_210 = std::sqrt(23625.0 / 32.0);
    const auto fs_15_8_22 = std::sqrt(2475.0 / 32.0);
    const auto fs_15_8_231 = std::sqrt(51975.0 / 64.0);
    const auto fs_15_8_462 = std::sqrt(51975.0 / 32.0);
    const auto fs_177_2431_3 = std::sqrt(93987.0 / 5909761.0);
    const auto fs_177_286_3 = std::sqrt(93987.0 / 81796.0);
    const auto fs_189_2431_5 = std::sqrt(178605.0 / 5909761.0);
    const auto fs_189_286_5 = std::sqrt(178605.0 / 81796.0);
    const auto fs_18_143_70 = std::sqrt(22680.0 / 20449.0);
    const auto fs_18_221_13 = std::sqrt(324.0 / 3757.0);
    const auto fs_18_221_14 = std::sqrt(4536.0 / 48841.0);
    const auto fs_18_221_3 = std::sqrt(972.0 / 48841.0);
    const auto fs_1_13_11 = std::sqrt(11.0 / 169.0);
    const auto fs_1_13_3 = std::sqrt(3.0 / 169.0);
    const auto fs_1_13_5 = std::sqrt(5.0 / 169.0);
    const auto fs_1_13_7 = std::sqrt(7.0 / 169.0);
    const auto fs_1_221_105 = std::sqrt(105.0 / 48841.0);
    const auto fs_1_221_14 = std::sqrt(14.0 / 48841.0);
    const auto fs_1_221_195 = std::sqrt(15.0 / 3757.0);
    const auto fs_1_221_231 = std::sqrt(231.0 / 48841.0);
    const auto fs_1_221_3003 = std::sqrt(231.0 / 3757.0);
    const auto fs_1_221_3094 = std::sqrt(14.0 / 221.0);
    const auto fs_1_221_455 = std::sqrt(35.0 / 3757.0);
    const auto fs_1_221_858 = std::sqrt(66.0 / 3757.0);
    const auto fs_1_221_910 = std::sqrt(70.0 / 3757.0);
    const auto fs_1_221_9282 = std::sqrt(42.0 / 221.0);
    const auto fs_1_26_10 = std::sqrt(5.0 / 338.0);
    const auto fs_1_26_22 = std::sqrt(11.0 / 338.0);
    const auto fs_1_39_14 = std::sqrt(14.0 / 1521.0);
    const auto fs_1_39_15 = std::sqrt(5.0 / 507.0);
    const auto fs_1_39_21 = std::sqrt(7.0 / 507.0);
    const auto fs_1_39_3 = std::sqrt(1.0 / 507.0);
    const auto fs_1_39_33 = std::sqrt(11.0 / 507.0);
    const auto fs_1_39_35 = std::sqrt(35.0 / 1521.0);
    const auto fs_1_39_42 = std::sqrt(14.0 / 507.0);
    const auto fs_1_39_6 = std::sqrt(2.0 / 507.0);
    const auto fs_1_442_10010 = std::sqrt(385.0 / 7514.0);
    const auto fs_1_442_14 = std::sqrt(7.0 / 97682.0);
    const auto fs_1_442_182 = std::sqrt(7.0 / 7514.0);
    const auto fs_1_442_2 = std::sqrt(1.0 / 97682.0);
    const auto fs_1_442_26 = std::sqrt(1.0 / 7514.0);
    const auto fs_1_442_462 = std::sqrt(231.0 / 97682.0);
    const auto fs_1_442_6006 = std::sqrt(231.0 / 7514.0);
    const auto fs_1_78_30 = std::sqrt(5.0 / 1014.0);
    const auto fs_1_78_6 = std::sqrt(1.0 / 1014.0);
    const auto fs_1_78_66 = std::sqrt(11.0 / 1014.0);
    const auto fs_20_13_2 = std::sqrt(800.0 / 169.0);
    const auto fs_20_221_11 = std::sqrt(4400.0 / 48841.0);
    const auto fs_20_429_14 = std::sqrt(5600.0 / 184041.0);
    const auto fs_210_2431_3 = std::sqrt(132300.0 / 5909761.0);
    const auto fs_21_221_6 = std::sqrt(2646.0 / 48841.0);
    const auto fs_21_26_6 = std::sqrt(1323.0 / 338.0);
    const auto fs_24_143_33 = std::sqrt(1728.0 / 1859.0);
    const auto fs_258_2431_3 = std::sqrt(199692.0 / 5909761.0);
    const auto fs_25_22_14 = std::sqrt(4375.0 / 242.0);
    const auto fs_25_429_14 = std::sqrt(8750.0 / 184041.0);
    const auto fs_25_4_14 = std::sqrt(4375.0 / 8.0);
    const auto fs_27_2431_165 = std::sqrt(10935.0 / 537251.0);
    const auto fs_27_286_165 = std::sqrt(10935.0 / 7436.0);
    const auto fs_2_221_1001 = std::sqrt(308.0 / 3757.0);
    const auto fs_2_221_1547 = std::sqrt(28.0 / 221.0);
    const auto fs_2_221_21 = std::sqrt(84.0 / 48841.0);
    const auto fs_2_221_231 = std::sqrt(924.0 / 48841.0);
    const auto fs_2_221_70 = std::sqrt(280.0 / 48841.0);
    const auto fs_2_221_715 = std::sqrt(220.0 / 3757.0);
    const auto fs_2_39_2 = std::sqrt(8.0 / 1521.0);
    const auto fs_2_39_21 = std::sqrt(28.0 / 507.0);
    const auto fs_2_39_7 = std::sqrt(28.0 / 1521.0);
    const auto fs_35_26_2 = std::sqrt(1225.0 / 338.0);
    const auto fs_35_8_2 = std::sqrt(1225.0 / 32.0);
    const auto fs_36_2431_70 = std::sqrt(90720.0 / 5909761.0);
    const auto fs_3_143_14 = std::sqrt(126.0 / 20449.0);
    const auto fs_3_187_42 = std::sqrt(378.0 / 34969.0);
    const auto fs_3_221_273 = std::sqrt(189.0 / 3757.0);
    const auto fs_3_221_30 = std::sqrt(270.0 / 48841.0);
    const auto fs_3_221_385 = std::sqrt(3465.0 / 48841.0);
    const auto fs_3_221_55 = std::sqrt(495.0 / 48841.0);
    const auto fs_3_22_42 = std::sqrt(189.0 / 242.0);
    const auto fs_3_26_273 = std::sqrt(189.0 / 52.0);
    const auto fs_3_26_30 = std::sqrt(135.0 / 338.0);
    const auto fs_3_442_1430 = std::sqrt(495.0 / 7514.0);
    const auto fs_3_442_2002 = std::sqrt(693.0 / 7514.0);
    const auto fs_3_442_26 = std::sqrt(9.0 / 7514.0);
    const auto fs_3_442_910 = std::sqrt(315.0 / 7514.0);
    const auto fs_3_4862_30030 = std::sqrt(945.0 / 82654.0);
    const auto fs_3_52_26 = std::sqrt(9.0 / 104.0);
    const auto fs_3_52_910 = std::sqrt(315.0 / 104.0);
    const auto fs_3_572_30030 = std::sqrt(945.0 / 1144.0);
    const auto fs_45_143_21 = std::sqrt(42525.0 / 20449.0);
    const auto fs_45_2431_110 = std::sqrt(20250.0 / 537251.0);
    const auto fs_45_2431_210 = std::sqrt(425250.0 / 5909761.0);
    const auto fs_45_286_110 = std::sqrt(10125.0 / 3718.0);
    const auto fs_45_286_210 = std::sqrt(212625.0 / 40898.0);
    const auto fs_45_374_2 = std::sqrt(2025.0 / 69938.0);
    const auto fs_45_44_2 = std::sqrt(2025.0 / 968.0);
    const auto fs_45_4862_286 = std::sqrt(2025.0 / 82654.0);
    const auto fs_45_4_7 = std::sqrt(14175.0 / 16.0);
    const auto fs_45_572_286 = std::sqrt(2025.0 / 1144.0);
    const auto fs_45_8_21 = std::sqrt(42525.0 / 64.0);
    const auto fs_48_143_15 = std::sqrt(34560.0 / 20449.0);
    const auto fs_48_2431_33 = std::sqrt(6912.0 / 537251.0);
    const auto fs_4_221_273 = std::sqrt(336.0 / 3757.0);
    const auto fs_4_221_91 = std::sqrt(112.0 / 3757.0);
    const auto fs_4_39_2 = std::sqrt(32.0 / 1521.0);
    const auto fs_573_4862_2 = std::sqrt(328329.0 / 11819522.0);
    const auto fs_573_572_2 = std::sqrt(328329.0 / 163592.0);
    const auto fs_57_4862_110 = std::sqrt(16245.0 / 1074502.0);
    const auto fs_57_572_110 = std::sqrt(16245.0 / 14872.0);
    const auto fs_5_11_21 = std::sqrt(525.0 / 121.0);
    const auto fs_5_13_14 = std::sqrt(350.0 / 169.0);
    const auto fs_5_13_15 = std::sqrt(375.0 / 169.0);
    const auto fs_5_13_21 = std::sqrt(525.0 / 169.0);
    const auto fs_5_13_3 = std::sqrt(75.0 / 169.0);
    const auto fs_5_13_33 = std::sqrt(825.0 / 169.0);
    const auto fs_5_13_35 = std::sqrt(875.0 / 169.0);
    const auto fs_5_13_42 = std::sqrt(1050.0 / 169.0);
    const auto fs_5_13_6 = std::sqrt(150.0 / 169.0);
    const auto fs_5_143_21 = std::sqrt(525.0 / 20449.0);
    const auto fs_5_1_14 = std::sqrt(350.0);
    const auto fs_5_1_2 = std::sqrt(50.0);
    const auto fs_5_221_143 = std::sqrt(275.0 / 3757.0);
    const auto fs_5_221_231 = std::sqrt(5775.0 / 48841.0);
    const auto fs_5_22_105 = std::sqrt(2625.0 / 484.0);
    const auto fs_5_22_14 = std::sqrt(175.0 / 242.0);
    const auto fs_5_22_210 = std::sqrt(2625.0 / 242.0);
    const auto fs_5_22_231 = std::sqrt(525.0 / 44.0);
    const auto fs_5_22_462 = std::sqrt(525.0 / 22.0);
    const auto fs_5_26_30 = std::sqrt(375.0 / 338.0);
    const auto fs_5_26_6 = std::sqrt(75.0 / 338.0);
    const auto fs_5_26_66 = std::sqrt(825.0 / 338.0);
    const auto fs_5_2_2 = std::sqrt(25.0 / 2.0);
    const auto fs_5_2_21 = std::sqrt(525.0 / 4.0);
    const auto fs_5_2_7 = std::sqrt(175.0 / 4.0);
    const auto fs_5_429_105 = std::sqrt(875.0 / 61347.0);
    const auto fs_5_429_14 = std::sqrt(350.0 / 184041.0);
    const auto fs_5_429_210 = std::sqrt(1750.0 / 61347.0);
    const auto fs_5_429_231 = std::sqrt(175.0 / 5577.0);
    const auto fs_5_429_462 = std::sqrt(350.0 / 5577.0);
    const auto fs_5_442_154 = std::sqrt(1925.0 / 97682.0);
    const auto fs_5_44_14 = std::sqrt(175.0 / 968.0);
    const auto fs_5_44_210 = std::sqrt(2625.0 / 968.0);
    const auto fs_5_4_105 = std::sqrt(2625.0 / 16.0);
    const auto fs_5_4_14 = std::sqrt(175.0 / 8.0);
    const auto fs_5_4_15 = std::sqrt(375.0 / 16.0);
    const auto fs_5_4_21 = std::sqrt(525.0 / 16.0);
    const auto fs_5_4_210 = std::sqrt(2625.0 / 8.0);
    const auto fs_5_4_231 = std::sqrt(5775.0 / 16.0);
    const auto fs_5_4_3 = std::sqrt(75.0 / 16.0);
    const auto fs_5_4_33 = std::sqrt(825.0 / 16.0);
    const auto fs_5_4_35 = std::sqrt(875.0 / 16.0);
    const auto fs_5_4_42 = std::sqrt(525.0 / 8.0);
    const auto fs_5_4_462 = std::sqrt(5775.0 / 8.0);
    const auto fs_5_4_6 = std::sqrt(75.0 / 8.0);
    const auto fs_5_858_14 = std::sqrt(175.0 / 368082.0);
    const auto fs_5_858_210 = std::sqrt(875.0 / 122694.0);
    const auto fs_5_8_14 = std::sqrt(175.0 / 32.0);
    const auto fs_5_8_210 = std::sqrt(2625.0 / 32.0);
    const auto fs_5_8_30 = std::sqrt(375.0 / 32.0);
    const auto fs_5_8_6 = std::sqrt(75.0 / 32.0);
    const auto fs_5_8_66 = std::sqrt(825.0 / 32.0);
    const auto fs_60_143_10 = std::sqrt(36000.0 / 20449.0);
    const auto fs_63_143_5 = std::sqrt(19845.0 / 20449.0);
    const auto fs_69_2431_5 = std::sqrt(23805.0 / 5909761.0);
    const auto fs_69_286_5 = std::sqrt(23805.0 / 81796.0);
    const auto fs_6_143_1430 = std::sqrt(360.0 / 143.0);
    const auto fs_6_221_66 = std::sqrt(2376.0 / 48841.0);
    const auto fs_6_2431_14 = std::sqrt(504.0 / 5909761.0);
    const auto fs_75_143_14 = std::sqrt(78750.0 / 20449.0);
    const auto fs_75_4862_66 = std::sqrt(16875.0 / 1074502.0);
    const auto fs_75_572_66 = std::sqrt(16875.0 / 14872.0);
    const auto fs_75_8_14 = std::sqrt(39375.0 / 32.0);
    const auto fs_7_221_66 = std::sqrt(3234.0 / 48841.0);
    const auto fs_7_78_2 = std::sqrt(49.0 / 3042.0);
    const auto fs_87_4862_22 = std::sqrt(7569.0 / 1074502.0);
    const auto fs_87_572_22 = std::sqrt(7569.0 / 14872.0);
    const auto fs_8_221_105 = std::sqrt(6720.0 / 48841.0);
    const auto fs_90_2431_21 = std::sqrt(170100.0 / 5909761.0);
    const auto fs_96_2431_15 = std::sqrt(138240.0 / 5909761.0);
    const auto fs_9_13_13 = std::sqrt(81.0 / 13.0);
    const auto fs_9_13_3 = std::sqrt(243.0 / 169.0);
    const auto fs_9_221_14 = std::sqrt(1134.0 / 48841.0);
    const auto fs_9_2431_10 = std::sqrt(810.0 / 5909761.0);
    const auto fs_9_2431_55 = std::sqrt(405.0 / 537251.0);
    const auto fs_9_286_10 = std::sqrt(405.0 / 40898.0);
    const auto fs_9_286_55 = std::sqrt(405.0 / 7436.0);

    // NOTE: the angular components are formed in 49 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph3_p2, ph3_p3, ph5_p2, ph5_p3, ph7_p2, ph7_p3, ph9_p2, ph9_p3, ph9_p8, ph9_p9, ab_2, pc_0, pc_1 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h3_p2 = ph3_p2[k];
        const auto h3_p3 = ph3_p3[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h7_p3 = ph7_p3[k];
        const auto h9_p2 = ph9_p2[k];
        const auto h9_p3 = ph9_p3[k];
        const auto h9_p8 = ph9_p8[k];
        const auto h9_p9 = ph9_p9[k];

        pc_0[k] = e_0 * fs_15_8_462 * h3_p3 + e_1 * fs_5_8_66 * h5_p3 - e_1 * fs_5_4_462 * r_2 * h3_p3 + e_2 * fs_9_286_55 * h7_p3 - e_2 * fs_5_26_66 * r_2 * h5_p3 + e_2 * fs_5_22_462 * r_4 * h3_p3 + e_3 * fs_1_442_2 * h9_p3 + e_3 * fs_1_221_9282 * h9_p9 - e_3 * fs_9_2431_55 * r_2 * h7_p3 + e_3 * fs_1_78_66 * r_4 * h5_p3 - e_3 * fs_5_429_462 * r_6 * h3_p3;

        pc_1[k] = e_0 * fs_15_8_231 * h3_p2 + e_1 * fs_5_4_33 * h5_p2 - e_1 * fs_5_4_231 * r_2 * h3_p2 + e_2 * fs_15_572_330 * h7_p2 - e_2 * fs_5_13_33 * r_2 * h5_p2 + e_2 * fs_5_22_231 * r_4 * h3_p2 + e_3 * fs_1_442_14 * h9_p2 + e_3 * fs_2_221_1547 * h9_p8 - e_3 * fs_15_4862_330 * r_2 * h7_p2 + e_3 * fs_1_39_33 * r_4 * h5_p2 - e_3 * fs_5_429_231 * r_6 * h3_p2;
    }

    // NOTE: the angular components are formed in 49 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph3_0, ph3_p1, ph5_0, ph5_p1, ph7_0, ph7_p1, ph7_p6, ph7_p7, ph9_0, ph9_p1, ph9_p6, ph9_p7, ab_2, pc_2, pc_3 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h3_0 = ph3_0[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h5_0 = ph5_0[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h7_0 = ph7_0[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h7_p6 = ph7_p6[k];
        const auto h7_p7 = ph7_p7[k];
        const auto h9_0 = ph9_0[k];
        const auto h9_p1 = ph9_p1[k];
        const auto h9_p6 = ph9_p6[k];
        const auto h9_p7 = ph9_p7[k];

        pc_2[k] = e_0 * fs_15_8_105 * h3_p1 + e_1 * fs_5_4_42 * h5_p1 - e_1 * fs_5_4_105 * r_2 * h3_p1 + e_2 * fs_135_572_10 * h7_p1 + e_2 * fs_3_572_30030 * h7_p7 - e_2 * fs_5_13_42 * r_2 * h5_p1 + e_2 * fs_5_22_105 * r_4 * h3_p1 + e_3 * fs_1_221_14 * h9_p1 + e_3 * fs_2_221_1001 * h9_p7 - e_3 * fs_135_4862_10 * r_2 * h7_p1 - e_3 * fs_3_4862_30030 * r_2 * h7_p7 + e_3 * fs_1_39_42 * r_4 * h5_p1 - e_3 * fs_5_429_105 * r_6 * h3_p1;

        pc_3[k] = e_0 * fs_15_4_21 * h3_0 + e_1 * fs_5_2_21 * h5_0 - e_1 * fs_5_2_21 * r_2 * h3_0 + e_2 * fs_45_143_21 * h7_0 + e_2 * fs_45_572_286 * h7_p6 - e_2 * fs_10_13_21 * r_2 * h5_0 + e_2 * fs_5_11_21 * r_4 * h3_0 + e_3 * fs_2_221_21 * h9_0 + e_3 * fs_1_442_10010 * h9_p6 - e_3 * fs_90_2431_21 * r_2 * h7_0 - e_3 * fs_45_4862_286 * r_2 * h7_p6 + e_3 * fs_2_39_21 * r_4 * h5_0 - e_3 * fs_10_429_21 * r_6 * h3_0;
    }

    // NOTE: the angular components are formed in 49 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph3_p1, ph5_p1, ph5_p5, ph7_p1, ph7_p5, ph9_p1, ph9_p5, ab_2, pc_4 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h3_p1 = ph3_p1[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h5_p5 = ph5_p5[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h7_p5 = ph7_p5[k];
        const auto h9_p1 = ph9_p1[k];
        const auto h9_p5 = ph9_p5[k];

        pc_4[k] = - e_0 * fs_15_8_14 * h3_p1 - e_1 * fs_5_4_35 * h5_p1 + e_1 * fs_5_8_6 * h5_p5 + e_1 * fs_5_4_14 * r_2 * h3_p1 - e_2 * fs_105_143_3 * h7_p1 + e_2 * fs_135_286_11 * h7_p5 + e_2 * fs_5_13_35 * r_2 * h5_p1 - e_2 * fs_5_26_6 * r_2 * h5_p5 - e_2 * fs_5_22_14 * r_4 * h3_p1 - e_3 * fs_1_221_105 * h9_p1 + e_3 * fs_1_442_6006 * h9_p5 + e_3 * fs_210_2431_3 * r_2 * h7_p1 - e_3 * fs_135_2431_11 * r_2 * h7_p5 - e_3 * fs_1_39_35 * r_4 * h5_p1 + e_3 * fs_1_78_6 * r_4 * h5_p5 + e_3 * fs_5_429_14 * r_6 * h3_p1;
    }

    // NOTE: the angular components are formed in 49 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph3_m3, ph3_p2, ph5_m3, ph5_p2, ph5_p4, ph7_m3, ph7_p2, ph7_p4, ph9_m3, ph9_p2, ph9_p4, ab_2, pc_5, pc_6 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

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

        pc_5[k] = e_0 * fs_15_16_14 * h3_p2 + e_1 * fs_35_8_2 * h5_p2 + e_1 * fs_5_4_6 * h5_p4 - e_1 * fs_5_8_14 * r_2 * h3_p2 + e_2 * fs_189_286_5 * h7_p2 + e_2 * fs_45_286_110 * h7_p4 - e_2 * fs_35_26_2 * r_2 * h5_p2 - e_2 * fs_5_13_6 * r_2 * h5_p4 + e_2 * fs_5_44_14 * r_4 * h3_p2 + e_3 * fs_1_221_231 * h9_p2 + e_3 * fs_1_221_858 * h9_p4 - e_3 * fs_189_2431_5 * r_2 * h7_p2 - e_3 * fs_45_2431_110 * r_2 * h7_p4 + e_3 * fs_7_78_2 * r_4 * h5_p2 + e_3 * fs_1_39_6 * r_4 * h5_p4 - e_3 * fs_5_858_14 * r_6 * h3_p2;

        pc_6[k] = - e_0 * f_15_8 * h3_m3 - e_1 * fs_5_2_7 * h5_m3 + e_1 * f_5_4 * r_2 * h3_m3 - e_2 * fs_45_286_210 * h7_m3 + e_2 * fs_10_13_7 * r_2 * h5_m3 - e_2 * f_5_22 * r_4 * h3_m3 - e_3 * fs_2_221_231 * h9_m3 + e_3 * fs_45_2431_210 * r_2 * h7_m3 - e_3 * fs_2_39_7 * r_4 * h5_m3 + e_3 * f_5_429 * r_6 * h3_m3;
    }

    // NOTE: the angular components are formed in 49 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph3_m2, ph5_m4, ph5_m2, ph7_m4, ph7_m2, ph9_m4, ph9_m2, ab_2, pc_7 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

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

        pc_7[k] = e_0 * fs_15_16_14 * h3_m2 - e_1 * fs_5_4_6 * h5_m4 + e_1 * fs_35_8_2 * h5_m2 - e_1 * fs_5_8_14 * r_2 * h3_m2 - e_2 * fs_45_286_110 * h7_m4 + e_2 * fs_189_286_5 * h7_m2 + e_2 * fs_5_13_6 * r_2 * h5_m4 - e_2 * fs_35_26_2 * r_2 * h5_m2 + e_2 * fs_5_44_14 * r_4 * h3_m2 - e_3 * fs_1_221_858 * h9_m4 + e_3 * fs_1_221_231 * h9_m2 + e_3 * fs_45_2431_110 * r_2 * h7_m4 - e_3 * fs_189_2431_5 * r_2 * h7_m2 - e_3 * fs_1_39_6 * r_4 * h5_m4 + e_3 * fs_7_78_2 * r_4 * h5_m2 - e_3 * fs_5_858_14 * r_6 * h3_m2;
    }

    // NOTE: the angular components are formed in 49 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph3_m1, ph5_m5, ph5_m1, ph7_m7, ph7_m6, ph7_m5, ph7_m1, ph9_m7, ph9_m6, ph9_m5, ph9_m1, ab_2, pc_8, pc_9, pc_10 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

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

        pc_8[k] = - e_0 * fs_15_8_14 * h3_m1 - e_1 * fs_5_8_6 * h5_m5 - e_1 * fs_5_4_35 * h5_m1 + e_1 * fs_5_4_14 * r_2 * h3_m1 - e_2 * fs_135_286_11 * h7_m5 - e_2 * fs_105_143_3 * h7_m1 + e_2 * fs_5_26_6 * r_2 * h5_m5 + e_2 * fs_5_13_35 * r_2 * h5_m1 - e_2 * fs_5_22_14 * r_4 * h3_m1 - e_3 * fs_1_442_6006 * h9_m5 - e_3 * fs_1_221_105 * h9_m1 + e_3 * fs_135_2431_11 * r_2 * h7_m5 + e_3 * fs_210_2431_3 * r_2 * h7_m1 - e_3 * fs_1_78_6 * r_4 * h5_m5 - e_3 * fs_1_39_35 * r_4 * h5_m1 + e_3 * fs_5_429_14 * r_6 * h3_m1;

        pc_9[k] = - e_2 * fs_45_572_286 * h7_m6 - e_3 * fs_1_442_10010 * h9_m6 + e_3 * fs_45_4862_286 * r_2 * h7_m6;

        pc_10[k] = - e_0 * fs_15_8_105 * h3_m1 - e_1 * fs_5_4_42 * h5_m1 + e_1 * fs_5_4_105 * r_2 * h3_m1 - e_2 * fs_3_572_30030 * h7_m7 - e_2 * fs_135_572_10 * h7_m1 + e_2 * fs_5_13_42 * r_2 * h5_m1 - e_2 * fs_5_22_105 * r_4 * h3_m1 - e_3 * fs_2_221_1001 * h9_m7 - e_3 * fs_1_221_14 * h9_m1 + e_3 * fs_3_4862_30030 * r_2 * h7_m7 + e_3 * fs_135_4862_10 * r_2 * h7_m1 - e_3 * fs_1_39_42 * r_4 * h5_m1 + e_3 * fs_5_429_105 * r_6 * h3_m1;
    }

    // NOTE: the angular components are formed in 49 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph3_m3, ph3_m2, ph5_m3, ph5_m2, ph7_m3, ph7_m2, ph9_m9, ph9_m8, ph9_m3, ph9_m2, ab_2, pc_11, pc_12 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h3_m3 = ph3_m3[k];
        const auto h3_m2 = ph3_m2[k];
        const auto h5_m3 = ph5_m3[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h7_m3 = ph7_m3[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h9_m9 = ph9_m9[k];
        const auto h9_m8 = ph9_m8[k];
        const auto h9_m3 = ph9_m3[k];
        const auto h9_m2 = ph9_m2[k];

        pc_11[k] = - e_0 * fs_15_8_231 * h3_m2 - e_1 * fs_5_4_33 * h5_m2 + e_1 * fs_5_4_231 * r_2 * h3_m2 - e_2 * fs_15_572_330 * h7_m2 + e_2 * fs_5_13_33 * r_2 * h5_m2 - e_2 * fs_5_22_231 * r_4 * h3_m2 - e_3 * fs_2_221_1547 * h9_m8 - e_3 * fs_1_442_14 * h9_m2 + e_3 * fs_15_4862_330 * r_2 * h7_m2 - e_3 * fs_1_39_33 * r_4 * h5_m2 + e_3 * fs_5_429_231 * r_6 * h3_m2;

        pc_12[k] = - e_0 * fs_15_8_462 * h3_m3 - e_1 * fs_5_8_66 * h5_m3 + e_1 * fs_5_4_462 * r_2 * h3_m3 - e_2 * fs_9_286_55 * h7_m3 + e_2 * fs_5_26_66 * r_2 * h5_m3 - e_2 * fs_5_22_462 * r_4 * h3_m3 - e_3 * fs_1_221_9282 * h9_m9 - e_3 * fs_1_442_2 * h9_m3 + e_3 * fs_9_2431_55 * r_2 * h7_m3 - e_3 * fs_1_78_66 * r_4 * h5_m3 + e_3 * fs_5_429_462 * r_6 * h3_m3;
    }

    // NOTE: the angular components are formed in 49 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph3_p3, ph5_p3, ph5_p4, ph7_p3, ph7_p4, ph7_p7, ph9_p3, ph9_p4, ph9_p7, ph9_p8, ab_2, pc_13, pc_14 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h3_p3 = ph3_p3[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h5_p4 = ph5_p4[k];
        const auto h7_p3 = ph7_p3[k];
        const auto h7_p4 = ph7_p4[k];
        const auto h7_p7 = ph7_p7[k];
        const auto h9_p3 = ph9_p3[k];
        const auto h9_p4 = ph9_p4[k];
        const auto h9_p7 = ph9_p7[k];
        const auto h9_p8 = ph9_p8[k];

        pc_13[k] = - e_1 * fs_15_8_22 * h5_p4 - e_2 * fs_3_26_30 * h7_p4 + e_2 * fs_15_26_22 * r_2 * h5_p4 - e_3 * fs_1_442_26 * h9_p4 + e_3 * fs_1_221_3094 * h9_p8 + e_3 * fs_3_221_30 * r_2 * h7_p4 - e_3 * fs_1_26_22 * r_4 * h5_p4;

        pc_14[k] = e_0 * fs_15_8_231 * h3_p3 - e_1 * fs_5_4_33 * h5_p3 - e_1 * fs_5_4_231 * r_2 * h3_p3 - e_2 * fs_57_572_110 * h7_p3 - e_2 * fs_3_52_910 * h7_p7 + e_2 * fs_5_13_33 * r_2 * h5_p3 + e_2 * fs_5_22_231 * r_4 * h3_p3 - e_3 * f_6_221 * h9_p3 + e_3 * fs_4_221_273 * h9_p7 + e_3 * fs_57_4862_110 * r_2 * h7_p3 + e_3 * fs_3_442_910 * r_2 * h7_p7 - e_3 * fs_1_39_33 * r_4 * h5_p3 - e_3 * fs_5_429_231 * r_6 * h3_p3;
    }

    // NOTE: the angular components are formed in 49 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph3_p1, ph3_p2, ph5_p2, ph5_p5, ph7_p1, ph7_p2, ph7_p5, ph7_p6, ph9_p1, ph9_p2, ph9_p5, ph9_p6, ab_2, pc_15, pc_16 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h3_p1 = ph3_p1[k];
        const auto h3_p2 = ph3_p2[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h5_p5 = ph5_p5[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h7_p5 = ph7_p5[k];
        const auto h7_p6 = ph7_p6[k];
        const auto h9_p1 = ph9_p1[k];
        const auto h9_p2 = ph9_p2[k];
        const auto h9_p5 = ph9_p5[k];
        const auto h9_p6 = ph9_p6[k];

        pc_15[k] = e_0 * fs_45_4_7 * h3_p2 - e_1 * f_15_4 * h5_p2 - e_1 * fs_15_2_7 * r_2 * h3_p2 - e_2 * fs_60_143_10 * h7_p2 - e_2 * fs_6_143_1430 * h7_p6 + e_2 * f_15_13 * r_2 * h5_p2 + e_2 * fs_15_11_7 * r_4 * h3_p2 - e_3 * fs_1_442_462 * h9_p2 + e_3 * fs_3_442_2002 * h9_p6 + e_3 * fs_120_2431_10 * r_2 * h7_p2 + e_3 * fs_12_2431_1430 * r_2 * h7_p6 - e_3 * f_1_13 * r_4 * h5_p2 - e_3 * fs_10_143_7 * r_6 * h3_p2;

        pc_16[k] = e_0 * fs_45_8_21 * h3_p1 - e_1 * f_15_4 * h5_p5 - e_1 * fs_15_4_21 * r_2 * h3_p1 - e_2 * fs_45_44_2 * h7_p1 - e_2 * fs_75_572_66 * h7_p5 + e_2 * f_15_13 * r_2 * h5_p5 + e_2 * fs_15_22_21 * r_4 * h3_p1 - e_3 * fs_2_221_70 * h9_p1 + e_3 * fs_2_221_1001 * h9_p5 + e_3 * fs_45_374_2 * r_2 * h7_p1 + e_3 * fs_75_4862_66 * r_2 * h7_p5 - e_3 * f_1_13 * r_4 * h5_p5 - e_3 * fs_5_143_21 * r_6 * h3_p1;
    }

    // NOTE: the angular components are formed in 49 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph3_0, ph5_0, ph5_p4, ph7_0, ph7_p4, ph9_0, ph9_p4, ab_2, pc_17 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h3_0 = ph3_0[k];
        const auto h5_0 = ph5_0[k];
        const auto h5_p4 = ph5_p4[k];
        const auto h7_0 = ph7_0[k];
        const auto h7_p4 = ph7_p4[k];
        const auto h9_0 = ph9_0[k];
        const auto h9_p4 = ph9_p4[k];

        pc_17[k] = e_0 * fs_15_2_14 * h3_0 + e_1 * fs_5_4_14 * h5_0 - e_1 * fs_15_8_10 * h5_p4 - e_1 * fs_5_1_14 * r_2 * h3_0 - e_2 * fs_75_143_14 * h7_0 - e_2 * fs_15_286_66 * h7_p4 - e_2 * fs_5_13_14 * r_2 * h5_0 + e_2 * fs_15_26_10 * r_2 * h5_p4 + e_2 * fs_10_11_14 * r_4 * h3_0 - e_3 * fs_9_221_14 * h9_0 + e_3 * fs_3_442_1430 * h9_p4 + e_3 * fs_150_2431_14 * r_2 * h7_0 + e_3 * fs_15_2431_66 * r_2 * h7_p4 + e_3 * fs_1_39_14 * r_4 * h5_0 - e_3 * fs_1_26_10 * r_4 * h5_p4 - e_3 * fs_20_429_14 * r_6 * h3_0;
    }

    // NOTE: the angular components are formed in 49 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph3_m2, ph3_p1, ph3_p3, ph5_m2, ph5_p1, ph5_p3, ph7_m2, ph7_p1, ph7_p3, ph9_m2, ph9_p1, ph9_p3, ab_2, pc_18, pc_19 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h3_m2 = ph3_m2[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h3_p3 = ph3_p3[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h7_p3 = ph7_p3[k];
        const auto h9_m2 = ph9_m2[k];
        const auto h9_p1 = ph9_p1[k];
        const auto h9_p3 = ph9_p3[k];

        pc_18[k] = - e_0 * fs_15_16_210 * h3_p1 - e_0 * fs_15_16_14 * h3_p3 - e_1 * fs_5_4_21 * h5_p1 - e_1 * fs_5_1_2 * h5_p3 + e_1 * fs_5_8_210 * r_2 * h3_p1 + e_1 * fs_5_8_14 * r_2 * h3_p3 + e_2 * fs_147_286_5 * h7_p1 + e_2 * fs_15_286_15 * h7_p3 + e_2 * fs_5_13_21 * r_2 * h5_p1 + e_2 * fs_20_13_2 * r_2 * h5_p3 - e_2 * fs_5_44_210 * r_4 * h3_p1 - e_2 * fs_5_44_14 * r_4 * h3_p3 + e_3 * fs_12_221_7 * h9_p1 + e_3 * fs_6_221_66 * h9_p3 - e_3 * fs_147_2431_5 * r_2 * h7_p1 - e_3 * fs_15_2431_15 * r_2 * h7_p3 - e_3 * fs_1_39_21 * r_4 * h5_p1 - e_3 * fs_4_39_2 * r_4 * h5_p3 + e_3 * fs_5_858_210 * r_6 * h3_p1 + e_3 * fs_5_858_14 * r_6 * h3_p3;

        pc_19[k] = e_0 * f_45_4 * h3_m2 + e_1 * fs_15_4_7 * h5_m2 - e_1 * f_15_2 * r_2 * h3_m2 - e_2 * fs_18_143_70 * h7_m2 - e_2 * fs_15_13_7 * r_2 * h5_m2 + e_2 * f_15_11 * r_4 * h3_m2 - e_3 * fs_7_221_66 * h9_m2 + e_3 * fs_36_2431_70 * r_2 * h7_m2 + e_3 * fs_1_13_7 * r_4 * h5_m2 - e_3 * f_10_143 * r_6 * h3_m2;
    }

    // NOTE: the angular components are formed in 49 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph3_m3, ph3_m1, ph5_m4, ph5_m3, ph5_m1, ph7_m4, ph7_m3, ph7_m1, ph9_m4, ph9_m3, ph9_m1, ab_2, pc_20, pc_21 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

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

        pc_20[k] = e_0 * fs_15_16_14 * h3_m3 - e_0 * fs_15_16_210 * h3_m1 + e_1 * fs_5_1_2 * h5_m3 - e_1 * fs_5_4_21 * h5_m1 - e_1 * fs_5_8_14 * r_2 * h3_m3 + e_1 * fs_5_8_210 * r_2 * h3_m1 - e_2 * fs_15_286_15 * h7_m3 + e_2 * fs_147_286_5 * h7_m1 - e_2 * fs_20_13_2 * r_2 * h5_m3 + e_2 * fs_5_13_21 * r_2 * h5_m1 + e_2 * fs_5_44_14 * r_4 * h3_m3 - e_2 * fs_5_44_210 * r_4 * h3_m1 - e_3 * fs_6_221_66 * h9_m3 + e_3 * fs_12_221_7 * h9_m1 + e_3 * fs_15_2431_15 * r_2 * h7_m3 - e_3 * fs_147_2431_5 * r_2 * h7_m1 + e_3 * fs_4_39_2 * r_4 * h5_m3 - e_3 * fs_1_39_21 * r_4 * h5_m1 - e_3 * fs_5_858_14 * r_6 * h3_m3 + e_3 * fs_5_858_210 * r_6 * h3_m1;

        pc_21[k] = e_1 * fs_15_8_10 * h5_m4 + e_2 * fs_15_286_66 * h7_m4 - e_2 * fs_15_26_10 * r_2 * h5_m4 - e_3 * fs_3_442_1430 * h9_m4 - e_3 * fs_15_2431_66 * r_2 * h7_m4 + e_3 * fs_1_26_10 * r_4 * h5_m4;
    }

    // NOTE: the angular components are formed in 49 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph3_m2, ph3_m1, ph5_m5, ph5_m2, ph7_m6, ph7_m5, ph7_m2, ph7_m1, ph9_m6, ph9_m5, ph9_m2, ph9_m1, ab_2, pc_22, pc_23 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h3_m2 = ph3_m2[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h5_m5 = ph5_m5[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h7_m6 = ph7_m6[k];
        const auto h7_m5 = ph7_m5[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h7_m1 = ph7_m1[k];
        const auto h9_m6 = ph9_m6[k];
        const auto h9_m5 = ph9_m5[k];
        const auto h9_m2 = ph9_m2[k];
        const auto h9_m1 = ph9_m1[k];

        pc_22[k] = - e_0 * fs_45_8_21 * h3_m1 + e_1 * f_15_4 * h5_m5 + e_1 * fs_15_4_21 * r_2 * h3_m1 + e_2 * fs_75_572_66 * h7_m5 + e_2 * fs_45_44_2 * h7_m1 - e_2 * f_15_13 * r_2 * h5_m5 - e_2 * fs_15_22_21 * r_4 * h3_m1 - e_3 * fs_2_221_1001 * h9_m5 + e_3 * fs_2_221_70 * h9_m1 - e_3 * fs_75_4862_66 * r_2 * h7_m5 - e_3 * fs_45_374_2 * r_2 * h7_m1 + e_3 * f_1_13 * r_4 * h5_m5 + e_3 * fs_5_143_21 * r_6 * h3_m1;

        pc_23[k] = - e_0 * fs_45_4_7 * h3_m2 + e_1 * f_15_4 * h5_m2 + e_1 * fs_15_2_7 * r_2 * h3_m2 + e_2 * fs_6_143_1430 * h7_m6 + e_2 * fs_60_143_10 * h7_m2 - e_2 * f_15_13 * r_2 * h5_m2 - e_2 * fs_15_11_7 * r_4 * h3_m2 - e_3 * fs_3_442_2002 * h9_m6 + e_3 * fs_1_442_462 * h9_m2 - e_3 * fs_12_2431_1430 * r_2 * h7_m6 - e_3 * fs_120_2431_10 * r_2 * h7_m2 + e_3 * f_1_13 * r_4 * h5_m2 + e_3 * fs_10_143_7 * r_6 * h3_m2;
    }

    // NOTE: the angular components are formed in 49 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph3_m3, ph5_m4, ph5_m3, ph7_m7, ph7_m4, ph7_m3, ph9_m8, ph9_m7, ph9_m4, ph9_m3, ab_2, pc_24, pc_25 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h3_m3 = ph3_m3[k];
        const auto h5_m4 = ph5_m4[k];
        const auto h5_m3 = ph5_m3[k];
        const auto h7_m7 = ph7_m7[k];
        const auto h7_m4 = ph7_m4[k];
        const auto h7_m3 = ph7_m3[k];
        const auto h9_m8 = ph9_m8[k];
        const auto h9_m7 = ph9_m7[k];
        const auto h9_m4 = ph9_m4[k];
        const auto h9_m3 = ph9_m3[k];

        pc_24[k] = - e_0 * fs_15_8_231 * h3_m3 + e_1 * fs_5_4_33 * h5_m3 + e_1 * fs_5_4_231 * r_2 * h3_m3 + e_2 * fs_3_52_910 * h7_m7 + e_2 * fs_57_572_110 * h7_m3 - e_2 * fs_5_13_33 * r_2 * h5_m3 - e_2 * fs_5_22_231 * r_4 * h3_m3 - e_3 * fs_4_221_273 * h9_m7 + e_3 * f_6_221 * h9_m3 - e_3 * fs_3_442_910 * r_2 * h7_m7 - e_3 * fs_57_4862_110 * r_2 * h7_m3 + e_3 * fs_1_39_33 * r_4 * h5_m3 + e_3 * fs_5_429_231 * r_6 * h3_m3;

        pc_25[k] = e_1 * fs_15_8_22 * h5_m4 + e_2 * fs_3_26_30 * h7_m4 - e_2 * fs_15_26_22 * r_2 * h5_m4 - e_3 * fs_1_221_3094 * h9_m8 + e_3 * fs_1_442_26 * h9_m4 - e_3 * fs_3_221_30 * r_2 * h7_m4 + e_3 * fs_1_26_22 * r_4 * h5_m4;
    }

    // NOTE: the angular components are formed in 49 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, ph5_p5, ph7_p4, ph7_p5, ph7_p6, ph7_p7, ph9_p4, ph9_p5, ph9_p6, ph9_p7, ab_2, pc_26, pc_27 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

        const auto h5_p5 = ph5_p5[k];
        const auto h7_p4 = ph7_p4[k];
        const auto h7_p5 = ph7_p5[k];
        const auto h7_p6 = ph7_p6[k];
        const auto h7_p7 = ph7_p7[k];
        const auto h9_p4 = ph9_p4[k];
        const auto h9_p5 = ph9_p5[k];
        const auto h9_p6 = ph9_p6[k];
        const auto h9_p7 = ph9_p7[k];

        pc_26[k] = e_1 * fs_15_8_22 * h5_p5 + e_2 * fs_9_13_3 * h7_p5 + e_2 * fs_3_26_273 * h7_p7 - e_2 * fs_15_26_22 * r_2 * h5_p5 + e_3 * fs_1_442_182 * h9_p5 + e_3 * fs_1_221_910 * h9_p7 - e_3 * fs_18_221_3 * r_2 * h7_p5 - e_3 * fs_3_221_273 * r_2 * h7_p7 + e_3 * fs_1_26_22 * r_4 * h5_p5;

        pc_27[k] = e_2 * f_3_2 * h7_p4 + e_2 * fs_3_52_26 * h7_p6 + e_3 * fs_1_221_195 * h9_p4 + e_3 * fs_3_442_910 * h9_p6 - e_3 * f_3_17 * r_2 * h7_p4 - e_3 * fs_3_442_26 * r_2 * h7_p6;
    }

    // NOTE: the angular components are formed in 49 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph3_p3, ph5_p3, ph5_p5, ph7_p3, ph7_p5, ph9_p3, ph9_p5, ab_2, pc_28 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h3_p3 = ph3_p3[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h5_p5 = ph5_p5[k];
        const auto h7_p3 = ph7_p3[k];
        const auto h7_p5 = ph7_p5[k];
        const auto h9_p3 = ph9_p3[k];
        const auto h9_p5 = ph9_p5[k];

        pc_28[k] = e_0 * fs_15_8_105 * h3_p3 - e_1 * fs_5_4_15 * h5_p3 + e_1 * fs_15_4_3 * h5_p5 - e_1 * fs_5_4_105 * r_2 * h3_p3 + e_2 * fs_573_572_2 * h7_p3 - e_2 * fs_87_572_22 * h7_p5 + e_2 * fs_5_13_15 * r_2 * h5_p3 - e_2 * fs_15_13_3 * r_2 * h5_p5 + e_2 * fs_5_22_105 * r_4 * h3_p3 + e_3 * fs_3_221_55 * h9_p3 + e_3 * fs_1_221_3003 * h9_p5 - e_3 * fs_573_4862_2 * r_2 * h7_p3 + e_3 * fs_87_4862_22 * r_2 * h7_p5 - e_3 * fs_1_39_15 * r_4 * h5_p3 + e_3 * fs_1_13_3 * r_4 * h5_p5 - e_3 * fs_5_429_105 * r_6 * h3_p3;
    }

    // NOTE: the angular components are formed in 49 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph3_p2, ph5_p2, ph5_p4, ph7_p2, ph7_p4, ph9_p2, ph9_p4, ab_2, pc_29 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

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

        pc_29[k] = e_0 * fs_45_8_21 * h3_p2 - e_1 * fs_15_4_3 * h5_p2 + e_1 * f_15_2 * h5_p4 - e_1 * fs_15_4_21 * r_2 * h3_p2 + e_2 * fs_111_572_30 * h7_p2 - e_2 * fs_27_286_165 * h7_p4 + e_2 * fs_15_13_3 * r_2 * h5_p2 - e_2 * f_30_13 * r_2 * h5_p4 + e_2 * fs_15_22_21 * r_4 * h3_p2 + e_3 * fs_5_442_154 * h9_p2 + e_3 * fs_5_221_143 * h9_p4 - e_3 * fs_111_4862_30 * r_2 * h7_p2 + e_3 * fs_27_2431_165 * r_2 * h7_p4 - e_3 * fs_1_13_3 * r_4 * h5_p2 + e_3 * f_2_13 * r_4 * h5_p4 - e_3 * fs_5_143_21 * r_6 * h3_p2;
    }

    // NOTE: the angular components are formed in 49 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph3_p1, ph3_p3, ph5_p1, ph5_p3, ph7_p1, ph7_p3, ph9_p1, ph9_p3, ab_2, pc_30 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h3_p1 = ph3_p1[k];
        const auto h3_p3 = ph3_p3[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h7_p3 = ph7_p3[k];
        const auto h9_p1 = ph9_p1[k];
        const auto h9_p3 = ph9_p3[k];

        pc_30[k] = e_0 * fs_15_8_210 * h3_p1 + e_0 * fs_15_8_14 * h3_p3 - e_1 * fs_5_4_21 * h5_p1 + e_1 * fs_35_8_2 * h5_p3 - e_1 * fs_5_4_210 * r_2 * h3_p1 - e_1 * fs_5_4_14 * r_2 * h3_p3 + e_2 * fs_69_286_5 * h7_p1 - e_2 * fs_48_143_15 * h7_p3 + e_2 * fs_5_13_21 * r_2 * h5_p1 - e_2 * fs_35_26_2 * r_2 * h5_p3 + e_2 * fs_5_22_210 * r_4 * h3_p1 + e_2 * fs_5_22_14 * r_4 * h3_p3 + e_3 * fs_15_221_7 * h9_p1 + e_3 * fs_15_442_66 * h9_p3 - e_3 * fs_69_2431_5 * r_2 * h7_p1 + e_3 * fs_96_2431_15 * r_2 * h7_p3 - e_3 * fs_1_39_21 * r_4 * h5_p1 + e_3 * fs_7_78_2 * r_4 * h5_p3 - e_3 * fs_5_429_210 * r_6 * h3_p1 - e_3 * fs_5_429_14 * r_6 * h3_p3;
    }

    // NOTE: the angular components are formed in 49 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph3_m1, ph3_0, ph3_p2, ph5_0, ph5_p2, ph7_m1, ph7_0, ph7_p2, ph9_m1, ph9_0, ph9_p2, ab_2, pc_31, pc_32 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h3_m1 = ph3_m1[k];
        const auto h3_0 = ph3_0[k];
        const auto h3_p2 = ph3_p2[k];
        const auto h5_0 = ph5_0[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h7_m1 = ph7_m1[k];
        const auto h7_0 = ph7_0[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h9_m1 = ph9_m1[k];
        const auto h9_0 = ph9_0[k];
        const auto h9_p2 = ph9_p2[k];

        pc_31[k] = e_0 * fs_75_8_14 * h3_0 + e_0 * fs_15_16_210 * h3_p2 - e_1 * fs_5_4_14 * h5_0 + e_1 * fs_5_8_30 * h5_p2 - e_1 * fs_25_4_14 * r_2 * h3_0 - e_1 * fs_5_8_210 * r_2 * h3_p2 - e_2 * fs_3_143_14 * h7_0 - e_2 * fs_177_286_3 * h7_p2 + e_2 * fs_5_13_14 * r_2 * h5_0 - e_2 * fs_5_26_30 * r_2 * h5_p2 + e_2 * fs_25_22_14 * r_4 * h3_0 + e_2 * fs_5_44_210 * r_4 * h3_p2 + e_3 * fs_18_221_14 * h9_0 + e_3 * fs_3_221_385 * h9_p2 + e_3 * fs_6_2431_14 * r_2 * h7_0 + e_3 * fs_177_2431_3 * r_2 * h7_p2 - e_3 * fs_1_39_14 * r_4 * h5_0 + e_3 * fs_1_78_30 * r_4 * h5_p2 - e_3 * fs_25_429_14 * r_6 * h3_0 - e_3 * fs_5_858_210 * r_6 * h3_p2;

        pc_32[k] = - e_0 * f_225_8 * h3_m1 + e_1 * f_75_4 * r_2 * h3_m1 + e_2 * fs_3_22_42 * h7_m1 - e_2 * f_75_22 * r_4 * h3_m1 - e_3 * fs_14_221_30 * h9_m1 - e_3 * fs_3_187_42 * r_2 * h7_m1 + e_3 * f_25_143 * r_6 * h3_m1;
    }

    // NOTE: the angular components are formed in 49 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph3_m3, ph3_m2, ph3_m1, ph5_m3, ph5_m2, ph5_m1, ph7_m3, ph7_m2, ph7_m1, ph9_m3, ph9_m2, ph9_m1, ab_2, pc_33, pc_34 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h3_m3 = ph3_m3[k];
        const auto h3_m2 = ph3_m2[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h5_m3 = ph5_m3[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h7_m3 = ph7_m3[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h7_m1 = ph7_m1[k];
        const auto h9_m3 = ph9_m3[k];
        const auto h9_m2 = ph9_m2[k];
        const auto h9_m1 = ph9_m1[k];

        pc_33[k] = - e_0 * fs_15_16_210 * h3_m2 - e_1 * fs_5_8_30 * h5_m2 + e_1 * fs_5_8_210 * r_2 * h3_m2 + e_2 * fs_177_286_3 * h7_m2 + e_2 * fs_5_26_30 * r_2 * h5_m2 - e_2 * fs_5_44_210 * r_4 * h3_m2 - e_3 * fs_3_221_385 * h9_m2 - e_3 * fs_177_2431_3 * r_2 * h7_m2 - e_3 * fs_1_78_30 * r_4 * h5_m2 + e_3 * fs_5_858_210 * r_6 * h3_m2;

        pc_34[k] = - e_0 * fs_15_8_14 * h3_m3 - e_0 * fs_15_8_210 * h3_m1 - e_1 * fs_35_8_2 * h5_m3 + e_1 * fs_5_4_21 * h5_m1 + e_1 * fs_5_4_14 * r_2 * h3_m3 + e_1 * fs_5_4_210 * r_2 * h3_m1 + e_2 * fs_48_143_15 * h7_m3 - e_2 * fs_69_286_5 * h7_m1 + e_2 * fs_35_26_2 * r_2 * h5_m3 - e_2 * fs_5_13_21 * r_2 * h5_m1 - e_2 * fs_5_22_14 * r_4 * h3_m3 - e_2 * fs_5_22_210 * r_4 * h3_m1 - e_3 * fs_15_442_66 * h9_m3 - e_3 * fs_15_221_7 * h9_m1 - e_3 * fs_96_2431_15 * r_2 * h7_m3 + e_3 * fs_69_2431_5 * r_2 * h7_m1 - e_3 * fs_7_78_2 * r_4 * h5_m3 + e_3 * fs_1_39_21 * r_4 * h5_m1 + e_3 * fs_5_429_14 * r_6 * h3_m3 + e_3 * fs_5_429_210 * r_6 * h3_m1;
    }

    // NOTE: the angular components are formed in 49 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph3_m2, ph5_m4, ph5_m2, ph7_m4, ph7_m2, ph9_m4, ph9_m2, ab_2, pc_35 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

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

        pc_35[k] = - e_0 * fs_45_8_21 * h3_m2 - e_1 * f_15_2 * h5_m4 + e_1 * fs_15_4_3 * h5_m2 + e_1 * fs_15_4_21 * r_2 * h3_m2 + e_2 * fs_27_286_165 * h7_m4 - e_2 * fs_111_572_30 * h7_m2 + e_2 * f_30_13 * r_2 * h5_m4 - e_2 * fs_15_13_3 * r_2 * h5_m2 - e_2 * fs_15_22_21 * r_4 * h3_m2 - e_3 * fs_5_221_143 * h9_m4 - e_3 * fs_5_442_154 * h9_m2 - e_3 * fs_27_2431_165 * r_2 * h7_m4 + e_3 * fs_111_4862_30 * r_2 * h7_m2 - e_3 * f_2_13 * r_4 * h5_m4 + e_3 * fs_1_13_3 * r_4 * h5_m2 + e_3 * fs_5_143_21 * r_6 * h3_m2;
    }

    // NOTE: the angular components are formed in 49 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph3_m3, ph5_m5, ph5_m3, ph7_m6, ph7_m5, ph7_m4, ph7_m3, ph9_m6, ph9_m5, ph9_m4, ph9_m3, ab_2, pc_36, pc_37 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h3_m3 = ph3_m3[k];
        const auto h5_m5 = ph5_m5[k];
        const auto h5_m3 = ph5_m3[k];
        const auto h7_m6 = ph7_m6[k];
        const auto h7_m5 = ph7_m5[k];
        const auto h7_m4 = ph7_m4[k];
        const auto h7_m3 = ph7_m3[k];
        const auto h9_m6 = ph9_m6[k];
        const auto h9_m5 = ph9_m5[k];
        const auto h9_m4 = ph9_m4[k];
        const auto h9_m3 = ph9_m3[k];

        pc_36[k] = - e_0 * fs_15_8_105 * h3_m3 - e_1 * fs_15_4_3 * h5_m5 + e_1 * fs_5_4_15 * h5_m3 + e_1 * fs_5_4_105 * r_2 * h3_m3 + e_2 * fs_87_572_22 * h7_m5 - e_2 * fs_573_572_2 * h7_m3 + e_2 * fs_15_13_3 * r_2 * h5_m5 - e_2 * fs_5_13_15 * r_2 * h5_m3 - e_2 * fs_5_22_105 * r_4 * h3_m3 - e_3 * fs_1_221_3003 * h9_m5 - e_3 * fs_3_221_55 * h9_m3 - e_3 * fs_87_4862_22 * r_2 * h7_m5 + e_3 * fs_573_4862_2 * r_2 * h7_m3 - e_3 * fs_1_13_3 * r_4 * h5_m5 + e_3 * fs_1_39_15 * r_4 * h5_m3 + e_3 * fs_5_429_105 * r_6 * h3_m3;

        pc_37[k] = - e_2 * fs_3_52_26 * h7_m6 - e_2 * f_3_2 * h7_m4 - e_3 * fs_3_442_910 * h9_m6 - e_3 * fs_1_221_195 * h9_m4 + e_3 * fs_3_442_26 * r_2 * h7_m6 + e_3 * f_3_17 * r_2 * h7_m4;
    }

    // NOTE: the angular components are formed in 49 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, ph5_m5, ph5_m4, ph7_m7, ph7_m6, ph7_m5, ph7_m4, ph9_m7, ph9_m6, ph9_m5, ph9_m4, ab_2, pc_38, pc_39, pc_40, pc_41 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

        const auto h5_m5 = ph5_m5[k];
        const auto h5_m4 = ph5_m4[k];
        const auto h7_m7 = ph7_m7[k];
        const auto h7_m6 = ph7_m6[k];
        const auto h7_m5 = ph7_m5[k];
        const auto h7_m4 = ph7_m4[k];
        const auto h9_m7 = ph9_m7[k];
        const auto h9_m6 = ph9_m6[k];
        const auto h9_m5 = ph9_m5[k];
        const auto h9_m4 = ph9_m4[k];

        pc_38[k] = - e_1 * fs_15_8_22 * h5_m5 - e_2 * fs_3_26_273 * h7_m7 - e_2 * fs_9_13_3 * h7_m5 + e_2 * fs_15_26_22 * r_2 * h5_m5 - e_3 * fs_1_221_910 * h9_m7 - e_3 * fs_1_442_182 * h9_m5 + e_3 * fs_3_221_273 * r_2 * h7_m7 + e_3 * fs_18_221_3 * r_2 * h7_m5 - e_3 * fs_1_26_22 * r_4 * h5_m5;

        pc_39[k] = - e_2 * fs_9_13_13 * h7_m6 - e_3 * fs_1_221_455 * h9_m6 + e_3 * fs_18_221_13 * r_2 * h7_m6;

        pc_40[k] = e_1 * fs_15_4_11 * h5_m5 - e_2 * fs_21_26_6 * h7_m5 - e_2 * fs_15_13_11 * r_2 * h5_m5 - e_3 * fs_4_221_91 * h9_m5 + e_3 * fs_21_221_6 * r_2 * h7_m5 + e_3 * fs_1_13_11 * r_4 * h5_m5;

        pc_41[k] = e_1 * fs_15_4_5 * h5_m4 - e_2 * fs_24_143_33 * h7_m4 - e_2 * fs_15_13_5 * r_2 * h5_m4 - e_3 * fs_2_221_715 * h9_m4 + e_3 * fs_48_2431_33 * r_2 * h7_m4 + e_3 * fs_1_13_5 * r_4 * h5_m4;
    }

    // NOTE: the angular components are formed in 49 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph3_m3, ph3_m2, ph3_m1, ph5_m3, ph5_m2, ph5_m1, ph7_m3, ph7_m2, ph7_m1, ph9_m3, ph9_m2, ph9_m1, ab_2, pc_42, pc_43, pc_44 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h3_m3 = ph3_m3[k];
        const auto h3_m2 = ph3_m2[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h5_m3 = ph5_m3[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h7_m3 = ph7_m3[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h7_m1 = ph7_m1[k];
        const auto h9_m3 = ph9_m3[k];
        const auto h9_m2 = ph9_m2[k];
        const auto h9_m1 = ph9_m1[k];

        pc_42[k] = e_0 * fs_15_4_21 * h3_m3 + e_1 * fs_5_4_3 * h5_m3 - e_1 * fs_5_2_21 * r_2 * h3_m3 + e_2 * fs_9_286_10 * h7_m3 - e_2 * fs_5_13_3 * r_2 * h5_m3 + e_2 * fs_5_11_21 * r_4 * h3_m3 - e_3 * fs_20_221_11 * h9_m3 - e_3 * fs_9_2431_10 * r_2 * h7_m3 + e_3 * fs_1_39_3 * r_4 * h5_m3 - e_3 * fs_10_429_21 * r_6 * h3_m3;

        pc_43[k] = e_0 * fs_15_2_14 * h3_m2 - e_1 * fs_5_2_2 * h5_m2 - e_1 * fs_5_1_14 * r_2 * h3_m2 + e_2 * fs_63_143_5 * h7_m2 + e_2 * fs_10_13_2 * r_2 * h5_m2 + e_2 * fs_10_11_14 * r_4 * h3_m2 - e_3 * fs_5_221_231 * h9_m2 - e_3 * fs_126_2431_5 * r_2 * h7_m2 - e_3 * fs_2_39_2 * r_4 * h5_m2 - e_3 * fs_20_429_14 * r_6 * h3_m2;

        pc_44[k] = e_0 * fs_75_8_14 * h3_m1 - e_1 * fs_5_4_35 * h5_m1 - e_1 * fs_25_4_14 * r_2 * h3_m1 + e_2 * fs_129_143_3 * h7_m1 + e_2 * fs_5_13_35 * r_2 * h5_m1 + e_2 * fs_25_22_14 * r_4 * h3_m1 - e_3 * fs_8_221_105 * h9_m1 - e_3 * fs_258_2431_3 * r_2 * h7_m1 - e_3 * fs_1_39_35 * r_4 * h5_m1 - e_3 * fs_25_429_14 * r_6 * h3_m1;
    }

    // NOTE: the angular components are formed in 49 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph3_0, ph3_p1, ph3_p2, ph5_0, ph5_p1, ph5_p2, ph7_0, ph7_p1, ph7_p2, ph9_0, ph9_p1, ph9_p2, ab_2, pc_45, pc_46, pc_47 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h3_0 = ph3_0[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h3_p2 = ph3_p2[k];
        const auto h5_0 = ph5_0[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h7_0 = ph7_0[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h9_0 = ph9_0[k];
        const auto h9_p1 = ph9_p1[k];
        const auto h9_p2 = ph9_p2[k];

        pc_45[k] = e_0 * f_75_2 * h3_0 - e_1 * f_35_4 * h5_0 - e_1 * f_25_1 * r_2 * h3_0 + e_2 * f_252_143 * h7_0 + e_2 * f_35_13 * r_2 * h5_0 + e_2 * f_50_11 * r_4 * h3_0 - e_3 * f_84_221 * h9_0 - e_3 * f_504_2431 * r_2 * h7_0 - e_3 * f_7_39 * r_4 * h5_0 - e_3 * f_100_429 * r_6 * h3_0;

        pc_46[k] = e_0 * fs_75_8_14 * h3_p1 - e_1 * fs_5_4_35 * h5_p1 - e_1 * fs_25_4_14 * r_2 * h3_p1 + e_2 * fs_129_143_3 * h7_p1 + e_2 * fs_5_13_35 * r_2 * h5_p1 + e_2 * fs_25_22_14 * r_4 * h3_p1 - e_3 * fs_8_221_105 * h9_p1 - e_3 * fs_258_2431_3 * r_2 * h7_p1 - e_3 * fs_1_39_35 * r_4 * h5_p1 - e_3 * fs_25_429_14 * r_6 * h3_p1;

        pc_47[k] = e_0 * fs_15_2_14 * h3_p2 - e_1 * fs_5_2_2 * h5_p2 - e_1 * fs_5_1_14 * r_2 * h3_p2 + e_2 * fs_63_143_5 * h7_p2 + e_2 * fs_10_13_2 * r_2 * h5_p2 + e_2 * fs_10_11_14 * r_4 * h3_p2 - e_3 * fs_5_221_231 * h9_p2 - e_3 * fs_126_2431_5 * r_2 * h7_p2 - e_3 * fs_2_39_2 * r_4 * h5_p2 - e_3 * fs_20_429_14 * r_6 * h3_p2;
    }

    // NOTE: the angular components are formed in 49 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph3_p3, ph5_p3, ph5_p4, ph5_p5, ph7_p3, ph7_p4, ph7_p5, ph9_p3, ph9_p4, ph9_p5, ab_2, pc_48, pc_49, pc_50 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

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

        pc_48[k] = e_0 * fs_15_4_21 * h3_p3 + e_1 * fs_5_4_3 * h5_p3 - e_1 * fs_5_2_21 * r_2 * h3_p3 + e_2 * fs_9_286_10 * h7_p3 - e_2 * fs_5_13_3 * r_2 * h5_p3 + e_2 * fs_5_11_21 * r_4 * h3_p3 - e_3 * fs_20_221_11 * h9_p3 - e_3 * fs_9_2431_10 * r_2 * h7_p3 + e_3 * fs_1_39_3 * r_4 * h5_p3 - e_3 * fs_10_429_21 * r_6 * h3_p3;

        pc_49[k] = e_1 * fs_15_4_5 * h5_p4 - e_2 * fs_24_143_33 * h7_p4 - e_2 * fs_15_13_5 * r_2 * h5_p4 - e_3 * fs_2_221_715 * h9_p4 + e_3 * fs_48_2431_33 * r_2 * h7_p4 + e_3 * fs_1_13_5 * r_4 * h5_p4;

        pc_50[k] = e_1 * fs_15_4_11 * h5_p5 - e_2 * fs_21_26_6 * h7_p5 - e_2 * fs_15_13_11 * r_2 * h5_p5 - e_3 * fs_4_221_91 * h9_p5 + e_3 * fs_21_221_6 * r_2 * h7_p5 + e_3 * fs_1_13_11 * r_4 * h5_p5;
    }

    // NOTE: the angular components are formed in 49 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, ph5_m5, ph7_m7, ph7_m6, ph7_m5, ph7_m4, ph7_p6, ph9_m7, ph9_m6, ph9_m5, ph9_m4, ph9_p6, ab_2, pc_51, pc_52, pc_53 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

        const auto h5_m5 = ph5_m5[k];
        const auto h7_m7 = ph7_m7[k];
        const auto h7_m6 = ph7_m6[k];
        const auto h7_m5 = ph7_m5[k];
        const auto h7_m4 = ph7_m4[k];
        const auto h7_p6 = ph7_p6[k];
        const auto h9_m7 = ph9_m7[k];
        const auto h9_m6 = ph9_m6[k];
        const auto h9_m5 = ph9_m5[k];
        const auto h9_m4 = ph9_m4[k];
        const auto h9_p6 = ph9_p6[k];

        pc_51[k] = - e_2 * fs_9_13_13 * h7_p6 - e_3 * fs_1_221_455 * h9_p6 + e_3 * fs_18_221_13 * r_2 * h7_p6;

        pc_52[k] = e_1 * fs_15_8_22 * h5_m5 - e_2 * fs_3_26_273 * h7_m7 + e_2 * fs_9_13_3 * h7_m5 - e_2 * fs_15_26_22 * r_2 * h5_m5 - e_3 * fs_1_221_910 * h9_m7 + e_3 * fs_1_442_182 * h9_m5 + e_3 * fs_3_221_273 * r_2 * h7_m7 - e_3 * fs_18_221_3 * r_2 * h7_m5 + e_3 * fs_1_26_22 * r_4 * h5_m5;

        pc_53[k] = - e_2 * fs_3_52_26 * h7_m6 + e_2 * f_3_2 * h7_m4 - e_3 * fs_3_442_910 * h9_m6 + e_3 * fs_1_221_195 * h9_m4 + e_3 * fs_3_442_26 * r_2 * h7_m6 - e_3 * f_3_17 * r_2 * h7_m4;
    }

    // NOTE: the angular components are formed in 49 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph3_m3, ph5_m5, ph5_m3, ph7_m5, ph7_m3, ph9_m5, ph9_m3, ab_2, pc_54 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

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

        pc_54[k] = e_0 * fs_15_8_105 * h3_m3 - e_1 * fs_15_4_3 * h5_m5 - e_1 * fs_5_4_15 * h5_m3 - e_1 * fs_5_4_105 * r_2 * h3_m3 + e_2 * fs_87_572_22 * h7_m5 + e_2 * fs_573_572_2 * h7_m3 + e_2 * fs_15_13_3 * r_2 * h5_m5 + e_2 * fs_5_13_15 * r_2 * h5_m3 + e_2 * fs_5_22_105 * r_4 * h3_m3 - e_3 * fs_1_221_3003 * h9_m5 + e_3 * fs_3_221_55 * h9_m3 - e_3 * fs_87_4862_22 * r_2 * h7_m5 - e_3 * fs_573_4862_2 * r_2 * h7_m3 - e_3 * fs_1_13_3 * r_4 * h5_m5 - e_3 * fs_1_39_15 * r_4 * h5_m3 - e_3 * fs_5_429_105 * r_6 * h3_m3;
    }

    // NOTE: the angular components are formed in 49 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph3_m2, ph5_m4, ph5_m2, ph7_m4, ph7_m2, ph9_m4, ph9_m2, ab_2, pc_55 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

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

        pc_55[k] = e_0 * fs_45_8_21 * h3_m2 - e_1 * f_15_2 * h5_m4 - e_1 * fs_15_4_3 * h5_m2 - e_1 * fs_15_4_21 * r_2 * h3_m2 + e_2 * fs_27_286_165 * h7_m4 + e_2 * fs_111_572_30 * h7_m2 + e_2 * f_30_13 * r_2 * h5_m4 + e_2 * fs_15_13_3 * r_2 * h5_m2 + e_2 * fs_15_22_21 * r_4 * h3_m2 - e_3 * fs_5_221_143 * h9_m4 + e_3 * fs_5_442_154 * h9_m2 - e_3 * fs_27_2431_165 * r_2 * h7_m4 - e_3 * fs_111_4862_30 * r_2 * h7_m2 - e_3 * f_2_13 * r_4 * h5_m4 - e_3 * fs_1_13_3 * r_4 * h5_m2 - e_3 * fs_5_143_21 * r_6 * h3_m2;
    }

    // NOTE: the angular components are formed in 49 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph3_m3, ph3_m2, ph3_m1, ph5_m3, ph5_m2, ph5_m1, ph7_m3, ph7_m2, ph7_m1, ph9_m3, ph9_m2, ph9_m1, ab_2, pc_56, pc_57 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h3_m3 = ph3_m3[k];
        const auto h3_m2 = ph3_m2[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h5_m3 = ph5_m3[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h7_m3 = ph7_m3[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h7_m1 = ph7_m1[k];
        const auto h9_m3 = ph9_m3[k];
        const auto h9_m2 = ph9_m2[k];
        const auto h9_m1 = ph9_m1[k];

        pc_56[k] = - e_0 * fs_15_8_14 * h3_m3 + e_0 * fs_15_8_210 * h3_m1 - e_1 * fs_35_8_2 * h5_m3 - e_1 * fs_5_4_21 * h5_m1 + e_1 * fs_5_4_14 * r_2 * h3_m3 - e_1 * fs_5_4_210 * r_2 * h3_m1 + e_2 * fs_48_143_15 * h7_m3 + e_2 * fs_69_286_5 * h7_m1 + e_2 * fs_35_26_2 * r_2 * h5_m3 + e_2 * fs_5_13_21 * r_2 * h5_m1 - e_2 * fs_5_22_14 * r_4 * h3_m3 + e_2 * fs_5_22_210 * r_4 * h3_m1 - e_3 * fs_15_442_66 * h9_m3 + e_3 * fs_15_221_7 * h9_m1 - e_3 * fs_96_2431_15 * r_2 * h7_m3 - e_3 * fs_69_2431_5 * r_2 * h7_m1 - e_3 * fs_7_78_2 * r_4 * h5_m3 - e_3 * fs_1_39_21 * r_4 * h5_m1 + e_3 * fs_5_429_14 * r_6 * h3_m3 - e_3 * fs_5_429_210 * r_6 * h3_m1;

        pc_57[k] = - e_0 * fs_15_16_210 * h3_m2 - e_1 * fs_5_8_30 * h5_m2 + e_1 * fs_5_8_210 * r_2 * h3_m2 + e_2 * fs_177_286_3 * h7_m2 + e_2 * fs_5_26_30 * r_2 * h5_m2 - e_2 * fs_5_44_210 * r_4 * h3_m2 - e_3 * fs_3_221_385 * h9_m2 - e_3 * fs_177_2431_3 * r_2 * h7_m2 - e_3 * fs_1_78_30 * r_4 * h5_m2 + e_3 * fs_5_858_210 * r_6 * h3_m2;
    }

    // NOTE: the angular components are formed in 49 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph3_0, ph3_p1, ph3_p2, ph5_0, ph5_p2, ph7_0, ph7_p1, ph7_p2, ph9_0, ph9_p1, ph9_p2, ab_2, pc_58, pc_59 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h3_0 = ph3_0[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h3_p2 = ph3_p2[k];
        const auto h5_0 = ph5_0[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h7_0 = ph7_0[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h9_0 = ph9_0[k];
        const auto h9_p1 = ph9_p1[k];
        const auto h9_p2 = ph9_p2[k];

        pc_58[k] = - e_0 * f_225_8 * h3_p1 + e_1 * f_75_4 * r_2 * h3_p1 + e_2 * fs_3_22_42 * h7_p1 - e_2 * f_75_22 * r_4 * h3_p1 - e_3 * fs_14_221_30 * h9_p1 - e_3 * fs_3_187_42 * r_2 * h7_p1 + e_3 * f_25_143 * r_6 * h3_p1;

        pc_59[k] = e_0 * fs_75_8_14 * h3_0 - e_0 * fs_15_16_210 * h3_p2 - e_1 * fs_5_4_14 * h5_0 - e_1 * fs_5_8_30 * h5_p2 - e_1 * fs_25_4_14 * r_2 * h3_0 + e_1 * fs_5_8_210 * r_2 * h3_p2 - e_2 * fs_3_143_14 * h7_0 + e_2 * fs_177_286_3 * h7_p2 + e_2 * fs_5_13_14 * r_2 * h5_0 + e_2 * fs_5_26_30 * r_2 * h5_p2 + e_2 * fs_25_22_14 * r_4 * h3_0 - e_2 * fs_5_44_210 * r_4 * h3_p2 + e_3 * fs_18_221_14 * h9_0 - e_3 * fs_3_221_385 * h9_p2 + e_3 * fs_6_2431_14 * r_2 * h7_0 - e_3 * fs_177_2431_3 * r_2 * h7_p2 - e_3 * fs_1_39_14 * r_4 * h5_0 - e_3 * fs_1_78_30 * r_4 * h5_p2 - e_3 * fs_25_429_14 * r_6 * h3_0 + e_3 * fs_5_858_210 * r_6 * h3_p2;
    }

    // NOTE: the angular components are formed in 49 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph3_p1, ph3_p3, ph5_p1, ph5_p3, ph7_p1, ph7_p3, ph9_p1, ph9_p3, ab_2, pc_60 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h3_p1 = ph3_p1[k];
        const auto h3_p3 = ph3_p3[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h7_p3 = ph7_p3[k];
        const auto h9_p1 = ph9_p1[k];
        const auto h9_p3 = ph9_p3[k];

        pc_60[k] = e_0 * fs_15_8_210 * h3_p1 - e_0 * fs_15_8_14 * h3_p3 - e_1 * fs_5_4_21 * h5_p1 - e_1 * fs_35_8_2 * h5_p3 - e_1 * fs_5_4_210 * r_2 * h3_p1 + e_1 * fs_5_4_14 * r_2 * h3_p3 + e_2 * fs_69_286_5 * h7_p1 + e_2 * fs_48_143_15 * h7_p3 + e_2 * fs_5_13_21 * r_2 * h5_p1 + e_2 * fs_35_26_2 * r_2 * h5_p3 + e_2 * fs_5_22_210 * r_4 * h3_p1 - e_2 * fs_5_22_14 * r_4 * h3_p3 + e_3 * fs_15_221_7 * h9_p1 - e_3 * fs_15_442_66 * h9_p3 - e_3 * fs_69_2431_5 * r_2 * h7_p1 - e_3 * fs_96_2431_15 * r_2 * h7_p3 - e_3 * fs_1_39_21 * r_4 * h5_p1 - e_3 * fs_7_78_2 * r_4 * h5_p3 - e_3 * fs_5_429_210 * r_6 * h3_p1 + e_3 * fs_5_429_14 * r_6 * h3_p3;
    }

    // NOTE: the angular components are formed in 49 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph3_p2, ph5_p2, ph5_p4, ph7_p2, ph7_p4, ph9_p2, ph9_p4, ab_2, pc_61 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

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

        pc_61[k] = e_0 * fs_45_8_21 * h3_p2 - e_1 * fs_15_4_3 * h5_p2 - e_1 * f_15_2 * h5_p4 - e_1 * fs_15_4_21 * r_2 * h3_p2 + e_2 * fs_111_572_30 * h7_p2 + e_2 * fs_27_286_165 * h7_p4 + e_2 * fs_15_13_3 * r_2 * h5_p2 + e_2 * f_30_13 * r_2 * h5_p4 + e_2 * fs_15_22_21 * r_4 * h3_p2 + e_3 * fs_5_442_154 * h9_p2 - e_3 * fs_5_221_143 * h9_p4 - e_3 * fs_111_4862_30 * r_2 * h7_p2 - e_3 * fs_27_2431_165 * r_2 * h7_p4 - e_3 * fs_1_13_3 * r_4 * h5_p2 - e_3 * f_2_13 * r_4 * h5_p4 - e_3 * fs_5_143_21 * r_6 * h3_p2;
    }

    // NOTE: the angular components are formed in 49 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph3_p3, ph5_p3, ph5_p5, ph7_p3, ph7_p4, ph7_p5, ph7_p6, ph9_p3, ph9_p4, ph9_p5, ph9_p6, ab_2, pc_62, pc_63 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h3_p3 = ph3_p3[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h5_p5 = ph5_p5[k];
        const auto h7_p3 = ph7_p3[k];
        const auto h7_p4 = ph7_p4[k];
        const auto h7_p5 = ph7_p5[k];
        const auto h7_p6 = ph7_p6[k];
        const auto h9_p3 = ph9_p3[k];
        const auto h9_p4 = ph9_p4[k];
        const auto h9_p5 = ph9_p5[k];
        const auto h9_p6 = ph9_p6[k];

        pc_62[k] = e_0 * fs_15_8_105 * h3_p3 - e_1 * fs_5_4_15 * h5_p3 - e_1 * fs_15_4_3 * h5_p5 - e_1 * fs_5_4_105 * r_2 * h3_p3 + e_2 * fs_573_572_2 * h7_p3 + e_2 * fs_87_572_22 * h7_p5 + e_2 * fs_5_13_15 * r_2 * h5_p3 + e_2 * fs_15_13_3 * r_2 * h5_p5 + e_2 * fs_5_22_105 * r_4 * h3_p3 + e_3 * fs_3_221_55 * h9_p3 - e_3 * fs_1_221_3003 * h9_p5 - e_3 * fs_573_4862_2 * r_2 * h7_p3 - e_3 * fs_87_4862_22 * r_2 * h7_p5 - e_3 * fs_1_39_15 * r_4 * h5_p3 - e_3 * fs_1_13_3 * r_4 * h5_p5 - e_3 * fs_5_429_105 * r_6 * h3_p3;

        pc_63[k] = e_2 * f_3_2 * h7_p4 - e_2 * fs_3_52_26 * h7_p6 + e_3 * fs_1_221_195 * h9_p4 - e_3 * fs_3_442_910 * h9_p6 - e_3 * f_3_17 * r_2 * h7_p4 + e_3 * fs_3_442_26 * r_2 * h7_p6;
    }

    // NOTE: the angular components are formed in 49 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, ph5_m4, ph5_p5, ph7_m4, ph7_p5, ph7_p7, ph9_m8, ph9_m4, ph9_p5, ph9_p7, ab_2, pc_64, pc_65 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

        const auto h5_m4 = ph5_m4[k];
        const auto h5_p5 = ph5_p5[k];
        const auto h7_m4 = ph7_m4[k];
        const auto h7_p5 = ph7_p5[k];
        const auto h7_p7 = ph7_p7[k];
        const auto h9_m8 = ph9_m8[k];
        const auto h9_m4 = ph9_m4[k];
        const auto h9_p5 = ph9_p5[k];
        const auto h9_p7 = ph9_p7[k];

        pc_64[k] = e_1 * fs_15_8_22 * h5_p5 + e_2 * fs_9_13_3 * h7_p5 - e_2 * fs_3_26_273 * h7_p7 - e_2 * fs_15_26_22 * r_2 * h5_p5 + e_3 * fs_1_442_182 * h9_p5 - e_3 * fs_1_221_910 * h9_p7 - e_3 * fs_18_221_3 * r_2 * h7_p5 + e_3 * fs_3_221_273 * r_2 * h7_p7 + e_3 * fs_1_26_22 * r_4 * h5_p5;

        pc_65[k] = - e_1 * fs_15_8_22 * h5_m4 - e_2 * fs_3_26_30 * h7_m4 + e_2 * fs_15_26_22 * r_2 * h5_m4 - e_3 * fs_1_221_3094 * h9_m8 - e_3 * fs_1_442_26 * h9_m4 + e_3 * fs_3_221_30 * r_2 * h7_m4 - e_3 * fs_1_26_22 * r_4 * h5_m4;
    }

    // NOTE: the angular components are formed in 49 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph3_m3, ph3_m2, ph5_m3, ph5_m2, ph7_m7, ph7_m6, ph7_m3, ph7_m2, ph9_m7, ph9_m6, ph9_m3, ph9_m2, ab_2, pc_66, pc_67 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

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

        pc_66[k] = e_0 * fs_15_8_231 * h3_m3 - e_1 * fs_5_4_33 * h5_m3 - e_1 * fs_5_4_231 * r_2 * h3_m3 + e_2 * fs_3_52_910 * h7_m7 - e_2 * fs_57_572_110 * h7_m3 + e_2 * fs_5_13_33 * r_2 * h5_m3 + e_2 * fs_5_22_231 * r_4 * h3_m3 - e_3 * fs_4_221_273 * h9_m7 - e_3 * f_6_221 * h9_m3 - e_3 * fs_3_442_910 * r_2 * h7_m7 + e_3 * fs_57_4862_110 * r_2 * h7_m3 - e_3 * fs_1_39_33 * r_4 * h5_m3 - e_3 * fs_5_429_231 * r_6 * h3_m3;

        pc_67[k] = e_0 * fs_45_4_7 * h3_m2 - e_1 * f_15_4 * h5_m2 - e_1 * fs_15_2_7 * r_2 * h3_m2 + e_2 * fs_6_143_1430 * h7_m6 - e_2 * fs_60_143_10 * h7_m2 + e_2 * f_15_13 * r_2 * h5_m2 + e_2 * fs_15_11_7 * r_4 * h3_m2 - e_3 * fs_3_442_2002 * h9_m6 - e_3 * fs_1_442_462 * h9_m2 - e_3 * fs_12_2431_1430 * r_2 * h7_m6 + e_3 * fs_120_2431_10 * r_2 * h7_m2 - e_3 * f_1_13 * r_4 * h5_m2 - e_3 * fs_10_143_7 * r_6 * h3_m2;
    }

    // NOTE: the angular components are formed in 49 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph3_m1, ph5_m5, ph5_m4, ph7_m5, ph7_m4, ph7_m1, ph9_m5, ph9_m4, ph9_m1, ab_2, pc_68, pc_69 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h3_m1 = ph3_m1[k];
        const auto h5_m5 = ph5_m5[k];
        const auto h5_m4 = ph5_m4[k];
        const auto h7_m5 = ph7_m5[k];
        const auto h7_m4 = ph7_m4[k];
        const auto h7_m1 = ph7_m1[k];
        const auto h9_m5 = ph9_m5[k];
        const auto h9_m4 = ph9_m4[k];
        const auto h9_m1 = ph9_m1[k];

        pc_68[k] = e_0 * fs_45_8_21 * h3_m1 + e_1 * f_15_4 * h5_m5 - e_1 * fs_15_4_21 * r_2 * h3_m1 + e_2 * fs_75_572_66 * h7_m5 - e_2 * fs_45_44_2 * h7_m1 - e_2 * f_15_13 * r_2 * h5_m5 + e_2 * fs_15_22_21 * r_4 * h3_m1 - e_3 * fs_2_221_1001 * h9_m5 - e_3 * fs_2_221_70 * h9_m1 - e_3 * fs_75_4862_66 * r_2 * h7_m5 + e_3 * fs_45_374_2 * r_2 * h7_m1 + e_3 * f_1_13 * r_4 * h5_m5 - e_3 * fs_5_143_21 * r_6 * h3_m1;

        pc_69[k] = e_1 * fs_15_8_10 * h5_m4 + e_2 * fs_15_286_66 * h7_m4 - e_2 * fs_15_26_10 * r_2 * h5_m4 - e_3 * fs_3_442_1430 * h9_m4 - e_3 * fs_15_2431_66 * r_2 * h7_m4 + e_3 * fs_1_26_10 * r_4 * h5_m4;
    }

    // NOTE: the angular components are formed in 49 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph3_m3, ph3_m1, ph3_p2, ph5_m3, ph5_m1, ph5_p2, ph7_m3, ph7_m1, ph7_p2, ph9_m3, ph9_m1, ph9_p2, ab_2, pc_70, pc_71 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h3_m3 = ph3_m3[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h3_p2 = ph3_p2[k];
        const auto h5_m3 = ph5_m3[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h7_m3 = ph7_m3[k];
        const auto h7_m1 = ph7_m1[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h9_m3 = ph9_m3[k];
        const auto h9_m1 = ph9_m1[k];
        const auto h9_p2 = ph9_p2[k];

        pc_70[k] = e_0 * fs_15_16_14 * h3_m3 + e_0 * fs_15_16_210 * h3_m1 + e_1 * fs_5_1_2 * h5_m3 + e_1 * fs_5_4_21 * h5_m1 - e_1 * fs_5_8_14 * r_2 * h3_m3 - e_1 * fs_5_8_210 * r_2 * h3_m1 - e_2 * fs_15_286_15 * h7_m3 - e_2 * fs_147_286_5 * h7_m1 - e_2 * fs_20_13_2 * r_2 * h5_m3 - e_2 * fs_5_13_21 * r_2 * h5_m1 + e_2 * fs_5_44_14 * r_4 * h3_m3 + e_2 * fs_5_44_210 * r_4 * h3_m1 - e_3 * fs_6_221_66 * h9_m3 - e_3 * fs_12_221_7 * h9_m1 + e_3 * fs_15_2431_15 * r_2 * h7_m3 + e_3 * fs_147_2431_5 * r_2 * h7_m1 + e_3 * fs_4_39_2 * r_4 * h5_m3 + e_3 * fs_1_39_21 * r_4 * h5_m1 - e_3 * fs_5_858_14 * r_6 * h3_m3 - e_3 * fs_5_858_210 * r_6 * h3_m1;

        pc_71[k] = e_0 * f_45_4 * h3_p2 + e_1 * fs_15_4_7 * h5_p2 - e_1 * f_15_2 * r_2 * h3_p2 - e_2 * fs_18_143_70 * h7_p2 - e_2 * fs_15_13_7 * r_2 * h5_p2 + e_2 * f_15_11 * r_4 * h3_p2 - e_3 * fs_7_221_66 * h9_p2 + e_3 * fs_36_2431_70 * r_2 * h7_p2 + e_3 * fs_1_13_7 * r_4 * h5_p2 - e_3 * f_10_143 * r_6 * h3_p2;
    }

    // NOTE: the angular components are formed in 49 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph3_p1, ph3_p3, ph5_p1, ph5_p3, ph7_p1, ph7_p3, ph9_p1, ph9_p3, ab_2, pc_72 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h3_p1 = ph3_p1[k];
        const auto h3_p3 = ph3_p3[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h7_p3 = ph7_p3[k];
        const auto h9_p1 = ph9_p1[k];
        const auto h9_p3 = ph9_p3[k];

        pc_72[k] = - e_0 * fs_15_16_210 * h3_p1 + e_0 * fs_15_16_14 * h3_p3 - e_1 * fs_5_4_21 * h5_p1 + e_1 * fs_5_1_2 * h5_p3 + e_1 * fs_5_8_210 * r_2 * h3_p1 - e_1 * fs_5_8_14 * r_2 * h3_p3 + e_2 * fs_147_286_5 * h7_p1 - e_2 * fs_15_286_15 * h7_p3 + e_2 * fs_5_13_21 * r_2 * h5_p1 - e_2 * fs_20_13_2 * r_2 * h5_p3 - e_2 * fs_5_44_210 * r_4 * h3_p1 + e_2 * fs_5_44_14 * r_4 * h3_p3 + e_3 * fs_12_221_7 * h9_p1 - e_3 * fs_6_221_66 * h9_p3 - e_3 * fs_147_2431_5 * r_2 * h7_p1 + e_3 * fs_15_2431_15 * r_2 * h7_p3 - e_3 * fs_1_39_21 * r_4 * h5_p1 + e_3 * fs_4_39_2 * r_4 * h5_p3 + e_3 * fs_5_858_210 * r_6 * h3_p1 - e_3 * fs_5_858_14 * r_6 * h3_p3;
    }

    // NOTE: the angular components are formed in 49 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph3_0, ph3_p1, ph5_0, ph5_p4, ph5_p5, ph7_0, ph7_p1, ph7_p4, ph7_p5, ph9_0, ph9_p1, ph9_p4, ph9_p5, ab_2, pc_73, pc_74 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h3_0 = ph3_0[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h5_0 = ph5_0[k];
        const auto h5_p4 = ph5_p4[k];
        const auto h5_p5 = ph5_p5[k];
        const auto h7_0 = ph7_0[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h7_p4 = ph7_p4[k];
        const auto h7_p5 = ph7_p5[k];
        const auto h9_0 = ph9_0[k];
        const auto h9_p1 = ph9_p1[k];
        const auto h9_p4 = ph9_p4[k];
        const auto h9_p5 = ph9_p5[k];

        pc_73[k] = e_0 * fs_15_2_14 * h3_0 + e_1 * fs_5_4_14 * h5_0 + e_1 * fs_15_8_10 * h5_p4 - e_1 * fs_5_1_14 * r_2 * h3_0 - e_2 * fs_75_143_14 * h7_0 + e_2 * fs_15_286_66 * h7_p4 - e_2 * fs_5_13_14 * r_2 * h5_0 - e_2 * fs_15_26_10 * r_2 * h5_p4 + e_2 * fs_10_11_14 * r_4 * h3_0 - e_3 * fs_9_221_14 * h9_0 - e_3 * fs_3_442_1430 * h9_p4 + e_3 * fs_150_2431_14 * r_2 * h7_0 - e_3 * fs_15_2431_66 * r_2 * h7_p4 + e_3 * fs_1_39_14 * r_4 * h5_0 + e_3 * fs_1_26_10 * r_4 * h5_p4 - e_3 * fs_20_429_14 * r_6 * h3_0;

        pc_74[k] = e_0 * fs_45_8_21 * h3_p1 + e_1 * f_15_4 * h5_p5 - e_1 * fs_15_4_21 * r_2 * h3_p1 - e_2 * fs_45_44_2 * h7_p1 + e_2 * fs_75_572_66 * h7_p5 - e_2 * f_15_13 * r_2 * h5_p5 + e_2 * fs_15_22_21 * r_4 * h3_p1 - e_3 * fs_2_221_70 * h9_p1 - e_3 * fs_2_221_1001 * h9_p5 + e_3 * fs_45_374_2 * r_2 * h7_p1 - e_3 * fs_75_4862_66 * r_2 * h7_p5 + e_3 * f_1_13 * r_4 * h5_p5 - e_3 * fs_5_143_21 * r_6 * h3_p1;
    }

    // NOTE: the angular components are formed in 49 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph3_p2, ph3_p3, ph5_p2, ph5_p3, ph7_p2, ph7_p3, ph7_p6, ph7_p7, ph9_p2, ph9_p3, ph9_p6, ph9_p7, ab_2, pc_75, pc_76 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

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

        pc_75[k] = e_0 * fs_45_4_7 * h3_p2 - e_1 * f_15_4 * h5_p2 - e_1 * fs_15_2_7 * r_2 * h3_p2 - e_2 * fs_60_143_10 * h7_p2 + e_2 * fs_6_143_1430 * h7_p6 + e_2 * f_15_13 * r_2 * h5_p2 + e_2 * fs_15_11_7 * r_4 * h3_p2 - e_3 * fs_1_442_462 * h9_p2 - e_3 * fs_3_442_2002 * h9_p6 + e_3 * fs_120_2431_10 * r_2 * h7_p2 - e_3 * fs_12_2431_1430 * r_2 * h7_p6 - e_3 * f_1_13 * r_4 * h5_p2 - e_3 * fs_10_143_7 * r_6 * h3_p2;

        pc_76[k] = e_0 * fs_15_8_231 * h3_p3 - e_1 * fs_5_4_33 * h5_p3 - e_1 * fs_5_4_231 * r_2 * h3_p3 - e_2 * fs_57_572_110 * h7_p3 + e_2 * fs_3_52_910 * h7_p7 + e_2 * fs_5_13_33 * r_2 * h5_p3 + e_2 * fs_5_22_231 * r_4 * h3_p3 - e_3 * f_6_221 * h9_p3 - e_3 * fs_4_221_273 * h9_p7 + e_3 * fs_57_4862_110 * r_2 * h7_p3 - e_3 * fs_3_442_910 * r_2 * h7_p7 - e_3 * fs_1_39_33 * r_4 * h5_p3 - e_3 * fs_5_429_231 * r_6 * h3_p3;
    }

    // NOTE: the angular components are formed in 49 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph3_m3, ph5_m3, ph5_p4, ph7_m3, ph7_p4, ph9_m9, ph9_m3, ph9_p4, ph9_p8, ab_2, pc_77, pc_78 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h3_m3 = ph3_m3[k];
        const auto h5_m3 = ph5_m3[k];
        const auto h5_p4 = ph5_p4[k];
        const auto h7_m3 = ph7_m3[k];
        const auto h7_p4 = ph7_p4[k];
        const auto h9_m9 = ph9_m9[k];
        const auto h9_m3 = ph9_m3[k];
        const auto h9_p4 = ph9_p4[k];
        const auto h9_p8 = ph9_p8[k];

        pc_77[k] = - e_1 * fs_15_8_22 * h5_p4 - e_2 * fs_3_26_30 * h7_p4 + e_2 * fs_15_26_22 * r_2 * h5_p4 - e_3 * fs_1_442_26 * h9_p4 - e_3 * fs_1_221_3094 * h9_p8 + e_3 * fs_3_221_30 * r_2 * h7_p4 - e_3 * fs_1_26_22 * r_4 * h5_p4;

        pc_78[k] = e_0 * fs_15_8_462 * h3_m3 + e_1 * fs_5_8_66 * h5_m3 - e_1 * fs_5_4_462 * r_2 * h3_m3 + e_2 * fs_9_286_55 * h7_m3 - e_2 * fs_5_26_66 * r_2 * h5_m3 + e_2 * fs_5_22_462 * r_4 * h3_m3 - e_3 * fs_1_221_9282 * h9_m9 + e_3 * fs_1_442_2 * h9_m3 - e_3 * fs_9_2431_55 * r_2 * h7_m3 + e_3 * fs_1_78_66 * r_4 * h5_m3 - e_3 * fs_5_429_462 * r_6 * h3_m3;
    }

    // NOTE: the angular components are formed in 49 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph3_m2, ph3_m1, ph5_m2, ph5_m1, ph7_m7, ph7_m2, ph7_m1, ph9_m8, ph9_m7, ph9_m2, ph9_m1, ab_2, pc_79, pc_80 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

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

        pc_79[k] = e_0 * fs_15_8_231 * h3_m2 + e_1 * fs_5_4_33 * h5_m2 - e_1 * fs_5_4_231 * r_2 * h3_m2 + e_2 * fs_15_572_330 * h7_m2 - e_2 * fs_5_13_33 * r_2 * h5_m2 + e_2 * fs_5_22_231 * r_4 * h3_m2 - e_3 * fs_2_221_1547 * h9_m8 + e_3 * fs_1_442_14 * h9_m2 - e_3 * fs_15_4862_330 * r_2 * h7_m2 + e_3 * fs_1_39_33 * r_4 * h5_m2 - e_3 * fs_5_429_231 * r_6 * h3_m2;

        pc_80[k] = e_0 * fs_15_8_105 * h3_m1 + e_1 * fs_5_4_42 * h5_m1 - e_1 * fs_5_4_105 * r_2 * h3_m1 - e_2 * fs_3_572_30030 * h7_m7 + e_2 * fs_135_572_10 * h7_m1 - e_2 * fs_5_13_42 * r_2 * h5_m1 + e_2 * fs_5_22_105 * r_4 * h3_m1 - e_3 * fs_2_221_1001 * h9_m7 + e_3 * fs_1_221_14 * h9_m1 + e_3 * fs_3_4862_30030 * r_2 * h7_m7 - e_3 * fs_135_4862_10 * r_2 * h7_m1 + e_3 * fs_1_39_42 * r_4 * h5_m1 - e_3 * fs_5_429_105 * r_6 * h3_m1;
    }

    // NOTE: the angular components are formed in 49 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph3_m1, ph5_m5, ph5_m1, ph7_m6, ph7_m5, ph7_m1, ph9_m6, ph9_m5, ph9_m1, ab_2, pc_81, pc_82 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h3_m1 = ph3_m1[k];
        const auto h5_m5 = ph5_m5[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h7_m6 = ph7_m6[k];
        const auto h7_m5 = ph7_m5[k];
        const auto h7_m1 = ph7_m1[k];
        const auto h9_m6 = ph9_m6[k];
        const auto h9_m5 = ph9_m5[k];
        const auto h9_m1 = ph9_m1[k];

        pc_81[k] = - e_2 * fs_45_572_286 * h7_m6 - e_3 * fs_1_442_10010 * h9_m6 + e_3 * fs_45_4862_286 * r_2 * h7_m6;

        pc_82[k] = e_0 * fs_15_8_14 * h3_m1 - e_1 * fs_5_8_6 * h5_m5 + e_1 * fs_5_4_35 * h5_m1 - e_1 * fs_5_4_14 * r_2 * h3_m1 - e_2 * fs_135_286_11 * h7_m5 + e_2 * fs_105_143_3 * h7_m1 + e_2 * fs_5_26_6 * r_2 * h5_m5 - e_2 * fs_5_13_35 * r_2 * h5_m1 + e_2 * fs_5_22_14 * r_4 * h3_m1 - e_3 * fs_1_442_6006 * h9_m5 + e_3 * fs_1_221_105 * h9_m1 + e_3 * fs_135_2431_11 * r_2 * h7_m5 - e_3 * fs_210_2431_3 * r_2 * h7_m1 - e_3 * fs_1_78_6 * r_4 * h5_m5 + e_3 * fs_1_39_35 * r_4 * h5_m1 - e_3 * fs_5_429_14 * r_6 * h3_m1;
    }

    // NOTE: the angular components are formed in 49 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph3_m2, ph3_p3, ph5_m4, ph5_m2, ph5_p3, ph7_m4, ph7_m2, ph7_p3, ph9_m4, ph9_m2, ph9_p3, ab_2, pc_83, pc_84 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

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

        pc_83[k] = - e_0 * fs_15_16_14 * h3_m2 - e_1 * fs_5_4_6 * h5_m4 - e_1 * fs_35_8_2 * h5_m2 + e_1 * fs_5_8_14 * r_2 * h3_m2 - e_2 * fs_45_286_110 * h7_m4 - e_2 * fs_189_286_5 * h7_m2 + e_2 * fs_5_13_6 * r_2 * h5_m4 + e_2 * fs_35_26_2 * r_2 * h5_m2 - e_2 * fs_5_44_14 * r_4 * h3_m2 - e_3 * fs_1_221_858 * h9_m4 - e_3 * fs_1_221_231 * h9_m2 + e_3 * fs_45_2431_110 * r_2 * h7_m4 + e_3 * fs_189_2431_5 * r_2 * h7_m2 - e_3 * fs_1_39_6 * r_4 * h5_m4 - e_3 * fs_7_78_2 * r_4 * h5_m2 + e_3 * fs_5_858_14 * r_6 * h3_m2;

        pc_84[k] = - e_0 * f_15_8 * h3_p3 - e_1 * fs_5_2_7 * h5_p3 + e_1 * f_5_4 * r_2 * h3_p3 - e_2 * fs_45_286_210 * h7_p3 + e_2 * fs_10_13_7 * r_2 * h5_p3 - e_2 * f_5_22 * r_4 * h3_p3 - e_3 * fs_2_221_231 * h9_p3 + e_3 * fs_45_2431_210 * r_2 * h7_p3 - e_3 * fs_2_39_7 * r_4 * h5_p3 + e_3 * f_5_429 * r_6 * h3_p3;
    }

    // NOTE: the angular components are formed in 49 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph3_p2, ph5_p2, ph5_p4, ph7_p2, ph7_p4, ph9_p2, ph9_p4, ab_2, pc_85 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

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

        pc_85[k] = e_0 * fs_15_16_14 * h3_p2 + e_1 * fs_35_8_2 * h5_p2 - e_1 * fs_5_4_6 * h5_p4 - e_1 * fs_5_8_14 * r_2 * h3_p2 + e_2 * fs_189_286_5 * h7_p2 - e_2 * fs_45_286_110 * h7_p4 - e_2 * fs_35_26_2 * r_2 * h5_p2 + e_2 * fs_5_13_6 * r_2 * h5_p4 + e_2 * fs_5_44_14 * r_4 * h3_p2 + e_3 * fs_1_221_231 * h9_p2 - e_3 * fs_1_221_858 * h9_p4 - e_3 * fs_189_2431_5 * r_2 * h7_p2 + e_3 * fs_45_2431_110 * r_2 * h7_p4 + e_3 * fs_7_78_2 * r_4 * h5_p2 - e_3 * fs_1_39_6 * r_4 * h5_p4 - e_3 * fs_5_858_14 * r_6 * h3_p2;
    }

    // NOTE: the angular components are formed in 49 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph3_0, ph3_p1, ph5_0, ph5_p1, ph5_p5, ph7_0, ph7_p1, ph7_p5, ph7_p6, ph9_0, ph9_p1, ph9_p5, ph9_p6, ab_2, pc_86, pc_87 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h3_0 = ph3_0[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h5_0 = ph5_0[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h5_p5 = ph5_p5[k];
        const auto h7_0 = ph7_0[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h7_p5 = ph7_p5[k];
        const auto h7_p6 = ph7_p6[k];
        const auto h9_0 = ph9_0[k];
        const auto h9_p1 = ph9_p1[k];
        const auto h9_p5 = ph9_p5[k];
        const auto h9_p6 = ph9_p6[k];

        pc_86[k] = - e_0 * fs_15_8_14 * h3_p1 - e_1 * fs_5_4_35 * h5_p1 - e_1 * fs_5_8_6 * h5_p5 + e_1 * fs_5_4_14 * r_2 * h3_p1 - e_2 * fs_105_143_3 * h7_p1 - e_2 * fs_135_286_11 * h7_p5 + e_2 * fs_5_13_35 * r_2 * h5_p1 + e_2 * fs_5_26_6 * r_2 * h5_p5 - e_2 * fs_5_22_14 * r_4 * h3_p1 - e_3 * fs_1_221_105 * h9_p1 - e_3 * fs_1_442_6006 * h9_p5 + e_3 * fs_210_2431_3 * r_2 * h7_p1 + e_3 * fs_135_2431_11 * r_2 * h7_p5 - e_3 * fs_1_39_35 * r_4 * h5_p1 - e_3 * fs_1_78_6 * r_4 * h5_p5 + e_3 * fs_5_429_14 * r_6 * h3_p1;

        pc_87[k] = e_0 * fs_15_4_21 * h3_0 + e_1 * fs_5_2_21 * h5_0 - e_1 * fs_5_2_21 * r_2 * h3_0 + e_2 * fs_45_143_21 * h7_0 - e_2 * fs_45_572_286 * h7_p6 - e_2 * fs_10_13_21 * r_2 * h5_0 + e_2 * fs_5_11_21 * r_4 * h3_0 + e_3 * fs_2_221_21 * h9_0 - e_3 * fs_1_442_10010 * h9_p6 - e_3 * fs_90_2431_21 * r_2 * h7_0 + e_3 * fs_45_4862_286 * r_2 * h7_p6 + e_3 * fs_2_39_21 * r_4 * h5_0 - e_3 * fs_10_429_21 * r_6 * h3_0;
    }

    // NOTE: the angular components are formed in 49 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph3_p1, ph3_p2, ph5_p1, ph5_p2, ph7_p1, ph7_p2, ph7_p7, ph9_p1, ph9_p2, ph9_p7, ph9_p8, ab_2, pc_88, pc_89 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

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

        pc_88[k] = e_0 * fs_15_8_105 * h3_p1 + e_1 * fs_5_4_42 * h5_p1 - e_1 * fs_5_4_105 * r_2 * h3_p1 + e_2 * fs_135_572_10 * h7_p1 - e_2 * fs_3_572_30030 * h7_p7 - e_2 * fs_5_13_42 * r_2 * h5_p1 + e_2 * fs_5_22_105 * r_4 * h3_p1 + e_3 * fs_1_221_14 * h9_p1 - e_3 * fs_2_221_1001 * h9_p7 - e_3 * fs_135_4862_10 * r_2 * h7_p1 + e_3 * fs_3_4862_30030 * r_2 * h7_p7 + e_3 * fs_1_39_42 * r_4 * h5_p1 - e_3 * fs_5_429_105 * r_6 * h3_p1;

        pc_89[k] = e_0 * fs_15_8_231 * h3_p2 + e_1 * fs_5_4_33 * h5_p2 - e_1 * fs_5_4_231 * r_2 * h3_p2 + e_2 * fs_15_572_330 * h7_p2 - e_2 * fs_5_13_33 * r_2 * h5_p2 + e_2 * fs_5_22_231 * r_4 * h3_p2 + e_3 * fs_1_442_14 * h9_p2 - e_3 * fs_2_221_1547 * h9_p8 - e_3 * fs_15_4862_330 * r_2 * h7_p2 + e_3 * fs_1_39_33 * r_4 * h5_p2 - e_3 * fs_5_429_231 * r_6 * h3_p2;
    }

    // NOTE: the angular components are formed in 49 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph3_p3, ph5_p3, ph7_p3, ph9_p3, ph9_p9, ab_2, pc_90 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h3_p3 = ph3_p3[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h7_p3 = ph7_p3[k];
        const auto h9_p3 = ph9_p3[k];
        const auto h9_p9 = ph9_p9[k];

        pc_90[k] = e_0 * fs_15_8_462 * h3_p3 + e_1 * fs_5_8_66 * h5_p3 - e_1 * fs_5_4_462 * r_2 * h3_p3 + e_2 * fs_9_286_55 * h7_p3 - e_2 * fs_5_26_66 * r_2 * h5_p3 + e_2 * fs_5_22_462 * r_4 * h3_p3 + e_3 * fs_1_442_2 * h9_p3 - e_3 * fs_1_221_9282 * h9_p9 - e_3 * fs_9_2431_55 * r_2 * h7_p3 + e_3 * fs_1_78_66 * r_4 * h5_p3 - e_3 * fs_5_429_462 * r_6 * h3_p3;
    }

    // NOTE: the values of a combination of angular components are stored as one
    // row of nvalues columns, with the component on bra side running slowest. The
    // rows which the symmetry relates to an already formed one are copied from it.

    const size_t sources[91] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90};

    for (size_t n = 0; n < 91; n++)
    {
        if (sources[n] != n) std::copy(values + sources[n] * nvalues, values + (sources[n] + 1) * nvalues, values + n * nvalues);
    }
}

}  // namespace simdt2ceri
