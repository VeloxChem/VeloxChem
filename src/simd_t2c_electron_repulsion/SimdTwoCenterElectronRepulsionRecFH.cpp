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



#include "SimdTwoCenterElectronRepulsionRecFH.hpp"

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
compute_fh_electron_repulsion(double                         *values,
                              const size_t                    nvalues,
                              const CBasisFunction           &bra,
                              const CBasisFunction           &ket,
                              const std::vector<CSimdMatrix> &harmonics,
                              const CSimdMatrix              &coordinates) -> void
{
    if ((bra.get_angular_momentum() != 3) || (ket.get_angular_momentum() != 5))
    {
        errors::assertMsgCritical(
            false, std::string("SimdTwoCenterElectronRepulsionFunc.compute_fh_electron_repulsion: Basis functions must be of angular momenta three and five"));
    }

    if (harmonics.size() < 8)
    {
        errors::assertMsgCritical(
            false, std::string("SimdTwoCenterElectronRepulsionFunc.compute_fh_electron_repulsion: Harmonics must reach angular momentum 8"));
    }

    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("SimdTwoCenterElectronRepulsionFunc.compute_fh_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    auto boys = CSimdVariableMatrix(std::vector<size_t>(nprims, nvalues), 10);

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
    // call, which fills the orders 0 to 8 of every row. The terms read the
    // orders 5 to 8 alone, and the orders below them are formed on the
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

            const auto ff_0 = fbase * aexp * aexp / (fexp * fexp * fexp * fexp * fexp);

            const auto ff_1 = fbase * aexp * aexp * fmu / (fexp * fexp * fexp * fexp * fexp);

            const auto ff_2 = fbase * aexp * aexp * fmu * fmu / (fexp * fexp * fexp * fexp * fexp);

            const auto ff_3 = fbase * aexp * aexp * fmu * fmu * fmu / (fexp * fexp * fexp * fexp * fexp);

            const auto *bv_0 = boys.data(6, i * nprim_b + j);

            const auto *bv_1 = boys.data(7, i * nprim_b + j);

            const auto *bv_2 = boys.data(8, i * nprim_b + j);

            const auto *bv_3 = boys.data(9, i * nprim_b + j);

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

    // NOTE: the factors of the terms depend on the angular momenta alone, so they
    // are formed once for the whole matrix instead of once for every atom pair.

    const auto f_15_11 = 15.0 / 11.0;
    const auto f_15_22 = 15.0 / 22.0;
    const auto f_180_1001 = 180.0 / 1001.0;
    const auto f_180_77 = 180.0 / 77.0;
    const auto f_18_11 = 18.0 / 11.0;
    const auto f_18_143 = 18.0 / 143.0;
    const auto f_1_11 = 1.0 / 11.0;
    const auto f_21_11 = 21.0 / 11.0;
    const auto f_21_143 = 21.0 / 143.0;
    const auto f_21_4 = 21.0 / 4.0;
    const auto f_225_14 = 225.0 / 14.0;
    const auto f_25_7 = 25.0 / 7.0;
    const auto f_2_11 = 2.0 / 11.0;
    const auto f_35_22 = 35.0 / 22.0;
    const auto f_36_11 = 36.0 / 11.0;
    const auto f_36_143 = 36.0 / 143.0;
    const auto f_45_7 = 45.0 / 7.0;
    const auto f_50_231 = 50.0 / 231.0;
    const auto f_56_143 = 56.0 / 143.0;
    const auto f_75_4 = 75.0 / 4.0;
    const auto f_7_33 = 7.0 / 33.0;
    const auto f_9_1 = 9.0;
    const auto f_9_2 = 9.0 / 2.0;
    const auto fs_10_11_2 = std::sqrt(200.0 / 121.0);
    const auto fs_10_231_14 = std::sqrt(200.0 / 7623.0);
    const auto fs_10_231_3 = std::sqrt(100.0 / 17787.0);
    const auto fs_10_231_7 = std::sqrt(100.0 / 7623.0);
    const auto fs_10_7_5 = std::sqrt(500.0 / 49.0);
    const auto fs_120_1001_3 = std::sqrt(43200.0 / 1002001.0);
    const auto fs_120_77_3 = std::sqrt(43200.0 / 5929.0);
    const auto fs_129_154_2 = std::sqrt(16641.0 / 11858.0);
    const auto fs_129_2002_2 = std::sqrt(16641.0 / 2004002.0);
    const auto fs_129_56_2 = std::sqrt(16641.0 / 1568.0);
    const auto fs_12_1001_105 = std::sqrt(2160.0 / 143143.0);
    const auto fs_12_11_3 = std::sqrt(432.0 / 121.0);
    const auto fs_12_143_11 = std::sqrt(144.0 / 1859.0);
    const auto fs_12_143_3 = std::sqrt(432.0 / 20449.0);
    const auto fs_12_77_105 = std::sqrt(2160.0 / 847.0);
    const auto fs_12_7_7 = std::sqrt(144.0 / 7.0);
    const auto fs_135_28_10 = std::sqrt(91125.0 / 392.0);
    const auto fs_135_28_7 = std::sqrt(18225.0 / 112.0);
    const auto fs_14_143_10 = std::sqrt(1960.0 / 20449.0);
    const auto fs_14_143_15 = std::sqrt(2940.0 / 20449.0);
    const auto fs_15_1001_105 = std::sqrt(3375.0 / 143143.0);
    const auto fs_15_11_3 = std::sqrt(675.0 / 121.0);
    const auto fs_15_143_3 = std::sqrt(675.0 / 20449.0);
    const auto fs_15_143_7 = std::sqrt(1575.0 / 20449.0);
    const auto fs_15_14_10 = std::sqrt(1125.0 / 98.0);
    const auto fs_15_14_7 = std::sqrt(225.0 / 28.0);
    const auto fs_15_16_14 = std::sqrt(1575.0 / 128.0);
    const auto fs_15_16_2 = std::sqrt(225.0 / 128.0);
    const auto fs_15_16_210 = std::sqrt(23625.0 / 128.0);
    const auto fs_15_16_30 = std::sqrt(3375.0 / 128.0);
    const auto fs_15_22_11 = std::sqrt(225.0 / 44.0);
    const auto fs_15_22_3 = std::sqrt(675.0 / 484.0);
    const auto fs_15_22_5 = std::sqrt(1125.0 / 484.0);
    const auto fs_15_22_7 = std::sqrt(1575.0 / 484.0);
    const auto fs_15_286_6 = std::sqrt(675.0 / 40898.0);
    const auto fs_15_28_105 = std::sqrt(3375.0 / 112.0);
    const auto fs_15_2_5 = std::sqrt(1125.0 / 4.0);
    const auto fs_15_44_10 = std::sqrt(1125.0 / 968.0);
    const auto fs_15_44_22 = std::sqrt(225.0 / 88.0);
    const auto fs_15_4_14 = std::sqrt(1575.0 / 8.0);
    const auto fs_15_4_3 = std::sqrt(675.0 / 16.0);
    const auto fs_15_4_7 = std::sqrt(1575.0 / 16.0);
    const auto fs_15_77_105 = std::sqrt(3375.0 / 847.0);
    const auto fs_15_8_105 = std::sqrt(23625.0 / 64.0);
    const auto fs_15_8_35 = std::sqrt(7875.0 / 64.0);
    const auto fs_15_8_42 = std::sqrt(4725.0 / 32.0);
    const auto fs_15_8_5 = std::sqrt(1125.0 / 64.0);
    const auto fs_1_11_11 = std::sqrt(1.0 / 11.0);
    const auto fs_1_11_3 = std::sqrt(3.0 / 121.0);
    const auto fs_1_11_5 = std::sqrt(5.0 / 121.0);
    const auto fs_1_11_7 = std::sqrt(7.0 / 121.0);
    const auto fs_1_143_105 = std::sqrt(105.0 / 20449.0);
    const auto fs_1_143_1155 = std::sqrt(105.0 / 1859.0);
    const auto fs_1_143_2002 = std::sqrt(14.0 / 143.0);
    const auto fs_1_143_286 = std::sqrt(2.0 / 143.0);
    const auto fs_1_143_30 = std::sqrt(30.0 / 20449.0);
    const auto fs_1_143_33 = std::sqrt(3.0 / 1859.0);
    const auto fs_1_143_42 = std::sqrt(42.0 / 20449.0);
    const auto fs_1_143_462 = std::sqrt(42.0 / 1859.0);
    const auto fs_1_143_858 = std::sqrt(6.0 / 143.0);
    const auto fs_1_22_10 = std::sqrt(5.0 / 242.0);
    const auto fs_1_22_22 = std::sqrt(1.0 / 22.0);
    const auto fs_1_286_10010 = std::sqrt(35.0 / 286.0);
    const auto fs_1_286_14 = std::sqrt(7.0 / 40898.0);
    const auto fs_1_286_2 = std::sqrt(1.0 / 40898.0);
    const auto fs_1_286_2002 = std::sqrt(7.0 / 286.0);
    const auto fs_1_286_22 = std::sqrt(1.0 / 3718.0);
    const auto fs_1_286_4290 = std::sqrt(15.0 / 286.0);
    const auto fs_1_286_6006 = std::sqrt(21.0 / 286.0);
    const auto fs_1_33_14 = std::sqrt(14.0 / 1089.0);
    const auto fs_1_33_15 = std::sqrt(5.0 / 363.0);
    const auto fs_1_33_21 = std::sqrt(7.0 / 363.0);
    const auto fs_1_33_3 = std::sqrt(1.0 / 363.0);
    const auto fs_1_33_33 = std::sqrt(1.0 / 33.0);
    const auto fs_1_33_35 = std::sqrt(35.0 / 1089.0);
    const auto fs_1_33_42 = std::sqrt(14.0 / 363.0);
    const auto fs_1_33_6 = std::sqrt(2.0 / 363.0);
    const auto fs_1_66_30 = std::sqrt(5.0 / 726.0);
    const auto fs_1_66_6 = std::sqrt(1.0 / 726.0);
    const auto fs_1_66_66 = std::sqrt(1.0 / 66.0);
    const auto fs_20_143_6 = std::sqrt(2400.0 / 20449.0);
    const auto fs_20_231_5 = std::sqrt(2000.0 / 53361.0);
    const auto fs_21_143_6 = std::sqrt(2646.0 / 20449.0);
    const auto fs_225_28_2 = std::sqrt(50625.0 / 392.0);
    const auto fs_25_14_2 = std::sqrt(625.0 / 98.0);
    const auto fs_25_231_2 = std::sqrt(1250.0 / 53361.0);
    const auto fs_27_1001_21 = std::sqrt(2187.0 / 143143.0);
    const auto fs_27_1001_35 = std::sqrt(3645.0 / 143143.0);
    const auto fs_27_28_21 = std::sqrt(2187.0 / 112.0);
    const auto fs_27_28_35 = std::sqrt(3645.0 / 112.0);
    const auto fs_27_77_21 = std::sqrt(2187.0 / 847.0);
    const auto fs_27_77_35 = std::sqrt(3645.0 / 847.0);
    const auto fs_2_143_1001 = std::sqrt(28.0 / 143.0);
    const auto fs_2_143_7 = std::sqrt(28.0 / 20449.0);
    const auto fs_2_33_2 = std::sqrt(8.0 / 1089.0);
    const auto fs_2_33_21 = std::sqrt(28.0 / 363.0);
    const auto fs_2_33_7 = std::sqrt(28.0 / 1089.0);
    const auto fs_30_1001_7 = std::sqrt(900.0 / 143143.0);
    const auto fs_30_77_7 = std::sqrt(900.0 / 847.0);
    const auto fs_30_7_3 = std::sqrt(2700.0 / 49.0);
    const auto fs_33_28_10 = std::sqrt(5445.0 / 392.0);
    const auto fs_33_56_70 = std::sqrt(5445.0 / 224.0);
    const auto fs_35_44_2 = std::sqrt(1225.0 / 968.0);
    const auto fs_3_1001_105 = std::sqrt(135.0 / 143143.0);
    const auto fs_3_1001_15 = std::sqrt(135.0 / 1002001.0);
    const auto fs_3_1001_210 = std::sqrt(270.0 / 143143.0);
    const auto fs_3_11_3 = std::sqrt(27.0 / 121.0);
    const auto fs_3_143_165 = std::sqrt(135.0 / 1859.0);
    const auto fs_3_143_3 = std::sqrt(27.0 / 20449.0);
    const auto fs_3_143_70 = std::sqrt(630.0 / 20449.0);
    const auto fs_3_14_21 = std::sqrt(27.0 / 28.0);
    const auto fs_3_14_70 = std::sqrt(45.0 / 14.0);
    const auto fs_3_182_70 = std::sqrt(45.0 / 2366.0);
    const auto fs_3_1_3 = std::sqrt(27.0);
    const auto fs_3_286_42 = std::sqrt(189.0 / 40898.0);
    const auto fs_3_286_858 = std::sqrt(27.0 / 286.0);
    const auto fs_3_28_105 = std::sqrt(135.0 / 112.0);
    const auto fs_3_28_15 = std::sqrt(135.0 / 784.0);
    const auto fs_3_28_210 = std::sqrt(135.0 / 56.0);
    const auto fs_3_2_15 = std::sqrt(135.0 / 4.0);
    const auto fs_3_2_5 = std::sqrt(45.0 / 4.0);
    const auto fs_3_4_3 = std::sqrt(27.0 / 16.0);
    const auto fs_3_77_105 = std::sqrt(135.0 / 847.0);
    const auto fs_3_77_15 = std::sqrt(135.0 / 5929.0);
    const auto fs_3_77_210 = std::sqrt(270.0 / 847.0);
    const auto fs_3_7_10 = std::sqrt(90.0 / 49.0);
    const auto fs_3_7_105 = std::sqrt(135.0 / 7.0);
    const auto fs_3_91_10 = std::sqrt(90.0 / 8281.0);
    const auto fs_45_1001_7 = std::sqrt(2025.0 / 143143.0);
    const auto fs_45_14_14 = std::sqrt(2025.0 / 14.0);
    const auto fs_45_14_3 = std::sqrt(6075.0 / 196.0);
    const auto fs_45_14_7 = std::sqrt(2025.0 / 28.0);
    const auto fs_45_154_30 = std::sqrt(30375.0 / 11858.0);
    const auto fs_45_2002_30 = std::sqrt(30375.0 / 2004002.0);
    const auto fs_45_28_105 = std::sqrt(30375.0 / 112.0);
    const auto fs_45_28_35 = std::sqrt(10125.0 / 112.0);
    const auto fs_45_28_42 = std::sqrt(6075.0 / 56.0);
    const auto fs_45_28_5 = std::sqrt(10125.0 / 784.0);
    const auto fs_45_28_7 = std::sqrt(2025.0 / 112.0);
    const auto fs_45_56_14 = std::sqrt(2025.0 / 224.0);
    const auto fs_45_56_2 = std::sqrt(2025.0 / 1568.0);
    const auto fs_45_56_210 = std::sqrt(30375.0 / 224.0);
    const auto fs_45_56_30 = std::sqrt(30375.0 / 1568.0);
    const auto fs_45_77_7 = std::sqrt(2025.0 / 847.0);
    const auto fs_45_7_5 = std::sqrt(10125.0 / 49.0);
    const auto fs_45_8_10 = std::sqrt(10125.0 / 32.0);
    const auto fs_45_8_7 = std::sqrt(14175.0 / 64.0);
    const auto fs_48_1001_7 = std::sqrt(2304.0 / 143143.0);
    const auto fs_48_77_7 = std::sqrt(2304.0 / 847.0);
    const auto fs_4_143_55 = std::sqrt(80.0 / 1859.0);
    const auto fs_4_33_2 = std::sqrt(32.0 / 1089.0);
    const auto fs_57_1001_6 = std::sqrt(19494.0 / 1002001.0);
    const auto fs_57_28_6 = std::sqrt(9747.0 / 392.0);
    const auto fs_57_77_6 = std::sqrt(19494.0 / 5929.0);
    const auto fs_5_11_2 = std::sqrt(50.0 / 121.0);
    const auto fs_5_11_21 = std::sqrt(525.0 / 121.0);
    const auto fs_5_11_7 = std::sqrt(175.0 / 121.0);
    const auto fs_5_143_66 = std::sqrt(150.0 / 1859.0);
    const auto fs_5_14_105 = std::sqrt(375.0 / 28.0);
    const auto fs_5_14_35 = std::sqrt(125.0 / 28.0);
    const auto fs_5_14_42 = std::sqrt(75.0 / 14.0);
    const auto fs_5_14_5 = std::sqrt(125.0 / 196.0);
    const auto fs_5_22_14 = std::sqrt(175.0 / 242.0);
    const auto fs_5_22_15 = std::sqrt(375.0 / 484.0);
    const auto fs_5_22_21 = std::sqrt(525.0 / 484.0);
    const auto fs_5_22_3 = std::sqrt(75.0 / 484.0);
    const auto fs_5_22_33 = std::sqrt(75.0 / 44.0);
    const auto fs_5_22_35 = std::sqrt(875.0 / 484.0);
    const auto fs_5_22_42 = std::sqrt(525.0 / 242.0);
    const auto fs_5_22_6 = std::sqrt(75.0 / 242.0);
    const auto fs_5_231_105 = std::sqrt(125.0 / 2541.0);
    const auto fs_5_231_35 = std::sqrt(125.0 / 7623.0);
    const auto fs_5_231_42 = std::sqrt(50.0 / 2541.0);
    const auto fs_5_231_5 = std::sqrt(125.0 / 53361.0);
    const auto fs_5_286_22 = std::sqrt(25.0 / 3718.0);
    const auto fs_5_28_14 = std::sqrt(25.0 / 56.0);
    const auto fs_5_28_2 = std::sqrt(25.0 / 392.0);
    const auto fs_5_28_210 = std::sqrt(375.0 / 56.0);
    const auto fs_5_28_30 = std::sqrt(375.0 / 392.0);
    const auto fs_5_44_30 = std::sqrt(375.0 / 968.0);
    const auto fs_5_44_6 = std::sqrt(75.0 / 968.0);
    const auto fs_5_44_66 = std::sqrt(75.0 / 88.0);
    const auto fs_5_462_14 = std::sqrt(25.0 / 15246.0);
    const auto fs_5_462_2 = std::sqrt(25.0 / 106722.0);
    const auto fs_5_462_210 = std::sqrt(125.0 / 5082.0);
    const auto fs_5_462_30 = std::sqrt(125.0 / 35574.0);
    const auto fs_5_77_10 = std::sqrt(250.0 / 5929.0);
    const auto fs_5_77_7 = std::sqrt(25.0 / 847.0);
    const auto fs_5_7_14 = std::sqrt(50.0 / 7.0);
    const auto fs_5_7_3 = std::sqrt(75.0 / 49.0);
    const auto fs_5_7_7 = std::sqrt(25.0 / 7.0);
    const auto fs_6_1001_21 = std::sqrt(108.0 / 143143.0);
    const auto fs_6_11_15 = std::sqrt(540.0 / 121.0);
    const auto fs_6_11_5 = std::sqrt(180.0 / 121.0);
    const auto fs_6_143_15 = std::sqrt(540.0 / 20449.0);
    const auto fs_6_143_42 = std::sqrt(1512.0 / 20449.0);
    const auto fs_6_143_5 = std::sqrt(180.0 / 20449.0);
    const auto fs_6_77_21 = std::sqrt(108.0 / 847.0);
    const auto fs_75_154_10 = std::sqrt(28125.0 / 11858.0);
    const auto fs_75_2002_10 = std::sqrt(28125.0 / 2004002.0);
    const auto fs_75_56_10 = std::sqrt(28125.0 / 1568.0);
    const auto fs_75_8_2 = std::sqrt(5625.0 / 32.0);
    const auto fs_7_66_2 = std::sqrt(49.0 / 2178.0);
    const auto fs_8_143_7 = std::sqrt(448.0 / 20449.0);
    const auto fs_90_1001_7 = std::sqrt(8100.0 / 143143.0);
    const auto fs_90_77_7 = std::sqrt(8100.0 / 847.0);

    // NOTE: the angular components are formed in 36 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph2_p1, ph2_p2, ph4_p1, ph4_p2, ph6_p1, ph6_p2, ph8_p1, ph8_p2, ph8_p7, ph8_p8, ab_2, pc_0, pc_1 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h2_p1 = ph2_p1[k];
        const auto h2_p2 = ph2_p2[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h8_p2 = ph8_p2[k];
        const auto h8_p7 = ph8_p7[k];
        const auto h8_p8 = ph8_p8[k];

        pc_0[k] = e_0 * fs_15_8_105 * h2_p2 + e_1 * fs_45_28_7 * h4_p2 - e_1 * fs_45_28_105 * r_2 * h2_p2 + e_2 * fs_5_44_6 * h6_p2 - e_2 * fs_45_77_7 * r_2 * h4_p2 + e_2 * fs_5_14_105 * r_4 * h2_p2 + e_3 * fs_1_286_2 * h8_p2 + e_3 * fs_2_143_1001 * h8_p8 - e_3 * fs_1_66_6 * r_2 * h6_p2 + e_3 * fs_45_1001_7 * r_4 * h4_p2 - e_3 * fs_5_231_105 * r_6 * h2_p2;

        pc_1[k] = e_0 * fs_15_8_42 * h2_p1 + e_1 * fs_27_28_35 * h4_p1 - e_1 * fs_45_28_42 * r_2 * h2_p1 + e_2 * fs_5_22_6 * h6_p1 - e_2 * fs_27_77_35 * r_2 * h4_p1 + e_2 * fs_5_14_42 * r_4 * h2_p1 + e_3 * fs_1_286_14 * h8_p1 + e_3 * fs_1_286_10010 * h8_p7 - e_3 * fs_1_33_6 * r_2 * h6_p1 + e_3 * fs_27_1001_35 * r_4 * h4_p1 - e_3 * fs_5_231_42 * r_6 * h2_p1;
    }

    // NOTE: the angular components are formed in 36 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph2_0, ph2_p1, ph4_0, ph4_p1, ph6_0, ph6_p1, ph6_p5, ph6_p6, ph8_0, ph8_p1, ph8_p5, ph8_p6, ab_2, pc_2, pc_3 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

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

        pc_2[k] = e_0 * fs_15_4_7 * h2_0 + e_1 * fs_45_14_7 * h4_0 - e_1 * fs_45_14_7 * r_2 * h2_0 + e_2 * fs_5_11_7 * h6_0 + e_2 * fs_5_44_66 * h6_p6 - e_2 * fs_90_77_7 * r_2 * h4_0 + e_2 * fs_5_7_7 * r_4 * h2_0 + e_3 * fs_2_143_7 * h8_0 + e_3 * fs_1_286_6006 * h8_p6 - e_3 * fs_2_33_7 * r_2 * h6_0 - e_3 * fs_1_66_66 * r_2 * h6_p6 + e_3 * fs_90_1001_7 * r_4 * h4_0 - e_3 * fs_10_231_7 * r_6 * h2_0;

        pc_3[k] = - e_0 * fs_15_16_14 * h2_p1 - e_1 * fs_15_28_105 * h4_p1 + e_1 * fs_45_56_14 * r_2 * h2_p1 - e_2 * fs_35_44_2 * h6_p1 + e_2 * fs_5_22_33 * h6_p5 + e_2 * fs_15_77_105 * r_2 * h4_p1 - e_2 * fs_5_28_14 * r_4 * h2_p1 - e_3 * fs_1_143_42 * h8_p1 + e_3 * fs_1_143_858 * h8_p5 + e_3 * fs_7_66_2 * r_2 * h6_p1 - e_3 * fs_1_33_33 * r_2 * h6_p5 - e_3 * fs_15_1001_105 * r_4 * h4_p1 + e_3 * fs_5_462_14 * r_6 * h2_p1;
    }

    // NOTE: the angular components are formed in 36 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph2_p2, ph4_m3, ph4_p2, ph4_p4, ph6_m3, ph6_p2, ph6_p4, ph8_m3, ph8_p2, ph8_p4, ab_2, pc_4, pc_5 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

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

        pc_4[k] = e_0 * fs_15_16_2 * h2_p2 + e_1 * fs_45_56_30 * h4_p2 + e_1 * fs_3_28_210 * h4_p4 - e_1 * fs_45_56_2 * r_2 * h2_p2 + e_2 * fs_5_22_35 * h6_p2 + e_2 * fs_5_22_42 * h6_p4 - e_2 * fs_45_154_30 * r_2 * h4_p2 - e_2 * fs_3_77_210 * r_2 * h4_p4 + e_2 * fs_5_28_2 * r_4 * h2_p2 + e_3 * fs_1_143_105 * h8_p2 + e_3 * fs_1_143_462 * h8_p4 - e_3 * fs_1_33_35 * r_2 * h6_p2 - e_3 * fs_1_33_42 * r_2 * h6_p4 + e_3 * fs_45_2002_30 * r_4 * h4_p2 + e_3 * fs_3_1001_210 * r_4 * h4_p4 - e_3 * fs_5_462_2 * r_6 * h2_p2;

        pc_5[k] = - e_1 * fs_45_28_7 * h4_m3 - e_2 * fs_5_11_21 * h6_m3 + e_2 * fs_45_77_7 * r_2 * h4_m3 - e_3 * fs_1_143_462 * h8_m3 + e_3 * fs_2_33_21 * r_2 * h6_m3 - e_3 * fs_45_1001_7 * r_4 * h4_m3;
    }

    // NOTE: the angular components are formed in 36 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph2_m2, ph2_m1, ph4_m4, ph4_m2, ph4_m1, ph6_m5, ph6_m4, ph6_m2, ph6_m1, ph8_m5, ph8_m4, ph8_m2, ph8_m1, ab_2, pc_6, pc_7 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

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

        pc_6[k] = e_0 * fs_15_16_2 * h2_m2 - e_1 * fs_3_28_210 * h4_m4 + e_1 * fs_45_56_30 * h4_m2 - e_1 * fs_45_56_2 * r_2 * h2_m2 - e_2 * fs_5_22_42 * h6_m4 + e_2 * fs_5_22_35 * h6_m2 + e_2 * fs_3_77_210 * r_2 * h4_m4 - e_2 * fs_45_154_30 * r_2 * h4_m2 + e_2 * fs_5_28_2 * r_4 * h2_m2 - e_3 * fs_1_143_462 * h8_m4 + e_3 * fs_1_143_105 * h8_m2 + e_3 * fs_1_33_42 * r_2 * h6_m4 - e_3 * fs_1_33_35 * r_2 * h6_m2 - e_3 * fs_3_1001_210 * r_4 * h4_m4 + e_3 * fs_45_2002_30 * r_4 * h4_m2 - e_3 * fs_5_462_2 * r_6 * h2_m2;

        pc_7[k] = - e_0 * fs_15_16_14 * h2_m1 - e_1 * fs_15_28_105 * h4_m1 + e_1 * fs_45_56_14 * r_2 * h2_m1 - e_2 * fs_5_22_33 * h6_m5 - e_2 * fs_35_44_2 * h6_m1 + e_2 * fs_15_77_105 * r_2 * h4_m1 - e_2 * fs_5_28_14 * r_4 * h2_m1 - e_3 * fs_1_143_858 * h8_m5 - e_3 * fs_1_143_42 * h8_m1 + e_3 * fs_1_33_33 * r_2 * h6_m5 + e_3 * fs_7_66_2 * r_2 * h6_m1 - e_3 * fs_15_1001_105 * r_4 * h4_m1 + e_3 * fs_5_462_14 * r_6 * h2_m1;
    }

    // NOTE: the angular components are formed in 36 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph2_m2, ph2_m1, ph4_m2, ph4_m1, ph6_m6, ph6_m2, ph6_m1, ph8_m8, ph8_m7, ph8_m6, ph8_m2, ph8_m1, ab_2, pc_8, pc_9, pc_10 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

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

        pc_8[k] = - e_2 * fs_5_44_66 * h6_m6 - e_3 * fs_1_286_6006 * h8_m6 + e_3 * fs_1_66_66 * r_2 * h6_m6;

        pc_9[k] = - e_0 * fs_15_8_42 * h2_m1 - e_1 * fs_27_28_35 * h4_m1 + e_1 * fs_45_28_42 * r_2 * h2_m1 - e_2 * fs_5_22_6 * h6_m1 + e_2 * fs_27_77_35 * r_2 * h4_m1 - e_2 * fs_5_14_42 * r_4 * h2_m1 - e_3 * fs_1_286_10010 * h8_m7 - e_3 * fs_1_286_14 * h8_m1 + e_3 * fs_1_33_6 * r_2 * h6_m1 - e_3 * fs_27_1001_35 * r_4 * h4_m1 + e_3 * fs_5_231_42 * r_6 * h2_m1;

        pc_10[k] = - e_0 * fs_15_8_105 * h2_m2 - e_1 * fs_45_28_7 * h4_m2 + e_1 * fs_45_28_105 * r_2 * h2_m2 - e_2 * fs_5_44_6 * h6_m2 + e_2 * fs_45_77_7 * r_2 * h4_m2 - e_2 * fs_5_14_105 * r_4 * h2_m2 - e_3 * fs_2_143_1001 * h8_m8 - e_3 * fs_1_286_2 * h8_m2 + e_3 * fs_1_66_6 * r_2 * h6_m2 - e_3 * fs_45_1001_7 * r_4 * h4_m2 + e_3 * fs_5_231_105 * r_6 * h2_m2;
    }

    // NOTE: the angular components are formed in 36 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph2_p2, ph4_p2, ph4_p3, ph6_p2, ph6_p3, ph6_p6, ph8_p2, ph8_p3, ph8_p6, ph8_p7, ab_2, pc_11, pc_12 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h2_p2 = ph2_p2[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h4_p3 = ph4_p3[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h6_p6 = ph6_p6[k];
        const auto h8_p2 = ph8_p2[k];
        const auto h8_p3 = ph8_p3[k];
        const auto h8_p6 = ph8_p6[k];
        const auto h8_p7 = ph8_p7[k];

        pc_11[k] = - e_1 * fs_15_4_3 * h4_p3 - e_2 * f_15_22 * h6_p3 + e_2 * fs_15_11_3 * r_2 * h4_p3 - e_3 * fs_1_286_22 * h8_p3 + e_3 * fs_1_286_6006 * h8_p7 + e_3 * f_1_11 * r_2 * h6_p3 - e_3 * fs_15_143_3 * r_4 * h4_p3;

        pc_12[k] = e_0 * fs_45_8_7 * h2_p2 - e_1 * fs_3_7_105 * h4_p2 - e_1 * fs_135_28_7 * r_2 * h2_p2 - e_2 * fs_15_44_10 * h6_p2 - e_2 * fs_15_44_22 * h6_p6 + e_2 * fs_12_77_105 * r_2 * h4_p2 + e_2 * fs_15_14_7 * r_4 * h2_p2 - e_3 * fs_1_143_30 * h8_p2 + e_3 * fs_1_143_2002 * h8_p6 + e_3 * fs_1_22_10 * r_2 * h6_p2 + e_3 * fs_1_22_22 * r_2 * h6_p6 - e_3 * fs_12_1001_105 * r_4 * h4_p2 - e_3 * fs_5_77_7 * r_6 * h2_p2;
    }

    // NOTE: the angular components are formed in 36 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph2_0, ph2_p1, ph4_0, ph4_p1, ph4_p4, ph6_0, ph6_p1, ph6_p4, ph6_p5, ph8_0, ph8_p1, ph8_p4, ph8_p5, ab_2, pc_13, pc_14 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h2_0 = ph2_0[k];
        const auto h2_p1 = ph2_p1[k];
        const auto h4_0 = ph4_0[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h4_p4 = ph4_p4[k];
        const auto h6_0 = ph6_0[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h6_p4 = ph6_p4[k];
        const auto h6_p5 = ph6_p5[k];
        const auto h8_0 = ph8_0[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h8_p4 = ph8_p4[k];
        const auto h8_p5 = ph8_p5[k];

        pc_13[k] = e_0 * fs_15_4_14 * h2_p1 - e_1 * fs_3_28_105 * h4_p1 - e_1 * fs_45_14_14 * r_2 * h2_p1 - e_2 * fs_10_11_2 * h6_p1 - e_2 * fs_5_22_33 * h6_p5 + e_2 * fs_3_77_105 * r_2 * h4_p1 + e_2 * fs_5_7_14 * r_4 * h2_p1 - e_3 * fs_3_286_42 * h8_p1 + e_3 * fs_3_286_858 * h8_p5 + e_3 * fs_4_33_2 * r_2 * h6_p1 + e_3 * fs_1_33_33 * r_2 * h6_p5 - e_3 * fs_3_1001_105 * r_4 * h4_p1 - e_3 * fs_10_231_14 * r_6 * h2_p1;

        pc_14[k] = e_0 * fs_45_8_7 * h2_0 + e_1 * fs_15_14_7 * h4_0 - e_1 * fs_3_2_5 * h4_p4 - e_1 * fs_135_28_7 * r_2 * h2_0 - e_2 * fs_15_22_7 * h6_0 - e_2 * f_15_22 * h6_p4 - e_2 * fs_30_77_7 * r_2 * h4_0 + e_2 * fs_6_11_5 * r_2 * h4_p4 + e_2 * fs_15_14_7 * r_4 * h2_0 - e_3 * fs_8_143_7 * h8_0 + e_3 * fs_12_143_11 * h8_p4 + e_3 * fs_1_11_7 * r_2 * h6_0 + e_3 * f_1_11 * r_2 * h6_p4 + e_3 * fs_30_1001_7 * r_4 * h4_0 - e_3 * fs_6_143_5 * r_4 * h4_p4 - e_3 * fs_5_77_7 * r_6 * h2_0;
    }

    // NOTE: the angular components are formed in 36 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph2_m2, ph2_p1, ph4_m2, ph4_p1, ph4_p3, ph6_m2, ph6_p1, ph8_m2, ph8_p1, ph8_p3, ab_2, pc_15, pc_16 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h2_m2 = ph2_m2[k];
        const auto h2_p1 = ph2_p1[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h4_p3 = ph4_p3[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h8_m2 = ph8_m2[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h8_p3 = ph8_p3[k];

        pc_15[k] = - e_0 * fs_15_4_3 * h2_p1 - e_1 * fs_75_56_10 * h4_p1 - e_1 * fs_33_56_70 * h4_p3 + e_1 * fs_45_14_3 * r_2 * h2_p1 + e_2 * fs_5_22_21 * h6_p1 + e_2 * fs_75_154_10 * r_2 * h4_p1 + e_2 * fs_3_14_70 * r_2 * h4_p3 - e_2 * fs_5_7_3 * r_4 * h2_p1 + e_3 * f_21_143 * h8_p1 + e_3 * fs_1_143_1155 * h8_p3 - e_3 * fs_1_33_21 * r_2 * h6_p1 - e_3 * fs_75_2002_10 * r_4 * h4_p1 - e_3 * fs_3_182_70 * r_4 * h4_p3 + e_3 * fs_10_231_3 * r_6 * h2_p1;

        pc_16[k] = e_0 * fs_15_8_5 * h2_m2 + e_1 * fs_30_7_3 * h4_m2 - e_1 * fs_45_28_5 * r_2 * h2_m2 - e_2 * fs_5_22_14 * h6_m2 - e_2 * fs_120_77_3 * r_2 * h4_m2 + e_2 * fs_5_14_5 * r_4 * h2_m2 - e_3 * fs_6_143_42 * h8_m2 + e_3 * fs_1_33_14 * r_2 * h6_m2 + e_3 * fs_120_1001_3 * r_4 * h4_m2 - e_3 * fs_5_231_5 * r_6 * h2_m2;
    }

    // NOTE: the angular components are formed in 36 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph2_m1, ph4_m4, ph4_m3, ph4_m1, ph6_m5, ph6_m4, ph6_m1, ph8_m5, ph8_m4, ph8_m3, ph8_m1, ab_2, pc_17, pc_18, pc_19 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h2_m1 = ph2_m1[k];
        const auto h4_m4 = ph4_m4[k];
        const auto h4_m3 = ph4_m3[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h6_m5 = ph6_m5[k];
        const auto h6_m4 = ph6_m4[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h8_m5 = ph8_m5[k];
        const auto h8_m4 = ph8_m4[k];
        const auto h8_m3 = ph8_m3[k];
        const auto h8_m1 = ph8_m1[k];

        pc_17[k] = - e_0 * fs_15_4_3 * h2_m1 + e_1 * fs_33_56_70 * h4_m3 - e_1 * fs_75_56_10 * h4_m1 + e_1 * fs_45_14_3 * r_2 * h2_m1 + e_2 * fs_5_22_21 * h6_m1 - e_2 * fs_3_14_70 * r_2 * h4_m3 + e_2 * fs_75_154_10 * r_2 * h4_m1 - e_2 * fs_5_7_3 * r_4 * h2_m1 - e_3 * fs_1_143_1155 * h8_m3 + e_3 * f_21_143 * h8_m1 - e_3 * fs_1_33_21 * r_2 * h6_m1 + e_3 * fs_3_182_70 * r_4 * h4_m3 - e_3 * fs_75_2002_10 * r_4 * h4_m1 + e_3 * fs_10_231_3 * r_6 * h2_m1;

        pc_18[k] = e_1 * fs_3_2_5 * h4_m4 + e_2 * f_15_22 * h6_m4 - e_2 * fs_6_11_5 * r_2 * h4_m4 - e_3 * fs_12_143_11 * h8_m4 - e_3 * f_1_11 * r_2 * h6_m4 + e_3 * fs_6_143_5 * r_4 * h4_m4;

        pc_19[k] = - e_0 * fs_15_4_14 * h2_m1 + e_1 * fs_3_28_105 * h4_m1 + e_1 * fs_45_14_14 * r_2 * h2_m1 + e_2 * fs_5_22_33 * h6_m5 + e_2 * fs_10_11_2 * h6_m1 - e_2 * fs_3_77_105 * r_2 * h4_m1 - e_2 * fs_5_7_14 * r_4 * h2_m1 - e_3 * fs_3_286_858 * h8_m5 + e_3 * fs_3_286_42 * h8_m1 - e_3 * fs_1_33_33 * r_2 * h6_m5 - e_3 * fs_4_33_2 * r_2 * h6_m1 + e_3 * fs_3_1001_105 * r_4 * h4_m1 + e_3 * fs_10_231_14 * r_6 * h2_m1;
    }

    // NOTE: the angular components are formed in 36 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph2_m2, ph4_m3, ph4_m2, ph6_m6, ph6_m3, ph6_m2, ph8_m7, ph8_m6, ph8_m3, ph8_m2, ab_2, pc_20, pc_21 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h2_m2 = ph2_m2[k];
        const auto h4_m3 = ph4_m3[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h6_m6 = ph6_m6[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h8_m7 = ph8_m7[k];
        const auto h8_m6 = ph8_m6[k];
        const auto h8_m3 = ph8_m3[k];
        const auto h8_m2 = ph8_m2[k];

        pc_20[k] = - e_0 * fs_45_8_7 * h2_m2 + e_1 * fs_3_7_105 * h4_m2 + e_1 * fs_135_28_7 * r_2 * h2_m2 + e_2 * fs_15_44_22 * h6_m6 + e_2 * fs_15_44_10 * h6_m2 - e_2 * fs_12_77_105 * r_2 * h4_m2 - e_2 * fs_15_14_7 * r_4 * h2_m2 - e_3 * fs_1_143_2002 * h8_m6 + e_3 * fs_1_143_30 * h8_m2 - e_3 * fs_1_22_22 * r_2 * h6_m6 - e_3 * fs_1_22_10 * r_2 * h6_m2 + e_3 * fs_12_1001_105 * r_4 * h4_m2 + e_3 * fs_5_77_7 * r_6 * h2_m2;

        pc_21[k] = e_1 * fs_15_4_3 * h4_m3 + e_2 * f_15_22 * h6_m3 - e_2 * fs_15_11_3 * r_2 * h4_m3 - e_3 * fs_1_286_6006 * h8_m7 + e_3 * fs_1_286_22 * h8_m3 - e_3 * f_1_11 * r_2 * h6_m3 + e_3 * fs_15_143_3 * r_4 * h4_m3;
    }

    // NOTE: the angular components are formed in 36 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, ph4_p3, ph4_p4, ph6_p3, ph6_p4, ph6_p6, ph8_p3, ph8_p4, ph8_p5, ph8_p6, ab_2, pc_22, pc_23 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
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
        const auto h6_p6 = ph6_p6[k];
        const auto h8_p3 = ph8_p3[k];
        const auto h8_p4 = ph8_p4[k];
        const auto h8_p5 = ph8_p5[k];
        const auto h8_p6 = ph8_p6[k];

        pc_22[k] = e_1 * fs_3_2_15 * h4_p4 + e_2 * fs_15_22_3 * h6_p4 + e_2 * fs_15_44_22 * h6_p6 - e_2 * fs_6_11_15 * r_2 * h4_p4 + e_3 * fs_1_143_33 * h8_p4 + e_3 * fs_1_286_2002 * h8_p6 - e_3 * fs_1_11_3 * r_2 * h6_p4 - e_3 * fs_1_22_22 * r_2 * h6_p6 + e_3 * fs_6_143_15 * r_4 * h4_p4;

        pc_23[k] = - e_1 * fs_3_4_3 * h4_p3 + e_2 * f_15_11 * h6_p3 + e_2 * fs_3_11_3 * r_2 * h4_p3 + e_3 * fs_5_286_22 * h8_p3 + e_3 * fs_1_286_4290 * h8_p5 - e_3 * f_2_11 * r_2 * h6_p3 - e_3 * fs_3_143_3 * r_4 * h4_p3;
    }

    // NOTE: the angular components are formed in 36 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph2_p2, ph4_p2, ph4_p4, ph6_p2, ph6_p4, ph8_p2, ph8_p4, ab_2, pc_24 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h2_p2 = ph2_p2[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h4_p4 = ph4_p4[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h6_p4 = ph6_p4[k];
        const auto h8_p2 = ph8_p2[k];
        const auto h8_p4 = ph8_p4[k];

        pc_24[k] = e_0 * fs_15_8_35 * h2_p2 - e_1 * fs_27_28_21 * h4_p2 + e_1 * fs_3_1_3 * h4_p4 - e_1 * fs_45_28_35 * r_2 * h2_p2 + e_2 * fs_35_44_2 * h6_p2 - e_2 * fs_5_22_15 * h6_p4 + e_2 * fs_27_77_21 * r_2 * h4_p2 - e_2 * fs_12_11_3 * r_2 * h4_p4 + e_2 * fs_5_14_35 * r_4 * h2_p2 + e_3 * fs_15_286_6 * h8_p2 + e_3 * fs_3_143_165 * h8_p4 - e_3 * fs_7_66_2 * r_2 * h6_p2 + e_3 * fs_1_33_15 * r_2 * h6_p4 - e_3 * fs_27_1001_21 * r_4 * h4_p2 + e_3 * fs_12_143_3 * r_4 * h4_p4 - e_3 * fs_5_231_35 * r_6 * h2_p2;
    }

    // NOTE: the angular components are formed in 36 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph2_p1, ph4_p1, ph4_p3, ph6_p1, ph6_p3, ph8_p1, ph8_p3, ab_2, pc_25 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h2_p1 = ph2_p1[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h4_p3 = ph4_p3[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h8_p3 = ph8_p3[k];

        pc_25[k] = e_0 * fs_15_16_210 * h2_p1 - e_1 * fs_12_7_7 * h4_p1 + e_1 * f_21_4 * h4_p3 - e_1 * fs_45_56_210 * r_2 * h2_p1 + e_2 * fs_5_44_30 * h6_p1 - e_2 * fs_15_22_3 * h6_p3 + e_2 * fs_48_77_7 * r_2 * h4_p1 - e_2 * f_21_11 * r_2 * h4_p3 + e_2 * fs_5_28_210 * r_4 * h2_p1 + e_3 * fs_3_143_70 * h8_p1 + e_3 * fs_5_143_66 * h8_p3 - e_3 * fs_1_66_30 * r_2 * h6_p1 + e_3 * fs_1_11_3 * r_2 * h6_p3 - e_3 * fs_48_1001_7 * r_4 * h4_p1 + e_3 * f_21_143 * r_4 * h4_p3 - e_3 * fs_5_462_210 * r_6 * h2_p1;
    }

    // NOTE: the angular components are formed in 36 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph2_m1, ph2_0, ph2_p2, ph4_m1, ph4_0, ph4_p2, ph6_m1, ph6_p2, ph8_m1, ph8_0, ph8_p2, ab_2, pc_26, pc_27 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h2_m1 = ph2_m1[k];
        const auto h2_0 = ph2_0[k];
        const auto h2_p2 = ph2_p2[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h4_0 = ph4_0[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h8_m1 = ph8_m1[k];
        const auto h8_0 = ph8_0[k];
        const auto h8_p2 = ph8_p2[k];

        pc_26[k] = e_0 * fs_45_8_10 * h2_0 + e_0 * fs_15_16_30 * h2_p2 - e_1 * fs_33_28_10 * h4_0 + e_1 * fs_129_56_2 * h4_p2 - e_1 * fs_135_28_10 * r_2 * h2_0 - e_1 * fs_45_56_30 * r_2 * h2_p2 - e_2 * fs_5_22_21 * h6_p2 + e_2 * fs_3_7_10 * r_2 * h4_0 - e_2 * fs_129_154_2 * r_2 * h4_p2 + e_2 * fs_15_14_10 * r_4 * h2_0 + e_2 * fs_5_28_30 * r_4 * h2_p2 + e_3 * fs_14_143_10 * h8_0 + e_3 * fs_15_143_7 * h8_p2 + e_3 * fs_1_33_21 * r_2 * h6_p2 - e_3 * fs_3_91_10 * r_4 * h4_0 + e_3 * fs_129_2002_2 * r_4 * h4_p2 - e_3 * fs_5_77_10 * r_6 * h2_0 - e_3 * fs_5_462_30 * r_6 * h2_p2;

        pc_27[k] = - e_0 * fs_75_8_2 * h2_m1 - e_1 * fs_3_28_15 * h4_m1 + e_1 * fs_225_28_2 * r_2 * h2_m1 + e_2 * fs_5_22_14 * h6_m1 + e_2 * fs_3_77_15 * r_2 * h4_m1 - e_2 * fs_25_14_2 * r_4 * h2_m1 - e_3 * fs_21_143_6 * h8_m1 - e_3 * fs_1_33_14 * r_2 * h6_m1 - e_3 * fs_3_1001_15 * r_4 * h4_m1 + e_3 * fs_25_231_2 * r_6 * h2_m1;
    }

    // NOTE: the angular components are formed in 36 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph2_m2, ph2_m1, ph4_m3, ph4_m2, ph4_m1, ph6_m3, ph6_m2, ph6_m1, ph8_m3, ph8_m2, ph8_m1, ab_2, pc_28, pc_29 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

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

        pc_28[k] = - e_0 * fs_15_16_30 * h2_m2 - e_1 * fs_129_56_2 * h4_m2 + e_1 * fs_45_56_30 * r_2 * h2_m2 + e_2 * fs_5_22_21 * h6_m2 + e_2 * fs_129_154_2 * r_2 * h4_m2 - e_2 * fs_5_28_30 * r_4 * h2_m2 - e_3 * fs_15_143_7 * h8_m2 - e_3 * fs_1_33_21 * r_2 * h6_m2 - e_3 * fs_129_2002_2 * r_4 * h4_m2 + e_3 * fs_5_462_30 * r_6 * h2_m2;

        pc_29[k] = - e_0 * fs_15_16_210 * h2_m1 - e_1 * f_21_4 * h4_m3 + e_1 * fs_12_7_7 * h4_m1 + e_1 * fs_45_56_210 * r_2 * h2_m1 + e_2 * fs_15_22_3 * h6_m3 - e_2 * fs_5_44_30 * h6_m1 + e_2 * f_21_11 * r_2 * h4_m3 - e_2 * fs_48_77_7 * r_2 * h4_m1 - e_2 * fs_5_28_210 * r_4 * h2_m1 - e_3 * fs_5_143_66 * h8_m3 - e_3 * fs_3_143_70 * h8_m1 - e_3 * fs_1_11_3 * r_2 * h6_m3 + e_3 * fs_1_66_30 * r_2 * h6_m1 - e_3 * f_21_143 * r_4 * h4_m3 + e_3 * fs_48_1001_7 * r_4 * h4_m1 + e_3 * fs_5_462_210 * r_6 * h2_m1;
    }

    // NOTE: the angular components are formed in 36 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph2_m2, ph4_m4, ph4_m3, ph4_m2, ph6_m4, ph6_m3, ph6_m2, ph8_m5, ph8_m4, ph8_m3, ph8_m2, ab_2, pc_30, pc_31 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h2_m2 = ph2_m2[k];
        const auto h4_m4 = ph4_m4[k];
        const auto h4_m3 = ph4_m3[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h6_m4 = ph6_m4[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h8_m5 = ph8_m5[k];
        const auto h8_m4 = ph8_m4[k];
        const auto h8_m3 = ph8_m3[k];
        const auto h8_m2 = ph8_m2[k];

        pc_30[k] = - e_0 * fs_15_8_35 * h2_m2 - e_1 * fs_3_1_3 * h4_m4 + e_1 * fs_27_28_21 * h4_m2 + e_1 * fs_45_28_35 * r_2 * h2_m2 + e_2 * fs_5_22_15 * h6_m4 - e_2 * fs_35_44_2 * h6_m2 + e_2 * fs_12_11_3 * r_2 * h4_m4 - e_2 * fs_27_77_21 * r_2 * h4_m2 - e_2 * fs_5_14_35 * r_4 * h2_m2 - e_3 * fs_3_143_165 * h8_m4 - e_3 * fs_15_286_6 * h8_m2 - e_3 * fs_1_33_15 * r_2 * h6_m4 + e_3 * fs_7_66_2 * r_2 * h6_m2 - e_3 * fs_12_143_3 * r_4 * h4_m4 + e_3 * fs_27_1001_21 * r_4 * h4_m2 + e_3 * fs_5_231_35 * r_6 * h2_m2;

        pc_31[k] = e_1 * fs_3_4_3 * h4_m3 - e_2 * f_15_11 * h6_m3 - e_2 * fs_3_11_3 * r_2 * h4_m3 - e_3 * fs_1_286_4290 * h8_m5 - e_3 * fs_5_286_22 * h8_m3 + e_3 * f_2_11 * r_2 * h6_m3 + e_3 * fs_3_143_3 * r_4 * h4_m3;
    }

    // NOTE: the angular components are formed in 36 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, ph4_m4, ph4_m3, ph6_m6, ph6_m5, ph6_m4, ph6_m3, ph8_m6, ph8_m5, ph8_m4, ph8_m3, ab_2, pc_32, pc_33, pc_34, pc_35 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

        const auto h4_m4 = ph4_m4[k];
        const auto h4_m3 = ph4_m3[k];
        const auto h6_m6 = ph6_m6[k];
        const auto h6_m5 = ph6_m5[k];
        const auto h6_m4 = ph6_m4[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h8_m6 = ph8_m6[k];
        const auto h8_m5 = ph8_m5[k];
        const auto h8_m4 = ph8_m4[k];
        const auto h8_m3 = ph8_m3[k];

        pc_32[k] = - e_1 * fs_3_2_15 * h4_m4 - e_2 * fs_15_44_22 * h6_m6 - e_2 * fs_15_22_3 * h6_m4 + e_2 * fs_6_11_15 * r_2 * h4_m4 - e_3 * fs_1_286_2002 * h8_m6 - e_3 * fs_1_143_33 * h8_m4 + e_3 * fs_1_22_22 * r_2 * h6_m6 + e_3 * fs_1_11_3 * r_2 * h6_m4 - e_3 * fs_6_143_15 * r_4 * h4_m4;

        pc_33[k] = - e_2 * fs_15_22_11 * h6_m5 - e_3 * fs_1_143_286 * h8_m5 + e_3 * fs_1_11_11 * r_2 * h6_m5;

        pc_34[k] = e_1 * f_9_1 * h4_m4 - e_2 * fs_15_22_5 * h6_m4 - e_2 * f_36_11 * r_2 * h4_m4 - e_3 * fs_4_143_55 * h8_m4 + e_3 * fs_1_11_5 * r_2 * h6_m4 + e_3 * f_36_143 * r_4 * h4_m4;

        pc_35[k] = e_1 * f_9_2 * h4_m3 - e_2 * fs_5_22_3 * h6_m3 - e_2 * f_18_11 * r_2 * h4_m3 - e_3 * fs_5_143_66 * h8_m3 + e_3 * fs_1_33_3 * r_2 * h6_m3 + e_3 * f_18_143 * r_4 * h4_m3;
    }

    // NOTE: the angular components are formed in 36 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph2_m2, ph2_m1, ph2_0, ph4_m2, ph4_m1, ph4_0, ph6_m2, ph6_m1, ph6_0, ph8_m2, ph8_m1, ph8_0, ab_2, pc_36, pc_37, pc_38 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h2_m2 = ph2_m2[k];
        const auto h2_m1 = ph2_m1[k];
        const auto h2_0 = ph2_0[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h4_0 = ph4_0[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h6_0 = ph6_0[k];
        const auto h8_m2 = ph8_m2[k];
        const auto h8_m1 = ph8_m1[k];
        const auto h8_0 = ph8_0[k];

        pc_36[k] = e_0 * fs_15_8_35 * h2_m2 - e_1 * fs_3_14_21 * h4_m2 - e_1 * fs_45_28_35 * r_2 * h2_m2 + e_2 * fs_5_11_2 * h6_m2 + e_2 * fs_6_77_21 * r_2 * h4_m2 + e_2 * fs_5_14_35 * r_4 * h2_m2 - e_3 * fs_20_143_6 * h8_m2 - e_3 * fs_2_33_2 * r_2 * h6_m2 - e_3 * fs_6_1001_21 * r_4 * h4_m2 - e_3 * fs_5_231_35 * r_6 * h2_m2;

        pc_37[k] = e_0 * fs_15_2_5 * h2_m1 - e_1 * fs_57_28_6 * h4_m1 - e_1 * fs_45_7_5 * r_2 * h2_m1 + e_2 * fs_5_22_35 * h6_m1 + e_2 * fs_57_77_6 * r_2 * h4_m1 + e_2 * fs_10_7_5 * r_4 * h2_m1 - e_3 * fs_14_143_15 * h8_m1 - e_3 * fs_1_33_35 * r_2 * h6_m1 - e_3 * fs_57_1001_6 * r_4 * h4_m1 - e_3 * fs_20_231_5 * r_6 * h2_m1;

        pc_38[k] = e_0 * f_75_4 * h2_0 - e_1 * f_45_7 * h4_0 - e_1 * f_225_14 * r_2 * h2_0 + e_2 * f_35_22 * h6_0 + e_2 * f_180_77 * r_2 * h4_0 + e_2 * f_25_7 * r_4 * h2_0 - e_3 * f_56_143 * h8_0 - e_3 * f_7_33 * r_2 * h6_0 - e_3 * f_180_1001 * r_4 * h4_0 - e_3 * f_50_231 * r_6 * h2_0;
    }

    // NOTE: the angular components are formed in 36 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph2_p1, ph2_p2, ph4_p1, ph4_p2, ph4_p3, ph6_p1, ph6_p2, ph6_p3, ph8_p1, ph8_p2, ph8_p3, ab_2, pc_39, pc_40, pc_41 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h2_p1 = ph2_p1[k];
        const auto h2_p2 = ph2_p2[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h4_p3 = ph4_p3[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h8_p2 = ph8_p2[k];
        const auto h8_p3 = ph8_p3[k];

        pc_39[k] = e_0 * fs_15_2_5 * h2_p1 - e_1 * fs_57_28_6 * h4_p1 - e_1 * fs_45_7_5 * r_2 * h2_p1 + e_2 * fs_5_22_35 * h6_p1 + e_2 * fs_57_77_6 * r_2 * h4_p1 + e_2 * fs_10_7_5 * r_4 * h2_p1 - e_3 * fs_14_143_15 * h8_p1 - e_3 * fs_1_33_35 * r_2 * h6_p1 - e_3 * fs_57_1001_6 * r_4 * h4_p1 - e_3 * fs_20_231_5 * r_6 * h2_p1;

        pc_40[k] = e_0 * fs_15_8_35 * h2_p2 - e_1 * fs_3_14_21 * h4_p2 - e_1 * fs_45_28_35 * r_2 * h2_p2 + e_2 * fs_5_11_2 * h6_p2 + e_2 * fs_6_77_21 * r_2 * h4_p2 + e_2 * fs_5_14_35 * r_4 * h2_p2 - e_3 * fs_20_143_6 * h8_p2 - e_3 * fs_2_33_2 * r_2 * h6_p2 - e_3 * fs_6_1001_21 * r_4 * h4_p2 - e_3 * fs_5_231_35 * r_6 * h2_p2;

        pc_41[k] = e_1 * f_9_2 * h4_p3 - e_2 * fs_5_22_3 * h6_p3 - e_2 * f_18_11 * r_2 * h4_p3 - e_3 * fs_5_143_66 * h8_p3 + e_3 * fs_1_33_3 * r_2 * h6_p3 + e_3 * f_18_143 * r_4 * h4_p3;
    }

    // NOTE: the angular components are formed in 36 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, ph4_m4, ph4_p4, ph6_m6, ph6_m4, ph6_p4, ph6_p5, ph8_m6, ph8_m4, ph8_p4, ph8_p5, ab_2, pc_42, pc_43, pc_44 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

        const auto h4_m4 = ph4_m4[k];
        const auto h4_p4 = ph4_p4[k];
        const auto h6_m6 = ph6_m6[k];
        const auto h6_m4 = ph6_m4[k];
        const auto h6_p4 = ph6_p4[k];
        const auto h6_p5 = ph6_p5[k];
        const auto h8_m6 = ph8_m6[k];
        const auto h8_m4 = ph8_m4[k];
        const auto h8_p4 = ph8_p4[k];
        const auto h8_p5 = ph8_p5[k];

        pc_42[k] = e_1 * f_9_1 * h4_p4 - e_2 * fs_15_22_5 * h6_p4 - e_2 * f_36_11 * r_2 * h4_p4 - e_3 * fs_4_143_55 * h8_p4 + e_3 * fs_1_11_5 * r_2 * h6_p4 + e_3 * f_36_143 * r_4 * h4_p4;

        pc_43[k] = - e_2 * fs_15_22_11 * h6_p5 - e_3 * fs_1_143_286 * h8_p5 + e_3 * fs_1_11_11 * r_2 * h6_p5;

        pc_44[k] = e_1 * fs_3_2_15 * h4_m4 - e_2 * fs_15_44_22 * h6_m6 + e_2 * fs_15_22_3 * h6_m4 - e_2 * fs_6_11_15 * r_2 * h4_m4 - e_3 * fs_1_286_2002 * h8_m6 + e_3 * fs_1_143_33 * h8_m4 + e_3 * fs_1_22_22 * r_2 * h6_m6 - e_3 * fs_1_11_3 * r_2 * h6_m4 + e_3 * fs_6_143_15 * r_4 * h4_m4;
    }

    // NOTE: the angular components are formed in 36 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph2_m2, ph4_m4, ph4_m3, ph4_m2, ph6_m4, ph6_m3, ph6_m2, ph8_m5, ph8_m4, ph8_m3, ph8_m2, ab_2, pc_45, pc_46 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h2_m2 = ph2_m2[k];
        const auto h4_m4 = ph4_m4[k];
        const auto h4_m3 = ph4_m3[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h6_m4 = ph6_m4[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h8_m5 = ph8_m5[k];
        const auto h8_m4 = ph8_m4[k];
        const auto h8_m3 = ph8_m3[k];
        const auto h8_m2 = ph8_m2[k];

        pc_45[k] = - e_1 * fs_3_4_3 * h4_m3 + e_2 * f_15_11 * h6_m3 + e_2 * fs_3_11_3 * r_2 * h4_m3 - e_3 * fs_1_286_4290 * h8_m5 + e_3 * fs_5_286_22 * h8_m3 - e_3 * f_2_11 * r_2 * h6_m3 - e_3 * fs_3_143_3 * r_4 * h4_m3;

        pc_46[k] = e_0 * fs_15_8_35 * h2_m2 - e_1 * fs_3_1_3 * h4_m4 - e_1 * fs_27_28_21 * h4_m2 - e_1 * fs_45_28_35 * r_2 * h2_m2 + e_2 * fs_5_22_15 * h6_m4 + e_2 * fs_35_44_2 * h6_m2 + e_2 * fs_12_11_3 * r_2 * h4_m4 + e_2 * fs_27_77_21 * r_2 * h4_m2 + e_2 * fs_5_14_35 * r_4 * h2_m2 - e_3 * fs_3_143_165 * h8_m4 + e_3 * fs_15_286_6 * h8_m2 - e_3 * fs_1_33_15 * r_2 * h6_m4 - e_3 * fs_7_66_2 * r_2 * h6_m2 - e_3 * fs_12_143_3 * r_4 * h4_m4 - e_3 * fs_27_1001_21 * r_4 * h4_m2 - e_3 * fs_5_231_35 * r_6 * h2_m2;
    }

    // NOTE: the angular components are formed in 36 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph2_m2, ph2_m1, ph4_m3, ph4_m2, ph4_m1, ph6_m3, ph6_m2, ph6_m1, ph8_m3, ph8_m2, ph8_m1, ab_2, pc_47, pc_48 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

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

        pc_47[k] = e_0 * fs_15_16_210 * h2_m1 - e_1 * f_21_4 * h4_m3 - e_1 * fs_12_7_7 * h4_m1 - e_1 * fs_45_56_210 * r_2 * h2_m1 + e_2 * fs_15_22_3 * h6_m3 + e_2 * fs_5_44_30 * h6_m1 + e_2 * f_21_11 * r_2 * h4_m3 + e_2 * fs_48_77_7 * r_2 * h4_m1 + e_2 * fs_5_28_210 * r_4 * h2_m1 - e_3 * fs_5_143_66 * h8_m3 + e_3 * fs_3_143_70 * h8_m1 - e_3 * fs_1_11_3 * r_2 * h6_m3 - e_3 * fs_1_66_30 * r_2 * h6_m1 - e_3 * f_21_143 * r_4 * h4_m3 - e_3 * fs_48_1001_7 * r_4 * h4_m1 - e_3 * fs_5_462_210 * r_6 * h2_m1;

        pc_48[k] = - e_0 * fs_15_16_30 * h2_m2 - e_1 * fs_129_56_2 * h4_m2 + e_1 * fs_45_56_30 * r_2 * h2_m2 + e_2 * fs_5_22_21 * h6_m2 + e_2 * fs_129_154_2 * r_2 * h4_m2 - e_2 * fs_5_28_30 * r_4 * h2_m2 - e_3 * fs_15_143_7 * h8_m2 - e_3 * fs_1_33_21 * r_2 * h6_m2 - e_3 * fs_129_2002_2 * r_4 * h4_m2 + e_3 * fs_5_462_30 * r_6 * h2_m2;
    }

    // NOTE: the angular components are formed in 36 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph2_0, ph2_p1, ph2_p2, ph4_0, ph4_p1, ph4_p2, ph6_p1, ph6_p2, ph8_0, ph8_p1, ph8_p2, ab_2, pc_49, pc_50 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h2_0 = ph2_0[k];
        const auto h2_p1 = ph2_p1[k];
        const auto h2_p2 = ph2_p2[k];
        const auto h4_0 = ph4_0[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h8_0 = ph8_0[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h8_p2 = ph8_p2[k];

        pc_49[k] = - e_0 * fs_75_8_2 * h2_p1 - e_1 * fs_3_28_15 * h4_p1 + e_1 * fs_225_28_2 * r_2 * h2_p1 + e_2 * fs_5_22_14 * h6_p1 + e_2 * fs_3_77_15 * r_2 * h4_p1 - e_2 * fs_25_14_2 * r_4 * h2_p1 - e_3 * fs_21_143_6 * h8_p1 - e_3 * fs_1_33_14 * r_2 * h6_p1 - e_3 * fs_3_1001_15 * r_4 * h4_p1 + e_3 * fs_25_231_2 * r_6 * h2_p1;

        pc_50[k] = e_0 * fs_45_8_10 * h2_0 - e_0 * fs_15_16_30 * h2_p2 - e_1 * fs_33_28_10 * h4_0 - e_1 * fs_129_56_2 * h4_p2 - e_1 * fs_135_28_10 * r_2 * h2_0 + e_1 * fs_45_56_30 * r_2 * h2_p2 + e_2 * fs_5_22_21 * h6_p2 + e_2 * fs_3_7_10 * r_2 * h4_0 + e_2 * fs_129_154_2 * r_2 * h4_p2 + e_2 * fs_15_14_10 * r_4 * h2_0 - e_2 * fs_5_28_30 * r_4 * h2_p2 + e_3 * fs_14_143_10 * h8_0 - e_3 * fs_15_143_7 * h8_p2 - e_3 * fs_1_33_21 * r_2 * h6_p2 - e_3 * fs_3_91_10 * r_4 * h4_0 - e_3 * fs_129_2002_2 * r_4 * h4_p2 - e_3 * fs_5_77_10 * r_6 * h2_0 + e_3 * fs_5_462_30 * r_6 * h2_p2;
    }

    // NOTE: the angular components are formed in 36 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph2_p1, ph4_p1, ph4_p3, ph6_p1, ph6_p3, ph8_p1, ph8_p3, ab_2, pc_51 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h2_p1 = ph2_p1[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h4_p3 = ph4_p3[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h8_p3 = ph8_p3[k];

        pc_51[k] = e_0 * fs_15_16_210 * h2_p1 - e_1 * fs_12_7_7 * h4_p1 - e_1 * f_21_4 * h4_p3 - e_1 * fs_45_56_210 * r_2 * h2_p1 + e_2 * fs_5_44_30 * h6_p1 + e_2 * fs_15_22_3 * h6_p3 + e_2 * fs_48_77_7 * r_2 * h4_p1 + e_2 * f_21_11 * r_2 * h4_p3 + e_2 * fs_5_28_210 * r_4 * h2_p1 + e_3 * fs_3_143_70 * h8_p1 - e_3 * fs_5_143_66 * h8_p3 - e_3 * fs_1_66_30 * r_2 * h6_p1 - e_3 * fs_1_11_3 * r_2 * h6_p3 - e_3 * fs_48_1001_7 * r_4 * h4_p1 - e_3 * f_21_143 * r_4 * h4_p3 - e_3 * fs_5_462_210 * r_6 * h2_p1;
    }

    // NOTE: the angular components are formed in 36 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph2_p2, ph4_p2, ph4_p3, ph4_p4, ph6_p2, ph6_p3, ph6_p4, ph8_p2, ph8_p3, ph8_p4, ph8_p5, ab_2, pc_52, pc_53 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h2_p2 = ph2_p2[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h4_p3 = ph4_p3[k];
        const auto h4_p4 = ph4_p4[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h6_p4 = ph6_p4[k];
        const auto h8_p2 = ph8_p2[k];
        const auto h8_p3 = ph8_p3[k];
        const auto h8_p4 = ph8_p4[k];
        const auto h8_p5 = ph8_p5[k];

        pc_52[k] = e_0 * fs_15_8_35 * h2_p2 - e_1 * fs_27_28_21 * h4_p2 - e_1 * fs_3_1_3 * h4_p4 - e_1 * fs_45_28_35 * r_2 * h2_p2 + e_2 * fs_35_44_2 * h6_p2 + e_2 * fs_5_22_15 * h6_p4 + e_2 * fs_27_77_21 * r_2 * h4_p2 + e_2 * fs_12_11_3 * r_2 * h4_p4 + e_2 * fs_5_14_35 * r_4 * h2_p2 + e_3 * fs_15_286_6 * h8_p2 - e_3 * fs_3_143_165 * h8_p4 - e_3 * fs_7_66_2 * r_2 * h6_p2 - e_3 * fs_1_33_15 * r_2 * h6_p4 - e_3 * fs_27_1001_21 * r_4 * h4_p2 - e_3 * fs_12_143_3 * r_4 * h4_p4 - e_3 * fs_5_231_35 * r_6 * h2_p2;

        pc_53[k] = - e_1 * fs_3_4_3 * h4_p3 + e_2 * f_15_11 * h6_p3 + e_2 * fs_3_11_3 * r_2 * h4_p3 + e_3 * fs_5_286_22 * h8_p3 - e_3 * fs_1_286_4290 * h8_p5 - e_3 * f_2_11 * r_2 * h6_p3 - e_3 * fs_3_143_3 * r_4 * h4_p3;
    }

    // NOTE: the angular components are formed in 36 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, ph4_m3, ph4_p4, ph6_m3, ph6_p4, ph6_p6, ph8_m7, ph8_m3, ph8_p4, ph8_p6, ab_2, pc_54, pc_55 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

        const auto h4_m3 = ph4_m3[k];
        const auto h4_p4 = ph4_p4[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h6_p4 = ph6_p4[k];
        const auto h6_p6 = ph6_p6[k];
        const auto h8_m7 = ph8_m7[k];
        const auto h8_m3 = ph8_m3[k];
        const auto h8_p4 = ph8_p4[k];
        const auto h8_p6 = ph8_p6[k];

        pc_54[k] = e_1 * fs_3_2_15 * h4_p4 + e_2 * fs_15_22_3 * h6_p4 - e_2 * fs_15_44_22 * h6_p6 - e_2 * fs_6_11_15 * r_2 * h4_p4 + e_3 * fs_1_143_33 * h8_p4 - e_3 * fs_1_286_2002 * h8_p6 - e_3 * fs_1_11_3 * r_2 * h6_p4 + e_3 * fs_1_22_22 * r_2 * h6_p6 + e_3 * fs_6_143_15 * r_4 * h4_p4;

        pc_55[k] = - e_1 * fs_15_4_3 * h4_m3 - e_2 * f_15_22 * h6_m3 + e_2 * fs_15_11_3 * r_2 * h4_m3 - e_3 * fs_1_286_6006 * h8_m7 - e_3 * fs_1_286_22 * h8_m3 + e_3 * f_1_11 * r_2 * h6_m3 - e_3 * fs_15_143_3 * r_4 * h4_m3;
    }

    // NOTE: the angular components are formed in 36 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph2_m2, ph2_m1, ph4_m2, ph4_m1, ph6_m6, ph6_m5, ph6_m2, ph6_m1, ph8_m6, ph8_m5, ph8_m2, ph8_m1, ab_2, pc_56, pc_57 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h2_m2 = ph2_m2[k];
        const auto h2_m1 = ph2_m1[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h6_m6 = ph6_m6[k];
        const auto h6_m5 = ph6_m5[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h8_m6 = ph8_m6[k];
        const auto h8_m5 = ph8_m5[k];
        const auto h8_m2 = ph8_m2[k];
        const auto h8_m1 = ph8_m1[k];

        pc_56[k] = e_0 * fs_45_8_7 * h2_m2 - e_1 * fs_3_7_105 * h4_m2 - e_1 * fs_135_28_7 * r_2 * h2_m2 + e_2 * fs_15_44_22 * h6_m6 - e_2 * fs_15_44_10 * h6_m2 + e_2 * fs_12_77_105 * r_2 * h4_m2 + e_2 * fs_15_14_7 * r_4 * h2_m2 - e_3 * fs_1_143_2002 * h8_m6 - e_3 * fs_1_143_30 * h8_m2 - e_3 * fs_1_22_22 * r_2 * h6_m6 + e_3 * fs_1_22_10 * r_2 * h6_m2 - e_3 * fs_12_1001_105 * r_4 * h4_m2 - e_3 * fs_5_77_7 * r_6 * h2_m2;

        pc_57[k] = e_0 * fs_15_4_14 * h2_m1 - e_1 * fs_3_28_105 * h4_m1 - e_1 * fs_45_14_14 * r_2 * h2_m1 + e_2 * fs_5_22_33 * h6_m5 - e_2 * fs_10_11_2 * h6_m1 + e_2 * fs_3_77_105 * r_2 * h4_m1 + e_2 * fs_5_7_14 * r_4 * h2_m1 - e_3 * fs_3_286_858 * h8_m5 - e_3 * fs_3_286_42 * h8_m1 - e_3 * fs_1_33_33 * r_2 * h6_m5 + e_3 * fs_4_33_2 * r_2 * h6_m1 - e_3 * fs_3_1001_105 * r_4 * h4_m1 - e_3 * fs_10_231_14 * r_6 * h2_m1;
    }

    // NOTE: the angular components are formed in 36 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph2_m1, ph4_m4, ph4_m3, ph4_m1, ph6_m4, ph6_m1, ph8_m4, ph8_m3, ph8_m1, ab_2, pc_58, pc_59 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h2_m1 = ph2_m1[k];
        const auto h4_m4 = ph4_m4[k];
        const auto h4_m3 = ph4_m3[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h6_m4 = ph6_m4[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h8_m4 = ph8_m4[k];
        const auto h8_m3 = ph8_m3[k];
        const auto h8_m1 = ph8_m1[k];

        pc_58[k] = e_1 * fs_3_2_5 * h4_m4 + e_2 * f_15_22 * h6_m4 - e_2 * fs_6_11_5 * r_2 * h4_m4 - e_3 * fs_12_143_11 * h8_m4 - e_3 * f_1_11 * r_2 * h6_m4 + e_3 * fs_6_143_5 * r_4 * h4_m4;

        pc_59[k] = e_0 * fs_15_4_3 * h2_m1 + e_1 * fs_33_56_70 * h4_m3 + e_1 * fs_75_56_10 * h4_m1 - e_1 * fs_45_14_3 * r_2 * h2_m1 - e_2 * fs_5_22_21 * h6_m1 - e_2 * fs_3_14_70 * r_2 * h4_m3 - e_2 * fs_75_154_10 * r_2 * h4_m1 + e_2 * fs_5_7_3 * r_4 * h2_m1 - e_3 * fs_1_143_1155 * h8_m3 - e_3 * f_21_143 * h8_m1 + e_3 * fs_1_33_21 * r_2 * h6_m1 + e_3 * fs_3_182_70 * r_4 * h4_m3 + e_3 * fs_75_2002_10 * r_4 * h4_m1 - e_3 * fs_10_231_3 * r_6 * h2_m1;
    }

    // NOTE: the angular components are formed in 36 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph2_p1, ph2_p2, ph4_p1, ph4_p2, ph4_p3, ph6_p1, ph6_p2, ph8_p1, ph8_p2, ph8_p3, ab_2, pc_60, pc_61 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h2_p1 = ph2_p1[k];
        const auto h2_p2 = ph2_p2[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h4_p3 = ph4_p3[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h8_p2 = ph8_p2[k];
        const auto h8_p3 = ph8_p3[k];

        pc_60[k] = e_0 * fs_15_8_5 * h2_p2 + e_1 * fs_30_7_3 * h4_p2 - e_1 * fs_45_28_5 * r_2 * h2_p2 - e_2 * fs_5_22_14 * h6_p2 - e_2 * fs_120_77_3 * r_2 * h4_p2 + e_2 * fs_5_14_5 * r_4 * h2_p2 - e_3 * fs_6_143_42 * h8_p2 + e_3 * fs_1_33_14 * r_2 * h6_p2 + e_3 * fs_120_1001_3 * r_4 * h4_p2 - e_3 * fs_5_231_5 * r_6 * h2_p2;

        pc_61[k] = - e_0 * fs_15_4_3 * h2_p1 - e_1 * fs_75_56_10 * h4_p1 + e_1 * fs_33_56_70 * h4_p3 + e_1 * fs_45_14_3 * r_2 * h2_p1 + e_2 * fs_5_22_21 * h6_p1 + e_2 * fs_75_154_10 * r_2 * h4_p1 - e_2 * fs_3_14_70 * r_2 * h4_p3 - e_2 * fs_5_7_3 * r_4 * h2_p1 + e_3 * f_21_143 * h8_p1 - e_3 * fs_1_143_1155 * h8_p3 - e_3 * fs_1_33_21 * r_2 * h6_p1 - e_3 * fs_75_2002_10 * r_4 * h4_p1 + e_3 * fs_3_182_70 * r_4 * h4_p3 + e_3 * fs_10_231_3 * r_6 * h2_p1;
    }

    // NOTE: the angular components are formed in 36 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph2_0, ph2_p1, ph4_0, ph4_p1, ph4_p4, ph6_0, ph6_p1, ph6_p4, ph6_p5, ph8_0, ph8_p1, ph8_p4, ph8_p5, ab_2, pc_62, pc_63 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h2_0 = ph2_0[k];
        const auto h2_p1 = ph2_p1[k];
        const auto h4_0 = ph4_0[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h4_p4 = ph4_p4[k];
        const auto h6_0 = ph6_0[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h6_p4 = ph6_p4[k];
        const auto h6_p5 = ph6_p5[k];
        const auto h8_0 = ph8_0[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h8_p4 = ph8_p4[k];
        const auto h8_p5 = ph8_p5[k];

        pc_62[k] = e_0 * fs_45_8_7 * h2_0 + e_1 * fs_15_14_7 * h4_0 + e_1 * fs_3_2_5 * h4_p4 - e_1 * fs_135_28_7 * r_2 * h2_0 - e_2 * fs_15_22_7 * h6_0 + e_2 * f_15_22 * h6_p4 - e_2 * fs_30_77_7 * r_2 * h4_0 - e_2 * fs_6_11_5 * r_2 * h4_p4 + e_2 * fs_15_14_7 * r_4 * h2_0 - e_3 * fs_8_143_7 * h8_0 - e_3 * fs_12_143_11 * h8_p4 + e_3 * fs_1_11_7 * r_2 * h6_0 - e_3 * f_1_11 * r_2 * h6_p4 + e_3 * fs_30_1001_7 * r_4 * h4_0 + e_3 * fs_6_143_5 * r_4 * h4_p4 - e_3 * fs_5_77_7 * r_6 * h2_0;

        pc_63[k] = e_0 * fs_15_4_14 * h2_p1 - e_1 * fs_3_28_105 * h4_p1 - e_1 * fs_45_14_14 * r_2 * h2_p1 - e_2 * fs_10_11_2 * h6_p1 + e_2 * fs_5_22_33 * h6_p5 + e_2 * fs_3_77_105 * r_2 * h4_p1 + e_2 * fs_5_7_14 * r_4 * h2_p1 - e_3 * fs_3_286_42 * h8_p1 - e_3 * fs_3_286_858 * h8_p5 + e_3 * fs_4_33_2 * r_2 * h6_p1 - e_3 * fs_1_33_33 * r_2 * h6_p5 - e_3 * fs_3_1001_105 * r_4 * h4_p1 - e_3 * fs_10_231_14 * r_6 * h2_p1;
    }

    // NOTE: the angular components are formed in 36 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph2_p2, ph4_p2, ph4_p3, ph6_p2, ph6_p3, ph6_p6, ph8_p2, ph8_p3, ph8_p6, ph8_p7, ab_2, pc_64, pc_65 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h2_p2 = ph2_p2[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h4_p3 = ph4_p3[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h6_p6 = ph6_p6[k];
        const auto h8_p2 = ph8_p2[k];
        const auto h8_p3 = ph8_p3[k];
        const auto h8_p6 = ph8_p6[k];
        const auto h8_p7 = ph8_p7[k];

        pc_64[k] = e_0 * fs_45_8_7 * h2_p2 - e_1 * fs_3_7_105 * h4_p2 - e_1 * fs_135_28_7 * r_2 * h2_p2 - e_2 * fs_15_44_10 * h6_p2 + e_2 * fs_15_44_22 * h6_p6 + e_2 * fs_12_77_105 * r_2 * h4_p2 + e_2 * fs_15_14_7 * r_4 * h2_p2 - e_3 * fs_1_143_30 * h8_p2 - e_3 * fs_1_143_2002 * h8_p6 + e_3 * fs_1_22_10 * r_2 * h6_p2 - e_3 * fs_1_22_22 * r_2 * h6_p6 - e_3 * fs_12_1001_105 * r_4 * h4_p2 - e_3 * fs_5_77_7 * r_6 * h2_p2;

        pc_65[k] = - e_1 * fs_15_4_3 * h4_p3 - e_2 * f_15_22 * h6_p3 + e_2 * fs_15_11_3 * r_2 * h4_p3 - e_3 * fs_1_286_22 * h8_p3 - e_3 * fs_1_286_6006 * h8_p7 + e_3 * f_1_11 * r_2 * h6_p3 - e_3 * fs_15_143_3 * r_4 * h4_p3;
    }

    // NOTE: the angular components are formed in 36 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph2_m2, ph2_m1, ph4_m2, ph4_m1, ph6_m6, ph6_m2, ph6_m1, ph8_m8, ph8_m7, ph8_m6, ph8_m2, ph8_m1, ab_2, pc_66, pc_67, pc_68 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

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

        pc_66[k] = e_0 * fs_15_8_105 * h2_m2 + e_1 * fs_45_28_7 * h4_m2 - e_1 * fs_45_28_105 * r_2 * h2_m2 + e_2 * fs_5_44_6 * h6_m2 - e_2 * fs_45_77_7 * r_2 * h4_m2 + e_2 * fs_5_14_105 * r_4 * h2_m2 - e_3 * fs_2_143_1001 * h8_m8 + e_3 * fs_1_286_2 * h8_m2 - e_3 * fs_1_66_6 * r_2 * h6_m2 + e_3 * fs_45_1001_7 * r_4 * h4_m2 - e_3 * fs_5_231_105 * r_6 * h2_m2;

        pc_67[k] = e_0 * fs_15_8_42 * h2_m1 + e_1 * fs_27_28_35 * h4_m1 - e_1 * fs_45_28_42 * r_2 * h2_m1 + e_2 * fs_5_22_6 * h6_m1 - e_2 * fs_27_77_35 * r_2 * h4_m1 + e_2 * fs_5_14_42 * r_4 * h2_m1 - e_3 * fs_1_286_10010 * h8_m7 + e_3 * fs_1_286_14 * h8_m1 - e_3 * fs_1_33_6 * r_2 * h6_m1 + e_3 * fs_27_1001_35 * r_4 * h4_m1 - e_3 * fs_5_231_42 * r_6 * h2_m1;

        pc_68[k] = - e_2 * fs_5_44_66 * h6_m6 - e_3 * fs_1_286_6006 * h8_m6 + e_3 * fs_1_66_66 * r_2 * h6_m6;
    }

    // NOTE: the angular components are formed in 36 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph2_m2, ph2_m1, ph4_m4, ph4_m2, ph4_m1, ph6_m5, ph6_m4, ph6_m2, ph6_m1, ph8_m5, ph8_m4, ph8_m2, ph8_m1, ab_2, pc_69, pc_70 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

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

        pc_69[k] = e_0 * fs_15_16_14 * h2_m1 + e_1 * fs_15_28_105 * h4_m1 - e_1 * fs_45_56_14 * r_2 * h2_m1 - e_2 * fs_5_22_33 * h6_m5 + e_2 * fs_35_44_2 * h6_m1 - e_2 * fs_15_77_105 * r_2 * h4_m1 + e_2 * fs_5_28_14 * r_4 * h2_m1 - e_3 * fs_1_143_858 * h8_m5 + e_3 * fs_1_143_42 * h8_m1 + e_3 * fs_1_33_33 * r_2 * h6_m5 - e_3 * fs_7_66_2 * r_2 * h6_m1 + e_3 * fs_15_1001_105 * r_4 * h4_m1 - e_3 * fs_5_462_14 * r_6 * h2_m1;

        pc_70[k] = - e_0 * fs_15_16_2 * h2_m2 - e_1 * fs_3_28_210 * h4_m4 - e_1 * fs_45_56_30 * h4_m2 + e_1 * fs_45_56_2 * r_2 * h2_m2 - e_2 * fs_5_22_42 * h6_m4 - e_2 * fs_5_22_35 * h6_m2 + e_2 * fs_3_77_210 * r_2 * h4_m4 + e_2 * fs_45_154_30 * r_2 * h4_m2 - e_2 * fs_5_28_2 * r_4 * h2_m2 - e_3 * fs_1_143_462 * h8_m4 - e_3 * fs_1_143_105 * h8_m2 + e_3 * fs_1_33_42 * r_2 * h6_m4 + e_3 * fs_1_33_35 * r_2 * h6_m2 - e_3 * fs_3_1001_210 * r_4 * h4_m4 - e_3 * fs_45_2002_30 * r_4 * h4_m2 + e_3 * fs_5_462_2 * r_6 * h2_m2;
    }

    // NOTE: the angular components are formed in 36 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph2_p2, ph4_p2, ph4_p3, ph4_p4, ph6_p2, ph6_p3, ph6_p4, ph8_p2, ph8_p3, ph8_p4, ab_2, pc_71, pc_72 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h2_p2 = ph2_p2[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h4_p3 = ph4_p3[k];
        const auto h4_p4 = ph4_p4[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h6_p4 = ph6_p4[k];
        const auto h8_p2 = ph8_p2[k];
        const auto h8_p3 = ph8_p3[k];
        const auto h8_p4 = ph8_p4[k];

        pc_71[k] = - e_1 * fs_45_28_7 * h4_p3 - e_2 * fs_5_11_21 * h6_p3 + e_2 * fs_45_77_7 * r_2 * h4_p3 - e_3 * fs_1_143_462 * h8_p3 + e_3 * fs_2_33_21 * r_2 * h6_p3 - e_3 * fs_45_1001_7 * r_4 * h4_p3;

        pc_72[k] = e_0 * fs_15_16_2 * h2_p2 + e_1 * fs_45_56_30 * h4_p2 - e_1 * fs_3_28_210 * h4_p4 - e_1 * fs_45_56_2 * r_2 * h2_p2 + e_2 * fs_5_22_35 * h6_p2 - e_2 * fs_5_22_42 * h6_p4 - e_2 * fs_45_154_30 * r_2 * h4_p2 + e_2 * fs_3_77_210 * r_2 * h4_p4 + e_2 * fs_5_28_2 * r_4 * h2_p2 + e_3 * fs_1_143_105 * h8_p2 - e_3 * fs_1_143_462 * h8_p4 - e_3 * fs_1_33_35 * r_2 * h6_p2 + e_3 * fs_1_33_42 * r_2 * h6_p4 + e_3 * fs_45_2002_30 * r_4 * h4_p2 - e_3 * fs_3_1001_210 * r_4 * h4_p4 - e_3 * fs_5_462_2 * r_6 * h2_p2;
    }

    // NOTE: the angular components are formed in 36 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph2_0, ph2_p1, ph4_0, ph4_p1, ph6_0, ph6_p1, ph6_p5, ph6_p6, ph8_0, ph8_p1, ph8_p5, ph8_p6, ab_2, pc_73, pc_74 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

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

        pc_73[k] = - e_0 * fs_15_16_14 * h2_p1 - e_1 * fs_15_28_105 * h4_p1 + e_1 * fs_45_56_14 * r_2 * h2_p1 - e_2 * fs_35_44_2 * h6_p1 - e_2 * fs_5_22_33 * h6_p5 + e_2 * fs_15_77_105 * r_2 * h4_p1 - e_2 * fs_5_28_14 * r_4 * h2_p1 - e_3 * fs_1_143_42 * h8_p1 - e_3 * fs_1_143_858 * h8_p5 + e_3 * fs_7_66_2 * r_2 * h6_p1 + e_3 * fs_1_33_33 * r_2 * h6_p5 - e_3 * fs_15_1001_105 * r_4 * h4_p1 + e_3 * fs_5_462_14 * r_6 * h2_p1;

        pc_74[k] = e_0 * fs_15_4_7 * h2_0 + e_1 * fs_45_14_7 * h4_0 - e_1 * fs_45_14_7 * r_2 * h2_0 + e_2 * fs_5_11_7 * h6_0 - e_2 * fs_5_44_66 * h6_p6 - e_2 * fs_90_77_7 * r_2 * h4_0 + e_2 * fs_5_7_7 * r_4 * h2_0 + e_3 * fs_2_143_7 * h8_0 - e_3 * fs_1_286_6006 * h8_p6 - e_3 * fs_2_33_7 * r_2 * h6_0 + e_3 * fs_1_66_66 * r_2 * h6_p6 + e_3 * fs_90_1001_7 * r_4 * h4_0 - e_3 * fs_10_231_7 * r_6 * h2_0;
    }

    // NOTE: the angular components are formed in 36 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, ph2_p1, ph2_p2, ph4_p1, ph4_p2, ph6_p1, ph6_p2, ph8_p1, ph8_p2, ph8_p7, ph8_p8, ab_2, pc_75, pc_76 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h2_p1 = ph2_p1[k];
        const auto h2_p2 = ph2_p2[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h8_p2 = ph8_p2[k];
        const auto h8_p7 = ph8_p7[k];
        const auto h8_p8 = ph8_p8[k];

        pc_75[k] = e_0 * fs_15_8_42 * h2_p1 + e_1 * fs_27_28_35 * h4_p1 - e_1 * fs_45_28_42 * r_2 * h2_p1 + e_2 * fs_5_22_6 * h6_p1 - e_2 * fs_27_77_35 * r_2 * h4_p1 + e_2 * fs_5_14_42 * r_4 * h2_p1 + e_3 * fs_1_286_14 * h8_p1 - e_3 * fs_1_286_10010 * h8_p7 - e_3 * fs_1_33_6 * r_2 * h6_p1 + e_3 * fs_27_1001_35 * r_4 * h4_p1 - e_3 * fs_5_231_42 * r_6 * h2_p1;

        pc_76[k] = e_0 * fs_15_8_105 * h2_p2 + e_1 * fs_45_28_7 * h4_p2 - e_1 * fs_45_28_105 * r_2 * h2_p2 + e_2 * fs_5_44_6 * h6_p2 - e_2 * fs_45_77_7 * r_2 * h4_p2 + e_2 * fs_5_14_105 * r_4 * h2_p2 + e_3 * fs_1_286_2 * h8_p2 - e_3 * fs_2_143_1001 * h8_p8 - e_3 * fs_1_66_6 * r_2 * h6_p2 + e_3 * fs_45_1001_7 * r_4 * h4_p2 - e_3 * fs_5_231_105 * r_6 * h2_p2;
    }

    // NOTE: the values of a combination of angular components are stored as one
    // row of nvalues columns, with the component on bra side running slowest. The
    // rows which the symmetry relates to an already formed one are copied from it.

    const size_t sources[77] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76};

    for (size_t n = 0; n < 77; n++)
    {
        if (sources[n] != n) std::copy(values + sources[n] * nvalues, values + (sources[n] + 1) * nvalues, values + n * nvalues);
    }
}

}  // namespace simdt2ceri
