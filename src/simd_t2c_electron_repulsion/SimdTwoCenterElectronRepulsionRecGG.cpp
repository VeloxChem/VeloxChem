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



#include "SimdTwoCenterElectronRepulsionRecGG.hpp"

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
compute_gg_electron_repulsion(double                         *values,
                              const size_t                    nvalues,
                              const CBasisFunction           &bra,
                              const CBasisFunction           &ket,
                              const std::vector<CSimdMatrix> &harmonics,
                              const CSimdMatrix              &coordinates) -> void
{
    if ((bra.get_angular_momentum() != 4) || (ket.get_angular_momentum() != 4))
    {
        errors::assertMsgCritical(
            false, std::string("SimdTwoCenterElectronRepulsionFunc.compute_gg_electron_repulsion: Basis functions must be of angular momenta four and four"));
    }

    if (harmonics.size() < 8)
    {
        errors::assertMsgCritical(
            false, std::string("SimdTwoCenterElectronRepulsionFunc.compute_gg_electron_repulsion: Harmonics must reach angular momentum 8"));
    }

    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("SimdTwoCenterElectronRepulsionFunc.compute_gg_electron_repulsion: Number of values exceeds number of atom pairs"));
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
    // orders 4 to 8 alone, and the orders below them are formed on the
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

            const auto ff_0 = fbase / (fexp * fexp * fexp * fexp);

            const auto ff_1 = fbase * fmu / (fexp * fexp * fexp * fexp);

            const auto ff_2 = fbase * fmu * fmu / (fexp * fexp * fexp * fexp);

            const auto ff_3 = fbase * fmu * fmu * fmu / (fexp * fexp * fexp * fexp);

            const auto ff_4 = fbase * fmu * fmu * fmu * fmu / (fexp * fexp * fexp * fexp);

            const auto *bv_0 = boys.data(5, i * nprim_b + j);

            const auto *bv_1 = boys.data(6, i * nprim_b + j);

            const auto *bv_2 = boys.data(7, i * nprim_b + j);

            const auto *bv_3 = boys.data(8, i * nprim_b + j);

            const auto *bv_4 = boys.data(9, i * nprim_b + j);

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

    auto *pc_0 = values + 0 * nvalues;
    auto *pc_1 = values + 1 * nvalues;
    auto *pc_2 = values + 2 * nvalues;
    auto *pc_3 = values + 3 * nvalues;
    auto *pc_4 = values + 4 * nvalues;
    auto *pc_5 = values + 5 * nvalues;
    auto *pc_6 = values + 6 * nvalues;
    auto *pc_7 = values + 7 * nvalues;
    auto *pc_8 = values + 8 * nvalues;
    auto *pc_9 = values + 10 * nvalues;
    auto *pc_10 = values + 11 * nvalues;
    auto *pc_11 = values + 12 * nvalues;
    auto *pc_12 = values + 13 * nvalues;
    auto *pc_13 = values + 14 * nvalues;
    auto *pc_14 = values + 15 * nvalues;
    auto *pc_15 = values + 16 * nvalues;
    auto *pc_16 = values + 17 * nvalues;
    auto *pc_17 = values + 20 * nvalues;
    auto *pc_18 = values + 21 * nvalues;
    auto *pc_19 = values + 22 * nvalues;
    auto *pc_20 = values + 23 * nvalues;
    auto *pc_21 = values + 24 * nvalues;
    auto *pc_22 = values + 25 * nvalues;
    auto *pc_23 = values + 26 * nvalues;
    auto *pc_24 = values + 30 * nvalues;
    auto *pc_25 = values + 31 * nvalues;
    auto *pc_26 = values + 32 * nvalues;
    auto *pc_27 = values + 33 * nvalues;
    auto *pc_28 = values + 34 * nvalues;
    auto *pc_29 = values + 35 * nvalues;
    auto *pc_30 = values + 40 * nvalues;
    auto *pc_31 = values + 41 * nvalues;
    auto *pc_32 = values + 42 * nvalues;
    auto *pc_33 = values + 43 * nvalues;
    auto *pc_34 = values + 44 * nvalues;
    auto *pc_35 = values + 50 * nvalues;
    auto *pc_36 = values + 51 * nvalues;
    auto *pc_37 = values + 52 * nvalues;
    auto *pc_38 = values + 53 * nvalues;
    auto *pc_39 = values + 60 * nvalues;
    auto *pc_40 = values + 61 * nvalues;
    auto *pc_41 = values + 62 * nvalues;
    auto *pc_42 = values + 70 * nvalues;
    auto *pc_43 = values + 71 * nvalues;
    auto *pc_44 = values + 80 * nvalues;

    // NOTE: the factors of the terms depend on the angular momenta alone, so they
    // are formed once for the whole matrix instead of once for every atom pair.

    const auto f_100_693 = 100.0 / 693.0;
    const auto f_105_16 = 105.0 / 16.0;
    const auto f_10_3 = 10.0 / 3.0;
    const auto f_10_33 = 10.0 / 33.0;
    const auto f_15_1 = 15.0;
    const auto f_15_4 = 15.0 / 4.0;
    const auto f_162_1001 = 162.0 / 1001.0;
    const auto f_162_77 = 162.0 / 77.0;
    const auto f_17_99 = 17.0 / 99.0;
    const auto f_18_11 = 18.0 / 11.0;
    const auto f_18_143 = 18.0 / 143.0;
    const auto f_196_1287 = 196.0 / 1287.0;
    const auto f_1_33 = 1.0 / 33.0;
    const auto f_1_9 = 1.0 / 9.0;
    const auto f_1_99 = 1.0 / 99.0;
    const auto f_20_21 = 20.0 / 21.0;
    const auto f_20_99 = 20.0 / 99.0;
    const auto f_21_2 = 21.0 / 2.0;
    const auto f_255_28 = 255.0 / 28.0;
    const auto f_25_2 = 25.0 / 2.0;
    const auto f_27_11 = 27.0 / 11.0;
    const auto f_27_143 = 27.0 / 143.0;
    const auto f_27_4 = 27.0 / 4.0;
    const auto f_2_1 = 2.0;
    const auto f_2_9 = 2.0 / 9.0;
    const auto f_30_7 = 30.0 / 7.0;
    const auto f_35_2 = 35.0 / 2.0;
    const auto f_35_8 = 35.0 / 8.0;
    const auto f_392_1287 = 392.0 / 1287.0;
    const auto f_40_693 = 40.0 / 693.0;
    const auto f_490_1287 = 490.0 / 1287.0;
    const auto f_4_99 = 4.0 / 99.0;
    const auto f_50_21 = 50.0 / 21.0;
    const auto f_50_33 = 50.0 / 33.0;
    const auto f_56_1287 = 56.0 / 1287.0;
    const auto f_5_1 = 5.0;
    const auto f_5_22 = 5.0 / 22.0;
    const auto f_5_3 = 5.0 / 3.0;
    const auto f_5_6 = 5.0 / 6.0;
    const auto f_5_66 = 5.0 / 66.0;
    const auto f_5_99 = 5.0 / 99.0;
    const auto f_75_7 = 75.0 / 7.0;
    const auto f_7_1287 = 7.0 / 1287.0;
    const auto f_81_1001 = 81.0 / 1001.0;
    const auto f_81_14 = 81.0 / 14.0;
    const auto f_81_28 = 81.0 / 28.0;
    const auto f_81_77 = 81.0 / 77.0;
    const auto f_85_42 = 85.0 / 42.0;
    const auto f_85_66 = 85.0 / 66.0;
    const auto f_85_693 = 85.0 / 693.0;
    const auto f_85_8 = 85.0 / 8.0;
    const auto f_99_28 = 99.0 / 28.0;
    const auto f_9_2 = 9.0 / 2.0;
    const auto f_9_7 = 9.0 / 7.0;
    const auto f_9_91 = 9.0 / 91.0;
    const auto fs_10_11_5 = std::sqrt(500.0 / 121.0);
    const auto fs_10_231_15 = std::sqrt(500.0 / 17787.0);
    const auto fs_10_33_15 = std::sqrt(500.0 / 363.0);
    const auto fs_10_693_21 = std::sqrt(100.0 / 22869.0);
    const auto fs_10_99_3 = std::sqrt(100.0 / 3267.0);
    const auto fs_135_56_6 = std::sqrt(54675.0 / 1568.0);
    const auto fs_13_198_6 = std::sqrt(169.0 / 6534.0);
    const auto fs_14_1287_858 = std::sqrt(392.0 / 3861.0);
    const auto fs_14_429_10 = std::sqrt(1960.0 / 184041.0);
    const auto fs_14_429_70 = std::sqrt(13720.0 / 184041.0);
    const auto fs_14_429_77 = std::sqrt(1372.0 / 16731.0);
    const auto fs_15_14_21 = std::sqrt(675.0 / 28.0);
    const auto fs_15_28_30 = std::sqrt(3375.0 / 392.0);
    const auto fs_15_28_6 = std::sqrt(675.0 / 392.0);
    const auto fs_15_4_15 = std::sqrt(3375.0 / 16.0);
    const auto fs_15_4_6 = std::sqrt(675.0 / 8.0);
    const auto fs_15_8_21 = std::sqrt(4725.0 / 64.0);
    const auto fs_1_33_11 = std::sqrt(1.0 / 99.0);
    const auto fs_1_33_30 = std::sqrt(10.0 / 363.0);
    const auto fs_1_66_42 = std::sqrt(7.0 / 726.0);
    const auto fs_1_99_105 = std::sqrt(35.0 / 3267.0);
    const auto fs_1_99_210 = std::sqrt(70.0 / 3267.0);
    const auto fs_1_99_42 = std::sqrt(14.0 / 3267.0);
    const auto fs_1_99_462 = std::sqrt(14.0 / 297.0);
    const auto fs_25_1386_42 = std::sqrt(625.0 / 45738.0);
    const auto fs_25_16_42 = std::sqrt(13125.0 / 128.0);
    const auto fs_25_21_3 = std::sqrt(625.0 / 147.0);
    const auto fs_25_33_3 = std::sqrt(625.0 / 363.0);
    const auto fs_25_4_3 = std::sqrt(1875.0 / 16.0);
    const auto fs_25_84_42 = std::sqrt(625.0 / 168.0);
    const auto fs_27_1001_35 = std::sqrt(3645.0 / 143143.0);
    const auto fs_27_14_5 = std::sqrt(3645.0 / 196.0);
    const auto fs_27_28_35 = std::sqrt(3645.0 / 112.0);
    const auto fs_27_77_35 = std::sqrt(3645.0 / 847.0);
    const auto fs_28_429_11 = std::sqrt(784.0 / 16731.0);
    const auto fs_2_33_11 = std::sqrt(4.0 / 99.0);
    const auto fs_2_33_7 = std::sqrt(28.0 / 1089.0);
    const auto fs_2_99_30 = std::sqrt(40.0 / 3267.0);
    const auto fs_2_99_66 = std::sqrt(8.0 / 297.0);
    const auto fs_35_1287_66 = std::sqrt(2450.0 / 50193.0);
    const auto fs_35_429_14 = std::sqrt(17150.0 / 184041.0);
    const auto fs_35_8_6 = std::sqrt(3675.0 / 32.0);
    const auto fs_45_14_15 = std::sqrt(30375.0 / 196.0);
    const auto fs_45_16_6 = std::sqrt(6075.0 / 128.0);
    const auto fs_45_28_21 = std::sqrt(6075.0 / 112.0);
    const auto fs_49_429_10 = std::sqrt(24010.0 / 184041.0);
    const auto fs_49_429_2 = std::sqrt(4802.0 / 184041.0);
    const auto fs_4_33_5 = std::sqrt(80.0 / 1089.0);
    const auto fs_4_99_15 = std::sqrt(80.0 / 3267.0);
    const auto fs_50_693_3 = std::sqrt(2500.0 / 160083.0);
    const auto fs_54_1001_5 = std::sqrt(14580.0 / 1002001.0);
    const auto fs_54_77_5 = std::sqrt(14580.0 / 5929.0);
    const auto fs_5_11_11 = std::sqrt(25.0 / 11.0);
    const auto fs_5_11_7 = std::sqrt(175.0 / 121.0);
    const auto fs_5_14_21 = std::sqrt(75.0 / 28.0);
    const auto fs_5_154_6 = std::sqrt(75.0 / 11858.0);
    const auto fs_5_21_21 = std::sqrt(25.0 / 21.0);
    const auto fs_5_22_11 = std::sqrt(25.0 / 44.0);
    const auto fs_5_22_30 = std::sqrt(375.0 / 242.0);
    const auto fs_5_231_21 = std::sqrt(25.0 / 2541.0);
    const auto fs_5_33_30 = std::sqrt(250.0 / 363.0);
    const auto fs_5_33_66 = std::sqrt(50.0 / 33.0);
    const auto fs_5_42_30 = std::sqrt(125.0 / 294.0);
    const auto fs_5_44_42 = std::sqrt(525.0 / 968.0);
    const auto fs_5_4_21 = std::sqrt(525.0 / 16.0);
    const auto fs_5_66_105 = std::sqrt(875.0 / 1452.0);
    const auto fs_5_66_210 = std::sqrt(875.0 / 726.0);
    const auto fs_5_66_42 = std::sqrt(175.0 / 726.0);
    const auto fs_5_66_462 = std::sqrt(175.0 / 66.0);
    const auto fs_5_693_30 = std::sqrt(250.0 / 160083.0);
    const auto fs_5_6_6 = std::sqrt(25.0 / 6.0);
    const auto fs_5_7_15 = std::sqrt(375.0 / 49.0);
    const auto fs_5_8_30 = std::sqrt(375.0 / 32.0);
    const auto fs_5_99_6 = std::sqrt(50.0 / 3267.0);
    const auto fs_65_132_6 = std::sqrt(4225.0 / 2904.0);
    const auto fs_75_14_3 = std::sqrt(16875.0 / 196.0);
    const auto fs_75_56_42 = std::sqrt(16875.0 / 224.0);
    const auto fs_7_1287_2310 = std::sqrt(3430.0 / 50193.0);
    const auto fs_7_2574_330 = std::sqrt(245.0 / 100386.0);
    const auto fs_7_2574_6006 = std::sqrt(343.0 / 7722.0);
    const auto fs_7_429_14 = std::sqrt(686.0 / 184041.0);
    const auto fs_7_429_286 = std::sqrt(98.0 / 1287.0);
    const auto fs_7_429_55 = std::sqrt(245.0 / 16731.0);
    const auto fs_7_429_715 = std::sqrt(245.0 / 1287.0);
    const auto fs_7_858_10 = std::sqrt(245.0 / 368082.0);
    const auto fs_7_858_1430 = std::sqrt(245.0 / 2574.0);
    const auto fs_7_858_2 = std::sqrt(49.0 / 368082.0);
    const auto fs_7_858_286 = std::sqrt(49.0 / 2574.0);
    const auto fs_9_1001_35 = std::sqrt(405.0 / 143143.0);
    const auto fs_9_11_5 = std::sqrt(405.0 / 121.0);
    const auto fs_9_143_5 = std::sqrt(405.0 / 20449.0);
    const auto fs_9_28_35 = std::sqrt(405.0 / 112.0);
    const auto fs_9_4_5 = std::sqrt(405.0 / 16.0);
    const auto fs_9_77_35 = std::sqrt(405.0 / 847.0);

    // NOTE: the angular components are formed in 24 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_0, ph2_p1, ph4_0, ph4_p1, ph6_0, ph6_p1, ph8_0, ph8_p1, ph8_p7, ph8_p8, ab_2, pc_0, pc_1 : simd::cache_line_size())
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
        const auto h8_p7 = ph8_p7[k];
        const auto h8_p8 = ph8_p8[k];

        pc_0[k] = e_0 * f_105_16 + e_1 * f_35_2 * h2_0 - e_1 * f_35_2 * r_2 + e_2 * f_9_2 * h4_0 - e_2 * f_15_1 * r_2 * h2_0 + e_2 * f_21_2 * r_4 + e_3 * f_10_33 * h6_0 - e_3 * f_18_11 * r_2 * h4_0 + e_3 * f_10_3 * r_4 * h2_0 - e_3 * f_2_1 * r_6 + e_4 * f_7_1287 * h8_0 - e_4 * fs_7_429_715 * h8_p8 - e_4 * f_4_99 * r_2 * h6_0 + e_4 * f_18_143 * r_4 * h4_0 - e_4 * f_20_99 * r_6 * h2_0 + e_4 * f_1_9 * r_8;

        pc_1[k] = - e_1 * fs_35_8_6 * h2_p1 - e_2 * fs_9_4_5 * h4_p1 + e_2 * fs_15_4_6 * r_2 * h2_p1 - e_3 * fs_5_66_42 * h6_p1 + e_3 * fs_9_11_5 * r_2 * h4_p1 - e_3 * fs_5_6_6 * r_4 * h2_p1 - e_4 * fs_7_858_2 * h8_p1 - e_4 * fs_7_858_1430 * h8_p7 + e_4 * fs_1_99_42 * r_2 * h6_p1 - e_4 * fs_9_143_5 * r_4 * h4_p1 + e_4 * fs_5_99_6 * r_6 * h2_p1;
    }

    // NOTE: the angular components are formed in 24 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, ph2_p2, ph4_p2, ph4_p3, ph6_p2, ph6_p3, ph6_p5, ph6_p6, ph8_p2, ph8_p3, ph8_p5, ph8_p6, ab_2, pc_2, pc_3 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

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

        pc_2[k] = e_1 * fs_5_4_21 * h2_p2 + e_2 * fs_27_28_35 * h4_p2 - e_2 * fs_15_14_21 * r_2 * h2_p2 + e_3 * fs_5_33_30 * h6_p2 - e_3 * fs_5_33_66 * h6_p6 - e_3 * fs_27_77_35 * r_2 * h4_p2 + e_3 * fs_5_21_21 * r_4 * h2_p2 + e_4 * fs_7_858_10 * h8_p2 - e_4 * fs_7_2574_6006 * h8_p6 - e_4 * fs_2_99_30 * r_2 * h6_p2 + e_4 * fs_2_99_66 * r_2 * h6_p6 + e_4 * fs_27_1001_35 * r_4 * h4_p2 - e_4 * fs_10_693_21 * r_6 * h2_p2;

        pc_3[k] = - e_2 * fs_9_4_5 * h4_p3 - e_3 * fs_10_33_15 * h6_p3 - e_3 * fs_5_11_11 * h6_p5 + e_3 * fs_9_11_5 * r_2 * h4_p3 - e_4 * fs_7_2574_330 * h8_p3 - e_4 * fs_7_858_286 * h8_p5 + e_4 * fs_4_99_15 * r_2 * h6_p3 + e_4 * fs_2_33_11 * r_2 * h6_p5 - e_4 * fs_9_143_5 * r_4 * h4_p3;
    }

    // NOTE: the angular components are formed in 24 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_2, pe_3, pe_4, ph4_m4, ph4_m3, ph6_m5, ph6_m4, ph6_m3, ph8_m5, ph8_m4, ph8_m3, ab_2, pc_4, pc_5 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

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

        pc_4[k] = e_2 * f_9_2 * h4_m4 + e_3 * fs_10_11_5 * h6_m4 - e_3 * f_18_11 * r_2 * h4_m4 + e_4 * fs_7_429_55 * h8_m4 - e_4 * fs_4_33_5 * r_2 * h6_m4 + e_4 * f_18_143 * r_4 * h4_m4;

        pc_5[k] = - e_2 * fs_9_4_5 * h4_m3 + e_3 * fs_5_11_11 * h6_m5 - e_3 * fs_10_33_15 * h6_m3 + e_3 * fs_9_11_5 * r_2 * h4_m3 + e_4 * fs_7_858_286 * h8_m5 - e_4 * fs_7_2574_330 * h8_m3 - e_4 * fs_2_33_11 * r_2 * h6_m5 + e_4 * fs_4_99_15 * r_2 * h6_m3 - e_4 * fs_9_143_5 * r_4 * h4_m3;
    }

    // NOTE: the angular components are formed in 24 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, ph2_m2, ph2_m1, ph4_m2, ph4_m1, ph6_m6, ph6_m2, ph6_m1, ph8_m8, ph8_m7, ph8_m6, ph8_m2, ph8_m1, ab_2, pc_6, pc_7, pc_8 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

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

        pc_6[k] = e_1 * fs_5_4_21 * h2_m2 + e_2 * fs_27_28_35 * h4_m2 - e_2 * fs_15_14_21 * r_2 * h2_m2 + e_3 * fs_5_33_66 * h6_m6 + e_3 * fs_5_33_30 * h6_m2 - e_3 * fs_27_77_35 * r_2 * h4_m2 + e_3 * fs_5_21_21 * r_4 * h2_m2 + e_4 * fs_7_2574_6006 * h8_m6 + e_4 * fs_7_858_10 * h8_m2 - e_4 * fs_2_99_66 * r_2 * h6_m6 - e_4 * fs_2_99_30 * r_2 * h6_m2 + e_4 * fs_27_1001_35 * r_4 * h4_m2 - e_4 * fs_10_693_21 * r_6 * h2_m2;

        pc_7[k] = - e_1 * fs_35_8_6 * h2_m1 - e_2 * fs_9_4_5 * h4_m1 + e_2 * fs_15_4_6 * r_2 * h2_m1 - e_3 * fs_5_66_42 * h6_m1 + e_3 * fs_9_11_5 * r_2 * h4_m1 - e_3 * fs_5_6_6 * r_4 * h2_m1 + e_4 * fs_7_858_1430 * h8_m7 - e_4 * fs_7_858_2 * h8_m1 + e_4 * fs_1_99_42 * r_2 * h6_m1 - e_4 * fs_9_143_5 * r_4 * h4_m1 + e_4 * fs_5_99_6 * r_6 * h2_m1;

        pc_8[k] = e_4 * fs_7_429_715 * h8_m8;
    }

    // NOTE: the angular components are formed in 24 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_0, ph2_p1, ph4_0, ph4_p1, ph6_0, ph6_p1, ph6_p5, ph6_p6, ph8_0, ph8_p1, ph8_p5, ph8_p6, ab_2, pc_9, pc_10 : simd::cache_line_size())
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
        const auto h6_p5 = ph6_p5[k];
        const auto h6_p6 = ph6_p6[k];
        const auto h8_0 = ph8_0[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h8_p5 = ph8_p5[k];
        const auto h8_p6 = ph8_p6[k];

        pc_9[k] = e_0 * f_105_16 + e_1 * f_35_8 * h2_0 - e_1 * f_35_2 * r_2 - e_2 * f_27_4 * h4_0 - e_2 * f_15_4 * r_2 * h2_0 + e_2 * f_21_2 * r_4 - e_3 * f_85_66 * h6_0 + e_3 * fs_5_66_462 * h6_p6 + e_3 * f_27_11 * r_2 * h4_0 + e_3 * f_5_6 * r_4 * h2_0 - e_3 * f_2_1 * r_6 - e_4 * f_56_1287 * h8_0 - e_4 * fs_14_1287_858 * h8_p6 + e_4 * f_17_99 * r_2 * h6_0 - e_4 * fs_1_99_462 * r_2 * h6_p6 - e_4 * f_27_143 * r_4 * h4_0 - e_4 * f_5_99 * r_6 * h2_0 + e_4 * f_1_9 * r_8;

        pc_10[k] = - e_1 * fs_25_16_42 * h2_p1 + e_2 * fs_9_28_35 * h4_p1 + e_2 * fs_75_56_42 * r_2 * h2_p1 + e_3 * fs_65_132_6 * h6_p1 + e_3 * fs_5_22_11 * h6_p5 - e_3 * fs_9_77_35 * r_2 * h4_p1 - e_3 * fs_25_84_42 * r_4 * h2_p1 + e_4 * fs_7_429_14 * h8_p1 - e_4 * fs_7_429_286 * h8_p5 - e_4 * fs_13_198_6 * r_2 * h6_p1 - e_4 * fs_1_33_11 * r_2 * h6_p5 + e_4 * fs_9_1001_35 * r_4 * h4_p1 + e_4 * fs_25_1386_42 * r_6 * h2_p1;
    }

    // NOTE: the angular components are formed in 24 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, ph2_p2, ph4_m3, ph4_p2, ph4_p4, ph6_m3, ph6_p2, ph6_p4, ph8_m3, ph8_p2, ph8_p4, ab_2, pc_11, pc_12 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

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

        pc_11[k] = e_1 * fs_15_8_21 * h2_p2 + e_2 * fs_9_28_35 * h4_p2 + e_2 * fs_9_4_5 * h4_p4 - e_2 * fs_45_28_21 * r_2 * h2_p2 - e_3 * fs_5_22_30 * h6_p2 - e_3 * f_5_22 * h6_p4 - e_3 * fs_9_77_35 * r_2 * h4_p2 - e_3 * fs_9_11_5 * r_2 * h4_p4 + e_3 * fs_5_14_21 * r_4 * h2_p2 - e_4 * fs_14_429_10 * h8_p2 - e_4 * fs_28_429_11 * h8_p4 + e_4 * fs_1_33_30 * r_2 * h6_p2 + e_4 * f_1_33 * r_2 * h6_p4 + e_4 * fs_9_1001_35 * r_4 * h4_p2 + e_4 * fs_9_143_5 * r_4 * h4_p4 - e_4 * fs_5_231_21 * r_6 * h2_p2;

        pc_12[k] = - e_2 * f_27_4 * h4_m3 + e_3 * fs_25_33_3 * h6_m3 + e_3 * f_27_11 * r_2 * h4_m3 + e_4 * fs_35_1287_66 * h8_m3 - e_4 * fs_10_99_3 * r_2 * h6_m3 - e_4 * f_27_143 * r_4 * h4_m3;
    }

    // NOTE: the angular components are formed in 24 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, ph2_m2, ph2_m1, ph4_m4, ph4_m2, ph4_m1, ph6_m5, ph6_m4, ph6_m2, ph6_m1, ph8_m5, ph8_m4, ph8_m2, ph8_m1, ab_2, pc_13, pc_14 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

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

        pc_13[k] = e_1 * fs_15_8_21 * h2_m2 - e_2 * fs_9_4_5 * h4_m4 + e_2 * fs_9_28_35 * h4_m2 - e_2 * fs_45_28_21 * r_2 * h2_m2 + e_3 * f_5_22 * h6_m4 - e_3 * fs_5_22_30 * h6_m2 + e_3 * fs_9_11_5 * r_2 * h4_m4 - e_3 * fs_9_77_35 * r_2 * h4_m2 + e_3 * fs_5_14_21 * r_4 * h2_m2 + e_4 * fs_28_429_11 * h8_m4 - e_4 * fs_14_429_10 * h8_m2 - e_4 * f_1_33 * r_2 * h6_m4 + e_4 * fs_1_33_30 * r_2 * h6_m2 - e_4 * fs_9_143_5 * r_4 * h4_m4 + e_4 * fs_9_1001_35 * r_4 * h4_m2 - e_4 * fs_5_231_21 * r_6 * h2_m2;

        pc_14[k] = - e_1 * fs_25_16_42 * h2_m1 + e_2 * fs_9_28_35 * h4_m1 + e_2 * fs_75_56_42 * r_2 * h2_m1 - e_3 * fs_5_22_11 * h6_m5 + e_3 * fs_65_132_6 * h6_m1 - e_3 * fs_9_77_35 * r_2 * h4_m1 - e_3 * fs_25_84_42 * r_4 * h2_m1 + e_4 * fs_7_429_286 * h8_m5 + e_4 * fs_7_429_14 * h8_m1 + e_4 * fs_1_33_11 * r_2 * h6_m5 - e_4 * fs_13_198_6 * r_2 * h6_m1 + e_4 * fs_9_1001_35 * r_4 * h4_m1 + e_4 * fs_25_1386_42 * r_6 * h2_m1;
    }

    // NOTE: the angular components are formed in 24 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, ph2_m1, ph4_m1, ph6_m6, ph6_m1, ph8_m7, ph8_m6, ph8_m1, ab_2, pc_15, pc_16 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h2_m1 = ph2_m1[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h6_m6 = ph6_m6[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h8_m7 = ph8_m7[k];
        const auto h8_m6 = ph8_m6[k];
        const auto h8_m1 = ph8_m1[k];

        pc_15[k] = - e_3 * fs_5_66_462 * h6_m6 + e_4 * fs_14_1287_858 * h8_m6 + e_4 * fs_1_99_462 * r_2 * h6_m6;

        pc_16[k] = e_1 * fs_35_8_6 * h2_m1 + e_2 * fs_9_4_5 * h4_m1 - e_2 * fs_15_4_6 * r_2 * h2_m1 + e_3 * fs_5_66_42 * h6_m1 - e_3 * fs_9_11_5 * r_2 * h4_m1 + e_3 * fs_5_6_6 * r_4 * h2_m1 + e_4 * fs_7_858_1430 * h8_m7 + e_4 * fs_7_858_2 * h8_m1 - e_4 * fs_1_99_42 * r_2 * h6_m1 + e_4 * fs_9_143_5 * r_4 * h4_m1 - e_4 * fs_5_99_6 * r_6 * h2_m1;
    }

    // NOTE: the angular components are formed in 24 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_0, ph4_0, ph4_p4, ph6_0, ph6_p4, ph8_0, ph8_p4, ab_2, pc_17 : simd::cache_line_size())
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
        const auto h4_p4 = ph4_p4[k];
        const auto h6_0 = ph6_0[k];
        const auto h6_p4 = ph6_p4[k];
        const auto h8_0 = ph8_0[k];
        const auto h8_p4 = ph8_p4[k];

        pc_17[k] = e_0 * f_105_16 - e_1 * f_5_1 * h2_0 - e_1 * f_35_2 * r_2 - e_2 * f_99_28 * h4_0 - e_2 * fs_27_28_35 * h4_p4 + e_2 * f_30_7 * r_2 * h2_0 + e_2 * f_21_2 * r_4 + e_3 * f_5_3 * h6_0 + e_3 * fs_5_11_7 * h6_p4 + e_3 * f_9_7 * r_2 * h4_0 + e_3 * fs_27_77_35 * r_2 * h4_p4 - e_3 * f_20_21 * r_4 * h2_0 - e_3 * f_2_1 * r_6 + e_4 * f_196_1287 * h8_0 - e_4 * fs_14_429_77 * h8_p4 - e_4 * f_2_9 * r_2 * h6_0 - e_4 * fs_2_33_7 * r_2 * h6_p4 - e_4 * f_9_91 * r_4 * h4_0 - e_4 * fs_27_1001_35 * r_4 * h4_p4 + e_4 * f_40_693 * r_6 * h2_0 + e_4 * f_1_9 * r_8;
    }

    // NOTE: the angular components are formed in 24 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, ph2_m2, ph2_p1, ph4_m2, ph4_p1, ph4_p3, ph6_p1, ph6_p3, ph8_m2, ph8_p1, ph8_p3, ab_2, pc_18, pc_19 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

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

        pc_18[k] = - e_1 * fs_45_16_6 * h2_p1 + e_2 * fs_27_14_5 * h4_p1 - e_2 * fs_9_28_35 * h4_p3 + e_2 * fs_135_56_6 * r_2 * h2_p1 - e_3 * fs_5_44_42 * h6_p1 + e_3 * fs_5_66_105 * h6_p3 - e_3 * fs_54_77_5 * r_2 * h4_p1 + e_3 * fs_9_77_35 * r_2 * h4_p3 - e_3 * fs_15_28_6 * r_4 * h2_p1 - e_4 * fs_49_429_2 * h8_p1 - e_4 * fs_7_1287_2310 * h8_p3 + e_4 * fs_1_66_42 * r_2 * h6_p1 - e_4 * fs_1_99_105 * r_2 * h6_p3 + e_4 * fs_54_1001_5 * r_4 * h4_p1 - e_4 * fs_9_1001_35 * r_4 * h4_p3 + e_4 * fs_5_154_6 * r_6 * h2_p1;

        pc_19[k] = e_1 * fs_15_4_15 * h2_m2 - e_2 * f_99_28 * h4_m2 - e_2 * fs_45_14_15 * r_2 * h2_m2 + e_3 * f_9_7 * r_2 * h4_m2 + e_3 * fs_5_7_15 * r_4 * h2_m2 + e_4 * fs_35_429_14 * h8_m2 - e_4 * f_9_91 * r_4 * h4_m2 - e_4 * fs_10_231_15 * r_6 * h2_m2;
    }

    // NOTE: the angular components are formed in 24 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, ph2_m1, ph4_m4, ph4_m3, ph4_m1, ph6_m5, ph6_m4, ph6_m3, ph6_m1, ph8_m5, ph8_m4, ph8_m3, ph8_m1, ab_2, pc_20, pc_21, pc_22 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

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

        pc_20[k] = - e_1 * fs_45_16_6 * h2_m1 + e_2 * fs_9_28_35 * h4_m3 + e_2 * fs_27_14_5 * h4_m1 + e_2 * fs_135_56_6 * r_2 * h2_m1 - e_3 * fs_5_66_105 * h6_m3 - e_3 * fs_5_44_42 * h6_m1 - e_3 * fs_9_77_35 * r_2 * h4_m3 - e_3 * fs_54_77_5 * r_2 * h4_m1 - e_3 * fs_15_28_6 * r_4 * h2_m1 + e_4 * fs_7_1287_2310 * h8_m3 - e_4 * fs_49_429_2 * h8_m1 + e_4 * fs_1_99_105 * r_2 * h6_m3 + e_4 * fs_1_66_42 * r_2 * h6_m1 + e_4 * fs_9_1001_35 * r_4 * h4_m3 + e_4 * fs_54_1001_5 * r_4 * h4_m1 + e_4 * fs_5_154_6 * r_6 * h2_m1;

        pc_21[k] = e_2 * fs_27_28_35 * h4_m4 - e_3 * fs_5_11_7 * h6_m4 - e_3 * fs_27_77_35 * r_2 * h4_m4 + e_4 * fs_14_429_77 * h8_m4 + e_4 * fs_2_33_7 * r_2 * h6_m4 + e_4 * fs_27_1001_35 * r_4 * h4_m4;

        pc_22[k] = e_1 * fs_25_16_42 * h2_m1 - e_2 * fs_9_28_35 * h4_m1 - e_2 * fs_75_56_42 * r_2 * h2_m1 - e_3 * fs_5_22_11 * h6_m5 - e_3 * fs_65_132_6 * h6_m1 + e_3 * fs_9_77_35 * r_2 * h4_m1 + e_3 * fs_25_84_42 * r_4 * h2_m1 + e_4 * fs_7_429_286 * h8_m5 - e_4 * fs_7_429_14 * h8_m1 + e_4 * fs_1_33_11 * r_2 * h6_m5 + e_4 * fs_13_198_6 * r_2 * h6_m1 - e_4 * fs_9_1001_35 * r_4 * h4_m1 - e_4 * fs_25_1386_42 * r_6 * h2_m1;
    }

    // NOTE: the angular components are formed in 24 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, ph2_m2, ph4_m2, ph6_m6, ph6_m2, ph8_m6, ph8_m2, ab_2, pc_23 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h2_m2 = ph2_m2[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h6_m6 = ph6_m6[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h8_m6 = ph8_m6[k];
        const auto h8_m2 = ph8_m2[k];

        pc_23[k] = - e_1 * fs_5_4_21 * h2_m2 - e_2 * fs_27_28_35 * h4_m2 + e_2 * fs_15_14_21 * r_2 * h2_m2 + e_3 * fs_5_33_66 * h6_m6 - e_3 * fs_5_33_30 * h6_m2 + e_3 * fs_27_77_35 * r_2 * h4_m2 - e_3 * fs_5_21_21 * r_4 * h2_m2 + e_4 * fs_7_2574_6006 * h8_m6 - e_4 * fs_7_858_10 * h8_m2 - e_4 * fs_2_99_66 * r_2 * h6_m6 + e_4 * fs_2_99_30 * r_2 * h6_m2 - e_4 * fs_27_1001_35 * r_4 * h4_m2 + e_4 * fs_10_693_21 * r_6 * h2_m2;
    }

    // NOTE: the angular components are formed in 24 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_m1, ph2_0, ph2_p2, ph4_m1, ph4_0, ph4_p2, ph6_m1, ph6_0, ph6_p2, ph8_m1, ph8_0, ph8_p2, ab_2, pc_24, pc_25 : simd::cache_line_size())
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

        pc_24[k] = e_0 * f_105_16 - e_1 * f_85_8 * h2_0 + e_1 * fs_25_4_3 * h2_p2 - e_1 * f_35_2 * r_2 + e_2 * f_81_28 * h4_0 - e_2 * fs_27_14_5 * h4_p2 + e_2 * f_255_28 * r_2 * h2_0 - e_2 * fs_75_14_3 * r_2 * h2_p2 + e_2 * f_21_2 * r_4 + e_3 * f_5_66 * h6_0 + e_3 * fs_5_66_210 * h6_p2 - e_3 * f_81_77 * r_2 * h4_0 + e_3 * fs_54_77_5 * r_2 * h4_p2 - e_3 * f_85_42 * r_4 * h2_0 + e_3 * fs_25_21_3 * r_4 * h2_p2 - e_3 * f_2_1 * r_6 - e_4 * f_392_1287 * h8_0 - e_4 * fs_14_429_70 * h8_p2 - e_4 * f_1_99 * r_2 * h6_0 - e_4 * fs_1_99_210 * r_2 * h6_p2 + e_4 * f_81_1001 * r_4 * h4_0 - e_4 * fs_54_1001_5 * r_4 * h4_p2 + e_4 * f_85_693 * r_6 * h2_0 - e_4 * fs_50_693_3 * r_6 * h2_p2 + e_4 * f_1_9 * r_8;

        pc_25[k] = - e_1 * fs_5_8_30 * h2_m1 + e_2 * f_81_28 * h4_m1 + e_2 * fs_15_28_30 * r_2 * h2_m1 - e_3 * fs_5_66_210 * h6_m1 - e_3 * f_81_77 * r_2 * h4_m1 - e_3 * fs_5_42_30 * r_4 * h2_m1 + e_4 * fs_49_429_10 * h8_m1 + e_4 * fs_1_99_210 * r_2 * h6_m1 + e_4 * f_81_1001 * r_4 * h4_m1 + e_4 * fs_5_693_30 * r_6 * h2_m1;
    }

    // NOTE: the angular components are formed in 24 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, ph2_m2, ph2_m1, ph4_m3, ph4_m2, ph4_m1, ph6_m3, ph6_m2, ph6_m1, ph8_m3, ph8_m2, ph8_m1, ab_2, pc_26, pc_27 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

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

        pc_26[k] = - e_1 * fs_25_4_3 * h2_m2 + e_2 * fs_27_14_5 * h4_m2 + e_2 * fs_75_14_3 * r_2 * h2_m2 - e_3 * fs_5_66_210 * h6_m2 - e_3 * fs_54_77_5 * r_2 * h4_m2 - e_3 * fs_25_21_3 * r_4 * h2_m2 + e_4 * fs_14_429_70 * h8_m2 + e_4 * fs_1_99_210 * r_2 * h6_m2 + e_4 * fs_54_1001_5 * r_4 * h4_m2 + e_4 * fs_50_693_3 * r_6 * h2_m2;

        pc_27[k] = e_1 * fs_45_16_6 * h2_m1 + e_2 * fs_9_28_35 * h4_m3 - e_2 * fs_27_14_5 * h4_m1 - e_2 * fs_135_56_6 * r_2 * h2_m1 - e_3 * fs_5_66_105 * h6_m3 + e_3 * fs_5_44_42 * h6_m1 - e_3 * fs_9_77_35 * r_2 * h4_m3 + e_3 * fs_54_77_5 * r_2 * h4_m1 + e_3 * fs_15_28_6 * r_4 * h2_m1 + e_4 * fs_7_1287_2310 * h8_m3 + e_4 * fs_49_429_2 * h8_m1 + e_4 * fs_1_99_105 * r_2 * h6_m3 - e_4 * fs_1_66_42 * r_2 * h6_m1 + e_4 * fs_9_1001_35 * r_4 * h4_m3 - e_4 * fs_54_1001_5 * r_4 * h4_m1 - e_4 * fs_5_154_6 * r_6 * h2_m1;
    }

    // NOTE: the angular components are formed in 24 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, ph2_m2, ph4_m4, ph4_m3, ph4_m2, ph6_m5, ph6_m4, ph6_m3, ph6_m2, ph8_m5, ph8_m4, ph8_m3, ph8_m2, ab_2, pc_28, pc_29 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

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

        pc_28[k] = - e_1 * fs_15_8_21 * h2_m2 - e_2 * fs_9_4_5 * h4_m4 - e_2 * fs_9_28_35 * h4_m2 + e_2 * fs_45_28_21 * r_2 * h2_m2 + e_3 * f_5_22 * h6_m4 + e_3 * fs_5_22_30 * h6_m2 + e_3 * fs_9_11_5 * r_2 * h4_m4 + e_3 * fs_9_77_35 * r_2 * h4_m2 - e_3 * fs_5_14_21 * r_4 * h2_m2 + e_4 * fs_28_429_11 * h8_m4 + e_4 * fs_14_429_10 * h8_m2 - e_4 * f_1_33 * r_2 * h6_m4 - e_4 * fs_1_33_30 * r_2 * h6_m2 - e_4 * fs_9_143_5 * r_4 * h4_m4 - e_4 * fs_9_1001_35 * r_4 * h4_m2 + e_4 * fs_5_231_21 * r_6 * h2_m2;

        pc_29[k] = e_2 * fs_9_4_5 * h4_m3 + e_3 * fs_5_11_11 * h6_m5 + e_3 * fs_10_33_15 * h6_m3 - e_3 * fs_9_11_5 * r_2 * h4_m3 + e_4 * fs_7_858_286 * h8_m5 + e_4 * fs_7_2574_330 * h8_m3 - e_4 * fs_2_33_11 * r_2 * h6_m5 - e_4 * fs_4_99_15 * r_2 * h6_m3 + e_4 * fs_9_143_5 * r_4 * h4_m3;
    }

    // NOTE: the angular components are formed in 24 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_0, ph2_p1, ph2_p2, ph4_0, ph4_p1, ph4_p2, ph6_0, ph6_p1, ph8_0, ph8_p1, ph8_p2, ab_2, pc_30, pc_31, pc_32 : simd::cache_line_size())
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
        const auto h2_p2 = ph2_p2[k];
        const auto h4_0 = ph4_0[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h6_0 = ph6_0[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h8_0 = ph8_0[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h8_p2 = ph8_p2[k];

        pc_30[k] = e_0 * f_105_16 - e_1 * f_25_2 * h2_0 - e_1 * f_35_2 * r_2 + e_2 * f_81_14 * h4_0 + e_2 * f_75_7 * r_2 * h2_0 + e_2 * f_21_2 * r_4 - e_3 * f_50_33 * h6_0 - e_3 * f_162_77 * r_2 * h4_0 - e_3 * f_50_21 * r_4 * h2_0 - e_3 * f_2_1 * r_6 + e_4 * f_490_1287 * h8_0 + e_4 * f_20_99 * r_2 * h6_0 + e_4 * f_162_1001 * r_4 * h4_0 + e_4 * f_100_693 * r_6 * h2_0 + e_4 * f_1_9 * r_8;

        pc_31[k] = - e_1 * fs_5_8_30 * h2_p1 + e_2 * f_81_28 * h4_p1 + e_2 * fs_15_28_30 * r_2 * h2_p1 - e_3 * fs_5_66_210 * h6_p1 - e_3 * f_81_77 * r_2 * h4_p1 - e_3 * fs_5_42_30 * r_4 * h2_p1 + e_4 * fs_49_429_10 * h8_p1 + e_4 * fs_1_99_210 * r_2 * h6_p1 + e_4 * f_81_1001 * r_4 * h4_p1 + e_4 * fs_5_693_30 * r_6 * h2_p1;

        pc_32[k] = e_1 * fs_15_4_15 * h2_p2 - e_2 * f_99_28 * h4_p2 - e_2 * fs_45_14_15 * r_2 * h2_p2 + e_3 * f_9_7 * r_2 * h4_p2 + e_3 * fs_5_7_15 * r_4 * h2_p2 + e_4 * fs_35_429_14 * h8_p2 - e_4 * f_9_91 * r_4 * h4_p2 - e_4 * fs_10_231_15 * r_6 * h2_p2;
    }

    // NOTE: the angular components are formed in 24 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_2, pe_3, pe_4, ph4_p3, ph4_p4, ph6_p3, ph6_p4, ph8_p3, ph8_p4, ab_2, pc_33, pc_34 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

        const auto h4_p3 = ph4_p3[k];
        const auto h4_p4 = ph4_p4[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h6_p4 = ph6_p4[k];
        const auto h8_p3 = ph8_p3[k];
        const auto h8_p4 = ph8_p4[k];

        pc_33[k] = - e_2 * f_27_4 * h4_p3 + e_3 * fs_25_33_3 * h6_p3 + e_3 * f_27_11 * r_2 * h4_p3 + e_4 * fs_35_1287_66 * h8_p3 - e_4 * fs_10_99_3 * r_2 * h6_p3 - e_4 * f_27_143 * r_4 * h4_p3;

        pc_34[k] = e_2 * f_9_2 * h4_p4 + e_3 * fs_10_11_5 * h6_p4 - e_3 * f_18_11 * r_2 * h4_p4 + e_4 * fs_7_429_55 * h8_p4 - e_4 * fs_4_33_5 * r_2 * h6_p4 + e_4 * f_18_143 * r_4 * h4_p4;
    }

    // NOTE: the angular components are formed in 24 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_0, ph2_p2, ph4_0, ph4_p2, ph6_0, ph6_p2, ph8_0, ph8_p2, ab_2, pc_35 : simd::cache_line_size())
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

        pc_35[k] = e_0 * f_105_16 - e_1 * f_85_8 * h2_0 - e_1 * fs_25_4_3 * h2_p2 - e_1 * f_35_2 * r_2 + e_2 * f_81_28 * h4_0 + e_2 * fs_27_14_5 * h4_p2 + e_2 * f_255_28 * r_2 * h2_0 + e_2 * fs_75_14_3 * r_2 * h2_p2 + e_2 * f_21_2 * r_4 + e_3 * f_5_66 * h6_0 - e_3 * fs_5_66_210 * h6_p2 - e_3 * f_81_77 * r_2 * h4_0 - e_3 * fs_54_77_5 * r_2 * h4_p2 - e_3 * f_85_42 * r_4 * h2_0 - e_3 * fs_25_21_3 * r_4 * h2_p2 - e_3 * f_2_1 * r_6 - e_4 * f_392_1287 * h8_0 + e_4 * fs_14_429_70 * h8_p2 - e_4 * f_1_99 * r_2 * h6_0 + e_4 * fs_1_99_210 * r_2 * h6_p2 + e_4 * f_81_1001 * r_4 * h4_0 + e_4 * fs_54_1001_5 * r_4 * h4_p2 + e_4 * f_85_693 * r_6 * h2_0 + e_4 * fs_50_693_3 * r_6 * h2_p2 + e_4 * f_1_9 * r_8;
    }

    // NOTE: the angular components are formed in 24 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, ph2_p1, ph4_p1, ph4_p3, ph6_p1, ph6_p3, ph8_p1, ph8_p3, ab_2, pc_36 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

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

        pc_36[k] = - e_1 * fs_45_16_6 * h2_p1 + e_2 * fs_27_14_5 * h4_p1 + e_2 * fs_9_28_35 * h4_p3 + e_2 * fs_135_56_6 * r_2 * h2_p1 - e_3 * fs_5_44_42 * h6_p1 - e_3 * fs_5_66_105 * h6_p3 - e_3 * fs_54_77_5 * r_2 * h4_p1 - e_3 * fs_9_77_35 * r_2 * h4_p3 - e_3 * fs_15_28_6 * r_4 * h2_p1 - e_4 * fs_49_429_2 * h8_p1 + e_4 * fs_7_1287_2310 * h8_p3 + e_4 * fs_1_66_42 * r_2 * h6_p1 + e_4 * fs_1_99_105 * r_2 * h6_p3 + e_4 * fs_54_1001_5 * r_4 * h4_p1 + e_4 * fs_9_1001_35 * r_4 * h4_p3 + e_4 * fs_5_154_6 * r_6 * h2_p1;
    }

    // NOTE: the angular components are formed in 24 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, ph2_p2, ph4_p2, ph4_p3, ph4_p4, ph6_p2, ph6_p3, ph6_p4, ph6_p5, ph8_p2, ph8_p3, ph8_p4, ph8_p5, ab_2, pc_37, pc_38 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

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
        const auto h6_p5 = ph6_p5[k];
        const auto h8_p2 = ph8_p2[k];
        const auto h8_p3 = ph8_p3[k];
        const auto h8_p4 = ph8_p4[k];
        const auto h8_p5 = ph8_p5[k];

        pc_37[k] = e_1 * fs_15_8_21 * h2_p2 + e_2 * fs_9_28_35 * h4_p2 - e_2 * fs_9_4_5 * h4_p4 - e_2 * fs_45_28_21 * r_2 * h2_p2 - e_3 * fs_5_22_30 * h6_p2 + e_3 * f_5_22 * h6_p4 - e_3 * fs_9_77_35 * r_2 * h4_p2 + e_3 * fs_9_11_5 * r_2 * h4_p4 + e_3 * fs_5_14_21 * r_4 * h2_p2 - e_4 * fs_14_429_10 * h8_p2 + e_4 * fs_28_429_11 * h8_p4 + e_4 * fs_1_33_30 * r_2 * h6_p2 - e_4 * f_1_33 * r_2 * h6_p4 + e_4 * fs_9_1001_35 * r_4 * h4_p2 - e_4 * fs_9_143_5 * r_4 * h4_p4 - e_4 * fs_5_231_21 * r_6 * h2_p2;

        pc_38[k] = - e_2 * fs_9_4_5 * h4_p3 - e_3 * fs_10_33_15 * h6_p3 + e_3 * fs_5_11_11 * h6_p5 + e_3 * fs_9_11_5 * r_2 * h4_p3 - e_4 * fs_7_2574_330 * h8_p3 + e_4 * fs_7_858_286 * h8_p5 + e_4 * fs_4_99_15 * r_2 * h6_p3 - e_4 * fs_2_33_11 * r_2 * h6_p5 - e_4 * fs_9_143_5 * r_4 * h4_p3;
    }

    // NOTE: the angular components are formed in 24 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_0, ph4_0, ph4_p4, ph6_0, ph6_p4, ph8_0, ph8_p4, ab_2, pc_39 : simd::cache_line_size())
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
        const auto h4_p4 = ph4_p4[k];
        const auto h6_0 = ph6_0[k];
        const auto h6_p4 = ph6_p4[k];
        const auto h8_0 = ph8_0[k];
        const auto h8_p4 = ph8_p4[k];

        pc_39[k] = e_0 * f_105_16 - e_1 * f_5_1 * h2_0 - e_1 * f_35_2 * r_2 - e_2 * f_99_28 * h4_0 + e_2 * fs_27_28_35 * h4_p4 + e_2 * f_30_7 * r_2 * h2_0 + e_2 * f_21_2 * r_4 + e_3 * f_5_3 * h6_0 - e_3 * fs_5_11_7 * h6_p4 + e_3 * f_9_7 * r_2 * h4_0 - e_3 * fs_27_77_35 * r_2 * h4_p4 - e_3 * f_20_21 * r_4 * h2_0 - e_3 * f_2_1 * r_6 + e_4 * f_196_1287 * h8_0 + e_4 * fs_14_429_77 * h8_p4 - e_4 * f_2_9 * r_2 * h6_0 + e_4 * fs_2_33_7 * r_2 * h6_p4 - e_4 * f_9_91 * r_4 * h4_0 + e_4 * fs_27_1001_35 * r_4 * h4_p4 + e_4 * f_40_693 * r_6 * h2_0 + e_4 * f_1_9 * r_8;
    }

    // NOTE: the angular components are formed in 24 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, ph2_p1, ph2_p2, ph4_p1, ph4_p2, ph6_p1, ph6_p2, ph6_p5, ph6_p6, ph8_p1, ph8_p2, ph8_p5, ph8_p6, ab_2, pc_40, pc_41 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

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

        pc_40[k] = - e_1 * fs_25_16_42 * h2_p1 + e_2 * fs_9_28_35 * h4_p1 + e_2 * fs_75_56_42 * r_2 * h2_p1 + e_3 * fs_65_132_6 * h6_p1 - e_3 * fs_5_22_11 * h6_p5 - e_3 * fs_9_77_35 * r_2 * h4_p1 - e_3 * fs_25_84_42 * r_4 * h2_p1 + e_4 * fs_7_429_14 * h8_p1 + e_4 * fs_7_429_286 * h8_p5 - e_4 * fs_13_198_6 * r_2 * h6_p1 + e_4 * fs_1_33_11 * r_2 * h6_p5 + e_4 * fs_9_1001_35 * r_4 * h4_p1 + e_4 * fs_25_1386_42 * r_6 * h2_p1;

        pc_41[k] = e_1 * fs_5_4_21 * h2_p2 + e_2 * fs_27_28_35 * h4_p2 - e_2 * fs_15_14_21 * r_2 * h2_p2 + e_3 * fs_5_33_30 * h6_p2 + e_3 * fs_5_33_66 * h6_p6 - e_3 * fs_27_77_35 * r_2 * h4_p2 + e_3 * fs_5_21_21 * r_4 * h2_p2 + e_4 * fs_7_858_10 * h8_p2 + e_4 * fs_7_2574_6006 * h8_p6 - e_4 * fs_2_99_30 * r_2 * h6_p2 - e_4 * fs_2_99_66 * r_2 * h6_p6 + e_4 * fs_27_1001_35 * r_4 * h4_p2 - e_4 * fs_10_693_21 * r_6 * h2_p2;
    }

    // NOTE: the angular components are formed in 24 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_0, ph2_p1, ph4_0, ph4_p1, ph6_0, ph6_p1, ph6_p6, ph8_0, ph8_p1, ph8_p6, ph8_p7, ab_2, pc_42, pc_43 : simd::cache_line_size())
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
        const auto h6_p6 = ph6_p6[k];
        const auto h8_0 = ph8_0[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h8_p6 = ph8_p6[k];
        const auto h8_p7 = ph8_p7[k];

        pc_42[k] = e_0 * f_105_16 + e_1 * f_35_8 * h2_0 - e_1 * f_35_2 * r_2 - e_2 * f_27_4 * h4_0 - e_2 * f_15_4 * r_2 * h2_0 + e_2 * f_21_2 * r_4 - e_3 * f_85_66 * h6_0 - e_3 * fs_5_66_462 * h6_p6 + e_3 * f_27_11 * r_2 * h4_0 + e_3 * f_5_6 * r_4 * h2_0 - e_3 * f_2_1 * r_6 - e_4 * f_56_1287 * h8_0 + e_4 * fs_14_1287_858 * h8_p6 + e_4 * f_17_99 * r_2 * h6_0 + e_4 * fs_1_99_462 * r_2 * h6_p6 - e_4 * f_27_143 * r_4 * h4_0 - e_4 * f_5_99 * r_6 * h2_0 + e_4 * f_1_9 * r_8;

        pc_43[k] = - e_1 * fs_35_8_6 * h2_p1 - e_2 * fs_9_4_5 * h4_p1 + e_2 * fs_15_4_6 * r_2 * h2_p1 - e_3 * fs_5_66_42 * h6_p1 + e_3 * fs_9_11_5 * r_2 * h4_p1 - e_3 * fs_5_6_6 * r_4 * h2_p1 - e_4 * fs_7_858_2 * h8_p1 + e_4 * fs_7_858_1430 * h8_p7 + e_4 * fs_1_99_42 * r_2 * h6_p1 - e_4 * fs_9_143_5 * r_4 * h4_p1 + e_4 * fs_5_99_6 * r_6 * h2_p1;
    }

    // NOTE: the angular components are formed in 24 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, ph2_0, ph4_0, ph6_0, ph8_0, ph8_p8, ab_2, pc_44 : simd::cache_line_size())
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

        pc_44[k] = e_0 * f_105_16 + e_1 * f_35_2 * h2_0 - e_1 * f_35_2 * r_2 + e_2 * f_9_2 * h4_0 - e_2 * f_15_1 * r_2 * h2_0 + e_2 * f_21_2 * r_4 + e_3 * f_10_33 * h6_0 - e_3 * f_18_11 * r_2 * h4_0 + e_3 * f_10_3 * r_4 * h2_0 - e_3 * f_2_1 * r_6 + e_4 * f_7_1287 * h8_0 + e_4 * fs_7_429_715 * h8_p8 - e_4 * f_4_99 * r_2 * h6_0 + e_4 * f_18_143 * r_4 * h4_0 - e_4 * f_20_99 * r_6 * h2_0 + e_4 * f_1_9 * r_8;
    }

    // NOTE: the values of a combination of angular components are stored as one
    // row of nvalues columns, with the component on bra side running slowest. The
    // rows which the symmetry relates to an already formed one are copied from it.

    const size_t sources[81] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 1, 10, 11, 12, 13, 14, 15, 16, 17, 2, 11, 20, 21, 22, 23, 24, 25, 26, 3, 12, 21, 30, 31, 32, 33, 34, 35, 4, 13, 22, 31, 40, 41, 42, 43, 44, 5, 14, 23, 32, 41, 50, 51, 52, 53, 6, 15, 24, 33, 42, 51, 60, 61, 62, 7, 16, 25, 34, 43, 52, 61, 70, 71, 8, 17, 26, 35, 44, 53, 62, 71, 80};

    for (size_t n = 0; n < 81; n++)
    {
        if (sources[n] != n) std::copy(values + sources[n] * nvalues, values + (sources[n] + 1) * nvalues, values + n * nvalues);
    }
}

}  // namespace simdt2ceri
