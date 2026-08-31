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



#include "SimdTwoCenterElectronRepulsionRecHH.hpp"

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
compute_hh_electron_repulsion(double                         *values,
                              const size_t                    nvalues,
                              const CBasisFunction           &bra,
                              const CBasisFunction           &ket,
                              const std::vector<CSimdMatrix> &harmonics,
                              const CSimdMatrix              &coordinates) -> void
{
    if ((bra.get_angular_momentum() != 5) || (ket.get_angular_momentum() != 5))
    {
        errors::assertMsgCritical(
            false, std::string("SimdTwoCenterElectronRepulsionFunc.compute_hh_electron_repulsion: Basis functions must be of angular momenta five and five"));
    }

    if (harmonics.size() < 10)
    {
        errors::assertMsgCritical(
            false, std::string("SimdTwoCenterElectronRepulsionFunc.compute_hh_electron_repulsion: Harmonics must reach angular momentum 10"));
    }

    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("SimdTwoCenterElectronRepulsionFunc.compute_hh_electron_repulsion: Number of values exceeds number of atom pairs"));
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
    // orders 5 to 10 alone, and the orders below them are formed on the
    // way to them by the recursion the Boys function is evaluated with.

    simdfunc::compute_boys_function(boys);

    // NOTE: the buffer holds the contracted prefactor of each exponent factor of
    // the terms, as the integrals of the angular components are formed straight
    // into the values and are not written a second time. Every exponent factor is
    // used with one order of the Boys function alone, so the order follows from
    // the factor and one accumulator per factor suffices.

    auto buffer = CSimdMatrix(6, nvalues);

    auto *pe_0 = buffer.data(0);
    auto *pe_1 = buffer.data(1);
    auto *pe_2 = buffer.data(2);
    auto *pe_3 = buffer.data(3);
    auto *pe_4 = buffer.data(4);
    auto *pe_5 = buffer.data(5);

    std::fill(pe_0, pe_0 + nvalues, 0.0);
    std::fill(pe_1, pe_1 + nvalues, 0.0);
    std::fill(pe_2, pe_2 + nvalues, 0.0);
    std::fill(pe_3, pe_3 + nvalues, 0.0);
    std::fill(pe_4, pe_4 + nvalues, 0.0);
    std::fill(pe_5, pe_5 + nvalues, 0.0);

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

            const auto ff_0 = fbase / (fexp * fexp * fexp * fexp * fexp);

            const auto ff_1 = fbase * fmu / (fexp * fexp * fexp * fexp * fexp);

            const auto ff_2 = fbase * fmu * fmu / (fexp * fexp * fexp * fexp * fexp);

            const auto ff_3 = fbase * fmu * fmu * fmu / (fexp * fexp * fexp * fexp * fexp);

            const auto ff_4 = fbase * fmu * fmu * fmu * fmu / (fexp * fexp * fexp * fexp * fexp);

            const auto ff_5 = fbase * fmu * fmu * fmu * fmu * fmu / (fexp * fexp * fexp * fexp * fexp);

            const auto *bv_0 = boys.data(6, i * nprim_b + j);

            const auto *bv_1 = boys.data(7, i * nprim_b + j);

            const auto *bv_2 = boys.data(8, i * nprim_b + j);

            const auto *bv_3 = boys.data(9, i * nprim_b + j);

            const auto *bv_4 = boys.data(10, i * nprim_b + j);

            const auto *bv_5 = boys.data(11, i * nprim_b + j);

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, bv_0, bv_1, bv_2, bv_3, bv_4, bv_5 : simd::cache_line_size())
            for (size_t k = 0; k < nvalues; k++)
            {
                pe_0[k] += ff_0 * bv_0[k];
                pe_1[k] += ff_1 * bv_1[k];
                pe_2[k] += ff_2 * bv_2[k];
                pe_3[k] += ff_3 * bv_3[k];
                pe_4[k] += ff_4 * bv_4[k];
                pe_5[k] += ff_5 * bv_5[k];
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
    auto *pc_11 = values + 12 * nvalues;
    auto *pc_12 = values + 13 * nvalues;
    auto *pc_13 = values + 14 * nvalues;
    auto *pc_14 = values + 15 * nvalues;
    auto *pc_15 = values + 16 * nvalues;
    auto *pc_16 = values + 17 * nvalues;
    auto *pc_17 = values + 18 * nvalues;
    auto *pc_18 = values + 19 * nvalues;
    auto *pc_19 = values + 20 * nvalues;
    auto *pc_20 = values + 21 * nvalues;
    auto *pc_21 = values + 24 * nvalues;
    auto *pc_22 = values + 25 * nvalues;
    auto *pc_23 = values + 26 * nvalues;
    auto *pc_24 = values + 27 * nvalues;
    auto *pc_25 = values + 28 * nvalues;
    auto *pc_26 = values + 29 * nvalues;
    auto *pc_27 = values + 30 * nvalues;
    auto *pc_28 = values + 31 * nvalues;
    auto *pc_29 = values + 32 * nvalues;
    auto *pc_30 = values + 36 * nvalues;
    auto *pc_31 = values + 37 * nvalues;
    auto *pc_32 = values + 38 * nvalues;
    auto *pc_33 = values + 39 * nvalues;
    auto *pc_34 = values + 40 * nvalues;
    auto *pc_35 = values + 41 * nvalues;
    auto *pc_36 = values + 42 * nvalues;
    auto *pc_37 = values + 43 * nvalues;
    auto *pc_38 = values + 48 * nvalues;
    auto *pc_39 = values + 49 * nvalues;
    auto *pc_40 = values + 50 * nvalues;
    auto *pc_41 = values + 51 * nvalues;
    auto *pc_42 = values + 52 * nvalues;
    auto *pc_43 = values + 53 * nvalues;
    auto *pc_44 = values + 54 * nvalues;
    auto *pc_45 = values + 60 * nvalues;
    auto *pc_46 = values + 61 * nvalues;
    auto *pc_47 = values + 62 * nvalues;
    auto *pc_48 = values + 63 * nvalues;
    auto *pc_49 = values + 64 * nvalues;
    auto *pc_50 = values + 65 * nvalues;
    auto *pc_51 = values + 72 * nvalues;
    auto *pc_52 = values + 73 * nvalues;
    auto *pc_53 = values + 74 * nvalues;
    auto *pc_54 = values + 75 * nvalues;
    auto *pc_55 = values + 76 * nvalues;
    auto *pc_56 = values + 84 * nvalues;
    auto *pc_57 = values + 85 * nvalues;
    auto *pc_58 = values + 86 * nvalues;
    auto *pc_59 = values + 87 * nvalues;
    auto *pc_60 = values + 96 * nvalues;
    auto *pc_61 = values + 97 * nvalues;
    auto *pc_62 = values + 98 * nvalues;
    auto *pc_63 = values + 108 * nvalues;
    auto *pc_64 = values + 109 * nvalues;
    auto *pc_65 = values + 120 * nvalues;

    // NOTE: the factors of the terms depend on the angular momenta alone, so they
    // are formed once for the whole matrix instead of once for every atom pair.

    const auto f_100_11 = 100.0 / 11.0;
    const auto f_100_33 = 100.0 / 33.0;
    const auto f_105_16 = 105.0 / 16.0;
    const auto f_10_11 = 10.0 / 11.0;
    const auto f_10_143 = 10.0 / 143.0;
    const auto f_10_187 = 10.0 / 187.0;
    const auto f_10_33 = 10.0 / 33.0;
    const auto f_120_11 = 120.0 / 11.0;
    const auto f_12_143 = 12.0 / 143.0;
    const auto f_13230_46189 = 13230.0 / 46189.0;
    const auto f_135_11 = 135.0 / 11.0;
    const auto f_135_2 = 135.0 / 2.0;
    const auto f_135_286 = 135.0 / 286.0;
    const auto f_135_4 = 135.0 / 4.0;
    const auto f_135_44 = 135.0 / 44.0;
    const auto f_145_22 = 145.0 / 22.0;
    const auto f_1575_16 = 1575.0 / 16.0;
    const auto f_15876_46189 = 15876.0 / 46189.0;
    const auto f_15_1 = 15.0;
    const auto f_15_143 = 15.0 / 143.0;
    const auto f_15_2 = 15.0 / 2.0;
    const auto f_18_143 = 18.0 / 143.0;
    const auto f_1_11 = 1.0 / 11.0;
    const auto f_20_11 = 20.0 / 11.0;
    const auto f_217_2717 = 217.0 / 2717.0;
    const auto f_217_286 = 217.0 / 286.0;
    const auto f_225_2 = 225.0 / 2.0;
    const auto f_238_143 = 238.0 / 143.0;
    const auto f_245_143 = 245.0 / 143.0;
    const auto f_24_11 = 24.0 / 11.0;
    const auto f_24_187 = 24.0 / 187.0;
    const auto f_25_1 = 25.0;
    const auto f_25_143 = 25.0 / 143.0;
    const auto f_270_143 = 270.0 / 143.0;
    const auto f_2835_46189 = 2835.0 / 46189.0;
    const auto f_30_11 = 30.0 / 11.0;
    const auto f_315_4 = 315.0 / 4.0;
    const auto f_315_8 = 315.0 / 8.0;
    const auto f_32_11 = 32.0 / 11.0;
    const auto f_32_187 = 32.0 / 187.0;
    const auto f_35_2717 = 35.0 / 2717.0;
    const auto f_35_286 = 35.0 / 286.0;
    const auto f_3_143 = 3.0 / 143.0;
    const auto f_405_143 = 405.0 / 143.0;
    const auto f_405_22 = 405.0 / 22.0;
    const auto f_45_1 = 45.0;
    const auto f_45_2 = 45.0 / 2.0;
    const auto f_45_8 = 45.0 / 8.0;
    const auto f_476_2717 = 476.0 / 2717.0;
    const auto f_490_2717 = 490.0 / 2717.0;
    const auto f_49_143 = 49.0 / 143.0;
    const auto f_50_11 = 50.0 / 11.0;
    const auto f_50_429 = 50.0 / 429.0;
    const auto f_511_2717 = 511.0 / 2717.0;
    const auto f_511_286 = 511.0 / 286.0;
    const auto f_525_8 = 525.0 / 8.0;
    const auto f_58_33 = 58.0 / 33.0;
    const auto f_58_561 = 58.0 / 561.0;
    const auto f_5_2 = 5.0 / 2.0;
    const auto f_5_429 = 5.0 / 429.0;
    const auto f_630_46189 = 630.0 / 46189.0;
    const auto f_63_46189 = 63.0 / 46189.0;
    const auto f_7560_46189 = 7560.0 / 46189.0;
    const auto f_75_1 = 75.0;
    const auto f_75_2 = 75.0 / 2.0;
    const auto f_75_22 = 75.0 / 22.0;
    const auto f_80_33 = 80.0 / 33.0;
    const auto f_80_561 = 80.0 / 561.0;
    const auto f_8_11 = 8.0 / 11.0;
    const auto f_8_187 = 8.0 / 187.0;
    const auto f_90_11 = 90.0 / 11.0;
    const auto f_945_16 = 945.0 / 16.0;
    const auto f_945_32 = 945.0 / 32.0;
    const auto f_98_2717 = 98.0 / 2717.0;
    const auto fs_105_16_15 = std::sqrt(165375.0 / 256.0);
    const auto fs_105_16_21 = std::sqrt(231525.0 / 256.0);
    const auto fs_105_16_5 = std::sqrt(55125.0 / 256.0);
    const auto fs_105_4_6 = std::sqrt(33075.0 / 8.0);
    const auto fs_105_8_14 = std::sqrt(77175.0 / 32.0);
    const auto fs_105_8_35 = std::sqrt(385875.0 / 64.0);
    const auto fs_10_11_11 = std::sqrt(100.0 / 11.0);
    const auto fs_10_11_35 = std::sqrt(3500.0 / 121.0);
    const auto fs_10_11_7 = std::sqrt(700.0 / 121.0);
    const auto fs_10_187_11 = std::sqrt(100.0 / 3179.0);
    const auto fs_10_187_7 = std::sqrt(700.0 / 34969.0);
    const auto fs_10_33_15 = std::sqrt(500.0 / 363.0);
    const auto fs_10_33_21 = std::sqrt(700.0 / 363.0);
    const auto fs_10_33_42 = std::sqrt(1400.0 / 363.0);
    const auto fs_10_33_5 = std::sqrt(500.0 / 1089.0);
    const auto fs_10_429_14 = std::sqrt(1400.0 / 184041.0);
    const auto fs_10_429_35 = std::sqrt(3500.0 / 184041.0);
    const auto fs_10_561_42 = std::sqrt(1400.0 / 104907.0);
    const auto fs_125_22_3 = std::sqrt(46875.0 / 484.0);
    const auto fs_126_46189_1001 = std::sqrt(111132.0 / 14919047.0);
    const auto fs_126_46189_12155 = std::sqrt(79380.0 / 877591.0);
    const auto fs_126_46189_2431 = std::sqrt(15876.0 / 877591.0);
    const auto fs_126_46189_3003 = std::sqrt(333396.0 / 14919047.0);
    const auto fs_135_143_5 = std::sqrt(91125.0 / 20449.0);
    const auto fs_135_22_5 = std::sqrt(91125.0 / 484.0);
    const auto fs_135_286_15 = std::sqrt(273375.0 / 81796.0);
    const auto fs_135_286_21 = std::sqrt(382725.0 / 81796.0);
    const auto fs_135_286_30 = std::sqrt(273375.0 / 40898.0);
    const auto fs_135_286_35 = std::sqrt(637875.0 / 81796.0);
    const auto fs_135_286_6 = std::sqrt(54675.0 / 40898.0);
    const auto fs_135_44_15 = std::sqrt(273375.0 / 1936.0);
    const auto fs_135_44_21 = std::sqrt(382725.0 / 1936.0);
    const auto fs_135_44_30 = std::sqrt(273375.0 / 968.0);
    const auto fs_135_44_35 = std::sqrt(637875.0 / 1936.0);
    const auto fs_135_44_6 = std::sqrt(54675.0 / 968.0);
    const auto fs_135_572_10 = std::sqrt(91125.0 / 163592.0);
    const auto fs_135_572_70 = std::sqrt(637875.0 / 163592.0);
    const auto fs_135_88_10 = std::sqrt(91125.0 / 3872.0);
    const auto fs_135_88_70 = std::sqrt(637875.0 / 3872.0);
    const auto fs_13_33_6 = std::sqrt(338.0 / 363.0);
    const auto fs_13_561_6 = std::sqrt(338.0 / 104907.0);
    const auto fs_14_209_2 = std::sqrt(392.0 / 43681.0);
    const auto fs_14_2717_429 = std::sqrt(588.0 / 51623.0);
    const auto fs_14_2717_858 = std::sqrt(1176.0 / 51623.0);
    const auto fs_1512_46189_77 = std::sqrt(16003008.0 / 193947611.0);
    const auto fs_15_1_14 = std::sqrt(3150.0);
    const auto fs_15_1_35 = std::sqrt(7875.0);
    const auto fs_15_22_154 = std::sqrt(1575.0 / 22.0);
    const auto fs_15_22_30 = std::sqrt(3375.0 / 242.0);
    const auto fs_15_22_55 = std::sqrt(1125.0 / 44.0);
    const auto fs_15_22_70 = std::sqrt(7875.0 / 242.0);
    const auto fs_15_2_15 = std::sqrt(3375.0 / 4.0);
    const auto fs_15_2_21 = std::sqrt(4725.0 / 4.0);
    const auto fs_15_2_5 = std::sqrt(1125.0 / 4.0);
    const auto fs_15_44_10 = std::sqrt(1125.0 / 968.0);
    const auto fs_15_44_210 = std::sqrt(23625.0 / 968.0);
    const auto fs_15_4_30 = std::sqrt(3375.0 / 8.0);
    const auto fs_189_46189_2431 = std::sqrt(35721.0 / 877591.0);
    const auto fs_189_46189_3003 = std::sqrt(750141.0 / 14919047.0);
    const auto fs_189_46189_330 = std::sqrt(1071630.0 / 193947611.0);
    const auto fs_189_46189_4290 = std::sqrt(1071630.0 / 14919047.0);
    const auto fs_189_46189_770 = std::sqrt(2500470.0 / 193947611.0);
    const auto fs_189_92378_10010 = std::sqrt(1250235.0 / 29838094.0);
    const auto fs_189_92378_110 = std::sqrt(178605.0 / 387895222.0);
    const auto fs_1_11_10 = std::sqrt(10.0 / 121.0);
    const auto fs_1_11_210 = std::sqrt(210.0 / 121.0);
    const auto fs_1_187_10 = std::sqrt(10.0 / 34969.0);
    const auto fs_1_187_210 = std::sqrt(210.0 / 34969.0);
    const auto fs_1_3_42 = std::sqrt(14.0 / 3.0);
    const auto fs_1_51_42 = std::sqrt(14.0 / 867.0);
    const auto fs_20_33_14 = std::sqrt(5600.0 / 1089.0);
    const auto fs_20_33_2 = std::sqrt(800.0 / 1089.0);
    const auto fs_20_33_35 = std::sqrt(14000.0 / 1089.0);
    const auto fs_20_561_2 = std::sqrt(800.0 / 314721.0);
    const auto fs_21_143_70 = std::sqrt(30870.0 / 20449.0);
    const auto fs_21_143_77 = std::sqrt(3087.0 / 1859.0);
    const auto fs_21_2717_10 = std::sqrt(4410.0 / 7382089.0);
    const auto fs_21_2717_165 = std::sqrt(6615.0 / 671099.0);
    const auto fs_21_2717_715 = std::sqrt(2205.0 / 51623.0);
    const auto fs_21_286_10 = std::sqrt(2205.0 / 40898.0);
    const auto fs_21_286_165 = std::sqrt(6615.0 / 7436.0);
    const auto fs_21_286_715 = std::sqrt(2205.0 / 572.0);
    const auto fs_252_46189_1430 = std::sqrt(635040.0 / 14919047.0);
    const auto fs_252_46189_2145 = std::sqrt(952560.0 / 14919047.0);
    const auto fs_252_46189_55 = std::sqrt(317520.0 / 193947611.0);
    const auto fs_25_11_2 = std::sqrt(1250.0 / 121.0);
    const auto fs_25_22_42 = std::sqrt(13125.0 / 242.0);
    const auto fs_25_2_2 = std::sqrt(625.0 / 2.0);
    const auto fs_25_2_3 = std::sqrt(1875.0 / 4.0);
    const auto fs_25_429_2 = std::sqrt(1250.0 / 184041.0);
    const auto fs_25_429_3 = std::sqrt(625.0 / 61347.0);
    const auto fs_2646_46189_33 = std::sqrt(21003948.0 / 193947611.0);
    const auto fs_28_143_55 = std::sqrt(3920.0 / 1859.0);
    const auto fs_2_11_154 = std::sqrt(56.0 / 11.0);
    const auto fs_2_11_30 = std::sqrt(120.0 / 121.0);
    const auto fs_2_11_55 = std::sqrt(20.0 / 11.0);
    const auto fs_2_11_70 = std::sqrt(280.0 / 121.0);
    const auto fs_2_187_154 = std::sqrt(56.0 / 3179.0);
    const auto fs_2_187_30 = std::sqrt(120.0 / 34969.0);
    const auto fs_2_187_55 = std::sqrt(20.0 / 3179.0);
    const auto fs_2_187_70 = std::sqrt(280.0 / 34969.0);
    const auto fs_2_33_14 = std::sqrt(56.0 / 1089.0);
    const auto fs_2_33_231 = std::sqrt(28.0 / 33.0);
    const auto fs_2_3_5 = std::sqrt(20.0 / 9.0);
    const auto fs_2_51_5 = std::sqrt(20.0 / 2601.0);
    const auto fs_2_561_14 = std::sqrt(56.0 / 314721.0);
    const auto fs_2_561_231 = std::sqrt(28.0 / 9537.0);
    const auto fs_30_11_5 = std::sqrt(4500.0 / 121.0);
    const auto fs_30_11_7 = std::sqrt(6300.0 / 121.0);
    const auto fs_315_32_30 = std::sqrt(1488375.0 / 512.0);
    const auto fs_315_92378_2002 = std::sqrt(694575.0 / 29838094.0);
    const auto fs_329_5434_6 = std::sqrt(324723.0 / 14764178.0);
    const auto fs_329_572_6 = std::sqrt(324723.0 / 163592.0);
    const auto fs_32_33_3 = std::sqrt(1024.0 / 363.0);
    const auto fs_32_561_3 = std::sqrt(1024.0 / 104907.0);
    const auto fs_35_2717_14 = std::sqrt(17150.0 / 7382089.0);
    const auto fs_35_2717_143 = std::sqrt(1225.0 / 51623.0);
    const auto fs_35_2717_286 = std::sqrt(2450.0 / 51623.0);
    const auto fs_35_2717_6 = std::sqrt(7350.0 / 7382089.0);
    const auto fs_35_2717_66 = std::sqrt(7350.0 / 671099.0);
    const auto fs_35_286_14 = std::sqrt(8575.0 / 40898.0);
    const auto fs_35_286_143 = std::sqrt(1225.0 / 572.0);
    const auto fs_35_286_286 = std::sqrt(1225.0 / 286.0);
    const auto fs_35_286_6 = std::sqrt(3675.0 / 40898.0);
    const auto fs_35_286_66 = std::sqrt(3675.0 / 3718.0);
    const auto fs_35_33_6 = std::sqrt(2450.0 / 363.0);
    const auto fs_35_4_6 = std::sqrt(3675.0 / 8.0);
    const auto fs_35_5434_154 = std::sqrt(8575.0 / 1342198.0);
    const auto fs_35_5434_330 = std::sqrt(18375.0 / 1342198.0);
    const auto fs_35_5434_858 = std::sqrt(3675.0 / 103246.0);
    const auto fs_35_572_154 = std::sqrt(8575.0 / 14872.0);
    const auto fs_35_572_330 = std::sqrt(18375.0 / 14872.0);
    const auto fs_35_572_858 = std::sqrt(3675.0 / 1144.0);
    const auto fs_35_858_6 = std::sqrt(1225.0 / 122694.0);
    const auto fs_378_46189_385 = std::sqrt(5000940.0 / 193947611.0);
    const auto fs_3_11_66 = std::sqrt(54.0 / 11.0);
    const auto fs_3_143_15 = std::sqrt(135.0 / 20449.0);
    const auto fs_3_143_21 = std::sqrt(189.0 / 20449.0);
    const auto fs_3_143_30 = std::sqrt(270.0 / 20449.0);
    const auto fs_3_143_35 = std::sqrt(315.0 / 20449.0);
    const auto fs_3_143_6 = std::sqrt(54.0 / 20449.0);
    const auto fs_3_187_66 = std::sqrt(54.0 / 3179.0);
    const auto fs_3_286_10 = std::sqrt(45.0 / 40898.0);
    const auto fs_3_286_70 = std::sqrt(315.0 / 40898.0);
    const auto fs_40_11_3 = std::sqrt(4800.0 / 121.0);
    const auto fs_42_2717_70 = std::sqrt(123480.0 / 7382089.0);
    const auto fs_42_2717_77 = std::sqrt(12348.0 / 671099.0);
    const auto fs_45_16_10 = std::sqrt(10125.0 / 128.0);
    const auto fs_45_16_70 = std::sqrt(70875.0 / 128.0);
    const auto fs_45_44_66 = std::sqrt(6075.0 / 88.0);
    const auto fs_45_4_30 = std::sqrt(30375.0 / 8.0);
    const auto fs_45_4_5 = std::sqrt(10125.0 / 16.0);
    const auto fs_45_8_15 = std::sqrt(30375.0 / 64.0);
    const auto fs_45_8_21 = std::sqrt(42525.0 / 64.0);
    const auto fs_45_8_30 = std::sqrt(30375.0 / 32.0);
    const auto fs_45_8_35 = std::sqrt(70875.0 / 64.0);
    const auto fs_45_8_6 = std::sqrt(6075.0 / 32.0);
    const auto fs_49_143_15 = std::sqrt(36015.0 / 20449.0);
    const auto fs_49_5434_210 = std::sqrt(252105.0 / 14764178.0);
    const auto fs_49_572_210 = std::sqrt(252105.0 / 163592.0);
    const auto fs_4_33_210 = std::sqrt(1120.0 / 363.0);
    const auto fs_4_33_462 = std::sqrt(224.0 / 33.0);
    const auto fs_4_561_210 = std::sqrt(1120.0 / 104907.0);
    const auto fs_4_561_462 = std::sqrt(224.0 / 9537.0);
    const auto fs_50_33_2 = std::sqrt(5000.0 / 1089.0);
    const auto fs_50_33_3 = std::sqrt(2500.0 / 363.0);
    const auto fs_50_561_3 = std::sqrt(2500.0 / 104907.0);
    const auto fs_525_16_2 = std::sqrt(275625.0 / 128.0);
    const auto fs_525_16_3 = std::sqrt(826875.0 / 256.0);
    const auto fs_56_2717_55 = std::sqrt(15680.0 / 671099.0);
    const auto fs_5_11_210 = std::sqrt(5250.0 / 121.0);
    const auto fs_5_11_30 = std::sqrt(750.0 / 121.0);
    const auto fs_5_11_462 = std::sqrt(1050.0 / 11.0);
    const auto fs_5_187_30 = std::sqrt(750.0 / 34969.0);
    const auto fs_5_1_14 = std::sqrt(350.0);
    const auto fs_5_1_35 = std::sqrt(875.0);
    const auto fs_5_22_14 = std::sqrt(175.0 / 242.0);
    const auto fs_5_22_231 = std::sqrt(525.0 / 44.0);
    const auto fs_5_286_30 = std::sqrt(375.0 / 40898.0);
    const auto fs_5_2_15 = std::sqrt(375.0 / 4.0);
    const auto fs_5_2_21 = std::sqrt(525.0 / 4.0);
    const auto fs_5_2_5 = std::sqrt(125.0 / 4.0);
    const auto fs_5_429_15 = std::sqrt(125.0 / 61347.0);
    const auto fs_5_429_21 = std::sqrt(175.0 / 61347.0);
    const auto fs_5_429_5 = std::sqrt(125.0 / 184041.0);
    const auto fs_5_4_42 = std::sqrt(525.0 / 8.0);
    const auto fs_63_2717_7 = std::sqrt(27783.0 / 7382089.0);
    const auto fs_63_286_7 = std::sqrt(27783.0 / 81796.0);
    const auto fs_63_46189_143 = std::sqrt(3969.0 / 14919047.0);
    const auto fs_63_46189_3003 = std::sqrt(83349.0 / 14919047.0);
    const auto fs_63_46189_30030 = std::sqrt(833490.0 / 14919047.0);
    const auto fs_63_46189_33 = std::sqrt(11907.0 / 193947611.0);
    const auto fs_63_46189_36465 = std::sqrt(59535.0 / 877591.0);
    const auto fs_63_46189_46189 = std::sqrt(3969.0 / 46189.0);
    const auto fs_63_46189_92378 = std::sqrt(7938.0 / 46189.0);
    const auto fs_63_92378_10010 = std::sqrt(138915.0 / 29838094.0);
    const auto fs_63_92378_2002 = std::sqrt(27783.0 / 29838094.0);
    const auto fs_63_92378_22 = std::sqrt(3969.0 / 387895222.0);
    const auto fs_65_44_6 = std::sqrt(12675.0 / 968.0);
    const auto fs_6_143_5 = std::sqrt(180.0 / 20449.0);
    const auto fs_735_32_6 = std::sqrt(1620675.0 / 512.0);
    const auto fs_75_22_11 = std::sqrt(5625.0 / 44.0);
    const auto fs_75_22_7 = std::sqrt(39375.0 / 484.0);
    const auto fs_75_2_2 = std::sqrt(5625.0 / 2.0);
    const auto fs_75_2_3 = std::sqrt(16875.0 / 4.0);
    const auto fs_75_44_30 = std::sqrt(84375.0 / 968.0);
    const auto fs_7_11_2 = std::sqrt(98.0 / 121.0);
    const auto fs_7_143_429 = std::sqrt(147.0 / 143.0);
    const auto fs_7_143_858 = std::sqrt(294.0 / 143.0);
    const auto fs_7_209_15 = std::sqrt(735.0 / 43681.0);
    const auto fs_7_22_15 = std::sqrt(735.0 / 484.0);
    const auto fs_7_247_55 = std::sqrt(2695.0 / 61009.0);
    const auto fs_7_26_55 = std::sqrt(2695.0 / 676.0);
    const auto fs_7_2717_1430 = std::sqrt(490.0 / 51623.0);
    const auto fs_7_2717_5005 = std::sqrt(1715.0 / 51623.0);
    const auto fs_7_286_1430 = std::sqrt(245.0 / 286.0);
    const auto fs_7_286_5005 = std::sqrt(1715.0 / 572.0);
    const auto fs_7_5434_286 = std::sqrt(49.0 / 103246.0);
    const auto fs_7_5434_6006 = std::sqrt(1029.0 / 103246.0);
    const auto fs_7_5434_66 = std::sqrt(147.0 / 1342198.0);
    const auto fs_7_572_286 = std::sqrt(49.0 / 1144.0);
    const auto fs_7_572_6006 = std::sqrt(1029.0 / 1144.0);
    const auto fs_7_572_66 = std::sqrt(147.0 / 14872.0);
    const auto fs_882_46189_165 = std::sqrt(11668860.0 / 193947611.0);
    const auto fs_8_11_5 = std::sqrt(320.0 / 121.0);
    const auto fs_8_11_7 = std::sqrt(448.0 / 121.0);
    const auto fs_8_187_5 = std::sqrt(320.0 / 34969.0);
    const auto fs_8_187_7 = std::sqrt(448.0 / 34969.0);
    const auto fs_8_33_35 = std::sqrt(2240.0 / 1089.0);
    const auto fs_8_561_35 = std::sqrt(2240.0 / 314721.0);
    const auto fs_945_46189_143 = std::sqrt(893025.0 / 14919047.0);
    const auto fs_98_2717_15 = std::sqrt(144060.0 / 7382089.0);

    // NOTE: the angular components are formed in 53 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph2_0, ph4_0, ph6_0, ph8_0, ph10_0, ph10_p10, ab_2, pc_0 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
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
        const auto h10_0 = ph10_0[k];
        const auto h10_p10 = ph10_p10[k];

        pc_0[k] = e_0 * f_945_32 + e_1 * f_1575_16 * h2_0 - e_1 * f_1575_16 * r_2 + e_2 * f_135_4 * h4_0 - e_2 * f_225_2 * r_2 * h2_0 + e_2 * f_315_4 * r_4 + e_3 * f_75_22 * h6_0 - e_3 * f_405_22 * r_2 * h4_0 + e_3 * f_75_2 * r_4 * h2_0 - e_3 * f_45_2 * r_6 + e_4 * f_35_286 * h8_0 - e_4 * f_10_11 * r_2 * h6_0 + e_4 * f_405_143 * r_4 * h4_0 - e_4 * f_50_11 * r_6 * h2_0 + e_4 * f_5_2 * r_8 + e_5 * f_63_46189 * h10_0 + e_5 * fs_63_46189_92378 * h10_p10 - e_5 * f_35_2717 * r_2 * h8_0 + e_5 * f_10_187 * r_4 * h6_0 - e_5 * f_18_143 * r_6 * h4_0 + e_5 * f_25_143 * r_8 * h2_0 - e_5 * f_1_11 * r_10;
    }

    // NOTE: the angular components are formed in 53 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph2_p1, ph4_p1, ph6_p1, ph8_p1, ph10_p1, ph10_p9, ab_2, pc_1 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
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

        const auto h2_p1 = ph2_p1[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h10_p1 = ph10_p1[k];
        const auto h10_p9 = ph10_p9[k];

        pc_1[k] = - e_1 * fs_315_32_30 * h2_p1 - e_2 * f_135_4 * h4_p1 + e_2 * fs_45_4_30 * r_2 * h2_p1 - e_3 * fs_15_44_210 * h6_p1 + e_3 * f_405_22 * r_2 * h4_p1 - e_3 * fs_15_4_30 * r_4 * h2_p1 - e_4 * fs_21_286_10 * h8_p1 + e_4 * fs_1_11_210 * r_2 * h6_p1 - e_4 * f_405_143 * r_4 * h4_p1 + e_4 * fs_5_11_30 * r_6 * h2_p1 - e_5 * fs_63_92378_22 * h10_p1 + e_5 * fs_63_46189_46189 * h10_p9 + e_5 * fs_21_2717_10 * r_2 * h8_p1 - e_5 * fs_1_187_210 * r_4 * h6_p1 + e_5 * f_18_143 * r_6 * h4_p1 - e_5 * fs_5_286_30 * r_8 * h2_p1;
    }

    // NOTE: the angular components are formed in 53 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph2_p2, ph4_p2, ph6_p2, ph8_p2, ph8_p8, ph10_p2, ph10_p8, ab_2, pc_2 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
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

        const auto h2_p2 = ph2_p2[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h8_p2 = ph8_p2[k];
        const auto h8_p8 = ph8_p8[k];
        const auto h10_p2 = ph10_p2[k];
        const auto h10_p8 = ph10_p8[k];

        pc_2[k] = e_1 * fs_105_16_15 * h2_p2 + e_2 * f_135_4 * h4_p2 - e_2 * fs_15_2_15 * r_2 * h2_p2 + e_3 * fs_25_22_42 * h6_p2 - e_3 * f_405_22 * r_2 * h4_p2 + e_3 * fs_5_2_15 * r_4 * h2_p2 + e_4 * fs_35_286_14 * h8_p2 + e_4 * fs_35_286_143 * h8_p8 - e_4 * fs_10_33_42 * r_2 * h6_p2 + e_4 * f_405_143 * r_4 * h4_p2 - e_4 * fs_10_33_15 * r_6 * h2_p2 + e_5 * fs_63_46189_33 * h10_p2 + e_5 * fs_189_46189_2431 * h10_p8 - e_5 * fs_35_2717_14 * r_2 * h8_p2 - e_5 * fs_35_2717_143 * r_2 * h8_p8 + e_5 * fs_10_561_42 * r_4 * h6_p2 - e_5 * f_18_143 * r_6 * h4_p2 + e_5 * fs_5_429_15 * r_8 * h2_p2;
    }

    // NOTE: the angular components are formed in 53 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_2, pe_3, pe_4, pe_5, ph4_p3, ph4_p4, ph6_p3, ph6_p4, ph6_p6, ph8_p3, ph8_p4, ph8_p6, ph8_p7, ph10_p3, ph10_p4, ph10_p6, ph10_p7, ab_2, pc_3, pc_4 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

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

        pc_3[k] = - e_2 * fs_45_8_21 * h4_p3 - e_3 * fs_75_22_7 * h6_p3 + e_3 * fs_135_44_21 * r_2 * h4_p3 - e_4 * fs_35_572_154 * h8_p3 + e_4 * fs_35_572_858 * h8_p7 + e_4 * fs_10_11_7 * r_2 * h6_p3 - e_4 * fs_135_286_21 * r_4 * h4_p3 - e_5 * fs_63_46189_143 * h10_p3 + e_5 * fs_126_46189_2431 * h10_p7 + e_5 * fs_35_5434_154 * r_2 * h8_p3 - e_5 * fs_35_5434_858 * r_2 * h8_p7 - e_5 * fs_10_187_7 * r_4 * h6_p3 + e_5 * fs_3_143_21 * r_6 * h4_p3;

        pc_4[k] = e_2 * fs_45_8_6 * h4_p4 + e_3 * fs_75_44_30 * h6_p4 + e_3 * fs_15_22_55 * h6_p6 - e_3 * fs_135_44_6 * r_2 * h4_p4 + e_4 * fs_35_572_330 * h8_p4 + e_4 * fs_7_286_5005 * h8_p6 - e_4 * fs_5_11_30 * r_2 * h6_p4 - e_4 * fs_2_11_55 * r_2 * h6_p6 + e_4 * fs_135_286_6 * r_4 * h4_p4 + e_5 * fs_63_92378_2002 * h10_p4 + e_5 * fs_126_46189_1001 * h10_p6 - e_5 * fs_35_5434_330 * r_2 * h8_p4 - e_5 * fs_7_2717_5005 * r_2 * h8_p6 + e_5 * fs_5_187_30 * r_4 * h6_p4 + e_5 * fs_2_187_55 * r_4 * h6_p6 - e_5 * fs_3_143_6 * r_6 * h4_p4;
    }

    // NOTE: the angular components are formed in 53 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_2, pe_3, pe_4, pe_5, ph4_m4, ph6_m6, ph6_m5, ph6_m4, ph8_m6, ph8_m5, ph8_m4, ph10_m6, ph10_m5, ph10_m4, ab_2, pc_5, pc_6 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

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

        pc_5[k] = - e_3 * fs_75_22_11 * h6_m5 - e_4 * fs_35_286_286 * h8_m5 + e_4 * fs_10_11_11 * r_2 * h6_m5 - e_5 * fs_63_46189_3003 * h10_m5 + e_5 * fs_35_2717_286 * r_2 * h8_m5 - e_5 * fs_10_187_11 * r_4 * h6_m5;

        pc_6[k] = e_2 * fs_45_8_6 * h4_m4 - e_3 * fs_15_22_55 * h6_m6 + e_3 * fs_75_44_30 * h6_m4 - e_3 * fs_135_44_6 * r_2 * h4_m4 - e_4 * fs_7_286_5005 * h8_m6 + e_4 * fs_35_572_330 * h8_m4 + e_4 * fs_2_11_55 * r_2 * h6_m6 - e_4 * fs_5_11_30 * r_2 * h6_m4 + e_4 * fs_135_286_6 * r_4 * h4_m4 - e_5 * fs_126_46189_1001 * h10_m6 + e_5 * fs_63_92378_2002 * h10_m4 + e_5 * fs_7_2717_5005 * r_2 * h8_m6 - e_5 * fs_35_5434_330 * r_2 * h8_m4 - e_5 * fs_2_187_55 * r_4 * h6_m6 + e_5 * fs_5_187_30 * r_4 * h6_m4 - e_5 * fs_3_143_6 * r_6 * h4_m4;
    }

    // NOTE: the angular components are formed in 53 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_2, pe_3, pe_4, pe_5, ph4_m3, ph6_m3, ph8_m7, ph8_m3, ph10_m7, ph10_m3, ab_2, pc_7 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h4_m3 = ph4_m3[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h8_m7 = ph8_m7[k];
        const auto h8_m3 = ph8_m3[k];
        const auto h10_m7 = ph10_m7[k];
        const auto h10_m3 = ph10_m3[k];

        pc_7[k] = - e_2 * fs_45_8_21 * h4_m3 - e_3 * fs_75_22_7 * h6_m3 + e_3 * fs_135_44_21 * r_2 * h4_m3 - e_4 * fs_35_572_858 * h8_m7 - e_4 * fs_35_572_154 * h8_m3 + e_4 * fs_10_11_7 * r_2 * h6_m3 - e_4 * fs_135_286_21 * r_4 * h4_m3 - e_5 * fs_126_46189_2431 * h10_m7 - e_5 * fs_63_46189_143 * h10_m3 + e_5 * fs_35_5434_858 * r_2 * h8_m7 + e_5 * fs_35_5434_154 * r_2 * h8_m3 - e_5 * fs_10_187_7 * r_4 * h6_m3 + e_5 * fs_3_143_21 * r_6 * h4_m3;
    }

    // NOTE: the angular components are formed in 53 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph2_m2, ph4_m2, ph6_m2, ph8_m8, ph8_m2, ph10_m8, ph10_m2, ab_2, pc_8 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
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

        const auto h2_m2 = ph2_m2[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h8_m8 = ph8_m8[k];
        const auto h8_m2 = ph8_m2[k];
        const auto h10_m8 = ph10_m8[k];
        const auto h10_m2 = ph10_m2[k];

        pc_8[k] = e_1 * fs_105_16_15 * h2_m2 + e_2 * f_135_4 * h4_m2 - e_2 * fs_15_2_15 * r_2 * h2_m2 + e_3 * fs_25_22_42 * h6_m2 - e_3 * f_405_22 * r_2 * h4_m2 + e_3 * fs_5_2_15 * r_4 * h2_m2 - e_4 * fs_35_286_143 * h8_m8 + e_4 * fs_35_286_14 * h8_m2 - e_4 * fs_10_33_42 * r_2 * h6_m2 + e_4 * f_405_143 * r_4 * h4_m2 - e_4 * fs_10_33_15 * r_6 * h2_m2 - e_5 * fs_189_46189_2431 * h10_m8 + e_5 * fs_63_46189_33 * h10_m2 + e_5 * fs_35_2717_143 * r_2 * h8_m8 - e_5 * fs_35_2717_14 * r_2 * h8_m2 + e_5 * fs_10_561_42 * r_4 * h6_m2 - e_5 * f_18_143 * r_6 * h4_m2 + e_5 * fs_5_429_15 * r_8 * h2_m2;
    }

    // NOTE: the angular components are formed in 53 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph2_m1, ph4_m1, ph6_m1, ph8_m1, ph10_m10, ph10_m9, ph10_m1, ab_2, pc_9, pc_10 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
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

        const auto h2_m1 = ph2_m1[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h8_m1 = ph8_m1[k];
        const auto h10_m10 = ph10_m10[k];
        const auto h10_m9 = ph10_m9[k];
        const auto h10_m1 = ph10_m1[k];

        pc_9[k] = - e_1 * fs_315_32_30 * h2_m1 - e_2 * f_135_4 * h4_m1 + e_2 * fs_45_4_30 * r_2 * h2_m1 - e_3 * fs_15_44_210 * h6_m1 + e_3 * f_405_22 * r_2 * h4_m1 - e_3 * fs_15_4_30 * r_4 * h2_m1 - e_4 * fs_21_286_10 * h8_m1 + e_4 * fs_1_11_210 * r_2 * h6_m1 - e_4 * f_405_143 * r_4 * h4_m1 + e_4 * fs_5_11_30 * r_6 * h2_m1 - e_5 * fs_63_46189_46189 * h10_m9 - e_5 * fs_63_92378_22 * h10_m1 + e_5 * fs_21_2717_10 * r_2 * h8_m1 - e_5 * fs_1_187_210 * r_4 * h6_m1 + e_5 * f_18_143 * r_6 * h4_m1 - e_5 * fs_5_286_30 * r_8 * h2_m1;

        pc_10[k] = - e_5 * fs_63_46189_92378 * h10_m10;
    }

    // NOTE: the angular components are formed in 53 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph2_0, ph4_0, ph6_0, ph8_0, ph8_p8, ph10_0, ph10_p8, ab_2, pc_11 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
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

        pc_11[k] = e_0 * f_945_32 + e_1 * f_315_8 * h2_0 - e_1 * f_1575_16 * r_2 - e_2 * f_135_4 * h4_0 - e_2 * f_45_1 * r_2 * h2_0 + e_2 * f_315_4 * r_4 - e_3 * f_120_11 * h6_0 + e_3 * f_405_22 * r_2 * h4_0 + e_3 * f_15_1 * r_4 * h2_0 - e_3 * f_45_2 * r_6 - e_4 * f_217_286 * h8_0 - e_4 * fs_21_286_715 * h8_p8 + e_4 * f_32_11 * r_2 * h6_0 - e_4 * f_405_143 * r_4 * h4_0 - e_4 * f_20_11 * r_6 * h2_0 + e_4 * f_5_2 * r_8 - e_5 * f_630_46189 * h10_0 + e_5 * fs_126_46189_12155 * h10_p8 + e_5 * f_217_2717 * r_2 * h8_0 + e_5 * fs_21_2717_715 * r_2 * h8_p8 - e_5 * f_32_187 * r_4 * h6_0 + e_5 * f_18_143 * r_6 * h4_0 + e_5 * f_10_143 * r_8 * h2_0 - e_5 * f_1_11 * r_10;
    }

    // NOTE: the angular components are formed in 53 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph2_p1, ph6_p1, ph8_p1, ph8_p7, ph10_p1, ph10_p7, ab_2, pc_12 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
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

        const auto h2_p1 = ph2_p1[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h8_p7 = ph8_p7[k];
        const auto h10_p1 = ph10_p1[k];
        const auto h10_p7 = ph10_p7[k];

        pc_12[k] = - e_1 * fs_735_32_6 * h2_p1 + e_2 * fs_105_4_6 * r_2 * h2_p1 + e_3 * fs_5_4_42 * h6_p1 - e_3 * fs_35_4_6 * r_4 * h2_p1 + e_4 * fs_7_11_2 * h8_p1 - e_4 * fs_7_286_1430 * h8_p7 - e_4 * fs_1_3_42 * r_2 * h6_p1 + e_4 * fs_35_33_6 * r_6 * h2_p1 + e_5 * fs_189_92378_110 * h10_p1 + e_5 * fs_63_46189_36465 * h10_p7 - e_5 * fs_14_209_2 * r_2 * h8_p1 + e_5 * fs_7_2717_1430 * r_2 * h8_p7 + e_5 * fs_1_51_42 * r_4 * h6_p1 - e_5 * fs_35_858_6 * r_8 * h2_p1;
    }

    // NOTE: the angular components are formed in 53 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph2_p2, ph4_p2, ph6_p2, ph6_p6, ph8_p2, ph8_p6, ph10_p2, ph10_p6, ab_2, pc_13 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
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

        const auto h2_p2 = ph2_p2[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h6_p6 = ph6_p6[k];
        const auto h8_p2 = ph8_p2[k];
        const auto h8_p6 = ph8_p6[k];
        const auto h10_p2 = ph10_p2[k];
        const auto h10_p6 = ph10_p6[k];

        pc_13[k] = e_1 * f_315_8 * h2_p2 + e_2 * fs_45_8_15 * h4_p2 - e_2 * f_45_1 * r_2 * h2_p2 - e_3 * fs_15_22_70 * h6_p2 - e_3 * fs_15_22_154 * h6_p6 - e_3 * fs_135_44_15 * r_2 * h4_p2 + e_3 * f_15_1 * r_4 * h2_p2 - e_4 * fs_49_572_210 * h8_p2 + e_4 * fs_7_572_286 * h8_p6 + e_4 * fs_2_11_70 * r_2 * h6_p2 + e_4 * fs_2_11_154 * r_2 * h6_p6 + e_4 * fs_135_286_15 * r_4 * h4_p2 - e_4 * f_20_11 * r_6 * h2_p2 - e_5 * fs_252_46189_55 * h10_p2 + e_5 * fs_252_46189_1430 * h10_p6 + e_5 * fs_49_5434_210 * r_2 * h8_p2 - e_5 * fs_7_5434_286 * r_2 * h8_p6 - e_5 * fs_2_187_70 * r_4 * h6_p2 - e_5 * fs_2_187_154 * r_4 * h6_p6 - e_5 * fs_3_143_15 * r_6 * h4_p2 + e_5 * f_10_143 * r_8 * h2_p2;
    }

    // NOTE: the angular components are formed in 53 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_2, pe_3, pe_4, pe_5, ph4_m4, ph4_p3, ph6_m4, ph6_p3, ph6_p5, ph8_m4, ph8_p3, ph8_p5, ph10_m4, ph10_p3, ph10_p5, ab_2, pc_14, pc_15 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

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

        pc_14[k] = - e_2 * fs_45_8_30 * h4_p3 + e_3 * fs_15_44_10 * h6_p3 - e_3 * fs_45_44_66 * h6_p5 + e_3 * fs_135_44_30 * r_2 * h4_p3 + e_4 * fs_28_143_55 * h8_p3 + e_4 * fs_7_143_429 * h8_p5 - e_4 * fs_1_11_10 * r_2 * h6_p3 + e_4 * fs_3_11_66 * r_2 * h6_p5 - e_4 * fs_135_286_30 * r_4 * h4_p3 + e_5 * fs_63_92378_10010 * h10_p3 + e_5 * fs_315_92378_2002 * h10_p5 - e_5 * fs_56_2717_55 * r_2 * h8_p3 - e_5 * fs_14_2717_429 * r_2 * h8_p5 + e_5 * fs_1_187_10 * r_4 * h6_p3 - e_5 * fs_3_187_66 * r_4 * h6_p5 + e_5 * fs_3_143_30 * r_6 * h4_p3;

        pc_15[k] = e_2 * f_135_4 * h4_m4 + e_3 * fs_30_11_5 * h6_m4 - e_3 * f_405_22 * r_2 * h4_m4 - e_4 * fs_7_26_55 * h8_m4 - e_4 * fs_8_11_5 * r_2 * h6_m4 + e_4 * f_405_143 * r_4 * h4_m4 - e_5 * fs_126_46189_3003 * h10_m4 + e_5 * fs_7_247_55 * r_2 * h8_m4 + e_5 * fs_8_187_5 * r_4 * h6_m4 - e_5 * f_18_143 * r_6 * h4_m4;
    }

    // NOTE: the angular components are formed in 53 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_2, pe_3, pe_4, pe_5, ph4_m3, ph6_m5, ph6_m3, ph8_m5, ph8_m3, ph10_m5, ph10_m3, ab_2, pc_16 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

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

        pc_16[k] = - e_2 * fs_45_8_30 * h4_m3 + e_3 * fs_45_44_66 * h6_m5 + e_3 * fs_15_44_10 * h6_m3 + e_3 * fs_135_44_30 * r_2 * h4_m3 - e_4 * fs_7_143_429 * h8_m5 + e_4 * fs_28_143_55 * h8_m3 - e_4 * fs_3_11_66 * r_2 * h6_m5 - e_4 * fs_1_11_10 * r_2 * h6_m3 - e_4 * fs_135_286_30 * r_4 * h4_m3 - e_5 * fs_315_92378_2002 * h10_m5 + e_5 * fs_63_92378_10010 * h10_m3 + e_5 * fs_14_2717_429 * r_2 * h8_m5 - e_5 * fs_56_2717_55 * r_2 * h8_m3 + e_5 * fs_3_187_66 * r_4 * h6_m5 + e_5 * fs_1_187_10 * r_4 * h6_m3 + e_5 * fs_3_143_30 * r_6 * h4_m3;
    }

    // NOTE: the angular components are formed in 53 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph2_m2, ph4_m2, ph6_m6, ph6_m2, ph8_m6, ph8_m2, ph10_m6, ph10_m2, ab_2, pc_17 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
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

        const auto h2_m2 = ph2_m2[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h6_m6 = ph6_m6[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h8_m6 = ph8_m6[k];
        const auto h8_m2 = ph8_m2[k];
        const auto h10_m6 = ph10_m6[k];
        const auto h10_m2 = ph10_m2[k];

        pc_17[k] = e_1 * f_315_8 * h2_m2 + e_2 * fs_45_8_15 * h4_m2 - e_2 * f_45_1 * r_2 * h2_m2 + e_3 * fs_15_22_154 * h6_m6 - e_3 * fs_15_22_70 * h6_m2 - e_3 * fs_135_44_15 * r_2 * h4_m2 + e_3 * f_15_1 * r_4 * h2_m2 - e_4 * fs_7_572_286 * h8_m6 - e_4 * fs_49_572_210 * h8_m2 - e_4 * fs_2_11_154 * r_2 * h6_m6 + e_4 * fs_2_11_70 * r_2 * h6_m2 + e_4 * fs_135_286_15 * r_4 * h4_m2 - e_4 * f_20_11 * r_6 * h2_m2 - e_5 * fs_252_46189_1430 * h10_m6 - e_5 * fs_252_46189_55 * h10_m2 + e_5 * fs_7_5434_286 * r_2 * h8_m6 + e_5 * fs_49_5434_210 * r_2 * h8_m2 + e_5 * fs_2_187_154 * r_4 * h6_m6 - e_5 * fs_2_187_70 * r_4 * h6_m2 - e_5 * fs_3_143_15 * r_6 * h4_m2 + e_5 * f_10_143 * r_8 * h2_m2;
    }

    // NOTE: the angular components are formed in 53 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph2_m1, ph4_m1, ph6_m1, ph8_m8, ph8_m7, ph8_m1, ph10_m9, ph10_m8, ph10_m7, ph10_m1, ab_2, pc_18, pc_19, pc_20 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
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

        pc_18[k] = - e_1 * fs_735_32_6 * h2_m1 + e_2 * fs_105_4_6 * r_2 * h2_m1 + e_3 * fs_5_4_42 * h6_m1 - e_3 * fs_35_4_6 * r_4 * h2_m1 + e_4 * fs_7_286_1430 * h8_m7 + e_4 * fs_7_11_2 * h8_m1 - e_4 * fs_1_3_42 * r_2 * h6_m1 + e_4 * fs_35_33_6 * r_6 * h2_m1 - e_5 * fs_63_46189_36465 * h10_m7 + e_5 * fs_189_92378_110 * h10_m1 - e_5 * fs_7_2717_1430 * r_2 * h8_m7 - e_5 * fs_14_209_2 * r_2 * h8_m1 + e_5 * fs_1_51_42 * r_4 * h6_m1 - e_5 * fs_35_858_6 * r_8 * h2_m1;

        pc_19[k] = e_4 * fs_21_286_715 * h8_m8 - e_5 * fs_126_46189_12155 * h10_m8 - e_5 * fs_21_2717_715 * r_2 * h8_m8;

        pc_20[k] = e_1 * fs_315_32_30 * h2_m1 + e_2 * f_135_4 * h4_m1 - e_2 * fs_45_4_30 * r_2 * h2_m1 + e_3 * fs_15_44_210 * h6_m1 - e_3 * f_405_22 * r_2 * h4_m1 + e_3 * fs_15_4_30 * r_4 * h2_m1 + e_4 * fs_21_286_10 * h8_m1 - e_4 * fs_1_11_210 * r_2 * h6_m1 + e_4 * f_405_143 * r_4 * h4_m1 - e_4 * fs_5_11_30 * r_6 * h2_m1 - e_5 * fs_63_46189_46189 * h10_m9 + e_5 * fs_63_92378_22 * h10_m1 - e_5 * fs_21_2717_10 * r_2 * h8_m1 + e_5 * fs_1_187_210 * r_4 * h6_m1 - e_5 * f_18_143 * r_6 * h4_m1 + e_5 * fs_5_286_30 * r_8 * h2_m1;
    }

    // NOTE: the angular components are formed in 53 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph2_0, ph4_0, ph6_0, ph6_p6, ph8_0, ph8_p6, ph10_0, ph10_p6, ab_2, pc_21 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
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

        pc_21[k] = e_0 * f_945_32 - e_1 * f_105_16 * h2_0 - e_1 * f_1575_16 * r_2 - e_2 * f_135_4 * h4_0 + e_2 * f_15_2 * r_2 * h2_0 + e_2 * f_315_4 * r_4 + e_3 * f_145_22 * h6_0 + e_3 * fs_5_11_462 * h6_p6 + e_3 * f_405_22 * r_2 * h4_0 - e_3 * f_5_2 * r_4 * h2_0 - e_3 * f_45_2 * r_6 + e_4 * f_511_286 * h8_0 - e_4 * fs_7_143_858 * h8_p6 - e_4 * f_58_33 * r_2 * h6_0 - e_4 * fs_4_33_462 * r_2 * h6_p6 - e_4 * f_405_143 * r_4 * h4_0 + e_4 * f_10_33 * r_6 * h2_0 + e_4 * f_5_2 * r_8 + e_5 * f_2835_46189 * h10_0 + e_5 * fs_189_46189_4290 * h10_p6 - e_5 * f_511_2717 * r_2 * h8_0 + e_5 * fs_14_2717_858 * r_2 * h8_p6 + e_5 * f_58_561 * r_4 * h6_0 + e_5 * fs_4_561_462 * r_4 * h6_p6 + e_5 * f_18_143 * r_6 * h4_0 - e_5 * f_5_429 * r_8 * h2_0 - e_5 * f_1_11 * r_10;
    }

    // NOTE: the angular components are formed in 53 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph2_p1, ph4_p1, ph6_p1, ph6_p5, ph8_p1, ph8_p5, ph10_p1, ph10_p5, ab_2, pc_22 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
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

        const auto h2_p1 = ph2_p1[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h6_p5 = ph6_p5[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h8_p5 = ph8_p5[k];
        const auto h10_p1 = ph10_p1[k];
        const auto h10_p5 = ph10_p5[k];

        pc_22[k] = - e_1 * fs_525_16_2 * h2_p1 + e_2 * fs_45_8_15 * h4_p1 + e_2 * fs_75_2_2 * r_2 * h2_p1 + e_3 * fs_5_22_14 * h6_p1 + e_3 * fs_5_22_231 * h6_p5 - e_3 * fs_135_44_15 * r_2 * h4_p1 - e_3 * fs_25_2_2 * r_4 * h2_p1 - e_4 * fs_329_572_6 * h8_p1 - e_4 * fs_7_572_6006 * h8_p5 - e_4 * fs_2_33_14 * r_2 * h6_p1 - e_4 * fs_2_33_231 * r_2 * h6_p5 + e_4 * fs_135_286_15 * r_4 * h4_p1 + e_4 * fs_50_33_2 * r_6 * h2_p1 - e_5 * fs_189_46189_330 * h10_p1 + e_5 * fs_945_46189_143 * h10_p5 + e_5 * fs_329_5434_6 * r_2 * h8_p1 + e_5 * fs_7_5434_6006 * r_2 * h8_p5 + e_5 * fs_2_561_14 * r_4 * h6_p1 + e_5 * fs_2_561_231 * r_4 * h6_p5 - e_5 * fs_3_143_15 * r_6 * h4_p1 - e_5 * fs_25_429_2 * r_8 * h2_p1;
    }

    // NOTE: the angular components are formed in 53 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph2_p2, ph4_m3, ph4_p4, ph6_m3, ph6_p2, ph6_p4, ph8_m3, ph8_p2, ph8_p4, ph10_m3, ph10_p2, ph10_p4, ab_2, pc_23, pc_24 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
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

        pc_23[k] = e_1 * fs_105_8_14 * h2_p2 + e_2 * fs_45_8_30 * h4_p4 - e_2 * fs_15_1_14 * r_2 * h2_p2 - e_3 * fs_5_2_5 * h6_p2 - e_3 * fs_65_44_6 * h6_p4 - e_3 * fs_135_44_30 * r_2 * h4_p4 + e_3 * fs_5_1_14 * r_4 * h2_p2 + e_4 * fs_7_22_15 * h8_p2 - e_4 * fs_7_572_66 * h8_p4 + e_4 * fs_2_3_5 * r_2 * h6_p2 + e_4 * fs_13_33_6 * r_2 * h6_p4 + e_4 * fs_135_286_30 * r_4 * h4_p4 - e_4 * fs_20_33_14 * r_6 * h2_p2 + e_5 * fs_189_46189_770 * h10_p2 + e_5 * fs_189_92378_10010 * h10_p4 - e_5 * fs_7_209_15 * r_2 * h8_p2 + e_5 * fs_7_5434_66 * r_2 * h8_p4 - e_5 * fs_2_51_5 * r_4 * h6_p2 - e_5 * fs_13_561_6 * r_4 * h6_p4 - e_5 * fs_3_143_30 * r_6 * h4_p4 + e_5 * fs_10_429_14 * r_8 * h2_p2;

        pc_24[k] = - e_2 * f_135_4 * h4_m3 + e_3 * fs_125_22_3 * h6_m3 + e_3 * f_405_22 * r_2 * h4_m3 - e_4 * fs_35_286_66 * h8_m3 - e_4 * fs_50_33_3 * r_2 * h6_m3 - e_4 * f_405_143 * r_4 * h4_m3 - e_5 * fs_189_46189_3003 * h10_m3 + e_5 * fs_35_2717_66 * r_2 * h8_m3 + e_5 * fs_50_561_3 * r_4 * h6_m3 + e_5 * f_18_143 * r_6 * h4_m3;
    }

    // NOTE: the angular components are formed in 53 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph2_m2, ph4_m4, ph6_m4, ph6_m2, ph8_m4, ph8_m2, ph10_m4, ph10_m2, ab_2, pc_25 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
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

        const auto h2_m2 = ph2_m2[k];
        const auto h4_m4 = ph4_m4[k];
        const auto h6_m4 = ph6_m4[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h8_m4 = ph8_m4[k];
        const auto h8_m2 = ph8_m2[k];
        const auto h10_m4 = ph10_m4[k];
        const auto h10_m2 = ph10_m2[k];

        pc_25[k] = e_1 * fs_105_8_14 * h2_m2 - e_2 * fs_45_8_30 * h4_m4 - e_2 * fs_15_1_14 * r_2 * h2_m2 + e_3 * fs_65_44_6 * h6_m4 - e_3 * fs_5_2_5 * h6_m2 + e_3 * fs_135_44_30 * r_2 * h4_m4 + e_3 * fs_5_1_14 * r_4 * h2_m2 + e_4 * fs_7_572_66 * h8_m4 + e_4 * fs_7_22_15 * h8_m2 - e_4 * fs_13_33_6 * r_2 * h6_m4 + e_4 * fs_2_3_5 * r_2 * h6_m2 - e_4 * fs_135_286_30 * r_4 * h4_m4 - e_4 * fs_20_33_14 * r_6 * h2_m2 - e_5 * fs_189_92378_10010 * h10_m4 + e_5 * fs_189_46189_770 * h10_m2 - e_5 * fs_7_5434_66 * r_2 * h8_m4 - e_5 * fs_7_209_15 * r_2 * h8_m2 + e_5 * fs_13_561_6 * r_4 * h6_m4 - e_5 * fs_2_51_5 * r_4 * h6_m2 + e_5 * fs_3_143_30 * r_6 * h4_m4 + e_5 * fs_10_429_14 * r_8 * h2_m2;
    }

    // NOTE: the angular components are formed in 53 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph2_m1, ph4_m1, ph6_m6, ph6_m5, ph6_m1, ph8_m6, ph8_m5, ph8_m1, ph10_m6, ph10_m5, ph10_m1, ab_2, pc_26, pc_27 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
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

        pc_26[k] = - e_1 * fs_525_16_2 * h2_m1 + e_2 * fs_45_8_15 * h4_m1 + e_2 * fs_75_2_2 * r_2 * h2_m1 - e_3 * fs_5_22_231 * h6_m5 + e_3 * fs_5_22_14 * h6_m1 - e_3 * fs_135_44_15 * r_2 * h4_m1 - e_3 * fs_25_2_2 * r_4 * h2_m1 + e_4 * fs_7_572_6006 * h8_m5 - e_4 * fs_329_572_6 * h8_m1 + e_4 * fs_2_33_231 * r_2 * h6_m5 - e_4 * fs_2_33_14 * r_2 * h6_m1 + e_4 * fs_135_286_15 * r_4 * h4_m1 + e_4 * fs_50_33_2 * r_6 * h2_m1 - e_5 * fs_945_46189_143 * h10_m5 - e_5 * fs_189_46189_330 * h10_m1 - e_5 * fs_7_5434_6006 * r_2 * h8_m5 + e_5 * fs_329_5434_6 * r_2 * h8_m1 - e_5 * fs_2_561_231 * r_4 * h6_m5 + e_5 * fs_2_561_14 * r_4 * h6_m1 - e_5 * fs_3_143_15 * r_6 * h4_m1 - e_5 * fs_25_429_2 * r_8 * h2_m1;

        pc_27[k] = - e_3 * fs_5_11_462 * h6_m6 + e_4 * fs_7_143_858 * h8_m6 + e_4 * fs_4_33_462 * r_2 * h6_m6 - e_5 * fs_189_46189_4290 * h10_m6 - e_5 * fs_14_2717_858 * r_2 * h8_m6 - e_5 * fs_4_561_462 * r_4 * h6_m6;
    }

    // NOTE: the angular components are formed in 53 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph2_m1, ph6_m1, ph8_m7, ph8_m1, ph10_m7, ph10_m1, ab_2, pc_28 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
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

        const auto h2_m1 = ph2_m1[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h8_m7 = ph8_m7[k];
        const auto h8_m1 = ph8_m1[k];
        const auto h10_m7 = ph10_m7[k];
        const auto h10_m1 = ph10_m1[k];

        pc_28[k] = e_1 * fs_735_32_6 * h2_m1 - e_2 * fs_105_4_6 * r_2 * h2_m1 - e_3 * fs_5_4_42 * h6_m1 + e_3 * fs_35_4_6 * r_4 * h2_m1 + e_4 * fs_7_286_1430 * h8_m7 - e_4 * fs_7_11_2 * h8_m1 + e_4 * fs_1_3_42 * r_2 * h6_m1 - e_4 * fs_35_33_6 * r_6 * h2_m1 - e_5 * fs_63_46189_36465 * h10_m7 - e_5 * fs_189_92378_110 * h10_m1 - e_5 * fs_7_2717_1430 * r_2 * h8_m7 + e_5 * fs_14_209_2 * r_2 * h8_m1 - e_5 * fs_1_51_42 * r_4 * h6_m1 + e_5 * fs_35_858_6 * r_8 * h2_m1;
    }

    // NOTE: the angular components are formed in 53 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph2_m2, ph4_m2, ph6_m2, ph8_m8, ph8_m2, ph10_m8, ph10_m2, ab_2, pc_29 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
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

        const auto h2_m2 = ph2_m2[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h8_m8 = ph8_m8[k];
        const auto h8_m2 = ph8_m2[k];
        const auto h10_m8 = ph10_m8[k];
        const auto h10_m2 = ph10_m2[k];

        pc_29[k] = - e_1 * fs_105_16_15 * h2_m2 - e_2 * f_135_4 * h4_m2 + e_2 * fs_15_2_15 * r_2 * h2_m2 - e_3 * fs_25_22_42 * h6_m2 + e_3 * f_405_22 * r_2 * h4_m2 - e_3 * fs_5_2_15 * r_4 * h2_m2 - e_4 * fs_35_286_143 * h8_m8 - e_4 * fs_35_286_14 * h8_m2 + e_4 * fs_10_33_42 * r_2 * h6_m2 - e_4 * f_405_143 * r_4 * h4_m2 + e_4 * fs_10_33_15 * r_6 * h2_m2 - e_5 * fs_189_46189_2431 * h10_m8 - e_5 * fs_63_46189_33 * h10_m2 + e_5 * fs_35_2717_143 * r_2 * h8_m8 + e_5 * fs_35_2717_14 * r_2 * h8_m2 - e_5 * fs_10_561_42 * r_4 * h6_m2 + e_5 * f_18_143 * r_6 * h4_m2 - e_5 * fs_5_429_15 * r_8 * h2_m2;
    }

    // NOTE: the angular components are formed in 53 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph2_0, ph4_0, ph4_p4, ph6_0, ph6_p4, ph8_0, ph8_p4, ph10_0, ph10_p4, ab_2, pc_30 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
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
        const auto h4_p4 = ph4_p4[k];
        const auto h6_0 = ph6_0[k];
        const auto h6_p4 = ph6_p4[k];
        const auto h8_0 = ph8_0[k];
        const auto h8_p4 = ph8_p4[k];
        const auto h10_0 = ph10_0[k];
        const auto h10_p4 = ph10_p4[k];

        pc_30[k] = e_0 * f_945_32 - e_1 * f_315_8 * h2_0 - e_1 * f_1575_16 * r_2 - e_2 * f_45_8 * h4_0 - e_2 * fs_45_8_35 * h4_p4 + e_2 * f_45_1 * r_2 * h2_0 + e_2 * f_315_4 * r_4 + e_3 * f_90_11 * h6_0 + e_3 * fs_30_11_7 * h6_p4 + e_3 * f_135_44 * r_2 * h4_0 + e_3 * fs_135_44_35 * r_2 * h4_p4 - e_3 * f_15_1 * r_4 * h2_0 - e_3 * f_45_2 * r_6 - e_4 * f_238_143 * h8_0 - e_4 * fs_21_143_77 * h8_p4 - e_4 * f_24_11 * r_2 * h6_0 - e_4 * fs_8_11_7 * r_2 * h6_p4 - e_4 * f_135_286 * r_4 * h4_0 - e_4 * fs_135_286_35 * r_4 * h4_p4 + e_4 * f_20_11 * r_6 * h2_0 + e_4 * f_5_2 * r_8 - e_5 * f_7560_46189 * h10_0 + e_5 * fs_252_46189_2145 * h10_p4 + e_5 * f_476_2717 * r_2 * h8_0 + e_5 * fs_42_2717_77 * r_2 * h8_p4 + e_5 * f_24_187 * r_4 * h6_0 + e_5 * fs_8_187_7 * r_4 * h6_p4 + e_5 * f_3_143 * r_6 * h4_0 + e_5 * fs_3_143_35 * r_6 * h4_p4 - e_5 * f_10_143 * r_8 * h2_0 - e_5 * f_1_11 * r_10;
    }

    // NOTE: the angular components are formed in 53 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph2_p1, ph4_p1, ph4_p3, ph6_p1, ph6_p3, ph8_p1, ph8_p3, ph10_p1, ph10_p3, ab_2, pc_31 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
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

        const auto h2_p1 = ph2_p1[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h4_p3 = ph4_p3[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h8_p3 = ph8_p3[k];
        const auto h10_p1 = ph10_p1[k];
        const auto h10_p3 = ph10_p3[k];

        pc_31[k] = - e_1 * fs_105_16_21 * h2_p1 + e_2 * fs_45_16_70 * h4_p1 - e_2 * fs_45_16_10 * h4_p3 + e_2 * fs_15_2_21 * r_2 * h2_p1 - e_3 * fs_40_11_3 * h6_p1 + e_3 * fs_15_22_30 * h6_p3 - e_3 * fs_135_88_70 * r_2 * h4_p1 + e_3 * fs_135_88_10 * r_2 * h4_p3 - e_3 * fs_5_2_21 * r_4 * h2_p1 + e_4 * fs_63_286_7 * h8_p1 - e_4 * fs_21_286_165 * h8_p3 + e_4 * fs_32_33_3 * r_2 * h6_p1 - e_4 * fs_2_11_30 * r_2 * h6_p3 + e_4 * fs_135_572_70 * r_4 * h4_p1 - e_4 * fs_135_572_10 * r_4 * h4_p3 + e_4 * fs_10_33_21 * r_6 * h2_p1 + e_5 * fs_378_46189_385 * h10_p1 + e_5 * fs_63_46189_30030 * h10_p3 - e_5 * fs_63_2717_7 * r_2 * h8_p1 + e_5 * fs_21_2717_165 * r_2 * h8_p3 - e_5 * fs_32_561_3 * r_4 * h6_p1 + e_5 * fs_2_187_30 * r_4 * h6_p3 - e_5 * fs_3_286_70 * r_6 * h4_p1 + e_5 * fs_3_286_10 * r_6 * h4_p3 - e_5 * fs_5_429_21 * r_8 * h2_p1;
    }

    // NOTE: the angular components are formed in 53 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph2_m2, ph4_m2, ph6_m2, ph8_m2, ph10_m2, ab_2, pc_32 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
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

        const auto h2_m2 = ph2_m2[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h8_m2 = ph8_m2[k];
        const auto h10_m2 = ph10_m2[k];

        pc_32[k] = e_1 * fs_105_8_35 * h2_m2 - e_2 * fs_45_8_21 * h4_m2 - e_2 * fs_15_1_35 * r_2 * h2_m2 + e_3 * fs_25_11_2 * h6_m2 + e_3 * fs_135_44_21 * r_2 * h4_m2 + e_3 * fs_5_1_35 * r_4 * h2_m2 + e_4 * fs_35_286_6 * h8_m2 - e_4 * fs_20_33_2 * r_2 * h6_m2 - e_4 * fs_135_286_21 * r_4 * h4_m2 - e_4 * fs_20_33_35 * r_6 * h2_m2 - e_5 * fs_1512_46189_77 * h10_m2 - e_5 * fs_35_2717_6 * r_2 * h8_m2 + e_5 * fs_20_561_2 * r_4 * h6_m2 + e_5 * fs_3_143_21 * r_6 * h4_m2 + e_5 * fs_10_429_35 * r_8 * h2_m2;
    }

    // NOTE: the angular components are formed in 53 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph2_m1, ph4_m3, ph4_m1, ph6_m3, ph6_m1, ph8_m3, ph8_m1, ph10_m3, ph10_m1, ab_2, pc_33 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
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

        const auto h2_m1 = ph2_m1[k];
        const auto h4_m3 = ph4_m3[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h8_m3 = ph8_m3[k];
        const auto h8_m1 = ph8_m1[k];
        const auto h10_m3 = ph10_m3[k];
        const auto h10_m1 = ph10_m1[k];

        pc_33[k] = - e_1 * fs_105_16_21 * h2_m1 + e_2 * fs_45_16_10 * h4_m3 + e_2 * fs_45_16_70 * h4_m1 + e_2 * fs_15_2_21 * r_2 * h2_m1 - e_3 * fs_15_22_30 * h6_m3 - e_3 * fs_40_11_3 * h6_m1 - e_3 * fs_135_88_10 * r_2 * h4_m3 - e_3 * fs_135_88_70 * r_2 * h4_m1 - e_3 * fs_5_2_21 * r_4 * h2_m1 + e_4 * fs_21_286_165 * h8_m3 + e_4 * fs_63_286_7 * h8_m1 + e_4 * fs_2_11_30 * r_2 * h6_m3 + e_4 * fs_32_33_3 * r_2 * h6_m1 + e_4 * fs_135_572_10 * r_4 * h4_m3 + e_4 * fs_135_572_70 * r_4 * h4_m1 + e_4 * fs_10_33_21 * r_6 * h2_m1 - e_5 * fs_63_46189_30030 * h10_m3 + e_5 * fs_378_46189_385 * h10_m1 - e_5 * fs_21_2717_165 * r_2 * h8_m3 - e_5 * fs_63_2717_7 * r_2 * h8_m1 - e_5 * fs_2_187_30 * r_4 * h6_m3 - e_5 * fs_32_561_3 * r_4 * h6_m1 - e_5 * fs_3_286_10 * r_6 * h4_m3 - e_5 * fs_3_286_70 * r_6 * h4_m1 - e_5 * fs_5_429_21 * r_8 * h2_m1;
    }

    // NOTE: the angular components are formed in 53 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph2_m1, ph4_m4, ph4_m1, ph6_m5, ph6_m4, ph6_m1, ph8_m5, ph8_m4, ph8_m1, ph10_m5, ph10_m4, ph10_m1, ab_2, pc_34, pc_35 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
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

        pc_34[k] = e_2 * fs_45_8_35 * h4_m4 - e_3 * fs_30_11_7 * h6_m4 - e_3 * fs_135_44_35 * r_2 * h4_m4 + e_4 * fs_21_143_77 * h8_m4 + e_4 * fs_8_11_7 * r_2 * h6_m4 + e_4 * fs_135_286_35 * r_4 * h4_m4 - e_5 * fs_252_46189_2145 * h10_m4 - e_5 * fs_42_2717_77 * r_2 * h8_m4 - e_5 * fs_8_187_7 * r_4 * h6_m4 - e_5 * fs_3_143_35 * r_6 * h4_m4;

        pc_35[k] = e_1 * fs_525_16_2 * h2_m1 - e_2 * fs_45_8_15 * h4_m1 - e_2 * fs_75_2_2 * r_2 * h2_m1 - e_3 * fs_5_22_231 * h6_m5 - e_3 * fs_5_22_14 * h6_m1 + e_3 * fs_135_44_15 * r_2 * h4_m1 + e_3 * fs_25_2_2 * r_4 * h2_m1 + e_4 * fs_7_572_6006 * h8_m5 + e_4 * fs_329_572_6 * h8_m1 + e_4 * fs_2_33_231 * r_2 * h6_m5 + e_4 * fs_2_33_14 * r_2 * h6_m1 - e_4 * fs_135_286_15 * r_4 * h4_m1 - e_4 * fs_50_33_2 * r_6 * h2_m1 - e_5 * fs_945_46189_143 * h10_m5 + e_5 * fs_189_46189_330 * h10_m1 - e_5 * fs_7_5434_6006 * r_2 * h8_m5 - e_5 * fs_329_5434_6 * r_2 * h8_m1 - e_5 * fs_2_561_231 * r_4 * h6_m5 - e_5 * fs_2_561_14 * r_4 * h6_m1 + e_5 * fs_3_143_15 * r_6 * h4_m1 + e_5 * fs_25_429_2 * r_8 * h2_m1;
    }

    // NOTE: the angular components are formed in 53 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph2_m2, ph4_m2, ph6_m6, ph6_m2, ph8_m6, ph8_m2, ph10_m6, ph10_m2, ab_2, pc_36 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
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

        const auto h2_m2 = ph2_m2[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h6_m6 = ph6_m6[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h8_m6 = ph8_m6[k];
        const auto h8_m2 = ph8_m2[k];
        const auto h10_m6 = ph10_m6[k];
        const auto h10_m2 = ph10_m2[k];

        pc_36[k] = - e_1 * f_315_8 * h2_m2 - e_2 * fs_45_8_15 * h4_m2 + e_2 * f_45_1 * r_2 * h2_m2 + e_3 * fs_15_22_154 * h6_m6 + e_3 * fs_15_22_70 * h6_m2 + e_3 * fs_135_44_15 * r_2 * h4_m2 - e_3 * f_15_1 * r_4 * h2_m2 - e_4 * fs_7_572_286 * h8_m6 + e_4 * fs_49_572_210 * h8_m2 - e_4 * fs_2_11_154 * r_2 * h6_m6 - e_4 * fs_2_11_70 * r_2 * h6_m2 - e_4 * fs_135_286_15 * r_4 * h4_m2 + e_4 * f_20_11 * r_6 * h2_m2 - e_5 * fs_252_46189_1430 * h10_m6 + e_5 * fs_252_46189_55 * h10_m2 + e_5 * fs_7_5434_286 * r_2 * h8_m6 - e_5 * fs_49_5434_210 * r_2 * h8_m2 + e_5 * fs_2_187_154 * r_4 * h6_m6 + e_5 * fs_2_187_70 * r_4 * h6_m2 + e_5 * fs_3_143_15 * r_6 * h4_m2 - e_5 * f_10_143 * r_8 * h2_m2;
    }

    // NOTE: the angular components are formed in 53 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_2, pe_3, pe_4, pe_5, ph4_m3, ph6_m3, ph8_m7, ph8_m3, ph10_m7, ph10_m3, ab_2, pc_37 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h4_m3 = ph4_m3[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h8_m7 = ph8_m7[k];
        const auto h8_m3 = ph8_m3[k];
        const auto h10_m7 = ph10_m7[k];
        const auto h10_m3 = ph10_m3[k];

        pc_37[k] = e_2 * fs_45_8_21 * h4_m3 + e_3 * fs_75_22_7 * h6_m3 - e_3 * fs_135_44_21 * r_2 * h4_m3 - e_4 * fs_35_572_858 * h8_m7 + e_4 * fs_35_572_154 * h8_m3 - e_4 * fs_10_11_7 * r_2 * h6_m3 + e_4 * fs_135_286_21 * r_4 * h4_m3 - e_5 * fs_126_46189_2431 * h10_m7 + e_5 * fs_63_46189_143 * h10_m3 + e_5 * fs_35_5434_858 * r_2 * h8_m7 - e_5 * fs_35_5434_154 * r_2 * h8_m3 + e_5 * fs_10_187_7 * r_4 * h6_m3 - e_5 * fs_3_143_21 * r_6 * h4_m3;
    }

    // NOTE: the angular components are formed in 53 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph2_0, ph2_p2, ph4_0, ph4_p2, ph6_0, ph6_p2, ph8_0, ph8_p2, ph10_0, ph10_p2, ab_2, pc_38 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
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

        pc_38[k] = e_0 * f_945_32 - e_1 * f_945_16 * h2_0 + e_1 * fs_525_16_3 * h2_p2 - e_1 * f_1575_16 * r_2 + e_2 * f_45_2 * h4_0 - e_2 * fs_45_4_5 * h4_p2 + e_2 * f_135_2 * r_2 * h2_0 - e_2 * fs_75_2_3 * r_2 * h2_p2 + e_2 * f_315_4 * r_4 - e_3 * f_30_11 * h6_0 + e_3 * fs_5_11_210 * h6_p2 - e_3 * f_135_11 * r_2 * h4_0 + e_3 * fs_135_22_5 * r_2 * h4_p2 - e_3 * f_45_2 * r_4 * h2_0 + e_3 * fs_25_2_3 * r_4 * h2_p2 - e_3 * f_45_2 * r_6 - e_4 * f_49_143 * h8_0 - e_4 * fs_21_143_70 * h8_p2 + e_4 * f_8_11 * r_2 * h6_0 - e_4 * fs_4_33_210 * r_2 * h6_p2 + e_4 * f_270_143 * r_4 * h4_0 - e_4 * fs_135_143_5 * r_4 * h4_p2 + e_4 * f_30_11 * r_6 * h2_0 - e_4 * fs_50_33_3 * r_6 * h2_p2 + e_4 * f_5_2 * r_8 + e_5 * f_13230_46189 * h10_0 + e_5 * fs_882_46189_165 * h10_p2 + e_5 * f_98_2717 * r_2 * h8_0 + e_5 * fs_42_2717_70 * r_2 * h8_p2 - e_5 * f_8_187 * r_4 * h6_0 + e_5 * fs_4_561_210 * r_4 * h6_p2 - e_5 * f_12_143 * r_6 * h4_0 + e_5 * fs_6_143_5 * r_6 * h4_p2 - e_5 * f_15_143 * r_8 * h2_0 + e_5 * fs_25_429_3 * r_8 * h2_p2 - e_5 * f_1_11 * r_10;
    }

    // NOTE: the angular components are formed in 53 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph2_m2, ph2_m1, ph4_m2, ph4_m1, ph6_m2, ph6_m1, ph8_m2, ph8_m1, ph10_m2, ph10_m1, ab_2, pc_39, pc_40 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
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

        pc_39[k] = - e_1 * fs_105_16_5 * h2_m1 + e_2 * fs_45_8_6 * h4_m1 + e_2 * fs_15_2_5 * r_2 * h2_m1 - e_3 * fs_10_11_35 * h6_m1 - e_3 * fs_135_44_6 * r_2 * h4_m1 - e_3 * fs_5_2_5 * r_4 * h2_m1 + e_4 * fs_49_143_15 * h8_m1 + e_4 * fs_8_33_35 * r_2 * h6_m1 + e_4 * fs_135_286_6 * r_4 * h4_m1 + e_4 * fs_10_33_5 * r_6 * h2_m1 - e_5 * fs_2646_46189_33 * h10_m1 - e_5 * fs_98_2717_15 * r_2 * h8_m1 - e_5 * fs_8_561_35 * r_4 * h6_m1 - e_5 * fs_3_143_6 * r_6 * h4_m1 - e_5 * fs_5_429_5 * r_8 * h2_m1;

        pc_40[k] = - e_1 * fs_525_16_3 * h2_m2 + e_2 * fs_45_4_5 * h4_m2 + e_2 * fs_75_2_3 * r_2 * h2_m2 - e_3 * fs_5_11_210 * h6_m2 - e_3 * fs_135_22_5 * r_2 * h4_m2 - e_3 * fs_25_2_3 * r_4 * h2_m2 + e_4 * fs_21_143_70 * h8_m2 + e_4 * fs_4_33_210 * r_2 * h6_m2 + e_4 * fs_135_143_5 * r_4 * h4_m2 + e_4 * fs_50_33_3 * r_6 * h2_m2 - e_5 * fs_882_46189_165 * h10_m2 - e_5 * fs_42_2717_70 * r_2 * h8_m2 - e_5 * fs_4_561_210 * r_4 * h6_m2 - e_5 * fs_6_143_5 * r_6 * h4_m2 - e_5 * fs_25_429_3 * r_8 * h2_m2;
    }

    // NOTE: the angular components are formed in 53 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph2_m1, ph4_m3, ph4_m1, ph6_m3, ph6_m1, ph8_m3, ph8_m1, ph10_m3, ph10_m1, ab_2, pc_41 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
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

        const auto h2_m1 = ph2_m1[k];
        const auto h4_m3 = ph4_m3[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h8_m3 = ph8_m3[k];
        const auto h8_m1 = ph8_m1[k];
        const auto h10_m3 = ph10_m3[k];
        const auto h10_m1 = ph10_m1[k];

        pc_41[k] = e_1 * fs_105_16_21 * h2_m1 + e_2 * fs_45_16_10 * h4_m3 - e_2 * fs_45_16_70 * h4_m1 - e_2 * fs_15_2_21 * r_2 * h2_m1 - e_3 * fs_15_22_30 * h6_m3 + e_3 * fs_40_11_3 * h6_m1 - e_3 * fs_135_88_10 * r_2 * h4_m3 + e_3 * fs_135_88_70 * r_2 * h4_m1 + e_3 * fs_5_2_21 * r_4 * h2_m1 + e_4 * fs_21_286_165 * h8_m3 - e_4 * fs_63_286_7 * h8_m1 + e_4 * fs_2_11_30 * r_2 * h6_m3 - e_4 * fs_32_33_3 * r_2 * h6_m1 + e_4 * fs_135_572_10 * r_4 * h4_m3 - e_4 * fs_135_572_70 * r_4 * h4_m1 - e_4 * fs_10_33_21 * r_6 * h2_m1 - e_5 * fs_63_46189_30030 * h10_m3 - e_5 * fs_378_46189_385 * h10_m1 - e_5 * fs_21_2717_165 * r_2 * h8_m3 + e_5 * fs_63_2717_7 * r_2 * h8_m1 - e_5 * fs_2_187_30 * r_4 * h6_m3 + e_5 * fs_32_561_3 * r_4 * h6_m1 - e_5 * fs_3_286_10 * r_6 * h4_m3 + e_5 * fs_3_286_70 * r_6 * h4_m1 + e_5 * fs_5_429_21 * r_8 * h2_m1;
    }

    // NOTE: the angular components are formed in 53 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph2_m2, ph4_m4, ph6_m4, ph6_m2, ph8_m4, ph8_m2, ph10_m4, ph10_m2, ab_2, pc_42 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
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

        const auto h2_m2 = ph2_m2[k];
        const auto h4_m4 = ph4_m4[k];
        const auto h6_m4 = ph6_m4[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h8_m4 = ph8_m4[k];
        const auto h8_m2 = ph8_m2[k];
        const auto h10_m4 = ph10_m4[k];
        const auto h10_m2 = ph10_m2[k];

        pc_42[k] = - e_1 * fs_105_8_14 * h2_m2 - e_2 * fs_45_8_30 * h4_m4 + e_2 * fs_15_1_14 * r_2 * h2_m2 + e_3 * fs_65_44_6 * h6_m4 + e_3 * fs_5_2_5 * h6_m2 + e_3 * fs_135_44_30 * r_2 * h4_m4 - e_3 * fs_5_1_14 * r_4 * h2_m2 + e_4 * fs_7_572_66 * h8_m4 - e_4 * fs_7_22_15 * h8_m2 - e_4 * fs_13_33_6 * r_2 * h6_m4 - e_4 * fs_2_3_5 * r_2 * h6_m2 - e_4 * fs_135_286_30 * r_4 * h4_m4 + e_4 * fs_20_33_14 * r_6 * h2_m2 - e_5 * fs_189_92378_10010 * h10_m4 - e_5 * fs_189_46189_770 * h10_m2 - e_5 * fs_7_5434_66 * r_2 * h8_m4 + e_5 * fs_7_209_15 * r_2 * h8_m2 + e_5 * fs_13_561_6 * r_4 * h6_m4 + e_5 * fs_2_51_5 * r_4 * h6_m2 + e_5 * fs_3_143_30 * r_6 * h4_m4 - e_5 * fs_10_429_14 * r_8 * h2_m2;
    }

    // NOTE: the angular components are formed in 53 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_2, pe_3, pe_4, pe_5, ph4_m3, ph6_m5, ph6_m3, ph8_m5, ph8_m3, ph10_m5, ph10_m3, ab_2, pc_43 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

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

        pc_43[k] = e_2 * fs_45_8_30 * h4_m3 + e_3 * fs_45_44_66 * h6_m5 - e_3 * fs_15_44_10 * h6_m3 - e_3 * fs_135_44_30 * r_2 * h4_m3 - e_4 * fs_7_143_429 * h8_m5 - e_4 * fs_28_143_55 * h8_m3 - e_4 * fs_3_11_66 * r_2 * h6_m5 + e_4 * fs_1_11_10 * r_2 * h6_m3 + e_4 * fs_135_286_30 * r_4 * h4_m3 - e_5 * fs_315_92378_2002 * h10_m5 - e_5 * fs_63_92378_10010 * h10_m3 + e_5 * fs_14_2717_429 * r_2 * h8_m5 + e_5 * fs_56_2717_55 * r_2 * h8_m3 + e_5 * fs_3_187_66 * r_4 * h6_m5 - e_5 * fs_1_187_10 * r_4 * h6_m3 - e_5 * fs_3_143_30 * r_6 * h4_m3;
    }

    // NOTE: the angular components are formed in 53 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_2, pe_3, pe_4, pe_5, ph4_m4, ph6_m6, ph6_m4, ph8_m6, ph8_m4, ph10_m6, ph10_m4, ab_2, pc_44 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h4_m4 = ph4_m4[k];
        const auto h6_m6 = ph6_m6[k];
        const auto h6_m4 = ph6_m4[k];
        const auto h8_m6 = ph8_m6[k];
        const auto h8_m4 = ph8_m4[k];
        const auto h10_m6 = ph10_m6[k];
        const auto h10_m4 = ph10_m4[k];

        pc_44[k] = - e_2 * fs_45_8_6 * h4_m4 - e_3 * fs_15_22_55 * h6_m6 - e_3 * fs_75_44_30 * h6_m4 + e_3 * fs_135_44_6 * r_2 * h4_m4 - e_4 * fs_7_286_5005 * h8_m6 - e_4 * fs_35_572_330 * h8_m4 + e_4 * fs_2_11_55 * r_2 * h6_m6 + e_4 * fs_5_11_30 * r_2 * h6_m4 - e_4 * fs_135_286_6 * r_4 * h4_m4 - e_5 * fs_126_46189_1001 * h10_m6 - e_5 * fs_63_92378_2002 * h10_m4 + e_5 * fs_7_2717_5005 * r_2 * h8_m6 + e_5 * fs_35_5434_330 * r_2 * h8_m4 - e_5 * fs_2_187_55 * r_4 * h6_m6 - e_5 * fs_5_187_30 * r_4 * h6_m4 + e_5 * fs_3_143_6 * r_6 * h4_m4;
    }

    // NOTE: the angular components are formed in 53 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph2_0, ph2_p1, ph4_0, ph4_p1, ph6_0, ph6_p1, ph8_0, ph8_p1, ph10_0, ph10_p1, ab_2, pc_45, pc_46 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
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

        pc_45[k] = e_0 * f_945_32 - e_1 * f_525_8 * h2_0 - e_1 * f_1575_16 * r_2 + e_2 * f_135_4 * h4_0 + e_2 * f_75_1 * r_2 * h2_0 + e_2 * f_315_4 * r_4 - e_3 * f_100_11 * h6_0 - e_3 * f_405_22 * r_2 * h4_0 - e_3 * f_25_1 * r_4 * h2_0 - e_3 * f_45_2 * r_6 + e_4 * f_245_143 * h8_0 + e_4 * f_80_33 * r_2 * h6_0 + e_4 * f_405_143 * r_4 * h4_0 + e_4 * f_100_33 * r_6 * h2_0 + e_4 * f_5_2 * r_8 - e_5 * f_15876_46189 * h10_0 - e_5 * f_490_2717 * r_2 * h8_0 - e_5 * f_80_561 * r_4 * h6_0 - e_5 * f_18_143 * r_6 * h4_0 - e_5 * f_50_429 * r_8 * h2_0 - e_5 * f_1_11 * r_10;

        pc_46[k] = - e_1 * fs_105_16_5 * h2_p1 + e_2 * fs_45_8_6 * h4_p1 + e_2 * fs_15_2_5 * r_2 * h2_p1 - e_3 * fs_10_11_35 * h6_p1 - e_3 * fs_135_44_6 * r_2 * h4_p1 - e_3 * fs_5_2_5 * r_4 * h2_p1 + e_4 * fs_49_143_15 * h8_p1 + e_4 * fs_8_33_35 * r_2 * h6_p1 + e_4 * fs_135_286_6 * r_4 * h4_p1 + e_4 * fs_10_33_5 * r_6 * h2_p1 - e_5 * fs_2646_46189_33 * h10_p1 - e_5 * fs_98_2717_15 * r_2 * h8_p1 - e_5 * fs_8_561_35 * r_4 * h6_p1 - e_5 * fs_3_143_6 * r_6 * h4_p1 - e_5 * fs_5_429_5 * r_8 * h2_p1;
    }

    // NOTE: the angular components are formed in 53 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph2_p2, ph4_p2, ph4_p3, ph6_p2, ph6_p3, ph8_p2, ph8_p3, ph10_p2, ph10_p3, ab_2, pc_47, pc_48 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
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

        const auto h2_p2 = ph2_p2[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h4_p3 = ph4_p3[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h8_p2 = ph8_p2[k];
        const auto h8_p3 = ph8_p3[k];
        const auto h10_p2 = ph10_p2[k];
        const auto h10_p3 = ph10_p3[k];

        pc_47[k] = e_1 * fs_105_8_35 * h2_p2 - e_2 * fs_45_8_21 * h4_p2 - e_2 * fs_15_1_35 * r_2 * h2_p2 + e_3 * fs_25_11_2 * h6_p2 + e_3 * fs_135_44_21 * r_2 * h4_p2 + e_3 * fs_5_1_35 * r_4 * h2_p2 + e_4 * fs_35_286_6 * h8_p2 - e_4 * fs_20_33_2 * r_2 * h6_p2 - e_4 * fs_135_286_21 * r_4 * h4_p2 - e_4 * fs_20_33_35 * r_6 * h2_p2 - e_5 * fs_1512_46189_77 * h10_p2 - e_5 * fs_35_2717_6 * r_2 * h8_p2 + e_5 * fs_20_561_2 * r_4 * h6_p2 + e_5 * fs_3_143_21 * r_6 * h4_p2 + e_5 * fs_10_429_35 * r_8 * h2_p2;

        pc_48[k] = - e_2 * f_135_4 * h4_p3 + e_3 * fs_125_22_3 * h6_p3 + e_3 * f_405_22 * r_2 * h4_p3 - e_4 * fs_35_286_66 * h8_p3 - e_4 * fs_50_33_3 * r_2 * h6_p3 - e_4 * f_405_143 * r_4 * h4_p3 - e_5 * fs_189_46189_3003 * h10_p3 + e_5 * fs_35_2717_66 * r_2 * h8_p3 + e_5 * fs_50_561_3 * r_4 * h6_p3 + e_5 * f_18_143 * r_6 * h4_p3;
    }

    // NOTE: the angular components are formed in 53 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_2, pe_3, pe_4, pe_5, ph4_p4, ph6_p4, ph6_p5, ph8_p4, ph8_p5, ph10_p4, ph10_p5, ab_2, pc_49, pc_50 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h4_p4 = ph4_p4[k];
        const auto h6_p4 = ph6_p4[k];
        const auto h6_p5 = ph6_p5[k];
        const auto h8_p4 = ph8_p4[k];
        const auto h8_p5 = ph8_p5[k];
        const auto h10_p4 = ph10_p4[k];
        const auto h10_p5 = ph10_p5[k];

        pc_49[k] = e_2 * f_135_4 * h4_p4 + e_3 * fs_30_11_5 * h6_p4 - e_3 * f_405_22 * r_2 * h4_p4 - e_4 * fs_7_26_55 * h8_p4 - e_4 * fs_8_11_5 * r_2 * h6_p4 + e_4 * f_405_143 * r_4 * h4_p4 - e_5 * fs_126_46189_3003 * h10_p4 + e_5 * fs_7_247_55 * r_2 * h8_p4 + e_5 * fs_8_187_5 * r_4 * h6_p4 - e_5 * f_18_143 * r_6 * h4_p4;

        pc_50[k] = - e_3 * fs_75_22_11 * h6_p5 - e_4 * fs_35_286_286 * h8_p5 + e_4 * fs_10_11_11 * r_2 * h6_p5 - e_5 * fs_63_46189_3003 * h10_p5 + e_5 * fs_35_2717_286 * r_2 * h8_p5 - e_5 * fs_10_187_11 * r_4 * h6_p5;
    }

    // NOTE: the angular components are formed in 53 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph2_0, ph2_p2, ph4_0, ph4_p2, ph6_0, ph6_p2, ph8_0, ph8_p2, ph10_0, ph10_p2, ab_2, pc_51 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
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

        pc_51[k] = e_0 * f_945_32 - e_1 * f_945_16 * h2_0 - e_1 * fs_525_16_3 * h2_p2 - e_1 * f_1575_16 * r_2 + e_2 * f_45_2 * h4_0 + e_2 * fs_45_4_5 * h4_p2 + e_2 * f_135_2 * r_2 * h2_0 + e_2 * fs_75_2_3 * r_2 * h2_p2 + e_2 * f_315_4 * r_4 - e_3 * f_30_11 * h6_0 - e_3 * fs_5_11_210 * h6_p2 - e_3 * f_135_11 * r_2 * h4_0 - e_3 * fs_135_22_5 * r_2 * h4_p2 - e_3 * f_45_2 * r_4 * h2_0 - e_3 * fs_25_2_3 * r_4 * h2_p2 - e_3 * f_45_2 * r_6 - e_4 * f_49_143 * h8_0 + e_4 * fs_21_143_70 * h8_p2 + e_4 * f_8_11 * r_2 * h6_0 + e_4 * fs_4_33_210 * r_2 * h6_p2 + e_4 * f_270_143 * r_4 * h4_0 + e_4 * fs_135_143_5 * r_4 * h4_p2 + e_4 * f_30_11 * r_6 * h2_0 + e_4 * fs_50_33_3 * r_6 * h2_p2 + e_4 * f_5_2 * r_8 + e_5 * f_13230_46189 * h10_0 - e_5 * fs_882_46189_165 * h10_p2 + e_5 * f_98_2717 * r_2 * h8_0 - e_5 * fs_42_2717_70 * r_2 * h8_p2 - e_5 * f_8_187 * r_4 * h6_0 - e_5 * fs_4_561_210 * r_4 * h6_p2 - e_5 * f_12_143 * r_6 * h4_0 - e_5 * fs_6_143_5 * r_6 * h4_p2 - e_5 * f_15_143 * r_8 * h2_0 - e_5 * fs_25_429_3 * r_8 * h2_p2 - e_5 * f_1_11 * r_10;
    }

    // NOTE: the angular components are formed in 53 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph2_p1, ph4_p1, ph4_p3, ph6_p1, ph6_p3, ph8_p1, ph8_p3, ph10_p1, ph10_p3, ab_2, pc_52 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
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

        const auto h2_p1 = ph2_p1[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h4_p3 = ph4_p3[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h8_p3 = ph8_p3[k];
        const auto h10_p1 = ph10_p1[k];
        const auto h10_p3 = ph10_p3[k];

        pc_52[k] = - e_1 * fs_105_16_21 * h2_p1 + e_2 * fs_45_16_70 * h4_p1 + e_2 * fs_45_16_10 * h4_p3 + e_2 * fs_15_2_21 * r_2 * h2_p1 - e_3 * fs_40_11_3 * h6_p1 - e_3 * fs_15_22_30 * h6_p3 - e_3 * fs_135_88_70 * r_2 * h4_p1 - e_3 * fs_135_88_10 * r_2 * h4_p3 - e_3 * fs_5_2_21 * r_4 * h2_p1 + e_4 * fs_63_286_7 * h8_p1 + e_4 * fs_21_286_165 * h8_p3 + e_4 * fs_32_33_3 * r_2 * h6_p1 + e_4 * fs_2_11_30 * r_2 * h6_p3 + e_4 * fs_135_572_70 * r_4 * h4_p1 + e_4 * fs_135_572_10 * r_4 * h4_p3 + e_4 * fs_10_33_21 * r_6 * h2_p1 + e_5 * fs_378_46189_385 * h10_p1 - e_5 * fs_63_46189_30030 * h10_p3 - e_5 * fs_63_2717_7 * r_2 * h8_p1 - e_5 * fs_21_2717_165 * r_2 * h8_p3 - e_5 * fs_32_561_3 * r_4 * h6_p1 - e_5 * fs_2_187_30 * r_4 * h6_p3 - e_5 * fs_3_286_70 * r_6 * h4_p1 - e_5 * fs_3_286_10 * r_6 * h4_p3 - e_5 * fs_5_429_21 * r_8 * h2_p1;
    }

    // NOTE: the angular components are formed in 53 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph2_p2, ph4_p4, ph6_p2, ph6_p4, ph8_p2, ph8_p4, ph10_p2, ph10_p4, ab_2, pc_53 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
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

        const auto h2_p2 = ph2_p2[k];
        const auto h4_p4 = ph4_p4[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h6_p4 = ph6_p4[k];
        const auto h8_p2 = ph8_p2[k];
        const auto h8_p4 = ph8_p4[k];
        const auto h10_p2 = ph10_p2[k];
        const auto h10_p4 = ph10_p4[k];

        pc_53[k] = e_1 * fs_105_8_14 * h2_p2 - e_2 * fs_45_8_30 * h4_p4 - e_2 * fs_15_1_14 * r_2 * h2_p2 - e_3 * fs_5_2_5 * h6_p2 + e_3 * fs_65_44_6 * h6_p4 + e_3 * fs_135_44_30 * r_2 * h4_p4 + e_3 * fs_5_1_14 * r_4 * h2_p2 + e_4 * fs_7_22_15 * h8_p2 + e_4 * fs_7_572_66 * h8_p4 + e_4 * fs_2_3_5 * r_2 * h6_p2 - e_4 * fs_13_33_6 * r_2 * h6_p4 - e_4 * fs_135_286_30 * r_4 * h4_p4 - e_4 * fs_20_33_14 * r_6 * h2_p2 + e_5 * fs_189_46189_770 * h10_p2 - e_5 * fs_189_92378_10010 * h10_p4 - e_5 * fs_7_209_15 * r_2 * h8_p2 - e_5 * fs_7_5434_66 * r_2 * h8_p4 - e_5 * fs_2_51_5 * r_4 * h6_p2 + e_5 * fs_13_561_6 * r_4 * h6_p4 + e_5 * fs_3_143_30 * r_6 * h4_p4 + e_5 * fs_10_429_14 * r_8 * h2_p2;
    }

    // NOTE: the angular components are formed in 53 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_2, pe_3, pe_4, pe_5, ph4_p3, ph6_p3, ph6_p5, ph8_p3, ph8_p5, ph10_p3, ph10_p5, ab_2, pc_54 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

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

        pc_54[k] = - e_2 * fs_45_8_30 * h4_p3 + e_3 * fs_15_44_10 * h6_p3 + e_3 * fs_45_44_66 * h6_p5 + e_3 * fs_135_44_30 * r_2 * h4_p3 + e_4 * fs_28_143_55 * h8_p3 - e_4 * fs_7_143_429 * h8_p5 - e_4 * fs_1_11_10 * r_2 * h6_p3 - e_4 * fs_3_11_66 * r_2 * h6_p5 - e_4 * fs_135_286_30 * r_4 * h4_p3 + e_5 * fs_63_92378_10010 * h10_p3 - e_5 * fs_315_92378_2002 * h10_p5 - e_5 * fs_56_2717_55 * r_2 * h8_p3 + e_5 * fs_14_2717_429 * r_2 * h8_p5 + e_5 * fs_1_187_10 * r_4 * h6_p3 + e_5 * fs_3_187_66 * r_4 * h6_p5 + e_5 * fs_3_143_30 * r_6 * h4_p3;
    }

    // NOTE: the angular components are formed in 53 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_2, pe_3, pe_4, pe_5, ph4_p4, ph6_p4, ph6_p6, ph8_p4, ph8_p6, ph10_p4, ph10_p6, ab_2, pc_55 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h4_p4 = ph4_p4[k];
        const auto h6_p4 = ph6_p4[k];
        const auto h6_p6 = ph6_p6[k];
        const auto h8_p4 = ph8_p4[k];
        const auto h8_p6 = ph8_p6[k];
        const auto h10_p4 = ph10_p4[k];
        const auto h10_p6 = ph10_p6[k];

        pc_55[k] = e_2 * fs_45_8_6 * h4_p4 + e_3 * fs_75_44_30 * h6_p4 - e_3 * fs_15_22_55 * h6_p6 - e_3 * fs_135_44_6 * r_2 * h4_p4 + e_4 * fs_35_572_330 * h8_p4 - e_4 * fs_7_286_5005 * h8_p6 - e_4 * fs_5_11_30 * r_2 * h6_p4 + e_4 * fs_2_11_55 * r_2 * h6_p6 + e_4 * fs_135_286_6 * r_4 * h4_p4 + e_5 * fs_63_92378_2002 * h10_p4 - e_5 * fs_126_46189_1001 * h10_p6 - e_5 * fs_35_5434_330 * r_2 * h8_p4 + e_5 * fs_7_2717_5005 * r_2 * h8_p6 + e_5 * fs_5_187_30 * r_4 * h6_p4 - e_5 * fs_2_187_55 * r_4 * h6_p6 - e_5 * fs_3_143_6 * r_6 * h4_p4;
    }

    // NOTE: the angular components are formed in 53 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph2_0, ph4_0, ph4_p4, ph6_0, ph6_p4, ph8_0, ph8_p4, ph10_0, ph10_p4, ab_2, pc_56 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
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
        const auto h4_p4 = ph4_p4[k];
        const auto h6_0 = ph6_0[k];
        const auto h6_p4 = ph6_p4[k];
        const auto h8_0 = ph8_0[k];
        const auto h8_p4 = ph8_p4[k];
        const auto h10_0 = ph10_0[k];
        const auto h10_p4 = ph10_p4[k];

        pc_56[k] = e_0 * f_945_32 - e_1 * f_315_8 * h2_0 - e_1 * f_1575_16 * r_2 - e_2 * f_45_8 * h4_0 + e_2 * fs_45_8_35 * h4_p4 + e_2 * f_45_1 * r_2 * h2_0 + e_2 * f_315_4 * r_4 + e_3 * f_90_11 * h6_0 - e_3 * fs_30_11_7 * h6_p4 + e_3 * f_135_44 * r_2 * h4_0 - e_3 * fs_135_44_35 * r_2 * h4_p4 - e_3 * f_15_1 * r_4 * h2_0 - e_3 * f_45_2 * r_6 - e_4 * f_238_143 * h8_0 + e_4 * fs_21_143_77 * h8_p4 - e_4 * f_24_11 * r_2 * h6_0 + e_4 * fs_8_11_7 * r_2 * h6_p4 - e_4 * f_135_286 * r_4 * h4_0 + e_4 * fs_135_286_35 * r_4 * h4_p4 + e_4 * f_20_11 * r_6 * h2_0 + e_4 * f_5_2 * r_8 - e_5 * f_7560_46189 * h10_0 - e_5 * fs_252_46189_2145 * h10_p4 + e_5 * f_476_2717 * r_2 * h8_0 - e_5 * fs_42_2717_77 * r_2 * h8_p4 + e_5 * f_24_187 * r_4 * h6_0 - e_5 * fs_8_187_7 * r_4 * h6_p4 + e_5 * f_3_143 * r_6 * h4_0 - e_5 * fs_3_143_35 * r_6 * h4_p4 - e_5 * f_10_143 * r_8 * h2_0 - e_5 * f_1_11 * r_10;
    }

    // NOTE: the angular components are formed in 53 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph2_p1, ph4_p1, ph6_p1, ph6_p5, ph8_p1, ph8_p5, ph10_p1, ph10_p5, ab_2, pc_57 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
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

        const auto h2_p1 = ph2_p1[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h6_p5 = ph6_p5[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h8_p5 = ph8_p5[k];
        const auto h10_p1 = ph10_p1[k];
        const auto h10_p5 = ph10_p5[k];

        pc_57[k] = - e_1 * fs_525_16_2 * h2_p1 + e_2 * fs_45_8_15 * h4_p1 + e_2 * fs_75_2_2 * r_2 * h2_p1 + e_3 * fs_5_22_14 * h6_p1 - e_3 * fs_5_22_231 * h6_p5 - e_3 * fs_135_44_15 * r_2 * h4_p1 - e_3 * fs_25_2_2 * r_4 * h2_p1 - e_4 * fs_329_572_6 * h8_p1 + e_4 * fs_7_572_6006 * h8_p5 - e_4 * fs_2_33_14 * r_2 * h6_p1 + e_4 * fs_2_33_231 * r_2 * h6_p5 + e_4 * fs_135_286_15 * r_4 * h4_p1 + e_4 * fs_50_33_2 * r_6 * h2_p1 - e_5 * fs_189_46189_330 * h10_p1 - e_5 * fs_945_46189_143 * h10_p5 + e_5 * fs_329_5434_6 * r_2 * h8_p1 - e_5 * fs_7_5434_6006 * r_2 * h8_p5 + e_5 * fs_2_561_14 * r_4 * h6_p1 - e_5 * fs_2_561_231 * r_4 * h6_p5 - e_5 * fs_3_143_15 * r_6 * h4_p1 - e_5 * fs_25_429_2 * r_8 * h2_p1;
    }

    // NOTE: the angular components are formed in 53 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph2_p2, ph4_p2, ph6_p2, ph6_p6, ph8_p2, ph8_p6, ph10_p2, ph10_p6, ab_2, pc_58 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
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

        const auto h2_p2 = ph2_p2[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h6_p6 = ph6_p6[k];
        const auto h8_p2 = ph8_p2[k];
        const auto h8_p6 = ph8_p6[k];
        const auto h10_p2 = ph10_p2[k];
        const auto h10_p6 = ph10_p6[k];

        pc_58[k] = e_1 * f_315_8 * h2_p2 + e_2 * fs_45_8_15 * h4_p2 - e_2 * f_45_1 * r_2 * h2_p2 - e_3 * fs_15_22_70 * h6_p2 + e_3 * fs_15_22_154 * h6_p6 - e_3 * fs_135_44_15 * r_2 * h4_p2 + e_3 * f_15_1 * r_4 * h2_p2 - e_4 * fs_49_572_210 * h8_p2 - e_4 * fs_7_572_286 * h8_p6 + e_4 * fs_2_11_70 * r_2 * h6_p2 - e_4 * fs_2_11_154 * r_2 * h6_p6 + e_4 * fs_135_286_15 * r_4 * h4_p2 - e_4 * f_20_11 * r_6 * h2_p2 - e_5 * fs_252_46189_55 * h10_p2 - e_5 * fs_252_46189_1430 * h10_p6 + e_5 * fs_49_5434_210 * r_2 * h8_p2 + e_5 * fs_7_5434_286 * r_2 * h8_p6 - e_5 * fs_2_187_70 * r_4 * h6_p2 + e_5 * fs_2_187_154 * r_4 * h6_p6 - e_5 * fs_3_143_15 * r_6 * h4_p2 + e_5 * f_10_143 * r_8 * h2_p2;
    }

    // NOTE: the angular components are formed in 53 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_2, pe_3, pe_4, pe_5, ph4_p3, ph6_p3, ph8_p3, ph8_p7, ph10_p3, ph10_p7, ab_2, pc_59 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h4_p3 = ph4_p3[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h8_p3 = ph8_p3[k];
        const auto h8_p7 = ph8_p7[k];
        const auto h10_p3 = ph10_p3[k];
        const auto h10_p7 = ph10_p7[k];

        pc_59[k] = - e_2 * fs_45_8_21 * h4_p3 - e_3 * fs_75_22_7 * h6_p3 + e_3 * fs_135_44_21 * r_2 * h4_p3 - e_4 * fs_35_572_154 * h8_p3 - e_4 * fs_35_572_858 * h8_p7 + e_4 * fs_10_11_7 * r_2 * h6_p3 - e_4 * fs_135_286_21 * r_4 * h4_p3 - e_5 * fs_63_46189_143 * h10_p3 - e_5 * fs_126_46189_2431 * h10_p7 + e_5 * fs_35_5434_154 * r_2 * h8_p3 + e_5 * fs_35_5434_858 * r_2 * h8_p7 - e_5 * fs_10_187_7 * r_4 * h6_p3 + e_5 * fs_3_143_21 * r_6 * h4_p3;
    }

    // NOTE: the angular components are formed in 53 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph2_0, ph4_0, ph6_0, ph6_p6, ph8_0, ph8_p6, ph10_0, ph10_p6, ab_2, pc_60 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
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

        pc_60[k] = e_0 * f_945_32 - e_1 * f_105_16 * h2_0 - e_1 * f_1575_16 * r_2 - e_2 * f_135_4 * h4_0 + e_2 * f_15_2 * r_2 * h2_0 + e_2 * f_315_4 * r_4 + e_3 * f_145_22 * h6_0 - e_3 * fs_5_11_462 * h6_p6 + e_3 * f_405_22 * r_2 * h4_0 - e_3 * f_5_2 * r_4 * h2_0 - e_3 * f_45_2 * r_6 + e_4 * f_511_286 * h8_0 + e_4 * fs_7_143_858 * h8_p6 - e_4 * f_58_33 * r_2 * h6_0 + e_4 * fs_4_33_462 * r_2 * h6_p6 - e_4 * f_405_143 * r_4 * h4_0 + e_4 * f_10_33 * r_6 * h2_0 + e_4 * f_5_2 * r_8 + e_5 * f_2835_46189 * h10_0 - e_5 * fs_189_46189_4290 * h10_p6 - e_5 * f_511_2717 * r_2 * h8_0 - e_5 * fs_14_2717_858 * r_2 * h8_p6 + e_5 * f_58_561 * r_4 * h6_0 - e_5 * fs_4_561_462 * r_4 * h6_p6 + e_5 * f_18_143 * r_6 * h4_0 - e_5 * f_5_429 * r_8 * h2_0 - e_5 * f_1_11 * r_10;
    }

    // NOTE: the angular components are formed in 53 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph2_p1, ph6_p1, ph8_p1, ph8_p7, ph10_p1, ph10_p7, ab_2, pc_61 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
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

        const auto h2_p1 = ph2_p1[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h8_p7 = ph8_p7[k];
        const auto h10_p1 = ph10_p1[k];
        const auto h10_p7 = ph10_p7[k];

        pc_61[k] = - e_1 * fs_735_32_6 * h2_p1 + e_2 * fs_105_4_6 * r_2 * h2_p1 + e_3 * fs_5_4_42 * h6_p1 - e_3 * fs_35_4_6 * r_4 * h2_p1 + e_4 * fs_7_11_2 * h8_p1 + e_4 * fs_7_286_1430 * h8_p7 - e_4 * fs_1_3_42 * r_2 * h6_p1 + e_4 * fs_35_33_6 * r_6 * h2_p1 + e_5 * fs_189_92378_110 * h10_p1 - e_5 * fs_63_46189_36465 * h10_p7 - e_5 * fs_14_209_2 * r_2 * h8_p1 - e_5 * fs_7_2717_1430 * r_2 * h8_p7 + e_5 * fs_1_51_42 * r_4 * h6_p1 - e_5 * fs_35_858_6 * r_8 * h2_p1;
    }

    // NOTE: the angular components are formed in 53 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph2_p2, ph4_p2, ph6_p2, ph8_p2, ph8_p8, ph10_p2, ph10_p8, ab_2, pc_62 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
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

        const auto h2_p2 = ph2_p2[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h8_p2 = ph8_p2[k];
        const auto h8_p8 = ph8_p8[k];
        const auto h10_p2 = ph10_p2[k];
        const auto h10_p8 = ph10_p8[k];

        pc_62[k] = e_1 * fs_105_16_15 * h2_p2 + e_2 * f_135_4 * h4_p2 - e_2 * fs_15_2_15 * r_2 * h2_p2 + e_3 * fs_25_22_42 * h6_p2 - e_3 * f_405_22 * r_2 * h4_p2 + e_3 * fs_5_2_15 * r_4 * h2_p2 + e_4 * fs_35_286_14 * h8_p2 - e_4 * fs_35_286_143 * h8_p8 - e_4 * fs_10_33_42 * r_2 * h6_p2 + e_4 * f_405_143 * r_4 * h4_p2 - e_4 * fs_10_33_15 * r_6 * h2_p2 + e_5 * fs_63_46189_33 * h10_p2 - e_5 * fs_189_46189_2431 * h10_p8 - e_5 * fs_35_2717_14 * r_2 * h8_p2 + e_5 * fs_35_2717_143 * r_2 * h8_p8 + e_5 * fs_10_561_42 * r_4 * h6_p2 - e_5 * f_18_143 * r_6 * h4_p2 + e_5 * fs_5_429_15 * r_8 * h2_p2;
    }

    // NOTE: the angular components are formed in 53 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph2_0, ph4_0, ph6_0, ph8_0, ph8_p8, ph10_0, ph10_p8, ab_2, pc_63 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
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

        pc_63[k] = e_0 * f_945_32 + e_1 * f_315_8 * h2_0 - e_1 * f_1575_16 * r_2 - e_2 * f_135_4 * h4_0 - e_2 * f_45_1 * r_2 * h2_0 + e_2 * f_315_4 * r_4 - e_3 * f_120_11 * h6_0 + e_3 * f_405_22 * r_2 * h4_0 + e_3 * f_15_1 * r_4 * h2_0 - e_3 * f_45_2 * r_6 - e_4 * f_217_286 * h8_0 + e_4 * fs_21_286_715 * h8_p8 + e_4 * f_32_11 * r_2 * h6_0 - e_4 * f_405_143 * r_4 * h4_0 - e_4 * f_20_11 * r_6 * h2_0 + e_4 * f_5_2 * r_8 - e_5 * f_630_46189 * h10_0 - e_5 * fs_126_46189_12155 * h10_p8 + e_5 * f_217_2717 * r_2 * h8_0 - e_5 * fs_21_2717_715 * r_2 * h8_p8 - e_5 * f_32_187 * r_4 * h6_0 + e_5 * f_18_143 * r_6 * h4_0 + e_5 * f_10_143 * r_8 * h2_0 - e_5 * f_1_11 * r_10;
    }

    // NOTE: the angular components are formed in 53 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph2_p1, ph4_p1, ph6_p1, ph8_p1, ph10_p1, ph10_p9, ab_2, pc_64 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
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

        const auto h2_p1 = ph2_p1[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h10_p1 = ph10_p1[k];
        const auto h10_p9 = ph10_p9[k];

        pc_64[k] = - e_1 * fs_315_32_30 * h2_p1 - e_2 * f_135_4 * h4_p1 + e_2 * fs_45_4_30 * r_2 * h2_p1 - e_3 * fs_15_44_210 * h6_p1 + e_3 * f_405_22 * r_2 * h4_p1 - e_3 * fs_15_4_30 * r_4 * h2_p1 - e_4 * fs_21_286_10 * h8_p1 + e_4 * fs_1_11_210 * r_2 * h6_p1 - e_4 * f_405_143 * r_4 * h4_p1 + e_4 * fs_5_11_30 * r_6 * h2_p1 - e_5 * fs_63_92378_22 * h10_p1 - e_5 * fs_63_46189_46189 * h10_p9 + e_5 * fs_21_2717_10 * r_2 * h8_p1 - e_5 * fs_1_187_210 * r_4 * h6_p1 + e_5 * f_18_143 * r_6 * h4_p1 - e_5 * fs_5_286_30 * r_8 * h2_p1;
    }

    // NOTE: the angular components are formed in 53 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph2_0, ph4_0, ph6_0, ph8_0, ph10_0, ph10_p10, ab_2, pc_65 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
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
        const auto h10_0 = ph10_0[k];
        const auto h10_p10 = ph10_p10[k];

        pc_65[k] = e_0 * f_945_32 + e_1 * f_1575_16 * h2_0 - e_1 * f_1575_16 * r_2 + e_2 * f_135_4 * h4_0 - e_2 * f_225_2 * r_2 * h2_0 + e_2 * f_315_4 * r_4 + e_3 * f_75_22 * h6_0 - e_3 * f_405_22 * r_2 * h4_0 + e_3 * f_75_2 * r_4 * h2_0 - e_3 * f_45_2 * r_6 + e_4 * f_35_286 * h8_0 - e_4 * f_10_11 * r_2 * h6_0 + e_4 * f_405_143 * r_4 * h4_0 - e_4 * f_50_11 * r_6 * h2_0 + e_4 * f_5_2 * r_8 + e_5 * f_63_46189 * h10_0 - e_5 * fs_63_46189_92378 * h10_p10 - e_5 * f_35_2717 * r_2 * h8_0 + e_5 * f_10_187 * r_4 * h6_0 - e_5 * f_18_143 * r_6 * h4_0 + e_5 * f_25_143 * r_8 * h2_0 - e_5 * f_1_11 * r_10;
    }

    // NOTE: the values of a combination of angular components are stored as one
    // row of nvalues columns, with the component on bra side running slowest. The
    // rows which the symmetry relates to an already formed one are copied from it.

    const size_t sources[121] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 1, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 2, 13, 24, 25, 26, 27, 28, 29, 30, 31, 32, 3, 14, 25, 36, 37, 38, 39, 40, 41, 42, 43, 4, 15, 26, 37, 48, 49, 50, 51, 52, 53, 54, 5, 16, 27, 38, 49, 60, 61, 62, 63, 64, 65, 6, 17, 28, 39, 50, 61, 72, 73, 74, 75, 76, 7, 18, 29, 40, 51, 62, 73, 84, 85, 86, 87, 8, 19, 30, 41, 52, 63, 74, 85, 96, 97, 98, 9, 20, 31, 42, 53, 64, 75, 86, 97, 108, 109, 10, 21, 32, 43, 54, 65, 76, 87, 98, 109, 120};

    for (size_t n = 0; n < 121; n++)
    {
        if (sources[n] != n) std::copy(values + sources[n] * nvalues, values + (sources[n] + 1) * nvalues, values + n * nvalues);
    }
}

}  // namespace simdt2ceri
