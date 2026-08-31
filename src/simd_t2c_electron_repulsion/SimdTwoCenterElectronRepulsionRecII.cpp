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



#include "SimdTwoCenterElectronRepulsionRecII.hpp"

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
compute_ii_electron_repulsion(double                         *values,
                              const size_t                    nvalues,
                              const CBasisFunction           &bra,
                              const CBasisFunction           &ket,
                              const std::vector<CSimdMatrix> &harmonics,
                              const CSimdMatrix              &coordinates) -> void
{
    if ((bra.get_angular_momentum() != 6) || (ket.get_angular_momentum() != 6))
    {
        errors::assertMsgCritical(
            false, std::string("SimdTwoCenterElectronRepulsionFunc.compute_ii_electron_repulsion: Basis functions must be of angular momenta six and six"));
    }

    if (harmonics.size() < 12)
    {
        errors::assertMsgCritical(
            false, std::string("SimdTwoCenterElectronRepulsionFunc.compute_ii_electron_repulsion: Harmonics must reach angular momentum 12"));
    }

    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("SimdTwoCenterElectronRepulsionFunc.compute_ii_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    auto boys = CSimdVariableMatrix(std::vector<size_t>(nprims, nvalues), 14);

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
    // call, which fills the orders 0 to 12 of every row. The terms read the
    // orders 6 to 12 alone, and the orders below them are formed on the
    // way to them by the recursion the Boys function is evaluated with.

    simdfunc::compute_boys_function(boys);

    // NOTE: the buffer holds the contracted prefactor of each exponent factor of
    // the terms, as the integrals of the angular components are formed straight
    // into the values and are not written a second time. Every exponent factor is
    // used with one order of the Boys function alone, so the order follows from
    // the factor and one accumulator per factor suffices.

    auto buffer = CSimdMatrix(7, nvalues);

    auto *pe_0 = buffer.data(0);
    auto *pe_1 = buffer.data(1);
    auto *pe_2 = buffer.data(2);
    auto *pe_3 = buffer.data(3);
    auto *pe_4 = buffer.data(4);
    auto *pe_5 = buffer.data(5);
    auto *pe_6 = buffer.data(6);

    std::fill(pe_0, pe_0 + nvalues, 0.0);
    std::fill(pe_1, pe_1 + nvalues, 0.0);
    std::fill(pe_2, pe_2 + nvalues, 0.0);
    std::fill(pe_3, pe_3 + nvalues, 0.0);
    std::fill(pe_4, pe_4 + nvalues, 0.0);
    std::fill(pe_5, pe_5 + nvalues, 0.0);
    std::fill(pe_6, pe_6 + nvalues, 0.0);

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

            const auto ff_0 = fbase / (fexp * fexp * fexp * fexp * fexp * fexp);

            const auto ff_1 = fbase * fmu / (fexp * fexp * fexp * fexp * fexp * fexp);

            const auto ff_2 = fbase * fmu * fmu / (fexp * fexp * fexp * fexp * fexp * fexp);

            const auto ff_3 = fbase * fmu * fmu * fmu / (fexp * fexp * fexp * fexp * fexp * fexp);

            const auto ff_4 = fbase * fmu * fmu * fmu * fmu / (fexp * fexp * fexp * fexp * fexp * fexp);

            const auto ff_5 = fbase * fmu * fmu * fmu * fmu * fmu / (fexp * fexp * fexp * fexp * fexp * fexp);

            const auto ff_6 = fbase * fmu * fmu * fmu * fmu * fmu * fmu / (fexp * fexp * fexp * fexp * fexp * fexp);

            const auto *bv_0 = boys.data(7, i * nprim_b + j);

            const auto *bv_1 = boys.data(8, i * nprim_b + j);

            const auto *bv_2 = boys.data(9, i * nprim_b + j);

            const auto *bv_3 = boys.data(10, i * nprim_b + j);

            const auto *bv_4 = boys.data(11, i * nprim_b + j);

            const auto *bv_5 = boys.data(12, i * nprim_b + j);

            const auto *bv_6 = boys.data(13, i * nprim_b + j);

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, bv_0, bv_1, bv_2, bv_3, bv_4, bv_5, bv_6 : simd::cache_line_size())
            for (size_t k = 0; k < nvalues; k++)
            {
                pe_0[k] += ff_0 * bv_0[k];
                pe_1[k] += ff_1 * bv_1[k];
                pe_2[k] += ff_2 * bv_2[k];
                pe_3[k] += ff_3 * bv_3[k];
                pe_4[k] += ff_4 * bv_4[k];
                pe_5[k] += ff_5 * bv_5[k];
                pe_6[k] += ff_6 * bv_6[k];
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
    auto *pc_13 = values + 14 * nvalues;
    auto *pc_14 = values + 15 * nvalues;
    auto *pc_15 = values + 16 * nvalues;
    auto *pc_16 = values + 17 * nvalues;
    auto *pc_17 = values + 18 * nvalues;
    auto *pc_18 = values + 19 * nvalues;
    auto *pc_19 = values + 20 * nvalues;
    auto *pc_20 = values + 21 * nvalues;
    auto *pc_21 = values + 22 * nvalues;
    auto *pc_22 = values + 23 * nvalues;
    auto *pc_23 = values + 24 * nvalues;
    auto *pc_24 = values + 25 * nvalues;
    auto *pc_25 = values + 28 * nvalues;
    auto *pc_26 = values + 29 * nvalues;
    auto *pc_27 = values + 30 * nvalues;
    auto *pc_28 = values + 31 * nvalues;
    auto *pc_29 = values + 32 * nvalues;
    auto *pc_30 = values + 33 * nvalues;
    auto *pc_31 = values + 34 * nvalues;
    auto *pc_32 = values + 35 * nvalues;
    auto *pc_33 = values + 36 * nvalues;
    auto *pc_34 = values + 37 * nvalues;
    auto *pc_35 = values + 38 * nvalues;
    auto *pc_36 = values + 42 * nvalues;
    auto *pc_37 = values + 43 * nvalues;
    auto *pc_38 = values + 44 * nvalues;
    auto *pc_39 = values + 45 * nvalues;
    auto *pc_40 = values + 46 * nvalues;
    auto *pc_41 = values + 47 * nvalues;
    auto *pc_42 = values + 48 * nvalues;
    auto *pc_43 = values + 49 * nvalues;
    auto *pc_44 = values + 50 * nvalues;
    auto *pc_45 = values + 51 * nvalues;
    auto *pc_46 = values + 56 * nvalues;
    auto *pc_47 = values + 57 * nvalues;
    auto *pc_48 = values + 58 * nvalues;
    auto *pc_49 = values + 59 * nvalues;
    auto *pc_50 = values + 60 * nvalues;
    auto *pc_51 = values + 61 * nvalues;
    auto *pc_52 = values + 62 * nvalues;
    auto *pc_53 = values + 63 * nvalues;
    auto *pc_54 = values + 64 * nvalues;
    auto *pc_55 = values + 70 * nvalues;
    auto *pc_56 = values + 71 * nvalues;
    auto *pc_57 = values + 72 * nvalues;
    auto *pc_58 = values + 73 * nvalues;
    auto *pc_59 = values + 74 * nvalues;
    auto *pc_60 = values + 75 * nvalues;
    auto *pc_61 = values + 76 * nvalues;
    auto *pc_62 = values + 77 * nvalues;
    auto *pc_63 = values + 84 * nvalues;
    auto *pc_64 = values + 85 * nvalues;
    auto *pc_65 = values + 86 * nvalues;
    auto *pc_66 = values + 87 * nvalues;
    auto *pc_67 = values + 88 * nvalues;
    auto *pc_68 = values + 89 * nvalues;
    auto *pc_69 = values + 90 * nvalues;
    auto *pc_70 = values + 98 * nvalues;
    auto *pc_71 = values + 99 * nvalues;
    auto *pc_72 = values + 100 * nvalues;
    auto *pc_73 = values + 101 * nvalues;
    auto *pc_74 = values + 102 * nvalues;
    auto *pc_75 = values + 103 * nvalues;
    auto *pc_76 = values + 112 * nvalues;
    auto *pc_77 = values + 113 * nvalues;
    auto *pc_78 = values + 114 * nvalues;
    auto *pc_79 = values + 115 * nvalues;
    auto *pc_80 = values + 116 * nvalues;
    auto *pc_81 = values + 126 * nvalues;
    auto *pc_82 = values + 127 * nvalues;
    auto *pc_83 = values + 128 * nvalues;
    auto *pc_84 = values + 129 * nvalues;
    auto *pc_85 = values + 140 * nvalues;
    auto *pc_86 = values + 141 * nvalues;
    auto *pc_87 = values + 142 * nvalues;
    auto *pc_88 = values + 154 * nvalues;
    auto *pc_89 = values + 155 * nvalues;
    auto *pc_90 = values + 168 * nvalues;

    // NOTE: the factors of the terms depend on the angular momenta alone, so they
    // are formed once for the whole matrix instead of once for every atom pair.

    const auto f_10395_16 = 10395.0 / 16.0;
    const auto f_10395_32 = 10395.0 / 32.0;
    const auto f_10395_64 = 10395.0 / 64.0;
    const auto f_10458_96577 = 10458.0 / 96577.0;
    const auto f_1050_2717 = 1050.0 / 2717.0;
    const auto f_105_143 = 105.0 / 143.0;
    const auto f_105_22 = 105.0 / 22.0;
    const auto f_105_247 = 105.0 / 247.0;
    const auto f_105_52 = 105.0 / 52.0;
    const auto f_10_143 = 10.0 / 143.0;
    const auto f_10_209 = 10.0 / 209.0;
    const auto f_120_187 = 120.0 / 187.0;
    const auto f_1215_11 = 1215.0 / 11.0;
    const auto f_1215_26 = 1215.0 / 26.0;
    const auto f_1215_8 = 1215.0 / 8.0;
    const auto f_12285_32 = 12285.0 / 32.0;
    const auto f_13230_96577 = 13230.0 / 96577.0;
    const auto f_135_1 = 135.0;
    const auto f_135_26 = 135.0 / 26.0;
    const auto f_1440_11 = 1440.0 / 11.0;
    const auto f_1485_8 = 1485.0 / 8.0;
    const auto f_14_143 = 14.0 / 143.0;
    const auto f_150_11 = 150.0 / 11.0;
    const auto f_15876_96577 = 15876.0 / 96577.0;
    const auto f_159_2431 = 159.0 / 2431.0;
    const auto f_15_1 = 15.0;
    const auto f_162_2431 = 162.0 / 2431.0;
    const auto f_16335_96577 = 16335.0 / 96577.0;
    const auto f_165_4 = 165.0 / 4.0;
    const auto f_180_1 = 180.0;
    const auto f_1890_11 = 1890.0 / 11.0;
    const auto f_189_323 = 189.0 / 323.0;
    const auto f_189_4199 = 189.0 / 4199.0;
    const auto f_18_221 = 18.0 / 221.0;
    const auto f_192_2431 = 192.0 / 2431.0;
    const auto f_1995_572 = 1995.0 / 572.0;
    const auto f_1_11 = 1.0 / 11.0;
    const auto f_1_13 = 1.0 / 13.0;
    const auto f_200_3553 = 200.0 / 3553.0;
    const auto f_20_323 = 20.0 / 323.0;
    const auto f_210_209 = 210.0 / 209.0;
    const auto f_2160_11 = 2160.0 / 11.0;
    const auto f_2178_96577 = 2178.0 / 96577.0;
    const auto f_2385_16 = 2385.0 / 16.0;
    const auto f_2385_22 = 2385.0 / 22.0;
    const auto f_252_2431 = 252.0 / 2431.0;
    const auto f_25_247 = 25.0 / 247.0;
    const auto f_26136_96577 = 26136.0 / 96577.0;
    const auto f_270_1 = 270.0;
    const auto f_27_221 = 27.0 / 221.0;
    const auto f_288_2431 = 288.0 / 2431.0;
    const auto f_2_13 = 2.0 / 13.0;
    const auto f_2_143 = 2.0 / 143.0;
    const auto f_300_11 = 300.0 / 11.0;
    const auto f_300_187 = 300.0 / 187.0;
    const auto f_30492_96577 = 30492.0 / 96577.0;
    const auto f_3087_8398 = 3087.0 / 8398.0;
    const auto f_3087_96577 = 3087.0 / 96577.0;
    const auto f_30_17 = 30.0 / 17.0;
    const auto f_318_143 = 318.0 / 143.0;
    const auto f_3225_44 = 3225.0 / 44.0;
    const auto f_324_143 = 324.0 / 143.0;
    const auto f_3375_16 = 3375.0 / 16.0;
    const auto f_3375_8 = 3375.0 / 8.0;
    const auto f_33_96577 = 33.0 / 96577.0;
    const auto f_350_2717 = 350.0 / 2717.0;
    const auto f_355_2717 = 355.0 / 2717.0;
    const auto f_3645_143 = 3645.0 / 143.0;
    const auto f_3675_286 = 3675.0 / 286.0;
    const auto f_36_13 = 36.0 / 13.0;
    const auto f_375_11 = 375.0 / 11.0;
    const auto f_375_143 = 375.0 / 143.0;
    const auto f_375_2 = 375.0 / 2.0;
    const auto f_375_22 = 375.0 / 22.0;
    const auto f_375_286 = 375.0 / 286.0;
    const auto f_375_4 = 375.0 / 4.0;
    const auto f_378_7429 = 378.0 / 7429.0;
    const auto f_378_96577 = 378.0 / 96577.0;
    const auto f_384_143 = 384.0 / 143.0;
    const auto f_396_96577 = 396.0 / 96577.0;
    const auto f_3_1 = 3.0;
    const auto f_3_221 = 3.0 / 221.0;
    const auto f_400_3553 = 400.0 / 3553.0;
    const auto f_405_13 = 405.0 / 13.0;
    const auto f_405_2 = 405.0 / 2.0;
    const auto f_430_3553 = 430.0 / 3553.0;
    const auto f_4320_143 = 4320.0 / 143.0;
    const auto f_4455_16 = 4455.0 / 16.0;
    const auto f_445_2717 = 445.0 / 2717.0;
    const auto f_45_2 = 45.0 / 2.0;
    const auto f_4725_16 = 4725.0 / 16.0;
    const auto f_4725_32 = 4725.0 / 32.0;
    const auto f_4725_8 = 4725.0 / 8.0;
    const auto f_495_16 = 495.0 / 16.0;
    const auto f_495_2 = 495.0 / 2.0;
    const auto f_504_143 = 504.0 / 143.0;
    const auto f_50_2717 = 50.0 / 2717.0;
    const auto f_50_323 = 50.0 / 323.0;
    const auto f_5229_4199 = 5229.0 / 4199.0;
    const auto f_525_11 = 525.0 / 11.0;
    const auto f_525_143 = 525.0 / 143.0;
    const auto f_525_2 = 525.0 / 2.0;
    const auto f_525_247 = 525.0 / 247.0;
    const auto f_525_286 = 525.0 / 286.0;
    const auto f_525_52 = 525.0 / 52.0;
    const auto f_54_13 = 54.0 / 13.0;
    const auto f_5670_143 = 5670.0 / 143.0;
    const auto f_576_143 = 576.0 / 143.0;
    const auto f_5_143 = 5.0 / 143.0;
    const auto f_5_247 = 5.0 / 247.0;
    const auto f_600_187 = 600.0 / 187.0;
    const auto f_60_11 = 60.0 / 11.0;
    const auto f_645_187 = 645.0 / 187.0;
    const auto f_645_22 = 645.0 / 22.0;
    const auto f_6480_143 = 6480.0 / 143.0;
    const auto f_6615_16 = 6615.0 / 16.0;
    const auto f_6615_4199 = 6615.0 / 4199.0;
    const auto f_675_8 = 675.0 / 8.0;
    const auto f_6993_8398 = 6993.0 / 8398.0;
    const auto f_6993_96577 = 6993.0 / 96577.0;
    const auto f_6_13 = 6.0 / 13.0;
    const auto f_7155_286 = 7155.0 / 286.0;
    const auto f_7260_96577 = 7260.0 / 96577.0;
    const auto f_7350_2717 = 7350.0 / 2717.0;
    const auto f_7425_16 = 7425.0 / 16.0;
    const auto f_7425_8 = 7425.0 / 8.0;
    const auto f_7455_2717 = 7455.0 / 2717.0;
    const auto f_7455_572 = 7455.0 / 572.0;
    const auto f_750_11 = 750.0 / 11.0;
    const auto f_75_1 = 75.0;
    const auto f_75_11 = 75.0 / 11.0;
    const auto f_75_13 = 75.0 / 13.0;
    const auto f_75_143 = 75.0 / 143.0;
    const auto f_75_17 = 75.0 / 17.0;
    const auto f_75_2 = 75.0 / 2.0;
    const auto f_75_22 = 75.0 / 22.0;
    const auto f_75_26 = 75.0 / 26.0;
    const auto f_7938_4199 = 7938.0 / 4199.0;
    const auto f_80_3553 = 80.0 / 3553.0;
    const auto f_825_2 = 825.0 / 2.0;
    const auto f_825_4 = 825.0 / 4.0;
    const auto f_8775_16 = 8775.0 / 16.0;
    const auto f_9345_2717 = 9345.0 / 2717.0;
    const auto f_9345_572 = 9345.0 / 572.0;
    const auto f_945_16 = 945.0 / 16.0;
    const auto f_945_4 = 945.0 / 4.0;
    const auto f_945_442 = 945.0 / 442.0;
    const auto f_945_5083 = 945.0 / 5083.0;
    const auto f_975_22 = 975.0 / 22.0;
    const auto f_975_4 = 975.0 / 4.0;
    const auto fs_1035_32_6 = std::sqrt(3213675.0 / 512.0);
    const auto fs_1035_44_6 = std::sqrt(3213675.0 / 968.0);
    const auto fs_1050_2717_21 = std::sqrt(23152500.0 / 7382089.0);
    const auto fs_105_104_130 = std::sqrt(55125.0 / 416.0);
    const auto fs_105_104_182 = std::sqrt(77175.0 / 416.0);
    const auto fs_105_104_42 = std::sqrt(231525.0 / 5408.0);
    const auto fs_105_1144_770 = std::sqrt(385875.0 / 59488.0);
    const auto fs_105_143_21 = std::sqrt(231525.0 / 20449.0);
    const auto fs_105_143_7 = std::sqrt(77175.0 / 20449.0);
    const auto fs_105_247_13 = std::sqrt(11025.0 / 4693.0);
    const auto fs_105_247_21 = std::sqrt(231525.0 / 61009.0);
    const auto fs_105_247_3 = std::sqrt(33075.0 / 61009.0);
    const auto fs_105_247_35 = std::sqrt(385875.0 / 61009.0);
    const auto fs_105_247_39 = std::sqrt(33075.0 / 4693.0);
    const auto fs_105_247_91 = std::sqrt(77175.0 / 4693.0);
    const auto fs_105_26_7 = std::sqrt(77175.0 / 676.0);
    const auto fs_105_2717_1001 = std::sqrt(77175.0 / 51623.0);
    const auto fs_105_2717_11 = std::sqrt(11025.0 / 671099.0);
    const auto fs_105_2717_1155 = std::sqrt(1157625.0 / 671099.0);
    const auto fs_105_2717_3003 = std::sqrt(231525.0 / 51623.0);
    const auto fs_105_2717_385 = std::sqrt(385875.0 / 671099.0);
    const auto fs_105_2717_715 = std::sqrt(55125.0 / 51623.0);
    const auto fs_105_2717_858 = std::sqrt(66150.0 / 51623.0);
    const auto fs_105_286_858 = std::sqrt(33075.0 / 286.0);
    const auto fs_105_494_130 = std::sqrt(55125.0 / 9386.0);
    const auto fs_105_494_182 = std::sqrt(77175.0 / 9386.0);
    const auto fs_105_494_42 = std::sqrt(231525.0 / 122018.0);
    const auto fs_105_52_13 = std::sqrt(11025.0 / 208.0);
    const auto fs_105_52_21 = std::sqrt(231525.0 / 2704.0);
    const auto fs_105_52_3 = std::sqrt(33075.0 / 2704.0);
    const auto fs_105_52_35 = std::sqrt(385875.0 / 2704.0);
    const auto fs_105_52_39 = std::sqrt(33075.0 / 208.0);
    const auto fs_105_52_91 = std::sqrt(77175.0 / 208.0);
    const auto fs_105_5434_770 = std::sqrt(385875.0 / 1342198.0);
    const auto fs_105_572_1001 = std::sqrt(77175.0 / 2288.0);
    const auto fs_105_572_1155 = std::sqrt(1157625.0 / 29744.0);
    const auto fs_105_572_3003 = std::sqrt(231525.0 / 2288.0);
    const auto fs_105_572_385 = std::sqrt(385875.0 / 29744.0);
    const auto fs_105_572_715 = std::sqrt(55125.0 / 2288.0);
    const auto fs_105_572_858 = std::sqrt(33075.0 / 1144.0);
    const auto fs_1089_96577_195 = std::sqrt(17788815.0 / 717470533.0);
    const auto fs_10_247_7 = std::sqrt(700.0 / 61009.0);
    const auto fs_10_2717_858 = std::sqrt(600.0 / 51623.0);
    const auto fs_10_323_7 = std::sqrt(700.0 / 104329.0);
    const auto fs_10_3553_154 = std::sqrt(1400.0 / 1147619.0);
    const auto fs_1134_96577_221 = std::sqrt(1285956.0 / 42204149.0);
    const auto fs_120_3553_7 = std::sqrt(100800.0 / 12623809.0);
    const auto fs_1215_572_154 = std::sqrt(10333575.0 / 14872.0);
    const auto fs_1215_572_330 = std::sqrt(22143375.0 / 14872.0);
    const auto fs_1260_2717_22 = std::sqrt(3175200.0 / 671099.0);
    const auto fs_126_143_14 = std::sqrt(222264.0 / 20449.0);
    const auto fs_126_143_3 = std::sqrt(47628.0 / 20449.0);
    const auto fs_126_2431_3 = std::sqrt(47628.0 / 5909761.0);
    const auto fs_126_2431_5 = std::sqrt(79380.0 / 5909761.0);
    const auto fs_126_4199_2145 = std::sqrt(2619540.0 / 1356277.0);
    const auto fs_126_4199_3315 = std::sqrt(238140.0 / 79781.0);
    const auto fs_126_7429_65 = std::sqrt(1031940.0 / 55190041.0);
    const auto fs_126_96577_12155 = std::sqrt(873180.0 / 42204149.0);
    const auto fs_126_96577_12597 = std::sqrt(47628.0 / 2221271.0);
    const auto fs_126_96577_143 = std::sqrt(174636.0 / 717470533.0);
    const auto fs_126_96577_15015 = std::sqrt(18336780.0 / 717470533.0);
    const auto fs_12_143_210 = std::sqrt(30240.0 / 20449.0);
    const auto fs_1323_8398_39 = std::sqrt(5250987.0 / 5425108.0);
    const auto fs_1323_96577_39 = std::sqrt(5250987.0 / 717470533.0);
    const auto fs_132_96577_20995 = std::sqrt(87120.0 / 2221271.0);
    const auto fs_135_143_210 = std::sqrt(3827250.0 / 20449.0);
    const auto fs_135_16_55 = std::sqrt(1002375.0 / 256.0);
    const auto fs_135_16_77 = std::sqrt(1403325.0 / 256.0);
    const auto fs_135_22_55 = std::sqrt(91125.0 / 44.0);
    const auto fs_135_22_77 = std::sqrt(127575.0 / 44.0);
    const auto fs_135_26_42 = std::sqrt(382725.0 / 338.0);
    const auto fs_135_286_2310 = std::sqrt(1913625.0 / 3718.0);
    const auto fs_135_4_2 = std::sqrt(18225.0 / 8.0);
    const auto fs_135_4_30 = std::sqrt(273375.0 / 8.0);
    const auto fs_135_4_33 = std::sqrt(601425.0 / 16.0);
    const auto fs_135_5434_10 = std::sqrt(91125.0 / 14764178.0);
    const auto fs_1470_2717_22 = std::sqrt(4321800.0 / 671099.0);
    const auto fs_1485_32_30 = std::sqrt(33078375.0 / 512.0);
    const auto fs_15_11_210 = std::sqrt(47250.0 / 121.0);
    const auto fs_15_11_385 = std::sqrt(7875.0 / 11.0);
    const auto fs_15_11_462 = std::sqrt(9450.0 / 11.0);
    const auto fs_15_17_7 = std::sqrt(1575.0 / 289.0);
    const auto fs_15_187_154 = std::sqrt(3150.0 / 3179.0);
    const auto fs_15_22_154 = std::sqrt(1575.0 / 22.0);
    const auto fs_15_2431_210 = std::sqrt(47250.0 / 5909761.0);
    const auto fs_15_2717_462 = std::sqrt(9450.0 / 671099.0);
    const auto fs_15_2717_715 = std::sqrt(1125.0 / 51623.0);
    const auto fs_15_2_7 = std::sqrt(1575.0 / 4.0);
    const auto fs_15_3553_462 = std::sqrt(9450.0 / 1147619.0);
    const auto fs_15_5434_462 = std::sqrt(4725.0 / 1342198.0);
    const auto fs_165_96577_17017 = std::sqrt(2096325.0 / 42204149.0);
    const auto fs_168_143_5 = std::sqrt(141120.0 / 20449.0);
    const auto fs_1701_8398_55 = std::sqrt(159137055.0 / 70526404.0);
    const auto fs_1701_96577_55 = std::sqrt(159137055.0 / 9327116929.0);
    const auto fs_1755_32_10 = std::sqrt(15400125.0 / 512.0);
    const auto fs_1755_44_10 = std::sqrt(15400125.0 / 968.0);
    const auto fs_180_187_7 = std::sqrt(226800.0 / 34969.0);
    const auto fs_1815_193154_78 = std::sqrt(9882675.0 / 1434941066.0);
    const auto fs_1890_143_5 = std::sqrt(17860500.0 / 20449.0);
    const auto fs_189_16796_858 = std::sqrt(1178793.0 / 10850216.0);
    const auto fs_189_193154_858 = std::sqrt(1178793.0 / 1434941066.0);
    const auto fs_189_4199_1155 = std::sqrt(41257755.0 / 17631601.0);
    const auto fs_189_4199_154 = std::sqrt(5501034.0 / 17631601.0);
    const auto fs_189_4199_78 = std::sqrt(214326.0 / 1356277.0);
    const auto fs_189_8398_1001 = std::sqrt(2750517.0 / 5425108.0);
    const auto fs_189_8398_2431 = std::sqrt(392931.0 / 319124.0);
    const auto fs_189_8398_4290 = std::sqrt(5893965.0 / 2712554.0);
    const auto fs_189_8398_8398 = std::sqrt(35721.0 / 8398.0);
    const auto fs_189_96577_1001 = std::sqrt(2750517.0 / 717470533.0);
    const auto fs_189_96577_2431 = std::sqrt(392931.0 / 42204149.0);
    const auto fs_189_96577_4290 = std::sqrt(11787930.0 / 717470533.0);
    const auto fs_189_96577_8398 = std::sqrt(71442.0 / 2221271.0);
    const auto fs_18_143_55 = std::sqrt(1620.0 / 1859.0);
    const auto fs_18_143_77 = std::sqrt(2268.0 / 1859.0);
    const auto fs_1953_8398_26 = std::sqrt(3814209.0 / 2712554.0);
    const auto fs_1953_96577_26 = std::sqrt(7628418.0 / 717470533.0);
    const auto fs_1995_572_7 = std::sqrt(27860175.0 / 327184.0);
    const auto fs_1_143_22 = std::sqrt(2.0 / 1859.0);
    const auto fs_1_143_30 = std::sqrt(30.0 / 20449.0);
    const auto fs_1_143_55 = std::sqrt(5.0 / 1859.0);
    const auto fs_1_143_7 = std::sqrt(7.0 / 20449.0);
    const auto fs_2025_16_10 = std::sqrt(20503125.0 / 128.0);
    const auto fs_2025_32_66 = std::sqrt(135320625.0 / 512.0);
    const auto fs_20_2717_21 = std::sqrt(8400.0 / 7382089.0);
    const auto fs_20_3553_210 = std::sqrt(84000.0 / 12623809.0);
    const auto fs_20_3553_385 = std::sqrt(14000.0 / 1147619.0);
    const auto fs_20_3553_462 = std::sqrt(16800.0 / 1147619.0);
    const auto fs_210_247_7 = std::sqrt(308700.0 / 61009.0);
    const auto fs_210_2717_858 = std::sqrt(264600.0 / 51623.0);
    const auto fs_2178_96577_182 = std::sqrt(66411576.0 / 717470533.0);
    const auto fs_21_143_462 = std::sqrt(18522.0 / 1859.0);
    const auto fs_21_2431_7 = std::sqrt(3087.0 / 5909761.0);
    const auto fs_21_4862_462 = std::sqrt(9261.0 / 1074502.0);
    const auto fs_2205_16796_130 = std::sqrt(24310125.0 / 10850216.0);
    const auto fs_2205_193154_130 = std::sqrt(24310125.0 / 1434941066.0);
    const auto fs_2205_2717_11 = std::sqrt(4862025.0 / 671099.0);
    const auto fs_2205_572_11 = std::sqrt(4862025.0 / 29744.0);
    const auto fs_225_11_7 = std::sqrt(354375.0 / 121.0);
    const auto fs_225_16_210 = std::sqrt(5315625.0 / 128.0);
    const auto fs_225_22_10 = std::sqrt(253125.0 / 242.0);
    const auto fs_225_22_210 = std::sqrt(5315625.0 / 242.0);
    const auto fs_225_286_10 = std::sqrt(253125.0 / 40898.0);
    const auto fs_225_44_21 = std::sqrt(1063125.0 / 1936.0);
    const auto fs_225_44_66 = std::sqrt(151875.0 / 88.0);
    const auto fs_225_4_10 = std::sqrt(253125.0 / 8.0);
    const auto fs_225_572_66 = std::sqrt(151875.0 / 14872.0);
    const auto fs_225_88_462 = std::sqrt(1063125.0 / 352.0);
    const auto fs_225_8_66 = std::sqrt(1670625.0 / 32.0);
    const auto fs_252_143_3 = std::sqrt(190512.0 / 20449.0);
    const auto fs_252_143_5 = std::sqrt(317520.0 / 20449.0);
    const auto fs_252_96577_2145 = std::sqrt(10478160.0 / 717470533.0);
    const auto fs_252_96577_3315 = std::sqrt(952560.0 / 42204149.0);
    const auto fs_2583_16796_66 = std::sqrt(220172337.0 / 141052808.0);
    const auto fs_2583_193154_66 = std::sqrt(220172337.0 / 18654233858.0);
    const auto fs_25_3553_70 = std::sqrt(43750.0 / 12623809.0);
    const auto fs_264_96577_273 = std::sqrt(1463616.0 / 717470533.0);
    const auto fs_264_96577_5005 = std::sqrt(26832960.0 / 717470533.0);
    const auto fs_2709_16796_10 = std::sqrt(36693405.0 / 141052808.0);
    const auto fs_2709_193154_10 = std::sqrt(36693405.0 / 18654233858.0);
    const auto fs_2709_8398_22 = std::sqrt(80725491.0 / 35263202.0);
    const auto fs_2709_96577_22 = std::sqrt(161450982.0 / 9327116929.0);
    const auto fs_270_11_2 = std::sqrt(145800.0 / 121.0);
    const auto fs_270_11_33 = std::sqrt(218700.0 / 11.0);
    const auto fs_27_143_154 = std::sqrt(10206.0 / 1859.0);
    const auto fs_27_143_330 = std::sqrt(21870.0 / 1859.0);
    const auto fs_27_4862_154 = std::sqrt(5103.0 / 1074502.0);
    const auto fs_27_4862_330 = std::sqrt(10935.0 / 1074502.0);
    const auto fs_2835_1144_10 = std::sqrt(40186125.0 / 654368.0);
    const auto fs_2835_143_3 = std::sqrt(24111675.0 / 20449.0);
    const auto fs_2835_143_5 = std::sqrt(40186125.0 / 20449.0);
    const auto fs_2835_286_14 = std::sqrt(56260575.0 / 40898.0);
    const auto fs_2835_286_3 = std::sqrt(24111675.0 / 81796.0);
    const auto fs_2835_32_10 = std::sqrt(40186125.0 / 512.0);
    const auto fs_2835_5434_10 = std::sqrt(40186125.0 / 14764178.0);
    const auto fs_2835_64_66 = std::sqrt(265228425.0 / 2048.0);
    const auto fs_297_96577_5005 = std::sqrt(33960465.0 / 717470533.0);
    const auto fs_2_143_30 = std::sqrt(120.0 / 20449.0);
    const auto fs_2_143_70 = std::sqrt(280.0 / 20449.0);
    const auto fs_30_143_210 = std::sqrt(189000.0 / 20449.0);
    const auto fs_30_187_210 = std::sqrt(189000.0 / 34969.0);
    const auto fs_30_187_385 = std::sqrt(31500.0 / 3179.0);
    const auto fs_30_187_462 = std::sqrt(37800.0 / 3179.0);
    const auto fs_30_2717_70 = std::sqrt(63000.0 / 7382089.0);
    const auto fs_30_2717_77 = std::sqrt(6300.0 / 671099.0);
    const auto fs_30_3553_21 = std::sqrt(18900.0 / 12623809.0);
    const auto fs_3105_572_6 = std::sqrt(28923075.0 / 163592.0);
    const auto fs_315_1144_462 = std::sqrt(2083725.0 / 59488.0);
    const auto fs_315_11_35 = std::sqrt(3472875.0 / 121.0);
    const auto fs_315_143_22 = std::sqrt(198450.0 / 1859.0);
    const auto fs_315_16_7 = std::sqrt(694575.0 / 256.0);
    const auto fs_315_22_7 = std::sqrt(694575.0 / 484.0);
    const auto fs_315_2717_462 = std::sqrt(4167450.0 / 671099.0);
    const auto fs_315_2717_715 = std::sqrt(496125.0 / 51623.0);
    const auto fs_315_286_70 = std::sqrt(3472875.0 / 40898.0);
    const auto fs_315_286_77 = std::sqrt(694575.0 / 7436.0);
    const auto fs_315_32_462 = std::sqrt(22920975.0 / 512.0);
    const auto fs_315_44_462 = std::sqrt(2083725.0 / 88.0);
    const auto fs_315_4_5 = std::sqrt(496125.0 / 16.0);
    const auto fs_315_5434_462 = std::sqrt(2083725.0 / 1342198.0);
    const auto fs_315_572_462 = std::sqrt(2083725.0 / 14872.0);
    const auto fs_315_572_715 = std::sqrt(496125.0 / 2288.0);
    const auto fs_315_8398_858 = std::sqrt(3274425.0 / 2712554.0);
    const auto fs_315_8_35 = std::sqrt(3472875.0 / 64.0);
    const auto fs_315_96577_858 = std::sqrt(6548850.0 / 717470533.0);
    const auto fs_330_96577_4862 = std::sqrt(2395800.0 / 42204149.0);
    const auto fs_3375_16_3 = std::sqrt(34171875.0 / 256.0);
    const auto fs_33_193154_182 = std::sqrt(7623.0 / 1434941066.0);
    const auto fs_33_193154_26 = std::sqrt(1089.0 / 1434941066.0);
    const auto fs_33_193154_910 = std::sqrt(38115.0 / 1434941066.0);
    const auto fs_33_96577_1352078 = std::sqrt(15246.0 / 96577.0);
    const auto fs_33_96577_146965 = std::sqrt(38115.0 / 2221271.0);
    const auto fs_33_96577_25194 = std::sqrt(6534.0 / 2221271.0);
    const auto fs_33_96577_3094 = std::sqrt(15246.0 / 42204149.0);
    const auto fs_33_96577_323323 = std::sqrt(83853.0 / 2221271.0);
    const auto fs_33_96577_429 = std::sqrt(35937.0 / 717470533.0);
    const auto fs_33_96577_461890 = std::sqrt(119790.0 / 2221271.0);
    const auto fs_33_96577_62985 = std::sqrt(16335.0 / 2221271.0);
    const auto fs_33_96577_676039 = std::sqrt(7623.0 / 96577.0);
    const auto fs_33_96577_910 = std::sqrt(76230.0 / 717470533.0);
    const auto fs_35_2717_3 = std::sqrt(3675.0 / 7382089.0);
    const auto fs_363_96577_65 = std::sqrt(658845.0 / 717470533.0);
    const auto fs_36_2431_2 = std::sqrt(2592.0 / 5909761.0);
    const auto fs_36_2431_33 = std::sqrt(3888.0 / 537251.0);
    const auto fs_375_22_3 = std::sqrt(421875.0 / 484.0);
    const auto fs_375_286_3 = std::sqrt(421875.0 / 81796.0);
    const auto fs_375_4_3 = std::sqrt(421875.0 / 16.0);
    const auto fs_375_88_70 = std::sqrt(4921875.0 / 3872.0);
    const auto fs_378_4199_22 = std::sqrt(3143448.0 / 17631601.0);
    const auto fs_378_4199_455 = std::sqrt(5000940.0 / 1356277.0);
    const auto fs_378_96577_1155 = std::sqrt(165031020.0 / 9327116929.0);
    const auto fs_378_96577_154 = std::sqrt(22004136.0 / 9327116929.0);
    const auto fs_378_96577_78 = std::sqrt(857304.0 / 717470533.0);
    const auto fs_396_96577_3003 = std::sqrt(36224496.0 / 717470533.0);
    const auto fs_396_96577_663 = std::sqrt(470448.0 / 42204149.0);
    const auto fs_3_143_10 = std::sqrt(90.0 / 20449.0);
    const auto fs_3_221_42 = std::sqrt(378.0 / 48841.0);
    const auto fs_3_2431_2310 = std::sqrt(1890.0 / 537251.0);
    const auto fs_3_286_66 = std::sqrt(27.0 / 3718.0);
    const auto fs_405_286_55 = std::sqrt(820125.0 / 7436.0);
    const auto fs_405_286_77 = std::sqrt(1148175.0 / 7436.0);
    const auto fs_405_32_154 = std::sqrt(12629925.0 / 512.0);
    const auto fs_405_32_330 = std::sqrt(27064125.0 / 512.0);
    const auto fs_405_44_10 = std::sqrt(820125.0 / 968.0);
    const auto fs_405_44_154 = std::sqrt(1148175.0 / 88.0);
    const auto fs_405_44_330 = std::sqrt(2460375.0 / 88.0);
    const auto fs_405_52_30 = std::sqrt(2460375.0 / 1352.0);
    const auto fs_420_2717_21 = std::sqrt(3704400.0 / 7382089.0);
    const auto fs_42_143_7 = std::sqrt(12348.0 / 20449.0);
    const auto fs_42_2431_35 = std::sqrt(61740.0 / 5909761.0);
    const auto fs_441_4199_165 = std::sqrt(32089365.0 / 17631601.0);
    const auto fs_441_4199_39 = std::sqrt(583443.0 / 1356277.0);
    const auto fs_441_8398_429 = std::sqrt(6417873.0 / 5425108.0);
    const auto fs_441_8398_442 = std::sqrt(194481.0 / 159562.0);
    const auto fs_441_96577_429 = std::sqrt(6417873.0 / 717470533.0);
    const auto fs_441_96577_442 = std::sqrt(388962.0 / 42204149.0);
    const auto fs_45_11_210 = std::sqrt(425250.0 / 121.0);
    const auto fs_45_16_2310 = std::sqrt(2338875.0 / 128.0);
    const auto fs_45_187_21 = std::sqrt(42525.0 / 34969.0);
    const auto fs_45_22_21 = std::sqrt(42525.0 / 484.0);
    const auto fs_45_22_2310 = std::sqrt(212625.0 / 22.0);
    const auto fs_45_2717_30 = std::sqrt(60750.0 / 7382089.0);
    const auto fs_45_2_42 = std::sqrt(42525.0 / 2.0);
    const auto fs_45_374_462 = std::sqrt(42525.0 / 6358.0);
    const auto fs_45_44_462 = std::sqrt(42525.0 / 88.0);
    const auto fs_45_8_210 = std::sqrt(212625.0 / 32.0);
    const auto fs_462_96577_2145 = std::sqrt(35218260.0 / 717470533.0);
    const auto fs_462_96577_442 = std::sqrt(426888.0 / 42204149.0);
    const auto fs_4725_16_3 = std::sqrt(66976875.0 / 256.0);
    const auto fs_4725_32_10 = std::sqrt(111628125.0 / 512.0);
    const auto fs_4725_32_3 = std::sqrt(66976875.0 / 1024.0);
    const auto fs_495_16_42 = std::sqrt(5145525.0 / 128.0);
    const auto fs_495_96577_2002 = std::sqrt(37733850.0 / 717470533.0);
    const auto fs_50_2717_21 = std::sqrt(52500.0 / 7382089.0);
    const auto fs_525_22_3 = std::sqrt(826875.0 / 484.0);
    const auto fs_525_286_21 = std::sqrt(5788125.0 / 81796.0);
    const auto fs_525_286_3 = std::sqrt(826875.0 / 81796.0);
    const auto fs_525_44_10 = std::sqrt(1378125.0 / 968.0);
    const auto fs_525_4_3 = std::sqrt(826875.0 / 16.0);
    const auto fs_525_572_10 = std::sqrt(1378125.0 / 163592.0);
    const auto fs_525_8_10 = std::sqrt(1378125.0 / 32.0);
    const auto fs_567_4199_221 = std::sqrt(321489.0 / 79781.0);
    const auto fs_567_8398_10 = std::sqrt(1607445.0 / 35263202.0);
    const auto fs_567_96577_10 = std::sqrt(3214890.0 / 9327116929.0);
    const auto fs_5_143_3 = std::sqrt(75.0 / 20449.0);
    const auto fs_5_143_7 = std::sqrt(175.0 / 20449.0);
    const auto fs_5_247_13 = std::sqrt(25.0 / 4693.0);
    const auto fs_5_247_21 = std::sqrt(525.0 / 61009.0);
    const auto fs_5_247_3 = std::sqrt(75.0 / 61009.0);
    const auto fs_5_247_35 = std::sqrt(875.0 / 61009.0);
    const auto fs_5_247_39 = std::sqrt(75.0 / 4693.0);
    const auto fs_5_247_91 = std::sqrt(175.0 / 4693.0);
    const auto fs_5_2717_1001 = std::sqrt(175.0 / 51623.0);
    const auto fs_5_2717_1155 = std::sqrt(2625.0 / 671099.0);
    const auto fs_5_2717_3003 = std::sqrt(525.0 / 51623.0);
    const auto fs_5_2717_385 = std::sqrt(875.0 / 671099.0);
    const auto fs_5_2717_715 = std::sqrt(125.0 / 51623.0);
    const auto fs_5_2717_858 = std::sqrt(150.0 / 51623.0);
    const auto fs_5_494_130 = std::sqrt(125.0 / 9386.0);
    const auto fs_5_494_182 = std::sqrt(175.0 / 9386.0);
    const auto fs_5_494_42 = std::sqrt(525.0 / 122018.0);
    const auto fs_5_5434_770 = std::sqrt(875.0 / 1342198.0);
    const auto fs_60_2717_22 = std::sqrt(7200.0 / 671099.0);
    const auto fs_630_11_5 = std::sqrt(1984500.0 / 121.0);
    const auto fs_630_2717_70 = std::sqrt(27783000.0 / 7382089.0);
    const auto fs_630_2717_77 = std::sqrt(2778300.0 / 671099.0);
    const auto fs_63_2431_14 = std::sqrt(55566.0 / 5909761.0);
    const auto fs_63_2431_3 = std::sqrt(11907.0 / 5909761.0);
    const auto fs_63_323_65 = std::sqrt(257985.0 / 104329.0);
    const auto fs_63_4199_12155 = std::sqrt(218295.0 / 79781.0);
    const auto fs_63_4199_12597 = std::sqrt(11907.0 / 4199.0);
    const auto fs_63_4199_143 = std::sqrt(43659.0 / 1356277.0);
    const auto fs_63_4199_15015 = std::sqrt(4584195.0 / 1356277.0);
    const auto fs_63_442_273 = std::sqrt(83349.0 / 15028.0);
    const auto fs_63_442_66 = std::sqrt(130977.0 / 97682.0);
    const auto fs_63_5083_273 = std::sqrt(83349.0 / 1987453.0);
    const auto fs_63_5083_66 = std::sqrt(261954.0 / 25836889.0);
    const auto fs_63_8398_165 = std::sqrt(654885.0 / 70526404.0);
    const auto fs_63_8398_20995 = std::sqrt(19845.0 / 16796.0);
    const auto fs_63_8398_663 = std::sqrt(11907.0 / 319124.0);
    const auto fs_63_8398_92378 = std::sqrt(43659.0 / 8398.0);
    const auto fs_63_96577_165 = std::sqrt(654885.0 / 9327116929.0);
    const auto fs_63_96577_20995 = std::sqrt(19845.0 / 2221271.0);
    const auto fs_63_96577_663 = std::sqrt(11907.0 / 42204149.0);
    const auto fs_63_96577_92378 = std::sqrt(87318.0 / 2221271.0);
    const auto fs_6615_32_3 = std::sqrt(131274675.0 / 1024.0);
    const auto fs_6615_64_10 = std::sqrt(218791125.0 / 2048.0);
    const auto fs_66_96577_138567 = std::sqrt(143748.0 / 2221271.0);
    const auto fs_66_96577_176358 = std::sqrt(182952.0 / 2221271.0);
    const auto fs_66_96577_30030 = std::sqrt(10062360.0 / 717470533.0);
    const auto fs_66_96577_455 = std::sqrt(152460.0 / 717470533.0);
    const auto fs_66_96577_4641 = std::sqrt(91476.0 / 42204149.0);
    const auto fs_66_96577_51051 = std::sqrt(1006236.0 / 42204149.0);
    const auto fs_675_16_22 = std::sqrt(5011875.0 / 128.0);
    const auto fs_675_16_30 = std::sqrt(6834375.0 / 128.0);
    const auto fs_675_16_55 = std::sqrt(25059375.0 / 256.0);
    const auto fs_675_16_7 = std::sqrt(3189375.0 / 256.0);
    const auto fs_675_286_210 = std::sqrt(47840625.0 / 40898.0);
    const auto fs_675_8_30 = std::sqrt(6834375.0 / 32.0);
    const auto fs_675_8_70 = std::sqrt(15946875.0 / 32.0);
    const auto fs_693_16796_1430 = std::sqrt(26413695.0 / 10850216.0);
    const auto fs_693_16796_286 = std::sqrt(5282739.0 / 10850216.0);
    const auto fs_693_193154_1430 = std::sqrt(26413695.0 / 1434941066.0);
    const auto fs_693_193154_286 = std::sqrt(5282739.0 / 1434941066.0);
    const auto fs_693_96577_1430 = std::sqrt(52827390.0 / 717470533.0);
    const auto fs_69_143_6 = std::sqrt(28566.0 / 20449.0);
    const auto fs_69_4862_6 = std::sqrt(14283.0 / 11819522.0);
    const auto fs_6_13_42 = std::sqrt(1512.0 / 169.0);
    const auto fs_6_143_2310 = std::sqrt(7560.0 / 1859.0);
    const auto fs_6_2431_210 = std::sqrt(7560.0 / 5909761.0);
    const auto fs_70_2717_22 = std::sqrt(9800.0 / 671099.0);
    const auto fs_72_143_2 = std::sqrt(10368.0 / 20449.0);
    const auto fs_72_143_33 = std::sqrt(15552.0 / 1859.0);
    const auto fs_735_2717_3 = std::sqrt(1620675.0 / 7382089.0);
    const auto fs_735_286_22 = std::sqrt(540225.0 / 3718.0);
    const auto fs_735_572_3 = std::sqrt(1620675.0 / 327184.0);
    const auto fs_756_96577_22 = std::sqrt(12573792.0 / 9327116929.0);
    const auto fs_756_96577_455 = std::sqrt(20003760.0 / 717470533.0);
    const auto fs_75_11_30 = std::sqrt(168750.0 / 121.0);
    const auto fs_75_11_70 = std::sqrt(393750.0 / 121.0);
    const auto fs_75_143_30 = std::sqrt(168750.0 / 20449.0);
    const auto fs_75_143_70 = std::sqrt(393750.0 / 20449.0);
    const auto fs_75_22_210 = std::sqrt(590625.0 / 242.0);
    const auto fs_75_22_22 = std::sqrt(5625.0 / 22.0);
    const auto fs_75_22_30 = std::sqrt(84375.0 / 242.0);
    const auto fs_75_22_385 = std::sqrt(196875.0 / 44.0);
    const auto fs_75_22_462 = std::sqrt(118125.0 / 22.0);
    const auto fs_75_22_55 = std::sqrt(28125.0 / 44.0);
    const auto fs_75_22_7 = std::sqrt(39375.0 / 484.0);
    const auto fs_75_286_22 = std::sqrt(5625.0 / 3718.0);
    const auto fs_75_286_30 = std::sqrt(84375.0 / 40898.0);
    const auto fs_75_286_55 = std::sqrt(28125.0 / 7436.0);
    const auto fs_75_286_7 = std::sqrt(39375.0 / 81796.0);
    const auto fs_75_2_30 = std::sqrt(84375.0 / 2.0);
    const auto fs_75_2_70 = std::sqrt(196875.0 / 2.0);
    const auto fs_75_374_70 = std::sqrt(196875.0 / 69938.0);
    const auto fs_75_44_154 = std::sqrt(39375.0 / 88.0);
    const auto fs_75_44_70 = std::sqrt(196875.0 / 968.0);
    const auto fs_75_4_22 = std::sqrt(61875.0 / 8.0);
    const auto fs_75_4_30 = std::sqrt(84375.0 / 8.0);
    const auto fs_75_4_55 = std::sqrt(309375.0 / 16.0);
    const auto fs_75_4_7 = std::sqrt(39375.0 / 16.0);
    const auto fs_7_143_3 = std::sqrt(147.0 / 20449.0);
    const auto fs_7_286_10 = std::sqrt(245.0 / 40898.0);
    const auto fs_810_143_2 = std::sqrt(1312200.0 / 20449.0);
    const auto fs_810_143_33 = std::sqrt(1968300.0 / 1859.0);
    const auto fs_84_143_35 = std::sqrt(246960.0 / 20449.0);
    const auto fs_84_2431_5 = std::sqrt(35280.0 / 5909761.0);
    const auto fs_882_96577_165 = std::sqrt(128357460.0 / 9327116929.0);
    const auto fs_882_96577_39 = std::sqrt(2333772.0 / 717470533.0);
    const auto fs_90_11_7 = std::sqrt(56700.0 / 121.0);
    const auto fs_924_96577_286 = std::sqrt(18783072.0 / 717470533.0);
    const auto fs_945_11_3 = std::sqrt(2679075.0 / 121.0);
    const auto fs_945_11_5 = std::sqrt(4465125.0 / 121.0);
    const auto fs_945_143_35 = std::sqrt(31255875.0 / 20449.0);
    const auto fs_945_16_14 = std::sqrt(6251175.0 / 128.0);
    const auto fs_945_16_3 = std::sqrt(2679075.0 / 256.0);
    const auto fs_945_16_30 = std::sqrt(13395375.0 / 128.0);
    const auto fs_945_16_70 = std::sqrt(31255875.0 / 128.0);
    const auto fs_945_22_14 = std::sqrt(6251175.0 / 242.0);
    const auto fs_945_22_3 = std::sqrt(2679075.0 / 484.0);
    const auto fs_945_2717_30 = std::sqrt(26790750.0 / 7382089.0);
    const auto fs_945_286_7 = std::sqrt(6251175.0 / 81796.0);
    const auto fs_945_32_22 = std::sqrt(9823275.0 / 512.0);
    const auto fs_945_32_30 = std::sqrt(13395375.0 / 512.0);
    const auto fs_945_32_55 = std::sqrt(49116375.0 / 1024.0);
    const auto fs_945_32_7 = std::sqrt(6251175.0 / 1024.0);
    const auto fs_945_572_30 = std::sqrt(13395375.0 / 163592.0);
    const auto fs_945_572_462 = std::sqrt(18753525.0 / 14872.0);
    const auto fs_945_8_3 = std::sqrt(2679075.0 / 64.0);
    const auto fs_945_8_5 = std::sqrt(4465125.0 / 64.0);
    const auto fs_99_193154_10010 = std::sqrt(3773385.0 / 1434941066.0);
    const auto fs_99_193154_2730 = std::sqrt(1029105.0 / 1434941066.0);
    const auto fs_99_96577_20995 = std::sqrt(49005.0 / 2221271.0);
    const auto fs_99_96577_36465 = std::sqrt(1617165.0 / 42204149.0);
    const auto fs_99_96577_58786 = std::sqrt(137214.0 / 2221271.0);
    const auto fs_99_96577_6006 = std::sqrt(4528062.0 / 717470533.0);
    const auto fs_9_11_10 = std::sqrt(810.0 / 121.0);
    const auto fs_9_13_30 = std::sqrt(2430.0 / 169.0);
    const auto fs_9_2431_55 = std::sqrt(405.0 / 537251.0);
    const auto fs_9_2431_77 = std::sqrt(567.0 / 537251.0);
    const auto fs_9_374_10 = std::sqrt(405.0 / 69938.0);
    const auto fs_9_442_30 = std::sqrt(1215.0 / 97682.0);

    // NOTE: the angular components are formed in 84 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph2_0, ph4_0, ph6_0, ph8_0, ph10_0, ph12_0, ph12_p12, ab_2, pc_0 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
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
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;
        const auto r_12 = r_2 * r_10;

        const auto h2_0 = ph2_0[k];
        const auto h4_0 = ph4_0[k];
        const auto h6_0 = ph6_0[k];
        const auto h8_0 = ph8_0[k];
        const auto h10_0 = ph10_0[k];
        const auto h12_0 = ph12_0[k];
        const auto h12_p12 = ph12_p12[k];

        pc_0[k] = e_0 * f_10395_64 + e_1 * f_10395_16 * h2_0 - e_1 * f_10395_16 * r_2 + e_2 * f_4455_16 * h4_0 - e_2 * f_7425_8 * r_2 * h2_0 + e_2 * f_10395_16 * r_4 + e_3 * f_75_2 * h6_0 - e_3 * f_405_2 * r_2 * h4_0 + e_3 * f_825_2 * r_4 * h2_0 - e_3 * f_495_2 * r_6 + e_4 * f_105_52 * h8_0 - e_4 * f_15_1 * r_2 * h6_0 + e_4 * f_1215_26 * r_4 * h4_0 - e_4 * f_75_1 * r_6 * h2_0 + e_4 * f_165_4 * r_8 + e_5 * f_189_4199 * h10_0 - e_5 * f_105_247 * r_2 * h8_0 + e_5 * f_30_17 * r_4 * h6_0 - e_5 * f_54_13 * r_6 * h4_0 + e_5 * f_75_13 * r_8 * h2_0 - e_5 * f_3_1 * r_10 + e_6 * f_33_96577 * h12_0 - e_6 * fs_33_96577_1352078 * h12_p12 - e_6 * f_378_96577 * r_2 * h10_0 + e_6 * f_5_247 * r_4 * h8_0 - e_6 * f_20_323 * r_6 * h6_0 + e_6 * f_27_221 * r_8 * h4_0 - e_6 * f_2_13 * r_10 * h2_0 + e_6 * f_1_13 * r_12;
    }

    // NOTE: the angular components are formed in 84 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph2_p1, ph4_p1, ph6_p1, ph8_p1, ph10_p1, ph12_p1, ph12_p11, ab_2, pc_1 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];

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
        const auto h12_p1 = ph12_p1[k];
        const auto h12_p11 = ph12_p11[k];

        pc_1[k] = - e_1 * f_10395_32 * h2_p1 - e_2 * fs_1485_32_30 * h4_p1 + e_2 * f_7425_16 * r_2 * h2_p1 - e_3 * fs_75_4_7 * h6_p1 + e_3 * fs_135_4_30 * r_2 * h4_p1 - e_3 * f_825_4 * r_4 * h2_p1 - e_4 * fs_105_52_3 * h8_p1 + e_4 * fs_15_2_7 * r_2 * h6_p1 - e_4 * fs_405_52_30 * r_4 * h4_p1 + e_4 * f_75_2 * r_6 * h2_p1 - e_5 * fs_63_8398_165 * h10_p1 + e_5 * fs_105_247_3 * r_2 * h8_p1 - e_5 * fs_15_17_7 * r_4 * h6_p1 + e_5 * fs_9_13_30 * r_6 * h4_p1 - e_5 * f_75_26 * r_8 * h2_p1 - e_6 * fs_33_193154_26 * h12_p1 - e_6 * fs_33_96577_676039 * h12_p11 + e_6 * fs_63_96577_165 * r_2 * h10_p1 - e_6 * fs_5_247_3 * r_4 * h8_p1 + e_6 * fs_10_323_7 * r_6 * h6_p1 - e_6 * fs_9_442_30 * r_8 * h4_p1 + e_6 * f_1_13 * r_10 * h2_p1;
    }

    // NOTE: the angular components are formed in 84 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph2_p2, ph4_p2, ph6_p2, ph8_p2, ph10_p2, ph10_p10, ph12_p2, ph12_p10, ab_2, pc_2 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];

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
        const auto h12_p2 = ph12_p2[k];
        const auto h12_p10 = ph12_p10[k];

        pc_2[k] = e_1 * fs_945_32_22 * h2_p2 + e_2 * fs_405_32_330 * h4_p2 - e_2 * fs_675_16_22 * r_2 * h2_p2 + e_3 * fs_75_22_385 * h6_p2 - e_3 * fs_405_44_330 * r_2 * h4_p2 + e_3 * fs_75_4_22 * r_4 * h2_p2 + e_4 * fs_105_572_1155 * h8_p2 - e_4 * fs_15_11_385 * r_2 * h6_p2 + e_4 * fs_1215_572_330 * r_4 * h4_p2 - e_4 * fs_75_22_22 * r_6 * h2_p2 + e_5 * fs_567_8398_10 * h10_p2 - e_5 * fs_63_4199_12597 * h10_p10 - e_5 * fs_105_2717_1155 * r_2 * h8_p2 + e_5 * fs_30_187_385 * r_4 * h6_p2 - e_5 * fs_27_143_330 * r_6 * h4_p2 + e_5 * fs_75_286_22 * r_8 * h2_p2 + e_6 * fs_33_193154_182 * h12_p2 - e_6 * fs_33_96577_323323 * h12_p10 - e_6 * fs_567_96577_10 * r_2 * h10_p2 + e_6 * fs_126_96577_12597 * r_2 * h10_p10 + e_6 * fs_5_2717_1155 * r_4 * h8_p2 - e_6 * fs_20_3553_385 * r_6 * h6_p2 + e_6 * fs_27_4862_330 * r_8 * h4_p2 - e_6 * fs_1_143_22 * r_10 * h2_p2;
    }

    // NOTE: the angular components are formed in 84 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_2, pe_3, pe_4, pe_5, pe_6, ph4_p3, ph6_p3, ph8_p3, ph10_p3, ph10_p9, ph12_p3, ph12_p9, ab_2, pc_3 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

        const auto h4_p3 = ph4_p3[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h8_p3 = ph8_p3[k];
        const auto h10_p3 = ph10_p3[k];
        const auto h10_p9 = ph10_p9[k];
        const auto h12_p3 = ph12_p3[k];
        const auto h12_p9 = ph12_p9[k];

        pc_3[k] = - e_2 * fs_405_32_154 * h4_p3 - e_3 * fs_75_22_462 * h6_p3 + e_3 * fs_405_44_154 * r_2 * h4_p3 - e_4 * fs_105_52_21 * h8_p3 + e_4 * fs_15_11_462 * r_2 * h6_p3 - e_4 * fs_1215_572_154 * r_4 * h4_p3 - e_5 * fs_189_4199_78 * h10_p3 - e_5 * fs_189_8398_8398 * h10_p9 + e_5 * fs_105_247_21 * r_2 * h8_p3 - e_5 * fs_30_187_462 * r_4 * h6_p3 + e_5 * fs_27_143_154 * r_6 * h4_p3 - e_6 * fs_33_193154_910 * h12_p3 - e_6 * fs_33_96577_146965 * h12_p9 + e_6 * fs_378_96577_78 * r_2 * h10_p3 + e_6 * fs_189_96577_8398 * r_2 * h10_p9 - e_6 * fs_5_247_21 * r_4 * h8_p3 + e_6 * fs_20_3553_462 * r_6 * h6_p3 - e_6 * fs_27_4862_154 * r_8 * h4_p3;
    }

    // NOTE: the angular components are formed in 84 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_2, pe_3, pe_4, pe_5, pe_6, ph4_p4, ph6_p4, ph8_p4, ph8_p8, ph10_p4, ph10_p8, ph12_p4, ph12_p8, ab_2, pc_4 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

        const auto h4_p4 = ph4_p4[k];
        const auto h6_p4 = ph6_p4[k];
        const auto h8_p4 = ph8_p4[k];
        const auto h8_p8 = ph8_p8[k];
        const auto h10_p4 = ph10_p4[k];
        const auto h10_p8 = ph10_p8[k];
        const auto h12_p4 = ph12_p4[k];
        const auto h12_p8 = ph12_p8[k];

        pc_4[k] = e_2 * fs_135_16_77 * h4_p4 + e_3 * fs_75_22_385 * h6_p4 - e_3 * fs_135_22_77 * r_2 * h4_p4 + e_4 * fs_105_52_35 * h8_p4 - e_4 * fs_105_52_13 * h8_p8 - e_4 * fs_15_11_385 * r_2 * h6_p4 + e_4 * fs_405_286_77 * r_4 * h4_p4 + e_5 * fs_441_4199_39 * h10_p4 - e_5 * fs_567_4199_221 * h10_p8 - e_5 * fs_105_247_35 * r_2 * h8_p4 + e_5 * fs_105_247_13 * r_2 * h8_p8 + e_5 * fs_30_187_385 * r_4 * h6_p4 - e_5 * fs_18_143_77 * r_6 * h4_p4 + e_6 * fs_33_96577_910 * h12_p4 - e_6 * fs_33_96577_62985 * h12_p8 - e_6 * fs_882_96577_39 * r_2 * h10_p4 + e_6 * fs_1134_96577_221 * r_2 * h10_p8 + e_6 * fs_5_247_35 * r_4 * h8_p4 - e_6 * fs_5_247_13 * r_4 * h8_p8 - e_6 * fs_20_3553_385 * r_6 * h6_p4 + e_6 * fs_9_2431_77 * r_8 * h4_p4;
    }

    // NOTE: the angular components are formed in 84 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_3, pe_4, pe_5, pe_6, ph6_m6, ph6_p5, ph8_m6, ph8_p5, ph8_p7, ph10_m6, ph10_p5, ph10_p7, ph12_m6, ph12_p5, ph12_p7, ab_2, pc_5, pc_6 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

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

        pc_5[k] = - e_3 * fs_75_4_7 * h6_p5 - e_4 * fs_105_104_182 * h8_p5 - e_4 * fs_105_104_130 * h8_p7 + e_4 * fs_15_2_7 * r_2 * h6_p5 - e_5 * fs_1323_8398_39 * h10_p5 - e_5 * fs_126_4199_3315 * h10_p7 + e_5 * fs_105_494_182 * r_2 * h8_p5 + e_5 * fs_105_494_130 * r_2 * h8_p7 - e_5 * fs_15_17_7 * r_4 * h6_p5 - e_6 * fs_33_96577_3094 * h12_p5 - e_6 * fs_33_96577_25194 * h12_p7 + e_6 * fs_1323_96577_39 * r_2 * h10_p5 + e_6 * fs_252_96577_3315 * r_2 * h10_p7 - e_6 * fs_5_494_182 * r_4 * h8_p5 - e_6 * fs_5_494_130 * r_4 * h8_p7 + e_6 * fs_10_323_7 * r_6 * h6_p5;

        pc_6[k] = e_3 * f_75_2 * h6_m6 + e_4 * fs_105_52_91 * h8_m6 - e_4 * f_15_1 * r_2 * h6_m6 + e_5 * fs_378_4199_455 * h10_m6 - e_5 * fs_105_247_91 * r_2 * h8_m6 + e_5 * f_30_17 * r_4 * h6_m6 + e_6 * fs_66_96577_4641 * h12_m6 - e_6 * fs_756_96577_455 * r_2 * h10_m6 + e_6 * fs_5_247_91 * r_4 * h8_m6 - e_6 * f_20_323 * r_6 * h6_m6;
    }

    // NOTE: the angular components are formed in 84 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_3, pe_4, pe_5, pe_6, ph6_m5, ph8_m7, ph8_m5, ph10_m7, ph10_m5, ph12_m7, ph12_m5, ab_2, pc_7 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h6_m5 = ph6_m5[k];
        const auto h8_m7 = ph8_m7[k];
        const auto h8_m5 = ph8_m5[k];
        const auto h10_m7 = ph10_m7[k];
        const auto h10_m5 = ph10_m5[k];
        const auto h12_m7 = ph12_m7[k];
        const auto h12_m5 = ph12_m5[k];

        pc_7[k] = - e_3 * fs_75_4_7 * h6_m5 + e_4 * fs_105_104_130 * h8_m7 - e_4 * fs_105_104_182 * h8_m5 + e_4 * fs_15_2_7 * r_2 * h6_m5 + e_5 * fs_126_4199_3315 * h10_m7 - e_5 * fs_1323_8398_39 * h10_m5 - e_5 * fs_105_494_130 * r_2 * h8_m7 + e_5 * fs_105_494_182 * r_2 * h8_m5 - e_5 * fs_15_17_7 * r_4 * h6_m5 + e_6 * fs_33_96577_25194 * h12_m7 - e_6 * fs_33_96577_3094 * h12_m5 - e_6 * fs_252_96577_3315 * r_2 * h10_m7 + e_6 * fs_1323_96577_39 * r_2 * h10_m5 + e_6 * fs_5_494_130 * r_4 * h8_m7 - e_6 * fs_5_494_182 * r_4 * h8_m5 + e_6 * fs_10_323_7 * r_6 * h6_m5;
    }

    // NOTE: the angular components are formed in 84 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_2, pe_3, pe_4, pe_5, pe_6, ph4_m4, ph6_m4, ph8_m8, ph8_m4, ph10_m8, ph10_m4, ph12_m8, ph12_m4, ab_2, pc_8 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

        const auto h4_m4 = ph4_m4[k];
        const auto h6_m4 = ph6_m4[k];
        const auto h8_m8 = ph8_m8[k];
        const auto h8_m4 = ph8_m4[k];
        const auto h10_m8 = ph10_m8[k];
        const auto h10_m4 = ph10_m4[k];
        const auto h12_m8 = ph12_m8[k];
        const auto h12_m4 = ph12_m4[k];

        pc_8[k] = e_2 * fs_135_16_77 * h4_m4 + e_3 * fs_75_22_385 * h6_m4 - e_3 * fs_135_22_77 * r_2 * h4_m4 + e_4 * fs_105_52_13 * h8_m8 + e_4 * fs_105_52_35 * h8_m4 - e_4 * fs_15_11_385 * r_2 * h6_m4 + e_4 * fs_405_286_77 * r_4 * h4_m4 + e_5 * fs_567_4199_221 * h10_m8 + e_5 * fs_441_4199_39 * h10_m4 - e_5 * fs_105_247_13 * r_2 * h8_m8 - e_5 * fs_105_247_35 * r_2 * h8_m4 + e_5 * fs_30_187_385 * r_4 * h6_m4 - e_5 * fs_18_143_77 * r_6 * h4_m4 + e_6 * fs_33_96577_62985 * h12_m8 + e_6 * fs_33_96577_910 * h12_m4 - e_6 * fs_1134_96577_221 * r_2 * h10_m8 - e_6 * fs_882_96577_39 * r_2 * h10_m4 + e_6 * fs_5_247_13 * r_4 * h8_m8 + e_6 * fs_5_247_35 * r_4 * h8_m4 - e_6 * fs_20_3553_385 * r_6 * h6_m4 + e_6 * fs_9_2431_77 * r_8 * h4_m4;
    }

    // NOTE: the angular components are formed in 84 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_2, pe_3, pe_4, pe_5, pe_6, ph4_m3, ph6_m3, ph8_m3, ph10_m9, ph10_m3, ph12_m9, ph12_m3, ab_2, pc_9 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

        const auto h4_m3 = ph4_m3[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h8_m3 = ph8_m3[k];
        const auto h10_m9 = ph10_m9[k];
        const auto h10_m3 = ph10_m3[k];
        const auto h12_m9 = ph12_m9[k];
        const auto h12_m3 = ph12_m3[k];

        pc_9[k] = - e_2 * fs_405_32_154 * h4_m3 - e_3 * fs_75_22_462 * h6_m3 + e_3 * fs_405_44_154 * r_2 * h4_m3 - e_4 * fs_105_52_21 * h8_m3 + e_4 * fs_15_11_462 * r_2 * h6_m3 - e_4 * fs_1215_572_154 * r_4 * h4_m3 + e_5 * fs_189_8398_8398 * h10_m9 - e_5 * fs_189_4199_78 * h10_m3 + e_5 * fs_105_247_21 * r_2 * h8_m3 - e_5 * fs_30_187_462 * r_4 * h6_m3 + e_5 * fs_27_143_154 * r_6 * h4_m3 + e_6 * fs_33_96577_146965 * h12_m9 - e_6 * fs_33_193154_910 * h12_m3 - e_6 * fs_189_96577_8398 * r_2 * h10_m9 + e_6 * fs_378_96577_78 * r_2 * h10_m3 - e_6 * fs_5_247_21 * r_4 * h8_m3 + e_6 * fs_20_3553_462 * r_6 * h6_m3 - e_6 * fs_27_4862_154 * r_8 * h4_m3;
    }

    // NOTE: the angular components are formed in 84 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph2_m2, ph4_m2, ph6_m2, ph8_m2, ph10_m10, ph10_m2, ph12_m10, ph12_m2, ab_2, pc_10 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];

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
        const auto h12_m10 = ph12_m10[k];
        const auto h12_m2 = ph12_m2[k];

        pc_10[k] = e_1 * fs_945_32_22 * h2_m2 + e_2 * fs_405_32_330 * h4_m2 - e_2 * fs_675_16_22 * r_2 * h2_m2 + e_3 * fs_75_22_385 * h6_m2 - e_3 * fs_405_44_330 * r_2 * h4_m2 + e_3 * fs_75_4_22 * r_4 * h2_m2 + e_4 * fs_105_572_1155 * h8_m2 - e_4 * fs_15_11_385 * r_2 * h6_m2 + e_4 * fs_1215_572_330 * r_4 * h4_m2 - e_4 * fs_75_22_22 * r_6 * h2_m2 + e_5 * fs_63_4199_12597 * h10_m10 + e_5 * fs_567_8398_10 * h10_m2 - e_5 * fs_105_2717_1155 * r_2 * h8_m2 + e_5 * fs_30_187_385 * r_4 * h6_m2 - e_5 * fs_27_143_330 * r_6 * h4_m2 + e_5 * fs_75_286_22 * r_8 * h2_m2 + e_6 * fs_33_96577_323323 * h12_m10 + e_6 * fs_33_193154_182 * h12_m2 - e_6 * fs_126_96577_12597 * r_2 * h10_m10 - e_6 * fs_567_96577_10 * r_2 * h10_m2 + e_6 * fs_5_2717_1155 * r_4 * h8_m2 - e_6 * fs_20_3553_385 * r_6 * h6_m2 + e_6 * fs_27_4862_330 * r_8 * h4_m2 - e_6 * fs_1_143_22 * r_10 * h2_m2;
    }

    // NOTE: the angular components are formed in 84 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph2_m1, ph4_m1, ph6_m1, ph8_m1, ph10_m1, ph12_m12, ph12_m11, ph12_m1, ab_2, pc_11, pc_12 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;

        const auto h2_m1 = ph2_m1[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h8_m1 = ph8_m1[k];
        const auto h10_m1 = ph10_m1[k];
        const auto h12_m12 = ph12_m12[k];
        const auto h12_m11 = ph12_m11[k];
        const auto h12_m1 = ph12_m1[k];

        pc_11[k] = - e_1 * f_10395_32 * h2_m1 - e_2 * fs_1485_32_30 * h4_m1 + e_2 * f_7425_16 * r_2 * h2_m1 - e_3 * fs_75_4_7 * h6_m1 + e_3 * fs_135_4_30 * r_2 * h4_m1 - e_3 * f_825_4 * r_4 * h2_m1 - e_4 * fs_105_52_3 * h8_m1 + e_4 * fs_15_2_7 * r_2 * h6_m1 - e_4 * fs_405_52_30 * r_4 * h4_m1 + e_4 * f_75_2 * r_6 * h2_m1 - e_5 * fs_63_8398_165 * h10_m1 + e_5 * fs_105_247_3 * r_2 * h8_m1 - e_5 * fs_15_17_7 * r_4 * h6_m1 + e_5 * fs_9_13_30 * r_6 * h4_m1 - e_5 * f_75_26 * r_8 * h2_m1 + e_6 * fs_33_96577_676039 * h12_m11 - e_6 * fs_33_193154_26 * h12_m1 + e_6 * fs_63_96577_165 * r_2 * h10_m1 - e_6 * fs_5_247_3 * r_4 * h8_m1 + e_6 * fs_10_323_7 * r_6 * h6_m1 - e_6 * fs_9_442_30 * r_8 * h4_m1 + e_6 * f_1_13 * r_10 * h2_m1;

        pc_12[k] = e_6 * fs_33_96577_1352078 * h12_m12;
    }

    // NOTE: the angular components are formed in 84 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph2_0, ph4_0, ph6_0, ph8_0, ph10_0, ph10_p10, ph12_0, ph12_p10, ab_2, pc_13 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
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
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;
        const auto r_12 = r_2 * r_10;

        const auto h2_0 = ph2_0[k];
        const auto h4_0 = ph4_0[k];
        const auto h6_0 = ph6_0[k];
        const auto h8_0 = ph8_0[k];
        const auto h10_0 = ph10_0[k];
        const auto h10_p10 = ph10_p10[k];
        const auto h12_0 = ph12_0[k];
        const auto h12_p10 = ph12_p10[k];

        pc_13[k] = e_0 * f_10395_64 + e_1 * f_10395_32 * h2_0 - e_1 * f_10395_16 * r_2 - e_2 * f_1485_8 * h4_0 - e_2 * f_7425_16 * r_2 * h2_0 + e_2 * f_10395_16 * r_4 - e_3 * f_375_4 * h6_0 + e_3 * f_135_1 * r_2 * h4_0 + e_3 * f_825_4 * r_4 * h2_0 - e_3 * f_495_2 * r_6 - e_4 * f_525_52 * h8_0 + e_4 * f_75_2 * r_2 * h6_0 - e_4 * f_405_13 * r_4 * h4_0 - e_4 * f_75_2 * r_6 * h2_0 + e_4 * f_165_4 * r_8 - e_5 * f_3087_8398 * h10_0 + e_5 * fs_63_8398_92378 * h10_p10 + e_5 * f_525_247 * r_2 * h8_0 - e_5 * f_75_17 * r_4 * h6_0 + e_5 * f_36_13 * r_6 * h4_0 + e_5 * f_75_26 * r_8 * h2_0 - e_5 * f_3_1 * r_10 - e_6 * f_396_96577 * h12_0 - e_6 * fs_66_96577_176358 * h12_p10 + e_6 * f_3087_96577 * r_2 * h10_0 - e_6 * fs_63_96577_92378 * r_2 * h10_p10 - e_6 * f_25_247 * r_4 * h8_0 + e_6 * f_50_323 * r_6 * h6_0 - e_6 * f_18_221 * r_8 * h4_0 - e_6 * f_1_13 * r_10 * h2_0 + e_6 * f_1_13 * r_12;
    }

    // NOTE: the angular components are formed in 84 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph2_p1, ph4_p1, ph6_p1, ph8_p1, ph10_p1, ph10_p9, ph12_p1, ph12_p9, ab_2, pc_14 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];

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
        const auto h12_p1 = ph12_p1[k];
        const auto h12_p9 = ph12_p9[k];

        pc_14[k] = - e_1 * fs_2835_64_66 * h2_p1 - e_2 * fs_135_16_55 * h4_p1 + e_2 * fs_2025_32_66 * r_2 * h2_p1 + e_3 * fs_225_88_462 * h6_p1 + e_3 * fs_135_22_55 * r_2 * h4_p1 - e_3 * fs_225_8_66 * r_4 * h2_p1 + e_4 * fs_315_143_22 * h8_p1 - e_4 * fs_45_44_462 * r_2 * h6_p1 - e_4 * fs_405_286_55 * r_4 * h4_p1 + e_4 * fs_225_44_66 * r_6 * h2_p1 + e_5 * fs_2709_16796_10 * h10_p1 + e_5 * fs_63_8398_20995 * h10_p9 - e_5 * fs_1260_2717_22 * r_2 * h8_p1 + e_5 * fs_45_374_462 * r_4 * h6_p1 + e_5 * fs_18_143_55 * r_6 * h4_p1 - e_5 * fs_225_572_66 * r_8 * h2_p1 + e_6 * fs_33_96577_429 * h12_p1 - e_6 * fs_99_96577_58786 * h12_p9 - e_6 * fs_2709_193154_10 * r_2 * h10_p1 - e_6 * fs_63_96577_20995 * r_2 * h10_p9 + e_6 * fs_60_2717_22 * r_4 * h8_p1 - e_6 * fs_15_3553_462 * r_6 * h6_p1 - e_6 * fs_9_2431_55 * r_8 * h4_p1 + e_6 * fs_3_286_66 * r_10 * h2_p1;
    }

    // NOTE: the angular components are formed in 84 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph2_p2, ph4_p2, ph6_p2, ph8_p2, ph8_p8, ph10_p2, ph10_p8, ph12_p2, ph12_p8, ab_2, pc_15 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];

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
        const auto h12_p2 = ph12_p2[k];
        const auto h12_p8 = ph12_p8[k];

        pc_15[k] = e_1 * fs_945_32_55 * h2_p2 + e_2 * fs_135_4_33 * h4_p2 - e_2 * fs_675_16_55 * r_2 * h2_p2 - e_3 * fs_75_44_154 * h6_p2 - e_3 * fs_270_11_33 * r_2 * h4_p2 + e_3 * fs_75_4_55 * r_4 * h2_p2 - e_4 * fs_315_572_462 * h8_p2 + e_4 * fs_105_52_39 * h8_p8 + e_4 * fs_15_22_154 * r_2 * h6_p2 + e_4 * fs_810_143_33 * r_4 * h4_p2 - e_4 * fs_75_22_55 * r_6 * h2_p2 - e_5 * f_6993_8398 * h10_p2 - e_5 * fs_63_8398_663 * h10_p8 + e_5 * fs_315_2717_462 * r_2 * h8_p2 - e_5 * fs_105_247_39 * r_2 * h8_p8 - e_5 * fs_15_187_154 * r_4 * h6_p2 - e_5 * fs_72_143_33 * r_6 * h4_p2 + e_5 * fs_75_286_55 * r_8 * h2_p2 - e_6 * fs_66_96577_455 * h12_p2 - e_6 * fs_132_96577_20995 * h12_p8 + e_6 * f_6993_96577 * r_2 * h10_p2 + e_6 * fs_63_96577_663 * r_2 * h10_p8 - e_6 * fs_15_2717_462 * r_4 * h8_p2 + e_6 * fs_5_247_39 * r_4 * h8_p8 + e_6 * fs_10_3553_154 * r_6 * h6_p2 + e_6 * fs_36_2431_33 * r_8 * h4_p2 - e_6 * fs_1_143_55 * r_10 * h2_p2;
    }

    // NOTE: the angular components are formed in 84 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_2, pe_3, pe_4, pe_5, pe_6, ph4_p3, ph6_p3, ph8_p3, ph8_p7, ph10_p3, ph10_p7, ph12_p3, ph12_p7, ab_2, pc_16 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];

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
        const auto h12_p3 = ph12_p3[k];
        const auto h12_p7 = ph12_p7[k];

        pc_16[k] = - e_2 * fs_315_32_462 * h4_p3 - e_3 * fs_75_44_154 * h6_p3 + e_3 * fs_315_44_462 * r_2 * h4_p3 + e_4 * fs_105_26_7 * h8_p3 + e_4 * fs_105_52_39 * h8_p7 + e_4 * fs_15_22_154 * r_2 * h6_p3 - e_4 * fs_945_572_462 * r_4 * h4_p3 + e_5 * fs_1953_8398_26 * h10_p3 - e_5 * fs_441_8398_442 * h10_p7 - e_5 * fs_210_247_7 * r_2 * h8_p3 - e_5 * fs_105_247_39 * r_2 * h8_p7 - e_5 * fs_15_187_154 * r_4 * h6_p3 + e_5 * fs_21_143_462 * r_6 * h4_p3 + e_6 * fs_99_193154_2730 * h12_p3 - e_6 * fs_99_96577_20995 * h12_p7 - e_6 * fs_1953_96577_26 * r_2 * h10_p3 + e_6 * fs_441_96577_442 * r_2 * h10_p7 + e_6 * fs_10_247_7 * r_4 * h8_p3 + e_6 * fs_5_247_39 * r_4 * h8_p7 + e_6 * fs_10_3553_154 * r_6 * h6_p3 - e_6 * fs_21_4862_462 * r_8 * h4_p3;
    }

    // NOTE: the angular components are formed in 84 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_2, pe_3, pe_4, pe_5, pe_6, ph4_p4, ph6_m5, ph6_p4, ph6_p6, ph8_p4, ph8_p6, ph10_m5, ph10_p4, ph10_p6, ph12_m5, ph12_p4, ph12_p6, ab_2, pc_17, pc_18 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

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

        pc_17[k] = e_2 * fs_45_16_2310 * h4_p4 + e_3 * fs_225_88_462 * h6_p4 + e_3 * fs_75_4_7 * h6_p6 - e_3 * fs_45_22_2310 * r_2 * h4_p4 - e_4 * fs_105_104_42 * h8_p4 + e_4 * fs_105_52_13 * h8_p6 - e_4 * fs_45_44_462 * r_2 * h6_p4 - e_4 * fs_15_2_7 * r_2 * h6_p6 + e_4 * fs_135_286_2310 * r_4 * h4_p4 - e_5 * fs_2205_16796_130 * h10_p4 - e_5 * fs_63_323_65 * h10_p6 + e_5 * fs_105_494_42 * r_2 * h8_p4 - e_5 * fs_105_247_13 * r_2 * h8_p6 + e_5 * fs_45_374_462 * r_4 * h6_p4 + e_5 * fs_15_17_7 * r_4 * h6_p6 - e_5 * fs_6_143_2310 * r_6 * h4_p4 - e_6 * fs_264_96577_273 * h12_p4 - e_6 * fs_396_96577_663 * h12_p6 + e_6 * fs_2205_193154_130 * r_2 * h10_p4 + e_6 * fs_126_7429_65 * r_2 * h10_p6 - e_6 * fs_5_494_42 * r_4 * h8_p4 + e_6 * fs_5_247_13 * r_4 * h8_p6 - e_6 * fs_15_3553_462 * r_6 * h6_p4 - e_6 * fs_10_323_7 * r_6 * h6_p6 + e_6 * fs_3_2431_2310 * r_8 * h4_p4;

        pc_18[k] = - e_3 * f_375_4 * h6_m5 + e_4 * f_75_2 * r_2 * h6_m5 + e_5 * fs_63_442_273 * h10_m5 - e_5 * f_75_17 * r_4 * h6_m5 + e_6 * fs_462_96577_442 * h12_m5 - e_6 * fs_63_5083_273 * r_2 * h10_m5 + e_6 * f_50_323 * r_6 * h6_m5;
    }

    // NOTE: the angular components are formed in 84 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_2, pe_3, pe_4, pe_5, pe_6, ph4_m4, ph6_m6, ph6_m4, ph8_m6, ph8_m4, ph10_m6, ph10_m4, ph12_m6, ph12_m4, ab_2, pc_19 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];

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
        const auto h12_m6 = ph12_m6[k];
        const auto h12_m4 = ph12_m4[k];

        pc_19[k] = e_2 * fs_45_16_2310 * h4_m4 - e_3 * fs_75_4_7 * h6_m6 + e_3 * fs_225_88_462 * h6_m4 - e_3 * fs_45_22_2310 * r_2 * h4_m4 - e_4 * fs_105_52_13 * h8_m6 - e_4 * fs_105_104_42 * h8_m4 + e_4 * fs_15_2_7 * r_2 * h6_m6 - e_4 * fs_45_44_462 * r_2 * h6_m4 + e_4 * fs_135_286_2310 * r_4 * h4_m4 + e_5 * fs_63_323_65 * h10_m6 - e_5 * fs_2205_16796_130 * h10_m4 + e_5 * fs_105_247_13 * r_2 * h8_m6 + e_5 * fs_105_494_42 * r_2 * h8_m4 - e_5 * fs_15_17_7 * r_4 * h6_m6 + e_5 * fs_45_374_462 * r_4 * h6_m4 - e_5 * fs_6_143_2310 * r_6 * h4_m4 + e_6 * fs_396_96577_663 * h12_m6 - e_6 * fs_264_96577_273 * h12_m4 - e_6 * fs_126_7429_65 * r_2 * h10_m6 + e_6 * fs_2205_193154_130 * r_2 * h10_m4 - e_6 * fs_5_247_13 * r_4 * h8_m6 - e_6 * fs_5_494_42 * r_4 * h8_m4 + e_6 * fs_10_323_7 * r_6 * h6_m6 - e_6 * fs_15_3553_462 * r_6 * h6_m4 + e_6 * fs_3_2431_2310 * r_8 * h4_m4;
    }

    // NOTE: the angular components are formed in 84 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_2, pe_3, pe_4, pe_5, pe_6, ph4_m3, ph6_m3, ph8_m7, ph8_m3, ph10_m7, ph10_m3, ph12_m7, ph12_m3, ab_2, pc_20 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];

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
        const auto h12_m7 = ph12_m7[k];
        const auto h12_m3 = ph12_m3[k];

        pc_20[k] = - e_2 * fs_315_32_462 * h4_m3 - e_3 * fs_75_44_154 * h6_m3 + e_3 * fs_315_44_462 * r_2 * h4_m3 - e_4 * fs_105_52_39 * h8_m7 + e_4 * fs_105_26_7 * h8_m3 + e_4 * fs_15_22_154 * r_2 * h6_m3 - e_4 * fs_945_572_462 * r_4 * h4_m3 + e_5 * fs_441_8398_442 * h10_m7 + e_5 * fs_1953_8398_26 * h10_m3 + e_5 * fs_105_247_39 * r_2 * h8_m7 - e_5 * fs_210_247_7 * r_2 * h8_m3 - e_5 * fs_15_187_154 * r_4 * h6_m3 + e_5 * fs_21_143_462 * r_6 * h4_m3 + e_6 * fs_99_96577_20995 * h12_m7 + e_6 * fs_99_193154_2730 * h12_m3 - e_6 * fs_441_96577_442 * r_2 * h10_m7 - e_6 * fs_1953_96577_26 * r_2 * h10_m3 - e_6 * fs_5_247_39 * r_4 * h8_m7 + e_6 * fs_10_247_7 * r_4 * h8_m3 + e_6 * fs_10_3553_154 * r_6 * h6_m3 - e_6 * fs_21_4862_462 * r_8 * h4_m3;
    }

    // NOTE: the angular components are formed in 84 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph2_m2, ph4_m2, ph6_m2, ph8_m8, ph8_m2, ph10_m8, ph10_m2, ph12_m8, ph12_m2, ab_2, pc_21 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];

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
        const auto h12_m8 = ph12_m8[k];
        const auto h12_m2 = ph12_m2[k];

        pc_21[k] = e_1 * fs_945_32_55 * h2_m2 + e_2 * fs_135_4_33 * h4_m2 - e_2 * fs_675_16_55 * r_2 * h2_m2 - e_3 * fs_75_44_154 * h6_m2 - e_3 * fs_270_11_33 * r_2 * h4_m2 + e_3 * fs_75_4_55 * r_4 * h2_m2 - e_4 * fs_105_52_39 * h8_m8 - e_4 * fs_315_572_462 * h8_m2 + e_4 * fs_15_22_154 * r_2 * h6_m2 + e_4 * fs_810_143_33 * r_4 * h4_m2 - e_4 * fs_75_22_55 * r_6 * h2_m2 + e_5 * fs_63_8398_663 * h10_m8 - e_5 * f_6993_8398 * h10_m2 + e_5 * fs_105_247_39 * r_2 * h8_m8 + e_5 * fs_315_2717_462 * r_2 * h8_m2 - e_5 * fs_15_187_154 * r_4 * h6_m2 - e_5 * fs_72_143_33 * r_6 * h4_m2 + e_5 * fs_75_286_55 * r_8 * h2_m2 + e_6 * fs_132_96577_20995 * h12_m8 - e_6 * fs_66_96577_455 * h12_m2 - e_6 * fs_63_96577_663 * r_2 * h10_m8 + e_6 * f_6993_96577 * r_2 * h10_m2 - e_6 * fs_5_247_39 * r_4 * h8_m8 - e_6 * fs_15_2717_462 * r_4 * h8_m2 + e_6 * fs_10_3553_154 * r_6 * h6_m2 + e_6 * fs_36_2431_33 * r_8 * h4_m2 - e_6 * fs_1_143_55 * r_10 * h2_m2;
    }

    // NOTE: the angular components are formed in 84 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph2_m1, ph4_m1, ph6_m1, ph8_m1, ph10_m10, ph10_m9, ph10_m1, ph12_m10, ph12_m9, ph12_m1, ab_2, pc_22, pc_23 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;

        const auto h2_m1 = ph2_m1[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h8_m1 = ph8_m1[k];
        const auto h10_m10 = ph10_m10[k];
        const auto h10_m9 = ph10_m9[k];
        const auto h10_m1 = ph10_m1[k];
        const auto h12_m10 = ph12_m10[k];
        const auto h12_m9 = ph12_m9[k];
        const auto h12_m1 = ph12_m1[k];

        pc_22[k] = - e_1 * fs_2835_64_66 * h2_m1 - e_2 * fs_135_16_55 * h4_m1 + e_2 * fs_2025_32_66 * r_2 * h2_m1 + e_3 * fs_225_88_462 * h6_m1 + e_3 * fs_135_22_55 * r_2 * h4_m1 - e_3 * fs_225_8_66 * r_4 * h2_m1 + e_4 * fs_315_143_22 * h8_m1 - e_4 * fs_45_44_462 * r_2 * h6_m1 - e_4 * fs_405_286_55 * r_4 * h4_m1 + e_4 * fs_225_44_66 * r_6 * h2_m1 - e_5 * fs_63_8398_20995 * h10_m9 + e_5 * fs_2709_16796_10 * h10_m1 - e_5 * fs_1260_2717_22 * r_2 * h8_m1 + e_5 * fs_45_374_462 * r_4 * h6_m1 + e_5 * fs_18_143_55 * r_6 * h4_m1 - e_5 * fs_225_572_66 * r_8 * h2_m1 + e_6 * fs_99_96577_58786 * h12_m9 + e_6 * fs_33_96577_429 * h12_m1 + e_6 * fs_63_96577_20995 * r_2 * h10_m9 - e_6 * fs_2709_193154_10 * r_2 * h10_m1 + e_6 * fs_60_2717_22 * r_4 * h8_m1 - e_6 * fs_15_3553_462 * r_6 * h6_m1 - e_6 * fs_9_2431_55 * r_8 * h4_m1 + e_6 * fs_3_286_66 * r_10 * h2_m1;

        pc_23[k] = - e_5 * fs_63_8398_92378 * h10_m10 + e_6 * fs_66_96577_176358 * h12_m10 + e_6 * fs_63_96577_92378 * r_2 * h10_m10;
    }

    // NOTE: the angular components are formed in 84 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph2_m1, ph4_m1, ph6_m1, ph8_m1, ph10_m1, ph12_m11, ph12_m1, ab_2, pc_24 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;

        const auto h2_m1 = ph2_m1[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h8_m1 = ph8_m1[k];
        const auto h10_m1 = ph10_m1[k];
        const auto h12_m11 = ph12_m11[k];
        const auto h12_m1 = ph12_m1[k];

        pc_24[k] = e_1 * f_10395_32 * h2_m1 + e_2 * fs_1485_32_30 * h4_m1 - e_2 * f_7425_16 * r_2 * h2_m1 + e_3 * fs_75_4_7 * h6_m1 - e_3 * fs_135_4_30 * r_2 * h4_m1 + e_3 * f_825_4 * r_4 * h2_m1 + e_4 * fs_105_52_3 * h8_m1 - e_4 * fs_15_2_7 * r_2 * h6_m1 + e_4 * fs_405_52_30 * r_4 * h4_m1 - e_4 * f_75_2 * r_6 * h2_m1 + e_5 * fs_63_8398_165 * h10_m1 - e_5 * fs_105_247_3 * r_2 * h8_m1 + e_5 * fs_15_17_7 * r_4 * h6_m1 - e_5 * fs_9_13_30 * r_6 * h4_m1 + e_5 * f_75_26 * r_8 * h2_m1 + e_6 * fs_33_96577_676039 * h12_m11 + e_6 * fs_33_193154_26 * h12_m1 - e_6 * fs_63_96577_165 * r_2 * h10_m1 + e_6 * fs_5_247_3 * r_4 * h8_m1 - e_6 * fs_10_323_7 * r_6 * h6_m1 + e_6 * fs_9_442_30 * r_8 * h4_m1 - e_6 * f_1_13 * r_10 * h2_m1;
    }

    // NOTE: the angular components are formed in 84 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph2_0, ph4_0, ph6_0, ph8_0, ph8_p8, ph10_0, ph10_p8, ph12_0, ph12_p8, ab_2, pc_25 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
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
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;
        const auto r_12 = r_2 * r_10;

        const auto h2_0 = ph2_0[k];
        const auto h4_0 = ph4_0[k];
        const auto h6_0 = ph6_0[k];
        const auto h8_0 = ph8_0[k];
        const auto h8_p8 = ph8_p8[k];
        const auto h10_0 = ph10_0[k];
        const auto h10_p8 = ph10_p8[k];
        const auto h12_0 = ph12_0[k];
        const auto h12_p8 = ph12_p8[k];

        pc_25[k] = e_0 * f_10395_64 + e_1 * f_945_16 * h2_0 - e_1 * f_10395_16 * r_2 - e_2 * f_270_1 * h4_0 - e_2 * f_675_8 * r_2 * h2_0 + e_2 * f_10395_16 * r_4 + e_3 * f_150_11 * h6_0 + e_3 * f_2160_11 * r_2 * h4_0 + e_3 * f_75_2 * r_4 * h2_0 - e_3 * f_495_2 * r_6 + e_4 * f_9345_572 * h8_0 - e_4 * fs_315_572_715 * h8_p8 - e_4 * f_60_11 * r_2 * h6_0 - e_4 * f_6480_143 * r_4 * h4_0 - e_4 * f_75_11 * r_6 * h2_0 + e_4 * f_165_4 * r_8 + e_5 * f_5229_4199 * h10_0 + e_5 * fs_63_4199_12155 * h10_p8 - e_5 * f_9345_2717 * r_2 * h8_0 + e_5 * fs_315_2717_715 * r_2 * h8_p8 + e_5 * f_120_187 * r_4 * h6_0 + e_5 * f_576_143 * r_6 * h4_0 + e_5 * f_75_143 * r_8 * h2_0 - e_5 * f_3_1 * r_10 + e_6 * f_2178_96577 * h12_0 - e_6 * fs_66_96577_138567 * h12_p8 - e_6 * f_10458_96577 * r_2 * h10_0 - e_6 * fs_126_96577_12155 * r_2 * h10_p8 + e_6 * f_445_2717 * r_4 * h8_0 - e_6 * fs_15_2717_715 * r_4 * h8_p8 - e_6 * f_80_3553 * r_6 * h6_0 - e_6 * f_288_2431 * r_8 * h4_0 - e_6 * f_2_143 * r_10 * h2_0 + e_6 * f_1_13 * r_12;
    }

    // NOTE: the angular components are formed in 84 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph2_p1, ph4_p1, ph6_p1, ph8_p1, ph8_p7, ph10_p1, ph10_p7, ph12_p1, ph12_p7, ab_2, pc_26 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];

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
        const auto h12_p1 = ph12_p1[k];
        const auto h12_p7 = ph12_p7[k];

        pc_26[k] = - e_1 * fs_6615_64_10 * h2_p1 + e_2 * fs_945_16_3 * h4_p1 + e_2 * fs_4725_32_10 * r_2 * h2_p1 + e_3 * fs_375_88_70 * h6_p1 - e_3 * fs_945_22_3 * r_2 * h4_p1 - e_3 * fs_525_8_10 * r_4 * h2_p1 - e_4 * fs_945_572_30 * h8_p1 - e_4 * fs_105_572_858 * h8_p7 - e_4 * fs_75_44_70 * r_2 * h6_p1 + e_4 * fs_2835_286_3 * r_4 * h4_p1 + e_4 * fs_525_44_10 * r_6 * h2_p1 - e_5 * fs_2583_16796_66 * h10_p1 + e_5 * fs_189_8398_2431 * h10_p7 + e_5 * fs_945_2717_30 * r_2 * h8_p1 + e_5 * fs_105_2717_858 * r_2 * h8_p7 + e_5 * fs_75_374_70 * r_4 * h6_p1 - e_5 * fs_126_143_3 * r_6 * h4_p1 - e_5 * fs_525_572_10 * r_8 * h2_p1 - e_6 * fs_363_96577_65 * h12_p1 - e_6 * fs_33_96577_461890 * h12_p7 + e_6 * fs_2583_193154_66 * r_2 * h10_p1 - e_6 * fs_189_96577_2431 * r_2 * h10_p7 - e_6 * fs_45_2717_30 * r_4 * h8_p1 - e_6 * fs_5_2717_858 * r_4 * h8_p7 - e_6 * fs_25_3553_70 * r_6 * h6_p1 + e_6 * fs_63_2431_3 * r_8 * h4_p1 + e_6 * fs_7_286_10 * r_10 * h2_p1;
    }

    // NOTE: the angular components are formed in 84 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph2_p2, ph4_p2, ph6_p2, ph6_p6, ph8_p2, ph8_p6, ph10_p2, ph10_p6, ph12_p2, ph12_p6, ab_2, pc_27 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];

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
        const auto h12_p2 = ph12_p2[k];
        const auto h12_p6 = ph12_p6[k];

        pc_27[k] = e_1 * fs_2835_32_10 * h2_p2 + e_2 * fs_1035_32_6 * h4_p2 - e_2 * fs_2025_16_10 * r_2 * h2_p2 - e_3 * fs_225_11_7 * h6_p2 - e_3 * fs_75_22_385 * h6_p6 - e_3 * fs_1035_44_6 * r_2 * h4_p2 + e_3 * fs_225_4_10 * r_4 * h2_p2 + e_4 * fs_105_143_21 * h8_p2 + e_4 * fs_105_572_715 * h8_p6 + e_4 * fs_90_11_7 * r_2 * h6_p2 + e_4 * fs_15_11_385 * r_2 * h6_p6 + e_4 * fs_3105_572_6 * r_4 * h4_p2 - e_4 * fs_225_22_10 * r_6 * h2_p2 + e_5 * fs_2709_8398_22 * h10_p2 + e_5 * fs_63_4199_143 * h10_p6 - e_5 * fs_420_2717_21 * r_2 * h8_p2 - e_5 * fs_105_2717_715 * r_2 * h8_p6 - e_5 * fs_180_187_7 * r_4 * h6_p2 - e_5 * fs_30_187_385 * r_4 * h6_p6 - e_5 * fs_69_143_6 * r_6 * h4_p2 + e_5 * fs_225_286_10 * r_8 * h2_p2 + e_6 * fs_99_193154_10010 * h12_p2 - e_6 * fs_99_96577_36465 * h12_p6 - e_6 * fs_2709_96577_22 * r_2 * h10_p2 - e_6 * fs_126_96577_143 * r_2 * h10_p6 + e_6 * fs_20_2717_21 * r_4 * h8_p2 + e_6 * fs_5_2717_715 * r_4 * h8_p6 + e_6 * fs_120_3553_7 * r_6 * h6_p2 + e_6 * fs_20_3553_385 * r_6 * h6_p6 + e_6 * fs_69_4862_6 * r_8 * h4_p2 - e_6 * fs_3_143_10 * r_10 * h2_p2;
    }

    // NOTE: the angular components are formed in 84 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_2, pe_3, pe_4, pe_5, pe_6, ph4_p3, ph6_p3, ph6_p5, ph8_p3, ph8_p5, ph10_p3, ph10_p5, ph12_p3, ph12_p5, ab_2, pc_28 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];

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
        const auto h12_p3 = ph12_p3[k];
        const auto h12_p5 = ph12_p5[k];

        pc_28[k] = - e_2 * fs_225_16_210 * h4_p3 + e_3 * fs_375_88_70 * h6_p3 - e_3 * fs_225_88_462 * h6_p5 + e_3 * fs_225_22_210 * r_2 * h4_p3 + e_4 * fs_105_572_385 * h8_p3 + e_4 * fs_105_572_3003 * h8_p5 - e_4 * fs_75_44_70 * r_2 * h6_p3 + e_4 * fs_45_44_462 * r_2 * h6_p5 - e_4 * fs_675_286_210 * r_4 * h4_p3 - e_5 * fs_693_16796_1430 * h10_p3 - e_5 * fs_693_16796_286 * h10_p5 - e_5 * fs_105_2717_385 * r_2 * h8_p3 - e_5 * fs_105_2717_3003 * r_2 * h8_p5 + e_5 * fs_75_374_70 * r_4 * h6_p3 - e_5 * fs_45_374_462 * r_4 * h6_p5 + e_5 * fs_30_143_210 * r_6 * h4_p3 - e_6 * fs_99_96577_6006 * h12_p3 - e_6 * fs_66_96577_51051 * h12_p5 + e_6 * fs_693_193154_1430 * r_2 * h10_p3 + e_6 * fs_693_193154_286 * r_2 * h10_p5 + e_6 * fs_5_2717_385 * r_4 * h8_p3 + e_6 * fs_5_2717_3003 * r_4 * h8_p5 - e_6 * fs_25_3553_70 * r_6 * h6_p3 + e_6 * fs_15_3553_462 * r_6 * h6_p5 - e_6 * fs_15_2431_210 * r_8 * h4_p3;
    }

    // NOTE: the angular components are formed in 84 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_2, pe_3, pe_4, pe_5, pe_6, ph4_m4, ph6_m4, ph8_m4, ph10_m4, ph12_m4, ab_2, pc_29 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

        const auto h4_m4 = ph4_m4[k];
        const auto h6_m4 = ph6_m4[k];
        const auto h8_m4 = ph8_m4[k];
        const auto h10_m4 = ph10_m4[k];
        const auto h12_m4 = ph12_m4[k];

        pc_29[k] = e_2 * fs_945_8_5 * h4_m4 + e_3 * f_150_11 * h6_m4 - e_3 * fs_945_11_5 * r_2 * h4_m4 - e_4 * fs_2205_572_11 * h8_m4 - e_4 * f_60_11 * r_2 * h6_m4 + e_4 * fs_2835_143_5 * r_4 * h4_m4 + e_5 * fs_63_4199_15015 * h10_m4 + e_5 * fs_2205_2717_11 * r_2 * h8_m4 + e_5 * f_120_187 * r_4 * h6_m4 - e_5 * fs_252_143_5 * r_6 * h4_m4 + e_6 * fs_924_96577_286 * h12_m4 - e_6 * fs_126_96577_15015 * r_2 * h10_m4 - e_6 * fs_105_2717_11 * r_4 * h8_m4 - e_6 * f_80_3553 * r_6 * h6_m4 + e_6 * fs_126_2431_5 * r_8 * h4_m4;
    }

    // NOTE: the angular components are formed in 84 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_2, pe_3, pe_4, pe_5, pe_6, ph4_m3, ph6_m5, ph6_m3, ph8_m5, ph8_m3, ph10_m5, ph10_m3, ph12_m5, ph12_m3, ab_2, pc_30 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];

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
        const auto h12_m5 = ph12_m5[k];
        const auto h12_m3 = ph12_m3[k];

        pc_30[k] = - e_2 * fs_225_16_210 * h4_m3 + e_3 * fs_225_88_462 * h6_m5 + e_3 * fs_375_88_70 * h6_m3 + e_3 * fs_225_22_210 * r_2 * h4_m3 - e_4 * fs_105_572_3003 * h8_m5 + e_4 * fs_105_572_385 * h8_m3 - e_4 * fs_45_44_462 * r_2 * h6_m5 - e_4 * fs_75_44_70 * r_2 * h6_m3 - e_4 * fs_675_286_210 * r_4 * h4_m3 + e_5 * fs_693_16796_286 * h10_m5 - e_5 * fs_693_16796_1430 * h10_m3 + e_5 * fs_105_2717_3003 * r_2 * h8_m5 - e_5 * fs_105_2717_385 * r_2 * h8_m3 + e_5 * fs_45_374_462 * r_4 * h6_m5 + e_5 * fs_75_374_70 * r_4 * h6_m3 + e_5 * fs_30_143_210 * r_6 * h4_m3 + e_6 * fs_66_96577_51051 * h12_m5 - e_6 * fs_99_96577_6006 * h12_m3 - e_6 * fs_693_193154_286 * r_2 * h10_m5 + e_6 * fs_693_193154_1430 * r_2 * h10_m3 - e_6 * fs_5_2717_3003 * r_4 * h8_m5 + e_6 * fs_5_2717_385 * r_4 * h8_m3 - e_6 * fs_15_3553_462 * r_6 * h6_m5 - e_6 * fs_25_3553_70 * r_6 * h6_m3 - e_6 * fs_15_2431_210 * r_8 * h4_m3;
    }

    // NOTE: the angular components are formed in 84 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph2_m2, ph4_m2, ph6_m6, ph6_m2, ph8_m6, ph8_m2, ph10_m6, ph10_m2, ph12_m6, ph12_m2, ab_2, pc_31 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];

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
        const auto h12_m6 = ph12_m6[k];
        const auto h12_m2 = ph12_m2[k];

        pc_31[k] = e_1 * fs_2835_32_10 * h2_m2 + e_2 * fs_1035_32_6 * h4_m2 - e_2 * fs_2025_16_10 * r_2 * h2_m2 + e_3 * fs_75_22_385 * h6_m6 - e_3 * fs_225_11_7 * h6_m2 - e_3 * fs_1035_44_6 * r_2 * h4_m2 + e_3 * fs_225_4_10 * r_4 * h2_m2 - e_4 * fs_105_572_715 * h8_m6 + e_4 * fs_105_143_21 * h8_m2 - e_4 * fs_15_11_385 * r_2 * h6_m6 + e_4 * fs_90_11_7 * r_2 * h6_m2 + e_4 * fs_3105_572_6 * r_4 * h4_m2 - e_4 * fs_225_22_10 * r_6 * h2_m2 - e_5 * fs_63_4199_143 * h10_m6 + e_5 * fs_2709_8398_22 * h10_m2 + e_5 * fs_105_2717_715 * r_2 * h8_m6 - e_5 * fs_420_2717_21 * r_2 * h8_m2 + e_5 * fs_30_187_385 * r_4 * h6_m6 - e_5 * fs_180_187_7 * r_4 * h6_m2 - e_5 * fs_69_143_6 * r_6 * h4_m2 + e_5 * fs_225_286_10 * r_8 * h2_m2 + e_6 * fs_99_96577_36465 * h12_m6 + e_6 * fs_99_193154_10010 * h12_m2 + e_6 * fs_126_96577_143 * r_2 * h10_m6 - e_6 * fs_2709_96577_22 * r_2 * h10_m2 - e_6 * fs_5_2717_715 * r_4 * h8_m6 + e_6 * fs_20_2717_21 * r_4 * h8_m2 - e_6 * fs_20_3553_385 * r_6 * h6_m6 + e_6 * fs_120_3553_7 * r_6 * h6_m2 + e_6 * fs_69_4862_6 * r_8 * h4_m2 - e_6 * fs_3_143_10 * r_10 * h2_m2;
    }

    // NOTE: the angular components are formed in 84 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph2_m1, ph4_m1, ph6_m1, ph8_m7, ph8_m1, ph10_m7, ph10_m1, ph12_m7, ph12_m1, ab_2, pc_32 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];

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
        const auto h12_m7 = ph12_m7[k];
        const auto h12_m1 = ph12_m1[k];

        pc_32[k] = - e_1 * fs_6615_64_10 * h2_m1 + e_2 * fs_945_16_3 * h4_m1 + e_2 * fs_4725_32_10 * r_2 * h2_m1 + e_3 * fs_375_88_70 * h6_m1 - e_3 * fs_945_22_3 * r_2 * h4_m1 - e_3 * fs_525_8_10 * r_4 * h2_m1 + e_4 * fs_105_572_858 * h8_m7 - e_4 * fs_945_572_30 * h8_m1 - e_4 * fs_75_44_70 * r_2 * h6_m1 + e_4 * fs_2835_286_3 * r_4 * h4_m1 + e_4 * fs_525_44_10 * r_6 * h2_m1 - e_5 * fs_189_8398_2431 * h10_m7 - e_5 * fs_2583_16796_66 * h10_m1 - e_5 * fs_105_2717_858 * r_2 * h8_m7 + e_5 * fs_945_2717_30 * r_2 * h8_m1 + e_5 * fs_75_374_70 * r_4 * h6_m1 - e_5 * fs_126_143_3 * r_6 * h4_m1 - e_5 * fs_525_572_10 * r_8 * h2_m1 + e_6 * fs_33_96577_461890 * h12_m7 - e_6 * fs_363_96577_65 * h12_m1 + e_6 * fs_189_96577_2431 * r_2 * h10_m7 + e_6 * fs_2583_193154_66 * r_2 * h10_m1 + e_6 * fs_5_2717_858 * r_4 * h8_m7 - e_6 * fs_45_2717_30 * r_4 * h8_m1 - e_6 * fs_25_3553_70 * r_6 * h6_m1 + e_6 * fs_63_2431_3 * r_8 * h4_m1 + e_6 * fs_7_286_10 * r_10 * h2_m1;
    }

    // NOTE: the angular components are formed in 84 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph2_m1, ph4_m1, ph6_m1, ph8_m8, ph8_m1, ph10_m9, ph10_m8, ph10_m1, ph12_m9, ph12_m8, ph12_m1, ab_2, pc_33, pc_34 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;

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

        pc_33[k] = e_4 * fs_315_572_715 * h8_m8 - e_5 * fs_63_4199_12155 * h10_m8 - e_5 * fs_315_2717_715 * r_2 * h8_m8 + e_6 * fs_66_96577_138567 * h12_m8 + e_6 * fs_126_96577_12155 * r_2 * h10_m8 + e_6 * fs_15_2717_715 * r_4 * h8_m8;

        pc_34[k] = e_1 * fs_2835_64_66 * h2_m1 + e_2 * fs_135_16_55 * h4_m1 - e_2 * fs_2025_32_66 * r_2 * h2_m1 - e_3 * fs_225_88_462 * h6_m1 - e_3 * fs_135_22_55 * r_2 * h4_m1 + e_3 * fs_225_8_66 * r_4 * h2_m1 - e_4 * fs_315_143_22 * h8_m1 + e_4 * fs_45_44_462 * r_2 * h6_m1 + e_4 * fs_405_286_55 * r_4 * h4_m1 - e_4 * fs_225_44_66 * r_6 * h2_m1 - e_5 * fs_63_8398_20995 * h10_m9 - e_5 * fs_2709_16796_10 * h10_m1 + e_5 * fs_1260_2717_22 * r_2 * h8_m1 - e_5 * fs_45_374_462 * r_4 * h6_m1 - e_5 * fs_18_143_55 * r_6 * h4_m1 + e_5 * fs_225_572_66 * r_8 * h2_m1 + e_6 * fs_99_96577_58786 * h12_m9 - e_6 * fs_33_96577_429 * h12_m1 + e_6 * fs_63_96577_20995 * r_2 * h10_m9 + e_6 * fs_2709_193154_10 * r_2 * h10_m1 - e_6 * fs_60_2717_22 * r_4 * h8_m1 + e_6 * fs_15_3553_462 * r_6 * h6_m1 + e_6 * fs_9_2431_55 * r_8 * h4_m1 - e_6 * fs_3_286_66 * r_10 * h2_m1;
    }

    // NOTE: the angular components are formed in 84 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph2_m2, ph4_m2, ph6_m2, ph8_m2, ph10_m10, ph10_m2, ph12_m10, ph12_m2, ab_2, pc_35 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];

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
        const auto h12_m10 = ph12_m10[k];
        const auto h12_m2 = ph12_m2[k];

        pc_35[k] = - e_1 * fs_945_32_22 * h2_m2 - e_2 * fs_405_32_330 * h4_m2 + e_2 * fs_675_16_22 * r_2 * h2_m2 - e_3 * fs_75_22_385 * h6_m2 + e_3 * fs_405_44_330 * r_2 * h4_m2 - e_3 * fs_75_4_22 * r_4 * h2_m2 - e_4 * fs_105_572_1155 * h8_m2 + e_4 * fs_15_11_385 * r_2 * h6_m2 - e_4 * fs_1215_572_330 * r_4 * h4_m2 + e_4 * fs_75_22_22 * r_6 * h2_m2 + e_5 * fs_63_4199_12597 * h10_m10 - e_5 * fs_567_8398_10 * h10_m2 + e_5 * fs_105_2717_1155 * r_2 * h8_m2 - e_5 * fs_30_187_385 * r_4 * h6_m2 + e_5 * fs_27_143_330 * r_6 * h4_m2 - e_5 * fs_75_286_22 * r_8 * h2_m2 + e_6 * fs_33_96577_323323 * h12_m10 - e_6 * fs_33_193154_182 * h12_m2 - e_6 * fs_126_96577_12597 * r_2 * h10_m10 + e_6 * fs_567_96577_10 * r_2 * h10_m2 - e_6 * fs_5_2717_1155 * r_4 * h8_m2 + e_6 * fs_20_3553_385 * r_6 * h6_m2 - e_6 * fs_27_4862_330 * r_8 * h4_m2 + e_6 * fs_1_143_22 * r_10 * h2_m2;
    }

    // NOTE: the angular components are formed in 84 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph2_0, ph4_0, ph6_0, ph6_p6, ph8_0, ph8_p6, ph10_0, ph10_p6, ph12_0, ph12_p6, ab_2, pc_36 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
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
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;
        const auto r_12 = r_2 * r_10;

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

        pc_36[k] = e_0 * f_10395_64 - e_1 * f_4725_32 * h2_0 - e_1 * f_10395_16 * r_2 - e_2 * f_1215_8 * h4_0 + e_2 * f_3375_16 * r_2 * h2_0 + e_2 * f_10395_16 * r_4 + e_3 * f_3225_44 * h6_0 + e_3 * fs_75_22_462 * h6_p6 + e_3 * f_1215_11 * r_2 * h4_0 - e_3 * f_375_4 * r_4 * h2_0 - e_3 * f_495_2 * r_6 - e_4 * f_1995_572 * h8_0 - e_4 * fs_105_286_858 * h8_p6 - e_4 * f_645_22 * r_2 * h6_0 - e_4 * fs_15_11_462 * r_2 * h6_p6 - e_4 * f_3645_143 * r_4 * h4_0 + e_4 * f_375_22 * r_6 * h2_0 + e_4 * f_165_4 * r_8 - e_5 * f_945_442 * h10_0 + e_5 * fs_189_8398_4290 * h10_p6 + e_5 * f_105_143 * r_2 * h8_0 + e_5 * fs_210_2717_858 * r_2 * h8_p6 + e_5 * f_645_187 * r_4 * h6_0 + e_5 * fs_30_187_462 * r_4 * h6_p6 + e_5 * f_324_143 * r_6 * h4_0 - e_5 * f_375_286 * r_8 * h2_0 - e_5 * f_3_1 * r_10 - e_6 * f_7260_96577 * h12_0 - e_6 * fs_330_96577_4862 * h12_p6 + e_6 * f_945_5083 * r_2 * h10_0 - e_6 * fs_189_96577_4290 * r_2 * h10_p6 - e_6 * f_5_143 * r_4 * h8_0 - e_6 * fs_10_2717_858 * r_4 * h8_p6 - e_6 * f_430_3553 * r_6 * h6_0 - e_6 * fs_20_3553_462 * r_6 * h6_p6 - e_6 * f_162_2431 * r_8 * h4_0 + e_6 * f_5_143 * r_10 * h2_0 + e_6 * f_1_13 * r_12;
    }

    // NOTE: the angular components are formed in 84 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph2_p1, ph4_p1, ph6_p1, ph6_p5, ph8_p1, ph8_p5, ph10_p1, ph10_p5, ph12_p1, ph12_p5, ab_2, pc_37 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];

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
        const auto h12_p1 = ph12_p1[k];
        const auto h12_p5 = ph12_p5[k];

        pc_37[k] = - e_1 * fs_4725_32_3 * h2_p1 + e_2 * fs_1755_32_10 * h4_p1 + e_2 * fs_3375_16_3 * r_2 * h2_p1 - e_3 * fs_225_44_21 * h6_p1 + e_3 * fs_75_44_154 * h6_p5 - e_3 * fs_1755_44_10 * r_2 * h4_p1 - e_3 * fs_375_4_3 * r_4 * h2_p1 - e_4 * f_105_22 * h8_p1 - e_4 * fs_105_572_1001 * h8_p5 + e_4 * fs_45_22_21 * r_2 * h6_p1 - e_4 * fs_15_22_154 * r_2 * h6_p5 + e_4 * fs_405_44_10 * r_4 * h4_p1 + e_4 * fs_375_22_3 * r_6 * h2_p1 + e_5 * fs_1701_8398_55 * h10_p1 + e_5 * fs_315_8398_858 * h10_p5 + e_5 * f_210_209 * r_2 * h8_p1 + e_5 * fs_105_2717_1001 * r_2 * h8_p5 - e_5 * fs_45_187_21 * r_4 * h6_p1 + e_5 * fs_15_187_154 * r_4 * h6_p5 - e_5 * fs_9_11_10 * r_6 * h4_p1 - e_5 * fs_375_286_3 * r_8 * h2_p1 + e_6 * fs_1815_193154_78 * h12_p1 - e_6 * fs_165_96577_17017 * h12_p5 - e_6 * fs_1701_96577_55 * r_2 * h10_p1 - e_6 * fs_315_96577_858 * r_2 * h10_p5 - e_6 * f_10_209 * r_4 * h8_p1 - e_6 * fs_5_2717_1001 * r_4 * h8_p5 + e_6 * fs_30_3553_21 * r_6 * h6_p1 - e_6 * fs_10_3553_154 * r_6 * h6_p5 + e_6 * fs_9_374_10 * r_8 * h4_p1 + e_6 * fs_5_143_3 * r_10 * h2_p1;
    }

    // NOTE: the angular components are formed in 84 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph2_p2, ph4_p2, ph4_p4, ph6_p2, ph6_p4, ph8_p2, ph8_p4, ph10_p2, ph10_p4, ph12_p2, ph12_p4, ab_2, pc_38 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];

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
        const auto h12_p2 = ph12_p2[k];
        const auto h12_p4 = ph12_p4[k];

        pc_38[k] = e_1 * fs_945_16_30 * h2_p2 - e_2 * fs_135_4_2 * h4_p2 + e_2 * fs_945_16_14 * h4_p4 - e_2 * fs_675_8_30 * r_2 * h2_p2 - e_3 * fs_225_44_21 * h6_p2 - e_3 * fs_375_88_70 * h6_p4 + e_3 * fs_270_11_2 * r_2 * h4_p2 - e_3 * fs_945_22_14 * r_2 * h4_p4 + e_3 * fs_75_2_30 * r_4 * h2_p2 + e_4 * fs_1995_572_7 * h8_p2 + e_4 * fs_105_1144_770 * h8_p4 + e_4 * fs_45_22_21 * r_2 * h6_p2 + e_4 * fs_75_44_70 * r_2 * h6_p4 - e_4 * fs_810_143_2 * r_4 * h4_p2 + e_4 * fs_2835_286_14 * r_4 * h4_p4 - e_4 * fs_75_11_30 * r_6 * h2_p2 - e_5 * fs_63_442_66 * h10_p2 + e_5 * fs_189_16796_858 * h10_p4 - e_5 * fs_105_143_7 * r_2 * h8_p2 - e_5 * fs_105_5434_770 * r_2 * h8_p4 - e_5 * fs_45_187_21 * r_4 * h6_p2 - e_5 * fs_75_374_70 * r_4 * h6_p4 + e_5 * fs_72_143_2 * r_6 * h4_p2 - e_5 * fs_126_143_14 * r_6 * h4_p4 + e_5 * fs_75_143_30 * r_8 * h2_p2 - e_6 * fs_66_96577_30030 * h12_p2 - e_6 * fs_264_96577_5005 * h12_p4 + e_6 * fs_63_5083_66 * r_2 * h10_p2 - e_6 * fs_189_193154_858 * r_2 * h10_p4 + e_6 * fs_5_143_7 * r_4 * h8_p2 + e_6 * fs_5_5434_770 * r_4 * h8_p4 + e_6 * fs_30_3553_21 * r_6 * h6_p2 + e_6 * fs_25_3553_70 * r_6 * h6_p4 - e_6 * fs_36_2431_2 * r_8 * h4_p2 + e_6 * fs_63_2431_14 * r_8 * h4_p4 - e_6 * fs_2_143_30 * r_10 * h2_p2;
    }

    // NOTE: the angular components are formed in 84 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_2, pe_3, pe_4, pe_5, pe_6, ph4_m3, ph6_m3, ph8_m3, ph10_m3, ph12_m3, ab_2, pc_39 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

        const auto h4_m3 = ph4_m3[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h8_m3 = ph8_m3[k];
        const auto h10_m3 = ph10_m3[k];
        const auto h12_m3 = ph12_m3[k];

        pc_39[k] = - e_2 * fs_945_8_3 * h4_m3 + e_3 * f_3225_44 * h6_m3 + e_3 * fs_945_11_3 * r_2 * h4_m3 - e_4 * fs_735_286_22 * h8_m3 - e_4 * f_645_22 * r_2 * h6_m3 - e_4 * fs_2835_143_3 * r_4 * h4_m3 + e_5 * fs_189_8398_1001 * h10_m3 + e_5 * fs_1470_2717_22 * r_2 * h8_m3 + e_5 * f_645_187 * r_4 * h6_m3 + e_5 * fs_252_143_3 * r_6 * h4_m3 + e_6 * fs_462_96577_2145 * h12_m3 - e_6 * fs_189_96577_1001 * r_2 * h10_m3 - e_6 * fs_70_2717_22 * r_4 * h8_m3 - e_6 * f_430_3553 * r_6 * h6_m3 - e_6 * fs_126_2431_3 * r_8 * h4_m3;
    }

    // NOTE: the angular components are formed in 84 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph2_m2, ph4_m4, ph4_m2, ph6_m4, ph6_m2, ph8_m4, ph8_m2, ph10_m4, ph10_m2, ph12_m4, ph12_m2, ab_2, pc_40 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];

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
        const auto h12_m4 = ph12_m4[k];
        const auto h12_m2 = ph12_m2[k];

        pc_40[k] = e_1 * fs_945_16_30 * h2_m2 - e_2 * fs_945_16_14 * h4_m4 - e_2 * fs_135_4_2 * h4_m2 - e_2 * fs_675_8_30 * r_2 * h2_m2 + e_3 * fs_375_88_70 * h6_m4 - e_3 * fs_225_44_21 * h6_m2 + e_3 * fs_945_22_14 * r_2 * h4_m4 + e_3 * fs_270_11_2 * r_2 * h4_m2 + e_3 * fs_75_2_30 * r_4 * h2_m2 - e_4 * fs_105_1144_770 * h8_m4 + e_4 * fs_1995_572_7 * h8_m2 - e_4 * fs_75_44_70 * r_2 * h6_m4 + e_4 * fs_45_22_21 * r_2 * h6_m2 - e_4 * fs_2835_286_14 * r_4 * h4_m4 - e_4 * fs_810_143_2 * r_4 * h4_m2 - e_4 * fs_75_11_30 * r_6 * h2_m2 - e_5 * fs_189_16796_858 * h10_m4 - e_5 * fs_63_442_66 * h10_m2 + e_5 * fs_105_5434_770 * r_2 * h8_m4 - e_5 * fs_105_143_7 * r_2 * h8_m2 + e_5 * fs_75_374_70 * r_4 * h6_m4 - e_5 * fs_45_187_21 * r_4 * h6_m2 + e_5 * fs_126_143_14 * r_6 * h4_m4 + e_5 * fs_72_143_2 * r_6 * h4_m2 + e_5 * fs_75_143_30 * r_8 * h2_m2 + e_6 * fs_264_96577_5005 * h12_m4 - e_6 * fs_66_96577_30030 * h12_m2 + e_6 * fs_189_193154_858 * r_2 * h10_m4 + e_6 * fs_63_5083_66 * r_2 * h10_m2 - e_6 * fs_5_5434_770 * r_4 * h8_m4 + e_6 * fs_5_143_7 * r_4 * h8_m2 - e_6 * fs_25_3553_70 * r_6 * h6_m4 + e_6 * fs_30_3553_21 * r_6 * h6_m2 - e_6 * fs_63_2431_14 * r_8 * h4_m4 - e_6 * fs_36_2431_2 * r_8 * h4_m2 - e_6 * fs_2_143_30 * r_10 * h2_m2;
    }

    // NOTE: the angular components are formed in 84 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph2_m1, ph4_m1, ph6_m5, ph6_m1, ph8_m5, ph8_m1, ph10_m5, ph10_m1, ph12_m5, ph12_m1, ab_2, pc_41 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];

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
        const auto h12_m5 = ph12_m5[k];
        const auto h12_m1 = ph12_m1[k];

        pc_41[k] = - e_1 * fs_4725_32_3 * h2_m1 + e_2 * fs_1755_32_10 * h4_m1 + e_2 * fs_3375_16_3 * r_2 * h2_m1 - e_3 * fs_75_44_154 * h6_m5 - e_3 * fs_225_44_21 * h6_m1 - e_3 * fs_1755_44_10 * r_2 * h4_m1 - e_3 * fs_375_4_3 * r_4 * h2_m1 + e_4 * fs_105_572_1001 * h8_m5 - e_4 * f_105_22 * h8_m1 + e_4 * fs_15_22_154 * r_2 * h6_m5 + e_4 * fs_45_22_21 * r_2 * h6_m1 + e_4 * fs_405_44_10 * r_4 * h4_m1 + e_4 * fs_375_22_3 * r_6 * h2_m1 - e_5 * fs_315_8398_858 * h10_m5 + e_5 * fs_1701_8398_55 * h10_m1 - e_5 * fs_105_2717_1001 * r_2 * h8_m5 + e_5 * f_210_209 * r_2 * h8_m1 - e_5 * fs_15_187_154 * r_4 * h6_m5 - e_5 * fs_45_187_21 * r_4 * h6_m1 - e_5 * fs_9_11_10 * r_6 * h4_m1 - e_5 * fs_375_286_3 * r_8 * h2_m1 + e_6 * fs_165_96577_17017 * h12_m5 + e_6 * fs_1815_193154_78 * h12_m1 + e_6 * fs_315_96577_858 * r_2 * h10_m5 - e_6 * fs_1701_96577_55 * r_2 * h10_m1 + e_6 * fs_5_2717_1001 * r_4 * h8_m5 - e_6 * f_10_209 * r_4 * h8_m1 + e_6 * fs_10_3553_154 * r_6 * h6_m5 + e_6 * fs_30_3553_21 * r_6 * h6_m1 + e_6 * fs_9_374_10 * r_8 * h4_m1 + e_6 * fs_5_143_3 * r_10 * h2_m1;
    }

    // NOTE: the angular components are formed in 84 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_3, pe_4, pe_5, pe_6, ph6_m6, ph8_m6, ph10_m6, ph12_m6, ab_2, pc_42 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h6_m6 = ph6_m6[k];
        const auto h8_m6 = ph8_m6[k];
        const auto h10_m6 = ph10_m6[k];
        const auto h12_m6 = ph12_m6[k];

        pc_42[k] = - e_3 * fs_75_22_462 * h6_m6 + e_4 * fs_105_286_858 * h8_m6 + e_4 * fs_15_11_462 * r_2 * h6_m6 - e_5 * fs_189_8398_4290 * h10_m6 - e_5 * fs_210_2717_858 * r_2 * h8_m6 - e_5 * fs_30_187_462 * r_4 * h6_m6 + e_6 * fs_330_96577_4862 * h12_m6 + e_6 * fs_189_96577_4290 * r_2 * h10_m6 + e_6 * fs_10_2717_858 * r_4 * h8_m6 + e_6 * fs_20_3553_462 * r_6 * h6_m6;
    }

    // NOTE: the angular components are formed in 84 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph2_m1, ph4_m1, ph6_m1, ph8_m7, ph8_m1, ph10_m7, ph10_m1, ph12_m7, ph12_m1, ab_2, pc_43 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];

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
        const auto h12_m7 = ph12_m7[k];
        const auto h12_m1 = ph12_m1[k];

        pc_43[k] = e_1 * fs_6615_64_10 * h2_m1 - e_2 * fs_945_16_3 * h4_m1 - e_2 * fs_4725_32_10 * r_2 * h2_m1 - e_3 * fs_375_88_70 * h6_m1 + e_3 * fs_945_22_3 * r_2 * h4_m1 + e_3 * fs_525_8_10 * r_4 * h2_m1 + e_4 * fs_105_572_858 * h8_m7 + e_4 * fs_945_572_30 * h8_m1 + e_4 * fs_75_44_70 * r_2 * h6_m1 - e_4 * fs_2835_286_3 * r_4 * h4_m1 - e_4 * fs_525_44_10 * r_6 * h2_m1 - e_5 * fs_189_8398_2431 * h10_m7 + e_5 * fs_2583_16796_66 * h10_m1 - e_5 * fs_105_2717_858 * r_2 * h8_m7 - e_5 * fs_945_2717_30 * r_2 * h8_m1 - e_5 * fs_75_374_70 * r_4 * h6_m1 + e_5 * fs_126_143_3 * r_6 * h4_m1 + e_5 * fs_525_572_10 * r_8 * h2_m1 + e_6 * fs_33_96577_461890 * h12_m7 + e_6 * fs_363_96577_65 * h12_m1 + e_6 * fs_189_96577_2431 * r_2 * h10_m7 - e_6 * fs_2583_193154_66 * r_2 * h10_m1 + e_6 * fs_5_2717_858 * r_4 * h8_m7 + e_6 * fs_45_2717_30 * r_4 * h8_m1 + e_6 * fs_25_3553_70 * r_6 * h6_m1 - e_6 * fs_63_2431_3 * r_8 * h4_m1 - e_6 * fs_7_286_10 * r_10 * h2_m1;
    }

    // NOTE: the angular components are formed in 84 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph2_m2, ph4_m2, ph6_m2, ph8_m8, ph8_m2, ph10_m8, ph10_m2, ph12_m8, ph12_m2, ab_2, pc_44 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];

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
        const auto h12_m8 = ph12_m8[k];
        const auto h12_m2 = ph12_m2[k];

        pc_44[k] = - e_1 * fs_945_32_55 * h2_m2 - e_2 * fs_135_4_33 * h4_m2 + e_2 * fs_675_16_55 * r_2 * h2_m2 + e_3 * fs_75_44_154 * h6_m2 + e_3 * fs_270_11_33 * r_2 * h4_m2 - e_3 * fs_75_4_55 * r_4 * h2_m2 - e_4 * fs_105_52_39 * h8_m8 + e_4 * fs_315_572_462 * h8_m2 - e_4 * fs_15_22_154 * r_2 * h6_m2 - e_4 * fs_810_143_33 * r_4 * h4_m2 + e_4 * fs_75_22_55 * r_6 * h2_m2 + e_5 * fs_63_8398_663 * h10_m8 + e_5 * f_6993_8398 * h10_m2 + e_5 * fs_105_247_39 * r_2 * h8_m8 - e_5 * fs_315_2717_462 * r_2 * h8_m2 + e_5 * fs_15_187_154 * r_4 * h6_m2 + e_5 * fs_72_143_33 * r_6 * h4_m2 - e_5 * fs_75_286_55 * r_8 * h2_m2 + e_6 * fs_132_96577_20995 * h12_m8 + e_6 * fs_66_96577_455 * h12_m2 - e_6 * fs_63_96577_663 * r_2 * h10_m8 - e_6 * f_6993_96577 * r_2 * h10_m2 - e_6 * fs_5_247_39 * r_4 * h8_m8 + e_6 * fs_15_2717_462 * r_4 * h8_m2 - e_6 * fs_10_3553_154 * r_6 * h6_m2 - e_6 * fs_36_2431_33 * r_8 * h4_m2 + e_6 * fs_1_143_55 * r_10 * h2_m2;
    }

    // NOTE: the angular components are formed in 84 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_2, pe_3, pe_4, pe_5, pe_6, ph4_m3, ph6_m3, ph8_m3, ph10_m9, ph10_m3, ph12_m9, ph12_m3, ab_2, pc_45 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

        const auto h4_m3 = ph4_m3[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h8_m3 = ph8_m3[k];
        const auto h10_m9 = ph10_m9[k];
        const auto h10_m3 = ph10_m3[k];
        const auto h12_m9 = ph12_m9[k];
        const auto h12_m3 = ph12_m3[k];

        pc_45[k] = e_2 * fs_405_32_154 * h4_m3 + e_3 * fs_75_22_462 * h6_m3 - e_3 * fs_405_44_154 * r_2 * h4_m3 + e_4 * fs_105_52_21 * h8_m3 - e_4 * fs_15_11_462 * r_2 * h6_m3 + e_4 * fs_1215_572_154 * r_4 * h4_m3 + e_5 * fs_189_8398_8398 * h10_m9 + e_5 * fs_189_4199_78 * h10_m3 - e_5 * fs_105_247_21 * r_2 * h8_m3 + e_5 * fs_30_187_462 * r_4 * h6_m3 - e_5 * fs_27_143_154 * r_6 * h4_m3 + e_6 * fs_33_96577_146965 * h12_m9 + e_6 * fs_33_193154_910 * h12_m3 - e_6 * fs_189_96577_8398 * r_2 * h10_m9 - e_6 * fs_378_96577_78 * r_2 * h10_m3 + e_6 * fs_5_247_21 * r_4 * h8_m3 - e_6 * fs_20_3553_462 * r_6 * h6_m3 + e_6 * fs_27_4862_154 * r_8 * h4_m3;
    }

    // NOTE: the angular components are formed in 84 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph2_0, ph4_0, ph4_p4, ph6_0, ph6_p4, ph8_0, ph8_p4, ph10_0, ph10_p4, ph12_0, ph12_p4, ab_2, pc_46 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
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
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;
        const auto r_12 = r_2 * r_10;

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

        pc_46[k] = e_0 * f_10395_64 - e_1 * f_4725_16 * h2_0 - e_1 * f_10395_16 * r_2 + e_2 * f_495_16 * h4_0 - e_2 * fs_315_8_35 * h4_p4 + e_2 * f_3375_8 * r_2 * h2_0 + e_2 * f_10395_16 * r_4 + e_3 * f_75_2 * h6_0 + e_3 * fs_225_11_7 * h6_p4 - e_3 * f_45_2 * r_2 * h4_0 + e_3 * fs_315_11_35 * r_2 * h4_p4 - e_3 * f_375_2 * r_4 * h2_0 - e_3 * f_495_2 * r_6 - e_4 * f_7455_572 * h8_0 - e_4 * fs_315_286_77 * h8_p4 - e_4 * f_15_1 * r_2 * h6_0 - e_4 * fs_90_11_7 * r_2 * h6_p4 + e_4 * f_135_26 * r_4 * h4_0 - e_4 * fs_945_143_35 * r_4 * h4_p4 + e_4 * f_375_11 * r_6 * h2_0 + e_4 * f_165_4 * r_8 + e_5 * f_6615_4199 * h10_0 + e_5 * fs_126_4199_2145 * h10_p4 + e_5 * f_7455_2717 * r_2 * h8_0 + e_5 * fs_630_2717_77 * r_2 * h8_p4 + e_5 * f_30_17 * r_4 * h6_0 + e_5 * fs_180_187_7 * r_4 * h6_p4 - e_5 * f_6_13 * r_6 * h4_0 + e_5 * fs_84_143_35 * r_6 * h4_p4 - e_5 * f_375_143 * r_8 * h2_0 - e_5 * f_3_1 * r_10 + e_6 * f_16335_96577 * h12_0 - e_6 * fs_495_96577_2002 * h12_p4 - e_6 * f_13230_96577 * r_2 * h10_0 - e_6 * fs_252_96577_2145 * r_2 * h10_p4 - e_6 * f_355_2717 * r_4 * h8_0 - e_6 * fs_30_2717_77 * r_4 * h8_p4 - e_6 * f_20_323 * r_6 * h6_0 - e_6 * fs_120_3553_7 * r_6 * h6_p4 + e_6 * f_3_221 * r_8 * h4_0 - e_6 * fs_42_2431_35 * r_8 * h4_p4 + e_6 * f_10_143 * r_10 * h2_0 + e_6 * f_1_13 * r_12;
    }

    // NOTE: the angular components are formed in 84 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph2_p1, ph4_p1, ph4_p3, ph6_p1, ph6_p3, ph8_p1, ph8_p3, ph10_p1, ph10_p3, ph12_p1, ph12_p3, ab_2, pc_47 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];

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
        const auto h12_p1 = ph12_p1[k];
        const auto h12_p3 = ph12_p3[k];

        pc_47[k] = - e_1 * fs_945_32_30 * h2_p1 + e_2 * f_2385_16 * h4_p1 - e_2 * fs_315_16_7 * h4_p3 + e_2 * fs_675_16_30 * r_2 * h2_p1 - e_3 * fs_75_22_210 * h6_p1 + e_3 * fs_225_44_21 * h6_p3 - e_3 * f_2385_22 * r_2 * h4_p1 + e_3 * fs_315_22_7 * r_2 * h4_p3 - e_3 * fs_75_4_30 * r_4 * h2_p1 + e_4 * fs_2835_1144_10 * h8_p1 - e_4 * fs_315_1144_462 * h8_p3 + e_4 * fs_15_11_210 * r_2 * h6_p1 - e_4 * fs_45_22_21 * r_2 * h6_p3 + e_4 * f_7155_286 * r_4 * h4_p1 - e_4 * fs_945_286_7 * r_4 * h4_p3 + e_4 * fs_75_22_30 * r_6 * h2_p1 - e_5 * fs_378_4199_22 * h10_p1 + e_5 * fs_441_8398_429 * h10_p3 - e_5 * fs_2835_5434_10 * r_2 * h8_p1 + e_5 * fs_315_5434_462 * r_2 * h8_p3 - e_5 * fs_30_187_210 * r_4 * h6_p1 + e_5 * fs_45_187_21 * r_4 * h6_p3 - e_5 * f_318_143 * r_6 * h4_p1 + e_5 * fs_42_143_7 * r_6 * h4_p3 - e_5 * fs_75_286_30 * r_8 * h2_p1 - e_6 * fs_1089_96577_195 * h12_p1 - e_6 * fs_297_96577_5005 * h12_p3 + e_6 * fs_756_96577_22 * r_2 * h10_p1 - e_6 * fs_441_96577_429 * r_2 * h10_p3 + e_6 * fs_135_5434_10 * r_4 * h8_p1 - e_6 * fs_15_5434_462 * r_4 * h8_p3 + e_6 * fs_20_3553_210 * r_6 * h6_p1 - e_6 * fs_30_3553_21 * r_6 * h6_p3 + e_6 * f_159_2431 * r_8 * h4_p1 - e_6 * fs_21_2431_7 * r_8 * h4_p3 + e_6 * fs_1_143_30 * r_10 * h2_p1;
    }

    // NOTE: the angular components are formed in 84 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph2_m2, ph4_m2, ph6_m2, ph8_m2, ph10_m2, ph12_m2, ab_2, pc_48 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;

        const auto h2_m2 = ph2_m2[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h8_m2 = ph8_m2[k];
        const auto h10_m2 = ph10_m2[k];
        const auto h12_m2 = ph12_m2[k];

        pc_48[k] = e_1 * fs_945_16_70 * h2_m2 - e_2 * fs_495_16_42 * h4_m2 - e_2 * fs_675_8_70 * r_2 * h2_m2 + e_3 * f_75_2 * h6_m2 + e_3 * fs_45_2_42 * r_2 * h4_m2 + e_3 * fs_75_2_70 * r_4 * h2_m2 - e_4 * fs_735_572_3 * h8_m2 - e_4 * f_15_1 * r_2 * h6_m2 - e_4 * fs_135_26_42 * r_4 * h4_m2 - e_4 * fs_75_11_70 * r_6 * h2_m2 - e_5 * fs_189_4199_154 * h10_m2 + e_5 * fs_735_2717_3 * r_2 * h8_m2 + e_5 * f_30_17 * r_4 * h6_m2 + e_5 * fs_6_13_42 * r_6 * h4_m2 + e_5 * fs_75_143_70 * r_8 * h2_m2 + e_6 * fs_693_96577_1430 * h12_m2 + e_6 * fs_378_96577_154 * r_2 * h10_m2 - e_6 * fs_35_2717_3 * r_4 * h8_m2 - e_6 * f_20_323 * r_6 * h6_m2 - e_6 * fs_3_221_42 * r_8 * h4_m2 - e_6 * fs_2_143_70 * r_10 * h2_m2;
    }

    // NOTE: the angular components are formed in 84 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph2_m1, ph4_m3, ph4_m1, ph6_m3, ph6_m1, ph8_m3, ph8_m1, ph10_m3, ph10_m1, ph12_m3, ph12_m1, ab_2, pc_49 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];

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
        const auto h12_m3 = ph12_m3[k];
        const auto h12_m1 = ph12_m1[k];

        pc_49[k] = - e_1 * fs_945_32_30 * h2_m1 + e_2 * fs_315_16_7 * h4_m3 + e_2 * f_2385_16 * h4_m1 + e_2 * fs_675_16_30 * r_2 * h2_m1 - e_3 * fs_225_44_21 * h6_m3 - e_3 * fs_75_22_210 * h6_m1 - e_3 * fs_315_22_7 * r_2 * h4_m3 - e_3 * f_2385_22 * r_2 * h4_m1 - e_3 * fs_75_4_30 * r_4 * h2_m1 + e_4 * fs_315_1144_462 * h8_m3 + e_4 * fs_2835_1144_10 * h8_m1 + e_4 * fs_45_22_21 * r_2 * h6_m3 + e_4 * fs_15_11_210 * r_2 * h6_m1 + e_4 * fs_945_286_7 * r_4 * h4_m3 + e_4 * f_7155_286 * r_4 * h4_m1 + e_4 * fs_75_22_30 * r_6 * h2_m1 - e_5 * fs_441_8398_429 * h10_m3 - e_5 * fs_378_4199_22 * h10_m1 - e_5 * fs_315_5434_462 * r_2 * h8_m3 - e_5 * fs_2835_5434_10 * r_2 * h8_m1 - e_5 * fs_45_187_21 * r_4 * h6_m3 - e_5 * fs_30_187_210 * r_4 * h6_m1 - e_5 * fs_42_143_7 * r_6 * h4_m3 - e_5 * f_318_143 * r_6 * h4_m1 - e_5 * fs_75_286_30 * r_8 * h2_m1 + e_6 * fs_297_96577_5005 * h12_m3 - e_6 * fs_1089_96577_195 * h12_m1 + e_6 * fs_441_96577_429 * r_2 * h10_m3 + e_6 * fs_756_96577_22 * r_2 * h10_m1 + e_6 * fs_15_5434_462 * r_4 * h8_m3 + e_6 * fs_135_5434_10 * r_4 * h8_m1 + e_6 * fs_30_3553_21 * r_6 * h6_m3 + e_6 * fs_20_3553_210 * r_6 * h6_m1 + e_6 * fs_21_2431_7 * r_8 * h4_m3 + e_6 * f_159_2431 * r_8 * h4_m1 + e_6 * fs_1_143_30 * r_10 * h2_m1;
    }

    // NOTE: the angular components are formed in 84 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_2, pe_3, pe_4, pe_5, pe_6, ph4_m4, ph6_m4, ph8_m4, ph10_m4, ph12_m4, ab_2, pc_50 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

        const auto h4_m4 = ph4_m4[k];
        const auto h6_m4 = ph6_m4[k];
        const auto h8_m4 = ph8_m4[k];
        const auto h10_m4 = ph10_m4[k];
        const auto h12_m4 = ph12_m4[k];

        pc_50[k] = e_2 * fs_315_8_35 * h4_m4 - e_3 * fs_225_11_7 * h6_m4 - e_3 * fs_315_11_35 * r_2 * h4_m4 + e_4 * fs_315_286_77 * h8_m4 + e_4 * fs_90_11_7 * r_2 * h6_m4 + e_4 * fs_945_143_35 * r_4 * h4_m4 - e_5 * fs_126_4199_2145 * h10_m4 - e_5 * fs_630_2717_77 * r_2 * h8_m4 - e_5 * fs_180_187_7 * r_4 * h6_m4 - e_5 * fs_84_143_35 * r_6 * h4_m4 + e_6 * fs_495_96577_2002 * h12_m4 + e_6 * fs_252_96577_2145 * r_2 * h10_m4 + e_6 * fs_30_2717_77 * r_4 * h8_m4 + e_6 * fs_120_3553_7 * r_6 * h6_m4 + e_6 * fs_42_2431_35 * r_8 * h4_m4;
    }

    // NOTE: the angular components are formed in 84 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph2_m1, ph4_m1, ph6_m5, ph6_m1, ph8_m5, ph8_m1, ph10_m5, ph10_m1, ph12_m5, ph12_m1, ab_2, pc_51 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];

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
        const auto h12_m5 = ph12_m5[k];
        const auto h12_m1 = ph12_m1[k];

        pc_51[k] = e_1 * fs_4725_32_3 * h2_m1 - e_2 * fs_1755_32_10 * h4_m1 - e_2 * fs_3375_16_3 * r_2 * h2_m1 - e_3 * fs_75_44_154 * h6_m5 + e_3 * fs_225_44_21 * h6_m1 + e_3 * fs_1755_44_10 * r_2 * h4_m1 + e_3 * fs_375_4_3 * r_4 * h2_m1 + e_4 * fs_105_572_1001 * h8_m5 + e_4 * f_105_22 * h8_m1 + e_4 * fs_15_22_154 * r_2 * h6_m5 - e_4 * fs_45_22_21 * r_2 * h6_m1 - e_4 * fs_405_44_10 * r_4 * h4_m1 - e_4 * fs_375_22_3 * r_6 * h2_m1 - e_5 * fs_315_8398_858 * h10_m5 - e_5 * fs_1701_8398_55 * h10_m1 - e_5 * fs_105_2717_1001 * r_2 * h8_m5 - e_5 * f_210_209 * r_2 * h8_m1 - e_5 * fs_15_187_154 * r_4 * h6_m5 + e_5 * fs_45_187_21 * r_4 * h6_m1 + e_5 * fs_9_11_10 * r_6 * h4_m1 + e_5 * fs_375_286_3 * r_8 * h2_m1 + e_6 * fs_165_96577_17017 * h12_m5 - e_6 * fs_1815_193154_78 * h12_m1 + e_6 * fs_315_96577_858 * r_2 * h10_m5 + e_6 * fs_1701_96577_55 * r_2 * h10_m1 + e_6 * fs_5_2717_1001 * r_4 * h8_m5 + e_6 * f_10_209 * r_4 * h8_m1 + e_6 * fs_10_3553_154 * r_6 * h6_m5 - e_6 * fs_30_3553_21 * r_6 * h6_m1 - e_6 * fs_9_374_10 * r_8 * h4_m1 - e_6 * fs_5_143_3 * r_10 * h2_m1;
    }

    // NOTE: the angular components are formed in 84 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph2_m2, ph4_m2, ph6_m6, ph6_m2, ph8_m6, ph8_m2, ph10_m6, ph10_m2, ph12_m6, ph12_m2, ab_2, pc_52 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];

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
        const auto h12_m6 = ph12_m6[k];
        const auto h12_m2 = ph12_m2[k];

        pc_52[k] = - e_1 * fs_2835_32_10 * h2_m2 - e_2 * fs_1035_32_6 * h4_m2 + e_2 * fs_2025_16_10 * r_2 * h2_m2 + e_3 * fs_75_22_385 * h6_m6 + e_3 * fs_225_11_7 * h6_m2 + e_3 * fs_1035_44_6 * r_2 * h4_m2 - e_3 * fs_225_4_10 * r_4 * h2_m2 - e_4 * fs_105_572_715 * h8_m6 - e_4 * fs_105_143_21 * h8_m2 - e_4 * fs_15_11_385 * r_2 * h6_m6 - e_4 * fs_90_11_7 * r_2 * h6_m2 - e_4 * fs_3105_572_6 * r_4 * h4_m2 + e_4 * fs_225_22_10 * r_6 * h2_m2 - e_5 * fs_63_4199_143 * h10_m6 - e_5 * fs_2709_8398_22 * h10_m2 + e_5 * fs_105_2717_715 * r_2 * h8_m6 + e_5 * fs_420_2717_21 * r_2 * h8_m2 + e_5 * fs_30_187_385 * r_4 * h6_m6 + e_5 * fs_180_187_7 * r_4 * h6_m2 + e_5 * fs_69_143_6 * r_6 * h4_m2 - e_5 * fs_225_286_10 * r_8 * h2_m2 + e_6 * fs_99_96577_36465 * h12_m6 - e_6 * fs_99_193154_10010 * h12_m2 + e_6 * fs_126_96577_143 * r_2 * h10_m6 + e_6 * fs_2709_96577_22 * r_2 * h10_m2 - e_6 * fs_5_2717_715 * r_4 * h8_m6 - e_6 * fs_20_2717_21 * r_4 * h8_m2 - e_6 * fs_20_3553_385 * r_6 * h6_m6 - e_6 * fs_120_3553_7 * r_6 * h6_m2 - e_6 * fs_69_4862_6 * r_8 * h4_m2 + e_6 * fs_3_143_10 * r_10 * h2_m2;
    }

    // NOTE: the angular components are formed in 84 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_2, pe_3, pe_4, pe_5, pe_6, ph4_m3, ph6_m3, ph8_m7, ph8_m3, ph10_m7, ph10_m3, ph12_m7, ph12_m3, ab_2, pc_53 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];

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
        const auto h12_m7 = ph12_m7[k];
        const auto h12_m3 = ph12_m3[k];

        pc_53[k] = e_2 * fs_315_32_462 * h4_m3 + e_3 * fs_75_44_154 * h6_m3 - e_3 * fs_315_44_462 * r_2 * h4_m3 - e_4 * fs_105_52_39 * h8_m7 - e_4 * fs_105_26_7 * h8_m3 - e_4 * fs_15_22_154 * r_2 * h6_m3 + e_4 * fs_945_572_462 * r_4 * h4_m3 + e_5 * fs_441_8398_442 * h10_m7 - e_5 * fs_1953_8398_26 * h10_m3 + e_5 * fs_105_247_39 * r_2 * h8_m7 + e_5 * fs_210_247_7 * r_2 * h8_m3 + e_5 * fs_15_187_154 * r_4 * h6_m3 - e_5 * fs_21_143_462 * r_6 * h4_m3 + e_6 * fs_99_96577_20995 * h12_m7 - e_6 * fs_99_193154_2730 * h12_m3 - e_6 * fs_441_96577_442 * r_2 * h10_m7 + e_6 * fs_1953_96577_26 * r_2 * h10_m3 - e_6 * fs_5_247_39 * r_4 * h8_m7 - e_6 * fs_10_247_7 * r_4 * h8_m3 - e_6 * fs_10_3553_154 * r_6 * h6_m3 + e_6 * fs_21_4862_462 * r_8 * h4_m3;
    }

    // NOTE: the angular components are formed in 84 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_2, pe_3, pe_4, pe_5, pe_6, ph4_m4, ph6_m4, ph8_m8, ph8_m4, ph10_m8, ph10_m4, ph12_m8, ph12_m4, ab_2, pc_54 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

        const auto h4_m4 = ph4_m4[k];
        const auto h6_m4 = ph6_m4[k];
        const auto h8_m8 = ph8_m8[k];
        const auto h8_m4 = ph8_m4[k];
        const auto h10_m8 = ph10_m8[k];
        const auto h10_m4 = ph10_m4[k];
        const auto h12_m8 = ph12_m8[k];
        const auto h12_m4 = ph12_m4[k];

        pc_54[k] = - e_2 * fs_135_16_77 * h4_m4 - e_3 * fs_75_22_385 * h6_m4 + e_3 * fs_135_22_77 * r_2 * h4_m4 + e_4 * fs_105_52_13 * h8_m8 - e_4 * fs_105_52_35 * h8_m4 + e_4 * fs_15_11_385 * r_2 * h6_m4 - e_4 * fs_405_286_77 * r_4 * h4_m4 + e_5 * fs_567_4199_221 * h10_m8 - e_5 * fs_441_4199_39 * h10_m4 - e_5 * fs_105_247_13 * r_2 * h8_m8 + e_5 * fs_105_247_35 * r_2 * h8_m4 - e_5 * fs_30_187_385 * r_4 * h6_m4 + e_5 * fs_18_143_77 * r_6 * h4_m4 + e_6 * fs_33_96577_62985 * h12_m8 - e_6 * fs_33_96577_910 * h12_m4 - e_6 * fs_1134_96577_221 * r_2 * h10_m8 + e_6 * fs_882_96577_39 * r_2 * h10_m4 + e_6 * fs_5_247_13 * r_4 * h8_m8 - e_6 * fs_5_247_35 * r_4 * h8_m4 + e_6 * fs_20_3553_385 * r_6 * h6_m4 - e_6 * fs_9_2431_77 * r_8 * h4_m4;
    }

    // NOTE: the angular components are formed in 84 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph2_0, ph2_p2, ph4_0, ph4_p2, ph6_0, ph6_p2, ph8_0, ph8_p2, ph10_0, ph10_p2, ph12_0, ph12_p2, ab_2, pc_55 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
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
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;
        const auto r_12 = r_2 * r_10;

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

        pc_55[k] = e_0 * f_10395_64 - e_1 * f_12285_32 * h2_0 + e_1 * fs_6615_32_3 * h2_p2 - e_1 * f_10395_16 * r_2 + e_2 * f_180_1 * h4_0 - e_2 * fs_315_4_5 * h4_p2 + e_2 * f_8775_16 * r_2 * h2_0 - e_2 * fs_4725_16_3 * r_2 * h2_p2 + e_2 * f_10395_16 * r_4 - e_3 * f_375_11 * h6_0 + e_3 * fs_75_22_210 * h6_p2 - e_3 * f_1440_11 * r_2 * h4_0 + e_3 * fs_630_11_5 * r_2 * h4_p2 - e_3 * f_975_4 * r_4 * h2_0 + e_3 * fs_525_4_3 * r_4 * h2_p2 - e_3 * f_495_2 * r_6 + e_4 * f_525_286 * h8_0 - e_4 * fs_315_286_70 * h8_p2 + e_4 * f_150_11 * r_2 * h6_0 - e_4 * fs_15_11_210 * r_2 * h6_p2 + e_4 * f_4320_143 * r_4 * h4_0 - e_4 * fs_1890_143_5 * r_4 * h4_p2 + e_4 * f_975_22 * r_6 * h2_0 - e_4 * fs_525_22_3 * r_6 * h2_p2 + e_4 * f_165_4 * r_8 + e_5 * f_189_323 * h10_0 + e_5 * fs_441_4199_165 * h10_p2 - e_5 * f_1050_2717 * r_2 * h8_0 + e_5 * fs_630_2717_70 * r_2 * h8_p2 - e_5 * f_300_187 * r_4 * h6_0 + e_5 * fs_30_187_210 * r_4 * h6_p2 - e_5 * f_384_143 * r_6 * h4_0 + e_5 * fs_168_143_5 * r_6 * h4_p2 - e_5 * f_75_22 * r_8 * h2_0 + e_5 * fs_525_286_3 * r_8 * h2_p2 - e_5 * f_3_1 * r_10 - e_6 * f_26136_96577 * h12_0 - e_6 * fs_396_96577_3003 * h12_p2 - e_6 * f_378_7429 * r_2 * h10_0 - e_6 * fs_882_96577_165 * r_2 * h10_p2 + e_6 * f_50_2717 * r_4 * h8_0 - e_6 * fs_30_2717_70 * r_4 * h8_p2 + e_6 * f_200_3553 * r_6 * h6_0 - e_6 * fs_20_3553_210 * r_6 * h6_p2 + e_6 * f_192_2431 * r_8 * h4_0 - e_6 * fs_84_2431_5 * r_8 * h4_p2 + e_6 * f_1_11 * r_10 * h2_0 - e_6 * fs_7_143_3 * r_10 * h2_p2 + e_6 * f_1_13 * r_12;
    }

    // NOTE: the angular components are formed in 84 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph2_m1, ph4_m1, ph6_m1, ph8_m1, ph10_m1, ph12_m1, ab_2, pc_56 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;

        const auto h2_m1 = ph2_m1[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h8_m1 = ph8_m1[k];
        const auto h10_m1 = ph10_m1[k];
        const auto h12_m1 = ph12_m1[k];

        pc_56[k] = - e_1 * fs_945_32_7 * h2_m1 + e_2 * fs_45_8_210 * h4_m1 + e_2 * fs_675_16_7 * r_2 * h2_m1 - e_3 * f_375_11 * h6_m1 - e_3 * fs_45_11_210 * r_2 * h4_m1 - e_3 * fs_75_4_7 * r_4 * h2_m1 + e_4 * fs_525_286_21 * h8_m1 + e_4 * f_150_11 * r_2 * h6_m1 + e_4 * fs_135_143_210 * r_4 * h4_m1 + e_4 * fs_75_22_7 * r_6 * h2_m1 - e_5 * fs_189_4199_1155 * h10_m1 - e_5 * fs_1050_2717_21 * r_2 * h8_m1 - e_5 * f_300_187 * r_4 * h6_m1 - e_5 * fs_12_143_210 * r_6 * h4_m1 - e_5 * fs_75_286_7 * r_8 * h2_m1 + e_6 * fs_2178_96577_182 * h12_m1 + e_6 * fs_378_96577_1155 * r_2 * h10_m1 + e_6 * fs_50_2717_21 * r_4 * h8_m1 + e_6 * f_200_3553 * r_6 * h6_m1 + e_6 * fs_6_2431_210 * r_8 * h4_m1 + e_6 * fs_1_143_7 * r_10 * h2_m1;
    }

    // NOTE: the angular components are formed in 84 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph2_m2, ph4_m2, ph6_m2, ph8_m2, ph10_m2, ph12_m2, ab_2, pc_57 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;

        const auto h2_m2 = ph2_m2[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h8_m2 = ph8_m2[k];
        const auto h10_m2 = ph10_m2[k];
        const auto h12_m2 = ph12_m2[k];

        pc_57[k] = - e_1 * fs_6615_32_3 * h2_m2 + e_2 * fs_315_4_5 * h4_m2 + e_2 * fs_4725_16_3 * r_2 * h2_m2 - e_3 * fs_75_22_210 * h6_m2 - e_3 * fs_630_11_5 * r_2 * h4_m2 - e_3 * fs_525_4_3 * r_4 * h2_m2 + e_4 * fs_315_286_70 * h8_m2 + e_4 * fs_15_11_210 * r_2 * h6_m2 + e_4 * fs_1890_143_5 * r_4 * h4_m2 + e_4 * fs_525_22_3 * r_6 * h2_m2 - e_5 * fs_441_4199_165 * h10_m2 - e_5 * fs_630_2717_70 * r_2 * h8_m2 - e_5 * fs_30_187_210 * r_4 * h6_m2 - e_5 * fs_168_143_5 * r_6 * h4_m2 - e_5 * fs_525_286_3 * r_8 * h2_m2 + e_6 * fs_396_96577_3003 * h12_m2 + e_6 * fs_882_96577_165 * r_2 * h10_m2 + e_6 * fs_30_2717_70 * r_4 * h8_m2 + e_6 * fs_20_3553_210 * r_6 * h6_m2 + e_6 * fs_84_2431_5 * r_8 * h4_m2 + e_6 * fs_7_143_3 * r_10 * h2_m2;
    }

    // NOTE: the angular components are formed in 84 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph2_m1, ph4_m3, ph4_m1, ph6_m3, ph6_m1, ph8_m3, ph8_m1, ph10_m3, ph10_m1, ph12_m3, ph12_m1, ab_2, pc_58 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];

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
        const auto h12_m3 = ph12_m3[k];
        const auto h12_m1 = ph12_m1[k];

        pc_58[k] = e_1 * fs_945_32_30 * h2_m1 + e_2 * fs_315_16_7 * h4_m3 - e_2 * f_2385_16 * h4_m1 - e_2 * fs_675_16_30 * r_2 * h2_m1 - e_3 * fs_225_44_21 * h6_m3 + e_3 * fs_75_22_210 * h6_m1 - e_3 * fs_315_22_7 * r_2 * h4_m3 + e_3 * f_2385_22 * r_2 * h4_m1 + e_3 * fs_75_4_30 * r_4 * h2_m1 + e_4 * fs_315_1144_462 * h8_m3 - e_4 * fs_2835_1144_10 * h8_m1 + e_4 * fs_45_22_21 * r_2 * h6_m3 - e_4 * fs_15_11_210 * r_2 * h6_m1 + e_4 * fs_945_286_7 * r_4 * h4_m3 - e_4 * f_7155_286 * r_4 * h4_m1 - e_4 * fs_75_22_30 * r_6 * h2_m1 - e_5 * fs_441_8398_429 * h10_m3 + e_5 * fs_378_4199_22 * h10_m1 - e_5 * fs_315_5434_462 * r_2 * h8_m3 + e_5 * fs_2835_5434_10 * r_2 * h8_m1 - e_5 * fs_45_187_21 * r_4 * h6_m3 + e_5 * fs_30_187_210 * r_4 * h6_m1 - e_5 * fs_42_143_7 * r_6 * h4_m3 + e_5 * f_318_143 * r_6 * h4_m1 + e_5 * fs_75_286_30 * r_8 * h2_m1 + e_6 * fs_297_96577_5005 * h12_m3 + e_6 * fs_1089_96577_195 * h12_m1 + e_6 * fs_441_96577_429 * r_2 * h10_m3 - e_6 * fs_756_96577_22 * r_2 * h10_m1 + e_6 * fs_15_5434_462 * r_4 * h8_m3 - e_6 * fs_135_5434_10 * r_4 * h8_m1 + e_6 * fs_30_3553_21 * r_6 * h6_m3 - e_6 * fs_20_3553_210 * r_6 * h6_m1 + e_6 * fs_21_2431_7 * r_8 * h4_m3 - e_6 * f_159_2431 * r_8 * h4_m1 - e_6 * fs_1_143_30 * r_10 * h2_m1;
    }

    // NOTE: the angular components are formed in 84 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph2_m2, ph4_m4, ph4_m2, ph6_m4, ph6_m2, ph8_m4, ph8_m2, ph10_m4, ph10_m2, ph12_m4, ph12_m2, ab_2, pc_59 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];

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
        const auto h12_m4 = ph12_m4[k];
        const auto h12_m2 = ph12_m2[k];

        pc_59[k] = - e_1 * fs_945_16_30 * h2_m2 - e_2 * fs_945_16_14 * h4_m4 + e_2 * fs_135_4_2 * h4_m2 + e_2 * fs_675_8_30 * r_2 * h2_m2 + e_3 * fs_375_88_70 * h6_m4 + e_3 * fs_225_44_21 * h6_m2 + e_3 * fs_945_22_14 * r_2 * h4_m4 - e_3 * fs_270_11_2 * r_2 * h4_m2 - e_3 * fs_75_2_30 * r_4 * h2_m2 - e_4 * fs_105_1144_770 * h8_m4 - e_4 * fs_1995_572_7 * h8_m2 - e_4 * fs_75_44_70 * r_2 * h6_m4 - e_4 * fs_45_22_21 * r_2 * h6_m2 - e_4 * fs_2835_286_14 * r_4 * h4_m4 + e_4 * fs_810_143_2 * r_4 * h4_m2 + e_4 * fs_75_11_30 * r_6 * h2_m2 - e_5 * fs_189_16796_858 * h10_m4 + e_5 * fs_63_442_66 * h10_m2 + e_5 * fs_105_5434_770 * r_2 * h8_m4 + e_5 * fs_105_143_7 * r_2 * h8_m2 + e_5 * fs_75_374_70 * r_4 * h6_m4 + e_5 * fs_45_187_21 * r_4 * h6_m2 + e_5 * fs_126_143_14 * r_6 * h4_m4 - e_5 * fs_72_143_2 * r_6 * h4_m2 - e_5 * fs_75_143_30 * r_8 * h2_m2 + e_6 * fs_264_96577_5005 * h12_m4 + e_6 * fs_66_96577_30030 * h12_m2 + e_6 * fs_189_193154_858 * r_2 * h10_m4 - e_6 * fs_63_5083_66 * r_2 * h10_m2 - e_6 * fs_5_5434_770 * r_4 * h8_m4 - e_6 * fs_5_143_7 * r_4 * h8_m2 - e_6 * fs_25_3553_70 * r_6 * h6_m4 - e_6 * fs_30_3553_21 * r_6 * h6_m2 - e_6 * fs_63_2431_14 * r_8 * h4_m4 + e_6 * fs_36_2431_2 * r_8 * h4_m2 + e_6 * fs_2_143_30 * r_10 * h2_m2;
    }

    // NOTE: the angular components are formed in 84 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_2, pe_3, pe_4, pe_5, pe_6, ph4_m3, ph6_m5, ph6_m3, ph8_m5, ph8_m3, ph10_m5, ph10_m3, ph12_m5, ph12_m3, ab_2, pc_60 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];

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
        const auto h12_m5 = ph12_m5[k];
        const auto h12_m3 = ph12_m3[k];

        pc_60[k] = e_2 * fs_225_16_210 * h4_m3 + e_3 * fs_225_88_462 * h6_m5 - e_3 * fs_375_88_70 * h6_m3 - e_3 * fs_225_22_210 * r_2 * h4_m3 - e_4 * fs_105_572_3003 * h8_m5 - e_4 * fs_105_572_385 * h8_m3 - e_4 * fs_45_44_462 * r_2 * h6_m5 + e_4 * fs_75_44_70 * r_2 * h6_m3 + e_4 * fs_675_286_210 * r_4 * h4_m3 + e_5 * fs_693_16796_286 * h10_m5 + e_5 * fs_693_16796_1430 * h10_m3 + e_5 * fs_105_2717_3003 * r_2 * h8_m5 + e_5 * fs_105_2717_385 * r_2 * h8_m3 + e_5 * fs_45_374_462 * r_4 * h6_m5 - e_5 * fs_75_374_70 * r_4 * h6_m3 - e_5 * fs_30_143_210 * r_6 * h4_m3 + e_6 * fs_66_96577_51051 * h12_m5 + e_6 * fs_99_96577_6006 * h12_m3 - e_6 * fs_693_193154_286 * r_2 * h10_m5 - e_6 * fs_693_193154_1430 * r_2 * h10_m3 - e_6 * fs_5_2717_3003 * r_4 * h8_m5 - e_6 * fs_5_2717_385 * r_4 * h8_m3 - e_6 * fs_15_3553_462 * r_6 * h6_m5 + e_6 * fs_25_3553_70 * r_6 * h6_m3 + e_6 * fs_15_2431_210 * r_8 * h4_m3;
    }

    // NOTE: the angular components are formed in 84 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_2, pe_3, pe_4, pe_5, pe_6, ph4_m4, ph6_m6, ph6_m4, ph8_m6, ph8_m4, ph10_m6, ph10_m4, ph12_m6, ph12_m4, ab_2, pc_61 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];

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
        const auto h12_m6 = ph12_m6[k];
        const auto h12_m4 = ph12_m4[k];

        pc_61[k] = - e_2 * fs_45_16_2310 * h4_m4 - e_3 * fs_75_4_7 * h6_m6 - e_3 * fs_225_88_462 * h6_m4 + e_3 * fs_45_22_2310 * r_2 * h4_m4 - e_4 * fs_105_52_13 * h8_m6 + e_4 * fs_105_104_42 * h8_m4 + e_4 * fs_15_2_7 * r_2 * h6_m6 + e_4 * fs_45_44_462 * r_2 * h6_m4 - e_4 * fs_135_286_2310 * r_4 * h4_m4 + e_5 * fs_63_323_65 * h10_m6 + e_5 * fs_2205_16796_130 * h10_m4 + e_5 * fs_105_247_13 * r_2 * h8_m6 - e_5 * fs_105_494_42 * r_2 * h8_m4 - e_5 * fs_15_17_7 * r_4 * h6_m6 - e_5 * fs_45_374_462 * r_4 * h6_m4 + e_5 * fs_6_143_2310 * r_6 * h4_m4 + e_6 * fs_396_96577_663 * h12_m6 + e_6 * fs_264_96577_273 * h12_m4 - e_6 * fs_126_7429_65 * r_2 * h10_m6 - e_6 * fs_2205_193154_130 * r_2 * h10_m4 - e_6 * fs_5_247_13 * r_4 * h8_m6 + e_6 * fs_5_494_42 * r_4 * h8_m4 + e_6 * fs_10_323_7 * r_6 * h6_m6 + e_6 * fs_15_3553_462 * r_6 * h6_m4 - e_6 * fs_3_2431_2310 * r_8 * h4_m4;
    }

    // NOTE: the angular components are formed in 84 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_3, pe_4, pe_5, pe_6, ph6_m5, ph8_m7, ph8_m5, ph10_m7, ph10_m5, ph12_m7, ph12_m5, ab_2, pc_62 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h6_m5 = ph6_m5[k];
        const auto h8_m7 = ph8_m7[k];
        const auto h8_m5 = ph8_m5[k];
        const auto h10_m7 = ph10_m7[k];
        const auto h10_m5 = ph10_m5[k];
        const auto h12_m7 = ph12_m7[k];
        const auto h12_m5 = ph12_m5[k];

        pc_62[k] = e_3 * fs_75_4_7 * h6_m5 + e_4 * fs_105_104_130 * h8_m7 + e_4 * fs_105_104_182 * h8_m5 - e_4 * fs_15_2_7 * r_2 * h6_m5 + e_5 * fs_126_4199_3315 * h10_m7 + e_5 * fs_1323_8398_39 * h10_m5 - e_5 * fs_105_494_130 * r_2 * h8_m7 - e_5 * fs_105_494_182 * r_2 * h8_m5 + e_5 * fs_15_17_7 * r_4 * h6_m5 + e_6 * fs_33_96577_25194 * h12_m7 + e_6 * fs_33_96577_3094 * h12_m5 - e_6 * fs_252_96577_3315 * r_2 * h10_m7 - e_6 * fs_1323_96577_39 * r_2 * h10_m5 + e_6 * fs_5_494_130 * r_4 * h8_m7 + e_6 * fs_5_494_182 * r_4 * h8_m5 - e_6 * fs_10_323_7 * r_6 * h6_m5;
    }

    // NOTE: the angular components are formed in 84 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph2_0, ph4_0, ph6_0, ph8_0, ph10_0, ph12_0, ab_2, pc_63 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
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
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;
        const auto r_12 = r_2 * r_10;

        const auto h2_0 = ph2_0[k];
        const auto h4_0 = ph4_0[k];
        const auto h6_0 = ph6_0[k];
        const auto h8_0 = ph8_0[k];
        const auto h10_0 = ph10_0[k];
        const auto h12_0 = ph12_0[k];

        pc_63[k] = e_0 * f_10395_64 - e_1 * f_6615_16 * h2_0 - e_1 * f_10395_16 * r_2 + e_2 * f_945_4 * h4_0 + e_2 * f_4725_8 * r_2 * h2_0 + e_2 * f_10395_16 * r_4 - e_3 * f_750_11 * h6_0 - e_3 * f_1890_11 * r_2 * h4_0 - e_3 * f_525_2 * r_4 * h2_0 - e_3 * f_495_2 * r_6 + e_4 * f_3675_286 * h8_0 + e_4 * f_300_11 * r_2 * h6_0 + e_4 * f_5670_143 * r_4 * h4_0 + e_4 * f_525_11 * r_6 * h2_0 + e_4 * f_165_4 * r_8 - e_5 * f_7938_4199 * h10_0 - e_5 * f_7350_2717 * r_2 * h8_0 - e_5 * f_600_187 * r_4 * h6_0 - e_5 * f_504_143 * r_6 * h4_0 - e_5 * f_525_143 * r_8 * h2_0 - e_5 * f_3_1 * r_10 + e_6 * f_30492_96577 * h12_0 + e_6 * f_15876_96577 * r_2 * h10_0 + e_6 * f_350_2717 * r_4 * h8_0 + e_6 * f_400_3553 * r_6 * h6_0 + e_6 * f_252_2431 * r_8 * h4_0 + e_6 * f_14_143 * r_10 * h2_0 + e_6 * f_1_13 * r_12;
    }

    // NOTE: the angular components are formed in 84 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph2_p1, ph4_p1, ph6_p1, ph8_p1, ph10_p1, ph12_p1, ab_2, pc_64 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];

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
        const auto h12_p1 = ph12_p1[k];

        pc_64[k] = - e_1 * fs_945_32_7 * h2_p1 + e_2 * fs_45_8_210 * h4_p1 + e_2 * fs_675_16_7 * r_2 * h2_p1 - e_3 * f_375_11 * h6_p1 - e_3 * fs_45_11_210 * r_2 * h4_p1 - e_3 * fs_75_4_7 * r_4 * h2_p1 + e_4 * fs_525_286_21 * h8_p1 + e_4 * f_150_11 * r_2 * h6_p1 + e_4 * fs_135_143_210 * r_4 * h4_p1 + e_4 * fs_75_22_7 * r_6 * h2_p1 - e_5 * fs_189_4199_1155 * h10_p1 - e_5 * fs_1050_2717_21 * r_2 * h8_p1 - e_5 * f_300_187 * r_4 * h6_p1 - e_5 * fs_12_143_210 * r_6 * h4_p1 - e_5 * fs_75_286_7 * r_8 * h2_p1 + e_6 * fs_2178_96577_182 * h12_p1 + e_6 * fs_378_96577_1155 * r_2 * h10_p1 + e_6 * fs_50_2717_21 * r_4 * h8_p1 + e_6 * f_200_3553 * r_6 * h6_p1 + e_6 * fs_6_2431_210 * r_8 * h4_p1 + e_6 * fs_1_143_7 * r_10 * h2_p1;
    }

    // NOTE: the angular components are formed in 84 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph2_p2, ph4_p2, ph4_p3, ph6_p2, ph6_p3, ph8_p2, ph8_p3, ph10_p2, ph10_p3, ph12_p2, ph12_p3, ab_2, pc_65, pc_66 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];

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
        const auto h12_p2 = ph12_p2[k];
        const auto h12_p3 = ph12_p3[k];

        pc_65[k] = e_1 * fs_945_16_70 * h2_p2 - e_2 * fs_495_16_42 * h4_p2 - e_2 * fs_675_8_70 * r_2 * h2_p2 + e_3 * f_75_2 * h6_p2 + e_3 * fs_45_2_42 * r_2 * h4_p2 + e_3 * fs_75_2_70 * r_4 * h2_p2 - e_4 * fs_735_572_3 * h8_p2 - e_4 * f_15_1 * r_2 * h6_p2 - e_4 * fs_135_26_42 * r_4 * h4_p2 - e_4 * fs_75_11_70 * r_6 * h2_p2 - e_5 * fs_189_4199_154 * h10_p2 + e_5 * fs_735_2717_3 * r_2 * h8_p2 + e_5 * f_30_17 * r_4 * h6_p2 + e_5 * fs_6_13_42 * r_6 * h4_p2 + e_5 * fs_75_143_70 * r_8 * h2_p2 + e_6 * fs_693_96577_1430 * h12_p2 + e_6 * fs_378_96577_154 * r_2 * h10_p2 - e_6 * fs_35_2717_3 * r_4 * h8_p2 - e_6 * f_20_323 * r_6 * h6_p2 - e_6 * fs_3_221_42 * r_8 * h4_p2 - e_6 * fs_2_143_70 * r_10 * h2_p2;

        pc_66[k] = - e_2 * fs_945_8_3 * h4_p3 + e_3 * f_3225_44 * h6_p3 + e_3 * fs_945_11_3 * r_2 * h4_p3 - e_4 * fs_735_286_22 * h8_p3 - e_4 * f_645_22 * r_2 * h6_p3 - e_4 * fs_2835_143_3 * r_4 * h4_p3 + e_5 * fs_189_8398_1001 * h10_p3 + e_5 * fs_1470_2717_22 * r_2 * h8_p3 + e_5 * f_645_187 * r_4 * h6_p3 + e_5 * fs_252_143_3 * r_6 * h4_p3 + e_6 * fs_462_96577_2145 * h12_p3 - e_6 * fs_189_96577_1001 * r_2 * h10_p3 - e_6 * fs_70_2717_22 * r_4 * h8_p3 - e_6 * f_430_3553 * r_6 * h6_p3 - e_6 * fs_126_2431_3 * r_8 * h4_p3;
    }

    // NOTE: the angular components are formed in 84 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_2, pe_3, pe_4, pe_5, pe_6, ph4_p4, ph6_p4, ph6_p5, ph8_p4, ph10_p4, ph10_p5, ph12_p4, ph12_p5, ab_2, pc_67, pc_68 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

        const auto h4_p4 = ph4_p4[k];
        const auto h6_p4 = ph6_p4[k];
        const auto h6_p5 = ph6_p5[k];
        const auto h8_p4 = ph8_p4[k];
        const auto h10_p4 = ph10_p4[k];
        const auto h10_p5 = ph10_p5[k];
        const auto h12_p4 = ph12_p4[k];
        const auto h12_p5 = ph12_p5[k];

        pc_67[k] = e_2 * fs_945_8_5 * h4_p4 + e_3 * f_150_11 * h6_p4 - e_3 * fs_945_11_5 * r_2 * h4_p4 - e_4 * fs_2205_572_11 * h8_p4 - e_4 * f_60_11 * r_2 * h6_p4 + e_4 * fs_2835_143_5 * r_4 * h4_p4 + e_5 * fs_63_4199_15015 * h10_p4 + e_5 * fs_2205_2717_11 * r_2 * h8_p4 + e_5 * f_120_187 * r_4 * h6_p4 - e_5 * fs_252_143_5 * r_6 * h4_p4 + e_6 * fs_924_96577_286 * h12_p4 - e_6 * fs_126_96577_15015 * r_2 * h10_p4 - e_6 * fs_105_2717_11 * r_4 * h8_p4 - e_6 * f_80_3553 * r_6 * h6_p4 + e_6 * fs_126_2431_5 * r_8 * h4_p4;

        pc_68[k] = - e_3 * f_375_4 * h6_p5 + e_4 * f_75_2 * r_2 * h6_p5 + e_5 * fs_63_442_273 * h10_p5 - e_5 * f_75_17 * r_4 * h6_p5 + e_6 * fs_462_96577_442 * h12_p5 - e_6 * fs_63_5083_273 * r_2 * h10_p5 + e_6 * f_50_323 * r_6 * h6_p5;
    }

    // NOTE: the angular components are formed in 84 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_3, pe_4, pe_5, pe_6, ph6_p6, ph8_p6, ph10_p6, ph12_p6, ab_2, pc_69 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h6_p6 = ph6_p6[k];
        const auto h8_p6 = ph8_p6[k];
        const auto h10_p6 = ph10_p6[k];
        const auto h12_p6 = ph12_p6[k];

        pc_69[k] = e_3 * f_75_2 * h6_p6 + e_4 * fs_105_52_91 * h8_p6 - e_4 * f_15_1 * r_2 * h6_p6 + e_5 * fs_378_4199_455 * h10_p6 - e_5 * fs_105_247_91 * r_2 * h8_p6 + e_5 * f_30_17 * r_4 * h6_p6 + e_6 * fs_66_96577_4641 * h12_p6 - e_6 * fs_756_96577_455 * r_2 * h10_p6 + e_6 * fs_5_247_91 * r_4 * h8_p6 - e_6 * f_20_323 * r_6 * h6_p6;
    }

    // NOTE: the angular components are formed in 84 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph2_0, ph2_p2, ph4_0, ph4_p2, ph6_0, ph6_p2, ph8_0, ph8_p2, ph10_0, ph10_p2, ph12_0, ph12_p2, ab_2, pc_70 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
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
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;
        const auto r_12 = r_2 * r_10;

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

        pc_70[k] = e_0 * f_10395_64 - e_1 * f_12285_32 * h2_0 - e_1 * fs_6615_32_3 * h2_p2 - e_1 * f_10395_16 * r_2 + e_2 * f_180_1 * h4_0 + e_2 * fs_315_4_5 * h4_p2 + e_2 * f_8775_16 * r_2 * h2_0 + e_2 * fs_4725_16_3 * r_2 * h2_p2 + e_2 * f_10395_16 * r_4 - e_3 * f_375_11 * h6_0 - e_3 * fs_75_22_210 * h6_p2 - e_3 * f_1440_11 * r_2 * h4_0 - e_3 * fs_630_11_5 * r_2 * h4_p2 - e_3 * f_975_4 * r_4 * h2_0 - e_3 * fs_525_4_3 * r_4 * h2_p2 - e_3 * f_495_2 * r_6 + e_4 * f_525_286 * h8_0 + e_4 * fs_315_286_70 * h8_p2 + e_4 * f_150_11 * r_2 * h6_0 + e_4 * fs_15_11_210 * r_2 * h6_p2 + e_4 * f_4320_143 * r_4 * h4_0 + e_4 * fs_1890_143_5 * r_4 * h4_p2 + e_4 * f_975_22 * r_6 * h2_0 + e_4 * fs_525_22_3 * r_6 * h2_p2 + e_4 * f_165_4 * r_8 + e_5 * f_189_323 * h10_0 - e_5 * fs_441_4199_165 * h10_p2 - e_5 * f_1050_2717 * r_2 * h8_0 - e_5 * fs_630_2717_70 * r_2 * h8_p2 - e_5 * f_300_187 * r_4 * h6_0 - e_5 * fs_30_187_210 * r_4 * h6_p2 - e_5 * f_384_143 * r_6 * h4_0 - e_5 * fs_168_143_5 * r_6 * h4_p2 - e_5 * f_75_22 * r_8 * h2_0 - e_5 * fs_525_286_3 * r_8 * h2_p2 - e_5 * f_3_1 * r_10 - e_6 * f_26136_96577 * h12_0 + e_6 * fs_396_96577_3003 * h12_p2 - e_6 * f_378_7429 * r_2 * h10_0 + e_6 * fs_882_96577_165 * r_2 * h10_p2 + e_6 * f_50_2717 * r_4 * h8_0 + e_6 * fs_30_2717_70 * r_4 * h8_p2 + e_6 * f_200_3553 * r_6 * h6_0 + e_6 * fs_20_3553_210 * r_6 * h6_p2 + e_6 * f_192_2431 * r_8 * h4_0 + e_6 * fs_84_2431_5 * r_8 * h4_p2 + e_6 * f_1_11 * r_10 * h2_0 + e_6 * fs_7_143_3 * r_10 * h2_p2 + e_6 * f_1_13 * r_12;
    }

    // NOTE: the angular components are formed in 84 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph2_p1, ph4_p1, ph4_p3, ph6_p1, ph6_p3, ph8_p1, ph8_p3, ph10_p1, ph10_p3, ph12_p1, ph12_p3, ab_2, pc_71 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];

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
        const auto h12_p1 = ph12_p1[k];
        const auto h12_p3 = ph12_p3[k];

        pc_71[k] = - e_1 * fs_945_32_30 * h2_p1 + e_2 * f_2385_16 * h4_p1 + e_2 * fs_315_16_7 * h4_p3 + e_2 * fs_675_16_30 * r_2 * h2_p1 - e_3 * fs_75_22_210 * h6_p1 - e_3 * fs_225_44_21 * h6_p3 - e_3 * f_2385_22 * r_2 * h4_p1 - e_3 * fs_315_22_7 * r_2 * h4_p3 - e_3 * fs_75_4_30 * r_4 * h2_p1 + e_4 * fs_2835_1144_10 * h8_p1 + e_4 * fs_315_1144_462 * h8_p3 + e_4 * fs_15_11_210 * r_2 * h6_p1 + e_4 * fs_45_22_21 * r_2 * h6_p3 + e_4 * f_7155_286 * r_4 * h4_p1 + e_4 * fs_945_286_7 * r_4 * h4_p3 + e_4 * fs_75_22_30 * r_6 * h2_p1 - e_5 * fs_378_4199_22 * h10_p1 - e_5 * fs_441_8398_429 * h10_p3 - e_5 * fs_2835_5434_10 * r_2 * h8_p1 - e_5 * fs_315_5434_462 * r_2 * h8_p3 - e_5 * fs_30_187_210 * r_4 * h6_p1 - e_5 * fs_45_187_21 * r_4 * h6_p3 - e_5 * f_318_143 * r_6 * h4_p1 - e_5 * fs_42_143_7 * r_6 * h4_p3 - e_5 * fs_75_286_30 * r_8 * h2_p1 - e_6 * fs_1089_96577_195 * h12_p1 + e_6 * fs_297_96577_5005 * h12_p3 + e_6 * fs_756_96577_22 * r_2 * h10_p1 + e_6 * fs_441_96577_429 * r_2 * h10_p3 + e_6 * fs_135_5434_10 * r_4 * h8_p1 + e_6 * fs_15_5434_462 * r_4 * h8_p3 + e_6 * fs_20_3553_210 * r_6 * h6_p1 + e_6 * fs_30_3553_21 * r_6 * h6_p3 + e_6 * f_159_2431 * r_8 * h4_p1 + e_6 * fs_21_2431_7 * r_8 * h4_p3 + e_6 * fs_1_143_30 * r_10 * h2_p1;
    }

    // NOTE: the angular components are formed in 84 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph2_p2, ph4_p2, ph4_p4, ph6_p2, ph6_p4, ph8_p2, ph8_p4, ph10_p2, ph10_p4, ph12_p2, ph12_p4, ab_2, pc_72 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];

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
        const auto h12_p2 = ph12_p2[k];
        const auto h12_p4 = ph12_p4[k];

        pc_72[k] = e_1 * fs_945_16_30 * h2_p2 - e_2 * fs_135_4_2 * h4_p2 - e_2 * fs_945_16_14 * h4_p4 - e_2 * fs_675_8_30 * r_2 * h2_p2 - e_3 * fs_225_44_21 * h6_p2 + e_3 * fs_375_88_70 * h6_p4 + e_3 * fs_270_11_2 * r_2 * h4_p2 + e_3 * fs_945_22_14 * r_2 * h4_p4 + e_3 * fs_75_2_30 * r_4 * h2_p2 + e_4 * fs_1995_572_7 * h8_p2 - e_4 * fs_105_1144_770 * h8_p4 + e_4 * fs_45_22_21 * r_2 * h6_p2 - e_4 * fs_75_44_70 * r_2 * h6_p4 - e_4 * fs_810_143_2 * r_4 * h4_p2 - e_4 * fs_2835_286_14 * r_4 * h4_p4 - e_4 * fs_75_11_30 * r_6 * h2_p2 - e_5 * fs_63_442_66 * h10_p2 - e_5 * fs_189_16796_858 * h10_p4 - e_5 * fs_105_143_7 * r_2 * h8_p2 + e_5 * fs_105_5434_770 * r_2 * h8_p4 - e_5 * fs_45_187_21 * r_4 * h6_p2 + e_5 * fs_75_374_70 * r_4 * h6_p4 + e_5 * fs_72_143_2 * r_6 * h4_p2 + e_5 * fs_126_143_14 * r_6 * h4_p4 + e_5 * fs_75_143_30 * r_8 * h2_p2 - e_6 * fs_66_96577_30030 * h12_p2 + e_6 * fs_264_96577_5005 * h12_p4 + e_6 * fs_63_5083_66 * r_2 * h10_p2 + e_6 * fs_189_193154_858 * r_2 * h10_p4 + e_6 * fs_5_143_7 * r_4 * h8_p2 - e_6 * fs_5_5434_770 * r_4 * h8_p4 + e_6 * fs_30_3553_21 * r_6 * h6_p2 - e_6 * fs_25_3553_70 * r_6 * h6_p4 - e_6 * fs_36_2431_2 * r_8 * h4_p2 - e_6 * fs_63_2431_14 * r_8 * h4_p4 - e_6 * fs_2_143_30 * r_10 * h2_p2;
    }

    // NOTE: the angular components are formed in 84 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_2, pe_3, pe_4, pe_5, pe_6, ph4_p3, ph6_p3, ph6_p5, ph8_p3, ph8_p5, ph10_p3, ph10_p5, ph12_p3, ph12_p5, ab_2, pc_73 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];

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
        const auto h12_p3 = ph12_p3[k];
        const auto h12_p5 = ph12_p5[k];

        pc_73[k] = - e_2 * fs_225_16_210 * h4_p3 + e_3 * fs_375_88_70 * h6_p3 + e_3 * fs_225_88_462 * h6_p5 + e_3 * fs_225_22_210 * r_2 * h4_p3 + e_4 * fs_105_572_385 * h8_p3 - e_4 * fs_105_572_3003 * h8_p5 - e_4 * fs_75_44_70 * r_2 * h6_p3 - e_4 * fs_45_44_462 * r_2 * h6_p5 - e_4 * fs_675_286_210 * r_4 * h4_p3 - e_5 * fs_693_16796_1430 * h10_p3 + e_5 * fs_693_16796_286 * h10_p5 - e_5 * fs_105_2717_385 * r_2 * h8_p3 + e_5 * fs_105_2717_3003 * r_2 * h8_p5 + e_5 * fs_75_374_70 * r_4 * h6_p3 + e_5 * fs_45_374_462 * r_4 * h6_p5 + e_5 * fs_30_143_210 * r_6 * h4_p3 - e_6 * fs_99_96577_6006 * h12_p3 + e_6 * fs_66_96577_51051 * h12_p5 + e_6 * fs_693_193154_1430 * r_2 * h10_p3 - e_6 * fs_693_193154_286 * r_2 * h10_p5 + e_6 * fs_5_2717_385 * r_4 * h8_p3 - e_6 * fs_5_2717_3003 * r_4 * h8_p5 - e_6 * fs_25_3553_70 * r_6 * h6_p3 - e_6 * fs_15_3553_462 * r_6 * h6_p5 - e_6 * fs_15_2431_210 * r_8 * h4_p3;
    }

    // NOTE: the angular components are formed in 84 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_2, pe_3, pe_4, pe_5, pe_6, ph4_p4, ph6_p4, ph6_p6, ph8_p4, ph8_p6, ph10_p4, ph10_p6, ph12_p4, ph12_p6, ab_2, pc_74 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];

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
        const auto h12_p4 = ph12_p4[k];
        const auto h12_p6 = ph12_p6[k];

        pc_74[k] = e_2 * fs_45_16_2310 * h4_p4 + e_3 * fs_225_88_462 * h6_p4 - e_3 * fs_75_4_7 * h6_p6 - e_3 * fs_45_22_2310 * r_2 * h4_p4 - e_4 * fs_105_104_42 * h8_p4 - e_4 * fs_105_52_13 * h8_p6 - e_4 * fs_45_44_462 * r_2 * h6_p4 + e_4 * fs_15_2_7 * r_2 * h6_p6 + e_4 * fs_135_286_2310 * r_4 * h4_p4 - e_5 * fs_2205_16796_130 * h10_p4 + e_5 * fs_63_323_65 * h10_p6 + e_5 * fs_105_494_42 * r_2 * h8_p4 + e_5 * fs_105_247_13 * r_2 * h8_p6 + e_5 * fs_45_374_462 * r_4 * h6_p4 - e_5 * fs_15_17_7 * r_4 * h6_p6 - e_5 * fs_6_143_2310 * r_6 * h4_p4 - e_6 * fs_264_96577_273 * h12_p4 + e_6 * fs_396_96577_663 * h12_p6 + e_6 * fs_2205_193154_130 * r_2 * h10_p4 - e_6 * fs_126_7429_65 * r_2 * h10_p6 - e_6 * fs_5_494_42 * r_4 * h8_p4 - e_6 * fs_5_247_13 * r_4 * h8_p6 - e_6 * fs_15_3553_462 * r_6 * h6_p4 + e_6 * fs_10_323_7 * r_6 * h6_p6 + e_6 * fs_3_2431_2310 * r_8 * h4_p4;
    }

    // NOTE: the angular components are formed in 84 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_3, pe_4, pe_5, pe_6, ph6_p5, ph8_p5, ph8_p7, ph10_p5, ph10_p7, ph12_p5, ph12_p7, ab_2, pc_75 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h6_p5 = ph6_p5[k];
        const auto h8_p5 = ph8_p5[k];
        const auto h8_p7 = ph8_p7[k];
        const auto h10_p5 = ph10_p5[k];
        const auto h10_p7 = ph10_p7[k];
        const auto h12_p5 = ph12_p5[k];
        const auto h12_p7 = ph12_p7[k];

        pc_75[k] = - e_3 * fs_75_4_7 * h6_p5 - e_4 * fs_105_104_182 * h8_p5 + e_4 * fs_105_104_130 * h8_p7 + e_4 * fs_15_2_7 * r_2 * h6_p5 - e_5 * fs_1323_8398_39 * h10_p5 + e_5 * fs_126_4199_3315 * h10_p7 + e_5 * fs_105_494_182 * r_2 * h8_p5 - e_5 * fs_105_494_130 * r_2 * h8_p7 - e_5 * fs_15_17_7 * r_4 * h6_p5 - e_6 * fs_33_96577_3094 * h12_p5 + e_6 * fs_33_96577_25194 * h12_p7 + e_6 * fs_1323_96577_39 * r_2 * h10_p5 - e_6 * fs_252_96577_3315 * r_2 * h10_p7 - e_6 * fs_5_494_182 * r_4 * h8_p5 + e_6 * fs_5_494_130 * r_4 * h8_p7 + e_6 * fs_10_323_7 * r_6 * h6_p5;
    }

    // NOTE: the angular components are formed in 84 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph2_0, ph4_0, ph4_p4, ph6_0, ph6_p4, ph8_0, ph8_p4, ph10_0, ph10_p4, ph12_0, ph12_p4, ab_2, pc_76 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
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
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;
        const auto r_12 = r_2 * r_10;

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

        pc_76[k] = e_0 * f_10395_64 - e_1 * f_4725_16 * h2_0 - e_1 * f_10395_16 * r_2 + e_2 * f_495_16 * h4_0 + e_2 * fs_315_8_35 * h4_p4 + e_2 * f_3375_8 * r_2 * h2_0 + e_2 * f_10395_16 * r_4 + e_3 * f_75_2 * h6_0 - e_3 * fs_225_11_7 * h6_p4 - e_3 * f_45_2 * r_2 * h4_0 - e_3 * fs_315_11_35 * r_2 * h4_p4 - e_3 * f_375_2 * r_4 * h2_0 - e_3 * f_495_2 * r_6 - e_4 * f_7455_572 * h8_0 + e_4 * fs_315_286_77 * h8_p4 - e_4 * f_15_1 * r_2 * h6_0 + e_4 * fs_90_11_7 * r_2 * h6_p4 + e_4 * f_135_26 * r_4 * h4_0 + e_4 * fs_945_143_35 * r_4 * h4_p4 + e_4 * f_375_11 * r_6 * h2_0 + e_4 * f_165_4 * r_8 + e_5 * f_6615_4199 * h10_0 - e_5 * fs_126_4199_2145 * h10_p4 + e_5 * f_7455_2717 * r_2 * h8_0 - e_5 * fs_630_2717_77 * r_2 * h8_p4 + e_5 * f_30_17 * r_4 * h6_0 - e_5 * fs_180_187_7 * r_4 * h6_p4 - e_5 * f_6_13 * r_6 * h4_0 - e_5 * fs_84_143_35 * r_6 * h4_p4 - e_5 * f_375_143 * r_8 * h2_0 - e_5 * f_3_1 * r_10 + e_6 * f_16335_96577 * h12_0 + e_6 * fs_495_96577_2002 * h12_p4 - e_6 * f_13230_96577 * r_2 * h10_0 + e_6 * fs_252_96577_2145 * r_2 * h10_p4 - e_6 * f_355_2717 * r_4 * h8_0 + e_6 * fs_30_2717_77 * r_4 * h8_p4 - e_6 * f_20_323 * r_6 * h6_0 + e_6 * fs_120_3553_7 * r_6 * h6_p4 + e_6 * f_3_221 * r_8 * h4_0 + e_6 * fs_42_2431_35 * r_8 * h4_p4 + e_6 * f_10_143 * r_10 * h2_0 + e_6 * f_1_13 * r_12;
    }

    // NOTE: the angular components are formed in 84 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph2_p1, ph4_p1, ph6_p1, ph6_p5, ph8_p1, ph8_p5, ph10_p1, ph10_p5, ph12_p1, ph12_p5, ab_2, pc_77 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];

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
        const auto h12_p1 = ph12_p1[k];
        const auto h12_p5 = ph12_p5[k];

        pc_77[k] = - e_1 * fs_4725_32_3 * h2_p1 + e_2 * fs_1755_32_10 * h4_p1 + e_2 * fs_3375_16_3 * r_2 * h2_p1 - e_3 * fs_225_44_21 * h6_p1 - e_3 * fs_75_44_154 * h6_p5 - e_3 * fs_1755_44_10 * r_2 * h4_p1 - e_3 * fs_375_4_3 * r_4 * h2_p1 - e_4 * f_105_22 * h8_p1 + e_4 * fs_105_572_1001 * h8_p5 + e_4 * fs_45_22_21 * r_2 * h6_p1 + e_4 * fs_15_22_154 * r_2 * h6_p5 + e_4 * fs_405_44_10 * r_4 * h4_p1 + e_4 * fs_375_22_3 * r_6 * h2_p1 + e_5 * fs_1701_8398_55 * h10_p1 - e_5 * fs_315_8398_858 * h10_p5 + e_5 * f_210_209 * r_2 * h8_p1 - e_5 * fs_105_2717_1001 * r_2 * h8_p5 - e_5 * fs_45_187_21 * r_4 * h6_p1 - e_5 * fs_15_187_154 * r_4 * h6_p5 - e_5 * fs_9_11_10 * r_6 * h4_p1 - e_5 * fs_375_286_3 * r_8 * h2_p1 + e_6 * fs_1815_193154_78 * h12_p1 + e_6 * fs_165_96577_17017 * h12_p5 - e_6 * fs_1701_96577_55 * r_2 * h10_p1 + e_6 * fs_315_96577_858 * r_2 * h10_p5 - e_6 * f_10_209 * r_4 * h8_p1 + e_6 * fs_5_2717_1001 * r_4 * h8_p5 + e_6 * fs_30_3553_21 * r_6 * h6_p1 + e_6 * fs_10_3553_154 * r_6 * h6_p5 + e_6 * fs_9_374_10 * r_8 * h4_p1 + e_6 * fs_5_143_3 * r_10 * h2_p1;
    }

    // NOTE: the angular components are formed in 84 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph2_p2, ph4_p2, ph6_p2, ph6_p6, ph8_p2, ph8_p6, ph10_p2, ph10_p6, ph12_p2, ph12_p6, ab_2, pc_78 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];

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
        const auto h12_p2 = ph12_p2[k];
        const auto h12_p6 = ph12_p6[k];

        pc_78[k] = e_1 * fs_2835_32_10 * h2_p2 + e_2 * fs_1035_32_6 * h4_p2 - e_2 * fs_2025_16_10 * r_2 * h2_p2 - e_3 * fs_225_11_7 * h6_p2 + e_3 * fs_75_22_385 * h6_p6 - e_3 * fs_1035_44_6 * r_2 * h4_p2 + e_3 * fs_225_4_10 * r_4 * h2_p2 + e_4 * fs_105_143_21 * h8_p2 - e_4 * fs_105_572_715 * h8_p6 + e_4 * fs_90_11_7 * r_2 * h6_p2 - e_4 * fs_15_11_385 * r_2 * h6_p6 + e_4 * fs_3105_572_6 * r_4 * h4_p2 - e_4 * fs_225_22_10 * r_6 * h2_p2 + e_5 * fs_2709_8398_22 * h10_p2 - e_5 * fs_63_4199_143 * h10_p6 - e_5 * fs_420_2717_21 * r_2 * h8_p2 + e_5 * fs_105_2717_715 * r_2 * h8_p6 - e_5 * fs_180_187_7 * r_4 * h6_p2 + e_5 * fs_30_187_385 * r_4 * h6_p6 - e_5 * fs_69_143_6 * r_6 * h4_p2 + e_5 * fs_225_286_10 * r_8 * h2_p2 + e_6 * fs_99_193154_10010 * h12_p2 + e_6 * fs_99_96577_36465 * h12_p6 - e_6 * fs_2709_96577_22 * r_2 * h10_p2 + e_6 * fs_126_96577_143 * r_2 * h10_p6 + e_6 * fs_20_2717_21 * r_4 * h8_p2 - e_6 * fs_5_2717_715 * r_4 * h8_p6 + e_6 * fs_120_3553_7 * r_6 * h6_p2 - e_6 * fs_20_3553_385 * r_6 * h6_p6 + e_6 * fs_69_4862_6 * r_8 * h4_p2 - e_6 * fs_3_143_10 * r_10 * h2_p2;
    }

    // NOTE: the angular components are formed in 84 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_2, pe_3, pe_4, pe_5, pe_6, ph4_p3, ph6_p3, ph8_p3, ph8_p7, ph10_p3, ph10_p7, ph12_p3, ph12_p7, ab_2, pc_79 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];

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
        const auto h12_p3 = ph12_p3[k];
        const auto h12_p7 = ph12_p7[k];

        pc_79[k] = - e_2 * fs_315_32_462 * h4_p3 - e_3 * fs_75_44_154 * h6_p3 + e_3 * fs_315_44_462 * r_2 * h4_p3 + e_4 * fs_105_26_7 * h8_p3 - e_4 * fs_105_52_39 * h8_p7 + e_4 * fs_15_22_154 * r_2 * h6_p3 - e_4 * fs_945_572_462 * r_4 * h4_p3 + e_5 * fs_1953_8398_26 * h10_p3 + e_5 * fs_441_8398_442 * h10_p7 - e_5 * fs_210_247_7 * r_2 * h8_p3 + e_5 * fs_105_247_39 * r_2 * h8_p7 - e_5 * fs_15_187_154 * r_4 * h6_p3 + e_5 * fs_21_143_462 * r_6 * h4_p3 + e_6 * fs_99_193154_2730 * h12_p3 + e_6 * fs_99_96577_20995 * h12_p7 - e_6 * fs_1953_96577_26 * r_2 * h10_p3 - e_6 * fs_441_96577_442 * r_2 * h10_p7 + e_6 * fs_10_247_7 * r_4 * h8_p3 - e_6 * fs_5_247_39 * r_4 * h8_p7 + e_6 * fs_10_3553_154 * r_6 * h6_p3 - e_6 * fs_21_4862_462 * r_8 * h4_p3;
    }

    // NOTE: the angular components are formed in 84 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_2, pe_3, pe_4, pe_5, pe_6, ph4_p4, ph6_p4, ph8_p4, ph8_p8, ph10_p4, ph10_p8, ph12_p4, ph12_p8, ab_2, pc_80 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

        const auto h4_p4 = ph4_p4[k];
        const auto h6_p4 = ph6_p4[k];
        const auto h8_p4 = ph8_p4[k];
        const auto h8_p8 = ph8_p8[k];
        const auto h10_p4 = ph10_p4[k];
        const auto h10_p8 = ph10_p8[k];
        const auto h12_p4 = ph12_p4[k];
        const auto h12_p8 = ph12_p8[k];

        pc_80[k] = e_2 * fs_135_16_77 * h4_p4 + e_3 * fs_75_22_385 * h6_p4 - e_3 * fs_135_22_77 * r_2 * h4_p4 + e_4 * fs_105_52_35 * h8_p4 + e_4 * fs_105_52_13 * h8_p8 - e_4 * fs_15_11_385 * r_2 * h6_p4 + e_4 * fs_405_286_77 * r_4 * h4_p4 + e_5 * fs_441_4199_39 * h10_p4 + e_5 * fs_567_4199_221 * h10_p8 - e_5 * fs_105_247_35 * r_2 * h8_p4 - e_5 * fs_105_247_13 * r_2 * h8_p8 + e_5 * fs_30_187_385 * r_4 * h6_p4 - e_5 * fs_18_143_77 * r_6 * h4_p4 + e_6 * fs_33_96577_910 * h12_p4 + e_6 * fs_33_96577_62985 * h12_p8 - e_6 * fs_882_96577_39 * r_2 * h10_p4 - e_6 * fs_1134_96577_221 * r_2 * h10_p8 + e_6 * fs_5_247_35 * r_4 * h8_p4 + e_6 * fs_5_247_13 * r_4 * h8_p8 - e_6 * fs_20_3553_385 * r_6 * h6_p4 + e_6 * fs_9_2431_77 * r_8 * h4_p4;
    }

    // NOTE: the angular components are formed in 84 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph2_0, ph4_0, ph6_0, ph6_p6, ph8_0, ph8_p6, ph10_0, ph10_p6, ph12_0, ph12_p6, ab_2, pc_81 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
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
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;
        const auto r_12 = r_2 * r_10;

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

        pc_81[k] = e_0 * f_10395_64 - e_1 * f_4725_32 * h2_0 - e_1 * f_10395_16 * r_2 - e_2 * f_1215_8 * h4_0 + e_2 * f_3375_16 * r_2 * h2_0 + e_2 * f_10395_16 * r_4 + e_3 * f_3225_44 * h6_0 - e_3 * fs_75_22_462 * h6_p6 + e_3 * f_1215_11 * r_2 * h4_0 - e_3 * f_375_4 * r_4 * h2_0 - e_3 * f_495_2 * r_6 - e_4 * f_1995_572 * h8_0 + e_4 * fs_105_286_858 * h8_p6 - e_4 * f_645_22 * r_2 * h6_0 + e_4 * fs_15_11_462 * r_2 * h6_p6 - e_4 * f_3645_143 * r_4 * h4_0 + e_4 * f_375_22 * r_6 * h2_0 + e_4 * f_165_4 * r_8 - e_5 * f_945_442 * h10_0 - e_5 * fs_189_8398_4290 * h10_p6 + e_5 * f_105_143 * r_2 * h8_0 - e_5 * fs_210_2717_858 * r_2 * h8_p6 + e_5 * f_645_187 * r_4 * h6_0 - e_5 * fs_30_187_462 * r_4 * h6_p6 + e_5 * f_324_143 * r_6 * h4_0 - e_5 * f_375_286 * r_8 * h2_0 - e_5 * f_3_1 * r_10 - e_6 * f_7260_96577 * h12_0 + e_6 * fs_330_96577_4862 * h12_p6 + e_6 * f_945_5083 * r_2 * h10_0 + e_6 * fs_189_96577_4290 * r_2 * h10_p6 - e_6 * f_5_143 * r_4 * h8_0 + e_6 * fs_10_2717_858 * r_4 * h8_p6 - e_6 * f_430_3553 * r_6 * h6_0 + e_6 * fs_20_3553_462 * r_6 * h6_p6 - e_6 * f_162_2431 * r_8 * h4_0 + e_6 * f_5_143 * r_10 * h2_0 + e_6 * f_1_13 * r_12;
    }

    // NOTE: the angular components are formed in 84 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph2_p1, ph4_p1, ph6_p1, ph8_p1, ph8_p7, ph10_p1, ph10_p7, ph12_p1, ph12_p7, ab_2, pc_82 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];

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
        const auto h12_p1 = ph12_p1[k];
        const auto h12_p7 = ph12_p7[k];

        pc_82[k] = - e_1 * fs_6615_64_10 * h2_p1 + e_2 * fs_945_16_3 * h4_p1 + e_2 * fs_4725_32_10 * r_2 * h2_p1 + e_3 * fs_375_88_70 * h6_p1 - e_3 * fs_945_22_3 * r_2 * h4_p1 - e_3 * fs_525_8_10 * r_4 * h2_p1 - e_4 * fs_945_572_30 * h8_p1 + e_4 * fs_105_572_858 * h8_p7 - e_4 * fs_75_44_70 * r_2 * h6_p1 + e_4 * fs_2835_286_3 * r_4 * h4_p1 + e_4 * fs_525_44_10 * r_6 * h2_p1 - e_5 * fs_2583_16796_66 * h10_p1 - e_5 * fs_189_8398_2431 * h10_p7 + e_5 * fs_945_2717_30 * r_2 * h8_p1 - e_5 * fs_105_2717_858 * r_2 * h8_p7 + e_5 * fs_75_374_70 * r_4 * h6_p1 - e_5 * fs_126_143_3 * r_6 * h4_p1 - e_5 * fs_525_572_10 * r_8 * h2_p1 - e_6 * fs_363_96577_65 * h12_p1 + e_6 * fs_33_96577_461890 * h12_p7 + e_6 * fs_2583_193154_66 * r_2 * h10_p1 + e_6 * fs_189_96577_2431 * r_2 * h10_p7 - e_6 * fs_45_2717_30 * r_4 * h8_p1 + e_6 * fs_5_2717_858 * r_4 * h8_p7 - e_6 * fs_25_3553_70 * r_6 * h6_p1 + e_6 * fs_63_2431_3 * r_8 * h4_p1 + e_6 * fs_7_286_10 * r_10 * h2_p1;
    }

    // NOTE: the angular components are formed in 84 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph2_p2, ph4_p2, ph6_p2, ph8_p2, ph8_p8, ph10_p2, ph10_p8, ph12_p2, ph12_p8, ab_2, pc_83 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];

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
        const auto h12_p2 = ph12_p2[k];
        const auto h12_p8 = ph12_p8[k];

        pc_83[k] = e_1 * fs_945_32_55 * h2_p2 + e_2 * fs_135_4_33 * h4_p2 - e_2 * fs_675_16_55 * r_2 * h2_p2 - e_3 * fs_75_44_154 * h6_p2 - e_3 * fs_270_11_33 * r_2 * h4_p2 + e_3 * fs_75_4_55 * r_4 * h2_p2 - e_4 * fs_315_572_462 * h8_p2 - e_4 * fs_105_52_39 * h8_p8 + e_4 * fs_15_22_154 * r_2 * h6_p2 + e_4 * fs_810_143_33 * r_4 * h4_p2 - e_4 * fs_75_22_55 * r_6 * h2_p2 - e_5 * f_6993_8398 * h10_p2 + e_5 * fs_63_8398_663 * h10_p8 + e_5 * fs_315_2717_462 * r_2 * h8_p2 + e_5 * fs_105_247_39 * r_2 * h8_p8 - e_5 * fs_15_187_154 * r_4 * h6_p2 - e_5 * fs_72_143_33 * r_6 * h4_p2 + e_5 * fs_75_286_55 * r_8 * h2_p2 - e_6 * fs_66_96577_455 * h12_p2 + e_6 * fs_132_96577_20995 * h12_p8 + e_6 * f_6993_96577 * r_2 * h10_p2 - e_6 * fs_63_96577_663 * r_2 * h10_p8 - e_6 * fs_15_2717_462 * r_4 * h8_p2 - e_6 * fs_5_247_39 * r_4 * h8_p8 + e_6 * fs_10_3553_154 * r_6 * h6_p2 + e_6 * fs_36_2431_33 * r_8 * h4_p2 - e_6 * fs_1_143_55 * r_10 * h2_p2;
    }

    // NOTE: the angular components are formed in 84 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_2, pe_3, pe_4, pe_5, pe_6, ph4_p3, ph6_p3, ph8_p3, ph10_p3, ph10_p9, ph12_p3, ph12_p9, ab_2, pc_84 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;

        const auto h4_p3 = ph4_p3[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h8_p3 = ph8_p3[k];
        const auto h10_p3 = ph10_p3[k];
        const auto h10_p9 = ph10_p9[k];
        const auto h12_p3 = ph12_p3[k];
        const auto h12_p9 = ph12_p9[k];

        pc_84[k] = - e_2 * fs_405_32_154 * h4_p3 - e_3 * fs_75_22_462 * h6_p3 + e_3 * fs_405_44_154 * r_2 * h4_p3 - e_4 * fs_105_52_21 * h8_p3 + e_4 * fs_15_11_462 * r_2 * h6_p3 - e_4 * fs_1215_572_154 * r_4 * h4_p3 - e_5 * fs_189_4199_78 * h10_p3 + e_5 * fs_189_8398_8398 * h10_p9 + e_5 * fs_105_247_21 * r_2 * h8_p3 - e_5 * fs_30_187_462 * r_4 * h6_p3 + e_5 * fs_27_143_154 * r_6 * h4_p3 - e_6 * fs_33_193154_910 * h12_p3 + e_6 * fs_33_96577_146965 * h12_p9 + e_6 * fs_378_96577_78 * r_2 * h10_p3 - e_6 * fs_189_96577_8398 * r_2 * h10_p9 - e_6 * fs_5_247_21 * r_4 * h8_p3 + e_6 * fs_20_3553_462 * r_6 * h6_p3 - e_6 * fs_27_4862_154 * r_8 * h4_p3;
    }

    // NOTE: the angular components are formed in 84 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph2_0, ph4_0, ph6_0, ph8_0, ph8_p8, ph10_0, ph10_p8, ph12_0, ph12_p8, ab_2, pc_85 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
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
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;
        const auto r_12 = r_2 * r_10;

        const auto h2_0 = ph2_0[k];
        const auto h4_0 = ph4_0[k];
        const auto h6_0 = ph6_0[k];
        const auto h8_0 = ph8_0[k];
        const auto h8_p8 = ph8_p8[k];
        const auto h10_0 = ph10_0[k];
        const auto h10_p8 = ph10_p8[k];
        const auto h12_0 = ph12_0[k];
        const auto h12_p8 = ph12_p8[k];

        pc_85[k] = e_0 * f_10395_64 + e_1 * f_945_16 * h2_0 - e_1 * f_10395_16 * r_2 - e_2 * f_270_1 * h4_0 - e_2 * f_675_8 * r_2 * h2_0 + e_2 * f_10395_16 * r_4 + e_3 * f_150_11 * h6_0 + e_3 * f_2160_11 * r_2 * h4_0 + e_3 * f_75_2 * r_4 * h2_0 - e_3 * f_495_2 * r_6 + e_4 * f_9345_572 * h8_0 + e_4 * fs_315_572_715 * h8_p8 - e_4 * f_60_11 * r_2 * h6_0 - e_4 * f_6480_143 * r_4 * h4_0 - e_4 * f_75_11 * r_6 * h2_0 + e_4 * f_165_4 * r_8 + e_5 * f_5229_4199 * h10_0 - e_5 * fs_63_4199_12155 * h10_p8 - e_5 * f_9345_2717 * r_2 * h8_0 - e_5 * fs_315_2717_715 * r_2 * h8_p8 + e_5 * f_120_187 * r_4 * h6_0 + e_5 * f_576_143 * r_6 * h4_0 + e_5 * f_75_143 * r_8 * h2_0 - e_5 * f_3_1 * r_10 + e_6 * f_2178_96577 * h12_0 + e_6 * fs_66_96577_138567 * h12_p8 - e_6 * f_10458_96577 * r_2 * h10_0 + e_6 * fs_126_96577_12155 * r_2 * h10_p8 + e_6 * f_445_2717 * r_4 * h8_0 + e_6 * fs_15_2717_715 * r_4 * h8_p8 - e_6 * f_80_3553 * r_6 * h6_0 - e_6 * f_288_2431 * r_8 * h4_0 - e_6 * f_2_143 * r_10 * h2_0 + e_6 * f_1_13 * r_12;
    }

    // NOTE: the angular components are formed in 84 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph2_p1, ph4_p1, ph6_p1, ph8_p1, ph10_p1, ph10_p9, ph12_p1, ph12_p9, ab_2, pc_86 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];

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
        const auto h12_p1 = ph12_p1[k];
        const auto h12_p9 = ph12_p9[k];

        pc_86[k] = - e_1 * fs_2835_64_66 * h2_p1 - e_2 * fs_135_16_55 * h4_p1 + e_2 * fs_2025_32_66 * r_2 * h2_p1 + e_3 * fs_225_88_462 * h6_p1 + e_3 * fs_135_22_55 * r_2 * h4_p1 - e_3 * fs_225_8_66 * r_4 * h2_p1 + e_4 * fs_315_143_22 * h8_p1 - e_4 * fs_45_44_462 * r_2 * h6_p1 - e_4 * fs_405_286_55 * r_4 * h4_p1 + e_4 * fs_225_44_66 * r_6 * h2_p1 + e_5 * fs_2709_16796_10 * h10_p1 - e_5 * fs_63_8398_20995 * h10_p9 - e_5 * fs_1260_2717_22 * r_2 * h8_p1 + e_5 * fs_45_374_462 * r_4 * h6_p1 + e_5 * fs_18_143_55 * r_6 * h4_p1 - e_5 * fs_225_572_66 * r_8 * h2_p1 + e_6 * fs_33_96577_429 * h12_p1 + e_6 * fs_99_96577_58786 * h12_p9 - e_6 * fs_2709_193154_10 * r_2 * h10_p1 + e_6 * fs_63_96577_20995 * r_2 * h10_p9 + e_6 * fs_60_2717_22 * r_4 * h8_p1 - e_6 * fs_15_3553_462 * r_6 * h6_p1 - e_6 * fs_9_2431_55 * r_8 * h4_p1 + e_6 * fs_3_286_66 * r_10 * h2_p1;
    }

    // NOTE: the angular components are formed in 84 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph2_p2, ph4_p2, ph6_p2, ph8_p2, ph10_p2, ph10_p10, ph12_p2, ph12_p10, ab_2, pc_87 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];

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
        const auto h12_p2 = ph12_p2[k];
        const auto h12_p10 = ph12_p10[k];

        pc_87[k] = e_1 * fs_945_32_22 * h2_p2 + e_2 * fs_405_32_330 * h4_p2 - e_2 * fs_675_16_22 * r_2 * h2_p2 + e_3 * fs_75_22_385 * h6_p2 - e_3 * fs_405_44_330 * r_2 * h4_p2 + e_3 * fs_75_4_22 * r_4 * h2_p2 + e_4 * fs_105_572_1155 * h8_p2 - e_4 * fs_15_11_385 * r_2 * h6_p2 + e_4 * fs_1215_572_330 * r_4 * h4_p2 - e_4 * fs_75_22_22 * r_6 * h2_p2 + e_5 * fs_567_8398_10 * h10_p2 + e_5 * fs_63_4199_12597 * h10_p10 - e_5 * fs_105_2717_1155 * r_2 * h8_p2 + e_5 * fs_30_187_385 * r_4 * h6_p2 - e_5 * fs_27_143_330 * r_6 * h4_p2 + e_5 * fs_75_286_22 * r_8 * h2_p2 + e_6 * fs_33_193154_182 * h12_p2 + e_6 * fs_33_96577_323323 * h12_p10 - e_6 * fs_567_96577_10 * r_2 * h10_p2 - e_6 * fs_126_96577_12597 * r_2 * h10_p10 + e_6 * fs_5_2717_1155 * r_4 * h8_p2 - e_6 * fs_20_3553_385 * r_6 * h6_p2 + e_6 * fs_27_4862_330 * r_8 * h4_p2 - e_6 * fs_1_143_22 * r_10 * h2_p2;
    }

    // NOTE: the angular components are formed in 84 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph2_0, ph4_0, ph6_0, ph8_0, ph10_0, ph10_p10, ph12_0, ph12_p10, ab_2, pc_88 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
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
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;
        const auto r_12 = r_2 * r_10;

        const auto h2_0 = ph2_0[k];
        const auto h4_0 = ph4_0[k];
        const auto h6_0 = ph6_0[k];
        const auto h8_0 = ph8_0[k];
        const auto h10_0 = ph10_0[k];
        const auto h10_p10 = ph10_p10[k];
        const auto h12_0 = ph12_0[k];
        const auto h12_p10 = ph12_p10[k];

        pc_88[k] = e_0 * f_10395_64 + e_1 * f_10395_32 * h2_0 - e_1 * f_10395_16 * r_2 - e_2 * f_1485_8 * h4_0 - e_2 * f_7425_16 * r_2 * h2_0 + e_2 * f_10395_16 * r_4 - e_3 * f_375_4 * h6_0 + e_3 * f_135_1 * r_2 * h4_0 + e_3 * f_825_4 * r_4 * h2_0 - e_3 * f_495_2 * r_6 - e_4 * f_525_52 * h8_0 + e_4 * f_75_2 * r_2 * h6_0 - e_4 * f_405_13 * r_4 * h4_0 - e_4 * f_75_2 * r_6 * h2_0 + e_4 * f_165_4 * r_8 - e_5 * f_3087_8398 * h10_0 - e_5 * fs_63_8398_92378 * h10_p10 + e_5 * f_525_247 * r_2 * h8_0 - e_5 * f_75_17 * r_4 * h6_0 + e_5 * f_36_13 * r_6 * h4_0 + e_5 * f_75_26 * r_8 * h2_0 - e_5 * f_3_1 * r_10 - e_6 * f_396_96577 * h12_0 + e_6 * fs_66_96577_176358 * h12_p10 + e_6 * f_3087_96577 * r_2 * h10_0 + e_6 * fs_63_96577_92378 * r_2 * h10_p10 - e_6 * f_25_247 * r_4 * h8_0 + e_6 * f_50_323 * r_6 * h6_0 - e_6 * f_18_221 * r_8 * h4_0 - e_6 * f_1_13 * r_10 * h2_0 + e_6 * f_1_13 * r_12;
    }

    // NOTE: the angular components are formed in 84 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph2_p1, ph4_p1, ph6_p1, ph8_p1, ph10_p1, ph12_p1, ph12_p11, ab_2, pc_89 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];
        const auto e_6 = pe_6[k];

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
        const auto h12_p1 = ph12_p1[k];
        const auto h12_p11 = ph12_p11[k];

        pc_89[k] = - e_1 * f_10395_32 * h2_p1 - e_2 * fs_1485_32_30 * h4_p1 + e_2 * f_7425_16 * r_2 * h2_p1 - e_3 * fs_75_4_7 * h6_p1 + e_3 * fs_135_4_30 * r_2 * h4_p1 - e_3 * f_825_4 * r_4 * h2_p1 - e_4 * fs_105_52_3 * h8_p1 + e_4 * fs_15_2_7 * r_2 * h6_p1 - e_4 * fs_405_52_30 * r_4 * h4_p1 + e_4 * f_75_2 * r_6 * h2_p1 - e_5 * fs_63_8398_165 * h10_p1 + e_5 * fs_105_247_3 * r_2 * h8_p1 - e_5 * fs_15_17_7 * r_4 * h6_p1 + e_5 * fs_9_13_30 * r_6 * h4_p1 - e_5 * f_75_26 * r_8 * h2_p1 - e_6 * fs_33_193154_26 * h12_p1 + e_6 * fs_33_96577_676039 * h12_p11 + e_6 * fs_63_96577_165 * r_2 * h10_p1 - e_6 * fs_5_247_3 * r_4 * h8_p1 + e_6 * fs_10_323_7 * r_6 * h6_p1 - e_6 * fs_9_442_30 * r_8 * h4_p1 + e_6 * f_1_13 * r_10 * h2_p1;
    }

    // NOTE: the angular components are formed in 84 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, pe_6, ph2_0, ph4_0, ph6_0, ph8_0, ph10_0, ph12_0, ph12_p12, ab_2, pc_90 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
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
        const auto r_6 = r_2 * r_4;
        const auto r_8 = r_2 * r_6;
        const auto r_10 = r_2 * r_8;
        const auto r_12 = r_2 * r_10;

        const auto h2_0 = ph2_0[k];
        const auto h4_0 = ph4_0[k];
        const auto h6_0 = ph6_0[k];
        const auto h8_0 = ph8_0[k];
        const auto h10_0 = ph10_0[k];
        const auto h12_0 = ph12_0[k];
        const auto h12_p12 = ph12_p12[k];

        pc_90[k] = e_0 * f_10395_64 + e_1 * f_10395_16 * h2_0 - e_1 * f_10395_16 * r_2 + e_2 * f_4455_16 * h4_0 - e_2 * f_7425_8 * r_2 * h2_0 + e_2 * f_10395_16 * r_4 + e_3 * f_75_2 * h6_0 - e_3 * f_405_2 * r_2 * h4_0 + e_3 * f_825_2 * r_4 * h2_0 - e_3 * f_495_2 * r_6 + e_4 * f_105_52 * h8_0 - e_4 * f_15_1 * r_2 * h6_0 + e_4 * f_1215_26 * r_4 * h4_0 - e_4 * f_75_1 * r_6 * h2_0 + e_4 * f_165_4 * r_8 + e_5 * f_189_4199 * h10_0 - e_5 * f_105_247 * r_2 * h8_0 + e_5 * f_30_17 * r_4 * h6_0 - e_5 * f_54_13 * r_6 * h4_0 + e_5 * f_75_13 * r_8 * h2_0 - e_5 * f_3_1 * r_10 + e_6 * f_33_96577 * h12_0 + e_6 * fs_33_96577_1352078 * h12_p12 - e_6 * f_378_96577 * r_2 * h10_0 + e_6 * f_5_247 * r_4 * h8_0 - e_6 * f_20_323 * r_6 * h6_0 + e_6 * f_27_221 * r_8 * h4_0 - e_6 * f_2_13 * r_10 * h2_0 + e_6 * f_1_13 * r_12;
    }

    // NOTE: the values of a combination of angular components are stored as one
    // row of nvalues columns, with the component on bra side running slowest. The
    // rows which the symmetry relates to an already formed one are copied from it.

    const size_t sources[169] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 1, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 2, 15, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 3, 16, 29, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 4, 17, 30, 43, 56, 57, 58, 59, 60, 61, 62, 63, 64, 5, 18, 31, 44, 57, 70, 71, 72, 73, 74, 75, 76, 77, 6, 19, 32, 45, 58, 71, 84, 85, 86, 87, 88, 89, 90, 7, 20, 33, 46, 59, 72, 85, 98, 99, 100, 101, 102, 103, 8, 21, 34, 47, 60, 73, 86, 99, 112, 113, 114, 115, 116, 9, 22, 35, 48, 61, 74, 87, 100, 113, 126, 127, 128, 129, 10, 23, 36, 49, 62, 75, 88, 101, 114, 127, 140, 141, 142, 11, 24, 37, 50, 63, 76, 89, 102, 115, 128, 141, 154, 155, 12, 25, 38, 51, 64, 77, 90, 103, 116, 129, 142, 155, 168};

    for (size_t n = 0; n < 169; n++)
    {
        if (sources[n] != n) std::copy(values + sources[n] * nvalues, values + (sources[n] + 1) * nvalues, values + n * nvalues);
    }
}

}  // namespace simdt2ceri
