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



#include "SimdTwoCenterElectronRepulsionRecHI.hpp"

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
compute_hi_electron_repulsion(double                         *values,
                              const size_t                    nvalues,
                              const CBasisFunction           &bra,
                              const CBasisFunction           &ket,
                              const std::vector<CSimdMatrix> &harmonics,
                              const CSimdMatrix              &coordinates) -> void
{
    if ((bra.get_angular_momentum() != 5) || (ket.get_angular_momentum() != 6))
    {
        errors::assertMsgCritical(
            false, std::string("SimdTwoCenterElectronRepulsionFunc.compute_hi_electron_repulsion: Basis functions must be of angular momenta five and six"));
    }

    if (harmonics.size() < 11)
    {
        errors::assertMsgCritical(
            false, std::string("SimdTwoCenterElectronRepulsionFunc.compute_hi_electron_repulsion: Harmonics must reach angular momentum 11"));
    }

    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("SimdTwoCenterElectronRepulsionFunc.compute_hi_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    auto boys = CSimdVariableMatrix(std::vector<size_t>(nprims, nvalues), 13);

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
    // call, which fills the orders 0 to 11 of every row. The terms read the
    // orders 6 to 11 alone, and the orders below them are formed on the
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

            const auto ff_0 = fbase * aexp / (fexp * fexp * fexp * fexp * fexp * fexp);

            const auto ff_1 = fbase * aexp * fmu / (fexp * fexp * fexp * fexp * fexp * fexp);

            const auto ff_2 = fbase * aexp * fmu * fmu / (fexp * fexp * fexp * fexp * fexp * fexp);

            const auto ff_3 = fbase * aexp * fmu * fmu * fmu / (fexp * fexp * fexp * fexp * fexp * fexp);

            const auto ff_4 = fbase * aexp * fmu * fmu * fmu * fmu / (fexp * fexp * fexp * fexp * fexp * fexp);

            const auto ff_5 = fbase * aexp * fmu * fmu * fmu * fmu * fmu / (fexp * fexp * fexp * fexp * fexp * fexp);

            const auto *bv_0 = boys.data(7, i * nprim_b + j);

            const auto *bv_1 = boys.data(8, i * nprim_b + j);

            const auto *bv_2 = boys.data(9, i * nprim_b + j);

            const auto *bv_3 = boys.data(10, i * nprim_b + j);

            const auto *bv_4 = boys.data(11, i * nprim_b + j);

            const auto *bv_5 = boys.data(12, i * nprim_b + j);

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
    const auto *ph11_m11 = harmonics[10].data(0);
    const auto *ph11_m10 = harmonics[10].data(1);
    const auto *ph11_m9 = harmonics[10].data(2);
    const auto *ph11_m8 = harmonics[10].data(3);
    const auto *ph11_m7 = harmonics[10].data(4);
    const auto *ph11_m6 = harmonics[10].data(5);
    const auto *ph11_m5 = harmonics[10].data(6);
    const auto *ph11_m4 = harmonics[10].data(7);
    const auto *ph11_m3 = harmonics[10].data(8);
    const auto *ph11_m2 = harmonics[10].data(9);
    const auto *ph11_m1 = harmonics[10].data(10);
    const auto *ph11_0 = harmonics[10].data(11);
    const auto *ph11_p1 = harmonics[10].data(12);
    const auto *ph11_p2 = harmonics[10].data(13);
    const auto *ph11_p3 = harmonics[10].data(14);
    const auto *ph11_p4 = harmonics[10].data(15);
    const auto *ph11_p5 = harmonics[10].data(16);
    const auto *ph11_p6 = harmonics[10].data(17);
    const auto *ph11_p7 = harmonics[10].data(18);
    const auto *ph11_p8 = harmonics[10].data(19);
    const auto *ph11_p9 = harmonics[10].data(20);
    const auto *ph11_p10 = harmonics[10].data(21);
    const auto *ph11_p11 = harmonics[10].data(22);

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
    auto *pc_117 = values + 117 * nvalues;
    auto *pc_118 = values + 118 * nvalues;
    auto *pc_119 = values + 119 * nvalues;
    auto *pc_120 = values + 120 * nvalues;
    auto *pc_121 = values + 121 * nvalues;
    auto *pc_122 = values + 122 * nvalues;
    auto *pc_123 = values + 123 * nvalues;
    auto *pc_124 = values + 124 * nvalues;
    auto *pc_125 = values + 125 * nvalues;
    auto *pc_126 = values + 126 * nvalues;
    auto *pc_127 = values + 127 * nvalues;
    auto *pc_128 = values + 128 * nvalues;
    auto *pc_129 = values + 129 * nvalues;
    auto *pc_130 = values + 130 * nvalues;
    auto *pc_131 = values + 131 * nvalues;
    auto *pc_132 = values + 132 * nvalues;
    auto *pc_133 = values + 133 * nvalues;
    auto *pc_134 = values + 134 * nvalues;
    auto *pc_135 = values + 135 * nvalues;
    auto *pc_136 = values + 136 * nvalues;
    auto *pc_137 = values + 137 * nvalues;
    auto *pc_138 = values + 138 * nvalues;
    auto *pc_139 = values + 139 * nvalues;
    auto *pc_140 = values + 140 * nvalues;
    auto *pc_141 = values + 141 * nvalues;
    auto *pc_142 = values + 142 * nvalues;

    // NOTE: the factors of the terms depend on the angular momenta alone, so they
    // are formed once for the whole matrix instead of once for every atom pair.

    const auto f_105_1 = 105.0;
    const auto f_105_2 = 105.0 / 2.0;
    const auto f_10_221 = 10.0 / 221.0;
    const auto f_12_13 = 12.0 / 13.0;
    const auto f_1386_4199 = 1386.0 / 4199.0;
    const auto f_138_2431 = 138.0 / 2431.0;
    const auto f_1449_2431 = 1449.0 / 2431.0;
    const auto f_145_4 = 145.0 / 4.0;
    const auto f_14_143 = 14.0 / 143.0;
    const auto f_1575_143 = 1575.0 / 143.0;
    const auto f_15_1 = 15.0;
    const auto f_15_13 = 15.0 / 13.0;
    const auto f_18_143 = 18.0 / 143.0;
    const auto f_210_143 = 210.0 / 143.0;
    const auto f_2205_16 = 2205.0 / 16.0;
    const auto f_225_26 = 225.0 / 26.0;
    const auto f_245_2 = 245.0 / 2.0;
    const auto f_24_221 = 24.0 / 221.0;
    const auto f_270_13 = 270.0 / 13.0;
    const auto f_2835_16 = 2835.0 / 16.0;
    const auto f_2835_8 = 2835.0 / 8.0;
    const auto f_29_13 = 29.0 / 13.0;
    const auto f_300_13 = 300.0 / 13.0;
    const auto f_315_11 = 315.0 / 11.0;
    const auto f_315_22 = 315.0 / 22.0;
    const auto f_32_221 = 32.0 / 221.0;
    const auto f_360_13 = 360.0 / 13.0;
    const auto f_36_13 = 36.0 / 13.0;
    const auto f_405_2 = 405.0 / 2.0;
    const auto f_40_13 = 40.0 / 13.0;
    const auto f_420_143 = 420.0 / 143.0;
    const auto f_420_2431 = 420.0 / 2431.0;
    const auto f_435_26 = 435.0 / 26.0;
    const auto f_4410_2431 = 4410.0 / 2431.0;
    const auto f_45_1 = 45.0;
    const auto f_45_11 = 45.0 / 11.0;
    const auto f_48_13 = 48.0 / 13.0;
    const auto f_490_143 = 490.0 / 143.0;
    const auto f_49_429 = 49.0 / 429.0;
    const auto f_50_1 = 50.0;
    const auto f_58_663 = 58.0 / 663.0;
    const auto f_60_1 = 60.0;
    const auto f_6300_2431 = 6300.0 / 2431.0;
    const auto f_6300_46189 = 6300.0 / 46189.0;
    const auto f_735_22 = 735.0 / 22.0;
    const auto f_75_4 = 75.0 / 4.0;
    const auto f_7_143 = 7.0 / 143.0;
    const auto f_80_663 = 80.0 / 663.0;
    const auto f_8_221 = 8.0 / 221.0;
    const auto f_90_13 = 90.0 / 13.0;
    const auto f_945_16 = 945.0 / 16.0;
    const auto f_945_8 = 945.0 / 8.0;
    const auto fs_1029_4862_33 = std::sqrt(3176523.0 / 2149004.0);
    const auto fs_1035_1144_14 = std::sqrt(7498575.0 / 654368.0);
    const auto fs_1035_4862_14 = std::sqrt(7498575.0 / 11819522.0);
    const auto fs_1035_92378_14 = std::sqrt(7498575.0 / 4266847442.0);
    const auto fs_1050_2431_11 = std::sqrt(1102500.0 / 537251.0);
    const auto fs_1050_2431_14 = std::sqrt(15435000.0 / 5909761.0);
    const auto fs_1050_2431_35 = std::sqrt(38587500.0 / 5909761.0);
    const auto fs_1050_2431_66 = std::sqrt(6615000.0 / 537251.0);
    const auto fs_1050_46189_11 = std::sqrt(1102500.0 / 193947611.0);
    const auto fs_1050_46189_14 = std::sqrt(15435000.0 / 2133423721.0);
    const auto fs_1050_46189_35 = std::sqrt(38587500.0 / 2133423721.0);
    const auto fs_1050_46189_66 = std::sqrt(6615000.0 / 193947611.0);
    const auto fs_105_104_130 = std::sqrt(55125.0 / 416.0);
    const auto fs_105_11_2 = std::sqrt(22050.0 / 121.0);
    const auto fs_105_11_21 = std::sqrt(231525.0 / 121.0);
    const auto fs_105_11_7 = std::sqrt(77175.0 / 121.0);
    const auto fs_105_143_10 = std::sqrt(110250.0 / 20449.0);
    const auto fs_105_143_22 = std::sqrt(22050.0 / 1859.0);
    const auto fs_105_22_14 = std::sqrt(77175.0 / 242.0);
    const auto fs_105_22_15 = std::sqrt(165375.0 / 484.0);
    const auto fs_105_22_21 = std::sqrt(231525.0 / 484.0);
    const auto fs_105_22_3 = std::sqrt(33075.0 / 484.0);
    const auto fs_105_22_33 = std::sqrt(33075.0 / 44.0);
    const auto fs_105_22_35 = std::sqrt(385875.0 / 484.0);
    const auto fs_105_22_42 = std::sqrt(231525.0 / 242.0);
    const auto fs_105_22_6 = std::sqrt(33075.0 / 242.0);
    const auto fs_105_2431_231 = std::sqrt(231525.0 / 537251.0);
    const auto fs_105_2431_4290 = std::sqrt(330750.0 / 41327.0);
    const auto fs_105_2431_715 = std::sqrt(55125.0 / 41327.0);
    const auto fs_105_286_165 = std::sqrt(165375.0 / 7436.0);
    const auto fs_105_286_21 = std::sqrt(231525.0 / 81796.0);
    const auto fs_105_286_35 = std::sqrt(385875.0 / 81796.0);
    const auto fs_105_2_11 = std::sqrt(121275.0 / 4.0);
    const auto fs_105_2_3 = std::sqrt(33075.0 / 4.0);
    const auto fs_105_2_5 = std::sqrt(55125.0 / 4.0);
    const auto fs_105_2_7 = std::sqrt(77175.0 / 4.0);
    const auto fs_105_442_130 = std::sqrt(55125.0 / 7514.0);
    const auto fs_105_442_91 = std::sqrt(77175.0 / 15028.0);
    const auto fs_105_44_30 = std::sqrt(165375.0 / 968.0);
    const auto fs_105_44_6 = std::sqrt(33075.0 / 968.0);
    const auto fs_105_44_66 = std::sqrt(33075.0 / 88.0);
    const auto fs_105_46189_4290 = std::sqrt(330750.0 / 14919047.0);
    const auto fs_105_46189_715 = std::sqrt(55125.0 / 14919047.0);
    const auto fs_105_4862_7293 = std::sqrt(33075.0 / 9724.0);
    const auto fs_105_4862_858 = std::sqrt(33075.0 / 82654.0);
    const auto fs_105_4_10 = std::sqrt(55125.0 / 8.0);
    const auto fs_105_4_22 = std::sqrt(121275.0 / 8.0);
    const auto fs_105_572_4290 = std::sqrt(165375.0 / 1144.0);
    const auto fs_105_572_715 = std::sqrt(55125.0 / 2288.0);
    const auto fs_105_8398_130 = std::sqrt(55125.0 / 2712554.0);
    const auto fs_105_8398_286 = std::sqrt(121275.0 / 2712554.0);
    const auto fs_107_2431_10 = std::sqrt(114490.0 / 5909761.0);
    const auto fs_108_4199_154 = std::sqrt(1796256.0 / 17631601.0);
    const auto fs_10_13_2 = std::sqrt(200.0 / 169.0);
    const auto fs_10_221_11 = std::sqrt(1100.0 / 48841.0);
    const auto fs_10_221_7 = std::sqrt(700.0 / 48841.0);
    const auto fs_10_2431_231 = std::sqrt(2100.0 / 537251.0);
    const auto fs_10_663_42 = std::sqrt(1400.0 / 146523.0);
    const auto fs_1113_4862_33 = std::sqrt(3716307.0 / 2149004.0);
    const auto fs_1197_9724_154 = std::sqrt(10029663.0 / 4298008.0);
    const auto fs_11_13_5 = std::sqrt(605.0 / 169.0);
    const auto fs_11_221_2 = std::sqrt(242.0 / 48841.0);
    const auto fs_11_26_42 = std::sqrt(2541.0 / 338.0);
    const auto fs_11_663_42 = std::sqrt(1694.0 / 146523.0);
    const auto fs_1200_2431_33 = std::sqrt(4320000.0 / 537251.0);
    const auto fs_1200_46189_33 = std::sqrt(4320000.0 / 193947611.0);
    const auto fs_120_13_3 = std::sqrt(43200.0 / 169.0);
    const auto fs_120_2431_2310 = std::sqrt(3024000.0 / 537251.0);
    const auto fs_120_46189_2310 = std::sqrt(3024000.0 / 193947611.0);
    const auto fs_125_4_3 = std::sqrt(46875.0 / 16.0);
    const auto fs_1260_2431_30 = std::sqrt(47628000.0 / 5909761.0);
    const auto fs_1260_46189_30 = std::sqrt(47628000.0 / 2133423721.0);
    const auto fs_126_2431_70 = std::sqrt(1111320.0 / 5909761.0);
    const auto fs_126_4199_110 = std::sqrt(1746360.0 / 17631601.0);
    const auto fs_1275_286_3 = std::sqrt(4876875.0 / 81796.0);
    const auto fs_12_13_5 = std::sqrt(720.0 / 169.0);
    const auto fs_12_13_7 = std::sqrt(1008.0 / 169.0);
    const auto fs_12_143_2 = std::sqrt(288.0 / 20449.0);
    const auto fs_12_2431_70 = std::sqrt(10080.0 / 5909761.0);
    const auto fs_12_4199_2145 = std::sqrt(23760.0 / 1356277.0);
    const auto fs_12_4199_6006 = std::sqrt(66528.0 / 1356277.0);
    const auto fs_12_4199_91 = std::sqrt(1008.0 / 1356277.0);
    const auto fs_135_1_2 = std::sqrt(36450.0);
    const auto fs_135_2431_105 = std::sqrt(1913625.0 / 5909761.0);
    const auto fs_135_2431_385 = std::sqrt(637875.0 / 537251.0);
    const auto fs_135_2431_5 = std::sqrt(91125.0 / 5909761.0);
    const auto fs_135_2431_770 = std::sqrt(1275750.0 / 537251.0);
    const auto fs_135_2_5 = std::sqrt(91125.0 / 4.0);
    const auto fs_135_46189_105 = std::sqrt(1913625.0 / 2133423721.0);
    const auto fs_135_46189_385 = std::sqrt(637875.0 / 193947611.0);
    const auto fs_135_46189_5 = std::sqrt(91125.0 / 2133423721.0);
    const auto fs_135_46189_770 = std::sqrt(1275750.0 / 193947611.0);
    const auto fs_135_4_11 = std::sqrt(200475.0 / 16.0);
    const auto fs_135_4_14 = std::sqrt(127575.0 / 8.0);
    const auto fs_135_4_15 = std::sqrt(273375.0 / 16.0);
    const auto fs_135_4_21 = std::sqrt(382725.0 / 16.0);
    const auto fs_135_4_3 = std::sqrt(54675.0 / 16.0);
    const auto fs_135_4_33 = std::sqrt(601425.0 / 16.0);
    const auto fs_135_4_35 = std::sqrt(637875.0 / 16.0);
    const auto fs_135_4_5 = std::sqrt(91125.0 / 16.0);
    const auto fs_135_52_66 = std::sqrt(601425.0 / 1352.0);
    const auto fs_135_572_105 = std::sqrt(1913625.0 / 327184.0);
    const auto fs_135_572_385 = std::sqrt(637875.0 / 29744.0);
    const auto fs_135_572_5 = std::sqrt(91125.0 / 327184.0);
    const auto fs_135_572_770 = std::sqrt(637875.0 / 14872.0);
    const auto fs_135_8398_22 = std::sqrt(200475.0 / 35263202.0);
    const auto fs_135_8_110 = std::sqrt(1002375.0 / 32.0);
    const auto fs_135_8_2 = std::sqrt(18225.0 / 32.0);
    const auto fs_135_8_6 = std::sqrt(54675.0 / 32.0);
    const auto fs_140_143_2 = std::sqrt(39200.0 / 20449.0);
    const auto fs_140_143_21 = std::sqrt(411600.0 / 20449.0);
    const auto fs_140_143_7 = std::sqrt(137200.0 / 20449.0);
    const auto fs_1470_2431_7 = std::sqrt(15126300.0 / 5909761.0);
    const auto fs_1470_46189_7 = std::sqrt(15126300.0 / 2133423721.0);
    const auto fs_147_187_3 = std::sqrt(64827.0 / 34969.0);
    const auto fs_147_2431_715 = std::sqrt(108045.0 / 41327.0);
    const auto fs_147_4862_22 = std::sqrt(21609.0 / 1074502.0);
    const auto fs_147_884_78 = std::sqrt(64827.0 / 30056.0);
    const auto fs_147_9724_1430 = std::sqrt(108045.0 / 330616.0);
    const auto fs_14_187_3 = std::sqrt(588.0 / 34969.0);
    const auto fs_14_2431_715 = std::sqrt(980.0 / 41327.0);
    const auto fs_14_429_2 = std::sqrt(392.0 / 184041.0);
    const auto fs_14_429_21 = std::sqrt(1372.0 / 61347.0);
    const auto fs_14_429_7 = std::sqrt(1372.0 / 184041.0);
    const auto fs_150_143_3 = std::sqrt(67500.0 / 20449.0);
    const auto fs_150_2431_1001 = std::sqrt(157500.0 / 41327.0);
    const auto fs_150_2717_3 = std::sqrt(67500.0 / 7382089.0);
    const auto fs_150_46189_1001 = std::sqrt(157500.0 / 14919047.0);
    const auto fs_15_104_2730 = std::sqrt(23625.0 / 416.0);
    const auto fs_15_1144_10010 = std::sqrt(7875.0 / 4576.0);
    const auto fs_15_11_5 = std::sqrt(1125.0 / 121.0);
    const auto fs_15_13_11 = std::sqrt(2475.0 / 169.0);
    const auto fs_15_13_210 = std::sqrt(47250.0 / 169.0);
    const auto fs_15_13_462 = std::sqrt(103950.0 / 169.0);
    const auto fs_15_13_65 = std::sqrt(1125.0 / 13.0);
    const auto fs_15_13_7 = std::sqrt(1575.0 / 169.0);
    const auto fs_15_187_1001 = std::sqrt(20475.0 / 3179.0);
    const auto fs_15_1_5 = std::sqrt(1125.0);
    const auto fs_15_1_7 = std::sqrt(1575.0);
    const auto fs_15_221_770 = std::sqrt(173250.0 / 48841.0);
    const auto fs_15_22_11 = std::sqrt(225.0 / 44.0);
    const auto fs_15_22_14 = std::sqrt(1575.0 / 242.0);
    const auto fs_15_22_15 = std::sqrt(3375.0 / 484.0);
    const auto fs_15_22_21 = std::sqrt(4725.0 / 484.0);
    const auto fs_15_22_3 = std::sqrt(675.0 / 484.0);
    const auto fs_15_22_33 = std::sqrt(675.0 / 44.0);
    const auto fs_15_22_35 = std::sqrt(7875.0 / 484.0);
    const auto fs_15_22_5 = std::sqrt(1125.0 / 484.0);
    const auto fs_15_2431_1001 = std::sqrt(1575.0 / 41327.0);
    const auto fs_15_2431_11 = std::sqrt(225.0 / 537251.0);
    const auto fs_15_2431_154 = std::sqrt(3150.0 / 537251.0);
    const auto fs_15_2431_66 = std::sqrt(1350.0 / 537251.0);
    const auto fs_15_2431_715 = std::sqrt(1125.0 / 41327.0);
    const auto fs_15_2431_858 = std::sqrt(1350.0 / 41327.0);
    const auto fs_15_26_105 = std::sqrt(23625.0 / 676.0);
    const auto fs_15_26_14 = std::sqrt(1575.0 / 338.0);
    const auto fs_15_26_231 = std::sqrt(51975.0 / 676.0);
    const auto fs_15_26_30 = std::sqrt(3375.0 / 338.0);
    const auto fs_15_2_11 = std::sqrt(2475.0 / 4.0);
    const auto fs_15_2_14 = std::sqrt(1575.0 / 2.0);
    const auto fs_15_2_15 = std::sqrt(3375.0 / 4.0);
    const auto fs_15_2_21 = std::sqrt(4725.0 / 4.0);
    const auto fs_15_2_3 = std::sqrt(675.0 / 4.0);
    const auto fs_15_2_33 = std::sqrt(7425.0 / 4.0);
    const auto fs_15_2_35 = std::sqrt(7875.0 / 4.0);
    const auto fs_15_2_5 = std::sqrt(1125.0 / 4.0);
    const auto fs_15_3553_1001 = std::sqrt(20475.0 / 1147619.0);
    const auto fs_15_4199_2431 = std::sqrt(2475.0 / 79781.0);
    const auto fs_15_4199_4862 = std::sqrt(4950.0 / 79781.0);
    const auto fs_15_4199_770 = std::sqrt(173250.0 / 17631601.0);
    const auto fs_15_442_2730 = std::sqrt(23625.0 / 7514.0);
    const auto fs_15_44_1001 = std::sqrt(20475.0 / 176.0);
    const auto fs_15_44_110 = std::sqrt(1125.0 / 88.0);
    const auto fs_15_44_2 = std::sqrt(225.0 / 968.0);
    const auto fs_15_44_6 = std::sqrt(675.0 / 968.0);
    const auto fs_15_4862_10010 = std::sqrt(7875.0 / 82654.0);
    const auto fs_15_4_110 = std::sqrt(12375.0 / 8.0);
    const auto fs_15_4_154 = std::sqrt(17325.0 / 8.0);
    const auto fs_15_4_2 = std::sqrt(225.0 / 8.0);
    const auto fs_15_4_30 = std::sqrt(3375.0 / 8.0);
    const auto fs_15_4_55 = std::sqrt(12375.0 / 16.0);
    const auto fs_15_4_6 = std::sqrt(675.0 / 8.0);
    const auto fs_15_4_70 = std::sqrt(7875.0 / 8.0);
    const auto fs_15_52_770 = std::sqrt(86625.0 / 1352.0);
    const auto fs_15_8398_2730 = std::sqrt(23625.0 / 2712554.0);
    const auto fs_15_8_10 = std::sqrt(1125.0 / 32.0);
    const auto fs_15_8_210 = std::sqrt(23625.0 / 32.0);
    const auto fs_15_92378_10010 = std::sqrt(7875.0 / 29838094.0);
    const auto fs_165_26_5 = std::sqrt(136125.0 / 676.0);
    const auto fs_165_4199_3 = std::sqrt(81675.0 / 17631601.0);
    const auto fs_165_52_42 = std::sqrt(571725.0 / 1352.0);
    const auto fs_1665_2431_10 = std::sqrt(27722250.0 / 5909761.0);
    const auto fs_1665_46189_10 = std::sqrt(27722250.0 / 2133423721.0);
    const auto fs_1665_572_10 = std::sqrt(13861125.0 / 163592.0);
    const auto fs_16_13_3 = std::sqrt(768.0 / 169.0);
    const auto fs_1785_286_5 = std::sqrt(15931125.0 / 81796.0);
    const auto fs_1827_9724_30 = std::sqrt(50068935.0 / 47278088.0);
    const auto fs_1869_4862_15 = std::sqrt(52397415.0 / 23639044.0);
    const auto fs_1890_2431_3 = std::sqrt(10716300.0 / 5909761.0);
    const auto fs_1890_46189_3 = std::sqrt(10716300.0 / 2133423721.0);
    const auto fs_189_2431_35 = std::sqrt(1250235.0 / 5909761.0);
    const auto fs_189_442_7 = std::sqrt(250047.0 / 195364.0);
    const auto fs_189_4862_286 = std::sqrt(35721.0 / 82654.0);
    const auto fs_189_4862_462 = std::sqrt(750141.0 / 1074502.0);
    const auto fs_189_9724_4290 = std::sqrt(535815.0 / 330616.0);
    const auto fs_18_2431_35 = std::sqrt(11340.0 / 5909761.0);
    const auto fs_18_4199_143 = std::sqrt(3564.0 / 1356277.0);
    const auto fs_18_4199_3003 = std::sqrt(74844.0 / 1356277.0);
    const auto fs_198_4199_35 = std::sqrt(1372140.0 / 17631601.0);
    const auto fs_1_13_14 = std::sqrt(14.0 / 169.0);
    const auto fs_1_13_231 = std::sqrt(231.0 / 169.0);
    const auto fs_1_187_385 = std::sqrt(35.0 / 3179.0);
    const auto fs_1_221_10 = std::sqrt(10.0 / 48841.0);
    const auto fs_1_221_1326 = std::sqrt(6.0 / 221.0);
    const auto fs_1_221_182 = std::sqrt(14.0 / 3757.0);
    const auto fs_1_221_210 = std::sqrt(210.0 / 48841.0);
    const auto fs_1_221_221 = std::sqrt(1.0 / 221.0);
    const auto fs_1_221_455 = std::sqrt(35.0 / 3757.0);
    const auto fs_1_2431_165 = std::sqrt(15.0 / 537251.0);
    const auto fs_1_2431_2145 = std::sqrt(15.0 / 41327.0);
    const auto fs_1_2_6 = std::sqrt(3.0 / 2.0);
    const auto fs_1_4862_30030 = std::sqrt(105.0 / 82654.0);
    const auto fs_1_4862_6006 = std::sqrt(21.0 / 82654.0);
    const auto fs_1_51_6 = std::sqrt(2.0 / 867.0);
    const auto fs_20_1_3 = std::sqrt(1200.0);
    const auto fs_20_663_2 = std::sqrt(800.0 / 439569.0);
    const auto fs_210_11_2 = std::sqrt(88200.0 / 121.0);
    const auto fs_210_143_11 = std::sqrt(44100.0 / 1859.0);
    const auto fs_210_143_3 = std::sqrt(132300.0 / 20449.0);
    const auto fs_210_143_42 = std::sqrt(1852200.0 / 20449.0);
    const auto fs_210_143_5 = std::sqrt(220500.0 / 20449.0);
    const auto fs_210_143_7 = std::sqrt(308700.0 / 20449.0);
    const auto fs_210_2431_165 = std::sqrt(661500.0 / 537251.0);
    const auto fs_210_2431_21 = std::sqrt(926100.0 / 5909761.0);
    const auto fs_210_2431_35 = std::sqrt(1543500.0 / 5909761.0);
    const auto fs_210_2717_5 = std::sqrt(220500.0 / 7382089.0);
    const auto fs_210_46189_165 = std::sqrt(661500.0 / 193947611.0);
    const auto fs_210_46189_21 = std::sqrt(926100.0 / 2133423721.0);
    const auto fs_210_46189_35 = std::sqrt(1543500.0 / 2133423721.0);
    const auto fs_21_143_105 = std::sqrt(46305.0 / 20449.0);
    const auto fs_21_221_14 = std::sqrt(6174.0 / 48841.0);
    const auto fs_21_221_3 = std::sqrt(1323.0 / 48841.0);
    const auto fs_21_221_39 = std::sqrt(1323.0 / 3757.0);
    const auto fs_21_221_442 = std::sqrt(882.0 / 221.0);
    const auto fs_21_221_455 = std::sqrt(15435.0 / 3757.0);
    const auto fs_21_2431_165 = std::sqrt(6615.0 / 537251.0);
    const auto fs_21_2431_462 = std::sqrt(18522.0 / 537251.0);
    const auto fs_21_374_385 = std::sqrt(15435.0 / 12716.0);
    const auto fs_21_4199_2145 = std::sqrt(72765.0 / 1356277.0);
    const auto fs_21_4199_286 = std::sqrt(9702.0 / 1356277.0);
    const auto fs_21_4199_715 = std::sqrt(24255.0 / 1356277.0);
    const auto fs_21_442_1326 = std::sqrt(1323.0 / 442.0);
    const auto fs_21_442_182 = std::sqrt(3087.0 / 7514.0);
    const auto fs_21_442_221 = std::sqrt(441.0 / 884.0);
    const auto fs_21_442_455 = std::sqrt(15435.0 / 15028.0);
    const auto fs_21_4862_165 = std::sqrt(6615.0 / 2149004.0);
    const auto fs_21_4862_2145 = std::sqrt(6615.0 / 165308.0);
    const auto fs_21_8398_390 = std::sqrt(6615.0 / 2712554.0);
    const auto fs_21_9724_30030 = std::sqrt(46305.0 / 330616.0);
    const auto fs_21_9724_6006 = std::sqrt(9261.0 / 330616.0);
    const auto fs_2205_32_2 = std::sqrt(4862025.0 / 512.0);
    const auto fs_2247_4862_10 = std::sqrt(25245045.0 / 11819522.0);
    const auto fs_225_104_30 = std::sqrt(759375.0 / 5408.0);
    const auto fs_225_221_13 = std::sqrt(50625.0 / 3757.0);
    const auto fs_225_2431_385 = std::sqrt(1771875.0 / 537251.0);
    const auto fs_225_26_11 = std::sqrt(556875.0 / 676.0);
    const auto fs_225_26_7 = std::sqrt(354375.0 / 676.0);
    const auto fs_225_4199_13 = std::sqrt(50625.0 / 1356277.0);
    const auto fs_225_442_30 = std::sqrt(759375.0 / 97682.0);
    const auto fs_225_46189_385 = std::sqrt(1771875.0 / 193947611.0);
    const auto fs_225_52_13 = std::sqrt(50625.0 / 208.0);
    const auto fs_225_52_30 = std::sqrt(759375.0 / 1352.0);
    const auto fs_225_572_385 = std::sqrt(1771875.0 / 29744.0);
    const auto fs_225_8398_30 = std::sqrt(759375.0 / 35263202.0);
    const auto fs_22_663_5 = std::sqrt(2420.0 / 439569.0);
    const auto fs_231_442_2 = std::sqrt(53361.0 / 97682.0);
    const auto fs_240_221_3 = std::sqrt(172800.0 / 48841.0);
    const auto fs_240_2431_35 = std::sqrt(2016000.0 / 5909761.0);
    const auto fs_240_2431_385 = std::sqrt(2016000.0 / 537251.0);
    const auto fs_240_4199_3 = std::sqrt(172800.0 / 17631601.0);
    const auto fs_240_46189_35 = std::sqrt(2016000.0 / 2133423721.0);
    const auto fs_240_46189_385 = std::sqrt(2016000.0 / 193947611.0);
    const auto fs_245_143_2 = std::sqrt(120050.0 / 20449.0);
    const auto fs_245_4_2 = std::sqrt(60025.0 / 8.0);
    const auto fs_24_2431_154 = std::sqrt(8064.0 / 537251.0);
    const auto fs_24_2431_165 = std::sqrt(8640.0 / 537251.0);
    const auto fs_252_2431_154 = std::sqrt(889056.0 / 537251.0);
    const auto fs_252_2431_165 = std::sqrt(952560.0 / 537251.0);
    const auto fs_25_13_3 = std::sqrt(1875.0 / 169.0);
    const auto fs_25_2431_143 = std::sqrt(625.0 / 41327.0);
    const auto fs_25_2_2 = std::sqrt(625.0 / 2.0);
    const auto fs_25_4_42 = std::sqrt(13125.0 / 8.0);
    const auto fs_265_2431_3 = std::sqrt(210675.0 / 5909761.0);
    const auto fs_270_2431_2 = std::sqrt(145800.0 / 5909761.0);
    const auto fs_27_4199_2002 = std::sqrt(112266.0 / 1356277.0);
    const auto fs_27_4199_286 = std::sqrt(16038.0 / 1356277.0);
    const auto fs_280_143_2 = std::sqrt(156800.0 / 20449.0);
    const auto fs_2835_16_2 = std::sqrt(8037225.0 / 128.0);
    const auto fs_2835_16_3 = std::sqrt(24111675.0 / 256.0);
    const auto fs_2835_2431_2 = std::sqrt(16074450.0 / 5909761.0);
    const auto fs_2835_32_10 = std::sqrt(40186125.0 / 512.0);
    const auto fs_2835_32_2 = std::sqrt(8037225.0 / 512.0);
    const auto fs_2835_32_3 = std::sqrt(24111675.0 / 1024.0);
    const auto fs_2835_64_10 = std::sqrt(40186125.0 / 2048.0);
    const auto fs_28_429_2 = std::sqrt(1568.0 / 184041.0);
    const auto fs_29_2431_143 = std::sqrt(841.0 / 41327.0);
    const auto fs_2_13_210 = std::sqrt(840.0 / 169.0);
    const auto fs_2_13_462 = std::sqrt(1848.0 / 169.0);
    const auto fs_2_143_105 = std::sqrt(420.0 / 20449.0);
    const auto fs_2_221_14 = std::sqrt(56.0 / 48841.0);
    const auto fs_2_221_154 = std::sqrt(616.0 / 48841.0);
    const auto fs_2_221_3 = std::sqrt(12.0 / 48841.0);
    const auto fs_2_221_30 = std::sqrt(120.0 / 48841.0);
    const auto fs_2_221_39 = std::sqrt(12.0 / 3757.0);
    const auto fs_2_221_442 = std::sqrt(8.0 / 221.0);
    const auto fs_2_221_455 = std::sqrt(140.0 / 3757.0);
    const auto fs_2_221_55 = std::sqrt(220.0 / 48841.0);
    const auto fs_2_221_70 = std::sqrt(280.0 / 48841.0);
    const auto fs_2_2431_462 = std::sqrt(168.0 / 537251.0);
    const auto fs_2_663_14 = std::sqrt(56.0 / 439569.0);
    const auto fs_2_663_231 = std::sqrt(308.0 / 146523.0);
    const auto fs_300_143_33 = std::sqrt(270000.0 / 1859.0);
    const auto fs_3087_9724_22 = std::sqrt(9529569.0 / 4298008.0);
    const auto fs_30_11_2 = std::sqrt(1800.0 / 121.0);
    const auto fs_30_13_35 = std::sqrt(31500.0 / 169.0);
    const auto fs_30_143_2310 = std::sqrt(189000.0 / 1859.0);
    const auto fs_30_1_2 = std::sqrt(1800.0);
    const auto fs_30_221_105 = std::sqrt(94500.0 / 48841.0);
    const auto fs_30_4199_105 = std::sqrt(94500.0 / 17631601.0);
    const auto fs_30_4199_11 = std::sqrt(9900.0 / 17631601.0);
    const auto fs_3105_1144_6 = std::sqrt(28923075.0 / 654368.0);
    const auto fs_3105_4862_6 = std::sqrt(28923075.0 / 11819522.0);
    const auto fs_3105_92378_6 = std::sqrt(28923075.0 / 4266847442.0);
    const auto fs_315_143_30 = std::sqrt(2976750.0 / 20449.0);
    const auto fs_315_16_14 = std::sqrt(694575.0 / 128.0);
    const auto fs_315_16_15 = std::sqrt(1488375.0 / 256.0);
    const auto fs_315_16_21 = std::sqrt(2083725.0 / 256.0);
    const auto fs_315_16_3 = std::sqrt(297675.0 / 256.0);
    const auto fs_315_16_33 = std::sqrt(3274425.0 / 256.0);
    const auto fs_315_16_35 = std::sqrt(3472875.0 / 256.0);
    const auto fs_315_16_42 = std::sqrt(2083725.0 / 128.0);
    const auto fs_315_16_6 = std::sqrt(297675.0 / 128.0);
    const auto fs_315_22_11 = std::sqrt(99225.0 / 44.0);
    const auto fs_315_22_3 = std::sqrt(297675.0 / 484.0);
    const auto fs_315_22_5 = std::sqrt(496125.0 / 484.0);
    const auto fs_315_22_7 = std::sqrt(694575.0 / 484.0);
    const auto fs_315_2431_286 = std::sqrt(198450.0 / 41327.0);
    const auto fs_315_286_35 = std::sqrt(3472875.0 / 81796.0);
    const auto fs_315_286_70 = std::sqrt(3472875.0 / 40898.0);
    const auto fs_315_32_30 = std::sqrt(1488375.0 / 512.0);
    const auto fs_315_32_6 = std::sqrt(297675.0 / 512.0);
    const auto fs_315_32_66 = std::sqrt(3274425.0 / 512.0);
    const auto fs_315_44_10 = std::sqrt(496125.0 / 968.0);
    const auto fs_315_44_22 = std::sqrt(99225.0 / 88.0);
    const auto fs_315_46189_286 = std::sqrt(198450.0 / 14919047.0);
    const auto fs_315_4862_1001 = std::sqrt(694575.0 / 165308.0);
    const auto fs_315_4862_11 = std::sqrt(99225.0 / 2149004.0);
    const auto fs_315_4862_154 = std::sqrt(694575.0 / 1074502.0);
    const auto fs_315_4862_66 = std::sqrt(297675.0 / 1074502.0);
    const auto fs_315_4862_715 = std::sqrt(496125.0 / 165308.0);
    const auto fs_315_4862_858 = std::sqrt(297675.0 / 82654.0);
    const auto fs_315_4_2 = std::sqrt(99225.0 / 8.0);
    const auto fs_315_572_286 = std::sqrt(99225.0 / 1144.0);
    const auto fs_315_8_2 = std::sqrt(99225.0 / 32.0);
    const auto fs_315_8_21 = std::sqrt(2083725.0 / 64.0);
    const auto fs_315_8_7 = std::sqrt(694575.0 / 64.0);
    const auto fs_32_663_3 = std::sqrt(1024.0 / 146523.0);
    const auto fs_33_4199_5 = std::sqrt(5445.0 / 17631601.0);
    const auto fs_35_143_30 = std::sqrt(36750.0 / 20449.0);
    const auto fs_35_143_6 = std::sqrt(7350.0 / 20449.0);
    const auto fs_35_143_66 = std::sqrt(7350.0 / 1859.0);
    const auto fs_35_1_2 = std::sqrt(2450.0);
    const auto fs_35_1_21 = std::sqrt(25725.0);
    const auto fs_35_1_7 = std::sqrt(8575.0);
    const auto fs_35_2431_11 = std::sqrt(1225.0 / 537251.0);
    const auto fs_35_2_14 = std::sqrt(8575.0 / 2.0);
    const auto fs_35_2_15 = std::sqrt(18375.0 / 4.0);
    const auto fs_35_2_21 = std::sqrt(25725.0 / 4.0);
    const auto fs_35_2_3 = std::sqrt(3675.0 / 4.0);
    const auto fs_35_2_33 = std::sqrt(40425.0 / 4.0);
    const auto fs_35_2_35 = std::sqrt(42875.0 / 4.0);
    const auto fs_35_2_42 = std::sqrt(25725.0 / 2.0);
    const auto fs_35_2_6 = std::sqrt(3675.0 / 2.0);
    const auto fs_35_4_30 = std::sqrt(18375.0 / 8.0);
    const auto fs_35_4_6 = std::sqrt(3675.0 / 8.0);
    const auto fs_35_4_66 = std::sqrt(40425.0 / 8.0);
    const auto fs_367_4862_2 = std::sqrt(134689.0 / 11819522.0);
    const auto fs_36_4199_1001 = std::sqrt(99792.0 / 1356277.0);
    const auto fs_36_4199_182 = std::sqrt(18144.0 / 1356277.0);
    const auto fs_36_4199_330 = std::sqrt(427680.0 / 17631601.0);
    const auto fs_36_4199_442 = std::sqrt(2592.0 / 79781.0);
    const auto fs_375_1144_770 = std::sqrt(4921875.0 / 59488.0);
    const auto fs_375_26_3 = std::sqrt(421875.0 / 676.0);
    const auto fs_375_4862_770 = std::sqrt(4921875.0 / 1074502.0);
    const auto fs_375_92378_770 = std::sqrt(4921875.0 / 387895222.0);
    const auto fs_3_13_154 = std::sqrt(1386.0 / 169.0);
    const auto fs_3_13_30 = std::sqrt(270.0 / 169.0);
    const auto fs_3_13_55 = std::sqrt(495.0 / 169.0);
    const auto fs_3_13_70 = std::sqrt(630.0 / 169.0);
    const auto fs_3_143_11 = std::sqrt(9.0 / 1859.0);
    const auto fs_3_143_14 = std::sqrt(126.0 / 20449.0);
    const auto fs_3_143_15 = std::sqrt(135.0 / 20449.0);
    const auto fs_3_143_21 = std::sqrt(189.0 / 20449.0);
    const auto fs_3_143_3 = std::sqrt(27.0 / 20449.0);
    const auto fs_3_143_33 = std::sqrt(27.0 / 1859.0);
    const auto fs_3_143_35 = std::sqrt(315.0 / 20449.0);
    const auto fs_3_143_5 = std::sqrt(45.0 / 20449.0);
    const auto fs_3_221_221 = std::sqrt(9.0 / 221.0);
    const auto fs_3_221_66 = std::sqrt(594.0 / 48841.0);
    const auto fs_3_221_91 = std::sqrt(63.0 / 3757.0);
    const auto fs_3_2431_10010 = std::sqrt(630.0 / 41327.0);
    const auto fs_3_2431_12155 = std::sqrt(45.0 / 2431.0);
    const auto fs_3_26_10 = std::sqrt(45.0 / 338.0);
    const auto fs_3_26_210 = std::sqrt(945.0 / 338.0);
    const auto fs_3_286_110 = std::sqrt(45.0 / 3718.0);
    const auto fs_3_286_2 = std::sqrt(9.0 / 40898.0);
    const auto fs_3_286_6 = std::sqrt(27.0 / 40898.0);
    const auto fs_3_374_22 = std::sqrt(9.0 / 6358.0);
    const auto fs_3_4199_11 = std::sqrt(99.0 / 17631601.0);
    const auto fs_3_4199_138567 = std::sqrt(297.0 / 4199.0);
    const auto fs_3_4199_143 = std::sqrt(99.0 / 1356277.0);
    const auto fs_3_4199_146965 = std::sqrt(315.0 / 4199.0);
    const auto fs_3_4199_176358 = std::sqrt(378.0 / 4199.0);
    const auto fs_3_4199_25194 = std::sqrt(54.0 / 4199.0);
    const auto fs_3_4199_30030 = std::sqrt(20790.0 / 1356277.0);
    const auto fs_3_4199_323323 = std::sqrt(693.0 / 4199.0);
    const auto fs_3_4199_33 = std::sqrt(297.0 / 17631601.0);
    const auto fs_3_4199_46189 = std::sqrt(99.0 / 4199.0);
    const auto fs_3_4199_62985 = std::sqrt(135.0 / 4199.0);
    const auto fs_3_4199_910 = std::sqrt(630.0 / 1356277.0);
    const auto fs_3_4199_92378 = std::sqrt(198.0 / 4199.0);
    const auto fs_3_4199_9282 = std::sqrt(378.0 / 79781.0);
    const auto fs_3_4862_30030 = std::sqrt(945.0 / 82654.0);
    const auto fs_3_8398_182 = std::sqrt(63.0 / 2712554.0);
    const auto fs_3_8398_2 = std::sqrt(9.0 / 35263202.0);
    const auto fs_3_8398_2002 = std::sqrt(693.0 / 2712554.0);
    const auto fs_3_8398_26 = std::sqrt(9.0 / 2712554.0);
    const auto fs_3_8398_6006 = std::sqrt(2079.0 / 2712554.0);
    const auto fs_3_8398_910 = std::sqrt(315.0 / 2712554.0);
    const auto fs_405_1144_462 = std::sqrt(3444525.0 / 59488.0);
    const auto fs_405_286_30 = std::sqrt(2460375.0 / 40898.0);
    const auto fs_405_4862_462 = std::sqrt(3444525.0 / 1074502.0);
    const auto fs_405_4_2 = std::sqrt(164025.0 / 8.0);
    const auto fs_405_4_3 = std::sqrt(492075.0 / 16.0);
    const auto fs_405_8_10 = std::sqrt(820125.0 / 32.0);
    const auto fs_405_92378_462 = std::sqrt(3444525.0 / 387895222.0);
    const auto fs_41_2431_66 = std::sqrt(10086.0 / 537251.0);
    const auto fs_42_221_91 = std::sqrt(12348.0 / 3757.0);
    const auto fs_42_2431_143 = std::sqrt(1764.0 / 41327.0);
    const auto fs_42_4199_429 = std::sqrt(58212.0 / 1356277.0);
    const auto fs_4305_2431_2 = std::sqrt(37066050.0 / 5909761.0);
    const auto fs_4305_46189_2 = std::sqrt(37066050.0 / 2133423721.0);
    const auto fs_4305_572_2 = std::sqrt(18533025.0 / 163592.0);
    const auto fs_435_1144_462 = std::sqrt(3973725.0 / 59488.0);
    const auto fs_435_1144_70 = std::sqrt(6622875.0 / 654368.0);
    const auto fs_435_4862_462 = std::sqrt(3973725.0 / 1074502.0);
    const auto fs_435_4862_70 = std::sqrt(6622875.0 / 11819522.0);
    const auto fs_435_92378_462 = std::sqrt(3973725.0 / 387895222.0);
    const auto fs_435_92378_70 = std::sqrt(6622875.0 / 4266847442.0);
    const auto fs_441_2431_143 = std::sqrt(194481.0 / 41327.0);
    const auto fs_441_4862_165 = std::sqrt(2917215.0 / 2149004.0);
    const auto fs_45_187_70 = std::sqrt(141750.0 / 34969.0);
    const auto fs_45_22_2 = std::sqrt(2025.0 / 242.0);
    const auto fs_45_22_3 = std::sqrt(6075.0 / 484.0);
    const auto fs_45_2431_66 = std::sqrt(12150.0 / 537251.0);
    const auto fs_45_2431_77 = std::sqrt(14175.0 / 537251.0);
    const auto fs_45_26_154 = std::sqrt(155925.0 / 338.0);
    const auto fs_45_26_30 = std::sqrt(30375.0 / 338.0);
    const auto fs_45_26_55 = std::sqrt(111375.0 / 676.0);
    const auto fs_45_26_70 = std::sqrt(70875.0 / 338.0);
    const auto fs_45_286_2002 = std::sqrt(14175.0 / 286.0);
    const auto fs_45_2_2 = std::sqrt(2025.0 / 2.0);
    const auto fs_45_2_3 = std::sqrt(6075.0 / 4.0);
    const auto fs_45_3553_70 = std::sqrt(141750.0 / 12623809.0);
    const auto fs_45_4199_231 = std::sqrt(467775.0 / 17631601.0);
    const auto fs_45_4199_429 = std::sqrt(66825.0 / 1356277.0);
    const auto fs_45_44_10 = std::sqrt(10125.0 / 968.0);
    const auto fs_45_44_70 = std::sqrt(70875.0 / 968.0);
    const auto fs_45_46189_66 = std::sqrt(12150.0 / 193947611.0);
    const auto fs_45_46189_77 = std::sqrt(14175.0 / 193947611.0);
    const auto fs_45_4862_10 = std::sqrt(10125.0 / 11819522.0);
    const auto fs_45_4_10 = std::sqrt(10125.0 / 8.0);
    const auto fs_45_52_10 = std::sqrt(10125.0 / 1352.0);
    const auto fs_45_52_210 = std::sqrt(212625.0 / 1352.0);
    const auto fs_45_572_66 = std::sqrt(6075.0 / 14872.0);
    const auto fs_45_572_77 = std::sqrt(14175.0 / 29744.0);
    const auto fs_45_8398_2002 = std::sqrt(155925.0 / 2712554.0);
    const auto fs_45_8_66 = std::sqrt(66825.0 / 32.0);
    const auto fs_495_4199_2 = std::sqrt(490050.0 / 17631601.0);
    const auto fs_49_2431_33 = std::sqrt(7203.0 / 537251.0);
    const auto fs_49_858_2 = std::sqrt(2401.0 / 368082.0);
    const auto fs_4_13_35 = std::sqrt(560.0 / 169.0);
    const auto fs_4_221_91 = std::sqrt(112.0 / 3757.0);
    const auto fs_4_663_210 = std::sqrt(1120.0 / 146523.0);
    const auto fs_4_663_462 = std::sqrt(2464.0 / 146523.0);
    const auto fs_50_663_3 = std::sqrt(2500.0 / 146523.0);
    const auto fs_525_2431_154 = std::sqrt(3858750.0 / 537251.0);
    const auto fs_525_286_11 = std::sqrt(275625.0 / 7436.0);
    const auto fs_525_286_14 = std::sqrt(1929375.0 / 40898.0);
    const auto fs_525_286_35 = std::sqrt(9646875.0 / 81796.0);
    const auto fs_525_286_66 = std::sqrt(826875.0 / 3718.0);
    const auto fs_525_46189_154 = std::sqrt(3858750.0 / 193947611.0);
    const auto fs_525_4862_143 = std::sqrt(275625.0 / 165308.0);
    const auto fs_525_572_154 = std::sqrt(1929375.0 / 14872.0);
    const auto fs_53_2431_33 = std::sqrt(8427.0 / 537251.0);
    const auto fs_5505_1144_2 = std::sqrt(30305025.0 / 654368.0);
    const auto fs_5505_4862_2 = std::sqrt(30305025.0 / 11819522.0);
    const auto fs_5505_92378_2 = std::sqrt(30305025.0 / 4266847442.0);
    const auto fs_555_2431_55 = std::sqrt(1540125.0 / 537251.0);
    const auto fs_555_46189_55 = std::sqrt(1540125.0 / 193947611.0);
    const auto fs_555_572_55 = std::sqrt(1540125.0 / 29744.0);
    const auto fs_5565_4862_3 = std::sqrt(92907675.0 / 23639044.0);
    const auto fs_55_4_5 = std::sqrt(15125.0 / 16.0);
    const auto fs_55_8_42 = std::sqrt(63525.0 / 32.0);
    const auto fs_57_4862_154 = std::sqrt(22743.0 / 1074502.0);
    const auto fs_5_13_42 = std::sqrt(1050.0 / 169.0);
    const auto fs_5_1_35 = std::sqrt(875.0);
    const auto fs_5_221_30 = std::sqrt(750.0 / 48841.0);
    const auto fs_5_221_91 = std::sqrt(175.0 / 3757.0);
    const auto fs_5_2431_7293 = std::sqrt(75.0 / 2431.0);
    const auto fs_5_2431_858 = std::sqrt(150.0 / 41327.0);
    const auto fs_5_2_210 = std::sqrt(2625.0 / 2.0);
    const auto fs_5_2_462 = std::sqrt(5775.0 / 2.0);
    const auto fs_5_4_14 = std::sqrt(175.0 / 8.0);
    const auto fs_5_4_231 = std::sqrt(5775.0 / 16.0);
    const auto fs_609_4862_143 = std::sqrt(370881.0 / 165308.0);
    const auto fs_60_13_3 = std::sqrt(10800.0 / 169.0);
    const auto fs_60_143_35 = std::sqrt(126000.0 / 20449.0);
    const auto fs_60_143_385 = std::sqrt(126000.0 / 1859.0);
    const auto fs_60_221_65 = std::sqrt(18000.0 / 3757.0);
    const auto fs_60_4199_286 = std::sqrt(79200.0 / 1356277.0);
    const auto fs_60_4199_33 = std::sqrt(118800.0 / 17631601.0);
    const auto fs_60_4199_65 = std::sqrt(18000.0 / 1356277.0);
    const auto fs_630_2431_35 = std::sqrt(13891500.0 / 5909761.0);
    const auto fs_630_2431_70 = std::sqrt(27783000.0 / 5909761.0);
    const auto fs_630_46189_35 = std::sqrt(13891500.0 / 2133423721.0);
    const auto fs_630_46189_70 = std::sqrt(27783000.0 / 2133423721.0);
    const auto fs_63_143_5 = std::sqrt(19845.0 / 20449.0);
    const auto fs_63_2431_2431 = std::sqrt(3969.0 / 2431.0);
    const auto fs_63_2431_715 = std::sqrt(19845.0 / 41327.0);
    const auto fs_63_4199_143 = std::sqrt(43659.0 / 1356277.0);
    const auto fs_63_442_221 = std::sqrt(3969.0 / 884.0);
    const auto fs_63_442_91 = std::sqrt(27783.0 / 15028.0);
    const auto fs_63_4862_10010 = std::sqrt(138915.0 / 82654.0);
    const auto fs_63_4862_12155 = std::sqrt(19845.0 / 9724.0);
    const auto fs_63_748_22 = std::sqrt(3969.0 / 25432.0);
    const auto fs_63_9724_30030 = std::sqrt(416745.0 / 330616.0);
    const auto fs_65_8_6 = std::sqrt(12675.0 / 32.0);
    const auto fs_6_143_5 = std::sqrt(180.0 / 20449.0);
    const auto fs_6_2431_2431 = std::sqrt(36.0 / 2431.0);
    const auto fs_6_2431_715 = std::sqrt(180.0 / 41327.0);
    const auto fs_6_4199_15 = std::sqrt(540.0 / 17631601.0);
    const auto fs_6_4199_15015 = std::sqrt(41580.0 / 1356277.0);
    const auto fs_6_4199_1547 = std::sqrt(252.0 / 79781.0);
    const auto fs_6_4199_2002 = std::sqrt(5544.0 / 1356277.0);
    const auto fs_6_4199_2431 = std::sqrt(396.0 / 79781.0);
    const auto fs_6_4199_24310 = std::sqrt(3960.0 / 79781.0);
    const auto fs_6_4199_41990 = std::sqrt(360.0 / 4199.0);
    const auto fs_6_4199_7735 = std::sqrt(1260.0 / 79781.0);
    const auto fs_70_143_14 = std::sqrt(68600.0 / 20449.0);
    const auto fs_70_143_15 = std::sqrt(73500.0 / 20449.0);
    const auto fs_70_143_21 = std::sqrt(102900.0 / 20449.0);
    const auto fs_70_143_3 = std::sqrt(14700.0 / 20449.0);
    const auto fs_70_143_33 = std::sqrt(14700.0 / 1859.0);
    const auto fs_70_143_35 = std::sqrt(171500.0 / 20449.0);
    const auto fs_70_143_42 = std::sqrt(205800.0 / 20449.0);
    const auto fs_70_143_6 = std::sqrt(29400.0 / 20449.0);
    const auto fs_70_1_2 = std::sqrt(9800.0);
    const auto fs_735_2431_35 = std::sqrt(18907875.0 / 5909761.0);
    const auto fs_735_286_7 = std::sqrt(3781575.0 / 81796.0);
    const auto fs_735_44_2 = std::sqrt(540225.0 / 968.0);
    const auto fs_735_46189_35 = std::sqrt(18907875.0 / 2133423721.0);
    const auto fs_735_4862_11 = std::sqrt(540225.0 / 2149004.0);
    const auto fs_735_572_35 = std::sqrt(18907875.0 / 327184.0);
    const auto fs_75_104_210 = std::sqrt(590625.0 / 5408.0);
    const auto fs_75_13_2 = std::sqrt(11250.0 / 169.0);
    const auto fs_75_187_5 = std::sqrt(28125.0 / 34969.0);
    const auto fs_75_221_10 = std::sqrt(56250.0 / 48841.0);
    const auto fs_75_221_6 = std::sqrt(33750.0 / 48841.0);
    const auto fs_75_2431_231 = std::sqrt(118125.0 / 537251.0);
    const auto fs_75_26_42 = std::sqrt(118125.0 / 338.0);
    const auto fs_75_286_1001 = std::sqrt(39375.0 / 572.0);
    const auto fs_75_3553_5 = std::sqrt(28125.0 / 12623809.0);
    const auto fs_75_4199_10 = std::sqrt(56250.0 / 17631601.0);
    const auto fs_75_4199_6 = std::sqrt(33750.0 / 17631601.0);
    const auto fs_75_442_210 = std::sqrt(590625.0 / 97682.0);
    const auto fs_75_44_5 = std::sqrt(28125.0 / 1936.0);
    const auto fs_75_46189_231 = std::sqrt(118125.0 / 193947611.0);
    const auto fs_75_4_11 = std::sqrt(61875.0 / 16.0);
    const auto fs_75_4_7 = std::sqrt(39375.0 / 16.0);
    const auto fs_75_52_10 = std::sqrt(28125.0 / 1352.0);
    const auto fs_75_52_6 = std::sqrt(16875.0 / 1352.0);
    const auto fs_75_572_231 = std::sqrt(118125.0 / 29744.0);
    const auto fs_75_8398_210 = std::sqrt(590625.0 / 35263202.0);
    const auto fs_75_8_30 = std::sqrt(84375.0 / 32.0);
    const auto fs_7707_9724_2 = std::sqrt(59397849.0 / 47278088.0);
    const auto fs_7_143_11 = std::sqrt(49.0 / 1859.0);
    const auto fs_7_143_3 = std::sqrt(147.0 / 20449.0);
    const auto fs_7_143_5 = std::sqrt(245.0 / 20449.0);
    const auto fs_7_143_7 = std::sqrt(343.0 / 20449.0);
    const auto fs_7_286_10 = std::sqrt(245.0 / 40898.0);
    const auto fs_7_286_22 = std::sqrt(49.0 / 3718.0);
    const auto fs_7_429_14 = std::sqrt(686.0 / 184041.0);
    const auto fs_7_429_15 = std::sqrt(245.0 / 61347.0);
    const auto fs_7_429_21 = std::sqrt(343.0 / 61347.0);
    const auto fs_7_429_3 = std::sqrt(49.0 / 61347.0);
    const auto fs_7_429_33 = std::sqrt(49.0 / 5577.0);
    const auto fs_7_429_35 = std::sqrt(1715.0 / 184041.0);
    const auto fs_7_429_42 = std::sqrt(686.0 / 61347.0);
    const auto fs_7_429_6 = std::sqrt(98.0 / 61347.0);
    const auto fs_7_442_78 = std::sqrt(147.0 / 7514.0);
    const auto fs_7_4862_1430 = std::sqrt(245.0 / 82654.0);
    const auto fs_7_858_30 = std::sqrt(245.0 / 122694.0);
    const auto fs_7_858_6 = std::sqrt(49.0 / 122694.0);
    const auto fs_7_858_66 = std::sqrt(49.0 / 11154.0);
    const auto fs_810_2431_30 = std::sqrt(19683000.0 / 5909761.0);
    const auto fs_810_46189_30 = std::sqrt(19683000.0 / 2133423721.0);
    const auto fs_840_2431_42 = std::sqrt(29635200.0 / 5909761.0);
    const auto fs_840_46189_42 = std::sqrt(29635200.0 / 2133423721.0);
    const auto fs_84_2431_1430 = std::sqrt(70560.0 / 41327.0);
    const auto fs_84_4199_143 = std::sqrt(77616.0 / 1356277.0);
    const auto fs_861_4862_66 = std::sqrt(2223963.0 / 1074502.0);
    const auto fs_87_4862_30 = std::sqrt(113535.0 / 11819522.0);
    const auto fs_89_2431_15 = std::sqrt(118815.0 / 5909761.0);
    const auto fs_8_221_5 = std::sqrt(320.0 / 48841.0);
    const auto fs_8_221_7 = std::sqrt(448.0 / 48841.0);
    const auto fs_8_2431_1430 = std::sqrt(640.0 / 41327.0);
    const auto fs_8_663_35 = std::sqrt(2240.0 / 439569.0);
    const auto fs_90_13_5 = std::sqrt(40500.0 / 169.0);
    const auto fs_90_13_7 = std::sqrt(56700.0 / 169.0);
    const auto fs_90_2431_2002 = std::sqrt(113400.0 / 41327.0);
    const auto fs_90_46189_2002 = std::sqrt(113400.0 / 14919047.0);
    const auto fs_945_16_11 = std::sqrt(9823275.0 / 256.0);
    const auto fs_945_16_14 = std::sqrt(6251175.0 / 128.0);
    const auto fs_945_16_15 = std::sqrt(13395375.0 / 256.0);
    const auto fs_945_16_21 = std::sqrt(18753525.0 / 256.0);
    const auto fs_945_16_3 = std::sqrt(2679075.0 / 256.0);
    const auto fs_945_16_33 = std::sqrt(29469825.0 / 256.0);
    const auto fs_945_16_35 = std::sqrt(31255875.0 / 256.0);
    const auto fs_945_16_5 = std::sqrt(4465125.0 / 256.0);
    const auto fs_945_16_7 = std::sqrt(6251175.0 / 256.0);
    const auto fs_945_2431_10 = std::sqrt(8930250.0 / 5909761.0);
    const auto fs_945_286_3 = std::sqrt(2679075.0 / 81796.0);
    const auto fs_945_32_10 = std::sqrt(4465125.0 / 512.0);
    const auto fs_945_32_11 = std::sqrt(9823275.0 / 1024.0);
    const auto fs_945_32_110 = std::sqrt(49116375.0 / 512.0);
    const auto fs_945_32_14 = std::sqrt(6251175.0 / 512.0);
    const auto fs_945_32_15 = std::sqrt(13395375.0 / 1024.0);
    const auto fs_945_32_2 = std::sqrt(893025.0 / 512.0);
    const auto fs_945_32_21 = std::sqrt(18753525.0 / 1024.0);
    const auto fs_945_32_22 = std::sqrt(9823275.0 / 512.0);
    const auto fs_945_32_3 = std::sqrt(2679075.0 / 1024.0);
    const auto fs_945_32_33 = std::sqrt(29469825.0 / 1024.0);
    const auto fs_945_32_35 = std::sqrt(31255875.0 / 1024.0);
    const auto fs_945_32_5 = std::sqrt(4465125.0 / 1024.0);
    const auto fs_945_32_6 = std::sqrt(2679075.0 / 512.0);
    const auto fs_945_46189_10 = std::sqrt(8930250.0 / 2133423721.0);
    const auto fs_945_4_2 = std::sqrt(893025.0 / 8.0);
    const auto fs_945_572_10 = std::sqrt(4465125.0 / 163592.0);
    const auto fs_945_64_110 = std::sqrt(49116375.0 / 2048.0);
    const auto fs_945_64_2 = std::sqrt(893025.0 / 2048.0);
    const auto fs_945_64_6 = std::sqrt(2679075.0 / 2048.0);
    const auto fs_945_8_2 = std::sqrt(893025.0 / 32.0);
    const auto fs_945_8_5 = std::sqrt(4465125.0 / 64.0);
    const auto fs_945_9724_10 = std::sqrt(4465125.0 / 47278088.0);
    const auto fs_9_143_2 = std::sqrt(162.0 / 20449.0);
    const auto fs_9_143_3 = std::sqrt(243.0 / 20449.0);
    const auto fs_9_221_7 = std::sqrt(567.0 / 48841.0);
    const auto fs_9_2431_286 = std::sqrt(162.0 / 41327.0);
    const auto fs_9_2431_462 = std::sqrt(3402.0 / 537251.0);
    const auto fs_9_26_66 = std::sqrt(2673.0 / 338.0);
    const auto fs_9_286_10 = std::sqrt(405.0 / 40898.0);
    const auto fs_9_4199_12597 = std::sqrt(243.0 / 4199.0);
    const auto fs_9_4199_14586 = std::sqrt(5346.0 / 79781.0);
    const auto fs_9_4199_165 = std::sqrt(13365.0 / 17631601.0);
    const auto fs_9_4199_2431 = std::sqrt(891.0 / 79781.0);
    const auto fs_9_4199_39 = std::sqrt(243.0 / 1356277.0);
    const auto fs_9_4862_4290 = std::sqrt(1215.0 / 82654.0);
    const auto fs_9_8398_1430 = std::sqrt(4455.0 / 2712554.0);
    const auto fs_9_8398_30030 = std::sqrt(93555.0 / 2712554.0);

    // NOTE: the angular components are formed in 124 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph1_p1, ph3_p1, ph5_p1, ph7_p1, ph9_p1, ph11_p1, ph11_p11, ab_2, pc_0 : simd::cache_line_size())
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

        const auto h1_p1 = ph1_p1[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h9_p1 = ph9_p1[k];
        const auto h11_p1 = ph11_p1[k];
        const auto h11_p11 = ph11_p11[k];

        pc_0[k] = e_0 * fs_945_32_33 * h1_p1 + e_1 * fs_945_32_22 * h3_p1 - e_1 * fs_945_16_33 * r_2 * h1_p1 + e_2 * fs_15_4_55 * h5_p1 - e_2 * fs_105_4_22 * r_2 * h3_p1 + e_2 * fs_135_4_33 * r_4 * h1_p1 + e_3 * fs_75_572_231 * h7_p1 - e_3 * fs_45_26_55 * r_2 * h5_p1 + e_3 * fs_315_44_22 * r_4 * h3_p1 - e_3 * fs_15_2_33 * r_6 * h1_p1 + e_4 * fs_21_4862_165 * h9_p1 - e_4 * fs_75_2431_231 * r_2 * h7_p1 + e_4 * fs_3_13_55 * r_4 * h5_p1 - e_4 * fs_105_143_22 * r_6 * h3_p1 + e_4 * fs_15_22_33 * r_8 * h1_p1 + e_5 * fs_3_8398_2 * h11_p1 + e_5 * fs_3_4199_323323 * h11_p11 - e_5 * fs_1_2431_165 * r_2 * h9_p1 + e_5 * fs_75_46189_231 * r_4 * h7_p1 - e_5 * fs_2_221_55 * r_6 * h5_p1 + e_5 * fs_7_286_22 * r_8 * h3_p1 - e_5 * fs_3_143_33 * r_10 * h1_p1;
    }

    // NOTE: the angular components are formed in 124 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph1_0, ph3_0, ph5_0, ph7_0, ph9_0, ph11_0, ph11_p10, ab_2, pc_1 : simd::cache_line_size())
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

        const auto h1_0 = ph1_0[k];
        const auto h3_0 = ph3_0[k];
        const auto h5_0 = ph5_0[k];
        const auto h7_0 = ph7_0[k];
        const auto h9_0 = ph9_0[k];
        const auto h11_0 = ph11_0[k];
        const auto h11_p10 = ph11_p10[k];

        pc_1[k] = e_0 * fs_945_32_11 * h1_0 + e_1 * fs_945_16_11 * h3_0 - e_1 * fs_945_16_11 * r_2 * h1_0 + e_2 * fs_75_4_11 * h5_0 - e_2 * fs_105_2_11 * r_2 * h3_0 + e_2 * fs_135_4_11 * r_4 * h1_0 + e_3 * fs_525_286_11 * h7_0 - e_3 * fs_225_26_11 * r_2 * h5_0 + e_3 * fs_315_22_11 * r_4 * h3_0 - e_3 * fs_15_2_11 * r_6 * h1_0 + e_4 * fs_315_4862_11 * h9_0 - e_4 * fs_1050_2431_11 * r_2 * h7_0 + e_4 * fs_15_13_11 * r_4 * h5_0 - e_4 * fs_210_143_11 * r_6 * h3_0 + e_4 * fs_15_22_11 * r_8 * h1_0 + e_5 * fs_3_4199_11 * h11_0 + e_5 * fs_3_4199_176358 * h11_p10 - e_5 * fs_15_2431_11 * r_2 * h9_0 + e_5 * fs_1050_46189_11 * r_4 * h7_0 - e_5 * fs_10_221_11 * r_6 * h5_0 + e_5 * fs_7_143_11 * r_8 * h3_0 - e_5 * fs_3_143_11 * r_10 * h1_0;
    }

    // NOTE: the angular components are formed in 124 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph1_p1, ph3_p1, ph5_p1, ph7_p1, ph9_p1, ph9_p9, ph11_p1, ph11_p9, ab_2, pc_2 : simd::cache_line_size())
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

        const auto h1_p1 = ph1_p1[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h9_p1 = ph9_p1[k];
        const auto h9_p9 = ph9_p9[k];
        const auto h11_p1 = ph11_p1[k];
        const auto h11_p9 = ph11_p9[k];

        pc_2[k] = - e_0 * fs_945_64_2 * h1_p1 - e_1 * fs_945_16_3 * h3_p1 + e_1 * fs_945_32_2 * r_2 * h1_p1 - e_2 * fs_75_8_30 * h5_p1 + e_2 * fs_105_2_3 * r_2 * h3_p1 - e_2 * fs_135_8_2 * r_4 * h1_p1 - e_3 * fs_525_286_14 * h7_p1 + e_3 * fs_225_52_30 * r_2 * h5_p1 - e_3 * fs_315_22_3 * r_4 * h3_p1 + e_3 * fs_15_4_2 * r_6 * h1_p1 - e_4 * fs_945_9724_10 * h9_p1 + e_4 * fs_63_4862_12155 * h9_p9 + e_4 * fs_1050_2431_14 * r_2 * h7_p1 - e_4 * fs_15_26_30 * r_4 * h5_p1 + e_4 * fs_210_143_3 * r_6 * h3_p1 - e_4 * fs_15_44_2 * r_8 * h1_p1 - e_5 * fs_3_4199_33 * h11_p1 + e_5 * fs_3_4199_92378 * h11_p9 + e_5 * fs_45_4862_10 * r_2 * h9_p1 - e_5 * fs_3_2431_12155 * r_2 * h9_p9 - e_5 * fs_1050_46189_14 * r_4 * h7_p1 + e_5 * fs_5_221_30 * r_6 * h5_p1 - e_5 * fs_7_143_3 * r_8 * h3_p1 + e_5 * fs_3_286_2 * r_10 * h1_p1;
    }

    // NOTE: the angular components are formed in 124 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph3_p2, ph5_p2, ph7_p2, ph9_p2, ph9_p8, ph11_p2, ph11_p8, ab_2, pc_3 : simd::cache_line_size())
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

        const auto h3_p2 = ph3_p2[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h9_p2 = ph9_p2[k];
        const auto h9_p8 = ph9_p8[k];
        const auto h11_p2 = ph11_p2[k];
        const auto h11_p8 = ph11_p8[k];

        pc_3[k] = e_1 * f_945_16 * h3_p2 + e_2 * fs_75_4_7 * h5_p2 - e_2 * f_105_2 * r_2 * h3_p2 + e_3 * fs_315_286_70 * h7_p2 - e_3 * fs_225_26_7 * r_2 * h5_p2 + e_3 * f_315_22 * r_4 * h3_p2 + e_4 * fs_315_4862_66 * h9_p2 + e_4 * fs_105_4862_7293 * h9_p8 - e_4 * fs_630_2431_70 * r_2 * h7_p2 + e_4 * fs_15_13_7 * r_4 * h5_p2 - e_4 * f_210_143 * r_6 * h3_p2 + e_5 * fs_3_4199_143 * h11_p2 + e_5 * fs_3_4199_46189 * h11_p8 - e_5 * fs_15_2431_66 * r_2 * h9_p2 - e_5 * fs_5_2431_7293 * r_2 * h9_p8 + e_5 * fs_630_46189_70 * r_4 * h7_p2 - e_5 * fs_10_221_7 * r_6 * h5_p2 + e_5 * f_7_143 * r_8 * h3_p2;
    }

    // NOTE: the angular components are formed in 124 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph3_p3, ph5_p3, ph7_p3, ph7_p7, ph9_p3, ph9_p7, ph11_p3, ph11_p7, ab_2, pc_4 : simd::cache_line_size())
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

        const auto h3_p3 = ph3_p3[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h7_p3 = ph7_p3[k];
        const auto h7_p7 = ph7_p7[k];
        const auto h9_p3 = ph9_p3[k];
        const auto h9_p7 = ph9_p7[k];
        const auto h11_p3 = ph11_p3[k];
        const auto h11_p7 = ph11_p7[k];

        pc_4[k] = - e_1 * fs_315_32_6 * h3_p3 - e_2 * fs_25_4_42 * h5_p3 + e_2 * fs_35_4_6 * r_2 * h3_p3 - e_3 * fs_525_286_35 * h7_p3 + e_3 * fs_105_572_715 * h7_p7 + e_3 * fs_75_26_42 * r_2 * h5_p3 - e_3 * fs_105_44_6 * r_4 * h3_p3 - e_4 * fs_315_4862_154 * h9_p3 + e_4 * fs_315_4862_858 * h9_p7 + e_4 * fs_1050_2431_35 * r_2 * h7_p3 - e_4 * fs_105_2431_715 * r_2 * h7_p7 - e_4 * fs_5_13_42 * r_4 * h5_p3 + e_4 * fs_35_143_6 * r_6 * h3_p3 - e_5 * fs_3_8398_2002 * h11_p3 + e_5 * fs_9_4199_2431 * h11_p7 + e_5 * fs_15_2431_154 * r_2 * h9_p3 - e_5 * fs_15_2431_858 * r_2 * h9_p7 - e_5 * fs_1050_46189_35 * r_4 * h7_p3 + e_5 * fs_105_46189_715 * r_4 * h7_p7 + e_5 * fs_10_663_42 * r_6 * h5_p3 - e_5 * fs_7_858_6 * r_8 * h3_p3;
    }

    // NOTE: the angular components are formed in 124 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_2, pe_3, pe_4, pe_5, ph5_m5, ph5_p4, ph7_m5, ph7_p4, ph7_p6, ph9_m5, ph9_p4, ph9_p6, ph11_m5, ph11_p4, ph11_p6, ab_2, pc_5, pc_6 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h5_m5 = ph5_m5[k];
        const auto h5_p4 = ph5_p4[k];
        const auto h7_m5 = ph7_m5[k];
        const auto h7_p4 = ph7_p4[k];
        const auto h7_p6 = ph7_p6[k];
        const auto h9_m5 = ph9_m5[k];
        const auto h9_p4 = ph9_p4[k];
        const auto h9_p6 = ph9_p6[k];
        const auto h11_m5 = ph11_m5[k];
        const auto h11_p4 = ph11_p4[k];
        const auto h11_p6 = ph11_p6[k];

        pc_5[k] = e_2 * fs_15_8_210 * h5_p4 + e_3 * fs_525_572_154 * h7_p4 + e_3 * fs_75_286_1001 * h7_p6 - e_3 * fs_45_52_210 * r_2 * h5_p4 + e_4 * fs_63_9724_30030 * h9_p4 + e_4 * fs_315_4862_715 * h9_p6 - e_4 * fs_525_2431_154 * r_2 * h7_p4 - e_4 * fs_150_2431_1001 * r_2 * h7_p6 + e_4 * fs_3_26_210 * r_4 * h5_p4 + e_5 * fs_3_8398_6006 * h11_p4 + e_5 * fs_6_4199_2431 * h11_p6 - e_5 * fs_3_4862_30030 * r_2 * h9_p4 - e_5 * fs_15_2431_715 * r_2 * h9_p6 + e_5 * fs_525_46189_154 * r_4 * h7_p4 + e_5 * fs_150_46189_1001 * r_4 * h7_p6 - e_5 * fs_1_221_210 * r_6 * h5_p4;

        pc_6[k] = - e_2 * f_75_4 * h5_m5 - e_3 * fs_525_286_66 * h7_m5 + e_3 * f_225_26 * r_2 * h5_m5 - e_4 * fs_315_4862_1001 * h9_m5 + e_4 * fs_1050_2431_66 * r_2 * h7_m5 - e_4 * f_15_13 * r_4 * h5_m5 - e_5 * fs_6_4199_2002 * h11_m5 + e_5 * fs_15_2431_1001 * r_2 * h9_m5 - e_5 * fs_1050_46189_66 * r_4 * h7_m5 + e_5 * f_10_221 * r_6 * h5_m5;
    }

    // NOTE: the angular components are formed in 124 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_2, pe_3, pe_4, pe_5, ph5_m4, ph7_m6, ph7_m4, ph9_m6, ph9_m4, ph11_m6, ph11_m4, ab_2, pc_7 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h5_m4 = ph5_m4[k];
        const auto h7_m6 = ph7_m6[k];
        const auto h7_m4 = ph7_m4[k];
        const auto h9_m6 = ph9_m6[k];
        const auto h9_m4 = ph9_m4[k];
        const auto h11_m6 = ph11_m6[k];
        const auto h11_m4 = ph11_m4[k];

        pc_7[k] = e_2 * fs_15_8_210 * h5_m4 - e_3 * fs_75_286_1001 * h7_m6 + e_3 * fs_525_572_154 * h7_m4 - e_3 * fs_45_52_210 * r_2 * h5_m4 - e_4 * fs_315_4862_715 * h9_m6 + e_4 * fs_63_9724_30030 * h9_m4 + e_4 * fs_150_2431_1001 * r_2 * h7_m6 - e_4 * fs_525_2431_154 * r_2 * h7_m4 + e_4 * fs_3_26_210 * r_4 * h5_m4 - e_5 * fs_6_4199_2431 * h11_m6 + e_5 * fs_3_8398_6006 * h11_m4 + e_5 * fs_15_2431_715 * r_2 * h9_m6 - e_5 * fs_3_4862_30030 * r_2 * h9_m4 - e_5 * fs_150_46189_1001 * r_4 * h7_m6 + e_5 * fs_525_46189_154 * r_4 * h7_m4 - e_5 * fs_1_221_210 * r_6 * h5_m4;
    }

    // NOTE: the angular components are formed in 124 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph3_m3, ph5_m3, ph7_m7, ph7_m3, ph9_m7, ph9_m3, ph11_m7, ph11_m3, ab_2, pc_8 : simd::cache_line_size())
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

        const auto h3_m3 = ph3_m3[k];
        const auto h5_m3 = ph5_m3[k];
        const auto h7_m7 = ph7_m7[k];
        const auto h7_m3 = ph7_m3[k];
        const auto h9_m7 = ph9_m7[k];
        const auto h9_m3 = ph9_m3[k];
        const auto h11_m7 = ph11_m7[k];
        const auto h11_m3 = ph11_m3[k];

        pc_8[k] = - e_1 * fs_315_32_6 * h3_m3 - e_2 * fs_25_4_42 * h5_m3 + e_2 * fs_35_4_6 * r_2 * h3_m3 - e_3 * fs_105_572_715 * h7_m7 - e_3 * fs_525_286_35 * h7_m3 + e_3 * fs_75_26_42 * r_2 * h5_m3 - e_3 * fs_105_44_6 * r_4 * h3_m3 - e_4 * fs_315_4862_858 * h9_m7 - e_4 * fs_315_4862_154 * h9_m3 + e_4 * fs_105_2431_715 * r_2 * h7_m7 + e_4 * fs_1050_2431_35 * r_2 * h7_m3 - e_4 * fs_5_13_42 * r_4 * h5_m3 + e_4 * fs_35_143_6 * r_6 * h3_m3 - e_5 * fs_9_4199_2431 * h11_m7 - e_5 * fs_3_8398_2002 * h11_m3 + e_5 * fs_15_2431_858 * r_2 * h9_m7 + e_5 * fs_15_2431_154 * r_2 * h9_m3 - e_5 * fs_105_46189_715 * r_4 * h7_m7 - e_5 * fs_1050_46189_35 * r_4 * h7_m3 + e_5 * fs_10_663_42 * r_6 * h5_m3 - e_5 * fs_7_858_6 * r_8 * h3_m3;
    }

    // NOTE: the angular components are formed in 124 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph3_m2, ph5_m2, ph7_m2, ph9_m8, ph9_m2, ph11_m8, ph11_m2, ab_2, pc_9 : simd::cache_line_size())
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

        const auto h3_m2 = ph3_m2[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h9_m8 = ph9_m8[k];
        const auto h9_m2 = ph9_m2[k];
        const auto h11_m8 = ph11_m8[k];
        const auto h11_m2 = ph11_m2[k];

        pc_9[k] = e_1 * f_945_16 * h3_m2 + e_2 * fs_75_4_7 * h5_m2 - e_2 * f_105_2 * r_2 * h3_m2 + e_3 * fs_315_286_70 * h7_m2 - e_3 * fs_225_26_7 * r_2 * h5_m2 + e_3 * f_315_22 * r_4 * h3_m2 - e_4 * fs_105_4862_7293 * h9_m8 + e_4 * fs_315_4862_66 * h9_m2 - e_4 * fs_630_2431_70 * r_2 * h7_m2 + e_4 * fs_15_13_7 * r_4 * h5_m2 - e_4 * f_210_143 * r_6 * h3_m2 - e_5 * fs_3_4199_46189 * h11_m8 + e_5 * fs_3_4199_143 * h11_m2 + e_5 * fs_5_2431_7293 * r_2 * h9_m8 - e_5 * fs_15_2431_66 * r_2 * h9_m2 + e_5 * fs_630_46189_70 * r_4 * h7_m2 - e_5 * fs_10_221_7 * r_6 * h5_m2 + e_5 * f_7_143 * r_8 * h3_m2;
    }

    // NOTE: the angular components are formed in 124 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph1_m1, ph3_m1, ph5_m1, ph7_m1, ph9_m9, ph9_m1, ph11_m11, ph11_m10, ph11_m9, ph11_m1, ab_2, pc_10, pc_11, pc_12 : simd::cache_line_size())
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

        const auto h1_m1 = ph1_m1[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h7_m1 = ph7_m1[k];
        const auto h9_m9 = ph9_m9[k];
        const auto h9_m1 = ph9_m1[k];
        const auto h11_m11 = ph11_m11[k];
        const auto h11_m10 = ph11_m10[k];
        const auto h11_m9 = ph11_m9[k];
        const auto h11_m1 = ph11_m1[k];

        pc_10[k] = - e_0 * fs_945_64_2 * h1_m1 - e_1 * fs_945_16_3 * h3_m1 + e_1 * fs_945_32_2 * r_2 * h1_m1 - e_2 * fs_75_8_30 * h5_m1 + e_2 * fs_105_2_3 * r_2 * h3_m1 - e_2 * fs_135_8_2 * r_4 * h1_m1 - e_3 * fs_525_286_14 * h7_m1 + e_3 * fs_225_52_30 * r_2 * h5_m1 - e_3 * fs_315_22_3 * r_4 * h3_m1 + e_3 * fs_15_4_2 * r_6 * h1_m1 - e_4 * fs_63_4862_12155 * h9_m9 - e_4 * fs_945_9724_10 * h9_m1 + e_4 * fs_1050_2431_14 * r_2 * h7_m1 - e_4 * fs_15_26_30 * r_4 * h5_m1 + e_4 * fs_210_143_3 * r_6 * h3_m1 - e_4 * fs_15_44_2 * r_8 * h1_m1 - e_5 * fs_3_4199_92378 * h11_m9 - e_5 * fs_3_4199_33 * h11_m1 + e_5 * fs_3_2431_12155 * r_2 * h9_m9 + e_5 * fs_45_4862_10 * r_2 * h9_m1 - e_5 * fs_1050_46189_14 * r_4 * h7_m1 + e_5 * fs_5_221_30 * r_6 * h5_m1 - e_5 * fs_7_143_3 * r_8 * h3_m1 + e_5 * fs_3_286_2 * r_10 * h1_m1;

        pc_11[k] = - e_5 * fs_3_4199_176358 * h11_m10;

        pc_12[k] = - e_0 * fs_945_32_33 * h1_m1 - e_1 * fs_945_32_22 * h3_m1 + e_1 * fs_945_16_33 * r_2 * h1_m1 - e_2 * fs_15_4_55 * h5_m1 + e_2 * fs_105_4_22 * r_2 * h3_m1 - e_2 * fs_135_4_33 * r_4 * h1_m1 - e_3 * fs_75_572_231 * h7_m1 + e_3 * fs_45_26_55 * r_2 * h5_m1 - e_3 * fs_315_44_22 * r_4 * h3_m1 + e_3 * fs_15_2_33 * r_6 * h1_m1 - e_4 * fs_21_4862_165 * h9_m1 + e_4 * fs_75_2431_231 * r_2 * h7_m1 - e_4 * fs_3_13_55 * r_4 * h5_m1 + e_4 * fs_105_143_22 * r_6 * h3_m1 - e_4 * fs_15_22_33 * r_8 * h1_m1 - e_5 * fs_3_4199_323323 * h11_m11 - e_5 * fs_3_8398_2 * h11_m1 + e_5 * fs_1_2431_165 * r_2 * h9_m1 - e_5 * fs_75_46189_231 * r_4 * h7_m1 + e_5 * fs_2_221_55 * r_6 * h5_m1 - e_5 * fs_7_286_22 * r_8 * h3_m1 + e_5 * fs_3_143_33 * r_10 * h1_m1;
    }

    // NOTE: the angular components are formed in 124 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph3_p2, ph5_p2, ph7_p2, ph9_p2, ph11_p2, ph11_p10, ab_2, pc_13 : simd::cache_line_size())
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

        const auto h3_p2 = ph3_p2[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h9_p2 = ph9_p2[k];
        const auto h11_p2 = ph11_p2[k];
        const auto h11_p10 = ph11_p10[k];

        pc_13[k] = - e_1 * fs_945_32_22 * h3_p2 - e_2 * fs_15_4_154 * h5_p2 + e_2 * fs_105_4_22 * r_2 * h3_p2 - e_3 * fs_135_572_385 * h7_p2 + e_3 * fs_45_26_154 * r_2 * h5_p2 - e_3 * fs_315_44_22 * r_4 * h3_p2 - e_4 * fs_21_221_3 * h9_p2 + e_4 * fs_135_2431_385 * r_2 * h7_p2 - e_4 * fs_3_13_154 * r_4 * h5_p2 + e_4 * fs_105_143_22 * r_6 * h3_p2 - e_5 * fs_3_8398_26 * h11_p2 + e_5 * fs_3_4199_146965 * h11_p10 + e_5 * fs_2_221_3 * r_2 * h9_p2 - e_5 * fs_135_46189_385 * r_4 * h7_p2 + e_5 * fs_2_221_154 * r_6 * h5_p2 - e_5 * fs_7_286_22 * r_8 * h3_p2;
    }

    // NOTE: the angular components are formed in 124 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph1_p1, ph5_p1, ph7_p1, ph9_p1, ph9_p9, ph11_p1, ph11_p9, ab_2, pc_14 : simd::cache_line_size())
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

        const auto h1_p1 = ph1_p1[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h9_p1 = ph9_p1[k];
        const auto h9_p9 = ph9_p9[k];
        const auto h11_p1 = ph11_p1[k];
        const auto h11_p9 = ph11_p9[k];

        pc_14[k] = e_0 * fs_945_64_110 * h1_p1 - e_1 * fs_945_32_110 * r_2 * h1_p1 - e_2 * fs_45_8_66 * h5_p1 + e_2 * fs_135_8_110 * r_4 * h1_p1 - e_3 * fs_15_52_770 * h7_p1 + e_3 * fs_135_52_66 * r_2 * h5_p1 - e_3 * fs_15_4_110 * r_6 * h1_p1 - e_4 * fs_63_748_22 * h9_p1 - e_4 * fs_63_442_221 * h9_p9 + e_4 * fs_15_221_770 * r_2 * h7_p1 - e_4 * fs_9_26_66 * r_4 * h5_p1 + e_4 * fs_15_44_110 * r_8 * h1_p1 - e_5 * fs_6_4199_15 * h11_p1 + e_5 * fs_6_4199_41990 * h11_p9 + e_5 * fs_3_374_22 * r_2 * h9_p1 + e_5 * fs_3_221_221 * r_2 * h9_p9 - e_5 * fs_15_4199_770 * r_4 * h7_p1 + e_5 * fs_3_221_66 * r_6 * h5_p1 - e_5 * fs_3_286_110 * r_10 * h1_p1;
    }

    // NOTE: the angular components are formed in 124 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph1_0, ph3_0, ph5_0, ph7_0, ph9_0, ph9_p8, ph11_0, ph11_p8, ab_2, pc_15 : simd::cache_line_size())
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

        const auto h1_0 = ph1_0[k];
        const auto h3_0 = ph3_0[k];
        const auto h5_0 = ph5_0[k];
        const auto h7_0 = ph7_0[k];
        const auto h9_0 = ph9_0[k];
        const auto h9_p8 = ph9_p8[k];
        const auto h11_0 = ph11_0[k];
        const auto h11_p8 = ph11_p8[k];

        pc_15[k] = e_0 * fs_945_16_5 * h1_0 + e_1 * fs_945_16_5 * h3_0 - e_1 * fs_945_8_5 * r_2 * h1_0 - e_2 * fs_15_1_5 * h5_0 - e_2 * fs_105_2_5 * r_2 * h3_0 + e_2 * fs_135_2_5 * r_4 * h1_0 - e_3 * fs_1785_286_5 * h7_0 + e_3 * fs_90_13_5 * r_2 * h5_0 + e_3 * fs_315_22_5 * r_4 * h3_0 - e_3 * fs_15_1_5 * r_6 * h1_0 - e_4 * fs_63_143_5 * h9_0 - e_4 * fs_63_2431_2431 * h9_p8 + e_4 * fs_210_143_5 * r_2 * h7_0 - e_4 * fs_12_13_5 * r_4 * h5_0 - e_4 * fs_210_143_5 * r_6 * h3_0 + e_4 * fs_15_11_5 * r_8 * h1_0 - e_5 * fs_33_4199_5 * h11_0 + e_5 * fs_3_4199_138567 * h11_p8 + e_5 * fs_6_143_5 * r_2 * h9_0 + e_5 * fs_6_2431_2431 * r_2 * h9_p8 - e_5 * fs_210_2717_5 * r_4 * h7_0 + e_5 * fs_8_221_5 * r_6 * h5_0 + e_5 * fs_7_143_5 * r_8 * h3_0 - e_5 * fs_6_143_5 * r_10 * h1_0;
    }

    // NOTE: the angular components are formed in 124 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph1_p1, ph3_p1, ph5_p1, ph7_p1, ph7_p7, ph9_p1, ph9_p7, ph11_p1, ph11_p7, ab_2, pc_16 : simd::cache_line_size())
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

        const auto h1_p1 = ph1_p1[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h7_p7 = ph7_p7[k];
        const auto h9_p1 = ph9_p1[k];
        const auto h9_p7 = ph9_p7[k];
        const auto h11_p1 = ph11_p1[k];
        const auto h11_p7 = ph11_p7[k];

        pc_16[k] = - e_0 * fs_945_64_6 * h1_p1 - e_1 * f_945_8 * h3_p1 + e_1 * fs_945_32_6 * r_2 * h1_p1 - e_2 * fs_15_8_10 * h5_p1 + e_2 * f_105_1 * r_2 * h3_p1 - e_2 * fs_135_8_6 * r_4 * h1_p1 + e_3 * fs_210_143_42 * h7_p1 - e_3 * fs_315_572_286 * h7_p7 + e_3 * fs_45_52_10 * r_2 * h5_p1 - e_3 * f_315_11 * r_4 * h3_p1 + e_3 * fs_15_4_6 * r_6 * h1_p1 + e_4 * fs_1827_9724_30 * h9_p1 - e_4 * fs_21_4862_2145 * h9_p7 - e_4 * fs_840_2431_42 * r_2 * h7_p1 + e_4 * fs_315_2431_286 * r_2 * h7_p7 - e_4 * fs_3_26_10 * r_4 * h5_p1 + e_4 * f_420_143 * r_6 * h3_p1 - e_4 * fs_15_44_6 * r_8 * h1_p1 + e_5 * fs_30_4199_11 * h11_p1 + e_5 * fs_6_4199_24310 * h11_p7 - e_5 * fs_87_4862_30 * r_2 * h9_p1 + e_5 * fs_1_2431_2145 * r_2 * h9_p7 + e_5 * fs_840_46189_42 * r_4 * h7_p1 - e_5 * fs_315_46189_286 * r_4 * h7_p7 + e_5 * fs_1_221_10 * r_6 * h5_p1 - e_5 * f_14_143 * r_8 * h3_p1 + e_5 * fs_3_286_6 * r_10 * h1_p1;
    }

    // NOTE: the angular components are formed in 124 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph3_p2, ph5_p2, ph7_p2, ph7_p6, ph9_p2, ph9_p6, ph11_p2, ph11_p6, ab_2, pc_17 : simd::cache_line_size())
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

        const auto h3_p2 = ph3_p2[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h7_p6 = ph7_p6[k];
        const auto h9_p2 = ph9_p2[k];
        const auto h9_p6 = ph9_p6[k];
        const auto h11_p2 = ph11_p2[k];
        const auto h11_p6 = ph11_p6[k];

        pc_17[k] = e_1 * fs_945_32_10 * h3_p2 + e_2 * fs_15_4_70 * h5_p2 - e_2 * fs_105_4_10 * r_2 * h3_p2 - e_3 * fs_735_286_7 * h7_p2 - e_3 * fs_15_44_1001 * h7_p6 - e_3 * fs_45_26_70 * r_2 * h5_p2 + e_3 * fs_315_44_10 * r_4 * h3_p2 - e_4 * fs_252_2431_165 * h9_p2 + e_4 * fs_63_2431_715 * h9_p6 + e_4 * fs_1470_2431_7 * r_2 * h7_p2 + e_4 * fs_15_187_1001 * r_2 * h7_p6 + e_4 * fs_3_13_70 * r_4 * h5_p2 - e_4 * fs_105_143_10 * r_6 * h3_p2 - e_5 * fs_9_8398_1430 * h11_p2 + e_5 * fs_15_4199_2431 * h11_p6 + e_5 * fs_24_2431_165 * r_2 * h9_p2 - e_5 * fs_6_2431_715 * r_2 * h9_p6 - e_5 * fs_1470_46189_7 * r_4 * h7_p2 - e_5 * fs_15_3553_1001 * r_4 * h7_p6 - e_5 * fs_2_221_70 * r_6 * h5_p2 + e_5 * fs_7_286_10 * r_8 * h3_p2;
    }

    // NOTE: the angular components are formed in 124 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph3_p3, ph5_p3, ph5_p5, ph7_p3, ph7_p5, ph9_p3, ph9_p5, ph11_p3, ph11_p5, ab_2, pc_18 : simd::cache_line_size())
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

        const auto h3_p3 = ph3_p3[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h5_p5 = ph5_p5[k];
        const auto h7_p3 = ph7_p3[k];
        const auto h7_p5 = ph7_p5[k];
        const auto h9_p3 = ph9_p3[k];
        const auto h9_p5 = ph9_p5[k];
        const auto h11_p3 = ph11_p3[k];
        const auto h11_p5 = ph11_p5[k];

        pc_18[k] = - e_1 * fs_315_16_6 * h3_p3 - e_2 * fs_55_8_42 * h5_p3 - e_2 * fs_15_8_210 * h5_p5 + e_2 * fs_35_2_6 * r_2 * h3_p3 + e_3 * fs_105_286_35 * h7_p3 - e_3 * fs_60_143_385 * h7_p5 + e_3 * fs_165_52_42 * r_2 * h5_p3 + e_3 * fs_45_52_210 * r_2 * h5_p5 - e_3 * fs_105_22_6 * r_4 * h3_p3 + e_4 * fs_1197_9724_154 * h9_p3 + e_4 * fs_189_9724_4290 * h9_p5 - e_4 * fs_210_2431_35 * r_2 * h7_p3 + e_4 * fs_240_2431_385 * r_2 * h7_p5 - e_4 * fs_11_26_42 * r_4 * h5_p3 - e_4 * fs_3_26_210 * r_4 * h5_p5 + e_4 * fs_70_143_6 * r_6 * h3_p3 + e_5 * fs_6_4199_2002 * h11_p3 + e_5 * fs_12_4199_2145 * h11_p5 - e_5 * fs_57_4862_154 * r_2 * h9_p3 - e_5 * fs_9_4862_4290 * r_2 * h9_p5 + e_5 * fs_210_46189_35 * r_4 * h7_p3 - e_5 * fs_240_46189_385 * r_4 * h7_p5 + e_5 * fs_11_663_42 * r_6 * h5_p3 + e_5 * fs_1_221_210 * r_6 * h5_p5 - e_5 * fs_7_429_6 * r_8 * h3_p3;
    }

    // NOTE: the angular components are formed in 124 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_2, pe_3, pe_4, pe_5, ph5_m4, ph7_m4, ph9_m4, ph11_m4, ab_2, pc_19 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h5_m4 = ph5_m4[k];
        const auto h7_m4 = ph7_m4[k];
        const auto h9_m4 = ph9_m4[k];
        const auto h11_m4 = ph11_m4[k];

        pc_19[k] = e_2 * f_60_1 * h5_m4 + e_3 * fs_105_286_165 * h7_m4 - e_3 * f_360_13 * r_2 * h5_m4 - e_4 * fs_441_2431_143 * h9_m4 - e_4 * fs_210_2431_165 * r_2 * h7_m4 + e_4 * f_48_13 * r_4 * h5_m4 - e_5 * fs_21_4199_715 * h11_m4 + e_5 * fs_42_2431_143 * r_2 * h9_m4 + e_5 * fs_210_46189_165 * r_4 * h7_m4 - e_5 * f_32_221 * r_6 * h5_m4;
    }

    // NOTE: the angular components are formed in 124 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph3_m3, ph5_m5, ph5_m3, ph7_m5, ph7_m3, ph9_m5, ph9_m3, ph11_m5, ph11_m3, ab_2, pc_20 : simd::cache_line_size())
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

        const auto h3_m3 = ph3_m3[k];
        const auto h5_m5 = ph5_m5[k];
        const auto h5_m3 = ph5_m3[k];
        const auto h7_m5 = ph7_m5[k];
        const auto h7_m3 = ph7_m3[k];
        const auto h9_m5 = ph9_m5[k];
        const auto h9_m3 = ph9_m3[k];
        const auto h11_m5 = ph11_m5[k];
        const auto h11_m3 = ph11_m3[k];

        pc_20[k] = - e_1 * fs_315_16_6 * h3_m3 + e_2 * fs_15_8_210 * h5_m5 - e_2 * fs_55_8_42 * h5_m3 + e_2 * fs_35_2_6 * r_2 * h3_m3 + e_3 * fs_60_143_385 * h7_m5 + e_3 * fs_105_286_35 * h7_m3 - e_3 * fs_45_52_210 * r_2 * h5_m5 + e_3 * fs_165_52_42 * r_2 * h5_m3 - e_3 * fs_105_22_6 * r_4 * h3_m3 - e_4 * fs_189_9724_4290 * h9_m5 + e_4 * fs_1197_9724_154 * h9_m3 - e_4 * fs_240_2431_385 * r_2 * h7_m5 - e_4 * fs_210_2431_35 * r_2 * h7_m3 + e_4 * fs_3_26_210 * r_4 * h5_m5 - e_4 * fs_11_26_42 * r_4 * h5_m3 + e_4 * fs_70_143_6 * r_6 * h3_m3 - e_5 * fs_12_4199_2145 * h11_m5 + e_5 * fs_6_4199_2002 * h11_m3 + e_5 * fs_9_4862_4290 * r_2 * h9_m5 - e_5 * fs_57_4862_154 * r_2 * h9_m3 + e_5 * fs_240_46189_385 * r_4 * h7_m5 + e_5 * fs_210_46189_35 * r_4 * h7_m3 - e_5 * fs_1_221_210 * r_6 * h5_m5 + e_5 * fs_11_663_42 * r_6 * h5_m3 - e_5 * fs_7_429_6 * r_8 * h3_m3;
    }

    // NOTE: the angular components are formed in 124 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph3_m2, ph5_m2, ph7_m6, ph7_m2, ph9_m6, ph9_m2, ph11_m6, ph11_m2, ab_2, pc_21 : simd::cache_line_size())
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

        const auto h3_m2 = ph3_m2[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h7_m6 = ph7_m6[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h9_m6 = ph9_m6[k];
        const auto h9_m2 = ph9_m2[k];
        const auto h11_m6 = ph11_m6[k];
        const auto h11_m2 = ph11_m2[k];

        pc_21[k] = e_1 * fs_945_32_10 * h3_m2 + e_2 * fs_15_4_70 * h5_m2 - e_2 * fs_105_4_10 * r_2 * h3_m2 + e_3 * fs_15_44_1001 * h7_m6 - e_3 * fs_735_286_7 * h7_m2 - e_3 * fs_45_26_70 * r_2 * h5_m2 + e_3 * fs_315_44_10 * r_4 * h3_m2 - e_4 * fs_63_2431_715 * h9_m6 - e_4 * fs_252_2431_165 * h9_m2 - e_4 * fs_15_187_1001 * r_2 * h7_m6 + e_4 * fs_1470_2431_7 * r_2 * h7_m2 + e_4 * fs_3_13_70 * r_4 * h5_m2 - e_4 * fs_105_143_10 * r_6 * h3_m2 - e_5 * fs_15_4199_2431 * h11_m6 - e_5 * fs_9_8398_1430 * h11_m2 + e_5 * fs_6_2431_715 * r_2 * h9_m6 + e_5 * fs_24_2431_165 * r_2 * h9_m2 + e_5 * fs_15_3553_1001 * r_4 * h7_m6 - e_5 * fs_1470_46189_7 * r_4 * h7_m2 - e_5 * fs_2_221_70 * r_6 * h5_m2 + e_5 * fs_7_286_10 * r_8 * h3_m2;
    }

    // NOTE: the angular components are formed in 124 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph1_m1, ph3_m1, ph5_m1, ph7_m7, ph7_m1, ph9_m8, ph9_m7, ph9_m1, ph11_m8, ph11_m7, ph11_m1, ab_2, pc_22, pc_23 : simd::cache_line_size())
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

        const auto h1_m1 = ph1_m1[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h7_m7 = ph7_m7[k];
        const auto h7_m1 = ph7_m1[k];
        const auto h9_m8 = ph9_m8[k];
        const auto h9_m7 = ph9_m7[k];
        const auto h9_m1 = ph9_m1[k];
        const auto h11_m8 = ph11_m8[k];
        const auto h11_m7 = ph11_m7[k];
        const auto h11_m1 = ph11_m1[k];

        pc_22[k] = - e_0 * fs_945_64_6 * h1_m1 - e_1 * f_945_8 * h3_m1 + e_1 * fs_945_32_6 * r_2 * h1_m1 - e_2 * fs_15_8_10 * h5_m1 + e_2 * f_105_1 * r_2 * h3_m1 - e_2 * fs_135_8_6 * r_4 * h1_m1 + e_3 * fs_315_572_286 * h7_m7 + e_3 * fs_210_143_42 * h7_m1 + e_3 * fs_45_52_10 * r_2 * h5_m1 - e_3 * f_315_11 * r_4 * h3_m1 + e_3 * fs_15_4_6 * r_6 * h1_m1 + e_4 * fs_21_4862_2145 * h9_m7 + e_4 * fs_1827_9724_30 * h9_m1 - e_4 * fs_315_2431_286 * r_2 * h7_m7 - e_4 * fs_840_2431_42 * r_2 * h7_m1 - e_4 * fs_3_26_10 * r_4 * h5_m1 + e_4 * f_420_143 * r_6 * h3_m1 - e_4 * fs_15_44_6 * r_8 * h1_m1 - e_5 * fs_6_4199_24310 * h11_m7 + e_5 * fs_30_4199_11 * h11_m1 - e_5 * fs_1_2431_2145 * r_2 * h9_m7 - e_5 * fs_87_4862_30 * r_2 * h9_m1 + e_5 * fs_315_46189_286 * r_4 * h7_m7 + e_5 * fs_840_46189_42 * r_4 * h7_m1 + e_5 * fs_1_221_10 * r_6 * h5_m1 - e_5 * f_14_143 * r_8 * h3_m1 + e_5 * fs_3_286_6 * r_10 * h1_m1;

        pc_23[k] = e_4 * fs_63_2431_2431 * h9_m8 - e_5 * fs_3_4199_138567 * h11_m8 - e_5 * fs_6_2431_2431 * r_2 * h9_m8;
    }

    // NOTE: the angular components are formed in 124 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph1_m1, ph5_m1, ph7_m1, ph9_m9, ph9_m1, ph11_m9, ph11_m1, ab_2, pc_24 : simd::cache_line_size())
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

        const auto h1_m1 = ph1_m1[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h7_m1 = ph7_m1[k];
        const auto h9_m9 = ph9_m9[k];
        const auto h9_m1 = ph9_m1[k];
        const auto h11_m9 = ph11_m9[k];
        const auto h11_m1 = ph11_m1[k];

        pc_24[k] = - e_0 * fs_945_64_110 * h1_m1 + e_1 * fs_945_32_110 * r_2 * h1_m1 + e_2 * fs_45_8_66 * h5_m1 - e_2 * fs_135_8_110 * r_4 * h1_m1 + e_3 * fs_15_52_770 * h7_m1 - e_3 * fs_135_52_66 * r_2 * h5_m1 + e_3 * fs_15_4_110 * r_6 * h1_m1 + e_4 * fs_63_442_221 * h9_m9 + e_4 * fs_63_748_22 * h9_m1 - e_4 * fs_15_221_770 * r_2 * h7_m1 + e_4 * fs_9_26_66 * r_4 * h5_m1 - e_4 * fs_15_44_110 * r_8 * h1_m1 - e_5 * fs_6_4199_41990 * h11_m9 + e_5 * fs_6_4199_15 * h11_m1 - e_5 * fs_3_221_221 * r_2 * h9_m9 - e_5 * fs_3_374_22 * r_2 * h9_m1 + e_5 * fs_15_4199_770 * r_4 * h7_m1 - e_5 * fs_3_221_66 * r_6 * h5_m1 + e_5 * fs_3_286_110 * r_10 * h1_m1;
    }

    // NOTE: the angular components are formed in 124 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph3_m2, ph5_m2, ph7_m2, ph9_m2, ph11_m10, ph11_m2, ab_2, pc_25 : simd::cache_line_size())
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

        const auto h3_m2 = ph3_m2[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h9_m2 = ph9_m2[k];
        const auto h11_m10 = ph11_m10[k];
        const auto h11_m2 = ph11_m2[k];

        pc_25[k] = e_1 * fs_945_32_22 * h3_m2 + e_2 * fs_15_4_154 * h5_m2 - e_2 * fs_105_4_22 * r_2 * h3_m2 + e_3 * fs_135_572_385 * h7_m2 - e_3 * fs_45_26_154 * r_2 * h5_m2 + e_3 * fs_315_44_22 * r_4 * h3_m2 + e_4 * fs_21_221_3 * h9_m2 - e_4 * fs_135_2431_385 * r_2 * h7_m2 + e_4 * fs_3_13_154 * r_4 * h5_m2 - e_4 * fs_105_143_22 * r_6 * h3_m2 - e_5 * fs_3_4199_146965 * h11_m10 + e_5 * fs_3_8398_26 * h11_m2 - e_5 * fs_2_221_3 * r_2 * h9_m2 + e_5 * fs_135_46189_385 * r_4 * h7_m2 - e_5 * fs_2_221_154 * r_6 * h5_m2 + e_5 * fs_7_286_22 * r_8 * h3_m2;
    }

    // NOTE: the angular components are formed in 124 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph3_p3, ph5_p3, ph7_p3, ph9_p3, ph9_p9, ph11_p3, ph11_p9, ab_2, pc_26 : simd::cache_line_size())
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

        const auto h3_p3 = ph3_p3[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h7_p3 = ph7_p3[k];
        const auto h9_p3 = ph9_p3[k];
        const auto h9_p9 = ph9_p9[k];
        const auto h11_p3 = ph11_p3[k];
        const auto h11_p9 = ph11_p9[k];

        pc_26[k] = e_1 * fs_315_32_66 * h3_p3 + e_2 * fs_5_2_462 * h5_p3 - e_2 * fs_35_4_66 * r_2 * h3_p3 + e_3 * fs_225_572_385 * h7_p3 - e_3 * fs_15_13_462 * r_2 * h5_p3 + e_3 * fs_105_44_66 * r_4 * h3_p3 + e_4 * fs_21_221_14 * h9_p3 + e_4 * fs_21_442_1326 * h9_p9 - e_4 * fs_225_2431_385 * r_2 * h7_p3 + e_4 * fs_2_13_462 * r_4 * h5_p3 - e_4 * fs_35_143_66 * r_6 * h3_p3 + e_5 * fs_3_8398_182 * h11_p3 + e_5 * fs_3_4199_62985 * h11_p9 - e_5 * fs_2_221_14 * r_2 * h9_p3 - e_5 * fs_1_221_1326 * r_2 * h9_p9 + e_5 * fs_225_46189_385 * r_4 * h7_p3 - e_5 * fs_4_663_462 * r_6 * h5_p3 + e_5 * fs_7_858_66 * r_8 * h3_p3;
    }

    // NOTE: the angular components are formed in 124 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph3_p2, ph5_p2, ph7_p2, ph9_p2, ph9_p8, ph11_p2, ph11_p8, ab_2, pc_27 : simd::cache_line_size())
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

        const auto h3_p2 = ph3_p2[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h9_p2 = ph9_p2[k];
        const auto h9_p8 = ph9_p8[k];
        const auto h11_p2 = ph11_p2[k];
        const auto h11_p8 = ph11_p8[k];

        pc_27[k] = - e_1 * fs_315_16_33 * h3_p2 + e_2 * fs_5_4_231 * h5_p2 + e_2 * fs_35_2_33 * r_2 * h3_p2 + e_3 * fs_30_143_2310 * h7_p2 - e_3 * fs_15_26_231 * r_2 * h5_p2 - e_3 * fs_105_22_33 * r_4 * h3_p2 + e_4 * fs_231_442_2 * h9_p2 - e_4 * fs_21_442_221 * h9_p8 - e_4 * fs_120_2431_2310 * r_2 * h7_p2 + e_4 * fs_1_13_231 * r_4 * h5_p2 + e_4 * fs_70_143_33 * r_6 * h3_p2 + e_5 * fs_9_4199_39 * h11_p2 + e_5 * fs_9_4199_12597 * h11_p8 - e_5 * fs_11_221_2 * r_2 * h9_p2 + e_5 * fs_1_221_221 * r_2 * h9_p8 + e_5 * fs_120_46189_2310 * r_4 * h7_p2 - e_5 * fs_2_663_231 * r_6 * h5_p2 - e_5 * fs_7_429_33 * r_8 * h3_p2;
    }

    // NOTE: the angular components are formed in 124 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph1_p1, ph3_p1, ph5_p1, ph7_p1, ph7_p7, ph9_p1, ph9_p7, ph11_p1, ph11_p7, ab_2, pc_28 : simd::cache_line_size())
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

        const auto h1_p1 = ph1_p1[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h7_p7 = ph7_p7[k];
        const auto h9_p1 = ph9_p1[k];
        const auto h9_p7 = ph9_p7[k];
        const auto h11_p1 = ph11_p1[k];
        const auto h11_p7 = ph11_p7[k];

        pc_28[k] = e_0 * fs_2835_64_10 * h1_p1 - e_1 * fs_315_16_15 * h3_p1 - e_1 * fs_2835_32_10 * r_2 * h1_p1 - e_2 * fs_65_8_6 * h5_p1 + e_2 * fs_35_2_15 * r_2 * h3_p1 + e_2 * fs_405_8_10 * r_4 * h1_p1 + e_3 * fs_45_44_70 * h7_p1 + e_3 * fs_105_572_4290 * h7_p7 + e_3 * fs_15_4_6 * r_2 * h5_p1 - e_3 * fs_105_22_15 * r_4 * h3_p1 - e_3 * fs_45_4_10 * r_6 * h1_p1 + e_4 * fs_7707_9724_2 * h9_p1 - e_4 * fs_609_4862_143 * h9_p7 - e_4 * fs_45_187_70 * r_2 * h7_p1 - e_4 * fs_105_2431_4290 * r_2 * h7_p7 - e_4 * fs_1_2_6 * r_4 * h5_p1 + e_4 * fs_70_143_15 * r_6 * h3_p1 + e_4 * fs_45_44_10 * r_8 * h1_p1 + e_5 * fs_9_4199_165 * h11_p1 + e_5 * fs_9_4199_14586 * h11_p7 - e_5 * fs_367_4862_2 * r_2 * h9_p1 + e_5 * fs_29_2431_143 * r_2 * h9_p7 + e_5 * fs_45_3553_70 * r_4 * h7_p1 + e_5 * fs_105_46189_4290 * r_4 * h7_p7 + e_5 * fs_1_51_6 * r_6 * h5_p1 - e_5 * fs_7_429_15 * r_8 * h3_p1 - e_5 * fs_9_286_10 * r_10 * h1_p1;
    }

    // NOTE: the angular components are formed in 124 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph1_0, ph3_0, ph5_0, ph7_0, ph7_p6, ph9_0, ph9_p6, ph11_0, ph11_p6, ab_2, pc_29 : simd::cache_line_size())
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

        const auto h1_0 = ph1_0[k];
        const auto h3_0 = ph3_0[k];
        const auto h5_0 = ph5_0[k];
        const auto h7_0 = ph7_0[k];
        const auto h7_p6 = ph7_p6[k];
        const auto h9_0 = ph9_0[k];
        const auto h9_p6 = ph9_p6[k];
        const auto h11_0 = ph11_0[k];
        const auto h11_p6 = ph11_p6[k];

        pc_29[k] = e_0 * fs_2835_32_3 * h1_0 + e_1 * fs_315_16_3 * h3_0 - e_1 * fs_2835_16_3 * r_2 * h1_0 - e_2 * fs_125_4_3 * h5_0 - e_2 * fs_35_2_3 * r_2 * h3_0 + e_2 * fs_405_4_3 * r_4 * h1_0 + e_3 * fs_945_286_3 * h7_0 + e_3 * fs_45_286_2002 * h7_p6 + e_3 * fs_375_26_3 * r_2 * h5_0 + e_3 * fs_105_22_3 * r_4 * h3_0 - e_3 * fs_45_2_3 * r_6 * h1_0 + e_4 * fs_5565_4862_3 * h9_0 - e_4 * fs_84_2431_1430 * h9_p6 - e_4 * fs_1890_2431_3 * r_2 * h7_0 - e_4 * fs_90_2431_2002 * r_2 * h7_p6 - e_4 * fs_25_13_3 * r_4 * h5_0 - e_4 * fs_70_143_3 * r_6 * h3_0 + e_4 * fs_45_22_3 * r_8 * h1_0 + e_5 * fs_165_4199_3 * h11_0 + e_5 * fs_15_4199_4862 * h11_p6 - e_5 * fs_265_2431_3 * r_2 * h9_0 + e_5 * fs_8_2431_1430 * r_2 * h9_p6 + e_5 * fs_1890_46189_3 * r_4 * h7_0 + e_5 * fs_90_46189_2002 * r_4 * h7_p6 + e_5 * fs_50_663_3 * r_6 * h5_0 + e_5 * fs_7_429_3 * r_8 * h3_0 - e_5 * fs_9_143_3 * r_10 * h1_0;
    }

    // NOTE: the angular components are formed in 124 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph1_p1, ph3_p1, ph5_p1, ph5_p5, ph7_p1, ph7_p5, ph9_p1, ph9_p5, ph11_p1, ph11_p5, ab_2, pc_30 : simd::cache_line_size())
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

        const auto h1_p1 = ph1_p1[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h5_p5 = ph5_p5[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h7_p5 = ph7_p5[k];
        const auto h9_p1 = ph9_p1[k];
        const auto h9_p5 = ph9_p5[k];
        const auto h11_p1 = ph11_p1[k];
        const auto h11_p5 = ph11_p5[k];

        pc_30[k] = - e_0 * fs_945_32_3 * h1_p1 - e_1 * fs_2205_32_2 * h3_p1 + e_1 * fs_945_16_3 * r_2 * h1_p1 + e_2 * fs_55_4_5 * h5_p1 + e_2 * fs_25_4_42 * h5_p5 + e_2 * fs_245_4_2 * r_2 * h3_p1 - e_2 * fs_135_4_3 * r_4 * h1_p1 + e_3 * fs_105_286_21 * h7_p1 - e_3 * fs_45_572_77 * h7_p5 - e_3 * fs_165_26_5 * r_2 * h5_p1 - e_3 * fs_75_26_42 * r_2 * h5_p5 - e_3 * fs_735_44_2 * r_4 * h3_p1 + e_3 * fs_15_2_3 * r_6 * h1_p1 - e_4 * fs_1869_4862_15 * h9_p1 - e_4 * fs_105_4862_858 * h9_p5 - e_4 * fs_210_2431_21 * r_2 * h7_p1 + e_4 * fs_45_2431_77 * r_2 * h7_p5 + e_4 * fs_11_13_5 * r_4 * h5_p1 + e_4 * fs_5_13_42 * r_4 * h5_p5 + e_4 * fs_245_143_2 * r_6 * h3_p1 - e_4 * fs_15_22_3 * r_8 * h1_p1 - e_5 * fs_135_8398_22 * h11_p1 + e_5 * fs_45_4199_429 * h11_p5 + e_5 * fs_89_2431_15 * r_2 * h9_p1 + e_5 * fs_5_2431_858 * r_2 * h9_p5 + e_5 * fs_210_46189_21 * r_4 * h7_p1 - e_5 * fs_45_46189_77 * r_4 * h7_p5 - e_5 * fs_22_663_5 * r_6 * h5_p1 - e_5 * fs_10_663_42 * r_6 * h5_p5 - e_5 * fs_49_858_2 * r_8 * h3_p1 + e_5 * fs_3_143_3 * r_10 * h1_p1;
    }

    // NOTE: the angular components are formed in 124 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph3_p2, ph5_p2, ph5_p4, ph7_p2, ph7_p4, ph9_p2, ph9_p4, ph11_p2, ph11_p4, ab_2, pc_31 : simd::cache_line_size())
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

        const auto h3_p2 = ph3_p2[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h5_p4 = ph5_p4[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h7_p4 = ph7_p4[k];
        const auto h9_p2 = ph9_p2[k];
        const auto h9_p4 = ph9_p4[k];
        const auto h11_p2 = ph11_p2[k];
        const auto h11_p4 = ph11_p4[k];

        pc_31[k] = e_1 * fs_315_4_2 * h3_p2 - e_2 * fs_5_4_14 * h5_p2 + e_2 * fs_55_8_42 * h5_p4 - e_2 * fs_70_1_2 * r_2 * h3_p2 - e_3 * fs_315_286_35 * h7_p2 - e_3 * fs_135_572_770 * h7_p4 + e_3 * fs_15_26_14 * r_2 * h5_p2 - e_3 * fs_165_52_42 * r_2 * h5_p4 + e_3 * fs_210_11_2 * r_4 * h3_p2 + e_4 * fs_1113_4862_33 * h9_p2 + e_4 * fs_21_9724_6006 * h9_p4 + e_4 * fs_630_2431_35 * r_2 * h7_p2 + e_4 * fs_135_2431_770 * r_2 * h7_p4 - e_4 * fs_1_13_14 * r_4 * h5_p2 + e_4 * fs_11_26_42 * r_4 * h5_p4 - e_4 * fs_280_143_2 * r_6 * h3_p2 + e_5 * fs_27_4199_286 * h11_p2 + e_5 * fs_9_8398_30030 * h11_p4 - e_5 * fs_53_2431_33 * r_2 * h9_p2 - e_5 * fs_1_4862_6006 * r_2 * h9_p4 - e_5 * fs_630_46189_35 * r_4 * h7_p2 - e_5 * fs_135_46189_770 * r_4 * h7_p4 + e_5 * fs_2_663_14 * r_6 * h5_p2 - e_5 * fs_11_663_42 * r_6 * h5_p4 + e_5 * fs_28_429_2 * r_8 * h3_p2;
    }

    // NOTE: the angular components are formed in 124 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph3_m3, ph5_m3, ph7_m3, ph9_m3, ph11_m3, ab_2, pc_32 : simd::cache_line_size())
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

        const auto h3_m3 = ph3_m3[k];
        const auto h5_m3 = ph5_m3[k];
        const auto h7_m3 = ph7_m3[k];
        const auto h9_m3 = ph9_m3[k];
        const auto h11_m3 = ph11_m3[k];

        pc_32[k] = - e_1 * fs_315_8_7 * h3_m3 - e_2 * f_145_4 * h5_m3 + e_2 * fs_35_1_7 * r_2 * h3_m3 + e_3 * fs_315_143_30 * h7_m3 + e_3 * f_435_26 * r_2 * h5_m3 - e_3 * fs_105_11_7 * r_4 * h3_m3 - e_4 * fs_1029_4862_33 * h9_m3 - e_4 * fs_1260_2431_30 * r_2 * h7_m3 - e_4 * f_29_13 * r_4 * h5_m3 + e_4 * fs_140_143_7 * r_6 * h3_m3 - e_5 * fs_42_4199_429 * h11_m3 + e_5 * fs_49_2431_33 * r_2 * h9_m3 + e_5 * fs_1260_46189_30 * r_4 * h7_m3 + e_5 * f_58_663 * r_6 * h5_m3 - e_5 * fs_14_429_7 * r_8 * h3_m3;
    }

    // NOTE: the angular components are formed in 124 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph3_m2, ph5_m4, ph5_m2, ph7_m4, ph7_m2, ph9_m4, ph9_m2, ph11_m4, ph11_m2, ab_2, pc_33 : simd::cache_line_size())
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

        const auto h3_m2 = ph3_m2[k];
        const auto h5_m4 = ph5_m4[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h7_m4 = ph7_m4[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h9_m4 = ph9_m4[k];
        const auto h9_m2 = ph9_m2[k];
        const auto h11_m4 = ph11_m4[k];
        const auto h11_m2 = ph11_m2[k];

        pc_33[k] = e_1 * fs_315_4_2 * h3_m2 - e_2 * fs_55_8_42 * h5_m4 - e_2 * fs_5_4_14 * h5_m2 - e_2 * fs_70_1_2 * r_2 * h3_m2 + e_3 * fs_135_572_770 * h7_m4 - e_3 * fs_315_286_35 * h7_m2 + e_3 * fs_165_52_42 * r_2 * h5_m4 + e_3 * fs_15_26_14 * r_2 * h5_m2 + e_3 * fs_210_11_2 * r_4 * h3_m2 - e_4 * fs_21_9724_6006 * h9_m4 + e_4 * fs_1113_4862_33 * h9_m2 - e_4 * fs_135_2431_770 * r_2 * h7_m4 + e_4 * fs_630_2431_35 * r_2 * h7_m2 - e_4 * fs_11_26_42 * r_4 * h5_m4 - e_4 * fs_1_13_14 * r_4 * h5_m2 - e_4 * fs_280_143_2 * r_6 * h3_m2 - e_5 * fs_9_8398_30030 * h11_m4 + e_5 * fs_27_4199_286 * h11_m2 + e_5 * fs_1_4862_6006 * r_2 * h9_m4 - e_5 * fs_53_2431_33 * r_2 * h9_m2 + e_5 * fs_135_46189_770 * r_4 * h7_m4 - e_5 * fs_630_46189_35 * r_4 * h7_m2 + e_5 * fs_11_663_42 * r_6 * h5_m4 + e_5 * fs_2_663_14 * r_6 * h5_m2 + e_5 * fs_28_429_2 * r_8 * h3_m2;
    }

    // NOTE: the angular components are formed in 124 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph1_m1, ph3_m1, ph5_m5, ph5_m1, ph7_m5, ph7_m1, ph9_m5, ph9_m1, ph11_m5, ph11_m1, ab_2, pc_34 : simd::cache_line_size())
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

        const auto h1_m1 = ph1_m1[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h5_m5 = ph5_m5[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h7_m5 = ph7_m5[k];
        const auto h7_m1 = ph7_m1[k];
        const auto h9_m5 = ph9_m5[k];
        const auto h9_m1 = ph9_m1[k];
        const auto h11_m5 = ph11_m5[k];
        const auto h11_m1 = ph11_m1[k];

        pc_34[k] = - e_0 * fs_945_32_3 * h1_m1 - e_1 * fs_2205_32_2 * h3_m1 + e_1 * fs_945_16_3 * r_2 * h1_m1 - e_2 * fs_25_4_42 * h5_m5 + e_2 * fs_55_4_5 * h5_m1 + e_2 * fs_245_4_2 * r_2 * h3_m1 - e_2 * fs_135_4_3 * r_4 * h1_m1 + e_3 * fs_45_572_77 * h7_m5 + e_3 * fs_105_286_21 * h7_m1 + e_3 * fs_75_26_42 * r_2 * h5_m5 - e_3 * fs_165_26_5 * r_2 * h5_m1 - e_3 * fs_735_44_2 * r_4 * h3_m1 + e_3 * fs_15_2_3 * r_6 * h1_m1 + e_4 * fs_105_4862_858 * h9_m5 - e_4 * fs_1869_4862_15 * h9_m1 - e_4 * fs_45_2431_77 * r_2 * h7_m5 - e_4 * fs_210_2431_21 * r_2 * h7_m1 - e_4 * fs_5_13_42 * r_4 * h5_m5 + e_4 * fs_11_13_5 * r_4 * h5_m1 + e_4 * fs_245_143_2 * r_6 * h3_m1 - e_4 * fs_15_22_3 * r_8 * h1_m1 - e_5 * fs_45_4199_429 * h11_m5 - e_5 * fs_135_8398_22 * h11_m1 - e_5 * fs_5_2431_858 * r_2 * h9_m5 + e_5 * fs_89_2431_15 * r_2 * h9_m1 + e_5 * fs_45_46189_77 * r_4 * h7_m5 + e_5 * fs_210_46189_21 * r_4 * h7_m1 + e_5 * fs_10_663_42 * r_6 * h5_m5 - e_5 * fs_22_663_5 * r_6 * h5_m1 - e_5 * fs_49_858_2 * r_8 * h3_m1 + e_5 * fs_3_143_3 * r_10 * h1_m1;
    }

    // NOTE: the angular components are formed in 124 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_3, pe_4, pe_5, ph7_m6, ph9_m6, ph11_m6, ab_2, pc_35 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

        const auto h7_m6 = ph7_m6[k];
        const auto h9_m6 = ph9_m6[k];
        const auto h11_m6 = ph11_m6[k];

        pc_35[k] = - e_3 * fs_45_286_2002 * h7_m6 + e_4 * fs_84_2431_1430 * h9_m6 + e_4 * fs_90_2431_2002 * r_2 * h7_m6 - e_5 * fs_15_4199_4862 * h11_m6 - e_5 * fs_8_2431_1430 * r_2 * h9_m6 - e_5 * fs_90_46189_2002 * r_4 * h7_m6;
    }

    // NOTE: the angular components are formed in 124 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph1_m1, ph3_m1, ph5_m1, ph7_m7, ph7_m1, ph9_m7, ph9_m1, ph11_m7, ph11_m1, ab_2, pc_36 : simd::cache_line_size())
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

        const auto h1_m1 = ph1_m1[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h7_m7 = ph7_m7[k];
        const auto h7_m1 = ph7_m1[k];
        const auto h9_m7 = ph9_m7[k];
        const auto h9_m1 = ph9_m1[k];
        const auto h11_m7 = ph11_m7[k];
        const auto h11_m1 = ph11_m1[k];

        pc_36[k] = - e_0 * fs_2835_64_10 * h1_m1 + e_1 * fs_315_16_15 * h3_m1 + e_1 * fs_2835_32_10 * r_2 * h1_m1 + e_2 * fs_65_8_6 * h5_m1 - e_2 * fs_35_2_15 * r_2 * h3_m1 - e_2 * fs_405_8_10 * r_4 * h1_m1 - e_3 * fs_105_572_4290 * h7_m7 - e_3 * fs_45_44_70 * h7_m1 - e_3 * fs_15_4_6 * r_2 * h5_m1 + e_3 * fs_105_22_15 * r_4 * h3_m1 + e_3 * fs_45_4_10 * r_6 * h1_m1 + e_4 * fs_609_4862_143 * h9_m7 - e_4 * fs_7707_9724_2 * h9_m1 + e_4 * fs_105_2431_4290 * r_2 * h7_m7 + e_4 * fs_45_187_70 * r_2 * h7_m1 + e_4 * fs_1_2_6 * r_4 * h5_m1 - e_4 * fs_70_143_15 * r_6 * h3_m1 - e_4 * fs_45_44_10 * r_8 * h1_m1 - e_5 * fs_9_4199_14586 * h11_m7 - e_5 * fs_9_4199_165 * h11_m1 - e_5 * fs_29_2431_143 * r_2 * h9_m7 + e_5 * fs_367_4862_2 * r_2 * h9_m1 - e_5 * fs_105_46189_4290 * r_4 * h7_m7 - e_5 * fs_45_3553_70 * r_4 * h7_m1 - e_5 * fs_1_51_6 * r_6 * h5_m1 + e_5 * fs_7_429_15 * r_8 * h3_m1 + e_5 * fs_9_286_10 * r_10 * h1_m1;
    }

    // NOTE: the angular components are formed in 124 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph3_m2, ph5_m2, ph7_m2, ph9_m8, ph9_m2, ph11_m8, ph11_m2, ab_2, pc_37 : simd::cache_line_size())
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

        const auto h3_m2 = ph3_m2[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h9_m8 = ph9_m8[k];
        const auto h9_m2 = ph9_m2[k];
        const auto h11_m8 = ph11_m8[k];
        const auto h11_m2 = ph11_m2[k];

        pc_37[k] = e_1 * fs_315_16_33 * h3_m2 - e_2 * fs_5_4_231 * h5_m2 - e_2 * fs_35_2_33 * r_2 * h3_m2 - e_3 * fs_30_143_2310 * h7_m2 + e_3 * fs_15_26_231 * r_2 * h5_m2 + e_3 * fs_105_22_33 * r_4 * h3_m2 + e_4 * fs_21_442_221 * h9_m8 - e_4 * fs_231_442_2 * h9_m2 + e_4 * fs_120_2431_2310 * r_2 * h7_m2 - e_4 * fs_1_13_231 * r_4 * h5_m2 - e_4 * fs_70_143_33 * r_6 * h3_m2 - e_5 * fs_9_4199_12597 * h11_m8 - e_5 * fs_9_4199_39 * h11_m2 - e_5 * fs_1_221_221 * r_2 * h9_m8 + e_5 * fs_11_221_2 * r_2 * h9_m2 - e_5 * fs_120_46189_2310 * r_4 * h7_m2 + e_5 * fs_2_663_231 * r_6 * h5_m2 + e_5 * fs_7_429_33 * r_8 * h3_m2;
    }

    // NOTE: the angular components are formed in 124 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph3_m3, ph5_m3, ph7_m3, ph9_m9, ph9_m3, ph11_m9, ph11_m3, ab_2, pc_38 : simd::cache_line_size())
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

        const auto h3_m3 = ph3_m3[k];
        const auto h5_m3 = ph5_m3[k];
        const auto h7_m3 = ph7_m3[k];
        const auto h9_m9 = ph9_m9[k];
        const auto h9_m3 = ph9_m3[k];
        const auto h11_m9 = ph11_m9[k];
        const auto h11_m3 = ph11_m3[k];

        pc_38[k] = - e_1 * fs_315_32_66 * h3_m3 - e_2 * fs_5_2_462 * h5_m3 + e_2 * fs_35_4_66 * r_2 * h3_m3 - e_3 * fs_225_572_385 * h7_m3 + e_3 * fs_15_13_462 * r_2 * h5_m3 - e_3 * fs_105_44_66 * r_4 * h3_m3 - e_4 * fs_21_442_1326 * h9_m9 - e_4 * fs_21_221_14 * h9_m3 + e_4 * fs_225_2431_385 * r_2 * h7_m3 - e_4 * fs_2_13_462 * r_4 * h5_m3 + e_4 * fs_35_143_66 * r_6 * h3_m3 - e_5 * fs_3_4199_62985 * h11_m9 - e_5 * fs_3_8398_182 * h11_m3 + e_5 * fs_1_221_1326 * r_2 * h9_m9 + e_5 * fs_2_221_14 * r_2 * h9_m3 - e_5 * fs_225_46189_385 * r_4 * h7_m3 + e_5 * fs_4_663_462 * r_6 * h5_m3 - e_5 * fs_7_858_66 * r_8 * h3_m3;
    }

    // NOTE: the angular components are formed in 124 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_2, pe_3, pe_4, pe_5, ph5_p4, ph7_p4, ph9_p4, ph9_p8, ph11_p4, ph11_p8, ab_2, pc_39 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h5_p4 = ph5_p4[k];
        const auto h7_p4 = ph7_p4[k];
        const auto h9_p4 = ph9_p4[k];
        const auto h9_p8 = ph9_p8[k];
        const auto h11_p4 = ph11_p4[k];
        const auto h11_p8 = ph11_p8[k];

        pc_39[k] = - e_2 * fs_15_4_154 * h5_p4 - e_3 * fs_75_104_210 * h7_p4 + e_3 * fs_45_26_154 * r_2 * h5_p4 - e_4 * fs_21_442_182 * h9_p4 + e_4 * fs_21_221_442 * h9_p8 + e_4 * fs_75_442_210 * r_2 * h7_p4 - e_4 * fs_3_13_154 * r_4 * h5_p4 - e_5 * fs_3_8398_910 * h11_p4 + e_5 * fs_3_4199_25194 * h11_p8 + e_5 * fs_1_221_182 * r_2 * h9_p4 - e_5 * fs_2_221_442 * r_2 * h9_p8 - e_5 * fs_75_8398_210 * r_4 * h7_p4 + e_5 * fs_2_221_154 * r_6 * h5_p4;
    }

    // NOTE: the angular components are formed in 124 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph3_p3, ph5_p3, ph7_p3, ph7_p7, ph9_p3, ph9_p7, ph11_p3, ph11_p7, ab_2, pc_40 : simd::cache_line_size())
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

        const auto h3_p3 = ph3_p3[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h7_p3 = ph7_p3[k];
        const auto h7_p7 = ph7_p7[k];
        const auto h9_p3 = ph9_p3[k];
        const auto h9_p7 = ph9_p7[k];
        const auto h11_p3 = ph11_p3[k];
        const auto h11_p7 = ph11_p7[k];

        pc_40[k] = e_1 * fs_315_16_33 * h3_p3 + e_2 * fs_5_4_231 * h5_p3 - e_2 * fs_35_2_33 * r_2 * h3_p3 - e_3 * fs_375_1144_770 * h7_p3 - e_3 * fs_105_104_130 * h7_p7 - e_3 * fs_15_26_231 * r_2 * h5_p3 + e_3 * fs_105_22_33 * r_4 * h3_p3 - e_4 * fs_189_442_7 * h9_p3 + e_4 * fs_21_221_39 * h9_p7 + e_4 * fs_375_4862_770 * r_2 * h7_p3 + e_4 * fs_105_442_130 * r_2 * h7_p7 + e_4 * fs_1_13_231 * r_4 * h5_p3 - e_4 * fs_70_143_33 * r_6 * h3_p3 - e_5 * fs_12_4199_91 * h11_p3 + e_5 * fs_36_4199_442 * h11_p7 + e_5 * fs_9_221_7 * r_2 * h9_p3 - e_5 * fs_2_221_39 * r_2 * h9_p7 - e_5 * fs_375_92378_770 * r_4 * h7_p3 - e_5 * fs_105_8398_130 * r_4 * h7_p7 - e_5 * fs_2_663_231 * r_6 * h5_p3 + e_5 * fs_7_429_33 * r_8 * h3_p3;
    }

    // NOTE: the angular components are formed in 124 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph3_p2, ph5_p2, ph7_p2, ph7_p6, ph9_p2, ph9_p6, ph11_p2, ph11_p6, ab_2, pc_41 : simd::cache_line_size())
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

        const auto h3_p2 = ph3_p2[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h7_p6 = ph7_p6[k];
        const auto h9_p2 = ph9_p2[k];
        const auto h9_p6 = ph9_p6[k];
        const auto h11_p2 = ph11_p2[k];
        const auto h11_p6 = ph11_p6[k];

        pc_41[k] = - e_1 * f_945_16 * h3_p2 + e_2 * fs_15_1_7 * h5_p2 + e_2 * f_105_2 * r_2 * h3_p2 - e_3 * fs_435_1144_70 * h7_p2 + e_3 * fs_15_1144_10010 * h7_p6 - e_3 * fs_90_13_7 * r_2 * h5_p2 - e_3 * f_315_22 * r_4 * h3_p2 - e_4 * fs_861_4862_66 * h9_p2 - e_4 * fs_189_4862_286 * h9_p6 + e_4 * fs_435_4862_70 * r_2 * h7_p2 - e_4 * fs_15_4862_10010 * r_2 * h7_p6 + e_4 * fs_12_13_7 * r_4 * h5_p2 + e_4 * f_210_143 * r_6 * h3_p2 - e_5 * fs_18_4199_143 * h11_p2 + e_5 * fs_6_4199_24310 * h11_p6 + e_5 * fs_41_2431_66 * r_2 * h9_p2 + e_5 * fs_9_2431_286 * r_2 * h9_p6 - e_5 * fs_435_92378_70 * r_4 * h7_p2 + e_5 * fs_15_92378_10010 * r_4 * h7_p6 - e_5 * fs_8_221_7 * r_6 * h5_p2 - e_5 * f_7_143 * r_8 * h3_p2;
    }

    // NOTE: the angular components are formed in 124 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph1_p1, ph3_p1, ph5_p1, ph5_p5, ph7_p1, ph7_p5, ph9_p1, ph9_p5, ph11_p1, ph11_p5, ab_2, pc_42 : simd::cache_line_size())
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

        const auto h1_p1 = ph1_p1[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h5_p5 = ph5_p5[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h7_p5 = ph7_p5[k];
        const auto h9_p1 = ph9_p1[k];
        const auto h9_p5 = ph9_p5[k];
        const auto h11_p1 = ph11_p1[k];
        const auto h11_p5 = ph11_p5[k];

        pc_42[k] = e_0 * fs_2835_32_2 * h1_p1 - e_1 * fs_945_16_3 * h3_p1 - e_1 * fs_2835_16_2 * r_2 * h1_p1 + e_2 * fs_15_4_30 * h5_p1 - e_2 * fs_75_4_7 * h5_p5 + e_2 * fs_105_2_3 * r_2 * h3_p1 + e_2 * fs_405_4_2 * r_4 * h1_p1 + e_3 * fs_1035_1144_14 * h7_p1 + e_3 * fs_435_1144_462 * h7_p5 - e_3 * fs_45_26_30 * r_2 * h5_p1 + e_3 * fs_225_26_7 * r_2 * h5_p5 - e_3 * fs_315_22_3 * r_4 * h3_p1 - e_3 * fs_45_2_2 * r_6 * h1_p1 - e_4 * fs_2247_4862_10 * h9_p1 - e_4 * fs_525_4862_143 * h9_p5 - e_4 * fs_1035_4862_14 * r_2 * h7_p1 - e_4 * fs_435_4862_462 * r_2 * h7_p5 + e_4 * fs_3_13_30 * r_4 * h5_p1 - e_4 * fs_15_13_7 * r_4 * h5_p5 + e_4 * fs_210_143_3 * r_6 * h3_p1 + e_4 * fs_45_22_2 * r_8 * h1_p1 - e_5 * fs_60_4199_33 * h11_p1 + e_5 * fs_60_4199_286 * h11_p5 + e_5 * fs_107_2431_10 * r_2 * h9_p1 + e_5 * fs_25_2431_143 * r_2 * h9_p5 + e_5 * fs_1035_92378_14 * r_4 * h7_p1 + e_5 * fs_435_92378_462 * r_4 * h7_p5 - e_5 * fs_2_221_30 * r_6 * h5_p1 + e_5 * fs_10_221_7 * r_6 * h5_p5 - e_5 * fs_7_143_3 * r_8 * h3_p1 - e_5 * fs_9_143_2 * r_10 * h1_p1;
    }

    // NOTE: the angular components are formed in 124 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph1_0, ph3_0, ph5_0, ph5_p4, ph7_0, ph7_p4, ph9_0, ph9_p4, ph11_0, ph11_p4, ab_2, pc_43 : simd::cache_line_size())
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

        const auto h1_0 = ph1_0[k];
        const auto h3_0 = ph3_0[k];
        const auto h5_0 = ph5_0[k];
        const auto h5_p4 = ph5_p4[k];
        const auto h7_0 = ph7_0[k];
        const auto h7_p4 = ph7_p4[k];
        const auto h9_0 = ph9_0[k];
        const auto h9_p4 = ph9_p4[k];
        const auto h11_0 = ph11_0[k];
        const auto h11_p4 = ph11_p4[k];

        pc_43[k] = e_0 * fs_945_8_2 * h1_0 - e_1 * fs_315_8_2 * h3_0 - e_1 * fs_945_4_2 * r_2 * h1_0 - e_2 * fs_25_2_2 * h5_0 - e_2 * fs_15_4_70 * h5_p4 + e_2 * fs_35_1_2 * r_2 * h3_0 + e_2 * fs_135_1_2 * r_4 * h1_0 + e_3 * fs_4305_572_2 * h7_0 + e_3 * fs_405_1144_462 * h7_p4 + e_3 * fs_75_13_2 * r_2 * h5_0 + e_3 * fs_45_26_70 * r_2 * h5_p4 - e_3 * fs_105_11_2 * r_4 * h3_0 - e_3 * fs_30_1_2 * r_6 * h1_0 - e_4 * fs_2835_2431_2 * h9_0 - e_4 * fs_63_4862_10010 * h9_p4 - e_4 * fs_4305_2431_2 * r_2 * h7_0 - e_4 * fs_405_4862_462 * r_2 * h7_p4 - e_4 * fs_10_13_2 * r_4 * h5_0 - e_4 * fs_3_13_70 * r_4 * h5_p4 + e_4 * fs_140_143_2 * r_6 * h3_0 + e_4 * fs_30_11_2 * r_8 * h1_0 - e_5 * fs_495_4199_2 * h11_0 + e_5 * fs_45_8398_2002 * h11_p4 + e_5 * fs_270_2431_2 * r_2 * h9_0 + e_5 * fs_3_2431_10010 * r_2 * h9_p4 + e_5 * fs_4305_46189_2 * r_4 * h7_0 + e_5 * fs_405_92378_462 * r_4 * h7_p4 + e_5 * fs_20_663_2 * r_6 * h5_0 + e_5 * fs_2_221_70 * r_6 * h5_p4 - e_5 * fs_14_429_2 * r_8 * h3_0 - e_5 * fs_12_143_2 * r_10 * h1_0;
    }

    // NOTE: the angular components are formed in 124 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph1_p1, ph3_p1, ph3_p3, ph5_p1, ph5_p3, ph7_p1, ph7_p3, ph9_p1, ph9_p3, ph11_p1, ph11_p3, ab_2, pc_44 : simd::cache_line_size())
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

        const auto h1_p1 = ph1_p1[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h3_p3 = ph3_p3[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h7_p3 = ph7_p3[k];
        const auto h9_p1 = ph9_p1[k];
        const auto h9_p3 = ph9_p3[k];
        const auto h11_p1 = ph11_p1[k];
        const auto h11_p3 = ph11_p3[k];

        pc_44[k] = - e_0 * fs_945_32_5 * h1_p1 - e_1 * fs_315_32_30 * h3_p1 - e_1 * fs_2205_32_2 * h3_p3 + e_1 * fs_945_16_5 * r_2 * h1_p1 + e_2 * fs_20_1_3 * h5_p1 + e_2 * fs_5_4_14 * h5_p3 + e_2 * fs_35_4_30 * r_2 * h3_p1 + e_2 * fs_245_4_2 * r_2 * h3_p3 - e_2 * fs_135_4_5 * r_4 * h1_p1 - e_3 * fs_735_572_35 * h7_p1 + e_3 * fs_135_572_105 * h7_p3 - e_3 * fs_120_13_3 * r_2 * h5_p1 - e_3 * fs_15_26_14 * r_2 * h5_p3 - e_3 * fs_105_44_30 * r_4 * h3_p1 - e_3 * fs_735_44_2 * r_4 * h3_p3 + e_3 * fs_15_2_5 * r_6 * h1_p1 + e_4 * f_1449_2431 * h9_p1 - e_4 * fs_189_4862_462 * h9_p3 + e_4 * fs_735_2431_35 * r_2 * h7_p1 - e_4 * fs_135_2431_105 * r_2 * h7_p3 + e_4 * fs_16_13_3 * r_4 * h5_p1 + e_4 * fs_1_13_14 * r_4 * h5_p3 + e_4 * fs_35_143_30 * r_6 * h3_p1 + e_4 * fs_245_143_2 * r_6 * h3_p3 - e_4 * fs_15_22_5 * r_8 * h1_p1 + e_5 * fs_36_4199_330 * h11_p1 + e_5 * fs_12_4199_6006 * h11_p3 - e_5 * f_138_2431 * r_2 * h9_p1 + e_5 * fs_9_2431_462 * r_2 * h9_p3 - e_5 * fs_735_46189_35 * r_4 * h7_p1 + e_5 * fs_135_46189_105 * r_4 * h7_p3 - e_5 * fs_32_663_3 * r_6 * h5_p1 - e_5 * fs_2_663_14 * r_6 * h5_p3 - e_5 * fs_7_858_30 * r_8 * h3_p1 - e_5 * fs_49_858_2 * r_8 * h3_p3 + e_5 * fs_3_143_5 * r_10 * h1_p1;
    }

    // NOTE: the angular components are formed in 124 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph3_m2, ph5_m2, ph7_m2, ph9_m2, ph11_m2, ab_2, pc_45 : simd::cache_line_size())
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

        const auto h3_m2 = ph3_m2[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h9_m2 = ph9_m2[k];
        const auto h11_m2 = ph11_m2[k];

        pc_45[k] = e_1 * fs_945_16_7 * h3_m2 - e_2 * f_45_1 * h5_m2 - e_2 * fs_105_2_7 * r_2 * h3_m2 + e_3 * fs_945_572_10 * h7_m2 + e_3 * f_270_13 * r_2 * h5_m2 + e_3 * fs_315_22_7 * r_4 * h3_m2 + e_4 * fs_21_2431_462 * h9_m2 - e_4 * fs_945_2431_10 * r_2 * h7_m2 - e_4 * f_36_13 * r_4 * h5_m2 - e_4 * fs_210_143_7 * r_6 * h3_m2 - e_5 * fs_36_4199_1001 * h11_m2 - e_5 * fs_2_2431_462 * r_2 * h9_m2 + e_5 * fs_945_46189_10 * r_4 * h7_m2 + e_5 * f_24_221 * r_6 * h5_m2 + e_5 * fs_7_143_7 * r_8 * h3_m2;
    }

    // NOTE: the angular components are formed in 124 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph1_m1, ph3_m3, ph3_m1, ph5_m3, ph5_m1, ph7_m3, ph7_m1, ph9_m3, ph9_m1, ph11_m3, ph11_m1, ab_2, pc_46 : simd::cache_line_size())
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

        const auto h1_m1 = ph1_m1[k];
        const auto h3_m3 = ph3_m3[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h5_m3 = ph5_m3[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h7_m3 = ph7_m3[k];
        const auto h7_m1 = ph7_m1[k];
        const auto h9_m3 = ph9_m3[k];
        const auto h9_m1 = ph9_m1[k];
        const auto h11_m3 = ph11_m3[k];
        const auto h11_m1 = ph11_m1[k];

        pc_46[k] = - e_0 * fs_945_32_5 * h1_m1 + e_1 * fs_2205_32_2 * h3_m3 - e_1 * fs_315_32_30 * h3_m1 + e_1 * fs_945_16_5 * r_2 * h1_m1 - e_2 * fs_5_4_14 * h5_m3 + e_2 * fs_20_1_3 * h5_m1 - e_2 * fs_245_4_2 * r_2 * h3_m3 + e_2 * fs_35_4_30 * r_2 * h3_m1 - e_2 * fs_135_4_5 * r_4 * h1_m1 - e_3 * fs_135_572_105 * h7_m3 - e_3 * fs_735_572_35 * h7_m1 + e_3 * fs_15_26_14 * r_2 * h5_m3 - e_3 * fs_120_13_3 * r_2 * h5_m1 + e_3 * fs_735_44_2 * r_4 * h3_m3 - e_3 * fs_105_44_30 * r_4 * h3_m1 + e_3 * fs_15_2_5 * r_6 * h1_m1 + e_4 * fs_189_4862_462 * h9_m3 + e_4 * f_1449_2431 * h9_m1 + e_4 * fs_135_2431_105 * r_2 * h7_m3 + e_4 * fs_735_2431_35 * r_2 * h7_m1 - e_4 * fs_1_13_14 * r_4 * h5_m3 + e_4 * fs_16_13_3 * r_4 * h5_m1 - e_4 * fs_245_143_2 * r_6 * h3_m3 + e_4 * fs_35_143_30 * r_6 * h3_m1 - e_4 * fs_15_22_5 * r_8 * h1_m1 - e_5 * fs_12_4199_6006 * h11_m3 + e_5 * fs_36_4199_330 * h11_m1 - e_5 * fs_9_2431_462 * r_2 * h9_m3 - e_5 * f_138_2431 * r_2 * h9_m1 - e_5 * fs_135_46189_105 * r_4 * h7_m3 - e_5 * fs_735_46189_35 * r_4 * h7_m1 + e_5 * fs_2_663_14 * r_6 * h5_m3 - e_5 * fs_32_663_3 * r_6 * h5_m1 + e_5 * fs_49_858_2 * r_8 * h3_m3 - e_5 * fs_7_858_30 * r_8 * h3_m1 + e_5 * fs_3_143_5 * r_10 * h1_m1;
    }

    // NOTE: the angular components are formed in 124 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_2, pe_3, pe_4, pe_5, ph5_m4, ph7_m4, ph9_m4, ph11_m4, ab_2, pc_47 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h5_m4 = ph5_m4[k];
        const auto h7_m4 = ph7_m4[k];
        const auto h9_m4 = ph9_m4[k];
        const auto h11_m4 = ph11_m4[k];

        pc_47[k] = e_2 * fs_15_4_70 * h5_m4 - e_3 * fs_405_1144_462 * h7_m4 - e_3 * fs_45_26_70 * r_2 * h5_m4 + e_4 * fs_63_4862_10010 * h9_m4 + e_4 * fs_405_4862_462 * r_2 * h7_m4 + e_4 * fs_3_13_70 * r_4 * h5_m4 - e_5 * fs_45_8398_2002 * h11_m4 - e_5 * fs_3_2431_10010 * r_2 * h9_m4 - e_5 * fs_405_92378_462 * r_4 * h7_m4 - e_5 * fs_2_221_70 * r_6 * h5_m4;
    }

    // NOTE: the angular components are formed in 124 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph1_m1, ph3_m1, ph5_m5, ph5_m1, ph7_m5, ph7_m1, ph9_m5, ph9_m1, ph11_m5, ph11_m1, ab_2, pc_48 : simd::cache_line_size())
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

        const auto h1_m1 = ph1_m1[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h5_m5 = ph5_m5[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h7_m5 = ph7_m5[k];
        const auto h7_m1 = ph7_m1[k];
        const auto h9_m5 = ph9_m5[k];
        const auto h9_m1 = ph9_m1[k];
        const auto h11_m5 = ph11_m5[k];
        const auto h11_m1 = ph11_m1[k];

        pc_48[k] = - e_0 * fs_2835_32_2 * h1_m1 + e_1 * fs_945_16_3 * h3_m1 + e_1 * fs_2835_16_2 * r_2 * h1_m1 + e_2 * fs_75_4_7 * h5_m5 - e_2 * fs_15_4_30 * h5_m1 - e_2 * fs_105_2_3 * r_2 * h3_m1 - e_2 * fs_405_4_2 * r_4 * h1_m1 - e_3 * fs_435_1144_462 * h7_m5 - e_3 * fs_1035_1144_14 * h7_m1 - e_3 * fs_225_26_7 * r_2 * h5_m5 + e_3 * fs_45_26_30 * r_2 * h5_m1 + e_3 * fs_315_22_3 * r_4 * h3_m1 + e_3 * fs_45_2_2 * r_6 * h1_m1 + e_4 * fs_525_4862_143 * h9_m5 + e_4 * fs_2247_4862_10 * h9_m1 + e_4 * fs_435_4862_462 * r_2 * h7_m5 + e_4 * fs_1035_4862_14 * r_2 * h7_m1 + e_4 * fs_15_13_7 * r_4 * h5_m5 - e_4 * fs_3_13_30 * r_4 * h5_m1 - e_4 * fs_210_143_3 * r_6 * h3_m1 - e_4 * fs_45_22_2 * r_8 * h1_m1 - e_5 * fs_60_4199_286 * h11_m5 + e_5 * fs_60_4199_33 * h11_m1 - e_5 * fs_25_2431_143 * r_2 * h9_m5 - e_5 * fs_107_2431_10 * r_2 * h9_m1 - e_5 * fs_435_92378_462 * r_4 * h7_m5 - e_5 * fs_1035_92378_14 * r_4 * h7_m1 - e_5 * fs_10_221_7 * r_6 * h5_m5 + e_5 * fs_2_221_30 * r_6 * h5_m1 + e_5 * fs_7_143_3 * r_8 * h3_m1 + e_5 * fs_9_143_2 * r_10 * h1_m1;
    }

    // NOTE: the angular components are formed in 124 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph3_m2, ph5_m2, ph7_m6, ph7_m2, ph9_m6, ph9_m2, ph11_m6, ph11_m2, ab_2, pc_49 : simd::cache_line_size())
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

        const auto h3_m2 = ph3_m2[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h7_m6 = ph7_m6[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h9_m6 = ph9_m6[k];
        const auto h9_m2 = ph9_m2[k];
        const auto h11_m6 = ph11_m6[k];
        const auto h11_m2 = ph11_m2[k];

        pc_49[k] = e_1 * f_945_16 * h3_m2 - e_2 * fs_15_1_7 * h5_m2 - e_2 * f_105_2 * r_2 * h3_m2 - e_3 * fs_15_1144_10010 * h7_m6 + e_3 * fs_435_1144_70 * h7_m2 + e_3 * fs_90_13_7 * r_2 * h5_m2 + e_3 * f_315_22 * r_4 * h3_m2 + e_4 * fs_189_4862_286 * h9_m6 + e_4 * fs_861_4862_66 * h9_m2 + e_4 * fs_15_4862_10010 * r_2 * h7_m6 - e_4 * fs_435_4862_70 * r_2 * h7_m2 - e_4 * fs_12_13_7 * r_4 * h5_m2 - e_4 * f_210_143 * r_6 * h3_m2 - e_5 * fs_6_4199_24310 * h11_m6 + e_5 * fs_18_4199_143 * h11_m2 - e_5 * fs_9_2431_286 * r_2 * h9_m6 - e_5 * fs_41_2431_66 * r_2 * h9_m2 - e_5 * fs_15_92378_10010 * r_4 * h7_m6 + e_5 * fs_435_92378_70 * r_4 * h7_m2 + e_5 * fs_8_221_7 * r_6 * h5_m2 + e_5 * f_7_143 * r_8 * h3_m2;
    }

    // NOTE: the angular components are formed in 124 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph3_m3, ph5_m3, ph7_m7, ph7_m3, ph9_m7, ph9_m3, ph11_m7, ph11_m3, ab_2, pc_50 : simd::cache_line_size())
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

        const auto h3_m3 = ph3_m3[k];
        const auto h5_m3 = ph5_m3[k];
        const auto h7_m7 = ph7_m7[k];
        const auto h7_m3 = ph7_m3[k];
        const auto h9_m7 = ph9_m7[k];
        const auto h9_m3 = ph9_m3[k];
        const auto h11_m7 = ph11_m7[k];
        const auto h11_m3 = ph11_m3[k];

        pc_50[k] = - e_1 * fs_315_16_33 * h3_m3 - e_2 * fs_5_4_231 * h5_m3 + e_2 * fs_35_2_33 * r_2 * h3_m3 + e_3 * fs_105_104_130 * h7_m7 + e_3 * fs_375_1144_770 * h7_m3 + e_3 * fs_15_26_231 * r_2 * h5_m3 - e_3 * fs_105_22_33 * r_4 * h3_m3 - e_4 * fs_21_221_39 * h9_m7 + e_4 * fs_189_442_7 * h9_m3 - e_4 * fs_105_442_130 * r_2 * h7_m7 - e_4 * fs_375_4862_770 * r_2 * h7_m3 - e_4 * fs_1_13_231 * r_4 * h5_m3 + e_4 * fs_70_143_33 * r_6 * h3_m3 - e_5 * fs_36_4199_442 * h11_m7 + e_5 * fs_12_4199_91 * h11_m3 + e_5 * fs_2_221_39 * r_2 * h9_m7 - e_5 * fs_9_221_7 * r_2 * h9_m3 + e_5 * fs_105_8398_130 * r_4 * h7_m7 + e_5 * fs_375_92378_770 * r_4 * h7_m3 + e_5 * fs_2_663_231 * r_6 * h5_m3 - e_5 * fs_7_429_33 * r_8 * h3_m3;
    }

    // NOTE: the angular components are formed in 124 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_2, pe_3, pe_4, pe_5, ph5_m4, ph5_p5, ph7_m4, ph7_p5, ph7_p7, ph9_m8, ph9_m4, ph9_p5, ph9_p7, ph11_m8, ph11_m4, ph11_p5, ph11_p7, ab_2, pc_51, pc_52 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h5_m4 = ph5_m4[k];
        const auto h5_p5 = ph5_p5[k];
        const auto h7_m4 = ph7_m4[k];
        const auto h7_p5 = ph7_p5[k];
        const auto h7_p7 = ph7_p7[k];
        const auto h9_m8 = ph9_m8[k];
        const auto h9_m4 = ph9_m4[k];
        const auto h9_p5 = ph9_p5[k];
        const auto h9_p7 = ph9_p7[k];
        const auto h11_m8 = ph11_m8[k];
        const auto h11_m4 = ph11_m4[k];
        const auto h11_p5 = ph11_p5[k];
        const auto h11_p7 = ph11_p7[k];

        pc_51[k] = e_2 * fs_15_4_154 * h5_m4 + e_3 * fs_75_104_210 * h7_m4 - e_3 * fs_45_26_154 * r_2 * h5_m4 - e_4 * fs_21_221_442 * h9_m8 + e_4 * fs_21_442_182 * h9_m4 - e_4 * fs_75_442_210 * r_2 * h7_m4 + e_4 * fs_3_13_154 * r_4 * h5_m4 - e_5 * fs_3_4199_25194 * h11_m8 + e_5 * fs_3_8398_910 * h11_m4 + e_5 * fs_2_221_442 * r_2 * h9_m8 - e_5 * fs_1_221_182 * r_2 * h9_m4 + e_5 * fs_75_8398_210 * r_4 * h7_m4 - e_5 * fs_2_221_154 * r_6 * h5_m4;

        pc_52[k] = e_2 * fs_15_4_55 * h5_p5 + e_3 * fs_225_104_30 * h7_p5 + e_3 * fs_15_104_2730 * h7_p7 - e_3 * fs_45_26_55 * r_2 * h5_p5 + e_4 * fs_21_442_455 * h9_p5 + e_4 * fs_42_221_91 * h9_p7 - e_4 * fs_225_442_30 * r_2 * h7_p5 - e_4 * fs_15_442_2730 * r_2 * h7_p7 + e_4 * fs_3_13_55 * r_4 * h5_p5 + e_5 * fs_3_4199_910 * h11_p5 + e_5 * fs_3_4199_9282 * h11_p7 - e_5 * fs_1_221_455 * r_2 * h9_p5 - e_5 * fs_4_221_91 * r_2 * h9_p7 + e_5 * fs_225_8398_30 * r_4 * h7_p5 + e_5 * fs_15_8398_2730 * r_4 * h7_p7 - e_5 * fs_2_221_55 * r_6 * h5_p5;
    }

    // NOTE: the angular components are formed in 124 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_2, pe_3, pe_4, pe_5, ph5_p4, ph7_p4, ph7_p6, ph9_p4, ph9_p6, ph11_p4, ph11_p6, ab_2, pc_53 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h5_p4 = ph5_p4[k];
        const auto h7_p4 = ph7_p4[k];
        const auto h7_p6 = ph7_p6[k];
        const auto h9_p4 = ph9_p4[k];
        const auto h9_p6 = ph9_p6[k];
        const auto h11_p4 = ph11_p4[k];
        const auto h11_p6 = ph11_p6[k];

        pc_53[k] = - e_2 * fs_45_8_66 * h5_p4 + e_3 * fs_75_52_10 * h7_p4 - e_3 * fs_15_13_65 * h7_p6 + e_3 * fs_135_52_66 * r_2 * h5_p4 + e_4 * fs_147_884_78 * h9_p4 + e_4 * fs_63_442_91 * h9_p6 - e_4 * fs_75_221_10 * r_2 * h7_p4 + e_4 * fs_60_221_65 * r_2 * h7_p6 - e_4 * fs_9_26_66 * r_4 * h5_p4 + e_5 * fs_21_8398_390 * h11_p4 + e_5 * fs_6_4199_7735 * h11_p6 - e_5 * fs_7_442_78 * r_2 * h9_p4 - e_5 * fs_3_221_91 * r_2 * h9_p6 + e_5 * fs_75_4199_10 * r_4 * h7_p4 - e_5 * fs_60_4199_65 * r_4 * h7_p6 + e_5 * fs_3_221_66 * r_6 * h5_p4;
    }

    // NOTE: the angular components are formed in 124 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph3_p3, ph5_p3, ph5_p5, ph7_p3, ph7_p5, ph9_p3, ph9_p5, ph11_p3, ph11_p5, ab_2, pc_54 : simd::cache_line_size())
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

        const auto h3_p3 = ph3_p3[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h5_p5 = ph5_p5[k];
        const auto h7_p3 = ph7_p3[k];
        const auto h7_p5 = ph7_p5[k];
        const auto h9_p3 = ph9_p3[k];
        const auto h9_p5 = ph9_p5[k];
        const auto h11_p3 = ph11_p3[k];
        const auto h11_p5 = ph11_p5[k];

        pc_54[k] = e_1 * fs_315_16_42 * h3_p3 - e_2 * fs_65_8_6 * h5_p3 + e_2 * fs_75_8_30 * h5_p5 - e_2 * fs_35_2_42 * r_2 * h3_p3 - e_3 * fs_75_44_5 * h7_p3 - e_3 * fs_555_572_55 * h7_p5 + e_3 * fs_15_4_6 * r_2 * h5_p3 - e_3 * fs_225_52_30 * r_2 * h5_p5 + e_3 * fs_105_22_42 * r_4 * h3_p3 + e_4 * fs_3087_9724_22 * h9_p3 + e_4 * fs_21_9724_30030 * h9_p5 + e_4 * fs_75_187_5 * r_2 * h7_p3 + e_4 * fs_555_2431_55 * r_2 * h7_p5 - e_4 * fs_1_2_6 * r_4 * h5_p3 + e_4 * fs_15_26_30 * r_4 * h5_p5 - e_4 * fs_70_143_42 * r_6 * h3_p3 + e_5 * fs_21_4199_286 * h11_p3 + e_5 * fs_6_4199_15015 * h11_p5 - e_5 * fs_147_4862_22 * r_2 * h9_p3 - e_5 * fs_1_4862_30030 * r_2 * h9_p5 - e_5 * fs_75_3553_5 * r_4 * h7_p3 - e_5 * fs_555_46189_55 * r_4 * h7_p5 + e_5 * fs_1_51_6 * r_6 * h5_p3 - e_5 * fs_5_221_30 * r_6 * h5_p5 + e_5 * fs_7_429_42 * r_8 * h3_p3;
    }

    // NOTE: the angular components are formed in 124 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_2, pe_3, pe_4, pe_5, ph5_p2, ph5_p4, ph7_p2, ph7_p4, ph9_p2, ph9_p4, ph11_p2, ph11_p4, ab_2, pc_55 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h5_p2 = ph5_p2[k];
        const auto h5_p4 = ph5_p4[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h7_p4 = ph7_p4[k];
        const auto h9_p2 = ph9_p2[k];
        const auto h9_p4 = ph9_p4[k];
        const auto h11_p2 = ph11_p2[k];
        const auto h11_p4 = ph11_p4[k];

        pc_55[k] = e_2 * fs_15_4_30 * h5_p2 + e_2 * fs_15_8_10 * h5_p4 - e_3 * fs_60_13_3 * h7_p2 + e_3 * fs_45_572_66 * h7_p4 - e_3 * fs_45_26_30 * r_2 * h5_p2 - e_3 * fs_45_52_10 * r_2 * h5_p4 + e_4 * fs_21_374_385 * h9_p2 - e_4 * fs_147_9724_1430 * h9_p4 + e_4 * fs_240_221_3 * r_2 * h7_p2 - e_4 * fs_45_2431_66 * r_2 * h7_p4 + e_4 * fs_3_13_30 * r_4 * h5_p2 + e_4 * fs_3_26_10 * r_4 * h5_p4 + e_5 * fs_3_4199_30030 * h11_p2 + e_5 * fs_105_8398_286 * h11_p4 - e_5 * fs_1_187_385 * r_2 * h9_p2 + e_5 * fs_7_4862_1430 * r_2 * h9_p4 - e_5 * fs_240_4199_3 * r_4 * h7_p2 + e_5 * fs_45_46189_66 * r_4 * h7_p4 - e_5 * fs_2_221_30 * r_6 * h5_p2 - e_5 * fs_1_221_10 * r_6 * h5_p4;
    }

    // NOTE: the angular components are formed in 124 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph1_p1, ph3_p1, ph3_p3, ph5_p1, ph5_p3, ph7_p1, ph7_p3, ph9_p1, ph9_p3, ph11_p1, ph11_p3, ab_2, pc_56 : simd::cache_line_size())
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

        const auto h1_p1 = ph1_p1[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h3_p3 = ph3_p3[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h7_p3 = ph7_p3[k];
        const auto h9_p1 = ph9_p1[k];
        const auto h9_p3 = ph9_p3[k];
        const auto h11_p1 = ph11_p1[k];
        const auto h11_p3 = ph11_p3[k];

        pc_56[k] = e_0 * fs_945_32_14 * h1_p1 - e_1 * fs_315_16_21 * h3_p1 + e_1 * fs_315_16_35 * h3_p3 - e_1 * fs_945_16_14 * r_2 * h1_p1 + e_2 * fs_5_2_210 * h5_p1 - e_2 * fs_55_4_5 * h5_p3 + e_2 * fs_35_2_21 * r_2 * h3_p1 - e_2 * fs_35_2_35 * r_2 * h3_p3 + e_2 * fs_135_4_14 * r_4 * h1_p1 - e_3 * fs_5505_1144_2 * h7_p1 + e_3 * fs_3105_1144_6 * h7_p3 - e_3 * fs_15_13_210 * r_2 * h5_p1 + e_3 * fs_165_26_5 * r_2 * h5_p3 - e_3 * fs_105_22_21 * r_4 * h3_p1 + e_3 * fs_105_22_35 * r_4 * h3_p3 - e_3 * fs_15_2_14 * r_6 * h1_p1 + e_4 * fs_126_2431_70 * h9_p1 - e_4 * fs_441_4862_165 * h9_p3 + e_4 * fs_5505_4862_2 * r_2 * h7_p1 - e_4 * fs_3105_4862_6 * r_2 * h7_p3 + e_4 * fs_2_13_210 * r_4 * h5_p1 - e_4 * fs_11_13_5 * r_4 * h5_p3 + e_4 * fs_70_143_21 * r_6 * h3_p1 - e_4 * fs_70_143_35 * r_6 * h3_p3 + e_4 * fs_15_22_14 * r_8 * h1_p1 + e_5 * fs_45_4199_231 * h11_p1 + e_5 * fs_21_4199_2145 * h11_p3 - e_5 * fs_12_2431_70 * r_2 * h9_p1 + e_5 * fs_21_2431_165 * r_2 * h9_p3 - e_5 * fs_5505_92378_2 * r_4 * h7_p1 + e_5 * fs_3105_92378_6 * r_4 * h7_p3 - e_5 * fs_4_663_210 * r_6 * h5_p1 + e_5 * fs_22_663_5 * r_6 * h5_p3 - e_5 * fs_7_429_21 * r_8 * h3_p1 + e_5 * fs_7_429_35 * r_8 * h3_p3 - e_5 * fs_3_143_14 * r_10 * h1_p1;
    }

    // NOTE: the angular components are formed in 124 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph1_0, ph3_0, ph3_p2, ph5_0, ph5_p2, ph7_0, ph7_p2, ph9_0, ph9_p2, ph11_0, ph11_p2, ab_2, pc_57 : simd::cache_line_size())
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

        const auto h1_0 = ph1_0[k];
        const auto h3_0 = ph3_0[k];
        const auto h3_p2 = ph3_p2[k];
        const auto h5_0 = ph5_0[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h7_0 = ph7_0[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h9_0 = ph9_0[k];
        const auto h9_p2 = ph9_p2[k];
        const auto h11_0 = ph11_0[k];
        const auto h11_p2 = ph11_p2[k];

        pc_57[k] = e_0 * fs_945_32_35 * h1_0 - e_1 * fs_315_16_35 * h3_0 + e_1 * fs_315_16_21 * h3_p2 - e_1 * fs_945_16_35 * r_2 * h1_0 + e_2 * fs_5_1_35 * h5_0 - e_2 * fs_20_1_3 * h5_p2 + e_2 * fs_35_2_35 * r_2 * h3_0 - e_2 * fs_35_2_21 * r_2 * h3_p2 + e_2 * fs_135_4_35 * r_4 * h1_0 - e_3 * fs_60_143_35 * h7_0 + e_3 * fs_405_286_30 * h7_p2 - e_3 * fs_30_13_35 * r_2 * h5_0 + e_3 * fs_120_13_3 * r_2 * h5_p2 - e_3 * fs_105_22_35 * r_4 * h3_0 + e_3 * fs_105_22_21 * r_4 * h3_p2 - e_3 * fs_15_2_35 * r_6 * h1_0 - e_4 * fs_189_2431_35 * h9_0 - e_4 * fs_252_2431_154 * h9_p2 + e_4 * fs_240_2431_35 * r_2 * h7_0 - e_4 * fs_810_2431_30 * r_2 * h7_p2 + e_4 * fs_4_13_35 * r_4 * h5_0 - e_4 * fs_16_13_3 * r_4 * h5_p2 + e_4 * fs_70_143_35 * r_6 * h3_0 - e_4 * fs_70_143_21 * r_6 * h3_p2 + e_4 * fs_15_22_35 * r_8 * h1_0 + e_5 * fs_198_4199_35 * h11_0 + e_5 * fs_18_4199_3003 * h11_p2 + e_5 * fs_18_2431_35 * r_2 * h9_0 + e_5 * fs_24_2431_154 * r_2 * h9_p2 - e_5 * fs_240_46189_35 * r_4 * h7_0 + e_5 * fs_810_46189_30 * r_4 * h7_p2 - e_5 * fs_8_663_35 * r_6 * h5_0 + e_5 * fs_32_663_3 * r_6 * h5_p2 - e_5 * fs_7_429_35 * r_8 * h3_0 + e_5 * fs_7_429_21 * r_8 * h3_p2 - e_5 * fs_3_143_35 * r_10 * h1_0;
    }

    // NOTE: the angular components are formed in 124 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph1_m1, ph3_m2, ph5_m2, ph5_m1, ph7_m2, ph7_m1, ph9_m2, ph9_m1, ph11_m2, ph11_m1, ab_2, pc_58, pc_59 : simd::cache_line_size())
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

        const auto h1_m1 = ph1_m1[k];
        const auto h3_m2 = ph3_m2[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h7_m1 = ph7_m1[k];
        const auto h9_m2 = ph9_m2[k];
        const auto h9_m1 = ph9_m1[k];
        const auto h11_m2 = ph11_m2[k];
        const auto h11_m1 = ph11_m1[k];

        pc_58[k] = - e_0 * fs_945_32_15 * h1_m1 + e_1 * fs_945_16_15 * r_2 * h1_m1 + e_2 * f_15_1 * h5_m1 - e_2 * fs_135_4_15 * r_4 * h1_m1 - e_3 * fs_15_26_105 * h7_m1 - e_3 * f_90_13 * r_2 * h5_m1 + e_3 * fs_15_2_15 * r_6 * h1_m1 + e_4 * fs_147_187_3 * h9_m1 + e_4 * fs_30_221_105 * r_2 * h7_m1 + e_4 * f_12_13 * r_4 * h5_m1 - e_4 * fs_15_22_15 * r_8 * h1_m1 - e_5 * fs_126_4199_110 * h11_m1 - e_5 * fs_14_187_3 * r_2 * h9_m1 - e_5 * fs_30_4199_105 * r_4 * h7_m1 - e_5 * f_8_221 * r_6 * h5_m1 + e_5 * fs_3_143_15 * r_10 * h1_m1;

        pc_59[k] = - e_1 * fs_315_16_21 * h3_m2 + e_2 * fs_20_1_3 * h5_m2 + e_2 * fs_35_2_21 * r_2 * h3_m2 - e_3 * fs_405_286_30 * h7_m2 - e_3 * fs_120_13_3 * r_2 * h5_m2 - e_3 * fs_105_22_21 * r_4 * h3_m2 + e_4 * fs_252_2431_154 * h9_m2 + e_4 * fs_810_2431_30 * r_2 * h7_m2 + e_4 * fs_16_13_3 * r_4 * h5_m2 + e_4 * fs_70_143_21 * r_6 * h3_m2 - e_5 * fs_18_4199_3003 * h11_m2 - e_5 * fs_24_2431_154 * r_2 * h9_m2 - e_5 * fs_810_46189_30 * r_4 * h7_m2 - e_5 * fs_32_663_3 * r_6 * h5_m2 - e_5 * fs_7_429_21 * r_8 * h3_m2;
    }

    // NOTE: the angular components are formed in 124 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph1_m1, ph3_m3, ph3_m1, ph5_m3, ph5_m1, ph7_m3, ph7_m1, ph9_m3, ph9_m1, ph11_m3, ph11_m1, ab_2, pc_60 : simd::cache_line_size())
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

        const auto h1_m1 = ph1_m1[k];
        const auto h3_m3 = ph3_m3[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h5_m3 = ph5_m3[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h7_m3 = ph7_m3[k];
        const auto h7_m1 = ph7_m1[k];
        const auto h9_m3 = ph9_m3[k];
        const auto h9_m1 = ph9_m1[k];
        const auto h11_m3 = ph11_m3[k];
        const auto h11_m1 = ph11_m1[k];

        pc_60[k] = - e_0 * fs_945_32_14 * h1_m1 - e_1 * fs_315_16_35 * h3_m3 + e_1 * fs_315_16_21 * h3_m1 + e_1 * fs_945_16_14 * r_2 * h1_m1 + e_2 * fs_55_4_5 * h5_m3 - e_2 * fs_5_2_210 * h5_m1 + e_2 * fs_35_2_35 * r_2 * h3_m3 - e_2 * fs_35_2_21 * r_2 * h3_m1 - e_2 * fs_135_4_14 * r_4 * h1_m1 - e_3 * fs_3105_1144_6 * h7_m3 + e_3 * fs_5505_1144_2 * h7_m1 - e_3 * fs_165_26_5 * r_2 * h5_m3 + e_3 * fs_15_13_210 * r_2 * h5_m1 - e_3 * fs_105_22_35 * r_4 * h3_m3 + e_3 * fs_105_22_21 * r_4 * h3_m1 + e_3 * fs_15_2_14 * r_6 * h1_m1 + e_4 * fs_441_4862_165 * h9_m3 - e_4 * fs_126_2431_70 * h9_m1 + e_4 * fs_3105_4862_6 * r_2 * h7_m3 - e_4 * fs_5505_4862_2 * r_2 * h7_m1 + e_4 * fs_11_13_5 * r_4 * h5_m3 - e_4 * fs_2_13_210 * r_4 * h5_m1 + e_4 * fs_70_143_35 * r_6 * h3_m3 - e_4 * fs_70_143_21 * r_6 * h3_m1 - e_4 * fs_15_22_14 * r_8 * h1_m1 - e_5 * fs_21_4199_2145 * h11_m3 - e_5 * fs_45_4199_231 * h11_m1 - e_5 * fs_21_2431_165 * r_2 * h9_m3 + e_5 * fs_12_2431_70 * r_2 * h9_m1 - e_5 * fs_3105_92378_6 * r_4 * h7_m3 + e_5 * fs_5505_92378_2 * r_4 * h7_m1 - e_5 * fs_22_663_5 * r_6 * h5_m3 + e_5 * fs_4_663_210 * r_6 * h5_m1 - e_5 * fs_7_429_35 * r_8 * h3_m3 + e_5 * fs_7_429_21 * r_8 * h3_m1 + e_5 * fs_3_143_14 * r_10 * h1_m1;
    }

    // NOTE: the angular components are formed in 124 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_2, pe_3, pe_4, pe_5, ph5_m4, ph5_m2, ph7_m4, ph7_m2, ph9_m4, ph9_m2, ph11_m4, ph11_m2, ab_2, pc_61 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h5_m4 = ph5_m4[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h7_m4 = ph7_m4[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h9_m4 = ph9_m4[k];
        const auto h9_m2 = ph9_m2[k];
        const auto h11_m4 = ph11_m4[k];
        const auto h11_m2 = ph11_m2[k];

        pc_61[k] = - e_2 * fs_15_8_10 * h5_m4 - e_2 * fs_15_4_30 * h5_m2 - e_3 * fs_45_572_66 * h7_m4 + e_3 * fs_60_13_3 * h7_m2 + e_3 * fs_45_52_10 * r_2 * h5_m4 + e_3 * fs_45_26_30 * r_2 * h5_m2 + e_4 * fs_147_9724_1430 * h9_m4 - e_4 * fs_21_374_385 * h9_m2 + e_4 * fs_45_2431_66 * r_2 * h7_m4 - e_4 * fs_240_221_3 * r_2 * h7_m2 - e_4 * fs_3_26_10 * r_4 * h5_m4 - e_4 * fs_3_13_30 * r_4 * h5_m2 - e_5 * fs_105_8398_286 * h11_m4 - e_5 * fs_3_4199_30030 * h11_m2 - e_5 * fs_7_4862_1430 * r_2 * h9_m4 + e_5 * fs_1_187_385 * r_2 * h9_m2 - e_5 * fs_45_46189_66 * r_4 * h7_m4 + e_5 * fs_240_4199_3 * r_4 * h7_m2 + e_5 * fs_1_221_10 * r_6 * h5_m4 + e_5 * fs_2_221_30 * r_6 * h5_m2;
    }

    // NOTE: the angular components are formed in 124 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph3_m3, ph5_m5, ph5_m3, ph7_m5, ph7_m3, ph9_m5, ph9_m3, ph11_m5, ph11_m3, ab_2, pc_62 : simd::cache_line_size())
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

        const auto h3_m3 = ph3_m3[k];
        const auto h5_m5 = ph5_m5[k];
        const auto h5_m3 = ph5_m3[k];
        const auto h7_m5 = ph7_m5[k];
        const auto h7_m3 = ph7_m3[k];
        const auto h9_m5 = ph9_m5[k];
        const auto h9_m3 = ph9_m3[k];
        const auto h11_m5 = ph11_m5[k];
        const auto h11_m3 = ph11_m3[k];

        pc_62[k] = - e_1 * fs_315_16_42 * h3_m3 - e_2 * fs_75_8_30 * h5_m5 + e_2 * fs_65_8_6 * h5_m3 + e_2 * fs_35_2_42 * r_2 * h3_m3 + e_3 * fs_555_572_55 * h7_m5 + e_3 * fs_75_44_5 * h7_m3 + e_3 * fs_225_52_30 * r_2 * h5_m5 - e_3 * fs_15_4_6 * r_2 * h5_m3 - e_3 * fs_105_22_42 * r_4 * h3_m3 - e_4 * fs_21_9724_30030 * h9_m5 - e_4 * fs_3087_9724_22 * h9_m3 - e_4 * fs_555_2431_55 * r_2 * h7_m5 - e_4 * fs_75_187_5 * r_2 * h7_m3 - e_4 * fs_15_26_30 * r_4 * h5_m5 + e_4 * fs_1_2_6 * r_4 * h5_m3 + e_4 * fs_70_143_42 * r_6 * h3_m3 - e_5 * fs_6_4199_15015 * h11_m5 - e_5 * fs_21_4199_286 * h11_m3 + e_5 * fs_1_4862_30030 * r_2 * h9_m5 + e_5 * fs_147_4862_22 * r_2 * h9_m3 + e_5 * fs_555_46189_55 * r_4 * h7_m5 + e_5 * fs_75_3553_5 * r_4 * h7_m3 + e_5 * fs_5_221_30 * r_6 * h5_m5 - e_5 * fs_1_51_6 * r_6 * h5_m3 - e_5 * fs_7_429_42 * r_8 * h3_m3;
    }

    // NOTE: the angular components are formed in 124 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_2, pe_3, pe_4, pe_5, ph5_m4, ph7_m6, ph7_m4, ph9_m6, ph9_m4, ph11_m6, ph11_m4, ab_2, pc_63 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h5_m4 = ph5_m4[k];
        const auto h7_m6 = ph7_m6[k];
        const auto h7_m4 = ph7_m4[k];
        const auto h9_m6 = ph9_m6[k];
        const auto h9_m4 = ph9_m4[k];
        const auto h11_m6 = ph11_m6[k];
        const auto h11_m4 = ph11_m4[k];

        pc_63[k] = e_2 * fs_45_8_66 * h5_m4 + e_3 * fs_15_13_65 * h7_m6 - e_3 * fs_75_52_10 * h7_m4 - e_3 * fs_135_52_66 * r_2 * h5_m4 - e_4 * fs_63_442_91 * h9_m6 - e_4 * fs_147_884_78 * h9_m4 - e_4 * fs_60_221_65 * r_2 * h7_m6 + e_4 * fs_75_221_10 * r_2 * h7_m4 + e_4 * fs_9_26_66 * r_4 * h5_m4 - e_5 * fs_6_4199_7735 * h11_m6 - e_5 * fs_21_8398_390 * h11_m4 + e_5 * fs_3_221_91 * r_2 * h9_m6 + e_5 * fs_7_442_78 * r_2 * h9_m4 + e_5 * fs_60_4199_65 * r_4 * h7_m6 - e_5 * fs_75_4199_10 * r_4 * h7_m4 - e_5 * fs_3_221_66 * r_6 * h5_m4;
    }

    // NOTE: the angular components are formed in 124 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_2, pe_3, pe_4, pe_5, ph5_m5, ph7_m7, ph7_m6, ph7_m5, ph9_m7, ph9_m6, ph9_m5, ph11_m7, ph11_m6, ph11_m5, ab_2, pc_64, pc_65, pc_66 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h5_m5 = ph5_m5[k];
        const auto h7_m7 = ph7_m7[k];
        const auto h7_m6 = ph7_m6[k];
        const auto h7_m5 = ph7_m5[k];
        const auto h9_m7 = ph9_m7[k];
        const auto h9_m6 = ph9_m6[k];
        const auto h9_m5 = ph9_m5[k];
        const auto h11_m7 = ph11_m7[k];
        const auto h11_m6 = ph11_m6[k];
        const auto h11_m5 = ph11_m5[k];

        pc_64[k] = - e_2 * fs_15_4_55 * h5_m5 - e_3 * fs_15_104_2730 * h7_m7 - e_3 * fs_225_104_30 * h7_m5 + e_3 * fs_45_26_55 * r_2 * h5_m5 - e_4 * fs_42_221_91 * h9_m7 - e_4 * fs_21_442_455 * h9_m5 + e_4 * fs_15_442_2730 * r_2 * h7_m7 + e_4 * fs_225_442_30 * r_2 * h7_m5 - e_4 * fs_3_13_55 * r_4 * h5_m5 - e_5 * fs_3_4199_9282 * h11_m7 - e_5 * fs_3_4199_910 * h11_m5 + e_5 * fs_4_221_91 * r_2 * h9_m7 + e_5 * fs_1_221_455 * r_2 * h9_m5 - e_5 * fs_15_8398_2730 * r_4 * h7_m7 - e_5 * fs_225_8398_30 * r_4 * h7_m5 + e_5 * fs_2_221_55 * r_6 * h5_m5;

        pc_65[k] = - e_3 * fs_225_52_13 * h7_m6 - e_4 * fs_21_221_455 * h9_m6 + e_4 * fs_225_221_13 * r_2 * h7_m6 - e_5 * fs_6_4199_1547 * h11_m6 + e_5 * fs_2_221_455 * r_2 * h9_m6 - e_5 * fs_225_4199_13 * r_4 * h7_m6;

        pc_66[k] = e_2 * fs_75_4_11 * h5_m5 + e_3 * fs_75_52_6 * h7_m5 - e_3 * fs_225_26_11 * r_2 * h5_m5 - e_4 * fs_105_442_91 * h9_m5 - e_4 * fs_75_221_6 * r_2 * h7_m5 + e_4 * fs_15_13_11 * r_4 * h5_m5 - e_5 * fs_36_4199_182 * h11_m5 + e_5 * fs_5_221_91 * r_2 * h9_m5 + e_5 * fs_75_4199_6 * r_4 * h7_m5 - e_5 * fs_10_221_11 * r_6 * h5_m5;
    }

    // NOTE: the angular components are formed in 124 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph3_m3, ph5_m4, ph5_m3, ph7_m4, ph7_m3, ph9_m4, ph9_m3, ph11_m4, ph11_m3, ab_2, pc_67, pc_68 : simd::cache_line_size())
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

        const auto h3_m3 = ph3_m3[k];
        const auto h5_m4 = ph5_m4[k];
        const auto h5_m3 = ph5_m3[k];
        const auto h7_m4 = ph7_m4[k];
        const auto h7_m3 = ph7_m3[k];
        const auto h9_m4 = ph9_m4[k];
        const auto h9_m3 = ph9_m3[k];
        const auto h11_m4 = ph11_m4[k];
        const auto h11_m3 = ph11_m3[k];

        pc_67[k] = - e_2 * fs_15_1_5 * h5_m4 + e_3 * fs_300_143_33 * h7_m4 + e_3 * fs_90_13_5 * r_2 * h5_m4 - e_4 * fs_147_2431_715 * h9_m4 - e_4 * fs_1200_2431_33 * r_2 * h7_m4 - e_4 * fs_12_13_5 * r_4 * h5_m4 - e_5 * fs_63_4199_143 * h11_m4 + e_5 * fs_14_2431_715 * r_2 * h9_m4 + e_5 * fs_1200_46189_33 * r_4 * h7_m4 + e_5 * fs_8_221_5 * r_6 * h5_m4;

        pc_68[k] = e_1 * fs_315_8_21 * h3_m3 - e_2 * fs_125_4_3 * h5_m3 - e_2 * fs_35_1_21 * r_2 * h3_m3 + e_3 * fs_1665_572_10 * h7_m3 + e_3 * fs_375_26_3 * r_2 * h5_m3 + e_3 * fs_105_11_21 * r_4 * h3_m3 - e_4 * fs_735_4862_11 * h9_m3 - e_4 * fs_1665_2431_10 * r_2 * h7_m3 - e_4 * fs_25_13_3 * r_4 * h5_m3 - e_4 * fs_140_143_21 * r_6 * h3_m3 - e_5 * fs_84_4199_143 * h11_m3 + e_5 * fs_35_2431_11 * r_2 * h9_m3 + e_5 * fs_1665_46189_10 * r_4 * h7_m3 + e_5 * fs_50_663_3 * r_6 * h5_m3 + e_5 * fs_14_429_21 * r_8 * h3_m3;
    }

    // NOTE: the angular components are formed in 124 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph1_m1, ph3_m2, ph3_m1, ph5_m2, ph5_m1, ph7_m2, ph7_m1, ph9_m2, ph9_m1, ph11_m2, ph11_m1, ab_2, pc_69, pc_70 : simd::cache_line_size())
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

        const auto h1_m1 = ph1_m1[k];
        const auto h3_m2 = ph3_m2[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h7_m1 = ph7_m1[k];
        const auto h9_m2 = ph9_m2[k];
        const auto h9_m1 = ph9_m1[k];
        const auto h11_m2 = ph11_m2[k];
        const auto h11_m1 = ph11_m1[k];

        pc_69[k] = e_1 * fs_315_16_14 * h3_m2 - e_2 * fs_25_2_2 * h5_m2 - e_2 * fs_35_2_14 * r_2 * h3_m2 + e_3 * fs_135_572_5 * h7_m2 + e_3 * fs_75_13_2 * r_2 * h5_m2 + e_3 * fs_105_22_14 * r_4 * h3_m2 + e_4 * fs_105_2431_231 * h9_m2 - e_4 * fs_135_2431_5 * r_2 * h7_m2 - e_4 * fs_10_13_2 * r_4 * h5_m2 - e_4 * fs_70_143_14 * r_6 * h3_m2 - e_5 * fs_27_4199_2002 * h11_m2 - e_5 * fs_10_2431_231 * r_2 * h9_m2 + e_5 * fs_135_46189_5 * r_4 * h7_m2 + e_5 * fs_20_663_2 * r_6 * h5_m2 + e_5 * fs_7_429_14 * r_8 * h3_m2;

        pc_70[k] = e_0 * fs_945_32_21 * h1_m1 - e_1 * fs_315_16_14 * h3_m1 - e_1 * fs_945_16_21 * r_2 * h1_m1 + e_2 * fs_5_1_35 * h5_m1 + e_2 * fs_35_2_14 * r_2 * h3_m1 + e_2 * fs_135_4_21 * r_4 * h1_m1 - e_3 * fs_1275_286_3 * h7_m1 - e_3 * fs_30_13_35 * r_2 * h5_m1 - e_3 * fs_105_22_14 * r_4 * h3_m1 - e_3 * fs_15_2_21 * r_6 * h1_m1 + e_4 * fs_21_143_105 * h9_m1 + e_4 * fs_150_143_3 * r_2 * h7_m1 + e_4 * fs_4_13_35 * r_4 * h5_m1 + e_4 * fs_70_143_14 * r_6 * h3_m1 + e_4 * fs_15_22_21 * r_8 * h1_m1 - e_5 * fs_108_4199_154 * h11_m1 - e_5 * fs_2_143_105 * r_2 * h9_m1 - e_5 * fs_150_2717_3 * r_4 * h7_m1 - e_5 * fs_8_663_35 * r_6 * h5_m1 - e_5 * fs_7_429_14 * r_8 * h3_m1 - e_5 * fs_3_143_21 * r_10 * h1_m1;
    }

    // NOTE: the angular components are formed in 124 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph1_0, ph3_0, ph5_0, ph7_0, ph9_0, ph11_0, ab_2, pc_71 : simd::cache_line_size())
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

        const auto h1_0 = ph1_0[k];
        const auto h3_0 = ph3_0[k];
        const auto h5_0 = ph5_0[k];
        const auto h7_0 = ph7_0[k];
        const auto h9_0 = ph9_0[k];
        const auto h11_0 = ph11_0[k];

        pc_71[k] = e_0 * f_2835_16 * h1_0 - e_1 * f_2205_16 * h3_0 - e_1 * f_2835_8 * r_2 * h1_0 + e_2 * f_50_1 * h5_0 + e_2 * f_245_2 * r_2 * h3_0 + e_2 * f_405_2 * r_4 * h1_0 - e_3 * f_1575_143 * h7_0 - e_3 * f_300_13 * r_2 * h5_0 - e_3 * f_735_22 * r_4 * h3_0 - e_3 * f_45_1 * r_6 * h1_0 + e_4 * f_4410_2431 * h9_0 + e_4 * f_6300_2431 * r_2 * h7_0 + e_4 * f_40_13 * r_4 * h5_0 + e_4 * f_490_143 * r_6 * h3_0 + e_4 * f_45_11 * r_8 * h1_0 - e_5 * f_1386_4199 * h11_0 - e_5 * f_420_2431 * r_2 * h9_0 - e_5 * f_6300_46189 * r_4 * h7_0 - e_5 * f_80_663 * r_6 * h5_0 - e_5 * f_49_429 * r_8 * h3_0 - e_5 * f_18_143 * r_10 * h1_0;
    }

    // NOTE: the angular components are formed in 124 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph1_p1, ph3_p1, ph3_p2, ph5_p1, ph5_p2, ph7_p1, ph7_p2, ph9_p1, ph9_p2, ph11_p1, ph11_p2, ab_2, pc_72, pc_73 : simd::cache_line_size())
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

        const auto h1_p1 = ph1_p1[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h3_p2 = ph3_p2[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h9_p1 = ph9_p1[k];
        const auto h9_p2 = ph9_p2[k];
        const auto h11_p1 = ph11_p1[k];
        const auto h11_p2 = ph11_p2[k];

        pc_72[k] = e_0 * fs_945_32_21 * h1_p1 - e_1 * fs_315_16_14 * h3_p1 - e_1 * fs_945_16_21 * r_2 * h1_p1 + e_2 * fs_5_1_35 * h5_p1 + e_2 * fs_35_2_14 * r_2 * h3_p1 + e_2 * fs_135_4_21 * r_4 * h1_p1 - e_3 * fs_1275_286_3 * h7_p1 - e_3 * fs_30_13_35 * r_2 * h5_p1 - e_3 * fs_105_22_14 * r_4 * h3_p1 - e_3 * fs_15_2_21 * r_6 * h1_p1 + e_4 * fs_21_143_105 * h9_p1 + e_4 * fs_150_143_3 * r_2 * h7_p1 + e_4 * fs_4_13_35 * r_4 * h5_p1 + e_4 * fs_70_143_14 * r_6 * h3_p1 + e_4 * fs_15_22_21 * r_8 * h1_p1 - e_5 * fs_108_4199_154 * h11_p1 - e_5 * fs_2_143_105 * r_2 * h9_p1 - e_5 * fs_150_2717_3 * r_4 * h7_p1 - e_5 * fs_8_663_35 * r_6 * h5_p1 - e_5 * fs_7_429_14 * r_8 * h3_p1 - e_5 * fs_3_143_21 * r_10 * h1_p1;

        pc_73[k] = e_1 * fs_315_16_14 * h3_p2 - e_2 * fs_25_2_2 * h5_p2 - e_2 * fs_35_2_14 * r_2 * h3_p2 + e_3 * fs_135_572_5 * h7_p2 + e_3 * fs_75_13_2 * r_2 * h5_p2 + e_3 * fs_105_22_14 * r_4 * h3_p2 + e_4 * fs_105_2431_231 * h9_p2 - e_4 * fs_135_2431_5 * r_2 * h7_p2 - e_4 * fs_10_13_2 * r_4 * h5_p2 - e_4 * fs_70_143_14 * r_6 * h3_p2 - e_5 * fs_27_4199_2002 * h11_p2 - e_5 * fs_10_2431_231 * r_2 * h9_p2 + e_5 * fs_135_46189_5 * r_4 * h7_p2 + e_5 * fs_20_663_2 * r_6 * h5_p2 + e_5 * fs_7_429_14 * r_8 * h3_p2;
    }

    // NOTE: the angular components are formed in 124 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph3_p3, ph5_p3, ph5_p4, ph7_p3, ph7_p4, ph9_p3, ph9_p4, ph11_p3, ph11_p4, ab_2, pc_74, pc_75 : simd::cache_line_size())
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

        const auto h3_p3 = ph3_p3[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h5_p4 = ph5_p4[k];
        const auto h7_p3 = ph7_p3[k];
        const auto h7_p4 = ph7_p4[k];
        const auto h9_p3 = ph9_p3[k];
        const auto h9_p4 = ph9_p4[k];
        const auto h11_p3 = ph11_p3[k];
        const auto h11_p4 = ph11_p4[k];

        pc_74[k] = e_1 * fs_315_8_21 * h3_p3 - e_2 * fs_125_4_3 * h5_p3 - e_2 * fs_35_1_21 * r_2 * h3_p3 + e_3 * fs_1665_572_10 * h7_p3 + e_3 * fs_375_26_3 * r_2 * h5_p3 + e_3 * fs_105_11_21 * r_4 * h3_p3 - e_4 * fs_735_4862_11 * h9_p3 - e_4 * fs_1665_2431_10 * r_2 * h7_p3 - e_4 * fs_25_13_3 * r_4 * h5_p3 - e_4 * fs_140_143_21 * r_6 * h3_p3 - e_5 * fs_84_4199_143 * h11_p3 + e_5 * fs_35_2431_11 * r_2 * h9_p3 + e_5 * fs_1665_46189_10 * r_4 * h7_p3 + e_5 * fs_50_663_3 * r_6 * h5_p3 + e_5 * fs_14_429_21 * r_8 * h3_p3;

        pc_75[k] = - e_2 * fs_15_1_5 * h5_p4 + e_3 * fs_300_143_33 * h7_p4 + e_3 * fs_90_13_5 * r_2 * h5_p4 - e_4 * fs_147_2431_715 * h9_p4 - e_4 * fs_1200_2431_33 * r_2 * h7_p4 - e_4 * fs_12_13_5 * r_4 * h5_p4 - e_5 * fs_63_4199_143 * h11_p4 + e_5 * fs_14_2431_715 * r_2 * h9_p4 + e_5 * fs_1200_46189_33 * r_4 * h7_p4 + e_5 * fs_8_221_5 * r_6 * h5_p4;
    }

    // NOTE: the angular components are formed in 124 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_2, pe_3, pe_4, pe_5, ph5_p5, ph7_p5, ph7_p6, ph9_p5, ph9_p6, ph11_p5, ph11_p6, ab_2, pc_76, pc_77 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h5_p5 = ph5_p5[k];
        const auto h7_p5 = ph7_p5[k];
        const auto h7_p6 = ph7_p6[k];
        const auto h9_p5 = ph9_p5[k];
        const auto h9_p6 = ph9_p6[k];
        const auto h11_p5 = ph11_p5[k];
        const auto h11_p6 = ph11_p6[k];

        pc_76[k] = e_2 * fs_75_4_11 * h5_p5 + e_3 * fs_75_52_6 * h7_p5 - e_3 * fs_225_26_11 * r_2 * h5_p5 - e_4 * fs_105_442_91 * h9_p5 - e_4 * fs_75_221_6 * r_2 * h7_p5 + e_4 * fs_15_13_11 * r_4 * h5_p5 - e_5 * fs_36_4199_182 * h11_p5 + e_5 * fs_5_221_91 * r_2 * h9_p5 + e_5 * fs_75_4199_6 * r_4 * h7_p5 - e_5 * fs_10_221_11 * r_6 * h5_p5;

        pc_77[k] = - e_3 * fs_225_52_13 * h7_p6 - e_4 * fs_21_221_455 * h9_p6 + e_4 * fs_225_221_13 * r_2 * h7_p6 - e_5 * fs_6_4199_1547 * h11_p6 + e_5 * fs_2_221_455 * r_2 * h9_p6 - e_5 * fs_225_4199_13 * r_4 * h7_p6;
    }

    // NOTE: the angular components are formed in 124 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_2, pe_3, pe_4, pe_5, ph5_m5, ph7_m7, ph7_m5, ph9_m7, ph9_m5, ph11_m7, ph11_m5, ab_2, pc_78 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h5_m5 = ph5_m5[k];
        const auto h7_m7 = ph7_m7[k];
        const auto h7_m5 = ph7_m5[k];
        const auto h9_m7 = ph9_m7[k];
        const auto h9_m5 = ph9_m5[k];
        const auto h11_m7 = ph11_m7[k];
        const auto h11_m5 = ph11_m5[k];

        pc_78[k] = e_2 * fs_15_4_55 * h5_m5 - e_3 * fs_15_104_2730 * h7_m7 + e_3 * fs_225_104_30 * h7_m5 - e_3 * fs_45_26_55 * r_2 * h5_m5 - e_4 * fs_42_221_91 * h9_m7 + e_4 * fs_21_442_455 * h9_m5 + e_4 * fs_15_442_2730 * r_2 * h7_m7 - e_4 * fs_225_442_30 * r_2 * h7_m5 + e_4 * fs_3_13_55 * r_4 * h5_m5 - e_5 * fs_3_4199_9282 * h11_m7 + e_5 * fs_3_4199_910 * h11_m5 + e_5 * fs_4_221_91 * r_2 * h9_m7 - e_5 * fs_1_221_455 * r_2 * h9_m5 - e_5 * fs_15_8398_2730 * r_4 * h7_m7 + e_5 * fs_225_8398_30 * r_4 * h7_m5 - e_5 * fs_2_221_55 * r_6 * h5_m5;
    }

    // NOTE: the angular components are formed in 124 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_2, pe_3, pe_4, pe_5, ph5_m4, ph7_m6, ph7_m4, ph9_m6, ph9_m4, ph11_m6, ph11_m4, ab_2, pc_79 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h5_m4 = ph5_m4[k];
        const auto h7_m6 = ph7_m6[k];
        const auto h7_m4 = ph7_m4[k];
        const auto h9_m6 = ph9_m6[k];
        const auto h9_m4 = ph9_m4[k];
        const auto h11_m6 = ph11_m6[k];
        const auto h11_m4 = ph11_m4[k];

        pc_79[k] = - e_2 * fs_45_8_66 * h5_m4 + e_3 * fs_15_13_65 * h7_m6 + e_3 * fs_75_52_10 * h7_m4 + e_3 * fs_135_52_66 * r_2 * h5_m4 - e_4 * fs_63_442_91 * h9_m6 + e_4 * fs_147_884_78 * h9_m4 - e_4 * fs_60_221_65 * r_2 * h7_m6 - e_4 * fs_75_221_10 * r_2 * h7_m4 - e_4 * fs_9_26_66 * r_4 * h5_m4 - e_5 * fs_6_4199_7735 * h11_m6 + e_5 * fs_21_8398_390 * h11_m4 + e_5 * fs_3_221_91 * r_2 * h9_m6 - e_5 * fs_7_442_78 * r_2 * h9_m4 + e_5 * fs_60_4199_65 * r_4 * h7_m6 + e_5 * fs_75_4199_10 * r_4 * h7_m4 + e_5 * fs_3_221_66 * r_6 * h5_m4;
    }

    // NOTE: the angular components are formed in 124 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph3_m3, ph5_m5, ph5_m3, ph7_m5, ph7_m3, ph9_m5, ph9_m3, ph11_m5, ph11_m3, ab_2, pc_80 : simd::cache_line_size())
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

        const auto h3_m3 = ph3_m3[k];
        const auto h5_m5 = ph5_m5[k];
        const auto h5_m3 = ph5_m3[k];
        const auto h7_m5 = ph7_m5[k];
        const auto h7_m3 = ph7_m3[k];
        const auto h9_m5 = ph9_m5[k];
        const auto h9_m3 = ph9_m3[k];
        const auto h11_m5 = ph11_m5[k];
        const auto h11_m3 = ph11_m3[k];

        pc_80[k] = e_1 * fs_315_16_42 * h3_m3 - e_2 * fs_75_8_30 * h5_m5 - e_2 * fs_65_8_6 * h5_m3 - e_2 * fs_35_2_42 * r_2 * h3_m3 + e_3 * fs_555_572_55 * h7_m5 - e_3 * fs_75_44_5 * h7_m3 + e_3 * fs_225_52_30 * r_2 * h5_m5 + e_3 * fs_15_4_6 * r_2 * h5_m3 + e_3 * fs_105_22_42 * r_4 * h3_m3 - e_4 * fs_21_9724_30030 * h9_m5 + e_4 * fs_3087_9724_22 * h9_m3 - e_4 * fs_555_2431_55 * r_2 * h7_m5 + e_4 * fs_75_187_5 * r_2 * h7_m3 - e_4 * fs_15_26_30 * r_4 * h5_m5 - e_4 * fs_1_2_6 * r_4 * h5_m3 - e_4 * fs_70_143_42 * r_6 * h3_m3 - e_5 * fs_6_4199_15015 * h11_m5 + e_5 * fs_21_4199_286 * h11_m3 + e_5 * fs_1_4862_30030 * r_2 * h9_m5 - e_5 * fs_147_4862_22 * r_2 * h9_m3 + e_5 * fs_555_46189_55 * r_4 * h7_m5 - e_5 * fs_75_3553_5 * r_4 * h7_m3 + e_5 * fs_5_221_30 * r_6 * h5_m5 + e_5 * fs_1_51_6 * r_6 * h5_m3 + e_5 * fs_7_429_42 * r_8 * h3_m3;
    }

    // NOTE: the angular components are formed in 124 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_2, pe_3, pe_4, pe_5, ph5_m4, ph5_m2, ph7_m4, ph7_m2, ph9_m4, ph9_m2, ph11_m4, ph11_m2, ab_2, pc_81 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h5_m4 = ph5_m4[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h7_m4 = ph7_m4[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h9_m4 = ph9_m4[k];
        const auto h9_m2 = ph9_m2[k];
        const auto h11_m4 = ph11_m4[k];
        const auto h11_m2 = ph11_m2[k];

        pc_81[k] = - e_2 * fs_15_8_10 * h5_m4 + e_2 * fs_15_4_30 * h5_m2 - e_3 * fs_45_572_66 * h7_m4 - e_3 * fs_60_13_3 * h7_m2 + e_3 * fs_45_52_10 * r_2 * h5_m4 - e_3 * fs_45_26_30 * r_2 * h5_m2 + e_4 * fs_147_9724_1430 * h9_m4 + e_4 * fs_21_374_385 * h9_m2 + e_4 * fs_45_2431_66 * r_2 * h7_m4 + e_4 * fs_240_221_3 * r_2 * h7_m2 - e_4 * fs_3_26_10 * r_4 * h5_m4 + e_4 * fs_3_13_30 * r_4 * h5_m2 - e_5 * fs_105_8398_286 * h11_m4 + e_5 * fs_3_4199_30030 * h11_m2 - e_5 * fs_7_4862_1430 * r_2 * h9_m4 - e_5 * fs_1_187_385 * r_2 * h9_m2 - e_5 * fs_45_46189_66 * r_4 * h7_m4 - e_5 * fs_240_4199_3 * r_4 * h7_m2 + e_5 * fs_1_221_10 * r_6 * h5_m4 - e_5 * fs_2_221_30 * r_6 * h5_m2;
    }

    // NOTE: the angular components are formed in 124 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph1_m1, ph3_m3, ph3_m1, ph5_m3, ph5_m1, ph7_m3, ph7_m1, ph9_m3, ph9_m1, ph11_m3, ph11_m1, ab_2, pc_82 : simd::cache_line_size())
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

        const auto h1_m1 = ph1_m1[k];
        const auto h3_m3 = ph3_m3[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h5_m3 = ph5_m3[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h7_m3 = ph7_m3[k];
        const auto h7_m1 = ph7_m1[k];
        const auto h9_m3 = ph9_m3[k];
        const auto h9_m1 = ph9_m1[k];
        const auto h11_m3 = ph11_m3[k];
        const auto h11_m1 = ph11_m1[k];

        pc_82[k] = e_0 * fs_945_32_14 * h1_m1 - e_1 * fs_315_16_35 * h3_m3 - e_1 * fs_315_16_21 * h3_m1 - e_1 * fs_945_16_14 * r_2 * h1_m1 + e_2 * fs_55_4_5 * h5_m3 + e_2 * fs_5_2_210 * h5_m1 + e_2 * fs_35_2_35 * r_2 * h3_m3 + e_2 * fs_35_2_21 * r_2 * h3_m1 + e_2 * fs_135_4_14 * r_4 * h1_m1 - e_3 * fs_3105_1144_6 * h7_m3 - e_3 * fs_5505_1144_2 * h7_m1 - e_3 * fs_165_26_5 * r_2 * h5_m3 - e_3 * fs_15_13_210 * r_2 * h5_m1 - e_3 * fs_105_22_35 * r_4 * h3_m3 - e_3 * fs_105_22_21 * r_4 * h3_m1 - e_3 * fs_15_2_14 * r_6 * h1_m1 + e_4 * fs_441_4862_165 * h9_m3 + e_4 * fs_126_2431_70 * h9_m1 + e_4 * fs_3105_4862_6 * r_2 * h7_m3 + e_4 * fs_5505_4862_2 * r_2 * h7_m1 + e_4 * fs_11_13_5 * r_4 * h5_m3 + e_4 * fs_2_13_210 * r_4 * h5_m1 + e_4 * fs_70_143_35 * r_6 * h3_m3 + e_4 * fs_70_143_21 * r_6 * h3_m1 + e_4 * fs_15_22_14 * r_8 * h1_m1 - e_5 * fs_21_4199_2145 * h11_m3 + e_5 * fs_45_4199_231 * h11_m1 - e_5 * fs_21_2431_165 * r_2 * h9_m3 - e_5 * fs_12_2431_70 * r_2 * h9_m1 - e_5 * fs_3105_92378_6 * r_4 * h7_m3 - e_5 * fs_5505_92378_2 * r_4 * h7_m1 - e_5 * fs_22_663_5 * r_6 * h5_m3 - e_5 * fs_4_663_210 * r_6 * h5_m1 - e_5 * fs_7_429_35 * r_8 * h3_m3 - e_5 * fs_7_429_21 * r_8 * h3_m1 - e_5 * fs_3_143_14 * r_10 * h1_m1;
    }

    // NOTE: the angular components are formed in 124 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph1_p1, ph3_m2, ph5_m2, ph5_p1, ph7_m2, ph7_p1, ph9_m2, ph9_p1, ph11_m2, ph11_p1, ab_2, pc_83, pc_84 : simd::cache_line_size())
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

        const auto h1_p1 = ph1_p1[k];
        const auto h3_m2 = ph3_m2[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h9_m2 = ph9_m2[k];
        const auto h9_p1 = ph9_p1[k];
        const auto h11_m2 = ph11_m2[k];
        const auto h11_p1 = ph11_p1[k];

        pc_83[k] = - e_1 * fs_315_16_21 * h3_m2 + e_2 * fs_20_1_3 * h5_m2 + e_2 * fs_35_2_21 * r_2 * h3_m2 - e_3 * fs_405_286_30 * h7_m2 - e_3 * fs_120_13_3 * r_2 * h5_m2 - e_3 * fs_105_22_21 * r_4 * h3_m2 + e_4 * fs_252_2431_154 * h9_m2 + e_4 * fs_810_2431_30 * r_2 * h7_m2 + e_4 * fs_16_13_3 * r_4 * h5_m2 + e_4 * fs_70_143_21 * r_6 * h3_m2 - e_5 * fs_18_4199_3003 * h11_m2 - e_5 * fs_24_2431_154 * r_2 * h9_m2 - e_5 * fs_810_46189_30 * r_4 * h7_m2 - e_5 * fs_32_663_3 * r_6 * h5_m2 - e_5 * fs_7_429_21 * r_8 * h3_m2;

        pc_84[k] = - e_0 * fs_945_32_15 * h1_p1 + e_1 * fs_945_16_15 * r_2 * h1_p1 + e_2 * f_15_1 * h5_p1 - e_2 * fs_135_4_15 * r_4 * h1_p1 - e_3 * fs_15_26_105 * h7_p1 - e_3 * f_90_13 * r_2 * h5_p1 + e_3 * fs_15_2_15 * r_6 * h1_p1 + e_4 * fs_147_187_3 * h9_p1 + e_4 * fs_30_221_105 * r_2 * h7_p1 + e_4 * f_12_13 * r_4 * h5_p1 - e_4 * fs_15_22_15 * r_8 * h1_p1 - e_5 * fs_126_4199_110 * h11_p1 - e_5 * fs_14_187_3 * r_2 * h9_p1 - e_5 * fs_30_4199_105 * r_4 * h7_p1 - e_5 * f_8_221 * r_6 * h5_p1 + e_5 * fs_3_143_15 * r_10 * h1_p1;
    }

    // NOTE: the angular components are formed in 124 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph1_0, ph3_0, ph3_p2, ph5_0, ph5_p2, ph7_0, ph7_p2, ph9_0, ph9_p2, ph11_0, ph11_p2, ab_2, pc_85 : simd::cache_line_size())
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

        const auto h1_0 = ph1_0[k];
        const auto h3_0 = ph3_0[k];
        const auto h3_p2 = ph3_p2[k];
        const auto h5_0 = ph5_0[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h7_0 = ph7_0[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h9_0 = ph9_0[k];
        const auto h9_p2 = ph9_p2[k];
        const auto h11_0 = ph11_0[k];
        const auto h11_p2 = ph11_p2[k];

        pc_85[k] = e_0 * fs_945_32_35 * h1_0 - e_1 * fs_315_16_35 * h3_0 - e_1 * fs_315_16_21 * h3_p2 - e_1 * fs_945_16_35 * r_2 * h1_0 + e_2 * fs_5_1_35 * h5_0 + e_2 * fs_20_1_3 * h5_p2 + e_2 * fs_35_2_35 * r_2 * h3_0 + e_2 * fs_35_2_21 * r_2 * h3_p2 + e_2 * fs_135_4_35 * r_4 * h1_0 - e_3 * fs_60_143_35 * h7_0 - e_3 * fs_405_286_30 * h7_p2 - e_3 * fs_30_13_35 * r_2 * h5_0 - e_3 * fs_120_13_3 * r_2 * h5_p2 - e_3 * fs_105_22_35 * r_4 * h3_0 - e_3 * fs_105_22_21 * r_4 * h3_p2 - e_3 * fs_15_2_35 * r_6 * h1_0 - e_4 * fs_189_2431_35 * h9_0 + e_4 * fs_252_2431_154 * h9_p2 + e_4 * fs_240_2431_35 * r_2 * h7_0 + e_4 * fs_810_2431_30 * r_2 * h7_p2 + e_4 * fs_4_13_35 * r_4 * h5_0 + e_4 * fs_16_13_3 * r_4 * h5_p2 + e_4 * fs_70_143_35 * r_6 * h3_0 + e_4 * fs_70_143_21 * r_6 * h3_p2 + e_4 * fs_15_22_35 * r_8 * h1_0 + e_5 * fs_198_4199_35 * h11_0 - e_5 * fs_18_4199_3003 * h11_p2 + e_5 * fs_18_2431_35 * r_2 * h9_0 - e_5 * fs_24_2431_154 * r_2 * h9_p2 - e_5 * fs_240_46189_35 * r_4 * h7_0 - e_5 * fs_810_46189_30 * r_4 * h7_p2 - e_5 * fs_8_663_35 * r_6 * h5_0 - e_5 * fs_32_663_3 * r_6 * h5_p2 - e_5 * fs_7_429_35 * r_8 * h3_0 - e_5 * fs_7_429_21 * r_8 * h3_p2 - e_5 * fs_3_143_35 * r_10 * h1_0;
    }

    // NOTE: the angular components are formed in 124 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph1_p1, ph3_p1, ph3_p3, ph5_p1, ph5_p3, ph7_p1, ph7_p3, ph9_p1, ph9_p3, ph11_p1, ph11_p3, ab_2, pc_86 : simd::cache_line_size())
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

        const auto h1_p1 = ph1_p1[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h3_p3 = ph3_p3[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h7_p3 = ph7_p3[k];
        const auto h9_p1 = ph9_p1[k];
        const auto h9_p3 = ph9_p3[k];
        const auto h11_p1 = ph11_p1[k];
        const auto h11_p3 = ph11_p3[k];

        pc_86[k] = e_0 * fs_945_32_14 * h1_p1 - e_1 * fs_315_16_21 * h3_p1 - e_1 * fs_315_16_35 * h3_p3 - e_1 * fs_945_16_14 * r_2 * h1_p1 + e_2 * fs_5_2_210 * h5_p1 + e_2 * fs_55_4_5 * h5_p3 + e_2 * fs_35_2_21 * r_2 * h3_p1 + e_2 * fs_35_2_35 * r_2 * h3_p3 + e_2 * fs_135_4_14 * r_4 * h1_p1 - e_3 * fs_5505_1144_2 * h7_p1 - e_3 * fs_3105_1144_6 * h7_p3 - e_3 * fs_15_13_210 * r_2 * h5_p1 - e_3 * fs_165_26_5 * r_2 * h5_p3 - e_3 * fs_105_22_21 * r_4 * h3_p1 - e_3 * fs_105_22_35 * r_4 * h3_p3 - e_3 * fs_15_2_14 * r_6 * h1_p1 + e_4 * fs_126_2431_70 * h9_p1 + e_4 * fs_441_4862_165 * h9_p3 + e_4 * fs_5505_4862_2 * r_2 * h7_p1 + e_4 * fs_3105_4862_6 * r_2 * h7_p3 + e_4 * fs_2_13_210 * r_4 * h5_p1 + e_4 * fs_11_13_5 * r_4 * h5_p3 + e_4 * fs_70_143_21 * r_6 * h3_p1 + e_4 * fs_70_143_35 * r_6 * h3_p3 + e_4 * fs_15_22_14 * r_8 * h1_p1 + e_5 * fs_45_4199_231 * h11_p1 - e_5 * fs_21_4199_2145 * h11_p3 - e_5 * fs_12_2431_70 * r_2 * h9_p1 - e_5 * fs_21_2431_165 * r_2 * h9_p3 - e_5 * fs_5505_92378_2 * r_4 * h7_p1 - e_5 * fs_3105_92378_6 * r_4 * h7_p3 - e_5 * fs_4_663_210 * r_6 * h5_p1 - e_5 * fs_22_663_5 * r_6 * h5_p3 - e_5 * fs_7_429_21 * r_8 * h3_p1 - e_5 * fs_7_429_35 * r_8 * h3_p3 - e_5 * fs_3_143_14 * r_10 * h1_p1;
    }

    // NOTE: the angular components are formed in 124 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_2, pe_3, pe_4, pe_5, ph5_p2, ph5_p4, ph7_p2, ph7_p4, ph9_p2, ph9_p4, ph11_p2, ph11_p4, ab_2, pc_87 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h5_p2 = ph5_p2[k];
        const auto h5_p4 = ph5_p4[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h7_p4 = ph7_p4[k];
        const auto h9_p2 = ph9_p2[k];
        const auto h9_p4 = ph9_p4[k];
        const auto h11_p2 = ph11_p2[k];
        const auto h11_p4 = ph11_p4[k];

        pc_87[k] = e_2 * fs_15_4_30 * h5_p2 - e_2 * fs_15_8_10 * h5_p4 - e_3 * fs_60_13_3 * h7_p2 - e_3 * fs_45_572_66 * h7_p4 - e_3 * fs_45_26_30 * r_2 * h5_p2 + e_3 * fs_45_52_10 * r_2 * h5_p4 + e_4 * fs_21_374_385 * h9_p2 + e_4 * fs_147_9724_1430 * h9_p4 + e_4 * fs_240_221_3 * r_2 * h7_p2 + e_4 * fs_45_2431_66 * r_2 * h7_p4 + e_4 * fs_3_13_30 * r_4 * h5_p2 - e_4 * fs_3_26_10 * r_4 * h5_p4 + e_5 * fs_3_4199_30030 * h11_p2 - e_5 * fs_105_8398_286 * h11_p4 - e_5 * fs_1_187_385 * r_2 * h9_p2 - e_5 * fs_7_4862_1430 * r_2 * h9_p4 - e_5 * fs_240_4199_3 * r_4 * h7_p2 - e_5 * fs_45_46189_66 * r_4 * h7_p4 - e_5 * fs_2_221_30 * r_6 * h5_p2 + e_5 * fs_1_221_10 * r_6 * h5_p4;
    }

    // NOTE: the angular components are formed in 124 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph3_p3, ph5_p3, ph5_p5, ph7_p3, ph7_p5, ph9_p3, ph9_p5, ph11_p3, ph11_p5, ab_2, pc_88 : simd::cache_line_size())
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

        const auto h3_p3 = ph3_p3[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h5_p5 = ph5_p5[k];
        const auto h7_p3 = ph7_p3[k];
        const auto h7_p5 = ph7_p5[k];
        const auto h9_p3 = ph9_p3[k];
        const auto h9_p5 = ph9_p5[k];
        const auto h11_p3 = ph11_p3[k];
        const auto h11_p5 = ph11_p5[k];

        pc_88[k] = e_1 * fs_315_16_42 * h3_p3 - e_2 * fs_65_8_6 * h5_p3 - e_2 * fs_75_8_30 * h5_p5 - e_2 * fs_35_2_42 * r_2 * h3_p3 - e_3 * fs_75_44_5 * h7_p3 + e_3 * fs_555_572_55 * h7_p5 + e_3 * fs_15_4_6 * r_2 * h5_p3 + e_3 * fs_225_52_30 * r_2 * h5_p5 + e_3 * fs_105_22_42 * r_4 * h3_p3 + e_4 * fs_3087_9724_22 * h9_p3 - e_4 * fs_21_9724_30030 * h9_p5 + e_4 * fs_75_187_5 * r_2 * h7_p3 - e_4 * fs_555_2431_55 * r_2 * h7_p5 - e_4 * fs_1_2_6 * r_4 * h5_p3 - e_4 * fs_15_26_30 * r_4 * h5_p5 - e_4 * fs_70_143_42 * r_6 * h3_p3 + e_5 * fs_21_4199_286 * h11_p3 - e_5 * fs_6_4199_15015 * h11_p5 - e_5 * fs_147_4862_22 * r_2 * h9_p3 + e_5 * fs_1_4862_30030 * r_2 * h9_p5 - e_5 * fs_75_3553_5 * r_4 * h7_p3 + e_5 * fs_555_46189_55 * r_4 * h7_p5 + e_5 * fs_1_51_6 * r_6 * h5_p3 + e_5 * fs_5_221_30 * r_6 * h5_p5 + e_5 * fs_7_429_42 * r_8 * h3_p3;
    }

    // NOTE: the angular components are formed in 124 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_2, pe_3, pe_4, pe_5, ph5_p4, ph7_p4, ph7_p6, ph9_p4, ph9_p6, ph11_p4, ph11_p6, ab_2, pc_89 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h5_p4 = ph5_p4[k];
        const auto h7_p4 = ph7_p4[k];
        const auto h7_p6 = ph7_p6[k];
        const auto h9_p4 = ph9_p4[k];
        const auto h9_p6 = ph9_p6[k];
        const auto h11_p4 = ph11_p4[k];
        const auto h11_p6 = ph11_p6[k];

        pc_89[k] = - e_2 * fs_45_8_66 * h5_p4 + e_3 * fs_75_52_10 * h7_p4 + e_3 * fs_15_13_65 * h7_p6 + e_3 * fs_135_52_66 * r_2 * h5_p4 + e_4 * fs_147_884_78 * h9_p4 - e_4 * fs_63_442_91 * h9_p6 - e_4 * fs_75_221_10 * r_2 * h7_p4 - e_4 * fs_60_221_65 * r_2 * h7_p6 - e_4 * fs_9_26_66 * r_4 * h5_p4 + e_5 * fs_21_8398_390 * h11_p4 - e_5 * fs_6_4199_7735 * h11_p6 - e_5 * fs_7_442_78 * r_2 * h9_p4 + e_5 * fs_3_221_91 * r_2 * h9_p6 + e_5 * fs_75_4199_10 * r_4 * h7_p4 + e_5 * fs_60_4199_65 * r_4 * h7_p6 + e_5 * fs_3_221_66 * r_6 * h5_p4;
    }

    // NOTE: the angular components are formed in 124 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_2, pe_3, pe_4, pe_5, ph5_m4, ph5_p5, ph7_m4, ph7_p5, ph7_p7, ph9_m8, ph9_m4, ph9_p5, ph9_p7, ph11_m8, ph11_m4, ph11_p5, ph11_p7, ab_2, pc_90, pc_91 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h5_m4 = ph5_m4[k];
        const auto h5_p5 = ph5_p5[k];
        const auto h7_m4 = ph7_m4[k];
        const auto h7_p5 = ph7_p5[k];
        const auto h7_p7 = ph7_p7[k];
        const auto h9_m8 = ph9_m8[k];
        const auto h9_m4 = ph9_m4[k];
        const auto h9_p5 = ph9_p5[k];
        const auto h9_p7 = ph9_p7[k];
        const auto h11_m8 = ph11_m8[k];
        const auto h11_m4 = ph11_m4[k];
        const auto h11_p5 = ph11_p5[k];
        const auto h11_p7 = ph11_p7[k];

        pc_90[k] = e_2 * fs_15_4_55 * h5_p5 + e_3 * fs_225_104_30 * h7_p5 - e_3 * fs_15_104_2730 * h7_p7 - e_3 * fs_45_26_55 * r_2 * h5_p5 + e_4 * fs_21_442_455 * h9_p5 - e_4 * fs_42_221_91 * h9_p7 - e_4 * fs_225_442_30 * r_2 * h7_p5 + e_4 * fs_15_442_2730 * r_2 * h7_p7 + e_4 * fs_3_13_55 * r_4 * h5_p5 + e_5 * fs_3_4199_910 * h11_p5 - e_5 * fs_3_4199_9282 * h11_p7 - e_5 * fs_1_221_455 * r_2 * h9_p5 + e_5 * fs_4_221_91 * r_2 * h9_p7 + e_5 * fs_225_8398_30 * r_4 * h7_p5 - e_5 * fs_15_8398_2730 * r_4 * h7_p7 - e_5 * fs_2_221_55 * r_6 * h5_p5;

        pc_91[k] = - e_2 * fs_15_4_154 * h5_m4 - e_3 * fs_75_104_210 * h7_m4 + e_3 * fs_45_26_154 * r_2 * h5_m4 - e_4 * fs_21_221_442 * h9_m8 - e_4 * fs_21_442_182 * h9_m4 + e_4 * fs_75_442_210 * r_2 * h7_m4 - e_4 * fs_3_13_154 * r_4 * h5_m4 - e_5 * fs_3_4199_25194 * h11_m8 - e_5 * fs_3_8398_910 * h11_m4 + e_5 * fs_2_221_442 * r_2 * h9_m8 + e_5 * fs_1_221_182 * r_2 * h9_m4 - e_5 * fs_75_8398_210 * r_4 * h7_m4 + e_5 * fs_2_221_154 * r_6 * h5_m4;
    }

    // NOTE: the angular components are formed in 124 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph3_m3, ph5_m3, ph7_m7, ph7_m3, ph9_m7, ph9_m3, ph11_m7, ph11_m3, ab_2, pc_92 : simd::cache_line_size())
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

        const auto h3_m3 = ph3_m3[k];
        const auto h5_m3 = ph5_m3[k];
        const auto h7_m7 = ph7_m7[k];
        const auto h7_m3 = ph7_m3[k];
        const auto h9_m7 = ph9_m7[k];
        const auto h9_m3 = ph9_m3[k];
        const auto h11_m7 = ph11_m7[k];
        const auto h11_m3 = ph11_m3[k];

        pc_92[k] = e_1 * fs_315_16_33 * h3_m3 + e_2 * fs_5_4_231 * h5_m3 - e_2 * fs_35_2_33 * r_2 * h3_m3 + e_3 * fs_105_104_130 * h7_m7 - e_3 * fs_375_1144_770 * h7_m3 - e_3 * fs_15_26_231 * r_2 * h5_m3 + e_3 * fs_105_22_33 * r_4 * h3_m3 - e_4 * fs_21_221_39 * h9_m7 - e_4 * fs_189_442_7 * h9_m3 - e_4 * fs_105_442_130 * r_2 * h7_m7 + e_4 * fs_375_4862_770 * r_2 * h7_m3 + e_4 * fs_1_13_231 * r_4 * h5_m3 - e_4 * fs_70_143_33 * r_6 * h3_m3 - e_5 * fs_36_4199_442 * h11_m7 - e_5 * fs_12_4199_91 * h11_m3 + e_5 * fs_2_221_39 * r_2 * h9_m7 + e_5 * fs_9_221_7 * r_2 * h9_m3 + e_5 * fs_105_8398_130 * r_4 * h7_m7 - e_5 * fs_375_92378_770 * r_4 * h7_m3 - e_5 * fs_2_663_231 * r_6 * h5_m3 + e_5 * fs_7_429_33 * r_8 * h3_m3;
    }

    // NOTE: the angular components are formed in 124 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph3_m2, ph5_m2, ph7_m6, ph7_m2, ph9_m6, ph9_m2, ph11_m6, ph11_m2, ab_2, pc_93 : simd::cache_line_size())
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

        const auto h3_m2 = ph3_m2[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h7_m6 = ph7_m6[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h9_m6 = ph9_m6[k];
        const auto h9_m2 = ph9_m2[k];
        const auto h11_m6 = ph11_m6[k];
        const auto h11_m2 = ph11_m2[k];

        pc_93[k] = - e_1 * f_945_16 * h3_m2 + e_2 * fs_15_1_7 * h5_m2 + e_2 * f_105_2 * r_2 * h3_m2 - e_3 * fs_15_1144_10010 * h7_m6 - e_3 * fs_435_1144_70 * h7_m2 - e_3 * fs_90_13_7 * r_2 * h5_m2 - e_3 * f_315_22 * r_4 * h3_m2 + e_4 * fs_189_4862_286 * h9_m6 - e_4 * fs_861_4862_66 * h9_m2 + e_4 * fs_15_4862_10010 * r_2 * h7_m6 + e_4 * fs_435_4862_70 * r_2 * h7_m2 + e_4 * fs_12_13_7 * r_4 * h5_m2 + e_4 * f_210_143 * r_6 * h3_m2 - e_5 * fs_6_4199_24310 * h11_m6 - e_5 * fs_18_4199_143 * h11_m2 - e_5 * fs_9_2431_286 * r_2 * h9_m6 + e_5 * fs_41_2431_66 * r_2 * h9_m2 - e_5 * fs_15_92378_10010 * r_4 * h7_m6 - e_5 * fs_435_92378_70 * r_4 * h7_m2 - e_5 * fs_8_221_7 * r_6 * h5_m2 - e_5 * f_7_143 * r_8 * h3_m2;
    }

    // NOTE: the angular components are formed in 124 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph1_m1, ph3_m1, ph5_m5, ph5_m1, ph7_m5, ph7_m1, ph9_m5, ph9_m1, ph11_m5, ph11_m1, ab_2, pc_94 : simd::cache_line_size())
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

        const auto h1_m1 = ph1_m1[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h5_m5 = ph5_m5[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h7_m5 = ph7_m5[k];
        const auto h7_m1 = ph7_m1[k];
        const auto h9_m5 = ph9_m5[k];
        const auto h9_m1 = ph9_m1[k];
        const auto h11_m5 = ph11_m5[k];
        const auto h11_m1 = ph11_m1[k];

        pc_94[k] = e_0 * fs_2835_32_2 * h1_m1 - e_1 * fs_945_16_3 * h3_m1 - e_1 * fs_2835_16_2 * r_2 * h1_m1 + e_2 * fs_75_4_7 * h5_m5 + e_2 * fs_15_4_30 * h5_m1 + e_2 * fs_105_2_3 * r_2 * h3_m1 + e_2 * fs_405_4_2 * r_4 * h1_m1 - e_3 * fs_435_1144_462 * h7_m5 + e_3 * fs_1035_1144_14 * h7_m1 - e_3 * fs_225_26_7 * r_2 * h5_m5 - e_3 * fs_45_26_30 * r_2 * h5_m1 - e_3 * fs_315_22_3 * r_4 * h3_m1 - e_3 * fs_45_2_2 * r_6 * h1_m1 + e_4 * fs_525_4862_143 * h9_m5 - e_4 * fs_2247_4862_10 * h9_m1 + e_4 * fs_435_4862_462 * r_2 * h7_m5 - e_4 * fs_1035_4862_14 * r_2 * h7_m1 + e_4 * fs_15_13_7 * r_4 * h5_m5 + e_4 * fs_3_13_30 * r_4 * h5_m1 + e_4 * fs_210_143_3 * r_6 * h3_m1 + e_4 * fs_45_22_2 * r_8 * h1_m1 - e_5 * fs_60_4199_286 * h11_m5 - e_5 * fs_60_4199_33 * h11_m1 - e_5 * fs_25_2431_143 * r_2 * h9_m5 + e_5 * fs_107_2431_10 * r_2 * h9_m1 - e_5 * fs_435_92378_462 * r_4 * h7_m5 + e_5 * fs_1035_92378_14 * r_4 * h7_m1 - e_5 * fs_10_221_7 * r_6 * h5_m5 - e_5 * fs_2_221_30 * r_6 * h5_m1 - e_5 * fs_7_143_3 * r_8 * h3_m1 - e_5 * fs_9_143_2 * r_10 * h1_m1;
    }

    // NOTE: the angular components are formed in 124 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_2, pe_3, pe_4, pe_5, ph5_m4, ph7_m4, ph9_m4, ph11_m4, ab_2, pc_95 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h5_m4 = ph5_m4[k];
        const auto h7_m4 = ph7_m4[k];
        const auto h9_m4 = ph9_m4[k];
        const auto h11_m4 = ph11_m4[k];

        pc_95[k] = e_2 * fs_15_4_70 * h5_m4 - e_3 * fs_405_1144_462 * h7_m4 - e_3 * fs_45_26_70 * r_2 * h5_m4 + e_4 * fs_63_4862_10010 * h9_m4 + e_4 * fs_405_4862_462 * r_2 * h7_m4 + e_4 * fs_3_13_70 * r_4 * h5_m4 - e_5 * fs_45_8398_2002 * h11_m4 - e_5 * fs_3_2431_10010 * r_2 * h9_m4 - e_5 * fs_405_92378_462 * r_4 * h7_m4 - e_5 * fs_2_221_70 * r_6 * h5_m4;
    }

    // NOTE: the angular components are formed in 124 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph1_m1, ph3_m3, ph3_m1, ph5_m3, ph5_m1, ph7_m3, ph7_m1, ph9_m3, ph9_m1, ph11_m3, ph11_m1, ab_2, pc_96 : simd::cache_line_size())
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

        const auto h1_m1 = ph1_m1[k];
        const auto h3_m3 = ph3_m3[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h5_m3 = ph5_m3[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h7_m3 = ph7_m3[k];
        const auto h7_m1 = ph7_m1[k];
        const auto h9_m3 = ph9_m3[k];
        const auto h9_m1 = ph9_m1[k];
        const auto h11_m3 = ph11_m3[k];
        const auto h11_m1 = ph11_m1[k];

        pc_96[k] = e_0 * fs_945_32_5 * h1_m1 + e_1 * fs_2205_32_2 * h3_m3 + e_1 * fs_315_32_30 * h3_m1 - e_1 * fs_945_16_5 * r_2 * h1_m1 - e_2 * fs_5_4_14 * h5_m3 - e_2 * fs_20_1_3 * h5_m1 - e_2 * fs_245_4_2 * r_2 * h3_m3 - e_2 * fs_35_4_30 * r_2 * h3_m1 + e_2 * fs_135_4_5 * r_4 * h1_m1 - e_3 * fs_135_572_105 * h7_m3 + e_3 * fs_735_572_35 * h7_m1 + e_3 * fs_15_26_14 * r_2 * h5_m3 + e_3 * fs_120_13_3 * r_2 * h5_m1 + e_3 * fs_735_44_2 * r_4 * h3_m3 + e_3 * fs_105_44_30 * r_4 * h3_m1 - e_3 * fs_15_2_5 * r_6 * h1_m1 + e_4 * fs_189_4862_462 * h9_m3 - e_4 * f_1449_2431 * h9_m1 + e_4 * fs_135_2431_105 * r_2 * h7_m3 - e_4 * fs_735_2431_35 * r_2 * h7_m1 - e_4 * fs_1_13_14 * r_4 * h5_m3 - e_4 * fs_16_13_3 * r_4 * h5_m1 - e_4 * fs_245_143_2 * r_6 * h3_m3 - e_4 * fs_35_143_30 * r_6 * h3_m1 + e_4 * fs_15_22_5 * r_8 * h1_m1 - e_5 * fs_12_4199_6006 * h11_m3 - e_5 * fs_36_4199_330 * h11_m1 - e_5 * fs_9_2431_462 * r_2 * h9_m3 + e_5 * f_138_2431 * r_2 * h9_m1 - e_5 * fs_135_46189_105 * r_4 * h7_m3 + e_5 * fs_735_46189_35 * r_4 * h7_m1 + e_5 * fs_2_663_14 * r_6 * h5_m3 + e_5 * fs_32_663_3 * r_6 * h5_m1 + e_5 * fs_49_858_2 * r_8 * h3_m3 + e_5 * fs_7_858_30 * r_8 * h3_m1 - e_5 * fs_3_143_5 * r_10 * h1_m1;
    }

    // NOTE: the angular components are formed in 124 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph3_p2, ph5_p2, ph7_p2, ph9_p2, ph11_p2, ab_2, pc_97 : simd::cache_line_size())
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

        const auto h3_p2 = ph3_p2[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h9_p2 = ph9_p2[k];
        const auto h11_p2 = ph11_p2[k];

        pc_97[k] = e_1 * fs_945_16_7 * h3_p2 - e_2 * f_45_1 * h5_p2 - e_2 * fs_105_2_7 * r_2 * h3_p2 + e_3 * fs_945_572_10 * h7_p2 + e_3 * f_270_13 * r_2 * h5_p2 + e_3 * fs_315_22_7 * r_4 * h3_p2 + e_4 * fs_21_2431_462 * h9_p2 - e_4 * fs_945_2431_10 * r_2 * h7_p2 - e_4 * f_36_13 * r_4 * h5_p2 - e_4 * fs_210_143_7 * r_6 * h3_p2 - e_5 * fs_36_4199_1001 * h11_p2 - e_5 * fs_2_2431_462 * r_2 * h9_p2 + e_5 * fs_945_46189_10 * r_4 * h7_p2 + e_5 * f_24_221 * r_6 * h5_p2 + e_5 * fs_7_143_7 * r_8 * h3_p2;
    }

    // NOTE: the angular components are formed in 124 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph1_p1, ph3_p1, ph3_p3, ph5_p1, ph5_p3, ph7_p1, ph7_p3, ph9_p1, ph9_p3, ph11_p1, ph11_p3, ab_2, pc_98 : simd::cache_line_size())
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

        const auto h1_p1 = ph1_p1[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h3_p3 = ph3_p3[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h7_p3 = ph7_p3[k];
        const auto h9_p1 = ph9_p1[k];
        const auto h9_p3 = ph9_p3[k];
        const auto h11_p1 = ph11_p1[k];
        const auto h11_p3 = ph11_p3[k];

        pc_98[k] = - e_0 * fs_945_32_5 * h1_p1 - e_1 * fs_315_32_30 * h3_p1 + e_1 * fs_2205_32_2 * h3_p3 + e_1 * fs_945_16_5 * r_2 * h1_p1 + e_2 * fs_20_1_3 * h5_p1 - e_2 * fs_5_4_14 * h5_p3 + e_2 * fs_35_4_30 * r_2 * h3_p1 - e_2 * fs_245_4_2 * r_2 * h3_p3 - e_2 * fs_135_4_5 * r_4 * h1_p1 - e_3 * fs_735_572_35 * h7_p1 - e_3 * fs_135_572_105 * h7_p3 - e_3 * fs_120_13_3 * r_2 * h5_p1 + e_3 * fs_15_26_14 * r_2 * h5_p3 - e_3 * fs_105_44_30 * r_4 * h3_p1 + e_3 * fs_735_44_2 * r_4 * h3_p3 + e_3 * fs_15_2_5 * r_6 * h1_p1 + e_4 * f_1449_2431 * h9_p1 + e_4 * fs_189_4862_462 * h9_p3 + e_4 * fs_735_2431_35 * r_2 * h7_p1 + e_4 * fs_135_2431_105 * r_2 * h7_p3 + e_4 * fs_16_13_3 * r_4 * h5_p1 - e_4 * fs_1_13_14 * r_4 * h5_p3 + e_4 * fs_35_143_30 * r_6 * h3_p1 - e_4 * fs_245_143_2 * r_6 * h3_p3 - e_4 * fs_15_22_5 * r_8 * h1_p1 + e_5 * fs_36_4199_330 * h11_p1 - e_5 * fs_12_4199_6006 * h11_p3 - e_5 * f_138_2431 * r_2 * h9_p1 - e_5 * fs_9_2431_462 * r_2 * h9_p3 - e_5 * fs_735_46189_35 * r_4 * h7_p1 - e_5 * fs_135_46189_105 * r_4 * h7_p3 - e_5 * fs_32_663_3 * r_6 * h5_p1 + e_5 * fs_2_663_14 * r_6 * h5_p3 - e_5 * fs_7_858_30 * r_8 * h3_p1 + e_5 * fs_49_858_2 * r_8 * h3_p3 + e_5 * fs_3_143_5 * r_10 * h1_p1;
    }

    // NOTE: the angular components are formed in 124 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph1_0, ph3_0, ph5_0, ph5_p4, ph7_0, ph7_p4, ph9_0, ph9_p4, ph11_0, ph11_p4, ab_2, pc_99 : simd::cache_line_size())
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

        const auto h1_0 = ph1_0[k];
        const auto h3_0 = ph3_0[k];
        const auto h5_0 = ph5_0[k];
        const auto h5_p4 = ph5_p4[k];
        const auto h7_0 = ph7_0[k];
        const auto h7_p4 = ph7_p4[k];
        const auto h9_0 = ph9_0[k];
        const auto h9_p4 = ph9_p4[k];
        const auto h11_0 = ph11_0[k];
        const auto h11_p4 = ph11_p4[k];

        pc_99[k] = e_0 * fs_945_8_2 * h1_0 - e_1 * fs_315_8_2 * h3_0 - e_1 * fs_945_4_2 * r_2 * h1_0 - e_2 * fs_25_2_2 * h5_0 + e_2 * fs_15_4_70 * h5_p4 + e_2 * fs_35_1_2 * r_2 * h3_0 + e_2 * fs_135_1_2 * r_4 * h1_0 + e_3 * fs_4305_572_2 * h7_0 - e_3 * fs_405_1144_462 * h7_p4 + e_3 * fs_75_13_2 * r_2 * h5_0 - e_3 * fs_45_26_70 * r_2 * h5_p4 - e_3 * fs_105_11_2 * r_4 * h3_0 - e_3 * fs_30_1_2 * r_6 * h1_0 - e_4 * fs_2835_2431_2 * h9_0 + e_4 * fs_63_4862_10010 * h9_p4 - e_4 * fs_4305_2431_2 * r_2 * h7_0 + e_4 * fs_405_4862_462 * r_2 * h7_p4 - e_4 * fs_10_13_2 * r_4 * h5_0 + e_4 * fs_3_13_70 * r_4 * h5_p4 + e_4 * fs_140_143_2 * r_6 * h3_0 + e_4 * fs_30_11_2 * r_8 * h1_0 - e_5 * fs_495_4199_2 * h11_0 - e_5 * fs_45_8398_2002 * h11_p4 + e_5 * fs_270_2431_2 * r_2 * h9_0 - e_5 * fs_3_2431_10010 * r_2 * h9_p4 + e_5 * fs_4305_46189_2 * r_4 * h7_0 - e_5 * fs_405_92378_462 * r_4 * h7_p4 + e_5 * fs_20_663_2 * r_6 * h5_0 - e_5 * fs_2_221_70 * r_6 * h5_p4 - e_5 * fs_14_429_2 * r_8 * h3_0 - e_5 * fs_12_143_2 * r_10 * h1_0;
    }

    // NOTE: the angular components are formed in 124 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph1_p1, ph3_p1, ph5_p1, ph5_p5, ph7_p1, ph7_p5, ph9_p1, ph9_p5, ph11_p1, ph11_p5, ab_2, pc_100 : simd::cache_line_size())
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

        const auto h1_p1 = ph1_p1[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h5_p5 = ph5_p5[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h7_p5 = ph7_p5[k];
        const auto h9_p1 = ph9_p1[k];
        const auto h9_p5 = ph9_p5[k];
        const auto h11_p1 = ph11_p1[k];
        const auto h11_p5 = ph11_p5[k];

        pc_100[k] = e_0 * fs_2835_32_2 * h1_p1 - e_1 * fs_945_16_3 * h3_p1 - e_1 * fs_2835_16_2 * r_2 * h1_p1 + e_2 * fs_15_4_30 * h5_p1 + e_2 * fs_75_4_7 * h5_p5 + e_2 * fs_105_2_3 * r_2 * h3_p1 + e_2 * fs_405_4_2 * r_4 * h1_p1 + e_3 * fs_1035_1144_14 * h7_p1 - e_3 * fs_435_1144_462 * h7_p5 - e_3 * fs_45_26_30 * r_2 * h5_p1 - e_3 * fs_225_26_7 * r_2 * h5_p5 - e_3 * fs_315_22_3 * r_4 * h3_p1 - e_3 * fs_45_2_2 * r_6 * h1_p1 - e_4 * fs_2247_4862_10 * h9_p1 + e_4 * fs_525_4862_143 * h9_p5 - e_4 * fs_1035_4862_14 * r_2 * h7_p1 + e_4 * fs_435_4862_462 * r_2 * h7_p5 + e_4 * fs_3_13_30 * r_4 * h5_p1 + e_4 * fs_15_13_7 * r_4 * h5_p5 + e_4 * fs_210_143_3 * r_6 * h3_p1 + e_4 * fs_45_22_2 * r_8 * h1_p1 - e_5 * fs_60_4199_33 * h11_p1 - e_5 * fs_60_4199_286 * h11_p5 + e_5 * fs_107_2431_10 * r_2 * h9_p1 - e_5 * fs_25_2431_143 * r_2 * h9_p5 + e_5 * fs_1035_92378_14 * r_4 * h7_p1 - e_5 * fs_435_92378_462 * r_4 * h7_p5 - e_5 * fs_2_221_30 * r_6 * h5_p1 - e_5 * fs_10_221_7 * r_6 * h5_p5 - e_5 * fs_7_143_3 * r_8 * h3_p1 - e_5 * fs_9_143_2 * r_10 * h1_p1;
    }

    // NOTE: the angular components are formed in 124 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph3_p2, ph5_p2, ph7_p2, ph7_p6, ph9_p2, ph9_p6, ph11_p2, ph11_p6, ab_2, pc_101 : simd::cache_line_size())
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

        const auto h3_p2 = ph3_p2[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h7_p6 = ph7_p6[k];
        const auto h9_p2 = ph9_p2[k];
        const auto h9_p6 = ph9_p6[k];
        const auto h11_p2 = ph11_p2[k];
        const auto h11_p6 = ph11_p6[k];

        pc_101[k] = - e_1 * f_945_16 * h3_p2 + e_2 * fs_15_1_7 * h5_p2 + e_2 * f_105_2 * r_2 * h3_p2 - e_3 * fs_435_1144_70 * h7_p2 - e_3 * fs_15_1144_10010 * h7_p6 - e_3 * fs_90_13_7 * r_2 * h5_p2 - e_3 * f_315_22 * r_4 * h3_p2 - e_4 * fs_861_4862_66 * h9_p2 + e_4 * fs_189_4862_286 * h9_p6 + e_4 * fs_435_4862_70 * r_2 * h7_p2 + e_4 * fs_15_4862_10010 * r_2 * h7_p6 + e_4 * fs_12_13_7 * r_4 * h5_p2 + e_4 * f_210_143 * r_6 * h3_p2 - e_5 * fs_18_4199_143 * h11_p2 - e_5 * fs_6_4199_24310 * h11_p6 + e_5 * fs_41_2431_66 * r_2 * h9_p2 - e_5 * fs_9_2431_286 * r_2 * h9_p6 - e_5 * fs_435_92378_70 * r_4 * h7_p2 - e_5 * fs_15_92378_10010 * r_4 * h7_p6 - e_5 * fs_8_221_7 * r_6 * h5_p2 - e_5 * f_7_143 * r_8 * h3_p2;
    }

    // NOTE: the angular components are formed in 124 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph3_p3, ph5_p3, ph7_p3, ph7_p7, ph9_p3, ph9_p7, ph11_p3, ph11_p7, ab_2, pc_102 : simd::cache_line_size())
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

        const auto h3_p3 = ph3_p3[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h7_p3 = ph7_p3[k];
        const auto h7_p7 = ph7_p7[k];
        const auto h9_p3 = ph9_p3[k];
        const auto h9_p7 = ph9_p7[k];
        const auto h11_p3 = ph11_p3[k];
        const auto h11_p7 = ph11_p7[k];

        pc_102[k] = e_1 * fs_315_16_33 * h3_p3 + e_2 * fs_5_4_231 * h5_p3 - e_2 * fs_35_2_33 * r_2 * h3_p3 - e_3 * fs_375_1144_770 * h7_p3 + e_3 * fs_105_104_130 * h7_p7 - e_3 * fs_15_26_231 * r_2 * h5_p3 + e_3 * fs_105_22_33 * r_4 * h3_p3 - e_4 * fs_189_442_7 * h9_p3 - e_4 * fs_21_221_39 * h9_p7 + e_4 * fs_375_4862_770 * r_2 * h7_p3 - e_4 * fs_105_442_130 * r_2 * h7_p7 + e_4 * fs_1_13_231 * r_4 * h5_p3 - e_4 * fs_70_143_33 * r_6 * h3_p3 - e_5 * fs_12_4199_91 * h11_p3 - e_5 * fs_36_4199_442 * h11_p7 + e_5 * fs_9_221_7 * r_2 * h9_p3 + e_5 * fs_2_221_39 * r_2 * h9_p7 - e_5 * fs_375_92378_770 * r_4 * h7_p3 + e_5 * fs_105_8398_130 * r_4 * h7_p7 - e_5 * fs_2_663_231 * r_6 * h5_p3 + e_5 * fs_7_429_33 * r_8 * h3_p3;
    }

    // NOTE: the angular components are formed in 124 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_2, pe_3, pe_4, pe_5, ph5_p4, ph7_p4, ph9_p4, ph9_p8, ph11_p4, ph11_p8, ab_2, pc_103 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h5_p4 = ph5_p4[k];
        const auto h7_p4 = ph7_p4[k];
        const auto h9_p4 = ph9_p4[k];
        const auto h9_p8 = ph9_p8[k];
        const auto h11_p4 = ph11_p4[k];
        const auto h11_p8 = ph11_p8[k];

        pc_103[k] = - e_2 * fs_15_4_154 * h5_p4 - e_3 * fs_75_104_210 * h7_p4 + e_3 * fs_45_26_154 * r_2 * h5_p4 - e_4 * fs_21_442_182 * h9_p4 - e_4 * fs_21_221_442 * h9_p8 + e_4 * fs_75_442_210 * r_2 * h7_p4 - e_4 * fs_3_13_154 * r_4 * h5_p4 - e_5 * fs_3_8398_910 * h11_p4 - e_5 * fs_3_4199_25194 * h11_p8 + e_5 * fs_1_221_182 * r_2 * h9_p4 + e_5 * fs_2_221_442 * r_2 * h9_p8 - e_5 * fs_75_8398_210 * r_4 * h7_p4 + e_5 * fs_2_221_154 * r_6 * h5_p4;
    }

    // NOTE: the angular components are formed in 124 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph3_m3, ph5_m3, ph7_m3, ph9_m9, ph9_m3, ph11_m9, ph11_m3, ab_2, pc_104 : simd::cache_line_size())
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

        const auto h3_m3 = ph3_m3[k];
        const auto h5_m3 = ph5_m3[k];
        const auto h7_m3 = ph7_m3[k];
        const auto h9_m9 = ph9_m9[k];
        const auto h9_m3 = ph9_m3[k];
        const auto h11_m9 = ph11_m9[k];
        const auto h11_m3 = ph11_m3[k];

        pc_104[k] = e_1 * fs_315_32_66 * h3_m3 + e_2 * fs_5_2_462 * h5_m3 - e_2 * fs_35_4_66 * r_2 * h3_m3 + e_3 * fs_225_572_385 * h7_m3 - e_3 * fs_15_13_462 * r_2 * h5_m3 + e_3 * fs_105_44_66 * r_4 * h3_m3 - e_4 * fs_21_442_1326 * h9_m9 + e_4 * fs_21_221_14 * h9_m3 - e_4 * fs_225_2431_385 * r_2 * h7_m3 + e_4 * fs_2_13_462 * r_4 * h5_m3 - e_4 * fs_35_143_66 * r_6 * h3_m3 - e_5 * fs_3_4199_62985 * h11_m9 + e_5 * fs_3_8398_182 * h11_m3 + e_5 * fs_1_221_1326 * r_2 * h9_m9 - e_5 * fs_2_221_14 * r_2 * h9_m3 + e_5 * fs_225_46189_385 * r_4 * h7_m3 - e_5 * fs_4_663_462 * r_6 * h5_m3 + e_5 * fs_7_858_66 * r_8 * h3_m3;
    }

    // NOTE: the angular components are formed in 124 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph3_m2, ph5_m2, ph7_m2, ph9_m8, ph9_m2, ph11_m8, ph11_m2, ab_2, pc_105 : simd::cache_line_size())
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

        const auto h3_m2 = ph3_m2[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h9_m8 = ph9_m8[k];
        const auto h9_m2 = ph9_m2[k];
        const auto h11_m8 = ph11_m8[k];
        const auto h11_m2 = ph11_m2[k];

        pc_105[k] = - e_1 * fs_315_16_33 * h3_m2 + e_2 * fs_5_4_231 * h5_m2 + e_2 * fs_35_2_33 * r_2 * h3_m2 + e_3 * fs_30_143_2310 * h7_m2 - e_3 * fs_15_26_231 * r_2 * h5_m2 - e_3 * fs_105_22_33 * r_4 * h3_m2 + e_4 * fs_21_442_221 * h9_m8 + e_4 * fs_231_442_2 * h9_m2 - e_4 * fs_120_2431_2310 * r_2 * h7_m2 + e_4 * fs_1_13_231 * r_4 * h5_m2 + e_4 * fs_70_143_33 * r_6 * h3_m2 - e_5 * fs_9_4199_12597 * h11_m8 + e_5 * fs_9_4199_39 * h11_m2 - e_5 * fs_1_221_221 * r_2 * h9_m8 - e_5 * fs_11_221_2 * r_2 * h9_m2 + e_5 * fs_120_46189_2310 * r_4 * h7_m2 - e_5 * fs_2_663_231 * r_6 * h5_m2 - e_5 * fs_7_429_33 * r_8 * h3_m2;
    }

    // NOTE: the angular components are formed in 124 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph1_m1, ph3_m1, ph5_m1, ph7_m7, ph7_m1, ph9_m7, ph9_m1, ph11_m7, ph11_m1, ab_2, pc_106 : simd::cache_line_size())
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

        const auto h1_m1 = ph1_m1[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h7_m7 = ph7_m7[k];
        const auto h7_m1 = ph7_m1[k];
        const auto h9_m7 = ph9_m7[k];
        const auto h9_m1 = ph9_m1[k];
        const auto h11_m7 = ph11_m7[k];
        const auto h11_m1 = ph11_m1[k];

        pc_106[k] = e_0 * fs_2835_64_10 * h1_m1 - e_1 * fs_315_16_15 * h3_m1 - e_1 * fs_2835_32_10 * r_2 * h1_m1 - e_2 * fs_65_8_6 * h5_m1 + e_2 * fs_35_2_15 * r_2 * h3_m1 + e_2 * fs_405_8_10 * r_4 * h1_m1 - e_3 * fs_105_572_4290 * h7_m7 + e_3 * fs_45_44_70 * h7_m1 + e_3 * fs_15_4_6 * r_2 * h5_m1 - e_3 * fs_105_22_15 * r_4 * h3_m1 - e_3 * fs_45_4_10 * r_6 * h1_m1 + e_4 * fs_609_4862_143 * h9_m7 + e_4 * fs_7707_9724_2 * h9_m1 + e_4 * fs_105_2431_4290 * r_2 * h7_m7 - e_4 * fs_45_187_70 * r_2 * h7_m1 - e_4 * fs_1_2_6 * r_4 * h5_m1 + e_4 * fs_70_143_15 * r_6 * h3_m1 + e_4 * fs_45_44_10 * r_8 * h1_m1 - e_5 * fs_9_4199_14586 * h11_m7 + e_5 * fs_9_4199_165 * h11_m1 - e_5 * fs_29_2431_143 * r_2 * h9_m7 - e_5 * fs_367_4862_2 * r_2 * h9_m1 - e_5 * fs_105_46189_4290 * r_4 * h7_m7 + e_5 * fs_45_3553_70 * r_4 * h7_m1 + e_5 * fs_1_51_6 * r_6 * h5_m1 - e_5 * fs_7_429_15 * r_8 * h3_m1 - e_5 * fs_9_286_10 * r_10 * h1_m1;
    }

    // NOTE: the angular components are formed in 124 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_3, pe_4, pe_5, ph7_m6, ph9_m6, ph11_m6, ab_2, pc_107 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

        const auto h7_m6 = ph7_m6[k];
        const auto h9_m6 = ph9_m6[k];
        const auto h11_m6 = ph11_m6[k];

        pc_107[k] = - e_3 * fs_45_286_2002 * h7_m6 + e_4 * fs_84_2431_1430 * h9_m6 + e_4 * fs_90_2431_2002 * r_2 * h7_m6 - e_5 * fs_15_4199_4862 * h11_m6 - e_5 * fs_8_2431_1430 * r_2 * h9_m6 - e_5 * fs_90_46189_2002 * r_4 * h7_m6;
    }

    // NOTE: the angular components are formed in 124 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph1_m1, ph3_m1, ph5_m5, ph5_m1, ph7_m5, ph7_m1, ph9_m5, ph9_m1, ph11_m5, ph11_m1, ab_2, pc_108 : simd::cache_line_size())
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

        const auto h1_m1 = ph1_m1[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h5_m5 = ph5_m5[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h7_m5 = ph7_m5[k];
        const auto h7_m1 = ph7_m1[k];
        const auto h9_m5 = ph9_m5[k];
        const auto h9_m1 = ph9_m1[k];
        const auto h11_m5 = ph11_m5[k];
        const auto h11_m1 = ph11_m1[k];

        pc_108[k] = e_0 * fs_945_32_3 * h1_m1 + e_1 * fs_2205_32_2 * h3_m1 - e_1 * fs_945_16_3 * r_2 * h1_m1 - e_2 * fs_25_4_42 * h5_m5 - e_2 * fs_55_4_5 * h5_m1 - e_2 * fs_245_4_2 * r_2 * h3_m1 + e_2 * fs_135_4_3 * r_4 * h1_m1 + e_3 * fs_45_572_77 * h7_m5 - e_3 * fs_105_286_21 * h7_m1 + e_3 * fs_75_26_42 * r_2 * h5_m5 + e_3 * fs_165_26_5 * r_2 * h5_m1 + e_3 * fs_735_44_2 * r_4 * h3_m1 - e_3 * fs_15_2_3 * r_6 * h1_m1 + e_4 * fs_105_4862_858 * h9_m5 + e_4 * fs_1869_4862_15 * h9_m1 - e_4 * fs_45_2431_77 * r_2 * h7_m5 + e_4 * fs_210_2431_21 * r_2 * h7_m1 - e_4 * fs_5_13_42 * r_4 * h5_m5 - e_4 * fs_11_13_5 * r_4 * h5_m1 - e_4 * fs_245_143_2 * r_6 * h3_m1 + e_4 * fs_15_22_3 * r_8 * h1_m1 - e_5 * fs_45_4199_429 * h11_m5 + e_5 * fs_135_8398_22 * h11_m1 - e_5 * fs_5_2431_858 * r_2 * h9_m5 - e_5 * fs_89_2431_15 * r_2 * h9_m1 + e_5 * fs_45_46189_77 * r_4 * h7_m5 - e_5 * fs_210_46189_21 * r_4 * h7_m1 + e_5 * fs_10_663_42 * r_6 * h5_m5 + e_5 * fs_22_663_5 * r_6 * h5_m1 + e_5 * fs_49_858_2 * r_8 * h3_m1 - e_5 * fs_3_143_3 * r_10 * h1_m1;
    }

    // NOTE: the angular components are formed in 124 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph3_m2, ph5_m4, ph5_m2, ph7_m4, ph7_m2, ph9_m4, ph9_m2, ph11_m4, ph11_m2, ab_2, pc_109 : simd::cache_line_size())
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

        const auto h3_m2 = ph3_m2[k];
        const auto h5_m4 = ph5_m4[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h7_m4 = ph7_m4[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h9_m4 = ph9_m4[k];
        const auto h9_m2 = ph9_m2[k];
        const auto h11_m4 = ph11_m4[k];
        const auto h11_m2 = ph11_m2[k];

        pc_109[k] = - e_1 * fs_315_4_2 * h3_m2 - e_2 * fs_55_8_42 * h5_m4 + e_2 * fs_5_4_14 * h5_m2 + e_2 * fs_70_1_2 * r_2 * h3_m2 + e_3 * fs_135_572_770 * h7_m4 + e_3 * fs_315_286_35 * h7_m2 + e_3 * fs_165_52_42 * r_2 * h5_m4 - e_3 * fs_15_26_14 * r_2 * h5_m2 - e_3 * fs_210_11_2 * r_4 * h3_m2 - e_4 * fs_21_9724_6006 * h9_m4 - e_4 * fs_1113_4862_33 * h9_m2 - e_4 * fs_135_2431_770 * r_2 * h7_m4 - e_4 * fs_630_2431_35 * r_2 * h7_m2 - e_4 * fs_11_26_42 * r_4 * h5_m4 + e_4 * fs_1_13_14 * r_4 * h5_m2 + e_4 * fs_280_143_2 * r_6 * h3_m2 - e_5 * fs_9_8398_30030 * h11_m4 - e_5 * fs_27_4199_286 * h11_m2 + e_5 * fs_1_4862_6006 * r_2 * h9_m4 + e_5 * fs_53_2431_33 * r_2 * h9_m2 + e_5 * fs_135_46189_770 * r_4 * h7_m4 + e_5 * fs_630_46189_35 * r_4 * h7_m2 + e_5 * fs_11_663_42 * r_6 * h5_m4 - e_5 * fs_2_663_14 * r_6 * h5_m2 - e_5 * fs_28_429_2 * r_8 * h3_m2;
    }

    // NOTE: the angular components are formed in 124 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph3_p3, ph5_p3, ph7_p3, ph9_p3, ph11_p3, ab_2, pc_110 : simd::cache_line_size())
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

        const auto h3_p3 = ph3_p3[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h7_p3 = ph7_p3[k];
        const auto h9_p3 = ph9_p3[k];
        const auto h11_p3 = ph11_p3[k];

        pc_110[k] = - e_1 * fs_315_8_7 * h3_p3 - e_2 * f_145_4 * h5_p3 + e_2 * fs_35_1_7 * r_2 * h3_p3 + e_3 * fs_315_143_30 * h7_p3 + e_3 * f_435_26 * r_2 * h5_p3 - e_3 * fs_105_11_7 * r_4 * h3_p3 - e_4 * fs_1029_4862_33 * h9_p3 - e_4 * fs_1260_2431_30 * r_2 * h7_p3 - e_4 * f_29_13 * r_4 * h5_p3 + e_4 * fs_140_143_7 * r_6 * h3_p3 - e_5 * fs_42_4199_429 * h11_p3 + e_5 * fs_49_2431_33 * r_2 * h9_p3 + e_5 * fs_1260_46189_30 * r_4 * h7_p3 + e_5 * f_58_663 * r_6 * h5_p3 - e_5 * fs_14_429_7 * r_8 * h3_p3;
    }

    // NOTE: the angular components are formed in 124 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph3_p2, ph5_p2, ph5_p4, ph7_p2, ph7_p4, ph9_p2, ph9_p4, ph11_p2, ph11_p4, ab_2, pc_111 : simd::cache_line_size())
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

        const auto h3_p2 = ph3_p2[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h5_p4 = ph5_p4[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h7_p4 = ph7_p4[k];
        const auto h9_p2 = ph9_p2[k];
        const auto h9_p4 = ph9_p4[k];
        const auto h11_p2 = ph11_p2[k];
        const auto h11_p4 = ph11_p4[k];

        pc_111[k] = e_1 * fs_315_4_2 * h3_p2 - e_2 * fs_5_4_14 * h5_p2 - e_2 * fs_55_8_42 * h5_p4 - e_2 * fs_70_1_2 * r_2 * h3_p2 - e_3 * fs_315_286_35 * h7_p2 + e_3 * fs_135_572_770 * h7_p4 + e_3 * fs_15_26_14 * r_2 * h5_p2 + e_3 * fs_165_52_42 * r_2 * h5_p4 + e_3 * fs_210_11_2 * r_4 * h3_p2 + e_4 * fs_1113_4862_33 * h9_p2 - e_4 * fs_21_9724_6006 * h9_p4 + e_4 * fs_630_2431_35 * r_2 * h7_p2 - e_4 * fs_135_2431_770 * r_2 * h7_p4 - e_4 * fs_1_13_14 * r_4 * h5_p2 - e_4 * fs_11_26_42 * r_4 * h5_p4 - e_4 * fs_280_143_2 * r_6 * h3_p2 + e_5 * fs_27_4199_286 * h11_p2 - e_5 * fs_9_8398_30030 * h11_p4 - e_5 * fs_53_2431_33 * r_2 * h9_p2 + e_5 * fs_1_4862_6006 * r_2 * h9_p4 - e_5 * fs_630_46189_35 * r_4 * h7_p2 + e_5 * fs_135_46189_770 * r_4 * h7_p4 + e_5 * fs_2_663_14 * r_6 * h5_p2 + e_5 * fs_11_663_42 * r_6 * h5_p4 + e_5 * fs_28_429_2 * r_8 * h3_p2;
    }

    // NOTE: the angular components are formed in 124 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph1_p1, ph3_p1, ph5_p1, ph5_p5, ph7_p1, ph7_p5, ph9_p1, ph9_p5, ph11_p1, ph11_p5, ab_2, pc_112 : simd::cache_line_size())
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

        const auto h1_p1 = ph1_p1[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h5_p5 = ph5_p5[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h7_p5 = ph7_p5[k];
        const auto h9_p1 = ph9_p1[k];
        const auto h9_p5 = ph9_p5[k];
        const auto h11_p1 = ph11_p1[k];
        const auto h11_p5 = ph11_p5[k];

        pc_112[k] = - e_0 * fs_945_32_3 * h1_p1 - e_1 * fs_2205_32_2 * h3_p1 + e_1 * fs_945_16_3 * r_2 * h1_p1 + e_2 * fs_55_4_5 * h5_p1 - e_2 * fs_25_4_42 * h5_p5 + e_2 * fs_245_4_2 * r_2 * h3_p1 - e_2 * fs_135_4_3 * r_4 * h1_p1 + e_3 * fs_105_286_21 * h7_p1 + e_3 * fs_45_572_77 * h7_p5 - e_3 * fs_165_26_5 * r_2 * h5_p1 + e_3 * fs_75_26_42 * r_2 * h5_p5 - e_3 * fs_735_44_2 * r_4 * h3_p1 + e_3 * fs_15_2_3 * r_6 * h1_p1 - e_4 * fs_1869_4862_15 * h9_p1 + e_4 * fs_105_4862_858 * h9_p5 - e_4 * fs_210_2431_21 * r_2 * h7_p1 - e_4 * fs_45_2431_77 * r_2 * h7_p5 + e_4 * fs_11_13_5 * r_4 * h5_p1 - e_4 * fs_5_13_42 * r_4 * h5_p5 + e_4 * fs_245_143_2 * r_6 * h3_p1 - e_4 * fs_15_22_3 * r_8 * h1_p1 - e_5 * fs_135_8398_22 * h11_p1 - e_5 * fs_45_4199_429 * h11_p5 + e_5 * fs_89_2431_15 * r_2 * h9_p1 - e_5 * fs_5_2431_858 * r_2 * h9_p5 + e_5 * fs_210_46189_21 * r_4 * h7_p1 + e_5 * fs_45_46189_77 * r_4 * h7_p5 - e_5 * fs_22_663_5 * r_6 * h5_p1 + e_5 * fs_10_663_42 * r_6 * h5_p5 - e_5 * fs_49_858_2 * r_8 * h3_p1 + e_5 * fs_3_143_3 * r_10 * h1_p1;
    }

    // NOTE: the angular components are formed in 124 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph1_0, ph3_0, ph5_0, ph7_0, ph7_p6, ph9_0, ph9_p6, ph11_0, ph11_p6, ab_2, pc_113 : simd::cache_line_size())
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

        const auto h1_0 = ph1_0[k];
        const auto h3_0 = ph3_0[k];
        const auto h5_0 = ph5_0[k];
        const auto h7_0 = ph7_0[k];
        const auto h7_p6 = ph7_p6[k];
        const auto h9_0 = ph9_0[k];
        const auto h9_p6 = ph9_p6[k];
        const auto h11_0 = ph11_0[k];
        const auto h11_p6 = ph11_p6[k];

        pc_113[k] = e_0 * fs_2835_32_3 * h1_0 + e_1 * fs_315_16_3 * h3_0 - e_1 * fs_2835_16_3 * r_2 * h1_0 - e_2 * fs_125_4_3 * h5_0 - e_2 * fs_35_2_3 * r_2 * h3_0 + e_2 * fs_405_4_3 * r_4 * h1_0 + e_3 * fs_945_286_3 * h7_0 - e_3 * fs_45_286_2002 * h7_p6 + e_3 * fs_375_26_3 * r_2 * h5_0 + e_3 * fs_105_22_3 * r_4 * h3_0 - e_3 * fs_45_2_3 * r_6 * h1_0 + e_4 * fs_5565_4862_3 * h9_0 + e_4 * fs_84_2431_1430 * h9_p6 - e_4 * fs_1890_2431_3 * r_2 * h7_0 + e_4 * fs_90_2431_2002 * r_2 * h7_p6 - e_4 * fs_25_13_3 * r_4 * h5_0 - e_4 * fs_70_143_3 * r_6 * h3_0 + e_4 * fs_45_22_3 * r_8 * h1_0 + e_5 * fs_165_4199_3 * h11_0 - e_5 * fs_15_4199_4862 * h11_p6 - e_5 * fs_265_2431_3 * r_2 * h9_0 - e_5 * fs_8_2431_1430 * r_2 * h9_p6 + e_5 * fs_1890_46189_3 * r_4 * h7_0 - e_5 * fs_90_46189_2002 * r_4 * h7_p6 + e_5 * fs_50_663_3 * r_6 * h5_0 + e_5 * fs_7_429_3 * r_8 * h3_0 - e_5 * fs_9_143_3 * r_10 * h1_0;
    }

    // NOTE: the angular components are formed in 124 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph1_p1, ph3_p1, ph5_p1, ph7_p1, ph7_p7, ph9_p1, ph9_p7, ph11_p1, ph11_p7, ab_2, pc_114 : simd::cache_line_size())
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

        const auto h1_p1 = ph1_p1[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h7_p7 = ph7_p7[k];
        const auto h9_p1 = ph9_p1[k];
        const auto h9_p7 = ph9_p7[k];
        const auto h11_p1 = ph11_p1[k];
        const auto h11_p7 = ph11_p7[k];

        pc_114[k] = e_0 * fs_2835_64_10 * h1_p1 - e_1 * fs_315_16_15 * h3_p1 - e_1 * fs_2835_32_10 * r_2 * h1_p1 - e_2 * fs_65_8_6 * h5_p1 + e_2 * fs_35_2_15 * r_2 * h3_p1 + e_2 * fs_405_8_10 * r_4 * h1_p1 + e_3 * fs_45_44_70 * h7_p1 - e_3 * fs_105_572_4290 * h7_p7 + e_3 * fs_15_4_6 * r_2 * h5_p1 - e_3 * fs_105_22_15 * r_4 * h3_p1 - e_3 * fs_45_4_10 * r_6 * h1_p1 + e_4 * fs_7707_9724_2 * h9_p1 + e_4 * fs_609_4862_143 * h9_p7 - e_4 * fs_45_187_70 * r_2 * h7_p1 + e_4 * fs_105_2431_4290 * r_2 * h7_p7 - e_4 * fs_1_2_6 * r_4 * h5_p1 + e_4 * fs_70_143_15 * r_6 * h3_p1 + e_4 * fs_45_44_10 * r_8 * h1_p1 + e_5 * fs_9_4199_165 * h11_p1 - e_5 * fs_9_4199_14586 * h11_p7 - e_5 * fs_367_4862_2 * r_2 * h9_p1 - e_5 * fs_29_2431_143 * r_2 * h9_p7 + e_5 * fs_45_3553_70 * r_4 * h7_p1 - e_5 * fs_105_46189_4290 * r_4 * h7_p7 + e_5 * fs_1_51_6 * r_6 * h5_p1 - e_5 * fs_7_429_15 * r_8 * h3_p1 - e_5 * fs_9_286_10 * r_10 * h1_p1;
    }

    // NOTE: the angular components are formed in 124 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph3_p2, ph5_p2, ph7_p2, ph9_p2, ph9_p8, ph11_p2, ph11_p8, ab_2, pc_115 : simd::cache_line_size())
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

        const auto h3_p2 = ph3_p2[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h9_p2 = ph9_p2[k];
        const auto h9_p8 = ph9_p8[k];
        const auto h11_p2 = ph11_p2[k];
        const auto h11_p8 = ph11_p8[k];

        pc_115[k] = - e_1 * fs_315_16_33 * h3_p2 + e_2 * fs_5_4_231 * h5_p2 + e_2 * fs_35_2_33 * r_2 * h3_p2 + e_3 * fs_30_143_2310 * h7_p2 - e_3 * fs_15_26_231 * r_2 * h5_p2 - e_3 * fs_105_22_33 * r_4 * h3_p2 + e_4 * fs_231_442_2 * h9_p2 + e_4 * fs_21_442_221 * h9_p8 - e_4 * fs_120_2431_2310 * r_2 * h7_p2 + e_4 * fs_1_13_231 * r_4 * h5_p2 + e_4 * fs_70_143_33 * r_6 * h3_p2 + e_5 * fs_9_4199_39 * h11_p2 - e_5 * fs_9_4199_12597 * h11_p8 - e_5 * fs_11_221_2 * r_2 * h9_p2 - e_5 * fs_1_221_221 * r_2 * h9_p8 + e_5 * fs_120_46189_2310 * r_4 * h7_p2 - e_5 * fs_2_663_231 * r_6 * h5_p2 - e_5 * fs_7_429_33 * r_8 * h3_p2;
    }

    // NOTE: the angular components are formed in 124 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph3_p3, ph5_p3, ph7_p3, ph9_p3, ph9_p9, ph11_p3, ph11_p9, ab_2, pc_116 : simd::cache_line_size())
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

        const auto h3_p3 = ph3_p3[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h7_p3 = ph7_p3[k];
        const auto h9_p3 = ph9_p3[k];
        const auto h9_p9 = ph9_p9[k];
        const auto h11_p3 = ph11_p3[k];
        const auto h11_p9 = ph11_p9[k];

        pc_116[k] = e_1 * fs_315_32_66 * h3_p3 + e_2 * fs_5_2_462 * h5_p3 - e_2 * fs_35_4_66 * r_2 * h3_p3 + e_3 * fs_225_572_385 * h7_p3 - e_3 * fs_15_13_462 * r_2 * h5_p3 + e_3 * fs_105_44_66 * r_4 * h3_p3 + e_4 * fs_21_221_14 * h9_p3 - e_4 * fs_21_442_1326 * h9_p9 - e_4 * fs_225_2431_385 * r_2 * h7_p3 + e_4 * fs_2_13_462 * r_4 * h5_p3 - e_4 * fs_35_143_66 * r_6 * h3_p3 + e_5 * fs_3_8398_182 * h11_p3 - e_5 * fs_3_4199_62985 * h11_p9 - e_5 * fs_2_221_14 * r_2 * h9_p3 + e_5 * fs_1_221_1326 * r_2 * h9_p9 + e_5 * fs_225_46189_385 * r_4 * h7_p3 - e_5 * fs_4_663_462 * r_6 * h5_p3 + e_5 * fs_7_858_66 * r_8 * h3_p3;
    }

    // NOTE: the angular components are formed in 124 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph3_m2, ph5_m2, ph7_m2, ph9_m2, ph11_m10, ph11_m2, ab_2, pc_117 : simd::cache_line_size())
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

        const auto h3_m2 = ph3_m2[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h9_m2 = ph9_m2[k];
        const auto h11_m10 = ph11_m10[k];
        const auto h11_m2 = ph11_m2[k];

        pc_117[k] = - e_1 * fs_945_32_22 * h3_m2 - e_2 * fs_15_4_154 * h5_m2 + e_2 * fs_105_4_22 * r_2 * h3_m2 - e_3 * fs_135_572_385 * h7_m2 + e_3 * fs_45_26_154 * r_2 * h5_m2 - e_3 * fs_315_44_22 * r_4 * h3_m2 - e_4 * fs_21_221_3 * h9_m2 + e_4 * fs_135_2431_385 * r_2 * h7_m2 - e_4 * fs_3_13_154 * r_4 * h5_m2 + e_4 * fs_105_143_22 * r_6 * h3_m2 - e_5 * fs_3_4199_146965 * h11_m10 - e_5 * fs_3_8398_26 * h11_m2 + e_5 * fs_2_221_3 * r_2 * h9_m2 - e_5 * fs_135_46189_385 * r_4 * h7_m2 + e_5 * fs_2_221_154 * r_6 * h5_m2 - e_5 * fs_7_286_22 * r_8 * h3_m2;
    }

    // NOTE: the angular components are formed in 124 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph1_m1, ph5_m1, ph7_m1, ph9_m9, ph9_m8, ph9_m1, ph11_m9, ph11_m8, ph11_m1, ab_2, pc_118, pc_119 : simd::cache_line_size())
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

        const auto h1_m1 = ph1_m1[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h7_m1 = ph7_m1[k];
        const auto h9_m9 = ph9_m9[k];
        const auto h9_m8 = ph9_m8[k];
        const auto h9_m1 = ph9_m1[k];
        const auto h11_m9 = ph11_m9[k];
        const auto h11_m8 = ph11_m8[k];
        const auto h11_m1 = ph11_m1[k];

        pc_118[k] = e_0 * fs_945_64_110 * h1_m1 - e_1 * fs_945_32_110 * r_2 * h1_m1 - e_2 * fs_45_8_66 * h5_m1 + e_2 * fs_135_8_110 * r_4 * h1_m1 - e_3 * fs_15_52_770 * h7_m1 + e_3 * fs_135_52_66 * r_2 * h5_m1 - e_3 * fs_15_4_110 * r_6 * h1_m1 + e_4 * fs_63_442_221 * h9_m9 - e_4 * fs_63_748_22 * h9_m1 + e_4 * fs_15_221_770 * r_2 * h7_m1 - e_4 * fs_9_26_66 * r_4 * h5_m1 + e_4 * fs_15_44_110 * r_8 * h1_m1 - e_5 * fs_6_4199_41990 * h11_m9 - e_5 * fs_6_4199_15 * h11_m1 - e_5 * fs_3_221_221 * r_2 * h9_m9 + e_5 * fs_3_374_22 * r_2 * h9_m1 - e_5 * fs_15_4199_770 * r_4 * h7_m1 + e_5 * fs_3_221_66 * r_6 * h5_m1 - e_5 * fs_3_286_110 * r_10 * h1_m1;

        pc_119[k] = e_4 * fs_63_2431_2431 * h9_m8 - e_5 * fs_3_4199_138567 * h11_m8 - e_5 * fs_6_2431_2431 * r_2 * h9_m8;
    }

    // NOTE: the angular components are formed in 124 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph1_m1, ph3_m1, ph5_m1, ph7_m7, ph7_m1, ph9_m7, ph9_m1, ph11_m7, ph11_m1, ab_2, pc_120 : simd::cache_line_size())
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

        const auto h1_m1 = ph1_m1[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h7_m7 = ph7_m7[k];
        const auto h7_m1 = ph7_m1[k];
        const auto h9_m7 = ph9_m7[k];
        const auto h9_m1 = ph9_m1[k];
        const auto h11_m7 = ph11_m7[k];
        const auto h11_m1 = ph11_m1[k];

        pc_120[k] = e_0 * fs_945_64_6 * h1_m1 + e_1 * f_945_8 * h3_m1 - e_1 * fs_945_32_6 * r_2 * h1_m1 + e_2 * fs_15_8_10 * h5_m1 - e_2 * f_105_1 * r_2 * h3_m1 + e_2 * fs_135_8_6 * r_4 * h1_m1 + e_3 * fs_315_572_286 * h7_m7 - e_3 * fs_210_143_42 * h7_m1 - e_3 * fs_45_52_10 * r_2 * h5_m1 + e_3 * f_315_11 * r_4 * h3_m1 - e_3 * fs_15_4_6 * r_6 * h1_m1 + e_4 * fs_21_4862_2145 * h9_m7 - e_4 * fs_1827_9724_30 * h9_m1 - e_4 * fs_315_2431_286 * r_2 * h7_m7 + e_4 * fs_840_2431_42 * r_2 * h7_m1 + e_4 * fs_3_26_10 * r_4 * h5_m1 - e_4 * f_420_143 * r_6 * h3_m1 + e_4 * fs_15_44_6 * r_8 * h1_m1 - e_5 * fs_6_4199_24310 * h11_m7 - e_5 * fs_30_4199_11 * h11_m1 - e_5 * fs_1_2431_2145 * r_2 * h9_m7 + e_5 * fs_87_4862_30 * r_2 * h9_m1 + e_5 * fs_315_46189_286 * r_4 * h7_m7 - e_5 * fs_840_46189_42 * r_4 * h7_m1 - e_5 * fs_1_221_10 * r_6 * h5_m1 + e_5 * f_14_143 * r_8 * h3_m1 - e_5 * fs_3_286_6 * r_10 * h1_m1;
    }

    // NOTE: the angular components are formed in 124 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph3_m2, ph5_m2, ph7_m6, ph7_m2, ph9_m6, ph9_m2, ph11_m6, ph11_m2, ab_2, pc_121 : simd::cache_line_size())
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

        const auto h3_m2 = ph3_m2[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h7_m6 = ph7_m6[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h9_m6 = ph9_m6[k];
        const auto h9_m2 = ph9_m2[k];
        const auto h11_m6 = ph11_m6[k];
        const auto h11_m2 = ph11_m2[k];

        pc_121[k] = - e_1 * fs_945_32_10 * h3_m2 - e_2 * fs_15_4_70 * h5_m2 + e_2 * fs_105_4_10 * r_2 * h3_m2 + e_3 * fs_15_44_1001 * h7_m6 + e_3 * fs_735_286_7 * h7_m2 + e_3 * fs_45_26_70 * r_2 * h5_m2 - e_3 * fs_315_44_10 * r_4 * h3_m2 - e_4 * fs_63_2431_715 * h9_m6 + e_4 * fs_252_2431_165 * h9_m2 - e_4 * fs_15_187_1001 * r_2 * h7_m6 - e_4 * fs_1470_2431_7 * r_2 * h7_m2 - e_4 * fs_3_13_70 * r_4 * h5_m2 + e_4 * fs_105_143_10 * r_6 * h3_m2 - e_5 * fs_15_4199_2431 * h11_m6 + e_5 * fs_9_8398_1430 * h11_m2 + e_5 * fs_6_2431_715 * r_2 * h9_m6 - e_5 * fs_24_2431_165 * r_2 * h9_m2 + e_5 * fs_15_3553_1001 * r_4 * h7_m6 + e_5 * fs_1470_46189_7 * r_4 * h7_m2 + e_5 * fs_2_221_70 * r_6 * h5_m2 - e_5 * fs_7_286_10 * r_8 * h3_m2;
    }

    // NOTE: the angular components are formed in 124 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph3_m3, ph5_m5, ph5_m3, ph7_m5, ph7_m3, ph9_m5, ph9_m3, ph11_m5, ph11_m3, ab_2, pc_122 : simd::cache_line_size())
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

        const auto h3_m3 = ph3_m3[k];
        const auto h5_m5 = ph5_m5[k];
        const auto h5_m3 = ph5_m3[k];
        const auto h7_m5 = ph7_m5[k];
        const auto h7_m3 = ph7_m3[k];
        const auto h9_m5 = ph9_m5[k];
        const auto h9_m3 = ph9_m3[k];
        const auto h11_m5 = ph11_m5[k];
        const auto h11_m3 = ph11_m3[k];

        pc_122[k] = e_1 * fs_315_16_6 * h3_m3 + e_2 * fs_15_8_210 * h5_m5 + e_2 * fs_55_8_42 * h5_m3 - e_2 * fs_35_2_6 * r_2 * h3_m3 + e_3 * fs_60_143_385 * h7_m5 - e_3 * fs_105_286_35 * h7_m3 - e_3 * fs_45_52_210 * r_2 * h5_m5 - e_3 * fs_165_52_42 * r_2 * h5_m3 + e_3 * fs_105_22_6 * r_4 * h3_m3 - e_4 * fs_189_9724_4290 * h9_m5 - e_4 * fs_1197_9724_154 * h9_m3 - e_4 * fs_240_2431_385 * r_2 * h7_m5 + e_4 * fs_210_2431_35 * r_2 * h7_m3 + e_4 * fs_3_26_210 * r_4 * h5_m5 + e_4 * fs_11_26_42 * r_4 * h5_m3 - e_4 * fs_70_143_6 * r_6 * h3_m3 - e_5 * fs_12_4199_2145 * h11_m5 - e_5 * fs_6_4199_2002 * h11_m3 + e_5 * fs_9_4862_4290 * r_2 * h9_m5 + e_5 * fs_57_4862_154 * r_2 * h9_m3 + e_5 * fs_240_46189_385 * r_4 * h7_m5 - e_5 * fs_210_46189_35 * r_4 * h7_m3 - e_5 * fs_1_221_210 * r_6 * h5_m5 - e_5 * fs_11_663_42 * r_6 * h5_m3 + e_5 * fs_7_429_6 * r_8 * h3_m3;
    }

    // NOTE: the angular components are formed in 124 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_2, pe_3, pe_4, pe_5, ph5_p4, ph7_p4, ph9_p4, ph11_p4, ab_2, pc_123 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h5_p4 = ph5_p4[k];
        const auto h7_p4 = ph7_p4[k];
        const auto h9_p4 = ph9_p4[k];
        const auto h11_p4 = ph11_p4[k];

        pc_123[k] = e_2 * f_60_1 * h5_p4 + e_3 * fs_105_286_165 * h7_p4 - e_3 * f_360_13 * r_2 * h5_p4 - e_4 * fs_441_2431_143 * h9_p4 - e_4 * fs_210_2431_165 * r_2 * h7_p4 + e_4 * f_48_13 * r_4 * h5_p4 - e_5 * fs_21_4199_715 * h11_p4 + e_5 * fs_42_2431_143 * r_2 * h9_p4 + e_5 * fs_210_46189_165 * r_4 * h7_p4 - e_5 * f_32_221 * r_6 * h5_p4;
    }

    // NOTE: the angular components are formed in 124 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph3_p3, ph5_p3, ph5_p5, ph7_p3, ph7_p5, ph9_p3, ph9_p5, ph11_p3, ph11_p5, ab_2, pc_124 : simd::cache_line_size())
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

        const auto h3_p3 = ph3_p3[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h5_p5 = ph5_p5[k];
        const auto h7_p3 = ph7_p3[k];
        const auto h7_p5 = ph7_p5[k];
        const auto h9_p3 = ph9_p3[k];
        const auto h9_p5 = ph9_p5[k];
        const auto h11_p3 = ph11_p3[k];
        const auto h11_p5 = ph11_p5[k];

        pc_124[k] = - e_1 * fs_315_16_6 * h3_p3 - e_2 * fs_55_8_42 * h5_p3 + e_2 * fs_15_8_210 * h5_p5 + e_2 * fs_35_2_6 * r_2 * h3_p3 + e_3 * fs_105_286_35 * h7_p3 + e_3 * fs_60_143_385 * h7_p5 + e_3 * fs_165_52_42 * r_2 * h5_p3 - e_3 * fs_45_52_210 * r_2 * h5_p5 - e_3 * fs_105_22_6 * r_4 * h3_p3 + e_4 * fs_1197_9724_154 * h9_p3 - e_4 * fs_189_9724_4290 * h9_p5 - e_4 * fs_210_2431_35 * r_2 * h7_p3 - e_4 * fs_240_2431_385 * r_2 * h7_p5 - e_4 * fs_11_26_42 * r_4 * h5_p3 + e_4 * fs_3_26_210 * r_4 * h5_p5 + e_4 * fs_70_143_6 * r_6 * h3_p3 + e_5 * fs_6_4199_2002 * h11_p3 - e_5 * fs_12_4199_2145 * h11_p5 - e_5 * fs_57_4862_154 * r_2 * h9_p3 + e_5 * fs_9_4862_4290 * r_2 * h9_p5 + e_5 * fs_210_46189_35 * r_4 * h7_p3 + e_5 * fs_240_46189_385 * r_4 * h7_p5 + e_5 * fs_11_663_42 * r_6 * h5_p3 - e_5 * fs_1_221_210 * r_6 * h5_p5 - e_5 * fs_7_429_6 * r_8 * h3_p3;
    }

    // NOTE: the angular components are formed in 124 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph3_p2, ph5_p2, ph7_p2, ph7_p6, ph9_p2, ph9_p6, ph11_p2, ph11_p6, ab_2, pc_125 : simd::cache_line_size())
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

        const auto h3_p2 = ph3_p2[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h7_p6 = ph7_p6[k];
        const auto h9_p2 = ph9_p2[k];
        const auto h9_p6 = ph9_p6[k];
        const auto h11_p2 = ph11_p2[k];
        const auto h11_p6 = ph11_p6[k];

        pc_125[k] = e_1 * fs_945_32_10 * h3_p2 + e_2 * fs_15_4_70 * h5_p2 - e_2 * fs_105_4_10 * r_2 * h3_p2 - e_3 * fs_735_286_7 * h7_p2 + e_3 * fs_15_44_1001 * h7_p6 - e_3 * fs_45_26_70 * r_2 * h5_p2 + e_3 * fs_315_44_10 * r_4 * h3_p2 - e_4 * fs_252_2431_165 * h9_p2 - e_4 * fs_63_2431_715 * h9_p6 + e_4 * fs_1470_2431_7 * r_2 * h7_p2 - e_4 * fs_15_187_1001 * r_2 * h7_p6 + e_4 * fs_3_13_70 * r_4 * h5_p2 - e_4 * fs_105_143_10 * r_6 * h3_p2 - e_5 * fs_9_8398_1430 * h11_p2 - e_5 * fs_15_4199_2431 * h11_p6 + e_5 * fs_24_2431_165 * r_2 * h9_p2 + e_5 * fs_6_2431_715 * r_2 * h9_p6 - e_5 * fs_1470_46189_7 * r_4 * h7_p2 + e_5 * fs_15_3553_1001 * r_4 * h7_p6 - e_5 * fs_2_221_70 * r_6 * h5_p2 + e_5 * fs_7_286_10 * r_8 * h3_p2;
    }

    // NOTE: the angular components are formed in 124 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph1_p1, ph3_p1, ph5_p1, ph7_p1, ph7_p7, ph9_p1, ph9_p7, ph11_p1, ph11_p7, ab_2, pc_126 : simd::cache_line_size())
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

        const auto h1_p1 = ph1_p1[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h7_p7 = ph7_p7[k];
        const auto h9_p1 = ph9_p1[k];
        const auto h9_p7 = ph9_p7[k];
        const auto h11_p1 = ph11_p1[k];
        const auto h11_p7 = ph11_p7[k];

        pc_126[k] = - e_0 * fs_945_64_6 * h1_p1 - e_1 * f_945_8 * h3_p1 + e_1 * fs_945_32_6 * r_2 * h1_p1 - e_2 * fs_15_8_10 * h5_p1 + e_2 * f_105_1 * r_2 * h3_p1 - e_2 * fs_135_8_6 * r_4 * h1_p1 + e_3 * fs_210_143_42 * h7_p1 + e_3 * fs_315_572_286 * h7_p7 + e_3 * fs_45_52_10 * r_2 * h5_p1 - e_3 * f_315_11 * r_4 * h3_p1 + e_3 * fs_15_4_6 * r_6 * h1_p1 + e_4 * fs_1827_9724_30 * h9_p1 + e_4 * fs_21_4862_2145 * h9_p7 - e_4 * fs_840_2431_42 * r_2 * h7_p1 - e_4 * fs_315_2431_286 * r_2 * h7_p7 - e_4 * fs_3_26_10 * r_4 * h5_p1 + e_4 * f_420_143 * r_6 * h3_p1 - e_4 * fs_15_44_6 * r_8 * h1_p1 + e_5 * fs_30_4199_11 * h11_p1 - e_5 * fs_6_4199_24310 * h11_p7 - e_5 * fs_87_4862_30 * r_2 * h9_p1 - e_5 * fs_1_2431_2145 * r_2 * h9_p7 + e_5 * fs_840_46189_42 * r_4 * h7_p1 + e_5 * fs_315_46189_286 * r_4 * h7_p7 + e_5 * fs_1_221_10 * r_6 * h5_p1 - e_5 * f_14_143 * r_8 * h3_p1 + e_5 * fs_3_286_6 * r_10 * h1_p1;
    }

    // NOTE: the angular components are formed in 124 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph1_0, ph3_0, ph5_0, ph7_0, ph9_0, ph9_p8, ph11_0, ph11_p8, ab_2, pc_127 : simd::cache_line_size())
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

        const auto h1_0 = ph1_0[k];
        const auto h3_0 = ph3_0[k];
        const auto h5_0 = ph5_0[k];
        const auto h7_0 = ph7_0[k];
        const auto h9_0 = ph9_0[k];
        const auto h9_p8 = ph9_p8[k];
        const auto h11_0 = ph11_0[k];
        const auto h11_p8 = ph11_p8[k];

        pc_127[k] = e_0 * fs_945_16_5 * h1_0 + e_1 * fs_945_16_5 * h3_0 - e_1 * fs_945_8_5 * r_2 * h1_0 - e_2 * fs_15_1_5 * h5_0 - e_2 * fs_105_2_5 * r_2 * h3_0 + e_2 * fs_135_2_5 * r_4 * h1_0 - e_3 * fs_1785_286_5 * h7_0 + e_3 * fs_90_13_5 * r_2 * h5_0 + e_3 * fs_315_22_5 * r_4 * h3_0 - e_3 * fs_15_1_5 * r_6 * h1_0 - e_4 * fs_63_143_5 * h9_0 + e_4 * fs_63_2431_2431 * h9_p8 + e_4 * fs_210_143_5 * r_2 * h7_0 - e_4 * fs_12_13_5 * r_4 * h5_0 - e_4 * fs_210_143_5 * r_6 * h3_0 + e_4 * fs_15_11_5 * r_8 * h1_0 - e_5 * fs_33_4199_5 * h11_0 - e_5 * fs_3_4199_138567 * h11_p8 + e_5 * fs_6_143_5 * r_2 * h9_0 - e_5 * fs_6_2431_2431 * r_2 * h9_p8 - e_5 * fs_210_2717_5 * r_4 * h7_0 + e_5 * fs_8_221_5 * r_6 * h5_0 + e_5 * fs_7_143_5 * r_8 * h3_0 - e_5 * fs_6_143_5 * r_10 * h1_0;
    }

    // NOTE: the angular components are formed in 124 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph1_p1, ph5_p1, ph7_p1, ph9_p1, ph9_p9, ph11_p1, ph11_p9, ab_2, pc_128 : simd::cache_line_size())
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

        const auto h1_p1 = ph1_p1[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h9_p1 = ph9_p1[k];
        const auto h9_p9 = ph9_p9[k];
        const auto h11_p1 = ph11_p1[k];
        const auto h11_p9 = ph11_p9[k];

        pc_128[k] = e_0 * fs_945_64_110 * h1_p1 - e_1 * fs_945_32_110 * r_2 * h1_p1 - e_2 * fs_45_8_66 * h5_p1 + e_2 * fs_135_8_110 * r_4 * h1_p1 - e_3 * fs_15_52_770 * h7_p1 + e_3 * fs_135_52_66 * r_2 * h5_p1 - e_3 * fs_15_4_110 * r_6 * h1_p1 - e_4 * fs_63_748_22 * h9_p1 + e_4 * fs_63_442_221 * h9_p9 + e_4 * fs_15_221_770 * r_2 * h7_p1 - e_4 * fs_9_26_66 * r_4 * h5_p1 + e_4 * fs_15_44_110 * r_8 * h1_p1 - e_5 * fs_6_4199_15 * h11_p1 - e_5 * fs_6_4199_41990 * h11_p9 + e_5 * fs_3_374_22 * r_2 * h9_p1 - e_5 * fs_3_221_221 * r_2 * h9_p9 - e_5 * fs_15_4199_770 * r_4 * h7_p1 + e_5 * fs_3_221_66 * r_6 * h5_p1 - e_5 * fs_3_286_110 * r_10 * h1_p1;
    }

    // NOTE: the angular components are formed in 124 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph3_p2, ph5_p2, ph7_p2, ph9_p2, ph11_p2, ph11_p10, ab_2, pc_129 : simd::cache_line_size())
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

        const auto h3_p2 = ph3_p2[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h9_p2 = ph9_p2[k];
        const auto h11_p2 = ph11_p2[k];
        const auto h11_p10 = ph11_p10[k];

        pc_129[k] = - e_1 * fs_945_32_22 * h3_p2 - e_2 * fs_15_4_154 * h5_p2 + e_2 * fs_105_4_22 * r_2 * h3_p2 - e_3 * fs_135_572_385 * h7_p2 + e_3 * fs_45_26_154 * r_2 * h5_p2 - e_3 * fs_315_44_22 * r_4 * h3_p2 - e_4 * fs_21_221_3 * h9_p2 + e_4 * fs_135_2431_385 * r_2 * h7_p2 - e_4 * fs_3_13_154 * r_4 * h5_p2 + e_4 * fs_105_143_22 * r_6 * h3_p2 - e_5 * fs_3_8398_26 * h11_p2 - e_5 * fs_3_4199_146965 * h11_p10 + e_5 * fs_2_221_3 * r_2 * h9_p2 - e_5 * fs_135_46189_385 * r_4 * h7_p2 + e_5 * fs_2_221_154 * r_6 * h5_p2 - e_5 * fs_7_286_22 * r_8 * h3_p2;
    }

    // NOTE: the angular components are formed in 124 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph1_m1, ph3_m1, ph5_m1, ph7_m1, ph9_m9, ph9_m1, ph11_m11, ph11_m10, ph11_m9, ph11_m1, ab_2, pc_130, pc_131, pc_132 : simd::cache_line_size())
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

        const auto h1_m1 = ph1_m1[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h7_m1 = ph7_m1[k];
        const auto h9_m9 = ph9_m9[k];
        const auto h9_m1 = ph9_m1[k];
        const auto h11_m11 = ph11_m11[k];
        const auto h11_m10 = ph11_m10[k];
        const auto h11_m9 = ph11_m9[k];
        const auto h11_m1 = ph11_m1[k];

        pc_130[k] = e_0 * fs_945_32_33 * h1_m1 + e_1 * fs_945_32_22 * h3_m1 - e_1 * fs_945_16_33 * r_2 * h1_m1 + e_2 * fs_15_4_55 * h5_m1 - e_2 * fs_105_4_22 * r_2 * h3_m1 + e_2 * fs_135_4_33 * r_4 * h1_m1 + e_3 * fs_75_572_231 * h7_m1 - e_3 * fs_45_26_55 * r_2 * h5_m1 + e_3 * fs_315_44_22 * r_4 * h3_m1 - e_3 * fs_15_2_33 * r_6 * h1_m1 + e_4 * fs_21_4862_165 * h9_m1 - e_4 * fs_75_2431_231 * r_2 * h7_m1 + e_4 * fs_3_13_55 * r_4 * h5_m1 - e_4 * fs_105_143_22 * r_6 * h3_m1 + e_4 * fs_15_22_33 * r_8 * h1_m1 - e_5 * fs_3_4199_323323 * h11_m11 + e_5 * fs_3_8398_2 * h11_m1 - e_5 * fs_1_2431_165 * r_2 * h9_m1 + e_5 * fs_75_46189_231 * r_4 * h7_m1 - e_5 * fs_2_221_55 * r_6 * h5_m1 + e_5 * fs_7_286_22 * r_8 * h3_m1 - e_5 * fs_3_143_33 * r_10 * h1_m1;

        pc_131[k] = - e_5 * fs_3_4199_176358 * h11_m10;

        pc_132[k] = e_0 * fs_945_64_2 * h1_m1 + e_1 * fs_945_16_3 * h3_m1 - e_1 * fs_945_32_2 * r_2 * h1_m1 + e_2 * fs_75_8_30 * h5_m1 - e_2 * fs_105_2_3 * r_2 * h3_m1 + e_2 * fs_135_8_2 * r_4 * h1_m1 + e_3 * fs_525_286_14 * h7_m1 - e_3 * fs_225_52_30 * r_2 * h5_m1 + e_3 * fs_315_22_3 * r_4 * h3_m1 - e_3 * fs_15_4_2 * r_6 * h1_m1 - e_4 * fs_63_4862_12155 * h9_m9 + e_4 * fs_945_9724_10 * h9_m1 - e_4 * fs_1050_2431_14 * r_2 * h7_m1 + e_4 * fs_15_26_30 * r_4 * h5_m1 - e_4 * fs_210_143_3 * r_6 * h3_m1 + e_4 * fs_15_44_2 * r_8 * h1_m1 - e_5 * fs_3_4199_92378 * h11_m9 + e_5 * fs_3_4199_33 * h11_m1 + e_5 * fs_3_2431_12155 * r_2 * h9_m9 - e_5 * fs_45_4862_10 * r_2 * h9_m1 + e_5 * fs_1050_46189_14 * r_4 * h7_m1 - e_5 * fs_5_221_30 * r_6 * h5_m1 + e_5 * fs_7_143_3 * r_8 * h3_m1 - e_5 * fs_3_286_2 * r_10 * h1_m1;
    }

    // NOTE: the angular components are formed in 124 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph3_m2, ph5_m2, ph7_m2, ph9_m8, ph9_m2, ph11_m8, ph11_m2, ab_2, pc_133 : simd::cache_line_size())
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

        const auto h3_m2 = ph3_m2[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h9_m8 = ph9_m8[k];
        const auto h9_m2 = ph9_m2[k];
        const auto h11_m8 = ph11_m8[k];
        const auto h11_m2 = ph11_m2[k];

        pc_133[k] = - e_1 * f_945_16 * h3_m2 - e_2 * fs_75_4_7 * h5_m2 + e_2 * f_105_2 * r_2 * h3_m2 - e_3 * fs_315_286_70 * h7_m2 + e_3 * fs_225_26_7 * r_2 * h5_m2 - e_3 * f_315_22 * r_4 * h3_m2 - e_4 * fs_105_4862_7293 * h9_m8 - e_4 * fs_315_4862_66 * h9_m2 + e_4 * fs_630_2431_70 * r_2 * h7_m2 - e_4 * fs_15_13_7 * r_4 * h5_m2 + e_4 * f_210_143 * r_6 * h3_m2 - e_5 * fs_3_4199_46189 * h11_m8 - e_5 * fs_3_4199_143 * h11_m2 + e_5 * fs_5_2431_7293 * r_2 * h9_m8 + e_5 * fs_15_2431_66 * r_2 * h9_m2 - e_5 * fs_630_46189_70 * r_4 * h7_m2 + e_5 * fs_10_221_7 * r_6 * h5_m2 - e_5 * f_7_143 * r_8 * h3_m2;
    }

    // NOTE: the angular components are formed in 124 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph3_m3, ph5_m3, ph7_m7, ph7_m3, ph9_m7, ph9_m3, ph11_m7, ph11_m3, ab_2, pc_134 : simd::cache_line_size())
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

        const auto h3_m3 = ph3_m3[k];
        const auto h5_m3 = ph5_m3[k];
        const auto h7_m7 = ph7_m7[k];
        const auto h7_m3 = ph7_m3[k];
        const auto h9_m7 = ph9_m7[k];
        const auto h9_m3 = ph9_m3[k];
        const auto h11_m7 = ph11_m7[k];
        const auto h11_m3 = ph11_m3[k];

        pc_134[k] = e_1 * fs_315_32_6 * h3_m3 + e_2 * fs_25_4_42 * h5_m3 - e_2 * fs_35_4_6 * r_2 * h3_m3 - e_3 * fs_105_572_715 * h7_m7 + e_3 * fs_525_286_35 * h7_m3 - e_3 * fs_75_26_42 * r_2 * h5_m3 + e_3 * fs_105_44_6 * r_4 * h3_m3 - e_4 * fs_315_4862_858 * h9_m7 + e_4 * fs_315_4862_154 * h9_m3 + e_4 * fs_105_2431_715 * r_2 * h7_m7 - e_4 * fs_1050_2431_35 * r_2 * h7_m3 + e_4 * fs_5_13_42 * r_4 * h5_m3 - e_4 * fs_35_143_6 * r_6 * h3_m3 - e_5 * fs_9_4199_2431 * h11_m7 + e_5 * fs_3_8398_2002 * h11_m3 + e_5 * fs_15_2431_858 * r_2 * h9_m7 - e_5 * fs_15_2431_154 * r_2 * h9_m3 - e_5 * fs_105_46189_715 * r_4 * h7_m7 + e_5 * fs_1050_46189_35 * r_4 * h7_m3 - e_5 * fs_10_663_42 * r_6 * h5_m3 + e_5 * fs_7_858_6 * r_8 * h3_m3;
    }

    // NOTE: the angular components are formed in 124 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_2, pe_3, pe_4, pe_5, ph5_m4, ph5_p5, ph7_m6, ph7_m4, ph7_p5, ph9_m6, ph9_m4, ph9_p5, ph11_m6, ph11_m4, ph11_p5, ab_2, pc_135, pc_136 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h5_m4 = ph5_m4[k];
        const auto h5_p5 = ph5_p5[k];
        const auto h7_m6 = ph7_m6[k];
        const auto h7_m4 = ph7_m4[k];
        const auto h7_p5 = ph7_p5[k];
        const auto h9_m6 = ph9_m6[k];
        const auto h9_m4 = ph9_m4[k];
        const auto h9_p5 = ph9_p5[k];
        const auto h11_m6 = ph11_m6[k];
        const auto h11_m4 = ph11_m4[k];
        const auto h11_p5 = ph11_p5[k];

        pc_135[k] = - e_2 * fs_15_8_210 * h5_m4 - e_3 * fs_75_286_1001 * h7_m6 - e_3 * fs_525_572_154 * h7_m4 + e_3 * fs_45_52_210 * r_2 * h5_m4 - e_4 * fs_315_4862_715 * h9_m6 - e_4 * fs_63_9724_30030 * h9_m4 + e_4 * fs_150_2431_1001 * r_2 * h7_m6 + e_4 * fs_525_2431_154 * r_2 * h7_m4 - e_4 * fs_3_26_210 * r_4 * h5_m4 - e_5 * fs_6_4199_2431 * h11_m6 - e_5 * fs_3_8398_6006 * h11_m4 + e_5 * fs_15_2431_715 * r_2 * h9_m6 + e_5 * fs_3_4862_30030 * r_2 * h9_m4 - e_5 * fs_150_46189_1001 * r_4 * h7_m6 - e_5 * fs_525_46189_154 * r_4 * h7_m4 + e_5 * fs_1_221_210 * r_6 * h5_m4;

        pc_136[k] = - e_2 * f_75_4 * h5_p5 - e_3 * fs_525_286_66 * h7_p5 + e_3 * f_225_26 * r_2 * h5_p5 - e_4 * fs_315_4862_1001 * h9_p5 + e_4 * fs_1050_2431_66 * r_2 * h7_p5 - e_4 * f_15_13 * r_4 * h5_p5 - e_5 * fs_6_4199_2002 * h11_p5 + e_5 * fs_15_2431_1001 * r_2 * h9_p5 - e_5 * fs_1050_46189_66 * r_4 * h7_p5 + e_5 * f_10_221 * r_6 * h5_p5;
    }

    // NOTE: the angular components are formed in 124 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_2, pe_3, pe_4, pe_5, ph5_p4, ph7_p4, ph7_p6, ph9_p4, ph9_p6, ph11_p4, ph11_p6, ab_2, pc_137 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_2 = pe_2[k];
        const auto e_3 = pe_3[k];
        const auto e_4 = pe_4[k];
        const auto e_5 = pe_5[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;
        const auto r_6 = r_2 * r_4;

        const auto h5_p4 = ph5_p4[k];
        const auto h7_p4 = ph7_p4[k];
        const auto h7_p6 = ph7_p6[k];
        const auto h9_p4 = ph9_p4[k];
        const auto h9_p6 = ph9_p6[k];
        const auto h11_p4 = ph11_p4[k];
        const auto h11_p6 = ph11_p6[k];

        pc_137[k] = e_2 * fs_15_8_210 * h5_p4 + e_3 * fs_525_572_154 * h7_p4 - e_3 * fs_75_286_1001 * h7_p6 - e_3 * fs_45_52_210 * r_2 * h5_p4 + e_4 * fs_63_9724_30030 * h9_p4 - e_4 * fs_315_4862_715 * h9_p6 - e_4 * fs_525_2431_154 * r_2 * h7_p4 + e_4 * fs_150_2431_1001 * r_2 * h7_p6 + e_4 * fs_3_26_210 * r_4 * h5_p4 + e_5 * fs_3_8398_6006 * h11_p4 - e_5 * fs_6_4199_2431 * h11_p6 - e_5 * fs_3_4862_30030 * r_2 * h9_p4 + e_5 * fs_15_2431_715 * r_2 * h9_p6 + e_5 * fs_525_46189_154 * r_4 * h7_p4 - e_5 * fs_150_46189_1001 * r_4 * h7_p6 - e_5 * fs_1_221_210 * r_6 * h5_p4;
    }

    // NOTE: the angular components are formed in 124 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph3_p3, ph5_p3, ph7_p3, ph7_p7, ph9_p3, ph9_p7, ph11_p3, ph11_p7, ab_2, pc_138 : simd::cache_line_size())
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

        const auto h3_p3 = ph3_p3[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h7_p3 = ph7_p3[k];
        const auto h7_p7 = ph7_p7[k];
        const auto h9_p3 = ph9_p3[k];
        const auto h9_p7 = ph9_p7[k];
        const auto h11_p3 = ph11_p3[k];
        const auto h11_p7 = ph11_p7[k];

        pc_138[k] = - e_1 * fs_315_32_6 * h3_p3 - e_2 * fs_25_4_42 * h5_p3 + e_2 * fs_35_4_6 * r_2 * h3_p3 - e_3 * fs_525_286_35 * h7_p3 - e_3 * fs_105_572_715 * h7_p7 + e_3 * fs_75_26_42 * r_2 * h5_p3 - e_3 * fs_105_44_6 * r_4 * h3_p3 - e_4 * fs_315_4862_154 * h9_p3 - e_4 * fs_315_4862_858 * h9_p7 + e_4 * fs_1050_2431_35 * r_2 * h7_p3 + e_4 * fs_105_2431_715 * r_2 * h7_p7 - e_4 * fs_5_13_42 * r_4 * h5_p3 + e_4 * fs_35_143_6 * r_6 * h3_p3 - e_5 * fs_3_8398_2002 * h11_p3 - e_5 * fs_9_4199_2431 * h11_p7 + e_5 * fs_15_2431_154 * r_2 * h9_p3 + e_5 * fs_15_2431_858 * r_2 * h9_p7 - e_5 * fs_1050_46189_35 * r_4 * h7_p3 - e_5 * fs_105_46189_715 * r_4 * h7_p7 + e_5 * fs_10_663_42 * r_6 * h5_p3 - e_5 * fs_7_858_6 * r_8 * h3_p3;
    }

    // NOTE: the angular components are formed in 124 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_1, pe_2, pe_3, pe_4, pe_5, ph3_p2, ph5_p2, ph7_p2, ph9_p2, ph9_p8, ph11_p2, ph11_p8, ab_2, pc_139 : simd::cache_line_size())
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

        const auto h3_p2 = ph3_p2[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h9_p2 = ph9_p2[k];
        const auto h9_p8 = ph9_p8[k];
        const auto h11_p2 = ph11_p2[k];
        const auto h11_p8 = ph11_p8[k];

        pc_139[k] = e_1 * f_945_16 * h3_p2 + e_2 * fs_75_4_7 * h5_p2 - e_2 * f_105_2 * r_2 * h3_p2 + e_3 * fs_315_286_70 * h7_p2 - e_3 * fs_225_26_7 * r_2 * h5_p2 + e_3 * f_315_22 * r_4 * h3_p2 + e_4 * fs_315_4862_66 * h9_p2 - e_4 * fs_105_4862_7293 * h9_p8 - e_4 * fs_630_2431_70 * r_2 * h7_p2 + e_4 * fs_15_13_7 * r_4 * h5_p2 - e_4 * f_210_143 * r_6 * h3_p2 + e_5 * fs_3_4199_143 * h11_p2 - e_5 * fs_3_4199_46189 * h11_p8 - e_5 * fs_15_2431_66 * r_2 * h9_p2 + e_5 * fs_5_2431_7293 * r_2 * h9_p8 + e_5 * fs_630_46189_70 * r_4 * h7_p2 - e_5 * fs_10_221_7 * r_6 * h5_p2 + e_5 * f_7_143 * r_8 * h3_p2;
    }

    // NOTE: the angular components are formed in 124 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph1_p1, ph3_p1, ph5_p1, ph7_p1, ph9_p1, ph9_p9, ph11_p1, ph11_p9, ab_2, pc_140 : simd::cache_line_size())
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

        const auto h1_p1 = ph1_p1[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h9_p1 = ph9_p1[k];
        const auto h9_p9 = ph9_p9[k];
        const auto h11_p1 = ph11_p1[k];
        const auto h11_p9 = ph11_p9[k];

        pc_140[k] = - e_0 * fs_945_64_2 * h1_p1 - e_1 * fs_945_16_3 * h3_p1 + e_1 * fs_945_32_2 * r_2 * h1_p1 - e_2 * fs_75_8_30 * h5_p1 + e_2 * fs_105_2_3 * r_2 * h3_p1 - e_2 * fs_135_8_2 * r_4 * h1_p1 - e_3 * fs_525_286_14 * h7_p1 + e_3 * fs_225_52_30 * r_2 * h5_p1 - e_3 * fs_315_22_3 * r_4 * h3_p1 + e_3 * fs_15_4_2 * r_6 * h1_p1 - e_4 * fs_945_9724_10 * h9_p1 - e_4 * fs_63_4862_12155 * h9_p9 + e_4 * fs_1050_2431_14 * r_2 * h7_p1 - e_4 * fs_15_26_30 * r_4 * h5_p1 + e_4 * fs_210_143_3 * r_6 * h3_p1 - e_4 * fs_15_44_2 * r_8 * h1_p1 - e_5 * fs_3_4199_33 * h11_p1 - e_5 * fs_3_4199_92378 * h11_p9 + e_5 * fs_45_4862_10 * r_2 * h9_p1 + e_5 * fs_3_2431_12155 * r_2 * h9_p9 - e_5 * fs_1050_46189_14 * r_4 * h7_p1 + e_5 * fs_5_221_30 * r_6 * h5_p1 - e_5 * fs_7_143_3 * r_8 * h3_p1 + e_5 * fs_3_286_2 * r_10 * h1_p1;
    }

    // NOTE: the angular components are formed in 124 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph1_0, ph3_0, ph5_0, ph7_0, ph9_0, ph11_0, ph11_p10, ab_2, pc_141 : simd::cache_line_size())
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

        const auto h1_0 = ph1_0[k];
        const auto h3_0 = ph3_0[k];
        const auto h5_0 = ph5_0[k];
        const auto h7_0 = ph7_0[k];
        const auto h9_0 = ph9_0[k];
        const auto h11_0 = ph11_0[k];
        const auto h11_p10 = ph11_p10[k];

        pc_141[k] = e_0 * fs_945_32_11 * h1_0 + e_1 * fs_945_16_11 * h3_0 - e_1 * fs_945_16_11 * r_2 * h1_0 + e_2 * fs_75_4_11 * h5_0 - e_2 * fs_105_2_11 * r_2 * h3_0 + e_2 * fs_135_4_11 * r_4 * h1_0 + e_3 * fs_525_286_11 * h7_0 - e_3 * fs_225_26_11 * r_2 * h5_0 + e_3 * fs_315_22_11 * r_4 * h3_0 - e_3 * fs_15_2_11 * r_6 * h1_0 + e_4 * fs_315_4862_11 * h9_0 - e_4 * fs_1050_2431_11 * r_2 * h7_0 + e_4 * fs_15_13_11 * r_4 * h5_0 - e_4 * fs_210_143_11 * r_6 * h3_0 + e_4 * fs_15_22_11 * r_8 * h1_0 + e_5 * fs_3_4199_11 * h11_0 - e_5 * fs_3_4199_176358 * h11_p10 - e_5 * fs_15_2431_11 * r_2 * h9_0 + e_5 * fs_1050_46189_11 * r_4 * h7_0 - e_5 * fs_10_221_11 * r_6 * h5_0 + e_5 * fs_7_143_11 * r_8 * h3_0 - e_5 * fs_3_143_11 * r_10 * h1_0;
    }

    // NOTE: the angular components are formed in 124 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, pe_3, pe_4, pe_5, ph1_p1, ph3_p1, ph5_p1, ph7_p1, ph9_p1, ph11_p1, ph11_p11, ab_2, pc_142 : simd::cache_line_size())
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

        const auto h1_p1 = ph1_p1[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h9_p1 = ph9_p1[k];
        const auto h11_p1 = ph11_p1[k];
        const auto h11_p11 = ph11_p11[k];

        pc_142[k] = e_0 * fs_945_32_33 * h1_p1 + e_1 * fs_945_32_22 * h3_p1 - e_1 * fs_945_16_33 * r_2 * h1_p1 + e_2 * fs_15_4_55 * h5_p1 - e_2 * fs_105_4_22 * r_2 * h3_p1 + e_2 * fs_135_4_33 * r_4 * h1_p1 + e_3 * fs_75_572_231 * h7_p1 - e_3 * fs_45_26_55 * r_2 * h5_p1 + e_3 * fs_315_44_22 * r_4 * h3_p1 - e_3 * fs_15_2_33 * r_6 * h1_p1 + e_4 * fs_21_4862_165 * h9_p1 - e_4 * fs_75_2431_231 * r_2 * h7_p1 + e_4 * fs_3_13_55 * r_4 * h5_p1 - e_4 * fs_105_143_22 * r_6 * h3_p1 + e_4 * fs_15_22_33 * r_8 * h1_p1 + e_5 * fs_3_8398_2 * h11_p1 - e_5 * fs_3_4199_323323 * h11_p11 - e_5 * fs_1_2431_165 * r_2 * h9_p1 + e_5 * fs_75_46189_231 * r_4 * h7_p1 - e_5 * fs_2_221_55 * r_6 * h5_p1 + e_5 * fs_7_286_22 * r_8 * h3_p1 - e_5 * fs_3_143_33 * r_10 * h1_p1;
    }

    // NOTE: the values of a combination of angular components are stored as one
    // row of nvalues columns, with the component on bra side running slowest. The
    // rows which the symmetry relates to an already formed one are copied from it.

    const size_t sources[143] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142};

    for (size_t n = 0; n < 143; n++)
    {
        if (sources[n] != n) std::copy(values + sources[n] * nvalues, values + (sources[n] + 1) * nvalues, values + n * nvalues);
    }
}

}  // namespace simdt2ceri
