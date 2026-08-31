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



#include "SimdTwoCenterElectronRepulsionRecDI.hpp"

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
compute_di_electron_repulsion(double                         *values,
                              const size_t                    nvalues,
                              const CBasisFunction           &bra,
                              const CBasisFunction           &ket,
                              const std::vector<CSimdMatrix> &harmonics,
                              const CSimdMatrix              &coordinates) -> void
{
    if ((bra.get_angular_momentum() != 2) || (ket.get_angular_momentum() != 6))
    {
        errors::assertMsgCritical(
            false, std::string("SimdTwoCenterElectronRepulsionFunc.compute_di_electron_repulsion: Basis functions must be of angular momenta two and six"));
    }

    if (harmonics.size() < 8)
    {
        errors::assertMsgCritical(
            false, std::string("SimdTwoCenterElectronRepulsionFunc.compute_di_electron_repulsion: Harmonics must reach angular momentum 8"));
    }

    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("SimdTwoCenterElectronRepulsionFunc.compute_di_electron_repulsion: Number of values exceeds number of atom pairs"));
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
    // orders 6 to 8 alone, and the orders below them are formed on the
    // way to them by the recursion the Boys function is evaluated with.

    simdfunc::compute_boys_function(boys);

    // NOTE: the buffer holds the contracted prefactor of each exponent factor of
    // the terms, as the integrals of the angular components are formed straight
    // into the values and are not written a second time. Every exponent factor is
    // used with one order of the Boys function alone, so the order follows from
    // the factor and one accumulator per factor suffices.

    auto buffer = CSimdMatrix(3, nvalues);

    auto *pe_0 = buffer.data(0);
    auto *pe_1 = buffer.data(1);
    auto *pe_2 = buffer.data(2);

    std::fill(pe_0, pe_0 + nvalues, 0.0);
    std::fill(pe_1, pe_1 + nvalues, 0.0);
    std::fill(pe_2, pe_2 + nvalues, 0.0);

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

            const auto ff_0 = fbase * aexp * aexp * aexp * aexp / (fexp * fexp * fexp * fexp * fexp * fexp);

            const auto ff_1 = fbase * aexp * aexp * aexp * aexp * fmu / (fexp * fexp * fexp * fexp * fexp * fexp);

            const auto ff_2 = fbase * aexp * aexp * aexp * aexp * fmu * fmu / (fexp * fexp * fexp * fexp * fexp * fexp);

            const auto *bv_0 = boys.data(7, i * nprim_b + j);

            const auto *bv_1 = boys.data(8, i * nprim_b + j);

            const auto *bv_2 = boys.data(9, i * nprim_b + j);

#pragma omp simd aligned(pe_0, pe_1, pe_2, bv_0, bv_1, bv_2 : simd::cache_line_size())
            for (size_t k = 0; k < nvalues; k++)
            {
                pe_0[k] += ff_0 * bv_0[k];
                pe_1[k] += ff_1 * bv_1[k];
                pe_2[k] += ff_2 * bv_2[k];
            }
        }
    }

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

    // NOTE: the factors of the terms depend on the angular momenta alone, so they
    // are formed once for the whole matrix instead of once for every atom pair.

    const auto f_12_11 = 12.0 / 11.0;
    const auto f_12_143 = 12.0 / 143.0;
    const auto f_13_55 = 13.0 / 55.0;
    const auto f_14_55 = 14.0 / 55.0;
    const auto f_15_11 = 15.0 / 11.0;
    const auto f_15_22 = 15.0 / 22.0;
    const auto f_1_11 = 1.0 / 11.0;
    const auto f_1_5 = 1.0 / 5.0;
    const auto f_21_11 = 21.0 / 11.0;
    const auto f_28_65 = 28.0 / 65.0;
    const auto f_2_11 = 2.0 / 11.0;
    const auto f_2_13 = 2.0 / 13.0;
    const auto f_2_5 = 2.0 / 5.0;
    const auto f_2_55 = 2.0 / 55.0;
    const auto f_39_22 = 39.0 / 22.0;
    const auto f_3_1 = 3.0;
    const auto f_3_11 = 3.0 / 11.0;
    const auto f_3_2 = 3.0 / 2.0;
    const auto f_45_11 = 45.0 / 11.0;
    const auto f_45_143 = 45.0 / 143.0;
    const auto f_45_4 = 45.0 / 4.0;
    const auto fs_12_11_7 = std::sqrt(1008.0 / 121.0);
    const auto fs_12_143_7 = std::sqrt(1008.0 / 20449.0);
    const auto fs_14_65_3 = std::sqrt(588.0 / 4225.0);
    const auto fs_15_11_7 = std::sqrt(1575.0 / 121.0);
    const auto fs_15_143_7 = std::sqrt(1575.0 / 20449.0);
    const auto fs_15_22_3 = std::sqrt(675.0 / 484.0);
    const auto fs_15_4_7 = std::sqrt(1575.0 / 16.0);
    const auto fs_18_11_3 = std::sqrt(972.0 / 121.0);
    const auto fs_18_143_3 = std::sqrt(972.0 / 20449.0);
    const auto fs_1_11_3 = std::sqrt(3.0 / 121.0);
    const auto fs_1_130_10 = std::sqrt(1.0 / 1690.0);
    const auto fs_1_130_1430 = std::sqrt(11.0 / 130.0);
    const auto fs_1_130_2 = std::sqrt(1.0 / 8450.0);
    const auto fs_1_130_2002 = std::sqrt(77.0 / 650.0);
    const auto fs_1_130_26 = std::sqrt(1.0 / 650.0);
    const auto fs_1_130_2730 = std::sqrt(21.0 / 130.0);
    const auto fs_1_130_30 = std::sqrt(3.0 / 1690.0);
    const auto fs_1_130_70 = std::sqrt(7.0 / 1690.0);
    const auto fs_1_130_910 = std::sqrt(7.0 / 130.0);
    const auto fs_1_13_22 = std::sqrt(22.0 / 169.0);
    const auto fs_1_26_66 = std::sqrt(33.0 / 338.0);
    const auto fs_1_55_22 = std::sqrt(2.0 / 275.0);
    const auto fs_1_55_30 = std::sqrt(6.0 / 605.0);
    const auto fs_1_55_55 = std::sqrt(1.0 / 55.0);
    const auto fs_1_55_7 = std::sqrt(7.0 / 3025.0);
    const auto fs_1_65_165 = std::sqrt(33.0 / 845.0);
    const auto fs_1_65_210 = std::sqrt(42.0 / 845.0);
    const auto fs_1_65_429 = std::sqrt(33.0 / 325.0);
    const auto fs_1_65_55 = std::sqrt(11.0 / 845.0);
    const auto fs_1_65_70 = std::sqrt(14.0 / 845.0);
    const auto fs_1_65_91 = std::sqrt(7.0 / 325.0);
    const auto fs_1_65_910 = std::sqrt(14.0 / 65.0);
    const auto fs_21_22_3 = std::sqrt(1323.0 / 484.0);
    const auto fs_21_44_10 = std::sqrt(2205.0 / 968.0);
    const auto fs_2_55_30 = std::sqrt(24.0 / 605.0);
    const auto fs_2_55_70 = std::sqrt(56.0 / 605.0);
    const auto fs_2_65_110 = std::sqrt(88.0 / 845.0);
    const auto fs_2_65_6 = std::sqrt(24.0 / 4225.0);
    const auto fs_2_65_91 = std::sqrt(28.0 / 325.0);
    const auto fs_3_110_66 = std::sqrt(27.0 / 550.0);
    const auto fs_3_11_105 = std::sqrt(945.0 / 121.0);
    const auto fs_3_11_15 = std::sqrt(135.0 / 121.0);
    const auto fs_3_11_165 = std::sqrt(135.0 / 11.0);
    const auto fs_3_11_210 = std::sqrt(1890.0 / 121.0);
    const auto fs_3_11_30 = std::sqrt(270.0 / 121.0);
    const auto fs_3_11_35 = std::sqrt(315.0 / 121.0);
    const auto fs_3_11_70 = std::sqrt(630.0 / 121.0);
    const auto fs_3_130_110 = std::sqrt(99.0 / 1690.0);
    const auto fs_3_130_70 = std::sqrt(63.0 / 1690.0);
    const auto fs_3_13_3 = std::sqrt(27.0 / 169.0);
    const auto fs_3_143_105 = std::sqrt(945.0 / 20449.0);
    const auto fs_3_143_15 = std::sqrt(135.0 / 20449.0);
    const auto fs_3_143_165 = std::sqrt(135.0 / 1859.0);
    const auto fs_3_143_210 = std::sqrt(1890.0 / 20449.0);
    const auto fs_3_143_35 = std::sqrt(315.0 / 20449.0);
    const auto fs_3_143_70 = std::sqrt(630.0 / 20449.0);
    const auto fs_3_1_7 = std::sqrt(63.0);
    const auto fs_3_22_10 = std::sqrt(45.0 / 242.0);
    const auto fs_3_22_2 = std::sqrt(9.0 / 242.0);
    const auto fs_3_22_22 = std::sqrt(9.0 / 22.0);
    const auto fs_3_22_30 = std::sqrt(135.0 / 242.0);
    const auto fs_3_22_330 = std::sqrt(135.0 / 22.0);
    const auto fs_3_22_55 = std::sqrt(45.0 / 44.0);
    const auto fs_3_22_7 = std::sqrt(63.0 / 484.0);
    const auto fs_3_22_70 = std::sqrt(315.0 / 242.0);
    const auto fs_3_286_10 = std::sqrt(45.0 / 40898.0);
    const auto fs_3_286_2 = std::sqrt(9.0 / 40898.0);
    const auto fs_3_286_330 = std::sqrt(135.0 / 3718.0);
    const auto fs_3_286_70 = std::sqrt(315.0 / 40898.0);
    const auto fs_3_2_30 = std::sqrt(135.0 / 2.0);
    const auto fs_3_2_42 = std::sqrt(189.0 / 2.0);
    const auto fs_3_4_105 = std::sqrt(945.0 / 16.0);
    const auto fs_3_4_15 = std::sqrt(135.0 / 16.0);
    const auto fs_3_4_165 = std::sqrt(1485.0 / 16.0);
    const auto fs_3_4_210 = std::sqrt(945.0 / 8.0);
    const auto fs_3_4_35 = std::sqrt(315.0 / 16.0);
    const auto fs_3_4_70 = std::sqrt(315.0 / 8.0);
    const auto fs_3_55_10 = std::sqrt(18.0 / 605.0);
    const auto fs_3_65_26 = std::sqrt(18.0 / 325.0);
    const auto fs_3_65_7 = std::sqrt(63.0 / 4225.0);
    const auto fs_3_8_10 = std::sqrt(45.0 / 32.0);
    const auto fs_3_8_2 = std::sqrt(9.0 / 32.0);
    const auto fs_3_8_330 = std::sqrt(1485.0 / 32.0);
    const auto fs_3_8_70 = std::sqrt(315.0 / 32.0);
    const auto fs_6_11_30 = std::sqrt(1080.0 / 121.0);
    const auto fs_6_11_42 = std::sqrt(1512.0 / 121.0);
    const auto fs_6_143_30 = std::sqrt(1080.0 / 20449.0);
    const auto fs_6_143_42 = std::sqrt(1512.0 / 20449.0);
    const auto fs_6_65_10 = std::sqrt(72.0 / 845.0);
    const auto fs_6_65_11 = std::sqrt(396.0 / 4225.0);
    const auto fs_6_65_21 = std::sqrt(756.0 / 4225.0);
    const auto fs_7_110_10 = std::sqrt(49.0 / 1210.0);
    const auto fs_7_55_3 = std::sqrt(147.0 / 3025.0);
    const auto fs_8_65_7 = std::sqrt(448.0 / 4225.0);
    const auto fs_9_11_14 = std::sqrt(1134.0 / 121.0);
    const auto fs_9_11_5 = std::sqrt(405.0 / 121.0);
    const auto fs_9_11_7 = std::sqrt(567.0 / 121.0);
    const auto fs_9_143_14 = std::sqrt(1134.0 / 20449.0);
    const auto fs_9_143_5 = std::sqrt(405.0 / 20449.0);
    const auto fs_9_143_7 = std::sqrt(567.0 / 20449.0);
    const auto fs_9_22_10 = std::sqrt(405.0 / 242.0);
    const auto fs_9_22_110 = std::sqrt(405.0 / 22.0);
    const auto fs_9_22_2 = std::sqrt(81.0 / 242.0);
    const auto fs_9_286_110 = std::sqrt(405.0 / 3718.0);
    const auto fs_9_286_2 = std::sqrt(81.0 / 40898.0);
    const auto fs_9_2_3 = std::sqrt(243.0 / 4.0);
    const auto fs_9_44_66 = std::sqrt(243.0 / 88.0);
    const auto fs_9_4_14 = std::sqrt(567.0 / 8.0);
    const auto fs_9_4_5 = std::sqrt(405.0 / 16.0);
    const auto fs_9_4_7 = std::sqrt(567.0 / 16.0);
    const auto fs_9_8_110 = std::sqrt(4455.0 / 32.0);
    const auto fs_9_8_2 = std::sqrt(81.0 / 32.0);

    // NOTE: the angular components are formed in 23 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph4_p2, ph4_p3, ph4_p4, ph6_p2, ph6_p3, ph6_p4, ph6_p6, ph8_p2, ph8_p3, ph8_p4, ph8_p6, ph8_p7, ph8_p8, ab_2, pc_0, pc_1, pc_2 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

        const auto h4_p2 = ph4_p2[k];
        const auto h4_p3 = ph4_p3[k];
        const auto h4_p4 = ph4_p4[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h6_p4 = ph6_p4[k];
        const auto h6_p6 = ph6_p6[k];
        const auto h8_p2 = ph8_p2[k];
        const auto h8_p3 = ph8_p3[k];
        const auto h8_p4 = ph8_p4[k];
        const auto h8_p6 = ph8_p6[k];
        const auto h8_p7 = ph8_p7[k];
        const auto h8_p8 = ph8_p8[k];

        pc_0[k] = e_0 * fs_9_8_110 * h4_p4 + e_1 * fs_3_22_22 * h6_p4 - e_1 * fs_9_22_110 * r_2 * h4_p4 + e_2 * fs_1_130_2 * h8_p4 - e_2 * fs_1_65_910 * h8_p8 - e_2 * fs_1_55_22 * r_2 * h6_p4 + e_2 * fs_9_286_110 * r_4 * h4_p4;

        pc_1[k] = e_0 * fs_3_4_165 * h4_p3 + e_1 * fs_3_22_55 * h6_p3 - e_1 * fs_3_11_165 * r_2 * h4_p3 + e_2 * fs_1_130_10 * h8_p3 - e_2 * fs_1_130_2730 * h8_p7 - e_2 * fs_1_55_55 * r_2 * h6_p3 + e_2 * fs_3_143_165 * r_4 * h4_p3;

        pc_2[k] = e_0 * fs_3_4_105 * h4_p2 + e_1 * fs_9_22_10 * h6_p2 - e_1 * fs_3_22_22 * h6_p6 - e_1 * fs_3_11_105 * r_2 * h4_p2 + e_2 * fs_1_130_30 * h8_p2 - e_2 * fs_1_130_2002 * h8_p6 - e_2 * fs_3_55_10 * r_2 * h6_p2 + e_2 * fs_1_55_22 * r_2 * h6_p6 + e_2 * fs_3_143_105 * r_4 * h4_p2;
    }

    // NOTE: the angular components are formed in 23 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph4_0, ph4_p1, ph4_p4, ph6_0, ph6_p1, ph6_p4, ph6_p5, ph8_0, ph8_p1, ph8_p4, ph8_p5, ab_2, pc_3, pc_4 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

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

        pc_3[k] = e_0 * fs_9_4_7 * h4_p1 + e_1 * fs_3_11_30 * h6_p1 - e_1 * fs_3_22_55 * h6_p5 - e_1 * fs_9_11_7 * r_2 * h4_p1 + e_2 * fs_1_130_70 * h8_p1 - e_2 * fs_1_130_1430 * h8_p5 - e_2 * fs_2_55_30 * r_2 * h6_p1 + e_2 * fs_1_55_55 * r_2 * h6_p5 + e_2 * fs_9_143_7 * r_4 * h4_p1;

        pc_4[k] = e_0 * fs_3_4_70 * h4_0 - e_0 * fs_3_8_2 * h4_p4 + e_1 * fs_3_11_70 * h6_0 - e_1 * fs_9_22_10 * h6_p4 - e_1 * fs_3_11_70 * r_2 * h4_0 + e_1 * fs_3_22_2 * r_2 * h4_p4 + e_2 * fs_1_65_70 * h8_0 - e_2 * fs_3_130_110 * h8_p4 - e_2 * fs_2_55_70 * r_2 * h6_0 + e_2 * fs_3_55_10 * r_2 * h6_p4 + e_2 * fs_3_143_70 * r_4 * h4_0 - e_2 * fs_3_286_2 * r_4 * h4_p4;
    }

    // NOTE: the angular components are formed in 23 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph4_m2, ph4_p1, ph4_p3, ph6_m2, ph6_p1, ph6_p3, ph8_m2, ph8_p1, ph8_p3, ab_2, pc_5, pc_6 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

        const auto h4_m2 = ph4_m2[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h4_p3 = ph4_p3[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h8_m2 = ph8_m2[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h8_p3 = ph8_p3[k];

        pc_5[k] = - e_0 * fs_3_8_70 * h4_p1 - e_0 * fs_3_8_10 * h4_p3 - e_1 * fs_21_22_3 * h6_p1 - e_1 * fs_3_11_30 * h6_p3 + e_1 * fs_3_22_70 * r_2 * h4_p1 + e_1 * fs_3_22_10 * r_2 * h4_p3 - e_2 * fs_3_65_7 * h8_p1 - e_2 * fs_1_65_165 * h8_p3 + e_2 * fs_7_55_3 * r_2 * h6_p1 + e_2 * fs_2_55_30 * r_2 * h6_p3 - e_2 * fs_3_286_70 * r_4 * h4_p1 - e_2 * fs_3_286_10 * r_4 * h4_p3;

        pc_6[k] = e_0 * fs_3_4_15 * h4_m2 + e_1 * fs_3_11_70 * h6_m2 - e_1 * fs_3_11_15 * r_2 * h4_m2 + e_2 * fs_1_65_210 * h8_m2 - e_2 * fs_2_55_70 * r_2 * h6_m2 + e_2 * fs_3_143_15 * r_4 * h4_m2;
    }

    // NOTE: the angular components are formed in 23 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph4_m4, ph4_m3, ph4_m1, ph6_m5, ph6_m4, ph6_m3, ph6_m1, ph8_m5, ph8_m4, ph8_m3, ph8_m1, ab_2, pc_7, pc_8, pc_9 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

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

        pc_7[k] = e_0 * fs_3_8_10 * h4_m3 - e_0 * fs_3_8_70 * h4_m1 + e_1 * fs_3_11_30 * h6_m3 - e_1 * fs_21_22_3 * h6_m1 - e_1 * fs_3_22_10 * r_2 * h4_m3 + e_1 * fs_3_22_70 * r_2 * h4_m1 + e_2 * fs_1_65_165 * h8_m3 - e_2 * fs_3_65_7 * h8_m1 - e_2 * fs_2_55_30 * r_2 * h6_m3 + e_2 * fs_7_55_3 * r_2 * h6_m1 + e_2 * fs_3_286_10 * r_4 * h4_m3 - e_2 * fs_3_286_70 * r_4 * h4_m1;

        pc_8[k] = e_0 * fs_3_8_2 * h4_m4 + e_1 * fs_9_22_10 * h6_m4 - e_1 * fs_3_22_2 * r_2 * h4_m4 + e_2 * fs_3_130_110 * h8_m4 - e_2 * fs_3_55_10 * r_2 * h6_m4 + e_2 * fs_3_286_2 * r_4 * h4_m4;

        pc_9[k] = - e_0 * fs_9_4_7 * h4_m1 + e_1 * fs_3_22_55 * h6_m5 - e_1 * fs_3_11_30 * h6_m1 + e_1 * fs_9_11_7 * r_2 * h4_m1 + e_2 * fs_1_130_1430 * h8_m5 - e_2 * fs_1_130_70 * h8_m1 - e_2 * fs_1_55_55 * r_2 * h6_m5 + e_2 * fs_2_55_30 * r_2 * h6_m1 - e_2 * fs_9_143_7 * r_4 * h4_m1;
    }

    // NOTE: the angular components are formed in 23 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph4_m4, ph4_m3, ph4_m2, ph6_m6, ph6_m4, ph6_m3, ph6_m2, ph8_m8, ph8_m7, ph8_m6, ph8_m4, ph8_m3, ph8_m2, ab_2, pc_10, pc_11, pc_12 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

        const auto h4_m4 = ph4_m4[k];
        const auto h4_m3 = ph4_m3[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h6_m6 = ph6_m6[k];
        const auto h6_m4 = ph6_m4[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h8_m8 = ph8_m8[k];
        const auto h8_m7 = ph8_m7[k];
        const auto h8_m6 = ph8_m6[k];
        const auto h8_m4 = ph8_m4[k];
        const auto h8_m3 = ph8_m3[k];
        const auto h8_m2 = ph8_m2[k];

        pc_10[k] = - e_0 * fs_3_4_105 * h4_m2 + e_1 * fs_3_22_22 * h6_m6 - e_1 * fs_9_22_10 * h6_m2 + e_1 * fs_3_11_105 * r_2 * h4_m2 + e_2 * fs_1_130_2002 * h8_m6 - e_2 * fs_1_130_30 * h8_m2 - e_2 * fs_1_55_22 * r_2 * h6_m6 + e_2 * fs_3_55_10 * r_2 * h6_m2 - e_2 * fs_3_143_105 * r_4 * h4_m2;

        pc_11[k] = - e_0 * fs_3_4_165 * h4_m3 - e_1 * fs_3_22_55 * h6_m3 + e_1 * fs_3_11_165 * r_2 * h4_m3 + e_2 * fs_1_130_2730 * h8_m7 - e_2 * fs_1_130_10 * h8_m3 + e_2 * fs_1_55_55 * r_2 * h6_m3 - e_2 * fs_3_143_165 * r_4 * h4_m3;

        pc_12[k] = - e_0 * fs_9_8_110 * h4_m4 - e_1 * fs_3_22_22 * h6_m4 + e_1 * fs_9_22_110 * r_2 * h4_m4 + e_2 * fs_1_65_910 * h8_m8 - e_2 * fs_1_130_2 * h8_m4 + e_2 * fs_1_55_22 * r_2 * h6_m4 - e_2 * fs_9_286_110 * r_4 * h4_m4;
    }

    // NOTE: the angular components are formed in 23 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph4_p3, ph4_p4, ph6_p3, ph6_p4, ph6_p5, ph6_p6, ph8_p3, ph8_p4, ph8_p5, ph8_p6, ph8_p7, ab_2, pc_13, pc_14, pc_15 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

        const auto h4_p3 = ph4_p3[k];
        const auto h4_p4 = ph4_p4[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h6_p4 = ph6_p4[k];
        const auto h6_p5 = ph6_p5[k];
        const auto h6_p6 = ph6_p6[k];
        const auto h8_p3 = ph8_p3[k];
        const auto h8_p4 = ph8_p4[k];
        const auto h8_p5 = ph8_p5[k];
        const auto h8_p6 = ph8_p6[k];
        const auto h8_p7 = ph8_p7[k];

        pc_13[k] = - e_1 * f_3_2 * h6_p5 - e_2 * fs_1_130_26 * h8_p5 - e_2 * fs_1_130_910 * h8_p7 + e_2 * f_1_5 * r_2 * h6_p5;

        pc_14[k] = e_0 * fs_3_8_330 * h4_p4 - e_1 * fs_9_44_66 * h6_p4 + e_1 * f_3_2 * h6_p6 - e_1 * fs_3_22_330 * r_2 * h4_p4 - e_2 * fs_2_65_6 * h8_p4 - e_2 * fs_2_65_91 * h8_p6 + e_2 * fs_3_110_66 * r_2 * h6_p4 - e_2 * f_1_5 * r_2 * h6_p6 + e_2 * fs_3_286_330 * r_4 * h4_p4;

        pc_15[k] = e_0 * fs_3_2_30 * h4_p3 - e_1 * fs_21_44_10 * h6_p3 + e_1 * fs_9_44_66 * h6_p5 - e_1 * fs_6_11_30 * r_2 * h4_p3 - e_2 * fs_1_65_55 * h8_p3 - e_2 * fs_1_65_429 * h8_p5 + e_2 * fs_7_110_10 * r_2 * h6_p3 - e_2 * fs_3_110_66 * r_2 * h6_p5 + e_2 * fs_6_143_30 * r_4 * h4_p3;
    }

    // NOTE: the angular components are formed in 23 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph4_p1, ph4_p2, ph4_p3, ph4_p4, ph6_p1, ph6_p2, ph6_p3, ph6_p4, ph8_p1, ph8_p2, ph8_p3, ph8_p4, ab_2, pc_16, pc_17 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

        const auto h4_p1 = ph4_p1[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h4_p3 = ph4_p3[k];
        const auto h4_p4 = ph4_p4[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h6_p4 = ph6_p4[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h8_p2 = ph8_p2[k];
        const auto h8_p3 = ph8_p3[k];
        const auto h8_p4 = ph8_p4[k];

        pc_16[k] = e_0 * fs_9_4_14 * h4_p2 + e_0 * fs_9_8_2 * h4_p4 - e_1 * fs_15_22_3 * h6_p2 + e_1 * fs_21_44_10 * h6_p4 - e_1 * fs_9_11_14 * r_2 * h4_p2 - e_1 * fs_9_22_2 * r_2 * h4_p4 - e_2 * f_2_13 * h8_p2 - e_2 * fs_2_65_110 * h8_p4 + e_2 * fs_1_11_3 * r_2 * h6_p2 - e_2 * fs_7_110_10 * r_2 * h6_p4 + e_2 * fs_9_143_14 * r_4 * h4_p2 + e_2 * fs_9_286_2 * r_4 * h4_p4;

        pc_17[k] = e_0 * fs_3_1_7 * h4_p1 + e_0 * f_3_1 * h4_p3 - e_1 * fs_3_22_30 * h6_p1 + e_1 * fs_15_22_3 * h6_p3 - e_1 * fs_12_11_7 * r_2 * h4_p1 - e_1 * f_12_11 * r_2 * h4_p3 - e_2 * fs_3_130_70 * h8_p1 - e_2 * fs_1_26_66 * h8_p3 + e_2 * fs_1_55_30 * r_2 * h6_p1 - e_2 * fs_1_11_3 * r_2 * h6_p3 + e_2 * fs_12_143_7 * r_4 * h4_p1 + e_2 * f_12_143 * r_4 * h4_p3;
    }

    // NOTE: the angular components are formed in 23 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph4_m2, ph4_m1, ph4_0, ph4_p2, ph6_m2, ph6_m1, ph6_0, ph6_p2, ph8_m2, ph8_m1, ph8_0, ph8_p2, ab_2, pc_18, pc_19, pc_20 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

        const auto h4_m2 = ph4_m2[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h4_0 = ph4_0[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h6_0 = ph6_0[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h8_m2 = ph8_m2[k];
        const auto h8_m1 = ph8_m1[k];
        const auto h8_0 = ph8_0[k];
        const auto h8_p2 = ph8_p2[k];

        pc_18[k] = e_0 * fs_15_4_7 * h4_0 + e_0 * fs_3_4_35 * h4_p2 - e_1 * fs_3_22_7 * h6_0 + e_1 * fs_3_22_30 * h6_p2 - e_1 * fs_15_11_7 * r_2 * h4_0 - e_1 * fs_3_11_35 * r_2 * h4_p2 - e_2 * fs_8_65_7 * h8_0 - e_2 * fs_6_65_10 * h8_p2 + e_2 * fs_1_55_7 * r_2 * h6_0 - e_2 * fs_1_55_30 * r_2 * h6_p2 + e_2 * fs_15_143_7 * r_4 * h4_0 + e_2 * fs_3_143_35 * r_4 * h4_p2;

        pc_19[k] = - e_0 * fs_3_2_30 * h4_m1 - e_1 * fs_3_22_7 * h6_m1 + e_1 * fs_6_11_30 * r_2 * h4_m1 + e_2 * fs_14_65_3 * h8_m1 + e_2 * fs_1_55_7 * r_2 * h6_m1 - e_2 * fs_6_143_30 * r_4 * h4_m1;

        pc_20[k] = - e_0 * fs_3_4_35 * h4_m2 - e_1 * fs_3_22_30 * h6_m2 + e_1 * fs_3_11_35 * r_2 * h4_m2 + e_2 * fs_6_65_10 * h8_m2 + e_2 * fs_1_55_30 * r_2 * h6_m2 - e_2 * fs_3_143_35 * r_4 * h4_m2;
    }

    // NOTE: the angular components are formed in 23 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph4_m4, ph4_m3, ph4_m2, ph4_m1, ph6_m4, ph6_m3, ph6_m2, ph6_m1, ph8_m4, ph8_m3, ph8_m2, ph8_m1, ab_2, pc_21, pc_22 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

        const auto h4_m4 = ph4_m4[k];
        const auto h4_m3 = ph4_m3[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h6_m4 = ph6_m4[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h8_m4 = ph8_m4[k];
        const auto h8_m3 = ph8_m3[k];
        const auto h8_m2 = ph8_m2[k];
        const auto h8_m1 = ph8_m1[k];

        pc_21[k] = - e_0 * f_3_1 * h4_m3 - e_0 * fs_3_1_7 * h4_m1 - e_1 * fs_15_22_3 * h6_m3 + e_1 * fs_3_22_30 * h6_m1 + e_1 * f_12_11 * r_2 * h4_m3 + e_1 * fs_12_11_7 * r_2 * h4_m1 + e_2 * fs_1_26_66 * h8_m3 + e_2 * fs_3_130_70 * h8_m1 + e_2 * fs_1_11_3 * r_2 * h6_m3 - e_2 * fs_1_55_30 * r_2 * h6_m1 - e_2 * f_12_143 * r_4 * h4_m3 - e_2 * fs_12_143_7 * r_4 * h4_m1;

        pc_22[k] = - e_0 * fs_9_8_2 * h4_m4 - e_0 * fs_9_4_14 * h4_m2 - e_1 * fs_21_44_10 * h6_m4 + e_1 * fs_15_22_3 * h6_m2 + e_1 * fs_9_22_2 * r_2 * h4_m4 + e_1 * fs_9_11_14 * r_2 * h4_m2 + e_2 * fs_2_65_110 * h8_m4 + e_2 * f_2_13 * h8_m2 + e_2 * fs_7_110_10 * r_2 * h6_m4 - e_2 * fs_1_11_3 * r_2 * h6_m2 - e_2 * fs_9_286_2 * r_4 * h4_m4 - e_2 * fs_9_143_14 * r_4 * h4_m2;
    }

    // NOTE: the angular components are formed in 23 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph4_m4, ph4_m3, ph6_m6, ph6_m5, ph6_m4, ph6_m3, ph8_m7, ph8_m6, ph8_m5, ph8_m4, ph8_m3, ab_2, pc_23, pc_24, pc_25, pc_26, pc_27 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

        const auto h4_m4 = ph4_m4[k];
        const auto h4_m3 = ph4_m3[k];
        const auto h6_m6 = ph6_m6[k];
        const auto h6_m5 = ph6_m5[k];
        const auto h6_m4 = ph6_m4[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h8_m7 = ph8_m7[k];
        const auto h8_m6 = ph8_m6[k];
        const auto h8_m5 = ph8_m5[k];
        const auto h8_m4 = ph8_m4[k];
        const auto h8_m3 = ph8_m3[k];

        pc_23[k] = - e_0 * fs_3_2_30 * h4_m3 - e_1 * fs_9_44_66 * h6_m5 + e_1 * fs_21_44_10 * h6_m3 + e_1 * fs_6_11_30 * r_2 * h4_m3 + e_2 * fs_1_65_429 * h8_m5 + e_2 * fs_1_65_55 * h8_m3 + e_2 * fs_3_110_66 * r_2 * h6_m5 - e_2 * fs_7_110_10 * r_2 * h6_m3 - e_2 * fs_6_143_30 * r_4 * h4_m3;

        pc_24[k] = - e_0 * fs_3_8_330 * h4_m4 - e_1 * f_3_2 * h6_m6 + e_1 * fs_9_44_66 * h6_m4 + e_1 * fs_3_22_330 * r_2 * h4_m4 + e_2 * fs_2_65_91 * h8_m6 + e_2 * fs_2_65_6 * h8_m4 + e_2 * f_1_5 * r_2 * h6_m6 - e_2 * fs_3_110_66 * r_2 * h6_m4 - e_2 * fs_3_286_330 * r_4 * h4_m4;

        pc_25[k] = e_1 * f_3_2 * h6_m5 + e_2 * fs_1_130_910 * h8_m7 + e_2 * fs_1_130_26 * h8_m5 - e_2 * f_1_5 * r_2 * h6_m5;

        pc_26[k] = e_1 * f_3_1 * h6_m6 + e_2 * fs_1_65_91 * h8_m6 - e_2 * f_2_5 * r_2 * h6_m6;

        pc_27[k] = e_1 * f_3_2 * h6_m5 + e_2 * fs_3_65_26 * h8_m5 - e_2 * f_1_5 * r_2 * h6_m5;
    }

    // NOTE: the angular components are formed in 23 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph4_m4, ph4_m3, ph4_m2, ph4_m1, ph6_m4, ph6_m3, ph6_m2, ph6_m1, ph8_m4, ph8_m3, ph8_m2, ph8_m1, ab_2, pc_28, pc_29, pc_30, pc_31 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

        const auto h4_m4 = ph4_m4[k];
        const auto h4_m3 = ph4_m3[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h6_m4 = ph6_m4[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h8_m4 = ph8_m4[k];
        const auto h8_m3 = ph8_m3[k];
        const auto h8_m2 = ph8_m2[k];
        const auto h8_m1 = ph8_m1[k];

        pc_28[k] = e_0 * fs_9_4_5 * h4_m4 + e_1 * f_3_11 * h6_m4 - e_1 * fs_9_11_5 * r_2 * h4_m4 + e_2 * fs_6_65_11 * h8_m4 - e_2 * f_2_55 * r_2 * h6_m4 + e_2 * fs_9_143_5 * r_4 * h4_m4;

        pc_29[k] = e_0 * fs_9_2_3 * h4_m3 - e_1 * f_15_22 * h6_m3 - e_1 * fs_18_11_3 * r_2 * h4_m3 + e_2 * fs_1_13_22 * h8_m3 + e_2 * f_1_11 * r_2 * h6_m3 + e_2 * fs_18_143_3 * r_4 * h4_m3;

        pc_30[k] = e_0 * fs_3_2_42 * h4_m2 - e_1 * f_15_11 * h6_m2 - e_1 * fs_6_11_42 * r_2 * h4_m2 + e_2 * fs_3_13_3 * h8_m2 + e_2 * f_2_11 * r_2 * h6_m2 + e_2 * fs_6_143_42 * r_4 * h4_m2;

        pc_31[k] = e_0 * fs_3_4_210 * h4_m1 - e_1 * f_39_22 * h6_m1 - e_1 * fs_3_11_210 * r_2 * h4_m1 + e_2 * fs_6_65_21 * h8_m1 + e_2 * f_13_55 * r_2 * h6_m1 + e_2 * fs_3_143_210 * r_4 * h4_m1;
    }

    // NOTE: the angular components are formed in 23 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph4_0, ph4_p1, ph4_p2, ph4_p3, ph6_0, ph6_p1, ph6_p2, ph6_p3, ph8_0, ph8_p1, ph8_p2, ph8_p3, ab_2, pc_32, pc_33, pc_34, pc_35 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

        const auto h4_0 = ph4_0[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h4_p3 = ph4_p3[k];
        const auto h6_0 = ph6_0[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h8_0 = ph8_0[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h8_p2 = ph8_p2[k];
        const auto h8_p3 = ph8_p3[k];

        pc_32[k] = e_0 * f_45_4 * h4_0 - e_1 * f_21_11 * h6_0 - e_1 * f_45_11 * r_2 * h4_0 + e_2 * f_28_65 * h8_0 + e_2 * f_14_55 * r_2 * h6_0 + e_2 * f_45_143 * r_4 * h4_0;

        pc_33[k] = e_0 * fs_3_4_210 * h4_p1 - e_1 * f_39_22 * h6_p1 - e_1 * fs_3_11_210 * r_2 * h4_p1 + e_2 * fs_6_65_21 * h8_p1 + e_2 * f_13_55 * r_2 * h6_p1 + e_2 * fs_3_143_210 * r_4 * h4_p1;

        pc_34[k] = e_0 * fs_3_2_42 * h4_p2 - e_1 * f_15_11 * h6_p2 - e_1 * fs_6_11_42 * r_2 * h4_p2 + e_2 * fs_3_13_3 * h8_p2 + e_2 * f_2_11 * r_2 * h6_p2 + e_2 * fs_6_143_42 * r_4 * h4_p2;

        pc_35[k] = e_0 * fs_9_2_3 * h4_p3 - e_1 * f_15_22 * h6_p3 - e_1 * fs_18_11_3 * r_2 * h4_p3 + e_2 * fs_1_13_22 * h8_p3 + e_2 * f_1_11 * r_2 * h6_p3 + e_2 * fs_18_143_3 * r_4 * h4_p3;
    }

    // NOTE: the angular components are formed in 23 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph4_p4, ph6_m5, ph6_p4, ph6_p5, ph6_p6, ph8_m7, ph8_m5, ph8_p4, ph8_p5, ph8_p6, ab_2, pc_36, pc_37, pc_38, pc_39 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

        const auto h4_p4 = ph4_p4[k];
        const auto h6_m5 = ph6_m5[k];
        const auto h6_p4 = ph6_p4[k];
        const auto h6_p5 = ph6_p5[k];
        const auto h6_p6 = ph6_p6[k];
        const auto h8_m7 = ph8_m7[k];
        const auto h8_m5 = ph8_m5[k];
        const auto h8_p4 = ph8_p4[k];
        const auto h8_p5 = ph8_p5[k];
        const auto h8_p6 = ph8_p6[k];

        pc_36[k] = e_0 * fs_9_4_5 * h4_p4 + e_1 * f_3_11 * h6_p4 - e_1 * fs_9_11_5 * r_2 * h4_p4 + e_2 * fs_6_65_11 * h8_p4 - e_2 * f_2_55 * r_2 * h6_p4 + e_2 * fs_9_143_5 * r_4 * h4_p4;

        pc_37[k] = e_1 * f_3_2 * h6_p5 + e_2 * fs_3_65_26 * h8_p5 - e_2 * f_1_5 * r_2 * h6_p5;

        pc_38[k] = e_1 * f_3_1 * h6_p6 + e_2 * fs_1_65_91 * h8_p6 - e_2 * f_2_5 * r_2 * h6_p6;

        pc_39[k] = - e_1 * f_3_2 * h6_m5 + e_2 * fs_1_130_910 * h8_m7 - e_2 * fs_1_130_26 * h8_m5 + e_2 * f_1_5 * r_2 * h6_m5;
    }

    // NOTE: the angular components are formed in 23 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph4_m4, ph4_m3, ph4_m2, ph6_m6, ph6_m5, ph6_m4, ph6_m3, ph6_m2, ph8_m6, ph8_m5, ph8_m4, ph8_m3, ph8_m2, ab_2, pc_40, pc_41, pc_42 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

        const auto h4_m4 = ph4_m4[k];
        const auto h4_m3 = ph4_m3[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h6_m6 = ph6_m6[k];
        const auto h6_m5 = ph6_m5[k];
        const auto h6_m4 = ph6_m4[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h8_m6 = ph8_m6[k];
        const auto h8_m5 = ph8_m5[k];
        const auto h8_m4 = ph8_m4[k];
        const auto h8_m3 = ph8_m3[k];
        const auto h8_m2 = ph8_m2[k];

        pc_40[k] = e_0 * fs_3_8_330 * h4_m4 - e_1 * f_3_2 * h6_m6 - e_1 * fs_9_44_66 * h6_m4 - e_1 * fs_3_22_330 * r_2 * h4_m4 + e_2 * fs_2_65_91 * h8_m6 - e_2 * fs_2_65_6 * h8_m4 + e_2 * f_1_5 * r_2 * h6_m6 + e_2 * fs_3_110_66 * r_2 * h6_m4 + e_2 * fs_3_286_330 * r_4 * h4_m4;

        pc_41[k] = e_0 * fs_3_2_30 * h4_m3 - e_1 * fs_9_44_66 * h6_m5 - e_1 * fs_21_44_10 * h6_m3 - e_1 * fs_6_11_30 * r_2 * h4_m3 + e_2 * fs_1_65_429 * h8_m5 - e_2 * fs_1_65_55 * h8_m3 + e_2 * fs_3_110_66 * r_2 * h6_m5 + e_2 * fs_7_110_10 * r_2 * h6_m3 + e_2 * fs_6_143_30 * r_4 * h4_m3;

        pc_42[k] = - e_0 * fs_9_8_2 * h4_m4 + e_0 * fs_9_4_14 * h4_m2 - e_1 * fs_21_44_10 * h6_m4 - e_1 * fs_15_22_3 * h6_m2 + e_1 * fs_9_22_2 * r_2 * h4_m4 - e_1 * fs_9_11_14 * r_2 * h4_m2 + e_2 * fs_2_65_110 * h8_m4 - e_2 * f_2_13 * h8_m2 + e_2 * fs_7_110_10 * r_2 * h6_m4 + e_2 * fs_1_11_3 * r_2 * h6_m2 - e_2 * fs_9_286_2 * r_4 * h4_m4 + e_2 * fs_9_143_14 * r_4 * h4_m2;
    }

    // NOTE: the angular components are formed in 23 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph4_m3, ph4_m2, ph4_m1, ph4_p1, ph6_m3, ph6_m2, ph6_m1, ph6_p1, ph8_m3, ph8_m2, ph8_m1, ph8_p1, ab_2, pc_43, pc_44, pc_45 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

        const auto h4_m3 = ph4_m3[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h8_m3 = ph8_m3[k];
        const auto h8_m2 = ph8_m2[k];
        const auto h8_m1 = ph8_m1[k];
        const auto h8_p1 = ph8_p1[k];

        pc_43[k] = - e_0 * f_3_1 * h4_m3 + e_0 * fs_3_1_7 * h4_m1 - e_1 * fs_15_22_3 * h6_m3 - e_1 * fs_3_22_30 * h6_m1 + e_1 * f_12_11 * r_2 * h4_m3 - e_1 * fs_12_11_7 * r_2 * h4_m1 + e_2 * fs_1_26_66 * h8_m3 - e_2 * fs_3_130_70 * h8_m1 + e_2 * fs_1_11_3 * r_2 * h6_m3 + e_2 * fs_1_55_30 * r_2 * h6_m1 - e_2 * f_12_143 * r_4 * h4_m3 + e_2 * fs_12_143_7 * r_4 * h4_m1;

        pc_44[k] = - e_0 * fs_3_4_35 * h4_m2 - e_1 * fs_3_22_30 * h6_m2 + e_1 * fs_3_11_35 * r_2 * h4_m2 + e_2 * fs_6_65_10 * h8_m2 + e_2 * fs_1_55_30 * r_2 * h6_m2 - e_2 * fs_3_143_35 * r_4 * h4_m2;

        pc_45[k] = - e_0 * fs_3_2_30 * h4_p1 - e_1 * fs_3_22_7 * h6_p1 + e_1 * fs_6_11_30 * r_2 * h4_p1 + e_2 * fs_14_65_3 * h8_p1 + e_2 * fs_1_55_7 * r_2 * h6_p1 - e_2 * fs_6_143_30 * r_4 * h4_p1;
    }

    // NOTE: the angular components are formed in 23 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph4_0, ph4_p1, ph4_p2, ph4_p3, ph6_0, ph6_p1, ph6_p2, ph6_p3, ph8_0, ph8_p1, ph8_p2, ph8_p3, ab_2, pc_46, pc_47 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

        const auto h4_0 = ph4_0[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h4_p3 = ph4_p3[k];
        const auto h6_0 = ph6_0[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h8_0 = ph8_0[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h8_p2 = ph8_p2[k];
        const auto h8_p3 = ph8_p3[k];

        pc_46[k] = e_0 * fs_15_4_7 * h4_0 - e_0 * fs_3_4_35 * h4_p2 - e_1 * fs_3_22_7 * h6_0 - e_1 * fs_3_22_30 * h6_p2 - e_1 * fs_15_11_7 * r_2 * h4_0 + e_1 * fs_3_11_35 * r_2 * h4_p2 - e_2 * fs_8_65_7 * h8_0 + e_2 * fs_6_65_10 * h8_p2 + e_2 * fs_1_55_7 * r_2 * h6_0 + e_2 * fs_1_55_30 * r_2 * h6_p2 + e_2 * fs_15_143_7 * r_4 * h4_0 - e_2 * fs_3_143_35 * r_4 * h4_p2;

        pc_47[k] = e_0 * fs_3_1_7 * h4_p1 - e_0 * f_3_1 * h4_p3 - e_1 * fs_3_22_30 * h6_p1 - e_1 * fs_15_22_3 * h6_p3 - e_1 * fs_12_11_7 * r_2 * h4_p1 + e_1 * f_12_11 * r_2 * h4_p3 - e_2 * fs_3_130_70 * h8_p1 + e_2 * fs_1_26_66 * h8_p3 + e_2 * fs_1_55_30 * r_2 * h6_p1 + e_2 * fs_1_11_3 * r_2 * h6_p3 + e_2 * fs_12_143_7 * r_4 * h4_p1 - e_2 * f_12_143 * r_4 * h4_p3;
    }

    // NOTE: the angular components are formed in 23 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph4_p2, ph4_p3, ph4_p4, ph6_p2, ph6_p3, ph6_p4, ph6_p5, ph6_p6, ph8_p2, ph8_p3, ph8_p4, ph8_p5, ph8_p6, ab_2, pc_48, pc_49, pc_50 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

        const auto h4_p2 = ph4_p2[k];
        const auto h4_p3 = ph4_p3[k];
        const auto h4_p4 = ph4_p4[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h6_p4 = ph6_p4[k];
        const auto h6_p5 = ph6_p5[k];
        const auto h6_p6 = ph6_p6[k];
        const auto h8_p2 = ph8_p2[k];
        const auto h8_p3 = ph8_p3[k];
        const auto h8_p4 = ph8_p4[k];
        const auto h8_p5 = ph8_p5[k];
        const auto h8_p6 = ph8_p6[k];

        pc_48[k] = e_0 * fs_9_4_14 * h4_p2 - e_0 * fs_9_8_2 * h4_p4 - e_1 * fs_15_22_3 * h6_p2 - e_1 * fs_21_44_10 * h6_p4 - e_1 * fs_9_11_14 * r_2 * h4_p2 + e_1 * fs_9_22_2 * r_2 * h4_p4 - e_2 * f_2_13 * h8_p2 + e_2 * fs_2_65_110 * h8_p4 + e_2 * fs_1_11_3 * r_2 * h6_p2 + e_2 * fs_7_110_10 * r_2 * h6_p4 + e_2 * fs_9_143_14 * r_4 * h4_p2 - e_2 * fs_9_286_2 * r_4 * h4_p4;

        pc_49[k] = e_0 * fs_3_2_30 * h4_p3 - e_1 * fs_21_44_10 * h6_p3 - e_1 * fs_9_44_66 * h6_p5 - e_1 * fs_6_11_30 * r_2 * h4_p3 - e_2 * fs_1_65_55 * h8_p3 + e_2 * fs_1_65_429 * h8_p5 + e_2 * fs_7_110_10 * r_2 * h6_p3 + e_2 * fs_3_110_66 * r_2 * h6_p5 + e_2 * fs_6_143_30 * r_4 * h4_p3;

        pc_50[k] = e_0 * fs_3_8_330 * h4_p4 - e_1 * fs_9_44_66 * h6_p4 - e_1 * f_3_2 * h6_p6 - e_1 * fs_3_22_330 * r_2 * h4_p4 - e_2 * fs_2_65_6 * h8_p4 + e_2 * fs_2_65_91 * h8_p6 + e_2 * fs_3_110_66 * r_2 * h6_p4 + e_2 * f_1_5 * r_2 * h6_p6 + e_2 * fs_3_286_330 * r_4 * h4_p4;
    }

    // NOTE: the angular components are formed in 23 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph4_m4, ph4_m3, ph6_m4, ph6_m3, ph6_p5, ph8_m8, ph8_m7, ph8_m4, ph8_m3, ph8_p5, ph8_p7, ab_2, pc_51, pc_52, pc_53 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

        const auto h4_m4 = ph4_m4[k];
        const auto h4_m3 = ph4_m3[k];
        const auto h6_m4 = ph6_m4[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h6_p5 = ph6_p5[k];
        const auto h8_m8 = ph8_m8[k];
        const auto h8_m7 = ph8_m7[k];
        const auto h8_m4 = ph8_m4[k];
        const auto h8_m3 = ph8_m3[k];
        const auto h8_p5 = ph8_p5[k];
        const auto h8_p7 = ph8_p7[k];

        pc_51[k] = - e_1 * f_3_2 * h6_p5 - e_2 * fs_1_130_26 * h8_p5 + e_2 * fs_1_130_910 * h8_p7 + e_2 * f_1_5 * r_2 * h6_p5;

        pc_52[k] = e_0 * fs_9_8_110 * h4_m4 + e_1 * fs_3_22_22 * h6_m4 - e_1 * fs_9_22_110 * r_2 * h4_m4 + e_2 * fs_1_65_910 * h8_m8 + e_2 * fs_1_130_2 * h8_m4 - e_2 * fs_1_55_22 * r_2 * h6_m4 + e_2 * fs_9_286_110 * r_4 * h4_m4;

        pc_53[k] = e_0 * fs_3_4_165 * h4_m3 + e_1 * fs_3_22_55 * h6_m3 - e_1 * fs_3_11_165 * r_2 * h4_m3 + e_2 * fs_1_130_2730 * h8_m7 + e_2 * fs_1_130_10 * h8_m3 - e_2 * fs_1_55_55 * r_2 * h6_m3 + e_2 * fs_3_143_165 * r_4 * h4_m3;
    }

    // NOTE: the angular components are formed in 23 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph4_m4, ph4_m2, ph4_m1, ph6_m6, ph6_m5, ph6_m4, ph6_m2, ph6_m1, ph8_m6, ph8_m5, ph8_m4, ph8_m2, ph8_m1, ab_2, pc_54, pc_55, pc_56 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

        const auto h4_m4 = ph4_m4[k];
        const auto h4_m2 = ph4_m2[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h6_m6 = ph6_m6[k];
        const auto h6_m5 = ph6_m5[k];
        const auto h6_m4 = ph6_m4[k];
        const auto h6_m2 = ph6_m2[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h8_m6 = ph8_m6[k];
        const auto h8_m5 = ph8_m5[k];
        const auto h8_m4 = ph8_m4[k];
        const auto h8_m2 = ph8_m2[k];
        const auto h8_m1 = ph8_m1[k];

        pc_54[k] = e_0 * fs_3_4_105 * h4_m2 + e_1 * fs_3_22_22 * h6_m6 + e_1 * fs_9_22_10 * h6_m2 - e_1 * fs_3_11_105 * r_2 * h4_m2 + e_2 * fs_1_130_2002 * h8_m6 + e_2 * fs_1_130_30 * h8_m2 - e_2 * fs_1_55_22 * r_2 * h6_m6 - e_2 * fs_3_55_10 * r_2 * h6_m2 + e_2 * fs_3_143_105 * r_4 * h4_m2;

        pc_55[k] = e_0 * fs_9_4_7 * h4_m1 + e_1 * fs_3_22_55 * h6_m5 + e_1 * fs_3_11_30 * h6_m1 - e_1 * fs_9_11_7 * r_2 * h4_m1 + e_2 * fs_1_130_1430 * h8_m5 + e_2 * fs_1_130_70 * h8_m1 - e_2 * fs_1_55_55 * r_2 * h6_m5 - e_2 * fs_2_55_30 * r_2 * h6_m1 + e_2 * fs_9_143_7 * r_4 * h4_m1;

        pc_56[k] = e_0 * fs_3_8_2 * h4_m4 + e_1 * fs_9_22_10 * h6_m4 - e_1 * fs_3_22_2 * r_2 * h4_m4 + e_2 * fs_3_130_110 * h8_m4 - e_2 * fs_3_55_10 * r_2 * h6_m4 + e_2 * fs_3_286_2 * r_4 * h4_m4;
    }

    // NOTE: the angular components are formed in 23 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph4_m3, ph4_m1, ph4_p2, ph6_m3, ph6_m1, ph6_p2, ph8_m3, ph8_m1, ph8_p2, ab_2, pc_57, pc_58 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

        const auto h4_m3 = ph4_m3[k];
        const auto h4_m1 = ph4_m1[k];
        const auto h4_p2 = ph4_p2[k];
        const auto h6_m3 = ph6_m3[k];
        const auto h6_m1 = ph6_m1[k];
        const auto h6_p2 = ph6_p2[k];
        const auto h8_m3 = ph8_m3[k];
        const auto h8_m1 = ph8_m1[k];
        const auto h8_p2 = ph8_p2[k];

        pc_57[k] = e_0 * fs_3_8_10 * h4_m3 + e_0 * fs_3_8_70 * h4_m1 + e_1 * fs_3_11_30 * h6_m3 + e_1 * fs_21_22_3 * h6_m1 - e_1 * fs_3_22_10 * r_2 * h4_m3 - e_1 * fs_3_22_70 * r_2 * h4_m1 + e_2 * fs_1_65_165 * h8_m3 + e_2 * fs_3_65_7 * h8_m1 - e_2 * fs_2_55_30 * r_2 * h6_m3 - e_2 * fs_7_55_3 * r_2 * h6_m1 + e_2 * fs_3_286_10 * r_4 * h4_m3 + e_2 * fs_3_286_70 * r_4 * h4_m1;

        pc_58[k] = e_0 * fs_3_4_15 * h4_p2 + e_1 * fs_3_11_70 * h6_p2 - e_1 * fs_3_11_15 * r_2 * h4_p2 + e_2 * fs_1_65_210 * h8_p2 - e_2 * fs_2_55_70 * r_2 * h6_p2 + e_2 * fs_3_143_15 * r_4 * h4_p2;
    }

    // NOTE: the angular components are formed in 23 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph4_0, ph4_p1, ph4_p3, ph4_p4, ph6_0, ph6_p1, ph6_p3, ph6_p4, ph8_0, ph8_p1, ph8_p3, ph8_p4, ab_2, pc_59, pc_60 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

        const auto h4_0 = ph4_0[k];
        const auto h4_p1 = ph4_p1[k];
        const auto h4_p3 = ph4_p3[k];
        const auto h4_p4 = ph4_p4[k];
        const auto h6_0 = ph6_0[k];
        const auto h6_p1 = ph6_p1[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h6_p4 = ph6_p4[k];
        const auto h8_0 = ph8_0[k];
        const auto h8_p1 = ph8_p1[k];
        const auto h8_p3 = ph8_p3[k];
        const auto h8_p4 = ph8_p4[k];

        pc_59[k] = - e_0 * fs_3_8_70 * h4_p1 + e_0 * fs_3_8_10 * h4_p3 - e_1 * fs_21_22_3 * h6_p1 + e_1 * fs_3_11_30 * h6_p3 + e_1 * fs_3_22_70 * r_2 * h4_p1 - e_1 * fs_3_22_10 * r_2 * h4_p3 - e_2 * fs_3_65_7 * h8_p1 + e_2 * fs_1_65_165 * h8_p3 + e_2 * fs_7_55_3 * r_2 * h6_p1 - e_2 * fs_2_55_30 * r_2 * h6_p3 - e_2 * fs_3_286_70 * r_4 * h4_p1 + e_2 * fs_3_286_10 * r_4 * h4_p3;

        pc_60[k] = e_0 * fs_3_4_70 * h4_0 + e_0 * fs_3_8_2 * h4_p4 + e_1 * fs_3_11_70 * h6_0 + e_1 * fs_9_22_10 * h6_p4 - e_1 * fs_3_11_70 * r_2 * h4_0 - e_1 * fs_3_22_2 * r_2 * h4_p4 + e_2 * fs_1_65_70 * h8_0 + e_2 * fs_3_130_110 * h8_p4 - e_2 * fs_2_55_70 * r_2 * h6_0 - e_2 * fs_3_55_10 * r_2 * h6_p4 + e_2 * fs_3_143_70 * r_4 * h4_0 + e_2 * fs_3_286_2 * r_4 * h4_p4;
    }

    // NOTE: the angular components are formed in 23 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph4_p1, ph4_p2, ph6_p1, ph6_p2, ph6_p5, ph6_p6, ph8_p1, ph8_p2, ph8_p5, ph8_p6, ab_2, pc_61, pc_62 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

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

        pc_61[k] = e_0 * fs_9_4_7 * h4_p1 + e_1 * fs_3_11_30 * h6_p1 + e_1 * fs_3_22_55 * h6_p5 - e_1 * fs_9_11_7 * r_2 * h4_p1 + e_2 * fs_1_130_70 * h8_p1 + e_2 * fs_1_130_1430 * h8_p5 - e_2 * fs_2_55_30 * r_2 * h6_p1 - e_2 * fs_1_55_55 * r_2 * h6_p5 + e_2 * fs_9_143_7 * r_4 * h4_p1;

        pc_62[k] = e_0 * fs_3_4_105 * h4_p2 + e_1 * fs_9_22_10 * h6_p2 + e_1 * fs_3_22_22 * h6_p6 - e_1 * fs_3_11_105 * r_2 * h4_p2 + e_2 * fs_1_130_30 * h8_p2 + e_2 * fs_1_130_2002 * h8_p6 - e_2 * fs_3_55_10 * r_2 * h6_p2 - e_2 * fs_1_55_22 * r_2 * h6_p6 + e_2 * fs_3_143_105 * r_4 * h4_p2;
    }

    // NOTE: the angular components are formed in 23 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph4_p3, ph4_p4, ph6_p3, ph6_p4, ph8_p3, ph8_p4, ph8_p7, ph8_p8, ab_2, pc_63, pc_64 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

        const auto h4_p3 = ph4_p3[k];
        const auto h4_p4 = ph4_p4[k];
        const auto h6_p3 = ph6_p3[k];
        const auto h6_p4 = ph6_p4[k];
        const auto h8_p3 = ph8_p3[k];
        const auto h8_p4 = ph8_p4[k];
        const auto h8_p7 = ph8_p7[k];
        const auto h8_p8 = ph8_p8[k];

        pc_63[k] = e_0 * fs_3_4_165 * h4_p3 + e_1 * fs_3_22_55 * h6_p3 - e_1 * fs_3_11_165 * r_2 * h4_p3 + e_2 * fs_1_130_10 * h8_p3 + e_2 * fs_1_130_2730 * h8_p7 - e_2 * fs_1_55_55 * r_2 * h6_p3 + e_2 * fs_3_143_165 * r_4 * h4_p3;

        pc_64[k] = e_0 * fs_9_8_110 * h4_p4 + e_1 * fs_3_22_22 * h6_p4 - e_1 * fs_9_22_110 * r_2 * h4_p4 + e_2 * fs_1_130_2 * h8_p4 + e_2 * fs_1_65_910 * h8_p8 - e_2 * fs_1_55_22 * r_2 * h6_p4 + e_2 * fs_9_286_110 * r_4 * h4_p4;
    }

    // NOTE: the values of a combination of angular components are stored as one
    // row of nvalues columns, with the component on bra side running slowest. The
    // rows which the symmetry relates to an already formed one are copied from it.

    const size_t sources[65] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64};

    for (size_t n = 0; n < 65; n++)
    {
        if (sources[n] != n) std::copy(values + sources[n] * nvalues, values + (sources[n] + 1) * nvalues, values + n * nvalues);
    }
}

}  // namespace simdt2ceri
