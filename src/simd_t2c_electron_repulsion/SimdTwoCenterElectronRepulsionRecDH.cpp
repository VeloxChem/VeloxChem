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



#include "SimdTwoCenterElectronRepulsionRecDH.hpp"

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
compute_dh_electron_repulsion(double                         *values,
                              const size_t                    nvalues,
                              const CBasisFunction           &bra,
                              const CBasisFunction           &ket,
                              const std::vector<CSimdMatrix> &harmonics,
                              const CSimdMatrix              &coordinates) -> void
{
    if ((bra.get_angular_momentum() != 2) || (ket.get_angular_momentum() != 5))
    {
        errors::assertMsgCritical(
            false, std::string("SimdTwoCenterElectronRepulsionFunc.compute_dh_electron_repulsion: Basis functions must be of angular momenta two and five"));
    }

    if (harmonics.size() < 7)
    {
        errors::assertMsgCritical(
            false, std::string("SimdTwoCenterElectronRepulsionFunc.compute_dh_electron_repulsion: Harmonics must reach angular momentum 7"));
    }

    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("SimdTwoCenterElectronRepulsionFunc.compute_dh_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    auto boys = CSimdVariableMatrix(std::vector<size_t>(nprims, nvalues), 9);

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
    // call, which fills the orders 0 to 7 of every row. The terms read the
    // orders 5 to 7 alone, and the orders below them are formed on the
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

            const auto ff_0 = fbase * aexp * aexp * aexp / (fexp * fexp * fexp * fexp * fexp);

            const auto ff_1 = fbase * aexp * aexp * aexp * fmu / (fexp * fexp * fexp * fexp * fexp);

            const auto ff_2 = fbase * aexp * aexp * aexp * fmu * fmu / (fexp * fexp * fexp * fexp * fexp);

            const auto *bv_0 = boys.data(6, i * nprim_b + j);

            const auto *bv_1 = boys.data(7, i * nprim_b + j);

            const auto *bv_2 = boys.data(8, i * nprim_b + j);

#pragma omp simd aligned(pe_0, pe_1, pe_2, bv_0, bv_1, bv_2 : simd::cache_line_size())
            for (size_t k = 0; k < nvalues; k++)
            {
                pe_0[k] += ff_0 * bv_0[k];
                pe_1[k] += ff_1 * bv_1[k];
                pe_2[k] += ff_2 * bv_2[k];
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

    // NOTE: the factors of the terms depend on the angular momenta alone, so they
    // are formed once for the whole matrix instead of once for every atom pair.

    const auto f_10_3 = 10.0 / 3.0;
    const auto f_10_33 = 10.0 / 33.0;
    const auto f_10_39 = 10.0 / 39.0;
    const auto f_15_2 = 15.0 / 2.0;
    const auto f_1_39 = 1.0 / 39.0;
    const auto f_1_6 = 1.0 / 6.0;
    const auto f_2_13 = 2.0 / 13.0;
    const auto f_3_13 = 3.0 / 13.0;
    const auto f_3_2 = 3.0 / 2.0;
    const auto f_5_13 = 5.0 / 13.0;
    const auto f_5_2 = 5.0 / 2.0;
    const auto f_5_3 = 5.0 / 3.0;
    const auto f_63_143 = 63.0 / 143.0;
    const auto fs_12_143_15 = std::sqrt(2160.0 / 20449.0);
    const auto fs_12_143_21 = std::sqrt(3024.0 / 20449.0);
    const auto fs_12_143_5 = std::sqrt(720.0 / 20449.0);
    const auto fs_15_4_2 = std::sqrt(225.0 / 8.0);
    const auto fs_18_143_10 = std::sqrt(3240.0 / 20449.0);
    const auto fs_1_11_10 = std::sqrt(10.0 / 121.0);
    const auto fs_1_11_7 = std::sqrt(7.0 / 121.0);
    const auto fs_1_1_10 = std::sqrt(10.0);
    const auto fs_1_1_7 = std::sqrt(7.0);
    const auto fs_1_26_30 = std::sqrt(15.0 / 338.0);
    const auto fs_1_33_105 = std::sqrt(35.0 / 363.0);
    const auto fs_1_33_35 = std::sqrt(35.0 / 1089.0);
    const auto fs_1_33_42 = std::sqrt(14.0 / 363.0);
    const auto fs_1_33_5 = std::sqrt(5.0 / 1089.0);
    const auto fs_1_39_15 = std::sqrt(5.0 / 507.0);
    const auto fs_1_39_21 = std::sqrt(7.0 / 507.0);
    const auto fs_1_39_5 = std::sqrt(5.0 / 1521.0);
    const auto fs_1_3_105 = std::sqrt(35.0 / 3.0);
    const auto fs_1_3_14 = std::sqrt(14.0 / 9.0);
    const auto fs_1_3_35 = std::sqrt(35.0 / 9.0);
    const auto fs_1_3_42 = std::sqrt(14.0 / 3.0);
    const auto fs_1_3_5 = std::sqrt(5.0 / 9.0);
    const auto fs_1_4_30 = std::sqrt(15.0 / 8.0);
    const auto fs_1_66_14 = std::sqrt(7.0 / 2178.0);
    const auto fs_1_66_2 = std::sqrt(1.0 / 2178.0);
    const auto fs_1_66_210 = std::sqrt(35.0 / 726.0);
    const auto fs_1_66_30 = std::sqrt(5.0 / 726.0);
    const auto fs_1_6_14 = std::sqrt(7.0 / 18.0);
    const auto fs_1_6_15 = std::sqrt(5.0 / 12.0);
    const auto fs_1_6_2 = std::sqrt(1.0 / 18.0);
    const auto fs_1_6_21 = std::sqrt(7.0 / 12.0);
    const auto fs_1_6_210 = std::sqrt(35.0 / 6.0);
    const auto fs_1_6_30 = std::sqrt(5.0 / 6.0);
    const auto fs_1_6_5 = std::sqrt(5.0 / 36.0);
    const auto fs_21_143_5 = std::sqrt(2205.0 / 20449.0);
    const auto fs_2_33_14 = std::sqrt(56.0 / 1089.0);
    const auto fs_2_33_3 = std::sqrt(4.0 / 363.0);
    const auto fs_2_33_7 = std::sqrt(28.0 / 1089.0);
    const auto fs_2_39_14 = std::sqrt(56.0 / 1521.0);
    const auto fs_2_39_35 = std::sqrt(140.0 / 1521.0);
    const auto fs_2_3_14 = std::sqrt(56.0 / 9.0);
    const auto fs_2_3_3 = std::sqrt(4.0 / 3.0);
    const auto fs_2_3_7 = std::sqrt(28.0 / 9.0);
    const auto fs_3_143_105 = std::sqrt(945.0 / 20449.0);
    const auto fs_3_143_143 = std::sqrt(9.0 / 143.0);
    const auto fs_3_143_165 = std::sqrt(135.0 / 1859.0);
    const auto fs_3_143_210 = std::sqrt(1890.0 / 20449.0);
    const auto fs_3_143_35 = std::sqrt(315.0 / 20449.0);
    const auto fs_3_143_66 = std::sqrt(54.0 / 1859.0);
    const auto fs_3_1_5 = std::sqrt(45.0);
    const auto fs_3_286_10 = std::sqrt(45.0 / 40898.0);
    const auto fs_3_286_1430 = std::sqrt(45.0 / 286.0);
    const auto fs_3_286_2 = std::sqrt(9.0 / 40898.0);
    const auto fs_3_286_2002 = std::sqrt(63.0 / 286.0);
    const auto fs_3_286_22 = std::sqrt(9.0 / 3718.0);
    const auto fs_3_286_30 = std::sqrt(135.0 / 40898.0);
    const auto fs_3_2_14 = std::sqrt(63.0 / 2.0);
    const auto fs_3_2_3 = std::sqrt(27.0 / 4.0);
    const auto fs_3_2_7 = std::sqrt(63.0 / 4.0);
    const auto fs_3_4_105 = std::sqrt(945.0 / 16.0);
    const auto fs_3_4_35 = std::sqrt(315.0 / 16.0);
    const auto fs_3_4_42 = std::sqrt(189.0 / 8.0);
    const auto fs_3_4_5 = std::sqrt(45.0 / 16.0);
    const auto fs_3_8_14 = std::sqrt(63.0 / 32.0);
    const auto fs_3_8_2 = std::sqrt(9.0 / 32.0);
    const auto fs_3_8_210 = std::sqrt(945.0 / 32.0);
    const auto fs_3_8_30 = std::sqrt(135.0 / 32.0);
    const auto fs_4_33_5 = std::sqrt(80.0 / 1089.0);
    const auto fs_4_3_5 = std::sqrt(80.0 / 9.0);
    const auto fs_5_33_2 = std::sqrt(50.0 / 1089.0);
    const auto fs_5_39_2 = std::sqrt(50.0 / 1521.0);
    const auto fs_5_39_3 = std::sqrt(25.0 / 507.0);
    const auto fs_5_3_2 = std::sqrt(50.0 / 9.0);
    const auto fs_5_6_2 = std::sqrt(25.0 / 18.0);
    const auto fs_5_6_3 = std::sqrt(25.0 / 12.0);
    const auto fs_6_143_105 = std::sqrt(3780.0 / 20449.0);
    const auto fs_6_143_5 = std::sqrt(180.0 / 20449.0);
    const auto fs_6_143_55 = std::sqrt(180.0 / 1859.0);
    const auto fs_7_12_6 = std::sqrt(49.0 / 24.0);
    const auto fs_7_78_6 = std::sqrt(49.0 / 1014.0);
    const auto fs_9_143_14 = std::sqrt(1134.0 / 20449.0);
    const auto fs_9_143_30 = std::sqrt(2430.0 / 20449.0);
    const auto fs_9_143_5 = std::sqrt(405.0 / 20449.0);
    const auto fs_9_286_110 = std::sqrt(405.0 / 3718.0);
    const auto fs_9_4_10 = std::sqrt(405.0 / 8.0);
    const auto fs_9_4_7 = std::sqrt(567.0 / 16.0);

    // NOTE: the angular components are formed in 18 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph3_p1, ph3_p2, ph3_p3, ph5_p1, ph5_p2, ph5_p3, ph5_p5, ph7_p1, ph7_p2, ph7_p3, ph7_p5, ph7_p6, ph7_p7, ab_2, pc_0, pc_1, pc_2 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

        const auto h3_p1 = ph3_p1[k];
        const auto h3_p2 = ph3_p2[k];
        const auto h3_p3 = ph3_p3[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h5_p5 = ph5_p5[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h7_p3 = ph7_p3[k];
        const auto h7_p5 = ph7_p5[k];
        const auto h7_p6 = ph7_p6[k];
        const auto h7_p7 = ph7_p7[k];

        pc_0[k] = e_0 * fs_3_4_105 * h3_p3 + e_1 * fs_1_6_15 * h5_p3 - e_1 * fs_1_3_105 * r_2 * h3_p3 + e_2 * fs_3_286_2 * h7_p3 - e_2 * fs_3_286_2002 * h7_p7 - e_2 * fs_1_39_15 * r_2 * h5_p3 + e_2 * fs_1_33_105 * r_4 * h3_p3;

        pc_1[k] = e_0 * fs_9_4_7 * h3_p2 + e_1 * h5_p2 - e_1 * fs_1_1_7 * r_2 * h3_p2 + e_2 * fs_3_286_10 * h7_p2 - e_2 * fs_3_286_1430 * h7_p6 - e_2 * f_2_13 * r_2 * h5_p2 + e_2 * fs_1_11_7 * r_4 * h3_p2;

        pc_2[k] = e_0 * fs_3_4_35 * h3_p1 + e_1 * fs_1_3_14 * h5_p1 - e_1 * fs_1_6_15 * h5_p5 - e_1 * fs_1_3_35 * r_2 * h3_p1 + e_2 * fs_3_286_30 * h7_p1 - e_2 * fs_9_286_110 * h7_p5 - e_2 * fs_2_39_14 * r_2 * h5_p1 + e_2 * fs_1_39_15 * r_2 * h5_p5 + e_2 * fs_1_33_35 * r_4 * h3_p1;
    }

    // NOTE: the angular components are formed in 18 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph3_0, ph3_p1, ph3_p3, ph5_0, ph5_p1, ph5_p3, ph5_p4, ph7_0, ph7_p1, ph7_p3, ph7_p4, ab_2, pc_3, pc_4 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

        const auto h3_0 = ph3_0[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h3_p3 = ph3_p3[k];
        const auto h5_0 = ph5_0[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h5_p4 = ph5_p4[k];
        const auto h7_0 = ph7_0[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h7_p3 = ph7_p3[k];
        const auto h7_p4 = ph7_p4[k];

        pc_3[k] = e_0 * fs_3_4_35 * h3_0 + e_1 * fs_1_3_35 * h5_0 - e_1 * h5_p4 - e_1 * fs_1_3_35 * r_2 * h3_0 + e_2 * fs_3_143_35 * h7_0 - e_2 * fs_3_143_165 * h7_p4 - e_2 * fs_2_39_35 * r_2 * h5_0 + e_2 * f_2_13 * r_2 * h5_p4 + e_2 * fs_1_33_35 * r_4 * h3_0;

        pc_4[k] = - e_0 * fs_3_8_30 * h3_p1 - e_0 * fs_3_8_2 * h3_p3 - e_1 * fs_5_6_3 * h5_p1 - e_1 * fs_1_3_14 * h5_p3 + e_1 * fs_1_6_30 * r_2 * h3_p1 + e_1 * fs_1_6_2 * r_2 * h3_p3 - e_2 * fs_3_143_35 * h7_p1 - e_2 * fs_3_143_105 * h7_p3 + e_2 * fs_5_39_3 * r_2 * h5_p1 + e_2 * fs_2_39_14 * r_2 * h5_p3 - e_2 * fs_1_66_30 * r_4 * h3_p1 - e_2 * fs_1_66_2 * r_4 * h3_p3;
    }

    // NOTE: the angular components are formed in 18 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph3_m3, ph3_m2, ph3_m1, ph5_m4, ph5_m3, ph5_m2, ph5_m1, ph7_m4, ph7_m3, ph7_m2, ph7_m1, ab_2, pc_5, pc_6, pc_7 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

        const auto h3_m3 = ph3_m3[k];
        const auto h3_m2 = ph3_m2[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h5_m4 = ph5_m4[k];
        const auto h5_m3 = ph5_m3[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h7_m4 = ph7_m4[k];
        const auto h7_m3 = ph7_m3[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h7_m1 = ph7_m1[k];

        pc_5[k] = e_0 * fs_3_4_5 * h3_m2 + e_1 * fs_1_3_35 * h5_m2 - e_1 * fs_1_3_5 * r_2 * h3_m2 + e_2 * fs_9_143_14 * h7_m2 - e_2 * fs_2_39_35 * r_2 * h5_m2 + e_2 * fs_1_33_5 * r_4 * h3_m2;

        pc_6[k] = e_0 * fs_3_8_2 * h3_m3 - e_0 * fs_3_8_30 * h3_m1 + e_1 * fs_1_3_14 * h5_m3 - e_1 * fs_5_6_3 * h5_m1 - e_1 * fs_1_6_2 * r_2 * h3_m3 + e_1 * fs_1_6_30 * r_2 * h3_m1 + e_2 * fs_3_143_105 * h7_m3 - e_2 * fs_3_143_35 * h7_m1 - e_2 * fs_2_39_14 * r_2 * h5_m3 + e_2 * fs_5_39_3 * r_2 * h5_m1 + e_2 * fs_1_66_2 * r_4 * h3_m3 - e_2 * fs_1_66_30 * r_4 * h3_m1;

        pc_7[k] = e_1 * h5_m4 + e_2 * fs_3_143_165 * h7_m4 - e_2 * f_2_13 * r_2 * h5_m4;
    }

    // NOTE: the angular components are formed in 18 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph3_m3, ph3_m2, ph3_m1, ph5_m5, ph5_m3, ph5_m2, ph5_m1, ph7_m7, ph7_m6, ph7_m5, ph7_m3, ph7_m2, ph7_m1, ab_2, pc_8, pc_9, pc_10 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

        const auto h3_m3 = ph3_m3[k];
        const auto h3_m2 = ph3_m2[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h5_m5 = ph5_m5[k];
        const auto h5_m3 = ph5_m3[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h7_m7 = ph7_m7[k];
        const auto h7_m6 = ph7_m6[k];
        const auto h7_m5 = ph7_m5[k];
        const auto h7_m3 = ph7_m3[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h7_m1 = ph7_m1[k];

        pc_8[k] = - e_0 * fs_3_4_35 * h3_m1 + e_1 * fs_1_6_15 * h5_m5 - e_1 * fs_1_3_14 * h5_m1 + e_1 * fs_1_3_35 * r_2 * h3_m1 + e_2 * fs_9_286_110 * h7_m5 - e_2 * fs_3_286_30 * h7_m1 - e_2 * fs_1_39_15 * r_2 * h5_m5 + e_2 * fs_2_39_14 * r_2 * h5_m1 - e_2 * fs_1_33_35 * r_4 * h3_m1;

        pc_9[k] = - e_0 * fs_9_4_7 * h3_m2 - e_1 * h5_m2 + e_1 * fs_1_1_7 * r_2 * h3_m2 + e_2 * fs_3_286_1430 * h7_m6 - e_2 * fs_3_286_10 * h7_m2 + e_2 * f_2_13 * r_2 * h5_m2 - e_2 * fs_1_11_7 * r_4 * h3_m2;

        pc_10[k] = - e_0 * fs_3_4_105 * h3_m3 - e_1 * fs_1_6_15 * h5_m3 + e_1 * fs_1_3_105 * r_2 * h3_m3 + e_2 * fs_3_286_2002 * h7_m7 - e_2 * fs_3_286_2 * h7_m3 + e_2 * fs_1_39_15 * r_2 * h5_m3 - e_2 * fs_1_33_105 * r_4 * h3_m3;
    }

    // NOTE: the angular components are formed in 18 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph3_p2, ph3_p3, ph5_p2, ph5_p3, ph5_p4, ph5_p5, ph7_p2, ph7_p3, ph7_p4, ph7_p5, ph7_p6, ab_2, pc_11, pc_12, pc_13 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

        const auto h3_p2 = ph3_p2[k];
        const auto h3_p3 = ph3_p3[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h5_p4 = ph5_p4[k];
        const auto h5_p5 = ph5_p5[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h7_p3 = ph7_p3[k];
        const auto h7_p4 = ph7_p4[k];
        const auto h7_p5 = ph7_p5[k];
        const auto h7_p6 = ph7_p6[k];

        pc_11[k] = - e_1 * fs_1_4_30 * h5_p4 - e_2 * fs_3_286_22 * h7_p4 - e_2 * fs_3_143_143 * h7_p6 + e_2 * fs_1_26_30 * r_2 * h5_p4;

        pc_12[k] = e_0 * fs_3_4_42 * h3_p3 - e_1 * fs_7_12_6 * h5_p3 + e_1 * fs_1_4_30 * h5_p5 - e_1 * fs_1_3_42 * r_2 * h3_p3 - e_2 * fs_6_143_5 * h7_p3 - e_2 * fs_6_143_55 * h7_p5 + e_2 * fs_7_78_6 * r_2 * h5_p3 - e_2 * fs_1_26_30 * r_2 * h5_p5 + e_2 * fs_1_33_42 * r_4 * h3_p3;

        pc_13[k] = e_0 * fs_3_2_14 * h3_p2 - e_1 * fs_5_6_2 * h5_p2 + e_1 * fs_7_12_6 * h5_p4 - e_1 * fs_2_3_14 * r_2 * h3_p2 - e_2 * fs_9_143_5 * h7_p2 - e_2 * fs_9_286_110 * h7_p4 + e_2 * fs_5_39_2 * r_2 * h5_p2 - e_2 * fs_7_78_6 * r_2 * h5_p4 + e_2 * fs_2_33_14 * r_4 * h3_p2;
    }

    // NOTE: the angular components are formed in 18 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph3_0, ph3_p1, ph3_p2, ph3_p3, ph5_0, ph5_p1, ph5_p2, ph5_p3, ph7_0, ph7_p1, ph7_p2, ph7_p3, ab_2, pc_14, pc_15 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

        const auto h3_0 = ph3_0[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h3_p2 = ph3_p2[k];
        const auto h3_p3 = ph3_p3[k];
        const auto h5_0 = ph5_0[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h7_0 = ph7_0[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h7_p3 = ph7_p3[k];

        pc_14[k] = e_0 * fs_3_8_210 * h3_p1 + e_0 * fs_3_8_14 * h3_p3 - e_1 * fs_1_6_21 * h5_p1 + e_1 * fs_5_6_2 * h5_p3 - e_1 * fs_1_6_210 * r_2 * h3_p1 - e_1 * fs_1_6_14 * r_2 * h3_p3 - e_2 * fs_12_143_5 * h7_p1 - e_2 * fs_12_143_15 * h7_p3 + e_2 * fs_1_39_21 * r_2 * h5_p1 - e_2 * fs_5_39_2 * r_2 * h5_p3 + e_2 * fs_1_66_210 * r_4 * h3_p1 + e_2 * fs_1_66_14 * r_4 * h3_p3;

        pc_15[k] = e_0 * fs_3_1_5 * h3_0 + e_0 * fs_3_2_3 * h3_p2 - e_1 * fs_1_6_5 * h5_0 + e_1 * fs_1_6_21 * h5_p2 - e_1 * fs_4_3_5 * r_2 * h3_0 - e_1 * fs_2_3_3 * r_2 * h3_p2 - e_2 * fs_21_143_5 * h7_0 - e_2 * fs_3_143_210 * h7_p2 + e_2 * fs_1_39_5 * r_2 * h5_0 - e_2 * fs_1_39_21 * r_2 * h5_p2 + e_2 * fs_4_33_5 * r_4 * h3_0 + e_2 * fs_2_33_3 * r_4 * h3_p2;
    }

    // NOTE: the angular components are formed in 18 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph3_m3, ph3_m2, ph3_m1, ph5_m4, ph5_m3, ph5_m2, ph5_m1, ph7_m4, ph7_m3, ph7_m2, ph7_m1, ab_2, pc_16, pc_17, pc_18, pc_19 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

        const auto h3_m3 = ph3_m3[k];
        const auto h3_m2 = ph3_m2[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h5_m4 = ph5_m4[k];
        const auto h5_m3 = ph5_m3[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h7_m4 = ph7_m4[k];
        const auto h7_m3 = ph7_m3[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h7_m1 = ph7_m1[k];

        pc_16[k] = - e_0 * fs_15_4_2 * h3_m1 - e_1 * fs_1_6_5 * h5_m1 + e_1 * fs_5_3_2 * r_2 * h3_m1 + e_2 * fs_12_143_21 * h7_m1 + e_2 * fs_1_39_5 * r_2 * h5_m1 - e_2 * fs_5_33_2 * r_4 * h3_m1;

        pc_17[k] = - e_0 * fs_3_2_3 * h3_m2 - e_1 * fs_1_6_21 * h5_m2 + e_1 * fs_2_3_3 * r_2 * h3_m2 + e_2 * fs_3_143_210 * h7_m2 + e_2 * fs_1_39_21 * r_2 * h5_m2 - e_2 * fs_2_33_3 * r_4 * h3_m2;

        pc_18[k] = - e_0 * fs_3_8_14 * h3_m3 - e_0 * fs_3_8_210 * h3_m1 - e_1 * fs_5_6_2 * h5_m3 + e_1 * fs_1_6_21 * h5_m1 + e_1 * fs_1_6_14 * r_2 * h3_m3 + e_1 * fs_1_6_210 * r_2 * h3_m1 + e_2 * fs_12_143_15 * h7_m3 + e_2 * fs_12_143_5 * h7_m1 + e_2 * fs_5_39_2 * r_2 * h5_m3 - e_2 * fs_1_39_21 * r_2 * h5_m1 - e_2 * fs_1_66_14 * r_4 * h3_m3 - e_2 * fs_1_66_210 * r_4 * h3_m1;

        pc_19[k] = - e_0 * fs_3_2_14 * h3_m2 - e_1 * fs_7_12_6 * h5_m4 + e_1 * fs_5_6_2 * h5_m2 + e_1 * fs_2_3_14 * r_2 * h3_m2 + e_2 * fs_9_286_110 * h7_m4 + e_2 * fs_9_143_5 * h7_m2 + e_2 * fs_7_78_6 * r_2 * h5_m4 - e_2 * fs_5_39_2 * r_2 * h5_m2 - e_2 * fs_2_33_14 * r_4 * h3_m2;
    }

    // NOTE: the angular components are formed in 18 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph3_m3, ph5_m5, ph5_m4, ph5_m3, ph7_m6, ph7_m5, ph7_m4, ph7_m3, ab_2, pc_20, pc_21, pc_22, pc_23, pc_24 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

        const auto h3_m3 = ph3_m3[k];
        const auto h5_m5 = ph5_m5[k];
        const auto h5_m4 = ph5_m4[k];
        const auto h5_m3 = ph5_m3[k];
        const auto h7_m6 = ph7_m6[k];
        const auto h7_m5 = ph7_m5[k];
        const auto h7_m4 = ph7_m4[k];
        const auto h7_m3 = ph7_m3[k];

        pc_20[k] = - e_0 * fs_3_4_42 * h3_m3 - e_1 * fs_1_4_30 * h5_m5 + e_1 * fs_7_12_6 * h5_m3 + e_1 * fs_1_3_42 * r_2 * h3_m3 + e_2 * fs_6_143_55 * h7_m5 + e_2 * fs_6_143_5 * h7_m3 + e_2 * fs_1_26_30 * r_2 * h5_m5 - e_2 * fs_7_78_6 * r_2 * h5_m3 - e_2 * fs_1_33_42 * r_4 * h3_m3;

        pc_21[k] = e_1 * fs_1_4_30 * h5_m4 + e_2 * fs_3_143_143 * h7_m6 + e_2 * fs_3_286_22 * h7_m4 - e_2 * fs_1_26_30 * r_2 * h5_m4;

        pc_22[k] = e_1 * f_5_2 * h5_m5 + e_2 * fs_3_143_66 * h7_m5 - e_2 * f_5_13 * r_2 * h5_m5;

        pc_23[k] = e_1 * h5_m4 + e_2 * fs_3_143_165 * h7_m4 - e_2 * f_2_13 * r_2 * h5_m4;

        pc_24[k] = e_0 * fs_3_2_7 * h3_m3 - e_1 * f_1_6 * h5_m3 - e_1 * fs_2_3_7 * r_2 * h3_m3 + e_2 * fs_9_143_30 * h7_m3 + e_2 * f_1_39 * r_2 * h5_m3 + e_2 * fs_2_33_7 * r_4 * h3_m3;
    }

    // NOTE: the angular components are formed in 18 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph3_m2, ph3_m1, ph3_0, ph3_p1, ph5_m2, ph5_m1, ph5_0, ph5_p1, ph7_m2, ph7_m1, ph7_0, ph7_p1, ab_2, pc_25, pc_26, pc_27, pc_28 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

        const auto h3_m2 = ph3_m2[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h3_0 = ph3_0[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h5_0 = ph5_0[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h7_m1 = ph7_m1[k];
        const auto h7_0 = ph7_0[k];
        const auto h7_p1 = ph7_p1[k];

        pc_25[k] = e_0 * fs_9_4_7 * h3_m2 - e_1 * h5_m2 - e_1 * fs_1_1_7 * r_2 * h3_m2 + e_2 * fs_18_143_10 * h7_m2 + e_2 * f_2_13 * r_2 * h5_m2 + e_2 * fs_1_11_7 * r_4 * h3_m2;

        pc_26[k] = e_0 * fs_9_4_10 * h3_m1 - e_1 * f_3_2 * h5_m1 - e_1 * fs_1_1_10 * r_2 * h3_m1 + e_2 * fs_6_143_105 * h7_m1 + e_2 * f_3_13 * r_2 * h5_m1 + e_2 * fs_1_11_10 * r_4 * h3_m1;

        pc_27[k] = e_0 * f_15_2 * h3_0 - e_1 * f_5_3 * h5_0 - e_1 * f_10_3 * r_2 * h3_0 + e_2 * f_63_143 * h7_0 + e_2 * f_10_39 * r_2 * h5_0 + e_2 * f_10_33 * r_4 * h3_0;

        pc_28[k] = e_0 * fs_9_4_10 * h3_p1 - e_1 * f_3_2 * h5_p1 - e_1 * fs_1_1_10 * r_2 * h3_p1 + e_2 * fs_6_143_105 * h7_p1 + e_2 * f_3_13 * r_2 * h5_p1 + e_2 * fs_1_11_10 * r_4 * h3_p1;
    }

    // NOTE: the angular components are formed in 18 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph3_p2, ph3_p3, ph5_p2, ph5_p3, ph5_p4, ph5_p5, ph7_p2, ph7_p3, ph7_p4, ph7_p5, ab_2, pc_29, pc_30, pc_31, pc_32 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

        const auto h3_p2 = ph3_p2[k];
        const auto h3_p3 = ph3_p3[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h5_p4 = ph5_p4[k];
        const auto h5_p5 = ph5_p5[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h7_p3 = ph7_p3[k];
        const auto h7_p4 = ph7_p4[k];
        const auto h7_p5 = ph7_p5[k];

        pc_29[k] = e_0 * fs_9_4_7 * h3_p2 - e_1 * h5_p2 - e_1 * fs_1_1_7 * r_2 * h3_p2 + e_2 * fs_18_143_10 * h7_p2 + e_2 * f_2_13 * r_2 * h5_p2 + e_2 * fs_1_11_7 * r_4 * h3_p2;

        pc_30[k] = e_0 * fs_3_2_7 * h3_p3 - e_1 * f_1_6 * h5_p3 - e_1 * fs_2_3_7 * r_2 * h3_p3 + e_2 * fs_9_143_30 * h7_p3 + e_2 * f_1_39 * r_2 * h5_p3 + e_2 * fs_2_33_7 * r_4 * h3_p3;

        pc_31[k] = e_1 * h5_p4 + e_2 * fs_3_143_165 * h7_p4 - e_2 * f_2_13 * r_2 * h5_p4;

        pc_32[k] = e_1 * f_5_2 * h5_p5 + e_2 * fs_3_143_66 * h7_p5 - e_2 * f_5_13 * r_2 * h5_p5;
    }

    // NOTE: the angular components are formed in 18 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph3_m3, ph3_m2, ph5_m5, ph5_m4, ph5_m3, ph5_m2, ph7_m6, ph7_m5, ph7_m4, ph7_m3, ph7_m2, ab_2, pc_33, pc_34, pc_35 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

        const auto h3_m3 = ph3_m3[k];
        const auto h3_m2 = ph3_m2[k];
        const auto h5_m5 = ph5_m5[k];
        const auto h5_m4 = ph5_m4[k];
        const auto h5_m3 = ph5_m3[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h7_m6 = ph7_m6[k];
        const auto h7_m5 = ph7_m5[k];
        const auto h7_m4 = ph7_m4[k];
        const auto h7_m3 = ph7_m3[k];
        const auto h7_m2 = ph7_m2[k];

        pc_33[k] = - e_1 * fs_1_4_30 * h5_m4 + e_2 * fs_3_143_143 * h7_m6 - e_2 * fs_3_286_22 * h7_m4 + e_2 * fs_1_26_30 * r_2 * h5_m4;

        pc_34[k] = e_0 * fs_3_4_42 * h3_m3 - e_1 * fs_1_4_30 * h5_m5 - e_1 * fs_7_12_6 * h5_m3 - e_1 * fs_1_3_42 * r_2 * h3_m3 + e_2 * fs_6_143_55 * h7_m5 - e_2 * fs_6_143_5 * h7_m3 + e_2 * fs_1_26_30 * r_2 * h5_m5 + e_2 * fs_7_78_6 * r_2 * h5_m3 + e_2 * fs_1_33_42 * r_4 * h3_m3;

        pc_35[k] = e_0 * fs_3_2_14 * h3_m2 - e_1 * fs_7_12_6 * h5_m4 - e_1 * fs_5_6_2 * h5_m2 - e_1 * fs_2_3_14 * r_2 * h3_m2 + e_2 * fs_9_286_110 * h7_m4 - e_2 * fs_9_143_5 * h7_m2 + e_2 * fs_7_78_6 * r_2 * h5_m4 + e_2 * fs_5_39_2 * r_2 * h5_m2 + e_2 * fs_2_33_14 * r_4 * h3_m2;
    }

    // NOTE: the angular components are formed in 18 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph3_m3, ph3_m2, ph3_m1, ph3_p1, ph5_m3, ph5_m2, ph5_m1, ph5_p1, ph7_m3, ph7_m2, ph7_m1, ph7_p1, ab_2, pc_36, pc_37, pc_38 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

        const auto h3_m3 = ph3_m3[k];
        const auto h3_m2 = ph3_m2[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h5_m3 = ph5_m3[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h7_m3 = ph7_m3[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h7_m1 = ph7_m1[k];
        const auto h7_p1 = ph7_p1[k];

        pc_36[k] = - e_0 * fs_3_8_14 * h3_m3 + e_0 * fs_3_8_210 * h3_m1 - e_1 * fs_5_6_2 * h5_m3 - e_1 * fs_1_6_21 * h5_m1 + e_1 * fs_1_6_14 * r_2 * h3_m3 - e_1 * fs_1_6_210 * r_2 * h3_m1 + e_2 * fs_12_143_15 * h7_m3 - e_2 * fs_12_143_5 * h7_m1 + e_2 * fs_5_39_2 * r_2 * h5_m3 + e_2 * fs_1_39_21 * r_2 * h5_m1 - e_2 * fs_1_66_14 * r_4 * h3_m3 + e_2 * fs_1_66_210 * r_4 * h3_m1;

        pc_37[k] = - e_0 * fs_3_2_3 * h3_m2 - e_1 * fs_1_6_21 * h5_m2 + e_1 * fs_2_3_3 * r_2 * h3_m2 + e_2 * fs_3_143_210 * h7_m2 + e_2 * fs_1_39_21 * r_2 * h5_m2 - e_2 * fs_2_33_3 * r_4 * h3_m2;

        pc_38[k] = - e_0 * fs_15_4_2 * h3_p1 - e_1 * fs_1_6_5 * h5_p1 + e_1 * fs_5_3_2 * r_2 * h3_p1 + e_2 * fs_12_143_21 * h7_p1 + e_2 * fs_1_39_5 * r_2 * h5_p1 - e_2 * fs_5_33_2 * r_4 * h3_p1;
    }

    // NOTE: the angular components are formed in 18 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph3_0, ph3_p1, ph3_p2, ph3_p3, ph5_0, ph5_p1, ph5_p2, ph5_p3, ph7_0, ph7_p1, ph7_p2, ph7_p3, ab_2, pc_39, pc_40 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

        const auto h3_0 = ph3_0[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h3_p2 = ph3_p2[k];
        const auto h3_p3 = ph3_p3[k];
        const auto h5_0 = ph5_0[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h7_0 = ph7_0[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h7_p3 = ph7_p3[k];

        pc_39[k] = e_0 * fs_3_1_5 * h3_0 - e_0 * fs_3_2_3 * h3_p2 - e_1 * fs_1_6_5 * h5_0 - e_1 * fs_1_6_21 * h5_p2 - e_1 * fs_4_3_5 * r_2 * h3_0 + e_1 * fs_2_3_3 * r_2 * h3_p2 - e_2 * fs_21_143_5 * h7_0 + e_2 * fs_3_143_210 * h7_p2 + e_2 * fs_1_39_5 * r_2 * h5_0 + e_2 * fs_1_39_21 * r_2 * h5_p2 + e_2 * fs_4_33_5 * r_4 * h3_0 - e_2 * fs_2_33_3 * r_4 * h3_p2;

        pc_40[k] = e_0 * fs_3_8_210 * h3_p1 - e_0 * fs_3_8_14 * h3_p3 - e_1 * fs_1_6_21 * h5_p1 - e_1 * fs_5_6_2 * h5_p3 - e_1 * fs_1_6_210 * r_2 * h3_p1 + e_1 * fs_1_6_14 * r_2 * h3_p3 - e_2 * fs_12_143_5 * h7_p1 + e_2 * fs_12_143_15 * h7_p3 + e_2 * fs_1_39_21 * r_2 * h5_p1 + e_2 * fs_5_39_2 * r_2 * h5_p3 + e_2 * fs_1_66_210 * r_4 * h3_p1 - e_2 * fs_1_66_14 * r_4 * h3_p3;
    }

    // NOTE: the angular components are formed in 18 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph3_p2, ph3_p3, ph5_p2, ph5_p3, ph5_p4, ph5_p5, ph7_p2, ph7_p3, ph7_p4, ph7_p5, ph7_p6, ab_2, pc_41, pc_42, pc_43 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

        const auto h3_p2 = ph3_p2[k];
        const auto h3_p3 = ph3_p3[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h5_p4 = ph5_p4[k];
        const auto h5_p5 = ph5_p5[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h7_p3 = ph7_p3[k];
        const auto h7_p4 = ph7_p4[k];
        const auto h7_p5 = ph7_p5[k];
        const auto h7_p6 = ph7_p6[k];

        pc_41[k] = e_0 * fs_3_2_14 * h3_p2 - e_1 * fs_5_6_2 * h5_p2 - e_1 * fs_7_12_6 * h5_p4 - e_1 * fs_2_3_14 * r_2 * h3_p2 - e_2 * fs_9_143_5 * h7_p2 + e_2 * fs_9_286_110 * h7_p4 + e_2 * fs_5_39_2 * r_2 * h5_p2 + e_2 * fs_7_78_6 * r_2 * h5_p4 + e_2 * fs_2_33_14 * r_4 * h3_p2;

        pc_42[k] = e_0 * fs_3_4_42 * h3_p3 - e_1 * fs_7_12_6 * h5_p3 - e_1 * fs_1_4_30 * h5_p5 - e_1 * fs_1_3_42 * r_2 * h3_p3 - e_2 * fs_6_143_5 * h7_p3 + e_2 * fs_6_143_55 * h7_p5 + e_2 * fs_7_78_6 * r_2 * h5_p3 + e_2 * fs_1_26_30 * r_2 * h5_p5 + e_2 * fs_1_33_42 * r_4 * h3_p3;

        pc_43[k] = - e_1 * fs_1_4_30 * h5_p4 - e_2 * fs_3_286_22 * h7_p4 + e_2 * fs_3_143_143 * h7_p6 + e_2 * fs_1_26_30 * r_2 * h5_p4;
    }

    // NOTE: the angular components are formed in 18 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph3_m3, ph3_m2, ph3_m1, ph5_m5, ph5_m3, ph5_m2, ph5_m1, ph7_m7, ph7_m6, ph7_m5, ph7_m3, ph7_m2, ph7_m1, ab_2, pc_44, pc_45, pc_46 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

        const auto h3_m3 = ph3_m3[k];
        const auto h3_m2 = ph3_m2[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h5_m5 = ph5_m5[k];
        const auto h5_m3 = ph5_m3[k];
        const auto h5_m2 = ph5_m2[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h7_m7 = ph7_m7[k];
        const auto h7_m6 = ph7_m6[k];
        const auto h7_m5 = ph7_m5[k];
        const auto h7_m3 = ph7_m3[k];
        const auto h7_m2 = ph7_m2[k];
        const auto h7_m1 = ph7_m1[k];

        pc_44[k] = e_0 * fs_3_4_105 * h3_m3 + e_1 * fs_1_6_15 * h5_m3 - e_1 * fs_1_3_105 * r_2 * h3_m3 + e_2 * fs_3_286_2002 * h7_m7 + e_2 * fs_3_286_2 * h7_m3 - e_2 * fs_1_39_15 * r_2 * h5_m3 + e_2 * fs_1_33_105 * r_4 * h3_m3;

        pc_45[k] = e_0 * fs_9_4_7 * h3_m2 + e_1 * h5_m2 - e_1 * fs_1_1_7 * r_2 * h3_m2 + e_2 * fs_3_286_1430 * h7_m6 + e_2 * fs_3_286_10 * h7_m2 - e_2 * f_2_13 * r_2 * h5_m2 + e_2 * fs_1_11_7 * r_4 * h3_m2;

        pc_46[k] = e_0 * fs_3_4_35 * h3_m1 + e_1 * fs_1_6_15 * h5_m5 + e_1 * fs_1_3_14 * h5_m1 - e_1 * fs_1_3_35 * r_2 * h3_m1 + e_2 * fs_9_286_110 * h7_m5 + e_2 * fs_3_286_30 * h7_m1 - e_2 * fs_1_39_15 * r_2 * h5_m5 - e_2 * fs_2_39_14 * r_2 * h5_m1 + e_2 * fs_1_33_35 * r_4 * h3_m1;
    }

    // NOTE: the angular components are formed in 18 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph3_m3, ph3_m1, ph3_p2, ph5_m4, ph5_m3, ph5_m1, ph5_p2, ph7_m4, ph7_m3, ph7_m1, ph7_p2, ab_2, pc_47, pc_48, pc_49 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

        const auto h3_m3 = ph3_m3[k];
        const auto h3_m1 = ph3_m1[k];
        const auto h3_p2 = ph3_p2[k];
        const auto h5_m4 = ph5_m4[k];
        const auto h5_m3 = ph5_m3[k];
        const auto h5_m1 = ph5_m1[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h7_m4 = ph7_m4[k];
        const auto h7_m3 = ph7_m3[k];
        const auto h7_m1 = ph7_m1[k];
        const auto h7_p2 = ph7_p2[k];

        pc_47[k] = e_1 * h5_m4 + e_2 * fs_3_143_165 * h7_m4 - e_2 * f_2_13 * r_2 * h5_m4;

        pc_48[k] = e_0 * fs_3_8_2 * h3_m3 + e_0 * fs_3_8_30 * h3_m1 + e_1 * fs_1_3_14 * h5_m3 + e_1 * fs_5_6_3 * h5_m1 - e_1 * fs_1_6_2 * r_2 * h3_m3 - e_1 * fs_1_6_30 * r_2 * h3_m1 + e_2 * fs_3_143_105 * h7_m3 + e_2 * fs_3_143_35 * h7_m1 - e_2 * fs_2_39_14 * r_2 * h5_m3 - e_2 * fs_5_39_3 * r_2 * h5_m1 + e_2 * fs_1_66_2 * r_4 * h3_m3 + e_2 * fs_1_66_30 * r_4 * h3_m1;

        pc_49[k] = e_0 * fs_3_4_5 * h3_p2 + e_1 * fs_1_3_35 * h5_p2 - e_1 * fs_1_3_5 * r_2 * h3_p2 + e_2 * fs_9_143_14 * h7_p2 - e_2 * fs_2_39_35 * r_2 * h5_p2 + e_2 * fs_1_33_5 * r_4 * h3_p2;
    }

    // NOTE: the angular components are formed in 18 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph3_0, ph3_p1, ph3_p3, ph5_0, ph5_p1, ph5_p3, ph5_p4, ph5_p5, ph7_0, ph7_p1, ph7_p3, ph7_p4, ph7_p5, ab_2, pc_50, pc_51, pc_52 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

        const auto h3_0 = ph3_0[k];
        const auto h3_p1 = ph3_p1[k];
        const auto h3_p3 = ph3_p3[k];
        const auto h5_0 = ph5_0[k];
        const auto h5_p1 = ph5_p1[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h5_p4 = ph5_p4[k];
        const auto h5_p5 = ph5_p5[k];
        const auto h7_0 = ph7_0[k];
        const auto h7_p1 = ph7_p1[k];
        const auto h7_p3 = ph7_p3[k];
        const auto h7_p4 = ph7_p4[k];
        const auto h7_p5 = ph7_p5[k];

        pc_50[k] = - e_0 * fs_3_8_30 * h3_p1 + e_0 * fs_3_8_2 * h3_p3 - e_1 * fs_5_6_3 * h5_p1 + e_1 * fs_1_3_14 * h5_p3 + e_1 * fs_1_6_30 * r_2 * h3_p1 - e_1 * fs_1_6_2 * r_2 * h3_p3 - e_2 * fs_3_143_35 * h7_p1 + e_2 * fs_3_143_105 * h7_p3 + e_2 * fs_5_39_3 * r_2 * h5_p1 - e_2 * fs_2_39_14 * r_2 * h5_p3 - e_2 * fs_1_66_30 * r_4 * h3_p1 + e_2 * fs_1_66_2 * r_4 * h3_p3;

        pc_51[k] = e_0 * fs_3_4_35 * h3_0 + e_1 * fs_1_3_35 * h5_0 + e_1 * h5_p4 - e_1 * fs_1_3_35 * r_2 * h3_0 + e_2 * fs_3_143_35 * h7_0 + e_2 * fs_3_143_165 * h7_p4 - e_2 * fs_2_39_35 * r_2 * h5_0 - e_2 * f_2_13 * r_2 * h5_p4 + e_2 * fs_1_33_35 * r_4 * h3_0;

        pc_52[k] = e_0 * fs_3_4_35 * h3_p1 + e_1 * fs_1_3_14 * h5_p1 + e_1 * fs_1_6_15 * h5_p5 - e_1 * fs_1_3_35 * r_2 * h3_p1 + e_2 * fs_3_286_30 * h7_p1 + e_2 * fs_9_286_110 * h7_p5 - e_2 * fs_2_39_14 * r_2 * h5_p1 - e_2 * fs_1_39_15 * r_2 * h5_p5 + e_2 * fs_1_33_35 * r_4 * h3_p1;
    }

    // NOTE: the angular components are formed in 18 loops, grouped so that
    // the values a loop loads stay within the vector registers of the machine.
    // Only the prefactors, the harmonics and the squared distance are loaded by
    // more than one loop, every other value once.

#pragma omp simd aligned(pe_0, pe_1, pe_2, ph3_p2, ph3_p3, ph5_p2, ph5_p3, ph7_p2, ph7_p3, ph7_p6, ph7_p7, ab_2, pc_53, pc_54 : simd::cache_line_size())
    for (size_t k = 0; k < nvalues; k++)
    {
        const auto e_0 = pe_0[k];
        const auto e_1 = pe_1[k];
        const auto e_2 = pe_2[k];

        const auto r_2 = ab_2[k];
        const auto r_4 = r_2 * r_2;

        const auto h3_p2 = ph3_p2[k];
        const auto h3_p3 = ph3_p3[k];
        const auto h5_p2 = ph5_p2[k];
        const auto h5_p3 = ph5_p3[k];
        const auto h7_p2 = ph7_p2[k];
        const auto h7_p3 = ph7_p3[k];
        const auto h7_p6 = ph7_p6[k];
        const auto h7_p7 = ph7_p7[k];

        pc_53[k] = e_0 * fs_9_4_7 * h3_p2 + e_1 * h5_p2 - e_1 * fs_1_1_7 * r_2 * h3_p2 + e_2 * fs_3_286_10 * h7_p2 + e_2 * fs_3_286_1430 * h7_p6 - e_2 * f_2_13 * r_2 * h5_p2 + e_2 * fs_1_11_7 * r_4 * h3_p2;

        pc_54[k] = e_0 * fs_3_4_105 * h3_p3 + e_1 * fs_1_6_15 * h5_p3 - e_1 * fs_1_3_105 * r_2 * h3_p3 + e_2 * fs_3_286_2 * h7_p3 + e_2 * fs_3_286_2002 * h7_p7 - e_2 * fs_1_39_15 * r_2 * h5_p3 + e_2 * fs_1_33_105 * r_4 * h3_p3;
    }

    // NOTE: the values of a combination of angular components are stored as one
    // row of nvalues columns, with the component on bra side running slowest. The
    // rows which the symmetry relates to an already formed one are copied from it.

    const size_t sources[55] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54};

    for (size_t n = 0; n < 55; n++)
    {
        if (sources[n] != n) std::copy(values + sources[n] * nvalues, values + (sources[n] + 1) * nvalues, values + n * nvalues);
    }
}

}  // namespace simdt2ceri
